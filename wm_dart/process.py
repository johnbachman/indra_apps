import os
import glob
import copy
import json
import yaml
import logging
from indra.sources import eidos, hume
from indra.sources.eidos import migration_table_processor
from indra.tools.live_curation import Corpus
from indra.tools import assemble_corpus as ac
from indra.sources.eidos.reader import EidosReader
from indra.belief.wm_scorer import get_eidos_scorer
from indra.preassembler.custom_preassembly import *
from indra.statements import Event, Influence, Association
from indra.preassembler.hierarchy_manager import YamlHierarchyManager
from indra.preassembler.make_eidos_hume_ontologies import eidos_ont_url, \
    load_yaml_from_url, rdf_graph_from_yaml


logger = logging.getLogger()
data_path = os.path.join(os.path.expanduser('~'), 'data', 'wm', 'dart')


def load_eidos():
    logger.info('Loading Eidos statements')
    fnames = glob.glob(os.path.join(data_path, 'eidos/jsonldDir/*.jsonld'))

    stmts = []
    for fname in fnames:
        doc_id = os.path.basename(fname).split('.')[0]
        ep = eidos.process_json_file(fname)
        fix_provenance(ep.statements, doc_id)
        stmts += ep.statements
    logger.info(f'Loaded {len(stmts)} statements from Eidos')
    return stmts


def load_hume():
    logger.info('Loading Hume statements')
    fnames = glob.glob(os.path.join(data_path,
                                    'hume/wm_dart.082919.v3.json-ld'))

    stmts = []
    for fname in fnames:
        hp = hume.process_jsonld_file(fname)
        stmts += hp.statements
    logger.info(f'Loaded {len(stmts)} statements from Hume')
    return stmts


def load_migration_spreadsheets():
    fname = os.path.join(os.pardir, 'wm_migration',
                         'grounded CAG links - New Ontology.xlsx')
    stmts = migration_table_processor.process_workbook(fname)
    return stmts


def fix_provenance(stmts, doc_id):
    """Move the document identifiers in evidences."""
    for stmt in stmts:
        for ev in stmt.evidence:
            prov = ev.annotations['provenance'][0]['document']
            prov['@id'] = doc_id


def remove_namespaces(stmts, namespaces):
    """Remove unnecessary namespaces from Concept grounding."""
    for stmt in stmts:
        for agent in stmt.agent_list():
            for namespace in namespaces:
                if namespace in copy.deepcopy(agent.db_refs):
                    agent.db_refs.pop(namespace, None)


def remove_raw_grounding(stmts):
    """Remove the raw_grounding annotation to decrease output size."""
    for stmt in stmts:
        for ev in stmt.evidence:
            if not ev.annotations:
                continue
            agents = ev.annotations.get('agents')
            if not agents:
                continue
            if 'raw_grounding' in agents:
                agents.pop('raw_grounding', None)


def get_events(stmts):
    """Return a list of all standalone events from a list of statements."""
    events = []
    for stmt in stmts:
        stmt = copy.deepcopy(stmt)
        if isinstance(stmt, Influence):
            for member in [stmt.subj, stmt.obj]:
                member.evidence = stmt.evidence[:]
                # Remove the context since it may be for the other member
                for ev in member.evidence:
                    ev.context = None
                events.append(member)
        elif isinstance(stmt, Association):
            for member in stmt.members:
                member.evidence = stmt.evidence[:]
                # Remove the context since it may be for the other member
                for ev in member.evidence:
                    ev.context = None
                events.append(member)
        elif isinstance(stmt, Event):
            events.append(stmt)
    return events


def get_non_events(stmts):
    """Return a list of statements that aren't Events"""
    return [st for st in stmts if not isinstance(st, Event)]


def check_event_context(events):
    for event in events:
        if not event.context and event.evidence[0].context:
            assert False, ('Event context issue', event, event.evidence)
        ej = event.to_json()
        if 'context' not in ej and 'context' in ej['evidence'][0]:
            assert False, ('Event context issue', event, event.evidence)


def reground_stmts(stmts):
    ont_manager = _make_un_ontology()
    eidos_reader = EidosReader()
    # Send the latest ontology and list of concept texts to Eidos
    yaml_str = yaml.dump(ont_manager.yaml_root)
    concepts = []
    for stmt in stmts:
        for concept in stmt.agent_list():
            concept_txt = concept.db_refs.get('TEXT')
            concepts.append(concept_txt)
    groundings = eidos_reader.reground_texts(concepts, yaml_str)
    # Update the corpus with new groundings
    idx = 0
    for stmt in stmts:
        for concept in stmt.agent_list():
            concept.db_refs['UN'] = groundings[idx]
            idx += 1
    return stmts


def _make_un_ontology():
    return YamlHierarchyManager(load_yaml_from_url(eidos_ont_url),
                                rdf_graph_from_yaml, True)


if __name__ == '__main__':
    mig_stmts = load_migration_spreadsheets()
    eidos_stmts = load_eidos()
    hume_stmts = load_hume()
    stmts = eidos_stmts + hume_stmts + mig_stmts
    reground_stmts(stmts)
    remove_namespaces(stmts, ['WHO', 'MITRE12', 'WM'])

    events = get_events(stmts)
    check_event_context(events)
    non_events = get_non_events(stmts)
    scorer = get_eidos_scorer()

    funs = {
        'grounding': (None, None),
        'location': (location_matches, location_refinement),
        'location_and_time': (location_time_matches,
                              location_time_refinement)
    }

    for key, (matches_fun, refinement_fun) in funs.items():
        assembled_non_events = ac.run_preassembly(non_events,
                                                  belief_scorer=scorer,
                                                  matches_fun=matches_fun,
                                                  refinement_fun=refinement_fun)
        assembled_events = ac.run_preassembly(events, belief_scorer=scorer,
                                              matches_fun=matches_fun,
                                              refinement_fun=refinement_fun)
        check_event_context(assembled_events)
        assembled_stmts = assembled_non_events + assembled_events
        remove_raw_grounding(assembled_stmts)
        corpus = Corpus(assembled_stmts, raw_statements=stmts)
        corpus.s3_put('dart-20190910-stmts-%s' % key)
        sj = stmts_to_json(assembled_stmts, matches_fun=matches_fun)
        with open(os.path.join(data_path,
                  'dart-20190910-stmts-%s.json' % key), 'w') as fh:
            json.dump(sj, fh, indent=1)
