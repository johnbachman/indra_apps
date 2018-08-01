import json
import pickle
from indra.util import write_unicode_csv
from indra.assemblers import PysbAssembler, EnglishAssembler, CyJSAssembler
from indra.explanation.model_checker import ModelChecker
import indra.tools.assemble_corpus as ac
from assemble_pysb import set_context, add_observables
from collections import defaultdict
from indra.assemblers import pysb_assembler as pa

from indra.statements import Dephosphorylation

import itertools
import process_data
import make_stmts_for_checking as make_stmts

from indra.util.graph import signed_edges_to_nodes, signed_node

def get_input_rules(stmt, mc):
    # Collect all input rules for this drug
    ag = stmt.agent_list()[0]
    enz_mps = list(pa.grounded_monomer_patterns(mc.model, ag))
    rules = set()
    for mp in enz_mps:
        rules |= mc._get_input_rules(mp)
    return rules

def dump_edge_file(im, filename):
    weight = 0.99
    with open(filename, 'wt') as f:
        for u, v in im.edges():
            f.write('%s\t%s\t%s\t%s\n' % (u, v, weight, 'D'))

def dump_prize_file(prize_dict, filename):
    with open(filename, 'wt') as f:
        for node, prize in prize_dict.items():
            f.write('%s\t%s\n' % (str(node), prize))



if __name__ == '__main__':
    print("Processing data")

    data = process_data.read_data(process_data.data_file)
    data_genes = process_data.get_all_gene_names(data)
    ab_map = process_data.get_antibody_map(data)

    print('Loading data statements.')
    data_stmts, data_values = make_stmts.run(dec_thresh=0.8, inc_thresh=1.2)
    all_data_stmts = [values.values() for values in data_stmts.values()]
    all_data_stmts = itertools.chain.from_iterable(all_data_stmts)
    all_data_stmts = list(itertools.chain.from_iterable(all_data_stmts))

    agent_obs = list(itertools.chain.from_iterable(ab_map.values()))
    # Here we need to cross-reference the antbody map with the data values
    agent_data = {}
    for drug_name, values in data_values.items():
        agent_data[drug_name] = {}
        for ab_name, value in values.items():
            agents = ab_map[ab_name]
            for agent in agents:
                agent_data[drug_name][agent] = value

    base_stmts = ac.load_statements('output/korkut_model_pysb_before_pa.pkl')
    for st in base_stmts:
        st.uuid = str(st.uuid)

    """
    # Merge the sources of statements
    # stmts = manual_stmts + base_stmts
    stmts = base_stmts
    #stmts = manual_stmts

    # Assemble model
    pa = PysbAssembler()
    pa.add_statements(stmts)
    model = pa.make_model()

    with open('korkut_pysb.pkl', 'wb') as f:
        print("Pickling PySB model")
        pickle.dump(pa.model, f)
    """

    with open('output/korkut_pysb_model.pkl', 'rb') as f:
        print("Unpickling PySB model")
        model = pickle.load(f)
    stmts = base_stmts

    mc = ModelChecker(model, all_data_stmts, agent_obs)
    mc.prune_influence_map()

    # NOTE NOTE NOTE
    signed = True
    # NOTE NOTE NOTE

    signed_im = signed_edges_to_nodes(mc.get_im())

    # Add edges between all drugs and the relevant input rules
    for drug_name, ab_dict in data_stmts.items():
        input_rules = set()
        first_stmt = True
        for ab, stmt_list in ab_dict.items():
            for stmt in stmt_list:
                if first_stmt:
                    input_rules |= get_input_rules(stmt, mc)
                    # Add a dummy node to the input rules
                    if signed:
                        edges = [(drug_name, signed_node(r, -1))
                                 for r in input_rules]
                        signed_im.add_edges_from(edges)
                    else:
                        edges = [(drug_name, r) for r in input_rules]
                        mc.get_im().add_edges_from(edges)
                    first_stmt = False
    if signed:
        dump_edge_file(signed_im, 'im_edges_signed.tsv')
    else:
        dump_edge_file(mc.get_im(), 'im_edges.tsv')

    # Now compile prize files for each drug from the observables
    for drug_name, ab_dict in data_stmts.items():
        agent_values = agent_data[drug_name]
        # Compile prize file from antibody list and values
        prizes = []
        obs_prizes = defaultdict(lambda: 0)
        for ab, stmt_list in ab_dict.items():
            value = data_values[drug_name][ab]
            prize_value = abs(value - 1) * 10
            if signed:
                prize_sign = -1 if value < 1 else 1
            for stmt in stmt_list:
                obs_names = mc.stmt_to_obs.get(stmt)
                if obs_names:
                    for obs_name in obs_names:
                        if signed:
                            obs_prizes[signed_node(obs_name, prize_sign)] = \
                              max(obs_prizes[signed_node(obs_name, prize_sign)],
                                  prize_value)
                        else:
                            obs_prizes[obs_name] = \
                                    max(obs_prizes[obs_name], prize_value)
                else:
                    print("No observables for statement %s" % stmt)
        # Add the drug target with a large prize
        obs_prizes[drug_name] = 100
        dump_prize_file(obs_prizes, '%s_prizes.tsv' % drug_name)

# Next step:
# Create signed digraph from the original digraph

# Filter digraph against data to remove nodes that conflict with the sign;
# Add negative prizes to penalize nodes that conflict with the sign, or;
# Simply filter nodes that conflict with the sign

# In both cases this will involve finding

# Add negative prizes to nodes that are measured but unchanged
# This would involve finding
"""
write_unicode_csv('model_check_results.csv', results)
path_stmts = get_path_stmts(results, model, base_stmts)
path_genes = get_path_genes(path_stmts)
#make_english_output(results, model, base_stmts)
#make_cyjs_network(results, model, base_stmts)
paths_json = export_json(results, model, base_stmts)
"""
