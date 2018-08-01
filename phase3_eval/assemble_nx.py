from indra.util import _require_python3
from os.path import join as pjoin
from indra.assemblers import GraphAssembler
import indra.tools.assemble_corpus as ac
import networkx as nx
import pickle

def assemble_nx(stmts, out_file_prefix, network_type):
    """Return a CX assembler."""
    stmts = ac.filter_belief(stmts, 0.95)
    stmts = ac.filter_top_level(stmts)

    if network_type == 'direct':
        stmts = ac.filter_direct(stmts)

    out_file = '%s_%s.pkl' % (out_file_prefix, network_type)

    ga = GraphAssembler()
    ga.add_statements(stmts)
    ga.make_model()
    nx_graph = nx.DiGraph(ga.graph)
    with open(out_file, 'wb') as f:
        pickle.dump(nx_graph, f)
