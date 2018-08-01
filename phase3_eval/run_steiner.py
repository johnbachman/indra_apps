import networkx as nx
from indra.explanation.steiner_approx import \
                    MemoizedShortestPath, directed_steiner_approx
import make_stmts_for_checking as make_stmts
import itertools
import pickle

def draw(g, filename):
    ag = nx.nx_agraph.to_agraph(g)
    ag.draw(filename, prog='dot')

def load_sif():
    edges = []
    with open('output/korkut_model_pysb_pysb.sif', 'rt') as f:
        for line in f.readlines():
            u, pol, v = line.strip().split(' ')
            edges.append((u, v))
    g = nx.DiGraph()
    g.add_edges_from(edges)
    return g

def load_im():
    with open('influence_map.pkl', 'rb') as f:
        im = pickle.load(f)
    # Flatten the multidigraph
    edges = set()
    for edge in im.edges():
        edges.add(edge)
    dg = nx.DiGraph()
    dg.add_edges_from(edges)
    return dg

def load_nx():
    with open('output/korkut_full_nx_direct.pkl', 'rb') as f:
        g = pickle.load(f)
    return g

# Add edges linking families to individual genes

"""
print("Loading influence map")
g = load_im()
braf_edges = []
for node in g.nodes():
    if node.startswith('BRAF'):
        braf_edges.append(('BRAF', node))
g.add_edges_from(braf_edges)
"""

def add_parent_node(g, parent, children):
    edges = [(parent, child) for child in children]
    g.add_edges_from(edges)
    return g

"""

# Set the parameters
dist = MemoizedShortestPath(g)
num_terminals_to_connect = 3
max_levels = 5

#g = add_parent_node(g, 'AKT', ['AKT1', 'AKT2', 'AKT3'])
#root = 'AKT'
#terminals = ['GSK3B', 'RPS6']

#g = add_parent_node(g, 'SFK', ['SRC', 'BLK', 'FYN'])
#root = 'BLK'
#terminals = ['EIF4EBP1', 'CHEK2']

root = 'EGFR'
terminals = ['KRAS', 'RPS6', 'MAPK1']

tree = directed_steiner_approx(g, terminals, num_terminals_to_connect, root,
                               max_levels, dist)

print("Getting shortest paths")
tree_edges = []
for u, v in tree.edges:
    sp = nx.shortest_path(g, u, v)
    for i, node in enumerate(sp[:-1]):
        tree_edges.append((node, sp[i+1]))

steiner_graph = nx.DiGraph()
steiner_graph.add_edges_from(tree_edges)

draw(steiner_graph, 'steiner_graph.pdf')
"""
if __name__ == '__main__':
    g = load_nx()
    with open('edges.tsv', 'wt') as f:
        for u, v in g.edges():
            f.write('%s\t%s\t%s\t%s\n' % (str(u), str(v), '1', 'D'))

