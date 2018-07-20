import dendropy as dy
import sys
import pandas as pd
from optparse import OptionParser
from abc import ABC, abstractmethod
import heapq
import math

# Glossary
# OLS: Ordinary Least Squares
# FM: Fitch-Margoliash --Least Squares with Variance scaling
# BE: Beyer --Least Squares with Standard Deviation scaling

class Algorithm(ABC):
    def __init__(self, tree):
        self.tree = tree

    @abstractmethod
    def placement_per_edge(self):
        pass

    @abstractmethod
    def error_per_edge(self, edge):
        pass

    def placement(self, query_name, selection_name):
        if selection_name == "HYBRID":
            sm = heapq.nsmallest(math.floor(math.log2(len(self.tree.leaf_nodes()))),
                                 tree.postorder_edge_iter(filter_fn = lambda e: e.head_node != tree.seed_node),
                                 key = lambda e: self.error_per_edge(e))
            placed_edge = min(sm, key = lambda e: e.x_1)
        elif selection_name == "ME":
            placed_edge = min(tree.postorder_edge_iter(filter_fn = lambda e: e.head_node != tree.seed_node),
                              key=lambda e: e.x_1)
        else: # selection_name == "MLSE"
            placed_edge = min(tree.postorder_edge_iter(filter_fn = lambda e: e.head_node != tree.seed_node),
                              key=lambda e: self.error_per_edge(e))
        insert(placed_edge, query_name)
        return self.tree

class OLS(Algorithm):

    # computes all alternative OLS placements
    def placement_per_edge(self):
        for e in tree.postorder_edge_iter():
            if e.head_node == tree.seed_node:
                continue
            a_11 = e.R + e.S
            a_12 = e.R - e.S
            a_21 = a_12
            a_22 = a_11
            c_1 = e.RD + e.SD - e.length * e.S - e.Rd - e.Sd
            c_2 = e.RD - e.SD + e.length * e.S - e.Rd + e.Sd
            x_1, x_2 = solve2_2(a_11, a_12, a_21, a_22, c_1, c_2)
            e.x_1 = x_1
            e.x_2 = x_2

    # computes OLS error (Q value) for a given edge
    def error_per_edge(self, edge):
        A = edge.RD2 + edge.SD2
        B = 2 * (edge.x_1 + edge.x_2) * edge.Rd + \
            2 * (edge.length + edge.x_1 - edge.x_2) * edge.Sd
        C = (edge.x_1 + edge.x_2) ** 2 * edge.R + \
            (edge.length + edge.x_1 - edge.x_2) ** 2 * edge.S
        D = -2 * (edge.x_1 + edge.x_2) * edge.RD - \
            2 * (edge.length + edge.x_1 - edge.x_2) * edge.SD
        E = -2 * edge.RDd - 2 * edge.SDd
        F = edge.Rd2 + edge.Sd2
        return A + B + C + D + E + F

class FM(Algorithm):

    # computes all alternative Fitch-Margoliash placements
    def placement_per_edge(self):
        for e in tree.postorder_edge_iter():
            if e.head_node == tree.seed_node:
                continue
            a_11 = e.R1_D2 + e.S1_D2
            a_12 = e.R1_D2 - e.S1_D2
            a_21 = a_12
            a_22 = a_11
            c_1 = e.R1_D + e.S1_D - e.length * e.S1_D2 - e.Rd_D2 - e.Sd_D2
            c_2 = e.R1_D - e.S1_D + e.length * e.S1_D2 - e.Rd_D2 + e.Sd_D2
            x_1, x_2 = solve2_2(a_11, a_12, a_21, a_22, c_1, c_2)
            e.x_1 = x_1
            e.x_2 = x_2

    # computes Fitch-Margoliash error (Q value) for a given edge
    def error_per_edge(self, edge):
        A = edge.R + edge.S
        B = 2 * (edge.x_1 + edge.x_2) * edge.Rd_D2 + \
            2 * (edge.length + edge.x_1 - edge.x_2) * edge.Sd_D2
        C = (edge.x_1 + edge.x_2) ** 2 * edge.R1_D2 + \
            (edge.length + edge.x_1 - edge.x_2) ** 2 * edge.S1_D2
        D = -2 * (edge.x_1 + edge.x_2) * edge.R1_D - \
            2 * (edge.length + edge.x_1 - edge.x_2) * edge.S1_D
        E = -2 * edge.Rd_D - 2 * edge.Sd_D
        F = edge.Rd2_D2 + edge.Sd2_D2
        return A + B + C + D + E + F

class BE(Algorithm):

    # computes all alternative BE placements
    def placement_per_edge(self):
        for e in tree.postorder_edge_iter():
            if e.head_node == tree.seed_node:
                continue
            a_11 = e.R1_D + e.S1_D
            a_12 = e.R1_D - e.S1_D
            a_21 = a_12
            a_22 = a_11
            c_1 = e.R + e.S - e.length * e.S1_D - e.Rd_D - e.Sd_D
            c_2 = e.R - e.S + e.length * e.S1_D - e.Rd_D + e.Sd_D
            x_1, x_2 = solve2_2(a_11, a_12, a_21, a_22, c_1, c_2)
            e.x_1 = x_1
            e.x_2 = x_2

    # computes BE error (Q value) for a given edge
    def error_per_edge(self, edge):
        A = edge.RD + edge.SD
        B = 2 * (edge.x_1 + edge.x_2) * edge.Rd_D + \
            2 * (edge.length + edge.x_1 - edge.x_2) * edge.Sd_D
        C = (edge.x_1 + edge.x_2) ** 2 * edge.R1_D + \
            (edge.length + edge.x_1 - edge.x_2) ** 2 * edge.S1_D
        D = -2 * (edge.x_1 + edge.x_2) * edge.R - \
            2 * (edge.length + edge.x_1 - edge.x_2) * edge.S
        E = -2 * edge.Rd - 2 * edge.Sd
        F = edge.Rd2_D + edge.Sd2_D
        return A + B + C + D + E + F

# inserts query into reference tree at given edge.

def insert(placed_edge, query_name):
    tailn = placed_edge.tail_node
    headn = placed_edge.head_node
    tailn.remove_child(headn)
    nn = dy.Node()
    nn.add_child(headn)
    qry = dy.Node(taxon=dy.Taxon(query_name))
    nn.add_child(qry)
    qry.edge_length = placed_edge.x_1
    tailn.add_child(nn)
    if placed_edge.head_node in list(master_edge.head_node.ancestor_iter()) or master_edge == placed_edge:
        nn.edge_length = placed_edge.length - max(placed_edge.x_2, 0)
        headn.edge_length = max(placed_edge.x_2, 0)
    else:
        nn.edge_length = max(placed_edge.x_2, 0)
        headn.edge_length = placed_edge.length - max(placed_edge.x_2, 0)

# solves two by two Ax=c linear system
def solve2_2(a_11, a_12, a_21, a_22, c_1, c_2):
    det = 1 / (a_11 * a_22 - a_12 * a_21)
    assert det is not 0
    x_1 = (a_22 * c_1 - a_12 * c_2) * det
    x_2 = (- a_21 * c_1 + a_11 * c_2) * det
    return (x_1, x_2)

def dfs_S_values(edge, downstream):
    if downstream.is_leaf():
        edge.S = 1
        edge.SDd = 0
        edge.Sd = 0
        edge.Sd_D = 0
        edge.Sd_D2 = 0
        edge.Sd2 = 0
        edge.Sd2_D = 0
        edge.Sd2_D2 = 0
        edge.SD = obs_dist[downstream.taxon.label]
        edge.SD2 = edge.SD * edge.SD
        edge.S1_D = 1/edge.SD
        edge.S1_D2 = 1/(edge.SD*edge.SD)


    else:
        inc = list(downstream.incident_edges())
        inc = filter(lambda x: x.length != None and x != edge, inc)
        edge.S, edge.SDd, edge.Sd, edge.Sd_D, edge.Sd_D2, edge.Sd2, edge.Sd2_D, edge.Sd2_D2, edge.SD2, edge.SD, edge.S1_D,  edge.S1_D2 = 12*[0]
        for d1  in inc:
            d1tips = [d1.head_node, d1.tail_node]
            d1tips.remove(downstream)
            if len(d1tips) == 0:
                print("error")
            dfs_S_values(d1, d1tips[0])
            edge.S += d1.S
            edge.SDd += d1.length * d1.SD + d1.SDd
            edge.Sd += d1.S * d1.length + d1.Sd
            edge.Sd_D += d1.length * d1.S1_D + d1.Sd_D
            edge.Sd_D2 += d1.length * d1.S1_D2 + d1.Sd_D2
            edge.Sd2 += d1.S * d1.length * d1.length + d1.Sd2 + 2 * d1.length * d1.Sd
            edge.Sd2_D += d1.S1_D * d1.length * d1.length + d1.Sd2_D + 2 * d1.length * d1.Sd_D
            edge.Sd2_D2 += d1.S1_D2 * d1.length * d1.length + d1.Sd2_D2 + 2 * d1.length * d1.Sd_D2
            edge.SD2 += d1.SD2
            edge.SD += d1.SD
            edge.S1_D += d1.S1_D
            edge.S1_D2 += d1.S1_D2


def dfs_R_values(edge, u1, upstream, downstream):
    if upstream.is_leaf():
        edge.R = 1
        edge.RDd = 0
        edge.Rd = 0
        edge.Rd_D = 0
        edge.Rd_D2 = 0
        edge.Rd2 = 0
        edge.Rd2_D = 0
        edge.Rd2_D2 = 0
        edge.RD = obs_dist[upstream.taxon.label]
        edge.RD2 = edge.RD * edge.RD
        edge.R1_D = 1/edge.RD
        edge.R1_D2 = 1/(edge.RD*edge.RD)

    else:
        edge.R = u1.R
        edge.RDd = u1.RD * u1.length + u1.RDd
        edge.Rd = u1.R * u1.length + u1.Rd
        edge.Rd_D = u1.length * u1.R1_D + u1.Rd_D
        edge.Rd_D2 = u1.length * u1.R1_D2 + u1.Rd_D2
        edge.Rd2 = u1.R * u1.length * u1.length + u1.Rd2 + 2 * u1.length * u1.Rd
        edge.Rd2_D = u1.R1_D * u1.length * u1.length + u1.Rd2_D + 2 * u1.length * u1.Rd_D
        edge.Rd2_D2 = u1.R1_D2 * u1.length * u1.length + u1.Rd2_D2 + 2 * u1.length * u1.Rd_D2
        edge.RD2 = u1.RD2
        edge.RD = u1.RD
        edge.R1_D = u1.R1_D
        edge.R1_D2 = u1.R1_D2
        inc = list(upstream.incident_edges())
        u2 = list(filter(lambda x: x.length != None and x != edge and x != u1, inc))
        if len(u2) == 1:
            edge.R += u2[0].S
            edge.RDd += u2[0].SD * u2[0].length + u2[0].SDd
            edge.Rd += u2[0].S * u2[0].length + u2[0].Sd
            edge.Rd_D += u2[0].length * u2[0].S1_D + u2[0].Sd_D
            edge.Rd_D2 += u2[0].length * u2[0].S1_D2 + u2[0].Sd_D2
            edge.Rd2 += u2[0].S * u2[0].length * u2[0].length + u2[0].Sd2 + 2 * u2[0].length * u2[0].Sd
            edge.Rd2_D += u2[0].S1_D * u2[0].length * u2[0].length + u2[0].Sd2_D + 2 * u2[0].length * u2[0].Sd_D
            edge.Rd2_D2 += u2[0].S1_D2 * u2[0].length * u2[0].length + u2[0].Sd2_D2 + 2 * u2[0].length * u2[0].Sd_D2
            edge.RD2 += u2[0].SD2
            edge.RD += u2[0].SD
            edge.R1_D += u2[0].S1_D
            edge.R1_D2 += u2[0].S1_D2

    if not downstream.is_leaf():
        inc = list(downstream.incident_edges())
        inc = filter(lambda x: x.length != None and x != edge, inc)
        for d1 in inc:
            d1tips = [d1.head_node, d1.tail_node]
            d1tips.remove(downstream)
            dfs_R_values(d1, edge, downstream, d1tips[0])


if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-t", "--tree", dest="tree_fp",
                      help="path to the reference tree", metavar="FILE")
    parser.add_option("-d", "--distances", dest="dist_fp",
                      help="path to the table of observed distances", metavar="FILE")
    parser.add_option("-q", "--query", dest="query_name",
                      help="name of the query taxon", metavar="NAME")
    parser.add_option("-a", "--algo", dest="algo_name", default = "OLS",
                      help="name of the algorithm (OLS, FM, or BE)", metavar="ALGO")
    parser.add_option("-s", "--selection", dest="selection_name", default = "MLSE",
                      help="name of the placement selection criteria (MLSE, ME, or HYBRID", metavar="CRITERIA")


    (options, args) = parser.parse_args()
    tree_fp = options.tree_fp
    dist_fp = options.dist_fp
    algo_name = options.algo_name
    query_name = options.query_name
    selection_name = options.selection_name
    hybrid = options.hybrid

    df = pd.read_csv(dist_fp, sep="\s+", names = ["taxa", "dist"] , dtype = {"taxa": object}, header = None)
    obs_dist = pd.Series(df.dist.values,index=df.taxa).to_dict()
    f = open(tree_fp)
    tree = dy.Tree.get_from_string(f.readline(), schema='newick')
    master_edge = next(tree.postorder_edge_iter())
    for l in tree.leaf_node_iter():
        if l.taxon.label not in obs_dist:
            raise ValueError('Taxon {} should be in distances table.'.format(l.taxon.label))
        if obs_dist[l.taxon.label] == 0:
            l.edge.x_1 = 0
            l.edge.x_2 = 0
            insert(l.edge, query_name)
            tree.write(file=sys.stdout, schema="newick")
            exit(0)
    dfs_S_values(master_edge, master_edge.tail_node)
    dfs_R_values(master_edge, None, master_edge.head_node, master_edge.tail_node)
    if algo_name == "BE":
        alg = BE(tree)
    elif algo_name == "FM":
        alg = FM(tree)
    else:
        alg = OLS(tree)
    alg.placement_per_edge()
    output_tree = alg.placement(query_name, selection_name)
    output_tree.write(file = sys.stdout, schema = "newick")
