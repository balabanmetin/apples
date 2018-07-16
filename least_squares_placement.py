import dendropy as dy
import sys
import pandas as pd
from optparse import OptionParser
from abc import ABC, abstractmethod

# Glossary
# OLS: Ordinary Least Squares
# FM: Fitch-Margoliash --Least Squares with Variance scaling
# BE: Beyer --Least Squares with Standard Deviation scaling


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
        edge.Sd = 0
        edge.Sd2 = 0
        edge.SD = obs_dist[downstream.taxon.label]
        edge.SD2 = edge.SD * edge.SD
        edge.SDd = 0
        edge.S1_D = 1/edge.SD
        edge.Sd_D = 0
        edge.S1_D2 = 1/(edge.SD*edge.SD)
        edge.Sd_D2 = 0
        edge.Sd2_D2 = 0



    else:
        inc = list(downstream.incident_edges())
        inc = filter(lambda x: x.length != None and x != edge, inc)
        edge.S, edge.Sd, edge.Sd2, edge.SD, edge.SD2, edge.SDd, edge.S1_D, edge.Sd_D, edge.S1_D2, edge.Sd_D2, edge.Sd2_D2 = 11*[0]
        for d1  in inc:
            d1tips = [d1.head_node, d1.tail_node]
            d1tips.remove(downstream)
            if len(d1tips) == 0:
                print("error")
            dfs_S_values(d1, d1tips[0])
            edge.S += d1.S
            edge.Sd += d1.S * d1.length + d1.Sd
            edge.Sd2 += d1.S * d1.length * d1.length + d1.Sd2 + 2 * d1.length * d1.Sd
            edge.SD += d1.SD
            edge.SD2 += d1.SD2
            edge.SDd += d1.length * d1.SD + d1.SDd
            edge.S1_D += d1.S1_D
            edge.Sd_D += d1.length * d1.S1_D + d1.Sd_D
            edge.S1_D2 += d1.S1_D2
            edge.Sd_D2 += d1.length * d1.S1_D2 + d1.Sd_D2
            edge.Sd2_D2 += d1.S1_D2 * d1.length * d1.length + d1.Sd2_D2 + 2 * d1.length * d1.Sd_D2

def dfs_R_values(edge, u1, upstream, downstream):
    if upstream.is_leaf():
        edge.R = 1
        edge.Rd = 0
        edge.Rd2 = 0
        edge.RD = obs_dist[upstream.taxon.label]
        edge.RD2 = edge.RD * edge.RD
        edge.RDd = 0
        edge.R1_D = 1/edge.RD
        edge.Rd_D = 0
        edge.R1_D2 = 1/(edge.RD*edge.RD)
        edge.Rd_D2 = 0
        edge.Rd2_D2 = 0

    else:
        edge.R = u1.R
        edge.Rd = u1.R * u1.length + u1.Rd
        edge.Rd2 = u1.R * u1.length * u1.length + u1.Rd2 + 2 * u1.length * u1.Rd
        edge.RD = u1.RD
        edge.RD2 = u1.RD2
        edge.RDd = u1.RD * u1.length + u1.RDd
        edge.R1_D = u1.R1_D
        edge.Rd_D = u1.length * u1.R1_D + u1.Rd_D
        edge.R1_D2 = u1.R1_D2
        edge.Rd_D2 = u1.length * u1.R1_D2 + u1.Rd_D2
        edge.Rd2_D2 = u1.R1_D2 * u1.length * u1.length + u1.Rd2_D2 + 2 * u1.length * u1.Rd_D2
        inc = list(upstream.incident_edges())
        u2 = list(filter(lambda x: x.length != None and x != edge and x != u1, inc))
        if len(u2) == 1:
            edge.R += u2[0].S
            edge.Rd += u2[0].S * u2[0].length + u2[0].Sd
            edge.Rd2 += u2[0].S * u2[0].length * u2[0].length + u2[0].Sd2 + 2 * u2[0].length * u2[0].Sd
            edge.RD += u2[0].SD
            edge.RD2 += u2[0].SD2
            edge.RDd += u2[0].SD * u2[0].length + u2[0].SDd
            edge.R1_D += u2[0].S1_D
            edge.Rd_D += u2[0].length * u2[0].S1_D + u2[0].Sd_D
            edge.R1_D2 += u2[0].S1_D2
            edge.Rd_D2 += u2[0].length * u2[0].S1_D2 + u2[0].Sd_D2
            edge.Rd2_D2 += u2[0].S1_D2 * u2[0].length * u2[0].length + u2[0].Sd2_D2 + 2 * u2[0].length * u2[0].Sd_D2

    if not downstream.is_leaf():
        inc = list(downstream.incident_edges())
        inc = filter(lambda x: x.length != None and x != edge, inc)
        for d1 in inc:
            d1tips = [d1.head_node, d1.tail_node]
            d1tips.remove(downstream)
            dfs_R_values(d1, edge, downstream, d1tips[0])

# computes all alternative OLS placements
def ols_parameters_per_edge(tree):
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
def ols_error(edge):
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

# computes all alternative Fitch-Margoliash placements
def fm_parameters_per_edge(tree):
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
def fm_error(edge):
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


def placement(tree, name, minimum_evolution):
    min = 2000000000
    placed_edge = None
    for e in tree.postorder_edge_iter():
        if e.head_node == tree.seed_node:
            continue
        if minimum_evolution:
            err = e.x_1
        else:
            err = ols_error(e)
        if err < min:
            min = err
            placed_edge = e
    tailn = placed_edge.tail_node
    headn = placed_edge.head_node
    tailn.remove_child(headn)
    nn = dy.Node()
    nn.add_child(headn)
    qry = dy.Node(taxon=dy.Taxon(name))
    nn.add_child(qry)
    qry.edge_length = placed_edge.x_1
    tailn.add_child(nn)
    if placed_edge.head_node in list(master_edge.head_node.ancestor_iter()) or master_edge == placed_edge:
        nn.edge_length = placed_edge.length - max(placed_edge.x_2,0)
        headn.edge_length = max(placed_edge.x_2,0)
    else:
        nn.edge_length = max(placed_edge.x_2,0)
        headn.edge_length = placed_edge.length - max(placed_edge.x_2,0)
    #print(placed_edge.x_1, placed_edge.x_2)


if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-t", "--tree", dest="tree_fp",
                      help="path to the reference tree", metavar="FILE")
    parser.add_option("-d", "--distances", dest="dist_fp",
                      help="path to the table of observed distances", metavar="FILE")
    parser.add_option("-q", "--query", dest="query_name",
                      help="name of the query taxon", metavar="NAME")
    parser.add_option("-a", "--algo", dest="algo_name",
                      help="name of the algorithm (OLS, or FM)", metavar="ALGO")
    parser.add_option("-m", "--me", action="store_true", dest="minimum_evolution", default=False,
                      help="performs the placement based on minimum evolution principle")

    (options, args) = parser.parse_args()
    tree_fp = options.tree_fp
    dist_fp = options.dist_fp
    algo_name = options.algo_name
    query_name = options.query_name
    minimum_evolution = options.minimum_evolution

    df = pd.read_csv(dist_fp, sep="\s+", names = ["taxa", "dist"] , dtype = {"taxa": object}, header = None)
    obs_dist = pd.Series(df.dist.values,index=df.taxa).to_dict()
    f = open(tree_fp)
    tree = dy.Tree.get_from_string(f.readline(), schema='newick')
    master_edge = next(tree.postorder_edge_iter())
    dfs_S_values(master_edge, master_edge.tail_node)
    dfs_R_values(master_edge, None, master_edge.head_node, master_edge.tail_node)
    if algo_name == "OLS":
        ols_parameters_per_edge(tree)
    else:
        fm_parameters_per_edge(tree)
    placement(tree, query_name, algo_name)
    tree.write(file = sys.stdout, schema = "newick")
