import dendropy as dy
import sys
import pandas as pd


def dfs_S_values(edge, downstream):
    if downstream.is_leaf():
        edge.S = 1
        edge.Sd = 0
        edge.Sd2 = 0
        edge.SD = obs_dist[downstream.taxon.label]
        edge.SD2 = edge.SD * edge.SD
        edge.SDd = 0


    else:
        inc = list(downstream.incident_edges())
        inc = filter(lambda x: x.length != None and x != edge, inc)
        d1, d2 = inc
        d1tips = [d1.head_node, d1.tail_node]
        d1tips.remove(downstream)
        d2tips = [d2.head_node, d2.tail_node]
        d2tips.remove(downstream)
        if len(d1tips) == 0 or len(d2tips) == 0:
            print("error")
        dfs_S_values(d1, d1tips[0])
        dfs_S_values(d2, d2tips[0])
        edge.S = d1.S + d2.S
        edge.Sd = d1.S * d1.length + d2.S * d2.length + d1.Sd + d2.Sd
        edge.Sd2 = d1.S * d1.length * d1.length + d2.S * d2.length * d2.length + d1.Sd2 + d2.Sd2 \
                   + 2 * d1.length * d1.Sd + 2 * d2.length * d2.Sd
        edge.SD = d1.SD + d2.SD
        edge.SD2 = d1.SD2 + d2.SD2
        edge.SDd = d1.length * d1.SD + d2.length * d2.SD + d1.SDd + d2.SDd


def dfs_R_values(edge, u1, upstream, downstream):
    if upstream.is_leaf():
        edge.R = 1
        edge.Rd = 0
        edge.Rd2 = 0
        edge.RD = obs_dist[upstream.taxon.label]
        edge.RD2 = edge.RD * edge.RD
        edge.RDd = 0

    else:
        inc = list(upstream.incident_edges())
        u2 = list(filter(lambda x: x.length != None and x != edge and x != u1, inc))
        assert len(u2) is 1
        edge.R = u1.R + u2[0].S
        edge.Rd = u2[0].S * u2[0].length + u1.R * u1.length + u1.Rd + u2[0].Sd
        edge.Rd2 = u2[0].S * u2[0].length * u2[0].length + u1.R * u1.length * u1.length + u1.Rd2 + u2[0].Sd2 \
                   + 2 * u2[0].length * u2[0].Sd + 2 * u1.length * u1.Rd
        edge.RD = u1.RD + u2[0].SD
        edge.RD2 = u1.RD2 + u2[0].SD2
        edge.RDd = u2[0].SD * u2[0].length + u1.RD * u1.length + u1.RDd + u2[0].SDd

    if not downstream.is_leaf():
        inc = list(downstream.incident_edges())
        inc = filter(lambda x: x.length != None and x != edge, inc)
        d1, d2 = inc
        d1tips = [d1.head_node, d1.tail_node]
        d1tips.remove(downstream)
        d2tips = [d2.head_node, d2.tail_node]
        d2tips.remove(downstream)
        dfs_R_values(d1, edge, downstream, d1tips[0])
        dfs_R_values(d2, edge, downstream, d2tips[0])


def solve2_2(a_11, a_12, a_21, a_22, c_1, c_2):
    det = 1 / (a_11 * a_22 - a_12 * a_21)
    assert det is not 0
    x_1 = (a_22 * c_1 - a_12 * c_2) * det
    x_2 = (- a_21 * c_1 + a_11 * c_2) * det
    return (x_1, x_2)


def ols_parameters_per_edge(tree):
    for e in tree.postorder_edge_iter():
        if e.head_node == tree.seed_node:
            continue
        a_11 = n
        a_12 = -n + 2 * e.R
        a_21 = a_12
        a_22 = a_11
        c_1 = Dsum - e.length * (n - e.R) - e.Rd - e.Sd
        c_2 = e.RD - e.SD + e.length * (n - e.R) - e.Rd + e.Sd
        x_1, x_2 = solve2_2(a_11, a_12, a_21, a_22, c_1, c_2)
        e.x_1 = x_1
        e.x_2 = x_2


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


def placement(tree, name, minimum_evolution):
    min = 999999999999
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
    print(placed_edge.x_1, placed_edge.x_2)


if __name__ == "__main__":
    fname = sys.argv[1]
    df = pd.read_csv(sys.argv[2], sep="\s+", names = ["taxa", "dist"] , header = None)
    obs_dist = pd.Series(df.dist.values,index=df.taxa).to_dict()
    f = open(fname)
    tree = dy.Tree.get_from_string(f.readline(), schema='newick')
    Dsum = sum(obs_dist.values())
    n = len(obs_dist)
    master_edge = next(tree.postorder_edge_iter())
    dfs_S_values(master_edge, master_edge.tail_node)
    dfs_R_values(master_edge, None, master_edge.head_node, master_edge.tail_node)
    ols_parameters_per_edge(tree)
    placement(tree, sys.argv[2].split("/")[-1].split(".")[0], 1)
    tree.write(file = sys.stdout, schema = "newick")
