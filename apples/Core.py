
class Core:
    def __init__(self, tree):
        self.tree = tree
        tree.master_edge = next(tree.postorder_edge_iter())
        master_edge = self.tree.master_edge
        self.tree_S_values(master_edge, master_edge.tail_node)
        self.tree_R_values(master_edge, None, master_edge.head_node, master_edge.tail_node)


    def dp(self, obs_dist):
        master_edge = self.tree.master_edge
        self.observed_S_values(master_edge, master_edge.tail_node, obs_dist)
        self.observed_R_values(master_edge, None, master_edge.head_node, master_edge.tail_node,obs_dist)

    def observed_S_values(self, edge, downstream, obs_dist):
        if downstream.is_leaf():
            edge.SDd = 0
            edge.Sd_D = 0
            edge.Sd_D2 = 0
            edge.Sd2_D = 0
            edge.Sd2_D2 = 0
            edge.SD = obs_dist[downstream.taxon.label]
            edge.SD2 = edge.SD * edge.SD
            edge.S1_D = 1 / edge.SD
            edge.S1_D2 = 1 / (edge.SD * edge.SD)


        else:
            inc = list(downstream.incident_edges())
            inc = filter(lambda e: e.head_node != self.tree.seed_node and e != edge, inc)
            edge.SDd, edge.Sd_D, edge.Sd_D2, edge.Sd2_D, edge.Sd2_D2, edge.SD2, edge.SD, edge.S1_D, edge.S1_D2 = 9 * [0]
            for d1 in inc:
                d1tips = [d1.head_node, d1.tail_node]
                d1tips.remove(downstream)
                if len(d1tips) == 0:
                    print("error")
                self.observed_S_values(d1, d1tips[0], obs_dist)
                edge.SDd += d1.length * d1.SD + d1.SDd
                edge.Sd_D += d1.length * d1.S1_D + d1.Sd_D
                edge.Sd_D2 += d1.length * d1.S1_D2 + d1.Sd_D2
                edge.Sd2_D += d1.S1_D * d1.length * d1.length + d1.Sd2_D + 2 * d1.length * d1.Sd_D
                edge.Sd2_D2 += d1.S1_D2 * d1.length * d1.length + d1.Sd2_D2 + 2 * d1.length * d1.Sd_D2
                edge.SD2 += d1.SD2
                edge.SD += d1.SD
                edge.S1_D += d1.S1_D
                edge.S1_D2 += d1.S1_D2

    def observed_R_values(self, edge, u1, upstream, downstream, obs_dist):
        if upstream.is_leaf():
            edge.RDd = 0
            edge.Rd_D = 0
            edge.Rd_D2 = 0
            edge.Rd2_D = 0
            edge.Rd2_D2 = 0
            edge.RD = obs_dist[upstream.taxon.label]
            edge.RD2 = edge.RD * edge.RD
            edge.R1_D = 1 / edge.RD
            edge.R1_D2 = 1 / (edge.RD * edge.RD)

        else:
            edge.RDd = u1.RD * u1.length + u1.RDd
            edge.Rd_D = u1.length * u1.R1_D + u1.Rd_D
            edge.Rd_D2 = u1.length * u1.R1_D2 + u1.Rd_D2
            edge.Rd2_D = u1.R1_D * u1.length * u1.length + u1.Rd2_D + 2 * u1.length * u1.Rd_D
            edge.Rd2_D2 = u1.R1_D2 * u1.length * u1.length + u1.Rd2_D2 + 2 * u1.length * u1.Rd_D2
            edge.RD2 = u1.RD2
            edge.RD = u1.RD
            edge.R1_D = u1.R1_D
            edge.R1_D2 = u1.R1_D2
            inc = list(upstream.incident_edges())
            u2 = list(filter(lambda e: e.head_node != self.tree.seed_node and e != u1 and e != edge, inc))
            if len(u2) == 1:
                edge.RDd += u2[0].SD * u2[0].length + u2[0].SDd
                edge.Rd_D += u2[0].length * u2[0].S1_D + u2[0].Sd_D
                edge.Rd_D2 += u2[0].length * u2[0].S1_D2 + u2[0].Sd_D2
                edge.Rd2_D += u2[0].S1_D * u2[0].length * u2[0].length + u2[0].Sd2_D + 2 * u2[0].length * u2[0].Sd_D
                edge.Rd2_D2 += u2[0].S1_D2 * u2[0].length * u2[0].length + u2[0].Sd2_D2 + 2 * u2[0].length * u2[0].Sd_D2
                edge.RD2 += u2[0].SD2
                edge.RD += u2[0].SD
                edge.R1_D += u2[0].S1_D
                edge.R1_D2 += u2[0].S1_D2

        if not downstream.is_leaf():
            inc = list(downstream.incident_edges())
            inc = filter(lambda e: e.head_node != self.tree.seed_node and e != edge, inc)
            for d1 in inc:
                d1tips = [d1.head_node, d1.tail_node]
                d1tips.remove(downstream)
                self.observed_R_values(d1, edge, downstream, d1tips[0], obs_dist)

    def tree_S_values(self, edge, downstream):
        if downstream.is_leaf():
            edge.S = 1
            edge.Sd = 0
            edge.Sd2 = 0

        else:
            inc = list(downstream.incident_edges())
            inc = filter(lambda e: e.head_node != self.tree.seed_node and e != edge, inc)
            edge.S, edge.Sd, edge.Sd2 = 3 * [0]
            for d1 in inc:
                d1tips = [d1.head_node, d1.tail_node]
                d1tips.remove(downstream)
                if len(d1tips) == 0:
                    print("error")
                self.tree_S_values(d1, d1tips[0])
                edge.S += d1.S
                edge.Sd += d1.S * d1.length + d1.Sd
                edge.Sd2 += d1.S * d1.length * d1.length + d1.Sd2 + 2 * d1.length * d1.Sd

    def tree_R_values(self, edge, u1, upstream, downstream):
        if upstream.is_leaf():
            edge.R = 1
            edge.Rd = 0
            edge.Rd2 = 0


        else:
            edge.R = u1.R
            edge.Rd = u1.R * u1.length + u1.Rd
            edge.Rd2 = u1.R * u1.length * u1.length + u1.Rd2 + 2 * u1.length * u1.Rd
            inc = list(upstream.incident_edges())
            u2 = list(filter(lambda e: e.head_node != self.tree.seed_node and e != u1 and e != edge, inc))
            if len(u2) == 1:
                edge.R += u2[0].S
                edge.Rd += u2[0].S * u2[0].length + u2[0].Sd
                edge.Rd2 += u2[0].S * u2[0].length * u2[0].length + u2[0].Sd2 + 2 * u2[0].length * u2[0].Sd

        if not downstream.is_leaf():
            inc = list(downstream.incident_edges())
            inc = filter(lambda e: e.head_node != self.tree.seed_node and e != edge, inc)
            for d1 in inc:
                d1tips = [d1.head_node, d1.tail_node]
                d1tips.remove(downstream)
                self.tree_R_values(d1, edge, downstream, d1tips[0])