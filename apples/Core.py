class Core:
    def __init__(self, tree):
        self.tree = tree

    def init(self):
        self.tree_S_values()
        self.tree_R_values()

    def dp(self, obs_dist):
        self.observed_S_values(obs_dist)
        self.observed_R_values()

    def dp_frag(self, subtree):
        self.all_S_values(subtree)
        self.all_R_values(subtree)

    def observed_S_values(self, obs_dist):
        for node in self.tree.traverse_postorder():
            if node.is_leaf():
                node.SDd = 0
                node.Sd_D = 0
                node.Sd_D2 = 0
                node.Sd2_D = 0
                node.Sd2_D2 = 0
                node.SD = obs_dist[node.label]
                node.SD2 = node.SD * node.SD
                node.S1_D = 1.0 / node.SD
                node.S1_D2 = 1.0 / (node.SD * node.SD)

            else:
                node.SDd, node.Sd_D, node.Sd_D2, node.Sd2_D, node.Sd2_D2, node.SD2, node.SD, node.S1_D, node.S1_D2 = 9 * [
                    0]
                for child in node.children:
                    node.SDd += child.edge_length * child.SD + child.SDd
                    node.Sd_D += child.edge_length * child.S1_D + child.Sd_D
                    node.Sd_D2 += child.edge_length * child.S1_D2 + child.Sd_D2
                    node.Sd2_D += child.S1_D * child.edge_length * child.edge_length + child.Sd2_D + 2 * child.edge_length * child.Sd_D
                    node.Sd2_D2 += child.S1_D2 * child.edge_length * child.edge_length + child.Sd2_D2 + 2 * child.edge_length * child.Sd_D2
                    node.SD2 += child.SD2
                    node.SD += child.SD
                    node.S1_D += child.S1_D
                    node.S1_D2 += child.S1_D2

    def observed_R_values(self):
        for node in self.tree.traverse_preorder():
            if node == self.tree.root:
                continue
            node.RDd, node.Rd_D, node.Rd_D2, node.Rd2_D, node.Rd2_D2, node.RD2, node.RD, node.R1_D, node.R1_D2 = 9 * [0]
            for sibling in node.parent.children:
                if sibling != node:
                    node.RDd += sibling.SD * sibling.edge_length + sibling.SDd
                    node.Rd_D += sibling.edge_length * sibling.S1_D + sibling.Sd_D
                    node.Rd_D2 += sibling.edge_length * sibling.S1_D2 + sibling.Sd_D2
                    node.Rd2_D += sibling.S1_D * sibling.edge_length * sibling.edge_length + sibling.Sd2_D + 2 * sibling.edge_length * sibling.Sd_D
                    node.Rd2_D2 += sibling.S1_D2 * sibling.edge_length * sibling.edge_length + sibling.Sd2_D2 + 2 * sibling.edge_length * sibling.Sd_D2
                    node.RD2 += sibling.SD2
                    node.RD += sibling.SD
                    node.R1_D += sibling.S1_D
                    node.R1_D2 += sibling.S1_D2
            if node.parent != self.tree.root:
                node.RDd += node.parent.RD * node.parent.edge_length + node.parent.RDd
                node.Rd_D += node.parent.edge_length * node.parent.R1_D + node.parent.Rd_D
                node.Rd_D2 += node.parent.edge_length * node.parent.R1_D2 + node.parent.Rd_D2
                node.Rd2_D += node.parent.R1_D * node.parent.edge_length * node.parent.edge_length + node.parent.Rd2_D + 2 * node.parent.edge_length * node.parent.Rd_D
                node.Rd2_D2 += node.parent.R1_D2 * node.parent.edge_length * node.parent.edge_length + node.parent.Rd2_D2 + 2 * node.parent.edge_length * node.parent.Rd_D2
                node.RD2 += node.parent.RD2
                node.RD += node.parent.RD
                node.R1_D += node.parent.R1_D
                node.R1_D2 += node.parent.R1_D2

    def tree_S_values(self):
        for node in self.tree.traverse_postorder():
            if node.is_leaf():
                node.S = 1
                node.Sd = 0
                node.Sd2 = 0

            else:
                node.S, node.Sd, node.Sd2 = 3 * [0]
                for child in node.children:
                    node.S += child.S
                    node.Sd += child.S * child.edge_length + child.Sd
                    node.Sd2 += child.S * child.edge_length * child.edge_length + child.Sd2 + 2 * child.edge_length * child.Sd

    def tree_R_values(self):
        for node in self.tree.traverse_preorder():
            if node == self.tree.root:
                continue
            node.R, node.Rd, node.Rd2 = 3 * [0]
            for sibling in node.parent.children:
                if sibling != node:
                    node.R += sibling.S
                    node.Rd += sibling.S * sibling.edge_length + sibling.Sd
                    node.Rd2 += sibling.S * sibling.edge_length * sibling.edge_length + sibling.Sd2 + \
                                2 * sibling.edge_length * sibling.Sd
            if node.parent != self.tree.root:
                node.R += node.parent.R
                node.Rd += node.parent.R * node.parent.edge_length + node.parent.Rd
                node.Rd2 += node.parent.R * node.parent.edge_length * node.parent.edge_length + node.parent.Rd2 + \
                            2 * node.parent.edge_length * node.parent.Rd

    def all_S_values(self, subtree):
        obs_dist = subtree.obs_dist
        for node in filter(lambda x: x.valid, subtree.traverse_postorder()):
            if node.is_leaf():
                node.BS = 1
                node.BSd = 0
                node.BSd2 = 0
                node.BSDd = 0
                node.BSD = obs_dist[node.label]
                node.BSD2 = node.BSD * node.BSD

                node.S = 1
                node.Sd = 0
                node.Sd2 = 0
                node.SDd = 0
                node.Sd_D = 0
                node.Sd_D2 = 0
                node.Sd2_D = 0
                node.Sd2_D2 = 0
                node.SD = obs_dist[node.label]
                node.SD2 = node.SD * node.SD
                node.S1_D = 1.0 / node.SD
                node.S1_D2 = 1.0 / (node.SD * node.SD)

            else:
                node.SDd, node.Sd_D, node.Sd_D2, node.Sd2_D, node.Sd2_D2, node.SD2, node.SD, node.S1_D, node.S1_D2 = 9 * [
                    0]
                node.S, node.Sd, node.Sd2 = 3 * [0]
                node.BS, node.BSd, node.BSd2 = 3 * [0]
                node.BSDd, node.BSD2, node.BSD = 3 * [0]
                bmecoef = 1 / len(list(filter(lambda x: x.valid, node.children)))
                for child in filter(lambda x: x.valid, node.children):
                    node.BS += bmecoef * child.BS
                    node.BSd += bmecoef * (child.BS * child.edge_length + child.BSd)
                    node.BSd2 += bmecoef * (
                            child.BS * child.edge_length * child.edge_length + child.BSd2 + 2 * child.edge_length * child.BSd)
                    node.BSDd += bmecoef * (child.edge_length * child.BSD + child.BSDd)
                    node.BSD2 += bmecoef * (child.BSD2)
                    node.BSD += bmecoef * (child.BSD)

                    node.S += child.S
                    node.Sd += child.S * child.edge_length + child.Sd
                    node.Sd2 += child.S * child.edge_length * child.edge_length + child.Sd2 + 2 * child.edge_length * child.Sd
                    node.SDd += child.edge_length * child.SD + child.SDd
                    node.Sd_D += child.edge_length * child.S1_D + child.Sd_D
                    node.Sd_D2 += child.edge_length * child.S1_D2 + child.Sd_D2
                    node.Sd2_D += child.S1_D * child.edge_length * child.edge_length + child.Sd2_D + 2 * child.edge_length * child.Sd_D
                    node.Sd2_D2 += child.S1_D2 * child.edge_length * child.edge_length + child.Sd2_D2 + 2 * child.edge_length * child.Sd_D2
                    node.SD2 += child.SD2
                    node.SD += child.SD
                    node.S1_D += child.S1_D
                    node.S1_D2 += child.S1_D2

    def all_R_values(self, subtree):
        for node in filter(lambda x: x.valid, subtree.traverse_preorder()):
            node.RDd, node.Rd_D, node.Rd_D2, node.Rd2_D, node.Rd2_D2, node.RD2, node.RD, node.R1_D, node.R1_D2 = 9 * [0]
            node.R, node.Rd, node.Rd2 = 3 * [0]
            node.BR, node.BRd, node.BRd2 = 3 * [0]
            node.BRDd, node.BRD2, node.BRD = 3 * [0]
            nonrootnode = 1 if node.parent != self.tree.root else 0
            bmecoef = 1 / (nonrootnode + len(list(filter(lambda x: x.valid and x != node, node.parent.children))))
            for sibling in filter(lambda x: x.valid and x != node, node.parent.children):
                node.BR += bmecoef * sibling.BS
                node.BRd += bmecoef * (sibling.BS * sibling.edge_length + sibling.BSd)
                node.BRd2 += bmecoef * (sibling.BS * sibling.edge_length * sibling.edge_length + sibling.BSd2 +
                                        2 * sibling.edge_length * sibling.BSd)
                node.BRDd += bmecoef * (sibling.BSD * sibling.edge_length + sibling.BSDd)
                node.BRD2 += bmecoef * sibling.BSD2
                node.BRD += bmecoef * sibling.BSD

                node.R += sibling.S
                node.Rd += sibling.S * sibling.edge_length + sibling.Sd
                node.Rd2 += sibling.S * sibling.edge_length * sibling.edge_length + sibling.Sd2 + \
                            2 * sibling.edge_length * sibling.Sd
                node.RDd += sibling.SD * sibling.edge_length + sibling.SDd
                node.Rd_D += sibling.edge_length * sibling.S1_D + sibling.Sd_D
                node.Rd_D2 += sibling.edge_length * sibling.S1_D2 + sibling.Sd_D2
                node.Rd2_D += sibling.S1_D * sibling.edge_length * sibling.edge_length + sibling.Sd2_D + 2 * sibling.edge_length * sibling.Sd_D
                node.Rd2_D2 += sibling.S1_D2 * sibling.edge_length * sibling.edge_length + sibling.Sd2_D2 + 2 * sibling.edge_length * sibling.Sd_D2
                node.RD2 += sibling.SD2
                node.RD += sibling.SD
                node.R1_D += sibling.S1_D
                node.R1_D2 += sibling.S1_D2
            if node.parent != subtree.root and node.parent.valid:
                node.BR += bmecoef * node.parent.BR
                node.BRd += bmecoef * (node.parent.BR * node.parent.edge_length + node.parent.BRd)
                node.BRd2 += bmecoef * (
                            node.parent.BR * node.parent.edge_length * node.parent.edge_length + node.parent.BRd2 +
                            2 * node.parent.edge_length * node.parent.BRd)
                node.BRDd += bmecoef * (node.parent.BRD * node.parent.edge_length + node.parent.BRDd)
                node.BRD2 += bmecoef * node.parent.BRD2
                node.BRD += bmecoef * node.parent.BRD

                node.R += node.parent.R
                node.Rd += node.parent.R * node.parent.edge_length + node.parent.Rd
                node.Rd2 += node.parent.R * node.parent.edge_length * node.parent.edge_length + node.parent.Rd2 + \
                            2 * node.parent.edge_length * node.parent.Rd
                node.RDd += node.parent.RD * node.parent.edge_length + node.parent.RDd
                node.Rd_D += node.parent.edge_length * node.parent.R1_D + node.parent.Rd_D
                node.Rd_D2 += node.parent.edge_length * node.parent.R1_D2 + node.parent.Rd_D2
                node.Rd2_D += node.parent.R1_D * node.parent.edge_length * node.parent.edge_length + node.parent.Rd2_D + 2 * node.parent.edge_length * node.parent.Rd_D
                node.Rd2_D2 += node.parent.R1_D2 * node.parent.edge_length * node.parent.edge_length + node.parent.Rd2_D2 + 2 * node.parent.edge_length * node.parent.Rd_D2
                node.RD2 += node.parent.RD2
                node.RD += node.parent.RD
                node.R1_D += node.parent.R1_D
                node.R1_D2 += node.parent.R1_D2
