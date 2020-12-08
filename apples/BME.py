from apples.Algorithm import Algorithm
from apples import util


class BME(Algorithm):

    def all_S_values(self):
        subtree = self.subtree
        obs_dist = subtree.obs_dist
        for node in filter(lambda x: x.valid, subtree.traverse_postorder()):
            if node.is_leaf():
                node.BS = 1
                node.BSd = 0
                node.BSd2 = 0
                node.BSDd = 0
                node.BSD = obs_dist[node.label]
                node.BSD2 = node.BSD * node.BSD
            else:
                node.BS, node.BSd, node.BSd2, node.BSDd, node.BSD2, node.BSD = 6 * [0]
                bmecoef = 1 / len(list(filter(lambda x: x.valid, node.children)))
                for child in filter(lambda x: x.valid, node.children):
                    node.BS += bmecoef * child.BS
                    node.BSd += bmecoef * (child.BS * child.edge_length + child.BSd)
                    node.BSd2 += bmecoef * (
                            child.BS * child.edge_length * child.edge_length + child.BSd2 + 2 * child.edge_length * child.BSd)
                    node.BSDd += bmecoef * (child.edge_length * child.BSD + child.BSDd)
                    node.BSD2 += bmecoef * child.BSD2
                    node.BSD += bmecoef * child.BSD

    def all_R_values(self):
        subtree = self.subtree
        for node in filter(lambda x: x.valid, subtree.traverse_preorder()):
            node.BR, node.BRd, node.BRd2, node.BRDd, node.BRD2, node.BRD = 6 * [0]
            nonrootnode = 1 if node.parent != subtree.root else 0
            bmecoef = 1 / (nonrootnode + len(list(filter(lambda x: x.valid and x != node, node.parent.children))))
            for sibling in filter(lambda x: x.valid and x != node, node.parent.children):
                node.BR += bmecoef * sibling.BS
                node.BRd += bmecoef * (sibling.BS * sibling.edge_length + sibling.BSd)
                node.BRd2 += bmecoef * (sibling.BS * sibling.edge_length * sibling.edge_length + sibling.BSd2 +
                                        2 * sibling.edge_length * sibling.BSd)
                node.BRDd += bmecoef * (sibling.BSD * sibling.edge_length + sibling.BSDd)
                node.BRD2 += bmecoef * sibling.BSD2
                node.BRD += bmecoef * sibling.BSD

            if node.parent != subtree.root and node.parent.valid:
                node.BR += bmecoef * node.parent.BR
                node.BRd += bmecoef * (node.parent.BR * node.parent.edge_length + node.parent.BRd)
                node.BRd2 += bmecoef * (
                            node.parent.BR * node.parent.edge_length * node.parent.edge_length + node.parent.BRd2 +
                            2 * node.parent.edge_length * node.parent.BRd)
                node.BRDd += bmecoef * (node.parent.BRD * node.parent.edge_length + node.parent.BRDd)
                node.BRD2 += bmecoef * node.parent.BRD2
                node.BRD += bmecoef * node.parent.BRD

    # computes all alternative Balanced Minimum evolution placements
    def placement_per_edge(self, negative_branch):
        for node in filter(lambda x: x.valid, self.subtree.traverse_postorder()):
            a_11 = node.BR + node.BS
            a_12 = node.BR - node.BS
            a_21 = a_12
            a_22 = a_11
            c_1 = node.BRD + node.BSD - node.edge_length * node.BS - node.BRd - node.BSd
            c_2 = node.BRD - node.BSD + node.edge_length * node.BS - node.BRd + node.BSd
            util.solve2_2(node, a_11, a_12, a_21, a_22, c_1, c_2, negative_branch)

    # computes BME error (Q value) for a given edge
    @staticmethod
    def error_per_edge(node):
        assert node.valid
        A = node.BRD2 + node.BSD2
        B = 2 * (node.x_1 + node.x_2) * node.BRd + \
            2 * (node.edge_length + node.x_1 - node.x_2) * node.BSd
        C = (node.x_1 + node.x_2) ** 2 * node.BR + \
            (node.edge_length + node.x_1 - node.x_2) ** 2 * node.BS
        D = -2 * (node.x_1 + node.x_2) * node.BRD - \
            2 * (node.edge_length + node.x_1 - node.x_2) * node.BSD
        E = -2 * node.BRDd - 2 * node.BSDd
        F = node.BRd2 + node.BSd2
        return A + B + C + D + E + F



