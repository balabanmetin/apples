from apples.Algorithm import Algorithm
from apples import util


class BE(Algorithm):

    def all_S_values(self):
        subtree = self.subtree
        obs_dist = subtree.obs_dist
        for node in filter(lambda x: x.valid, subtree.traverse_postorder()):
            if node.is_leaf():
                node.S = 1
                node.Sd = 0
                node.Sd_D = 0
                node.Sd2_D = 0
                node.SD = obs_dist[node.label]
                node.S1_D = 1.0 / node.SD

            else:
                node.S, node.Sd, node.Sd_D, node.Sd2_D, node.SD, node.S1_D = 6 * [0]
                for child in filter(lambda x: x.valid, node.children):
                    node.S += child.S
                    node.Sd += child.S * child.edge_length + child.Sd
                    node.Sd_D += child.edge_length * child.S1_D + child.Sd_D
                    node.Sd2_D += child.S1_D * child.edge_length * child.edge_length + child.Sd2_D + 2 * child.edge_length * child.Sd_D
                    node.SD += child.SD
                    node.S1_D += child.S1_D

    def all_R_values(self):
        subtree = self.subtree
        for node in filter(lambda x: x.valid, subtree.traverse_preorder()):
            node.R, node.Rd, node.Rd_D, node.Rd2_D, node.RD, node.R1_D = 6 * [0]
            for sibling in filter(lambda x: x.valid and x != node, node.parent.children):
                node.R += sibling.S
                node.Rd += sibling.S * sibling.edge_length + sibling.Sd
                node.Rd_D += sibling.edge_length * sibling.S1_D + sibling.Sd_D
                node.Rd2_D += sibling.S1_D * sibling.edge_length * sibling.edge_length + sibling.Sd2_D + 2 * sibling.edge_length * sibling.Sd_D
                node.RD += sibling.SD
                node.R1_D += sibling.S1_D
            if node.parent != subtree.root and node.parent.valid:
                node.R += node.parent.R
                node.Rd += node.parent.R * node.parent.edge_length + node.parent.Rd
                node.Rd_D += node.parent.edge_length * node.parent.R1_D + node.parent.Rd_D
                node.Rd2_D += node.parent.R1_D * node.parent.edge_length * node.parent.edge_length + node.parent.Rd2_D + 2 * node.parent.edge_length * node.parent.Rd_D
                node.RD += node.parent.RD
                node.R1_D += node.parent.R1_D

    # computes all alternative BE placements
    def placement_per_edge(self, negative_branch):
        for node in filter(lambda x: x.valid, self.subtree.traverse_postorder()):
            a_11 = node.R1_D + node.S1_D
            a_12 = node.R1_D - node.S1_D
            a_21 = a_12
            a_22 = a_11
            c_1 = node.R + node.S - node.edge_length * node.S1_D - node.Rd_D - node.Sd_D
            c_2 = node.R - node.S + node.edge_length * node.S1_D - node.Rd_D + node.Sd_D
            util.solve2_2(node, a_11, a_12, a_21, a_22, c_1, c_2, negative_branch)

    # computes BE error (Q value) for a given edge
    @staticmethod
    def error_per_edge(node):
        assert node.valid
        A = node.RD + node.SD
        B = 2 * (node.x_1 + node.x_2) * node.Rd_D + \
            2 * (node.edge_length + node.x_1 - node.x_2) * node.Sd_D
        C = (node.x_1 + node.x_2) ** 2 * node.R1_D + \
            (node.edge_length + node.x_1 - node.x_2) ** 2 * node.S1_D
        D = -2 * (node.x_1 + node.x_2) * node.R - \
            2 * (node.edge_length + node.x_1 - node.x_2) * node.S
        E = -2 * node.Rd - 2 * node.Sd
        F = node.Rd2_D + node.Sd2_D
        return A + B + C + D + E + F