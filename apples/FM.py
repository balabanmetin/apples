from apples.Algorithm import Algorithm
from apples import util


class FM(Algorithm):

    def all_S_values(self):
        subtree = self.subtree
        obs_dist = subtree.obs_dist
        for node in filter(lambda x: x.valid, subtree.traverse_postorder()):
            if node.is_leaf():
                node.S = 1
                node.Sd_D = 0
                node.Sd_D2 = 0
                node.Sd2_D2 = 0
                node.S1_D = 1.0 / obs_dist[node.label]
                node.S1_D2 = 1.0 / (obs_dist[node.label] * obs_dist[node.label])

            else:
                node.S, node.Sd_D, node.Sd_D2, node.Sd2_D2, node.S1_D, node.S1_D2 = 6 * [0]
                for child in filter(lambda x: x.valid, node.children):
                    node.S += child.S
                    node.Sd_D += child.edge_length * child.S1_D + child.Sd_D
                    node.Sd_D2 += child.edge_length * child.S1_D2 + child.Sd_D2
                    node.Sd2_D2 += child.S1_D2 * child.edge_length * child.edge_length + child.Sd2_D2 + 2 * child.edge_length * child.Sd_D2
                    node.S1_D += child.S1_D
                    node.S1_D2 += child.S1_D2

    def all_R_values(self):
        subtree = self.subtree
        for node in filter(lambda x: x.valid, subtree.traverse_preorder()):
            node.R, node.Rd_D, node.Rd_D2, node.Rd2_D2, node.R1_D, node.R1_D2 = 6 * [0]
            for sibling in filter(lambda x: x.valid and x != node, node.parent.children):
                node.R += sibling.S
                node.Rd_D += sibling.edge_length * sibling.S1_D + sibling.Sd_D
                node.Rd_D2 += sibling.edge_length * sibling.S1_D2 + sibling.Sd_D2
                node.Rd2_D2 += sibling.S1_D2 * sibling.edge_length * sibling.edge_length + sibling.Sd2_D2 + 2 * sibling.edge_length * sibling.Sd_D2
                node.R1_D += sibling.S1_D
                node.R1_D2 += sibling.S1_D2
            if node.parent != subtree.root and node.parent.valid:
                node.R += node.parent.R
                node.Rd_D += node.parent.edge_length * node.parent.R1_D + node.parent.Rd_D
                node.Rd_D2 += node.parent.edge_length * node.parent.R1_D2 + node.parent.Rd_D2
                node.Rd2_D2 += node.parent.R1_D2 * node.parent.edge_length * node.parent.edge_length + node.parent.Rd2_D2 + 2 * node.parent.edge_length * node.parent.Rd_D2
                node.R1_D += node.parent.R1_D
                node.R1_D2 += node.parent.R1_D2


    # computes all alternative Fitch-Margoliash placements
    def placement_per_edge(self, negative_branch):
        for node in filter(lambda x: x.valid, self.subtree.traverse_postorder()):
            a_11 = node.R1_D2 + node.S1_D2
            a_12 = node.R1_D2 - node.S1_D2
            a_21 = a_12
            a_22 = a_11
            c_1 = node.R1_D + node.S1_D - node.edge_length * node.S1_D2 - node.Rd_D2 - node.Sd_D2
            c_2 = node.R1_D - node.S1_D + node.edge_length * node.S1_D2 - node.Rd_D2 + node.Sd_D2
            util.solve2_2(node, a_11, a_12, a_21, a_22, c_1, c_2, negative_branch)

    # computes Fitch-Margoliash error (Q value) for a given edge
    @staticmethod
    def error_per_edge(node):
        assert node.valid
        A = node.R + node.S
        B = 2 * (node.x_1 + node.x_2) * node.Rd_D2 + \
            2 * (node.edge_length + node.x_1 - node.x_2) * node.Sd_D2
        C = (node.x_1 + node.x_2) ** 2 * node.R1_D2 + \
            (node.edge_length + node.x_1 - node.x_2) ** 2 * node.S1_D2
        D = -2 * (node.x_1 + node.x_2) * node.R1_D - \
            2 * (node.edge_length + node.x_1 - node.x_2) * node.S1_D
        E = -2 * node.Rd_D - 2 * node.Sd_D
        F = node.Rd2_D2 + node.Sd2_D2
        return A + B + C + D + E + F