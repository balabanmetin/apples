from apples.Algorithm import Algorithm
from apples import util


# Glossary
# OLS: Ordinary Least Squares
# FM: Fitch-Margoliash --Least Squares with \delta^2 scaling
# BE: Beyer --Least Squares with \delta scaling

class OLS(Algorithm):

    def all_S_values(self):
        subtree = self.subtree
        obs_dist = subtree.obs_dist
        for node in filter(lambda x: x.valid, subtree.traverse_postorder()):
            if node.is_leaf():
                node.S = 1
                node.Sd = 0
                node.Sd2 = 0
                node.SDd = 0
                node.SD = obs_dist[node.label]
                node.SD2 = node.SD * node.SD

            else:
                node.SDd, node.SD2, node.SD, node.S, node.Sd, node.Sd2 = 6 * [0]
                for child in filter(lambda x: x.valid, node.children):
                    node.S += child.S
                    node.Sd += child.S * child.edge_length + child.Sd
                    node.Sd2 += child.S * child.edge_length * child.edge_length + child.Sd2 +\
                                2 * child.edge_length * child.Sd
                    node.SDd += child.edge_length * child.SD + child.SDd
                    node.SD2 += child.SD2
                    node.SD += child.SD

    def all_R_values(self):
        subtree = self.subtree
        for node in filter(lambda x: x.valid, subtree.traverse_preorder()):
            node.RDd, node.RD2, node.RD, node.R, node.Rd, node.Rd2 = 6 * [0]
            for sibling in filter(lambda x: x.valid and x != node, node.parent.children):
                node.R += sibling.S
                node.Rd += sibling.S * sibling.edge_length + sibling.Sd
                node.Rd2 += sibling.S * sibling.edge_length * sibling.edge_length + sibling.Sd2 + \
                            2 * sibling.edge_length * sibling.Sd
                node.RDd += sibling.SD * sibling.edge_length + sibling.SDd
                node.RD2 += sibling.SD2
                node.RD += sibling.SD
            if node.parent != subtree.root and node.parent.valid:
                node.R += node.parent.R
                node.Rd += node.parent.R * node.parent.edge_length + node.parent.Rd
                node.Rd2 += node.parent.R * node.parent.edge_length * node.parent.edge_length + node.parent.Rd2 + \
                            2 * node.parent.edge_length * node.parent.Rd
                node.RDd += node.parent.RD * node.parent.edge_length + node.parent.RDd
                node.RD2 += node.parent.RD2
                node.RD += node.parent.RD

    # computes all alternative OLS placements
    def placement_per_edge(self, negative_branch):
        for node in filter(lambda x: x.valid, self.subtree.traverse_postorder()):
            a_11 = node.R + node.S
            a_12 = node.R - node.S
            a_21 = a_12
            a_22 = a_11
            c_1 = node.RD + node.SD - node.edge_length * node.S - node.Rd - node.Sd
            c_2 = node.RD - node.SD + node.edge_length * node.S - node.Rd + node.Sd
            util.solve2_2(node, a_11, a_12, a_21, a_22, c_1, c_2, negative_branch)

    # computes OLS error (Q value) for a given edge
    @staticmethod
    def error_per_edge(node):
        assert node.valid
        A = node.RD2 + node.SD2
        B = 2 * (node.x_1 + node.x_2) * node.Rd + \
            2 * (node.edge_length + node.x_1 - node.x_2) * node.Sd
        C = (node.x_1 + node.x_2) ** 2 * node.R + \
            (node.edge_length + node.x_1 - node.x_2) ** 2 * node.S
        D = -2 * (node.x_1 + node.x_2) * node.RD - \
            2 * (node.edge_length + node.x_1 - node.x_2) * node.SD
        E = -2 * node.RDd - 2 * node.SDd
        F = node.Rd2 + node.Sd2
        return A + B + C + D + E + F