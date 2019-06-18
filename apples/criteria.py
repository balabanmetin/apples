from apples.Algorithm import Algorithm
from apples import util


# Glossary
# OLS: Ordinary Least Squares
# FM: Fitch-Margoliash --Least Squares with \delta^2 scaling
# BE: Beyer --Least Squares with \delta scaling



class OLS(Algorithm):

    # computes all alternative OLS placements
    def placement_per_edge(self, negative_branch):
        for node in filter(lambda x: x.valid, self.tree.traverse_postorder()):
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


class FM(Algorithm):

    # computes all alternative Fitch-Margoliash placements
    def placement_per_edge(self, negative_branch):
        for node in filter(lambda x: x.valid, self.tree.traverse_postorder()):
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


class BE(Algorithm):

    # computes all alternative BE placements
    def placement_per_edge(self, negative_branch):
        for node in filter(lambda x: x.valid, self.tree.traverse_postorder()):
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
