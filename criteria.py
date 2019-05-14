from Algorithm import Algorithm
import util


# Glossary
# OLS: Ordinary Least Squares
# FM: Fitch-Margoliash --Least Squares with \delta^2 scaling
# BE: Beyer --Least Squares with \delta scaling



class OLS(Algorithm):

    # computes all alternative OLS placements
    def placement_per_edge(self, negative_branch):
        for e in self.tree.postorder_edge_iter():
            if e.head_node == self.tree.seed_node:
                continue
            a_11 = e.R + e.S
            a_12 = e.R - e.S
            a_21 = a_12
            a_22 = a_11
            c_1 = e.RD + e.SD - e.length * e.S - e.Rd - e.Sd
            c_2 = e.RD - e.SD + e.length * e.S - e.Rd + e.Sd
            util.solve2_2(e, a_11, a_12, a_21, a_22, c_1, c_2, negative_branch)

    # computes OLS error (Q value) for a given edge
    @staticmethod
    def error_per_edge(edge):
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
    def placement_per_edge(self, negative_branch):
        for e in self.tree.postorder_edge_iter():
            if e.head_node == self.tree.seed_node:
                continue
            a_11 = e.R1_D2 + e.S1_D2
            a_12 = e.R1_D2 - e.S1_D2
            a_21 = a_12
            a_22 = a_11
            c_1 = e.R1_D + e.S1_D - e.length * e.S1_D2 - e.Rd_D2 - e.Sd_D2
            c_2 = e.R1_D - e.S1_D + e.length * e.S1_D2 - e.Rd_D2 + e.Sd_D2
            util.solve2_2(e, a_11, a_12, a_21, a_22, c_1, c_2, negative_branch)

    # computes Fitch-Margoliash error (Q value) for a given edge
    @staticmethod
    def error_per_edge(edge):
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
    def placement_per_edge(self, negative_branch):
        for e in self.tree.postorder_edge_iter():
            if e.head_node == self.tree.seed_node:
                continue
            a_11 = e.R1_D + e.S1_D
            a_12 = e.R1_D - e.S1_D
            a_21 = a_12
            a_22 = a_11
            c_1 = e.R + e.S - e.length * e.S1_D - e.Rd_D - e.Sd_D
            c_2 = e.R - e.S + e.length * e.S1_D - e.Rd_D + e.Sd_D
            util.solve2_2(e, a_11, a_12, a_21, a_22, c_1, c_2, negative_branch)

    # computes BE error (Q value) for a given edge
    @staticmethod
    def error_per_edge(edge):
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
