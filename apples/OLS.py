from apples.Algorithm import Algorithm
from apples import util


# Glossary
# OLS: Ordinary Least Squares
# FM: Fitch-Margoliash --Least Squares with \delta^2 scaling
# BE: Beyer --Least Squares with \delta scaling


class OLS(Algorithm):
    def all_S_values(self):
        """
        Calculate and assign the S, Sd, Sd2, SDd, SD, and SD2 values
        (left-hand side dynamic programming values) for each node in the subtree.
        It is computing certain summary statistics (S values) for each node in a tree,
        which will be used to compute the least squares error for the edge.
        These S values are calculated based on the node's siblings and parent,
        and they are used in dynamic programming on the tree.
        The calculations involve updating the S values for each node based on the node's siblings and parent,
        and the code id traversing the tree in a specific order to perform these calculations
        """
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
                    node.Sd2 += (
                        child.S * child.edge_length * child.edge_length + child.Sd2 + 2 * child.edge_length * child.Sd
                    )
                    node.SDd += child.edge_length * child.SD + child.SDd
                    node.SD2 += child.SD2
                    node.SD += child.SD

    def all_R_values(self):
        """
        Calculate the R values (right-hand side dynamic programming values) for each node in the subtree.
        It is computing certain summary statistics (R values) for each node in a tree,
        which will be used to compute the least squares error for the edge.
        These R values are calculated based on the node's siblings and parent,
        and they are used in dynamic programming on the tree.
        The calculations involve updating the R values for each node based on the node's siblings and parent,
        and the code id traversing the tree in a specific order to perform these calculations
        """
        subtree = self.subtree
        for node in filter(lambda x: x.valid, subtree.traverse_preorder()):
            node.RDd, node.RD2, node.RD, node.R, node.Rd, node.Rd2 = 6 * [0]
            for sibling in filter(lambda x: x.valid and x != node, node.parent.children):
                node.R += sibling.S
                node.Rd += sibling.S * sibling.edge_length + sibling.Sd
                node.Rd2 += (
                    sibling.S * sibling.edge_length * sibling.edge_length
                    + sibling.Sd2
                    + 2 * sibling.edge_length * sibling.Sd
                )
                node.RDd += sibling.SD * sibling.edge_length + sibling.SDd
                node.RD2 += sibling.SD2
                node.RD += sibling.SD
            if node.parent != subtree.root and node.parent.valid:
                node.R += node.parent.R
                node.Rd += node.parent.R * node.parent.edge_length + node.parent.Rd
                node.Rd2 += (
                    node.parent.R * node.parent.edge_length * node.parent.edge_length
                    + node.parent.Rd2
                    + 2 * node.parent.edge_length * node.parent.Rd
                )
                node.RDd += node.parent.RD * node.parent.edge_length + node.parent.RDd
                node.RD2 += node.parent.RD2
                node.RD += node.parent.RD

    # computes all alternative OLS placements
    def placement_per_edge(self, negative_branch):
        """
        Calculate placement per edge for the given the flag for the negative branches allowed or not.
        It calculates placement per edge for a given subtree by solving a 2x2 linear system for each node (edge)
        in the subtree. It calculates two values (implicitly), pendant and distal edge lengths for each edge.
        The values are stored in the node attributes.
        """
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
        """
        Calculate the error per edge based on the given node attributes and return the result.
        The error is defined as the sum of the following terms:
        A = RD2 + SD2
        B = 2 * (x_1 + x_2) * Rd + 2 * (edge_length + x_1 - x_2) * Sd
        C = (x_1 + x_2)^2 * R + (edge_length + x_1 - x_2)^2 * S
        D = -2 * (x_1 + x_2) * RD - 2 * (edge_length + x_1 - x_2) * SD
        E = -2 * RDd - 2 * SDd
        F = Rd2 + Sd2

        This is derived by writing the least squares objective explicitly.

        Parameters:
        - node: the node object containing attributes such as valid, RD2, SD2, x_1, x_2, Rd, Sd,
        edge_length, R, S, RD, SD, RDd, SDd, Rd2, and Sd2

        Returns:
        - The calculated error per edge
        """
        assert node.valid
        A = node.RD2 + node.SD2
        B = 2 * (node.x_1 + node.x_2) * node.Rd + 2 * (node.edge_length + node.x_1 - node.x_2) * node.Sd
        C = (node.x_1 + node.x_2) ** 2 * node.R + (node.edge_length + node.x_1 - node.x_2) ** 2 * node.S
        D = -2 * (node.x_1 + node.x_2) * node.RD - 2 * (node.edge_length + node.x_1 - node.x_2) * node.SD
        E = -2 * node.RDd - 2 * node.SDd
        F = node.Rd2 + node.Sd2
        return A + B + C + D + E + F
