from abc import ABC, abstractmethod
import heapq
import math


class Algorithm(ABC):
    def __init__(self, subtree):
        """
        Initialize an Algorithm instance with a given subtree object.

        Args:
            subtree: The subtree to which the algorithm will be applied.
        """
        self.subtree = subtree

    def dp_frag(self):
        """
        Perform dynamic programming to calculate all S and R values.
        """
        self.all_S_values()
        self.all_R_values()

    @abstractmethod
    def all_S_values(self):
        """
        Abstract method to calculate all S values.
        Must be implemented by subclasses.
        """
        pass

    @abstractmethod
    def all_R_values(self):
        """
        Abstract method to calculate all R values.
        Must be implemented by subclasses.
        """
        pass

    @abstractmethod
    def placement_per_edge(self, negative_branch):
        """
        Abstract method to calculate the best pendant and distal edge lengths for all branches.

        Args:
            negative_branch: Whether negative branches are allowed or not.
        """
        pass

    @staticmethod
    def error_per_edge(edge):
        """
        Static method to calculate the weighted least error for the best placement of an edge.

        Args:
            edge: The edge for which to calculate the error.

        Returns:
            The error associated with the given edge.
        """
        pass

    def placement(self, selection_name):
        """
        Determine the best placement based on a selection strategy.

        Args:
            selection_name: A string representing the selection strategy.
                            Can be 'HYBRID', 'ME', or 'MLSE'.

        Returns:
            A tuple containing placement information and a flag indicating
            potential misplacement.
        """
        valids = filter(lambda x: x.valid, self.subtree.traverse_postorder())

        if selection_name == 'HYBRID':
            sm = heapq.nsmallest(
                math.floor(math.log2(self.subtree.num_nodes)), valids, key=lambda e: self.error_per_edge(e)
            )
            placed_edge = min(sm, key=lambda e: e.x_1)
            pendant = placed_edge.x_1
            relative_distal = placed_edge.x_2

        elif selection_name == 'ME':
            placed_edge = min(valids, key=lambda e: e.x_1)
            pendant = placed_edge.x_1
            relative_distal = placed_edge.x_2
        else:  # selection_name == "MLSE"
            placed_edge = min(valids, key=lambda e: self.error_per_edge(e))
            pendant = placed_edge.x_1
            relative_distal = placed_edge.x_2

        e_per_edge = self.error_per_edge(placed_edge)
        if pendant == 0 and e_per_edge > 0 and (relative_distal == 0 or relative_distal == placed_edge.edge_length):
            potential_misplacement_flag = 1
        else:
            potential_misplacement_flag = 0
        return (
            [placed_edge.edge_index, e_per_edge, 1, placed_edge.edge_length - relative_distal, pendant],
            potential_misplacement_flag,
        )
