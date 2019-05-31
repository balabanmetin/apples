
from abc import ABC, abstractmethod
import heapq
import math
import util

class Algorithm(ABC):
    def __init__(self, tree):
        self.tree = tree

    @abstractmethod
    def placement_per_edge(self, negative_branch):
        pass

    @staticmethod
    def error_per_edge(edge):
        pass

    def placement(self, selection_name):

        def find_distal_pendant(placed_edge, x_1, x_2):
            pendant = x_1
            if placed_edge.head_node in list(
                    self.tree.master_edge.head_node.ancestor_iter()) or self.tree.master_edge == placed_edge:
                distal = max(x_2, 0)
            else:
                distal = placed_edge.length - max(x_2, 0)
            return [placed_edge.edge_index, 0, 1, distal, pendant]

        if selection_name == "HYBRID":
            sm = heapq.nsmallest(math.floor(math.log2(len(self.tree.leaf_nodes()))),
                                 self.tree.postorder_edge_iter(filter_fn=lambda e: e.head_node != self.tree.seed_node),
                                 key=lambda e: self.error_per_edge(e))
            placed_edge = min(sm, key=lambda e: e.x_1_neg)
            pendant = placed_edge.x_1_neg
            relative_distal = placed_edge.x_2_neg

        elif selection_name == "ME":
            placed_edge = min(self.tree.postorder_edge_iter(filter_fn=lambda e: e.head_node != self.tree.seed_node),
            key=lambda e: e.x_1_neg)
            pendant = placed_edge.x_1_neg
            relative_distal = placed_edge.x_2_neg
        else:  # selection_name == "MLSE"
            placed_edge = min(self.tree.postorder_edge_iter(filter_fn=lambda e: e.head_node != self.tree.seed_node),
                              key=lambda e: self.error_per_edge(e))
            pendant = placed_edge.x_1
            relative_distal = placed_edge.x_2

        return find_distal_pendant(placed_edge, pendant, relative_distal)
