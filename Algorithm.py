
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

    def placement(self, query_name, selection_name):
        if selection_name == "HYBRID":
            sm = heapq.nsmallest(math.floor(math.log2(len(self.tree.leaf_nodes()))),
                                 self.tree.postorder_edge_iter(filter_fn=lambda e: e.head_node != self.tree.seed_node),
                                 key=lambda e: self.error_per_edge(e))
            placed_edge = min(sm, key=lambda e: e.x_1_neg)
            util.insert(self.tree, placed_edge, query_name, placed_edge.x_1_neg, placed_edge.x_2_neg)
        elif selection_name == "ME":
            placed_edge = min(self.tree.postorder_edge_iter(filter_fn=lambda e: e.head_node != self.tree.seed_node),
                              key=lambda e: e.x_1_neg)
            util.insert(self.tree, placed_edge, query_name, placed_edge.x_1_neg, placed_edge.x_2_neg)
        else:  # selection_name == "MLSE"
            placed_edge = min(self.tree.postorder_edge_iter(filter_fn=lambda e: e.head_node != self.tree.seed_node),
                              key=lambda e: self.error_per_edge(e))
            util.insert(self.tree, placed_edge, query_name, placed_edge.x_1, placed_edge.x_2)
        return self.tree