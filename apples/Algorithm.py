
from abc import ABC, abstractmethod
import heapq
import math


class Algorithm(ABC):
    def __init__(self, subtree):
        self.subtree = subtree

    def dp_frag(self):
        self.all_S_values()
        self.all_R_values()

    @abstractmethod
    def all_S_values(self):
        pass

    @abstractmethod
    def all_R_values(self):
        pass

    @abstractmethod
    def placement_per_edge(self, negative_branch):
        pass

    @staticmethod
    def error_per_edge(edge):
        pass

    def placement(self, selection_name):

        valids = filter(lambda x: x.valid, self.subtree.traverse_postorder())

        if selection_name == "HYBRID":
            sm = heapq.nsmallest(math.floor(math.log2(self.subtree.num_nodes)),
                                 valids,
                                 key=lambda e: self.error_per_edge(e))
            placed_edge = min(sm, key=lambda e: e.x_1)
            pendant = placed_edge.x_1
            relative_distal = placed_edge.x_2

        elif selection_name == "ME":
            placed_edge = min(valids, key=lambda e: e.x_1)
            pendant = placed_edge.x_1
            relative_distal = placed_edge.x_2
        else:  # selection_name == "MLSE"
            placed_edge = min(valids, key=lambda e: self.error_per_edge(e))
            pendant = placed_edge.x_1
            relative_distal = placed_edge.x_2

        return [placed_edge.edge_index, self.error_per_edge(placed_edge),
                1, placed_edge.edge_length - relative_distal, pendant]
