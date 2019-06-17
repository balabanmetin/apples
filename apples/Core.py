import numpy as np
from apples import util
from scipy.special import binom
import heapq
import math

class Core:
    def __init__(self, tree, k):
        self.k = k
        self.tree = tree


    def dp(self, obs_dist):
        self.dpm = np.zeros((2, self.tree.num_nodes(), 3, 3)) # S-R, nodes, tree,  obs
        self.xs = np.zeros((self.tree.num_nodes(), 2, 2)) #  nodes, neg,  x
        for node in self.tree.traverse_postorder(): #b=0
            for a in range(3):
                if node.is_leaf():
                    delta = obs_dist[node.label]
                    self.dpm[0,node.edge_index,0,a] = delta ** (a + self.k)
                else:
                    if node != self.tree.root:
                        for child in node.children:
                            self.dpm[0,node.edge_index,0,a] += self.dpm[0,child.edge_index,0,a]
        for node in self.tree.traverse_postorder():
            for b in range(1,3):
                for a in range(3):
                    if node.is_leaf():
                        self.dpm[0,node.edge_index,b,a] = 0
                    else:
                        if node != self.tree.root:
                            for child in node.children:
                                for j in range(b+1):
                                    self.dpm[0,node.edge_index,b,a] += ((child.edge_length) ** j) * binom(b, j) * \
                                                                          self.dpm[0,child.edge_index,b - j,a]

        for node in self.tree.traverse_preorder(): #TODO make algorithm linear time for star-tree
            for b in range(0,3):
                for a in range(3):
                    if node != self.tree.root:
                        for j in range(b + 1):
                            for sibling in node.parent.children:
                                if sibling != node:
                                    self.dpm[1,node.edge_index,b,a] += ((sibling.edge_length) ** j) * binom(b, j) * \
                                                                          self.dpm[0,sibling.edge_index,b - j,a]
                            if node.parent != self.tree.root:
                                self.dpm[1,node.edge_index,b,a] += ((node.parent.edge_length) ** j) * binom(b, j) * \
                                                                      self.dpm[1,node.parent.edge_index,b - j,a]

    def placement_per_edge(self, negative_branch):
        dpm = self.dpm
        tree = self.tree
        for node in tree.traverse_postorder():
            if node == tree.root:
                continue
            i = node.edge_index
            a_11 = dpm[1,i,0,0] + dpm[0,i,0,0]
            a_12 = dpm[1,i,0,0] - dpm[0,i,0,0]
            a_21 = a_12
            a_22 = a_11
            c_1 = dpm[1,i,0,1] + dpm[0,i,0,1] - node.edge_length * dpm[0,i,0,0] - dpm[1,i,1,0] - dpm[0,i,1,0]
            c_2 = dpm[1,i,0,1] - dpm[0,i,0,1] + node.edge_length * dpm[0,i,0,0] - dpm[1,i,1,0] + dpm[0,i,1,0]
            self.xs[i] = util.solve2_2(i, a_11, a_12, a_21, a_22, c_1, c_2, negative_branch)

    def error_per_edge(self, node):
        dpm = self.dpm
        xs = self.xs
        i = node.edge_index
        return dpm[1,i,0,2] + dpm[0,i,0,2] + dpm[1,i,2,0] + dpm[0,i,2,0] + \
               2 * (xs[i][0][0] + xs[i][0][1]) * (dpm[1, i, 1, 0] - dpm[1, i, 0, 1]) + \
               2 * (xs[i][0][0] + node.edge_length - xs[i][0][1]) * (dpm[0, i, 1, 0] - dpm[0, i, 0, 1]) + \
               (xs[i][0][0] + xs[i][0][1])**2 * dpm[1, i, 0, 0] + (xs[i][0][0] + node.edge_length - xs[i][0][1])**2 * dpm[0, i, 0, 0] - \
               2 * dpm[0, i, 1, 1] - 2 * dpm[1, i, 1, 1]

            #x = np.dot(np.linalg.inv(A),C)
            #print(x)

    def placement(self, selection_name):
        valids = filter(lambda x : x != self.tree.root, self.tree.traverse_postorder())
        if selection_name == "HYBRID":
            sm = heapq.nsmallest(math.floor(math.log2(len(self.tree.num_nodes(internal=False)))),
                                 valids,
                                 key=lambda e: self.error_per_edge(e))
            placed_edge = min(sm, key=lambda e: self.xs[e.edge_index][1][0])
            pendant, relative_distal = self.xs[placed_edge.edge_index][1]

        elif selection_name == "ME":
            placed_edge = min(valids, key=lambda e: self.xs[e.edge_index][1][0])
            pendant, relative_distal = self.xs[placed_edge.edge_index][1]
        else:  # selection_name == "MLSE"
            placed_edge = min(valids, key=lambda e: self.error_per_edge(e))
            pendant, relative_distal = self.xs[placed_edge.edge_index][1]
        return [placed_edge.edge_index, 0, 1, placed_edge.edge_length - relative_distal, pendant]
