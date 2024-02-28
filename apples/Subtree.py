from collections import deque

from apples.PrioritySet import PrioritySet


class Subtree:
    def __init__(self, obs_dist, name_to_node_map):
        """
        Initializes the subtree with the given observed distances and name to node map.
         It validates the minumum number of edges required to select the subtree.

        Parameters:
            obs_dist (type): The observation distribution.
            name_to_node_map (type): The map of node names to node objects.

        Returns:
            None
        """
        self.root = self.validate_edges(obs_dist, name_to_node_map, True)
        self.obs_dist = obs_dist
        self.name_to_node_map = name_to_node_map

    def validate_edges(self, obs_dist, name_to_node_map, valid):
        """
        Validates edges based on the given observed distances. It finds the subtree
        spun by the nodes in the gived observed distances map. It returns the root.
        It can be used to select the subtree or unselect it by toggling valid flag.
        """
        ps = PrioritySet()
        for k, v in obs_dist.items():
            if k in name_to_node_map:
                ps.add(
                    name_to_node_map[k], -name_to_node_map[k].level
                )  # adding minus to give priority to the largest level

        count = 0
        while len(ps) > 1:
            x = ps.get()
            x.valid = valid
            count += 1
            ps.add(x.parent, -x.parent.level)
        self.num_nodes = count
        return ps.get()

    def traverse_preorder(self):
        """
        Traverse the tree in preorder fashion, yielding each node as it is visited.
        """
        s = deque()
        s.append(self.root)
        while len(s) != 0:
            n = s.pop()
            s.extend([c for c in n.children if c.valid])
            yield n

    def traverse_postorder(self):
        """
        Traverse the tree in postorder fashion and yield each node.
        """
        s1 = deque()
        s2 = deque()
        s1.append(self.root)
        while len(s1) != 0:
            n = s1.pop()
            s2.append(n)
            s1.extend([c for c in n.children if c.valid])

        while len(s2) != 0:
            n = s2.pop()
            yield n

    def unroll_changes(self):
        """
        Method to unroll changes, and validate edges of the observation distribution and name to node map.
        """
        self.validate_edges(self.obs_dist, self.name_to_node_map, False)
