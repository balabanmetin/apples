from collections import deque

from apples.PrioritySet import PrioritySet


class Subtree:

    def __init__(self, obs_dist, name_to_node_map):
        self.root = self.validate_edges(obs_dist, name_to_node_map, True)
        self.obs_dist = obs_dist
        self.name_to_node_map = name_to_node_map

    def validate_edges(self, obs_dist, name_to_node_map, valid):
        ps = PrioritySet()
        for k, v in obs_dist.items():
            if k in name_to_node_map:
                ps.add(name_to_node_map[k],
                       -name_to_node_map[k].level)  # adding minus to give priority to the largest level

        count = 0
        while len(ps) > 1:
            x = ps.get()
            x.valid = valid
            count += 1
            ps.add(x.parent, -x.parent.level)
        self.num_nodes = count
        return ps.get()

    def traverse_preorder(self):
        s = deque()
        s.append(self.root)
        while len(s) != 0:
            n = s.pop()
            s.extend([c for c in n.children if c.valid])
            yield n

    def traverse_postorder(self):
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
        self.validate_edges(self.obs_dist, self.name_to_node_map, False)