from collections import deque

from apples.PrioritySet import PrioritySet


class Subtree:

    def __init__(self, obs_dist, name_to_node_map, query_name):
        self.query_name = query_name
        self.root = self.validate_edges(obs_dist, name_to_node_map, True)
        self.obs_dist = obs_dist
        self.name_to_node_map = name_to_node_map

    def edge_contribution_counter(self):
        contr_map = dict()
        for k, v in self.obs_dist.items():
            if k in self.name_to_node_map:
                edge_count = 0
                current = self.name_to_node_map[k]
                while True:
                    edge_count += 1
                    if current.parent.visitcounts > 1:
                        break
                    current = current.parent
                if current.parent == self.root:
                    def descend(parent,child):
                        if child.visitcounts > 1:
                            return 1
                        validgrandchild = [c for c in child.children if c.valid][0]
                        return 1 + descend(child, validgrandchild)
                    sibling = [c for c in self.root.children if c != current and c.valid][0]
                    edge_count += descend(self.root, sibling)
                contr_map[k] = edge_count
        return contr_map

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
            if not valid:
                x.parent.visitcounts = 0
            else:
                if hasattr(x.parent, "visitcounts"):
                    x.parent.visitcounts += 1
                else:
                    x.parent.visitcounts = 1
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
