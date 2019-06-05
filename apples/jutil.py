
def join_jplace(lst):
    if len(lst) == 1:
        return lst[0]
    else:
        result = lst[0]
        for i in range(1,len(lst)):
            result["placements"] = result["placements"] + lst[i]["placements"]
        return result

def extended_newick(tree):
    """Newick printing algorithm is based on treeswift"""

    if tree.seed_node.edge_length is None:
        suffix = ''
    elif isinstance(tree.seed_node.edge_length, int):
        suffix = ':%d' % tree.seed_node.edge_length
    elif isinstance(tree.seed_node.edge_length, float) and tree.seed_node.edge_length.is_integer():
        suffix = ':%d' % int(tree.seed_node.edge_length)
    else:
        suffix = ':%s' % str(tree.seed_node.edge_length)
    strng = _nodeprint(tree.seed_node)
    count = tree.seed_node.edge.edge_index
    if tree.is_rooted:
        return '[&R] %s%s{%d};' % (strng, suffix, count)
    else:
        return '%s%s{%d};' % (strng, suffix, count)


def _nodeprint(root):
    node_to_str = dict()
    counter=0

    for node in root.postorder_iter():
        node.edge.edge_index = counter
        counter += 1

    for node in root.postorder_iter():
        if node.is_leaf():
            if node.taxon.label is None:
                node_to_str[node] = ''
            else:
                node_to_str[node] = str(node.taxon.label)
        else:
            out = ['(']
            for c in node.child_node_iter():
                out.append(node_to_str[c])
                if c.edge_length is not None:
                    if isinstance(c.edge_length, int):
                        l_str = str(c.edge_length)
                    elif isinstance(c.edge_length, float) and c.edge_length.is_integer():
                        l_str = str(int(c.edge_length))
                    else:
                        l_str = str(c.edge_length)
                    out.append(':%s' % l_str)
                out.append('{%d}' % c.edge.edge_index)
                out.append(',')
                del node_to_str[c]
            out.pop()  # trailing comma
            out.append(')')
            if node.label is not None:
                out.append(str(node.label))
            node_to_str[node] = ''.join(out)
    return node_to_str[root]
