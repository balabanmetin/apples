
def join_jplace(lst):
    result = lst[0]
    if len(lst) == 1:
        if result["placements"][0]["p"][0][0] == -1:
            result["placements"] = []
    else:
        for i in range(1,len(lst)):
            if lst[i]["placements"][0]["p"][0][0] != -1:
                result["placements"] = result["placements"] + lst[i]["placements"]
    return result

def extended_newick(tree):
    """Newick printing algorithm is based on treeswift"""

    # if tree.root.edge_length is None:
    #     suffix = ''
    # elif isinstance(tree.root.edge_length, int):
    #     suffix = ':%d' % tree.root.edge_length
    # elif isinstance(tree.root.edge_length, float) and tree.root.edge_length.is_integer():
    #     suffix = ':%d' % int(tree.root.edge_length)
    # else:
    #     suffix = ':%s' % str(tree.root.edge_length)
    suffix = ''
    strng = _nodeprint(tree.root)
    if tree.is_rooted:
        return '[&R] %s%s;' % (strng, suffix)
    else:
        return '%s%s;' % (strng, suffix)


def _nodeprint(root):
    node_to_str = dict()

    for node in root.traverse_postorder():
        if node.is_leaf():
            if node.label is None:
                node_to_str[node] = ''
            else:
                node_to_str[node] = str(node.label)
        else:
            out = ['(']
            for c in node.children:
                out.append(node_to_str[c])
                if c.edge_length is not None:
                    if isinstance(c.edge_length, int):
                        l_str = str(c.edge_length)
                    elif isinstance(c.edge_length, float) and c.edge_length.is_integer():
                        l_str = str(int(c.edge_length))
                    else:
                        l_str = str(c.edge_length)
                    out.append(':%s' % l_str)
                out.append('{%d}' % c.edge_index)
                out.append(',')
                del node_to_str[c]
            out.pop()  # trailing comma
            out.append(')')
            if node.label is not None:
                out.append(str(node.label))
            node_to_str[node] = ''.join(out)
    return node_to_str[root]
