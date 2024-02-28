# solves two by two Ax=c linear system and returns the optimal values defined by
# constraints. This is thanks to the convexity of (weighted) least squared error
from collections import deque


def solve2_2(node, a_11, a_12, a_21, a_22, c_1, c_2, negative_branch):
    """
    Solves a system of linear equations for a 2x2 matrix and updates the
    node's attributes. It solves a system of two linear equations with two
    variables using the Cramer's rule. It calculates the values of x1 and x2
    based on the input coefficients and constants, and assigns the results to
    the node attributes. If negative branches are not allowed, it uses
    convexity properties to find the non-negative branch lengths that minimize
    the error.

    Args:
    - node: the node object to update
    - a_11, a_12, a_21, a_22: coefficients of the linear equations
    - c_1, c_2: constants of the linear equations
    - negative_branch: a boolean indicating whether to use the negative branch

    Returns:
    - None
    """
    edge_length = node.edge_length
    det = 1 / (a_11 * a_22 - a_12 * a_21)
    assert det != 0
    x_1_neg = (a_22 * c_1 - a_12 * c_2) * det
    x_2_neg = (-a_21 * c_1 + a_11 * c_2) * det
    x_1 = x_1_neg
    x_2 = x_2_neg
    if not negative_branch:
        if x_1_neg < 0 and x_2_neg < 0:
            x_1 = 0
            x_2 = 0
        elif x_1_neg > 0 and x_2_neg < 0:
            x_1 = max(c_1 * 1.0 / a_11, 0)
            x_2 = 0
        elif x_1_neg < 0 and 0 <= x_2_neg and x_2_neg <= edge_length:
            x_1 = 0
            x_2 = min(max(c_2 * 1.0 / a_22, 0), edge_length)
        elif x_1_neg < 0 and x_2_neg > edge_length:
            x_1 = 0
            x_2 = edge_length
        elif x_1_neg > 0 and x_2_neg > edge_length:
            x_1 = max((c_1 * 1.0 - a_12 * edge_length) / a_11, 0)
            x_2 = edge_length
        else:
            x_1 = x_1
            x_2 = x_2
    node.x_1 = x_1
    node.x_2 = x_2
    node.x_1_neg = x_1_neg
    node.x_2_neg = x_2_neg


def index_edges(tree):
    """
    Function to index the edges of a tree. It takes a tree as input. It initializes
    a counter to 0 and then iterates through the nodes of the tree in postorder.
    For each node, it assigns the current counter value to the edge_index attribute
    of the node, increments the counter, and sets the valid attribute of the node
    to False.
    """
    counter = 0
    for node in tree.traverse_postorder():
        node.edge_index = counter
        counter += 1
        node.valid = False


def set_levels(tree):
    """
    Sets the level attribute for each node in the tree based on their distance from the root node.
    Parameters:
        tree: the tree data structure to operate on
    Return:
        None
    """
    root = tree.root
    root.level = 0
    q = deque()
    q.append(root)
    while len(q) != 0:
        n = q.popleft()
        for c in n.children:
            c.level = n.level + 1
        q.extend(n.children)
