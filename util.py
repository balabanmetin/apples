

# solves two by two Ax=c linear system and returns the optimal values defined by
# constraints. This is thanks to the convexity of (weighted) least squared error
def solve2_2(e, a_11, a_12, a_21, a_22, c_1, c_2, negative_branch):
    edge_length = e.length
    det = 1 / (a_11 * a_22 - a_12 * a_21)
    assert det is not 0
    x_1_neg = (a_22 * c_1 - a_12 * c_2) * det
    x_2_neg = (- a_21 * c_1 + a_11 * c_2) * det
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
    e.x_1 = x_1
    e.x_2 = x_2
    e.x_1_neg = x_1_neg
    e.x_2_neg = x_2_neg