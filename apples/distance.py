import numpy as np

def jc69(a,b):
    mm = 0
    valid = 0
    for c1,c2 in zip(a,b):
        if c1 == "-" or c2 == "-":
            continue
        else:
            valid+=1
            if c1 != c2:
                mm += 1
    if not valid :
        return (5.0)
    p = mm*1.0/valid
    if (abs(p) - np.finfo(float).eps < 0):
        return 0.0

    loc = 1 - (4 * p / 3)

    if (0 >= loc):
        return (5.0)

    return -0.75 * (np.log(loc))
