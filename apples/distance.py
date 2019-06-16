import numpy as np


def jc69(a2,b2):
    nondash=np.logical_and(a2!=b'-',b2!=b'-')
    valid = np.count_nonzero(nondash)
    if not valid:
        return 5.0
    p = np.count_nonzero(np.logical_and(a2!=b2, nondash))*1.0/valid
    if (p - np.finfo(float).eps < 0):
        return 0.0
    else:
        loc = 1 - (4 * p / 3)
        if (0 >= loc):
            return  5.0
        else:
            return -0.75*np.log(loc)
