import itertools
import numpy as np
import multiprocessing as mp
import time


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
    p = mm*1.0/valid
    if (abs(p) - np.finfo(float).eps < 0):
        return 0.0

    loc = 1 - (4 * p / 3)

    if (0 >= loc):
        return (5.0)

    return -0.75 * (np.log(loc))

def calc_mp(indices, data, num_thread):
    # construct pool
    if num_thread:
        pool = mp.Pool(num_thread)
    else:
        pool = mp.Pool(mp.cpu_count())

    # we are going to populate the matrix; organize all the inputs; then map them
    matrix = [[0] * len(data) for i in range(len(indices))]
    args = [(data[i_a], data[i_b]) for i_a, i_b in itertools.product(indices, range(len(data)))]
    #time_start = time.time()
    results = pool.starmap(jc69, args)
    #print("TOOK {}".format(time.time() - time_start))

    # unpack the results into the matrix
    for i_tuple, result in zip([(i_com, i_b) for i_com, i_b in list(itertools.product(enumerate(indices), range(len(data))))], results):
        # unpack
        i_com, i_b = i_tuple
        i_ind, i_a = i_com
        # set it in the matrix
        matrix[i_ind][i_b] = result

    return matrix