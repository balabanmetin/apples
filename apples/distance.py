import numpy as np
np.seterr(all='raise')

#  Retrieved from belvu
#  https://github.com/richarddurbin/acedb/blob/86d8c1c92d8a2c58cd090e85a1e17771612bfcb9/w9/belvu.c
#  BLOSUM62 930809
#  Matrix made by matblas from blosum62.iij
#  * column uses minimum score
#  BLOSUM Clustered Scoring Matrix in 1/2 Bit Units
#  Blocks Database = /data/blocks_5.0/blocks.dat
#  Cluster Percentage: >= 62
#  Entropy =   0.6979, Expected =  -0.5209
#  Note: to use with a2b[], always subtract 1 from the values !!!!
#  A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   b   z   x
BLOSUM62 = np.array([
   4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0,
  -1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3,
  -2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3,
  -2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3,
   0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1,
  -1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2,
  -1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2,
   0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3,
  -2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3,
  -1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3,
  -1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1,
  -1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2,
  -1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1,
  -2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1,
  -1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2,
   1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2,
   0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0,
  -3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3,
  -2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1,
   0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4,
 ])


NA=0
a2i = np. array([
    NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
    NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
    NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
    NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
    NA, 0,NA, 4, 3, 6,13, 7, 8, 9,NA,11,10,12, 2,NA,
    14, 5, 1,15,16,NA,19,17,NA,18,NA,NA,NA,NA,NA,NA,
    NA, 0,NA, 4, 3, 6,13, 7, 8, 9,NA,11,10,12, 2,NA,
    14, 5, 1,15,16,NA,19,17,NA,18,NA,NA,NA,NA,NA,NA,

    NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
    NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
    NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
    NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
    NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
    NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
    NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
    NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA
  ])


def scoredist(a2, b2):
    nondash = np.logical_and(a2 != b'-', b2 != b'-')
    valid = np.count_nonzero(nondash)
    if not valid or valid / len(nondash) < 0.001:
        return -1.0  # treat as missing data
    a_as_int = a2i[a2.view(np.uint8)]
    b_as_int = a2i[b2.view(np.uint8)]
    aa_ind = 21 * a_as_int
    bb_ind = 21 * b_as_int
    ab_ind = 20 * a_as_int + b_as_int
    expect = -0.5209 * valid
    aa_tot = np.sum(np.dot(nondash, BLOSUM62[aa_ind]))
    bb_tot = np.sum(np.dot(nondash, BLOSUM62[bb_ind]))
    ab_tot = np.sum(np.dot(nondash, BLOSUM62[ab_ind]))
    maxsc = (aa_tot + bb_tot) / 2.0
    od = (ab_tot - expect) / (maxsc - expect)
    if 0 >= od:
        return -1.0  # treat as missing data
    # if od < 0.05:
    #     od = 0.05  # Limit to 300 PAM;  len==0 if no overlap */
    if od > 1.0:
        od = 1.0

    cd = -np.log(od)
    cd = cd * 1.337  # Magic scaling factor optimized for Dayhoff data
    # if cd > 3.0:
    #     cd = 3.0
    return cd


def jc69(a2, b2):
    nondash = np.logical_and(a2 != b'-', b2 != b'-')
    valid = np.count_nonzero(nondash)
    if not valid or valid / len(nondash) < 0.001:
        return -1.0  # treat as missing data
    p = np.count_nonzero(np.logical_and(a2 != b2, nondash)) * 1.0 / valid
    if (p - np.finfo(float).eps < 0):
        return 0.0
    else:
        loc = 1 - (4 * p / 3)
        if 0 >= loc:
            return -1.0  # treat as missing data
        else:
            return -0.75 * np.log(loc)