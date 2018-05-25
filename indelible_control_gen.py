import numpy as np
import scipy as sci
import sys
from optparse import OptionParser


if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-t", "--tree", dest="tree_fp",
                      help="file containing all gene trees, one tree per line",
                      metavar="FILE")
    parser.add_option("-o", "--output", dest = "output_fp",
                      help="path to the output file", metavar="FILE")

    parser.add_option("-l", "--length", dest = "seq_length",
                      help = "sequence length", metavar = "INT")

    (options, args) = parser.parse_args()
    tree_fp = options.tree_fp
    output_fp = options.output_fp
    seq_length = options.seq_length
    
    with open(tree_fp) as f:
        trees = f.readlines()

    with open(output_fp, "w") as f:
        f.write("[TYPE] NUCLEOTIDE 1\n")
        digits = len(str(len(trees)))
        for i in range(len(trees)):
            pi = np.random.dirichlet([36,26,28,32])
            gtr = np.random.dirichlet([600,300,500,500,600,150])
            alpha = np.random.exponential(1.2)
            if alpha < 0.1:
                alpha = 0
            f.write("[MODEL] GTR" + str(i).zfill(digits) + "\n")
            f.write("\t [submodel] GTR")
            for j in gtr:
                f.write(" " + "{:.6f}".format(j))
            f.write("\n\t [statefreq]")
            for j in pi:
                f.write(" " + "{:.6f}".format(j))
            f.write("\n\t [rates] 0 " + "{:.6f}".format(alpha) + " 0 \n")
        f.write("[SETTINGS]\n\t[randomseed] 2478\n\t[fileperrep] FALSE\n")
        for i in range(len(trees)):
            f.write("[TREE] T" + str(i).zfill(digits) + " " + trees[i])
        for i in range(len(trees)):
            f.write("[PARTITIONS] T" + str(i).zfill(digits) + " " + "[T"
                    + str(i).zfill(digits) + " GTR" + str(i).zfill(digits)
                    + " " + seq_length + "]\n")
        f.write("[EVOLVE]\n")
        for i in range(len(trees)):
            f.write("T" + str(i).zfill(digits) + " 1 " + str(i).zfill(digits)
                    + "\n")
        
        


