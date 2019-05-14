#!/usr/bin/env python3

import dendropy as dy
import subprocess
import sys
from sys import platform as _platform
from optparse import OptionParser
import re
import os.path
import tempfile
import util
from core import Core
from criteria import OLS, FM, BE



if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-t", "--tree", dest="tree_fp",
                      help="path to the reference tree", metavar="FILE")
    parser.add_option("-d", "--distances", dest="dist_fp",
                      help="path to the table of observed distances", metavar="FILE")
    parser.add_option("-a", "--aln", dest="aln_fp",
                      help="path to the alignment file (PHYLIP), containing reference and query sequence",
                      metavar="FILE")
    parser.add_option("-m", "--method", dest="method_name", default="FM",
                      help="name of the weighted least squares method (OLS, FM, or BE)", metavar="METHOD")
    parser.add_option("-s", "--selection", dest="selection_name", default="MLSE",
                      help="name of the placement selection criteria (MLSE, ME, or HYBRID", metavar="CRITERIA")
    parser.add_option("-n", "--negative", dest="negative_branch", action='store_true',
                      help="relaxes positivity constraint on new branch lengths, i.e. allows negative branch lengths")
    parser.add_option("-p", "--protein", dest="protein_seqs", action='store_true',
                      help="input sequences are protein sequences")
    parser.add_option("-P", "--placement", dest="placement_set_size", default="-1",
                      help="number of query sequences. First P sequences in the input alignment must be \
                           queries, and the rest must be backbone sequences", metavar="NUMBER")

    (options, args) = parser.parse_args()
    tree_fp = options.tree_fp
    dist_fp = options.dist_fp
    aln_fp = options.aln_fp
    method_name = options.method_name
    selection_name = options.selection_name
    negative_branch = options.negative_branch
    protein_seqs = options.protein_seqs
    placement_set_size = options.placement_set_size


    f = open(tree_fp)
    tree_string = f.readline()
    f.close()

    if aln_fp:
        if dist_fp:
            raise ValueError('Input should be either an alignment or a distance matrix, but not both!')
        if protein_seqs:
            datatype='p'
        else:
            datatype='d'
        if placement_set_size == "-1":
            raise ValueError('Number of queries is not specified!')
        else:
            query_size = "-Q {} ".format(placement_set_size)
        nldef = tempfile.NamedTemporaryFile(delete=True, mode='w+t')
        if _platform == "darwin":
            fastme_exec = os.path.join(os.path.dirname(__file__), 'tools/fastme-darwin64')
        elif _platform == "linux" or _platform == "linux2":
            fastme_exec = os.path.join(os.path.dirname(__file__), 'tools/fastme-linux64')
        elif _platform == "win32" or _platform == "win64" or _platform == "msys":
            fastme_exec = os.path.join(os.path.dirname(__file__), 'tools/fastme-win64.exe')
        else:
            # Unrecognised system
            raise ValueError('Your system {} is not supported yet.' % _platform)

        dist_fp = tempfile.NamedTemporaryFile(delete=True, mode='w+t').name

        s = [fastme_exec, "-c", "-{}J".format(datatype), "-i", aln_fp, "-O", dist_fp, "-Q", placement_set_size]
        subprocess.call(s, stdout = nldef, stderr = nldef)


    tbl = open(dist_fp)
    tags = list(re.split("\s+", tbl.readline().rstrip()))[1:]

    for line in tbl.readlines():
        dists = list(re.split("\s+", line.strip()))
        query_name = dists[0]
        obs_dist = dict(zip(tags, map(float, dists[1:])))
        tree = dy.Tree.get_from_string(tree_string, schema='newick', preserve_underscores=True)
        tree.master_edge = next(tree.postorder_edge_iter())

        for l in tree.leaf_node_iter():
            if l.taxon.label not in obs_dist:
                raise ValueError('Taxon {} should be in distances table.'.format(l.taxon.label))
            if obs_dist[l.taxon.label] == 0:
                util.insert(tree, l.edge, query_name, 0, 0)
                tree.write(file=sys.stdout, schema="newick")
                exit(1)

        core = Core(obs_dist, tree)
        core.dp()

        if method_name == "BE":
            alg = BE(core.tree)
        elif method_name == "FM":
            alg = FM(core.tree)
        else:
            alg = OLS(core.tree)
        alg.placement_per_edge(negative_branch)
        output_tree = alg.placement(query_name, selection_name)
        output_tree.write(file=sys.stdout, schema="newick")

    tbl.close()
