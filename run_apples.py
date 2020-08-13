#!/usr/bin/env python3

from optparse import OptionParser
import re
from apples.Core import Core
from apples.readfq import readfq
from apples import util
from apples.runquery import runquery
from apples.Reference import *
import multiprocessing as mp
from apples.jutil import extended_newick, join_jplace
import sys
import json
import numpy as np
import treeswift as ts



if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-t", "--tree", dest="tree_fp",
                      help="path to the reference tree", metavar="FILE")
    parser.add_option("-d", "--distances", dest="dist_fp",
                      help="path to the table of observed distances", metavar="FILE")
    parser.add_option("-o", "--output", dest="output_fp",
                      help="path for the output jplace file",
                      metavar="FILE")
    parser.add_option("-s", "--ref", dest="ref_fp",
                      help="path to the reference alignment file (FASTA), containing reference sequences",
                      metavar="FILE")
    parser.add_option("-x", "--extendedref", dest="extended_ref_fp",
                      help="path to the extened reference alignment file (FASTA), "
                           "containing reference and query sequences",
                      metavar="FILE")
    parser.add_option("-q", "--query", dest="query_fp",
                      help="path to the query alignment file (FASTA), containing query sequences",
                      metavar="FILE")
    parser.add_option("-m", "--method", dest="method_name", default="FM",
                      help="name of the weighted least squares method (OLS, FM, or BE)", metavar="METHOD")
    parser.add_option("-c", "--criterion", dest="criterion_name", default="MLSE",
                      help="name of the placement selection criterion (MLSE, ME, or HYBRID", metavar="CRITERIA")
    parser.add_option("-n", "--negative", dest="negative_branch", action='store_true',
                      help="relaxes positivity constraint on new branch lengths, i.e. allows negative branch lengths")
    # parser.add_option("-p", "--protein", dest="protein_seqs", action='store_true',
    #                  help="input sequences are protein sequences")
    parser.add_option("-T", "--threads", dest="num_thread", default="0",
                      help="number of cores used in placement. "
                           "0 to use all cores in the running machine", metavar="NUMBER")
    parser.add_option("-f", "--filter", dest="filt_threshold", default="5",
                      help="ignores distances higher than the given threshold. "
                           "Use when long distances have a high bias or variance.", metavar="NUMBER")
    parser.add_option("-b", "--base", dest="base_observation_threshold", default="0",
                      help="minimum number of observations kept for "
                           "each query ignoring the filter threshold.", metavar="NUMBER")

    (options, args) = parser.parse_args()

    # protein_seqs = options.protein_seqs
    options.num_thread = int(options.num_thread)
    options.filt_threshold = float(options.filt_threshold)
    options.base_observation_threshold = float(options.base_observation_threshold)

    if not options.num_thread:
        options.num_thread = mp.cpu_count()

    if options.ref_fp:
        if options.dist_fp:
            raise ValueError('Input should be either an alignment or a distance matrix, but not both!')

        querytags = []
        queryseqs = []
        num_query = 0

        reference = Reduced_reference(options.ref_fp, options.tree_fp,
                                      options.filt_threshold, options.base_observation_threshold)

        if options.query_fp and options.extended_ref_fp:
            raise ValueError('Input should be either an extended alignment or a query alignment, but not both!')
        if options.query_fp:
            f = open(options.query_fp)
            for name, seq, qual in readfq(f):
                querytags.append(name)
                queryseqs.append(np.frombuffer(seq.upper().encode(), dtype='S1'))
            f.close()
            num_query = len(querytags)
        else:
            f = open(options.extended_ref_fp)
            setreftags = set(reftags)
            translation = str.maketrans('abcdefghijklmnopqrstuvwxyz', '-' * 26)
            for name, seq, qual in readfq(f):
                if name not in setreftags:
                    querytags.append(name)
                    queryseqs.append(np.frombuffer(seq.translate(translation).encode(), dtype='S1'))
            num_query = len(querytags)
            f.close()

        def set_queries(querytags, queryseqs):
            for querytags, queryseqs in zip(querytags, queryseqs):
                yield (treecore, treecore_frag, options, querytags, queryseqs, None, reference)


        queries = set_queries(querytags, queryseqs)

    else:
        f = open(options.dist_fp)


        def read_dismat(f):
            tags = list(re.split("\s+", f.readline().rstrip()))[1:]
            for line in f.readlines():
                dists = list(re.split("\s+", line.strip()))
                query_name = dists[0]
                obs_dist = dict(zip(tags, map(float, dists[1:])))
                yield (treecore, treecore_frag, options, query_name, None, obs_dist, None)


        queries = read_dismat(f)

    f = open(options.tree_fp)
    tree_string = f.readline()
    f.close()

    first_read_tree = ts.read_tree(tree_string, schema='newick')
    util.index_edges(first_read_tree)
    extended_newick_string = extended_newick(first_read_tree)
    treecore = Core(first_read_tree)
    treecore.init()
    second_read_tree = ts.read_tree(tree_string, schema='newick')
    util.index_edges(second_read_tree)
    treecore_frag = Core(second_read_tree)

    pool = mp.Pool(options.num_thread)
    results = pool.starmap(runquery, queries)

    result = join_jplace(results)
    result["tree"] = extended_newick_string
    result["metadata"] = {"invocation": " ".join(sys.argv)}
    result["fields"] = ["edge_num", "likelihood", "like_weight_ratio", "distal_length", "pendant_length"]
    result["version"] = 3

    if options.output_fp:
        f = open(options.output_fp, "w")
    else:
        f = sys.stdout
    f.write(json.dumps(result, sort_keys=True, indent=4))
    f.close()
