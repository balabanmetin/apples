#!/usr/bin/env python3

from optparse import OptionParser
import re
from apples.Core import Core
from apples.fasta2dic import fasta2dic
from apples import util
from apples.Reference import Reduced_reference
import multiprocessing as mp
from apples.jutil import extended_newick, join_jplace
import sys
import json
import treeswift as ts
from apples.criteria import OLS, FM, BE
from sys import platform as _platform


def runquery(query_name, query_seq, obs_dist):
    jplace = dict()
    jplace["placements"] = [{"p": [[0, 0, 1, 0, 0]], "n": [query_name]}]

    if not obs_dist:
        obs_dist = reference.get_obs_dist(query_seq, query_name)
    else:
        tx = 0
        for k, v in sorted(obs_dist.items(), key=lambda kv: kv[1]):
            if v == -1:
                continue
            tx += 1
            if tx > options.base_observation_threshold and v > options.filt_threshold:
                obs_dist[k] = -1

    for l in treecore.tree.traverse_postorder(internal=False):
        if l.label not in obs_dist:
            raise ValueError('Taxon {} should be in distances table.'.format(l.label))
        if obs_dist[l.label] == 0:
            jplace["placements"][0]["p"][0][0] = l.edge_index
            return jplace

    # number of non-infinity distances
    tx = len([v for k, v in obs_dist.items() if v > -1])

    if tx <= 2:
        sys.stderr.write('Taxon {} cannot be placed. At least three non-infinity distances '
                         'should be observed to place a taxon. '
                         'Consequently, this taxon is ignored (no output).\n'.format(query_name))
        jplace["placements"][0]["p"][0][0] = -1
        return jplace

    if -1 not in obs_dist.values():
        tc = treecore
        tc.dp(obs_dist)
    else:
        tc = treecore_frag
        tc.validate_edges(obs_dist)
        tc.dp_frag(obs_dist)

    if options.method_name == "BE":
        alg = BE(tc.tree)
    elif options.method_name == "FM":
        alg = FM(tc.tree)
    else:
        alg = OLS(tc.tree)
    alg.placement_per_edge(options.negative_branch)
    jplace["placements"][0]["p"] = [alg.placement(options.criterion_name, query_name)]
    return jplace


if __name__ == "__main__":
    mp.set_start_method('fork')
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
    parser.add_option("-p", "--protein", dest="protein_seqs", action='store_true', default = False,
                      help="input sequences are protein sequences")
    parser.add_option("-T", "--threads", dest="num_thread", default="0",
                      help="number of cores used in placement. "
                           "0 to use all cores in the running machine", metavar="NUMBER")
    parser.add_option("-f", "--filter", dest="filt_threshold", default="5",
                      help="ignores distances higher than the given threshold. "
                           "Use when long distances have a high bias or variance.", metavar="NUMBER")
    parser.add_option("-b", "--base", dest="base_observation_threshold", default="0",
                      help="minimum number of observations kept for "
                           "each query ignoring the filter threshold.", metavar="NUMBER")
    parser.add_option("-X", "--mask", dest="mask_lowconfidence", action='store_true', default = False,
                      help="masks low confidence characters in the alignments indicated by lowercase characters "
                           "output by softwares like SEPP.")

    (options, args) = parser.parse_args()

    options.num_thread = int(options.num_thread)
    options.filt_threshold = float(options.filt_threshold)
    options.base_observation_threshold = float(options.base_observation_threshold)

    if not options.num_thread:
        options.num_thread = mp.cpu_count()

    if options.ref_fp:
        if options.dist_fp:
            raise ValueError('Input should be either an alignment or a distance matrix, but not both!')

        reference = Reduced_reference(options.ref_fp, options.protein_seqs, options.tree_fp,
                                      options.filt_threshold, options.base_observation_threshold)

        if options.query_fp and options.extended_ref_fp:
            raise ValueError('Input should be either an extended alignment or a query alignment, but not both!')
        if options.query_fp:
            query_dict = fasta2dic(options.query_fp, options.protein_seqs, options.mask_lowconfidence)
        else:  # must be extended reference
            extended_dict = fasta2dic(options.extended_ref_fp, options.protein_seqs, options.mask_lowconfidence)
            query_dict = {key: value for key, value in extended_dict.items() if key not in reference.refs}

        def set_queries(query_dict):
            for query_name, query_seq in query_dict.items():
                yield (query_name, query_seq, None)

        queries = set_queries(query_dict)

    else:
        f = open(options.dist_fp)

        def read_dismat(f):
            tags = list(re.split("\s+", f.readline().rstrip()))[1:]
            for line in f.readlines():
                dists = list(re.split("\s+", line.strip()))
                query_name = dists[0]
                obs_dist = dict(zip(tags, map(float, dists[1:])))
                yield (query_name, None, obs_dist)


        queries = read_dismat(f)
        f.close()

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

    if _platform == "win32" or _platform == "win64" or _platform == "msys":
        # if windows, multithreading is not supported until either
        # processes can be forked in windows or apples works with spawn.
        results = list(map(lambda a: runquery(a[0], *a[1:]), queries))
    else:
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
