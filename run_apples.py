#!/usr/bin/env python3

from optparse import OptionParser
import re
from apples.Core import Core
from apples.criteria import OLS, FM, BE
from apples import readfq, distance, util
import multiprocessing as mp
from apples.jutil import extended_newick, join_jplace
import sys
import json
import numpy as np
import treeswift as ts




def runquery(query_name, query_seq, obs_dist):
    jplace=dict()
    jplace["placements"] = [{"p":[[0,0,1,0,0]] ,"n":[query_name]}]
    if not obs_dist:
        obs_dist = {query_name: 0}
        for tagr,seqr in zip(reftags, refseqs):
            obs_dist[tagr] = distance.jc69(query_seq, seqr)

    for l in treecore.tree.traverse_postorder(internal=False):
        if l.label not in obs_dist:
            raise ValueError('Taxon {} should be in distances table.'.format(l.label))
        if obs_dist[l.label] == 0:
            jplace["placements"][0]["p"][0][0] = l.edge_index
            return jplace

    if -1 not in obs_dist.values():
        tc = treecore
        tc.dp(obs_dist)
    else:
        tc = treecore_frag
        tc.validate_edges(obs_dist)
        tc.dp_frag(obs_dist)

    if method_name == "BE":
        alg = BE(tc.tree)
    elif method_name == "FM":
        alg = FM(tc.tree)
    else:
        alg = OLS(tc.tree)
    alg.placement_per_edge(negative_branch)
    jplace["placements"][0]["p"] = [alg.placement(criterion_name)]
    return jplace


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
                      help="path to the extened reference alignment file (FASTA), containing reference and query sequences",
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
    #parser.add_option("-p", "--protein", dest="protein_seqs", action='store_true',
    #                  help="input sequences are protein sequences")
    parser.add_option("-T", "--threads", dest="num_thread", default="0",
                      help="number of cores used in placement. 0 to use all cores in the running machine", metavar="NUMBER")

    (options, args) = parser.parse_args()
    tree_fp = options.tree_fp
    dist_fp = options.dist_fp
    ref_fp = options.ref_fp
    output_fp = options.output_fp
    query_fp = options.query_fp
    extended_ref_fp = options.extended_ref_fp
    method_name = options.method_name
    criterion_name = options.criterion_name
    negative_branch = options.negative_branch
    #protein_seqs = options.protein_seqs
    num_thread = int(options.num_thread)
    if not num_thread:
        num_thread = mp.cpu_count()

    if ref_fp:
        if dist_fp:
            raise ValueError('Input should be either an alignment or a distance matrix, but not both!')

        reftags = []
        refseqs = []
        querytags = []
        queryseqs = []
        num_query = 0

        f = open(ref_fp)
        for name, seq, qual in readfq.readfq(f):
            reftags.append(name)
            refseqs.append(np.frombuffer(seq.upper().encode(), dtype='S1'))
        f.close()

        if query_fp and extended_ref_fp:
            raise ValueError('Input should be either an extended alignment or a query alignment, but not both!')
        if query_fp:
            f = open(query_fp)
            for name, seq, qual in readfq.readfq(f):
                querytags.append(name)
                queryseqs.append(np.frombuffer(seq.upper().encode(), dtype='S1'))
            f.close()
            num_query = len(querytags)
        else:
            f = open(extended_ref_fp)
            setreftags = set(reftags)
            translation = str.maketrans('abcdefghijklmnopqrstuvwxyz', '-'*26)
            for name, seq, qual in readfq.readfq(f):
                if name not in setreftags:
                    querytags.append(name)
                    queryseqs.append(np.frombuffer(seq.translate(translation).encode(), dtype='S1'))
            num_query = len(querytags)
            f.close()

        queries = zip(querytags, queryseqs, num_query*[None])

    else:
        f = open(dist_fp)
        def read_dismat(f):
            tags = list(re.split("\s+", f.readline().rstrip()))[1:]
            for line in f.readlines():
                dists = list(re.split("\s+", line.strip()))
                query_name = dists[0]
                obs_dist = dict(zip(tags, map(float, dists[1:])))
                yield (query_name, None, obs_dist)
        queries = read_dismat(f)

    f = open(tree_fp)
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

    pool = mp.Pool(num_thread)
    results = pool.starmap(runquery, queries)

    result = join_jplace(results)
    result["tree"] = extended_newick_string
    result["metadata"] = {"invocation":" ".join(sys.argv)}
    result["fields"] = ["edge_num", "likelihood", "like_weight_ratio", "distal_length", "pendant_length"]
    result["version"] = 3

    if output_fp:
        f = open(output_fp, "w")
    else:
        f = sys.stdout
    f.write(json.dumps(result, sort_keys=True, indent=4))
    f.close()
