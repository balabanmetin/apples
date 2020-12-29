#!/usr/bin/env python3

import re
from apples.PoolQueryWorker import PoolQueryWorker
from apples.fasta2dic import fasta2dic
from apples import util
from apples.Reference import Reduced_reference
from apples.options_run import options_config
import multiprocessing as mp
from apples.jutil import extended_newick, join_jplace
import sys
import json
import treeswift as ts
from sys import platform as _platform
from apples.reestimate_backbone import reestimate_backbone
import time
import logging


if __name__ == "__main__":
    mp.set_start_method('fork')
    startb = time.time()
    options, args = options_config()
    logging.info("[%s] Options are parsed." % time.strftime("%H:%M:%S"))

    if options.reestimate_backbone:  # reestimate backbone branch lengths
        reestimate_backbone(options)

    if options.ref_fp:
        start = time.time()
        reference = Reduced_reference(options.ref_fp, options.protein_seqs, options.tree_fp,
                                      options.filt_threshold)
        reference.set_baseobs(options.base_observation_threshold)
        logging.info(
            "[%s] Reduced reference is computed in %.3f seconds." % (time.strftime("%H:%M:%S"),(time.time() - start)))

        start = time.time()
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
        logging.info(
            "[%s] Query sequences are prepared in %.3f seconds." % (time.strftime("%H:%M:%S"), (time.time() - start)))
    else:
        options.reestimate_backbone = False
        reference = None

        def read_dismat(f):
            tags = list(re.split("\s+", f.readline().rstrip()))[1:]
            for line in f.readlines():
                dists = list(re.split("\s+", line.strip()))
                query_name = dists[0]
                obs_dist = dict(zip(tags, map(float, dists[1:])))
                yield (query_name, None, obs_dist)

        start = time.time()
        f = open(options.dist_fp)
        queries = read_dismat(f)
        logging.info(
            "[%s] Query sequences are prepared in %.3f seconds." % (time.strftime("%H:%M:%S"), (time.time() - start)))

    start = time.time()
    f = open(options.tree_fp)
    orig_tree_string = f.readline()
    f.close()
    logging.info(
        "[%s] Read tree string in %.3f seconds." % (time.strftime("%H:%M:%S"), (time.time() - start)))

    start = time.time()
    first_read_tree = ts.read_tree(options.tree_fp, schema='newick')
    logging.info(
        "[%s] Tree is parsed in %.3f seconds." % (time.strftime("%H:%M:%S"), (time.time() - start)))
    start = time.time()
    util.index_edges(first_read_tree)
    util.set_levels(first_read_tree)

    # create a dictionary where keys are leaf labels and values are
    # pendant edge index for that leaf
    name_to_node_map = {}
    for l in first_read_tree.traverse_postorder(internal=False):
        name_to_node_map[l.label] = l

    extended_newick_string = extended_newick(first_read_tree)
    logging.info(
        "[%s] Tree preprocessing is completed in %.3f seconds." % (time.strftime("%H:%M:%S"), (time.time() - start)))

    startq = time.time()
    queryworker = PoolQueryWorker()
    queryworker.set_class_attributes(reference, options, name_to_node_map)
    if _platform == "win32" or _platform == "win64" or _platform == "msys":
        # if windows, multithreading is not supported until either
        # processes can be forked in windows or apples works with spawn.
        results = list(map(lambda a: queryworker.runquery(a[0], *a[1:]), queries))   # a.k.a starmap
    else:
        pool = mp.Pool(options.num_thread)
        results = pool.starmap(queryworker.runquery, queries)
    logging.info(
        "[%s] Processed all queries in %.3f seconds." % (time.strftime("%H:%M:%S"), (time.time() - startq)))
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
    f.write("\n")
    f.close()
    logging.warning(
        "[%s] APPLES finished in %.3f seconds." % (time.strftime("%H:%M:%S"), (time.time() - startb)))
