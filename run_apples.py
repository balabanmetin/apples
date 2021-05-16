#!/usr/bin/env python3
import re
from apples.PoolQueryWorker import PoolQueryWorker
from apples.fasta2dic import fasta2dic
from apples.Reference import ReducedReference
from apples.OptionsRun import options_config
import multiprocessing as mp
from apples.jutil import join_jplace
import sys
import json
from sys import platform as _platform

from apples.prepareTree import prepareTree
import time
import logging
import pickle


if __name__ == "__main__":
    mp.set_start_method('fork')
    startb = time.time()
    options, args = options_config()
    logging.info("[%s] Options are parsed." % time.strftime("%H:%M:%S"))

    if options.database_fp:
        # unpickle tree from database
        start = time.time()
        fdtb = open(options.database_fp, "rb")
        up = pickle.Unpickler(fdtb)
        first_read_tree = up.load()
        name_to_node_map = up.load()
        extended_newick_string = up.load()
        logging.info(
            "[%s] Tree is loaded from APPLES database in %.3f seconds." % (
                time.strftime("%H:%M:%S"), (time.time() - start)))

    if options.tree_fp:
        first_read_tree, name_to_node_map, extended_newick_string = prepareTree(options)

    if options.dist_fp:
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
    else:
        if options.ref_fp:
            start = time.time()
            reference = ReducedReference(options.ref_fp, options.protein_seqs, options.tree_fp,
                                         options.filt_threshold, options.num_thread)
            logging.info(
                "[%s] Reduced reference is computed in %.3f seconds." % (
                    time.strftime("%H:%M:%S"), (time.time() - start)))
        else:  # options.database_fp
            start = time.time()
            reference = up.load()
            fdtb.close()
            logging.info(
                "[%s] Reduced reference is loaded from APPLES database in %.3f seconds." % (
                    time.strftime("%H:%M:%S"), (time.time() - start)))

        reference.set_baseobs(options.base_observation_threshold)
        start = time.time()
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
