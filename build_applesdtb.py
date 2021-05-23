#!/usr/bin/env python3
from apples.Reference import ReducedReference
from apples.OptionsBuild import options_config
import multiprocessing as mp
import time
import logging
from apples.prepareTree import prepareTree
import pickle

if __name__ == "__main__":
    mp.set_start_method('fork')
    startb = time.time()
    options, args = options_config()
    logging.info("[%s] Options are parsed." % time.strftime("%H:%M:%S"))

    first_read_tree, name_to_node_map, extended_newick_string = prepareTree(options)
    start = time.time()
    reference = ReducedReference(options.ref_fp, options.protein_seqs, options.tree_fp,
                                 options.filt_threshold, options.num_thread)
    logging.info(
        "[%s] Reduced reference is computed in %.3f seconds." % (time.strftime("%H:%M:%S"), (time.time() - start)))

    with open(options.output_fp, "wb") as f:
        p = pickle.Pickler(f)
        p.dump(first_read_tree)
        p.dump(name_to_node_map)
        p.dump(extended_newick_string)
        p.dump(reference)

    logging.warning(
        "[%s] APPLES database is built in %.3f seconds." % (time.strftime("%H:%M:%S"), (time.time() - startb)))
