from apples import util
from apples.jutil import extended_newick
from apples.reestimateBackbone import reestimate_backbone
import treeswift as ts
import time
import logging


def prepareTree(options):
    if options.reestimate_backbone:  # reestimate backbone branch lengths
        reestimate_backbone(options)

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
    return first_read_tree, name_to_node_map, extended_newick_string