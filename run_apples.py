#!/usr/bin/env python3

import re
from apples.Core import Core
from apples.Subtree import Subtree
from apples.fasta2dic import fasta2dic
from apples import util
from apples.Reference import Reduced_reference
from apples.options import options_config
import multiprocessing as mp
from apples.jutil import extended_newick, join_jplace
import sys
import json
import treeswift as ts
from apples.criteria import OLS, FM, BE
from sys import platform as _platform
import tempfile
from subprocess import Popen, PIPE
import pkg_resources


def runquery(query_name, query_seq, obs_dist):
    jplace = dict()
    jplace["placements"] = [{"p": [[0, 0, 1, 0, 0]], "n": [query_name]}]

    if not obs_dist:
        obs_dist = reference.get_obs_dist(query_seq, query_name)
    else:
        def valid_dists(obs_dist, name_to_node_map):
            tx = 0
            for k, v in sorted(obs_dist.items(), key=lambda kv: kv[1]):
                if v < 0 or k not in name_to_node_map:
                    continue
                tx += 1
                if tx > options.base_observation_threshold and v > options.filt_threshold:
                    break
                yield (k, v)
        obs_dist = dict(valid_dists(obs_dist, name_to_node_map))

    if len(obs_dist) <= 2:
        sys.stderr.write('Taxon {} cannot be placed. At least three non-infinity distances '
                         'should be observed to place a taxon. '
                         'Consequently, this taxon is ignored (no output).\n'.format(query_name))
        jplace["placements"][0]["p"][0][0] = -1
        return jplace

    for k, v in obs_dist.items():
        if v == 0 and k != query_name:
            jplace["placements"][0]["p"][0][0] = name_to_node_map[k].edge_index
            return jplace

    tc = treecore_frag
    subtree = Subtree(obs_dist, name_to_node_map)
    tc.dp_frag(subtree)

    if options.method_name == "BE":
        alg = BE(tc.tree)
    elif options.method_name == "FM":
        alg = FM(tc.tree)
    else:
        alg = OLS(tc.tree)
    alg.placement_per_edge(options.negative_branch)
    jplace["placements"][0]["p"] = [alg.placement(options.criterion_name)]
    subtree.unroll_changes()
    return jplace


if __name__ == "__main__":
    mp.set_start_method('fork')

    parser = options_config()
    (options, args) = parser.parse_args()

    options.num_thread = int(options.num_thread)
    options.filt_threshold = float(options.filt_threshold)
    options.base_observation_threshold = float(options.base_observation_threshold)
    options.reestimate_backbone = not options.disable_reestimation

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
        options.reestimate_backbone = False

        def read_dismat(f):
            tags = list(re.split("\s+", f.readline().rstrip()))[1:]
            for line in f.readlines():
                dists = list(re.split("\s+", line.strip()))
                query_name = dists[0]
                obs_dist = dict(zip(tags, map(float, dists[1:])))
                yield (query_name, None, obs_dist)

        f = open(options.dist_fp)
        queries = read_dismat(f)

    f = open(options.tree_fp)
    orig_tree_string = f.readline()
    f.close()

    if options.reestimate_backbone:  # reestimate backbone branch lengths
        assert options.ref_fp
        orig_branch_tree = ts.read_tree(orig_tree_string, schema='newick')
        orig_branch_tree.resolve_polytomies()
        orig_branch_resolved_fp = tempfile.NamedTemporaryFile(delete=True, mode='w+t').name
        orig_branch_tree.write_tree_newick(orig_branch_resolved_fp)

        if _platform == "darwin":
            fasttree_exec = pkg_resources.resource_filename('apples', "tools/FastTree-darwin")
        elif _platform == "linux" or _platform == "linux2":
            fasttree_exec = pkg_resources.resource_filename('apples', "tools/FastTree-linux")
        elif _platform == "win32" or _platform == "win64" or _platform == "msys":
            fasttree_exec = pkg_resources.resource_filename('apples', "tools/FastTree.exe")
        else:
            # Unrecognised system
            raise ValueError('Your system {} is not supported yet.' % _platform)

        bb_fp = tempfile.NamedTemporaryFile(delete=True, mode='w+t')
        fasttree_log = tempfile.NamedTemporaryFile(delete=True, mode='w+t').name

        s = [fasttree_exec, "-nosupport", "-nome", "-noml", "-log", fasttree_log,
             "-intree", orig_branch_resolved_fp]
        if not options.protein_seqs:
            s.append("-nt")
        with open(options.ref_fp, "r") as rf:
            with Popen(s, stdout=PIPE, stdin=rf, stderr=sys.stderr) as p:
                tree_string = p.stdout.read().decode('utf-8')
                print(tree_string)
    else:
        tree_string = orig_tree_string

    first_read_tree = ts.read_tree(tree_string, schema='newick')
    util.index_edges(first_read_tree)
    util.set_levels(first_read_tree)

    # create a dictionary where keys are leaf labels and values are
    # pendant edge index for that leaf
    name_to_node_map = {}
    for l in first_read_tree.traverse_postorder(internal=False):
        name_to_node_map[l.label] = l

    extended_newick_string = extended_newick(first_read_tree)
    util.index_edges(first_read_tree)
    util.set_levels(first_read_tree)
    treecore_frag = Core(first_read_tree)

    if _platform == "win32" or _platform == "win64" or _platform == "msys":
        # if windows, multithreading is not supported until either
        # processes can be forked in windows or apples works with spawn.
        results = list(map(lambda a: runquery(a[0], *a[1:]), queries))   # a.k.a starmap
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
