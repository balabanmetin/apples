import sys
import treeswift as ts
from sys import platform as _platform
import tempfile
from subprocess import Popen
import pkg_resources
import time
import logging


def reestimate_backbone(options):
    assert options.ref_fp
    start = time.time()
    orig_branch_tree = ts.read_tree(options.tree_fp, schema='newick')
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

    bb_fp = tempfile.NamedTemporaryFile(delete=False, mode='w+t')
    fasttree_log = tempfile.NamedTemporaryFile(delete=False, mode='w+t').name
    logging.info("FastTree log file is located here: %s" % fasttree_log)

    s = [fasttree_exec, "-nosupport", "-nome", "-noml", "-log", fasttree_log,
         "-intree", orig_branch_resolved_fp]
    if not options.protein_seqs:
        s.append("-nt")
    with open(options.ref_fp, "r") as rf:
        with Popen(s, stdout=bb_fp, stdin=rf, stderr=sys.stderr):
            options.tree_fp = bb_fp.name
            # tree_string = p.stdout.read().decode('utf-8')
            # print(tree_string)
    logging.info(
        "[%s] Reestimated branch lengths in %.3f seconds." % (time.strftime("%H:%M:%S"), (time.time() - start)))
