import sys
import treeswift as ts
from sys import platform as _platform
import tempfile
from subprocess import Popen, PIPE
import pkg_resources
import time
import logging


def print_ident(tree):
    """
    Function to print the length of each node's postorder traversal.
    """
    for i in tree.root.children:
        print(
            len(list(i.traverse_postorder())),
        )
    print('')


def reestimate_backbone(options):
    """
    This code snippet is a Python function that reestimates the branch lengths of a phylogenetic tree.
    It reads in a tree from a file, resolves polytomies, checks for branch length completeness,
    runs a tree estimation program (FastTree) based on the input sequences and fixed tree topology,
    and then updates the branch lengths of the input tree based on the re-estimated tree.
    Finally, it logs the time taken for the reestimation process.
    Internally, FastTree computes branch lengths using distance-based branch lenghts, either JC69 or BLOSUM45.
    """
    assert options.ref_fp
    start = time.time()
    orig_branch_tree = ts.read_tree(options.tree_fp, schema='newick')
    if len(orig_branch_tree.root.children) > 2:  # 3
        rooted = False
    else:
        rooted = True
    orig_branch_tree.suppress_unifurcations()
    if len(orig_branch_tree.root.children) > 3:
        # polytomy at the root
        orig_branch_tree.resolve_polytomies()
    else:
        # root node is ok, resolve the other nodes
        for i in orig_branch_tree.root.children:
            i.resolve_polytomies()
    all_branches_have_length = True
    for n in orig_branch_tree.traverse_postorder(internal=True, leaves=True):
        if not n.is_root() and n.edge_length is None:
            all_branches_have_length = False
            break

    if rooted and all_branches_have_length:
        left, right = orig_branch_tree.root.children
        if left.children:
            thetwo = [next(c.traverse_postorder(internal=False)) for c in left.children]
            theone = [next(right.traverse_postorder(internal=False))]
            lengthtwoside = left.edge_length
            lengthoneside = right.edge_length
        else:
            thetwo = [next(c.traverse_postorder(internal=False)) for c in right.children]
            theone = [next(left.traverse_postorder(internal=False))]
            lengthtwoside = right.edge_length
            lengthoneside = left.edge_length

    orig_branch_resolved_fp = tempfile.NamedTemporaryFile(delete=True, mode='w+t').name
    orig_branch_tree.write_tree_newick(orig_branch_resolved_fp)

    if _platform == 'darwin':
        fasttree_exec = pkg_resources.resource_filename('apples', 'tools/FastTree-darwin')
    elif _platform == 'linux' or _platform == 'linux2':
        fasttree_exec = pkg_resources.resource_filename('apples', 'tools/FastTree-linux')
    elif _platform == 'win32' or _platform == 'win64' or _platform == 'msys':
        fasttree_exec = pkg_resources.resource_filename('apples', 'tools/FastTree.exe')
    else:
        # Unrecognised system
        raise ValueError('Your system {} is not supported yet.' % _platform)

    bb_fp = tempfile.NamedTemporaryFile(delete=False, mode='w+t')
    fasttree_log = tempfile.NamedTemporaryFile(delete=False, mode='w+t').name
    logging.info('FastTree log file is located here: %s' % fasttree_log)

    s = [fasttree_exec, '-nosupport', '-nome', '-noml', '-log', fasttree_log, '-intree', orig_branch_resolved_fp]
    if not options.protein_seqs:
        s.append('-nt')

    # match the rooting of the fasttree output to the original input tree.
    with open(options.ref_fp, 'r') as rf:
        with Popen(s, stdout=PIPE, stdin=rf, stderr=sys.stderr) as p:
            # options.tree_fp = bb_fp.name
            tree_string = p.stdout.read().decode('utf-8')
            if rooted and all_branches_have_length:
                ft = ts.read_tree_newick(tree_string)
                for n in ft.traverse_postorder(internal=False):
                    if n.label == theone[0].label:
                        theone_inft = n
                        break
                ft.reroot(theone_inft)
                mrca = ft.mrca([n.label for n in thetwo])
                mrca_edge_length = mrca.edge_length
                ft.reroot(mrca, length=mrca_edge_length / 2)
                if lengthtwoside + lengthoneside > 0:
                    for i in range(2):
                        if ft.root.children[i] == mrca:
                            ft.root.children[i].edge_length = (
                                mrca_edge_length * lengthtwoside / (lengthtwoside + lengthoneside)
                            )
                            ft.root.children[1 - i].edge_length = (
                                mrca_edge_length * lengthoneside / (lengthtwoside + lengthoneside)
                            )
                ft.is_rooted = False
                tree_string = str(ft)
            with open(bb_fp.name, 'w') as ntree:
                ntree.write(tree_string.strip())
                ntree.write('\n')
                options.tree_fp = bb_fp.name
    logging.info(
        '[%s] Reestimated branch lengths in %.3f seconds.' % (time.strftime('%H:%M:%S'), (time.time() - start))
    )
