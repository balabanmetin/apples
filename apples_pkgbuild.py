from apples.Reference import Reduced_reference
from apples.options_pkgbuild import options_config
import multiprocessing as mp
import json
import time
import logging
from apples.reestimate_backbone import reestimate_backbone

if __name__ == "__main__":
    mp.set_start_method('fork')
    startb = time.time()
    options, args = options_config()
    logging.info("[%s] Options are parsed." % time.strftime("%H:%M:%S"))

    if options.reestimate_backbone:  # reestimate backbone branch lengths
        reestimate_backbone(options)

    start = time.time()
    reference = Reduced_reference(options.ref_fp, options.protein_seqs, options.tree_fp,
                                  options.filt_threshold)
    reference.representatives = [(i.tobytes().decode('utf-8'),j) for i,j in reference.representatives]
    reference.refs = {i: j.tobytes().decode('utf-8') for i,j in reference.refs.items()}
    logging.info(
        "[%s] Reduced reference is computed in %.3f seconds." % (time.strftime("%H:%M:%S"), (time.time() - start)))


    pkg = dict()

    pkg["tree"] = open(options.tree_fp).readline().strip()
    pkg["sequences"] = reference.refs
    pkg["representatives"] = reference.representatives
    pkg["threshold"] = reference.threshold
    with open(options.output_fp, "w") as f:
        json.dump(pkg, f, sort_keys=True, indent=4)

