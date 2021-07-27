from optparse import OptionParser
from multiprocessing import cpu_count
import logging
from apples.version import __version__

class OptionsBasic(OptionParser):
    def __init__(self, output_filetype):
        super().__init__()
        self.add_option("-t", "--tree", dest="tree_fp",
                        help="path to the reference tree", metavar="FILE")
        self.add_option("-o", "--output", dest="output_fp",
                        help="path for the output %s file" % output_filetype,
                        metavar="FILE")
        self.add_option("-s", "--ref", dest="ref_fp",
                        help="path to the reference alignment file (FASTA), containing reference sequences",
                        metavar="FILE")
        self.add_option("-p", "--protein", dest="protein_seqs", action='store_true', default=False,
                        help="input sequences are protein sequences")
        self.add_option("-T", "--threads", dest="num_thread", type=int, default=0,
                        help="number of cores used in placement. "
                             "0 to use all cores in the running machine", metavar="NUMBER")
        self.add_option("-f", "--filter", dest="filt_threshold", type=float, default=0.2,
                        help="ignores distances higher than the given threshold. "
                             "Improves accuracy when long distances have a high bias or variance.", metavar="NUMBER")
        self.add_option("-D", "--disable-reestimation", dest="disable_reestimation", action='store_true',
                        default=False,
                        help="disables minimum evolution branch length reestimation of the backbone tree. "
                             "This option has no effect if input alignment is not provided.")
        self.add_option("--debug", dest="debug_mode", action='store_true', default=False,
                        help="Enables debug mode.")
        self.add_option("-v", "--version", dest="print_version", action='store_true', default=False,
                        help="print APPLES version number. ")

    def parse(self):
        (options, args) = self.parse_args()
        if options.print_version:
            print("APPLES version " + __version__, flush=True)
            exit(0)

        options.reestimate_backbone = not options.disable_reestimation

        if options.debug_mode:
            root = logging.getLogger()
            root.setLevel(logging.DEBUG)
        if not options.num_thread:
            options.num_thread = cpu_count()
        return options, args
