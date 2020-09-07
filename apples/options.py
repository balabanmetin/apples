from optparse import OptionParser


def options_config():
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
                      help="path to the extened reference alignment file (FASTA), "
                           "containing reference and query sequences",
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
    parser.add_option("-p", "--protein", dest="protein_seqs", action='store_true', default=False,
                      help="input sequences are protein sequences")
    parser.add_option("-T", "--threads", dest="num_thread", default="0",
                      help="number of cores used in placement. "
                           "0 to use all cores in the running machine", metavar="NUMBER")
    parser.add_option("-f", "--filter", dest="filt_threshold", default="5",
                      help="ignores distances higher than the given threshold. "
                           "Use when long distances have a high bias or variance.", metavar="NUMBER")
    parser.add_option("-b", "--base", dest="base_observation_threshold", default="0",
                      help="minimum number of observations kept for "
                           "each query ignoring the filter threshold.", metavar="NUMBER")
    parser.add_option("-X", "--mask", dest="mask_lowconfidence", action='store_true', default=False,
                      help="masks low confidence characters in the alignments indicated by lowercase characters "
                           "output by softwares like SEPP.")
    parser.add_option("-D", "--disable-reestimation", dest="disable_reestimation", action='store_true', default=False,
                      help="disables minimum evolution branch length reestimation of the backbone tree. "
                           "This option has no effect if input alignment is not provided.")
    return parser
