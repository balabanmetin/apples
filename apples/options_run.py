from apples.options_basic import options_basic


def options_config():
    parser = options_basic("jplace")
    parser.add_option("-d", "--distances", dest="dist_fp",
                      help="path to the table of observed distances", metavar="FILE")
    parser.add_option("-x", "--extendedref", dest="extended_ref_fp",
                      help="path to the extened reference alignment file (FASTA), "
                           "containing reference and query sequences",
                      metavar="FILE")
    parser.add_option("-q", "--query", dest="query_fp",
                      help="path to the query alignment file (FASTA), containing query sequences",
                      metavar="FILE")
    parser.add_option("-m", "--method", dest="method_name", default="FM",
                      help="name of the weighted least squares method (OLS, FM, BME, or BE)", metavar="METHOD")
    parser.add_option("-c", "--criterion", dest="criterion_name", default="MLSE",
                      help="name of the placement selection criterion (MLSE, ME, or HYBRID", metavar="CRITERIA")
    parser.add_option("-n", "--negative", dest="negative_branch", action='store_true',
                      help="relaxes positivity constraint on new branch lengths, i.e. allows negative branch lengths")
    parser.add_option("-b", "--base", dest="base_observation_threshold", type=int, default=25,
                      help="minimum number of observations kept for "
                           "each query ignoring the filter threshold.", metavar="NUMBER")
    parser.add_option("-X", "--mask", dest="mask_lowconfidence", action='store_true', default=False,
                      help="masks low confidence characters in the alignments indicated by lowercase characters "
                           "output by softwares like SEPP.")

    (options, args) = parser.parse()
    return options, args
