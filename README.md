


------------------------------------
Summary
------------------------------------
APPLES stands for `Accurate Phylogenetic Placement with LEast Squares` and addresses the problem of phylogenetic placement of DNA and protein sequences into an already existing reference tree.  APPLES is a command-line tool and it can run on **Linux, Mac OSX, and Windows**.


------------------------------------
Publication
------------------------------------

* Metin Balaban, Shahab Sarmashghi, and Siavash Mirarab. “APPLES: Scalable Distance-Based Phylogenetic Placement with or without Alignments.” Systematic Biology 69, no. 3 (2020): 566–78. [https://doi.org/10.1093/sysbio/syz063](https://doi.org/10.1093/sysbio/syz063)


------------------------------------
Requirements
------------------------------------
1. Python: Version >= 3.0

------------------------------------
Installation on Linux, Mac OSX, or Windows
------------------------------------

Install APPLES using the following command in the command-line:

`pip install apples`


---------------------------------------------
Getting Started with APPLES
---------------------------------------------

For listing all options, run the following command:

`run_apples.py -h`

---------------------------------------------
Input & Output Specification
---------------------------------------------

Input reference (backbone) tree must be in newick format. APPLES can perform placements based on FASTA alignments of nucleotide or amino acid sequences or a distance table.

APPLES input can be either

* An alignment of query sequences to the backbone. In this case, by default, hamming distances will be computed using pairwise comparison and will be then corrected using then JC69 model. If the input is amino acid alignment, APPLES uses scoredist model instead.
* A distance matrix computed using other tools. 

#### Input an alignment 
APPLES require a reference alignment and a query alignment. All species in the backbone tree must have a corresponding sequence in the reference alignment. You can find an example reference alignment and query alignment for ten query sequences under [data/ref.fa](data/ref.fa) and [data/query.fa](data/query.fa) respectively. 

You can run APPLES with the following command on the example input nucleotide alignment dataset:

`run_apples.py -s data/ref.fa -q data/query.fa -t data/backbone.nwk`

If input sequences are amino acid, use the flag `-p`.

#### Input a distance matrix
The format for distance matrix is a tab delimited csv file with column and row headers. Rows should represent query sequences and columns should represent reference sequences. You can find an example distance matrix for ten query sequences under [data/dist.mat](data/dist.mat). 
You can run APPLES on the example distance matrix by running the following command:

`run_apples.py -d data/dist.mat -t data/backbone.nwk`

#### Output
Output is a jplace file containing placement results for all queries. For more information about jplace files, please refer to Matsen et. al. (2012) [https://doi.org/10.1371/journal.pone.0031009](https://doi.org/10.1371/journal.pone.0031009). The output file can be specified using `-o` command. When output file is not specified, the result will be printed to the standard output.

### ! IMPORTANT NOTE !

Backbone tree provided to APPLES has to have its branch lengths estimated using a distance based method such as minimum evolution. This is a requirement for getting good results. We recommend [FastTree2](http://www.microbesonline.org/fasttree/) for re-estimating branch lengths if the backbone tree is estimated using Maximum Likelihood based methods (e.g. RAxML, PASTA). We support re-estimation of branch lengths within APPLES; in fact by default, FastTree is run prior to placement using the following command.

`FastTreeMP -nosupport -nt -nome -noml -log tree.log -intree backbone.nwk < ref.fa > minimum_evo_backbone.nwk`

Then perform placement on the new tree:

`run_apples.py -s ref.fa -q query.fa -t minimum_evo_backbone.nwk`

If you already reestimated the backbone branch lengths, you can skip this step in APPLES using the option "-D" for speedup.


------------------------------------
Detailed tutorial on how to run APPLES on various datasets
------------------------------------

Please refer to the tutorial below for detailed examples of usage in alignment-based and alignment-free settings.

[https://github.com/smirarab/tutorials/blob/master/Skmer-APPLES-tutorial.md](https://github.com/smirarab/tutorials/blob/master/Skmer-APPLES-tutorial.md)

-------------
CHANGELOG
-------------
2.0.10
* APPLES-2 raises a warning when a query sequence is already present in the backbone. 
APPLES-2 places the query anyway (after adding "-query" suffix on the name), ignoring the record in the backbone.
 
2.0.9
* Amino acid distances are now computed using the same BLOSUM45 matrix used in FastTree2.
 
2.0.8
* Pairwise distance overlap requirement parameter can now be changed by the user.