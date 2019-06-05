


------------------------------------
Summary
------------------------------------
APPLES stands for `Accurate Phylogenetic Placement with LEast Squares` and addresses the problem of phylogenetic placement of DNA and protein sequences into an already existing reference tree.  APPLES is a command-line tool and it can run on **Linux, Mac OSX, and Windows**.


------------------------------------
Preprint
------------------------------------

[APPLES: Distance-based Phylogenetic Placement for Scalable and Assembly-free Sample Identification,
Metin Balaban, Shahab Sarmashghi, Siavash Mirarab,
bioRxiv 475566; doi: https://doi.org/10.1101/475566](https://doi.org/10.1101/475566)

------------------------------------
Requirements
------------------------------------
1. Python: Version >= 3.0
2. Dendropy: Version >= 4.0.0

------------------------------------
Installation on Linux, Mac OSX, or Windows
------------------------------------

Install Dendropy using the following command in the command-line:

`pip3 install --user dendropy`

Then, clone the repository using the following command:

`git clone https://github.com/balabanmetin/apples.git`

Once the repository is downloaded, make [apples.py](run_apples.py) executable.


---------------------------------------------
Getting Started with APPLES
---------------------------------------------

For listing all options, run the following command:

`apples.py -h`

---------------------------------------------
Input & Output Specification
---------------------------------------------

Input reference (backbone) tree must be in newick format. APPLES can perform placements based on FASTA alignments of __nucleotide sequences__ (functionality for amino acid sequences will be added soon) or a distance table.
#### Input an alignment 
APPLES require a reference alignment and a query alignment. All species in the backbone tree must have a corresponding sequence in the reference alignment. You can find an example reference alignment and query alignment for ten query sequences under [data/ref.fa](data/ref.fa) and [data/query.fa](data/query.fa) respectively.

You can run APPLES with the following command on the example input alignment dataset:

`apples.py -s data/ref.fa -q data/query.fa -t data/backbone.nwk`

#### Input a distance matrix
The format for distance matrix is a tab delimited csv file with column and row headers. Rows should represent query sequences and columns should represent reference sequences. You can find an example distance matrix for ten query sequences under [data/dist.mat](data/dist.mat). 
You can run APPLES on the example distance matrix by running the following command:

`apples.py -d data/dist.mat -t data/backbone.nwk`

#### Output
Output is a jplace file containing placement results for all queries. For more information about jplace files, please refer to Matsen et. al. (2012) [https://doi.org/10.1371/journal.pone.0031009](https://doi.org/10.1371/journal.pone.0031009). The output file can be specified using `-o` command. When output file is not specified, the result will be printed to the standard output.

