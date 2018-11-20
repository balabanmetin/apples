------------------------------------
Summary
------------------------------------
APPLES stands for `Accurate Phylogenetic Placement with LEast Squares` and addresses the problem of phylogenetic placement of DNA and protein sequences into an already existing reference tree.

Requirements
-------------------
1. Python: Version >= 3.0
2. Dendropy: Version >= 4.0.0

---------------------------------------------
Getting Started with APPLES
---------------------------------------------
You can get started running APPLES with the following command on the example dataset:

`python apples.py -a data/aln.phy -t backbone.nwk -P 10`

For more options, run the following command:

`python apples.py -h`

---------------------------------------------
Input & Output Specification
---------------------------------------------

Input reference (backbone) tree must be in newick format. APPLES can perform placements based on PHYLIP alignments xor a distance table.
#### Input an alignment 
The first `k` entry of the input alignment must be the `k` query sequences that will be placed. All species in the backbone tree must have a corresponding sequence in the alignment.

#### Input a distance matrix
The format for distance matrix is a tab delimited csv file with column and row headers. You can find an example distance matrix for ten query sequences under `data/dist.mat`.
You can run APPLES on the example distance matrix by running the following command:

`python apples.py -d data/dist.mat -t data/backbone.nwk -P 10`

#### Output
Output is multiple lines of newick trees, i-th line is corresponding to the placement of the i-th entry in alignment/distance matrix.
