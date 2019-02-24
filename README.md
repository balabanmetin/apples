------------------------------------
Summary
------------------------------------
APPLES stands for `Accurate Phylogenetic Placement with LEast Squares` and addresses the problem of phylogenetic placement of DNA and protein sequences into an already existing reference tree. APPLES is a command line tool and it can run on **Linux, Mac OSX, and Windows**.

Requirements
-------------------
1. Python: Version >= 3.0
2. Dendropy: Version >= 4.0.0

------------------------------------
Installation on Linux, Mac OSX, or Windows
------------------------------------

Install Dendropy using the following command in the command-line:

`pip3 install --user dendropy`

Then, clone the repository using the following command:

`git clone https://github.com/balabanmetin/apples.git`

Finally, change the working directory:

`cd apples`

---------------------------------------------
Getting Started with APPLES
---------------------------------------------

For listing all options, run the following command:

`python3 apples.py -h`

---------------------------------------------
Input & Output Specification
---------------------------------------------

Input reference (backbone) tree must be in newick format. APPLES can perform placements based on PHYLIP alignments or a distance table.
#### Input an alignment 
The first `k` entry of the input alignment must be the `k` query sequences that will be placed. All species in the backbone tree must have a corresponding sequence in the alignment. You can find an example distance matrix for ten query sequences under [data/aln.phy](data/aln.phy).

You can run APPLES with the following command on the example input alignment dataset:

`python3 apples.py -a data/aln.phy -t data/backbone.nwk -P 10`

#### Input a distance matrix
The format for distance matrix is a tab delimited csv file with column and row headers. You can find an example distance matrix for ten query sequences under [data/dist.mat](data/dist.mat).
You can run APPLES on the example distance matrix by running the following command:

`python3 apples.py -d data/dist.mat -t data/backbone.nwk -P 10`

#### Output
Output is multiple lines of newick trees, i-th line is corresponding to the placement of the i-th entry in alignment/distance matrix.

---------------------------------------------
Acknowledgement
---------------------------------------------

APPLES uses my personal fork of FastME 2.0 software, developed by Lefort V., Desper R., Gascuel O.. I do not own or maintain the FastME 2.0 software, all credit goes to actual authors.
The original software can be found via the following link:

http://www.atgc-montpellier.fr/fastme/

