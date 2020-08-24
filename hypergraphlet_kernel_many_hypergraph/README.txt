Hypergraph Kernel Framework version 1.0 (Adapted from Vladimir Vacic graphlet kernel code)

Jose Lugo-Martinez, jlugomar@andrew.cmu.edu
Computational Biology Department
School of Computer Science
Carnegie Mellon University

This is beta version 1.0 of the hypergraph kernels software for vertex classification, 
edge classification and link prediction. This version provides a framework for users to 
run a variety of hypergraph-based kernel methods on (sparse) vertex- and egde-labeled 
hypergraphs. Kernels implemented include: cumulative random walk, standard random walk, 
standard hypergraphlet kernel, label substitutions hypergraphlet kernel, edge indels 
hypergraphlet kernel and edit distance hypergraphlet kernel.

Please direct all comments and bug reports of this version to jlugomar@andrew.cmu.edu

A detailed description of each method and discussion of parameters can be found in:

Lugo-Martinez et al. "Classification in biological networks with hypergraphlet kernels", 
2017, arXiv:1703.04823.

Please cite the arXiv paper if you use or modify this software.

------------------------------------------------------------------------
COMPILATION
------------------------------------------------------------------------

To compile the hypergraph kernel framework, type "make" on the command prompt. 
Program file "run_hypergraphkernel" will be generated. Running the binary with -h switch 
will list all the command line options.


------------------------------------------------------------------------
PROGRAM OPTIONS
------------------------------------------------------------------------

Usage: run_hyperkernel -p FILE -n FILE -g HG_PATH -l L_PATH -e E_PATH -t TYPE -z TASK -[k|s] OUTPUT [...]
Options:

  -h         Displays this message.

  -t TYPE    Kernel type (0-Cumulative Random Walk, 1-Standard Random Walk, 2-Standard Hypergraphlet, 3-Label Substitutions Hypergraphlet, 4-Hyperedge Indels Hypergraphlet, 5-Edit Distance Hypergraphlet)
             Defaults to standard graphlet.

  -z TASK    Task type (0-Vertex classification, 1-Edge classification, 2-Link prediction)
             Defaults to vertex classification.

  -p FILE    List of positive (vertices) examples.
  -n FILE    List of negative (vertices) examples.
  -g HG_PATH Path to hypergraph files.
  -l L_PATH  Path to vertex labels file for input hypergraph files.
  -e E_PATH  Path to edge labels file for input hypergraph files.

  -N         Normalize the kernel matrix.
             Defaults to false.

  -k KERNEL  Output file for the kernel matrix in standard output.
  -s SPARSE  Output file for the sparse attribute matrix (SVML).
             Defaults to KERNEL.

  -I STEPS   Number of steps. (Needed for Random Walk Hypergraph kernels)
             Defaults to 100,000 steps.

  -R RESTART Restart probability. (Needed for Random Walk Hypergraph kernels)
             Defaults to 0.3

  -S SIMVMAT Similarity matrix for weighting each possible vertex label substitution. (Needed for Label Substitutions and Edit Distance Hypergraph kernels)
             Defaults to uniform weights (i.e. all label mismatches are equally weighted to 1).

  -P SIMEMAT Similarity matrix for weighting each possible hyperedge label substitution. (Needed for Label Substitutions, Edge Indels and Edit Distance Hypergraph kernels)
             Defaults to uniform weights (i.e. all label mismatches are equally weighted to 1).

  -K VLABMIS Fraction of nodes in the n-hypergraphlet allowed to have vertex label mismatches. (Needed for Label Substitutions and Edit Distance Hypergraph kernels)
             For example, K=0.34 allows 1-vertex label mismatch for 3- and 4-hypergraphlets, whereas K=0.26 allows 1-vertex label mismatch for 4-hypergraphlets only.
             Defaults to 0.0 (i.e. no vertex label mismatches allowed).

  -M ELABMIS Total number of hyperedge label mismacthes allowed between graphlets. (Needed for Label Substitutions, Edge Indels and Edit Distance Hypergraph kernels)
             Defaults to 0.

  -m EDGMIS  Total number of hyperedge insertions or deletions allowed between hypergraphlets. (Needed for Edge Indels and Edit Distance Hypergraph Kernels)
             Defaults to 1.

  -V ALPHA   Vertex labels alphabet over problem statement is defined (i.e. all possible labels for a vertex). (Needed for Label Substitutions and Edit Distance Hypergraph Kernels)

  -E ALPHA   Hyperedge labels alphabet over problem statement is defined (i.e. all possible labels for an hyperedge). (Needed for Label Substitutions, Edge Indels and Edit Distance Hypergraph kernels)

  -c LABELS  Output file for each example class label.

  -v         Verbose (prints progress messages).


------------------------------------------------------------------------
EXAMPLE
------------------------------------------------------------------------

The examples subdirectory contains two subfolders. Each subfolder contains 
files with positives and negatives vertices within a vertex- and 
edge-labeled hypergraph. For example subfolder example1 contains two files 
example1.pos and example1.neg, respectively. These files have two fields which 
indicate the hypergraph file and corresponding vertex and hyperedges files, 
as well as vertex of interest within each input hypergraph file.

Preparing a kernel matrix from hypergraph files is one step process:

run_hypergraphkernel is used to generate the user-specified hypergraph-based 
kernel method. There are several output options (see config.h), out of which 
SVM^Light format is the easiest to use, because it can readily be used with 
SVM^Light (see http://svmlight.joachims.org). In summary, this is a space-separated
file with the first field equal to 1 (for positives) or -1 (for negatives),
and all other non-zero entries are feature_key:feature_value pairs. 

A sample "run_example.sh" shell script is provided as an illustrative example 
of the process and parameters required to run each kernel method.


------------------------------------------------------------------------
COMMENT REGARDING LARGE DATASETS
------------------------------------------------------------------------

Using the SVM^Light format directly on very large datasets is not 
efficient, because the dot-product between feature vectors has to be
unnecessarily recomputed over and over again. There is an option to 
precompute the kernel matrix and save it (for different output options, 
see config.h). However, using it with SVM^Light is slightly more complicated, 
because the reader for this file format has to be custom coded and according 
to the SVM^Light copyright agreement we are not allowed to distribute modified
SVM^Light code.

