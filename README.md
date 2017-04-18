# PBSuite

This package was originally developed by Adam English to utilize PacBio long reads for gap filling and analysis of structural variants.

# jelly2

In this copy of the package, I am stripping out the functions for analysis of structural variants and re-factoring the gap filling module to reflect the latest developments in PacBio technology and Python conventions.

Thus the focus of this work is to improve the gap filling and scaffolding capabilities of the jelly module within PBSuite. The goal is to yield an updated and improved version of jelly for finishing de novo genome assemblies.

# Dependencies

* Python (v2.7)
* BLASR (https://github.com/PacificBiosciences/blasr)
* NetworkX (https://github.com/networkx/networkx)
* Pysam (https://github.com/pysam-developers/pysam)
* Minimap (https://github.com/lh3/minimap)
* Miniasm (https://github.com/lh3/miniasm)
* Racon (https://github.com/isovic/racon)

# Design Philosophy

BLASR is used to conduct the initial alignment of PacBio reads to the draft genome sequences, finding reads that span gaps. A NetworkX graph is constructed to store the gap information, as well as add new gaps detected from the alignment that were not in the original assembly. Gap sequences are assembled with Minimap (overlap), Miniasm (layout), and Racon (consensus). Assembled gap sequences are then inserted into the original assembly and the new sequences are written into a Fasta file.
