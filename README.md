# Jelly2

This project started out as a hack of [PBJelly](https://github.com/dbrowneup/PBSuite) in an attempt to improve the efficiency of the code. I ended up rebuilding the core pipeline, focusing it down into four stages: Setup, Support, Assembly, and Placement. The stages are managed by a single driver script, `Jelly2.py`, making this program fairly straightforward to use. To see the usage, run:

```
$ python Jelly2.py --help
```

Make sure that the `src` directory is in your `PYTHONPATH` variable, so the driver script can find the modules. Though this repository was created on May 8, 2017, it contains the git tracking history of the above-linked PBJelly repository in which this project began.

# Dependencies

* Python (v2.7)
* [BLASR](https://github.com/PacificBiosciences/blasr)
* [pysam](https://github.com/pysam-developers/pysam)
* [pyfaidx](https://github.com/mdshw5/pyfaidx)
* [Minimap](https://github.com/lh3/minimap)
* [Miniasm](https://github.com/lh3/miniasm)
* [Racon](https://github.com/isovic/racon)

# Design Philosophy

To find PacBio reads that span the gaps predicted in the scaffolds, gap-flanking sequences are extracted from the assembly. Separate Fasta files are created for the left flanks and the right flanks. The PacBio reads are then aligned against each flank file with BLASR. Supporting reads are those that align to each flank, spanning the predicted gap size within a certain margin of error. If there are enough support reads for a gap, the reads will be assembled with Minimap and Miniasm. A consensus will be determined with Minimap and Racon. If the assembled gap sequence passes quality control, it will be placed into the scaffold, filling the gap.

# Contact

Any questions can be addressed to: dbrowne@tamu.edu
