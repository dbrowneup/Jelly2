#!/usr/bin/env python

import argparse

from Setup import Setup
from Support import Support
from Assembly import Assembly
from Placement import Placement

def usage():
    return """
    python Jelly2.py [options] <Scaffolds.fasta> <PacBio_Subreads.bam>
    
    Jelly2 will take the input scaffolds and subreads, and process them
    with four stages: setup, support, assembly, and placement.
    
    The setup stage takes the input scaffolds, identifies gaps, and extracts
    gap-flanking sequences. The left and right sequences flanking gaps are written 
    to new, separate Fasta files, for the alignment of PacBio reads with BLASR. 
    By aligning only to flanking sequences, we save time and increase specificity.
    Performing alignments on the separated flanks allows easy identification of PacBio
    reads that bridge gaps by finding reads that align properly to each flank.
    
    Fasta header output format:
        $ cat Scaffold_gaps.L.fa
        >Scaffold_1.gap.1.L
        ATGC
    
        $ cat Scaffold_gaps.R.fa
        >Scaffold_1.gap.1.R
        ATGC
    
    A gap table describing the scaffolds is written in standard BED format:
        Scaffold_1 gap_start gap_end
    
    The support stage performs the alignment with BLASR, loads the alignments and flanks 
    into memory with pysam and pyfaidx respectively. If there are enough reads supporting 
    a gap, the reads that span the gap are written to disk for assembly. The assembly stage
    goes through each gap and assembles the reads with Minimap, Miniasm, and Racon. The placement
    stage takes gap sequences that were assembled into a single contig and writes a new Fasta file 
    with the gap-filled scaffolds.
    """

def main():
    parser = argparse.ArgumentParser(description='Determine PacBio read support for gaps in scaffolds', usage=usage())
    # Main positional arguments
    parser.add_argument('scaffolds', action='store', help='The input scaffolds in Fasta format')
    parser.add_argument('subreads', action='store', help='The PacBio subreads in BAM format')
    # Arguments for Setup
    setup_args = parser.add_argument_group('Setup')
    setup_args.add_argument('-n', '--min_gap', dest='min_gap', type=int, default=25, \
        help='Minimum number of consecutive Ns to be considered a gap, default=25')
    setup_args.add_argument('-x', '--max_gap', dest='max_gap', type=int, default=1000000, \
        help='Maximum number of consecutive Ns to be considered a gap, default=Inf')
    setup_args.add_argument('-f', '--flank_size', dest='flank_size', type=int, default=1000, \
        help='Number of extracted bases flanking gaps and scaffold ends, default=1000')
    # Arguments for Support
    support_args = parser.add_argument_group('Support')
    support_args.add_argument('-b', '--blasr', dest='blasr', type=str, \
        help='Parameters to pass to BLASR', default='')
    support_args.add_argument('-d', '--min_reads', dest='min_reads', type=int, \
        help='The minimum number of reads required to support a gap', default=5)
    support_args.add_argument('-w', '--wiggle', dest='wiggle', type=int, \
        help='The percent deviation allowed from predicted gap size', default=0.5)
    # Arguments for Assembly
    assembly_args = parser.add_argument_group('Assembly')
    assembly_args.add_argument('-t', '--threads', dest='threads', type=int, \
        help='Number of threads to use for Minimap', default=1)
    assembly_args.add_argument('-m', '--minimap', dest='minimap', \
        help='Parameters to pass to Minimap', default='-Sw5 -L100 -m0')
    assembly_args.add_argument('-a', '--miniasm', dest='miniasm', \
        help='Parameters to pass to Miniasm', default='')
    assembly_args.add_argument('-r', '--racon', dest='racon', \
        help='Parameters to pass to Racon', default='')
#    # Arguments for Placement
#    placement_args = parser.add_argument_group('Placement')
    # Parse the arguments
    args = parser.parse_args()
    # Run Setup
    setup = Setup()
    setup.run(args)
    # Run Support
    support = Support()
    support.mapping(args)
    support.find_support(args)
    # Run Assembly
    assembly = Assembly()
    assembly.assemble_gaps(args)
    # Run Placement
    placement = Placement(args)
    placement.fill_gaps()



if __name__ == '__main__':
    main()