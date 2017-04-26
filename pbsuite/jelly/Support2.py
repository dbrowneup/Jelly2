#!/usr/bin/env python

import os
import argparse
import subprocess
import itertools as it

import pysam
import networkx as nx
from pyfaidx import Fasta

def usage():
    return """
    USAGE: python Support2.py [<options>] <Scaffolds.fa> <gapInfo.bed>

    This program will load the reference Fasta containing the scaffold sequences,
    the BAM alignment files of PacBio reads to flanks, and the BED table containing
    the gap coordinates. The alignments are then scanned for each gap to identify
    PacBio reads that bridge the gap. If insufficient reads are found to support
    gap filling, the gap will remain unfilled. Links between scaffold ends with
    good support will be treated as gaps and filled, joining the scaffolds. The
    output of this program is a directory with sub-directories for each supported 
    gap in the assembly. Contained in each sub-directory are the reads that bridge
    that gap in FastQ format. These will subsequently be assembled with minimap/miniasm
    and polished with minimap/racon.
    """


"""
This Class is intended to store gap information, including the gap coordinates,
the flanking sequences, and the bridging PacBio reads
"""
class GapGraph():
	def __init__(self):
		self.G = nx.DiGraph()


"""
This Class will hold the scaffold information, including the inter-gap sequences,
the set of GapGraph objects, and edges between scaffold ends
"""
class ScaffoldGraph():
	def __init__(self):
		self.G = nx.DiGraph()


"""
This Class will hold the methods for parsing information from the BAM alignments,
the BED table of gap coordinates, and the reference Fasta file
"""
class ReadParser():
    def __init__(self, basename):
    	self.endsL = pysam.AlignmentFile(basename+'_ends.L.bam', 'rb')
    	self.endsR = pysam.AlignmentFile(basename+'_ends.R.bam', 'rb')
    	self.gapsL = pysam.AlignmentFile(basename+'_gaps.L.bam', 'rb')
    	self.gapsR = pysam.AlignmentFile(basename+'_gaps.R.bam', 'rb')


def main():
	# Parse arguments for program
    parser = argparse.ArgumentParser(description='Determine PacBio read support for gaps in scaffolds', usage=usage())
    parser.add_argument('scaffolds', action='store', help='The input scaffolds in Fasta format')
    parser.add_argument('gap_info', action='store', help='The input gap info in BED format')
    parser.add_argument('-m', '--min_reads', dest='min_reads', type=int, default=5, \
    	help='The minimum number of reads required to support a gap')
    parser.add_argument('-w', '--wiggle', dest='wiggle', type=int, default=0.5, \
    	help='The percentage of deviation allowed from predicted gap size')
    args = parser.parse_args()
    # Load BAM alignments, gap BED table, reference Fasta file
    basename = '.'.join(self.options.scaffolds.split('.')[:-1])
    reads = ReadParser(basename)
    gaps_L = open(args.gap_info, 'r').read().split('\n')
    gaps_L = [x.split('\t') for x in gaps_L]
    gaps_D = {x[0]: [(x[1], x[2])] if x[0] not in gaps_D else gaps_D[x[0]].append((x[1], x[2])) for x in gaps_L}
    ref = Fasta(args.scaffolds)
    # Iterate through scaffolds and determine support
    os.mkdir('Gap_Support')
    for scaf in ref:
    	# Continue if there are no gaps in scaffold
    	try:
            gaps = gaps_D[str(scaf.name)]
        except KeyError:
            continue
        # Iterate through gaps and check support
        for i, gap in enumerate(gaps):
        	gap_size = gap[1] - gap[0]
            readsL = [L for L in reads.endsL.fetch(str(scaf.name)+'.gap.'+str(i)+'.L')])
            readsR = [R for R in reads.endsR.fetch(str(scaf.name)+'.gap.'+str(i)+'.R')])
            # Determine the number of supporting reads, store alignments in tuple
            support = [(L, R) if L.query_name == R.query_name for L, R in it.product(readsL, readsR)]
            if len(support) < args.min_reads:
            	continue
            # Iterate through supporting reads and measure wiggle
            fastq = list()
            for L, R in support:
            	read_span = R.query_alignment_start - L.query_alignment_end
            	if (gap_size - gap_size * args.wiggle) < read_span < (gap_size + gap_size * args.wiggle):
            		sequence = L.query_sequence[L.query_alignment_end:R.query_alignment_start]
            		quality = L.query_qualities[L.query_alignment_end:R.query_alignment_start]
            		fastq.append('@'+str(L.query_name)+'\n'+str(sequence)+'\n+\n'+str(quality)+'\n')
            # Continue if too many reads fail wiggle-check
            if len(fastq) < args.min_reads:
                continue
            # Create sub-directory and write FastQ output
            path = 'Gap_Support/'+str(scaf.name)+'.gap.'+str(i)
            os.mkdir(path)
            with open(path+'/reads.fq', 'a') as output:
                for line in fastq:
                	output.write(line)



if __name__ == '__main__':
	main()