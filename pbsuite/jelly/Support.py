#!/usr/bin/env python

import os
import argparse
import subprocess
import itertools as it

import pysam
from pyfaidx import Fasta

def usage():
    return """
    USAGE: python Support2.py [options] <Scaffolds.fa> <gapInfo.bed>

    This program will load the reference Fasta containing the scaffold sequences,
    the BAM alignment files of PacBio reads to flanks, and the BED table containing
    the gap coordinates. The alignments are then scanned for each gap to identify
    PacBio reads that bridge the gap. If insufficient reads are found to support
    gap filling, the gap will remain unfilled. The output of this program is a directory 
    with sub-directories for each supported gap in the assembly. Contained in each 
    sub-directory are the reads that bridge that gap in FastQ format. These will 
    subsequently be assembled and polished with Minimap/Miniasm/Racon.
    """

class Support():
    def __init__(self):
        # Parse arguments for program
        parser = argparse.ArgumentParser(description='Determine PacBio read support for gaps in scaffolds', usage=usage())
        parser.add_argument('scaffolds', action='store', help='The input scaffolds in Fasta format')
        parser.add_argument('gap_info', action='store', help='The input gap info in BED format')
        parser.add_argument('-m', '--min_reads', dest='min_reads', type=int, default=5, \
            help='The minimum number of reads required to support a gap')
        parser.add_argument('-w', '--wiggle', dest='wiggle', type=int, default=0.5, \
            help='The percentage of deviation allowed from predicted gap size')
        args = parser.parse_args()

    def mapping(self, reads, scaffoldName, blasr_params):
        # Run the BLASR mapping jobs
        basename = '.'.join(scaffoldName.split('.')[:-1])
        gapsL = {"map": {"reads": reads, "flanks": basename+'_gaps.L.fa', "param": blasr_params}, "out": basename+"_gaps.L.bam"}
        gapsR = {"map": {"reads": reads, "flanks": basename+'_gaps.R.fa', "param": blasr_params}, "out": basename+"_gaps.R.bam"}
        mappingTemplate = Template("blasr ${reads} ${flanks} --bam --hitPolicy allbest ${param}")
        mappingJobs = [gapsL, gapsR]
        for job in mappingJobs:
            with open(job["out"]) as output:
                p1 = subprocess.Popen(mappingTemplate.substitute(job["map"]).split(' '), stdout=subprocess.PIPE)
                p2 = subprocess.Popen('samtools view -F4 -b'.split(' '), stdin=p1.stdout, stdout=output)
                p2.communicate()
        # Index the BAM alignment files
        gapsL = {"aligns": basename+"_gaps.L.bam"}
        gapsR = {"aligns": basename+"_gaps.R.bam"}
        indexingTemplate = Template("samtools index ${aligns}")
        indexingJobs = [gapsL, gapsR]
        for job in indexingJobs:
            subprocess.call(indexingTemplate.substitute(job).split(' '))
    
    def find_support(self):
        # Load BAM alignments, gap BED table, reference Fasta file
        basename = '.'.join(self.scaffolds.split('.')[:-1])
        gapsL = pysam.AlignmentFile(basename+'_gaps.L.bam', 'rb')
        gapsR = pysam.AlignmentFile(basename+'_gaps.R.bam', 'rb')
        gap_list = open(args.gap_info, 'r').read().split('\n')[:-1]
        gap_list = [x.split('\t') for x in gap_list]
        gap_dict = {x[0]: [(x[1], x[2])] if x[0] not in gap_dict else gap_dict[x[0]].append((x[1], x[2])) for x in gap_list}
        ref = Fasta(args.scaffolds)
        # Iterate through scaffolds and determine support
        os.mkdir('Gap_Support')
        for scaf in ref:
        	# Continue if there are no gaps in scaffold
        	try:
                gaps = gap_dict[str(scaf.name)]
            except KeyError:
                continue
            # Iterate through gaps and check support
            for i, gap in enumerate(gaps):
            	gap_size = gap[1] - gap[0]
                readsL = [L for L in gapsL.fetch(str(scaf.name)+'.gap.'+str(i+1)+'.L')])
                readsR = [R for R in gapsR.fetch(str(scaf.name)+'.gap.'+str(i+1)+'.R')])
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
                path = 'Gap_Support/'+str(scaf.name)+'.gap.'+str(i+1)
                os.mkdir(path)
                with open(path+'/reads.fq', 'a') as output:
                    for read in fastq:
                    	output.write(read)
