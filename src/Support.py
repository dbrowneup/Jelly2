#!/usr/bin/env python

import os
import subprocess
import itertools as it
from string import Template

import pysam
from pyfaidx import Fasta

class Support():
    def __init__(self):
        return
    
    def mapping(self, args):
        basename = '.'.join(args.scaffolds.split('.')[:-1])
        # Run the BLASR mapping jobs
        mapL = {"reads": args.subreads, "flanks": basename+'_gaps.L.fa', "threads": args.threads, "param": args.blasr, "out": "aligned_gaps.L.bam"}
        mapR = {"reads": args.subreads, "flanks": basename+'_gaps.R.fa', "threads": args.threads, "param": args.blasr, "out": "aligned_gaps.R.bam"}
        mappingTemplate = Template("blasr ${reads} ${flanks} --nproc ${threads} --bam --out ${out} --hitPolicy allbest ${param}")
        for job in [mapL, mapR]:
            subprocess.call(mappingTemplate.substitute(job).split(' '))
    
    def sorting(self, args):
        # Sort the BAM alignment files
        pysam.sort("-@", str(args.threads), "-o", "sorted_gaps.L.bam", "aligned_gaps.L.bam")
        pysam.sort("-@", str(args.threads), "-o", "sorted_gaps.R.bam", "aligned_gaps.R.bam")
    
    def indexing(self, args):
        # Index the BAM alignment files
        indexL = {"aligns": "sorted_gaps.L.bam"}
        indexR = {"aligns": "sorted_gaps.R.bam"}
        indexingTemplate = Template("samtools index ${aligns}")
        for job in [indexL, indexR]:
            subprocess.call(indexingTemplate.substitute(job).split(' '))
    
    def find_support(self, args):
        # Load BAM alignments, gap BED table, reference Fasta file
        basename = '.'.join(args.scaffolds.split('.')[:-1])
        gapsL = pysam.AlignmentFile('sorted_gaps.L.bam', 'rb')
        gapsR = pysam.AlignmentFile('sorted_gaps.R.bam', 'rb')
        gap_list = open(basename+'_gapInfo.bed', 'r').read().split('\n')[:-1]
        gap_list = [x.split('\t') for x in gap_list]
        print "Number of gaps loaded in gap_list:", len(gap_list)
        gap_dict = dict()
        gap_dict = {x[0]: [] for x in gap_list if x[0] not in gap_dict}
        for x in gap_list:
            gap_dict[x[0]].append((int(x[1]), int(x[2])))
        ref = Fasta(args.scaffolds)
        # Iterate through scaffolds and determine support
        os.mkdir('Gap_Support')
        supported_gaps = list()
        for scaf in ref:
            # Continue if there are no gaps in scaffold
            try:
                gaps = gap_dict[str(scaf.name)]
                print str(len(gaps)), "gaps loaded for scaffold", str(scaf.name)
            except KeyError:
                print "No gaps loaded for scaffold", str(scaf.name)
                continue
            # Iterate through gaps and check support
            for i, gap in enumerate(gaps):
                gap_size = gap[1] - gap[0]
                readsL = [L for L in gapsL.fetch(str(scaf.name)+'.gap.'+str(i+1)+'.L')]
                readsR = [R for R in gapsR.fetch(str(scaf.name)+'.gap.'+str(i+1)+'.R')]
                # Determine the number of supporting reads, store alignments in list of tuples
                support = [(L, R) for L, R in it.product(readsL, readsR) if L.query_name == R.query_name]
                print "Scaffold:", str(scaf.name), "Gap:", str(i+1), "Size:", str(gap_size), "Support:", len(support)
                if len(support) < args.min_reads:
                    continue
                # Iterate through supporting reads and measure wiggle
                fasta = list()
                for L, R in support:
                    read_span = R.query_alignment_start - L.query_alignment_end
                    print "Read Span:", str(read_span)
                    if (gap_size - gap_size * args.wiggle) < read_span < (gap_size + gap_size * args.wiggle):
                        sequence = L.query_sequence[L.query_alignment_end:R.query_alignment_start]
#                        quality = L.query_qualities[L.query_alignment_end:R.query_alignment_start]
                        fasta.append('>'+str(L.query_name)+'\n'+str(sequence)+'\n')
                # Continue if too many reads fail wiggle-check
                if len(fastq) < args.min_reads:
                    continue
                # Create sub-directory and write FastQ output
                gap_name = str(scaf.name)+'.gap.'+str(i+1)
                supported_gaps.append((gap_name, str(len(fasta))))
                path = 'Gap_Support/'+gap_name
                os.mkdir(path)
                with open(path+'/reads.fa', 'a') as output:
                    for read in fasta:
                        output.write(read)
        with open('Supported_Gaps.txt', 'w') as output:
            for gap_name, support in supported_gaps:
                output.write(gap_name+'\t'+support+'\n')
