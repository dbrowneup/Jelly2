#!/usr/bin/env python

import argparse

from pyfaidx import Fasta

class Placement():
    def __init__(self, scaffolds, gap_table):
        self.gap = self.load_gaps(gap_table)
        self.ref = Fasta(scaffolds)

    def load_gaps(self, gap_table):
    	# Open gap table and read file into list
        gap_L = open(gap_table).read().split('\n')[:-1]
        gap_L = [x.split('\t') for x in gap_L]
        # Parse table into dictionary where scaffold names are keys
        gap_D = {x[0]: [[x[1], x[2]]] if x[0] not in gap_D else gap_D[x[0]].append([x[1], x[2]]) for x in gap_L}
        # Load assembled gap sequences
        for scaf, gaps in gap_D.iteritems():
            for i, g in enumerate(gaps):
                try:
                    seq = Fasta('Gap_Support/'+str(scaf)+'.gap.'+str(i+1)+'/consensus.fa'
                    if sum(1 for _ in seq) == 1:
                        gap_D[scaf][i].append(str(seq))
                    else:
                    	continue
                except IOError:
                    continue
