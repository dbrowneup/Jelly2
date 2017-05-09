#!/usr/bin/env python

from pyfaidx import Fasta

class Placement():
    def __init__(self, args):
        self.basename = '.'.join(args.scaffolds.split('.')[:-1])
        self.gap = self.load_gaps(self.basename+'_gapInfo.bed')
        self.ref = Fasta(args.scaffolds)

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
                    seq = [str(f) for f in Fasta('Gap_Support/'+str(scaf)+'.gap.'+str(i+1)+'/consensus.fa')]
                    if len(seq) == 1:
                        gap_D[scaf][i].append(str(seq[0]))
                    else:
                        continue
                except IOError:
                    continue
        return gap_D

    def fill_gaps(self):
        with open(self.basename+'_GapFilled.fa', 'w') as out:
            for scaf in self.ref:
                # Try to load gaps, but if scaffold has no gaps, write to output and continue
                try:
                    gaps = self.gap[str(scaf.name)]
                except KeyError:
                    out.write('>'+str(scaf.name)+'\n'+str(scaf)+'\n')
                    continue
                # If gaps were loaded, parse and write to output
                last_index = 0
                scaffold = list()
                for i, g in enumerate(gaps):
                    # Load scaffold sequence preceding gap
                    scaffold.append(scaf[last_index:g[0]])
                    # Load fill sequence if gap was assembled, else fill with N
                    try:
                        scaffold.append(g[2])
                    except IndexError:
                        scaffold.append('N' * (g[1] - g[0]))
                    # Move marker to end of current gap
                    last_index = g[1]
                else:
                    # Load final stretch of scaffold sequence past last gap
                    scaffold.append(scaf[last_index:])
                # Write scaffold to output
                out.write('>'+str(scaf.name)+'\n'+''.join(scaffold)+'\n')
