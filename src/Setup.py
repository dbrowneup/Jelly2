#!/usr/bin/env python

import os
import re
import sys

from pyfaidx import Fasta

class Setup():
    def __init__(self):
        return

    def run(self, args):
        # Open Fastaq and BED output files
        basename = '.'.join(args.scaffolds.split('.')[:-1])
        gapsL = open(basename+'_gaps.L.fa', 'w')
        gapsR = open(basename+'_gaps.R.fa', 'w')
        gapTableOut = open(basename+'_gapInfo.bed', 'w')
        # Load the reference Fasta file
        try:
            reference = Fasta(args.scaffolds)
        except IOError:
            sys.exit("Cannot open reference Fasta, please check file integrity")
        # Implementing flank extraction procedure
        for scaf in reference:
            # Determine gap coordinates in scaffold
            seq = str(scaf)
            gapCoords = []
            query = "[^Nn]([Nn]{%d,%d})[^Nn]"
            param = (args.min_gap, args.max_gap)
            for gap in re.finditer(query % param, seq):
                gapCoords.append((gap.start(), gap.end()))
            # Continue if no gaps are found
            if len(gapCoords) == 0:
                continue
            # Consolidate gaps that are too close (indicates low quality regions)
            i = 0
            while i < len(gapCoords)-1:
                if gapCoords[i+1][0] - gapCoords[i][1] < args.flank_size:
                    gapCoords[i+1][0] = gapCoords[i][0]
                    del gapCoords[i]
                else:
                    i += 1
            # Iterate through gaps and write flanks to fasta, coordinates to bed
            for i, coords in enumerate(gapCoords):
                gapStart, gapEnd = coords
                gapTableOut.write("%s\t%i\t%i\n" % (str(scaf.name), gapStart, gapEnd))
                flankL = seq[gapStart-args.flank_size:gapStart]
                flankR = seq[gapEnd:gapEnd+args.flank_size]
                gapsL.write('>'+str(scaf.name)+'.gap.'+str(i+1)+'.L\n'+flankL+'\n')
                gapsR.write('>'+str(scaf.name)+'.gap.'+str(i+1)+'.R\n'+flankR+'\n')
        # Close shop
        gapsL.close()
        gapsR.close()
        gapTableOut.close()
