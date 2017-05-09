#!/usr/bin/env python

import os
import re
import sys
import argparse

from pyfaidx import Fasta

def usage():
    return """
    python Setup.py [options] <inputScaffolding.fasta>

    Take the input scaffolds, identify gaps, and extract flanking regions.
    The sequences flanking gaps, as well as the scaffold ends, are written
    to new Fasta files, for the alignment of PacBio reads. By aligning only
    to flanking sequences, we save time and increase specificity.
    """

class Setup():
    def __init__(self):
        parser = argparse.ArgumentParser(description='process input for jelly2 pipeline', usage=usage())
        parser.add_argument('scaffolds', action='store', help='The input scaffolds in Fasta format')
        parser.add_argument('-g', '--gapOutput', dest='gapOutput', \
            help="Create the table for gapInformation", default=None)
        parser.add_argument("--minGap", dest="minGap", type=int, default=25, \
            help="Minimum number of consecutive Ns to be considered a gap, default=25")
        parser.add_argument("--maxGap", dest="maxGap", type=int, default=1000000, \
            help="Maximum number of consecutive Ns to be considered a gap, default=Inf")
        parser.add_argument('--flankSize', dest='flankSize', type=int, default=1000, \
            help='Number of extracted bases flanking gaps and scaffold ends, default=1000')
        self.options = parser.parse_args()
        if not os.path.isfile(self.options.scaffolds):
            sys.exit("Error! Scaffold File is not a file / does not exist")
        if not self.options.scaffolds.endswith(".fasta") and not self.options.scaffolds.endswith(".fa"):
            sys.exit("Reference must end in extension .fasta or .fa! Please rename it.")

    def run(self):
        """
        Output stream is 4 separate Fasta files containing left and right flanking
        sequences for scaffolds and gaps within them. PacBio reads can then be aligned
        to each Fasta file separately, selecting the best hit for each read. The alignments
        can then be loaded into memory and reads that bridge gaps can be identified.

        Fasta header output format:
            Scaffold_ends.L.fa
            >Scaffold_1.end.L
            ATGC

            Scaffold_ends.R.fa
            >Scaffold_1.end.R
            ATGC

            Scaffold_gaps.L.fa
            >Scaffold_1.gap.1.L
            ATGC

            Scaffold_gaps.R.fa
            >Scaffold_1.gap.1.R
            ATGC
        
        Flanking sequences of X bp are stored for each scaffold and for each gap.
        The gap table is written in standard BED format:
            Scaffold_1 gapStart gapEnd
        """
        # Open Fasta output files
        basename = '.'.join(self.options.scaffolds.split('.')[:-1])
#        endsL = open(basename+'_ends.L.fa', 'w')
#        endsR = open(basename+'_ends.R.fa', 'w')
        gapsL = open(basename+'_gaps.L.fa', 'w')
        gapsR = open(basename+'_gaps.R.fa', 'w')
        # Open gap BED table output
        try:
            gapTableOut = open(self.options.gapOutput,'w')
        except Exception:
            gapTableOut = open(basename+'_gapInfo.bed', 'w')
        # Load the reference Fasta file
        try:
            reference = Fasta(self.options.scaffolds)
        except Exception:
            sys.exit("Cannot open reference Fasta, please check file integrity")
        # Implementing flank extraction procedure
        for scaf in reference:
            # Extract scaffold end sequences
            seq = str(scaf)
#            endsL.write(">"+str(scaf.name)+'.end.L\n'+seq[:self.options.flankSize]+'\n')
#            endsR.write(">"+str(scaf.name)+'.end.R\n'+seq[-self.options.flankSize:]+'\n')
            # Determine gap coordinates in scaffold
            gapCoords = []
            query = "[^Nn]([Nn]{%d,%d})[^Nn]"
            param = (self.options.minGap, self.options.maxGap)
            for gap in re.finditer(query % param, seq):
                gapCoords.append((gap.start()+1, gap.end()-1))
            # Continue if no gaps are found
            if len(gapCoords) == 0:
                continue
            # Consolidate gaps that are too close -- indicating low quality regions.
            i = 0
            while i < len(gapCoords)-1:
                if gapCoords[i+1][0] - gapCoords[i][1] < self.options.flankSize:
                    gapCoords[i+1][0] = gapCoords[i][0]
                    del gapCoords[i]
                else:
                    i += 1
            # Iterate through gaps and write flanks to fasta, coordinates to bed
            for i, coords in enumerate(gapCoords):
                gapStart, gapEnd = coords
                gapTableOut.write("%s\t%i\t%i\n" % (str(scaff.name), gapStart, gapEnd))
                flankL = seq[gapStart-self.options.flankSize:gapStart]
                flankR = seq[gapEnd:gapEnd+self.options.flankSize]
                gapsL.write('>'+str(scaf.name)+'.gap.'+str(i+1)+'.L\n'+flankL+'\n')
                gapsR.write('>'+str(scaf.name)+'.gap.'+str(i+1)+'.R\n'+flankR+'\n')
        # Close shop
#        endsL.close()
#        endsR.close()
        gapsL.close()
        gapsR.close()
        gapTableOut.close()

if __name__ == '__main__':
    setup = Setup()
    setup.run()
