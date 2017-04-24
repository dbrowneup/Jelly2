#!/usr/bin/env python

import os
import re
import logging
import argparse

from pyfaidx import Fasta

from pbsuite.utils.setupLogging import *

def USAGE():
    return """
    python Setup.py [<options>] <inputScaffolding.fasta>

    Take the input scaffolds, identify gaps, and extract flanking regions.
    The sequences flanking gaps, as well as the scaffold ends, are written
    to a new Fasta file, for the alignment of PacBio reads. By aligning only
    to sequences flanking gaps, we save time and increase specificity.
    """

class Setup():
    def __init__(self):
        parser = argparse.ArgumentParser(description='process input for jelly2 pipeline', usage=USAGE())
        parser.add_argument('scaffolds', action='store', help='The input scaffolds in Fasta format')
        parser.add_argument('-g', '--gapOutput', dest='gapOutput', \
            help="Create the table for gapInformation", default=None)
        parser.add_argument("--debug",action="store_true", help="Increases verbosity of logging" )
        parser.add_argument("--minGap", dest="minGap", type="int", default=25, \
            help="Minimum number of consecutive Ns to be considered a gap. default=25")
        parser.add_argument("--maxGap", dest="maxGap", type="int", default=1000000, \
            help="Maximum number of consecutive Ns to be considered a gap default=Inf")
        parser.add_argument('--flankSize', dest='flankSize', type='int', default=1000, \
            help='Number of extracted bases flanking gaps and scaffold ends, default=1000')
        self.options = parser.parse_arguments()
        if not os.path.isfile(self.options.scaffolds):
            parser.error("Error! Scaffold File is not a file / does not exist")
        if not self.options.scaffolds.endswith(".fasta") or not self.options.scaffolds.endswith(".fa"):
            parser.error("Reference must end in extension .fasta or .fa! Please rename it.")
        setupLogging(self.options.debug)

    def run(self):
        """
        Fasta header output format:
            >Scaffold_1.end.L
            >Scaffold_1.end.R
            >Scaffold_1.gap.1.L
            >Scaffold_1.gap.1.R
        
        Flanking sequences of n bp are stored for each scaffold and for each gap.
        The gap table is written in standard BED format:
            Scaffold_1 gapStart gapEnd
        """
        # Fasta reference output
        flankName = '.'.join(self.options.scaffolds.split('.')[:-1])+'.flanks.fasta'
        flankOutput = open(flankName, 'w')      
        # Gaps output
        if self.options.gapOutput is not None:
            gapTableOut = open(self.options.gapOutput,'w')
        # Implementing new flank extraction procedure
        reference = Fasta(self.options.scaffolds)
        for scaf in reference:
            # Extract scaffold end sequences
            seq = str(scaf)
            if 'N' not in seq[:self.options.flankSize]:
                flankOutput.write(">"+str(scaf.name)+'.end.L\n'+seq[:self.options.flankSize]+'\n')
            if 'N' not in seq[-self.options.flankSize:]:
                flankOutput.write(">"+str(scaf.name)+'.end.R\n'+seq[-self.options.flankSize:]+'\n')
            # Determine gap coordinates in scaffold
            gapCoords = []
            query = "[^Nn]([Nn]{%d,%d})[^Nn]"
            param = (self.options.minGap, self.options.maxGap)
            for gap in re.finditer(query % param, seq):
                gapCoords.append((gap.start()+1, gap.end()-1))
            # Continue if no gaps are found
            if len(gapCoords) == 0:
                gapTableOut.write("\t".join([str(scaf.name), 'na', 'na'])+'\n')
                logging.debug("Scaffold %s is empty" % str(scaf.name))
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
                if self.options.gapOutput is not None:
                    gapTableOut.write("%s\t%i\t%i" % (scaff.name, gapStart, gapEnd))
                flankL = seq[gapStart-self.options.flankSize:gapStart]
                flankR = seq[gapEnd:gapEnd+self.options.flankSize]
                flankOutput.write('>'+str(scaf.name)+'.gap.'+str(i+1)+'.L\n'+flankL+'\n')
                flankOutput.write('>'+str(scaf.name)+'.gap.'+str(i+1)+'.R\n'+flankR+'\n')
        # Close shop
        flankOutput.close()
        if self.options.gapOutput is not None:
            gapTableOut.close()
        logging.info("Finished!")

if __name__ == '__main__':
    setup = Setup()
    setup.run()
