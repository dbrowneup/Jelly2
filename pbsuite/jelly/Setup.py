#!/usr/bin/env python

import re
import os
import sys
import logging
import argparse
import subprocess
#from optparse import OptionParser

from pyfaidx import Fasta

from pbsuite.utils.setupLogging import *
from pbsuite.utils.FileHandlers import FastaFile, Gap, wrap, qwrap
from pbsuite.utils.CommandRunner import exe

def USAGE():
    return """
    python Setup.py [<options>] <inputScaffolding.fasta>

    Take the input scaffolds, identify gaps, and extract flanking regions.
    The sequences flanking gaps, as well as the scaffold ends, are written
    to a new Fasta file, for the alignment of PacBio reads. By aligning only
    to sequences flanking gaps, we save time and increase specificity.
    """

refParser = re.compile("(.*)\|(ref\d{7})/?(\d+)?$")

class Setup():
    def __init__(self):
        parser = argparse.ArgumentParser(description='process input for jelly2 pipeline', usage=USAGE())
        parser.add_argument('scaffolds', action='store', help='The input scaffolds in Fasta format')
        parser.add_argument('-g', '--gapOutput', action='store_true', \
            help="Create the table for gapInformation")
        parser.add_argument('-i', '--index', dest='index', action='store_true', default=False, \
            help="Create the .sa index for faster blasr alignments")
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
#        qn = self.options.scaffolds[:self.options.scaffolds.rindex('.fasta')]
#        qualInputName = qn+".qual"
#        if not os.path.isfile(qualInputName):
#            self.qualInput = None
#        else:
#            self.qualInput = qualInputName
        self.qualInput = None

    def run(self):
        """
        Fasta header output format:
            >Scaffold_1.flank.L
            >Scaffold_1.flank.R
            >Scaffold_1.gap1.L
            >Scaffold_1.gap1.R
        
        Flanking sequences of n bp are stored for each scaffold and for each gap.
        The gap table is written in standard BED format:
            Scaffold_1 gapStart gapEnd
        """
        #Fasta Ref Output
        scaffTempName = '.'.join(self.options.scaffolds.split('.')[:-1])+'.flanks.fasta'
        scaffOutput = open(scaffTempName, 'w')     
#        #Qual Ref Output
#        if self.qualInput is not None:
#            qualTempName= self.qualInput+".tempQual"
#            qualOutput = open(qualTempName, 'w')   
        #Gaps Output
        if self.options.gapOutput:
            gapTableOut = open('.'.join(self.options.scaffolds.split('.')[:-1])+'.gaps.bed','w')
#        else:
#            gapTableOut = False
#        logging.info("Creating reference sequence index names and identifying gaps")
#        refTemplate = "ref%07d"
#        refId = 1
#        #Read References
#        reference = FastaFile(self.options.scaffolds)
#        if self.qualInput is not None:
#            qualReference = QualFile(self.qualInput)
        # Implementing new flank extraction procedure
        reference = Fasta(self.options.scaffolds)
        for scaf in reference:
            # Extract scaffold end sequences
            seq = str(scaf)
            if 'N' not in seq[:self.options.flankSize]:
                scaffOutput.write(">"+str(scaf.name)+'.flank.L\n'+seq[:self.options.flankSize]+'\n')
            if 'N' not in seq[-self.options.flankSize:]:
                scaffOutput.write(">"+str(scaf.name)+'.flank.R\n'+seq[-self.options.flankSize:]+'\n')
            # Determine gap coordinates in scaffold
            gapCoords = []
            query = "[^Nn]([Nn]{%d,%d})[^Nn]"
            param = (self.options.minGap, self.options.maxGap)
            for gap in re.finditer(query % param, seq):
                gapCoords.append((gap.start() + 1, gap.end() - 1))
            # Continue if no gaps are found
            if len(gapCoords) == 0:
                gapTableOut.write("\t".join([scaf.name, 'na', 'na'])+'\n')
                logging.debug("Scaffold %s is empty" % scaffName)
                continue
            # Consolidate gaps that are too close -- indicating low quality regions.
            i = 0
            while i < len(gapCoords)-1:
                if gapCoords[i+1][0] - gapCoords[i][1] < self.options.flankSize:
                    gapCoords[i+1][0] = gapCoords[i][0]
                    del gapCoords[i]
                else:
                    i += 1
            # Iterate through gaps and write to table
            for i, coords in enumerate(gapCoords):
                gapStart, gapEnd = coords
                if self.options.gapOutput:
                    gapTableOut.write("%s\t%i\t%i" % (scaff.name, gapStart, gapEnd))
                flankL = seq[gapStart-1000:gapStart]
                flankR = seq[gapEnd:gapEnd+1000]
                scaffOutput.write('>'+str(scaff.name)+'.gap'+str(i+1)+'.L\n'+flankL+'\n')
                scaffOutput.write('>'+str(scaff.name)+'.gap'+str(i+1)+'.R\n'+flankR+'\n')

#        for key in reference:
#            scaffIndex = refTemplate % refId
#            scaffName = key.replace(' ','_')
#            refId += 1
#            scaffName = scaffName + "|" + scaffIndex
#            scaffOutput.write(">"+scaffName+"\n"+wrap(reference[key])+"\n")
#            if self.qualInput is not None:
#                qualOutput.write(">"+scaffName+"\n"+qwrap(qualReference[key])+"\n")
#            gapCoords = []
#            for gap in re.finditer("[^Nn]([Nn]{%d,%s})[^Nn]" % (self.options.minGap, self.options.maxGap), reference[key]):
#                gapCoords.append([gap.start() + 1, gap.end() - 1])
#            if len(gapCoords) == 0: #no Gaps
#                gapTableOut.write("\t".join([scaffName, 'na', 'na', scaffIndex+"_0_0", '3'])+'\n')
#                logging.debug("Scaffold %s is empty" % scaffName)
#                continue
#            #Consolidate gaps that are too close -- indicating LQ regions.
#            i = 0
#            while i < len(gapCoords)-1:
#                if gapCoords[i+1][0] - gapCoords[i][1] < 25:
#                    gapCoords[i+1][0] = gapCoords[i][0]
#                    del gapCoords[i]
#                else:
#                    i += 1
#            prevEnd = 0#Contig Start Tracking
#            idx = 0
#            #Make the first gap
#            prevEnd = gapCoords[0][1]
#            gapCoords[0][1]-gapCoords[0][0]
#            flag = Gap.BEGIN
#            if len(gapCoords) == 1:
#                flag += Gap.END
#            if gapTableOut:
#                gapTableOut.write("%s\t%i\t%i\t%s_%i_%i\t%d\n" \
#                    % (scaffName, gapCoords[0][0], gapCoords[0][1], scaffIndex, idx, idx+1, flag))
#            #Now Go Through the rest of the gaps
#            for i in range(1, len(gapCoords)):
#                idx += 1
#                prevEnd = gapCoords[i][1]
#                gapCoords[i][1]-gapCoords[i][0]
#                if gapTableOut:
#                    if i == len(gapCoords)-1:
#                        flag = Gap.END
#                    else:
#                        flag = 0
#                    gapTableOut.write("%s\t%i\t%i\t%s_%i_%i\t%d\n" \
#                        % (scaffName, gapCoords[i][0], gapCoords[i][1], scaffIndex, idx, idx+1, flag))
        #Close shop
        scaffOutput.close()
        os.rename(self.options.scaffolds, self.options.scaffolds+".original")
        os.rename(scaffTempName, self.options.scaffolds)
#        if self.qualInput is not None:
#            qualOutput.close()
#            os.rename(self.qualInput, self.qualInput+".original")
#            os.rename(qualTempName, self.qualInput)
        if gapTableOut:
            gapTableOut.close()
#        if self.options.index:
#            logging.info("Creating .sa indexes for references")
#            r, o, e = exe("sawriter %s.sa %s" % (self.options.scaffolds, self.options.scaffolds))
#            if r != 0:
#                logging.error("sawriter returned %d" % r)
#                logging.error("Ensure it's in your path")
#                exit(1)
#            logging.debug(str(o) + ' ' + str(e))       
        logging.info("Finished!")

if __name__ == '__main__':
    setup = Setup()
    setup.run()
