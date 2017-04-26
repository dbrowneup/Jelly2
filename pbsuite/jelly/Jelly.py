#!/usr/bin/env python
"""
document the protocol
"""

import os
import sys
import time
import logging
import argparse
import subprocess
from string import Template
from xml.etree import ElementTree
from glob import glob

STAGES = ["setup", "mapping", "support", "extraction", "assembly", "output"]

from pbsuite.jelly import Stages
from pbsuite.utils.setupLogging import *
from pbsuite.utils.CommandRunner import * 

def usage():
    return """
    USAGE: Jelly.py <stage> <protocol.xml> [-x "<options for stage>"]
    
    Jelly is the driver for each stage of the 
    reference genome upgrade. 
    
    <stage> is one of
        setup
        mapping
        support
        extraction
        assembly
        output
    <protocol.xml> contains the information about the
    data and parameters Jelly will run. See README.txt
    or the documentation for the Protocol's format.
    """

class JellyProtocol():
    """
    Independent method to parse protocols
    """
    def __init__(self, fileName):
        self.protocolName = os.path.abspath(fileName)
        self.parseProtocol()
    
    def parseProtocol(self):
        # Load protocol into memory with ElementTree
        try:
            protocol = ElementTree.parse(self.protocolName)
        except Exception:
            logging.error(("Invalid Protocol"))
            sys.exit(1)
        root = protocol.getroot()
        refNode = root.find("reference")
        # Check for correct protocol elements
        if refNode is None:
            logging.error("Protocol doesn't have <reference> element.")
            sys.exit(1)
        # Check for reference fasta file
        self.scaffoldName = refNode.text
        if not os.path.exists(self.scaffoldName):
            logging.error("Reference %s Does Not Exist!" % (self.scaffoldName))
            sys.exit(1)
        self.referenceNameBase = self.scaffoldName[:self.scaffoldName.rindex(".fasta")]
        self.gapTable = self.referenceNameBase+".gapInfo.bed"
        # Load base directory into memory
        inputNode = root.find("input")
        if inputNode.attrib.has_key("baseDir"):
            self.baseDir = inputNode.attrib["baseDir"]
        else:
            self.baseDir = ""
        # Check integrity of each input
        self.inputs = []
        for part in inputNode:
            fullPath = os.path.join(self.baseDir, part.text)
            if not os.path.exists(fullPath):
                logging.error("Input %s doesn't exit! Exiting" % (fullPath))
                exit(0)
            if not (fullPath.lower().endswith('.fasta') or fullPath.lower().endswith('.fastq')):
                logging.error("Input %s needs to end with .fasta or .fastq! Exiting" % (fullPath))
                exit(0)
            self.inputs.append(fullPath)
        if len(self.inputs) == 0:
            logging.error("Protocol doesn't specifiy any input inputs!")
            sys.exit(1)
        # Assign output directory
        outputNode = root.find("outputDir")
        if outputNode is None:
            logging.warning("Output directory not specified. Using pwd. Hope you're cool with that...")
            self.outDir = os.getcwd()
        else:
            self.outDir = outputNode.text
        # Load command execution mechanism from protocol
        self.runXML = root.find("cluster")
        # Assign any BLASR parameters from protocol
        blasrNode = root.find("blasr")
        if blasrNode is None:
            logging.warning("No blasr parameters!?")
            self.blasrParams = ""
        else:
            self.blasrParams = blasrNode.text

class JellyRunner():
    """
    Take a JellyProtocol and loads in the variables in a way that 
    JellyRunner can use it easily!
    """
    def __init__(self):
        # Parse arguments: $ python Jelly.py <stage> <protocol> [-x "<options>"]
        parser = argparse.ArgumentParser(description='jelly2 genome improvement tool', usage=usage())
        parser.add_argument('stage', help="the stage you're trying to run", action='store')
        parser.add_argument('protocol', help='XML-formatted protocol file', action='store')
        parser.add_argument('-x', help="options to pass into the stage you're running",
                            dest='extras', metavar='<options>', default='')
        parser.add_argument('--debug', help='enable debugging',
                            dest='debug', action='store_true')
        self.options = parser.parse_args()
        self.executeStage = self.options.stage
        self.protocolName = os.path.abspath(self.options.protocol)
        # Modulate logging, print citation, parse protocol
        setupLogging(self.options.debug)
        sys.stderr.write("""
             Please Cite: English, Adam C., Stephen Richards, Yi Han, Min Wang,
                 Vanesa Vee, Jiaxin Qu, Xiang Qin, et al. "Mind the
                 Gap: Upgrading Genomes with Pacific Biosciences RS
                 Long-Read Sequencing Technology." PLoS ONE 7, no. 11
                 (November 21, 2012): e47768.
                 doi:10.1371/journal.pone.0047768.\n\n""")
        self.parseProtocol()
        
    def parseProtocol(self):
        self.protocol = JellyProtocol(self.protocolName)
        self.parseCluster(self.protocol.runXML)
        
    def parseCluster(self, xmlNode):
        # Set up the command running mechanism based on protocol input
        if xmlNode is None:
            self.runCmd = CommandRunner()
        else:
            command = xmlNode.find("command")
            if command is None:
                logging.error(("""You're trying to use a cluster template, but you didn't 
                                specify the template. Please read the documentation."""))
                sys.exit(1)
            nJobs = xmlNode.find("nJobs")
            if nJobs is None or nJobs.text == '0':
                logging.warning(("Not specifying number of jobs may overload clusters."))
                nJobs = 0
            else:
                nJobs = int(nJobs.text)
            cmdTemplate = command.text
            self.runCmd = CommandRunner(cmdTemplate, nJobs)

    def run(self):
        logging.info("Executing Stage: %s" % self.executeStage)
        if self.options.debug:
            Stages.DEBUG = "--debug"
        # Setup before a stage and run
        if self.executeStage == "setup":
            wDir = os.path.dirname(self.protocol.scaffoldName)
            myCommands = [Stages.setup(self.protocol.scaffoldName, self.protocol.gapTable, self.options.extras)]

        elif self.executeStage == "mapping":
            wDir = os.path.join(self.protocol.outDir, "mapping")
            try:
                os.mkdir(wDir)
            except OSError:
                logging.debug("%s already exists" % wDir)
            if not os.path.exists(wDir):
                logging.warning("%s was not created. Check write permissions!" % wDir)
                sys.exit(1)
            myCommands = Stages.mapping(self.protocol.inputs, wDir, self.protocol.scaffoldName, \
                self.protocol.referenceIndex, self.protocol.blasrParams, self.options.extras)

        elif self.executeStage == "support":
            wDir = os.path.join(self.protocol.outDir, "support")
            try:
                os.mkdir(wDir)
            except OSError:
                logging.debug("%s already exists" % wDir)
            if not os.path.exists(wDir):
                logging.warning("%s was not created. Check write permissions!" % wDir)
                sys.exit(1)
            myCommands = Stages.support(self.protocol.outDir, self.protocol.gapTable, wDir, self.options.extras)

        elif self.executeStage == "extraction":
            wDir = os.path.join(self.protocol.outDir, "assembly")
            try:
                os.mkdir(wDir)
            except OSError:
                logging.debug("%s already exists" % wDir)
            if not os.path.exists(wDir):
                logging.warning("%s was not created. Check write permissions!" % wDir)
                sys.exit(1)
            myCommands = [Stages.extraction(wDir, self.protocolName, self.options.extras)]

        elif self.executeStage == "assembly":
            wDir = os.path.join(self.protocol.outDir, "assembly")
            myCommands = Stages.assembly(wDir, self.protocol.gapTable, self.options.extras)

        elif self.executeStage == "output":
            wDir = os.path.join(self.protocol.outDir, "assembly")
            myCommands = [Stages.collection(self.protocol.outDir, self.protocol, self.options.extras)]
            
        logging.debug("CommandRunner Returned: " + str(self.runCmd(myCommands, wDir, self.executeStage)))
        logging.info("Finished %s Stage: %s" % (self.runCmd.runType, self.executeStage))
        
if __name__ == '__main__':
    prog = JellyRunner()
    prog.run()
