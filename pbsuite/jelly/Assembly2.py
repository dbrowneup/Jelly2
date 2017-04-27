#!/usr/bin/env python

import argparse
import subprocess
from string import Template

def usage():
    return """
    USAGE: python Assembly2.py [options] <gapName>

    This script is a wrapper around Minimap/Miniasm/Racon for 
    assembling the reads and generating consensus.
    """

def main():
	# Parser the arguments
    parser = argparse.ArgumentParser(description='Assemble and polish gap-supporting PacBio reads', usage=usage())
    parser.add_argument('gapName', action='store', help='The PacBio reads to be assembled')
    parser.add_argument('-t', '--threads', dest='threads', type=int \
    	help='Number of threads to use for Minimap', default=1)
    parser.add_argument('-m', '--minimap', dest='minimap', \
    	help='Parameters to pass to Minimap', default='-Sw5 -L100 -m0')
    parser.add_argument('-a', '--miniasm', dest='miniasm', \
    	help='Parameters to pass to Miniasm', default='')
    parser.add_argument('-r', '--racon', dest='racon', \
    	help='Parameters to pass to Racon', default='')
    args = parser.parse_args()
    # Map reads against reads
    reads_path = 'Gap_Support/'+args.gapName+'/reads.fq'
    align_path = 'Gap_Support/'+args.gapName+'/reads.paf.gz'
    mappingTemplate = Template("minimap ${params} ${target} ${query}")
    mappingKey = {"params": args.minimap, "target": reads_path, "query": reads_path}
    with open(align_path, 'w') as output:
        p1 = subprocess.Popen(mappingTemplate.substitute(mappingKey).split(' '), stdout=subprocess.PIPE)
        p2 = subprocess.Popen('gzip -1'.split(' '), stdin=p1.stdout, stdout=output)
        p2.communicate()
    # Assemble reads based on mapping
    graph_path = 'Gap_Support/'+args.gapName+'/assembly.gfa'
    fasta_path = 'Gap_Support/'+args.gapName+'/assembly.fa'
    assemblyTemplate = Template("miniasm ${params} -f ${reads} ${aligns}")
    assemblyKey = {"params", args.miniasm, "reads": reads_path, "aligns": align_path}
    with open(graph_path, 'w') as output:
        p1 = subprocess.Popen(assemblyTemplate.substitute(assemblyKey).split(' '), stdout=output)
        p1.communicate()
    with open(fasta_path, 'w') as output:
        p1 = subprocess.Popen(str("""awk '/^S/{print(">"$2"\n"$3)}' """+graph_path).split(' '), stdout=output)
        p1.communicate()
    # Map reads against assembly
    mappingKey = {"params": "", "target": fasta_path, "query": reads_path}
    align_path = 'Gap_Support/'+args.gapName+'/consensus.paf'
    with open(align_path, 'w') as output:
        p1 = subprocess.Popen(mappingTemplate.substitute(mappingKey).split(' '), stdout=output)
        p1.communicate()
    # Generate consensus
    final_path = 'Gap_Support/'_args.gapName+'/consensus.fa'
    consensusTemplate = Template('racon ${params} ${reads} ${aligns} ${graph} ${final}')
    consensusKey = {"params": args.racon, "reads": reads_path, "aligns": align_path, "graph": graph_path, "final": final_path}
    p1  = subprocess.Popen(consensusTemplate.substitute(consensusKey).split(' '))
    p1.communicate()

if __name__ = '__main__':
    main()