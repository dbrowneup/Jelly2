#!/usr/bin/env python

import argparse
import subprocess
from string import Template


class Assembly():
    def __init__(self):
        return
    
    def assemble_gaps(self, args):
        gaps = open('Supported_Gaps.txt', 'r').read().split('\n')[:-1]
        gaps = [x.split('\t') for x in gaps]
        for gap in gaps:
            # Map reads against reads
            reads_path = 'Gap_Support/'+gap[0]+'/reads.fq'
            align_path = 'Gap_Support/'+gap[0]+'/reads.paf.gz'
            mappingTemplate = Template("minimap ${params} ${target} ${query}")
            mappingKey = {"params": args.minimap, "target": reads_path, "query": reads_path}
            with open(align_path, 'w') as output:
                p1 = subprocess.Popen(mappingTemplate.substitute(mappingKey).split(' '), stdout=subprocess.PIPE)
                p2 = subprocess.Popen('gzip -1'.split(' '), stdin=p1.stdout, stdout=output)
                p2.communicate()
            # Assemble reads based on mapping
            graph_path = 'Gap_Support/'+gap[0]+'/assembly.gfa'
            fasta_path = 'Gap_Support/'+gap[0]+'/assembly.fa'
            assemblyTemplate = Template("miniasm ${params} -f ${reads} ${aligns}")
            assemblyKey = {"params": args.miniasm, "reads": reads_path, "aligns": align_path}
            with open(graph_path, 'w') as output:
                p3 = subprocess.Popen(assemblyTemplate.substitute(assemblyKey).split(' '), stdout=output)
                p3.communicate()
            with open(fasta_path, 'w') as output:
                p4 = subprocess.Popen(str("""awk '/^S/{print(">"$2"\n"$3)}' """+graph_path).split(' '), stdout=output)
                p4.communicate()
            # Map reads against assembly
            mappingKey = {"params": "", "target": fasta_path, "query": reads_path}
            align_path = 'Gap_Support/'+gap[0]+'/consensus.paf'
            with open(align_path, 'w') as output:
                p5 = subprocess.Popen(mappingTemplate.substitute(mappingKey).split(' '), stdout=output)
                p5.communicate()
            # Generate consensus
            final_path = 'Gap_Support/'+gap[0]+'/consensus.fa'
            consensusTemplate = Template('racon ${params} ${reads} ${aligns} ${graph} ${final}')
            consensusKey = {"params": args.racon, "reads": reads_path, "aligns": align_path, "graph": graph_path, "final": final_path}
            p6 = subprocess.Popen(consensusTemplate.substitute(consensusKey).split(' '))
            p6.communicate()
