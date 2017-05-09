#!/usr/bin/env python

from Setup import Setup
from Support import Support
from Assembly import Assembly
from Placement import Placement

def main():
    parser = argparse.ArgumentParser(description='Determine PacBio read support for gaps in scaffolds', usage=usage())
    parser.add_argument('scaffolds', action='store', help='The input scaffolds in Fasta format')
    parser.add_argument('gap_info', action='store', help='The input gap info in BED format')
    parser.add_argument('-m', '--min_reads', dest='min_reads', type=int, default=5, \
        help='The minimum number of reads required to support a gap')
    parser.add_argument('-w', '--wiggle', dest='wiggle', type=int, default=0.5, \
        help='The percentage of deviation allowed from predicted gap size')
    args = parser.parse_args()
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



if __name__ == '__main__':
    main()