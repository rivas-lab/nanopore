#!/usr/bin/env python2

import fileinput, argparse, sys
import pysam
from utils import *

_README_ = '''
-------------------------------------------------------------------------
Extract reads in BAM file based on mismatch ratio
-------------------------------------------------------------------------
'''

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=_README_)
    parser.add_argument('ref', metavar='r',
                        help='reference file (hg19.fa)')
    parser.add_argument('-i', metavar='i',
                        default = None,
                        help='input sam file (without header lines)')
    parser.add_argument('-o', metavar='o',
                        default = None,
                        help='output file')
    parser.add_argument('-e', metavar='e',
                        default = None,
                        help='error sam file (without header lines)')
    parser.add_argument('-q', metavar='q', type=int,
                        default = -1,
                        help='Base call Q-score threshold')
    parser.add_argument('-g', metavar='g',
                        default = '_',
                        help="gap char (default: '_')")

    
    args = parser.parse_args()
    
    main_extract(sam_f = args.i, 
                 ref = pysam.FastaFile(args.ref),
                 gap_char = args.g,
                 qval_thr = args.q,                     
                 min_len = 10000, 
                 max_mismatch_rate = 0.1,
                 outfile = args.o,
                 errfile = args.e)

if __name__ == "__main__":
    main()
