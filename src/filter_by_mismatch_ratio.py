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
    parser.add_argument('-i', metavar='i', type = file, default = sys.stdin,
                        help='input sam file (without header lines)')
    parser.add_argument('-o', metavar='o', type = argparse.FileType('w'),
                        default = sys.stdout,
                        help='output file')
    parser.add_argument('-e', metavar='e', type = argparse.FileType('w'),
                        default = sys.stderr,
                        help='error sam file (without header lines)')
    parser.add_argument('-l', metavar='l', type=int,
                        default = 10000,
                        help='minimum length (default: 10000)')
    parser.add_argument('-m', metavar='m', type=float,
                        default = 0.1,
                        help='maximum mismatch ratio (default: 0.1)')    
    args = parser.parse_args()
    
    main_extract(in_f = args.i,
                 out = args.o,
                 err = args.e,
                 ref = pysam.FastaFile(args.ref),
                 min_len = args.l, 
                 max_mismatch_rate = args.m)

if __name__ == "__main__":
    main()
