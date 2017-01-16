#!/usr/bin/env python2

import fileinput, argparse, sys
import pysam
from utils import *

_README_ = '''
-------------------------------------------------------------------------
Show SNPs
-------------------------------------------------------------------------
'''

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=_README_)
    parser.add_argument('ref', metavar='r',
                        help='reference file (hg19.fa)')    
    parser.add_argument('-i', metavar='i', type = file, default = sys.stdin,
                        help='input sam file (without header lines)')
    parser.add_argument('-o', metavar='o', type = argparse.FileType('w'), default = sys.stdout,
                        help='output file')
    parser.add_argument('-e', metavar='e', type = argparse.FileType('w'), default = sys.stderr,
                        help='error sam file (without header lines)')
    parser.add_argument('-s', metavar='s',
                        default = None,
                        help='vcf file (dbSNP)')    
    parser.add_argument('--validation', metavar='V',
                        default = None,
                        help='vcf file (for validation)')        
    parser.add_argument('-q', metavar='q', type=int, default = -1,
                        help='Base call Q-score threshold')
    parser.add_argument('-g', metavar='g',
                        default = '_',
                        help="gap char (default: '_')")
    parser.add_argument('--showname', action = 'store_const',
                        const=True, default = False,
                        help= 'print name of the sequence') 

    
    args = parser.parse_args()
    
    main_dump_snps(in_f = args.i, out = args.o, err = args.e, 
                   ref = pysam.FastaFile(args.ref),
                   gap_char = args.g,
                   qval_thr = args.q,                     
                   showname = args.showname, 
                   vcf_f = args.s,
                   validation_f = args.validation)

if __name__ == "__main__":
    main()
