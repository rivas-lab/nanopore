#!/usr/bin/env python3

_README_ = '''
-------------------------------------------------------------------------
count match/mismatch/indel (python3 implementaion)
Author: Yosuke Tanigawa (ytanigaw@stanford.edu)
-------------------------------------------------------------------------
'''

import numpy as np
import argparse

class stats():
    def __int__(self, match, mismatch, insertion, deletion):
        self.data = np.array([match, mismatch, insertion, deletion])
    def match(self):
        return(self.data[0])
    def mismatch(self):
        return(self.data[1])
    def insertion(self):
        return(self.data[2])
    def deletion(self):
        return(self.data[3])    
    def __str__(self):
        return(np.array_str(self.data))

class aln():
    def __init__(self, a_filed, alignment):
        self.a_filed = a_filed
        self.src = []
        self.start = []
        self.size = []
        self.strand = []
        self.srcSize = []
        self.text = []
        self.other = []
        for entry in alignment:
            if(entry.startswith('s')):
                elements = entry.split()
                self.src.append(elements[1])
                self.start.append(elements[2])
                self.size.append(elements[3])
                self.strand.append(elements[4])
                self.srcSize.append(elements[5])
                self.text.append(elements[6].upper())
            else:
                self.other.append(entry)
        self.num_seq = len(self.src)
    def count_stats_pair(self, ref, target):
        insertion = 0
        deletion = 0
        match = 0
        mismatch = 0
        for position in range(len(self.text[ref])):
            if(self.text[ref] == '-'):
                insertion += 1
            elif(self.text[target] == '-'):
                deletion += 1
            elif(self.text[ref] == self.text[target]):
                match += 1
            else:
                mismatch += 1
        stat = stats(match, mismatch, insertion, deletion)
        return(stat.data)
    def count_stats(self, ref = 0, targets = None):
        if(targets == None):
            return([self.count_stats_pair(ref, j) for j in range(self.num_seq)])
        else:
            return([self.count_stats_pair(ref, j) for j in targets])
    def __str__(self):
        strs = []
        for i in range(self.num_seq):
            strs.append(self.text[i])
        return('\n'.join(strs))
    
    
def count_mismatches_main(maf_file, out_file, debug):    
    alignment_block_str = []
    for l in maf_file:
        if(not l.startswith('#') and len(l.strip()) > 0):
            if(l.startswith('a')):
                if(len(alignment_block_str) > 0):
                    alignment = aln(alignment_block_str[0], alignment_block_str[1:])
                    print(alignment.count_stats())
                    alignment_block = []
                    alignment_block_str.append(l.strip())
            else:
                alignment_block_str.append(l.strip())            
    alignment = aln(alignment_block_str[0], alignment_block_str[1:])
    print(alignment.count_stats())
        
        
def main():
    # This function serves as a parser
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=_README_)
    parser.add_argument('i', metavar='i', type=argparse.FileType('r'),
                        help='input file (maf format)')
    parser.add_argument('-o', metavar='o', type=argparse.FileType('w'),
                        default = './count.out', 
                        help='output file name (default = ./count.out)')
    parser.add_argument('--debug', action = 'store_const',
                        const=True, default = False,
                        help='debug mode (***)')
    parser.add_argument('--version', action='version', version='%(prog)s 2016-10-27')
    
    args = parser.parse_args()

    count_mismatches_main(args.i, args.o, args.debug)
    
if __name__ == "__main__":
    main()