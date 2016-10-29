#!/usr/bin/env python3

_README_ = '''
-------------------------------------------------------------------------
count match/mismatch/indel (python3 implementaion)
Author: Yosuke Tanigawa (ytanigaw@stanford.edu)
-------------------------------------------------------------------------
'''

import argparse

class maf_item():
    # This class can store an alignment in a maf file
    # 
    # For example,
    #   a score=27 EG2=4.7e+04 E=2.6e-05
    #   s humanMito 2170 145 + 16571 AGTAGGCCTAAAAGCAGCCACCAATTAAGAAAGCGTT...
    #   s fuguMito  1648 142 + 16447 AGTAGGCTTAGAAGCAGCCACCA--CAAGAAAGCGTT...
    # This data can be stored as follows
    #  entry = maf_item('score=27 EG2=4.7e+04 E=2.6e-05')
    #  entry.append('humanMito', 2170, 145, '+', 16571, 'AGTAGGCCTAAAAGCAGCCACCAATTAAGAAAGCGTT...')
    #  entry.append('fuguMito',  1648, 142, '+', 16447, 'AGTAGGCTTAGAAGCAGCCACCA--CAAGAAAGCGTT...')
    # http://last.cbrc.jp/doc/last-tutorial.html
    def __init__(self, a):
        self._a = a
        self._name    = []
        self._name_hash = {}
        self._start   = []
        self._alnSize = []
        self._strand  = []
        self._seqSize = []
        self._alignment = []
        self.other_fields = {}
    def __str__(self):
        strs = []
        strs.append(self._a + " ({0} sequences)".format(self.size()))
        return('\n'.join(strs))
    def size(self):
        return(len(self._name))
    def append(self, l):
        name, start, alnSize, strand, seqSize, alignment = l[0], l[1], l[2], l[3], l[4], l[5]
        # add a sequence into alignment
        self._name_hash[name] = self.size()
        self._name.append(name)
        self._start.append(start)
        self._alnSize.append(alnSize)
        self._strand.append(strand)
        self._seqSize.append(seqSize)
        self._alignment.append(alignment.upper())
    def add_field(self, key, val):
        # add some additional information about alignment
        self.other_fields[key] = val
        
    def count_sub(self, ref_id = 0, tar_id = 1):
        # count match/mismatch/insertion/deletion
        if(self.size() <= max(ref_id, tar_id)):
            return([-1, -1, -1, -1])
        else:
            match = 0
            mismatch = 0
            insertion = 0
            deletion = 0
            seq_ref = self.get(ref_id)[5]
            seq_tar = self.get(tar_id)[5]
            for position in range(len(seq_ref)):
                if(seq_ref[position] == '-'):
                    insertion += 1
                elif(seq_tar[position] == '-'):
                    deletion += 1
                elif(seq_ref[position] == seq_tar[position]):
                    match += 1
                else:
                    mismatch += 1
            return([match, mismatch, insertion, deletion])
    def count(self, ref_id = 0, tar_id = 1):
        count_stats = []
        count_stats += self.count_sub(ref_id, tar_id)
        count_stats += self.get(ref_id)[:-1]
        count_stats += self.get(tar_id)[:-1]
        return(count_stats)
    def get_id(self, name):
        return(self._name_hash[name])
    def get_by_name(self, name):
        return(self.get(self.get_id(name)))
    def get(self, i):
        # get a component of the alignment
        return([self._name[i], self._start[i], self._alnSize[i], 
                self._strand[i], self._seqSize[i], self._alignment[i]])
    def dump(self):
        # show the alignment
        strs = []
        strs.append(self._a + " ({0} sequences)".format(self.size()))
        strs.append('\t'.join(['name', 'start', 'alnSize', 'strand', 'seqSize', 'alignment']))
        strs += ['\t'.join([str(x) for x in self.get(i)]) for i in range(self.size())]
        return('\n'.join(strs))

def count_main(maf_file, out_file, debug):
    entry = None
    strs = []
    strs.append('\t'.join(['match', 'mismatch', 'insertion', 'deletion',
                           'ref_name', 'ref_start', 'ref_alnSize', 'ref_strand', 'ref_seqSize',
                           'tar_name', 'tar_start', 'tar_alnSize', 'tar_strand', 'tar_seqSize']))

    for l in maf_file:
        if((not l.startswith('#')) and len(l.strip()) > 0):
            if(l.startswith('a')):
                if(entry != None):
                    strs.append('\t'.join([str(x) for x in entry.count(ref_id=0, tar_id=1)]))
                entry = maf_item(l[2:])
            elif(l.startswith('s')):
                entry.append(l.split()[1:])
            else:
                splitted_line = l.split()
                entry.add_field(splitted_line[0], ' '.join(splitted_line[1:]))
    if(entry != None):
        strs.append('\t'.join([str(x) for x in entry.count(ref_id=0, tar_id=1)]))

    if(out_file != None):
        with open(out_file, 'w') as f:
            f.write('\n'.join(strs))
    else:
        print('\n'.join(strs))

def main():
    # This function serves as a parser
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=_README_)
    parser.add_argument('-i', metavar='i', type=argparse.FileType('r'),
                        #default = None,
                        help='input file (maf format)')
    parser.add_argument('-o', metavar='o', # type=argparse.FileType('w'),
                        default = None, 
                        help='output file name (default = stdout)')
    parser.add_argument('--debug', action = 'store_const',
                        const=True, default = False,
                        help='debug mode (***)')
    parser.add_argument('--version', action='version', version='%(prog)s 2016-10-28')
    
    args = parser.parse_args()

    count_main(args.i, args.o, args.debug)
    
if __name__ == "__main__":
    main()
    
