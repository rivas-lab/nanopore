#!/usr/bin/env python2

import fileinput, argparse
import pysam

_README_ = '''
-------------------------------------------------------------------------
count match/mismatches from sam file
-------------------------------------------------------------------------
'''

def parse_cigar(cigar):
    '''
    parse CIGAR string and return them as list of 2-tupples and
    count the length of reference sequence corresponding to the query
    '''
    parsed = []
    start = 0
    i = 0
    ref_len = 0
    while(i < len(cigar)):
        while('0' <= cigar[i] <= '9'):
            i = i + 1
        segment_length = int(cigar[start:i])
        segment_type = cigar[i]
        if(segment_type != 'I'):
            ref_len = ref_len + segment_length
        parsed.append((segment_length, segment_type))
        i = i + 1
        start = i
    return(parsed, ref_len)

def retrieve_alignment(seq, ref, qual, cigar_list, 
                       gap_char = '_', qval_thr = 20):
    '''
    Based on parsed CIGAR string and read(seq) and reference(ref),
    reconstruct an alignment and count match/mismatches
    '''
    aln_seq = []
    aln_ref = []
    aln_chr = []
    aln_qual = []
    snps = []
    ptr_seq = 0
    ptr_ref = 0
    counts_high_q = {'M':0, 'I':0, 'D':0, 'N':0, 'S':0, 'H':0, 'P':0, '=':0, 'X':0}
    counts_low_q = {'M':0, 'I':0, 'D':0, 'N':0, 'S':0, 'H':0, 'P':0, '=':0, 'X':0}
    q_thr = str(unichr(33 + qval_thr))
    for i in xrange(len(cigar_list)):
        if(cigar_list[i][1] == 'M'):
            if(qual != "*"):
                aln_qual.append(qual[ptr_seq : ptr_seq + cigar_list[i][0]])
            for j in xrange(cigar_list[i][0]):
                aln_seq.append(seq[ptr_seq])
                aln_ref.append(ref[ptr_ref])
                if(seq[ptr_seq].upper() == ref[ptr_ref].upper()):
                    aln_chr.append('=')
                    if(qual != "*" and q_thr > qual[ptr_seq]):
                        counts_low_q['='] = counts_low_q['='] + 1
                    else:
                        counts_high_q['='] = counts_high_q['='] + 1
                else:
                    aln_chr.append('X')
                    if(qual != "*" and q_thr > qual[ptr_seq]):
                        counts_low_q['X'] = counts_low_q['X'] + 1                
                    else:
                        counts_high_q['X'] = counts_high_q['X'] + 1
                        snps.append((ptr_ref, ref[ptr_ref], seq[ptr_seq]))
                ptr_seq = ptr_seq + 1
                ptr_ref = ptr_ref + 1
        elif(cigar_list[i][1] == 'I'):
            if(qual != "*"):
                aln_qual.append(qual[ptr_seq : ptr_seq + cigar_list[i][0]])
            aln_seq.append(seq[ptr_seq : ptr_seq + cigar_list[i][0]])
            for j in xrange(cigar_list[i][0]):
                if(qual != "*" and q_thr > qual[ptr_seq]):
                    counts_low_q['I'] = counts_low_q['I'] + 1
                else:
                    counts_high_q['I'] = counts_high_q['I'] + 1                    
                ptr_seq = ptr_seq + 1
            aln_ref.append(gap_char * cigar_list[i][0])    
            aln_chr.append('I' * cigar_list[i][0])    
        elif(cigar_list[i][1] == 'D'):
            if(qual != "*"):
                aln_qual.append(' ' * cigar_list[i][0])
            aln_seq.append(gap_char * cigar_list[i][0])
            aln_ref.append(ref[ptr_ref : ptr_ref + cigar_list[i][0]])
            ptr_ref = ptr_ref + cigar_list[i][0]
            aln_chr.append('D' * cigar_list[i][0])    
            counts_high_q['D'] = counts_high_q['D'] + cigar_list[i][0]                                           
        else:
            if(qual != "*"):
                aln_qual.append(qual[ptr_seq : ptr_seq + cigar_list[i][0]])
            aln_seq.append(seq[ptr_seq : ptr_seq + cigar_list[i][0]])
            for j in xrange(cigar_list[i][0]):
                if(qual != "*" and q_thr > qual[ptr_seq]):
                    counts_low_q[cigar_list[i][1]] = counts_low_q[cigar_list[i][1]] + 1
                else:                    
                    counts_high_q[cigar_list[i][1]] = counts_high_q[cigar_list[i][1]] + 1
                ptr_seq = ptr_seq + 1            
            aln_ref.append(ref[ptr_ref : ptr_ref + cigar_list[i][0]])
            ptr_ref = ptr_ref + cigar_list[i][0]
            aln_chr.append(cigar_list[i][1] * cigar_list[i][0])
    return(''.join(aln_seq), ''.join(aln_ref), ''.join(aln_chr), ''.join(aln_qual),
           ptr_ref, counts_high_q, counts_low_q, snps)    

def process_entry(e, reference, gap_char = '_', qval_thr = 20):
    '''
    process one line in a sam file and returns an object of sam_entry class
    '''
    cigar_list, ref_len = parse_cigar(e[5])
    (aln_seq, aln_ref, aln_chr, aln_qual, ptr_ref, counts_high_q, counts_low_q, snps) = \
    retrieve_alignment(seq = e[9],
                       ref = reference.fetch(reference = e[2],
                                             start = int(e[3]) - 1, 
                                             end = int(e[3]) - 1 + ref_len),
                       qual = e[10],
                       cigar_list = cigar_list, 
                       gap_char = gap_char)
    return(sam_entry(e[2], int(e[3]), int(e[3]) + ptr_ref, len(e[9]), 
                     counts_high_q, counts_low_q, qval_thr, 
                     aln_seq, aln_ref, aln_chr, aln_qual, snps, e))

class sam_entry:
    def __init__(self, name, aln_start, aln_end, seq_len,
                 counts_high_q, counts_low_q, qscore_t,
                 aln_seq, aln_ref, aln_chr, aln_qual, snps, raw):
        self.name = name
        self.aln_start = aln_start
        self.aln_end = aln_end
        self.seq_len = seq_len
        self.counts_high_q = counts_high_q
        self.counts_low_q = counts_low_q
        self.qscore_t = qscore_t
        self.aln_seq = aln_seq
        self.aln_ref = aln_ref
        self.aln_chr = aln_chr
        self.aln_qual = aln_qual
        self.snps = snps
        self.raw = raw
    def format_aln(self, width = None, qual = False):
        strs = []
        if(width == None):
            width = len(self.aln_chr)
        for batch in xrange(int((len(self.aln_chr) + width - 1) / width)):            
            strs.append(self.aln_seq[batch * width : (batch + 1) * width])
            strs.append(self.aln_chr[batch * width : (batch + 1) * width])
            strs.append(self.aln_ref[batch * width : (batch + 1) * width])
            if(qual):
                strs.append(self.aln_qual[batch * width : (batch + 1) * width])            
            strs.append('')
        return('\n'.join(strs))
    def format_snps(self):
        strs = []
        strs.append('\t'.join(['pos', 'ref', 'read']))
        strs.append('\n'.join('\t'.join(l) for l in 
                              [[str(i[0] + self.aln_start), i[1], i[2]] for i in self.snps]))
        return('\n'.join(strs))
    
    def counts(self, c):
        return(self.counts_high_q[c] + self.counts_low_q[c])        
    def format_count(self, qscore = False):
        if(qscore):
            return('\t'.join(sum([[i[0], i[1]] for i in 
                                  zip(["{:6d}".format(self.counts_high_q[c]) 
                                       for c in ['=', 'X', 'I', 'D', 'N', 'S', 'H', 'P']], 
                                      ["{:6d}".format(self.counts_low_q[c]) 
                                       for c in ['=', 'X', 'I', 'D', 'N', 'S', 'H', 'P']])
                                 ], [])))            
        else:
            return('\t'.join(["{:6d}".format(self.counts(c)) 
                              for c in ['=', 'X', 'I', 'D', 'N', 'S', 'H', 'P']]))
               
def format_counts_head(qscore_t = -1):
    strs = []
    if(qscore_t >= 0):
        strs.append('\t'.join(["{:>6}\t".format(c) 
                    for c in ['=', 'X', 'I', 'D', 'N', 'S', 'H', 'P']]))
        strs.append('\t'.join(['>={:2d}'.format(qscore_t), '<{:2d}'.format(qscore_t)] * 8))
    else:
        strs.append('\t'.join(["{:>6}".format(c) 
                               for c in ['=', 'X', 'I', 'D', 'N', 'S', 'H', 'P']]))        
    return('\n'.join(strs))            

def show_counts_main(sam_f, ref, gap_char = '_', qval_thr = -1):
    with open(sam_f, 'r') as f:
        for line in f:
            entry = line.strip().split()
            data = process_entry(entry, ref, gap_char, qval_thr)
            print data.format_count(qval_thr >= 0)
                
def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=_README_)
    parser.add_argument('ref', metavar='r',
                        help='reference file (hg19.fa)')
    parser.add_argument('-i', metavar='i',
                        default = None,
                        help='input sam file (without header lines)')
    parser.add_argument('-q', metavar='q', type=int,
                        default = -1,
                        help='Base call Q-score threshold')
    parser.add_argument('-g', metavar='g',
                        default = '_',
                        help="gap char (default: '_')")

    
    args = parser.parse_args()

    print format_counts_head()
    show_counts_main(sam_f = args.i, 
                     ref = pysam.FastaFile(args.ref),
                     gap_char = '_',
                     qval_thr = args.q)

if __name__ == "__main__":
    main()
    
