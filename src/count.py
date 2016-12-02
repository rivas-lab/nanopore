#!/usr/bin/env python2

import pysam

_README_ = '''
-------------------------------------------------------------------------
count match/mismatches from sam file
-------------------------------------------------------------------------
'''

def read_sam(sam_f, sep = None):
    '''
    read a sam file. assumes there is no header component.
    '''
    with open(sam_f, 'r') as f:
        if(sep != None):
            entries = [line.strip().split(sep) for line in f]
        else:
            entries = [line.strip().split() for line in f]
    return(entries)

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
            aln_qual.append(qual[ptr_seq : ptr_seq + cigar_list[i][0]])
            for j in xrange(cigar_list[i][0]):
                aln_seq.append(seq[ptr_seq])
                aln_ref.append(ref[ptr_ref])
                if(seq[ptr_seq].upper() == ref[ptr_ref].upper()):
                    aln_chr.append('=')
                    if(q_thr <= qual[ptr_seq]):
                        counts_high_q['='] = counts_high_q['='] + 1
                    else:
                        counts_low_q['='] = counts_low_q['='] + 1
                else:
                    aln_chr.append('X')
                    if(q_thr <= qual[ptr_seq]):
                        counts_high_q['X'] = counts_high_q['X'] + 1
                        snps.append((ptr_ref, ref[ptr_ref], seq[ptr_seq]))
                    else:
                        counts_low_q['X'] = counts_low_q['X'] + 1                
                ptr_seq = ptr_seq + 1
                ptr_ref = ptr_ref + 1
        elif(cigar_list[i][1] == 'I'):
            aln_qual.append(qual[ptr_seq : ptr_seq + cigar_list[i][0]])
            aln_seq.append(seq[ptr_seq : ptr_seq + cigar_list[i][0]])
            for j in xrange(cigar_list[i][0]):
                if(q_thr <= qual[ptr_seq]):
                    counts_high_q['I'] = counts_high_q['I'] + 1
                else:
                    counts_low_q['I'] = counts_low_q['I'] + 1
                ptr_seq = ptr_seq + 1
            aln_ref.append(gap_char * cigar_list[i][0])    
            aln_chr.append('I' * cigar_list[i][0])    
        elif(cigar_list[i][1] == 'D'):    
            aln_qual.append(' ' * cigar_list[i][0])
            aln_seq.append(gap_char * cigar_list[i][0])
            aln_ref.append(ref[ptr_ref : ptr_ref + cigar_list[i][0]])
            ptr_ref = ptr_ref + cigar_list[i][0]
            aln_chr.append('D' * cigar_list[i][0])    
            counts_high_q['D'] = counts_high_q['D'] + cigar_list[i][0]                                           
        else:
            aln_qual.append(qual[ptr_seq : ptr_seq + cigar_list[i][0]])
            aln_seq.append(seq[ptr_seq : ptr_seq + cigar_list[i][0]])
            for j in xrange(cigar_list[i][0]):
                if(q_thr <= qual[ptr_seq]):
                    counts_high_q[cigar_list[i][1]] = counts_high_q[cigar_list[i][1]] + 1
                else:
                    counts_low_q[cigar_list[i][1]] = counts_low_q[cigar_list[i][1]] + 1
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

def format_counts(list, qscore_t = None):
    strs = []
    if(qscore_t != None):
        strs.append('\t'.join(["{:>6}\t".format(c) 
                    for c in ['=', 'X', 'I', 'D', 'N', 'S', 'H', 'P']]))
        strs.append('\t'.join(['>={:2d}'.format(qscore_t), '<{:2d}'.format(qscore_t)] * 8))
        return('\n'.join(strs + [i.format_count(True) for i in list]))
    else:
        strs.append('\t'.join(["{:>6}".format(c) 
                               for c in ['=', 'X', 'I', 'D', 'N', 'S', 'H', 'P']]))
        return('\n'.join(strs + [i.format_count(False) for i in list]))

def get_counts_main(sam_f, reference, gap_char = '_', qval_thr = 20):
    entries = read_sam(sam_f)
    data = [process_entry(e, reference, gap_char, qval_thr) for e in entries]
    return(data)

def main():
    sam_f = '/home/ytanigaw/data/nanopore/20161008_wgs_caucasian_48hr.10k.bwa.mapq60.20kb.chr11.sam.body'
    ref_f = '/share/PI/mrivas/data/hg19/hg19.fa'
    hg19 = pysam.FastaFile(ref_f)
    
    data = get_counts_main(sam_f, hg19, gap_char = '_', qval_thr = 20)
    
    print format_counts(data, qscore_t = 20)
    
if __name__ == "__main__":
    main()
    