#!/usr/bin/env python2

import fileinput, argparse, sys, subprocess
import pysam

_README_ = '''
-------------------------------------------------------------------------
Utilty file
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
    (aln_seq, aln_ref, aln_chr, aln_qual, ptr_ref, 
     counts_high_q, counts_low_q, snps) = \
    retrieve_alignment(seq = e[9],
                       ref = reference.fetch(reference = e[2],
                                             start = int(e[3]) - 1, 
                                             end = int(e[3]) - 1 + ref_len),
                       qual = e[10],
                       cigar_list = cigar_list, 
                       gap_char = gap_char, 
                       qval_thr = qval_thr)
    return(sam_entry(e[0], e[2], int(e[3]), int(e[3]) + ptr_ref, len(e[9]), 
                     counts_high_q, counts_low_q, qval_thr, 
                     aln_seq, aln_ref, aln_chr, aln_qual, snps, e))

class tabix_lookup:
    def __init__(self, vcf_f, rname, pos1, pos2):    
        '''
        look up vcf file with tabix and store the results
        '''
        try:
            query_str = '{rname}:{pos1}-{pos2}'.format(rname = rname,
                                                       pos1  = pos1, 
                                                       pos2  = pos2)
            tabix_res = subprocess.check_output(["tabix", vcf_f, query_str]).split('\n')
        except:
            print "Unexpected error:", sys.exc_info()[0]
            raise
        else:            
            if(len(tabix_res) > 0):
                # store the results table as a nested list                
                self.data = [x.split() for x in tabix_res if len(x) > 0]
            else:
                self.data = []
    def has_hit(self):
        return(len(self.data) > 0)
    def get_exact_match(self, snp):
        exact_match = None
        if(len(self.data) > 0):
            for entry in self.data:
                if(len(entry) > 0 and 
                   entry[0]         == snp.rname and
                   int(entry[1])    == snp.pos and
                   entry[3].upper() == snp.ref.upper() and
                   snp.seq.upper() in [alt.upper() for alt in entry[4].split(',')]):
                    exact_match = entry
                    break
        return(exact_match)

class snp:
    def __init__(self, rname, pos, ref, seq, base_call_q = -1):
        self.rname = rname
        self.pos = pos
        self.ref = ref
        self.seq = seq
        self.base_call_q = base_call_q
        self.var = None
        self.dbsnp = None
        self.dbsnp_tried = False
        self.validation = None
        self.is_valid = None
        
    def dbsnp_lookup(self, dbsnp_f, validation_f = None):
        self.dbsnp = tabix_lookup(dbsnp_f, self.rname, self.pos, self.pos)
        if(self.dbsnp.has_hit()):
            dbsnp_match = self.dbsnp.get_exact_match(self)
            if(dbsnp_match is not None):
                self.var = dbsnp_match[2]
                if(validation_f is not None):
                    self.validation = tabix_lookup(validation_f, self.rname, self.pos, self.pos)
                    validation_match = self.validation.get_exact_match(self)
                    if(validation_match is not None and 
                       self.seq in [alt.upper() for alt in validation_match[4].split(',')]):
                        self.is_valid = True
                    else:
                        self.is_valid = False
            else:
                self.var = '!'
        self.dbsnp_tried = True        
    def has_hit_on_dbsnp(self, dbsnp_f, vcf = None):
        if(not self.dbsnp_tried):
            self.dbsnp_lookup(dbsnp_f, vcf)        
        return(self.dbsnp.has_hit())
    def has_var_id(self, dbsnp_f, vcf = None):
        if(not self.dbsnp_tried):
            self.dbsnp_lookup(dbsnp_f, vcf)
        return(self.var != None and len(self.var) > 1)
    def is_validated(self, dbsnp_f, vcf = None):
        if(not self.dbsnp_tried):
            self.dbsnp_lookup(dbsnp_f, vcf)
        return(self.is_valid)        
    def __str__(self, full = False):
        if(self.dbsnp.has_hit()):            
            return(','.join([str(x) for x in 
                              ['{rname}:{pos}'.format(rname = self.rname, pos = self.pos), 
                               self.ref, self.seq, self.var, self.is_valid]]))
        else:
            return(','.join([str(x) for x in 
                              ['{rname}:{pos}'.format(rname = self.rname, pos = self.pos), 
                               self.ref, self.seq, '*', self.is_valid]]))
                                 
class sam_entry:
    def __init__(self, qname, rname, aln_start, aln_end, seq_len,
                 counts_high_q, counts_low_q, qscore_t,
                 aln_seq, aln_ref, aln_chr, aln_qual, snps, raw):
        self.qname = qname
        self.rname = rname        
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
        self.snps = [snp(rname, s[0] + aln_start, s[1], s[2]) for s in snps]
        self.raw = raw
    def dbsnp_lookup(self, vcf_f, validation_f):
        for s in self.snps:
            s.dbsnp_lookup(vcf_f, validation_f)
    def snp_nums(self, vcf_f, validation_f = None):
        total_num = len(self.snps)
        has_hit_on_dbsnp = sum([int(s.has_hit_on_dbsnp(vcf_f, validation_f)) 
                                for s in self.snps])
        has_var_id       = sum([int(s.has_var_id(vcf_f, validation_f)) 
                                for s in self.snps])
        has_validated    = sum([int(s.is_validated(vcf_f, validation_f) == True) 
                                for s in self.snps])
        return([total_num, has_hit_on_dbsnp, has_var_id, has_validated])
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
    def format_snps(self, vcf_f = None, validation_f = None):
        strs = []
        if(vcf_f is not None):
            self.dbsnp_lookup(vcf_f, validation_f)
        strs.append(';'.join([str(s) for s in self.snps]))
        return('\n'.join(strs))
    
    def counts(self, c):
        return(self.counts_high_q[c] + self.counts_low_q[c])  
    def mismatch_rate(self):
        return(1.0 * self.counts('X') / (self.counts('=') + self.counts('X')))
    def format_count(self, qscore = False, showname = False):
        if(qscore):
            return('\t'.join(([self.name] if showname else []) +
                             sum([[i[0], i[1]] for i in 
                                  zip(["{:6d}".format(self.counts_high_q[c]) 
                                       for c in ['=', 'X', 'I', 'D', 'N', 'S', 'H', 'P']], 
                                      ["{:6d}".format(self.counts_low_q[c]) 
                                       for c in ['=', 'X', 'I', 'D', 'N', 'S', 'H', 'P']])
                                 ], [])))            
        else:
            return('\t'.join(([self.name] if showname else []) + 
                             ["{:6d}".format(self.counts(c)) 
                              for c in ['=', 'X', 'I', 'D', 'N', 'S', 'H', 'P']]))
               
def format_counts_head(qscore_t = -1, showname = False):
    strs = []
    if(qscore_t >= 0):
        strs.append('\t'.join((["name"] if showname else []) + 
                              ["{:>6}\t".format(c) 
                    for c in ['=', 'X', 'I', 'D', 'N', 'S', 'H', 'P']]))
        strs.append('\t'.join(([""] if showname else []) + ['>={:2d}'.format(qscore_t), '<{:2d}'.format(qscore_t)] * 8))
    else:
        strs.append('\t'.join((["name"] if showname else []) + 
                              ["{:>6}".format(c) 
                               for c in ['=', 'X', 'I', 'D', 'N', 'S', 'H', 'P']]))        
    return('\n'.join(strs))            

##########################################################
def main_counts(sam_f, ref, gap_char = '_', qval_thr = -1, showname = False, 
                outfile = None, errfile = None):
    if outfile is None:
        out = sys.stdout
    else:
        out = open(outfile, 'w')

    if errfile is None:
        err = sys.stderr
    else:
        err = open(errfile, 'w')
       
    out.write(format_counts_head(qscore_t = qval_thr, showname = showname) + '\n')
    with open(sam_f, 'r') as f:
        for line in f:
            entry = line.strip().split()
            if(((int(entry[1]) >> 11) % 2) != 0):
                err.write("split alignment\n")
            else:
              # if it is not a supplementary alignment
              try:
                  data = process_entry(entry, ref, gap_char, qval_thr)
              except IndexError as e:
                  sys.stderr.write("Index Error\n")
                  err.write(line)
              except:
                  print "Unexpected error:", sys.exc_info()[0]
                  raise
              else:            
                out.write(data.format_count(qval_thr >= 0, showname) + '\n')
            


def main_extract(in_f, out, err, ref, 
                 min_len, max_mismatch_rate):
    for line in in_f:
        entry = line.strip().split()
        if(((int(entry[1]) >> 11) % 2) == 0):
            # if it is not a supplementary alignment
            try:
                data = process_entry(entry, ref)
            except IndexError as e:
                sys.stderr.write("IndexError\n")
                err.write('\t'.join(["IndexError", line]))
            except ValueError as e:
                sys.stderr.write("ValueError\n")
                err.write('\t'.join(["ValueError", line]))                
            except:
                sys.stderr.write("UnexpectedError\n")
                sys.stderr.write("continuing...\n")
                err.write('\t'.join(["UnexpectedError", line]))                
                #print "Unexpected error:", sys.exc_info()[0]
                #raise
            else:
                if(data.seq_len >= min_len and
                   data.mismatch_rate() <= max_mismatch_rate):
                    out.write(line)

def main_dump_snps(in_f, out, err, ref, 
                   gap_char = '_', qval_thr = -1, showname = False, 
                   vcf_f = None, validation_f = None):     
    snpstr_head = 'snps([<pos>,<ref>,<seq>,<varid>,<validated>;]+)'
    if(showname):
        out.write('\t'.join(['name', '#SNPs', '#SNPs_with_hits_to_dbSNP', '#SNPs_with_var_id', '#SNPs_with_var_id(validated)', snpstr_head]) + '\n')
    else:
        out.write('\t'.join(['#SNPs', '#SNPs_with_hits_to_dbSNP', '#SNPs_with_var_id', '#SNPs_with_var_id(validated)', snpstr_head]) + '\n')
    for line in in_f:
        entry = line.strip().split()
        if(((int(entry[1]) >> 11) % 2) == 0):
            # if it is not a supplementary alignment
            try:
                data = process_entry(entry, ref, gap_char, qval_thr)
            except IndexError as e:
                sys.stderr.write("IndexError\n")
                err.write('\t'.join(["IndexError", line]))
            except ValueError as e:
                sys.stderr.write("ValueError\n")
                err.write('\t'.join(["ValueError", line]))                
            except:
                sys.stderr.write("UnexpectedError\n")
                sys.stderr.write("continuing...\n")
                err.write('\t'.join(["UnexpectedError", line]))                
                #print "Unexpected error:", sys.exc_info()[0]
                #raise
            else:
                snpstr = str(data.format_snps(vcf_f, validation_f))
                snp_nums = data.snp_nums(vcf_f, validation_f)                  
                if(showname):
                    out.write('\t'.join([data.qname] + [str(x) for x in snp_nums] + [snpstr]) + '\n')
                else:
                    out.write('\t'.join([str(x) for x in snp_nums] + [snpstr]) + '\n')
