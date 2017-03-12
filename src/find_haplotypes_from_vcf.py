import numpy as np
import pandas as pd
import itertools as it
import logging
import sys
import collections
import pysam

import misc

class VCF:
    '''Class for VCF file
    (very crappy implementatoin)    
    '''
    def __init__(self, parsed_line):        
        self.CHROM  = parsed_line[0]
        self.POS    = int(parsed_line[1])
        self.ID     = parsed_line[2]
        self.REF    = parsed_line[3]
        self.ALT    = parsed_line[4]
        self.QUAL   = parsed_line[5]
        self.FILTER = parsed_line[6]
        self.INFO   = parsed_line[7]
        self.GF     = parsed_line[8]
        self.GT     = parsed_line[9]
        self.raw    = parsed_line
    def is_phased(self):
        return('|' in self.GT)
    def get_GT(self):
        return(np.array([int(x) for x in self.GT.split('|')]))
    def __str__(self):
        return('\t'.join(self.raw))
        
def tabix_lookup(tbx, reference_name, start, end):
    '''Fetch VCF files and return as a list of VCF() objects
    '''
    fetched = tbx.fetch(reference_name, start, end,
                        parser=pysam.asTuple())
    return([VCF(x) for x in fetched])

def find_haplotypes_from_vcf(tbx, bim, chrom, reference_start, reference_end):
    '''read vcf file and return haplotypes as nd-array
    '''
    
    # prepare arrays to write haplotype
    index_l, index_r = bim.find_index_interval(chrom, 
                                               reference_start, 
                                               reference_end)
    hap_ary = np.zeros((index_r - index_l, 2), dtype = np.int32) - 9

    # obtain variant info from vcf file
    variants = tabix_lookup(tbx, 
                            'chr{}'.format(chrom), 
                            reference_start,
                            reference_end)
        
    for v in variants:
        index = bim.find_index(chrom, v.POS)
        if(bim.get_bp(chrom)[index] == v.POS):
            if(v.is_phased):
#                logger.debug('{} {} {} {} {}'.format(v.REF, v.ALT, index-index_l, 
#                                                     bim.get_pri(chrom)[index],
#                                                     bim.get_sec(chrom)[index]))
                if(v.REF == bim.get_pri(chrom)[index]):
                    for allele in np.where(v.get_GT() == 0)[0]:
                        hap_ary[index - index_l][allele] = 0
                elif(v.REF == bim.get_sec(chrom)[index]):
                    for allele in np.where(v.get_GT() == 0)[0]:
                        hap_ary[index- index_l][allele] = 1                    

                if(v.ALT == bim.get_pri(chrom)[index]):
                    for allele in np.where(v.get_GT() == 1)[0]:
                        hap_ary[index - index_l][allele] = 0
                elif(v.ALT == bim.get_sec(chrom)[index]):
                    for allele in np.where(v.get_GT() == 1)[0]:
                        hap_ary[index- index_l][allele] = 1

            else:
                logger.warning('Unphased vcf file has not supported yet')
    return(np.transpose(hap_ary))
