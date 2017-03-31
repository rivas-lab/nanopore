import numpy as np
import pandas as pd
import itertools as it
import logging
import sys
import collections
import pysam
import bisect

import misc
from misc import Mismatch
from misc import Nucleotide

class Haplotype:
    def __init__(self, chromosome, index_l, index_r, hap_ary):
        self.chromosome = chromosome
        self.index_l = index_l
        self.index_r = index_r
        self.hap_ary = hap_ary
        self.hash_key = hash((self.chromosome,
                              self.index_l,
                              self.index_r,
                              self.get_bitstr()))
        
    def get_bitstr(self):
        return(''.join([str(x) if x in [0, 1] else '?' 
                        for x in self.hap_ary]))
    
    def is_on_the_same_region(self, other):
        return(self.chromosome == other.chromosome and
               self.index_l == other.index_l and
               self.index_r == other.index_r)
    
    def dist(self, other):
        '''Hamming distance between two haplotypes
        '''
        if(self.is_on_the_same_region(other)):
            return(np.sum(self.hap_ary != other.hap_ary))
        else:
            return(-1)
    
    def __hash__(self):
        return(self.hash_key)
    
    def __eq__(self, other):
        return(self.__hash__() == other.__hash__())

    def __ne__(self, other):        
        return(not(self == other))
    
    def __str__(self):
        return('Hap_{}:{}-{}_{}'.format(self.chromosome, 
                                    self.index_l, 
                                    self.index_r, 
                                    self.get_bitstr()))