import numpy as np
import pandas as pd
import itertools as it
import logging
import sys
import collections
import pysam
import bisect

import pgenlib as pg

import misc
from misc import Mismatch
from misc import Nucleotide
from haplotype import Haplotype

class population_reference():
    def __init__(self, params):
        self.params  = params
        self.pgens   = {}
        
    def load_pgen(self, chromosome):
        self.pgens[str(chromosome)] =\
            pg.PgenReader(misc.get_pgen_file_name(chromosome=chromosome,
                                                  param_obj = self.params))   

    def is_loaded(self, chromosome):
        return(str(chromosome) in self.pgens)
        
    def get_pgen(self, chromosome):
        if(not self.is_loaded(chromosome)):
            self.load_pgen(chromosome)
        return(self.pgens[str(chromosome)])
    
    def read_haplotype(self, chromosome, index_l, index_r):
        pgen = self.get_pgen(chromosome)
        sample_ct = pgen.get_raw_sample_ct()
        alleles_list = np.zeros((index_r - index_l, 2 * sample_ct),
                                dtype = np.int32)
        pgen.read_alleles_range(index_l, index_r, alleles_list)
        return([Haplotype(chromosome, index_l, index_r, alleles_list[:, i]) 
                for i in range(alleles_list.shape[1])])