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

class bim_file():
    def __init__(self, params):
        self.params  = params
        self.raw_tables  = {}
        self.id = {}
        self.morgan = {}
        self.bp = {}
        self.pri = {}
        self.sec = {}
        
    def load_raw_table(self, chromosome):
        self.raw_tables[str(chromosome)] = \
            pd.read_csv(misc.get_bim_file_name(chromosome=chromosome,
                                               param_obj = self.params),
                        sep = '\t', 
                        names = ['chr', 'id', 'morgan', 'bp', 'pri', 'sec'])            
        self.id[    str(chromosome)] = np.array(self.raw_tables[str(chromosome)].ix[:, 1])            
        self.morgan[str(chromosome)] = np.array(self.raw_tables[str(chromosome)].ix[:, 2])            
        self.bp[    str(chromosome)] = np.array([int(x) for x in self.raw_tables[str(chromosome)].ix[:, 3]])
        self.pri[   str(chromosome)] = np.array(self.raw_tables[str(chromosome)].ix[:, 4])
        self.sec[   str(chromosome)] = np.array(self.raw_tables[str(chromosome)].ix[:, 5])            
        
    def is_loaded(self, chromosome):
        return(str(chromosome) in self.raw_tables)
        
    def get_raw_table(self, chromosome):
        if(not self.is_loaded(chromosome)):
            self.load_raw_table(chromosome)
        return(self.raw_tables[str(chromosome)])  
    
    def get_id(self, chromosome):
        if(not self.is_loaded(chromosome)):
            self.load_raw_table(chromosome)
        return(self.id[str(chromosome)])  

    def get_morgan(self, chromosome):
        if(not self.is_loaded(chromosome)):
            self.load_raw_table(chromosome)
        return(self.morgan[str(chromosome)])  
    
    def get_bp(self, chromosome):
        if(not self.is_loaded(chromosome)):
            self.load_raw_table(chromosome)
        return(self.bp[str(chromosome)])  

    def get_pri(self, chromosome):
        if(not self.is_loaded(chromosome)):
            self.load_raw_table(chromosome)
        return(self.pri[str(chromosome)])  

    def get_sec(self, chromosome):
        if(not self.is_loaded(chromosome)):
            self.load_raw_table(chromosome)
        return(self.sec[str(chromosome)])  
    
    def find_index(self, chromosome, position):
        bp = self.get_bp(chromosome)
        return(min(len(bp) - 1,
                   max(0,
                       bisect.bisect_left(bp, position))))

    def find_index_list(self, chromosome, positions):
        return(np.array([self.find_index(chromosome, position)
                         for position in positions]))
    
    def find_index_interval(self, chromosome, pos_l, pos_r):
        '''Find semi-open interval of indicies of SNPs that 
           overlaps with a mapped fragment spanning [pos_l, pos_r)
        '''
        index_l = self.find_index(chromosome, pos_l)
        index_r = self.find_index(chromosome, pos_r)
        return((index_l, index_r))
    
    def find_index_exact(self, chromosome, position):
        '''
        Find index on bim file by chromosom and position
        If there is an exact hit, return the index on bim
        else, return None to indicate it is NOT polymorphic region
        '''
        return(index 
               if (self.get_bp(chromosome)[self.find_index(chromosome, 
                                                           position)] 
                   == position)
               else None)

    def find_index_exact_list(self, chromosome, positions):
        return(np.array([self.find_index_exact(chromosome, position)
                         for position in positions]))        