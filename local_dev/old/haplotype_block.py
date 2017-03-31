import numpy as np
import itertools as it
import logging
import sys
import collections
import pysam
import enum

import misc
from misc import Mismatch
from misc import Nucleotide

class haplotype_block:
    def __init__(self, chromosome, idx_l, idx_r, bim, pgen):
        self.chromosome = chromosome
        self.idx_l = idx_l
        self.idx_r = idx_r        
        self.bim = bim
        self.bim_idx_l, self.bim_idx_r = \
        bim.find_index_interval(chromosome, idx_l, idx_r)
        
        # read_pgen_file
        alleles = pgen.read_alleles(chromosome, self.bim_idx_l, self.bim_idx_r)
        
        self.alleles = alleles # we don't need it. keeping for debug purpose
        
        self.pop_size = alleles.shape[1]
        
        cnt = collections.Counter([tuple(alleles[:, i]) 
                                   for i in range(self.pop_size)])
        
        # number of unique haplotypes
        self.n_uniq_haps = len(cnt.keys())
        
        # find unique haplotypes and their occurance in population
        cnt_sorted = cnt.most_common(self.n_uniq_haps)
        self.haps = np.array([cnt_sorted[i][0] for i
                              in range(self.n_uniq_haps)], dtype = np.int8)
        self.cnts = np.array([cnt_sorted[i][1] for i 
                              in range(self.n_uniq_haps)], dtype = np.int)
        
        # log p(data | haplotype) = \sum_{read \in data} p(read | haplotype)
        self.log_likelihood = np.zeros(self.n_uniq_haps, dtype = np.float)
        self.n_pos_total = 0

    def log_prior(self):
        return(np.log(self.cnts) - np.log(self.pop_size))
    
    def log_posterior(self):
        log_potential = np.log(self.cnts) + self.log_likelihood
        log_part = np.log(np.sum(np.exp(log_potential)))
        return(log_potential - log_part)

    def update_model_read_specific(self, read, bim):
        # take intersection of read and the haplotype block
        bim_min = max(self.bim_idx_l, read.get_bim_interval(bim)[0])
        bim_max = min(self.bim_idx_r, read.get_bim_interval(bim)[1])
        
        # number of positions we consider in this update
        n_pos = max(0, bim_max - bim_min)
        
        # mismatch between reference and read
        mm_ref_vs_read = \
        read.get_mismatches_polymorphic_pos(bim, 
                                            bim.get_bp(self.chromosome)[bim_min], 
                                            bim.get_bp(self.chromosome)[bim_max])
        # polymorphic positions
        poly_pos = bim.get_bp(self.chromosome)[bim_min:bim_max]
        
        # For each haplotype, find set of 
        mm_ref_vs_haps = [poly_pos[self.haps[i][bim_min - self.bim_idx_l :
                                                bim_max - self.bim_idx_l] == 1]
                          for i in range(self.n_uniq_haps)]
        
        # number of mismatches between read and haplotype
        #   need to fix later
        x = np.array([len(set(mm_ref_vs_read) ^ set(mm_ref_vs_haps[i]))
                      for i in range(self.n_uniq_haps)])
                
        e = read.non_polymorphic_error_rate(bim)
        
        self.log_likelihood = self.log_likelihood + np.log(e) * x + np.log(1 - e) * (n_pos - x)
        self.n_pos_total = self.n_pos_total + n_pos
    
    def map_value(self):
        return(np.exp(max(self.log_posterior())))
    
    def map_index(self):
        return(np.where(self.log_posterior() >= max(self.log_posterior()))[0])
    
    def plot(self, filename = None):
        if(filename is not None):
            logger.info('Writing image to {}'.format(filename))
        misc.make_hist(x = np.exp(self.log_posterior()),
                  title = ' '.join(['Haplotype probability distribution',
                                    'after {:.2} observations'.format(1.0 * self.n_pos_total / 
                                                                      (self.bim_idx_r - self.bim_idx_l)),
                                    '\n({} unique haplotypes in'.format(self.n_uniq_haps),
                                    'chr{}:{}-{} ({}SNPs))'.format(self.chromosome,
                                                                   self.idx_l,
                                                                   self.idx_r, 
                                                                   self.bim_idx_r - self.bim_idx_l
                                                                  )]),
                  xlabel = ' '.join(['Haplotype probability',
                                     '(population size: {} haplotypes)'.format(self.pop_size)]),
                  ylabel = 'Frequency of haplotype probabilities', 
                  filename = None)    
    