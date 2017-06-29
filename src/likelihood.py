from __future__ import print_function

import sys
import os
import logging
from logging.config import dictConfig
import numpy as np
import pandas as pd
import itertools as it
import collections as cl
import bisect

from numba import jit

logging_config = dict(
    version = 1,
    formatters = {
        'f': {'format':
              '%(asctime)s %(name)-12s %(levelname)-8s %(message)s'}
        },
    handlers = {
        'h': {'class': 'logging.StreamHandler',
              'formatter': 'f',
              'level': logging.DEBUG}
        },
    root = {
        'handlers': ['h'],
        'level': logging.DEBUG,
        },
)
dictConfig(logging_config)

def read_prior_cnts(block_df, prior_count_dir):
    
    logger_cnt = logging.getLogger('read_prior_cnts')    
    logger_cnt.info(
        'reading prior counts from {}'.format(prior_count_dir)
    )        
    
    prior_cnt_keys = [None] * len(block_df)
    prior_cnt_vals = [None] * len(block_df)
    for block_id in range(len(block_df)):
        if(block_id % 100 == 0):
            logger_cnt.info(
                'reading block {} of {}'.format(block_id, len(block_df))
            )    
        cnt = np.load('{}/{}.npz'.format(prior_count_dir, block_id))
        prior_cnt_keys[block_id] = cnt['keys']
        prior_cnt_vals[block_id] = cnt['vals']
        
    logger_cnt.info(
        'prior counts is loaded on memory'
    )    
        
    return prior_cnt_keys, prior_cnt_vals

@jit
def find_block_intervals_from_hap(block_df, hapStart, hapEnd):
    """Align two cordinates (one on read, the other on block_df)"""
    hap_len = hapEnd - hapStart
    hap_to_assign_len = hap_len 
    start_block = bisect.bisect_right(block_df.bimStart.as_matrix(), hapStart) - 1
    current_block = start_block
    
    res = dict([])
    
    while(hap_to_assign_len > 0):
        if(start_block == current_block):
            current_block_SNP_start = hapStart - block_df.bimStart[current_block] 
        else:
            current_block_SNP_start = 0
        if(hap_to_assign_len < block_df.n_SNPs[current_block]):
            current_block_SNP_end = hap_to_assign_len 
            res[current_block] = ((hap_len - hap_to_assign_len, 
                                   hap_len - hap_to_assign_len + current_block_SNP_end - current_block_SNP_start),
                                  (current_block_SNP_start, current_block_SNP_end))
            hap_to_assign_len = 0
        else:
            current_block_SNP_end = block_df.n_SNPs[current_block] 
            res[current_block] = ((hap_len - hap_to_assign_len, 
                                   hap_len - hap_to_assign_len + current_block_SNP_end - current_block_SNP_start),
                                  (current_block_SNP_start, current_block_SNP_end))
            hap_to_assign_len = hap_to_assign_len - (current_block_SNP_end - current_block_SNP_start)  
            current_block = current_block + 1    
    
    return res

def find_ld_blocks(line, block_df):
    logger_find_ld_blocks = logging.getLogger('find_ld_blocks')
    logger_find_ld_blocks.debug('passing {} {}'.format(line[6], line[7]))
    return find_block_intervals_from_hap(block_df, int(line[6]), int(line[7]))

def find_read_specific_error_rate(line):
    return float(line[5])

def find_observation(line):
    obs1d = np.array([1 - int(x) for x in list(line[-1])], dtype=np.int8)
    return obs1d.reshape((1, len(obs1d)))

@jit
def compute_likelihood_from_a_read(line, block_df, prior_cnt_keys):
    """Given one line of the input hap file, compute and update the log-likelihood
    """
    
    logger_read_ll = logging.getLogger('compute_likelihood_from_a_read')
    
    obs = find_observation(line)
    err_rate = find_read_specific_error_rate(line)
    
    likelihood_of_read = dict([])
    
    logger_read_ll.debug('parsed ld blocks are {}'.format(str(find_ld_blocks(line, block_df))))
    
    for block_id, v in find_ld_blocks(line, block_df).items():
        logger_read_ll.debug('The LD we are working on is: {} {}'.format(block_id, v))
        
        obs_s, obs_e = v[0]
        block_s, block_e = v[1]    
                
        n_mismatch = np.sum(
            prior_cnt_keys[block_id][:, block_s:block_e] != obs[:, obs_s:obs_e], axis = 1
        )
        n_match = (obs_e - obs_s) - n_mismatch
        
        likelihood_of_read[block_id] = n_mismatch * np.log(err_rate) + n_match * np.log(1 - err_rate)

    return likelihood_of_read

def update_likelihood_dict(likelihood_of_read, log_likelihood):   
    for block_id, ll in likelihood_of_read.items():
        if(block_id in log_likelihood):
            log_likelihood[block_id] = log_likelihood[block_id] + ll
        else:
            log_likelihood[block_id] = ll
    return log_likelihood

def update_ll(hap_f, block_df, prior_cnt_keys, log_likelihood = dict([])):

    logger_ll = logging.getLogger('update_ll')
        
    with open(hap_f, 'r') as f:
        for line_num, raw_line in enumerate(f):
            line  = raw_line.rstrip().split('\t')
            
#            if(line_num % 100 == 0):            
            if True:
                logger_ll.info(
                    'processing read {} {}'.format(line_num, line[3])
                )                                

            likelihood_of_read = compute_likelihood_from_a_read(line, block_df, prior_cnt_keys)            
            log_likelihood = update_likelihood_dict(likelihood_of_read, log_likelihood)
    return log_likelihood

def main():
    read_dir='/oak/stanford/groups/mrivas/public_data/nanopore-wgs-consortium/rel3/hg19/chr20'
    hap_basename='rel3.chr20.12500.10k-chr20impv1-keep-maf0005-snv-biallelic-geno01-hwe1e-10.hap'
    hap_f='{}/{}'.format(read_dir, hap_basename)

    data_dir='/oak/stanford/groups/mrivas/users/ytanigaw/nanopore-data'

    block_tsv_f='{}/{}'.format(
        data_dir,
        'chr20impv1-keep-maf0005-snv-biallelic-geno01-hwe1e-10-block-stronglow050-stronghigh083-infofrac10.tsv'
    )

    prior_count_dir='{}/{}'.format(data_dir, 'prior_count')
    log_likelihood_dir='{}/{}/{}'.format(data_dir, 'log_likelihood', hap_basename[:-4])
    if not os.path.exists(log_likelihood_dir):
        os.makedirs(log_likelihood_dir)    

    block_df = pd.read_csv(block_tsv_f, sep='\t')

    prior_cnt_keys, prior_cnt_vals = read_prior_cnts(block_df, prior_count_dir)

    log_likelihood = update_ll(hap_f, block_df, prior_cnt_keys)

    for block_id, ll in log_likelihood.items():
        np.savez(
            '{}/ll{}.npz'.format(log_likelihood_dir, block_id), 
            ll = ll
        )

if __name__ == "__main__":
    main()
    
