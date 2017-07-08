from __future__ import print_function

_README_ = '''
-------------------------------------------------------------------------
Compute the log-likelihood of reads given population reference and reads

Author: Yosuke Tanigawa (ytanigaw@stanford.edu)
Date: 2017/07/07
-------------------------------------------------------------------------
'''

import sys
import os
import errno
import argparse
import logging
from logging.config import dictConfig
import numpy as np
import pandas as pd
import itertools as it
import collections as cl
import bisect

from numba import jit

from read_npz import *

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
        'level': logging.INFO,
        },
)
dictConfig(logging_config)


@jit
def find_block_intervals_from_hap(bins_df, hapStart, hapEnd):
    """Align two cordinates (one on read, the other lon bins_df)"""
    
    hap_len = hapEnd - hapStart
    hap_to_assign_len = hap_len 
    start_block = bisect.bisect_right(bins_df.bimIdStart.as_matrix(), hapStart) - 1
    current_block = start_block

    res = dict([])
    
    while(hap_to_assign_len > 0):
        
        if(start_block == current_block):
            current_block_SNP_start = hapStart - bins_df.bimIdStart[current_block] 
        else:
            current_block_SNP_start = 0
            
        if(hap_to_assign_len < bins_df.nSNPs[current_block] - current_block_SNP_start):
            current_block_SNP_end = hap_to_assign_len + current_block_SNP_start
            
            res[current_block] = ((hap_len - hap_to_assign_len, 
                                   hap_len - hap_to_assign_len + current_block_SNP_end - current_block_SNP_start),
                                  (current_block_SNP_start, current_block_SNP_end))
            hap_to_assign_len = 0
        else:
            current_block_SNP_end = bins_df.nSNPs[current_block] 
            
            res[current_block] = ((hap_len - hap_to_assign_len, 
                                   hap_len - hap_to_assign_len + current_block_SNP_end - current_block_SNP_start),
                                  (current_block_SNP_start, current_block_SNP_end))
            hap_to_assign_len = hap_to_assign_len - (current_block_SNP_end - current_block_SNP_start)  
            current_block = current_block + 1    
    
    return res

def find_ld_blocks(line, bins_df):
    logger_find_ld_blocks = logging.getLogger('find_ld_blocks')
    logger_find_ld_blocks.debug('passing {} {}'.format(line[6], line[7]))
    return find_block_intervals_from_hap(bins_df, int(line[6]), int(line[7]))

def find_read_specific_error_rate(line):
    return float(line[5])

def find_observation(line):
    obs1d = np.array([1 - int(x) for x in list(line[-1])], dtype=np.int8)
    return obs1d.reshape((1, len(obs1d)))

@jit
def compute_likelihood_from_a_read(line, bins_df, hapkey):
    """Given one line of the input hap file, compute and update the log-likelihood
    """
    
    logger_read_ll = logging.getLogger('compute_likelihood_from_a_read')
    
    logger_read_ll.debug('line = {}'.format(line))  
    
    obs = find_observation(line)
    err_rate = find_read_specific_error_rate(line)
    
    likelihood_of_read = dict([])
    
    logger_read_ll.debug('parsed ld blocks are {}'.format(str(find_ld_blocks(line, bins_df))))
    
    for block_id, v in find_ld_blocks(line, bins_df).items():
        logger_read_ll.debug('The LD we are working on is: {} {}'.format(block_id, v))
                
        obs_s, obs_e = v[0]
        block_s, block_e = v[1]    

        n_mismatch = np.sum(
            hapkey[block_id][:, block_s:block_e] != obs[:, obs_s:obs_e], axis = 1
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

def update_ll(hap_f, bins_df, hapkey, log_likelihood = dict([])):

    logger_ll = logging.getLogger('update_ll')
        
    with open(hap_f, 'r') as f:
        for line_num, raw_line in enumerate(f):
            line  = raw_line.rstrip().split('\t')
            
            if(line_num % 100 == 0):            
#            if True:
                logger_ll.info(
                    'processing read {} {}'.format(line_num, line[3])
                )                                

            likelihood_of_read = compute_likelihood_from_a_read(line, bins_df, hapkey)            
            log_likelihood = update_likelihood_dict(likelihood_of_read, log_likelihood)
    return log_likelihood


def log_likelihood_main(hap_f, hapkey_dir, bins_df, log_likelihood_dir):
    hapkey = read_hapkey(hapkey_dir, bins_df)
    log_likelihood = update_ll(hap_f, bins_df, hapkey)
    
    for block_id, ll in log_likelihood.items():
        np.savez(
            '{}/{}.npz'.format(log_likelihood_dir, block_id),
            log_likelihood = ll
        )          
        
    
def main():    
    logger_main = logging.getLogger('main')    
    
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=_README_
    )
        
    parser.add_argument('--hap', metavar='r', required=True,
                        help='[input] hap file (from reads, use bam_to_hap to generate hap file from bam file)')
    parser.add_argument('--key', metavar='k', required=True,
                        help='[input] haplotype keys dir') 
    parser.add_argument('--bins', metavar='b', required=True,
                        help='[input] bins.tsv file')
    parser.add_argument('--ll', metavar='l', required=True,
                        help='[output] log-likelihood dir')
    
    args = parser.parse_args()

    # check if the hap file exists
    if not os.path.isfile(args.hap):
        raise FileNotFoundError(
            errno.ENOENT, os.strerror(errno.ENOENT), args.hap
        )        

    logger_main.info('running log_likelihood')
    logger_main.info('  hap  : {}'.format(args.hap))
    logger_main.info('  keys : {}'.format(args.key))    
    logger_main.info('  bins : {}'.format(args.bins))
    logger_main.info('  ll   : {}'.format(args.ll))    
                
    # read the bin definition tsv file
    bins_df = pd.read_csv(args.bins, sep='\t')
                          
    # create the output directories if not exist    
    if not os.path.exists(args.ll):
        os.makedirs(args.ll)
    
    log_likelihood_main(
        hap_f = args.hap, 
        hapkey_dir = args.key, 
        bins_df = bins_df,
        log_likelihood_dir = args.ll
    )


if __name__ == "__main__":
    main()
    
