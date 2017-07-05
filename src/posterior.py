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
        'level': logging.INFO,
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

def read_log_likelihood(block_df, log_likelihood_dir):
    
    logger_ll = logging.getLogger('read_log_likelihood')
    logger_ll.info(
        'reading read_log_likelihood from {}'.format(log_likelihood_dir)
    )        
    
    skipped_blocks = []
    
    log_likelihood = dict([])
    for block_id in range(len(block_df)):
        if(block_id % 100 == 0):
            logger_ll.info(
                'reading block {} of {}'.format(block_id, len(block_df))
            )    
        npz_file = '{}/ll{}.npz'.format(log_likelihood_dir, block_id)
        if not os.path.isfile(npz_file):
            skipped_blocks.append(block_id)
        else:
            ll = np.load(npz_file)
            log_likelihood[block_id] = ll['ll']

    if(len(skipped_blocks) > 0):
        logger_ll.info(
            'skipped blocks are {}'.format(skipped_blocks)
        )
        
    
    logger_ll.info(
        'log likelihood is loaded on memory'
    )    
        
    return log_likelihood

def compute_log_posterior_sub(prior_cnt_vals, log_likelihood, block_id):
    prior_log_prob = np.log(
        1.0 * prior_cnt_vals[block_id] / np.sum(prior_cnt_vals[block_id])
    )

    log_joint_prob = prior_log_prob + log_likelihood[block_id]

    log_partition = np.log(np.sum(np.exp(log_joint_prob)))

    log_posterior = log_joint_prob - log_partition    

    return log_posterior

def compute_log_posterior(prior_cnt_vals, log_likelihood):
    logger_clp = logging.getLogger('compute_log_posterior')
    logger_clp.info(
        'computing log posterior probabilities'
    )            
    log_posterior = dict([])
    for block_id in sorted(log_likelihood.keys()):
        if(block_id % 100 == 0):
            logger_clp.info(
                'processing block {}'.format(block_id)
            )            
        log_posterior[block_id] = compute_log_posterior_sub(prior_cnt_vals, log_likelihood, block_id)

    logger_clp.info(
        'posterior computation done!'
    )    
        
    return log_posterior

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
    posterior_dir='{}/{}/{}'.format(data_dir, 'posterior', hap_basename[:-4])
    if not os.path.exists(posterior_dir):
        os.makedirs(posterior_dir)
    
    block_df = pd.read_csv(block_tsv_f, sep='\t')

    prior_cnt_keys, prior_cnt_vals = read_prior_cnts(block_df, prior_count_dir)

    log_likelihood = read_log_likelihood(block_df, log_likelihood_dir)

    log_posterior = compute_log_posterior(prior_cnt_vals, log_likelihood)    
    
    for block_id, lp in log_posterior.items():
        np.savez(
            '{}/{}.npz'.format(posterior_dir, block_id), 
            lp = lp
        )

if __name__ == "__main__":
    main()
    
