from __future__ import print_function

_README_ = '''
-------------------------------------------------------------------------
Compute the log posterior of haplotypes given prior counts and log-likelihood

Author: Yosuke Tanigawa (ytanigaw@stanford.edu)
Date: 2017/07/07
-------------------------------------------------------------------------
'''

import sys
import os
import argparse
import logging
from logging.config import dictConfig
import numpy as np
import pandas as pd
import itertools as it
import collections as cl
import bisect
from scipy.misc import logsumexp

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


def compute_log_posterior_sub(hapcnt, log_likelihood, block_id):
    prior_log_prob = np.log(
        1.0 * hapcnt[block_id] / np.sum(hapcnt[block_id])
    )

    log_joint_prob = prior_log_prob + log_likelihood[block_id]
    
    log_partition = logsumexp(log_joint_prob)

    log_posterior = log_joint_prob - log_partition    

    return log_posterior

def compute_log_posterior(hapcnt, log_likelihood):
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
        log_posterior[block_id] = compute_log_posterior_sub(hapcnt, log_likelihood, block_id)

    logger_clp.info(
        'posterior computation done!'
    )    
        
    return log_posterior


def log_posterior_main(log_likelihood_dir, hapcnt_dir, bins_df, log_posterior_dir):
    hapcnt = read_hapcnt(hapcnt_dir, bins_df)
    log_likelihood = read_log_likelihood(log_likelihood_dir, bins_df)
    
    log_posterior = compute_log_posterior(hapcnt, log_likelihood)    
    
    for block_id, lp in log_posterior.items():
        np.savez(
            '{}/{}.npz'.format(log_posterior_dir, block_id), 
            log_posterior = lp
        )
        
    
def main():    
    logger_main = logging.getLogger('main')    
    
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=_README_
    )
    
    parser.add_argument('--ll', metavar='l', required=True,
                        help='[input] log-likelihood dir')        
    parser.add_argument('--cnt', metavar='c', required=True,
                        help='[input] haplotype counts dir')
    parser.add_argument('--bins', metavar='b', required=True,
                        help='[input] bins.tsv file')
    parser.add_argument('--lp', metavar='p', required=True,
                        help='[output] log-posterior dir')
    
    args = parser.parse_args()

    logger_main.info('running log_posterior')
    logger_main.info('  ll   : {}'.format(args.ll))
    logger_main.info('  cnt  : {}'.format(args.cnt))    
    logger_main.info('  bins : {}'.format(args.bins))
    logger_main.info('  lp   : {}'.format(args.lp))    
                
    # read the bin definition tsv file
    bins_df = pd.read_csv(args.bins, sep='\t')
                          
    # create the output directories if not exist    
    if not os.path.exists(args.lp):
        os.makedirs(args.lp)
    
    log_posterior_main(
        log_likelihood_dir = args.ll, 
        hapcnt_dir = args.cnt,
        bins_df = bins_df,
        log_posterior_dir = args.lp
    )

    
if __name__ == "__main__":
    main()
    
