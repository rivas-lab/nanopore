from __future__ import print_function

_README_ = '''
-------------------------------------------------------------------------
Write the results of the inference to a pgen file

Author: Yosuke Tanigawa (ytanigaw@stanford.edu)
Date: 2017/07/11
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

import pgenlib as pg

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
        'level': logging.DEBUG,
        },
)
dictConfig(logging_config)


def argmax2(ary):
    return (-ary).argpartition(2)[:2]


def find_hap_index_block(log_posterior_block, homo_threshold_log = np.log(0.9)):
    logger_find_hap_index_block = logging.getLogger('find_hap_index_block')

    top2 = argmax2(log_posterior_block)

    logger_find_hap_index_block.debug([np.exp(x) for x in log_posterior_block[top2]])

    if(log_posterior_block[top2[0]] >= homo_threshold_log):
        return np.array([top2[0]] * 2)
    else:
        return top2

    
def find_hap_block(hapkey, block_id, log_posterior, homo_threshold_log = np.log(0.9)):
    return np.array(
        hapkey[block_id][
            find_hap_index_block(log_posterior[block_id], homo_threshold_log)
        ].transpose(),
        dtype=np.int32
    )    

def pgenout_main(bins_df, hapkey_dir, log_posterior_dir, pgen_out, homo_threshold = 0.9):
    hapkey = read_hapkey(hapkey_dir, bins_df)
    log_posterior = read_log_posterior(log_posterior_dir, bins_df)
    
    with pg.PgenWriter(pgen_out, 1, bins_df.nSNPs.sum(), True) as pgw:
        for block_id in range(len(bins_df)):
            pgw.append_alleles_batch(
                find_hap_block(
                    hapkey, block_id, log_posterior, 
                    homo_threshold_log = np.log(homo_threshold)
                ).copy(order='C')
            )
    
def main():    
    logger_main = logging.getLogger('main')    
    
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=_README_
    )
    
    parser.add_argument('--bins', metavar='b', required=True,
                        help='[input] bins.tsv file')        
    parser.add_argument('--key', metavar='k', required=True,
                        help='[input] haplotype keys dir') 
    parser.add_argument('--lp', metavar='p', required=True,
                        help='[input] log-posterior dir')
    parser.add_argument('--out', metavar='o', required=True,
                        help='[output] output pgen file')
    parser.add_argument('--homo', metavar='H', type=float, default=0.9,
                        help='posterior probability threshold for homozygosity')

    
    args = parser.parse_args()

    logger_main.info('running log_likelihood')
    logger_main.info('  bins : {}'.format(args.bins))    
    logger_main.info('  keys : {}'.format(args.key))    
    logger_main.info('  lp   : {}'.format(args.lp))    
    logger_main.info('  homo : {}'.format(args.homo))
    logger_main.info('  out  : {}'.format(args.out))
                
    # read the bin definition tsv file
    bins_df = pd.read_csv(args.bins, sep='\t')
                                  
    pgenout_main(
        bins_df = bins_df, 
        hapkey_dir = args.key, 
        log_posterior_dir = args.lp,
        pgen_out = args.out, 
        homo_threshold = args.homo
    )


if __name__ == "__main__":
    main()
    
