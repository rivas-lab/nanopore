from __future__ import print_function

_README_ = '''
-------------------------------------------------------------------------
Find the prior count of the population reference

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
import multiprocessing

import pgenlib as pg

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


def read_alleles_block(pgen_f, bins_df, block_id):
    """wrapper function of pgenlib.PgenReader.read_alleles_range for a LD block"""    
    bim_s = bins_df.bimIdStart[block_id]
    bim_e = bins_df.bimIdEnd[block_id]
    
    with pg.PgenReader(pgen_f) as pgr:
        buf_ndary = np.zeros(
            (bim_e - bim_s, pgr.get_raw_sample_ct() * 2), 
            dtype=np.int32
        )
        pgr.read_alleles_range(bim_s, bim_e, buf_ndary)
        
    return buf_ndary


def filter_missing_alleles(geno):
    return geno[:, np.sum(geno == -9, axis = 0) == 0]


def count_alleles_freq(geno):
    tuples = [tuple(geno[:, hap]) for hap in range(geno.shape[1])]
    return cl.Counter(tuples)


def prior_count_bin(pgen_f, bins_df, block_id):
    geno = filter_missing_alleles(
        read_alleles_block(pgen_f, bins_df, block_id)
    )        
    return count_alleles_freq(geno), geno.shape


def prior_count_main(pgen_f, bins_df, hapkey_dir, prior_count_dir):
    logger_cnt_m = logging.getLogger('prior_count_main')
    
    logger_cnt_m.info(
        '\t'.join(['block_id', 'shape(raw_genotypes)', 
                   'alleles_without_missing_data', 'num_of_uniq_alleles'])
    )
    
    for block_id in range(len(bins_df)):
        cnt, geno_shape = prior_count_bin(pgen_f, bins_df, block_id)        
        logger_cnt_m.info(
            '{:3d}\t{}\t{}\t{}'.format(block_id, geno_shape, sum(cnt.values()), len(cnt))
        )        
        
        np.savez(
            '{}/{}.npz'.format(hapkey_dir, block_id),
            hapkey = np.array(cnt.keys(), dtype=np.int8),
        )                
        
        np.savez(
            '{}/{}.npz'.format(prior_count_dir, block_id),
            prior_count = np.array(cnt.values())
        )                    
    

def main():    
    logger_main = logging.getLogger('main')    
    
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=_README_
    )
        
    parser.add_argument('--pgen', metavar='p', required=True,
                        help='[input] pgen file')
    parser.add_argument('--bins', metavar='b', required=True,
                        help='[input] bins.tsv file')
    parser.add_argument('--key', metavar='k', required=True,
                        help='[output] haplotype keys dir')
    parser.add_argument('--cnt', metavar='c', required=True,
                        help='[output] haplotype counts dir')
    
    args = parser.parse_args()

    # check if the pgen file exists
    if not os.path.isfile(args.pgen):
        raise FileNotFoundError(
            errno.ENOENT, os.strerror(errno.ENOENT), args.pgen
        )        

    logger_main.info('running prior_count')
    logger_main.info('  pgen : {}'.format(args.pgen))
    logger_main.info('  bins : {}'.format(args.bins))
    logger_main.info('  keys : {}'.format(args.key))
    logger_main.info('  cnts : {}'.format(args.cnt))    
                
    # read the bin definition tsv file
    bins_df = pd.read_csv(args.bins, sep='\t')
                          
    # create the output directories if not exist
    for outDir in [args.key, args.cnt]:
        if not os.path.exists(outDir):
            os.makedirs(outDir)
            
    prior_count_main(args.pgen, bins_df, args.key, args.cnt)

    
if __name__ == '__main__':
    main()
