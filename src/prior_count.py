from __future__ import print_function

import sys
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


def read_alleles_block(pgen_f, block_df, block_id):
    """wrapper function of pgen.PgenReader.read_alleles_range for a LD block"""
    bim_interval = block_df.bim_interval[block_id]
    with pg.PgenReader(pgen_f) as pgr:
        buf_ndary = np.zeros(
            (bim_interval[1] - bim_interval[0], pgr.get_raw_sample_ct() * 2), 
            dtype=np.int32
        )
        pgr.read_alleles_range(bim_interval[0], bim_interval[1], buf_ndary)
    return buf_ndary


def filter_missing_alleles(geno):
    return geno[:, np.sum(geno == -9, axis = 0) == 0]


def count_alleles_freq(geno):
    tuples = [tuple(1 - geno[:, hap])
              for hap in range(geno.shape[1])]
    return cl.Counter(tuples)


def prior_count(pgen_f, block_df, block_id):
    logger_cnt = logging.getLogger('prior_count')
    cnt = count_alleles_freq(
        filter_missing_alleles(
            read_alleles_block(pgen_f, block_df, block_id)
        )
    )    
    logger_cnt.info(
        'block {:3d} {} bytes'.format(block_id, sys.getsizeof(cnt))
    )    
    return cnt

def main():
    data_dir='/oak/stanford/groups/mrivas/users/ytanigaw/nanopore-data/'
    geno_bed_log='chr20impv1-keep-maf0005-snv-biallelic-geno01-hwe1e-10.log'
    block_log='chr20impv1-keep-maf0005-snv-biallelic-geno01-hwe1e-10-block-stronglow050-stronghigh083-infofrac10.log'
    pgen_log='chr20impv1-keep-maf0005-snv-biallelic-geno01-hwe1e-10-pg.log'
    pgen_f = '{}{}.pgen'.format(data_dir, pgen_log[:-4])
    block_bed_f = '{}{}.bed'.format(data_dir, block_log[:-4])
    
    prior_dir='{}prior_count'.format(data_dir)    
    
    block_bed = pd.read_csv(block_bed_f, sep='\t', names=['chrom', 'chromStart', 'chromEnd', 'name'])
    block_bed['bim_interval'] = block_bed.name.map(lambda x: [int(pos) for pos in x.split(':')])

    for block_id in range(len(block_bed)):
        cnt = prior_count(pgen_f, block_bed, block_id)
        np.savez(
            '{}/{}.npz'.format(prior_dir, block_id),
            keys = np.array(cnt.keys(), dtype=np.int8),
            vals = np.array(cnt.values())
        )                
    
if __name__ == '__main__':
    main()
    
    