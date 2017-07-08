import numpy as np
import os
import logging
from logging.config import dictConfig

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

def read_generic(npz_dir, bins_df, npz_key = None,
                 logger_name = 'read_generic'):
    
    variable_name = 'npz' if (npz_key is None) else npz_key
    
    logger_this = logging.getLogger(logger_name)    
    logger_this.info(
        'reading {} from {}'.format(variable_name, npz_dir)
    )        
    
    data = dict([])
    skipped_blocks = set([])
    
    for block_id in range(len(bins_df)):
        
        if(block_id % 100 == 0):
            logger_this.info(
                'reading block {} of {}'.format(block_id, len(bins_df))
            )    
            
        npz_file = '{}/{}.npz'.format(npz_dir, block_id)
        
        if not os.path.isfile(npz_file):
            skipped_blocks.add(block_id)
        else:                      
            data_point = np.load(npz_file)

            if(npz_key is None):
                npz_key = data_point.files[0]

            data[block_id] = data_point[npz_key]
    
    if(len(skipped_blocks) > 0):
        logger_this.info(
            'skipped blocks are {}'.format(sorted(skipped_blocks))
        )        
        
        
    logger_this.info(
        '{} is loaded on memory'.format(variable_name)
    )    
        
    return data

def read_hapkey(hapkey_dir, bins_df):
    return read_generic(npz_dir = hapkey_dir, bins_df = bins_df, 
                        logger_name = 'read_hapkey')

def read_hapcnt(hapcnt_dir, bins_df):
    return read_generic(npz_dir = hapcnt_dir, bins_df = bins_df, 
                        logger_name = 'read_hapcnt')

def read_log_likelihood(log_likelihood_dir, bins_df):
    return read_generic(npz_dir = log_likelihood_dir, bins_df = bins_df, 
                        logger_name = 'read_log_likelihood')

def read_log_posterior(log_posterior_dir, bins_df):
    return read_generic(npz_dir = log_posterior_dir, bins_df = bins_df, 
                        logger_name = 'read_log_posterior')
