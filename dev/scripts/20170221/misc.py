import numpy as np
import itertools as it
import logging
import sys
import re
import collections
import enum

@enum.unique
class Nucleotide(enum.Enum):
    N = -1
    A = 0
    C = 1
    G = 2
    T = 3
    
#class Mismatch():
#    def __init__(self, reference_position, reference, read, quality):
#        self.reference_position = reference_position
#        self.reference = reference
#        self.read      = read
#        self.quality   = quality
Mismatch = collections.namedtuple('Mismatch', 'reference_position reference read quality')   


def get_ref_file_name(filename_template, chromosome, extension):
    '''Convert filename template to actual file name
    
    filename_template: full path to reference file with meta variables${CHR} ${EXT}
      ${CHR} will be replaced with chromosome
      ${EXT} will be replaced with extension (fa, pgen)
    chromosome:        chromosome
    extension:         extension
    '''    
    return(re.sub(r'\${CHR}(.*)\${EXT}', r'{}\1{}'.format(str(chromosome),
                                                          str(extension)),
                  filename_template))

def get_fasta_file_name(param_obj, chromosome):
    return(get_ref_file_name(param_obj['ref_genome_template'], 
                             chromosome, 
                             param_obj['fasta_ext']))

def get_pgen_file_name(param_obj, chromosome):
    return(get_ref_file_name(param_obj['ref_population_template'], 
                             chromosome, 
                             'pgen'))

def get_bim_file_name(param_obj, chromosome):
    return(get_ref_file_name(param_obj['ref_population_template'], 
                             chromosome, 
                             'bim'))
