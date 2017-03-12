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

class read:
    '''class to manage mapped read    
    '''
    def __init__(self, aligned_segment, reference):
        # This copy will be removed
        self.aligned_segment = aligned_segment
        
        # find all mismatches
        self.find_mismatches(aligned_segment, reference)
        
        # copy useful info
        self.query_name      = aligned_segment.query_name     
        self.reference_name  = aligned_segment.reference_name
        self.reference_start = aligned_segment.reference_start
        self.reference_end   = aligned_segment.reference_end
        self.length          = aligned_segment.reference_length
        self.mapping_quality = aligned_segment.mapping_quality


    def get_mismatches(self, quality_threshod = 14):
        return([i for i in self.mismatches 
                if i.quality >= quality_threshod and
                   i.reference is not Nucleotide['N']])

    def n_mismatch(self, quality_threshod = 14):
        return(len(self.get_mismatches(quality_threshod)))
    
    def n_match(self, quality_threshod = 14):
        return(self.length - self.n_mismatch(quality_threshod))
    
        
    def find_mismatches(self, aligned_segment, reference):
        '''find all mismatches
        [args]
        pysam.AlignedSegment aligned_segment: mapped fragment
        pysam.FastxFile      reference:       reference sequence
        '''
        aligned_pairs = np.array(aligned_segment.get_aligned_pairs(matches_only=True, with_seq=True))
        
        # fetch corresponding reference sequence
        reference_str = reference.fetch(reference = aligned_segment.reference_name,
                                        start     = aligned_segment.reference_start, 
                                        end       = aligned_segment.reference_end).upper()
        
        # obtain nucleotide letters on both read and reference
        read_letters = np.array([aligned_segment.query_sequence[int(read_position)].upper()
                                 for read_position in aligned_pairs[:, 0]])
        ref_letters  = np.array([reference_str[int(ref_position) - 
                                               aligned_segment.reference_start] 
                                 for ref_position in aligned_pairs[:, 1]])

        # enumerate all the mismatches by comparing nucleotide letters
        self.mismatches = [Mismatch(reference_position = \
                                      int(aligned_pairs[mismatch_pos_on_pairs][1]),
                                    reference = \
                                      Nucleotide[ref_letters[mismatch_pos_on_pairs]],
                                    read      = \
                                      Nucleotide[read_letters[mismatch_pos_on_pairs]],
                                    quality   = \
                                      aligned_segment.query_qualities[int(aligned_pairs[mismatch_pos_on_pairs][0])])
                           for mismatch_pos_on_pairs 
                           in np.where(read_letters != ref_letters)[0]]

    def __str__(self):
        return('\t'.join(['{}:{}-{}'.format(self.reference_name, 
                                            self.reference_start,
                                            self.reference_end),
                          str(self.n_match()),
                          str(self.n_mismatch()),
                          str(self.mapping_quality),
                          self.query_name]))