_README_ = '''
-------------------------------------------------------------------------
Convert BAM file to hap file

hap file is a tab-delimited file with the following fields:
1. chrom - The name of the chromosome (e.g. chr3, chrY, chr2_random)
2. chromStart - The starting position of the feature in the chromosome or scaffold. The first base in a chromosome is numbered 0.
3. chromEnd - The ending position of the feature in the chromosome or scaffold. The chromEnd base is not included in the display of the feature.
4. name - Defines the name of the line.
5. MAPQ - MAPping Quality.
6. ReadErr - Read specific error rate. The ratio of mismatches to the reference sequence in the non-polymorphic sites
7. bimStart - Starting position of the read as an index on the given bim file
8. bimEnd - Ending position of the read as an index on the given bim file
9. HAP - haplotype string in the interval [bimStart, bimEnd). 1 indicates major allele and 0 indicates minor allele.

Author: Yosuke Tanigawa (ytanigaw@stanford.edu)
Date: 2017/04/09
-------------------------------------------------------------------------
'''

import sys
import argparse
import numpy as np

import pysam

from Position import Position
from FastaFile import FastaFile
from BimFile import BimFile
from Read import Read
from Nucleotide import Nucleotide


def filter_mismatches(read, bim, quality_threshod=14):
    """
    Given a read and a bim file (and quality threshold),
    identify list of SNPs
      - is a mismatch
      - base call quality value >= quality threshold
      - the position is on bim file
      - read is minor allele in the bim file
      - reference is major allele in the bim file
    """
    # Corresponding id interval (on bim file) is
    read_bim_l, read_bim_r = bim.find_index_interval(query_l = read.l, query_r = read.r)
    # Get polymorphic SNP positions (written on bim)
    read_snp_pos = bim.get_bp(read.chr)[read_bim_l : read_bim_r]
    # Get a candidate of SNPs
    #  - Is a mismatch
    #  - Base call quality value >= quality_threshod
    #  - The position is on bim file
    read_poly_mm = read.get_mismatches_on_polymorphic_sites(read_snp_pos, quality_threshod)
    # take the position of candidate SNPs
    read_poly_mm_pos = np.array([x.reference_position for x in read_poly_mm])
    # convert to bim indeces
    read_poly_mm_bim_id = bim.find_index_list([Position(read.chr, ref_pos) for ref_pos in read_poly_mm_pos])
    # minor allele on bim file
    read_bim_a1 = [Nucleotide[bim.get_allele_1(read.chr)[id]] for id in read_poly_mm_bim_id]
    # major allele on bim file
    read_bim_a2 = [Nucleotide[bim.get_allele_2(read.chr)[id]] for id in read_poly_mm_bim_id]
    # return list of mismatches
    return [read_poly_mm[i] for i in range(len(read_poly_mm)) if
            read_poly_mm[i].reference == read_bim_a2[i] and
            read_poly_mm[i].read == read_bim_a1[i]]

def find_haplotype(read, read_poly_mm, bim):
    """Given list of polymorphic SNPs on a read and corresponding interval on bim file,
    construct a haplotype on the read
    """
    # bim id interval (on bim file)
    read_bim_l, read_bim_r = bim.find_index_interval(query_l=read.l, query_r=read.r)
    read_poly_mm_pos = np.array([x.reference_position for x in read_poly_mm])
    read_poly_mm_bim_id = bim.find_index_list([Position(read.chr, ref_pos)
                                               for ref_pos in read_poly_mm_pos])
    read_hap = np.array([0 if i in set(read_poly_mm_bim_id) else 1
                         for i in range(read_bim_l, read_bim_r)])
    return read_hap


def process_read(read, bim, quality_threshod=14):
    """Given a read, find bim interval and haplotype representation
    """
    # bim id interval (on bim file)
    read_bim_l, read_bim_r = bim.find_index_interval(query_l=read.l, query_r=read.r)
    # polymorphic SNPs on reads
    read_poly_mm = filter_mismatches(read, bim, quality_threshod)
    # Construct the haplotype of the read
    read_hap = find_haplotype(read, read_poly_mm, bim)
    # convert to string representation
    read_hap_str = ''.join([str(x) for x in read_hap])
    # compute read specific error rate
    read_non_poly_mismatch = read.n_mismatch() - len(read_poly_mm)
    # read specific error rate with Jeffreys prior
    read_error_rate = 1.0 * (1 + read_non_poly_mismatch) / (1 + read.query_alignment_length - (read_bim_r - read_bim_l))
    return (read_bim_l, read_bim_r, "{:e}".format(read_error_rate), read_hap_str)


def bamToHap(fasta, bim, bam, quality_threshod = 14):
    for x in bam:
        r = Read(x, fasta)
        try:
            read_bim_l, read_bim_r, read_error_rate, read_hap_str = process_read(r, bim, quality_threshod)
        except KeyError as e:
            print >> sys.stderr, r.query_name, "KeyError", e
        except Exception as e:
            print >> sys.stderr, r.query_name, type(e), e.args, e
        else:
            # name, chr, start, end, MAPQ, bim_l, bim_r, hapstr
            read_output = (
            r.reference_name, r.l.get_pos(), r.r.get_pos(),
            r.query_name,
            r.mapping_quality, read_error_rate,
            read_bim_l, read_bim_r,
            read_hap_str)
            print '\t'.join([str(x) for x in read_output])


def bamToHap_main():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=_README_)

    parser.add_argument('--fasta_ext', metavar='F',
                        default='fa',
                        help='extension of fasta file (default: fa)')

    parser.add_argument('-q', metavar='q', type=int,
                        default=14,
                        help='quality threshold (default: 14)')

    parser.add_argument('--fasta', metavar='f', required=True,
                        help='path for fasta file')

    parser.add_argument('--bim', metavar='b', required=True,
                        help='path for bim file')

    parser.add_argument('bam', metavar='r',
                        help='path for bam file (read)')

    args = parser.parse_args()

    fasta = FastaFile(args.fasta, fast_ext=args.fasta_ext)
    bim = BimFile(args.bim)
    bam = pysam.AlignmentFile(filename=args.bam, mode='rb')

    bamToHap(fasta, bim, bam, quality_threshod = args.q)


if __name__ == '__main__':
    bamToHap_main()
