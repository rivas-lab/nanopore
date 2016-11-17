#!/usr/bin/python

"""
Read BAM files to obtain coverage information for selected regions.

Author: Peter Humburg (2011) based on original code by Manuel Rivas
"""

###################################################################################################

import array
import sys,re
import os
import pysam
import SAMpileuphelper
import time
import fastafile
import bamfileutils
import chaplotype
import variant
import samtoolsWrapper
import logging

from optparse import OptionParser
from SAMpileuphelper import ascii_list

from multiprocessing import Pool

## set up logging
logger = logging.getLogger('Log')

###################################################################################################

def requiredArgs(options):
    return ['pif', 'tgf', 'ref', 'mqthr', 'bqthr', 'outputdir', 'ncpu', 'mode','is_paired','rmdup','getindels']

###################################################################################################

def inputFiles(options):
    return ['pif', 'tgf', 'ref']

###################################################################################################

class Pileup:
    """
    Class to handle fetch callbacks from samtools.
    This will accumulate information about all alignments in a given region,
    counting the number of times each base is observed at each position in the region.
    """
    cover = None
    start = -1
    end = -1
    mqthr = 1
    bqthr = 22
    base = {'A':0, 'C':1, 'G':2, 'T':3, 'a':0, 'c':1, 'g':2, 't':3}
    is_exome = False
    rmdup = False
    paired = False

    def __init__(self, start, end, mqthr, bqthr, mode,paired,rmdup):
        self.start = int(start)
        self.end = int(end)
        self.mqthr = int(mqthr)
        self.bqthr = int(bqthr)
        self.cover = [[0 for col in range(12)] for row in range(end - start + 1)]
        if mode == 'exome': self.is_exome = True
        if paired: self.paired = True
        if rmdup: self.rmdup = True

    def __getitem__(self, key):
        return self.cover[key]

    def __len__(self):
        return len(self.cover)

    def __call__(self, align):
        ## skip low quality alignments
        if align.mapq < self.mqthr: return

        ## skip marked duplicates when analysing exome data
        if align.is_qcfail: return
        if self.paired and not align.is_proper_pair: return
        if self.rmdup and align.is_duplicate: return
        if self.is_exome and align.is_duplicate: return

        ## determine index of first read base that lies inside target region
        ## in read and region coordinates
        offset = align.pos - self.start
        qIdx = max(0, -offset)
        rIdx = max(0, offset)

        ## parse cigar string
        cigar = align.cigar
        cigar.reverse()
        op = cigar.pop()
        cigarLen = op[1]


        ## discard hard clipped portion of read

        while op[0] == 5 or op[0] == 4:
           op = cigar.pop()
           cigarLen = op[1]

           
       ## find first CIGAR operation that overlaps with target region

        while cigarLen < qIdx:
            ## insertion before start of target region
            if op[0] == 1:
                qIdx += op[1]
            ## deletion before start of target region
            if op[0] == 2:
                qIdx -= op[1]
                cigarLen -= op[1]
            op = cigar.pop()
            cigarLen += op[1]

        cIdx = qIdx - cigarLen + op[1]
        while True:
            ## process bases covered by current CIGAR operation
            ## alignment match
            if op[0] == 0:
                qual = ascii_list(align.qqual)
                seq = align.query
                while cIdx < op[1] and qIdx < len(seq) and rIdx < len(self.cover):
                    if qual[qIdx] >= self.bqthr and seq[qIdx] != "N":
                        idx = self.base[seq[qIdx]]
                        if align.is_reverse:
                            idx += 6
                        self.cover[rIdx][idx] += 1
                    qIdx += 1
                    rIdx += 1
                    cIdx += 1
            ## insertion
            if op[0] == 1:
                seq = align.query
                idx = 5
                if align.is_reverse:
                    idx += 6
                self.cover[rIdx][idx] += 1
                qIdx += op[1]
            ## deletion
            elif op[0] == 2:
                seq = align.query
                idx = 4
                if align.is_reverse:
                    idx += 6
                self.cover[rIdx][idx] += 1
                rIdx += op[1]
                if(qIdx == -offset):
                    qIdx -= 1

            ## skipped region in reference
            elif op[0] == 3:
                qIdx -= op[1]

            if len(cigar) > 0 and rIdx < len(self.cover):
                op = cigar.pop()
                cIdx = 0
            else:
                break


###################################################################################################

def parsePIF(file):
    """
    Read the pool info file and return list of BAM file names.
    """
    poolfile = open(file,'r')
    pools = poolfile.readlines()
    poolfile.close()
    bamfiles = []
    for line in pools[1:]:
        line = line.rstrip()
        line = line.split()
        bamfiles.append(line[0])
    return bamfiles

###################################################################################################

def parseTGF(file):
    """
    Read the target information file and return list with chromosome, start and end positions.
    Target regions are sorted internally and overlapping regions are merged.
    Raises a ValueError if TGF file exists but is empty.
    """
    target = open(file, 'r').readlines()
    ## check that at least one region is defined
    if len(target) <= 1: raise ValueError('No target regions defined in TGF file')
    
    ## ensure target regions are sorted
    target = target[1:]
    for i in range(len(target)):
        target[i] = target[i].rstrip()
        target[i] = target[i].split()
    ## sort regions by chromosome name, start and end
    target.sort(key=lambda region: int(region[3]))
    target.sort(key=lambda region: int(region[2]))
    target.sort(key=lambda region: region[1])
    
    chrom = [target[0][1]]
    start = [int(target[0][2])]
    end = [int(target[0][3])]
        
    for line in target[1:]:   
        ## check for overlapping regions
        if line[1] != chrom[-1] or int(line[2]) > end[-1]:
            chrom.append(line[1])
            start.append(int(line[2]) -1)
            end.append(int(line[3]))
        else:
            end[-1] = int(line[3])
    return chrom,start,end

###################################################################################################

def formatOutput(chrom, pos, ref, cover, indelList):
    """
    Format a single line of output.
    """
    idx = {'A':0, 'C':1, 'G':2, 'T':3, 'a':0, 'c':1, 'g':2, 't':3}
    bases = ['A', 'C', 'G', 'T']

    if ref == 'N':
        coverSum = [ cover[i] + cover[i+6] for i in range(4)]
        ref = bases[coverSum.index(max(coverSum))]

    out = chrom + ' ' + str(pos) + ' ' + chrom + ':' + str(pos) + ' ' + ref + ' ' + str(idx[ref])
    
    for i in xrange(len(cover)):
        # dFwd
        if i == 4:

            out += ' '
            count = 0
            if indelList is not None:
                for indel in indelList:
                    if len(indel.removed) > 0:
                        count += indel.nFwdReads
            out += str(count)

        # iFwd
        elif i == 5:

            out += ' '
            count = 0
            if indelList is not None:
                for indel in indelList:
                    if len(indel.added) > 0:
                        count += indel.nFwdReads
            out += str(count)

        # dRev
        elif i == 10:

            out += ' '
            count = 0
            if indelList is not None:
                for indel in indelList:
                    if len(indel.removed) > 0:
                        count += indel.nRevReads
            out += str(count)

        # iRev
        elif i == 11:

            out += ' '
            count = 0
            if indelList is not None:
                for indel in indelList:
                    if len(indel.added) > 0:
                        count += indel.nRevReads
            out += str(count)

        else:
            out += ' ' + str(cover[i])

    # Now write out the insertion/deletion sequences in the same order as the coverage 
    # values
    out += ' '

    if indelList is not None:
        for indel in indelList:
            if len(indel.removed) > 0:
                out += str(indel.removed) + ":" + str(indel.nFwdReads) + "," + str(indel.nRevReads) + ";"
        if out[-1] == ",":
            out = out[:-1]
    if(out[-1] == ";"):
        out = out[:-1]

    out += ' '

    if indelList is not None:
        for indel in indelList:
            if len(indel.added) > 0:
                out += str(indel.added) + ":" + str(indel.nFwdReads) + "," + str(indel.nRevReads) + ";"
        if out[-1] == ";":
            out = out[:-1]
    if(out[-1] == ";"):
        out = out[:-1]

    out += '\n'
    return out

###################################################################################################

def getIndelsByPosition(chrom, start, end, bamFile, minMapQual, minReadQual, ref):
    """
    Return a hash of lists of indel candidates, keyed by position. For each key/position, store
    a list of all indels seen at that position, after left-normalisation.
    """
    indels = {}
    refFile = fastafile.FastaFile(ref, ref + ".fai")
    rawIndelCandidates = variant.VariantCandidateGenerator((chrom,start,end), [bamFile], refFile, int(minMapQual), 3, int(minReadQual), 10000000, 100, 0, 1).Candidates()
    lastIndel = None

    # Loop through indels, normalising each one, and merging where appropriate.
    for rawIndel in rawIndelCandidates:
        normalisedIndel = chaplotype.convertVariantToStandardFormat(rawIndel, refFile, 100)

        if lastIndel is None:
            lastIndel = normalisedIndel
        elif normalisedIndel == lastIndel:
            lastIndel.nSupportingReads += normalisedIndel.nSupportingReads
        else:
            indelPos = lastIndel.refPos

            try:
                indels[indelPos].append(lastIndel)
            except KeyError:
                indels[indelPos] = [lastIndel]

            lastIndel = normalisedIndel

    # Make sure to catch last one
    if lastIndel is not None:
        indelPos = lastIndel.refPos
        try:
            indels[indelPos].append(lastIndel)
        except KeyError:
            indels[indelPos] = [lastIndel]

    return indels

###################################################################################################

def processBAM(file, ref, chrom, start, end, outdir, bqthr, mqthr, mode,rmdup, paired,getindels):
    """
    Process a single BAM file to produce coverage file.
    """
    #TODO: support for creation of .bes files and plots. Maybe better to create these as part of the error model module
    filebamname = os.path.basename(file)
    fileoutput = os.path.join(outdir, str(filebamname) + '.' + 'pileup.' + str(bqthr) +'thresholded.coverage')
    outf = open(fileoutput,'w')
    outf.write('chr offset loc ref_base ref_idx afwd ' + 'cfwd gfwd tfwd dfwd ifwd arev crev grev trev drev irev dseq iseq\n')
    bam = samtoolsWrapper.Samfile(file, "rb")
    pysamBam = pysam.Samfile(file, "rb")
    # reference sequence
    pysamRefFile = pysam.Fastafile(ref)

    for thisChrom,thisStart,thisEnd in zip(chrom,start,end):

        # get coverage for each region

        indelHash = getIndelsByPosition(thisChrom,thisStart,thisEnd,bam,mqthr,bqthr,ref)

        cover = Pileup(thisStart, thisEnd-1, mqthr, bqthr, mode,paired,rmdup)
        pysamBam.fetch(thisChrom, thisStart, thisEnd, callback=cover)
        seq = pysamRefFile.fetch(thisChrom, thisStart, thisEnd)

        for j in xrange(len(cover)):
            if thisStart+j in indelHash:
                outf.write(formatOutput(thisChrom, thisStart+j+1, seq[j], cover[j], indelHash[thisStart+j]))
            elif not getindels:
                outf.write(formatOutput(thisChrom, thisStart+j+1, seq[j], cover[j], None))

    outf.close()
    bam.close()
    pysamBam.close()
    pysamRefFile.close()

###################################################################################################
            

def readBAMs(options):
    """
    Read all BAM files and create coverage files.
    Several parallel jobs will be generated if ncpu > 1.
    """
    ## list of bam files
    bamfiles = parsePIF(options.pif)

    ## get target regions
    chrom,start,end = parseTGF(options.tgf)

    logger.info("Using " + str(options.ncpu) + " cores")
    ncpus = int(float(options.ncpu))
    start_time = time.time()
    
    if ncpus > 1:
        pool = Pool(processes=int(options.ncpu))
         
        jobs = []
        for file in bamfiles:
            jobitem = pool.apply_async(processBAM,(file,options.ref,chrom,start,end,options.outputdir,options.bqthr,options.mqthr,options.mode,options.rmdup,options.is_paired,options.getindels))
            jobs.append(jobitem)
        for result in jobs:
            result.get()

    else:
        for file in bamfiles:
            processBAM(file,options.ref,chrom, start, end,options.outputdir,options.bqthr,options.mqthr, options.mode,options.rmdup,options.is_paired,options.getindels)
    logger.info("Time elapsed: %.3f seconds", time.time() - start_time)

###################################################################################################

def main(options):
    """
    Top-level function, called from sysygy wrapper.
    """    
    readBAMs(options)

###################################################################################################

if __name__ == "__main__":
    raise Exeption("This module should not be called directly. Use 'syzygy [options] --module ReadBam' instead.")
