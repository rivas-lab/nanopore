# nanopore
Tools for analysis of long-read sequencing data.


# Directory Structure
- src: source files
- dev: scratch partition for development

# Method
## pipline
- FAST5
  - extract sequence information with portools, etc.
- FASTQ.gz
  - `src/fastq_filter.sh`
  - filter by read length
- FASTQ.gz (filtered)
  - map to the reference genome with bwa 
- BAM (mapped)
  - `src/bam_filter.sh`
  - filter by the following
    - mapping quality (MAPQ)
    - mapped fragment length
    - drop secondary alignment
- BAM (filtered)
  - count match, mismatch to estimate mismatch ratio
  - drop fragments with $\geq 10%$ mismatch ratio
  - query mismatches to SNP data base (dbSNP)
  - drop fragments with $\leq 10$ SNPs
  - [optional] query putative SNPs to Platinum genome
- BAM (informative) + SNP summary file
