#!/bin/bash
#SBATCH --job-name=bwa.10k
#SBATCH   --output=bwa.10k.out
#SBATCH    --error=bwa.10k.err
#################
#time you think you need; default is 2 hours
#format could be dd-hh:mm:ss, hh:mm:ss, mm:ss, or mm
#SBATCH --time=2-0:00:00
#################
#quality of service; think of it as job priority
#SBATCH --qos=normal
# We are submitting to the dev partition, there are several on sherlock: normal, gpu, owners, hns, bigmem (jobs requiring >64Gigs RAM)
#
#SBATCH -p normal
#################
#number of nodes you are requesting
#SBATCH --nodes=1
#SBATCH --mem=63000
#################

ref="${PI_HOME}/data/hg19/hg19.fa"
data_head="${PI_HOME}/data/nanopore/20161008_wgs_caucasian_48hr"

bwa mem -x ont2d ${ref} ${data_head}.10k.fq > ${data_head}.10k.bwa.sam
