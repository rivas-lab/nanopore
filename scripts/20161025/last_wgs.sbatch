#!/bin/bash
#SBATCH --job-name=last_wgs
#SBATCH   --output=last_wgs.out
#SBATCH    --error=last_wgs.err
#################
#time you think you need; default is 2 hours
#format could be dd-hh:mm:ss, hh:mm:ss, mm:ss, or mm
#SBATCH --time=4:00:00
#################
#quality of service; think of it as job priority
#SBATCH --qos=normal
# We are submitting to the dev partition, there are several on sherlock: normal, gpu, owners, hns, bigmem (jobs requiring >64Gigs RAM)
#
#SBATCH -p normal
#################
#number of nodes you are requesting
#SBATCH --nodes=1
#################
#memory per node; default is 4000 MB per CPU
#SBATCH --mem=16000
#################

data_head="${PI_HOME}/data/nanopore/20161008_wgs_caucasian_48hr"

ml load python/2.7.5
ml load hdf5/1.8.16

poretools fastq ${data_head} > ${data_head}.fq
	  
lastal -Q1 -D100 \
       ${PI_HOME}/data/hg19/hg19.lastdb \
       ${data_head}.fq | last-split -c0  > ${data_head}.maf
