#!/bin/bash
#SBATCH --job-name=last_cDNA
#SBATCH   --output=last_cDNA.out
#SBATCH    --error=last_cDNA.err
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
#SBATCH --mem=64000
#################

data_head="${PI_HOME}/data/nanopore/20161006_minion_human_cDNA"

#ml load python/2.7.5
#ml load hdf5/1.8.16

#poretools fastq ${data_head} > ${data_head}.fq
	  
lastal -Q1 -D100 \
       ${PI_HOME}/data/hg19/hg19.lastdb \
       ${data_head}.fq | last-split -c0  > ${data_head}.maf
