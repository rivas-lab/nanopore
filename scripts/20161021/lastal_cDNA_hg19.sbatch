#!/bin/bash
#SBATCH --job-name=lastal_cDNA_hg19
#SBATCH   --output=lastal_cDNA_hg19.out
#SBATCH    --error=lastal_cDNA_hg19.err
#################
#time you think you need; default is 2 hours
#format could be dd-hh:mm:ss, hh:mm:ss, mm:ss, or mm
#SBATCH --time=2:00:00
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

lastal -Q1 -D100 \
       ${PI_HOME}/data/hg19/hg19.lastdb \
       ${HOME}/data/cDNA-sample.fq \
       | \
       last-split \
	   -c0 \
	   > ${HOME}/data/cDNA-sample.maf
