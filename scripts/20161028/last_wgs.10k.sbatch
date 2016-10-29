#!/bin/bash
#SBATCH --job-name=last_wgs.10k
#SBATCH   --output=last_wgs.10k.out
#SBATCH    --error=last_wgs.10k.err
#################
#time you think you need; default is 2 hours
#format could be dd-hh:mm:ss, hh:mm:ss, mm:ss, or mm
#SBATCH --time=12:00:00
#################
#quality of service; think of it as job priority
#SBATCH --qos=normal
# We are submitting to the dev partition, there are several on sherlock: normal, gpu, owners, hns, bigmem (jobs requiring >64Gigs RAM)
#
#SBATCH -p normal
#################
#number of nodes you are requesting
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8000
#################

data_head="${PI_HOME}/data/nanopore/20161008_wgs_caucasian_48hr"
	  
lastal -P 8 -Q1 -D100 \
       ${PI_HOME}/data/hg19/hg19.lastdb \
       ${data_head}.10k.fq | last-split -c0  > ${data_head}.10k.maf
