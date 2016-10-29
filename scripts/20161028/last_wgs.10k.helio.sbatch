#!/bin/bash
#SBATCH --job-name=last_wgs.10k.helio
#SBATCH   --output=last_wgs.10k.helio.out
#SBATCH    --error=last_wgs.10k.helio.err
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
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=15G
#################

data_head="${PI_HOME}/data/nanopore/20161008_wgs_caucasian_48hr"

parallel-fastq-j8 "lastal -Q1 -D100 -s2 -T0 -a1 ${PI_HOME}/data/hg19/hg19.lastdb | last-split -c0" \
		  < ${data_head}.10k.fq > ${data_head}.10k.helio.maf
