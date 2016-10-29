#!/bin/bash
#SBATCH --job-name=last_split.helio.10k
#SBATCH   --output=last_split.helio.10k.out
#SBATCH    --error=last_split.helio.10k.err
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

data_head="${PI_HOME}/data/nanopore/20161008_wgs_caucasian_48hr"
	  
cat ${data_head}.10k.helio.no-split.maf |
    last-split -c0 > ${data_head}.10k.helio.split.maf
