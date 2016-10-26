#!/bin/bash
#SBATCH --job-name=length_dist_wgs
#SBATCH   --output=length_dist_wgs.out
#SBATCH    --error=length_dist_wgs.err
#################
#time you think you need; default is 2 hours
#format could be dd-hh:mm:ss, hh:mm:ss, mm:ss, or mm
#SBATCH --time=0:30:00
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
#SBATCH --mem=4000
#################

data_head="${PI_HOME}/data/nanopore/20161008_wgs_caucasian_48hr"

awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' ${data_head}.fq |sort -nr > ${data_head}.fq.length
