#!/bin/bash
#SBATCH --job-name=count.plane.no-split.10k
#SBATCH   --output=count.plane.no-split.10k.out
#SBATCH    --error=count.plane.no-split.10k.err
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
#SBATCH --mem=8000
#################

data_head="${PI_HOME}/data/nanopore/20161008_wgs_caucasian_48hr"

ml load anaconda/anaconda3

${HOME}/projects/nanopore/src/count.py \
       -i ${data_head}.10k.plane.no-split.maf \
       -o ${data_head}.10k.plane.no-split.count
