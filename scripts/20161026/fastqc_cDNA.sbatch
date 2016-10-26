#!/bin/bash
#SBATCH --job-name=fastqc_cDNA
#SBATCH   --output=fastqc_cDNA.out
#SBATCH    --error=fastqc_cDNA.err
#################
#time you think you need; default is 2 hours
#format could be dd-hh:mm:ss, hh:mm:ss, mm:ss, or mm
#SBATCH --time=0:20:00
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
#SBATCH --cpus-per-task=16
#################

data_head="${PI_HOME}/data/nanopore/20161006_minion_human_cDNA"
mkdir -p "${data_head}.fastqc"
fastqc --nano -t 16 -o "${data_head}.fastqc" "${data_head}/"
