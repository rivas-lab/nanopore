#!/bin/bash
#
#all commands that start with SBATCH contain commands that are just used by SLURM for scheduling
#################
#set a job name
#SBATCH --job-name=count-fast5
#################
#a file for job output, you can check job progress
#SBATCH --output=20161008_wgs_caucasian_48hr-count-fast5.out
#################
# a file for errors from the job
#SBATCH --error=20161008_wgs_caucasian_48hr-count-fast5.err
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
#SBATCH --mem=4000
#################

find ${HOME}/data/nanopore/20161008_wgs_caucasian_48hr/ -name "*.fast5"| wc -l
