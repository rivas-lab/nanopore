#!/bin/bash
#SBATCH --job-name=fastqc_wgs
#SBATCH   --output=fastqc_wgs.out
#SBATCH    --error=fastqc_wgs.err
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

mkdir -p "${data_head}.fastqc"
java -Xms8g -Xmx16g -classpath /share/PI/mrivas/bin/FastQC/ uk.ac.babraham.FastQC.FastQCApplication -o "${data_head}.fastqc" ${data_head}.fq
