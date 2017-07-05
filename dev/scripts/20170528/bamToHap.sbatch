#!/bin/bash

# test bamToHap.py script on cluster

#SBATCH --job-name=bamToHap
#SBATCH   --output=bamToHap.%j.out
#SBATCH    --error=bamToHap.%j.err
#SBATCH --time=6:00:00
#SBATCH --qos=normal
#SBATCH -p normal
#SBATCH --nodes=1
#SBATCH --cores=1
#SBATCH --mem=4000
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ytanigaw@stanford.edu
#################
cat $0 > bamToHap.${SLURM_JOBID}.sh
#ml load anaconda
source activate pgenlib
#if [ ! -d ${LOCAL_SCRATCH} ]; then mkdir -p ${LOCAL_SCRATCH}; fi
#SCRATCH_TODAY=$SCRATCH/$(date +%Y%m%d)
#if [ ! -e ${SCRATCH_TODAY} ]; then mkdir -p ${SCRATCH_TODAY}; fi
#################

fasta=/oak/stanford/groups/mrivas/public_data/hg19/chr
bim=/oak/stanford/groups/mrivas/private_data/ukbb/pgen/chr20impv1-keep-maf001-snv-biallelic-pgen.bim
bam=/oak/stanford/groups/mrivas/public_data/nanopore-wgs-consortium/rel3/hg19/chr20/rel3.chr20.12500.10k.bam

hap=${bam%.bam}-ukbb-keep-maf001-snv-biallelic-pgen.hap
err=${hap%.hap}.err

src=$HOME/projects/nanopore/src/bamToHap/

python2 $src/bamToHap.py --fasta $fasta --bim $bim $bam >$hap 2>$err


