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
ml load anaconda
source activate pgenlib
#if [ ! -d ${LOCAL_SCRATCH} ]; then mkdir -p ${LOCAL_SCRATCH}; fi
#SCRATCH_TODAY=$SCRATCH/$(date +%Y%m%d)
#if [ ! -e ${SCRATCH_TODAY} ]; then mkdir -p ${SCRATCH_TODAY}; fi
#################

fasta=$PI_HOME/data/hg19/chr
bim=$PI_HOME/data/1000genomes/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes-pgen.bim
bam=$PI_HOME/data/nanopore-wgs-consortium/rel3/hg19/chr20/rel3.chr20.12500.10k.head.bam

hap=${bam%.bam}.hap
err=${bam%.bam}.err

cd $HOME/projects/nanopore/src/bamToHap/

python2 bamToHap.py --fasta $fasta --bim $bim $bam >$hap 2>$err


