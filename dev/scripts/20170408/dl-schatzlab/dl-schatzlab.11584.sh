#!/bin/bash

# download 
#
#  NA12878 Human Reference on Oxford Nanopore MinION
#    Read Alignments and Structural Variant Calls
#

#SBATCH --job-name=dl-schatzlab
#SBATCH   --output=dl-schatzlab.%j.out
#SBATCH    --error=dl-schatzlab.%j.err
#SBATCH --time=7-0:00:00
#SBATCH --qos=normal
#SBATCH -p mrivas
#SBATCH --nodes=1
#SBATCH --cores=1
#SBATCH --mem=2000
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ytanigaw@stanford.edu
#################
cat $0 > dl-schatzlab.${SLURM_JOBID}.sh
#ml load anaconda
#source activate pgenlib
#if [ ! -d ${LOCAL_SCRATCH} ]; then mkdir -p ${LOCAL_SCRATCH}; fi
#SCRATCH_TODAY=$SCRATCH/$(date +%Y%m%d)
#if [ ! -e ${SCRATCH_TODAY} ]; then mkdir -p ${SCRATCH_TODAY}; fi
#################

target=/scratch/PI/mrivas/data/nanopore-wgs-consortium-schatzlab
src=http://labshare.cshl.edu/shares/schatzlab/www-data/fsedlaze/NA12878_nanopore/version1.0/ngm_Nanopore_human_ngmlr-0.2.3_mapped.bam

cd $target

wget ${src}.sniffles1kb_auto_noalts.vcf_summary
wget ${src}.sniffles1kb_auto_noalts.vcf_summary.md5sum

wget ${src}.sniffles1kb_auto_noalts.vcf.gz
wget ${src}.sniffles1kb_auto_noalts.vcf.gz.md5sum

wget ${src}
wget ${src}.md5sum

