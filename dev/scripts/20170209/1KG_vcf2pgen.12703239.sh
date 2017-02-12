#!/bin/bash
#SBATCH --job-name=1KG_vcf2pgen
#SBATCH   --output=1KG_vcf2pgen.%j.out
#SBATCH    --error=1KG_vcf2pgen.%j.err
#SBATCH --time=12:00:00
#SBATCH --qos=normal
#SBATCH -p normal
#SBATCH --nodes=1
#SBATCH --mem=16000
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ytanigaw@stanford.edu
#################
cat $0 > 1KG_vcf2pgen.${SLURM_JOBID}.sh
ml load anaconda
source activate pgenlib
if [ ! -d ${LOCAL_SCRATCH} ]; then mkdir -p ${LOCAL_SCRATCH}; fi
SCRATCH_TODAY=$SCRATCH/$(date +%Y%m%d)
if [ ! -e ${SCRATCH_TODAY} ]; then mkdir -p ${SCRATCH_TODAY}; fi
#################

DATA_DIR=$PI_HOME/data/1000genomes

MEMORY=16000
THREADS=8

for CHR in `seq 1 19` 21 22; do
    VCF=$DATA_DIR/ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
    PGEN=$DATA_DIR/ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes-pgen

    echo $VCF >&2
    echo $PGEN >&2

    ${PI_HOME}/bin/plink2-20170208 --version >&2

    ${PI_HOME}/bin/plink2-20170208 \
	      --memory ${MEMORY} \
	      --threads ${THREADS} \
	      --vcf ${VCF} \
	      --make-bpgen \
	      --out ${PGEN}
done