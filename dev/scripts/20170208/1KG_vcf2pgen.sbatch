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
#TEMP_DIR=$SCRATCH_TODAY

MEMORY=16000
THREADS=8
CHR=20

VCF=$DATA_DIR/ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
#TEMP=$TEMP_DIR/ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes-temp
PGEN=$DATA_DIR/ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes-pgen

echo $VCF >&2
#echo $TEMP >&2
echo $PGEN >&2

#${PI_HOME}/bin/plink --version >&2

#${PI_HOME}/bin/plink \
#    --memory ${MEMORY} \
#    --threads ${THREADS} \
#    --vcf ${VCF} \
#    --make-bed \
#    --out ${TEMP}

${PI_HOME}/bin/plink2 --version >&2

${PI_HOME}/bin/plink2 \
    --memory ${MEMORY} \
    --threads ${THREADS} \
    --vcf ${VCF} \
    --make-pgen \
    --out ${PGEN}

#${PI_HOME}/bin/plink2 \
#    --memory ${MEMORY} \
#    --threads ${THREADS} \
#    --bfile ${TEMP} \
#    --make-bpgen \
#    --out ${PGEN}
