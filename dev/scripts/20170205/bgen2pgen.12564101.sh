#!/bin/bash
#SBATCH --job-name=bgen2pgen
#SBATCH   --output=bgen2pgen.%j.out
#SBATCH    --error=bgen2pgen.%j.err
#SBATCH --time=1:00:00
#SBATCH --qos=normal
#SBATCH -p dev
#SBATCH --nodes=1
#SBATCH --mem=8000
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ytanigaw@stanford.edu
#################
cat $0 > bgen2pgen.${SLURM_JOBID}.sh
ml load anaconda
source activate pgenlib
if [ ! -d ${LOCAL_SCRATCH} ]; then mkdir -p ${LOCAL_SCRATCH}; fi
SCRATCH_TODAY=$SCRATCH/$(date +%Y%m%d)
if [ ! -e ${SCRATCH_TODAY} ]; then mkdir -p ${SCRATCH_TODAY}; fi
#################

DATA_DIR=$PI_HOME/ukbb/download
TEMP_DIR=$SCRATCH_TODAY

BGEN=$DATA_DIR/chr20impv1.bgen
PGEN=$DATA_DIR/chr20impv1-pgen
SAMPLE=$DATA_DIR/impv1.sample
TEMP=$TEMP_DIR/chr20impv1-temp

echo $BGEN >&2
echo $SAMPLE >&2
echo $TEMP >&2
echo $PGEN >&2

${PI_HOME}/bin/plink --version >&2

${PI_HOME}/bin/plink \
    --memory 8000 \
    --bgen ${BGEN} \
    --sample ${SAMPLE} \
    --out ${TEMP}

${PI_HOME}/bin/plink2 --version >&2

${PI_HOME}/bin/plink2 \
    --memory 8000 \
    --bfile ${TEMP} \
    --make-bpgen \
    --out ${PGEN}
