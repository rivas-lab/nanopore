#!/bin/bash

#SBATCH --job-name=WTCHG-compress
#SBATCH   --output=WTCHG-compress.%j.out
#SBATCH    --error=WTCHG-compress.%j.err
#SBATCH --time=2-0:00:00
#SBATCH --qos=normal
#SBATCH -p mrivas
#SBATCH --nodes=1
#SBATCH --mem=24000
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ytanigaw@stanford.edu
#################
cat $0 > WTCHG-compress.${SLURM_JOBID}.sh

TMP=$(mktemp -d)
cp /share/PI/mrivas/data/nanopore-WTCHG/WTCHG-Gplc-NA12878-minion-R9.4-v0.3.tar $TMP/
cd $TMP
gzip -9 WTCHG-Gplc-NA12878-minion-R9.4-v0.3.tar
rm /share/PI/mrivas/data/nanopore-WTCHG/WTCHG-Gplc-NA12878-minion-R9.4-v0.3.tar
cp $TMP/WTCHG-Gplc-NA12878-minion-R9.4-v0.3.tar.gz /share/PI/mrivas/data/nanopore-WTCHG/
rm -rf $TMP
#################
