#!/bin/bash
#SBATCH --job-name=nanopore-wgs.dl
#SBATCH   --output=nanopore-wgs.dl.out
#SBATCH    --error=nanopore-wgs.dl.err
#SBATCH --time=0-10:00:00
#SBATCH --qos=normal
#SBATCH -p normal
#SBATCH --nodes=1
#SBATCH --mem=4000
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ytanigaw@stanford.edu
#################

DIR=${PI_HOME}/data/nanopore-wgs-consortium
cd ${DIR}
wget http://s3.climb.ac.uk/nanopore-human-wgs/152889212-FAB45271.fastq.gz
