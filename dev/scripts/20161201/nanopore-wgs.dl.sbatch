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

LOC="http://s3.climb.ac.uk/nanopore-human-wgs"

DIR=${PI_HOME}/data/nanopore-wgs-consortium
CUTOFF=25000

PROG="/home/ytanigaw/projects/nanopore/src/filter_by_length_fastq.sh"

FILE="3775529215-FAB49914"

while read FILE; do
    echo $FILE
    if [ -e ${DIR}/${FILE}.${CUTOFF}.fq.gz ]; then
	echo "[skip] $FILE"
    else
	echo $FILE
	wget ${LOC}/${FILE}.fastq.gz -O - \
	    | zcat \
	    | ${PROG} ${CUTOFF} \
	    | gzip -9 > ${DIR}/${FILE}.${CUTOFF}.fq.gz
    fi
done < <(wget -O - https://raw.githubusercontent.com/nanopore-wgs-consortium/NA12878/master/README.md |grep "http:" |awk '{print $13}'|sed 's/.*wgs\/\(.*\).fastq.gz.*/\1/g')
