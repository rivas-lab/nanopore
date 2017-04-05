#!/bin/bash

# mapping (perform bwa mem and samtools sort)

#SBATCH --job-name=bwa
#SBATCH   --output=bwa.%j.out
#SBATCH    --error=bwa.%j.err
#SBATCH --time=2-0:00:00
#SBATCH --qos=normal
#SBATCH -p mrivas
#SBATCH --nodes=1
#SBATCH --mem=40000
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ytanigaw@stanford.edu
#################
cat $0 > bwa.${SLURM_JOBID}.sh
#ml load anaconda
#source activate pgenlib
if [ ! -d ${LOCAL_SCRATCH} ]; then mkdir -p ${LOCAL_SCRATCH}; fi
SCRATCH_TODAY=$SCRATCH/$(date +%Y%m%d)
if [ ! -e ${SCRATCH_TODAY} ]; then mkdir -p ${SCRATCH_TODAY}; fi
#################

# project config
project=nanopore
project_root=${HOME}/projects/${project}
scratch_root=${PI_SCRATCH}/projects/${project}

# parameters
verbose=1
threads=20
memory=40000

# reference file
ref=${PI_HOME}/data/hg19/hg19.fa

# input data

data_prefix=${PI_HOME}/data/nanopore-wgs-consortium/rel1/rel1.25k
fastq=${data_prefix}.fq.gz

# intermediate file
bam_all_tmp=${data_prefix}.tmp.bam

# sorted mapped bam file
bam_all=${data_prefix}.bam

cd $LOCAL_SCRATCH
pwd >&2

# map to reference genome
if [ ! -f ${bam_all} ]; then
    if [ "${verbose}" -eq 1 ]; then
	echo "${fastq} --> ${bam_all}" >&2
    fi

    if [ ! -f ${ref}.sa ]; then
	if [ "${verbose}" -eq 1 ]; then
	    echo "indexing ${ref} with bwa index" >&2
	fi
	/share/PI/mrivas/bin/bwa index ${ref}
    fi

    if [ ! -f ${bam_all_tmp} ]; then
	zcat ${fastq} \
	    | /share/PI/mrivas/bin/bwa mem -x ont2d -t ${threads} $ref - \
	    | /share/PI/mrivas/bin/samtools view -Sb - > ${bam_all_tmp}
    fi

    if [ ! -f ${bam_all} ]; then
	if [ "${verbose}" -eq 1 ]; then
	    echo "sorting ${bam_all_tmp}" >&2
	fi

	/share/PI/mrivas/bin/samtools sort \
				      -l 9 \
				      -@ ${threads} \
				      -m ${memory}M \
				      -o ${bam_all} \
				      ${bam_all_tmp}
    fi

    if [ "${verbose}" -eq 1 ]; then
	echo "indexing ${bam_all}" >&2
    fi

    /share/PI/mrivas/bin/samtools index ${bam_all}

    echo "please delete ${bam_all_tmp}"

elif [ ${verbose} -eq 1 ]; then
    echo "${bam_all} already exists" >&2
fi
