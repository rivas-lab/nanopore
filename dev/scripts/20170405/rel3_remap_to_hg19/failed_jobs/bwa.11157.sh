#!/bin/bash

# mapping (perform bwa mem and samtools sort)
#   Given a bam file mapped to GRCh38,
#    - extract fastq
#    - map to hg19
#    - sort & index

#SBATCH --job-name=bwa
#SBATCH   --output=bwa.%j.out
#SBATCH    --error=bwa.%j.err
#SBATCH --time=2-0:00:00
#SBATCH --qos=normal
#SBATCH -p mrivas
#SBATCH --nodes=1
#SBATCH --mem=20000
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
threads=10
memory=20000

cd $SCRATCH_TODAY
echo "current location (SCRATCH_TODAY) is:" >&2
pwd >&2

for CHR in $(seq 8 14) ; do

    # reference file
    ref=${PI_HOME}/data/hg19/chr${CHR}.fa

    # input data
    source_data_loc=${PI_SCRATCH}/data/nanopore-wgs-consortium-rel3
    fastq=${source_data_loc}/chr${CHR}.sorted.bam

    target_data_loc=${PI_SCRATCH}/data/nanopore-wgs-consortium-rel3-remap
    data_prefix=${target_data_loc}/chr${CHR}

    # intermediate file
    bam_all_tmp=${SCRATCH_TODAY}/chr${CHR}.tmp.bam

    # sorted mapped bam file
    bam_all=${data_prefix}.bam

    # map to reference genome
    if [ ! -f ${bam_all} ]; then
	if [ ${verbose} -eq 1 ]; then
	    echo "${fastq} --> ${bam_all}" >&2
	fi

	if [ ! -f ${ref}.sa ]; then
	    if [ ${verbose} -eq 1 ]; then
		echo "indexing ${ref} with bwa index" >&2
	    fi
	    /share/PI/mrivas/bin/bwa index ${ref}
	fi

	if [ ! -f ${bam_all_tmp} ]; then
	    /share/PI/mrivas/bin/samtools fastq ${fastq} \
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

	if [ ${verbose} -eq 1 ]; then
	    echo "indexing ${bam_all}" >&2
	fi

	/share/PI/mrivas/bin/samtools index ${bam_all}

	echo "please delete ${bam_all_tmp}" >&2
	echo "written to ${bam_all}"

    elif [ ${verbose} -eq 1 ]; then
	echo "${bam_all} already exists" >&2
    fi
done
