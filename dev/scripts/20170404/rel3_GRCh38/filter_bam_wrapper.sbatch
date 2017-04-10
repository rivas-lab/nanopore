#!/bin/bash

# filter BAM file by mapped fragment length and FLAG 

#SBATCH --job-name=filter_bam_wrapper
#SBATCH   --output=filter_bam_wrapper.%j.out
#SBATCH    --error=filter_bam_wrapper.%j.err
#SBATCH --time=2-0:00:00
#SBATCH --qos=normal
#SBATCH -p mrivas
#SBATCH --nodes=1
#SBATCH --mem=40000
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ytanigaw@stanford.edu
#################
cat $0 > filter_bam_wrapper.${SLURM_JOBID}.sh
#ml load anaconda
#source activate pgenlib
if [ ! -d ${LOCAL_SCRATCH} ]; then mkdir -p ${LOCAL_SCRATCH}; fi
SCRATCH_TODAY=$SCRATCH/$(date +%Y%m%d)
if [ ! -e ${SCRATCH_TODAY} ]; then mkdir -p ${SCRATCH_TODAY}; fi

# parameters
verbose=1
threads=20
memory=40000
bam_cutoff_length=10
bam_mapq_min=50

# project config
project=nanopore
project_root=${HOME}/projects/${project}
scratch_root=${PI_SCRATCH}/projects/${project}



bam_original_path=/scratch/PI/mrivas/data/nanopore-wgs-consortium-rel3

for CHR in $(seq 1 22) X Y M ; do
    bam_original_name=chr${CHR}.sorted.bam
    bam_original=${bam_original_path}/${bam_original_name}

    tmp_sam_file=${SCRATCH_TODAY}/${SLURM_JOBID}.sam
    tmp_bam_file=${tmp_sam_file%.sam}.bam

    bam_filtered_name=${bam_original_name%.bam}.${bam_cutoff_length}k.bam
    bam_filtered=${bam_original_path}/${bam_filtered_name}

    if [ ! -f ${bam_filtered} ]; then
	if [ "${verbose}" -eq 1 ]; then
	    echo "${bam_original} --> ${bam_filtered}" >&2
	fi

	# filter BAM file and write to a sam file

	/share/PI/mrivas/bin/samtools view -H ${bam_original} > ${tmp_sam_file}
	/share/PI/mrivas/bin/samtools view -q ${bam_mapq_min} ${bam_original} \
	    | /share/PI/mrivas/bin/parallel --no-notice -k --jobs=${threads} --pipe  \
					    "${project_root}/src/bam_filter.sh -c ${bam_cutoff_length} -k" \
					    >> ${tmp_sam_file}


	# convert to BAM
	/share/PI/mrivas/bin/samtools view -Sb ${tmp_sam_file} > ${tmp_bam_file}


	# sort

	if [ "${verbose}" -eq 1 ]; then
	    echo "sorting ${tmp_bam_file}" >&2
	fi

	/share/PI/mrivas/bin/samtools sort \
				      -l 9 \
				      -@ ${threads} \
				      -m ${memory}M \
				      -o ${bam_filtered} \
				      ${tmp_bam_file}

	/share/PI/mrivas/bin/samtools index ${bam_filtered}

    elif [ ${verbose} -eq 1 ]; then
	echo "${bam_filtered} already exists" >&2
    fi
done
