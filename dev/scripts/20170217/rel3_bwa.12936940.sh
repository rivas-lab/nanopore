#!/bin/bash
#SBATCH --job-name=rel3_bwa
#SBATCH   --output=rel3_bwa.%j.out
#SBATCH    --error=rel3_bwa.%j.err
#SBATCH --time=2-0:00:00
#SBATCH --qos=normal
#SBATCH -p normal
#SBATCH --nodes=1
#SBATCH --mem=24000
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ytanigaw@stanford.edu
#################
cat $0 > rel3_bwa.${SLURM_JOBID}.sh
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
threads=12
memory=24000
bam_cutoff_length=10
bam_mapq_min=50
chr=20

# reference file
ref=${PI_HOME}/data/hg19/chr${chr}.fa

# input data
data_prefix=${PI_SCRATCH}/data/NA12878-fast5/poretools_fastq.12894489.geq12500
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
    
    zcat ${fastq} \
	| /share/PI/mrivas/bin/bwa mem -x ont2d -t ${threads} $ref - \				   
	| samtools view -Sb - > ${bam_all_tmp}

    if [ "${verbose}" -eq 1 ]; then 
	echo "sorting ${bam_all_tmp}" >&2
    fi

    /share/PI/mrivas/bin/samtools sort \
	     -f -l 9 \
	     -@ ${threads} \
	     -m ${memory}M \
	     ${bam_all_tmp} \
	     ${bam_all}

    if [ "${verbose}" -eq 1 ]; then 
	echo "indexing ${bam_all}" >&2
    fi

    /share/PI/mrivas/bin/samtools index ${bam_all} 
#    rm ${bam_all_tmp}
elif [ ${verbose} -eq 1 ]; then
    echo "${bam_all} already exists" >&2
fi
