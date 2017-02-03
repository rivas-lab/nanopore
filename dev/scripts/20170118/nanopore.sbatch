#!/bin/bash
#SBATCH --job-name=nanopore
#SBATCH   --output=nanopore.%j.out
#SBATCH    --error=nanopore.%j.err
#SBATCH --time=12:00:00
#SBATCH --qos=normal
#SBATCH -p normal
#SBATCH --nodes=1
#SBATCH --mem=16000
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ytanigaw@stanford.edu
#################
cat $0 > nanopore.${SLURM_JOBID}.sh
ml load anaconda
source activate pgenlib
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
fastq_cutoff_length=0
bam_cutoff_length=10
bam_mapq_min=50

# input data
dataset_name=wtchg-chr22
fastqs=${PI_HOME}/data/nanopore-wgs-consortium/*.fq.gz
bam=${PI_HOME}/data/WTCHG/v0.3/NA12878-minion-R9.4-pass-v0.3_chr22.bam

# reference file
hg19=${PI_HOME}/data/hg19/hg19.fa

# intermediate files
fastq_filtered=${scratch_root}/${dataset_name}.${fastq_cutoff_length}k.fq.gz
bam_unfiltered=${scratch_root}/${dataset_name}.${fastq_cutoff_length}k.unsorted.unfiltered.bam
sam_unsorted=${scratch_root}/${dataset_name}.${fastq_cutoff_length}k.unsorted.sam
bam_unsorted=${scratch_root}/${dataset_name}.${fastq_cutoff_length}k.unsorted.bam
bam_sorted=${scratch_root}/${dataset_name}.${fastq_cutoff_length}k.bam

## filter fastq file by its length
#if [ ! -f ${fastq_filtered} ]; then
#    if [ "${verbose}" -eq 1 ]; then 
#	echo "${fastqs} --> ${fastq_filtered}" >&2
#    fi
#    zcat ${fastqs} | head -n100 \
#	| parallel --no-notice -k --jobs=${threads} --pipe -N4 \
#		   "${project_root}/src/fastq_filter.sh -c ${fastq_cutoff_length} -k" \
#	| gzip -9 - > ${fastq_filtered}	
#elif [ ${verbose} -eq 1 ]; then
#    echo "${fastq_filtered} already exists" >&2
#fi

# map to reference genome
# filter by map quality
#if [ ! -f ${bam_unfiltered} ]; then
#    if [ "${verbose}" -eq 1 ]; then 
#	echo "${fastq_filtered} --> ${bam_unfiltered}" >&2
#    fi
#    zcat ${fastq_filtered} \
#	| bwa mem -x ont2d -t ${threads} $hg19 - \
#	| samtools view -Sb - > ${bam_unfiltered}
#elif [ ${verbose} -eq 1 ]; then
#    echo "${bam_unfiltered} already exists" >&2
#fi

# use the provided bam file
cp ${bam} ${bam_unfiltered}

# filter by map quality score, mapped fragment length and
# bitwise FLAG (drop supplementary alignmnet)
if [ ! -f ${bam_unsorted} ]; then
    if [ "${verbose}" -eq 1 ]; then 
	echo "${bam_unfiltered} --> ${bam_unsorted}" >&2
    fi
    samtools view -H ${bam_unfiltered} > ${sam_unsorted}
    samtools view -q ${bam_mapq_min} ${bam_unfiltered} \
	| parallel --no-notice -k --jobs=${threads} --pipe  \
		   "${project_root}/src/bam_filter.sh -c ${bam_cutoff_length} -k" \
	>> ${sam_unsorted}
    samtools view -Sb ${sam_unsorted} > ${bam_unsorted}
    rm ${sam_unsorted}
elif [ ${verbose} -eq 1 ]; then
    echo "${bam_unsorted} already exists" >&2
fi

