#!/bin/bash
#SBATCH --job-name=bowtie
#SBATCH   --output=bowtie.%j.out
#SBATCH    --error=bowtie.%j.err
#SBATCH --time=1:00:00
#SBATCH --qos=normal
#SBATCH -p normal
#SBATCH --nodes=1
#SBATCH --mem=16000
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ytanigaw@stanford.edu
#################
cat $0 > bowtie.$SLURM_JOBID.sh
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
fastq_cutoff_length=30
bam_mapq_min=50

# input data
fastqs=$PI_HOME/data/nanopore-wgs-consortium/*.fq.gz

# reference file
hg19=${PI_HOME}/data/hg19/hg19.fa

# intermediate files
fastq_filtered=${scratch_root}/ont.${fastq_cutoff_length}k.fq.gz
bam_unfiltered=${scratch_root}/ont.${fastq_cutoff_length}k.unsorted.unfiltered.bam
bam_unsorted=${scratch_root}/ont.${fastq_cutoff_length}k.unsorted.bam
bam_sorted=${scratch_root}/ont.${fastq_cutoff_length}k.bam

# filter fastq file by its length
if [ ! -f ${fastq_filtered} ]; then
    if [ "${verbose}" -eq 1 ]; then 
	echo "${fastqs} --> ${fastq_filtered}" >&2
    fi
    zcat $fastqs | head -n100 \
	| parallel --no-notice -k --jobs=${threads} --pipe -N4 \
		   "${project_root}/src/fastq_filter.sh -c ${fastq_cutoff_length} -k" \
	| gzip -9 - > ${fastq_filtered}	
elif [ ${verbose} -eq 1 ]; then
    echo "${fastq_filtered} already exists" >&2
fi

# map to reference genome
# filter by map quality
if [ ! -f ${bam_unfiltered} ]; then
    if [ "${verbose}" -eq 1 ]; then 
	echo "${fastq_filtered} --> ${bam_unfiltered}" >&2
    fi
    zcat ${fastq_filtered} \
	| bwa mem -x ont2d -t ${threads} $hg19 - \
	| samtools view -Sb -q ${bam_mapq_min} - > ${bam_unfiltered}
elif [ ${verbose} -eq 1 ]; then
    echo "${bam_unfiltered} already exists" >&2
fi

# filter by mapped fragment length and 
