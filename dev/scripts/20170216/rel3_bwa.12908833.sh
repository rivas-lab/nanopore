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

# intermediate file
bam_long_tmp=${data_prefix}.geq${bam_cutoff_length}k.tmp.bam

# sorted mapped bam file with long mapped fragments
bam_long=${data_prefix}.geq${bam_cutoff_length}k.bam

# map to reference genome
if [ ! -f ${bam_all} ]; then
    if [ "${verbose}" -eq 1 ]; then 
	echo "${fastq_filtered} --> ${bam_all}" >&2
    fi

    if [ ! -f ${ref}.sa ]; then
	if [ "${verbose}" -eq 1 ]; then 
	    echo "indexing ${ref} with bwa index" >&2
	fi
	bwa index ${ref}
    fi
    
    zcat ${fastq} \
	| bwa mem -x ont2d -t ${threads} $ref - \
	| samtools view -Sb - > ${bam_all_tmp}
    samtools sort \
	     -l 9 \
	     -@ ${threads} -m ${memory} \
	     ${bam_all_tmp} > ${bam_all}
    samtools index ${bam_all}
    rm ${bam_all_tmp}
elif [ ${verbose} -eq 1 ]; then
    echo "${bam_all} already exists" >&2
fi

# filter by map quality score, mapped fragment length and
# bitwise FLAG (drop supplementary alignmnet)
if [ ! -f ${bam_long} ]; then
    if [ "${verbose}" -eq 1 ]; then 
	echo "${bam_all} --> ${bam_long}" >&2
    fi
    cat <(samtools view -H ${bam_all} > ${sam_long}) \
	<(samtools view -q ${bam_mapq_min} ${bam_all} \	   
    		| parallel --no-notice -k --jobs=${threads} --pipe  \
			   "${project_root}/src/bam_filter.sh -c ${bam_cutoff_length} -k") \
	| samtools view -Sb - > ${bam_long_tmp}		      

    samtools sort \
	     -l 9 \
	     -@ ${threads} -m ${memory} \
	     ${bam_long_tmp} > ${bam_long}
    samtools index ${bam_long}
    rm ${bam_long_tmp}
   
elif [ ${verbose} -eq 1 ]; then
    echo "${bam_long} already exists" >&2
fi
