#!/bin/bash
set -beEu -o pipefail

inBam="../../../../public_data/input/read/chr20.sorted.bam"
outSam="../../../../sandbox/data_prep/read/remap_from_hg38_to_hg19/chr20.sorted.split"

bash $(dirname $(readlink -m $0))/$(basename ${0%-exec.sh}).sh $inBam $outSam
