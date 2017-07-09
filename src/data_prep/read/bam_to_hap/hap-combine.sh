#!/bin/bash
set -beEu -o pipefail

export MODULEPATH=$PI_HOME/.modules:$MODULEPATH
ml load ytanigaw-sh2

ary_lookup="$nanopore/src/data_prep/read/bam_to_hap/array-lookup.tsv"
out="$nanopore/public_data/intermediate/read/chr20.sorted"

echo "[$0] combining hap files to ${out}.hap"
parallel --jobs 10 cat {} :::: <( cat $ary_lookup | awk '{print $2}' | sed -e 's/bam$/hap/g' ) > ${out}.hap
echo "[$0] combining err files to ${out}.err"
parallel --jobs 10 cat {} :::: <( cat $ary_lookup | awk '{print $2}' | sed -e 's/bam$/err/g' ) > ${out}.err

