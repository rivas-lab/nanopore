#!/bin/bash
set -beEu -o pipefail

fasta=/oak/stanford/groups/mrivas/public_data/hg19/chr
bim=../../../../private_data/intermediate/population_ref/chr20-alleles.bim 

ary_lookup="/oak/stanford/groups/mrivas/users/ytanigaw/repos/rivas-lab/nanopore/src/data_prep/read/bam_to_hap/array-lookup.tsv"

parallel --jobs 10 bash ${0%.sh}-sub.sh $fasta $bim {} :::: <( cat $ary_lookup | awk '{print $2}' )

