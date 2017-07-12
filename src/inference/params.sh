#!/bin/sh
set -beEu -o pipefail

homo_thr=0.9

# input
## read
hap_f="../../public_data/intermediate/read/chr20.sorted.hap"

## population reference
pgen_f="../../public_data/intermediate/population_ref/chr20-alleles.pgen"
bins_f="../../public_data/intermediate/population_ref/chr20-bins.tsv"

## intermediate files
inference_d="../../public_data/intermediate/inference"
hapkey_d="$inference_d/hapkey/chr20"
hapcnt_d="$inference_d/hapcnt/chr20"
log_likelihood_d="$inference_d/log_likelihood/chr20/chr20.sorted"
log_posterior_d="$inference_d/log_posterior/chr20/chr20.sorted"

# results file
pgen_out="../../private_data/output/$(basename ${hap_f%.hap})-$(basename ${pgen_f%.pgen}).pgen"
