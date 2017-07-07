#!/bin/bash
set -beEu -o pipefail

ary_lookup="/oak/stanford/groups/mrivas/users/ytanigaw/repos/rivas-lab/nanopore/src/data_prep/read/select_long_fragments/splithg19.tsv"

parallel --jobs 10 bash select_long_fragments-sub-parallel.sh {} :::: <( cat $ary_lookup | awk '{print $2}' )

