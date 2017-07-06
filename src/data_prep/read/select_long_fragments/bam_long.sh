#!/bin/bash
set -beEu -o pipefail

parallel --jobs=12 bash $(readlink -m ${0%.sh})-sub.sh :::: <( find $(dirname $(readlink -m $0))/bam/ -name '*.bam' )

