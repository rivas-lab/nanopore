#!/bin/bash
set -beEu -o pipefail

. $(dirname $(readlink -m $0))/params.sh

python $(basename ${0%.sh}).py --bins $bins_f --key $hapkey_d --lp $log_posterior_d --out $pgen_out --homo $homo_thr 2>&1 | tee $(basename ${0%.sh}).out

