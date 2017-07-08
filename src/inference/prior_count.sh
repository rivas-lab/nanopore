#!/bin/bash
set -beEu -o pipefail

. $(dirname $(readlink -m $0))/params.sh

python $(basename ${0%.sh}).py --pgen $pgen_f --bins $bins_f --cnt $hapcnt_d --key $hapkey_d  2>&1 | tee $(basename ${0%.sh}).out

