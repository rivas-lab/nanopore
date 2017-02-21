#!/bin/sh

source activate pgenlib

which python
python --version

python2 <<EOF

import numpy as np
import pgenlib as pg

filename = './test.pgen'
sample_ct = 1
variant_ct = 3
nonref_flags = True

allele_codes = np.zeros(sample_ct * 2, dtype = np.int32)

pgwriter = pg.PgenWriter(filename = filename, 
                         sample_ct = sample_ct,
                         variant_ct = variant_ct, 
                         nonref_flags = nonref_flags, 
                         hardcall_phase_present = True)

for i in range(variant_ct):
    pgwriter.append_alleles(allele_codes = allele_codes,
                            all_phased = True)

#    if we pass as our args as positional args, it works perfectly.
#    pgwriter.append_alleles(allele_codes, True)

pgwriter.close()

EOF
