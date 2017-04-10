import numpy as np
import pgenlib as pg


class PgenFile:
    """Class to load pgenfile
    """
    def __init__(self, filename):
        self.filename = filename
        self.pgen = {}

    def load(self, chromosome):
        self.pgen[self.uniq_key(chromosome)] = \
            pg.PgenReader(self.filename)

    def uniq_key(self, chromosome):
        return str(chromosome)

    def is_loaded(self, chromosome):
        return self.uniq_key(chromosome) in self.pgen

    def get_pgen(self, chromosome):
        if not self.is_loaded(chromosome):
            self.load(chromosome)
        return self.pgen[self.uniq_key(chromosome)]

    def read_alleles(self, chromosome, index_l, index_r):
        pgen = self.get_pgen(chromosome)
        sample_ct = pgen.get_raw_sample_ct()
        alleles_list = np.zeros((index_r - index_l, 2 * sample_ct),
                                dtype=np.int32)
        pgen.read_alleles_range(index_l, index_r, alleles_list)
        return alleles_list
