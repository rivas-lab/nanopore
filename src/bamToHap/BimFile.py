import numpy as np
import pandas as pd
import bisect


class BimFile:
    """A class to manage collection of bim files
    """
    def __init__(self, filename):
        self.filename = filename
        self.raw_tables = {}
        self.id = {}
        self.morgan = {}
        self.bp = {}
        self.allele_1 = {}
        self.allele_2 = {}

    def load_raw_table(self, chromosome, index_base=1):
        self.raw_tables[self.uniq_key(chromosome)] = \
            pd.read_csv(self.filename, sep='\t',
                        names=['chr', 'id', 'morgan', 'bp', 'allele_1', 'allele_2'])

        self.id[self.uniq_key(chromosome)] = \
            np.array(self.raw_tables[self.uniq_key(chromosome)].ix[:, 1])

        self.morgan[self.uniq_key(chromosome)] = \
            np.array(self.raw_tables[self.uniq_key(chromosome)].ix[:, 2])

        self.bp[self.uniq_key(chromosome)] = \
            np.array([int(x) - index_base for x
                      in self.raw_tables[self.uniq_key(chromosome)].ix[:, 3]])

        self.allele_1[self.uniq_key(chromosome)] = \
            np.array(self.raw_tables[self.uniq_key(chromosome)].ix[:, 4])

        self.allele_2[self.uniq_key(chromosome)] = \
            np.array(self.raw_tables[self.uniq_key(chromosome)].ix[:, 5])

    def is_loaded(self, chromosome):
        return self.uniq_key(chromosome) in self.raw_tables

    def uniq_key(self, chromosome):
        return str(chromosome)

    def get_raw_table(self, chromosome):
        if not self.is_loaded(chromosome):
            self.load_raw_table(chromosome)
        return self.raw_tables[self.uniq_key(chromosome)]

    def get_id(self, chromosome):
        if not self.is_loaded(chromosome):
            self.load_raw_table(chromosome)
        return self.id[self.uniq_key(chromosome)]

    def get_morgan(self, chromosome):
        if not self.is_loaded(chromosome):
            self.load_raw_table(chromosome)
        return self.morgan[self.uniq_key(chromosome)]

    def get_bp(self, chromosome):
        if not self.is_loaded(chromosome):
            self.load_raw_table(chromosome)
        return self.bp[self.uniq_key(chromosome)]

    def get_allele_1(self, chromosome):
        if not self.is_loaded(chromosome):
            self.load_raw_table(chromosome)
        return self.allele_1[self.uniq_key(chromosome)]

    def get_allele_2(self, chromosome):
        if not self.is_loaded(chromosome):
            self.load_raw_table(chromosome)
        return self.allele_2[self.uniq_key(chromosome)]

    def find_index(self, query):
        bp = self.get_bp(str(query.get_chrom()))
        return (min(len(bp) - 1,
                    max(0,
                        bisect.bisect_left(bp, query.get_pos()))))

    def find_index_list(self, queries):
        return (np.array([self.find_index(query)
                          for query in queries]))

    def find_index_interval(self, query_l, query_r):
        """Find semi-open interval of indicies of SNPs that
           overlaps with a mapped fragment spanning [pos_l, pos_r)
        """
        index_l = self.find_index(query_l)
        index_r = self.find_index(query_r)
        return index_l, index_r

    def find_index_exact(self, query):
        """Find index on bim file by chromosom and position
        If there is an exact hit, return the index on bim
        else, return None to indicate it is NOT polymorphic region
        """
        index = self.find_index(query)
        return (index
                if (self.get_bp(query.get_chrom())[index] == query.get_pos())
                else None)

    def find_index_exact_list(self, queries):
        return (np.array([self.find_index_exact(query)
                          for query in queries]))
