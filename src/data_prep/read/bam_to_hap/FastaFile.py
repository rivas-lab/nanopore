import pysam


class FastaFile:
    """A class to load and manage collection of FASTA files
    """
    def __init__(self, filename, fast_ext='fa'):
        self.filename = filename
        self.fast_ext = fast_ext
        self.seq = {}

    def load(self, chromosome):
        self.seq[self.uniq_key(chromosome)] = \
            pysam.FastaFile('{}{}.{}'.format(self.filename, chromosome, self.fast_ext))

    def uniq_key(self, chromosome):
        return str(chromosome)

    def is_loaded(self, chromosome):
        return self.uniq_key(chromosome) in self.seq

    def get_seq(self, chromosome):
        if not self.is_loaded(chromosome):
            self.load(chromosome)
        return self.seq[self.uniq_key(chromosome)]
