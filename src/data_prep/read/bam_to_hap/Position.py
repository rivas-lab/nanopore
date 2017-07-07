class Position:
    """Class to store genomic position (contig name + position on the contig)
    """

    def __init__(self, chrom, pos):
        if isinstance(chrom, int):
            self.chrom = str(chrom)
        elif chrom.startswith('chr'):
            self.chrom = str(chrom[3:])
        else:
            self.chrom = str(chrom)

        if isinstance(pos, int):
            self.pos = pos
        else:
            self.pos = int(pos)

    def get_pos(self):
        return self.pos

    def get_chrom(self):
        return self.chrom

    def __str__(self):
        return '{}:{}'.format(self.chrom, self.pos)

    def __eq__(self, other):
        return ((self.chrom == other.chrom)
                and (self.pos == other.pos))

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash((self.chrom, self.pos))
