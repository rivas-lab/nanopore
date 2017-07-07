import enum

@enum.unique
class Nucleotide(enum.Enum):
    N = -1
    A = 0
    C = 1
    G = 2
    T = 3
