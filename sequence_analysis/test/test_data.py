from openalea.sequence_analysis.data import path
from os.path import join as pj


def test_read_data():
    from openalea.sequence_analysis.sequences import Sequences
    seq = Sequences(pj(path , 'sequences_tutorial.dat'))
