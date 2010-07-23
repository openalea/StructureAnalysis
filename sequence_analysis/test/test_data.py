from openalea.sequence_analysis.data import path



def test_read_data():
    from openalea.sequence_analysis.sequences import Sequences
    seq = Sequences(path + 'sequences_tutorial.dat')
