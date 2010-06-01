from openalea.sequence_analysis import data

def test():
    path = data.path
    seq = data.data_tutorial.seq1

def test_read_data():
    path = data.path
    from openalea.sequence_analysis.sequences import Sequences
    seq = Sequences(path + 'sequences_tutorial.dat')
