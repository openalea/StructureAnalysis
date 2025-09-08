from tools import runTestClass, robust_path as get_shared_data


def test_get_shared_data():
    from openalea.sequence_analysis.sequences import Sequences
    seq = Sequences(str(get_shared_data('wij1.seq')))
    assert(seq)

if __name__ == "__main__":
    test_get_shared_data()