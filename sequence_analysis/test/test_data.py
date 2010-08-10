from openalea.sequence_analysis import get_shared_data


def test_get_shared_data():
    from openalea.sequence_analysis.sequences import Sequences
    seq = Sequences(get_shared_data('wij1.seq'))
