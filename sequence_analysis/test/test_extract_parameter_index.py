""" extract parameter functions


"""
__revision__ = "$Id: test_extract_parameter_index.py 9401 2010-08-10 12:24:59Z cokelaer $"


from openalea.sequence_analysis.data_transform import IndexParameterExtract
from openalea.sequence_analysis.sequences import Sequences
from openalea.sequence_analysis import get_shared_data

def test1():
    """FIXME markovian_sequences call"""
    seq69 = Sequences(get_shared_data("pin_laricio_7x.seq"))
    a = IndexParameterExtract(seq69, 1929)
    b = IndexParameterExtract(seq69, 1929, 1994)
    c = seq69.index_parameter_extract(1929, -1).markovian_sequences()
    d = seq69.index_parameter_extract(1929, 1994).markovian_sequences()
    assert str(a) == str(c)
    assert str(b) == str(d)


if __name__ == "__main__":
    test1()

