""" extract parameter functions


"""
__revision__ = "$Id$"


from openalea.sequence_analysis.data_transform import IndexParameterExtract
from openalea.sequence_analysis.sequences import Sequences
from tools import runTestClass, robust_path as get_shared_data

def test1():
    """FIXME markovian_sequences call"""
    seq69 = Sequences(str(get_shared_data("pin_laricio_7x.seq")))
    a = IndexParameterExtract(seq69, 1929)
    b = IndexParameterExtract(seq69, 1929, 1994)
    c = seq69.index_parameter_extract(1929, -1).markovian_sequences()
    d = seq69.index_parameter_extract(1929, 1994).markovian_sequences()
    assert str(a) == str(c)
    assert str(b) == str(d)


if __name__ == "__main__":
    test1()

