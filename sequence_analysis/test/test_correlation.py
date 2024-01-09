""" Tests on ComputeAutoCorrelation, ComputeParialAutoCorrelation,
ComputewhiteNoiseCorrelation

.. author:: Thomas Cokelaer, Thomas.Cokelaer@inria.fr

"""
__revision__ = "$Id: test_correlation.py 9401 2010-08-10 12:24:59Z cokelaer $"


from openalea.sequence_analysis.sequences import Sequences
from openalea.sequence_analysis.data_transform import MovingAverage, VariableScaling
from openalea.sequence_analysis.correlation import ComputeCorrelation
from openalea.sequence_analysis.correlation import type_dict
from openalea.sequence_analysis.correlation import ComputeWhiteNoiseCorrelation
from openalea.sequence_analysis.correlation import ComputePartialAutoCorrelation

from openalea.stat_tool.distribution import Distribution
from tools import runTestClass
from openalea.sequence_analysis import get_shared_data

class Data():

    def __init__(self):

        self.sequence = self.create_sequence_data()
        self.type_map = type_dict

    def create_sequence_data(self):

        seq66 = Sequences(get_shared_data( "laricio_date66.seq"))
        seq69 = MovingAverage(VariableScaling(seq66, 3, 100),
                          Distribution("B", 0, 6, 0.5), BeginEnd=True,
                          Output="Residual")
        return seq69

def CorrelationData(index=1):
    """Returns a correlation

    index from 1 to 3"""
    seq66 = Sequences(get_shared_data( "laricio_date66.seq"))
    ret = ComputeCorrelation(seq66, index)
    return ret




class TestComputeCorrelation(Data):

    def __init__(self):
        Data.__init__(self)
        self.variable = 2
    def compute_correlation_type(self, variable, type, MaxLag=10,
                                 Normalization="Exact"):
        seq = self.sequence
        cf = ComputeCorrelation(seq, variable,
                               Type=type, MaxLag=MaxLag, Normalization=Normalization)
        assert cf.type == self.type_map[type]
        return cf

    def test_correlation_no_optional_arguments(self):
        seq = self.sequence
        cf = ComputeCorrelation(seq, self.variable)

    def test_spearman(self):
        seq = Sequences(get_shared_data( "laricio_date66.seq"))
        ComputeCorrelation(seq, 1, Type="Spearman")
        ComputeCorrelation(seq, 1, 2,Type="Spearman")
        try:
            dummy = 3
            ComputeCorrelation(seq, 1, 2, dummy, Type="Spearman")
        except:
            assert True

    def test_pearson(self):
        # used by test_ComputewhiteNoiseCorrelation
        return self.compute_correlation_type(self.variable, "Pearson", 10)

    def test_norm1(self):
        self.compute_correlation_type(self.variable, "Pearson", 10, "Approximated")

    def test_norm2(self):
        self.compute_correlation_type(self.variable, "Pearson", 10, "Exact")

    def test_norm3(self):
        try:
            self.compute_correlation_type(self.variable, "Pearson", 10, "Typolabel")
            assert False
        except Exception:
            assert True



class TestComputeWhiteNoiseCorrelation(TestComputeCorrelation):

    def __init__(self):
        TestComputeCorrelation.__init__(self)
        self.correlation = self.test_pearson()

    def test_filter(self):
        data = self.correlation
        ComputeWhiteNoiseCorrelation(data, [1, 1, 1])

    def test_order(self):
        data = self.correlation
        ComputeWhiteNoiseCorrelation(data, 1)

    def test_distribution(self):
        data = self.correlation
        ComputeWhiteNoiseCorrelation(data , Distribution("BINOMIAL", 0,4,0.5))

class TestComputePartialAutoCorrelation(Data):

    def __init__(self):
        Data.__init__(self)

    def test_compute_partial_auto_correlation(self):
        ComputePartialAutoCorrelation(self.sequence, 2, MaxLag=5)
        ComputePartialAutoCorrelation(self.sequence, MaxLag=5)


if __name__ == "__main__":
    runTestClass(TestComputeCorrelation())
    runTestClass(TestComputePartialAutoCorrelation())
    runTestClass(TestComputeWhiteNoiseCorrelation())
