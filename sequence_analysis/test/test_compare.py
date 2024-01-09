"""tests on comparison methods

.. author:: Thomas Cokelaer, Thomas.Cokelaer@inria.fr

.. todo:: systematic tests
"""
__revision__ = "$Id: test_compare.py 9401 2010-08-10 12:24:59Z cokelaer $"

from openalea.sequence_analysis.sequences import Sequences
from openalea.sequence_analysis.hidden_semi_markov import HiddenSemiMarkov
from openalea.sequence_analysis.compare import Compare
from openalea.stat_tool.vectors import VectorDistance, Vectors
from openalea.sequence_analysis.compare import Compare
from openalea.stat_tool.data_transform import ValueSelect, ExtractHistogram
from openalea.sequence_analysis.estimate import Estimate
from openalea.sequence_analysis.data_transform import Thresholding

from openalea.sequence_analysis import get_shared_data
from tools import runTestClass


class _Compare():
    """
    a main class to test the Compare function on different type of
    data structure.
    """
    def __init__(self):
        self.data = None

    def create_data(self):
        raise NotImplemented

    def test_compare(self):
        raise NotImplemented



class Test_Compare_Histograms(_Compare):
    def __init__(self):
        _Compare.__init__(self)
        self.data = self.create_data()

    def create_data(self):
        seq0 = Sequences(get_shared_data( "chene_sessile_15pa.seq"))
        vec10 = Vectors(seq0)
        vec95 = ValueSelect(vec10, 1, 95)
        vec96 = ValueSelect(vec10, 1, 96)
        vec97 = ValueSelect(vec10, 1, 97)
        return [vec95, vec96, vec97]

    def test_compare(self):
        seq = self.data
        res = Compare(ExtractHistogram(seq[0], 2),
                      ExtractHistogram(seq[1], 2),
                      ExtractHistogram(seq[2], 2), "N")
        assert res


class Test_Compare_Sequences(_Compare):

    def __init__(self):
        _Compare.__init__(self)
        self.data = self.create_data()

    def create_data(self):
        seq = Sequences(get_shared_data( 'dupreziana_a1.seq'))
        return seq

    def test_compare(self):
        seq = self.data
        matrix20 = Compare(seq)
        assert matrix20


class Test_Compare_Sequences_VectorDistance(_Compare):

    def __init__(self):
        _Compare.__init__(self)
        self.data = self.create_data()

    def create_data(self):
        seq = Sequences(get_shared_data( 'dupreziana_a1.seq'))
        return seq

    def test_compare(self):
        seq = self.data
        matrix20 = Compare(seq, VectorDistance("N", "N"))
        assert matrix20


class Test_Compare_Vectors_VectorDistance(_Compare):

    def __init__(self):
        _Compare.__init__(self)
        self.data = self.create_data()

    def create_data(self):
        data = Vectors([[1, 2, 3], [1, 3, 1], [4, 5, 6]])
        return data

    def test_compare(self):
        data = self.data
        a =  Compare(data, VectorDistance("N", "N", "N"))

        assert a.nb_row == 3
        assert a.nb_column == 3
        


class Test_Compare_hsmc_with_sequences(_Compare):

    def __init__(self):
        _Compare.__init__(self)
        self.data = self.create_data()

    def create_data(self):
        hsmc0 = HiddenSemiMarkov(get_shared_data( "belren1.hsc"))
        hsmc1 = HiddenSemiMarkov(get_shared_data("elstar1.hsc"))
        seq0 = Sequences(get_shared_data( "belren1.seq"))
        seq1 = Sequences(get_shared_data( "elstar1.seq"))
        data0 = Estimate(seq0, "HIDDEN_SEMI-MARKOV", hsmc0)
        data1 = Estimate(seq1, "HIDDEN_SEMI-MARKOV", hsmc1)
        return [seq0, seq1, data0, data1]

    def test_compare(self):
        data = self.data
        matrix20 = Compare(Thresholding(data[2], MinProbability=0.001), data[0],
                           Thresholding(data[3], MinProbability=0.001), data[1],
                           10000)
        assert matrix20


if __name__ ==  "__main__":
    runTestClass(Test_Compare_Histograms())
    runTestClass(Test_Compare_Sequences())
    runTestClass(Test_Compare_Sequences_VectorDistance())
    runTestClass(Test_Compare_Vectors_VectorDistance())
    runTestClass(Test_Compare_hsmc_with_sequences())
