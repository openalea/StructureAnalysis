"""comparison test

:Author: Thomas Cokelaer, Thomas.Cokelaer@inria.fr

histo done
vector done

.. todo:: sequence and markov
"""
__version__ = "$Id: test_comparison.py 9380 2010-08-06 17:53:45Z cokelaer $"


import os

from openalea.stat_tool.histogram import Histogram
from openalea.stat_tool.vectors import Vectors, VectorDistance
from openalea.stat_tool.data_transform import SelectVariable
from openalea.stat_tool.comparison import Compare, ComparisonTest

from tools import runTestClass

from openalea.stat_tool import get_shared_data
class TestHisto():
    """a simple unittest class"""

    def __init__(self):
        self.meri1 = Histogram(get_shared_data("meri1.his"))
        self.meri2 = Histogram(get_shared_data("meri1.his"))
        self.meri3 = Histogram(get_shared_data("meri1.his"))
    def test_comparisontest(self):
        meri1 = self.meri1
        meri2 = self.meri2

        assert ComparisonTest("F", meri1, meri2) == meri1.f_comparison(meri2)
        assert ComparisonTest("T", meri1, meri2) == meri1.t_comparison(meri2)
        assert ComparisonTest("W", meri1, meri2) == meri1.wmw_comparison(meri2)

    def test_comparison_histo(self):
        # check both the long and short argument (O, S, N)
        meri1 = self.meri1
        meri2 = self.meri2
        meri3 = self.meri3

        c1 = Compare(meri1, meri2, meri3, 'N')
        c2 = Compare(meri1, meri2, meri3, 'O')
        c3 = Compare(meri1, meri2, meri3, 'S')

        c1_long = Compare(meri1, meri2, meri3, 'NUMERIC')
        c2_long = Compare(meri1, meri2, meri3, 'ORDINAL')
        c3_long = Compare(meri1, meri2, meri3, 'SYMBOLIC')

        assert c1 == c1_long
        assert c2 == c2_long
        assert c3 == c3_long

        assert meri1.compare_histo(meri2, meri3, "N") == c1
        assert meri1.compare_histo(meri2, meri3, "S") == c3
        assert meri1.compare_histo(meri2, meri3, "O") == c2

        assert meri1.compare_histo(meri2, meri3, "N") == c1_long
        assert meri1.compare_histo(meri2, meri3, "S") == c3_long
        assert meri1.compare_histo(meri2, meri3, "O") == c2_long

    def test_comparison_wrong_argument_1(self):
        meri1 = self.meri1
        try:
            Compare("N")
            assert False
        except NotImplementedError:
            assert True

    def test_comparison_wrong_argument_2(self):
        meri1 = self.meri1
        meri2 = self.meri2
        try:
            meri1.compare_histo(meri1, meri2, "F")
            assert False
        except KeyError:
            assert True

    def test_comparisontest_wrong_argument(self):
        meri1 = self.meri1
        meri2 = self.meri2
        try:
            ComparisonTest("N", meri1, meri2)
            assert False
        except TypeError:
            assert True

    def test_comparison_histo_filename(self):
        meri1 = self.meri1
        meri2 = self.meri2
        meri3 = self.meri3

        _c1 = Compare(meri1, meri2, meri3, 'N', Filename='result.dat')
        os.remove('result.dat')
        _c1 = Compare(meri1, meri2, meri3, 'N', Filename='result.dat', \
                     Format='ASCII')
        os.remove('result.dat')
        _c1 = Compare(meri1, meri2, meri3, 'N', Filename='result.dat', \
                     Format='SpreadSheet')
        os.remove('result.dat')
        try:
            _c1 = Compare(meri1, meri2, meri3, 'N', Filename='result.dat', \
                         Format='badname')
            assert False
        except ValueError:
            assert True


class TestVectors():

    def __init__(self):
        pass

    def test_compare_vectors(self):

        vec10 = Vectors(get_shared_data("chene_sessile.vec"))
        vec15 = SelectVariable(vec10, [1, 3, 6], Mode="Reject")
        assert vec15

        matrix10 = Compare(vec15, VectorDistance("N", "N", "N"))
        assert matrix10
        assert str(vec15.compare(VectorDistance("N", "N", "N"),
                                 True)) == str(matrix10)



if __name__ == "__main__":
    runTestClass(TestHisto())
    runTestClass(TestVectors())
