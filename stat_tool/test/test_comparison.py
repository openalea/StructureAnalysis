"""comparison test"""
__revision__ = "$Id: $"

from openalea.stat_tool.histogram import Histogram
from openalea.stat_tool.vectors import Vectors, VectorDistance
from openalea.stat_tool.data_transform import SelectVariable
from openalea.stat_tool.comparison import Compare, ComparisonTest


class Test:
    """a simple unittest class"""

    def test_comparisontest(self):

        meri1 = Histogram("meri1.his")
        meri2 = Histogram("meri2.his")

        assert ComparisonTest("F", meri1, meri2)
        assert ComparisonTest("T", meri1, meri2)
        assert ComparisonTest("W", meri1, meri2)

    def test_comparison(self):

        meri1 = Histogram("meri1.his")
        meri2 = Histogram("meri2.his")
        meri3 = Histogram("meri3.his")

        assert Compare(meri1, meri2, meri3, 'N')
        assert Compare(meri1, meri2, meri3, 'O')
        assert Compare(meri1, meri2, meri3, 'S')

    def test_compare_vectors(self):
        
        #from vectors import Vectors, VectorDistance
        #from data_transform import SelectVariable

        vec10 = Vectors("chene_sessile.vec")
        vec15 = SelectVariable(vec10, [1, 3, 6], Mode="Reject")
        assert vec15

        matrix10 = Compare(vec15, VectorDistance("N", "N", "N"))
        assert matrix10

    def test_compare_sequence(self):
        pass
        
    def test_compare_markov(self):
        pass



if __name__=="__main__":
    # perform all the test in the class Test (unit tests)
    test = Test()
    for method in dir(test):
        if method.startswith('_'):
            continue
        if callable(getattr(test, method)):
            getattr(test, method)()
        else:
            print 'skipping'
    # and functional tests.    


