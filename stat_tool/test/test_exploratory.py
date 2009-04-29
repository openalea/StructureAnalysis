"""Functional test"""
__revision__ = "$Id$"


from openalea.stat_tool.histogram import Histogram
from openalea.stat_tool.estimate import Estimate
from openalea.stat_tool.comparison import Compare, ComparisonTest
from openalea.stat_tool.data_transform import Merge

class Test:
    
    def test_aml_1(self):
        """ Test AML Functions """

        meri1 = Histogram("data/meri1.his")
        meri2 = Histogram("data/meri2.his")
        meri3 = Histogram("data/meri3.his")
        meri4 = Histogram("data/meri4.his")
        meri5 = Histogram("data/meri5.his")
        
        # Plot(meri1, meri2, meri3, meri4, meri5)
        
        Compare(meri1, meri2, meri3, meri4, meri5, "N")
        
        ComparisonTest("F", meri1, meri2)
        ComparisonTest("T", meri1, meri2)
        ComparisonTest("W", meri1, meri2)
        
        ComparisonTest("F", meri1, meri3)
        ComparisonTest("T", meri1, meri3)
        ComparisonTest("W", meri1, meri3)
        
        # estimation of a mixture of two distributions assuming a first 
        # sub-population of GUs made only of a preformed part and a second 
        # sub-population made of both a preformed part and a neoformed part
        
        _mixt1 = Estimate(meri2, "MIXTURE", "B", "B")
        
        meri = Merge(meri1, meri2, meri3, meri4, meri5)
        
        # model selection approach: estimation of both the mixture parameters 
        # and the number of components 
        
        _mixt2 = Estimate(meri, "MIXTURE", "B", "B", "B", "B",  
                          NbComponent="Estimated")
        # mixt2 = Estimate(meri, "MIXTURE", "NB", "NB")
        # Plot(ExtractDistribution(mixt2, "Mixture"))
        # Display(mixt2)
        
        peup1 = Histogram("data/peup1.his")
        peup2 = Histogram("data/peup2.his")
        peup3 = Histogram("data/peup3.his")
        peup4 = Histogram("data/peup4.his")
        peup5 = Histogram("data/peup5.his")
        peup6 = Histogram("data/peup6.his")
        
        _mixt10 = Estimate(peup2, "MIXTURE", "B", "NB", "NB", "NB", 
                           NbComponent="Estimated")
        
        peup = Merge(peup1, peup2, peup3, peup4, peup5, peup6)
        
        mixt11 = Estimate(peup, "MIXTURE", "B", "NB", "NB", "NB", 
                          NbComponent="Estimated")
        assert mixt11
        mixt11 = Estimate(peup, "MIXTURE", "B", "NB")
        assert mixt11
        
        
