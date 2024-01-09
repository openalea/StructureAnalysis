"""vectors functional test from stat_tool_test.aml

#######################################################################
#
#  vectors : contingency table, one-way variance analysis,
#  linear regression or nonparametric regression (loess smoother)
#
#  Oak trunk annual shoots
#
#  INDEX_PARAMETER : TIME  (year of growth - 1995, 1996, 1997)
#  VARIABLE 1 : length of the annual shoot (cm)
#  VARIABLE 2 : diameter of the annual shoot (1/10 de mm)
#  VARIABLE 3 : number of cycles
#  VARIABLE 4 : number of nodes
#  VARIABLE 5 : number de branches
#
#######################################################################
"""
__version__ = "$Id: test_vectors_functional.py 7886 2010-02-09 07:49:53Z cokelaer $"

from openalea.stat_tool.vectors import Vectors, \
    VarianceAnalysis, ContingencyTable
from openalea.stat_tool.data_transform import ValueSelect
from openalea.stat_tool.regression import Regression
from openalea.stat_tool.data_transform import ExtractHistogram, \
     SelectVariable, ExtractDistribution
from openalea.stat_tool.comparison import Compare
from openalea.stat_tool.estimate import Estimate
from openalea.stat_tool.output import Display, Plot

def test():
    vec10 = Vectors("data/chene_sessile.vec")
    Plot(vec10)
    # plot of the pointwise averages
    Plot(Regression(vec10, "MovingAverage", 1, 2, [1]))
    
    vec95 = ValueSelect(vec10, 1, 1995)    
    vec96 = ValueSelect(vec10, 1, 1996)
    vec97 = ValueSelect(vec10, 1, 1997)
    
    VarianceAnalysis(vec10, 1, 2, "N")
    Compare(ExtractHistogram(vec95, 2), ExtractHistogram(vec96, 2), \
            ExtractHistogram(vec95, 2), "N")
    
    
    Plot(ExtractHistogram(vec95, 2), ExtractHistogram(vec96, 2), \
         ExtractHistogram(vec97, 2))
    
    ContingencyTable(vec10, 1, 4)
    
    
    # one-way variance analysis based on ranks
    VarianceAnalysis(vec10, 1, 4, "O")
    Compare(ExtractHistogram(vec95, 4), ExtractHistogram(vec96, 4), \
            ExtractHistogram(vec95, 4), "O")
    
    # looks like it is not plotted
    Plot(ExtractHistogram(vec95, 4), ExtractHistogram(vec96, 4),
        ExtractHistogram(vec97, 4))
    Plot(ExtractHistogram(vec95, 5), ExtractHistogram(vec96, 5),
        ExtractHistogram(vec97, 5))
    Plot(ExtractHistogram(vec95, 6), ExtractHistogram(vec96, 6),
        ExtractHistogram(vec97, 6))
    
    
    vec11 = ValueSelect(vec10, 4, 1)
    vec12 = ValueSelect(vec10, 4, 2)
    vec13 = ValueSelect(vec10, 4, 3, 4)
    
    Plot(ExtractHistogram(vec11, 2), ExtractHistogram(vec12, 2),
            ExtractHistogram(vec13, 2))
    Plot(ExtractHistogram(vec11, 5), ExtractHistogram(vec12, 5),
        ExtractHistogram(vec13, 5))
    
    mixt20 = Estimate(ExtractHistogram(vec10, 2), \
                      "MIXTURE", "NB", "NB", "NB", "NB", \
                      NbComponent="Estimated")
    Display(mixt20)
    
    Plot(mixt20)
    Plot(ExtractDistribution(mixt20, "Mixture"))
    
    _mixt21 = Estimate(ExtractHistogram(vec10, 5), \
                       "MIXTURE", "NB", "NB", "NB", "NB", \
                       NbComponent="Estimated")
    
    vec9596 = ValueSelect(vec10, 1, 1995, 1996)
    
    Plot(ExtractHistogram(ValueSelect(vec9596, 4, 1), 6), \
         ExtractHistogram(ValueSelect(vec9596, 4, 2), 6), \
         ExtractHistogram(ValueSelect(vec9596, 4, 3, 4), 6))
    
    # linear regression
    regress10 = Regression(vec10, "Linear", 5, 2)
    Display(regress10)
    Plot(regress10)
    
    # nonparametric regression (loess smoother)
    
    _regress11 = Regression(vec10, "NearestNeighbors", 5, 2, 0.3)
    _regress12 = Regression(vec9596, "Linear", 5, 6)
    _regress13 = Regression(vec9596, "NearestNeighbors", 5, 6, 0.5)
    
    _vec15 = SelectVariable(vec10, [1, 3, 6], Mode="Reject")
    

if __name__ == "__main__":
    test1()
