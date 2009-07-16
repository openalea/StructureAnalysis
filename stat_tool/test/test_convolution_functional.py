"""Convolution tests

.. author:: Thomas Cokelaer, Thomas.Cokelaer@inria.fr

"""
__revision__ = "$Id: test_convolution.py 6609 2009-07-06 09:31:24Z cokelaer $"


from openalea.stat_tool import *

from tools import interface




def test1():

    convol1 = Convolution(Distribution("B", 0, 10, 0.5),
        Distribution("NB", 0, 10, 0.5))
    convol1 = Convolution('data/convolution1.conv')
  
    convol_histo1 = Simulate(convol1, 200)

    histo20 = ExtractHistogram(convol_histo1, "Elementary", 1)

    convol2 = Estimate(convol_histo1, "CONVOLUTION", 
                       ExtractDistribution(convol1, "Elementary", 1),
                       MinInfBound=0)

    histo21 = ExtractHistogram(ExtractData(convol2), 'Elementary', 1)
    histo22 = ToHistogram(ExtractDistribution(convol2, 'Elementary', 1))

    Plot(histo21, histo22)

    histo_b2 = Histogram("data/nothofagus_antarctica_bud_2.his")
    histo_s2 = Histogram("data/nothofagus_antarctica_shoot_2.his")

    # Estimator="Likelihood" (default) / "PenalizedLikelihood" / "Parametric"
    # Si Estimator="PenalizedLikelihood", options supplementaires possibles
    # Penalty="FirstDifference" / "SecondDifference" (default) / "Entropy", Weight,
    # Outside="Zero" (default) / "Continuation" (cf. stat_funs4.cpp).

    # "NP" or "NON_PARAMETRIC" (for estimation only)

    convol30 = Estimate(Shift(histo_s2, 1), "CONVOLUTION", 
                        Estimate(histo_b2, "NP"), NbIteration=500)
    convol31 = Estimate(Shift(histo_s2, 1), "CONVOLUTION", 
                        Estimate(histo_b2, "NP"), NbIteration=100,
                        Estimator="PenalizedLikelihood", Weight=0.5)
    convol32 = Estimate(Shift(histo_s2, 1), "CONVOLUTION",
                       Estimate(histo_b2, "NP"),
                        Estimator="Parametric")

    # FIX Cast 
    Plot(convol30, convol31, convol32)
    Display(convol31)
    #if DISABLE_PLOT==False:
    Plot(convol31)
    Plot(ExtractDistribution(convol31, "Convolution"))
    Save(convol31, "data/nothofagus_antartica_2.xls", Format="SpreadSheet")


if __name__=="__main__":
    test1()
