"""Convolution functional test extracted from original version of
stat_tool_test.aml"""
__version__ = "$Id: test_convolution_functional.py 8552 2010-03-25 15:06:07Z cokelaer $"


from openalea.stat_tool import Convolution, Histogram, Distribution
from openalea.stat_tool import Simulate, ExtractHistogram, ToHistogram
from openalea.stat_tool import Estimate, ExtractData, ExtractDistribution
from openalea.stat_tool import Plot, Shift, Display, Save


def test():
    convol2 = Convolution(Distribution("B", 0, 10, 0.5),
                          Distribution("NB", 0, 10, 0.5))
    Plot(convol2, Title='convol2')

    convol1 = Convolution('data/convolution1.conv')
    Plot(convol1, Title='Convol1 (from file data/convolution.conv)')

    Plot(convol1, convol2, Title='The two input convolutions')

    convol_histo1 = Simulate(convol1, 200)
    Plot(convol_histo1, Title='convol_histo1 from simulate')

    histo20 = ExtractHistogram(convol_histo1, "Elementary", 1)
    Plot(histo20)

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
    # Penalty="FirstDifference" / "SecondDifference" (default) / "Entropy",
    # Weight,
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


    Plot(convol30, convol31, convol32)
    Display(convol31)

    Plot(convol31)
    Plot(ExtractDistribution(convol31, "Convolution"))
    Save(convol31, "data/nothofagus_antartica_2.xls", Format="SpreadSheet")


if __name__ == "__main__":
    test()
