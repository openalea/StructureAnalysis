"""Compound tests"""
__revision__ = "$Id: test_compound.py 6609 2009-07-06 09:31:24Z cokelaer $"

from openalea.stat_tool import *

from tools import interface


    
def test1():
    """tests on compound extract from stat_tool_test.aml
    """

    cdist1 = Compound("data/compound1.cd")

    chisto1 = Simulate(cdist1, 200)

    histo30 = ExtractHistogram(chisto1, "Sum")

    cdist2 = Estimate(chisto1, "COMPOUND", ExtractDistribution(cdist1, "Elementary"), "Sum", MinInfBound=0)

    histo31 = ExtractHistogram(ExtractData(cdist2), "Sum")
    histo32 = ToHistogram(ExtractDistribution(cdist2, "Sum"))
    Plot(histo31, histo32)

    peup1 = Histogram("data/peup1.his")
    mixt4 = Estimate(peup1, "MIXTURE", "B", "NB")
    histo33 = ToHistogram(ExtractDistribution(mixt4, "Component", 2))
    histo34 = Shift(histo33, -11)
    Plot(histo34)

    cdist3 = Estimate(histo34, "COMPOUND", Distribution("B", 0, 1, 0.7), "Sum")
    Plot(cdist3)
    cdist3 = Estimate(histo34, "COMPOUND", Distribution("B", 0, 1, 0.7), "Sum", MinInfBound=0)


if __name__=="__main__":
    test1()
