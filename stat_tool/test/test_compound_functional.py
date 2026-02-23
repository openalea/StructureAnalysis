"""Compound functional test extracted from original version
of stat_tool_test.aml


"""
__version__ = "$Id$"

import os

import numpy as np

try:
    from .tools import interface
    from .tools import robust_path as get_shared_data
except ImportError:
    from tools import interface
    from tools import robust_path as get_shared_data


from openalea.stat_tool import Compound, Histogram, Distribution

from openalea.stat_tool import (
    Simulate, ExtractHistogram, ToHistogram,
    Estimate, ExtractData, ExtractDistribution,
    Plot, Shift
)

# read data file
def test():
    # cdist1 = Compound("data/compound1.cd")
    # Plot(cdist1)
    #
    # chisto1 = Simulate(cdist1, 200)
    # Plot(chisto1)
    #
    # histo30 = ExtractHistogram(chisto1, "Sum")
    # Plot(histo30)
    #
    # histo30e = ExtractHistogram(chisto1, "Elementary")
    # Plot(histo30e)
    #
    #
    # cdist2 = Estimate(chisto1, "COMPOUND",
    #                   ExtractDistribution(cdist1, "Elementary"),
    #                   "Sum", MinInfBound=0)
    #
    # histo31 = ExtractHistogram(ExtractData(cdist2), "Sum")
    # histo32 = ToHistogram(ExtractDistribution(cdist2, "Sum"))
    # Plot(histo31, histo32)

    peup1 = Histogram(get_shared_data("peup1.his"))
    mixt4 = Estimate(peup1, "MIXTURE", "B", "NB")
    histo33 = ToHistogram(ExtractDistribution(mixt4, "Component", 2))
    histo33array = np.array(histo33)
    bad_shift_val = int(np.min(np.where(histo33array > 0))) + 1
    try:
        histo34 = Shift(histo33, -bad_shift_val)
    except Exception:
        pass
    else:
        raise ValueError("Failed to detect bad shift value")
    histo34 = Shift(histo33, -4)
    Plot(histo34)

    cdist3 = Estimate(histo34, "COMPOUND",
                      Distribution("B", 0, 1, 0.7), "Sum")
    Plot(cdist3)
    cdist3 = Estimate(histo34, "COMPOUND",
                       Distribution("B", 0, 1, 0.7), "Sum", MinInfBound=0)

    Plot(cdist3)


if __name__ == "__main__":
    test()
