""" Test moving average

.. author:: Thomas Cokelaer, Thomas.Cokelaer@inria.fr

.. todo:: to be done
"""
__revision__ = "$Id$"

from openalea.sequence_analysis.data_transform import MovingAverage
from openalea.stat_tool.distribution import Distribution
from data import seq70
from tools import runTestClass


class Test():
    def __init__(self):
        self.data = seq70

    def test_distribution(self):
        seq70 = self.data
        MovingAverage(seq70, Distribution("B", 0, 16, 0.5), BeginEnd=True)
        MovingAverage(seq70, Distribution("B", 0, 16, 0.5), BeginEnd=True, Output="Residual")

    def test_frequencies(self):
        seq70 = self.data
        MovingAverage(seq70, [1, 1, 1], BeginEnd=True)
        MovingAverage(seq70, [1, 1, 1], BeginEnd=True, Output="Residual")

    def test_filter(self):
        """test not yet implemented"""
        pass

if __name__ == "__main__":
    runTestClass(Test())