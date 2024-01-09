""" Test moving average

.. author:: Thomas Cokelaer, Thomas.Cokelaer@inria.fr

.. todo:: to be done
"""
__revision__ = "$Id: test_moving_average.py 9885 2010-11-06 18:19:34Z cokelaer $"

from openalea.sequence_analysis import *
from openalea.stat_tool.distribution import Distribution
from tools import runTestClass

seq = Sequences(get_shared_data("pin_laricio_7x.seq"))
seq70 = Cluster(seq, "Step", 1, 10)


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
