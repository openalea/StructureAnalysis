""" Test moving average

.. author:: Thomas Cokelaer, Thomas.Cokelaer@inria.fr

.. todo:: to be done
"""
__revision__ = "$Id:  $"

from openalea.sequence_analysis.data_transform import MovingAverage
from openalea.stat_tool.distribution import Distribution
from data import seq70

def test_distribution():
    MovingAverage(seq70, Distribution("B", 0, 16, 0.5), BeginEnd=True)
    MovingAverage(seq70, Distribution("B", 0, 16, 0.5), BeginEnd=True, Output="Residual")

def test_frequencies():
    MovingAverage(seq70, [1, 1, 1], BeginEnd=True)
    MovingAverage(seq70, [1, 1, 1], BeginEnd=True, Output="Residual")

def test_filter():
    """test not yet implemented"""
    pass
