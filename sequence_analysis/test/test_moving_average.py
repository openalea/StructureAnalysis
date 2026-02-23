"""Test moving average

.. author:: Thomas Cokelaer, Thomas.Cokelaer@inria.fr

.. todo:: to be done
"""

__revision__ = "$Id$"

import pytest

from openalea.sequence_analysis import Cluster, Distribution, MovingAverage, Sequences

from .tools import robust_path as get_shared_data


@pytest.fixture
def create_data_moving_average():
    seq = Sequences(get_shared_data("pin_laricio_7x.seq"))
    return Cluster(seq, "Step", 1, 10)


class Test:
    def test_distribution(self, create_data_moving_average):
        seq70 = create_data_moving_average
        MovingAverage(seq70, Distribution("B", 0, 16, 0.5), BeginEnd=True)
        MovingAverage(
            seq70, Distribution("B", 0, 16, 0.5), BeginEnd=True, Output="Residual"
        )

    def test_frequencies(self, create_data_moving_average):
        seq70 = create_data_moving_average
        MovingAverage(seq70, [1, 1, 1], BeginEnd=True)
        MovingAverage(seq70, [1, 1, 1], BeginEnd=True, Output="Residual")

    def test_filter(self):
        """test not yet implemented"""
        pass
