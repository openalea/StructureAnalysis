"""Difference tests

.. author:: Thomas Cokelaer, Thomas.Cokelaer@inria.fr


"""

__revision__ = "$Id$"

import pytest

from openalea.sequence_analysis.data_transform import Difference
from openalea.sequence_analysis import Sequences
from .tools import robust_path as get_shared_data


@pytest.fixture
def build_seq1():
    return Sequences(str(get_shared_data("sequences1.seq")))


@pytest.fixture
def build_seqn():
    return Sequences(str(get_shared_data("sequences2.seq")))


def test_difference1(build_seq1):
    """difference test to finalise"""
    data = build_seq1
    res = Difference(data, 1)
    assert str(res) == str(data.difference(1, False))
    assert res.cumul_length == 50


def test_difference1_first_element(build_seq1):
    """difference test to finalise"""
    data = build_seq1
    res = Difference(data, 1, True)
    assert str(res) == str(data.difference(1, True))
    assert res.cumul_length == 52


def test_differencen(build_seqn):
    """difference test to finalise"""
    data = build_seqn
    res = Difference(data, Variable=1)
    assert str(res) == str(data.difference(1, False))
    assert res.cumul_length == 23
