"""tests on comparison methods

.. author:: Thomas Cokelaer, Thomas.Cokelaer@inria.fr

.. todo:: systematic tests
"""

__revision__ = "$Id$"

import pytest

from openalea.sequence_analysis.sequences import Sequences
from openalea.sequence_analysis.hidden_semi_markov import HiddenSemiMarkov
from openalea.sequence_analysis.compare import Compare
from openalea.stat_tool.vectors import VectorDistance, Vectors
from openalea.sequence_analysis.compare import Compare
from openalea.stat_tool.data_transform import ValueSelect, ExtractHistogram
from openalea.sequence_analysis.estimate import Estimate
from openalea.sequence_analysis.data_transform import Thresholding

from .tools import runTestClass, robust_path as get_shared_data


class _Compare:
    """
    a main class to test the Compare function on different type of
    data structure.
    """

    def __init__(self):
        self.data = None

    def create_data(self):
        raise NotImplemented

    def test_compare(self):
        raise NotImplemented


@pytest.fixture
def create_data_compare_histogram():
    seq0 = Sequences(str(get_shared_data("chene_sessile_15pa.seq")))
    vec10 = Vectors(seq0)
    vec95 = ValueSelect(vec10, 1, 95)
    vec96 = ValueSelect(vec10, 1, 96)
    vec97 = ValueSelect(vec10, 1, 97)
    return [vec95, vec96, vec97]


def test_compare_compare_histogram(create_data_compare_histogram):
    seq = create_data_compare_histogram
    res = Compare(
        ExtractHistogram(seq[0], 2),
        ExtractHistogram(seq[1], 2),
        ExtractHistogram(seq[2], 2),
        "N",
    )
    assert res


@pytest.fixture
def create_data_sequences():
    return Sequences(str(get_shared_data("dupreziana_a1.seq")))


def test_compare_sequences(create_data_sequences):
    seq = create_data_sequences
    matrix20 = Compare(seq)
    assert matrix20


@pytest.fixture
def create_data_vectordistance():
    return Sequences(str(get_shared_data("dupreziana_a1.seq")))


def test_compare_vectordistance(create_data_vectordistance):
    seq = create_data_vectordistance
    matrix20 = Compare(seq, VectorDistance("N", "N"))
    assert matrix20


@pytest.fixture
def create_data_vectors_vectordistance():
    return Vectors([[1, 2, 3], [1, 3, 1], [4, 5, 6]])


def test_compare_vectors_vectordistance(create_data_vectors_vectordistance):
    data = create_data_vectors_vectordistance
    a = Compare(data, VectorDistance("N", "N", "N"))
    assert a.nb_row == 3
    assert a.nb_column == 3


@pytest.fixture
def create_data_compare_hsmc_with_sequences():
    hsmc0 = HiddenSemiMarkov(str(get_shared_data("belren1.hsc")))
    hsmc1 = HiddenSemiMarkov(str(get_shared_data("elstar1.hsc")))
    seq0 = Sequences(str(get_shared_data("belren1.seq")))
    seq1 = Sequences(str(get_shared_data("elstar1.seq")))
    data0 = Estimate(seq0, "HIDDEN_SEMI-MARKOV", hsmc0)
    data1 = Estimate(seq1, "HIDDEN_SEMI-MARKOV", hsmc1)
    return [seq0, seq1, data0, data1]


def test_compare_compare_hsmc_with_sequences(create_data_compare_hsmc_with_sequences):
    data = create_data_compare_hsmc_with_sequences
    matrix20 = Compare(
        Thresholding(data[2], MinProbability=0.001),
        data[0],
        Thresholding(data[3], MinProbability=0.001),
        data[1],
        10000,
    )
    assert matrix20


class Test_Compare_hsmc_with_sequences(_Compare):
    def __init__(self):
        _Compare.__init__(self)
        self.data = self.create_data()

    def create_data(self):
        hsmc0 = HiddenSemiMarkov(str(get_shared_data("belren1.hsc")))
        hsmc1 = HiddenSemiMarkov(str(get_shared_data("elstar1.hsc")))
        seq0 = Sequences(str(get_shared_data("belren1.seq")))
        seq1 = Sequences(str(get_shared_data("elstar1.seq")))
        data0 = Estimate(seq0, "HIDDEN_SEMI-MARKOV", hsmc0)
        data1 = Estimate(seq1, "HIDDEN_SEMI-MARKOV", hsmc1)
        return [seq0, seq1, data0, data1]

    def test_compare(self):
        data = self.data
        matrix20 = Compare(
            Thresholding(data[2], MinProbability=0.001),
            data[0],
            Thresholding(data[3], MinProbability=0.001),
            data[1],
            10000,
        )
        assert matrix20
