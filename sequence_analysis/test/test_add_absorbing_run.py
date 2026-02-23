"""tests on the method AddAbsorbingRun

.. author:: Thomas Cokelaer, Thomas.Cokelaer@inria.fr

.. todo:: markov case ?
"""

__revision__ = "$Id$"

import pytest

from dataclasses import dataclass
from typing import Any

from openalea.sequence_analysis.sequences import Sequences
from openalea.sequence_analysis.semi_markov import SemiMarkov
from openalea.sequence_analysis.data_transform import AddAbsorbingRun
from openalea.sequence_analysis import get_shared_data

MAX_RUN_LENGTH = 20  # hardcoded values in CPP code


@dataclass
class AbsorbingData:
    data: Any
    max_length: float
    max_run_length: float


@pytest.fixture(params=["raw", "sequences", "semimarkov"])
def AddAbsorbingRunData(request):
    if request.param == "raw":
        return AbsorbingData(None, -1, 20)
    elif request.param == "sequences":
        return Sequences(str(get_shared_data("sequences1.seq")))
    elif request.param == "semimarkov":
        markov = SemiMarkov(str(get_shared_data("test_semi_markov.dat")))
        return markov.simulation_nb_elements(1, 1000, True)


class TestAddAbsorbingRun:
    """
    a main class to test the AddAbsorbingrun function on different type of
    data structure.

    """

    def test_max_length(self, AddAbsorbingRunData):
        seq = AddAbsorbingRunData.data
        assert seq.max_length == self.max_length

    def test_boost_versus_module(self, AddAbsorbingRunData):
        seq = AddAbsorbingRunData.data
        sequence_length = -1
        run_length = 6

        boost = seq.add_absorbing_run(sequence_length, run_length)
        module1 = AddAbsorbingRun(
            seq, SequenceLength=sequence_length, RunLength=run_length
        )
        module2 = AddAbsorbingRun(seq, RunLength=run_length)

        assert str(module1) == str(boost)
        assert str(module2) == str(boost)

    def test_no_arguments(self, AddAbsorbingRunData):
        seq = AddAbsorbingRunData.data
        assert AddAbsorbingRun(seq)

    def test_wrong_run_length(self, AddAbsorbingRunData):
        seq = AddAbsorbingRunData.data
        try:
            # second arguments must be less than MAX_RU_LENGTH
            _res = AddAbsorbingRun(seq, -1, AddAbsorbingRunData["MAX_RUN_LENGTH"] + 1)
            assert False
        except Exception:
            assert True

    def test_wrong_sequence_length(self, AddAbsorbingRunData):
        seq = AddAbsorbingRunData.data
        try:
            # second arguments must be less than MAX_RU_LENGTH
            _res = AddAbsorbingRun(seq, AddAbsorbingRunData["max_length"] - 1, -1)
            assert False
        except Exception:
            assert True


@pytest.fixture
def seq():
    return Sequences(str(get_shared_data("sequences1.seq")))


@pytest.fixture
def semi_markov():
    markov = SemiMarkov(str(get_shared_data("test_semi_markov.dat")))
    return markov.simulation_nb_elements(1, 1000, True)


def test_max_length_seq(seq):
    assert seq.max_length == 30


def test_boost_versus_module_seq(seq):
    sequence_length = -1
    run_length = 6

    boost = seq.add_absorbing_run(sequence_length, run_length)
    module1 = AddAbsorbingRun(seq, SequenceLength=sequence_length, RunLength=run_length)
    module2 = AddAbsorbingRun(seq, RunLength=run_length)

    assert str(module1) == str(boost)
    assert str(module2) == str(boost)


def test_no_arguments_seq(seq):
    assert AddAbsorbingRun(seq)


def test_wrong_run_length_seq(seq):
    try:
        # second arguments must be less than MAX_RU_LENGTH
        _res = AddAbsorbingRun(seq, -1, MAX_RUN_LENGTH + 1)
        assert False
    except Exception:
        assert True


def test_wrong_sequence_length_seq(seq):
    try:
        # second arguments must be less than MAX_RU_LENGTH
        max_length = 30
        _res = AddAbsorbingRun(seq, max_length - 1, -1)
        assert False
    except Exception:
        assert True
