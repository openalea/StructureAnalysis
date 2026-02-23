"""Test renewal data structure

.. author:: Thomas Cokelaer, Thomas.Cokelaer@inria.fr

"""

__revision__ = "$Id$"

import pytest

from openalea.sequence_analysis import _sequence_analysis as sa
from openalea.sequence_analysis.hidden_variable_order_markov import (
    HiddenVariableOrderMarkov,
)
from openalea.sequence_analysis.hidden_semi_markov import HiddenSemiMarkov
from openalea.sequence_analysis.renewal import Renewal
from .tools import runTestClass, robust_path as get_shared_data

N = 10


@pytest.fixture
def create_data_hidden_semi_markov():
    return HiddenSemiMarkov(str(get_shared_data("test_hidden_semi_markov.dat")))


# SEMI MARKOV case
def test_semi_markov_iterator(create_data_hidden_semi_markov):
    hsm = create_data_hidden_semi_markov
    smi = sa._SemiMarkovIterator(hsm)
    sim = smi.simulation(N, True)


def hsm_iterator(fn):
    hsm = HiddenSemiMarkov(fn)
    it = sa._SemiMarkovIterator(hsm)
    return it


def test_semi_markov_iterator2(create_data_hidden_semi_markov):
    fn = create_data_hidden_semi_markov
    smi = hsm_iterator(fn)
    sim = smi.simulation(N, True)


# VARIABLE ORDER MARKOV case


@pytest.fixture
def create_data_variable_order_markov():
    return str(get_shared_data("dupreziana21.hc"))


def vom_iterator(fn):
    vom = HiddenVariableOrderMarkov(fn)
    it = sa._VariableOrderMarkovIterator(vom)
    return it


def test_variable_order_markov_iterator(create_data_variable_order_markov):
    vom = HiddenVariableOrderMarkov(create_data_variable_order_markov)
    smi = sa._VariableOrderMarkovIterator(vom)
    sim = smi.simulation(N, True)


def test_variable_order_markov_iterator2(create_data_variable_order_markov):
    fn = create_data_variable_order_markov
    smi = vom_iterator(fn)
    sim = smi.simulation(N, True)


# RENEWAL case
def renewal_iterator(fn):
    ren = Renewal(fn)
    it = sa._RenewalIterator(ren)
    return it


def _test_renewal_iterator2():
    """to be fixed"""
    fn = path + "abri13.ren"
    smi = renewal_iterator(fn)
    sim = smi.simulation(N, True)


def _test_renewal_iterator():
    """to be fixed"""
    ren = Renewal(path + "abri13.ren")
    print((type(ren)))
    smi = sa._RenewalIterator(ren)
    sim = smi.simulation(N, True)
