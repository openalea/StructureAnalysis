""" Test renewal data structure

.. author:: Thomas Cokelaer, Thomas.Cokelaer@inria.fr

"""
__revision__ = "$Id:  $"

from os import sep
from openalea.sequence_analysis import _sequence_analysis as sa
from openalea.sequence_analysis.hidden_variable_order_markov import *
from openalea.sequence_analysis.hidden_semi_markov import *


def test_semi_markov_iterator():
    hsm = HiddenSemiMarkov('data'+sep+'hidden_semi_markov.dat')
    smi = sa._Semi_markov_iterator(hsm)
    sim = smi.simulation(10, True)

def hsm_iterator(fn):
    hsm = HiddenSemiMarkov(fn)
    it = sa._Semi_markov_iterator(hsm)
    return it

def test_semi_markov_iterator2():
    fn = 'data/hidden_semi_markov.dat'
    smi = hsm_iterator(fn)
    sim = smi.simulation(10, True)

def test_variable_order_markov_iterator():
    vom = HiddenVariableOrderMarkov('data'+sep+'dupreziana21.hc')
    smi = sa._Variable_order_markov_iterator(vom)
    sim = smi.simulation(10, True)

def test_renewal_iterator():
    """to be done. needs renewal data """
#    vom = HiddenVariableOrderMarkov('data'+sep+'dupreziana21.hc')
#    smi = sa._Variable_order_markov_iterator(vom)
#    sim = smi.simulation(10, True)
