""" Test renewal data structure

.. author:: Thomas Cokelaer, Thomas.Cokelaer@inria.fr

"""
__revision__ = "$Id: test_iterator.py 9885 2010-11-06 18:19:34Z cokelaer $"

from openalea.sequence_analysis import _sequence_analysis as sa
from openalea.sequence_analysis.hidden_variable_order_markov import *
from openalea.sequence_analysis.hidden_semi_markov import *
from openalea.sequence_analysis.renewal import *
from openalea.sequence_analysis import get_shared_data

N = 10
import os
# SEMI MARKOV case
def test_semi_markov_iterator():
    hsm = HiddenSemiMarkov(get_shared_data('test_hidden_semi_markov.dat'))
    smi = sa._SemiMarkovIterator(hsm)
    sim = smi.simulation(N, True)

def hsm_iterator(fn):
    hsm = HiddenSemiMarkov(fn)
    it = sa._SemiMarkovIterator(hsm)
    return it

def test_semi_markov_iterator2():
    fn = get_shared_data('test_hidden_semi_markov.dat')
    smi = hsm_iterator(fn)
    sim = smi.simulation(N, True)

# VARIABLE ORDER MARKOV case
def vom_iterator(fn):

    vom = HiddenVariableOrderMarkov(fn)
    it = sa._VariableOrderMarkovIterator(vom)
    return it

def test_variable_order_markov_iterator():
    vom = HiddenVariableOrderMarkov(get_shared_data('dupreziana21.hc'))
    smi = sa._VariableOrderMarkovIterator(vom)
    sim = smi.simulation(N, True)

def test_variable_order_markov_iterator2():
    fn = get_shared_data('dupreziana21.hc')
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
    print type(ren)
    smi = sa._RenewalIterator(ren)
    sim = smi.simulation(N, True)



if __name__ == "__main__":

    test_semi_markov_iterator()
    test_semi_markov_iterator2()

    test_variable_order_markov_iterator()
    test_variable_order_markov_iterator2()

    #test_renewal_iterator()
    #test_renewal_iterator2()

