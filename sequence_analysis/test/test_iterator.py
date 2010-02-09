""" Test renewal data structure

.. author:: Thomas Cokelaer, Thomas.Cokelaer@inria.fr

"""
__revision__ = "$Id$"

from os import sep
from openalea.sequence_analysis import _sequence_analysis as sa
from openalea.sequence_analysis.hidden_variable_order_markov import *
from openalea.sequence_analysis.hidden_semi_markov import *
from openalea.sequence_analysis.renewal import *


N = 10

# SEMI MARKOV case
def test_semi_markov_iterator():
    hsm = HiddenSemiMarkov('data'+sep+'hidden_semi_markov.dat')
    smi = sa._Semi_markov_iterator(hsm)
    sim = smi.simulation(N, True)

def hsm_iterator(fn):
    hsm = HiddenSemiMarkov(fn)
    it = sa._Semi_markov_iterator(hsm)
    return it

def test_semi_markov_iterator2():
    fn = 'data'+sep+'hidden_semi_markov.dat'
    smi = hsm_iterator(fn)  
    sim = smi.simulation(N, True)

# VARIABLE ORDER MARKOV case
def vom_iterator(fn):
    
    vom = HiddenVariableOrderMarkov(fn)
    it = sa._Variable_order_markov_iterator(vom)
    return it

def test_variable_order_markov_iterator():
    vom = HiddenVariableOrderMarkov('data'+sep+'dupreziana21.hc')
    smi = sa._Variable_order_markov_iterator(vom)
    sim = smi.simulation(N, True)

def test_variable_order_markov_iterator2():
    fn = 'data'+sep+'dupreziana21.hc'
    smi = vom_iterator(fn)
    sim = smi.simulation(N, True)

# RENEWAL case
def renewal_iterator(fn):
    ren = Renewal(fn)
    it = sa._Renewal_iterator(ren)
    return it

def _test_renewal_iterator2():
    """to be fixed"""
    fn = "data" + sep + "abri13.ren"
    smi = renewal_iterator(fn)
    sim = smi.simulation(N, True)

def _test_renewal_iterator():
    """to be fixed"""
    ren = Renewal("data" + sep + "abri13.ren")
    print type(ren)
    smi = sa._Renewal_iterator(ren)
    sim = smi.simulation(N, True)



if __name__ == "__main__":
    
    test_semi_markov_iterator()
    test_semi_markov_iterator2()
    
    test_variable_order_markov_iterator()
    test_variable_order_markov_iterator2()
    
    #test_renewal_iterator()
    #test_renewal_iterator2()
    