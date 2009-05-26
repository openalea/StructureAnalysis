"""tests on the method estimate

.. author:: Thomas Cokelaer, Thomas.Cokelaer@inria.fr

.. todo:: finalise
"""
__revision__ = "$Id: $"

from openalea.sequence_analysis.sequences import Sequences
from openalea.stat_tool.vectors import Vectors
from openalea.sequence_analysis.hidden_variable_order_markov import HiddenVariableOrderMarkov
from openalea.sequence_analysis.estimate import Estimate
from openalea.stat_tool.data_transform import ExtractHistogram

from data import *

   
class Test_Estimate_Histogram():
    def __init__(self):
        self.data = self.create_data()

    def create_data(self):
        seq0 = Sequences("data/chene_sessile_15pa.seq")
        vec10 = Vectors(seq0)
        return vec10
    
    def test_estimate_mixture(self):
        mixt20 = Estimate(ExtractHistogram(self.data, 2), 
                          "MIXTURE", "NB", "NB", "NB", "NB", 
                          NbComponent="Estimated")
        assert mixt20.nb_component() == 2

    def test_estimate_mixture2(self):
        mixt20 = Estimate(ExtractHistogram(self.data, 5), 
                          "MIXTURE", "NB", "NB", "NB", "NB", 
                          NbComponent="Estimated")
        assert mixt20.nb_component() == 3

class Test_Estimate_VARIABLE_ORDER_MARKOV():
    
    def __init__(self):
        self.sequence = seq10
        self.type = "VARIABLE_ORDER_MARKOV"
        
    def test_estimate(self):
        mc10 = Estimate(self.sequence, self.type, "Ordinary",
                         MaxOrder=5, GlobalInitialTransition=True)

    def test_estimate1(self):        
        mc11 = Estimate(self.sequence , self.type, "Ordinary",
                        MaxOrder=5, GlobalInitialTransition=False)

    def test_estimate2(self):
        mc12 = Estimate(self.sequence, self.type, "Ordinary",
                         Algorithm="LocalBIC", Threshold=10.,
                         MaxOrder=5, GlobalInitialTransition=False,
                         GlobalSample=False)
    
    def test_estimate3(self):
        mc13 = Estimate(self.sequence, self.type, "Ordinary",
                        Algorithm="Context", Threshold=1.,
                        MaxOrder=5, GlobalInitialTransition=False, 
                        GlobalSample=False)


class Test_Estimate_VARIABLE_ORDER_MARKOV_from_markovian():
    def test_estimate(self):
        mc11 = Estimate(seq10 , "VARIABLE_ORDER_MARKOV", "Ordinary",
                        MaxOrder=5, GlobalInitialTransition=False)
        mc2 = Estimate(seq2, "VARIABLE_ORDER_MARKOV", 
                       mc11, GlobalInitialTransition=False)

class Test_Estimate_HIDDEN_VARIABLE_ORDER_MARKOV():
    def test_estimate(self):    
        from data import hvom_sample, seq1
        hmc_estimated = Estimate(seq1, "HIDDEN_VARIABLE_ORDER_MARKOV", hvom_sample,
                         GlobalInitialTransition=True, NbIteration=80)
        assert hmc_estimated
        
class Test_Estimate_time_events():
    """test not yet implemented"""
    pass

class Test_Estimate_tops():
    """tests not yet implemented"""
    pass
    
