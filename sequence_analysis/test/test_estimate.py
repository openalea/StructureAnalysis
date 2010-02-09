"""tests on the method estimate

.. author:: Thomas Cokelaer, Thomas.Cokelaer@inria.fr

.. todo:: finalise
"""
__revision__ = "$Id$"

from openalea.sequence_analysis.sequences import Sequences
from openalea.stat_tool.vectors import Vectors
from openalea.sequence_analysis.hidden_variable_order_markov import \
    HiddenVariableOrderMarkov
from openalea.sequence_analysis.hidden_semi_markov import \
    HiddenSemiMarkov
    
from openalea.sequence_analysis.estimate import Estimate
from openalea.stat_tool.data_transform import ExtractHistogram

from data import *
from tools import runTestClass
from test_tops import TopsData
from test_hidden_semi_markov import HiddenSemiMarkovData
from test_semi_markov import SemiMarkovData

   
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
        assert mixt20.nb_component == 2

    def test_estimate_mixture2(self):
        mixt20 = Estimate(ExtractHistogram(self.data, 5), 
                          "MIXTURE", "NB", "NB", "NB", "NB", 
                          NbComponent="Estimated")
        assert mixt20.nb_component == 3

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
    def test_estimate4(self):
        for Algorithm in ["CTM_BIC", "CTM_KT", "Context"]:
            mc13 = Estimate(self.sequence, self.type, "Ordinary",
                        Algorithm=Algorithm,
                        MaxOrder=5, GlobalInitialTransition=False, 
                        GlobalSample=False)
    def test_estimate_error1(self):
        """test that Estimator and Algorith=CTM_KT are incompatible"""
        try:
            mc13 = Estimate(self.sequence, self.type, "Ordinary",
                        Algorithm="CTM_KT", Estimator="Laplace",
                        MaxOrder=5, GlobalInitialTransition=False, 
                        GlobalSample=False)
            assert False
        except:
            assert True
    
            
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
        
class Test_Estimate_HIDDEN_SEMI_MARKOV():
    
    def __init__(self):
        self.data = HiddenSemiMarkovData()
        self.sequence = Sequences("data/wij1.seq")

        
    def test_estimate(self):
        seq = self.sequence
        # data is a hsm class
        Estimate(seq, "HIDDEN_SEMI-MARKOV", self.data)


class Test_Estimate_SEMI_MARKOV():
    
    def __init__(self):
        #self.data = SemiMarkovData()
        self.sequence = Sequences("data/wij1.seq")
        
    def _test_estimate(self):
        seq = self.sequence
        # data is a hsm class
        Estimate(seq, "SEMI-MARKOV", "Ordinary")
        
        
class Test_Estimate_time_events():
    """test not yet implemented"""
    pass

class Test_Estimate_tops():
    """tests not yet implemented"""
    
    def __init__(self):
        self.data = TopsData()
        
    def test_estimate(self):
        Estimate(self.data, MinPosition=1, MaxPosition=10)
    


if __name__ == "__main__":
    runTestClass(Test_Estimate_HIDDEN_VARIABLE_ORDER_MARKOV())
    runTestClass(Test_Estimate_VARIABLE_ORDER_MARKOV_from_markovian())
    runTestClass(Test_Estimate_VARIABLE_ORDER_MARKOV())
    runTestClass(Test_Estimate_Histogram())
    runTestClass(Test_Estimate_tops())
    runTestClass(Test_Estimate_HIDDEN_SEMI_MARKOV())
    runTestClass(Test_Estimate_SEMI_MARKOV())
