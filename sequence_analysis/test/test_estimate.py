"""tests on the method estimate

.. author:: Thomas Cokelaer, Thomas.Cokelaer@inria.fr

.. todo:: finalise
"""
__revision__ = "$Id: test_estimate.py 9885 2010-11-06 18:19:34Z cokelaer $"

from openalea.stat_tool.vectors import Vectors
from openalea.stat_tool.data_transform import ExtractHistogram

from tools import runTestClass
from test_tops import TopsData
from test_hidden_semi_markov import HiddenSemiMarkovData
from test_semi_markov import SemiMarkovData

from openalea.sequence_analysis import *

_seq1 = Sequences(get_shared_data('dupreziana_20a2.seq'))

seq2 = RemoveRun(_seq1, 1, 0, "End")
seq3 = Sequences(get_shared_data('dupreziana_40a2.seq'))
seq4_0 = RemoveRun(seq3, 2, 0, "End")
seq4 = SegmentationExtract(seq4_0, 1, 2)
seq5 = Sequences(get_shared_data('dupreziana_60a2.seq'))
seq6_0 = RemoveRun(seq5, 2, 0, "End")
seq6 = LengthSelect(SegmentationExtract(seq6_0, 1, 2), 1, Mode="Reject")
seq7 = Sequences(get_shared_data('dupreziana_80a2.seq'))
seq8_0 = RemoveRun(seq7, 2, 0, "End")
seq8 = SegmentationExtract(seq8_0, 1, 2)
seq10 = Merge(seq2, seq4, seq6, seq8)




class Test_Estimate_Histogram():
    def __init__(self):
        self.data = self.create_data()

    def create_data(self):
        seq0 = Sequences(get_shared_data( "chene_sessile_15pa.seq"))
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
        seq1 = Sequences(get_shared_data('sequences1.seq'))
        hvom_sample = HiddenVariableOrderMarkov(get_shared_data("dupreziana21.hc"))
        hmc_estimated = Estimate(seq1, "HIDDEN_VARIABLE_ORDER_MARKOV", hvom_sample,
                         GlobalInitialTransition=True, NbIteration=80)
        assert hmc_estimated

class Test_Estimate_HIDDEN_SEMI_MARKOV():

    def __init__(self):
        self.data = HiddenSemiMarkovData()
        self.sequence = Sequences(get_shared_data( "wij1.seq"))


    def test_estimate(self):
        seq = self.sequence
        # data is a hsm class
        Estimate(seq, "HIDDEN_SEMI-MARKOV", self.data)


class Test_Estimate_SEMI_MARKOV():

    def __init__(self):
        #self.data = SemiMarkovData()
        self.sequence = Sequences(get_shared_data( "wij1.seq"))

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
