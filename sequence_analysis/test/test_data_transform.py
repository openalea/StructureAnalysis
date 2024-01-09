"""various function tests

.. author:: Thomas Cokelaer, Thomas.Cokelaer@inria.fr

"""
__revision__ = "$Id: test_data_transform.py 9885 2010-11-06 18:19:34Z cokelaer $"


from openalea.sequence_analysis import get_shared_data
from openalea.sequence_analysis.data_transform import *
from openalea.sequence_analysis.sequences import Sequences
from openalea.stat_tool import Mixture
from openalea.stat_tool import Distribution
from openalea.stat_tool import Simulate
from openalea.stat_tool import Merge
from openalea.stat_tool import SelectStep
from openalea.stat_tool import Plot
from openalea.stat_tool import Display
from openalea.stat_tool import Vectors
from tools import runTestClass

class TestRemoveRun():

    def __init__(self):
        self.data = Sequences(get_shared_data("sequences1.seq"))

    def _test_sequences_1(self):

        seq1 = self.data
        seq2 = seq1.remove_run(1, 0,"e",2)
        seq3 = RemoveRun(seq1,1, 0,"e", MaxLength=2)
        assert str(seq3)==str(seq2)

    def test_incorrect_value(self):
        seq1 = self.data
        try:
            seq1.remove_run(1,3,'e',10)
            assert False
        except:
            assert True
        try:
            seq1.remove_run(1,-1,'e',10)
            assert False
        except:
            assert True

    def test_incorrect_variable(self):
        seq1 = self.data
        try:
            seq1.remove_run(0,2,'e',10)
            assert False
        except:
            assert True

    def _test_sequences_2(self):
        seq1 = Sequences(get_shared_data("sequences2.seq"))
        seq2 = seq1.remove_run(1, 0,"e",2)
        seq3 = RemoveRun(seq1,1, 0,"e",2)
        assert str(seq3)==str(seq2)


def test_markov_data():
    """not implemented"""
    pass

def test_semi_markov_data():
    """not implemented"""
    pass

def test_discrete_sequences():
    """not implemented"""
    pass

def test_compute_state_sequence():
    from openalea.sequence_analysis import HiddenSemiMarkov
    seq = Sequences(get_shared_data( "wij1.seq"))
    hsmc0 = HiddenSemiMarkov(get_shared_data(  "wij1.hsc"))
    ComputeStateSequences(seq, hsmc0, Algorithm="ForwardBackward", Characteristics=True)


def test_transition_count():
    seq = Sequences(get_shared_data( "wij1.seq"))
    TransitionCount(seq, 5, Begin=True, Estimator="MaximumLikelihood",
                    Filename = "ASCII")


def test_merge():

    mixt1 = Mixture(0.6, Distribution("B", 2, 18, 0.5),
                    0.4, Distribution("NB", 10, 10, 0.5))

    mixt_histo1 = Simulate(mixt1, 200)

    histo10 = mixt_histo1.extract_component(1)
    histo11 = mixt_histo1.extract_component(2)

    histo12 = Merge(histo10, histo11)

    assert histo12


def test_select_step():
    """
    #########################################################################
    #
    #  Well-log data; used in Fearnhead and Clifford "On-line Inference for 
    #  Hidden Markov Models via Particle Filters". Measurements of Nuclear-response
    #  of a well-bore over time. Data from O Ruanaidh, J. J. K. and 
    #  Fitzgerald, W. J. (1996). "Numerical Bayesion Methods Applied to Signal 
    #  Processing". New York: Springer.
    #
    #########################################################################
    """
    seq1 = Sequences(get_shared_data( "well_log_filtered.seq"))
    Plot(seq1, ViewPoint="Data")
    Plot(seq1)

    SelectStep(seq1, 1000)
    Plot(seq1)

    #Display(seq1, 1, 17, "Gaussian", ViewPoint="SegmentProfile", NbSegmentation=5)
    Plot(seq1, 1, 17, "Gaussian", ViewPoint="SegmentProfile")

    # seq20 = Segmentation(seq1, 1, 20, "Gaussian")
    # seq40 = Segmentation(seq1, 1, 40, "Gaussian")

    # seq20 = Segmentation(seq1, 1, 20, "Mean")
    # seq40 = Segmentation(seq1, 1, 40, "Mean")

    # seq16 = Segmentation(seq1, 1, 16, "Gaussian", NbSegment->"Fixed")

    vec1 = Vectors(seq1)
    Plot(vec1)

    SelectStep(vec1, 1000)
    Plot(vec1)

if __name__ == "__main__":
    runTestClass(TestMerge())
    runTestClass(TestRemoveRun())
    test_compute_state_sequence()
