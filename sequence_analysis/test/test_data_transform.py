"""various function tests

.. author:: Thomas Cokelaer, Thomas.Cokelaer@inria.fr

"""
__revision__ = "$Id:  $"



from openalea.sequence_analysis.data_transform import *
from openalea.sequence_analysis.sequences import Sequences
from openalea.stat_tool import Mixture
from openalea.stat_tool import Distribution
from openalea.stat_tool import Simulate
from openalea.stat_tool import Merge

class TestRemoveRun():

    def __init__(self):
        self.data = Sequences("data/sequences1.seq")
        
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
        seq1 = Sequences("data/sequences2.seq")
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



class TestMerge():
    
    def __init__(self):
        pass
        #data.__init__(self)
        
    def test_merge(self):
        
        mixt1 = Mixture(0.6, Distribution("B", 2, 18, 0.5), 
                        0.4, Distribution("NB", 10, 10, 0.5))

        mixt_histo1 = Simulate(mixt1, 200)

        histo10 = mixt_histo1.extract_component(1)
        histo11 = mixt_histo1.extract_component(2)

        histo12 = Merge(histo10, histo11)

        assert histo12
