"""various function tests

.. author:: Thomas Cokelaer, Thomas.Cokelaer@inria.fr

"""
__revision__ = "$Id:  $"



from openalea.sequence_analysis.data_transform import *
from openalea.sequence_analysis.sequences import Sequences

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
