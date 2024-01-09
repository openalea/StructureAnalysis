"""tests on the method AddAbsorbingRun

.. author:: Thomas Cokelaer, Thomas.Cokelaer@inria.fr 

.. todo:: markov case ? 
"""
__revision__ = "$Id: test_add_absorbing_run.py 9885 2010-11-06 18:19:34Z cokelaer $"

from openalea.sequence_analysis.sequences import Sequences
from openalea.sequence_analysis.semi_markov import SemiMarkov
from openalea.sequence_analysis.data_transform import AddAbsorbingRun
from openalea.sequence_analysis import get_shared_data

class _AddAbsorbingRun():
    """
    a main class to test the AddAbsorbingrun function on different type of 
    data structure.
    
    """

    def __init__(self):
        self.data = None 
        self.max_length = -1
        self.MAX_RUN_LENGTH = 20 # hardcoded values in CPP code
    
    def test_max_length(self):
        seq = self.data
        assert seq.max_length == self.max_length
        
        
    def test_boost_versus_module(self):
        seq = self.data
        sequence_length = -1
        run_length = 6
        
        boost = seq.add_absorbing_run(sequence_length, run_length)
        module1 = AddAbsorbingRun(seq, 
                                  SequenceLength=sequence_length,
                                  RunLength=run_length)
        module2 = AddAbsorbingRun(seq, RunLength=run_length)

        assert str(module1)==str(boost)
        assert str(module2)==str(boost)
        
    def test_no_arguments(self):
        seq = self.data
        assert AddAbsorbingRun(seq)
    
    def test_wrong_run_length(self):
        seq = self.data       
        try:
            #second arguments must be less than MAX_RU_LENGTH
            _res = AddAbsorbingRun(seq, -1, self.MAX_RUN_LENGTH + 1)
            assert False
        except Exception:
            assert True

    def test_wrong_sequence_length(self):
        seq = self.data       
        try:
            #second arguments must be less than MAX_RU_LENGTH
            _res = AddAbsorbingRun(seq, self.max_length -1, -1)
            assert False
        except Exception:
            assert True

class Test_AddAbsorbingRun_Sequences(_AddAbsorbingRun):
    """sequences case"""
    def __init__(self):
        _AddAbsorbingRun.__init__(self)
        self.data = self.create_data()
        self.max_length = 30

    def create_data(self):
        seq = Sequences(get_shared_data('sequences1.seq'))
        return seq
    
class Test_AddAbsorbingRun_SemiMarkov(_AddAbsorbingRun):
    """semi markov case"""
    def __init__(self):
        _AddAbsorbingRun.__init__(self)
        self.data = self.create_data()
        self.max_length = 1000

    def create_data(self):
        markov = SemiMarkov(get_shared_data('test_semi_markov.dat') )
        semi_markov_data = markov.simulation_nb_elements(1, 1000, True)
        return semi_markov_data
    
