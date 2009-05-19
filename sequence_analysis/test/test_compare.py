"""tests on the method AddAbsorbingRun

.. author:: Thomas Cokelaer, Thomas.Cokelaer@inria.fr 
"""
__revision__ = "$Id: $"

from openalea.sequence_analysis.sequences import Sequences
from openalea.sequence_analysis.semi_markov import SemiMarkov
from openalea.sequence_analysis.compare import Compare
from openalea.stat_tool.vectors import VectorDistance


class _Compare():
    """
    a main class to test the AddAbsorbingrun function on different type of 
    data structure.
    """

    def __init__(self):
        self.data = None 
        
    

class Test_AddAbsorbingRun_Sequences(_Compare):
    """sequences case"""
    def __init__(self):
        _Compare.__init__(self)
        self.data = self.create_data()
        self.max_length = 30

    def create_data(self):
        seq = Sequences('data/dupreziana_a1.seq')
        return seq

    def test_compare_vector_distance(self):
        seq = self.data
        matrix20 = Compare(seq, VectorDistance("N", "N"))
