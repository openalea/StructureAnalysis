""" Cluster tests

.. todo:: check the AddVariable option (sequences) and sequences cases
"""
__revision__ = "$Id: test_cluster.py 6325 2009-04-29 16:20:55Z cokelaer $"


from openalea.sequence_analysis.data_transform import ComputeSelfTransition
from data import seqn, seq1




def test_ComputeSelfTransition():
    ComputeSelfTransition(seq1)
    ComputeSelfTransition(seqn)

def test_ComputeSelfTransition_order():
    """not implemented"""
    ComputeSelfTransition(seqn, Order=2)

