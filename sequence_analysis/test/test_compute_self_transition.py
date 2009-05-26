""" ComputeSelfTransition tests

.. author:: Thomas Cokelaer, Thomas.Cokelaer@inria.fr
"""
__revision__ = "$Id: $"


from openalea.sequence_analysis.data_transform import ComputeSelfTransition
from data import seqn, seq1




def test_ComputeSelfTransition():
    ComputeSelfTransition(seq1)
    ComputeSelfTransition(seqn)

def test_ComputeSelfTransition_order():
    """not implemented see export_markovian_sequences code
    the order arguments is protected..."""
    ComputeSelfTransition(seqn, Order=2)

