""" ComputeSelfTransition tests

.. author:: Thomas Cokelaer, Thomas.Cokelaer@inria.fr
"""
__revision__ = "$Id$"


from openalea.sequence_analysis.data_transform import ComputeSelfTransition
from openalea.sequence_analysis import Sequences, get_shared_data

seq1 = Sequences(get_shared_data("sequences1.seq"))
seqn = Sequences(get_shared_data("sequences2.seq"))





def test_ComputeSelfTransition():
    ComputeSelfTransition(seq1)
    ComputeSelfTransition(seqn)

def test_ComputeSelfTransition_order():
    """not implemented see export_markovian_sequences code
    the order arguments is protected..."""
    ComputeSelfTransition(seqn, Order=2)

