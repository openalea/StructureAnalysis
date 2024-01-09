""" Difference tests

.. author:: Thomas Cokelaer, Thomas.Cokelaer@inria.fr

 
"""
__revision__ = "$Id: test_difference.py 9885 2010-11-06 18:19:34Z cokelaer $"


from openalea.sequence_analysis.data_transform import Difference
from openalea.sequence_analysis import Sequences, get_shared_data

seq1 = Sequences(get_shared_data("sequences1.seq"))
seqn = Sequences(get_shared_data("sequences2.seq"))


def test_difference1():
    """difference test to finalise"""
    data = seq1
    res = Difference(data,1)
    assert str(res)==str(data.difference(1, False))
    assert res.cumul_length == 50

def test_difference1_first_element():
    """difference test to finalise"""
    data = seq1
    res = Difference(data,1, True)
    assert str(res)==str(data.difference(1, True))
    assert res.cumul_length == 52

def test_differencen():
    """difference test to finalise"""
    data = seqn
    res = Difference(data,Variable=1)
    assert str(res)==str(data.difference(1, False))
    assert res.cumul_length == 23

