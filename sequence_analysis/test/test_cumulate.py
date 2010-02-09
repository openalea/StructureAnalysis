""" Cumulate tests

.. author:: Thomas Cokelaer, Thomas.Cokelaer@inria.fr

.. todo : in general, variable index starts at 1 when calling Cumulate!
Do we want to start at 0 ? Since later on, python calls will 
start the index at 0 ?
 
"""
__revision__ = "$Id$"


from openalea.sequence_analysis.data_transform import Cumulate
from data import seqn
from data import seq1




def test_cumulate1():
    data = seq1
    a = Cumulate(data)
    b = data.cumulate(1).markovian_sequences()
   
    assert str(a) == str(b)
    assert a.get_max_value(0) == 29
    assert b.get_max_value(0) == 29
    
    assert data.get_max_value(0) == 2


def test_cumulaten():
    for var in range(1, seqn.nb_variable+1):
        assert str(seqn.cumulate(var).markovian_sequences()) == str(Cumulate(seqn, var))


if __name__ == "__main__":
    test_cumulate1()
    test_cumulaten()


