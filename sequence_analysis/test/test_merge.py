"""unitary or functional tests on Merge.
 See also test_semi_markov, test_time_events and so on
 
 .. author:: Thomas Cokelaer, Thomas.Cokelaer@inria.fr
"""
__revision__ = "$Id: test_merge.py 7884 2010-02-09 07:42:43Z cokelaer $"


from openalea.stat_tool.data_transform import Merge 

from test_correlation import CorrelationData
from test_tops import TopsData
from test_semi_markov import SemiMarkovData


def test_merge_histo():
    """see stat_tool/test/test_data_transform"""
    pass


def test_merge_time_events():
    """see test_time_events_functional"""
    pass


def test_merge_renewal_data():
    """See test_renewal_functional"""
    pass


def test_merge_sequences():
    from test_sequences import Test as sequences_data
    sequences = sequences_data()
    data = sequences.build_data()
    assert Merge(data, data)
    

def test_merge_vom_data():
    """test not yet implemented"""
    pass


def _test_merge_semi_markov_data():
    sm1 = SemiMarkovData()
    sm2 = SemiMarkovData()
    sm = Merge(sm1, sm2)
    
    #todo this plot does not work right now ?
    #sm.plot()


def test_merge_nonhomogenesous_markov_data():
    """test not yet implemented"""
    pass


def test_merge_tops():
    t1 = TopsData()
    t2 = TopsData()
    t = Merge(t1, t2)
    t.plot()
    
def test_merge_correlation():
    c1 = CorrelationData(1)
    c2 = CorrelationData(2)
    c3 = CorrelationData(3)
    c = Merge(c1,c2,c3)
    c_bis = c1.merge([c2,c3])
    assert str(c)==str(c_bis)
    c.plot()


if __name__ == "__main__":
    test_merge_tops()
    test_merge_correlation()
    #test_merge_semi_markov_data()
    test_merge_sequences()

