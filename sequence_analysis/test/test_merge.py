"""unitary or functional tests on Merge.
 See also test_semi_markov, test_time_events and so on
 
 .. author:: Thomas Cokelaer, Thomas.Cokelaer@inria.fr
"""
__revision__ = "$Id:  $"


from openalea.stat_tool.data_transform import Merge 
#from test_sequences import Test as sequences_data


def test_merge_histo():
    """test not yet implemented"""
    pass


def test_merge_time_events():
    """test not yet implemented"""
    pass


def test_merge_renewal_data():
    """test not yet implemented"""
    pass


def test_merge_sequences():
    from test_sequences import Test as sequences_data
    sequences = sequences_data()
    data = sequences.build_data()
    assert Merge(data, data)
    

def test_merge_vom_data():
    """test not yet implemented"""
    pass


def test_merge_semi_markov_data():
    """test not yet implemented"""
    pass


def test_merge_nonhomogenesous_markov_data():
    """test not yet implemented"""
    pass


def test_merge_tops():
    """test not yet implemented"""
    pass


def test_correlation():
    """test not yet implemented"""
    pass





