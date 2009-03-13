__doc__ = """ Simulation functions """
#__docformat__ = "Restructuredtext"


def Simulate(obj, *args):
    """ Generation of a random sample from a distribution """

    return obj.simulate(*args)



################################################################################
from openalea.stat_tool import get_test_file

class Test:
    
    def test_simulate_mixture(self):
        
        from mixture import Mixture
        m = Mixture(get_test_file("mixture1.mixt"))
        
        s1 = Simulate(m, 1000)
        assert s1
        

