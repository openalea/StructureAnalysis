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
        
        s1.plot()
        
    def test_simulate_convolution(self):
        
        from convolution import Convolution
        c = Convolution(get_test_file("convolution1.conv"))
        
        s1 = Simulate(c, 1000)
        assert s1
        
        s1.plot()
        

