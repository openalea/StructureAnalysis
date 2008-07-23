__doc__ = """ Simulation functions """
__docformat__ = "Restructuredtext"


def Simulate(obj, *args):
    """ Generation of a random saple from a distribution """

    return obj.simulate(*args)



################################################################################

class Test:
    
    def test_simulate_mixture(self):
        
        from mixture import Mixture
        m = Mixture("../../../test/mixture1.mixt")
        
        s1 = Simulate(m, 1000)
        assert s1
        

