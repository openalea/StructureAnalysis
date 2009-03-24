"""simulate tests"""
__revision__ = "$Id: $"

from openalea.stat_tool.mixture import Mixture
from openalea.stat_tool.convolution import Convolution
from openalea.stat_tool.simulate import Simulate


class Test:
    
    def test_simulate_mixture(self):
        
        m = Mixture("mixture1.mixt")
        s1 = Simulate(m, 1000)
        assert s1
        s1.plot()
        
    def test_simulate_convolution(self):
        
        c = Convolution("convolution1.conv")
        s1 = Simulate(c, 1000)
        assert s1
        s1.plot()
        

if __name__=="__main__":
    # perform all the test in the class Test (unit tests)
    test = Test()
    for method in dir(test):
        if method.startswith('_'):
            continue
        if callable(getattr(test, method)):
            getattr(test, method)()
        else:
            print 'skipping'
    # and functional tests.    
    
