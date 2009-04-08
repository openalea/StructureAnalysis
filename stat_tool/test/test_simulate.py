"""simulate tests"""
__revision__ = "$Id$"

from openalea.stat_tool.mixture import Mixture
from openalea.stat_tool.convolution import Convolution
from openalea.stat_tool.compound import Compound
from openalea.stat_tool.simulate import Simulate
from openalea.stat_tool.distribution import Distribution


class Test:
    
    def __init__(self):
        pass
    
    def test_simulate_mixture(self):
        m = Mixture("mixture1.mixt")
        s1 = Simulate(m, 1000)
        assert s1
    
    def test_simulate_convolution(self):
        c = Convolution("convolution1.conv")
        s1 = Simulate(c, 1000)
        assert s1
    
    def test_simulate_compound(self):
        c = Compound("compound1.cd")
        s1 = Simulate(c, 1000)
        assert s1
    
    def test_simulate_distribution(self):
        c = Distribution("distribution1.dist")
        s1 = Simulate(c, 1000)
        assert s1

