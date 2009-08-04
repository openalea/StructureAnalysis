"""simulate tests"""
__version__ = "$Id$"

from openalea.stat_tool.mixture import Mixture
from openalea.stat_tool.convolution import Convolution
from openalea.stat_tool.compound import Compound
from openalea.stat_tool.simulate import Simulate
from openalea.stat_tool.distribution import Distribution


def test_simulate_mixture():
    m = Mixture("data/mixture1.mixt")
    s1 = Simulate(m, 1000)
    assert s1

    
def test_simulate_convolution():
    c = Convolution("data/convolution1.conv")
    s1 = Simulate(c, 1000)
    assert s1

    
def test_simulate_compound():
    c = Compound("data/compound1.cd")
    s1 = Simulate(c, 1000)
    assert s1

    
def test_simulate_distribution():
    c = Distribution("data/distribution1.dist")
    s1 = Simulate(c, 1000)
    assert s1


if __name__ == "__main__":
    test_simulate_mixture()
    test_simulate_distribution()
    test_simulate_convolution()
    test_simulate_compound()
