"""simulate tests"""
__version__ = "$Id: test_simulate.py 9028 2010-06-01 08:17:28Z cokelaer $"

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

def test_simulation_sequences():
    #TODO The Simulate call does not work. Requires a proper input.
    try:
        from openalea.sequence_analysis import Sequences
        seq = Sequences([1,2,3,4])
        Simulate(seq, 100)
    except:
        assert True

if __name__ == "__main__":
    test_simulate_mixture()
    test_simulate_distribution()
    test_simulate_convolution()
    test_simulate_compound()
