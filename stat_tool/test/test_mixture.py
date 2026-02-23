"""mixture tests"""

__version__ = "$Id$"

try:
    from .tools import interface
    from .tools import robust_path as get_shared_data
except ImportError:
    from tools import interface
    from tools import robust_path as get_shared_data

from openalea.stat_tool.data_transform import ExtractDistribution
from openalea.stat_tool.distribution import Binomial
from openalea.stat_tool.estimate import Estimate
from openalea.stat_tool.histogram import Histogram
from openalea.stat_tool.mixture import Mixture

import pytest

@pytest.fixture
def data():
    d1 = Binomial(0, 12, 0.1)
    d2 = Binomial(0, 12, 0.5)
    d3 = Binomial(0, 12, 0.8)

    mixt = Mixture(0.1, d1, 0.2, d2, 0.7, d3)
    assert mixt.nb_component == 3
    return mixt

@pytest.fixture
def myi(data):
    return interface(data, "data/mixture1.mixt", Mixture)



def test_constructor_from_file(myi):
    myi.constructor_from_file()

def test_constructor_from_file_failure(myi):
    myi.constructor_from_file_failure()

def test_constructor_from_dists_failure():
    d1 = Binomial(0, 12, 0.1)
    d2 = Binomial(0, 12, 0.5)
    try:
        _mixt = Mixture(0.1, d1, d2)
        assert False
    except TypeError:
        assert True

def test_print(myi):
    myi.print_data()

def test_display(myi):
    myi.display()
    myi.display_versus_ascii_write()
    myi.display_versus_str()

def test_len(data):
    c = data
    assert len(c) == 3

def test_plot(myi):
    myi.plot()

def test_save(myi):
    myi.save()

def test_plot_write(myi):
    myi.plot_write()

def test_file_ascii_write(myi):
    myi.file_ascii_write()

def test_spreadsheet_write(myi):
    myi.spreadsheet_write()

def test_simulate(myi):
    sim = myi.simulate()
    sim.plot()

def test_estimate(myi):
    sim = myi.simulate()
    # 3 Binomial distribution to match th original data
    est = Estimate(sim, "Mixture", "B", "B", "B")
    est.plot()

def test_extract(data):
    """run and test the extract methods"""

    m = data

    assert m.extract_weight() == ExtractDistribution(m, "Weight")
    assert m.extract_mixture() == ExtractDistribution(m, "Mixture")

    assert ExtractDistribution(m, "Component", 1) == Binomial(0, 12, 0.1)
    assert ExtractDistribution(m, "Component", 2) == Binomial(0, 12, 0.5)
    assert ExtractDistribution(m, "Component", 3) == Binomial(0, 12, 0.8)

    assert m.extract_component(1) == Binomial(0, 12, 0.1)
    assert m.extract_component(2) == Binomial(0, 12, 0.5)
    assert m.extract_component(3) == Binomial(0, 12, 0.8)

def test_extract_data():
    """run and test the extract_data methods"""

    h = Histogram(str(get_shared_data("meri2.his")))
    m = h.estimate_DiscreteMixture("B", "NB")

    d = m.extract_data()
    assert d

def test_truncate(data):
    s = data
    _res = s.truncate(4)

