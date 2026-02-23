"""histogram tests"""


try:
    from .tools import interface
    from .tools import robust_path as get_shared_data
except ImportError:
    from tools import interface
    from tools import robust_path as get_shared_data

import openalea.stat_tool as stat_tool
from openalea.stat_tool.distribution import (
    Binomial,
    Distribution,
    ToDistribution,
    ToHistogram,
)
from openalea.stat_tool.histogram import Histogram

from openalea.stat_tool.distribution import set_seed

import pytest

@pytest.fixture
def data():
    v = Histogram([0, 1, 2, 3])
    return v

@pytest.fixture
def myi(data):
    return interface(data, "data/peup1.his", Histogram)


def test_constructor_from_list():
    v = Histogram([0, 1, 2, 3])
    assert v

def test_constructor_from_integers():
    v = Histogram(0, 1, 2, 3)
    assert v

def test_constructor_from_file(myi):
    myi.constructor_from_file()

def test_constructor_from_file_failure(myi):
    myi.constructor_from_file_failure()

def test_print(myi):
    myi.print_data()

def test_display(myi):
    myi.display()
    myi.display_versus_ascii_write()
    myi.display_versus_str()

def test_len(data):
    v = data
    assert len(v) == 4

def test_plot(myi):
    myi.plot()

def test_plot_null_variance():
    """Test plotting frequency distribution
    with null variance"""
    v = Histogram([0] * 100)
    v.plot()

def test_save(myi):
    myi.save()

def test_plot_write(myi):
    myi.plot_write()

def test_file_ascii_write(myi):
    myi.file_ascii_write()

def test_spreadsheet_write(myi):
    myi.spreadsheet_write()

def test_survival_plot_write(myi):
    myi.survival_plot_write()

def test_survival_spreadsheet_write(myi):
    myi.survival_spreadsheet_write()

def test_extract_data():
    """test extract_data"""
    h = Histogram(str(get_shared_data("meri1.his")))

    e = h.estimate_nonparametric()

    assert e.extract_data()

def test_survival_ascii_write(data):
    """test display"""
    h = data
    h.survival_ascii_write()

def test_container():
    """container / iterator"""
    h = Histogram(str(get_shared_data("meri1.his")))

    assert h[0] == 0
    assert h[10] == 1

def test_to_histogram():
    """Test the ToHistogram function"""

    d = Distribution("NEGATIVE_BINOMIAL", 0, 1, 0.5)
    set_seed(0)
    h = d.simulate(1000)
    d2 = ToDistribution(h)
    assert h and d2

    h2 = ToHistogram(d2)
    assert h2
    assert h == h2

def test_extract_model():
    set_seed(0)
    d = Binomial(0, 10, 0.5)
    d == d.simulate(1000).extract_model()


def test_extract_vec():
    """Extract histogram from Vectors"""

    vs = stat_tool.Vectors("data/cvectors.vec")
    h = vs.extract(4)
    assert h
    print(h)

if __name__ == "__main__":

    def data():
        v = Histogram([0, 1, 2, 3])
        return v            
    
    myi = interface(data, "data/peup1.his", Histogram)
    myi.data = data()

    test_save(myi)
