"""Convolution tests

:Author: Thomas Cokelaer, Thomas.Cokelaer@inria.fr

"""

__version__ = "$Id$"

from .tools import interface, runTestClass

from openalea.stat_tool import Convolution, Estimate
from openalea.stat_tool.data_transform import ExtractDistribution
from openalea.stat_tool.distribution import Binomial, NegativeBinomial

from openalea.stat_tool.distribution import set_seed

import pytest

@pytest.fixture
def data():
    d1 = Binomial(0, 10, 0.5)
    d2 = NegativeBinomial(0, 1, 0.1)
    conv = Convolution(d1, d2)
    return conv

@pytest.fixture
def myi(data):
    return interface(data, "data/convolution1.conv", Convolution)

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

def test_len(myi):
    c = myi.data
    assert len(c) == 2
    assert len(c) == c.nb_distribution()

def test_plot(myi):
    d = myi.data
    d.plot()

def test_save(myi):
    myi.save()

def test_plot_write(myi):
    myi.plot_write()

def test_file_ascii_write(myi):
    myi.file_ascii_write()

def test_spreadsheet_write(myi):
    myi.spreadsheet_write()

def test_simulate(myi):
    myi.simulate()

def test_extract(data):
    """run and test the extract methods"""

    m = data
    assert m.extract_convolution() == ExtractDistribution(m, "Convolution")
    assert m.extract(1) == Binomial(0, 10, 0.5)
    assert m.extract(2) == NegativeBinomial(0, 1, 0.1)
    assert ExtractDistribution(m, "Elementary", 1) == Binomial(0, 10, 0.5)
    assert ExtractDistribution(m, "Elementary", 2) == NegativeBinomial(0, 1, 0.1)

def test_extract_data(data):
    """run and test the extract_data methods"""
    set_seed(0)
    m = data
    s = m.simulate(1000)
    e = s.estimate_convolution(Binomial(0, 10, 0.5), Estimator="Parametric")
    d = e.extract_data()
    assert d
    eprime = Estimate(s, "CONVOLUTION", Binomial(0, 10, 0.5))

    # todo: find robust assert ?
    print(eprime.get_mean)
    assert eprime

def test_truncate(data):
    s = data
    _res = s.truncate(4)

