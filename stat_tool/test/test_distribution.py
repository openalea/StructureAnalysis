# -*- coding: utf-8 -*-
"""Distribution tests"""

__version__ = "$Id$"

from .tools import interface

from openalea.stat_tool._stat_tool import _DiscreteParametricModel

from openalea.stat_tool import Estimate
from openalea.stat_tool.distribution import (
    Binomial,
    Distribution,
    Multinomial,
    NegativeBinomial,
    Poisson,
    ToDistribution,
    ToHistogram,
    Uniform,
    set_seed,
)
from openalea.stat_tool.histogram import Histogram

import pytest

@pytest.fixture
def data():
    d1 = Binomial(0, 10, 0.5)
    d2 = Distribution("BINOMIAL", 0, 10, 0.5)
    assert d1 == d2
    return d1

@pytest.fixture
def myi(data):
    return interface(data, "data/distribution1.dist", Distribution)



def test_constructor_from_file(myi):
    myi.constructor_from_file()

def test_constructor_from_file_failure(myi):
    myi.constructor_from_file_failure()

def test_constructors():
    h = Histogram([1, 2, 3, 4, 5, 6, 1, 2, 3])
    assert h

    # from histogram
    dist = Distribution(h)
    assert dist

    # from parametric model
    pm = _DiscreteParametricModel(h)
    dist = Distribution(pm)
    assert dist

def test_constructor_failure():
    try:
        _h = Distribution("Whatever", 1, 1)
        assert False
    except KeyError:
        assert True

def test_print(myi):
    myi.print_data()

def test_display(myi):
    myi.display()
    myi.display_versus_ascii_write()
    myi.display_versus_str()

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
    myi.simulate()

def test_survival_ascii_write(myi):
    myi.survival_ascii_write()

def test_survival_spreadsheet_write(myi):
    myi.survival_spreadsheet_write()

def test_extract_data(myi):
    s = myi.simulate()
    e = s.estimate_parametric("B")
    d = e.extract_data()
    assert d
    _eprime = Estimate(s, "Binomial")

def test_truncate(data):
    # todo:find a test
    s = data
    _res = s.truncate(4)


class TestDistribution:
    """Test the distribution (Uniform, Binomial, ...)

    test the sup_bound, inf_bound, probability, parameter,ident

    test the ToDistribution and ToHistogram
    """
    def test_to_histogram(self):
        set_seed(0)
        d = Distribution("NEGATIVE_BINOMIAL", 0, 1, 0.5)
        h = d.simulate(1000)
        d2 = ToDistribution(h)
        assert h and d2

        h2 = ToHistogram(d2)
        assert h2
        assert h == h2

    def test_uniform(self):
        set_seed(0)
        d = Distribution("UNIFORM", 0, 10)
        s = d.simulate(1000)
        assert list(s)
        import numpy as np

        assert np.sum(s) == 1000
        assert d.get_sup_bound == 10
        assert d.get_inf_bound == 0
        assert d.get_probability == -1
        assert d.get_parameter == -1
        from openalea.stat_tool.enums import distribution_identifier_type as dist_type

        assert dist_type["UNIFORM"] == d.get_ident()

        d = Uniform(0, 10)
        s = d.simulate(1000)
        assert list(s)
        assert np.sum(s) == 1000

        m = d.simulate(1000).extract_model()
        assert isinstance(m, _DiscreteParametricModel)

    def test_binomial(self):
        set_seed(0)
        d = Distribution("BINOMIAL", 0, 10, 0.5)
        assert list(d.simulate(1000))
        assert d.get_sup_bound == 10
        assert d.get_inf_bound == 0
        assert d.get_probability == 0.5
        assert d.get_parameter == -1
        assert d.get_ident() == 1

        d = Binomial(0, 10, 0.5)
        assert list(d.simulate(1000))

        m = d.simulate(1000).extract_model()
        assert isinstance(m, _DiscreteParametricModel)

    def test_poisson(self):
        set_seed(0)
        d = Distribution("POISSON", 0, 2)
        assert list(d.simulate(1000))
        assert d.get_sup_bound == -1
        assert d.get_inf_bound == 0
        assert d.get_probability == -1
        assert d.get_parameter == 2
        assert d.get_ident() == 2

        d = Poisson(0, 2)
        assert list(d.simulate(1000))

        m = d.simulate(1000).extract_model()
        assert isinstance(m, _DiscreteParametricModel)

    def test_neg_binomial(self):
        set_seed(0)
        d = Distribution("NEGATIVE_BINOMIAL", 0, 1, 0.5)
        assert list(d.simulate(1000))
        assert d.get_sup_bound == -1
        assert d.get_inf_bound == 0
        assert d.get_probability == 0.5
        assert d.get_parameter == 1
        assert d.get_ident() == 3
        d = NegativeBinomial(0, 1, 0.5)
        assert list(d.simulate(1000))

        m = d.simulate(1000).extract_model()
        assert isinstance(m, _DiscreteParametricModel)

    def test_multimodial(self):
        """multinomial not yet implemented"""
        try:
            Distribution("MULTINOMIAL", 10)
            assert False
        except:
            assert True
        try:
            d = Multinomial()
            assert False
        except:
            assert True

    def test_simulation(self):
        set_seed(0)
        """simulate a vector of realizations"""
        d = Uniform(0, 1)
        try:
            import numpy

            l = numpy.zeros(100000)
            for i in range(len(l)):
                l[i] = d.simulation()
            assert (0 < sum(l)) & (sum(l) < len(l))
        except:
            pass

    def test_getters(self):
        dist = Distribution(Histogram([1, 1, 1, 2, 2, 2, 3, 3, 3]))
        assert dist.get_mean == 2
        assert 0.3333 < dist.get_max < 0.3334
        assert dist.get_mean == 2
        assert 0.6666 < dist.get_variance < 0.6667

    def test_set_seed(self):
        """ "Setting simulation seed"""

        d = Distribution("POISSON", 0, 2.5)
        set_seed(0)
        # Simulate and get histograms
        s1 = str(list(d.simulate(200)))
        s2 = str(list(d.simulate(200)))
        assert s1 != s2
        set_seed(0)
        s3 = str(list(d.simulate(200)))
        assert s1 == s3
        #return s1, s2, s3

