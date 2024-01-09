# -*- coding: utf-8 -*-
"""Distribution tests"""
__version__ = "$Id: test_distribution.py 9869 2010-11-04 17:31:41Z dbarbeau $"


from openalea.stat_tool import _stat_tool
from openalea.stat_tool.distribution import Distribution, Uniform, Binomial
from openalea.stat_tool.distribution import NegativeBinomial, Poisson
from openalea.stat_tool.distribution import Multinomial
from openalea.stat_tool.distribution import ToHistogram, ToDistribution
from openalea.stat_tool.histogram import Histogram
from openalea.stat_tool import Estimate



from openalea.stat_tool._stat_tool import _DiscreteParametricModel

from tools import interface
from tools import runTestClass


class Test(interface):
    """a simple unittest class


    """
    def __init__(self):
        interface.__init__(self,
                           self.build_data(),
                           "data/distribution1.dist",
                           Distribution)

    def build_data(self):
        d1 = Binomial(0, 10, 0.5)
        d2 = Distribution("BINOMIAL", 0, 10, 0.5)
        assert d1 == d2
        return d1

    def test_constructor_from_file(self):
        self.constructor_from_file()

    def test_constructor_from_file_failure(self):
        self.constructor_from_file_failure()

    def test_constructors(self):
        h = Histogram([1, 2, 3, 4, 5, 6, 1, 2, 3])
        assert h

        # from histogram
        dist = Distribution(h)
        assert dist

        #from parametric model
        pm = _DiscreteParametricModel(h)
        dist = Distribution(pm)
        assert dist

    def test_constructor_failure(self):
        try:
            _h = Distribution("Whatever", 1, 1)
            assert False
        except KeyError:
            assert True

    def test_print(self):
        self.print_data()

    def test_display(self):
        self.display()
        self.display_versus_ascii_write()
        self.display_versus_str()

    def test_len(self):
        """not implemented; irrelevant?"""
        pass

    def test_plot(self):
        self.plot()

    def test_save(self):
        self.save()

    def test_plot_write(self):
        self.plot_write()

    def test_file_ascii_write(self):
        self.file_ascii_write()

    def test_spreadsheet_write(self):
        self.spreadsheet_write()

    def test_simulate(self):
        self.simulate()

    def test_survival_ascii_write(self):
        self.survival_ascii_write()

    def test_survival_spreadsheet_write(self):
        self.survival_spreadsheet_write()

    def test_extract(self):
        pass

    def test_extract_data(self):
        s = self.simulate()
        e = s.estimate_parametric("B")
        d = e.extract_data()
        assert d
        _eprime = Estimate(s, "Binomial")

    def test_truncate(self):
        # todo:find a test
        s = self.data
        _res = s.truncate(4)


class TestDistribution():
    """Test the distribution (Uniform, Binomial, ...)

    test the sup_bound, inf_bound, probability, parameter,ident

    test the ToDistribution and ToHistogram
    """
    def __init__(self):
        pass

    def test_to_histogram(self):

        d = Distribution("NEGATIVE_BINOMIAL", 0, 1, 0.5)
        h = d.simulate(1000)
        d2 = ToDistribution(h)
        assert h and d2

        h2 = ToHistogram(d2)
        assert h2
        assert h == h2

    def test_uniform(self):

        d = Distribution("UNIFORM", 0, 10)
        assert list(d.simulate(1000))
        assert d.get_sup_bound == 10
        assert d.get_inf_bound == 0
        assert d.get_probability == -1
        assert d.get_parameter == -1
        assert d.get_ident == 4

        d = Uniform(0, 10)
        assert list(d.simulate(1000))

        m = d.simulate(1000).extract_model()
        assert isinstance(m, _DiscreteParametricModel)

    def test_binomial(self):

        d = Distribution("BINOMIAL", 0, 10, 0.5)
        assert list(d.simulate(1000))
        assert d.get_sup_bound == 10
        assert d.get_inf_bound == 0
        assert d.get_probability == 0.5
        assert d.get_parameter == -1
        assert d.get_ident == 1

        d = Binomial(0, 10, 0.5)
        assert list(d.simulate(1000))

        m = d.simulate(1000).extract_model()
        assert isinstance(m, _DiscreteParametricModel)

    def test_poisson(self):

        d = Distribution("POISSON", 0, 2)
        assert list(d.simulate(1000))
        assert d.get_sup_bound == -1
        assert d.get_inf_bound == 0
        assert d.get_probability == -1
        assert d.get_parameter == 2
        assert d.get_ident == 2

        d = Poisson(0, 2)
        assert list(d.simulate(1000))

        m = d.simulate(1000).extract_model()
        assert isinstance(m, _DiscreteParametricModel)

    def test_neg_binomial(self):

        d = Distribution("NEGATIVE_BINOMIAL", 0, 1, 0.5)
        assert list(d.simulate(1000))
        assert d.get_sup_bound == -1
        assert d.get_inf_bound == 0
        assert d.get_probability == 0.5
        assert d.get_parameter == 1
        assert d.get_ident == 3
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
        """simulate a vector of realizations"""
        d = Uniform(0, 1)
        try:
            import numpy
            l = numpy.zeros(100000)
            for i in range(len(l)):
                l[i] = d.simulation()
            assert((0 < sum(l)) & (sum(l) < len(l)))
        except:
            pass
    def test_getters(self):
        dist = Distribution(Histogram([1, 1, 1, 2, 2, 2, 3, 3, 3]))
        assert dist.get_mean == 2
        assert 0.3333 < dist.get_max < 0.3334
        assert dist.get_mean == 2
        assert 0.6666 < dist.get_variance < 0.6667






if __name__ == "__main__":
    runTestClass(Test())
    runTestClass(TestDistribution())
