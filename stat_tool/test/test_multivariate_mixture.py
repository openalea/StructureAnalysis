# -*- coding: utf-8 -*-
"""Multivariate mixture tests
#
#
#########################################################################
"""
__version__ = "$Id$"


from openalea.stat_tool.plot import DISABLE_PLOT

from openalea.stat_tool.multivariate_mixture import _MultivariateMixture
from openalea.stat_tool.distribution import Binomial
from openalea.stat_tool.distribution import Poisson

from tools import interface
from tools import runTestClass


class Test(interface):
    """a simple unittest class"""

    def __init__(self):
        interface.__init__(self,
                           self.build_data(),
                           "data/mixture_mv1.mixt",
                           _MultivariateMixture)

    def build_data(self):
        d11 = Binomial(0, 12, 0.1)
        d12 = Binomial(0, 12, 0.5)
        d13 = Binomial(0, 12, 0.8)

        d21 = Poisson(0, 18.0)
        d22 = Poisson(0, 5.0)
        d23 = Poisson(0, .20)

        data = _MultivariateMixture([0.1, 0.2, 0.7], [[d11, d21], [d12, d22], [d13, d23]])
        assert data.nb_component == 3
        assert data.nb_variable == 2
        return data

    def test_empty(self):
        """skip test_empty

        because there is an empty constructor in VectorDistance; """
        pass

    def _test_constructor_from_file(self):
        #raise error (proba non equal to 1) when nosetests used from parent directory.
        self.constructor_from_file()

    def test_constructor_from_file_failure(self):
        self.constructor_from_file_failure()

    def test_print(self):
        self.print_data()

    def test_display(self):
        self.display()
        self.display_versus_ascii_write()
        self.display_versus_str()

    def test_len(self):
        c = self.data
        assert len(c) == 3

    def test_plot(self):
        if DISABLE_PLOT == False:
            self.data.plot(1)
         #   self.data.plot(2)

    def _test_save(self):
        self.save()

    def test_plot_write(self):
        self.plot_write()

    def test_file_ascii_write(self):
        self.file_ascii_write()

    def test_spreadsheet_write(self):
        self.spreadsheet_write()

    def test_simulate(self):
        self.simulate()

    def test_extract(self):
        pass

    def test_extract_data(self):
        pass

    def test_simulate2(self):

        d11 = Binomial(0, 12, 0.1)
        d12 = Binomial(0, 12, 0.5)
        d13 = Binomial(0, 12, 0.8)
        d21 = Poisson(0, 18.0)
        d22 = Poisson(0, 5.0)
        d23 = Poisson(0, .20)
        m = _MultivariateMixture([0.1, 0.2, 0.7], [[d11, d21], [d12, d22], [d13, d23]])
        v = m.simulate(5000)
        assert v

        m_estim_model = v.mixture_estimation(m, 100, [True, True])
        assert m_estim_model
        m_estim_nbcomp = v.mixture_estimation(2)
        assert m_estim_nbcomp

    def _test_permutation(self):
        data1 = self.data

        data2 = data1.state_permutation([0, 2, 1])
        _data3 = data2.state_permutation([0, 1, 2])

        assert str(data1)==str(data2)


if __name__ == "__main__":
    test = Test()
    runTestClass(test)
