"""Convolution tests

:Author: Thomas Cokelaer, Thomas.Cokelaer@inria.fr

"""
__version__ = "$Id: test_convolution.py 8625 2010-04-13 11:34:37Z cokelaer $"

from openalea.stat_tool import Convolution
from openalea.stat_tool.distribution import Binomial, NegativeBinomial
from openalea.stat_tool.data_transform import ExtractDistribution
from openalea.stat_tool import Estimate

from tools import interface
from tools import runTestClass


class Test(interface):
    """a simple unittest class

    See tools to get documentation of the following methods
    """

    def __init__(self):
        interface.__init__(self,
            self.build_data(),
            "data/convolution1.conv",
            Convolution)

    def build_data(self):
        d1 = Binomial(0, 10, 0.5)
        d2 = NegativeBinomial(0, 1, 0.1)
        conv = Convolution(d1, d2, d1, d2)
        conv = Convolution(d1, d2)
        return conv

    def test_constructor_from_file(self):
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
        assert len(c) == 2
        assert len(c) == c.nb_distribution()

    def test_plot(self):
        d = self.data
        d.plot()


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

    def test_extract(self):
        """run and test the extract methods"""

        m = self.data
        assert m.extract_convolution() == ExtractDistribution(m, "Convolution")
        assert m.extract(1) == Binomial(0, 10, 0.5)
        assert m.extract(2) == NegativeBinomial(0, 1, 0.1)
        assert ExtractDistribution(m, "Elementary", 1) == Binomial(0, 10, 0.5)
        assert ExtractDistribution(m, "Elementary", 2) == \
            NegativeBinomial(0, 1, 0.1)

    def test_extract_data(self):
        """run and test the extract_data methods"""

        m = self.data
        s = m.simulate(1000)
        e = s.estimate_convolution(Binomial(0, 10, 0.5), Estimator="Parametric")
        d = e.extract_data()
        assert d
        eprime = Estimate(s, "CONVOLUTION", Binomial(0, 10, 0.5))

        # todo: find robust assert ?
        print eprime.get_mean

    def test_truncate(self):
        s = self.data
        _res = s.truncate(4)



if __name__ == "__main__":
    runTestClass(Test())

