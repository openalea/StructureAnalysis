"""histogram tests"""
__version__ = "$Id: test_histogram.py 9397 2010-08-10 10:50:07Z cokelaer $"

from openalea.stat_tool.distribution import ToHistogram, Binomial
from openalea.stat_tool.distribution import ToDistribution
from openalea.stat_tool.histogram import Histogram
from openalea.stat_tool.distribution import Distribution

from tools import interface
from tools import runTestClass

from openalea.stat_tool import get_shared_data

class Test(interface):
    """a simple unittest class

    """
    def __init__(self):
        interface.__init__(self,
                           self.build_data(),
                           get_shared_data("peup1.his"),
                           Histogram)

    def build_data(self):
        v = Histogram([0, 1, 2, 3])
        assert v
        return v

    def test_constructor_from_list(self):
        v = Histogram([0, 1, 2, 3])
        assert v

    def test_constructor_from_integers(self):
        v = Histogram(0, 1, 2, 3)
        assert v

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
        v = self.data
        assert len(v) == 4

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

    def test_survival_plot_write(self):
        self.survival_plot_write()

    def test_survival_spreadsheet_write(self):
        self.survival_spreadsheet_write()

    def test_simulate(self):
        pass

    def test_extract(self):
        pass

    def test_extract_data(self):
        """todo : check if this test makes sense"""
        h = Histogram(get_shared_data("meri1.his"))

        e = h.estimate_nonparametric()

        assert e

    def test_survival_ascii_write(self):
        """ test display"""
        h = self.data
        h.survival_ascii_write()

    def test_container(self):
        """ container / iterator"""
        h = Histogram(get_shared_data("meri1.his"))

        assert h[0] == 0
        assert h[10] == 1

    def test_to_histogram(self):
        """Test the ToHistogram function"""

        d = Distribution("NEGATIVE_BINOMIAL", 0, 1, 0.5)
        h = d.simulate(1000)
        d2 = ToDistribution(h)
        assert h and d2

        h2 = ToHistogram(d2)
        assert h2
        assert h == h2

    def test_extract_model(self):
        d = Binomial(0, 10, 0.5)
        d == d.simulate(1000).extract_model()


def test_extract_vec():
    """Extract histogram from Vectors"""
    import openalea.stat_tool as stat_tool
    vs = stat_tool.Vectors("data/cvectors.vec")
    h = vs.extract(4)
    assert h
    print h


if __name__ == "__main__":
    #runTestClass(Test())
    test_extract_vec()

#test_extract_vec()
