"""mixture tests"""
__version__ = "$Id: test_mixture.py 9397 2010-08-10 10:50:07Z cokelaer $"


from openalea.stat_tool.mixture import Mixture
from openalea.stat_tool.data_transform import ExtractDistribution
from openalea.stat_tool.histogram import Histogram
from openalea.stat_tool.distribution import Binomial
from openalea.stat_tool.estimate import Estimate

from tools import interface
from tools import runTestClass

from openalea.stat_tool import get_shared_data

class Test(interface):
    """a simple unittest class"""

    def __init__(self):
        interface.__init__(self,
                           self.build_data(),
                           "data/mixture1.mixt",
                           Mixture)

    def build_data(self):
        d1 = Binomial(0, 12, 0.1)
        d2 = Binomial(0, 12, 0.5)
        d3 = Binomial(0, 12, 0.8)

        mixt = Mixture(0.1, d1, 0.2, d2, 0.7, d3)
        assert mixt.nb_component == 3
        return mixt

    def test_constructor_from_file(self):
        self.constructor_from_file()

    def test_constructor_from_file_failure(self):
        self.constructor_from_file_failure()

    def test_constructor_from_dists_failure(self):
        d1 = Binomial(0, 12, 0.1)
        d2 = Binomial(0, 12, 0.5)
        try:
            _mixt = Mixture(0.1, d1, d2)
            assert False
        except TypeError:
            assert True

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
        sim = self.simulate()
        sim.plot()

    def test_estimate(self):
        sim = self.simulate()
        # 3 Binomial distribution to match th original data
        est = Estimate(sim, "Mixture", "B", "B", "B")
        est.plot()


    def test_extract(self):
        """run and test the extract methods"""

        m = self.data

        assert m.extract_weight() == ExtractDistribution(m, "Weight")
        assert m.extract_mixture() == ExtractDistribution(m, "Mixture")

        assert ExtractDistribution(m, "Component", 1) == Binomial(0, 12, 0.1)
        assert ExtractDistribution(m, "Component", 2) == Binomial(0, 12, 0.5)
        assert ExtractDistribution(m, "Component", 3) == Binomial(0, 12, 0.8)

        assert m.extract_component(1) == Binomial(0, 12, 0.1)
        assert m.extract_component(2) == Binomial(0, 12, 0.5)
        assert m.extract_component(3) == Binomial(0, 12, 0.8)

    def test_extract_data(self):
        """run and test the extract_data methods"""

        h = Histogram(get_shared_data( "meri2.his"))
        m = h.estimate_mixture("B", "NB")

        d = m.extract_data()
        assert d

    def test_truncate(self):
        s = self.data
        _res = s.truncate(4)



if __name__ == "__main__":
    runTestClass(Test())

