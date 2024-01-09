"""Compound tests

:Author: Thomas Cokelaer, Thomas.Cokelaer@inria.fr

"""
__version__ = "$Id: test_compound.py 8625 2010-04-13 11:34:37Z cokelaer $"

from openalea.stat_tool.compound import Compound
from openalea.stat_tool.data_transform import ExtractDistribution
from openalea.stat_tool.distribution import Binomial, NegativeBinomial
from openalea.stat_tool.estimate import Estimate

from tools import interface
from tools import runTestClass


class Test(interface):
    """a simple unittest class

    See tools to get documentation of the following methods
    """

    def __init__(self):
        interface.__init__(self,
            self.build_data(),
            "data/compound1.cd",
            Compound)

    def build_data(self):
        d1 = Binomial(2, 5, 0.5)
        d2 = NegativeBinomial(0, 2, 0.5)
        comp = Compound(d1, d2)
        return comp

    def test_constructor_from_file(self):
        self.constructor_from_file()

    def test_constructor_from_file_failure(self):
        self.constructor_from_file_failure()

    def test_constructor_from_compound(self):
        # to be removed
        #compound1 = self.data
        #compound2 = Compound(compound1)
        #assert str(compound1) == str(compound2)
        pass

    def test_constructor_from_dists_and_threshold(self):
        compound1 = self.data
        compound2 = Compound(Binomial(2, 5, 0.5),
                             NegativeBinomial(0, 2, 0.5),
                             Threshold=0.99999)
        assert str(compound1) == str(compound2)

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
        sim = self.simulate()
        sim.plot()

    def test_extract(self):
        """run and test the extract methods"""
        m = self.data
        assert m.extract_compound() == ExtractDistribution(m, "Compound")
        assert m.extract_sum() == Binomial(2, 5, 0.5)
        assert m.extract_sum() == ExtractDistribution(m, "Sum")
        assert m.extract_elementary() == NegativeBinomial(0, 2, 0.5)
        assert m.extract_elementary() == ExtractDistribution(m, "Elementary")

    def test_extract_data(self):
        """todo : check if this test makes sense"""

        s = self.simulate()
        #e = Estimate(s, "Compound",  Binomial(2, 5, 0.5), "Sum")
        d = s.extract_sum()
        assert d
        _eprime = Estimate(s, "COMPOUND", Binomial(0, 10, 0.5), "Sum")

    def test_truncate(self):
        s = self.data
        _res = s.truncate(4)


if __name__ == "__main__":
    runTestClass(Test())
