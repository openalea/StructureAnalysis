"""vectors tests"""
__version__ = "$Id: test_vectors.py 9380 2010-08-06 17:53:45Z cokelaer $"

from openalea.stat_tool.vectors import Vectors
from openalea.stat_tool.vectors import VarianceAnalysis
from openalea.stat_tool.vectors import ContingencyTable
from openalea.stat_tool.vectors import ComputeRankCorrelation

from tools import interface
from tools import runTestClass

from openalea.stat_tool import get_shared_data

class Test(interface):
    """a simple unittest class


    .. todo:: possible issue with save where Format has to be set to "Data".
        in other words, Format=Data is not taken as the dafault value.

    """

    def __init__(self):
        interface.__init__(self,
                           self.build_data(),
                           "data/vectors.vec",
                           Vectors)

        self.vec10 = self.build_data_2()

    def build_data(self):
        v = Vectors([[1, 2, 3], [1, 3, 1]])

        assert 2 == v.nb_vector
        assert 3 == v.nb_variable
        assert [1, 2] == v.get_identifiers()
        assert v
        return v

    def build_data_2(self):
        return Vectors(get_shared_data("chene_sessile.vec"))

    def test_identifiers(self):
        v = Vectors([[1, 2, 3], [4, 5, 6], [7, 8, 9]], Identifiers=[1, 2, 4])
        assert v.get_identifiers() == [1, 2, 4]

    def test_types(self):

        v2 = Vectors([[1,2.,3.],[1,5.,1.]], Identifiers=[1,2])


    def test_constructor_from_file(self):
        self.constructor_from_file()

    def test_constructor_from_file_failure(self):
        self.constructor_from_file_failure()

    def test_constructor_one_variable(self):

        v = Vectors([1,2,3])
        assert v.nb_variable == 3

    def test_constructor_identifiers_failure(self):
        try:
            # should be Vectors([[1,2,3], Identifiers=[1])
            v = Vectors([1,2,3], Identifiers=[1,2])
            assert False
        except ValueError:
            assert True

    def test_print(self):
        self.print_data()

    def test_display(self):
        self.display()
        self.display_versus_ascii_write()
        self.display_versus_str()

    def test_vectors_pylist(self):
        """test vector constructor from list"""
        v = Vectors([[0, 1, 2, 3], [4, 5, 6, 7]])
        assert v
        v = Vectors([[1.2, 0.34], [1.2, 0.34]])
        assert v
        v = Vectors([[1]])
        assert v
        v = Vectors([[0, 1, 2, 3], [4, 5, 6, 7], [1, 2, 3, 4]])
        assert v

        try:
            v = Vectors([[0, 1, 2, 3], [4, 5, 6, 7], [1, 2, 3]])
            assert False
        except:
            assert True

        try:
            v = Vectors([[0, 1, 2, 3], [4, 5, 6, 7], [1.2, 2, 3]])
            assert False
        except:
            assert True



    def test_len(self):
        v = self.data
        assert len(v) == 2
        assert len(v) == v.nb_vector

    def test_plot(self):
        self.plot()

    def test_save(self):
        self.save(Format="Data")

    def test_plot_write(self):
        self.plot_write()

    def test_file_ascii_write(self):
        self.file_ascii_write()

    def test_file_ascii_data_write(self):
        self.file_ascii_data_write()

    def test_spreadsheet_write(self):
        self.spreadsheet_write()

    def test_simulate(self):
        pass

    def test_extract(self):
        """run and test the extract methods"""
        m = self.data
        m.extract(1)

    def test_extract_data(self):
        """not relevant"""
        pass

    def test_vectors_container(self):
        """vector container : len"""
        v = Vectors([[0, 1, 2, 3], [4, 5, 6 , 7]])
        assert len(v) == 2

        for i in v:
            assert len(i) == 4

        assert v[0] == range(4)
        assert v[1][1] == 5

    def test_variance_analysis(self):
        """test VarianceAnalysis"""
        # todo: finalise and make method variance_analysis robust.
        vec10 = self.vec10
        va = VarianceAnalysis(vec10, 1, 4, "O")
        assert vec10.variance_analysis(1, 4, 1, "whatever", "whatever") == \
            str(va)
    
        try:
            va = VarianceAnalysis(vec10, 1, 4, "DUMMY")
            assert False
        except:
            assert True

    def test_contingency_table(self):
        """test contingency table"""
        vec10 = self.vec10
        ct = ContingencyTable(vec10, 1, 4)
        assert ct and str(ct)

        ct2 = vec10.contingency_table(1, 4, "what", "what")
        assert ct == ct2

    def test_rank_computation(self):
        vec10 = self.vec10
        ComputeRankCorrelation(vec10, Type="Kendall", FileName="test")
        #ComputeRankCorrelation(vec10, Type="Spearman", FileName="test")

if __name__ == "__main__":
    test = Test()
    runTestClass(test)
