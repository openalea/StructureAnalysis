""" distance matrix tests


:Author: Thomas Cokelaer, Thomas.Cokelaer@inria.fr

"""
__version__ = "$Id: test_distance_matrix.py 9380 2010-08-06 17:53:45Z cokelaer $"


from openalea.stat_tool.vectors import Vectors, VectorDistance
from openalea.stat_tool.comparison import Compare
from openalea.stat_tool.data_transform import SelectVariable
from openalea.stat_tool.cluster import Clustering, ToDistanceMatrix

from tools import interface
from tools import runTestClass
from openalea.stat_tool import get_shared_data

class Test(interface):

    def __init__(self):
        interface.__init__(self,
                           self.build_data(),
                           "data/distribution1.dist",
                           ToDistanceMatrix)

    def build_data(self):
        vec10 = Vectors(get_shared_data("chene_sessile.vec"))
        vec15 = SelectVariable(vec10, [1, 3, 6], Mode="Reject")
        matrix10 = Compare(vec15, VectorDistance("N", "N", "N"))
        c1 = Clustering(matrix10, "Partition", 3, Prototypes=[1, 3, 12],
                        Algorithm="Divisive")
        return ToDistanceMatrix(c1)

    def test_len(self):
        pass

    def test_empty(self):
        pass
        #self.empty()

    def test_constructor_from_file(self):
        pass
        #self.constructor_from_file()

    def test_constructor_from_file_failure(self):
        pass
        #self.constructor_from_file_failure()

    def test_print(self):
        self.print_data()

    def test_display(self):
        self.display()
        self.display_versus_ascii_write()
        self.display_versus_str()

    def test_plot(self):
        self.plot()

    def test_save(self):
        pass
        #self.save()

    def test_extract(self):
        pass

    def test_extract_data(self):
        pass

    def test_symmetrize(self):
        data = self.data
        assert data.symmetrize()

    def test_unnormalize(self):
        data = self.data
        assert data.symmetrize().unnormalize()

    def test_select_individual(self):
        data = self.data
        keep_false = data.select_individual([1], keep=False)
        keep_true = data.select_individual([1], keep=True)

        assert keep_false.nb_row == 2
        assert keep_false.nb_column == 2

        assert keep_true.nb_row == 1
        assert keep_true.nb_column == 1

    def test_get_distance(self):
        data = self.data
        assert data.get_distance(0, 0)

    def test_get_length(self):
        data = self.data
        assert data.get_length(0, 0)

    def test_get_substitution(self):
        data = self.data
        # no substituion computed for vectors 
        assert data.get_substitution_distance(0, 0) == -1

    def test_get_insertion(self):
        data = self.data
        # no insertion computed for vectors 
        assert data.get_insertion_distance(0, 0) == -1

    def test_get_transposition(self):
        data = self.data
        # no transposition computed for vectors 
        assert data.get_transposition_distance(0, 0) == -1

    def test_get_deletion(self):
        data = self.data
        # no deletion computed for vectors 
        assert data.get_deletion_distance(0, 0) == -1




    def test_get_length_outside_range(self):
        data = self.data
        try:
            data.get_length(-1, 0)
            assert False
        except TypeError:
            assert True



if __name__ == "__main__":
    runTestClass(Test())

