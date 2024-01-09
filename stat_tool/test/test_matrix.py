"""vectors tests"""
__version__ = "$Id: test_matrix.py 9380 2010-08-06 17:53:45Z cokelaer $"

from openalea.stat_tool.vectors import Vectors, VectorDistance
from openalea.stat_tool.data_transform import SelectVariable
from openalea.stat_tool.comparison import Compare

from tools import interface
from tools import runTestClass
from openalea.stat_tool import get_shared_data

class Test(interface):
    """a simple unittest class for distance matrix"""
    def __init__(self):
        interface.__init__(self,
                           self.build_data(),
                           None,
                           Vectors)
        # init expect the 4th argument to be provided.
        # vectors is therefore passed as dummy structure
    def build_data(self):
        vec10 = Vectors(get_shared_data("chene_sessile.vec"))
        vec15 = SelectVariable(vec10, [1, 3, 6], Mode="Reject")
        matrix10 = Compare(vec15, VectorDistance("N", "N", "N"))
        assert 138 == matrix10.nb_row
        assert 138 == matrix10.nb_column
        return matrix10

    def test_empty(self):
        """no construtor"""
        pass

    def test_constructor_from_file(self):
        """no construtor"""
        pass

    def test_constructor_from_file_failure(self):
        """no construtor"""
        pass

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

    def _test_save(self):
        pass
        #self.save()

    def test_plot_write(self):
        self.plot_write()

    def test_file_ascii_write(self):
        self.file_ascii_write()

    def test_spreadsheet_write(self):
        self.spreadsheet_write()

    def test_simulate(self):
        """no simulate method"""
        pass

    def test_extract(self):
        """no such method"""
        pass

    def test_extract_data(self):
        """no such method"""
        pass


 
if __name__ == "__main__":
    runTestClass(Test())  
