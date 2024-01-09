"""vectors tests

:Author: Thomas Cokelaer, Thomas.Cokelaer@inria.fr

"""
__version__ = "$Id: test_vector_distance.py 8938 2010-05-19 15:15:06Z cokelaer $"


from openalea.stat_tool.vectors import VectorDistance

from tools import interface
from tools import runTestClass


class Test(interface):
    """a simple unittest class"""

    def __init__(self):
        interface.__init__(self,
                           self.build_data(),
                           "data/vector_distance.vd",
                           VectorDistance)
    def build_data(self):
        v = VectorDistance('N', 'O', 'S')
        assert v.distance_type == 0
        assert v
        return v

    def test_empty(self):
        """skip test_empty
        because there is an empty constructor in VectorDistance; """
        pass

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
        assert len(v) == 3
        assert len(v) == v.nb_variable

    def test_plot(self):
        pass

    def test_save(self):
        self.save()

    def _test_plot_write(self):
        self.plot_write()

    def test_file_ascii_write(self):
        self.file_ascii_write()

    def _test_spreadsheet_write(self):
        #not implemented
        self.spreadsheet_write()

    def test_simulate(self):
        pass

    def test_extract(self):
        pass

    def test_extract_data(self):
        pass

    def test_vector_distance(self):
        """ test vector distance constructors"""
        v = VectorDistance('NUMERIC', 'ORDINAL', 'SYMBOLIC')
        assert v and len(v) == 3

        v = VectorDistance(2.3, 'N', 4, 'O', 6, 'S')
        assert v and len(v) == 3

        # wrapper should work. needs python interface ?
        #v = VectorDistance( (2.3, 'N'),  (4, 'O'), (6, 'S'))
        #assert v and len(v) == 3

        v = VectorDistance(2.3, 'N', 4, 'O', 6, 'S', Distance="QUADRATIC")
        assert v and len(v) == 3

        v = VectorDistance('NUMERIC', 'ORDINAL', 'SYMBOLIC',
                           Distance="QUADRATIC")
        assert v and len(v) == 3

        assert str(VectorDistance('N', 'O', 'S'))

    def test_constructor_failure(self):
        try:
            v = VectorDistance('N', 'DUMMY')
            assert False
        except:
            assert True

if __name__ == "__main__":
    runTestClass(Test())
