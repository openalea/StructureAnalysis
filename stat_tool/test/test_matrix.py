"""vectors tests"""
__revision__ = "$Id$"

from openalea.stat_tool.vectors import Vectors, VectorDistance
from openalea.stat_tool.data_transform import SelectVariable
from openalea.stat_tool.comparison import Compare

from tools import interface


class Test(interface):
    """a simple unittest class

    ** from the interface**
    Constructors from file to be done
    ascii_write                 ok
    display                     ok    
    file_ascii_write            ok
    plot                        no                  
    save                        ok
    plot_print                  ok
    plot_write                  ok
    spreadsheet_write           ok
    **others**
    str                         ok
    len                         no 
    old_plot                    ok
    hierarchical_clustering     ok
    partitioning_clusters       ok
    partitioning_prototype      ok
    get_nb_row                  ok
    get_nb_column               ok

    """
    def __init__(self):
        self.data = self.build_data()
        self.filename = None
        self.structure = None

    def build_data(self):
        vec10 = Vectors("chene_sessile.vec")
        vec15 = SelectVariable(vec10, [1, 3, 6], Mode="Reject")
        matrix10 = Compare(vec15, VectorDistance("N", "N", "N"))
        assert 138 == matrix10.get_nb_row()
        assert 138 == matrix10.get_nb_column()
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

    def test_vector_distance(self):
        v = VectorDistance('N', 'O', 'S') 
        assert v and len(v) == 3

        v = VectorDistance('NUMERIC', 'ORDINAL', 'SYMBOLIC')
        assert v and len(v) == 3

        v = VectorDistance(2.3, 'N', 4, 'O', 6, 'S')
        assert v and len(v) == 3

        v = VectorDistance( (2.3, 'N'),  (4, 'O'), (6, 'S'))
        assert v and len(v) == 3

        v = VectorDistance(2.3, 'N', 4, 'O', 6, 'S', distance = "QUADRATIC")
        assert v and len(v) == 3

        v = VectorDistance('NUMERIC', 'ORDINAL', 'SYMBOLIC', \
                           Distance="QUADRATIC")
        assert v and len(v) == 3

        assert str(VectorDistance('N', 'O', 'S'))


