"""vectors tests"""
__revision__ = "$Id: $"

from openalea.stat_tool.vectors import Vectors, VectorDistance, \
    VarianceAnalysis, ContingencyTable
from openalea.stat_tool.data_transform import ValueSelect
from openalea.stat_tool.regression import Regression
from openalea.stat_tool.data_transform import ExtractHistogram, \
    SelectIndividual, SelectVariable, ExtractDistribution
from openalea.stat_tool.comparison import Compare
from openalea.stat_tool.estimate import Estimate
from openalea.stat_tool.cluster import Clustering
from openalea.stat_tool.output import Display, Plot, Save

from openalea.stat_tool.plot import DISABLE_PLOT

from tools import interface

class Test(interface):
    """a simple unittest class
    
    Integration Test 
    ================
    
    * 'ok' means works and testedPerform test on 
    * 'works' means that the output has b=not been tested yet
    
    ========================    ==================================
    ** from the interface**
    ascii_write                 ok
    display                     ok
    extract_data                nothing to be done
    file_ascii_write            ok
    get_plotable                what is it for ?     
    plot                        Not working                       
    save                        ok
    plot_print                  ok
    simulate                    ok
    plot_write                  ok
    spreadsheet_write           ok
    survival_ascii_write        ok
    survival_spreadsheet_write  ok
   
    str                         ok
    len                         ok 
    """

    def __init__(self):
        self.data = self.build_data()
        self.filename = "vector_distance.vd"
        self.structure = VectorDistance
 
        
    def build_data(self):
        v = VectorDistance('N', 'O', 'S') 
        assert v 
        return v
        
    def _test_empty(self):
        self.empty()
        
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

    def _test_plot(self):
        """run plotting routines """
        # todo: does not produce anythinh but expected ?
#        c = self.test_build_vectors()
#        if DISABLE_PLOT==False:
#            c.plot()
        pass
    
    def test_save(self):
        self.save()

    def _test_plot_write(self):
        self.plot_write()

    def test_file_ascii_write(self):
        self.file_ascii_write()
      
    def _test_spreadsheet_write(self):
        self.spreadsheet_write()
    
    def test_simulate(self):
        pass
        
    def test_extract(self):
        pass
    
    def test_extract_data(self):
        pass
    
    def test_vector_distance(self):
        """ test vector distance"""
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


