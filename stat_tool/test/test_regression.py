""" Cluster tests"""
__revision__ = "$Id$"


from openalea.stat_tool.regression import Regression
from openalea.stat_tool.vectors import Vectors

from tools import interface

class Test(interface):
    """a simple unittest class
    
    Integration test 
    ================
    
    * 'ok' means works and tested 
    
    ========================    ==================================
    ** from the interface**
    ascii_write(False)          ok
    display                     ok
    file_ascii_write            ok
    old_plot                    ok
    plot                        ok
    plot_print()                ok
    plot_write                  ok  
    save                        ok
    spreadsheet_write           ok
    file_ascii_write            ok
    str                         ok
    """
      
    def __init__(self):
        interface.__init__(self, self.build_data(), None, Regression)
    
    def build_data(self):
        v = Vectors([[0, 0], [1, 1], [2, 2], [3, 3]])
        r1 = Regression(v, "Linear", 1, 2)
        return r1
    
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
        self.save(skip_reading=True)
                
    def test_plot_write(self):
        self.plot_write()

    def test_file_ascii_write(self):
        self.file_ascii_write()
      
    def test_spreadsheet_write(self):
        self.spreadsheet_write()

    def test_extract(self):
        pass
    
    def test_extract_data(self):
        pass

    def test_linear_regression(self):

        v = Vectors([[0, 0], [1, 1], [2, 2], [3, 3]])
        r1 = Regression(v, "Linear", 1, 2)
        r = v.linear_regression(1, 2)

        assert r
        assert r1
        assert str(r) == str(r1)

    def test_moving_average(self):
        
        v = Vectors([[0, 0], [1, 1], [2, 2], [3, 3]])

        # Test algorithm
        try:
            r = v.moving_average_regression(1, 2, [1, ], 'n') 
            assert False
        except:
            assert True

        r1 = Regression(v, "MovingAverage" , 1, 2, [1, ])
        r = v.moving_average_regression(1, 2, [1, ], 'a') 
        assert r
        assert r1
        assert str(r)==str(r1)
       
    def test_nearest_neighbours(self):
        
        v = Vectors([[0, 0], [1, 1], [2, 2], [3, 3]])

        r1 = Regression(v, "NearestNeighbours", 1, 2, 1, Weighting=False)
        r = v.nearest_neighbours_regression(1, 2, 1., False) 
        assert r
        assert r1
        assert str(r) == str(r1)

    def test_badtype(self):

        v = Vectors([[0, 0], [1, 1], [2, 2], [3, 3]])

        try:
            Regression(v, "N", 1, 2, [1, ])
            assert False
        except TypeError:
            assert True


