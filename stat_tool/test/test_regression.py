""" Cluster tests"""
__revision__ = "$Id$"


from openalea.stat_tool.regression import Regression
from openalea.stat_tool._stat_tool import _Regression_kernel as RegressionKernel
from openalea.stat_tool.vectors import Vectors

from tools import interface

class TestRegression(interface):
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
        self.v = Vectors([[0, 0], [1, 1], [2, 2], [3, 3]])
        r1 = Regression(self.v, "Linear", 1, 2)

        assert r1.get_nb_vector() == 4
        return r1
   
    def test_get_residuals(self):
        for i in range(0,self.data.get_nb_vector()):
            assert self.data.get_residual(i) == 0

    def tst_get_vectors(self):
        v = self.data.get_vectors()
        for i in range(0, self.data.get_nb_vector()):
            assert v[i] == self.v[i]

    def test_print(self):
        self.print_data()
        
    def test_display(self):
        self.display()
        self.display_versus_ascii_write()
        self.display_versus_str()
        
    def test_len(self):
        """not implemented; irrelevant?"""
        assert self.data.get_nb_vector() == 4
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
        #get back the original vectors
        v = self.v
        #and its regression
        r1 = self.data
        #compare with the direct usage of linear regression
        r = v.linear_regression(1, 2)

        assert r
        assert r1
        assert str(r) == str(r1)

    def test_moving_average(self):
        
        v = self.v
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
        
        v = self.v

        r1 = Regression(v, "NearestNeighbours", 1, 2, 1, Weighting=False)
        r = v.nearest_neighbours_regression(1, 2, 1., False) 
        assert r
        assert r1
        assert str(r) == str(r1)

    def test_badtype(self):

        v = self.v

        try:
            Regression(v, "N", 1, 2, [1, ])
            assert False
        except TypeError:
            assert True


class TestRegressionKernel():

    def __init__(self):
        self.data =  RegressionKernel(4, 0, 10)

    def test_max_value(self):
        assert self.data.get_max_value() == 10
    
    def test_min_value(self):
        assert self.data.get_min_value() == 0
    
    def test_get_ident(self):
        assert self.data.get_ident() == 4

    def others_to_be_done(self):
        #there are other methods that need to be tested with an 
        #appropriate examples:
        #get_point
        #get_regression_df
        #get_residual_df  
        #get_nb_parameter
        #get_parameter
        pass



    
