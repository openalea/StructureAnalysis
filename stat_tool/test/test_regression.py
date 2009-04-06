""" Cluster tests"""
__revision__ = "$Id: $"


from openalea.stat_tool.regression import Regression
from openalea.stat_tool.vectors import Vectors


class Test:
    """a simple unittest class
    
    Integration test 
    ================
    
    * 'ok' means works and testedPerform test on 
    * 'works' means that the output has b=not been tested yet
    
    ========================    ==================================
    ** from the interface**
    ascii_write(False)          ok
    display                     ok
    file_ascii_write            ok
    get_plotable                ?
    old_plot                    ?
    plot                        gnuplot
    plot_print()                ok
    plot_write                  ok  
    save                        ok
    spreadsheet_write           ok
    file_ascii_write            ok
    
    str                 ok
    """
    
    def __init__(self):
        self.data = None
        self.filename = None
        self.Structure = Regression
    
    
    
    
    def test_linear_regression(self):

        v = Vectors([[0,0], [1,1], [2,2], [3,3]])
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



if __name__=="__main__":
    # perform all the test in the class Test (unit tests)
    test = Test()
    for method in dir(test):
        if method.startswith('_'):
            continue
        if callable(getattr(test, method)):
            getattr(test, method)()
        else:
            print 'skipping'
    # and functional tests.    




