"""histogram tests"""
__revision__ = "$Id: $"

from openalea.stat_tool.distribution import ToHistogram, Binomial
from openalea.stat_tool.distribution import ToDistribution
from openalea.stat_tool.histogram import Histogram
from openalea.stat_tool.distribution import Distribution
from openalea.stat_tool.output  import Display, Save

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
    plot                        ok                      
    save                        ok
    plot_print                  ok
    simulate                    ok
    plot_write                  ok
    spreadsheet_write           ok
    survival_ascii_write        ok
    survival_spreadsheet_write  ok

    **others**
    extract_model               ok

    ** see test_cluster**
    cluster_information         ok
    cluster_limit               ok
    cluster_step                ok
    transcode                   ok

    **comparison**
    compare                     ???
    compare_histo               ok
    t_comparison                ok
    wmw_comparison              ok
    f_comparison                ok
    
    
    survival_get_plotable        ???
    
    **see test_estimate**
    estimate_compound            ok
    estimate_convolution         ok     
    estimate_mixture             ok 
    estimate_nonparametric       ok
    estimate_parametric          ok
    compound_estimation          equivalent to estimate_mixture and not tested
    convolution_estimation       equivalent to estimate_mixture and not tested
    mixture_estimation           equivalent to estimate_mixture and tested
    parametric_estimation        equivalent to estimate_mixture and not tested
      
    **see data_transform**
    fit                            ok
    merge                          ok
    shift                          ok
    value_select                   ok
    """
    def __init__(self):
        self.data = self.build_data()
        self.filename = "peup1.his"
        self.structure = Histogram
        
    def build_data(self):
        v = Histogram([0, 1, 2, 3])
        assert v
        return v
      
    def test_empty(self):
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
        assert len(v) == 4

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
        
    def test_plot_write(self):
        self.plot_write()

    def test_file_ascii_write(self):
        self.file_ascii_write()
      
    def test_spreadsheet_write(self):
        self.spreadsheet_write()
    
    def test_simulate(self):
        pass
     
    def test_extract(self):
        pass
    
    def test_extract_data(self):
        """todo : check if this test makes sense"""
        h = Histogram("meri1.his")
        e = h.estimate_nonparametric()
        
        assert e
    
    def test_survival_ascii_write(self):
        """ test display"""
        h = self.data
        h.survival_ascii_write()
        
    def test_container(self):
        """ container / iterator"""
        h = Histogram("meri1.his")

        assert h[0] == 0
        assert h[10] == 1

    def test_to_histogram(self):
        """Test the ToHistogram function"""
        
        d = Distribution("NEGATIVE_BINOMIAL", 0, 1, 0.5)
        h = d.simulate(1000)
        d2 = ToDistribution(h)
        assert h and d2

        h2 = ToHistogram(d2)
        assert h2
        assert h == h2

    def test_extract_model(self):    
        d = Binomial(0,10,0.5)
        d == d.simulate(1000).extract_model()
        


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


