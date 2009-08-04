"""Mixture tests
#
#  Frequency distributions
#
#  Objective: 
    Analyzing the number of nodes of growth units in selected architectural
    position considering the respective roles of preformation and neoformation,
#
#  Methods: comparison tests, one-way variance analysis,
#           estimation of finite mixture of distributions.
#
#  Wild cherry tree: number of nodes per growth unit (GU)
#
#  Data: Dominique Fournier
#
#  meri1.his: order 1,
#  meri1.his: order 2,
#  meri1.his: order 3, GU 1,
#  meri1.his: order 3, GU 2,
#  meri5.his: short shoots.
#
#
#  Poplar: number of nodes per growth unit
#
#  Data: Yves Caraglio and Herve Rey
#
#  peup1.his: order 2,
#  peup2.his: order 3,
#  peup3.his: order 4,
#  peup4.his: order 5,
#  peup5.his: order 3, GU 4,
#  peup6.his: order 3, acrotony.
#
#########################################################################
"""
__version__ = "$Id$"


from openalea.stat_tool.plot import DISABLE_PLOT

from openalea.stat_tool.mvmixture import _MvMixture
from openalea.stat_tool.distribution import Binomial
from openalea.stat_tool.distribution import Poisson

from tools import interface
from tools import runTestClass


class Test(interface):
    """a simple unittest class"""
   
    def __init__(self):
        interface.__init__(self,
                           self.build_data(),
                           "data/mixture_mv1.mixt",
                           _MvMixture)
        
    def build_data(self):
        d11 = Binomial(0, 12, 0.1)
        d12 = Binomial(0, 12, 0.5)
        d13 = Binomial(0, 12, 0.8)

        d21 = Poisson(0, 18.0)
        d22 = Poisson(0, 5.0)
        d23 = Poisson(0, .20)
        
        data = _MvMixture([0.1, 0.2, 0.7], [[d11, d21], [d12, d22], [d13, d23]])
        assert data.nb_component == 3
        assert data.nb_variable == 2
        return data
       
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
        c = self.data
        assert len(c) == 3

    def test_plot(self):
        if DISABLE_PLOT == False:    
            self.data.plot(1)
         #   self.data.plot(2)

    def _test_save(self):
        self.save()
                
    def test_plot_write(self):
        self.plot_write()

    def test_file_ascii_write(self):
        self.file_ascii_write()
      
    def test_spreadsheet_write(self):
        self.spreadsheet_write()
    
    def test_simulate(self):
        self.simulate()     
        
    def test_extract(self):
        pass
    
    def test_extract_data(self):
        pass

    def test_simulate2(self):
        
        dummy = 1
        d11 = Binomial(0, 12, 0.1)
        d12 = Binomial(0, 12, 0.5)
        d13 = Binomial(0, 12, 0.8)
        d21 = Poisson(0, 18.0)
        d22 = Poisson(0, 5.0)
        d23 = Poisson(0, .20)
        m = _MvMixture([0.1, 0.2, 0.7], [[d11, d21], [d12, d22], [d13, d23]])
        v = m.simulate(5000)
        assert v

        m_estim_model = v.mixture_estimation(m, 100, [True, True])
        assert m_estim_model
        m_estim_nbcomp = v.mixture_estimation(2)
        assert m_estim_nbcomp

    def _test_permutation(self):
        data1 = self.data
        
        data2 = data1.state_permutation([0, 2, 1])
        _data3 = data2.state_permutation([0, 1, 2])
        
        assert str(data1)==str(data2)


if __name__ == "__main__":
    runTestClass(Test())
