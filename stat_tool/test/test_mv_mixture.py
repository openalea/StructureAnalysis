"""Mixture tests
#
#  Frequency distributions
#
#  Objective: Analyzing the number of nodes of growth units in selected architectural
#             position considering the respective roles of preformation and neoformation,
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
__revision__ = "$Id$"

import os

from openalea.stat_tool.plot import DISABLE_PLOT

from openalea.stat_tool.mixture import Mixture, _MvMixture
from openalea.stat_tool.data_transform import ExtractDistribution
from openalea.stat_tool.histogram import Histogram
from openalea.stat_tool.distribution import Uniform, Binomial, Poisson
from openalea.stat_tool.estimate import Estimate
from openalea.stat_tool.comparison import Compare, ComparisonTest

from tools import interface
 

class Test(interface):
    """a simple unittest class
    
    Integration test 
    ================
    
    * 'ok' means works and testedPerform test on 
    
    ========================    ==================================
    ** from the interface**
    ascii_write                 ok
    display                     ok    
    extract_data                ok
    file_ascii_write            ok
    get_plotable                ?     
    plot                        ok  but only gnuplot implemented                     
    save                        ok
    plot_print                  ok
    simulate                    ok
    plot_write                  ok
    
    extraxt_mixture
    extract_data        
    extract_weight              ok but does it make sense?
    
    str                         ok
    len                         ok
    nb_component                ok
    old_plot                    ok   
    state_permutation           ok
    
    _criteria                   need a test
    _is_parametric              need a test
    
    save_backup                 same as save
    ========================    ==================================    
    """
   
   
    def __init__(self):
        interface.__init__(self,
                           self.build_data(),
                           "mixture_mv1.mixt",
                           _MvMixture)
        
    def build_data(self):
        d11 = Binomial(0, 12, 0.1)
        d12 = Binomial(0, 12, 0.5)
        d13 = Binomial(0, 12, 0.8)

        d21 = Poisson(0, 18.0)
        d22 = Poisson(0, 5.0)
        d23 = Poisson(0, .20)
        
        data = _MvMixture([0.1, 0.2, 0.7], [[d11, d21], [d12, d22], [d13, d23]])
        assert data.nb_component()==3
        assert data.nb_variable()==2
        return data
       
    def _test_empty(self):
        """there is an empty constructor in MvMixture!"""
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
        
        data2 = data1.state_permutation([0,2,1])
        data3 = data2.state_permutation([0,1,2])
        
        assert str(data1)==str(data2)
