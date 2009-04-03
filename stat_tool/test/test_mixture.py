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
__revision__ = "$Id: $"

import os

from openalea.stat_tool.plot import DISABLE_PLOT

from openalea.stat_tool.mixture import Mixture, _MvMixture
from openalea.stat_tool.data_transform import ExtractDistribution
from openalea.stat_tool.histogram import Histogram
from openalea.stat_tool.distribution import Uniform, Binomial, Poisson
from openalea.stat_tool.estimate import Estimate
from openalea.stat_tool.comparison import Compare, ComparisonTest

from openalea.stat_tool.output import plot
from openalea.stat_tool.data_transform import Merge, Shift, ExtractData
from openalea.stat_tool.simulate import Simulate
import openalea.stat_tool.distribution as distribution 
from openalea.stat_tool.cluster import Cluster

from tools import interface

class Test(interface):
    """a simple unittest class
    
    Integration test 
    ================
    
    * 'ok' means works and testedPerform test on 
    * 'works' means that the output has b=not been tested yet
    
    ========================    ==================================
    ** from the interface**
    ascii_write                 ok
    display                     ok    
    extract_data                ok
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
    extract_mixture             ok
    extratc_component           ok
    extra_weight                ok
    str                         ok
    len                         not relevant
    nb_component                ok
    old_plot                    works   
    ========================    ==================================    
    """

    def __init__(self):
        self.data = self.build_data()
        self.filename = "mixture1.mixt"
        self.structure = Mixture
        
    def build_data(self):
        d1 = Binomial(0, 12, 0.1)
        d2 = Binomial(0, 12, 0.5)
        d3 = Binomial(0, 12, 0.8)
        
        data = Mixture(0.1, d1, 0.2, d2, 0.7, d3)
        assert data.nb_component() == 3
        return data
            
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
        c = self.data
        assert len(c) == 3

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
    
    def test_simulate(self):
        self.simulate()     
        
    def test_extract(self):
        """run and test the extract methods"""

        m = self.data

        assert m.extract_weight() == ExtractDistribution(m, "Weight")
        assert m.extract_mixture() == ExtractDistribution(m, "Mixture")

        assert ExtractDistribution(m, "Component", 1) == Binomial(0, 12, 0.1)
        assert ExtractDistribution(m, "Component", 2) == Binomial(0, 12, 0.5)
        assert ExtractDistribution(m, "Component", 3) == Binomial(0, 12, 0.8)

        assert m.extract_component(1) == Binomial(0, 12, 0.1)
        assert m.extract_component(2) == Binomial(0, 12, 0.5)
        assert m.extract_component(3) == Binomial(0, 12, 0.8)

    def test_extract_data(self):
        """run and test the extract_data methods""" 

        h = Histogram("meri2.his")
        m = h.estimate_mixture(["B", "NB"])

        d = m.extract_data()
        assert d



def test1():
    plot.DISABLE_PLOT = DISABLE_PLOT
 
    meri1 = Histogram("meri1.his")
    meri2 = Histogram("meri2.his")
    meri3 = Histogram("meri3.his")
    meri4 = Histogram("meri4.his")
    meri5 = Histogram("meri5.his")

    #Plot(meri1, meri2, meri3, meri4, meri5)
    Compare(meri1, meri2, meri3, meri4, meri5, "N")


    ComparisonTest("F", meri1, meri2)
    ComparisonTest("T", meri1, meri2)
    ComparisonTest("W", meri1, meri2)

    ComparisonTest("F", meri1, meri3)
    ComparisonTest("T", meri1, meri3)
    ComparisonTest("W", meri1, meri3)

    # estimation of a mixture of two distributions assuming a first sub-population of GUs
    # made only of a preformed part and a second sub-population made of both a preformed part
    # and a neoformed part
    
    _mixt1 = Estimate(meri2, "MIXTURE", "B", "B")
    
    meri = Merge(meri1, meri2, meri3, meri4, meri5)
    
    # model selection approach: estimation of both the mixture parameters and
    # the number of components 
    
    mixt2 = Estimate(meri, "MIXTURE", "B", "B", "B", "B",  NbComponent="Estimated")
    # mixt2 = Estimate(meri, "MIXTURE", "NB", "NB")
    # Plot(ExtractDistribution(mixt2, "Mixture"))
    # Display(mixt2)
    
    _mixt_data = ExtractData(mixt2)
    
    
    dist5 = Estimate(meri5, "BINOMIAL")
    # Display(dist5, Detail->2)
    # Plot(dist5)
    
    histo5 = Simulate(dist5, 100)
    # Display(histo5, Detail->2)
    # Plot(histo5)
    
    
    peup1 = Histogram("peup1.his")
    peup2 = Histogram("peup2.his")
    peup3 = Histogram("peup3.his")
    peup4 = Histogram("peup4.his")
    peup5 = Histogram("peup5.his")
    peup6 = Histogram("peup6.his")
    
    mixt10 = Estimate(peup2, "MIXTURE", "B", "NB", "NB", "NB", NbComponent="Estimated")
    
    peup = Merge(peup1, peup2, peup3, peup4, peup5, peup6)
    
    histo1 = Shift(peup, -1)
    histo2 = Cluster(peup, "Information", 0.8)
    histo3 = Cluster(peup, "Step", 10)
    histo4 = Cluster(peup, "Limit", [13, 24])
    # Display(histo4, Detail->2)
    # Plot(histo4)


    mixt11 = Estimate(peup, "MIXTURE", "B", "NB", "NB", "NB", NbComponent="Estimated")
    # mixt11 = Estimate(peup, "MIXTURE", "B", "NB")

    d11 = distribution.Binomial(0, 12, 0.1)
    d12 = distribution.Binomial(2, 13, 0.6)
    d13 = distribution.Binomial(3, 15, 0.9)
    
    d21 = distribution.Poisson(0, 25.0)
    d22 = distribution.Poisson(0, 5.0)
    d23 = distribution.Poisson(0, 0.2)

    m = _MvMixture([0.1, 0.2, 0.7], [[d11, d21], [d12, d22], [d13, d23]])
    print m

    #m2 = _MvMixture("mixture_mv1.mixt")
    #print m2

    #print "Egalite des melanges construits par liste ",\
    #  "de distributions et par fichiers : ", str(str(m)==str(m2))

    #m = _MvMixture("mixture_mv_nonparam.mixt")
   # print m

    #print "Simulation de melanges multivaries : "
    #v = m.simulate(5000)
    #print v

    #m.plot(variable=1, Title="Simulated mixture")

    #print "Estimation de melanges multivaries ", \
    #    "d'apres un modele initial : "
    #m_estim_model = v.mixture_estimation(m, 100,  [True, True])
      
    #m_estim_model.plot(variable = 1, Title="Estimated mixture")
    
    #print "Estimation de melanges multivaries ", \
    #    "d'apres un nombre de composantes : "
        
    #m_estim_nbcomp = v.mixture_estimation(3, 100, [True, True])
    
    #m_estim_nbcomp.plot(variable = 1, Title="Estimated mixture")


