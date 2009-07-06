"""Convolution tests

.. author:: Thomas Cokelaer, Thomas.Cokelaer@inria.fr

"""
__revision__ = "$Id$"


from openalea.stat_tool import Convolution, Distribution
from openalea.stat_tool.distribution import Binomial, NegativeBinomial
from openalea.stat_tool.data_transform import ExtractDistribution
from openalea.stat_tool.plot import DISABLE_PLOT
from openalea.stat_tool.output import Display, Save
from openalea.stat_tool import Estimate
from openalea.stat_tool import Simulate

from tools import interface

class Test(interface):
    """a simple unittest class
    
    """
 
    def __init__(self):
        interface.__init__(self, 
            self.build_data(), 
            "data/convolution1.conv", 
            Convolution)
    
    def build_data(self):
        d1 = Binomial(0, 10, 0.5)
        d2 = NegativeBinomial(0, 1, 0.1)
        conv = Convolution(d1, d2) 
        return conv
    
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
        assert len(c) == 2
        assert len(c) == c.nb_distribution()

    def test_plot(self):
        d = self.data
        d.plot()
        

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
        assert m.extract_convolution() == ExtractDistribution(m, "Convolution")
        assert m.extract(1) == Binomial(0, 10, 0.5)
        assert m.extract(2) == NegativeBinomial(0, 1, 0.1)
        assert ExtractDistribution(m, "Elementary", 1) == Binomial(0, 10, 0.5)
        assert ExtractDistribution(m, "Elementary", 2) == NegativeBinomial(0, 1, 0.1)

    def test_extract_data(self):
        """run and test the extract_data methods"""

        m = self.data
        s = m.simulate(1000)
        e = s.estimate_convolution(Binomial(0, 10, 0.5), Estimator="Parametric")
        d = e.extract_data()
        assert d        
        eprime = Estimate(s, "CONVOLUTION", Binomial(0, 10, 0.5))
        # todo: check that e and eprime are similar



def test1():
    """to be cleaned"""
    from openalea.stat_tool import ExtractHistogram, ExtractData
    from openalea.stat_tool import ToHistogram
    from openalea.stat_tool import Histogram, Shift, Display, Save, Plot

    _convol0 = Convolution(Distribution("B", 0, 10, 0.5), 
                           Distribution("NB", 0, 10, 0.5))
    convol1 = Convolution("data/convolution1.conv")

    convol_histo1 = Simulate(convol1, 200)

    _histo20 = ExtractHistogram(convol_histo1, "Elementary", 1)

    convol2 = Estimate(convol_histo1, "CONVOLUTION", 
                       ExtractDistribution(convol1, "Elementary", 1),
                       MinInfBound=0)

    _histo21 = ExtractHistogram(ExtractData(convol2), 'Elementary', 1)
    _histo22 = ToHistogram(ExtractDistribution(convol2, 'Elementary', 1))

    histo_b2 = Histogram("data/nothofagus_antarctica_bud_2.his")
    histo_s2 = Histogram("data/nothofagus_antarctica_shoot_2.his")

    # Estimator="Likelihood" (default) / "PenalizedLikelihood" / "Parametric"
    # Si Estimator="PenalizedLikelihood", options supplementaires possibles
    # Penalty="FirstDifference" / "SecondDifference" (default) / "Entropy", Weight,
    # Outside="Zero" (default) / "Continuation" (cf. stat_funs4.cpp).

    # "NP" or "NON_PARAMETRIC" (for estimation only)

    _convol30 = Estimate(Shift(histo_s2, 1), "CONVOLUTION", 
                        Estimate(histo_b2, "NP"), NbIteration=500)
    convol31 = Estimate(Shift(histo_s2, 1), "CONVOLUTION", 
                        Estimate(histo_b2, "NP"), NbIteration=100,
                        Estimator="PenalizedLikelihood", Weight=0.5)
    _convol32 = Estimate(Shift(histo_s2, 1), "CONVOLUTION",
                        Estimate(histo_b2, "NP"),
                        Estimator="Parametric")
    Display(convol31)
    #if DISABLE_PLOT==False:
        #Plot(convol31)
        #Plot(ExtractDistribution(convol31, "Convolution"))
    Save(convol31, "data/nothofagus_antartica_2.xls", Format="SpreadSheet")


