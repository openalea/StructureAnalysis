"""Convolution tests"""
__revision__ = "$Id: $"

from openalea.stat_tool import Convolution, Distribution
from openalea.stat_tool.distribution import Binomial, NegativeBinomial
from openalea.stat_tool.data_transform import ExtractDistribution


class Test:
    """a simple unittest class"""
 
    def test_empty(self):
        """Test that empty constructor fails"""
        try:
            _m = Convolution()
            assert False
        except Exception:
            assert True

    def test_constructor_fromfile(self):
        """run constructor with filename argument"""
        c = Convolution("convolution1.conv")
        assert c
        
        return c

    def test_constructor_fromfile_failure(self):
        """run constructor with filename argument"""
        try:
            _h = Convolution("convolution1.con")
            assert False
        except Exception:
            assert True

    def test_build_convolution(self):
        """run constructor with filename two distributions as arguments"""

        d1 = Binomial(0, 10, 0.5)
        d2 = NegativeBinomial(0, 1, 0.1)

        m = Convolution(d1, d2)
        
        assert m
        return m

    def test_plot(self):
        """run plotting routines """
        m = self.test_build_convolution()
        m.plot()

        assert str(m)
        m.display()

    def test_simulation(self):
        """Test the simulate method"""

        m = self.test_build_convolution()
        s = m.simulate(1000)

        assert len(s) == 1000
        assert s.nb_histogram() == 2
        assert str(s)

    def test_extract(self):
        """run and test the extract methods"""
      
        m = self.test_build_convolution()
        assert m.nb_distribution() == 2

        assert m.extract_convolution() == ExtractDistribution(m, "Convolution")

        assert m.extract_elementary(1) == Binomial(0, 10, 0.5)
        assert m.extract_elementary(2) == NegativeBinomial(0, 1, 0.1)

        assert ExtractDistribution(m, "Elementary", 1) == Binomial(0, 10, 0.5)
        assert ExtractDistribution(m, "Elementary", 2) == NegativeBinomial(0, 1, 0.1)

    def test_extract_data(self):
        """run and test the extract_data methods"""
        #from openalea.stat_tool.distribution import Binomial

        m = self.test_build_convolution()
        s = m.simulate(1000)

        m = s.estimate_convolution(Binomial(0, 10, 0.5), Estimator="Parametric")

        d = m.extract_data()
        assert d



def test1():
    from openalea.stat_tool import Estimate, Simulate
    from openalea.stat_tool import ExtractHistogram, ExtractData
    from openalea.stat_tool import ToHistogram
    from openalea.stat_tool import Histogram, Shift, Display, Save
  
    convol1 = Convolution(Distribution("B", 0, 10, 0.5), Distribution("NB", 0, 10, 0.5))
    convol1 = Convolution("convolution1.conv")
    
    convol_histo1 = Simulate(convol1, 200)
    
    histo20 = ExtractHistogram(convol_histo1, "Elementary", 1)
    
    convol2 = Estimate(convol_histo1, "CONVOLUTION", 
                       ExtractDistribution(convol1, "Elementary", 1), MinInfBound=0)
    
    histo21 = ExtractHistogram(ExtractData(convol2), 'Elementary', 1)
    histo22 = ToHistogram(ExtractDistribution(convol2,'Elementary', 1))
    
    histo_b2 = Histogram("nothofagus_antarctica_bud_2.his")
    histo_s2 = Histogram("nothofagus_antarctica_shoot_2.his")
    
    # Estimator="Likelihood" (default) / "PenalizedLikelihood" / "Parametric"
    # Si Estimator="PenalizedLikelihood", options supplementaires possibles
    # Penalty="FirstDifference" / "SecondDifference" (default) / "Entropy", Weight,
    # Outside="Zero" (default) / "Continuation" (cf. stat_funs4.cpp).
    
    # "NP" or "NON_PARAMETRIC" (for estimation only)
    
    convol30 = Estimate(Shift(histo_s2, 1), "CONVOLUTION", Estimate(histo_b2, "NP"), 
                        NbIteration=500)
    convol31 = Estimate(Shift(histo_s2, 1), "CONVOLUTION", Estimate(histo_b2, "NP"), 
                        NbIteration=100, Estimator="PenalizedLikelihood", Weight=0.5)
    convol32 = Estimate(Shift(histo_s2, 1), "CONVOLUTION", Estimate(histo_b2, "NP"), 
                        Estimator="Parametric")
    Display(convol31)
    # Plot(convol31)
    # Plot(ExtractDistribution(convol31, "Convolution"))
    Save(convol31, "nothofagus_antartica_2.xld", Format="SpreadSheet")




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
    test1()
