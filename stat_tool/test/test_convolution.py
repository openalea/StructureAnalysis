"""Convolution tests"""
__revision__ = "$Id: $"


from openalea.stat_tool import Convolution, Distribution
from openalea.stat_tool.distribution import Binomial, NegativeBinomial
from openalea.stat_tool.data_transform import ExtractDistribution
from openalea.stat_tool.plot import DISABLE_PLOT
from openalea.stat_tool.output import Display, Save
from openalea.stat_tool import Estimate
from openalea.stat_tool import Simulate


class TestUnit:
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
        """run constructor with wrong filename argument"""
        try:
            _h = Convolution("convolution1.con")
            assert False
        except Exception:
            assert True

    def test_build_convolution(self):
        """build convolution"""

        d1 = Binomial(0, 10, 0.5)
        d2 = NegativeBinomial(0, 1, 0.1)

        c = Convolution(d1, d2)

        assert c
        return c

    def test_print(self):
        """test that print command exists"""
        c = self.test_build_convolution()
        print c

    def test_display(self):
        """check that .display and Display calls are equivalent"""
        c = self.test_build_convolution()
        c.display() == c.ascii_write(False)
        s = str(c)
        assert c.display() == s
        assert c.display()==Display(c)

    def str(self):
        self.test_display()

    def test_ascii_write(self):
        self.test_display()

    def test_len(self):
        c = self.test_build_convolution()
        
        assert len(c) == 2
        assert len(c) == c.nb_distribution()

    def test_plot(self):
        """run plotting routines """
        c = self.test_build_convolution()
        if DISABLE_PLOT==False:
            c.plot()
 
    def test_save(self):
        c1 = self.test_build_convolution()
        c1.save('test1.dat')
        Save(c1, 'test2.dat')
        
        c1_read = Convolution('test1.dat')
        c2_read = Convolution('test2.dat')
        
        assert c1 and c1_read and c2_read
        assert len(c1)==len(c1_read) and len(c2_read)
        
        assert str(c1_read) == str(c2_read)
        
    def test_plot_write(self):
        c = self.test_build_convolution()
        c.plot_write('test', 'title')

    def test_file_ascii_write(self):
        d = self.test_build_convolution()
        d.file_ascii_write('test.dat', True)

    def test_spreadsheet_write(self):
        d = self.test_build_convolution()
        d.spreadsheet_write('test.dat')
            
    def test_simulation(self):
        """Test the simulate method"""

        c = self.test_build_convolution()
        s = c.simulate(1000)
        
        assert s.nb_histogram() == 2
        assert str(s)
        
        Simulate(c, 1000)

    def test_extract(self):
        """run and test the extract methods"""

        m = self.test_build_convolution()

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

        e = s.estimate_convolution(Binomial(0, 10, 0.5), Estimator="Parametric")

        d = e.extract_data()
        assert d
        
        eprime = Estimate(s, "CONVOLUTION", Binomial(0, 10, 0.5))
        # todo: check that e and eprime are similar





def test1():
    from openalea.stat_tool import ExtractHistogram, ExtractData
    from openalea.stat_tool import ToHistogram
    from openalea.stat_tool import Histogram, Shift, Display, Save, Plot

    _convol0 = Convolution(Distribution("B", 0, 10, 0.5), 
                           Distribution("NB", 0, 10, 0.5))
    convol1 = Convolution("convolution1.conv")

    convol_histo1 = Simulate(convol1, 200)

    _histo20 = ExtractHistogram(convol_histo1, "Elementary", 1)

    convol2 = Estimate(convol_histo1, "CONVOLUTION", 
                       ExtractDistribution(convol1, "Elementary", 1),
                       MinInfBound=0)

    _histo21 = ExtractHistogram(ExtractData(convol2), 'Elementary', 1)
    _histo22 = ToHistogram(ExtractDistribution(convol2, 'Elementary', 1))

    histo_b2 = Histogram("nothofagus_antarctica_bud_2.his")
    histo_s2 = Histogram("nothofagus_antarctica_shoot_2.his")

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
    if DISABLE_PLOT==False:
        Plot(convol31)
        Plot(ExtractDistribution(convol31, "Convolution"))
    Save(convol31, "nothofagus_antartica_2.xls", Format="SpreadSheet")




if __name__=="__main__":
    # perform all the test in the class Test (unit tests)
    test = TestUnit()
    for method in dir(test):
        if method.startswith('_'):
            continue
        if callable(getattr(test, method)):
            getattr(test, method)()
        else:
            print 'skipping'
    # and functional tests.
    test1()
