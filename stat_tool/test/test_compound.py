"""Compound tests"""
__revision__ = "$Id: $"


from openalea.stat_tool.compound import Compound
from openalea.stat_tool.data_transform import ExtractDistribution
from openalea.stat_tool.distribution import Uniform


class Test:
    """a simple unittest class"""
    
    def test_empty(self):
        """Test that empty constructor fails"""
        try:
            m = Compound()
            assert False
        except TypeError:
            assert True

    def test_file(self):
        """run constructor with filename argument"""
        c = Compound("compound1.cd")
        assert c

    def test_build_compound(self):
        """run constructor with  two distributions as arguments"""
        
        d1 = Uniform(0, 10)
        d2 = Uniform(10, 20)
        
        m = Compound(d1, d2)
        
        assert m
        return m

    def test_plot(self):        
        """run plotting routines """
        m = self.test_build_compound()
        m.plot()

        assert str(m)
        m.display()

    def test_simulation(self):
        """Test the simulate method"""
        
        m = self.test_build_compound()
        s = m.simulate(1000)

        assert len(s) == 1000
        assert str(s)

    def test_extract(self):
        """run and test the extract methods"""
        m = self.test_build_compound()

        assert m.extract_compound() == ExtractDistribution(m, "Compound")

        assert m.extract_sum() == Uniform(0, 10)
        assert m.extract_sum() == ExtractDistribution(m, "Sum")

        assert m.extract_elementary() == Uniform(10, 20)
        assert m.extract_elementary() == ExtractDistribution(m, "Elementary")

    def test_extract_data(self):
        """todo : check if this test makes sense"""
        from openalea.stat_tool.distribution import Binomial

        m = self.test_build_compound()
        s = m.simulate(1000)

        m = s.estimate_compound(Binomial(0, 10, 0.5))

        d = m.extract_data()
        assert d
    
from openalea.stat_tool.estimate import Estimate
from openalea.stat_tool.simulate import Simulate
from openalea.stat_tool.data_transform import ExtractHistogram, ExtractData, Shift
from openalea.stat_tool.histogram import Histogram
from openalea.stat_tool.distribution import ToHistogram

def test1():
    
    cdist1 = Compound("compound1.cd")

    chisto1 = Simulate(cdist1, 200)

    histo30 = ExtractHistogram(chisto1, "Sum")

    #cdist2 = Estimate(chisto1, "COMPOUND", ExtractDistribution(cdist1, "Elementary"), "Sum", MinInfBound=0)

    cdist2 = Estimate(chisto1, "COMPOUND", ExtractDistribution(cdist1, "Elementary"), ExtractDistribution(cdist1, "Sum"), MinInfBound=0)
    
    histo31 = ExtractHistogram(ExtractData(cdist2), "Sum")
    histo32 = ToHistogram(ExtractDistribution(cdist2, "Sum"))
    
    peup1 = Histogram("peup1.his")
    mixt4 = Estimate(peup1, "MIXTURE", "B", "NB")
    histo33 = ToHistogram(ExtractDistribution(mixt4, "Component", 2))
    histo34 = Shift(histo33, -11)

    #cdist3 = Estimate(histo34, "COMPOUND", Distribution("B", 0, 1, 0.7), "Sum")
    #cdist4 = Estimate(histo34, "COMPOUND", Distribution("B", 0, 1, 0.7), "Sum", MinInfBound=0)



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
