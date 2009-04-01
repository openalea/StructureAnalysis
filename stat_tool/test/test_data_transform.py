"""test data_transform"""
__revision__ = "$Id: $"

from openalea.stat_tool import vectors, histogram, distribution, mixture, convolution, simulate
from openalea.stat_tool.vectors import Vectors
from openalea.stat_tool.histogram import Histogram
from openalea.stat_tool.distribution import Distribution, Uniform, Binomial, NegativeBinomial
from openalea.stat_tool.mixture import Mixture
from openalea.stat_tool.convolution import Convolution
from openalea.stat_tool.simulate import Simulate
from openalea.stat_tool.compound import Compound

from openalea.stat_tool.data_transform import Shift, Merge, Fit, ValueSelect, SelectVariable, SelectIndividual, MergeVariable, ExtractDistribution, ExtractHistogram, ExtractData

class TestShift:

    def test_shift_compound_data(self):
        
        d1 = Binomial(0, 10, 0.5)
        d2 = NegativeBinomial(0, 1, 0.1)
        comp = Compound(d1, d2)
        comp_data = comp.simulate(1000)
        assert comp_data.shift(20) == Shift(comp_data, 20)
       
    def test_shift_convolution_data(self):
        
        d1 = Binomial(0, 10, 0.5)
        d2 = NegativeBinomial(0, 1, 0.1)
        conv = Convolution(d1, d2)
        conv_data = conv.simulate(1000)
        assert conv_data.shift(20) == Shift(conv_data, 20) 
    
    def test_shift_mixture_data(self):
        
        d1 = Binomial(0, 10, 0.5)
        d2 = NegativeBinomial(0, 1, 0.1)
        mixt = Mixture(0.1, d1, 0.4, d2)
        mixt_data = mixt.simulate(1000)
        assert mixt_data.shift(20) == Shift(mixt_data, 20)
         
    def test_shift_histo(self):
        # equivalent to test Distribution data
        h = Histogram("meri2.his")
        assert h.shift(2)

        assert Shift(h, 2) == h.shift(2)

    def test_shift_vector(self):

        vn = Vectors([[0.2, 1., 2, 3], [4.2, 5, 6, 7]])
        v1 = Vectors([[1.]])

        assert Shift(v1, 2)
        assert Shift(vn, 1, 2)
        
        assert str(Shift(vn, 1,2)) == str(vn.shift(1,2))
        
    def test_shift_sequence(self):
        pass
        #raise NotImplementedError()

class TestFit:
    def test_fit(self):

        meri5 = Histogram("meri5.his")
        dist5 = Fit(meri5, Distribution("B", 0, 10, 0.437879))
        assert dist5


class TestSelect:
    def test_value_select_float(self):
    
        v = Vectors([[0.1, 0.3, 4.2], 
                     [0.5, 2.3, 1.2],
                     [4.5, 6.3, 3.2],
                     ])

        v1b = ValueSelect(v, 1, 0.2, 2.0, Mode="keep")
        v2 =  ValueSelect(v, 2, 1.0, 6.0, Mode="keep")
        v3 =  ValueSelect(v, 3, 1.0, 2.0, Mode="keep")
        
        assert v and v1b and v2 and v3
        print len(v1b)
        assert len(v1b) == 2
        
    def test_value_select_int(self):

        v = Vectors([[1, 3, 4], 
                     [4, 12, 2],
                     [8, 7, 3],
                     ])

        v1b =  ValueSelect(v, 1, 2, 8, Mode="Reject")

        assert v1b

        try:
            v1 =  ValueSelect(v, 1, 2, 3, Mode="Keep")
            assert False
        except Exception:
            # Empty sample
            assert True
    
    def test_select_variable(self):

        a = [[1, 3, 4], 
             [4, 12, 2],
             [8, 7, 3],]

        v = Vectors(a)

        for i in range(3):
            v1 =  SelectVariable(v, i+1, Mode="Keep")
            assert v1
            assert len(v1) == 3
            assert len(v1[0]) == 1

            for j in range(3):
                assert v1[j][0] == a[j][i]


    def test_select_individual(self):
        
        a = [[1, 3, 4], 
             [4, 12, 2],
             [8, 7, 3],]

        v = Vectors(a)
        selection = SelectIndividual(v, [1,2], Mode="Keep")
        
        assert len(selection) == 2
        print selection[0] == [1, 3, 4]


class TestExtract:

    def test_build_mixture(self):


        d1 = Uniform(0, 10)
        d2 = Uniform(10, 20)
        d3 = Uniform(20, 30)

        m = Mixture(0.1, d1, 0.2, d2, 0.7, d3)
        assert m
        return m


    def test_build_convolution(self):

        d1 = Binomial(0, 10, 0.5)
        d2 = NegativeBinomial(0, 1, 0.1)

        m = Convolution(d1, d2)
        assert m
        return m


    def test_extract_data(self):

        h = Histogram("meri2.his")
        mixt = h.estimate_mixture(["B", "NB"])

        assert ExtractData(mixt)
        #assert Convolution().extract_data()
        #assert Compound().extract_data()


    def test_extract_mixture_distribution(self):

        mixt = self.test_build_mixture()

        assert ExtractDistribution(mixt, "Weight")
        assert ExtractDistribution(mixt, "Mixture")
        assert ExtractDistribution(mixt, "Component", 1)
        assert ExtractDistribution(mixt, "Component", 2)
        assert ExtractDistribution(mixt, "Component", 3)

        try:
            ExtractDistribution(mixt, "Component", 0)
            assert False
        except: # Bas distrubition index
            assert True


    def test_extract_convolution_distribution(self):

        convol = self.test_build_convolution()

        assert ExtractDistribution(convol, "Convolution")
        assert ExtractDistribution(convol, "Elementary", 1)
        assert ExtractDistribution(convol, "Elementary", 2)


    def test_extract_histogram(self):

        h = Histogram("meri2.his")
        mixt = h.estimate_mixture(["B", "NB"])

        assert ExtractHistogram(mixt, "Weight")
        assert ExtractHistogram(mixt, "Mixture")
        assert ExtractHistogram(mixt, "Component", 1)
        assert ExtractHistogram(mixt, "Component", 2)

        try:
            ExtractHistogram(mixt, "Component", 3)
            assert False
        except: # Bas distrubition index
            assert True


class TestMerge:
    
    def test_merge(self):
        
        mixt1 = Mixture(0.6, Distribution("B", 2, 18, 0.5), 
                        0.4, Distribution("NB", 10, 10, 0.5))

        mixt_histo1 = Simulate(mixt1, 200)

        histo10 = mixt_histo1.extract_component(1)
        histo11 = mixt_histo1.extract_component(2)

        histo12 = Merge(histo10, histo11)

        assert histo12

    def test_merge_histo(self):
        meri1 = Histogram("meri1.his")
        meri2 = Histogram("meri2.his")
        meri3 = Histogram("meri3.his")
        meri4 = Histogram("meri4.his")
        meri5 = Histogram("meri5.his")

        meri = Merge(meri1, meri2, meri3, meri4, meri5)
        assert meri

    def test_merge_vectors(self):
        
        a = [[1, 3, 4], 
             [4, 12, 2],
             [8, 7, 3],]

        b = [[2, 78, 45], 
             [6, 2, 122],
             [3, 4, 31],]


        v1 = Vectors(a)
        v2 = Vectors(b)

        v = Merge(v1, v2)
        assert v

    def test_mergevariable(self):
        a = [[1, 3, 4], 
             [4, 12, 2],
             [8, 7, 3],]

        v = Vectors(a)

        v1 =  SelectVariable(v, 1, Mode="Keep")
        v2 =  SelectVariable(v, 2, Mode="Keep")
        v3 =  SelectVariable(v, 3, Mode="Keep")

        merged = MergeVariable(v1,v2,v3)

        for i in range(3):
            for j in range(3):
                assert merged[i][j] == v[i][j]
        

    
if __name__=="__main__":
    # perform all the test in the class Test (unit tests)
    testshift = TestShift()
    testmerge = TestMerge()
    testextract = TestExtract()
    testselect = TestSelect()
    testfit = TestFit()
    
    for test in [testshift, testmerge, testextract, testselect]:
        for method in dir(test):
            if method.startswith('_'):
                continue
            if callable(getattr(test, method)):
                print 'testing %s ' % method
                getattr(test, method)()
            else:
                print 'skipping'
                    
    # and functional tests.    




