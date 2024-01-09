"""test data_transform methods and compare them with class members methods

:Author: Thomas Cokelaer, Thomas.Cokelaer@inria.fr

.. todo:: to be clean
"""
__version__ = "$Id: test_data_transform.py 9397 2010-08-10 10:50:07Z cokelaer $"

from openalea.stat_tool.vectors import Vectors
from openalea.stat_tool.histogram import Histogram
from openalea.stat_tool.distribution import Distribution, Uniform, Binomial, \
    NegativeBinomial
from openalea.stat_tool.mixture import Mixture
from openalea.stat_tool.convolution import Convolution
from openalea.stat_tool.simulate import Simulate
from openalea.stat_tool.compound import Compound
from openalea.stat_tool.output import Plot
from openalea.stat_tool.data_transform import Shift, Merge, Fit, ValueSelect, \
    SelectVariable, SelectIndividual, MergeVariable, ExtractDistribution, \
    ExtractHistogram, ExtractData, SelectStep

from tools import runTestClass

from openalea.stat_tool import get_shared_data

class data():
    """Create some common data sets for the classes
    shift, fit, merge, ...
    """

    def __init__(self):
        self.N = 1000

    def comp_data(self):
        d1 = Binomial(0, 10, 0.5)
        d2 = NegativeBinomial(0, 1, 0.1)
        comp = Compound(d1, d2)
        comp_data = comp.simulate(self.N)
        return comp_data

    def conv_data(self):
        d1 = Binomial(0, 10, 0.5)
        d2 = NegativeBinomial(0, 1, 0.1)
        conv = Convolution(d1, d2)
        conv_data = conv.simulate(self.N)
        return conv_data

    def mixt_data(self):
        d1 = Binomial(0, 10, 0.5)
        d2 = NegativeBinomial(0, 1, 0.1)
        mixt = Mixture(0.1, d1, 0.4, d2)
        mixt_data = mixt.simulate(self.N)
        return mixt_data

    def hist_data(self):
        h = Histogram(get_shared_data( "meri2.his"))
        return h

    def int_vector_data(self):
        a = [[1, 3, 4],
             [4, 12, 2],
             [8, 7, 3],]
        v = Vectors(a)
        return v

    def float_vector_data(self):
        v = Vectors([[0.1, 0.3, 4.2],
                     [0.5, 2.3, 1.2],
                     [4.5, 6.3, 3.2],
                     ])
        return v

    def mixt(self):
        d1 = Uniform(0, 10)
        d2 = Uniform(10, 20)
        d3 = Uniform(20, 30)
        m = Mixture(0.1, d1, 0.2, d2, 0.7, d3)
        return m

    def comp(self):
        d1 = Binomial(2, 5, 0.5)
        d2 = NegativeBinomial(0, 2, 0.5)
        c = Compound(d1, d2)
        return c

    def conv(self):
        d1 = Binomial(0, 10, 0.5)
        d2 = NegativeBinomial(0, 1, 0.1)
        m = Convolution(d1, d2)
        return m

class TestSelectStep(data):

    def __init__(self):
        data.__init__(self)

    def test_select_step_vectors(self):
        data = self.int_vector_data
        try:
            SelectStep(data,1, 100)
            assert False
        except:
            assert True


class TestShift(data):

    def __init__(self):
        data.__init__(self)

    def test_shift_compound_data(self):
        comp_data = self.comp_data()
        assert comp_data.shift(20) == Shift(comp_data, 20)

    def test_shift_convolution_data(self):
        conv_data = self.conv_data()
        assert conv_data.shift(20) == Shift(conv_data, 20)

    def test_shift_mixture_data(self):
        mixt_data = self.mixt_data()
        assert mixt_data.shift(20) == Shift(mixt_data, 20)

    def test_shift_histo(self):
        # equivalent to test Distribution data
        h = data.hist_data(self)
        assert h.shift(2)
        assert Shift(h, 2) == h.shift(2)

    def test_shift_vector(self):
        vn = Vectors([[0.2, 1., 2, 3],
                      [4.2, 5, 6, 7]])
        v1 = Vectors([[1.]])
        assert Shift(v1, 2)
        assert Shift(vn, 1, 2)
        assert str(Shift(vn, 1, 2)) == str(vn.shift(1, 2))


class TestFit:

    def __init__(self):
        pass

    def test_fit_histogram(self):
        meri5 = Histogram(get_shared_data( "meri5.his"))
        dist1 = Fit(meri5, Distribution("B", 0, 10, 0.437879))
        dist2 = meri5.fit(Distribution("B", 0, 10, 0.437879))
        assert str(dist1)==str(dist2)


class TestSelectHist:

    def __init__(self):
        pass
    def test_value_select_float(self):
        meri1 = Histogram(get_shared_data("meri1.his"))
        # note keep=False is equivalent to Mode=keep,is this correct ?
        assert str(ValueSelect(meri1, 0, 10, Mode="Keep"))==\
            str(meri1.value_select( min=0, max=10, keep=True))


class TestValueSelect(data):

    def __init__(self):
        data.__init__(self)

    def test_value_select_float(self):

        v = self.float_vector_data()

        v1b = ValueSelect(v, 1, 0.2, 2.0, Mode="Keep")
        v2 =  ValueSelect(v, 2, 1.0, 6.0, Mode="Keep")
        v3 =  ValueSelect(v, 3, 1.0, 2.0, Mode="Keep")

        assert v and v1b and v2 and v3
        assert len(v1b) == 1

        assert str(ValueSelect(v, 1, 0.2, 2, Mode="Keep")) == \
                str(v.value_select(1, 0.2, 2, keep=True))
        assert str(ValueSelect(v, 2, 1.0, 6.0, Mode="Keep")) == \
                str(v.value_select(2, 1.0, 6.0, keep=True))
        assert str(ValueSelect(v, 3, 1.0, 2.0, Mode="Keep")) == \
                str(v.value_select(3, 1.0, 2., keep=True))

    def test_value_select_int(self):

        v = self.int_vector_data()
        v1b =  ValueSelect(v, 1, 2, 8, Mode="Reject")
        assert v1b

        try:
            _v1 =  ValueSelect(v, 1, 2, 3, Mode="Keep")
            assert False
        except Exception:
            # Empty sample
            assert True

    def test_compound_data(self):
        comp_data = self.comp_data()
        assert comp_data.value_select(1, 20, True)

    def test_convolution_data(self):
        conv_data = self.conv_data()
        assert conv_data.value_select(1, 20, True)

    def test_mixture_data(self):
        mixt_data = self.mixt_data()
        assert mixt_data.value_select(1, 20, True)

    def test_histo_data(self):
        hist_data = self.hist_data()
        assert hist_data.value_select(1, 20, True)

class TestSelectVariable(data):

    def __init__(self):
        data.__init__(self)

    def test_vector(self):

        v = self.int_vector_data()

        for i in range(3):
            v1 =  SelectVariable(v, i+1, Mode="Keep")
            assert v1
            assert len(v1) == 3
            assert len(v1[0]) == 1
            assert str(v.select_variable([i+1], keep=True)) == \
                str(SelectVariable(v, i+1, Mode="Keep"))
            for j in range(3):
                assert v1[j][0] == v[j][i]


class TestSelectIndividual(data):

    def __init__(self):
        data.__init__(self)

    def test_vector_integer(self):

        v = self.int_vector_data()

        selection = SelectIndividual(v, [1, 2], Mode="Keep")
        selection2 = v.select_individual([1, 2], keep=True)

        assert str(selection)==str(selection2)

        assert len(selection) == 2

    def test_distance_matrix(self):
        """not implemented - distance matrix """
        pass


class TestExtract():
    """
    Extract("Compound") see test_compound
    Extract("Convolution") see test_convolution
    Extract("Mixture") see test_mixture
    Sum, elementary, component, weight see test_compound, and so on.
    """
    def __init__(self):
        pass



class TestExtractData():
    """
    See other test file.
    """

    def __init__(self):
        pass

    def test_histo_extract_data(self):

        h = Histogram(get_shared_data("meri2.his"))
        mixt = h.estimate_mixture(["B", "NB"])
        assert ExtractData(mixt)
        assert mixt.extract_data() == ExtractData(mixt)


class TestExtractDistribution(data):

    def __init__(self):
        data.__init__(self)

    def test_mixture(self):
        mixt = self.mixt()

        assert ExtractDistribution(mixt, "Weight")
        assert ExtractDistribution(mixt, "Weight") == mixt.extract_weight()
        assert ExtractDistribution(mixt, "Mixture") == mixt.extract_mixture()
        assert ExtractDistribution(mixt, "Component", 1) == \
            mixt.extract_component(1)
        assert ExtractDistribution(mixt, "Component", 2) == \
            mixt.extract_component(2)
        assert ExtractDistribution(mixt, "Component", 3) == \
            mixt.extract_component(3)

        try:
            ExtractDistribution(mixt, "Component", 0)
            assert False
        except: # Bas distrubition index
            assert True

    def test_convolution(self):
        convol = self.conv()

        assert ExtractDistribution(convol, "Convolution") == \
            convol.extract_convolution()
        assert ExtractDistribution(convol, "Elementary", 1) == \
            convol.extract_elementary(1)
        assert ExtractDistribution(convol, "Elementary", 2) == \
            convol.extract_elementary(2)

    def test_compound(self):
        comp = self.comp()

        assert ExtractDistribution(comp, "Compound") == \
            comp.extract_compound()
        assert ExtractDistribution(comp, "Elementary") == \
            comp.extract_elementary()
        assert ExtractDistribution(comp, "Sum") == \
            comp.extract_sum()



class TestExtractHistogram:

    def __init__(self):
        pass

    def test_mixture(self):

        h = Histogram(get_shared_data("meri2.his"))
        mixt = h.estimate_mixture(["B", "NB"])

        assert ExtractHistogram(mixt, "Weight") == \
            mixt.extract_weight()
        assert ExtractHistogram(mixt, "Mixture") == \
            mixt.extract_mixture()
        assert ExtractHistogram(mixt, "Component", 1) == \
            mixt.extract_component(1)
        assert ExtractHistogram(mixt, "Component", 2) == \
            mixt.extract_component(2)

        try:
            ExtractHistogram(mixt, "Component", 3)
            assert False
        except: # Bas distrubition index
            assert True

    def test_convolution(self):
        convol = Convolution("data/convolution1.conv")
        convol_histo = Simulate(convol, 200)
        _histo = ExtractHistogram(convol_histo, "Elementary", 1)
        _histo = ExtractHistogram(convol_histo, "Convolution")

    def test_compound(self):
        comp  = Compound("data/compound1.cd")
        comp_histo = Simulate(comp, 200)
        _histo = ExtractHistogram(comp_histo, "Sum")
        _histo = ExtractHistogram(comp_histo, "Elementary")

    def test_vectors(self):
        v = Vectors([[1, 2, 3, 4, 5, 6, 7]])
        ExtractHistogram(v, 1)
        v = Vectors([[1, 2], [3, 4]])
        ExtractHistogram(v, 1)


class TestMerge(data):

    def __init__(self):
        data.__init__(self)

    def test_merge(self):

        mixt1 = Mixture(0.6, Distribution("B", 2, 18, 0.5),
                        0.4, Distribution("NB", 10, 10, 0.5))

        mixt_histo1 = Simulate(mixt1, 200)

        histo10 = mixt_histo1.extract_component(1)
        histo11 = mixt_histo1.extract_component(2)

        histo12 = Merge(histo10, histo11)

        assert histo12
        Plot(histo12)

    def test_merge_histo(self):
        meri1 = Histogram(get_shared_data( "meri1.his"))
        meri2 = Histogram(get_shared_data( "meri2.his"))
        meri3 = Histogram(get_shared_data( "meri3.his"))
        meri4 = Histogram(get_shared_data( "meri4.his"))
        meri5 = Histogram(get_shared_data( "meri5.his"))
       

        meri = Merge(meri1, meri2, meri3, meri4, meri5)
        assert meri
        meri_bis = meri1.merge([meri2, meri3, meri4, meri5])
        assert meri_bis
        assert str(meri)==str(meri_bis)
        Plot(meri)

    def test_merge_vectors(self):

        v1 = self.int_vector_data()
        b = [[2, 78, 45],
             [6, 2, 122],
             [3, 4, 31],]
        v2 = Vectors(b)

        v = Merge(v1, v2)
        assert v

        a = v1.merge([v2])
        b = v2.merge([v1])
        assert str(a)==str(b)

        assert str(a)==str(v)
        Plot(v)


class TestMergeVariable(data):

    def __init__(self):
        data.__init__(self)

    def test_vector(self):
        v = self.int_vector_data()

        v1 =  SelectVariable(v, 1, Mode="keep")
        v2 =  SelectVariable(v, 2, Mode="Keep")
        v3 =  SelectVariable(v, 3, Mode="Keep")

        merged = MergeVariable(v1, v2, v3)
        merged2 = v1.merge_variable([v2, v3], 1)
        for i in range(3):
            for j in range(3):
                assert merged[i][j] == v[i][j]

        assert str(merged) == str(merged2)



if __name__ == "__main__":
    runTestClass(TestSelectHist())
    runTestClass(TestFit())
    runTestClass(TestShift())
    runTestClass(TestExtractHistogram())
    runTestClass(TestMerge())
    runTestClass(TestMergeVariable())
    runTestClass(TestValueSelect())
    runTestClass(TestSelectVariable())
    runTestClass(TestSelectIndividual())
    runTestClass(TestExtract())
    runTestClass(TestExtractData())
    runTestClass(TestExtractDistribution())
