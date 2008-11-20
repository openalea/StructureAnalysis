__doc__ = """Comparison functions"""
__docformat__ = "restructuredtext"



import _stat_tool
from distribution import Distribution
from histogram import Histogram



def compare_histo(histo, *args, **kargs):
    """    
    Comparison of frequency distributions.

    Usage
    -----
      * ``compare_histo(histo1, histo2, ..., type, FileName="result", Format="ASCII")``
      
    Parameters
    ----------
      * histo1, histo2,... (histogram, mixture_data, convolution_data, compound_data),
      * type (string): variable type ("NUMERIC" ("N"), "ORDINAL" ("O") or "SYMBOLIC" ("S")).
      
    Keywords
    --------
      * FileName (string): name of the result file,
      * Format (string): format of the result file: "ASCII" (default format) or "SpreadSheet". 
      This optional argument can only be used in conjunction with the optional argument FileName. 

    Return
    ------
      The comparison result.
    """
	
    type_dict = {
	    "O" : _stat_tool.ORDINAL,
	    "ORDINAL" : _stat_tool.ORDINAL,
	    "N" : _stat_tool.NUMERIC,
	    "NUMERIC" : _stat_tool.NUMERIC,
	    "S" : _stat_tool.SYMBOLIC,
	    "SYMBOLICL" : _stat_tool.SYMBOLIC,
	    }

    type = args[-1]
    try:
       type = type_dict[type.upper()]
    except KeyError:
       print "Type not reconized."
       return

    base = histo
    histos = args[0:-1]

    filename = kargs.get("Filename", None)
    format = kargs.get("Format", "a")

    ret = base.compare(histos, type, filename, format)
    return ret


_stat_tool._Histogram.compare_histo = compare_histo	    



def compare_vectors(vec, vector_distance):
    """Comparison of vectors.

    Usage
    -----
       * compare_vectors(vec, vector_distance) 
 
    Parameters
    ----------
      * vec (_Vectors),
      * vector_distance (_VectorDistance). 

    Return
    ------
     An object of type _DistanceMatrix is returned. 

    Background
    ----------
    The type _VectorDistance implements standardization procedures. 
    The objective of standardization is to avoid the dependence on 
    the variable type (chosen among symbolic, ordinal, numeric and circular) 
    and, for numeric variables, on the choice of the measurement units 
    by converting the original variables to unitless variables. 

    See Also
    --------
     `VectorDistance`, `Clustering`.  
     """

    return vec.compare(vector_distance)


def compare_seq(seq, *args, **kargs):
    raise NotImplementedError()


def compare_markov(mc, *args, **kargs):
    raise NotImplementedError()



def Compare(arg1, *args, **kargs):
    """ 
    Comparison function for

    For more details see :
      `compare_histo` : Histograms comparison
      `compare_vectors` : Vectors comparison
      `compare_seq` : Sequences comparison
      `compare_markov` : Markovian models comparison
    """

    p1 = arg1
    
    if(isinstance(p1, _stat_tool._Histogram)):
        return compare_histo(arg1, *args, **kargs)
    
    elif(isinstance(p1, _stat_tool._Vectors)):
        return compare_vectors(arg1, *args, **kargs)
    
    elif(isinstance(p1, _stat_tool._Sequence)):
        return compare_seq(arg1, *args, **kargs)

    elif(isinstance(p1, _stat_tool._Markovian)):
        return compare_markov(arg1, *args, **kargs)
    



def ComparisonTest(type, histo1, histo2):
    """
    Test of comparaison of frequency distributions.
    
    Usage
    -----
      * ``ComparisonTest(type, histo1, histo2)``

    Parameters
    ----------
       * type(string) : type of test "F" (Fisher-Snedecor), "T" (Student) 
         or "W" (Wilcoxon-Mann-Whitney)
       * histo1, histo2 (Histogram, MixtureData, ConvolutionData, CompoundData)

    Return
    ------
       A string containing the result of the tests

    Background
    ----------
       The objective is to compare two independent random samples in order to decide 
    if they have been drawn from the same population or not.
    In the case of samples from normal populations, the Fisher-Snedecor ("F") test 
    enables to test is the two variances are not significantly different. The normal 
    distribution assumption should be checked for instance by the exam of the shape 
    coefficients (skewness and kurtosis coefficients). The test statistic is:
    

    where  and  are the means of the samples.
    
       The Fisher-Snedecor variable with degrees of freedom and degrees of freedom can 
    be interpreted as the ratio of the variance estimators of the two samples. In practice,  
    the larger estimated variance is put at the denominator. Hence . The critical region is 
    of the form (one-sided test).
    
       In the case of samples from normal populations with equal variances, the Student ("T") 
    test enables to test if the two means are not significantly different. The test statistic is: 

    
       The critical region is of the form (two-sided test). For sufficiently large sample 
    sizes, this test of sample mean comparison can be used for samples from non-normal 
    populations with unequal variances. This test is said to be robust.
    """

    type = type.lower()
    type_dict = { 
	    "f": "f_comparison",
	    "t": "t_comparison",
	    "w": "wmw_comparison",
	    }

    if(not type_dict.has_key(type)):
	    raise TypeError()

    func = getattr(histo1, type_dict[type])
    ret = func(histo2)
    return ret



############# Tests ############################################################
from openalea.stat_tool import get_test_file

class Test:
    def test_comparisontest(self):

        meri1 = Histogram(get_test_file("meri1.his"))
        meri2 = Histogram(get_test_file("meri2.his"))

        assert ComparisonTest("F", meri1, meri2)
        assert ComparisonTest("T", meri1, meri2)
        assert ComparisonTest("W", meri1, meri2)


    def test_comparison(self):

        meri1 = Histogram(get_test_file("meri1.his"))
        meri2 = Histogram(get_test_file("meri2.his"))
        meri3 = Histogram(get_test_file("meri3.his"))

        assert Compare(meri1, meri2, meri3, 'N')
        assert Compare(meri1, meri2, meri3, 'O')
        assert Compare(meri1, meri2, meri3, 'S')


    def test_compare_vectors(self):
        
        from vectors import Vectors, VectorDistance
        from data_transform import SelectVariable

        vec10 = Vectors(get_test_file("chene_sessile.vec"))
        vec15 = SelectVariable(vec10, [1, 3, 6], Mode="Reject")
        assert vec15

        matrix10 = Compare(vec15, VectorDistance("N", "N", "N"))
        assert matrix10


    def test_compare_sequence(self):
        pass
        

    def test_compare_markov(self):
        pass


