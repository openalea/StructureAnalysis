"""Comparison functions"""
__revision__ = "$Id$"

import _stat_tool
#from distribution import Distribution
from histogram import Histogram


def compare_histo(histo, *args, **kargs):
    """    
    Comparison of frequency distributions.
  
    :Parameters:
      * `histo1`, `histo2`,... (histogram, mixture_data, convolution_data, compound_data),
      * type (string): variable type ("NUMERIC" ("N"), "ORDINAL" ("O") or "SYMBOLIC" ("S")).
      
    :Keywords:
      * FileName (string): name of the result file,
      * Format (string): format of the result file: "ASCII" (default format) or "SpreadSheet". 
        This optional argument can only be used in conjunction with the optional argument FileName. 

    :Returns:
        The comparison result.
      
    :Examples:
    .. doctest::
        :options: +SKIP
    
        >>> compare_histo(histo1, histo2, ..., type, FileName="result", 
        ... Format="ASCII")
          
    """
    type_dict = {
	    "O" : _stat_tool.ORDINAL,
	    "ORDINAL" : _stat_tool.ORDINAL,
	    "N" : _stat_tool.NUMERIC,
	    "NUMERIC" : _stat_tool.NUMERIC,
	    "S" : _stat_tool.SYMBOLIC,
	    "SYMBOLICL" : _stat_tool.SYMBOLIC,
	    }

    _type = args[-1]
    try:
        _type = type_dict[_type.upper()]
    except KeyError:
        print "Type not recognized."
        return

    base = histo
    histos = args[0:-1]

    filename = kargs.get("Filename", None)
    format = kargs.get("Format", "a")

    ret = base.compare(histos, _type, filename, format)
    return ret


_stat_tool._Histogram.compare_histo = compare_histo	    


def compare_vectors(vec, vector_distance):
    """Comparison of vectors.

    The type _VectorDistance implements standardization procedures. 
    The objective of standardization is to avoid the dependence on 
    the variable type (chosen among symbolic, ordinal, numeric and circular) 
    and, for numeric variables, on the choice of the measurement units 
    by converting the original variables to dimensionless variables. 
 
    :Parameters:
      - `vec` (_Vectors) 
      - `vector_distance` (_VectorDistance) 

    :Returns:
        An object of type _DistanceMatrix is returned. 

    :Examples:
    
    .. doctest::
        :options: +SKIP
    
        >>> compare_vectors(vec, vector_distance) 
    
    .. seealso::
        `VectorDistance`, `Clustering`.  
     """

    return vec.compare(vector_distance)


def compare_seq(seq, *args, **kargs):
    """not implemented"""
    raise NotImplementedError()


def compare_markov(mc, *args, **kargs):
    """not implemented"""
    raise NotImplementedError()



def Compare(arg1, *args, **kargs):
    """ 
    Comparison functions factory
    
    :Parameters:
      - `arg1` - should be in: 
          * `compare_histo` : Histograms comparison
          * `compare_vectors` : Vectors comparison
          * `compare_seq` : Sequences comparison
          * `compare_markov` : Markovian models comparison
          
    .. todo:: Get the AMAPMod documentation
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
    r"""
    Test of comparaison of frequency distributions.
    
    The objective is to compare two independent random samples in order to decide 
    if they have been drawn from the same population or not.
    In the case of samples from normal populations, the Fisher-Snedecor ("F") test 
    enables to test is the two variances are not significantly different. The normal 
    distribution assumption should be checked for instance by the exam of the shape 
    coefficients (skewness and kurtosis coefficients). The test statistic is:
    
    .. math::
        F_{n_1-1,n_2-1} = \frac
            {
            \frac{\displaystyle\sum_{i=1}^{n_1}\left( x_{1i}-m_1 \right)^2}{n_1-1}
            }
            {
            \frac{\displaystyle\sum_{i=1}^{n_2}\left( x_{2i}-m_2 \right)^2}{n_2-1}
            }
    
    where :math:`m_1` and :math:`m_2` are the means of the samples.
    
    The Fisher-Snedecor variable :math:`F_{n_1-1,n_2-1}` with :math:`n_1-1` degrees
    of freedom and :math:`n_2-1` degrees of freedom can 
    be interpreted as the ratio of the variance estimators of the two samples.
    In practice, the larger estimated variance is put at the denominator. Hence 
    :math:`F_{n_1-1,n_2-1} \geq 1` . The critical region is of the form
    :math:`F_{n_1-1,n_2-1} > f` (one-sided test).
    
    In the case of samples from normal populations with equal variances,
    the Student ("T") test enables to test if the two means are not significantly
    different. The test statistic is:
    
    .. math::
        T_{n_1+n_2 - 2} = \frac{m_1 - m_2}{
        \sqrt{\left(
            \displaystyle\sum_{i=1}^{n_1}\left( x_{1i}-m_1 \right)^2{n_1-1}
            +
            \displaystyle\sum_{i=1}^{n_2}\left( x_{2i}-m_1 \right)^2{n_2-1}
            \right)
            \left( \frac{1}{n_1} + \frac{1}{n_2}\right)
        }        
        } \sqrt{n_1 + n_2 - 2}

    The critical region is of the form :math:`\left| T_{n_1+n_2-2}\right| > t` 
    (two-sided test). For sufficiently large sample 
    sizes, this test of sample mean comparison can be used for samples from non-normal 
    populations with unequal variances. This test is said to be robust.
    
    The Wilcoxon-Mann-Whitney ("W") test is a distribution-free test relying on
    the homogeneity of the ranking of the two sample (ranks of one sample should
    not cluster at either or both ends of the range). It can be seen as the
    non-parametric analog of the Student's t test and can be applied to compare
    two sets of observations measures on an interval scale when it is supposed
    that the data are non-normally distributed, or to compare two sets of
    observations measured on an ordinal scale.

    :Parameters:
       * type(string) : type of test "F" (Fisher-Snedecor), "T" (Student) 
         or "W" (Wilcoxon-Mann-Whitney)
       * histo1, histo2 (Histogram, MixtureData, ConvolutionData, CompoundData)

    :Returns:
       A string containing the result of the tests
       
    :Examples:

    .. doctest::
        :options: +SKIP

        >>> ComparisonTest(type, histo1, histo2)
    
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


    
