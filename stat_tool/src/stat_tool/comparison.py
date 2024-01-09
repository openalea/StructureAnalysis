#!/usr/bin/env python
#-*- coding: utf-8 -*-
"""Comparison

.. topic:: comparison.py summary

    A module dedicated to Comparison tests

    :Code status: mature
    :Documentation status: to be completed
    :Author: Thomas Cokelaer <Thomas.Cokelaer@sophia.inria.fr>

    :Revision: $Id: comparison.py 15183 2013-11-06 10:35:50Z jbdurand $
    
"""
__version__ = "$Id: comparison.py 15183 2013-11-06 10:35:50Z jbdurand $"

import error

from enums import variable_type
from enums import format_type


from openalea.stat_tool._stat_tool import \
    _CompoundData,\
    _DiscreteDistributionData,\
    _DiscreteMixtureData,\
    _ConvolutionData,\
    _Vectors,\
    _VectorDistance,\
    _FrequencyDistribution

__all__ = ['Compare', 'ComparisonTest']


def compare_histo(histo, *args, **kargs):
    """Comparison of frequency distributions.

    :Parameters:
      * `histo1`, `histo2`, ... (histogram, mixture_data, convolution_data, compound_data),
      * `type` (string): variable type ("NUMERIC" ("N"), "ORDINAL" ("O") or "SYMBOLIC" ("S")).

    :Keywords:
      - FileName (string) : name of the result file
      - Format (string) : format of the result file: "ASCII" (default format) or "SpreadSheet".
        This optional argument can only be used in conjunction with the optional argument FileName.

    :Returns:
      The comparison result.

    :Examples:

    .. doctest::
        :options: +SKIP

        >>> compare_histo(histo1, histo2, ..., type, FileName="result",
        ... Format="ASCII")

    .. seealso::
        :func:`~openalea.stat_tool.comparison.Compare`

    """
    utype = args[-1]
    if utype not in variable_type.keys():
        raise KeyError("%s not found. Allowed keys are %s"
                       % (utype, variable_type.keys()))


    utype = variable_type[args[-1]]

    error.CheckType([histo],
                        [[_DiscreteDistributionData, _DiscreteMixtureData,
                          _ConvolutionData, _CompoundData]])

    histos = args[0:-1]
    for h in histos:
        error.CheckType([h],
                        [[_DiscreteDistributionData, _DiscreteMixtureData,
                          _ConvolutionData, _CompoundData]])
    filename = kargs.get("Filename", None)
    format = error.ParseKargs(kargs, "Format", "ASCII",
                                  possible=format_type)

    ret = histo.compare(histos, utype, filename, format)

    return ret


_FrequencyDistribution.compare_histo = compare_histo	


def compare_vectors(vec, vector_distance, Standardization=True):
    """Comparison of vectors.

    The type _VectorDistance implements standardization procedures.
    The objective of standardization is to avoid the dependence on
    the variable type (chosen among symbolic, ordinal, numeric and circular)
    and, for numeric variables, on the choice of the measurement units
    by converting the original variables to dimensionless variables.

    :Parameters:
     - `vec` (_Vectors) : test
     - `vector_distance` (_VectorDistance) : test

    :Returns:
      An object of type _DistanceMatrix is returned.

    :Examples:

    .. doctest::
        :options: +SKIP

        >>> compare_vectors(vec, vector_distance)

    .. seealso::
        :func:`~openalea.stat_tool.vectors.VectorDistance`,
        :func:`~openalea.stat_tool.cluster.Clustering`,
        :func:`~openalea.stat_tool.comparison.Compare`
     """
    error.CheckType([vec, vector_distance], [_Vectors, _VectorDistance])
    error.CheckType([Standardization], [bool])

    return vec.compare(vector_distance, Standardization)


def Compare(arg1, *args, **kargs):
    """Comparison functions factory

    :Parameters:
      - `arg1` should be in :
          * `compare_histo` : Histograms comparison
          * `compare_vectors` : Vectors comparison
          * `compare_seq` : Sequences comparison
          * `compare_markov` : Markovian models comparison

    .. seealso::
        :func:`~openalea.stat_tool.comparison.compare_histo`,
        :func:`~openalea.stat_tool.comparison.compare_vectors`
        :func:`~openalea.stat_tool.comparison.compare_seq`
        :func:`~openalea.stat_tool.comparison.compare_markov`

    .. todo:: Get the AMAPMod documentation


    """

    p1 = arg1

    if isinstance(p1, _Vectors):
        ret = compare_vectors(arg1, *args, **kargs)
    elif isinstance(p1, _FrequencyDistribution):
        ret = compare_histo(arg1, *args, **kargs)
    else:
        raise NotImplementedError("First argument must be either Vectors or FrequencyDistribution")


    return ret


def ComparisonTest(utype, histo1, histo2):
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
    error.CheckType([histo1, histo2],
                        [[_DiscreteDistributionData, _DiscreteMixtureData,
                          _ConvolutionData, _CompoundData]]*2)

    utype = utype.lower()
    #todo: move this dict to enumerate.py ?
    type_dict = {
	    "f": "f_comparison",
	    "t": "t_comparison",
	    "w": "wmw_comparison",
	    }

    if not type_dict.has_key(utype):
        raise TypeError("to be done")

    func = getattr(histo1, type_dict[utype])
    ret = func(histo2)

    return ret



