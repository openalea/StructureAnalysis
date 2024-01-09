#!/usr/bin/env python
#-*- coding: utf-8 -*-
"""

.. topic:: correlation.py summary

    A module dedicated to Correlation functions

    :Code status: mature
    :Documentation status: to be completed
    :Author: Thomas Cokelaer <Thomas.Cokelaer@sophia.inria.fr>

    :Revision: $Id: correlation.py 15827 2014-02-19 16:17:45Z jbdurand $
"""
__version__ = "$Id: correlation.py 15827 2014-02-19 16:17:45Z jbdurand $"


import openalea.stat_tool.interface as interface
import openalea.stat_tool.error as error
from openalea.stat_tool._stat_tool import I_DEFAULT
from openalea.stat_tool._stat_tool import _DiscreteParametricModel, _DiscreteMixture, \
    _Convolution, _Compound

# from openalea.stat_tool.mixture import Mixture
    
from openalea.sequence_analysis._sequence_analysis import _Correlation
from openalea.sequence_analysis._sequence_analysis import _NonHomogeneousMarkovData
from openalea.sequence_analysis._sequence_analysis import _HiddenVariableOrderMarkov
from openalea.sequence_analysis._sequence_analysis import _VariableOrderMarkovData
from openalea.sequence_analysis._sequence_analysis import _SemiMarkovData
from openalea.sequence_analysis._sequence_analysis import _Sequences
from openalea.sequence_analysis._sequence_analysis import _MarkovianSequences
from openalea.sequence_analysis._sequence_analysis import _VariableOrderMarkov


from openalea.sequence_analysis.enums import type_dict, norm_type
from openalea.sequence_analysis.data_transform import _check_nb_variable
from openalea.stat_tool import error

from openalea.sequence_analysis._sequence_analysis import MAX_LAG

__all__ = ['_Correlation',
           'ComputeCorrelation',
           'ComputeAutoCorrelation',
           'ComputeWhiteNoiseCorrelation',
           'ComputePartialAutoCorrelation']


# Extend dynamically class
interface.extend_class( _Correlation, interface.StatInterface)


def ComputeCorrelation(obj, *args, **kargs):
    """Computation of sample autocorrelation or cross-correlation functions.

    :Examples:

    .. doctest::
        :options: +SKIP
        
        >>> ComputeCorrelation(seq1, MaxLag=10, Type="Spearman", Normalization="Exact")
        >>> ComputeCorrelation(seqn, variable, MaxLag=10, Type="Spearman", Normalization="Exact")
        >>> ComputeCorrelation(seqn, variable1, variable2, MaxLag=10, Type="Spearman", Normalization="Exact")

    :Arguments:

    * seq1 (sequences, discrete_sequences, markov_data, semi-markov_data):
      univariate sequences,
    * seqn (sequences, discrete_sequences, markov_data, semi-markov_data):
      multivariate sequences,
    * variable (int): variable index (computation of a sample autocorrelation function).
    * variable1, variable2 (int): variable indices (computation of a sample
      cross-correlation function).

    :Optional Arguments:

    * Type (string): type of correlation coefficient: "Pearson" (linear
      correlation coefficient - default value), "Spearman" or "Kendall" (rank
      correlation coefficients).
    * MaxLag (int): maximum lag. A default value is computed from the sequence
      length distribution,
    * Normalization (STRING): normalization of the correlation coefficients:
      "Approximated" (the default - usual convention for time series analysis)
      or "Exact", (highly recommended for sample of short sequences). This
      optional argument can only be used if the optional argument Type is set
      at "Pearson" or "Spearman".

    :Returned Object:

    If variable, or variable1 and variable2 are valid indices of variables (and
    are different if two indices are given) and if 0 <= MaxLag < (maximum length
    of sequences), then an object of type correlation is returned, otherwise no
    object is returned.

    :Background:

    In the univariate case or if only variable is given, a sample
    autocorrelation function is computed. If variable1 and variable2 are given,
    a sample cross-correlation function is computed.

    .. seealso::

        :func:`~openalea.sequence_analysis.correlation.ComputePartialAutoCorrelation`,
        :func:`~openalea.sequence_analysis.correlation.ComputeWhiteNoiseCorrelation`

"""

    error.CheckType([obj], [[_Sequences, _MarkovianSequences,
                            _VariableOrderMarkovData, _SemiMarkovData,
                            _NonHomogeneousMarkovData]])

    if obj.nb_variable == 1:
        variable1 = 1
        variable2 = 1
    else:
        error.CheckType([args[0]], [int])
        #todo: check that variable1 <= nb_variable and > 0
        variable1 = args[0]
        if len(args) == 1:
            variable2 = variable1
        elif len(args) == 2:
            #todo: check that variable1 <= nb_variable and > 0
            error.CheckType([args[1]], [int])
            variable2 = args[1]
        else:
            raise TypeError("1 or 2  non-optional arguments required")


    max_lag = error.ParseKargs(kargs, "MaxLag", I_DEFAULT)
    itype = error.ParseKargs(kargs, "Type", "Pearson", type_dict)
    normalization = error.ParseKargs(kargs, "Normalization", "Exact", norm_type)
    IndividualMean = error.ParseKargs(kargs, "IndividualMean", False)


    #if normalization_option and ((type == SPEARMAN2) or (type == KENDALL)):
    #    raise Exception

    #if individual_mean_option and (type != PEARSON):
    #    raise Exception

    # check argument validity.
    return obj.correlation_computation(variable1, variable2,
                                   itype, max_lag, normalization,
                                   IndividualMean)



def ComputeAutoCorrelation(obj, *args, **kargs):
    """
    ComputeAutoCorrelation
    """
    error.CheckType([obj], [[_VariableOrderMarkov,
                             _HiddenVariableOrderMarkov,
                             _VariableOrderMarkovData]])

    error.CheckArgumentsLength(args, 1, 2)

    if len(args) == 1:
        variable = 1
        value = args[0]
    elif len(args) == 2:
        variable = args[0]
        value = args[1]

    #_check_nb_variable(obj, variable)
    max_lag = error.ParseKargs(kargs, "MaxLag", MAX_LAG)

    error.CheckType([variable, value, max_lag], [int, int, int])

    if len(args) == 1:
        return obj.state_autocorrelation_computation(value, max_lag)
    elif len(args)==2:
        return obj.output_autocorrelation_computation(variable, value, max_lag)



def ComputeWhiteNoiseCorrelation(obj, itype):
    """ComputeWhiteNoiseCorrelation

    Computation of the sample autocorrelation or cross-correlation function induced on a white noise sequence by filtering.

    :Usage:

    .. doctest::
        :options: +SKIP
        
        >>> ComputeWhiteNoiseAutoCorrelation(cf, order)
        >>> ComputeWhiteNoiseAutoCorrelation(cf, filter)
        >>> ComputeWhiteNoiseAutoCorrelation(cf, frequencies)
        >>> ComputeWhiteNoiseAutoCorrelation(cf, dist)

    :Arguments:

    * cf (correlation): sample autocorrelation or cross-correlation function (in the Pearsons sense),
    * order (int): order of differencing,
    * filter (array(real)): filter values on a half width i.e. from one extremity to the central value (with the constraint that filteri + filterm = 1),
    * frequencies (array(int)): frequencies defining the filter,
    * dist (distribution, mixture, convolution, compound): symmetric distribution whose size of the support is even defining the filter (for instance Distribution("BINOMIAL", 0, 4, 0.5)),

    :Returned Object:

    No object is returned.

    :Background:

    The application of linear filters for trend removal induces an autocorrelation structure. The effect of a given linear filter on the autocorrelation structure of the residual sequence can be roughly described as follows: the number of non-zero induced autocorrelation coefficients increase with the width of the filter while their numerical magnitudes decrease.

    .. seealso::  :func:`~openalea.sequence_analysis.correlation.ComputeCorrelation`.
    """

    error.CheckType([obj], [_Correlation])
    error.CheckType([itype], [[int, list, _DiscreteParametricModel,
                              _DiscreteMixture, _Convolution, _Compound ]])

    if isinstance(itype, int):
        obj.white_noise_correlation_order(itype)
    elif isinstance(itype, list):
        obj.white_noise_correlation_filter(itype)
    else:
        obj.white_noise_correlation_dist(itype)

    return obj


def ComputePartialAutoCorrelation(obj, *args, **kargs):
    """ComputePartialAutoCorrelation

    Computation of sample partial autocorrelation functions.

    :Usage:

    .. doctest::
        :options: +SKIP
        
        >>> ComputePartialAutoCorrelation(seq1, MaxLag=10, Type="Kendall")
        >>> ComputePartialAutoCorrelation(seqn, variable, MaxLag=10, Type="Kendall")

    :Arguments:

    * seq1 (**sequences**, discrete_sequences, markov_data, semi-markov_data):
      univariate sequences,
    * seqn (sequences, discrete_sequences, markov_data, semi-markov_data):
      multivariate sequences,
    * variable (int): variable index.

    :Optional Arguments:

    * MaxLag (int): maximum lag. A default value is computed from the sequence
      length distribution,
    * Type (string): type of correlation coefficient: "Pearson" (linear
      correlation coefficient - the default) or "Kendall" (rank correlation coefficient).

    :Returned Object:

    If variable is a valid variable index and if 1 <= MaxLag < (maximum length
    of sequences), an object of type correlation is returned, otherwise no object
    is returned.

    :Background:

    The partial autocorrelation coefficient at lag k measures the correlation
    between :math:`x_i` and :math:`x_{t+k}` not accounted for by
    :math:`x_{t+1}, ..., x_{t+k-1}`  (or after adjusting for the effects of
    :math:`x_{t+1}, ..., x_{t+k-1}`).

    .. seealso::

        :func:`~openalea.sequence_analysis.correlation.ComputeCorrelation`
    """
    error.CheckType([obj], [[_Sequences, _MarkovianSequences,
                             _VariableOrderMarkovData, _SemiMarkovData,
                             _NonHomogeneousMarkovData]])

    error.CheckArgumentsLength(args, 0, 1)


    if len(args) == 0:
        variable = 1
    else:
        variable = args[0]

    max_lag = error.ParseKargs(kargs, "MaxLag", MAX_LAG)
    Type = error.ParseKargs(kargs, "Type", "Pearson", type_dict)

    error.CheckType([variable, max_lag], [int, int])
    _check_nb_variable(obj, variable)

    #todo check that  Type is Pearson or Kendall
    return obj.partial_autocorrelation_computation(variable, Type, max_lag)

