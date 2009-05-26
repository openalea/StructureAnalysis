"""Correlation methods

.. author:: Thomas Cokelaer, Thomas.Cokelaer@inria.fr

"""
__revision__ = "$Id: $"



import openalea.stat_tool.interface as interface
import _sequence_analysis
import tools

from _sequence_analysis import _Correlation

__all__ = ['_Correlation', 
           'ComputeCorrelation',
           'ComputeAutoCorrelation',
           'ComputeWhiteNoiseCorrelation',
           'ComputePartialAutoCorrelation']


# Extend dynamically class
interface.extend_class( _Correlation, interface.StatInterface)

type_dict = {"Pearson": 0,
             "Spearman": 1,
             "Kendall": 2,
             "Spearman2": 3}

norm_type = {"Approximated": 0,
                 "Exact": 1, }

def ComputeCorrelation(obj, *args, **kargs):
    """Computation of sample autocorrelation or cross-correlation functions.

    :Examples:

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
    if obj.nb_variable == 1:
        variable1 = 1
        variable2 = 1
    else:
        variable1 = args[0]
        if len(args)==1:
            variable2 = variable1
        elif len(args)==2:
            variable2 = args[1]
        else:
            raise TypeError("1 or 2  non-optional arguments required")
    
    
    max_lag = tools.__parse_kargs__(kargs, "MaxLag", -1)
    itype = tools.__parse_kargs__(kargs, "Type", "Pearson", type_dict) 
    normalization = tools.__parse_kargs__(kargs, "Normalization", "Exact", 
                                          norm_type)

    return obj.correlation_computation(variable1, variable2,
                                   itype, max_lag, normalization)



def ComputeAutoCorrelation(obj, *args, **kargs):
    """
    ComputeAutoCorrelation
    """

    if len(args) == 1:
        variable = 1
        value = args[0]
    elif len(args) == 2:
        variable = args[0]
        value = args[1]

    MAXLAG = 100
    max_lag = tools.__parse_kargs__(kargs, "MaxLag", MAXLAG)


    if len(args) == 1:
        if isinstance(obj, _sequence_analysis._Variable_order_markov):
            return obj.state_autocorrelation_computation(
                value, max_lag)
    elif len(args)==2:
        return obj.output_autocorrelation_computation(
                variable, value, max_lag)



def ComputeWhiteNoiseCorrelation(obj, itype):
    """ComputeWhiteNoiseCorrelation

    Computation of the sample autocorrelation or cross-correlation function induced on a white noise sequence by filtering.

    :Usage:

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

    if isinstance(itype, int):
        obj.white_noise_correlation_order(itype)
    elif isinstance(itype, list):
        obj.white_noise_correlation_filter(itype)
    else:
        try:
            obj.white_noise_correlation_dist(itype)
        except TypeError:
            raise TypeError("second argument must be either an integer, a list or a Distribution type")

    return obj


def ComputePartialAutoCorrelation(obj, *args, **kargs):
    """ComputePartialAutoCorrelation

    Computation of sample partial autocorrelation functions.

    :Usage:

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


    #if obj.nb_variable==1 and len(args)==1:
    if len(args) == 1:
        variable = 1
        value = args[0]
    #elif obj.nb_variable!=1 and len(args)==2:
    elif len(args)==2:
        variable = args[0]
        value = args[1]

    MAXLAG = 100
    max_lag = tools.__parse_kargs__(kargs, "MaxLag", MAXLAG)
    Type = tools.__parse_kargs__(kargs, "Type", "Pearson", type_dict) 

    if len(args) == 1:
        return obj.partial_autocorrelation_computation(variable, Type, max_lag)

