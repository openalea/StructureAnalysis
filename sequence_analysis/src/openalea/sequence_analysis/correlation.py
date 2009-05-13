"""Tops"""
__revision__ = "$Id: vectors.py 6217 2009-04-08 12:40:15Z cokelaer $"

import os
import openalea.stat_tool.interface as interface
import _sequence_analysis

from _sequence_analysis import _Correlation

__all__ = ['_Correlation','ComputeCorrelation', 'ComputeAutoCorrelation']


# Extend dynamically class
interface.extend_class( _Correlation, interface.StatInterface)

type_dict = {"Pearson": 0,
                 "Spearman": 1,
                 "Kendall":2,
                 "Spearman2":3}

def ComputeCorrelation(obj, *args, **kargs):
    """Computation of sample autocorrelation or cross-correlation functions.
   
    :Examples:
        >>>	ComputeCorrelation(seq1, MaxLag->10, Type->"Spearman", Normalization->"Exact")
	    >>> ComputeCorrelation(seqn, variable, MaxLag->10, Type->"Spearman", Normalization->"Exact")
        >>>	ComputeCorrelation(seqn, variable1, variable2, MaxLag->10, Type->"Spearman", Normalization->"Exact")	

    :Arguments:
    	seq1 (sequences, discrete_sequences, markov_data, semi-markov_data): univariate sequences,
	    seqn (sequences, discrete_sequences, markov_data, semi-markov_data): multivariate sequences,
    	variable (int): variable index (computation of a sample autocorrelation function).
	    variable1, variable2 (int): variable indices (computation of a sample cross-correlation function).
	
    :Optional Arguments:
    	Type (string): type of correlation coefficient: "Pearson" (linear correlation coefficient - default value), "Spearman" or "Kendall" (rank correlation coefficients).
	    MaxLag (int): maximum lag. A default value is computed from the sequence length distribution,
	    Normalization (STRING): normalization of the correlation coefficients: "Approximated" (the default - usual convention for time series analysis) or "Exact", (highly recommended for sample of short sequences). This optional argument can only be used if the optional argument Type is set at "Pearson" or "Spearman".
	
    :Returned Object:
    	If variable, or variable1 and variable2 are valid indices of variables (and are different if two indices are given) 
        and if 0 < MaxLag < (maximum length of sequences), an object of type correlation is returned, otherwise no object is returned.	
  
    :Background:
    	In the univariate case or if only variable is given, a sample autocorrelation function is computed. If variable1 and variable2 are given, a sample cross-correlation function is computed.
	
    .. seealso::
    	:class:`ComputePartialAutoCorrelation`
	    :class:`ComputeWhiteNoiseAutoCorrelation`
"""
    if obj.nb_variable==1:
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
    #check validit of the arguments
    max_lag = kargs.get("MaxLag", -1) #todo set default values
    user_type = kargs.get("Type", "Pearson")
    user_normalization = kargs.get("Normalization", "Exact")
    #todo check it is valid, i.e. in Pearson, SpearMan,Kendall,Spearman2
    
   
    
    norm_type = {"Approximated": 0,
                 "Exact": 1, }
    
    itype = type_dict[user_type]
    normalization = norm_type[user_normalization]
    
    return obj.correlation_computation(variable1, variable2, 
                                   itype, max_lag, normalization)
    
    

def ComputeAutoCorrelation(obj, *args, **kargs):
    """
    ComputeAutoCorrelation
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
    max_lag = kargs.get("MaxLag", MAXLAG) 
    
    
    if len(args) == 1:
        if isinstance(obj, _sequence_analysis._Variable_order_markov):
            return obj.state_autocorrelation_computation(
                value, max_lag)
    elif len(args)==2:
        return obj.output_autocorrelation_computation(
                variable, value, max_lag)

    

def ComputewhiteNoiseCorrelation(obj, itype):
     
    if isinstance(itype, int):
        obj.white_noise_correlation_order(itype)
    elif isinstance(itype, list):
        obj.white_noise_correlation_filter(itype)
    else:
        try:
            obj.white_noise_correlation_dist(itype)
        except TypeError:
            raise TypeError("second argument must be either an integer, a list or a Distribution type")
        
        
        
def ComputePartialCorrelation(obj, *args, **kargs):
    """
    ComputeAutoCorrelation
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
    max_lag = kargs.get("MaxLag", MAXLAG) 
   
    Type = kargs.get("Type", "Pearson")
    Type = type_dict[Type]
    
    if len(args) == 1:
        return obj.partial_autocorrelation_computation(
                variable, Type, max_lag)
            