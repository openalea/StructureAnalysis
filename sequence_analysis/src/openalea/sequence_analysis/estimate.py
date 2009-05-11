""" Estimation functions """
__revision__ = "$Id: estimate.py 6167 2009-04-01 16:41:14Z cokelaer $"

import sys,os

import openalea.stat_tool._stat_tool as _stat_tool
import _sequence_analysis
import sequences

# to be checked or improve. Maybe the c++ code could be more explicit, i.e.,
# e switched to elementary and so on.

stochastic_process_type = {
    'Ordinary': 'o',
    'Equilibrium' : 'e'
    }


algorithm = {
  'CTM_BIC':_sequence_analysis.Algorithm.CTM_BIC , # // algorithme Context Tree Maximizing/BIC
  'CTM_KT': _sequence_analysis.Algorithm.CTM_KT,   # // algorithme Context Tree Maximizing/Krichevsky-Trofimov
  'LOCAL_BIC': _sequence_analysis.Algorithm.LOCAL_BIC , # // algorithme d'elagage recursif/BIC 
  'CONTEXT': _sequence_analysis.Algorithm.CONTEXT # // algorithme Context
};


estimator = {
             
  'MaximumLikelihood' : _sequence_analysis.Estimator.MAXIMUM_LIKELIHOOD ,
  'Laplace' :_sequence_analysis.Estimator.LAPLACE ,
  'AdaptativeLaplace' :_sequence_analysis.Estimator.ADAPTATIVE_LAPLACE ,
  'UniformSubset' :_sequence_analysis.Estimator.UNIFORM_SUBSET ,
  'UniformCardinality':_sequence_analysis.Estimator.UNIFORM_CARDINALITY
 
}

likelihood_penalty_type = {
    'AIC': _stat_tool.AIC,
    'AICc': _stat_tool.AICc,
    'BIC': _stat_tool.BIC,
    'BICc': _stat_tool.BICc,
    'ICL' : _stat_tool.ICL,
    'ICLc': _stat_tool.ICLc,
    }



def estimate_markovian_sequences(obj, type, *args, **kargs):
    """

    """

    #fct_map = {
    #    "VARIABLE_ORDER_MARKOV" : variable_order_markov_estimation,
     #   }

    ORDER = 8
    LOCAL_BIC_THRESHOLD = 10

    Order = kargs.get("Order", 1)
    MaxOrder = kargs.get("MaxOrder", ORDER) 
    MinOrder = kargs.get("MinOrder", 0)
    Algorithm = kargs.get("Algorithm", _sequence_analysis.Algorithm.LOCAL_BIC)
    Threshold = kargs.get("Threshold", LOCAL_BIC_THRESHOLD)
    Estimator = kargs.get("Estimator", _sequence_analysis.Estimator.LAPLACE)
    GlobalInitialTransition = kargs.get("GlobalInitialTransition", True)
    GlobalSample = kargs.get("GlobalSample", True)
    CountingFlag = kargs.get("CountingFlag", True)
    Penalty = kargs.get("Penalty", likelihood_penalty_type['BIC'])
    # check type
 
    
    if Order :
        order_estimation = False
        
    
    if isinstance(args[0], str):
        try:
            if (args[0]):
                Type = stochastic_process_type[args[0]]
        except KeyError:
            raise AttributeError("Bad type. Possible types are %s" 
                                 % (str(stochastic_process_type.keys())))
         
        if order_estimation is True:
             markov = obj.variable_order_markov_estimation1( 
                Type, MinOrder, MaxOrder, Algorithm, Threshold, Estimator , 
                  GlobalInitialTransition , GlobalSample , CountingFlag);
        else:
              markov = obj.variable_order_markov_estimation2(
                    Type, MaxOrder, GlobalInitialTransition, CountingFlag)
            
    elif isinstance(args[0], _sequence_analysis._Variable_order_markov):
        markov = obj.variable_order_markov_estimation3(args[0], 
                      GlobalInitialTransition, CountingFlag)
    

    elif isinstance(args[0], list):
        symbol = args[0]
        markov = obj.lumpability_estimation(symbol, Penalty,
                                         Order, CountingFlag)
   
    else:
        raise KeyError("jfjf")
    
    return markov

