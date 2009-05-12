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

markovian_algorithms = {
'Forward':_stat_tool.RestorationAlgorithm.FORWARD,                    
'EM':_stat_tool.RestorationAlgorithm.FORWARD_BACKWARD,           
'MCEM':_stat_tool.RestorationAlgorithm.FORWARD_BACKWARD_SAMPLING,  
'ForwardBackwardSampling':_stat_tool.RestorationAlgorithm.FORWARD_DYNAMIC_PROGRAMMING,
'GeneralizedViterbi':_stat_tool.RestorationAlgorithm.GENERALIZED_VITERBI,
'Gibbs':_stat_tool.RestorationAlgorithm.GIBBS_SAMPLING,             
'NoComputation':_stat_tool.RestorationAlgorithm.NO_COMPUTATION,             
'Viterbi':_stat_tool.RestorationAlgorithm.VITERBI,
}
                       

algorithm = {
  'CTM_BIC':_sequence_analysis.Algorithm.CTM_BIC , # // algorithme Context Tree Maximizing/BIC
  'CTM_KT': _sequence_analysis.Algorithm.CTM_KT,   # // algorithme Context Tree Maximizing/Krichevsky-Trofimov
  'LocalBIC': _sequence_analysis.Algorithm.LOCAL_BIC , # // algorithme d'elagage recursif/BIC 
  'Context': _sequence_analysis.Algorithm.CONTEXT # // algorithme Context
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

def estimate_hidden_variable_order_markov(obj, *args, **kargs):
    """
    
    """
    #default values
    MIN_NB_STATE_SEQUENCE = 1
    MAX_NB_STATE_SEQUENCE = 10  
    NB_STATE_SEQUENCE_PARAMETER = 1.
    
    GlobalInitialTransition = kargs.get("GlobalInitialTransition", True)
    NbIteration = kargs.get("NbIteration", 80)
    CountingFlag = kargs.get("CountingFlag", True)
    StateSequence = kargs.get("StateSequence", True)
    Parameter = kargs.get("Parameter", NB_STATE_SEQUENCE_PARAMETER)
    MinNbSequence = kargs.get("MinNbSequence", MIN_NB_STATE_SEQUENCE)
    MaxNbSequence = kargs.get("MaxNbSequence", MAX_NB_STATE_SEQUENCE)
    Algorithm = kargs.get("Algorithm", 'EM')
    
    
    #check that a proper algorithm has been chosen
    if Algorithm not in markovian_algorithms.keys():
        raise KeyError("Algorithm must be in %s " % markovian_algorithms.keys())
    
    if not isinstance(args[0], _sequence_analysis._Hidden_variable_order_markov):
        raise TypeError('First argument must be a hidden_variable_order_markov')
    
    if Algorithm == 'EM':
        hmarkov = obj.hidden_variable_order_markov_estimation(
                args[0], GlobalInitialTransition,
                CountingFlag, StateSequence, NbIteration)

    elif Algorithm == 'MCEM':
        hmarkov = obj.hidden_variable_order_markov_stochastic_estimation(
                        args[0], GlobalInitialTransition, MinNbSequence,
                        MaxNbSequence, Parameter, CountingFlag,
                        StateSequence, NbIteration)

    return hmarkov
      
       



def estimate_variable_order_markov(obj, *args, **kargs):
    """

    """

    #fct_map = {
    #    "VARIABLE_ORDER_MARKOV" : variable_order_markov_estimation,
     #   }

    ORDER = 8
    LOCAL_BIC_THRESHOLD = 10

    Order = kargs.get("Order", None)
    MaxOrder = kargs.get("MaxOrder", ORDER) 
    MinOrder = kargs.get("MinOrder", 0)
    Algorithm = kargs.get("Algorithm", "LocalBIC")
    Threshold = kargs.get("Threshold", LOCAL_BIC_THRESHOLD)
    Estimator = kargs.get("Estimator", "Laplace")
    GlobalInitialTransition = kargs.get("GlobalInitialTransition", True)
    GlobalSample = kargs.get("GlobalSample", True)
    CountingFlag = kargs.get("CountingFlag", True)
    Penalty = kargs.get("Penalty","BIC")
   
    #convert to integers 
    try:
        Estimator = estimator[Estimator]
    except:
        pass
    try:
        Algorithm = algorithm[Algorithm]
    except:
        pass
    try:
        Penalty = likelihood_penalty_type[Penalty]
    except:
        pass
    
#    Estimator = estimator[Estimator]
    #Algorithm = algorithm[Algorithm]        
    #Penalty = likelihood_penalty_type[Penalty]
    
    
       
    if isinstance(args[0], str):
        order_estimation = True
        if Order is not None:
            order_estimation = False
        try:
            if (args[0]):
                Type = stochastic_process_type[args[0]]
        except KeyError:
            raise AttributeError("Bad type. Possible types are %s" 
                                 % (str(stochastic_process_type.keys())))
         
        if order_estimation is True:
            print "esitmation1" 

            
            markov = obj.variable_order_markov_estimation1( 
                Type, MinOrder, MaxOrder, Algorithm, Threshold, Estimator , 
                  GlobalInitialTransition , GlobalSample , CountingFlag);
        else:
            print "esitmation2"
            markov = obj.variable_order_markov_estimation2(
                    Type, MaxOrder, GlobalInitialTransition, CountingFlag)
            
    elif isinstance(args[0], _sequence_analysis._Variable_order_markov):
        print "esitmation3"
        markov = obj.variable_order_markov_estimation3(args[0], 
                      GlobalInitialTransition, CountingFlag)
    

    elif isinstance(args[0], list):
        symbol = args[0]
        markov = obj.lumpability_estimation(symbol, Penalty,
                                         Order, CountingFlag)
   
    else:
        raise KeyError("jfjf")
    
    return markov



def Estimate(obj, itype, *args, **kargs):
    """
    
    """
    fct_map = {
        "VARIABLE_ORDER_MARKOV": "estimate_variable_order_markov",
        "HIDDEN_VARIABLE_ORDER_MARKOV": "estimate_hidden_variable_order_markov"
        }
    
    #fct = getattr(obj, "estimate_%s" % fct_map[type] )
    
    
    if itype in fct_map.keys():
        if itype == "VARIABLE_ORDER_MARKOV":    
            return estimate_variable_order_markov(obj, *args, **kargs)
        elif itype == "HIDDEN_VARIABLE_ORDER_MARKOV":
            return estimate_hidden_variable_order_markov(obj, *args, **kargs)
    else:
        from openalea.stat_tool.estimate import Estimate as HistoEstimate
        return HistoEstimate(obj, itype, *args, **kargs)
    
    
    
    