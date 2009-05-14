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

sub_markovian_algorithms = {                    
    'EM':_stat_tool.RestorationAlgorithm.FORWARD_BACKWARD,           
    'MCEM':_stat_tool.RestorationAlgorithm.FORWARD_BACKWARD_SAMPLING,  
}               

algorithm = {
    'CTM_BIC':_sequence_analysis.Algorithm.CTM_BIC , # // algorithme Context Tree Maximizing/BIC
    'CTM_KT': _sequence_analysis.Algorithm.CTM_KT,   # // algorithme Context Tree Maximizing/Krichevsky-Trofimov
    'LocalBIC': _sequence_analysis.Algorithm.LOCAL_BIC , # // algorithme d'elagage recursif/BIC 
    'Context': _sequence_analysis.Algorithm.CONTEXT # // algorithme Context
}


estimator = {
    'MaximumLikelihood' : _sequence_analysis.Estimator.MAXIMUM_LIKELIHOOD ,
    'Laplace' :_sequence_analysis.Estimator.LAPLACE ,
    'AdaptativeLaplace' :_sequence_analysis.Estimator.ADAPTATIVE_LAPLACE ,
    'UniformSubset' :_sequence_analysis.Estimator.UNIFORM_SUBSET ,
    'UniformCardinality':_sequence_analysis.Estimator.UNIFORM_CARDINALITY
 
}

likelihood_penalty_type = {
    'AIC': _stat_tool.LikelihoodPenaltyType.AIC,
    'AICc': _stat_tool.LikelihoodPenaltyType.AICc,
    'BIC': _stat_tool.LikelihoodPenaltyType.BIC,
    'BICc': _stat_tool.LikelihoodPenaltyType.BICc,
    'ICL' : _stat_tool.LikelihoodPenaltyType.ICL,
    'ICLc': _stat_tool.LikelihoodPenaltyType.ICLc,
    }
#todo add this enumerate in boost_python 

COMPUTED = 0
ESTIMATED = 1
ONE_STEP_LATE = 2

mean_computation_map = {
    "Computed" : COMPUTED,
    "Estimated" : ESTIMATED,
    "OneStepLate" : ONE_STEP_LATE        
}

#todo add this enumerate in boost_python 
PARTIAL_LIKELIHOOD =0
COMPLETE_LIKELIHOOD =1
KAPLAN_MEIER =2

estimator_hidden_semi_markov = {
    "CompleteLikelihood" : COMPLETE_LIKELIHOOD,
    "PartialLikelihood" :  PARTIAL_LIKELIHOOD,
    "KaplanMeier" : KAPLAN_MEIER
}


def __parse_kargs__(kargs, key, default=None, map=None):
    """
    convert the key (string) into appropriate enumerate value
    map is a required dictionary 
    """
    user_choice = kargs.get(key, default)
    try:
        return map[user_choice]
    except KeyError:    
        raise KeyError("Wrong choice for %. Possible choices are %s " % (key, map.keys()))


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
      
       
def estimate_hidden_semi_markov(obj, *args, **kargs):
    """
    >>> hsmc21 = Estimate(seq21, "HIDDEN_SEMI-MARKOV", hsmc0)
    
    """

    """  nb_required , 
         
                
         mean_computation = COMPUTED
         
    """
    
    MIN_NB_STATE_SEQUENCE = 1
    MAX_NB_STATE_SEQUENCE = 10  
    NB_STATE_SEQUENCE_PARAMETER = 1
   
    Algorithm = kargs.get("Algorithm", "EM") #FORWARD_BACKWARD
    
    _AlgorithmCheck = __parse_kargs__(kargs, "Algorithm", default='EM',
                                map=sub_markovian_algorithms)
    Estimator = __parse_kargs__(kargs, "Estimator",
                                default='CompleteLikelihood',
                                map=estimator_hidden_semi_markov)
    MeanComputation = __parse_kargs__(kargs, "OccupancyMean",
                                      default='Computed',
                                      map=mean_computation_map)
        
    StateSequence = kargs.get("StateSequence", True)
    Counting = kargs.get("Counting", True)
    InitialOccupancyMean = kargs.get("InitialOccupancyMean", None)
    OccupancyMean = kargs.get("InitialOccupancyMean", None)
    NbIteration = kargs.get("NbIteration", -1)
    MinNbSequence = kargs.get("MinNbSequence", MIN_NB_STATE_SEQUENCE)
    MaxNbSequence = kargs.get("MaxNbSequence", MAX_NB_STATE_SEQUENCE)
    Parameter = kargs.get("Parameter",NB_STATE_SEQUENCE_PARAMETER )
     

#todo check this ?? not really needed for now      
#if ((algorithm != FORWARD_BACKWARD_SAMPLING) && (min_nb_state_sequence_option)) {
      #genAMLError(ERRORMSG(FORBIDDEN_OPTION_ss) , "Estimate" , "MinNbStateSequence");
#if ((algorithm != FORWARD_BACKWARD_SAMPLING) && (max_nb_state_sequence_option)) {
    #genAMLError(ERRORMSG(FORBIDDEN_OPTION_ss) , "Estimate" , "MaxNbStateSequence");
#if ((algorithm != FORWARD_BACKWARD_SAMPLING) && (parameter_option)) {
     #genAMLError(ERRORMSG(FORBIDDEN_OPTION_ss) , "Estimate" , "Parameter");

   
    if isinstance(args[0], str) or isinstance(args[0], int):  
        #note : do we really need the second isinstance ? 
        #this is the case in the AML code but later there is an error raised 
        #if args[0] is not a string... 
        
        if args[0] == "Ordinary":
            Type = 'o'
        elif args[0] == "Equilibrium":
            Type = 'e'
        else:
            raise AttributeError("type must be Ordinary or Equilibrium")
        
        if not isinstance(args[1], int):
            raise AttributeError("nb_state nnot provided. expected an integer after the type.")
        
        NbState = args[1]
        
        if Type == 'o':
            if args[2] not in ["LeftRight", "Irreducible"]:
                raise AttributeError("third arguments must be a strin : LeftRight or Irreducible.")
            if args[2] == "LeftRight":
                LeftRight = True
        elif Type == 'e':
            LeftRight = False
                 
      #todo traise error if (((type != 'e') || (estimator == PARTIAL_LIKELIHOOD) || (algorithm != FORWARD_BACKWARD)) &&(occupancy_mean_option)) {
    
        if Algorithm == "EM":
            hsmarkov = obj.hidden_semi_markov_estimation_model( Type , NbState ,
                         LeftRight , Estimator , Counting , StateSequence ,
                         OccupancyMean ,NbIteration , MeanComputation)

        elif Algorithm == "MCEM": #FORWARD_BACKWARD_SAMPLING : 
            #if Estimator == KAPLAN_MEIER:
            #    Estimator = COMPLETE_LIKELIHOOD

            hsmarkov = obj.hidden_semi_markov_stochastic_estimation_model(
                Type, NbState, LeftRight, MinNbSequence, MaxNbSequence,
                Parameter, Estimator, Counting, StateSequence,
                OccupancyMean, NbIteration)
        

    
    if isinstance(args[0], _sequence_analysis._Hidden_semi_markov):
        hsmarkov = args[0] 
        if Algorithm == 'EM':
            return obj.hidden_semi_markov_estimation(hsmarkov,
                                Estimator, Counting, StateSequence,
                                NbIteration, MeanComputation)
        elif Algorithm == 'MCEM':
            return obj.hidden_semi_markov_stochastic_estimation(hsmarkov,
                            MinNbSequence, MaxNbSequence,
                            Parameter, Estimator, Counting,
                            StateSequence, NbIteration)



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
       
    if isinstance(args[0], str):
        order_estimation = True
        if Order is not None:
            order_estimation = False
            MaxOrder = Order
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
        print "estimation4"
        symbol = args[0]
        markov = obj.lumpability_estimation(symbol, Penalty,
                                         Order, CountingFlag)
   
    else:
        raise KeyError("jfjf")
    
    return markov



def Estimate(obj, itype, *args, **kargs):
    """
    
    """
    """    fct_map = {
        "VARIABLE_ORDER_MARKOV": "estimate_variable_order_markov",
        "HIDDEN_VARIABLE_ORDER_MARKOV": "estimate_hidden_variable_order_markov",
        "HIDDEN_SEMI-MARKOV": "estimate_hidden_semi_markov"
        }
    """
    #fct = getattr(obj, "estimate_%s" % fct_map[type] )
    
    
    #if itype in fct_map.keys():
    if itype == "VARIABLE_ORDER_MARKOV":    
        return estimate_variable_order_markov(obj, *args, **kargs)
    elif itype == "HIDDEN_VARIABLE_ORDER_MARKOV":
        return estimate_hidden_variable_order_markov(obj, *args, **kargs)
    elif itype == "HIDDEN_SEMI-MARKOV":
        return estimate_hidden_semi_markov(obj, *args, **kargs)
    else:
        from openalea.stat_tool.estimate import Estimate as HistoEstimate
        return HistoEstimate(obj, itype, *args, **kargs)
    
    
    
    