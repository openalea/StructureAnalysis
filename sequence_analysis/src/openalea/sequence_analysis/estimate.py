""" Estimation functions """
__revision__ = "$Id:  $"

import sys,os

import openalea.stat_tool._stat_tool as _stat_tool
import _sequence_analysis
import sequences

__all__ = ['Estimate']
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

def _estimate_non_homogeneous_markov(obj, *args, **kargs):
 
    STAT_LINEAR = 0
    STAT_LOGISTIC = 1
    STAT_MONOMOLECULAR =2
    STAT_NONPARAMETRIC = 3

    ident_map = {
                 "VOID" :-1,
                 "MONOMOLECULAR":STAT_MONOMOLECULAR,
                 "LOGISTIC":STAT_LOGISTIC,
                 }
    CountingFlag = kargs.get("CountingFlag", True)
    
    ident = [ident_map[x] for x in args]
    
    return obj.nonhomogeneous_markov_estimation(ident, CountingFlag);


def _estimate_hidden_variable_order_markov(obj, *args, **kargs):
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
      
      
def _estimate_semi_markov(obj, *args, **kargs):
 
    PARTIAL_LIKELIHOOD = 0 
    COMPLETE_LIKELIHOOD  =1
    KAPLAN_MEIER=2
    estimator_semi_markov_map = {
        "PartialLikelihood": PARTIAL_LIKELIHOOD,
        "CompleteLikelihood": COMPLETE_LIKELIHOOD,
        "KaplanMeier": KAPLAN_MEIER           
                             }
    if isinstance(args[0], str):  
        
        if args[0] == "Ordinary":
            Type = 'o'
        elif args[0] == "Equilibrium":
            Type = 'e'
        else:
            raise AttributeError("type must be Ordinary or Equilibrium")
    else:
        raise AttributeError("type must be Ordinary or Equilibrium")
        
    NbIteration = kargs.get("NbIteration", -1)
    Counting = kargs.get("Counting", True)
    Estimator = kargs.get("Estimator", COMPLETE_LIKELIHOOD)
    
    Estimator = estimator_semi_markov_map[Estimator]
    MeanComputation = __parse_kargs__(kargs, "OccupancyMean",
                                      default='Computed',
                                      map=mean_computation_map)
    return obj.semi_markov_estimation(Type , Estimator , Counting,
                                        NbIteration , MeanComputation)


       
def _estimate_hidden_semi_markov(obj, *args, **kargs):
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
    OccupancyMean = kargs.get("InitialOccupancyMean", -1)
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
            print 'hidden semi markov 1'
            hsmarkov = obj.hidden_semi_markov_estimation_model( Type , NbState ,
                         LeftRight , Estimator , Counting , StateSequence ,
                         OccupancyMean ,NbIteration , MeanComputation)
            return hsmarkov

        elif Algorithm == "MCEM": #FORWARD_BACKWARD_SAMPLING : 
            #if Estimator == KAPLAN_MEIER:
            #    Estimator = COMPLETE_LIKELIHOOD
            print 'hidden semi markov 2'
            hsmarkov = obj.hidden_semi_markov_stochastic_estimation_model(
                Type, NbState, LeftRight, MinNbSequence, MaxNbSequence,
                Parameter, Estimator, Counting, StateSequence,
                OccupancyMean, NbIteration)
            return hsmarkov
        

    
    elif isinstance(args[0], _sequence_analysis._Hidden_semi_markov):
        
        hsmarkov = args[0] 
        if Algorithm == 'EM':
            print 'hidden semi markov 3'
            output = obj.hidden_semi_markov_estimation(hsmarkov,
                                Estimator, Counting, StateSequence,
                                NbIteration, MeanComputation)
            return output
        elif Algorithm == 'MCEM':
            print 'hidden semi markov 4'
            return obj.hidden_semi_markov_stochastic_estimation(hsmarkov,
                            MinNbSequence, MaxNbSequence,
                            Parameter, Estimator, Counting,
                            StateSequence, NbIteration)
    else:
        raise TypeError("should not reach this part of the code.check the usage")
    
    raise TypeError("should not reach this part of the code. check the usage")
    

def _estimate_variable_order_markov(obj, *args, **kargs):
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
       

    
    
    if len(args)>0 and isinstance(args[0], str):
        print "Type"
        
        
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

def _estimate_top(obj, *args, **kargs):
    """
    
    """
    MinPosition = kargs.get("MinPosition",1)
    MaxPosition = kargs.get("MaxPosition", obj.max_position)
    Neighbourhood = kargs.get("Neighbourhood", 1)
    if kargs.get("Neighborhood", None):
        Neighbourhood = kargs.get("Neighborhood", None)
    EqualProbability = kargs.get("EqualProbability", False)
    
    return obj.estimation(MinPosition, MaxPosition, Neighbourhood, 
                          EqualProbability)
     
    return None



def _estimate_dispatch(obj, itype, *args, **kargs):
    """
    
    """
    """    fct_map = {
        "VARIABLE_ORDER_MARKOV": "estimate_variable_order_markov",
        "HIDDEN_VARIABLE_ORDER_MARKOV": "estimate_hidden_variable_order_markov",
        "HIDDEN_SEMI-MARKOV": "estimate_hidden_semi_markov"
        }
    """
    #fct = getattr(obj, "estimate_%s" % fct_map[type] )
    fct_map_distribution = [
        "NONPARAMETRIC" ,
        "NP" ,
        "B" ,
        "BINOMIAL" ,
        "P" ,
        "POISSON" ,
        "NB" ,
        "NEGATIVE_BINOMIAL" ,
        "U" ,
        "UNIFORM" ,
        "MIXTURE" ,
        "CONVOLUTION" ,
        "COMPOUND",
        ]
    fct_map = ["VARIABLE_ORDER_MARKOV",
               "HIDDEN_VARIABLE_ORDER_MARKOV",
               "HIDDEN_SEMI-MARKOV",
               "SEMI-MARKOV",
               "NON-HOMOGENEOUS_MARKOV",
               "MARKOV"
               ]
    if (itype not in fct_map_distribution) and (itype not in fct_map):
        raise KeyError("Valid type are %s or %s" 
                       % (str(fct_map),
                          str(fct_map_distribution)))

        
        
    if itype == "VARIABLE_ORDER_MARKOV":    
        return _estimate_variable_order_markov(obj, *args, **kargs)
    elif itype == "MARKOV":    
        return _estimate_variable_order_markov(obj, *args, **kargs)
    elif itype == "HIDDEN_VARIABLE_ORDER_MARKOV":
        return _estimate_hidden_variable_order_markov(obj, *args, **kargs)
    elif itype == "HIDDEN_SEMI-MARKOV":
        return _estimate_hidden_semi_markov(obj, *args, **kargs)
    elif itype == "SEMI-MARKOV":
        return _estimate_semi_markov(obj, *args, **kargs)
    elif itype == "NON-HOMOGENEOUS_MARKOV":
        return _estimate_non_homogeneous_markov(obj, *args, **kargs)
    else:
        from openalea.stat_tool.estimate import Estimate as HistoEstimate
        return HistoEstimate(obj, itype, *args, **kargs)
    
    
    
def Estimate(obj, *args, **kargs):
    """Estimate
    
    * Estimation of distributions.
    * Estimation of 'top' parameters.
    * Estimation of a renewal process from count data.
    * Estimation of (hidden) Markovian models.

    :Usage:
    
    >>> Estimate(histo, "NON-PARAMETRIC")
    >>> Estimate(histo, "NB", MinInfBound->1, InfBoundStatus->"Fixed")
    >>> Estimate(histo, "MIXTURE", "B", dist,..., MinInfBound->1, InfBoundStatus->"Fixed", 
        DistInfBoundStatus->"Fixed") 
    >>> Estimate(histo, "MIXTURE", "B", "NB",..., MinInfBound->1, InfBoundStatus->"Fixed", 
        DistInfBoundStatus->"Fixed", NbComponent->"Estimated", Penalty->"AIC") 
    >>> Estimate(histo, "CONVOLUTION", dist,MinInfBound->1, Parametric->False) 
    >>> Estimate(histo, "CONVOLUTION", dist,InitialDistribution->initial_dist, Parametric->False) 
    >>> Estimate(histo, "COMPOUND", dist, unknown, Parametric->False, MinInfBound->0) 
    >>> Estimate(histo, "COMPOUND", dist, unknown, InitialDistribution->initial_dist, Parametric->False) 
    
    >>> Estimate(top, MinPosition->1, MaxPosition->5, Neighbourhood->2,    EqualProba->True)
    
    >>> Estimate(timev, type, NbIteration->10,Parametric->True)
    >>> Estimate(timev, type, InitialInterEvent->initial_dist,    NbIteration->10, Parametric->True)      
    
    >>> Estimate(seq, "MARKOV", Order->2, Counting->False)
    >>> Estimate(seq, "MARKOV", MaxOrder->3, Penalty->"AIC", Counting->False) 
    >>> Estimate(seq, "MARKOV", states, Penalty->"AIC", Order->2, Counting->False) 
    >>> Estimate(seq, "NON-HOMOGENEOUS_MARKOV", MONOMOLECULAR, VOID, Counting->False) 
    >>> Estimate(seq, "SEMI-MARKOV", Counting->False)
    >>> Estimate(seq, "HIDDEN_MARKOV", nb_state, structure, SelfTransition->0.9, NbIteration->10, 
        StateSequences->"Viterbi", Counting->False) 
    >>> Estimate(seq, "HIDDEN_MARKOV", hmc, Algorithm->"Viterbi",
        NbIteration->10, Order->2, Counting->False) 
    >>> Estimate(seq, "HIDDEN_MARKOV", "NbState", min_nb_state,
        max_nb_state, Penalty->"AIC", Order->2, Counting->False) 
    >>> Estimate(seq, "HIDDEN_MARKOV", "NbState", hmc, state,
        max_nb_state, Penalty->"AIC", SelfTransition->0.9, Counting->False) 
    >>> Estimate(seq, "HIDDEN_SEMI-MARKOV", nb_state, structure,
        OccupancyMean->20, NbIteration->10, Estimator->"PartialLikelihood",
        StateSequences->"Viterbi", Counting->False) 
    >>> Estimate(seq, "HIDDEN_SEMI-MARKOV", hsmc, Algorithm->"Viterbi", NbIteration->10, Counting->False) 
    
    :Arguments:
    
    * histo (histogram, mixture_data, convolution_data, compound_data),
    * dist (distribution, mixture, convolution, compound),
    * unknown (string): type of unknown distribution: "Sum" or "Elementary".
    * top (tops),
    * timev (time_events, renewal_data),
    * type (string): type or renewal process: "Ordinary" or "Equilibrium".
    * seq (discrete_sequences, markov_data, semi-markov_data),
    * states, ... (array(int)):  new states corresponding to a partition of the
      original state space, 
    * hmc (hidden_markov),
    * structure (string): structural properties of  the underlying Markov chain: 
      "Irreductible" or "LeftRight" (i.e. a succession of transient states and a 
      final absorbing state),
    * nb_state (int): number of states with 2 <= nb_state <=  15,
    * min_nb_state (int): minimum number of states,
    * max_nb_state (int): maximum number of states with 2 <= min_nb_state < max_nb_state <= 15
      or (number of states of the initial hidden Markov chain hmc) < max_nb_state<= 15.
    * state (int): state to be duplicated,
    * hsmc (hidden_semi-markov).

    :Optional Arguments:
    
    **distribution case**
    
    * MinInfBound (int): lower bound to the range of possible values (0 - default
      value - or 1). This optional argument cannot be used in conjunction with the 
      optional argument InitialDistribution.
    * InfBoundStatus (string): shifting or not of the distribution: "Free" (default
      value) or "Fixed". This optional argument cannot be used if the second mandatory 
      argument giving the model type is "NON-PARAMETRIC" ("NP").
    * DistInfBoundStatus (string): shifting or not of the subsequent components of
      the mixture: "Free" (default value) or "Fixed". This optional argument can 
      only be used if the second mandatory argument giving the distribution type is "MIXTURE".
    * NbComponent (string): estimation of the number of components of the mixture:
      "Fixed" (default value) or "Estimated". This optional argument can only be 
      used if the second mandatory argument giving the distribution type is "MIXTURE".
      the number of estimated components is comprised between 1 and a maximum number
      which is given by the number of specified parametric distributions in the 
      mandatory arguments (all of these distributions are assumed to be unknown).
    * Penalty (string): type of penalty function for model selection: "AIC" 
      (Akaike Information Criterion), "AICc" (corrected Akaike Information Criterion 
      - default value) or "BIC" (Bayesian Information Criterion). This optional 
      argument can only be used if the second mandatory argument giving the distribution 
      type is "MIXTURE" and if the optional argument NbComponent is set at "Estimated".
    * Parametric (bool): reestimation of a discrete nonparametric or parametric 
      distribution (default value: True). This optional argument can only be used if
      the second mandatory argument giving the distribution type is "CONVOLUTION" 
      or "COMPOUND".
    * InitialDistribution (distribution, mixture, convolution, compound): initial
      distribution for the EM deconvolution-type algorithm. This optional argument
      can only be used if the second mandatory argument giving the distribution type
      is "CONVOLUTION" or "COMPOUND". This optional argument cannot be used in 
      conjunction with the optional argument MinInfBound.
      
    .. note:: the optional arguments MinInfBound and InitialDistribution are mutually exclusive.
    
    
    **top case**
    
    * MinPosition (int): lower position taken into account for the estimation of 'top' parameters.
    * MaxPosition (int): upper position taken into account for the estimation of 'top' parameters.
    * Neighbourhood (int): neighbourhood taken into account for the estimation of 'top' parameters.
    * EqualProba (bool): growth probabilities of the parent shoot and of the offspring shoots equal or not (default value: False).
    
    **renewal case**
    
    * InitialInterEvent (distribution, mixture, convolution, compound): initial inter-event distribution for the EM algorithm.
    * NbIteration (int): number of iterations of the EM algorithm.
    * Parametric (bool): reestimation of a discrete nonparametric or parametric distribution (default value: False).
    
    **markovian case**
    
    * Counting (bool): computation of counting distributions (default value: True).
    * Order (int): Markov chain order (default value: 1). This optional argument can only be used if the second mandatory argument giving the model type is "MARKOV", "NON-HOMOGENEOUS_MARKOV" or "HIDDEN_MARKOV".
    * MaxOrder (int): maximum order of the Markov chain (default value: 4). This optional argument can only be used if the second mandatory argument giving the model type is "MARKOV".
    * Penalty (string): type of penalty function for model selection: "AIC" (Akaike Information Criterion), "AICc" (corrected Akaike Information Criterion) or "BIC" (Bayesian Information Criterion). This optional argument can only be used if the second mandatory argument giving the model type is "MARKOV" (default value: "BIC") and if the optional argument MaxOrder is set or else, if a new set of states is given (defining a partition of the original state space) or else, if the second mandatory argument giving the model type is "HIDDEN_MARKOV" and the third "NbState" (default value: "AICc").
    * Algorithm (string): type of algorithm: "ForwardBackward" (the default) or "Viterbi". This optional argument can only be used if the second mandatory argument giving the model type is "HIDDEN_MARKOV" or "HIDDEN_SEMI-MARKOV".
    * NbIteration (int): number of iterations of the estimation algorithm.
    * SelfTransition (real): self-transition probability. This optional argument can only be used if the second mandatory argument giving the model type is "HIDDEN_MARKOV" and if the initial model used in the iterative estimation procedure (EM algorithm) is only specified by its number of states, its structural properties and eventually its order.
    * OccupancyMean (int/real): average state occupancy. This optional argument can only be used if the second mandatory argument giving the model type is "HIDDEN_SEMI-MARKOV" and if the initial model used in the iterative estimation procedure (EM algorithm) is only specified by its number of states and its structural properties.
    * Estimator (string): type of estimator: "CompleteLikelihood" (the default) or "PartialLikelihood". In this latter case, the contribution of the time spent in the last visited state is not taken into account inthe estimation of the state occupancy distributions. This optional argument can only be used if the second mandatory argument giving the model type is "HIDDEN_SEMI-MARKOV" and the optional argument Algorithm is set at "ForwardBackward". 
    * StateSequences (string): Computation of the optimal state sequences: no computation (the default), "ForwardBackward" or "Viterbi". This optional argument can only be used if the second mandatory argument giving the model type is "HIDDEN_MARKOV" or "HIDDEN_SEMI-MARKOV" and if the optional argument Algorithm is not set at "Viterbi".
    
    
    :Returned Object:
    
    **distribution case**
    
    In case of success of the estimation procedure, the type of the returned object 
    (chosen among distribution, mixture, convolution or compound) is given by the 
    second mandatory argument. Otherwise no object is returned. The returned object
    of type distribution, mixture, convolution or compound contains both the estimated 
    distribution and the data used in the estimation procedure. In the case of
    mixtures, convolutions, or compound (or stopped-sum) distributions, the returned 
    object contains pseudo-data computed as a byproduct of the EM algorithm which 
    can be extracted by the function ExtractData.
    
    **top case**
    
    In case of success of the estimation procedure, an object of type top_parameters 
    is returned, otherwise no object is returned. The returned object of type top_parameters
    contains both the estimated model and the data used for the estimation.    

    **renewal case**
    
    In case of success of the estimation procedure, an object of type renewal is
    returned, otherwise no object is returned. The returned object of type renewal
    contains both the estimated renewal process and the count data used in the 
    estimation procedure.
    
    **markovian case**
    
    In case of success of the estimation procedure, the type of the returned object 
    (chosen among markov, semi-markov, hidden_markov, hidden_semi-markov) is given
    by the second mandatory argument. Otherwise no object is returned. If the 
    second mandatory argument is "NON-HOMOGENEOUS_MARKOV", in case of success
    of the estimation procedure, the returned object is of type markov. If the
    second mandatory argument is "NON-HOMOGENEOUS_MARKOV", the subsequent 
    arguments chosen among "VOID" (homogeneous state), "MONOMOLECULAR" or 
    "LOGISTIC", specify the evolution of the self-transition probabilities 
    as a function of the index parameter. The returned object of type markov,
    semi-markov, hidden_markov or hidden_semi-markov contains both the estimated 
    distribution and the data used in the estimation procedure. In the case of
    the estimation of a hidden Markov chain or a hidden semi-Markov chain, 
    the returned object contains pseudo-data (optimal state sequences 
    corresponding to the observed sequences used in the estimation procedure)
    computed as a byproduct of the EM algorithm which can be extracted by the
    function ExtractData.


    :Background:
    
    The aim of the model of 'tops' is to related the growth of offspring shoots 
    to the growth of their parent shoot in the case of immediate branching. A 
    model of 'tops' is defined by three parameters, namely the growth probability
    of the parent shoot, the growth probability of the offspring shoots (both in
    the sense of Bernoulli processes) and the growth rhythm ratio offspring 
    shoots / parent shoot.

    :Description (markovian case):
    
    In the case of hidden Markovian models (second mandatory argument "HIDDEN_MARKOV"
    or "HIDDEN_SEMI-MARKOV"), either the forward-backward algorithm or the Viterbi
    algorithm can be used for estimation. The Viterbi algorithm should only be
    used for the estimation of hidden Markovian models based on an underlying 
    "left-right" Markov chain (i.e. constituted of a succession of transient states 
    and a final absorbing state). Hence, in this case, the model structure is
    implicitly "LeftRight" and should not be given as argument (only the number of 
    states should be given as argument). Since the optimal state sequences are computed
    by the Viterbi algorithm, the optional argument StateSequences cannot be used if 
    the optional argument Algorithm is set at "Viterbi".

    .. seealso::
        
        :func:`~openalea.stat_tool.data_transform.ExtractData`, 
        :func:`~openalea.stat_tool.data_transform.ExtractDistribution`.
        :func:`~openalea.sequence_analysis.data_transform.AddAbsorbingRun`,
        :func:`ModelSelectionTest`.

    """
    if isinstance(args[0], str):
        return _estimate_dispatch(obj, args[0], *args[1:], **kargs)
    else:
        return _estimate_top(obj, *args, **kargs) 
    
