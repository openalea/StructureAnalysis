"""common enumerate

:Author: Thomas Cokelaer <Thomas.Cokelaer@inria.fr>

"""
__version__ = "$Id$"

from openalea.stat_tool._stat_tool import *

import openalea.stat_tool._stat_tool as _stat_tool
from openalea.sequence_analysis import _sequence_analysis as _sequence_analysis

import openalea.stat_tool.enumerate as enumerate_st



#todo: get rid if type_dict, rename the enumerate properly

type_dict = enumerate_st.pearson_type



norm_type = {"Approximated": _sequence_analysis.NormType.APPROXIMATED,
             "Exact": _sequence_analysis.NormType.EXACT }

# not complete !! see MarkovianSequnece
markovian_sequence_type = {
    "Length": _sequence_analysis.MarkovianSequenceType.LENGTH,
    "NbRun": _sequence_analysis.MarkovianSequenceType.NB_RUN,
    "NbOccurrence": _sequence_analysis.MarkovianSequenceType.NB_OCCURRENCE,
    "FirstOccurrence": \
        _sequence_analysis.MarkovianSequenceType.FIRST_OCCURRENCE,
    "Mean": _sequence_analysis.MarkovianSequenceType.SEQUENCE_MEAN,
    "Cumul": _sequence_analysis.MarkovianSequenceType.SEQUENCE_CUMUL
        }



#todo: get rid if mode_type, rename the enumerate properly
mode_type = enumerate_st.round_type

INTER_EVENT = _sequence_analysis.RenewalType.INTER_EVENT
WITHIN_OBSERVATION_PERIOD =  \
    _sequence_analysis.RenewalType.WITHIN_OBSERVATION_PERIOD
LENGTH_BIAS =  _sequence_analysis.RenewalType.LENGTH_BIAS
BACKWARD_RECURRENCE_TIME = \
    _sequence_analysis.RenewalType.BACKWARD_RECURRENCE_TIME
FORWARD_RECURRENCE_TIME = _sequence_analysis.RenewalType.FORWARD_RECURRENCE_TIME
NB_EVENT = _sequence_analysis.RenewalType.NB_EVENT
MIXTURE = _sequence_analysis.RenewalType.MIXTURE
     
renewal_nb_event_map = { 
    "InterEvent":   _sequence_analysis.RenewalType.INTER_EVENT,
    "LengthBias":   _sequence_analysis.RenewalType.LENGTH_BIAS,
    "Backward":     _sequence_analysis.RenewalType.BACKWARD_RECURRENCE_TIME,
    "Forward":      _sequence_analysis.RenewalType.FORWARD_RECURRENCE_TIME,
    "Mixture":      _sequence_analysis.RenewalType.MIXTURE,
    "Within":       _sequence_analysis.RenewalType.WITHIN_OBSERVATION_PERIOD,
    "NbEvent":      _sequence_analysis.RenewalType.NB_EVENT,
                        }

sub_func_map = {
    "Sequence": _sequence_analysis.OutputType.SEQUENCE,
    "SubtractionResidual": _sequence_analysis.OutputType.SUBTRACTION_RESIDUAL ,
    "DivisionResidual":_sequence_analysis.OutputType.DIVISION_RESIDUAL,
    "Residual": _sequence_analysis.OutputType.SUBTRACTION_RESIDUAL,

    }
func_map = {
    "Sequence": _sequence_analysis.OutputType.SEQUENCE,
    "Trend": _sequence_analysis.OutputType.TREND,
    "SubtractionResidual": _sequence_analysis.OutputType.SUBTRACTION_RESIDUAL ,
    "Residual": _sequence_analysis.OutputType.SUBTRACTION_RESIDUAL,
    "DivisionResidual":_sequence_analysis.OutputType.DIVISION_RESIDUAL,
    "StandardizedResidual": _sequence_analysis.OutputType.STANDARDIZED_RESIDUAL
    }
#todo: simplify
output_sequence_map = func_map
output_map = func_map 

estimator_map = {
    "MaximumLikelihood": _sequence_analysis.Estimator.MAXIMUM_LIKELIHOOD,
    "Laplace": _sequence_analysis.Estimator.LAPLACE,
    "AdaptativeLaplace": _sequence_analysis.Estimator.ADAPTATIVE_LAPLACE,
    "UniformSubset": _sequence_analysis.Estimator.UNIFORM_SUBSET,
    "UniformCardinality": _sequence_analysis.Estimator.UNIFORM_CARDINALITY
    }
estimator = estimator_map


SEQUENCE = 0
SUBTRACTION_RESIDUAL = 1
STANDARDIZED_RESIDUAL = 2
    
seq_map = {"Observation":0, 
               "FirstOccurrence":1,
               "Recurrence":2, 
               "Sojourn":3,
               "NbRun":6,
               "NbOccurrence":7,
               "Forward": -1}
    
_INITIAL_RUN = 4
_FINAL_RUN = 5
_LENGTH = 8
_SEQUENCE_CUMUL = 9
_SEQUENCE_MEAN = 10
   
histogram_type = {
                  "FinalRun":5,
                  "InitialRun":4,
                  }

model_type = {
              "Multinomial": _sequence_analysis.ChangeType.MULTINOMIAL_CHANGE,
              "Poisson": _sequence_analysis.ChangeType.POISSON_CHANGE,
              "Ordinal": _sequence_analysis.ChangeType.ORDINAL_GAUSSIAN_CHANGE,
              "Gaussian": _sequence_analysis.ChangeType.GAUSSIAN_CHANGE,        
              "Mean": _sequence_analysis.ChangeType.MEAN_CHANGE,
              "Variance": _sequence_analysis.ChangeType.VARIANCE_CHANGE,
              "MeanVariance": _sequence_analysis.ChangeType.MEAN_VARIANCE_CHANGE
              }

# todo:: ckean the 2 next enums ?
type_map = {
            "INT": _sequence_analysis.INT_VALUE, 
            "REAL" : _sequence_analysis.REAL_VALUE,
            "STATE": _sequence_analysis.STATE,
            "NB_INTERNODE":  _sequence_analysis.NB_INTERNODE,
            "AUXILIARY":  _sequence_analysis.AUXILIARY,
            }
    
index_parameter_type_map = {
        "IMPLICIT_TYPE": _sequence_analysis.IMPLICIT_TYPE,
        "TIME": _sequence_analysis.TIME,
        "TIME_INTERVAL": _sequence_analysis.TIME_INTERVAL,
        "POSITION": _sequence_analysis.POSITION,
        "POSITION_INTERVAL": _sequence_analysis.POSITION_INTERVAL
        }

 
stochastic_process_type = {
    'Ordinary': 'o',
    'Equilibrium' : 'e'
    }

markovian_algorithms = {
    'Forward':_stat_tool.RestorationAlgorithm.FORWARD,                    
    'EM':_stat_tool.RestorationAlgorithm.FORWARD_BACKWARD,           
    'MCEM':_stat_tool.RestorationAlgorithm.FORWARD_BACKWARD_SAMPLING,  
    'ForwardBackwardSampling': \
        _stat_tool.RestorationAlgorithm.FORWARD_DYNAMIC_PROGRAMMING,
    'GeneralizedViterbi': \
        _stat_tool.RestorationAlgorithm.GENERALIZED_VITERBI,
    'Gibbs':_stat_tool.RestorationAlgorithm.GIBBS_SAMPLING,             
    'NoComputation':_stat_tool.RestorationAlgorithm.NO_COMPUTATION,             
    'Viterbi':_stat_tool.RestorationAlgorithm.VITERBI,
    }

sub_markovian_algorithms = {                    
    'EM':_stat_tool.RestorationAlgorithm.FORWARD_BACKWARD,           
    'MCEM':_stat_tool.RestorationAlgorithm.FORWARD_BACKWARD_SAMPLING,  
    }               


sub_markovian_algorithms_2 = {
    'Forward':_stat_tool.RestorationAlgorithm.FORWARD,                    
    'Viterbi':_stat_tool.RestorationAlgorithm.VITERBI,
    }

  
algorithm = {
    'CTM_BIC': _sequence_analysis.Algorithm.CTM_BIC , 
    'BIC': _sequence_analysis.Algorithm.CTM_BIC ,
    'CTM_KT': _sequence_analysis.Algorithm.CTM_KT,   
    'KT': _sequence_analysis.Algorithm.CTM_KT,
    'LocalBIC': _sequence_analysis.Algorithm.LOCAL_BIC, 
    'Context': _sequence_analysis.Algorithm.CONTEXT
}



likelihood_penalty_type = {
    'AIC': _stat_tool.LikelihoodPenaltyType.AIC,
    'AICc': _stat_tool.LikelihoodPenaltyType.AICc,
    'BIC': _stat_tool.LikelihoodPenaltyType.BIC,
    'BICc': _stat_tool.LikelihoodPenaltyType.BICc,
    'ICL' : _stat_tool.LikelihoodPenaltyType.ICL,
    'ICLc': _stat_tool.LikelihoodPenaltyType.ICLc,
    }


mean_computation_map = {
    "Computed" : _stat_tool.ClusterType.COMPUTED,
    "Estimated" : _stat_tool.ClusterType.ESTIMATED,
    "OneStepLate" : _stat_tool.ClusterType.ONE_STEP_LATE        
    }

estimator_semi_markov_type = {
    "CompleteLikelihood" : _stat_tool.EstimatorHSMType.COMPLETE_LIKELIHOOD,
    "PartialLikelihood" :  _stat_tool.EstimatorHSMType.PARTIAL_LIKELIHOOD,
    "KaplanMeier" : _stat_tool.EstimatorHSMType.KAPLAN_MEIER
    }


ident_map = {
             "VOID" : _stat_tool.I_DEFAULT,
             "LINEAR": _stat_tool.RegressionType.STAT_LINEAR,
             "MONOMOLECULAR": _stat_tool.RegressionType.STAT_MONOMOLECULAR,
             "LOGISTIC": _stat_tool.RegressionType.STAT_LOGISTIC,
             "NONPARAMETRIC": _stat_tool.RegressionType.STAT_NONPARAMETRIC,

             }

output_type = { 
               "ChangePoint" : _sequence_analysis.ChangePointType.CHANGE_POINT,
               "Segment" : _sequence_analysis.ChangePointType.SEGMENT
               }

nb_segment_map = {
                  "Fixed": False,
                  "Estimated":True 
                  }


begin_aligned_map = {
                  "Aligned": False,
                  "Free":True 
                  }
end_aligned_map = begin_aligned_map

from openalea.sequence_analysis._sequence_analysis import \
    _Sequences, _MarkovianSequences, _VariableOrderMarkovData, \
    _SemiMarkovData, _NonHomogeneousMarkovData, _VariableOrderMarkov,\
    _HiddenVariableOrderMarkov, _HiddenSemiMarkov, _SemiMarkov
    
sequence_alignment_first_arg = [_Sequences,
                                _MarkovianSequences,
                                _VariableOrderMarkovData,
                                _SemiMarkovData,
                                _NonHomogeneousMarkovData]
    
markov_model_comparison_first_arg = \
    [_VariableOrderMarkov,
     _HiddenVariableOrderMarkov,
     _HiddenSemiMarkov,
     _SemiMarkov]
    
markov_model_for_sequences_first_arg = \
    [_MarkovianSequences,
     _VariableOrderMarkovData,
     _SemiMarkovData,
     _NonHomogeneousMarkovData]
   
markov_model_for_sequences_second_arg = \
    [_VariableOrderMarkov,
     _SemiMarkov,
     _HiddenVariableOrderMarkov,
     _HiddenSemiMarkov]

ms_vomd_smd_list = [_MarkovianSequences,
                    _VariableOrderMarkovData,
                    _SemiMarkovData]

# todo: use those simpler naming convention: 
ms_vomd_smd_nhmd = markov_model_for_sequences_first_arg

output_sequence = {
                   "DistanceMatrix": 'm',
                   "Sequences": 's'
                   }

indel_cost_map = {"Adaptative": _sequence_analysis.IndelCost.ADAPTATIVE, 
                     "Fixed": _sequence_analysis.IndelCost.FIXED} 