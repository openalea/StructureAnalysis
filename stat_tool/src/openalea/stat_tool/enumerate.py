"""common enumerate

:Author: Thomas Cokelaer <Thomas.Cokelaer@inria.fr>

"""
__version__ = "$Id$"

import _stat_tool


# map to enumerate in boost python
criterion_type = {
    "NearestNeighbor":  _stat_tool.CriterionType.NEAREST_NEIGHBOR,
    "FarthestNeighbor": _stat_tool.CriterionType.FARTHEST_NEIGHBOR,
    "Averaging":        _stat_tool.CriterionType.AVERAGING,
    }

distance_type = {
    "ABSOLUTE_VALUE":   _stat_tool.DistanceType.ABSOLUTE_VALUE,
    "QUADRATIC":        _stat_tool.DistanceType.QUADRATIC
    }

fitting_type = {
    "STANDARD_NORMAL": _stat_tool.FittingType.STANDARD_NORMAL,
    "CHI2": _stat_tool.FittingType.CHI2,
    "FISHER": _stat_tool.FittingType.FISHER,
    "STUDENT": _stat_tool.FittingType.STUDENT
    }

round_type = {
    "FLOOR":_stat_tool.RoundType.FLOOR,
    "ROUND":_stat_tool.RoundType.ROUND,
    "CEIL":_stat_tool.RoundType.CEIL
    }

distribution_identifier_type = {
    "B":                _stat_tool.DistributionIdentifierType.BINOMIAL,
    "BINOMIAL":         _stat_tool.DistributionIdentifierType.BINOMIAL,
    "P":                _stat_tool.DistributionIdentifierType.POISSON,
    "POISSON":          _stat_tool.DistributionIdentifierType.POISSON,
    "NB":               _stat_tool.DistributionIdentifierType.NEGATIVE_BINOMIAL,
    "NEGATIVE_BINOMIAL":_stat_tool.DistributionIdentifierType.NEGATIVE_BINOMIAL,
    "U": _stat_tool.DistributionIdentifierType.UNIFORM,
    "UNIFORM": _stat_tool.DistributionIdentifierType.UNIFORM,
    "M": _stat_tool.DistributionIdentifierType.MULTINOMIAL,
    "MULTINOMIAL": _stat_tool.DistributionIdentifierType.MULTINOMIAL,
    }

variable_type = {
    "O" : _stat_tool.VariableType.ORDINAL,
    "ORDINAL" : _stat_tool.VariableType.ORDINAL,
    "N" : _stat_tool.VariableType.NUMERIC,
    "NUMERIC" : _stat_tool.VariableType.NUMERIC,
    "S" : _stat_tool.VariableType.SYMBOLIC,
    "SYMBOLIC" : _stat_tool.VariableType.SYMBOLIC,
    "C" : _stat_tool.VariableType.CIRCULAR,
    "CIRCULAR" : _stat_tool.VariableType.CIRCULAR,
    }

pearson_type = {
    "PEARSON": _stat_tool.PearsonType.PEARSON,
    "SPEARMAN": _stat_tool.PearsonType.SPEARMAN,
    "KENDALL": _stat_tool.PearsonType.KENDALL,
    "SPEARMAN2": _stat_tool.PearsonType.SPEARMAN2
    }

smoothing_penalty_type = {
    'FirstDifference': _stat_tool.FIRST_DIFFERENCE,
    'SecondDifference': _stat_tool.SECOND_DIFFERENCE,
    'Entropy': _stat_tool.ENTROPY,
    }

estimator_type = {
    "Zero":                 _stat_tool.EstimatorType.ZERO,
    "Likelihood":           _stat_tool.EstimatorType.LIKELIHOOD,
    "PenalizedLikelihood":  _stat_tool.EstimatorType.PENALIZED_LIKELIHOOD,
    "Parametric": _stat_tool.EstimatorType.PARAMETRIC_REGULARIZATION,
    "ParametricRegularization": 
        _stat_tool.EstimatorType.PARAMETRIC_REGULARIZATION
    }

outside_type = {
    "Zero": _stat_tool.ZERO,
    "Continuation": _stat_tool.CONTINUATION,
    }

#  enum_<stat_tool::wrap_util::UniqueInt<3, 11> >("TOBEDONEType")
#  .value("PARTIAL_LIKELIHOOD", PARTIAL_LIKELIHOOD)
#  .value("COMPLETE_LIKELIHOOD", COMPLETE_LIKELIHOOD)
#  .value("KAPLAN_MEIER", KAPLAN_MEIER)
#  .export_values()
#  ;


#  enum_<stat_tool::wrap_util::UniqueInt<3, 12> >("ClusterType")
#  .value("COMPUTED", COMPUTED)
#  .value("ESTIMATED", ESTIMATED)
#  .value("ONE_STEP_LATE", ONE_STEP_LATE)
#  .export_values()
#  ;

likelihood_penalty_type = {
    'AIC': _stat_tool.AIC,
    'AICc': _stat_tool.AICc,
    'BIC': _stat_tool.BIC,
    'BICc': _stat_tool.BICc,
    'ICL' : _stat_tool.ICL,
    'ICLc': _stat_tool.ICLc,
    }

algorithm_type = {
     "Agglomerative":   _stat_tool.AlgoType.AGGLOMERATIVE,
     "Divisive":        _stat_tool.AlgoType.DIVISIVE,
     "Ordering":        _stat_tool.AlgoType.ORDERING,
     }




# used by Clustering
format_type = { \
    "ASCII" :'a',
    "SpreadSheet": 's',
    "" : 'n'
    }

cluster_type = \
    {
     "Step":        "cluster_step",
     "Limit":       "cluster_limit",
     "Information": "cluster_information"
     }

compound_type = {
    'Elementary': 'e',
    'Sum': 's',
    }

variance_type = {
    "N": _stat_tool.NUMERIC,
    "NUMERIC" : _stat_tool.NUMERIC,
    "O": _stat_tool.ORDINAL,
    "ORDINAL":  _stat_tool.ORDINAL,
    }

keep_type = {
             "keep": "Keep",
             "Keep": "Keep",
             "reject": "Reject",
             "Reject": "Reject",
             "Rejected": "Reject"
             }
