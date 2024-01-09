#!/usr/bin/env python
#-*- coding: utf-8 -*-
"""Estimate functions

.. topic:: estimate.py summary

    A module dedicated to estimation functionalities

    :Code status: mature
    :Documentation status: to be completed
    :Author: Thomas Cokelaer <Thomas.Cokelaer@sophia.inria.fr>

    :Revision: $Id: estimate.py 15183 2013-11-06 10:35:50Z jbdurand $
    
.. warning:: sequence analysis package also contains an estimate module and
    function
"""

__version__ = "$Id: estimate.py 15183 2013-11-06 10:35:50Z jbdurand $"

import error
import interface

from enums import likelihood_penalty_type
from enums import smoothing_penalty_type
from enums import outside_type
from enums import compound_type
from enums import estimator_type
from enums import distribution_identifier_type as dist_type

from openalea.stat_tool._stat_tool import _DiscreteParametricModel
from openalea.stat_tool._stat_tool import _DiscreteParametric
from openalea.stat_tool._stat_tool import _Compound
from openalea.stat_tool._stat_tool import _Convolution
from openalea.stat_tool._stat_tool import _Distribution
from openalea.stat_tool._stat_tool import _DiscreteMixture
from openalea.stat_tool._stat_tool import _FrequencyDistribution
from openalea.stat_tool._stat_tool import LikelihoodPenaltyType

__all__ = ["Estimate", "EstimateFunctions"]

class EstimateFunctions(object):
    """
    Class containing histogram estimation functions
    This class must not be used alone, but through an histogram object
    """

    def estimate_nonparametric(histo):
        """
        Estimate a non parametric distribution

        :Parameters:
          * histo (histogram, mixture_data, convolution_data, compound_data)

        :Usage:

        .. doctest::
            :options: +SKIP

            >>> Estimate(histo, "NON-PARAMETRIC")
            >>> estimate_nonparametric(histo)

        """
        return  _DiscreteParametricModel(histo)


    def estimate_parametric(histo, ident,
                            MinInfBound=0,
                            InfBoundStatus="Free"):
        """ Estimate a parametric discrete distribution (binomial,
        Poisson or negative binomial distribution with an additional shift
        parameter)

        :Parameters:
          * histo (histogram, mixture_data, convolution_data, compound_data),
          * ident ("BINOMIAL", "POISSON", "NEGATIVE_BINOMIAL", "UNIFORM")
          * MinInfBound (int): lower bound to the range of possible values (0 - default value - or 1).
          * InfBoundStatus (string): shifting or not of the distribution:
                                      "Free" (default value) or "Fixed". T

        :Usage:

        .. doctest::
            :options: +SKIP

            >>> estimate_parametric(histo, ident, MinInfBound=0, InfBoundStatus="Free")
            >>> Estimate(histo, "NB", MinInfBound=1, InfBoundStatus="Fixed")

        """

        error.CheckType([ident, MinInfBound, InfBoundStatus],
                        [str, int, str])

        flag = bool(InfBoundStatus == "Free")

        try:
            ident_id = dist_type[ident]
        except KeyError:
            raise KeyError("Valid type are %s" % (str(dist_type.keys())))


        return histo.parametric_estimation(ident_id, MinInfBound, flag)




    def estimate_DiscreteMixture(histo, *args, **kargs):


        """ Estimate a finite  mixture of discrete distributions


        :Parameters:

          * histo (histogram, mixture_data, convolution_data, compound_data),
          * distributions (list) : a list of distribution object
                                   or distribution label(string) : 'B', 'NB', 'U', 'P', ...
          * unknown (string): type of unknown distribution: "Sum" or "Elementary".

        :Keywords:

          * MinInfBound (int): lower bound to the range of possible values (0 -default- or 1). \
                               This optional argument cannot be used in conjunction \
                               with the optional argument InitialDistribution.
          * InfBoundStatus (string): shifting or not of the distribution: "Free" (default value) or "Fixed".
          * DistInfBoundStatus (string): shifting or not of the subsequent components of \
                                         the mixture: "Free" (default value) or "Fixed".
          * NbComponent (string): estimation of the number of components of the mixture: \
                                  "Fixed" (default value) or "Estimated". Le number of estimated \
                                  components is comprised between\
                                  1 and a maximum number which is given by the number of specified \
                                  parametric distributions in the mandatory arguments \
                                  (all of these distributions are assumed to be unknown).
          * Penalty (string): type of Penalty function for model selection: \
                              "AIC" (Akaike Information Criterion), \
                              "AICc" (corrected Akaike Information Criterion) \
                              "BIC" (Bayesian Information Criterion - default value). \
                              "BICc" (corrected Bayesian Information Criterion). \

                              This optional argument can only be used if the optional argument
                              NbComponent is set at "Estimated".

        :Examples:

        .. doctest::
            :options: +SKIP

            >>> estimate_DiscreteMixture(histo, "MIXTURE", "B", dist,...,,
                             MinInfBound=1, InfBoundStatus="Fixed",
                             DistInfBoundStatus="Fixed")
            >>> estimate_DiscreteMixture(histo, "MIXTURE", "B", "NB",...,,
                               MinInfBound=1, InfBoundStatus="Fixed",
                               DistInfBoundStatus="Fixed",
                               NbComponent="Estimated", Penalty="AIC")
            >>> Estimate(histo, "MIXTURE", "B", dist, MinInfBound=1, InfBoundStatus="Fixed",
                    DistInfBoundStatus="Fixed")
            >>> Estimate(histo, "MIXTURE", "B", "NB",
                    MinInfBound=1, InfBoundStatus="Fixed",
                    DistInfBoundStatus="Fixed",
                    NbComponent="Estimated", Penalty="AIC")


        """
        #alias

        #error.CheckArgumentsLength(args, 1, 1)

        # get user arguments
        # list of distributions can be either a list or several arguments
        # e.g.: estimate_DiscreteMixture(["B","B"]) or estimate_DiscreteMixture("B", "B")
        if len(args)==1 and type(args[0])==list:
            distributions = args[0]
        else:
            distributions = list(args)

        InfBoundStatus = kargs.get("InfBoundStatus","Free")
        DistInfBoundStatus = kargs.get("DistInfBoundStatus", "Free")
        NbComponent = kargs.get("NbComponent", "Fixed")


        MinInfBound = kargs.get("MinInfBound", 0)
        Penalty = error.ParseKargs(kargs, "Penalty", "AIC",
                                likelihood_penalty_type)

        #should be before the conversion to booleans
        error.CheckType([MinInfBound, InfBoundStatus, DistInfBoundStatus,
                         NbComponent, Penalty],
                        [int, str, str, str, LikelihoodPenaltyType])


        # transform into boolean when needed
        InfBoundStatus = bool(InfBoundStatus == "Free")
        DistInfBoundStatus = bool(DistInfBoundStatus == "Free")
        NbComponent = bool(NbComponent == "Estimated")


        estimate = [] # list of bool
        pcomponent = [] # list of distribution
        ident = [] # list of distribution identifier

        # Parse list of distribution that could be defined by a distribution,
        # compound, mixture, convolution or simplya string such as "B",
        # "Mixture", ...

        for dist in distributions:

            if isinstance(dist, str):
                dist_authorised = ["BINOMIAL", "B", "POISSON",
                                   "P", "NB", "NEGATIVE_BINOMIAL"]
                if dist not in dist_authorised:
                    raise ValueError("""If distribution is a string, then it
                        must be in %s. You provided %s"""
                        % (dist_authorised, dist))
                #todo: check that poisson is allowed

                pcomponent.append(_DiscreteParametric(0, dist_type[dist]))
                ident.append(dist_type[dist])
                estimate.append(True)
            elif isinstance(dist, _DiscreteParametricModel):
                pcomponent.append(_DiscreteParametric(dist))
                ident.append(None)
                estimate.append(False)
            elif type(dist) in [_DiscreteMixture, _Convolution, _Compound]:
                pcomponent.append(_Distribution(dist))
                ident.append(None)
                estimate.append(False)
            else:
                raise ValueError("""In the case of a MIXTURE estimation,
                argument related to distributions must be either string, or
                Distribution, Mixture, Convolution, Compound. %s provided"""
                % dist)

        # check parameters
        if not NbComponent and Penalty:
            raise TypeError("""
            Penalty can only be used with NbComponent set to 'Estimated'""")

        if not NbComponent: # "FIXED"
            imixt = _DiscreteMixture(pcomponent)
            ret = histo.mixture_estimation1(imixt, estimate, MinInfBound,
                                            InfBoundStatus, DistInfBoundStatus)

            return ret
        else:  # "ESTIMATED"
            ret = histo.mixture_estimation2(ident, MinInfBound, InfBoundStatus,
                                            DistInfBoundStatus, Penalty)
            return ret

    def estimate_compound(histo, *args, **kargs):
        """estimate a compound


        :Usage:

        .. doctest::
            :options: +SKIP

            >>> Estimate(histo, "COMPOUND", dist, unknown,
                    Parametric=False, MinInfBound=0)
                    Estimate(histo, "COMPOUND", dist, unknown,
                    InitialDistribution=initial_dist, Parametric=False)
        """

        if len(args)<2:
            raise ValueError("expect at least three arguments")

        known_distribution = args[0]
        ##if isinstance(known_distribution, _DiscreteParametricModel):
        #    known_distribution = _DiscreteParametric(known_distribution)
        #elif type(known_distribution) in [_DiscreteMixture, _Convolution, _Compound]:
        #    known_distribution = _Distribution(known_distribution)
        #else:
        #    raise TypeError("""
        #    argument "known_distribution" must be of type _DiscreteMixture,
        #     _COnvolution, _Compound or _DiscreteParametricModel""")

        Type = args[1]
        error.CheckType([Type], [str])

        Weight = kargs.get("Weight", -1)
        NbIteration = kargs.get("NbIteration", -1)
        InitialDistribution = kargs.get("InitialDistribution", None)
        MinInfBound = kargs.get("MinInfBound", 0)

        Estimator = error.ParseKargs(kargs, "Estimator", "Likelihood",
                                estimator_type)
        Penalty = error.ParseKargs(kargs, "Penalty", "SecondDifference",
                                smoothing_penalty_type)
        Outside = error.ParseKargs(kargs, "Outside", "Zero", outside_type)

        if MinInfBound and InitialDistribution:
            raise ValueError("""MinInfBound and InitialDistribution cannot be
                             used together.""")
        #if Estimator != _stat_tool.PENALIZED_LIKELIHOOD:
        #    if Penalty or Weight or Outside:
        #        raise ValueError("""Estimator cannot be used with O
        #            utside or Weight or Penalty option""")

	    #The second argument is either a string (e.g.,"Sum") or an unknown
        #distribution.
        try:
            if Type:
                Type = compound_type[Type]
        except KeyError:
            raise AttributeError("Bad type. Possible types are %s"
                                 % (str(compound_type.keys())))

        #The second argument is either a string (e.g.,"Sum") or an unknown
        #distribution.
        unknown_distribution = None

        if InitialDistribution:
            unknown_distribution = InitialDistribution
            if isinstance(unknown_distribution, _Distribution):
                unknown_distribution = _DiscreteParametric(unknown_distribution)
            elif type(unknown_distribution) in \
                [_DiscreteMixture, _Convolution, _Compound]:
                unknown_distribution = _Distribution(unknown_distribution)
            else:
                raise TypeError("""
                    argument "known_distribution" must be of type
                     _DiscreteMixture, _COnvolution, _Compound or _DiscreteParametricModel""")
            if Type == 's':

                return histo.compound_estimation1(
                    unknown_distribution, known_distribution, Type,
                    Estimator, NbIteration, Weight, Penalty, Outside)
            elif Type == 'e':

                return histo.compound_estimation1(
                           known_distribution, unknown_distribution, Type,
                           Estimator, NbIteration, Weight, Penalty, Outside)
            else:
                raise KeyError("should not enter here.")	
        else:
            return histo.compound_estimation2(
                            known_distribution, Type, MinInfBound,  Estimator,
                            NbIteration, Weight, Penalty, Outside)


    def estimate_convolution(histo, *args, **kargs):
        """ Estimate a convolution

        :Usage:

        .. doctest::
            :options: +SKIP

            >>> Estimate(histo, "CONVOLUTION", dist,
                    MinInfBound=1, Parametric=False)
                    Estimate(histo, "CONVOLUTION", dist,
                    InitialDistribution=initial_dist, Parametric=False)

        """

        if len(args)==0:
            raise ValueError("expect at least two argument")
        known_distribution = args[0]
        Weight = kargs.get("Weight", -1)
        NbIteration = kargs.get("NbIteration", -1)
        InitialDistribution = kargs.get("InitialDistribution", None)
        MinInfBound = kargs.get("MinInfBound", 0)

        Estimator = error.ParseKargs(kargs, "Estimator", "Likelihood",
                                estimator_type)
        Penalty = error.ParseKargs(kargs, "Penalty", "SecondDifference",
                                smoothing_penalty_type)
        Outside = error.ParseKargs(kargs, "Outside", "Zero", outside_type)

        if isinstance(known_distribution, _DiscreteParametricModel):
            known_distribution = _DiscreteParametric(known_distribution)
        elif type(known_distribution) in [_DiscreteMixture, _Convolution, _Compound]:
            known_distribution = _Distribution(known_distribution)
        else:
            raise TypeError("""
            argument "known_distribution" must be of type _DiscreteMixture, _COnvolution,
            _Compound or _DiscreteParametricModel""")

        if InitialDistribution:
            return histo.convolution_estimation1(known_distribution,
                                                InitialDistribution, Estimator,
                                                NbIteration, Weight, Penalty,
                                                Outside)

        else :
            return histo.convolution_estimation2(known_distribution,
                                                 MinInfBound,
                                                 Estimator, NbIteration, Weight,
                                                 Penalty, Outside)


# Extend _Histogram class
_FrequencyDistribution = interface.extend_class( _FrequencyDistribution, EstimateFunctions)




def Estimate(histo, itype, *args, **kargs):
    """Estimate function

    This function is a dispatcher to several estimate functions depending on
    the first argument and the type.

    :param obj: the input object (may be histogram, sequence, compound, ...)
    :param itype: string.

    .. seealso::
        :func:`~openalea.stat_tool.estimate.EstimateFunctions.estimate_nonparametric`,
        :func:`~openalea.stat_tool.estimate.EstimateFunctions.estimate_parametric`,
        :func:`~openalea.stat_tool.estimate.EstimateFunctions.estimate_DiscreteMixture`,
        :func:`~openalea.stat_tool.estimate.EstimateFunctions.estimate_convolution`,
        :func:`~openalea.stat_tool.estimate.EstimateFunctions.estimate_compound`,
    """

    fct_map = {
        "NONPARAMETRIC" : _FrequencyDistribution.estimate_nonparametric,
        "NP" : _FrequencyDistribution.estimate_nonparametric,
        "B" :  _FrequencyDistribution.estimate_parametric,
        "BINOMIAL" :  _FrequencyDistribution.estimate_parametric,
        "P" :  _FrequencyDistribution.estimate_parametric,
        "POISSON" :  _FrequencyDistribution.estimate_parametric,
        "NB" :  _FrequencyDistribution.estimate_parametric,
        "NEGATIVE_BINOMIAL" :  _FrequencyDistribution.estimate_parametric,
        "U" : _FrequencyDistribution.estimate_parametric,
        "UNIFORM" : _FrequencyDistribution.estimate_parametric,
        "MIXTURE" : _FrequencyDistribution.estimate_DiscreteMixture,
        "CONVOLUTION" : _FrequencyDistribution.estimate_convolution,
        "COMPOUND": _FrequencyDistribution.estimate_compound,
        }

    Type = itype.upper()

    # sequence analysis case
    if Type not in fct_map.keys():
        try:
            from openalea.sequence_analysis.estimate import Estimate \
                as SeqEstimate
        except:
            raise ImportError("Could not import sequence_analysis")
        return SeqEstimate(histo, itype, *args, **kargs)
    else:
        fct = fct_map[Type]

        if fct == _FrequencyDistribution.estimate_parametric:
            return fct(histo, Type, *args, **kargs)
        elif fct == _FrequencyDistribution.estimate_DiscreteMixture:
            return fct(histo, *args, **kargs)
        else:
            return fct(histo, *args, **kargs)

