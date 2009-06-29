""" Estimation functions """
__revision__ = "$Id$"

import sys,os
sys.path.append(os.path.abspath("."))

import _stat_tool
import interface
import distribution

# to be checked or improve. Maybe the c++ code could be more explicit, i.e.,
# e switched to elementary and so on.
compound_type = {
    'Elementary': 'e',
    'Sum': 's',
    }

likelihood_penalty_type = {
    'AIC': _stat_tool.AIC,
    'AICc': _stat_tool.AICc,
    'BIC': _stat_tool.BIC,
    'BICc': _stat_tool.BICc,
    'ICL' : _stat_tool.ICL,
    'ICLc': _stat_tool.ICLc,
    }

smoothing_penalty_type = {
    'FirstDifference': _stat_tool.FIRST_DIFFERENCE,
    'SecondDifference': _stat_tool.SECOND_DIFFERENCE,
    'Entropy': _stat_tool.ENTROPY,
    }

outside_type = {
    "Zero": _stat_tool.ZERO,
    "Continuation": _stat_tool.CONTINUATION,
    }


estimator_type = {
    "Likelihood": _stat_tool.LIKELIHOOD,
    "PenalizedLikelihood": _stat_tool.PENALIZED_LIKELIHOOD,
    "Parametric": _stat_tool.PARAMETRIC_REGULARIZATION,
    }


class EstimateFunctions(object):
    """
    Class containing histogram estimation functions
    This class must not be used alone, but through an histogram object
    """

    def estimate_nonparametric( histo):
        """
        Estimate a non parametric distribution

        :Parameters:
          * histo (histogram, mixture_data, convolution_data, compound_data)

        :Examples:

        .. doctest::
            :options: +SKIP

            >>> estimate_nonparametric(histo)
        """
        return  _stat_tool._ParametricModel(histo)

    def estimate_parametric(histo, ident, MinInfBound=0, 
                            InfBoundStatus="Free"):
        """ Estimate a parametric discrete distribution (binomial, Poisson or negative binomial distribution with an additional shit parameter)

        :Parameters:
          * histo (histogram, mixture_data, convolution_data, compound_data),
          * ident ("BINOMIAL", "POISSON", "NEGATIVE_BINOMIAL", "UNIFORM")
          * MinInfBound (int): lower bound to the range of possible values (0 - default value - or 1).
          * InfBoundStatus (string): shifting or not of the distribution:
                                      "Free" (default value) or "Fixed". T

        :Examples:

        .. doctest::
            :options: +SKIP

            >>> estimate_parametric(histo, ident, MinInfBound=0, InfBoundStatus="Free")

        """

        flag = bool(InfBoundStatus == "Free")

        map_ident = distribution.distribution_type

        try:
            ident_id = map_ident[ident]
        except KeyError:
            raise KeyError("Valid type are %s"%(str(map_ident.keys())))

        return histo.parametric_estimation(ident_id, MinInfBound, flag)




    def estimate_mixture(histo, distributions,
                         MinInfBound=0, InfBoundStatus="Free",
                         DistInfBoundStatus="Free",
                         NbComponent = "Fixed", Penalty='AIC'):
        """ Estimate a finite  mixture of discrete distributions


        :Parameters:

          * histo (histogram, mixture_data, convolution_data, compound_data),
          * distributions (list) : a list of distribution object
                                   or distribution label(string) : 'B', 'NB', 'U', 'P', ...

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

            >>> estimate_mixture(histo, ("MIXTURE", "B", dist,...,),
                             MinInfBound=1, InfBoundStatus="Fixed",
                             DistInfBoundStatus="Fixed")
            >>> estimate_mixture(histo, ("MIXTURE", "B", "NB",...,),
                               MinInfBound=1, InfBoundStatus="Fixed",
                               DistInfBoundStatus="Fixed",
                               NbComponent="Estimated", Penalty="AIC")

        """

        estimate = [] # list of bool
        pcomponent = [] # list of distribution (parametric)
        ident = [] # list of distribution identifier

        # Parse list of distribution
        for d in distributions:

            estimate.append(False)
            pcomponent.append(None)
            ident.append(None)

            try:
                type = distribution.get_distribution_type(d,
                                                          [_stat_tool.BINOMIAL,
                                                           _stat_tool.POISSON,
                                                           _stat_tool.NEGATIVE_BINOMIAL,]
                                                          )
                # todo: why _stat_tool.UNIFORM is not included ?

                estimate[-1] = True
                ident[-1] = type
                pcomponent[-1] = _stat_tool._Parametric(0, type)

            except AttributeError:
                if(not isinstance(d, str)):
                    pcomponent[-1] = _stat_tool._Parametric(d)
                else:
                    raise

        flag = bool(InfBoundStatus == "Free")
        component_flag = bool(DistInfBoundStatus == "Free")
        nb_component_estimation = bool(NbComponent == "Estimated")

        try:
            if(Penalty):
                Penalty = likelihood_penalty_type[Penalty]

        except KeyError:
            raise AttributeError("Bad penalty. Possible penalty are %s"%(
                    str(likelihood_penalty_type.keys())))

        # check parameters
        if(not nb_component_estimation and Penalty):
            raise TypeError("Penalty can only be used with NbComponent set to 'Estimated'")

        if(not nb_component_estimation): # "FIXED"
            imixt = _stat_tool._Mixture(pcomponent)
            return histo.mixture_estimation(imixt, estimate, MinInfBound,
                                            flag, component_flag)

        else:  # "ESTIMATED"
            return histo.mixture_estimation(ident, MinInfBound, flag,
                                            component_flag, Penalty) 

    def estimate_compound(histo, known_distribution, Type,
                          InitialDistribution=None,
                          MinInfBound=0,  Weight=-1. , NbIteration=-1,
                          Penalty="SecondDifference", Outside="Zero",
                          Estimator="Likelihood"):
        """estimate a compound"""


 
        try:
            if (Type):
                Type = compound_type[Type]
        except KeyError:
            raise AttributeError("Bad type. Possible types are %s"
                                 % (str(compound_type.keys())))

        try:
            if(Estimator):
                Estimator = estimator_type[Estimator]
        except KeyError:
            raise AttributeError("Bad estimator. Possible estimator are %s"
                                 % (str(estimator_type.keys())))

        try:
            if(Penalty):
                Penalty = smoothing_penalty_type[Penalty]
        except KeyError:
            raise AttributeError("Bad penalty. Possible penalty are %s"
                                 % (str(smoothing_penalty_type.keys())))

        try:
            if(Outside):
                Outside = outside_type[Outside]
        except KeyError:
            raise AttributeError("Bad side effect management type: should be %s"
                                 % (str(outside_type.keys())))


        #print InitialDistribution
        #print known_distribution

        if (InitialDistribution):
            raise("to be checked carefully")
            print InitialDistribution
            return histo.compound_estimation1(
                            known_distribution, unknown_distribution, Type,
                            Estimator, NbIteration, Weight, Penalty, Outside)

        else:
            return histo.compound_estimation2(
                            known_distribution, Type, MinInfBound,  Estimator,
                            NbIteration, Weight, Penalty, Outside)


    def estimate_convolution( histo, known_distribution,
                             InitialDistribution=None,
                             MinInfBound=0, Estimator="Likelihood",
                             NbIteration=-1, Weight=-1.,
                             Penalty="SecondDifference", Outside="Zero"):

        """ Estimate a convolution


        """

        try:
            if(Estimator):
                Estimator = estimator_type[Estimator]

        except KeyError:
            raise AttributeError("Bad estimator. Possible estimator are %s"
                                 % (str(estimator_type.keys())))

        try:
            if(Penalty):
                Penalty = smoothing_penalty_type[Penalty]

        except KeyError:
            raise AttributeError("Bad penalty. Possible penalty are %s"
                                 % (str(smoothing_penalty_type.keys())))

        try:
            if(Outside):
                Outside = outside_type[Outside]

        except KeyError:
            raise AttributeError("Bad side effect management type: should be %s"
                                 % (str(outside_type.keys())))


        if (InitialDistribution) :
            return histo.convolution_estimation(known_distribution, 
                                                InitialDistribution, Estimator,
                                                NbIteration, Weight, Penalty, 
                                                Outside)

        else :
            return histo.convolution_estimation(known_distribution, MinInfBound,
                                                Estimator, NbIteration, Weight,
                                                Penalty, Outside)


# Extend _Histogram class
_Histogram = interface.extend_class( _stat_tool._Histogram, EstimateFunctions)




def Estimate(histo, itype, *args, **kargs):
    """
    Estimation function for AML compatibility

    .. seealso::
        :func:`~openalea.stat_tool.estimate.EstimateFunctions.estimation_nonparametric`,
        :func:`~openalea.stat_tool.estimate.estimation_parametric`,
        :func:`~openalea.stat_tool.estimate.estimation_mixture`,
        :func:`~openalea.stat_tool.estimate.estimation_convolution`,
        :func:`~openalea.stat_tool.estimate.estimation_compound`,
    """

    fct_map = {
        "NONPARAMETRIC" : _Histogram.estimate_nonparametric,
        "NP" : _Histogram.estimate_nonparametric,
        "B" :  _Histogram.estimate_parametric,
        "BINOMIAL" :  _Histogram.estimate_parametric,
        "P" :  _Histogram.estimate_parametric,
        "POISSON" :  _Histogram.estimate_parametric,
        "NB" :  _Histogram.estimate_parametric,
        "NEGATIVE_BINOMIAL" :  _Histogram.estimate_parametric,
        "U" : _Histogram.estimate_parametric,
        "UNIFORM" : _Histogram.estimate_parametric,
        "MIXTURE" : _Histogram.estimate_mixture,
        "CONVOLUTION" : _Histogram.estimate_convolution,
        "COMPOUND": _Histogram.estimate_compound,
        }

    itype = itype.upper()

    try:
        fct = fct_map[itype]
        if (fct == _Histogram.estimate_parametric):
            return fct(histo, itype, *args, **kargs)
        elif (fct == _Histogram.estimate_mixture):
            return fct(histo, args, **kargs)
#       elif (fct == _Histogram.estimate_compound):
#            return fct(histo, *args, **kargs)
        else:
            return fct(histo, *args, **kargs)

    except KeyError:
        raise KeyError("Valid type are %s" % (str(fct_map.keys())))




