#!/usr/bin/env python
#-*- coding: utf-8 -*-
"""Distributions

.. topic:: distribution.py summary

    Provides standard distributions

    :Code status: mature
    :Documentation status: to be completed
    :Author: Thomas Cokelaer <Thomas.Cokelaer@sophia.inria.fr>

    :Revision: $Id: distribution.py 17945 2015-04-07 07:44:42Z pradal $
"""
__version__ = "$Id: distribution.py 17945 2015-04-07 07:44:42Z pradal $"

import interface
import error

from openalea.stat_tool._stat_tool import *

from openalea.stat_tool._stat_tool import _DiscreteParametricModel
from openalea.stat_tool._stat_tool import _DiscreteDistributionData
from openalea.stat_tool._stat_tool import _Distribution
from openalea.stat_tool._stat_tool import I_DEFAULT
from openalea.stat_tool._stat_tool import D_DEFAULT
from openalea.stat_tool._stat_tool import D_INF
from openalea.stat_tool._stat_tool import MAX_DIFF_BOUND
from openalea.stat_tool._stat_tool import MAX_MEAN
from openalea.stat_tool._stat_tool import VariableType
from openalea.stat_tool._stat_tool import VariableTypeBis
from openalea.stat_tool._stat_tool import RestorationAlgorithm


from enums import distribution_identifier_type

__all__ = ["_Distribution",
           "_DiscreteParametricModel",
           "Distribution",
           "Binomial",
           "Poisson",
           "Uniform",
           "NegativeBinomial",
           "Multinomial",
           "ToHistogram",
           "ToDistribution"]

def Distribution(utype, *args):
    """
    Construction of a parametric discrete distribution (either binomial,
    Poisson, negative binomial or uniform) from the name and the parameters
    of the distribution or from an ASCII file.

    A supplementary shift parameter (argument inf_bound) is required with
    respect to the usual definitions of these discrete distributions.
    Constraints over parameters are given in the file syntax corresponding
    to the type distribution(cf. File Syntax).

    :Parameters:
      * `inf_bound` (int) : lower bound to the range of possible values
        (shift parameter),
      * `sup_bound` (int) : upper bound to the range of possible values \
      (only relevant for binomial or uniform distributions),
      * `param` (int, real) : parameter of either the Poisson distribution or \
      the negative binomial distribution.
      * `proba` (int, float) : probability of success \
      (only relevant for binomial or negative binomial distributions),
      * `file_name` (string).

      .. note:: the names of the parametric discrete distributions can be
        summarized by their first letters:

        * "B" ("BINOMIAL"),
        * "P" ("POISSON"),
        * "NB" ("NEGATIVE_BINOMIAL"),
        * "U" ("UNIFORM"),
        * "M" ("MULTINOMIAL"),


    :Returns:
        If the construction succeeds, an object of type distribution is
        returned, otherwise no object is returned.

    :Examples:

    .. doctest::
        :options: +SKIP

        >>> Distribution("BINOMIAL", inf_bound, sup_bound, proba)
        >>> Distribution("POISSON", inf_bound, param)
        >>> Distribution("NEGATIVE_BINOMIAL", inf_bound, param, proba)
        >>> Distribution("UNIFORM", inf_bound, sup_bound)
        >>> Distribution(file_name)

    .. seealso::
        :func:`~openalea.stat_tool.output.Save`,
        :func:`~openalea.stat_tool.estimate.Estimate`
        :func:`~openalea.stat_tool.simulate.Simulate`.
    """
    # Constructor from Filename or Histogram or parametricmodel
    if(len(args) == 0):
        error.CheckType([utype],
                        [[str, _DiscreteDistributionData, _DiscreteParametricModel]],
                        arg_id=[1])
        result =  _DiscreteParametricModel(utype)
    # from parameters
    if len(args)>0:
        error.CheckArgumentsLength(args, 1)
        if utype in ["B",  "BINOMIAL"]:
            result = Binomial(*args)
        elif utype in ["P", "POISSON"]:
            result = Poisson(*args)
        elif utype in ["M", "MULTINOMIAL"]:
            raise NotImplementedError("Multinomial not yet implemented")
        elif utype in ["NB", "NEGATIVE_BINOMIAL"]:
            result = NegativeBinomial(*args)
        elif utype in ["U", "UNIFORM"]:
            result = Uniform(*args)
        else:
            raise KeyError(" %s not found. Allowed keys are %s"
                           % (utype, distribution_identifier_type.keys()))

    return result


def Binomial(inf_bound, sup_bound=I_DEFAULT, \
             proba=D_DEFAULT):
    """
    Construction of a binomial distribution

    :param float inf_bound: lower bound to the range of possible values    (shift parameter)
    :param float sup_bound: upper bound to the range of possilbe values
    :param float proba: probability of `success`

    .. plot::
        :width: 50%
        :include-source:

        from openalea.stat_tool.distribution import Binomial
        b = Binomial(0,10,0.5)
        b.plot(legend_size=8)

    """
    # todo: seg fault when passing -1 as first arguments if there
    # is no assert here below
    # memory leak ?
    # todo:  returns error if ((inf_bound < min_inf_bound) ||
    # (inf_bound > MAX_INF_BOUND)) {

    error.CheckType([inf_bound, sup_bound, proba], [int, int, [int, float]])
    assert inf_bound >= 0
    assert inf_bound < sup_bound
    assert (sup_bound - inf_bound) <= MAX_DIFF_BOUND
    assert proba <= 1. and proba > 0

    param = D_DEFAULT

    return _DiscreteParametricModel(BINOMIAL,
        inf_bound, sup_bound, param, proba)


def Poisson(inf_bound, param=D_DEFAULT):
    """
    Construction of a poisson distribution

    :Parameters:
      * `inf_bound` (int) : lower bound to the range of possible values
        (shift parameter)
      * `param` (int, float) : parameter of the Poisson distribution
    
    .. plot::
        :width: 50%
        :include-source:

        from openalea.stat_tool.distribution import Poisson
        b = Poisson(1,1.5)
        b.plot(legend_size=8)
    """

    assert inf_bound >= 0
    assert param > inf_bound
    assert param > 0. and param < MAX_MEAN

    sup_bound = I_DEFAULT
    proba = D_DEFAULT

    return _DiscreteParametricModel(POISSON, \
        inf_bound, sup_bound, param, proba)


def NegativeBinomial(inf_bound, param=D_DEFAULT, \
                     proba=D_DEFAULT):
    """
    Construction of a negative binomial distribution

    :Parameters:
      * inf_bound (int) : lower bound to the range of possible values
        (shift parameter)
      * param (int, float) : parameter of the Poisson distribution
      * proba (int, float) : probability of 'success'

    .. plot::
        :width: 50%
        :include-source:

        from openalea.stat_tool.distribution import NegativeBinomial
        b = NegativeBinomial(1,2 ,.5)
        b.plot(legend_size=12)

    """

    # Check parameters
    assert inf_bound >= 0
    assert param > 0
    assert proba <= 1. and proba > 0
    assert (param * (1. - proba) / proba) <= MAX_MEAN

    sup_bound = I_DEFAULT

    return _DiscreteParametricModel(NEGATIVE_BINOMIAL, \
        inf_bound, sup_bound, param, proba)


def Uniform(inf_bound, sup_bound=I_DEFAULT):
    """
    Construction of a uniform distribution

    :Parameters:
      * inf_bound (int) : lower bound to the range of possible values
        (shift parameter)
      * sup_bound (int) : upper bound to the range of possilbe values

    .. plot::
        :width: 50%
        :include-source:

        from openalea.stat_tool.distribution import Uniform
        b = Uniform(1,10)
        b.plot(legend_size=8)

    """

    # Check parameters
    #todo:  returns error if ((inf_bound < min_inf_bound) ||
    # (inf_bound > MAX_INF_BOUND)) {
    assert inf_bound >= 0
    assert sup_bound >= 0
    assert inf_bound <= sup_bound
    assert (sup_bound - inf_bound) < MAX_DIFF_BOUND

    param = D_DEFAULT
    proba = D_DEFAULT
    cumul_threshold = CUMUL_THRESHOLD
    return _DiscreteParametricModel(UNIFORM, \
        inf_bound, sup_bound, param, proba, cumul_threshold)


def Multinomial():
    """to be done"""
    raise NotImplementedError("Multinomial not yet implemented")

# Extend _DiscreteParametricModel
interface.extend_class( _DiscreteParametricModel, interface.StatInterface)

# Cast Functions


def ToDistribution(histo):
    """ Cast an object of type `_DiscreteDistributionData` into an object
    of type `_Distribution`.

    :Parameters:
      * `histo` (DiscreteDistributionData)

    :Returns:
        If the object histo contains a 'model' part, an object
        of type `_Distribution` is returned, otherwise no object
        is returned.

    :Examples:

    .. doctest::
        :options: +SKIP

        >>> ToDistribution(histo)

    .. seealso::
        :func:`~openalea.stat_tool.distribution.ToHistogram`
    """
    return histo.extract_model()


def ToHistogram(dist):
    """Cast an object of type `_Distribution` into an object of
    type `_DiscreteDistributionData`.

    :Parameters:
      * dist (distribution).

    :Returns:
        If the object dist contains a 'data' part, an object of
        type `_DiscreteDistributionData` is returned, otherwise no object
        is returned.

    :Examples:

    .. doctest::
        :options: +SKIP

        >>> ToHistogram(dist)

    .. seealso::
        :func:`~openalea.stat_tool.distribution.ToDistribution`
    """
    return dist.extract_data()


