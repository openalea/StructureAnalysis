""" Distributions module

:Author: Thomas Cokelaer <Thomas.Cokelaer@inria.fr>
 """
__version__ = "$Id$"

import interface
import error

import _stat_tool

from _stat_tool import _DiscreteParametricModel
from _stat_tool import _DistributionData
from _stat_tool import _Distribution

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
                        [[str, _DistributionData, _DiscreteParametricModel]],
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


def Binomial(inf_bound, sup_bound=_stat_tool.I_DEFAULT, \
             proba=_stat_tool.D_DEFAULT):
    """
    Construction of a binomial distribution

    :Parameters:
      * `inf_bound` (int) : lower bound to the range of possible values
        (shift parameter)
      * `sup_bound` (int) : upper bound to the range of possilbe values
      * `proba` (int, float) : probability of 'success'
    """
    # todo: seg fault when passing -1 as first arguments if there
    # is no assert here below
    # memory leak ?
    # todo:  returns error if ((inf_bound < min_inf_bound) ||
    # (inf_bound > MAX_INF_BOUND)) {

    error.CheckType([inf_bound, sup_bound, proba], [int, int, [int, float]])
    assert inf_bound >= 0
    assert inf_bound < sup_bound
    assert (sup_bound - inf_bound) <= _stat_tool.MAX_DIFF_BOUND
    assert proba <= 1. and proba > 0

    param = _stat_tool.D_DEFAULT

    return _DiscreteParametricModel(_stat_tool.BINOMIAL, \
        inf_bound, sup_bound, param, proba)


def Poisson(inf_bound, param=_stat_tool.D_DEFAULT):
    """
    Construction of a poisson distribution

    :Parameters:
      * `inf_bound` (int) : lower bound to the range of possible values
        (shift parameter)
      * `param` (int, float) : parameter of the Poisson distribution
    """

    assert inf_bound >= 0
    assert param > inf_bound
    assert param > 0. and param < _stat_tool.MAX_MEAN

    sup_bound = _stat_tool.I_DEFAULT
    proba = _stat_tool.D_DEFAULT

    return _DiscreteParametricModel(_stat_tool.POISSON, \
        inf_bound, sup_bound, param, proba)


def NegativeBinomial(inf_bound, param=_stat_tool.D_DEFAULT, \
                     proba=_stat_tool.D_DEFAULT):
    """
    Construction of a negative binomial distribution

    :Parameters:
      * inf_bound (int) : lower bound to the range of possible values
        (shift parameter)
      * param (int, float) : parameter of the Poisson distribution
      * proba (int, float) : probability of 'success'
    """

    # Check parameters
    assert inf_bound >= 0
    assert param > 0
    assert proba <= 1. and proba > 0
    assert (param * (1. - proba) / proba) <= _stat_tool.MAX_MEAN

    sup_bound = _stat_tool.I_DEFAULT

    return _DiscreteParametricModel(_stat_tool.NEGATIVE_BINOMIAL, \
        inf_bound, sup_bound, param, proba)


def Uniform(inf_bound, sup_bound=_stat_tool.I_DEFAULT):
    """
    Construction of a uniform distribution

    :Parameters:
      * inf_bound (int) : lower bound to the range of possible values
        (shift parameter)
      * sup_bound (int) : upper bound to the range of possilbe values
    """

    # Check parameters
    #todo:  returns error if ((inf_bound < min_inf_bound) ||
    # (inf_bound > MAX_INF_BOUND)) {
    assert inf_bound >= 0
    assert sup_bound >= 0
    assert inf_bound <= sup_bound
    assert (sup_bound - inf_bound) < _stat_tool.MAX_DIFF_BOUND

    param = _stat_tool.D_DEFAULT
    proba = _stat_tool.D_DEFAULT
    cumul_threshold = _stat_tool.CUMUL_THRESHOLD
    return _DiscreteParametricModel(_stat_tool.UNIFORM, \
        inf_bound, sup_bound, param, proba, cumul_threshold)


def Multinomial():
    """to be done"""
    raise NotImplementedError("Multinomial not yet implemented")

# Extend _DiscreteParametricModel
interface.extend_class( _DiscreteParametricModel, interface.StatInterface)

# Cast Functions


def ToDistribution(histo):
    """ Cast an object of type `_DistributionData` into an object
    of type `_Distribution`.

    :Parameters:
      * `histo` (DistributionData)

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
    type `_DistributionData`.

    :Parameters:
      * dist (distribution).

    :Returns:
        If the object dist contains a 'data' part, an object of
        type `_DistributionData` is returned, otherwise no object
        is returned.

    :Examples:

    .. doctest::
        :options: +SKIP

        >>> ToHistogram(dist)

    .. seealso::
        :func:`~openalea.stat_tool.distribution.ToDistribution`
    """
    return dist.extract_data()



