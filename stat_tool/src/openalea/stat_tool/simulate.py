#!/usr/bin/env python
#-*- coding: utf-8 -*-
"""Simulate functions

.. topic:: simulate.py summary

    A module dedicated to simulate functionalities

    :Code status: mature
    :Documentation status: to be completed
    :Author: Thomas Cokelaer <Thomas.Cokelaer@sophia.inria.fr>

"""

__all__ = ["Simulate"]


from . import error
from .enums import *


def Simulate(obj, *args):
    """Generation of a random sample from a distribution.

    :Parameters:

      * `dist` (distribution),
      * `mixt` (mixture)
      * `convol` (convolution)
      * `compound` (compound),
      * `size` (int): sample size.

    :Returns:

        If the first argument is of type distribution and if 0 < size < 1000000,
        an object of type HISTOGRAM is returned, otherwise no object is returned.
        If the first argument is of type mixture and if 0 < size < 1000000, an
        object of type mixture_data is returned, otherwise no object is returned.
        If the first argument is of type convolution and if 0 < size < 1000000, an
        object of type convolution_data is returned, otherwise no object is returned.
        If the first argument is of type compound and if 0 < size < 1000000, an
        object of type compound_data is returned, otherwise no object is returned.
        The returned object of type HISTOGRAM, mixture_data, convolution_data or
        compound_data contains both the simulated sample and the model used for
        simulation.

    :Example:

    .. doctest::
        :options: +SKIP

        >>> Simulate(dist, size)
        >>> Simulate(mixt, size)
        >>> Simulate(convol, size)
        >>> Simulate(compound, size)

    :See Also:

        Distribution,
        Mixture,
        Convolution,
        Compound,
        ExtractHistogram.
    """
    error.CheckArgumentsLength(args, 1, 1)
    return obj.simulate(args[0])
