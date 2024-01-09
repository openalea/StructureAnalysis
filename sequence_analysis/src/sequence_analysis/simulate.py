#!/usr/bin/env python
#-*- coding: utf-8 -*-
"""

.. topic:: simulate.py summary

    A module dedicated to Simulate functionalities

    :Code status: mature
    :Documentation status: to be completed
    :Author: Thomas Cokelaer <Thomas.Cokelaer@sophia.inria.fr>

    :Revision: $Id: simulate.py 15827 2014-02-19 16:17:45Z jbdurand $
    
"""
__version__ = "$Id: simulate.py 15827 2014-02-19 16:17:45Z jbdurand $"


from openalea.stat_tool.simulate import Simulate as SimulateDistribution

from openalea.stat_tool._stat_tool import _FrequencyDistribution
from openalea.stat_tool._stat_tool import _DiscreteMixtureData
from openalea.stat_tool._stat_tool import _ConvolutionData
from openalea.stat_tool._stat_tool import _CompoundData

from openalea.sequence_analysis._sequence_analysis import srand
from openalea.sequence_analysis._sequence_analysis import _Renewal
# from openalea.sequence_analysis._sequence_analysis import _TopParameters
from openalea.sequence_analysis._sequence_analysis import _TimeEvents
from openalea.sequence_analysis._sequence_analysis import \
    _VariableOrderMarkov,\
    _HiddenVariableOrderMarkov,\
    _SemiMarkov,\
    _NonHomogeneousMarkov,\
    _HiddenSemiMarkov
from openalea.sequence_analysis._sequence_analysis import \
    _MarkovianSequences,\
    _VariableOrderMarkovData,\
    _SemiMarkovData,\
    _NonHomogeneousMarkovData

from openalea.stat_tool import error
from enums import stochastic_process_type
from openalea.stat_tool.enums import bool_type


def Simulate(obj, *args, **kargs):
    """Simulate

    * Generation of a random sample from a distribution.
    * Generation of a random sample from a renewal process.
    * Generation of a sample of sequences from a (hidden) Markovian process.
    * Generation of 'tops' from 'top' parameters.


    :Usage:

    .. doctest::
        :options: +SKIP
        
        >>> Simulate(dist, size)
        >>> Simulate(mixt, size)
        >>> Simulate(convol, size)
        >>> Simulate(compound, size)
        >>> Simulate(renew, type, time_histo)
        >>> Simulate(renew, type, size, time)
        >>> Simulate(renew, type, size, timev)
        >>> Simulate(markov, length_histo)
        >>> Simulate(markov, size, length)
        >>> Simulate(markov, size, seq)
        >>> Simulate(semi_markov, length_histo, Counting=False)
        >>> Simulate(semi_markov, size, length, Counting=False)
        >>> Simulate(semi_markov, size, seq, Counting=False)
        >>> Simulate(top_param, size, nb_trial, NbAxillary=2)

    :Arguments (distribution case):

    * dist (distribution),
    * mixt (mixture),
    * convol (convolution),
    * compound (compound),
    * size (int): sample size.

    :Arguments (renewal case):

    * renew (renewal),
    * type (string): type of renewal process: "Ordinary" or "Equilibriun",
    * time_histo (histogram, mixture_data, convolution_data, compound_data): 
      frequency distribution of the length of the observation period,
    * size (int): sample size,
    * time (int): length of the observation period,
    * timev (time_events, renewal_data),

    :Arguments (markovian case):

    * markov (markov, hidden_markov),
    * semi_markov (semi-markov, hidden_semi-markov),
    * length_histo (histogram, mixture_data, convolution_data, compound_data):
      frequency distribution of sequence lengths,
    * size (int): sample size,
    * length (int): sequence length,
    * seq (discrete_sequences, markov_data, semi-markov_data),

    :Arguments (top case):

    * top_param (top_parameters),
    * size (int): sample size,
    * nb_trial (int): number of Bernoulli trials for the growth of the parent shoot,

    :Optional Arguments:

    * Counting (bool): computation of counting distributions (default value: True).
    * NbAxillary (int): number of offspring shoots per node (default value: 1,
      should be < 4).

    :Returned Object:

    * If the first argument is of type distribution and if 0 < size < 1000000, 
      an object of type HISTOGRAM is returned, otherwise no object is returned.
    * If the first argument is of type mixture and if 0 < size < 1000000, an 
      object of type mixture_data is returned, otherwise no object is returned.
    * If the first argument is of type convolution and if 0 < size < 1000000,
      an object of type convolution_data is returned, otherwise no object is returned.
    * If the first argument is of type compound and if 0 < size < 1000000, an 
      object of type compound_data is returned, otherwise no object is returned.
    * The returned object of type HISTOGRAM, mixture_data, convolution_data or
      compound_data contains both the simulated sample and the model used for 
      simulation.
    * If 0 < (sample size) < 1000000, if (minimum length of the observation
      period) > 2 and if (maximum length of the observation period) XX 1000,
      an object of type renewal_data is returned, otherwise no object is
      returned. The returned object contains both the simulated sample (not
      only count data but also time sequences) and the model used for simulation.
    * If the first argument is of type markov or hidden_markov, if 0 < (sample size)
      < 100000, if (minimum sequence length) < 2, if (maximum sequence length) <
      1000 and if (cumulative sequence length) < 1000000, an object of type
      markov_data is returned, otherwise no object is returned.
    * If the first argument is of type semi-markov or hidden_semi-markov, if
      0 < (sample size) < 100000, if (minimum sequence length) < 2, if (maximum
      sequence length) < 1000 and if (cumulative sequence length) < 1000000, an
      object of type semi-markov_data is returned, otherwise no object is returned.
    * The returned object contains both the simulated sequences and the model
      used for simulation.
    * If 0 < size < 100000 and if 0 < nb_trial < 1000, an object of type tops
      is returned, otherwise no object is returned. The returned object contains
      both the simulated 'tops' and the model used for simulation.


    :Description:

    If the fourth argument is of type time_events or renewal_data, the simulated
    sample has the same distribution of length of the observation period than 
    the original sample given by this fourth argument. This simulation mode is
    particularly useful to study the effects of length biasing.

    If the third argument is of type discrete_sequences, markov_data or 
    semi-markov_data, the simulated sequences has the same length distribution
    than the original sample given by this third argument. This simulation 
    mode is particularly useful to study the effects of length biasing.

    .. seealso::

        :func:`~openalea.stat_tool.distribution.Distribution`,
        :func:`~openalea.stat_tool.mixture.Mixture`,
        :func:`~openalea.stat_tool.convolution.Convolution`,
        :func:`~openalea.stat_tool.compound.Compound`,
        :func:`~openalea.stat_tool.data_transform.ExtractHistogram`.
        :func:`~openalea.sequence_analysis.semi_markov.SemiMarkov`,
        :func:`~openalea.sequence_analysis.hidden_semi_markov.HiddenSemiMarkov`,
        :func:`~openalea.sequence_analysis.renewal.Renewal`,
        :func:`~openalea.sequence_analysis.top_parameters.TopParameters`,
    """

    _valid_dists = [_FrequencyDistribution, _DiscreteMixtureData, \
                              _ConvolutionData, _CompoundData]

    # standard distribution case
    if len(args) == 1  and isinstance(args[0], int):
        return SimulateDistribution(obj, args[0])
    # top parameters case
    elif isinstance(obj, _TopParameters):
        error.CheckArgumentsLength(args, 2, 2)
        NbAxillary = kargs.get("NbAxillary", 1)
        error.CheckType([args[0], args[1], NbAxillary], [int, int, int])
        return obj.simulate(args[0], args[1], NbAxillary)

    # Renewal case
    elif isinstance(obj, _Renewal):
        # check that args[0] is a correct key
        Type = error.CheckDictKeys(args[0], stochastic_process_type)



        error.CheckType([args[1]], [[_FrequencyDistribution, _DiscreteMixtureData, \
                                     _ConvolutionData, _CompoundData, int]])

        if type(args[1]) in _valid_dists:
            ret = obj.simulation_histogram(Type, args[1])
        elif isinstance(args[1], int):
            error.CheckType([args[2]], [[int, _TimeEvents]])
            if isinstance(args[2], int):
                ret = obj.simulation_nb_elements(Type, args[1], args[2])
            else:
                ret =  obj.simulation_time_events(Type, args[1], args[2])
        return ret
    # other cases (Markovian, semi_markov, hidden_semi_markov and so on
    else:


        error.CheckType([obj], [[_VariableOrderMarkov, _SemiMarkov,
                                 _HiddenVariableOrderMarkov,
                                 _NonHomogeneousMarkov, _HiddenSemiMarkov]])

        Counting = error.ParseKargs(kargs, "Counting", True, bool_type)
        if type(obj) not in [_SemiMarkov, _HiddenSemiMarkov]:
            Counting = True

        #order of the if statements is important ! Keep it that way
        if isinstance(args[0], int) and isinstance(args[1], int):
            ret = obj.simulation_nb_elements(args[0], args[1], Counting)

        #here the second arguments is data structure such as Sequences
        elif isinstance(args[0], int):

            error.CheckType([args[1]], [[_MarkovianSequences,
                                         _VariableOrderMarkovData,
                                         _SemiMarkovData,
                                         _NonHomogeneousMarkovData]])
            ret = obj.simulation_markovian_sequences(args[0],
                                                      args[1], Counting)

        # first argument is a compound_data or equivalent
        else:
            error.CheckType([args[0]], [_valid_dists])
            ret =  obj.simulation_histogram(args[0], Counting, False)

        return ret
