#!/usr/bin/env python
#-*- coding: utf-8 -*-
"""Data transformation functions

.. topic:: data_transform.py summary

    A module dedicated to Output functions (plot, displya, save)

    :Code status: mature
    :Documentation status: to be completed
    :Author: Samuel Dufour-Kowalski <samuel.dufour@sophia.inria.fr>
        Thomas Cokelaer <Thomas.Cokelaer@sophia.inria.fr>

    :Revision: $Id: data_transform.py 15183 2013-11-06 10:35:50Z jbdurand $

.. testsetup:: *

    from openalea.stat_tool import *
    
.. warning:: sequence analysis package also contain a data transform module
"""
__revision__ = "$Id: data_transform.py 15183 2013-11-06 10:35:50Z jbdurand $"

import error
from openalea.stat_tool._stat_tool import \
     _DiscreteMixture, \
     _DiscreteMixtureData, \
     _Compound, \
     _CompoundData, \
     _Convolution, \
     _ConvolutionData,\
     _Vectors

from enums import keep_type



def Merge(obj, *args):
    """
    Merging of objects of the same 'data' type or merging of sample correlation functions.

    :Parameters:

      * histo1, histo2, ... (_Histogram, _DiscreteMixtureData, _ConvolutionData, _CompoundData),
      * vec1, vec2, ... (_Vectors),
      * timev1, timev2, ... (_TimeEvents, _RenewalData),
      * seq1, seq2, ... (_Sequences, _DiscreteSequences, _MarkovData, _SemiMarkovData),
      * top1, top2, ... (_Tops),
      * correl1, correl2, ... (_Correlation).

    :Returns:

        If the arguments are of type _Histogram, _DiscreteMixtureData, _ConvolutionData, _CompoundData
        an object of type _Histogram is returned.

        If the arguments are of type _Vectors and if the vectors have the same number of variables,
        an object of type vectors is returned, otherwise no object is returned.

        If the arguments are of type _TimeEvents, _RenewalData, an object of type
        _TimeEvents is returned.

        If the arguments are of type _Sequences, _DiscreteSequences, _MarkovData, _SemiMarkovData
        and if the sequences have the same number of variables, an object of type _Sequences
        is returned.

        If the arguments are of type _Tops, an object of type _Tops is returned.
        If the arguments are of type correlation, an object of type correlation is returned.

    :Examples:

    .. doctest::
        :options: +SKIP

        >>> Merge(histo1, histo2,...)
        >>> Merge(vec1, vec2,...)
        >>> Merge(timev1, timev2,...)
        >>> Merge(seq1, seq2,...)
        >>> Merge(discrete_seq1, discrete_seq2,...)
        >>> Merge(top1, top2,...)
        >>> Merge(correl1, correl2,...)

    .. seealso::
        :func:`~openalea.stat_tool.cluster.Cluster`,
        :func:`~openalea.stat_tool.data_transform.Shift`,
        :func:`~openalea.stat_tool.cluster.Transcode`,
        :func:`~openalea.stat_tool.data_transform.ValueSelect`,
        :func:`~openalea.stat_tool.data_transform.MergeVariable`,
        :func:`~openalea.stat_tool.data_transform.SelectIndividual`,
        :func:`~openalea.stat_tool.data_transform.SelectVariable`,
        :func:`~openalea.stat_tool.cluster.NbEventSelect`,
        :func:`~openalea.stat_tool.cluster.TimeScaling`,
        :func:`~openalea.stat_tool.cluster.TimeSelect`,
        :func:`~openalea.stat_tool.cluster.AddAbsorbingRun`,
        :func:`~openalea.stat_tool.cluster.Cumulate`,
        :func:`~openalea.stat_tool.cluster.Difference`,
        :func:`~openalea.stat_tool.cluster.IndexExtract`,
        :func:`~openalea.stat_tool.cluster.LengthSelect`,
        :func:`~openalea.stat_tool.cluster.MovingAverage`,
        :func:`~openalea.stat_tool.cluster.RecurrenceTimeSequences`,
        :func:`~openalea.stat_tool.cluster.RemoveRun`,
        :func:`~openalea.stat_tool.cluster.Reverse`,
        :func:`~openalea.stat_tool.cluster.SegmentationExtract`,
        :func:`~openalea.stat_tool.cluster.VariableScaling`,
        :func:`~openalea.stat_tool.cluster.RemoveApicalInternodes`
    """

    ret =  obj.merge(list(args))
    if hasattr(ret, 'markovian_sequences'):
        try:
            newret = ret.markovian_sequences()
            return newret
        except:
            import warnings
            warnings.warn("""markovian_sequences did not work. Check that this
             is normal behaviour""")

    return ret


def MergeVariable(obj, *args, **kargs):
    """
    Merging of variables.

    :Parameters:

        * vec1, vec2, ... (_Vectors),
        * seq1, seq2, ... (_Sequences, _DiscreteSequences, _MarkovData, _SemiMarkovData).

    :Keywords:

      * RefSample (int): reference sample to define individual identifiers
        (the default: no reference sample).

    :Returns:

        If the arguments are of type _Vectors and if the number of vectors is the same
        for each sample, an object of type _Vectors is returned.

        If the arguments are of type _Sequences, _DiscreteSequences, _MarkovData, _SemiMarkovData,
        if all the variables are of type STATE, and if the number and the lengths of sequences
        are the same for each sample, an object of type _Sequences or _DiscreteSequences is returned.

        The returned object is of type _DiscreteSequences if all the variables are of type STATE,
        if the possible values for each variable are consecutive from 0 and if the number of
        possible values for each variable is < 15.

    :Examples:

    .. doctest::
        :options: +SKIP

        >>> MergeVariable(histo1, histo2)
        >>> MergeVariable(vec1, vec2,..., RefSample=2)
        >>> MergeVariable(seq1, seq2,..., RefSample=2)

    .. seealso::
        :func:`~openalea.stat_tool.cluster.Cluster`,
        :func:`~openalea.stat_tool.data_transform.Shift`,
        :func:`~openalea.stat_tool.cluster.Transcode`,
        :func:`~openalea.stat_tool.data_transform.ValueSelect`,
        :func:`~openalea.stat_tool.data_transform.Merge`,
        :func:`~openalea.stat_tool.data_transform.SelectIndividual`,
        :func:`~openalea.stat_tool.data_transform.SelectVariable`,
        :func:`~openalea.stat_tool.cluster.AddAbsorbingRun`,
        :func:`~openalea.stat_tool.cluster.Cumulate`,
        :func:`~openalea.stat_tool.cluster.Difference`,
        :func:`~openalea.stat_tool.cluster.IndexExtract`,
        :func:`~openalea.stat_tool.cluster.LengthSelect`,
        :func:`~openalea.stat_tool.cluster.MovingAverage`,
        :func:`~openalea.stat_tool.cluster.RecurrenceTimeSequences`,
        :func:`~openalea.stat_tool.cluster.RemoveRun`,
        :func:`~openalea.stat_tool.cluster.Reverse`,
        :func:`~openalea.stat_tool.cluster.SegmentationExtract`,
        :func:`~openalea.stat_tool.cluster.VariableScaling`,
    """
    #todo:manage the marjkovian_sequences conversion if input
    # is made os Sequences
    arg1 = args[0]
    for arg in args:
        error.CheckType( [arg], [type(arg1)])

    RefSample = kargs.get("RefSample", -1)
    error.CheckType([RefSample],  [int])

    return obj.merge_variable(list(args), RefSample)


def ExtractData(model):
    """Extraction of the 'data' part of an object of type 'model'.

    This function enables to extract the 'data' part of an object of type 'model'
    when the estimation of model parameters from data gives rise to the construction
    of pseudo-data. This situation is notably exemplified by the computation of
    optimal state sequences from estimated hidden Markovian processes (optional
    argument StateSequences of the function Estimate set at "ForwardBackward" or
    "Viterbi").

    :Parameters:
      * mixt (_DiscreteMixture),
      * convol (_Convolution),
      * compound (_Compound),
      * hmc (_HiddenMarkov),
      * hsmc (_HiddenSemiMarkov).

    :Returns:
      - If mixt contains a 'data' part, an object of type `_DiscreteMixtureData` is returned.
      - If convol contains a 'data' part, an object of type `_ConvolutionData` is returned.
      - If compound contains a 'data' part, an object of type `_CompoundData` is returned.
      - If hmc contains a 'data' part, an object of type `_MarkovData` is returned.
      - If hsmc contains a 'data' part, an object of type `_SemiMarkovData` is returned.

    :Examples:

    .. doctest::
        :options: +SKIP

        >>> ExtractData(mixt)
        >>> ExtractData(convol)
        >>> ExtractData(compound)
        >>> ExtractData(hmc)
        >>> ExtractData(hsmc)

    .. seealso::
        :func:`~openalea.stat_tool.estimate.Estimate`
    """

    return model.extract_data()


def ExtractDistribution(model, *args):
    """
    Extraction of a distribution from an object of type 'model'.

    :Parameters:

      * mixt (_DiscreteMixture),
      * mixt_type (string): type of distribution: "Weight" or "Mixture",
      * index (int): index of the elementary distribution,
      * convol (_Convolution),
      * compound (_Compound),
      * compound_type (string): type of distribution: "Sum", "Elementary" or "Compound",
      * renew (renewal),
      * renew_type (string): type of distribution: "InterEvent", "Backward", "Forward", "LengthBias" or "Mixture",
      * time (int): observation period,
      * markov (markov, semi-markov, hidden_markov, hidden_semi-markov),
      * markov_type (string): type of distribution: "Observation", "FirstOccurrence", "Recurrence", "Sojourn", "NbRun" or "NbOccurrence",
      * state (int): state,
      * variable (int): variable index,
      * output (int): output,
      * top_param (top_parameters),
      * position (int): position.

    :Returns:

        If the arguments (mixt_type, index, compound_type, renew_type, time,
        markov_type, state, variable, output, position) defined an existing
        distribution, an object of type `_Distribution` is returned.

    :Examples:

    .. doctest::
        :options: +SKIP

        >>> ExtractDistribution(mixt, mixt_type)
        >>> ExtractDistribution(mixt, "Component", index)
        >>> ExtractDistribution(convol, "Elementary", index)
        >>> ExtractDistribution(convol, "Convolution")
        >>> ExtractDistribution(compound, compound_type)
        >>> ExtractDistribution(renew, renew_type)
        >>> ExtractDistribution(renew, "NbEvent", time)
        >>> ExtractDistribution(markov, markov_type, state)
        >>> ExtractDistribution(markov, markov_type, variable, output)
        >>> ExtractDistribution(top_param, position)

    .. seealso::
        :func:`~openalea.stat_tool.output.Plot`,
        :func:`~openalea.stat_tool.data_transform.Fit`,
        :func:`~openalea.stat_tool.simulate.Simulate`.
    """

    return Extract(model, *args)


def ExtractHistogram(data, *args, **kargs):
    """
    Extraction of a frequency distribution from an object of type 'data'.

    :Parameters:

      * mixt_histo (_DiscreteMixtureData),
      * mixt_type (string): type of distribution: "Weight" or "Mixture",
      * index (int): index of the elementary distribution,
      * convol_histo (_ConvolutionData),
      * compound_histo (_CompoundData),
      * compound_type (string): type of distribution: "Sum", "Elementary" or "Compound",
      * vec1 (_Vectors) : values,
      * vecn (_Vectors) : vectors,
      * variable (int) : variable index
      * timev (_TimeEvents, _RenewalData)
      * timev_type (string):
      * time  (int)  : observation period

    :Returns:

        If the arguments (mixt_type, index, compound_type, renew_type, time,
        markov_type, state, variable, output, position) defined an existing frequency
        distribution, an object of type `_Histogram` is returned.

    :Examples:

    .. doctest::
        :options: +SKIP

        >>> ExtractHistogram(mixt_histo, mixt_type)
        >>> ExtractHistogram(mixt_histo, "Component", index)
        >>> ExtractHistogram(mixt_histo, "Mixture")
        >>> ExtractHistogram(mixt_histo, "Weight")
        >>> ExtractHistogram(convol_histo, "Elementary", index)
        >>> ExtractHistogram(convol_histo, "Convolution")
        >>> ExtractHistogram(compound_histo, compound_type)
        >>> ExtractHistogram(compound_histo, "Sum")
        >>> ExtractHistogram(compound_histo, "Elementary")
        >>> ExtractHistogram(vec1)
        >>> ExtractHistogram(vecn, variable)
        >>> ExtractHistogram(renewval_data, renew_type)
        >>> ExtractHistogram(timev, timev_type)
        >>> ExtractHistogram(timev, "NbEvent", time)
        >>> ExtractHistogram(seq, "Length")
        >>> ExtractHistogram(seq, "Value")
        >>> ExtractHistogram(seq, "Value", variable)
        >>> ExtractHistogram(discrete_seq1, seq_type, value)
        >>> ExtractHistogram(discrete_seqn, seq_type, variable, value)
        >>> ExtractHistogram(simul_seq1, "Observation", value)
        >>> ExtractHistogram(simul_seq1, "Observation", variable, value)
        >>> ExtractHistogram(tops, "Main")
        >>> ExtractHistogram(top, "NbAxillary", position)


    .. seealso::
        :func:`~openalea.stat_tool.output.Plot`,
        :func:`~openalea.stat_tool.data_transform.Fit`,
        :func:`~openalea.stat_tool.simulate.Simulate`.
    """
    return Extract(data, *args, **kargs)


def Extract(obj, *args, **kargs):
    """
    Common method to redirect extract function call
    See`ExtractHistogram` or `ExtractDistribution`
    """

    ret = None

    if type(obj) in [_DiscreteMixture, _DiscreteMixtureData]:
        assert len(args)>=1
        error.CheckType([args[0]], [str])
        if args[0] == 'Mixture':
            assert len(args) == 1
            ret = obj.extract_DiscreteMixture()
        elif args[0] == 'Component':
            assert len(args) == 2
            error.CheckType([args[1]], [int])
            ret = obj.extract_component(args[1])
        elif args[0] == 'Weight':
            assert len(args) == 1
            ret = obj.extract_weight()
        else:
            raise ValueError("Excepted Component, Weight or Mixture")
    elif type(obj) in [_Convolution, _ConvolutionData]:
        assert len(args)>=1
        error.CheckType([args[0]], [[str, int]])
        if args[0] == 'Elementary' or isinstance(args[0], int):
            if len(args) == 1:
                error.CheckType([args[0]], [int])
                ret = obj.extract_elementary(args[0])
            elif len(args)==2:
                error.CheckType([args[0], args[1]], [str, int])
                ret = obj.extract_elementary(args[1])
        elif args[0] == 'Convolution':
            error.CheckType([args[0]], [[str, int]])
            ret = obj.extract_convolution()
        else:
            raise ValueError("Excepted \"Elementaty\", or index")

    elif type(obj) in [_Compound, _CompoundData]:
        assert len(args) == 1
        if args[0] == 'Sum':
            ret = obj.extract_sum()
        elif args[0] == 'Elementary':
            ret = obj.extract_elementary()
        elif args[0] == 'Compound':
            ret = obj.extract_compound()
        else:
            raise ValueError("Excepted Sum, Elementary or Compound")

    elif isinstance(obj, _Vectors):
        # _Vectors with one variable

        try:
            nb_var = obj.nb_variable
            if (nb_var>1):
                try:
                    variable = args[0]
                except IndexError:
                    raise TypeError("""Extract with vectors object need 1
                     arguments (variable) if nb variable>1""")
            else:
                variable = 1

            return obj.extract(variable)

        except AttributeError:
            raise ValueError("unknown issue while extracting vectors")
    else:
        # related to Top, Renewal, Markov , ...
        try:
            from openalea.sequence_analysis.data_transform import Extract \
                as newExtract
            ret = newExtract(obj, *args, **kargs)
        except ValueError:
            pass

    return ret


def SelectVariable(obj, variables, Mode="Keep"):
    """
    Selection of variables.

    :Parameters:

      * vec (vectors),
      * seq (sequences, discrete_sequences, markov_data, semi-markov_data),
      * variable (int): variable index.
      * variables (array(int)): variable indices.

    :Keywords:

      * Mode (string): conservation or rejection of the selected variables: "Keep" (default) or "Reject".

    :Returns:

      If either variable or variables[1], ..., variables[n] are valid indices of variables,
      an object of type vectors (respectively sequences or discrete_sequences) is returned,
      otherwise no object is returned. In the case of a first argument of type sequences,
      the returned object is of type discrete_sequences if all the variables are of type STATE,
      if the possible values for each variable are consecutive from 0 and if the number of
      possible values for each variable is < 15.

    :Examples:

    .. doctest::
        :options: +SKIP

        >>> SelectVariable(vec, variable, Mode="Reject")
        >>> SelectVariable(vec, variables, Mode="Reject")
        >>> SelectVariable(seq, variable, Mode="Reject")
        >>> SelectVariable(seq, variables, Mode="Reject")

    .. seealso::
        `AddAbsorbingRun`,
        :func:`~openalea.stat_tool.cluster.Cluster`,
        :func:`~openalea.stat_tool.cumulate.Cumulate`,
        `Difference`,
        `IndexExtract`,
        `LengthSelect`,
        :func:`~openalea.stat_tool.data_transform.Merge`,
        :func:`~openalea.stat_tool.data_transform.MergeVariable`,
        `MovingAverage`,
        `RecurrenceTimeSequences`,
        `RemoveRun`,
        `Reverse`,
        :func:`~openalea.stat_tool.data_transform.SelectIndividual`,
        :func:`~openalea.stat_tool.data_transform.Shift`,
        :func:`~openalea.stat_tool.cluster.Transcode`,
        :func:`~openalea.stat_tool.data_transform.ValueSelect`,
        `SegmentationExtract`,
        `VariableScaling`.
    """

    error.CheckType([variables, Mode] , [[int, list], str])
    #todo: check that Mode is in ["Keep", "Reject"]

    keep = bool(Mode == "Keep" or Mode == "keep")

    if isinstance(variables, int):
        variables = [variables]

    ret = obj.select_variable(variables, keep)

    return ret

def SelectIndividual(obj, identifiers, Mode="Keep"):
    """
    Selection of vectors, sequences, tops or patterns (in a dissimilarity matrix).

    :Parameters:

      * vec (vectors),
      * seq (sequences, discrete_sequences, markov_data, semi-markov_data),
      * top (tops),
      * dist_matrix (distance_matrix),
      * identifiers (array(int)): identifiers.

    :Keywords:

       Mode (string): conservation or rejection of the selected individuals: "Keep" (default) or "Reject".

    :Returns:

        If identifiers[1], ..., identifiers[n] are valid identifiers of vectors (respectively
        sequences, tops or patterns compared in a dissimilarity matrix), an object of type vectors
        (respectively sequences or discrete_sequences, tops or distance_matrix) is returned,
        otherwise no object is returned. In the case of a first argument of type sequences,
        discrete_sequences, markov_data, semi-markov_data, the returned object is of type
        discrete_sequences if all the variables are of type STATE, if the possible values for
        each variable are consecutive from 0 and if the number of possible values for each variable
        is < 15.

    :Examples:

    .. doctest::
        :options: +SKIP

        >>> SelectIndividual(vec, identifiers, Mode="Reject")
        >>> SelectIndividual(seq, identifiers, Mode="Reject")
        >>> SelectIndividual(top, identifiers, Mode="Reject")
        >>> SelectIndividual(dist_matrix, identifiers, Mode="Reject")

    .. seealso::
        :func:`~openalea.stat_tool.cluster.Cluster`,
        :func:`~openalea.stat_tool.data_transform.Merge`,
        :func:`~openalea.stat_tool.data_transform.Shift`,
        :func:`~openalea.stat_tool.cluster.Transcode`,
        :func:`~openalea.stat_tool.data_transform.ValueSelect`,
        :func:`~openalea.stat_tool.data_transform.MergeVariable`,
        :func:`~openalea.stat_tool.data_transform.SelectVariable`
        `AddAbsorbingRun`,
        `Cumulate`,
        `Difference`,
        `IndexExtract`,
        `LengthSelect`,
        `MovingAverage`,
        `RecurrenceTimeSequences`,
        `RemoveSeries`,
        `Reverse`,
        `SegmentationExtract`,
        `VariableScaling`,
        `RemoveApicalInternodes`,
        `Symmetrize`.

    """
    error.CheckType([identifiers, Mode] , [list, str])

    #todo: CHECK THAT Mode is in ["Keep", "Reject"]
    keep = bool(Mode == "Keep" or Mode == "keep")

    ret = None
    try:
        ret = obj.select_individual(identifiers, keep)
    except:
        raise Exception("Could not run extract_data on the input variable. ")

    #if ret:
    #    try:
            # if obj is a sequence, returns markovian_sequences
    #        return ret.markovian_sequences()
    #    except AttributeError:
    #        return ret
    #else:
    #    raise Exception("Must not enter here")
    # the code above prevent tests to succeed.
    return ret


def ValueSelect(obj, *args, **kargs):
    """Selection of individuals according to the values taken by a variable

    :Parameters:

      * histo (histogram, mixture_data, convolution_data, compound_data),
      * value (int): value,
      * min_value (int): minimum value,
      * max_value (int): maximum value,
      * vec1 (vectors): values,
      * vecn (vectors): vectors,
      * variable (int): variable index,
      * seq1 (sequences, discrete_sequences, markov_data, semi-markov_data): univariate sequences,
      * seqn (sequences, discrete_sequences, markov_data, semi-markov_data): multivariate sequences.

    :Keywords:

      * Mode (string): conservation or rejection of selected individuals: "Keep" (the default) or "Reject".

    :Returns:

        If either value 0 or if 0 < min_value < max_value and if the range of values defined either
        by value or by min_value and max_value enables to select individuals, an object of
        type HISTOGRAM is returned (respectively vectors, sequences or discrete_sequences),
        otherwise no object is returned. In the case of a first argument of type sequences,
        discrete_sequences, markov_data or semi-markov_data, the returned object is of type
        discrete_sequences if all the variables are of type STATE, if the possible values for each
        variable are consecutive from 0 and if the number of possible values for each variable is < 15.

    :Examples:

    .. doctest::
        :options: +SKIP

        >>> ValueSelect(histo, value, Mode="Reject")
        >>> ValueSelect(histo, min_value, max_value, Mode="Reject")
        >>> ValueSelect(vec1, value, Mode="Reject")
        >>> ValueSelect(vec1, min_value, max_value, Mode="Reject")
        >>> ValueSelect(vecn, variable, value, Mode="Reject")
        >>> ValueSelect(vecn, variable, min_value, max_value, Mode="Reject")
        >>> ValueSelect(seq1, value, Mode="Reject")
        >>> ValueSelect(seq1, min_value, max_value, Mode="Reject")
        >>> ValueSelect(seqn, variable, value, Mode="Reject")
        >>> ValueSelect(seqn, variable, min_value, max_value, Mode="Reject")

    .. seealso::
        :func:`~openalea.stat_tool.cluster.Cluster`,
        :func:`~openalea.stat_tool.data_transform.Merge`,
        :func:`~openalea.stat_tool.data_transform.Shift`,
        :func:`~openalea.stat_tool.data_transform.Transcode`,
        :func:`~openalea.stat_tool.data_transform.SelectIndividual`,
        :func:`~openalea.stat_tool.data_transform.MergeVariable`,
        :func:`~openalea.stat_tool.data_transform.SelectVariable`
        Cumulate`
        Difference`
        IndexExtract`
        LengthSelect`,
        MovingAverage`,
        RecurrenceTimeSequences`
        RemoveRun`,
        Reverse`,
        SegmentationExtract`,
        VariableScaling`.
    """
    error.CheckArgumentsLength(args, 1, 3)
    Mode = error.ParseKargs(kargs, "Mode", "Keep", keep_type)
    #keep = bool(Mode == "Keep" or Mode == "keep")
    keep = bool(Mode == "Keep")
    # Test for vectors
    try:
        nb_variable = obj.nb_variable
    except AttributeError:
        nb_variable = 0


    if len(args) == 3 :
        variable, umin, umax = args

    elif len(args) == 2:
        # 2 cases (min_value, max_value) or (variable, value)
        if nb_variable:
            variable, umin = args
            umax = umin
        else:
            umin, umax = args

    elif len(args) == 1:
        value = args[0]
        error.CheckType([value], [[int, tuple, list]])
        if isinstance(value, tuple) and len(value) == 2:
            umin, umax = value
        else:
            umin = umax = value

    if(nb_variable):    # Vectors, sequences
        return obj.value_select(variable, umin, umax, keep)
    else:
        return obj.value_select(umin, umax, keep)


def Shift(obj, *args):
    """
    Shifting of values

    :Parameters:

      * histo (histogram, mixture_data, convolution_data, compound_data),
      * param (int): shifting parameter,
      * vec1 (vectors): values,
      * vecn (vectors): vectors,
      * variable (int): variable index,
      * seq1 (sequences): univariate sequences,
      * seqn (sequences): multivariate sequences.

    :Returns:

      If the shifting makes that the lower bound to the possible values is positive, an
      object of type HISTOGRAM (respectively _Vectors, _Sequences) is returned. In the
      case of a first argument of type sequences, the returned object is of type
      discrete_sequences if all the variables are of type STATE, if the possible values
      for each variable are consecutive from 0 and if the number of possible values for
      each variable is 15.

    :Examples:

    .. doctest::
        :options: +SKIP

        >>> Shift(histo, param)
        >>> Shift(vec1, param)
        >>> Shift(vecn, variable, param)
        >>> Shift(seq1, param)
        >>> Shift(seqn, variable, param)

    .. seealso::
        :func:`~openalea.stat_tool.cluster.Cluster`,
        :func:`~openalea.stat_tool.data_transform.Merge`,
        :func:`~openalea.stat_tool.data_transform.Transcode`,
        :func:`~openalea.stat_tool.data_transform.SelectIndividual`,
        :func:`~openalea.stat_tool.data_transform.MergeVariable`,
        :func:`~openalea.stat_tool.data_transform.SelectVariable`
        :func:`~openalea.stat_tool.data_transform.AddAbsorbingRun`,
        :func:`~openalea.stat_tool.data_transform.Cumulate`,
        :func:`~openalea.stat_tool.data_transform.Difference`,
        :func:`~openalea.stat_tool.data_transform.Lengthselect`,
        :func:`~openalea.stat_tool.data_transform.MovingAverage`,
        :func:`~openalea.stat_tool.data_transform.IndexExtract`,
        :func:`~openalea.stat_tool.data_transform.RecurrenceTimeSequences`,
        :func:`~openalea.stat_tool.data_transform.RemoveRun`,
        :func:`~openalea.stat_tool.data_transform.Reverse`,
        :func:`~openalea.stat_tool.data_transform.SegmentationExtract`,
        :func:`~openalea.stat_tool.data_transform.ValueSelect`,
        :func:`~openalea.stat_tool.data_transform.VariableScaling`.
    """
    error.CheckArgumentsLength(args, 1, 3)

    try:
        nb_variable = obj.nb_variable
    except AttributeError:
        nb_variable = 0

    if nb_variable == 1:
        param = args[0]
        ret = obj.shift(1, param)
    elif nb_variable > 1:
        variable = args[0]
        param = args[1]
        ret = obj.shift(variable, param)
    else:
        param = args[0]
        ret = obj.shift(param)

    return ret


def Fit(histo, dist):
    r"""
    Fit of a frequency distribution by a theoretical distribution.

    The result is displayed in the shell window (characteristics of the frequency
    and theoretical distributions, log-likelihood of the data for the theoretical distribution,
    information - maximum possible log-likelihood of the data -, c2 goodness of fit test).

    The difference between the information measure and the log-likelihood is the
    Kullback-Leibler divergence from the observed distribution to the theoretical distribution.
    It is also one-half the deviance of the theoretical distribution.

    Assume that a sample of size n is generated by a given random variable.
    The statistic measures the random deviation between the observed frequencies fi and the
    theoretical frequencies npi:

    .. math::

        D^2 = \sum_{i=0}^k \frac{\left(f_i - n p_i\right)^2}{n p_i}
        \textrm{with}  \sum_{i=0}^k f_i = n

    If each theoretical frequency npi is greater than a given threshold (between 1 and 5
    according to the authors), has a c2 with k - 1 degrees of freedom.

    :Parameters:

      * histo (histogram, mixture_data, convolution_data, compound_data),
      * dist (distribution, mixture, convolution, compound).

    :Returns:

       Distribution


    :Examples:

    .. doctest::
        :options: +SKIP

        >>> Fit(histo, dist)

    .. todo:: documenation: get back the latex equations

    """

    return histo.fit(dist)



def Unnormalize(obj):
    """Unnormalize

    :Parameters:

        * dist_matrix (distance_matrix).

     :Returns:

        An object of type distance_matrix is returned.

    :Examples:

    .. doctest::
        :options: +SKIP

        >>>  Unnormalize(dist_matrix)

    """
    return obj.unnormalize()







def Symmetrize(obj):
    """

    .. todo:: documentation
    """
    ret = None
    try:
        ret = obj.symmetrize()
    except ValueError:
        pass
    if ret:
        return ret

def TruncateDistribution(obj, variable):
    """
    .. todo:: to be tests

    """
    return obj.truncate(variable)



def SelectStep(obj, *args):
    """Change the internal step of a vector or a sequence

    :param obj: the vector or sequence objet
    :param argument 1: the new step

    :Example:

    .. doctest::
        :options: +SKIP

        >>> seq = Sequences([])
        >>> SelectStep(seq, 100)
        >>> Plot(seq)

    .. todo:: shall we move this function to sequence_analysis package?

    """
    error.CheckArgumentsLength(args, 1, 2)

    try:
        nb_variable = obj.nb_variable
    except AttributeError:
        raise TypeError("object has no nb_variable. Check that it is a Vector or Sequence")


    if len(args)==2:
        variable, step = args
        error.CheckType([step], [[int, float]])
        error.CheckType([variable], [[int]])
    elif len(args)==1 and nb_variable==1:
        variable = 1
        step = args[0]
        error.CheckType([step], [[int, float]])
    else:
        if nb_variable!=1:
            raise SyntaxError("Wrong number of arguments. The number of variable is greater than 1 (%s) therefore you must provide a variable and a step like in SelectStep(object, 1, 100)" % nb_variable)
        else:
            raise ValueError("UnknownError")

    #obj.get_marginal_histogram(variable)
    ret = obj.select_step(variable, step)
    return ret
