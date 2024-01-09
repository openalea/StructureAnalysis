#!/usr/bin/env python
#-*- coding: utf-8 -*-
"""

.. topic:: data_transform.py summary

    A module dedicated to Variable Order Markov objects

    :Code status: mature
    :Documentation status: to be completed
    :Author: Thomas Cokelaer <Thomas.Cokelaer@sophia.inria.fr>

    :Revision: $Id: data_transform.py 15827 2014-02-19 16:17:45Z jbdurand $
    
"""
__version__ = "$Id: data_transform.py 15827 2014-02-19 16:17:45Z jbdurand $"

from openalea.stat_tool.error import *

from openalea.stat_tool._stat_tool import _Compound
from openalea.stat_tool._stat_tool import _DiscreteMixture
from openalea.stat_tool._stat_tool import _Convolution
from openalea.stat_tool._stat_tool import _DiscreteParametricModel
from openalea.stat_tool._stat_tool import _Vectors
from openalea.stat_tool._stat_tool import I_DEFAULT

from _sequence_analysis import _MarkovianSequences
from _sequence_analysis import _Sequences
from _sequence_analysis import _SemiMarkov
from _sequence_analysis import _SemiMarkovData
# from _sequence_analysis import _TopParameters
from _sequence_analysis import _VariableOrderMarkov
from _sequence_analysis import _VariableOrderMarkovData
from _sequence_analysis import _NonHomogeneousMarkov
from _sequence_analysis import _NonHomogeneousMarkovData
from _sequence_analysis import _Renewal
from _sequence_analysis import _HiddenVariableOrderMarkov
from _sequence_analysis import _HiddenSemiMarkov
# from _sequence_analysis import _Tops
from _sequence_analysis import _RenewalData
from _sequence_analysis import _TimeEvents

from openalea.stat_tool._stat_tool import _VectorDistance
from openalea.stat_tool import error
from openalea.stat_tool import I_DEFAULT

from openalea.sequence_analysis._sequence_analysis import MEAN_CHANGE
from openalea.sequence_analysis._sequence_analysis import MEAN_VARIANCE_CHANGE
from openalea.sequence_analysis._sequence_analysis import NB_EVENT
from openalea.stat_tool._stat_tool import MIN_PROBABILITY

from enums import markovian_sequence_type, output_map, histogram_type
from enums import mode_type, func_map, estimator_map, model_type, seq_map
from enums import renewal_nb_event_map, sub_func_map, nb_segment_map
from enums import output_type


from openalea.stat_tool.enums import keep_type, bool_type
from openalea.stat_tool.enums import format_type




def _check_nb_variable(obj, variable):
    """checks that nb_variable of an object is compatible with user variable


    :param obj: an object that must contain a nb_variable method
    :param variable: a variable

    :returns: variable

    if variable is None and obj.nb_variable =1 , returns 1
    else returns same value as input.
    """
    try:
        _dummy = obj.nb_variable
    except:
        raise TypeError("nb_variable method not found. Check your input")

    if obj.nb_variable == 1:
        if variable == None:
            variable = 1
        elif variable != 1:
            raise TypeError("only one variable available. Set variable to 1")

    if variable == None:
        raise TypeError("""variable argument must be set because number
        of variable is larger than 1""")

    if variable > obj.nb_variable:
        raise TypeError("variable larger than max number of variables %s" %
                        obj.nb_variable)

    return variable


def RemoveRun(obj, *args, **kargs):
    """RemoveRun

    Removal of the first or last run of a given value (for a given variable) in a sequence.

    :Usage:

    .. doctest::
        :options: +SKIP
        
        >>> RemoveRun(seq1, value, position, MaxLength=4)
        >>> RemoveRun(seqn, variable, value, position, MaxLength=4)

    :Arguments:

    * seq1 (sequences, discrete_sequences, markov_data, semi-markov_data): univariate sequences,
    * seqn (sequences, discrete_sequences, markov_data, semi-markov_data): multivariate sequences,
    * variable (int): variable index, (optional if only one variable
    * value (int): value,
    * position (string): position of the removed run: "Begin" (or "b") or "End" (or "e").

    :Optional Arguments:

    * MaxLength (int): maximum length of the removed run (default behaviour: the entire run is removed).

    :Returned Object:

    * If variable is a valid index of variable, if value is a possible value and
      if MaxLength > 0, an object of type sequences or discrete_sequences is
      returned, otherwise no object is returned. The returned object is of type
      discrete_sequences if all the variables are of type STATE, if the possible
      values for each variable are consecutive from 0 and if the number of
      possible values for each variable is <= 15.

    :Examples:

    .. doctest::
        :options: +SKIP
    
        >>> RemoveRun(seq1, 0, "End")
        >>> RemoveRun(seq5, 2, 0, "End")

    .. seealso::

        :func:`AddAbsorbingRun`,
        :func:`~openalea.stat_tool.cluster.Cluster`,
        :func:`Cumulate`,
        :func:`Difference`,
        :func:`IndexParameterExtract`,
        :func:`LengthSelect`,
        :func:`~openalea.stat_tool.data_transform.Merge`,
        :func:`~openalea.stat_tool.data_transform.MergeVariable`,
        :func:`MovingAverage`,
        :func:`RecurrenceTimeSequences`,
        :func:`Reverse`,
        :func:`SegmentationExtract`,
        :func:`~openalea.stat_tool.data_transform.SelectIndividual`,
        :func:`~openalea.stat_tool.data_transform.SelectVariable`,
        :func:`~openalea.stat_tool.data_transform.Shift`,
        :func:`~openalea.stat_tool.cluster.Transcode`,
        :func:`~openalea.stat_tool.data_transform.ValueSelect`,
        :func:`VariableScaling`.
    """
    error.CheckArgumentsLength(args, 2)
    MaxLength = kargs.get("MaxLength", I_DEFAULT) # max length will be used
    error.CheckType([MaxLength], [int])
    error.CheckType([obj, MaxLength], [[_Sequences, _MarkovianSequences,
               _VariableOrderMarkovData, _SemiMarkovData,
               _NonHomogeneousMarkovData], int])

    if len(args) == 3:
        variable = args[0]
        value = args[1]
        position = args[2]
        error.CheckType([variable, value, position], [int, int, str])

    elif len(args) == 2:
        variable = 1
        value = args[0]
        position = args[1]
        error.CheckType([variable, value, position], [int, int, str])

    if position in ['End', 'e']:
        position = 'e'
    elif position == ['b', 'Begin']:
        position = 'b'
    else:
        raise TypeError("position must be 'End' or 'e' or 'Begin' or 'b'")


    sequence = obj.remove_run(variable, value, position, MaxLength)

    if hasattr(sequence, 'markovian_sequences'):
        return sequence.markovian_sequences()
    else:
        return sequence


def ExtractVectors(obj, key, *args):
    """ExtractVectors

    Extraction of vectors from global characteristics of sequences (length or counting characteristics).

    :Usage:

    .. doctest::
        :options: +SKIP
    
        >>> ExtractVectors(seq, "Length")
        >>> ExtractVectors(seq1, "NbRun", value)
        >>> ExtractVectors(seq1, "NbOccurrence", value)
        >>> ExtractVectors(seqn, "NbRun", variable, value)
        >>> ExtractVectors(seqn, "NbOccurrence", variable, value)

    :Arguments:

    * seq (sequences, discrete_sequences, markov_data, semi-markov_data),
    * seq1 (discrete_sequences, markov_data, semi-markov_data): univariate sequences,
    * seqn (discrete_sequences, markov_data, semi-markov_data): multivariate sequences,
    * value (int): value,
    * variable (int): variable index.

    :Returned Object:

    An object of type vectors is returned.

    :Description:

    The type of global characteristic is given by a key word chosen among "Length", "NbRun" or "NbOccurrence". In the case of counting characteristics, the key word "NbRun" or "NbOccurrence" should be followed by a variable index in the case of multivariate sequences, and by the value of interest.

    .. seealso::

        :func:`~openalea.stat_tool.data_transform.ExtractHistogram`,
        :func:`~openalea.stat_tool.output.Plot`,
        :func:`~openalea.stat_tool.data_transform.MergeVariable`,
        :func:`~openalea.sequence_analysis.compare.Compare` (renewal process),
        :func:`~openalea.stat_tool.vectors.ContingencyTable`,,
        :func:`~openalea.stat_tool.regression.Regression`,
        :func:`~openalea.stat_tool.vectors.VarianceAnalysis`.
    """
    error.CheckArgumentsLength(args, 0, 2)
    error.CheckType([obj, key], [[_Sequences, _MarkovianSequences,
                             _VariableOrderMarkovData, _SemiMarkovData,
                             _NonHomogeneousMarkovData], str])


    error.CheckDictKeys(key, markovian_sequence_type)

    if key == "Length":
        variable = -1
        value = -1
    elif key in ["Cumul", "Mean"]:
        value = -1
        variable = obj.nb_variable
    elif key in ["NbRun", "FirstOccurrence", "NbOccurrence"]:
        if len(args)==2:
            variable = args[0]
            value = args[1]
        elif len(args)==1 and obj.nb_variable==1:
            value = args[0]
            variable = obj.nb_variable # should be equal to 1
            assert variable == 1

    error.CheckType([variable, value], [int, int])

    rkey = markovian_sequence_type[key]
    sequence = obj.extract_vectors(rkey, variable, value)

    return sequence





def SegmentationExtract(obj, variable, values , **kargs):
    """SegmentationExtract

    Extraction by segmentation of sub-sequences.

    :Usage:

    .. doctest::
        :options: +SKIP
    
        >>> SegmentationExtract(seqn, variable, value, Mode="Reject")
        >>> SegmentationExtract(seqn, variable, values, Mode="Reject")

    :Arguments:

    * seqn (sequences, discrete_sequences, markov_data, semi-markov_data):
      multivariate sequences,
    * variable (int): variable index,
    * value (int): value,
    * values (ARRAY(int)): values.

    :Optional Arguments:

    Mode (string): conservation or rejection of the selected sub-sequences:
    "Keep" (the default) or "Reject".

    :Returned Object:

    If all the variables are of type STATE, if variable is a valid index of
    variable and if either value or values[1], ..., values[n] are possible values,
    an object of type sequences or discrete_sequences is returned, otherwise no
    object is returned. The returned object is of type discrete_sequences if
    all the variables are of type STATE, if the possible values for each variable
    are consecutive from 0 and if the number of possible values for each
    variable is <= 15.

    :Description:

    Sub-sequences corresponding to run of value or values[1], ..., values[n] are
    extracted (or to all the possible values except value or
    values[1], ..., values[n]) are extracted.

    .. seealso::

        :func:`AddAbsorbingRun`,
        :func:`~openalea.stat_tool.cluster.Cluster`,
        :func:`Cumulate`,
        :func:`Difference`,
        :func:`IndexParameterExtract`,
        :func:`LengthSelect`,
        :func:`~openalea.stat_tool.data_transform.Merge`,
        :func:`~openalea.stat_tool.data_transform.MergeVariable`,
        :func:`MovingAverage`,
        :func:`RecurrenceTimeSequences`,
        :func:`RemoveRun`,
        :func:`Reverse`,
        :func:`~openalea.stat_tool.data_transform.SelectIndividual`,
        :func:`~openalea.stat_tool.data_transform.SelectVariable`,
        :func:`~openalea.stat_tool.data_transform.Shift`,
        :func:`~openalea.stat_tool.cluster.Transcode`,
        :func:`~openalea.stat_tool.data_transform.ValueSelect`,
        :func:`VariableScaling`.
    """

    mode = error.ParseKargs(kargs, "Mode", "Keep", keep_type)
    mode = bool(mode == "Keep" or mode == "keep")

    error.CheckType([obj, variable, values], [[_Sequences, _MarkovianSequences,
                            _VariableOrderMarkovData, _SemiMarkovData,
                            _NonHomogeneousMarkovData], int, [int, list]])

    if type(values) == int:
        sequence = obj.segmentation_extract(variable, [values], mode)
    elif type(values) == list:
        sequence = obj.segmentation_extract(variable, values, mode)
    else:
        raise TypeError("values must be an int or a list")

    if hasattr(sequence, 'markovian_sequences'):
        return sequence.markovian_sequences()
    else:
        return sequence


def LengthSelect(obj, minLength, *args, **kargs):
    """LengthSelect

    Selection of sequences according to a length criterion.

    :Usage:

    .. doctest::
        :options: +SKIP
    
        >>> LengthSelect(seq, length, Mode="Reject")
        >>> LengthSelect(seq, min_length, max_length, Mode="Reject")

    :Arguments:

    * seq (sequences, discrete_sequences, markov_data, semi-markov_data),
    * length (int): length,
    * min_length (int): minimum length,
    * max_length (int): maximum length.

    :Optional Arguments:

    Mode (string): conservation or rejection of the selected sequences: "Keep"
    (the default) or "Reject".

    :Returned Object:

    If length > 0 or if 0 < min_length < max_length and if the range of lengths
    defined either by length or by min_length and max_length enables to select
    sequences, an object of type sequences or discrete_sequences is returned,
    otherwise no object is returned. The returned object is of type
    discrete_sequences if all the variables are of type STATE, if the
    possible values for each variable are consecutive from 0 and if the
    number of possible values for each variable is < 15.

    .. seealso::

        :func:`AddAbsorbingRun`,
        :func:`~openalea.stat_tool.cluster.Cluster`,
        :func:`Cumulate`,
        :func:`Difference`,
        :func:`IndexParameterExtract`,
        :func:`~openalea.stat_tool.data_transform.Merge`,
        :func:`~openalea.stat_tool.data_transform.MergeVariable`,
        :func:`MovingAverage`,
        :func:`RecurrenceTimeSequences`,
        :func:`RemoveRun`,
        :func:`Reverse`,
        :func:`SegmentationExtract`,
        :func:`~openalea.stat_tool.data_transform.SelectIndividual`,
        :func:`~openalea.stat_tool.data_transform.SelectVariable`,
        :func:`~openalea.stat_tool.data_transform.Shift`,
        :func:`~openalea.stat_tool.cluster.Transcode`,
        :func:`~openalea.stat_tool.data_transform.ValueSelect`,
        :func:`VariableScaling`.
    """
    mode = error.ParseKargs(kargs, "Mode", "Keep", keep_type)
    mode = bool(mode == "Keep" or mode == "keep")

    maxLength = None
    if len(args) == 0:
        maxLength = minLength
    elif len(args)==1:
        maxLength = args[0]

    error.CheckType([obj], [[_Sequences, _MarkovianSequences,
                            _VariableOrderMarkovData, _SemiMarkovData,
                            _NonHomogeneousMarkovData]])

    error.CheckType([minLength, maxLength], [int, int])

    sequence =  obj.length_select(minLength, maxLength, mode)

    if hasattr(sequence, 'markovian_sequences'):
        return sequence.markovian_sequences()
    else:
        return sequence


def RecurrenceTimeSequences(obj, *args):
    """RecurrenceTimeSequences

    Computation of recurrence time sequences for a given value (and for a given variable).

    :Usage:

    .. doctest::
        :options: +SKIP
    
        >>> RecurrenceTimeSequences(seq1, value)
        >>> RecurrenceTimeSequences(seqn, variable, value)

    :Arguments:

    * seq1 (sequences, discrete_sequences, markov_data, semi-markov_data): univariate sequences,
    * seqn (sequences, discrete_sequences, markov_data, semi-markov_data): multivariate sequences,
    * variable (int): variable index,
    * value (int): value.

    :Returned Object:

    If the selected variable is of type STATE and if value is a possible value, an object of type sequences or discrete_sequences is returned, otherwise no object is returned. The returned object is of type discrete_sequences if all the variables are of type STATE, if the possible values for each variable are consecutive from 0 and if the number of possible values for each variable is <= 15.

    .. seealso::

        :func:`AddAbsorbingRun`,
        :func:`~openalea.stat_tool.cluster.Cluster`,
        :func:`Cumulate`,
        :func:`Difference`,
        :func:`IndexParameterExtract`,
        :func:`LengthSelect`,
        :func:`~openalea.stat_tool.data_transform.Merge`,
        :func:`~openalea.stat_tool.data_transform.MergeVariable`,
        :func:`MovingAverage`,
        :func:`RemoveRun`,
        :func:`Reverse`,
        :func:`SegmentationExtract`,
        :func:`~openalea.stat_tool.data_transform.SelectIndividual`,
        :func:`~openalea.stat_tool.data_transform.SelectVariable`,
        :func:`~openalea.stat_tool.data_transform.Shift`,
        :func:`~openalea.stat_tool.cluster.Transcode`,
        :func:`~openalea.stat_tool.data_transform.ValueSelect`,
        :func:`VariableScaling`.
    """
    error.CheckType([obj], [[_Sequences, _MarkovianSequences,
                            _VariableOrderMarkovData, _SemiMarkovData,
                            _NonHomogeneousMarkovData]])
    if len(args)==1:
        variable = 1
        value = args[0]
    elif len(args)==2:
        variable = args[0]
        value = args[1]

    error.CheckType([variable, value], [int, int])

    sequence = obj.recurrence_time_sequences(variable, value)

    if hasattr(sequence, 'markovian_sequences'):
        return sequence.markovian_sequences()
    else:
        return sequence

def WordCount(obj, *args, **kargs):
    """WordCount

    .. note::

        * args = variable, value,word_length
        * kargs = begin state, end_state, min_frequency

    .. todo:: the documentation...
    """
    error.CheckArgumentsLength(args, 1, 2)
    error.CheckType([obj], [[_MarkovianSequences, _VariableOrderMarkovData,
                            _SemiMarkovData, _NonHomogeneousMarkovData]])


    if len(args) == 1:
        variable = obj.nb_variable
        assert variable == 1
        word_length = args[0]
    elif len(args) == 2:
        variable = args[0]
        assert variable > 0
        word_length = args[1]
    error.CheckType([variable, word_length], [int, int])

    begin_state = error.ParseKargs(kargs, "BeginState", I_DEFAULT)
    end_state = error.ParseKargs(kargs, "EndState", I_DEFAULT)
    min_frequency = error.ParseKargs(kargs, "MinFrequency", 1)

    error.CheckType([begin_state, end_state, min_frequency], [int, int, int])

    return obj.word_count(variable, word_length, begin_state,
                          end_state, min_frequency)


def AddAbsorbingRun(obj, **kargs):
    """AddAbsorbingRun

    Addition of a run of absorbing vectors at the end of sequences.

    :Usage:

    .. doctest::
        :options: +SKIP
    
        >>> seq = Sequences([1,2,3,5])
        >>> AddAbsorbingRun(seq, RunLength=20)
        >>> AddAbsorbingRun(seq, SequenceLength=22, RunLength=20)

    :Arguments:

    * seq (distance_sequences, markov_data, semi-markov_data)

    :Optional Arguments:

    * SequenceLength (int): length of the sequences. A default value is
      computed from the maximum sequence length. must be less than
      max_length + 20.
    * RunLength (int): length of the runs. A default value is computed from
      the maximum sequence length. must be less than 20.

    :Returns:

    An object of type discrete_sequences is returned.

    .. seealso::

        :func:`~openalea.stat_tool.cluster.Cluster`,
        :func:`Cumulate`,
        :func:`IndexParameterExtract`,
        :func:`LengthSelect`,
        :func:`~openalea.stat_tool.data_transform.Merge`,
        :func:`~openalea.stat_tool.data_transform.MergeVariable`,
        :func:`MovingAverage`,
        :func:`RecurrenceTimeSequences`,
        :func:`~openalea.sequence_analysis.data_transform.Reverse`,
        :func:`~openalea.sequence_analysis.data_transform.SegmentationExtract`,
        :func:`~openalea.stat_tool.data_transform.SelectIndividual`,
        :func:`~openalea.stat_tool.data_transform.SelectVariable`,
        :func:`~openalea.stat_tool.data_transform.Shift`,
        :func:`~openalea.stat_tool.cluster.Transcode`,
        :func:`~openalea.stat_tool.data_transform.ValueSelect`,
        :func:`~openalea.sequence_analysis.data_transform.VariableScaling`.

    """
    error.CheckType([obj], [[_MarkovianSequences, _VariableOrderMarkovData,
                            _SemiMarkovData, _NonHomogeneousMarkovData,
                            _Tops]])

    SequenceLength = error.ParseKargs(kargs, "SequenceLength", I_DEFAULT)
    RunLength = error.ParseKargs(kargs, "RunLength", I_DEFAULT)

    #todo: replace this by sequence_analysis constant
    MAX_ABSORBING_RUN_LENGTH = 20

    try:
        return obj.add_absorbing_run(SequenceLength, RunLength)
    except:

        max_max_length = obj.max_length + MAX_ABSORBING_RUN_LENGTH
        raise ValueError("abnormal termination. "
                         "Check the arguments."
                         "RunLength must be positive and  less than %s"
                         % MAX_ABSORBING_RUN_LENGTH,
                         "SequenceLength must be positive, larger than "
                         "max_length (%s) and less than max_length +"
                         " MAX_RUN_LENGTH (%s)"
                         % (obj.max_length, max_max_length))

def Reverse(obj):
    """Reverse

    Reversing of sequences or 'tops'.

    :Usage:

    .. doctest::
        :options: +SKIP
        
        >>> Reverse(seq)
        >>> Reverse(discrete_seq)
        >>> Reverse(top)

    :Arguments:

    * seq (sequences),
    * discrete_seq (discrete_sequences, markov_data, semi-markov_data),
    * top (tops).

    :Returned Object:

    If the argument is of type sequences, an object of type sequences is returned. If the argument is of type discrete_sequences, markov_data, semi-markov_data, an object of type discrete_sequences is returned. If the argument is of type tops, an object of type tops is returned.

    .. seealso::

        :func:`AddAbsorbingRun`,
        :func:`~openalea.stat_tool.cluster.Cluster`,
        :func:`Cumulate`,
        :func:`Difference`,
        :func:`IndexParameterExtract`,
        :func:`LengthSelect`,
        :func:`~openalea.stat_tool.data_transform.Merge`,
        :func:`~openalea.stat_tool.data_transform.MergeVariable`,
        :func:`MovingAverage`,
        :func:`RecurrenceTimeSequences`,
        :func:`RemoveRun`,
        :func:`~openalea.sequence_analysis.tops.RemoveApicalInternodes`,
        :func:`SegmentationExtract`,
        :func:`~openalea.stat_tool.data_transform.SelectIndividual`,
        :func:`~openalea.stat_tool.data_transform.SelectVariable`,
        :func:`~openalea.stat_tool.data_transform.Shift`,
        :func:`~openalea.stat_tool.cluster.Transcode`,
        :func:`~openalea.stat_tool.data_transform.ValueSelect`,
        :func:`VariableScaling`.
    """

    error.CheckType([obj], [[_MarkovianSequences, _VariableOrderMarkovData,
                            _SemiMarkovData, _NonHomogeneousMarkovData,
                            _Tops]])

    ret = obj.reverse()

    if hasattr(ret, 'markovian_sequences'):
        return ret.markovian_sequences()
    else:
        return ret


def Thresholding(obj, MinProbability=MIN_PROBABILITY):
    """

    .. todo:: documentation
    """
    error.CheckType([MinProbability], [[int, float]])
    error.CheckType([obj], [[_VariableOrderMarkov, _HiddenSemiMarkov,
                            _HiddenVariableOrderMarkov, _SemiMarkov]])

    return obj.thresholding(MinProbability)


def Cumulate(obj, Variable=I_DEFAULT):
    """Cumulate

    Sum of successive values along sequences.

    :Usage:

    .. doctest::
        :options: +SKIP
    
        >>> Cumulate(seq, Variable=1)

    :Arguments:

    * seq (sequences, discrete_sequences, markov_data, semi-markov_data).

    :Optional Arguments:

    * Variable (int): variable index.

    :Returned Object:

    The returned object is of type sequences or discrete_sequences. The returned object is of type discrete_sequences if all the variables are of type STATE, if the possible values for each variable are consecutive from 0 and if the number of possible values for each variable is < 15.

    :Background:

    Cumulate is the inverse function of Difference with the optional argument FirstValue set at True.

    .. seealso::

        :func:`AddAbsorbingRun`,
        :func:`~openalea.stat_tool.cluster.Cluster`,
        :func:`Difference`,
        :func:`IndexParameterSelect`,
        :func:`LengthSelect`,
        :func:`~openalea.stat_tool.data_transform.Merge`,
        :func:`~openalea.stat_tool.data_transform.MergeVariable`,
        :func:`MovingAverage`,
        :func:`RecurrenceTimeSequences`,
        :func:`RemoveSeries`,
        :func:`Reverse`,
        :func:`SegmentationExtract`,
        :func:`~openalea.stat_tool.data_transform.SelectIndividual`,
        :func:`~openalea.stat_tool.data_transform.SelectVariable`,
        :func:`~openalea.stat_tool.data_transform.Shift`,
        :func:`~openalea.stat_tool.cluster.Transcode`,
        :func:`~openalea.stat_tool.data_transform.ValueSelect`,
        :func:`VariableScaling`.
    """

    error.CheckType([obj, Variable], [[_Sequences, _MarkovianSequences,
                                       _VariableOrderMarkovData,
                                       _SemiMarkovData,
                                       _NonHomogeneousMarkovData], int])

    ret = obj.cumulate(Variable)


    if hasattr(ret, 'markovian_sequences'):
        return ret.markovian_sequences()
    else:
        return ret



def Difference(obj, Variable=I_DEFAULT, FirstElement=False):
    """Difference

    First-order differencing of sequences.

    :Usage:
    
    .. doctest::
        :options: +SKIP
        
        Difference(seq, Variable=1, FirstElement=True)

    :Arguments:

    * seq (sequences, discrete_sequences, markov_data, semi-markov_data).

    :Optional Arguments:

    * Variable (int): variable index,
    * FirstElement (bool): first element of sequences kept or not (default value: False).

    :Returned Object:

    The returned object is of type sequences or discrete_sequences. The returned object is of type discrete_sequences if all the variables are of type STATE, if the possible values for each variable are consecutive from 0 and if the number of possible values for each variable is < 15.

    :Background:

    If the first element of sequences are kept (FirstValue=True), Difference is the inverse function of Cumulate.

    .. seealso::

        :func:`AddAbsorbingRun`,
        :func:`~openalea.stat_tool.cluster.Cluster`,
        :func:`Cumulate`,
        :func:`IndexParameterExtract`,
        :func:`LengthSelect`,
        :func:`~openalea.stat_tool.data_transform.Merge`,
        :func:`~openalea.stat_tool.data_transform.MergeVariable`,
        :func:`MovingAverage`,
        :func:`RecurrenceTimeSequences`,
        :func:`RemoveSeries`,
        :func:`Reverse`,
        :func:`SegmentationExtract`,
        :func:`~openalea.stat_tool.data_transform.SelectIndividual`,
        :func:`~openalea.stat_tool.data_transform.SelectVariable`,
        :func:`~openalea.stat_tool.data_transform.Shift`,
        :func:`~openalea.stat_tool.cluster.Transcode`,
        :func:`~openalea.stat_tool.data_transform.ValueSelect`,
        :func:`VariableScaling`.
    """
    error.CheckType([obj, Variable, FirstElement],
                    [[_Sequences, _MarkovianSequences,
                      _VariableOrderMarkovData, _SemiMarkovData,
                      _NonHomogeneousMarkovData], int, bool])

    ret = obj.difference(Variable, FirstElement)

    #todo: markovian_sequences does nto seem to work
    #if hasattr(ret, 'markovian_sequences'):
    #    return ret.markovian_sequences()
    #else:
    #    return ret
    return ret


def IndexParameterExtract(obj, minIndex, MaxIndex=I_DEFAULT):
    """IndexExtract

    Extraction of sub-sequences corresponding to a range of index parameters.

    :Usage:
    
    .. doctest::
        :options: +SKIP

        >>> IndexParameterExtract(seq, min_index, MaxIndex=40)
        >>> IndexParameterExtract(seq, min_index, MaxIndex)
        
    :Arguments:

    * seq (sequences, discrete_sequences, markov_data, semi-markov_data),
    * min_index (int): minimum index parameter.

    :Optional Arguments:

    * MaxIndex (int): maximum index parameter (default behaviour: the end of sequences is kept).

    :Returned Object:

    If 0 < min_index < (maximum index parameter if the optional argument MaxIndex is set) < (maximum length of sequences), the returned object is of type sequences or discrete_sequences, otherwise no object is returned. The returned object is of type discrete_sequences if all the variables are of type STATE, if the possible values for each variable are consecutive from 0 and if the number of possible values for each variable is < 15.

    .. seealso::

        :func:`AddAbsorbingRun`,
        :func:`~openalea.stat_tool.cluster.Cluster`,
        :func:`Cumulate`,
        :func:`Difference`,
        :func:`LengthSelect`,
        :func:`~openalea.stat_tool.data_transform.Merge`,
        :func:`~openalea.stat_tool.data_transform.MergeVariable`,
        :func:`MovingAverage`,
        :func:`RecurrenceTimeSequences`,
        :func:`RemoveRun`,
        :func:`Reverse`,
        :func:`SegmentationExtract`,
        :func:`~openalea.stat_tool.data_transform.SelectIndividual`,
        :func:`~openalea.stat_tool.data_transform.SelectVariable`,
        :func:`~openalea.stat_tool.data_transform.Shift`,
        :func:`~openalea.stat_tool.cluster.Transcode`,
        :func:`~openalea.stat_tool.data_transform.ValueSelect`,
        :func:`VariableScaling`.
    """
    error.CheckType([obj, minIndex], [[_Sequences, _MarkovianSequences,
                             _VariableOrderMarkovData, _SemiMarkovData,
                             _NonHomogeneousMarkovData], int])


    ret = obj.index_parameter_extract(minIndex, MaxIndex)

    try:
        if hasattr(ret, "markovian_sequences"):
            return ret.markovian_sequences()
    except:
        return ret

def IndexParameterSelect(obj, minIndex, *args, **kargs):
    """IndexExtract

    Extraction of sub-sequences corresponding to a range of index parameters.

    :Usage:
    
    .. doctest::
        :options: +SKIP

        >>> IndexParameterSelect(seq, min_index)
        >>> IndexParameterSelect(seq, min_index, max_index)
        >>> IndexParameterSelect(seq, min_index, max_index, Mode="Keep")

    :Arguments:

    * seq (sequences, discrete_sequences, markov_data, semi-markov_data),
    * min_index (int): minimum index parameter.

    :Optional Arguments:

    * max_index (int): maximum index parameter (default behaviour: the end of 
      sequences is kept).
    * Mode: Keep or Reject
    
    :Returned Object:

    If 0 < min_index < (maximum index parameter if the optional argument 
    MaxIndex is set) < (maximum length of sequences), the returned object is of
    type sequences or discrete_sequences, otherwise no object is returned. The
    returned object is of type discrete_sequences if all the variables are of 
    type STATE, if the possible values for each variable are consecutive from 0 
    and if the number of possible values for each variable is < 15.

    .. seealso::

        :func:`AddAbsorbingRun`,
        :func:`~openalea.stat_tool.cluster.Cluster`,
        :func:`Cumulate`,
        :func:`Difference`,
        :func:`LengthSelect`,
        :func:`~openalea.stat_tool.data_transform.Merge`,
        :func:`~openalea.stat_tool.data_transform.MergeVariable`,
        :func:`MovingAverage`,
        :func:`RecurrenceTimeSequences`,
        :func:`RemoveRun`,
        :func:`Reverse`,
        :func:`SegmentationExtract`,
        :func:`~openalea.stat_tool.data_transform.SelectIndividual`,
        :func:`~openalea.stat_tool.data_transform.SelectVariable`,
        :func:`~openalea.stat_tool.data_transform.Shift`,
        :func:`~openalea.stat_tool.cluster.Transcode`,
        :func:`~openalea.stat_tool.data_transform.ValueSelect`,
        :func:`VariableScaling`.
    """
    error.CheckArgumentsLength(args, 0, 1)
    error.CheckType([obj], [[_Sequences, _MarkovianSequences,
                             _VariableOrderMarkovData, _SemiMarkovData,
                             _NonHomogeneousMarkovData]])

    mode = error.ParseKargs(kargs, "Mode", "Keep", keep_type)
    mode = bool(mode == "Keep" or mode == "keep")


    if len(args) == 0:
        maxIndex = minIndex
    elif len(args)==1:
        maxIndex = args[0]

    error.CheckType([minIndex, maxIndex], [int, int])

    ret = obj.index_parameter_select(minIndex, maxIndex, mode)

    if hasattr(ret, "markovian_sequences"):
        return ret.markovian_sequences()
    else:
        return ret


def ComputeStateSequences(obj, data, **kargs):
    """ComputeStateSequences

    Computation of the optimal state sequences corresponding to the observed sequences using a hidden Markov chain or a hidden semi-Markov chain.

    :Usage:
    
    .. doctest::
        :options: +SKIP

        >>> ComputeStateSequences(seq, hmc, Algorithm="ForwardBackward", Characteristics=True)
        >>> ComputeStateSequences(seq, hsmc, Algorithm="ForwardBackward", Characteristics=True)

    :Arguments:

    * hc (hidden_markov),
    * hsmc (hidden_semi-markov),
    * seq (discrete_sequences, markov_data, semi-markov_data).

    :Optional Arguments:

    * Algorithm (string): type of algorithm: "Viterbi" (the default) or "ForwardBackward".
    * Characteristics (bool): computation of the characteristic distributions of the model taking into account the lengths of the segmented sequences (default value: False).

    :Returned Object:

    If the second mandatory argument is of type hidden_markov, an object of type markov_data is returned. If the second mandatory argument is of type hidden_semi-markov, an object of type semi-markov_data is returned. The returned object contains both the sequences (including the optimal state sequences as a supplementary variable) and the model.

    :Background:

    This restoration of the state sequence is either performed by a dynamic
    programming algorithm referred to as Viterbi algorithm which maximizes
    the joint probability of the state and output sequence  (global criterion)
    or by a forward-backward algorithm which chooses state j at time t to maximize  (local criterion).

    .. seealso::

        :func:`~openalea.stat_tool.data_transform.ExtractHistogram`.

    """
    error.CheckType([obj], [[_MarkovianSequences,
                             _VariableOrderMarkovData, _SemiMarkovData]])

    Characteristics = error.ParseKargs(kargs, "Characteristics",
                                       True, bool_type)

    error.CheckType([data], [[_HiddenVariableOrderMarkov,
                             _HiddenSemiMarkov]])
    ret = data.state_sequence_computation(obj, Characteristics)

    return ret

def MovingAverage(obj, itype, *args, **kargs):
    """MovingAverage

    Extraction of trends or residuals using a symmetric smoothing filter.

    :Usage:
    
    .. doctest::
        :options: +SKIP

        >>> MovingAverage(seq, filter, Variable=1, BeginEnd=True, Output="Residual", FileName="filtered_sequences")
        >>> MovingAverage(seq, frequencies, Variable=1, BeginEnd=True,  Output="Residual", FileName="filtered_sequences")
        >>> MovingAverage(seq, dist, Variable=1, BeginEnd=True,  Output="Residual", FileName="filtered_sequences")

    :Arguments:

    * seq (sequences, discrete_sequences, markov_data, semi-markov_data)
    * filter (array(real)): filter values on the half width i.e. from one extremity to the central value (with the constraintfilter@i + filter@m = 1),
    * frequencies (array(int)): frequencies defining the filter,
    * dist (distribution, mixture, convolution, compound): symmetric distribution, whose size of the support is even, defining the filter (for instance Distribution("BINOMIAL", 0, 4, 0.5)),

    :Optional Arguments:

    * Variable (int): variable index,
    * BeginEnd (bool): begin and end of sequences filtered or suppresses (default value: False),
    * Output (string): "Trend" (the default), "Residual" or "DivisionResidual",
    * FileName (string): name of file of non-rounded filtered sequences.

    :Background:

    Consider a symmetric smoothing filter of half width q applied to a sequence of length t. Whenever a symmetric smoothing filter is chosen, there is likely to be an end-effect problem for t<q or t>t - q - 1. We chose to apply the following solution to the first and the last q terms of the sequences: we define Xt:=X0 for t<0 and Xt:=Xt-1 for t>t-1.

    :Returned Object:

    An object of type sequences or discrete_sequences is returned. An object of type discrete_sequences is returned if all the variables are of type STATE, if the possible values for each variable are consecutive from 0 and if the number of possible values for each variable is <= 15.

    .. seealso::

        :func:`AddAbsorbingRun`,
        :func:`~openalea.stat_tool.cluster.Cluster`,
        :func:`Cumulate`,
        :func:`Difference`,
        :func:`IndexParameterExtract`,
        :func:`IndexParameterSelect`,
        :func:`~openalea.stat_tool.data_transform.Merge`,
        :func:`~openalea.stat_tool.data_transform.MergeVariable`,
        :func:`RecurrenceTimeSequences`,
        :func:`RemoveRun`,
        :func:`Reverse`,
        :func:`SegmentationExtract`,
        :func:`~openalea.stat_tool.data_transform.SelectIndividual`,
        :func:`~openalea.stat_tool.data_transform.SelectVariable`,
        :func:`~openalea.stat_tool.data_transform.Shift`,
        :func:`~openalea.stat_tool.cluster.Transcode`,
        :func:`~openalea.stat_tool.data_transform.ValueSelect`,
        :func:`VariableScaling`.


    .. todo:: seq79 = MovingAverage(seq70, [1, 1, 1], BeginEnd=True, Output="Residual")
        returns wrong results. comapre with aml

    """

    error.CheckArgumentsLength(args, 0)
    error.CheckType([obj], [[_Sequences, _MarkovianSequences,
                             _VariableOrderMarkovData, _SemiMarkovData,
                             _NonHomogeneousMarkovData]])
    error.CheckType([itype], [[list, _DiscreteParametricModel, _DiscreteMixture,
                               _Convolution, _Compound]])
    Variable = error.ParseKargs(kargs, "Variable", I_DEFAULT)
    Output = error.ParseKargs(kargs, "Output", "Trend", func_map)
    BeginEnd = error.ParseKargs(kargs, "BeginEnd", False,
                                {"False":False, "True":True,
                                 False:False, True:True})

    if isinstance(itype, list):
        #todo build a 2*N+1 filter here or inside moving_average function ?
        return obj.moving_average(itype, Variable, BeginEnd, Output)
    else: #distribution
        return obj.moving_average_from_distribution(itype, Variable,
                                                    BeginEnd, Output)


def ComputeSelfTransition(obj, Order=1):
    """ComputeSelfTransition

    Computation of the self-transition probabilities as a function of the 
    index parameter from discrete sequences.

    :Usage:
    
    .. doctest::
        :options: +SKIP

        >>> ComputeSelfTransition(seq, Order=2)

    :Arguments:

    * seq (discrete_sequences, markov_data, semi-markov_data).

    :Optional Arguments:

    * Order (int): Markov chain order (default value: 1).

    :Returned Object:

    No object returned.
    """
    #todo: Order is not used
    error.CheckType([obj], [[_MarkovianSequences,
                             _VariableOrderMarkovData, _SemiMarkovData,
                             _NonHomogeneousMarkovData]])

    obj.self_transition_computation()


def TransitionCount(obj, max_order, **kargs):

    """.. todo:: documentation"""
    error.CheckType([obj], [[_MarkovianSequences,
                             _VariableOrderMarkovData, _SemiMarkovData,
                             _NonHomogeneousMarkovData]])
    begin = error.ParseKargs(kargs, "begin", False, bool_type)
    estimator = error.ParseKargs(kargs, "Estimator", "Maximum_Likelihood",
                                 estimator_map)
    error.CheckType([max_order, begin], [int, bool])

    #todo: check that format_type =ascii/spreadheet is expected. Just the name
    #may be sufficient
    filename = error.ParseKargs(kargs, "Filename", None, format_type)

    obj.transition_count(max_order, begin, estimator, filename)


def Cross(obj):
    """.. todo:: documentation"""
    error.CheckType([obj], [[_Sequences, _MarkovianSequences,
                             _VariableOrderMarkovData, _SemiMarkovData,
                             _NonHomogeneousMarkovData]])
    sequence =  obj.cross()
    if hasattr(sequence, 'markovian_sequences'):
        return sequence.markovian_sequences()
    else:
        return sequence


def PointwiseAverage(obj, **kargs):
    """.. todo:: documentation"""


    error.CheckType([obj], [[_Sequences, _MarkovianSequences,
                             _VariableOrderMarkovData, _SemiMarkovData,
                             _NonHomogeneousMarkovData]])

    StandardDeviation = error.ParseKargs(kargs, "StandardDeviation", False)
    Output = error.ParseKargs(kargs, "Output", "Sequence", output_map)
    Filename = kargs.get("Filename", None)
    Format = error.ParseKargs(kargs, "Format", 'ASCII', format_type)

    return obj.pointwise_average(StandardDeviation, Output, Filename, Format)


def ConsecutiveValues(obj, *args, **kargs):
    """ConsecutiveValues"""
    error.CheckType([obj], [[_MarkovianSequences,
                             _VariableOrderMarkovData, _SemiMarkovData,
                             _NonHomogeneousMarkovData]])
    AddVariable = kargs.get("AddVariable", False)

    if len(args) == 1:
        variable = args[0]
    else:
        variable = 1

    variable = _check_nb_variable(obj, variable)
    error.CheckType([variable, AddVariable], [int, bool])
    sequence = obj.consecutive_values(variable, AddVariable)

    return sequence



def Round(obj,  **kargs):
    """Round"""

    Variable = kargs.get("Variable", I_DEFAULT)
    Mode = kargs.get("Mode", "Round")
    Mode = mode_type[Mode]

    ret =  obj.round(Variable, Mode)
    if hasattr(ret, 'markovian_sequences'):
        return ret.markovian_sequences()
    else:
        return ret


def TimeScaling(obj, scaling_factor=0):
    """TimeScaling

    Change of the time unit of data of type {time interval between two observation dates, number of events occurring between these two observation dates}.

    :Usage:
    .. doctest::
        :options: +SKIP

        >>> TimeScaling(timev, scaling_factor)

    :Arguments:

    * timev (time_events, renewal_data),
    * scaling_factor (int): scaling factor.

    :Returned Object:

    If scaling_factor > 0, an object of type time_events is returned, otherwise no object is returned.

    .. seealso::

        :func:`~openalea.stat_tool.data_transform.Merge`,
        :func:`~openalea.sequence_analysis.time_events.NbEventSelect`,
        :func:`TimeSelect`.
    """
    error.CheckType([obj, scaling_factor],
                    [[ _TimeEvents, _RenewalData], int])
    return obj.time_scaling(scaling_factor)

def TimeSelect(obj, *args):
    """TimeSelect

    Selection of data item of type {time interval between two observation dates, number of events occurring between these two observation dates} according to a length of the observation period criterion.

    :Usage:
    
    .. doctest::
        :options: +SKIP

        >>> TimeSelect(timev, time)
        >>> TimeSelect(timev, min_time, max_time)

    :Arguments:

    * timev (time_events, renewal_data),
    * time (int): time interval between two observation dates,
    * min_time (int): minimum time interval between two observation dates,
    * max_time (int): maximum time interval between two observation dates.

    :Returned Object:

    If either time > 0 or if 0 < min_time < max_time and if the range of lengths of the observation period defined either by time or by min_time and max_time enables to select data items of type {time interval between two observation dates, number of events}, an object of type time_events is returned, otherwise no object is returned.

    .. seealso::

        :func:`~openalea.stat_tool.data_transform.Merge`,
        :func:`~openalea.sequence_analysis.time_events.NbEventSelect`,
        :func:`TimeScaling`.
    """
    error.CheckType([obj], [[_TimeEvents, _RenewalData]])
    error.CheckArgumentsLength(args, 1, 2)

    if len(args) == 2:
        min_time = args[0]
        max_time = args[1]
    elif len(args) == 1:
        min_time = args[0]
        max_time = min_time

    error.CheckType([min_time, max_time], [int, int])

    return obj.time_select(min_time, max_time)


def Segmentation(obj, *args, **kargs):
    """Segmentation

    :Usage:

    .. doctest::
        :options: +SKIP

        >>> Segmentation(seq, 1,2, VectorDistance("N"), Output="Segment")
    
    """
    error.CheckArgumentsLength(args, 2)
    error.CheckType([obj], [[_Sequences, _MarkovianSequences,
                             _VariableOrderMarkovData, _SemiMarkovData]])

    nb_variable = obj.nb_variable

    error.CheckType([args[0]], [[int, list]])

    # Segmentation from list and strings
    # >>> Segmentation(seq80, [6, 5, 5, 6, 4, 4], "Mean", Output="Residual")
    # >>> Segmentation(seq80, [5, 5, 6, 5, 5, 4], "Gaussian")
    if isinstance(args[0], list):
        nb_segments = obj.nb_sequence
        nb_segments = args[0]
        #segments = args[0]
        Model = __get_model__(args[1:], nb_variable)
        Output = error.ParseKargs(kargs, "Output", "Sequence", sub_func_map)
        identifier = I_DEFAULT

        seq = obj.segmentation_array(nb_segments, Model, identifier , Output)


    # Segmentation(seq, 1,2, VectorDistance("N"), Output="Segment")
    # elif type(args[2]) == _VectorDistance:
    #     error.CheckType([args[0], args[1]], [int, int])
    #     Output = error.ParseKargs(kargs, "Output", "Sequence", output_type)
    #     seq = obj.segmentation_vector_distance(args[0], args[1], args[2],
    #                                            Output)

    # Segmenation when args[0] is int
    # >>> Segmentation(seq80, 2, 6, "Gaussian", NbChangePoint="Fixed",
    # ... Output="Residual")
    elif isinstance(args[0], int):
        error.CheckType([args[1]], [[int, list]])
        Model = __get_model__(args[2:], nb_variable)
        # >>> Segmentation(seq80, 2, 10, "Gaussian")
        if type(args[1]) == int:
            #error.CheckType([args[2]], [[int]])
            nb_segment_estimation = error.ParseKargs(kargs, "NbSegment",
                                                     "Estimated",
                                                     nb_segment_map)

            Output = error.ParseKargs(kargs, "Output", "Sequence", sub_func_map)

            if kargs.get("NbSegment") and kargs.get("Output"):
                raise ValueError("NbSegment and Output option incompatible")

            if nb_segment_estimation is False:
                if args[1] == 1:
                    seq = obj.segmentation_int_int(args[0], args[1], Model,
                                                   Output)
                else:
                    nb_segment = args[0]
                    seq = obj.segmentation_array(nb_segment, Model,
                                                args[0], Output)
            else:
                seq = obj.segmentation_model(args[0], args[1], Model)

        else:
            # type(args[1]) == list
            # >>> Segmentation(seq70, 1, [1935, 1961, 1972, 1990],
            #    "Gaussian", "Gaussian"),

            #todo: to finish
            error.CheckType([args[0], args[1]], [int, list])

            Model = __get_model__(args[2:], nb_variable)
            Output = error.ParseKargs(kargs, "Output", "Sequence", sub_func_map)

            change_point = args[1]

            nb_segment = len(args[1]) + 1

            seq = obj.segmentation_change_point(args[0] , nb_segment ,
                            change_point , Model , Output)
    return seq


def SojournTimeSequences(obj, Variable = 1):
    """.. todo:: documentation"""
    error.CheckType([obj], [[_Sequences, _MarkovianSequences,
                             _VariableOrderMarkovData, _SemiMarkovData,
                             _NonHomogeneousMarkovData]])

    #todo:  check that variable is valid
    ret =  obj.sojourn_time_sequences(Variable)
    if hasattr(ret, 'markovian_sequences'):
        return ret.markovian_sequences()
    else:
        return ret

def VariableScaling(obj, *args):
    """VariableScaling

    Change of the unit of a variable.

    :Usage:
    .. doctest::
        :options: +SKIP

        >>> VariableScaling(seq1, scaling_factor)
        >>> VariableScaling(seqn, variable, scaling_factor)

    :Arguments:

    * seq1 (sequences, discrete_sequences, markov_data, semi-markov_data): univariate sequences,
    * seqn (sequences, discrete_sequences, markov_data, semi-markov_data): multivariate sequences.
    * variable (int): variable index,
    * scaling_factor (int): scaling factor.

    :Returned Object:

    If scaling_factor > 0, an object of type sequences is returned, otherwise no object is returned.

    :Background:

    The function VariableScaling is mainly useful as a preprocessing when one wants to study the correlation structure of residual sequences. This function enables by an appropriate scaling to control the rounding of the residual sequences and hence to obtain exact sample correlation functions.

    .. seealso::

        :func:`AddAbsorbingRun`,
        :func:`~openalea.stat_tool.cluster.Cluster`,
        :func:`Cumulate`,
        :func:`Difference`,
        :func:`IndexParameterExtract`,
        :func:`LengthSelect`,
        :func:`~openalea.stat_tool.data_transform.Merge`,
        :func:`~openalea.stat_tool.data_transform.MergeVariable`,
        :func:`MovingAverage`,
        :func:`RecurrenceTimeSequences`,
        :func:`RemoveRun`,
        :func:`Reverse`,
        :func:`SegmentationExtract`,
        :func:`~openalea.stat_tool.data_transform.SelectIndividual`,
        :func:`~openalea.stat_tool.data_transform.SelectVariable`,
        :func:`~openalea.stat_tool.data_transform.Shift`,
        :func:`~openalea.stat_tool.cluster.Transcode`,
        :func:`~openalea.stat_tool.data_transform.ValueSelect`.
    """
    error.CheckType([obj], [[_Vectors, _Sequences, _MarkovianSequences,
                             _VariableOrderMarkovData, _SemiMarkovData,
                             _NonHomogeneousMarkovData]])

    if len(args) == 2:
        variable = args[0]
        scaling = args[1]
    elif len(args) == 1:
        variable = 1
        scaling = args[0]

    _check_nb_variable(obj, variable)

    error.CheckType([variable, scaling], [int, int])
    return obj.scaling(variable, scaling)

def __get_model__(data, nb_variable):
    """ """

    Model = []
    for i in range(0, nb_variable):
        Model.append(model_type[data[i]])

        if (i == 0) and (model_type[data[i]] == MEAN_CHANGE) or \
        (model_type[data[i]] == MEAN_VARIANCE_CHANGE):

            for _j in range(1, nb_variable):
                Model.append(model_type[i])
            break #todo check that it break the main for i in range() loop

    return Model

def vec2list(vector):
    """
    transform a vector from stat_tool vector into a python list
    """
    output = []
    for i in range(0, vector.nb_vector):
        output.append(vector[i][0])
    return output


def Extract(obj, *args, **kargs):
    """
    .. todo:: check the case elif type(obj) in [_Sequences, _MarkovianSequences,
        _VariableOrderMarkovData, _SemiMarkovData, _NonHomogeneousMarkovData]
    """

    if isinstance(obj, _Tops):
        raise NotImplemented("Extract Tops case not yet implemented")
    elif isinstance(obj, _TopParameters):
        position = args[0]
        error.CheckType(position, [int])
        ret = obj.extract(position)

    elif type(obj) in [_TimeEvents, _RenewalData]:
        if len(args) >=1:
            key = args[0]
            if key == "NbEvent":
                error.CheckType([args[1]], [int])
                ret = obj.extract(renewal_nb_event_map[key], args[1])

            elif key == "ObservationTime":
                ret = obj.get_htime()

            elif key in renewal_nb_event_map.keys():

                ret = obj.extract(renewal_nb_event_map[key], I_DEFAULT)
        else:
            raise NotImplemented("Not yet implemented")

    # renewal/timeevents case
    elif type(obj) in  [_Renewal]:
        key = args[0]
        error.CheckType([key], [str])
        error.CheckDictKeys(key, renewal_nb_event_map)


        if key == "NbEvent":
            error.CheckArgumentsLength(args, 2, 2)
            time = args[1]
            ret =  obj.extract(renewal_nb_event_map[key], time)
        else:
            time = 0
            ret = obj.extract(renewal_nb_event_map[key], time)


    # todo: bug ? here, we need to specify the _sequence_analysis module
    # before _MarkovianSequences otherwise even though obj is
    # Markovian_Seqiences,
    # it is not found. Found with stat_tool/test/stat_tool_test.py
    elif type(obj) in [_Sequences,
                       _MarkovianSequences,
                       _VariableOrderMarkovData,
                       _SemiMarkovData,
                       _NonHomogeneousMarkovData]:
        error.CheckType([args[0]], [str])
        if args[0] == 'Value':
            if obj.nb_variable == 1:
                error.CheckArgumentsLength(args, 1, 1)
                #todo: the following call does not work
                #error.CheckArgumentsLength(args, 2, 2)
                variable = 1
            else:
                error.CheckArgumentsLength(args, 1, 2)
                #todo: the following call does not work
                #error.CheckArgumentsLength(args, 3, 3)
                error.CheckType([args[1]], [int])
                variable = args[1]

            ret = obj.extract_value(variable)

        elif args[0] == 'Length':
            #todo: the following call does not work
            error.CheckArgumentsLength(args, 1, 1)
            ret = obj.extract_length()
        else:
            ident = args[0]
            if ident not in markovian_sequence_type.keys():
                print 'Wrong type %s. Use one of %s ' % (ident, markovian_sequence_type.keys())
                raise ValueError('error in extract histogram to be done')
            if obj.nb_variable == 1:
                variable = 1
                value = args[1]
                error.CheckType([args[1]], [int])
                ret = obj.extract(markovian_sequence_type[ident], variable, value)
            else:
                variable = args[1]
                value = args[2]
                error.CheckType([args[1], args[2]], [[int], [int]])
                ret = obj.extract(markovian_sequence_type[ident], variable, value)


    elif type(obj) in [ _VariableOrderMarkov,
                       _HiddenVariableOrderMarkov,
                       _SemiMarkov,
                       _NonHomogeneousMarkov,
                       _HiddenSemiMarkov]:
        error.CheckType([args[0]], [str])
        key = seq_map[args[0]]
        # forward is just for semi_markov and hidden+semi_markov
        # add a CheckType. In addition, args[2] must be int ?
        if args[0] == "Forward":

            error.CheckArgumentsLength(args, 2, 3)
            error.CheckType([args[1]], [int])
            error.CheckType([args[0]], [[_SemiMarkov, _HiddenSemiMarkov]])
            HistogramType = kargs.get("HistogramType", "FinalRun")
            HistogramType = histogram_type[HistogramType]
            ret = obj.extract(seq_map[args[0]], args[1], HistogramType)

        elif len(args) == 2:
            # todo: ident cannot be Observation
            variable = 1
            value = args[1]
            error.CheckType([value], [int])
            ret = obj.extract(seq_map[args[0]], variable, value)

        elif len(args)==3:
            #todo: check validity of variable ?
            error.CheckType([args[1], args[2]], [int, int])
            variable = args[1]
            value = args[2]
            print seq_map[args[0]]
            ret = obj.extract(seq_map[args[0]], variable, value)
        else:
            raise KeyError("expect only 1 or 2 arguments after the type")

    else:
        print 'switch to openalea.stat_tool'
        # related to Top, Renewal, Markov that are in sequence
        try:
            from openalea.stat_tool.data_transform import Extract as newExtract
            ret = newExtract(obj, *args, **kargs)
        except:
            pass
    #todo check if ret not None
    if ret is None:
        raise Exception("return object cannot not be None")
    else:
        return ret


def BuildAuxiliaryVariable(obj):
    error.CheckType([obj], [[ _VariableOrderMarkovData, _SemiMarkovData ]])

    try:
        ret =  obj.build_auxiliary_variable()
    except:
        ValueError("build_axiliary_variable failed with the given object.Check the types")
    return ret
