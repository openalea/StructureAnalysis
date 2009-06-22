"""Data transform methods

.. author:: Thomas Cokelaer, Thomas.Cokelaer@inria.fr
uthor:/

"""
__revision__ = "$Id: $"

import _sequence_analysis
from openalea.stat_tool._stat_tool import _VectorDistance


def __get_mode__(kargs):
    """ Return True if kargs has "keep" for the "mode" key """

    mode = kargs.get("mode", None)
    if(not mode):
        mode = kargs.get("Mode", None)
    if(mode == "Keep" or mode == "keep"):
        keep = True
    else : keep = False

    return keep


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
        raise TypeError("variable argument must be set because number of variable is larger than 1")
    
    if variable > obj.nb_variable:
        raise TypeError("variable larger than max number of variables %s" % 
                        obj.nb_variable)
        
    return variable
        
        
def RemoveRun(obj, *args, **kargs):
    """RemoveRun

    Removal of the first or last run of a given value (for a given variable) in a sequence.

    :Usage:

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
    
    MaxLength = kargs.get("MaxLength", -1) # max length will be used
    if len(args) == 3:
        variable = args[0]
        value = args[1]
        position = args[2]
    elif len(args) == 2:
        variable = None
        value = args[0]
        position = args[1]
    else:
        raise TypeError("expect 2 or 3 arguments and 1 optional argument (MaxLength)")

    # position argument
    # remove_run method takes as its first argument a char in ['e','b'] that 
    # corresponds to BEGIN and END. The user can use both notation so we need
    # to convert the position to the expected character if neeeded.
    if position == 'End':
        position = 'e'
    if position == 'Begin':
        position = 'b'
    if position not in ['b', 'e']:
        raise TypeError("position must be 'e' for end or 'b' for begin")
    if value is None:
        raise TypeError("value must be provided")

   # variable argument
    variable = _check_nb_variable(obj, variable)
    
    # value must be provided
    if value is None:
        raise TypeError("value must be provided")
    
        
    sequence = obj.remove_run(variable, value, position, MaxLength)
    return sequence.markovian_sequences()


def ExtractVectors(obj, key, *args):
    """ExtractVectors

    Extraction of vectors from global characteristics of sequences (length or counting characteristics).
    
    :Usage:

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
    
    type_dict = {
        "Length"            : _sequence_analysis.LENGTH,
        "NbRun"             : _sequence_analysis.NB_RUN,
        "NbOccurrence"      : _sequence_analysis.NB_OCCURRENCE,
        "FirstOccurrence"   : _sequence_analysis.FIRST_OCCURRENCE,
        "Mean"              : _sequence_analysis.SEQUENCE_MEAN,
        "Cumul"             : _sequence_analysis.SEQUENCE_CUMUL
        }        

    if key not in type_dict.keys():
        raise TypeError("key must be in %s" % type_dict.keys())
     
    
    if obj.nb_variable == 1: 
        variable = 1
    else:
        variable = obj.nb_variable
        
    if key == "Length":
        sequence = obj.extract_vectors(type_dict[key], -1, -1)
    elif key in ["Cumul", "Mean"]:
        sequence = obj.extract_vectors(type_dict[key], variable, -1)
    elif len(args)==2:
        variable = args[0]
        value = args[1]
        sequence = obj.extract_vectors(type_dict[key], variable, value)
    elif len(args)==1 and obj.nb_variable==1:
        value = args[0]
        sequence = obj.extract_vectors(type_dict[key], 1, value)
    else:
        raise TypeError("nb_variable =1 so varaible must be 1")
    
    #sequence = obj.extract_vectors(key, variable, value)
    return sequence
#.markovian_sequences()

   


def SegmentationExtract(obj, variable, values , Mode="Keep"):
    """SegmentationExtract 
    
    Extraction by segmentation of sub-sequences.
  
    :Usage:

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

    keep = bool(Mode == "Keep" or Mode == "keep")


    # Test if variables is a list
    try:
        _v = values[0]
    except TypeError:
        values = [values, ]

    sequence = obj.segmentation_extract(variable, list(values), keep)
    return sequence.markovian_sequences()

def LengthSelect(obj, minLength, *args, **kargs):
    """LengthSelect
    
    Selection of sequences according to a length criterion.
    
    :Usage:

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
    mode = __get_mode__(kargs)
    if len(args) == 0:
        maxLength = minLength
    elif len(args)==1:
        maxLength = args[0]
    else:
        raise KeyError("one or two arguments required and one optional arguments (Mode). see usage")

    sequence =  obj.length_select(minLength, maxLength, mode)
    return sequence.markovian_sequences()


def RecurrenceTimeSequences(obj, *args):
    """RecurrenceTimeSequences

    Computation of recurrence time sequences for a given value (and for a given variable).

    :Usage:

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
    if len(args)==1:
        variable = 1
        value = args[0]
    elif len(args)==2:
        variable = args[0]
        value = args[1]
     
    sequence = obj.recurrence_time_sequences(variable, value)
    return sequence.markovian_sequences()

def WordCount(obj, *args, **kargs):
    """WordCount

    .. note::

        * args = variable, value,word_length
        * kargs = begin state, end_state, min_frequency
    
    .. todo:: the documentation...
    """
    
    if len(args)==1:
        variable = 1
        word_length = args[0]
    elif len(args)==2:
        variable = args[0]
        word_length = args[1]
    begin_state = kargs.get("BeginState", -1)
    end_state = kargs.get("EndState", -1)
    min_frequency = kargs.get("MinFrequency", 1)
    
    return obj.word_count(variable, word_length, begin_state, 
                          end_state, min_frequency)
    
   
def AddAbsorbingRun(obj, SequenceLength=-1, RunLength=-1):
    """AddAbsorbingRun
    
    Addition of a run of absorbing vectors at the end of sequences.

    :Usage:
    
    >>> AddAbsorbingRun(seq, RunLength=20)
    >>> AddAbsorbingRun(seq, SequenceLength=30, RunLength=20)
      
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
    
    .. todo:: wrap the constant values
    """                                  

    MAX_ABSORBING_RUN_LENGTH = 20
        
    try:
        return obj.add_absorbing_run(SequenceLength, RunLength)
    except:
        
        max_max_length = obj.max_length + MAX_ABSORBING_RUN_LENGTH
        raise ValueError("abnormal termination. " 
                         "Check the arguments."
                         "RunLength must be positive and  less than %s" % MAX_ABSORBING_RUN_LENGTH,
                         "SequenceLength must be positive, larger than max_length (%s) and less than max_length + MAX_RUN_LENGTH (%s)" 
                         % (obj.max_length, max_max_length))
        
def Reverse(obj):
    """Reverse

    Reversing of sequences or 'tops'.
    
    :Usage:
    
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
    
    ret = obj.reverse()
    
    try:
        return ret.markovian_sequences()
    except TypeError:
        return ret

    
def Thresholding(obj, MinProbability=1e-5):
    """
    
    .. todo:: documentation
    """
    return obj.thresholding(MinProbability)
    
    
def Cumulate(obj, Variable=1):
    """Cumulate
    
    Sum of successive values along sequences.
  
    :Usage:

    Cumulate(seq, Variable=1)    
    
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
    return obj.cumulate(Variable)
    

def Difference(obj, Variable=-1, FirstElement=False):
    """Difference
    
    First-order differencing of sequences.
  
    :Usage:

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
    return obj.difference(Variable, FirstElement)
    
    
def IndexParameterExtract(obj, minIndex, maxIndex=40):    
    """IndexExtract
    
    Extraction of sub-sequences corresponding to a range of index parameters.
  
    :Usage:
    
    >>> IndexParameterExtract(seq, min_index, MaxIndex=40)    
  
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

    return obj.index_parameter_extract(minIndex, maxIndex)
    
def IndexParameterSelect(obj, minIndex, *args, **kargs):    
    """IndexExtract
    
    Extraction of sub-sequences corresponding to a range of index parameters.

    :Usage:

    >>> IndexParameterSelect(seq, min_index, MaxIndex=40)    

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

    mode = __get_mode__(kargs)
    if len(args) == 0:
        maxIndex = minIndex
    elif len(args)==1:
        maxIndex = args[0]
    else:
        raise KeyError("one or two arguments required and one optional arguments (Mode). see usage")

    return obj.index_parameter_select(minIndex, maxIndex, mode)
    
    
def ComputeStateSequences(obj, data, Characteristics=True):
    """ComputeStateSequences
    
    Computation of the optimal state sequences corresponding to the observed sequences using a hidden Markov chain or a hidden semi-Markov chain.
  
    :Usage:
    
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

    if isinstance(obj, (_sequence_analysis._Markovian_sequences, 
                        _sequence_analysis._Variable_order_markov_data,
                        _sequence_analysis._Semi_markov_data)) and \
            isinstance(data, (_sequence_analysis._Hidden_variable_order_markov, 
                           _sequence_analysis._Hidden_semi_markov
                           )):
            try:            
                return data.state_sequence_computation(obj, Characteristics)
            except:
                raise Exception('call to state_sequence failed')
    else:
        raise Exception("first or second argument has the wrong type")
    

def MovingAverage(obj, itype, Variable=1, BeginEnd=False, Output="Trend"):
    """MovingAverage

    Extraction of trends or residuals using a symmetric smoothing filter.

    :Usage:
    
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
    func_map = {
                 "Sequence":0,
                 "Trend":1,
                 "SubtractionResidual":2,
                 "Residual":2,
                 "DivisionResidual":3
                 }
    if Output not in func_map.keys():
        raise KeyError("output choice must be in %s" % func_map.keys())
    
    if isinstance(itype, list):
        #todo build a 2*N+1 filter here or inside moving_average function ?  
        return obj.moving_average(itype, Variable, BeginEnd, func_map[Output])
    else: #distribution
        return obj.moving_average_from_distribution(itype, Variable, 
                                                    BeginEnd, func_map[Output])
    
    
def ComputeSelfTransition(obj, Order=1):
    """ComputeSelfTransition
    
    Computation of the self-transition probabilities as a function of the index parameter from discrete sequences.
    
    :Usage:

    >>> ComputeSelfTransition(seq, Order=2)    
    
    :Arguments:
    
    * seq (discrete_sequences, markov_data, semi-markov_data).
    
    :Optional Arguments:
    
    * Order (int): Markov chain order (default value: 1).
    
    :Returned Object:
    
    No object returned.
    """
    obj.self_transition_computation()
    
    
def TransitionCount(obj, max_order, begin=False, 
                    estimator="MAXIMUM_LIKELIHOOD", filename=None):
    """.. todo:: documentation"""
    estimator_map = {
                     "MaximumLikelihood": 0, # MAXIMUM_LIKELIHOOD;
                     "Laplace": 1, # LAPLACE;
                     "AdaptativeLaplace": 2, # ADAPTATIVE_LAPLACE;
                     "UniformSubset": 3, # UNIFORM_SUBSET;
                     "UniformCardinality" :4 #estimator = UNIFORM_CARDINALITY;
                     }
    estimator_int = estimator_map[estimator]
    obj.transition_count(max_order, begin, estimator_int, filename)
    
def Cross(obj):    
    """.. todo:: documentation"""
    return obj.cross()
    
    
def PointwiseAverage(obj, *args, **kargs):
    """.. todo:: documentation"""

    SEQUENCE = 0
    SUBTRACTION_RESIDUAL = 1
    STANDARDIZED_RESIDUAL = 2
    
    output_map = {
       "Sequence": SEQUENCE,
       "SubtractionResidual": SUBTRACTION_RESIDUAL,
       "Residual": SUBTRACTION_RESIDUAL,
       "StandardizedResidual":  STANDARDIZED_RESIDUAL
       }
            
    StandardDeviation = kargs.get("StandardDeviation", False)
    Output = kargs.get("Output", "Sequence")
    Output = output_map[Output]
    Filename = kargs.get("Filename", None)
    Format = kargs.get("Format", 'a')
    
    return obj.pointwise_average(StandardDeviation, Output,
                                 Filename, Format)

    
def ConsecutiveValues(obj, *args, **kargs):    
    """ConsecutiveValues"""
        
    AddVariable = kargs.get("AddVariable", False)

    if len(args)==1:
        variable = args[0]
    else:
        variable = 1
     
    variable = _check_nb_variable(obj, variable)
    
    sequence = obj.consecutive_values(variable, AddVariable)
    
    return sequence.markovian_sequences()

def Round(obj,  **kargs):
    """Round"""

    CEIL =  2
    FLOOR = 0
    ROUND = 1
    
    mode_map = {
            "Floor":FLOOR,
            "Ceil":CEIL,
            "Round":ROUND
            }
    Variable = kargs.get("Variable", -1)
    Mode = kargs.get("Mode", "Round")
    Mode = mode_map[Mode]
    
    return obj.round(Variable, Mode) 
    
    
def TimeScaling(obj, scaling_factor=0):
    """TimeScaling
    
    Change of the time unit of data of type {time interval between two observation dates, number of events occurring between these two observation dates}.

    :Usage:

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

    return obj.time_scaling(scaling_factor)    
   
def TimeSelect(obj, *args):
    """TimeSelect
    
    Selection of data item of type {time interval between two observation dates, number of events occurring between these two observation dates} according to a length of the observation period criterion.
  
    :Usage:

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

    if len(args) == 2:
        min_time = args[0]
        max_time = args[1]
    elif len(args) == 1:
        min_time = args[0]
        max_time = min_time
        
    return obj.time_select(min_time, max_time)


def Segmentation(obj, *args, **kargs):
    """Segmentation

    :Usage:
    
    >>> Segmentation(seq, 1,2, VectorDistance("N"), Output="Segment")
    >>> S
    """
    SEQUENCE = 0
    #TREND = 1
    SUBTRACTION_RESIDUAL = 2
    DIVISION_RESIDUAL = 3
    STANDARDIZED_RESIDUAL = 4
    
    output_sequence_map = {
       "Sequence": SEQUENCE,
       "SubtractionResidual": SUBTRACTION_RESIDUAL,
       "Residual": SUBTRACTION_RESIDUAL,
       "StandardizedResidual":  STANDARDIZED_RESIDUAL,
       "DivisionResidual": DIVISION_RESIDUAL
       }
            
    CHANGE_POINT = 0
    SEGMENT = 1
    
    output_type = { 
            "ChangePoint" : CHANGE_POINT,
             "Segment" : SEGMENT
          }
    
    nb_variable = obj.nb_variable

    #from an array, and a string
    if isinstance(args[0], list) and isinstance(args[1], str):
        nb_sequence = obj.nb_sequence
        nb_segment = args[0]
        Model = __get_model__(args[1:], nb_variable)
        
        Output = kargs.get("Output", "Sequence")
        Output = output_sequence_map[Output] 
        identifier = -1
        
        seq = obj.segmentation_array(nb_segment, Model, identifier , Output)
    #from 2 ints, a vectordistance and the output 
    elif isinstance(args[0], int) and  isinstance(args[1], int) \
        and isinstance(args[2], _VectorDistance):
        Output = kargs.get("Output", "Segment")
        Output = output_type[Output] 
        seq = obj.segmentation_vector_distance(args[0], 
                                               args[1], args[2], Output)
    #from 2 ints, and model_type
    elif isinstance(args[0], int) and isinstance(args[1], int)  \
        and isinstance(args[2], str):
        
        nb_sequence = obj.nb_sequence
        
        Model = __get_model__(args[2:], nb_variable)

        Output = kargs.get("Output", "Sequence")
        Output = output_sequence_map[Output]
        
        nb_segment_map = {
                          "Fixed": False,
                          "Estimated":True 
                          }
        NbSegment = kargs.get("NbSegment", "Estimated")
        nb_segment_estimation =  nb_segment_map[NbSegment]
        
        if nb_segment_estimation is False:
            if args[1] == 1 :
                #.. todo:: which call ?????????????????
                seq = obj.segmentation_int_int(args[0] , args[1] ,
                                        Model, Output)
            else:   
                nb_segment = args[0]
                seq = obj.segmentation_array(nb_segment, Model,
                                                args[0], Output)
        else:
            # correct a priori given the aml example seq80
            
            seq = obj.segmentation_model(args[0], args[1], Model)
        
    #from int and a list
    elif isinstance(args[0], int) and isinstance(args[1], list):
        Model = __get_model__(args[2:], nb_variable)
                
        change_point = args[1]
        Output = kargs.get("Output", "Sequence")
        Output = output_sequence_map[Output]
        nb_segment = len(args[1])+1 #.. todo:: to be checked

        seq = obj.segmentation_change_point(args[0] , nb_segment ,
                            change_point , Model , Output)
    return seq
 
 
def SojournTimeSequences(obj, Variable = 1):
    """.. todo:: documentation"""
    
    return obj.sojourn_time_sequences(Variable).markovian_sequences()

def VariableScaling(obj, *args): 
    """VariableScaling
    
    Change of the unit of a variable.
  
    :Usage:
    
    * VariableScaling(seq1, scaling_factor)
    * VariableScaling(seqn, variable, scaling_factor)    
  
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

    if len(args) == 2:
        variable = args[0]
        scaling = args[1]
    elif len(args) == 1:
        variable = 1
        scaling = args[0]
        
    if obj.nb_variable == 1 and variable != 1:
        raise TypeError("nb_variable is 1 but your provided a different value. Check the arguments.")
    return obj.scaling(variable, scaling)
 
def __get_model__(data, nb_variable): 
    """ """
    MULTINOMIAL_CHANGE = 0
    POISSON_CHANGE = 1
    ORDINAL_GAUSSIAN_CHANGE = 2
    GAUSSIAN_CHANGE = 3
    MEAN_CHANGE = 4
    VARIANCE_CHANGE = 5
    MEAN_VARIANCE_CHANGE= 6

    model_type = {
                      "Multinomial":MULTINOMIAL_CHANGE,
                      "Poisson":POISSON_CHANGE,
                      "Ordinal": ORDINAL_GAUSSIAN_CHANGE,
                      "Gaussian": GAUSSIAN_CHANGE,        
                      "Mean": MEAN_CHANGE,
                      "Variance": VARIANCE_CHANGE,
                      "MeanVariance": MEAN_VARIANCE_CHANGE
                      }
    Model = []
    for i in range(0, nb_variable):
        Model.append(model_type[data[i]])

        if (i == 0) and (model_type[data[i]] == MEAN_CHANGE) or \
        (model_type[data[i]] == MEAN_VARIANCE_CHANGE):
        
            for j in range(1, nb_variable):
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
