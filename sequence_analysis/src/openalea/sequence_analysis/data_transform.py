"""data transform"""
__revision__ = "$Id: vectors.py 6217 2009-04-08 12:40:15Z cokelaer $"


#import os
import openalea.stat_tool.interface as interface
import _sequence_analysis


#import openalea.sequence_analysis._sequence_analysis as _sequence_analysis


def __get_mode__(kargs):
    """ Return True if kargs has "keep" for the "mode" key """

    mode = kargs.get("mode", None)
    if(not mode):
        mode = kargs.get("Mode", None)
    if(mode == "Keep" or mode == "keep") : keep = True
    else : keep = False

    return keep


def check_nb_variable(obj, variable):
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

	RemoveRun(seq1, value, position, MaxLength->4)
	RemoveRun(seqn, variable, value, position, MaxLength->4)	

    :Arguments:
	    * seq1 (sequences, discrete_sequences, markov_data, semi-markov_data): univariate sequences,
	    * seqn (sequences, discrete_sequences, markov_data, semi-markov_data): multivariate sequences,
	    * variable (int): variable index, (optional if only one variable 
	    * value (int): value,
	    * position (string): position of the removed run: "Begin" (or "b") or 
	      "End" (or "e").
	
    :Optional Arguments:
    	* MaxLength (int): maximum length of the removed run (default behaviour: 
    	  the entire run is removed).
	
    :Returned Object:

	If variable is a valid index of variable, if value is a possible value and
	 if MaxLength > 0, an object of type sequences or discrete_sequences is 
	 returned, otherwise no object is returned. The returned object is of type
	  discrete_sequences if all the variables are of type STATE, if the possible 
	  values for each variable are consecutive from 0 and if the number of 
	  possible values for each variable is <= 15.

    :Examples:
    
    >>> RemoveRun(seq1, 0, "End") 	
    >>> RemoveRun(seq5, 2, 0, "End")

    :See Also:
	
	AddAbsorbingRun,
	Cluster, 
	Cumulate, 
	Difference, 
	IndexExtract, 
	LengthSelect, 
	Merge, 
	MergeVariable, 
	MovingAverage, 
	RecurrenceTimeSequences, 
	Reverse, 
	SegmentationExtract, 
	SelectIndividual, 
	SelectVariable, 
	Shift, 
	Transcode, 
	ValueSelect,
	VariableScaling.
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
    if position == 'End': position = 'e'
    if position == 'Begin': position = 'b'
    if position not in ['b','e']:
        raise TypeError("position must be 'e' for end or 'b' for begin")
    if value is None:
        raise TypeError("value must be provided")

   # variable argument
    variable = check_nb_variable(obj, variable)
    
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
    
    ExtractHistogram, 
    Plot, 
    MergeVariable, 
    Compare (renewal process), 
    ContingencyTable, 
    Regression, 
    Varianceanalysis.
    """
    
    type_dict = {
        "Length"            : _sequence_analysis.LENGTH,
        "NbRun"             : _sequence_analysis.NB_RUN,
        "NbOccurrence"      : _sequence_analysis.NB_OCCURRENCE,
        "FirstOccurrence"   : _sequence_analysis.FIRST_OCCURRENCE,
        "Mean"              : _sequence_analysis.SEQUENCE_MEAN,
        "Cumul"             : _sequence_analysis.SEQUENCE_CUMUL
        }        

    #.value("OBSERVATION",OBSERVATION)
    #.value("RECURRENCE_TIME",RECURRENCE_TIME)
    #.value("SOJOURN_TIME",SOJOURN_TIME)
    #.value("INITIAL_RUN",INITIAL_RUN)
    #.value("FINAL_RUN",FINAL_RUN)
    
    if key not in type_dict.keys():
        raise TypeError("key must be in %s" % type_dict.keys())
    
    key = type_dict[key]
    
    variable = -1 
    
    
    if len(args)==2:
        variable = args[0]
        value = args[1]
        sequence = obj.extract_vectors(key, variable, value)
    elif len(args)==1 and obj.nb_variable==1:
        value = args[0]
        sequence = obj.extract_vectors(key, 1,value)
    elif len(args)==0:
        # why the overloading in export_markovian does not work ? 
        #here variable is 1 and whatever value given, the results is the same
        sequence = obj.extract_vectors(key,1,0)
    else:
        raise TypeError("nb_variable =1 so varaible must be 1")
    
    #sequence = obj.extract_vectors(key, variable, value)
    return sequence
#.markovian_sequences()

   


def SegmentationExtract(obj, variable, values , Mode="Keep"):
    """SegmentationExtract 
    
    Extraction by segmentation of sub-sequences.
  
    :Usage
    SegmentationExtract(seqn, variable, value, Mode->"Reject")
    SegmentationExtract(seqn, variable, values, Mode->"Reject")    
  
    :Arguments:
    seqn (sequences, discrete_sequences, markov_data, semi-markov_data): multivariate sequences,
    variable (int): variable index,
    value (int): value,
    values (ARRAY(int)): values.
    
    :Optional Arguments: 
    Mode (string): conservation or rejection of the selected sub-sequences: "Keep" (the default) or "Reject".
    
    :Returned Object:
    If all the variables are of type STATE, if variable is a valid index of variable and if either value or values@1, ..., values@n are possible values, an object of type sequences or discrete_sequences is returned, otherwise no object is returned. The returned object is of type discrete_sequences if all the variables are of type STATE, if the possible values for each variable are consecutive from 0 and if the number of possible values for each variable is <= 15.
    
    :Description:
    Sub-sequences corresponding to run of value or values@1, ..., values@n are extracted (or to all the possible values except value or values@1, ..., values@n) are extracted.
    
    .. seealso::
    
    AddAbsorbingRun,
    Cluster, 
    Cumulate, 
    Difference, 
    IndexExtract, 
    LengthSelect, 
    Merge, 
    MergeVariable, 
    MovingAverage, 
    RecurrenceTimeSequences, 
    RemoveRun, 
    Reverse, 
    SelectIndividual, 
    SelectVariable, 
    Shift, 
    Transcode, 
    ValueSelect,
    VariableScaling.
"""

    keep = bool(Mode == "Keep" or Mode == "keep")


    # Test if variables is a list
    try:
        v = values[0]
    except TypeError:
        values = [values,]

    sequence = obj.segmentation_extract(variable, list(values), keep)
    return sequence.markovian_sequences()

def LengthSelect(obj, minLength, *args, **kargs):
    """LengthSelect
    
    Selection of sequences according to a length criterion.
    
    :Usage:
    >>> LengthSelect(seq, length, Mode->"Reject")
    >>> LengthSelect(seq, min_length, max_length, Mode->"Reject")    
  
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
    
    AddAbsorbingRun,
    Cluster, 
    Cumulate, 
    Difference, 
    IndexExtract, 
    Merge, 
    MergeVariable, 
    MovingAverage, 
    RecurrenceTimeSequences, 
    RemoveRun, 
    Reverse, 
    SegmentationExtract, 
    SelectIndividual, 
    SelectVariable, 
    Shift, 
    Transcode, 
    ValueSelect,
    VariableScaling.
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
    Usage
    RecurrenceTimeSequences(seq1, value)
    RecurrenceTimeSequences(seqn, variable, value)    
    Arguments
    seq1 (sequences, discrete_sequences, markov_data, semi-markov_data): univariate sequences,
    seqn (sequences, discrete_sequences, markov_data, semi-markov_data): multivariate sequences,
    variable (int): variable index,
    value (int): value.
    
    Returned Object
    If the selected variable is of type STATE and if value is a possible value, an object of type sequences or discrete_sequences is returned, otherwise no object is returned. The returned object is of type discrete_sequences if all the variables are of type STATE, if the possible values for each variable are consecutive from 0 and if the number of possible values for each variable is <= 15.
    
    See Also
    
    AddAbsorbingRun,
    Cluster, 
    Cumulate, 
    Difference, 
    IndexExtract, 
    LengthSelect, 
    Merge, 
    MergeVariable, 
    MovingAverage, 
    RemoveRun, 
    Reverse, 
    SegmentationExtract, 
    SelectIndividual, 
    SelectVariable, 
    Shift, 
    Transcode, 
    ValueSelect,
    VariableScaling.
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
    """
    *args = variable, value,word_length
    **kargs = begin state, end_state, min_frequency
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

      Usage
    
      AddAbsorbingRun(seq, Length->50)
    
      Arguments
    seq (distance_sequences, markov_data, semi-markov_data)
    
      Optional Arguments 
    Length (int): length of the sequences. A default value is computer from the maximum sequence length.
    
      Returned Object
    An object of type discrete_sequences is returned.

    See Also
    
    Cluster, 
    Cumulate, 
    IndexExtract,
    LengthSelect,
    Merge,
    MergeVariable,
    MovingAverage,
    RecurrenceTimeSequences,
    Reverse,
    SegmentationExtract,
    SelectIndividual,
    SlectVariable,
    Shift,
    Transcode,
    ValueSelect,
    VariableScaling.
    """
    #todo possibly : set default Length to -1 and add SequenceLength or RunLength optional arguments
    return obj.add_absorbing_run(SequenceLength, RunLength)

def Reverse(obj):
    """
    Reverse
    Reversing of sequences or 'tops'.
    
    :Usage:
    
    Reverse(seq)
    Reverse(discrete_seq)
    
    Reverse(top)    
      
    :Arguments:
    
    seq (sequences),
    discrete_seq (discrete_sequences, markov_data, semi-markov_data),
    
    top (tops).
    
  Returned Object
    If the argument is of type sequences, an object of type sequences is returned. If the argument is of type discrete_sequences, markov_data, semi-markov_data, an object of type discrete_sequences is returned. If the argument is of type tops, an object of type tops is returned.
    
  See Also
    
    AddAbsorbingRun,
    Cluster, 
    Cumulate, 
    Difference, 
    IndexExtract, 
    LengthSelect, 
    Merge, 
    MergeVariable, 
    MovingAverage, 
    RecurrenceTimeSequences, 
    RemoveRun, 
    RemoveApicalInternodes, 
    SegmentationExtract, 
    SelectIndividual, 
    SelectVariable, 
    Shift, 
    Transcode, 
    ValueSelect,
    VariableScaling.
    """
    
    ret = obj.reverse()
    
    try:
        return ret.markovian_sequences()
    except TypeError:
        return ret

    
def Thresholding(obj, MinProbability=1e-5):
    
    return obj.thresholding(MinProbability)
    
    
def Cumulate(obj, Variable=1):
    #todo check that it is 
    #SEQUENCES :, MARKOVIAN_SEQUENCES,VARIABLE_ORDER_MARKOV_DATA, SEMI_MARKOV_DATA,NONHOMOGENEOUS_MARKOV_DATA :    
    return obj.cumulate(Variable)
    
def Difference(obj, Variable=1, FirstElement=False):
    #SEQUENCES :, MARKOVIAN_SEQUENCES,VARIABLE_ORDER_MARKOV_DATA, SEMI_MARKOV_DATA,NONHOMOGENEOUS_MARKOV_DATA :
    return obj.difference(Variable, FirstElement)
    
    
def IndexParameterExtract(obj, minIndex, maxIndex=40):    
    """cf IndexExtract"""
    return obj.index_parameter_extract(minIndex, maxIndex)
    
def IndexParameterSelect(obj, minIndex, *args, **kargs):    
    """cf IndexExtract"""
    mode = __get_mode__(kargs)
    if len(args) == 0:
        maxIndex = minIndex
    elif len(args)==1:
        maxIndex = args[0]
    else:
        raise KeyError("one or two arguments required and one optional arguments (Mode). see usage")

    return obj.index_parameter_select(minIndex, maxIndex, mode)
    
    
def ComputeStateSequences(obj, data, Characteristics=True):
    
    return data.compute_state_sequences(obj, Characteristics)
    

def MovingAverage(obj, itype, Variable=1, BeginEnd=False, output="Trend"):
    """
    itype is list of vlaues or a distribution
    
    put Moving Average doc here 
    """
    func_map = {
                 "Sequence":0,
                 "Trend":1,
                 "SubtractionResidual":2,
                 "Residual":2,
                 "DivisionResidual":3
                 }
    if output not in func_map.keys():
        raise KeyError("output choice must be in %s" % func_map.keys())
    
    if isinstance(itype, list):
        #todo build a 2*N+1 filter here or inside moving_average function ?  
        return obj.moving_average(itype, Variable, BeginEnd, func_map[output])
    else: #distribution
        return obj.moving_average(itype, Variable, BeginEnd, func_map[output])
        
def ComputeSelfTransition(obj):
    """ """
    obj.self_transition_computation()
    
def TransitionCount(obj, max_order, begin=False, 
                    estimator="MAXIMUM_LIKELIHOOD", filename=None):
    
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
    return obj.cross();
    
    
def PointwiseAverage(obj, *args, **kargs):
    
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
                                 Filename, Format);

    
def ConsecutiveValues(obj, *args, **kargs):    
        
        
    StandardDeviation = kargs.get("AddVariable", False)

    if len(args)==1:
        variable = variable
    else:
        variable = 1
     
    variable = check_nb_variable(obj, variable)
    
    sequence = obj.consecutive_values(variable, AddVariable)
    
    return sequence.markovian_sequences()

    
    
    
def vec2list(vector):
    """
    transform a vector from stat_tool vector into a python list
    """
    output = []
    for i in range(0, vector.nb_vector):
        output.append(vector[i][0])
    return output