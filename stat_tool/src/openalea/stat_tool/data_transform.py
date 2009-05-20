__doc__ = """Data transformation functions"""
__revision__ = "$Id$"

def Merge(obj, *args):    
    """    
    Merging of objects of the same 'data' type or merging of sample correlation functions.

    :Parameters:
    
      * histo1, histo2, ... (_Histogram, _MixtureData, _ConvolutionData, _CompoundData),
      * vec1, vec2, ... (_Vectors),
      * timev1, timev2, ... (_TimeEvents, _RenewalData),
      * seq1, seq2, ... (_Sequences, _DiscreteSequences, _MarkovData, _SemiMarkovData),
      * top1, top2, ... (_Tops),
      * correl1, correl2, ... (_Correlation). 

    :Returns:
    
        If the arguments are of type _Histogram, _MixtureData, _ConvolutionData, _CompoundData     
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
    
    return obj.merge(list(args))



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

    RefSample = kargs.get("RefSample", -1)
    
    return obj.merge_variable(list(args), RefSample)



################################################################################

def ExtractData(model):
    """Extraction of the 'data' part of an object of type 'model'.

    This function enables to extract the 'data' part of an object of type 'model' 
    when the estimation of model parameters from data gives rise to the construction 
    of pseudo-data. This situation is notably exemplified by the computation of 
    optimal state sequences from estimated hidden Markovian processes (optional 
    argument StateSequences of the function Estimate set at "ForwardBackward" or 
    "Viterbi"). 

    :Parameters:
      * mixt (_Mixture),
      * convol (_Convolution),
      * compound (_Compound),
      * hmc (_HiddenMarkov),
      * hsmc (_HiddenSemiMarkov).

    :Returns:
      - If mixt contains a 'data' part, an object of type `_MixtureData` is returned.
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

      * mixt (_Mixture),
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
    
      * mixt_histo (_MixtureData),
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


    #arg 0 is misture : Wiegth, Component, Mixture
    #arg0 is convol: Elementary and convolution
    #arg0 is compound : compoubd
    #arg0 is markov: recurrence, sojourn, ...
     
    func_map = { 
        "Weight" : "extract_weight",
        "Component"  : "extract_component",
        "Mixture" : "extract_mixture",
        "Elementary" : "extract_elementary",
        "Convolution" : "extract_convolution",
        "Sum" : "extract_sum",
        "Compound" : "extract_compound",
        "Value": "extract_value",
        "Length":"extract_length"
        }
  
    #todo add a enumerate fro m boost python here
    INTER_EVENT = 0 
    WITHIN_OBSERVATION_PERIOD =1
    LENGTH_BIAS =2
    BACKWARD_RECURRENCE_TIME =3
    FORWARD_RECURRENCE_TIME =4
    NB_EVENT =5
    MIXTURE =6

     
    renewal_nb_event_map = { 
            "InterEvent" : INTER_EVENT,
            "LengthBias" : LENGTH_BIAS,
            "Backward" : BACKWARD_RECURRENCE_TIME,
            "Forward" : FORWARD_RECURRENCE_TIME,
            "Mixture" : MIXTURE,
            "Within": WITHIN_OBSERVATION_PERIOD,
            
            }
    
    NbEvent = kargs.get("NbEvent", NB_EVENT)  
    
    #todo renewal and time events case here 
    
    #top parameters
    if len(args)==1 and isinstance(args[0], int):
        position = args[0]
        return obj.extract(position)
    # vectors case (only 1 compulsary argument that is an int (variable)
    elif (not isinstance(args[0], str)):
        # _Vectors with one variable
        try:
            nb_var = obj.nb_variable
            if (nb_var>1):
                try:
                    variable = args[0]                    
                except IndexError:
                    raise TypeError("Extract with vectors object need 1 arguments (variable) if nb variable>1")
            else:
                variable = 1
                
            return obj.extract(variable)

        except AttributeError:
            raise TypeError("Expect an extract command as first argument." + \
                                "Possible command are : %s"%(str(func_map.keys())))
    else:
    
        # Others cases
        key = args[0]
        
        seq_map = {"Observation":0, 
               "FirstOccurrence":1,
               "Recurrence":2, 
               "Sojourn":3,
               "NbRun":6,
               "NbOccurrence":7,
               "Forward": -1}
        
        
        if key in seq_map:
            f = getattr(obj, "extract")
        else:
            try:
                func_name = func_map[key]
                f = getattr(obj, func_name)
            except KeyError:
            
                raise AttributeError("Object has no attribute '%s'"%(key))


        #case sequences with Value
        if key=="Value":
            if obj.nb_variable == 1:
                variable=1
            else:
                variable = args[1]
            return obj.extract_value(variable) 
        if key=="Length":
            return obj.extract_length()
        #case extract related to semimarkov, vom, hsom, ....

        INITIAL_RUN= 4
        FINAL_RUN= 5
        LENGTH= 8
        SEQUENCE_CUMUL= 9
        SEQUENCE_MEAN= 10
   
        histogram_type = {
                      "FinalRun":5,
                      "InitialRun":4,
                      } 
    
        if key in seq_map.keys():
            if key == "Forward":
                histogram_type = "FinalRun"
                HistogramType = kargs.get("HistogramType", histogram_type)
                HistogramType = histogram_type[HistogramType]
                return f(seq_map[args[0]], args[1], HistogramType)
                    
            if len(args) == 2:
                variable = 1
                value = args[1]
            elif len(args)==3: 
                variable = args[1]
                value = args[2]
            else:
                raise KeyError("expect only 1 or 2 arguments after the type")
                
            return f(seq_map[args[0]], variable, value)
            
    
        
        return f(*args[1:])


################################################################################



# Utility function

def __get_mode__(kargs):
    """ Return True if kargs has "keep" for the "mode" key """

    mode = kargs.get("mode", None)
    if(not mode): 
        mode = kargs.get("Mode", True)
    if(mode == "Keep" or mode == "keep") : keep = True
    else : keep = False

    return keep




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

    keep = bool(Mode == "Keep" or Mode == "keep")
    

    if isinstance(variables, int):
        variables = [variables]
        
    return obj.select_variable(variables, keep)
    
        


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

    keep = bool(Mode == "Keep")

    ret = obj.select_individual(identifiers, keep)
    try:
        # if obj is a sequence, returns markovian_sequences
        return ret.markovian_sequences()
    except:
        return ret
        

def ValueSelect(obj, *args, **kargs):
    """ 
    Selection of individuals according to the values taken by a variable

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
    Mode = kargs.get("Mode", "Keep")
    keep = bool(Mode == "Keep" or Mode == "keep")

    # Test for vectors
    try:
        nb_variable = obj.nb_variable      
    except AttributeError:
        nb_variable = 0

    # Parse args
    l = len(args)
    
    if l == 3 :
        variable, min, max = args

    elif l == 2:
        # 2 cases (min_value, max_value) or (variable, value)
        if(nb_variable):
            variable, min = args
            max = min
        else:
           min, max = args

    elif(l == 1):
        value = args[0]
        if(isinstance(value, tuple) and len(value) == 2):
            min, max = value
        else:
            min = max = value

    # Check variable
    #if(nb_variable and variable > len(obj)):
        #raise ValueError("Variable is greater than object size")

    
    if(nb_variable):    # Vectors, sequences
        return obj.value_select(variable, min, max, keep)
    else:
        return obj.value_select(min, max, keep)




################################################################################


    

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
    
    try:
        # Try vector.shift
        nb_var = obj.nb_variable

        if(nb_var > 1):

            try:
                variable = args[0] 
                param = args[1]
            except IndexError:
                raise TypeError("Shift with vectors object need 2 arguments (variable, param)")
            
        else:
            variable = 1 
            param = args[0]

        return obj.shift(variable, param)


    except AttributeError:
        return obj.shift(*args)




def Fit(histo, dist):
    """
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

           with       

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

    """

    return histo.fit(dist)


    
def Unnormalize(obj):
    
    ret = None
    try:
        ret = obj.unnormalize()
    except:
        pass
    if ret:
        return ret    


    
def Symmetrize(obj):
    
    ret = None
    try:
        ret = obj.symmetrize()
    except:
        pass
    if ret:
        return ret    

def TruncateDistribution(obj, variable):
    """
    .. todo:: to be tests
    
    """
    return obj.truncate(variable)
