__doc__ = """ Data transformation functions """
__docformat__ = " restructuredtext "


def Merge(obj, *args):    
    """    
    Merging of objects of the same 'data' type or merging of sample correlation functions.

    Usage
    -----
      * ``Merge(histo1, histo2,...)``
      * ``Merge(vec1, vec2,...)``
      * ``Merge(timev1, timev2,...)``
      * ``Merge(seq1, seq2,...)``
      * ``Merge(discrete_seq1, discrete_seq2,...)``
      * ``Merge(top1, top2,...)``
      * ``Merge(correl1, correl2,...)``

    Parameters
    ----------
      * histo1, histo2, ... (_Histogram, _MixtureData, _ConvolutionData, _CompoundData),
      * vec1, vec2, ... (_Vectors),
      * timev1, timev2, ... (_TimeEvents, _RenewalData),
      * seq1, seq2, ... (_Sequences, _DiscreteSequences, _MarkovData, _SemiMarkovData),
      * top1, top2, ... (_Tops),
      * correl1, correl2, ... (_Correlation). 

    Return
    ------
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

    See Also
    --------
    `Cluster`, `Shift`, `Transcode`, `ValueSelect`, `MergeVariable`, `SelectIndividual`, 
    `SelectVariable`, `NbEventSelect`, `TimeScaling`, `TimeSelect`, `AddAbsorbingRun`, 
    `Cumulate`, `Difference`, `IndexExtract`, `LengthSelect`, `MovingAverage`, 
    `RecurrenceTimeSequences`, `RemoveRun`, `Reverse`, `SegmentationExtract`, `VariableScaling`, 
    `RemoveApicalInternodes` """
    
    return obj.merge(list(args))



def MergeVariable(obj, *args, **kargs):    
    """     
    Merging of variables.

    Usage
    -----
      * MergeVariable ``(vec1, vec2,..., RefSample=2)``
      * MergeVariable ``(seq1, seq2,..., RefSample=2)`` 

    Parameters
    ----------
    vec1, vec2, ... (_Vectors),
    seq1, seq2, ... (_Sequences, _DiscreteSequences, _MarkovData, _SemiMarkovData). 

    Keyworkds
    ---------
      * RefSample (int): reference sample to define individual identifiers 
    (the default: no reference sample).

    Return
    ------
    If the arguments are of type _Vectors and if the number of vectors is the same 
    for each sample, an object of type _Vectors is returned.

    If the arguments are of type _Sequences, _DiscreteSequences, _MarkovData, _SemiMarkovData,
    if all the variables are of type STATE, and if the number and the lengths of sequences 
    are the same for each sample, an object of type _Sequences or _DiscreteSequences is returned.
    
    The returned object is of type _DiscreteSequences if all the variables are of type STATE, 
    if the possible values for each variable are consecutive from 0 and if the number of 
    possible values for each variable is < 15. 

    See Also
    --------
    `Merge`, `Cluster`, `Shift`, `Transcode`, `ValueSelect`, `SelectIndividual`, `SelectVariable`, 
    `AddAbsorbingRun`, `Cumulate`, `Difference`, `IndexExtract`, `LengthSelect`, 
    `MovingAverage`, `RecurrenceTimeSequences`, `RemoveRun`, `Reverse`, 
    `SegmentationExtract`, `VariableScaling`.
    """

    RefSample = kargs.get("RefSample", -1)
    
    return obj.merge_variable(list(args), RefSample)



################################################################################

def ExtractData(model):
    """    
    Extraction of the 'data' part of an object of type 'model'.
    
    Usage
    -----
      * ``ExtractData(mixt)``
      * ``ExtractData(convol)``
      * ``ExtractData(compound)``

      * ``ExtractData(hmc)``
      * ``ExtractData(hsmc)`` 

    Parameters
    ----------
      * mixt (_Mixture),
      * convol (_Convolution),
      * compound (_Compound),

      * hmc (_HiddenMarkov),
      * hsmc (_HiddenSemiMarkov). 

    Return
    ------
    If mixt contains a 'data' part, an object of type `_MixtureData` is returned.
    If convol contains a 'data' part, an object of type `_ConvolutionData is returned.
    If compound contains a 'data' part, an object of type `_CompoundData is returned.
    If hmc contains a 'data' part, an object of type `_MarkovData` is returned.
    If hsmc contains a 'data' part, an object of type `_SemiMarkovData` is returned.

    Description
    -----------
    This function enables to extract the 'data' part of an object of type 'model' 
    when the estimation of model parameters from data gives rise to the construction 
    of pseudo-data. This situation is notably exemplified by the computation of 
    optimal state sequences from estimated hidden Markovian processes (optional 
    argument StateSequences of the function Estimate set at "ForwardBackward" or 
    "Viterbi"). 

    See Also
    --------
    `Estimate`
    """

    return model.extract_data()



def ExtractDistribution(model, *args):
    """    
    Extraction of a distribution from an object of type 'model'.

    Usage
    -----
      * ``ExtractDistribution(mixt, mixt_type)``
      * ``ExtractDistribution(mixt, "Component", index)``
      * ``ExtractDistribution(convol, "Elementary", index)``
      * ``ExtractDistribution(convol, "Convolution")``
      * ``ExtractDistribution(compound, compound_type)``

      * ``ExtractDistribution(renew, renew_type)``
      * ``ExtractDistribution(renew, "NbEvent", time)``

      * ``ExtractDistribution(markov, markov_type, state)``
      * ``ExtractDistribution(markov, markov_type, variable, output)``

      * ``ExtractDistribution(top_param, position)`` 
      
    Parameters
    ----------
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

    Return
    ------
    If the arguments (mixt_type, index, compound_type, renew_type, time, 
    markov_type, state, variable, output, position) defined an existing 
    distribution, an object of type `_Distribution` is returned. 

    See Also
    --------
    `Plot`, `Fit`, `Simulate`. 
    """
    
    return Extract(model, *args)
    

def ExtractHistogram(data, *args):
    """  
    Extraction of a frequency distribution from an object of type 'data'.

    Usage
    -----
      * ``ExtractHistogram(mixt_histo, mixt_type)``
      * ``ExtractHistogram(mixt_histo, "Component", index)``
      * ``ExtractHistogram(convol_histo, "Elementary", index)``
      * ``ExtractHistogram(convol_histo, "Convolution")``
      * ``ExtractHistogram(compound_histo, compound_type)``

      * ``ExtractHistogram(vec1)``
      * ``ExtractHistogram(vecn, variable)``

      * ``ExtractHistogram(renewval_data, renew_type)``
      * ``ExtractHistogram(timev, timev_type)``
      * ``ExtractHistogram(timev, "NbEvent", time)``

      * ``ExtractHistogram(seq, "Length")``
      * ``ExtractHistogram(seq, "Value")``
      * ``ExtractHistogram(seq, "Value", variable)``
      * ``ExtractHistogram(discrete_seq1, seq_type, value)``
      * ``ExtractHistogram(discrete_seqn, seq_type, variable, value)``
      * ``ExtractHistogram(simul_seq1, "Observation", value)``
      * ``ExtractHistogram(simul_seq1, "Observation", variable, value)``

      * ``ExtractHistogram(tops, "Main")``
      * ``ExtractHistogram(top, "NbAxillary", position)``
      
      
    Parameters
    ----------
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

      ...

    Return
    ------
    If the arguments (mixt_type, index, compound_type, renew_type, time, 
    markov_type, state, variable, output, position) defined an existing frequency
    distribution, an object of type `_Histogram` is returned. 

    See Also
    --------
    `Plot`, `Fit`, `Simulate`. 
    """
    
    return Extract(data, *args)


def Extract(obj, *args):
    """ 
    Common method to redirect extract function call
    See`ExtractHistogram` or `ExtractDistribution`
    """

    func_map = { 
        "Weight" : "extract_weight",
        "Component"  : "extract_component",
        "Mixture" : "extract_mixture",
        "Elementary" : "extract_elementary",
        "Convolution" : "extract_convolution",
        "Sum" : "extract_sum",
        "Compound" : "extract_compound",
        }


    # Test without "string" command (ex for _Vecotors)
    if(len(args) == 0 or not isinstance(args[0], str)):

        # _Vectors with one variable
        try:
            if(obj.get_nb_variable() == 1):
                args.append(1)
        except:
            pass

        try:
            return obj.extract(*args)

        except AttributeError:
            raise TypeError("Expect an extract command as first argument." + \
                                "Possible command are : %s"%(str(func_map.keys())))

    # Others cases
    key = args[0]
           
    try:
        func_name = func_map[key]
        f = getattr(obj, func_name)

    except KeyError:
        raise AttributeError("Object has no attribute '%s'"%(key))

    except AttributeError:
        raise

    # Call function
    return f(*args[1:])


################################################################################



# Utility function

def __get_mode__(kargs):
    """ Return True if kargs has "keep" for the "mode" key """

    mode = kargs.get("mode", None)
    if(not mode): kargs.get("Mode", None)
    if(mode == "Keep") : keep = True
    else : keep = False

    return keep




def SelectVariable(obj, variables, Mode="Keep"):
    """ 
    Selection of variables.

    Usage
    -----

      * ``SelectVariable(vec, variable, Mode="Reject")``
      * ``SelectVariable(vec, variables, Mode="Reject")``

      * ``SelectVariable(seq, variable, Mode="Reject")``
      * ``SelectVariable(seq, variables, Mode="Reject")``

    Parameters
    ----------
      * vec (vectors),

      * seq (sequences, discrete_sequences, markov_data, semi-markov_data),

      * variable (int): variable index.
      * variables (array(int)): variable indices. 

    Keywords
    --------
      * Mode (string): conservation or rejection of the selected variables: "Keep" (default) or "Reject". 

    Return
    ------
      If either variable or variables[1], ..., variables[n] are valid indices of variables, 
      an object of type vectors (respectively sequences or discrete_sequences) is returned, 
      otherwise no object is returned. In the case of a first argument of type sequences, 
      the returned object is of type discrete_sequences if all the variables are of type STATE, 
      if the possible values for each variable are consecutive from 0 and if the number of 
      possible values for each variable is < 15. 

    See Also
    --------
      `AddAbsorbingRun`, `Cluster`, `Cumulate`, `Difference`, `IndexExtract`, `LengthSelect`, 
      `Merge`, `MergeVariable`, `MovingAverage`, `RecurrenceTimeSequences`, `RemoveRun`, 
      `Reverse`, `SelectIndividual`, `Shift`, `Transcode`, `ValueSelect`, 
      `SegmentationExtract`, `VariableScaling`. 
    """

    keep = bool(Mode == "Keep")

    # Test if variables is a list
    try:
        v = variables[0]
    except TypeError:
        variables = [variables,]

    return obj.select_variable(list(variables), keep)


def SelectIndividual(obj, identifiers, Mode="Keep"):
    """    
    Selection of vectors, sequences, tops or patterns (in a dissimilarity matrix).

    Usage
    -----
      * ``SelectIndividual(vec, identifiers, Mode="Reject")``
      * ``SelectIndividual(seq, identifiers, Mode="Reject")``
      * ``SelectIndividual(top, identifiers, Mode="Reject")``
      * ``SelectIndividual(dist_matrix, identifiers, Mode="Reject")`` 

    Parameters
    ----------
      vec (vectors),

      seq (sequences, discrete_sequences, markov_data, semi-markov_data),

      top (tops),

      dist_matrix (distance_matrix),
 
      identifiers (array(int)): identifiers. 

    Keywords
    --------
       Mode (string): conservation or rejection of the selected individuals: "Keep" (default) or "Reject". 

    Returned
    --------
    If identifiers[1], ..., identifiers[n] are valid identifiers of vectors (respectively 
    sequences, tops or patterns compared in a dissimilarity matrix), an object of type vectors 
    (respectively sequences or discrete_sequences, tops or distance_matrix) is returned, 
    otherwise no object is returned. In the case of a first argument of type sequences, 
    discrete_sequences, markov_data, semi-markov_data, the returned object is of type 
    discrete_sequences if all the variables are of type STATE, if the possible values for 
    each variable are consecutive from 0 and if the number of possible values for each variable 
    is < 15. 

    See Also
    --------
    `Cluster`, `Merge`, `Shift`, `Transcode`, `ValueSelect`, `MergeVariable`, `SelectVariable`, 
    `AddAbsorbingRun`, `Cumulate`, `Difference`, `IndexExtract`, `LengthSelect`, 
    `MovingAverage`, `RecurrenceTimeSequences`, `RemoveSeries`, `Reverse`, `SegmentationExtract`, 
    `VariableScaling`, `RemoveApicalInternodes`, `Symmetrize`. 
    """

    keep = bool(Mode == "Keep")

    return obj.select_individual(identifiers, keep)

    

    
def ValueSelect(obj, *args, **kargs):
    """ 
    Selection of individuals according to the values taken by a variable

    Usage
    -----
      * ``ValueSelect(histo, value, Mode="Reject")``
      * ``ValueSelect(histo, min_value, max_value, Mode="Reject")``

      * ``ValueSelect(vec1, value, Mode="Reject")``
      * ``ValueSelect(vec1, min_value, max_value, Mode="Reject")``
      * ``ValueSelect(vecn, variable, value, Mode="Reject")``
      * ``ValueSelect(vecn, variable, min_value, max_value, Mode="Reject")``

      * ``ValueSelect(seq1, value, Mode="Reject")``
      * ``ValueSelect(seq1, min_value, max_value, Mode="Reject")``
      * ``ValueSelect(seqn, variable, value, Mode="Reject")``
      * ``ValueSelect(seqn, variable, min_value, max_value, Mode="Reject")``

    Parameters
    ----------
      * histo (histogram, mixture_data, convolution_data, compound_data),

      * value (int): value,
      * min_value (int): minimum value,
      * max_value (int): maximum value,

      * vec1 (vectors): values,
      * vecn (vectors): vectors,
      * variable (int): variable index,

      * seq1 (sequences, discrete_sequences, markov_data, semi-markov_data): univariate sequences,
      * seqn (sequences, discrete_sequences, markov_data, semi-markov_data): multivariate sequences. 
    Keywords
    --------
      * Mode (string): conservation or rejection of selected individuals: "Keep" (the default) or "Reject". 
   
    Return
    ------
    
    If either value 0 or if 0 < min_value < max_value and if the range of values defined either 
    by value or by min_value and max_value enables to select individuals, an object of 
    type HISTOGRAM is returned (respectively vectors, sequences or discrete_sequences), 
    otherwise no object is returned. In the case of a first argument of type sequences, 
    discrete_sequences, markov_data or semi-markov_data, the returned object is of type 
    discrete_sequences if all the variables are of type STATE, if the possible values for each 
    variable are consecutive from 0 and if the number of possible values for each variable is < 15. 

    See Also
    --------
    `Cluster`, `Merge`, `Shift`, `Transcode`, `SelectIndividual`, `MergeVariable`, 
    `SelectVariable`, `Cumulate`, `Difference`, `IndexExtract`, `LengthSelect`, 
    `MovingAverage`, `RecurrenceTimeSequences`, `RemoveRun`, `Reverse`, 
    `SegmentationExtract`, `VariableScaling`. 
    """

    keep = __get_mode__(kargs)

    # Test for vectors
    try:
        nb_variable = obj.get_nb_variable()
        variable = 1
    except AttributeError:
        nb_variable = 0


    # Parse args
    l = len(args)
    
    if(l == 3):
        variable, min, max = args

    elif(l == 2):
        # 2 cases (min_value, max_value) or (variable, value)
        if(nb_variable):
            variable, min = args
            max = min
        else:
           min, max = args

    elif(l == 1):
        v = args[0]
        if(isinstance(v, tuple) and len(v) == 2):
            min, max = v
        else:
            min = max = v

    # Check variable
    if(nb_variable and variable > len(obj)):
        raise ValueError("Variable is greater than object size")

    
    if(nb_variable):    # Vectors
        return obj.value_select(variable, min, max, keep)
    else:
        return obj.value_select(min, max, keep)




################################################################################


    

def Shift(obj, *args):
    """ 
    Shifting of values
    
    Usage
    -----
      * ``Shift(histo, param)``
      * ``Shift(vec1, param)``
      * ``Shift(vecn, variable, param)``
      * ``Shift(seq1, param)``
      * ``Shift(seqn, variable, param)``

    Arguments
    ---------
      * histo (histogram, mixture_data, convolution_data, compound_data),

      * param (int): shifting parameter,

      * vec1 (vectors): values,
      * vecn (vectors): vectors,
      * variable (int): variable index,

      * seq1 (sequences): univariate sequences,
      * seqn (sequences): multivariate sequences. 

    Return
    ------
      If the shifting makes that the lower bound to the possible values is positive, an 
      object of type HISTOGRAM (respectively _Vectors, _Sequences) is returned. In the 
      case of a first argument of type sequences, the returned object is of type 
      discrete_sequences if all the variables are of type STATE, if the possible values 
      for each variable are consecutive from 0 and if the number of possible values for 
      each variable is 15. 

    See Also
    --------
      `Cluster`, `Merge`, `Transcode`, `MergeVariable`, `SelectIndividual`, `SelectVariable`, 
      `AddAbsorbingRun`, `Cumulate`, `Difference`, `Lengthselect`, `MovingAverage`, `IndexExtract`, 
      `RecurrenceTimeSequences`, `RemoveRun`, `Reverse`, `SegmentationExtract`, `ValueSelect`, 
      `VariableScaling`.
    """
    
    try:
        # Try vector.shift
        nb_var = obj.get_nb_variable()

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

    Usage
    -----
      * ``Fit(histo, dist)``

    Arguments
    ---------
      * histo (histogram, mixture_data, convolution_data, compound_data),
      * dist (distribution, mixture, convolution, compound). 

    Return
    ------
       Distribution

    Description
    -----------
       The result is displayed in the shell window (characteristics of the frequency 
       and theoretical distributions, log-likelihood of the data for the theoretical distribution, 
       information - maximum possible log-likelihood of the data -, c2 goodness of fit test). 

    Background
    ----------
       The difference between the information measure and the log-likelihood is the 
       Kullback-Leibler divergence from the observed distribution to the theoretical distribution. 
       It is also one-half the deviance of the theoretical distribution.

       Assume that a sample of size n is generated by a given random variable. 
       The statistic measures the random deviation between the observed frequencies fi and the 
       theoretical frequencies npi:

    	   with   	

       If each theoretical frequency npi is greater than a given threshold (between 1 and 5 
       according to the authors), has a c2 with k - 1 degrees of freedom. 
    """

    return histo.fit(dist)


################################################################################


################################################################################
from openalea.stat_tool import get_test_prefix


class TestShift:

    def test_shift_histo(self):

        from histogram import Histogram 

        h = Histogram(get_test_prefix() + "meri2.his")
        assert h.shift(2)

        assert Shift(h, 2) == h.shift(2)


    def test_shift_vector(self):

        from vectors import Vectors

        vn = Vectors([[0.2, 1., 2, 3], [4.2, 5, 6, 7]])
        v1 = Vectors([[1.]])

        assert Shift(v1, 2)
        assert Shift(vn, 1, 2)



    def test_shift_sequence(self):
        raise NotImplementedError()


    def test_fit(self):

        from distribution import Distribution
        from histogram import Histogram

        meri5 = Histogram(get_test_prefix() + "meri5.his")
        dist5 = Fit(meri5, Distribution("B", 0, 10, 0.437879))
        assert dist5



class TestSelect:
    def test_value_select_float(self):
    
        from openalea.stat_tool import Vectors

        v = Vectors([[0.1, 0.3, 4.2], 
                     [0.5, 2.3, 1.2],
                     [4.5, 6.3, 3.2],
                     ])

        v1b = ValueSelect(v, 1, 0.2, 2.0, Mode="keep")
        v2 =  ValueSelect(v, 2, 1.0, 6.0, Mode="keep")
        v3 =  ValueSelect(v, 3, 1.0, 2.0, Mode="keep")
        
        assert v and v1b and v2 and v3
        print len(v1b)
        assert len(v1b) == 2
        


    def test_value_select_int(self):

        from openalea.stat_tool import Vectors

        v = Vectors([[1, 3, 4], 
                     [4, 12, 2],
                     [8, 7, 3],
                     ])

        v1b =  ValueSelect(v, 1, 2, 8, Mode="Reject")

        assert v1b

        try:
            v1 =  ValueSelect(v, 1, 2, 3, Mode="Keep")
            assert False
        except Exception:
            # Empty sample
            assert True
    


    def test_select_variable(self):

        from openalea.stat_tool import Vectors

        a = [[1, 3, 4], 
             [4, 12, 2],
             [8, 7, 3],]

        v = Vectors(a)

        for i in range(3):
            v1 =  SelectVariable(v, i+1, Mode="Keep")
            assert v1
            assert len(v1) == 3
            assert len(v1[0]) == 1

            for j in range(3):
                assert v1[j][0] == a[j][i]


    def test_select_individual(self):
        
        from openalea.stat_tool import Vectors

        a = [[1, 3, 4], 
             [4, 12, 2],
             [8, 7, 3],]

        v = Vectors(a)
        selection = SelectIndividual(v, [1,2], Mode="Keep")
        
        assert len(selection) == 2
        print selection[0] == [1, 3, 4]


class TestExtract:

    def test_build_mixture(self):

        from distribution import Uniform
        from mixture import Mixture

        d1 = Uniform(0, 10)
        d2 = Uniform(10, 20)
        d3 = Uniform(20, 30)

        m = Mixture(0.1, d1, 0.2, d2, 0.7, d3)
        assert m
        return m


    def test_build_convolution(self):

        from distribution import Binomial
        from distribution import NegativeBinomial
        from convolution import Convolution

        d1 = Binomial(0, 10, 0.5)
        d2 = NegativeBinomial(0, 1, 0.1)

        m = Convolution(d1, d2)
        assert m
        return m


    def test_extract_data(self):

        from histogram import Histogram

        h = Histogram(get_test_prefix() + "meri2.his")
        mixt = h.estimate_mixture(["B", "NB"])

        assert ExtractData(mixt)
        #assert Convolution().extract_data()
        #assert Compound().extract_data()


    def test_extract_mixture_distribution(self):

        mixt = self.test_build_mixture()

        assert ExtractDistribution(mixt, "Weight")
        assert ExtractDistribution(mixt, "Mixture")
        assert ExtractDistribution(mixt, "Component", 1)
        assert ExtractDistribution(mixt, "Component", 2)
        assert ExtractDistribution(mixt, "Component", 3)

        try:
            ExtractDistribution(mixt, "Component", 0)
            assert False
        except: # Bas distrubition index
            assert True


    def test_extract_convolution_distribution(self):

        convol = self.test_build_convolution()

        assert ExtractDistribution(convol, "Convolution")
        assert ExtractDistribution(convol, "Elementary", 1)
        assert ExtractDistribution(convol, "Elementary", 2)


    def test_extract_histogram(self):

        from histogram import Histogram

        h = Histogram(get_test_prefix() + "meri2.his")
        mixt = h.estimate_mixture(["B", "NB"])

        assert ExtractHistogram(mixt, "Weight")
        assert ExtractHistogram(mixt, "Mixture")
        assert ExtractHistogram(mixt, "Component", 1)
        assert ExtractHistogram(mixt, "Component", 2)

        try:
            ExtractHistogram(mixt, "Component", 3)
            assert False
        except: # Bas distrubition index
            assert True


class TestMerge:
    
    def test_merge(self):
        
        from mixture import Mixture
        from distribution import Distribution
        from simulate import Simulate

        mixt1 = Mixture(0.6, Distribution("B", 2, 18, 0.5), 
                        0.4, Distribution("NB", 10, 10, 0.5))

        mixt_histo1 = Simulate(mixt1, 200)

        histo10 = mixt_histo1.extract_component(1)
        histo11 = mixt_histo1.extract_component(2)

        histo12 = Merge(histo10, histo11)

        assert histo12

    def test_merge_histo(self):
        from histogram import Histogram
        meri1 = Histogram(get_test_prefix() + "meri1.his")
        meri2 = Histogram(get_test_prefix() + "meri2.his")
        meri3 = Histogram(get_test_prefix() + "meri3.his")
        meri4 = Histogram(get_test_prefix() + "meri4.his")
        meri5 = Histogram(get_test_prefix() + "meri5.his")

        meri = Merge(meri1, meri2, meri3, meri4, meri5)
        assert meri


    def test_merge_vectors(self):
        
        from vectors import Vectors
        a = [[1, 3, 4], 
             [4, 12, 2],
             [8, 7, 3],]

        b = [[2, 78, 45], 
             [6, 2, 122],
             [3, 4, 31],]


        v1 = Vectors(a)
        v2 = Vectors(b)

        v = Merge(v1, v2)
        assert v


    def test_mergevariable(self):
        
        from vectors import Vectors
        a = [[1, 3, 4], 
             [4, 12, 2],
             [8, 7, 3],]

        v = Vectors(a)

        v1 =  SelectVariable(v, 1, Mode="Keep")
        v2 =  SelectVariable(v, 2, Mode="Keep")
        v3 =  SelectVariable(v, 3, Mode="Keep")

        merged = MergeVariable(v1,v2,v3)

        for i in range(3):
            for j in range(3):
                assert merged[i][j] == v[i][j]
        

    


