#!/usr/bin/env python
#-*- coding: utf-8 -*-
"""Sequences modules

.. module:: sequences
    :synopsis: a odule dedicated to Sequences

.. topic:: summary

    A module dedicated to Sequences

    :Code status: mature
    :Documentation status: to be completed
    :Author: Thomas Cokelaer <Thomas.Cokelaer@sophia.inria.fr>
    :Revision: $Id$
    :Usage: >>> from openalea.sequence_analysis.sequences import *

"""
__version__ = "$Id$"

import os
import openalea.stat_tool.interface as interface
from openalea.sequence_analysis._sequence_analysis import (
    _Sequences, 
    _MarkovianSequences, 
    _RenewalData,
    _MarkovianSequences,
    _VariableOrderMarkovData,
    _SemiMarkovData,
    _NonHomogeneousMarkovData,
)


#import _sequence_analysis
from openalea.stat_tool import error



from .enums import index_parameter_type_map, ms_vomd_smd_nhmd
from openalea.stat_tool.enums import bool_type

__all__ = ['Sequences',
           '_Sequences',
           '_MarkovianSequences',
           'LumpabilityTest',
           'RemoveIndexParameter',
           'TransformPosition',
           'SaveMTG',
           'ComputeInitialRun',
           'IndexParameterType',
           'ComputeInitialRun',
           'Split']


# Extend dynamically class
interface.extend_class( _Sequences, interface.StatInterface)
interface.extend_class( _MarkovianSequences, interface.StatInterface)

# Add methods to _Vectors



def LumpabilityTest(obj, *args, **kargs):
    """.. todo:: documenation"""

    error.CheckArgumentsLength(args, 1, 1)
    error.CheckType([obj], [[_MarkovianSequences, _VariableOrderMarkovData,
                            _SemiMarkovData, _NonHomogeneousMarkovData]])


    symbol = args[0]
    Order = kargs.get("Order", 1)

    error.CheckType([symbol, Order], [list, int])

    ret = obj.lumpability_test(symbol, Order)

    if ret is False:
        raise TypeError("warning: false status returned by lumpability test")


def RemoveIndexParameter(obj):
    """.. todo:: documenation

    input can be sequence, markovian_sequences,
    nonhomogeneous_markov, variable_order_markov"""

    error.CheckType([obj], [[_Sequences, _MarkovianSequences,
                             _VariableOrderMarkovData, _SemiMarkovData,
                             _NonHomogeneousMarkovData]])
    if isinstance(obj, _Sequences):
        return obj.remove_index_parameter().markovian_sequences()
    else:
        return obj.remove_index_parameter()

def TransformPosition(obj, step=None):
    """.. todo:: documenation

    input is a sequence only"""
    error.CheckType([obj, step], [_Sequences, int])

    ret = obj.transform_position(step)
    if hasattr(obj, 'markovian_sequences'):
        return obj.markovian_sequences()
    else:
        return ret

def ComputeInitialRun(obj):
    """.. todo:: documenation

    input can be sequence, markovian_sequences, nonhomogeneous_markov,
    variable_order_markov
    """
    error.CheckType([obj], [[_MarkovianSequences, _VariableOrderMarkovData,
                             _SemiMarkovData, _NonHomogeneousMarkovData]])
    return obj.initial_run_computation()


def Split(obj, step):
    """.. todo:: documentaiton

    input markovian
    """
    error.CheckType([obj], [[_MarkovianSequences, _VariableOrderMarkovData,
                            _SemiMarkovData, _NonHomogeneousMarkovData]])
    error.CheckType([step], [int])

    return obj.split(step)

def SaveMTG(obj, Filename=None, Type=None):
    """
    Save sequence in MTG format.
    :param Type : list of "S" or "N" characters for Symbolic or Numeric
    """
    from openalea.stat_tool.enums import sub_variable_type

    error.CheckType([obj], [ms_vomd_smd_nhmd])
    error.CheckType([Filename, Type], [str, list])

    type = []
    for pstr in Type:
        type.append(sub_variable_type[pstr])

    obj.mtg_write(Filename, type)



def Sequences(obj, **kargs):
    """Construction of a set of sequences from multidimensional arrays
    of integers, from data generated by a renewal process or from an
    ASCII file.

    The data structure of type array(array(array(int))) should be
    constituted at the most internal level of arrays of constant size. If the
    optional argument IndexParameter is set at "Position" or "Time", the data
    structure of type array(array(array(int))) is constituted at the most
    internal level of arrays of size 1 + n (index parameter, n variables
    attached to the explicit index parameter). If the optional argument
    IndexParameter is set at "Position", only the index parameter of the
    last array of size 1 + n is considered and the first component of successive
    elementary arrays (representing the index parameter) should be
    ncreasing. If the optional argument IndexParameter is set at "Time", the
    first component of successive elementary arrays should be strictly
    increasing.


    :Parameters:

    * array1 (array(array(int))): input data for univariate sequences
    * arrayn (array(array(array(int)))): input data for multivariate sequences,
    * timev (renewal_data),
    * file_name (string).

    :Optional Parameters:

    * Identifiers (array(int)): explicit identifiers of sequences. This 
      optional argument can only be used if the first argument is of 
      type array(array(int / array(int))).
    * VertexIdentifiers (array(array(int))): explicit identifiers of vectors. 

    * IndexParameter (string): type of the explicit index parameter: "Position"
      or "Time" (the default: implicit discrete index parameter starting at 0). 
      This optional argument can only be used if the first argument is of type 
      array(array(int / array(int))).
    
    .. todo:: IndexParameterType

    :Returns:

    If the construction succeeds, an object of type sequences or 
    discrete_sequences is returned, otherwise no object is returned. The 
    returned object is of type discrete_sequences if all the variables are of 
    type STATE, if the possible values for each variable are consecutive from 0 
    and if the number of possible values for each variable is <= 15.

    :Examples:

    .. doctest::
        
        >>> # Single univariate sequence case (array1). 
        >>> seq1 = Sequences([1, 2, 3], Identifiers=[8])
        >>> seq1.nb_sequence
        1
        >>> seq1.nb_variable
        1
        >>> # General case arrayn
        >>> seq = Sequences([ 
        ...    [[1,2],[3,4]], 
        ...    [[21,22],[23,24]], 
        ...    [[31,32],[33,34], [35,36] ]], 
        ...    Identifiers = [1,8,12],
        ...    VertexIdentifiers = [[1,2],[3,4],[5,6,7]])
        >>> seq.nb_sequence
        3
        >>> seq.nb_variable
        2
        >>> seq.max_length
        3

    .. doctest::
        :options: +SKIP
        
            >>> Sequences(timev)
            >>> Sequences(file_name)
    
    .. seealso::

       :class:`~openalea.stat_tool.output.Save`,
       :func:`~openalea.sequence_analysis.data_transform.AddAbsorbingRun`,
       :func:`~openalea.stat_tool.cluster.Cluster`,
       :func:`~openalea.sequence_analysis.data_transform.Cumulate`,
       :func:`~openalea.sequence_analysis.data_transform.Difference`,
       :func:`~openalea.sequence_analysis.data_transform.IndexParameterExtract`,
       :func:`~openalea.sequence_analysis.data_transform.LengthSelect`,
       :func:`~openalea.stat_tool.data_transform.Merge`,
       :func:`~openalea.stat_tool.data_transform.MergeVariable`,
       :func:`~openalea.sequence_analysis.data_transform.MovingAverage`,
       :func:`~openalea.sequence_analysis.data_transform.RecurrenceTimeSequences`,
       :func:`~openalea.sequence_analysis.data_transform.RemoveRun`,
       :func:`~openalea.sequence_analysis.data_transform.Reverse`,
       :func:`~openalea.sequence_analysis.data_transform.SegmentationExtract`,
       :func:`~openalea.stat_tool.data_transform.SelectIndividual`,
       :func:`~openalea.stat_tool.data_transform.SelectVariable`,
       :func:`~openalea.stat_tool.data_transform.Shift`,
       :func:`~openalea.stat_tool.cluster.Transcode`,
       :func:`~openalea.stat_tool.data_transform.ValueSelect`,
       :func:`~openalea.sequence_analysis.data_transform.VariableScaling`.
       :func:`~openalea.stat_tool.data_transform.ExtractHistogram`,
       :func:`~openalea.sequence_analysis.data_transform.ExtractVectors`,
       :func:`~openalea.sequence_analysis.correlation.ComputeCorrelation`,
       :func:`~openalea.sequence_analysis.correlation.ComputePartialAutoCorrelation`,
       :func:`~openalea.sequence_analysis.data_transform.ComputeSelfTransition`,
       :func:`~openalea.sequence_analysis.compare.Compare`,
       :func:`~openalea.sequence_analysis.estimate.Estimate`,
       :func:`~openalea.sequence_analysis.data_transform.ComputeStateSequences`,
       :func:`~openalea.sequence_analysis.simulate.Simulate`.


    """
    import numpy

    sequence = None

    error.CheckType([obj], [[str, _RenewalData, list]])

    if isinstance(obj, str):
        filename = obj
        if os.path.isfile(filename):
            OldFormat = error.ParseKargs(kargs, "OldFormat", False, bool_type)
            sequence = _Sequences(filename, OldFormat)
        else:
            raise IOError("bad file name %s" % filename)
        if hasattr(sequence, 'markovian_sequences'):
            try:
                sequence = sequence.markovian_sequences()
            except Exception:
                pass
        try:
            sequence.nb_sequence
        except ValueError:
            raise ValueError("File read but issue while parsing. Returned sequence is not valid")
        return sequence
    elif isinstance(obj, _RenewalData):
        sequence = _Sequences(obj)
        if hasattr(sequence, 'markovian_sequences'):
            try:
                sequence = sequence.markovian_sequences()
            except Exception:
                pass



        return sequence

    # otherwise, we switch to a list constructor that requires a list of seqs
    # transform input into array of arrays of arrays 
    # case 1: general case where input = [[[1,2],[3,4]],[[1,2],[3,4], [5,6]]] nothing to do
    # case 2: univariate single sequence, input = [1,2,3,4,5,6] so it is [[[1],[2],[3],...]]
    # case 3: univariate sequences input = [[1,2],[3,4],[5,6,7]] (i.e, different vector sizes)
    # case 4: multivariate sequence input = [[1,2],[3,4],[5,6]]

    Verbose = error.ParseKargs(kargs, "Verbose", False)

    Univariate = error.ParseKargs(kargs, "Univariate", False)



    if type(obj)==list:
        first_sequence = obj[0]
        if (type(first_sequence) in [int, float]):
            obj = [[[x] for x in obj]]
            if Verbose:print('this is a single univariate sequence')
        elif type(first_sequence)==list:
            #either a single multivariate sequence ot general case of several sequences multivariates
            if type(first_sequence[0]) == list:
                if Verbose: print('this is the general case, nothing to do')
            elif type(first_sequence[0]) in [int, float]:
                lengths = numpy.array([len(x) for x in obj])
                if lengths.var()==0:
                    if Verbose:print('this is the ambiguous case')

                    if lengths[0]<5 and Univariate==False:
                        if Verbose:print('this is 1 single multivariate sequence')
                        obj = [obj]
                    else:
                        if Verbose:print('this is univariate sequences')
                        res = []
                        for x in obj:
                            res.append([[y] for y in x])
                        obj = res
                else:
                    if Verbose:print('this is univariate sequences')
                    res = []
                    for x in obj:
                        res.append([[y] for y in x])
                    obj = res
            else:
                print(SyntaxError('wrong syntax for input object'))


    # 0 for int, 1 for float. By default all variables are int
    #now, we loop over all sequences and sequences and if a variable 
    # is found to be float, then the type is float.
    # once a float is found, there is no need to carry on the current variable
    InputTypes = [0] * len(obj[0][0])
    nb_variables = len(obj[0][0])
    for seq in obj:
        for vec in seq:
            for index, var in enumerate(vec):
                assert type(var) in [int, float], "wrong types var=%s and its type is %s" % (var, type(var))
                if type(var)==float:
                    InputTypes[index]=1



    from openalea.sequence_analysis._sequence_analysis import TIME, POSITION, \
        IMPLICIT_TYPE
    #error.CheckArgumentsLength(args, 1, 1)


    IndexParameterType = error.ParseKargs(kargs, "IndexParameterType", "IMPLICIT_TYPE", index_parameter_type_map)
    IndexParameter = error.ParseKargs(kargs, "IndexParameter",  [])
    Identifiers = error.ParseKargs(kargs, "Identifiers", [])
    VertexIdentifiers = error.ParseKargs(kargs, "VertexIdentifiers", [])

 
    # build up a list of unique identifiers if none is provided
    lengths=[]
    for seq in obj:
        lengths.append(len(seq))
    # all values must be positive strictly
    if len(Identifiers)>0:
        assert len([x for x in Identifiers if x<=0]) == 0
    else:
        #create a standard identifiers list [0,1,2,....]
        for i, seq in enumerate(obj):
            Identifiers.append(i)

    # build up a list of unique vertex identifiers if none is provided
    if len(VertexIdentifiers)>0:
        assert len([x for x in VertexIdentifiers if x<=0]) == 0
    else:
        #create a standard identifiers list [0,1,2,....] for each sequences ?
        index = 0
        for i, seq in enumerate(obj):
            VertexIdentifiers.append([])
            for vec in  seq:
                VertexIdentifiers[i].append(index)
                index+=1

    # check unicity of vertex identifiers
    idents = []
    for seq in VertexIdentifiers:
        for ident in seq:
            idents.append(ident)
    assert len(set(idents)) == len(idents), "ERROR, VertexIdentifiers must be made of unique identifiers (for each vector)"

    # check unicity of identifiers
    idents = []
    for ident in Identifiers:
        idents.append(ident)
    assert len(set(idents)) == len(idents), "ERROR, Identifiers must be made of unique identifiers (for each sequence)"

    if len(IndexParameter)==0:
        index = 0
        for i, seq in enumerate(obj):
            IndexParameter.append([])
            for vec in seq:
                IndexParameter[i].append(index)
                index+=1
            if IndexParameterType==POSITION:
                IndexParameter[i].append(index)
                index+=1

    for i, seq in enumerate(obj):
        #print len(seq), len(IndexParameter)
        if IndexParameterType==POSITION:
            assert len(seq)==len(IndexParameter[i])-1, "ERROR, wrong IndexParameterLength. When ParameterType=POSITION, ParameterIndex length must be equla to the sequence length +1"
        else:
            assert len(seq)==len(IndexParameter[i]), "ERROR, wrong IndexParameterLength. ParameterIndex length must be equal to the sequence length."
   #todo check that indesparameter length is correct (length of vectors +1 if position)


    valid_param = [POSITION, TIME, IMPLICIT_TYPE]
    if IndexParameterType not in valid_param:
        raise ValueError("""IndexParameter can be only %s if first
            argument is a list""" % valid_param)
    sequence = _Sequences(obj, Identifiers, VertexIdentifiers, IndexParameter, InputTypes, IndexParameterType)

    if hasattr(sequence, 'markovian_sequences'):
        try:
            sequence = sequence.markovian_sequences()
        except Exception:
            pass

    return sequence









def IndexParameterType(obj):
    """

    input can be sequence, markovian_sequences,
    nonhomogeneous_markov, variable_order_markov

    .. doctest::
        :options: +SKIP
        
        >>> obj.index_parameter_type
        3
        >>> IndexParameterType(obj)
        openalea.sequence_analysis._sequence_analysis.IndexParameterType.POSITON
    """

    error.CheckType([obj], [[_Sequences, _MarkovianSequences,
                             _VariableOrderMarkovData, _SemiMarkovData,
                             _NonHomogeneousMarkovData]])
    from openalea.sequence_analysis.enums import index_parameter_type_map

    type = obj.index_parameter_type
    for key, value in index_parameter_type_map.items():
        if value==type:
            return key
