#!/usr/bin/env python
#-*- coding: utf-8 -*-
"""Vectors function and class

.. topic:: vectors.py summary

    A module dedicated to Vectors

    :Code status: mature
    :Documentation status: to be completed
    :Author: Thomas Cokelaer <Thomas.Cokelaer@sophia.inria.fr>

    :Revision: $Id: vectors.py 9404 2010-08-10 14:45:08Z cokelaer $
    
"""
__version__ = "$Id: vectors.py 9404 2010-08-10 14:45:08Z cokelaer $"

import interface
import error

from openalea.stat_tool.multivariate_mixture import _MultivariateMixture

from openalea.stat_tool._stat_tool import _Vectors
from openalea.stat_tool._stat_tool import _VectorDistance
from openalea.stat_tool._stat_tool import I_DEFAULT
from openalea.stat_tool.enums import format_type
from openalea.stat_tool.enums import variable_type
from openalea.stat_tool.enums import variance_type
from openalea.stat_tool.enums import distance_type

__all__ = ['Vectors',
           '_Vectors',
           'VectorDistance',
           '_VectorDistance',
           'ContingencyTable',
           'VarianceAnalysis',
           'ComputeRankCorrelation']


############### VectorDistance #################################################




# Add methods to _Vectors

def _Vectors_mixture_estimation(self, model,
                                nb_iteration=I_DEFAULT,
                                force_param=None):
    """Estimate a mixture from _Vectors given initial model or number of
    components, the maximal number of iterations and a flag for using parametric
    observation distributions or not, within a given family
    """
    if force_param is None:
        force_param = []

    error.CheckType([nb_iteration, force_param], [int, list])

    # model is a MultivariateMixture class
    error.CheckType([model], [[int, _MultivariateMixture]])
    if type(model) == int:
        return _Vectors.mixture_estimation_nb_component(self, model,
                                            nb_iteration, force_param)
    else:
        return _Vectors.mixture_estimation_model(self, model,
                                            nb_iteration, force_param)

_Vectors.mixture_estimation = _Vectors_mixture_estimation


def Vectors(*args, **kargs):
    """
    Construction of a set of vectors from a multidimensional array,
    from a set of sequences or from an ASCII file.

    The data structure of type list(list(int)) should be constituted at the
    most internal level of arrays of constant size.

    :Parameters:
      - `list` (list(list(int))) :
      - `seq` (sequences, discrete_sequences, markov_data, semi-markov_data)
      - `file_name` (string) :

    :Keywords:
      - Identifiers (list(int)): explicit identifiers of vectors.
        This optional argument can only be used if the first mandatory argument is of
        type list(list(int)).
      - IndexVariable (bool): taking into account of the implicit index parameter as
        a supplementary variable (default value: False). This optional argument can
        only be used if the first mandatory argument is of type `sequences`,
        `discrete_sequences`, `markov_data` or `semi-markov_data`.

    :Returns:
       If the construction succeeds, an object of type vectors is returned,
       otherwise no object is returned.

    :Examples:

    .. doctest::
        :options: +SKIP

        >>> Vectors(list, Identifiers=[1, 8, 12])
        >>> Vectors(seq, IndexVariable=True)
        >>> Vectors(file_name)

    .. seealso::
        :func:`~openalea.stat_tool.output.Save`,
        :func:`~openalea.stat_tool.data_transform.ExtractHistogram`,
        :func:`~openalea.stat_tool.cluster.Cluster`,
        :func:`~openalea.stat_tool.data_transform.Merge`,
        :func:`~openalea.stat_tool.data_transform.MergeVariable`,
        :func:`~openalea.stat_tool.data_transform.SelectIndividual`,
        :func:`~openalea.stat_tool.data_transform.SelectVariable`,
        :func:`~openalea.stat_tool.data_transform.Shift`,
        :func:`~openalea.stat_tool.cluster.Transcode`,
        :func:`~openalea.stat_tool.data_transform.ValueSelect`,
        :func:`~openalea.stat_tool.comparison.Compare`,
        :func:`~openalea.stat_tool.comparison.ComputeRankCorrelation`,
        :func:`~openalea.stat_tool.comparison.ContingencyTable`,
        :func:`~openalea.stat_tool.comparison.Regression`,
        :func:`~openalea.stat_tool.comparison.VarianceAnalysis`
    """
    error.CheckArgumentsLength(args, 1, 1)
    error.CheckKargs(kargs, possible_kargs = ["Identifiers", 
                                              "IndexVariable"])

    obj = args[0]
    ret = None

    import openalea.core.path

    if isinstance(obj, str):
        # constructor from a filename
        ret = _Vectors(args[0])
    elif isinstance(obj, openalea.core.path.path):
        # constructor from a path
        ret = _Vectors(str(args[0]))
    elif isinstance(obj, list):
        # Normal usage is Vectors([ [1,2,3],  [1,2,3], [4,5,6]])
        # If only one variable is requited, then Normal usage is
        # Vectors([ [1,2,3] ]). Yet, to simplify usage, if there is only
        # one variable, the followin if allows us to use Vectors([1,2,3])
        if type(obj[0])!=list:
            obj = [obj]



        # 0 for int, 1 for float. By default all variables are int
        #now, we loop over all sequences and sequences and if a variable 
        # is found to be float, then the type is float.
        # once a float is found, there is no need to carry on the current variable
        InputTypes = [0] * len(obj[0])
        nb_variables = len(obj[0])
        for vec in obj:
            for index, var in enumerate(vec):
                assert type(var) in [int, float], "wrong types var=%s and its type is %s" % (var, type(var))
                if type(var)==float:
                    InputTypes[index]=1


        # from a list and an optional argument

        # first, get the Identifiers and check its type
        identifiers = error.ParseKargs(kargs, "Identifiers")
        if identifiers:
            error.CheckType([identifiers], [[list]], variable_pos=[2])

            if len(identifiers) != len(obj):
                raise ValueError("""Identifiers must be a list,
                which size equals vectors's length""")
            #iif InputTypes: 
            ret = _Vectors(obj, identifiers, InputTypes)
            #else:
            #    ret = _Vectors(obj, identifiers)
        else:

            #create a standard identifiers list [0,1,2,....] for each sequences ?
            identifiers = []
            for i, vec in enumerate(obj):
                identifiers.append(i+1)
    
            print identifiers
            #if InputTypes:
            ret = _Vectors(obj, identifiers, InputTypes)
            #else:
            #    ret = _Vectors(obj, [])
    else:
        # from a sequence
        index_variable = error.ParseKargs(kargs, "IndexVariable", False,
                                          [True, False])
        error.CheckType([index_variable], [bool], variable_pos=[2])
        ret = obj.build_vectors(index_variable)


    return ret


interface.extend_class( _Vectors, interface.StatInterface)



def VectorDistance(*args, **kargs):
    """
    Construction of an object of type vector_distance from types (and eventually weights)
    of variables or from an ASCII file.

    The type _VectorDistance implements standardization procedures. The objective of
    standardization is to avoid the dependence on the variable type
    (chosen among symbolic, ordinal, numeric and circular) and, for numeric variables,
    on the choice of the measurement units by converting the original variables to
    unitless variables.

    :Parameters:
      - `type1`, `type2`, ... (string):
        variable types ("NUMERIC" ("N"), "ORDINAL" ("O") or "SYMBOLIC" ("S")),
      - `weight1`, `weight2`, ... (float): weights of variables,
      - `file_name` (string).

    :Keywords:

      * Distance (string): distance type: "ABSOLUTE_VALUE" (default) or "QUADRATIC".
        This optional argument is only relevant in the multivariate case.

    :Returns:
        If the construction succeeds, an object of type vector_distance is returned.

    :Examples:

    .. doctest::
        :options: +SKIP

        >>> VectorDistance(type1, type2,..., Distance="QUADRATIC")
        >>> VectorDistance(weight1, type1, weight2, type2,..., Distance="QUADRATIC")
        >>> VectorDistance(file_name)

    .. seealso::
        :func:`~openalea.stat_tool.comparison.Compare`
    """

    error_arguments = ["",
                       """If first argument is a number, following
    argument must be in ["N", "O", "S"]. Check documentation by typing
    VectorDistance? .""",
    ""]

    distance = error.ParseKargs(kargs, "Distance", "ABSOLUTE_VALUE",
                                distance_type)


    # Case VectorDistance("O", "N", "S")
    if args[0] in variable_type.keys():
        # check that all following arguments (if any) are correct
        types = []
        for arg, index in zip(args, range(0, len(args))):
            # check that the arguments are correct
            if arg not in variable_type.keys():
                raise ValueError(error_arguments[1])
            else:
                types.append(variable_type[arg])
        # assign a uniform weights since none were provided
        weights = [1./len(types) for _elem in types]

        return _VectorDistance(types, weights, distance)
    # same as above but with weights VectorDistance(0.5, "N", 0.5, "S")
    if isinstance(args[0], int) or isinstance(args[0], float):
        types = list(args[1:len(args):2])
        weights = list(args[0:len(args):2])
        assert len(types)==len(weights)

        # check that types are strings
        error.CheckType(types, [str]*len(types))
        # check that weights are integer or floats
        error.CheckType(weights, [[int, float]]*len(weights))

        # convert to vector_distance_type
        for arg, index in zip(types, range(0, len(types))):
            types[index] = variable_type[types[index]]

        return _VectorDistance(types, weights, distance)
    # filename case
    elif isinstance(args[0], str) and len(args)==1 and \
            args[0] not in variable_type.keys():
        return _VectorDistance(args[0])


# Extend dynamically class
interface.extend_class( _VectorDistance, interface.StatInterface)


def VarianceAnalysis(*args, **kargs):
    """
    One-way variance analysis.

    :Examples:

    .. doctest::
        :options: +SKIP

        >>> VarianceAnalysis(vec, class_variable, response_variable,
        ... type, FileName="result", Format="SpreadSheet")

    :Parameters:

      * vec (_Vectors),
      * class_variable (int): index of the class or group variable,
      * response_variable (int): index of the response variable,
      * type (string): type of the response variable ("NUMERIC" ("N") or
        "ORDINAL" ("O")).

    :Keywords:

      * FileName (string): name of the result file,
      * Format (string): format of the result file: "ASCII" (default format)
        or "SpreadSheet". This optional argument can only be used in conjunction with the optional argument FileName.

    :Returns:

        The variance analysis result as a string

    """
    error.CheckArgumentsLength(args, 4, 4)
    error.CheckKargs(kargs, possible_kargs = ["FileName", "Format"])

    #kargs
    filename = error.ParseKargs(kargs, "FileName", default="result")
    format = error.ParseKargs(kargs, "Format", default="ASCII",
                              possible=format_type)

    #args
    vec = args[0]
    class_variable = args[1]
    response_variable = args[2]
    utype = args[3]
    error.CheckType([vec, class_variable, response_variable, utype],
                    [_Vectors, int, int, str])

    try:
        utype = variance_type[args[3]]
    except KeyError:
        raise KeyError("Possible type are : " + str(variance_type.keys()))


    return vec.variance_analysis(class_variable, response_variable, utype,
                                 filename, format)



def ContingencyTable(*args, **kargs):
    """
    Computation of a contingency table.

    :Parameters:

      * vec (_Vectors),
      * variable1, variable2 (int): variable indices,

    :Keywords:

      * FileName (string): name of the result file,
      * Format (string): format of the result file: "ASCII" (default format) or "SpreadSheet".
        This optional argument can only be used in conjunction with the optional argument FileName.

    :Returns:

        The contingency table result as a string

    :Examples:

    .. doctest::
        :options: +SKIP

        >>> ContingencyTable(vec, variable1, variable2, FileName="result", Format="SpreadSheet")


    """
    error.CheckArgumentsLength(args, 3, 3)
    error.CheckKargs(kargs, possible_kargs = ["FileName", "Format"])

    #kargs
    filename = error.ParseKargs(kargs, "FileName", default="result")
    format = error.ParseKargs(kargs, "Format", default="ASCII",
                              possible=format_type)

    #args
    vec = args[0]
    variable1 = args[1]
    variable2 = args[2]
    error.CheckType([vec, variable1, variable2], [_Vectors, int, int])

    return vec.contingency_table(variable1, variable2, filename, format)


def ComputeRankCorrelation(*args, **kargs):
    """ComputeRankCorrelation

    Computation of the rank correlation matrix.

    :Usage:

    >>> vec = Vectors([1,2,3,4,5,4,3,2,1])
    >>> ComputeRankCorrelation(vec, Type="Spearman", FileName='')

    :Arguments:

    * vec (vectors).

    :Optional Arguments:

    * Type (string): type of rank correlation coefficient:
      "Spearman" (the default) or "Kendall".

    :Returned Object:

    No object returned.
    """

    func_map = {
            "Spearman": 0,
            "Kendall": 1
            }

    error.CheckArgumentsLength(args, 1, 1)
    error.CheckKargs(kargs, possible_kargs = ["Type", "FileName"])

    #kargs
    utype = error.ParseKargs(kargs, "Type", default="Spearman",
                             possible=func_map)
    filename = error.ParseKargs(kargs, "FileName", default=None)

    #args
    vec = args[0]

    error.CheckType([vec], [_Vectors])

    _a = vec.rank_correlation_computation(utype, filename)




interface.extend_class( _Vectors, interface.StatInterface)

