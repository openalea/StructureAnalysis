#!/usr/bin/env python
#-*- coding: utf-8 -*-
"""Cluster functions and classes

.. topic:: cluster.py summary

    A module dedicated to Clustering

    :Code status: mature
    :Documentation status: to be completed
    :Author: Samuel Dufour-Kowalski <samuel.dufour@sophia.inria.fr>
        Thomas Cokelaer <Thomas.Cokelaer@sophia.inria.fr>

    :Revision: $Id: cluster.py 15183 2013-11-06 10:35:50Z jbdurand $

"""
__version__ = "$Id: cluster.py 15183 2013-11-06 10:35:50Z jbdurand $"

import error
import interface

from openalea.stat_tool._stat_tool import _DistanceMatrix
from openalea.stat_tool._stat_tool import _Cluster
from openalea.stat_tool._stat_tool import _Dendrogram

from enums import format_type
from enums import criterion_type
from enums import algorithm_type
from enums import cluster_type

__all__ = [
     "_DistanceMatrix",
     "_Cluster",
     "_Dendrogram",
     "Cluster",
     "Transcode",
     "Clustering",
     "ToDistanceMatrix", ]


# Extend classes dynamically
interface.extend_class(_DistanceMatrix, interface.StatInterface)
interface.extend_class(_Cluster, interface.StatInterface)
interface.extend_class(_Dendrogram, interface.StatInterface)


def Cluster(obj, utype, *args, **kargs):
    """Clustering of values.

    In the case of the clustering of values of a frequency distribution on the
    basis of an information measure criterion (argument `Information`), both the
    information measure ratio and the selected optimal step are given in the
    shell window.

    The clustering mode `Step` (and its variant `Information`) is naturally
    adapted to numeric variables while the clustering mode `Limit` applies to
    both symbolic (nominal) and numeric variables. In the case of a symbolic
    variable, the function `Cluster` with the mode `Limit` can be seen as a
    dedicated interface of the more general function `Transcode`.

    :Parameters:

      * `histo` (`_FrequencyDistribution`, `_DiscreteMixtureData`, `_ConvolutionData`, `_CompoundData`),
      * `step` (int) - step for the clustering of values
      * `information_ratio` (float) - proportion of the information measure of \
        the original sample for determining the clustering step,
      * `limits` (list(int)) - first values corresponding to the new classes \
        classes 1, ..., nb_class - 1. By convention, the first value corresponding \
        to the first class is 0,
      * `vec1` (`_Vector`) - values,
      * `vecn` (`_Vectors`) - vectors,
      * `variable` (int) - variable index,
      * `seq1` (`_Sequences`) - univariate sequences,
      * `seqn` (`_Sequences`) - multivariate sequences,
      * `discrete_seq1` (`_DiscreteSequences`, `_Markov`, `_SemiMarkovData`) -
        discrete univariate sequences,
      * `discrete_seqn` (`_DiscreteSequences`, `_Markov`, `_SemiMarkovData`) -
        discrete multivariate sequences.

    :Keywords:

      * `AddVariable` (bool) : addition (instead of simple replacement) of the variable
        corresponding to the clustering of values (default value: False).
        This optional argument can only be used if the first argument is of
        type `_DiscreteSequences`, `_Markov` or `_SemiMarkovData`. The addition
        of the clustered variable is particularly useful if one wants to evaluate
        a lumpability hypothesis.

    :Returns:

      * If `step` > 0, or if 0 <  `information_ratio` <  1, or if 0 < limits[1]
        < limits[2] < ... < limits[nb_class - 1] < (maximum possible value of histo),
        an object of type _FrequencyDistribution is returned.
      * If variable is a valid index of a variable and if `step` > 0, or
        if 0 < limits[1] < limits[2] < ... < limits[nb_class - 1] < (maximum possible
        value taken by the selected variable of `vec1` or `vecn`), an object of type
        `_Vectors` is returned.
      * If variable is a valid index of a variable of type STATE and if `step` > 0, or \
        if 0 < limits[1] < limits[2] < ... < limits[nb_class - 1] < (maximum
        possible value taken by the selected variable of `seq1`, `seqn`, `discrete_seq1`
        or `discrete_seqn`), an object of type `_Sequences` or `_DiscreteSequences`
        is returned.
      * In the case of a first argument of type `_Sequences`, an object of type
        `_DiscreteSequences` is returned if all the variables are of type STATE,
        if the possible values taken by each variable are consecutive from 0 and
        if the number of possible values for each variable is < 15.

    :Examples:

    .. doctest::
        :options: +SKIP

        >>> Cluster(histo, "Step", step)
        >>> Cluster(histo, "Information", information_ratio)
        >>> Cluster(histo, "Limit", limits)
        >>> Cluster(vec1, "Step", step)
        >>> Cluster(vecn, "Step", variable, step)
        >>> Cluster(vec1, "Limit", limits)
        >>> Cluster(vecn, "Limit", variable, limits)
        >>> Cluster(seq1, "Step", step)
        >>> Cluster(seqn, "Step", variable, step)
        >>> Cluster(discrete_seq1, "Step", step, AddVariable=True)
        >>> Cluster(discrete_seqn, "Step", variable, step, AddVariable=True)
        >>> Cluster(seq1, "Limit", limits)
        >>> Cluster(seqn, "Limit", variable, limits)
        >>> Cluster(discrete_seq1, "Limit", limits, AddVariable=True)
        >>> Cluster(discrete_seqn, "Limit", variable, limits, AddVariable=True)

    .. seealso::
        :func:`~openalea.stat_tool.data_transform.Merge`,
        :func:`~openalea.stat_tool.data_transform.Shift`,
        :func:`~openalea.stat_tool.data_transform.ValueSelect`,
        :func:`~openalea.stat_tool.data_transform.MergeVariable`,
        :func:`~openalea.stat_tool.data_transform.SelectIndividual`,
        :func:`~openalea.stat_tool.data_transform.SelectVariable`,
        :func:`~openalea.stat_tool.cluster.Transcode`,
        :func:`~openalea.stat_tool.data_transform.AddAbsorbingRun`,
        :func:`~openalea.stat_tool.data_transform.Cumulate`,
        :func:`~openalea.stat_tool.data_transform.Difference`,
        :func:`~openalea.stat_tool.data_transform.IndexExtract`,
        :func:`~openalea.stat_tool.data_transform.LengthSelect`,
        :func:`~vplants.sequence_analysis.data_transform.MovingAverage`,
        :func:`~openalea.stat_tool.data_transform.RecurrenceTimeSequences`,
        :func:`~openalea.stat_tool.data_transform.Removerun`,
        :func:`~openalea.stat_tool.data_transform.Reverse`,
        :func:`~openalea.stat_tool.data_transform.SegmentationExtract`,
        :func:`~openalea.stat_tool.data_transform.VariableScaling`.
    """

    # fixme: what about the Mode in the Step case ?
    # check markovian_sequences call in Sequences
    AddVariable = error.ParseKargs(kargs, "AddVariable", False,
                                   possible=[False, True])


    error.CheckArgumentsLength(args, 1, 2)


    # search for the function name
    if hasattr(obj, cluster_type[utype]):
        func = getattr(obj, cluster_type[utype])
    else:
        raise KeyError("""Possible action are : 'Step', 'Information' or
        'Limit'. Information cannot be used with Vectors objects""")

    # check if nb_variable is available (vectors, sequences)
    if hasattr(obj, 'nb_variable'):
        nb_variable = obj.nb_variable
    else:
        nb_variable = 1

    #check types
    if nb_variable == 1:
        if len(args) == 1:
            if utype == "Step":
                error.CheckType([args[0]], [int])
            if utype == "Limit":
                error.CheckType([args[0]], [list])
            if utype == "Information":
                error.CheckType([args[0]], [[int, float]])
            try:
                ret = func(args[0]) # histogram case
            except:
                try:
                    ret = func(1, args[0]) # vector case
                except:
                    try:
                        ret = func(1, args[0], AddVariable) # sequences case
                    except:
                        pass
        else:
            raise ValueError("""Extra arguments provided 
            (to specify variable value ?). Consider removing it. 
            Be aware that nb_variable equals 1""")

    else:
        if len(args) == 2:
            if utype == "Step":
                error.CheckType([args[0]], [int])
                error.CheckType([args[1]], [[int, float]])
            if utype == "Limit":
                error.CheckType([args[0]], [int])
                error.CheckType([args[1]], [list])
            try:
                ret = func(*args)
            except:
                ret = func(args[0], args[1], AddVariable)
        else:
            raise ValueError("""Extra arguments provided 
            (to specify variable value ?).
            Consider removing it. Be aware that nb_variable equals 1""")




    if hasattr(ret, 'markovian_sequences'):
        ret = ret.markovian_sequences()

    return ret


def Transcode(obj, *args, **kargs):
    """
    Transcoding of values.

    The function `Cluster` with the mode "Limit" can be seen as a dedicated interface of
    the more general function Transcode.

    :Parameters:

      * `histo` (_FrequencyDistribution, _MixtureData, _ConvolutionData, _CompoundData),
      * `new_values` (array(int)) - new values replacing the old ones min, min + 1, ..., max.
      * `vec1` (_Vectors) - values,
      * `vecn` (_Vectors) - vectors,
      * `variable` (int) - variable index,
      * `seq1` (_Sequences) - univariate sequences,
      * `seqn` (_Ssequences) - multivariate sequences,
      * `discrete_seq1` (_DiscreteSequences, _MarkovData, _SemiMarkovData) - discrete univariate sequences,
      * `discrete_seqn` (_DiscreteSequences, _MarkovData, _SemiMarkovData) - discrete multivariate sequences.

    :Keywords:

      * AddVariable (bool): addition (instead of simple replacement) of the variable
        to which the transcoding is applied (default value: False). This optional argument
        can only be used if the first argument is of type (_DiscreteSequences, _MarkovData,
        _SemiMarkovData).

    :Returns:

        If the new values are in same number as the old values and are consecutive from 0,
        an object of type _FrequencyDistribution is returned (respectively _Vectors, _Sequences or
        _DiscreteSequences). In the case of a first argument of type _Sequences, the
        returned object is of type _DiscreteSequences if all the variables are of type STATE,
        if the possible values for each variable are consecutive from 0 and if the number of
        possible values for each variable is < 15.

    :Examples:

    .. doctest::
        :options: +SKIP

        >>> Transcode(histo, new_values)
        >>> Transcode(vec1, new_values)
        >>> Transcode(vecn, variable, new_values)
        >>> Transcode(seq1, new_values)
        >>> Transcode(seqn, variable, new_values)
        >>> Transcode(discrete_seq1, new_values, AddVariable=True)
        >>> Transcode(discrete_seqn, variable, new_values, AddVariable=True)

    .. seealso::

        :func:`~openalea.stat_tool.cluster.Clustering`,
        :func:`~openalea.stat_tool.data_transform.Merge`,
        :func:`~openalea.stat_tool.data_transform.Shift`,
        :func:`~openalea.stat_tool.data_transform.ValueSelect`,
        :func:`~openalea.stat_tool.data_transform.MergeVariable`,
        :func:`~openalea.stat_tool.data_transform.SelectIndividual`,
        :func:`~openalea.stat_tool.data_transform.SelectVariable`,
        :func:`~openalea.stat_tool.cumulate.Cumulate`,
        :func:`~openalea.stat_tool.data_transform.AddAbsorbingRun`,
        :func:`~openalea.stat_tool.cumulate.Cumulate`,
        :func:`~openalea.stat_tool.data_transform.Difference`,
        :func:`~openalea.stat_tool.data_transform.IndexExtract`,
        :func:`~openalea.stat_tool.data_transform.LengthSelect`,
        :func:`~openalea.stat_tool.data_transform.MovingAverage`,
        :func:`~openalea.stat_tool.data_transform.RecurrenceTimeSequences`,
        :func:`~openalea.stat_tool.data_transform.Removerun`,
        :func:`~openalea.stat_tool.data_transform.Reverse`,
        :func:`~openalea.stat_tool.data_transform.SegmentationExtract`,
        :func:`~openalea.stat_tool.data_transform.VariableScaling`.
    """
    AddVariable = error.ParseKargs(kargs, "AddVariable", False,
                                   possible=[False, True])

    myerror = "Arguments do not seem to be correct"

    if hasattr(obj, 'nb_variable'):# case sequence, vectors
        nb_variable = obj.nb_variable
        if len(args)==1 and nb_variable == 1:


            try:
                ret = obj.transcode(1, args[0], AddVariable)
            except:
                ret = obj.transcode(1, args[0])

        elif len(args)==2 and nb_variable!=1:
            try:
                ret = obj.transcode(args[0], args[1], AddVariable)
            except:
                ret = obj.transcode(args[0], args[1])

        else:
            raise ValueError(myerror)




    else:# case histogram and co
        nb_variable = None
        new_values = args[0]
        if len(args)>1:
            raise ValueError(myerror)
        ret  = obj.transcode(new_values)


    if ret is None:
        raise Exception("transcode function did not return anything...")
    else:
        #Sequence case to be converted to Markovian_sequences
        # todo: test and checks
        if hasattr(ret, 'markovian_sequences'):
            func = getattr(ret, 'markovian_sequences')
            ret = func()
        return ret




def Clustering(matrix, utype, *args, **kargs):
    """
    Application of clustering methods (either partitioning methods or hierarchical methods)
    to dissimilarity matrices between patterns.

    In the case where the composition of clusters is a priori fixed,
    the function Clustering simply performs an evaluation of the a priori fixed
    partition.

    :Parameters:
      * `dissimilarity_matrix` (distance_matrix) - dissimilarity matrix between patterns,
      * `nb_cluster` (int) - number of clusters,
      * `clusters` (list(list(int))) - cluster composition.

    :Keywords:
      * `Prototypes` (list(int)): cluster prototypes.
      * `Algorithm` (string): "Agglomerative", "Divisive" or "Ordering"
      * `Criterion` (string): "FarthestNeighbor" or "Averaging"
      * `Filename` (string): filename
      * `Format` (string) : "ASCII" or "SpreadSheet"

    :Returns:
        If the second mandatory argument is "Partitioning" and
        if 2 < nb_cluster < (number of patterns), an object of type clusters
        is returned

    :Examples:

    .. doctest::
        :options: +SKIP

        >>> Clustering(dissimilarity_matrix, "Partition", nb_cluster, Prototypes=[1, 3, 12])
        >>> Clustering(dissimilarity_matrix, "Partition", clusters)
        >>> Clustering(dissimilarity_matrix, "Hierarchy", Algorithm="Agglomerative")
        >>> Clustering(dissimilarity_matrix, "Hierarchy", Algorithm="Divisive")

    .. seealso::
        :func:`~openalea.stat_tool.data_transform.SelectIndividual`,
        `Symmetrize`,
        :func:`~openalea.stat_tool.comparison.Compare`,
        :func:`~openalea.stat_tool.cluster.ToDistanceMatrix`.

    .. note:: if type=Partition, Algorthim must be 1 (divisive) or 2 (ordering).

    .. note:: if type!=Divisive criterion must be provided
    """
    #TODO: check this case :
    #Clustering(dissimilarity_matrix, "Partition", clusters)

    error.CheckType([matrix], [_DistanceMatrix])

    Algorithm = error.ParseKargs(kargs, "Algorithm",
                                 default="Divisive",
                                 possible=algorithm_type)
    # Switch for each type of clustering
    # first the partition case
    if utype == "Partition":
        error.CheckArgumentsLength(args, 1, 1)
        error.CheckKargs(kargs,
                         ["Algorithm", "Prototypes", "Initialization"])
        Initialization = error.ParseKargs(kargs, "Initialization", 1,
                                          possible=[1, 2])

        if Algorithm == algorithm_type["Agglomerative"]:
            raise ValueError("""If partition is on, Algorithm cannot
                    be agglomerative""")


        if(isinstance(args[0], int)): #int case
            # if Prototypes is empty, the wrapping will send an
            # int * = 0 to the prototyping function, as expected
            Prototypes = kargs.get("Prototypes", [])
            nb_cluster = args[0]
            return matrix.partitioning_prototype(nb_cluster, Prototypes,
                                                 Initialization, Algorithm)
        elif isinstance(args[0], list): # array case
            #todo:: array of what kind of object?
            #need a test
            return matrix.partitioning_clusters(args[0])
        else:
            raise TypeError("""
            With Partition as second argument, the third one must be either
            an int or an array.""")

    elif utype == "Hierarchy":
        error.CheckKargs(kargs,
                        ["Algorithm",  "FileName", "Criterion", "Format"])

        Algorithm = error.ParseKargs(kargs, "Algorithm",
                                             default="Agglomerative",
                                             possible=algorithm_type)

        Criterion = error.ParseKargs(kargs, "Criterion", "Averaging",
                                     possible=criterion_type)

        # fixme: is it correct to set "" to the filename by defautl ?
        # if set to None, the prototype does not match
        filename = kargs.get("Filename", None)
        format = error.ParseKargs(kargs, "Format", "ASCII",
                                  possible=format_type)
        #check options
        if Algorithm != algorithm_type["Agglomerative"] and \
            kargs.get("Criterion"):

            raise ValueError("""
                In the Hierarchy case, if Algorithm is different from
                AGGLOMERATIVE, then Criterion cannot be used.""")
        return matrix.hierarchical_clustering(Algorithm, Criterion,
                                            filename, format)


    else:
        raise KeyError("Second argument must be 'Partitioning' or 'Hierarchy'")




def ToDistanceMatrix(distance_matrix):
    """
    Cast and object of type CLUSTER into an object of type DISTANCE_MATRIX.


    :Parameters:
      * distance_matrix

    :Returns:
        An object of type distance_matrix is returned.

    :Examples:

    .. doctest::
        :options: +SKIP

        >>> ToDistanceMatrix(distance_matrix)

    .. seealso::
        :func:`~openalea.stat_tool.cluster.Clustering`,

    """
    error.CheckType([distance_matrix], [[_Cluster, _DistanceMatrix]])

    try:
        return _DistanceMatrix(distance_matrix)
    except:
        raise TypeError("Input arguments must be of type Cluster")

