"""Cluster functions and classes"""
__revision__ = "$Id$"

import _stat_tool
import interface
from _stat_tool import _DistanceMatrix
from _stat_tool import _Clusters
from _stat_tool import _Dendrogram

__all__ = [
     "_DistanceMatrix",
     "_Clusters",
     "_Dendrogram",
     "Cluster",
     "Transcode",
     "Clustering",
     "ToDistanceMatrix",]


# Extend classes dynamically
interface.extend_class( _DistanceMatrix, interface.StatInterface)
interface.extend_class( _Clusters, interface.StatInterface)
interface.extend_class( _Dendrogram, interface.StatInterface)


def Cluster(obj, type, *args, **kargs):
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

      * `histo` (`_Histogram`, `_MixtureData`, `_ConvolutionData`, `_CompoundData`),
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
        an object of type _Histogram is returned.
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
        :func:`~openalea.stat_tool.data_transform.MovingAverage`,
        :func:`~openalea.stat_tool.data_transform.RecurrenceTimeSequences`,
        :func:`~openalea.stat_tool.data_transform.Removerun`, 
        :func:`~openalea.stat_tool.data_transform.Reverse`,
        :func:`~openalea.stat_tool.data_transform.SegmentationExtract`,
        :func:`~openalea.stat_tool.data_transform.VariableScaling`. 
    """
    
    AddVariable = kargs.get("AddVariable", False)

    type_map = { "Step": "cluster_step", 
                 "Limit": "cluster_limit",
                 "Information" : "cluster_information" }

    # Get param : step, ratio or limits
    try:
        param = args[0]
    except IndexError:
        raise TypeError("Cluster expect an argument for '%s'"%(type))

    # Map type with a function name
    try:
        func_name = type_map[type]
        func = getattr(obj, func_name)
    except KeyError:
        raise KeyError("Possible action are : 'Step' or 'Information' or 'Limit'")

    except AttributeError:
        raise AttributeError("Object doesn't support %s."%(func_name))


    # optinal arg : variable
    try:
        nb_variable = obj.get_nb_variable()
        if(nb_variable == 1):
            variable = 1

        return func(variable, param)

    except AttributeError: #no get_nb_variable
        return func(param)
    
    

def Transcode(obj, param1, param2=None, AddVariable=False):
    """    
    Transcoding of values.
    
    The function `Cluster` with the mode "Limit" can be seen as a dedicated interface of 
    the more general function Transcode.
    
    :Parameters:

      * `histo` (_Histogram, _MixtureData, _ConvolutionData, _CompoundData),
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
        an object of type _Histogram is returned (respectively _Vectors, _Sequences or 
        _DiscreteSequences). In the case of a first argument of type _Sequences, the 
        returned object is of type _DiscreteSequences if all the variables are of type STATE, 
        if the possible values for each variable are consecutive from 0 and if the number of 
        possible values for each variable is < 15. 

    :Examples:
    
      >>> Transcode(histo, new_values)
      >>> Transcode(vec1, new_values)
      >>> Transcode(vecn, variable, new_values)
      >>> Transcode(seq1, new_values)
      >>> Transcode(seqn, variable, new_values)
      >>> Transcode(discrete_seq1, new_values, AddVariable=True)
      >>> Transcode(discrete_seqn, variable, new_values, AddVariable=True)
    
    .. seealso::
    
        `Cluster`, `Merge`, `Shift`, `ValueSelect`, `MergeVariable`, `SelectIndividual`, 
        `SelectVariable`, `AddAbsorbingRun`, `Cumulate`, `Difference`, `IndexExtract`, 
        `LengthSelect`, `MovingAverage`, `RecurrenceTimeSequences`, `RemoveRun`, 
        `Reverse`, `SegmentationExtract`, `VariableScaling`. 
    """

    # Default variable value
    try:
        nb_variable = obj.get_nb_variable()
        if(nb_variable == 1 and not param2): 
            param2 = param1
            param1 = 1
    except AttributeError:
        pass

    # Manager Parameters
    params = [param1]
    if(param2):
        params.append(param2)
    if(AddVariable):
        params.append(AddVariable)

    return obj.transcode(*params)


def Clustering(matrix, type, *args, **kargs):
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
        if 2 < nb_cluster < (number of patterns), an object of type clusters is returned
        
    :Examples:
      >>> Clustering(dissimilarity_matrix, "Partition", nb_cluster, Prototypes=[1, 3, 12])
      >>> Clustering(dissimilarity_matrix, "Partition", clusters)
      >>> Clustering(dissimilarity_matrix, "Hierarchy", Algorithm="Agglomerative")
      >>> Clustering(dissimilarity_matrix, "Hierarchy", Algorithm="Divisive")
    
    .. seealso::
        `SelectIndividual`,
        `Symmetrize`,
        `Compare`,
        `ToDistanceMatrix`. 
    """
    
    format_map = { "ASCII" :'a',
                   "SpreadSheet": 's',
                   "" : 'n'
                   }

    criterion_map = {
        "FarthestNeighbor" : _stat_tool.FARTHEST_NEIGHBOR,
        "Averaging" : _stat_tool.AVERAGING,
        }

    algorithm_map = {
        "Agglomerative" : _stat_tool.AGGLOMERATIVE,
        "Divisive": _stat_tool.DIVISIVE,
        "Ordering": _stat_tool.ORDERING,
        }
 


    # Get Keywords (set default values)
    Prototypes = kargs.get("Prototypes", [])
    FileName = kargs.get("Filename", "")
    Algorithm = kargs.get("Algorithm", "Agglomerative")
    Criterion = kargs.get("Criterion", "Averaging")
    Initialization = kargs.get("Initialization", 1)
    Format = kargs.get("Format", "")

    
    # Swith for each type of clustering

    if(type == "Partition"):
        if(isinstance(args[0], int)):

            nb_cluster = args[0]
            # WARNING : in this case, algorithm is an int (default is 1)
            if(not isinstance(Algorithm, int)):
                if(kargs.has_key("Algorithm")):
                    raise TypeError("Algorithm must be an int")
                else:
                    Algorithm = 1

            if(len(Prototypes) != nb_cluster):
                raise TypeError("Prototypes must have %i values"%(nb_cluster))
            
            return matrix.partitioning_prototype(nb_cluster, Prototypes, 
                                                 Initialization, Algorithm)
        else:
            return matrix.partitioning_clusters(args[0])

    elif(type == "Hierarchy"):
        
        # Transform Keywords to value
        try:
            Algorithm = algorithm_map[Algorithm]
        except KeyError:
            raise KeyError("Invalid Algorithm. Possible values are : "
                           + str(algorithm_map.keys()))

        try:
            Criterion = criterion_map[Criterion]
        except KeyError:
            raise KeyError("Invalid Criterion. Possible values are : " 
                           + str(criterion_map.keys()))

        try:
            Format = format_map[Format]
        except KeyError:
            raise KeyError("Invalid Format. Possible values are : "
                           + str(format_map.keys()))


        return matrix.hierarchical_clustering(Algorithm, Criterion,
                                              FileName, Format)
        
        
    else:
        raise KeyError("Possible action are : 'Partitioning' or 'Hierarchy'")
                 
        


def ToDistanceMatrix(clusters):
    """    
    Cast and object of type CLUSTERS into an object of type DISTANCE_MATRIX.
 
    
    :Parameters:
      * clusters (clusters) 

    :Returns:    
        An object of type distance_matrix is returned. 

    :Examples:
        >>> ToDistanceMatrix(clusters)

    .. seealso::
        `Clustering`
    """
    
    return _DistanceMatrix(clusters)

################################################################################
from openalea.stat_tool import get_test_file


class Test:
    """Test class to test cluster function and classes"""
    
    def __init__(self):
        pass

    def test_cluster(self):

        from histogram import Histogram
        fagus = Histogram(get_test_file("fagus1.his"))

        histo2 = Cluster(fagus, "Step", 2)
        histo3 = Cluster(fagus, "Information", 0.8)
        histo4 = Cluster(fagus, "Limit", [2, 4, 6, 8, 10])

        assert histo2
        assert histo3
        assert histo4

    def test_transcode(self):

        from histogram import Histogram
        fagus = Histogram(get_test_file("fagus1.his"))

        histo5 = Transcode(fagus, [1, 2, 2, 3, 3, 4, 4, 5])
        
        assert histo5

    
    def test_transcode_err(self):

        from histogram import Histogram
        fagus = Histogram(get_test_file("fagus1.his"))

        try:
            _histo5 = Transcode(fagus, [1, 2, 2, 3, 3, 4, ])
            assert False
        except:
            assert True

    def test_clustering(self):        
        
        from vectors import Vectors, VectorDistance
        from comparison import Compare
        from data_transform import SelectVariable

        vec10 = Vectors(get_test_file("chene_sessile.vec"))
        vec15 = SelectVariable(vec10, [1, 3, 6], Mode="Reject")

        assert vec15

        matrix10 = Compare(vec15, VectorDistance("N", "N", "N"))

        c1 = Clustering(matrix10, "Partition", 3, Prototypes=[1, 3, 12])
        Clustering(matrix10, "Hierarchy", Algorithm="Agglomerative")
        Clustering(matrix10, "Hierarchy", Algorithm="Divisive")

        assert c1
        assert ToDistanceMatrix(c1)


