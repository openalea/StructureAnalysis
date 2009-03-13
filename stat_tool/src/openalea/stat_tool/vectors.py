__doc__ = "Vectors function and class"
__docformat__ = "restructuredtext"

import _stat_tool
import interface

from _stat_tool import _Vectors
from _stat_tool import _VectorDistance

__all__ = ['Vectors',
           '_Vectors',
           'VectorDistance',
           '_VectorDistance',
           'ContingencyTable',
           'VarianceAnalysis']


# Extend dynamically class
interface.extend_class( _Vectors, interface.StatInterface)
interface.extend_class( _VectorDistance, interface.StatInterface)

# Add methods to _Vectors

def _Vectors_mixture_estimation(self, model, nb_iteration=_stat_tool.I_DEFAULT,
                                force_param=[]):
    """Estimate a mixture from _Vectors given initial model or number of components,
    the maximal number of iterations and a flag for using parametric observation
    distributions or not, within a given family"""
    return _stat_tool._Vectors.mixture_estimation_wrap(self, model, nb_iteration, force_param);

_Vectors.mixture_estimation = _Vectors_mixture_estimation


def Vectors(*args):
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
        >>> Vectors(list, Identifiers=[1, 8, 12])
        >>> Vectors(seq, IndexVariable=True)
        >>> Vectors(file_name)
     
    .. seealso::
        `Save`, `ExtractHistogram`, `Cluster`, `Merge`, `MergeVariable`, `SelectIndividual`, 
        `SelectVariable`, `Shift`, `Transcode`, `ValueSelect`, `Compare`, 
        `ComputeRankCorrelation`, `ContingencyTable`, `Regression`, `VarianceAnalysis`
    """
    return _Vectors(*args)



############### VectorDistance #################################################

vector_distance_type = \
    {
    "S":  _stat_tool.SYMBOLIC,
    "SYMBOLIC" : _stat_tool.SYMBOLIC,
    "N": _stat_tool.NUMERIC, 
    "NUMERIC" : _stat_tool.NUMERIC,
    "O": _stat_tool.ORDINAL,
    "ORDINAL":  _stat_tool.ORDINAL,
    "C":  _stat_tool.CIRCULAR,
    "CIRCULAR" : _stat_tool.CIRCULAR,

    _stat_tool.SYMBOLIC: _stat_tool.SYMBOLIC,
    _stat_tool.NUMERIC : _stat_tool.NUMERIC,
    _stat_tool.ORDINAL : _stat_tool.ORDINAL,
    _stat_tool.CIRCULAR : _stat_tool.CIRCULAR,
    }


distance_type = \
    {
    "ABSOLUTE_VALUE" : _stat_tool.ABSOLUTE_VALUE,
    "QUADRATIC" : _stat_tool.QUADRATIC,

    _stat_tool.ABSOLUTE_VALUE: _stat_tool.ABSOLUTE_VALUE,
    _stat_tool.QUADRATIC : _stat_tool.QUADRATIC,
    }



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
        >>> VectorDistance(type1, type2,..., Distance="QUADRATIC")
        >>> VectorDistance(weight1, type1, weight2, type2,..., Distance="QUADRATIC")
        >>> VectorDistance(file_name)



    .. seealso::
        `Compare`
    """

    # filename
    if(len(args)==1 and isinstance(args[0], str)):
        filename = args[0]
        return _stat_tool._VectorDistance.__init__(self, filename)

    # Get keyword parameters
    distance = None
    if(not distance) : distance = kargs.get("Distance", None)
    if(not distance) : distance = _stat_tool.ABSOLUTE_VALUE
    else: distance = distance_type[distance]

    # Parse arguments
    weigths = []
    types = []
    
    cindex = 0
    while cindex < len(args):
        arg = args[cindex]

        # (weigth, type)
        if(isinstance(arg, float) or isinstance(arg, int) ):
            weigths.append(arg)
            types.append(vector_distance_type[args[cindex+1]])
            cindex +=1 

        elif(isinstance(arg, str)):
            weigths.append(0.)
            types.append(vector_distance_type[arg])
            
        elif(isinstance(arg, tuple) and len(arg) == 2):             
            weigths.append(arg[0])
            types.append(vector_distance_type[arg[1]])

        else:
            raise TypeError("Unexpected argument" + str(arg))

        cindex += 1

    return _stat_tool._VectorDistance(types, weigths, distance)
    

################################################################################

def VarianceAnalysis(vec, class_variable, response_variable, 
                     type, FileName="result", Format="ASCII"):
    """
    One-way variance analysis.
    
    :Examples:
    
        >>> VarianceAnalysis(vec, class_variable, response_variable, type, FileName="result", Format="SpreadSheet")
      
    :Parameters:

      * vec (_Vectors),
      * class_variable (int): index of the class or group variable,
      * response_variable (int): index of the response variable,
      * type (string): type of the response variable ("NUMERIC" ("N") or "ORDINAL" ("O")). 

    :Keywords:
    
      * FileName (string): name of the result file,
      * Format (string): format of the result file: "ASCII" (default format) or "SpreadSheet". This optional argument can only be used in conjunction with the optional argument FileName. 

    :Returns:
    
        The variance analysis result as a string

    """

    variance_type = {
        "N": _stat_tool.NUMERIC, 
        "NUMERIC" : _stat_tool.NUMERIC,
        "O": _stat_tool.ORDINAL,
        "ORDINAL":  _stat_tool.ORDINAL,
        }

    try:
        type = variance_type[type]
    except KeyError:
        raise KeyError("Possible type are : " + str(variance_type.keys()))

    return vec.variance_analysis(class_variable, response_variable, type, FileName, Format)



def ContingencyTable(vec, variable1, variable2, 
                     FileName="result", Format="ASCII"):
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
    
        >>> ContingencyTable(vec, variable1, variable2, FileName="result", Format="SpreadSheet")
      
  
    """

    return vec.contingency_table(variable1, variable2, FileName, Format)



################################################################################
from openalea.stat_tool import get_test_file

class Test:

    def test_vectors_pylist(self):

        v = Vectors([[0,1,2,3], [4,5,6,7]])
        assert v

        v = Vectors([[1.2, 0.34], [1.2, 0.34]])
        assert v

        v = Vectors([[1]])
        assert v

        v = Vectors([[0,1,2,3], [4,5,6,7], [1,2,3,4]])
        assert v

        try:
            v = Vectors([[0,1,2,3], [4,5,6,7], [1,2,3]])
            assert False
        except:
            assert True

        try:
            v = Vectors([[0,1,2,3], [4,5,6,7], [1.2,2,3]])
            assert False
        except:
            assert True


    def test_vectors_file(self):
        # File operations

        v1 = Vectors(get_test_file("vectors.vec"))
        v1.save(get_test_file("vectors2.vec"), ViewPoint="Data")

        v2 = Vectors(get_test_file("vectors2.vec"))

        assert v1 and v2
        assert len(v2) == len(v1)

        for i in xrange(len(v2)):
            assert v1[i] == v2[i]



    def test_vectors_display(self):
        # ASCII representation
        v = Vectors([[0,1,2,3], [4,5,6,7]])
        s = str(v)
        assert v.display() == s



    def test_vectors_container(self):

        v = Vectors([[0,1,2,3], [4,5,6,7]])
        assert len(v) == 2

        for i in v:
            assert len(i) == 4

        assert v[0] == range(4)
        assert v[1][1] == 5


    def test_vector_distance(self):

        v = VectorDistance('N', 'O', 'S') 
        assert v and len(v) == 3

        v = VectorDistance('NUMERIC', 'ORDINAL', 'SYMBOLIC')
        assert v and len(v) == 3

        v = VectorDistance(2.3, 'N', 4, 'O', 6, 'S')
        assert v and len(v) == 3

        v = VectorDistance( (2.3, 'N'),  (4, 'O'), (6, 'S'))
        assert v and len(v) == 3

        v = VectorDistance(2.3, 'N', 4, 'O', 6, 'S', distance = "QUADRATIC")
        assert v and len(v) == 3

        v = VectorDistance('NUMERIC', 'ORDINAL', 'SYMBOLIC', Distance="QUADRATIC")
        assert v and len(v) == 3

        assert str(VectorDistance('N', 'O', 'S'))


    def test_variance_analysis(self):
        vec10 = Vectors(get_test_file("chene_sessile.vec"))
        va = VarianceAnalysis(vec10, 1, 4, "O")
        assert va and str(va)


    def test_contingency_table(self):
        vec10 = Vectors(get_test_file("chene_sessile.vec"))
        ct = ContingencyTable(vec10, 1, 4)
        assert ct and str(ct)

        


