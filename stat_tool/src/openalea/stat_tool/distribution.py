__doc__ = """ Distributions """
__docformat__ = "restructuredtext"

import _stat_tool
import interface


# Public members
from _stat_tool import _Distribution 
from _stat_tool import _ParametricModel
from _stat_tool import DistributionIdentifier

__all__ = ["_Distribution", 
           "_ParametricModel",
           "Distribution",
           "Binomial",
           "Poisson",
           "Uniform",
           "NegativeBinomial",
           "ToHistogram",
           "ToDistribution",
           "distribution_type",
           "DistributionIdentifier"
           ]



distribution_type = \
    {
    "B":  _stat_tool.BINOMIAL,
    "BINOMIAL" : _stat_tool.BINOMIAL,
    "P": _stat_tool.POISSON, 
    "POISSON" : _stat_tool.POISSON,
    "NB": _stat_tool.NEGATIVE_BINOMIAL,
    "NEGATIVE_BINOMIAL":  _stat_tool.NEGATIVE_BINOMIAL,
    "U":  _stat_tool.UNIFORM,
    "UNIFORM" : _stat_tool.UNIFORM,

    _stat_tool.BINOMIAL: _stat_tool.BINOMIAL,
    _stat_tool.POISSON : _stat_tool.POISSON,
    _stat_tool.UNIFORM : _stat_tool.UNIFORM,
    _stat_tool.NEGATIVE_BINOMIAL : _stat_tool.NEGATIVE_BINOMIAL,
    }


def get_distribution_type(typeid, filter=None):
    """ Return a distribution type constante corresponding to type 
    typeid must be in filter list
    """

    # Convert distribution type 
    if(isinstance(typeid, str)):
        typeid = typeid.upper()

    try:
        typeid = distribution_type[typeid]
        if(filter and not typeid in filter):
            raise TypeError("Bad distribution type %s. Should be %s"%(typeid, str(filter)))
        return typeid
    except KeyError:
        error = "Invalid Distribution type, possible type are %s"%(str(distribution_type.keys()),)
        raise AttributeError(error)



def Distribution(type_or_filename, *args):
    """
    Construction of a parametric discrete distribution (either binomial, 
    Poisson, negative binomial or uniform) from the name and the parameters of the 
    distribution or from an ASCII file.

    Usage
    -----
      * ``Distribution("BINOMIAL", inf_bound, sup_bound, proba)``
      * ``Distribution("POISSON", inf_bound, param)``
      * ``Distribution("NEGATIVE_BINOMIAL", inf_bound, param, proba)``
      * ``Distribution("UNIFORM", inf_bound, sup_bound)``
      * ``Distribution(file_name)``

    Parameters
    ----------
      * inf_bound (int): lower bound to the range of possible values (shift parameter),
      * sup_bound (int): upper bound to the range of possible values \
      (only relevant for binomial or uniform distributions),
      * param (int, real): parameter of either the Poisson distribution or \
      the negative binomial distribution.
      * proba (int, float): probability of success \
      (only relevant for binomial or negative binomial distributions),
      * file_name (string). 
    
      Remark: the names of the parametric discrete distributions can be summarized \
    by their first letters: 

        * "B" ("BINOMIAL"), 
        * "P" ("POISSON"),  
        * "NB" ("NEGATIVE_BINOMIAL"), 
        * "U" ("UNIFORM"),

    
    Return
    ------
    If the construction succeeds, an object of type distribution is returned, 
    otherwise no object is returned. 

    Description
    -----------
    A supplementary shift parameter (argument inf_bound) is required with respect 
    to the usual definitions of these discrete distributions. Constraints over 
    parameters are given in the file syntax corresponding to the type distribution 
    (cf. File Syntax). 

    See Also
    --------
    `Save`, `Estimate`, `Simulate`.
    """
        
    typeid = type_or_filename

    # Filename
    if( len(args) == 0 ):
        filename = typeid
        return _ParametricModel(filename)

    # Convert distribution type 
    typeid = get_distribution_type(typeid)

    if(typeid == _stat_tool.BINOMIAL):
        return Binomial(*args)
    
    elif(typeid == _stat_tool.POISSON):
        return Poisson(*args)
    
    elif(typeid == _stat_tool.NEGATIVE_BINOMIAL):
        return NegativeBinomial(*args)
            
    elif(typeid == _stat_tool.UNIFORM):
        return Uniform(*args)



def Binomial(inf_bound, sup_bound=_stat_tool.I_DEFAULT, proba=_stat_tool.D_DEFAULT):
    """
    Construction of a binomial distribution
    
    Arguments:
      inf_bound (int) : lower bound to the range of possible values (shift parameter)
      sup_bound (int) : upper bound to the range of possilbe values 
      proba (int, float) :  probability of 'success' 
    """

    # Check parameters
    assert inf_bound < sup_bound
    assert (sup_bound - inf_bound) < _stat_tool.MAX_DIFF_BOUND
    assert proba <=1. and proba >=0

    param = _stat_tool.D_DEFAULT
    return _ParametricModel(_stat_tool.BINOMIAL, inf_bound, sup_bound, param, proba)

    

def Poisson(inf_bound, param=_stat_tool.D_DEFAULT):
    """
    Construction of a poisson distribution

    Arguments:
        inf_bound (int) : lower bound to the range of possible values (shift parameter)
        param (int, float) : parameter of the Poisson distribution
    """

    assert param > 0. and param < _stat_tool.MAX_MEAN

    sup_bound = _stat_tool.I_DEFAULT
    proba = _stat_tool.D_DEFAULT

    return _ParametricModel(_stat_tool.POISSON, inf_bound, sup_bound, param, proba)



def NegativeBinomial(inf_bound, param=_stat_tool.D_DEFAULT, proba=_stat_tool.D_DEFAULT):
    """
    Construction of a negative binomial distribution
    
    Arguments:
        inf_bound (int) : lower bound to the range of possible values (shift parameter)
        param (int, float) : parameter of the Poisson distribution
        proba (int, float) :  probability of 'success' 
    """

    # Check parameters
    assert proba <=1. and proba >=0
    assert (param * (1. - proba) / proba) <  _stat_tool.MAX_MEAN
        
    sup_bound = _stat_tool.I_DEFAULT

    return _ParametricModel(_stat_tool.NEGATIVE_BINOMIAL, inf_bound, sup_bound, param, proba)


def Uniform(inf_bound, sup_bound=_stat_tool.I_DEFAULT):
    """
    Construction of a uniform distribution

    Arguments:
        inf_bound (int) : lower bound to the range of possible values (shift parameter)
        sup_bound (int) : upper bound to the range of possilbe values 
    """

    # Check parameters

    assert inf_bound < sup_bound
    assert (sup_bound - inf_bound) < _stat_tool.MAX_DIFF_BOUND
    
    param = _stat_tool.D_DEFAULT
    proba = _stat_tool.D_DEFAULT
        
    return _ParametricModel(_stat_tool.UNIFORM, inf_bound, sup_bound, param, proba)




# Extend _ParametricModel
interface.extend_class( _stat_tool._ParametricModel, interface.StatInterface)
 


# Cast Functions

def ToDistribution(histo):
    """    
    Cast an object of type `_DistributionData` into an object of type `_Distribution`.

    Usage
    -----
      * ``ToDistribution(histo)``
    
    Parameters
    ----------
      * histo (DistributionData)

    Return
    ------
    If the object histo contains a 'model' part, an object of type `_Distribution` 
    is returned, otherwise no object is returned. 

    See Also
    --------
      `ToHistogram`
    """
    return histo.extract_model()


def ToHistogram(dist):
    """ 
    Cast an object of type `_Distribution` into an object of type `_DistributionData`.

    Usage
    -----
      * ``ToHistogram(dist)``
    
    Parameters
    ----------
      * dist (distribution). 
    
    Return
    ------
    If the object dist contains a 'data' part, an object of type `_DistributionData` 
    is returned, otherwise no object is returned. 

    See Also
    --------
    `ToDistribution`
    """
    return dist.extract_data()


############################# Test Distribution #################################

class TestDistribution:

    def test_constructor(self):
        pass


    def test_ascii(self):
        # ASCII representation
        #d = Distribution("../../../test/meri1.his")
        #assert d

        d = Distribution ("UNIFORM", 0, 10)

        s = str(d)
        assert d.display() == s


    def test_uniform(self):

        d = Distribution ("UNIFORM", 0, 10)
        assert list(d.simulate(1000))

        d = Uniform (0, 10)
        assert list(d.simulate(1000))

        m = d.simulate(1000).extract_model()
        assert isinstance(m, _stat_tool._ParametricModel)


    def test_binomial(self):

        d = Distribution ("BINOMIAL", 0, 10, 0.5)
        assert list(d.simulate(1000))

        d = Binomial(0, 10, 0.5)
        assert list(d.simulate(1000))

        m = d.simulate(1000).extract_model()
        assert isinstance(m, _stat_tool._ParametricModel)


    def test_poisson(self):

        d = Distribution ("POISSON", 0, 2)
        assert list(d.simulate(1000))

        d = Poisson(0, 2)
        assert list(d.simulate(1000))

        m = d.simulate(1000).extract_model()
        assert isinstance(m, _stat_tool._ParametricModel)


    def test_neg_binomial(self):
        d = Distribution ("NEGATIVE_BINOMIAL", 0, 1, 0.5)
        assert list(d.simulate(1000))

        d = NegativeBinomial (0, 1, 0.5)
        assert list(d.simulate(1000))

        m = d.simulate(1000).extract_model()
        assert isinstance(m, _stat_tool._ParametricModel)


    def test_plot(self):

        d = Distribution ("POISSON", 0, 2)
        d.plot()


    def test_fromnothing(self):

        # From nothing
        try:
            d = Distribution()
            assert False

        except Exception:
            assert True


    def test_to_histogram(self):

        d = Distribution ("NEGATIVE_BINOMIAL", 0, 1, 0.5)
        h  = d.simulate(1000)
        d2 = ToDistribution(h)
        assert h and d2

        h2 = ToHistogram(d2)
        assert h2
        assert h == h2
