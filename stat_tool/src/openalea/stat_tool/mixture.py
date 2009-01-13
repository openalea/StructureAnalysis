__docformat__ = "restructuredtext"
__doc__ = """ Mixture """

import os

import interface
import _stat_tool

from _stat_tool import _Mixture
from _stat_tool import _MixtureData
from _stat_tool import _MvMixture
from _stat_tool import _MvMixtureData

__all__ = ['_Mixture', 
           '_MixtureData',
           '_MvMixture',
           '_MvMixtureData',
           'Mixture',
           ]




def Mixture(*args):
    """
    Construction of a mixture of distributions from elementary distributions and 
    associated weights or from an ASCII file.

    Usage
    -----
      * ``Mixture(weight1, dist1, weight2, dist2,...)``
      * ``Mixture(filename)``

    Parameters
    ----------
      * weight1, weight2, ... (float): weights of each component. \
     These weights should sum to one (they constitute a discrete distribution).
      * dist1, dist2, ... (`_ParametricModel`, `_Mixture`, `_Convolution`, `_Compound`)\ 
      elementary distributions (or components).
      * filename (string). 

    Return
    ------
    If the construction succeeds, an object of type mixture is returned, otherwise no object 
    is returned. 

    Background
    ----------
    A mixture is a parametric model of classification where each elementary distribution or 
    component represents a class with its associated weight. 

    See Also
    --------
    `Save`, `Estimate`, `Simulate`.

    """

    if(len(args)==0) : 
        raise TypeError()

    # filename
    if(len(args)==1) :
        return _stat_tool._Mixture(args[0])

    # build list of weights and distributions
    else:
        nb_param = len(args)

        if((nb_param % 2) != 0) :
            raise TypeError("Number of parameters must be pair")

        weights = []
        dists = []
        
        for i in xrange(nb_param / 2):
            weights.append(args[i * 2])
            dists.append(args[i * 2 + 1])

        return _stat_tool._Mixture(weights, dists)
    


# Extend _Mixture
interface.extend_class( _stat_tool._Mixture, interface.StatInterface)

# Extend _MixtureData
interface.extend_class( _stat_tool._MixtureData, interface.StatInterface)

# Extend _MvMixture
interface.extend_class( _stat_tool._MvMixture, interface.StatInterface)

# Extend _MvMixtureData
interface.extend_class( _stat_tool._MvMixtureData, interface.StatInterface)

# Add methods to _MvMixture

def _MvMixture_old_plot(self, variable, Title=""):
    """Plot a given variable"""
    if ((variable < 0) or (variable >= self.nb_variable())):
        raise IndexError, "variable index out of range: "+str(variable)
    file_id = str(variable+1)
    if (not self._is_parametric(variable)):
        file_id += "0"
    interface.StatInterface.old_plot(self, Title=Title, Suffix=file_id)

def _MvMixture_get_plotable(self):
    """Return a plotable object (not yet implemented)"""
    return None

def _MvMixture_criteria(self):
    """Extract the value of each selection criterion"""
    disp = self.display()
    criteria = {}
    names = ["AIC", "AICc", "BIC", "BICc"]
    for name in names:
        f = disp.find("(" + name + "):")
        if (f != -1):
            pos = f + len(name) + 3
            i = disp.find("\n", pos)
            try:
                val = float(disp[pos:i])
            except ValueError:
                pass
            else:
                if str(val).upper() != "NAN":
                    criteria[name] = val
    return criteria

_stat_tool._MvMixture.save_backup = _stat_tool._MvMixture.save

def _MvMixture_save(self, file_name, format="ASCII", overwrite=False):
        """Save MvMixture object into a file.
        
        Argument file_name is a string designing the file name and path.
        String argument format must be "ASCII" or "SpreadSheet".
        Boolean argument overwrite is false is the file should not 
        be overwritten."""
        if not overwrite:
            try:
                f = file(file_name, 'r')
            except IOError:
                f = file(file_name, 'w+')
            else:
                msg = "File " + file_name + " already exists"
                raise IOError, msg
            f.close()
        import string
        if not (string.upper(format)=="ASCII" 
                or string.upper(format)=="SPREADSHEET"):
            msg = "unknown file format: " + str(format)
            raise ValueError, msg
        else:
            try:
                _stat_tool._MvMixture.save_backup(self, file_name, Detail=1, 
                                                  ViewPoint='', Format=format)
            except RuntimeError, error:
                raise FormatError, error

_stat_tool._MvMixture.state_permutation_backup = _stat_tool._MvMixture.state_permutation

def _MvMixture_state_permutation(self, perm):
  """Permutation of the states of self.
  perm[i]==j means that current state i will become new state j.
        
  Usage:  state_permutation(list)"""
  self.state_permutation_backup(perm)

_MvMixture.old_plot = _MvMixture_old_plot
_MvMixture.get_plotable = _MvMixture_get_plotable
_MvMixture._criteria = _MvMixture_criteria
_MvMixture.save = _MvMixture_save
_MvMixture.state_permutation = _MvMixture_state_permutation

# Add methods to _MvMixtureData

def _MvMixtureData_old_plot(self, variable, Title=""):
    """Plot a given variable"""
    if ((variable < 0) or (variable >= self.get_nb_variable())):
        raise IndexError, "variable index out of range: "+str(variable)
    file_id = str(variable+1)
    interface.StatInterface.old_plot(self, Title=Title, Suffix=file_id)

def _MvMixtureData_get_plotable(self):
    """Return a plotable object (not yet implemented)"""
    return None

_MvMixtureData.old_plot = _MvMixtureData_old_plot
_MvMixtureData.get_plotable = _MvMixtureData_get_plotable

########################## Test Mixture ########################################
from openalea.stat_tool import get_test_file

class Test:
    def test_emty(self):

        try:
            m = Mixture()
            assert False

        except Exception:
            assert True


    def test_file(self):

        m = Mixture(get_test_file("mixture1.mixt"))
        assert m


    def test_build_mixture(self):

        from distribution import Uniform


        d1 = Uniform(0, 10)
        d2 = Uniform(10, 20)
        d3 = Uniform(20, 30)

        m = Mixture(0.1, d1, 0.2, d2, 0.7, d3)
        assert m
        return m


    def __test_plot(self):

        m = self.test_build_mixture()
        #    m.plot()

        assert str(m)
        m.display()


    def test_simulation(self):

        m = self.test_build_mixture()
        s = m.simulate(10)

        assert s.nb_component() == 3
        assert str(s)


    def test_extract(self):

        from data_transform import ExtractDistribution
        from distribution import Uniform


        m = self.test_build_mixture()
        assert m.nb_component() == 3

        assert m.extract_weight() == ExtractDistribution(m, "Weight")

        assert m.extract_mixture() == ExtractDistribution(m, "Mixture")

        assert ExtractDistribution(m, "Component", 1) == Uniform(0, 10)
        assert ExtractDistribution(m, "Component", 2) == Uniform(10, 20)
        assert ExtractDistribution(m, "Component", 3) == Uniform(20, 30)

        assert m.extract_component(1) == Uniform(0, 10)
        assert m.extract_component(2) == Uniform(10, 20)
        assert m.extract_component(3) == Uniform(20, 30)


    def test_extract_data(self):

        from histogram import Histogram 

        h = Histogram(get_test_file("meri2.his"))
        m = h.estimate_mixture(["B", "NB"])

        d = m.extract_data()
        assert d


    def test_build_mv_mixture(self):

        from distribution import Binomial, Poisson
        from mixture import _MvMixture

        d11 = Binomial(0, 12, 0.1)
        d12 = Binomial(0, 12, 0.5)
        d13 = Binomial(0, 12, 0.8)

        d21 = Poisson(0, 18.0)
        d22 = Poisson(0, 5.0)
        d23 = Poisson(0, .20)
        
        m = _MvMixture([0.1, 0.2, 0.7], [[d11, d21], [d12, d22], [d13, d23]])
        assert m
        return m
    
    def test_mv_fromfile(self):

        from mixture import _MvMixture

        # From file
        m = _MvMixture(get_test_file("mixture_mv1.mixt"))
        assert m
    
        # File
        m.save("test_mv.mixt")

        m1 = _MvMixture(get_test_file("mixture_mv1.mixt"))
        m2 = _MvMixture("test_mv.mixt")
        assert m.nb_component() == m2.nb_component()
        assert str(m) == str(m2)
        
        os.remove("test_mv.mixt")    

        mnp = _MvMixture(get_test_file("mixture_mv_nonparam.mixt"))
        assert m

        try:
            h = _MvMixture(get_test_file("no_such_file.mixt"))
            assert False
        except Exception:
            assert True

    def test_simulate_estimate_mv_mixture(self):

        from distribution import Binomial, Poisson
        from mixture import _MvMixture
        
        d11 = Binomial(0, 12, 0.1)
        d12 = Binomial(2, 13, 0.6)
        d13 = Binomial(3, 15, 0.9)

        d21 = Poisson(0, 25.0)
        d22 = Poisson(0, 5.0)
        d23 = Poisson(0, 0.2)

        
        m = _MvMixture([0.1, 0.2, 0.7], [[d11, d21], [d12, d22], [d13, d23]])
        v = m.simulate(5000);
        assert v

        m_estim_model = v.mixture_estimation(m, 100,  [True, True]);
        assert m_estim_model
        m_estim_nbcomp = v.mixture_estimation(2);
        assert m_estim_nbcomp
