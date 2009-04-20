""" Mixture class"""
__revision__ = "$Id$"


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
           'Mixture', ]


def Mixture(*args):
    """Construction of a mixture of distributions from elementary distributions
    and associated weights or from an ASCII file.

    A mixture is a parametric model of classification where each elementary
    distribution or component represents a class with its associated weight.

    :Parameters:
      * `weight1`, `weight2`, ... (float) - weights of each component.
         These weights should sum to one (they constitute a discrete
         distribution).
      * `dist1`, `dist2`, ... (`_ParametricModel`, `_Mixture`, `_Convolution`,
        `_Compound`) elementary distributions (or components).
      * `filename` (string) -

    :Returns:
        If the construction succeeds, an object of type mixture is returned,
        otherwise no object is returned.

    :Examples:

    .. doctest::
        :options: +SKIP

        >>> Mixture(weight1, dist1, weight2, dist2,...)
        >>> Mixture(filename)

    .. seealso::
        :func:`~openalea.stat_tool.output.Save`,
        :func:`~openalea.stat_tool.estimate.Estimate`,
        :func:`~openalea.stat_tool.simulate.Simulate`.

    """

    if (len(args)==0):
        raise TypeError()

    # filename
    if (len(args)==1):
        return _stat_tool._Mixture(args[0])

    # build list of weights and distributions
    else:
        nb_param = len(args)

        if ((nb_param % 2) != 0):
            raise TypeError("Number of parameters must be pair")

        weights = []
        dists = []

        for i in xrange(nb_param / 2):
            weights.append(args[i * 2])
            dists.append(args[i * 2 + 1])

        return _stat_tool._Mixture(weights, dists)

# Extend _Mixture
interface.extend_class(_stat_tool._Mixture, interface.StatInterface)

# Extend _MixtureData
interface.extend_class(_stat_tool._MixtureData, interface.StatInterface)

# Extend _MvMixture
interface.extend_class(_stat_tool._MvMixture, interface.StatInterface)

# Extend _MvMixtureData
interface.extend_class(_stat_tool._MvMixtureData, interface.StatInterface)

# Add methods to _MvMixture


def _MvMixture_old_plot(self, variable, Title=""):
    """Plot a given variable"""
    if ((variable < 0) or (variable >= self.nb_variable())):
        raise IndexError, "variable index out of range: " + str(variable)
    file_id = str(variable + 1)
    if (not self._is_parametric(variable)):
        file_id += "0"
    interface.StatInterface.old_plot(self, Title=Title, Suffix=file_id)

def _MvMixture_old_print(self, variable, Title=""):
    """Print a given variable into .ps file"""
    if ((variable < 0) or (variable >= self.nb_variable())):
        raise IndexError, "variable index out of range: "+str(variable)
    file_id = str(variable+1)
    if (not self._is_parametric(variable)):
        file_id += "0"
    interface.StatInterface.plot_print(self, Title=Title, Suffix=file_id)

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
        be overwritten.
        """
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

    Usage:  state_permutation(list)
    """
    self.state_permutation_backup(perm)

_MvMixture.old_plot = _MvMixture_old_plot
_MvMixture.plot_print = _MvMixture_old_print
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
