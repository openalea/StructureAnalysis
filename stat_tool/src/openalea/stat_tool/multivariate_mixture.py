# -*- coding: utf-8 -*-
"""Mv Mixture class"""
__version__ = "$Id$"


import interface
import error


from _stat_tool import _MultivariateMixture
from _stat_tool import _MultivariateMixtureData

__all__ = ['_MultivariateMixture',
           '_MultivariateMixtureData']



# Extend _MultivariateMixture
interface.extend_class(_MultivariateMixture, interface.StatInterface)

# Extend _MultivariateMixtureData
interface.extend_class(_MultivariateMixtureData, interface.StatInterface)

# Add methods to _MultivariateMixture


def _MultivariateMixture_old_plot(self, variable, Title=""):
    """Plot a given variable"""
    if ((variable < 0) or (variable >= self.nb_variable)):
        raise IndexError, "variable index out of range: " + str(variable)
    file_id = str(variable + 1)
    if (not self._is_parametric(variable)):
        file_id += "0"
    interface.StatInterface.old_plot(self, Title=Title, Suffix=file_id)

def _MultivariateMixture_old_print(self, variable, Title=""):
    """Print a given variable into .ps file"""
    if ((variable < 0) or (variable >= self.nb_variable)):
        raise IndexError, "variable index out of range: "+str(variable)
    file_id = str(variable+1)
    if (not self._is_parametric(variable)):
        file_id += "0"
    interface.StatInterface.plot_print(self, Title=Title, Suffix=file_id)

def _MultivariateMixture_get_plotable(self):
    """Return a plotable object (not yet implemented)"""
    return None


def _MultivariateMixture_criteria(self):
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

_MultivariateMixture.save_backup = _MultivariateMixture.save


def _MultivariateMixture_save(self, file_name, format="ASCII", overwrite=False):
    """Save MultivariateMixture object into a file.

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

    if not (format.upper()=="ASCII"
           or format.upper()=="SPREADSHEET"):
        msg = "unknown file format: " + str(format)
        raise ValueError, msg
    else:
        try:
            _MultivariateMixture.save_backup(self, file_name, Detail=1,
                                                  ViewPoint='', Format=format)
        except RuntimeError, my_error:
            raise my_error

_MultivariateMixture.state_permutation_backup = _MultivariateMixture.state_permutation


def _MultivariateMixture_state_permutation(self, perm):
    """Permutation of the states of self.
    perm[i]==j means that current state i will become new state j.

    Usage:  state_permutation(list)
    """
    self.state_permutation_backup(perm)

_MultivariateMixture.old_plot = _MultivariateMixture_old_plot
_MultivariateMixture.plot_print = _MultivariateMixture_old_print
_MultivariateMixture.get_plotable = _MultivariateMixture_get_plotable
_MultivariateMixture._criteria = _MultivariateMixture_criteria
_MultivariateMixture.save = _MultivariateMixture_save
_MultivariateMixture.state_permutation = _MultivariateMixture_state_permutation

# Add methods to _MultivariateMixtureData


def _MultivariateMixtureData_old_plot(self, variable, Title=""):
    """Plot a given variable"""
    if ((variable < 0) or (variable >= self.nb_variable)):
        raise IndexError, "variable index out of range: "+str(variable)
    file_id = str(variable+1)
    interface.StatInterface.old_plot(self, Title=Title, Suffix=file_id)


def _MultivariateMixtureData_get_plotable(self):
    """Return a plotable object (not yet implemented)"""
    return None

_MultivariateMixtureData.old_plot = _MultivariateMixtureData_old_plot
_MultivariateMixtureData.get_plotable = _MultivariateMixtureData_get_plotable
