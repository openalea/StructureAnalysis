"""Mv Mixture class"""
__version__ = "$Id$"


import interface
import error


from _stat_tool import _MvMixture
from _stat_tool import _MvMixtureData

__all__ = ['_MvMixture',
           '_MvMixtureData']



# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# THW FOLLOWING HAS TO BE CHECK WITH JBD + tests 

# Extend _MvMixture
interface.extend_class(_MvMixture, interface.StatInterface)

# Extend _MvMixtureData
interface.extend_class(_MvMixtureData, interface.StatInterface)

# Add methods to _MvMixture


def _MvMixture_old_plot(self, variable, Title=""):
    """Plot a given variable"""
    if ((variable < 0) or (variable >= self.nb_variable)):
        raise IndexError, "variable index out of range: " + str(variable)
    file_id = str(variable + 1)
    if (not self._is_parametric(variable)):
        file_id += "0"
    interface.StatInterface.old_plot(self, Title=Title, Suffix=file_id)

def _MvMixture_old_print(self, variable, Title=""):
    """Print a given variable into .ps file"""
    if ((variable < 0) or (variable >= self.nb_variable)):
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

_MvMixture.save_backup = _MvMixture.save


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

    if not (format.upper()=="ASCII"
           or format.upper()=="SPREADSHEET"):
        msg = "unknown file format: " + str(format)
        raise ValueError, msg
    else:
        try:
            _MvMixture.save_backup(self, file_name, Detail=1,
                                                  ViewPoint='', Format=format)
        except RuntimeError, my_error:
            raise my_error

_MvMixture.state_permutation_backup = _MvMixture.state_permutation


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
