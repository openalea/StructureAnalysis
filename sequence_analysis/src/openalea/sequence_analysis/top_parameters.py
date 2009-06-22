"""Top parameters

author: Thomas Cokelaer, Thomas.Cokelaer@inria.fr

"""
__revision__ = "$Id: $"

import os
import openalea.stat_tool.interface as interface
from openalea.sequence_analysis._sequence_analysis import _Top_parameters

__all__ = ['Top_parameters',
           '_Top_parameters']


# Extend dynamically class
interface.extend_class( _Top_parameters, interface.StatInterface)

# Add methods to _Vectors


def Top_parameters(*args, **kargs):
    """TopParameters

    Construction of 'top' parameters from the three parameters or from an ASCII file.

    :Usage:

    >>> TopParameters(proba, axillary_proba, rhythm_ratio, MaxPosition=40)
    >>> TopParameters(file_name, MaxPosition=40)    


    :Arguments:

    * proba (int, real): growth probability of the parent shoot,
    * axillary_proba (int, real): growth probability of the offspring shoots,
    * rhythm_ratio (int, real): growth rhythm ratio offspring shoots / parent shoot,
    * file_name (string).
    
    :Optional Arguments:

    MaxPosition (int): maximum position for the computation of the distributions of the number of internodes of offspring shoots (default value: 20).
    
    :Returned Object:

    If the construction succeeds, an object of type top_parameters is returned, otherwise no object is returned.
    
    :Background:

    The aim of the model of 'tops' is to related the growth of offspring shoots to the growth of their parent shoot in the case of immediate branching. In the case where the arguments are the three 'top' parameters, the constraints over these parameters are described in the definition of the syntactic form of the type top_parameters (cf. File Syntax).
    
    .. seealso::

        :func:`~openalea.stat_tool.output.Save`, 
        :func:`~openalea.sequence_analysis.simulate.Simulate`.

    """
    
    if len(args)==1: 
        #filename case
        if isinstance(args[0], str):
            filename = args[0]
            if os.path.isfile(filename):
                return _Top_parameters(filename, True)
            else:
                raise IOError("bad file name")
        #sequences case
    else:
        raise TypeError("Expected a valid filename or a list of lists (e.g., [[1,0],[0,1]])")
    
    

 
    


