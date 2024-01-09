""" Interfaces for stat_tool objects 

.. note:: this is for developers usage only, not part of the user library
"""
__version__ = "$Id: interface.py 9033 2010-06-01 12:06:16Z cokelaer $"

from output import StatInterface



def extend_class(cls, *base_class):
    """ Extend boost python class

    :Parameters:

      * `cls` - the class to extend
      * `base_class` - the base class to extend

    :returns:
        the modified cls
    """

    b = list(cls.__bases__)
    for c in base_class:
        if (c not in b):
            b.append(c)
    cls.__bases__ = tuple(b)

    return cls
