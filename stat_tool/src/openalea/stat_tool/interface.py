__doc__ = """ Interfaces for stat_tool objects """
__docformat__ = "restructuredtext"



def extend_class(cls, *base_class):
    """ Extend boost python class 
    cls : the class to extend
    base_class : the base class to extend

    return : the modified cls
    """
    
    b = list(cls.__bases__)
    for c in base_class:
        b.append(c)
    cls.__bases__ = tuple(b)

    return cls




from output import StatInterface    

        
