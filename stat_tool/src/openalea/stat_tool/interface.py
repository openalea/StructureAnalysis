""" Interfaces for stat_tool objects """
#__docformat__ = "restructuredtext"



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


from output import StatInterface    

        
