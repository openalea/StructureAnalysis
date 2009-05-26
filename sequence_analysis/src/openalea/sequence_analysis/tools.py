""" common tools to all modules

.. author:: Thomas Cokelaer, Thomas.Cokelaer@inria.fr 
"""
__revision__ = "$Id:  $"




def __parse_kargs__(kargs, key, default=None, dict_map=None):
    """
    Check the presence of key in kargs. 
    
    
    :param kargs: a list of arguments and their values
    :param key:" a key to look for in the kargs variable
    :param default: the default value of the key if not found in kargs
    :param dict_map: if provided, returns the value of the key with respect to this dictionary
    
    :return: 
    
    * if dict_map is provided returns the value corresponding to key found in dict_map. 
    * if dict_map is none, returns the value of key found in kargs (typically int, bool)
    """
    
    
    user_choice = kargs.get(key, default)
    
    if dict_map:
        try:
            return dict_map[user_choice]
        except KeyError:    
            raise KeyError("Wrong choice for %. Possible choices are %s " % 
                           (key, dict_map.keys()))
    else:
        return user_choice
