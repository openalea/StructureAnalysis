""" Error class for stat_tool """
__version__ = "$Id$"


#from _stat_tool import _FormatError

#__all__ = ["_FormatError",]

arguments_labels = {1:'first', 
                    2:'second', 
                    3:'third', 
                    4:'fourth', 
                    5:'fifth'}



def CheckArgumentsLength(args, min_nargs=0, max_nargs=1e6, optional=False):
    """Check that the number of arguments is in the expected range.

    :param: min is the minimum number of arguments expected.
    :param: max is the maximum number of arguments expected (strict).
    
    if max is not provided, we assume that it can be any value
    """
    l = len(args)

    # todo: check that min_args is not none and raise expcetion otherwise
         
    if (l < min_nargs) or (l > max_nargs):
        if optional: 
            msg_optional = 'OPTIONAL'
        else:
            msg_optional = ''
        msg = """Expected number of %s arguments in the range [%s, %s]  
            but %s provided' """ % (msg_optional, min_nargs, max_nargs, l)
                
        raise Exception(msg)


def CheckOptionalArgumentsLength(args, min_nargs=None, max_nargs=None):
    """Check that the number of arguments is expected range.

    See CheckArgumentsLength function.
    """
    CheckArgumentsLength(args, min_nargs, max_nargs, optional=True)


def CheckType(variables, types, **kargs):
    """
    
    :param variable: a variable to test
    :param types: list of types or unique type (e.g., int, str)
    
    optional argument: variable_id corresponding to position of variable
    
    #CheckType(1, int, variable_index=0) NOT IMPLEMENTED 
    CheckType([1,'a'], [[int, float],str], variable_index=[0,1])

    """
    
    variable_id = kargs.get("variable_id", None)
    
    if type(variables)!=list:
        raise TypeError("first argument(variables) must be a list")
        
    if type(types)!=list:
        raise TypeError("second argument(variables) must be a list")

    if type(variable_id)!=list and variable_id is not None:
        raise TypeError("second argument(variables) must be a list")

        
    n = len(variables)
    #populate labels from 0 to N or with contents of variable_id 
    labels = []
    if variable_id is None:
        for index in range(1, n+1):
            labels.append(arguments_labels[index])
    else:
        for var in variable_id:
            if var in arguments_labels.keys():
                labels.append(arguments_labels[var])
            else:
                labels.append("unknown position (larger than 5)")
    
    if len(variables)!=len(types):
        raise ValueError('length of first and second arguments must be equal')


    for variable, mytype, label in zip(variables, types, labels):
        if type(mytype)!=list:
            mytype = [mytype]
        if type(variable) not in mytype:
            raise Exception(
                    """Expect %s argument's type to be in %s. 
                    \nFound %s instead.""" % (label, mytype, type(variable)))
    

def CheckDictKeys(key, udict):
    if key not in udict.keys():
        raise KeyError('Key %s not found. Possible choices are %s.' % 
                            (key, udict.keys()))


def ParseKargs(kargs, keyword, default=None, possible=None):
    """
    :param possible: a list of possible values that kargs[keyword] can take
    
    :Example:
        distance = error.ParseKargs(kargs, "Distance", "ABSOLUTE_VALUE", 
                                distance_type.keys())
        distance = distance_type[distance]
        
        or simply (if a dictionary is used as fourth argument)
        
        distance = error.ParseKargs(kargs, "Distance", "ABSOLUTE_VALUE", 
                                distance_type)
        
        
    """
    ret = kargs.get(keyword, default)
    if possible and type(possible)==list:
        if ret not in possible:
            raise ValueError("keyword **%s** can only take those values %s." 
                             % (keyword, possible))
            
    elif possible and type(possible)==dict:
        if ret not in possible.keys():
            raise ValueError("keyword **%s** can only take those values %s." 
                             % (keyword, possible.keys()))
        ret = possible[ret]
        
    return ret
            
    
    

def CheckKargs(kargs, possible_kargs, dicts=None):
    """Parse kargs argument and check that it belongs to possible_kargs
    
    :param kargs: list of arguments
    :param possible_kargs: list of values that can take kargs
    :param dicts: list of dictionaries possible for each kargs
    
    :Usage:
        Let us assume a function which prototypes is::
    
            def TestFunction(kargs):
            
        Then, we called it as follows::
            
            TestFunction(firstname="James", surname="Brown")
        
        This function only accept *firstname* in ["James", "Roger"] and
        *surname* in ["Brown", "Moore"].
        
        We would add the following code:: 
        
            allowed_firstmame = ["James", "Roger"]
            allowed_surname = ["Brown", "Moore"]
            CheckArgs(kargs, ["firstname", "surname"], 
                        allowed=[allowed_firstname, allowed_surname])
    """
    # check that number of arguments is correct
    CheckArgumentsLength(kargs, 0, len(possible_kargs))

    if dicts: 
        assert len(possible_kargs)==len(dicts)
        
    # check that argument names are correct
    for karg in kargs:
        if karg not in possible_kargs:
            raise KeyError('Argument %s not allowed (%s).', 
                                karg, possible_kargs)
    
    # check that values of the arguments are allowed
    if dicts:
        for karg, mydict in zip(possible_kargs, dicts):
            if kargs.get(karg) not in mydict.keys() and \
                kargs.get(karg) is not None:
                raise ValueError("""
                    Value of the Argument \"%s\" (given %s) not found in 
                    the list of allowed values (%s)""" 
                    % (karg, kargs.get(karg), mydict.keys()))
                
                
def _myreturn(ret, func_name, msg=None):
    """ specialized version of *return*
    
    Check if the object to be returned is None. 
    If so, returns the calling function (func_name) and an error message.
    
    """
    
    if msg == None:
        msg = "Function %s did not return anything. Check your arguments" \
            % func_name
    if ret:
        return ret
    else:
        raise Exception(msg)                
                    
                    
STAT_TOOL_ERROR_MSG_RETURN_NONE = "Function did not return anything. Check your arguments"
STAT_TOOL_NB_VARIABLE_ERROR = \
    """Extra arguments provided (to specify variable value ?). Consider removing 
    it. Be aware that nb_variable equals 1"""
