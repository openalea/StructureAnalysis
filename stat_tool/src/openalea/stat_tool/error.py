""" Error class for stat_tool """
__revision__ = "$Id$"


from _stat_tool import _FormatError

__all__ = ["_FormatError",]

arguments_labels = {1:'first', 2:'second', 3:'third', 4:'fourth', 5:'fifth'}


class StatToolError(Exception):
    def __init__(self, *args, **kargs):
        #fixme:  is this the correct syntax for the init method ?
        Exception.__init__(self, args, kargs)
        self.nargs = kargs.get("nargs", None)
        
        # error message
        self.msg = kargs.get("msg", None)
        if len(args)==1:
            self.msg = args[0]
        
        # arguments error
        self.min_nargs = kargs.get("min_nargs", None)
        self.max_nargs = kargs.get("max_nargs", None)
        self.nargs_given = kargs.get("nargs_given", None)

    def __str__(self):
        """switch to the error messages."""
        
        # number of arguments case
        if self.msg:
            msg = self.msg
        else:
            msg = 'Unknown error'
        
        return msg


def CheckArgumentsLength(args, min_nargs=None, max_nargs=1e6, optional=False):
    """Check that the number of arguments is expected range.

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
        msg = 'Expected number of %s arguments in the range [%s, %s] but %s provided' \
                % (msg_optional, min_nargs, max_nargs, l)
                
        raise StatToolError(min_nargs=min_nargs, max_nargs=max_nargs, 
                            nargs_given=l, msg=msg)


def CheckOptionalArgumentsLength(args, min_nargs=None, max_nargs=None):
    """Check that the number of arguments is expected range.

    See CheckArgumentsLength function.
    """
    CheckArgumentsLength(args, min_nargs, max_nargs, optional=True)


def CheckType(variable, *args, **kargs):
    """
    utype is a list of types expected
    """
    variable_id = kargs.get("variable_id", 1)

    try:
        label = arguments_labels[variable_id]
    except:
        label = ''
    
    if type(args[0])==list:
      if type(variable) not in args[0]:
            raise StatToolError('Expect %s argument to be a of type %s. Found %s instead.' % 
                            (label, args[0], type(variable)))
    elif type(variable) not in list(args):
        raise StatToolError('Expect %s argument to be a of type %s. Found %s instead.' % 
                            (label, list(args), type(variable)))

def DictWrongKeys(udict, key):
    raise StatToolError('Key %s not found. Possible choices are %s.' % 
                            (key, udict))
    