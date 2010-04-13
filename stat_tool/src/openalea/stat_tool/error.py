# -*- coding: utf-8 -*-
""" Error class to manage standard errors in stat_tool """
__version__ = "$Id$"


#from _stat_tool import _FormatError

#__all__ = ["_FormatError",]

arguments_labels = {1:'first',
                    2:'second',
                    3:'third',
                    4:'fourth',
                    5:'fifth',
                    6:'sixth',
                    7:'seventh',
                    8:'eighth',
                    9:'ninth',
                    10:'tenth'
                    }



def CheckArgumentsLength(args, min_nargs=0, max_nargs=1e6, optional=False):
    """Check that input's length (tuple) is within a given range

    :param args: a tuple containing user arguments
    :param min: minimum number of arguments expected.
    :param max: maximum number of arguments expected (strict).

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



def CheckType(variables, types, **kargs):
    """Check types of input variables

    :param variables: a list of variables to be checked
    :param types: list of types

    :optional argument:

    variable_id corresponding to the position of each variable.

    :example:

        >>> #CheckType(1, int, variable_index=0) NOT IMPLEMENTED
        >>> CheckType([1,'a'], [int,str], variable_id=[0,1])
        >>> CheckType([1,'a'], [[int, float],str], variable_id=[0,1])

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
            try:
                labels.append(arguments_labels[index])
            except:
                raise KeyError("""consider adding keys in the labels dictionary:
                 you've reached the limit""")
    else:
        for var in variable_id:
            if var in arguments_labels.keys():
                labels.append(arguments_labels[var])
            else:
                labels.append("unknown position (larger than 5)")

    if len(variables) != len(types):
        raise ValueError('length of first and second arguments must be equal')


    for variable, mytype, label in zip(variables, types, labels):
        if type(mytype)!=list:
            mytype = [mytype]
        if type(variable) not in mytype:
            raise Exception(
                    """Expect %s argument's type to be in %s.
                    \nFound %s instead.""" % (label, mytype, type(variable)))


def CheckDictKeys(key, udict):
    """check that a key is contained in a dictionary and raise error otherwise.

    .. seealso:: ParseKargs
    """
    if key not in udict.keys():
        raise KeyError('Key %s not found. Possible choices are %s.' %
                            (key, udict.keys()))
    else:
        return udict[key]


def ParseKargs(kargs, keyword, default=None, possible=None):
    """Parse and check presence of a key in a dictionary

    Enhanced version of kargs.get

    :param kargs: a dictionary
    :param keyword: a key to look for
    :param default: value to assigned to keyword if keyword not found in kargs
    :param possible: values that kargs[keyword] can take (either a list or dict)

    :Example:

        >>> distance = error.ParseKargs(kargs, "Distance", "ABSOLUTE_VALUE",
                                distance_type.keys())
        >>> distance = distance_type[distance]

    or simply (if a dictionary is used as fourth argument)

        >>>distance = error.ParseKargs(kargs, "Distance", "ABSOLUTE_VALUE",
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
    """Check that a list of keywords are present in kargs

    :param kargs: dictionary containing a set of keys
    :param possible_kargs: a list containing a keywords
    :param dicts: list of dictionaries possible for each kargs


    :Usage:

        if kargs is {"key1", None, "key2": True, "key3": "dummy"}

        If a function requires "key1" and "key2" to work properly,
        we will use:

            CheckKargs(kargs, ["key1", "key2"])


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
    #if dicts:
    #    for karg, mydict in zip(possible_kargs, dicts):
    #        if kargs.get(karg) not in mydict.keys() and \
    #            kargs.get(karg) is not None:
    #            raise ValueError("""
    #                Value of the Argument \"%s\" (given %s) not found in
    #                the list of allowed values (%s)"""
    #                % (karg, kargs.get(karg), mydict.keys()))



STAT_TOOL_NB_VARIABLE_ERROR = \
    """Extra arguments provided (to specify variable value ?).
    Consider removing it. Be aware that nb_variable equals 1"""

