#!/usr/bin/env python
#-*- coding: utf-8 -*-
"""Functions dedicated to check function and class arguments

.. topic:: error module summary

    :Code status: mature
    :Documentation status: mature
    :Author: Thomas Cokelaer <Thomas.Cokelaer@sophia.inria.fr>
    :Revision: $Id: error.py 15344 2013-12-04 08:53:04Z jbdurand $

.. testsetup:: *

    from openalea.stat_tool.error import *
"""
__version__ = "$Id: error.py 15344 2013-12-04 08:53:04Z jbdurand $"



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

__all__ = ['CheckArgumentsLength',
           'CheckType',
           'CheckDictKeys',
           'ParseKargs',
           'CheckKargs',
           'FormatError']



def CheckType(variables, types, **kargs):
    """Check types of input list of variables

    .. warning:: only list are supported.

    :param variables: a list of variables to be checked
    :param types: list of types
    :param variable_pos: optional argument that provides the position of
        each variable. Used to enhance output message. For instance, if only
        a specific arguments (let us say the third one) has to be checked,
        use variable_pos=[3] and in case of errorm, the error message will be
        'the third argument is incorrect'

    .. todo:: consider removing the optional argument that is hardly used.

    :Examples:

        >>> #CheckType(1, int, variable_index=0) NOT IMPLEMENTED
        >>> CheckType([1, 'a'], [int, str], variable_pos=[0,1])
        >>> CheckType([1, 'a'], [[int, float], str], variable_pos=[0,1])

    """
    variable_pos = kargs.get("variable_pos", None)

    if type(variables)!=list:
        raise TypeError("first argument(variables) must be a list")

    if type(types)!=list:
        raise TypeError("second argument(variables) must be a list")

    if type(variable_pos)!=list and variable_pos is not None:
        raise TypeError("variable_pos must be a list")


    n = len(variables)
    #populate labels from 0 to N or with contents of variable_pos
    labels = []
    if variable_pos is None:
        for index in range(1, n+1):
            try:
                labels.append(arguments_labels[index])
            except:
                raise KeyError("""consider adding keys in the labels dictionary:
                 you've reached the limit""")
    else:
        for var in variable_pos:
            if var in arguments_labels.keys():
                labels.append(arguments_labels[var])
            else:
                labels.append("unknown position (larger than 10)")

    if len(variables) != len(types):
        raise ValueError('length of first and second arguments must be equal')


    for variable, mytype, label in zip(variables, types, labels):
        if type(mytype)!=list:
            mytype = [mytype]
        if type(variable) not in mytype:
            raise Exception(
                    """Expect %s argument's type to be in %s.
                    \nFound %s instead.""" % (label, mytype, type(variable)))

def CheckValue(variables, values, **kargs):
    """Check values of input list of variables

    .. warning:: only list are supported.

    :param variables: a list of variables to be checked
    :param values: list of possible values
    :param variable_pos: optional argument that provides the position of
        each variable. Used to enhance output message. For instance, if only
        a specific arguments (let us say the third one) has to be checked,
        use variable_pos=[3] and in case of errorm, the error message will be
        'the third argument is incorrect'

    .. todo:: consider removing the optional argument that is hardly used.

    :Examples:

        >>> CheckValue([1, 'a'], [[0,1], ['a']], variable_pos=[1,2])
        >>> CheckValue([1, 'a'], [1, ['a', 'b']], variable_pos=[1,2])

    """
    variable_pos = kargs.get("variable_pos", None)

    if type(variables)!=list:
        raise TypeError("first argument(variables) must be a list")

    if type(values)!=list:
        raise TypeError("second argument(variables) must be a list")

    if type(variable_pos)!=list and variable_pos is not None:
        raise TypeError("variable_pos must be a list")


    n = len(variables)
    #populate labels from 0 to N or with contents of variable_pos
    labels = []
    if variable_pos is None:
        for index in range(1, n+1):
            try:
                labels.append(arguments_labels[index])
            except:
                raise KeyError("""consider adding keys in the labels dictionary:
                 you've reached the limit""")
    else:
        for var in variable_pos:
            if var in arguments_labels.keys():
                labels.append(arguments_labels[var])
            else:
                labels.append("unknown position (larger than 10)")

    if len(variables) != len(values):
        raise ValueError('length of first and second arguments must be equal')


    for variable, myval, label in zip(variables, values, labels):
        if type(myval)!=list:
            myval = [myval]
        if variable not in myval:
            raise Exception(
                    """Expect %s argument's value to be in %s.
                    \nFound %s instead.""" % (label, myval, variable))

def CheckArgumentsLength(args, min_nargs=0, max_nargs=32):
    """Check that the number of arguments is valid

    This function check that the number of arguments of the list/tuple of
    arguments is in the given range.

    Used by functions to check the validity of the list of arguments
    (usually denoted `*args`) provided by the user.

    :param args: a tuple containing user arguments
    :param min: minimum number of arguments expected.
    :param max: maximum number of arguments expected (strict) (default is 32).

    :Example:

        >>> args = ('a','b')
        >>> CheckArgumentsLength(args, 1, 2)
    """
    CheckType([min_nargs], [[int, float]])
    l = len(args)
    assert max_nargs >= min_nargs, \
        "max_nargs must be greater or equal to min_args"
    assert max_nargs >= 0 and min_nargs >= 0, \
        "max_nargs and min_margs must be striclty positive"
    if (l < min_nargs) or (l > max_nargs):
        msg = """Expected at least %s arguments and at most %s arguments but
            %s were provided' """ % (min_nargs, max_nargs, l)
        raise Exception(msg)


def CheckDictKeys(key, udict):
    """check that a key is contained in a dictionary and raise error otherwise.

    :param key:
    :param udict: a valid dictionary
    :returns: the value corresponding to the key

    :Example:

        >>> d = {'a': [1,2], 'b':[1,2]}
        >>> res = CheckDictKeys('a', d)

    """
    CheckType([key], [str])
    if key not in udict.keys():
        raise KeyError('Key %s not found. Possible choices are %s.' %
                            (key, udict.keys()))
    else:
        return udict[key]


def ParseKargs(kargs, keyword, default=None, possible=None):
    """Utility to parse and check `**kargs` optional arguments

    This is an improved version of kargs.get() to be used after the function
    definition.

    :param kargs: a dictionary
    :param keyword: a key to look for
    :param default: value to assigned to keyword if not found in kargs
    :param possible: possible values that kargs[keyword] can take
        (either a list or dict)

    :Example:

        >>> kargs = {'a':[1,2], 'verbose':True}
        >>> a = ParseKargs(kargs, "a", [1,1])
        >>> verbose = ParseKargs(kargs, "verbose", False, possible=[False, True])

    The fourth argument may be a dictionary (values are irrelevant)

        >>> mykeys = {False:None, True:None}
        >>> distance = ParseKargs(kargs, "verbose", False, possible=mykeys)
    """
    CheckType([kargs], [dict])
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




def CheckKargs(kargs, possible_kargs):
    """Check that a list of keywords are present in kargs

    :param kargs: a dictionary such as `**kargs`
    :param possible_kargs: a list of possible keywords

    :Example:

        >>> d = {'a':[1,1], 'b':[1,2]}
        >>> CheckKargs(d, ['a', 'b'])
    """
    # check that number of arguments is correct

    CheckArgumentsLength(kargs, 0, len(possible_kargs))

    # check that argument names are correct
    for karg in kargs:
        if karg not in possible_kargs:
            raise KeyError('Argument %s not allowed (possible args are %s).',
                                karg, possible_kargs)


class FormatError(Exception):
         """Exceptions related to the statistical modules."""

         def __init__(self, error=None):
             """Initialize a FormatError exception."""
             if error is None:
                 self.__error=""
             else:
                 self.__error=error

         def _error(self):
             return str(self.__error)

         def __str__(self):
             return str(self.__error)
