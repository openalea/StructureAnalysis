from openalea.stat_tool.error import *






def test_CheckArgumentsLength():
    args = ('a', 'b', 'c')
    CheckArgumentsLength(args, 1, 3)

    # try wrong values of max_nargs
    try:
        CheckArgumentsLength(args, 1, -1)
        assert False
    except:
        assert True

    # try wrong values of min_nargs
    try:
        CheckArgumentsLength(args, -1)
        assert False
    except:
        assert True

    # try wrong values of max_nargs
    try:
        CheckArgumentsLength(args, 3, 1)
        assert False
    except:
        assert True
    
    # try wrong values of max_nargs
    try:
        CheckArgumentsLength(args, 4, 7)
        assert False
    except:
        assert True

def test_CheckType():

    CheckType([1], [int])
    CheckType([1], [[int,float]])
    CheckType([1], [[int,float]], variable_id=[1])

    #first argument must be a list
    try:
        CheckType(1, [int])
        assert False
    except:
        assert True

    # arguments must be a list
    try:
        CheckType(1, int)
        assert False
    except:
        assert True

    # second argument ust be a list
    try:
        CheckType([1], int)
        assert False
    except:
        assert True
   
    # check that length are equal 
    try:
        CheckType([1], [int,int])
        assert False
    except:
        assert True
    
    #check that third arg is a list
    try:
        CheckType([1], [int], variable_id=1)
        assert False
    except:
        assert True

    #check that list length is less than 10
    try:
        CheckType([1,2,3,4,5,6,7,8,9,10,11], [int,int,int,int,int,int,int,int,int,int,int])
        assert False
    except:
        assert True
    
    #variable_pos too large (greater than 10)
    try:
        CheckType([12], [int], variable_pos=[12])
        assert False
    except:
        assert True
    #variable_pos type is wrong
    try:
        CheckType([12], [int], variable_pos=12)
        assert False
    except:
        assert True


def test_CheckDictKeys():

    CheckDictKeys('a', {'a':[1,2]})==[1,2]
    try:
        CheckDictKeys('b', {'a':[1,2]})
        assert False
    except:
        assert True


def test_ParseKargs():
    kargs = {'a':[1,2], 'verbose':True}
    ParseKargs(kargs, 'a', [1,1])
    ParseKargs(kargs, 'a', [1,1], possible=[[1,1], [1,2]])

    keys = {False:False, True:True}
    ParseKargs(kargs, 'verbose', False, possible=keys)

    try:
        ParseKargs('a', 'a')
        assert False
    except:
        assert True
    #possible keys from dictionary are not found
    try:
        ParseKargs(kargs, 'a', possible={'dummy':None})
        assert False
    except:
        assert True

    #possible keys from list are not found
    try:
        ParseKargs(kargs, 'a', possible=['dummy'])
        assert False
    except:
        assert True



def test_CheckKargs():

    d = {'a':[1,1], 'b':[1,2]}
    CheckKargs(d, ['a', 'b'])

    #incorrect input dictionary that do bot contain b
    try:
        CheckKargs(d, ['a'])
        assert False
    except:
        assert True
    try:
        CheckKargs(d, ['a','d'])
        assert False
    except:
        assert True


def test_class_error():

    fe = FormatError()
    fe._error()
    fe = FormatError(error="test")
    fe._error()
    print fe
