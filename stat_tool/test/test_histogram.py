from openalea.stat_tool.histogram import *

def test_fromnothing():

    # From noting
    try:
        h = Histogram()
        assert False

    except TypeError:
        assert True

def test_fromfile():

    # From file
    
    h = Histogram("peup1.his")
    assert h
    
    try:
        h = Histogram("peup1.hi")
        assert False

    except StatToolError:
        assert True

def test_fromlist():

    # From list
    
    h = Histogram([1,2,3])
    assert h

    try:
        h = Histogram([0, 'a', 'b', 3])
        assert False

    except TypeError:
        assert True

    
