__doc__ = """ Histogram functions and classes """
__docformat__ = "restructuredtext"

import sys,os
sys.path.append(os.path.abspath("."))

import _stat_tool
import interface

from _stat_tool import _DistributionData

# Extend _DistributionData class dynamically
interface.extend_class( _stat_tool._DistributionData, interface.StatInterface)

__all__ = ["_DistributionData",
           "Histogram",
           ]


def Histogram(arg):
    """
    Construction of a frequency distribution from an object of type list(int) or from
    an ascii file.
    
    Usage
    -----
      * ``Histogram(list)``
      * ``Histogram(filename)``

    Parameters
    ----------
      * list (list(int))
      * filename (string)

    Return
    ------
       If the construction succeeds, an object of type `_DistributionData` is returned, 
    
    Description
    -----------
    In the file syntax, the frequencies fi for each possible value i are given
    in two columns. In the case of an argument of type (list(int)), it is simply assumed 
    that each array element represents one data item.

    See Also
    --------
       `Save`, `Cluster`, `Merge`, `Shift`, `Transcode`, `ValueSelect`, `Compare`, `Estimate`

    """

    return _DistributionData(arg)








############################## Tests ###########################################
from openalea.stat_tool import get_test_file

# Test Histogram
class Test:

    def test_constructor(self):
    
        h = Histogram([0,1,2,3])
        assert h


    def test_ascii(self):
        # ASCII representation
        h = Histogram(get_test_file("meri1.his"))
        assert h

        s = str(h)
        assert h.display() == s

        h.ascii_write(True)
        h.survival_ascii_write()
    


    def __test_plot(self):
        # plot
        h = Histogram(get_test_file("meri1.his"))
        h.plot()


    def test_container(self):
        # container / iterator
        h = Histogram(get_test_file("meri1.his"))

        assert h[0] == 0
        assert h[10] == 1


    def test_fromnothing(self):

        # From nothing
        try:
            h = Histogram()
            assert False

        except TypeError:
            assert True


    def test_fromfile(self):

        # From file
        h = Histogram(get_test_file("meri1.his"))
        assert h
    
        # File
        h.save("test.his")

        h1 = Histogram(get_test_file("meri1.his"))
        h2 = Histogram("test.his")
        assert len(h) == len(h2)
        assert list(h) == list(h2)
        os.remove("test.his")
    
    
        try:
            h = Histogram(get_test_file("peup1.hi"))
            assert False
        except Exception:
            assert True
        

    def test_len(self):
        h = Histogram(range(10))
        assert len(h) == 10

    
    def test_fromlist(self):

        # From list
        h = Histogram([1,2,3])
        assert h

        try:
            h = Histogram([0, 'a', 'b', 3])
            assert False

        except TypeError:
            assert True


	
