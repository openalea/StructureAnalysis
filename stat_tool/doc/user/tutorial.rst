.. testsetup:: *
   
    from openalea.stat_tool import *


Tutorial
########


Histogram
=========

Let us suppose, we want to look at the file `meri5.his` in the `./test/` directory. It looks like::

    0    0
    1    0
    2    4
    3   14
    4   17
    5   18
    6   10
    7    3

We can read it as an histogram using the
:func:`~openalea.stat_tool.histogram.Histogram` function as follows:

.. filename with respect to the directory where sphinx is launch
.. doctest::

    >>> h = Histogram('./test/meri5.his')
   
The object `h` has a few methods, The `display` method :func:`~openalea.stat_tool.output.Display` returns information on the screen
   
.. doctest::

   >>> h.display()  

::

    histogram - sample size: 66
    mean: 4.37879   variance: 1.62354   standard deviation: 1.27418
    coefficient of skewness: 0.0727983   coefficient of kurtosis: -0.709664
    mean absolute deviation: 1.06841   coefficient of concentration: 0.161214
    information: -107.512 (-1.62897)


