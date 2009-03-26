.. testsetup:: *
   
    from openalea.stat_tool import *


Tutorial
########


Histogram
=========

Loading a file
--------------
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

display the content
-------------------   
The object `h` has a few methods, The `display` method :func:`~openalea.stat_tool.output.Display` returns information on the screen
   
.. doctest::

    >>> h.display()  
    'histogram - sample size: 66\nmean: 4.37879   variance: 1.62354   standard deviation: 1.27418\ncoefficient of skewness: 0.0727983   coefficient of kurtosis: -0.709664\nmean absolute deviation: 1.06841   coefficient of concentration: 0.161214\ninformation: -107.512 (-1.62897)\n'

Saving
------

.. doctest::

    h.save('whatever_name.his') #doctest: +SKIP


Vectors
=======

Similarly, you can upload a vector using::

.. doctest::
    :options: +SKIP

    >>> h = Vectors('./test/chene_sessile.vec')
    >>> h.display() #doctest: +SKIP



