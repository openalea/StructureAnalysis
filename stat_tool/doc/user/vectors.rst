.. define some aliases:
.. _vectors: syntax.html#type-vectors

.. define the setup for doctest:
.. testsetup:: *
   
    from openalea.stat_tool import *
    from pylab import savefig



Vectors
=======
Here is a brief description of the Vector type.

Constructor
-----------

Vectors can be generated either by loading an ascii file -- following the syntax given in vectors_ -- or with a list of values, that will be processed into a vectors.

So, if you have an ASCII file, you can load it using the :func:`~openalea.stat_tool.vectors.Vectors` function as follows:

.. filename with respect to the directory where sphinx is launch
.. doctest::

    >>> v1 = Vectors('./test/chene_sessile.vec')

Otherwise, if you have a list of values, provide it using python syntax as follows:

.. doctest::

    >>> v2 = Vectors([[1,2], [3,4]])

.. note:: Note the syntax, which is a list of lists

    >>> v2.get_nb_variable()
    2
    >>> v2.get_nb_vector()
    2
    >>> v2.get_identifiers()
    [1, 2]

You can then access to the vectors usin indexes::

    >>> v2[1]
    [1, 2]
    >>> v2[1][0]
    3

display
-------
The object `h` has a few methods, The `display` method :func:`~openalea.stat_tool.output.Display` returns information on the screen
   
.. doctest::

    >>> v2.display() #doctest: +SKIP
    >>> # equivalently
    >>> Display(v2)
    '2 vectors\n\n2 VARIABLES\n\nVARIABLE 1 : INT   (minimum value: 1, maximum value: 3)\n\nmarginal histogram - sample size: 2\nmean: 2   variance: 2   standard deviation: 1.41421\n\n   | marginal histogram\n0  0\n1  1\n2  0\n3  1\n\nVARIABLE 2 : INT   (minimum value: 2, maximum value: 4)\n\nmarginal histogram - sample size: 2\nmean: 3   variance: 2   standard deviation: 1.41421\n\n   | marginal histogram\n0  0\n1  0\n2  1\n3  0\n4  1\n\ncorrelation matrix\n\n   1  2\n1  1  1\n2  1  1\n\nreference t-value: 1e+37   reference critical probability: 0.05\nlimit correlation coefficient: 1\n\nreference t-value: 1e+37   reference critical probability: 0.01\nlimit correlation coefficient: 1\n'
    
A nicer layout can be obtained: by swithcing the \n character to a return carriage using the **print** command:

.. doctest::

    >>> print v2 #doctest: +SKIP
    2 vectors

    2 VARIABLES

    VARIABLE 1 : INT   (minimum value: 1, maximum value: 3)

    marginal histogram - sample size: 2
    mean: 2   variance: 2   standard deviation: 1.41421

       | marginal histogram
    0  0
    1  1
    2  0
    3  1

    VARIABLE 2 : INT   (minimum value: 2, maximum value: 4)
    ...

 
printing ASCII information (exhaustive output) on the screen or in a file:

.. doctest::

    >>> print v2.ascii_write(True) #doctest: +SKIP
    >>> print v2.file_ascii_write('output.dat', True) #doctest: +SKIP
    >>> print v2.save('output.dat') #doctest: +SKIP

The two last lines are equivalent.


save in a gnuplot file with plot_write method::

    >>> v1.plot_write('output', 'title')

clustering
-----------

.. doctest::
    :options: +SKIP

    >>> h1.cluster_information()
    >>> h1.cluster_limit([1,2])
    >>> h1.cluster_step()


