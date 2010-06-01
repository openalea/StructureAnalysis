.. testsetup:: *

    from openalea.stat_tool import *
    from openalea.sequence_analysis import *
    from openalea.sequence_analysis.data import path


Sequences tutorial
##################

.. note:: preliminary code::

    >>> from openalea.stat_tool import *
    >>> from openalea.sequence_analysis import *
    >>> from openalea.sequence_analysis.data import path



Let us start from the following file data, which is fulfills the sequences syntax (see 
the stat_tool user guide documentation for the syntax).

::

    1 VARIABLE

    VARIABLE 1 : VALUE

    0 0 0 0 0 3 3 0 3 3 0 0 4 1 4 4 0 0 0 2 1 0 3 0 0 0 0 2 0 3 3 0 1 4 3 0 0 0 0 0 0 0 0 0 4 0 4
    0 0 0 0 0 0 0 0 0 4 0 0 0 1 0 3 0 1 0 0 0 0 0 1 0
    0 0 0 0 0 0 0 0 0 3 0 0 0 0 0 0 0 0 4 1 3 4 4 0 0 0 4 4 0 3 0 0 0 0 1 0 2 0 4 4 0 0 0 0 4 0 4 4 0 4 4 0 4 4 0 4 0 0 0 0 0
    0 0 0 0 0 0 0 0 3 0 0 0 0 0 0 3 3 0 0 4 0 4 0 0 0 0 0 4 0 0 0 4 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 4 1 0 3 3 0 4 0 4 0 1 3 2 0 0
    0 0 0 0 0 0 0 0 4 4 0 4 0 3 0 0 0 4 0 0 0 0 0 4 0 1 0 0 0 0 0 0 2 1 1 4 2 0 0 0 0 0 0 0 4 4 0 0 4 0 4 0 0 0 0 4 4
    0 0 0 0 0 0 0 0 0 0 0 3 0 2 0 0 0 0 0 4 0 0 0 0 0 1 0 0 1 0 0 1 0 1 3 3 3 4 0 2 0 2 3 0 0 0 0 0
    0 0 3 0 0 0
    0 0 0 0

This file can be read using :func:`~openalea.sequence_analysis.sequences.Sequences` as follows

.. doctest::

    >>> seq = Sequences(path +  'sequences_tutorial.dat')


Now, you can start to introspect the object `seq`. 
For instance, you can obtain the number of variables, the maximum length among the sequences or the number of elements over all sequences:

.. doctest::

    >>> seq.nb_variable
    1
    >>> seq.nb_sequence
    8
    >>> seq.cumul_length
    310
    >>> seq.max_length
    62


Then, you may not be interested in all the sequences, but only the two first one. This can be done 
using :func:`~openalea.stat_tool.data_transform.SelectIndividual` from the **stat_tool package**:



.. doctest::

    >>> s1 = SelectIndividual(seq, [1])

.. note:: All Functions have an object equivalent but there are usually more difficult to use
   (not type or bound checks)

The object equivalent works as follows.

.. doctest::

    >>> s1 = seq.select_individual([1], True)

The extracted sequences can now be plotted:

.. doctest::

    >>> Plot(s1)

or displayed::

    >>> print Display(s1)


or introspect:

.. doctest::

    >>> s1.get_length(0)
    47


in order to access to the data (array of arrays):

.. doctest::

    >>> s1[0][0]
    [0]





.. plot:: pyplots/test.py
