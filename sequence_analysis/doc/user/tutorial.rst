.. testsetup:: *

    from openalea.stat_tool import *
    from openalea.sequence_analysis import *
    from openalea.sequence_analysis.data import path

Sequences tutorial (beginners)
##############################

.. topic:: Sequence tutorial on manipulating sequences

    Functions used:
        * data: Sequences, Vectors
        * output: Plot, Display
        * data exploration:

    :Authors: Thomas Cokelaer
    :Prerequisites: Sequences file syntax (see stat_tool package)

.. contents::



Import python package
======================


.. note:: preliminary code::

    >>> from openalea.stat_tool import *
    >>> from openalea.sequence_analysis import *
    >>> from openalea.sequence_analysis.data import path


Read some data
================

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

.. note:: min_length is not implemented but can be retrieved as follows::

    >>> min([len(s) for s in seq])

Then, you may not be interested in all the sequences, but only the two first one. This can be done 
using :func:`~openalea.stat_tool.data_transform.SelectIndividual` from the **stat_tool package**:

Select variable or individual
==============================

SelectIndividual
----------------
.. doctest::

    >>> s1 = SelectIndividual(seq, [1])

.. note:: All Functions have an object equivalent but there are usually more difficult to use
   (not type or bound checks)

The object equivalent works as follows.

.. doctest::

    >>> s1 = seq.select_individual([1], True)

The extracted sequences can now be displayed::

    >>> print Display(s1)

or introspect to check values:

.. doctest::

    >>> l = s1.get_length(0)
    >>> assert l==47


in order to access to the data (array of arrays) use indices as follows:

.. doctest::

    >>> s1[0][0]
    [0]

where the first index is the sequence number and the second one corresponds to the vector index.

.. note:: It is always possible to convert sequences into list or numpy array. Just be cautious with the indices.
    For instance to get the first sequence (here univariate)

        >>> import numpy
        >>> numpy.array(seq[0]).flatten()
        array([0, 0, 0, 0, 0, 3, 3, 0, 3, 3, 0, 0, 4, 1, 4, 4, 0, 0, 0, 2, 1, 0, 3, 0, 0, 0, 0, 2, 0, 3, 3, 0, 1, 4, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 4])


the `flatten` method allows to remove the list of list (for the univariate case this is quite useful).

SelectVariable
---------------

The example above provides univariate sequences so the :func:`openalea.stat_tool.data_transform.SelectVariable` is useless here but would work as follows::

    >>> SelectVariable(seq, 1)


Plotting
=========

Let us start to read the data again and use the  Plot function with a data ViewPoint::

    >>>    from openalea.sequence_analysis import *
    >>>    from openalea.sequence_analysis.data import path 
    >>>    seq = Sequences(path +  'sequences_tutorial.dat')
    >>>    Plot(seq, ViewPoint="Data")

.. plot::
    :width: 480px
    :height: 480px

    from openalea.sequence_analysis import *
    from openalea.sequence_analysis.data import path
    seq = Sequences(path +  'sequences_tutorial.dat')
    Plot(seq, ViewPoint="Data")

Then, we can look at an histogram of the values::

    >>> Plot(ExtractHistogram(seq, "Value"))

.. plot::
    :width: 480px
    :height: 480px

    from openalea.sequence_analysis import *
    from openalea.sequence_analysis.data import path
    seq = Sequences(path +  'sequences_tutorial.dat')
    Plot(ExtractHistogram(seq, "Value"))

and finally an histogram of the sequences length

    >>> Plot(ExtractHistogram(seq, "Length"))

.. plot::
    :width: 480px
    :height: 480px

    from openalea.sequence_analysis import *
    from openalea.sequence_analysis.data import path
    seq = Sequences(path +  'sequences_tutorial.dat')
    Plot(ExtractHistogram(seq, "Length"))






