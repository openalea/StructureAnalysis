.. testsetup:: *

    from openalea.stat_tool import *
    from openalea.sequence_analysis import *
    from openalea.sequence_analysis import shared_data_path as path
    from os.path import join as pj

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

Introduction
============

In this tutorial you will learn how to read Sequences stored in a file. At the end of this tutorial you should be able
to manipulate Sequences so as to extract relevant information. In addition, you should be able to display and plot 
Sequences.


In order to go through this tutorial, you will first need to import the stat_tool ans sequence_analysis packages

.. code-block:: python

    >>> from openalea.stat_tool import *
    >>> from openalea.sequence_analysis import *

Then, we will need some data to play with. We will use some shared data available within the sequence_analysis package by typing::

    >>> from openalea.sequence_analysis import shared_data_path as path
    >>> from os.path import join as pj

Read some data
================

We will use data that are stored as :func:`Sequences` (see the stat_tool user guide documentation for the syntax). The file that is used within this tutorial contains univariate sequences with one sequence per line.

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

    >>> seq = Sequences(pj(path ,  'sequences_tutorial.dat'))


Now, you can start to introspect the object `seq`. For instance, you can obtain the number of variables, the maximum length among the sequences or the number of elements over all sequences as follows

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
   (not type or bound checks) but can be used to introspect all possible methods.

The object equivalent works as follows.

.. doctest::

    >>> s2 = seq.select_individual([1], True)
    >>> assert s1[0]==s2[0]

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
        >>> numpy.array(seq[0]).flatten() #doctest: +SKIP
        array([0, 0, 0, 0, 0, 3, 3, 0, 3, 3, 0, 0, 4, 1, 4, 4, 0, 0, 0, 2, 1, 0, 3, 0, 0, 0, 0, 2, 0, 3, 3, 0, 1, 4, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 4])


the `flatten` method allows to remove the list of list (for the univariate case this is quite useful).

SelectVariable
---------------

The example above provides univariate sequences so the :func:`openalea.stat_tool.data_transform.SelectVariable` is useless here but would work as follows::

    >>> SelectVariable(seq, 1)

.. warning:: to be back compatible with old AML syntax, variable index starts at 1 (not zero).

Plotting
=========

Data viewpoint
----------------

If you  know what you are doing, you may use your favorite plotter and extract data 
by hand. However, stat_tool and sequence analysis provides
dedicated viewpoints to each type of data structures, which majes life easier. These viewpoints derive from the 
:func:`~openalea.stat_tool.output.Plot` function. By default the Plot function is quite verbose and will plot curves 
for each sequence and each variable. So, we should use a specific viewpoint. The "Data" viewpoint is a good starting point.

::

    >>>    Plot(seq, ViewPoint="Data")

.. plot::
    :width: 480px
    :height: 480px

    from openalea.sequence_analysis import *
    from openalea.sequence_analysis import shared_data_path
    from os.path import join as pj
    seq = Sequences(pj(shared_data_path, 'sequences_tutorial.dat'))
    Plot(seq, ViewPoint="Data")

As you can see, all sequences are plotted. You may select only a subsets and plot them as follows::

    >>>    Plot(SelectIndividual(seq, [1,2]), ViewPoint="Data")

.. plot::
    :width: 480px
    :height: 480px

    from openalea.sequence_analysis import *
    from openalea.sequence_analysis import shared_data_path
    from os.path import join as pj
    seq = Sequences(pj(shared_data_path, 'sequences_tutorial.dat'))
    Plot(SelectIndividual(seq, [1,2]), ViewPoint="Data")

Normal Viewpoint
-----------------

Normal viewpoint is done as follows::

    >>> Plot(seq)

We let the user look at the output by himself. You will see a lot of outputs specific to markovian sequences (Observation, sojourn time, ...). 

Now, if you do::
    
    >>> Plot(SelectIndividual(seq, [1,2]))

You won't get at all the same kind of output. The reason being that SelectIndividual does not return a Markovian Sequence. 

Sojourn and Markovian sequences specific viewpoints
----------------------------------------------------
    
    >>>    Plot(seq, "Sojourn")

.. plot::
    :width: 480px
    :height: 480px

    from openalea.sequence_analysis import *
    from openalea.sequence_analysis import shared_data_path
    from os.path import join as pj
    seq = Sequences(pj(shared_data_path, 'sequences_tutorial.dat'))
    Plot(seq,  "Sojourn")



Others
------

Then, we can look at an histogram of the values::

    >>> Plot(ExtractHistogram(seq, "Value"))

.. plot::
    :width: 480px
    :height: 480px

    from openalea.sequence_analysis import *
    from openalea.sequence_analysis import shared_data_path
    from os.path import join as pj
    seq = Sequences(pj(shared_data_path,  'sequences_tutorial.dat'))
    Plot(ExtractHistogram(seq, "Value"))

and finally an histogram of the sequences length::

    >>> Plot(ExtractHistogram(seq, "Length"))

.. plot::
    :width: 480px
    :height: 480px

    from openalea.sequence_analysis import *
    from openalea.sequence_analysis import shared_data_path
    from os.path import join as pj
    seq = Sequences(pj(shared_data_path,  'sequences_tutorial.dat'))
    Plot(ExtractHistogram(seq, "Length"))






