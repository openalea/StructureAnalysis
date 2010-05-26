Sequences
#########

Let us start from the following file that follows the sequences syntax (see stat_tool user guide documentation for the syntax).

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

you read the file as follows

.. code-block:: python

    >>> from openalea.stat_tool import *
    >>> from openalea.sequence_analysis import *
    >>> seq = Sequences('sequences.dat')
    

Then, you can introspect the object `seq`. For instance, you can obtain the number of variables, the maximum length among the sequences os the number of elements over all sequences::

    >>> seq.nb_variable
    1
    >>> seq.nb_sequence
    8
    >>> seq.cumul_length
    310
    >>> seq.max_length
    62


Then, you may not be interested in all the sequences, but only the two first one. Use `SelectIndividual` ::

    >>> s1 = SelectIndividual(seq, [1])

which is equivalent to the object syntax (which is more complicated since you have to provide a second argument)::
    
    >>> s1 = seq.select_individual([1], True)

so, in the following we only use the AML syntax with functions except if the function does not exist.

The extracted sequences can now be plot::

    >>> s1.plot()

or displayed::

    >>> print s1.display()
    

or introspect::

    >>> s1.get_length()
    47


in order to access to the data (array of arrays)::

    >>> data = s1[0][0]



.. plot:: pyplots/test.py
