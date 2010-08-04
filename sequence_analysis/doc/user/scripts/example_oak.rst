
Sequences advanced tutorial: Oak trunk annual shoots
=====================================================

.. topic:: Illustrative example of stat_tool and sequence_analysis packages

    Functions used:
        * data: Sequences, Vectors
        * output: Plot, Display
        * data exploration: ExtractHistogram, Cluster, Regression, ValueSelect, VarianceAnalysis, Compare, ContingencyTable, Estimate

    :Authors: Thomas Cokelaer
    :References: Oak trunk annual shoot AML example :download:`example_oak.py`
    :Data: Patrick Heuret

.. contents::

:Material:

  * data :download:`../../../share/data/chene_sessile_15pa.seq`
  * whole script :download:`example_oak.py`

In this example, we will explain some functionalities of sequence_analysis package by going through the python script `example_oak.py` step by step. First, we will start to explain what is a sequence data file and how to explore and visualize its contents.

.. _preliminaries:

Preliminaries: the data set
----------------------------

The file `chene_sessile_15pa.seq` is used as a data set. Its
extension `seq` indicates that this file contains a list of sequences. 

This file contains an header that provides the number of variables and their
types. Here there are 6 variables that are labelled as **STATE** variables. 
There are 46 sequences in total, and a sequence looks like::

    95 110 219 2 52 14 | 96 17 119 2 24 9 | 97 57 101 2 33 1

.. note:: Each element in the sequence is separated by the other by the pipe
    symbol, which is optional when the sequence is univariate. In the sequence
    above, we have 6 variables and each sequence is 3-elements long (one for
    each year 95, 96, 97).

The Sequences data structures allows to read a file and returns a sequences object that can be indexed like a normal list. The first index being the sequence index and the
second index being the vector index.

::

    >>> seq = Sequences('chene_sessile_15pa.seq')
    >>> seq[0]
    [[95, 110, 219, 2, 52, 14], [96, 17, 119, 2, 24, 9], [97, 57, 101, 2, 33, 1]]
    >>> seq[0][1]
    [96, 17, 119, 2, 24, 9]
    >>> seq[0,1]
    [96, 17, 119, 2, 24, 9]

.. note:: Note the syntax to access the second vector of the first sequence. First, 
   indices starts at zero like in C/python languages and the index syntax may be either   [i][j] or [i,j].


You may want to extract the vectors one by one as follows::

    >>> seq.nb_sequence
    46
    >>> vec = Vectors(seq)
    >>> vec.nb_vectors
    138
    >>> assert vec[3*10] == seq[10][0]

138 is simply the number of sequences times the number of element per sequences (3 here).

Read sequences in a file and plot a data viewpoint
---------------------------------------------------

Let us come back to the example. The first step consists in reading the data. We need to use use the :func:`~openalea.sequence_analysis.sequences.Sequences` function that will check the header and return a sequence object. Then, we can plot the results using a Data viewpoint (see :func:`openalea.stat_tool.output.Plot`).

.. literalinclude:: example_oak.py
    :lines: 1-26
    :linenos:

.. plot:: pyplots/example_oak_1.py
    :width: 480px
    :height: 480px

Looking at one variable
-----------------------
On the first line of the following code, we first extract a given variable (the third one) and then cluster it before plotting the histogram. Note that the extraction is made in such a way that each vector in each sequence is extracted in a single object. So you have 138 vectors in the **marginal3** variable.



.. literalinclude:: example_oak.py
    :lines: 27-32
    :linenos:

.. plot:: pyplots/example_oak_2.py
    :width: 480px
    :height: 480px

Regression and vectors
------------------------

Here below, we first get a list of vectors in place of the sequence (vectors and sequences are briefly cover in the :ref:`preliminaries` section.

The variable marginal3 could have been obtained as follows::


    >>> vec10 = Vectors(seq0)
    >>> marginal3 = ExtractHistogram(vec10, 3)

Now, coming back to the first (year) and second variable (growth length), we plot the average of the second variable versus the year as follows

.. literalinclude:: example_oak.py
    :lines: 35-39
    :linenos:

.. plot:: pyplots/example_oak_3.py
    :width: 480px
    :height: 480px

Comparison and variance analysis
-----------------------------------

Then, we want to look those data in more details by selecting the length year by year. So, we first select each year one by one (line 1 to 3 below). Here, :func:`openalea.stat_tool.data_transform.ValueSelect` selects the first variable when its values is 95, or 96 or 97 from which an histogram is extracted. All 3 histograms can be plotted together (line 11).


.. literalinclude:: example_oak.py
    :lines: 40-52
    :linenos:

.. plot:: pyplots/example_oak_4.py
    :width: 480px
    :height: 480px

The output of the variance analysis :func:`openalea.stat_tool.vectors.VarianceAnalysis` (line 5) is ::


    value                                 95          96          97
    sample size                           46          46          46
    mean                             58.4348     24.2174     45.7609
    variance                         438.873     88.5739     501.875
    standard deviation               20.9493     9.41137     22.4026
    mean absolute deviation          16.3724     7.41588     19.5227
    coefficient of concentration     0.19955    0.215518    0.277259
    coefficient of skewness       -0.0383523   0.0248062    0.155391
    coefficient of kurtosis        -0.219122  -0.0256427    -1.10015
    
    source of variation | degrees of freedom | sum of squares | mean square
    between samples    2  27532.2  13766.1
    within samples   135  46319.5  343.107
    total            137  73851.7  539.064
    
    F-test (2 degrees of freedom, 135 degrees of freedom)
    F-value: 40.1219   critical probability: 2.69378e-12
    reference F-value: 4.11997   reference critical probability: 0.05
    reference F-value: 6.27039   reference critical probability: 0.01
    

Comparison and ContingencyTable
-----------------------------------------------------

Here again, we will compare variables but looking at the ranks rather than the years.

We first look at the contingency table ::

    >>> ContingencyTable(vec, 1, 4)

::

    contingency table
    
          1    2    3    4
    95    1   41    3    1   46
    96   38    7    1    0   46
    97   18   27    1    0   46
         57   75    5    1  138
    
    deviation table
    
                1          2          3          4
    95        -18         16    1.33333   0.666667
    96         19        -18  -0.666667  -0.333333
    97         -1          2  -0.666667  -0.333333
    
    chi-square contribution table
    
                  1            2            3            4
    95     0.270397     0.162371    0.0169137    0.0211421
    96     0.301275     0.205501   0.00422842   0.00528553
    97  0.000834557   0.00253705   0.00422842   0.00528553
    
    chi-square test (6 degrees of freedom)
    chi-square value: 63.0653   critical probability: 7.9565e-11
    reference chi-square value: 12.5708   reference critical probability: 0.05
    reference chi-square value: 16.7967   reference critical probability: 0.01

and then plot some histograms.

.. literalinclude:: example_oak.py
    :lines: 52-63
    :linenos:

.. plot:: pyplots/example_oak_5.py
    :width: 480px
    :height: 480px



Estimate
-----------

Now, we extract an historam from the second variable and estimated a model that is a mixture of 4 negatibe binomial
distribution and plot the results either the mixture model or component by component

.. literalinclude:: example_oak.py
    :lines: 71-75
    :linenos:

.. plot:: pyplots/example_oak_6.py
    :width: 480px
    :height: 480px


Regression
-------------

.. literalinclude:: example_oak.py
    :lines: 80-83
    :linenos:

.. plot:: pyplots/example_oak_7.py
    :width: 480px
    :height: 480px









