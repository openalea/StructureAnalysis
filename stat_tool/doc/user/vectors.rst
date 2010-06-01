.. define some aliases:
.. _vectors: syntax.html#type-vectors
.. .. _histogram_tutorial: tutorial.html#histogram
.. _convolution_tutorial: tutorial.html#convolution

.. define the setup for doctest:
.. testsetup:: *
   
    from openalea.stat_tool import *
    import pylab
    from pylab import savefig




Vectors
=======
Let us finish with a third type of objects, that are vectors, which contains 
slighty more methods in addition to those already presented in the 
`histogram_tutorial_ <test>` and `Convolution <convolution_tutorial_>` sections

Constructor
-----------

Again, :class:`~openalea.stat_tool.vectors.Vectors` can be generated either by 
loading an ASCII or directly using python lists as follows:

.. filename with respect to the directory where sphinx is launch
.. doctest::

    >>> v1 = Vectors('./test/data/chene_sessile.vec') #doctest: +SKIP
    >>> v2 = Vectors([[1,2], [3,4]])

.. note:: Note the syntax, which is a list of lists

Then, you can access to various information using:

.. doctest::

    >>> v2.nb_variable
    2
    >>> v2.nb_vector
    2
    >>> v2.get_identifiers()
    [1, 2]

Finally, container are available and you can access to the data as follows 
(starting at 0):

    >>> v2[1]
    [3, 4]
    >>> v2[1][0]
    3

Display, Save, str()  methods are available as in the previous cases. 

However, there is no plotting routines available.
    
There are many more methods available, some of them are explained here below

VarianceAnalysis
----------------

Here is the usage of One-way variance analysis.
    
.. doctest::

    >>> print VarianceAnalysis(v2, 1,2,"O")
    value                          1   3
    sample size                    1   1
    mean                           2   4
    variance                       0   0
    standard deviation             0   0
    mean absolute deviation        0   0
    coefficient of concentration   1   1
    coefficient of skewness        0   0
    coefficient of kurtosis       -2  -2
    <BLANKLINE>
       | frequency distribution 1 | frequency distribution 3 | cumulative distribution 1 function | cumulative distribution 3 function
    0  0  0  0  0
    1  0  0  0  0
    2  1  0  1  0
    3     0     0
    4     1     1
    <BLANKLINE>
    Kruskal-Wallis test
    chi-square test (1 degree of freedom)
    chi-square value: 1   critical probability: 0.315013
    reference chi-square value: 3.74866   reference critical probability: 0.05
    <BLANKLINE>



Compare
-------

.. doctest::



    >>> print Compare(ExtractHistogram(v2, 1), ExtractHistogram(v2,2), "O")
    frequency distribution 1 - sample size: 2
    mean: 2   variance: 2   standard deviation: 1.41421
    coefficient of skewness: 0   coefficient of kurtosis: -2.5
    mean absolute deviation: 1   coefficient of concentration: 0.25
    information: -1.38629 (-0.693147)
    <BLANKLINE>
    frequency distribution 2 - sample size: 2
    mean: 3   variance: 2   standard deviation: 1.41421
    coefficient of skewness: 0   coefficient of kurtosis: -2.5
    mean absolute deviation: 1   coefficient of concentration: 0.166667
    information: -1.38629 (-0.693147)
    <BLANKLINE>
       | frequency distribution 1 | frequency distribution 2 | cumulative distribution 1 function | cumulative distribution 2 function
    0  0  0    0    0
    1  1  0  0.5    0
    2  0  1  0.5  0.5
    3  1  0    1  0.5
    4     1         1
    <BLANKLINE>
    dissimilarities between frequency distributions
    <BLANKLINE>
                | frequency distribution 1 | frequency distribution 2
    frequency distribution 1      0   0.5
    frequency distribution 2   -0.5     0
    <BLANKLINE>
    Kruskal-Wallis test
    chi-square test (1 degree of freedom)
    chi-square value: 0.6   critical probability: 0.448429
    reference chi-square value: 3.74866   reference critical probability: 0.05
    <BLANKLINE>

ContingencyTable
----------------

.. doctest::

    >>> print ContingencyTable(v2, 1, 2)
    contingency table
    <BLANKLINE>
       2  3  4
    1  1  0  0  1
    2  0  0  0  0
    3  0  0  1  1
       1  0  1  2
    <BLANKLINE>
    deviation table
    <BLANKLINE>
          2     3     4
    1   0.5     0  -0.5
    2     0     0     0
    3  -0.5     0   0.5
    <BLANKLINE>
    chi-square contribution table
    <BLANKLINE>
          2     3     4
    1  0.25     0  0.25
    2     0     0     0
    3  0.25     0  0.25
    <BLANKLINE>
    chi-square test (1 degree of freedom)
    chi-square value: 2   critical probability: 0.160475
    reference chi-square value: 3.74866   reference critical probability: 0.05
    <BLANKLINE>
             
