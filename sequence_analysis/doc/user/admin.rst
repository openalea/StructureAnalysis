Current Developments
####################


.. topic:: This page is not part of the documentation

    This is a worksheet to gather various information about current of future developments.

.. contents::



Test
====

There are many tests in sequence analysis package. They can be run using nosetests
to check the validity of the tests in the ./test directory. Another solution
is to check the test in the module's docstring if there are not flagged with 
#doctest:+SKIP keyword.


Using nosetests as follows::

    python setup.py nosetest

as of end of May 2010, the test coverage (with 921 tests) was:


=========================================================== ====== ======= ========
File                                                        Missed Skipped Percent
=========================================================== ====== ======= ========
openalea.sequence_analysis                                  19      19     100%
openalea.sequence_analysis.compare                          108     93      86%
openalea.sequence_analysis.correlation                       72     72     100%
openalea.sequence_analysis.data_transform                   429    317      73%
openalea.sequence_analysis.enums                             62     62     100%
openalea.sequence_analysis.estimate                         286    224      78%
openalea.sequence_analysis.hidden_semi_markov                20     20     100%
openalea.sequence_analysis.hidden_variable_order_markov      14     14     100%
openalea.sequence_analysis.nonhomogeneous_markov             15     15     100%
openalea.sequence_analysis.renewal                           68     68     100%
openalea.sequence_analysis.semi_markov                       15     15     100%
openalea.sequence_analysis.sequences                        165    139      84%
openalea.sequence_analysis.simulate                          46     43      93%
openalea.sequence_analysis.time_events                       42     34      80%
openalea.sequence_analysis.top_parameters                    33     33     100%
openalea.sequence_analysis.tops                              32     20      62%
openalea.sequence_analysis.variable_order_markov             18     18     100%
TOTAL                                                       1444    1206    83%
=========================================================== ====== ======= ========

In addition more than 10 tests are available using doctest. Use the following command to run those tests::

    python setup.py build_sphinx -b doctest








