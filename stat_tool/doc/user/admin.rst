
Current Developments
####################


.. topic:: This page is not part of the documentation

    This is a wokrsheet to gather various information about current of future development.

.. contents::

.. methodology (in french)



.. include:: ../../methodo.txt



TODO
====

* Check that tests run under windows
* Issue with the Histogram distribution in test_save:
  if h = Histogram('mixture1.mixt/')
  len(h) returns 76 but this seem to be the length of the original data set,
  not the histogram itself. Is this what we want ? 
* plotable inside MvMixture(Data) ? 
* consistency of property names: nb_vector or nb_vectors, nb_variable, nb_sequence, ...

Test
====

Using nosetests as follows::

    python setup.py nosetest

as of end of May 2010, the test coverage (with 891 tests) was:


======================================= ======= ====== ======= ========
File                                    Covered Missed Skipped Percent
======================================= ======= ====== ======= ========
openalea.stat_tool.plot                 51      175    121     22 %
openalea.stat_tool.multivariate_mixture 45      53     38      45 %
openalea.stat_tool.output               514     315    203     62 %
openalea.stat_tool.estimate             248     52     131     82 %
openalea.stat_tool.cluster              317     37     121     89 %
openalea.stat_tool.data_transform       605     52     210     92 %
openalea.stat_tool.regression           93      7      49      93 %
openalea.stat_tool.vectors              263     1      143     99 %
openalea.stat_tool                      24      0      25      100 %
openalea.stat_tool.comparison           174     0      71      100 %
openalea.stat_tool.compound             61      0      35      100 %
openalea.stat_tool.convolution          51      0      30      100 %
openalea.stat_tool.distribution         189     0      78      100 %
openalea.stat_tool.enums                170     0      54      100 %
openalea.stat_tool.error                166     0      69      100 %
openalea.stat_tool.histogram            47      0      28      100 %
openalea.stat_tool.interface            17      0      10      100 %
openalea.stat_tool.mixture              61      0      32      100 %
openalea.stat_tool.simulate             44      0      15      100 %
Total                                   3140    692    1463    81%
======================================= ======= ====== ======= ========

In addition more than 130 tests are available using doctest using the following command::

    python setup.py build_sphinx -b doctest


