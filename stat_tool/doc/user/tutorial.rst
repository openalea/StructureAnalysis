.. _stat_tool_tutorial:

.. define some aliases:
.. _histogram: syntax.html#type-histogram

.. define the setup for doctest:
.. testsetup:: *
   
    from openalea.stat_tool import *
    from pylab import savefig

.. contents::


Tutorial
########



.. .. include:: histogram.rst
.. .. include:: vectors.rst

.. toctree::
    :maxdepth: 1

    histogram.rst
    vectors.rst


Distribution, Histogram, Compound
=================================

.. doctest:: 

    # reads a set of histograms
    meri1 = Histogram("meri1.his")
    meri2 = Histogram("meri2.his")
    meri3 = Histogram("meri3.his")
    meri4 = Histogram("meri4.his")
    meri5 = Histogram("meri5.his")

    # misx them 
    mixt1 = Mixture("mixture1.mixt") #doctest: +SKIP
    mixt1 = Mixture(0.6, Distribution("B", 2, 18, 0.5), 0.4, Distribution("NB", 10, 10, 0.5)) #doctest: +SKIP

    # simulate
    mixt_histo1 = Simulate(mixt1, 200) #doctest: +SKIP

    mixt_histo1.plot() #doctest: +SKIP



.. doctest::

    >>> h = Histogram("meri2.his") #doctest: +SKIP 
    >>> m = h.estimate_mixture(["B", "NB"]) #doctest: +SKIP 
    >>> d = m.extract_data() #doctest: +SKIP
    >>> d.plot(show=False) #doctest: +SKIP
    >>> savefig('doc/user/tutorial_fig1.png')
    

.. figure:: tutorial_fig1.png

