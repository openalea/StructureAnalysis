.. _stat_tool_tutorial:

.. define the setup for doctest:
.. testsetup:: *

    from openalea.stat_tool import *
    from pylab import savefig


Tutorial
########

.. contents::

.. note:: most of the docstrings provided in this tutorial are tested with
    doctest and should work out of the box. However, you will need to import
    relevant modules.  For instance:

    >>> from openalea.stat_tool import *
    >>> from pylab import savefig



Object types
============

There are several types of object in the the stat_tool module, which are
associated with an ASCII format as well as a set of methods. A full description
of each of those objects and their methods are available in the
:ref:`Reference Guide <stat_tool_reference>`.

In this tutorial we will show how to use some of them (
:func:`~openalea.stat_tool.histogram.Histogram`,
:func:`~openalea.stat_tool.convolution.Convolution`, and
:func:`~openalea.stat_tool.vectors.Vectors`) together with their methods.

Because many methods are common to those different objects, this tutorial should
be sufficient to allow you to declare and manipulate other objects such as
:func:`~openalea.stat_tool.compound.Compound`,
:func:`~openalea.stat_tool.distribution.Distribution` and
:func:`~openalea.stat_tool.vectors.VectorDistance`.


.. include:: histogram.rst
.. include:: convolution.rst
.. include:: vectors.rst


.. just to prevent annoying warnings
.. htmlonly::

    .. toctree::
        :hidden:

        convolution.rst
        histogram.rst
        vectors.rst

.. todo:: regression, data_transform, comparison, cluster, estimate, simulate
