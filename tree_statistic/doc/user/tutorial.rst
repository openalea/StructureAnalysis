.. _tree_statistic_tutorial:

.. define the setup for doctest:
.. testsetup:: *
   
    from openalea.stat_tool import *
    from openalea.tree_statistic.trees import *
    from openalea.tree_statistic.etrees import *
    from openalea.tree_statistic.hmt import *
    import openalea
    from pylab import savefig


Tutorial
########

.. contents::

.. note:: most of the docstrings provided in this tutorial are tested with 
    doctest and should work out of the box. However, you will need to import 
    relevant modules.  For instance:
    
    >>> from openalea.stat_tool import *
    >>> from openalea.tree_statistic.trees import *
    >>> from openalea.tree_statistic.etrees import *
    >>> from openalea.tree_statistic.hmt import *
    >>> import openalea
    >>> from pylab import savefig

    
    

Object types
============

There are several types of object in the the tree_statistic module.
A full description of each of them is available in the
:ref:`Reference Guide <tree_statistic_reference>`.

In this tutorial we will show how to use some of them (
:func:`~openalea.tree_statistic.trees.Tree`, 
:func:`~openalea.tree_statistic.trees.Trees`, and
:func:`~openalea.tree_statistic.hmt.HiddenMarkovTree`), together with their methods.

.. include:: tree.rst
..
    .. include:: trees.rst
    .. include:: hmt.rst


.. just to prevent annoying warnings 
.. htmlonly::

    .. toctree::
        :hidden:
    
