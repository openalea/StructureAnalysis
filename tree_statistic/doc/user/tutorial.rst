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
    
    .. code-block:: python
        from openalea.stat_tool import *
        from openalea.tree_statistic.trees import *
        from openalea.tree_statistic.etrees import *
        from openalea.tree_statistic.hmt import *
        import openalea
        from pylab import savefig


.. .. todo:: uncomment the tutorial_object_types.rst in the tutorial.rst file once commited

.. include:: tutorial_object_types.rst

.. include:: tree.rst

..    .. include:: trees.rst
..    .. include:: hmt.rst


.. just to prevent annoying warnings 
.. htmlonly::

    .. toctree::
        :hidden:
    
