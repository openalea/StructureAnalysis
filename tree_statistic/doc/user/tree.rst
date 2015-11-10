.. define the setup for doctest:
.. testsetup:: *
   

    from openalea.stat_tool import *
    from openalea.tree_statistic.trees import *

TreeStructure
=============

Here is a brief description of 
the :class:`~openalea.tree_statistic.trees.TreeStructure`  class.

The *TreeStructure*
class is used for the representation of single trees
(as opposed to forests). No attributes can be attached to the vertices
(as opposed to the `Tree`_ class).

Constructor
-----------

Instances of class TreeStructure can be obtained either by simulation,
or by extraction from a `Tree`_ object.

Building a *TreeStructure* instance by simulation requires the distribution of the number
of children for each vertex to be specified; *e.g.* a discrete uniform
distribution on 0, 1, 2, 3. The tree is then generated using a branching
process (Galton-Watson process). A maximal depth and a maximal number
of children also have to be specified, and are used as stopping criterion
for the simulation algorithm.

.. filename with respect to the directory where sphinx is launch

.. doctest::

    >>> inf_bound = 0
    >>> sup_bound = 3
    >>> distrib = Uniform(inf_bound, sup_bound)
    >>> max_depth = 3
    >>> max_size = 10
    >>> R = TreeStructure(distrib, max_size, max_depth)

TreeValue
=========

The :class:`~openalea.tree_statistic.trees.TreeValue` 
class is used to represent the attributes attached
to each vertex of a `Tree`_ object. It is very much comparable
to a python list with a fixed number of elements, and containing
only objects of type 'int' and 'float'. As in python lists,
the first element of a *TreeValue* has index 0.

Constructor
-----------

Instances of class `TreeValue`_ can be obtained either by conversion from
a list into a `TreeValue`_ instance, or by extraction from a `Tree`_
object, through the method :ref:`Get <get-put>`.

.. doctest::

    >>> tv = TreeValue([1., 0])
    >>> print tv #doctest: +SKIP
    [1.0, 0]
    >>> print tv[0] #doctest: +SKIP
    [1.0]

Tree
====

Here is a brief description of the :class:`~openalea.tree_statistic.trees.Tree`
class.

The *Tree* class is used for the representation of single trees
(as opposed to forests). To each vertex of a *Tree* object, several
attributes (or variables) can be attached. The number of variables
must be the same for each vertex.

Constructor
-----------

Instances of class `Tree`_ can be obtained either by addition of variables
to a `TreeStructure`_ object, by extraction from a `Trees`_,
or from a :ref:`MTG` file (see section `Trees`_).
If this file contains more than a tree, only
the first tree of the file will be returned. 

To build a `Tree`_ object from a `TreeStructure`_ object, a default value
for the attributes has to be provided in the constructor.

.. doctest::

    >>> tv = [1., 0, 1, 2.]
    >>> T1 = Tree(tv, R)

Attributes
----------

The object **T1** has a few methods, among which some aims at printing
information on the screen.
The :meth:`~openalea.tree_statistic.trees.Tree.Attributes` method is one
of them. This methods prints the name of the attributes 
(also named variables). The default
name of attribute number `i` is `"Variable"+str(i)`. The attributes
can only be changed if the `Tree`_ instance was built from a MTG file
(see below).

The number of attributes can be accessed to, using 
the :meth:`~openalea.tree_statistic.trees.Tree.NbVariables()` method.
This number is the same for all the vertices of a given instance, 
and cannot be changed directly for this instance. 

The attributes have types, which can be either `INT_VALUE` (for integer 
attributes) or `REAL_VALUE` (for floating attributes). The type of 
variable `i` in instance **T1** is given by `T1.Type(i)`, and the list
of every types by `T1.Types()`. 
The :meth:`~openalea.tree_statistic.trees.Tree.NbInt` method returns
the number of attributes with type `INT_VALUE`, whereas the 
:meth:`~openalea.tree_statistic.trees.Tree.NbFloat` method returns
the number of attributes with type `REAL_VALUE`.

.. doctest::

    >>> print(T1.NbVariables()) #doctest +SKIP
    4
    >>> print(T1.NbInt())
    2
    >>> print(T1.NbFloat())
    2

Display
-------

The :meth:`~openalea.tree_statistic.trees.Tree.Display` method
provides an ASCII output of a `Tree`, as shown below:
   
.. doctest::

    >>> T1.Display() #doctest: +SKIP
    vids: [ Variable0, Variable1, Variable2, Variable3 ]
    0: [1.0, 0, 1, 2.0]
    |-1: [1.0, 0, 1, 2.0]+
    | |-4: [1.0, 0, 1, 2.0]+
    |
    |-2: [1.0, 0, 1, 2.0]+
    | |-5: [1.0, 0, 1, 2.0]+
    |
    |-3: [1.0, 0, 1, 2.0]+

The first line sums up the semantics of every further line, as follows.
The numbers preceding the colons *:* denote the vertex identifiers,
or `vids`.
The quantities between brackets are the values of the attributes
for each vertex.
Thus, the first line also contains the names of the variables.

.. note:: The *Display* method has the following optional arguments,
    which take boolean values: *vids*, *attributes* and *mtg_vids*
    (which must be equal to **True** if, respectively, the vids,
    the attributes and the corresponding vids in
    the MTG where the data come from, have to be
    printed). The *mtg_vids* argument is available only if the `Tree`_
    instance was built from a MTG.

The :func:`str` and the :func:`print` functions have the same
effects than .Display(vids=False)

.. doctest::
    
    >>> print(T1) #doctest: +SKIP
    [ Variable0, Variable1, Variable2, Variable3 ]
    [1.0, 0, 1, 2.0]
    |-[1.0, 0, 1, 2.0]+
    | |-[1.0, 0, 1, 2.0]+
    |
    |-[1.0, 0, 1, 2.0]+
    | |-[1.0, 0, 1, 2.0]+
    |
    |-[1.0, 0, 1, 2.0]+

Root vertex
-----------

Root vertex
-----------

The :func:`~openalea.tree_statistic.trees.Tree.Root` method
returns the vertex identifier of the root vertex, as illustrated
below:

.. doctest::

    >>> T1.Root() #doctest: +SKIP
    0

..  _get-put:

Changing the attribute values
-----------------------------

Now in our example, all attributes have the same value for
every vertex. To access the value of the attributes for a given vertex,
the :meth:`~openalea.tree_statistic.trees.Tree.Get` method
has to be called (using as argument the vid of a vertex).
The :meth:`~openalea.tree_statistic.trees.Tree.Put` method allows
these attributes to be changed (using as arguments the vid of
the considered vertex and the list of the new values).

.. doctest::

    >>> T1.Root() #doctest: +SKIP
    0
    >>> print(T1.Get(T1.Root())) #doctest: +SKIP
    [1.0, 0, 1, 2.0]
    >>> T1.Put(T1.Root(), [3.1, 5, 8, -2.2])
    >>> print(T1.Get(T1.Root())) #doctest: +SKIP
    [3.1000000000000001, 5, 8, -2.2000000000000002]

..  _tree-save:

Saving
------

Any `Tree`_ instance can be saved into a file, using the
:meth:`~source_openalea.tree_statistic.trees.Tree.Save` method:

.. doctest::
    :options: +SKIP
    
    >>> T1.Save('test.mtg', variable_names=["V1", "V2", "V3", "V4"])

This creates a file *test.mtg* where the tree structure is stored,
as well as the attributes, under the mtg format.

.. note:: In addition to the file name, *Save* takes two optional
    arguments: the boolean `overwrite` (= **True** if any existing
    file with the same name can be overwritten) and the list
    *variable_names* of the names of the variables, if they have
    to be renamed in the MTG file.

Then, you can construct a new `Tree`_ instance as follows:

.. doctest::
    :options: +SKIP
    
    >>> T2 = Tree('test.mtg')


Trees
=====

Here is a brief description of the :class:`~openalea.tree_statistic.trees.Trees`
class.

The *Trees* class is used for the representation of forests, *i.e.* collections
of `Tree`_ objects

Constructor
-----------

Instances of class `Trees`_ can be obtained either from a list of `Tree`_
objects, a `Trees`_ object (by copy), or from a MTG file.

To build an instance of `Trees`_ from a list of  `Tree`_ objects, its constructor
must be called using that list as argument.

.. doctest::

    >>> F = Trees([T1, T2])
    >>> F.Save('test_forest.mtg', variable_names=["Ln", "Fl", "Lf", "LgEn"])

.. note:: The *Save* method is used exactly the same way 
    as for `Tree`_ objects - see Section :ref:`Saving <tree-save>`.
 
To copy a `Trees`_ instance, its constructor must be called using 
that instance as argument.

.. doctest::

    >>> F2 = Trees(F)

To build an instance of `Trees`_ from a MTG file, the name of that file must be 
provided as argument of the constructor.

.. note:: In addition to the file name, this constructor has two optional
    arguments : a filter on the vertices (which is a boolean function), 
    a list of strings defining the attribute names, a list of functions 
    defining how the attributes are computed, and the considered scale of the MTG.
    The default value for the scale is the finest scale of the MTG (*i.e.*
    the scale with highest value).

.. doctest::

    >>> F3 = Trees('test_forest.mtg')

Attributes
----------

The names of the attributes for the objects **F1**, **F2** and **F3**
can be accessed to using the method :meth:`~openalea.tree_statistic.trees.Trees.Attributes`.

.. doctest::

    >>> print(F3.Attributes()) #doctest: +SKIP
    ['Ln', 'Fl', 'Lf', 'Lge']


Display
-------

The :meth:`~openalea.tree_statistic.trees.Trees.Display` method
provides an ASCII output of a `Tree`, as shown below:
   
.. doctest::

    >>> T1.Display() #doctest: +SKIP
    vids: [ Variable0, Variable1, Variable2, Variable3 ]
    0: [1.0, 0, 1, 2.0]
    |-1: [1.0, 0, 1, 2.0]+
    | |-4: [1.0, 0, 1, 2.0]+
    |
    |-2: [1.0, 0, 1, 2.0]+
    | |-5: [1.0, 0, 1, 2.0]+
    |
    |-3: [1.0, 0, 1, 2.0]+

The first line sums up the semantics of every further line, as follows.
The numbers preceding the colons *:* denote the vertex identifiers,
or `vids`.
The quantities between brackets are the values of the attributes
for each vertex.
Thus, the first line also contains the names of the variables.

.. note:: The *Display* method has the following optional arguments,
    which take boolean values: *vids*, *attributes* and *mtg_vids*
    (which must be equal to **True** if, respectively, the vids,
    the attributes and the corresponding vids in
    the MTG where the data come from, have to be
    printed). The *mtg_vids* argument is available only if the `Tree`_
    instance was built from a MTG.

The :func:`str` and the :func:`print` functions have the same
effects than .Display(vids=False)

.. todo:: Tree: EdgeType, Size, SelectSubTree
    
    Trees: Merge, SelectVariable, SelectIndividual, Shift, Transcode, Cluster, 
    Merge, MergeVariable, NbVariables, Plot

.. keep paragraph below
    However, it is possible
    to create other instances by adding or deleting attributes 
    (see :ref:`Merging and Selecting variables <tree-save>`.)

..  comment lines below
    Plotting
    --------

    old AML style

    .. doctest::
        :options: +SKIP
        
        h.old_plot()

    new style, either with GNUPLOT or MATPLOTLIB. By default, matplotlib is used if
    it is implemented:

    .. doctest::
        
        >>> clf()
        >>> h1.plot(show=False)
        >>> savefig('doc/user/stat_tool_histogram_plot.png')
        >>> # by default, the Plot routine uses matplolib (if available)
        >>> # but you can still use gnuplot 
        >>> plot.set_plotter(plot.gnuplot()) #doctest: +SKIP
        >>> # and come back to matplotlib later on
        >>> plot.set_plotter(plot.mtplotlib()) #doctest: +SKIP


    .. figure:: stat_tool_histogram_plot.png
        :width: 50%
        :align: center

    There are other methods related to GNUPLOT that we will not supported anymore
    in the future::

        >>> h1.plot_write('output', 'title')
        >>> h1.print_plot() # save gnuplot output in a postscript file

    Clustering
    ----------

    Histograms can be clustered. See :func:`~openalea.stat_tool.cluster.Cluster`

    .. doctest::
        :options: +SKIP

        >>> h1.cluster_information(0.5) 
        # equivalently
        >>> Cluster(h1, "Information", 0.5)
        >>> h1.cluster_limit([1,2])
        # equivalently
        >>> Cluster(h1, "Limit", [1,2])
        >>> h1.cluster_step(3)
        # equivalently
        >>> Cluster(h1, "Step", 3)
        
    .. warning:: Again, although the function is equivalent to the method, we 
        advice you to use the functions. See Display section for details.


    Merging
    -------

    the following examples illustrates the usage of the 
    :func:`~openalea.stat_tool.data_transform.Merge` function. See also 
    Figure :ref:`fig_merging` for the output plots.

    .. doctest::

        >>> # load two histograms
        >>> h1 = Histogram('./test/meri1.his')
        >>> clf(); h1.plot(show=False); savefig('doc/user/stat_tool_histogram_h1.png')
        >>> h5 = Histogram('./test/meri5.his')
        >>> clf(); h5.plot(show=False); savefig('doc/user/stat_tool_histogram_h5.png')

    The two original histograms are shown here below:

    +---------------------------------------+----------------------------------------+
    | .. image:: stat_tool_histogram_h1.png | .. image:: stat_tool_histogram_h5.png  |
    |     :width: 100%                      |     :width: 100%                       |
    +---------------------------------------+----------------------------------------+

    .. doctest::

        >>> a = Merge(h1,h5)
        >>> b= h1.merge([h5])
        >>> c = h5.merge([h1])
        >>> clf(); a.plot(show=False)
        >>> savefig('doc/user/stat_tool_histogram_merging.png')

    .. _fig_merging:
    .. figure:: stat_tool_histogram_merging.png
        :width: 50%
        :align: center

        **Figure: The merging of two histograms**




