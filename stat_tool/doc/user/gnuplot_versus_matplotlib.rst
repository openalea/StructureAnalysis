.. testsetup:: *

    from openalea.stat_tool import *
    import pylab
    from pylab import savefig, clf
    h1 = Histogram('./test/data/meri1.his')


Gnuplot and Matplotlib Plots
============================

Switching between  Gnuplot and/or Matplotlib
--------------------------------------------

::

    >>> h1 = Histogram('./test/data/meri1.his')

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


Saving in postscript with gnuplot
---------------------------------

There are other methods related to GNUPLOT (if implemented to output the results
in a file or ps file)

.. doctest::
    :options: +SKIP

    >>> h1.plot_write('output', 'title')
    >>> h1.plot_print() # save gnuplot output in a postscript file
