.. define some aliases:
.. _histogram: syntax.html#type-histogram

.. define the setup for doctest:
.. testsetup:: *
   
    from openalea.stat_tool import *
    import pylab
    from pylab import savefig, clf



Histogram
=========
Here is a brief description of the Histogram type.

Constructor
-----------

Histogram can be generated either by loading an ascii file -- following the syntax given in histogram_ -- or with a list of values, that will be processed into an histogram.

So, if you have an ASCII file, you can load it using the :func:`~openalea.stat_tool.histogram.Histogram` function as follows:

.. filename with respect to the directory where sphinx is launch
.. doctest::

    >>> h1 = Histogram('./test/meri5.his')

Otherwise, if you have a list of values, provide it using python syntax as follows:

.. doctest::

    >>> h2 = Histogram([1,2,2,3,4,4,4,5])

Now, you can use the methods bounded to the `Histogram` class. Many methods are also functions. To differentiate them, the first letter of the function is capitalized (see next example)


display
-------
The object `h` has a few methods, The `display` method :func:`~openalea.stat_tool.output.Display` returns information on the screen
   
.. doctest::

    >>> h1.display()
    'histogram - sample size: 66\nmean: 4.37879   variance: 1.62354   standard deviation: 1.27418\ncoefficient of skewness: 0.0727983   coefficient of kurtosis: -0.709664\nmean absolute deviation: 1.06841   coefficient of concentration: 0.161214\ninformation: -107.512 (-1.62897)\n'
    >>> Display(h1)
    'histogram - sample size: 66\nmean: 4.37879   variance: 1.62354   standard deviation: 1.27418\ncoefficient of skewness: 0.0727983   coefficient of kurtosis: -0.709664\nmean absolute deviation: 1.06841   coefficient of concentration: 0.161214\ninformation: -107.512 (-1.62897)\n'
    

printing ASCII information :

.. doctest::

    >>> print h1.ascii_write(True) #doctest: +SKIP
    >>> print h1.ascii_write(False) #doctest: +SKIP
    histogram - sample size: 66
    mean: 4.37879   variance: 1.62354   standard deviation: 1.27418
    coefficient of skewness: 0.0727983   coefficient of kurtosis: -0.709664
    mean absolute deviation: 1.06841   coefficient of concentration: 0.161214
    information: -107.512 (-1.62897)
    >>>

plotting
--------

old AML style

.. doctest::
    :options: +SKIP
    
    h.plot_print()

new style, either with GNUPLOT or MATPLOTLIB. By default, matplotlib is used:

.. doctest::
    
    >>> clf()
    >>> h1.plot(show=False)
    >>> savefig('doc/user/stat_tool_histogram_plot.png')
    >>> # by default, the Plot routine uses matplolib (if availabl)
    >>> # but you can still use gnuplot 
    >>> plot.set_plotter(plot.gnuplot()) #doctest: +SKIP
    >>> # and come back to matplotlib 
    >>> plot.set_plotter(plot.mtplotlib()) #doctest: +SKIP


.. figure:: stat_tool_histogram_plot.png
    :width: 50%
    :align: center

save in a gnuplot file with plot_write method::

    >>> h1.plot_write('output', 'title')


    >>> h1.print_plot() # save gnuplot output in a postscript file
clustering
-----------

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


Merging
-------

the following examples illustrates the usage of the :func:`¬openalea.stat_tool.data_transform.Merge` function. See also Figure :ref:`fig_merging` for the output plots.

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




.. h.estimate_convolution        h.plot_write h.estimate_mixture            h.save h.estimate_nonparametric      h.shift h.estimate_parametric         h.spreadsheet_write h.extract_model               h.survival_ascii_write h.f_comparison                h.survival_get_plotable h.file_ascii_write            h.survival_plot_write h.fit                         h.survival_spreadsheet_write h.get_plotable                h.t_comparison h.compare                     h.transcode h.compare_histo               h.mixture_estimation          h.value_select h.compound_estimation         h.old_plot                    h.wmw_comparison h.convolution_estimation      h.parametric_estimation        h.estimate_compound           h.plot_print                  
