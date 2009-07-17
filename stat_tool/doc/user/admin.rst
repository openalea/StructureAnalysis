
Current Developments
####################

.. contents::

.. methodology (in french)



.. include:: ../../methodo.txt



Work in progress
================

documentation
-------------

The sphinx documentation is now used to manage the documentation so as to 
include a reference guide as well as a User guide (and this administration
section).

All docstrings are currently compatible with the sphinx requirements (reST
 syntax).

aml syntax
----------
Syntax shoudl stick to English or American over all the functions.


Wrapper and validation 
----------------------

Here is a list of data structures together with the methods that are available, 
either through the boost_python interface or directly from python.


.. note:: In the following table, 
    
    * `X` means that the function does not work or is not finalised, 
    * `ok` means implemented, tested and working
    * `I` or `irrelevant` means means irrelevant: this functionality does not exist for this class


There are classes/modules related to data types like compound.py and modules
related to functions like cluster.py. First, let us look at the data type
 
**Data types**


In general there are two constructors from a filename or from distributions::

    Convolution(filename)
    Compound(filename)
    Mixture(filename)
    Distribution(filename)
    Histogram(filename)
    Vectors(filename)  

and

::

    Convolution(d1,d2)
    Compound(d1,d2)
    Mixture(0,1,d1,0,2,d2)
    Distribution(d1)
    Histogram([1,2,...])
    Vectors([[],[]])

All of those types have common methods that are summarized in the following 
table. Here are some aliaes being used in the table.

    compound        = 1;
    convolution     = 2;
    distribution    = 3;
    histogram       = 4;
    mixture         = 5;
    vector          = 6;
    matrix          = 7;
    mv_mixture      = 8;
    regression      = 9;
    vector_distance = 10
    

=========================== ======= ======= ======= ======= ======= ======= ======
Basic methods               1       2       3       4       5       6       7    
=========================== ======= ======= ======= ======= ======= ======= ======
**Constructors**
From filename               ok      ok      ok      ok      ok      ok      X
From distribution           ok      ok      ok      ok      ok      ok      X
**methods**
ascii_write                 ok      ok      ok      ok      ok      ok      ok
display==Display            ok      ok      ok      ok      ok      ok      ok
extract_data                ok      ok      ok      ok      ok      ok      ok
file_ascii_write            ok      ok      ok      ok      ok      ok      ok
plot==Plot                  ok      ok      ok      ok      ok      **I?**  **I?**
save==Save                  ok      ok      ok      ok      ok      ok      ok
plot_print                  ok      ok      ok      ok      ok      ok      ok
simulate==Simulate          ok      ok      ok      ok      ok      ok      ok
plot_write                  ok      ok      ok      ok      ok      ok      ok
spreadsheet_write           ok      ok      ok      ok      ok      ok      ok
str                         ok      ok      ok      ok      ok      ok      ok
len                         I       ok      I       ok      ok      ok      I
old_plot                    ok      ok      ok      ok      ok      ok      ok
=========================== ======= ======= ======= ======= ======= ======= ======

**Specific methods**


======================= ===========
Compound                status
======================= ===========
extract_compound        ok
extract_elementary      ok  
extract_sum             ok
======================= ===========

and 

======================= ===========
Convolution             status
======================= ===========
extract_convolution     ok
extract_elementary      ok  
======================= ===========

and

=========================== ===========
Distribution                status
=========================== ===========
survival_ascii_write        ok
survival_spreadsheet_write  ok
survival_plot_write         ok
ident                       ok            
probability                 ok
sup_bound                   ok
parameter                   ok
inf_bound                   ok     
=========================== ===========

and

=========================== ===========
Histogram                   status
=========================== ===========
spreadsheet_write           ok
survival_ascii_write        ok
survival_spreadsheet_write  ok
extract_model               ok
**see test_cluster**
cluster_information         ok
cluster_limit               ok
cluster_step                ok
transcode                   ok
**see test_comparison**
compare                     ???
compare_histo               ok
t_comparison                ok
wmw_comparison              ok
f_comparison                ok
**see test_estimate**
estimate_compound           ok
estimate_convolution        ok     
estimate_mixture            ok 
estimate_nonparametric      ok
estimate_parametric         ok
compound_estimation         X
convolution_estimation      X
mixture_estimation          X
parametric_estimation       X
**see data_transform**
fit                         ok
merge                       ok
shift                       ok
value_select                ok  
=========================== ===========

and

======================= ===========
Mixture                 status
======================= ===========
extract_weight          ok
extract_mixture         ok  
======================= ===========

and

=========================== ===========
Vectors                     status
=========================== ===========
survival_ascii_write        ok
survival_spreadsheet_write  ok
get_nb_vector               ok
get_nb_variable             ok
get_identifiers             ok
contingency_table           ok
variance_analysis           ok
extract                     seg fault e.g. str(v.extract(1))
**see test_data_transform**
merge                       ok
shift                       ok
merge_variable              ok      
select_individual           ok
select_variable             ok
value_select                ok
**see test_regression**
linear_regression           ok
moving_average_regression   ok
nearest_neighbours_regress  ok
**see test_cluster**
cluster_limit               ok
cluster_step                ok
transcode                   to be done
compare                     ok     
mixture_cluster             to be done
mixture_estimation          to be done
mixture_estimation_wrap     to be done
=========================== ===========

and

=========================== ===============
Matrix                      status
=========================== ===============
hierarchical_clustering     ok
partitioning_clusters       ok
partitioning_prototype      ok
get_nb_row                  ok
get_nb_column               ok
=========================== ===============

and

=========================== ==========================
others
=========================== ==========================
ExtractHistogram            ok
ExtractDistribution         ok
=========================== ==========================

Other methods
=============

=================================================   ===
**cluster** (histogram)
=================================================   ===
cluster_step==Cluster(step=...)                     ok
cluster_information==Cluster(information=...)       ok
cluster_limit==Cluster(limit=...)                   ok
cluster.ToDistanceMatrix
transcode==Transcode                                ok
cluster._Dendrogram       
cluster.Clustering
=================================================   ===


=============================================================================   ===
**regression**
=============================================================================   ===
Regression(v, "Linear")==v.linear_regression(...)                               ok
Regression(v, "MovingAverage")==v.moving_average_regression(...)                ok
Regression(v, "NearestNeighbours")==v.nearest_neighbours_regression(...)        ok
=============================================================================   ===

=================================================   ===
**comparison**
=================================================   ===
comparison.compare_seq       
comparison.compare_vectors   
comparison.Compare   
comparison.ComparisonTest    
comparison.compare_histo     
comparison.compare_markov  
comparison.Histogram         
=================================================   ===



TODO
====

* Check that tests run under windows
* Issue with the UNIFORM distribution while saving and loading a mixture of
  uniform distribution. See test_save in test_mixture.
* Issue with the Histogram distribution in test_save:
  if h = Histogram('mixture1.mixt/')
  len(h) returns 76 but this seem to be the length of the original data set,
  not the histogram itself. Is this what we want ? 
* VarianceAnalysis works but variance_analysis wrapping is not robust, not well
  documented, leads to crashses.
* plotable inside MvMixture(Data) ? 
* plotable in vectors
* check Display and Ploy functions (viewpoint, Details, ...)

DONE
====

* language for the syntax is American
* Comments all functional tests and add title to the Plot commands 
