

Administration documents
########################

.. contents::

.. methodology (in french)
.. include:: ../../methodo.txt


Work in progress
================

documentation
-------------

The sphinx documentation is now used to manage the documentation so as to include a reference guide as well as a User guide (and this administration section).

All docstrings are currently compatible with the sphinx requirements (reST syntax)

Wrapper and validation 
----------------------

Here is a list of data structures together with the methods that are available, either through the boost_python interface or directly from python.


.. note:: In the following table, 
    
    * `X` means that the function does not work or is not finalised, 
    * `ok` means implemented,tested and working
    * `I` or `irrelevant` means means irrelevant: this functionality does not exist for this class

Class family
++++++++++++

=========================== =========================== ======================= =========================== =============================== ========================
command/data structure      convolution                 compound                Mixture                     Distribution                    Histogram
=========================== =========================== ======================= =========================== =============================== ========================
**Constructors**
--------------------------------------------------------------------------------------------------------------------------------------------------------------------
constructor (from file)     c = Convolution(filename)   c = Compound(filename)  m = Mixture(filename)       d = Distribution(filename)      h = Histogram(filename)
constructor (from dist)     c = Convolution(d1,d2)      c = Compound(d1,d2)     m = Mixture(0,1,d1,0,2,d2)  d = Distribution(d1)            h = Histogram([1,2,...])
**display and IO**
--------------------------------------------------------------------------------------------------------------------------------------------------------------------
.display                    ok                          ok                      ok                          ok                              ok
Display(var)                ok                          ok                      ok                          ok                              ok
.save(filename)             ok                          ok                      ok                          ok                              ok
Save(filename)              ok                          ok                      ok                          ok                              ok
print var                   ok                          ok                      ok                          ok                              ok
var.ascii_write(True)       ok                          ok                      ok                          ok                              ok
var.file_ascii_write(True)  ok                          ok                      ok                          ok                              ok
var.spreadsheet_write(True) ok                          ok                      ok                          ok                              ok
var.survival_ascii_write()  I                           I                       I                           ok                              ok
str(var)                    ok                          ok                      ok                          ok                              ok
len(var)                    ok                          I?                      ok                          I                               ok
container: var[]            I                           I?                      I?                          I?                              ok
**Plotting routines**
--------------------------------------------------------------------------------------------------------------------------------------------------------------------
var.plot() #gnuplot         ok                          ok                      ok                          ok                              ok
var.plot() #matplotlib      ok                          runs but to be checked  ok                          ok                              ok
var.print_plot()            ok                          ok                      ok                          ok                              ok
var.plot_write()            ok                          ok                      ok                          ok                              ok
**family**
--------------------------------------------------------------------------------------------------------------------------------------------------------------------
var.simulate                ok                          ok                      ok                          ok                              I
Simulate(var, "",...)       ok                          ok                      ok                          ok                              I
var.estimate_<name>()       ok                          ok                      ok                          ok (name=paramtric)             ok (name=non_parametric)
Estimate(var, ""...)        ok                          ok                      I                           I                               I
var.extract_convolution()   ok                          I                       I                           I                               I
var.extract_compound()      I                           ok                      I                           I                               I
var.extract_mixture()       I                           I                       ok                          I                               I
var.extract_elementary()    ok                          ok                      I                           I                               I
var.extract_component       I                           I                       ok                          I                               I
var.extract_data()          ok                          ok                      ok                          TODO                            TODO
var.extract_weight()        I                           I                       ok                          I                               I
var.extract_sum()           I                           ok                      I?                          I                               I
var.extract_elementary()    ok                          ok                      I                           I                               I
--------------------------------------------------------------------------------------------------------------------------------------------------------------------
**specialized function**
ToHistogram                 I                           I                       I                           I                               ok
ToDistribution              I                           I                       I                           I                               ok
=========================== =========================== ======================= =========================== =============================== ========================



=========================== =========================== ======================= =========================== ===============================
command/data structure      convolutionData             compoundData            MixtureData                 DistributionData               
=========================== =========================== ======================= =========================== ===============================
**Data transformation**
-------------------------------------------------------------------------------------------------------------------------------------------
fit
merge
shift
cluster_information     
cluster_limit
get_plotable            
t_comparison
cluster_step 
transcode
compare
ms.value_select
compare_histo        
wmw_comparison
merge                 
f_comparison            
=========================== =========================== ======================= =========================== ===============================

**cluster**
cluster.ToDistanceMatrix
cluster.Transcode
cluster._Dendrogram       
cluster.Cluster           
cluster.Clustering

**regression**

comparison.compare_seq       
comparison.compare_vectors   
comparison.Compare           comparison.ComparisonTest    
comparison.compare_histo     
comparison.compare_markov    comparison.Histogram         

**vectors**


vectors.vector_distance_type
vectors.Vectors
vectors.ContingencyTable             
vectors.interface                    
vectors._Vectors_mixture_estimation
vectors.distance_type               
vectors.VarianceAnalysis             
vectors.VectorDistance             

v.ascii_data_write               v.mixture_cluster                
v.ascii_write                    v.get_identifiers                v.mixture_estimation             v.save
v.mixture_estimation_wrap        v.select_individual
v.cluster_limit                  v.get_nb_variable                v.select_variable
v.cluster_step                   v.get_nb_vector                  v.moving_average_regression      
v.compare                        v.get_plotable                   v.nearest_neighbours_regression  v.shift
v.contingency_table             __new__                        v.spreadsheet_write
v.old_plot                       
v.plot                           v.transcode
v.display                        v.plot_print                     v.value_select
v.linear_regression              v.plot_write                     v.variance_analysis
v.extract                        v.merge                         _
v.file_ascii_write               v.merge_variable                 





Function family
+++++++++++++++

=========== =========== ==========
Estimate    Simulate    Cumulate
=========== =========== ==========
in progress in progress inprogress
=========== =========== ==========
