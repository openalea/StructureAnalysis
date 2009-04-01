

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

=========================== =========================== ======================= =========================== =========================== ======================== =====================
command/data structure      convolution                 compound                Mixture                     Distribution                    Histogram            Vectors   
=========================== =========================== ======================= =========================== =========================== ======================== =====================
**Constructors**
--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
constructor (from file)     c = Convolution(filename)   c = Compound(filename)  m = Mixture(filename)       d = Distribution(filename)  h = Histogram(filename)  v = Vectors(filename)  
constructor (from dist)     c = Convolution(d1,d2)      c = Compound(d1,d2)     m = Mixture(0,1,d1,0,2,d2)  d = Distribution(d1)        h = Histogram([1,2,...]) v = Vectors([[],[]])
**display and IO**
--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
.display                    ok                          ok                      ok                          ok                          ok                       ok
Display(var)                ok                          ok                      ok                          ok                          ok                       ok
.save(filename)             ok                          ok                      ok                          ok                          ok                       ok
Save(filename)              ok                          ok                      ok                          ok                          ok                       ok
print var                   ok                          ok                      ok                          ok                          ok                       ok
var.ascii_write(True)       ok                          ok                      ok                          ok                          ok                       ok
var.file_ascii_write(True)  ok                          ok                      ok                          ok                          ok                       ok
var.spreadsheet_write(True) ok                          ok                      ok                          ok                          ok                       ok
var.survival_ascii_write()  I                           I                       I                           ok                          ok                       I
str(var)                    ok                          ok                      ok                          ok                          ok                       ok
len(var)                    ok                          I?                      ok                          I                           ok                       ok
container: var[]            I                           I?                      I?                          I?                          ok                       ok
**Plotting routines**
--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
var.plot() #gnuplot         ok                          ok                      ok                          ok                          ok                       X
var.plot() #matplotlib      ok                          runs but to be checked  ok                          ok                          ok                       X
var.print_plot()            ok                          ok                      ok                          ok                          ok
var.plot_write()            ok                          ok                      ok                          ok                          ok                       ok
**family**
--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
var.simulate                ok                          ok                      ok                          ok                          I                        I
Simulate(var, "",...)       ok                          ok                      ok                          ok                          I                        I
var.estimate_<name>()       ok                          ok                      ok                          ok (name=paramtric)         ok (name=non_parametric) I
Estimate(var, ""...)        ok                          ok                      I                           I                           I                        I
var.extract_convolution()   ok                          I                       I                           I                           I                        I  
var.extract_compound()      I                           ok                      I                           I                           I                        I
var.extract_mixture()       I                           I                       ok                          I                           I                        I 
var.extract_elementary()    ok                          ok                      I                           I                           I                        I
var.extract_component       I                           I                       ok                          I                           I                        I
var.extract_data()          ok                          ok                      ok                          TODO                        TODO                     I
var.extract_weight()        I                           I                       ok                          I                           I                        I
var.extract_sum()           I                           ok                      I?                          I                           I                        I
var.extract_elementary()    ok                          ok                      I                           I                           I                        I
--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
**specialized function**
ToHistogram                 I                           I                       I                           I                           ok
ToDistribution              I                           I                       I                           I                           ok
=========================== =========================== ======================= =========================== =========================== ======================== =====================



=========================== =========================== ======================= =========================== =============================== ========
command/data structure      convolutionData             compoundData            MixtureData                 DistributionData                Vectors
=========================== =========================== ======================= =========================== =============================== ========
**Data transformation**
----------------------------------------------------------------------------------------------------------------------------------------------------
.fit=Fit
.merge=Merge
.shift==Shift                ok                         ok                      ok                          ok                              ok
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
=========================== =========================== ======================= =========================== =============================== ========

**cluster**
cluster_step 
cluster_information     
cluster_limit
cluster.ToDistanceMatrix
cluster.Transcode
cluster._Dendrogram       
cluster.Cluster           
cluster.Clustering

**regression**

**comparison**
comparison.compare_seq       
comparison.compare_vectors   
comparison.Compare   
comparison.ComparisonTest    
comparison.compare_histo     
comparison.compare_markov  
comparison.Histogram         

**vectors**


.vector_distance_type
.ContingencyTable             ok v.contingency_table           documentation to finalise what are the extra input arguments (filename and format ?)
.interface                    
.distance_type               
.VarianceAnalysis            ok variance_analys XXXX see(issue) 
.VectorDistance            ok 

v.mixture_cluster                
v.get_identifiers     
v.mixture_estimation
v.mixture_estimation_wrap  
v.select_individual
v.cluster_limit     
v.get_nb_variable
v.select_variable
v.cluster_step     
v.get_nb_vector
v.moving_average_regression      
v.compare          
v.nearest_neighbours_regression  v.shift
v.spreadsheet_write
v.transcode
v.value_select
v.linear_regression  
v.extract                       
v.merge                         
v.file_ascii_write 
v.merge_variable                 




Current issues
==============

* Issue with the UNIFORM distribution while saving and loading a mixture of uniform distribution. See test_save in test_mixture.
* Issue with the Histogram distribution in test_save:
  if h = Histogram('mixture1.mixt/')
    len(h) returns 76 but this seem to be the length of the original data set, not the histogram itself. Is this what we want ? 
* VarianceAnalysis works but variance_analysis wrapping is not robust, not well documented, leads to crashses.
