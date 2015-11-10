.. _stat_tool_list:

.. note:: Functions/Classes with a link have been included into the python module. Tests have also been done. 

List of AML functions of the STAT module
########################################

List of AML functions from STAT: alphabetic order
=================================================


- :func:`~openalea.stat_tool.AddAbsorbingRun`
- :func:`~openalea.stat_tool.cluster.Cluster`
- :func:`~openalea.stat_tool.cluster.Clustering`
- :func:`~openalea.stat_tool.comparison.compare_histo` (distributions)
- :func:`~openalea.stat_tool.comparison.compare_vectors` (vectors)
- :func:`~openalea.stat_tool.comparison.compare_seq` (sequences)
- :func:`~openalea.stat_tool.comparison.compare_markov` (Markovian models) 
- :func:`~openalea.stat_tool.comparison.ComparisonTest`
- :func:`~openalea.stat_tool.compound.Compound` 
- `ComputeCorrelation` *in aml/src/cpp/stat_funs5.cpp* 
- `ComputePartialAutoCorrelation` 
- `ComputeRankCorrelation` 
- `ComputeSelfTransition` 
- `ComputeStateSequences` 
- `ComputeWhiteNoiseCorrelation` 
- :func:`~openalea.stat_tool.vectors.ContingencyTable` 
- :func:`~openalea.stat_tool.convolution.Convolution` 
- `Cumulate` 
- `Difference` 
- :func:`~openalea.stat_tool.output.Display` 
- :class:`~openalea.stat_tool.distribution.Distribution` 
- :class:`~openalea.stat_tool.estimate.Estimate` (distributions) 
- `Estimate` (renewal process) 
- `Estimate` (Markovian models) 
- `Estimate` ('top' parameters) 
- :func:`~openalea.stat_tool.data_transform.ExtractData` 
- :func:`~openalea.stat_tool.data_transform.ExtractDistribution`
- :func:`~openalea.stat_tool.data_transform.ExtractHistogram`
- `ExtracVectors`
- :func:`~openalea.stat_tool.data_transform.Fit`
- `HiddenMarkov` 
- `HiddenSemiMarkov` 
- :func:`~openalea.stat_tool.histogram.Histogram`
- `IndexSelect` 
- `LengthSelect` 
- `Load` 
- `Markov` 
- :func:`~openalea.stat_tool.data_transform.Merge`
- :func:`~openalea.stat_tool.data_transform.MergeVariable`
- :func:`~openalea.stat_tool.mixture.Mixture` 
- `ModelSelectionTest`
- `MovingAverage`
- `NbEventSelect`
- :func:`~openalea.stat_tool.output.Plot`, NewPlot 
- `RecurrenceTimeSequences`
- `Regression`
- `RemoveApicalInternodes`
- `RemoveRun`
- `Renewal`
- `Reverse`
- :func:`~openalea.stat_tool.output.Save`
- `SegmentationExtract`
- :func:`~openalea.stat_tool.data_transform.SelectIndividual`
- :func:`~openalea.stat_tool.data_transform.SelectVariable`
- `SemiMarkov`
- `Sequences`
- :func:`~openalea.stat_tool.data_transform.Shift`
- :class:`~openalea.stat_tool.simulate.Simulate` (distributions) 
- `Simulate` (renewal process) 
- `Simulate` (Markovian models) 
- `Simulate` ('topt' parameters) 
- `Symmetrize`
- `TimeEvents`
- `TimeScaling` 
- `TimeSelect`
- :func:`~openalea.stat_tool.cluster.ToDistanceMatrix`
- :func:`~openalea.stat_tool.distribution.ToDistribution`
- :func:`~openalea.stat_tool.distribution.ToHistogram`
- `TopParameters`
- `Tops`
- :func:`~openalea.stat_tool.cluster.Transcode`
- `TransformPosition` 
- :func:`~openalea.stat_tool.data_transform.ValueSelect`
- `VariableScaling`
- :func:`~openalea.stat_tool.vectors.VarianceAnalysis` 
- :func:`~openalea.stat_tool.vectors.VectorDistance`
- :class:`~openalea.stat_tool.vectors.Vectors`

List of AML functions from STAT: by category
============================================

Input/output functions
----------------------
- :func:`~openalea.stat_tool.compound.Compound` : construction d'un objet de type COMPOUND
- :func:`~openalea.stat_tool.convolution.Convolution`: CONVOLUTION constructor,
- Distribution: DISTRIBUTION constructor,
- HiddenMarkov: HIDDEN_MARKOV constructor,
- HiddenSemiMarkov: HIDDEN_SEMI-MARKOV constructor,
- :func:`~openalea.stat_tool.histogram.Histogram`: HISTOGRAM constructor,
- Markov: MARKOV constructor,
- Mixture: MIXTURE constructor,
- Renewal: RENEWAL constructor,
- SemiMarkov: SEMI-MARKOV constructor,
- Sequences: SEQUENCES constructor,
- TimeEvents: TIME_EVENTS constructor,
- TopParameters: TOP_PARAMETERS constructor,
- Tops: TOPS constructor,
- :func:`~openalea.stat_tool.vectors.VectorDistance`: VECTOR_DISTANCE constructor,
- Vectors: VECTORS, constructor,
- Load: restoration of an object saved as a binary file
- :func:`~openalea.stat_tool.output.Display`: ASCII output,
- :func:`~openalea.stat_tool.output.Plot`: graphical output,
- :func:`~openalea.stat_tool.output.Print`: ASCII print,
- :func:`~openalea.stat_tool.output.Save`: save in a file.

Functions of data manipulation:
-------------------------------

- :func:`~openalea.stat_tool.data_transform.Merge` merging of objects of the same 'data' type or merging of sample correlation functions,
- :class:`~openalea.stat_tool.cluster.Cluster`: clustering of values,
- :func:`~openalea.stat_tool.data_transform.Shift` shifting of values,
- :func:`~openalea.stat_tool.cluster.Transcode`: transcoding of values,
- :func:`~openalea.stat_tool.data_transform.SelectIndividual` selection of individuals,
- :func:`~openalea.stat_tool.data_transform.ValueSelect` selection of individuals according to the values taken by a variable.
- :func:`~openalea.stat_tool.data_transform.MergeVariable` merging of variables,
- :func:`~openalea.stat_tool.data_transform.SelectVariable` selection of variables.
  
set of count data of type {time interval between two observation dates, number of events occurring between these two observation dates}:

- NbEventSelect: selection of data item according to a number of events criterion,
- TimeScaling: change of the time unit,
- TimeSelect: selection of data item according to a length of the observation period criterion.

set of sequences:

- AddAbsorbingRun: addition of a run of absorbing vectors at the end of sequences,
- Cumulate: sum of successive values along sequences,
- Difference: first-order differencing of sequences,
- IndexExtract: extraction of sub-sequences corresponding to a range of index parameters,
- LengthSelect: selection of sequences according to a length criterion,
- MovingAverage: extraction of trends or residuals using a symmetric smoothing filter,
- RecurrenceTimeSequences: computation of recurrence time sequences for a given value,
- RemoveRun: removal of the first or last run of a given value (for a given variable) in a sequence,
- Reverse: reversing of sequences or 'tops',
- SegmentationExtract: extraction of sub-sequences by segmentation,
- VariableScaling: change of the unit of a variable.

set of 'tops':
  - RemoveApicalInternodes: removal of the apical internodes of the parent shoot of a 'top'.

dissimilarity matrix:
  - Symmetrize: symmetrization of a dissimilarity matrix.

Statistical functions:
----------------------
- :func:`~openalea.stat_tool.cluster.Clustering` application of clustering methods (either partitioning methods or hierarchical methods) to dissimilarity matrices between patterns,
- :func:`~openalea.stat_tool.comparison.Compare` comparison of frequency distributions, vectors, sequences, Markovian models for sequences or Markovian models,
- :func:`~openalea.stat_tool.comparison.ComparisonTest` test of comparison of frequency distributions,
- ComputeCorrelation: computation of sample autocorrelation or cross-correlation functions,
- ComputePartialAutoCorrelation: computation of sample partial autocorrelation functions,
- ComputeRankCorrelation: computation of a rank correlation matrix,
- ComputeStateSequences: computation of the optimal state sequences corresponding to the observed sequences using a hidden Markov chain or a hidden semi-Markov chain,
- ComputeWhiteNoiseAutoCorrelation: computation of the autocorrelation function induced on a white noise sequence by filtering,
- :func:`~openalea.stat_tool.vectors.ContingencyTable`: computation of a contingency table,
- :class:`~openalea.stat_tool.estimate.Estimate`: estimation of distributions, renewal processes, Markovian models or 'top' parametres from data sample,
- :func:`~openalea.stat_tool.data_transform.Fit` fit of a frequency distribution by a theoretical distribution,
- ModelSelectionTest: test for selecting the order of a Markov chain or an aggregation of states of a Markov chain,
- Regression: simple (either linear or nonparametric) regression,
- Simulate: generation of random samples from distributions, renewal processes, Markovian models or 'top' parametres,
- :func:`~openalea.stat_tool.vectors.VarianceAnalysis`: one-way variance analysis.

Miscellaneous functions
-----------------------
- ComputeSelfTransition: computation of the self-transition probabilities as a function of the index parameter from discrete sequences,
- :func:`~openalea.stat_tool.data_transform.ExtractData` extraction of the 'data' part of an object of type 'model',
- :func:`~openalea.stat_tool.data_transform.ExtractDistribution` extraction of a distribution from an object of type 'model',
- :func:`~openalea.stat_tool.data_transform.ExtractHistogram` extraction of a frequency distribution from an object of type 'data',
- ExtractVectors: extraction of vectors from global characteristics of sequences (length or counting characteristics),
- :func:`~openalea.stat_tool.cluster.ToDistanceMatrix` cast of an object of type CLUSTERS into an object of type DISTANCE-MATRIX
- :func:`~openalea.stat_tool.distribution.ToDistribution`: cast of an object of type HISTOGRAM into an object of type DISTRIBUTION
- :func:`~openalea.stat_tool.distribution.ToHistogram`: cast of an object of type DISTRIBUTION into an object of type HISTOGRAM
- TransformPosition: discretization of inter-position intervals. 

List by type
============

type clusters
-------------

* function returning an object of type CLUSTERS:
  - Load
  - :func:`~openalea.stat_tool.cluster.Clustering`  
* function taking as argument an object of type CLUSTERS:
  - :func:`~openalea.stat_tool.output.Display`
  - :func:`~openalea.stat_tool.output.Plot`
  - :func:`~openalea.stat_tool.output.Print`
  - :func:`~openalea.stat_tool.output.Save`
  - :func:`~openalea.stat_tool.cluster.ToDistanceMatrix`
