
# This file has been generated at Wed Aug 11 17:21:31 2010

from openalea.core import *


__name__ = 'vplants.stats'

__editable__ = True
__description__ = 'stat_tool and sequence analysis nodes.'
__license__ = 'CECILL-C'
__url__ = 'http://openalea.gforge.inria.fr/doc/vplants/stat_tool/doc/_build/html/contents.html'
__alias__ = []
__version__ = '1.2.0'
__authors__ = 'OpenAlea consortium'
__institutes__ = 'INRIA/CIRAD'
__icon__ = 'icon.png'


__all__ = ['stat_py_value_select', 'stat_py_convolution', 'stat_py_simulate_dist', 'stat_py_dist_poisson', 'stat_py_semimarkov', 'stat_py_segmentation_sample', 'stat_py_shiftn', 'stat_py_pointwise_average', 'stat_py_estimate_mixture', 'stat_PyObjectFromFile', 'stat_py_compute_correlation', 'stat_py_vectors', 'stat_py_cluster', 'stat_py_select_variable', 'stat_py_dist_binomial', 'stat_py_sequences', 'stat_py_compute_correlation_mult', 'stat_py_to_distribution', 'stat_py_hiddensemimarkov', 'stat_py_to_histogram', 'stat_py_compound', 'stat_py_regression', 'stat_py_plot_segprofile', 'stat_PyRenewal', 'stat_py_merge', 'stat_py_plot', 'data_Selector', 'stat_py_comparisontest', 'stat_py_estimate_compound', 'stat_py_cumulate', 'stat_py_estimate_dist', 'stat_py_extract_distribution', 'stat_py_shift', 'stat_py_dist_negativebinomial', 'stat_py_markov', 'stat_py_dist_uniform', 'stat_py_fit', 'stat_py_segmentation', 'stat_py_select_individual', 'stat_py_hiddenmarkov', 'stat_py_extract_histogram', 'stat_py_compare_vectors', 'stat_py_histogram', 'stat_PyRenewalAscii', 'stat_py_extract_data', 'stat_py_compare_frequency', 'stat_py_estimate_conv', 'stat_py_extract_histogram_sequences']



stat_py_value_select = Factory(name='ValueSelect',
                description='Extraction of vectors according to a a criteria upon a given variable',
                category='STAT',
                nodemodule='sequence_analysis_wralea.stat',
                nodeclass='py_value_select',
                inputs=[{'name': 'obj'}, {'interface': IInt, 'name': 'variable', 'value': 1}, {'interface': IInt, 'name': 'value', 'value': None}],
                outputs=None,
                widgetmodule=None,
                widgetclass=None,
               )




stat_py_convolution = Factory(name='Convolution',
                description='Construction of an object of type convolution from an elementary distribution.',
                category='STAT.FrequencyDistribution',
                nodemodule='sequence_analysis_wralea.stat',
                nodeclass='py_convolution',
                inputs=None,
                outputs=None,
                widgetmodule=None,
                widgetclass=None,
               )




stat_py_simulate_dist = Factory(name='Simulate(distribution)',
                description='Generation of a random sample from a distribution.',
                category='STAT.FrequencyDistribution',
                nodemodule='sequence_analysis_wralea.stat',
                nodeclass='py_simulate_dist',
                inputs=None,
                outputs=None,
                widgetmodule=None,
                widgetclass=None,
               )




stat_py_dist_poisson = Factory(name='Poisson',
                description='Construction of a poisson distribution.',
                category='STAT.FrequencyDistribution',
                nodemodule='sequence_analysis_wralea.stat',
                nodeclass='py_dist_poisson',
                inputs=None,
                outputs=None,
                widgetmodule=None,
                widgetclass=None,
               )




stat_py_semimarkov = Factory(name='SemiMarkov',
                description='Construction of a semi-Markov chain from an ASCII file.',
                category='STAT.Sequence',
                nodemodule='sequence_analysis_wralea.stat',
                nodeclass='py_semimarkov',
                inputs=({'interface': IFileStr, 'name': 'filename'}, {'interface': IInt, 'name': 'length', 'value': 20}, {'interface': IBool, 'name': 'counting', 'value': True}),
                outputs=None,
                widgetmodule=None,
                widgetclass=None,
               )




stat_py_segmentation_sample = Factory(name='SegmentationSample',
                description='Segmentation of a sample of sequences.',
                category='STAT.Sequence',
                nodemodule='sequence_analysis_wralea.stat',
                nodeclass='py_segmentation_sample',
                inputs=[{'name': 'seq'}, {'interface': ISequence, 'name': 'nb_segment'}, {'interface': IEnumStr(enum=['Multinomial', 'Ordinal', 'Poisson', 'Gaussian', 'Mean', 'Variance', 'MeanVariance']), 'name': 'available models'}, {'interface': IStr, 'name': 'model'}, {'interface': IEnumStr(enum=['Sequence', 'Residual']), 'name': 'Output', 'value': 'Sequence'}],
                outputs=None,
                widgetmodule=None,
                widgetclass=None,
               )




stat_py_shiftn = Factory(name='Shift multivariate',
                description='Shifting of values.',
                category='STAT',
                nodemodule='sequence_analysis_wralea.stat',
                nodeclass='py_shiftn',
                inputs=None,
                outputs=None,
                widgetmodule=None,
                widgetclass=None,
               )




stat_py_pointwise_average = Factory(name='PointwiseAverage',
                description='Pointwise averaging of a sequence.',
                category='STAT.Sequence',
                nodemodule='sequence_analysis_wralea.stat',
                nodeclass='py_pointwise_average',
                inputs=[{'name': 'seq'}, {'interface': IBool, 'name': 'StandardDeviation', 'value': False}, {'interface': IEnumStr(enum=['Sequence', 'Residual', 'StandardizedResidual']), 'name': 'Output', 'value': 'Sequence'}, {'interface': IDirStr, 'name': 'DirName'}, {'interface': IStr, 'name': 'FileName'}, {'interface': IEnumStr(enum=['ASCII', 'SpreadSheet']), 'name': 'Format', 'value': 'ASCII'}],
                outputs=None,
                widgetmodule=None,
                widgetclass=None,
               )




stat_py_estimate_mixture = Factory(name='EstimateMixture',
                description='Estimation of a finite mixture.',
                category='STAT.FrequencyDistribution',
                nodemodule='sequence_analysis_wralea.stat',
                nodeclass='py_estimate_mixture',
                inputs=({'name': 'histo'}, {'interface': IStr, 'name': 'components'}, {'interface': IInt(min=0, max=1, step=1), 'name': 'MinInfBound', 'value': 0}, {'interface': IEnumStr(enum=['Free', 'Fixed']), 'name': 'InfBoundStatus', 'value': 'Free'}, {'interface': IEnumStr(enum=['Free', 'Fixed']), 'name': 'DistInfBoundStatus', 'value': 'Free'}, {'interface': IEnumStr(enum=['Fixed', 'Estimated']), 'name': 'NbComponents', 'value': 'Fixed'}, {'interface': IEnumStr(enum=['AIC', 'AICc', 'BIC', 'BICc']), 'name': 'Penalty', 'value': 'BICc'}),
                outputs=None,
                widgetmodule=None,
                widgetclass=None,
               )




stat_PyObjectFromFile = Factory(name='Load_from_file',
                description='Construction of a statistical object from an ASCII file.',
                category='STAT',
                nodemodule='sequence_analysis_wralea.stat',
                nodeclass='PyObjectFromFile',
                inputs=None,
                outputs=None,
                widgetmodule=None,
                widgetclass=None,
               )




stat_py_compute_correlation = Factory(name='ComputeCorrelation',
                description='Computation of sample autocorrelation function.',
                category='STAT.Sequence',
                nodemodule='sequence_analysis_wralea.stat',
                nodeclass='py_compute_correlation',
                inputs=[{'name': 'seq'}, {'interface': IInt, 'name': 'MaxLag', 'value': -1}, {'interface': IEnumStr(enum=['Pearson', 'Spearman', 'Kendall']), 'name': 'Type', 'value': 'Pearson'}, {'interface': IEnumStr(enum=['Approximated', 'Exact']), 'name': 'Normalization', 'value': 'Exact'}],
                outputs=None,
                widgetmodule=None,
                widgetclass=None,
               )




stat_py_vectors = Factory(name='Vectors',
                description='Construction of a set of vectors from a multidimensional array, or from a set of sequences.',
                category='STAT',
                nodemodule='sequence_analysis_wralea.stat',
                nodeclass='py_vectors',
                inputs=None,
                outputs=None,
                widgetmodule=None,
                widgetclass=None,
               )




stat_py_cluster = Factory(name='Cluster',
                description='Clustering of values',
                category='STAT',
                nodemodule='sequence_analysis_wralea.stat',
                nodeclass='py_cluster',
                inputs=({'name': 'obj'}, {'interface': IEnumStr(enum=['Step', 'Information', 'Limit']), 'name': 'mode', 'value': 'Step'}, {'interface': IInt(min=1, max=16777216, step=1), 'name': 'variable_rank', 'value': None}, {'interface': IInt, 'name': 'step'}, {'interface': IFloat, 'name': 'information_ratio'}, {'interface': ISequence, 'name': 'limits'}),
                outputs=None,
                widgetmodule=None,
                widgetclass=None,
               )




stat_py_select_variable = Factory(name='SelectVariable',
                description='Selection of variables.',
                category='STAT',
                nodemodule='sequence_analysis_wralea.stat',
                nodeclass='py_select_variable',
                inputs=[{'name': 'multivariate'}, {'interface': ISequence, 'name': 'variables'}, {'interface': IEnumStr(enum=['Keep', 'Reject']), 'name': 'mode', 'value': 'Keep'}],
                outputs=None,
                widgetmodule=None,
                widgetclass=None,
               )




stat_py_dist_binomial = Factory(name='Binomial',
                description='Construction of a binomial distribution.',
                category='STAT.FrequencyDistribution',
                nodemodule='sequence_analysis_wralea.stat',
                nodeclass='py_dist_binomial',
                inputs=None,
                outputs=None,
                widgetmodule=None,
                widgetclass=None,
               )




stat_py_sequences = Factory(name='Sequences',
                description='',
                category='STAT.Sequence',
                nodemodule='sequence_analysis_wralea.stat',
                nodeclass='py_sequences',
                inputs=({'interface': ISequence, 'name': 'seq', 'value': []}, {'interface': ISequence, 'name': 'identifiers'}, {'interface': IEnumStr(enum=['Position', 'Time']), 'name': 'indexParameter', 'value': 'Position'}),
                outputs=({'interface': ISequence, 'name': 'sequences'},),
                widgetmodule=None,
                widgetclass=None,
               )




stat_py_compute_correlation_mult = Factory(name='ComputeCorrelation multivariate',
                description='Computation of sample autocorrelation  or cross-correlation functions.',
                category='STAT.Sequence',
                nodemodule='sequence_analysis_wralea.stat',
                nodeclass='py_compute_correlation_mult',
                inputs=[{'name': 'seq'}, {'interface': IBool, 'name': 'autocorrelation', 'value': True}, {'interface': IInt(min=1, max=16777216, step=1), 'name': 'variable1'}, {'interface': IInt(min=1, max=16777216, step=1), 'name': 'variable2'}, {'interface': IInt, 'name': 'MaxLag', 'value': -1}, {'interface': IEnumStr(enum=['Pearson', 'Spearman', 'Kendall']), 'name': 'Type', 'value': 'Pearson'}, {'interface': IEnumStr(enum=['Approximated', 'Exact']), 'name': 'Normalization', 'value': 'Exact'}],
                outputs=None,
                widgetmodule=None,
                widgetclass=None,
               )




stat_py_to_distribution = Factory(name='ToDistribution',
                description='Cast an object of type DistributionData into a Distribution.',
                category='STAT.Sequence',
                nodemodule='sequence_analysis_wralea.stat',
                nodeclass='py_to_distribution',
                inputs=None,
                outputs=None,
                widgetmodule=None,
                widgetclass=None,
               )




stat_py_hiddensemimarkov = Factory(name='HiddenSemiMarkov',
                description='Construction of an object of type hidden_semi-markov from an ASCII file.',
                category='STAT.Sequence',
                nodemodule='sequence_analysis_wralea.stat',
                nodeclass='py_hiddensemimarkov',
                inputs=({'interface': IFileStr, 'name': 'filename'}, {'interface': IInt, 'name': 'length', 'value': 20}, {'interface': IBool, 'name': 'counting', 'value': True}),
                outputs=None,
                widgetmodule=None,
                widgetclass=None,
               )




stat_py_to_histogram = Factory(name='ToHistogram',
                description='Cast an object of type Distribution into a DistributionData.',
                category='STAT.Sequence',
                nodemodule='sequence_analysis_wralea.stat',
                nodeclass='py_to_histogram',
                inputs=None,
                outputs=None,
                widgetmodule=None,
                widgetclass=None,
               )




stat_py_compound = Factory(name='Compound',
                description='Construction of an object of type compound from a sum     distribution and an elementary distribution.',
                category='STAT.FrequencyDistribution',
                nodemodule='sequence_analysis_wralea.stat',
                nodeclass='py_compound',
                inputs=None,
                outputs=None,
                widgetmodule=None,
                widgetclass=None,
               )




stat_py_regression = Factory(name='Regression',
                description='Simple Regression with a single explanatory variable.',
                category='STAT',
                nodemodule='sequence_analysis_wralea.stat',
                nodeclass='py_regression',
                inputs=[{'name': 'vec'}, {'interface': IEnumStr(enum=['Linear', 'MovingAverage', 'NearestNeighbors']), 'name': 'RegressionModel', 'value': 'Linear'}, {'interface': IInt, 'name': 'explanatoryVariable'}, {'interface': IInt, 'name': 'responseVariable'}, {'interface': ISequence, 'name': 'filter'}, {'interface': ISequence, 'name': 'frequencies'}, {'name': 'distribution'}, {'interface': IFloat, 'name': 'span'}, {'interface': IEnumStr(enum=['Averaging', 'LeastSquares']), 'name': 'Algorithm', 'value': 'Averaging'}, {'interface': IBool, 'name': 'Weighting', 'value': True}],
                outputs=None,
                widgetmodule=None,
                widgetclass=None,
               )




stat_py_plot_segprofile = Factory(name='PlotSegmentProfile (gnuplot)',
                description='Graphical output of segment profiles.',
                category='STAT.Sequence',
                nodemodule='sequence_analysis_wralea.stat',
                nodeclass='py_plot_segprofile',
                inputs=({'name': 'sequence'}, {'interface': IInt, 'name': 'individual'}, {'interface': IInt, 'name': 'nb_segment'}, {'interface': IStr, 'name': 'model'}, {'interface': IEnumStr(enum=['Segment', 'ChangePoint']), 'name': 'output', 'value': 'Segment'}),
                outputs=(),
                widgetmodule=None,
                widgetclass=None,
               )




stat_PyRenewal = Factory(name='Renewal',
                description='Construction of a (either ordinary or equilibrium) renewal process from an inter-event distribution.',
                category='STAT.Renewal',
                nodemodule='sequence_analysis_wralea.stat',
                nodeclass='PyRenewal',
                inputs=None,
                outputs=None,
                widgetmodule=None,
                widgetclass=None,
               )




stat_py_merge = Factory(name='Merge',
                description='Merge of objects of the same data type or merging of sample correlation functions.',
                category='STAT',
                nodemodule='sequence_analysis_wralea.stat',
                nodeclass='py_merge',
                inputs=({'interface': ISequence, 'name': 'data'},),
                outputs=None,
                widgetmodule=None,
                widgetclass=None,
               )




stat_py_plot = Factory(name='Plot (stat_tool)',
                description='Graphical output of an object of the STAT module using either GNUPLOT or matplotlib software.',
                category='Unclassified',
                nodemodule='sequence_analysis_wralea.stat',
                nodeclass='py_plot',
                inputs=({'name': 'obj'}, {'interface': IInt, 'name': 'fig_id'}, {'interface': IStr, 'name': 'title', 'value': ''}, {'interface': IEnumStr(enum=['Data', 'Survival', 'v']), 'name': 'viewpoint', 'value': 'v'}),
                outputs=(),
                widgetmodule=None,
                widgetclass=None,
               )




data_Selector = Factory(name='Selector',
                description='',
                category='Unclassified',
                nodemodule='data',
                nodeclass='Selector',
                inputs=[{'interface': IEnumStr, 'name': 'filename', 'value': None, 'desc': ''}],
                outputs=None,
                widgetmodule=None,
                widgetclass=None,
               )




stat_py_comparisontest = Factory(name='ComparisonTest',
                description='Test of comparison of frequency distributions.',
                category='STAT.FrequencyDistribution',
                nodemodule='sequence_analysis_wralea.stat',
                nodeclass='py_comparisontest',
                inputs=({'interface': IEnumStr(enum=['F', 'T', 'W']), 'name': 'Type', 'value': 'F'}, {'name': 'histo1'}, {'name': 'histo2'}),
                outputs=None,
                widgetmodule=None,
                widgetclass=None,
               )




stat_py_estimate_compound = Factory(name='EstimateCompound',
                description='Estimation of a compound of discrete distributions.',
                category='STAT.FrequencyDistribution',
                nodemodule='sequence_analysis_wralea.stat',
                nodeclass='py_estimate_compound',
                inputs=({'name': 'histo'}, {'name': 'dist'}, {'interface': IEnumStr(enum=['Sum', 'Elementary']), 'name': 'type'}, {'interface': IEnumStr(enum=['Parametric', 'PenalizedLikelihood', 'Zero', 'Likelihood', 'ParametricRegularization']), 'name': 'Estimator', 'value': 'Likelihood'}, {'name': 'InitialDistribution'}, {'interface': IInt(min=0, max=1, step=1), 'name': 'MinInfBound', 'value': 0}, {'interface': IInt, 'name': 'NbIteration', 'value': -1}, {'interface': IEnumStr(enum=['FirstDifference', 'Entropy', 'SecondDifference']), 'name': 'Penalty', 'value': 'SecondDifference'}, {'interface': IFloat, 'name': 'Weight', 'value': -1.0}, {'interface': IEnumStr(enum=['Continuation', 'Zero']), 'name': 'Outside', 'value': 'Zero'}),
                outputs=None,
                widgetmodule=None,
                widgetclass=None,
               )




stat_py_cumulate = Factory(name='Cumulate',
                description='Sum of successive values along sequences.',
                category='STAT.Sequence',
                nodemodule='sequence_analysis_wralea.stat',
                nodeclass='py_cumulate',
                inputs=None,
                outputs=None,
                widgetmodule=None,
                widgetclass=None,
               )




stat_py_estimate_dist = Factory(name='EstimateDistribution',
                description='Estimation of a discrete distribution.',
                category='STAT.FrequencyDistribution',
                nodemodule='sequence_analysis_wralea.stat',
                nodeclass='py_estimate_dist',
                inputs=({'name': 'histo'}, {'interface': IEnumStr(enum=['NONPARAMETRIC', 'BINOMIAL', 'POISSON', 'NEGATIVE_BINOMIAL']), 'name': 'distribution', 'value': 'NONPARAMETRIC'}, {'interface': IInt(min=0, max=1, step=1), 'name': 'MinInfBound', 'value': 0}, {'interface': IEnumStr(enum=['Free', 'Fixed']), 'name': 'InfBoundStatus', 'value': 'Free'}),
                outputs=None,
                widgetmodule=None,
                widgetclass=None,
               )




stat_py_extract_distribution = Factory(name='ExtractDistribution',
                description='Extraction of the distribution part of an object of type model.',
                category='STAT',
                nodemodule='sequence_analysis_wralea.stat',
                nodeclass='py_extract_distribution',
                inputs=[{'name': 'obj'}, {'interface': IEnumStr(enum=['Sum', 'Elementary', 'Component', 'Weight', 'Mixture', 'Convolution', 'Compound']), 'name': 'key', 'value': 'Sum'}, {'interface': IFloat, 'name': 'index', 'value': None}],
                outputs=None,
                widgetmodule=None,
                widgetclass=None,
               )




stat_py_shift = Factory(name='Shift',
                description='Shifting of values.',
                category='STAT',
                nodemodule='sequence_analysis_wralea.stat',
                nodeclass='py_shift',
                inputs=None,
                outputs=None,
                widgetmodule=None,
                widgetclass=None,
               )




stat_py_dist_negativebinomial = Factory(name='NegativeBinomial',
                description='Construction of a negative binomial distribution.',
                category='STAT.FrequencyDistribution',
                nodemodule='sequence_analysis_wralea.stat',
                nodeclass='py_dist_negativebinomial',
                inputs=None,
                outputs=None,
                widgetmodule=None,
                widgetclass=None,
               )




stat_py_markov = Factory(name='Markov',
                description='',
                category='STAT.Sequence',
                nodemodule='sequence_analysis_wralea.stat',
                nodeclass='py_markov',
                inputs=({'interface': IFileStr, 'name': 'filename'}, {'interface': IInt, 'name': 'length', 'value': 20}),
                outputs=None,
                widgetmodule=None,
                widgetclass=None,
               )




stat_py_dist_uniform = Factory(name='Uniform',
                description='Construction of a uniform distribution.',
                category='STAT.FrequencyDistribution',
                nodemodule='sequence_analysis_wralea.stat',
                nodeclass='py_dist_uniform',
                inputs=None,
                outputs=None,
                widgetmodule=None,
                widgetclass=None,
               )




stat_py_fit = Factory(name='Fit',
                authors='Thomas Cokelaer',
                description='Fit distribution to histogram.',
                category='STAT',
                nodemodule='sequence_analysis_wralea.stat',
                nodeclass='py_fit',
                inputs=None,
                outputs=None,
                widgetmodule=None,
                widgetclass=None,
               )




stat_py_segmentation = Factory(name='Segmentation',
                description='Segmentation of a sequence.',
                category='STAT.Sequence',
                nodemodule='sequence_analysis_wralea.stat',
                nodeclass='py_segmentation',
                inputs=[{'name': 'seq'}, {'interface': IInt, 'name': 'individual'}, {'interface': IInt, 'name': 'nb_segment'}, {'interface': ISequence, 'name': 'change points'}, {'interface': IEnumStr(enum=['Multinomial', 'Ordinal', 'Poisson', 'Gaussian', 'Mean', 'Variance', 'MeanVariance']), 'name': 'available models'}, {'interface': IStr, 'name': 'model'}, {'interface': IEnumStr(enum=['Fixed', 'Estimated']), 'name': 'NbSegment', 'value': 'Estimated'}, {'interface': IEnumStr(enum=['Sequence', 'Residual']), 'name': 'Output', 'value': 'Sequence'}],
                outputs=None,
                widgetmodule=None,
                widgetclass=None,
               )




stat_py_select_individual = Factory(name='SelectIndividual',
                description='Selection of individuals.',
                category='STAT',
                nodemodule='sequence_analysis_wralea.stat',
                nodeclass='py_select_individual',
                inputs=[{'name': 'data'}, {'interface': ISequence, 'name': 'individuals'}, {'interface': IEnumStr(enum=['Keep', 'Reject']), 'name': 'mode', 'value': 'Keep'}],
                outputs=None,
                widgetmodule=None,
                widgetclass=None,
               )




stat_py_hiddenmarkov = Factory(name='HiddenMarkov',
                description='Construction of an object of type hidden_markov from an ASCII file.',
                category='STAT.Sequence',
                nodemodule='sequence_analysis_wralea.stat',
                nodeclass='py_hiddenmarkov',
                inputs=({'interface': IFileStr, 'name': 'filename'}, {'interface': IInt, 'name': 'length', 'value': 20}),
                outputs=None,
                widgetmodule=None,
                widgetclass=None,
               )




stat_py_extract_histogram = Factory(name='ExtractHistogram',
                description='Extraction of the histogram part of an object of type model.',
                category='STAT',
                nodemodule='sequence_analysis_wralea.stat',
                nodeclass='py_extract_histogram',
                inputs=[{'name': 'obj'}, {'interface': IEnumStr(enum=['Sum', 'Elementary']), 'name': 'choice', 'value': 'Sum'}],
                outputs=None,
                widgetmodule=None,
                widgetclass=None,
               )




stat_py_compare_vectors = Factory(name='Compare(Vectors)',
                description='Comparison of vectors.',
                category='STAT.Vectors',
                nodemodule='sequence_analysis_wralea.stat',
                nodeclass='py_compare_vectors',
                inputs=None,
                outputs=None,
                widgetmodule=None,
                widgetclass=None,
               )




stat_py_histogram = Factory(name='Histogram',
                description='Histogram from a sequence of int.',
                category='STAT.FrequencyDistribution',
                nodemodule='sequence_analysis_wralea.stat',
                nodeclass='py_histogram',
                inputs=[{'interface': None, 'name': 'seq', 'value': None}],
                outputs=({'interface': None, 'name': 'out'},),
                widgetmodule=None,
                widgetclass=None,
               )




stat_PyRenewalAscii = Factory(name='Renewal_from_file',
                description='Construction of a (either ordinary or equilibrium) renewal process from an inter-event distribution.',
                category='STAT.Renewal',
                nodemodule='sequence_analysis_wralea.stat',
                nodeclass='PyRenewalAscii',
                inputs=None,
                outputs=None,
                widgetmodule=None,
                widgetclass=None,
               )




stat_py_extract_data = Factory(name='ExtractData',
                description='Extraction of the data part of an object of type model.',
                category='STAT.Sequence',
                nodemodule='sequence_analysis_wralea.stat',
                nodeclass='py_extract_data',
                inputs=None,
                outputs=None,
                widgetmodule=None,
                widgetclass=None,
               )




stat_py_compare_frequency = Factory(name='Compare(frequency)',
                description='Comparison of frequency distributions.',
                category='STAT.FrequencyDistribution',
                nodemodule='sequence_analysis_wralea.stat',
                nodeclass='py_compare_frequency',
                inputs=({'interface': ISequence, 'name': 'histos'}, {'interface': IEnumStr(enum=['Numeric', 'Ordinal', 'Symbolic']), 'name': 'Type', 'value': 'Numeric'}),
                outputs=None,
                widgetmodule=None,
                widgetclass=None,
               )




stat_py_estimate_conv = Factory(name='EstimateConvolution',
                description='Estimation of a convolution of discrete distributions.',
                category='STAT.FrequencyDistribution',
                nodemodule='sequence_analysis_wralea.stat',
                nodeclass='py_estimate_conv',
                inputs=({'name': 'histo'}, {'name': 'dist'}, {'interface': IEnumStr(enum=['Likelihood', 'PenalizedLikelihood', 'Parametric']), 'name': 'Estimator', 'value': 'Likelihood'}, {'name': 'InitialDistribution'}, {'interface': IInt(min=0, max=1, step=1), 'name': 'MinInfBound', 'value': 0}, {'interface': IInt, 'name': 'NbIteration', 'value': -1}, {'interface': IEnumStr(enum=['FirstDifference', 'SecondDifference', 'Entropy']), 'name': 'Penalty', 'value': 'SecondDifference'}, {'interface': IFloat, 'name': 'Weight', 'value': -1.0}, {'interface': IEnumStr(enum=['Zero', 'Continuation']), 'name': 'Outside', 'value': 'Zero'}),
                outputs=None,
                widgetmodule=None,
                widgetclass=None,
               )




stat_py_extract_histogram_sequences = Factory(name='ExtractHistogram(Sequences)',
                description='Extraction of the histogram part of an object of type model.',
                category='STAT',
                nodemodule='sequence_analysis_wralea.stat',
                nodeclass='py_extract_histogram_sequences',
                inputs=[{'name': 'obj'}, {'interface': IEnumStr(enum=['Sequences', 'Vectors']), 'name': 'input data', 'value': None}, {'interface': IEnumStr(enum=['Value', 'Length']), 'name': 'choice', 'value': 'Value'}, {'interface': IInt, 'name': 'variable', 'value': 1}],
                outputs=None,
                widgetmodule=None,
                widgetclass=None,
               )




