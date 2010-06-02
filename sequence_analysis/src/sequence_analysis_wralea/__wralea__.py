# -*- python -*-
#
#       aml package: AMAPmod package interface for the amlPy module
#
#       Copyright or (C) or Copr. 2006 INRIA - CIRAD - INRA  
#
#       File author(s): Christophe Pradal <christophe.prada@cirad.fr>
#
#       Distributed under the Cecill-C License.
#       See accompanying file LICENSE.txt or copy at
#           http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html
# 
#       OpenAlea WebSite : http://openalea.gforge.inria.fr
#


__doc__ = """
Wralea for AMAPmod.
"""

__license__ = "Cecill"
__revision__ = " $Id$"


from openalea.core import Factory, IFileStr, IInt, IBool, IFloat, \
    ISequence, IEnumStr, IStr, IDirStr
#from openalea.core.external import *


import stat 

__openalea_url = 'http://openalea.gforge.inria.fr'

__name__ = "vplants.stats"
__version__ = '0.7.0'
__license__ = 'CECILL-C'
__authors__ = 'OpenAlea consortium'
__institutes__ = 'INRIA/CIRAD'
__description__ = 'stat_tool nodes.'
__url__ = __openalea_url + '/doc/vplants/stat_tool/doc/_build/html/contents.html'
__editable__ = 'True'
__icon__ = 'icon.png'


__all__ = ['binomial',
           'cluster',
           'compare',
           'compare_vectors',
           'comparison_test',
           'compound',
           'compute_correlation',
           'compute_correlation_multivariate',
           'convolution',
           'cumulate',
           'estimate_compound',
           'estimate_convolution',
           'estimate_distribution',
           'estimate_mixture',
           'extract_data',
           'extract_distribution',
           'extract_histogram',
           'hidden_markov',
           'hidden_semi_markov',
           'histogram',
           'keep_type',
           'load_from_file',
           'markov',
           'merge',
           'negative_binomial',
           'plot',
           'plot_segment_profile',
           'point_wise_average',
           'poisson',
           'regression',
           'renewal',
           'renewal_from_file',
           'segmentation_sample',
           'segmentation',
           'select_individual',
           'select_variable',
           'semi_markov',
           'sequences',
           'shift',
           'shift_multivariate',
           'simulate',
           'stat',
           'to_distribution',
           'to_histogram',
           'uniform',
           'vectors']



# Input/output functions

compound = Factory(name="Compound",
             description="Construction of an object of type compound from a sum \
    distribution and an elementary distribution.",
             category="STAT.FrequencyDistribution",
             nodemodule="stat",
             nodeclass="py_compound",
             )



load_from_file = Factory(name="Load_from_file",
              description="Construction of a statistical object from an ASCII file.",
              category="STAT",
              nodemodule="stat",
              nodeclass="PyObjectFromFile",
              )



convolution = Factory(name="Convolution",

              description="Construction of an object of type convolution from an elementary distribution.",
              category="STAT.FrequencyDistribution",
              nodemodule="stat",
              nodeclass="py_convolution",
              )


binomial = Factory(name="Binomial",
              description="Construction of a binomial distribution.",
              category="STAT.FrequencyDistribution",
              nodemodule="stat",
              nodeclass="py_dist_binomial",
              )

#__all__.append(nf.name.lower())

poisson = Factory(name="Poisson",
              description="Construction of a poisson distribution.",
              category="STAT.FrequencyDistribution",
              nodemodule="stat",
              nodeclass="py_dist_poisson",
              )
#__all__.append(nf.name.lower())


negative_binomial = Factory(name="NegativeBinomial",
              description="Construction of a negative binomial distribution.",
              category="STAT.FrequencyDistribution",
              nodemodule="stat",
              nodeclass="py_dist_negativebinomial",
              )
#__all__.append(nf.name.lower())


uniform = Factory(name="Uniform",
              description="Construction of a uniform distribution.",
              category="STAT.FrequencyDistribution",
              nodemodule="stat",
              nodeclass="py_dist_uniform",
              )
#__all__.append('nf')

hidden_markov = Factory(name="HiddenMarkov",
              description="Construction of an object of type hidden_markov from an ASCII file.",
              category="STAT.Sequence",
              nodemodule="stat",
              nodeclass="py_hiddenmarkov",
              inputs=(dict(name='filename', interface=IFileStr),
                      dict(name='length', interface=IInt, value=20)),
              )
#__all__.append('nf')

hidden_semi_markov = Factory(name="HiddenSemiMarkov",
              description="Construction of an object of type hidden_semi-markov from an ASCII file.",
              category="STAT.Sequence",
              nodemodule="stat",
              nodeclass="py_hiddensemimarkov",
              inputs=(dict(name='filename', interface=IFileStr),
                      dict(name='length', interface=IInt, value=20),
                      dict(name='counting', interface=IBool, value=True)),
              )
#__all__.append('nf')

markov = Factory(name="Markov",
              description="",
              category="STAT.Sequence",
              nodemodule="stat",
              nodeclass="py_markov",
              inputs=(dict(name='filename', interface=IFileStr),
                      dict(name='length', interface=IInt, value=20)),
              )
#__all__.append('nf')

# TODO Mixture with a set of dist and proba

'''
'''

renewal = Factory(name="Renewal",
              description="Construction of a (either ordinary or equilibrium) renewal process from an inter-event distribution.",
              category="STAT.Renewal",
              nodemodule="stat",
              nodeclass="PyRenewal",
              )
#__all__.append('nf')

renewal_from_file = Factory(name="Renewal_from_file",
              description="Construction of a (either ordinary or equilibrium) renewal process from an inter-event distribution.",
              category="STAT.Renewal",
              nodemodule="stat",
              nodeclass="PyRenewalAscii",
              )
#__all__.append('nf')

semi_markov = Factory(name="SemiMarkov",
              description="Construction of a semi-Markov chain from an ASCII file.",
              category="STAT.Sequence",
              nodemodule="stat",
              nodeclass="py_semimarkov",
              inputs=(dict(name='filename', interface=IFileStr),
                      dict(name='length', interface=IInt, value=20),
                      dict(name='counting', interface=IBool, value=True)),
              )
#__all__.append('nf')
sequences = Factory(name="Sequences",
              description="",
              category="STAT.Sequence",
              nodemodule="stat",
              nodeclass="py_sequences",
              inputs=(dict(name='seq', interface=ISequence, value=[]),
                      dict(name="identifiers", interface=ISequence),
                      dict(name="indexParameter",
                           interface=IEnumStr(["Position", "Time"]),
                           value="Position")),
              outputs=(dict(name="sequences", interface=ISequence),)
              )
#__all__.append('nf')

'''
nf = Factory( name= "TimeEvents", 
              description= "", 
              category = "STAT", 
              nodemodule = "stat",
              nodeclass = "py_timeevents",
              )
package.add_factory( nf )

nf = Factory( name= "TopParameters", 
              description= "", 
              category = "STAT", 
              nodemodule = "stat",
              nodeclass = "py_topparameters",
              )
package.add_factory( nf )

nf = Factory( name= "Tops", 
              description= "", 
              category = "STAT", 
              nodemodule = "stat",
              nodeclass = "py_tops",
              )
package.add_factory( nf )

nf = Factory( name= "VectorDistance", 
              description= "", 
              category = "STAT", 
              nodemodule = "stat",
              nodeclass = "py_vectordistance",
              )
package.add_factory( nf )


nf = Factory( name= "", 
              description= "", 
              category = "STAT", 
              nodemodule = "stat",
              nodeclass = "py_",
              )
package.add_factory( nf )

'''
histogram = Factory(name="Histogram",
              description="Histogram from a sequence of int.",
              category="STAT.FrequencyDistribution",
              nodemodule="stat",
              nodeclass="py_histogram",
              inputs=(dict(name='sequence', interface=ISequence, value=[]),),
              )
#__all__.append('nf')




#function = 'Plot'
#cmd = 'from openalea.stat_tool import %s' % function
#exec(cmd)
#description = getattr(eval(function), '__doc__')

plot = Factory(name="Plot (stat_tool)",
              description="Graphical output of an object of the STAT module using either GNUPLOT or matplotlib software.",
              category="STAT",
              nodemodule="stat",
              nodeclass="py_plot",
              inputs=(dict(name='obj'),
                        dict(name='fig_id', interface=IInt),
                        dict(name='title', interface=IStr, value=''),
                        ),
              outputs=(),
              )
#__all__.append('nf')

plot_segment_profile = Factory(name="PlotSegmentProfile (gnuplot)",
              description="Graphical output of segment profiles.",
              category="STAT.Sequence",
              nodemodule="stat",
              nodeclass="py_plot_segprofile",
              inputs=(dict(name='sequence'),
                         dict(name='individual', interface=IInt),
                         dict(name='nb_segment', interface=IInt),
                         dict(name='model', interface=IStr),
                         dict(name='output', interface=IEnumStr(['Segment', 'ChangePoint']), value='Segment'),
                     ),
              outputs=(),)
#__all__.append('nf')



#//////////////////////////////////////////////////////////////////////////////
# Compare algorithms
compare = Factory(name="Compare(frequency)",
              description="Comparison of frequency distributions.",
              category="STAT.FrequencyDistribution",
              nodemodule="stat",
              nodeclass="py_compare_frequency",
              inputs=(dict(name='histos', interface=ISequence),
                      dict(name="Type",
                           interface=IEnumStr(["Numeric", "Ordinal", "Symbolic"]),
                           value="Numeric"),
                     ),
              )
#__all__.append('nf')

compare_vectors = Factory(name="Compare(Vectors)",
              description="Comparison of vectors.",
              category="STAT.Vectors",
              nodemodule="stat",
              nodeclass="py_compare_vectors",
              )
#__all__.append('nf')

#//////////////////////////////////////////////////////////////////////////////

comparison_test = Factory(name="ComparisonTest",
              description="Test of comparison of frequency distributions.",
              category="STAT.FrequencyDistribution",
              nodemodule="stat",
              nodeclass="py_comparisontest",
              inputs=(dict(name="Type",
                           interface=IEnumStr(["F", "T", "W"]),
                           value="F"),
                      dict(name='histo1'),
                      dict(name='histo2')),
              )
#__all__.append('nf')

#//////////////////////////////////////////////////////////////////////////////

estimate_distribution = Factory(name="EstimateDistribution",
              description="Estimation of a discrete distribution.",
              category="STAT.FrequencyDistribution",
              nodemodule="stat",
              nodeclass="py_estimate_dist",
              inputs=(
                      dict(name='histo'),
                      dict(name="distribution",
                           interface=IEnumStr(["NONPARAMETRIC",
                                                "BINOMIAL",
                                                "POISSON",
                                                "NEGATIVE_BINOMIAL"]),
                           value="NONPARAMETRIC"),
                      dict(name='MinInfBound', interface=IInt(0, 1), value=0),
                      dict(name='InfBoundStatus', interface=IEnumStr(['Free', 'Fixed']), value="Free"),
                     )
              )
#__all__.append('nf')

estimate_mixture = Factory(name="EstimateMixture",
              description="Estimation of a finite mixture.",
              category="STAT.FrequencyDistribution",
              nodemodule="stat",
              nodeclass="py_estimate_mixture",
              inputs=(
                      dict(name='histo'),
                      dict(name='components', interface=IStr),
                      dict(name='MinInfBound', interface=IInt(0, 1), value=0),
                      dict(name='InfBoundStatus', interface=IEnumStr(['Free', 'Fixed']), value="Free"),
                      dict(name='DistInfBoundStatus', interface=IEnumStr(['Free', 'Fixed']), value="Free"),
                      dict(name='NbComponents', interface=IEnumStr(['Fixed', 'Estimated']), value="Fixed"),
                      dict(name='Penalty', interface=IEnumStr(['AIC', 'AICc', 'BIC', 'BICc']), value='BICc'),
                     )
              )
#__all__.append('nf')

estimate_convolution = Factory(name="EstimateConvolution",
              description="Estimation of a convolution of discrete distributions.",
              category="STAT.FrequencyDistribution",
              nodemodule="stat",
              nodeclass="py_estimate_conv",
              inputs=(
                      dict(name='histo'),
                      dict(name='dist'),
                      dict(name='Estimator',
                           interface=IEnumStr(["Likelihood",
                                              "PenalizedLikelihood",
                                              "Parametric"]),
                           value="Likelihood"),
                      dict(name='InitialDistribution'),
                      dict(name='MinInfBound', interface=IInt(0, 1), value=0),
                      dict(name='NbIteration', interface=IInt, value= -1),
                      dict(name='Penalty', interface=IEnumStr(['FirstDifference', 'SecondDifference', 'Entropy']), value='SecondDifference'),
                      dict(name='Weight', interface=IFloat, value= -1.),
                      dict(name='Outside', interface=IEnumStr(['Zero', 'Continuation']), value='Zero'),
                     )
              )
#__all__.append('nf')


cmd = 'from openalea.stat_tool.enums import smoothing_penalty_type, estimator_type, outside_type'
exec(cmd)

estimate_compound = Factory(name="EstimateCompound",
             description="Estimation of a compound of discrete distributions.",
             category="STAT.FrequencyDistribution",
             nodemodule="stat",
             nodeclass="py_estimate_compound",
             inputs=(
                     dict(name='histo'),
                     dict(name='dist'),
                     dict(name='type', interface=IEnumStr(['Sum', 'Elementary'])),
                     dict(name='Estimator', interface=IEnumStr(estimator_type.keys()), value="Likelihood"),
                     dict(name='InitialDistribution'),
                     dict(name='MinInfBound', interface=IInt(0, 1), value=0),
                     dict(name='NbIteration', interface=IInt, value= -1),
                     dict(name='Penalty', interface=IEnumStr(smoothing_penalty_type.keys()), value='SecondDifference'),
                     dict(name='Weight', interface=IFloat, value= -1.),
                     dict(name='Outside', interface=IEnumStr(outside_type.keys()), value='Zero'),
                     )
              )
#__all__.append('nf')

'''    
nf = Factory( name= "Estimate(distribution)", 
              description= "Estimation of distributions.", 
              category = "STAT", 
              nodemodule = "stat",
              nodeclass = "py_estimate_distrib",
              inputs=(dict(name="histo"), 
                      dict(name="distribution", 
                           interface= IEnumStr(["NON-PARAMETRIC",
                                                "BINOMIAL",
                                                "POISSON",
                                                "NEGATIVE_BINOMIAL",
                                                "UNIFORM",
                                                "MIXTURE",
                                                "CONVOLUTION",
                                                "COMPOUND"]),
                           value="NON-PARAMETRIC"),
                      dict(name="mixtures",interface=IStr),
                      dict(name="unknow",interface=IEnumStr(["","Sum","Elementary"]),value=""),
                      dict(name="MinInfBound",interface=IInt(0,1), value = 0),
                      dict(name="InfBoundStatus",interface=IEnumStr(["Free","Fixed"]), value="Free"),
                      dict(name="DistInfBoundStatus",interface=IEnumStr(["Fixed","Estimated"]), value="Fixed"),
                      dict(name="NbComponents",interface=IEnumStr(["Fixed","Estimated"]), value="Fixed"),
                     dict(name="Penalty",interface=IEnumStr(["AIC","AICc","BIC"]), value="AICc"),
                      dict(name="Parametric",interface=IBool, value=True),
                      dict(name="InitialDistribution"),
                     )
            )
package.add_factory( nf )
'''
#//////////////////////////////////////////////////////////////////////////////

merge = Factory(name="Merge",
              description="Merge of objects of the same data type or merging of sample correlation functions.",
              category="STAT",
              nodemodule="stat",
              nodeclass="py_merge",
              inputs=(dict(name='data', interface=ISequence),),
              )
#__all__.append('nf')


extract_histogram = Factory(name="ExtractHistogram",
              description="Extraction of the histogram part of an object of type model.",
              category="STAT",
              nodemodule="stat",
              nodeclass="py_extract_histogram",
              inputs=(
                      [dict(name='obj'),
                      dict(name='choice', interface=IEnumStr(['Sum', 'Elementary']), value='Sum'),
                      ]
                     )
              )
#__all__.append('nf')

extract_distribution = Factory(name="ExtractDistribution",
             description="Extraction of the distribution part of an object of type model.",
             category="STAT",
             nodemodule="stat",
             nodeclass="py_extract_distribution",
             inputs=(
                     [dict(name='obj'),
                      dict(name='key', interface=IEnumStr(
                        ['Sum', 'Elementary', 'Component', 'Weight',
                          'Mixture', 'Convolution', 'Compound']),
                           value='Sum'),
                      dict(name='index', interface=IFloat, value=None),
                      ]
                      )
              )
#__all__.append('nf')




shift = Factory(name="Shift",
              description="Shifting of values.",
              category="STAT",
              nodemodule="stat",
              nodeclass="py_shift",
              )
#__all__.append('nf')

shift_multivariate = Factory(name="Shift multivariate",
              description="Shifting of values.",
              category="STAT",
              nodemodule="stat",
              nodeclass="py_shiftn",
              )
#__all__.append('nf')

cluster = Factory(name="Cluster",
              description="Clustering of values",
              category="STAT",
              nodemodule="stat",
              nodeclass="py_cluster",
              inputs=(dict(name='obj'),
                       dict(name='mode',
                            interface=IEnumStr(["Step", "Information", "Limit"]),
                            value="Step"),
                       dict(name='variable_rank', interface=IInt(min=1), value=None),
                       dict(name='step', interface=IInt),
                       dict(name='information_ratio', interface=IFloat),
                       dict(name='limits', interface=ISequence),
                      ),
              )
#__all__.append('nf')


simulate = Factory(name="Simulate(distribution)",
              description="Generation of a random sample from a distribution.",
              category="STAT.FrequencyDistribution",
              nodemodule="stat",
              nodeclass="py_simulate_dist",
              )
#__all__.append('nf')

select_variable = Factory(name="SelectVariable",
              description="Selection of variables.",
              category="STAT",
              nodemodule="stat",
              nodeclass="py_select_variable",
              inputs=[dict(name='multivariate'),
                        dict(name='variables', interface=ISequence),
                        dict(name='mode', interface=IEnumStr(['Keep', 'Reject']), value='Keep')],
            )
#__all__.append('nf')

select_individual = Factory(name="SelectIndividual",
              description="Selection of individuals.",
              category="STAT",
              nodemodule="stat",
              nodeclass="py_select_individual",
              inputs=[dict(name='data'),
                        dict(name='individuals', interface=ISequence),
                        dict(name='mode', interface=IEnumStr(['Keep', 'Reject']), value='Keep')],
            )
#__all__.append('nf')

segmentation = Factory(name="Segmentation",
              description="Segmentation of a sequence.",
              category="STAT.Sequence",
              nodemodule="stat",
              nodeclass="py_segmentation",
              inputs=[dict(name='seq'),
                        dict(name='individual', interface=IInt),
                        dict(name='nb_segment', interface=IInt),
                        dict(name='change points', interface=ISequence),
                        dict(name='available models',
                              interface=IEnumStr(['Multinomial', 'Ordinal', 'Poisson', 'Gaussian', 'Mean', 'Variance', 'MeanVariance']),),
                        dict(name='model', interface=IStr),
                        dict(name='NbSegment', interface=IEnumStr(['Fixed', 'Estimated']), value='Estimated'),
                        dict(name='Output', interface=IEnumStr(['Sequence', 'Residual']), value='Sequence')],
            )
#__all__.append('nf')

segmentation_sample = Factory(name="SegmentationSample",
              description="Segmentation of a sample of sequences.",
              category="STAT.Sequence",
              nodemodule="stat",
              nodeclass="py_segmentation_sample",
              inputs=[dict(name='seq'),
                        dict(name='nb_segment', interface=ISequence),
                        dict(name='available models', interface=IEnumStr(['Multinomial', 'Ordinal', 'Poisson', 'Gaussian', 'Mean', 'Variance', 'MeanVariance']),),
                        dict(name='model', interface=IStr),
                        dict(name='Output', interface=IEnumStr(['Sequence', 'Residual']), value='Sequence')],
            )
#__all__.append('nf')

compute_correlation = Factory(name="ComputeCorrelation",
              description="Computation of sample autocorrelation function.",
              category="STAT.Sequence",
              nodemodule="stat",
              nodeclass="py_compute_correlation",
              inputs=[dict(name='seq'),
                        dict(name='MaxLag', interface=IInt(), value= -1),
                        dict(name='Type', interface=IEnumStr(['Pearson', 'Spearman', 'Kendall']), value='Pearson'),
                        dict(name='Normalization', interface=IEnumStr(['Approximated', 'Exact']), value='Exact')],
            )
#__all__.append('nf')

compute_correlation_multivariate = Factory(\
    name="ComputeCorrelation multivariate",
    description="Computation of sample autocorrelation  or cross-correlation functions.",
    category="STAT.Sequence",
              nodemodule="stat",
              nodeclass="py_compute_correlation_mult",
              inputs=[dict(name='seq'),
                        dict(name='autocorrelation', interface=IBool, value=True),
                        dict(name='variable1', interface=IInt(min=1)),
                        dict(name='variable2', interface=IInt(min=1)),
                        dict(name='MaxLag', interface=IInt(), value= -1),
                        dict(name='Type', interface=IEnumStr(['Pearson', 'Spearman', 'Kendall']), value='Pearson'),
                        dict(name='Normalization', interface=IEnumStr(['Approximated', 'Exact']), value='Exact')],
            )

#__all__.append('nf')

point_wise_average = Factory(name="PointwiseAverage",
              description="Pointwise averaging of a sequence.",
              category="STAT.Sequence",
              nodemodule="stat",
              nodeclass="py_pointwise_average",
              inputs=[dict(name='seq',),
                            dict(name='StandardDeviation', interface=IBool, value=False),
                        dict(name='Output', interface=IEnumStr(['Sequence', 'Residual', 'StandardizedResidual']), value='Sequence'),
                        dict(name='DirName', interface=IDirStr),
                        dict(name='FileName', interface=IStr),
                        dict(name='Format', interface=IEnumStr(['ASCII', 'SpreadSheet']), value='ASCII'), ],
            )

#__all__.append('nf')

vectors = Factory(name="Vectors",
              description="Construction of a set of vectors from a multidimensional array, or from a set of sequences.",
              category="STAT",
              nodemodule="stat",
              nodeclass="py_vectors",
              )
#__all__.append('nf')

regression = Factory(name="Regression",
              description="Simple Regression with a single explanatory variable.",
              category="STAT",
              nodemodule="stat",
              nodeclass="py_regression",
              inputs=[ dict(name='vec'),
                         dict(name='RegressionModel',
                              interface=IEnumStr([ "Linear", "MovingAverage", "NearestNeighbors"]), value="Linear"),
                         dict(name='explanatoryVariable', interface=IInt),
                         dict(name='responseVariable', interface=IInt),
                         dict(name='filter', interface=ISequence),
                        dict(name='frequencies', interface=ISequence),
                        dict(name='distribution'),
                        dict(name='span', interface=IFloat),
                        dict(name='Algorithm', interface=IEnumStr(['Averaging', 'LeastSquares']), value='Averaging'),
                        dict(name='Weighting', interface=IBool, value=True),
 ],
              )
#__all__.append('nf')

cumulate = Factory(name="Cumulate",
              description="Sum of successive values along sequences.",
              category="STAT.Sequence",
              nodemodule="stat",
              nodeclass="py_cumulate")

#__all__.append('nf')

to_histogram = Factory(name="ToHistogram",
              description="Cast an object of type Distribution into a DistributionData.",
              category="STAT.Sequence",
              nodemodule="stat",
              nodeclass="py_to_histogram")

#__all__.append(nf.name.lower())


to_distribution = Factory(name="ToDistribution",
              description="Cast an object of type DistributionData into a Distribution.",
              category="STAT.Sequence",
              nodemodule="stat",
              nodeclass="py_to_distribution")

#__all__.append(nf.name.lower())



extract_data = Factory(name="ExtractData",
             description="Extraction of the data part of an object of type model.",
             category="STAT.Sequence",
             nodemodule="stat",
             nodeclass="py_extract_data")

#__all__.append(nf.name.lower())


