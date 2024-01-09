"""
#########################################################################
#
#  Objective: analyzing the apical growth and the branching of Corsican pines,
#
#  Methods: extraction of trends and residuals by linear filtering
#           (differencing and weighted moving average), computation of
#           sample autocorrelation functions, change-point detection models and
#           hidden semi-Markov chains.
#
#  70-year-old Corsican pines.
#
#  Data: Celine Meredieu and Yves Caraglio
#
#  INDEX_PARAMETER : TIME (year of growth),
#
#  VARIABLE 1 : length of the annual shoot (cm),
#  VARIABLE 2 : number of  branches per annual shoot / tier.
#
#########################################################################
"""
__revision__ = "$Id: functional3.py 9401 2010-08-10 12:24:59Z cokelaer $"


import os
from openalea.sequence_analysis import *
from openalea.sequence_analysis.estimate import  Estimate
from openalea.sequence_analysis.compare import  Compare

seq69 = Sequences(get_shared_data( "pin_laricio_7x.seq"))
seq70 = Cluster(seq69, "Step", 1, 10)
#seq70 = IndexParameterExtract(Cluster(seq69, "Step", 2, 10), 1927, MaxIndex=1992)

seq2 = SelectVariable(seq70, 1)
Plot(seq2, 2, 5, "Gaussian", ViewPoint="SegmentProfile")

Plot(seq70, ViewPoint="Data")
Plot(Cumulate(seq70), ViewPoint="Data")

vec70 = Vectors(seq70)
Plot(Regression(vec70, "MovingAverage", 1, 2, [1]))
Plot(Regression(vec70, "MovingAverage", 1, 3, [1]))

vec71 =  Vectors(SelectIndividual(seq70, [1, 2, 3]))
Plot(Regression(vec71, "MovingAverage", 1, 2, [1]))
Plot(Regression(vec71, "MovingAverage", 1, 3, [1]))

seq71 = PointwiseAverage(SelectIndividual(seq70, [1, 2, 3]), StandardDeviation=True)
Plot(SelectIndividual(seq71, [0, 1, 2, 3]), ViewPoint="Data")
Plot(SelectIndividual(seq71, [0, 4]), ViewPoint="Data")

seq72 = PointwiseAverage(SelectIndividual(seq70, [1, 2, 3]), Output="Residual")
seq72 = PointwiseAverage(SelectIndividual(seq70, [1, 2, 3]), Output="StandardizedResidual")
Plot(SelectIndividual(seq72, [1, 2, 3]), ViewPoint="Data")
Plot(SelectIndividual(Cumulate(seq72), [1, 2, 3]), ViewPoint="Data")

vec73 =  Vectors(SelectIndividual(seq70, [4, 5, 6]))
Plot(Regression(vec73, "MovingAverage", 1, 2, [1]))
Plot(Regression(vec73, "MovingAverage", 1, 3, [1]))

seq73 = PointwiseAverage(SelectIndividual(seq70, [4, 5, 6]), StandardDeviation=True)
Plot(SelectIndividual(seq73, [3, 4, 5, 6]), ViewPoint="Data")
Plot(SelectIndividual(seq73, [3, 7]), ViewPoint="Data")

seq74 = PointwiseAverage(SelectIndividual(seq70, [4, 5, 6]), Output="Residual")
seq74 = PointwiseAverage(SelectIndividual(seq70, [4, 5, 6]), Output="StandardizedResidual")
Plot(SelectIndividual(seq74, [4, 5, 6]), ViewPoint="Data")
Plot(SelectIndividual(Cumulate(seq74), [4, 5, 6]), ViewPoint="Data")


matrix70 = Compare(seq70, VectorDistance("N", "N"), IndelFactor=1., End="Free")
matrix70 = Compare(seq70, VectorDistance("N", "N"), IndelFactor=1.)
Display(Clustering(matrix70, "Partition", 3))
Clustering(matrix70, "Hierarchy")

# extraction of trends (slowly varying component) and residuals (rapidly varying component)
# by symmetric smoothing filters and computation of sample autocorrelation functions from residuals

seq75 = Difference(seq70)
acf11 = Merge(ComputeCorrelation(seq75, 1, MaxLag=10), ComputeCorrelation(seq75, 2, MaxLag=10))
ComputeWhiteNoiseCorrelation(acf11, 1)
Plot(acf11)

# symmetric smoothing filters of half-width 3

filter1 = Convolution(Distribution("B", 0, 6, 0.2), Distribution("B", 0, 6, 0.8))
filter2 = Convolution(Distribution("B", 0, 4, 0.2), Distribution("B", 0, 4, 0.5), Distribution("B", 0, 4, 0.8))
filter3 = Convolution(Distribution("U", 0, 2), Distribution("U", 0, 2), Distribution("U", 0, 2), Distribution("U", 0, 2), Distribution("U", 0, 2), Distribution("U", 0, 2))
filter4 = Convolution(Distribution("U", 0, 3), Distribution("U", 0, 3), Distribution("U", 0, 3), Distribution("U", 0, 3))
filter5 = Convolution(Distribution("U", 0, 4), Distribution("U", 0, 4), Distribution("U", 0, 4))
filter6 = Convolution(Distribution("U", 0, 6), Distribution("U", 0, 6))
Plot(filter1, filter2, Distribution("B", 0, 12, 0.5), filter3, filter4, filter5, filter6, Distribution("U", 0, 12))

seq76 = MovingAverage(seq70, Distribution("B", 0, 16, 0.5), BeginEnd=True)

seq77 = MovingAverage(seq70, Distribution("B", 0, 16, 0.5), BeginEnd=True, Output="Residual")
acf12 = Merge(ComputeCorrelation(seq73, 1, MaxLag=10), ComputeCorrelation(seq73, 2, MaxLag=10))
ComputeWhiteNoiseCorrelation(acf12, Distribution("B", 0, 6, 0.5))
Plot(acf12)

seq78 = MovingAverage(seq70, [1, 1, 1], BeginEnd=True)
seq79 = MovingAverage(seq70, [1, 1, 1], BeginEnd=True, Output="Residual")
acf13 = Merge(ComputeCorrelation(seq75, 1, MaxLag=10), ComputeCorrelation(seq75, 2, MaxLag=10))
ComputeWhiteNoiseCorrelation(acf13, [1, 1, 1])
Plot(acf13)

# multiple change-point models

# "Gaussian", "Multinomial", "Poisson", "Ordinal", "Variance",
# "Mean" (multivariate specification), "MeanVariance" (multivariate specification)

seq80 = SelectVariable(seq70, 1)
Display(seq80, 2, 5, "Gaussian", ViewPoint="SegmentProfile", NbSegmentation=5)
Plot(seq80, 2, 5, "Gaussian", ViewPoint="SegmentProfile")
Plot(seq80, 2, 5, "Gaussian", ViewPoint="SegmentProfile", Output="ChangePoint")

seq81 = Segmentation(seq80, 2, 10, "Gaussian")
Plot(seq81, ViewPoint="Data")

seq82 = Segmentation(seq80, [5, 6, 6, 4, 4, 4], "Gaussian")
seq83 = SelectVariable(seq82, 1, Mode="Reject")
Plot(SelectIndividual(seq83, [1, 2]), ViewPoint="Data")
Plot(SelectIndividual(seq83, [4, 5]), ViewPoint="Data")

# multivariate segmentation

Display(seq70, 5, 4, "Gaussian", "Gaussian", ViewPoint="SegmentProfile", NbSegmentation=5)
Plot(seq70, 5, 4, "Gaussian", "Gaussian", ViewPoint="SegmentProfile")
Plot(seq70, 5, 4, "Gaussian", "Gaussian", ViewPoint="SegmentProfile", Output="ChangePoint")

# estimation of a hidden semi-Markov chain

hmc60 = HiddenSemiMarkov(get_shared_data( "pin_laricio_6.hsc"))
hmc6 = Estimate(seq70, "HIDDEN_SEMI-MARKOV", hmc60)

hsmc60 = HiddenSemiMarkov(get_shared_data( "pin_laricio_6.hsc"))
hsmc6 = Estimate(seq70, "HIDDEN_SEMI-MARKOV", hsmc60)
hsmc61 = Estimate(seq70, "HIDDEN_SEMI-MARKOV", "Ordinary", 6, "LeftRight")

Plot(ExtractDistribution(hsmc6, "Observation", 1, 0), ExtractDistribution(hsmc6, "Observation", 1, 1), ExtractDistribution(hsmc6, "Observation", 1, 2), ExtractDistribution(hsmc6, "Observation", 1, 3), ExtractDistribution(hsmc6, "Observation", 1, 4), ExtractDistribution(hsmc6, "Observation", 1, 5))
Plot(ExtractDistribution(hsmc6, "Observation", 2, 0), ExtractDistribution(hsmc6, "Observation", 2, 1), ExtractDistribution(hsmc6, "Observation", 2, 2), ExtractDistribution(hsmc6, "Observation", 2, 3), ExtractDistribution(hsmc6, "Observation", 2, 4), ExtractDistribution(hsmc6, "Observation", 2, 5))

# 1, 3, 5
Plot(hsmc6, 5, ViewPoint="StateProfile")

seq60 = ExtractData(hsmc6)
Display(seq60, ViewPoint="Data", Format="Line")

seq61 = SojournTimeSequences(seq60, 1)
Display(seq61, ViewPoint="Data", Format="Line")





mixt61 = Mixture(21. / 406., ExtractDistribution(hsmc6, "Observation", 1, 0), 29. / 406., ExtractDistribution(hsmc6, "Observation", 1, 1), 140. / 406., ExtractDistribution(hsmc6, "Observation", 1, 2), 87. / 406., ExtractDistribution(hsmc6, "Observation", 1, 3), 71. / 406., ExtractDistribution(hsmc6, "Observation", 1, 4), 58. / 406., ExtractDistribution(hsmc6, "Observation", 1, 5))
mixt61 = Mixture(0.0497296, ExtractDistribution(hsmc6, "Observation", 1, 0), 0.0750034, ExtractDistribution(hsmc6, "Observation", 1, 1), 0.3416, ExtractDistribution(hsmc6, "Observation", 1, 2), 0.207806, ExtractDistribution(hsmc6, "Observation", 1, 3), 0.159426, ExtractDistribution(hsmc6, "Observation", 1, 4), 0.166434, ExtractDistribution(hsmc6, "Observation", 1, 5))
Plot(Fit(ExtractHistogram(seq70, "Value", 1), ExtractDistribution(mixt61, "Mixture")))

mixt62 = Mixture(21. / 406., ExtractDistribution(hsmc6, "Observation", 2, 0), 29. / 406., ExtractDistribution(hsmc6, "Observation", 2, 1), 140. / 406., ExtractDistribution(hsmc6, "Observation", 2, 2), 87. / 406., ExtractDistribution(hsmc6, "Observation", 2, 3), 71. / 406., ExtractDistribution(hsmc6, "Observation", 2, 4), 58. / 406., ExtractDistribution(hsmc6, "Observation", 2, 5))
mixt62 = Mixture(0.0497296, ExtractDistribution(hsmc6, "Observation", 2, 0), 0.0750034, ExtractDistribution(hsmc6, "Observation", 2, 1), 0.3416, ExtractDistribution(hsmc6, "Observation", 2, 2), 0.207806, ExtractDistribution(hsmc6, "Observation", 2, 3), 0.159426, ExtractDistribution(hsmc6, "Observation", 2, 4), 0.166434, ExtractDistribution(hsmc6, "Observation", 2, 5))
Plot(Fit(ExtractHistogram(seq70, "Value", 2), ExtractDistribution(mixt62, "Mixture")))


# comparason with the segmentations deduced from the 6-state hidden semi-Markov chain

seq46 = Merge( Segmentation(seq70, 1, [1935, 1961, 1972, 1990], "Gaussian", "Gaussian"),
Segmentation(seq70, 2, [1932, 1936, 1961, 1984, 1986], "Gaussian", "Gaussian"),
Segmentation(seq70, 3, [1932, 1949, 1971, 1985, 1990], "Gaussian", "Gaussian"),
Segmentation(seq70, 4, [1930, 1953, 1963, 1977], "Gaussian", "Gaussian"),
Segmentation(seq70, 5, [1931, 1963, 1975, 1991], "Gaussian", "Gaussian"),
Segmentation(seq70, 6, [1931, 1943, 1960, 1976], "Gaussian", "Gaussian"))

seq47 = SelectVariable(seq46, [3, 5])
Plot(seq47, ViewPoint="Data")

seq48 = Segmentation(seq80, [5, 5, 6, 5, 5, 4], "Gaussian")
seq49 = SelectVariable(seq48, [3])
# these two lines works together
# seq47 = SelectVariable(seq46, [3])
#Plot(Merge(SelectIndividual(seq47, [1]), SelectIndividual(seq49, [1])), ViewPoint=Data)


# analyse des residus

seq50 = Segmentation(seq80, [5, 6, 6, 4, 4, 4], "Gaussian", Output="Residual")
seq51 = PointwiseAverage(seq50, StandardDeviation=True)
Plot(SelectIndividual(seq51, [0, 7]), ViewPoint="Data")
Plot(Regression(Vectors(seq50), "MovingAverage", 1, 3, [1]))

acf50 = ComputeCorrelation(seq50, 1, MaxLag=10)
Plot(acf50)
seq51 = Merge(SegmentationExtract(seq50, 1, 0), SegmentationExtract(seq50, 1, 1), SegmentationExtract(seq50, 1, 2), SegmentationExtract(seq50, 1, 3), SegmentationExtract(seq50, 1, 4), SegmentationExtract(seq50, 1, 5))
acf51 = ComputeCorrelation(seq51, MaxLag=10)
Plot(acf51)

seq52 = Segmentation(seq80, 2, 6, "Gaussian", NbChangePoint="Fixed", Output="Residual")


seq55 = Segmentation(seq80, [6, 5, 5, 6, 4, 4], "Mean", Output="Residual")
seq56 = PointwiseAverage(seq55, StandardDeviation=True)
Plot(SelectIndividual(seq56, [0, 7]), ViewPoint="Data")
Plot(Regression(Vectors(seq55), "MovingAverage", 1, 3, [1]))

seq57 = Segmentation(seq80, [5, 5, 5, 4, 4, 4], "Mean")
seq58 = Segmentation(seq80, [5, 5, 5, 4, 4, 4], "Gaussian")
Display(MergeVariable(SelectVariable(seq57, 1), seq58), ViewPoint="Data", Format="Line")

