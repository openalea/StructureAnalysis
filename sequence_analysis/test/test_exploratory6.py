"""
#########################################################################
#
#  Objective: analyzing the apical growth and the branching of Corsican pines,
#
#  Methods: extraction of trends and residuals by linear filtering
#           (differencing and weighted moving average), computation of
#           sample autocorrelation functions.
#
#  75-year-old Corsican pines (66 years measured).
#
#  Data: Celine Meredieu and Yves Caraglio
#
#  VARIABLE 1 : year of growth (explicit index parameter),
#  VARIABLE 2 : length of the annual shoot (mm),
#  VARIABLE 3 : number of  branches per annual shoot / tier.
#
#########################################################################
"""
__revision__ = "$Id: test_exploratory6.py 9401 2010-08-10 12:24:59Z cokelaer $"


from openalea.sequence_analysis import *
from openalea.sequence_analysis.compare import Compare as Compare
from openalea.sequence_analysis import get_shared_data

seq66 = Sequences(get_shared_data("laricio_date66.seq"))
Plot(seq66, ViewPoint="Data")
#Plot(Cumulate(seq66), ViewPoint="Data")

vec66 = Vectors(seq66)
regress66_1 = Regression(vec66, "MovingAverage", 1, 2, [1])
Plot(regress66_1)
regress66_2 = Regression(vec66, "MovingAverage", 1, 3, [1])

regress66_23 = Regression(vec66, "NearestNeighbours", 2, 3, 0.3)
Display(regress66_23)
Plot(regress66_23)

vec70 =  Vectors(SelectIndividual(seq66, [1, 2, 3]))
regress70_1 = Regression(vec70, "MovingAverage", 1, 2, [1])
Plot(regress70_1)
regress70_2 = Regression(vec70, "MovingAverage", 1, 3, [1])

vec71 =  Vectors(SelectIndividual(seq66, [4, 5, 6]))
regress71_1 = Regression(vec71, "MovingAverage", 1, 2, [1])
Plot(regress71_1)
regress71_2 = Regression(vec71, "MovingAverage", 1, 3, [1])

matrix66 = Compare(SelectVariable(seq66, 1, Mode="Reject"), VectorDistance("N", "N"))
Display(Clustering(matrix66, "Partition", 3))
Clustering(matrix66, "Hierarchy")

# extraction of trends (slowly varying component) and residuals (rapidly varying component)
# by symmetric smoothing filters and computation of sample autocorrelation functions from residuals

seq67 = Difference(seq66)
acf11 = Merge(ComputeCorrelation(seq67, 2, MaxLag=10),\
               ComputeCorrelation(seq67, 3, MaxLag=10))
acf11 = Merge(ComputeCorrelation(seq67, 2, MaxLag=10, Normalization="Exact"),\
              ComputeCorrelation(seq67, 3, MaxLag=10, Normalization="Exact"))
ComputeWhiteNoiseCorrelation(acf11, 1)
Plot(acf11)

# symmetric smoothing filters of half-width 3

filter1 = Convolution(Distribution("B", 0, 3, 0.2), Distribution("B", 0, 3, 0.8))
filter2 = Convolution(Distribution("B", 0, 2, 0.2), Distribution("B", 0, 2, 0.5), Distribution("B", 0, 2, 0.8))
filter3 = Convolution(Distribution("U", 0, 2), Distribution("U", 0, 2), Distribution("U", 0, 2))
filter4 = Convolution(Distribution("U", 0, 3), Distribution("U", 0, 3))
Plot(filter1, filter2, Distribution("B", 0, 6, 0.5), filter3, filter4, Distribution("U", 0, 6))

seq68 = MovingAverage(seq66, Distribution("B", 0, 6, 0.5), BeginEnd=True)

seq69 = MovingAverage(VariableScaling(seq66, 3, 100), Distribution("B", 0, 6, 0.5), BeginEnd=True, Output="Residual")
acf12 = Merge(ComputeCorrelation(seq69, 2, MaxLag=10),\
               ComputeCorrelation(seq69, 3, MaxLag=10))
acf12 = Merge(ComputeCorrelation(seq69, 2, MaxLag=10, Normalization="Exact"),\
              ComputeCorrelation(seq69, 3, MaxLag=10, Normalization="Exact"))
ComputeWhiteNoiseCorrelation(acf12, Distribution("B", 0, 6, 0.5))
Plot(acf12)

seq70 = MovingAverage(seq66, [1, 1, 1], BeginEnd=True)
seq71 = MovingAverage(VariableScaling(seq66, 3, 100), [1, 1, 1], BeginEnd=True, Output="Residual")
acf13 = Merge(ComputeCorrelation(seq71, 2, MaxLag=10),\
               ComputeCorrelation(seq71, 3, MaxLag=10))
acf13 = Merge(ComputeCorrelation(seq71, 2, MaxLag=10, Normalization="Exact"),\
              ComputeCorrelation(seq71, 3, MaxLag=10, Normalization="Exact"))
ComputeWhiteNoiseCorrelation(acf13, [1, 1, 1])
Plot(acf13)

seq80 = Sequences(get_shared_data( "laricio_position66.seq"), OldFormat=True)

#Plot(Cumulate(seq80), ViewPoint="Data")
