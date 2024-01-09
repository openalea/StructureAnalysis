"""
#########################################################################
#
#  Oak trunk annual shoots
#
#  Data: Patrick Heuret
#
#  VARIABLE 1 : year of growth (95, 96, 97) (index parameter of sequences)
#  VARIABLE 2 : length of the annual shoot (cm)
#  VARIABLE 3 : diameter of the annual shoot (1/10 de mm)
#  VARIABLE 4 : number of cycles
#  VARIABLE 5 : number of nodes
#  VARIABLE 6 : number de branches
#
#########################################################################
"""
__revision__ = "$Id: test_exploratory.py 9401 2010-08-10 12:24:59Z cokelaer $"

from openalea.sequence_analysis import *
from openalea.sequence_analysis.estimate import  Estimate
from openalea.sequence_analysis import get_shared_data

seq0 = Sequences(get_shared_data("chene_sessile_15pa.seq"))
Plot(seq0, ViewPoint="Data")

# change of unit for the variable diameter of the annual shoot

marginal3 = ExtractHistogram(seq0, "Value", 3)
Plot(Cluster(marginal3, "Information", 0.75))
Plot(Cluster(marginal3, "Information", 0.61))
Plot(Cluster(marginal3, "Step", 10))

vec10 = Vectors(seq0)

# plot of the average sequence
Plot(Regression(vec10, "MovingAverage", 1, 2, [1]))

vec95 = ValueSelect(vec10, 1, 95)
vec96 = ValueSelect(vec10, 1, 96)
vec97 = ValueSelect(vec10, 1, 97)

VarianceAnalysis(vec10, 1, 2, "N")


print type(ExtractHistogram(vec95, 2))

Compare(ExtractHistogram(vec95, 2), ExtractHistogram(vec96, 2), ExtractHistogram(vec97, 2), "N")
Plot(ExtractHistogram(vec95, 2), ExtractHistogram(vec96, 2), ExtractHistogram(vec97, 2))

ContingencyTable(vec10, 1, 4)

# one-way variance analysis based on ranks

VarianceAnalysis(vec10, 1, 4, "O")
Compare(ExtractHistogram(vec95, 4), ExtractHistogram(vec96, 4), ExtractHistogram(vec97, 4), "O")
Plot(ExtractHistogram(vec95, 4), ExtractHistogram(vec96, 4), ExtractHistogram(vec97, 4))

Plot(ExtractHistogram(vec95, 5), ExtractHistogram(vec96, 5), ExtractHistogram(vec97, 5))
Plot(ExtractHistogram(vec95, 6), ExtractHistogram(vec96, 6), ExtractHistogram(vec97, 6))

vec11 = ValueSelect(vec10, 4, 1)
vec12 = ValueSelect(vec10, 4, 2)
vec13 = ValueSelect(vec10, 4, 3, 4)

Plot(ExtractHistogram(vec11, 2), ExtractHistogram(vec12, 2), ExtractHistogram(vec13, 2))
Plot(ExtractHistogram(vec11, 5), ExtractHistogram(vec12, 5), ExtractHistogram(vec13, 5))

mixt20 = Estimate(ExtractHistogram(vec10, 2), "MIXTURE", "NB", "NB", "NB", "NB", NbComponent="Estimated")
Display(mixt20)
Plot(mixt20)
Plot(ExtractDistribution(mixt20, "Mixture"))

mixt21 = Estimate(ExtractHistogram(vec10, 5), "MIXTURE", "NB", "NB", "NB", "NB", NbComponent="Estimated")

vec9596 = ValueSelect(vec10, 1, 95, 96)
Plot(ExtractHistogram(ValueSelect(vec9596, 4, 1), 6), ExtractHistogram(ValueSelect(vec9596, 4, 2), 6), ExtractHistogram(ValueSelect(vec9596, 4, 3, 4), 6))

regress10 = Regression(vec10, "Linear", 5, 2)
Display(regress10)
Plot(regress10)

# nonparametric regression (loess smoother)

regress11 = Regression(vec10, "NearestNeighbours",  5, 2, 0.3)

regress12 = Regression(vec9596, "Linear", 5, 6)
regress13 = Regression(vec9596, "NearestNeighbours", 5, 6, 0.5)

vec15 = SelectVariable(vec10, [1, 3, 6], Mode="Reject")

# computation of a distance matrix using a standardization procedure

matrix10 = Compare(vec15, VectorDistance("N", "N", "N"))

# clustering using a partitioning method

Display(Clustering(matrix10, "Partition", 2))

vec151 = SelectIndividual(vec10,  [69, 48, 41, 44, 32, 47, 81, 95, 11, 36, 75, 108, 56, 83, 38, 98, 113, 134, 110, 101, 77, 35, 74, 80, 50, 24, 89, 128, 5, 45, 8, 116, 119, 132, 61, 78, 53, 29, 131, 65, 90, 96, 104, 20, 86, 66, 42, 68, 125, 14, 23, 54, 33, 26, 71, 129, 102, 51, 70, 111, 138, 19, 127, 62, 117, 137, 2, 28, 17])
vec152 = SelectIndividual(vec10, [100, 13, 133, 105, 72, 9, 93, 109, 30, 115, 63, 7, 55, 37, 15, 114, 106, 46, 73, 18, 3, 87, 58, 43, 60, 76, 52, 6, 39, 31, 12, 99, 121, 123, 22, 79, 94, 88, 21, 97, 25, 40, 57, 136, 67, 49, 10, 4, 120, 92, 27, 91, 64, 124, 16, 130, 84, 107, 126, 103, 122, 112, 59, 1, 82, 34, 135, 118, 85])
Plot(ExtractHistogram(vec151, 4), ExtractHistogram(vec152, 4))

matrix11 = Compare(vec15, VectorDistance("N", "O", "N"))

vec16 = SelectVariable(vec9596, [1, 3], Mode="Reject")
matrix12 = Compare(vec16, VectorDistance("N", "N", "N", "N"))
matrix13 = Compare(vec16, VectorDistance("N", "O", "N", "N"))


