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
__revision__ = "$Id: test_exploratory.py 8676 2010-04-20 15:28:46Z cokelaer $"

from openalea.sequence_analysis import *
from openalea.sequence_analysis.estimate import  Estimate
from openalea.sequence_analysis.data import path

seq0 = Sequences(path +"chene_sessile_15pa.seq")
#Plot(seq0, ViewPoint="Data")


# change of unit for the variable diameter of the annual shoot

marginal3 = ExtractHistogram(seq0, "Value", 3)
#Plot(Cluster(marginal3, "Information", 0.75))
#Plot(Cluster(marginal3, "Information", 0.61))
#Plot(Cluster(marginal3, "Step", 10))

vec10 = Vectors(seq0)

# plot of the average sequence
#Plot(Regression(vec10, "MovingAverage", 1, 2, [1]))

vec95 = ValueSelect(vec10, 1, 95)
vec96 = ValueSelect(vec10, 1, 96)
vec97 = ValueSelect(vec10, 1, 97)

VarianceAnalysis(vec10, 1, 2, "N")


print type(ExtractHistogram(vec95, 2))

Compare(ExtractHistogram(vec95, 2), ExtractHistogram(vec96, 2), ExtractHistogram(vec97, 2), "N")
#Plot(ExtractHistogram(vec95, 2), ExtractHistogram(vec96, 2), ExtractHistogram(vec97, 2))

ContingencyTable(vec10, 1, 4)

# one-way variance analysis based on ranks

VarianceAnalysis(vec10, 1, 4, "O")
Compare(ExtractHistogram(vec95, 4), ExtractHistogram(vec96, 4), ExtractHistogram(vec97, 4), "O")
#Plot(ExtractHistogram(vec95, 4), ExtractHistogram(vec96, 4), ExtractHistogram(vec97, 4))

#Plot(ExtractHistogram(vec95, 5), ExtractHistogram(vec96, 5), ExtractHistogram(vec97, 5))
#Plot(ExtractHistogram(vec95, 6), ExtractHistogram(vec96, 6), ExtractHistogram(vec97, 6))

vec11 = ValueSelect(vec10, 4, 1)
vec12 = ValueSelect(vec10, 4, 2)
vec13 = ValueSelect(vec10, 4, 3, 4)

#Plot(ExtractHistogram(vec11, 2), ExtractHistogram(vec12, 2), ExtractHistogram(vec13, 2))
#Plot(ExtractHistogram(vec11, 5), ExtractHistogram(vec12, 5), ExtractHistogram(vec13, 5))

mixt20 = Estimate(ExtractHistogram(vec10, 2), "MIXTURE", "NB", "NB", "NB", "NB", NbComponent="Estimated")
Display(mixt20)
Plot(mixt20)
Plot(ExtractDistribution(mixt20, "Mixture"))

mixt21 = Estimate(ExtractHistogram(vec10, 5), "MIXTURE", "NB", "NB", "NB", "NB", NbComponent="Estimated")




