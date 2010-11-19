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
from openalea.sequence_analysis import get_shared_data as path
from os.path import join as pj
seq0 = Sequences(pj(path ,"chene_sessile_15pa.seq"))
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


#print type(ExtractHistogram(vec95, 2))

Compare(ExtractHistogram(vec95, 2), ExtractHistogram(vec96, 2), ExtractHistogram(vec97, 2), "N")
Plot(ExtractHistogram(vec95, 2), ExtractHistogram(vec96, 2), ExtractHistogram(vec97, 2))



