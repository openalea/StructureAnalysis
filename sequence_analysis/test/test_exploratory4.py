"""
#########################################################################
#
#  Objective: analyzing jointly the phyllotaxis and the number of offspring shoots
#             in sequences whose index parameter is the node rank,
#
#  Methods: extraction of characteristic distributions
#           (intensity, recurrence time and sojourn time distributions),
#           sequence alignment.
#
#  Cypress dupreziana: order 1 axes
#
#  Data: Daniel Barthelemy
#
#  VARIABLE 1 : phyllotaxis (number of leaves per node),
#  VARIABLE 2 : number of offspring shoots per node.
#
#########################################################################
"""
__revision__ = "$Id: test_exploratory4.py 9401 2010-08-10 12:24:59Z cokelaer $"

from openalea.sequence_analysis import *
from openalea.sequence_analysis.compare import Compare as Compare
from openalea.sequence_analysis import get_shared_data


seq1 = Sequences(get_shared_data( "dupreziana_a1.seq"))

Display(seq1, ViewPoint="Data", Format="Line")

vec20 = MergeVariable(ExtractVectors(seq1, "NbOccurrence", 1, 3), ExtractVectors(seq1, "Length"))
Display(vec20)
#todo 
#Plot(vec20,1)
Plot(vec20)

seq2 = Shift(seq1, 1, -3)

seq3 = SegmentationExtract(seq1, 1, 3)
seq4 = SegmentationExtract(seq1, 1, 4)
Plot(ExtractHistogram(seq3, "Recurrence", 1), ExtractHistogram(seq4, "Recurrence", 1))
Plot(ExtractHistogram(seq3, "Sojourn", 1), ExtractHistogram(seq4, "Sojourn", 1))
Plot(ExtractHistogram(seq3, "Recurrence", 2), ExtractHistogram(seq4, "Recurrence", 2))
Plot(ExtractHistogram(seq3, "Sojourn", 2), ExtractHistogram(seq4, "Sojourn", 2))

Plot(ExtractHistogram(seq4, "Recurrence", 2), ExtractHistogram(seq4, "Recurrence", 3), ExtractHistogram(seq4, "Recurrence", 4))

matrix20 = Compare(seq1, VectorDistance("N", "N"))
Plot(matrix20)
Display(Clustering(matrix20, "Partition", 2))
Clustering(matrix20, "Hierarchy")

matrix21 = Compare(seq1, VectorDistance("O", "O"))
Plot(matrix21)
Display(Clustering(matrix21, "Partition", 2))
Clustering(matrix21, "Hierarchy")

seq11 = SelectIndividual(seq2, [18, 9, 10, 31, 6, 14, 29, 16, 1, 12, 5, 7, 25, 22, 17, 30, 13, 4, 21, 27, 20, 24])
seq12 = SelectIndividual(seq2, [28, 19, 32, 23, 26, 11, 3, 15, 8, 33, 2])
#todo
#Plot(ExtractHistogram(seq11, "FirstOccurrence", 1, 0), ExtractHistogram(seq12, "FirstOccurrence", 1, 0))
ComparisonTest("W", ExtractHistogram(seq11, "Length"), ExtractHistogram(seq12, "Length"))
Plot(ExtractHistogram(seq11, "Length"), ExtractHistogram(seq12, "Length"))

#todo
#Plot(seq2, "Intensity", 1)

#Plot(seq2, "Intensity", 2)
seq5 = RemoveRun(seq2, 2, 0, "End")
Plot(seq5, "Intensity", 1)

