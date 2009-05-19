#########################################################################
#
#  Objective: identifying branching and axillary flowering zones and repeated patterns
#             in sequences whose index parameter is the node rank,
#
#  Methods: extraction of characteristic distributions (intensity,
#           sojourn time distributions), sequence alignment,
#           clustering methods (either partitioning or hierarchical methods).
#
#  1st annual shoot of apple tree trunks
#
#  Data: Evelyne Costes
#
#  VARIABLE 1 : type of axillary production (0: latent bud, 1: one-year-delayed short shoot,
#               2: one-year-delayed long shoot, 3: one-year-delayed flowering shoot,
#               4: immediate shoot).
#
#########################################################################
from openalea.sequence_analysis import *
seq20 = Sequences("data/belren1.seq")
seq21 = Sequences("data/elstar1.seq")
seq22 = Sequences("data/fuji1.seq")
seq23 = Sequences("data/gala1.seq")
seq24 = Sequences("data/granny1.seq")
seq25 = Sequences("data/reinet1.seq")

# Display(seq25, ViewPoint="Data")
# Plot(seq25, "Intensity")
# Plot(seq25, "Sojourn")

seq26 = Reverse(seq25)
# Plot(seq26, "Intensity")
# Plot(seq26, "FirstOccurrence")

# Sojourn time (run length) distributions

seq30 = Merge(seq20, seq21, seq22, seq23, seq24, seq25)
# Plot(seq30, "Sojourn")
# Plot(ExtractHistogram(seq30, "Sojourn", 1), ExtractHistogram(seq30, "Sojourn", 2), ExtractHistogram(seq30, "Sojourn", 3),  ExtractHistogram(seq30, "Sojourn", 4))

#todo
#mc30 = Estimate(seq30, "MARKOV", MaxOrder=4)

# Plot(mc30, "Sojourn")
# Display(Estimate(seq30, "MARKOV"))

seq31 = Cluster(seq30, "Limit", [1, 4])
#todo
#mc31 = Estimate(seq31, "MARKOV", Order=2)

# Plot(mc31, "Sojourn")
# Display(Estimate(seq31, "MARKOV"))

# comparison of sequences by dynamic programming algorithms

seq32 = Merge(seq20, seq25)
matrix30 = Compare(seq32)
matrix31 = Compare(seq32, VectorDistance("S"))
matrix32 = Compare(seq32, VectorDistance("S"), Transposition=True)
matrix33 = Compare(seq32, VectorDistance("data/align1.a"), Transposition=True)

#todo
#Display(Clustering(matrix33, "Partition", 2))
#Clustering(matrix33, "Hierarchy")

Compare(seq25, TestSequence=9, RefSequence=1)
Compare(seq25, VectorDistance("S"), TestSequence=9, RefSequence=1)
Compare(seq25, VectorDistance("S"), TestSequence=9, RefSequence=1, Transposition=True)
