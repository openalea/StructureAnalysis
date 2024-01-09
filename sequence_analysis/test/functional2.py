""" functional tests

"""
__revision__ = "$Id: functional2.py 9885 2010-11-06 18:19:34Z cokelaer $"



from openalea.sequence_analysis import *
from openalea.sequence_analysis.estimate import  Estimate
from openalea.sequence_analysis.compare import  Compare
from openalea.sequence_analysis import get_shared_data


seq20 = Sequences(get_shared_data("belren1.seq"))
seq21 = Sequences(get_shared_data("elstar1.seq"))
seq22 = Sequences(get_shared_data("fuji1.seq"))
seq23 = Sequences(get_shared_data("gala1.seq"))
seq24 = Sequences(get_shared_data("granny1.seq"))
seq25 = Sequences(get_shared_data("reinet1.seq"))
seq26 = Sequences(get_shared_data("wij1.seq"))

Display(seq25, ViewPoint="Data")
Plot(seq25, "Intensity")
Plot(seq25, "Sojourn")

seq26 = Reverse(seq25)
Plot(seq26, "Intensity")
Plot(seq26, "FirstOccurrence")

# Sojourn time (run length) distributions

seq30 = Merge(seq20, seq21, seq22, seq23, seq24, seq25)
Plot(seq30, "Sojourn")
Plot(ExtractHistogram(seq30, "Sojourn", 1), ExtractHistogram(seq30, "Sojourn", 2), ExtractHistogram(seq30, "Sojourn", 3), ExtractHistogram(seq30, "Sojourn", 4))

mc30 = Estimate(seq30, "VARIABLE_ORDER_MARKOV", "Ordinary", MaxOrder=4, GlobalInitialTransition=False)
mc30 = Estimate(seq30, "VARIABLE_ORDER_MARKOV", "Ordinary", MaxOrder=4, Algorithm="BIC", GlobalInitialTransition=False)
#todo empty plot
#Plot(mc30, "Sojourn")
Display(Estimate(seq30, "VARIABLE_ORDER_MARKOV", "Ordinary", Order=1))
Display(Estimate(seq30, "VARIABLE_ORDER_MARKOV", "Ordinary", Order=2, GlobalInitialTransition=False))

seq31 = Cluster(seq30, "Limit", [1, 4])
mc31 = Estimate(seq30, "VARIABLE_ORDER_MARKOV", "Ordinary", MaxOrder=4, GlobalInitialTransition=False)
mc31 = Estimate(seq31, "VARIABLE_ORDER_MARKOV", "Ordinary", Order=2, GlobalInitialTransition=False)
Plot(mc31, "Sojourn")
Display(Estimate(seq31, "VARIABLE_ORDER_MARKOV", "Ordinary", Order=1))

# comparison of sequences by dynamic programming algorithms

seq32 = Merge(seq20, seq25)
matrix30 = Compare(seq32)
matrix31 = Compare(seq32, VectorDistance("S"))
matrix32 = Compare(seq32, VectorDistance("S"), Transposition=True)
matrix33 = Compare(seq32, VectorDistance(get_shared_data("test_align1.a")), Transposition=True)

Display(Clustering(matrix33, "Partition", 2))
Clustering(matrix33, "Hierarchy", Algorithm="Agglomerative")
Clustering(matrix33, "Hierarchy", Algorithm="Divisive")

# multiple alignment

seq33 = Compare(SelectIndividual(seq25, [10, 11, 12, 14, 15]), VectorDistance("S"), Output="Sequences", Algorithm="Agglomerative")
seq34 = Compare(SelectIndividual(seq25, [10, 11, 12, 14, 15]), VectorDistance("S"), Output="Sequences", Algorithm="Divisive")
seq35 = Compare(SelectIndividual(seq25, [10, 11, 12, 14, 15]), VectorDistance("S"), Output="Sequences", Algorithm="Ordering")

Compare(seq25, TestSequence=9, RefSequence=1)
Compare(seq25, VectorDistance("S"), TestSequence=9, RefSequence=1)
Compare(seq25, VectorDistance("S"), TestSequence=9, RefSequence=1, Transposition=True)

# multiple change-point models

Display(seq25, 14, 6, "Multinomial", ViewPoint="SegmentProfile")
Display(seq25, 14, 6, "Multinomial", ViewPoint="SegmentProfile", Output="ChangePoint")
Plot(seq25, 14, 6, "Multinomial", ViewPoint="SegmentProfile")
Plot(seq25, 14, 6, "Multinomial", ViewPoint="SegmentProfile", Output="ChangePoint")
# hidden semi-Markov chains

hsmc0 = HiddenSemiMarkov(get_shared_data("belren1.hsc"))
hsmc20 = Estimate(seq20, "HIDDEN_SEMI-MARKOV", hsmc0)

hsmc0 = HiddenSemiMarkov(get_shared_data("elstar1.hsc"))
hsmc21 = Estimate(seq21, "HIDDEN_SEMI-MARKOV", hsmc0)

hsmc0 = HiddenSemiMarkov(get_shared_data("fuji1.hsc"))
hsmc22 = Estimate(seq22, "HIDDEN_SEMI-MARKOV", hsmc0)

hsmc0 = HiddenSemiMarkov(get_shared_data("gala1.hsc"))
hsmc23 = Estimate(seq23, "HIDDEN_SEMI-MARKOV", hsmc0)

hsmc0 = HiddenSemiMarkov(get_shared_data("granny1.hsc"))
hsmc24 = Estimate(seq24, "HIDDEN_SEMI-MARKOV", hsmc0)

hsmc0 = HiddenSemiMarkov(get_shared_data("reinet1.hsc"))
hsmc25 = Estimate(seq25, "HIDDEN_SEMI-MARKOV", hsmc0)

Display(hsmc25)
Plot(hsmc25, "Intensity", 1)
Plot(hsmc25, "FirstOccurrence", 1)
Plot(hsmc25, "Counting", 1)

# state
Plot(hsmc25, "Intensity")
Plot(hsmc25, "Sojourn")
# observed
Plot(hsmc25, "Sojourn",1)

Plot(hsmc25, 1, ViewPoint="StateProfile")
Plot(hsmc25, 1, ViewPoint="StateProfile", Output='InState')
Plot(hsmc25, 1, ViewPoint="StateProfile", Output='OutState')

seq25_1 = ExtractData(hsmc25)
Display(seq25_1, ViewPoint="Data", Format="Line")

hsmc0 = HiddenSemiMarkov(get_shared_data("wij1.hsc"))
hsmc26 = Estimate(seq26, "HIDDEN_SEMI-MARKOV", hsmc0)

# model comparison

#Thresholding(hsmc20, MinProbability=0.001)
#Thresholding(hsmc21, MinProbability=0.001)
#Thresholding(hsmc22, MinProbability=0.001)
#Thresholding(hsmc23, MinProbability=0.001)
#Thresholding(hsmc24, MinProbability=0.001)
#Thresholding(hsmc25, MinProbability=0.001)
#Thresholding(hsmc26, MinProbability=0.001)


#matrix20 = Compare(Thresholding(hsmc22, MinProbability=0.001), seq22, 10000)

#matrix20 = Compare(Thresholding(hsmc20, MinProbability=0.001), seq20, Thresholding(hsmc21, MinProbability=0.001), seq21, Thresholding(hsmc22, MinProbability=0.001), seq22, Thresholding(hsmc24, MinProbability=0.001), seq24, Thresholding(hsmc25, MinProbability=0.001), seq25, Thresholding(hsmc26, MinProbability=0.001), seq26, 10000)

#TODO unstable the line above works, the line below does not
#matrix20 = Compare(Thresholding(hsmc20, MinProbability=0.001), seq20, Thresholding(hsmc21, MinProbability=0.001), seq21, Thresholding(hsmc22, MinProbability=0.001), seq22, Thresholding(hsmc23, MinProbability=0.001), seq23, Thresholding(hsmc24, MinProbability=0.001), seq24, Thresholding(hsmc25, MinProbability=0.001), seq25, Thresholding(hsmc26, MinProbability=0.001), seq26, 10000, FileName="ASCII/cultivar1_models.txt")

# may be slow 
#matrix21 = Compare(Thresholding(hsmc20, MinProbability=0.001), Thresholding(hsmc21, MinProbability=0.001), Thresholding(hsmc22, MinProbability=0.001), Thresholding(hsmc22, MinProbability=0.001), Thresholding(hsmc24, MinProbability=0.001), Thresholding(hsmc25, MinProbability=0.001), Thresholding(hsmc26, MinProbability=0.001), 100, 90)
#matrix21 = Compare(Thresholding(hsmc20, MinProbability=0.001), Thresholding(hsmc21, MinProbability=0.001), Thresholding(hsmc22, MinProbability=0.001), Thresholding(hsmc22, MinProbability=0.001), Thresholding(hsmc24, MinProbability=0.001), Thresholding(hsmc25, MinProbability=0.001), Thresholding(hsmc26, MinProbability=0.001), 100, 90, FileName="ASCII/cultivar1_models_90.txt")


#Plot(matrix20)
