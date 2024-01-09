""" functional tests


.. todo:: to be done
"""
__revision__ = "$Id: functional1.py 9401 2010-08-10 12:24:59Z cokelaer $"

import os
from openalea.stat_tool import *
from openalea.sequence_analysis import *
from openalea.sequence_analysis import get_shared_data

seq1 = Sequences(get_shared_data( 'dupreziana_20a2.seq'))   # correct
seq2 = RemoveRun(seq1, 1, 0, "End")              # correct

histo21 = ExtractHistogram(seq2, "Recurrence", 1)  # correct
histo22 = ExtractHistogram(seq2, "Recurrence", 2)  # correct

seq3 = Sequences(get_shared_data( 'dupreziana_40a2.seq'))   #correct
seq4_0 = RemoveRun(seq3, 2, 0, "End")            #correct
seq4 = SegmentationExtract(seq4_0, 1, 2)         #correct


seq5 = Sequences(get_shared_data( 'dupreziana_60a2.seq'))   #correct
seq6_0 = RemoveRun(seq5, 2, 0, "End")            #correct
seq6 = LengthSelect(SegmentationExtract(seq6_0, 1, 2), 1, Mode="Reject") #correct


seq7 = Sequences(get_shared_data( 'dupreziana_80a2.seq'))   #correct
seq8_0 = RemoveRun(seq7, 2, 0, "End")            #correct
seq8 = SegmentationExtract(seq8_0, 1, 2)         #correct


seq10 = Merge(seq2, seq4, seq6, seq8)

seq10_1 = RecurrenceTimeSequences(seq10,1)
seq10_2 = RecurrenceTimeSequences(seq10,2)


vec10 = MergeVariable(ExtractVectors(seq10, "Length"),ExtractVectors(seq10, "NbOccurrence",1,1),ExtractVectors(seq10, "NbOccurrence",1,2),ExtractVectors(seq10, "Cumul"))

seq11 = Transcode(seq10, [0, 1, 0])
seq12 = Transcode(seq10, [0, 0, 1])

acf1 = Merge(ComputeCorrelation(seq11, MaxLag=15), ComputeCorrelation(seq12, MaxLag=15))

acf2 = Merge(ComputeCorrelation(seq11, Type="Spearman", MaxLag=15),
             ComputeCorrelation(seq12, Type="Spearman", MaxLag=15))
acf3 = Merge(ComputeCorrelation(seq11, Type="Kendall", MaxLag=15),
             ComputeCorrelation(seq12, Type="Kendall", MaxLag=15))



WordCount(seq10, 3, BeginState=1, EndState=1, MinFrequency=10)
WordCount(seq10, 4, BeginState=2, EndState=2)
WordCount(seq10, 4, BeginState=2, EndState=1)
 


mc10 = Estimate(seq10, "VARIABLE_ORDER_MARKOV", "Ordinary", MaxOrder=5, GlobalInitialTransition=True)

mc11 = Estimate(seq10, "VARIABLE_ORDER_MARKOV", "Ordinary", MaxOrder=5, GlobalInitialTransition=False)

Plot(mc11, "Intensity")

mc12 = Estimate(seq10, "VARIABLE_ORDER_MARKOV", "Ordinary", Algorithm="LocalBIC", Threshold=10., MaxOrder=5, GlobalInitialTransition=False, GlobalSample=False)
mc13 = Estimate(seq10, "VARIABLE_ORDER_MARKOV", "Ordinary", Algorithm="Context", Threshold=1., MaxOrder=5, GlobalInitialTransition=False, GlobalSample=False)

acf11 = ComputeAutoCorrelation(mc11, 1, MaxLag=20)
acf12 = ComputeAutoCorrelation(mc11, 2, MaxLag=20)

mc2 = Estimate(seq2, "VARIABLE_ORDER_MARKOV", mc11, GlobalInitialTransition=False)
mc4 = Estimate(seq4, "VARIABLE_ORDER_MARKOV", mc11, GlobalInitialTransition=False)
mc6 = Estimate(seq6, "VARIABLE_ORDER_MARKOV", mc11, GlobalInitialTransition=False)
mc8 = Estimate(seq8, "VARIABLE_ORDER_MARKOV", mc11, GlobalInitialTransition=False)


#TODO compare functions crashes sometimes
#matrix1 = Compare(Thresholding(mc2, MinProbability=0.001), seq10, Thresholding(mc4, MinProbability=0.001), seq10, Thresholding(mc6, MinProbability=0.001), seq10, Thresholding(mc8, MinProbability=0.001), seq10, 10000)
#matrix2 = Compare(Thresholding(mc2, MinProbability=0.001), seq2, Thresholding(mc4, MinProbability=0.001), seq4,	Thresholding(mc6, MinProbability=0.001), seq6,	Thresholding(mc8, MinProbability=0.001), seq8, 10000)


#Compare(seq10, Thresholding(mc2, MinProbability=0.001), Thresholding(mc4, MinProbability=.001), Thresholding(mc6, MinProbability=0.001), Thresholding(mc8, MinProbability=0.001))



# test #
hmc9 = HiddenVariableOrderMarkov(get_shared_data( "dupreziana21.hc"))
hmc10 = Estimate(seq10, "HIDDEN_VARIABLE_ORDER_MARKOV", hmc9, GlobalInitialTransition=True, NbIteration=80)
hmc11 = Estimate(seq10, "HIDDEN_VARIABLE_ORDER_MARKOV", hmc9, GlobalInitialTransition=False, NbIteration=80)


acf21 = ComputeAutoCorrelation(hmc11, 1, 1, MaxLag=20)
acf22 = ComputeAutoCorrelation(hmc11, 1, 2, MaxLag=20)

seq15 = Simulate(hmc11, 10000, seq10)

