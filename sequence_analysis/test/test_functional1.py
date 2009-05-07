import os
from openalea.sequence_analysis.data_transform import *
from openalea.sequence_analysis.sequences import *
from openalea.stat_tool import *


path = 'test' + os.sep + 'data' +os.sep
path = 'data' + os.sep

seq1 = Sequences(path + 'dupreziana_20a2.seq')   # correct
seq2 = RemoveRun(seq1, 1, 0, "End")              # correct

histo21 = ExtractHistogram(seq2, "Recurrence", 1)  # correct
histo22 = ExtractHistogram(seq2, "Recurrence", 2)  # correct

seq3 = Sequences(path + 'dupreziana_40a2.seq')   #correct
seq4_0 = RemoveRun(seq3, 2, 0, "End")            #correct
seq4 = SegmentationExtract(seq4_0, 1, 2)         #correct


seq5 = Sequences(path + 'dupreziana_60a2.seq')   #correct
seq6_0 = RemoveRun(seq5, 2, 0, "End")            #correct
seq6 = LengthSelect(SegmentationExtract(seq6_0, 1, 2), 1, Mode="Reject") #correct


seq7 = Sequences(path + 'dupreziana_80a2.seq')   #correct
seq8_0 = RemoveRun(seq7, 2, 0, "End")            #correct
seq8 = SegmentationExtract(seq8_0, 1, 2)         #correct


seq10 = Merge(seq2, seq4, seq6, seq8)

seq10_1 = RecurrenceTimeSequences(seq10,1)
seq10_2 = RecurrenceTimeSequences(seq10,2)


vec10 = MergeVariable(ExtractVectors(seq10, "Length"),ExtractVectors(seq10, "NbOccurrence",1,1),ExtractVectors(seq10, "NbOccurrence",1,2),ExtractVectors(seq10, "Cumul"))


#to be done
#seq11 = Transcore(seq10, [0,1,0])
#seq12 = Transcore(seq10, [0,0,1])
#ComputeCorrlection(seq10, type=Spearman, maxlag=15) Perason by default, spearman or kendall as well

#WordCount ?
#Estimate(seq10, "Variable_ORDER_MARKOV","ordinary", maxorder=5, globalInitialtransition=False)

#ComputeAutoCorrelation# 

#Comapre
#Thresholding
#hiddentVaraibleOrderMarkov Estimate, and so on
#simulate, fit
