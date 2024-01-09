"""
#########################################################################
#
#  Sequence samples
#
#########################################################################
#
#  Objective: identifying repeated patterns in branching sequences
#             whose index parameter is the node rank,
#
#  Methods: extraction of characteristic distributions
#           (intensity, recurrence time distributions),
#           computation of sample autocorrelation functions,
#           estimation of the order of a Markov chain.
#
#  Data: Daniel Barthelemy
#
#  Cypress dupreziana: order 2 axes in selected positions (20, 40, 60 et 80 nodes from the apex)
#
#  VARIABLE 1 : number of offspring shoots per node.
#
#########################################################################
"""
__revision__ = "$Id: test_exploratory2.py 9885 2010-11-06 18:19:34Z cokelaer $"

from openalea.sequence_analysis import *
from openalea.sequence_analysis.estimate import Estimate as Estimate

from openalea.sequence_analysis import get_shared_data

def test_exploratory():
    seq19 = Sequences(get_shared_data( "dupreziana_20a2.seq"))
    seq20 = RemoveRun(seq19, 0, "End")
    histo201 = ExtractHistogram(seq20, "Recurrence", 1)
    histo202 = ExtractHistogram(seq20, "Recurrence", 2)

    seq38 = Sequences(get_shared_data( "dupreziana_40a2.seq"))
    seq39 = RemoveRun(seq38, 2, 0, "End")
    seq40 = SegmentationExtract(seq39, 1, 2)
    histo401 = ExtractHistogram(seq40, "Recurrence", 1)
    histo402 = ExtractHistogram(seq40, "Recurrence", 2)

    seq58 = Sequences(get_shared_data( "dupreziana_60a2.seq"))
    seq59 = RemoveRun(seq58, 2, 0, "End")
    seq60 =  LengthSelect(SegmentationExtract(seq59, 1, 2), 1, Mode="Reject")
    histo601 = ExtractHistogram(seq60, "Recurrence", 1)
    histo602 = ExtractHistogram(seq60, "Recurrence", 2)

    seq78 = Sequences(get_shared_data( "dupreziana_80a2.seq"))
    seq79 = RemoveRun(seq78, 2, 0, "End")
    seq80 = SegmentationExtract(seq79, 1, 2)
    histo801 = ExtractHistogram(seq80, "Recurrence", 1)
    histo802 = ExtractHistogram(seq80, "Recurrence", 2)


    Plot(histo201, histo401, histo601, histo801)
    Plot(histo202, histo402, histo602, histo802)
    Plot(ExtractHistogram(seq20, "Length"), ExtractHistogram(seq40, "Length"), ExtractHistogram(seq60, "Length"), ExtractHistogram(seq80, "Length"))
    
    seq10 = Merge(seq20, seq40, seq60, seq80)
    Display(seq10, ViewPoint="Data")
    #Plot(seq10, "Intensity")
    #Plot(seq10, "Recurrence")
    #Plot(seq10, "Sojourn")

    # plot of a sample Spearman (rank based) autocorrelation function

    Plot(ComputeCorrelation(seq10, Type="Spearman", MaxLag=15, Normalization="Exact"))

    seq11 = Transcode(seq10, [0, 1, 0])
    seq12 = Transcode(seq10, [0, 0, 1])
    acf1 = Merge(ComputeCorrelation(seq11, MaxLag=15, Normalization="Exact"),\
             ComputeCorrelation(seq12, MaxLag=15, Normalization="Exact"))
    Plot(acf1)
    Display(acf1)

    acf2 = Merge(ComputeCorrelation(seq11, Type="Spearman", MaxLag=15, Normalization="Exact"),\
             ComputeCorrelation(seq12, Type="Spearman", MaxLag=15, Normalization="Exact"))
    acf3 = Merge(ComputeCorrelation(seq11, Type="Kendall", MaxLag=15),\
             ComputeCorrelation(seq12, Type="Kendall", MaxLag=15))

# model selection approach: estimation of both the parameters (initial probabilities and
# transition probabilities) and the order (memory length) of a Markov chain

#todo
#mc10 = Estimate(seq10, "MARKOV", MaxOrder=4)
#Plot(mc10, "Intensity")
#Plot(mc10, "Recurrence")

