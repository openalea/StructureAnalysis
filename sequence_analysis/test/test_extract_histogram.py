""" Test extract histogram

.. author:: Thomas Cokelaer, Thomas.Cokelaer@inria.fr

.. todo:: to be done
"""
__revision__ = "$Id: test_extract_histogram.py 9885 2010-11-06 18:19:34Z cokelaer $"


from openalea.stat_tool.data_transform import ValueSelect, ExtractHistogram
from openalea.sequence_analysis import *

seq = Sequences(get_shared_data("pin_laricio_7x.seq"))
seq_cluster = Cluster(seq, "Step", 1, 10)
_seq1 = Sequences(get_shared_data('dupreziana_20a2.seq'))
seq2 = RemoveRun(_seq1, 1, 0, "End")
seq3 = Sequences(get_shared_data('dupreziana_40a2.seq'))
seq4_0 = RemoveRun(seq3, 2, 0, "End")
seq4 = SegmentationExtract(seq4_0, 1, 2)
seq5 = Sequences(get_shared_data('dupreziana_60a2.seq'))
seq6_0 = RemoveRun(seq5, 2, 0, "End")
seq6 = LengthSelect(SegmentationExtract(seq6_0, 1, 2), 1, Mode="Reject")
seq7 = Sequences(get_shared_data('dupreziana_80a2.seq'))
seq8_0 = RemoveRun(seq7, 2, 0, "End")
seq8 = SegmentationExtract(seq8_0, 1, 2)
seq10 = Merge(seq2, seq4, seq6, seq8)
seq11 = Transcode(seq10, [0, 1, 0])

seq0 = Sequences(get_shared_data("chene_sessile_15pa.seq"))
vec10 = Vectors(seq0)
vec95 = ValueSelect(vec10, 1, 95)
seq20 = Sequences(get_shared_data("belren1.seq"))

def test_vectors():
    ExtractHistogram(vec95, 2)

def test_sequences():
    ExtractHistogram(seq20, "Recurrence", 1)
    ExtractHistogram(seq20, "Recurrence", 2)
    ExtractHistogram(seq20, "Length")


    ExtractHistogram(seq11, "FirstOccurrence", 1, 0)


    ExtractHistogram(seq0, "Value", 3)

    ExtractHistogram(seq_cluster, "Value", 1)
    ExtractHistogram(seq_cluster, "Value", 2)

def test_time_events():
    """not implemented"""

def test_renewal():
    """not implemented"""

