"""common data files

.. author:: Thomas Cokelaer, Thomas.Cokelaer@inria.fr
"""

import os
from openalea.stat_tool.vectors import Vectors
from openalea.stat_tool.data_transform import Merge
from openalea.stat_tool import *
from openalea.sequence_analysis import *
from os.path import join as pj

from openalea.sequence_analysis.data import path

vecn = Vectors([[1,2,3],[4,5,6], [7,8,9]])
vec1 = Vectors([[1],[4], [7]])
seq1 = Sequences('data' + os.sep + 'sequences1.seq')
seqn = Sequences('data' + os.sep + 'sequences2.seq')
hmc = HiddenSemiMarkov('data' + os.sep  + 'hidden_semi_markov.dat')

hvom_sample = HiddenVariableOrderMarkov(pj(path, "dupreziana21.hc"))



files = {
"belren1.seq": pj(path, "belren1.seq"),
"elstar1.seq":pj(path, "elstar1.seq"),
"fuji1.seq":pj(path, "fuji1.seq"),
"gala1.seq":pj(path, "gala1.seq"),
"granny1.seq":pj(path, "granny1.seq"),
"reinet1.seq":pj(path, "reinet1.seq"),
"dupreziana_20a2.seq": pj(path, "dupreziana_20a2.seq"),
"dupreziana_40a2.seq": pj(path, "dupreziana_40a2.seq"),
"dupreziana_60a2.seq": pj(path, "dupreziana_60a2.seq"),
"dupreziana_80a2.seq": pj(path, "dupreziana_80a2.seq"),
"dupreziana21.hc": pj(path, "dupreziana21.hc"),
"chene_sessile_15pa.seq": pj(path, "chene_sessile_15pa.seq"),
"pin_laricio_7x.seq":pj(path , "pin_laricio_7x.seq"),
"pin_laricio_6.hsc":pj(path , "pin_laricio_6.hsc"),
"abri13.ren":pj(path, 'abri13.ren')

}


_seq1 = Sequences(pj(path , 'dupreziana_20a2.seq'))
seq2 = RemoveRun(_seq1, 1, 0, "End")

seq3 = Sequences(pj(path , 'dupreziana_40a2.seq'))
seq4_0 = RemoveRun(seq3, 2, 0, "End")
seq4 = SegmentationExtract(seq4_0, 1, 2)

seq5 = Sequences(pj(path , 'dupreziana_60a2.seq'))
seq6_0 = RemoveRun(seq5, 2, 0, "End")
seq6 = LengthSelect(SegmentationExtract(seq6_0, 1, 2), 1, Mode="Reject")

seq7 = Sequences(pj(path , 'dupreziana_80a2.seq'))
seq8_0 = RemoveRun(seq7, 2, 0, "End")
seq8 = SegmentationExtract(seq8_0, 1, 2)

seq10 = Merge(seq2, seq4, seq6, seq8)
seq11 = Transcode(seq10, [0, 1, 0])


seq0 = Sequences(pj(path , "chene_sessile_15pa.seq"))
vec10 = Vectors(seq0)
vec95 = ValueSelect(vec10, 1, 95)


seq20 = Sequences(pj(path , "belren1.seq"))



_seq69 = Sequences(pj(path , "pin_laricio_7x.seq"))
seq70 = Cluster(_seq69, "Step", 1, 10)
hsmc60 = HiddenSemiMarkov(pj(path , 'pin_laricio_6.hsc'))
hsmc6 = Estimate(seq70, "HIDDEN_SEMI-MARKOV", hsmc60)
