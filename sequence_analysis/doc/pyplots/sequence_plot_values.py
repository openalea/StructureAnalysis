from openalea.sequence_analysis import *
from openalea.sequence_analysis.data import path
seq = Sequences(path +  'sequences_tutorial.dat')
Plot(seq, ViewPoint="Data")

Plot(ExtractHistogram(seq, "Value"))

Plot(ExtractHistogram(seq, "Length"))

