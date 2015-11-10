from openalea.sequence_analysis import *
from openalea.sequence_analysis.estimate import  Estimate
from openalea.sequence_analysis import get_shared_data
seq0 = Sequences(get_shared_data("chene_sessile_15pa.seq"))

Plot(seq0, ViewPoint="Data")
