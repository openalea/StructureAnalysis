from openalea.sequence_analysis import *
from openalea.sequence_analysis.estimate import  Estimate
from openalea.sequence_analysis import shared_data_path as path
from os.path import join as pj
seq0 = Sequences(pj(path ,"chene_sessile_15pa.seq"))

Plot(seq0, ViewPoint="Data")
