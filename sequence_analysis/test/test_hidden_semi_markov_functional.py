# -*- coding: utf-8 -*-
"""tests on mv_mixture"""
__version__ = "$Id$"

from openalea.stat_tool import _stat_tool
from openalea.sequence_analysis import _sequence_analysis
from openalea.sequence_analysis.hidden_semi_markov import HiddenSemiMarkov
from openalea.sequence_analysis.simulate import Simulate
from openalea.sequence_analysis.data_transform import Thresholding

from openalea.stat_tool.data_transform import *
from openalea.stat_tool.cluster import Cluster
from openalea.stat_tool.cluster import Transcode, Cluster

import openalea.stat_tool.plot #import DISABLE_PLOT
openalea.stat_tool.plot.DISABLE_PLOT = True

from tools import interface
from tools import runTestClass, robust_path as get_shared_data

import os

from openalea.stat_tool.output import plot, Plot

from openalea.stat_tool.plot import get_plotter, mplotlib 

def test1():

    from openalea.sequence_analysis import Estimate
    
    # TODO: find model with more separated states
    hsm = HiddenSemiMarkov(str(get_shared_data('test_hidden_semi_markov.dat')))

    # Simulate nb_seq with length seq_length
    nb_seq = 30
    seq_length = 100
    seq = hsm.simulation_nb_sequences(nb_seq, seq_length, True)
    assert(len(seq) == nb_seq)
    assert(len(seq[0]) == seq_length)
    # NB: hsm has 1 output process but simulation includes the hidden state
    assert(len(seq[0][1]) == 2)

    nb_states = 8
    seq_estim = seq.select_variable([2], True)

    hsmc_est = Estimate(seq_estim, "HIDDEN_SEMI-MARKOV", "Ordinary", nb_states, "LeftRight", Nbiteration=300)   
    print(hsmc_est.display())

    # TODO: find adequate error message in 
    # hsmc_est = Estimate(seq, "HIDDEN_SEMI-MARKOV", "Ordinary", nb_states, "LeftRight", Nbiteration=300)
    plotter = mplotlib()
    hsmc_est.plot("Intensity", 1)

if __name__ == "__main__":
    test1()
