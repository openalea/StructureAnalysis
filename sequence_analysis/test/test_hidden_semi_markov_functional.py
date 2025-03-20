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
# openalea.stat_tool.plot.DISABLE_PLOT = True
from openalea.stat_tool.plot import DISABLE_PLOT
# DISABLE_PLOT = False
DISABLE_PLOT = True

from tools import interface
from tools import runTestClass, robust_path as get_shared_data

import os

from openalea.stat_tool.output import plot, Plot

plot.DISABLE_PLOT = DISABLE_PLOT

from openalea.stat_tool.plot import get_plotter, mplotlib 

from openalea.stat_tool.distribution import set_seed

def test1():
    
    set_seed(0)
    from openalea.sequence_analysis import Estimate
    
    # TODO: find model with more separated states
    hsm = HiddenSemiMarkov(str(get_shared_data('test_hidden_semi_markov.dat')))
    # seg fault
    # hsmc_est.plot("Intensity", 1) 
    
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
    # hsmc_est.plot("Intensity", 1)

def test2():

    from openalea.sequence_analysis import Estimate
    from openalea.sequence_analysis import seq_map
    set_seed(0)
    
    # TODO: find model with more separated states
    hsm = HiddenSemiMarkov(str(get_shared_data('test_hidden_semi_markov_param.dat')))
    # seg fault
    # hsmc_est.plot("Intensity", 1) 
        
    # Simulate nb_seq with length seq_length
    nb_seq = 30
    seq_length = 100
    seq = hsm.simulation_nb_sequences(nb_seq, seq_length, True)
    assert(len(seq) == nb_seq)
    assert(len(seq[0]) == seq_length)
    # NB: hsm has 1 output process but simulation includes the hidden state
    assert(len(seq[0][1]) == 2)

    nb_states = 3
    seq_estim = seq.select_variable([2], True)
    print(seq_estim.display())

    # TODO: why are two states the same?
    # hsmc_est = Estimate(seq_estim, "HIDDEN_SEMI-MARKOV", "Ordinary", nb_states, "Irreducible", Nbiteration=300)   
    # print(hsmc_est.display())
    #
    # # TODO: find adequate error message in 
    # plotter = mplotlib()
    # # hsmc_est.plot("Intensity", 1)
    #
    # hsmc_est.extract_histogram(1,1).plot()
    # hsmc_est.extract(seq_map['Observation'],1,1).plot(Title="Observation distribution for state 1")
    # hsmc_est.extract(seq_map['Sojourn'],0,0).plot(Title="Sojourn distribution for state 0")

    hsmc_est_file = Estimate(seq_estim, "HIDDEN_SEMI-MARKOV", hsm, Nbiteration=300)   
    print(hsmc_est_file.display())
    
if __name__ == "__main__":
    # test1()
    test2()
