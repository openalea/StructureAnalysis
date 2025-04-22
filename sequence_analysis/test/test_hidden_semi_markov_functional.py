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
DISABLE_PLOT = False
# DISABLE_PLOT = True

from tools import interface
from tools import runTestClass, robust_path as get_shared_data

import os

from openalea.stat_tool.output import plot, Plot

plot.DISABLE_PLOT = DISABLE_PLOT

from openalea.stat_tool.plot import get_plotter, mplotlib 

from openalea.stat_tool.distribution import set_seed

def test1():
    """Estimate HSMC with nonparametric emission distributions"""
    set_seed(0)

    # TODO: find model with more separated states
    hsm = HiddenSemiMarkov(str(get_shared_data('test_hidden_semi_markov.dat')))
    # seg fault
    hsm.plot("Intensity", 1)  
    hsm.plot("Observation", 1)
    hsm.plot("Counting", 1)
    hsm.plot("Recurrence", 1)
    hsm.plot("Sojourn", 1)
    try:
        hsm.plot("InitialRun", 1)
    except:
        pass
    else:
        raise RuntimeError("Failed to raise error")
    try:
        hsm.plot("FinalRun", 1)
    except:
        pass
    else:
        raise RuntimeError("Failed to raise error")            
    hsm.plot("FirstOccurrence", 1)
    hsm.plot()

    # Simulate nb_seq with length seq_length
    nb_seq = 30
    seq_length = 100
    seq = hsm.simulation_nb_sequences(nb_seq, seq_length, True)
    assert(len(seq) == nb_seq)
    assert(len(seq[0]) == seq_length)
    # NB: hsm has 1 output process but simulation includes the hidden state
    assert(len(seq[0][1]) == 2)

    print(seq[0])
    obs = seq.select_variable([1], keep=False)

    obs.plot("Intensity", 0)

    from openalea.sequence_analysis import Estimate
    nb_states = 8

    hsmc_est = Estimate(obs, "HIDDEN_SEMI-MARKOV", "Ordinary", nb_states, "LeftRight", Nbiteration=300)   
    print(hsmc_est.display())

    # TODO: find adequate error message in 
    # hsmc_est = Estimate(seq, "HIDDEN_SEMI-MARKOV", "Ordinary", nb_states, "LeftRight", Nbiteration=300)
    plotter = mplotlib()
    # seg fault
    hsmc_est.plot("Intensity", 1)
    hsmc_est.plot("Observation", 1)
    hsmc_est.plot("Counting", 1)
    hsmc_est.plot("Recurrence", 1)
    hsmc_est.plot("Sojourn", 1)    

    from openalea.sequence_analysis import seq_map
    hsmc_est.extract(seq_map['Observation'],1,1).plot(Title="Observation distribution for state 1")
    hsmc_est.extract(seq_map['Sojourn'],0,0).plot(Title="Sojourn distribution for state 0")


def test2():
    """Estimate HSMC with parametric emission distributions"""

    from openalea.sequence_analysis import Estimate
    from openalea.sequence_analysis import seq_map
    
        
    hsm = HiddenSemiMarkov(str(get_shared_data('test_hidden_semi_markov_param.dat')))
    
    # TODO: find adequate error message
    # hsm.plot("Intensity", 1)
    hsm.plot("Intensity", 0)   
    hsm.plot("Observation", 1) 
        
    # Simulate nb_seq with length seq_length
    nb_seq = 30
    seq_length = 500
    seq = hsm.simulation_nb_sequences(nb_seq, seq_length, True)
    assert(len(seq) == nb_seq)
    assert(len(seq[0]) == seq_length)
    # NB: hsm has 1 output process but simulation includes the hidden state
    assert(len(seq[0][1]) == 2)

    nb_states = 3
    obs = seq.select_variable([2], True)

    set_seed(0)
    hsmc_est = Estimate(obs, "HIDDEN_SEMI-MARKOV", "Ordinary", nb_states, "Irreducible", Nbiteration=300)   
    print(hsmc_est.display())
    
    hsmc_est.plot("Observation", 1)
    
    from openalea.sequence_analysis import seq_map
    hsmc_est.extract(seq_map['Observation'],1,0).plot(Title="Observation distribution for state 0")
    hsmc_est.extract(seq_map['Observation'],1,1).plot(Title="Observation distribution for state 1")
    hsmc_est.extract(seq_map['Observation'],1,2).plot(Title="Observation distribution for state 2")
    
if __name__ == "__main__":
    test1()
    test2()
