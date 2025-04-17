# -*- coding: utf-8 -*-
"""tests on mv_mixture"""
__version__ = "$Id$"

from openalea.stat_tool import _stat_tool
from openalea.sequence_analysis import _sequence_analysis
from openalea.sequence_analysis.hidden_semi_markov import HiddenSemiMarkov
from openalea.sequence_analysis.simulate import Simulate
from openalea.sequence_analysis.data_transform import Thresholding
from openalea.sequence_analysis.sequences import Sequences, IndexParameterType

from openalea.stat_tool.data_transform import *
from openalea.stat_tool.cluster import Cluster
from openalea.stat_tool.cluster import Transcode, Cluster

import openalea.stat_tool.plot #import DISABLE_PLOT
# openalea.stat_tool.plot.DISABLE_PLOT = True
from openalea.stat_tool.plot import DISABLE_PLOT
DISABLE_PLOT = False

import os

from openalea.stat_tool.output import plot, Plot

plot.DISABLE_PLOT = DISABLE_PLOT

from openalea.stat_tool.plot import get_plotter, mplotlib 
from openalea.stat_tool.distribution import set_seed

def test1():
    
    from pathlib import Path
    from openalea.sequence_analysis import _MarkovianSequences
    
    data_path = Path(openalea.sequence_analysis.__path__[0])
    data_path = str(Path.joinpath(data_path.parent.parent.parent.absolute(), "share","data"))   
    model_file = "switching_lmm_irred.hsc"
    
    hsm = HiddenSemiMarkov(data_path + os.sep + model_file)
    print(hsm.display())
    from openalea.sequence_analysis import Simulate
    nb_seq = 30
    seq_length = 100
    set_seed(0)
    seq = hsm.simulation_nb_sequences(nb_seq, seq_length, True)
    assert(len(seq) == nb_seq)
    assert(len(seq[0]) == seq_length)
    # NB: hsm has 1 output process but simulation includes the hidden state
    assert(len(seq[0][1]) == 2)
    seqm = _MarkovianSequences(seq)
    from openalea.sequence_analysis import Estimate
    hsmd = hsm.simulation_nb_sequences(nb_seq, seq_length, True)
    seq_estim = hsmd.select_variable([2], True)
    from openalea.stat_tool import NegativeBinomial
    d = NegativeBinomial(1, 10, 0.5)
    index = []
    for u in range(nb_seq):
        indexut = []
        for t  in range(seq_length):
            indexut.append([d.simulation()])
        index.append(indexut)
    seq_index = Sequences(index)    
    seq_index_markov = _MarkovianSequences(seq_index)
    seq_merge = seq_estim.merge_variable([seq_index_markov], 1)
    try:
        seq_index = seq_merge.set_variable_as_index_parameter(2, "TIME")
    except:
        pass
    else:
        raise error.FormatError("Failed to raise exception for bad index values")
    
    index = []
    for u in range(nb_seq):
        indexut = []
        indexut.append([d.simulation()])
        for t  in range(1, seq_length):
            indexut.append([indexut[-1][-1]+d.simulation()])
        index.append(indexut)
    seq_index = Sequences(index)    
    seq_index_markov = _MarkovianSequences(seq_index)
    seq_merge = seq_estim.merge_variable([seq_index_markov], 1)
    seq_index = seq_merge.set_variable_as_index_parameter(2, "TIME")
    seq_index_markov2 = _MarkovianSequences(seq_index)
    assert IndexParameterType(seq_index_markov2)=='TIME'
    for u in range(nb_seq):
        for t  in range(0, seq_length):
            assert(int(seq_index_markov2.get_index_parameter(u,t)) == index[u][t][0])
            
    # Resimulate with new index
    hsmd = hsm.semi_markov_switching_lm_simulation(1, seq_index_markov2, openalea.stat_tool.I_DEFAULT, True);
    assert(not(hsmd is None))
    
    # Plot index and simulated values per state
    values = {}
    for u in range(hsmd.nb_sequence):
        for t in range(len(hsmd[u])):
            state = hsmd[u][t][0]
            if state in values.keys():
                values[state][0] += [hsmd.get_index_parameter(u,t)] # predictor
                # values[state][0] += [index[u][t][0]] # predictor
                values[state][1] += [hsmd[u][t][1]] # simulated value
            else:
                values[state] = [[],[]]
    nb_states = len(values.keys())
    import matplotlib.pyplot as plt       
    import numpy as np

    fig, subfigs = plt.subplots(int(np.ceil(nb_states / 2)), 2) 
    for k in values.keys():    
        subfigs[int(k/2)][k%2].plot(values[k][0],values[k][1],'o')
    plt.show()
    seq_estim = hsmd.select_variable([2], True);

    hsmc_est_file = Estimate(seq_estim, "HIDDEN_SEMI-MARKOV", hsm, Nbiteration=300)   
    

if __name__ == "__main__":
   test1()
