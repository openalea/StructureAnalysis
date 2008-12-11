# comparison of the outputs for the Sequences and Trees classes
import sys, os
import openalea.stat_tool as stat_tool
import openalea.tree_statistic.trees as trees
import openalea.tree_statistic.hmt as hmt
import openalea.aml as aml
nb_children=1
size=100
ident=stat_tool.DistributionIdentifier.UNIFORM
parameter=stat_tool.D_DEFAULT
distrib= stat_tool._ParametricModel(ident, 1, 1, parameter, 1)
# name of 
hmotrefpath= "./hmot_np.hmt";
hmotinitpath= "./hmot_np_init.hmt";
hmcinitpath= "./hmc_init.hmc";
HMT=hmt.HiddenMarkovTree(hmotrefpath)
HMTI=hmt.HiddenMarkovTree(hmotinitpath)
HMCI=aml.HiddenSemiMarkov(hmcinitpath)
T=HMT.Simulate(2, size, nb_children)
# delete state variable
T=T.SelectVariable(1, "Keep")
# build sequence from T
S=T.BuildSequences()

# Print Sequences object
print S
# Print Trees object
print T

# Display sequence object
aml.Display(S)
# Display forest object
T.Display()

# Display sequence object and data
aml.Display(S, ViewPoint="Data")
# Display forest object and data
T.Display(ViewPoint="Data")

# Estimate a hidden Markov tree
hmt=T.Estimate("HIDDEN_MARKOV_TREE", HMTI, 10)
# Estimate a hidden Markov chain
hmc=aml.Estimate(S, "HIDDEN_SEMI-MARKOV", HMCI, NbIteration=10)

# Print an estimated HMC
print hmc
# Print an estimated HMT
print hmt

# Display an estimated HMT
hmt.Display()
# Display an estimated HMC
aml.Display(hmc)

# Plot an estimated HMT
# hmt.Plot()
# Plot an estimated HMC
