# -*- coding: utf-8 -*-
# a test for the class hmt.Hmt: estimation - syntax
import sys, os
import openalea.stat_tool as stat_tool
import openalea.tree_statistic.trees as trees, openalea.tree_statistic.hmt as hmt
inf_bound = 1
sup_bound = 3
distrib = stat_tool.Uniform(inf_bound, sup_bound)
# Distribution used for the number of children and the tree attributes
file_name = "hmot_np_2s.hmt"
# read a HMT from a file
print 'A hidden Markov tree H read from file "', file_name, '"'
H = hmt.HiddenMarkovIndOutTree(file_name)
H2O = hmt.HiddenMarkovIndOutTree("hmot_np_2o.hmt")
sample_size = 2
tree_size = 30
nb_children = 2
H.Display()
# simulate labels from this HMT
print "Label simulation using H (1st tree): "
T = H.Simulate(sample_size, tree_size, nb_children)
print T.Tree(0)
# delete the state variable
print "Delete the state variable and the smoothed probabilities."
T = T.SelectVariable(1, "Keep")
print "State variable deleted: "
print T.Tree(0)
# parameter estimation from the simulated trees
# initialization from a model
print "Parameter estimation from the simulated trees" + \
" using the initial model specification:"
EH = T.Estimate("HIDDEN_MARKOV_TREE", H, 20)
EH.Display()
# checking exceptions raised by Estimate
print "Parameter estimation using a initial HMT "+ \
"with bad number of output processes:"
try:
    EH = T.Estimate("HIDDEN_MARKOV_TREE", H2O, 20)
except trees.StatTreeError, f:
    print f
else:
    print "Failed to raise exception for bad number of output processes"
print "Parameter estimation with bad type of algorithm:"
try:
    EH = T.Estimate("HIDDEN_MARKOV_TREE", H, 20, Algorithm="Saem")
except ValueError, f:
    print f
else:
    print "Failed to raise exception for bad type of algorithm"
print "Parameter estimation with bad value for Saem exponent:"
try:
    EH = T.Estimate("HIDDEN_MARKOV_TREE", H, 20, Algorithm="GibbsSampling", Saem=2.)
except trees.StatTreeError, f:
    print f
else:
    print "Failed to raise exception for bad value of Saem"
EH.Display()
# parameter estimation from the simulated trees
# initialization from a self-transition probability
print "Parameter estimation from the simulated trees" + \
" using the self-transition probability specification:"
EH = T.Estimate("HIDDEN_MARKOV_TREE", 2, "Irreductible", 0.8, 20)
EH.Display()
print "Parameter estimation from the simulated trees" + \
" using self-transition specification and EM 'a la Gibbs':"
EH = T.Estimate("HIDDEN_MARKOV_TREE", 2, "Irreducible", 0.8, 20, Algorithm="GibbsSampling",
              Saem=0.)
EH.Display()
print "Parameter estimation from the simulated trees" + \
" using self-transition specification and CEM:"
EH = T.Estimate("HIDDEN_MARKOV_TREE", 2, "Irreducible", 0.8, 20, Algorithm="Viterbi",
              Saem=0.)
EH.Display()
print "Parameter estimation from the simulated trees" + \
" using self-transition specification and SAEM:"
EH = T.Estimate("HIDDEN_MARKOV_TREE", 2, "Irreducible", 0.8, 20, Algorithm="ForwardBackwardSampling",
              Saem=0.5)
EH.Display()
# Extract the trees from model
print "Extract the data part of the initial (file) HMT, with no data:"
try:
    HMTD = H.ExtractData()
except trees.StatTreeError, f:
    print f
else:
    print "Failed to raise exception for no data part in HiddenMarkovIndOutTree"
    print HMTD
print "Extract the data part of the estimated HMT"
HMTD = EH.ExtractData()
if (str(HMTD)!=str(T)):
    print "HMTD:"
    print HMTD
    print "T:"
    print T
else:
    print "Result is identical with the simulated data"
    HMTD.Display()
# Extract the state trees
print "Computation of the state trees"
S = T.ComputeStateTrees(H, "Viterbi")
# Copy the state trees
CPS = hmt.HiddenMarkovTreeData(S)
S.Display(Detail=1)
HISTO = S.ExtractHistogram("Value", 1)
print "Marginal distribution for variable 1:"
print HISTO
print "Observation histogram for variable 1 and state 1:"
OBS1 = S.ExtractHistogram("Observation", 1, 1)
print OBS1
# EH.Plot("Observation", variable=0)
# S.Plot("Observation", variable=0)
# State profiles
EH.Display(ViewPoint="StateProfile", TreeId=1, NbStateTrees=4)
# compute entropy values
print EH.EntropyComputation()
print EH.EntropyComputation(0)
# parameter estimation from the simulated trees
# initialization from a self-transition probability
# forcing parametric estimation
print "Parameter estimation from the simulated trees" + \
" using the self-transition probability specification" + \
" forcing parametric estimation:"
EH = T.Estimate("HIDDEN_MARKOV_TREE", 2, "Irreductible", 0.999, 20,
              ForceParametric=[True])
EH.Display()
# Re-estimate on segmented tree
print "Estimate an HMT on state tree: "
SS = S.SelectVariable([0])
SS.ToIntType(0)
EH2 = SS.Estimate("HIDDEN_MARKOV_TREE", 2, "Irreductible", 0.999, 20)
EH2.Display()

# HMT with discrete parametric observations
print "Estimation of an HMT with parametric discrete observations"
HP = hmt.HiddenMarkovIndOutTree("hmot_param.hmt")
TP = HP.Simulate(sample_size*10, tree_size, nb_children)
TP.ToIntType(0)
EHP = TP.Estimate("HIDDEN_MARKOV_TREE", 2, "Irreductible", 0.999, 20)
SP = TP.ComputeStateTrees(EHP, "Viterbi")
OBSP1 = SP.ExtractHistogram("Value", 3)
#EHP.Display()

