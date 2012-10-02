# -*- coding: utf-8 -*-
# a test for the class hmt.Hmt: estimation accurracy
import sys, os
import openalea.stat_tool as stat_tool
import openalea.tree_statistic.trees as trees, openalea.tree_statistic.hmt as hmt
inf_bound=1
sup_bound=3
probability= 0.6
ident=stat_tool.DistributionIdentifierType.UNIFORM
parameter=stat_tool.D_DEFAULT
distrib= stat_tool.distribution._DiscreteParametricModel(ident, inf_bound, sup_bound, parameter, probability)
# Distribution used for the number of children and the tree attributes
H=hmt.HiddenMarkovIndOutTree("hmot.hmt")
sample_size=200
tree_size=30
nb_children=2
print 'A hidden Markov tree H read from file "hmot.hmt":'
H.Display()
T=H.Simulate(sample_size, tree_size, nb_children)
print "Label simulation using H: "
T.Display()
print "Delete the state variable"
T=T.SelectVariable([1, 2])
T.Display()
print "Parameter estimation from the simulated trees " + \
"using initial model specification:"
EH=T.Estimate("HIDDEN_MARKOV_TREE", H, 100)
EH.Display()
# raise UserWarning, "Warning: need to be debugged"
print "Parameter estimation from the simulated trees " + \
"using state number and self transition specifications:"
EH=T.Estimate("HIDDEN_MARKOV_TREE", 2, "IRREDUCIBLE", 0.99, 100)
EH.Display()
