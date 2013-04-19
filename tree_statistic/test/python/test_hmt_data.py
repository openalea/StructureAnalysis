# -*- coding: utf-8 -*-
# a test for the class hmt.Hmt: estimation - syntax
import sys, os
import openalea.stat_tool as stat_tool
from openalea.stat_tool.plot import gnuplot, set_plotter
# set_plotter(gnuplot())
import openalea.tree_statistic.trees as trees, openalea.tree_statistic.hmt as hmt
inf_bound=0
sup_bound=3
probability= 0.6
ident=stat_tool.DistributionIdentifierType.UNIFORM
parameter=stat_tool.D_DEFAULT
distrib= stat_tool.distribution._DiscreteParametricModel(ident, inf_bound, sup_bound, parameter, probability)
# Distribution used for the number of children and the tree attributes
file_name="hmot_np_2s.hmt"
# read a non parametric HMT from a file
print 'A hidden Markov tree H read from file "', file_name, '"'
H=hmt.HiddenMarkovIndOutTree(file_name)
sample_size=10
tree_size=30
nb_children=3
print H
# simulate labels from this HMT
print "Label simulation using H: "
T=H.Simulate(sample_size, tree_size, nb_children)
HISTO1=T.ExtractHistogram("Value", 1)
print "Marginal distribution for variable 1:"
print HISTO1
print "Observation histogram for variable 1 and state 0:"
OBS1=T.ExtractHistogram("Observation", 1, 0)
OBS1.display(Detail=2)
# OBS1.plot()
file_name="hmot_param.hmt"
# read a non parametric HMT from a file
print 'A hidden Markov tree H read from file "', file_name, '"'
H=hmt.HiddenMarkovIndOutTree(file_name)
sample_size=20
tree_size=30
nb_children=3
print H
# simulate labels from this HMT
print "Label simulation using H: "
T2=H.Simulate(sample_size, tree_size, nb_children)
T2.Display()
HISTO2=T2.ExtractHistogram("Value", 1)
print "Marginal distribution for variable 1:"
print HISTO2
print "Observation histogram for variable 1 and state 1:"
OBS2=T2.ExtractHistogram("Observation", 1, 1)
print OBS2.display()
# OBS2.plot()
# T2.Plot("Observation", 1)
# Copy constructor of HiddenMarkovIndOutTreeData
CPT2 = hmt.HiddenMarkovTreeData(T2)
OBS3 = CPT2.ExtractHistogram("Observation", 1, 1)

ECPT2 = T.ExtractMarkov()

HIST3 = T2.ExtractHistogram("Value", variable=2)
print HIST3.display(Detail=2)
HIST3.plot()
# T2.Plot("FirstOccurrenceRoot", variable=0)
T2.MPlot("FirstOccurrenceRoot", variable=0)

