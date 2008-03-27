# a test for the class hmt.Hmt: constructor and basic methods
import sys, os
import openalea.stat_tool as stat_tool
import openalea.tree_statistic.trees as trees, openalea.tree_statistic.hmt as hmt
inf_bound=1
sup_bound=3
probability= 0.6
ident=stat_tool.DistributionIdentifier.UNIFORM
parameter=stat_tool.D_DEFAULT
distrib= stat_tool.Parametric(ident, inf_bound, sup_bound, parameter, probability)
# Distribution used for the number of children and the tree attributes
print distrib
max_depth=3
max_size=10
nbtrees=4
# defining a set of trees
tree_list=[]
tv=trees.TreeValue([1., 0])
R=trees.TreeStructure(distrib, max_size, max_depth)
tmp_tree=trees.Tree(tv, R)
n=3
tree_list.append(trees.Tree(tmp_tree))
while n < nbtrees:
    n=n+1
    R.Simulate(distrib, max_size, max_depth)
    tmp_tree=trees.Tree(tv, R)
    tree_list.append(trees.Tree(tmp_tree))    
distrib_list=[]
for i in range(tmp_tree.NbInt()):
    distrib_list.append(distrib)
    
for n in range(len(tree_list)):
    tree_list[n].Simulate(distrib_list)
# initializing a Trees object
T_ind=trees.Trees(tree_list)
H=hmt.HiddenMarkovTree("hmot.hmt")
sample_size=2
tree_size=20
nb_children=3
print 'A hidden Markov tree H read from file "hmot.hmt":'
print H
H.Save("hmot_ascii.hmt", "ASCII", True)
I=hmt.HiddenMarkovTree("hmot_ascii.hmt")
I.Save("hmot_ascii.hmt", "ASCII", True)
J=hmt.HiddenMarkovTree("hmot_ascii.hmt")
if str(I)==str(J):
    print "load o (save o load) == load"
else:
    print "load o (save o load) != load"
    print I
    print J
print "Saving H as a spreadsheet"
I.Save("hmot_spreadsheet.hmt", "Spreadsheet", True)
print "Simulation of the labels using H, tree size and number of children:"
T=H.Simulate(sample_size, tree_size, nb_children)
# print T
print "Simulation of the labels using H, tree size and a Trees object:"
T=H.Simulate(sample_size, T_ind)
Histo_size=T_ind.ExtractHistogram("Size")
Histo_nb_children=T_ind.ExtractHistogram("NbChildren")
T=H.Simulate(Histo_size, Histo_nb_children)
print "Simulated labels:"
print T
for t in range(T.NbTrees()):
    T.Tree(t).Display()
print "Type of the attributes:", T.Types()
print "Copy a hmt"
HC=hmt.HiddenMarkovTree(H)
print HC