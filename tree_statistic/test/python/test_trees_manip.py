# -*- coding: utf-8 -*-
# a test for the class trees.Trees: tree manipulation (merge, cluster, etc.)
import sys, os
import openalea.stat_tool as stat_tool
import openalea.tree_statistic.trees as trees
inf_bound = 1
sup_bound = 3
probability = 0.6
children_distrib = stat_tool.Uniform(inf_bound, sup_bound)
# distribution of the children
attributes_distrib = stat_tool.Uniform(0, 10)
# distribution of the attributes
max_depth = 3
max_size = 10
nbtrees = 2
# defining a set of trees
tree_list = []
# simulation of the structure
R = trees.TreeStructure(children_distrib, max_size, max_depth)
tmp_tree = trees.Tree([1., 0], R)
n = 1
tree_list.append(trees.Tree(tmp_tree))
while n < nbtrees:
    n = n+1
    R.Simulate(children_distrib, max_size, max_depth)
    tmp_tree = trees.Tree([1., 0], R)
    tree_list.append(trees.Tree(tmp_tree))

distrib_list = []
# simulation of the labels
for i in range(tmp_tree.NbInt()):
    distrib_list.append(attributes_distrib)

for n in range(len(tree_list)):
    tree_list[n].Simulate(distrib_list)

# initialize a Trees object
T = trees.Trees(tree_list)
print T
# T.Plot()
# T.Display()
# extract histograms
print "Histogram of the tree size: "
histo = T.ExtractHistogram("Size")
print(histo.display(Detail=2))
print "Histogram of the number of children: "
histo = T.ExtractHistogram("NbChildren")
print(histo.display(Detail=2))
print "Histogram of marginal distribution: "
histo = T.ExtractHistogram("Value", 1)
print(histo.display(Detail=2))
# save trees
T.Save("fake_mtg_forest.txt", True)
# merge trees
# several (2) sets of trees
trees_list = []
t = 0
while t < 2:
    t = t+1
    tree_list=[]
    n = 0
    while n < 2:
        n = n+1
        R.Simulate(attributes_distrib, max_size, max_depth)
        tmp_tree = trees.Tree([1., 0], R)
        tree_list.append(trees.Tree(tmp_tree))
    for n in range(len(tree_list)):
        tree_list[n].Simulate(distrib_list)
    trees_list.append(trees.Trees(tree_list))

V = T.Merge(trees_list)
print "Reference tree: "
for t in range(T.NbTrees()):
    print "Tree number ", t, ": "
    print T.Tree(t)

print "Number of other trees: 4 (2 Trees objects)"
print "Number of trees merged into V: ", V.NbTrees()
equal = True
if str(V.Tree(V.NbTrees()-1))!=str(trees_list[1].Tree(n)):
    msg="Last merged trees do not match: \n"
    msg+=str(V.Tree(V.NbTrees()-1))+"\n"
    msg+=str(trees_list[1].Tree(n))+"\n"
else:
    print "Last merged trees do match!"
# merge variables
# several (2) sets of trees
trees_list=[]
for n in range(2):
    tree_list=[]
    for t in range(T.NbTrees()):
        structure = trees.TreeStructure(T.Tree(t))
        tmp_tree = trees.Tree([0], structure)
        tree_list.append(trees.Tree(tmp_tree))
        tree_list[t].Simulate(distrib_list)
    trees_list.append(trees.Trees(tree_list))

V = T.MergeVariable(trees_list)
print "Reference tree: "
print T.Tree(0)
print "Processes to be merged:"
for tl in trees_list:
    print tl.Tree(0)
print "Merged variables:"
print V.Tree(0)
rootval = eval(str(T.Tree(0).Get(0)))
for tl in trees_list:
    rootval+=(eval(str(tl.Tree(0).Get(0))))
if str(V.Tree(0).Get(0))!=str(rootval):
    msg="Merged variables do not match for tree 0 at root node: \n"
    msg+=str(V.Tree(0).Get(0))+"\n"
    msg+=str(rootval)+"\n"
    raise Warning, msg
else:
    print "Merged variables do match for tree 0 at root node!"

# check exceptions raised by ExtractHistogram
try:
    T.ExtractHistogram("BadFeature")
except ValueError, v:
    print v
else:
    print "Failed to raise exception for bad feature histogram."
try:
    T.ExtractHistogram("Value", 0)
except trees.StatTreeError, e:
    print e
else:
    print "Failed to raise exception for bad variable for marginal histogram."
try:
    T.ExtractHistogram("Value", 2)
except IndexError, e:
    print e
else:
    print "Failed to raise exception for bad variable for marginal histogram."
# check exceptions raised by Merge
try:
    T.Merge([])
except UserWarning, w:
    print w
else:
    print "Failed to raise exception for empty tree list."
n = 0
tree_list=[]
while n < 2:
    n = n+1
    R.Simulate(attributes_distrib, max_size, max_depth)
    tmp_tree = trees.Tree([0., 0, 0], R)
    tree_list.append(trees.Tree(tmp_tree))

trees_list=[]
trees_list.append(trees.Trees(tree_list))
try:
    V = T.Merge(trees_list)
except trees.StatTreeError, e:
    print e
else:
    print "Failed to raise exception for bad number of variables."
n = 0
tree_list=[]
while n < 2:
    n = n+1
    R.Simulate(attributes_distrib, max_size, max_depth)
    tmp_tree = trees.Tree([0, 0], R)
    tree_list.append(trees.Tree(tmp_tree))
trees_list=[]
trees_list.append(trees.Trees(tree_list))
try:
    V = T.Merge(trees_list)
except trees.StatTreeError, e:
    print e
else:
    print "Failed to raise exception for bad variable type."
try:
    V = T.Merge([0, w])
except TypeError, e:
    print e
else:
    print "Failed to raise exception for bad argument type."
# check exceptions raised by MergeVariable
try:
    M = T.MergeVariable([])
except UserWarning, w:
    print w
else:
    print "Failed to raise exception for empty tree list."
    print "Result: ", M
    print M.Tree(0)
trees_list=[]
t = 0
while t < 2:
    t = t+1
    tree_list=[]
    n = 0
    while n < 2:
        n = n+1
        R.Simulate(attributes_distrib, max_size, max_depth)
        tmp_tree = trees.Tree([1., 0], R)
        tree_list.append(trees.Tree(tmp_tree))
    for n in range(len(tree_list)):
        tree_list[n].Simulate(distrib_list)
    trees_list.append(trees.Trees(tree_list))
try:
    V = T.MergeVariable(trees_list)
except trees.StatTreeError, e:
    print e
else:
    print "Failed to raise exception for bad structure."
try:
    V = T.MergeVariable([0, T])
except TypeError, e:
    print e
else:
    print "Failed to raise exception for bad argument type."
# cluster variables:
# mode "Limit"
limits=[5, 10]
print "Cluster variable 1 in T with limits ", limits
V = T.Cluster("Limit", 1, limits)
for t in range(V.NbTrees()):
    print "Tree number ", t, ": "
    print V.Tree(t)
# mode "Step"
step = 3
print "Cluster variable 1 in T with step ", step
V = T.Cluster("Step", 1, step)
for t in range(V.NbTrees()):
    print "Tree number ", t, ": "
    print V.Tree(t)
limits=[]
for c in range((10/step)):
    limits.append(step*(c+1))
limits.append(11)
print limits
V2 = T.Cluster("Limit", 1, limits)
equal = True
for t in range(V.NbTrees()):
    if str(V.Tree(t))!=str(V2.Tree(t)):
        equal = False
        print "Modes 'Limit' and 'Step' do not match on tree ", t
        print V.Tree(V.NbTrees()-1)
        print trees_list[1].Tree(n)
if equal:
    print "Modes 'Limit' and 'Step' do match!"
# check exceptions raised by Cluster
try:
    T.Cluster("Limit", 1, [0, 1])
except trees.StatTreeError, e:
    print e
else:
    print "Failed to raise exception for non positive limit."
try:
    T.Cluster("Limit", 1, range(12))
except trees.StatTreeError, e:
    print e
else:
    print "Failed to raise exception for number of clusters too high."
try:
    V = T.Cluster("Limit", 1, [5])
except trees.StatTreeError, e:
    print e
else:
    print "Failed to raise exception for number of clusters too low."
    print V
try:
    T.Cluster("Limit", 0, limits)
except trees.StatTreeError, e:
    print e
else:
    print "Failed to raise exception for bad variable."
try:
    T.Cluster("Step", 1, 0)
except trees.StatTreeError, e:
    print e
else:
    print "Failed to raise exception for non positive step."
try:
    T.Cluster("Equal", 1, 12)
except ValueError, e:
    print e
else:
    print "Failed to raise exception for incorrect mode."
try:
    T.Cluster("Step", 1, [5, 9])
except TypeError, e:
    print e
else:
    print "Failed to raise exception for incorrect step type."
try:
    T.Cluster("Limit", 1, [5, V2])
except TypeError, e:
    print e
else:
    print "Failed to raise exception for incorrect limit type."
try:
    T.Cluster("Limit", 1, V2)
except TypeError, e:
    print e
else:
    print "Failed to raise exception for incorrect limit type."
try:
    T.Cluster("Limit", 1, 3)
except TypeError, e:
    print e
else:
    print "Failed to raise exception for incorrect limit type."
try:
    T.Cluster("Step", -1, 3)
except IndexError, e:
    print e
else:
    "Failed to raise exception for bad variable index."
# Transcode variables:
# mode "Limit"
new_values=[0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5]
print "Transcode variable 1 in T with new values ", new_values
V = T.Transcode(1, new_values)
for t in range(V.NbTrees()):
    print "Tree number ", t, ": "
    print V.Tree(t)
# check exceptions raised by Transcode
try:
    T.Transcode(1, [0, 1])
except ValueError, e:
    print e
else:
    print "Failed to raise exception for bad number of classes."
try:
    V = T.Transcode(1, range(0, 5) + range(6,13))
except trees.StatTreeError, e:
    print e
else:
    print "Failed to raise exception for non-consecutive symbols."
    for t in range(V.NbTrees()):
        print "Tree number ", t, ": "
        print V.Tree(t)
try:
    V = T.Transcode(1, [])
except ValueError, e:
    print e
else:
    print "Failed to raise exception for nonpositive number of classes."
    print V
try:
    T.Transcode(0, range(0, 3))
except IndexError, e:
    print e
else:
    print "Failed to raise exception for bad variable."
try:
    T.Transcode(1, "Step")
except TypeError, e:
    print e
else:
    print "Failed to raise exception for incorrect array type."
try:
    T.Transcode(1, range(0,10)+ [T])
except TypeError, e:
    print e
else:
    print "Failed to raise exception for incorrect symbol type."
try:
    T.Transcode( -1, range(0, 11))
except IndexError, e:
    print e
else:
    "Failed to raise exception for bad variable index."
# select variables:
# mode "Keep"
variables=[1]
print "Select variable 1 in T"
V = T.SelectVariable(variables)
print V
for t in range(V.NbTrees()):
    print "Tree number ", t, ": "
    print V.Tree(t)
# mode "Reject"
variables = 1
print "Reject variable 1 in T"
V = T.SelectVariable(variables, "Reject")
for t in range(V.NbTrees()):
    print "Tree number ", t, ": "
    print V.Tree(t)
# check exceptions raised by SelectVariable
try:
    T.SelectVariable([0, 1], "Reject")
except trees.StatTreeError, e:
    print e
else:
    print "Failed to raise exception for rejecting all variables."
# shift
# integral variable
p = 2
print "Shift tree with param " + str(p) +" using variable 1"
V = T.Shift(variable=1, param=p)
for t in range(V.NbTrees()):
    print "Tree number ", t, ": "
    print V.Tree(t)
p = 0.5
# floating variable
print "Shift tree with param " + str(p) +" using variable 0"
V = T.Shift(variable=0, param=p)
for t in range(V.NbTrees()):
    print "Tree number ", t, ": "
    print V.Tree(t)

# Extract vectors and sequences
print "Extract non-redundant sequences"
print str(T.Size()) + " vertices"
Seq = T.BuildPySequences(False)
print str(Seq.build_vectors(True).nb_vector) + " items in sequence"
print "Extract vectors"
Vec = T.BuildVectors()
print str(Vec.nb_vector) + " vectors"
print "Sequences :"
for i in range(len(Seq)):
    print Seq[i]




