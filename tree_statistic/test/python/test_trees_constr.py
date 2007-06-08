# a test for the class trees.Trees: constructor and basic methods
import sys
import os
AMLLIBDIR=os.getenv("HOME")+"/devlp/AMAPmod/build-linux/lib"
sys.path.append(AMLLIBDIR)
map(lambda x:sys.path.append(os.getenv("AMAPDIR")+"/STAT_TREES/Python/"+x),
    ["Int_fl_containers/", "Trees/", "Hmt/"])
sys.path.append(os.getenv("AMAPDIR")+"/STAT/Python/")
import trees
import stat_tools
inf_bound=1
sup_bound=3
probability= 0.6
ident=stat_tools.DistributionIdentifier.UNIFORM
parameter=stat_tools.D_DEFAULT
distrib= stat_tools.Parametric(ident, inf_bound, sup_bound, parameter, probability)
print "Distribution used for the number of children and the tree attributes:"
print distrib
max_depth=3
max_size=10
nbtrees=40
# define a set of trees
tree_list=[]
tv=[1., 0, 1, 2.] # trees.TreeValue([1., 0])
R=trees.TreeStructure(distrib, max_size, max_depth)
tmp_tree=trees.Tree(tv, R)
n=1
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
# initialize a Trees object
T=trees.Trees(tree_list)
print "A set of trees with variable types", T.Types(), \
      " (initialized by ", tv, "):"
print T
print "Display Trees object with detail level 2:"
T.Display(Detail=2)
print "Display Trees object with data and detail level 1:"
T.Display(ViewPoint="DATA", Detail=1)
T.Plot()
print "Number of trees (supposedly ", nbtrees, "): ", T.NbTrees()
print "Comparison of trees contained in Trees " \
      "and those used for initialization: "
equal=True
differences=[]
for n in range(len(tree_list)):
    if str(T.Tree(n))!=str(tree_list[n]):
        equal=False
        differences.append(n)
if not equal:
    print "Differences found: ", differences
else:
    print "Match!"
    # print "Trees.Tree(", n, "):"
    # print T.Tree(n)
    # print "Tree number ", n, ":"
    # print tree_list[n]
    
# trees with one vertex
print "Test using trees of size one"
omax_size=1
onbtrees=2
# define a set of trees
otree_list=[]
otv=[1., 0, 1, 2.] # trees.TreeValue([1., 0])
OR=trees.TreeStructure(distrib, omax_size, max_depth)
otmp_tree=trees.Tree(otv, OR)
on=1
otree_list.append(trees.Tree(otmp_tree))
while on < onbtrees:
    on=on+1
    OR.Simulate(distrib, omax_size, max_depth)
    otmp_tree=trees.Tree(otv, OR)
    otree_list.append(trees.Tree(otmp_tree))    
for on in range(len(otree_list)):
    otree_list[on].Simulate(distrib_list)
# initialize a Trees object
OT=trees.Trees(otree_list)
print "A set of trees with variable types", T.Types(), \
      " (initialized by ", tv, "):"
print OT
# copy constructor
print "Test of the constructor by copy: "
T2=trees.Trees(T)
print T2
# constructor from a MTG
print "Read Trees object from MTG file 'sample_mtg_forest.txt':"
T=trees.Trees("sample_mtg_forest.txt")
print "Read ", T.NbTrees(), " trees."
print T
print "Print 1st tree of Trees object:"
print T.Tree(0).Display()
import amlPy
amlPy.setmode(1)
# reading a Tree from a MTG with filter and custom attributes
print "Read tree from MTG file 'sample_mtg_forest.txt' using a filter", \
      "and custom attributes:"
f=lambda x: amlPy.Feature(x, "Length")*3*(amlPy.Feature(x, "Diam")/2)**2
# FM=amlPy.MTG("sample_mtg_forest.txt")
# f=lambda x: amlPy.Feature(x, "Diam")
# print [amlPy.Feature(2, "Length"), amlPy.Feature(2, "Diam")]
# print f(5)
filter=lambda x: x < 6
# filter=lambda x: True
attributes=["anything"]
T=trees.Trees("sample_mtg_forest.txt", filter, attributes, [f], scale=2)
print "Read ", T.NbTrees(), " trees."
print T
print "Print 1st tree of Trees object:"
print T.Tree(0).Display()
# checking exceptions
# invalid file name
try:
    T=trees.Trees("no_such_file.t\/t")
except IOError, e:
    print e
else:
    print "Failed to raise exception for invalid MTG file name"
    raise e
# inconsistent attribute number
try:
    T=trees.Trees("sample_mtg_forest.txt", filter, [], [f], scale=2)
except ValueError, v:
    print v
else:
    print "Failed to raise exception for inconsistent attribute number"
# bad attribute name
try:
    T=trees.Trees("sample_mtg_forest.txt", filter, [f], [f], scale=2)
except TypeError, t:
    print t
else:
    print "Failed to raise exception for bad attribute name"
# bad attribute function
try:
    T=trees.Trees("sample_mtg_forest.txt", filter, attributes, attributes, 
                  scale=2)
except TypeError, t:
    print t
else:
    print "Failed to raise exception for bad attribute function"
# bad attribute type
try:
    T=trees.Trees("sample_mtg_forest.txt", filter, attributes, 
    [lambda x: amlPy.Descendants(x)], scale=2)
except TypeError, t:
    print t
else:
    print "Failed to raise exception for bad attribute type"
# bad attribute name
try:
    T=trees.Trees("sample_mtg_forest.txt", filter, [f], [f], scale=2)
except TypeError, t:
    print t
else:
    print "Failed to raise exception for bad attribute name"
# bad filtering function
try:
    T=trees.Trees("sample_mtg_forest.txt", attributes, attributes, [f], scale=2)
except TypeError, t:
    print t
else:
    print "Failed to raise exception for bad attribute function"
# filtering function does not filter descendants
try:
    T=trees.Trees("sample_mtg_forest.txt", lambda x: x != 2, attributes, [f], 
                  scale=2)
except IndexError, i:
    print i
else:
    print "Failed to raise exception for filter not filtering descendants"
try:
    T=trees.Trees("sample_mtg_forest.txt", filter, attributes, [f], scale=1)
except TypeError, t:
    print t

# T.Plot(ViewPoint="FirstOccurrenceRoot")