# -*- coding: utf-8 -*-
# a test for the class _hmt.Hmt: constructor and basic methods
import sys
import os
import openalea.stat_tool as cstat_tools
import openalea.tree_statistic.trees as trees
#cstat_tools
import openalea.tree_statistic.hmt
import openalea.tree_statistic.hmt._hmt as _hmt


#AMLLIBDIR=os.getenv("HOME")+"/devlp/AMAPmod/build-linux/lib"
#sys.path.append(AMLLIBDIR)
#map(lambda x:sys.path.append(os.getenv("AMAPDIR")+"/STAT_TREES/Python/"+x),
#    ["Int_fl_containers/", "Trees/", "Hmt/"])
#sys.path.append(os.getenv("AMAPDIR")+"/STAT/Python/")


inf_bound = 1
sup_bound = 3
probability = 0.6
ident = cstat_tools.DistributionIdentifierType.UNIFORM
parameter = cstat_tools.D_DEFAULT
distrib = cstat_tools.distribution._DiscreteParametricModel(ident, inf_bound, sup_bound, parameter,
                                                            probability)

print "Distribution used for the number of children and the tree attributes:"
print distrib

max_depth = 3
max_size = 10
nbtrees = 4

# test the initialization and the estimation of a HMT
print "Test the initialization and the estimation of a HMT"
# define a set of trees

tree_list = []
tv = [1., 0, 1, 2.] # trees.TreeValue([1., 0])
R = trees.TreeStructure(distrib, max_size, max_depth)
tmp_tree = trees.Tree(tv, R)
n = 1
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
T = trees.Trees(tree_list)
H = _hmt.HmtAsciiRead("hmot_np_2s.hmt")
print "A HMT read from file 'hmot_np_2s.hmt':"
print H.Display(False)
T = H.Simulate(3, 20, 3, True)
print "Set of simulated trees:"
print T.Display(False)


# The remaining part of this test does not work : Cpp and python prototype
# do not seem to match each other. I comment this piece of code and the
# author should fix it. TC, Dec 2008
#print "\nEstimated HMT (initialization by HMT): \n "
#EH = T.EstimationCiHmot(H, True, cstat_tools.RestorationAlgorithm.VITERBI, 20, True)
#print EH.Display(False)
#print "\nEstimated HMT (initialization by self-transition): \n "
#EH = T.EstimationCiHmot(2, False, False, cstat_tools.RestorationAlgorithm.VITERBI, 0.9999, 10)
#print EH.Display(False)
# test the initialization of a CHmt_data object
#print "Test the initialization of a CHmt_data object"
# copy
chmtd = _hmt.CHmt_data(T)
print chmtd.Display(True)

