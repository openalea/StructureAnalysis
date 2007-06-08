# a test for the class hmt.Hmt: estimation accurracy
import sys, os
AMLLIBDIR=os.getenv("HOME")+"/devlp/AMAPmod/build-linux/lib"
sys.path.append(AMLLIBDIR)
map(lambda x:sys.path.append(os.getenv("AMAPDIR")+"/STAT_TREES/Python/"+x),
    ["Int_fl_containers/", "Trees/", "Hmt/"])
sys.path.append(os.getenv("AMAPDIR")+"/STAT/Python/")
import stat_tools, trees, hmt
inf_bound=1
sup_bound=3
probability= 0.6
ident=stat_tools.DistributionIdentifier.UNIFORM
parameter=stat_tools.D_DEFAULT
distrib= stat_tools.Parametric(ident, inf_bound, sup_bound, parameter, probability)
# Distribution used for the number of children and the tree attributes
H=hmt.HiddenMarkovTree("hmot.hmt")
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
EH=T.Estimate("HIDDEN_MARKOV_TREE", 3, "LEFTRIGHT", 0.99, 100)
EH.Display()