# a test for the class hmt.Hmt: estimation - syntax
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
file_name="hmot_np_2s.hmt"
# read a non parametric HMT from a file
print 'A hidden Markov tree H read from file "', file_name, '"'
H=hmt.HiddenMarkovTree(file_name)
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
OBS1.Display(Detail=2)
# OBS1.Plot()
file_name="hmot_param.hmt"
# read a non parametric HMT from a file
print 'A hidden Markov tree H read from file "', file_name, '"'
H=hmt.HiddenMarkovTree(file_name)
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
OBS2.Display()
# OBS2.Plot()