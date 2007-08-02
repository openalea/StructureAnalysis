# a test for the class hmt.Hmt: estimation - syntax
import sys, os
AMLLIBDIR=os.getenv("HOME")+"/devlp/AMAPmod/build-linux/lib"
sys.path.append(AMLLIBDIR)
map(lambda x:sys.path.append(os.getenv("AMAPDIR")+"/STAT_TREES/Python/"+x),
    ["Int_fl_containers/", "Trees/", "Hmt/"])
sys.path.append(os.getenv("AMAPDIR")+"/STAT_TOOL/Python/")
import stat_tools, trees, hmt
inf_bound=1
sup_bound=3
probability= 0.6
ident=stat_tools.DistributionIdentifier.UNIFORM
parameter=stat_tools.D_DEFAULT
distrib= stat_tools.Parametric(ident, inf_bound, sup_bound, parameter, probability)
# Distribution used for the number of children and the tree attributes
file_name="hmot_np_2s.hmt"
# read a HMT from a file
print 'A hidden Markov tree H read from file "', file_name, '"'
H=hmt.HiddenMarkovTree(file_name)
H2O=hmt.HiddenMarkovTree("hmot_np_2o.hmt")
sample_size=2
tree_size=30
nb_children=2
H.Display()
# simulate labels from this HMT
print "Label simulation using H (1st tree): "
T=H.Simulate(sample_size, tree_size, nb_children)
print T.Tree(0)
# delete the state variable
print "Delete the state variable and the smoothed probabilities."
T=T.SelectVariable(1, "Keep")
print "State variable deleted: "
print T.Tree(0)
# parameter estimation from the simulated trees
# initialization from a model
print "Parameter estimation from the simulated trees" + \
" using the initial model specification:"
EH=T.Estimate("HIDDEN_MARKOV_TREE", H, 20)
EH.Display()
# checking exceptions raised by Estimate
print "Parameter estimation using a initial HMT "+ \
"with bad number of output processes:"
try:
    EH=T.Estimate("HIDDEN_MARKOV_TREE", H2O, 20)
except stat_tools.FormatError, f:
    print f
else:
    print "Failed to raise exception for bad number of output processes"
EH.Display()
# parameter estimation from the simulated trees
# initialization from a self-transition probability
print "Parameter estimation from the simulated trees" + \
" using the self-transition probability specification:"
EH=T.Estimate("HIDDEN_MARKOV_TREE", 2, "Irreductible", 0.999, 20)
EH.Display()
# Extract the trees from model
print "Extract the data part of the initial (file) HMT"
try:
    HMTD=H.ExtractData()
except stat_tools.FormatError, f:
    print f
else:
    print "Failed to raise exception for no data part in HiddenMarkovTree"
    print HMTD
print "Extract the data part of the estimated HMT"
HMTD=EH.ExtractData()
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
S=T.ComputeStateTrees(H, "Viterbi")
# Copy the state trees
CPS=hmt.HiddenMarkovTreeData(S)
S.Display(Detail=1)
HISTO=S.ExtractHistogram("Value", 1)
print "Marginal distribution for variable 1:"
print HISTO
print "Observation histogram for variable 1 and state 1:"
OBS1=S.ExtractHistogram("Observation", 1, 1)
print OBS1
# State profiles
EH.Display(ViewPoint="StateProfile", TreeId=1, NbStateTrees=4)
# parameter estimation from the simulated trees
# initialization from a self-transition probability
# forcing parametric estimation
print "Parameter estimation from the simulated trees" + \
" using the self-transition probability specification" + \
" forcing parametric estimation:"
EH=T.Estimate("HIDDEN_MARKOV_TREE", 2, "Irreductible", 0.999, 20, 
              ForceParametric=[True])
EH.Display()
# Re-estimate on segmented tree
print "Estimate an HMT on state tree: "
SS=S.SelectVariable([0])
EH2=SS.Estimate("HIDDEN_MARKOV_TREE", 2, "Irreductible", 0.999, 20)
EH2.Display()
