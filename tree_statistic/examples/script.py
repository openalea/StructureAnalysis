# -*- coding: utf-8 -*-
# create a MarkovOutTreeData with types
# {Short / Long} x {Vegetative / Floral}
import sys, os, datetime
import openalea.tree_statistic.trees as trees
import openalea.tree_statistic.trees as trees, openalea.tree_statistic.mot as mot
import openalea.stat_tool as stat_tool
import openalea.tree_statistic.dgdistributions as dgd
import matplotlib.pyplot as plt

#from stat_tool.plot import gnuplot, mplotlib
#stat_tool.plot.PLOTTER = mplotlib()
#stat_tool.plot.PLOTTER = gnuplot()

import openalea.aml as amlPy
from amlPy import *

T = trees.Trees("Data/mtg_fuji_42_47.txt")
# 2 pommiers a l'echelle de la pousse annuelle
# variable 0 : nombre d'entre-noeus
# variable 1 : etat (0: long veg. 1: long fl. 2: court veg. 3: court fl.
# variable 2 : nombre de cycles de la pousse annuelle
# variable 3 : descendance censuree (cas en 1999)

# recupere le dictionnaire des sommets censures
censored = []
for tid in range(T.NbTrees()):
    censored += [{}]
    TA = T.Tree(tid)
    # identificateurs 
    vids = TA.TreeVertex().values()
    for v in vids:
        # variables observees
        val = TA.Get(v)
        censored[tid][v] = (val[3] > 0)

TASStates = T.SelectVariable([1])
TAS = mot.MarkovOutTreeData(T, StateVariable=1)
TAS.SetCensoredDict(censored)
#TAS.Plot(ViewPoint="StateProcess")
#TAS.Plot(ViewPoint="StateProcess")
#GP = TAS.ExtractGenerationProcess()
# TASStates = mot.MarkovOutTreeData(TASStates, StateVariable=0)

# Tailles d'arbres : 2197  2317

#EMBP = TAS.Estimate("MARKOV_OUT_TREE", 0, 0, 0, [], True, [], 0, True)

# Configuration de descendants etat parent 0

TAS._UpdateMarkovReestimation(0, 0, 0, [], ParentDependent=True)


V0 = TAS.ExtractVectorsGenerationProcess([0])

Data0 = dgd.DiscreteGraphicalData(4,V0)
# Print multivariate histogram
Data0.Display()
# Parametric graphical model estimated (using conditional entropies)
# with specified number of edges
M = Data0.ParametricEstimation(dgd.Algorithms.ML, number_of_edges = 0)
# Draw graph
M.Graph()
plt.show()

# Estimation of a sequence of probabilistic graphical models
# with increasing number of edges
MLs0 = Data0.DistributionEstimation(dgd.Algorithms.ML, all = True)
Data0.LogLikelihood(MLs0)
plt.show() # Draw log-likelihood as a function of the number of edges
Data0.PenalizedLogLikelihood(MLs0, dgd.Criterions.BIC)
plt.show()

# Best model with 6 edges
# M6 = Data0.ParametricEstimation(dgd.Algorithms.ML, number_of_edges = 6)


V1 = TAS.ExtractVectorsGenerationProcess([1])
Data1 = dgd.DiscreteGraphicalData(4,V1)
MLs1 = Data1.DistributionEstimation(dgd.Algorithms.ML, all = True)
Data1.LogLikelihood(MLs1)
plt.show()
Data1.PenalizedLogLikelihood(MLs1, dgd.Criterions.BIC)
plt.show()

V2 = TAS.ExtractVectorsGenerationProcess([2])
Data2 = dgd.DiscreteGraphicalData(4,V2)
MLs2 = Data2.DistributionEstimation(dgd.Algorithms.ML, all = True)
Data2.LogLikelihood(MLs2)
plt.show()
Data2.PenalizedLogLikelihood(MLs2, dgd.Criterions.BIC)
plt.show()

V3 = TAS.ExtractVectorsGenerationProcess([3])
Data3 = dgd.DiscreteGraphicalData(4,V3)
MLs3 = Data3.DistributionEstimation(dgd.Algorithms.ML, all = True)
Data3.LogLikelihood(MLs3)
plt.show()
Data3.PenalizedLogLikelihood(MLs3, dgd.Criterions.BIC)
plt.show()

# M3 = Data0.ParametricEstimation(dgd.Algorithms.ML, number_of_edges = 3)

