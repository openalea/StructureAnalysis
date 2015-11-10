"""
#########################################################################
#
#  discrete distributions/histograms, comparison of histograms/frequency
#  distributions
#
#  beech, Wild cherry tree: number of nodes per growth unit (GU)
#
#  meri1.his: order 1,
#  meri1.his: order 2,
#  meri1.his: order 3, GU 1,
#  meri1.his: order 3, GU 2,
#  meri5.his: short shoots.
#
#########################################################################
"""
from openalea.stat_tool.cluster import Cluster, Clustering
from openalea.stat_tool import VectorDistance, ComparisonTest
from openalea.stat_tool import Distribution, ExtractDistribution
from openalea.stat_tool import SelectIndividual
from openalea.stat_tool import Simulate, ExtractHistogram, Vectors
from openalea.stat_tool import Estimate, Plot, Fit, Merge
from openalea.stat_tool import Histogram, SelectVariable, ValueSelect
from openalea.stat_tool import Transcode, Compare, Regression, Shift
from openalea.stat_tool import ExtractData, Mixture, ToHistogram 
from openalea.stat_tool import Compound, Convolution, VarianceAnalysis
from openalea.stat_tool import ContingencyTable

from openalea.sequence_analysis import Sequences

from openalea.stat_tool import get_shared_data


dist0 = Distribution("NEGATIVE_BINOMIAL", 0, 1, 0.3)
dist0 = Distribution("data//distribution1.dist")

dist1 = Distribution("B", 0, 10, 0.3)
dist1 = Distribution("NB", 0, 3.5, 0.3)

# computation of survivior function and hasard rate

#Display(dist1, ViewPoint="Survival")
Plot(dist1, ViewPoint="Survival")

# simulation/estimation

histo1 = Simulate(dist1, 200)

#Display(histo1, ViewPoint="Survival")
Plot(histo1, ViewPoint="Survival")

# "B" "BINOMIAL", "P" "POISSON", "NB" "NEGATIVE_BINOMIAL"
# InfBoundStatus="Free" (default) / "Fixed"
# MinInfBound=0 (default) / 1

dist2 = Estimate(histo1, "NEGATIVE_BINOMIAL", MinInfBound=0,
                 InfBoundStatus="Fixed")

fagus = Histogram("data//fagus1.his")

# transformation of histograms, extraction/filter

histo2 = Cluster(fagus, "Step", 2)
histo3 = Cluster(fagus, "Information", 0.8)
histo4 = Cluster(fagus, "Limit", [2, 4, 6, 8, 10])
histo5 = Transcode(fagus, [1, 2, 2, 3, 3, 4, 4, 5])
#Display(histo5, Detail=2)

histo7 = Shift(fagus, -2)

histo8 = ValueSelect(fagus, 2, 8)

dist3 = Estimate(fagus, "B")

# comparison of histograms
meri1 = Histogram(get_shared_data( "meri1.his"))
meri2 = Histogram(get_shared_data( "meri2.his"))
meri3 = Histogram(get_shared_data( "meri3.his"))
meri4 = Histogram(get_shared_data( "meri4.his"))
meri5 = Histogram(get_shared_data( "meri5.his"))


# Compare(meri1, meri2, meri3, meri4, meri5, "N", FileName="ASCII/meri.cmp")
Compare(meri1, meri2, meri3, meri4, meri5, "O")

ComparisonTest("F", meri1, meri2)
ComparisonTest("T", meri1, meri2)
ComparisonTest("W", meri1, meri2)

# fit of a known distribution to an  histogram

dist5 = Fit(meri5, Distribution("B", 0, 10, 0.437879))
#Display(dist5, Detail=2)
Plot(dist5)


#########################################################################
#
#  finite mixture of discrete distributions
#
#########################################################################

mixt1 = Mixture("data//mixture1.mixt")
mixt1 = Mixture(0.6, Distribution("B", 2, 18, 0.5), 0.4,
                Distribution("NB", 10, 10, 0.5))

mixt_histo1 = Simulate(mixt1, 200)

# extraction of histograms/frequency distributions corresponding to 
# a given mixture component
# (i.e. elementary distributions which are combined by mixture)

histo10 = ExtractHistogram(mixt_histo1, "Component", 1)
histo11 = ExtractHistogram(mixt_histo1, "Component", 2)
histo12 = Merge(histo10, histo11)
histo13 = ExtractHistogram(mixt_histo1, "Weight")

# estimation

mixt2 = Estimate(mixt_histo1, "MIXTURE", "B", "NB", MinInfBound=0,
                 InfBoundStatus="Fixed", DistInfBoundStatus="Fixed")

mixt_histo2 = ExtractData(mixt2)

histo14 = ExtractHistogram(ExtractData(mixt2), "Component", 1)
histo15 = ToHistogram(ExtractDistribution(mixt2, "Component", 1))

# estimation and selection of the number of components

mixt3 = Estimate(meri1, "MIXTURE", Distribution("B", 6, 7, 0.5), "B")

# Penalty="AIC"/ "AICc" / "BIC" / "BICc" (default),  NbComponent="Estimated"
 
meri = Merge(meri1, meri2, meri3, meri4, meri5)

mixt2 = Estimate(meri, "MIXTURE", "B", "B", "B", "B", 
                 NbComponent="Estimated", Penalty="BIC")
#Display(mixt2, Detail=2)
dist_mixt = ExtractDistribution(mixt2, "Mixture")
Plot(dist_mixt)


#########################################################################
#
#  convolution of discrete distributions
#
#########################################################################

convol1 = Convolution(Distribution("B", 0, 10, 0.5),
                      Distribution("NB", 0, 10, 0.5))
convol1 = Convolution("data//convolution1.conv")

convol_histo1 = Simulate(convol1, 200)

# is it the correct syntax ???
histo20 = ExtractHistogram(convol_histo1, 1)

convol2 = Estimate(convol_histo1, "CONVOLUTION",
                   ExtractDistribution(convol1, 1), MinInfBound=0)

histo21 = ExtractHistogram(ExtractData(convol2), 1)
histo22 = ToHistogram(ExtractDistribution(convol2, 1))

histo_b2 = Histogram("data//nothofagus_antarctica_bud_2.his")
histo_s2 = Histogram("data//nothofagus_antarctica_shoot_2.his")

# Estimator="Likelihood" (default) / "PenalizedLikelihood" / "Parametric"
# Si Estimator="PenalizedLikelihood", options supplementaires possibles
# Penalty="FirstDifference" / "SecondDifference" (default) / "Entropy", Weight,
# Outside="Zero" (default) / "Continuation" (cf. stat_funs4.cpp).

# "NP" or "NON_PARAMETRIC" (for estimation only)

convol30 = Estimate(Shift(histo_s2, 1), "CONVOLUTION",
                    Estimate(histo_b2, "NP"), NbIteration=500)
convol31 = Estimate(Shift(histo_s2, 1), "CONVOLUTION",
                    Estimate(histo_b2, "NP"), NbIteration=100,
                    Estimator="PenalizedLikelihood", Weight=0.5)
convol32 = Estimate(Shift(histo_s2, 1), "CONVOLUTION",
                    Estimate(histo_b2, "NP"), Estimator="Parametric")
#Display(convol31)
Plot(convol31)
Plot(ExtractDistribution(convol31, "Convolution"))
#Save(convol31, "Spreadsheet/nothofagus_antartica_2.xld", Format="SpreadSheet")


#########################################################################
#
#  compound distribution : not publicly distributed in the Python version
#
#########################################################################

cdist1 = Compound("data//compound1.cd")

chisto1 = Simulate(cdist1, 200)

histo30 = ExtractHistogram(chisto1, "Sum")

cdist2 = Estimate(chisto1, "COMPOUND", 
                  ExtractDistribution(cdist1, "Elementary"), "Sum",
                  MinInfBound=0)

histo31 = ExtractHistogram(ExtractData(cdist2), "Sum")
histo32 = ToHistogram(ExtractDistribution(cdist2, "Sum"))
    
peup1 = Histogram(get_shared_data( "peup1.his"))
mixt4 = Estimate(peup1, "MIXTURE", "B", "NB")
histo33 = ToHistogram(ExtractDistribution(mixt4, "Component", 2))
histo34 = Shift(histo33, -11)

cdist3 = Estimate(histo34, "COMPOUND", Distribution("B", 0, 1, 0.7), "Sum")
cdist3 = Estimate(histo34, "COMPOUND", Distribution("B", 0, 1, 0.7), "Sum",
                  MinInfBound=0)


#########################################################################
#
#  vectors : contingency table, one-way variance analysis,
#  linear regression or nonparametric regression (loess smoother)
#
#  Oak trunk annual shoots
#
#  INDEX_PARAMETER : TIME  (year of growth - 1995, 1996, 1997)
#  VARIABLE 1 : length of the annual shoot (cm)
#  VARIABLE 2 : diameter of the annual shoot (1/10 de mm)
#  VARIABLE 3 : number of cycles
#  VARIABLE 4 : number of nodes
#  VARIABLE 5 : number de branches
#
#########################################################################


seq0 = Sequences("data/chene_sessile_15pa.seq")
#Plot(seq0, ViewPoint="Data")

# change of unit for the variable diameter of the annual shoot

marginal2 = ExtractHistogram(seq0, "Value", 2)
print marginal2
Plot(Cluster(marginal2, "Information", 0.75))
Plot(Cluster(marginal2, "Information", 0.61))
Plot(Cluster(marginal2, "Step", 10))

vec10 = Vectors(seq0)

# plot of the pointwise averages
# Plot(Regression(vec10, "MovingAverage", 1, 2, [1]))

vec95 = ValueSelect(vec10, 1, 1995)
vec96 = ValueSelect(vec10, 1, 1996)
vec97 = ValueSelect(vec10, 1, 1997)

VarianceAnalysis(vec10, 1, 2, "N")
Compare(ExtractHistogram(vec95, 2), ExtractHistogram(vec96, 2),
        ExtractHistogram(vec97, 2), "N")
# Plot(ExtractHistogram(vec95, 2), ExtractHistogram(vec96, 2),
#     ExtractHistogram(vec97, 2))

ContingencyTable(vec10, 1, 4)

# one-way variance analysis based on ranks

VarianceAnalysis(vec10, 1, 4, "O")
Compare(ExtractHistogram(vec95, 4), ExtractHistogram(vec96, 4),
        ExtractHistogram(vec97, 4), "O")
# Plot(ExtractHistogram(vec95, 4), ExtractHistogram(vec96, 4),
# ExtractHistogram(vec97, 4))

# Plot(ExtractHistogram(vec95, 5), ExtractHistogram(vec96, 5),
#     ExtractHistogram(vec97, 5))
# Plot(ExtractHistogram(vec95, 6), ExtractHistogram(vec96, 6),
#     ExtractHistogram(vec97, 6))

vec11 = ValueSelect(vec10, 4, 1)
vec12 = ValueSelect(vec10, 4, 2)
vec13 = ValueSelect(vec10, 4, 3, 4)

#Plot(ExtractHistogram(vec11, 2), ExtractHistogram(vec12, 2),
#     ExtractHistogram(vec13, 2))
#Plot(ExtractHistogram(vec11, 5), ExtractHistogram(vec12, 5),
#     ExtractHistogram(vec13, 5))

mixt20 = Estimate(ExtractHistogram(vec10, 2),
                  "MIXTURE", "NB", "NB", "NB", "NB", NbComponent="Estimated")
# Display(mixt20)
Plot(mixt20)
Plot(ExtractDistribution(mixt20, "Mixture"))

mixt21 = Estimate(ExtractHistogram(vec10, 5),
                  "MIXTURE", "NB", "NB", "NB", "NB", NbComponent="Estimated")

vec9596 = ValueSelect(vec10, 1, 1995, 1996)
Plot(ExtractHistogram(ValueSelect(vec9596, 4, 1), 6),
     ExtractHistogram(ValueSelect(vec9596, 4, 2), 6),
     ExtractHistogram(ValueSelect(vec9596, 4, 3, 4), 6))

# linear regression

regress10 = Regression(vec10, "Linear", 5, 2)
# Display(regress10)
Plot(regress10)

# nonparametric regression (loess smoother)

regress11 = Regression(vec10, "NearestNeighbors",  5, 2, 0.3)

regress12 = Regression(vec9596, "Linear", 5, 6)
regress13 = Regression(vec9596, "NearestNeighbors", 5, 6, 0.5)

vec15 = SelectVariable(vec10, [1, 3, 6], Mode="Reject")


#########################################################################
#
#  Distance matrix and clustering (partitioning or hierarchical methods) 
#
#########################################################################

# computation of a distance matrix using a standardization procedure

matrix10 = Compare(vec15, VectorDistance("N", "N", "N"))

# clustering using a partitioning method

# Display(Clustering(matrix10, "Partition", 2))

vec151 = SelectIndividual(vec10,  
                          [69, 48, 41, 44, 32, 47, 81, 95, 11, 36, 75, 108, 56,
                           83, 38, 98, 113, 134, 110, 101, 77, 35, 74, 80, 50,
                           24, 89, 128, 5, 45, 8, 116, 119, 132, 61, 78, 53,
                           29, 131, 65, 90, 96, 104, 20, 86, 66, 42, 68,
                           125, 14, 23, 54, 33, 26, 71, 129, 102, 51, 70,
                           111, 138, 19, 127, 62, 117, 137, 2, 28, 17])

vec152 = SelectIndividual(vec10, [100, 13, 133, 105, 72, 9, 93, 109, 30, 115,
                                  63, 7, 55, 37, 15, 114, 106, 46, 73, 18, 3,
                                  87, 58, 43, 60, 76, 52, 6, 39, 31, 12, 99,
                                  121, 123, 22, 79, 94, 88, 21, 97, 25, 40,
                                  57, 136, 67, 49, 10, 4, 120, 92, 27, 91,
                                  64, 124, 16, 130, 84, 107, 126, 103, 122, 
                                  112, 59, 1, 82, 34, 135, 118, 85])
Plot(ExtractHistogram(vec151, 4), ExtractHistogram(vec152, 4))

matrix11 = Compare(vec15, VectorDistance("N", "O", "N"))

Clustering(matrix10, "Hierarchy", Algorithm="Agglomerative")
Clustering(matrix10, "Hierarchy", Algorithm="Divisive")

vec16 = SelectVariable(vec9596, [1, 3], Mode="Reject")
matrix12 = Compare(vec16, VectorDistance("N", "N", "N", "N"))
matrix13 = Compare(vec16, VectorDistance("N", "O", "N", "N"))


