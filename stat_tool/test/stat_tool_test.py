#########################################################################
#
#  discrete distributions, comparison of frequency distributions
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

import structure_analysis.stat_tool as st
from structure_analysis.stat_tool.distribution import Distribution
from structure_analysis.stat_tool.frequency_distribution import FrequencyDistribution
from structure_analysis.stat_tool.discrete_mixture import Mixture
from structure_analysis.stat_tool.estimate import Estimate

dist0 = Distribution("NEGATIVE_BINOMIAL", 0, 1, 0.3)
# print(dist0)
# print(dist0.ascii_write(True)) # Display(dist0, Detail->2)

dist1 = Distribution("data/distribution1.dist")

dist2 = Distribution("B", 0, 10, 0.3)
dist3 = Distribution("NB", 0, 3.5, 0.3)
# print (dist3.survival_ascii_write(False)) # Display(dist1, ViewPoint->"Survival")

#dist3.plot()

histo1 = dist1.simulation(200)
#histo1.plot()

dist2 = Estimate(histo1, "NEGATIVE_BINOMIAL", MinInfBound=0, InfBoundStatus="Fixed")


fagus = FrequencyDistribution("data/fagus1.his")

## transformation of frequency distributions, extraction/filter

histo2 = Cluster(fagus, "Step", 2)
histo3 = Cluster(fagus, "Information", 0.8)
histo4 = Cluster(fagus, "Limit", [2, 4, 6, 8, 10])
histo5 = Transcode(fagus, [1, 2, 2, 3, 3, 4, 4, 5])
# Display(histo5, Detail->2)
# Plot(fagus, histo2, histo3, histo4, histo5)

histo7 = Shift(fagus, -2)

histo8 = ValueSelect(fagus, 2, 8)

dist3 = Estimate(fagus, "B")
# Display(dist3, Detail->2)
# Plot(dist3)

## comparison of frequency distributions

meri1 = FrequencyDistribution("data/meri1.his")
meri2 = FrequencyDistribution("data/meri2.his")
meri3 = FrequencyDistribution("data/meri3.his")
meri4 = FrequencyDistribution("data/meri4.his")
meri5 = FrequencyDistribution("data/meri5.his")

# Compare(meri1, meri2, meri3, meri4, meri5, "Nu", FileName->"ASCII/meri.cmp")
Compare(meri1, meri2, meri3, meri4, meri5, "O")

ComparisonTest("F", meri1, meri2)
ComparisonTest("T", meri1, meri2)
ComparisonTest("W", meri1, meri2)

# fit of a known distribution to a frequency distribution

dist5 = Fit(meri5, Distribution("B", 0, 10, 0.437879))
# Display(dist5, Detail->2)
# Plot(dist5)


#########################################################################
#
#  finite mixture of discrete distributions
#
#########################################################################
DiscreteMixture = Mixture

mixt1 = DiscreteMixture("data//mixture1.mixt")
mixt1 = DiscreteMixture(0.6, Distribution("B", 2, 18, 0.5),
                0.4, Distribution("NB", 10, 10, 0.5))

mixt_histo1 = mixt1.simulation(200)

## extraction of frequency distributions corresponding to a given mixture component
## (i.e. elementary distributions which are combined by mixture)
#########################################################################

histo10 = ExtractFrequencyDistribution(mixt_histo1, "Component", 1)
histo11 = ExtractFrequencyDistribution(mixt_histo1, "Component", 2)
histo12 = Merge(histo10, histo11)
histo13 = ExtractFrequencyDistribution(mixt_histo1, "Weight")
# Plot(histo10, histo11, histo12)

## estimation
#########################################################################

mixt2 = Estimate(mixt_histo1, "DISCRETE_MIXTURE", "B", "NB", MinInfBound->0, InfBoundStatus->"Fixed", DistInfBoundStatus->"Fixed")
# Display(mixt2, Detail->2)
# Plot(mixt2)
mixt_histo2 = ExtractData(mixt2)

histo14 = ExtractFrequencyDistribution(ExtractData(mixt2), "Component", 1)
histo15 = ToFrequencyDistribution(ExtractDistribution(mixt2, "Component", 1))

## estimation with a known component
#########################################################################

mixt3 = Estimate(meri1, "DISCRETE_MIXTURE", Distribution("B", 6, 7, 0.5), "B")
# Display(mixt3, Detail->2)
# Plot(mixt3)

## estimation and selection of the number of components
#########################################################################

# NbComponent->"Fixed" (default) / "Estimated"
# Criterion->"AIC"/ "AICc" / "BIC" / "BICc" (default), valid option if NbComponent->"Estimated"

meri = Merge(meri1, meri2, meri3, meri4, meri5)

mixt4 = Estimate(meri, "DISCRETE_MIXTURE", "B", "B", "B", "B",
    NbComponent="Estimated", Criterion="BIC")
# Display(mixt4, Detail->2)
# Plot(mixt4)
mixt_histo4 = ExtractDistribution(mixt4, "Mixture")
# Plot(mixt_histo4)

#########################################################################
# TODO
#########################################################################

#########################################################################
#
#  convolution of discrete distributions
#
#########################################################################
