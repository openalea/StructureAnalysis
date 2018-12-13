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

meri1 = FrequencyDistribution("data/meri1.his")
meri2 = FrequencyDistribution("data/meri2.his")
meri3 = FrequencyDistribution("data/meri3.his")
meri4 = FrequencyDistribution("data/meri4.his")
meri5 = FrequencyDistribution("data/meri5.his")



#########################################################################
#
#  finite mixture of discrete distributions
#
#########################################################################

mixt1 = Mixture("data//mixture1.mixt")
mixt1 = Mixture(0.6, Distribution("B", 2, 18, 0.5),
                0.4, Distribution("NB", 10, 10, 0.5))

mixt_histo1 = mixt1.simulation(200)

