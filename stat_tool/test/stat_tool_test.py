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
from structure_analysis.stat_tool.frequency_distribution import Histogram
from structure_analysis.stat_tool.discrete_mixture import Mixture

dist0 = Distribution("NEGATIVE_BINOMIAL", 0, 1, 0.3)
dist1 = Distribution("data/distribution1.dist")

dist2 = Distribution("B", 0, 10, 0.3)
dist3 = Distribution("NB", 0, 3.5, 0.3)

#dist3.plot()

histo1 = dist1.simulation(200)

fagus = Histogram("data/fagus1.his")

meri1 = Histogram("data/meri1.his")
meri2 = Histogram("data/meri2.his")
meri3 = Histogram("data/meri3.his")
meri4 = Histogram("data/meri4.his")
meri5 = Histogram("data/meri5.his")



#########################################################################
#
#  finite mixture of discrete distributions
#
#########################################################################

mixt1 = Mixture("data//mixture1.mixt")
mixt1 = Mixture(0.6, Distribution("B", 2, 18, 0.5),
                0.4, Distribution("NB", 10, 10, 0.5))

mixt_histo1 = mixt1.simulation(200)

