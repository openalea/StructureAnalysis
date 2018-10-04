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

dist0 = Distribution("NEGATIVE_BINOMIAL", 0, 1, 0.3)
dist1 = Distribution("data/distribution1.dist")

dist2 = Distribution("B", 0, 10, 0.3)
dist3 = Distribution("NB", 0, 3.5, 0.3)

#dist3.plot()

