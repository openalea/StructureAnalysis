"""
Simple test to illustrate usage of GNUPlot on Mixture data
"""
# import stat_tool modules
from openalea.stat_tool import *
from openalea.stat_tool.plot  import *

# Download histograms
meri1 = Histogram("../../test/data/meri1.his")
meri2 = Histogram("../../test/data/meri2.his")
meri3 = Histogram("../../test/data/meri3.his")
meri4 = Histogram("../../test/data/meri4.his")
meri5 = Histogram("../../test/data/meri5.his")

#Merge the histograms
meri = Merge(meri1, meri2, meri3, meri4, meri5)
# Estimate
mixt2 = Estimate(meri, "MIXTURE", "B", "B", "B", "B",  NbComponent="Estimated")

# plot using gnuplot 
plotter = gnuplot()
plotable = mixt2.get_plotable()
a = plotter.plot(plotable, "dummy")
# save plots in postscript format
mixt2.plot_print()

#matplotlib by default but not yet implemented
#Plot(mixt2)    


