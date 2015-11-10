"""
Simple test to illustrate usage of GNUPlot on Mixture data
"""
# import stat_tool modules
from openalea.stat_tool import *
from openalea.stat_tool.plot  import *

# Download histograms
meri1 = Histogram(get_shared_data("meri1.his"))
meri2 = Histogram(get_shared_data("meri2.his"))
meri3 = Histogram(get_shared_data("meri3.his"))
meri4 = Histogram(get_shared_data("meri4.his"))
meri5 = Histogram(get_shared_data("meri5.his"))

#Merge the histograms
meri = Merge(meri1, meri2, meri3, meri4, meri5)
Plot(meri)
# Estimate
mixt2 = Estimate(meri, "MIXTURE", "B", "B", "B", "B",  NbComponent="Estimated")

import matplotlib as mpl
mpl.rc('figure', figsize=(8,8))
Plot(mixt2, fontsize=8) 
