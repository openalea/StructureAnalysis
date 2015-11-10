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
meri_a = Merge(meri1, meri2, meri3, meri4, meri5)
print meri_a.display()

meri_b = meri1.merge([meri2, meri3, meri4, meri5])
print meri_b.display()


