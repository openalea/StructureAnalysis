import stat_tool
import os
%reload_ext stat_tool.mplotlib
%reload_ext stat_tool.aml

meri1 = stat_tool.Histogram(str(os.path.join("..", "share", "data", "meri1.his")))
meri2 = stat_tool.Histogram(str(os.path.join("..", "share", "data", "meri2.his")))
meri3 = stat_tool.Histogram(str(os.path.join("..", "share", "data", "meri3.his")))
meri4 = stat_tool.Histogram(str(os.path.join("..", "share", "data", "meri4.his")))
meri5 = stat_tool.Histogram(str(os.path.join("..", "share", "data", "meri5.his")))
meri0 = stat_tool.Merge(meri1, meri2, meri3, meri4, meri5)

mixt = stat_tool.MixtureEstimation(meri0, 4, "BINOMIAL", display=True)

from matplotlib import pyplot
fig = mixt.plot()
pyplot.show()
