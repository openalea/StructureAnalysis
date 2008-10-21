
#########################################################################
#########################################################################
#
#  Frequency distributions
#
#  Objective: Analyzing the number of nodes of growth units in selected architectural
#             position considering the respective roles of preformation and neoformation,
#
#  Methods: comparison tests, one-way variance analysis,
#           estimation of finite mixture of distributions.
#
#  Wild cherry tree: number of nodes per growth unit (GU)
#
#  Data: Dominique Fournier
#
#  meri1.his: order 1,
#  meri1.his: order 2,
#  meri1.his: order 3, GU 1,
#  meri1.his: order 3, GU 2,
#  meri5.his: short shoots.
#
#
#  Poplar: number of nodes per growth unit
#
#  Data: Yves Caraglio and Herve Rey
#
#  peup1.his: order 2,
#  peup2.his: order 3,
#  peup3.his: order 4,
#  peup4.his: order 5,
#  peup5.his: order 3, GU 4,
#  peup6.his: order 3, acrotony.
#
#########################################################################

#EchoOn()
from openalea.stat_tool import *

plot.DISABLE_PLOT = False
 
meri1 = Histogram("meri1.his")
meri2 = Histogram("meri2.his")
meri3 = Histogram("meri3.his")
meri4 = Histogram("meri4.his")
meri5 = Histogram("meri5.his")

Plot(meri1, meri2, meri3, meri4, meri5)
Compare(meri1, meri2, meri3, meri4, meri5, "N")


ComparisonTest("F", meri1, meri2)
ComparisonTest("T", meri1, meri2)
ComparisonTest("W", meri1, meri2)

ComparisonTest("F", meri1, meri3)
ComparisonTest("T", meri1, meri3)
ComparisonTest("W", meri1, meri3)

# estimation of a mixture of two distributions assuming a first sub-population of GUs
# made only of a preformed part and a second sub-population made of both a preformed part
# and a neoformed part

mixt1 = Estimate(meri2, "MIXTURE", "B", "B")

meri = Merge(meri1, meri2, meri3, meri4, meri5)

# model selection approach: estimation of both the mixture parameters and
# the number of components 

mixt2 = Estimate(meri, "MIXTURE", "B", "B", "B", "B",  NbComponent="Estimated")
# mixt2 = Estimate(meri, "MIXTURE", "NB", "NB")
# Plot(ExtractDistribution(mixt2, "Mixture"))
# Display(mixt2)

mixt_data = ExtractData(mixt2)


dist5 = Estimate(meri5, "BINOMIAL")
# Display(dist5, Detail->2)
# Plot(dist5)

histo5 = Simulate(dist5, 100)
# Display(histo5, Detail->2)
# Plot(histo5)


peup1 = Histogram("peup1.his")
peup2 = Histogram("peup2.his")
peup3 = Histogram("peup3.his")
peup4 = Histogram("peup4.his")
peup5 = Histogram("peup5.his")
peup6 = Histogram("peup6.his")

mixt10 = Estimate(peup2, "MIXTURE", "B", "NB", "NB", "NB", NbComponent="Estimated")

peup = Merge(peup1, peup2, peup3, peup4, peup5, peup6)

histo1 = Shift(peup, -1)
histo2 = Cluster(peup, "Information", 0.8)
histo3 = Cluster(peup, "Step", 10)
histo4 = Cluster(peup, "Limit", [13, 24])
# Display(histo4, Detail->2)
# Plot(histo4)


mixt11 = Estimate(peup, "MIXTURE", "B", "NB", "NB", "NB", NbComponent="Estimated")
# mixt11 = Estimate(peup, "MIXTURE", "B", "NB")

d11 = distribution.Binomial(2, 12, 0.1)
d12 = distribution.Binomial(0, 10, 0.5)
d13 = distribution.Binomial(3, 10, 0.8)

d21 = distribution.Poisson(2, 8.0)
d22 = distribution.Poisson(4, 5.0)
d23 = distribution.Poisson(0, 2.0)

m = mixture._MvMixture([0.1, 0.2, 0.7], [[d11, d21], [d12, d22], [d13, d23]])
print m

m2 = mixture._MvMixture("mixture_mv1.mixt")
print m2

print "Egalite des melanges construits par liste ",\
      "de distributions et par fichiers : ", str(str(m)==str(m2))

m = mixture._MvMixture("mixture_mv_nonparam.mixt")
print m

print "Simulation de melanges multivaries : "
v = m.simulate(300);
print v

# m.old_plot(Title="Simulated mixture")
# m.old_plot(variable=1, Title="")
