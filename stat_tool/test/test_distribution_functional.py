"""Mixture tests from exploratory.aml
#
#  Frequency distributions
#
#  Objective: Analyzing the number of nodes of growth units in
#    selected architectural position considering the respective roles of
#    preformation and neoformation,
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
"""
__version__ = "$Id: test_distribution_functional.py 9397 2010-08-10 10:50:07Z cokelaer $"


from openalea.stat_tool.plot import DISABLE_PLOT


from openalea.stat_tool.data_transform import ExtractDistribution
from openalea.stat_tool.histogram import Histogram
from openalea.stat_tool.estimate import Estimate
from openalea.stat_tool.comparison import Compare, ComparisonTest

from openalea.stat_tool.output import plot, Plot, Display
from openalea.stat_tool.data_transform import Merge, Shift, ExtractData
from openalea.stat_tool.simulate import Simulate
from openalea.stat_tool.cluster import Cluster

from openalea.stat_tool import get_shared_data

plot.DISABLE_PLOT = DISABLE_PLOT


def test():
    meri1 = Histogram(get_shared_data("meri1.his"))
    meri2 = Histogram(get_shared_data("meri2.his"))
    meri3 = Histogram(get_shared_data("meri3.his"))
    meri4 = Histogram(get_shared_data("meri4.his"))
    meri5 = Histogram(get_shared_data("meri5.his"))


    Plot(meri1, meri2, meri3, meri4, meri5)
    Compare(meri1, meri2, meri3, meri4, meri5, "N")


    ComparisonTest("F", meri1, meri2)
    ComparisonTest("T", meri1, meri2)
    ComparisonTest("W", meri1, meri2)

    ComparisonTest("F", meri1, meri3)
    ComparisonTest("T", meri1, meri3)
    ComparisonTest("W", meri1, meri3)

    # Estimation of a mixture of two distributions assuming a first
    # sub-population of GUs made only of a preformed part and a second
    # sub-population made of both a preformed part and a neoformed part


    _mixt1 = Estimate(meri2, "MIXTURE", "B", "B")

    meri = Merge(meri1, meri2, meri3, meri4, meri5)

    #model selection approach: estimation of both the mixture parameters and
    # the number of components

    mixt2 = Estimate(meri, "MIXTURE", "B", "B", "B", "B",
                     NbComponent="Estimated")
    mixt2 = Estimate(meri, "MIXTURE", "NB", "NB")
    Plot(mixt2)
    Plot(ExtractDistribution(mixt2, "Mixture"))

    print type(ExtractDistribution(mixt2, "Component", 1))


    Plot(ExtractDistribution(mixt2, "Component", 1),
         ExtractDistribution(mixt2, "Component", 2))
    Display(mixt2)

    _mixt_data = ExtractData(mixt2)

    dist5 = Estimate(meri5, "BINOMIAL")
    Display(dist5, Detail=2)
    Plot(dist5)

    histo5 = Simulate(dist5, 100)
    Display(histo5, Detail=2)
    Plot(histo5)


    peup1 = Histogram(get_shared_data( "peup1.his"))
    peup2 = Histogram(get_shared_data( "peup2.his"))
    peup3 = Histogram(get_shared_data( "peup3.his"))
    peup4 = Histogram(get_shared_data( "peup4.his"))
    peup5 = Histogram(get_shared_data( "peup5.his"))
    peup6 = Histogram(get_shared_data( "peup6.his"))

    _mixt10 = Estimate(peup2, "MIXTURE", "B", "NB", "NB", "NB",
                       NbComponent="Estimated")

    peup = Merge(peup1, peup2, peup3, peup4, peup5, peup6)

    _histo1 = Shift(peup, -1)
    _histo2 = Cluster(peup, "Information", 0.8)
    _histo3 = Cluster(peup, "Step", 10)
    histo4 = Cluster(peup, "Limit", [13, 24])
    Display(histo4, Detail=2)
    Plot(histo4)


    _mixt11 = Estimate(peup, "MIXTURE", "B", "NB", "NB", "NB",
                       NbComponent="Estimated")
    _mixt11 = Estimate(peup, "MIXTURE", "B", "NB")



if __name__ == "__main__":
    test()
