""" extract from exploratoty.aml

#########################################################################
#
#  Multivariate samples
#
#  Objective: investigating the relations between different variables
#             measured on annual shoots,
#
#  Methods: contingency table, one-way variance analysis,
#           linear/nonparametric regression, estimation of
#           finite mixture of distributions, clustering methods.
#
#########################################################################
#
#  beech annual shoots
#
#  Data: Eric Nicolini
#
#  VARIABLE 1 : orientation (
                0: horizontal,
                1: oblique,
                2: vertical,
               3: pruned 3,
                4: absent)
#  VARIABLE 2 : number of cycles
#  VARIABLE 3 : apex death (
                1: monocyclic alive,
                2: monocyclic dead,
                3: polycyclic alive,
                4: polycyclic spring dead,
                5: polycyclic summer dead,
#               6: polycyclic spring-summer dead)
#  VARIABLE 4 : number of nodes
#  VARIABLE 5 : length of the annual shoot (cm)
#
#########################################################################
"""

from openalea.stat_tool import Vectors
from openalea.stat_tool import ExtractHistogram
from openalea.stat_tool import Compare
from openalea.stat_tool import Estimate
from openalea.stat_tool import Plot
from openalea.stat_tool import Display
from openalea.stat_tool import ExtractDistribution
from openalea.stat_tool import ValueSelect
from openalea.stat_tool import VarianceAnalysis
from openalea.stat_tool import Regression
from openalea.stat_tool import VectorDistance
from openalea.stat_tool import Transcode

try:
    from openalea.sequence_analysis import Sequences
except ImportError:
    Sequences = None

def test():
    if Sequences:
        pass
    else:
        return
    vec1 = Vectors(Sequences("../examples/Sample/Sequences/hetre.seq"))

    Plot(ExtractHistogram(ValueSelect(vec1, 1, 0), 4),
         ExtractHistogram(ValueSelect(vec1, 1, 1), 4),
         ExtractHistogram(ValueSelect(vec1, 1, 2), 4),
         ExtractHistogram(ValueSelect(vec1, 1, 3), 4),
         ExtractHistogram(ValueSelect(vec1, 1, 4), 4))

    VarianceAnalysis(vec1, 2, 4, "N")
    Plot(ExtractHistogram(ValueSelect(vec1, 2, 1), 4),
         ExtractHistogram(ValueSelect(vec1, 2, 2), 4))

    mixt10 = Estimate(ExtractHistogram(vec1, 4), "MIXTURE", "B", "B", "B", "B",
                      NbComponent="Estimated")

    Plot(Estimate(ExtractHistogram(ValueSelect(vec1, 2, 1), 4), "NP"),
         ExtractDistribution(mixt10, "Component", 1))
    Plot(Estimate(ExtractHistogram(ValueSelect(vec1, 2, 2, 3), 4), "NP"),
         ExtractDistribution(mixt10, "Component", 2))

    Plot(ExtractHistogram(ValueSelect(vec1, 2, 1), 5),
         ExtractHistogram(ValueSelect(vec1, 2, 2), 5))

    Plot(ExtractHistogram(ValueSelect(vec1, 3, 1), 4),
         ExtractHistogram(ValueSelect(vec1, 3, 2), 4))
    Plot(ExtractHistogram(ValueSelect(vec1, 3, 3), 4),
         ExtractHistogram(ValueSelect(vec1, 3, 4, 6), 4))

    regress1 = Regression(vec1, "Linear", 4, 5)
    Display(regress1)
    Plot(regress1)

    regress2 = Regression(vec1, "NearestNeighbors", 4, 5, 0.1)
    Plot(regress2)
    regress2 = Regression(vec1, "NearestNeighbors", 4, 5, 0.2)

    matrix1 = Compare(Transcode(vec1, 3, [1, 0, 1, 0, 0, 0]),
                      VectorDistance("S", "O", "S", "N", "N"))


if __name__ == "__main__":
    test()
