"""test extracted from stat_tool_test

#########################################################################
#
#  Distance matrix and clustering (partitioning or hierarchical methods) 
#
#########################################################################

"""


from openalea.stat_tool import Vectors
from openalea.stat_tool import Compare
from openalea.stat_tool import SelectVariable
from openalea.stat_tool import ValueSelect
from openalea.stat_tool import VectorDistance
from openalea.stat_tool import Plot
from openalea.stat_tool import Clustering
from openalea.stat_tool import SelectIndividual
from openalea.stat_tool import ExtractHistogram
from openalea.stat_tool import Display

from openalea.stat_tool import get_shared_data

def test():
    vec10 = Vectors(get_shared_data('chene_sessile.vec'))
    vec15 = SelectVariable(vec10, [1, 3, 6], Mode="Reject")
    vec9596 = ValueSelect(vec10, 1, 1995, 1996)

    # computation of a distance matrix using a standardization procedure

    matrix10 = Compare(vec15, VectorDistance("N", "N", "N"))
    Plot(matrix10)
    # clustering using a partitioning method

    cluster = Clustering(matrix10, "Partition", 5)
    Plot(cluster)
    Display(Clustering(matrix10, "Partition", 2))

    vec151 = SelectIndividual(vec10,  \
        [69, 48, 41, 44, 32, 47, 81, 95, 11, 36, 75, 108, 56, 83, 38, 98
        , 113, 134, 110, 101, 77, 35, 74, 80, 50, 24, 89, 128, 5, 45, 8,
        116, 119, 132, 61, 78, 53, 29, 131, 65, 90, 96, 104, 20, 86, 66,
        42, 68, 125, 14, 23, 54, 33, 26, 71, 129, 102, 51, 70, 111, 138,
        19, 127, 62, 117, 137, 2, 28, 17])
    vec152 = SelectIndividual(vec10, \
        [100, 13, 133, 105, 72, 9, 93, 109, 30, 115, 63, 7, 55, 37, 15,
        114, 106, 46, 73, 18, 3, 87, 58, 43, 60, 76, 52, 6, 39, 31, 12, 
        99, 121, 123, 22, 79, 94, 88, 21, 97, 25, 40, 57, 136, 67, 49, 10, 
        4, 120, 92, 27, 91, 64, 124, 16, 130, 84, 107, 126, 103, 122, 112, 
        59, 1, 82, 34, 135, 118, 85])

    Plot(ExtractHistogram(vec151, 4), ExtractHistogram(vec152, 4))

    matrix11 = Compare(vec15, VectorDistance("N", "O", "N"))

    Clustering(matrix10, "Hierarchy", Algorithm="Agglomerative")
    Clustering(matrix10, "Hierarchy", Algorithm="Divisive")

    vec16 = SelectVariable(vec9596, [1, 3], Mode="Reject")
    matrix12 = Compare(vec16, VectorDistance("N", "N", "N", "N"))
    matrix13 = Compare(vec16, VectorDistance("N", "O", "N", "N"))

