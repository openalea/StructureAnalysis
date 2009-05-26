"""extract_distribution tests

.. author:: Thomas Cokelaer, Thomas.Cokelaer@inria.fr

"""
__revision__ = "$Id:  $"


from openalea.stat_tool.data_transform import ExtractDistribution
from data import *



def test_hidden_semi_markov():
    for process in range(1, hsmc6.nb_output_process):
        for state in range(0,6):
            ExtractDistribution(hsmc6, "Observation", process, state)


def test_renewal():
    """tests not yet implemented"""
def test_semi_markov():
    """tests not yet implemented"""
def test_top_param():
    """tests not yet implemented"""
