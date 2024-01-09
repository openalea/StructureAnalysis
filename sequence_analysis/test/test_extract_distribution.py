"""extract_distribution tests

.. author:: Thomas Cokelaer, Thomas.Cokelaer@inria.fr

"""
__revision__ = "$Id: test_extract_distribution.py 9885 2010-11-06 18:19:34Z cokelaer $"


from openalea.stat_tool.data_transform import ExtractDistribution
from openalea.sequence_analysis import *



def test_hidden_semi_markov():
    seq = Sequences(get_shared_data("pin_laricio_7x.seq"))
    seq_cluster = Cluster(seq, "Step", 1, 10)

    hsmc60 = HiddenSemiMarkov(get_shared_data('pin_laricio_6.hsc'))
    hsmc6 = Estimate(seq_cluster, "HIDDEN_SEMI-MARKOV", hsmc60)

    for process in range(1, hsmc6.nb_output_process):
        for state in range(0,6):
            ExtractDistribution(hsmc6, "Observation", process, state)


def test_renewal():
    """tests not yet implemented"""
def test_semi_markov():
    """tests not yet implemented"""
def test_top_param():
    """tests not yet implemented"""
