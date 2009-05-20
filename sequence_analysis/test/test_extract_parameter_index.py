""" Cluster tests

.. todo:: check the AddVariable option (sequences) and sequences cases
"""
__revision__ = "$Id: test_cluster.py 6325 2009-04-29 16:20:55Z cokelaer $"


from openalea.sequence_analysis.data_transform import IndexParameterExtract
from openalea.sequence_analysis.sequences import Sequences



def test1():
    seq69 = Sequences("data//pin_laricio_7x.seq")
    IndexParameterExtract(seq69, 1928, 1929)

