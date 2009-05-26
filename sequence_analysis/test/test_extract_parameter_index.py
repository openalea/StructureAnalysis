""" extract parameter functions


"""
__revision__ = "$Id: $"


from openalea.sequence_analysis.data_transform import IndexParameterExtract
from openalea.sequence_analysis.sequences import Sequences



def test1():
    seq69 = Sequences("data//pin_laricio_7x.seq")
    IndexParameterExtract(seq69, 1928, 1929)

