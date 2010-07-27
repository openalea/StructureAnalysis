"""
provides data samples to play with in test and documentation
"""

import os
import openalea.sequence_analysis as sa


path = sa.__path__[0] + (os.sep + '..') * 3 \
    + os.sep + 'share' + os.sep + 'data' +os.sep
path = os.path.realpath(path) + os.sep

class _data_tutorial():
    def __init__(self):
        from openalea.sequence_analysis import Sequences
        self.seq1 = Sequences([1,2,3,4,5,6,7,8,9,10])


data_tutorial = _data_tutorial()

