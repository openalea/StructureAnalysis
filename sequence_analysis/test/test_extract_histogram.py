""" Test extract histogram

.. author:: Thomas Cokelaer, Thomas.Cokelaer@inria.fr

.. todo:: to be done
"""
__revision__ = "$Id$"


from openalea.stat_tool.data_transform import ValueSelect, ExtractHistogram

from data import *


def test_vectors():
    ExtractHistogram(vec95, 2)

def test_sequences():
    ExtractHistogram(seq20, "Recurrence", 1)
    ExtractHistogram(seq20, "Recurrence", 2)
    ExtractHistogram(seq20, "Length")


    ExtractHistogram(seq11, "FirstOccurrence", 1, 0)


    ExtractHistogram(seq0, "Value", 3)

    ExtractHistogram(seq70, "Value", 1)
    ExtractHistogram(seq70, "Value", 2)

def test_time_events():
    """not implemented"""

def test_renewal():
    """not implemented"""

