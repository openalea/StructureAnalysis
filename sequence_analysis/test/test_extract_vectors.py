""" Test extract vectors

.. author:: Thomas Cokelaer, Thomas.Cokelaer@inria.fr

"""
__revision__ = "$Id$"


from data import seq20 
from openalea.sequence_analysis.data_transform import ExtractVectors
from tools import runTestClass


class Test():
    def __init__(self):
        self.data = seq20
        
    def test_extract_vectors_length(self):
        data = self.data
        res = ExtractVectors(data, "Length") 

    def test_extract_vectors_nb_run(self):
        data = self.data
        res = ExtractVectors(data, "NbRun", 1)
        res = ExtractVectors(data, "NbRun", 1, 1)

    def test_extract_vectors_occurrence(self):
        data = self.data
        res = ExtractVectors(data, "NbOccurrence", 1)
        res = ExtractVectors(data, "NbOccurrence", 1, 1)




if __name__ == "__main__":
    runTestClass(Test())
