""" Test extract vectors

.. author:: Thomas Cokelaer, Thomas.Cokelaer@inria.fr

"""
__revision__ = "$Id: test_extract_vectors.py 9885 2010-11-06 18:19:34Z cokelaer $"


from openalea.sequence_analysis import *
from tools import runTestClass


seq0 = Sequences(get_shared_data("chene_sessile_15pa.seq"))
vec10 = Vectors(seq0)
vec95 = ValueSelect(vec10, 1, 95)
seq20 = Sequences(get_shared_data("belren1.seq"))


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
