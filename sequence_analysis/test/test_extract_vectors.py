
from data import seq20 as data
from openalea.sequence_analysis.data_transform import ExtractVectors

def test_extract_vectors_length():
    res = ExtractVectors(data, "Length") 


def test_extract_vectors_nb_run():
    res = ExtractVectors(data, "NbRun", 1)
    res = ExtractVectors(data, "NbRun", 1, 1)

def test_extract_vectors_occurrence():
    res = ExtractVectors(data, "NbOccurrence", 1)
    res = ExtractVectors(data, "NbOccurrence", 1, 1)
