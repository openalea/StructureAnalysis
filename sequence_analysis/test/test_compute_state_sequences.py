


from data import seq1, hmc
from openalea.sequence_analysis.data_transform import ComputeStateSequences



def test_ComputeStateSequences():
    ComputeStateSequences(seq1, hmc)
