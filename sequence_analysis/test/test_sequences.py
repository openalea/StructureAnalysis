"""Sequences tests

.. author:: Thomas Cokelaer, Thomas.Cokelaer@inria.fr

"""

__revision__ = "$Id$"


from openalea.stat_tool import _stat_tool
import pytest
from openalea.sequence_analysis import _sequence_analysis
from openalea.sequence_analysis.sequences import Sequences, Split, SaveMTG
from openalea.sequence_analysis.data_transform import (
    Cumulate,
    Difference,
    ExtractVectors,
    IndexParameterExtract,
    IndexParameterSelect,
    RecurrenceTimeSequences,
)

from openalea.stat_tool.data_transform import *
from openalea.stat_tool.cluster import Cluster
from openalea.stat_tool.cluster import Transcode, Cluster

from .tools import interface
from .tools import runTestClass

from openalea.sequence_analysis.sequences import Sequences, IndexParameterType
from .tools import robust_path as get_shared_data


@pytest.fixture
def build_data():
    """todo: check identifier output. should be a list"""
    # build a list of 2 sequences with a variable that should be identical
    # to sequences1.seq
    data = Sequences(
        [
            [
                1,
                0,
                0,
                0,
                1,
                1,
                2,
                0,
                2,
                2,
                2,
                1,
                1,
                0,
                1,
                0,
                1,
                1,
                1,
                1,
                0,
                1,
                1,
                1,
                0,
                1,
                2,
                2,
                2,
                1,
            ],
            [0, 0, 0, 1, 1, 0, 2, 0, 2, 2, 2, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0],
        ]
    )
    assert data

    assert data.nb_sequence == 2
    assert data.nb_variable == 1
    assert data.cumul_length == 52
    assert data.max_length == 30

    assert [0, 1] == data.get_identifiers()

    return data


@pytest.fixture
def build_seqn():
    return Sequences([[[1, 1, 1], [12, 12, 12]], [[2, 2, 2], [22, 23, 24]]])


@pytest.fixture
def build_seq1():
    return Sequences([[1, 1, 1], [2, 2, 2]])


@pytest.fixture
def build_seq_realn():
    return Sequences(
        [
            [[1.5, 1.5, 1.5], [12.5, 12.5, 12.5]],
            [[2.5, 2.5, 2.5], [22.5, 23.5, 24.5]],
        ]
    )


def build_seq_wrong_identifiers():
    try:
        s = Sequences([[1, 1, 1], [2, 2, 2]], Identifiers=[-1, 1])
        assert False
    except:
        assert True


def test_len(build_data):
    seq = build_data
    assert len(seq) == 2
    assert len(seq) == seq.nb_sequence


def test_extract(build_seqn):
    # todo
    seqn = build_seqn
    assert seqn.extract_value(1)
    assert seqn.extract_value(2)


def test_index_parameter_type():
    seq1 = Sequences([[1.0, 1, 1], [2.0, 2, 2.0]])
    assert IndexParameterType(seq1) == "IMPLICIT_TYPE"
    seq1 = Sequences([[1.0, 1, 1], [2.0, 2, 2.0]], IndexParameterType="TIME")
    assert IndexParameterType(seq1) == "TIME"
    seq1 = Sequences([[1.0, 1, 1], [2.0, 2, 2.0]], IndexParameterType="POSITION")
    assert IndexParameterType(seq1) == "POSITION"


def test_constructors():
    # heterogeneous or homogeneous type
    seq1 = Sequences([1, 2, 3, 4])
    seq1 = Sequences([1, 2, 3, 4.0])
    assert seq1.nb_sequence == 1
    assert seq1.nb_variable == 1

    # single sequence multivariate
    seq2 = Sequences([[1, 2], [3, 4], [5, 6]])
    assert seq2.nb_sequence == 1
    assert seq2.nb_variable == 2
    # ambiguous case (length>5)
    seq2 = Sequences([[1, 2, 3, 4, 5, 6], [3, 4, 3, 4, 5, 6], [5, 6, 4, 5, 6, 7]])
    assert seq2.nb_sequence == 3
    assert seq2.nb_variable == 1

    # univariates sequences
    seq3 = Sequences([[1, 2], [3, 4], [5, 6, 7]])
    assert seq3.nb_sequence == 3
    assert seq3.nb_variable == 1

    # general case
    seq4 = Sequences(
        [[[1, 2], [3, 4]], [[21, 22], [23, 24]], [[31, 32], [33, 34], [35, 36]]]
    )
    assert seq4.nb_sequence == 3
    assert seq4.nb_variable == 2

    seq4 = Sequences(
        [[[1, 2], [3, 4]], [[21, 22], [23, 24]], [[31, 32], [33, 34], [35, 36]]],
        VertexIdentifiers=[[1, 2], [3, 4], [5, 6, 7]],
        Identifiers=[1, 2, 3],
    )

    seq4 = Sequences(
        [[[1, 2], [3, 4]], [[21, 22], [23, 24]], [[31, 32], [33, 34], [35, 36]]],
        IndexParameterType="POSITION",
        IndexParameter=[[0, 1, 10], [2, 3, 11], [4, 5, 6, 12]],
    )

    seq4 = Sequences(
        [[[1, 2], [3, 4]], [[21, 22], [23, 24]], [[31, 32], [33, 34], [35, 36]]],
        IndexParameterType="TIME",
        IndexParameter=[[0, 1], [2, 3], [4, 5, 6]],
    )


def test_constructor_one_sequence():
    # two sequences with 2 variables (int)
    s = Sequences([[1, 1, 1], [2, 2, 2]])
    assert s
    # two sequences with 2 variables (real)
    s = Sequences([[1.0, 1.0, 1.0], [2.0, 2.0, 2.0]])
    assert s
    # two sequences with 2 variables mix of (real and int)
    # works because the first one is float so all others are assume to be float as well
    s = Sequences([[1.0, 1.0, 1.0], [2.0, 2.0, 2]])
    assert s
    # here it fails because the first number is int but others may be float
    try:
        s = Sequences([[1, 1.0, 1.0], [2.0, 2.0, 2]])
        assert False
    except:
        assert True


def test_constructor_two_sequences():
    # two sequences with 2 variables (int)
    s = Sequences([[[1, 1, 1], [12, 12, 12]], [[2, 2, 2], [22, 23, 24]]])
    assert s
    s = Sequences(
        [[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]], [[2.0, 2.0, 2.0], [2.0, 2.0, 2.0]]]
    )
    assert s


def test_container(build_seq1, build_seqn):
    # first index is the sequence and second index is the variable
    s = build_seqn
    assert s[0, 0] == [1, 1, 1]
    assert s[0, 1] == [12, 12, 12]
    assert s[1, 0] == [2, 2, 2]
    assert s[1, 1] == [22, 23, 24]
    assert s[1, 1][1] == 23
    assert len(s) == 2

    s = build_seq1
    assert s[0, 0] == [1, 1, 1]
    assert s[0, 1] == [2, 2, 2]
    assert s[0, 0][0] == 1
    assert len(s) == 1


def test_value_select(build_seqn):
    "test_value_select implemented but need to be checked"
    seqn = build_seqn
    a = seqn.value_select(1, 1, 2, True)
    assert a
    assert str(ValueSelect(seqn, 1, 1, 2)) == str(seqn.value_select(1, 1, 2, True))


def test_select_variable_int(build_seqn):
    "test_select_variable_int implemented but need to be checked (index issue)"
    # !!!!!!!! NEED to CHECK THE INDEX 0, 1 , ... or 1,2,....
    # what about identifiers ?
    # Variable seems to start at 1 not 0
    s = build_seqn
    # select variable 1
    select = s.select_variable([1], keep=True)
    assert select[0, 0] == [1]
    assert select[1, 0] == [2]


def test_select_variable_real(build_seq_realn):
    "test_select_variable_real implemented but need to be checked (index issue)"
    # !!!!!!!! NEED to CHECK THE INDEX 0, 1 , ... or 1,2,....
    # what about identifiers ?
    # Variable seems to start at 1 not 0
    s = build_seq_realn
    # select variable 1
    select = s.select_variable([1], keep=True)
    assert select[0, 0] == [1.5]
    assert select[1, 0] == [2.5]


def test_select_individual(build_seqn):
    # select one or several sequences
    s = build_seqn

    # select all
    select = s.select_individual([0, 1], keep=True)
    assert s.display() == select.markovian_sequences().display()

    # select first sequence only
    select = s.select_individual([0], keep=True)
    assert select[0, 0] == [1, 1, 1]
    assert select[0, 1] == [12, 12, 12]
    try:
        select[1, 0]
        assert False
    except:
        assert True

    # select second sequence only
    select = s.select_individual([1], keep=True)
    assert select[0, 0] == [2, 2, 2]
    assert select[0, 1] == [22, 23, 24]
    try:
        select[1, 0]
        assert False
    except:
        assert True


def test_shift_seqn(build_seqn):
    s = build_seqn
    shifted = s.shift(1, 2)
    assert shifted[0, 0] == [3, 1, 1]
    assert shifted[0, 1] == [14, 12, 12]
    assert shifted[1, 0] == [4, 2, 2]
    assert shifted[1, 1] == [24, 23, 24]


def test_shift_seq1(build_seq1):
    s = build_seq1
    shifted = s.shift(1, 2)
    assert shifted[0, 0] == [3, 1, 1]
    assert shifted[0, 1] == [4, 2, 2]


def test_threshold_seq1(build_seqn):
    s = build_seqn
    thresholded = s.thresholding(1, 10, "ABOVE")
    for x in thresholded:
        for v in x:
            assert v[0] <= 10


def test_threshold_seq():
    s = Sequences(
        [
            [[1.01, 1.07], [2.01, 2.07], [1.99, 1.07], [2.41, 2.07]],
            [[1.97, 1.07], [1.98, 2.07], [1.99, 1.07], [2.00, 2.07]],
        ]
    )
    thresholded = s.thresholding(1, 1.99, "ABOVE")
    for x in thresholded:
        for v in x:
            assert v[0] <= 1.99
    thresholded = s.thresholding(1, 1.99, "BELOW")
    for x in thresholded:
        for v in x:
            assert v[0] >= 1.99

    assert s


def test_merge(build_seqn, build_seq1):
    s1 = build_seqn
    s2 = build_seqn
    s3 = build_seq1

    sall = s1.merge([s2])
    assert sall.nb_sequence == 4
    assert sall.nb_variable == 3

    sall = s1.merge([s3])
    assert sall.nb_sequence == 3
    assert sall.nb_variable == 3


def test_merge_and_Merge(build_seqn):
    s1 = build_seqn
    s2 = build_seqn

    a = s1.merge([s2])
    b = s2.merge([s1])
    v = Merge(s1, s2)

    assert str(a) == str(b)
    assert str(a) == str(v)


def test_imcompatible_merge(build_seqn, build_seq1):
    s1 = build_seqn
    s3 = build_seq1
    try:
        s1.merge(s3)
        assert False
    except:
        assert True


def test_merge_variable(build_seqn, build_seq1):
    s1 = build_seqn
    s2 = build_seqn
    s3 = build_seq1

    sall = s1.merge_variable([s2], 1)  # why 1 ? same result with 2 !
    assert sall.nb_sequence == 2
    assert sall.nb_variable == 6

    # sall =  s1.merge_variable([s3],1)
    # assert sall.nb_sequence == 2
    # assert sall.nb_variable == 3


def test_merge_variable_and_MergeVariable(build_seqn, build_seq1):
    s1 = build_seqn
    s2 = build_seqn
    s3 = build_seq1

    a = s1.merge_variable([s2], 1)
    b = s2.merge_variable([s1], 1)
    v = MergeVariable(s1, s2)

    assert str(a) == str(b)
    assert str(a) == str(v)


def test_cluster_step():
    seq1 = Sequences([[1, 2, 3], [1, 3, 1], [4, 5, 6]])
    assert str(Cluster(seq1, "Step", 1, 2)) == str(seq1.cluster_step(1, 2, True))
    seqn = Sequences([[[1, 2, 3], [1, 3, 1]], [[4, 5, 6], [7, 8, 9]]])
    assert str(Cluster(seqn, "Step", 1, 2)) == str(seqn.cluster_step(1, 2, True))


def test_cluster_limit():
    seq1 = Sequences([[1, 2, 3], [1, 3, 1], [4, 5, 6]])
    assert str(Cluster(seq1, "Limit", 1, [2])) == str(seq1.cluster_limit(1, [2], True))
    seqn = Sequences([[[1, 2, 3], [1, 3, 1]], [[4, 5, 6], [7, 8, 9]]])
    assert str(Cluster(seqn, "Limit", 1, [2, 4, 6])) == str(
        seqn.cluster_limit(1, [2, 4, 6], True)
    )


def test_transcode(build_seqn):
    """This functionality need to be checked.

    See also the vector case!"""
    seq = build_seqn
    assert str(
        seq.transcode(
            1,
            [0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0o1, 1, 1, 1, 1, 0, 0],
            False,
        )
    ) == str(
        Transcode(
            seq,
            1,
            [0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0o1, 1, 1, 1, 1, 0, 0],
        )
    )


def test_reverse(build_seqn):
    """reverse to be checked. seems to give same output as input"""
    s = build_seqn
    s.reverse()


def test_max_length(build_data):
    s = build_data
    assert s.max_length == 30


def test_get_max_value(build_data):
    s = build_data
    assert s.get_max_value(0) == 2


def test_get_min_value(build_data):
    s = build_data
    assert s.get_min_value(0) == 0


def test_get_length(build_data):
    s = build_data
    assert s.get_length(0) == 30
    assert s.get_length(1) == 22


def test_difference(build_data):
    data = build_data
    assert str(Difference(data, 1)) == str(data.difference(1, False))
    res = Difference(data, 1)
    assert res.cumul_length == 50


def test_cumulate(build_data):
    # see also test_cumulate for more tests
    s = build_data
    res = Cumulate(s)
    assert res.cumul_length == 52


def test_extract_vectors(build_data):
    """see test_extract_vectors"""
    ExtractVectors(build_data, "Length")


def _test_index_parameter_extract(build_data):
    """fixme: markovian_sequences should be in wrapper ?"""
    aml = build_data.index_parameter_extract(0, 29).markovian_sequences()
    mod = IndexParameterExtract(build_data, 0, MaxIndex=29)
    assert str(aml) == str(mod)


def test_index_parameter_select():
    """test to be done"""
    pass


def test_recurrence_time_sequences(build_data):
    aml = build_data.recurrence_time_sequences(1, 1)
    mod = RecurrenceTimeSequences(build_data, 1, 1)
    assert str(aml.markovian_sequences()) == str(mod)


def test_remove_run():
    """test to be done"""
    pass


def test_transform_position():
    """test to be done"""
    pass


def test_segmentation_extract():
    """test to be done"""
    pass


def test_variable_scaling():
    """test to be done"""
    pass


def test_remove_index_parameter():
    """test to be done"""
    pass


def test_write_mtg(build_data):
    import os

    build_data.mtg_write("test.mtg", [1, 2])
    os.remove("test.mtg")


def test_SaveMTG(build_data):
    SaveMTG(build_data, Filename="test.mtg", Type=["N"])


def test_split():
    # markovian sequences
    data = Sequences(str(get_shared_data("vanille_m.seq")))
    Split(data, 2)


def test_initial_run():
    from openalea.sequence_analysis import ComputeInitialRun

    # markovian sequences
    data = Sequences(str(get_shared_data("vanille_m.seq")))
    ComputeInitialRun(data)


"""

seq.segmentation_extract
seq.pointwise_average
seq.cross
seq.get_index_parameter_type
seq.sojourn_time_sequences
seq.remove_index_parameter
seq.round
"""
