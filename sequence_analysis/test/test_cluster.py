"""Cluster tests

.. author:: Thomas Cokelaer, Thomas.Cokelaer@inria.fr

.. todo:: check the AddVariable option (sequences) and sequences cases
"""

__revision__ = "$Id$"

import os

import pytest
from openalea.sequence_analysis.sequences import Sequences

# from openalea.sequence_analysis.semi_markov import SemiMarkov
from openalea.stat_tool.cluster import Cluster
from openalea.stat_tool.histogram import Histogram
from openalea.stat_tool.convolution import Convolution
from openalea.stat_tool.compound import Compound
from openalea.stat_tool.vectors import Vectors

from .tools import runTestClass, robust_path as get_shared_data


@pytest.fixture
def create_data_vectorsn():
    return Vectors([[1, 2, 3], [1, 3, 1], [4, 5, 6]])


def test_cluster_step_vectorsn(create_data_vectorsn):
    data = create_data_vectorsn
    cluster1 = data.cluster_step(1, 2)
    cluster2 = Cluster(data, "Step", 1, 2)
    assert str(cluster1) == str(cluster2)


def test_cluster_limit_vectorsn(create_data_vectorsn):
    data = create_data_vectorsn
    cluster1 = data.cluster_limit(1, [2, 4, 6])
    cluster2 = Cluster(data, "Limit", 1, [2, 4, 6])
    assert str(cluster1) == str(cluster2)


@pytest.fixture
def create_data_vector1():
    return Vectors([[1], [1], [4]])


def test_cluster_step_vector1(create_data_vector1):
    data = create_data_vector1
    cluster1 = data.cluster_step(1, 2)
    cluster2 = Cluster(data, "Step", 2)
    assert str(cluster1) == str(cluster2)


def test_cluster_limit_vector1(create_data_vector1):
    data = create_data_vector1
    cluster1 = data.cluster_limit(1, [2, 4, 6])
    cluster2 = Cluster(data, "Limit", [2, 4, 6])
    assert str(cluster1) == str(cluster2)


@pytest.fixture
def create_data_sequences1():
    return Sequences(str(get_shared_data("sequences1.seq")))


def test_cluster_step_sequences1(create_data_sequences1):
    data = create_data_sequences1
    mode = False
    cluster1 = data.cluster_step(1, 2, mode)
    cluster2 = Cluster(data, "Step", 2)
    assert str(cluster1) == str(cluster2)


def test_cluster_limit_sequences1(create_data_sequences1):
    data = create_data_sequences1
    print(data.nb_variable)
    cluster1 = data.cluster_limit(1, [2], False)
    cluster2 = Cluster(data, "Limit", [2], AddVariable=False)
    assert str(cluster1) == str(cluster2)


@pytest.fixture
def create_data_sequencen():
    return Sequences(str(get_shared_data("sequences2.seq")))


def test_cluster_step_sequencen(create_data_sequencen):
    data = create_data_sequencen
    mode = True
    cluster1 = data.cluster_step(1, 2, mode)
    cluster2 = Cluster(data, "Step", 1, 2)
    assert str(cluster1) == str(cluster2)


def test_cluster_limit_sequencen(create_data_sequencen):
    data = create_data_sequencen
    cluster1 = data.cluster_limit(1, [2], True)
    cluster2 = Cluster(data, "Limit", 1, [2])
    assert str(cluster1) == str(cluster2)
