"""Testing class model."""
from ema import EyeMovementData
from ema import Model
from ema import TEST_PATH
import numpy as np
import os
import pytest


@pytest.fixture
def model_init_by_default():
    return Model(EyeMovementData(), k=5)


@pytest.fixture
def model_init_with_hsmc_file():
    return Model(EyeMovementData(), init_hsmc_file=os.path.join(TEST_PATH, 'input_model.hsmc'))


@pytest.fixture
def model_init_random():
    return Model(EyeMovementData(), random_init=True, k=5, n_init=2,
                 n_random_seq=30, n_iter_init=10)


def test_model_init_by_default_constructor(model_init_by_default):
    assert model_init_by_default._output_process_name == 'READMODE'
    assert model_init_by_default._model_type == 'HIDDEN_SEMI-MARKOV_CHAIN'
    assert model_init_by_default._k == 5
    assert not model_init_by_default._random_init
    assert model_init_by_default._n_init is None
    assert model_init_by_default._n_iter_init is None
    assert len(model_init_by_default.parameters.initial_probabilities.
               initial_probabilities) == 5
    assert model_init_by_default.parameters.transition_probabilities.\
        transition_probabilities.shape == (5, 5)
    assert len(model_init_by_default.parameters.occupancy_distributions) == 5
    assert len(model_init_by_default.parameters.output_processes) == 1
    assert len(model_init_by_default.parameters.output_processes[0].
               observation_distributions) == 5
    assert len(model_init_by_default.parameters.output_processes[0].
               observation_distributions[0].outputs) == 5
    assert len(model_init_by_default.parameters.output_processes[0].
               observation_distributions[1].outputs) == 5
    assert len(model_init_by_default.parameters.output_processes[0].
               observation_distributions[2].outputs) == 5
    assert len(model_init_by_default.parameters.output_processes[0].
              observation_distributions[3].outputs) == 5
    assert len(model_init_by_default.parameters.output_processes[0].
               observation_distributions[4].outputs) == 5


def test_model_init_with_hsmc_file(model_init_with_hsmc_file):
    assert model_init_with_hsmc_file._output_process_name == 'READMODE'
    assert model_init_with_hsmc_file._model_type == 'HIDDEN_SEMI-MARKOV_CHAIN'
    assert model_init_with_hsmc_file._k == 5
    assert not model_init_with_hsmc_file._random_init
    assert model_init_with_hsmc_file._n_init is None
    assert model_init_with_hsmc_file._n_iter_init is None
    assert len(model_init_with_hsmc_file.parameters.initial_probabilities.
               initial_probabilities) == 5
    assert model_init_with_hsmc_file.parameters.transition_probabilities.\
        transition_probabilities.shape == (5, 5)
    assert len(model_init_with_hsmc_file.parameters.occupancy_distributions) == 5
    assert len(model_init_with_hsmc_file.parameters.output_processes) == 1
    assert len(model_init_with_hsmc_file.parameters.output_processes[0].
               observation_distributions) == 5
    assert len(model_init_with_hsmc_file.parameters.output_processes[0].
               observation_distributions[0].outputs) == 5
    assert len(model_init_with_hsmc_file.parameters.output_processes[0].
               observation_distributions[1].outputs) == 5
    assert len(model_init_with_hsmc_file.parameters.output_processes[0].
               observation_distributions[2].outputs) == 5
    assert len(model_init_with_hsmc_file.parameters.output_processes[0].
              observation_distributions[3].outputs) == 5
    assert len(model_init_with_hsmc_file.parameters.output_processes[0].
               observation_distributions[4].outputs) == 5


def test_model_init_random(model_init_random):
    assert model_init_random._output_process_name == 'READMODE'
    assert model_init_random._model_type == 'HIDDEN_SEMI-MARKOV_CHAIN'
    assert model_init_random._k == 5
    assert model_init_random._random_init
    assert model_init_random._n_init == 2
    assert model_init_random._n_iter_init == 10
    assert len(model_init_random.parameters.initial_probabilities.
               initial_probabilities) == 5
    assert model_init_random.parameters.transition_probabilities.\
        transition_probabilities.shape == (5, 5)
    assert len(model_init_random.parameters.occupancy_distributions) == 5
    assert len(model_init_random.parameters.output_processes) == 1
    assert len(model_init_random.parameters.output_processes[0].
               observation_distributions) == 5

    """ Not necessarily of length 5 (may have 0s)
    assert len(model_init_random.parameters.output_processes[0].
               observation_distributions[0].outputs) == 5
    assert len(model_init_random.parameters.output_processes[0].
               observation_distributions[1].outputs) == 5
    assert len(model_init_random.parameters.output_processes[0].
               observation_distributions[2].outputs) == 5
    assert len(model_init_random.parameters.output_processes[0].
              observation_distributions[3].outputs) == 5
    assert len(model_init_random.parameters.output_processes[0].
               observation_distributions[4].outputs) == 5
    """



def test_update_hsmc_file(model_init_by_default):
    model_init_by_default.update_hsmc_file()


def test_update_parameters(model_init_by_default):
    model_init_by_default.update_parameters()


def test_iterate_EM(model_init_by_default):
    n_iter = model_init_by_default._n_iter
    model_init_by_default.iterate_em(10)
    n_iter += 10
    assert model_init_by_default._n_iter == n_iter


def test_update_restored_data(model_init_by_default):
    model_init_by_default.update_restored_data()


def test_print_hsmc_file(model_init_by_default):
    model_init_by_default.print_hsmc_file(verbose=False)


def test_generate_random_sequences(model_init_by_default):
    assert len(model_init_by_default.generate_random_sequences(10)) == 10


def test_secure_probabilities_sum(model_init_by_default):
    Model.secure_probabilities_sum(model_init_by_default._hsmc_file)
    np.testing.assert_almost_equal(
        np.sum(model_init_by_default.parameters.initial_probabilities.\
               initial_probabilities), 1)
    np.testing.assert_almost_equal(
        np.sum(model_init_by_default.parameters.transition_probabilities.\
               transition_probabilities, 1), np.array([1, 1, 1, 1, 1]))
