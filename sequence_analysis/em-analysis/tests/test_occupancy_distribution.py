"""Testing class OccupancyDistribution."""
from ema import EyeMovementData
from ema import Model
from ema import OccupancyDistribution
from ema import TEST_PATH
import numpy as np
import os
import pytest


@pytest.fixture
def model():
    return Model(EyeMovementData(), init_hsmc_file=os.path.join(TEST_PATH, 'input_model.hsmc'))


@pytest.fixture
def absorbing_state():
    return OccupancyDistribution(None, state_number=0, state_type='ABSORBING')


@pytest.fixture
def binomial_state():
    return OccupancyDistribution(
        None, state_number=1, state_type='RECURRENT',
        distribution_name='BINOMIAL', lower_bound=1, upper_bound=5,
        probability=0.3)


@pytest.fixture
def negative_binomial_state():
    return OccupancyDistribution(
        None, state_number=2, state_type='RECURRENT',
        distribution_name='NEGATIVE_BINOMIAL', lower_bound=2, parameter=5,
        probability=0.1)


@pytest.fixture
def poisson_state():
    return OccupancyDistribution(
        None, state_number=3, state_type='RECURRENT',
        distribution_name='POISSON', lower_bound=5, parameter=4)


def test_absorbing_state_constructor(absorbing_state):
    assert absorbing_state._state_number == 0
    assert absorbing_state._state_type == 'ABSORBING'
    assert absorbing_state._distribution_name is None
    assert absorbing_state._lower_bound is None
    assert absorbing_state._upper_bound is None
    assert absorbing_state._parameter is None
    assert absorbing_state._probability is None


def test_binomial_state_constructor(binomial_state):
    assert binomial_state._state_number == 1
    assert binomial_state._state_type == 'RECURRENT'
    assert binomial_state._distribution_name == 'BINOMIAL'
    assert binomial_state._lower_bound == 1
    assert binomial_state._upper_bound == 5
    assert binomial_state._parameter is None
    assert binomial_state._probability == 0.3


def test_negative_binomial_state_constructor(negative_binomial_state):
    assert negative_binomial_state._state_number == 2
    assert negative_binomial_state._state_type == 'RECURRENT'
    assert negative_binomial_state._distribution_name == 'NEGATIVE_BINOMIAL'
    assert negative_binomial_state._lower_bound == 2
    assert negative_binomial_state._upper_bound is None
    assert negative_binomial_state._parameter == 5
    assert negative_binomial_state._probability == 0.1


def test_poisson_state_constructor(poisson_state):
    assert poisson_state._state_number == 3
    assert poisson_state._state_type == 'RECURRENT'
    assert poisson_state._distribution_name == 'POISSON'
    assert poisson_state._lower_bound == 5
    assert poisson_state._upper_bound is None
    assert poisson_state._parameter == 4
    assert poisson_state._probability is None


def test_absorbing_state_getter(absorbing_state):
    assert absorbing_state.state_number == 0
    assert absorbing_state.state_type == 'ABSORBING'
    assert absorbing_state.distribution_name is None
    assert absorbing_state.lower_bound is None
    assert absorbing_state.upper_bound is None
    assert absorbing_state.parameter is None
    assert absorbing_state.probability is None


def test_binomial_state_getter(binomial_state):
    assert binomial_state.state_number == 1
    assert binomial_state.state_type == 'RECURRENT'
    assert binomial_state.distribution_name == 'BINOMIAL'
    assert binomial_state.lower_bound == 1
    assert binomial_state.upper_bound == 5
    assert binomial_state.parameter is None
    assert binomial_state.probability == 0.3


def test_negative_binomial_state_getter(negative_binomial_state):
    assert negative_binomial_state.state_number == 2
    assert negative_binomial_state.state_type == 'RECURRENT'
    assert negative_binomial_state.distribution_name == 'NEGATIVE_BINOMIAL'
    assert negative_binomial_state.lower_bound == 2
    assert negative_binomial_state.upper_bound is None
    assert negative_binomial_state.parameter == 5
    assert negative_binomial_state.probability == 0.1


def test_poisson_state_getter(poisson_state):
    assert poisson_state.state_number == 3
    assert poisson_state.state_type == 'RECURRENT'
    assert poisson_state.distribution_name == 'POISSON'
    assert poisson_state.lower_bound == 5
    assert poisson_state.upper_bound is None
    assert poisson_state.parameter == 4
    assert poisson_state.probability is None


def test_absorbing_state_setter(model):
    model.parameters.occupancy_distributions[0].state_type = 'ABSORBING'
    transition_probabilities_line = np.zeros(model.parameters.occupancy_distributions[0]._model._k)
    transition_probabilities_line[0] = 1
    assert model.parameters.occupancy_distributions[0].distribution_name is None
    assert model.parameters.occupancy_distributions[0].lower_bound is None
    assert model.parameters.occupancy_distributions[0].upper_bound is None
    assert model.parameters.occupancy_distributions[0].parameter is None
    assert model.parameters.occupancy_distributions[0].probability is None
    assert (model.parameters.transition_probabilities.transition_probabilities[0, :] ==
            transition_probabilities_line).all()


def test_recurrent_state_setter(model):
    model.parameters.occupancy_distributions[0].state_type = 'RECURRENT'
    model.parameters.occupancy_distributions[0].distribution_name = 'BINOMIAL'
    model.parameters.occupancy_distributions[0].lower_bound = 1
    model.parameters.occupancy_distributions[0].upper_bound = 20
    model.parameters.occupancy_distributions[0].parameter = None
    model.parameters.occupancy_distributions[0].probability = 0.2
    assert model.parameters.occupancy_distributions[0].state_type == 'RECURRENT'
    assert model.parameters.occupancy_distributions[0].distribution_name == 'BINOMIAL'
    assert model.parameters.occupancy_distributions[0].lower_bound == 1
    assert model.parameters.occupancy_distributions[0].upper_bound == 20
    assert model.parameters.occupancy_distributions[0].parameter is None
    assert model.parameters.occupancy_distributions[0].probability == 0.2
