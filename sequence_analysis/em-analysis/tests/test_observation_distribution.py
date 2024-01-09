"""Testing class ObservationDistribution."""
from ema import EyeMovementData
from ema import Model
from ema.observation_distribution import ObservationDistribution
from ema import TEST_PATH
import os
import pytest


@pytest.fixture
def model():
    return Model(EyeMovementData(), init_hsmc_file=os.path.join(TEST_PATH, 'input_model.hsmc'))


@pytest.fixture
def observation_distribution():
    return ObservationDistribution(None, 0, [0.2, 0.8])


def test_observation_distribution_constructor(observation_distribution):
    assert observation_distribution._observation_distribution_number == 0
    assert observation_distribution._outputs == [0.2, 0.8]


def test_observation_distribution_getter(observation_distribution):
    assert observation_distribution.observation_distribution_number == 0
    assert observation_distribution.outputs == [0.2, 0.8]


def test_observation_distribution_setter(model):
    model.parameters.output_processes[0].observation_distributions[0].outputs \
        = [0.1, 0.9, 0., 0., 0.]
    assert model.parameters.output_processes[0].observation_distributions[0]. \
        outputs == [0.1, 0.9, 0., 0., 0.]
