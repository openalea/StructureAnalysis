"""Testing class OutputProcess."""
from ema import EyeMovementData
from ema import Model
from ema.observation_distribution import ObservationDistribution
from ema import OutputProcess
from ema import TEST_PATH
import os
import pytest


@pytest.fixture
def model():
    return Model(EyeMovementData(), init_hsmc_file=os.path.join(TEST_PATH, 'input_model.hsmc'))


@pytest.fixture
def output_process():
    return OutputProcess(
        None, 1, 'CATEGORICAL',
        [ObservationDistribution(None, 1, [0.2, 0.8]),
         ObservationDistribution(None, 2, [0.1, 0.9])])


def test_output_process_constructor(output_process):
    assert(output_process._output_process_number == 1)
    assert(output_process._output_process_type == 'CATEGORICAL')
    assert(len(output_process._observation_distributions) == 2)


def test_output_process_getter(output_process):
    assert(output_process.output_process_number == 1)
    assert(output_process.output_process_type == 'CATEGORICAL')
    assert(len(output_process._observation_distributions) == 2)


""" Does not make sense
def test_output_process_setter(model):
    model.parameters.output_processes[0].output_process_number = 2
    assert model.parameters.output_processes[0].output_process_number == 2
"""