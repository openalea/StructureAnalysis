"""Simulate semi Markovian data
"""

from ema import EyeMovementData
from ema import Model
from ema import SequentialData
from openalea.sequence_analysis import HiddenSemiMarkov
import random

random.seed(1)

input_model_file = "./input_model2.hsmc"
output_model_file = "./output_model2.hsmc"

hsmm = HiddenSemiMarkov(input_model_file)

simulated_data = hsmm.simulation_markovian_sequences(
    1000, EyeMovementData().get_input_sequence, False)

text_reading_list = []
for text_reading in simulated_data:
    fixation_list = []
    for fixation in text_reading:
        fixation_list.append(fixation[0])
    text_reading_list.append(fixation_list)
simulated_data = text_reading_list

model = Model(SequentialData(simulated_data), random_init=True, k=5, n_init=2, n_random_seq=30, n_iter_init=10)
model.hsmm.save(output_model_file)
