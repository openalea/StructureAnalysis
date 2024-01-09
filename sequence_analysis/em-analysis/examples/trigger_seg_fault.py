from ema import EyeMovementData
from ema import Model
from ema import SequentialData
from ema import MODELS_PATH
import numpy as np
from numpy.random import dirichlet
from openalea.sequence_analysis import HiddenSemiMarkov
import os
import random
import tempfile

random.seed(1)
np.random.seed(2)

# seed 8 file error
# seed 11 seg fault

def generate_random_transition_matrix(size_of_hidden_state_space):
    A = np.zeros(size_of_hidden_state_space**2).reshape(size_of_hidden_state_space, size_of_hidden_state_space)
    A_temp = dirichlet(np.ones(size_of_hidden_state_space-1), size_of_hidden_state_space)
    for i in xrange(0, size_of_hidden_state_space):
        diag_pass = False
        for j in xrange(0, size_of_hidden_state_space):
            if i == j:
                diag_pass = True
            if i != j:
                if diag_pass:
                    A[i, j] = A_temp[i, j-1]
                else:
                    A[i, j] = A_temp[i, j]
    return A


def generate_random_sojourn_distributions(size_of_hidden_state_space):
    distribution_list = ['geometric', 'binomial', 'poisson', 'negative binomial']
    sojourn_distribution_list = []
    distribution = np.random.randint(0, len(distribution_list), size_of_hidden_state_space)
    for i in xrange(0, size_of_hidden_state_space):
        if distribution_list[distribution[i]] == 'geometric':
            probability = np.random.uniform(0, 1, 1)[0]
            sojourn_distribution_list.append(['geometric', probability])
        elif distribution_list[distribution[i]] == 'binomial':
            upper_bound = np.random.randint(6, 12, 1)[0]
            probability = np.random.uniform(0, 1, 1)[0]
            sojourn_distribution_list.append(['binomial', upper_bound, probability])
        elif distribution_list[distribution[i]] == 'poisson':
            parameter = np.random.uniform(1, 10, 1)[0]
            sojourn_distribution_list.append(['poisson', parameter])
        elif distribution_list[distribution[i]] == 'negative binomial':
            parameter = np.random.randint(1, 5, 1)[0]
            probability = np.random.uniform(0, 1, 1)[0]
            sojourn_distribution_list.append(['negative binomial', parameter, probability])
    return sojourn_distribution_list


def write_input_file(pi, A, B, sojourn_distribution_list):
    tmp_path = tempfile.mkdtemp(prefix='ema-tmp-dir-')
    model_id = tmp_path[-6:]
    hsmc_file = os.path.join(MODELS_PATH, model_id + '.hsmc')
    k = len(pi)
    l = B.shape[1]
    with open(hsmc_file, 'w') as f:
        f.write('HIDDEN_SEMI-MARKOV_CHAIN\n\n')
        f.write(str(k) + ' STATES\n\n')
        f.write('INITIAL_PROBABILITIES\n')
        for i in xrange(0, k):
            f.write(str(pi[i]))
            if i != k - 1:
                f.write('     ')
        f.write('\n\n')
        f.write('TRANSITION_PROBABILITIES\n')
        for i in xrange(0, k):
            for j in xrange(0, k):
                f.write(str(A[i, j]))
                if j != k - 1:
                    f.write('     ')
            f.write('\n')
        f.write('\n\n')
        for i in xrange(0, k):
            f.write('STATE ' + str(i) + ' OCCUPANCY_DISTRIBUTION\n')
            if sojourn_distribution_list[i][0] == 'geometric':
                f.write('NEGATIVE_BINOMIAL INF_BOUND : 1 PARAMETER : 1 PROBABILITY : '
                        + str(sojourn_distribution_list[i][1]) + '\n\n')
            elif sojourn_distribution_list[i][0] == 'binomial':
                f.write('BINOMIAL INF_BOUND : 1 SUP_BOUND : '
                        + str(sojourn_distribution_list[i][1])
                        + ' PROBABILITY : ' + str(sojourn_distribution_list[i][2]) + '\n\n')
            elif sojourn_distribution_list[i][0] == 'poisson':
                f.write('POISSON INF_BOUND : 1 PARAMETER : '
                        + str(sojourn_distribution_list[i][1]) + '\n\n')
            elif sojourn_distribution_list[i][0] == 'negative binomial':
                f.write('NEGATIVE_BINOMIAL INF_BOUND : 1 PARAMETER : '
                        + str(sojourn_distribution_list[i][1])
                        + ' PROBABILITY : ' + str(sojourn_distribution_list[i][2]) + '\n\n')
        f.write('\n')
        f.write('1 OUTPUT_PROCESS\n\n')
        f.write('OUTPUT_PROCESS 1 : CATEGORICAL\n\n')
        for i in xrange(0, k):
            f.write('STATE ' + str(i) + ' OBSERVATION_DISTRIBUTION\n')
            for j in xrange(0, l):
                f.write('OUTPUT ' + str(j) + ' : ' + str(B[i, j]) + '\n')
            f.write('\n')
        return hsmc_file


size_of_hidden_state_space = 4
size_of_observation_state_space = 5

number_of_simulations = 5

m1_criterion_list = []
m2_criterion_list = []

for i in xrange(0, number_of_simulations):
    pi = dirichlet(np.ones(size_of_hidden_state_space), 1)
    A = generate_random_transition_matrix(size_of_hidden_state_space)
    B = dirichlet(np.ones(size_of_observation_state_space), size_of_hidden_state_space)
    sojourn_distribution_list = generate_random_sojourn_distributions(size_of_hidden_state_space)
    input_file = write_input_file(pi[0], A, B, sojourn_distribution_list)
    print input_file
    hsmm = HiddenSemiMarkov(input_file)
    simulated_data = hsmm.simulation_markovian_sequences(100, EyeMovementData().get_input_sequence, False)
    text_reading_list = []
    for text_reading in simulated_data:
        fixation_list = []
        for fixation in text_reading:
            fixation_list.append(fixation[0])
        text_reading_list.append(fixation_list)
    simulated_data = text_reading_list
    m1 = Model(SequentialData(simulated_data), k=4)
    m1.iterate_em(100)
    m1_criterion_list.append(m1.criterion)
    m2 = Model(SequentialData(simulated_data), random_init=True, k=4, n_init=3, n_random_seq=30, n_iter_init=50)
    m2.iterate_em(50)
    m2_criterion_list.append(m2.criterion)