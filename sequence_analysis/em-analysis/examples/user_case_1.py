"""Testing a user case."""
from ema import EyeMovementData
from ema import Model

# import data
data = EyeMovementData()

# initialize a model with random initializations for EM
model = Model(data, random_init=True, k=5, n_init=2, n_random_seq=40,
              n_iter_init=20)

# check out the criterion
print(model.criterion)

# iterate EM for convergence
model.iterate_em(2000)

# check out the criterion
print(model.criterion)

# manually check model parameters
model.print_hsmc_file(verbose=False)
# the 4th state seems to be an absorbing state according to the transition
# matrix and the occupancy distribution
model.parameters.transition_probabilities.transition_probabilities = (
    [0, 0, 0, 0, 1], 4)
# or equivalently
# model.parameters.occupancy_distributions[4].state_type = 'ABSORBING'

# iterate EM for new convergence
model.iterate_em(1000)

# check out the criterion
print(model.criterion)
