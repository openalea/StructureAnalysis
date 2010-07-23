HIDDEN_MARKOV_CHAIN

3 STATES

INITIAL_PROBABILITIES
0.8       0.1       0.1       

TRANSITION_PROBABILITIES      # memory
# 0.352885  0.547115  0.1            0   non-terminal
# 0.381192  0.423546  0.195262     0 0   non-terminal
0.496324  0.403676  0.1        0 0 0   #     terminal
0.367958  0.397887  0.234155   1 0 0   #     terminal
0.1       0.459574  0.440426   2 0 0   #     terminal
0.256309  0.643691  0.1          1 0   #     terminal
0.680294  0.219706  0.1          2 0   #     terminal
# 0.7594    0.1406    0.1            1   non-terminal
0.782464  0.117536  0.1          0 1   #     terminal
0.625846  0.274154  0.1          1 1   #     terminal
0.72      0.18      0.1          2 1   #     terminal
0.8       0.1       0.1            2   #     terminal

1 OUTPUT_PROCESS

OUTPUT_PROCESS 1 : NONPARAMETRIC

STATE 0 OBSERVATION_DISTRIBUTION
OUTPUT 0 : 0.4
OUTPUT 1 : 0.3
OUTPUT 2 : 0.3

# OUTPUT 0 : 1.0

STATE 1 OBSERVATION_DISTRIBUTION
OUTPUT 0 : 0.3
OUTPUT 1 : 0.4
OUTPUT 2 : 0.3

# OUTPUT 0 : 0.4
# OUTPUT 1 : 0.6

STATE 2 OBSERVATION_DISTRIBUTION
OUTPUT 0 : 0.3
OUTPUT 1 : 0.3
OUTPUT 2 : 0.4
