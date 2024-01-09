"""Testing a user case."""
from ema import EyeMovementData
from ema import Indicator
from ema import Model

# import data
data = EyeMovementData()

# model with auto-initilization (non random)
model = Model(data, k=5)

# iterate EM for convergence
model.iterate_em(1000)

# plot indicators
fdur = Indicator(model, 'FDUR')
fdur.boxplot()

sacamp = Indicator(model, 'SACAMP')
sacamp.boxplot()

winc = Indicator(model, 'WINC')
winc.boxplot()

cosinst = Indicator(model, 'COSINST')
cosinst.boxplot()
