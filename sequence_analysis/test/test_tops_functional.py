from openalea.sequence_analysis import *

param1 = TopParameters(0.6, 0.6, 1.2)
param1 = TopParameters(get_shared_data("test_param1.p"), MaxPosition=20)
Plot(param1)
Display(param1)

top1 = Simulate(param1, 200, 30)
Plot(top1)
Display(top1)

top2 = Simulate(param1, 200, 20)

top5 = Merge(top1,  top2)

param2 = Estimate(top1)
Plot(param2)
Display(param2)

param4 = Estimate(top1, MinPosition=1, MaxPosition=10)
param5 = Estimate(top1, MinPosition=11, MaxPosition=20)
