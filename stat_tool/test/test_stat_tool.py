# import openalea.stat_tool as stat_tool

# inf_bound=0
# sup_bound=3
# probability= 0.6
# sample_size= 10
# sample= []
# val= [0]
# ident=stat_tool.DistributionIdentifier.UNIFORM
# parameter=stat_tool.D_DEFAULT
# distrib= stat_tool.Parametric(ident, inf_bound, sup_bound, 
#                                parameter, probability)
# print "A uniform distribution:"
# print distrib
# print "Simulation from this distribution: "
# print distrib.simulation()
# for i in range(sample_size):
#     val[0]= distrib.simulation()
#     sample= sample+val
# print sample
# print "Histogram build from this distribution: "
# distrib_data=stat_tool.Histogram(sample)
# print distrib_data
# distrib_data.Display(Detail=2)
# # distrib_data.Plot()
