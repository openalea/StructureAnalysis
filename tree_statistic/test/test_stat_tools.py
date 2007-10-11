import openalea.stat_tool
stat_tool=openalea.stat_tool

inf_bound=0
sup_bound=3
probability= 0.6
sample_size= 10
sample= []
val= [0]
ident=stat_tools.DistributionIdentifier.UNIFORM
parameter=stat_tools.D_DEFAULT
distrib= stat_tools.Parametric(ident, inf_bound, sup_bound, 
                               parameter, probability)
print "A uniform distribution:"
print distrib
print "Simulation from this distribution: "
print distrib.simulation()
for i in range(sample_size):
    val[0]= distrib.simulation()
    sample= sample+val
print sample
print "Histogram build from this distribution: "
distrib_data=stat_tools.Histogram(sample)
print distrib_data
distrib_data.Display(Detail=2)
# distrib_data.Plot()
