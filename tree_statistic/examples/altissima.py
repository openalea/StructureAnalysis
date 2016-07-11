import matplotlib.pyplot as plt
import openalea.stat_tool.vectors as vectors
import openalea.tree_statistic.dgdistributions as dgd

# Work on Altissima data
Altissima = vectors.Vectors('./Data/altissima.vec')
l = list(Altissima)
rl = list()
rlt = list()
for i in l:
	rlt = list()
	for j in range(0, len(i)):
		if(i[j] > 0):
			i[j] = 1
	for j in range(0,15):
		if(i[j] > 0):
			rlt.append(1)
		else:
			rlt.append(0)
	rl.append(rlt)



data = dgd.DiscreteGraphicalData(15, rl)
data.Display()

data.DistributionEstimation(dgd.Algorithms.ML, number_of_edges = 3)
data.ParametricEstimation(dgd.Algorithms.ML, number_of_edges = 3)

# Non-parametric estimation of DiscreteGraphicalDistributions (all distributions that an algorithm car inferred are saved !)
CLTs = data.DistributionEstimation(dgd.Algorithms.CLT, all = True)
CLRs = data.DistributionEstimation(dgd.Algorithms.CLR, all = True)
MRNETs = data.DistributionEstimation(dgd.Algorithms.MRNET, all = True)
Relevances = data.DistributionEstimation(dgd.Algorithms.Relevance, all = True)
MLs = data.DistributionEstimation(dgd.Algorithms.ML, all = True)

# Plot loglikelihoods (for each algorithm plot all the log-likelihood of inferred distributions)
data.LogLikelihood(CLTs, label = 'CLT algorithm')
data.LogLikelihood(CLRs, label = 'CLR algorithm')
data.LogLikelihood(MRNETs, label = 'MRNET algorithm')
data.LogLikelihood(Relevances, label = 'Relevance algorithm')
data.LogLikelihood(MLs, label = 'ML algorithm')
plt.legend()
plt.show()

# Plot BIC (for each algorithm plot all the BIC of inferred distributions)
data.PenalizedLogLikelihood(CLTs, dgd.Criterions.BIC, label = 'CLT algorithm')
data.PenalizedLogLikelihood(CLRs, dgd.Criterions.BIC, label = 'CLR algorithm')
data.PenalizedLogLikelihood(MRNETs, dgd.Criterions.BIC, label = 'MRNET algorithm')
data.PenalizedLogLikelihood(Relevances, dgd.Criterions.BIC, label = 'Relevance algorithm')
data.PenalizedLogLikelihood(MLs, dgd.Criterions.BIC, label = 'ML algorithm')
plt.legend()
plt.show()

# Plot AIC (for each algorithm plot all the AIC of inferred distributions)
data.PenalizedLogLikelihood(CLTs, dgd.Criterions.AIC, label = 'CLT algorithm')
data.PenalizedLogLikelihood(CLRs, dgd.Criterions.AIC, label = 'CLR algorithm')
data.PenalizedLogLikelihood(MRNETs, dgd.Criterions.AIC, label = 'MRNET algorithm')
data.PenalizedLogLikelihood(Relevances, dgd.Criterions.AIC, label = 'Relevance algorithm')
data.PenalizedLogLikelihood(MLs, dgd.Criterions.AIC, label = 'ML algorithm')
plt.legend()
plt.show()

# Get the best inferred DiscreteGraphicalDistribution according to BIC criterion
BIC_MLs = data.GetPenalizedLogLikelihood(MLs, dgd.Criterions.BIC)
best_ML = MLs[BIC_MLs.index(max(BIC_MLs))]
best_ML.Graph()
plt.show()
best_ML.Formula()
