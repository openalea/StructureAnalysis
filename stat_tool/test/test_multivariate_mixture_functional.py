# -*- coding: utf-8 -*-
"""tests on mv_mixture"""
__version__ = "$Id$"


from openalea.stat_tool.multivariate_mixture import _MultivariateMixture
from openalea.stat_tool.distribution import Binomial, Poisson
from openalea.stat_tool.distribution import set_seed

from openalea.stat_tool.output import plot, Plot

from openalea.stat_tool.plot import get_plotter, mplotlib 

def test1():

    set_seed(4)
    plotter = mplotlib()
     
    d11 = Binomial(0, 12, 0.1)
    d12 = Binomial(0, 12, 0.6)
    d13 = Binomial(0, 12, 0.9)
    
    d21 = Poisson(0, 25.0)
    d22 = Poisson(0, 5.0)
    d23 = Poisson(0, 0.2)

    m = _MultivariateMixture([0.1, 0.2, 0.7], [[d11, d21], [d12, d22], [d13, d23]])
    print(m)

    Plot(m, Title="Mixture model used for simulation: ")

    print("Simulate multivariate mixture: ")
    v = m.simulate(500000)
    print(v)

    # TODO: matplotlib output for ._MultivariateMixtureData
    assert(v.get_plotable())
    Plot(v, Title="Simulated mixture: ")
    

    print("Estimate multivariate mixture")
    #    "d'apres un modele initial : "
    m_estim_model = v.mixture_estimation(m, 100,  [True, True])

    extracted_mixture = m_estim_model.extract_mixture(1)
    extracted_mixture.old_plot(variable=1, Title="Marginal distribution")
    # TODO: check why clustering is so bad (or histograms are so far from distributions)
    Plot(m_estim_model, Title="Estimated mixture")
    print(m_estim_model.display(True))

    print("Estimation de melanges multivaries ", 
        "d'apres un nombre de composantes : ")
        
    m_estim_nbcomp = v.mixture_estimation(2, 100, [True, True])
    
    m_estim_nbcomp.plot(variable = 1, Title="Estimated mixture")

    clust_entropy = m_estim_nbcomp.cluster_data(v , True)
    clust_plain = m_estim_nbcomp.cluster_data(v , False)


if __name__ == "__main__":
    test1()
