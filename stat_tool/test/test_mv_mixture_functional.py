"""tests on mv_mixture"""
__version__ = "$Id$"


from openalea.stat_tool.mvmixture import _MvMixture
from openalea.stat_tool.distribution import Binomial, Poisson

from openalea.stat_tool.output import Plot


def test1():


    #ORIGINAL AML stops HERE 
     
    d11 = Binomial(0, 12, 0.1)
    d12 = Binomial(2, 13, 0.6)
    d13 = Binomial(3, 15, 0.9)
    
    d21 = Poisson(0, 25.0)
    d22 = Poisson(0, 5.0)
    d23 = Poisson(0, 0.2)

    m = _MvMixture([0.1, 0.2, 0.7], [[d11, d21], [d12, d22], [d13, d23]])
    print m

    #m2 = _MvMixture("mixture_mv1.mixt")
    #print m2

    #print "Egalite des melanges construits par liste ",\
    #  "de distributions et par fichiers : ", str(str(m)==str(m2))

    #m = _MvMixture("mixture_mv_nonparam.mixt")
   # print m

    print "Simulation de melanges multivaries : "
    v = m.simulate(5000)
    print v

    Plot(m, variable=1, Title="Simulated mixture")

    print "Estimation de melanges multivaries ", \
    #    "d'apres un modele initial : "
    m_estim_model = v.mixture_estimation(m, 100,  [True, True])
      
    Plot(m_estim_model, variable = 1, Title="Estimated mixture")
    
    print "Estimation de melanges multivaries ", \
        "d'apres un nombre de composantes : "
        
    m_estim_nbcomp = v.mixture_estimation(3, 100, [True, True])
    
    m_estim_nbcomp.plot(variable = 1, Title="Estimated mixture")


if __name__ == "__main__":
    test1()
