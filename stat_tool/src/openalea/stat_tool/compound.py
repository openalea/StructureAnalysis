__doc__ = """ Compound """
__docformat__ = "restructuredtext"

import interface
import _stat_tool


#from _stat_tool import _Compound
#from _stat_tool import _CompoundData

# __all__ = ['Compound',
#            '_Compound',
#            '_CompoundData',
#            ]



def Compound(*args):
    """
    Construction of a compound of distributions from a sum distribution and an elementary distribution 
    or from an ASCII file.
    
    Usage:
      Compound(sum_dist, dist)
      Compound(filename)

    Arguments:
      sum_dist, dist (DISTRIBUTION, MIXTURE, CONVOLUTION, COMPOUND)
      filename(STRING)

    Description:

        """

    if((len(args)==0) or (len(args)>2)) : 
        raise TypeError("Bad number of arguments")

    # filename
    if(len(args)==1) :
        return _stat_tool._Compound(args[0])

    # build list of distributions
    if(len(args)==2) :
        return _stat_tool._Compound(args[0], args[1])


# Extend _Compound
#interface.extend_class( _stat_tool._Compound, interface.StatInterface)


# Extend _CompoundData
# interface.extend_class( _stat_tool._CompoundData, interface.StatInterface)



########################## Test Compound ########################################
from openalea.stat_tool import get_test_file

#class Test:
 #    def test_emty(self):
#         try:
#             m = Compound()
#             assert False

#         except TypeError:
#             assert True

#     def test_file(self):
#         c = Compound(get_test_file("compound1.cd"))
#         assert c


#     def test_build_compound(self):
#         from distribution import Uniform

#         d1 = Uniform(0,10)
#         d2 = Uniform(10,20)

#         m = Compound(d1, d2)
#         assert m
#         return m


#     def test_plot(self):

#         m = test_build_compound()
#         m.plot()

#         assert str(m)
#         m.display()


#     def test_simulation(self):

#         m = test_build_compound()
#         s = m.simulate(1000)

#         assert len(s) == 1000
#         assert str(s)


#     def test_extract(self):
#         from extract import ExtractDistribution
#         from distribution import Uniform

#         m = test_build_compound()

#         assert m.extract_compound() == ExtractDistribution(m, "Compound")

#         assert m.extract_sum() == Uniform(0,10)
#         assert m.extract_sum() == ExtractDistribution(m, "Sum")

#         assert m.extract_elementary() == Uniform(10,20)
#         assert m.extract_elementary() == ExtractDistribution(m, "Elementary")


#     def test_extract_data(self):

#         assert False    
