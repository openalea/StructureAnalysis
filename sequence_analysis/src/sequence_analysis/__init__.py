"""Sequence Analysis init file"""
__revision__ = "$Id: __init__.py 9407 2010-08-10 16:22:46Z cokelaer $"

from openalea.stat_tool import *
#import openalea.stat_tool._stat_tool

from openalea.deploy.shared_data import get_shared_data_path
from os.path import join as pj

def get_shared_data(file):
    import openalea.sequence_analysis
    shared_data_path = get_shared_data_path(openalea.sequence_analysis.__path__)
    return pj(shared_data_path, file)





import openalea.stat_tool.interface as interface


from correlation import *
from simulate import *
from compare import *


from time_events import *
from sequences import *
from hidden_semi_markov import *
from hidden_variable_order_markov import *
from semi_markov import *
from data_transform  import *

from estimate import *
from nonhomogeneous_markov import *
from renewal import *
from variable_order_markov import *
from distance_matrix import *
