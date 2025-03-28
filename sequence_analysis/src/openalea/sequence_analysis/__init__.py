"""Sequence Analysis init file"""
__revision__ = "$Id$"

from openalea.stat_tool import *
#import openalea.stat_tool._stat_tool

from openalea.deploy.shared_data import get_shared_data_path
from os.path import join as pj

def get_shared_data(file):
    import openalea.sequence_analysis
    shared_data_path = get_shared_data_path(openalea.sequence_analysis.__path__)
    return pj(shared_data_path, file)


import openalea.stat_tool.interface as interface


from .correlation import *
from .simulate import *
from .compare import *

from .time_events import *
# from top_parameters import *
# from tops import *
from .sequences import *
from .hidden_semi_markov import *
from .hidden_variable_order_markov import *
from .semi_markov import *
from .data_transform  import *

from .estimate import *
from .nonhomogeneous_markov import *
from .renewal import *
from .variable_order_markov import *
from .distance_matrix import *
from .enums_seq import *
