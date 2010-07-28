"""Sequence Analysis init file"""
__revision__ = "$Id$"

from openalea.stat_tool import *
from openalea.stat_tool._stat_tool import *
#must be after from openalea_stat_tool import *
from openalea.deploy.shared_data import get_shared_data_path
shared_data_path = get_shared_data_path(__path__)

import _sequence_analysis

import openalea.stat_tool.interface as interface


from correlation import *
from simulate import *
from compare import *


from time_events import *
from top_parameters import *
from tops import *
from sequences import *
from hidden_semi_markov import *
from hidden_variable_order_markov import *
from semi_markov import *
from data_transform  import *

from estimate import *
from nonhomogeneous_markov import *
from renewal import *
from variable_order_markov import *

# must be after from openaeal.stat_tool import *
from openalea.deploy.metainfo import read_metainfo
metadata = read_metainfo(__path__[0] + '/../../../metainfo.ini', verbose=False)
__version__ = metadata['version']
