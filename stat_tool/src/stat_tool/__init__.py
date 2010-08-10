# -*- coding: utf-8 -*-
# -*- python -*-
#
#       vplants.stat_tool
#
#       Copyright 2008 INRIA - CIRAD
#
#       File author(s):
#           Samuel Dufour-Kowalski <samuel.dufour@sophia.inria.fr>
#           Yann Guedon <yann.guedon@cirad.fr>
#           Thomas Cokelaer <Thomas.Cokelaer@sophia.inria.fr>
#
#       Distributed under the GPL 2.0 License.
#       See accompanying file LICENSE.txt or copy at
#           http://www.gnu.org/licenses/gpl-2.0.txt
#
"""Stat_Tool init file"""
__revision__ = "$Id$"


from openalea.deploy.shared_data import get_shared_data_path
from os.path import join as pj

def get_shared_data(file):
    import openalea.stat_tool
    shared_data_path = get_shared_data_path(openalea.stat_tool.__path__)
    return pj(shared_data_path, file)

import _stat_tool
from distribution import *
from histogram import *
from mixture import *
from multivariate_mixture import *
from compound import *
from convolution import *
from vectors import *

from estimate import *
from simulate import *
from comparison import *
from regression import *

from output import  *
from data_transform import *
from cluster import *

# Constant
from _stat_tool import I_DEFAULT
from _stat_tool import D_DEFAULT
from _stat_tool import D_INF
from _stat_tool import MAX_DIFF_BOUND
from _stat_tool import MAX_MEAN
from _stat_tool import VariableType
from _stat_tool import VariableTypeBis
from _stat_tool import RestorationAlgorithm

# Error
from _stat_tool import StatError
