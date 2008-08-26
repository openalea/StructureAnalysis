# -*- python -*-
#
#       vplants.stat_tool
#
#       Copyright 2008 INRIA - CIRAD
#
#       File author(s): 
#           Samuel Dufour-Kowalski <samuel.dufour@sophia.inria.fr>
#           Yann Guedon <yann.guedon@cirad.fr>
#
#       Distributed under the GPL 2.0 License.
#       See accompanying file LICENSE.txt or copy at
#           http://www.gnu.org/licenses/gpl-2.0.txt
# 

__doc__="""
Stat_Tool
"""

# Test prefix
test_prefix = "../../../test/"

def get_test_prefix():
    global test_prefix
    return test_prefix


def set_test_prefix(prefix=""):
    global test_prefix
    test_prefix = prefix


from distribution import *
from histogram import *
from mixture import *
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

# Exception
from error import StatToolError
