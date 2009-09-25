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
#
#       Distributed under the GPL 2.0 License.
#       See accompanying file LICENSE.txt or copy at
#           http://www.gnu.org/licenses/gpl-2.0.txt
# 
"""Stat_Tool"""

# Test prefix

def get_test_file(filename=""):
    """
    look for data in different directories
    """
 
    path1 = os.path.abspath("../../../test/")
    path2 = os.path.abspath("./test/")
    path3 = os.path.abspath("./")

    if os.path.isdir(path1):
        return os.path.join(path1, filename)
    elif os.path.isdir(path2):
        return os.path.join(path2, filename)
    else:
        return os.path.join(path3, filename)



import _stat_tool
from distribution import *
from histogram import *
from mixture import *
from mvmixture import *
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

from error import FormatError as FormatError

# Constant
from _stat_tool import I_DEFAULT
from _stat_tool import D_DEFAULT
from _stat_tool import D_INF
from _stat_tool import MAX_DIFF_BOUND
from _stat_tool import MAX_MEAN
from _stat_tool import VariableType
from _stat_tool import VariableTypeBis
from _stat_tool import RestorationAlgorithm

#
