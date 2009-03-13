# -*- python -*-
#
#       vplants.stat_tool
#
#       Copyright 2006-2007 INRIA - CIRAD - INRA  
#
#       File author(s): Samuel Dufour-Kowalski <samuel.dufour@sophia.inria.fr>
#
#       Distributed under the GPL 2.0 License.
#       See accompanying file LICENSE.txt or copy at
#           http://www.gnu.org/licenses/gpl-2.0.txt
# 
#       OpenAlea WebSite : http://openalea.gforge.inria.fr
#
"""Error classes"""

__license__= "GPL2.0"
__revision__=" $Id$ "


class StatToolError(Exception):
    def __init__(self, arg):
        Exception.__init__(self, str(arg))


