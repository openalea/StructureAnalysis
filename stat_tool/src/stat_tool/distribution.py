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

__doc__="""
Histogram classes
"""

__license__= "GPL2.0"
__revision__=" $Id: sceneobject.py 559 2007-05-25 12:25:30Z dufourko $ "



from _stat_tool import _Format_error
from _stat_tool import _Distribution_data, _histogram_ascii_read

from error import StatToolError



class Histogram(object):
    """ Histogram class """

    def __init__(self, arg):
        """
        Constructor
          @param arg : call correct constructor depending of the type
        """

        # Factory map
        fmap = { str : self.from_file,
                 list : self.from_list,
                 }

        try:
            f = fmap[type(arg)]
            f(arg)
            
        except KeyError:
            raise StatToolError("Bad Argument %s", arg)


    def from_file(self, filename):
        """ Build histogram from file """

        ferror = _Format_error()
        self.__wrapped = _histogram_ascii_read(ferror, filename)

        if(not self.__wrapped):
            raise StatToolError("Bad Argument %s", filename)
            


    def from_list(self, array):
        """ Build histogram from a list of int """

        self.__wrapped = _Distribution_data(array)


    def plot(self):
        """ Plot function """
        pass


    def display(self):
        """ Display function """
        pass


    def save(self, filename):
        """ File output """
        pass
        
    

    
