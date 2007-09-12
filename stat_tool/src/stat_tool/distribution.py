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


class StatModel(object):
    """ Interface of Stat_Tool classes """

    def display(self, **args):
        raise NotImplementedError()


    def plot(self, **args):
        raise NotImplementedError()


    def save(self, filename, **args):
        raise NotImplementedError()



class DistributionData(StatModel):
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


    def __str__(self):
        return str(self.__wrapped)


    def from_file(self, filename):
        """ Build histogram from file """

        ferror = _Format_error()
        self.__wrapped = _histogram_ascii_read(ferror, filename)

        if(not self.__wrapped):
            raise StatToolError(ferror)


    def from_list(self, array):
        """
        Build histogram from a list of int
        Raise a TypeError Exception if list if badly constructed
        """

        self.__wrapped = _Distribution_data(array)

        if(not self.__wrapped):
            raise StatToolError("Bad Argument %s"%(str(array),))


    def plot(self, **args):
        """ Plot function """
        pass


    def display(self, **args):
        """ Display function """

        try:
            survival = args['ViewPoint'] is "Survival"
        except:
            survival = False

            
        if(survival):
            self.__wrapped.ascii
        


    def save(self, filename, **args):
        """ File output """
        pass
        


class ParametricModel(StatModel):
    """ Parametric model """
    pass

    
