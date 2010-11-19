#!/usr/bin/env python
#-*- coding: utf-8 -*-
"""

.. topic:: distance_matrix.py summary

    A module dedicated to distance matrix plotting

    :Code status: mature
    :Documentation status: to be completed
    :Author: Thomas Cokelaer <Thomas.Cokelaer@sophia.inria.fr>

    :Revision: $Id: sequences.py 9311 2010-07-26 16:39:44Z cokelaer $

"""
__version__ = "$Id: sequences.py 9311 2010-07-26 16:39:44Z cokelaer $"

import numpy

__all__ = ['ImprovedDistanceMatrix']

class ImprovedDistanceMatrix(object):
    """temporary class to manipulate and plot Distance Matrix

    A distance matrix is generated by functions `Compare`, which is
    not associated to specific viewpoint except an average distance over row
    or column. 

    This class provides more complex plotting routines such as histogram or images.

    
    .. plot::
        :include-source:
        :width: 440px
        :height: 360px
        
        from openalea.stat_tool import *
        from openalea.sequence_analysis import *
        from os.path import join as pj
        seq = Sequences(get_shared_data('chene_sessile_15pa.seq'))
        distance_matrix = Compare(seq, VectorDistance("O","O","O","O", "O", "O"))
        Plot(distance_matrix)
        
        
        

    """

    def __init__(self, dm):
        """
        :param dm: a distance matrix


        """
        self.distance_matrix = dm
        self.distance = numpy.array(numpy.zeros((dm.nb_row,dm.nb_column)))
        for i in range(0,dm.nb_row):
            for j in range(0,dm.nb_column):
                self.distance[i,j] = dm.get_distance(i,j)
        self.length = numpy.array(numpy.zeros((dm.nb_row,dm.nb_column)))
        for i in range(0,dm.nb_row):
            for j in range(0,dm.nb_column):
                self.length[i,j] = dm.get_length(i,j)

    def show_distance(self, fig=1, norm=True, **kargs):
        """Plot the distance values between two sequences for all
        the sequences. Normalizes distance by the length.
        
        :param norm: True by default. If not, distances are not divided by the length
        
        .. plot::
            :include-source:
        
            from openalea.stat_tool import *
            from openalea.sequence_analysis import *
            from os.path import join as pj
            seq = Sequences(get_shared_data('chene_sessile_15pa.seq'))
            distance_matrix = Compare(seq, VectorDistance("O","O","O","O", "O", "O"))
            imp_dm = ImprovedDistanceMatrix(distance_matrix)
            imp_dm.show_distance()
        """
        import pylab
        pylab.figure(fig)
        pylab.clf()
        if norm==True:
            pylab.imshow(self.distance/self.length)
        else:
            pylab.imshow(self.distance)
        pylab.xlabel("Sequence #")
        pylab.ylabel("Sequence #")
        pylab.title("Distance Between Sequences #")
        pylab.colorbar()
        pylab.show()

    def hist(self, minv=None, maxv=None, bins=10, norm=True):
        """Plot histogram of normalised distances
        
        :param norm: True by default. If not, distances are not divided by the length
        :param bin: number of histogram bins (10)
        :param minv: remove values below minv (None)
        :param maxv: remove values above maxv (None)
        
        .. plot::
            :include-source:
        
        
            from openalea.stat_tool import *
            from openalea.sequence_analysis import *
            from os.path import join as pj
            seq = Sequences(get_shared_data_path('chene_sessile_15pa.seq'))
            distance_matrix = Compare(seq, VectorDistance("O","O","O","O", "O", "O"))
            imp_dm = ImprovedDistanceMatrix(distance_matrix)
            imp_dm.hist()
        """
        import pylab
        import matplotlib.mlab as mlab
        
        index = numpy.where(self.length>0)
        data = self.distance[index].flatten()
        length = self.length[index].flatten()
        
        if norm==True:
            data = data/length

        #length is not required anymore
        del length
        if minv!=None:
            data = data[numpy.where(data>minv)[0]]
            
        if maxv!=None:
            data = data[numpy.where(data<maxv)[0]]

        n, bins, patches = pylab.hist(data, bins, normed=1, facecolor='green', alpha=0.75)

        # add a 'best fit' line
        mu = numpy.mean(data)
        sigma = numpy.std(data)
        y = mlab.normpdf( bins, mu, sigma)
        l = pylab.plot(bins, y, 'r--', linewidth=1)

        pylab.show()


    def cumulated_distance(self, verbose=False):
        c = numpy.zeros(self.distance_matrix.nb_row)

        for i in range(0,self.distance_matrix.nb_row):
            index = numpy.where(self.length[i,:]>0)
            c[i] = numpy.mean(self.distance[i,index]/self.length[i,index])

        i = numpy.argsort(c)
        if verbose==True:
            for x in i:
                print x+1, c[x]
        return c

