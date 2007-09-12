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
Main interface
"""

__license__= "GPL2.0"
__revision__=" $Id: sceneobject.py 559 2007-05-25 12:25:30Z dufourko $ "

from error import *
from distribution import *


# API

def Histogram(arg):
    """
    Histogram:
    =========
    Construction of a frequency distribution from an object of type array(int) or from an ASCII file.
    Usage
    -----
    Histogram(array)

    Histogram(file_name)
    
    Arguments
    ---------
    array (array(int)),

    file_name (string).
    
    Returned Object
    ---------------
    If the construction succeeds, an object of type HISTOGRAM is returned, otherwise no object is returned.
    
    Description
    -----------
    In the file syntax, the frequencies fi for each possible value i are given in two columns. In the case of an argument of type (array(int)), it is simply assumed that each array element represents one data item. 

    See Also
    --------
    Save, Cluster, Merge, Shift, Transcode, ValueSelect, Compare (distributions), Estimate (distributions).
    """
    return DistributionData(arg)





def Display(obj, ViewPoint=None, Detail=1, Format="Column"):
    """
Display:
=======
    
    ASCII output of an object of the STAT module
    (Print: ASCII printing of an object of the STAT module.
    Print takes the same arguments as Display).

Usage
-----
    Display(obj, Detail=2)

    Display(vec, ViewPoint="Data", Detail=2)
    Display(seq, ViewPoint="Data", Format="Line", Detail=2)

    Display(dist, ViewPoint="Survival")
    Display(histo, ViewPoint="Survival")

    Display(hmc, identifier, ViewPoint="StateProfile")
    Display(hsmc, identifier, ViewPoint="StateProfile")
    
Arguments
---------
    obj: STAT module object,

    vec (vectors),
    seq (sequences, discrete_sequences, markov_data, semi-markov_data, tops),

    dist (distribution, mixture, convolution, compound),
    histo (histogram, mixture_data, convolution_data, compound_data),

    hmc (hidden_markov),
    hsmc (hidden_semi-markov),
    identifier (int): identifier of a sequence. 

Optional Arguments
------------------
    * ViewPoint (string): point of view on the object ("Survival" or "Data" or "StateProfile").
     This optional argument can be set at
          -> "Data" only if the first argument is of type vectors, sequences, discrete_sequences, markov_data, semi-markov_data or tops,
          -> "Survival" only if the first argument is of type distribution, mixture, convolution, compound, histogram, mixture_data, convolution_data or compound_data
          -> "StateProfile" only if the first argument is of type hidden_markov or hidden_semi-markov.

    * Detail (int): level of detail: 1 (default value) or 2. This optional argument cannot be used if the optional argument ViewPoint is set at "Survival" or "StateProfile".
       
    * Format (string): format of sequences (only relevant for multivariate sequences): "Column" (default value) or "Line". This optional argument can only be used if the optional argument ViewPoint is set at "Data", and hence, if the first argument is of type sequences, discrete_sequences, markov_data, semi-markov_data or tops.
       
Description
-----------
    ASCII output of sets of sequences or tops (ViewPoint="Data"): the format "Column" corresponds to the ASCII file syntax for objects of type sequences or tops. For a given value of the index parameter, the different variables are successively displayed. With the format "Line", the univariate sequence for each variable are displayed on consecutive lines. In the case of univariate sequences, the two formats give the same output.

    ASCII output of a (frequency) distribution and the associate hazard or survival rates (ViewPoint="Survival"): It is assumed that the (frequency) distribution represents lifetime and the hazard or survival rates are deduced from this lifetime distribution.

    ASCII output of the state profile given by the smoothed probabilities as a function of the index parameter t computed from the parameters of a hidden Markovian model for the sequence (ViewPoint="StateProfile").
    
Returned Object
---------------
    No object returned.
    
See Also
--------
    Plot, Save.
"""

    obj.display(ViewPoint, Detail, Format)



def Save(obj, filename, ViewPoint="Data", Detail=1, Format="ASCII"):
    """
Save
====
    Saving of an object of the STAT module in a file.
    
Usage
-----
    Save(obj, file_name, Format="ASCII", Detail=2)

    Save(histo, file_name, ViewPoint="Data")
    Save(vec, file_name, ViewPoint="Data", Detail=2)
    Save(timev, file_name, ViewPoint="Data")
    Save(seq, file_name, ViewPoint="Data", Format="Line",
    Detail=2)

    Save(dist, file_name, ViewPoint="Survival",
    Format="SpreadSheet")
    Save(histo, file_name, ViewPoint="Survival",
    Format="SpreadSheet")

    Save(hmc, ViewPoint="StateProfile", Sequence=1,
    Format="SpreadSheet")
    Save(hsmc, ViewPoint="StateProfile", Sequence=1,
    Format="SpreadSheet")
    
Arguments
---------
    obj: object of the STAT module (except objects of type vector_distance),
    file_name (string),

    histo (histogram, mixture_data, convolution_data, compound_data),
    vec (vectors),
    timev (time_events, renewal_data),
    seq (sequences, discrete_sequences, markov_data, semi-markov_data, tops).

    dist (distribution, mixture, convolution, compound),

    hmc (hidden_markov),
    hsmc (hidden_semi-markov).
    
Optional Arguments
------------------
    ViewPoint (string): point of view on the object ("Data" or "Survival" or "StateProfile"). This optional argument can be set at "Data" only if the first argument is of type vectors, sequences, discrete_sequences, markov_data, semi-markov_data or tops, can be set at "Survival" only if the first argument is of type distribution, mixture, convolution, compound, histogram, mixture_data, convolution_data or compound_data and can be set at "StateProfile" only if the first argument is of type hidden_markov or hidden_semi-markov.
    Detail (int): level of detail: 1 (default value) or 2. This optional argument can only be used if the optional argument ViewPoint is not set, or if the optional argument ViewPoint is set at "Data" and if the first mandatory argument is of type vectors, sequences, discrete_sequences, markov_data, semi-markov_data or tops.
    Format (string): file format: "ASCII" (default format), "Binary" or "SpreadSheet". These file formats cannot be specified if the optional argument ViewPoint is set at "Data". The optional argument Format can only be set at "Binary" if the optional argument ViewPoint is not set.
    Format (string): format of sequences (only relevant for multivariate sequences): "Column" (the default) or "Line". This optional argument can only be used if the optional argument ViewPoint is set at "Data", and hence, if the first mandatory argument is of type sequences, discrete_sequences, markov_data, semi-markov_data or tops.
    Sequence (int): identifier of a sequence. This optional argument can only be used if the optional argument ViewPoint is set at "StateProfile", and hence, if the first mandatory argument is of type hidden_markov or hidden_semi-markov. 

Description
-----------
    Saving of sets of sequences or 'tops' (ViewPoint="Data"): the format "Column" corresponds to the ASCII file syntax for objects of type sequences or tops. For a given value of the index parameter, the different variables are successively written. With the format "Line", the univariate sequence for each variable are written on consecutive lines. In the case of univariate sequences, the two formats give the same file.

    Saving of a (frequency) distribution and the associate hazard or survival rates (ViewPoint="Survival"): It is assumed that the (frequency) distribution represents lifetime and the hazard or survival rates are deduced from this lifetime distribution.
    Saving of the state profile given by the smoothed probabilities as a function of the index parameter t computed from the parameters of a hidden Markovian model for the sequence (ViewPoint="StateProfile"). 
Background
    The persistence mechanism is implemented by the Save function with the optional argument Format set at "Binary" for saving and by the function Load for restoration.
    
Returned Object
---------------
    No object returned.
    
See Also
--------
    Display, Plot, Compound, Convolution, Distribution, HiddenMarkov, HiddenSemiMarkov, Histogram, Markov, Mixture, Renewal, SemiMarkov, Sequences, TimeEvents, Tops, TopParameters
"""

    obj.save(filename, ViewPoint, Detail, Format)


def Plot(obj, **args):
    """
Plot
====
    Graphical output of an object of the STAT module using the GNUPLOT software.
    
Usage
-----
    Plot(obj1, Title="Distribution")
    w1 = NewPlot(obj1, Title="Distribution")
    Plot(Window=w1)
    w1 = NewPlot()

    Plot(vec1, Title="Values")
    Plot()
    Plot(vecn, variable, Title="Vectors")
    Plot(variable)

    Plot(obj2, type, Title="Sequences")
    Plot(type)
    Plot(obj3, type, variable, Title="Multivariate sequences")
    Plot(type, variable)

    Plot(dist1, dist2,..., Title="Family of distributions")
    Plot(histo1, histo2,...,
    Title="Family of frequency distributions")

    Plot(seq, ViewPoint="Data")

    Plot(dist, ViewPoint="Survival", Title="Survival rates")
    Plot(histo, ViewPoint="Survival", Title="Survival rates")

    Plot(hmc, identifier, ViewPoint="StateProfile",
    Title="Smoothed probabilities")
    Plot(hsmc, identifier, ViewPoint="StateProfile",
    Title="Smoothed probabilities")
    
Arguments
---------
    obj1 (distribution, mixture, convolution, compound, histogram, mixture_data, convolution_data, compound_data, renewal, time_events, renewal_data, sequences, distance_matrix, top_parameters, tops),
    vec1 (vectors): values,
    vecn (vectors): vectors,
    variable (int): variable index,

    obj2 (markov, semi-markov, hidden_markov, hidden_semi-markov, discrete_sequences, markov_data, semi-markov_data): Markovian model for discrete univariate sequences or discrete univariate sequences,
    obj3: (markov, semi-markov, hidden_markov, hidden_semi-markov, discrete_sequences, markov_data, semi-markov_data): Markovian model for discrete multivariate sequences or discrete multivariate sequences,
    type (string): type of graphical outputs in the case of Markovian models or sequences: "SelfTransition", "Observation", "Intensity", "FirstOccurrence", "Recurrence", "Sojourn" or "Counting",

    dist1, dist2, ... (distribution, mixture, convolution, compound),
    histo1, histo2, ... (histogram, mixture_data, convolution_data, compound_data),

    seq (sequences, discrete_sequences, markov_data, semi-markov_data, tops),

    dist (distribution, mixture, convolution, compound),
    histo (histogram, mixture_data, convolution_data, compound_data),

    hmc (hidden_markov),
    hsmc (hidden_semi-markov),
    identifier (int): identifier of a sequence.
    
Optional Arguments
------------------
    Window (window): window (the default: last used window). This optional argument cannot be used with the function NewPlot.
    ViewPoint (string): point of view on the object ("Data" or "Survival" or "StateProfile"). This optional argument can be set at "Data" only if the first mandatory argument is of type sequences, discrete_sequences, markov_data, semi-markov_data or tops, can be set at "Survival" only if the first mandatory argument is of type distribution, mixture, convolution, compound, histogram, mixture_data, convolution_data or compound_data and can be set at "StateProfile" only if the first mandatory argument is of type hidden_markov or hidden_semi-markov.
    Title (string): graphic title (the default: no title). 
Returned Object
    No object returned by the function Plot. An object of type window is returned by the function NewPlot.
    
Description
-----------
    In the case of Markovian models or sequences, the graphical outputs are grouped as follows:
    "SelfTransition": self-transition probability as a function of the index parameter (non-homogeneous Markov chain),
    "Observation": observation distributions attached to each state of the underlying (semi-)Markov chain (lumped processes or hidden Markovian processes),
    "Intensity": (empirical) probabilities of states/outputs as a function of the index parameter,
    "FirstOccurrence": (frequency) distributions of the time-up to the first occurrence of a state/output (or first-passage time in a state/output distributions),
    "Recurrence" (frequency) distributions of the recurrence time in a state/output,
    "Sojourn": (frequency) distributions of the sojourn time in a state/output (or state/output occupancy distributions). For the frequency distributions extracted from sequences, the sojourn times in the last visited states which are considered as censored are isolated.
    "Counting": counting (frequency) distributions (either distributions of the number of runs (or clumps) of a state/output per sequence or distributions of the number of occurrences of a state/output per sequence).
Background
    Graphical output of a (frequency) distribution and the associate hazard or survival rates (ViewPoint="Survival"): It is assumed that the (frequency) distribution represents lifetime and the hazard or survival rates are deduced from this lifetime distribution.
    Graphical output of the state profile given by the smoothed probabilities as a function of the index parameter t computed from the parameters of a hidden Markovian model for the sequence (ViewPoint="StateProfile"). 

See Also
--------
    Display, Save.
"""

    obj.plot(**args)
