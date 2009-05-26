"""Simulate

.. author:: Thomas Cokelaer, Thomas.Cokelaer@inria.fr 
"""
__revision__ = "$Id:  $"

import openalea.stat_tool._stat_tool as _stat_tool
import _sequence_analysis

from openalea.stat_tool.simulate import Simulate as SimulateDistribution
 
def Simulate(obj, *args, **kargs):
    """Simulate
    
    * Generation of a random sample from a distribution.
    * Generation of a random sample from a renewal process.
    * Generation of a sample of sequences from a (hidden) Markovian process.
    * Generation of 'tops' from 'top' parameters.

    
    :Usage:
    
    >>> Simulate(dist, size)
    >>> Simulate(mixt, size)
    >>> Simulate(convol, size)
    >>> Simulate(compound, size)
    
    >>> Simulate(renew, type, time_histo)
    >>> Simulate(renew, type, size, time)
    >>> Simulate(renew, type, size, timev)
    
    >>> Simulate(markov, length_histo)
    >>> Simulate(markov, size, length)
    >>> Simulate(markov, size, seq)
    >>> Simulate(semi_markov, length_histo, Counting=False)
    >>> Simulate(semi_markov, size, length, Counting=False)
    >>> Simulate(semi_markov, size, seq, Counting=False)
    
    >>> Simulate(top_param, size, nb_trial, NbAxillary=2)    

    :Arguments (distribution case):
    
    * dist (distribution),
    * mixt (mixture),
    * convol (convolution),
    * compound (compound),
    * size (int): sample size.
    
    :Arguments (renewal case):
    
    * renew (renewal),
    * type (string): type of renewal process: "Ordinary" or "Equilibriun",
    * time_histo (histogram, mixture_data, convolution_data, compound_data): frequency distribution of the length of the observation period,
    * size (int): sample size,
    * time (int): length of the observation period,
    * timev (time_events, renewal_data),
    
    :Arguments (markovian case):
    
    * markov (markov, hidden_markov),
    * semi_markov (semi-markov, hidden_semi-markov),
    * length_histo (histogram, mixture_data, convolution_data, compound_data): frequency distribution of sequence lengths,
    * size (int): sample size,
    * length (int): sequence length,
    * seq (discrete_sequences, markov_data, semi-markov_data),
    
    :Arguments (top case)

    * top_param (top_parameters),
    * size (int): sample size,
    * nb_trial (int): number of Bernoulli trials for the growth of the parent shoot,
    
    :Optional Arguments: 
    
    * Counting (bool): computation of counting distributions (default value: True).
    * NbAxillary (int): number of offspring shoots per node (default value: 1, should be < 4).

    :Returned Object:
    
    * If the first argument is of type distribution and if 0 < size < 1000000, an object of type HISTOGRAM is returned, otherwise no object is returned.
    * If the first argument is of type mixture and if 0 < size < 1000000, an object of type mixture_data is returned, otherwise no object is returned.
    * If the first argument is of type convolution and if 0 < size < 1000000, an object of type convolution_data is returned, otherwise no object is returned.
    * If the first argument is of type compound and if 0 < size < 1000000, an object of type compound_data is returned, otherwise no object is returned.
    * The returned object of type HISTOGRAM, mixture_data, convolution_data or compound_data contains both the simulated sample and the model used for simulation.
    * If 0 < (sample size) < 1000000, if (minimum length of the observation period) > 2 and if (maximum length of the observation period) XX 1000, an object of type renewal_data is returned, otherwise no object is returned. The returned object contains both the simulated sample (not only count data but also time sequences) and the model used for simulation.
    * If the first argument is of type markov or hidden_markov, if 0 < (sample size)
      < 100000, if (minimum sequence length) < 2, if (maximum sequence length) < 
      1000 and if (cumulative sequence length) < 1000000, an object of type 
      markov_data is returned, otherwise no object is returned.
    * If the first argument is of type semi-markov or hidden_semi-markov, if 
      0 < (sample size) < 100000, if (minimum sequence length) < 2, if (maximum
      sequence length) < 1000 and if (cumulative sequence length) < 1000000, an
      object of type semi-markov_data is returned, otherwise no object is returned.
    * The returned object contains both the simulated sequences and the model 
      used for simulation.
    * If 0 < size < 100000 and if 0 < nb_trial < 1000, an object of type tops
      is returned, otherwise no object is returned. The returned object contains
      both the simulated 'tops' and the model used for simulation.

    
    :Description:
    
    If the fourth argument is of type time_events or renewal_data, the simulated sample has the same distribution of length of the observation period than the original sample given by this fourth argument. This simulation mode is particularly useful to study the effects of length biasing.

   If the third argument is of type discrete_sequences, markov_data or semi-markov_data, the simulated sequences has the same length distribution than the original sample given by this third argument. This simulation mode is particularly useful to study the effects of length biasing.

    .. seealso::
    
        :func:`~openalea.stat_tool.distribution.Distribution`,
        :func:`openalea.stat_tool.mixture.Mixture`,
        :func:`openalea.stat_tool.convolution.Convolution`,
        :func:`openalea.stat_tool.compound.Compound`,
        :func:`openalea.stat_tool.data_transform.ExtractHistogram`.
        :func:`~openalea.sequence_analysis.Markov`,
        :func:`~openalea.sequence_analysis.semi_markov.SemiMarkov`,
        :func:`~openalea.sequence_analysis.hidden_markov.HiddenMarkov`,
        :func:`~openalea.sequence_analysis.hidden_semi_markov.HiddenSemiMarkov`,
        :func:`~openalea.sequence_analysis.renewal.Renewal`,
        :func:`~openalea.sequence_analysis.top_parameters.TopParameters`,    
    """
    
    # switch input obj argument: 
    
    # standard distribution case
    if len(args)==1  and isinstance(args[0], int):
        return SimulateDistribution(obj, args[0])
    # top parameters case
    elif isinstance(obj, _sequence_analysis._Top_parameters):
        try:
            arg1 = args[0]
            arg2 = args[1]
        except TypeError:
            raise TypeError("request two arguments with top_parameters simulation")
    
        NbAxillary = kargs.get("NbAxillary", 1)
        
        if isinstance(arg1,int) and isinstance(arg2,int):    
            return obj.simulation(arg1, arg2, NbAxillary)
        else:
            raise TypeError("With top_parameters simulation, second and third arguments must be integers")
    # Renewal case
    elif isinstance(obj, _sequence_analysis._Renewal):
        itype = args[0]
        if itype == 'Ordinary':
            Type = 'o'
        elif itype == 'Equilibrium':
            Type = 'e'
        else:
            raise KeyError("first argument must be Equilibrium or Ordinary")  
        
        if isinstance(args[1], int) and isinstance(args[2], int):
              return obj.simulation_nb_elements(Type, args[1], args[2])
        elif isinstance(args[1], int):
              return obj.simulation_time_events(Type, args[1], args[2])
        else:
            return obj.simulation_histogram(Type, args[1])
    # other cases (Markovian, semi_markov, hidden_semi_markov and so on
    else:
        CountingFlag = kargs.get("CountingFlag", True)
        
        #order of the if statements is important ! Keep it that way
        if isinstance(args[0], int) and isinstance(args[1], int):
            if isinstance(obj, _sequence_analysis._Semi_markov) or isinstance(obj, _sequence_analysis._Hidden_semi_markov):
                return obj.simulation_nb_elements(args[0], args[1], CountingFlag)
            else:
                return obj.simulation_nb_elements(args[0], args[1], True)
        #here the second arguments is data structure such as Sequences
        elif isinstance(args[0], int):
            if isinstance(obj, _sequence_analysis._Semi_markov) or isinstance(obj, _sequence_analysis._Hidden_semi_markov):
                return obj.simulation_markovian_sequences(args[0], args[1], CountingFlag)
            else:
                return obj.simulation_markovian_sequences(args[0], args[1], True)
        # first argument is a compound_data or equivalent
        else:
            if isinstance(obj, _sequence_analysis._Semi_markov) or isinstance(obj, _sequence_analysis._Hidden_semi_markov):
                return obj.simulation_histogram(args[0], CountingFlag, False)
            else:
                return obj.simulation_histogram(args[0], True, False)
            
                
        
         
    
    
        
   
    
    
    
