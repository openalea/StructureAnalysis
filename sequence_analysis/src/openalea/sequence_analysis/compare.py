"""
compare tools for sequences
"""
import openalea.sequence_analysis._sequence_analysis as _sequence_analysis
import openalea.stat_tool._stat_tool as _stat_tool
#from distribution import Distribution
from openalea.stat_tool.comparison import compare_histo, compare_vectors

def compare_markovian_sequences_and_variable_order_markov(obj, *args, **kargs):

    filename = kargs.get("Filename", None)
    markov_list = []
    
    for arg in args:
        if (isinstance(arg, _sequence_analysis._Variable_order_markov)):
            markov_list.append(arg)

    if(isinstance(obj, _sequence_analysis._Markovian_sequences)):
        return obj.comparison_variable_order_markov(markov_list, filename)


def compare_variable_order_markov_and_sequences(obj, *args, **kargs):
    """not implemented"""

    variable_order_markov = []
    markov_sequence = []
    nb_seq = None
    filename = kargs.get("Filename", None)
    
    if(isinstance(obj, _sequence_analysis._Variable_order_markov)):
        for arg in args:
            if (isinstance(arg, _sequence_analysis._Variable_order_markov)):
                variable_order_markov.append(arg)
            elif (isinstance(arg, _sequence_analysis._Markovian_sequences)):
                markov_sequence.append(arg)
            elif (isinstance(arg, int)):
                  nb_seq = arg
            else:
                raise KeyError("arguments must be alternance of Variable_order_markov and Markovian_sequences %s provided", type(arg))
                
        
        if nb_seq is None:
            raise KeyError("number of sequences must be provided")
        return obj.divergence_computation(variable_order_markov, 
                                          markov_sequence, nb_seq, filename)


def compare_markov(mc, *args, **kargs):
    """not implemented"""
    raise NotImplementedError()



def Compare(arg1, *args, **kargs):
    """Comparison functions factory

    :Parameters:
      - `arg1` should be in : 
          * `compare_histo` : Histograms comparison
          * `compare_vectors` : Vectors comparison
          * `compare_seq` : Sequences comparison
          * `compare_markov` : Markovian models comparison

    .. seealso::
        :func:`~openalea.stat_tool.comparison.compare_histo`,
        :func:`~openalea.stat_tool.comparison.compare_vectors`
        :func:`~openalea.stat_tool.comparison.compare_seq`
        :func:`~openalea.stat_tool.comparison.compare_markov`

    .. todo:: Get the AMAPMod documentation
    
    
    """

    p1 = arg1
    
    if(isinstance(p1, _stat_tool._Histogram)):
        return compare_histo(arg1, *args, **kargs)
    
    elif(isinstance(p1, _stat_tool._Vectors)):
        return compare_vectors(arg1, *args, **kargs)
    
    elif(isinstance(p1, _sequence_analysis._Variable_order_markov)):
        return compare_variable_order_markov_and_sequences(arg1, *args, **kargs)

    elif(isinstance(p1, _sequence_analysis._Markovian_sequences) and
         isinstance(args[0], _sequence_analysis._Variable_order_markov)):
        return compare_markovian_sequences_and_variable_order_markov(arg1, *args, **kargs)
    
