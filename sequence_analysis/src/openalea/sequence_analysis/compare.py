"""
compare tools for sequences
"""


import openalea.stat_tool._stat_tool as _stat_tool
from openalea.stat_tool.comparison import compare_histo, compare_vectors
import _sequence_analysis

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


def compare_sequences(seq, *args, **kargs):
    #int indel_cost = ADAPTATIVE , algorithm = AGGLOMERATIVE;
    #double indel_factor , transposition_factor = TRANSPOSITION_FACTOR;
    INDEL_FACTOR_1 = 0.51
    INDEL_FACTOR_N = 0.51
    TRANSPOSITION_FACTOR = 0.
    INDEL_DISTANCE = 1.0
    
    RefSequence = kargs.get("RefSequence", -1)
    TestSequence = kargs.get("TestSequence", -1)
    Begin = kargs.get("Begin", "Aligned")
    End = kargs.get("End", "Aligned")
    FileName = kargs.get("FileName", None)
    Format = kargs.get("Format", 'a') # a for ASCII
    AlignmentFormat = kargs.get("AlignmentFormat", "a") # a for ASCII
    AlignmentFileName = kargs.get("AlignmentFileName", None)
    IndelCost = kargs.get("IndelCost", "Adaptative")
    IndelFactor = kargs.get("IndelFactor", INDEL_FACTOR_1)
    TranspositionFlag = kargs.get("Transposition", False)
    TranspositionFactor = kargs.get("TranspositionFactor", TRANSPOSITION_FACTOR)
    
    begin_end_map = {"Aligned":False, "Free":True} 
    indelcost_map = {"Adaptative":_sequence_analysis.IndelCost.ADAPTATIVE, 
                     "Fixed": _sequence_analysis.IndelCost.FIXED}    
    try:
        Begin = begin_end_map[Begin]
        End = begin_end_map[End]
    except KeyError:
        raise KeyError("wrong Begin or End argument. Use one of %" % begin_end_map.keys())
    
    try:
        IndelCost = indelcost_map[IndelCost]
    except KeyError:
        raise KeyError("wrong Begin or End argument. Use one of %" % indelcost_map.keys())
        
    print IndelCost
    
    if len(args)==1:
        print args[0]
        if isinstance(args[0], _stat_tool._VectorDistance):
            dist_matrix = seq.alignment_vector_distance(args[0], RefSequence, TestSequence, Begin, End,
                                    IndelCost, IndelFactor, TranspositionFlag,
                                    TranspositionFactor, FileName, Format,
                                    AlignmentFileName, AlignmentFormat)
    else:
        dist_matrix = seq.alignment(RefSequence, TestSequence, Begin ,
                                  End , FileName , Format , AlignmentFileName,
                                  AlignmentFormat)
    return dist_matrix
        

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
    
    if (isinstance(p1, _stat_tool._Histogram)):
        return compare_histo(arg1, *args, **kargs)
    
    elif (isinstance(p1, _stat_tool._Vectors)):
        return compare_vectors(arg1, *args, **kargs)
    
    elif (isinstance(p1, _sequence_analysis._Variable_order_markov)):
        return compare_variable_order_markov_and_sequences(arg1, *args, **kargs)

    elif (isinstance(p1, _sequence_analysis._Markovian_sequences)):
         
          if len(args)== 0:
              return compare_sequences(arg1, *args, **kargs)
          if len(args) > 0:
              if (isinstance(args[0], _sequence_analysis._Variable_order_markov)):
                  return compare_markovian_sequences_and_variable_order_markov(arg1, *args, **kargs)
              elif  (isinstance(args[0], _stat_tool._VectorDistance)):
                  return compare_sequences(arg1, *args, **kargs)
    
    
    
