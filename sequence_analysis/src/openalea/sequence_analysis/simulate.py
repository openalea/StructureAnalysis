        

#import openalea.sequence_analysis._sequence_analysis as _sequence_analysis
import openalea.stat_tool._stat_tool as _stat_tool
import _sequence_analysis

from openalea.stat_tool.simulate import Simulate as SimulateDistribution
 
def Simulate(obj, *args, **kargs):
    """
    
    """
    
    # switch input obj argument: 
    
    # standard distribution case
    if len(args)==1 and len(kargs==0) and isinstance(args[0], int):
        return SimulateDistribution(obj, args[0])
    # top parameters case
    elif isinstance(obj, _sequence_analysis._Top_parameters):
        try:
            arg1 = args[0]
            arg2 = args[1]
        except TypeError:
            raise TypeError("request two arguments With top_parameters simulation")
    
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
                return obj.simulation_nb_sequences(args[0], args[1], CountingFlag)
            else:
                return obj.simulation_nb_sequences(args[0], args[1], True)
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
            
                
        
         
    
    
        
   
    
    
    