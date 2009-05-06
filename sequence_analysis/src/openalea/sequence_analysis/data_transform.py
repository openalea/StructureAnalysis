"""data transform"""
__revision__ = "$Id: vectors.py 6217 2009-04-08 12:40:15Z cokelaer $"


#def MergeVariable(obj, *args, **kargs):
# RefSample = kargs.get("RefSample", -1)

#    return obj.merge_variable(list(args), RefSample) 
def RemoveRun(obj, variable=1, value=None, position=None, MaxLength=4):
    """RemoveRun

	Removal of the first or last run of a given value (for a given variable) in a sequence.

    :Usage:

	RemoveRun(seq1, value, position, MaxLength->4)
	RemoveRun(seqn, variable, value, position, MaxLength->4)	

    :Arguments:
	    seq1 (sequences, discrete_sequences, markov_data, semi-markov_data): univariate sequences,
	    seqn (sequences, discrete_sequences, markov_data, semi-markov_data): multivariate sequences,
	    variable (int): variable index,
	    value (int): value,
	    position (string): position of the removed run: "Begin" or "End".
	
    :Optional Arguments:
    	MaxLength (int): maximum length of the removed run (default behaviour: the entire run is removed).
	
    :Returned Object:

	If variable is a valid index of variable, if value is a possible value and if MaxLength > 0, an object of type sequences or discrete_sequences is returned, otherwise no object is returned. The returned object is of type discrete_sequences if all the variables are of type STATE, if the possible values for each variable are consecutive from 0 and if the number of possible values for each variable is <= 15.

    :Examples:
    >>> RemoveRun(seq1, 0, "End") 	
    >>> RemoveRun(seq5, 2, 0, "End")

    :See Also:
	
	AddAbsorbingRun,
	Cluster, 
	Cumulate, 
	Difference, 
	IndexExtract, 
	LengthSelect, 
	Merge, 
	MergeVariable, 
	MovingAverage, 
	RecurrenceTimeSequences, 
	Reverse, 
	SegmentationExtract, 
	SelectIndividual, 
	SelectVariable, 
	Shift, 
	Transcode, 
	ValueSelect,
	VariableScaling.
    """

    # arguments : indice de la variable,
    #             valeur, position ('b' : begin, 'e' : end),
    #             longueur maximum de la serie supprimee.
    
    if position == 'End': position = 'e'
    if position == 'Begin': position = 'b'
    
    if position not in ['b','e']:
        raise "position must be 'e' for end or 'b' for begin"
    return obj.remove_run(variable, value, position, MaxLength)

