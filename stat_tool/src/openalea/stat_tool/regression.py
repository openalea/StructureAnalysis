"""Regression"""
__revision__ = "$Id$"

import interface
import _stat_tool

from _stat_tool import _Regression

__all__ = ['_Regression',
           'Regression']

def Regression(vec, type, 
               explanatory_variable, response_variable, 
               param=None, 
               Algorithm="Averaging", Weighting=False):
    """
    Simple regression (with a single explanatory variable).

    :Parameters:
    
      * vec : vectors
        vectors
      * type : string
        `"Linear"` or `"MovingAverage"` or `"NearestNeighbors"`
      * explanatory_variable : int
        index of the explanatory variable
      * response_variable : int 
        index of the response variable
      * filter : list of float
        filter values on the half width i.e. from one extremity to the central value
        (with the constraint filter[i] + filter[m] = 1),
      * frequencies : list of float
        frequencies defining the filter,
      * dist : distribution, mixture, convolution, compound 
        symmetric distribution, whose size of the support is even, defining the filter 
        (for instance Distribution("BINOMIAL",0,4,0.5)),
      * span : float
        proportion of individuals in each neighbourhood. 

    :Keywords:
    
      * Algorithm : string
          - `"Averaging"` (default) 
          - `"LeastSquares"`
            
        This optional argument can only be used if the second mandatory argument
        specifying the regression type is "MovingAverage".

      * Weighting : bool 
        weighting or not of the neighbors according to their distance to the 
        computed point (default value: True). This optional argument can only
        be used if the second mandatory argument specifying the regression type 
        is "NearestNeighbors". 

    :Returns:
    
        An object of type regression is returned.

    :Examples:

    .. doctest::
        :options: +SKIP
    
        >>> Regression(vec, "Linear", explanatory_variable, response_variable)
        >>> Regression(vec, "MovingAverage", explanatory_variable, response_variable, filter, Algorithm="LeastSquares")
        >>> Regression(vec, "MovingAverage", explanatory_variable, response_variable, frequencies, Algorithm="LeastSquares")
        >>> Regression(vec, "MovingAverage", explanatory_variable, response_variable, dist, Algorithm="LeastSquares")
        >>> Regression(vec, "NearestNeighbors", explanatory_variable, response_variable, span, Weighting=False)
    
    .. seealso::
        :func:`~openalea.stat_tool.output.Plot`

    """

    algo_map = { 'Averaging' : 'a',
                 'LeastSquares' : 's',
                 }

    try:
        Algorithm = algo_map[Algorithm]
    except KeyError:
        raise KeyError("Bad Algorithm. Possible algorithms are %s"%(algo_map.keys(),))
        

    if(type == "Linear"):
        return vec.linear_regression(explanatory_variable, response_variable, )

    elif(type == "MovingAverage"):
        # param is a list of float or a Distribution
        if isinstance(param[0], int):
            return vec.moving_average_regression_values(explanatory_variable, 
                            response_variable, param, Algorithm)
        else:
            return vec.moving_average_regression_distribution(explanatory_variable, 
                            response_variable, param, Algorithm)
            
            

    elif(type == "NearestNeighbors"):
        span = param
        return vec.nearest_neighbours_regression(explanatory_variable, response_variable, 
                                                 float(span), bool(Weighting))

    else:
        raise TypeError("Bad Regression type. Must be 'Linear' or " +
                        "'MovingAverage' or 'NearestNeighbors'" )



# Extend _Regression class dynamically
interface.extend_class( _stat_tool._Regression, interface.StatInterface)




