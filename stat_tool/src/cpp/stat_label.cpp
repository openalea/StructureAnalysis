/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       V-Plants: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2017 CIRAD/INRA/Inria Virtual Plants
 *
 *       File author(s): Yann Guedon (yann.guedon@cirad.fr)
 *
 *       $Source$
 *       $Id$
 *
 *       Forum for V-Plants developers:
 *
 *  ----------------------------------------------------------------------------
 *
 *                      GNU General Public Licence
 *
 *       This program is free software; you can redistribute it and/or
 *       modify it under the terms of the GNU General Public License as
 *       published by the Free Software Foundation; either version 2 of
 *       the License, or (at your option) any later version.
 *
 *       This program is distributed in the hope that it will be useful,
 *       but WITHOUT ANY WARRANTY; without even the implied warranty of
 *       MERCHANTABILITY or FITNESS For A PARTICULAR PURPOSE. See the
 *       GNU General Public License for more details.
 *
 *       You should have received a copy of the GNU General Public
 *       License along with this program; see the file COPYING. If not,
 *       write to the Free Software Foundation, Inc., 59
 *       Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 *
 *  ----------------------------------------------------------------------------
 */



namespace stat_tool {


/****************************************************************
 *
 *  keywords (file format)
 */


const char *STAT_word[] = {
  "INF_BOUND" ,
  "SUP_BOUND" ,
  "PARAMETER" ,
  "PROBABILITY" ,
  "NO_SEGMENT" ,
  "SEQUENCE_LENGTH" ,

  "SHAPE" ,
  "SCALE" ,
  "ZERO_PROBABILITY" ,
  "MEAN" ,
  "STANDARD_DEVIATION" ,
  "MEAN_DIRECTION" ,
  "CONCENTRATION" ,
  "INTERCEPT" ,
  "SLOPE" ,
  "AUTOREGRESSIVE_COEFF" ,

  "MIXTURE" ,
  "CONVOLUTION" ,
  "DISTRIBUTIONS" ,
  "DISTRIBUTION" ,
  "WEIGHT" ,
  "WEIGHTS" ,
  "COMPONENTS" ,
  "COMPONENT" ,

  "COMPOUND_DISTRIBUTION" ,
  "SUM_DISTRIBUTION" ,
  "ELEMENTARY_DISTRIBUTION" ,

  "VARIABLE" ,
  "VARIABLES" ,

  "DISTANCE" ,
  "CATEGORIES" ,
  "PERIOD" ,

  "STATES" ,
  "INITIAL_PROBABILITIES" ,
  "TRANSITION_PROBABILITIES" ,

  "STATE" ,
  "FUNCTION" ,

  "OUTPUT_PROCESS" ,
  "OUTPUT_PROCESSES" ,

  "NONPARAMETRIC" ,   // for ascending compatibility

  "CATEGORICAL" ,

  "PARAMETRIC" ,   // for ascending compatibility

  "DISCRETE_PARAMETRIC" ,
  "CONTINUOUS_PARAMETRIC" ,
  "OBSERVATION_DISTRIBUTION" ,
  "OBSERVATION_MODEL" ,

  "OBSERVATION_PROBABILITIES" ,   // for ascending compatibility

  "OUTPUT"
};


const char *STAT_discrete_distribution_word[] = {
  "CATEGORICAL" ,
  "BINOMIAL" ,
  "POISSON" ,
  "NEGATIVE_BINOMIAL" ,
  "POISSON_GEOMETRIC" ,
  "UNIFORM" ,
  "PRIOR_SEGMENT_LENGTH"
};


const char *STAT_discrete_distribution_letter[] = {
  "C" ,
  "B" ,
  "P" ,
  "NB" ,
  "PG" ,
  "U" ,
  "PSL"
};


const char *STAT_continuous_distribution_word[] = {
  "GAMMA" ,
  "INVERSE_GAUSSIAN" ,
  "GAUSSIAN" ,
  "VON_MISES" ,
  "ZERO_INFLATED_GAMMA" ,
  "LINEAR_MODEL" ,
  "AUTOREGRESSIVE_MODEL"
};


const char *STAT_continuous_distribution_letter[] = {
  "Ga" ,
  "IG" ,
  "G" ,
  "VM" ,
  "ZIGa" ,
  "LM" ,
  "AR"
};


const char *STAT_variable_type_word[] = {
  "NOMINAL" ,
  "ORDINAL" ,
  "NUMERIC" ,
  "CIRCULAR"
};


const char *STAT_variable_type_letter[] = {
  "No" ,
  "O" ,
  "Nu" ,
  "C"
};


const char *STAT_distance_word[] = {
  "ABSOLUTE_VALUE" ,
  "QUADRATIC"
};


const char *STAT_criterion_word[] = {
  "AIC" ,
  "AICc" ,
  "BIC" ,
  "BICc" ,
  "ICL" ,
  "ICLc" ,
  "mBIC" ,
  "log-likelihood slope" ,
  "dimension jump" ,
  "segmentation log-likelihood slope" ,
  "segmentation dimension jump"
};


const char *STAT_function_word[] = {
  "LINEAR" ,
  "LOGISTIC" ,
  "MONOMOLECULAR"
};


const char *STAT_variable_word[] = {
  "INT" ,
  "REAL" ,
  "STATE" ,
  "VALUE" ,         // for ascending compatibility
//  "NB_INTERNODE" ,
  "AUXILIARY"
};



/****************************************************************
 *
 *  Labels
 */


const char *STAT_label[] = {
  "line" ,
  "word" ,

  "<Return>: continue." ,
  "End." ,

  "model" ,

  "mean" ,
  "median" ,
  "mode" ,
  "variance" ,
  "standard deviation" ,
  "quantile" ,
  "lower quartile" ,
  "upper quartile" ,
  "mean absolute deviation" ,
  "coefficient of variation" ,
  "variance/mean ratio" ,
  "coefficient of skewness" ,
  "coefficient of kurtosis" ,
  "information" ,
  "coefficient of concentration" ,
  "smoothness" ,

  "mean direction" ,
  "mean resultant length" ,
  "circular standard deviation" ,

  "mean confidence interval" ,

  "one-sided" ,
  "two-sided" ,
  "chi-square test" ,
  "F-test" ,
  "t-test" ,
  "Wilcoxon-Mann-Whitney test" ,
  "Kruskal-Wallis test" ,
  "degree of freedom" ,
  "degrees of freedom" ,
  "standard normal value" ,
  "chi-square value" ,
  "F-value" ,
  "t-value" ,
  "critical probability" ,
  "reference" ,
  "test" ,
  "P(X1 < X2)" ,
  "P(X1 = X2)" ,
  "P(X1 > X2)" ,

  "fit" ,
  "log-likelihood" ,
  "normalized" ,
  "maximum possible log-likelihood" ,
  "deviance" ,
  "log-likelihood for the optimal classification" ,
  "maximum possible log-likelihood for the optimal classification" ,
  "classification entropy" ,
  "free parameter" ,
  "free parameters" ,
  "penalyzed log-likelihood" ,
  "iteration" ,
  "iterations" ,

  "distribution" ,
  "distributions" ,
  "unproper" ,
  "complementary probability" ,
  "cumulative" ,
  "survivor" ,
  "function" ,
  "matching" ,
  "concentration" ,
  "curve" ,

  "mixture" ,
  "weight" ,
  "component" ,
  "categorical" ,
  "discrete parametric" ,
  "distances between components" ,
  "posterior assignment probability" ,
  "assignment entropy" ,  
  "convolution" ,
  "compound" ,
  "sum" ,
  "elementary" ,
  "death probability" ,
  "survival probability" ,
  "frequency" ,

  "backward" ,
  "forward" ,
  "inter-event time within the observation period" ,
  "sojourn time" ,

  "frequency distribution" ,
  "frequency distributions" ,
  "histogram" ,
  "bin width" ,
  "sample" ,
  "sample size" ,
  "dissimilarities" ,

  "information ratio" ,
  "clustering step" ,
  "category" ,

  "vector" ,
  "vectors" ,
  "number of vectors" ,
  "variable" ,
  "identifier" ,
  "minimum value" ,
  "maximum value" ,
  "marginal" ,
  "variance-covariance matrix" ,
  "correlation matrix" ,
  "Spearman rank correlation matrix" ,
  "Kendall rank correlation matrix" ,
  "limit correlation coefficient" ,
  "Spearman limit rank correlation coefficient" ,
  "Kendall limit rank correlation coefficient" ,
  "sup norm distance" ,
  "overlap" ,

  "contingency table" ,
  "deviation table" ,
  "chi-square contribution table" ,

  "analysis of variance" ,
  "source of variation" ,
  "between samples" ,
  "within samples" ,
  "total" ,
  "sum of squares" ,
  "mean square" ,

  "intercept" ,
  "slope" ,
  "explanatory variable" ,
  "response variable" ,
  "correlation coefficient" ,
  "R-squared" ,
//  "regression variation / total variation" ,
  "coefficient of determination" ,
  "95% confidence limits for a null autoregressive coefficient" ,
  "residual" ,
  "standardized residual" ,
  "linear" ,
  "linear model" ,
  "nonparametric" ,
  "regression" ,
  "estimation" ,

  "distance" ,
  "cumulated distance per" ,
  "length" ,
  "deletion" ,
  "insertion" ,
  "insertion/deletion" ,
  "match" ,
  "substitution" ,
  "transposition" ,
  "rate" ,

  "matrix" ,
  "number of rows" ,
  "number of columns" ,

  "verified symmetry rate" ,
  "verified triangle inequality rate" ,

  "cluster" ,
  "clusters" ,
  "neighbor" ,
  "most distant" ,
  "diameter" ,
  "separation" ,
  "pattern level" ,
  "isolated" ,
  "non-isolated" ,
  "within" ,
  "between" ,
  "ratio" ,
  "prototype" ,
  "step" ,
  "child" ,
  "composition" ,
  "child cluster distance scale coefficient" ,
  "diameter scale coefficient" ,
  "dendrogram scale" ,

  "class" ,
  "state" ,
  "states" ,
  "transient" ,
  "recurrent" ,
  "absorbing" ,
  "memory" ,
  "stationary probabilities" ,

  "value" ,
  "values" ,
  "output" ,
  "output process" ,
  "observation" ,
  "observation probability matrix" ,
  "theoretical" ,
  "restoration" ,
  "weights" ,
  "distances between observation distributions" ,
  "distances between observation distributions for consecutive states" ,
  "q-q plot"
};



/****************************************************************
 *
 *  Error messages for lexical analysis of files
 */


const char *STAT_parsing[] = {
  "bad keyword" ,
  "format error" ,

  "bad distribution name" ,
  "bad parameter name" ,
  "bad parameter index" ,
  "bad number of parameters" ,
  "bad separator" ,
  "bad parameter value" ,
  "sum of probabilities different to 1" ,

  "bad data type" ,
  "bad number of tokens per line" ,
  "empty sample" ,

  "values not ordered" ,
  "value too large" ,

  "bad number of distributions" ,
  "bad distribution index" ,
  "bad weight value" ,

  "bad number of components" ,
  "bad component index" ,

  "bad number of variables" ,
  "bad variable index" ,
  "bad variable type" ,

  "bad number of categories" ,
  "bad local distance" ,
  "triangle inequality not verified" ,
  "bad period value" ,

  "bad number of states" ,
  "bad order" ,
  "bad initial probability" ,
  "bad transition probability" ,
  "Markov chain structure not allowed" ,
  "irreducible" ,

  "bad state index" ,
  "bad output index" ,
  "bad observation probability" ,
  "bad number of observation distributions" ,
  "non-consecutive outputs" ,
  "overlap of observation distributions" ,
  "non-overlap of observation distributions" ,
  "bad observation distribution type" ,
  "bad number of output processes" ,
  "bad output process index" ,
  "vector does not define a permutation"  // STATR_NO_PERMUTATION
};



/****************************************************************
 *
 *  Error messages
 */


const char *STAT_error[] = {
  "bad file name" ,
  "bad file prefix" ,
  "empty sample" ,
  "bad size of the sample" ,
  "bad sizes of the samples" ,
  "non-existing distribution" ,
  "non-existing frequency distribution" ,
  "no data associated with the model" ,

  "too many distributions" ,
  "too many histograms" ,
  "null variance: plot not possible" ,

  "distribution range of values incompatible with frequency distribution range of values" ,
  "bad minimum inferior bound" ,
  "estimation failure" ,

  "number of complete intervals too small" ,
  "complete time interval: bad minimum value" ,
  "forward recurrence time interval: bad minimum value" ,
  "no event observation period: bad minimum value" ,
  "bad duration mean estimation method" ,
  "initial inter-event distribution support incompatible with data" ,

  "bad distribution index" ,
  "bad frequency distribution index" ,
  "bad number of distributions" ,
  "bad weight initialization step" ,
  "inferior bound: distribution flag incompatible with mixture flag" ,

  "not present" ,
  "bad output process index" ,
  "bad output process type" ,
  "optimal assignment not in the data" ,

  "bad number of values per variable" ,
  "bad number of output processes" ,
  "bad number of outputs" ,
  "bad minimum value: should be positive" ,

  "bad minimum value of known distribution" ,
  "bad mean of known distribution" ,
  "bad number of iterations" ,
  "bad number of assignments" ,
  "bad penalty weight" ,

  "bad clustering step" ,
  "bad number of classes" ,
  "bad shift value" ,
  "bad threshold value" ,
  "bad rounded value" ,
  "smaller than" ,
  "greater than" ,
  "not allowed" ,
  "bad number of categories" ,
  "non-consecutive categories" ,
  "missing category" ,
  "bad cluster limit" ,
  "bad information ratio: should be between 0 and 1" ,
  "null reference information: variance should be strictly positive" ,
  "bad scaling coefficient" ,
  "bad value" ,
  "bad minimum value" ,
  "bad maximum value" ,

  "marginal frequency distribution not built" ,
  "marginal histogram not built" ,
  "bad histogram bin width" ,
  "bad histogram min value" ,
  "bad" ,
  "bad number of vectors" ,
  "bad number of variables" ,
  "bad number of selected variables" ,
  "bad number of summed variables" ,
  "bad variable type" ,
  "bad variable index" ,
  "bad variable indices: should be different" ,
  "already used" ,
  "already selected" ,
  "bad sample index" ,
  "rank correlation coefficient computation not possible" ,
  "should be shifted or scaled down" ,
  "number of possible values incompatible with period" ,
  "bad number of possible values" ,
  "missing value" ,
  "least-square algorithm not allowed" ,
  "regression failure" ,

  "bad matrix dimensions" ,
  "only infinite distances" ,
  "matrix already symmetrical" ,
  "unsymmetrical matrix" ,
  "dissimilarity measures not normalized by lengths" ,

  "symmetry not verified" ,
  "triangle inequality not verified" ,

  "bad matrix structure" ,
  "a square matrix with dimension greater than 1" ,
  "only clusters with a single element" ,
  "bad number of clusters" ,
  "bad pattern type" ,
  "bad prototype identifier" ,
  "bad number of" ,

  "initial probability equal to 0: should be stricly positive" ,

  "odd" ,
  "non-symmetrical distribution" ,
  "unproper distribution"
};


};  // namespace stat_tool
