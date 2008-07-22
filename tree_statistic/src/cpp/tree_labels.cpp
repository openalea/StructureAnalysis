/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       AMAPmod: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2002 UMR Cirad/Inra Modelisation des Plantes
 *
 *       File author(s): Y. Guedon (yann.guedon@cirad.fr)
 *
 *       $Source: /usr/cvsmaster/AMAPmod/src/STAT_TREES/src/tree_labels.cpp,v $
 *       $Id: tree_labels.cpp 3186 2007-05-25 15:10:30Z dufourko $
 *
 *       Forum for AMAPmod developers: amldevlp@cirad.fr
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



/****************************************************************
 *
 *  Mots cles (format des fichiers) :
 */

namespace Stat_trees
{

const char *STAT_TREES_word[] = {
  "HIDDEN_MARKOV_TREE" ,               // STATW_HIDDEN_MARKOV_TREE
  "EQUILIBRIUM_HIDDEN_MARKOV_TREE" ,   // STATW_EQUILIBRIUM_HIDDEN_MARKOV_TREE
  "HIDDEN_MARKOV_OUT_TREE" ,           // STATW_HIDDEN_MARKOV_OUT_TREE
  "EQUILIBRIUM_HIDDEN_MARKOV_OUT_TREE" // STATW_EQUILIBRIUM_HIDDEN_MARKOV_OUT_TREE
};


const char *STAT_TREES_type[] = {
  "INTEGER VALUE" ,
  "REAL VALUE" ,
  "STATE" ,
  "TIME" ,
  "TIME INTERVAL" ,
  "POSITION" ,
  "POSITION INTERVAL" ,
  "NB INTERNODE"
};


/****************************************************************
 *
 *  Labels :
 */


const char *STAT_TREES_label[] = {
  "log-likelihood of the state tree" ,  // STATL_STATE_TREE_LIKELIHOOD
  "log-likelihood of the state trees" , //  STATL_STATE_TREES_LIKELIHOOD
  "log-likelihood of the observed trees" , //  STATL_OBSERVED_TREES_LIKELIHOOD
  "information of the trees in the iid case" , // STATL_TREES_IID_INFORMATION

  "smoothed" , // STATL_SMOOTHED
  "observed" , //  STATL_OBSERVED


  "probability of no-occurrence of state" , // STATL_STATE_NO_OCCURRENCE
  "probability of no-occurrence of output" , //   STATL_OUTPUT_NO_OCCURRENCE
  "absorption probability of state" , // STATL_STATE_ABSORPTION
  "absorption probability of output" , // STATL_OUTPUT_ABSORPTION

  "path length (starting from root) up to the first occurrence of state" , //STATL_STATE_FIRST_OCCURRENCE_ROOT
  "path length (starting from root) up to the first occurrence of output" , // STATL_OUTPUT_FIRST_OCCURRENCE_ROOT
  "path length (starting from terminal node) up to the first occurrence of state" , // STATL_STATE_FIRST_OCCURRENCE_LEAVES
  "path length (starting from terminal nodes) up to the first occurrence of output" , // STATL_OUTPUT_FIRST_OCCURRENCE_LEAVES
  "sojourn size" , // STATL_SOJOURN_SIZE
  "mixture of " , // STATL_MIXTURE_OF
  "number of zones of state" ,
  "number of zones of output" ,
  "number of occurrences of state" , // STATL_STATE_NB_OCCURRENCE
  "number of occurrences of output" , // STATL_OUTPUT_NB_OCCURRENCE
  "per tree" ,
  "per size" ,
  "state probabilities" , // STATL_STATE_PROBABILITY
  "posterior state probabilities" ,
  "posterior in state probabilities" ,
  "posterior out state probabilities" ,
  "conditional entropy" ,
  "marginal entropy" ,
  "sum of marginal entropies" ,
  "partial state tree entropy" ,
  "state tree entropy" ,
  "Gini index" ,
  "upper bound" ,
  "number of state trees" ,
  "maximum posterior state probabilities" ,
  "maximum posterior in state probabilities" ,
  "maximum posterior out state probabilities" ,
  "ambiguity" , // STATL_AMBIGUITY

  "tree" ,
  "trees" ,
  "vertex" ,

  "tree size" ,
  "cumulative size" ,
  "number of children" ,
  "cumulative number of children",
  "time" // STATL_TIME
};



/****************************************************************
 *
 *  Messages d'erreur pour l'analyse des fichiers :
 */





/****************************************************************
 *
 *  Messages d'erreur de traitement :
 */


const char *STAT_TREES_error[] = {
  "bad number of states" , // SEQR_NB_STATE
  "bad number of values per variable" , // STATR_VARIABLE_NB_VALUE
  "bad number of output processes" ,  // STATR_NB_OUTPUT_PROCESS
  "bad number of integer output processes" ,  // STATR_NB_INT_OUTPUT_PROCESS
  "bad number of floating output processes" ,  // STATR_NB_REAL_OUTPUT_PROCESS
  "bad output process index" ,  // STATR_OUTPUT_PROCESS_INDEX
  "bad number of outputs" , // STATR_NB_OUTPUT
  "bad self-transition probability" , // STATR_SELF_TRANSITION
  "bad number of state trees", // STATR_NB_STATE_TREES
  "bad number of trees" , // STATR_NB_TREES
  "bad tree identifier" ,
  "number of children too high" ,
  "number of children too low" ,
  "bad number of children" , // STATR_TREE_NB_CHILDREN
  "bad tree size" , // STATR_TREE_SIZE
  "bad vertex identifier" ,
  "child IDs do not match" ,
  "tree size too small" ,
  "tree size too big" ,
  "tree cumulative size too big" ,
  "tree cumulative number of children too high" ,
  "failure in the computation of the optimal state trees" ,

  "not present" , // STATR_NOT_PRESENT

  "vector does not define a permutation" , // STATR_NO_PERMUTATION

  "bad exponent for saem" , // STATR_SAEM_EXP
  "bad estimation algorithm" , // STATR_EM_ALGORITHM

  "bad variable type" , // STAR_VARIABLE_TYPE
  "variable 1: bad variable type" , // STAR_VARIABLE_1_TYPE
  "bad minimum tree size" , // STATR_MIN_TREE_SIZE
  "bad maximum tree size" , // STATR_MAX_TREE_SIZE
  "bad number of selected values" , // STATR_NB_SELECTED_VALUE

  "state trees not in the data",
  "characteristics not computed" , // STATR_CHARACTERISTICS_NOT_COMPUTED
  "non-existing characteristic distribution" // STATR_NON_EXISTING_CHARACTERISTIC_DISTRIBUTION
};

}; // end namespace
