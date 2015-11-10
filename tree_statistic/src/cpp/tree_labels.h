/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       V-Plants: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2010 CIRAD/INRIA Virtual Plants
 *
 *       File author(s): J.-B. Durand (jean-baptiste.durand@imag.fr)
 *
 *       $Source$
 *       $Id: tree_labels.h 3186 2007-05-25 15:10:30Z dufourko $
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



#ifndef TREE_LABELS_H
#define TREE_LABELS_H

namespace Stat_trees
{

/****************************************************************
 *
 *  Keyword identifiers (file format):
 */


enum {
  TREESTATW_HIDDEN_MARKOV_TREE ,
  TREESTATW_EQUILIBRIUM_HIDDEN_MARKOV_TREE ,
  TREESTATW_HIDDEN_MARKOV_IND_OUT_TREE ,
  TREESTATW_EQUILIBRIUM_HIDDEN_MARKOV_IND_OUT_TREE ,
  TREESTATW_MARKOV_OUT_TREE ,
  TREESTATW_GENERATION_DISTRIBUTION ,
  TREESTATW_GENERATION ,
  TREESTATW_GENERATION_PROCESS ,
  TREESTATW_GENERATION_PROCESSES ,
  TREESTATW_PARENT_STATE ,
  TREESTATW_FACTOR ,
  TREESTATW_EXTERNAL_FACTOR ,
  TREESTATW_FACTORS ,
  TREESTATW_EXTERNAL_FACTORS ,
  TREESTATW_CHILD ,
  TREESTATW_CHILDREN ,
  TREESTATW_ORDERED_CHILD ,
  TREESTATW_ORDERED_CHILDREN ,
  TREESTATW_VALUE ,
  TREESTATW_VARIABLE_ORDER_MARKOV_CHAIN ,
  TREESTATW_VARIABLE_ORDER_MARKOV_CHAINS
};

extern const char *STAT_TREES_word[];
extern const char *STAT_TREES_type[];

/****************************************************************
 *
 *  Multivariate distribution identifiers:
 */

extern const char *STAT_TREES_multivariate_distribution_word[];

/****************************************************************
 *
 *  Labels:
 */

enum {
  TREESTATL_DISCRETE_MULTIVARIATE_DISTRIBUTION
};

extern const char *STAT_MULTIVARIATE_label[];

/****************************************************************
 *
 *  Keywords:
 */

enum {
  TREESTATW_COMPOUNDING_DISTRIBUTION
};


extern const char *STAT_MULTIVARIATE_PARSING_word[];




enum {
  TREESTATL_STATE_TREE_LIKELIHOOD ,
  TREESTATL_STATE_TREES_LIKELIHOOD ,
  TREESTATL_OBSERVED_TREES_LIKELIHOOD ,
  TREESTATL_TREES_IID_INFORMATION ,

  TREESTATL_SMOOTHED ,
  TREESTATL_OBSERVED ,

  TREESTATL_STATE_NO_OCCURRENCE ,
  TREESTATL_OUTPUT_NO_OCCURRENCE ,
  TREESTATL_STATE_ABSORPTION ,
  TREESTATL_OUTPUT_ABSORPTION ,

  TREESTATL_STATE_FIRST_OCCURRENCE_ROOT ,
  TREESTATL_OUTPUT_FIRST_OCCURRENCE_ROOT ,
  TREESTATL_STATE_FIRST_OCCURRENCE_LEAVES ,
  TREESTATL_OUTPUT_FIRST_OCCURRENCE_LEAVES ,
  TREESTATL_SOJOURN_SIZE ,
  TREESTATL_MIXTURE_OF ,
  TREESTATL_STATE_NB_ZONES ,
  TREESTATL_OUTPUT_NB_ZONES ,
  TREESTATL_STATE_NB_OCCURRENCES ,
  TREESTATL_OUTPUT_NB_OCCURRENCES ,
  TREESTATL_PER_TREE ,
  TREESTATL_PER_SIZE ,

  TREESTATL_STATE_PROBABILITY ,
  TREESTATL_OBSERVATION_DISTRIBUTION_DISTANCE ,
  TREESTATL_POSTERIOR_STATE_PROBABILITY ,
  TREESTATL_POSTERIOR_IN_STATE_PROBABILITY ,
  TREESTATL_POSTERIOR_OUT_STATE_PROBABILITY ,
  TREESTATL_CONDITIONAL_ENTROPY ,
  TREESTATL_MARGINAL_ENTROPY ,
  TREESTATL_MARGINAL_ENTROPY_SUM ,
  TREESTATL_PARTIAL_STATE_TREE_ENTROPY ,
  TREESTATL_STATE_TREE_ENTROPY ,
  TREESTATL_GINI_INDEX ,
  TREESTATL_UPPER_BOUND ,
  TREESTATL_NB_STATE_TREE ,
  TREESTATL_MAX_POSTERIOR_STATE_PROBABILITY ,
  TREESTATL_MAX_POSTERIOR_IN_STATE_PROBABILITY ,
  TREESTATL_MAX_POSTERIOR_OUT_STATE_PROBABILITY ,
  TREESTATL_AMBIGUITY ,

  TREESTATL_TREE ,
  TREESTATL_TREES ,
  TREESTATL_VERTEX ,

  TREESTATL_TREE_SIZE ,
  TREESTATL_CUMULATIVE_SIZE ,
  TREESTATL_TREE_CHILDREN ,
  TREESTATL_CUMULATIVE_CHILDREN ,
  TREESTATL_TIME ,

  TREESTATL_GENERATION ,
  TREESTATL_GENERATION_PROCESS ,
  TREESTATL_PARENT
};


extern const char *STAT_TREES_label[];


/****************************************************************
 *
 *  Identifiers for parsing error messages:
 */

enum {
  TREESTATP_BAD_PROBABILITY ,
  TREESTATP_BAD_ONE_MINUS_PROBABILITIES ,
  TREESTATP_ORDERED_CHILDREN ,
  TREESTATP_NB_VOMC ,
  TREESTATP_NB_GENERATION_PROCESS ,
  TREESTATP_NB_CHILDREN_BRANCHING ,
  TREESTATP_NB_FACTORS ,
  TREESTATP_FACTOR_VALUE ,
  TREESTATP_NB_FACTOR_VALUES
};

extern const char *STAT_TREES_parsing[];

/****************************************************************
 *
 *  Identifiers for other error messages:
 */


enum {
  TREESTATR_NB_STATE ,
  TREESTATR_VARIABLE_NB_VALUE ,
//  TREESTATR_NB_OUTPUT_PROCESS ,
  TREESTATR_NB_INT_OUTPUT_PROCESS ,
  TREESTATR_NB_REAL_OUTPUT_PROCESS ,
  TREESTATR_OUTPUT_PROCESS_INDEX ,
//  TREESTATR_NB_OUTPUT ,
  TREESTATR_SELF_TRANSITION ,
  TREESTATR_NB_STATE_TREES ,
  TREESTATR_NB_TREES ,
  TREESTATR_TREE_IDENTIFIER ,
  TREESTATR_WIDE_TREE ,
  TREESTATR_NARROW_TREE ,
  TREESTATR_TREE_NB_CHILDREN ,
  TREESTATR_TREE_SIZE ,
  TREESTATR_VERTEX_ID ,
  TREESTATR_CHILD_ID ,
  TREESTATR_SMALL_TREE_SIZE ,
  TREESTATR_BIG_TREE_SIZE ,
  TREESTATR_TREE_CUMULATIVE_SIZE ,
  TREESTATR_TREE_CUMULATIVE_CHILDREN ,
  TREESTATR_STATE_TREE_COMPUTATION_FAILURE ,

  TREESTATR_NOT_PRESENT ,

  TREESTATR_NO_MODEL ,

  TREESTATR_NO_PERMUTATION ,

  TREESTATR_SAEM_EXP ,
  TREESTATR_EM_ALGORITHM ,

  TREESTATR_VARIABLE_TYPE ,
  TREESTATR_VARIABLE_1_TYPE ,
  TREESTATR_MIN_TREE_SIZE ,
  TREESTATR_MAX_TREE_SIZE ,
  TREESTATR_NB_SELECTED_VALUE ,

  TREESTATR_STATE_TREES ,
  TREESTATR_CHARACTERISTICS_NOT_COMPUTED ,
  TREESTATR_NON_EXISTING_CHARACTERISTIC_DISTRIBUTION ,
  TREESTATR_BAD_FACTORS ,
  TREESTATR_NON_EXISTING_GENERATION_PROCESS
};


extern const char *STAT_TREES_error[];

};// end namespace

#endif
