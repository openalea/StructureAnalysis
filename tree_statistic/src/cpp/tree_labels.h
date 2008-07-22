/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       AMAPmod: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2002 UMR Cirad/Inra Modelisation des Plantes
 *
 *       File author(s): Y. Guedon (yann.guedon@cirad.fr)
 *
 *       $Source: /usr/cvsmaster/AMAPmod/src/STAT_TREES/src/tree_labels.h,v $
 *       $Id: tree_labels.h 3186 2007-05-25 15:10:30Z dufourko $
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



#ifndef TREE_LABELS_H
#define TREE_LABELS_H

namespace Stat_trees
{

/****************************************************************
 *
 *  Identificateurs des mots cles (format des fichiers) :
 */


enum {
  STATW_HIDDEN_MARKOV_TREE ,
  STATW_EQUILIBRIUM_HIDDEN_MARKOV_TREE ,
  STATW_HIDDEN_MARKOV_OUT_TREE ,
  STATW_EQUILIBRIUM_HIDDEN_MARKOV_OUT_TREE
};

extern const char *STAT_TREES_word[];
extern const char *STAT_TREES_type[];


/****************************************************************
 *
 *  Identificateurs des labels :
 */


enum {
  STATL_STATE_TREE_LIKELIHOOD ,
  STATL_STATE_TREES_LIKELIHOOD ,
  STATL_OBSERVED_TREES_LIKELIHOOD ,
  STATL_TREES_IID_INFORMATION ,

  STATL_SMOOTHED ,
  STATL_OBSERVED ,

  STATL_STATE_NO_OCCURRENCE ,
  STATL_OUTPUT_NO_OCCURRENCE ,
  STATL_STATE_ABSORPTION ,
  STATL_OUTPUT_ABSORPTION ,

  STATL_STATE_FIRST_OCCURRENCE_ROOT ,
  STATL_OUTPUT_FIRST_OCCURRENCE_ROOT ,
  STATL_STATE_FIRST_OCCURRENCE_LEAVES ,
  STATL_OUTPUT_FIRST_OCCURRENCE_LEAVES ,
  STATL_SOJOURN_SIZE ,
  STATL_MIXTURE_OF ,
  STATL_STATE_NB_ZONES ,
  STATL_OUTPUT_NB_ZONES ,
  STATL_STATE_NB_OCCURRENCES ,
  STATL_OUTPUT_NB_OCCURRENCES ,
  STATL_PER_TREE ,
  STATL_PER_SIZE ,

  STATL_STATE_PROBABILITY ,
  STATL_POSTERIOR_STATE_PROBABILITY ,
  STATL_POSTERIOR_IN_STATE_PROBABILITY ,
  STATL_POSTERIOR_OUT_STATE_PROBABILITY ,
  STATL_CONDITIONAL_ENTROPY ,
  STATL_MARGINAL_ENTROPY ,
  STATL_MARGINAL_ENTROPY_SUM ,
  STATL_PARTIAL_STATE_TREE_ENTROPY ,
  STATL_STATE_TREE_ENTROPY ,
  STATL_GINI_INDEX ,
  STATL_UPPER_BOUND ,
  STATL_NB_STATE_TREE ,
  STATL_MAX_POSTERIOR_STATE_PROBABILITY ,
  STATL_MAX_POSTERIOR_IN_STATE_PROBABILITY ,
  STATL_MAX_POSTERIOR_OUT_STATE_PROBABILITY ,
  STATL_AMBIGUITY ,

  STATL_TREE ,
  STATL_TREES ,
  STATL_VERTEX ,

  STATL_TREE_SIZE ,
  STATL_CUMULATIVE_SIZE ,
  STATL_TREE_CHILDREN ,
  STATL_CUMULATIVE_CHILDREN ,
  STATL_TIME
};


extern const char *STAT_TREES_label[];


/****************************************************************
 *
 *  Identificateurs des messages d'erreur pour l'analyse des fichiers :
 */



/****************************************************************
 *
 *  Identificateurs des messages d'erreur de traitement :
 */


enum {
  STATR_NB_STATE ,
  STATR_VARIABLE_NB_VALUE ,
  STATR_NB_OUTPUT_PROCESS ,
  STATR_NB_INT_OUTPUT_PROCESS ,
  STATR_NB_REAL_OUTPUT_PROCESS ,
  STATR_OUTPUT_PROCESS_INDEX ,
  STATR_NB_OUTPUT ,
  STATR_SELF_TRANSITION ,
  STATR_NB_STATE_TREES ,
  STATR_NB_TREES ,
  STATR_TREE_IDENTIFIER ,
  STATR_WIDE_TREE ,
  STATR_NARROW_TREE ,
  STATR_TREE_NB_CHILDREN ,
  STATR_TREE_SIZE ,
  STATR_VERTEX_ID ,
  STATR_CHILD_ID ,
  STATR_SMALL_TREE_SIZE ,
  STATR_BIG_TREE_SIZE ,
  STATR_TREE_CUMULATIVE_SIZE ,
  STATR_TREE_CUMULATIVE_CHILDREN ,
  STATR_STATE_TREE_COMPUTATION_FAILURE ,

  STATR_NOT_PRESENT ,

  STATR_NO_PERMUTATION ,

  STATR_SAEM_EXP ,
  STATR_EM_ALGORITHM ,

  STATR_VARIABLE_TYPE ,
  STATR_VARIABLE_1_TYPE ,
  STATR_MIN_TREE_SIZE ,
  STATR_MAX_TREE_SIZE ,
  STATR_NB_SELECTED_VALUE ,

  STATR_STATE_TREES ,
  STATR_CHARACTERISTICS_NOT_COMPUTED ,
  STATR_NON_EXISTING_CHARACTERISTIC_DISTRIBUTION
};


extern const char *STAT_TREES_error[];

};// end namespace

#endif
