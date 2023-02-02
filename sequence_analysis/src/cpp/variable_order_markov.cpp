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



#include <math.h>

#include <string>
#include <sstream>
#include <iomanip>

#include <boost/tokenizer.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/math/distributions/normal.hpp>

#include "stat_tool/stat_label.h"

#include "variable_order_markov.h"
#include "sequence_label.h"

using namespace std;
using namespace boost;
using namespace boost::math;
using namespace stat_tool;


namespace sequence_analysis {



/*--------------------------------------------------------------*/
/**
 *  \brief Default constructor of the VariableOrderMarkovChain class.
 */
/*--------------------------------------------------------------*/

VariableOrderMarkovChain::VariableOrderMarkovChain()

{
  max_order = 0;

  memo_type = NULL;
  order = NULL;
  state = NULL;

  parent = NULL;
  child = NULL;

  next = NULL;
  nb_memory = NULL;
  previous = NULL;

  state_process = NULL;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Constructor of the VariableOrderMarkovChain class.
 *
 *  \param[in] itype     process type (ORDINARY/EQUILIBRIUM),
 *  \param[in] inb_state number of states,
 *  \param[in] inb_row   number of memories.
 */
/*--------------------------------------------------------------*/

VariableOrderMarkovChain::VariableOrderMarkovChain(process_type itype , int inb_state , int inb_row)
:Chain(itype , inb_state , inb_row , true)

{
  int i;

  max_order = 0;

  memo_type = new memory_type[nb_row];
  order = new int[nb_row];
  state = new int*[nb_row];
  parent = new int[nb_row];
  child = new int*[nb_row];

  for (i = 0;i < nb_row;i++) {
    state[i] = NULL;
    child[i] = NULL;
  }

  next = NULL;
  nb_memory = NULL;
  previous = NULL;

  state_process = new CategoricalSequenceProcess(nb_state , nb_state);
}


/*--------------------------------------------------------------*/
/**
 *  \brief Constructor of the VariableOrderMarkovChain class.
 *
 *  \param[in] itype      process type (ORDINARY/EQUILIBRIUM),
 *  \param[in] inb_state  number of states,
 *  \param[in] inb_row    number of memories,
 *  \param[in] imax_order maximum order.
 */
/*--------------------------------------------------------------*/

VariableOrderMarkovChain::VariableOrderMarkovChain(process_type itype , int inb_state ,
                                                   int inb_row , int imax_order)
:Chain(itype , inb_state , inb_row , true)

{
  int i;


  max_order = imax_order;

  memo_type = new memory_type[nb_row];
  order = new int[nb_row];
  state = new int*[nb_row];
  parent = new int[nb_row];
  child = new int*[nb_row];

  for (i = 0;i < nb_row;i++) {
    state[i] = new int[max_order];
    child[i] = NULL;
  }

  next = NULL;
  nb_memory = NULL;
  previous = NULL;

  state_process = new CategoricalSequenceProcess(nb_state , nb_state);
}


/*--------------------------------------------------------------*/
/**
 *  \brief Construction of a VariableOrderMarkovChain object of fixed order.
 *
 *  \param[in] itype     process type (ORDINARY/EQUILIBRIUM),
 *  \param[in] inb_state number of states,
 *  \param[in] iorder    order,
 *  \param[in] init_flag flag initialization.
 */
/*--------------------------------------------------------------*/

VariableOrderMarkovChain::VariableOrderMarkovChain(process_type itype , int inb_state ,
                                                   int iorder , bool init_flag)
:Chain(itype , inb_state , (int)(pow((double)inb_state , iorder + 1) - 1) / (inb_state - 1) , init_flag)

{
  int i , j;

  max_order = iorder;

  memo_type = new memory_type[nb_row];
  order = new int[nb_row];
  state = new int*[nb_row];
  parent = new int[nb_row];
  child = new int*[nb_row];

  // root (zero order)

  memo_type[0] = NON_TERMINAL;
  order[0] = 0;
  state[0] = NULL;
  parent[0] = -1;
  child[0] = new int[nb_state];

  for (i = 1;i < nb_row;i++) {

    // case increase of the order

    if (order[i - 1] < max_order) {
      order[i] = order[i - 1] + 1;
      if (order[i] < max_order) {
        memo_type[i] = NON_TERMINAL;
      }
      else {
        memo_type[i] = TERMINAL;
      }

      state[i] = new int[order[i]];
      for (j = 0;j < order[i - 1];j++) {
        state[i][j] = state[i - 1][j];
      }
      state[i][order[i] - 1] = 0;

      parent[i] = i - 1;
      child[i - 1][0] = i;
    }

    else {

      // case stable (maximum) order

      if (state[i - 1][order[i - 1] - 1] < nb_state - 1) {
        memo_type[i] = TERMINAL;
        order[i] = max_order;
        state[i] = new int[order[i]];
        for (j = 0;j < order[i] - 1;j++) {
          state[i][j] = state[i - 1][j];
        }
        state[i][order[i] - 1] = state[i - 1][order[i] - 1] + 1;

        parent[i] = i - state[i][order[i] - 1] - 1;
        child[parent[i]][state[i][order[i] - 1]] = i;
      }

      // case decrease of the order

      else {
        memo_type[i] = NON_TERMINAL;

        for (j = order[i - 1] - 2;j >= 0;j--) {
          if (state[i - 1][j] != nb_state - 1) {
            break;
          }
        }
        order[i] = j + 1;

        state[i] = new int[order[i]];
        for (j = 0;j < order[i] - 1;j++) {
          state[i][j] = state[i - 1][j];
        }
        state[i][order[i] - 1] = state[i - 1][order[i] - 1] + 1;

        // search for the parent vertex

        find_parent_memory(i);
      }
    }

    if (memo_type[i] == NON_TERMINAL) {
      child[i] = new int[nb_state];
    }
    else {
      child[i] = NULL;
    }
  }

  // computation of the transitions between terminal memories

  next = NULL;
  nb_memory = NULL;
  previous = NULL;

  build_memory_transition();

  state_process = new CategoricalSequenceProcess(nb_state , nb_state , false);
}


/*--------------------------------------------------------------*/
/**
 *  \brief Completion of the memory tree.
 *
 *  \param[in] markov reference on a VariableOrderMarkovChain object.
 */
/*--------------------------------------------------------------*/

void VariableOrderMarkovChain::memory_tree_completion(const VariableOrderMarkovChain &markov)

{
  bool prefix;
  int i , j , k , m;
  int bnb_memory , border , *markov_next , *completion_next;
  VariableOrderMarkovChain *completion;


  bnb_memory = markov.nb_row;
  for (i = 1;i < markov.nb_row;i++) {
    if (markov.order[i] == markov.max_order) {
      bnb_memory--;
    }
  }

  completion = new VariableOrderMarkovChain(markov.type , markov.nb_state ,
                                            (int)(pow((double)markov.nb_state , markov.max_order) - 1) / (markov.nb_state - 1) - bnb_memory);
  completion->nb_row = 0;

  for (i = 1;i < markov.nb_row;i++) {
    if (markov.order[i] > 1) {

      // search for the vertex corresponding to the longer proper prefix

      prefix = false;

      for (j = 0;j < markov.nb_row;j++) {
        if (markov.order[j] == markov.order[i] - 1) {
          for (k = 0;k < markov.order[j];k++) {
            if (markov.state[j][k] != markov.state[i][k + 1]) {
              break;
            }
          }

          if (k == markov.order[j]) {
            prefix = true;
            break;
          }
        }
      }

      if (!prefix) {
        for (j = 0;j < completion->nb_row;j++) {
          if (completion->order[j] == markov.order[i] - 1) {
            for (k = 0;k < completion->order[j];k++) {
              if (completion->state[j][k] != markov.state[i][k + 1]) {
                break;
              }
            }

            if (k == completion->order[j]) {
              prefix = true;
              break;
            }
          }
        }

        if (!prefix) {

          // construction of the vertex corresponding to the longer proper prefix

          completion->order[completion->nb_row] = markov.order[i] - 1;
          completion->state[completion->nb_row] = new int[completion->order[completion->nb_row]];
          for (j = 0;j < completion->order[completion->nb_row];j++) {
            completion->state[completion->nb_row][j] = markov.state[i][j + 1];
          }
          for (j = 0;j < markov.nb_state;j++) {
            completion->transition[completion->nb_row][j] = markov.transition[i][j];
          }
          (completion->nb_row)++;

          for (j = completion->order[completion->nb_row - 1] - 1;j >= 2;j--) {
            prefix = false;

            for (k = 0;k < markov.nb_row;k++) {
              if (markov.order[k] == j) {
                for (m = 0;m < markov.order[k];m++) {
                  if (markov.state[k][m] != completion->state[completion->nb_row - 1][m + 1]) {
                    break;
                  }
                }

                if (m == markov.order[k]) {
                  prefix = true;
                  break;
                }
              }
            }

            if (!prefix) {
              for (k = 0;k < completion->nb_row - 1;k++) {
                if (completion->order[k] == j) {
                  for (m = 0;m < completion->order[k];m++) {
                    if (completion->state[k][m] != completion->state[completion->nb_row - 1][m + 1]) {
                      break;
                    }
                  }

                  if (m == completion->order[k]) {
                    prefix = true;
                    break;
                  }
                }
              }

              if (!prefix) {

                // construction of the vertex corresponding of the longer proper prefix

                completion->order[completion->nb_row] = j;
                completion->state[completion->nb_row] = new int[completion->order[completion->nb_row]];
                for (k = 0;k < completion->order[completion->nb_row];k++) {
                  completion->state[completion->nb_row][k] = completion->state[completion->nb_row - 1][k + 1];
                }
                for (k = 0;k < markov.nb_state;k++) {
                  completion->transition[completion->nb_row][k] = completion->transition[completion->nb_row - 1][k];
                }
                (completion->nb_row)++;
              }

              else {
                break;
              }
            }
          }
        }
      }
    }
  }

# ifdef DEBUG
  {
    cout << "\n";
    for (i = 0;i < completion->nb_row;i++) {
      cout << completion->order[i] << " | ";
      for (j = completion->order[i] - 1;j >= 0;j--) {
        cout << completion->state[i][j] << " ";
      }
      cout << endl;
    }
  }
# endif

  // search for the following vertex in the ordering of the completed memory tree 

  markov_next = new int[markov.nb_row];
  completion_next = new int[completion->nb_row];

  for (i = 0;i < markov.nb_row;i++) {
    markov_next[i] = I_DEFAULT;

    if ((markov.memo_type[i] == TERMINAL) && (markov.order[i] < markov.max_order - 1)) {

      // search for the 1st child (state 0) of the terminal vertex

      for (j = 0;j < completion->nb_row;j++) {
        if (completion->order[j] == markov.order[i] + 1) {
          for (k = 0;k < markov.order[i];k++) {
            if (completion->state[j][k] != markov.state[i][k]) {
              break;
            }
          }

          if ((k == markov.order[i]) && (completion->state[j][completion->order[j] - 1] == 0)) {
            markov_next[i] = j;
            break;
          }
        }
      }
    }
  }

  for (i = 0;i < completion->nb_row;i++) {
    completion_next[i] = I_DEFAULT;

    // search for the 1st child (state 0) of the built vertex

    if (completion->order[i] < markov.max_order - 1) {
      for (j = 0;j < completion->nb_row;j++) {
        if (completion->order[j] == completion->order[i] + 1) {
          for (k = 0;k < completion->order[i];k++) {
            if (completion->state[j][k] != completion->state[i][k]) {
              break;
            }
          }

          if ((k == completion->order[i]) && (completion->state[j][completion->order[j] - 1] == 0)) {
            completion_next[i] = j;
            break;
          }
        }
      }
    }

    if ((completion->order[i] == markov.max_order - 1) || (j == completion->nb_row)) {
      if (completion->state[i][completion->order[i] - 1] < markov.nb_state - 1) {

        // search for siblings (following states) of the built vertex

        for (j = 0;j < completion->nb_row;j++) {
          if (completion->order[j] == completion->order[i]) {
            for (k = 0;k < completion->order[i] - 1;k++) {
              if (completion->state[j][k] != completion->state[i][k]) {
                break;
              }
            }

            if ((k == completion->order[i] - 1) &&
                (completion->state[j][completion->order[j] - 1] == completion->state[i][completion->order[i] - 1] + 1)) {
              completion_next[i] = j;
              break;
            }
          }
        }
      }

      else {
        for (j = completion->order[i] - 2;j >= 0;j--) {
          if (completion->state[i][j] != markov.nb_state - 1) {
            break;
          }
        }

        border = j + 1;
        for (j = 0;j < completion->nb_row;j++) {
          if (completion->order[j] == border) {
            for (k = 0;k < border - 1;k++) {
              if (completion->state[j][k] != completion->state[i][k]) {
                break;
              }
            }

            if ((k == border - 1) &&
                (completion->state[j][border - 1] == completion->state[i][border - 1] + 1)) {
              completion_next[i] = j;
              break;
            }
          }
        }
      }
    }
  }

# ifdef DEBUG
  {
    cout << "\n";
    for (i = 0;i < markov.nb_row;i++) {
      for (j = markov.order[i] - 1;j >= 0;j--) {
        cout << markov.state[i][j] << " ";
      }
      if (markov_next[i] != I_DEFAULT) {
        cout << " -> ";
        for (j = completion->order[markov_next[i]] - 1;j >= 0;j--) {
          cout << completion->state[markov_next[i]][j] << " ";
        }
      }
      cout << endl;
    }

    cout << "\n";
    for (i = 0;i < completion->nb_row;i++) {
      for (j = completion->order[i] - 1;j >= 0;j--) {
        cout << completion->state[i][j] << " ";
      }
      if (completion_next[i] != I_DEFAULT) {
        cout << " -> ";
        for (j = completion->order[completion_next[i]] - 1;j >= 0;j--) {
          cout << completion->state[completion_next[i]][j] << " ";
        }
      }
      cout << endl;
    }
  }
# endif

  // copy of parameters

  type = markov.type;
  nb_state = markov.nb_state;
  nb_row = markov.nb_row + completion->nb_row;

  if (markov.nb_component > 0) {
    accessibility = new bool*[nb_state];
    for (i = 0;i < nb_state;i++) {
      accessibility[i] = new bool[nb_state];
      for (j = 0;j < nb_state;j++) {
        accessibility[i][j] = markov.accessibility[i][j];
      }
    }

    nb_component = markov.nb_component;
    component_nb_state = new int[nb_component];
    component = new int*[nb_component];

    for (i = 0;i < nb_component;i++) {
      component_nb_state[i] = markov.component_nb_state[i];
      component[i] = new int[component_nb_state[i]];
      for (j = 0;j < component_nb_state[i];j++) {
        component[i][j] = markov.component[i][j];
      }
    }

    stype = new state_type[nb_state];
    for (i = 0;i < nb_state;i++) {
      stype[i] = markov.stype[i];
    }
  }

  else {
    accessibility = NULL;
    nb_component = 0;
    component_nb_state = NULL;
    component = NULL;
    stype = NULL;
  }

  initial = new double[type == ORDINARY ? nb_state : nb_row];

  if (type == ORDINARY) {
    for (i = 0;i < nb_state;i++) {
      initial[i] = markov.initial[i];
    }
  }

  transition = new double*[nb_row];
  for (i = 0;i < nb_row;i++) {
    transition[i] = new double[nb_state];
  }

  cumul_initial = NULL;
  cumul_transition = NULL;

  max_order = markov.max_order;

  memo_type = new memory_type[nb_row];
  order = new int[nb_row];
  state = new int*[nb_row];
  parent = new int[nb_row];
  child = new int*[nb_row];

  next = NULL;
  nb_memory = NULL;
  previous = NULL;

  // insertion of the memories in the out-tree

  i = 0;
  for (j = 0;j < markov.nb_row;j++) {
    for (k = 0;k < nb_state;k++) {
      transition[i][k] = markov.transition[j][k];
    }

    memo_type[i] = markov.memo_type[j];

    order[i] = markov.order[j];
    state[i] = new int[order[i]];
    for (k = 0;k < order[i];k++) {
      state[i][k] = markov.state[j][k];
    }

    parent[i] = markov.parent[j];

    if ((markov.memo_type[j] == NON_TERMINAL) || (markov_next[j] != I_DEFAULT)) {
      child[i] = new int[nb_state];

      if (markov.memo_type[j] == NON_TERMINAL) {
        for (k = 0;k < nb_state;k++) {
          child[i][k] = markov.child[j][k];
        }
      }
    }

    else {
      child[i] = NULL;
    }

    i++;

    if (markov_next[j] != I_DEFAULT) {
      k = markov_next[j];
      bnb_memory = i;

      do {
        memo_type[i] = COMPLETION;

        order[i] = completion->order[k];
        state[i] = new int[order[i]];
        for (m = 0;m < order[i];m++) {
          state[i][m] = completion->state[k][m];
        }

        find_parent_memory(i);

        for (m = 0;m < nb_state;m++) {
          transition[i][m] = transition[parent[i]][m];
        }

        if ((completion_next[k] != I_DEFAULT) &&
            (completion->order[completion_next[k]] == completion->order[k] + 1)) {
          child[i] = new int[nb_state];
        }
        else {
          child[i] = NULL;
        }

        i++;
        k = completion_next[k];
      }
      while (k != I_DEFAULT);

      // update of the relationships parent/children

      for (k = 0;k < bnb_memory;k++) {
        if (memo_type[k] == NON_TERMINAL) {
          for (m = 0;m < nb_state;m++) {
            if (child[k][m] >= bnb_memory) {
              child[k][m] += i - bnb_memory;
            }
          }
        }
      }

      for (k = j + 1;k < markov.nb_row;k++) {
        if (markov.parent[k] >= bnb_memory) {
          markov.parent[k] += i - bnb_memory;
        }
        if (markov.memo_type[k] == NON_TERMINAL) {
          for (m = 0;m < nb_state;m++) {
            if (markov.child[k][m] >= bnb_memory) {
              markov.child[k][m] += i - bnb_memory;
            }
          }
        }
      }
    }
  }

# ifdef DEBUG
  {
    cout << "\n";
    cout << "Suffix free? " << (check_free_suffix() ? "True" : "False") << endl;

    for (i = 0;i < nb_row;i++) {
      cout << i << "  ";
      for (j = max_order - 1;j >= order[i];j--) {
        cout << "  ";
      }
      for (j = order[i] - 1;j >= 0;j--) {
        cout << state[i][j] << " ";
      }

      switch (memo_type[i]) {
      case NON_TERMINAL :
        cout << "  " << SEQ_label[SEQL_NON_TERMINAL];
        break;
      case TERMINAL :
        cout << "      " << SEQ_label[SEQL_TERMINAL];
        break;
      case COMPLETION :
        cout << "    " << SEQ_label[SEQL_COMPLETION];
        break;
      }

      cout << " | " << parent[i];

      if (child[i]) {
        cout << " |";
        for (j = 0;j < nb_state;j++) {
          cout << " " << child[i][j];
        }
      }

      else {
        cout << "  ";
        for (j = 0;j < nb_state;j++) {
          cout << "  ";
        }
      }

      cout << endl;
    }
  }
# endif

  delete completion;

  delete [] markov_next;
  delete [] completion_next;

  build_memory_transition();
  build_previous_memory();
}


/*--------------------------------------------------------------*/
/**
 *  \brief Construction of a VariableOrderMarkovChain object with completion of
 *         the memory tree, computation of the transition distributions for the non-terminal memories
 *         (ordinary process) or computation of the stationary distribution (equilibrium process).
 *
 *  \param[in] markov reference on a VariableOrderMarkovChain object.
 */
/*--------------------------------------------------------------*/

void VariableOrderMarkovChain::build(const VariableOrderMarkovChain &markov)

{
  int i;
  int nb_terminal;


  memory_tree_completion(markov);

  switch (type) {

  case ORDINARY : {
    non_terminal_transition_probability_computation();
    break;
  }

  case EQUILIBRIUM : {
    nb_terminal = (nb_row - 1) * (nb_state - 1) / nb_state + 1;

    for (i = 1;i < nb_row;i++) {
      if (!child[i]) {
        initial[i] = 1. / (double)nb_terminal;
      }
      else {
        initial[i] = 0.;
      }
    }

    initial_probability_computation();
    break;
  }
  }

  state_process = new CategoricalSequenceProcess(nb_state , nb_state);
}


/*--------------------------------------------------------------*/
/**
 *  \brief Copy of a VariableOrderMarkovChain object.
 *
 *  \param[in] markov reference on a VariableOrderMarkovChain object.
 */
/*--------------------------------------------------------------*/

void VariableOrderMarkovChain::copy(const VariableOrderMarkovChain &markov)

{
  int i , j;


  memo_type = new memory_type[nb_row];
  for (i = 0;i < nb_row;i++) {
    memo_type[i] = markov.memo_type[i];
  }

  order = new int[nb_row];
  for (i = 0;i < nb_row;i++) {
    order[i] = markov.order[i];
  }

  max_order = markov.max_order;

  state = new int*[nb_row];
  for (i = 0;i < nb_row;i++) {
    state[i] = new int[order[i]];
    for (j = 0;j < order[i];j++) {
      state[i][j] = markov.state[i][j];
    }
  }

  parent = new int[nb_row];
  for (i = 0;i < nb_row;i++) {
    parent[i] = markov.parent[i];
  }

  child = new int*[nb_row];
  for (i = 0;i < nb_row;i++) {
    if (markov.child[i]) {
      child[i] = new int[nb_state];
      for (j = 0;j < nb_state;j++) {
        child[i][j] = markov.child[i][j];
      }
    }
    else {
      child[i] = NULL;
    }
  }

  if (markov.next) {
    next = new int*[nb_row];
    for (i = 0;i < nb_row;i++) {
      if (markov.next[i]) {
        next[i] = new int[nb_state];
        for (j = 0;j < nb_state;j++) {
          next[i][j] = markov.next[i][j];
        }
      }
      else {
        next[i] = NULL;
      }
    }
  }
  else {
    next = NULL;
  }

  if (markov.nb_memory) {
    nb_memory = new int[nb_row];
    for (i = 0;i < nb_row;i++) {
      nb_memory[i] = markov.nb_memory[i];
    }
  }
  else {
    nb_memory = NULL;
  }

  if (markov.previous) {
    previous = new int*[nb_row];
    for (i = 0;i < nb_row;i++) {
      if (markov.previous[i]) {
        previous[i] = new int[nb_memory[i]];
        for (j = 0;j < nb_memory[i];j++) {
          previous[i][j] = markov.previous[i][j];
        }
      }
      else {
        previous[i] = NULL;
      }
    }
  }
  else {
    previous = NULL;
  }

  state_process = new CategoricalSequenceProcess(*(markov.state_process));
}


/*--------------------------------------------------------------*/
/**
 *  \brief Destruction of the data members of a VariableOrderMarkovChain object.
 */
/*--------------------------------------------------------------*/

void VariableOrderMarkovChain::remove()

{
  int i;


  delete [] memo_type;
  delete [] order;

  if (state) {
    for (i = 0;i < nb_row;i++) {
      delete [] state[i];
    }
    delete [] state;
  }

  delete [] parent;

  if (child) {
    for (i = 0;i < nb_row;i++) {
      delete [] child[i];
    }
    delete [] child;
  }

  if (next) {
    for (i = 1;i < nb_row;i++) {
      delete [] next[i];
    }
    delete [] next;
  }

  delete [] nb_memory;

  if (previous) {
    for (i = 1;i < nb_row;i++) {
      delete [] previous[i];
    }
    delete [] previous;
  }

  delete state_process;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Destructor of the VariableOrderMarkovChain class.
 */
/*--------------------------------------------------------------*/

VariableOrderMarkovChain::~VariableOrderMarkovChain()

{
  remove();
}


/*--------------------------------------------------------------*/
/**
 *  \brief Assignment operator of the VariableOrderMarkovChain class.
 *
 *  \param[in] markov reference on a VariableOrderMarkovChain object.
 *
 *  \return           VariableOrderMarkovChain object.
 */
/*--------------------------------------------------------------*/

VariableOrderMarkovChain& VariableOrderMarkovChain::operator=(const VariableOrderMarkovChain &markov)

{
  if (&markov != this) {
    remove();
    Chain::remove();

    Chain::copy(markov);
    copy(markov);
  }

  return *this;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Search for the parent memory.
 *
 *  \param[in] index memory index.
 */
/*--------------------------------------------------------------*/

void VariableOrderMarkovChain::find_parent_memory(int index)

{
  int i;


  for (i = index - 1;i >= 0;i--) {
    if (order[i] == order[index] - 1) {
      parent[index] = i;
      child[i][state[index][order[index] - 1]] = index;
      break;
    }
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Construction of the transitions between memories.
 */
/*--------------------------------------------------------------*/

void VariableOrderMarkovChain::build_memory_transition()

{
  if (!next) {
    int i , j , k;
    int bnb_memory;


#   ifdef DEBUG
    cout << "\n";
#   endif

    next = new int*[nb_row];
    next[0] = NULL;

    for (i = 1;i < nb_row;i++) {
      if ((type == ORDINARY) || (!child[i])) {
        next[i] = new int[nb_state];

#       ifdef DEBUG
        for (j = order[i] - 1;j >= 0;j--) {
          cout << state[i][j] << " ";
        }
#       endif

        bnb_memory = 0;
        for (j = 1;j < nb_row;j++) {
          if (((child[i]) && (child[j]) && (order[j] == order[i] + 1)) ||
              ((!child[j]) && (order[j] <= order[i] + 1))) {
            for (k = 0;k < order[j] - 1;k++) {
              if (state[j][k + 1] != state[i][k]) {
                break;
              }
            }

            if ((order[j] == 1) || (k == order[j] - 1)) {
//              if ((memo_type[i] == NON_TERMINAL) || (transition[i][state[j][0]] > 0.)) {
                next[i][state[j][0]] = j;
/*              }
              else {
                next[i][state[j][0]] = I_DEFAULT;
              } */

#             ifdef DEBUG
              cout << "| ";
              for (k = order[j] - 1;k >= 0;k--) {
                cout << state[j][k] << " ";
              }
#             endif

              bnb_memory++;
              if (bnb_memory == nb_state) {
                break;
              }
            }
          }
        }

#       ifdef DEBUG
        cout << endl;
#       endif

      }

      else {
        next[i] = NULL;
      }
    }
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Extraction of the previous memories.
 */
/*--------------------------------------------------------------*/

void VariableOrderMarkovChain::build_previous_memory()

{
  if ((next) && (!nb_memory) && (!previous)) {
    int i , j;
    int *buffer;


    nb_memory = new int[nb_row];
    previous = new int*[nb_row];
    nb_memory[0] = 0;
    previous[0] = NULL;

    buffer = new int[nb_row - 1];
    for (i = 1;i < nb_row;i++) {
      nb_memory[i] = 0;

//      if (next[i]) {
      if ((type == ORDINARY) || (!child[i])) {
        for (j = 1;j < nb_row;j++) {
          if ((next[j]) && (next[j][state[i][0]] == i)) {
            buffer[nb_memory[i]] = j;
            nb_memory[i]++;
          }
        }
      }

      if (nb_memory[i] > 0) {
        previous[i] = new int[nb_memory[i]];
        for (j = 0;j < nb_memory[i];j++) {
          previous[i][j] = buffer[j];
        }
      }
      else {
        previous[i] = NULL;
      }
    }

    delete [] buffer;
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Checking that the set of terminal memories is suffix-free.
 *
 *  \return suffix-free or not.
 */
/*--------------------------------------------------------------*/

bool VariableOrderMarkovChain::check_free_suffix() const

{
  bool free_suffix = true;
  int i , j , k;


  for (i = 1;i < nb_row;i++) {
    if ((!child[i]) && (order[i] >= 2)) {
      for (j = 1;j < nb_row;j++) {
        if ((!child[j]) && (order[j] < order[i])) {
          for (k = 0;k < order[j];k++) {
            if (state[i][k] != state[j][k]) {
              break;
            }
          }

          if (k == order[j]) {
            free_suffix = false;
            break;
          }
        }
      }

      if (j < nb_row) {
        break;
      }
    }
  }

  return free_suffix;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Construction of the matrix of possible transitions between states
 *         (adjacency matrix of the graph of possible transitions).
 *
 *  \return matrix of possible transitions.
 */
/*--------------------------------------------------------------*/

bool** VariableOrderMarkovChain::logic_transition_computation() const

{
  bool **logic_transition;
  int i , j;
  double **order1_transition;


  logic_transition = new bool*[nb_state];
  for (i = 0;i < nb_state;i++) {
    logic_transition[i] = new bool[nb_state];
    logic_transition[i][i] = false;
  }

  order1_transition = new double*[nb_state];
  for (i = 0;i < nb_state;i++) {
    order1_transition[i] = new double[nb_state];
    for (j = 0;j < nb_state;j++) {
      order1_transition[i][j] = 0.;
    }
  }

  for (i = 1;i < nb_row;i++) {
    if (memo_type[i] == TERMINAL) {
      for (j = 0;j < nb_state;j++) {
        order1_transition[state[i][0]][j] += transition[i][j];
      }
    }
  }

  for (i = 0;i < nb_state;i++) {
    for (j = 0;j < nb_state;j++) {
      if (j != i) {
        logic_transition[i][j] = (order1_transition[i][j] == 0. ? false : true);
      }
    }
  }

  for (i = 0;i < nb_state;i++) {
    delete [] order1_transition[i];
  }
  delete [] order1_transition;

  return logic_transition;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Extraction of the variable-order Markov chain classes
 *         (transient/recurrent/absorbing) from the state accessibility.
 */
/*--------------------------------------------------------------*/

void VariableOrderMarkovChain::component_computation()

{
  bool **logic_transition;
  int i;


  logic_transition = logic_transition_computation();

  Chain::component_computation(logic_transition);

  for (i = 0;i < nb_state;i++) {
    delete [] logic_transition[i];
  }
  delete [] logic_transition;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Construction of the non-terminal memories.
 */
/*--------------------------------------------------------------*/

void VariableOrderMarkovChain::build_non_terminal()

{
  int i , j , k , m;
  int nb_non_terminal , nb_terminal , bnb_memory;


  nb_non_terminal = (nb_row - 1) / nb_state;
  nb_terminal = (nb_row - 1) * (nb_state - 1) / nb_state + 1;
//  nb_terminal = nb_non_terminal * (nb_state - 1) + 1;

  parent[0] = -1;

  i = 0;
  for (j = nb_non_terminal;j < nb_non_terminal + nb_terminal;j++) {
    for (k = order[j] - 1;k >= 0;k--) {
      if (state[j][k] != 0) {
        break;
      }
    }
    bnb_memory = order[j] - 1 - k;

    // insertion of the non-terminal memories

    for (k = order[j] - bnb_memory;k < order[j];k++) {
      memo_type[i] = NON_TERMINAL;

      for (m = 0;m < nb_state;m++) {
        transition[i][m] = 0.;
      }

      order[i] = k;
      for (m = 0;m < order[i];m++) {
        state[i][m] = state[j][m];
      }

      child[i] = new int[nb_state];
      i++;
    }

    // copy of the terminal memory

    memo_type[i] = TERMINAL;
    if (i < j) {
      for (k = 0;k < nb_state;k++) {
        transition[i][k] = transition[j][k];
      }
      order[i] = order[j];
      for (k = 0;k < order[i];k++) {
        state[i][k] = state[j][k];
      }
    }

    // parent-children relationships

    for (k = 0;k < bnb_memory;k++) {
      parent[i - k] = i - k - 1;
      child[i - k - 1][state[i - k][order[i - k] - 1]] = i - k;
    }

    find_parent_memory(i - bnb_memory);
    i++;
  }

# ifdef DEBUG
  {
    cout << "\n";
    cout << "Suffix free? " << (check_free_suffix() ? "True" : "False") << endl;

    for (i = 0;i < nb_row;i++) {
      for (j = max_order - 1;j >= order[i];j--) {
        cout << "  ";
      }
      for (j = order[i] - 1;j >= 0;j--) {
        cout << state[i][j] << " ";
      }

      switch (memo_type[i]) {
      case NON_TERMINAL :
        cout << "  " << SEQ_label[SEQL_NON_TERMINAL];
        break;
      case TERMINAL :
        cout << "      " << SEQ_label[SEQL_TERMINAL];
        break;
      }

      cout << " | " << parent[i];

      if (child[i]) {
        cout << " |";
        for (j = 0;j < nb_state;j++) {
          cout << " " << child[i][j];
        }
      }

      else {
        cout << "  ";
        for (j = 0;j < nb_state;j++) {
          cout << "  ";
        }
      }

      cout << endl;
    }
  }
# endif

}


/*--------------------------------------------------------------*/
/**
 *  \brief Application of a threshold on the variable-order Markov chain parameters.
 *
 *  \param[in] min_probability minimum probability.
 */
/*--------------------------------------------------------------*/

void VariableOrderMarkovChain::thresholding(double min_probability)

{
  bool stop;
  int i , j;
  int nb_correction;
  double norm;


  if (min_probability > THRESHOLDING_FACTOR / (double)nb_state) {
    min_probability = THRESHOLDING_FACTOR / (double)nb_state;
  }

  if (type == ORDINARY) {
    do {
      stop = true;
      nb_correction = 0;
      norm = 0.;

      for (i = 0;i < nb_state;i++) {
        if (initial[i] <= min_probability) {
          nb_correction++;
          initial[i] = min_probability;
        }
        else {
          norm += initial[i];
        }
      }

      if (nb_correction > 0) {
        for (i = 0;i < nb_state;i++) {
          if (initial[i] > min_probability) {
            initial[i] *= (1. - nb_correction * min_probability) / norm;
            if (initial[i] < min_probability) {
              stop = false;
            }
          }
        }
      }
    }
    while (!stop);
  }

  for (i = 1;i < nb_row;i++) {
    if ((memo_type[i] == TERMINAL) || ((type == ORDINARY) &&
         (memo_type[i] == NON_TERMINAL))) {
      do {
        stop = true;
        nb_correction = 0;
        norm = 0.;

        for (j = 0;j < nb_state;j++) {
          if ((accessibility[state[i][0]][j]) && (transition[i][j] <= min_probability)) {
            nb_correction++;
            transition[i][j] = min_probability;
          }
          else {
            norm += transition[i][j];
          }
        }

        if (nb_correction > 0) {
          for (j = 0;j < nb_state;j++) {
            if (transition[i][j] > min_probability) {
              transition[i][j] *= (1. - nb_correction * min_probability) / norm;
              if (transition[i][j] < min_probability) {
                stop = false;
              }
            }
          }
        }
      }
      while (!stop);
    }

    else if (memo_type[i] == COMPLETION) {
      for (j = 0;j < nb_state;j++) {
        transition[i][j] = transition[parent[i]][j];
      }
    }
  }

/*  if (accessibility) {
    for (i = 0;i < nb_state;i++) {
      delete [] accessibility[i];
    }
    delete [] accessibility;
  }
  accessibility = NULL;

  delete [] component_nb_state;

  if (component) {
    for (i = 0;i < nb_component;i++) {
      delete [] component[i];
    }
    delete [] component;
  }
  nb_component = 0;

  component_computation(); */

  if (next) {
    for (i = 1;i < nb_row;i++) {
      delete [] next[i];
    }
    delete [] next;
  }
  next = NULL;

  delete [] nb_memory;
  nb_memory = NULL;

  if (previous) {
    for (i = 1;i < nb_row;i++) {
      delete [] previous[i];
    }
    delete [] previous;
  }
  previous = NULL;

  build_memory_transition();
  build_previous_memory();
}


/*--------------------------------------------------------------*/
/**
 *  \brief Determination of the maximum memory order.
 */
/*--------------------------------------------------------------*/

void VariableOrderMarkovChain::max_order_computation()

{
  int i;


  max_order = 0;
  for (i = 0;i < nb_row;i++) {
    if ((memo_type[i] == TERMINAL) && (order[i] > max_order)) {
      max_order = order[i];
    }
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the number of parameters of a VariableOrderMarkovChain object.
 *
 *  \param[in] min_probability minimum probability.
 *
 *  \return                    number of parameters.
 */
/*--------------------------------------------------------------*/

int VariableOrderMarkovChain::nb_parameter_computation(double min_probability) const

{
  int i , j;
  int nb_parameter = 0;


  // particular case order 0

  if (max_order == 1) {
    for (i = 0;i < nb_state - 1;i++) {
      for (j = 2;j <= nb_state;j++) {
        if (transition[j][i] != transition[1][i]) {
          break;
        }
      }

      if (j <= nb_state) {
        break;
      }
    }

    if (i == nb_state - 1) {
      for (i = 0;i < nb_state;i++) {
        if (transition[0][i] > min_probability) {
          nb_parameter++;
        }
      }

      nb_parameter--;
    }
  }

  if (nb_parameter == 0) {
    for (i = 1;i < nb_row;i++) {
      if (memo_type[i] == TERMINAL) {
//      if ((memo_type[i] == TERMINAL) || ((type == ORDINARY) &&
//           (memo_type[i] == NON_TERMINAL))) {
        for (j = 0;j < nb_state;j++) {
          if (transition[i][j] > min_probability) {
            nb_parameter++;
          }
        }

        nb_parameter--;
      }
    }
  }

  return nb_parameter;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the number of transient parameters corresponding to
 *         the transition distributions attached to the non-terminal memories of
 *         an ordinary variable-order Markov chain.
 *
 *  \param[in] min_probability minimum probability.
 *
 *  \return                    number of transient parameters.
 */
/*--------------------------------------------------------------*/

int VariableOrderMarkovChain::nb_transient_parameter_computation(double min_probability) const

{
  int i , j;
  int nb_parameter = 0;


  if (type == ORDINARY) {
    for (i = 1;i < nb_row;i++) {
      if (memo_type[i] == NON_TERMINAL) {
        for (j = 0;j < nb_state;j++) {
          if (transition[i][j] > min_probability) {
            nb_parameter++;
          }
        }

        nb_parameter--;
      }
    }
  }

  return nb_parameter;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Analysis of the format of a VariableOrderMarkovChain object.
 *
 *  \param[in] error   reference on a StatError object,
 *  \param[in] in_file stream,
 *  \param[in] line    reference on the file line index,
 *  \param[in] type    process type (ORDINARY/EQUILIBRIUM).
 *
 *  \return            VariableOrderMarkovChain object.
 */
/*--------------------------------------------------------------*/

VariableOrderMarkovChain* VariableOrderMarkovChain::parsing(StatError &error , ifstream &in_file ,
                                                            int &line , process_type type)

{
  string buffer;
  size_t position;
  typedef tokenizer<char_separator<char>> tokenizer;
  char_separator<char> separator(" \t");
  streampos transition_line;
  bool status = true , lstatus , increase , **logic_transition;
  int i , j;
  int read_line , tline , value , nb_state = 0 , order , previous_order , max_order = 0 , buff ,
      nb_terminal , nb_non_terminal , memory , state[ORDER] , previous_state[ORDER];
  double proba , cumul , *initial;
  VariableOrderMarkovChain *markov;


  markov = NULL;

  // analysis of the line defining the number of states

  while (getline(in_file , buffer)) {
    line++;

#   ifdef DEBUG
    cout << line << "  " << buffer << endl;
#   endif

    position = buffer.find('#');
    if (position != string::npos) {
      buffer.erase(position);
    }
    i = 0;

    tokenizer tok_buffer(buffer , separator);

    for (tokenizer::iterator token = tok_buffer.begin();token != tok_buffer.end();token++) {
      switch (i) {

      // test number of states

      case 0 : {
        lstatus = true;

/*        try {
          value = stoi(*token);   in C++ 11
        }
        catch(invalid_argument &arg) {
          lstatus = false;
        } */
        value = atoi(token->c_str());

        if (lstatus) {
          if ((value < 2) || (value > NB_STATE)) {
            lstatus = false;
          }
          else {
            nb_state = value;
          }
        }

        if (!lstatus) {
          status = false;
          error.update(STAT_parsing[STATP_NB_STATE] , line , i + 1);
        }
        break;
      }

      // test STATES keyword

      case 1 : {
        if (*token != STAT_word[STATW_STATES]) {
          status = false;
          error.correction_update(STAT_parsing[STATP_KEYWORD] , STAT_word[STATW_STATES] , line , i + 1);
        }
        break;
      }
      }

      i++;
    }

    if (i > 0) {
      if (i != 2) {
        status = false;
        error.update(STAT_parsing[STATP_FORMAT] , line);
      }
      break;
    }
  }

  if (nb_state == 0) {
    status = false;
    error.update(STAT_parsing[STATP_FORMAT] , line);
  }

  if (status) {
    initial = new double[nb_state];

    // 1st pass: search for the number of terminal memories and the maximum order,
    // analysis of the initial probabilities, transition probabilities and memories

    read_line = 0;
    while (getline(in_file , buffer)) {
      line++;

#     ifdef DEBUG
      cout << line << "  " << buffer << endl;
#     endif

      position = buffer.find('#');
      if (position != string::npos) {
        buffer.erase(position);
      }
      i = 0;

      tokenizer tok_buffer(buffer , separator);

      if ((read_line == 0) || ((type == ORDINARY) && (read_line == 2))) {
        for (tokenizer::iterator token = tok_buffer.begin();token != tok_buffer.end();token++) {

          // test INITIAL_PROBABILITIES/TRANSITION_PROBABILITIES keyword

          if (i == 0) {
            if ((type == ORDINARY) && (read_line == 0)) {
              if (*token != STAT_word[STATW_INITIAL_PROBABILITIES]) {
                status = false;
                error.correction_update(STAT_parsing[STATP_KEYWORD] , STAT_word[STATW_INITIAL_PROBABILITIES] , line);
              }
            }

            else {
              if (*token != STAT_word[STATW_TRANSITION_PROBABILITIES]) {
                status = false;
                error.correction_update(STAT_parsing[STATP_KEYWORD] , STAT_word[STATW_TRANSITION_PROBABILITIES] , line);
              }
            }
          }

          i++;
        }

        if (i > 0) {
          if (((type == ORDINARY) && (read_line == 2)) || ((type == EQUILIBRIUM) && (read_line == 0))) {
            transition_line = in_file.tellg();
            tline = line;
          }

          if (i != 1) {
            status = false;
            error.update(STAT_parsing[STATP_FORMAT] , line);
          }
        }
      }

      else {
        cumul = 0.;

        for (tokenizer::iterator token = tok_buffer.begin();token != tok_buffer.end();token++) {
          if (i < nb_state) {
            lstatus = true;

/*            try {
              proba = stod(*token);   in C++ 11
            }
            catch (invalid_argument &arg) {
              lstatus = false;
            } */
            proba = atof(token->c_str());

            if (lstatus) {
              if ((proba < 0.) || (proba > 1. - cumul + DOUBLE_ERROR)) {
                lstatus = false;
              }

              else {
                cumul += proba;
                if ((type == ORDINARY) && (read_line == 1)) {
                  initial[i] = proba;
                }
              }
            }

            if (!lstatus) {
              status = false;
              if ((type == ORDINARY) && (read_line == 1)) {
                error.update(STAT_parsing[STATP_INITIAL_PROBABILITY] , line , i + 1);
              }
              else {
                error.update(STAT_parsing[STATP_TRANSITION_PROBABILITY] , line , i + 1);
              }
            }
          }

          else if ((type == EQUILIBRIUM) || (read_line >= 3)) {
            lstatus = true;

/*            try {
              value = stoi(*token);   in C++ 11
            }
            catch(invalid_argument &arg) {
              lstatus = false;
            } */
            value = atoi(token->c_str());

            if (lstatus) {
              if ((value < 0) || (value >= nb_state)) {
                lstatus = false;
              }
              else if (i - nb_state < ORDER) {
                state[i - nb_state] = value;
              }
            }

            if (!lstatus) {
              status = false;
              error.update(SEQ_parsing[SEQP_STATE] , line , i + 1);
            }
          }

          i++;
        }

        if (i > 0) {
          if ((type == ORDINARY) && (read_line == 1) && (i != nb_state)) {
            status = false;
            error.update(STAT_parsing[STATP_FORMAT] , line);
          }

          if (cumul < 1. - DOUBLE_ERROR) {
            status = false;
            error.update(STAT_parsing[STATP_PROBABILITY_SUM] , line);
          }

          if (((type == ORDINARY) && (read_line >= 3)) || ((type == EQUILIBRIUM) && (read_line >= 1))) {
            if (i <= nb_state) {
              status = false;
              error.update(STAT_parsing[STATP_FORMAT] , line);
            }

            else {
              order = i - nb_state;

              if (order > ORDER) {
                status = false;
                error.update(STAT_parsing[STATP_ORDER] , line);
              }

              else if (order > max_order) {
                max_order = order;
              }

              if (status) {
                for (j = 0;j < order / 2;j++) {
                  buff = state[j];
                  state[j] = state[order - j - 1];
                  state[order - j - 1] = buff;
                }

                if (read_line - (type == ORDINARY ? 3 : 1) == 0) {
                  for (j = 0;j < order;j++) {
                    if (state[j] != 0) {
                      status = false;
                      error.update(SEQ_parsing[SEQP_STATE] , line , nb_state + order - j);
                    }
                  }
                }

                else {

                  // checking of the memory succession

                  increase = true;
                  for (j = MIN(previous_order , order) - 1;j >= 0;j--) {
                    if (increase) {
                      if (previous_state[j] < nb_state - 1) {
                        if (state[j] != previous_state[j] + 1) {
                          status = false;
                          error.update(SEQ_parsing[SEQP_STATE] , line , nb_state + order - j);
                        }
                        increase = false;
                      }

                      else {
                        if (state[j] != 0) {
                          status = false;
                          error.update(SEQ_parsing[SEQP_STATE] , line , nb_state + order - j);
                        }
                      }
                    }

                    else {
                      if (state[j] != previous_state[j]) {
                        status = false;
                        error.update(SEQ_parsing[SEQP_STATE] , line , nb_state + order - j);
                      }
                    }
                  }

                  // case increase of the memory length or stable memory length

                  if (order >= previous_order) {
                    for (j = order - 1;j >= previous_order;j--) {
                      if (state[j] != 0) {
                        status = false;
                        error.update(SEQ_parsing[SEQP_STATE] , line , nb_state + order - j);
                      }
                    }
                  }

                  // case decrease of the memory length

                  else {
                    for (j = order;j < previous_order;j++) {
                      if (previous_state[j] != nb_state - 1) {
                        status = false;
                        error.update(SEQ_parsing[SEQP_STATE] , line);
                      }
                    }
                  }
                }

                // search for the last memory

                for (j = 0;j < order;j++) {
                  if (state[j] != nb_state - 1) {
                    break;
                  }
                }

                if (j == order) {
                  read_line++;
                  break;
                }

                else {
                  for (j = 0;j < order;j++) {
                    previous_state[j] = state[j];
                  }
                  previous_order = order;
                }
              }
            }
          }
        }
      }

      if (i > 0) {
        read_line++;
      }
    }

    // checking of the number of memories

    nb_terminal = read_line - (type == ORDINARY ? 3 : 1);
    if ((nb_state > 2) && (nb_terminal % (nb_state - 1) != 1)) {
      status = false;
      error.update(SEQ_parsing[SEQP_NB_MEMORY] , line);
    }

    if (status) {
      in_file.clear();
      in_file.seekg(transition_line);

      nb_non_terminal = (nb_terminal - 1) / (nb_state - 1);
      markov = new VariableOrderMarkovChain(type , nb_state , nb_non_terminal + nb_terminal , max_order);

      if (type == ORDINARY) {
        for (i = 0;i < nb_state;i++) {
          markov->initial[i] = initial[i];
        }
      }

      // 2nd pass: reading of the transition probabilities and the memories

      memory = nb_non_terminal;
      line = tline;
      while (getline(in_file , buffer)) {
        line++;

#       ifdef DEBUG
        cout << line << "  " << buffer << endl;
#       endif

        position = buffer.find('#');
        if (position != string::npos) {
          buffer.erase(position);
        }
        i = 0;

        tokenizer tok_buffer(buffer , separator);

        for (tokenizer::iterator token = tok_buffer.begin();token != tok_buffer.end();token++) {
          if (i < nb_state) {
//            markov->transition[memory][i] = stod(*token);   in C++ 11
            markov->transition[memory][i] = atof(token->c_str());
          }

          else {
//            state[i - nb_state] = stoi(*token);   in C++ 11
            state[i - nb_state] = atoi(token->c_str());
          }

          i++;
        }

        if (i > 0) {
          markov->order[memory] = i - nb_state;
          for (j = 0;j < markov->order[memory];j++) {
            markov->state[memory][j] = state[markov->order[memory] - j - 1];
          }

          // search for the last memory

          for (j = 0;j < markov->order[memory];j++) {
            if (markov->state[memory][j] != nb_state - 1) {
              break;
            }
          }

          if (j == order) {
            break;
          }
          else {
            memory++;
          }
        }
      }

      markov->build_non_terminal();

      // test accessible states

      logic_transition = markov->logic_transition_computation();
      status = markov->strongly_connected_component_research(error , logic_transition);

      if (status) {

        // test irreducibility in the equilibrium process case

        markov->Chain::component_computation(logic_transition);
        if ((type == EQUILIBRIUM) && (markov->nb_component > 1)) {
          status = false;
          error.correction_update(STAT_parsing[STATP_CHAIN_STRUCTURE] , STAT_parsing[STATP_IRREDUCIBLE]);
        }
      }

      for (i = 0;i < nb_state;i++) {
        delete [] logic_transition[i];
      }
      delete [] logic_transition;

      if (!status) {
        delete markov;
        markov = NULL;
      }
      else {
        markov->max_order_computation();
      }
    }

    delete [] initial;
  }

  return markov;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of the memory out-tree.
 *
 *  \param[in,out] os        stream,
 *  \param[in]     file_flag flag file.
 */
/*--------------------------------------------------------------*/

ostream& VariableOrderMarkovChain::ascii_memory_tree_print(ostream &os , bool file_flag) const

{
  int i , j , k;
  int bnb_memory , width = column_width(nb_state);
  ios_base::fmtflags format_flags;


  format_flags = os.setf(ios::left , ios::adjustfield);

  os << "\n";
  if (file_flag) {
    os << "# ";
  }
  os << SEQ_label[SEQL_MEMORY_TREE] << endl;

  os << "\n";
  for (i = 0;i < nb_row;i++) {
    if (!child[i]) {
      if (file_flag) {
        os << "# ";
      }

      for (j = order[i] - 1;j >= 1;j--) {
        if (state[i][j] != 0) {
          break;
        }
      }
      bnb_memory = order[i] - j;

      for (j = 0;j < order[i] - bnb_memory;j++) {
        if (state[i][j] < nb_state - 1) {
          os << "|";
        }
        else {
          os << " ";
        }

        for (k = 0;k < (j + 1) * (width + 1) + 1;k++) {
          os << " ";
        }
      }

      os << "|";
      for (j = order[i] - bnb_memory;j < order[i];j++) {
        for (k = 0;k < 3;k++) {
          os << "_";
        }
        for (k = j;k >= 0;k--) {
          os << setw(width) << state[i][k];
          if (k > 0) {
            os << " ";
          }
        }
      }

      os << endl;
    }
  }

  os.setf(format_flags , ios::adjustfield);

# ifdef MESSAGE
  os << "\n";
  for (i = 0;i < nb_row;i++) {
    if (child[i]) {
      if (file_flag) {
        os << "# ";
      }

      for (j = max_order - 1;j >= order[i];j--) {
        os << "  ";
      }
      for (j = order[i] - 1;j >= 0;j--) {
        os << state[i][j] << " ";
      }

      for (j = 0;j < nb_state;j++) {
        os << "  ";
        for (k = max_order - 1;k >= order[child[i][j]];k--) {
          os << "  ";
        }
        for (k = order[child[i][j]] - 1;k >= 0;k--) {
          os << state[child[i][j]][k] << " ";
        }
      }

      if (memo_type[i] == COMPLETION) {
        os << "    " << SEQ_label[SEQL_COMPLETION];
      }

      os << endl;
    }
  }
# endif

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of the memory prefix in-tree.
 *
 *  \param[in,out] os        stream,
 *  \param[in]     file_flag flag file.
 */
/*--------------------------------------------------------------*/

ostream& VariableOrderMarkovChain::ascii_transition_tree_print(ostream &os , bool file_flag) const

{
  int i , j , k;
  int min_order , nb_root , memory , width = column_width(nb_state) , *nb_next_memory ,
      *root , *nb_leaf_memory , *nb_drawn_next_memory , *nb_drawn_leaf_memory;
  ios_base::fmtflags format_flags;


  // computation of the number of following memories of higher length 

  nb_next_memory = new int[nb_row];
  root = new int[nb_row];
  nb_leaf_memory = new int[nb_row];
  nb_drawn_next_memory = new int[nb_row];
  nb_drawn_leaf_memory = new int[nb_row];

  min_order = max_order;

  for (i = 1;i < nb_row;i++) {
    nb_next_memory[i] = 0;

    if ((type == ORDINARY) || (!child[i])) {
      for (j = 0;j < nb_state;j++) {
        if (order[next[i][j]] == order[i] + 1) {
          nb_next_memory[i]++;
        }
      }

      if (nb_next_memory[i] == 0) {
        nb_leaf_memory[i] = 1;
      }
      else {
        nb_leaf_memory[i] = 0;
      }

      if (order[i] < min_order) {
        min_order = order[i];
      }
    }
  }

  // computation of the roots

  nb_root = 0;
  for (i = 1;i < nb_row;i++) {
    if ((type == ORDINARY) || (!child[i])) {
      for (j = 0;j < nb_memory[i];j++) {
        if (order[previous[i][j]] == order[i] - 1) {
          break;
        }
      }

      if (j == nb_memory[i]) {
        root[nb_root++] = i;
      }
    }
  }

  // computation of the number of leaf memories

  for (i = max_order;i >= min_order + 1;i--) {
    for (j = 1;j < nb_row;j++) {
      if ((order[j] == i) && (nb_leaf_memory[j] > 0) &&
          ((type == ORDINARY) || (!child[j]))) {
        for (k = 0;k < nb_memory[j];k++) {
          if ((order[previous[j][k]] == order[j] - 1) &&
              ((type == ORDINARY) || (!child[previous[j][k]]))) {
            nb_leaf_memory[previous[j][k]] += nb_leaf_memory[j];
            break;
          }
        }
      }
    }
  }

# ifdef DEBUG
  cout << "\n";
  for (i = 1;i < nb_row;i++) {
    if ((type == ORDINARY) || (!child[i])) {
      for (j = max_order - 1;j >= order[i];j--) {
        cout << "  ";
      }
      for (j = order[i] - 1;j >= 0;j--) {
        cout << state[i][j] << " ";
      }

      cout << "   "  << nb_next_memory[i] << " | " << nb_leaf_memory[i];
      for (j = 0;j < nb_root;j++) {
        if (root[j] == i) {
          cout << "   root";
          break;
        }
      }
      cout << endl;
    }
  }
# endif

  for (i = 1;i < nb_row;i++) {
    if ((type == ORDINARY) || (!child[i])) {
      nb_drawn_next_memory[i] = nb_next_memory[i];
      nb_drawn_leaf_memory[i] = nb_leaf_memory[i];
    }
  }

  format_flags = os.setf(ios::left , ios::adjustfield);

  os << "\n";
  if (file_flag) {
    os << "# ";
  }
  os << SEQ_label[SEQL_TRANSITION_TREE] << endl;

  for (i = 0;i < nb_root;i++) {
    os << "\n";
    for (j = 0;j < nb_leaf_memory[root[i]];j++) {
      if (file_flag) {
        os << "# ";
      }

      memory = root[i];
      for (;;) {
        if (!child[memory]) {
          if (nb_drawn_leaf_memory[memory] == nb_leaf_memory[memory]) {
            if (memory != root[i]) {
              for (k = 0;k < 3;k++) {
                os << "_";
              }
            }

            for (k = order[memory] - 1;k >= 0;k--) {
              os << setw(width) << state[memory][k];
              if (k > 0) {
                os << " ";
              }
            }
          }

          else {
            if (memory != root[i]) {
              for (k = 0;k < 3;k++) {
                os << " ";
              }
            }

            for (k = 0;k < order[memory] * (width + 1) - 2;k++) {
              os << " ";
            }

            if (nb_drawn_next_memory[memory] > 0) {
              os << "|";
            }
            else {
              os << " ";
            }
          }
        }

        else {
          if (memory != root[i]) {
            for (k = 0;k < 3;k++) {
              os << " ";
            }
          }

          for (k = 0;k < order[memory] * (width + 1) - 1;k++) {
            os << " ";
          }
        }

        nb_drawn_leaf_memory[memory]--;

        if (nb_next_memory[memory] == 0) {
          os << endl;
          break;
        }

        for (k = 0;k < nb_state;k++) {
          if ((order[next[memory][k]] == order[memory] + 1) &&
              (nb_drawn_leaf_memory[next[memory][k]] > 0)) {
            if (nb_drawn_leaf_memory[next[memory][k]] == nb_leaf_memory[next[memory][k]]) {
              nb_drawn_next_memory[memory]--;
            }
            memory = next[memory][k];
            break;
          }
        }
      }
    }
  }

  os.setf(format_flags , ios::adjustfield);

  delete [] nb_next_memory;
  delete [] root;
  delete [] nb_leaf_memory;
  delete [] nb_drawn_next_memory;
  delete [] nb_drawn_leaf_memory;

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of the variable-order Markov chain parameters.
 *
 *  \param[in,out] os        stream,
 *  \param[in]     file_flag flag file.
 */
/*--------------------------------------------------------------*/

ostream& VariableOrderMarkovChain::ascii_print(ostream &os , bool file_flag) const

{
  int i , j , k;
  int buff , width;
  double *stationary_probability;
  ios_base::fmtflags format_flags;


  format_flags = os.setf(ios::left , ios::adjustfield);

  os << "\n" << nb_state << " " << STAT_word[STATW_STATES] << endl;

  // computation of the column width

  width = column_width((type == ORDINARY ? nb_state : nb_row) , initial);

  for (i = 1;i < nb_row;i++) {
    buff = column_width(nb_state , transition[i]);
    if (buff > width) {
      width = buff;
    }
  }
  width += ASCII_SPACE;

  os << "\n";

  switch (type) {

  case ORDINARY : {
    os << STAT_word[STATW_INITIAL_PROBABILITIES] << endl;
    for (i = 0;i < nb_state;i++) {
      os << setw(width) << initial[i];
    }
    os << endl;
    break;
  }

  case EQUILIBRIUM : {

    // computation of the stationary probabilities corresponding to the memories added by completion

    stationary_probability = new double[nb_row];

    for (i = 1;i < nb_row;i++) {
      stationary_probability[i] = initial[i];
    }
    for (i = nb_row - 1;i >= 1;i--) {
      if (memo_type[i] == COMPLETION) {
        stationary_probability[parent[i]] += stationary_probability[i];
      }
    }

    if (file_flag) {
      os << "# ";
    }
    os << STAT_label[STATL_STATIONARY_PROBABILITIES] << endl;

    for (i = 1;i < nb_row;i++) {
      if (memo_type[i] == TERMINAL) {
        if (file_flag) {
          os << "# ";
        }
        os << setw(width) << stationary_probability[i];

        os << " ";
        for (j = max_order - 1;j >= order[i];j--) {
          os << "  ";
        }
        for (j = order[i] - 1;j >= 0;j--) {
          os << state[i][j] << " ";
        }
        os << endl;
      }
    }

    delete [] stationary_probability;
    break;
  }
  }

  os << "\n" << STAT_word[STATW_TRANSITION_PROBABILITIES] << "      ";
  if (file_flag) {
    os << "# ";
  }
  os << STAT_label[STATL_MEMORY] << endl;

  for (i = 1;i < nb_row;i++) {
    if ((memo_type[i] != TERMINAL) && (file_flag)) {
      os << "# ";
    }

    for (j = 0;j < nb_state;j++) {
      os << setw(width) << transition[i][j];
    }

    os << " ";
    for (j = max_order - 1;j >= order[i];j--) {
      os << "  ";
    }
    for (j = order[i] - 1;j >= 0;j--) {
      os << state[i][j] << " ";
    }

    if ((memo_type[i] == TERMINAL) && (file_flag)) {
      os << "# ";
    }

    switch (memo_type[i]) {

    case NON_TERMINAL : {
      os << "  " << SEQ_label[SEQL_NON_TERMINAL];
      break;
    }

    case TERMINAL : {
      os << "      " << SEQ_label[SEQL_TERMINAL];
      if (child[i]) {
        os << " (" << SEQ_label[SEQL_COMPLETED] << ")";
      }
      break;
    }

    case COMPLETION : {
      os << "    " << SEQ_label[SEQL_COMPLETION] << " ("
         << (child[i] ? SEQ_label[SEQL_NON_TERMINAL] : SEQ_label[SEQL_TERMINAL]) << ")";
      break;
    }
    }

    os << endl;
  }

  if (nb_component > 0) {
    for (i = 0;i < nb_component;i++) {
      os << "\n";
      if (file_flag) {
        os << "# ";
      }

      switch (stype[component[i][0]]) {
      case TRANSIENT :
        os << STAT_label[STATL_TRANSIENT] << " ";
        break;
      default :
        os << STAT_label[STATL_RECURRENT] << " ";
        break;
      }
      os << STAT_label[STATL_CLASS] << ": " << STAT_label[component_nb_state[i] == 1 ? STATL_STATE : STATL_STATES];

      for (j = 0;j < component_nb_state[i];j++) {
        os << " " << component[i][j];
      }

      if (stype[component[i][0]] == ABSORBING) {
        os << " (" << STAT_label[STATL_ABSORBING] << " " << STAT_label[STATL_STATE] << ")";
      }
    }
    os << endl;
  }

  ascii_memory_tree_print(os , file_flag);
  ascii_transition_tree_print(os , file_flag);

  os << "\n";
  if (file_flag) {
    os << "# ";
  }
  os << SEQ_label[SEQL_MEMORY_TRANSITION_MATRIX] << endl;

  os << "\n";
  for (i = 1;i < nb_row;i++) {
    if ((type == ORDINARY) || (!child[i])) {
      if (file_flag) {
        os << "# ";
      }

      for (j = max_order - 1;j >= order[i];j--) {
        os << "  ";
      }
      for (j = order[i] - 1;j >= 0;j--) {
        os << state[i][j] << " ";
      }

      for (j = 0;j < nb_state;j++) {
        os << "  ";
        for (k = max_order - 1;k >= order[next[i][j]];k--) {
          os << "  ";
        }
        for (k = order[next[i][j]] - 1;k >= 0;k--) {
          os << state[next[i][j]][k] << " ";
        }
      }

      switch (memo_type[i]) {

      case NON_TERMINAL : {
        os << "  " << SEQ_label[SEQL_NON_TERMINAL];
        break;
      }

      case TERMINAL : {
        os << "      " << SEQ_label[SEQL_TERMINAL];
        if (child[i]) {
          os << " (" << SEQ_label[SEQL_COMPLETED] << ")";
        }
        break;
      }

      case COMPLETION : {
        os << "    " << SEQ_label[SEQL_COMPLETION] << " ("
           << (child[i] ? SEQ_label[SEQL_NON_TERMINAL] : SEQ_label[SEQL_TERMINAL]) << ")";
        break;
      }
      }

      os << endl;
    }
  }

# ifdef DEBUG
  if (previous) {
    os << "\n";
    for (i = 1;i < nb_row;i++) {
      if ((type == ORDINARY) || (!child[i])) {
        for (j = max_order - 1;j >= order[i];j--) {
          os << "  ";
        }
        for (j = order[i] - 1;j >= 0;j--) {
          os << state[i][j] << " ";
        }

        for (j = 0;j < nb_memory[i];j++) {
          os << "  ";
          for (k = max_order - 1;k >= order[previous[i][j]];k--) {
            os << "  ";
          }
          for (k = order[previous[i][j]] - 1;k >= 0;k--) {
            os << state[previous[i][j]][k] << " ";
          }
        }
        os << endl;
      }
    }
  }
# endif

  os.setf(format_flags , ios::adjustfield);

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of the variable-order Markov chain parameters at the spreadsheet format.
 *
 *  \param[in,out] os stream.
 */
/*--------------------------------------------------------------*/

ostream& VariableOrderMarkovChain::spreadsheet_print(ostream &os) const

{
  int i , j , k;
  double *stationary_probability;


  os << "\n" << nb_state << "\t" << STAT_word[STATW_STATES] << endl;

  switch (type) {

  case ORDINARY : {
    os << "\n" << STAT_word[STATW_INITIAL_PROBABILITIES] << endl;
    for (i = 0;i < nb_state;i++) {
      os << initial[i] << "\t";
    }
    os << endl;
    break;
  }

  case EQUILIBRIUM : {

    // computation of the stationary probabilities corresponding to the memories added by completion

    stationary_probability = new double[nb_row];

    for (i = 1;i < nb_row;i++) {
      stationary_probability[i] = initial[i];
    }
    for (i = nb_row - 1;i >= 1;i--) {
      if (memo_type[i] == COMPLETION) {
        stationary_probability[parent[i]] += stationary_probability[i];
      }
    }

    os << "\n" << STAT_label[STATL_STATIONARY_PROBABILITIES] << endl;
    for (i = 1;i < nb_row;i++) {
      if (memo_type[i] == TERMINAL) {
        os << stationary_probability[i] << "\t\t";
        for (j = order[i] - 1;j >= 0;j--) {
          os << state[i][j] << " ";
        }
        os << endl;
      }
    }

    delete [] stationary_probability;
    break;
  }
  }

  os << "\n" << STAT_word[STATW_TRANSITION_PROBABILITIES]
     << "\t\t" << STAT_label[STATL_MEMORY] << endl;

  for (i = 1;i < nb_row;i++) {
    for (j = 0;j < nb_state;j++) {
      os << transition[i][j] << "\t";
    }

    os << "\t";
    for (j = order[i] - 1;j >= 0;j--) {
      os << state[i][j] << " ";
    }

    os << "\t";
    switch (memo_type[i]) {

    case NON_TERMINAL : {
      os << SEQ_label[SEQL_NON_TERMINAL];
      break;
    }

    case TERMINAL : {
      os << SEQ_label[SEQL_TERMINAL];
      if (child[i]) {
        os << "\t" << SEQ_label[SEQL_COMPLETED];
      }
      break;
    }

    case COMPLETION : {
      os << SEQ_label[SEQL_COMPLETION] << "\t"
         << (child[i] ? SEQ_label[SEQL_NON_TERMINAL] : SEQ_label[SEQL_TERMINAL]);
      break;
    }
    }

    os << endl;
  }

  if (nb_component > 0) {
    for (i = 0;i < nb_component;i++) {
      switch (stype[component[i][0]]) {
      case TRANSIENT :
        os << "\n"<< STAT_label[STATL_TRANSIENT] << " ";
        break;
      default :
        os << "\n" << STAT_label[STATL_RECURRENT] << " ";
        break;
      }
      os << STAT_label[STATL_CLASS] << "\t" << STAT_label[component_nb_state[i] == 1 ? STATL_STATE : STATL_STATES];

      for (j = 0;j < component_nb_state[i];j++) {
        os << "\t" << component[i][j];
      }

      if (stype[component[i][0]] == ABSORBING) {
        os << " (" << STAT_label[STATL_ABSORBING] << " " << STAT_label[STATL_STATE] << ")";
      }
    }
    os << endl;
  }

  os << SEQ_label[SEQL_MEMORY_TRANSITION_MATRIX] << endl;

  os << "\n";
  for (i = 1;i < nb_row;i++) {
    if ((type == ORDINARY) || (!child[i])) {
      for (j = order[i] - 1;j >= 0;j--) {
        os << state[i][j] << " ";
      }

      os << "\t";
      for (j = 0;j < nb_state;j++) {
        os << "\t";
        for (k = order[next[i][j]] - 1;k >= 0;k--) {
          os << state[next[i][j]][k] << " ";
        }
      }

      os << "\t";
      switch (memo_type[i]) {

      case NON_TERMINAL : {
        os << SEQ_label[SEQL_NON_TERMINAL];
        break;
      }

      case TERMINAL : {
        os << SEQ_label[SEQL_TERMINAL];
        if (child[i]) {
          os << "\t" << SEQ_label[SEQL_COMPLETED];
        }
        break;
      }

      case COMPLETION : {
        os << SEQ_label[SEQL_COMPLETION] << "\t"
           << (child[i] ? SEQ_label[SEQL_NON_TERMINAL] : SEQ_label[SEQL_TERMINAL]);
        break;
      }
      }

      os << endl;
    }
  }

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Default constructor of the VariableOrderMarkov class.
 */
/*--------------------------------------------------------------*/

VariableOrderMarkov::VariableOrderMarkov()

{
  nb_iterator = 0;
  markov_data = NULL;

  nb_output_process = 0;
  categorical_process = NULL;
  discrete_parametric_process = NULL;
  continuous_parametric_process = NULL;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Constructor of the VariableOrderMarkov class.
 *
 *  \param[in] itype     process type (ORDINARY/EQUILIBRIUM),
 *  \param[in] inb_state number of states,
 *  \param[in] inb_row   number of memories.
 */
/*--------------------------------------------------------------*/

VariableOrderMarkov::VariableOrderMarkov(process_type itype , int inb_state , int inb_row)
:VariableOrderMarkovChain(itype , inb_state , inb_row)

{
  nb_iterator = 0;
  markov_data = NULL;

  nb_output_process = 0;
  categorical_process = NULL;
  discrete_parametric_process = NULL;
  continuous_parametric_process = NULL;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Constructor of the VariableOrderMarkov class.
 *
 *  \param[in] itype      process type (ORDINARY/EQUILIBRIUM),
 *  \param[in] inb_state  number of states,
 *  \param[in] inb_row    number of memories,
 *  \param[in] imax_order maximum order.
 */
/*--------------------------------------------------------------*/

VariableOrderMarkov::VariableOrderMarkov(process_type itype , int inb_state ,
                                         int inb_row , int imax_order)
:VariableOrderMarkovChain(itype , inb_state , inb_row , imax_order)

{
  nb_iterator = 0;
  markov_data = NULL;

  nb_output_process = 0;
  categorical_process = NULL;
  discrete_parametric_process = NULL;
  continuous_parametric_process = NULL;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Construction of a VariableOrderMarkov object of fixed order.
 *
 *  \param[in] itype              process type (ORDINARY/EQUILIBRIUM),
 *  \param[in] inb_state          number of states,
 *  \param[in] iorder             order,
 *  \param[in] init_flag          flag initialization,
 *  \param[in  inb_output_process number of observation processes,
 *  \param[in  nb_value           number of observed values for the categorical observation process.
 */
/*--------------------------------------------------------------*/

VariableOrderMarkov::VariableOrderMarkov(process_type itype , int inb_state ,
                                         int iorder , bool init_flag ,
                                         int inb_output_process , int nb_value)
:VariableOrderMarkovChain(itype , inb_state , iorder , init_flag)

{
  nb_iterator = 0;
  markov_data = NULL;

  nb_output_process = inb_output_process;

  if (nb_output_process == 1) {
    categorical_process = new CategoricalSequenceProcess*[nb_output_process];
    categorical_process[0] = new CategoricalSequenceProcess(nb_state , nb_value , true);
  }
  else {
    categorical_process = NULL;
  }

  discrete_parametric_process = NULL;
  continuous_parametric_process = NULL;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Constructor of the VariableOrderMarkov class with completion of
 *         the memory tree.
 *
 *  \param[in] markov             reference on a VariableOrderMarkov object,
 *  \param[in] inb_output_process number of observation processes,
 *  \param[in] nb_value           number of observed values for the categorical observation process.
 */
/*--------------------------------------------------------------*/

VariableOrderMarkov::VariableOrderMarkov(const VariableOrderMarkov &markov ,
                                         int inb_output_process , int nb_value)

{
  memory_tree_completion(markov);

  nb_iterator = 0;

  if (markov.markov_data) {
    markov_data = new VariableOrderMarkovData(*(markov.markov_data) , false);
  }
  else {
    markov_data = NULL;
  }

  nb_output_process = inb_output_process;

  state_process = new CategoricalSequenceProcess(nb_state , nb_state , false);

  if (nb_output_process == 1) {
    categorical_process = new CategoricalSequenceProcess*[nb_output_process];
    categorical_process[1] = new CategoricalSequenceProcess(nb_state , nb_value , true);
  }
  else {
    categorical_process = NULL;
  }

  discrete_parametric_process = NULL;
  continuous_parametric_process = NULL;
}


/*--------------------------------------------------------------*/
/*
 *  \brief Constructor of the VariableOrderMarkov class with completion of
 *         the memory tree.
 *
 *  \param[in] markov             reference on a VariableOrderMarkov object,
 *  \param[in] inb_output_process number of observation processes,
 *  \param[in] nb_value           number of observed values for each observation process.
 */
/*--------------------------------------------------------------*/

/* VariableOrderMarkov::VariableOrderMarkov(const VariableOrderMarkov &markov ,
                                         int inb_output_process , int *nb_value)

{
  int i;


  memory_tree_completion(markov);

  nb_iterator = 0;

  if (markov.markov_data) {
    markov_data = new VariableOrderMarkovData(*(markov.markov_data) , false);
  }
  else {
    markov_data = NULL;
  }

  nb_output_process = inb_output_process;

  categorical_process = new CategoricalSequenceProcess*[nb_output_process];
  discrete_parametric_process = new DiscreteParametricProcess*[nb_output_process];
  continuous_parametric_process = new ContinuousParametricProcess*[nb_output_process];

  for (i = 0;i < nb_output_process;i++) {
    if (nb_value[i] == I_DEFAULT) {
      categorical_process[i] = NULL;
      discrete_parametric_process[i] = NULL;
      continuous_parametric_process[i] = new ContinuousParametricProcess(nb_state);
    }

    else if (nb_value[i] <= NB_OUTPUT) {
      categorical_process[i] = new CategoricalSequenceProcess(nb_state , nb_value[i] , true);
      discrete_parametric_process[i] = NULL;
      continuous_parametric_process[i] = NULL;
    }

    else {
      categorical_process[i] = NULL;
      discrete_parametric_process[i] = new DiscreteParametricProcess(nb_state , (int)(nb_value[i] * SAMPLE_NB_VALUE_COEFF));
      continuous_parametric_process[i] = NULL;
    }
  }
} */


/*--------------------------------------------------------------*/
/**
 *  \brief Constructor of the VariableOrderMarkov class with completion of
 *         the memory tree.
 *
 *  \param[in] pmarkov      pointer on a VariableOrderMarkovChain object,
 *  \param[in] pobservation pointer on a CategoricalProcess object,
 *  \param[in] length       sequence length.
 */
/*--------------------------------------------------------------*/

VariableOrderMarkov::VariableOrderMarkov(const VariableOrderMarkovChain *pmarkov ,
                                         const CategoricalProcess *pobservation , int length)

{
  build(*pmarkov);

  nb_iterator = 0;
  markov_data = NULL;

  nb_output_process = (pobservation ? 1 : 0);

  if (nb_output_process == 1) {
    categorical_process = new CategoricalSequenceProcess*[nb_output_process];
    categorical_process[0] = new CategoricalSequenceProcess(*pobservation);
  }
  else {
    categorical_process = NULL;
  }

  discrete_parametric_process = NULL;
  continuous_parametric_process = NULL;

  characteristic_computation(length , true);
}


/*--------------------------------------------------------------*/
/**
 *  \brief Copy of a VariableOrderMarkov object.
 *
 *  \param[in] markov    reference on a VariableOrderMarkov object.
 *  \param[in] data_flag flag copy of the included VariableOrderMarkovData object.
 */
/*--------------------------------------------------------------*/

void VariableOrderMarkov::copy(const VariableOrderMarkov &markov , bool data_flag)

{
  int i , j;


  nb_iterator = 0;

  if ((data_flag) && (markov.markov_data)) {
    markov_data = new VariableOrderMarkovData(*(markov.markov_data) , false);
  }
  else {
    markov_data = NULL;
  }

  nb_output_process = markov.nb_output_process;

  if (markov.categorical_process) {
    categorical_process = new CategoricalSequenceProcess*[nb_output_process];

    for (i = 0;i < nb_output_process;i++) {
      if (markov.categorical_process[i]) {
        categorical_process[i] = new CategoricalSequenceProcess(*(markov.categorical_process[i]));
      }
      else {
        categorical_process[i] = NULL;
      }
    }
  }

  else {
    categorical_process = NULL;
  }

  if (markov.discrete_parametric_process) {
    discrete_parametric_process = new DiscreteParametricProcess*[nb_output_process];

    for (i = 0;i < nb_output_process;i++) {
      if (markov.discrete_parametric_process[i]) {
        discrete_parametric_process[i] = new DiscreteParametricProcess(*(markov.discrete_parametric_process[i]));
      }
      else {
        discrete_parametric_process[i] = NULL;
      }
    }
  }

  else {
    discrete_parametric_process = NULL;
  }

  if (markov.continuous_parametric_process) {
    continuous_parametric_process = new ContinuousParametricProcess*[nb_output_process];

    for (i = 0;i < nb_output_process;i++) {
      if (markov.continuous_parametric_process[i]) {
        continuous_parametric_process[i] = new ContinuousParametricProcess(*(markov.continuous_parametric_process[i]));
      }
      else {
        continuous_parametric_process[i] = NULL;
      }
    }
  }

  else {
    continuous_parametric_process = NULL;
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Destruction of the data members of a VariableOrderMarkov object.
 */
/*--------------------------------------------------------------*/

void VariableOrderMarkov::remove()

{
  int i;


  delete markov_data;

  if (categorical_process) {
    for (i = 0;i < nb_output_process;i++) {
      delete categorical_process[i];
    }
    delete [] categorical_process;
  }

  if (discrete_parametric_process) {
    for (i = 0;i < nb_output_process;i++) {
      delete discrete_parametric_process[i];
    }
    delete [] discrete_parametric_process;
  }

  if (continuous_parametric_process) {
    for (i = 0;i < nb_output_process;i++) {
      delete continuous_parametric_process[i];
    }
    delete [] continuous_parametric_process;
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Destructor of the VariableOrderMarkov class.
 */
/*--------------------------------------------------------------*/

VariableOrderMarkov::~VariableOrderMarkov()

{
  remove();
}


/*--------------------------------------------------------------*/
/**
 *  \brief Destruction of a VariableOrderMarkov object taking account of
 *         the number of iterators pointing to it.
 */
/*--------------------------------------------------------------*/

void VariableOrderMarkov::conditional_delete()

{
  if (nb_iterator == 0) {
    delete this;
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Assignment operator of the VariableOrderMarkov class.
 *
 *  \param[in] markov reference on a VariableOrderMarkov object.
 *
 *  \return           VariableOrderMarkov object.
 */
/*--------------------------------------------------------------*/

VariableOrderMarkov& VariableOrderMarkov::operator=(const VariableOrderMarkov &markov)

{
  if ((&markov != this) && (nb_iterator == 0)) {
    remove();
    VariableOrderMarkovChain::remove();
    Chain::remove();

    Chain::copy(markov);
    VariableOrderMarkovChain::copy(markov);
    copy(markov);
  }

  return *this;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Extraction of a distribution.
 *
 *  \param[in] error     reference on a StatError object,
 *  \param[in] dist_type distribution type,
 *  \param[in] variable  variable index,
 *  \param[in] value     state or observation.
 *
 *  \return              DiscreteParametricModel object.
 */
/*--------------------------------------------------------------*/

DiscreteParametricModel* VariableOrderMarkov::extract(StatError &error , process_distribution dist_type ,
                                                      int variable , int value) const

{
  bool status = true;
  int hvariable;
  Distribution *pdist;
  DiscreteParametric *pparam;
  DiscreteParametricModel *dist;
  FrequencyDistribution *phisto;
  CategoricalSequenceProcess *process;


  dist = NULL;
  error.init();

  pdist = NULL;
  pparam = NULL;

  if (dist_type == OBSERVATION) {
    if ((variable < 1) || (variable > nb_output_process)) {
      status = false;
      error.update(STAT_error[STATR_OUTPUT_PROCESS_INDEX]);
    }

    else {
      if ((value < 0) || (value >= nb_state)) {
        status = false;
        ostringstream error_message;
        error_message << STAT_label[STATL_STATE] << " " << value << " "
                      << STAT_error[STATR_NOT_PRESENT];
        error.update((error_message.str()).c_str());
      }

      else {
        if (categorical_process[variable - 1]) {
          pdist = categorical_process[variable - 1]->observation[value];
        }
        else if (discrete_parametric_process[variable - 1]) {
          pparam = discrete_parametric_process[variable - 1]->observation[value];
        }
        else {
          status = false;
          ostringstream correction_message;
          correction_message << STAT_label[STATL_CATEGORICAL] << " or "
                             << STAT_label[STATL_DISCRETE_PARAMETRIC];
          error.correction_update(STAT_error[STATR_OUTPUT_PROCESS_TYPE] , (correction_message.str()).c_str());
        }
      }
    }
  }

  else {
    if ((variable < 0) || (variable > nb_output_process)) {
      status = false;
      error.update(STAT_error[STATR_OUTPUT_PROCESS_INDEX]);
    }

    else {
      if (variable == 0) {
        process = state_process;
      }

      else {
        process = categorical_process[variable - 1];

        if (!process) {
          status = false;
          ostringstream error_message;
          error_message << STAT_label[STATL_VARIABLE] << " " << variable << ": " 
                        << SEQ_error[SEQR_CHARACTERISTICS_NOT_COMPUTED];
          error.update((error_message.str()).c_str());
        }
      }

      if ((process) && ((value < 0) || (value >= process->nb_value))) {
        status = false;
        ostringstream error_message;
        error_message << STAT_label[variable == 0 ? STATL_STATE : STATL_OUTPUT] << " "
                      << value << " " << STAT_error[STATR_NOT_PRESENT];
        error.update((error_message.str()).c_str());
      }

      if (status) {
        switch (dist_type) {
        case FIRST_OCCURRENCE :
          pdist = process->first_occurrence[value];
          break;
        case RECURRENCE_TIME :
          pdist = process->recurrence_time[value];
          break;
        case SOJOURN_TIME :
          pparam = process->sojourn_time[value];
          break;
        case NB_RUN :
          pdist = process->nb_run[value];
          break;
        case NB_OCCURRENCE :
          pdist = process->nb_occurrence[value];
          break;
        }

        if ((!pdist) && (!pparam)) {
          status = false;
          error.update(SEQ_error[SEQR_NON_EXISTING_CHARACTERISTIC_DISTRIBUTION]);
        }
      }
    }
  }

  if (status) {
    phisto = NULL;

    if (markov_data) {
      switch (markov_data->type[0]) {
      case STATE :
        hvariable = variable;
        break;
      case INT_VALUE :
        hvariable = variable - 1;
        break;
      }

      if (hvariable >= 0) {
        switch (dist_type) {

        case OBSERVATION : {
          if ((markov_data->observation_distribution) &&
              (markov_data->observation_distribution[hvariable])) {
            phisto = markov_data->observation_distribution[hvariable][value];
          }
          break;
        }

        case FIRST_OCCURRENCE : {
          phisto = markov_data->characteristics[hvariable]->first_occurrence[value];
          break;
        }

        case RECURRENCE_TIME : {
          if (markov_data->characteristics[hvariable]->recurrence_time[value]->nb_element > 0) {
            phisto = markov_data->characteristics[hvariable]->recurrence_time[value];
          }
          break;
        }

        case SOJOURN_TIME : {
          if (markov_data->characteristics[hvariable]->sojourn_time[value]->nb_element > 0) {
            phisto = markov_data->characteristics[hvariable]->sojourn_time[value];
          }
          break;
        }

        case NB_RUN : {
          phisto = markov_data->characteristics[hvariable]->nb_run[value];
          break;
        }

        case NB_OCCURRENCE : {
          phisto = markov_data->characteristics[hvariable]->nb_occurrence[value];
          break;
        }
        }
      }
    }

    if (pdist) {
      dist = new DiscreteParametricModel(*pdist , phisto);
    }
    else if (pparam) {
      dist = new DiscreteParametricModel(*pparam , phisto);
    }
  }

  return dist;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Extraction of the VariableOrderMarkovData object included in
 *         a VariableOrderMarkov object.
 *
 *  \param[in] error reference on a StatError object.
 *
 *  \return          VariableOrderMarkovData object.
 */
/*--------------------------------------------------------------*/

VariableOrderMarkovData* VariableOrderMarkov::extract_data(StatError &error) const

{
  bool status = true;
  VariableOrderMarkovData *seq;


  seq = NULL;
  error.init();

  if (!markov_data) {
    status = false;
    error.update(STAT_error[STATR_NO_DATA]);
  }
  else if (nb_output_process + 1 != markov_data->nb_variable) {
    status = false;
    error.update(SEQ_error[SEQR_STATE_SEQUENCES]);
  }

  if (status) {
    seq = new VariableOrderMarkovData(*markov_data);
    seq->markov = new VariableOrderMarkov(*this , false);
  }

  return seq;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Application of a threshold on the probability parameters of a variable-order Markov chain.
 *
 *  \param[in] min_probability minimum probability.
 *
 *  \return                    VariableOrderMarkov object.
 */
/*--------------------------------------------------------------*/

VariableOrderMarkov* VariableOrderMarkov::thresholding(double min_probability) const

{
  int i;
  VariableOrderMarkov *markov;


  markov = new VariableOrderMarkov(*this , false);
  markov->VariableOrderMarkovChain::thresholding(min_probability);

  for (i = 0;i < markov->nb_output_process;i++) {
    if (markov->categorical_process[i]) {
      markov->categorical_process[i]->thresholding(min_probability);
    }
  }

  return markov;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Construction of a VariableOrderMarkov object from a file.
 *
 *  \param[in] error  reference on a StatError object,
 *  \param[in] path   file path,
 *  \param[in] length sequence length.
 *
 *  \return           VariableOrderMarkov object.
 */
/*--------------------------------------------------------------*/

VariableOrderMarkov* VariableOrderMarkov::ascii_read(StatError &error ,
                                                     const string path , int length)

{
  string buffer;
  size_t position;
  typedef tokenizer<char_separator<char>> tokenizer;
  char_separator<char> separator(" \t");
  process_type type = DEFAULT_TYPE;
  bool status;
  int i;
  int line;
  const VariableOrderMarkovChain *imarkov;
  const CategoricalProcess *observation;
  VariableOrderMarkov *markov;
  ifstream in_file(path.c_str());


  markov = NULL;
  error.init();

  if (!in_file) {
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {
    status = true;
    line = 0;

    if (length < 2) {
      status = false;
      error.update(SEQ_error[SEQR_SHORT_SEQUENCE_LENGTH]);
    }
    if (length > MAX_LENGTH) {
      status = false;
      error.update(SEQ_error[SEQR_LONG_SEQUENCE_LENGTH]);
    }

    while (getline(in_file , buffer)) {
      line++;

#     ifdef DEBUG
      cout << line << "  " << buffer << endl;
#     endif

      position = buffer.find('#');
      if (position != string::npos) {
        buffer.erase(position);
      }
      i = 0;

      tokenizer tok_buffer(buffer , separator);

      for (tokenizer::iterator token = tok_buffer.begin();token != tok_buffer.end();token++) {

        // test (EQUILIBRIUM_)MARKOV_CHAIN keyword

        if (i == 0) {
          if (*token == SEQ_word[SEQW_MARKOV_CHAIN]) {
            type = ORDINARY;
          }
          else if (*token == SEQ_word[SEQW_EQUILIBRIUM_MARKOV_CHAIN]) {
            type = EQUILIBRIUM;
          }
          else {
            status = false;
            ostringstream correction_message;
            correction_message << SEQ_word[SEQW_MARKOV_CHAIN] << " or "
                               << SEQ_word[SEQW_EQUILIBRIUM_MARKOV_CHAIN];
            error.correction_update(STAT_parsing[STATP_KEYWORD] ,
                                    (correction_message.str()).c_str() , line);
          }
        }

        i++;
      }

      if (i > 0) {
        if (i != 1) {
          status = false;
          error.update(STAT_parsing[STATP_FORMAT] , line);
        }
        break;
      }
    }

    if (type != DEFAULT_TYPE) {

      // analysis of the format and reading of the variable-order Markov chain

      imarkov = VariableOrderMarkovChain::parsing(error , in_file , line , type);

      if (imarkov) {

        // analysis of the format and reading of the categorical observation distributions

        observation = NULL;

        while (getline(in_file , buffer)) {
          line++;

#         ifdef DEBUG
          cout << line << "  " << buffer << endl;
#         endif

          position = buffer.find('#');
          if (position != string::npos) {
            buffer.erase(position);
          }
          i = 0;

          tokenizer tok_buffer(buffer , separator);

          for (tokenizer::iterator token = tok_buffer.begin();token != tok_buffer.end();token++) {

            // test OUTPUT_PROCESS keyword

            if (i == 0) {
              if (*token != STAT_word[STATW_OUTPUT_PROCESS]) {
                status = false;
                error.correction_update(STAT_parsing[STATP_KEYWORD] , STAT_word[STATW_OUTPUT_PROCESS] , line);
              }
            }

            i++;
          }

          if (i > 0) {
            if (i != 1) {
              status = false;
              error.update(STAT_parsing[STATP_FORMAT] , line);
            }

            observation = CategoricalProcess::parsing(error , in_file , line ,
                                                      ((Chain*)imarkov)->nb_state ,
                                                      HIDDEN_MARKOV , false);
            if (!observation) {
              status = false;
            }

            break;
          }
        }

        while (getline(in_file , buffer)) {
          line++;

#         ifdef DEBUG
          cout << line << "  " << buffer << endl;
#         endif

          position = buffer.find('#');
          if (position != string::npos) {
            buffer.erase(position);
          }
          if (!(trim_right_copy_if(buffer , is_any_of(" \t")).empty())) {
            status = false;
            error.update(STAT_parsing[STATP_FORMAT] , line);
          }
        }

        if (status) {
          markov = new VariableOrderMarkov(imarkov , observation , length);

#         ifdef DEBUG
          imarkov->ascii_memory_tree_print(cout);
          markov->ascii_write(cout);
#         endif

        }

        delete imarkov;
        delete observation;
      }
    }
  }

  return markov;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing on a single line of a VariableOrderMarkov object.
 *
 *  \param[in,out] os stream.
 */
/*--------------------------------------------------------------*/

ostream& VariableOrderMarkov::line_write(ostream &os) const

{
  os << nb_state << " " << STAT_word[STATW_STATES] << "   "
     << SEQ_label[SEQL_MAX_ORDER] << " " << max_order;

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a VariableOrderMarkov object and the associated data structure.
 *
 *  \param[in,out] os         stream,
 *  \param[in]     seq        pointer on a VariableOrderMarkovData object,
 *  \param[in]     exhaustive flag detail level,
 *  \param[in]     file_flag  flag file,
 *  \param[in]     hidden     flag hidden model.
 */
/*--------------------------------------------------------------*/

ostream& VariableOrderMarkov::ascii_write(ostream &os , const VariableOrderMarkovData *seq ,
                                          bool exhaustive , bool file_flag , bool hidden) const

{
  bool **logic_transition;
  int i , j , k;
  int buff , variable , max_memory_count , *memory_count , width[2];
  double standard_normal_value , half_confidence_interval , **distance;
  FrequencyDistribution *marginal_dist = NULL , **observation_dist = NULL;
  Histogram *marginal_histo = NULL , **observation_histo = NULL;
  SequenceCharacteristics *characteristics = NULL;
  ios_base::fmtflags format_flags;


  format_flags = os.setf(ios::left , ios::adjustfield);

  if (hidden) {
    switch (type) {
    case ORDINARY :
      os << SEQ_word[SEQW_HIDDEN_MARKOV_CHAIN] << endl;
      break;
    case EQUILIBRIUM :
      os << SEQ_word[SEQW_EQUILIBRIUM_HIDDEN_MARKOV_CHAIN] << endl;
      break;
    }
  }

  else {
    switch (type) {
    case ORDINARY :
      os << SEQ_word[SEQW_MARKOV_CHAIN] << endl;
      break;
    case EQUILIBRIUM :
      os << SEQ_word[SEQW_EQUILIBRIUM_MARKOV_CHAIN] << endl;
      break;
    }
  }

  // writing of the variable-order Markov chain parameters

  ascii_print(os , file_flag);

//  if ((nb_component == 1) && (seq) && ((!hidden) || (seq->type[0] == STATE))) {
  if ((seq) && ((!hidden) || (seq->type[0] == STATE))) {
    normal dist;
    standard_normal_value = quantile(complement(dist , 0.025));

    if (!hidden) {
      width[0] = 0;
      for (i = 1;i < nb_row;i++) {
        if ((memo_type[i] == TERMINAL) || ((type == ORDINARY) &&
             (memo_type[i] == NON_TERMINAL))) {
          buff = column_width(nb_state , transition[i]);
          if (buff > width[0]) {
            width[0] = buff;
          }
        }
      }
      width[0]++;
    }

    memory_count = new int[nb_row];
    max_memory_count = 0;
    for (i = 1;i < nb_row;i++) {
//      if (memo_type[i] == TERMINAL) {
      if ((memo_type[i] == TERMINAL) || ((type == ORDINARY) &&
           (memo_type[i] == NON_TERMINAL)) || (hidden)) {
        memory_count[i] = 0;
        for (j = 0;j < nb_state;j++) {
          memory_count[i] += seq->chain_data->transition[i][j];
        }

        if (memory_count[i] > max_memory_count) {
          max_memory_count = memory_count[i];
        }
      }
    }
    width[1] = column_width(max_memory_count) + ASCII_SPACE;

    os << "\n";
    if (file_flag) {
      os << "# ";
    }
    if (hidden) {
      os << SEQ_label[SEQL_TRANSITION_COUNTS] << endl;
    }
    else {
      os << SEQ_label[SEQL_TRANSITION_PROBABILITIY_CONFIDENCE_INTERVAL] << endl;
    }

    os << "\n";
    for (i = 1;i < nb_row;i++) {
//      if (memo_type[i] == TERMINAL) {
      if ((memo_type[i] == TERMINAL) || ((type == ORDINARY) &&
           (memo_type[i] == NON_TERMINAL)) || (hidden)) {
        if (memory_count[i] > 0.) {
          if (file_flag) {
            os << "# ";
          }

          if (hidden) {
            for (j = 0;j < nb_state;j++) {
              os << setw(width[1]) << seq->chain_data->transition[i][j];
            }
            os << "  ";
          }

          else {
            if (memo_type[i] == TERMINAL) {
              for (j = 0;j < nb_state;j++) {
                if ((transition[i][j] > 0.) && (transition[i][j] < 1.)) {
                  half_confidence_interval = standard_normal_value *
                                             sqrt(transition[i][j] * (1. - transition[i][j]) / memory_count[i]);
                  os << setw(width[0]) << MAX(transition[i][j] - half_confidence_interval , 0.)
                     << setw(width[0]) << MIN(transition[i][j] + half_confidence_interval , 1.)
                     << "| ";
                }
                else {
                  os << setw(width[0]) << " "
                     << setw(width[0]) << " "
                     << "| ";
                }
              }
            }

            else {
              for (j = 0;j < nb_state;j++) {
                os << setw(width[0]) << " "
                   << setw(width[0]) << " "
                   << "| ";
              }
            }
          }

          os << setw(width[1]) << memory_count[i] << "  ";

          for (j = max_order - 1;j >= order[i];j--) {
            os << "  ";
          }
          for (j = order[i] - 1;j >= 0;j--) {
            os << state[i][j] << " ";
          }

          os << endl;
        }
      }
    }

    delete [] memory_count;
  }

  if ((seq) && (seq->type[0] == STATE)) {
    characteristics = seq->characteristics[0];
  }
  else {
    characteristics = NULL;
  }

  state_process->ascii_print(os , 0 , NULL , NULL , characteristics ,
                             exhaustive , file_flag);

  if (hidden) {
    for (i = 0;i < nb_output_process;i++) {
      if (discrete_parametric_process[i]) {
        if (discrete_parametric_process[i]->weight) {
          width[0] = column_width(nb_state , discrete_parametric_process[i]->weight->mass);
        }
        else {
          width[0] = 0;
        }
        if (discrete_parametric_process[i]->restoration_weight) {
          buff = column_width(nb_state , discrete_parametric_process[i]->restoration_weight->mass);
          if (buff > width[0]) {
            width[0] = buff;
          }
        }
        width[0]++;

        if (discrete_parametric_process[i]->weight) {
          os << "\n";
          if (file_flag) {
            os << "# ";
          }
          os << STAT_label[STATL_THEORETICAL] << " " << SEQ_label[SEQL_STATE_PROBABILITY] << ": ";

          for (j = 0;j < nb_state;j++) {
            os << setw(width[0]) << discrete_parametric_process[i]->weight->mass[j];
          }
          os << endl;
        }

        if (discrete_parametric_process[i]->restoration_weight) {
          if (file_flag) {
            os << "# ";
          }
          os << STAT_label[STATL_RESTORATION] << " " << SEQ_label[SEQL_STATE_PROBABILITY] << ": ";

          for (j = 0;j < nb_state;j++) {
            os << setw(width[0]) << discrete_parametric_process[i]->restoration_weight->mass[j];
          }
          os << endl;
        }
        break;
      }

      else if (continuous_parametric_process[i]) {
        if (continuous_parametric_process[i]->weight) {
          width[0] = column_width(nb_state , continuous_parametric_process[i]->weight->mass);
        }
        else {
          width[0] = 0;
        }
        if (continuous_parametric_process[i]->restoration_weight) {
          buff = column_width(nb_state , continuous_parametric_process[i]->restoration_weight->mass);
          if (buff > width[0]) {
            width[0] = buff;
          }
        }
        width[0]++;

        if (continuous_parametric_process[i]->weight) {
          os << "\n";
          if (file_flag) {
            os << "# ";
          }
          os << STAT_label[STATL_THEORETICAL] << " " << SEQ_label[SEQL_STATE_PROBABILITY] << ": ";

          for (j = 0;j < nb_state;j++) {
            os << setw(width[0]) << continuous_parametric_process[i]->weight->mass[j];
          }
          os << endl;
        }

        if (continuous_parametric_process[i]->restoration_weight) {
          if (file_flag) {
            os << "# ";
          }
          os << STAT_label[STATL_RESTORATION] << " " << SEQ_label[SEQL_STATE_PROBABILITY] << ": ";

          for (j = 0;j < nb_state;j++) {
            os << setw(width[0]) << continuous_parametric_process[i]->restoration_weight->mass[j];
          }
          os << endl;
        }
        break;
      }
    }

    os << "\n" << nb_output_process << " "
       << STAT_word[nb_output_process == 1 ? STATW_OUTPUT_PROCESS : STATW_OUTPUT_PROCESSES] << endl;
  }

  // writing of the distributions associated with each observation process

  if (hidden) {
    logic_transition = logic_transition_computation();

    distance = new double*[nb_state];
    for (i = 0;i < nb_state;i++) {
      distance[i] = new double[nb_state];
    }
  }

  for (i = 0;i < nb_output_process;i++) {
    os << "\n" << STAT_word[STATW_OUTPUT_PROCESS];

    if (hidden) {
      os << " " << i + 1;

      if (categorical_process[i]) {
        os << " : " << STAT_word[STATW_CATEGORICAL];
      }
      else if (discrete_parametric_process[i]) {
        os << " : " << STAT_word[STATW_DISCRETE_PARAMETRIC];
      }
      else {
        os << " : " << STAT_word[STATW_CONTINUOUS_PARAMETRIC];
      }
    }
    os << endl;

    if (seq) {
      switch (seq->type[0]) {
      case STATE :
        variable = i + 1;
        break;
      default :
        variable = i;
        break;
      }

      if (seq->observation_distribution) {
        observation_dist = seq->observation_distribution[variable];
      }
      marginal_dist = seq->marginal_distribution[variable];

      if (seq->observation_histogram) {
        observation_histo = seq->observation_histogram[variable];
      }
      marginal_histo = seq->marginal_histogram[variable];

      characteristics = seq->characteristics[variable];
    }

    if (categorical_process[i]) {
      categorical_process[i]->ascii_print(os , i + 1 , observation_dist , marginal_dist ,
                                          characteristics , exhaustive , file_flag);

      if (hidden) {
        for (j = 0;j < nb_state;j++) {
          distance[j][j] = 0.;

          for (k = j + 1;k < nb_state;k++) {
            if ((logic_transition[j][k]) || (logic_transition[k][j])) {
              distance[j][k] = categorical_process[i]->observation[j]->overlap_distance_computation(*(categorical_process[i]->observation[k]));
            }
            else {
              distance[j][k] = 1.;
            }

            distance[k][j] = distance[j][k];
          }
        }
      }
    }

    else if (discrete_parametric_process[i]) {
      discrete_parametric_process[i]->ascii_print(os , observation_dist , marginal_dist ,
                                                  exhaustive , file_flag);

      if (hidden) {
        for (j = 0;j < nb_state;j++) {
          distance[j][j] = 0.;

          for (k = j + 1;k < nb_state;k++) {
            if ((logic_transition[j][k]) || (logic_transition[k][j])) {
              distance[j][k] = discrete_parametric_process[i]->observation[j]->sup_norm_distance_computation(*(discrete_parametric_process[i]->observation[k]));
            }
            else {
              distance[j][k] = 1.;
            }

            distance[k][j] = distance[j][k];
          }
        }
      }
    }

    else {
      continuous_parametric_process[i]->ascii_print(os , observation_histo , observation_dist ,
                                                    marginal_histo , marginal_dist ,
                                                    exhaustive , file_flag);

      if (hidden) {
        for (j = 0;j < nb_state;j++) {
          distance[j][j] = 0.;

          for (k = j + 1;k < nb_state;k++) {
            if ((logic_transition[j][k]) || (logic_transition[k][j])) {
              distance[j][k] = continuous_parametric_process[i]->observation[j]->sup_norm_distance_computation(*(continuous_parametric_process[i]->observation[k]));
            }
            else {
              distance[j][k] = 1.;
            }

            distance[k][j] = distance[j][k];
          }
        }
      }
    }

    if (hidden) {
      width[0] = column_width(nb_state , distance[0]);
      for (j = 1;j < nb_state;j++) {
        buff = column_width(nb_state , distance[j]);
        if (buff > width[0]) {
          width[0] = buff;
        }
      }
      width[0] += ASCII_SPACE;

      os.setf(ios::left , ios::adjustfield);

      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << STAT_label[STATL_CONSECUTIVE_STATE_OBSERVATION_DISTRIBUTION_DISTANCE] << endl;

      for (j = 0;j < nb_state;j++) {
        if (file_flag) {
          os << "# ";
        }
        for (k = 0;k < nb_state;k++) {
          if ((k != j) && (logic_transition[j][k])) {
            os << setw(width[0]) << distance[j][k];
          }
          else {
            os << setw(width[0]) << "_";
          }
        }
        os << endl;
      }
    }
  }

  if (hidden) {
    for (i = 0;i < nb_state;i++) {
      delete [] logic_transition[i];
    }
    delete [] logic_transition;

    for (i = 0;i < nb_state;i++) {
      delete [] distance[i];
    }
    delete [] distance;
  }

  if (seq) {
    int nb_parameter = nb_parameter_computation(hidden ? MIN_PROBABILITY : 0.) , nb_transient_parameter;
    double information;


    // writing of the sequence length frequency distribution

    os << "\n";
    if (file_flag) {
      os << "# ";
    }
    os << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";
    seq->length_distribution->ascii_characteristic_print(os , false , file_flag);

    if (exhaustive) {
      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << "   | " << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << endl;
      seq->length_distribution->ascii_print(os , file_flag);
    }

    os << "\n";
    if (file_flag) {
      os << "# ";
    }
    os << SEQ_label[SEQL_CUMUL_LENGTH] << ": " << seq->cumul_length << endl;

    // writing of the information quantity of the observed sequences in the i.i.d. case

    for (i = 0;i < seq->nb_variable;i++) {
      if (seq->type[i] == REAL_VALUE) {
        break;
      }
    }

    if (i == seq->nb_variable) {
      information = seq->iid_information_computation();

      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << SEQ_label[SEQL_IID_INFORMATION] << ": " << information << " ("
         << information / seq->cumul_length << ")" << endl;
    }

    // writing of the (penalized) log-likelihoods of the model for sequences

    if (hidden) {
      if (seq->restoration_likelihood != D_INF) {
        os << "\n";
        if (file_flag) {
          os << "# ";
        }
        os << SEQ_label[SEQL_STATE_SEQUENCES_LIKELIHOOD] << ": " << seq->restoration_likelihood << "   ("
           << STAT_label[STATL_NORMALIZED] << ": " << seq->restoration_likelihood / seq->cumul_length << ")" << endl;
      }

      if (seq->sample_entropy != D_DEFAULT) {
        os << "\n";
        if (file_flag) {
          os << "# ";
        }
        os << SEQ_label[SEQL_STATE_SEQUENCE_ENTROPY] << ": " << seq->sample_entropy << "   ("
           << STAT_label[STATL_NORMALIZED] << ": " << seq->sample_entropy / seq->cumul_length << ")" << endl;
      }

      if (seq->likelihood != D_INF) {
        os << "\n";
        if (file_flag) {
          os << "# ";
        }
        os << SEQ_label[SEQL_OBSERVED_SEQUENCES_LIKELIHOOD] << ": " << seq->likelihood << "   ("
           << STAT_label[STATL_NORMALIZED] << ": " << seq->likelihood / seq->cumul_length << ")" << endl;
      }
    }

    else {
      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << STAT_label[STATL_LIKELIHOOD] << ": " << seq->likelihood << "   ("
         << STAT_label[STATL_NORMALIZED] << ": " << seq->likelihood / seq->cumul_length << ")" << endl;
    }

    if (seq->likelihood != D_INF) {
      if (type == ORDINARY) {
        nb_transient_parameter = nb_transient_parameter_computation(hidden ? MIN_PROBABILITY : 0.);

        os << "\n";
        if (file_flag) {
          os << "# ";
        }
        os << nb_transient_parameter << " "
           << SEQ_label[nb_transient_parameter == 1 ? SEQL_FREE_TRANSIENT_PARAMETER : SEQL_FREE_TRANSIENT_PARAMETERS] << endl;
      }

      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << nb_parameter << " " << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS]
         << "   2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[AIC] << "): "
         << 2 * (seq->likelihood - nb_parameter) << endl;

      if ((type == ORDINARY) && (nb_transient_parameter > 0)) {
        if (file_flag) {
          os << "# ";
        }
        os << nb_transient_parameter + nb_parameter << " " << STAT_label[STATL_FREE_PARAMETERS]
           << "   2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[AIC] << "): "
           << 2 * (seq->likelihood - nb_transient_parameter - nb_parameter) << endl;
      }

      if (nb_parameter < seq->cumul_length - 1) {
        os << "\n";
        if (file_flag) {
          os << "# ";
        }
        os << nb_parameter << " " << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS]
           << "   2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[AICc] << "): "
           << 2 * (seq->likelihood - (double)(nb_parameter * seq->cumul_length) /
               (double)(seq->cumul_length - nb_parameter - 1)) << endl;
      }

      if ((type == ORDINARY) && (nb_transient_parameter > 0) &&
          (nb_transient_parameter + nb_parameter < seq->cumul_length - 1)) {
        if (file_flag) {
          os << "# ";
        }
        os << nb_transient_parameter + nb_parameter << " " << STAT_label[STATL_FREE_PARAMETERS]
           << "   2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[AICc] << "): "
           << 2 * (seq->likelihood - (double)((nb_transient_parameter + nb_parameter) * seq->cumul_length) /
               (double)(seq->cumul_length - nb_transient_parameter - nb_parameter - 1)) << endl;
      }

      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << nb_parameter << " " << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS]
         << "   2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[BIC] << "): "
         << 2 * seq->likelihood - nb_parameter * log((double)seq->cumul_length) << endl;

      if ((type == ORDINARY) && (nb_transient_parameter > 0)) {
        if (file_flag) {
          os << "# ";
        }
        os << nb_transient_parameter + nb_parameter << " " << STAT_label[STATL_FREE_PARAMETERS]
           << "   2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[BIC] << "): "
           << 2 * seq->likelihood - (nb_transient_parameter + nb_parameter) * log((double)seq->cumul_length) << endl;
      }

      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << nb_parameter + (type == ORDINARY ? nb_transient_parameter : 0) << " " << STAT_label[STATL_FREE_PARAMETERS]
         << "   2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[BICc] << "): "
         << 2 * seq->likelihood - penalty_computation(hidden , (hidden ? MIN_PROBABILITY : 0.)) << endl;
    }

    if ((hidden) && (seq->likelihood != D_INF)) {
      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << nb_parameter << " " << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS]
         << "   2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[ICL] << "): "
         << 2 * (seq->likelihood - seq->sample_entropy) - nb_parameter * log((double)seq->cumul_length) << endl;

      if ((type == ORDINARY) && (nb_transient_parameter > 0)) {
        if (file_flag) {
          os << "# ";
        }
        os << nb_transient_parameter + nb_parameter << " " << STAT_label[STATL_FREE_PARAMETERS]
           << "   2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[ICL] << "): "
           << 2 * (seq->likelihood - seq->sample_entropy) - (nb_transient_parameter + nb_parameter) * log((double)seq->cumul_length) << endl;
      }

      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << nb_parameter + (type == ORDINARY ? nb_transient_parameter : 0) << " " << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS]
         << "   2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[ICLc] << "): "
         << 2 * (seq->likelihood - seq->sample_entropy) - penalty_computation(hidden , MIN_PROBABILITY) << endl;
    }
  }

  os.setf(format_flags , ios::adjustfield);

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a VariableOrderMarkov object.
 *
 *  \param[in,out] os         stream,
 *  \param[in]     exhaustive flag detail level.
 */
/*--------------------------------------------------------------*/

ostream& VariableOrderMarkov::ascii_write(ostream &os , bool exhaustive) const

{
  return ascii_write(os , markov_data , exhaustive , false);
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a VariableOrderMarkov object in a file.
 *
 *  \param[in] error      reference on a StatError object,
 *  \param[in] path       file path,
 *  \param[in] exhaustive flag detail level.
 *
 *  \return               error status.
 */
/*--------------------------------------------------------------*/

bool VariableOrderMarkov::ascii_write(StatError &error , const string path ,
                                      bool exhaustive) const

{
  bool status;
  ofstream out_file(path.c_str());


  error.init();

  if (!out_file) {
    status = false;
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {
    status = true;
    ascii_write(out_file , markov_data , exhaustive , true);
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a VariableOrderMarkov object and the associated data structure
 *         in a file at the spreadsheet format.
 *
 *  \param[in,out] os     stream,
 *  \param[in]     seq    pointer on a VariableOrderMarkovData object,
 *  \param[in]     hidden flag hidden model.
 */
/*--------------------------------------------------------------*/

ostream& VariableOrderMarkov::spreadsheet_write(ostream &os ,
                                                const VariableOrderMarkovData *seq ,
                                                bool hidden) const

{
  bool **logic_transition;
  int i , j , k;
  int variable;
  double **distance;
  FrequencyDistribution *marginal_dist = NULL , **observation_dist = NULL;
  Histogram *marginal_histo = NULL , **observation_histo = NULL;
  SequenceCharacteristics *characteristics = NULL;


  if (hidden) {
    switch (type) {
    case ORDINARY :
      os << SEQ_word[SEQW_HIDDEN_MARKOV_CHAIN] << endl;
      break;
    case EQUILIBRIUM :
      os << SEQ_word[SEQW_EQUILIBRIUM_HIDDEN_MARKOV_CHAIN] << endl;
      break;
    }
  }

  else {
    switch (type) {
    case ORDINARY :
      os << SEQ_word[SEQW_MARKOV_CHAIN] << endl;
      break;
    case EQUILIBRIUM :
      os << SEQ_word[SEQW_EQUILIBRIUM_MARKOV_CHAIN] << endl;
      break;
    }
  }

  // writing of the variable-order Markov chain parameters

  spreadsheet_print(os);

  if ((seq) && (seq->type[0] == STATE)) {
    characteristics = seq->characteristics[0];
  }
  else {
    characteristics = NULL;
  }

  state_process->spreadsheet_print(os , 0 , NULL , NULL , characteristics);

  // writing of the distributions associated with each observation process

  if (hidden) {
    os << "\n" << nb_output_process << "\t"
       << STAT_word[nb_output_process == 1 ? STATW_OUTPUT_PROCESS : STATW_OUTPUT_PROCESSES] << endl;
  }

  if (hidden) {
    logic_transition = logic_transition_computation();

    distance = new double*[nb_state];
    for (i = 0;i < nb_state;i++) {
      distance[i] = new double[nb_state];
    }
  }

  for (i = 0;i < nb_output_process;i++) {
    os << "\n" << STAT_word[STATW_OUTPUT_PROCESS];

    if (hidden) {
      os << "\t" << i + 1;

      if (categorical_process[i]) {
        os << "\t" << STAT_word[STATW_CATEGORICAL];
      }
      else if (discrete_parametric_process[i]) {
        os << "\t" << STAT_word[STATW_DISCRETE_PARAMETRIC];
      }
      else {
        os << "\t" << STAT_word[STATW_CONTINUOUS_PARAMETRIC];
      }
    }
    os << endl;

    if (seq) {
      switch (seq->type[0]) {
      case STATE :
        variable = i + 1;
        break;
      default :
        variable = i;
        break;
      }

      if (seq->observation_distribution) {
        observation_dist = seq->observation_distribution[variable];
      }
      marginal_dist = seq->marginal_distribution[variable];

      if (seq->observation_histogram) {
        observation_histo = seq->observation_histogram[variable];
      }
      marginal_histo = seq->marginal_histogram[variable];

      characteristics = seq->characteristics[variable];
    }

    if (categorical_process[i]) {
      categorical_process[i]->spreadsheet_print(os , i + 1 , observation_dist , marginal_dist ,
                                                characteristics);

      if (hidden) {
        for (j = 0;j < nb_state;j++) {
          for (k = j + 1;k < nb_state;k++) {
            if ((logic_transition[j][k]) || (logic_transition[k][j])) {
              distance[j][k] = categorical_process[i]->observation[j]->overlap_distance_computation(*(categorical_process[i]->observation[k]));
              distance[k][j] = distance[j][k];
            }
          }
        }
      }
    }

    else if (discrete_parametric_process[i]) {
      discrete_parametric_process[i]->spreadsheet_print(os , observation_dist , marginal_dist);

      if (hidden) {
        for (j = 0;j < nb_state;j++) {
          for (k = j + 1;k < nb_state;k++) {
            if ((logic_transition[j][k]) || (logic_transition[k][j])) {
              distance[j][k] = discrete_parametric_process[i]->observation[j]->sup_norm_distance_computation(*(discrete_parametric_process[i]->observation[k]));
              distance[k][j] = distance[j][k];
            }
          }
        }
      }
    }

    else {
      continuous_parametric_process[i]->spreadsheet_print(os , observation_histo , observation_dist ,
                                                          marginal_histo , marginal_dist);

      if (hidden) {
        for (j = 0;j < nb_state;j++) {
          for (k = j + 1;k < nb_state;k++) {
            if ((logic_transition[j][k]) || (logic_transition[k][j])) {
              distance[j][k] = continuous_parametric_process[i]->observation[j]->sup_norm_distance_computation(*(continuous_parametric_process[i]->observation[k]));
              distance[k][j] = distance[j][k];
            }
          }
        }
      }
    }

    if (hidden) {
      os << "\n" << STAT_label[STATL_CONSECUTIVE_STATE_OBSERVATION_DISTRIBUTION_DISTANCE] << endl;

      for (j = 0;j < nb_state;j++) {
        for (k = 0;k < nb_state;k++) {
          if ((k != j) && (logic_transition[j][k])) {
            os << distance[j][k];
          }
          os << "\t";
        }
        os << endl;
      }
    }
  }

  if (hidden) {
    for (i = 0;i < nb_state;i++) {
      delete [] logic_transition[i];
    }
    delete [] logic_transition;

   for (i = 0;i < nb_state;i++) {
     delete [] distance[i];
    }
    delete [] distance;
  }

  if (seq) {
    int nb_parameter = nb_parameter_computation(hidden ? MIN_PROBABILITY : 0.) , nb_transient_parameter;
    double information;


    // writing of the sequence length frequency distribution

    os << "\n" << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\t";
    seq->length_distribution->spreadsheet_characteristic_print(os);

    os << "\n\t" << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << endl;
    seq->length_distribution->spreadsheet_print(os);

    os << "\n" << SEQ_label[SEQL_CUMUL_LENGTH] << "\t" << seq->cumul_length << endl;

    // writing of the information quantity of the observed sequences in the i.i.d. case

    for (i = 0;i < seq->nb_variable;i++) {
      if (seq->type[i] == REAL_VALUE) {
        break;
      }
    }

    if (i == seq->nb_variable) {
      information = seq->iid_information_computation();

      os << "\n" << SEQ_label[SEQL_IID_INFORMATION] << "\t" << information << "\t"
         << information / seq->cumul_length << endl;
    }

    // writing of the (penalized) log-likelihoods of the model for sequences

    if (hidden) {
      if (seq->restoration_likelihood != D_INF) {
        os << "\n" << SEQ_label[SEQL_STATE_SEQUENCES_LIKELIHOOD] << "\t" << seq->restoration_likelihood << "\t"
           << STAT_label[STATL_NORMALIZED] << "\t" << seq->restoration_likelihood / seq->cumul_length << endl;
      }

      if (seq->sample_entropy != D_DEFAULT) {
        os << "\n" << SEQ_label[SEQL_STATE_SEQUENCE_ENTROPY] << "\t" << seq->sample_entropy << "\t"
           << STAT_label[STATL_NORMALIZED] << "\t" << seq->sample_entropy / seq->cumul_length << endl;
      }

      if (seq->likelihood != D_INF) {
        os << "\n" << SEQ_label[SEQL_OBSERVED_SEQUENCES_LIKELIHOOD] << "\t" << seq->likelihood << "\t"
           << STAT_label[STATL_NORMALIZED] << "\t" << seq->likelihood / seq->cumul_length << endl;
      }
    }

    else {
      os << "\n" << STAT_label[STATL_LIKELIHOOD] << "\t" << seq->likelihood << "\t"
         << STAT_label[STATL_NORMALIZED] << "\t" << seq->likelihood / seq->cumul_length << endl;
    }

    if (seq->likelihood != D_INF) {
      if (type == ORDINARY) {
        nb_transient_parameter = nb_transient_parameter_computation(hidden ? MIN_PROBABILITY : 0.);

        os << "\n" << nb_transient_parameter << "\t"
           << SEQ_label[nb_transient_parameter == 1 ? SEQL_FREE_TRANSIENT_PARAMETER : SEQL_FREE_TRANSIENT_PARAMETERS] << endl;
      }

      os << "\n" << nb_parameter << "\t" << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS] << "\t"
         << "2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[AIC] << ")\t"
         << 2 * (seq->likelihood - nb_parameter) << endl;
      if ((type == ORDINARY) && (nb_transient_parameter > 0)) {
        os << nb_transient_parameter + nb_parameter << "\t" << STAT_label[STATL_FREE_PARAMETERS] << "\t"
           << "2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[AIC] << ")\t"
           << 2 * (seq->likelihood - nb_transient_parameter - nb_parameter) << endl;
      }

      if (nb_parameter < seq->cumul_length - 1) {
        os << "\n" << nb_parameter << "\t" << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS] << "\t"
           << "2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[AICc] << ")\t"
           << 2 * (seq->likelihood - (double)(nb_parameter * seq->cumul_length) /
               (double)(seq->cumul_length - nb_parameter - 1)) << endl;
      }
      if ((type == ORDINARY) && (nb_transient_parameter > 0) &&
          (nb_transient_parameter + nb_parameter < seq->cumul_length - 1)) {
        os << nb_transient_parameter + nb_parameter << "\t" << STAT_label[STATL_FREE_PARAMETERS] << "\t"
           << "2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[AICc] << ")\t"
           << 2 * (seq->likelihood - (double)((nb_transient_parameter + nb_parameter) * seq->cumul_length) /
               (double)(seq->cumul_length - nb_transient_parameter - nb_parameter - 1)) << endl;
      }

      os << "\n" << nb_parameter << "\t" << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS] << "\t"
         << "2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[BIC] << ")\t"
         << 2 * seq->likelihood - nb_parameter * log((double)seq->cumul_length) << endl;
      if ((type == ORDINARY) && (nb_transient_parameter > 0)) {
        os << nb_transient_parameter + nb_parameter << "\t" << STAT_label[STATL_FREE_PARAMETERS] << "\t"
           << "2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[BIC] << ")\t"
           << 2 * seq->likelihood - (nb_transient_parameter + nb_parameter) * log((double)seq->cumul_length) << endl;
      }

      os << "\n" << nb_parameter + (type == ORDINARY ? nb_transient_parameter : 0) << "\t" << STAT_label[STATL_FREE_PARAMETERS] << "\t"
         << "2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[BICc] << ")\t"
         << 2 * seq->likelihood - penalty_computation(hidden , (hidden ? MIN_PROBABILITY : 0.)) << endl;
    }

    if ((hidden) && (seq->likelihood != D_INF)) {
      os << "\n" << nb_parameter << "\t" << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS] << "\t"
         << "2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[ICL] << ")\t"
         << 2 * (seq->likelihood - seq->sample_entropy) - nb_parameter * log((double)seq->cumul_length) << endl;
      if ((type == ORDINARY) && (nb_transient_parameter > 0)) {
        os << nb_transient_parameter + nb_parameter << "\t" << STAT_label[STATL_FREE_PARAMETERS] << "\t"
           << "2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[ICL] << ")\t"
           << 2 * (seq->likelihood - seq->sample_entropy) - (nb_transient_parameter + nb_parameter) * log((double)seq->cumul_length) << endl;
      }

      os << "\n" << nb_parameter << "\t" << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS] << "\t"
         << "2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[ICLc] << ")\t"
         << 2 * (seq->likelihood - seq->sample_entropy) - penalty_computation(hidden , MIN_PROBABILITY) << endl;
    }
  }

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a VariableOrderMarkov object in a file at the spreadsheet format.
 *
 *  \param[in] error reference on a StatError object,
 *  \param[in] path  file path.
 *
 *  \return          error status.
 */
/*--------------------------------------------------------------*/

bool VariableOrderMarkov::spreadsheet_write(StatError &error , const string path) const

{
  bool status;
  ofstream out_file(path.c_str());


  error.init();

  if (!out_file) {
    status = false;
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {
    status = true;
    spreadsheet_write(out_file , markov_data);
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Plot of a VariableOrderMarkov object and the associated data structure
 *         using Gnuplot.
 *
 *  \param[in] prefix file prefix,
 *  \param[in] title  figure title,
 *  \param[in] seq    pointer on a VariableOrderMarkovData object.
 *
 *  \return           error status.
 */
/*--------------------------------------------------------------*/

bool VariableOrderMarkov::plot_write(const char *prefix , const char *title ,
                                     const VariableOrderMarkovData *seq) const

{
  bool status;
  int i;
  int variable , nb_value = I_DEFAULT;
  double *empirical_cdf[2];
  FrequencyDistribution *length_distribution = NULL , *marginal_dist = NULL , **observation_dist = NULL;
  Histogram *marginal_histo = NULL , **observation_histo = NULL;
  SequenceCharacteristics *characteristics = NULL;


  if ((seq) && (seq->type[0] == STATE)) {
    characteristics = seq->characteristics[0];
    length_distribution = seq->length_distribution;
  }
  else {
    characteristics = NULL;
  }

  status = state_process->plot_print(prefix , title , 0 , NULL , NULL ,
                                     characteristics , length_distribution);

  if (status) {
    if (seq) {
      length_distribution = seq->length_distribution;
    }

    for (i = 0;i < nb_output_process;i++) {
      if (seq) {
        switch (seq->type[0]) {
        case STATE :
          variable = i + 1;
          break;
        default :
          variable = i;
          break;
        }

        if (seq->observation_distribution) {
          observation_dist = seq->observation_distribution[variable];
        }
        marginal_dist = seq->marginal_distribution[variable];

        if (seq->observation_histogram) {
          observation_histo = seq->observation_histogram[variable];
        }
        marginal_histo = seq->marginal_histogram[variable];

        characteristics = seq->characteristics[variable];

        if (continuous_parametric_process[i]) {
          nb_value = seq->cumulative_distribution_function_computation(variable , empirical_cdf);
        }
      }

      if (categorical_process[i]) {
        categorical_process[i]->plot_print(prefix , title , i + 1 , observation_dist ,
                                           marginal_dist , characteristics ,
                                           length_distribution);
      }
      else if (discrete_parametric_process[i]) {
        discrete_parametric_process[i]->plot_print(prefix , title , i + 1 , observation_dist ,
                                                   marginal_dist);
      }
      else {
        continuous_parametric_process[i]->plot_print(prefix , title , i + 1 ,
                                                     observation_histo , observation_dist ,
                                                     marginal_histo , marginal_dist ,
                                                     nb_value , (seq ? empirical_cdf : NULL));
        if (seq) {
          delete [] empirical_cdf[0];
          delete [] empirical_cdf[1];
        }
      }
    }
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Plot of a VariableOrderMarkov object using Gnuplot.
 *
 *  \param[in] error  reference on a StatError object,
 *  \param[in] prefix file prefix,
 *  \param[in] title  figure title.
 *
 *  \return           error status.
 */
/*--------------------------------------------------------------*/

bool VariableOrderMarkov::plot_write(StatError &error , const char *prefix ,
                                     const char *title) const

{
  bool status = plot_write(prefix , title , markov_data);

  error.init();

  if (!status) {
    error.update(STAT_error[STATR_FILE_PREFIX]);
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Plot of a VariableOrderMarkov object and the associated data structure.
 *
 *  \param[in] seq pointer on a VariableOrderMarkovData object.
 *
 *  \return        MultiPlotSet object.
 */
/*--------------------------------------------------------------*/

MultiPlotSet* VariableOrderMarkov::get_plotable(const VariableOrderMarkovData *seq) const

{
  int i , j;
  int nb_plot_set , index_length , index , variable;
  FrequencyDistribution *length_distribution = NULL , *marginal_dist = NULL , **observation_dist = NULL;
  Histogram *marginal_histo = NULL , **observation_histo = NULL;
  SequenceCharacteristics *characteristics = NULL;
  MultiPlotSet *plot_set;


  if ((seq) && (seq->type[0] == STATE)) {
    characteristics = seq->characteristics[0];
  }
  else {
    characteristics = NULL;
  }

  // computation of the number of plots

  nb_plot_set = 0;

  if ((state_process->index_value) || (characteristics)) {
    nb_plot_set++;

    if (characteristics) {
      index_length = characteristics->index_value->plot_length_computation();

      if (characteristics->index_value->frequency[index_length - 1] < MAX_FREQUENCY) {
        nb_plot_set++;
      }
      nb_plot_set++;
    }
  }

  if ((state_process->first_occurrence) || (characteristics)) {
    for (i = 0;i < nb_state;i++) {
      if ((state_process->first_occurrence) &&
          (state_process->first_occurrence[i])) {
        nb_plot_set++;
      }
      else if ((characteristics) && (i < characteristics->nb_value) &&
               (characteristics->first_occurrence[i]->nb_element > 0)) {
        nb_plot_set++;
      }
    }
  }

  if ((state_process->recurrence_time) || (characteristics)) {
    for (i = 0;i < nb_state;i++) {
      if ((state_process->recurrence_time) &&
          (state_process->recurrence_time[i])) {
        nb_plot_set++;
      }
      else if ((characteristics) && (i < characteristics->nb_value) &&
               (characteristics->recurrence_time[i]->nb_element > 0)) {
        nb_plot_set++;
      }
    }
  }

  if ((state_process->sojourn_time) || (characteristics)) {
    for (i = 0;i < nb_state;i++) {
      if ((state_process->sojourn_time) &&
          (state_process->sojourn_time[i])) {
        nb_plot_set++;
      }
      else if ((characteristics) && (i < characteristics->nb_value) &&
               (characteristics->sojourn_time[i]->nb_element > 0)) {
        nb_plot_set++;
      }

      if ((characteristics) && (i < characteristics->nb_value) &&
          (characteristics->initial_run) &&
          (characteristics->initial_run[i]->nb_element > 0)) {
        nb_plot_set++;
      }

      if ((characteristics) && (i < characteristics->nb_value) &&
          (characteristics->final_run[i]->nb_element > 0)) {
        nb_plot_set++;
      }
    }
  }

  if ((state_process->nb_run) || (state_process->nb_occurrence) ||
      ((characteristics) && (characteristics->nb_run) && (characteristics->nb_occurrence))) {
    for (i = 0;i < nb_state;i++) {
      if (state_process->nb_run) {
        nb_plot_set++;
      }
      else if ((characteristics) && (i < characteristics->nb_value) &&
               (characteristics->nb_run) && (characteristics->nb_run[i]->nb_element > 0)) {
        nb_plot_set++;
      }

      if (state_process->nb_occurrence) {
        nb_plot_set++;
      }
      else if ((characteristics) && (i < characteristics->nb_value) &&
               (characteristics->nb_occurrence) &&
               (characteristics->nb_occurrence[i]->nb_element > 0)) {
        nb_plot_set++;
      }
    }

    if ((characteristics) && (characteristics->nb_run) && (characteristics->nb_occurrence)) {
      nb_plot_set++;
    }
  }

  for (i = 0;i < nb_output_process;i++) {
    if (seq) {
      switch (seq->type[0]) {
      case STATE :
        variable = i + 1;
        break;
      default :
        variable = i;
        break;
      }

      characteristics = seq->characteristics[variable];
    }

    if (categorical_process[i]) {
      if ((categorical_process[i]->index_value) || (characteristics)) {
        nb_plot_set++;

        if (characteristics) {
          index_length = characteristics->index_value->plot_length_computation();

          if (characteristics->index_value->frequency[index_length - 1] < MAX_FREQUENCY) {
            nb_plot_set++;
          }
          nb_plot_set++;
        }
      }

      if ((categorical_process[i]->first_occurrence) || (characteristics)) {
        for (j = 0;j < categorical_process[i]->nb_value;j++) {
          if ((categorical_process[i]->first_occurrence) &&
              (categorical_process[i]->first_occurrence[j])) {
            nb_plot_set++;
          }
          else if ((characteristics) && (j < characteristics->nb_value) &&
                   (characteristics->first_occurrence[j]->nb_element > 0)) {
            nb_plot_set++;
          }
        }
      }

      if ((categorical_process[i]->recurrence_time) || (characteristics)) {
        for (j = 0;j < categorical_process[i]->nb_value;j++) {
          if ((categorical_process[i]->recurrence_time) &&
              (categorical_process[i]->recurrence_time[j])) {
            nb_plot_set++;
          }
          else if ((characteristics) && (i < characteristics->nb_value) &&
                   (characteristics->recurrence_time[j]->nb_element > 0)) {
            nb_plot_set++;
          }
        }
      }

      if ((categorical_process[i]->sojourn_time) || (characteristics)) {
        for (j = 0;j < categorical_process[i]->nb_value;j++) {
          if ((categorical_process[i]->sojourn_time) &&
              (categorical_process[i]->sojourn_time[j])) {
            nb_plot_set++;
          }
          else if ((characteristics) && (i < characteristics->nb_value) &&
                   (characteristics->sojourn_time[j]->nb_element > 0)) {
            nb_plot_set++;
          }

/*          if ((characteristics) && (j < characteristics->nb_value) &&
              (characteristics->initial_run) &&
              (characteristics->initial_run[j]->nb_element > 0)) {
            nb_plot_set++;
          } */

          if ((characteristics) && (j < characteristics->nb_value) &&
              (characteristics->final_run[j]->nb_element > 0)) {
            nb_plot_set++;
          }
        }
      }

      if ((categorical_process[i]->nb_run) || (categorical_process[i]->nb_occurrence) ||
          ((characteristics) && (characteristics->nb_run) && (characteristics->nb_occurrence))) {
        for (j = 0;j < categorical_process[i]->nb_value;j++) {
          if (categorical_process[i]->nb_run) {
            nb_plot_set++;
          }
          else if ((characteristics) && (j < characteristics->nb_value) &&
                   (characteristics->nb_run) && (characteristics->nb_run[j]->nb_element > 0)) {
            nb_plot_set++;
          }

          if (categorical_process[i]->nb_occurrence) {
            nb_plot_set++;
          }
          else if ((characteristics) && (j < characteristics->nb_value) &&
                   (characteristics->nb_occurrence) &&
                   (characteristics->nb_occurrence[j]->nb_element > 0)) {
            nb_plot_set++;
          }
        }

        if ((characteristics) && (characteristics->nb_run) && (characteristics->nb_occurrence)) {
          nb_plot_set++;
        }
      }
    }

    if ((seq->observation_distribution) || (seq->observation_histogram)) {
      nb_plot_set += nb_state;
    }
    else {
      nb_plot_set++;
    }

    if ((categorical_process[i]) && (seq->marginal_distribution[variable])) {
      if ((categorical_process[i]->weight) &&
          (categorical_process[i]->mixture)) {
        nb_plot_set++;
      }
      if ((categorical_process[i]->restoration_weight) &&
          (categorical_process[i]->restoration_mixture)) {
        nb_plot_set++;
      }
    }

    if ((discrete_parametric_process[i]) && (seq->marginal_distribution[variable])) {
      if ((discrete_parametric_process[i]->weight) &&
          (discrete_parametric_process[i]->mixture)) {
        nb_plot_set += 2;
      }
      if ((discrete_parametric_process[i]->restoration_weight) &&
          (discrete_parametric_process[i]->restoration_mixture)) {
        nb_plot_set += 2;
      }
    }

    if ((continuous_parametric_process[i]) && ((seq->marginal_histogram[variable]) ||
         (seq->marginal_distribution[variable]))) {
      if (continuous_parametric_process[i]->weight) {
        nb_plot_set += 2;
      }
      if (continuous_parametric_process[i]->restoration_weight) {
        nb_plot_set += 2;
      }
    }
  }

  plot_set = new MultiPlotSet(nb_plot_set , nb_output_process + 1);
  plot_set->border = "15 lw 0";

  if ((seq) && (seq->type[0] == STATE)) {
    characteristics = seq->characteristics[0];
    length_distribution = seq->length_distribution;
  }
  else {
    characteristics = NULL;
  }

  index = 0;
  plot_set->variable_nb_viewpoint[0] = 0;
  state_process->plotable_write(*plot_set , index , 0 , NULL , NULL , characteristics ,
                                length_distribution);

  if (seq) {
    length_distribution = seq->length_distribution;
  }

  for (i = 0;i < nb_output_process;i++) {
    if (seq) {
      switch (seq->type[0]) {
      case STATE :
        variable = i + 1;
        break;
      default :
        variable = i;
        break;
      }

      if (seq->observation_distribution) {
        observation_dist = seq->observation_distribution[variable];
      }
      marginal_dist = seq->marginal_distribution[variable];

      if (seq->observation_histogram) {
        observation_histo = seq->observation_histogram[variable];
      }
      marginal_histo = seq->marginal_histogram[variable];

      characteristics = seq->characteristics[variable];
    }

    if (categorical_process[i]) {
      plot_set->variable_nb_viewpoint[i] = 0;
      categorical_process[i]->plotable_write(*plot_set , index , i + 1 , observation_dist ,
                                             marginal_dist , characteristics ,
                                             length_distribution);
    }
    else if (discrete_parametric_process[i]){
      discrete_parametric_process[i]->plotable_write(*plot_set , index , i + 1 , observation_dist ,
                                                     marginal_dist);
    }
    else {
      continuous_parametric_process[i]->plotable_write(*plot_set , index , i + 1 ,
                                                       observation_histo , observation_dist ,
                                                       marginal_histo , marginal_dist);
    }
  }

  return plot_set;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Plot of a VariableOrderMarkov object.
 *
 *  \return MultiPlotSet object.
 */
/*--------------------------------------------------------------*/

MultiPlotSet* VariableOrderMarkov::get_plotable() const

{
  return get_plotable(markov_data);
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the number of parameters of a VariableOrderMarkov object.
 *
 *  \param[in] min_probability minimum probability.
 *
 *  \return                    number of parameters.
 */
/*--------------------------------------------------------------*/

int VariableOrderMarkov::nb_parameter_computation(double min_probability) const

{
  int i;
  int nb_parameter = VariableOrderMarkovChain::nb_parameter_computation(min_probability);


  for (i = 0;i < nb_output_process;i++) {
    if (categorical_process[i]) {
      nb_parameter += categorical_process[i]->nb_parameter_computation(min_probability);
    }
    else if (discrete_parametric_process[i]) {
      nb_parameter += discrete_parametric_process[i]->nb_parameter_computation();
    }
    else {
      nb_parameter += continuous_parametric_process[i]->nb_parameter_computation();
    }
  }

  return nb_parameter;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of an adaptative penalty.
 *
 *  \param[in] hidden          flag hidden model,
 *  \param[in] min_probability minimum probability.
 *
 *  \return                    adaptative penalty.
 */
/*--------------------------------------------------------------*/

double VariableOrderMarkov::penalty_computation(bool hidden , double min_probability) const

{
  int i , j , k;
  int nb_parameter , sample_size;
  double sum , *state_marginal , *memory;
  double penalty = 0.;


  if (markov_data) {
    if (hidden) {
      switch (type) {

      case ORDINARY : {
        memory = memory_computation();

        sum = 0.;
        for (i = 1;i < nb_row;i++) {
//          if (memo_type[i] == TERMINAL) {
            sum += memory[i];
//          }
        }
        for (i = 1;i < nb_row;i++) {
//          if (memo_type[i] == TERMINAL) {
            memory[i] /= sum;
//          }
        }
        break;
      }

      case EQUILIBRIUM : {
        memory = new double[nb_row];
        for (i = 1;i < nb_row;i++) {
          memory[i] = initial[i];
        }
        break;
      }
      }

      state_marginal = new double[nb_state];
      for (i = 0;i < nb_state;i++) {
        state_marginal[i] = 0.;
      }
      for (i = 1;i < nb_row;i++) {
//        if (memo_type[i] == TERMINAL) {
        state_marginal[state[i][0]] += memory[i];
//        }
      }
    }

    for (i = 1;i < nb_row;i++) {
//      if (memo_type[i] == TERMINAL) {
      if ((memo_type[i] == TERMINAL) || ((type == ORDINARY) &&
           (memo_type[i] == NON_TERMINAL))) {
        nb_parameter = 0;
        if (!hidden) {
          sample_size = 0;
        }
        for (j = 0;j < nb_state;j++) {
          if (transition[i][j] > min_probability) {
            nb_parameter++;
            if (!hidden) {
              sample_size += markov_data->chain_data->transition[i][j];
            }
          }
        }

        nb_parameter--;

        if (nb_parameter > 0) {
          if (hidden) {
            if (memory[i] > 0.) {
              penalty += nb_parameter * log(memory[i] * markov_data->cumul_length);
            }
          }
          else {
            if (sample_size > 0) {
              penalty += nb_parameter * log((double)sample_size);
            }
          }
        }
      }
    }

    for (i = 0;i < nb_output_process;i++) {
      if (categorical_process[i]) {
        for (j = 0;j < nb_state;j++) {
          nb_parameter = 0;
          for (k = 0;k < categorical_process[i]->nb_value;k++) {
            if (categorical_process[i]->observation[j]->mass[k] > min_probability) {
              nb_parameter++;
            }
          }

          nb_parameter--;

          if (nb_parameter > 0) {
            if (hidden) {
              penalty += nb_parameter * log(state_marginal[j] * markov_data->cumul_length);
            }
            else {
              penalty += nb_parameter *
                         log((double)markov_data->marginal_distribution[0]->frequency[j]);
            }
          }
        }
      }

      else if (discrete_parametric_process[i]) {
        for (j = 0;j < nb_state;j++) {
          nb_parameter = discrete_parametric_process[i]->observation[j]->nb_parameter_computation();

          if (hidden) {
            penalty += nb_parameter * log(state_marginal[j] * markov_data->cumul_length);
          }
          else {
            penalty += nb_parameter *
                       log((double)markov_data->marginal_distribution[0]->frequency[j]);
          }
        }
      }

      else {
        for (j = 0;j < nb_state;j++) {
          nb_parameter = continuous_parametric_process[i]->observation[j]->nb_parameter_computation();

          if (hidden) {
            penalty += nb_parameter * log(state_marginal[j] * markov_data->cumul_length);
          }
          else {
            penalty += nb_parameter *
                       log((double)markov_data->marginal_distribution[0]->frequency[j]);
          }
        }
      }
    }

    if (hidden) {
      delete [] memory;
      delete [] state_marginal;
    }
  }

  return penalty;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Default constructor of the VariableOrderMarkovData class.
 */
/*--------------------------------------------------------------*/

VariableOrderMarkovData::VariableOrderMarkovData()

{
  markov = NULL;
  chain_data = NULL;

  likelihood = D_INF;
  restoration_likelihood = D_INF;
  sample_entropy = D_DEFAULT;

  posterior_probability = NULL;
  entropy = NULL;
  nb_state_sequence = NULL;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Constructor of the VariableOrderMarkovData class.
 *
 *  \param[in] ilength_distribution sequence length frequency distribution,
 *  \param[in] inb_variable         number of variables,
 *  \param[in] itype                variable types,
 *  \param[in] init_flag            flag initialization.
 */
/*--------------------------------------------------------------*/

VariableOrderMarkovData::VariableOrderMarkovData(const FrequencyDistribution &ilength_distribution ,
                                                 int inb_variable , variable_nature *itype , bool init_flag)
:MarkovianSequences(ilength_distribution , inb_variable , itype , init_flag)

{
  markov = NULL;
  chain_data = NULL;

  likelihood = D_INF;
  restoration_likelihood = D_INF;
  sample_entropy = D_DEFAULT;

  posterior_probability = NULL;
  entropy = NULL;
  nb_state_sequence = NULL;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Construction of a VariableOrderMarkovData object from
 *         a MarkovianSequences object adding a state variable.
 *
 *  \param[in] seq reference on a MarkovianSequences object.
 */
/*--------------------------------------------------------------*/

VariableOrderMarkovData::VariableOrderMarkovData(const MarkovianSequences &seq)
:MarkovianSequences(seq , ADD_STATE_VARIABLE , UNCHANGED)

{
  markov = NULL;
  chain_data = NULL;

  likelihood = D_INF;
  restoration_likelihood = D_INF;
  sample_entropy = D_DEFAULT;

  posterior_probability = NULL;
  entropy = NULL;
  nb_state_sequence = NULL;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Construction of a VariableOrderMarkovData object from
 *         a MarkovianSequences object.
 *
 *  \param[in] seq              reference on a MarkovianSequences object,
 *  \param[in] transform        type of transform (SEQUENCE_COPY/ADD_STATE_VARIABLE),
 *  \param[in] initial_run_flag addition/removing of the initial run length frequency distributions.
 */
/*--------------------------------------------------------------*/

VariableOrderMarkovData::VariableOrderMarkovData(const MarkovianSequences &seq ,
                                                 sequence_transformation transform , bool initial_run_flag)
:MarkovianSequences(seq , transform , (initial_run_flag ? ADD_INITIAL_RUN : REMOVE_INITIAL_RUN))

{
  markov = NULL;
  chain_data = NULL;

  likelihood = D_INF;
  restoration_likelihood = D_INF;
  sample_entropy = D_DEFAULT;

  posterior_probability = NULL;
  entropy = NULL;
  nb_state_sequence = NULL;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Copy of a VariableOrderMarkovData object.
 *
 *  \param[in] seq        reference on a VariableOrderMarkovData object,
 *  \param[in] model_flag flag copy of the included VariableOrderMarkov object.
 */
/*--------------------------------------------------------------*/

void VariableOrderMarkovData::copy(const VariableOrderMarkovData &seq ,
                                   bool model_flag)

{
  int i;


  if ((model_flag) && (seq.markov)) {
    markov = new VariableOrderMarkov(*(seq.markov) , false);
  }
  else {
    markov = NULL;
  }

  if (seq.chain_data) {
    chain_data = new VariableOrderMarkovChainData(*(seq.chain_data));
  }
  else {
    chain_data = NULL;
  }

  likelihood = seq.likelihood;
  restoration_likelihood = seq.restoration_likelihood;
  sample_entropy = seq.sample_entropy;

  if (seq.posterior_probability) {
    posterior_probability = new double[nb_sequence];
    for (i = 0;i < nb_sequence;i++) {
      posterior_probability[i] = seq.posterior_probability[i];
    }
  }
  else {
    posterior_probability = NULL;
  }

  if (seq.entropy) {
    entropy = new double[nb_sequence];
    for (i = 0;i < nb_sequence;i++) {
      entropy[i] = seq.entropy[i];
    }
  }
  else {
    entropy = NULL;
  }

  if (seq.nb_state_sequence) {
    nb_state_sequence = new double[nb_sequence];
    for (i = 0;i < nb_sequence;i++) {
      nb_state_sequence[i] = seq.nb_state_sequence[i];
    }
  }
  else {
    nb_state_sequence = NULL;
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Destructor of the VariableOrderMarkovData class.
 */
/*--------------------------------------------------------------*/

VariableOrderMarkovData::~VariableOrderMarkovData()

{
  delete markov;
  delete chain_data;

  delete [] posterior_probability;
  delete [] entropy;
  delete [] nb_state_sequence;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Assignment operator of the VariableOrderMarkovData class.
 *
 *  \param[in] seq reference on a VariableOrderMarkovData object.
 *
 *  \return        VariableOrderMarkovData object.
 */
/*--------------------------------------------------------------*/

VariableOrderMarkovData& VariableOrderMarkovData::operator=(const VariableOrderMarkovData &seq)

{
  if (&seq != this) {
    delete markov;
    delete chain_data;

    delete [] posterior_probability;
    delete [] entropy;
    delete [] nb_state_sequence;

    remove();
    Sequences::remove();

    Sequences::copy(seq);
    MarkovianSequences::copy(seq);
    copy(seq);
  }

  return *this;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Extraction of a frequency distribution.
 *
 *  \param[in] error      reference on a StatError object,
 *  \param[in] histo_type frequency distribution type,
 *  \param[in] variable   variable index,
 *  \param[in] value      state or observation.
 *
 *  \return               DiscreteDistributionData object.
 */
/*--------------------------------------------------------------*/

DiscreteDistributionData* VariableOrderMarkovData::extract(StatError &error , process_distribution histo_type ,
                                                           int variable , int value) const

{
  bool status = true;
  Distribution *pdist;
  DiscreteParametric *pparam;
  FrequencyDistribution *phisto;
  DiscreteDistributionData *histo;
  CategoricalSequenceProcess *process;


  histo = NULL;
  error.init();

  if (histo_type == OBSERVATION) {
    if ((variable < 2) || (variable > nb_variable)) {
      status = false;
      error.update(STAT_error[STATR_VARIABLE_INDEX]);
    }

    else {
      variable--;

      if ((value < 0) || (value >= marginal_distribution[0]->nb_value)) {
        status = false;
        ostringstream error_message;
        error_message << STAT_label[STATL_STATE] << " " << value << " "
                      << STAT_error[STATR_NOT_PRESENT];
        error.update((error_message.str()).c_str());
      }

      else {
        if (!observation_distribution[variable]) {
          status = false;
          ostringstream error_message;
          error_message << STAT_label[STATL_VARIABLE] << " " << variable + 1 << ": "
                        << STAT_error[STATR_VARIABLE_TYPE];
          error.correction_update((error_message.str()).c_str() , STAT_variable_word[INT_VALUE]);
        }

        else {
          phisto = observation_distribution[variable][value];

          if (phisto->nb_element == 0) {
            status = false;
            error.update(STAT_error[STATR_EMPTY_SAMPLE]);
          }
        }
      }
    }
  }

  else {
    if ((variable < 1) || (variable > nb_variable)) {
      status = false;
      error.update(STAT_error[STATR_VARIABLE_INDEX]);
    }

    else {
      variable--;

      if (!characteristics[variable]) {
        status = false;
        ostringstream error_message;
        error_message << STAT_label[STATL_VARIABLE] << " " << variable + 1 << ": "
                      << SEQ_error[SEQR_CHARACTERISTICS_NOT_COMPUTED];
        error.update((error_message.str()).c_str());
      }

      else if ((value < 0) || (value >= marginal_distribution[variable]->nb_value)) {
        status = false;
        ostringstream error_message;
        error_message << STAT_label[variable == 0 ? STATL_STATE : STATL_OUTPUT] << " "
                      << value << " " << STAT_error[STATR_NOT_PRESENT];
        error.update((error_message.str()).c_str());
      }

      else {
        switch (histo_type) {
        case FIRST_OCCURRENCE :
          phisto = characteristics[variable]->first_occurrence[value];
          break;
        case RECURRENCE_TIME :
          phisto = characteristics[variable]->recurrence_time[value];
          break;
        case SOJOURN_TIME :
          phisto = characteristics[variable]->sojourn_time[value];
          break;
        case FINAL_RUN :
          phisto = characteristics[variable]->final_run[value];
          break;
        case NB_RUN :
          phisto = characteristics[variable]->nb_run[value];
          break;
        case NB_OCCURRENCE :
          phisto = characteristics[variable]->nb_occurrence[value];
          break;
        }

        if (phisto->nb_element == 0) {
          status = false;
          error.update(STAT_error[STATR_EMPTY_SAMPLE]);
        }
      }
    }
  }

  if (status) {
    if (variable == 0) {
      process = markov->state_process;
    }
    else {
      process = markov->categorical_process[variable - 1];
    }

    pdist = NULL;
    pparam = NULL;

    switch (histo_type) {

    case OBSERVATION : {
      if (markov->categorical_process[variable - 1]) {
        pdist = markov->categorical_process[variable - 1]->observation[value];
      }
      else if (markov->discrete_parametric_process[variable - 1]) {
        pparam = markov->discrete_parametric_process[variable - 1]->observation[value];
      }
      break;
    }

    case FIRST_OCCURRENCE : {
      pdist = process->first_occurrence[value];
      break;
    }

    case RECURRENCE_TIME : {
      pdist = process->recurrence_time[value];
      break;
    }

    case SOJOURN_TIME : {
      pparam = process->sojourn_time[value];
      break;
    }

    case NB_RUN : {
      pdist = process->nb_run[value];
      break;
    }

    case NB_OCCURRENCE : {
      pdist = process->nb_occurrence[value];
      break;
    }
    }

    if (pdist) {
      histo = new DiscreteDistributionData(*phisto , pdist);
    }
    else {
      histo = new DiscreteDistributionData(*phisto , pparam);
    }
  }

  return histo;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Copy of a VariableOrderMarkovData object transforming
 *         the implicit index parameters in explicit index parameters.
 *
 *  \param[in] error reference on a StatError object.
 *
 *  \return          VariableOrderMarkovData object.
 */
/*--------------------------------------------------------------*/

VariableOrderMarkovData* VariableOrderMarkovData::explicit_index_parameter(StatError &error) const

{
  VariableOrderMarkovData *seq;


  error.init();

  if (index_parameter) {
    seq = NULL;
    error.update(SEQ_error[SEQR_INDEX_PARAMETER_TYPE]);
  }
  else {
    seq = new VariableOrderMarkovData(*this , true , EXPLICIT_INDEX_PARAMETER);
  }

  return seq;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Removing of the index parameters.
 *
 *  \param[in] error reference on a StatError object.
 *
 *  \return          VariableOrderMarkovData object.
 */
/*--------------------------------------------------------------*/

VariableOrderMarkovData* VariableOrderMarkovData::remove_index_parameter(StatError &error) const

{
  VariableOrderMarkovData *seq;


  error.init();

  if (!index_parameter) {
    seq = NULL;
    error.update(SEQ_error[SEQR_INDEX_PARAMETER_TYPE]);
  }
  else {
    seq = new VariableOrderMarkovData(*this , true , REMOVE_INDEX_PARAMETER);
  }

  return seq;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Construction of the auxiliary variables corresponding to
 *         the restored state sequences.
 *
 *  \param[in] error reference on a StatError object.
 *
 *  \return          VariableOrderMarkovData object.
 */
/*--------------------------------------------------------------*/

MarkovianSequences* VariableOrderMarkovData::build_auxiliary_variable(StatError &error) const

{
  bool status = true;
  int i;
  MarkovianSequences *seq;


  seq = NULL;
  error.init();

  if (type[0] != STATE) {
    status = false;
    ostringstream error_message;
    error_message << STAT_label[STATL_VARIABLE] << " 1: "
                  << STAT_error[STATR_VARIABLE_TYPE];
    error.correction_update((error_message.str()).c_str() , STAT_variable_word[STATE]);
  }

  for (i = 0;i < markov->nb_output_process;i++) {
    if (((markov->discrete_parametric_process) && (markov->discrete_parametric_process[i])) ||
        ((markov->continuous_parametric_process) && (markov->continuous_parametric_process[i]))) {
      break;
    }
  }

  if (i == markov->nb_output_process) {
    status = false;
    error.update(SEQ_error[SEQR_PARAMETRIC_PROCESS]);
  }

  if (status) {
    seq = MarkovianSequences::build_auxiliary_variable(markov->discrete_parametric_process ,
                                                       markov->continuous_parametric_process);
  }

  return seq;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Building of residual sequences on the basis of restored state sequences.
 *
 *  \param[in] error reference on a StatError object.
 *
 *  \return          MarkovianSequences object.
 */
/*--------------------------------------------------------------*/

MarkovianSequences* VariableOrderMarkovData::residual_sequences(StatError &error) const

{
  MarkovianSequences *seq;


  error.init();

  if (type[0] != STATE) {
    seq = NULL;

    ostringstream error_message;
    error_message << STAT_label[STATL_VARIABLE] << " 1: "
                  << STAT_error[STATR_VARIABLE_TYPE];
    error.correction_update((error_message.str()).c_str() , STAT_variable_word[STATE]);
  }

  else {
    seq = MarkovianSequences::residual_sequences(markov->categorical_process ,
                                                 markov->discrete_parametric_process ,
                                                 markov->continuous_parametric_process);
  }

  return seq;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a VariableOrderMarkovData object.
 *
 *  \param[in,out] os         stream,
 *  \param[in]     exhaustive flag detail level.
 */
/*--------------------------------------------------------------*/

ostream& VariableOrderMarkovData::ascii_write(ostream &os , bool exhaustive) const

{
  if (markov) {
    markov->ascii_write(os , this , exhaustive , false ,
                        CategoricalSequenceProcess::test_hidden(markov->nb_output_process , markov->categorical_process));
  }

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a VariableOrderMarkovData object in a file.
 *
 *  \param[in] error      reference on a StatError object,
 *  \param[in] path       file path,
 *  \param[in] exhaustive flag detail level.
 *
 *  \return               error status.
 */
/*--------------------------------------------------------------*/

bool VariableOrderMarkovData::ascii_write(StatError &error , const string path ,
                                          bool exhaustive) const

{
  bool status = false;


  if (markov) {
    ofstream out_file(path.c_str());

    error.init();

    if (!out_file) {
      status = false;
      error.update(STAT_error[STATR_FILE_NAME]);
    }

    else {
      status = true;
      markov->ascii_write(out_file , this , exhaustive , true ,
                          CategoricalSequenceProcess::test_hidden(markov->nb_output_process , markov->categorical_process));
    }
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a VariableOrderMarkovData object.
 *
 *  \param[in,out] os         stream,
 *  \param[in]     format     format (line/column),
 *  \param[in]     exhaustive flag detail level.
 */
/*--------------------------------------------------------------*/

ostream& VariableOrderMarkovData::ascii_data_write(ostream &os , output_sequence_format format ,
                                                   bool exhaustive) const

{
  MarkovianSequences::ascii_write(os , exhaustive , false);
  ascii_print(os , format , false , posterior_probability , entropy , nb_state_sequence);

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a VariableOrderMarkovData object in a file.
 *
 *  \param[in] error      reference on a StatError object,
 *  \param[in] path       file path,
 *  \param[in] format     format (line/column),
 *  \param[in] exhaustive flag detail level.
 *
 *  \return               error status.
 */
/*--------------------------------------------------------------*/

bool VariableOrderMarkovData::ascii_data_write(StatError &error , const string path ,
                                               output_sequence_format format , bool exhaustive) const

{
  bool status = false;
  ofstream out_file(path.c_str());


  error.init();

  if (!out_file) {
    status = false;
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {
    status = true;
    if (format != 'a') {
      MarkovianSequences::ascii_write(out_file , exhaustive , true);
    }
    ascii_print(out_file , format , true , posterior_probability , entropy , nb_state_sequence);
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a VariableOrderMarkovData object in a file at the spreadsheet format.
 *
 *  \param[in] error reference on a StatError object,
 *  \param[in] path  file path.
 *
 *  \return          error status.
 */
/*--------------------------------------------------------------*/

bool VariableOrderMarkovData::spreadsheet_write(StatError &error , const string path) const

{
  bool status = false;


  if (markov) {
    ofstream out_file(path.c_str());

    error.init();

    if (!out_file) {
      status = false;
      error.update(STAT_error[STATR_FILE_NAME]);
    }

    else {
      status = true;
      markov->spreadsheet_write(out_file , this ,
                                CategoricalSequenceProcess::test_hidden(markov->nb_output_process , markov->categorical_process));
    }
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Plot of a VariableOrderMarkovData object using Gnuplot.
 *
 *  \param[in] error  reference on a StatError object,
 *  \param[in] prefix file prefix,
 *  \param[in] title  figure title.
 *
 *  \return           error status.
 */
/*--------------------------------------------------------------*/

bool VariableOrderMarkovData::plot_write(StatError &error , const char *prefix ,
                                         const char *title) const

{
  bool status = false;


  if (markov) {
    status = markov->plot_write(prefix , title , this);

    error.init();

    if (!status) {
      error.update(STAT_error[STATR_FILE_PREFIX]);
    }
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Plot of a VariableOrderMarkovData object.
 *
 *  \return MultiPlotSet object.
 */
/*--------------------------------------------------------------*/

MultiPlotSet* VariableOrderMarkovData::get_plotable() const

{
  MultiPlotSet *plot_set;


  if (markov) {
    plot_set = markov->get_plotable(this);
  }
  else {
    plot_set = NULL;
  }

  return plot_set;
}


};  // namespace sequence_analysis
