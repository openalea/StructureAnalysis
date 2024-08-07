/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       StructureAnalysis: Identifying patterns in plant architecture and development
 *
 *       Copyright 1995-2018 CIRAD AGAP
 *
 *       File author(s): Yann Guedon (yann.guedon@cirad.fr)
 *
 *       $Source$
 *       $Id$
 *
 *       Forum for StructureAnalysis developers:
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



#include <sstream>

#include "variable_order_markov.h"
#include "sequence_label.h"

using namespace std;
using namespace stat_tool;


namespace sequence_analysis {



/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the state probabilities as a function of
 *         the index parameter for a variable-order Markov chain.
 */
/*--------------------------------------------------------------*/

void VariableOrderMarkovChain::index_state_distribution()

{
  int i , j , k;
  double *memory , *previous_memory;
  Curves *index_state;


  index_state = state_process->index_value;

  // initialization of the probabilities of the memories and the states

  memory = new double[nb_row];
  previous_memory = new double[nb_row];

  switch (type) {

  case ORDINARY : {
    for (i = 1;i < nb_row;i++) {
      if (order[i] == 1) {
        index_state->point[state[i][0]][0] = initial[state[i][0]];
        memory[i] = initial[state[i][0]];
      }
      else {
        memory[i] = 0.;
      }
    }
    break;
  }

  case EQUILIBRIUM : {
    for (i = 0;i < nb_state;i++) {
      index_state->point[i][0] = 0.;
    }

    for (i = 1;i < nb_row;i++) {
      if (!child[i]) {
        index_state->point[state[i][0]][0] += initial[i];
        memory[i] = initial[i];
      }
      else {
        memory[i] = 0.;
      }
    }
    break;
  }
  }

  // computation of the state probabilities as a function of the index parameter

  for (i = 1;i < index_state->length;i++) {

    // update of the probabilities of the memories

    for (j = 1;j < nb_row;j++) {
      previous_memory[j] = memory[j];
    }

    // computation of the probabilities of the memories

    for (j = 1;j < nb_row;j++) {
      memory[j] = 0.;
      for (k = 0;k < nb_memory[j];k++) {
        memory[j] += transition[previous[j][k]][state[j][0]] * previous_memory[previous[j][k]];
      }
    }

    // computation of the state probabilities

    for (j = 0;j < nb_state;j++) {
      index_state->point[j][i] = 0.;
    }
    for (j = 1;j < nb_row;j++) {
      index_state->point[state[j][0]][i] += memory[j];
    }
  }

  delete [] memory;
  delete [] previous_memory;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the probabilities of the memories for an ordinary variable-order Markov chain
 *         taking account of the sequence length distribution.
 *
 *  \return memory probabilities.
 */
/*--------------------------------------------------------------*/

double* VariableOrderMarkovChain::memory_computation() const

{
  int i , j , k;
  double *average_memory , *memory , *previous_memory;


  average_memory = new double[nb_row];
  for (i = 1;i < nb_row;i++) {
    average_memory[i] = 0.;
  }

  // initialization of the probabilities of the memories

  memory = new double[nb_row];
  previous_memory = new double[nb_row];

  for (i = 1;i < nb_row;i++) {
    if (order[i] == 1) {
//      average_memory[i] += initial[state[i][0]];
      memory[i] = initial[state[i][0]];
    }
    else {
      memory[i] = 0.;
    }
  }

  // computation of the probabilities of the memories as a function of the index parameter

  for (i = 1;i < state_process->length->nb_value - 2;i++) {

    // update of the probabilities of the memories

    for (j = 1;j < nb_row;j++) {
      previous_memory[j] = memory[j];
    }

    // computation of the probabilities of the memories

    for (j = 1;j < nb_row;j++) {
      memory[j] = 0.;
      for (k = 0;k < nb_memory[j];k++) {
        memory[j] += transition[previous[j][k]][state[j][0]] * previous_memory[previous[j][k]];
      }
    }

    // accumulation of the probabilities of the memories

    for (j = 1;j < nb_row;j++) {
      average_memory[j] += memory[j] * (1. - state_process->length->cumul[i]);
    }
  }

  delete [] memory;
  delete [] previous_memory;

  return average_memory;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the probability of not visiting a state
 *         for an ordinary variable-order Markov chain.
 *
 *  \param[in] istate    state,
 *  \param[in] increment threshold on the sum of the probabilities of the memories.
 */
/*--------------------------------------------------------------*/

void VariableOrderMarkovChain::state_no_occurrence_probability(int istate , double increment)

{
  int i;

  for (i = 0;i < nb_state;i++) {
    if ((i != istate) && (!accessibility[i][istate])) {
      break;
    }
  }

  if (i < nb_state) {
    int j , k;
    double memory_sum , *memory , *previous_memory ,
           &no_occurrence = state_process->no_occurrence[istate];


    // initialization of the probabilities of the memories

    memory = new double[nb_row];
    previous_memory = new double[nb_row];

    memory_sum = 0.;
    no_occurrence = 0.;

    for (i = 1;i < nb_row;i++) {
      if (order[i] == 1) {
        if (state[i][0] != istate) {
          if (accessibility[state[i][0]][istate]) {
            memory[i] = initial[state[i][0]];
            memory_sum += memory[i];
          }
          else {
            memory[i] = 0.;
            no_occurrence += initial[state[i][0]];
          }
        }

        else {
          memory[i] = 0.;
        }
      }

      else {
        memory[i] = 0.;
      }
    }

    i = 1;

    while ((memory_sum > increment) || (i < (nb_state - 1) * max_order)) {

      // update of the probabilities of the memories

      for (j = 1;j < nb_row;j++) {
        previous_memory[j] = memory[j];
      }

      // computation of the probabilities of the memories and update of
      // the probability of not visiting the selected state

      memory_sum = 0.;

      for (j = 1;j < nb_row;j++) {
        memory[j] = 0.;

        for (k = 1;k < MIN(i , order[j]);k++) {
          if ((state[j][k] == istate) || (!accessibility[state[j][k]][istate])) {
            break;
          }
        }

        if ((k == MIN(i , order[j])) && (state[j][0] != istate)) {
          if (accessibility[state[j][0]][istate]) {
            for (k = 0;k < nb_memory[j];k++) {
              if (state[previous[j][k]][0] != istate) {
                memory[j] += transition[previous[j][k]][state[j][0]] * previous_memory[previous[j][k]];
              }
            }

            memory_sum += memory[j];
          }

          else {
            for (k = 0;k < nb_memory[j];k++) {
              if (state[previous[j][k]][0] != istate) {
                no_occurrence += transition[previous[j][k]][state[j][0]] * previous_memory[previous[j][k]];
              }
            }
          }
        }
      }

      i++;
    }

    delete [] memory;
    delete [] previous_memory;
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the distribution of the time to the 1st occurrence of
 *         a state for a variable-order Markov chain.
 *
 *  \param[in] istate          state,
 *  \param[in] min_nb_value    minimum number of values,
 *  \param[in] cumul_threshold threshold on the cumulative distribution function.
 */
/*--------------------------------------------------------------*/

void VariableOrderMarkovChain::state_first_occurrence_distribution(int istate , int min_nb_value ,
                                                                   double cumul_threshold)

{
  int i , j , k;
  double *memory , *previous_memory , *pmass , *pcumul;
  Distribution *first_occurrence;


  first_occurrence = state_process->first_occurrence[istate];
  first_occurrence->complement = state_process->no_occurrence[istate];

  pmass = first_occurrence->mass;
  pcumul = first_occurrence->cumul;

  // initialization of the probabilities of the memories

  memory = new double[nb_row];
  previous_memory = new double[nb_row];

  switch (type) {

  case ORDINARY : {
    for (i = 1;i < nb_row;i++) {
      if (order[i] == 1) {
        if (state[i][0] != istate) {
          memory[i] = initial[state[i][0]];
        }
        else {
          memory[i] = 0.;
          *pmass = initial[state[i][0]];
        }
      }

      else {
        memory[i] = 0.;
      }
    }
    break;
  }

  case EQUILIBRIUM : {
    *pmass = 0.;

    for (i = 1;i < nb_row;i++) {
      if (!child[i]) {
        if (state[i][0] != istate) {
          memory[i] = initial[i];
        }
        else {
          memory[i] = 0.;
          *pmass += initial[i];
        }
      }

      else {
        memory[i] = 0.;
      }
    }
    break;
  }
  }

  *pcumul = *pmass;

  i = 1;

  while (((*pcumul < cumul_threshold - first_occurrence->complement) || (i < min_nb_value)) &&
         (i < first_occurrence->alloc_nb_value)) {

    // update of the probabilities of the memories

    for (j = 1;j < nb_row;j++) {
      previous_memory[j] = memory[j];
    }

    // computation of the probabilities of the memories and the current probability mass

    *++pmass = 0.;

    for (j = 1;j < nb_row;j++) {
      memory[j] = 0.;

      for (k = 1;k < MIN(i , order[j]);k++) {
        if (state[j][k] == istate) {
          break;
        }
      }

      if (k == MIN(i , order[j])) {
        if (state[j][0] != istate) {
          for (k = 0;k < nb_memory[j];k++) {
            if (state[previous[j][k]][0] != istate) {
              memory[j] += transition[previous[j][k]][state[j][0]] * previous_memory[previous[j][k]];
            }
          }
        }

        else {
          for (k = 0;k < nb_memory[j];k++) {
            if (state[previous[j][k]][0] != istate) {
              *pmass += transition[previous[j][k]][state[j][0]] * previous_memory[previous[j][k]];
            }
          }
        }
      }
    }

    // update of the cumulative distribution function

    pcumul++;
    *pcumul = *(pcumul - 1) + *pmass;
    i++;
  }

  first_occurrence->nb_value = i;

# ifdef DEBUG
  if (first_occurrence->complement > 0.) {
    cout << "\n" << SEQ_label[SEQL_STATE_NO_OCCURRENCE] << " " << istate << " : "
         << first_occurrence->complement << " | "
         << 1. - first_occurrence->cumul[first_occurrence->nb_value - 1] << endl;
  }
# endif

  first_occurrence->offset_computation();
  first_occurrence->max_computation();
  first_occurrence->mean_computation();
  first_occurrence->variance_computation();

  delete [] memory;
  delete [] previous_memory;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the probability of leaving definitively a state
 *         for an ordinary variable-order Markov chain.
 *
 *  \param[in] imemory   memory distribution,
 *  \param[in] istate    state,
 *  \param[in] increment threshold on the sum of the probabilities of the memories.
 */
/*--------------------------------------------------------------*/

void VariableOrderMarkovChain::state_leave_probability(const double *imemory , int istate ,
                                                       double increment)

{
  if (stype[istate] == TRANSIENT) {
    int i , j , k;
    double memory_sum , *memory , *previous_memory ,
           &leave = state_process->leave[istate];


    memory = new double[nb_row];
    previous_memory = new double[nb_row];

    // initialization of the probabilities of the memories

    memory_sum = 0.;

    for (i = 1;i < nb_row;i++) {
      if (state[i][0] == istate) {
        memory[i] = imemory[i];
        if (order[i] == 1) {
          memory[i] += initial[state[i][0]];
        }
        memory_sum += memory[i];
      }

      else {
        memory[i] = 0.;
      }
    }

    for (i = 1;i < nb_row;i++) {
      if (state[i][0] == istate) {
        memory[i] /= memory_sum;
      }
    }

    leave = 0.;
    i = 1;

    do {

      // update of the probabilities of the memories

      for (j = 1;j < nb_row;j++) {
        previous_memory[j] = memory[j];
      }

      // computation of the probabilities of the memories and update of
      // the probability of leaving definitively the selected state

      memory_sum = 0.;

      for (j = 1;j < nb_row;j++) {
        memory[j] = 0.;

        for (k = 1;k < MIN(i , order[j]);k++) {
          if ((state[j][k] == istate) || (!accessibility[state[j][k]][istate])) {
            break;
          }
        }

        if ((((k == i) && (i < order[j]) && (state[j][k] == istate)) ||
             (k == order[j])) && (state[j][0] != istate)) {
          if (accessibility[state[j][0]][istate]) {
            for (k = 0;k < nb_memory[j];k++) {
              memory[j] += transition[previous[j][k]][state[j][0]] * previous_memory[previous[j][k]];
            }
            memory_sum += memory[j];
          }

          else {
            for (k = 0;k < nb_memory[j];k++) {
              leave += transition[previous[j][k]][state[j][0]] * previous_memory[previous[j][k]];
            }
          }
        }
      }

      i++;
    }
    while ((memory_sum > increment) || (i < (nb_state - 1) * max_order));

    delete [] memory;
    delete [] previous_memory;
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the distribution of the recurrence time in a state
 *         for a variable-order Markov chain.
 *
 *  \param[in] imemory         memory distribution,
 *  \param[in] istate          state,
 *  \param[in] min_nb_value    minimum number of values,
 *  \param[in] cumul_threshold threshold on the cumulative distribution function.
 */
/*--------------------------------------------------------------*/

void VariableOrderMarkovChain::state_recurrence_time_distribution(const double *imemory , int istate ,
                                                                  int min_nb_value , double cumul_threshold)

{
  int i , j , k;
  double sum , *memory , *previous_memory , *pmass , *pcumul;
  Distribution *recurrence_time;


  recurrence_time = state_process->recurrence_time[istate];
  recurrence_time->complement = state_process->leave[istate];

  pmass = recurrence_time->mass;
  pcumul = recurrence_time->cumul;
  *pmass = 0.;
  *pcumul = 0.;

  memory = new double[nb_row];
  previous_memory = new double[nb_row];

  // initialization of the probabilities of the memories

  sum = 0.;

  for (i = 1;i < nb_row;i++) {
    if (state[i][0] == istate) {
      memory[i] = imemory[i];
      if ((type == ORDINARY) && (order[i] == 1)) {
        memory[i] += initial[state[i][0]];
      }
      sum += memory[i];
    }

    else {
      memory[i] = 0.;
    }
  }

  for (i = 1;i < nb_row;i++) {
    if (state[i][0] == istate) {
      memory[i] /= sum;
    }
  }

  i = 1;

  do {

    // update of the probabilities of the memories

    for (j = 1;j < nb_row;j++) {
      previous_memory[j] = memory[j];
    }

    // computation of the probabilities of the memories and the current probability mass

    *++pmass = 0.;

    for (j = 1;j < nb_row;j++) {
      memory[j] = 0.;

      for (k = 1;k < MIN(i , order[j]);k++) {
        if (state[j][k] == istate) {
          break;
        }
      }

      if (((k == i) && (i < order[j]) && (state[j][k] == istate)) || (k == order[j])) {
        if (state[j][0] != istate) {
          for (k = 0;k < nb_memory[j];k++) {
            memory[j] += transition[previous[j][k]][state[j][0]] * previous_memory[previous[j][k]];
          }
        }

        else {
          for (k = 0;k < nb_memory[j];k++) {
            *pmass += transition[previous[j][k]][state[j][0]] * previous_memory[previous[j][k]];
          }
        }
      }
    }

    // update of the cumulative distribution function

    pcumul++;
    *pcumul = *(pcumul - 1) + *pmass;
    i++;
  }
  while (((*pcumul < cumul_threshold - recurrence_time->complement) || (i < min_nb_value)) &&
         (i < recurrence_time->alloc_nb_value));

  recurrence_time->nb_value = i;
  recurrence_time->nb_value_computation();

  if (recurrence_time->nb_value > 0) {
    recurrence_time->offset_computation();
    recurrence_time->max_computation();
    recurrence_time->mean_computation();
    recurrence_time->variance_computation();
  }

  else {
    delete state_process->recurrence_time[istate];
    state_process->recurrence_time[istate] = NULL;
  }

  delete [] memory;
  delete [] previous_memory;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the distribution of the sojourn time in a state
 *         for a variable-order Markov chain.
 *
 *  \param[in] imemory         memory distribution,
 *  \param[in] istate          state,
 *  \param[in] min_nb_value    minimum number of values,
 *  \param[in] cumul_threshold threshold on the cumulative distribution function.
 */
/*--------------------------------------------------------------*/

void VariableOrderMarkovChain::state_sojourn_time_distribution(const double *imemory , int istate ,
                                                               int min_nb_value , double cumul_threshold)

{
  int i , j , k;
  int self_index;
  double sum , *memory , *previous_memory , *pmass , *pcumul;
  DiscreteParametric *sojourn_time;


  sojourn_time = state_process->sojourn_time[istate];

  pmass = sojourn_time->mass;
  pcumul = sojourn_time->cumul;
  *pmass = 0.;
  *pcumul = 0.;

  memory = new double[nb_row];
  previous_memory = new double[nb_row];

  // initialization of the probabilities of the memories

  sum = 0.;

  for (i = 1;i < nb_row;i++) {
    if ((state[i][0] == istate) && ((order[i] == 1) ||
         ((order[i] > 1) && (state[i][1] != istate)))) {
      memory[i] = imemory[i];
      if ((type == ORDINARY) && (order[i] == 1)) {
        memory[i] += initial[state[i][0]];
      }
      sum += memory[i];
    }

    else {
      memory[i] = 0.;
    }
  }

  for (i = 1;i < nb_row;i++) {
    if ((state[i][0] == istate) && ((order[i] == 1) ||
         ((order[i] > 1) && (state[i][1] != istate)))) {
      memory[i] /= sum;
    }
  }

  // sojourn time < maximum order of the Markov chain

  for (i = 1;i < max_order;i++) {

    // update of the probabilities of the memories

    for (j = 1;j < nb_row;j++) {
      previous_memory[j] = memory[j];
    }

    // computation of the probabilities of the memories and the current probability mass

    *++pmass = 0.;

    for (j = 1;j < nb_row;j++) {
      memory[j] = 0.;

      for (k = 0;k <= MIN(i , order[j] - 1);k++) {
        if (state[j][k] != istate) {
          break;
        }
      }

      if (((k == i + 1) && (i < order[j] - 1) && (state[j][k] != istate)) || (k == order[j])) {
        for (k = 0;k < nb_memory[j];k++) {
          memory[j] += transition[previous[j][k]][state[j][0]] * previous_memory[previous[j][k]];
          *pmass += (1. - transition[previous[j][k]][state[j][0]]) * previous_memory[previous[j][k]];
        }
      }
    }

    // update of the cumulative distribution function

    pcumul++;
    *pcumul = *(pcumul - 1) + *pmass;
  }

  // computation of the probability masses of the geometric tail

  for (j = 1;j < nb_row;j++) {
    if (!child[j]) {
      for (k = 0;k < order[j];k++) {
        if (state[j][k] != istate) {
          break;
        }
      }

      if (k == order[j]) {
        self_index = j;
        break;
      }
    }
  }

  while (((*pcumul < cumul_threshold) || (i < min_nb_value)) &&
         (i < sojourn_time->alloc_nb_value)) {
    *++pmass = memory[self_index] * (1. - transition[self_index][istate]);
    memory[self_index] *= transition[self_index][istate];
    pcumul++;
    *pcumul = *(pcumul - 1) + *pmass;
    i++;
  }
  sojourn_time->nb_value = i;

  while (*pmass-- == 0.) {
    (sojourn_time->nb_value)--;
  }

  sojourn_time->offset_computation();
  sojourn_time->max_computation();
  sojourn_time->mean_computation();
  sojourn_time->variance_computation();

  delete [] memory;
  delete [] previous_memory;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the observation probabilities as a function of
 *         the index parameter for a hidden variable-order Markov chain.
 *
 *  \param[in] variable observation process index.
 */
/*--------------------------------------------------------------*/

void VariableOrderMarkov::index_output_distribution(int variable)

{
  int i , j , k;
  Curves *index_state , *index_value;


  index_value = categorical_process[variable]->index_value;

  // computation of the state probabilities

  if (!(state_process->index_value)) {
    state_process->index_value = new Curves(nb_state , index_value->length);
    index_state_distribution();
  }
  index_state = state_process->index_value;

  // incorporation of the observation probabilities

  for (i = 0;i < index_value->length;i++) {
    for (j = 0;j < categorical_process[variable]->nb_value;j++) {
      index_value->point[j][i] = 0.;
      for (k = 0;k < nb_state;k++) {
        index_value->point[j][i] += categorical_process[variable]->observation[k]->mass[j] *
                                    index_state->point[k][i];
      }
    }
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the probability of not observing a value for
 *         for a hidden ordinary variable-order Markov chain.
 *
 *  \param[in] variable  observation process index,
 *  \param[in] output    observation,
 *  \param[in] increment threshold on the sum of the probabilities of the memories.
 */
/*--------------------------------------------------------------*/

void VariableOrderMarkov::output_no_occurrence_probability(int variable , int output ,
                                                           double increment)

{
  bool status = false , *output_accessibility;
  int i , j , k;
  double memory_sum , sum , *observation , *memory , *previous_memory ,
         &no_occurrence = categorical_process[variable]->no_occurrence[output];


  observation = new double[nb_state];
  for (i = 0;i < nb_state;i++) {
    observation[i] = categorical_process[variable]->observation[i]->mass[output];
  }

  // computation of the accessibility of the selected observation from a given state

  output_accessibility = new bool[nb_state];

  for (i = 0;i < nb_state;i++) {
    output_accessibility[i] = false;

    for (j = 0;j < nb_state;j++) {
      if (j == i) {
        if (observation[j] > 0.) {
          output_accessibility[i] = true;
          break;
        }
      }

      else {
        if ((accessibility[i][j]) && (observation[j] > 0.)) {
          output_accessibility[i] = true;
          break;
        }
      }
    }

    if (!output_accessibility[i]) {
      status = true;
    }
  }

  if (status) {

    // initialization of the probabilities of the memories

    memory = new double[nb_row];
    previous_memory = new double[nb_row];

    memory_sum = 0.;
    no_occurrence = 0.;

    for (i = 1;i < nb_row;i++) {
      if (order[i] == 1) {
        if (output_accessibility[state[i][0]]) {
          memory[i] = (1. - observation[state[i][0]]) * initial[state[i][0]];
          memory_sum += memory[i];
        }
        else {
          memory[i] = 0.;
          no_occurrence += initial[state[i][0]];
        }
      }

      else {
        memory[i] = 0.;
      }
    }

    i = 1;

    while ((memory_sum > increment) || (i < nb_state * max_order)) {

      // update of the probabilities of the memories

      for (j = 1;j < nb_row;j++) {
        previous_memory[j] = memory[j];
      }

      // computation of the probabilities of the memories and update of
      // the probability of not observing the selected observation

      memory_sum = 0.;

      for (j = 1;j < nb_row;j++) {
        sum = 0.;
        for (k = 0;k < nb_memory[j];k++) {
          sum += transition[previous[j][k]][state[j][0]] * previous_memory[previous[j][k]];
        }

        if (output_accessibility[state[j][0]]) {
          memory[j] = (1. - observation[state[j][0]]) * sum;
          memory_sum += memory[j];
        }
        else {
          memory[j] = 0.;
          no_occurrence += sum;
        }
      }

      i++;
    }

    delete [] memory;
    delete [] previous_memory;
  }

  delete [] observation;
  delete [] output_accessibility;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the distribution of the time to the 1st occurrence of
 *         a categorical observation for a hidden variable-order Markov chain.
 *
 *  \param[in] variable        observation process index,
 *  \param[in] output          observation,
 *  \param[in] min_nb_value    minimum number of values,
 *  \param[in] cumul_threshold threshold on the cumulative distribution function.
 */
/*--------------------------------------------------------------*/

void VariableOrderMarkov::output_first_occurrence_distribution(int variable , int output ,
                                                               int min_nb_value ,
                                                               double cumul_threshold)

{
  int i , j , k;
  double sum , *observation , *memory , *previous_memory , *pmass , *pcumul;
  Distribution *first_occurrence;


  first_occurrence = categorical_process[variable]->first_occurrence[output];
  first_occurrence->complement = categorical_process[variable]->no_occurrence[output];

  pmass = first_occurrence->mass;
  pcumul = first_occurrence->cumul;

  observation = new double[nb_state];
  for (i = 0;i < nb_state;i++) {
    observation[i] = categorical_process[variable]->observation[i]->mass[output];
  }

  // initialization of the probabilities of the memories

  memory = new double[nb_row];
  previous_memory = new double[nb_row];

  *pmass = 0.;

  switch (type) {

  case ORDINARY : {
    for (i = 1;i < nb_row;i++) {
      if (order[i] == 1) {
        memory[i] = (1. - observation[state[i][0]]) * initial[state[i][0]];
        *pmass += observation[state[i][0]] * initial[state[i][0]];
      }
      else {
        memory[i] = 0.;
      }
    }
    break;
  }

  case EQUILIBRIUM : {
    for (i = 1;i < nb_row;i++) {
      if (!child[i]) {
        memory[i] = (1. - observation[state[i][0]]) * initial[i];
        *pmass += observation[state[i][0]] * initial[i];
      }
      else {
        memory[i] = 0.;
      }
    }
    break;
  }
  }

  *pcumul = *pmass;

  i = 1;

  while (((*pcumul < cumul_threshold - first_occurrence->complement) || (i < min_nb_value)) &&
         (i < first_occurrence->alloc_nb_value)) {

    // update of the probabilities of the memories

    for (j = 1;j < nb_row;j++) {
      previous_memory[j] = memory[j];
    }

    // computation of the probabilities of the memories and the current probability mass

    *++pmass = 0.;

    for (j = 1;j < nb_row;j++) {
      sum = 0.;
      for (k = 0;k < nb_memory[j];k++) {
        sum += transition[previous[j][k]][state[j][0]] * previous_memory[previous[j][k]];
      }

      memory[j] = (1. - observation[state[j][0]]) * sum;
      *pmass += observation[state[j][0]] * sum;
    }

    // update of the cumulative distribution function

    pcumul++;
    *pcumul = *(pcumul - 1) + *pmass;
    i++;
  }

  first_occurrence->nb_value = i;

  first_occurrence->offset_computation();
  first_occurrence->max_computation();
  first_occurrence->mean_computation();
  first_occurrence->variance_computation();

  delete [] observation;
  delete [] memory;
  delete [] previous_memory;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the probability of leaving definitively a categorical observation
 *         for a hidden ordinary variable-order Markov chain.
 *
 *  \param[in] imemory   memory distribution,
 *  \param[in] variable  observation process index,
 *  \param[in] output    observation,
 *  \param[in] increment threshold on the sum of the probabilities of the memories.
 */
/*--------------------------------------------------------------*/

void VariableOrderMarkov::output_leave_probability(const double *imemory , int variable ,
                                                   int output , double increment)

{
  bool status = false , *output_accessibility;
  int i , j , k;
  double memory_sum , sum , *observation , *memory , *previous_memory ,
         &leave = categorical_process[variable]->leave[output];


  observation = new double[nb_state];
  for (i = 0;i < nb_state;i++) {
    observation[i] = categorical_process[variable]->observation[i]->mass[output];
  }

  // computation of the accessibility of the selected observation from a given state

  output_accessibility = new bool[nb_state];

  for (i = 0;i < nb_state;i++) {
    output_accessibility[i] = false;

    for (j = 0;j < nb_state;j++) {
      if (j == i) {
        if (observation[j] > 0.) {
          output_accessibility[i] = true;
          break;
        }
      }

      else {
        if ((accessibility[i][j]) && (observation[j] > 0.)) {
          output_accessibility[i] = true;
          break;
        }
      }
    }

    if (!output_accessibility[i]) {
      status = true;
    }
  }

  if (status) {
    memory = new double[nb_row];
    previous_memory = new double[nb_row];

    // initialization of the probabilities of the memories

    memory_sum = 0.;

    for (i = 1;i < nb_row;i++) {
      memory[i] = imemory[i];
      if (order[i] == 1) {
        memory[i] += initial[state[i][0]];
      }
      memory[i] *= observation[state[i][0]];

      memory_sum += memory[i];
    }

    for (i = 1;i < nb_row;i++) {
      memory[i] /= memory_sum;
    }

    leave = 0.;
    i = 1;

    do {

      // update of the probabilities of the memories

      for (j = 1;j < nb_row;j++) {
        previous_memory[j] = memory[j];
      }

      // computation of the probabilities of the memories and update of
      // the probability of leaving definitively the selected observation

      memory_sum = 0.;

      for (j = 1;j < nb_row;j++) {
        memory[j] = 0.;

        if (observation[state[j][0]] < 1.) {
          sum = 0.;
          for (k = 0;k < nb_memory[j];k++) {
            sum += transition[previous[j][k]][state[j][0]] * previous_memory[previous[j][k]];
          }

          if (output_accessibility[state[j][0]]) {
            memory[j] = (1. - observation[state[j][0]]) * sum;
            memory_sum += memory[j];
          }
          else {
            leave += sum;
          }
        }
      }

      i++;
    }
    while ((memory_sum > increment) || (i < nb_state * max_order));

    delete [] memory;
    delete [] previous_memory;
  }

  delete [] observation;
  delete [] output_accessibility;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the distribution of the recurrence time in a categorical observation
 *         for a hidden variable-order Markov chain.
 *
 *  \param[in] imemory         memory distribution,
 *  \param[in] variable        observation process index,
 *  \param[in] output          observation,
 *  \param[in] min_nb_value    minimum number of values,
 *  \param[in] cumul_threshold threshold on the cumulative distribution function.
 */
/*--------------------------------------------------------------*/

void VariableOrderMarkov::output_recurrence_time_distribution(const double *imemory , int variable ,
                                                              int output , int min_nb_value ,
                                                              double cumul_threshold)

{
  int i , j , k;
  double sum , *observation , *memory , *previous_memory , *pmass , *pcumul;
  Distribution *recurrence_time;


  recurrence_time = categorical_process[variable]->recurrence_time[output];
  recurrence_time->complement = categorical_process[variable]->leave[output];

  pmass = recurrence_time->mass;
  pcumul = recurrence_time->cumul;
  *pmass = 0.;
  *pcumul = 0.;

  observation = new double[nb_state];
  for (i = 0;i < nb_state;i++) {
    observation[i] = categorical_process[variable]->observation[i]->mass[output];
  }

  memory = new double[nb_row];
  previous_memory = new double[nb_row];

  // initialization of the probabilities of the memories

  sum = 0.;

  for (i = 1;i < nb_row;i++) {
    memory[i] = imemory[i];
    if ((type == ORDINARY) && (order[i] == 1)) {
      memory[i] += initial[state[i][0]];
    }
    memory[i] *= observation[state[i][0]];

    sum += memory[i];
  }

  for (i = 1;i < nb_row;i++) {
    memory[i] /= sum;
  }

  i = 1;

  do {

    // update of the probabilities of the memories

    for (j = 1;j < nb_row;j++) {
      previous_memory[j] = memory[j];
    }

    // computation of the probabilities of the memories and the current probability mass

    *++pmass = 0.;

    for (j = 1;j < nb_row;j++) {
      sum = 0.;
      for (k = 0;k < nb_memory[j];k++) {
        sum += transition[previous[j][k]][state[j][0]] * previous_memory[previous[j][k]];
      }

      memory[j] = (1. - observation[state[j][0]]) * sum;
      *pmass += observation[state[j][0]] * sum;
    }

    // update of the cumulative distribution function

    pcumul++;
    *pcumul = *(pcumul - 1) + *pmass;
    i++;
  }
  while (((*pcumul < cumul_threshold - recurrence_time->complement) || (i < min_nb_value)) &&
         (i < recurrence_time->alloc_nb_value));

  recurrence_time->nb_value = i;
  recurrence_time->nb_value_computation();

  if (recurrence_time->nb_value > 0) {
    recurrence_time->offset_computation();
    recurrence_time->max_computation();
    recurrence_time->mean_computation();
    recurrence_time->variance_computation();
  }

  else {
    delete categorical_process[variable]->recurrence_time[output];
    categorical_process[variable]->recurrence_time[output] = NULL;
  }

  delete [] observation;
  delete [] memory;
  delete [] previous_memory;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the distribution of the sojourn time in a categorical observation
 *         for a hidden variable-order Markov chain.
 *
 *  \param[in] imemory         memory distribution,
 *  \param[in] variable        observation process index,
 *  \param[in] output          observation,
 *  \param[in] min_nb_value    minimum number of values,
 *  \param[in] cumul_threshold threshold on the cumulative distribution function.
 */
/*--------------------------------------------------------------*/

void VariableOrderMarkov::output_sojourn_time_distribution(const double *imemory , int variable ,
                                                           int output , int min_nb_value ,
                                                           double cumul_threshold)

{
  int i , j , k;
  double sum , *observation , *memory , *previous_memory , *pmass , *pcumul ,
         &absorption = categorical_process[variable]->absorption[output];
  DiscreteParametric *sojourn_time;


  sojourn_time = categorical_process[variable]->sojourn_time[output];

  pmass = sojourn_time->mass;
  pcumul = sojourn_time->cumul;
  *pmass = 0.;
  *pcumul = 0.;

  observation = new double[nb_state];
  for (i = 0;i < nb_state;i++) {
    observation[i] = categorical_process[variable]->observation[i]->mass[output];
  }

  memory = new double[nb_row];
  previous_memory = new double[nb_row];

  // initialization of the probabilities of the memories

  sum = 0.;

  for (i = 1;i < nb_row;i++) {
    if (order[i] == 1) {
      memory[i] = 0.;
      for (j = 0;j < nb_memory[i];j++) {
        memory[i] += (1. - observation[state[previous[i][j]][0]]) * transition[previous[i][j]][state[i][0]] *
                     imemory[previous[i][j]];
      }

      if (type == ORDINARY) {
        memory[i] += initial[state[i][0]];
        for (j = 0;j < nb_memory[i];j++) {
          memory[i] += (1. - observation[state[previous[i][j]][0]]) * transition[previous[i][j]][state[i][0]] *
                       initial[state[previous[i][j]][0]];
        }
      }

      memory[i] *= observation[state[i][0]];
    }

    else {
      memory[i] = (1. - observation[state[i][1]]) * observation[state[i][0]] * imemory[i];
    }

    sum += memory[i];
  }

  for (i = 1;i < nb_row;i++) {
    memory[i] /= sum;
  }

  i = 1;

  do {

    // update of the probabilities of the memories

    for (j = 1;j < nb_row;j++) {
      previous_memory[j] = memory[j];
    }

    // computation of the probabilities of the memories and the current probability mass

    absorption = 0.;
    *++pmass = 0.;

    for (j = 1;j < nb_row;j++) {
      sum = 0.;
      for (k = 0;k < nb_memory[j];k++) {
        sum += transition[previous[j][k]][state[j][0]] * previous_memory[previous[j][k]];
      }

      if ((stype[state[j][0]] == ABSORBING) && (observation[state[j][0]] == 1.)) {
        absorption += sum;
      }

      memory[j] = observation[state[j][0]] * sum;
      *pmass += (1. - observation[state[j][0]]) * sum;
    }

    // update of the cumulative distribution function

    pcumul++;
    *pcumul = *(pcumul - 1) + *pmass;
    i++;
  }
  while (((*pcumul < cumul_threshold - absorption) || (i < min_nb_value)) &&
         (i < sojourn_time->alloc_nb_value));

  if (*pcumul == 0.) {
    absorption = 1.;
    delete categorical_process[variable]->sojourn_time[output];
    categorical_process[variable]->sojourn_time[output] = NULL;
  }

  else {
    sojourn_time->nb_value = i;
    sojourn_time->complement = absorption;

#   ifdef DEBUG
    if (absorption > 0.) {
      cout << "\n" << SEQ_label[SEQL_OUTPUT_ABSORPTION] << " " << output << " : "
           << absorption << " | " << 1. - sojourn_time->cumul[sojourn_time->nb_value - 1] << endl;
    }
#   endif

    sojourn_time->offset_computation();
    sojourn_time->max_computation();
    sojourn_time->mean_computation();
    sojourn_time->variance_computation();
  }

  delete [] observation;
  delete [] memory;
  delete [] previous_memory;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the autocorrelation function for a state
 *         of a variable-order Markov chain (binarized state process).
 *
 *  \param[in] error   reference on a StatError object,
 *  \param[in] istate  state,
 *  \param[in] max_lag maximum lag,
 *  \param[in] seq     pointer on a VariableOrderMarkovData object.
 *
 *  \return            Correlation object.
 */
/*--------------------------------------------------------------*/

Correlation* VariableOrderMarkovChain::state_autocorrelation_computation(StatError &error ,
                                                                         int istate , int max_lag ,
                                                                         const MarkovianSequences *seq) const

{
  bool status = true;
  int i , j , k;
  int *category;
  double sum , norm , mean , *average_memory , *memory , *previous_memory , *ppoint;
  Correlation *correl;
  MarkovianSequences *binary_seq;


  correl = NULL;
  error.init();

/*  if (nb_component > 1) {
    status = false;
    error.correction_update(STAT_parsing[STATP_CHAIN_STRUCTURE] , STAT_parsing[STATP_IRREDUCIBLE]);
  } */

  if ((istate < 0) || (istate >= nb_state)) {
    status = false;
    ostringstream error_message;
    error_message << STAT_label[STATL_STATE] << " " << istate << " "
                  << STAT_error[STATR_NOT_PRESENT];
    error.update((error_message.str()).c_str());
  }

  else {
    for (i = 0;i < nb_component;i++) {
      if ((component_nb_state[i] == 1) && (component[i][0] == istate)) {
        status = false;
        error.update(SEQ_error[SEQR_SINGLE_STATE_COMPONENT]);
        break;
      }
    }
  }

  if ((max_lag < max_order) || (max_lag > MAX_LAG)) {
    status = false;
    error.update(SEQ_error[SEQR_MAX_LAG]);
  }

  if (status) {
    if ((seq) && (seq->type[0] == STATE)) {
      correl = new Correlation(2 , max_lag + 1 , true , PEARSON);
    }
    else {
      correl = new Correlation(1 , max_lag + 1 , false , PEARSON);
    }

    i = 0;
    if ((seq) && (seq->type[0] == STATE)) {
      correl->variable_type[i] = OBSERVED_STATE;
      correl->variable1[i++] = istate;
    }
    correl->variable_type[i] = THEORETICAL_STATE;
    correl->variable1[i] = istate;

    switch (type) {

    case ORDINARY : {
      average_memory = memory_computation();
      break;
    }

    case EQUILIBRIUM : {
      average_memory = new double[nb_row];
      for (i = 1;i < nb_row;i++) {
        average_memory[i] = initial[i];
      }
      break;
    }
    }

    ppoint = correl->point[((seq) && (seq->type[0] == STATE)) ? 1 : 0];
    *ppoint = 1.;

    memory = new double[nb_row];
    previous_memory = new double[nb_row];

    // initialization of the probabilities of the memories

    sum = 0.;
    norm = 0.;

    for (i = 1;i < nb_row;i++) {
      if (state[i][0] == istate) {
        memory[i] = average_memory[i];
        if ((type == ORDINARY) && (order[i] == 1)) {
          memory[i] += initial[state[i][0]];
        }
        sum += memory[i];
      }

      else {
        memory[i] = 0.;

        norm += average_memory[i];
        if ((type == ORDINARY) && (order[i] == 1)) {
          norm += initial[state[i][0]];
        }
      }
    }

    for (i = 1;i < nb_row;i++) {
      if (state[i][0] == istate) {
        memory[i] /= sum;
      }
    }

    mean = sum / (sum + norm);

    for (i = 1;i <= max_lag;i++) {

      // update of the probabilities of the memories

      for (j = 1;j < nb_row;j++) {
        previous_memory[j] = memory[j];
      }

      // computation of the probabilities of the memories and the current autocorrelation coefficient

      *++ppoint = 0.;

      for (j = 1;j < nb_row;j++) {
        memory[j] = 0.;

        if (((i < order[j]) && (state[j][i] == istate)) || (i >= order[j])) {
          for (k = 0;k < nb_memory[j];k++) {
            memory[j] += transition[previous[j][k]][state[j][0]] * previous_memory[previous[j][k]];
          }

          if (state[j][0] == istate) {
            *ppoint += memory[j];
          }
        }
      }

      *ppoint = (*ppoint - mean) / (1. - mean);
    }

    delete [] average_memory;
    delete [] memory;
    delete [] previous_memory;

    if ((seq) && (seq->type[0] == STATE)) {
      category = new int[nb_state];
      for (i = 0;i < nb_state;i++) {
        category[i] = 0;
      }
      category[istate] = 1;

      binary_seq = seq->transcode(error , 1 , category);
      binary_seq->correlation_computation(*correl , 0 , 0 , EXACT);
      delete [] category;
      delete binary_seq;
    }
  }

  return correl;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the autocorrelation function for a state
 *         of a variable-order Markov chain (binarized state process).
 *
 *  \param[in] error   reference on a StatError object,
 *  \param[in] istate  state,
 *  \param[in] max_lag maximum lag.
 *
 *  \return            Correlation object.
 */
/*--------------------------------------------------------------*/

Correlation* VariableOrderMarkov::state_autocorrelation_computation(StatError &error ,
                                                                    int istate , int max_lag) const

{
  Correlation *correl;


  correl = VariableOrderMarkovChain::state_autocorrelation_computation(error , istate , max_lag , markov_data);

  return correl;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the autocorrelation function for a state
 *         of a variable-order Markov chain (binarized state process).
 *
 *  \param[in] error   reference on a StatError object,
 *  \param[in] istate  state,
 *  \param[in] max_lag maximum lag.
 *
 *  \return            Correlation object.
 */
/*--------------------------------------------------------------*/

Correlation* VariableOrderMarkovData::state_autocorrelation_computation(StatError &error ,
                                                                        int istate , int max_lag) const

{
  Correlation *correl;


  correl = markov->VariableOrderMarkovChain::state_autocorrelation_computation(error , istate , max_lag , this);

  return correl;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the autocorrelation function for a categorical observation
 *         of a hidden variable-order Markov chain (binarized observation process).
 *
 *  \param[in] error    reference on a StatError object,
 *  \param[in] variable observation process index
 *  \param[in] output   observation,
 *  \param[in] max_lag  maximum lag,
 *  \param[in] seq      pointer on a VariableOrderMarkovData object.
 *
 *  \return             Correlation object.
 */
/*--------------------------------------------------------------*/

Correlation* VariableOrderMarkov::output_autocorrelation_computation(StatError &error , int variable ,
                                                                     int output , int max_lag ,
                                                                     const VariableOrderMarkovData *seq) const

{
  bool status = true;
  int i , j , k;
  int seq_variable , *category;
  double sum , norm , mean , *average_memory , *observation , *memory ,
         *previous_memory , *ppoint;
  Correlation *correl;
  MarkovianSequences *binary_seq;


  correl = NULL;
  error.init();

/*  if (nb_component > 1) {
    status = false;
    error.correction_update(STAT_parsing[STATP_CHAIN_STRUCTURE] , STAT_parsing[STATP_IRREDUCIBLE]);
  } */

  if (nb_output_process == 0) {
    status = false;
    error.update(STAT_error[STATR_NB_OUTPUT_PROCESS]);
  }

  else {
    if ((variable < 1) || (variable > nb_output_process) || (!categorical_process[variable - 1])) {
      status = false;
      error.update(STAT_error[STATR_OUTPUT_PROCESS_INDEX]);
    }

    else {
      variable--;

      if ((output < 0) || (output >= categorical_process[variable]->nb_value)) {
        status = false;
        ostringstream error_message;
        error_message << STAT_label[STATL_OUTPUT] << " " << output << " "
                      << STAT_error[STATR_NOT_PRESENT];
        error.update((error_message.str()).c_str());
      }
    }
  }

  if ((max_lag < max_order) || (max_lag > MAX_LAG)) {
    status = false;
    error.update(SEQ_error[SEQR_MAX_LAG]);
  }

  if (status) {
    if (seq) {
      correl = new Correlation(2 , max_lag + 1 , true , PEARSON);
    }
    else {
      correl = new Correlation(1 , max_lag + 1 , false , PEARSON);
    }

    i = 0;
    if (seq) {
      correl->variable_type[i] = OBSERVED_OUTPUT;
      correl->variable1[i++] = output;
    }
    correl->variable_type[i] = THEORETICAL_OUTPUT;
    correl->variable1[i] = output;

    switch (type) {

    case ORDINARY : {
      average_memory = memory_computation();
      break;
    }

    case EQUILIBRIUM : {
      average_memory = new double[nb_row];
      for (i = 1;i < nb_row;i++) {
        average_memory[i] = initial[i];
      }
      break;
    }
    }

    ppoint = correl->point[seq ? 1 : 0];
    *ppoint = 1.;

    observation = new double[nb_state];
    for (i = 0;i < nb_state;i++) {
      observation[i] = categorical_process[variable]->observation[i]->mass[output];
    }

    memory = new double[nb_row];
    previous_memory = new double[nb_row];

    // initialization of the probabilities of the memories

    sum = 0.;
    norm = 0.;

    for (i = 1;i < nb_row;i++) {
      memory[i] = average_memory[i];
      if ((type == ORDINARY) && (order[i] == 1)) {
        memory[i] += initial[state[i][0]];
      }
      norm += memory[i];

      memory[i] *= observation[state[i][0]];
      sum += memory[i];
    }

    for (i = 1;i < nb_row;i++) {
      memory[i] /= sum;
    }

    mean = sum / norm;

    for (i = 1;i <= max_lag;i++) {

      // update of the probabilities of the memories

      for (j = 1;j < nb_row;j++) {
        previous_memory[j] = memory[j];
      }

      // computation of the probabilities of the memories and the current autocorrelation coefficient

      *++ppoint = 0.;

      for (j = 1;j < nb_row;j++) {
        sum = 0.;
        for (k = 0;k < nb_memory[j];k++) {
          sum += transition[previous[j][k]][state[j][0]] * previous_memory[previous[j][k]];
        }

        memory[j] = sum;
        *ppoint += observation[state[j][0]] * sum;
      }

      *ppoint = (*ppoint - mean) / (1. - mean);
    }

    delete [] average_memory;
    delete [] observation;
    delete [] memory;
    delete [] previous_memory;

    if (seq) {
      switch (seq->type[0]) {
      case INT_VALUE :
        seq_variable = variable - 1;
        break;
      case STATE :
        seq_variable = variable;
        break;
      }

      category = new int[categorical_process[variable]->nb_value];
      for (i = 0;i < categorical_process[variable]->nb_value;i++) {
        category[i] = 0;
      }
      category[output] = 1;

      binary_seq = seq->transcode(error , seq_variable + 1 , category);
      binary_seq->correlation_computation(*correl , seq_variable , seq_variable , EXACT);
      delete [] category;
      delete binary_seq;
    }
  }

  return correl;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the autocorrelation function for a categorical observation
 *         of a hidden variable-order Markov chain (binarized observation process).
 *
 *  \param[in] error    reference on a StatError object,
 *  \param[in] variable observation process index,
 *  \param[in] output   observation,
 *  \param[in] max_lag  maximum lag.
 *
 *  \return             Correlation object.
 */
/*--------------------------------------------------------------*/

Correlation* VariableOrderMarkov::output_autocorrelation_computation(StatError &error ,
                                                                     int variable , int output ,
                                                                     int max_lag) const

{
  Correlation *correl;


  correl = output_autocorrelation_computation(error , variable , output , max_lag , markov_data);

  return correl;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the autocorrelation function for a categorical observation
 *         of a hidden variable-order Markov chain (binarized observation process).
 *
 *  \param[in] error    reference on a StatError object,
 *  \param[in] variable observation process index,
 *  \param[in] output   observation,
 *  \param[in] max_lag  maximum lag.
 *
 *  \return             Correlation object.
 */
/*--------------------------------------------------------------*/

Correlation* VariableOrderMarkovData::output_autocorrelation_computation(StatError &error ,
                                                                         int variable , int output ,
                                                                         int max_lag) const

{
  Correlation *correl;


  correl = markov->output_autocorrelation_computation(error , variable , output , max_lag , this);

  return correl;
}


};  // namespace sequence_analysis
