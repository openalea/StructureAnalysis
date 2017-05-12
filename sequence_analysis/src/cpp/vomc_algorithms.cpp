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

#include <boost/math/distributions/normal.hpp>

#include "stat_tool/stat_label.h"

#include "variable_order_markov.h"
#include "sequence_label.h"

using namespace std;
using namespace boost::math;
using namespace stat_tool;


namespace sequence_analysis {



/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the transition probabilities corresponding to the non-terminal memories
 *         for an ordinary variable-order Markov chain.
 */
/*--------------------------------------------------------------*/

void VariableOrderMarkovChain::non_terminal_transition_probability_computation()

{
  int i , j , k;
  int nb_terminal;
  double sum , *memory , *previous_memory;


  for (i = 1;i < nb_row;i++) {
    if (memo_type[i] == NON_TERMINAL) {
      sum = 0.;
      for (j = 0;j < nb_state;j++) {
        sum += transition[i][j];
      }
      break;
    }
  }

  if ((i < nb_row) && (sum == 0.)) {

    // initialization of the probabilities of the memories

    memory = new double[nb_row];
    previous_memory = new double[nb_row];

    for (i = 1;i < nb_row;i++) {
      if (order[i] == 1) {
        memory[i] = initial[state[i][0]];
      }
      else {
        memory[i] = 0.;
      }
    }

    // computation of the probabilities of each memory as a function of the index parameter

    for (i = 1;i < nb_row;i++) {
      if (memo_type[i] == NON_TERMINAL) {
        for (j = 0;j < nb_state;j++) {
          transition[i][j] = 1. / (double)nb_state;
        }
      }
    }

    nb_terminal = (nb_row - 1) * (nb_state - 1) / nb_state + 1;
    i = 1;

    do {

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

      // computation of the sum of the absolute differences of the probabilities of the memories

      sum = 0.;
      for (j = 1;j < nb_row;j++) {
        sum += fabs(memory[j] - previous_memory[j]);
      }

      i++;
    }
    while (((i < max_order) || (sum / nb_terminal > STATIONARY_PROBABILITY_THRESHOLD)) &&
           (i < STATIONARY_PROBABILITY_LENGTH));

#   ifdef DEBUG
    cout << "LENGTH: " << i << endl;
#   endif

    // extraction of the transition probabilities corresponding to the non-terminal memories

    for (i = nb_row - 1;i >= 1;i--) {
      if (child[i]) {
        memory[i] = 0.;

        if (memo_type[i] == NON_TERMINAL) {
          for (j = 0;j < nb_state;j++) {
            transition[i][j] = 0.;
          }
        }

        for (j = 0;j < nb_state;j++) {
          memory[i] += memory[child[i][j]];

          if (memo_type[i] == NON_TERMINAL) {
            for (k = 0;k < nb_state;k++) {
              transition[i][k] += transition[child[i][j]][k] * memory[child[i][j]];
            }
          }
        }

        if (memo_type[i] == NON_TERMINAL) {
          for (j = 0;j < nb_state;j++) {
            transition[i][j] /= memory[i];
          }
        }
      }
    }

    delete [] memory;
    delete [] previous_memory;
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the stationary distribution for an equilibrium variable-order Markov chain.
 */
/*--------------------------------------------------------------*/

void VariableOrderMarkovChain::initial_probability_computation()

{
  int i , j , k;
  int nb_terminal;
  double sum , *memory , *previous_memory;


  // initialization of the probabilities of the memories

  memory = new double[nb_row];
  previous_memory = new double[nb_row];

  for (i = 1;i < nb_row;i++) {
    if (!child[i]) {
      memory[i] = initial[i];
    }
    else {
      memory[i] = 0.;
    }
  }

  // computation of the probabilities of each memory as a function of the index parameter

  nb_terminal = (nb_row - 1) * (nb_state - 1) / nb_state + 1;
  i = 1;

  do {

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

    // computation of the sum of the absolute differences of the probabilities of the memories

    sum = 0.;
    for (j = 1;j < nb_row;j++) {
      sum += fabs(memory[j] - previous_memory[j]);
    }

    i++;
  }
  while ((i < max_order) || (sum / nb_terminal > STATIONARY_PROBABILITY_THRESHOLD));
//         && (i < STATIONARY_PROBABILITY_LENGTH));

# ifdef DEBUG
  cout << "LENGTH: " << i << endl;
# endif

  initial[0] = 0.;
  for (i = 1;i < nb_row;i++) {
    initial[i] = memory[i];
  }

  delete [] memory;
  delete [] previous_memory;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the log-likelihood of a variable-order Markov chain for sequences.
 *
 *  \param[in] seq   reference on a MarkovianSequences object,
 *  \param[in] index sequence index.
 *
 *  \return          log-likelihood.
 */
/*--------------------------------------------------------------*/

double VariableOrderMarkov::likelihood_computation(const MarkovianSequences &seq , int index) const

{
  int i , j , k;
  int nb_value , memory , start , length , *pstate , **pioutput;
  double likelihood = 0. , proba , **proutput;


  // checking of the compatibility of the model with the data

  if (nb_output_process + 1 == seq.nb_variable) {
    if (state_process->nb_value < seq.marginal_distribution[0]->nb_value) {
      likelihood = D_INF;
    }

    for (i = 0;i < nb_output_process;i++) {
      if (((categorical_process[i]) || (discrete_parametric_process[i])) &&
          (seq.marginal_distribution[i])) {
        if (categorical_process[i]) {
          nb_value = categorical_process[i]->nb_value;
        }
        else {
          nb_value = discrete_parametric_process[i]->nb_value;
        }

        if (nb_value < seq.marginal_distribution[i + 1]->nb_value) {
          likelihood = D_INF;
          break;
        }
      }
    }
  }

  else {
    likelihood = D_INF;
  }

  if (likelihood != D_INF) {
    if (nb_output_process > 0) {
      pioutput = new int*[nb_output_process];
      proutput = new double*[nb_output_process];
    }

    for (i = 0;i < seq.nb_sequence;i++) {
      if ((index == I_DEFAULT) || (index == i)) {
        switch (type) {

        case ORDINARY : {
          pstate = seq.int_sequence[i][0];
          proba = initial[*pstate];
          memory = child[0][*pstate];
          start = 1;
          break;
        }

        case EQUILIBRIUM : {
          if (max_order <= seq.length[i]) {
            for (j = 1;j < nb_row;j++) {
              if (!child[j]) {
                pstate = seq.int_sequence[i][0] + max_order;

                for (k = 0;k < order[j];k++) {
                  if (*--pstate != state[j][k]) {
                    break;
                  }
                }

                if (k == order[j]) {
                  proba = initial[j];
                  memory = j;
                  pstate = seq.int_sequence[i][0] + max_order - 1;
                  start = max_order;
                  break;
                }
              }
            }
          }

          else {
            likelihood = D_INF;
          }
          break;
        }
        }

        if (likelihood == D_INF) {
          break;
        }

        if (proba > 0.) {
          likelihood += log(proba);
        }
        else {
          likelihood = D_INF;
          break;
        }

        for (j = 0;j < nb_output_process;j++) {
          switch (seq.type[j + 1]) {
          case INT_VALUE :
            pioutput[j] = seq.int_sequence[i][j + 1] + start - 1;
            break;
          case REAL_VALUE :
            proutput[j] = seq.real_sequence[i][j + 1] + start - 1;
            break;
          }

          if (categorical_process[j]) {
            proba = categorical_process[j]->observation[*pstate]->mass[*pioutput[j]];
          }

          else if (discrete_parametric_process[j]) {
            proba = discrete_parametric_process[j]->observation[*pstate]->mass[*pioutput[j]];
          }

          else {
            switch (seq.type[j + 1]) {
            case INT_VALUE :
              proba = continuous_parametric_process[j]->observation[*pstate]->mass_computation(*pioutput[j] - seq.min_interval[j + 1] / 2 , *pioutput[j] + seq.min_interval[j + 1] / 2);
              break;
            case REAL_VALUE :
              proba = continuous_parametric_process[j]->observation[*pstate]->mass_computation(*proutput[j] - seq.min_interval[j + 1] / 2 , *proutput[j] + seq.min_interval[j + 1] / 2);
               break;
            }
          }

          if (proba > 0.) {
            likelihood += log(proba);
          }
          else {
            likelihood = D_INF;
            break;
          }
        }

        for (j = start;j < seq.length[i];j++) {
          proba = transition[memory][*++pstate];
          if (proba > 0.) {
            likelihood += log(proba);
          }
          else {
            likelihood = D_INF;
            break;
          }

          for (k = 0;k < nb_output_process;k++) {
            if (categorical_process[k]) {
              proba = categorical_process[k]->observation[*pstate]->mass[*++pioutput[k]];
            }

            else if (discrete_parametric_process[k]) {
              proba = discrete_parametric_process[k]->observation[*pstate]->mass[*++pioutput[k]];
            }

            else {
              switch (seq.type[k + 1]) {
              case INT_VALUE :
                pioutput[k]++;
                proba = continuous_parametric_process[k]->observation[*pstate]->mass_computation(*pioutput[k] - seq.min_interval[k + 1] / 2 , *pioutput[k] + seq.min_interval[k + 1] / 2);
                break;
              case REAL_VALUE :
                proutput[k]++;
                proba = continuous_parametric_process[k]->observation[*pstate]->mass_computation(*proutput[k] - seq.min_interval[k + 1] / 2 , *proutput[k] + seq.min_interval[k + 1] / 2);
                break;
              }
            }

            if (proba > 0.) {
              likelihood += log(proba);
            }
            else {
              likelihood = D_INF;
              break;
            }
          }

          if (likelihood == D_INF) {
            break;
          }

          memory = next[memory][*pstate];
        }

        if (likelihood == D_INF) {
          break;
        }
      }
    }

    if ((likelihood != D_INF) && (type == EQUILIBRIUM)) {
      length = 0;
      for (i = 0;i < seq.nb_sequence;i++) {
        if ((index == I_DEFAULT) || (index == i)) {
          length += MIN(max_order - 1 , seq.length[i]);
        }
      }

      if (index == i) {
        likelihood = likelihood * seq.length[i] / (seq.length[i] - length);
      }
      else {
        likelihood = likelihood * seq.cumul_length / (seq.cumul_length - length);
      }
    }

    if (nb_output_process > 0) {
      delete [] pioutput;
      delete [] proutput;
    }
  }

  return likelihood;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the log-likelihood of initial states and transitions
 *         for a variable-order Markov chain.
 *
 *  \param[in] chain_data reference on a VariableOrderMarkovChainData object.
 *
 *  \return               log-likelihood.
 */
/*--------------------------------------------------------------*/

double VariableOrderMarkov::likelihood_computation(const VariableOrderMarkovChainData &chain_data) const

{
  int i , j;
  double likelihood;


  if ((chain_data.nb_state != nb_state) || (chain_data.nb_row != nb_row)) {
    likelihood = D_INF;
  }

  else {
    likelihood = 0.;

    for (i = 0;i < (type == ORDINARY ? nb_state : nb_row);i++) {
      if (chain_data.initial[i] > 0) {
        if (initial[i] > 0.) {
          likelihood += chain_data.initial[i] * log(initial[i]);
        }
        else {
          likelihood = D_INF;
          break;
        }
      }
    }

    if (likelihood != D_INF) {
      for (i = 1;i < nb_row;i++) {
        if ((memo_type[i] == TERMINAL) || ((type == ORDINARY) &&
             (memo_type[i] == NON_TERMINAL))) {
          for (j = 0;j < nb_state;j++) {
            if (chain_data.transition[i][j] > 0) {
              if (transition[i][j] > 0.) {
                likelihood += chain_data.transition[i][j] * log(transition[i][j]);
              }
              else {
                likelihood = D_INF;
                break;
              }
            }
          }

          if (likelihood == D_INF) {
            break;
          }
        }
      }
    }
  }

  return likelihood;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Correction of the observation frequency distributions.
 *
 *  \param[in] corrected_observation reference on FrequencyDistribution objects,
 *  \param[in] variable              variable index,
 *  \param[in] start                 start.
 */
/*--------------------------------------------------------------*/

void VariableOrderMarkovData::observation_frequency_distribution_correction(FrequencyDistribution **corrected_observation ,
                                                                            int variable , int start) const

{
  int i , j;
  int *pstate , *poutput;


  for (i = 0;i < nb_sequence;i++) {
    pstate = int_sequence[i][0];
    poutput = int_sequence[i][variable];
    for (j = 0;j < MIN(start , length[i]);j++) {
      (corrected_observation[*pstate++]->frequency[*poutput++])--;
    }
  }

  // computation of the characteristics of the observation frequency distributions

  for (i = 0;i < marginal_distribution[0]->nb_value;i++) {
    corrected_observation[i]->nb_value_computation();
    corrected_observation[i]->offset_computation();
    corrected_observation[i]->nb_element_computation();
    corrected_observation[i]->max_computation();

    if (!characteristics[variable]) {
      corrected_observation[i]->mean_computation();
      corrected_observation[i]->variance_computation();
    }
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the log-likelihood of a variable-order Markov chain for sequences.
 *
 *  \param[in] seq reference on a VariableOrderMarkovData object.
 *
 *  \return        log-likelihood.
 */
/*--------------------------------------------------------------*/

double VariableOrderMarkov::likelihood_computation(const VariableOrderMarkovData &seq) const

{
  int i , j;
  int nb_value , length;
  double buff , likelihood = 0.;
  FrequencyDistribution **observation;


  // checking of the compatibility of the model with the data

  if (nb_output_process + 1 == seq.nb_variable) {
    if ((!(seq.marginal_distribution[0])) || (nb_state < seq.marginal_distribution[0]->nb_value)) {
      likelihood = D_INF;
    }

    for (i = 0;i < nb_output_process;i++) {
      if ((categorical_process[i]) || (discrete_parametric_process[i])) {
        if (categorical_process[i]) {
          nb_value = categorical_process[i]->nb_value;
        }
        else {
          nb_value = discrete_parametric_process[i]->nb_value;
        }

        if (nb_value < seq.marginal_distribution[i + 1]->nb_value) {
          likelihood = D_INF;
          break;
        }
      }

      else if (!(seq.marginal_distribution[i + 1])) {
        likelihood = D_INF;
        break;
      }
    }
  }

  else {
    likelihood = D_INF;
  }

  if (likelihood != D_INF) {
    likelihood = likelihood_computation(*(seq.chain_data));

    if ((likelihood != D_INF) && (nb_output_process > 0)) {
      observation = new FrequencyDistribution*[nb_state];

      for (i = 0;i < nb_output_process;i++) {
        switch (type) {

        case ORDINARY : {
          for (j = 0;j < nb_state;j++) {
            observation[j] = seq.observation_distribution[i + 1][j];
          }
          break;
        }

        case EQUILIBRIUM : {
          for (j = 0;j < nb_state;j++) {
            observation[j] = new FrequencyDistribution(*(seq.observation_distribution[i + 1][j]));
          }
          break;
        }
        }

        if (type == EQUILIBRIUM) {
          seq.observation_frequency_distribution_correction(observation , i , max_order - 1);
        }

        if (categorical_process[i]) {
          for (j = 0;j < nb_state;j++) {
            buff = categorical_process[i]->observation[j]->likelihood_computation(*observation[j]);

            if (buff != D_INF) {
              likelihood += buff;
            }
            else {
              likelihood = D_INF;
              break;
            }
          }
        }

        else if (discrete_parametric_process[i]) {
          for (j = 0;j < nb_state;j++) {
            buff = discrete_parametric_process[i]->observation[j]->likelihood_computation(*observation[j]);

            if (buff != D_INF) {
              likelihood += buff;
            }
            else {
              likelihood = D_INF;
              break;
            }
          }
        }

        else {
          for (j = 0;j < nb_state;j++) {
            buff = continuous_parametric_process[i]->observation[j]->likelihood_computation(*observation[j] , (int)seq.min_interval[i]);

            if (buff != D_INF) {
              likelihood += buff;
            }
            else {
              likelihood = D_INF;
              break;
            }
          }
        }

        if (type == EQUILIBRIUM) {
          for (j = 0;j < nb_state;j++) {
            delete observation[j];
          }
        }

        if (likelihood == D_INF) {
          break;
        }
      }

      delete [] observation;
    }

    if ((likelihood != D_INF) && (type == EQUILIBRIUM)) {
      length = 0;
      for (i = 0;i < seq.nb_sequence;i++) {
        length += MIN(max_order - 1 , seq.length[i]);
      }

      likelihood = likelihood * seq.cumul_length / (seq.cumul_length - length);
    }
  }

  return likelihood;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Counting of initial states and transitions.
 *
 *  \param[in] chain_data   reference on a VariableOrderMarkovChainData object,
 *  \param[in] markov       reference on a VariableOrderMarkovChain object,
 *  \param[in] begin        flag for taking account of the beginning of sequences,
 *  \param[in] non_terminal flag for cumulating on the non-terminal memories.
 */
/*--------------------------------------------------------------*/

void MarkovianSequences::transition_count_computation(const VariableOrderMarkovChainData &chain_data ,
                                                      const VariableOrderMarkovChain &markov ,
                                                      bool begin , bool non_terminal) const

{
  int i , j , k;
  int memory , start , *pstate;


  for (i = 0;i < (markov.type == ORDINARY ? chain_data.nb_state : chain_data.nb_row);i++) {
    chain_data.initial[i] = 0;
  }

  for (i = 0;i < chain_data.nb_row;i++) {
    for (j = 0;j < chain_data.nb_state;j++) {
      chain_data.transition[i][j] = 0;
    }
  }

  // extraction of initial states and transitions

  for (i = 0;i < nb_sequence;i++) {
    switch (markov.type) {

    case ORDINARY : {
      pstate = int_sequence[i][0];
      (chain_data.initial[*pstate])++;
      memory = markov.child[0][*pstate];
      start = 1;
      break;
    }

    case EQUILIBRIUM : {
      if (markov.max_order <= length[i]) {
        for (j = 1;j < chain_data.nb_row;j++) {
          if (!markov.child[j]) {
            pstate = int_sequence[i][0] + markov.max_order;

            for (k = 0;k < markov.order[j];k++) {
              if (*--pstate != markov.state[j][k]) {
                break;
              }
            }

            if (k == markov.order[j]) {
              (chain_data.initial[j])++;
              memory = j;
              pstate = int_sequence[i][0] + markov.max_order - 1;
              start = markov.max_order;
              break;
            }
          }
        }
      }
      break;
    }
    }

    for (j = start;j < length[i];j++) {
      pstate++;
      if ((begin) || (!markov.child[memory])) {
        (chain_data.transition[memory][*pstate])++;
      }
      memory = markov.next[memory][*pstate];
    }
  }

  // extraction of the transition counts corresponding to non-terminal memories

  for (i = chain_data.nb_row - 1;i >= 1;i--) {
    if ((markov.memo_type[i] == COMPLETION) || (non_terminal)) {
      for (j = 0;j < chain_data.nb_state;j++) {
        chain_data.transition[markov.parent[i]][j] += chain_data.transition[i][j];
      }
    }
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Construction of initial state and transition counts.
 *
 *  \param[in] markov       reference on a VariableOrderMarkovChain object,
 *  \param[in] begin        flag for taking account of the beginning of sequences,
 *  \param[in] non_terminal flag for cumulating on the non-terminal memories.
 */
/*--------------------------------------------------------------*/

void VariableOrderMarkovData::build_transition_count(const VariableOrderMarkovChain &markov ,
                                                     bool begin , bool non_terminal)

{
  chain_data = new VariableOrderMarkovChainData(markov.type , markov.nb_state , markov.nb_row);
  transition_count_computation(*chain_data , markov , begin , non_terminal);
}


/*--------------------------------------------------------------*/
/**
 *  \brief Estimation of a variable-order Markov chain from the initial state and transition counts.
 *
 *  \param[in] markov        reference on a VariableOrderMarkovChain object,
 *  \param[in] non_terminal  flag for estiming the transition probabilities on the non-terminal memories,
 *  \param[in] estimator     estimator (maximum likelihood, Laplace, adaptative Laplace),
 *  \param[in] laplace_coeff Laplace estimator coefficient.
 */
/*--------------------------------------------------------------*/

void VariableOrderMarkovChainData::estimation(VariableOrderMarkovChain &markov , bool non_terminal ,
                                              transition_estimator estimator , double laplace_coeff) const

{
  int i , j;
  int sum , nb_parameter;


  // estimation of the initial probabilities

  if (markov.type == ORDINARY) {
    sum = 0;
    for (i = 0;i < nb_state;i++) {
      sum += initial[i];
    }

    for (i = 0;i < nb_state;i++) {
      markov.initial[i] = (double)initial[i] / (double)sum;
    }
  }

  // estimation of the transition probabilities

  for (i = 0;i < nb_row;i++) {
    if ((markov.memo_type[i] == TERMINAL) || ((markov.memo_type[i] == NON_TERMINAL) &&
         ((markov.type == ORDINARY) || (non_terminal)))) {
      sum = 0;
      for (j = 0;j < nb_state;j++) {
        sum += transition[i][j];
      }

      if ((estimator == ADAPTATIVE_LAPLACE) || (estimator == UNIFORM_SUBSET) ||
          (estimator == UNIFORM_CARDINALITY)) {
        nb_parameter = 0;
        for (j = 0;j < nb_state;j++) {
          if (transition[i][j] > 0) {
            nb_parameter++;
          }
        }
      }

      switch (estimator) {

      case MAXIMUM_LIKELIHOOD : {
        if (sum > 0) {
          for (j = 0;j < nb_state;j++) {
            markov.transition[i][j] = (double)transition[i][j] / (double)sum;
          }
        }

        else if (i > 0) {
          for (j = 0;j < markov.state[i][0];j++) {
            markov.transition[i][j] = 0.;
          }
          markov.transition[i][markov.state[i][0]] = 1.;
          for (j = markov.state[i][0] + 1;j < nb_state;j++) {
            markov.transition[i][j] = 0.;
          }
        }
        break;
      }

      // Laplace estimator

      case LAPLACE : {
        for (j = 0;j < nb_state;j++) {
          markov.transition[i][j] = (double)(transition[i][j] + laplace_coeff) /
                                    (double)(sum + nb_state * laplace_coeff);
        }
        break;
      }

      // adaptative Laplace estimator (Vert, 2001)

      case ADAPTATIVE_LAPLACE : {
        if (sum > 0) {
          for (j = 0;j < nb_state;j++) {
            markov.transition[i][j] = ((double)transition[i][j] + (double)nb_parameter / (double)nb_state) /
                                      (double)(sum + nb_parameter);
          }
        }

        else {
          for (j = 0;j < nb_state;j++) {
            markov.transition[i][j] = 1 / (double)nb_state;
          }
        }
        break;
      }

      // Ristad (1995) estimator 1

      case UNIFORM_SUBSET : {
        for (j = 0;j < nb_state;j++) {
          if (transition[i][j] > 0) {
/*            markov.transition[i][j] = ((double)(transition[i][j] + 1) * (double)(sum + 1 - nb_parameter)) /
                                      ((double)(sum + nb_parameter) * (double)(sum + 1 - nb_parameter) +
                                       nb_parameter * (nb_state - nb_parameter)); */
            markov.transition[i][j] = ((double)transition[i][j] * (double)(sum + nb_parameter) *
                                       (double)(sum + 1 - nb_parameter)) / ((double)sum * ((double)(sum + nb_parameter) *
                                        (double)(sum + 1 - nb_parameter) + nb_parameter * (nb_state - nb_parameter)));
          }

          else {
            markov.transition[i][j] = ((double)nb_parameter) / ((double)(sum + nb_parameter) *
                                       (double)(sum + 1 - nb_parameter) + nb_parameter * (nb_state - nb_parameter));
          }
        }
        break;
      }

      // Ristad (1995) estimator 2

      case UNIFORM_CARDINALITY : {
        if (nb_parameter == nb_state) {
          for (j = 0;j < nb_state;j++) {
//            markov.transition[i][j] = (double)(transition[i][j] + 1) / (double)(sum + nb_state);
            markov.transition[i][j] = (double)transition[i][j] / (double)sum;
          }
        }

        else {
          for (j = 0;j < nb_state;j++) {
            if (transition[i][j] > 0) {
/*              markov.transition[i][j] = ((double)(transition[i][j] + 1) * (double)(sum + 1 - nb_parameter)) /
                                        ((double)(sum * sum + sum + 2 * nb_parameter)); */
              markov.transition[i][j] = ((double)transition[i][j] * ((double)sum * (double)(sum + 1) +
                                          nb_parameter * (1 - nb_parameter))) /
                                        ((double)sum * (double)(sum * sum + sum + 2 * nb_parameter));
            }
            else {
              markov.transition[i][j] = (double)(nb_parameter * (nb_parameter + 1)) /
                                        ((nb_state - nb_parameter) * (double)(sum * sum + sum + 2 * nb_parameter));
            }
          }
        }
        break;
      }
      }
    }

    else if (markov.memo_type[i] == COMPLETION) {
      for (j = 0;j < nb_state;j++) {
        markov.transition[i][j] = markov.transition[markov.parent[i]][j];
      }
    }

    else {
      for (j = 0;j < nb_state;j++) {
        markov.transition[i][j] = 0.;
      }
    }
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Estimation of a zero-order Markov chain.
 *
 *  \param[in] markov reference on a VariableOrderMarkov object.
 */
/*--------------------------------------------------------------*/

void VariableOrderMarkovData::order0_estimation(VariableOrderMarkov &markov) const

{
  int i , j;
//  int sum;


/*  sum = 0;
  for (i = 0;i < chain_data->nb_state;i++) {
    sum += chain_data->initial[i] + chain_data->transition[0][i];
  } */

  for (i = 0;i < chain_data->nb_state;i++) {
    markov.initial[i] = (double)marginal_distribution[0]->frequency[i] /
                        (double)marginal_distribution[0]->nb_element;
//    markov.initial[i] = (double)(chain_data->initial[i] + chain_data->transition[0][i]) /
//                        (double)sum;
    for (j = 0;j <= chain_data->nb_state;j++) {
      markov.transition[j][i] = markov.initial[i];
    }
  }

# ifdef DEBUG
  cout << "\n";
  for (i = 0;i < chain_data->nb_state;i++) {
    cout << STAT_label[STATL_STATE] << " " << i << ": " << markov.initial[i] << " | "
         << (double)(chain_data->initial[i] + chain_data->transition[0][i]) / (double)cumul_length << endl;
  }
# endif

}


/*--------------------------------------------------------------*/
/**
 *  \brief Estimation of a variable-order Markov chain.
 *
 *  \param[in] error                     reference on a StatError object,
 *  \param[in] display                   flag for displaying estimation intermediate results,
 *  \param[in] itype                     process type (ORDINARY/EQUILIBRIUM),
 *  \param[in] min_order                 minimum order of the variable-order Markov chain,
 *  \param[in] max_order                 maximum order of the variable-order Markov chain,
 *  \param[in] algorithm                 algorithm (CTM_BIC/CTM_KT/LOCAL_BIC/CONTEXT),
 *  \param[in] threshold                 threshold on the memory pruning,
 *  \param[in] estimator                 estimator (maximum likelihood, Laplace, adaptative Laplace),
 *  \param[in] global_initial_transition type of estimation of the initial transition probabilities (ordinary process case),
 *  \param[in] global_sample             type of management of the sample size
 *                                       (for the LOCAL_BIC or CONTEXT algorithm),
 *  \param[in] counting_flag             flag on the computation of the counting distributions.
 *
 *  \return                              VariableOrderMarkov object.
 */
/*--------------------------------------------------------------*/

VariableOrderMarkov* MarkovianSequences::variable_order_markov_estimation(StatError &error , bool display ,
                                                                          process_type itype , int min_order , int max_order ,
                                                                          memory_tree_selection algorithm , double threshold ,
                                                                          transition_estimator estimator , bool global_initial_transition ,
                                                                          bool global_sample , bool counting_flag) const

{
  bool status = true , order0 , *active_memory , *selected_memory;
  int i , j , k;
  int sample_size , length_nb_sequence , nb_row , state , nb_terminal , *memory_count ,
      *nb_parameter , *diff_nb_parameter;
  double num , denom , max_likelihood , *memory_likelihood , *diff_likelihood;
  VariableOrderMarkov *markov , *completed_markov;
  VariableOrderMarkovChainData *chain_data;
  VariableOrderMarkovData *seq;


  completed_markov = NULL;
  error.init();

  if ((type[0] != INT_VALUE) && (type[0] != STATE)) {
    status = false;
    ostringstream correction_message;
    correction_message << STAT_variable_word[INT_VALUE] << " or " << STAT_variable_word[STATE];
    error.correction_update(STAT_error[STATR_VARIABLE_TYPE] , (correction_message.str()).c_str());
  }

  else {
    if ((marginal_distribution[0]->nb_value < 2) ||
        (marginal_distribution[0]->nb_value > NB_STATE)) {
      status = false;
      error.update(SEQ_error[SEQR_NB_STATE]);
    }

    else if (!characteristics[0]) {
      for (i = 0;i < marginal_distribution[0]->nb_value;i++) {
        if (marginal_distribution[0]->frequency[i] == 0) {
          status = false;
          ostringstream error_message;
          error_message << SEQ_error[SEQR_MISSING_STATE] << " " << i;
          error.update((error_message.str()).c_str());
        }
      }
    }

    if ((min_order < 0) || (min_order >= max_order)) {
      status = false;
      error.update(SEQ_error[SEQR_MIN_ORDER]);
    }
    if ((max_order <= min_order) || (max_order > ORDER)) {
      status = false;
      error.update(SEQ_error[SEQR_MAX_ORDER]);
    }
    else {
      if ((int)pow((double)marginal_distribution[0]->nb_value , max_order + 1) > NB_PARAMETER) {
        status = false;
        error.update(SEQ_error[SEQR_NB_PARAMETER]);
      }
    }
  }

  if (nb_variable > 1) {
    if (nb_variable > 2) {
      status = false;
      error.correction_update(STAT_error[STATR_NB_VARIABLE] , "1 or 2");
    }

    if ((type[1] != INT_VALUE) && (type[1] != STATE)) {
      status = false;
      ostringstream error_message , correction_message;
      error_message << STAT_label[STATL_VARIABLE] << " " << 2 << ": "
                    << STAT_error[STATR_VARIABLE_TYPE];
      correction_message << STAT_variable_word[INT_VALUE] << " or "
                         << STAT_variable_word[STATE];
      error.correction_update((error_message.str()).c_str() , (correction_message.str()).c_str());
    }

    else {
      if (test_hidden(1)) {
        status = false;
        ostringstream error_message;
        error_message << STAT_label[STATL_VARIABLE] << " " << 2 << ": "
                      << SEQ_error[SEQR_OVERLAP];
        error.update((error_message.str()).c_str());
      }

      if (!characteristics[1]) {
        for (i = 0;i < marginal_distribution[1]->nb_value;i++) {
          if (marginal_distribution[1]->frequency[i] == 0) {
            status = false;
            ostringstream error_message;
            error_message << STAT_label[STATL_VARIABLE] << " " << 2 << ": "
                          << STAT_error[STATR_MISSING_VALUE] << " " << i;
            error.update((error_message.str()).c_str());
          }
        }
      }
    }
  }

  if (status) {
    if (display) {
      length_nb_sequence = nb_sequence;

      sample_size = cumul_length;
      cout << "\n" << STAT_label[STATL_SAMPLE_SIZE] << ":";
      for (i = 0;i <= MIN((int)::round(log((double)cumul_length) / log((double)marginal_distribution[0]->nb_value)) , max_length - 2);i++) {
        cout << " " << sample_size;
        sample_size -= length_nb_sequence;
        length_nb_sequence -= length_distribution->frequency[i + 1];
      }
      cout << endl;

      cout << SEQ_label[SEQL_RECOMMENDED_MAX_ORDER] << ": "
           << MIN((int)::round(log((double)cumul_length) / log((double)marginal_distribution[0]->nb_value)) , max_length - 2) << endl;

/*      if ((algorithm == CONTEXT) && (threshold == CONTEXT_THRESHOLD)) {
        Test test(CHI2);
        test.df1 = marginal_distribution[0]->nb_value - 1;
        test.critical_probability = 0.05;
        test.chi2_value_computation();

        threshold = test.value;

        cout << "\n" << SEQ_label[SEQL_PRUNING_THRESHOLD] << ": " << threshold << endl;
      } */
    }

    markov = new VariableOrderMarkov(itype , marginal_distribution[0]->nb_value , max_order , true);

    // counting of transitions and computation of the corresponding number of free parameters

    chain_data = new VariableOrderMarkovChainData(markov->type , markov->nb_state , markov->nb_row);
    transition_count_computation(*chain_data , *markov , false , true);
    chain_data->estimation(*markov , true , estimator);

    memory_count = new int[markov->nb_row];
    memory_likelihood = new double[markov->nb_row];
    nb_parameter = new int[markov->nb_row];
    diff_likelihood = new double[markov->nb_row];
    diff_nb_parameter = new int[markov->nb_row];

    for (i = 0;i < markov->nb_row;i++) {
      memory_count[i] = 0;
      for (j = 0;j < markov->nb_state;j++) {
        memory_count[i] += chain_data->transition[i][j];
      }
    }

    // BIC-type estimator

    if (algorithm != CTM_KT) {
      for (i = 0;i < markov->nb_row;i++) {
        memory_likelihood[i] = 0.;
        if (memory_count[i] > 0) {
          for (j = 0;j < markov->nb_state;j++) {
            if (chain_data->transition[i][j] > 0) {
              memory_likelihood[i] += chain_data->transition[i][j] * log(markov->transition[i][j]);
            }
          }
        }
      }
    }

    // Krichevsky-Trofimov estimator

    else {
      for (i = 0;i < markov->nb_row;i++) {
        memory_likelihood[i] = 0.;
        if (memory_count[i] > 0) {
          denom = memory_count[i] - 1. + (double)markov->nb_state / 2.;
          for (j = 0;j < markov->nb_state;j++) {
            if (chain_data->transition[i][j] > 0) {
              num = chain_data->transition[i][j] - 0.5;
              for (k = 0;k < chain_data->transition[i][j];k++) {
                memory_likelihood[i] += log(num) - log(denom);
                num--;
                denom--;
              }
            }
          }
        }
      }
    }

    for (i = 0;i < markov->nb_row;i++) {
      nb_parameter[i] = -1;
      for (j = 0;j < markov->nb_state;j++) {
        if (chain_data->transition[i][j] > 0) {
//        if (markov->transition[i][j] > 0.) {
          nb_parameter[i]++;
        }
      }

      if (nb_parameter[i] == -1) {
        nb_parameter[i] = 0;
      }
      else if ((algorithm == CTM_BIC) || (algorithm == LOCAL_BIC)) {
        memory_likelihood[i] = 2 * memory_likelihood[i] -
                               nb_parameter[i] * log((double)memory_count[0]);
      }
    }

#   ifdef DEBUG
    for (i = 0;i < markov->nb_row;i++) {
      for (j = markov->max_order - 1;j >= markov->order[i];j--) {
        cout << "  ";
      }
      for (j = markov->order[i] - 1;j >= 0;j--) {
        cout << markov->state[i][j] << " ";
      }
      cout << "   ";
      for (j = 0;j < markov->nb_state;j++) {
        cout << chain_data->transition[i][j] << "  ";
      }
      cout << " | ";
      for (j = 0;j < markov->nb_state;j++) {
        cout << markov->transition[i][j] << "  ";
      }
      cout << "   " << memory_likelihood[i] << "    " << nb_parameter[i] << endl;
    }
#   endif

    // pruning of the memory tree

    if (display) {
      if ((algorithm == CONTEXT) && (global_sample)) {
        cout << "\n" << SEQ_label[SEQL_PRUNING_THRESHOLD] << ": "
             << threshold * log((double)memory_count[0]) << endl;
      }

      cout << "\n";
    }

    if ((algorithm == CTM_BIC) || (algorithm == CTM_KT)) {
      active_memory = new bool[markov->nb_row];
      for (i = 0;i < markov->nb_row;i++) {
        active_memory[i] = false;
      }

      selected_memory = new bool[markov->nb_row];
      selected_memory[0] = true;
      for (i = 1;i < markov->nb_row;i++) {
        if (markov->order[i] <= min_order) {
          selected_memory[i] = true;
        }
        else {
          selected_memory[i] = false;
        }
      }

      // computation of the maximum BIC or Krichevsky-Trofimov estimator for each memory sub-tree

      for (i = markov->nb_row - 1;i >= 0;i--) {
        if ((markov->memo_type[i] == NON_TERMINAL) && (nb_parameter[i] > 0) &&
            (memory_count[i] >= MEMORY_MIN_COUNT)) {
          max_likelihood = 0.;
          diff_nb_parameter[i] = -nb_parameter[i];
          for (j = 0;j < markov->nb_state;j++) {
            max_likelihood += memory_likelihood[markov->child[i][j]];
            diff_nb_parameter[i] += nb_parameter[markov->child[i][j]];
          }

          diff_likelihood[i] = max_likelihood - memory_likelihood[i];

          if (max_likelihood > memory_likelihood[i]) {
            memory_likelihood[i] = max_likelihood;
            active_memory[i] = true;
          }
        }

        else if (display) {
          diff_likelihood[i] = 0.;
          diff_nb_parameter[i] = 0;
        }
      }

      // construction by chaining of the memory tree

      nb_row = 1;
      for (i = 0;i < markov->nb_row;i++) {
        if ((markov->order[i] >= min_order) && (active_memory[i]) && (selected_memory[i]) &&
            (diff_nb_parameter[i] >= 0) && (diff_likelihood[i] >= threshold)) {
//            (diff_likelihood[i] >= threshold)) {
          markov->memo_type[i] = NON_TERMINAL;
          for (j = 0;j < markov->nb_state;j++) {
            selected_memory[markov->child[i][j]] = true;
            markov->memo_type[markov->child[i][j]] = TERMINAL;
          }
          nb_row += markov->nb_state;
        }
      }

      if (nb_row == 1) {
        order0 = true;
        for (i = 0;i < markov->nb_state;i++) {
          selected_memory[markov->child[0][i]] = true;
          markov->memo_type[markov->child[0][i]] = TERMINAL;
        }
        nb_row += markov->nb_state;
      }

      else {
        order0 = false;
      }

      if (display) {
        for (i = 0;i < markov->nb_row;i++) {
          if ((nb_parameter[i] > 0) && (memory_count[i] >= MEMORY_MIN_COUNT)) {
            for (j = markov->max_order - 1;j >= markov->order[i];j--) {
              cout << "  ";
            }
            for (j = markov->order[i] - 1;j >= 0;j--) {
              cout << markov->state[i][j] << " ";
            }

            cout << "   " << diff_likelihood[i] << "   " << diff_nb_parameter[i] << " | " << memory_count[i]
                 << " | " << active_memory[i] << "   " << selected_memory[i] << endl;
          }
        }
      }

      for (i = markov->nb_row - 1;i >= 0;i--) {
        if (selected_memory[i]) {
          if ((markov->memo_type[i] == TERMINAL) && (markov->child[i])) {
            delete [] markov->child[i];
            markov->child[i] = NULL;
          }
        }

        else {
          markov->memo_type[i] = PRUNED;
          delete [] markov->state[i];
          markov->state[i] = NULL;
          if (markov->child[i]) {
            delete [] markov->child[i];
            markov->child[i] = NULL;
          }
        }
      }
    }

    else {
      order0 = false;
      nb_row = markov->nb_row;

      for (i = markov->nb_row - 1;i >= 0;i--) {
        if ((markov->memo_type[i] == NON_TERMINAL) && (markov->order[i] >= min_order)) {
          for (j = 0;j < markov->nb_state;j++) {
            if (markov->memo_type[markov->child[i][j]] != TERMINAL) {
              break;
            }
          }

          if (j == markov->nb_state) {
            if ((nb_parameter[i] > 0) && (memory_count[i] >= MEMORY_MIN_COUNT)) {
/*              if (algorithm == LOCAL_BIC) {
                if ((nb_parameter[i] > 0) && (memory_count[i] >= MEMORY_MIN_COUNT)) {
                  diff_likelihood[i] = -memory_likelihood[i];
                  diff_nb_parameter[i] = -nb_parameter[i];
                  for (j = 0;j < markov->nb_state;j++) {
                    diff_likelihood[i] += memory_likelihood[markov->child[i][j]];
                    diff_nb_parameter[i] += nb_parameter[markov->child[i][j]];
                  }
                }
              } */

              switch (algorithm) {
              case LOCAL_BIC :
                diff_likelihood[i] = 0.;
                diff_nb_parameter[i] = -nb_parameter[i];
                break;
              case CONTEXT :
                max_likelihood = D_INF;
                break;
              }

              for (j = 0;j < markov->nb_state;j++) {
                if (algorithm == CONTEXT) {
                  diff_likelihood[i] = 0.;
                }
                for (k = 0;k < markov->nb_state;k++) {
                  if (chain_data->transition[markov->child[i][j]][k] > 0) {
                    diff_likelihood[i] += chain_data->transition[markov->child[i][j]][k] *
                                          log(markov->transition[markov->child[i][j]][k] /
                                              markov->transition[i][k]);
                  }
                }

                if ((algorithm == CONTEXT) && (diff_likelihood[i] > max_likelihood)) {
                  max_likelihood = diff_likelihood[i];

                  if (display) {
                    state = j;
                  }
                }

                if (algorithm == LOCAL_BIC) {
                  diff_nb_parameter[i] += nb_parameter[markov->child[i][j]];
                }
              }

              if (algorithm == LOCAL_BIC) {
                diff_likelihood[i] = 2 * diff_likelihood[i] - diff_nb_parameter[i] *
                                     log((double)memory_count[global_sample ? 0 : i]);
              }

//              if ((display) && (diff_likelihood >= threshold)) {
              if (display) {
                for (j = markov->max_order - 1;j >= markov->order[i];j--) {
                  cout << "  ";
                }
                for (j = markov->order[i] - 1;j >= 0;j--) {
                  cout << markov->state[i][j] << " ";
                }

                switch (algorithm) {

                case LOCAL_BIC : {
                  cout << "   " << diff_likelihood[i] << "   " << diff_nb_parameter[i]
                       << "   " << memory_count[i] << endl;
                  break;
                }

                case CONTEXT : {
                  cout << "   " << 2 * max_likelihood << "   " << state;
                  if (!global_sample) {
                    cout << "   " << threshold * log((double)memory_count[i]);
                  }
                  cout << endl;
                  break;
                }
                }
              }
            }

            if ((nb_parameter[i] == 0) || (memory_count[i] < MEMORY_MIN_COUNT) ||
                ((algorithm == LOCAL_BIC) && ((diff_nb_parameter[i] < 0) || (diff_likelihood[i] < threshold))) ||
                ((algorithm == CONTEXT) &&
                 (2 * max_likelihood < threshold * log((double)memory_count[global_sample ? 0 : i])))) {
//                 (2 * max_likelihood < threshold))) {
              if (i > 0) {
                markov->memo_type[i] = TERMINAL;

                for (j = 0;j < markov->nb_state;j++) {
                  markov->memo_type[markov->child[i][j]] = PRUNED;
                  delete [] markov->state[markov->child[i][j]];
                  markov->state[markov->child[i][j]] = NULL;
                }
                delete [] markov->child[i];
                markov->child[i] = NULL;

                nb_row -= markov->nb_state;
              }

              else {
                order0 = true;
              }
            }
          }
        }
      }
    }

    // copy of the conserved memories

    i = 1;
    for (j = 1;j < markov->nb_row;j++) {
      if (markov->memo_type[j] != PRUNED) {
        if (i != j) {
          for (k = 0;k < markov->nb_state;k++) {
            markov->transition[i][k] = markov->transition[j][k];
          }

          markov->memo_type[i] = markov->memo_type[j];
          markov->order[i] = markov->order[j];

          delete [] markov->state[i];
          markov->state[i] = new int[markov->order[i]];
          for (k = 0;k < markov->order[i];k++) {
            markov->state[i][k] = markov->state[j][k];
          }

          for (k = i - 1;k >= 0;k--) {
            if ((markov->child[k]) &&
                (markov->child[k][markov->state[i][markov->order[i] - 1]] == j)) {
              markov->child[k][markov->state[i][markov->order[i] - 1]] = i;
              markov->parent[i] = k;
              break;
            }
          }

          if (markov->child[j]) {
            if (!markov->child[i]) {
              markov->child[i] = new int[markov->nb_state];
            }
            for (k = 0;k < markov->nb_state;k++) {
              markov->child[i][k] = markov->child[j][k];
            }
          }

          else if (markov->child[i]) {
            delete [] markov->child[i];
            markov->child[i] = NULL;
          }
        }

        i++;
      }
    }

    for (i = 1;i < markov->nb_row;i++) {
      delete [] markov->next[i];
    }
    delete [] markov->next;
    markov->next = NULL;

    markov->nb_row = nb_row;
    markov->max_order_computation();

    delete chain_data;

    delete [] memory_count;
    delete [] memory_likelihood;
    delete [] nb_parameter;
    delete [] diff_likelihood;
    delete [] diff_nb_parameter;

    if ((algorithm == CTM_BIC) || (algorithm == CTM_KT)) {
      delete [] active_memory;
      delete [] selected_memory;
    }

    completed_markov = new VariableOrderMarkov(*markov , nb_variable - 1 ,
                                               (nb_variable == 2 ? marginal_distribution[1]->nb_value : 0));
    delete markov;

    completed_markov->markov_data = new VariableOrderMarkovData(*this , SEQUENCE_COPY ,
                                                                (completed_markov->type == EQUILIBRIUM ? true : false));

    seq = completed_markov->markov_data;
    seq->state_variable_init();

    if (order0) {
      seq->build_transition_count(*completed_markov , true , true);
      seq->order0_estimation(*completed_markov);
    }

    else {
      seq->build_transition_count(*completed_markov , true ,
                                  (((completed_markov->type == ORDINARY) && (global_initial_transition)) ? true : false));
      seq->chain_data->estimation(*completed_markov);

      if ((completed_markov->type == ORDINARY) && (global_initial_transition)) {
        for (i = 1;i < completed_markov->nb_row;i++) {
          if (completed_markov->memo_type[i] == NON_TERMINAL) {
            for (j = 0;j < completed_markov->nb_state;j++) {
              for (k = 0;k < completed_markov->nb_state;k++) {
                seq->chain_data->transition[i][k] -= seq->chain_data->transition[completed_markov->child[i][j]][k];
              }
            }
          }
        }
      }
    }

    if (completed_markov->type == EQUILIBRIUM) {
      nb_terminal = (completed_markov->nb_row - 1) * (completed_markov->nb_state - 1) /
                    completed_markov->nb_state + 1;

      for (i = 1;i < completed_markov->nb_row;i++) {
        if (!completed_markov->child[i]) {
          completed_markov->initial[i] = 1. / (double)nb_terminal;
        }
        else {
          completed_markov->initial[i] = 0.;
        }
      }

      completed_markov->initial_probability_computation();
    }

    // estimation of the categorical observation distributions

    if (completed_markov->nb_output_process == 1) {
      seq->build_observation_frequency_distribution(completed_markov->nb_state);

      for (i = 0;i < completed_markov->nb_state;i++) {
        seq->observation_distribution[1][i]->distribution_estimation(completed_markov->categorical_process[0]->observation[i]);
      }
    }

#   ifdef DEBUG
    for (i = 1;i < completed_markov->nb_row;i++) {
      for (j = completed_markov->max_order - 1;j >= completed_markov->order[i];j--) {
        cout << "  ";
      }
      for (j = completed_markov->order[i] - 1;j >= 0;j--) {
        cout << completed_markov->state[i][j] << " ";
      }
      cout << "   ";
      for (j = 0;j < completed_markov->nb_state;j++) {
        cout << seq->chain_data->transition[i][j] << "  ";
      }
      cout << " | ";
      for (j = 0;j < completed_markov->nb_state;j++) {
        cout << completed_markov->transition[i][j] << "  ";
      }
      cout << endl;
    }
#   endif

    // computation of the log-likelihood and the characteristic distributions of the model

    seq->likelihood = completed_markov->likelihood_computation(*seq);

    if (display) {
      cout << "\n" << STAT_label[STATL_LIKELIHOOD] << ": " << seq->likelihood
           << " | " << completed_markov->likelihood_computation(*seq , I_DEFAULT) << endl;
    }

    if (seq->likelihood == D_INF) {
      delete completed_markov;
      completed_markov = NULL;
      error.update(STAT_error[STATR_ESTIMATION_FAILURE]);
    }

    else {
      completed_markov->component_computation();
      completed_markov->characteristic_computation(*seq , counting_flag , I_DEFAULT , false);
    }
  }

  return completed_markov;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Estimation of a variable-order Markov chain.
 *
 *  \param[in] error                     reference on a StatError object,
 *  \param[in] imarkov                   reference on a VariableOrderMarkov object,
 *  \param[in] global_initial_transition type of estimation of the initial transition probabilities (ordinary process case),
 *  \param[in] counting_flag             flag on the computation of the counting distributions.
 *
 *  \return                              VariableOrderMarkov object.
 */
/*--------------------------------------------------------------*/

VariableOrderMarkov* MarkovianSequences::variable_order_markov_estimation(StatError &error ,
                                                                          const VariableOrderMarkov &imarkov,
                                                                          bool global_initial_transition ,
                                                                          bool counting_flag) const

{
  bool status = true;
  int i , j , k;
  int nb_terminal;
  VariableOrderMarkov *markov;
  VariableOrderMarkovData *seq;


  markov = NULL;
  error.init();

  if ((type[0] != INT_VALUE) && (type[0] != STATE)) {
    status = false;
    ostringstream correction_message;
    correction_message << STAT_variable_word[INT_VALUE] << " or " << STAT_variable_word[STATE];
    error.correction_update(STAT_error[STATR_VARIABLE_TYPE] , (correction_message.str()).c_str());
  }

  else {
    if ((marginal_distribution[0]->nb_value < 2) ||
        (marginal_distribution[0]->nb_value > NB_STATE)) {
      status = false;
      error.update(SEQ_error[SEQR_NB_STATE]);
    }

    else if (!characteristics[0]) {
      for (i = 0;i < marginal_distribution[0]->nb_value;i++) {
        if (marginal_distribution[0]->frequency[i] == 0) {
          status = false;
          ostringstream error_message;
          error_message << SEQ_error[SEQR_MISSING_STATE] << " " << i;
          error.update((error_message.str()).c_str());
        }
      }
    }
  }

  if (nb_variable > 1) {
    if (nb_variable > 2) {
      status = false;
      error.correction_update(STAT_error[STATR_NB_VARIABLE] , "1 or 2");
    }

    if ((type[1] != INT_VALUE) && (type[1] != STATE)) {
      status = false;
      ostringstream error_message , correction_message;
      error_message << STAT_label[STATL_VARIABLE] << " " << 2 << ": "
                    << STAT_error[STATR_VARIABLE_TYPE];
      correction_message << STAT_variable_word[INT_VALUE] << " or "
                         << STAT_variable_word[STATE];
      error.correction_update((error_message.str()).c_str() , (correction_message.str()).c_str());
    }

    else {
      if (test_hidden(1)) {
        status = false;
        ostringstream error_message;
        error_message << STAT_label[STATL_VARIABLE] << " " << 2 << ": "
                      << SEQ_error[SEQR_OVERLAP];
        error.update((error_message.str()).c_str());
      }

      if (marginal_distribution[1]->nb_value > NB_STATE) {
        status = false;
        error.update(STAT_error[STATR_NB_OUTPUT]);
      }

/*      if (!characteristics[1]) {
        for (i = 0;i < marginal_distribution[1]->nb_value;i++) {
          if (marginal_distribution[1]->frequency[i] == 0) {
            status = false;
            ostringstream error_message;
            error_message << STAT_label[STATL_VARIABLE] << " " << 2 << ": "
                          << STAT_error[STATR_MISSING_VALUE] << " " << i;
            error.update((error_message.str()).c_str());
          }
        }
      } */
    }
  }

  if (imarkov.nb_output_process != nb_variable - 1) {
    status = false;
    error.update(STAT_error[STATR_NB_OUTPUT_PROCESS]);
  }

  if (status) {
    markov = new VariableOrderMarkov(imarkov , false);

    markov->markov_data = new VariableOrderMarkovData(*this , SEQUENCE_COPY , (markov->type == EQUILIBRIUM ? true : false));

    seq = markov->markov_data;
    seq->state_variable_init();
    seq->build_transition_count(*markov , true ,
                               (((markov->type == ORDINARY) && (global_initial_transition)) ? true : false));
    seq->chain_data->estimation(*markov);

    if ((markov->type == ORDINARY) && (global_initial_transition)) {
      for (i = 1;i < markov->nb_row;i++) {
        if (markov->memo_type[i] == NON_TERMINAL) {
          for (j = 0;j < markov->nb_state;j++) {
            for (k = 0;k < markov->nb_state;k++) {
              seq->chain_data->transition[i][k] -= seq->chain_data->transition[markov->child[i][j]][k];
            }
          }
        }
      }
    }

    if (markov->type == EQUILIBRIUM) {
      nb_terminal = (markov->nb_row - 1) * (markov->nb_state - 1) / markov->nb_state + 1;

      for (i = 1;i < markov->nb_row;i++) {
        if (!markov->child[i]) {
          markov->initial[i] = 1. / (double)nb_terminal;
        }
        else {
          markov->initial[i] = 0.;
        }
      }

      markov->initial_probability_computation();
    }

    // estimation of the categorical observation distributions

    if (markov->nb_output_process == 1) {
      seq->build_observation_frequency_distribution(markov->nb_state);

      for (i = 0;i < markov->nb_state;i++) {
        seq->observation_distribution[1][i]->distribution_estimation(markov->categorical_process[0]->observation[i]);
      }
    }

    // computation of the log-likelihood and the characteristic distributions of the model

    seq->likelihood = markov->likelihood_computation(*seq);

#   ifdef MESSAGE
    cout << "\n" << STAT_label[STATL_LIKELIHOOD] << ": " << seq->likelihood
         << " | " << markov->likelihood_computation(*seq , I_DEFAULT) << endl;
#   endif

    if (seq->likelihood == D_INF) {
      delete markov;
      markov = NULL;
      error.update(STAT_error[STATR_ESTIMATION_FAILURE]);
    }

    else {
      markov->component_computation();
      markov->characteristic_computation(*seq , counting_flag , I_DEFAULT , false);
    }
  }

  return markov;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Estimation of a fixed-order  Markov chain.
 *
 *  \param[in] error                     reference on a StatError object,
 *  \param[in] type                      process type (ORDINARY/EQUILIBRIUM),
 *  \param[in] order                     Markov chain order,
 *  \param[in] global_initial_transition type of estimation of the initial transition probabilities (ordinary process case),
 *  \param[in] counting_flag             flag on the computation of the counting distributions.
 *
 *  \return                              VariableOrderMarkov object.
 */
/*--------------------------------------------------------------*/

VariableOrderMarkov* MarkovianSequences::variable_order_markov_estimation(StatError &error ,
                                                                          process_type type , int order ,
                                                                          bool global_initial_transition ,
                                                                          bool counting_flag) const

{
  bool status = true;
  VariableOrderMarkov *imarkov , *markov;


  markov = NULL;
  error.init();

  if ((order < 1) || (order > ORDER)) {
    status = false;
    error.update(SEQ_error[SEQR_ORDER]);
  }
  else {
    if ((int)pow((double)marginal_distribution[0]->nb_value , order + 1) > NB_PARAMETER) {
      status = false;
      error.update(SEQ_error[SEQR_NB_PARAMETER]);
    }
  }

  if (status) {
    imarkov = new VariableOrderMarkov(type , marginal_distribution[0]->nb_value , order , true ,
                                      nb_variable - 1 , (nb_variable == 2 ? marginal_distribution[1]->nb_value : 0));
    imarkov->build_previous_memory();

    markov = variable_order_markov_estimation(error , *imarkov , global_initial_transition ,
                                              counting_flag);
    delete imarkov;
  }

  return markov;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of transition counts.
 *
 *  \param[in,out] os    stream,
 *  \param[in]     begin flag for taking account of the beginning of sequences.
 */
/*--------------------------------------------------------------*/

ostream& VariableOrderMarkov::transition_count_ascii_write(ostream &os , bool begin) const

{
  bool *bic_memory , *kt_memory;
  int i , j , k;
  int buff , max_memory_count , row , initial_count , *memory_count , *max_state ,
      *nb_parameter , *diff_nb_parameter , width[3];
  double standard_normal_value , half_confidence_interval , diff , max_abs_diff , child_likelihood ,
         child_krichevsky_trofimov , num , denom , *diff_count , *initial_likelihood ,
         *memory_likelihood , *max_likelihood , *krichevsky_trofimov , *max_krichevsky_trofimov ,
         **confidence_limit , **transition_likelihood , **diff_likelihood;
  ios_base::fmtflags format_flags;


  format_flags = os.setf(ios::right , ios::adjustfield);

  memory_count = new int[nb_row];

  if (begin) {
    initial_count = 0;
    for (i = 0;i < nb_state;i++) {
      initial_count += markov_data->chain_data->initial[i];
    }
    width[0] = column_width(initial_count);

    max_memory_count = 0;
  }

  else {
    width[0] = 0;
    row = 0;
  }

  for (i = (begin ? 1 : 0);i < nb_row;i++) {
    memory_count[i] = 0;
    for (j = 0;j < nb_state;j++) {
      memory_count[i] += markov_data->chain_data->transition[i][j];
    }

    if ((begin) && (memory_count[i] > max_memory_count)) {
      max_memory_count = memory_count[i];
      row = i;
    }
  }

  buff = column_width(memory_count[row]);
  if (buff > width[0]) {
    width[0] = buff;
  }
  width[0] += ASCII_SPACE;

  // Peres-Shields fluctuation estimator

  diff_count = new double[nb_row];
  max_state = new int[nb_row];

  diff_count[0] = memory_count[row];

  for (i = 1;i < nb_row;i++) {
    diff_count[i] = memory_count[row];
    if (((!begin) || (order[i] > 1)) && (memory_count[i] > 0)) {
      max_abs_diff = 0.;
      for (j = 0;j < nb_state;j++) {
        diff = fabs(transition[parent[i]][j] * memory_count[i] -
                    markov_data->chain_data->transition[i][j]);
//        if (diff > diff_count[i]) {
//          diff_count[i] = diff;
        if (diff > max_abs_diff) {
          max_abs_diff = diff;
          diff_count[i] = transition[parent[i]][j] * memory_count[i] -
                          markov_data->chain_data->transition[i][j];
          max_state[i] = j;
        }
      }
    }
  }
  width[1] = column_width(nb_row - 1 , diff_count + 1) + ASCII_SPACE;

  if (begin) {
    os << "\n" << SEQ_label[SEQL_INITIAL_COUNTS] << endl;

    for (i = 0;i < max_order;i++) {
      os << "  ";
    }
    os << "  ";

    for (i = 0;i < nb_state;i++) {
      os << setw(width[0]) << i;
    }
    os << endl;

    for (i = 0;i < max_order;i++) {
      os << "  ";
    }
    os << "  ";

    for (i = 0;i < nb_state;i++) {
      os << setw(width[0]) << markov_data->chain_data->initial[i];
    }
    os << "  " << setw(width[0]) << initial_count << endl;
  }

//  os << "\nthreshold: " << pow((double)memory_count[0] , 0.75) << endl;

  os << "\n" << SEQ_label[SEQL_TRANSITION_COUNTS] << endl;

  for (i = 0;i < max_order;i++) {
    os << "  ";
  }
  os << "  ";

  for (i = 0;i < nb_state;i++) {
    os << setw(width[0]) << i;
  }
  os << setw(width[0]) << " "
     << "    " << SEQ_label[SEQL_MAX_TRANSITION_COUNT_DIFFERENCE];

  for (i = (begin ? 1 : 0);i <= max_order;i++) {
    os << "\n";
    for (j = 0;j < nb_row;j++) {
      if (order[j] == i) {
        for (k = max_order - 1;k >= order[j];k--) {
          os << "  ";
        }
        for (k = order[j] - 1;k >= 0;k--) {
          os << state[j][k] << " ";
        }
        os << "  ";

        for (k = 0;k < nb_state;k++) {
          os << setw(width[0]) << markov_data->chain_data->transition[j][k];
        }
        os << "  " << setw(width[0]) << memory_count[j];

        if (diff_count[j] != memory_count[row]) {
          os << setw(width[1]) << diff_count[j] << " (" << max_state[j] << ")";
        }
        os << endl;
      }
    }
  }
  os << endl;

  // computation of column widths

  width[1] = (begin ? column_width(nb_state , initial) : 0);

  for (i = (begin ? 1 : 0);i < nb_row;i++) {
    buff = column_width(nb_state , transition[i]);
    if (buff > width[1]) {
      width[1] = buff;
    }
  }
  width[1] += ASCII_SPACE;

  if (begin) {
    os << "\n" << STAT_word[STATW_INITIAL_PROBABILITIES] << endl;

    for (i = 0;i < max_order;i++) {
      os << "  ";
    }
    os << "  ";

    for (i = 0;i < nb_state;i++) {
      os << setw(width[1]) << i;
    }
    os << endl;

    for (i = 0;i < max_order;i++) {
      os << "  ";
    }
    os << "  ";

    for (i = 0;i < nb_state;i++) {
      os << setw(width[1]) << initial[i];
    }
    os << endl;
  }

  os << "\n" << STAT_word[STATW_TRANSITION_PROBABILITIES] << endl;

  for (i = 0;i < max_order;i++) {
    os << "  ";
  }
  os << "  ";

  for (i = 0;i < nb_state;i++) {
    os << setw(width[1]) << i;
  }

  for (i = (begin ? 1 : 0);i <= max_order;i++) {
    os << "\n";
    for (j = 0;j < nb_row;j++) {
      if ((order[j] == i) && (memory_count[j] > 0)) {
        for (k = max_order - 1;k >= order[j];k--) {
          os << "  ";
        }
        for (k = order[j] - 1;k >= 0;k--) {
          os << state[j][k] << " ";
        }
        os << "  ";

        for (k = 0;k < nb_state;k++) {
          os << setw(width[1]) << transition[j][k];
        }
        os << endl;
      }
    }
  }

  // computation of confidence intervals on the transition probabilities

  confidence_limit = new double*[nb_row];
  for (i = (begin ? 1 : 0);i < nb_row;i++) {
    if (memory_count[i] > 0) {
      confidence_limit[i] = new double[nb_state * 2];
      for (j = 0;j < nb_state * 2;j++) {
        confidence_limit[i][j] = 0.;
      }
    }
  }

  normal dist;
  standard_normal_value = quantile(complement(dist , 0.025));

  for (i = (begin ? 1 : 0);i < nb_row;i++) {
    if (memory_count[i] > 0) {
      for (j = 0;j < nb_state;j++) {
        if ((transition[i][j] > 0.) && (transition[i][j] < 1.)) {
          half_confidence_interval = standard_normal_value *
                                     sqrt(transition[i][j] * (1. - transition[i][j]) / memory_count[i]);
          confidence_limit[i][2 * j] = MAX(transition[i][j] - half_confidence_interval , 0.);
          confidence_limit[i][2 * j + 1] = MIN(transition[i][j] + half_confidence_interval , 1.);
        }
      }
    }
  }

  // computation of column widths

  width[1] = 0;
  for (i = (begin ? 1 : 0);i < nb_row;i++) {
    if (memory_count[i] > 0) {
      buff = column_width(2 * nb_state , confidence_limit[i]);
      if (buff > width[1]) {
        width[1] = buff;
      }
    }
  }
  width[1] += ASCII_SPACE;

  os << "\n" << SEQ_label[SEQL_TRANSITION_PROBABILITIY_CONFIDENCE_INTERVAL] << endl;

  for (i = 0;i < max_order;i++) {
    os << "  ";
  }
  os << "  ";

  for (i = 0;i < nb_state;i++) {
    os << setw(width[1]) << i
       << setw(width[1]) << " ";
  }
  os << "  " << SEQ_label[SEQL_COUNT];

  for (i = (begin ? 1 : 0);i <= max_order;i++) {
    os << "\n";
    for (j = 0;j < nb_row;j++) {
      if ((order[j] == i) && (memory_count[j] > 0)) {
        for (k = max_order - 1;k >= order[j];k--) {
          os << "  ";
        }
        for (k = order[j] - 1;k >= 0;k--) {
          os << state[j][k] << " ";
        }
        os << "  ";

        for (k = 0;k < nb_state;k++) {
          if ((transition[j][k] > 0.) && (transition[j][k] < 1.)) {
            os << setw(width[1]) << confidence_limit[j][2 * k]
               << setw(width[1]) << confidence_limit[j][2 * k + 1];
          }
          else {
            os << setw(width[1]) << " "
               << setw(width[1]) << " ";
          }
        }
        os << "  " << setw(width[0]) << memory_count[j] << endl;
      }
    }
  }

  // computation of log-likelihoods

  if (begin) {
    initial_likelihood = new double[nb_state + 1];

    initial_likelihood[nb_state] = 0.;
    for (i = 0;i < nb_state;i++) {
      if (markov_data->chain_data->initial[i] > 0) {
        initial_likelihood[i] = markov_data->chain_data->initial[i] * log(initial[i]);
        initial_likelihood[nb_state] += initial_likelihood[i];
      }
      else {
        initial_likelihood[i] = 0.;
      }
    }
  }

  transition_likelihood = new double*[nb_row];
  memory_likelihood = new double[nb_row];
  krichevsky_trofimov = new double[nb_row];
  nb_parameter = new int[nb_row];

  diff_likelihood = new double*[2];
  diff_likelihood[0] = new double[nb_row];
  diff_likelihood[1] = new double[nb_row];

  max_likelihood = new double[nb_row];
  bic_memory = new bool[nb_row];
  max_krichevsky_trofimov = new double[nb_row];
  kt_memory = new bool[nb_row];
  diff_nb_parameter = new int[nb_row];

  if (begin) {
    memory_likelihood[0] = 0.;
    krichevsky_trofimov[0] = 0.;
    nb_parameter[0] = 0;

    diff_likelihood[0][0] = 0;
    diff_likelihood[1][0] = 0;
    max_likelihood[0] = 0;
    max_krichevsky_trofimov[0] = 0;
  }

  for (i = (begin ? 1 : 0);i < nb_row;i++) {
    transition_likelihood[i] = new double[nb_state];

    memory_likelihood[i] = 0.;
    nb_parameter[i] = 0;
    for (j = 0;j < nb_state;j++) {
      if (markov_data->chain_data->transition[i][j] > 0) {
        nb_parameter[i]++;
        transition_likelihood[i][j] = markov_data->chain_data->transition[i][j] * log(transition[i][j]);
        memory_likelihood[i] += transition_likelihood[i][j];
      }
      else {
        transition_likelihood[i][j] = 0.;
      }
    }

    if (nb_parameter[i] > 0) {
      nb_parameter[i]--;
    }

    // Krichevsky-Trofimov estimator

    krichevsky_trofimov[i] = 0.;
//    krichevsky_trofimov[i] = 1.;
    if (memory_count[i] > 0) {
      denom = memory_count[i] - 1. + (double)nb_state / 2.;
      for (j = 0;j < nb_state;j++) {
        if (markov_data->chain_data->transition[i][j] > 0) {
          num = markov_data->chain_data->transition[i][j] - 0.5;
          for (k = 0;k < markov_data->chain_data->transition[i][j];k++) {
//          krichevsky_trofimov[i] *= num / denom;
            krichevsky_trofimov[i] += log(num) - log(denom);
            num--;
            denom--;
          }
        }
      }
    }
  }

  for (i = nb_row - 1;i >= (begin ? 1 : 0);i--) {
    max_likelihood[i] = 2 * memory_likelihood[i] - nb_parameter[i] * log((double)memory_count[0]);
    max_krichevsky_trofimov[i] = krichevsky_trofimov[i];

    if (order[i] == max_order) {
      diff_likelihood[0][i] = 0.;
      diff_likelihood[1][i] = 0.;
    }

    else {
      diff = -memory_likelihood[i];
      child_likelihood = 0.;
      child_krichevsky_trofimov = 0.;
      diff_nb_parameter[i] = -nb_parameter[i];
      for (j = 0;j < nb_state;j++) {
        diff += memory_likelihood[child[i][j]];
        child_likelihood += max_likelihood[child[i][j]];
        child_krichevsky_trofimov += max_krichevsky_trofimov[child[i][j]];
        diff_nb_parameter[i] += nb_parameter[child[i][j]];
      }
      diff_likelihood[0][i] = 2 * diff - diff_nb_parameter[i] * log((double)memory_count[0]);
      diff_likelihood[1][i] = 2 * diff - diff_nb_parameter[i] * log((double)memory_count[i]);

      if (child_likelihood > max_likelihood[i]) {
        max_likelihood[i] = child_likelihood;
        bic_memory[i] = true;
      }
      else {
        bic_memory[i] = false;
      }

      if (child_krichevsky_trofimov > max_krichevsky_trofimov[i]) {
        max_krichevsky_trofimov[i] = child_krichevsky_trofimov;
        kt_memory[i] = true;
      }
      else {
        kt_memory[i] = false;
      }
    }
  }

  // computation of column widths

  if (begin) {
    width[1] = column_width(nb_state + 1 , initial_likelihood);
  }
  else {
    width[1] = 0;
  }

  for (i = (begin ? 1 : 0);i < nb_row;i++) {
    buff = column_width(nb_state , transition_likelihood[i]);
    if (buff > width[1]) {
      width[1] = buff;
    }
  }

  buff = column_width(nb_row , memory_likelihood);
  if (buff > width[1]) {
    width[1] = buff;
  }
  buff = column_width(nb_row , krichevsky_trofimov);
  if (buff > width[1]) {
    width[1] = buff;
  }

  buff = column_width(nb_row , diff_likelihood[0]);
  if (buff > width[1]) {
    width[1] = buff;
  }
  buff = column_width(nb_row , diff_likelihood[1]);
  if (buff > width[1]) {
    width[1] = buff;
  }
  buff = column_width(nb_row , max_likelihood);
  if (buff > width[1]) {
    width[1] = buff;
  }
  buff = column_width(nb_row , max_krichevsky_trofimov);
  if (buff > width[1]) {
    width[1] = buff;
  }
  width[1] += ASCII_SPACE;

  width[2] = column_width(nb_state - 1) + ASCII_SPACE;

  os << "\n\n" << SEQ_label[SEQL_LIKELIHOODS] << endl;

  if (begin) {
    for (i = 0;i < max_order;i++) {
      os << "  ";
    }
    os << "  ";

    for (i = 0;i < nb_state;i++) {
      os << setw(width[1]) << i;
    }
    os << endl;

    for (i = 0;i < max_order;i++) {
      os << "  ";
    }
    os << "  ";

    for (i = 0;i < nb_state;i++) {
      os << setw(width[1]) << initial_likelihood[i];
    }
    os << "  " << setw(width[1]) << initial_likelihood[nb_state] << "\n\n";
  }

  for (i = 0;i < max_order;i++) {
    os << "  ";
  }
  os << "  ";

  for (i = 0;i < nb_state;i++) {
    os << setw(width[1]) << i;
  }
  os << setw(width[1]) << " "
     << "    " << SEQ_label[SEQL_COUNT]
     << "    " << STAT_label[STATL_FREE_PARAMETERS];
  if (!begin) {
/*    os << "    " << STAT_criterion_word[BIC] << "    " << STAT_criterion_word[BICc]
       << "    " << SEQ_label[SEQL_KRICHEVSKY_TROFIMOV]; */

    os << "    " << SEQ_label[SEQL_DELTA] << " " << STAT_label[STATL_FREE_PARAMETERS]
       << "    " << SEQ_label[SEQL_DELTA] << " " << STAT_criterion_word[BIC]
       << "    " << SEQ_label[SEQL_DELTA] << " " << STAT_criterion_word[BICc]
       << "    " << SEQ_label[SEQL_DELTA] << " " << STAT_criterion_word[BIC]
       << "    " << SEQ_label[SEQL_DELTA] << " " << SEQ_label[SEQL_KRICHEVSKY_TROFIMOV];
  }
  os << endl;

  for (i = (begin ? 1 : 0);i <= max_order;i++) {
    os << "\n";
    for (j = 0;j < nb_row;j++) {
      if (order[j] == i) {
        for (k = max_order - 1;k >= order[j];k--) {
          os << "  ";
        }
        for (k = order[j] - 1;k >= 0;k--) {
          os << state[j][k] << " ";
        }
        os << "  ";

        for (k = 0;k < nb_state;k++) {
          os << setw(width[1]) << transition_likelihood[j][k];
        }

        os << "  " << setw(width[1]) << memory_likelihood[j]
           << "  " << setw(width[0]) << memory_count[j]
           << "  " << setw(width[2]) << nb_parameter[j];

/*        if (!begin) {
          os << "  " << setw(width[1]) << 2 * memory_likelihood[j] -
             nb_parameter[j] * log((double)memory_count[0]);

          if (memory_count[parent[j]] > 0) {
            os << "  " << setw(width[1]) << 2 * memory_likelihood[j] -
               nb_parameter[j] * log((double)memory_count[parent[j]]);
          }
          else {
            os << "  " << setw(width[1]) << 0;
          }

          os << "  " << setw(width[1]) << krichevsky_trofimov[j];
        } */

        if ((!begin) && (order[j] < max_order)) {
          os << "  " << setw(width[2]) << diff_nb_parameter[j]
             << "  " << setw(width[1]) << diff_likelihood[0][j];
          if (memory_count[j] > 0) {
            os << setw(width[1]) << diff_likelihood[1][j];
          }
          else {
            os << setw(width[1]) << 0;
          }
          os << "  " << setw(width[1]) << max_likelihood[j] << "  " << bic_memory[j]
             << setw(width[1]) << max_krichevsky_trofimov[j] << "  " << kt_memory[j];
        }

        os << endl;
      }
    }
  }

  delete [] diff_count;
  delete [] max_state;

  for (i = (begin ? 1 : 0);i < nb_row;i++) {
    if (memory_count[i] > 0) {
      delete [] confidence_limit[i];
    }
  }
  delete [] confidence_limit;

  delete [] memory_count;

  if (begin) {
    delete [] initial_likelihood;
  }

  for (i = (begin ? 1 : 0);i < nb_row;i++) {
    delete [] transition_likelihood[i];
  }
  delete [] transition_likelihood;

  delete [] memory_likelihood;
  delete [] krichevsky_trofimov;
  delete [] nb_parameter;

  delete [] diff_likelihood[0];
  delete [] diff_likelihood[1];
  delete [] diff_likelihood;

  delete [] max_likelihood;
  delete [] bic_memory;
  delete [] max_krichevsky_trofimov;
  delete [] kt_memory;
  delete [] diff_nb_parameter;

  os.setf(format_flags , ios::adjustfield);

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Counting of transitions for successive orders.
 *
 *  \param[in] error     reference on a StatError object,
 *  \param[in] display   flag for displaying transition counts,
 *  \param[in] max_order maximum order,
 *  \param[in] begin     flag for taking account of the beginning of sequences,
 *  \param[in] estimator estimator (maximum likelihood, Laplace, adaptative Laplace),
 *  \param[in] path      file path.
 *
 *  \return              error status.
 */
/*--------------------------------------------------------------*/

bool MarkovianSequences::transition_count(StatError &error , bool display , int max_order ,
                                          bool begin , transition_estimator estimator ,
                                          const string path) const

{
  bool status = true;
  int i;
  VariableOrderMarkov *markov;
  VariableOrderMarkovData *seq;


  error.init();

  if (nb_variable > 1) {
    status = false;
    error.correction_update(STAT_error[STATR_NB_VARIABLE] , 1);
  }

  if ((type[0] != INT_VALUE) && (type[0] != STATE)) {
    status = false;
    ostringstream correction_message;
    correction_message << STAT_variable_word[INT_VALUE] << " or " << STAT_variable_word[STATE];
    error.correction_update(STAT_error[STATR_VARIABLE_TYPE] , (correction_message.str()).c_str());
  }

  else {
    if ((marginal_distribution[0]->nb_value < 2) ||
        (marginal_distribution[0]->nb_value > NB_STATE)) {
      status = false;
      error.update(SEQ_error[SEQR_NB_STATE]);
    }

    else if (!characteristics[0]) {
      for (i = 0;i < marginal_distribution[0]->nb_value;i++) {
        if (marginal_distribution[0]->frequency[i] == 0) {
          status = false;
          ostringstream error_message;
          error_message << SEQ_error[SEQR_MISSING_STATE] << " " << i;
          error.update((error_message.str()).c_str());
        }
      }
    }
  }

  if ((max_order < 1) || (max_order > ORDER)) {
    status = false;
    error.update(SEQ_error[SEQR_ORDER]);
  }

  if (status) {
    markov = new VariableOrderMarkov(ORDINARY , marginal_distribution[0]->nb_value , max_order , true);
    markov->build_previous_memory();

    markov->markov_data = new VariableOrderMarkovData(*this , SEQUENCE_COPY , false);

    seq = markov->markov_data;
    seq->state_variable_init();
    seq->build_transition_count(*markov , begin , !begin);
    seq->chain_data->estimation(*markov , true , estimator);

    if (display) {
      markov->transition_count_ascii_write(cout , begin);
    }

    if (!path.empty()) {
      ofstream out_file(path.c_str());

      if (!out_file) {
        status = false;
        error.update(STAT_error[STATR_FILE_NAME]);
      }

      else {
        status = true;
        markov->transition_count_ascii_write(out_file , begin);
      }
    }

    delete markov;
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of the results of a comparison of models for a sample of sequences.
 *
 *  \param[in,out] os         stream,
 *  \param[in]     nb_model   number of models,
 *  \param[in]     likelihood log-likelihoods,
 *  \param[in]     label      model label,
 *  \param[in]     exhaustive flag detail level,
 *  \param[in]     algorithm  type of algorithm (NO_LATENT_STRUCTURE/FORWARD/VITERBI).
 */
/*--------------------------------------------------------------*/

ostream& MarkovianSequences::likelihood_write(ostream &os , int nb_model , double **likelihood ,
                                              const char *label , bool exhaustive ,
                                              latent_structure_algorithm algorithm) const

{
  bool *status;
  int i , j , k , m;
  int buff , model , min , width[3] , *rank_cumul , **rank;
  double max_likelihood , likelihood_cumul;
  ios_base::fmtflags format_flags;


  format_flags = os.setf(ios::right , ios::adjustfield);

  // computation of column widths

  if (exhaustive) {
    width[0] = column_width(nb_sequence);

    width[1] = 0;
    for (i = 0;i < nb_sequence;i++) {
      buff = column_width(nb_model , likelihood[i]);
      if (buff > width[1]) {
        width[1] = buff;
      }

/*      for (j = 0;j < nb_model;j++) {
        buff = column_width(1 , likelihood[i] + j ,
                            (likelihood[i][j] == D_INF ? 1. : 1. / length[i]));
        if (buff > width[1]) {
          width[1] = buff;
        }
      } */
    }
    width[1] += ASCII_SPACE;

    width[2] = column_width(nb_model) + ASCII_SPACE;

    // writing of the matrix of log-likelihoods of each model for each sequence

    os << "             ";
    for (i = 0;i < nb_model;i++) {
      os << " | " << label << " " << i + 1;
    }

    if (nb_model == 2) {
      os << " | " << SEQ_label[SEQL_LIKELIHOOD_RATIO];
    }
    os << endl;

    for (i = 0;i < nb_sequence;i++) {
      os << SEQ_label[SEQL_SEQUENCE] << " ";
      os << setw(width[0]) << i + 1 << ":";
      for (j = 0;j < nb_model;j++) {
//        os << setw(width[1]) << (likelihood[i][j] == D_INF ? likelihood[i][j] : likelihood[i][j] / length[i]);
        os << setw(width[1]) << likelihood[i][j];
      }

      max_likelihood = likelihood[i][0];
      model = 0;
      for (j = 1;j < nb_model;j++) {
        if (likelihood[i][j] > max_likelihood) {
          model = j;
          max_likelihood = likelihood[i][j];
        }
      }
      os << setw(width[2]) << model + 1;

      if (nb_model == 2) {
        if (likelihood[i][0] > likelihood[i][1]) {
          os << "   " << exp(likelihood[i][1] - likelihood[i][0]);
        }
        else {
          os << "   " << exp(likelihood[i][0] - likelihood[i][1]);
        }
      }
      os << endl;
    }
    os << endl;
  }

  // extraction of model ranks for each sequence

  rank = new int*[nb_model];
  for (i = 0;i < nb_model;i++) {
    rank[i] = new int[nb_model + 1];
    for (j = 0;j < nb_model;j++) {
      rank[i][j] = 0;
    }
  }

  status = new bool[nb_model];

  for (i = 0;i < nb_sequence;i++) {
    for (j = 0;j < nb_model;j++) {
      status[j] = true;
    }

    j = 0;
    while (j < nb_model) {
      max_likelihood = 2 * D_INF;
      for (k = 0;k < nb_model;k++) {
        if ((status[k]) && (likelihood[i][k] > max_likelihood)) {
          model = k;
          max_likelihood = likelihood[i][k];
        }
      }
      status[model] = false;
      rank[model][j]++;

      m = 1;
      for (k = 0;k < nb_model;k++) {
        if ((status[k]) && (likelihood[i][k] == max_likelihood)) {
          status[k] = false;
          rank[k][j]++;
          m++;
        }
      }

      j += m;
    }
  }

  // extraction of model ranks for all the sequences

  rank_cumul = new int[nb_model];

  for (i = 0;i < nb_model;i++) {
    status[i] = true;
    rank_cumul[i] = 0;
    for (j = 0;j < nb_model;j++) {
      rank_cumul[i] += rank[i][j] * (j + 1);
    }
  }

  i = 0;
  while (i < nb_model) {
    min = nb_model * nb_sequence + 1;
    for (j = 0;j < nb_model;j++) {
      if ((status[j]) && (rank_cumul[j] < min)) {
        model = j;
        min = rank_cumul[j];
      }
    }
    status[model] = false;
    rank[model][nb_model] = i + 1;

    k = 1;
    for (j = 0;j < nb_model;j++) {
      if ((status[j]) && (rank_cumul[j] == min)) {
        status[j] = false;
        rank[j][nb_model] = i + 1;
        k++;
      }
    }

    i += k;
  }

  width[0] = column_width(nb_model);
  width[1] = column_width(nb_sequence) + ASCII_SPACE;

  for (i = 0;i < nb_model;i++) {
    os << label << " ";
    os << setw(width[0]) << i + 1 << ":";
    for (j = 0;j < nb_model;j++) {
      os << setw(width[1]) << rank[i][j];
    }
    os << "  (" << rank[i][nb_model] << ")";

    // computation of the log-likelihood of a model for the sequences

    if (exhaustive) {
      likelihood_cumul = 0.;
      for (j = 0;j < nb_sequence;j++) {
        if (likelihood[j][i] != D_INF) {
          likelihood_cumul += likelihood[j][i];
        }
        else {
          likelihood_cumul = D_INF;
          break;
        }
      }

      switch (algorithm) {
      case NO_LATENT_STRUCTURE :
        os << " | " << STAT_label[STATL_LIKELIHOOD];
        break;
      case FORWARD :
        os << " | " << SEQ_label[SEQL_OBSERVED_SEQUENCES_LIKELIHOOD];
        break;
      case VITERBI :
        os << " | " << SEQ_label[SEQL_STATE_SEQUENCES_LIKELIHOOD];
        break;
      }
/*      if (likelihood_cumul != D_INF) {
        os << ": " << likelihood_cumul / cumul_length;
      }
      else { */
        os << ": " << likelihood_cumul;
//      }
    }
    os << endl;
  }
  os << endl;

  for (i = 0;i < nb_model;i++) {
    delete [] rank[i];
  }
  delete [] rank;

  delete [] status;
  delete [] rank_cumul;

  os.setf(format_flags , ios::adjustfield);

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of the results of a comparison of models for a sample of sequences in a file.
 *
 *  \param[in] error      reference on a StatError object,
 *  \param[in] path       file path,
 *  \param[in] nb_model   number of models,
 *  \param[in] likelihood log-likelihoods,
 *  \param[in] label      model label,
 *  \param[in] algorithm  type of algorithm (NO_LATENT_STRUCTURE/FORWARD/VITERBI).
 *
 *  \return               error status.
 */
/*--------------------------------------------------------------*/

bool MarkovianSequences::likelihood_write(StatError &error , const string path ,
                                          int nb_model , double **likelihood , const char *label ,
                                          latent_structure_algorithm algorithm) const

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
    likelihood_write(out_file , nb_model , likelihood , label , true , algorithm);
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Comparison of variable-order Markov chains for a sample of sequences.
 *
 *  \param[in] error    reference on a StatError object,
 *  \param[in] display  flag for displaying the results of model comparison,
 *  \param[in] nb_model number of variable-order Markov chains,
 *  \param[in] imarkov  pointer on VariableOrderMarkov objects,
 *  \param[in] path     file path.
 *
 *  \return             error status.
 */
/*--------------------------------------------------------------*/

bool MarkovianSequences::comparison(StatError &error , bool display , int nb_model ,
                                    const VariableOrderMarkov **imarkov ,
                                    const string path) const

{
  bool status = true;
  int i , j;
  double **likelihood;


  error.init();

  if ((type[0] != INT_VALUE) && (type[0] != STATE)) {
    status = false;
    ostringstream correction_message;
    correction_message << STAT_variable_word[INT_VALUE] << " or " << STAT_variable_word[STATE];
    error.correction_update(STAT_error[STATR_VARIABLE_TYPE] , (correction_message.str()).c_str());
  }

  else if (!characteristics[0]) {
    for (i = 0;i < marginal_distribution[0]->nb_value;i++) {
      if (marginal_distribution[0]->frequency[i] == 0) {
        status = false;
        ostringstream error_message;
        error_message << SEQ_error[SEQR_MISSING_STATE] << " " << i;
        error.update((error_message.str()).c_str());
      }
    }
  }

  if (nb_variable > 1) {
    if (nb_variable > 2) {
      status = false;
      error.correction_update(STAT_error[STATR_NB_VARIABLE] , "1 or 2");
    }

    if ((type[1] != INT_VALUE) && (type[1] != STATE)) {
      status = false;
      ostringstream error_message , correction_message;
      error_message << STAT_label[STATL_VARIABLE] << " " << 2 << ": "
                    << STAT_error[STATR_VARIABLE_TYPE];
      correction_message << STAT_variable_word[INT_VALUE] << " or "
                         << STAT_variable_word[STATE];
      error.correction_update((error_message.str()).c_str() , (correction_message.str()).c_str());
    }

    else {
      if (test_hidden(1)) {
        status = false;
        ostringstream error_message;
        error_message << STAT_label[STATL_VARIABLE] << " " << 2 << ": "
                      << SEQ_error[SEQR_OVERLAP];
        error.update((error_message.str()).c_str());
      }

      if (!characteristics[1]) {
        for (i = 0;i < marginal_distribution[1]->nb_value;i++) {
          if (marginal_distribution[1]->frequency[i] == 0) {
            status = false;
            ostringstream error_message;
            error_message << STAT_label[STATL_VARIABLE] << " " << 2 << ": "
                          << STAT_error[STATR_MISSING_VALUE] << " " << i;
            error.update((error_message.str()).c_str());
          }
        }
      }
    }
  }

  for (i = 0;i < nb_model;i++) {
    if (imarkov[i]->nb_output_process + 1 != nb_variable) {
      status = false;
      ostringstream error_message;
      error_message << SEQ_label[SEQL_MARKOV_CHAIN] << " " << i + 1 << ": "
                    << STAT_error[STATR_NB_OUTPUT_PROCESS];
      error.update((error_message.str()).c_str());
    }

    else {
      if (imarkov[i]->state_process->nb_value < marginal_distribution[0]->nb_value) {
        status = false;
        ostringstream error_message;
        error_message << SEQ_label[SEQL_MARKOV_CHAIN] << " " << i + 1 << ": "
                      << SEQ_error[SEQR_NB_STATE];
        error.update((error_message.str()).c_str());
      }

      if (nb_variable == 2) {
        if (imarkov[i]->categorical_process[0]->nb_value < marginal_distribution[1]->nb_value) {
          status = false;
          ostringstream error_message;
          error_message << SEQ_label[SEQL_MARKOV_CHAIN] << " " << i + 1 << ": "
                        << STAT_error[STATR_NB_OUTPUT];
          error.update((error_message.str()).c_str());
        }
      }
    }
  }

  if (status) {
    likelihood = new double*[nb_sequence];
    for (i = 0;i < nb_sequence;i++) {
      likelihood[i] = new double[nb_model];
    }

    // for each sequence, computation of the log-likelihood for each model

    for (i = 0;i < nb_sequence;i++) {
      for (j = 0;j < nb_model;j++) {
        likelihood[i][j] = imarkov[j]->likelihood_computation(*this , i);
      }
    }

    if (display) {
      likelihood_write(cout , nb_model , likelihood , SEQ_label[SEQL_MARKOV_CHAIN] , true);
    }
    if (!path.empty()) {
      status = likelihood_write(error , path , nb_model , likelihood , SEQ_label[SEQL_MARKOV_CHAIN]);
    }

    for (i = 0;i < nb_sequence;i++) {
      delete [] likelihood[i];
    }
    delete [] likelihood;
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Simulation using a variable-order Markov chain.
 *
 *  \param[in] error               reference on a StatError object,
 *  \param[in] length_distribution sequence length frequency distribution,
 *  \param[in] counting_flag       flag on the computation of the counting distributions,
 *  \param[in] divergence_flag     flag on the computation of a Kullback-Leibler divergence.
 *
 *  \return                        VariableOrderMarkovData.
 */
/*--------------------------------------------------------------*/

VariableOrderMarkovData* VariableOrderMarkov::simulation(StatError &error ,
                                                         const FrequencyDistribution &length_distribution ,
                                                         bool counting_flag ,
                                                         bool divergence_flag) const

{
  bool status = true , hidden;
  int i , j , k;
  int memory , cumul_length , *decimal_scale , *pstate , **pioutput;
  variable_nature *itype;
  double buff , min_location , likelihood , **proutput;
  Distribution *weight , *restoration_weight;
  VariableOrderMarkov *markov;
  VariableOrderMarkovData *seq;


  seq = NULL;
  error.init();

  if ((length_distribution.nb_element < 1) || (length_distribution.nb_element > NB_SEQUENCE)) {
    status = false;
    error.update(SEQ_error[SEQR_NB_SEQUENCE]);
  }
  if (length_distribution.offset < 2) {
    status = false;
    error.update(SEQ_error[SEQR_SHORT_SEQUENCE_LENGTH]);
  }
  if (length_distribution.nb_value - 1 > MAX_LENGTH) {
    status = false;
    error.update(SEQ_error[SEQR_LONG_SEQUENCE_LENGTH]);
  }

  if (status) {
    cumul_length = 0;
    for (i = length_distribution.offset;i < length_distribution.nb_value;i++) {
      cumul_length += i * length_distribution.frequency[i];
    }

    if (cumul_length > CUMUL_LENGTH) {
      status = false;
      error.update(SEQ_error[SEQR_CUMUL_SEQUENCE_LENGTH]);
    }
  }

  if (status) {
    if (length_distribution.nb_value - 1 > COUNTING_MAX_LENGTH) {
      counting_flag = false;
    }
    hidden = CategoricalSequenceProcess::test_hidden(nb_output_process , categorical_process);

    // initializations

    itype = new variable_nature[nb_output_process + 1];

    itype[0] = STATE;
    for (i = 0;i < nb_output_process;i++) {
      if (!continuous_parametric_process[i]) {
        itype[i + 1] = INT_VALUE;
      }
      else {
        itype[i + 1] = REAL_VALUE;
      }
    }

    seq = new VariableOrderMarkovData(length_distribution , nb_output_process + 1 , itype);
    delete [] itype;

    seq->markov = new VariableOrderMarkov(*this , false);

    markov = seq->markov;
    markov->create_cumul();
    markov->cumul_computation();

    if (markov->nb_output_process > 0) {
      pioutput = new int*[markov->nb_output_process];
      proutput = new double*[markov->nb_output_process];

      decimal_scale = new int[markov->nb_output_process];

      for (i = 0;i < markov->nb_output_process;i++) {
        if (markov->continuous_parametric_process[i]) {
          switch (markov->continuous_parametric_process[i]->ident) {

          case GAUSSIAN : {
            min_location = fabs(markov->continuous_parametric_process[i]->observation[0]->location);
            for (j = 1;j < markov->nb_state;j++) {
              buff = fabs(markov->continuous_parametric_process[i]->observation[j]->location);
              if (buff < min_location) {
                min_location = buff;
              }
            }

            buff = (int)ceil(log(min_location) / log(10));
            if (buff < GAUSSIAN_MAX_NB_DECIMAL) {
              decimal_scale[i] = pow(10 , (GAUSSIAN_MAX_NB_DECIMAL - buff));
            }
            else {
              decimal_scale[i] = 1;
            }
 
#           ifdef MESSAGE
            cout << "\nScale: " << i + 1 << " " << decimal_scale[i] << endl;
#           endif

            break;
          }

          case VON_MISES : {
            switch (markov->continuous_parametric_process[i]->unit) {
            case DEGREE :
              decimal_scale[i] = DEGREE_DECIMAL_SCALE;
              break;
            case RADIAN :
              decimal_scale[i] = RADIAN_DECIMAL_SCALE;
              break;
            }

            for (j = 0;j < markov->nb_state;j++) {
              markov->continuous_parametric_process[i]->observation[j]->von_mises_cumul_computation();
            }
            break;
          }
          }
        }
      }
    }

    for (i = 0;i < seq->nb_sequence;i++) {
      pstate = seq->int_sequence[i][0];

      for (j = 0;j < markov->nb_output_process;j++) {
        switch (seq->type[j + 1]) {
        case INT_VALUE :
          pioutput[j] = seq->int_sequence[i][j + 1];
          break;
        case REAL_VALUE :
          proutput[j] = seq->real_sequence[i][j + 1];
          break;
        }
      }

      switch (markov->type) {
      case ORDINARY :
        *pstate = cumul_method(markov->nb_state , markov->cumul_initial);
        memory = markov->child[0][*pstate];
        break;
      case EQUILIBRIUM :
        memory = cumul_method(markov->nb_row , markov->cumul_initial);
        *pstate = markov->state[memory][0];
        break;
      }

      for (j = 0;j < markov->nb_output_process;j++) {
        if (markov->categorical_process[j]) {
          *pioutput[j] = markov->categorical_process[j]->observation[*pstate]->simulation();
        }
        else if (markov->discrete_parametric_process[j]) {
          *pioutput[j] = markov->discrete_parametric_process[j]->observation[*pstate]->simulation();
        }
        else {
          *proutput[j] = round(markov->continuous_parametric_process[j]->observation[*pstate]->simulation() * decimal_scale[j]) / decimal_scale[j];
        }
      }

      for (j = 1;j < seq->length[i];j++) {
        *++pstate = cumul_method(markov->nb_state , markov->cumul_transition[memory]);

        for (k = 0;k < markov->nb_output_process;k++) {
          if (markov->categorical_process[k]) {
            *++pioutput[k] = markov->categorical_process[k]->observation[*pstate]->simulation();
          }
          else if (markov->discrete_parametric_process[k]) {
            *++pioutput[k] = markov->discrete_parametric_process[k]->observation[*pstate]->simulation();
          }
          else {
            *++proutput[k] = round(markov->continuous_parametric_process[k]->observation[*pstate]->simulation() * decimal_scale[k]) / decimal_scale[k];
          }
        }

        memory = markov->next[memory][*pstate];
      }
    }

    markov->remove_cumul();

    if (markov->nb_output_process > 0) {
      delete [] pioutput;
      delete [] proutput;

      delete [] decimal_scale;

      for (i = 0;i < markov->nb_output_process;i++) {
        if ((markov->continuous_parametric_process[i]) &&
            (markov->continuous_parametric_process[i]->ident == VON_MISES)) {
          for (j = 0;j < markov->nb_state;j++) {
            delete [] markov->continuous_parametric_process[i]->observation[j]->cumul;
            markov->continuous_parametric_process[i]->observation[j]->cumul = NULL;
          }
        }
      }
    }

    // computation of the characteristics of the generated sequences

    seq->min_value[0] = 0;
    seq->max_value[0] = nb_state - 1;
    seq->build_marginal_frequency_distribution(0);

    for (i = 1;i < seq->nb_variable;i++) {
      seq->min_value_computation(i);
      seq->max_value_computation(i);

      seq->build_marginal_frequency_distribution(i);
      seq->min_interval_computation(i);
    }

    seq->build_transition_count(*markov);
    seq->build_observation_frequency_distribution(nb_state);
    seq->build_observation_histogram(nb_state);
    seq->build_characteristic();

/*    if ((seq->max_value[0] < nb_state - 1) || (!(seq->characteristics[0]))) {
      delete seq;
      seq = NULL;
      error.update(SEQ_error[SEQR_STATES_NOT_REPRESENTED]);
    }

    else if (!divergence_flag) { */
    if (!divergence_flag) {
      markov->characteristic_computation(*seq , counting_flag);

      // computation of the log-likelihood of the model for the generated sequences

      likelihood = markov->likelihood_computation(*seq);

      if (likelihood == D_INF) {
        likelihood = markov->likelihood_computation(*seq , I_DEFAULT);
      }

#     ifdef DEBUG
      else {
        cout << "\n" << STAT_label[STATL_LIKELIHOOD] << ": " << likelihood
             << " | " << markov->likelihood_computation(*seq , I_DEFAULT) << endl;
      }
#     endif

      if (hidden) {
        seq->restoration_likelihood = likelihood;
      }
      else {
        seq->likelihood = likelihood;
      }

      // computation of the mixtures of observation distributions (theoretical weights and weights deduced from the restoration)

      if (hidden) {
        weight = markov->state_process->weight_computation();
        restoration_weight = seq->weight_computation();

        for (i = 0;i < markov->nb_output_process;i++) {
          if (markov->categorical_process[i]) {
            delete markov->categorical_process[i]->weight;
            delete markov->categorical_process[i]->mixture;
            markov->categorical_process[i]->weight = new Distribution(*weight);
            markov->categorical_process[i]->mixture = markov->categorical_process[i]->mixture_computation(markov->categorical_process[i]->weight);

            delete markov->categorical_process[i]->restoration_weight;
            delete markov->categorical_process[i]->restoration_mixture;
            markov->categorical_process[i]->restoration_weight = new Distribution(*restoration_weight);
            markov->categorical_process[i]->restoration_mixture = markov->categorical_process[i]->mixture_computation(markov->categorical_process[i]->restoration_weight);
          }

          else if (markov->discrete_parametric_process[i]) {
            delete markov->discrete_parametric_process[i]->weight;
            delete markov->discrete_parametric_process[i]->mixture;
            markov->discrete_parametric_process[i]->weight = new Distribution(*weight);
            markov->discrete_parametric_process[i]->mixture = markov->discrete_parametric_process[i]->mixture_computation(markov->discrete_parametric_process[i]->weight);

            delete markov->discrete_parametric_process[i]->restoration_weight;
            delete markov->discrete_parametric_process[i]->restoration_mixture;
            markov->discrete_parametric_process[i]->restoration_weight = new Distribution(*restoration_weight);
            markov->discrete_parametric_process[i]->restoration_mixture = markov->discrete_parametric_process[i]->mixture_computation(markov->discrete_parametric_process[i]->restoration_weight);
          }

          else if (markov->continuous_parametric_process[i]) {
            delete markov->continuous_parametric_process[i]->weight;
            markov->continuous_parametric_process[i]->weight = new Distribution(*weight);

            delete markov->continuous_parametric_process[i]->restoration_weight;
            markov->continuous_parametric_process[i]->restoration_weight = new Distribution(*restoration_weight);
          }
        }

        delete weight;
        delete restoration_weight;
      }
    }
  }

  return seq;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Simulation using a variable-order Markov chain.
 *
 *  \param[in] error         reference on a StatError object,
 *  \param[in] nb_sequence   number of sequences,
 *  \param[in] length        sequence length,
 *  \param[in] counting_flag flag on the computation of the counting distributions.
 *
 *  \return                  VariableOrderMarkovData object.
 */
/*--------------------------------------------------------------*/

VariableOrderMarkovData* VariableOrderMarkov::simulation(StatError &error ,
                                                         int nb_sequence , int length ,
                                                         bool counting_flag) const

{
  bool status = true;
  VariableOrderMarkovData *seq;


  seq = NULL;
  error.init();

  if ((nb_sequence < 1) || (nb_sequence > NB_SEQUENCE)) {
    status = false;
    error.update(SEQ_error[SEQR_NB_SEQUENCE]);
  }
  if (length < 2) {
    status = false;
    error.update(SEQ_error[SEQR_SHORT_SEQUENCE_LENGTH]);
  }
  if (length > MAX_LENGTH) {
    status = false;
    error.update(SEQ_error[SEQR_LONG_SEQUENCE_LENGTH]);
  }

  if (status) {
    FrequencyDistribution length_distribution(length + 1);

    length_distribution.nb_element = nb_sequence;
    length_distribution.offset = length;
    length_distribution.max = nb_sequence;
    length_distribution.mean = length;
    length_distribution.variance = 0.;
    length_distribution.frequency[length] = nb_sequence;

    seq = simulation(error , length_distribution , counting_flag);
  }

  return seq;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Simulation using a variable-order Markov chain.
 *
 *  \param[in] error         reference on a StatError object,
 *  \param[in] nb_sequence   number of sequences,
 *  \param[in] iseq          reference on a MarkovianSequences object,
 *  \param[in] counting_flag flag on the computation of the counting distributions.
 *
 *  \return                  VariableOrderMarkovData object.
 */
/*--------------------------------------------------------------*/

VariableOrderMarkovData* VariableOrderMarkov::simulation(StatError &error , int nb_sequence ,
                                                         const MarkovianSequences &iseq ,
                                                         bool counting_flag) const

{
  FrequencyDistribution *length_distribution;
  VariableOrderMarkovData *seq;


  error.init();

  if ((nb_sequence < 1) || (nb_sequence > NB_SEQUENCE)) {
    seq = NULL;
    error.update(SEQ_error[SEQR_NB_SEQUENCE]);
  }

  else {
    length_distribution = iseq.length_distribution->frequency_scale(nb_sequence);

    seq = simulation(error , *length_distribution , counting_flag);
    delete length_distribution;
  }

  return seq;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of Kullback-Leibler divergences between variable-order Markov chains.
 *
 *  \param[in] error               reference on a StatError object,
 *  \param[in] display             flag for displaying the matrix of pairwise distances between models,
 *  \param[in] nb_model            number of variable-order Markov chains,
 *  \param[in] imarkov             pointer on VariableOrderMarkov objects,
 *  \param[in] length_distribution sequence length frequency distribution,
 *  \param[in] path                file path.
 *
 *  \return                        DistanceMatrix object.
 */
/*--------------------------------------------------------------*/

DistanceMatrix* VariableOrderMarkov::divergence_computation(StatError &error , bool display , int nb_model ,
                                                            const VariableOrderMarkov **imarkov ,
                                                            FrequencyDistribution **length_distribution ,
                                                            const string path) const

{
  bool status = true , lstatus;
  int i , j , k;
  int cumul_length , nb_failure;
  double **likelihood;
  long double divergence;
  const VariableOrderMarkov **markov;
  MarkovianSequences *iseq , *seq;
  VariableOrderMarkovData *simul_seq;
  DistanceMatrix *dist_matrix;
  ofstream *out_file;


  dist_matrix = NULL;
  error.init();

  for (i = 0;i < nb_model - 1;i++) {
    if (imarkov[i]->type != type) {
      status = false;
      ostringstream error_message;
      error_message << SEQ_label[SEQL_MARKOV_CHAIN] << " " << i + 2 << ": "
                    << SEQ_error[SEQR_MODEL_TYPE];
      error.update((error_message.str()).c_str());
    }

    if (imarkov[i]->nb_output_process == nb_output_process) {
      if (imarkov[i]->nb_state != nb_state) {
        status = false;
        ostringstream error_message;
        error_message << SEQ_label[SEQL_MARKOV_CHAIN] << " " << i + 2 << ": "
                      << SEQ_error[SEQR_NB_STATE];
        error.update((error_message.str()).c_str());
      }

      if (nb_output_process == 1) {
        if (imarkov[i]->categorical_process[0]->nb_value != categorical_process[0]->nb_value) {
          status = false;
          ostringstream error_message;
          error_message << SEQ_label[SEQL_MARKOV_CHAIN] << " " << i + 2 << ": "
                        << STAT_error[STATR_NB_OUTPUT];
          error.update((error_message.str()).c_str());
        }
      }
    }

    else if ((nb_output_process == 0) && (imarkov[i]->nb_output_process == 1)) {
      if (imarkov[i]->categorical_process[0]->nb_value != nb_state) {
        status = false;
        ostringstream error_message;
        error_message << SEQ_label[SEQL_MARKOV_CHAIN] << " " << i + 2 << ": "
                      << STAT_error[STATR_NB_OUTPUT];
        error.update((error_message.str()).c_str());
      }
    }

    else {  // if ((nb_output_process == 1) && (imarkov[i]->nb_output_process == 0))
      if (imarkov[i]->nb_state != categorical_process[0]->nb_value) {
        status = false;
        ostringstream error_message;
        error_message << SEQ_label[SEQL_MARKOV_CHAIN] << " " << i + 2 << ": "
                      << SEQ_error[SEQR_NB_STATE];
        error.update((error_message.str()).c_str());
      }
    }
  }

  for (i = 0;i < nb_model;i++) {
    lstatus = true;

    if ((length_distribution[i]->nb_element < 1) || (length_distribution[i]->nb_element > NB_SEQUENCE)) {
      lstatus = false;
      ostringstream error_message;
      error_message << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " "
                    << i + 1 << ": "  << SEQ_error[SEQR_NB_SEQUENCE];
      error.update((error_message.str()).c_str());
    }
    if (length_distribution[i]->offset < 2) {
      lstatus = false;
      ostringstream error_message;
      error_message << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " "
                    << i + 1 << ": "  << SEQ_error[SEQR_SHORT_SEQUENCE_LENGTH];
      error.update((error_message.str()).c_str());
    }
    if (length_distribution[i]->nb_value - 1 > MAX_LENGTH) {
      lstatus = false;
      ostringstream error_message;
      error_message << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " "
                    << i + 1 << ": "  << SEQ_error[SEQR_LONG_SEQUENCE_LENGTH];
      error.update((error_message.str()).c_str());
    }

    if (!lstatus) {
      status = false;
    }

    else {
      cumul_length = 0;
      for (j = length_distribution[i]->offset;j < length_distribution[i]->nb_value;j++) {
        cumul_length += j * length_distribution[i]->frequency[j];
      }

      if (cumul_length > CUMUL_LENGTH) {
        status = false;
        ostringstream error_message;
        error_message << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " "
                      << i + 1 << ": "  << SEQ_error[SEQR_CUMUL_SEQUENCE_LENGTH];
        error.update((error_message.str()).c_str());
      }
    }
  }

  if (status) {
    out_file = NULL;

    if (!path.empty()) {
      out_file = new ofstream(path.c_str());

      if (!out_file) {
        error.update(STAT_error[STATR_FILE_NAME]);
        if (display) {
          cout << error;
        }
      }
    }

    markov = new const VariableOrderMarkov*[nb_model];

    markov[0] = this;
    for (i = 1;i < nb_model;i++) {
      markov[i] = imarkov[i - 1];
    }

    dist_matrix = new DistanceMatrix(nb_model , SEQ_label[SEQL_MARKOV_CHAIN]);

    for (i = 0;i < nb_model;i++) {

      // generation of a sample of sequences using a variable-order Markov chain

      simul_seq = markov[i]->simulation(error , *length_distribution[i] , false , true);

      likelihood = new double*[simul_seq->nb_sequence];
      for (j = 0;j < simul_seq->nb_sequence;j++) {
        likelihood[j] = new double[nb_model];
      }

      for (j = 0;j < simul_seq->nb_sequence;j++) {
        likelihood[j][i] = markov[i]->likelihood_computation(*simul_seq , j);

        if ((display) && (likelihood[j][i] == D_INF)) {
          cout << "\nERROR - " << SEQ_error[SEQR_REFERENCE_MODEL] << ": " << i + 1 << endl;
        }
      }

      if (markov[i]->nb_output_process == 1) {
        iseq = simul_seq->remove_variable_1();
      }
      else {
        iseq = simul_seq;
      }

      // computation of the log-likelihood of each variable-order Markov chain for the sample of sequences

      for (j = 0;j < nb_model;j++) {
        if (j != i) {
          if (markov[j]->nb_output_process == 1) {
            seq = iseq->transcode(error , markov[j]->categorical_process[0]);
          }
          else {
            seq = iseq;
          }

          divergence = 0.;
          cumul_length = 0;
          nb_failure = 0;

          for (k = 0;k < seq->nb_sequence;k++) {
            likelihood[k][j] = markov[j]->likelihood_computation(*seq , k);

//            if (divergence != -D_INF) {
              if (likelihood[k][j] != D_INF) {
                divergence += likelihood[k][i] - likelihood[k][j];
                cumul_length += seq->length[k];
              }
              else {
                nb_failure++;
//                divergence = -D_INF;
              }
//            }
          }

          if ((display) && (nb_failure > 0)) {
            cout << "\nWARNING - " << SEQ_error[SEQR_REFERENCE_MODEL] << ": " << i + 1 << ", "
                 << SEQ_error[SEQR_TARGET_MODEL] << ": " << j + 1 << " - "
                 << SEQ_error[SEQR_DIVERGENCE_NB_FAILURE] << ": " << nb_failure << endl;
          }

//          if (divergence != -D_INF) {
            dist_matrix->update(i + 1 , j + 1 , divergence , cumul_length);
//          }

          if (markov[j]->nb_output_process == 1) {
            delete seq;
          }
        }
      }

      if (display) {
        cout << SEQ_label[SEQL_MARKOV_CHAIN] << " " << i + 1 << ": " << simul_seq->nb_sequence << " "
             << SEQ_label[SEQL_SIMULATED] << " " << SEQ_label[simul_seq->nb_sequence == 1 ? SEQL_SEQUENCE : SEQL_SEQUENCES] << endl;
        simul_seq->likelihood_write(cout , nb_model , likelihood , SEQ_label[SEQL_MARKOV_CHAIN]);
      }
      if (out_file) {
        *out_file << SEQ_label[SEQL_MARKOV_CHAIN] << " " << i + 1 << ": " << simul_seq->nb_sequence << " "
                  << SEQ_label[SEQL_SIMULATED] << " " << SEQ_label[simul_seq->nb_sequence == 1 ? SEQL_SEQUENCE : SEQL_SEQUENCES] << endl;
        simul_seq->likelihood_write(*out_file , nb_model , likelihood , SEQ_label[SEQL_MARKOV_CHAIN]);
      }

      for (j = 0;j < simul_seq->nb_sequence;j++) {
        delete [] likelihood[j];
      }
      delete [] likelihood;

      if (markov[i]->nb_output_process == 1) {
        delete iseq;
      }
      delete simul_seq;
    }

    if (out_file) {
      out_file->close();
      delete out_file;
    }

    delete markov;
  }

  return dist_matrix;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of Kullback-Leibler divergences between variable-order Markov chains.
 *
 *  \param[in] error       reference on a StatError object,
 *  \param[in] display     flag for displaying the matrix of pairwise distances between models,
 *  \param[in] nb_model    number of variable-order Markov chains,
 *  \param[in] markov      pointer on VariableOrderMarkov objects,
 *  \param[in] nb_sequence number of sequences,
 *  \param[in] length      sequence length,
 *  \param[in] path        file path.
 *
 *  \return                DistanceMatrix object.
 */
/*--------------------------------------------------------------*/

DistanceMatrix* VariableOrderMarkov::divergence_computation(StatError &error , bool display , int nb_model ,
                                                            const VariableOrderMarkov **markov ,
                                                            int nb_sequence , int length , const string path) const

{
  bool status = true;
  int i;
  FrequencyDistribution **length_distribution;
  DistanceMatrix *dist_matrix;


  dist_matrix = NULL;
  error.init();

  if ((nb_sequence < 1) || (nb_sequence > NB_SEQUENCE)) {
    status = false;
    error.update(SEQ_error[SEQR_NB_SEQUENCE]);
  }
  if (length < 2) {
    status = false;
    error.update(SEQ_error[SEQR_SHORT_SEQUENCE_LENGTH]);
  }
  if (length > MAX_LENGTH) {
    status = false;
    error.update(SEQ_error[SEQR_LONG_SEQUENCE_LENGTH]);
  }

  if (status) {
    length_distribution = new FrequencyDistribution*[nb_model];

    length_distribution[0] = new FrequencyDistribution(length + 1);

    length_distribution[0]->nb_element = nb_sequence;
    length_distribution[0]->offset = length;
    length_distribution[0]->max = nb_sequence;
    length_distribution[0]->mean = length;
    length_distribution[0]->variance = 0.;
    length_distribution[0]->frequency[length] = nb_sequence;

    for (i = 1;i < nb_model;i++) {
      length_distribution[i] = new FrequencyDistribution(*length_distribution[0]);
    }

    dist_matrix = divergence_computation(error , display , nb_model , markov , length_distribution , path);

    for (i = 0;i < nb_model;i++) {
      delete length_distribution[i];
    }
    delete [] length_distribution;
  }

  return dist_matrix;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of Kullback-Leibler divergences between variable-order Markov chains.
 *
 *  \param[in] error       reference on a StatError object,
 *  \param[in] display     flag for displaying the matrix of pairwise distances between models,
 *  \param[in] nb_model    number of variable-order Markov chains,
 *  \param[in] markov      pointer on VariableOrderMarkov objects,
 *  \param[in] nb_sequence number of generated sequences,
 *  \param[in] seq         pointer on MarkovianSequences objects,
 *  \param[in] path        file path.
 *
 *  \return                DistanceMatrix object.
 */
/*--------------------------------------------------------------*/

DistanceMatrix* VariableOrderMarkov::divergence_computation(StatError &error , bool display , int nb_model ,
                                                            const VariableOrderMarkov **markov ,
                                                            int nb_sequence , const MarkovianSequences **seq ,
                                                            const string path) const

{
  int i;
  FrequencyDistribution **length_distribution;
  DistanceMatrix *dist_matrix;


  error.init();

  if ((nb_sequence < 1) || (nb_sequence > NB_SEQUENCE)) {
    dist_matrix = NULL;
    error.update(SEQ_error[SEQR_NB_SEQUENCE]);
  }

  else {
    length_distribution = new FrequencyDistribution*[nb_model];
    for (i = 0;i < nb_model;i++) {
      length_distribution[i] = seq[i]->length_distribution->frequency_scale(nb_sequence);
    }

    dist_matrix = divergence_computation(error , display , nb_model , markov , length_distribution , path);

    for (i = 0;i < nb_model;i++) {
      delete length_distribution[i];
    }
    delete [] length_distribution;
  }

  return dist_matrix;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Constructor of the VariableOrderMarkovIterator class.
 *
 *  \param[in] imarkov pointer on a VariableOrderMarkov object.
 */
/*--------------------------------------------------------------*/

VariableOrderMarkovIterator::VariableOrderMarkovIterator(VariableOrderMarkov *imarkov)

{
  markov = imarkov;
  (markov->nb_iterator)++;

  if ((!(markov->cumul_initial)) || (!(markov->cumul_transition))) {
    markov->create_cumul();
    markov->cumul_computation();
  }

  memory = I_DEFAULT;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Copy of a VariableOrderMarkovIterator object.
 *
 *  \param[in] iter reference on a VariableOrderMarkovIterator object.
 */
/*--------------------------------------------------------------*/

void VariableOrderMarkovIterator::copy(const VariableOrderMarkovIterator &iter)

{
  markov = iter.markov;
  (markov->nb_iterator)++;

  memory = iter.memory;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Destructor of the VariableOrderMarkovIterator class.
 */
/*--------------------------------------------------------------*/

VariableOrderMarkovIterator::~VariableOrderMarkovIterator()

{
  (markov->nb_iterator)--;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Assignment operator of the VariableOrderMarkovIterator class.
 *
 *  \param[in] iter reference on a VariableOrderMarkovIterator object.
 *
 *  \return         VariableOrderMarkovIterator object.
 */
/*--------------------------------------------------------------*/

VariableOrderMarkovIterator& VariableOrderMarkovIterator::operator=(const VariableOrderMarkovIterator &iter)

{
  if (&iter != this) {
    (markov->nb_iterator)--;
    copy(iter);
  }

  return *this;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Simulation using a variable-order Markov chain.
 *
 *  \param[in] int_seq        sequence,
 *  \param[in] length         sequence length,
 *  \param[in] initialization flag initialization.
 *
 *  \return                   error status.
 */
/*--------------------------------------------------------------*/

bool VariableOrderMarkovIterator::simulation(int **int_seq , int length , bool initialization)

{
  bool status;


  if ((memory == I_DEFAULT) && (!initialization)) {
    status = false;
  }

  else {
    int i , j;
    int offset = 0 , *pstate , **pioutput;
//    double **proutput;


    status = true;

    if (markov->nb_output_process > 0) {
      pioutput = new int*[markov->nb_output_process];
//      proutput = new double*[markov->nb_output_process];
    }

    pstate = int_seq[0];
    for (i = 0;i < markov->nb_output_process;i++) {
/*      switch (type[i + 1]) {
      case INT_VALUE : */
        pioutput[i] = int_seq[i + 1];
/*        break;
      case REAL_VALUE :
        proutput[i] = real_seq[i + 1];
        break;
      } */
    }

    if (initialization) {
      switch (markov->type) {
      case ORDINARY :
        *pstate = cumul_method(markov->nb_state , markov->cumul_initial);
        memory = markov->child[0][*pstate];
        break;
      case EQUILIBRIUM :
        memory = cumul_method(markov->nb_row , markov->cumul_initial);
        *pstate = markov->state[memory][0];
        break;
      }

      for (i = 0;i < markov->nb_output_process;i++) {
        if (markov->categorical_process[i]) {
          *pioutput[i]++ = markov->categorical_process[i]->observation[*pstate]->simulation();
        }
        else if (markov->discrete_parametric_process[i]) {
          *pioutput[i]++ = markov->discrete_parametric_process[i]->observation[*pstate]->simulation();
        }
        else {
//          *proutput[i]++ = markov->continuous_parametric_process[i]->observation[*pstate]->simulation();
        }
      }

      pstate++;
      offset++;
    }

    for (i = offset;i < length;i++) {
      *pstate = cumul_method(markov->nb_state , markov->cumul_transition[memory]);

      for (j = 0;j < markov->nb_output_process;j++) {
        if (markov->categorical_process[j]) {
          *pioutput[j]++ = markov->categorical_process[j]->observation[*pstate]->simulation();
        }
        else if (markov->categorical_process[j]) {
          *pioutput[j]++ = markov->discrete_parametric_process[j]->observation[*pstate]->simulation();
        }
        else {
//          *proutput[j]++ = markov->continuous_parametric_process[j]->observation[*pstate]->simulation();
        }
      }

      memory = markov->next[memory][*pstate++];
    }

    if (markov->nb_output_process > 0) {
      delete [] pioutput;
//      delete [] proutput;
    }
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Simulation using a variable-order Markov chain.
 *
 *  \param[in] length         sequence length,
 *  \param[in] initialization flag initialization.
 *
 *  \return                   generated sequence.
 */
/*--------------------------------------------------------------*/

int** VariableOrderMarkovIterator::simulation(int length , bool initialization)

{
  int i;
  int **int_seq;


  if ((memory == I_DEFAULT) && (!initialization)) {
    int_seq = NULL;
  }

  else {
    int_seq = new int*[markov->nb_output_process + 1];
    for (i = 0;i <= markov->nb_output_process;i++) {
      int_seq[i] = new int[length];
    }

    simulation(int_seq , length , initialization);
  }

  return int_seq;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Correction of the log-likelihood of a variable-order Markov chain for sequences.
 *
 *  \param[in] seq reference on a VariableOrderMarkovData object.
 *
 *  \return          log-likelihood.
 */
/*--------------------------------------------------------------*/

double VariableOrderMarkov::likelihood_correction(const VariableOrderMarkovData &seq) const

{
  int i;
  double correction;


  correction = 0.;

  for (i = 0;i < nb_state;i++) {
    if (seq.chain_data->initial[i] > 0) {
      correction += seq.chain_data->initial[i] * log(initial[i]);
    }
  }

  if (nb_output_process > 0) {
    for (i = 0;i < seq.nb_sequence;i++) {
      correction += log(categorical_process[0]->observation[seq.int_sequence[i][0][0]]->mass[seq.int_sequence[i][1][0]]);
    }
  }

  return correction;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Estimation of a lumped Markov chain.
 *
 *  \param[in] error         reference on a StatError object,
 *  \param[in] display       flag for displaying the results of the estimation,
 *  \param[in] category      transcoding table,
 *  \param[in] criterion     model selection criterion (AIC(c)/BIC),
 *  \param[in] order         Markov chain order,
 *  \param[in] counting_flag flag on the computation of the counting distributions.
 *
 *  \return                  VariableOrderMarkov object.
 */
/*--------------------------------------------------------------*/

VariableOrderMarkov* MarkovianSequences::lumpability_estimation(StatError &error , bool display , int *category ,
                                                                model_selection_criterion criterion ,
                                                                int order , bool counting_flag) const

{
  bool status = true , *presence;
  int i;
  int max_category , nb_state[2] , nb_parameter[2];
  double penalty , max_likelihood , likelihood[2] , penalized_likelihood[2];
  VariableOrderMarkov *markov , *lumped_markov;
  MarkovianSequences *seq;


  markov = NULL;
  error.init();

  if (nb_variable > 1) {
    status = false;
    error.correction_update(STAT_error[STATR_NB_VARIABLE] , 1);
  }

  if ((type[0] != INT_VALUE) && (type[0] != STATE)) {
    status = false;
    ostringstream correction_message;
    correction_message << STAT_variable_word[INT_VALUE] << " or " << STAT_variable_word[STATE];
    error.correction_update(STAT_error[STATR_VARIABLE_TYPE] , (correction_message.str()).c_str());
  }

  else {
    max_category = 0;
    for (i = 0;i < marginal_distribution[0]->nb_value;i++) {
      if ((category[i] < 0) || (category[i] >= marginal_distribution[0]->nb_value - 1)) {
        status = false;
        ostringstream error_message;
        error_message << STAT_label[STATL_CATEGORY] << " " << category[i] << " "
                      << STAT_error[STATR_NOT_ALLOWED];
        error.update((error_message.str()).c_str());
      }
      else if (category[i] > max_category) {
        max_category = category[i];
      }
    }

    if (max_category == 0) {
      status = false;
      error.update(STAT_error[STATR_NB_CATEGORY]);
    }

    if (status) {
      presence = new bool[max_category + 1];
      for (i = 0;i <= max_category;i++) {
        presence[i] = false;
      }

      for (i = 0;i < marginal_distribution[0]->nb_value;i++) {
        presence[category[i]] = true;
      }

      for (i = 0;i <= max_category;i++) {
        if (!presence[i]) {
          status = false;
          ostringstream error_message;
          error_message << STAT_error[STATR_MISSING_CATEGORY] << " " << i;
          error.update((error_message.str()).c_str());
        }
      }

      delete [] presence;
    }
  }

  if ((order < 1) || (order > ORDER)) {
    status = false;
    error.update(SEQ_error[SEQR_ORDER]);
  }

  if (status) {
//    markov = variable_order_markov_estimation(error , type , order , true , false);
    markov = variable_order_markov_estimation(error , ORDINARY , order , true , false);

    if (markov) {

      // computation of the compensation term

      switch (criterion) {
      case AIC :
        penalty = 1.;
        break;
      case BIC :
        penalty = 0.5 * log((double)cumul_length);
        break;
      }

      nb_state[1] = markov->nb_state;
      nb_parameter[1] = markov->nb_parameter_computation();
      likelihood[1] = markov->markov_data->likelihood -
                      markov->likelihood_correction(*(markov->markov_data));

      if (criterion == AICc) {
        if (nb_parameter[1] < cumul_length - 1) {
          penalized_likelihood[1] = likelihood[1] - (double)(nb_parameter[1] * cumul_length) /
                                    (double)(cumul_length - nb_parameter[1] - 1);
        }
        else {
          penalized_likelihood[1] = D_INF;
        }
      }

      else {
        penalized_likelihood[1] = likelihood[1] - nb_parameter[1] * penalty;
      }

      max_likelihood = penalized_likelihood[1];

      seq = transcode(error , 1 , category , true);

//      lumped_markov = seq->variable_order_markov_estimation(error , type , order , true , false);
      lumped_markov = seq->variable_order_markov_estimation(error , ORDINARY , order , true , false);

      if (lumped_markov) {
        if (display) {
          int j , k;
          int nb_output , sum , lumped_nb_parameter , *pstate , *poutput , *pfrequency ,
              ***observation_data;
          double lumped_likelihood , lumped_penalized_likelihood;


          // 2nd lumpability property (two-state- i.e. transition-dependent observation probabilities)

          observation_data = new int**[seq->marginal_distribution[0]->nb_value];
          for (i = 0;i < seq->marginal_distribution[0]->nb_value;i++) {
            observation_data[i] = new int*[seq->marginal_distribution[1]->nb_value];
            for (j = 0;j < seq->marginal_distribution[1]->nb_value;j++) {
              observation_data[i][j] = new int[seq->marginal_distribution[1]->nb_value];
            }
            for (j = 0;j < seq->marginal_distribution[0]->nb_value;j++) {
              pfrequency = observation_data[i][j];
              for (k = 0;k < seq->marginal_distribution[1]->nb_value;k++) {
                *pfrequency++ = 0;
              }
            }
          }

          // accumulation of the observation frequencies

          for (i = 0;i < seq->nb_sequence;i++) {
            pstate = seq->int_sequence[i][0] + 1;
            poutput = seq->int_sequence[i][1] + 1;
            for (j = 1;j < seq->length[i];j++) {
              observation_data[*(pstate - 1)][*pstate][*poutput]++;
              pstate++;
              poutput++;
            }
          }

          // estimation of the observation probabilities, computation of the log-likelihood and
          // of the number of free parameters

          lumped_nb_parameter = lumped_markov->nb_parameter_computation() -
                                lumped_markov->categorical_process[0]->nb_parameter_computation(0.);
          lumped_likelihood = lumped_markov->likelihood_computation(*(lumped_markov->markov_data->chain_data));

          if (criterion == AICc) {
            if (lumped_nb_parameter < cumul_length - 1) {
              lumped_penalized_likelihood = lumped_likelihood - (double)(lumped_nb_parameter * cumul_length) /
                                            (double)(cumul_length - lumped_nb_parameter - 1);
            }
            else {
              lumped_penalized_likelihood = D_INF;
            }
          }

          else {
            lumped_penalized_likelihood = lumped_likelihood - lumped_nb_parameter * penalty;
          }

          cout << "\n" << lumped_markov->nb_state << " " << STAT_label[STATL_STATES]
               << "   2 * " << SEQ_label[SEQL_MARKOV_CHAIN] << " " << STAT_label[STATL_LIKELIHOOD] << ": "
               << 2 * lumped_likelihood << "   " << lumped_nb_parameter << " "
               << STAT_label[lumped_nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS]
               << "   2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("
               << STAT_criterion_word[criterion] << "): " << 2 * lumped_penalized_likelihood << endl;

          cout << "\n" << STAT_word[STATW_OBSERVATION_PROBABILITIES] << endl;

          for (i = 0;i < seq->marginal_distribution[0]->nb_value;i++) {
            for (j = 0;j < seq->marginal_distribution[0]->nb_value;j++) {
              nb_output = 0;
              sum = 0;
              pfrequency = observation_data[i][j];
              for (k = 0;k < seq->marginal_distribution[1]->nb_value;k++) {
                if (*pfrequency > 0) {
                  nb_output++;
                  sum += *pfrequency;
                }
                pfrequency++;
              }

              if (nb_output > 1) {
                cout << i << " -> " << j << " : ";

                lumped_nb_parameter += (nb_output - 1);

                pfrequency = observation_data[i][j];
                for (k = 0;k < seq->marginal_distribution[1]->nb_value;k++) {
                  if (*pfrequency > 0) {
                    cout << k << " (" << (double)*pfrequency / (double)sum << ") | ";

                    lumped_likelihood += *pfrequency * log((double)*pfrequency / (double)sum);
                  }
                  pfrequency++;
                }

                cout << endl;
              }
            }
          }

          if (criterion == AICc) {
            if (lumped_nb_parameter < cumul_length - 1) {
              lumped_penalized_likelihood = lumped_likelihood - (double)(lumped_nb_parameter * cumul_length) /
                                            (double)(cumul_length - lumped_nb_parameter - 1);
            }
            else {
              lumped_penalized_likelihood = D_INF;
            }
          }

          else {
            lumped_penalized_likelihood = lumped_likelihood - lumped_nb_parameter * penalty;
          }

          cout << "\n" << lumped_markov->nb_state << " " << STAT_label[STATL_STATES]
               << "   2 * " << STAT_label[STATL_LIKELIHOOD] << ": " << 2 * lumped_likelihood << "   "
               << lumped_nb_parameter << " " << STAT_label[lumped_nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS]
               << "   2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("
               << STAT_criterion_word[criterion] << "): " << 2 * lumped_penalized_likelihood << endl;

          // 3rd lumpability property (output-state-dependent observation probabilities)

          for (i = 0;i < seq->marginal_distribution[0]->nb_value;i++) {
            for (j = 0;j < seq->marginal_distribution[1]->nb_value;j++) {
              pfrequency = observation_data[i][j];
              for (k = 0;k < seq->marginal_distribution[1]->nb_value;k++) {
                *pfrequency++ = 0;
              }
            }
          }

          // accumulation of the observation frequencies

          for (i = 0;i < seq->nb_sequence;i++) {
            pstate = seq->int_sequence[i][0] + 1;
            poutput = seq->int_sequence[i][1] + 1;
            for (j = 1;j < seq->length[i];j++) {
              observation_data[*pstate][*(poutput - 1)][*poutput]++;
              pstate++;
              poutput++;
            }
          }

          // estimation of the observation probabilities, computation of the log-likelihood and
          // of the number of free parameters

          lumped_nb_parameter = lumped_markov->nb_parameter_computation() -
                                lumped_markov->categorical_process[0]->nb_parameter_computation(0.);
          lumped_likelihood = lumped_markov->likelihood_computation(*(lumped_markov->markov_data->chain_data));

          cout << "\n" << STAT_word[STATW_OBSERVATION_PROBABILITIES] << endl;

          for (i = 0;i < seq->marginal_distribution[0]->nb_value;i++) {
            for (j = 0;j < seq->marginal_distribution[1]->nb_value;j++) {
              nb_output = 0;
              sum = 0;
              pfrequency = observation_data[i][j];
              for (k = 0;k < seq->marginal_distribution[1]->nb_value;k++) {
                if (*pfrequency > 0) {
                  nb_output++;
                  sum += *pfrequency;
                }
                pfrequency++;
              }

              if (nb_output > 1) {
                cout << j << ", " << i << " : ";

                lumped_nb_parameter += (nb_output - 1);

                pfrequency = observation_data[i][j];
                for (k = 0;k < seq->marginal_distribution[1]->nb_value;k++) {
                  if (*pfrequency > 0) {
                    cout << k << " (" << (double)*pfrequency / (double)sum << ") | ";

                    lumped_likelihood += *pfrequency * log((double)*pfrequency / (double)sum);
                  }
                  pfrequency++;
                }

                cout << endl;
              }
            }
          }

          if (criterion == AICc) {
            if (lumped_nb_parameter < cumul_length - 1) {
              lumped_penalized_likelihood = lumped_likelihood - (double)(lumped_nb_parameter * cumul_length) /
                                            (double)(cumul_length - lumped_nb_parameter - 1);
            }
            else {
              lumped_penalized_likelihood = D_INF;
            }
          }

          else {
            lumped_penalized_likelihood = lumped_likelihood - lumped_nb_parameter * penalty;
          }

          cout << "\n" << lumped_markov->nb_state << " " << STAT_label[STATL_STATES]
               << "   2 * " << STAT_label[STATL_LIKELIHOOD] << ": " << 2 * lumped_likelihood << "   "
               << lumped_nb_parameter << " " << STAT_label[lumped_nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS]
               << "   2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("
               << STAT_criterion_word[criterion] << "): " << 2 * lumped_penalized_likelihood << endl;

          for (i = 0;i < seq->marginal_distribution[0]->nb_value;i++) {
            for (j = 0;j < seq->marginal_distribution[1]->nb_value;j++) {
              delete [] observation_data[i][j];
            }
            delete [] observation_data[i];
          }
          delete [] observation_data;
        }

        nb_state[0] = lumped_markov->nb_state;
        nb_parameter[0] = lumped_markov->nb_parameter_computation();
        likelihood[0] = lumped_markov->markov_data->likelihood -
                        lumped_markov->likelihood_correction(*(lumped_markov->markov_data));

        if (criterion == AICc) {
          if (nb_parameter[0] < cumul_length - 1) {
            penalized_likelihood[0] = likelihood[0] - (double)(nb_parameter[0] * cumul_length) /
                                      (double)(cumul_length - nb_parameter[0] - 1);
          }
          else {
            penalized_likelihood[0] = D_INF;
          }
        }

        else {
          penalized_likelihood[0] = likelihood[0] - nb_parameter[0] * penalty;
        }

#       ifdef DEBUG
        if (penalized_likelihood[0] > max_likelihood) {
          markov->ascii_write(os);
        }
        else {
          lumped_markov->ascii_write(os);
        }
#       endif

        if (penalized_likelihood[0] > max_likelihood) {
          max_likelihood = penalized_likelihood[0];
          delete markov;
          markov = lumped_markov;
        }
        else {
          delete lumped_markov;
        }

#       ifdef DEBUG
        lumpability_test(error , category , cout , order);
#       endif

#       ifdef MESSAGE
        if (display) {
/*          double norm = 0. , weight[2];

          for (i = 0;i < 2;i++) {
            weight[i] = exp(penalized_likelihood[i] - max_likelihood);
            norm += weight[i];
          } */

          for (i = 0;i < 2;i++) {
            cout << "\n" << nb_state[i] << " " << STAT_label[STATL_STATES]
                 << "   2 * " << STAT_label[STATL_LIKELIHOOD] << ": " << 2 * likelihood[i] << "   "
                 << nb_parameter[i] << " " << STAT_label[nb_parameter[i] == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS]
                 << "   2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("
                 << STAT_criterion_word[criterion] << "): " << 2 * penalized_likelihood[i] << endl;
//                 << "   " << STAT_label[STATL_WEIGHT] << ": " << weight[i] / norm << endl;
          }
        }
#       endif

      }

      delete seq;
    }

    // computation of the characteristic distributions of the model

    markov->component_computation();
    markov->characteristic_computation(*(markov->markov_data) , counting_flag , I_DEFAULT , false);
  }

  return markov;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Test of state lumpability for a Markov chain.
 *
 *  \param[in] error    reference on a StatError object,
 *  \param[in] display  flag for displaying the test results,
 *  \param[in] category transcoding table,
 *  \param[in] order    Markov chain order.
 *
 *  \return             error status.
 */
/*--------------------------------------------------------------*/

bool MarkovianSequences::lumpability_test(StatError &error , bool display ,
                                          int *category , int order) const

{
  bool status = true , *presence;
  int i , j , k;
  int max_category , df , sum , *ftransition;
  double value , var1 , var2;
  Test *test;
  VariableOrderMarkov *markov , *lumped_markov;
  MarkovianSequences *seq;


  error.init();

  if (nb_variable > 1) {
    status = false;
    error.correction_update(STAT_error[STATR_NB_VARIABLE] , 1);
  }

  if ((type[0] != INT_VALUE) && (type[0] != STATE)) {
    status = false;
    ostringstream correction_message;
    correction_message << STAT_variable_word[INT_VALUE] << " or " << STAT_variable_word[STATE];
    error.correction_update(STAT_error[STATR_VARIABLE_TYPE] , (correction_message.str()).c_str());
  }

  else {
    if ((marginal_distribution[0]->nb_value < 2) ||
        (marginal_distribution[0]->nb_value > NB_STATE)) {
      status = false;
      error.update(SEQ_error[SEQR_NB_STATE]);
    }

    else if (!characteristics[0]) {
      for (i = 0;i < marginal_distribution[0]->nb_value;i++) {
        if (marginal_distribution[0]->frequency[i] == 0) {
          status = false;
          ostringstream error_message;
          error_message << SEQ_error[SEQR_MISSING_STATE] << " " << i;
          error.update((error_message.str()).c_str());
        }
      }
    }

    max_category = 0;
    for (i = 0;i < marginal_distribution[0]->nb_value;i++) {
      if ((category[i] < 0) || (category[i] >= marginal_distribution[0]->nb_value - 1)) {
        status = false;
        ostringstream error_message;
        error_message << STAT_label[STATL_CATEGORY] << " " << category[i] << " "
                      << STAT_error[STATR_NOT_ALLOWED];
        error.update((error_message.str()).c_str());
      }
      else if (category[i] > max_category) {
        max_category = category[i];
      }
    }

    if (max_category == 0) {
      status = false;
      error.update(STAT_error[STATR_NB_CATEGORY]);
    }

    if (status) {
      presence = new bool[max_category + 1];
      for (i = 0;i <= max_category;i++) {
        presence[i] = false;
      }

      for (i = 0;i < marginal_distribution[0]->nb_value;i++) {
        presence[category[i]] = true;
      }

      for (i = 0;i <= max_category;i++) {
        if (!presence[i]) {
          status = false;
          ostringstream error_message;
          error_message << STAT_error[STATR_MISSING_CATEGORY] << " " << i;
          error.update((error_message.str()).c_str());
        }
      }

      delete [] presence;
    }
  }

  if ((order < 1) || (order > ORDER)) {
    status = false;
    error.update(SEQ_error[SEQR_ORDER]);
  }

  if (status) {
//    markov = variable_order_markov_estimation(error , type , order , true , false);
    markov = variable_order_markov_estimation(error , ORDINARY , order , true , false);

    seq = transcode(error , 1 , category , true);

//    lumped_markov = seq->variable_order_markov_estimation(error , type , order , true , false);
    lumped_markov = seq->variable_order_markov_estimation(error , ORDINARY , order , true , false);

    df = markov->nb_parameter_computation() - lumped_markov->nb_parameter_computation();

    value = 0.;

    for (i = 1;i < markov->nb_row;i++) {
      if (markov->memo_type[i] == TERMINAL) {
        ftransition = markov->markov_data->chain_data->transition[i];
        sum = 0;
        for (j = 0;j < markov->nb_state;j++) {
          sum += *ftransition++;
        }

        if (sum > 0) {
          for (j = 1;j < lumped_markov->nb_row;j++) {
            if ((lumped_markov->memo_type[j] == TERMINAL) &&
                (lumped_markov->order[j] == markov->order[i])) {
              for (k = 0;k < lumped_markov->order[j];k++) {
                if (lumped_markov->state[j][k] != category[markov->state[i][k]]) {
                  break;
                }
              }

              if (k == lumped_markov->order[j]) {
                ftransition = markov->markov_data->chain_data->transition[i];
                for (k = 0;k < markov->nb_state;k++) {
                  var1 = (double)sum * lumped_markov->categorical_process[0]->observation[category[k]]->mass[k] *
                         lumped_markov->transition[j][category[k]];
                  if (var1 > 0.) {
                    var2 = *ftransition - var1;
                    value += var2 * var2 / var1;
                  }
                  ftransition++;
                }
                break;
              }
            }
          }
        }
      }
    }

    test = new Test(CHI2 , true , df , I_DEFAULT , value);

    test->chi2_critical_probability_computation();

#   ifdef MESSAGE
    cout << *test;
#   endif

    delete test;

    value = 2 * (markov->markov_data->likelihood - markov->likelihood_correction(*(markov->markov_data)) -
            (lumped_markov->markov_data->likelihood - lumped_markov->likelihood_correction(*(lumped_markov->markov_data))));

    test = new Test(CHI2 , true , df , I_DEFAULT , value);

    test->chi2_critical_probability_computation();

#   ifdef MESSAGE
    cout << "\n" << SEQ_label[SEQL_LIKELIHOOD_RATIO_TEST] << "\n" << *test;
#   endif

    delete test;

    if (display) {
      int k;
      int nb_output , sum , lumped_nb_parameter , *pstate , *poutput , *pfrequency ,
          ***observation_data;
      double lumped_likelihood , *pproba , ***observation_proba;


      // 2eme lumpability property (two-state- i.e. transition-dependent observation probabilities)

      observation_data = new int**[seq->marginal_distribution[0]->nb_value];
      observation_proba = new double**[seq->marginal_distribution[0]->nb_value];
      for (i = 0;i < seq->marginal_distribution[0]->nb_value;i++) {
        observation_data[i] = new int*[seq->marginal_distribution[1]->nb_value];
        observation_proba[i] = new double*[seq->marginal_distribution[1]->nb_value];
        for (j = 0;j < seq->marginal_distribution[1]->nb_value;j++) {
          observation_data[i][j] = new int[seq->marginal_distribution[1]->nb_value];
          observation_proba[i][j] = new double[seq->marginal_distribution[1]->nb_value];
        }
        for (j = 0;j < seq->marginal_distribution[0]->nb_value;j++) {
          pfrequency = observation_data[i][j];
          for (k = 0;k < seq->marginal_distribution[1]->nb_value;k++) {
            *pfrequency++ = 0;
          }
        }
      }

      // accumulation of the observation frequencies

      for (i = 0;i < seq->nb_sequence;i++) {
        pstate = seq->int_sequence[i][0] + 1;
        poutput = seq->int_sequence[i][1] + 1;
        for (j = 1;j < seq->length[i];j++) {
          observation_data[*(pstate - 1)][*pstate][*poutput]++;
          pstate++;
          poutput++;
        }
      }

      // estimation of the observation probabilities, computation of the log-likelihood and
      // of the number of free parameters

      lumped_nb_parameter = lumped_markov->nb_parameter_computation() -
                            lumped_markov->categorical_process[0]->nb_parameter_computation(0.);
      lumped_likelihood = lumped_markov->likelihood_computation(*(lumped_markov->markov_data->chain_data));

      for (i = 0;i < seq->marginal_distribution[0]->nb_value;i++) {
        for (j = 0;j < seq->marginal_distribution[0]->nb_value;j++) {
          nb_output = 0;
          sum = 0;
          pfrequency = observation_data[i][j];
          for (k = 0;k < seq->marginal_distribution[1]->nb_value;k++) {
            if (*pfrequency > 0) {
              nb_output++;
              sum += *pfrequency;
            }
            pfrequency++;
          }

          if (nb_output > 1) {
            lumped_nb_parameter += (nb_output - 1);

            pfrequency = observation_data[i][j];
            for (k = 0;k < seq->marginal_distribution[1]->nb_value;k++) {
              if (*pfrequency > 0) {
                lumped_likelihood += *pfrequency * log((double)*pfrequency / (double)sum);
              }
              pfrequency++;
            }
          }

          if (sum > 0) {
            pproba = observation_proba[i][j];
            pfrequency = observation_data[i][j];
            for (k = 0;k < seq->marginal_distribution[1]->nb_value;k++) {
              *pproba++ = (double)*pfrequency++ / (double)sum;
            }
          }
        }
      }

      df = markov->nb_parameter_computation() - lumped_nb_parameter;

      value = 0.;

      for (i = 1;i < markov->nb_row;i++) {
        if (markov->memo_type[i] == TERMINAL) {
          ftransition = markov->markov_data->chain_data->transition[i];
          sum = 0;
          for (j = 0;j < markov->nb_state;j++) {
            sum += *ftransition++;
          }

          if (sum > 0) {
            for (j = 1;j < lumped_markov->nb_row;j++) {
              if ((lumped_markov->memo_type[j] == TERMINAL) &&
                  (lumped_markov->order[j] == markov->order[i])) {
                for (k = 0;k < lumped_markov->order[j];k++) {
                  if (lumped_markov->state[j][k] != category[markov->state[i][k]]) {
                    break;
                  }
                }

                if (k == lumped_markov->order[j]) {
                  ftransition = markov->markov_data->chain_data->transition[i];
                  for (k = 0;k < markov->nb_state;k++) {
                    var1 = (double)sum * observation_proba[category[markov->state[i][0]]][category[k]][k] *
                           lumped_markov->transition[j][category[k]];
                    if (var1 > 0.) {
                      var2 = *ftransition - var1;
                      value += var2 * var2 / var1;
                    }
                    ftransition++;
                  }
                  break;
                }
              }
            }
          }
        }
      }

      test = new Test(CHI2 , true , df , I_DEFAULT , value);

      test->chi2_critical_probability_computation();

      cout << "\n" << *test;

      delete test;

      value = 2 * (markov->markov_data->likelihood - markov->likelihood_correction(*(markov->markov_data)) - lumped_likelihood);

      test = new Test(CHI2 , true , df , I_DEFAULT , value);

      test->chi2_critical_probability_computation();

      cout << "\n" << SEQ_label[SEQL_LIKELIHOOD_RATIO_TEST] << "\n" << *test;

      delete test;

      // 3rd lumpability property (output-state-dependent observation probabilities)

      for (i = 0;i < seq->marginal_distribution[0]->nb_value;i++) {
        for (j = 0;j < seq->marginal_distribution[1]->nb_value;j++) {
          pfrequency = observation_data[i][j];
          for (k = 0;k < seq->marginal_distribution[1]->nb_value;k++) {
            *pfrequency++ = 0;
          }
        }
      }

      // accumulation of the observation frequencies

      for (i = 0;i < seq->nb_sequence;i++) {
        pstate = seq->int_sequence[i][0] + 1;
        poutput = seq->int_sequence[i][1] + 1;
        for (j = 1;j < seq->length[i];j++) {
          observation_data[*pstate][*(poutput - 1)][*poutput]++;
          pstate++;
          poutput++;
        }
      }

      // estimation of the observation probabilities, computation of the log-likelihood and
      // of the number of free parameters

      lumped_nb_parameter = lumped_markov->nb_parameter_computation() -
                            lumped_markov->categorical_process[0]->nb_parameter_computation(0.);
      lumped_likelihood = lumped_markov->likelihood_computation(*(lumped_markov->markov_data->chain_data));

      for (i = 0;i < seq->marginal_distribution[0]->nb_value;i++) {
        for (j = 0;j < seq->marginal_distribution[1]->nb_value;j++) {
          nb_output = 0;
          sum = 0;
          pfrequency = observation_data[i][j];
          for (k = 0;k < seq->marginal_distribution[1]->nb_value;k++) {
            if (*pfrequency > 0) {
              nb_output++;
              sum += *pfrequency;
            }
            pfrequency++;
          }

          if (nb_output > 1) {
            lumped_nb_parameter += (nb_output - 1);

            pfrequency = observation_data[i][j];
            for (k = 0;k < seq->marginal_distribution[1]->nb_value;k++) {
              if (*pfrequency > 0) {
                lumped_likelihood += *pfrequency * log((double)*pfrequency / (double)sum);
              }
              pfrequency++;
            }
          }

          if (sum > 0) {
            pproba = observation_proba[i][j];
            pfrequency = observation_data[i][j];
            for (k = 0;k < seq->marginal_distribution[1]->nb_value;k++) {
              *pproba++ = (double)*pfrequency++ / (double)sum;
            }
          }
        }
      }

      df = markov->nb_parameter_computation() - lumped_nb_parameter;

      value = 0.;

      for (i = 1;i < markov->nb_row;i++) {
        if (markov->memo_type[i] == TERMINAL) {
          ftransition = markov->markov_data->chain_data->transition[i];
          sum = 0;
          for (j = 0;j < markov->nb_state;j++) {
            sum += *ftransition++;
          }

          if (sum > 0) {
            for (j = 1;j < lumped_markov->nb_row;j++) {
              if ((lumped_markov->memo_type[j] == TERMINAL) &&
                  (lumped_markov->order[j] == markov->order[i])) {
                for (k = 0;k < lumped_markov->order[j];k++) {
                  if (lumped_markov->state[j][k] != category[markov->state[i][k]]) {
                    break;
                  }
                }

                if (k == lumped_markov->order[j]) {
                  ftransition = markov->markov_data->chain_data->transition[i];
                  for (k = 0;k < markov->nb_state;k++) {
                    var1 = (double)sum * observation_proba[category[k]][markov->state[i][0]][k] *
                           lumped_markov->transition[j][category[k]];
                    if (var1 > 0.) {
                      var2 = *ftransition - var1;
                      value += var2 * var2 / var1;
                    }
                    ftransition++;
                  }
                  break;
                }
              }
            }
          }
        }
      }

      test = new Test(CHI2 , true , df , I_DEFAULT , value);

      test->chi2_critical_probability_computation();

      cout << "\n" << *test;

      delete test;

      value = 2 * (markov->markov_data->likelihood - markov->likelihood_correction(*(markov->markov_data)) - lumped_likelihood);

      test = new Test(CHI2 , true , df , I_DEFAULT , value);

      test->chi2_critical_probability_computation();

      cout << "\n" << SEQ_label[SEQL_LIKELIHOOD_RATIO_TEST] << "\n" << *test;

      delete test;

      for (i = 0;i < seq->marginal_distribution[0]->nb_value;i++) {
        for (j = 0;j < seq->marginal_distribution[1]->nb_value;j++) {
          delete [] observation_data[i][j];
          delete [] observation_proba[i][j];
        }
        delete [] observation_data[i];
        delete [] observation_proba[i];
      }
      delete [] observation_data;
      delete [] observation_proba;
    }

    delete markov;
    delete seq;
    delete lumped_markov;
  }

  return status;
}


};  // namespace sequence_analysis
