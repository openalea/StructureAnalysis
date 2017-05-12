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

#include <boost/math/distributions/normal.hpp>

#include "stat_tool/stat_label.h"

#include "semi_markov.h"
#include "sequence_label.h"

using namespace std;
using namespace boost::math;
using namespace stat_tool;


namespace sequence_analysis {



/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the stationary distribution for an equilibrium semi-Markov chain.
 */
/*--------------------------------------------------------------*/

void SemiMarkovChain::initial_probability_computation()

{
  int i , j , k;
  double sum , *state , *state_out , **state_in;
  DiscreteParametric *occupancy;


  state = new double[nb_state];
  state_out = new double[nb_state];

  state_in = new double*[STATIONARY_PROBABILITY_LENGTH];
  for (i = 0;i < STATIONARY_PROBABILITY_LENGTH;i++) {
    state_in[i] = new double[nb_state];
  }

  i = 0;

  do {
    if (i > 0) {
      sum = 0.;
    }

    for (j = 0;j < nb_state;j++) {
      if (i > 0) {
        sum += fabs(state_in[i - 1][j] - state_out[j]);
      }

      switch (sojourn_type[j]) {

      // case semi-Markovian state

      case SEMI_MARKOVIAN : {
        if (i == 0) {
          state[j] = initial[j];
        }
        else {
          state[j] += state_in[i - 1][j] - state_out[j];
        }

        occupancy = state_process->sojourn_time[j];
        state_out[j] = 0.;

        for (k = 1;k <= MIN(i + 1 , occupancy->nb_value - 1);k++) {
          if (k < i + 1) {
            state_out[j] += occupancy->mass[k] * state_in[i - k][j];
          }
          else {
            state_out[j] += forward[j]->mass[k] * initial[j];
          }
        }
        break;
      }

      // case Markovian state

      case MARKOVIAN : {
        if (i == 0) {
          state_out[j] = initial[j];
        }
        else {
          state_out[j] = state_in[i - 1][j];
        }
        break;
      }
      }
    }

    for (j = 0;j < nb_state;j++) {
      state_in[i][j] = 0.;
      for (k = 0;k < nb_state;k++) {
        state_in[i][j] += transition[k][j] * state_out[k];
      }
    }

#   ifdef DEBUG
//    if ((i > 0) && (i % 100 == 0)) {
      cout << i << "  ";
      for (j = 0;j < nb_state;j++) {
        cout << state[j] << " ";
      }
      cout << " | " << sum / nb_state << endl;
//    }
#   endif

    i++;
  }
  while (((i == 1) || (sum / nb_state > STATIONARY_PROBABILITY_THRESHOLD)) &&
         (i < STATIONARY_PROBABILITY_LENGTH));

# ifdef DEBUG
  cout << "\n" << SEQ_label[SEQL_LENGTH] << ": "  << i << endl;
# endif

  for (j = 0;j < nb_state;j++) {
    switch (sojourn_type[j]) {

    // case semi-Markovian state

    case SEMI_MARKOVIAN :
      initial[j] = state_in[i - 1][j] - state_out[j] + state[j];
//      initial[j] = state[j];
      break;

    // case Markovian state

    case MARKOVIAN :
      initial[j] = state_in[i - 1][j];
//      initial[j] = state_out[j];
      break;
    }
  }

  // renormalization for taking account of the thresholds applied on
  // the cumulative state occupancy distribution functions

  sum = 0.;
  for (i = 0;i < nb_state;i++) {
    sum += initial[i];
  }
  for (i = 0;i < nb_state;i++) {
    initial[i] /= sum;
  }

  delete [] state;
  delete [] state_out;

  for (i = 0;i < STATIONARY_PROBABILITY_LENGTH;i++) {
    delete [] state_in[i];
  }
  delete [] state_in;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the log-likelihood of a semi-Markov chain for sequences.
 *
 *  \param[in] seq   reference on a MarkovianSequences object,
 *  \param[in] index sequence index.
 *
 *  \return          log-likelihood.
 */
/*--------------------------------------------------------------*/

double SemiMarkov::likelihood_computation(const MarkovianSequences &seq , int index) const

{
  int i , j , k , m;
  int nb_value , occupancy , *pstate , **pioutput;
  double likelihood = 0. , proba , residual , **proutput;


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
        pstate = seq.int_sequence[i][0];

        proba = initial[*pstate];
        if (proba > 0.) {
          likelihood += log(proba);
        }
        else {
          likelihood = D_INF;
          break;
        }

        if (nb_output_process > 0) {
          for (j = 0;j < nb_output_process;j++) {
            switch (seq.type[j + 1]) {
            case INT_VALUE :
              pioutput[j] = seq.int_sequence[i][j + 1];
              break;
            case REAL_VALUE :
              proutput[j] = seq.real_sequence[i][j + 1];
              break;
            }
          }
        }

        j = 0;
        do {
          if (j > 0) {
            pstate++;

            proba = transition[*(pstate - 1)][*pstate];
            if (proba > 0.) {
              likelihood += log(proba);
            }
            else {
              likelihood = D_INF;
              break;
            }
          }

          if (transition[*pstate][*pstate] < 1.) {
            occupancy = 1;

            if (sojourn_type[*pstate] == SEMI_MARKOVIAN) {
              while ((j + occupancy < seq.length[i]) && (*(pstate + 1) == *pstate)) {
                occupancy++;
                pstate++;
              }

              proba = 0.;
              if ((type == EQUILIBRIUM) && (j == occupancy)) {
                if (occupancy < forward[*pstate]->nb_value) {
                  if (j + occupancy < seq.length[i]) {
                    proba = forward[*pstate]->mass[occupancy];
                  }
                  else {
                    proba = (1. - forward[*pstate]->cumul[occupancy - 1]);
                  }
                }
              }

              else {
                if (occupancy < state_process->sojourn_time[*pstate]->nb_value) {
                  if (j + occupancy < seq.length[i]) {
                    proba = state_process->sojourn_time[*pstate]->mass[occupancy];
                  }
                  else {
                    proba = (1. - state_process->sojourn_time[*pstate]->cumul[occupancy - 1]);
                  }
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
          }

          else {
            occupancy = seq.length[i] - j;
          }

          if (nb_output_process > 0) {
            for (k = j;k < j + occupancy;k++) {
              for (m = 0;m < nb_output_process;m++) {
                if (categorical_process[m]) {
                  proba = categorical_process[m]->observation[*pstate]->mass[*pioutput[m]];
                }

                else if (discrete_parametric_process[m]) {
                  proba = discrete_parametric_process[m]->observation[*pstate]->mass[*pioutput[m]];
                }

                else {
                  if (((continuous_parametric_process[m]->ident == GAMMA) ||
                       (continuous_parametric_process[m]->ident == ZERO_INFLATED_GAMMA)) && (seq.min_value[m + 1] < seq.min_interval[m + 1] / 2)) {
                    switch (seq.type[m + 1]) {
                    case INT_VALUE :
                      proba = continuous_parametric_process[m]->observation[*pstate]->mass_computation(*pioutput[m] , *pioutput[m] + seq.min_interval[m + 1]);
                      break;
                    case REAL_VALUE :
                      proba = continuous_parametric_process[m]->observation[*pstate]->mass_computation(*proutput[m] , *proutput[m] + seq.min_interval[m + 1]);
                      break;
                    }
                  }

                  else if (continuous_parametric_process[m]->ident == LINEAR_MODEL) {
                    switch (seq.type[m + 1]) {
                    case INT_VALUE :
                      residual = *pioutput[m] - (continuous_parametric_process[m]->observation[*pstate]->intercept +
                                  continuous_parametric_process[m]->observation[*pstate]->slope *
                                  (seq.index_param_type == IMPLICIT_TYPE ? k : seq.index_parameter[i][k]));
                      break;
                    case REAL_VALUE :
                      residual = *proutput[m] - (continuous_parametric_process[m]->observation[*pstate]->intercept +
                                  continuous_parametric_process[m]->observation[*pstate]->slope *
                                  (seq.index_param_type == IMPLICIT_TYPE ? k : seq.index_parameter[i][k]));
                      break;
                    }

                    proba = continuous_parametric_process[m]->observation[*pstate]->mass_computation(residual - seq.min_interval[m + 1] / 2 , residual + seq.min_interval[m + 1] / 2);
                  }

                  else if (continuous_parametric_process[m]->ident == AUTOREGRESSIVE_MODEL) {
                    if (k == 0) {
                      switch (seq.type[m + 1]) {
                      case INT_VALUE :
                        residual = *pioutput[m] - continuous_parametric_process[m]->observation[*pstate]->location;
                        break;
                      case REAL_VALUE :
                        residual = *proutput[m] - continuous_parametric_process[m]->observation[*pstate]->location;
                        break;
                      }
                    }

                    else {
                      switch (seq.type[m + 1]) {
                      case INT_VALUE :
                        residual = *pioutput[m] - (continuous_parametric_process[m]->observation[*pstate]->location +
                                    continuous_parametric_process[m]->observation[*pstate]->autoregressive_coeff *
                                    (*(pioutput[m] - 1) - continuous_parametric_process[m]->observation[*pstate]->location));
                        break;
                      case REAL_VALUE :
                        residual = *proutput[m] - (continuous_parametric_process[m]->observation[*pstate]->location +
                                    continuous_parametric_process[m]->observation[*pstate]->autoregressive_coeff *
                                    (*(proutput[m] - 1) - continuous_parametric_process[m]->observation[*pstate]->location));
                        break;
                      }
                    }

                    proba = continuous_parametric_process[m]->observation[*pstate]->mass_computation(residual - seq.min_interval[m + 1] / 2 , residual + seq.min_interval[m + 1] / 2);
                  }

                  else {
                    switch (seq.type[m + 1]) {
                    case INT_VALUE :
                      proba = continuous_parametric_process[m]->observation[*pstate]->mass_computation(*pioutput[m] - seq.min_interval[m + 1] / 2 , *pioutput[m] + seq.min_interval[m + 1] / 2);
                      break;
                    case REAL_VALUE :
                      proba = continuous_parametric_process[m]->observation[*pstate]->mass_computation(*proutput[m] - seq.min_interval[m + 1] / 2 , *proutput[m] + seq.min_interval[m + 1] / 2);
                      break;
                    }
                  }
                }

                switch (seq.type[m + 1]) {
                case INT_VALUE :
                  pioutput[m]++;
                  break;
                case REAL_VALUE :
                  proutput[m]++;
                  break;
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
            }

            if (likelihood == D_INF) {
              break;
            }
          }

          j += occupancy;
        }
        while (j < seq.length[i]);

        if (likelihood == D_INF) {
          break;
        }
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
 *  \brief Computation of the log-likelihood of a semi-Markov chain for sequences.
 *
 *  \param[in] seq reference on a SemiMarkovData object.
 *
 *  \return        log-likelihood.
 */
/*--------------------------------------------------------------*/

double SemiMarkov::likelihood_computation(const SemiMarkovData &seq) const

{
  int i , j;
  int nb_value;
  double buff , likelihood = 0.;
  FrequencyDistribution **initial_run , **final_run , **single_run;


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

      else if ((continuous_parametric_process[i]->ident == LINEAR_MODEL) ||
               (continuous_parametric_process[i]->ident == AUTOREGRESSIVE_MODEL) || (!(seq.marginal_distribution[i + 1]))) {
        likelihood = D_INF;
        break;
      }
    }
  }

  else {
    likelihood = D_INF;
  }

  if (likelihood != D_INF) {
    likelihood = Chain::likelihood_computation(*(seq.chain_data));

    if (likelihood != D_INF) {
      if (type == EQUILIBRIUM) {

        // construction of the censored sojourn time frequency distributions

        initial_run = new FrequencyDistribution*[seq.marginal_distribution[0]->nb_value];
        for (i = 0;i < seq.marginal_distribution[0]->nb_value;i++) {
          initial_run[i] = new FrequencyDistribution(seq.max_length);
        }

        final_run = new FrequencyDistribution*[seq.marginal_distribution[0]->nb_value];
        for (i = 0;i < seq.marginal_distribution[0]->nb_value;i++) {
          final_run[i] = new FrequencyDistribution(seq.max_length);
        }

        single_run = new FrequencyDistribution*[seq.marginal_distribution[0]->nb_value];
        for (i = 0;i < seq.marginal_distribution[0]->nb_value;i++) {
          single_run[i] = new FrequencyDistribution(seq.max_length + 1);
        }

        // update of the censored sojourn time frequency distributions

        seq.censored_sojourn_time_frequency_distribution_computation(initial_run , final_run , single_run);
      }

      for (i = 0;i < nb_state;i++) {
        if (sojourn_type[i] == SEMI_MARKOVIAN) {
          buff = state_process->sojourn_time[i]->likelihood_computation(*(seq.characteristics[0]->sojourn_time[i]));

          if (buff != D_INF) {
            likelihood += buff;
          }
          else {
            likelihood = D_INF;
            break;
          }

          switch (type) {

          case ORDINARY : {
            buff = state_process->sojourn_time[i]->survivor_likelihood_computation(*(seq.characteristics[0]->final_run[i]));

            if (buff != D_INF) {
              likelihood += buff;
            }
            else {
              likelihood = D_INF;
            }
            break;
          }

          case EQUILIBRIUM : {
            buff = state_process->sojourn_time[i]->survivor_likelihood_computation(*(final_run[i]));

            if (buff != D_INF) {
              likelihood += buff;
              buff = forward[i]->likelihood_computation(*(initial_run[i]));

              if (buff != D_INF) {
                likelihood += buff;
                buff = forward[i]->survivor_likelihood_computation(*(single_run[i]));

                if (buff != D_INF) {
                  likelihood += buff;
                }
                else {
                  likelihood = D_INF;
                }
              }

              else {
                likelihood = D_INF;
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
        }
      }

      if (type == EQUILIBRIUM) {
        for (i = 0;i < seq.marginal_distribution[0]->nb_value;i++) {
          delete initial_run[i];
        }
        delete [] initial_run;

        for (i = 0;i < seq.marginal_distribution[0]->nb_value;i++) {
          delete final_run[i];
        }
        delete [] final_run;

        for (i = 0;i < seq.marginal_distribution[0]->nb_value;i++) {
          delete single_run[i];
        }
        delete [] single_run;
      }
    }

    if (likelihood != D_INF) {
      for (i = 0;i < nb_output_process;i++) {
        if (categorical_process[i]) {
          for (j = 0;j < nb_state;j++) {
            buff = categorical_process[i]->observation[j]->likelihood_computation(*(seq.observation_distribution[i + 1][j]));

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
            buff = discrete_parametric_process[i]->observation[j]->likelihood_computation(*(seq.observation_distribution[i + 1][j]));

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
            buff = continuous_parametric_process[i]->observation[j]->likelihood_computation(*(seq.observation_distribution[i + 1][j]) ,
                                                                                            (int)seq.min_interval[i]);

            if (buff != D_INF) {
              likelihood += buff;
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

  return likelihood;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Counting of initial states and transitions.
 *
 *  \param[in] chain_data reference on a ChainData object,
 *  \param[in] smarkov    flags on the self-transition probabilities.
 */
/*--------------------------------------------------------------*/

void MarkovianSequences::transition_count_computation(const ChainData &chain_data ,
                                                      const SemiMarkov *smarkov) const

{
  int i , j;
  int *pstate;


  for (i = 0;i < chain_data.nb_state;i++) {
    chain_data.initial[i] = 0;
  }

  for (i = 0;i < chain_data.nb_row;i++) {
    for (j = 0;j < chain_data.nb_state;j++) {
      chain_data.transition[i][j] = 0;
    }
  }

  // extraction of initial states and transitions

  for (i = 0;i < nb_sequence;i++) {
    pstate = int_sequence[i][0];
    (chain_data.initial[*pstate])++;

    for (j = 1;j < length[i];j++) {
      pstate++;
      (chain_data.transition[*(pstate - 1)][*pstate])++;
    }
  }

  if (smarkov) {
    for (i = 0;i < chain_data.nb_state;i++) {
      if (smarkov->sojourn_type[i] == SEMI_MARKOVIAN) {
        chain_data.transition[i][i] = 0;
      }
    }
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Construction of the initial state and transition counts.
 *
 *  \param[in] smarkov flags on the self-transition probabilities.
 */
/*--------------------------------------------------------------*/

void SemiMarkovData::build_transition_count(const SemiMarkov *smarkov)

{
  chain_data = new ChainData(ORDINARY , marginal_distribution[0]->nb_value ,
                             marginal_distribution[0]->nb_value);
  transition_count_computation(*chain_data , smarkov);
}


/*--------------------------------------------------------------*/
/**
 *  \brief Estimation of a semi-Markov chain.
 *
 *  \param[in] error          reference on a StatError object,
 *  \param[in] display        flag for displaying estimation intermediate results,
 *  \param[in] itype          process type (ORDINARY/EQUILIBRIUM),
 *  \param[in] estimator      estimator type for the reestimation of the state occupancy distribution
 *                            (complete or partial likelihood),
 *  \param[in] counting_flag  flag on the computation of the counting distributions,
 *  \param[in] nb_iter        number of iterations,
 *  \param[in] mean_estimator method for the computation of the state occupancy
 *                            distribution mean (equilibrium semi-Markov chain).
 *
 *  \return                   SemiMarkov object.
 */
/*--------------------------------------------------------------*/

SemiMarkov* MarkovianSequences::semi_markov_estimation(StatError &error , bool display , process_type itype ,
                                                       censoring_estimator estimator , bool counting_flag , int nb_iter ,
                                                       duration_distribution_mean_estimator mean_estimator) const

{
  bool status = true;
  int i , j;
  int nb_likelihood_decrease , *occupancy_survivor , *censored_occupancy_survivor , nb_value[1];
  double likelihood , previous_likelihood , hlikelihood , occupancy_mean;
  DiscreteParametric *occupancy;
  Forward *forward;
  Reestimation<double> *occupancy_reestim , *length_bias_reestim;
  SemiMarkov *smarkov;
  SemiMarkovData *seq;
  FrequencyDistribution *complete_run , *censored_run , *pfinal_run , *hreestim ,
                        **initial_run , **final_run , **single_run;
  const FrequencyDistribution *prun[3];


  smarkov = NULL;
  error.init();

  if ((type[0] != INT_VALUE) && (type[0] != STATE)) {
    status = false;
    ostringstream correction_message;
    correction_message << STAT_variable_word[INT_VALUE] << " or " << STAT_variable_word[STATE];
    error.correction_update(STAT_error[STATR_VARIABLE_TYPE] , (correction_message.str()).c_str());
  }

  else {
    if ((marginal_distribution[0]->nb_value < 2) ||
        (marginal_distribution[0]->nb_value > NB_OUTPUT)) {
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

  if (status) {
    if (max_length > COUNTING_MAX_LENGTH) {
      counting_flag = false;
    }

    if (itype == EQUILIBRIUM) {

      // construction of the censored sojourn time frequency distributions

      initial_run = new FrequencyDistribution*[marginal_distribution[0]->nb_value];
      for (i = 0;i < marginal_distribution[0]->nb_value;i++) {
        initial_run[i] = new FrequencyDistribution(max_length);
      }

      final_run = new FrequencyDistribution*[marginal_distribution[0]->nb_value];
      for (i = 0;i < marginal_distribution[0]->nb_value;i++) {
        final_run[i] = new FrequencyDistribution(max_length);
      }

      single_run = new FrequencyDistribution*[marginal_distribution[0]->nb_value];
      for (i = 0;i < marginal_distribution[0]->nb_value;i++) {
        single_run[i] = new FrequencyDistribution(max_length + 1);
      }

      // update of the censored sojourn time frequency distributions

      censored_sojourn_time_frequency_distribution_computation(initial_run , final_run , single_run);
    }

    if (nb_variable == 2) {
      nb_value[0] = marginal_distribution[1]->nb_value;
    }

    smarkov = new SemiMarkov(itype , marginal_distribution[0]->nb_value ,
                             nb_variable - 1 , nb_value);
    smarkov->semi_markov_data = new SemiMarkovData(*this , SEQUENCE_COPY ,
                                                   (itype == EQUILIBRIUM ? true : false));

    seq = smarkov->semi_markov_data;
    seq->state_variable_init();
    seq->build_transition_count();

    for (i = 0;i < smarkov->nb_state;i++) {
      if ((seq->characteristics[0]->sojourn_time[i]->nb_element > 0) ||
          ((itype == EQUILIBRIUM) && (initial_run[i]->nb_element > 0))) {
        seq->chain_data->transition[i][i] = 0;
      }
    }

    // estimation of the Markov chain parameters

    seq->chain_data->estimation(*smarkov);
    smarkov->component_computation();

    if ((itype == EQUILIBRIUM) && (smarkov->nb_component > 1)) {
      delete smarkov;
      smarkov = NULL;
      error.correction_update(STAT_parsing[STATP_CHAIN_STRUCTURE] , STAT_parsing[STATP_IRREDUCIBLE]);
    }

    else {

      // estimation of the state occupancy distributions

      if (estimator != PARTIAL_LIKELIHOOD) {
        occupancy_survivor = new int[max_length];
        censored_occupancy_survivor = new int[max_length + 1];
        occupancy_reestim = new Reestimation<double>(max_length + 1);
      }

      smarkov->sojourn_type = new state_sojourn_type[smarkov->nb_state];
      smarkov->state_process->absorption = new double[smarkov->nb_state];
      smarkov->state_process->sojourn_time = new DiscreteParametric*[smarkov->nb_state];
      smarkov->forward = new Forward*[smarkov->nb_state];
      for (i = 0;i < smarkov->nb_state;i++) {
        smarkov->forward[i] = NULL;
        smarkov->state_process->sojourn_time[i] = NULL;
      }

      for (i = 0;i < smarkov->nb_state;i++) {
        if ((seq->characteristics[0]->sojourn_time[i]->nb_element == 0) &&
            ((itype == ORDINARY) || (initial_run[i]->nb_element == 0))) {
          smarkov->sojourn_type[i] = MARKOVIAN;
          smarkov->state_process->absorption[i] = 1.;
        }

        else {
          smarkov->state_process->absorption[i] = 0.;

          if ((itype == EQUILIBRIUM) && (seq->characteristics[0]->sojourn_time[i]->nb_element == 0)) {
            occupancy = NULL;
          }

          else if ((estimator != PARTIAL_LIKELIHOOD) && (itype == EQUILIBRIUM) && 
                   ((initial_run[i]->nb_element > 0) || (single_run[i]->nb_element > 0))) {

            //  initialization of the state occupancy distribution

            prun[0] = seq->characteristics[0]->sojourn_time[i];
            prun[1] = seq->characteristics[0]->sojourn_time[i];
            complete_run = new FrequencyDistribution(2 , prun);

//            prun[0] = seq->initial_run[0][i];
//            prun[1] = final_run[i];

//            prun[0] = initial_run[i];
//            prun[1] = seq->characteristics[0]->final_run[i];
//            censored_run = new FrequencyDistribution(2 , prun);

            prun[0] = initial_run[i];
            prun[1] = final_run[i];
            prun[2] = single_run[i];
            censored_run = new FrequencyDistribution(3 , prun);

            complete_run->state_occupancy_estimation(censored_run , occupancy_reestim ,
                                                     occupancy_survivor ,
                                                     censored_occupancy_survivor);
            delete complete_run;
            delete censored_run;

            if ((estimator == KAPLAN_MEIER) && (single_run[i]->nb_element == 0)) {
              occupancy = occupancy_reestim->type_parametric_estimation(1 , true , OCCUPANCY_THRESHOLD);
            }

            else {
              occupancy = new DiscreteParametric(NEGATIVE_BINOMIAL , 1 , I_DEFAULT , 1. ,
                                                 (occupancy_reestim->mean > 1. ? 1. / occupancy_reestim->mean : 0.99) ,
                                                 OCCUPANCY_THRESHOLD);
              occupancy->init(CATEGORICAL , I_DEFAULT , I_DEFAULT , D_DEFAULT , D_DEFAULT);
              forward = new Forward(*occupancy , occupancy->alloc_nb_value);

              delete occupancy_reestim;
              occupancy_reestim = new Reestimation<double>(MAX(occupancy->alloc_nb_value , max_length + 1));
              length_bias_reestim = new Reestimation<double>(MAX(occupancy->alloc_nb_value , max_length + 1));

              likelihood = D_INF;
              j = 0;

              do {
                j++;

                // computation of the reestimation quantities of the state occupancy distribution

                occupancy->expectation_step(*(seq->characteristics[0]->sojourn_time[i]) ,
                                            *(initial_run[i]) , *(final_run[i]) ,
                                            *(single_run[i]) , occupancy_reestim ,
                                            length_bias_reestim , j);

                switch (mean_estimator) {
                case COMPUTED :
                  occupancy_mean = interval_bisection(occupancy_reestim , length_bias_reestim);
                  break;
                case ONE_STEP_LATE :
                  occupancy_mean = occupancy->mean;
                  break;
                }

                occupancy_reestim->equilibrium_process_estimation(length_bias_reestim , occupancy ,
                                                                  occupancy_mean);
                forward->computation(*occupancy);

                previous_likelihood = likelihood;
                likelihood = occupancy->state_occupancy_likelihood_computation(*forward , *(seq->characteristics[0]->sojourn_time[i]) ,
                                                                               *(initial_run[i]) , *(final_run[i]) ,
                                                                               *(single_run[i]));

                if ((display) && ((j < 10) || ((j < 100) && (j % 10 == 0)) || ((j < 1000) && (j % 100 == 0)) || (j % 1000 == 0))) {
                  cout << STAT_label[STATL_ITERATION] << " " << j << "   "
                       << STAT_label[STATL_LIKELIHOOD] << ": " << likelihood << "   "
                       << STAT_label[STATL_SMOOTHNESS] << ": " << occupancy->second_difference_norm_computation() << endl;
                }
              }
              while ((likelihood != D_INF) && (((nb_iter == I_DEFAULT) && (j < OCCUPANCY_NB_ITER) &&
                       ((likelihood - previous_likelihood) / -likelihood > OCCUPANCY_LIKELIHOOD_DIFF)) ||
                      ((nb_iter != I_DEFAULT) && (j < nb_iter))));

              if (likelihood != D_INF) {
                if (display) {
                  cout << "\n" << STAT_label[STATL_STATE] << " " << i << " - "
                       << j << " " << STAT_label[STATL_ITERATIONS] << "   "
                       << STAT_label[STATL_LIKELIHOOD] << ": " << likelihood << "   "
                       << STAT_label[STATL_SMOOTHNESS] << ": " << occupancy->second_difference_norm_computation() << "\n" << endl;
                }

                hreestim = new FrequencyDistribution(MAX(occupancy->alloc_nb_value , max_length + 1));

                likelihood = D_INF;
                nb_likelihood_decrease = 0;

                j = 0;
                do {
                  j++;

                  // computation of the reestimation quantities of the state occupancy distribution

                  occupancy->expectation_step(*(seq->characteristics[0]->sojourn_time[i]) ,
                                              *(initial_run[i]) , *(final_run[i]) ,
                                              *(single_run[i]) , occupancy_reestim ,
                                              length_bias_reestim , j , true , mean_estimator);

                  hreestim->update(occupancy_reestim , (int)(occupancy_reestim->nb_element *
                                   MAX(sqrt(occupancy_reestim->variance) , 1.) * OCCUPANCY_COEFF));
                  hlikelihood = hreestim->Reestimation<int>::type_parametric_estimation(occupancy , 1 , true ,
                                                                                        OCCUPANCY_THRESHOLD);

                  if (hlikelihood == D_INF) {
                    likelihood = D_INF;
                  }

                  else {
                    occupancy->computation(hreestim->nb_value , OCCUPANCY_THRESHOLD);
                    forward->copy(*occupancy);
                    forward->computation(*occupancy);

                    previous_likelihood = likelihood;
                    likelihood = occupancy->state_occupancy_likelihood_computation(*forward , *(seq->characteristics[0]->sojourn_time[i]) ,
                                                                                   *(initial_run[i]) , *(final_run[i]) ,
                                                                                   *(single_run[i]));

                    if (likelihood < previous_likelihood) {
                      nb_likelihood_decrease++;
                    }
                    else {
                      nb_likelihood_decrease = 0;
                    }

                    if ((display) && ((j < 10) || ((j < 100) && (j % 10 == 0)) || ((j < 1000) && (j % 100 == 0)) || (j % 1000 == 0))) {
                      cout << STAT_label[STATL_ITERATION] << " " << j << "   "
                           << STAT_label[STATL_LIKELIHOOD] << ": " << likelihood << "   "
                           << STAT_label[STATL_SMOOTHNESS] << ": " << occupancy->second_difference_norm_computation() << endl;
                    }
                  }
                }
                while ((likelihood != D_INF) && (j < OCCUPANCY_NB_ITER) &&
                       (((likelihood - previous_likelihood) / -likelihood > OCCUPANCY_LIKELIHOOD_DIFF) ||
                        (hlikelihood == D_INF) || (nb_likelihood_decrease == 1)));

                delete hreestim;

                if (likelihood != D_INF) {
                  if (display) {
                    cout << "\n" << STAT_label[STATL_STATE] << " " << i << " - "
                         << j << " " << STAT_label[STATL_ITERATIONS] << "   "
                         << STAT_label[STATL_LIKELIHOOD] << ": " << likelihood << "   "
                         << STAT_label[STATL_SMOOTHNESS] << ": " << occupancy->second_difference_norm_computation() << "\n" << endl;
                  }
                }

                else {
                  delete occupancy;
                  occupancy = NULL;
                }
              }

              else {
                delete occupancy;
                occupancy = NULL;
              }

              delete forward;
              delete length_bias_reestim;
            }
          }

          else if ((estimator != PARTIAL_LIKELIHOOD) && (((itype == ORDINARY) && (seq->characteristics[0]->final_run[i]->nb_element > 0) &&
                     (seq->characteristics[0]->final_run[i]->nb_value > seq->characteristics[0]->sojourn_time[i]->nb_value)) || ((itype == EQUILIBRIUM) &&
                     (final_run[i]->nb_element > 0) && (final_run[i]->nb_value > seq->characteristics[0]->sojourn_time[i]->nb_value)))) {
            switch (itype) {
            case ORDINARY :
              pfinal_run = seq->characteristics[0]->final_run[i];
              break;
            case EQUILIBRIUM :
              pfinal_run = final_run[i];
              break;
            }

            //  initialization of the state occupancy distribution

            seq->characteristics[0]->sojourn_time[i]->state_occupancy_estimation(pfinal_run , occupancy_reestim ,
                                                                                 occupancy_survivor ,
                                                                                 censored_occupancy_survivor);
            occupancy = new DiscreteParametric(NEGATIVE_BINOMIAL , 1 , I_DEFAULT , 1. ,
                                               1. / occupancy_reestim->mean , OCCUPANCY_THRESHOLD);
            occupancy->init(CATEGORICAL , I_DEFAULT , I_DEFAULT , D_DEFAULT , D_DEFAULT);

            delete occupancy_reestim;
            occupancy_reestim = new Reestimation<double>(MAX(occupancy->alloc_nb_value , max_length + 1));

            likelihood = D_INF;
            j = 0;

            do {
              j++;

              // computation of the reestimation quantities of the state occupancy distribution

              occupancy->expectation_step(*(seq->characteristics[0]->sojourn_time[i]) ,
                                          *pfinal_run , occupancy_reestim , j);
              occupancy_reestim->distribution_estimation(occupancy);

              previous_likelihood = likelihood;
              likelihood = occupancy->state_occupancy_likelihood_computation(*(seq->characteristics[0]->sojourn_time[i]) ,
                                                                             *pfinal_run);

              if ((display) && ((j < 10) || ((j < 100) && (j % 10 == 0)) || ((j < 1000) && (j % 100 == 0)) || (j % 1000 == 0))) {
                cout << STAT_label[STATL_ITERATION] << " " << j << "   "
                     << STAT_label[STATL_LIKELIHOOD] << ": " << likelihood << "   "
                     << STAT_label[STATL_SMOOTHNESS] << ": " << occupancy->second_difference_norm_computation() << endl;
              }
            }
            while ((likelihood != D_INF) && (((nb_iter == I_DEFAULT) && (j < OCCUPANCY_NB_ITER) &&
                     ((likelihood - previous_likelihood) / -likelihood > OCCUPANCY_LIKELIHOOD_DIFF)) ||
                    ((nb_iter != I_DEFAULT) && (j < nb_iter))));

            if (likelihood != D_INF) {
              if (display) {
                cout << "\n" << STAT_label[STATL_STATE] << " " << i << " - "
                     << j << " " << STAT_label[STATL_ITERATIONS] << "   "
                     << STAT_label[STATL_LIKELIHOOD] << ": " << likelihood << "   "
                     << STAT_label[STATL_SMOOTHNESS] << ": " << occupancy->second_difference_norm_computation() << "\n" << endl;
              }

              hreestim = new FrequencyDistribution(MAX(occupancy->alloc_nb_value , max_length + 1));

              likelihood = D_INF;
              nb_likelihood_decrease = 0;

              j = 0;
              do {
                j++;

                // computation of the reestimation quantities of the state occupancy distribution

                occupancy->expectation_step(*(seq->characteristics[0]->sojourn_time[i]) ,
                                            *pfinal_run , occupancy_reestim , j);

                hreestim->update(occupancy_reestim , (int)(occupancy_reestim->nb_element *
                                 MAX(sqrt(occupancy_reestim->variance) , 1.) * OCCUPANCY_COEFF));
                hlikelihood = hreestim->Reestimation<int>::type_parametric_estimation(occupancy , 1 , true ,
                                                                                      OCCUPANCY_THRESHOLD);

                if (hlikelihood == D_INF) {
                  likelihood = D_INF;
                }

                else {
                  occupancy->computation(hreestim->nb_value , OCCUPANCY_THRESHOLD);

                  previous_likelihood = likelihood;
                  likelihood = occupancy->state_occupancy_likelihood_computation(*(seq->characteristics[0]->sojourn_time[i]) ,
                                                                                 *pfinal_run);

                  if (likelihood < previous_likelihood) {
                    nb_likelihood_decrease++;
                  }
                  else {
                    nb_likelihood_decrease = 0;
                  }

                  if ((display) && ((j < 10) || ((j < 100) && (j % 10 == 0)) || ((j < 1000) && (j % 100 == 0)) || (j % 1000 == 0))) {
                    cout << STAT_label[STATL_ITERATION] << " " << j << "   "
                         << STAT_label[STATL_LIKELIHOOD] << ": " << likelihood << "   "
                         << STAT_label[STATL_SMOOTHNESS] << ": " << occupancy->second_difference_norm_computation() << endl;
                  }
                }
              }
              while ((likelihood != D_INF) && (j < OCCUPANCY_NB_ITER) &&
                     (((likelihood - previous_likelihood) / -likelihood > OCCUPANCY_LIKELIHOOD_DIFF) ||
                      (hlikelihood == D_INF) || (nb_likelihood_decrease == 1)));

              delete hreestim;

              if (likelihood != D_INF) {
                if (display) {
                  cout << "\n" << STAT_label[STATL_STATE] << " " << i << " - "
                       << j << " " << STAT_label[STATL_ITERATIONS] << "   "
                       << STAT_label[STATL_LIKELIHOOD] << ": " << likelihood << "   "
                       << STAT_label[STATL_SMOOTHNESS] << ": " << occupancy->second_difference_norm_computation() << "\n" << endl;
                }
              }

              else {
                delete occupancy;
                occupancy = NULL;
              }
            }

            else {
              delete occupancy;
              occupancy = NULL;
            }
          }

          else if ((estimator != PARTIAL_LIKELIHOOD) && (((itype == ORDINARY) && (seq->characteristics[0]->final_run[i]->nb_element > 0)) ||
                    ((itype == EQUILIBRIUM) && (final_run[i]->nb_element > 0)))) {
            seq->characteristics[0]->sojourn_time[i]->state_occupancy_estimation((itype == ORDINARY ? seq->characteristics[0]->final_run[i] : final_run[i]) ,
                                                                                 occupancy_reestim , occupancy_survivor ,
                                                                                 censored_occupancy_survivor);

            occupancy = occupancy_reestim->type_parametric_estimation(1 , true , OCCUPANCY_THRESHOLD);

/*          occupancy = new DiscreteParametric(occupancy_reestim->nb_value);
            occupancy_reestim->distribution_estimation(occupancy); */
          }

          else {
            occupancy = seq->characteristics[0]->sojourn_time[i]->Reestimation<int>::type_parametric_estimation(1 , true ,
                                                                                                                OCCUPANCY_THRESHOLD);
          }

          if (occupancy) {
            if (occupancy->mean == 1.) {
              smarkov->sojourn_type[i] = MARKOVIAN;
            }

            else {
              smarkov->sojourn_type[i] = SEMI_MARKOVIAN;
              smarkov->state_process->sojourn_time[i] = new DiscreteParametric(*occupancy);
              if (smarkov->stype[i] == RECURRENT) {
                smarkov->forward[i] = new Forward(*(smarkov->state_process->sojourn_time[i]));
              }
            }

            delete occupancy;
          }

          else {
            delete smarkov;
            smarkov = NULL;

            ostringstream error_message;
            error_message << STAT_label[STATL_STATE] << " " << i << " "
                          << SEQ_label[SEQL_OCCUPANCY_DISTRIBUTION] << " "
                          << STAT_error[STATR_ESTIMATION_FAILURE];
            error.update((error_message.str()).c_str());
            break;
          }
        }
      }

      if (estimator != PARTIAL_LIKELIHOOD) {
        delete [] occupancy_survivor;
        delete [] censored_occupancy_survivor;
        delete occupancy_reestim;
      }
    }

    if (itype == EQUILIBRIUM) {
      for (i = 0;i < marginal_distribution[0]->nb_value;i++) {
        delete initial_run[i];
      }
      delete [] initial_run;

      for (i = 0;i < marginal_distribution[0]->nb_value;i++) {
        delete final_run[i];
      }
      delete [] final_run;

      for (i = 0;i < marginal_distribution[0]->nb_value;i++) {
        delete single_run[i];
      }
      delete [] single_run;
    }

    if (smarkov) {
      if (itype == EQUILIBRIUM) {
        for (i = 0;i < smarkov->nb_state;i++) {
          smarkov->initial[i] = 1. / (double)smarkov->nb_state;
        }
        smarkov->initial_probability_computation();
      }

      // estimation of categorical observation distributions

      if (smarkov->nb_output_process == 1) {
        seq->build_observation_frequency_distribution(smarkov->nb_state);

        for (i = 0;i < smarkov->nb_state;i++) {
          seq->observation_distribution[1][i]->distribution_estimation(smarkov->categorical_process[0]->observation[i]);
        }
      }

      // computation of the log-likelihood and the characteristic distributions of the estimated semi-Markov chain

      seq->likelihood = smarkov->likelihood_computation(*seq);

#     ifdef DEBUG
      cout << "\n" << STAT_label[STATL_LIKELIHOOD] << ": " << seq->likelihood << " | "
           << smarkov->likelihood_computation(*seq , I_DEFAULT) << endl;
#     endif

      if (seq->likelihood == D_INF) {
        delete smarkov;
        smarkov = NULL;
        error.update(STAT_error[STATR_ESTIMATION_FAILURE]);
      }

      else {
        smarkov->characteristic_computation(*seq , counting_flag , I_DEFAULT , false);
      }
    }
  }

  return(smarkov);
}


/*--------------------------------------------------------------*/
/**
 *  \brief Comparison of semi-Markov chains for a sample of sequences.
 *
 *  \param[in] error    reference on a StatError object,
 *  \param[in] display  flag for displaying the results of model comparison,
 *  \param[in] nb_model number of semi-Markov chains,
 *  \param[in] ismarkov pointer on SemiMarkov objects,
 *  \param[in] path     file path.
 *
 *  \return             error status.
 */
/*--------------------------------------------------------------*/

bool MarkovianSequences::comparison(StatError &error , bool display , int nb_model ,
                                    const SemiMarkov **ismarkov , const string path) const

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
    if (ismarkov[i]->nb_output_process + 1 != nb_variable) {
      status = false;
      ostringstream error_message;
      error_message << SEQ_label[SEQL_SEMI_MARKOV_CHAIN] << " " << i + 1 << ": "
                    << STAT_error[STATR_NB_OUTPUT_PROCESS];
      error.update((error_message.str()).c_str());
    }

    else {
      if (ismarkov[i]->state_process->nb_value < marginal_distribution[0]->nb_value) {
        status = false;
        ostringstream error_message;
        error_message << SEQ_label[SEQL_SEMI_MARKOV_CHAIN] << " " << i + 1 << ": "
                      << SEQ_error[SEQR_NB_STATE];
        error.update((error_message.str()).c_str());
      }

      if (nb_variable == 2) {
        if (ismarkov[i]->categorical_process[0]->nb_value < marginal_distribution[1]->nb_value) {
          status = false;
          ostringstream error_message;
          error_message << SEQ_label[SEQL_SEMI_MARKOV_CHAIN] << " " << i + 1 << ": "
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
        likelihood[i][j] = ismarkov[j]->likelihood_computation(*this , i);
      }
    }

    if (display) {
      likelihood_write(cout , nb_model , likelihood , SEQ_label[SEQL_SEMI_MARKOV_CHAIN] , true);
    }
    if (!path.empty()) {
      status = likelihood_write(error , path , nb_model , likelihood , SEQ_label[SEQL_SEMI_MARKOV_CHAIN]);
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
 *  \brief Simulation using a semi-Markov chain.
 *
 *  \param[in] error               reference on a StatError object,
 *  \param[in] length_distribution sequence length frequency distribution,
 *  \param[in] counting_flag       flag on the computation of the counting distributions,
 *  \param[in] divergence_flag     flag on the computation of a Kullback-Leibler divergence.
 *
 *  \return                        SemiMarkovData object.
 */
/*--------------------------------------------------------------*/

SemiMarkovData* SemiMarkov::simulation(StatError &error , const FrequencyDistribution &length_distribution ,
                                       bool counting_flag , bool divergence_flag) const

{
  bool status = true , hidden;
  int i , j , k , m;
  int cumul_length , occupancy , *decimal_scale , *pstate , **pioutput;
  variable_nature *itype;
  double buff , min_location , likelihood , **proutput;
  Distribution *weight , *restoration_weight;
  SemiMarkov *smarkov;
  SemiMarkovData *seq;


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

    seq = new SemiMarkovData(length_distribution , nb_output_process + 1 , itype);
    delete [] itype;

    seq->semi_markov = new SemiMarkov(*this , false);

    smarkov = seq->semi_markov;
    smarkov->create_cumul();
    smarkov->cumul_computation();

    if (smarkov->nb_output_process > 0) {
      pioutput = new int*[smarkov->nb_output_process];
      proutput = new double*[smarkov->nb_output_process];

      decimal_scale = new int[smarkov->nb_output_process];

      for (i = 0;i < smarkov->nb_output_process;i++) {
        if (smarkov->continuous_parametric_process[i]) {
          switch (smarkov->continuous_parametric_process[i]->ident) {

          case GAMMA : {
            min_location = smarkov->continuous_parametric_process[i]->observation[0]->location * smarkov->continuous_parametric_process[i]->observation[0]->dispersion;
            for (j = 1;j < smarkov->nb_state;j++) {
              buff = smarkov->continuous_parametric_process[i]->observation[j]->location * smarkov->continuous_parametric_process[i]->observation[0]->dispersion;
              if (buff < min_location) {
                min_location = buff;
              }
            }

            buff = (int)ceil(log(min_location) / log(10));
            if (buff < GAMMA_MAX_NB_DECIMAL) {
              decimal_scale[i] = pow(10 , (GAMMA_MAX_NB_DECIMAL - buff));
            }
            else {
              decimal_scale[i] = 1;
            }

#           ifdef MESSAGE
            cout << "\nScale: " << i + 1 << " " << decimal_scale[i] << endl;
#           endif

            break;
          }

          case GAUSSIAN : {
            min_location = fabs(smarkov->continuous_parametric_process[i]->observation[0]->location);
            for (j = 1;j < smarkov->nb_state;j++) {
              buff = fabs(smarkov->continuous_parametric_process[i]->observation[j]->location);
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
            switch (smarkov->continuous_parametric_process[i]->unit) {
            case DEGREE :
              decimal_scale[i] = DEGREE_DECIMAL_SCALE;
              break;
            case RADIAN :
              decimal_scale[i] = RADIAN_DECIMAL_SCALE;
              break;
            }

            for (j = 0;j < smarkov->nb_state;j++) {
              smarkov->continuous_parametric_process[i]->observation[j]->von_mises_cumul_computation();
            }
            break;
          }
          }
        }
      }
    }

    for (i = 0;i < seq->nb_sequence;i++) {
      pstate = seq->int_sequence[i][0];
      *pstate = cumul_method(smarkov->nb_state , smarkov->cumul_initial);

      for (j = 0;j < smarkov->nb_output_process;j++) {
        switch (seq->type[j + 1]) {
        case INT_VALUE :
          pioutput[j] = seq->int_sequence[i][j + 1];
          break;
        case REAL_VALUE :
          proutput[j] = seq->real_sequence[i][j + 1];
          break;
        }
      }

      j = 0;
      do {
        if (j > 0) {
          pstate++;
          *pstate = cumul_method(smarkov->nb_state , smarkov->cumul_transition[*(pstate - 1)]);
        }

        switch (smarkov->sojourn_type[*pstate]) {

        case SEMI_MARKOVIAN : {
          if ((smarkov->type == EQUILIBRIUM) && (j == 0)) {
            occupancy = smarkov->forward[*pstate]->simulation();
          }
          else {
            occupancy = smarkov->state_process->sojourn_time[*pstate]->simulation();
          }

          if (j + occupancy > seq->length[i]) {
            occupancy = seq->length[i] - j;
          }
          break;
        }

        case MARKOVIAN : {
          if (smarkov->transition[*pstate][*pstate] < 1.) {
            occupancy = 1;
          }
          else {
            occupancy = seq->length[i] - j;
          }
          break;
        }
        }

        for (k = 1;k < occupancy;k++) {
          pstate++;
          *pstate = *(pstate - 1);
        }

        for (k = j;k < j + occupancy;k++) {
          for (m = 0;m < smarkov->nb_output_process;m++) {
            if (smarkov->categorical_process[m]) {
              *pioutput[m] = smarkov->categorical_process[m]->observation[*pstate]->simulation();
            }

            else if (smarkov->discrete_parametric_process[m]) {
              *pioutput[m] = smarkov->discrete_parametric_process[m]->observation[*pstate]->simulation();
            }

            else {
              if (smarkov->continuous_parametric_process[m]->ident == LINEAR_MODEL) {
                *proutput[m] = smarkov->continuous_parametric_process[m]->observation[*pstate]->intercept +
                               smarkov->continuous_parametric_process[m]->observation[*pstate]->slope * k +
//                               round(smarkov->continuous_parametric_process[m]->observation[*pstate]->simulation() * decimal_scale[m]) / decimal_scale[m];
                               smarkov->continuous_parametric_process[m]->observation[*pstate]->simulation();
              }

              else if (smarkov->continuous_parametric_process[m]->ident == AUTOREGRESSIVE_MODEL) {
                if (k == 0) {
                  *proutput[m] = smarkov->continuous_parametric_process[m]->observation[*pstate]->location +
//                                 round(smarkov->continuous_parametric_process[m]->observation[*pstate]->simulation() * decimal_scale[m]) / decimal_scale[m];
                                 smarkov->continuous_parametric_process[m]->observation[*pstate]->simulation();
                }
                else {
                  *proutput[m] = smarkov->continuous_parametric_process[m]->observation[*pstate]->location +
                                 smarkov->continuous_parametric_process[m]->observation[*pstate]->autoregressive_coeff *
                                 (*(proutput[m] - 1) - smarkov->continuous_parametric_process[m]->observation[*pstate]->location) +
//                                 round(smarkov->continuous_parametric_process[m]->observation[*pstate]->simulation() * decimal_scale[m]) / decimal_scale[m];
                                 smarkov->continuous_parametric_process[m]->observation[*pstate]->simulation();
                }
              }

              else {
                *proutput[m] = round(smarkov->continuous_parametric_process[m]->observation[*pstate]->simulation() * decimal_scale[m]) / decimal_scale[m];
              }
            }

            switch (seq->type[m + 1]) {
            case INT_VALUE :
              pioutput[m]++;
              break;
            case REAL_VALUE :
              proutput[m]++;
              break;
            }
          }
        }

        j += occupancy;
      }
      while (j < seq->length[i]);
    }

    smarkov->remove_cumul();

    if (smarkov->nb_output_process > 0) {
      delete [] pioutput;
      delete [] proutput;

      delete [] decimal_scale;

      for (i = 0;i < smarkov->nb_output_process;i++) {
        if ((smarkov->continuous_parametric_process[i]) &&
            (smarkov->continuous_parametric_process[i]->ident == VON_MISES)) {
          for (j = 0;j < smarkov->nb_state;j++) {
            delete [] smarkov->continuous_parametric_process[i]->observation[j]->cumul;
            smarkov->continuous_parametric_process[i]->observation[j]->cumul = NULL;
          }
        }
      }
    }

    // extraction of the characteristics of the generated sequences

    seq->min_value[0] = 0;
    seq->max_value[0] = nb_state - 1;
    seq->build_marginal_frequency_distribution(0);

    for (i = 1;i < seq->nb_variable;i++) {
      seq->min_value_computation(i);
      seq->max_value_computation(i);

      seq->build_marginal_frequency_distribution(i);
      seq->min_interval_computation(i);
    }

    seq->build_transition_count(smarkov);
    seq->build_observation_frequency_distribution(nb_state);
    seq->build_observation_histogram(nb_state);
    seq->build_characteristic(I_DEFAULT , true , (type == EQUILIBRIUM ? true : false));

/*    if ((seq->max_value[0] < nb_state - 1) || (!(seq->characteristics[0]))) {
      delete seq;
      seq = NULL;
      error.update(SEQ_error[SEQR_STATES_NOT_REPRESENTED]);
    }

    else if (!divergence_flag) { */
    if (!divergence_flag) {
      smarkov->characteristic_computation(*seq , counting_flag);

      // computation of the log-likelihood of the model for the generated sequences

      likelihood = smarkov->likelihood_computation(*seq);

      if (likelihood == D_INF) {
        likelihood = smarkov->likelihood_computation(*seq , I_DEFAULT);
      }

#     ifdef DEBUG
      else {
        cout << "\n" << STAT_label[STATL_LIKELIHOOD] << ": " << likelihood
             << " | " << smarkov->likelihood_computation(*seq , I_DEFAULT) << endl;
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
        weight = NULL;
        restoration_weight = NULL;

        for (i = 0;i < smarkov->nb_output_process;i++) {
          if ((smarkov->categorical_process[i]) || (smarkov->discrete_parametric_process[i]) ||
              ((smarkov->continuous_parametric_process[i]) &&
               (smarkov->continuous_parametric_process[i]->ident != LINEAR_MODEL))) {
            weight = smarkov->state_process->weight_computation();
            restoration_weight = seq->weight_computation();
            break;
          }
        }

        for (i = 0;i < smarkov->nb_output_process;i++) {
          if (smarkov->categorical_process[i]) {
            delete smarkov->categorical_process[i]->weight;
            delete smarkov->categorical_process[i]->mixture;
            smarkov->categorical_process[i]->weight = new Distribution(*weight);
            smarkov->categorical_process[i]->mixture = smarkov->categorical_process[i]->mixture_computation(smarkov->categorical_process[i]->weight);
            delete smarkov->categorical_process[i]->restoration_weight;
            delete smarkov->categorical_process[i]->restoration_mixture;
            smarkov->categorical_process[i]->restoration_weight = new Distribution(*restoration_weight);
            smarkov->categorical_process[i]->restoration_mixture = smarkov->categorical_process[i]->mixture_computation(smarkov->categorical_process[i]->restoration_weight);
          }

          else if (smarkov->discrete_parametric_process[i]) {
            delete smarkov->discrete_parametric_process[i]->weight;
            delete smarkov->discrete_parametric_process[i]->mixture;
            smarkov->discrete_parametric_process[i]->weight = new Distribution(*weight);
            smarkov->discrete_parametric_process[i]->mixture = smarkov->discrete_parametric_process[i]->mixture_computation(smarkov->discrete_parametric_process[i]->weight);

            delete smarkov->discrete_parametric_process[i]->restoration_weight;
            delete smarkov->discrete_parametric_process[i]->restoration_mixture;
            smarkov->discrete_parametric_process[i]->restoration_weight = new Distribution(*restoration_weight);
            smarkov->discrete_parametric_process[i]->restoration_mixture = smarkov->discrete_parametric_process[i]->mixture_computation(smarkov->discrete_parametric_process[i]->restoration_weight);
          }

          else if ((smarkov->continuous_parametric_process[i]) &&
                   (smarkov->continuous_parametric_process[i]->ident != LINEAR_MODEL)) {
            delete smarkov->continuous_parametric_process[i]->weight;
            smarkov->continuous_parametric_process[i]->weight = new Distribution(*weight);

            delete smarkov->continuous_parametric_process[i]->restoration_weight;
            smarkov->continuous_parametric_process[i]->restoration_weight = new Distribution(*restoration_weight);
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
 *  \brief Simulation using a semi-Markov chain.
 *
 *  \param[in] error         reference on a StatError object,
 *  \param[in] nb_sequence   number of sequences,
 *  \param[in] length        sequence length,
 *  \param[in] counting_flag flag on the computation of the counting distributions.
 *
 *  \return                  SemiMarkovData object.
 */
/*--------------------------------------------------------------*/

SemiMarkovData* SemiMarkov::simulation(StatError &error , int nb_sequence ,
                                       int length , bool counting_flag) const

{
  bool status = true;
  SemiMarkovData *seq;


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
 *  \brief Simulation using a semi-Markov chain.
 *
 *  \param[in] error         reference on a StatError object,
 *  \param[in] nb_sequence   number of sequences,
 *  \param[in] iseq          reference on a MarkovianSequences object,
 *  \param[in] counting_flag flag on the computation of the counting distributions.
 *
 *  \return                  SemiMarkovData object.
 */
/*--------------------------------------------------------------*/

SemiMarkovData* SemiMarkov::simulation(StatError &error , int nb_sequence ,
                                       const MarkovianSequences &iseq , bool counting_flag) const

{
  FrequencyDistribution *length_distribution;
  SemiMarkovData *seq;


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
 *  \brief Computation of Kullback-Leibler divergences between semi-Markov chains.
 *
 *  \param[in] error               reference on a StatError object,
 *  \param[in] display             flag for displaying the matrix of pairwise distances between models,
 *  \param[in] nb_model            number of semi-Markov chains,
 *  \param[in] ismarkov            pointer on SemiMarkov objects,
 *  \param[in] length_distribution sequence length frequency distribution,
 *  \param[in] path                file path.
 *
 *  \return                        DistanceMatrix object.
 */
/*--------------------------------------------------------------*/

DistanceMatrix* SemiMarkov::divergence_computation(StatError &error , bool display ,
                                                   int nb_model , const SemiMarkov **ismarkov ,
                                                   FrequencyDistribution **length_distribution ,
                                                   const string path) const

{
  bool status = true , lstatus;
  int i , j , k;
  int cumul_length , nb_failure;
  double **likelihood;
  long double divergence;
  const SemiMarkov **smarkov;
  MarkovianSequences *iseq , *seq;
  SemiMarkovData *simul_seq;
  DistanceMatrix *dist_matrix;
  ofstream *out_file;


  dist_matrix = NULL;
  error.init();

  for (i = 0;i < nb_model - 1;i++) {
    if (ismarkov[i]->type != type) {
      status = false;
      ostringstream error_message;
      error_message << SEQ_label[SEQL_SEMI_MARKOV_CHAIN] << " " << i + 2 << ": "
                    << SEQ_error[SEQR_MODEL_TYPE];
      error.update((error_message.str()).c_str());
    }

    if (ismarkov[i]->nb_output_process == nb_output_process) {
      if (ismarkov[i]->nb_state != nb_state) {
        status = false;
        ostringstream error_message;
        error_message << SEQ_label[SEQL_SEMI_MARKOV_CHAIN] << " " << i + 2 << ": "
                      << SEQ_error[SEQR_NB_STATE];
        error.update((error_message.str()).c_str());
      }

      if (nb_output_process == 1) {
        if (ismarkov[i]->categorical_process[0]->nb_value != categorical_process[0]->nb_value) {
          status = false;
          ostringstream error_message;
          error_message << SEQ_label[SEQL_SEMI_MARKOV_CHAIN] << " " << i + 2 << ": "
                        << STAT_error[STATR_NB_OUTPUT];
          error.update((error_message.str()).c_str());
        }
      }
    }

    else if ((nb_output_process == 0) && (ismarkov[i]->nb_output_process == 1)) {
      if (ismarkov[i]->categorical_process[0]->nb_value != nb_state) {
        status = false;
        ostringstream error_message;
        error_message << SEQ_label[SEQL_SEMI_MARKOV_CHAIN] << " " << i + 2 << ": "
                      << STAT_error[STATR_NB_OUTPUT];
        error.update((error_message.str()).c_str());
      }
    }

    else {  //  if ((nb_output_process == 1) && (ismarkov[i]->nb_output_process == 0))
      if (ismarkov[i]->nb_state != categorical_process[0]->nb_value) {
        status = false;
        ostringstream error_message;
        error_message << SEQ_label[SEQL_SEMI_MARKOV_CHAIN] << " " << i + 2 << ": "
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

    smarkov = new const SemiMarkov*[nb_model];

    smarkov[0] = this;
    for (i = 1;i < nb_model;i++) {
      smarkov[i] = ismarkov[i - 1];
    }

    dist_matrix = new DistanceMatrix(nb_model , SEQ_label[SEQL_SEMI_MARKOV_CHAIN]);

    for (i = 0;i < nb_model;i++) {

      // generation of a sample of sequences using a semi-Markov chain

      simul_seq = smarkov[i]->simulation(error , *length_distribution[i] , false , true);

      likelihood = new double*[simul_seq->nb_sequence];
      for (j = 0;j < simul_seq->nb_sequence;j++) {
        likelihood[j] = new double[nb_model];
      }

      for (j = 0;j < simul_seq->nb_sequence;j++) {
        likelihood[j][i] = smarkov[i]->likelihood_computation(*simul_seq , j);

        if ((display) && (likelihood[j][i] == D_INF)) {
          cout << "\nERROR - " << SEQ_error[SEQR_REFERENCE_MODEL] << ": " << i + 1 << endl;
        }
      }

      if (smarkov[i]->nb_output_process == 1) {
        iseq = simul_seq->remove_variable_1();
      }
      else {
        iseq = simul_seq;
      }

      // computation of the log-likelihood of each semi-Markov chain for the sample of sequences

      for (j = 0;j < nb_model;j++) {
        if (j != i) {
          if (smarkov[j]->nb_output_process == 1) {
            seq = iseq->transcode(error , smarkov[j]->categorical_process[0]);
          }
          else {
            seq = iseq;
          }

          divergence = 0.;
          cumul_length = 0;
          nb_failure = 0;

          for (k = 0;k < seq->nb_sequence;k++) {
            likelihood[k][j] = smarkov[j]->likelihood_computation(*seq , k);

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

          if (smarkov[j]->nb_output_process == 1) {
            delete seq;
          }
        }
      }

      if (display) {
        cout << SEQ_label[SEQL_SEMI_MARKOV_CHAIN] << " " << i + 1 << ": " << simul_seq->nb_sequence << " "
             << SEQ_label[SEQL_SIMULATED] << " " << SEQ_label[simul_seq->nb_sequence == 1 ? SEQL_SEQUENCE : SEQL_SEQUENCES] << endl;
        simul_seq->likelihood_write(cout , nb_model , likelihood , SEQ_label[SEQL_SEMI_MARKOV_CHAIN]);
      }
      if (out_file) {
        *out_file << SEQ_label[SEQL_SEMI_MARKOV_CHAIN] << " " << i + 1 << ": " << simul_seq->nb_sequence << " "
                  << SEQ_label[SEQL_SIMULATED] << " " << SEQ_label[simul_seq->nb_sequence == 1 ? SEQL_SEQUENCE : SEQL_SEQUENCES] << endl;
        simul_seq->likelihood_write(*out_file , nb_model , likelihood , SEQ_label[SEQL_SEMI_MARKOV_CHAIN]);
      }

      for (j = 0;j < simul_seq->nb_sequence;j++) {
        delete [] likelihood[j];
      }
      delete [] likelihood;

      if (smarkov[i]->nb_output_process == 1) {
        delete iseq;
      }
      delete simul_seq;
    }

    if (out_file) {
      out_file->close();
      delete out_file;
    }

    delete smarkov;
  }

  return dist_matrix;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of Kullback-Leibler divergences between semi-Markov chains.
 *
 *  \param[in] error       reference on a StatError object,
 *  \param[in] display     flag for displaying the matrix of pairwise distances between models,
 *  \param[in] nb_model    number of semi-Markov chains,
 *  \param[in] smarkov     pointer on SemiMarkov objects,
 *  \param[in] nb_sequence number of generated sequences,
 *  \param[in] length      sequence length,
 *  \param[in] path        file path.
 *
 *  \return                DistanceMatrix object.
 */
/*--------------------------------------------------------------*/

DistanceMatrix* SemiMarkov::divergence_computation(StatError &error , bool display ,
                                                   int nb_model , const SemiMarkov **smarkov ,
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

    dist_matrix = divergence_computation(error , display , nb_model , smarkov , length_distribution , path);

    for (i = 0;i < nb_model;i++) {
      delete length_distribution[i];
    }
    delete [] length_distribution;
  }

  return dist_matrix;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of Kullback-Leibler divergences between semi-Markov chains.
 *
 *  \param[in] error       reference on a StatError object,
 *  \param[in] display     flag for displaying the matrix of pairwise distances between models,
 *  \param[in] nb_model    number of semi-Markov chains,
 *  \param[in] smarkov     pointer on SemiMarkov objects,
 *  \param[in] nb_sequence number of generated sequences,
 *  \param[in] seq         pointer on MarkovianSequences objects,
 *  \param[in] path        file path.
 *
 *  \return                DistanceMatrix object.
 */
/*--------------------------------------------------------------*/

DistanceMatrix* SemiMarkov::divergence_computation(StatError &error , bool display ,
                                                   int nb_model , const SemiMarkov **smarkov ,
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

    dist_matrix = divergence_computation(error , display , nb_model , smarkov , length_distribution , path);

    for (i = 0;i < nb_model;i++) {
      delete length_distribution[i];
    }
    delete [] length_distribution;
  }

  return dist_matrix;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Constructor of the SemiMarkovIterator class.
 *
 *  \param[in] ismarkov pointer on a SemiMarkov object.
 */
/*--------------------------------------------------------------*/

SemiMarkovIterator::SemiMarkovIterator(SemiMarkov *ismarkov)

{
  semi_markov = ismarkov;
  (semi_markov->nb_iterator)++;

  if ((!(semi_markov->cumul_initial)) || (!(semi_markov->cumul_transition))) {
    semi_markov->create_cumul();
    semi_markov->cumul_computation();
  }

  state = I_DEFAULT;
  occupancy = 0;
  counter = I_DEFAULT;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Copy of a SemiMarkovIterator object.
 *
 *  \param[in] iter reference on a SemiMarkovIterator object.
 */
/*--------------------------------------------------------------*/

void SemiMarkovIterator::copy(const SemiMarkovIterator &iter)

{
  semi_markov = iter.semi_markov;
  (semi_markov->nb_iterator)++;

  state = iter.state;
  occupancy = iter.occupancy;
  counter = iter.counter;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Destructor of the SemiMarkovIterator class.
 */
/*--------------------------------------------------------------*/

SemiMarkovIterator::~SemiMarkovIterator()

{
  (semi_markov->nb_iterator)--;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Assignment operator of the SemiMarkovIterator class.
 *
 *  \param[in] iter reference on a SemiMarkovIterator object.
 *
 *  \return         SemiMarkovIterator object.
 */
/*--------------------------------------------------------------*/

SemiMarkovIterator& SemiMarkovIterator::operator=(const SemiMarkovIterator &iter)

{
  if (&iter != this) {
    (semi_markov->nb_iterator)--;
    copy(iter);
  }

  return *this;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Simulation using a semi-Markov chain.
 *
 *  \param[in] int_seq        sequence,
 *  \param[in] length         sequence length,
 *  \param[in] initialization flag initialization.
 *
 *  \return                   error status.
 */
/*--------------------------------------------------------------*/

bool SemiMarkovIterator::simulation(int **int_seq , int length , bool initialization)

{
  bool status;


  if ((state == I_DEFAULT) && (!initialization)) {
    status = false;
  }

  else {
    int i , j;
    int *pstate , **pioutput;
//    double **proutput;


    status = true;

    if (semi_markov->nb_output_process > 0) {
      pioutput = new int*[semi_markov->nb_output_process];
//      proutput = new double*[semi_markov->nb_output_process];
    }

    if (initialization) {
      state = cumul_method(semi_markov->nb_state , semi_markov->cumul_initial);

      switch (semi_markov->sojourn_type[state]) {

      case SEMI_MARKOVIAN : {
        switch (semi_markov->type) {
        case ORDINARY :
          occupancy = semi_markov->state_process->sojourn_time[state]->simulation();
          break;
        case EQUILIBRIUM :
          occupancy = semi_markov->forward[state]->simulation();
          break;
        }
        break;
      }

      case MARKOVIAN : {
        if (semi_markov->transition[state][state] < 1.) {
          occupancy = 1;
        }
        break;
      }
      }

      counter = 0;
    }

    pstate = int_seq[0];
    for (i = 0;i < semi_markov->nb_output_process;i++) {
/*      switch (type[i + 1]) {
      case INT_VALUE : */
      pioutput[i] = int_seq[i + 1];
/*        break;
      case REAL_VALUE :
        proutput[i] = real_seq[i + 1];
        break;
      } */
    }

    for (i = 0;i < length;i++) {
      counter++;
      *pstate++ = state;

      for (j = 0;j < semi_markov->nb_output_process;j++) {
        if (semi_markov->categorical_process[j]) {
          *pioutput[j]++ = semi_markov->categorical_process[j]->observation[state]->simulation();
        }
        else if (semi_markov->discrete_parametric_process[j]) {
          *pioutput[j]++ = semi_markov->discrete_parametric_process[j]->observation[state]->simulation();
        }
        else {
//          *proutput[j]++ = semi_markov->continuous_parametric_process[j]->observation[state]->simulation();
        }
      }

      if ((semi_markov->transition[state][state] < 1.) && (counter == occupancy)) {
        state = cumul_method(semi_markov->nb_state , semi_markov->cumul_transition[state]);

        switch (semi_markov->sojourn_type[state]) {
        case SEMI_MARKOVIAN :
          occupancy = semi_markov->state_process->sojourn_time[state]->simulation();
          break;
        case MARKOVIAN :
          occupancy = 1;
          break;
        }

        counter = 0;
      }
    }

    if (semi_markov->nb_output_process > 0) {
      delete [] pioutput;
//      delete [] proutput;
    }
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Simulation using a semi-Markov chain.
 *
 *  \param[in] length         sequence length,
 *  \param[in] initialization flag initialization.
 *
 *  \return                   generated sequence.
 */
/*--------------------------------------------------------------*/

int** SemiMarkovIterator::simulation(int length , bool initialization)

{
  int i;
  int **int_seq;


  if ((state == I_DEFAULT) && (!initialization)) {
    int_seq = NULL;
  }

  else {
    int_seq = new int*[semi_markov->nb_output_process + 1];
    for (i = 0;i <= semi_markov->nb_output_process;i++) {
      int_seq[i] = new int[length];
    }

    simulation(int_seq , length , initialization);
  }

  return int_seq;
}


};  // namespace sequence_analysis
