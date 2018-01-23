/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       StructureAnalysis: Exploring and Analyzing Plant Architecture
 *
 *       Copyright 1995-2018 CIRAD AGAP
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
#include <sstream>

#include "stat_tool/stat_label.h"

#include "stat_tool/distribution_reestimation.hpp"   // problem compiler C++ Windows

#include "hidden_semi_markov.h"
#include "sequence_label.h"

using namespace std;
using namespace boost::math;
using namespace stat_tool;


namespace sequence_analysis {



/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the log-likelihood of a hidden semi-Markov chain for sequences
 *         using the forward algorithm.
 *
 *  \param[in] seq                   reference on a MarkovianSequences object,
 *  \param[in] posterior_probability pointer on the posterior probabilities of the most probable state sequences,
 *  \param[in] index                 sequence index.
 *
 *  \return                          log-likelihood.
 */
/*--------------------------------------------------------------*/

double HiddenSemiMarkov::likelihood_computation(const MarkovianSequences &seq ,
                                                double *posterior_probability , int index) const

{
  int i , j , k , m;
  int nb_value , length , **pioutput;
  double likelihood = 0. , seq_likelihood , obs_product , residual , **observation ,
         *norm , *state_norm , *forward1 , **state_in , **proutput;
  DiscreteParametric *occupancy;


  // checking of the compatibility of the model with the data

  if (nb_output_process == seq.nb_variable) {
    for (i = 0;i < nb_output_process;i++) {
      if ((categorical_process[i]) || (discrete_parametric_process[i])) {
        if (categorical_process[i]) {
          nb_value = categorical_process[i]->nb_value;
        }
        else {
          nb_value = discrete_parametric_process[i]->nb_value;
        }

        if (nb_value < seq.marginal_distribution[i]->nb_value) {
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

    // initializations

    length = (index == I_DEFAULT ? seq.max_length : seq.length[index]);

    observation = new double*[length];
    for (i = 0;i < length;i++) {
      observation[i] = new double[nb_state];
    }

    norm = new double[length];
    state_norm = new double[nb_state];
    forward1 = new double[nb_state];

    state_in = new double*[length - 1];
    for (i = 0;i < length - 1;i++) {
      state_in[i] = new double[nb_state];
    }

    pioutput = new int*[seq.nb_variable];
    proutput = new double*[seq.nb_variable];

    for (i = 0;i < seq.nb_sequence;i++) {
      if ((index == I_DEFAULT) || (index == i)) {
        for (j = 0;j < seq.nb_variable;j++) {
          switch (seq.type[j]) {
          case INT_VALUE :
            pioutput[j] = seq.int_sequence[i][j];
            break;
          case REAL_VALUE :
            proutput[j] = seq.real_sequence[i][j];
            break;
          }
        }
        seq_likelihood = 0.;

        for (j = 0;j < seq.length[i];j++) {
          norm[j] = 0.;

          for (k = 0;k < nb_state;k++) {

            // computation of the observation probabilities

            observation[j][k] = 1.;
            for (m = 0;m < nb_output_process;m++) {
              if (categorical_process[m]) {
                observation[j][k] *= categorical_process[m]->observation[k]->mass[*pioutput[m]];
              }

              else if (discrete_parametric_process[m]) {
                observation[j][k] *= discrete_parametric_process[m]->observation[k]->mass[*pioutput[m]];
              }

              else {
                if (((continuous_parametric_process[m]->ident == GAMMA) ||
                     (continuous_parametric_process[m]->ident == ZERO_INFLATED_GAMMA)) && (seq.min_value[m] < seq.min_interval[m] / 2)) {
                  switch (seq.type[m]) {
                  case INT_VALUE :
                    observation[j][k] *= continuous_parametric_process[m]->observation[k]->mass_computation(*pioutput[m] , *pioutput[m] + seq.min_interval[m]);
                    break;
                  case REAL_VALUE :
                    observation[j][k] *= continuous_parametric_process[m]->observation[k]->mass_computation(*proutput[m] , *proutput[m] + seq.min_interval[m]);
                    break;
                  }
                }

                else if (continuous_parametric_process[m]->ident == LINEAR_MODEL) {
                  switch (seq.type[m]) {
                  case INT_VALUE :
                    residual = *pioutput[m] - (continuous_parametric_process[m]->observation[k]->intercept +
                                continuous_parametric_process[m]->observation[k]->slope *
                                (seq.index_param_type == IMPLICIT_TYPE ? j : seq.index_parameter[i][j]));
                    break;
                  case REAL_VALUE :
                    residual = *proutput[m] - (continuous_parametric_process[m]->observation[k]->intercept +
                                continuous_parametric_process[m]->observation[k]->slope *
                                (seq.index_param_type == IMPLICIT_TYPE ? j : seq.index_parameter[i][j]));
                    break;
                  }

                  observation[j][k] *= continuous_parametric_process[m]->observation[k]->mass_computation(residual , residual);
                }

                else if (continuous_parametric_process[m]->ident == AUTOREGRESSIVE_MODEL) {
                  if (j == 0) {
                    switch (seq.type[m]) {
                    case INT_VALUE :
                      residual = *pioutput[m] - continuous_parametric_process[m]->observation[k]->location;
                      break;
                    case REAL_VALUE :
                      residual = *proutput[m] - continuous_parametric_process[m]->observation[k]->location;
                      break;
                    }
                  }

                  else {
                    switch (seq.type[m]) {
                    case INT_VALUE :
                      residual = *pioutput[m] - (continuous_parametric_process[m]->observation[k]->location +
                                  continuous_parametric_process[m]->observation[k]->autoregressive_coeff *
                                  (*(pioutput[m] - 1) - continuous_parametric_process[m]->observation[k]->location));
                      break;
                    case REAL_VALUE :
                      residual = *proutput[m] - (continuous_parametric_process[m]->observation[k]->location +
                                  continuous_parametric_process[m]->observation[k]->autoregressive_coeff *
                                  (*(proutput[m] - 1) - continuous_parametric_process[m]->observation[k]->location));
                      break;
                    }
                  }

                  observation[j][k] *= continuous_parametric_process[m]->observation[k]->mass_computation(residual , residual);
                }

                else {
                  switch (seq.type[m]) {
                  case INT_VALUE :
                    observation[j][k] *= continuous_parametric_process[m]->observation[k]->mass_computation(*pioutput[m] - seq.min_interval[m] / 2 , *pioutput[m] + seq.min_interval[m] / 2);
                    break;
                  case REAL_VALUE :
                    observation[j][k] *= continuous_parametric_process[m]->observation[k]->mass_computation(*proutput[m] - seq.min_interval[m] / 2 , *proutput[m] + seq.min_interval[m] / 2);
                    break;
                  }
                }
              }
            }

            switch (sojourn_type[k]) {

            // case semi-Markovian state

            case SEMI_MARKOVIAN : {
              if (j == 0) {
                state_norm[k] = initial[k];
              }
              else {
                state_norm[k] += state_in[j - 1][k] - forward1[k];
              }
              state_norm[k] *= observation[j][k];

              norm[j] += state_norm[k];
              break;
            }

            // case Markovian state

            case MARKOVIAN : {
              if (j == 0) {
                forward1[k] = initial[k];
              }
              else {
                forward1[k] = state_in[j - 1][k];
              }
              forward1[k] *= observation[j][k];

              norm[j] += forward1[k];
              break;
            }
            }
          }

          if (norm[j] > 0.) {
            for (k = 0;k < nb_state;k++) {
              switch (sojourn_type[k]) {
              case SEMI_MARKOVIAN :
                state_norm[k] /= norm[j];
                break;
              case MARKOVIAN :
                forward1[k] /= norm[j];
                break;
              }
            }

            seq_likelihood += log(norm[j]);
          }

          else {
            seq_likelihood = D_INF;
            break;
          }

          for (k = 0;k < nb_state;k++) {

            // case semi-Markovian state

            if (sojourn_type[k] == SEMI_MARKOVIAN) {
              occupancy = state_process->sojourn_time[k];
              obs_product = 1.;
              forward1[k] = 0.;

              if (j < seq.length[i] - 1) {
                for (m = 1;m <= MIN(j + 1 , occupancy->nb_value - 1);m++) {
                  obs_product *= observation[j - m + 1][k] / norm[j - m + 1];
                  if (obs_product == 0.) {
                    break;
                  }

                  if (m < j + 1) {
                    forward1[k] += obs_product * occupancy->mass[m] * state_in[j - m][k];
                  }

                  else {
                    switch (type) {
                    case ORDINARY :
                      forward1[k] += obs_product * occupancy->mass[m] * initial[k];
                      break;
                    case EQUILIBRIUM :
                      forward1[k] += obs_product * forward[k]->mass[m] * initial[k];
                      break;
                    }
                  }
                }
              }

              else {
                for (m = 1;m <= MIN(j + 1 , occupancy->nb_value - 1);m++) {
                  obs_product *= observation[j - m + 1][k] / norm[j - m + 1];
                  if (obs_product == 0.) {
                    break;
                  }

                  if (m < j + 1) {
                    forward1[k] += obs_product * (1. - occupancy->cumul[m - 1]) * state_in[j - m][k];
                  }

                  else {
                    switch (type) {
                    case ORDINARY :
                      forward1[k] += obs_product * (1. - occupancy->cumul[m - 1]) * initial[k];
                      break;
                    case EQUILIBRIUM :
                      forward1[k] += obs_product * (1. - forward[k]->cumul[m - 1]) * initial[k];
                      break;
                    }
                  }
                }
              }
            }
          }

          if (j < seq.length[i] - 1) {
            for (k = 0;k < nb_state;k++) {
              state_in[j][k] = 0.;
              for (m = 0;m < nb_state;m++) {
                state_in[j][k] += transition[m][k] * forward1[m];
              }
            }
          }

          for (k = 0;k < seq.nb_variable;k++) {
            switch (seq.type[k]) {
            case INT_VALUE :
              pioutput[k]++;
              break;
            case REAL_VALUE :
              proutput[k]++;
              break;
            }
          }
        }

        if (seq_likelihood != D_INF) {
          likelihood += seq_likelihood;
          if (posterior_probability) {
            posterior_probability[i] = exp(posterior_probability[i] - seq_likelihood);
          }
        }

        else {
          likelihood = D_INF;
          break;
        }
      }
    }

    for (i = 0;i < length;i++) {
      delete [] observation[i];
    }
    delete [] observation;

    delete [] norm;
    delete [] state_norm;
    delete [] forward1;

    for (i = 0;i < length - 1;i++) {
      delete [] state_in[i];
    }
    delete [] state_in;

    delete [] pioutput;
    delete [] proutput;
  }

  return likelihood;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Estimation of a hidden semi-Markov chain using the EM algorithm.
 *
 *  \param[in] error             reference on a StatError object,
 *  \param[in] display           flag for displaying estimation intermediate results,
 *  \param[in] ihsmarkov         initial hidden semi-Markov chain,
 *  \param[in] poisson_geometric flag on the estimation of Poisson geometric state occupancy distributions,
 *  \param[in] common_dispersion flag common dispersion parameter (continuous observation processes),
 *  \param[in] estimator         estimator type for the reestimation of the state occupancy distributions
 *                               (complete or partial likelihood),
 *  \param[in] counting_flag     flag on the computation of the counting distributions,
 *  \param[in] state_sequence    flag on the computation of the restored state sequences,
 *  \param[in] nb_iter           number of iterations,
 *  \param[in] mean_estimator    method for the computation of the state occupancy
 *                               distribution mean (equilibrium semi-Markov chain).
 *
 *  \return                      HiddenSemiMarkov object.
 */
/*--------------------------------------------------------------*/

HiddenSemiMarkov* MarkovianSequences::hidden_semi_markov_estimation(StatError &error , bool display ,
                                                                    const HiddenSemiMarkov &ihsmarkov ,
                                                                    bool poisson_geometric , bool common_dispersion ,
                                                                    censoring_estimator estimator , bool counting_flag ,
                                                                    bool state_sequence , int nb_iter ,
                                                                    duration_distribution_mean_estimator mean_estimator) const

{
  bool status;
  int i , j , k , m , n;
  int max_nb_value , iter , nb_likelihood_decrease , offset , nb_value , *occupancy_nb_value ,
      *censored_occupancy_nb_value , **pioutput;
  double likelihood = D_INF , previous_likelihood , occupancy_likelihood , observation_likelihood ,
         min_likelihood , obs_product , residual , buff , sum , occupancy_mean , **observation ,
         *norm , *state_norm , **forward1 , **state_in , *backward , **backward1 , *auxiliary ,
         *ofrequency , *lfrequency , *occupancy_survivor , *censored_occupancy_survivor , diff ,
         variance , **mean_direction , global_mean_direction , concentration , ***state_sequence_count ,
         **proutput;
  double *complete_occupancy_weight , *censored_occupancy_weight;
  Distribution *weight;
  DiscreteParametric *occupancy;
  ChainReestimation<double> *chain_reestim;
  Reestimation<double> **occupancy_reestim , **length_bias_reestim , **censored_occupancy_reestim ,
                       ***observation_reestim;
  FrequencyDistribution *hoccupancy , *hobservation;
  HiddenSemiMarkov *hsmarkov;
  SemiMarkovData *seq;

# ifdef DEBUG
  double test[NB_STATE][4];
# endif


  hsmarkov = NULL;
  error.init();

  // test number of values for each variable

  status = false;
  for (i = 0;i < nb_variable;i++) {
    if (max_value[i] > min_value[i]) {
      status = true;
      break;
    }
  }

  if (!status) {
    error.update(STAT_error[STATR_VARIABLE_NB_VALUE]);
  }

  for (i = 0;i < nb_variable;i++) {
    if ((type[i] != INT_VALUE) && (type[i] != REAL_VALUE) && (type[i] != STATE)) {
      status = false;
      ostringstream error_message , correction_message;
      error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                    << STAT_error[STATR_VARIABLE_TYPE];
      correction_message << STAT_variable_word[INT_VALUE] << " or "
                         << STAT_variable_word[REAL_VALUE];
      error.correction_update((error_message.str()).c_str() , (correction_message.str()).c_str());
    }
  }

  if (ihsmarkov.nb_output_process != nb_variable) {
    status = false;
    error.update(STAT_error[STATR_NB_OUTPUT_PROCESS]);
  }

  else {
    for (i = 0;i < nb_variable;i++) {
      if ((ihsmarkov.categorical_process[i]) || (ihsmarkov.discrete_parametric_process[i])) {
        if (type[i] == REAL_VALUE) {
          status = false;
          ostringstream error_message;
          error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                        << STAT_error[STATR_VARIABLE_TYPE];
          error.correction_update((error_message.str()).c_str() , STAT_variable_word[INT_VALUE]);
        }

        else {
          if (min_value[i] < 0) {
            status = false;
            ostringstream error_message;
            error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                          << STAT_error[STATR_POSITIVE_MIN_VALUE];
            error.update((error_message.str()).c_str());
          }

          if (!marginal_distribution[i]) {
            status = false;
            ostringstream error_message;
            error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                          << STAT_error[STATR_MARGINAL_FREQUENCY_DISTRIBUTION];
            error.update((error_message.str()).c_str());
          }

          else {
            if (((ihsmarkov.categorical_process[i]) &&
                 (ihsmarkov.categorical_process[i]->nb_value != marginal_distribution[i]->nb_value)) ||
                ((ihsmarkov.discrete_parametric_process[i]) &&
                 (ihsmarkov.discrete_parametric_process[i]->nb_value < marginal_distribution[i]->nb_value))) {
              status = false;
              ostringstream error_message;
              error_message << STAT_label[STATL_OUTPUT_PROCESS] << " " << i + 1 << ": "
                            << STAT_error[STATR_NB_OUTPUT];
              error.update((error_message.str()).c_str());
            }

            else if ((ihsmarkov.categorical_process[i]) && (!characteristics[i])) {
              for (j = 0;j < marginal_distribution[i]->nb_value;j++) {
                if (marginal_distribution[i]->frequency[j] == 0) {
                  status = false;
                  ostringstream error_message;
                  error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                                << STAT_error[STATR_MISSING_VALUE] << " " << j;
                  error.update((error_message.str()).c_str());
                }
              }
            }
          }
        }
      }

      else if ((ihsmarkov.continuous_parametric_process[i]) &&
               (ihsmarkov.continuous_parametric_process[i]->ident == LINEAR_MODEL) &&
               (ihsmarkov.nb_component < ihsmarkov.nb_state)) {
        status = false;
        error.update(SEQ_error[SEQR_MODEL_STRUCTURE]);
      }
    }
  }

  if ((nb_iter != I_DEFAULT) && (nb_iter < 1)) {
    status = false;
    error.update(STAT_error[STATR_NB_ITERATION]);
  }

  if (status) {
    if (max_length > COUNTING_MAX_LENGTH) {
      counting_flag = false;
    }

    // construction of the hidden semi-Markov chain

    hsmarkov = new HiddenSemiMarkov(ihsmarkov , false , (int)(max_length * SAMPLE_NB_VALUE_COEFF));

    if (hsmarkov->type == EQUILIBRIUM) {
      for (i = 0;i < hsmarkov->nb_state;i++) {
        hsmarkov->initial[i] = 1. / (double)hsmarkov->nb_state;
      }
    }

    if (common_dispersion) {
      for (i = 0;i < hsmarkov->nb_output_process;i++) {
        if (hsmarkov->continuous_parametric_process[i]) {
          hsmarkov->continuous_parametric_process[i]->tied_dispersion = true;
        }
      }
    }

#   ifdef DEBUG
    cout << *hsmarkov;
#   endif

    // construction of the data structures of the algorithm

    observation = new double*[max_length];
    for (i = 0;i < max_length;i++) {
      observation[i] = new double[hsmarkov->nb_state];
    }

    norm = new double[max_length];
    state_norm = new double[hsmarkov->nb_state];

    forward1 = new double*[max_length];
    for (i = 0;i < max_length;i++) {
      forward1[i] = new double[hsmarkov->nb_state];
    }

    state_in = new double*[max_length - 1];
    for (i = 0;i < max_length - 1;i++) {
      state_in[i] = new double[hsmarkov->nb_state];
    }

    backward = new double[hsmarkov->nb_state];

    backward1 = new double*[max_length];
    for (i = 0;i < max_length;i++) {
      backward1[i] = new double[hsmarkov->nb_state];
    }

    auxiliary = new double[hsmarkov->nb_state];

    chain_reestim = new ChainReestimation<double>(hsmarkov->type , hsmarkov->nb_state , hsmarkov->nb_state);

    occupancy_nb_value = new int[hsmarkov->nb_state];
    occupancy_reestim = new Reestimation<double>*[hsmarkov->nb_state];
    if (hsmarkov->type == EQUILIBRIUM) {
      length_bias_reestim = new Reestimation<double>*[hsmarkov->nb_state];
    }

    for (i = 0;i < hsmarkov->nb_state;i++) {
      switch (hsmarkov->sojourn_type[i]) {

      case SEMI_MARKOVIAN : {
        if (estimator == COMPLETE_LIKELIHOOD) {
          occupancy_nb_value[i] = hsmarkov->state_process->sojourn_time[i]->alloc_nb_value;
        }
        else {
          occupancy_nb_value[i] = MIN(hsmarkov->state_process->sojourn_time[i]->alloc_nb_value ,
                                      max_length);
        }

        occupancy_reestim[i] = new Reestimation<double>(occupancy_nb_value[i]);
        if (hsmarkov->type == EQUILIBRIUM) {
          length_bias_reestim[i] = new Reestimation<double>(occupancy_nb_value[i]);
        }
        break;
      }

      case MARKOVIAN : {
        occupancy_reestim[i] = NULL;
        if (hsmarkov->type == EQUILIBRIUM) {
          length_bias_reestim[i] = NULL;
        }
        break;
      }
      }
    }

    max_nb_value = 0;
    for (i = 0;i < hsmarkov->nb_state;i++) {
      if ((hsmarkov->sojourn_type[i] == SEMI_MARKOVIAN) && (occupancy_nb_value[i] > max_nb_value)) {
        max_nb_value = occupancy_nb_value[i];
      }
    }

    if (estimator == KAPLAN_MEIER) {
      censored_occupancy_nb_value = new int[hsmarkov->nb_state];
      censored_occupancy_reestim = new Reestimation<double>*[hsmarkov->nb_state];
      for (i = 0;i < hsmarkov->nb_state;i++) {
        switch (hsmarkov->sojourn_type[i]) {
        case SEMI_MARKOVIAN :
          censored_occupancy_nb_value[i] = MIN(hsmarkov->state_process->sojourn_time[i]->alloc_nb_value ,
                                               max_length + 1);
          censored_occupancy_reestim[i] = new Reestimation<double>(censored_occupancy_nb_value[i]);
          break;
        case MARKOVIAN :
          censored_occupancy_reestim[i] = NULL;
          break;
        }
      }

      occupancy_survivor = new double[max_nb_value];
      censored_occupancy_survivor = new double[max_nb_value + 1];
    }

    hoccupancy = new FrequencyDistribution(max_nb_value);

    if ((display) && (hsmarkov->type == ORDINARY)) {
      complete_occupancy_weight = new double[hsmarkov->nb_state];
      censored_occupancy_weight = new double[hsmarkov->nb_state];
    }

    observation_reestim = new Reestimation<double>**[hsmarkov->nb_output_process];
    for (i = 0;i < hsmarkov->nb_output_process;i++) {
      if ((marginal_distribution[i]) && ((!(hsmarkov->continuous_parametric_process[i])) ||
           ((hsmarkov->continuous_parametric_process[i]->ident != LINEAR_MODEL) &&
            (hsmarkov->continuous_parametric_process[i]->ident != AUTOREGRESSIVE_MODEL)))) {
        observation_reestim[i] = new Reestimation<double>*[hsmarkov->nb_state];
        for (j = 0;j < hsmarkov->nb_state;j++) {
          observation_reestim[i][j] = new Reestimation<double>(marginal_distribution[i]->nb_value);
        }
      }

      else {
        observation_reestim[i] = NULL;
      }
    }

    max_nb_value = 0;
    for (i = 0;i < hsmarkov->nb_output_process;i++) {
      if ((hsmarkov->discrete_parametric_process[i]) &&
          (max_nb_value < marginal_distribution[i]->nb_value)) {
        max_nb_value = marginal_distribution[i]->nb_value;
      }
    }

    if (max_nb_value > 0) {
      hobservation = new FrequencyDistribution(max_nb_value);
    }
    else {
      hobservation = NULL;
    }

    for (i = 0;i < hsmarkov->nb_output_process;i++) {
      if ((hsmarkov->continuous_parametric_process[i]) &&
          (hsmarkov->continuous_parametric_process[i]->ident == VON_MISES)) {
        break;
      }
    }

    if (i < hsmarkov->nb_output_process) {
      mean_direction = new double*[hsmarkov->nb_state];
      for (i = 0;i < hsmarkov->nb_state;i++) {
        mean_direction[i] = new double[4];
      }
    }
    else {
      mean_direction = NULL;
    }

    for (i = 0;i < hsmarkov->nb_output_process;i++) {
      if ((!marginal_distribution[i]) || ((hsmarkov->continuous_parametric_process[i]) &&
           ((hsmarkov->continuous_parametric_process[i]->ident == LINEAR_MODEL) ||
            (hsmarkov->continuous_parametric_process[i]->ident == AUTOREGRESSIVE_MODEL)))) {
        break;
      }
    }

    if (i < hsmarkov->nb_output_process) {
      state_sequence_count = new double**[nb_sequence];
      for (i = 0;i < nb_sequence;i++) {
        state_sequence_count[i] = new double*[length[i]];
        for (j = 0;j < length[i];j++) {
          state_sequence_count[i][j] = new double[hsmarkov->nb_state];
        }
      }
    }
    else {
      state_sequence_count = NULL;
    }

    pioutput = new int*[nb_variable];
    proutput = new double*[nb_variable];

    iter = 0;
    nb_likelihood_decrease = 0;

    do {
      iter++;
      previous_likelihood = likelihood;
      likelihood = 0.;

      // initialization of the reestimation quantities

      chain_reestim->init();

      for (i = 0;i < hsmarkov->nb_state;i++) {
        if (hsmarkov->sojourn_type[i] == SEMI_MARKOVIAN) {
          for (j = 0;j < occupancy_nb_value[i];j++) {
            occupancy_reestim[i]->frequency[j] = 0.;
          }

          if (hsmarkov->type == EQUILIBRIUM) {
            for (j = 0;j < occupancy_nb_value[i];j++) {
              length_bias_reestim[i]->frequency[j] = 0.;
            }
          }

          if (estimator == KAPLAN_MEIER) {
            for (j = 0;j < censored_occupancy_nb_value[i];j++) {
              censored_occupancy_reestim[i]->frequency[j] = 0.;
            }
          }

          if ((display) && (hsmarkov->type == ORDINARY)) {
            complete_occupancy_weight[i] = 0.;
            censored_occupancy_weight[i] = 0.;
          }
        }
      }

      for (i = 0;i < hsmarkov->nb_output_process;i++) {
        if (observation_reestim[i]) {
          for (j = 0;j < hsmarkov->nb_state;j++) {
            for (k = 0;k < marginal_distribution[i]->nb_value;k++) {
              observation_reestim[i][j]->frequency[k] = 0.;
            }
          }
        }
      }

      if (state_sequence_count) {
        for (i = 0;i < nb_sequence;i++) {
          for (j = 0;j < length[i];j++) {
            for (k = 0;k < hsmarkov->nb_state;k++) {
              state_sequence_count[i][j][k] = 0.;
            }
          }
        }
      }

#     ifdef DEBUG
      for (i = 0;i < hsmarkov->nb_state;i++) {
        for (j = 0;j < 4;j++) {
          test[i][j] = 0.;
        }
      }
#     endif

      for (i = 0;i < nb_sequence;i++) {
        for (j = 0;j < nb_variable;j++) {
          switch (type[j]) {
          case INT_VALUE :
            pioutput[j] = int_sequence[i][j];
            break;
          case REAL_VALUE :
            proutput[j] = real_sequence[i][j];
            break;
          }
        }

        // forward recurrence

        for (j = 0;j < length[i];j++) {
          norm[j] = 0.;

          for (k = 0;k < hsmarkov->nb_state;k++) {

            // computation of the observation probabilities

            observation[j][k] = 1.;
            for (m = 0;m < hsmarkov->nb_output_process;m++) {
              if (hsmarkov->categorical_process[m]) {
                observation[j][k] *= hsmarkov->categorical_process[m]->observation[k]->mass[*pioutput[m]];
              }

              else if (hsmarkov->discrete_parametric_process[m]) {
                observation[j][k] *= hsmarkov->discrete_parametric_process[m]->observation[k]->mass[*pioutput[m]];
              }

              else {
                if (((hsmarkov->continuous_parametric_process[m]->ident == GAMMA) ||
                     (hsmarkov->continuous_parametric_process[m]->ident == ZERO_INFLATED_GAMMA)) && (min_value[m] < min_interval[m] / 2)) {
                  switch (type[m]) {
                  case INT_VALUE :
                    observation[j][k] *= hsmarkov->continuous_parametric_process[m]->observation[k]->mass_computation(*pioutput[m] , *pioutput[m] + min_interval[m]);
                    break;
                  case REAL_VALUE :
                    observation[j][k] *= hsmarkov->continuous_parametric_process[m]->observation[k]->mass_computation(*proutput[m] , *proutput[m] + min_interval[m]);
                    break;
                  }
                }

                else if (hsmarkov->continuous_parametric_process[m]->ident == LINEAR_MODEL) {
                  switch (type[m]) {
                  case INT_VALUE :
                    residual = *pioutput[m] - (hsmarkov->continuous_parametric_process[m]->observation[k]->intercept +
                                hsmarkov->continuous_parametric_process[m]->observation[k]->slope *
                                (index_param_type == IMPLICIT_TYPE ? j : index_parameter[i][j]));
                    break;
                  case REAL_VALUE :
                    residual = *proutput[m] - (hsmarkov->continuous_parametric_process[m]->observation[k]->intercept +
                                hsmarkov->continuous_parametric_process[m]->observation[k]->slope *
                                (index_param_type == IMPLICIT_TYPE ? j : index_parameter[i][j]));
                    break;
                  }

                  observation[j][k] *= hsmarkov->continuous_parametric_process[m]->observation[k]->mass_computation(residual , residual);

#                 ifdef DEBUG
                  cout << STAT_label[STATL_STATE] << " " << k << "  " << SEQ_label[SEQL_SEQUENCE] << " " << i << "  "
                       << SEQ_label[SEQL_INDEX] << " " << j << ": " << residual << " "
                       << hsmarkov->continuous_parametric_process[m]->observation[k]->mass_computation(residual , residual) << endl;
#                 endif

                }

                else if (hsmarkov->continuous_parametric_process[m]->ident == AUTOREGRESSIVE_MODEL) {
                  if (j == 0) {
                    switch (type[m]) {
                    case INT_VALUE :
                      residual = *pioutput[m] - hsmarkov->continuous_parametric_process[m]->observation[k]->location;
                      break;
                    case REAL_VALUE :
                      residual = *proutput[m] - hsmarkov->continuous_parametric_process[m]->observation[k]->location;
                      break;
                    }
                  }

                  else {
                    switch (type[m]) {
                    case INT_VALUE :
                      residual = *pioutput[m] - (hsmarkov->continuous_parametric_process[m]->observation[k]->location +
                                  hsmarkov->continuous_parametric_process[m]->observation[k]->autoregressive_coeff *
                                  (*(pioutput[m] - 1) - hsmarkov->continuous_parametric_process[m]->observation[k]->location));
                      break;
                    case REAL_VALUE :
                      residual = *proutput[m] - (hsmarkov->continuous_parametric_process[m]->observation[k]->location +
                                  hsmarkov->continuous_parametric_process[m]->observation[k]->autoregressive_coeff *
                                  (*(proutput[m] - 1) - hsmarkov->continuous_parametric_process[m]->observation[k]->location));
                      break;
                    }
                  }

                  observation[j][k] *= hsmarkov->continuous_parametric_process[m]->observation[k]->mass_computation(residual , residual);
                }

                else {
                  switch (type[m]) {
                  case INT_VALUE :
                    observation[j][k] *= hsmarkov->continuous_parametric_process[m]->observation[k]->mass_computation(*pioutput[m] - min_interval[m] / 2 , *pioutput[m] + min_interval[m] / 2);
                    break;
                  case REAL_VALUE :
                    observation[j][k] *= hsmarkov->continuous_parametric_process[m]->observation[k]->mass_computation(*proutput[m] - min_interval[m] / 2 , *proutput[m] + min_interval[m] / 2);
                    break;
                  }
                }
              }
            }

            switch (hsmarkov->sojourn_type[k]) {

            // case semi-Markovian state

            case SEMI_MARKOVIAN : {
              if (j == 0) {
                state_norm[k] = hsmarkov->initial[k];
              }
              else {
                state_norm[k] += state_in[j - 1][k] - forward1[j - 1][k];
              }
              state_norm[k] *= observation[j][k];

              norm[j] += state_norm[k];
              break;
            }

            // case Markovian state

            case MARKOVIAN : {
              if (j == 0) {
                forward1[j][k] = hsmarkov->initial[k];
              }
              else {
                forward1[j][k] = state_in[j - 1][k];
              }
              forward1[j][k] *= observation[j][k];

              norm[j] += forward1[j][k];
              break;
            }
            }
          }

          if (norm[j] > 0.) {
            for (k = 0;k < hsmarkov->nb_state;k++) {
              switch (hsmarkov->sojourn_type[k]) {
              case SEMI_MARKOVIAN :
                state_norm[k] /= norm[j];
                break;
              case MARKOVIAN :
                forward1[j][k] /= norm[j];
                break;
              }
            }

            likelihood += log(norm[j]);
          }

          else {
            likelihood = D_INF;
            break;
          }

          for (k = 0;k < hsmarkov->nb_state;k++) {

            // case semi-Markovian state

            if (hsmarkov->sojourn_type[k] == SEMI_MARKOVIAN) {
              occupancy = hsmarkov->state_process->sojourn_time[k];
              obs_product = 1.;
              forward1[j][k] = 0.;

              if (j < length[i] - 1) {
                for (m = 1;m <= MIN(j + 1 , occupancy->nb_value - 1);m++) {
                  obs_product *= observation[j - m + 1][k] / norm[j - m + 1];
                  if (obs_product == 0.) {
                    break;
                  }

                  if (m < j + 1) {
                    forward1[j][k] += obs_product * occupancy->mass[m] * state_in[j - m][k];
                  }

                  else {
                    switch (hsmarkov->type) {
                    case ORDINARY :
                      forward1[j][k] += obs_product * occupancy->mass[m] * hsmarkov->initial[k];
                      break;
                    case EQUILIBRIUM :
                      forward1[j][k] += obs_product * hsmarkov->forward[k]->mass[m] * hsmarkov->initial[k];
                      break;
                    }
                  }
                }
              }

              else {
                for (m = 1;m <= MIN(j + 1 , occupancy->nb_value - 1);m++) {
                  obs_product *= observation[j - m + 1][k] / norm[j - m + 1];
                  if (obs_product == 0.) {
                    break;
                  }

                  if (m < j + 1) {
                    forward1[j][k] += obs_product * (1. - occupancy->cumul[m - 1]) * state_in[j - m][k];
                  }

                  else {
                    switch (hsmarkov->type) {
                    case ORDINARY :
                      forward1[j][k] += obs_product * (1. - occupancy->cumul[m - 1]) *
                                        hsmarkov->initial[k];
                      break;
                    case EQUILIBRIUM :
                      forward1[j][k] += obs_product * (1. - hsmarkov->forward[k]->cumul[m - 1]) *
                                        hsmarkov->initial[k];
                      break;
                    }
                  }
                }
              }
            }
          }

          if (j < length[i] - 1) {
            for (k = 0;k < hsmarkov->nb_state;k++) {
              state_in[j][k] = 0.;
              for (m = 0;m < hsmarkov->nb_state;m++) {
                state_in[j][k] += hsmarkov->transition[m][k] * forward1[j][m];
              }
            }
          }

          for (k = 0;k < nb_variable;k++) {
            switch (type[k]) {
            case INT_VALUE :
              pioutput[k]++;
              break;
            case REAL_VALUE :
              proutput[k]++;
              break;
            }
          }
        }

        if (likelihood == D_INF) {
          break;
        }

#       ifdef DEBUG
        for (j = 0;j < length[i];j++) {
          cout << j << " : ";
          for (k = 0;k < hsmarkov->nb_state;k++) {
            cout << forward1[j][k] << " ";
//            cout << observation[j][k] << " ";
          }
          cout << endl;
        }
        cout << endl;
#       endif

        // backward recurrence

        for (j = 0;j < nb_variable;j++) {
          if (type[j] == INT_VALUE) {
            pioutput[j]--;
          }
        }

        j = length[i] - 1;
        for (k = 0;k < hsmarkov->nb_state;k++) {
          backward[k] = forward1[j][k];
          backward1[j][k] = backward[k];

          // accumulation of the reestimation quantities of the observation distributions

          for (m = 0;m < hsmarkov->nb_output_process;m++) {
            if (observation_reestim[m]) {
              observation_reestim[m][k]->frequency[*pioutput[m]] += backward[k];
            }
          }

          if (state_sequence_count) {
            state_sequence_count[i][j][k] += backward[k];
          }
        }

        for (j = length[i] - 2;j >= 0;j--) {
          for (k = 0;k < nb_variable;k++) {
            if (type[k] == INT_VALUE) {
              pioutput[k]--;
            }
          }

          for (k = 0;k < hsmarkov->nb_state;k++) {
            auxiliary[k] = 0.;

            switch (hsmarkov->sojourn_type[k]) {

            // case semi-Markovian state

            case SEMI_MARKOVIAN : {
              occupancy = hsmarkov->state_process->sojourn_time[k];
              obs_product = 1.;

              for (m = 1;m < MIN(length[i] - j , occupancy->nb_value);m++) {
                obs_product *= observation[j + m][k] / norm[j + m];
                if (obs_product == 0.) {
                  break;
                }

                if (backward1[j + m][k] > 0.) {
//                if (forward1[j + m][k] > 0.) {
                  if (m < length[i] - j - 1) {
                    buff = backward1[j + m][k] * obs_product * occupancy->mass[m] /
                           forward1[j + m][k];

                    // accumulation of the reestimation quantities of the state occupancy distributions

                    occupancy_reestim[k]->frequency[m] += buff * state_in[j][k];
                  }

                  else {
                    buff = obs_product * (1. - occupancy->cumul[m - 1]);

                    // accumulation of the reestimation quantities of the state occupancy distributions

                    switch (estimator) {

                    case COMPLETE_LIKELIHOOD : {
                      for (n = m;n < occupancy->nb_value;n++) {
                        occupancy_reestim[k]->frequency[n] += obs_product * occupancy->mass[n] *
                                                              state_in[j][k];
                      }
                      break;
                    }

                    case KAPLAN_MEIER : {
                      censored_occupancy_reestim[k]->frequency[m] += buff * state_in[j][k];
                      break;
                    }
                    }
                  }

                  auxiliary[k] += buff;
                }
              }
              break;
            }

            // case Markovian state

            case MARKOVIAN : {
              if (backward1[j + 1][k] > 0.) {
//              if (forward1[j + 1][k] > 0.) {
                auxiliary[k] = backward1[j + 1][k] / state_in[j][k];

/*                auxiliary[k] = backward1[j + 1][k] * observation[j + 1][k] /
                               (forward1[j + 1][k] * norm[j + 1]); */

              }
              break;
            }
            }
          }

          for (k = 0;k < hsmarkov->nb_state;k++) {
            backward1[j][k] = 0.;

            for (m = 0;m < hsmarkov->nb_state;m++) {
              buff = auxiliary[m] * hsmarkov->transition[k][m] * forward1[j][k];
              backward1[j][k] += buff;

              // accumulation of the reestimation quantities of the transition probabilities

              chain_reestim->transition[k][m] += buff;
            }

            switch (hsmarkov->sojourn_type[k]) {

            // case semi-Markovian state

            case SEMI_MARKOVIAN : {
              backward[k] = backward[k] + backward1[j][k] - auxiliary[k] * state_in[j][k];
              if (backward[k] < 0.) {
                backward[k] = 0.;
              }
              if (backward[k] > 1.) {
                backward[k] = 1.;
              }
              break;
            }

            // case Markovian state

            case MARKOVIAN : {
              backward[k] = backward1[j][k];
              break;
            }
            }

            // accumulation of the reestimation quantities of the observation distributions

            for (m = 0;m < hsmarkov->nb_output_process;m++) {
              if (observation_reestim[m]) {
                observation_reestim[m][k]->frequency[*pioutput[m]] += backward[k];
              }
            }
          }

          if (state_sequence_count) {
            for (k = 0;k < hsmarkov->nb_state;k++) {
              state_sequence_count[i][j][k] += backward[k];
            }
          }
        }

        // accumulation of the reestimation quantities of the initial probabilities

        if (hsmarkov->type == ORDINARY) {
          for (j = 0;j < hsmarkov->nb_state;j++) {
            chain_reestim->initial[j] += backward[j];
          }
        }

        // accumulation of the reestimation quantities of the initial state occupancy distributions

        if ((hsmarkov->type == ORDINARY) || (estimator == COMPLETE_LIKELIHOOD)) {
          for (j = 0;j < hsmarkov->nb_state;j++) {
            if ((hsmarkov->sojourn_type[j] == SEMI_MARKOVIAN) && (hsmarkov->initial[j] > 0.)) {
              occupancy = hsmarkov->state_process->sojourn_time[j];
              obs_product = 1.;
              if (hsmarkov->type == EQUILIBRIUM) {
                sum = 0.;
              }

              for (k = 1;k < MIN(length[i] + 1 , occupancy->nb_value);k++) {
                obs_product *= observation[k - 1][j] / norm[k - 1];
                if (obs_product == 0.) {
                  break;
                }

                if (backward1[k - 1][j] > 0.) {
//                if (forward1[k - 1][j] > 0.) {
                  if (k < length[i]) {
                    switch (hsmarkov->type) {
                    case ORDINARY :
                      occupancy_reestim[j]->frequency[k] += backward1[k - 1][j] * obs_product *
                                                            occupancy->mass[k] * hsmarkov->initial[j] /
                                                            forward1[k - 1][j];
                      break;
                    case EQUILIBRIUM :
                      sum += backward1[k - 1][j] * obs_product / forward1[k - 1][j];
                      length_bias_reestim[j]->frequency[k] += sum * occupancy->mass[k] * hsmarkov->initial[j] /
                                                              occupancy->mean;
                      break;
                    }
                  }

                  else {
                    switch (estimator) {

                    case COMPLETE_LIKELIHOOD : {
                      for (m = k;m < occupancy->nb_value;m++) {
                        switch (hsmarkov->type) {
                        case ORDINARY :
                          occupancy_reestim[j]->frequency[m] += obs_product * occupancy->mass[m] *
                                                                hsmarkov->initial[j];
                          break;
                        case EQUILIBRIUM :
                          length_bias_reestim[j]->frequency[m] += (sum + obs_product * (m + 1 - k)) * occupancy->mass[m] *
                                                                  hsmarkov->initial[j] / occupancy->mean;
                          break;
                        }
                      }
                      break;
                    }

                    case KAPLAN_MEIER : {
                      censored_occupancy_reestim[j]->frequency[k] += obs_product *
                                                                     (1. - occupancy->cumul[k - 1]) *
                                                                     hsmarkov->initial[j];
                      break;
                    }
                    }
                  }
                }
              }
            }
          }
        }

#       ifdef DEBUG
        for (j = length[i] - 1;j >= 0;j--) {
          cout << j << " : ";
          double sum = 0.;
          for (k = 0;k < hsmarkov->nb_state;k++) {
            sum += backward[j][k];
            cout << backward[j][k];
            if ((hsmarkov->sojourn_type[k] == SEMI_MARKOVIAN) && (j < length[i] - 1)){
              cout << " (" << backward1[j][k] << ") ";
            }
          }
          cout << "| " << sum << endl;

          for (k = 0;k < hsmarkov->nb_state;k++) {
            if (hsmarkov->sojourn_type[k] == SEMI_MARKOVIAN) {
              if (j < length[i] - 1) {
                test[k][0] += backward1[j][k];
                test[k][1] += auxiliary[k] * state_in[j][k];
              }
              else {
                test[k][2] += backward[j][k];
              }
              if (j == 0) {
                test[k][3] += backward[j][k];
              }
            }
          }
        }
#       endif

        if ((display) && (hsmarkov->type == ORDINARY)) {
          for (j = 0;j < hsmarkov->nb_state;j++) {
            if (hsmarkov->sojourn_type[j] == SEMI_MARKOVIAN) {
              for (k = 0;k < length[i] - 1;k++) {
                complete_occupancy_weight[j] += backward1[k][j];
              }
              censored_occupancy_weight[j] += backward1[length[i] - 1][j];
            }
          }
        }
      }

      if (likelihood != D_INF) {
        if (likelihood < previous_likelihood) {
          nb_likelihood_decrease++;
        }
        else {
          nb_likelihood_decrease = 0;
        }

        // reestimation of the initial probabilities

        if (hsmarkov->type == ORDINARY) {
          reestimation(hsmarkov->nb_state , chain_reestim->initial ,
                       hsmarkov->initial , MIN_PROBABILITY , false);
        }

        // reestimation of the transition probabilities

        for (i = 0;i < hsmarkov->nb_state;i++) {
          reestimation(hsmarkov->nb_state , chain_reestim->transition[i] ,
                       hsmarkov->transition[i] , MIN_PROBABILITY , false);
        }

        // reestimation of the state occupancy distributions

        min_likelihood = 0.;

        for (i = 0;i < hsmarkov->nb_state;i++) {
          if (hsmarkov->sojourn_type[i] == SEMI_MARKOVIAN) {
            occupancy = hsmarkov->state_process->sojourn_time[i];

            if (estimator == KAPLAN_MEIER) {
              occupancy_reestim[i]->nb_value_computation();
              occupancy_reestim[i]->offset_computation();
              occupancy_reestim[i]->nb_element_computation();

              censored_occupancy_reestim[i]->nb_value_computation();
              censored_occupancy_reestim[i]->offset_computation();
              censored_occupancy_reestim[i]->nb_element_computation();

              if (censored_occupancy_reestim[i]->nb_element > 0.) {

#               ifdef DEBUG
                cout << "\n" << STAT_label[STATL_STATE] << " " << i << " (" << test[i][2]
                     << " | " << censored_occupancy_reestim[i]->nb_element << ") - ";

                occupancy_reestim[i]->max_computation();
                occupancy_reestim[i]->mean_computation();
                occupancy_reestim[i]->variance_computation();

                occupancy_reestim[i]->ascii_characteristic_print(cout);
#               endif

                occupancy_reestim[i]->state_occupancy_estimation(censored_occupancy_reestim[i] ,
                                                                 occupancy_reestim[i] ,
                                                                 occupancy_survivor ,
                                                                 censored_occupancy_survivor , false);
              }
            }

#           ifdef DEBUG
            cout << STAT_label[STATL_STATE] << " " << i << " (";
#           endif

            if ((hsmarkov->type == ORDINARY) || (estimator == PARTIAL_LIKELIHOOD)) {
              occupancy_reestim[i]->nb_value_computation();
              occupancy_reestim[i]->offset_computation();
              occupancy_reestim[i]->nb_element_computation();
              occupancy_reestim[i]->max_computation();
              occupancy_reestim[i]->mean_computation();
              occupancy_reestim[i]->variance_computation();

#             ifdef DEBUG
              if (hsmarkov->type == ORDINARY) {
                switch (estimator) {
                case COMPLETE_LIKELIHOOD :
                  cout << test[i][0] + test[i][2] << " | " << test[i][1] + test[i][3];
                  break;
                case PARTIAL_LIKELIHOOD :
                  cout << test[i][0];
                  break;
                }
                cout << " | " << occupancy_reestim[i]->nb_element << ") - ";
                occupancy_reestim[i]->ascii_characteristic_print(cout);
              }
#             endif

            }

            else {
              offset = 1;
              nb_value = occupancy_nb_value[i];

              ofrequency = occupancy_reestim[i]->frequency + occupancy_nb_value[i];
              lfrequency = length_bias_reestim[i]->frequency + occupancy_nb_value[i];
              while ((*--ofrequency == 0) && (*--lfrequency == 0) && (nb_value > 2)) {
                nb_value--;
              }
              occupancy_reestim[i]->nb_value = nb_value;
              length_bias_reestim[i]->nb_value = nb_value;

              ofrequency = occupancy_reestim[i]->frequency + offset;
              lfrequency = length_bias_reestim[i]->frequency + offset;
              while ((*ofrequency++ == 0) && (*lfrequency++ == 0) && (offset < nb_value - 1)) {
                offset++;
              }
              occupancy_reestim[i]->offset = offset;
              length_bias_reestim[i]->offset = offset;

              occupancy_reestim[i]->nb_element_computation();
              length_bias_reestim[i]->nb_element_computation();

#             ifdef DEBUG
              occupancy_reestim[i]->max_computation();
              occupancy_reestim[i]->mean_computation();
              occupancy_reestim[i]->variance_computation();

              cout << test[i][1] << " | " << occupancy_reestim[i]->nb_element << ") - ";
              occupancy_reestim[i]->ascii_characteristic_print(cout);

              length_bias_reestim[i]->max_computation();
              length_bias_reestim[i]->mean_computation();
              length_bias_reestim[i]->variance_computation();

              cout << STAT_label[STATL_STATE] << " " << i << " (" << test[i][3] << " | "
                   << length_bias_reestim[i]->nb_element << ") - ";
              length_bias_reestim[i]->ascii_characteristic_print(cout);
#             endif

              switch (mean_estimator) {
              case COMPUTED :
                occupancy_mean = interval_bisection(occupancy_reestim[i] , length_bias_reestim[i]);
                break;
              case ONE_STEP_LATE :
                occupancy_mean = occupancy->mean;
                break;
              }

#             ifdef DEBUG
              cout << STAT_label[STATL_SOJOURN_TIME] << " " << STAT_label[STATL_MEAN] << ": "
                   << occupancy_mean << endl;
#             endif

              occupancy_reestim[i]->equilibrium_process_combination(length_bias_reestim[i] , occupancy_mean);

#             ifdef DEBUG
              cout << test[i][0] + test[i][2] << " | " << test[i][1] + test[i][3] << " | "
                   << occupancy_reestim[i]->nb_element << ") - ";
              occupancy_reestim[i]->ascii_characteristic_print(cout);
#             endif
            }

            hoccupancy->update(occupancy_reestim[i] ,
                               MAX((int)(occupancy_reestim[i]->nb_element *
                                         MAX(sqrt(occupancy_reestim[i]->variance) , 1.) * OCCUPANCY_COEFF) , MIN_NB_ELEMENT));
            if (iter <= EXPLORATION_NB_ITER) {
              occupancy_likelihood = hoccupancy->Reestimation<int>::parametric_estimation(occupancy , 1 , true ,
                                                                                          OCCUPANCY_THRESHOLD , poisson_geometric);
            }
            else {
              occupancy_likelihood = hoccupancy->Reestimation<int>::type_parametric_estimation(occupancy , 1 , true ,
                                                                                               OCCUPANCY_THRESHOLD , poisson_geometric);
            }

            if (occupancy_likelihood == D_INF) {
              min_likelihood = D_INF;
            }
            else {
              occupancy->computation(hoccupancy->nb_value , OCCUPANCY_THRESHOLD);
              if (hsmarkov->type == EQUILIBRIUM) {
                hsmarkov->forward[i]->copy(*occupancy);
                hsmarkov->forward[i]->computation(*occupancy);
              }
            }

#           ifdef DEBUG
            cout << STAT_word[STATW_STATE] << " " << i << " " << SEQ_word[SEQW_OCCUPANCY_DISTRIBUTION] << endl;
            occupancy->ascii_print(cout);
#           endif

          }
        }

        if (hsmarkov->type == EQUILIBRIUM) {
          hsmarkov->initial_probability_computation();
        }

        // reestimation of the observation distributions

        for (i = 0;i < hsmarkov->nb_output_process;i++) {
          if (hsmarkov->categorical_process[i]) {
            for (j = 0;j < hsmarkov->nb_state;j++) {
              reestimation(marginal_distribution[i]->nb_value , observation_reestim[i][j]->frequency ,
                           hsmarkov->categorical_process[i]->observation[j]->mass ,
                           MIN_PROBABILITY , false);
            }
          }

          else if (observation_reestim[i]) {
            for (j = 0;j < hsmarkov->nb_state;j++) {
              observation_reestim[i][j]->nb_value_computation();
              observation_reestim[i][j]->offset_computation();
              observation_reestim[i][j]->nb_element_computation();
              observation_reestim[i][j]->max_computation();
              if ((hsmarkov->discrete_parametric_process[i]) ||
                  (hsmarkov->continuous_parametric_process[i]->ident != ZERO_INFLATED_GAMMA)) {
                observation_reestim[i][j]->mean_computation();
                observation_reestim[i][j]->variance_computation(true);
//                observation_reestim[i][j]->variance_computation();
              }
            }

            if (hsmarkov->discrete_parametric_process[i]) {
              for (j = 0;j < hsmarkov->nb_state;j++) {
                hobservation->update(observation_reestim[i][j] ,
                                     MAX((int)(observation_reestim[i][j]->nb_element *
                                               MAX(sqrt(observation_reestim[i][j]->variance) , 1.) * OBSERVATION_COEFF) , MIN_NB_ELEMENT));
                observation_likelihood = hobservation->Reestimation<int>::type_parametric_estimation(hsmarkov->discrete_parametric_process[i]->observation[j] ,
                                                                                                     0 , true , OBSERVATION_THRESHOLD);

                if (observation_likelihood == D_INF) {
                  min_likelihood = D_INF;
                }
                else {
                  hsmarkov->discrete_parametric_process[i]->observation[j]->computation(marginal_distribution[i]->nb_value ,
                                                                                        OBSERVATION_THRESHOLD);

                  if (hsmarkov->discrete_parametric_process[i]->observation[j]->ident == BINOMIAL) {
                    for (k = hsmarkov->discrete_parametric_process[i]->observation[j]->nb_value;k < marginal_distribution[i]->nb_value;k++) {
                      hsmarkov->discrete_parametric_process[i]->observation[j]->mass[k] = 0.;
                    }
                  }
                }
              }
            }

            else {
              switch (hsmarkov->continuous_parametric_process[i]->ident) {

              case GAMMA : {
                for (j = 0;j < hsmarkov->nb_state;j++) {
                  observation_reestim[i][j]->gamma_estimation(hsmarkov->continuous_parametric_process[i]->observation[j] , iter);
                }
                break;
              }

              case ZERO_INFLATED_GAMMA : {
                for (j = 0;j < hsmarkov->nb_state;j++) {
                  observation_reestim[i][j]->zero_inflated_gamma_estimation(hsmarkov->continuous_parametric_process[i]->observation[j] , iter);
                }
                break;
              }

              case GAUSSIAN : {
                for (j = 0;j < hsmarkov->nb_state;j++) {
                  hsmarkov->continuous_parametric_process[i]->observation[j]->location = observation_reestim[i][j]->mean;
                }

                if (common_dispersion) {
                  variance = 0.;
                  buff = 0.;

                  for (j = 0;j < hsmarkov->nb_state;j++) {
                    for (k = observation_reestim[i][j]->offset;k < observation_reestim[i][j]->nb_value;k++) {
                      diff = k - observation_reestim[i][j]->mean;
                      variance += observation_reestim[i][j]->frequency[k] * diff * diff;
                    }

                    buff += observation_reestim[i][j]->nb_element;
                  }

                  variance /= buff;
//                  variance /= (buff - 1);

                  for (j = 0;j < hsmarkov->nb_state;j++) {
                    hsmarkov->continuous_parametric_process[i]->observation[j]->dispersion = sqrt(variance);
                  }
                }

                else {
                  for (j = 0;j < hsmarkov->nb_state;j++) {
                    hsmarkov->continuous_parametric_process[i]->observation[j]->dispersion = sqrt(observation_reestim[i][j]->variance);
                    if (hsmarkov->continuous_parametric_process[i]->observation[j]->dispersion /
                        hsmarkov->continuous_parametric_process[i]->observation[j]->location < GAUSSIAN_MIN_VARIATION_COEFF) {
                      hsmarkov->continuous_parametric_process[i]->observation[j]->dispersion = hsmarkov->continuous_parametric_process[i]->observation[j]->location * GAUSSIAN_MIN_VARIATION_COEFF;
                    }
                  }
                }

                break;
              }

              case VON_MISES : {
                for (j = 0;j < hsmarkov->nb_state;j++) {
                  observation_reestim[i][j]->mean_direction_computation(mean_direction[j]);
                  hsmarkov->continuous_parametric_process[i]->observation[j]->location = mean_direction[j][3];
                }

                if (common_dispersion) {
                  global_mean_direction = 0.;
                  buff = 0.;

                  for (j = 0;j < hsmarkov->nb_state;j++) {
                    global_mean_direction += observation_reestim[i][j]->nb_element * mean_direction[j][2];
                    buff += observation_reestim[i][j]->nb_element;
                  }
                  concentration = von_mises_concentration_computation(global_mean_direction / buff);

                  for (j = 0;j < hsmarkov->nb_state;j++) {
                    hsmarkov->continuous_parametric_process[i]->observation[j]->dispersion = concentration;
                  }
                }

                else {
                  for (j = 0;j < hsmarkov->nb_state;j++) {
                    hsmarkov->continuous_parametric_process[i]->observation[j]->dispersion = von_mises_concentration_computation(mean_direction[j][2]);
                  }
                }
                break;
              }
              }
            }
          }

          else {
            switch (hsmarkov->continuous_parametric_process[i]->ident) {
            case GAMMA :
              gamma_estimation(state_sequence_count , i ,
                               hsmarkov->continuous_parametric_process[i] , iter);
              break;
            case ZERO_INFLATED_GAMMA :
              zero_inflated_gamma_estimation(state_sequence_count , i ,
                                             hsmarkov->continuous_parametric_process[i] , iter);
              break;
            case GAUSSIAN :
              gaussian_estimation(state_sequence_count , i ,
                                  hsmarkov->continuous_parametric_process[i]);
              break;
            case VON_MISES :
              von_mises_estimation(state_sequence_count , i ,
                                   hsmarkov->continuous_parametric_process[i]);
              break;
            case LINEAR_MODEL :
              linear_model_estimation(state_sequence_count , i ,
                                      hsmarkov->continuous_parametric_process[i]);
              break;
            case AUTOREGRESSIVE_MODEL :
              autoregressive_model_estimation(state_sequence_count , i ,
                                              hsmarkov->continuous_parametric_process[i]);
              break;
            }
          }
        }
      }

      if (display) {
        cout << STAT_label[STATL_ITERATION] << " " << iter << "   "
             << SEQ_label[SEQL_OBSERVED_SEQUENCES_LIKELIHOOD] << ": " << likelihood << endl;
      }

#     ifdef DEBUG
      if (iter % 5 == 0) {
        cout << *hsmarkov;
      }
#     endif

    }
    while ((likelihood != D_INF) && (((nb_iter == I_DEFAULT) && (iter < SEMI_MARKOV_NB_ITER) &&
             (((likelihood - previous_likelihood) / -likelihood > SEMI_MARKOV_LIKELIHOOD_DIFF) ||
              (min_likelihood == D_INF) || (nb_likelihood_decrease == 1))) ||
            ((nb_iter != I_DEFAULT) && (iter < nb_iter))));

    if (likelihood != D_INF) {
      if (display) {
        cout << "\n" << iter << " " << STAT_label[STATL_ITERATIONS] << endl;

        if (hsmarkov->type == ORDINARY) {
          cout << "\n" << SEQ_label[SEQL_OCCUPANCY_WEIGHTS] << endl;
          for (i = 0;i < hsmarkov->nb_state;i++) {
            if (hsmarkov->sojourn_type[i] == SEMI_MARKOVIAN) {
              cout << STAT_label[STATL_STATE] << " " << i << ": " << complete_occupancy_weight[i] << ", "
                   << censored_occupancy_weight[i];
              if ((complete_occupancy_weight[i] > 0.) && (censored_occupancy_weight[i] > 0.)) {
                cout <<"  (" << complete_occupancy_weight[i] / (complete_occupancy_weight[i] + censored_occupancy_weight[i]) << ", "
                     << censored_occupancy_weight[i] / (complete_occupancy_weight[i] + censored_occupancy_weight[i]) << ")";
              }
              cout << endl;
            }
          }
        }
      }

      // reestimation of the initial probabilities

      if (hsmarkov->type == ORDINARY) {
        reestimation(hsmarkov->nb_state , chain_reestim->initial ,
                     hsmarkov->initial , MIN_PROBABILITY , true);
      }

      // reestimation of the transition probabilities

      for (i = 0;i < hsmarkov->nb_state;i++) {
        reestimation(hsmarkov->nb_state , chain_reestim->transition[i] ,
                     hsmarkov->transition[i] , MIN_PROBABILITY , true);
      }

      if (hsmarkov->type == EQUILIBRIUM) {
        hsmarkov->initial_probability_computation();
      }

      for (i = 0;i < hsmarkov->nb_state;i++) {
        if ((hsmarkov->sojourn_type[i] == SEMI_MARKOVIAN) &&
            (hsmarkov->state_process->sojourn_time[i]->mean == 1.)) {
          hsmarkov->sojourn_type[i] = MARKOVIAN;
          delete hsmarkov->state_process->sojourn_time[i];
          hsmarkov->state_process->sojourn_time[i] = NULL;
          delete hsmarkov->forward[i];
          hsmarkov->forward[i] = NULL;
        }
      }

      // reestimation of the categorical observation distributions

      for (i = 0;i < hsmarkov->nb_output_process;i++) {
        if (hsmarkov->categorical_process[i]) {
          for (j = 0;j < hsmarkov->nb_state;j++) {
            reestimation(marginal_distribution[i]->nb_value , observation_reestim[i][j]->frequency ,
                         hsmarkov->categorical_process[i]->observation[j]->mass ,
                         MIN_PROBABILITY , true);
          }
        }

        else if (hsmarkov->discrete_parametric_process[i]) {
          hsmarkov->discrete_parametric_process[i]->nb_value_computation();
        }
      }
    }

    // destruction of the data structures of the algorithm

    for (i = 0;i < max_length;i++) {
      delete [] observation[i];
    }
    delete [] observation;

    delete [] norm;
    delete [] state_norm;

    for (i = 0;i < max_length;i++) {
      delete [] forward1[i];
    }
    delete [] forward1;

    for (i = 0;i < max_length - 1;i++) {
      delete [] state_in[i];
    }
    delete [] state_in;

    delete [] backward;

    for (i = 0;i < max_length;i++) {
      delete [] backward1[i];
    }
    delete [] backward1;

    delete [] auxiliary;

    delete chain_reestim;

    delete [] occupancy_nb_value;

    for (i = 0;i < hsmarkov->nb_state;i++) {
      delete occupancy_reestim[i];
    }
    delete [] occupancy_reestim;

    if (hsmarkov->type == EQUILIBRIUM) {
      for (i = 0;i < hsmarkov->nb_state;i++) {
        delete length_bias_reestim[i];
      }
      delete [] length_bias_reestim;
    }

    if (estimator == KAPLAN_MEIER) {
      delete [] censored_occupancy_nb_value;

      for (i = 0;i < hsmarkov->nb_state;i++) {
        delete censored_occupancy_reestim[i];
      }
      delete [] censored_occupancy_reestim;

      delete [] occupancy_survivor;
      delete [] censored_occupancy_survivor;
    }

    delete hoccupancy;

    if ((display) && (hsmarkov->type == ORDINARY)) {
      delete [] complete_occupancy_weight;
      delete [] censored_occupancy_weight;
    }

    for (i = 0;i < hsmarkov->nb_output_process;i++) {
      if (observation_reestim[i]) {
        for (j = 0;j < hsmarkov->nb_state;j++) {
          delete observation_reestim[i][j];
        }
        delete [] observation_reestim[i];
      }
    }
    delete [] observation_reestim;

    delete hobservation;

    if (mean_direction) {
      for (i = 0;i < hsmarkov->nb_state;i++) {
        delete [] mean_direction[i];
      }
      delete [] mean_direction;
    }

    if (state_sequence_count) {
      for (i = 0;i < nb_sequence;i++) {
        for (j = 0;j < length[i];j++) {
          delete [] state_sequence_count[i][j];
        }
        delete [] state_sequence_count[i];
      }
      delete [] state_sequence_count;
    }

    delete [] pioutput;
    delete [] proutput;

    if (likelihood == D_INF) {
      delete hsmarkov;
      hsmarkov = NULL;
      error.update(STAT_error[STATR_ESTIMATION_FAILURE]);
    }

    else {
      if (state_sequence) {
        hsmarkov->semi_markov_data = new SemiMarkovData(*this , ADD_STATE_VARIABLE ,
                                                        (hsmarkov->type == EQUILIBRIUM ? true : false));
        seq = hsmarkov->semi_markov_data;

        for (i = 0;i < hsmarkov->nb_output_process;i++) {
          if (((hsmarkov->discrete_parametric_process[i]) || (hsmarkov->continuous_parametric_process[i])) &&
              (seq->characteristics[i + 1])) {
            delete seq->characteristics[i + 1];
            seq->characteristics[i + 1] = NULL;
          }
        }

        hsmarkov->forward_backward(*seq);

        hsmarkov->create_cumul();
        hsmarkov->log_computation();
        hsmarkov->viterbi(*seq);
        hsmarkov->remove_cumul();

        seq->min_value_computation(0);
        seq->max_value_computation(0);
        seq->build_marginal_frequency_distribution(0);
        seq->build_characteristic(0 , true , (hsmarkov->type == EQUILIBRIUM ? true : false));

        seq->build_transition_count(hsmarkov);
        seq->build_observation_frequency_distribution(hsmarkov->nb_state);
        seq->build_observation_histogram(hsmarkov->nb_state);

        // computation of the state occupancy distributions

        for (i = 0;i < hsmarkov->nb_state;i++) {
          if (hsmarkov->sojourn_type[i] == SEMI_MARKOVIAN) {
            hsmarkov->state_process->sojourn_time[i]->computation((seq->characteristics[0] ? seq->characteristics[0]->sojourn_time[i]->nb_value : 1) ,
                                                                  OCCUPANCY_THRESHOLD);
            if (hsmarkov->stype[i] == RECURRENT) {
              if (hsmarkov->type == ORDINARY) {
                hsmarkov->forward[i]->copy(*(hsmarkov->state_process->sojourn_time[i]));
              }
              hsmarkov->forward[i]->computation(*(hsmarkov->state_process->sojourn_time[i]));
            }
          }
        }

        // computation of the mixtures of observation distributions (weights deduced from the restoration)

        weight = NULL;

        for (i = 0;i < hsmarkov->nb_output_process;i++) {
          if ((hsmarkov->categorical_process[i]) || (hsmarkov->discrete_parametric_process[i]) ||
              ((hsmarkov->continuous_parametric_process[i]) &&
               (hsmarkov->continuous_parametric_process[i]->ident != LINEAR_MODEL) &&
               (hsmarkov->continuous_parametric_process[i]->ident != AUTOREGRESSIVE_MODEL))) {
            weight = seq->weight_computation();
            break;
          }
        }

        for (i = 0;i < hsmarkov->nb_output_process;i++) {
          if (hsmarkov->categorical_process[i]) {
            hsmarkov->categorical_process[i]->restoration_weight = new Distribution(*weight);
            hsmarkov->categorical_process[i]->restoration_mixture = hsmarkov->categorical_process[i]->mixture_computation(hsmarkov->categorical_process[i]->restoration_weight);
          }

          else if (hsmarkov->discrete_parametric_process[i]) {
            for (j = 0;j < hsmarkov->nb_state;j++) {
              hsmarkov->discrete_parametric_process[i]->observation[j]->cumul_computation();
//              hsmarkov->discrete_parametric_process[i]->observation[j]->computation(seq->observation_distribution[i + 1][j]->nb_value ,
//                                                                                    OBSERVATION_THRESHOLD);
            }

            hsmarkov->discrete_parametric_process[i]->restoration_weight = new Distribution(*weight);
            hsmarkov->discrete_parametric_process[i]->restoration_mixture = hsmarkov->discrete_parametric_process[i]->mixture_computation(hsmarkov->discrete_parametric_process[i]->restoration_weight);
          }

          else if ((hsmarkov->continuous_parametric_process[i]) &&
                   (hsmarkov->continuous_parametric_process[i]->ident != LINEAR_MODEL) &&
                   (hsmarkov->continuous_parametric_process[i]->ident != AUTOREGRESSIVE_MODEL)) {
            hsmarkov->continuous_parametric_process[i]->restoration_weight = new Distribution(*weight);
          }
        }

        delete weight;

        if ((display) && (seq->characteristics[0])) {
          cout << "\n" << SEQ_label[SEQL_STATE_SEQUENCES_LIKELIHOOD] << ": " << seq->restoration_likelihood;

          for (i = 0;i < nb_variable;i++) {
            if (type[i] == REAL_VALUE) {
              break;
            }
          }
          if (i == nb_variable) {
            cout << " | " << hsmarkov->SemiMarkov::likelihood_computation(*seq);
          }
          cout << endl;
        }
      }

      else {
        if (hsmarkov->type == ORDINARY) {
          for (i = 0;i < hsmarkov->nb_state;i++) {
            if ((hsmarkov->sojourn_type[i] == SEMI_MARKOVIAN) && (hsmarkov->stype[i] == RECURRENT)) {
              hsmarkov->forward[i]->copy(*(hsmarkov->state_process->sojourn_time[i]));
              hsmarkov->forward[i]->computation(*(hsmarkov->state_process->sojourn_time[i]));
            }
          }
        }

        hsmarkov->semi_markov_data = new SemiMarkovData(*this , SEQUENCE_COPY ,
                                                        (hsmarkov->type == EQUILIBRIUM ? true : false));
        seq = hsmarkov->semi_markov_data;
        if (seq->type[0] == STATE) {
          seq->state_variable_init(INT_VALUE);
        }

        for (i = 0;i < hsmarkov->nb_output_process;i++) {
          if (((hsmarkov->discrete_parametric_process[i]) || (hsmarkov->continuous_parametric_process[i])) &&
              (seq->characteristics[i])) {
            delete seq->characteristics[i];
            seq->characteristics[i] = NULL;
          }
        }
      }

      for (i = 0;i < hsmarkov->nb_output_process;i++) {
        if (hsmarkov->categorical_process[i]) {
          for (j = 0;j < hsmarkov->nb_state;j++) {
            hsmarkov->categorical_process[i]->observation[j]->cumul_computation();

            hsmarkov->categorical_process[i]->observation[j]->max_computation();
//            hsmarkov->categorical_process[i]->observation[j]->mean_computation();
//            hsmarkov->categorical_process[i]->observation[j]->variance_computation();
          }
        }
      }

      // computation of the log-likelihood and the characteristic distributions of the model

      seq->likelihood = hsmarkov->likelihood_computation(*this , seq->posterior_probability);

      hsmarkov->component_computation();
      hsmarkov->characteristic_computation(*seq , counting_flag , I_DEFAULT , false);

      // computation of the mixtures of observation distributions (theoretical weights)

      weight = NULL;

      for (i = 0;i < hsmarkov->nb_output_process;i++) {
        if ((hsmarkov->categorical_process[i]) || (hsmarkov->discrete_parametric_process[i]) ||
            ((hsmarkov->continuous_parametric_process[i]) &&
             (hsmarkov->continuous_parametric_process[i]->ident != LINEAR_MODEL) &&
             (hsmarkov->continuous_parametric_process[i]->ident != AUTOREGRESSIVE_MODEL))) {
          switch (hsmarkov->type) {
          case ORDINARY :
            weight = hsmarkov->state_process->weight_computation();
            break;
          case EQUILIBRIUM :
            weight = new Distribution(hsmarkov->nb_state , hsmarkov->initial);
            break;
          }
          break;
        }
      }

      for (i = 0;i < hsmarkov->nb_output_process;i++) {
        if (hsmarkov->categorical_process[i]) {
          hsmarkov->categorical_process[i]->weight = new Distribution(*weight);
          hsmarkov->categorical_process[i]->mixture = hsmarkov->categorical_process[i]->mixture_computation(hsmarkov->categorical_process[i]->weight);
        }

        else if (hsmarkov->discrete_parametric_process[i]) {
          hsmarkov->discrete_parametric_process[i]->weight = new Distribution(*weight);
          hsmarkov->discrete_parametric_process[i]->mixture = hsmarkov->discrete_parametric_process[i]->mixture_computation(hsmarkov->discrete_parametric_process[i]->weight);
        }

        else if ((hsmarkov->continuous_parametric_process[i]) &&
                 (hsmarkov->continuous_parametric_process[i]->ident != LINEAR_MODEL) &&
                 (hsmarkov->continuous_parametric_process[i]->ident != AUTOREGRESSIVE_MODEL)) {
          hsmarkov->continuous_parametric_process[i]->weight = new Distribution(*weight);
        }
      }

      delete weight;

      if ((display) && (state_sequence) && (seq->nb_sequence <= POSTERIOR_PROBABILITY_NB_SEQUENCE)) {
        int *pstate;

        cout << "\n" << SEQ_label[SEQL_POSTERIOR_STATE_SEQUENCE_PROBABILITY] << endl;
        for (i = 0;i < seq->nb_sequence;i++) {
          cout << SEQ_label[SEQL_SEQUENCE] << " " << seq->identifier[i] << ": "
               << seq->posterior_probability[i];

          if (hsmarkov->nb_component == hsmarkov->nb_state) {
            cout << " | " << SEQ_label[SEQL_STATE_BEGIN] << ": ";

            pstate = seq->int_sequence[i][0] + 1;
            if (seq->index_parameter) {
              for (j = 1;j < seq->length[i];j++) {
                if (*pstate != *(pstate - 1)) {
                  cout << seq->index_parameter[i][j] << ", ";
                }
                pstate++;
              }
            }

            else {
              for (j = 1;j < seq->length[i];j++) {
                if (*pstate != *(pstate - 1)) {
                  cout << j << ", ";
                }
                pstate++;
              }
            }
          }

          cout << endl;
        }
      }
    }
  }

  return hsmarkov;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Estimation of a hidden semi-Markov chain using the EM algorithm.
 *
 *  \param[in] error             reference on a StatError object,
 *  \param[in] display           flag for displaying estimation intermediate results,
 *  \param[in] itype             process type (ORDINARY/EQUILIBRIUM),
 *  \param[in] nb_state          number of states,
 *  \param[in] left_right        flag on the Markov chain structure,
 *  \param[in] occupancy_mean    mean state occupancy,
 *  \param[in] poisson_geometric flag on the estimation of Poisson geometric state occupancy distributions,
 *  \param[in] common_dispersion flag common dispersion parameter (continuous observation processes),
 *  \param[in] estimator         estimator type for the reestimation of the state occupancy distributions
 *                               (complete or partial likelihood),
 *  \param[in] counting_flag     flag on the computation of the counting distributions,
 *  \param[in] state_sequence    flag on the computation of the restored state sequences,
 *  \param[in] nb_iter           number of iterations,
 *  \param[in] mean_estimator    method for the computation of the state occupancy
 *                               distribution mean (equilibrium semi-Markov chain).
 *
 *  \return                      HiddenSemiMarkov object.
 */
/*--------------------------------------------------------------*/

HiddenSemiMarkov* MarkovianSequences::hidden_semi_markov_estimation(StatError &error , bool display ,
                                                                    process_type itype , int nb_state ,
                                                                    bool left_right , double occupancy_mean ,
                                                                    bool poisson_geometric , bool common_dispersion ,
                                                                    censoring_estimator estimator , bool counting_flag ,
                                                                    bool state_sequence , int nb_iter ,
                                                                    duration_distribution_mean_estimator mean_estimator) const

{
  bool status = true;
  int i;
  int nb_value[SEQUENCE_NB_VARIABLE];
  double proba , mean , variance;
  HiddenSemiMarkov *ihsmarkov , *hsmarkov;


  hsmarkov = NULL;
  error.init();

  if ((nb_state < 2) || (nb_state > NB_STATE)) {
    status = false;
    error.update(SEQ_error[SEQR_NB_STATE]);
  }
  if ((occupancy_mean != D_DEFAULT) && (occupancy_mean <= 1.)) {
    status = false;
    error.update(SEQ_error[SEQR_OCCUPANCY]);
  }

  for (i = 0;i < nb_variable;i++) {
    if ((type[i] != INT_VALUE) && (type[i] != REAL_VALUE) && (type[i] != STATE)) {
      status = false;
      ostringstream error_message , correction_message;
      error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                    << STAT_error[STATR_VARIABLE_TYPE];
      correction_message << STAT_variable_word[INT_VALUE] << " or "
                         << STAT_variable_word[REAL_VALUE];
      error.correction_update((error_message.str()).c_str() , (correction_message.str()).c_str());
    }
  }

  if (status) {
    for (i = 0;i < nb_variable;i++) {
      if (marginal_distribution[i]) {
        nb_value[i] = marginal_distribution[i]->nb_value;
      }
      else {
        nb_value[i] = I_DEFAULT;
      }
    }

    ihsmarkov = new HiddenSemiMarkov(itype , nb_state , nb_variable , nb_value);

    // initialization of the Markov chain parameters

    ihsmarkov->init(left_right , 0.);

    // initialization of the state occupancy distributions

    if (occupancy_mean == D_DEFAULT) {
      occupancy_mean = MAX(length_distribution->mean , OCCUPANCY_MEAN);
    }

    ihsmarkov->sojourn_type = new state_sojourn_type[nb_state];
    ihsmarkov->state_process->absorption = new double[nb_state];
    ihsmarkov->state_process->sojourn_time = new DiscreteParametric*[nb_state];
    ihsmarkov->forward = new Forward*[nb_state];

    for (i = 0;i < nb_state;i++) {
      if (ihsmarkov->stype[i] != ABSORBING) {
        ihsmarkov->sojourn_type[i] = SEMI_MARKOVIAN;
        ihsmarkov->state_process->absorption[i] = 0.;
        proba = 1. / occupancy_mean;
        ihsmarkov->state_process->sojourn_time[i] = new DiscreteParametric(NEGATIVE_BINOMIAL , 1 ,
                                                                           I_DEFAULT , 1. , proba ,
                                                                           OCCUPANCY_THRESHOLD);

        if (ihsmarkov->stype[i] == RECURRENT) {
          ihsmarkov->forward[i] = new Forward(*(ihsmarkov->state_process->sojourn_time[i]) ,
                                              ihsmarkov->state_process->sojourn_time[i]->alloc_nb_value);
        }
        else {
          ihsmarkov->forward[i] = NULL;
        }
      }

      else {
        ihsmarkov->sojourn_type[i] = MARKOVIAN;
        ihsmarkov->state_process->absorption[i] = 1.;
        ihsmarkov->state_process->sojourn_time[i] = NULL;
        ihsmarkov->forward[i] = NULL;
      }
    }

    // initialization of the observation distributions

    for (i = 0;i < ihsmarkov->nb_output_process;i++) {
      if (ihsmarkov->categorical_process[i]) {
        ihsmarkov->categorical_process[i]->init();
      }

      else if (ihsmarkov->discrete_parametric_process[i]) {
        ihsmarkov->discrete_parametric_process[i]->init();
      }

      else {
        mean = mean_computation(i);
        variance = variance_computation(i , mean);

        ihsmarkov->continuous_parametric_process[i]->init(GAUSSIAN , min_value[i] , max_value[i] ,
                                                          mean , variance);
      }
    }

    hsmarkov = hidden_semi_markov_estimation(error , display , *ihsmarkov , poisson_geometric ,
                                             common_dispersion , estimator , counting_flag ,
                                             state_sequence , nb_iter , mean_estimator);
    delete ihsmarkov;
  }

  return hsmarkov;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Estimation of a hidden semi-Markov chain using the MCEM algorithm.
 *
 *  \param[in] error                  reference on a StatError object,
 *  \param[in] display                flag for displaying estimation intermediate results,
 *  \param[in] ihsmarkov              initial hidden semi-Markov chain,
 *  \param[in] poisson_geometric      flag on the estimation of Poisson geometric state occupancy distributions,
 *  \param[in] common_dispersion      flag common dispersion parameter (continuous observation processes),
 *  \param[in] min_nb_state_sequence  minimum number of generated sequences,
 *  \param[in] max_nb_state_sequence  maximum number of generated sequences,
 *  \param[in] parameter              parameter for defining the number of generated sequences,
 *  \param[in] estimator              estimator type for the reestimation of the state occupancy distributions
 *  \param[in]                        (complete or partial likelihood),
 *  \param[in] counting_flag          flag on the computation of the counting distributions,
 *  \param[in] state_sequence         flag on the computation of the restored state sequences,
 *  \param[in] nb_iter                number of iterations.
 *
 *  \return                           HiddenSemiMarkov object.
 */
/*--------------------------------------------------------------*/

HiddenSemiMarkov* MarkovianSequences::hidden_semi_markov_stochastic_estimation(StatError &error , bool display ,
                                                                               const HiddenSemiMarkov &ihsmarkov ,
                                                                               bool poisson_geometric , bool common_dispersion ,
                                                                               int min_nb_state_sequence ,
                                                                               int max_nb_state_sequence , double parameter ,
                                                                               censoring_estimator estimator ,
                                                                               bool counting_flag , bool state_sequence ,
                                                                               int nb_iter) const

{
  bool status;
  int i , j , k , m , n;
  int max_nb_value , iter , nb_state_sequence , state_occupancy , nb_likelihood_decrease ,
      *occupancy_nb_value , *state_seq , *pstate , ***state_sequence_count , nb_element , **pioutput;
  double likelihood = D_INF , previous_likelihood , occupancy_likelihood , observation_likelihood ,
         min_likelihood , obs_product , residual , **observation , *norm , *state_norm , **forward1 ,
         **state_in , *backward , *cumul_backward , *occupancy_survivor , *censored_occupancy_survivor ,
         diff , variance , **mean_direction , concentration , global_mean_direction , **proutput;
  Distribution *weight;
  DiscreteParametric *occupancy;
  ChainReestimation<double> *chain_reestim;
  Reestimation<double> *bcomplete_run , *censored_run , **complete_run , **final_run , **initial_run ,
                       **single_run , ***observation_reestim;
  HiddenSemiMarkov *hsmarkov;
  SemiMarkovData *seq;
  const Reestimation<double> *prun[3];

# ifdef DEBUG
  double sum;
# endif


  hsmarkov = NULL;
  error.init();

  // test number of values for each variable

  status = false;
  for (i = 0;i < nb_variable;i++) {
    if (max_value[i] > min_value[i]) {
      status = true;
      break;
    }
  }

  if (!status) {
    error.update(STAT_error[STATR_VARIABLE_NB_VALUE]);
  }

  for (i = 0;i < nb_variable;i++) {
    if ((type[i] != INT_VALUE) && (type[i] != REAL_VALUE) && (type[i] != STATE)) {
      status = false;
      ostringstream error_message , correction_message;
      error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                    << STAT_error[STATR_VARIABLE_TYPE];
      correction_message << STAT_variable_word[INT_VALUE] << " or "
                         << STAT_variable_word[REAL_VALUE];
      error.correction_update((error_message.str()).c_str() , (correction_message.str()).c_str());
    }
  }

  if (ihsmarkov.nb_output_process != nb_variable) {
    status = false;
    error.update(STAT_error[STATR_NB_OUTPUT_PROCESS]);
  }

  else {
    for (i = 0;i < nb_variable;i++) {
      if ((ihsmarkov.categorical_process[i]) || (ihsmarkov.discrete_parametric_process[i])) {
        if (type[i] == REAL_VALUE) {
          status = false;
          ostringstream error_message;
          error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                        << STAT_error[STATR_VARIABLE_TYPE];
          error.correction_update((error_message.str()).c_str() , STAT_variable_word[INT_VALUE]);
        }

        else {
          if (min_value[i] < 0) {
            status = false;
            ostringstream error_message;
            error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                          << STAT_error[STATR_POSITIVE_MIN_VALUE];
            error.update((error_message.str()).c_str());
          }

          if (!marginal_distribution[i]) {
            status = false;
            ostringstream error_message;
            error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                          << STAT_error[STATR_MARGINAL_FREQUENCY_DISTRIBUTION];
            error.update((error_message.str()).c_str());
          }

          else {
            if (((ihsmarkov.categorical_process[i]) &&
                 (ihsmarkov.categorical_process[i]->nb_value != marginal_distribution[i]->nb_value)) ||
                ((ihsmarkov.discrete_parametric_process[i]) &&
                 (ihsmarkov.discrete_parametric_process[i]->nb_value < marginal_distribution[i]->nb_value))) {
              status = false;
              ostringstream error_message;
              error_message << STAT_label[STATL_OUTPUT_PROCESS] << " " << i + 1 << ": "
                            << STAT_error[STATR_NB_OUTPUT];
              error.update((error_message.str()).c_str());
            }

            else if ((ihsmarkov.categorical_process[i]) && (!characteristics[i])) {
              for (j = 0;j < marginal_distribution[i]->nb_value;j++) {
                if (marginal_distribution[i]->frequency[j] == 0) {
                  status = false;
                  ostringstream error_message;
                  error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                                << STAT_error[STATR_MISSING_VALUE] << " " << j;
                  error.update((error_message.str()).c_str());
                }
              }
            }
          }
        }
      }

      else if ((ihsmarkov.continuous_parametric_process[i]) &&
               (ihsmarkov.continuous_parametric_process[i]->ident == LINEAR_MODEL) &&
               (ihsmarkov.nb_component < ihsmarkov.nb_state)) {
        status = false;
        error.update(SEQ_error[SEQR_MODEL_STRUCTURE]);
      }
    }
  }

  if ((min_nb_state_sequence < 1) || (min_nb_state_sequence > max_nb_state_sequence)) {
    status = false;
    error.update(SEQ_error[SEQR_MIN_NB_STATE_SEQUENCE]);
  }

  if ((nb_iter != I_DEFAULT) && (nb_iter < 1)) {
    status = false;
    error.update(STAT_error[STATR_NB_ITERATION]);
  }

  if (status) {
    if (max_length > COUNTING_MAX_LENGTH) {
      counting_flag = false;
    }

    // construction of the hidden semi-Markov chain

    hsmarkov = new HiddenSemiMarkov(ihsmarkov , false , (int)(max_length * SAMPLE_NB_VALUE_COEFF));

    if (hsmarkov->type == EQUILIBRIUM) {
      for (i = 0;i < hsmarkov->nb_state;i++) {
        hsmarkov->initial[i] = 1. / (double)hsmarkov->nb_state;
      }
    }

    if (common_dispersion) {
      for (i = 0;i < hsmarkov->nb_output_process;i++) {
        if (hsmarkov->continuous_parametric_process[i]) {
          hsmarkov->continuous_parametric_process[i]->tied_dispersion = true;
        }
      }
    }

#   ifdef DEBUG
    cout << *hsmarkov;
#   endif

    // construction of the data structures of the algorithm

    observation = new double*[max_length];
    for (i = 0;i < max_length;i++) {
      observation[i] = new double[hsmarkov->nb_state];
    }

    norm = new double[max_length];
    state_norm = new double[hsmarkov->nb_state];

    forward1 = new double*[max_length];
    for (i = 0;i < max_length;i++) {
      forward1[i] = new double[hsmarkov->nb_state];
    }

    state_in = new double*[max_length - 1];
    for (i = 0;i < max_length - 1;i++) {
      state_in[i] = new double[hsmarkov->nb_state];
    }

    backward = new double[max_length + 1];
    cumul_backward = new double[max_length + 1];

    state_seq = new int[max_length];

    chain_reestim = new ChainReestimation<double>(hsmarkov->type , hsmarkov->nb_state , hsmarkov->nb_state);

    occupancy_nb_value = new int[hsmarkov->nb_state];
    complete_run = new Reestimation<double>*[hsmarkov->nb_state];
    final_run = new Reestimation<double>*[hsmarkov->nb_state];
    if (hsmarkov->type == EQUILIBRIUM) {
      initial_run = new Reestimation<double>*[hsmarkov->nb_state];
      single_run = new Reestimation<double>*[hsmarkov->nb_state];
    }

    for (i = 0;i < hsmarkov->nb_state;i++) {
      switch (hsmarkov->sojourn_type[i]) {

      case SEMI_MARKOVIAN : {
        occupancy_nb_value[i] = MIN(hsmarkov->state_process->sojourn_time[i]->alloc_nb_value ,
                                    max_length + 1);

        complete_run[i] = new Reestimation<double>(occupancy_nb_value[i]);
        final_run[i] = new Reestimation<double>(occupancy_nb_value[i]);
        if (hsmarkov->type == EQUILIBRIUM) {
          initial_run[i] = new Reestimation<double>(occupancy_nb_value[i]);
          single_run[i] = new Reestimation<double>(occupancy_nb_value[i]);
        }
        break;
      }

      case MARKOVIAN : {
        complete_run[i] = NULL;
        final_run[i] = NULL;
        if (hsmarkov->type == EQUILIBRIUM) {
          initial_run[i] = NULL;
          single_run[i] = NULL;
        }
        break;
      }
      }
    }

    max_nb_value = 0;
    for (i = 0;i < hsmarkov->nb_state;i++) {
      if ((hsmarkov->sojourn_type[i] == SEMI_MARKOVIAN) && (occupancy_nb_value[i] > max_nb_value)) {
        max_nb_value = occupancy_nb_value[i];
      }
    }

    if (estimator != PARTIAL_LIKELIHOOD) {
      occupancy_survivor = new double[max_nb_value];
      censored_occupancy_survivor = new double[max_nb_value + 1];
    }

    observation_reestim = new Reestimation<double>**[hsmarkov->nb_output_process];
    for (i = 0;i < hsmarkov->nb_output_process;i++) {
      if ((marginal_distribution[i]) && ((!(hsmarkov->continuous_parametric_process[i])) ||
           ((hsmarkov->continuous_parametric_process[i]->ident != LINEAR_MODEL) &&
            (hsmarkov->continuous_parametric_process[i]->ident != AUTOREGRESSIVE_MODEL)))) {
        observation_reestim[i] = new Reestimation<double>*[hsmarkov->nb_state];
        for (j = 0;j < hsmarkov->nb_state;j++) {
          observation_reestim[i][j] = new Reestimation<double>(marginal_distribution[i]->nb_value);
        }
      }

      else {
        observation_reestim[i] = NULL;
      }
    }

    for (i = 0;i < hsmarkov->nb_output_process;i++) {
      if ((hsmarkov->continuous_parametric_process[i]) &&
          (hsmarkov->continuous_parametric_process[i]->ident == VON_MISES)) {
        break;
      }
    }

    if (i < hsmarkov->nb_output_process) {
      mean_direction = new double*[hsmarkov->nb_state];
      for (i = 0;i < hsmarkov->nb_state;i++) {
        mean_direction[i] = new double[4];
      }
    }
    else {
      mean_direction = NULL;
    }

    for (i = 0;i < hsmarkov->nb_output_process;i++) {
      if ((!marginal_distribution[i]) || ((hsmarkov->continuous_parametric_process[i]) &&
           ((hsmarkov->continuous_parametric_process[i]->ident == LINEAR_MODEL) ||
            (hsmarkov->continuous_parametric_process[i]->ident == AUTOREGRESSIVE_MODEL)))) {
        break;
      }
    }

    if (i < hsmarkov->nb_output_process) {
      state_sequence_count = new int**[nb_sequence];
      for (i = 0;i < nb_sequence;i++) {
        state_sequence_count[i] = new int*[length[i]];
        for (j = 0;j < length[i];j++) {
          state_sequence_count[i][j] = new int[hsmarkov->nb_state];
        }
      }
    }
    else {
      state_sequence_count = NULL;
    }

    for (i = 0;i < hsmarkov->nb_output_process;i++) {
      if ((!marginal_distribution[i]) || ((hsmarkov->continuous_parametric_process[i]) &&
           (hsmarkov->continuous_parametric_process[i]->ident == AUTOREGRESSIVE_MODEL))) {
        break;
      }
    }

    pioutput = new int*[nb_variable];
    proutput = new double*[nb_variable];

    iter = 0;
    nb_likelihood_decrease = 0;

    do {
      previous_likelihood = likelihood;
      likelihood = 0.;

      // computation of the number of generated state sequences

      if (min_nb_state_sequence + (int)::round(parameter * iter) < max_nb_state_sequence) {
        nb_state_sequence = min_nb_state_sequence + (int)::round(parameter * iter);
      }
      else {
        nb_state_sequence = max_nb_state_sequence;
      }

/*      nb_state_sequence = max_nb_state_sequence - (int)::round((max_nb_state_sequence - min_nb_state_sequence) *
                          exp(-parameter * iter)); */

      iter++;

      // initialization of the reestimation quantities

      chain_reestim->init();

      for (i = 0;i < hsmarkov->nb_state;i++) {
        if (hsmarkov->sojourn_type[i] == SEMI_MARKOVIAN) {
          for (j = 0;j < occupancy_nb_value[i];j++) {
            complete_run[i]->frequency[j] = 0.;
          }

          for (j = 0;j < occupancy_nb_value[i];j++) {
            final_run[i]->frequency[j] = 0.;
          }

          if (hsmarkov->type == EQUILIBRIUM) {
            for (j = 0;j < occupancy_nb_value[i];j++) {
              initial_run[i]->frequency[j] = 0.;
            }

            for (j = 0;j < occupancy_nb_value[i];j++) {
              single_run[i]->frequency[j] = 0.;
            }
          }
        }
      }

      for (i = 0;i < hsmarkov->nb_output_process;i++) {
        if (observation_reestim[i]) {
          for (j = 0;j < hsmarkov->nb_state;j++) {
            for (k = 0;k < marginal_distribution[i]->nb_value;k++) {
              observation_reestim[i][j]->frequency[k] = 0.;
            }
          }
        }
      }

      if (state_sequence_count) {
        for (i = 0;i < nb_sequence;i++) {
          for (j = 0;j < length[i];j++) {
            for (k = 0;k < hsmarkov->nb_state;k++) {
              state_sequence_count[i][j][k] = 0;
            }
          }
        }
      }

      for (i = 0;i < nb_sequence;i++) {

        // forward recurrence

        for (j = 0;j < nb_variable;j++) {
          switch (type[j]) {
          case INT_VALUE :
            pioutput[j] = int_sequence[i][j];
            break;
          case REAL_VALUE :
            proutput[j] = real_sequence[i][j];
            break;
          }
        }

        for (j = 0;j < length[i];j++) {
          norm[j] = 0.;

          for (k = 0;k < hsmarkov->nb_state;k++) {

            // computation of the observation probabilities

            observation[j][k] = 1.;
            for (m = 0;m < hsmarkov->nb_output_process;m++) {
              if (hsmarkov->categorical_process[m]) {
                observation[j][k] *= hsmarkov->categorical_process[m]->observation[k]->mass[*pioutput[m]];
              }

              else if (hsmarkov->discrete_parametric_process[m]) {
                observation[j][k] *= hsmarkov->discrete_parametric_process[m]->observation[k]->mass[*pioutput[m]];
              }

              else {
                if (((hsmarkov->continuous_parametric_process[m]->ident == GAMMA) ||
                     (hsmarkov->continuous_parametric_process[m]->ident == ZERO_INFLATED_GAMMA)) && (min_value[m] < min_interval[m] / 2)) {
                  switch (type[m]) {
                  case INT_VALUE :
                    observation[j][k] *= hsmarkov->continuous_parametric_process[m]->observation[k]->mass_computation(*pioutput[m] , *pioutput[m] + min_interval[m]);
                    break;
                  case REAL_VALUE :
                    observation[j][k] *= hsmarkov->continuous_parametric_process[m]->observation[k]->mass_computation(*proutput[m] , *proutput[m] + min_interval[m]);
                    break;
                  }
                }

                else if (hsmarkov->continuous_parametric_process[m]->ident == LINEAR_MODEL) {
                  switch (type[m]) {
                  case INT_VALUE :
                    residual = *pioutput[m] - (hsmarkov->continuous_parametric_process[m]->observation[k]->intercept +
                                hsmarkov->continuous_parametric_process[m]->observation[k]->slope *
                                (index_param_type == IMPLICIT_TYPE ? j : index_parameter[i][j]));
                    break;
                  case REAL_VALUE :
                    residual = *proutput[m] - (hsmarkov->continuous_parametric_process[m]->observation[k]->intercept +
                                hsmarkov->continuous_parametric_process[m]->observation[k]->slope *
                                (index_param_type == IMPLICIT_TYPE ? j : index_parameter[i][j]));
                    break;
                  }

                  observation[j][k] *= hsmarkov->continuous_parametric_process[m]->observation[k]->mass_computation(residual , residual);
                }

                else if (hsmarkov->continuous_parametric_process[m]->ident == AUTOREGRESSIVE_MODEL) {
                  if (j == 0) {
                    switch (type[m]) {
                    case INT_VALUE :
                      residual = *pioutput[m] - hsmarkov->continuous_parametric_process[m]->observation[k]->location;
                      break;
                    case REAL_VALUE :
                      residual = *proutput[m] - hsmarkov->continuous_parametric_process[m]->observation[k]->location;
                      break;
                    }
                  }

                  else {
                    switch (type[m]) {
                    case INT_VALUE :
                      residual = *pioutput[m] - (hsmarkov->continuous_parametric_process[m]->observation[k]->location +
                                  hsmarkov->continuous_parametric_process[m]->observation[k]->autoregressive_coeff *
                                  (*(pioutput[m] - 1) - hsmarkov->continuous_parametric_process[m]->observation[k]->location));
                      break;
                    case REAL_VALUE :
                      residual = *proutput[m] - (hsmarkov->continuous_parametric_process[m]->observation[k]->location +
                                  hsmarkov->continuous_parametric_process[m]->observation[k]->autoregressive_coeff *
                                  (*(proutput[m] - 1) - hsmarkov->continuous_parametric_process[m]->observation[k]->location));
                      break;
                    }
                  }

                  observation[j][k] *= hsmarkov->continuous_parametric_process[m]->observation[k]->mass_computation(residual , residual);
                }

                else {
                  switch (type[m]) {
                  case INT_VALUE :
                    observation[j][k] *= hsmarkov->continuous_parametric_process[m]->observation[k]->mass_computation(*pioutput[m] - min_interval[m] / 2 , *pioutput[m] + min_interval[m] / 2);
                    break;
                  case REAL_VALUE :
                    observation[j][k] *= hsmarkov->continuous_parametric_process[m]->observation[k]->mass_computation(*proutput[m] - min_interval[m] / 2 , *proutput[m] + min_interval[m] / 2);
                    break;
                  }
                }
              }
            }

            switch (hsmarkov->sojourn_type[k]) {

            // case semi-Markovian state

            case SEMI_MARKOVIAN : {
              if (j == 0) {
                state_norm[k] = hsmarkov->initial[k];
              }
              else {
                state_norm[k] += state_in[j - 1][k] - forward1[j - 1][k];
              }
              state_norm[k] *= observation[j][k];

              norm[j] += state_norm[k];
              break;
            }

            // case Markovian state

            case MARKOVIAN : {
              if (j == 0) {
                forward1[j][k] = hsmarkov->initial[k];
              }
              else {
                forward1[j][k] = state_in[j - 1][k];
              }
              forward1[j][k] *= observation[j][k];

              norm[j] += forward1[j][k];
              break;
            }
            }
          }

          if (norm[j] > 0.) {
            for (k = 0;k < hsmarkov->nb_state;k++) {
              switch (hsmarkov->sojourn_type[k]) {
              case SEMI_MARKOVIAN :
                state_norm[k] /= norm[j];
                break;
              case MARKOVIAN :
                forward1[j][k] /= norm[j];
                break;
              }
            }

            likelihood += log(norm[j]);
          }

          else {
            likelihood = D_INF;
            break;
          }

          for (k = 0;k < hsmarkov->nb_state;k++) {

            // case semi-Markovian state

            if (hsmarkov->sojourn_type[k] == SEMI_MARKOVIAN) {
              occupancy = hsmarkov->state_process->sojourn_time[k];
              obs_product = 1.;
              forward1[j][k] = 0.;

              if (j < length[i] - 1) {
                for (m = 1;m <= MIN(j + 1 , occupancy->nb_value - 1);m++) {
                  obs_product *= observation[j - m + 1][k] / norm[j - m + 1];
                  if (obs_product == 0.) {
                    break;
                  }

                  if (m < j + 1) {
                    forward1[j][k] += obs_product * occupancy->mass[m] * state_in[j - m][k];
                  }

                  else {
                    switch (hsmarkov->type) {
                    case ORDINARY :
                      forward1[j][k] += obs_product * occupancy->mass[m] * hsmarkov->initial[k];
                      break;
                    case EQUILIBRIUM :
                      forward1[j][k] += obs_product * hsmarkov->forward[k]->mass[m] * hsmarkov->initial[k];
                      break;
                    }
                  }
                }
              }

              else {
                for (m = 1;m <= MIN(j + 1 , occupancy->nb_value - 1);m++) {
                  obs_product *= observation[j - m + 1][k] / norm[j - m + 1];
                  if (obs_product == 0.) {
                    break;
                  }

                  if (m < j + 1) {
                    forward1[j][k] += obs_product * (1. - occupancy->cumul[m - 1]) * state_in[j - m][k];
                  }

                  else {
                    switch (hsmarkov->type) {
                    case ORDINARY :
                      forward1[j][k] += obs_product * (1. - occupancy->cumul[m - 1]) *
                                        hsmarkov->initial[k];
                      break;
                    case EQUILIBRIUM :
                      forward1[j][k] += obs_product * (1. - hsmarkov->forward[k]->cumul[m - 1]) *
                                        hsmarkov->initial[k];
                      break;
                    }
                  }
                }
              }
            }
          }

          if (j < length[i] - 1) {
            for (k = 0;k < hsmarkov->nb_state;k++) {
              state_in[j][k] = 0.;
              for (m = 0;m < hsmarkov->nb_state;m++) {
                state_in[j][k] += hsmarkov->transition[m][k] * forward1[j][m];
              }
            }
          }

          for (k = 0;k < nb_variable;k++) {
            switch (type[k]) {
            case INT_VALUE :
              pioutput[k]++;
              break;
            case REAL_VALUE :
              proutput[k]++;
              break;
            }
          }
        }

        if (likelihood == D_INF) {
          break;
        }

#       ifdef DEBUG
        for (j = 0;j < length[i];j++) {
          cout << j << " : ";
          for (k = 0;k < hsmarkov->nb_state;k++) {
            cout << forward1[j][k] << " ";
          }
          cout << endl;
        }
        cout << endl;
#       endif

        // backward passes

        for (j = 0;j < nb_state_sequence;j++) {
          k = length[i] - 1;
          pstate = state_seq + k;
          for (m = 0;m < nb_variable;m++) {
            if (type[m] == INT_VALUE) {
              pioutput[m] = int_sequence[i][m] + k;
            }
          }

          cumul_computation(hsmarkov->nb_state , forward1[k] , cumul_backward);
          *pstate = cumul_method(hsmarkov->nb_state , cumul_backward);

          // accumulation of the reestimation quantities of the observation distributions

          for (m = 0;m < hsmarkov->nb_output_process;m++) {
            if (observation_reestim[m]) {
              (observation_reestim[m][*pstate]->frequency[*pioutput[m]])++;
            }
          }

          if (state_sequence_count) {
            (state_sequence_count[i][k][*pstate])++;
          }

          do {

            // case semi-Markovian state

            if (hsmarkov->sojourn_type[*pstate] == SEMI_MARKOVIAN) {
              occupancy = hsmarkov->state_process->sojourn_time[*pstate];
              obs_product = 1.;

              if (k < length[i] - 1) {
                for (m = 1;m <= MIN(k + 1 , occupancy->nb_value - 1);m++) {
                  obs_product *= observation[k - m + 1][*pstate] / norm[k - m + 1];
                  if (obs_product == 0.) {
                    break;
                  }

                  if (m < k + 1) {
                    backward[m] = obs_product * occupancy->mass[m] * state_in[k - m][*pstate] /
                                  forward1[k][*pstate];
                  }

                  else {
                    switch (hsmarkov->type) {
                    case ORDINARY :
                      backward[m] = obs_product * occupancy->mass[m] * hsmarkov->initial[*pstate] /
                                    forward1[k][*pstate];
                      break;
                    case EQUILIBRIUM :
                      backward[m] = obs_product * hsmarkov->forward[*pstate]->mass[m] * hsmarkov->initial[*pstate] /
                                    forward1[k][*pstate];
                      break;
                    }
                  }
                }
              }

              else {
                for (m = 1;m <= MIN(k + 1 , occupancy->nb_value - 1);m++) {
                  obs_product *= observation[k - m + 1][*pstate] / norm[k - m + 1];
                  if (obs_product == 0.) {
                    break;
                  }

                  if (m < k + 1) {
                    backward[m] = obs_product * (1. - occupancy->cumul[m - 1]) * state_in[k - m][*pstate] /
                                  forward1[k][*pstate];
                  }

                  else {
                    switch (hsmarkov->type) {
                    case ORDINARY :
                      backward[m] = obs_product * (1. - occupancy->cumul[m - 1]) *
                                    hsmarkov->initial[*pstate] / forward1[k][*pstate];
                      break;
                    case EQUILIBRIUM :
                      backward[m] = obs_product * (1. - hsmarkov->forward[*pstate]->cumul[m - 1]) *
                                    hsmarkov->initial[*pstate] / forward1[k][*pstate];
                      break;
                    }
                  }
                }
              }

              cumul_computation(m - 1 , backward + 1 , cumul_backward);
              state_occupancy = 1 + cumul_method(m - 1 , cumul_backward);

#             ifdef DEBUG
              sum = 0.;
              for (n = 1;n < m;n++) {
                sum += backward[n];
              }
              if ((sum < 1. - DOUBLE_ERROR) || (sum > 1. + DOUBLE_ERROR)) {
                cout << "\nERROR: " << k << " " << sum << endl;
              }
#             endif

              // accumulation of the reestimation quantities of the state occupancy distributions

              if (k < length[i] - 1) {
                if (state_occupancy < k + 1) {
                  (complete_run[*pstate]->frequency[state_occupancy])++;
                }

                else {
                  switch (hsmarkov->type) {
                  case ORDINARY :
                    (complete_run[*pstate]->frequency[state_occupancy])++;
                    break;
                  case EQUILIBRIUM :
                    (initial_run[*pstate]->frequency[state_occupancy])++;
                    break;
                  }
                }
              }

              else {
                if (state_occupancy < k + 1) {
                  (final_run[*pstate]->frequency[state_occupancy])++;
                }

                else {
                  switch (hsmarkov->type) {
                  case ORDINARY :
                    (final_run[*pstate]->frequency[state_occupancy])++;
                    break;
                  case EQUILIBRIUM :
                    (single_run[*pstate]->frequency[state_occupancy])++;
                    break;
                  }
                }
              }

              for (m = 1;m < state_occupancy;m++) {
                pstate--;
                *pstate = *(pstate + 1);

                // accumulation of the reestimation quantities of the observation distributions

                for (n = 0;n < hsmarkov->nb_output_process;n++) {
                  if (observation_reestim[n]) {
                    (observation_reestim[n][*pstate]->frequency[*--pioutput[n]])++;
                  }
                }

                if (state_sequence_count) {
                  (state_sequence_count[i][k - m + 1][*pstate])++;
                }
              }
              k -= (state_occupancy - 1);

              if (k == 0) {
                break;
              }
            }

            k--;
            for (m = 0;m < hsmarkov->nb_state;m++) {
              backward[m] = hsmarkov->transition[m][*pstate] * forward1[k][m] / state_in[k][*pstate];
            }
            cumul_computation(hsmarkov->nb_state , backward , cumul_backward);
            *--pstate = cumul_method(hsmarkov->nb_state , cumul_backward);

            // accumulation of the reestimation quantities of the transition probabilities and
            // the observation distributions

            (chain_reestim->transition[*pstate][*(pstate + 1)])++;

            for (m = 0;m < hsmarkov->nb_output_process;m++) {
              if (observation_reestim[m]) {
                (observation_reestim[m][*pstate]->frequency[*--pioutput[m]])++;
              }
            }

            if (state_sequence_count) {
              (state_sequence_count[i][k][*pstate])++;
            }

#           ifdef DEBUG
            sum = 0.;
            for (m = 0;m < nb_state;m++) {
              sum += backward[m];
            }
            if ((sum < 1. - DOUBLE_ERROR) || (sum > 1. + DOUBLE_ERROR)) {
              cout << "\nERROR: " << k << " " << sum << endl;
            }
#           endif

          }
          while (k > 0);

          // accumulation of the reestimation quantities of the initial probabilities

          if (hsmarkov->type == ORDINARY) {
            (chain_reestim->initial[*pstate])++;
          }
        }
      }

      if (likelihood != D_INF) {
        if (likelihood < previous_likelihood) {
          nb_likelihood_decrease++;
        }
        else {
          nb_likelihood_decrease = 0;
        }

        // reestimation of the initial probabilities

        if (hsmarkov->type == ORDINARY) {
          reestimation(hsmarkov->nb_state , chain_reestim->initial ,
                       hsmarkov->initial , MIN_PROBABILITY , false);
        }

        // reestimation of the transition probabilities

        for (i = 0;i < hsmarkov->nb_state;i++) {
          reestimation(hsmarkov->nb_state , chain_reestim->transition[i] ,
                       hsmarkov->transition[i] , MIN_PROBABILITY , false);
        }

        // reestimation of the state occupancy distributions

        min_likelihood = 0.;

        for (i = 0;i < hsmarkov->nb_state;i++) {
          if (hsmarkov->sojourn_type[i] == SEMI_MARKOVIAN) {
            occupancy = hsmarkov->state_process->sojourn_time[i];

            complete_run[i]->nb_value_computation();
            complete_run[i]->offset_computation();
            complete_run[i]->nb_element_computation();

#           ifdef DEBUG
            cout << "\n" << STAT_label[STATL_STATE] << " " << i << " ";

            complete_run[i]->max_computation();
            complete_run[i]->mean_computation();
            complete_run[i]->variance_computation();

            complete_run[i]->print(cout);
#           endif

            if ((iter > STOCHASTIC_EXPLORATION_NB_ITER) && (estimator == COMPLETE_LIKELIHOOD)) {
              final_run[i]->nb_value_computation();
              final_run[i]->offset_computation();
              final_run[i]->nb_element_computation();

              switch (hsmarkov->type) {
              case ORDINARY : {
                if (final_run[i]->nb_element > 0.) {
                  complete_run[i]->state_occupancy_estimation(final_run[i] , complete_run[i] ,
                                                              occupancy_survivor ,
                                                              censored_occupancy_survivor , false);
                }
                break;
              }

              case EQUILIBRIUM : {
                initial_run[i]->nb_value_computation();
                initial_run[i]->offset_computation();
                initial_run[i]->nb_element_computation();

                single_run[i]->nb_value_computation();
                single_run[i]->offset_computation();
                single_run[i]->nb_element_computation();

                prun[0] = complete_run[i];
                prun[1] = complete_run[i];
                bcomplete_run = new Reestimation<double>(2 , prun);

                prun[0] = initial_run[i];
                prun[1] = final_run[i];
                prun[2] = single_run[i];
                censored_run = new Reestimation<double>(3 , prun);

#               ifdef DEBUG
                censored_run->print(cout);
#               endif

                bcomplete_run->state_occupancy_estimation(censored_run , complete_run[i] ,
                                                          occupancy_survivor ,
                                                          censored_occupancy_survivor , false);
                delete bcomplete_run;
                delete censored_run;
                break;
              }
              }

              if ((hsmarkov->type == EQUILIBRIUM) || (final_run[i]->nb_element > 0.)) {
                complete_run[i]->nb_value_computation();
                complete_run[i]->offset_computation();
                complete_run[i]->nb_element_computation();

#               ifdef DEBUG
                complete_run[i]->max_computation();
                complete_run[i]->mean_computation();
                complete_run[i]->variance_computation();

                complete_run[i]->print(cout);
#               endif

              }
            }

            complete_run[i]->max_computation();
            complete_run[i]->mean_computation();
            complete_run[i]->variance_computation();

            if (iter <= EXPLORATION_NB_ITER) {
              occupancy_likelihood = complete_run[i]->parametric_estimation(occupancy , 1 , true ,
                                                                            OCCUPANCY_THRESHOLD , poisson_geometric);
            }
            else {
              occupancy_likelihood = complete_run[i]->type_parametric_estimation(occupancy , 1 , true ,
                                                                                 OCCUPANCY_THRESHOLD , poisson_geometric);
            }

#           ifdef DEBUG
            if (i == 1) {
              occupancy->print(cout);
            }
#           endif

            if (occupancy_likelihood == D_INF) {
              min_likelihood = D_INF;
            }
            else {
              occupancy->computation(complete_run[i]->nb_value , OCCUPANCY_THRESHOLD);
              if (hsmarkov->type == EQUILIBRIUM) {
                hsmarkov->forward[i]->copy(*occupancy);
                hsmarkov->forward[i]->computation(*occupancy);
              }
            }

#           ifdef DEBUG
            cout << STAT_word[STATW_STATE] << " " << i << " " << SEQ_word[SEQW_OCCUPANCY_DISTRIBUTION] << endl;
            occupancy->ascii_print(cout);
#           endif

          }
        }

        if (hsmarkov->type == EQUILIBRIUM) {
          hsmarkov->initial_probability_computation();
        }

        // reestimation of the observation distributions

        for (i = 0;i < hsmarkov->nb_output_process;i++) {
          if (hsmarkov->categorical_process[i]) {
            for (j = 0;j < hsmarkov->nb_state;j++) {
              reestimation(marginal_distribution[i]->nb_value , observation_reestim[i][j]->frequency ,
                           hsmarkov->categorical_process[i]->observation[j]->mass ,
                           MIN_PROBABILITY , false);
            }
          }

          else if (observation_reestim[i]) {
            for (j = 0;j < hsmarkov->nb_state;j++) {
              observation_reestim[i][j]->nb_value_computation();
              observation_reestim[i][j]->offset_computation();
              observation_reestim[i][j]->nb_element_computation();
              observation_reestim[i][j]->max_computation();
              if ((hsmarkov->discrete_parametric_process[i]) ||
                  (hsmarkov->continuous_parametric_process[i]->ident != ZERO_INFLATED_GAMMA)) {
                observation_reestim[i][j]->mean_computation();
//                observation_reestim[i][j]->variance_computation();
                observation_reestim[i][j]->variance_computation(true);
              }
            }

            if (hsmarkov->discrete_parametric_process[i]) {
              for (j = 0;j < hsmarkov->nb_state;j++) {
                observation_likelihood = observation_reestim[i][j]->type_parametric_estimation(hsmarkov->discrete_parametric_process[i]->observation[j] ,
                                                                                               0 , true , OBSERVATION_THRESHOLD);

                if (observation_likelihood == D_INF) {
                  min_likelihood = D_INF;
                }
                else {
                  hsmarkov->discrete_parametric_process[i]->observation[j]->computation(marginal_distribution[i]->nb_value ,
                                                                                        OBSERVATION_THRESHOLD);

                  if (hsmarkov->discrete_parametric_process[i]->observation[j]->ident == BINOMIAL) {
                    for (k = hsmarkov->discrete_parametric_process[i]->observation[j]->nb_value;k < marginal_distribution[i]->nb_value;k++) {
                      hsmarkov->discrete_parametric_process[i]->observation[j]->mass[k] = 0.;
                    }
                  }
                }
              }
            }

            else {
              switch (hsmarkov->continuous_parametric_process[i]->ident) {

              case GAMMA : {
                for (j = 0;j < hsmarkov->nb_state;j++) {
                  observation_reestim[i][j]->gamma_estimation(hsmarkov->continuous_parametric_process[i]->observation[j] , iter);
                }
                break;
              }

              case ZERO_INFLATED_GAMMA : {
                for (j = 0;j < hsmarkov->nb_state;j++) {
                  observation_reestim[i][j]->zero_inflated_gamma_estimation(hsmarkov->continuous_parametric_process[i]->observation[j] , iter);
                }
                break;
              }

              case GAUSSIAN : {
                for (j = 0;j < hsmarkov->nb_state;j++) {
                  hsmarkov->continuous_parametric_process[i]->observation[j]->location = observation_reestim[i][j]->mean;
                }

                if (common_dispersion) {
                  variance = 0.;
                  nb_element = 0;

                  for (j = 0;j < hsmarkov->nb_state;j++) {
                    for (k = observation_reestim[i][j]->offset;k < observation_reestim[i][j]->nb_value;k++) {
                      diff = k - observation_reestim[i][j]->mean;
                      variance += observation_reestim[i][j]->frequency[k] * diff * diff;
                    }

                    nb_element += observation_reestim[i][j]->nb_element;
                  }

                  variance /= nb_element;
//                  variance /= (nb_element - 1);

                  for (j = 0;j < hsmarkov->nb_state;j++) {
                    hsmarkov->continuous_parametric_process[i]->observation[j]->dispersion = sqrt(variance);
                  }
                }

                else {
                  for (j = 0;j < hsmarkov->nb_state;j++) {
                    hsmarkov->continuous_parametric_process[i]->observation[j]->dispersion = sqrt(observation_reestim[i][j]->variance);
                    if (hsmarkov->continuous_parametric_process[i]->observation[j]->dispersion /
                        hsmarkov->continuous_parametric_process[i]->observation[j]->location < GAUSSIAN_MIN_VARIATION_COEFF) {
                      hsmarkov->continuous_parametric_process[i]->observation[j]->dispersion = hsmarkov->continuous_parametric_process[i]->observation[j]->location * GAUSSIAN_MIN_VARIATION_COEFF;
                    }
                  }
                }
                break;
              }

              case VON_MISES : {
                for (j = 0;j < hsmarkov->nb_state;j++) {
                  observation_reestim[i][j]->mean_direction_computation(mean_direction[j]);
                  hsmarkov->continuous_parametric_process[i]->observation[j]->location = mean_direction[j][3];
                }

                if (common_dispersion) {
                  global_mean_direction = 0.;
                  nb_element = 0;

                  for (j = 0;j < hsmarkov->nb_state;j++) {
                    global_mean_direction += observation_reestim[i][j]->nb_element * mean_direction[j][2];
                    nb_element += observation_reestim[i][j]->nb_element;
                  }
                  concentration = von_mises_concentration_computation(global_mean_direction / nb_element);

                  for (j = 0;j < hsmarkov->nb_state;j++) {
                    hsmarkov->continuous_parametric_process[i]->observation[j]->dispersion = concentration;
                  }
                }

                else {
                  for (j = 0;j < hsmarkov->nb_state;j++) {
                    hsmarkov->continuous_parametric_process[i]->observation[j]->dispersion = von_mises_concentration_computation(mean_direction[j][2]);
                  }
                }
                break;
              }
              }
            }
          }

          else {
            switch (hsmarkov->continuous_parametric_process[i]->ident) {
            case GAMMA :
              gamma_estimation(state_sequence_count , i ,
                               hsmarkov->continuous_parametric_process[i] , iter);
              break;
            case ZERO_INFLATED_GAMMA :
              zero_inflated_gamma_estimation(state_sequence_count , i ,
                                             hsmarkov->continuous_parametric_process[i] , iter);
              break;
            case GAUSSIAN :
              gaussian_estimation(state_sequence_count , i ,
                                  hsmarkov->continuous_parametric_process[i]);
              break;
            case VON_MISES :
              von_mises_estimation(state_sequence_count , i ,
                                   hsmarkov->continuous_parametric_process[i]);
              break;
            case LINEAR_MODEL :
              linear_model_estimation(state_sequence_count , i ,
                                      hsmarkov->continuous_parametric_process[i]);
              break;
            case AUTOREGRESSIVE_MODEL :
              autoregressive_model_estimation(state_sequence_count , i ,
                                              hsmarkov->continuous_parametric_process[i]);
              break;
            }
          }
        }
      }

      if (display) {
        cout << STAT_label[STATL_ITERATION] << " " << iter << "   "
             << SEQ_label[SEQL_OBSERVED_SEQUENCES_LIKELIHOOD] << ": " << likelihood
             << "   (" << nb_state_sequence << ")" << endl;
      }

#     ifdef DEBUG
      if (iter % 5 == 0) {
        cout << *hsmarkov;
      }
#     endif

    }
    while ((likelihood != D_INF) && ((iter <= STOCHASTIC_EXPLORATION_NB_ITER + 2) ||
            ((nb_iter == I_DEFAULT) && (iter < SEMI_MARKOV_NB_ITER) &&
             (((likelihood - previous_likelihood) / -likelihood > SEMI_MARKOV_LIKELIHOOD_DIFF) ||
              (min_likelihood == D_INF) || (nb_likelihood_decrease == 1))) ||
            ((nb_iter != I_DEFAULT) && (iter < nb_iter))));

    if (likelihood != D_INF) {
      if (display) {
        cout << "\n" << iter << " " << STAT_label[STATL_ITERATIONS] << endl;

        if (hsmarkov->type == EQUILIBRIUM) {
          for (i = 0;i < hsmarkov->nb_state;i++) {
            if (single_run[i]->nb_element > 0) {
              cout << "\n" << SEQ_label[SEQL_BIASED] << " " << STAT_label[STATL_STATE] << " " << i
                   << " " << SEQ_label[SEQL_OCCUPANCY] << " "  << STAT_label[STATL_DISTRIBUTION] << endl;
            }
          }
        }
      }

      // reestimation of the initial probabilities

      if (hsmarkov->type == ORDINARY) {
        reestimation(hsmarkov->nb_state , chain_reestim->initial ,
                     hsmarkov->initial , MIN_PROBABILITY , true);
      }

      // reestimation of the transition probabilities

      for (i = 0;i < hsmarkov->nb_state;i++) {
        reestimation(hsmarkov->nb_state , chain_reestim->transition[i] ,
                     hsmarkov->transition[i] , MIN_PROBABILITY , true);
      }

      if (hsmarkov->type == EQUILIBRIUM) {
        hsmarkov->initial_probability_computation();
      }

      for (i = 0;i < hsmarkov->nb_state;i++) {
        if ((hsmarkov->sojourn_type[i] == SEMI_MARKOVIAN) &&
            (hsmarkov->state_process->sojourn_time[i]->mean == 1.)) {
          hsmarkov->sojourn_type[i] = MARKOVIAN;
          delete hsmarkov->state_process->sojourn_time[i];
          hsmarkov->state_process->sojourn_time[i] = NULL;
          delete hsmarkov->forward[i];
          hsmarkov->forward[i] = NULL;
        }
      }

      // reestimation of the categorical observation distributions

      for (i = 0;i < hsmarkov->nb_output_process;i++) {
        if (hsmarkov->categorical_process[i]) {
          for (j = 0;j < hsmarkov->nb_state;j++) {
            reestimation(marginal_distribution[i]->nb_value , observation_reestim[i][j]->frequency ,
                         hsmarkov->categorical_process[i]->observation[j]->mass ,
                         MIN_PROBABILITY , true);
          }
        }

        else if (hsmarkov->discrete_parametric_process[i]) {
          hsmarkov->discrete_parametric_process[i]->nb_value_computation();
        }
      }
    }

    // destruction of the data structures of the algorithm

    for (i = 0;i < max_length;i++) {
      delete [] observation[i];
    }
    delete [] observation;

    delete [] norm;
    delete [] state_norm;

    for (i = 0;i < max_length;i++) {
      delete [] forward1[i];
    }
    delete [] forward1;

    for (i = 0;i < max_length - 1;i++) {
      delete [] state_in[i];
    }
    delete [] state_in;

    delete [] backward;
    delete [] cumul_backward;

    delete [] state_seq;

    delete chain_reestim;

    for (i = 0;i < hsmarkov->nb_state;i++) {
      delete complete_run[i];
    }
    delete [] complete_run;

    for (i = 0;i < hsmarkov->nb_state;i++) {
      delete final_run[i];
    }
    delete [] final_run;

    if (hsmarkov->type == EQUILIBRIUM) {
      for (i = 0;i < hsmarkov->nb_state;i++) {
        delete initial_run[i];
      }
      delete [] initial_run;

      for (i = 0;i < hsmarkov->nb_state;i++) {
        delete single_run[i];
      }
      delete [] single_run;
    }

    delete [] occupancy_nb_value;

    if (estimator != PARTIAL_LIKELIHOOD) {
      delete [] occupancy_survivor;
      delete [] censored_occupancy_survivor;
    }

    for (i = 0;i < hsmarkov->nb_output_process;i++) {
      if (observation_reestim[i]) {
        for (j = 0;j < hsmarkov->nb_state;j++) {
          delete observation_reestim[i][j];
        }
        delete [] observation_reestim[i];
      }
    }
    delete [] observation_reestim;

    if (mean_direction) {
      for (i = 0;i < hsmarkov->nb_state;i++) {
        delete [] mean_direction[i];
      }
      delete [] mean_direction;
    }

    if (state_sequence_count) {
      for (i = 0;i < nb_sequence;i++) {
        for (j = 0;j < length[i];j++) {
          delete [] state_sequence_count[i][j];
        }
        delete [] state_sequence_count[i];
      }
      delete [] state_sequence_count;
    }

    delete [] pioutput;
    delete [] proutput;

    if (likelihood == D_INF) {
      delete hsmarkov;
      hsmarkov = NULL;
      error.update(STAT_error[STATR_ESTIMATION_FAILURE]);
    }

    else {
      if (state_sequence) {
        hsmarkov->semi_markov_data = new SemiMarkovData(*this , ADD_STATE_VARIABLE ,
                                                        (hsmarkov->type == EQUILIBRIUM ? true : false));
        seq = hsmarkov->semi_markov_data;

        for (i = 0;i < hsmarkov->nb_output_process;i++) {
          if (((hsmarkov->discrete_parametric_process[i]) || (hsmarkov->continuous_parametric_process[i])) &&
              (seq->characteristics[i + 1])) {
            delete seq->characteristics[i + 1];
            seq->characteristics[i + 1] = NULL;
          }
        }

        hsmarkov->forward_backward(*seq);

        hsmarkov->create_cumul();
        hsmarkov->log_computation();
        hsmarkov->viterbi(*seq);
        hsmarkov->remove_cumul();

        seq->min_value_computation(0);
        seq->max_value_computation(0);
        seq->build_marginal_frequency_distribution(0);
        seq->build_characteristic(0 , true , (hsmarkov->type == EQUILIBRIUM ? true : false));

        seq->build_transition_count(hsmarkov);
        seq->build_observation_frequency_distribution(hsmarkov->nb_state);
        seq->build_observation_histogram(hsmarkov->nb_state);

        // computation of the state occupancy distributions

        for (i = 0;i < hsmarkov->nb_state;i++) {
          if (hsmarkov->sojourn_type[i] == SEMI_MARKOVIAN) {
            hsmarkov->state_process->sojourn_time[i]->computation((seq->characteristics[0] ? seq->characteristics[0]->sojourn_time[i]->nb_value : 1) ,
                                                                  OCCUPANCY_THRESHOLD);
            if (hsmarkov->stype[i] == RECURRENT) {
              if (hsmarkov->type == ORDINARY) {
                hsmarkov->forward[i]->copy(*(hsmarkov->state_process->sojourn_time[i]));
              }
              hsmarkov->forward[i]->computation(*(hsmarkov->state_process->sojourn_time[i]));
            }
          }
        }

        // computation of the mixtures of observation distributions (weights deduced from the restoration)

        weight = NULL;

        for (i = 0;i < hsmarkov->nb_output_process;i++) {
          if ((hsmarkov->categorical_process[i]) || (hsmarkov->discrete_parametric_process[i]) ||
              ((hsmarkov->continuous_parametric_process[i]) &&
               (hsmarkov->continuous_parametric_process[i]->ident != LINEAR_MODEL) &&
               (hsmarkov->continuous_parametric_process[i]->ident != AUTOREGRESSIVE_MODEL))) {
            weight = seq->weight_computation();
            break;
          }
        }

        for (i = 0;i < hsmarkov->nb_output_process;i++) {
          if (hsmarkov->categorical_process[i]) {
            hsmarkov->categorical_process[i]->restoration_weight = new Distribution(*weight);
            hsmarkov->categorical_process[i]->restoration_mixture = hsmarkov->categorical_process[i]->mixture_computation(hsmarkov->categorical_process[i]->restoration_weight);
          }

          else if (hsmarkov->discrete_parametric_process[i]) {
            for (j = 0;j < hsmarkov->nb_state;j++) {
              hsmarkov->discrete_parametric_process[i]->observation[j]->cumul_computation();
            }

            hsmarkov->discrete_parametric_process[i]->restoration_weight = new Distribution(*weight);
            hsmarkov->discrete_parametric_process[i]->restoration_mixture = hsmarkov->discrete_parametric_process[i]->mixture_computation(hsmarkov->discrete_parametric_process[i]->restoration_weight);
          }

          else if ((hsmarkov->continuous_parametric_process[i]) &&
                   (hsmarkov->continuous_parametric_process[i]->ident != LINEAR_MODEL) &&
                   (hsmarkov->continuous_parametric_process[i]->ident != AUTOREGRESSIVE_MODEL)) {
            hsmarkov->continuous_parametric_process[i]->restoration_weight = new Distribution(*weight);
          }
        }

        delete weight;

        if ((display) && (seq->characteristics[0])) {
          cout << "\n" << SEQ_label[SEQL_STATE_SEQUENCES_LIKELIHOOD] << ": " << seq->restoration_likelihood;

          for (i = 0;i < nb_variable;i++) {
            if (type[i] == REAL_VALUE) {
              break;
            }
          }
          if (i == nb_variable) {
            cout << " | " << hsmarkov->SemiMarkov::likelihood_computation(*seq);
          }
          cout << endl;
        }
      }

      else {
        if (hsmarkov->type == ORDINARY) {
          for (i = 0;i < hsmarkov->nb_state;i++) {
            if ((hsmarkov->sojourn_type[i] == SEMI_MARKOVIAN) && (hsmarkov->stype[i] == RECURRENT)) {
              hsmarkov->forward[i]->copy(*(hsmarkov->state_process->sojourn_time[i]));
              hsmarkov->forward[i]->computation(*(hsmarkov->state_process->sojourn_time[i]));
            }
          }
        }

        hsmarkov->semi_markov_data = new SemiMarkovData(*this , SEQUENCE_COPY ,
                                                        (hsmarkov->type == EQUILIBRIUM ? true : false));
        seq = hsmarkov->semi_markov_data;
        if (seq->type[0] == STATE) {
          seq->state_variable_init(INT_VALUE);
        }

        for (i = 0;i < hsmarkov->nb_output_process;i++) {
          if (((hsmarkov->discrete_parametric_process[i]) || (hsmarkov->continuous_parametric_process[i])) &&
              (seq->characteristics[i])) {
            delete seq->characteristics[i];
            seq->characteristics[i] = NULL;
          }
        }
      }

      for (i = 0;i < hsmarkov->nb_output_process;i++) {
        if (hsmarkov->categorical_process[i]) {
          for (j = 0;j < hsmarkov->nb_state;j++) {
            hsmarkov->categorical_process[i]->observation[j]->cumul_computation();

            hsmarkov->categorical_process[i]->observation[j]->max_computation();
//            hsmarkov->categorical_process[i]->observation[j]->mean_computation();
//            hsmarkov->categorical_process[i]->observation[j]->variance_computation();
          }
        }
      }

      // computation of the log-likelihood and the characteristic distributions of the model

      seq->likelihood = hsmarkov->likelihood_computation(*this , seq->posterior_probability);

      hsmarkov->component_computation();
      hsmarkov->characteristic_computation(*seq , counting_flag , I_DEFAULT , false);

      // computation of the mixtures of observation distributions (theoretical weights)

      weight = NULL;

      for (i = 0;i < hsmarkov->nb_output_process;i++) {
        if ((hsmarkov->categorical_process[i]) || (hsmarkov->discrete_parametric_process[i]) ||
            ((hsmarkov->continuous_parametric_process[i]) &&
             (hsmarkov->continuous_parametric_process[i]->ident != LINEAR_MODEL) &&
             (hsmarkov->continuous_parametric_process[i]->ident != AUTOREGRESSIVE_MODEL))) {
          switch (hsmarkov->type) {
          case ORDINARY :
            weight = hsmarkov->state_process->weight_computation();
            break;
          case EQUILIBRIUM :
            weight = new Distribution(hsmarkov->nb_state , hsmarkov->initial);
            break;
          }
          break;
        }
      }

      for (i = 0;i < hsmarkov->nb_output_process;i++) {
        if (hsmarkov->categorical_process[i]) {
          hsmarkov->categorical_process[i]->weight = new Distribution(*weight);
          hsmarkov->categorical_process[i]->mixture = hsmarkov->categorical_process[i]->mixture_computation(hsmarkov->categorical_process[i]->weight);
        }

        else if (hsmarkov->discrete_parametric_process[i]) {
          hsmarkov->discrete_parametric_process[i]->weight = new Distribution(*weight);
          hsmarkov->discrete_parametric_process[i]->mixture = hsmarkov->discrete_parametric_process[i]->mixture_computation(hsmarkov->discrete_parametric_process[i]->weight);
        }

        else if ((hsmarkov->continuous_parametric_process[i]) &&
                 (hsmarkov->continuous_parametric_process[i]->ident != LINEAR_MODEL) &&
                 (hsmarkov->continuous_parametric_process[i]->ident != AUTOREGRESSIVE_MODEL)) {
          hsmarkov->continuous_parametric_process[i]->weight = new Distribution(*weight);
        }
      }

      delete weight;

      // update of the sample sizes for the computation of the confidence intervals on the slopes and
      // the correlation coefficients (linear observation model) and the autoregressive coefficients (autoregressive models)

      for (i = 0;i < hsmarkov->nb_output_process;i++) {
        if ((hsmarkov->continuous_parametric_process[i]->ident == LINEAR_MODEL) ||
            (hsmarkov->continuous_parametric_process[i]->ident == AUTOREGRESSIVE_MODEL)) {
          for (j = 0;j < hsmarkov->nb_state;j++) {
            hsmarkov->continuous_parametric_process[i]->observation[j]->sample_size /= nb_state_sequence;
          }
        }
      }

      if ((display) && (state_sequence) && (seq->nb_sequence <= POSTERIOR_PROBABILITY_NB_SEQUENCE)) {
        cout << "\n" << SEQ_label[SEQL_POSTERIOR_STATE_SEQUENCE_PROBABILITY] << endl;
        for (i = 0;i < seq->nb_sequence;i++) {
          cout << SEQ_label[SEQL_SEQUENCE] << " " << seq->identifier[i] << ": "
               << seq->posterior_probability[i];

          if (hsmarkov->nb_component == hsmarkov->nb_state) {
            cout << " | " << SEQ_label[SEQL_STATE_BEGIN] << ": ";

            pstate = seq->int_sequence[i][0] + 1;
            if (seq->index_parameter) {
              for (j = 1;j < seq->length[i];j++) {
                if (*pstate != *(pstate - 1)) {
                  cout << seq->index_parameter[i][j] << ", ";
                }
                pstate++;
              }
            }

            else {
              for (j = 1;j < seq->length[i];j++) {
                if (*pstate != *(pstate - 1)) {
                  cout << j << ", ";
                }
                pstate++;
              }
            }
          }

          cout << endl;
        }
      }
    }
  }

  return hsmarkov;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Estimation of a hidden semi-Markov chain using the MCEM algorithm.
 *
 *  \param[in] error                 reference on a StatError object,
 *  \param[in] display               flag for displaying estimation intermediate results,
 *  \param[in] itype                 process type (ORDINARY/EQUILIBRIUM),
 *  \param[in] nb_state              number of states,
 *  \param[in] left_right            flag on the Markov chain structure,
 *  \param[in] occupancy_mean        mean state occupancy,
 *  \param[in] poisson_geometric     flag on the estimation of Poisson geometric state occupancy distributions,
 *  \param[in] common_dispersion     flag common dispersion parameter (continuous observation processes),
 *  \param[in] min_nb_state_sequence minimum number of generated sequences,
 *  \param[in] max_nb_state_sequence maximum number of generated sequences,
 *  \param[in] parameter             parameter for defining the number of generated sequences,
 *  \param[in] estimator             estimator type for the reestimation of the state occupancy distributions
 *                                   (complete or partial likelihood),
 *  \param[in] counting_flag         flag on the computation of the counting distributions,
 *  \param[in] state_sequence        flag on the computation of the restored state sequences,
 *  \param[in] nb_iter               number of iterations.
 *
 *  \return                          HiddenSemiMarkov object.
 */
/*--------------------------------------------------------------*/

HiddenSemiMarkov* MarkovianSequences::hidden_semi_markov_stochastic_estimation(StatError &error , bool display ,
                                                                               process_type itype , int nb_state ,
                                                                               bool left_right , double occupancy_mean ,
                                                                               bool poisson_geometric , bool common_dispersion ,
                                                                               int min_nb_state_sequence ,
                                                                               int max_nb_state_sequence , double parameter ,
                                                                               censoring_estimator estimator ,
                                                                               bool counting_flag , bool state_sequence ,
                                                                               int nb_iter) const

{
  bool status = true;
  int i;
  int nb_value[SEQUENCE_NB_VARIABLE];
  double proba , mean , variance;
  HiddenSemiMarkov *ihsmarkov , *hsmarkov;


  hsmarkov = NULL;
  error.init();

  if ((nb_state < 2) || (nb_state > NB_STATE)) {
    status = false;
    error.update(SEQ_error[SEQR_NB_STATE]);
  }
  if ((occupancy_mean != D_DEFAULT) && (occupancy_mean <= 1.)) {
    status = false;
    error.update(SEQ_error[SEQR_OCCUPANCY]);
  }

  for (i = 0;i < nb_variable;i++) {
    if ((type[i] != INT_VALUE) && (type[i] != REAL_VALUE) && (type[i] != STATE)) {
      status = false;
      ostringstream error_message , correction_message;
      error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                    << STAT_error[STATR_VARIABLE_TYPE];
      correction_message << STAT_variable_word[INT_VALUE] << " or "
                         << STAT_variable_word[REAL_VALUE];
      error.correction_update((error_message.str()).c_str() , (correction_message.str()).c_str());
    }
  }

  if (status) {
    for (i = 0;i < nb_variable;i++) {
      if (marginal_distribution[i]) {
        nb_value[i] = marginal_distribution[i]->nb_value;
      }
      else {
        nb_value[i] = I_DEFAULT;
      }
    }

    ihsmarkov = new HiddenSemiMarkov(itype , nb_state , nb_variable , nb_value);

    // initialization of the Markov chain parameters

    ihsmarkov->init(left_right , 0.);

    // initialization of the state occupancy distributions

    if (occupancy_mean == D_DEFAULT) {
      occupancy_mean = MAX(length_distribution->mean , OCCUPANCY_MEAN);
    }

    ihsmarkov->sojourn_type = new state_sojourn_type[nb_state];
    ihsmarkov->state_process->absorption = new double[nb_state];
    ihsmarkov->state_process->sojourn_time = new DiscreteParametric*[nb_state];
    ihsmarkov->forward = new Forward*[nb_state];

    for (i = 0;i < nb_state;i++) {
      if (ihsmarkov->stype[i] != ABSORBING) {
        ihsmarkov->sojourn_type[i] = SEMI_MARKOVIAN;
        ihsmarkov->state_process->absorption[i] = 0.;
        proba = 1. / occupancy_mean;
        ihsmarkov->state_process->sojourn_time[i] = new DiscreteParametric(NEGATIVE_BINOMIAL , 1 ,
                                                                           I_DEFAULT , 1. , proba ,
                                                                           OCCUPANCY_THRESHOLD);

        if (ihsmarkov->stype[i] == RECURRENT) {
          ihsmarkov->forward[i] = new Forward(*(ihsmarkov->state_process->sojourn_time[i]) ,
                                              ihsmarkov->state_process->sojourn_time[i]->alloc_nb_value);
        }
        else {
          ihsmarkov->forward[i] = NULL;
        }
      }

      else {
        ihsmarkov->sojourn_type[i] = MARKOVIAN;
        ihsmarkov->state_process->absorption[i] = 1.;
        ihsmarkov->state_process->sojourn_time[i] = NULL;
        ihsmarkov->forward[i] = NULL;
      }
    }

    // initialization of the observation distributions

    for (i = 0;i < ihsmarkov->nb_output_process;i++) {
      if (ihsmarkov->categorical_process[i]) {
        ihsmarkov->categorical_process[i]->init();
      }

      else if (ihsmarkov->discrete_parametric_process[i]) {
        ihsmarkov->discrete_parametric_process[i]->init();
      }

      else {
        mean = mean_computation(i);
        variance = variance_computation(i , mean);

        ihsmarkov->continuous_parametric_process[i]->init(GAUSSIAN , min_value[i] , max_value[i] ,
                                                          mean , variance);
      }
    }

    hsmarkov = hidden_semi_markov_stochastic_estimation(error , display , *ihsmarkov , poisson_geometric ,
                                                        common_dispersion , min_nb_state_sequence ,
                                                        max_nb_state_sequence , parameter , estimator ,
                                                        counting_flag , state_sequence , nb_iter);
    delete ihsmarkov;
  }

  return hsmarkov;
}


};  // namespace sequence_analysis
