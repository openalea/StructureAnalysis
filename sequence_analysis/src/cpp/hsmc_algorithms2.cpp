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



#include <math.h>

#include <string>
#include <sstream>

#include "stat_tool/stat_label.h"

#include "hidden_semi_markov.h"
#include "sequence_label.h"

using namespace std;
using namespace stat_tool;


namespace sequence_analysis {



/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the state sequence entropies using the forward-backward algorithm.
 *
 *  \param[in] seq reference on a SemiMarkovData object,
 */
/*--------------------------------------------------------------*/

void HiddenSemiMarkov::forward_backward(SemiMarkovData &seq) const

{
  bool posterior_state_probability_flag;
  int i , j , k , m , n;
  int **pioutput;
  double seq_likelihood , obs_product , residual , buff , sum , **observation ,
         *norm , *state_norm , **forward1 , **state_in , *transition_predicted ,
         *occupancy_predicted , **state_entropy , **predicted_entropy , **proutput;
  DiscreteParametric *occupancy;

# ifdef MESSAGE
  double entropy , **backward , **backward1 , *auxiliary , *occupancy_auxiliary ,
         **transition_entropy , **occupancy_entropy;
# endif


  // initializations

  seq.entropy = new double[seq.nb_sequence];
  seq.nb_state_sequence = new double[seq.nb_sequence];

  posterior_state_probability_flag = parallel_initial_state();
  if (posterior_state_probability_flag) {
    seq.posterior_state_probability = new double[seq.nb_sequence];
  }

  observation = new double*[seq.max_length];
  for (i = 0;i < seq.max_length;i++) {
    observation[i] = new double[nb_state];
  }

  norm = new double[seq.max_length];
  state_norm = new double[nb_state];

  forward1 = new double*[seq.max_length];
  for (i = 0;i < seq.max_length;i++) {
    forward1[i] = new double[nb_state];
  }

  state_in = new double*[seq.max_length - 1];
  for (i = 0;i < seq.max_length - 1;i++) {
    state_in[i] = new double[nb_state];
  }

  transition_predicted = new double[nb_state];
  occupancy_predicted = new double[seq.max_length + 1];

  state_entropy = new double*[seq.max_length];
  for (i = 0;i < seq.max_length;i++) {
    state_entropy[i] = new double[nb_state];
  }

  predicted_entropy = new double*[seq.max_length];
  for (i = 0;i < seq.max_length;i++) {
    predicted_entropy[i] = new double[nb_state];
  }

# ifdef MESSAGE
  backward = new double*[seq.max_length];
  for (i = 0;i < seq.max_length;i++) {
    backward[i] = new double[nb_state];
  }

  backward1 = new double*[seq.max_length];
  for (i = 0;i < seq.max_length;i++) {
    backward1[i] = new double[nb_state];
  }

  auxiliary = new double[nb_state];
  occupancy_auxiliary = new double[seq.max_length + 1];

  transition_entropy = new double*[nb_state];
  for (i = 0;i < nb_state;i++) {
    transition_entropy[i] = new double[nb_state];
  }

  occupancy_entropy = new double*[nb_state];
  for (i = 0;i < nb_state;i++) {
    switch (sojourn_type[i]) {
    case SEMI_MARKOVIAN :
      occupancy = state_process->sojourn_time[i];
      occupancy_entropy[i] = new double[MIN(seq.max_length , occupancy->nb_value)];
      break;
    case MARKOVIAN :
      occupancy_entropy[i] = NULL;
      break;
    }
  }
# endif

  pioutput = new int*[nb_output_process];
  proutput = new double*[nb_output_process];

  seq.sample_entropy = 0.;

  for (i = 0;i < seq.nb_sequence;i++) {
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

    // forward recurrence

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
                 (continuous_parametric_process[m]->ident == ZERO_INFLATED_GAMMA)) && (seq.min_value[m + 1] < seq.min_interval[m + 1] / 2)) {
              switch (seq.type[m + 1]) {
              case INT_VALUE :
                observation[j][k] *= continuous_parametric_process[m]->observation[k]->mass_computation(*pioutput[m] , *pioutput[m] + seq.min_interval[m + 1]);
                break;
              case REAL_VALUE :
                observation[j][k] *= continuous_parametric_process[m]->observation[k]->mass_computation(*proutput[m] , *proutput[m] + seq.min_interval[m + 1]);
                break;
              }
            }

            else if (continuous_parametric_process[m]->ident == LINEAR_MODEL) {
              switch (seq.type[m + 1]) {
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
                switch (seq.type[m + 1]) {
                case INT_VALUE :
                  residual = *pioutput[m] - continuous_parametric_process[m]->observation[k]->location;
                  break;
                case REAL_VALUE :
                  residual = *proutput[m] - continuous_parametric_process[m]->observation[k]->location;
                  break;
                }
              }

              else {
                switch (seq.type[m + 1]) {
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
              switch (seq.type[m + 1]) {
              case INT_VALUE :
                observation[j][k] *= continuous_parametric_process[m]->observation[k]->mass_computation(*pioutput[m] - seq.min_interval[m + 1] / 2 , *pioutput[m] + seq.min_interval[m + 1] / 2);
                break;
              case REAL_VALUE :
                observation[j][k] *= continuous_parametric_process[m]->observation[k]->mass_computation(*proutput[m] - seq.min_interval[m + 1] / 2 , *proutput[m] + seq.min_interval[m + 1] / 2);
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
            state_norm[k] += state_in[j - 1][k] - forward1[j - 1][k];
          }
          state_norm[k] *= observation[j][k];

          norm[j] += state_norm[k];
          break;
        }

        // case Markovian state

        case MARKOVIAN : {
          if (j == 0) {
            forward1[j][k] = initial[k];
            state_entropy[j][k] = 0.;
          }
          else {
            forward1[j][k] = state_in[j - 1][k];
            state_entropy[j][k] = predicted_entropy[j - 1][k];
          }
          forward1[j][k] *= observation[j][k];

          norm[j] += forward1[j][k];
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
            forward1[j][k] /= norm[j];
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
          forward1[j][k] = 0.;

          if (j < seq.length[i] - 1) {
            for (m = 1;m <= MIN(j + 1 , occupancy->nb_value - 1);m++) {
              obs_product *= observation[j - m + 1][k] / norm[j - m + 1];
              if (obs_product == 0.) {
                break;
              }

              if (m < j + 1) {
                occupancy_predicted[m] = obs_product * occupancy->mass[m] * state_in[j - m][k];
//                forward1[j][k] += obs_product * occupancy->mass[m] * state_in[j - m][k];
              }

              else {
                switch (type) {
                case ORDINARY :
                  occupancy_predicted[m] = obs_product * occupancy->mass[m] * initial[k];
//                  forward1[j][k] += obs_product * occupancy->mass[m] * initial[k];
                  break;
                case EQUILIBRIUM :
                  occupancy_predicted[m] = obs_product * forward[k]->mass[m] * initial[k];
//                  forward1[j][k] += obs_product * forward[k]->mass[m] * initial[k];
                  break;
                }
              }

              forward1[j][k] += occupancy_predicted[m];
            }
          }

          else {
            for (m = 1;m <= MIN(j + 1 , occupancy->nb_value - 1);m++) {
              obs_product *= observation[j - m + 1][k] / norm[j - m + 1];
              if (obs_product == 0.) {
                break;
              }

              if (m < j + 1) {
                occupancy_predicted[m] = obs_product * (1. - occupancy->cumul[m - 1]) * state_in[j - m][k];
//                forward1[j][k] += obs_product * (1. - occupancy->cumul[m - 1]) * state_in[j - m][k];
              }

              else {
                switch (type) {
                case ORDINARY :
                  occupancy_predicted[m] = obs_product * (1. - occupancy->cumul[m - 1]) * initial[k];
//                  forward1[j][k] += obs_product * (1. - occupancy->cumul[m - 1]) * initial[k];
                  break;
                case EQUILIBRIUM :
                  occupancy_predicted[m] = obs_product * (1. - forward[k]->cumul[m - 1]) * initial[k];
//                  forward1[j][k] += obs_product * (1. - forward[k]->cumul[m - 1]) * initial[k];
                  break;
                }
              }

              forward1[j][k] += occupancy_predicted[m];
            }
          }

          state_entropy[j][k] = 0.;

          if (forward1[j][k] > 0.) {
            for (n = 1;n < m;n++) {
              buff = occupancy_predicted[n] / forward1[j][k];
              if (buff > 0.) {
                if (n < j + 1) {
                  state_entropy[j][k] += buff * (predicted_entropy[j - n][k] - log(buff));
                }
                else {
                  state_entropy[j][k] -= buff * log(buff);
                }
              }
            }

            if (state_entropy[j][k] < 0.) {
              state_entropy[j][k] = 0.;
            }
          }
        }
      }

      if (j < seq.length[i] - 1) {
        for (k = 0;k < nb_state;k++) {
          state_in[j][k] = 0.;
          for (m = 0;m < nb_state;m++) {
            transition_predicted[m] = transition[m][k] * forward1[j][m];
            state_in[j][k] += transition_predicted[m];
//            state_in[j][k] += transition[m][k] * forward1[j][m];
          }

          predicted_entropy[j][k] = 0.;

          if (state_in[j][k] > 0.) {
            for (m = 0;m < nb_state;m++) {
              buff = transition_predicted[m] / state_in[j][k];
              if (buff > 0.) {
                predicted_entropy[j][k] += buff * (state_entropy[j][m] - log(buff));
              }
            }

            if (predicted_entropy[j][k] < 0.) {
              predicted_entropy[j][k] = 0.;
            }
          }
        }
      }

      for (k = 0;k < nb_output_process;k++) {
        switch (seq.type[k + 1]) {
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
      seq.entropy[i] = 0.;
      j = seq.length[i] - 1;
      for (k = 0;k < nb_state;k++) {
        if (forward1[j][k] > 0.) {
          seq.entropy[i] += forward1[j][k] * (state_entropy[j][k] - log(forward1[j][k]));
        }
      }
      seq.sample_entropy += seq.entropy[i];

/*      for (j = 0;j < nb_state;j++) {
        if (sojourn_type[j] == SEMI_MARKOVIAN) {
          for (k = 0;k < seq.length[i];k++) {
            state_entropy[k][j] = 0.;
          }
        }
        } */

      // backward recurrence

#     ifdef MESSAGE
      entropy = 0.;

      for (j = 0;j < nb_output_process;j++) {
        switch (seq.type[j + 1]) {
        case INT_VALUE :
          pioutput[j]--;
          break;
        case REAL_VALUE :
          proutput[j]--;
          break;
        }
      }

      for (j = 0;j < nb_state;j++) {
        for (k = 0;k < nb_state;k++) {
          transition_entropy[j][k] = 0.;
        }
      }

      for (j = 0;j < nb_state;j++) {
        if (sojourn_type[j] == SEMI_MARKOVIAN) {
          occupancy = state_process->sojourn_time[j];
          for (k = occupancy->offset;k < MIN(seq.length[i] , occupancy->nb_value);k++) {
            occupancy_entropy[j][k] = 0.;
          }
        }
      }

      j = seq.length[i] - 1;
      for (k = 0;k < nb_state;k++) {
        backward[j][k] = forward1[j][k];
        backward1[j][k] = backward[j][k];

        if (backward[j][k] > 0.) {
          for (m = 0;m < nb_output_process;m++) {
            if (categorical_process[m]) {
              if (categorical_process[m]->observation[k]->mass[*pioutput[m]] > 0.) {
                entropy -= backward[j][k] * log(categorical_process[m]->observation[k]->mass[*pioutput[m]]);
              }
            }

            else if (discrete_parametric_process[m]) {
              if (discrete_parametric_process[m]->observation[k]->mass[*pioutput[m]] > 0.) {
                entropy -= backward[j][k] * log(discrete_parametric_process[m]->observation[k]->mass[*pioutput[m]]);
              }
            }

            else {
              if (((continuous_parametric_process[m]->ident == GAMMA) ||
                   (continuous_parametric_process[m]->ident == ZERO_INFLATED_GAMMA)) && (seq.min_value[m + 1] < seq.min_interval[m + 1] / 2)) {
                switch (seq.type[m + 1]) {
                case INT_VALUE :
                  entropy -= backward[j][k] * log(continuous_parametric_process[m]->observation[k]->mass_computation(*pioutput[m] , *pioutput[m] + seq.min_interval[m + 1]));
                  break;
                case REAL_VALUE :
                  entropy -= backward[j][k] * log(continuous_parametric_process[m]->observation[k]->mass_computation(*proutput[m] , *proutput[m] + seq.min_interval[m + 1]));
                  break;
                }
              }

              else if (continuous_parametric_process[m]->ident == LINEAR_MODEL) {
                switch (seq.type[m + 1]) {
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

                entropy -= backward[j][k] * log(continuous_parametric_process[m]->observation[k]->mass_computation(residual , residual));
              }

              else if (continuous_parametric_process[m]->ident == AUTOREGRESSIVE_MODEL) {
                if (j == 0) {
                  switch (seq.type[m + 1]) {
                  case INT_VALUE :
                    residual = *pioutput[m] - continuous_parametric_process[m]->observation[k]->location;
                    break;
                  case REAL_VALUE :
                    residual = *proutput[m] - continuous_parametric_process[m]->observation[k]->location;
                    break;
                  }
                }

                else {
                  switch (seq.type[m + 1]) {
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

                entropy -= backward[j][k] * log(continuous_parametric_process[m]->observation[k]->mass_computation(residual , residual));
              }

              else {
                switch (seq.type[m + 1]) {
                case INT_VALUE :
                  entropy -= backward[j][k] * log(continuous_parametric_process[m]->observation[k]->mass_computation(*pioutput[m] - seq.min_interval[m + 1] / 2 , *pioutput[m] + seq.min_interval[m + 1] / 2));
                  break;
                case REAL_VALUE :
                  entropy -= backward[j][k] * log(continuous_parametric_process[m]->observation[k]->mass_computation(*proutput[m] - seq.min_interval[m + 1] / 2 , *proutput[m] + seq.min_interval[m + 1] / 2));
                  break;
                }
              }
            }
          }
        }
      }

      for (j = seq.length[i] - 2;j >= 0;j--) {
        for (k = 0;k < nb_output_process;k++) {
          switch (seq.type[k + 1]) {
          case INT_VALUE :
            pioutput[k]--;
            break;
          case REAL_VALUE :
            proutput[k]--;
            break;
          }
        }

        for (k = 0;k < nb_state;k++) {
          auxiliary[k] = 0.;

          switch (sojourn_type[k]) {

          // case semi-Markovian state

          case SEMI_MARKOVIAN : {
            occupancy = state_process->sojourn_time[k];
            obs_product = 1.;

            for (m = 1;m < MIN(seq.length[i] - j , occupancy->nb_value);m++) {
              obs_product *= observation[j + m][k] / norm[j + m];
              if (obs_product == 0.) {
                break;
              }

              occupancy_auxiliary[m] = 0.;

              if (backward1[j + m][k] > 0.) {
//              if (forward1[j + m][k] > 0.) {
                if (m < seq.length[i] - j - 1) {
                  buff = backward1[j + m][k] * obs_product * occupancy->mass[m] /
                         forward1[j + m][k];
                  occupancy_auxiliary[m] = buff * state_in[j][k];
                  occupancy_entropy[k][m] += occupancy_auxiliary[m];

/*                  if (occupancy->mass[m] > 0.) {
                    entropy -= occupancy_auxiliary[m] * log(occupancy->mass[m]);
                  } */
                }

                else {
                  buff = obs_product * (1. - occupancy->cumul[m - 1]);
                  occupancy_auxiliary[m] = buff * state_in[j][k];
                  if (occupancy->cumul[m - 1] < 1.) {
                    entropy -= occupancy_auxiliary[m] * log(1. - occupancy->cumul[m - 1]);
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
//            if (forward1[j + 1][k] > 0.) {
              auxiliary[k] = backward1[j + 1][k] / state_in[j][k];

/*              auxiliary[k] = backward1[j + 1][k] * observation[j + 1][k] /
                             (forward1[j + 1][k] * norm[j + 1]); */
            }
            break;
          }
          }
        }

        for (k = 0;k < nb_state;k++) {
          backward1[j][k] = 0.;

          for (m = 0;m < nb_state;m++) {
            buff = auxiliary[m] * transition[k][m] * forward1[j][k];
            backward1[j][k] += buff;
            transition_entropy[k][m] += buff;

/*            if (transition[k][m] > 0.) {
              entropy -= buff * log(transition[k][m]);
            } */
          }

          switch (sojourn_type[k]) {

          // case semi-Markovian state

          case SEMI_MARKOVIAN : {
            backward[j][k] = backward[j + 1][k] + backward1[j][k] - auxiliary[k] * state_in[j][k];
            if (backward[j][k] < 0.) {
              backward[j][k] = 0.;
            }
            if (backward[j][k] > 1.) {
              backward[j][k] = 1.;
            }
            break;
          }

          // case Markovian state

          case MARKOVIAN : {
            backward[j][k] = backward1[j][k];
            break;
          }
          }

          if (backward[j][k] > 0.) {
            for (m = 0;m < nb_output_process;m++) {
              if (categorical_process[m]) {
                if (categorical_process[m]->observation[k]->mass[*pioutput[m]] > 0.) {
                  entropy -= backward[j][k] * log(categorical_process[m]->observation[k]->mass[*pioutput[m]]);
                }
              }

              else if (discrete_parametric_process[m]) {
                if (discrete_parametric_process[m]->observation[k]->mass[*pioutput[m]] > 0.) {
                  entropy -= backward[j][k] * log(discrete_parametric_process[m]->observation[k]->mass[*pioutput[m]]);
                }
              }

              else {
                if (((continuous_parametric_process[m]->ident == GAMMA) ||
                     (continuous_parametric_process[m]->ident == ZERO_INFLATED_GAMMA)) && (seq.min_value[m + 1] < seq.min_interval[m + 1] / 2)) {
                  switch (seq.type[m + 1]) {
                  case INT_VALUE :
                    entropy -= backward[j][k] * log(continuous_parametric_process[m]->observation[k]->mass_computation(*pioutput[m] , *pioutput[m] + seq.min_interval[m + 1]));
                    break;
                  case REAL_VALUE :
                    entropy -= backward[j][k] * log(continuous_parametric_process[m]->observation[k]->mass_computation(*proutput[m] , *proutput[m] + seq.min_interval[m + 1]));
                    break;
                  }
                }

                else if (continuous_parametric_process[m]->ident == LINEAR_MODEL) {
                  switch (seq.type[m + 1]) {
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

                  entropy -= backward[j][k] * log(continuous_parametric_process[m]->observation[k]->mass_computation(residual , residual));
                }

                else if (continuous_parametric_process[m]->ident == AUTOREGRESSIVE_MODEL) {
                  if (j == 0) {
                    switch (seq.type[m + 1]) {
                    case INT_VALUE :
                      residual = *pioutput[m] - continuous_parametric_process[m]->observation[k]->location;
                      break;
                    case REAL_VALUE :
                      residual = *proutput[m] - continuous_parametric_process[m]->observation[k]->location;
                      break;
                    }
                  }

                  else {
                    switch (seq.type[m + 1]) {
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

                  entropy -= backward[j][k] * log(continuous_parametric_process[m]->observation[k]->mass_computation(residual , residual));
                }

                else {
                  switch (seq.type[m + 1]) {
                  case INT_VALUE :
                    entropy -= backward[j][k] * log(continuous_parametric_process[m]->observation[k]->mass_computation(*pioutput[m] - seq.min_interval[m + 1] / 2 , *pioutput[m] + seq.min_interval[m + 1] / 2));
                    break;
                  case REAL_VALUE :
                    entropy -= backward[j][k] * log(continuous_parametric_process[m]->observation[k]->mass_computation(*proutput[m] - seq.min_interval[m + 1] / 2 , *proutput[m] + seq.min_interval[m + 1] / 2));
                    break;
                  }
                }
              }
            }
          }
        }
      }

      if (posterior_state_probability_flag) {
        seq.posterior_state_probability[i] = 0.;
        for (j = 0;j < nb_state;j++) {
          if (backward[0][j] > seq.posterior_state_probability[i]) {
            seq.posterior_state_probability[i] = backward[0][j];
          }
        }
      }

      for (j = 0;j < nb_state;j++) {
        if (initial[j] > 0.) {
          entropy -= backward[0][j] * log(initial[j]);
        }
      }

      for (j = 0;j < nb_state;j++) {
        for (k = 0;k < nb_state;k++) {
          if (transition[j][k] > 0.) {
            entropy -= transition_entropy[j][k] * log(transition[j][k]);
          }
        }
      }

      for (j = 0;j < nb_state;j++) {
        if (sojourn_type[j] == SEMI_MARKOVIAN) {
          occupancy = state_process->sojourn_time[j];

          if (initial[j] > 0.) {
            obs_product = 1.;

#           ifdef DEBUG
            backward0[j] = 0.;
#           endif

            for (k = 1;k < MIN(seq.length[i] + 1 , occupancy->nb_value);k++) {
              obs_product *= observation[k - 1][j] / norm[k - 1];
              if (obs_product == 0.) {
                break;
              }

              occupancy_auxiliary[k] = 0.;

              if (backward1[k - 1][j] > 0.) {
//              if (forward1[k - 1][j] > 0.) {
                if (k < seq.length[i]) {
                  switch (type) {

                  case ORDINARY : {
                    occupancy_auxiliary[k] = backward1[k - 1][j] * obs_product * occupancy->mass[k] *
                                             initial[j] / forward1[k - 1][j];
                    occupancy_entropy[j][k] += occupancy_auxiliary[k];

/*                    if (occupancy->mass[k] > 0.) {
                      entropy -= occupancy_auxiliary[k] * log(occupancy->mass[k]);
                    } */
                    break;
                  }

                  case EQUILIBRIUM : {
                    occupancy_auxiliary[k] = backward1[k - 1][j] * obs_product * forward[j]->mass[k] *
                                             initial[j] / forward1[k - 1][j];
                    if (forward[j]->mass[k] > 0.) {
                      entropy -= occupancy_auxiliary[k] * log(forward[j]->mass[k]);
                    }
                    break;
                  }
                  }
                }

                else {
                  switch (type) {

                  case ORDINARY : {
                    occupancy_auxiliary[k] = obs_product * (1. - occupancy->cumul[k - 1]) * initial[j];
                    if (occupancy->cumul[k - 1] < 1.) {
                      entropy -= occupancy_auxiliary[k] * log(1. - occupancy->cumul[k - 1]);
                    }
                    break;
                  }

                  case EQUILIBRIUM : {
                    occupancy_auxiliary[k] = obs_product * (1. - forward[j]->cumul[k - 1]) * initial[j];
                    if (forward[j]->cumul[k - 1] < 1.) {
                      entropy -= occupancy_auxiliary[k] * log(1. - forward[j]->cumul[k - 1]);
                    }
                    break;
                  }
                  }
                }

#               ifdef DEBUG
                backward0[j] += occupancy_auxiliary[k];
#               endif

              }
            }

#           ifdef DEBUG
            cout << j << " " << backward[0][j] << " " << backward0[j] << endl;
#           endif
          }

          for (k = occupancy->offset;k < MIN(seq.length[i] , occupancy->nb_value);k++) {
            if (occupancy->mass[k] > 0.) {
              entropy -= occupancy_entropy[j][k] * log(occupancy->mass[k]);
            }
          }
        }
      }

      entropy += seq_likelihood;

      if ((entropy < seq.entropy[i] - DOUBLE_ERROR) || (entropy > seq.entropy[i] + DOUBLE_ERROR)) {
        cout << "\nERROR: " << i << " " << seq.entropy[i] << " " << entropy << endl;
      }
#     endif

      // computation of the number of state sequences

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

      // forward recurrence

      for (j = 0;j < seq.length[i];j++) {
        for (k = 0;k < nb_state;k++) {

          // computation of the indicator functions of the observation probabilities

          if (observation[j][k] > 0.) {
            observation[j][k] = 1.;
          }

          forward1[j][k] = 0.;

          switch (sojourn_type[k]) {

          // case semi-Markovian state

          case SEMI_MARKOVIAN : {
            occupancy = state_process->sojourn_time[k];

            if (j < seq.length[i] - 1) {
              for (m = 1;m <= MIN(j + 1 , occupancy->nb_value - 1);m++) {
                if (observation[j - m + 1][k] == 0.) {
                  break;
                }

                if (m < j + 1) {
                  if (occupancy->mass[m] > 0.) {
                    forward1[j][k] += state_in[j - m][k];
                  }
                }

                else {
                  if (initial[k] > 0.) {
                    switch (type) {

                    case ORDINARY : {
                      if (occupancy->mass[m] > 0.) {
                        forward1[j][k]++;
                      }
                      break;
                    }

                    case EQUILIBRIUM : {
                      if (forward[k]->mass[m] > 0.) {
                        forward1[j][k]++;
                      }
                      break;
                    }
                    }
                  }
                }
              }
            }

            else {
              for (m = 1;m <= MIN(j + 1 , occupancy->nb_value - 1);m++) {
                if (observation[j - m + 1][k] == 0.) {
                  break;
                }

                if (m < j + 1) {
                  if (1. - occupancy->cumul[m - 1] > 0.) {
                    forward1[j][k] += state_in[j - m][k];
                  }
                }

                else {
                  if (initial[k] > 0.) {
                    switch (type) {

                    case ORDINARY : {
                      if (1. - occupancy->cumul[m - 1] > 0.) {
                        forward1[j][k]++;
                      }
                      break;
                    }

                    case EQUILIBRIUM : {
                      if (1. - forward[k]->cumul[m - 1] > 0.) {
                        forward1[j][k]++;
                      }
                      break;
                    }
                    }
                  }
                }
              }
            }
            break;
          }

          // case Markovian state

          case MARKOVIAN : {
            if (observation[j][k] == 1.) {
              if (j == 0) {
                if (initial[k] > 0.) {
                  forward1[j][k] = 1.;
                }
              }
              else {
                forward1[j][k] = state_in[j - 1][k];
              }
            }
            break;
          }
          }
        }

        if (j < seq.length[i] - 1) {
          for (k = 0;k < nb_state;k++) {
            state_in[j][k] = 0.;
            for (m = 0;m < nb_state;m++) {
              if (transition[m][k] > 0.) {
                state_in[j][k] += forward1[j][m];
              }
            }
          }
        }

        for (k = 0;k < nb_output_process;k++) {
          switch (seq.type[k + 1]) {
          case INT_VALUE :
            pioutput[k]++;
            break;
          case REAL_VALUE :
            proutput[k]++;
            break;
          }
        }
      }

      seq.nb_state_sequence[i] = 0.;
      j = seq.length[i] - 1;
      for (k = 0;k < nb_state;k++) {
        seq.nb_state_sequence[i] += forward1[j][k];
      }
    }
  }

  for (i = 0;i < seq.max_length;i++) {
    delete [] observation[i];
  }
  delete [] observation;

  delete [] norm;
  delete [] state_norm;

  for (i = 0;i < seq.max_length;i++) {
    delete [] forward1[i];
  }
  delete [] forward1;

  for (i = 0;i < seq.max_length - 1;i++) {
    delete [] state_in[i];
  }
  delete [] state_in;

  delete [] transition_predicted;
  delete [] occupancy_predicted;

  for (i = 0;i < seq.max_length;i++) {
    delete [] state_entropy[i];
  }
  delete [] state_entropy;

  for (i = 0;i < seq.max_length;i++) {
    delete [] predicted_entropy[i];
  }
  delete [] predicted_entropy;

# ifdef MESSAGE
  for (i = 0;i < seq.max_length;i++) {
    delete [] backward[i];
  }
  delete [] backward;

  for (i = 0;i < seq.max_length;i++) {
    delete [] backward1[i];
  }
  delete [] backward1;

  delete [] auxiliary;
  delete [] occupancy_auxiliary;

  for (i = 0;i < nb_state;i++) {
    delete [] transition_entropy[i];
  }
  delete [] transition_entropy;

  for (i = 0;i < nb_state;i++) {
    delete [] occupancy_entropy[i];
  }
  delete [] occupancy_entropy;
# endif

  delete [] pioutput;
  delete [] proutput;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of state and entropy profiles using the forward-backward algorithm.
 *
 *  \param[in] seq                  reference on a MarkovianSequences object,
 *  \param[in] index                sequence index,
 *  \param[in] os                   stream,
 *  \param[in] plot_set             pointer on a MultiPlotSet object,
 *  \param[in] output               output type,
 *  \param[in] format               output format (ASCII/SPREADSHEET/GNUPLOT/PLOT),
 *  \param[in] max_marginal_entropy reference on the maximum marginal entropy,
 *  \param[in] entropy1             reference on the entropy (for the plots).
 *
 *  \return                         log-likelihood for the observed sequence.
 */
/*--------------------------------------------------------------*/

double HiddenSemiMarkov::forward_backward(MarkovianSequences &seq , int index , ostream *os ,
                                          MultiPlotSet *plot_set , state_profile output ,
                                          output_format format , double &max_marginal_entropy ,
                                          double &entropy1) const

{
  int i , j , k , m;
  int *pstate , **pioutput;
  double seq_likelihood , state_seq_likelihood , obs_product , residual , entropy2 , buff , sum ,
         backward_max , **observation , *norm , *state_norm , **forward1 , **state_in ,
         **backward , **backward1 , *auxiliary , *occupancy_auxiliary , **backward_output ,
         *transition_predicted , *occupancy_predicted , **state_entropy , **predicted_entropy ,
         **transition_entropy , **occupancy_entropy , *partial_entropy , *conditional_entropy ,
         *marginal_entropy , **proutput;
  DiscreteParametric *occupancy;


  // initializations

  observation = new double*[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    observation[i] = new double[nb_state];
  }

  norm = new double[seq.length[index]];
  state_norm = new double[nb_state];

  forward1 = new double*[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    forward1[i] = new double[nb_state];
  }

  state_in = new double*[seq.length[index] - 1];
  for (i = 0;i < seq.length[index] - 1;i++) {
    state_in[i] = new double[nb_state];
  }

  backward = new double*[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    backward[i] = new double[nb_state];
  }

  backward1 = new double*[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    backward1[i] = new double[nb_state];
  }

  auxiliary = new double[nb_state];
  occupancy_auxiliary = new double[seq.length[index] + 1];

  if (output == SSTATE) {
    backward_output = backward;
  }
  else {
    backward_output = new double*[seq.length[index]];
    for (i = 0;i < seq.length[index];i++) {
      backward_output[i] = new double[nb_state];
    }
  }

  transition_predicted = new double[nb_state];
  occupancy_predicted = new double[seq.length[index] + 1];

  state_entropy = new double*[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    state_entropy[i] = new double[nb_state];
  }

  predicted_entropy = new double*[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    predicted_entropy[i] = new double[nb_state];
  }

  transition_entropy = new double*[nb_state];
  for (i = 0;i < nb_state;i++) {
    transition_entropy[i] = new double[nb_state];
  }

  occupancy_entropy = new double*[nb_state];
  for (i = 0;i < nb_state;i++) {
    switch (sojourn_type[i]) {
    case SEMI_MARKOVIAN :
      occupancy = state_process->sojourn_time[i];
      occupancy_entropy[i] = new double[MIN(seq.length[index] , occupancy->nb_value)];
      break;
    case MARKOVIAN :
      occupancy_entropy[i] = NULL;
      break;
    }
  }

  partial_entropy = new double[seq.length[index]];
  conditional_entropy = new double[seq.length[index]];
  marginal_entropy = new double[seq.length[index]];

# ifdef DEBUG
  double *backward0;

  backward0 = new double[nb_state];
# endif

  pioutput = new int*[nb_output_process];
  proutput = new double*[nb_output_process];

  for (i = 0;i < nb_output_process;i++) {
    switch (seq.type[i + 1]) {
    case INT_VALUE :
      pioutput[i] = seq.int_sequence[index][i + 1];
      break;
    case REAL_VALUE :
      proutput[i] = seq.real_sequence[index][i + 1];
      break;
    }
  }

  // forward recurrence

  seq_likelihood = 0.;
  for (i = 0;i < seq.length[index];i++) {
    norm[i] = 0.;

    for (j = 0;j < nb_state;j++) {

      // computation of the observation probabilities

      observation[i][j] = 1.;
      for (k = 0;k < nb_output_process;k++) {
        if (categorical_process[k]) {
          observation[i][j] *= categorical_process[k]->observation[j]->mass[*pioutput[k]];
        }

        else if (discrete_parametric_process[k]) {
          observation[i][j] *= discrete_parametric_process[k]->observation[j]->mass[*pioutput[k]];
        }

        else {
          if (((continuous_parametric_process[k]->ident == GAMMA) ||
               (continuous_parametric_process[k]->ident == ZERO_INFLATED_GAMMA)) && (seq.min_value[k + 1] < seq.min_interval[k + 1] / 2)) {
            switch (seq.type[k + 1]) {
            case INT_VALUE :
              observation[i][j] *= continuous_parametric_process[k]->observation[j]->mass_computation(*pioutput[k] , *pioutput[k] + seq.min_interval[k + 1]);
              break;
            case REAL_VALUE :
              observation[i][j] *= continuous_parametric_process[k]->observation[j]->mass_computation(*proutput[k] , *proutput[k] + seq.min_interval[k + 1]);
              break;
            }
          }

          else if (continuous_parametric_process[k]->ident == LINEAR_MODEL) {
            switch (seq.type[k + 1]) {
            case INT_VALUE :
              residual = *pioutput[k] - (continuous_parametric_process[k]->observation[j]->intercept +
                          continuous_parametric_process[k]->observation[j]->slope *
                          (seq.index_param_type == IMPLICIT_TYPE ? i : seq.index_parameter[index][i]));
              break;
            case REAL_VALUE :
              residual = *proutput[k] - (continuous_parametric_process[k]->observation[j]->intercept +
                          continuous_parametric_process[k]->observation[j]->slope *
                          (seq.index_param_type == IMPLICIT_TYPE ? i : seq.index_parameter[index][i]));
              break;
            }

            observation[i][j] *= continuous_parametric_process[k]->observation[j]->mass_computation(residual , residual);
          }

          else if (continuous_parametric_process[k]->ident == AUTOREGRESSIVE_MODEL) {
            if (i == 0) {
              switch (seq.type[k + 1]) {
              case INT_VALUE :
                residual = *pioutput[k] - continuous_parametric_process[k]->observation[j]->location;
                break;
              case REAL_VALUE :
                residual = *proutput[k] - continuous_parametric_process[k]->observation[j]->location;
                break;
              }
            }

            else {
              switch (seq.type[k + 1]) {
              case INT_VALUE :
                residual = *pioutput[k] - (continuous_parametric_process[k]->observation[j]->location +
                            continuous_parametric_process[k]->observation[j]->autoregressive_coeff *
                            (*(pioutput[k] - 1) - continuous_parametric_process[k]->observation[j]->location));
                break;
              case REAL_VALUE :
                residual = *proutput[k] - (continuous_parametric_process[k]->observation[j]->location +
                            continuous_parametric_process[k]->observation[j]->autoregressive_coeff *
                            (*(proutput[k] - 1) - continuous_parametric_process[k]->observation[j]->location));
                break;
              }
            }

            observation[i][j] *= continuous_parametric_process[k]->observation[j]->mass_computation(residual , residual);
          }

          else {
            switch (seq.type[k + 1]) {
            case INT_VALUE :
              observation[i][j] *= continuous_parametric_process[k]->observation[j]->mass_computation(*pioutput[k] - seq.min_interval[k + 1] / 2 , *pioutput[k] + seq.min_interval[k + 1] / 2);
              break;
            case REAL_VALUE :
              observation[i][j] *= continuous_parametric_process[k]->observation[j]->mass_computation(*proutput[k] - seq.min_interval[k + 1] / 2 , *proutput[k] + seq.min_interval[k + 1] / 2);
              break;
            }
          }
        }
      }

      switch (sojourn_type[j]) {

      // case semi-Markovian state 

      case SEMI_MARKOVIAN : {
        if (i == 0) {
          state_norm[j] = initial[j];
        }
        else {
          state_norm[j] += state_in[i - 1][j] - forward1[i - 1][j];
        }
        state_norm[j] *= observation[i][j];

        norm[i] += state_norm[j];
        break;
      }

      // case Markovian state

      case MARKOVIAN : {
        if (i == 0) {
          forward1[i][j] = initial[j];
          state_entropy[i][j] = 0.;
        }
        else {
          forward1[i][j] = state_in[i - 1][j];
          state_entropy[i][j] = predicted_entropy[i - 1][j];
        }
        forward1[i][j] *= observation[i][j];

        norm[i] += forward1[i][j];
        break;
      }
      }
    }

    if (norm[i] > 0.) {
      for (j = 0;j < nb_state;j++) {
        switch (sojourn_type[j]) {
        case SEMI_MARKOVIAN :
          state_norm[j] /= norm[i];
          break;
        case MARKOVIAN :
          forward1[i][j] /= norm[i];
          break;
        }
      }

      seq_likelihood += log(norm[i]);
    }

    else {
      seq_likelihood = D_INF;
      break;
    }

    for (j = 0;j < nb_state;j++) {

      // case semi-Markovian state

      if (sojourn_type[j] == SEMI_MARKOVIAN) {
        occupancy = state_process->sojourn_time[j];
        obs_product = 1.;
        forward1[i][j] = 0.;

        if (i < seq.length[index] - 1) {
          for (k = 1;k <= MIN(i + 1 , occupancy->nb_value - 1);k++) {
            obs_product *= observation[i - k + 1][j] / norm[i - k + 1];
            if (obs_product == 0.) {
              break;
            }

            if (k < i + 1) {
              occupancy_predicted[k] = obs_product * occupancy->mass[k] * state_in[i - k][j];
//              forward1[i][j] += obs_product * occupancy->mass[k] * state_in[i - k][j];
            }

            else {
              switch (type) {
              case ORDINARY :
                occupancy_predicted[k] = obs_product * occupancy->mass[k] * initial[j];
//                forward1[i][j] += obs_product * occupancy->mass[k] * initial[j];
                break;
              case EQUILIBRIUM :
                occupancy_predicted[k] = obs_product * forward[j]->mass[k] * initial[j];
//                forward1[i][j] += obs_product * forward[j]->mass[k] * initial[j];
                break;
              }
            }

            forward1[i][j] += occupancy_predicted[k];
          }
        }

        else {
          for (k = 1;k <= MIN(i + 1 , occupancy->nb_value - 1);k++) {
            obs_product *= observation[i - k + 1][j] / norm[i - k + 1];
            if (obs_product == 0.) {
              break;
            }

            if (k < i + 1) {
              occupancy_predicted[k] = obs_product * (1. - occupancy->cumul[k - 1]) * state_in[i - k][j];
//              forward1[i][j] += obs_product * (1. - occupancy->cumul[k - 1]) * state_in[i - k][j];
            }

            else {
              switch (type) {
              case ORDINARY :
                occupancy_predicted[k] = obs_product * (1. - occupancy->cumul[k - 1]) * initial[j];
//                forward1[i][j] += obs_product * (1. - occupancy->cumul[k - 1]) * initial[j];
                break;
              case EQUILIBRIUM :
                occupancy_predicted[k] = obs_product * (1. - forward[j]->cumul[k - 1]) * initial[j];
//                forward1[i][j] += obs_product * (1. - forward[j]->cumul[k - 1]) * initial[j];
                break;
              }
            }

            forward1[i][j] += occupancy_predicted[k];
          }
        }

        state_entropy[i][j] = 0.;

        if (forward1[i][j] > 0.) {
          for (m = 1;m < k;m++) {
            buff = occupancy_predicted[m] / forward1[i][j];
            if (buff > 0.) {
              if (m < i + 1) {
                state_entropy[i][j] += buff * (predicted_entropy[i - m][j] - log(buff));
              }
              else {
                state_entropy[i][j] -= buff * log(buff);
              }
            }
          }

          if (state_entropy[i][j] < 0.) {
            state_entropy[i][j] = 0.;
          }
        }
      }
    }

    if (i < seq.length[index] - 1) {
      for (j = 0;j < nb_state;j++) {
        state_in[i][j] = 0.;
        for (k = 0;k < nb_state;k++) {
          transition_predicted[k] = transition[k][j] * forward1[i][k];
          state_in[i][j] += transition_predicted[k];
//          state_in[i][j] += transition[k][j] * forward1[i][k];
        }

        predicted_entropy[i][j] = 0.;

        if (state_in[i][j] > 0.) {
          for (k = 0;k < nb_state;k++) {
            buff = transition_predicted[k] / state_in[i][j];
            if (buff > 0.) {
              predicted_entropy[i][j] += buff * (state_entropy[i][k] - log(buff));
            }
          }

          if (predicted_entropy[i][j] < 0.) {
            predicted_entropy[i][j] = 0.;
          }
        }
      }
    }

    for (j = 0;j < nb_output_process;j++) {
      switch (seq.type[j + 1]) {
      case INT_VALUE :
        pioutput[j]++;
        break;
      case REAL_VALUE :
        proutput[j]++;
        break;
      }
    }
  }

  if (seq_likelihood != D_INF) {
    entropy1 = 0.;
    i = seq.length[index] - 1;
    for (j = 0;j < nb_state;j++) {
      if (forward1[i][j] > 0.) {
        entropy1 += forward1[i][j] * (state_entropy[i][j] - log(forward1[i][j]));
      }
    }

    for (i = 0;i < nb_state;i++) {
      if (sojourn_type[i] == SEMI_MARKOVIAN) {
        for (j = 0;j < seq.length[index];j++) {
          state_entropy[j][i] = 0.;
        }
      }
    }

    // backward recurrence

    for (i = 0;i < nb_output_process;i++) {
      switch (seq.type[i + 1]) {
      case INT_VALUE :
        pioutput[i]--;
        break;
      case REAL_VALUE :
        proutput[i]--;
        break;
      }
    }

    for (i = 0;i < nb_state;i++) {
      for (j = 0;j < nb_state;j++) {
        transition_entropy[i][j] = 0.;
      }
    }

    for (i = 0;i < nb_state;i++) {
      if (sojourn_type[i] == SEMI_MARKOVIAN) {
        occupancy = state_process->sojourn_time[i];
        for (j = occupancy->offset;j < MIN(seq.length[index] , occupancy->nb_value);j++) {
          occupancy_entropy[i][j] = 0.;
        }
      }
    }

    entropy2 = 0.;

    i = seq.length[index] - 1;
    for (j = 0;j < nb_state;j++) {
      backward[i][j] = forward1[i][j];
      backward1[i][j] = backward[i][j];

      if (output == OUT_STATE) {
        backward_output[i][j] = backward[i][j];
      }

      if (backward[i][j] > 0.) {
        for (k = 0;k < nb_output_process;k++) {
          if (categorical_process[k]) {
            if (categorical_process[k]->observation[j]->mass[*pioutput[k]] > 0.) {
              entropy2 -= backward[i][j] * log(categorical_process[k]->observation[j]->mass[*pioutput[k]]);
            }
          }

          else if (discrete_parametric_process[k]) {
            if (discrete_parametric_process[k]->observation[j]->mass[*pioutput[k]] > 0.) {
              entropy2 -= backward[i][j] * log(discrete_parametric_process[k]->observation[j]->mass[*pioutput[k]]);
            }
          }

          else {
            if (((continuous_parametric_process[k]->ident == GAMMA) ||
                 (continuous_parametric_process[k]->ident == ZERO_INFLATED_GAMMA)) && (seq.min_value[k + 1] < seq.min_interval[k + 1] / 2)) {
              switch (seq.type[k + 1]) {
              case INT_VALUE :
                entropy2 -= backward[i][j] * log(continuous_parametric_process[k]->observation[j]->mass_computation(*pioutput[k] , *pioutput[k] + seq.min_interval[k + 1]));
                break;
              case REAL_VALUE :
                entropy2 -= backward[i][j] * log(continuous_parametric_process[k]->observation[j]->mass_computation(*proutput[k] , *proutput[k] + seq.min_interval[k + 1]));
                break;
              }
            }

            else if (continuous_parametric_process[k]->ident == LINEAR_MODEL) {
              switch (seq.type[k + 1]) {
              case INT_VALUE :
                residual = *pioutput[k] - (continuous_parametric_process[k]->observation[j]->intercept +
                            continuous_parametric_process[k]->observation[j]->slope *
                            (seq.index_param_type == IMPLICIT_TYPE ? i : seq.index_parameter[index][i]));
                break;
              case REAL_VALUE :
                residual = *proutput[k] - (continuous_parametric_process[k]->observation[j]->intercept +
                            continuous_parametric_process[k]->observation[j]->slope *
                            (seq.index_param_type == IMPLICIT_TYPE ? i : seq.index_parameter[index][i]));
                break;
              }

              entropy2 -= backward[i][j] * log(continuous_parametric_process[k]->observation[j]->mass_computation(residual , residual));
            }

            else if (continuous_parametric_process[k]->ident == AUTOREGRESSIVE_MODEL) {
              if (i == 0) {
                switch (seq.type[k + 1]) {
                case INT_VALUE :
                  residual = *pioutput[k] - continuous_parametric_process[k]->observation[j]->location;
                  break;
                case REAL_VALUE :
                  residual = *proutput[k] - continuous_parametric_process[k]->observation[j]->location;
                  break;
                }
              }

              else {
                switch (seq.type[k + 1]) {
                case INT_VALUE :
                  residual = *pioutput[k] - (continuous_parametric_process[k]->observation[j]->location +
                              continuous_parametric_process[k]->observation[j]->autoregressive_coeff *
                              (*(pioutput[k] - 1) - continuous_parametric_process[k]->observation[j]->location));
                  break;
                case REAL_VALUE :
                  residual = *proutput[k] - (continuous_parametric_process[k]->observation[j]->location +
                              continuous_parametric_process[k]->observation[j]->autoregressive_coeff *
                              (*(proutput[k] - 1) - continuous_parametric_process[k]->observation[j]->location));
                  break;
                }
              }

              entropy2 -= backward[i][j] * log(continuous_parametric_process[k]->observation[j]->mass_computation(residual , residual));
            }

            else {
              switch (seq.type[k + 1]) {
              case INT_VALUE :
                entropy2 -= backward[i][j] * log(continuous_parametric_process[k]->observation[j]->mass_computation(*pioutput[k] - seq.min_interval[k + 1] / 2 , *pioutput[k] + seq.min_interval[k + 1] / 2));
                break;
              case REAL_VALUE :
                entropy2 -= backward[i][j] * log(continuous_parametric_process[k]->observation[j]->mass_computation(*proutput[k] - seq.min_interval[k + 1] / 2 , *proutput[k] + seq.min_interval[k + 1] / 2));
                break;
              }
            }
          }
        }
      }
    }

    for (i = seq.length[index] - 2;i >= 0;i--) {
      for (j = 0;j < nb_output_process;j++) {
        switch (seq.type[j + 1]) {
        case INT_VALUE :
          pioutput[j]--;
          break;
        case REAL_VALUE :
          proutput[j]--;
          break;
        }
      }

      for (j = 0;j < nb_state;j++) {
        auxiliary[j] = 0.;

        switch (sojourn_type[j]) {

        // case semi-Markovian state

        case SEMI_MARKOVIAN : {
          occupancy = state_process->sojourn_time[j];
          obs_product = 1.;

          for (k = 1;k < MIN(seq.length[index] - i , occupancy->nb_value);k++) {
            obs_product *= observation[i + k][j] / norm[i + k];
            if (obs_product == 0.) {
              break;
            }

            occupancy_auxiliary[k] = 0.;

            if (backward1[i + k][j] > 0.) {
//            if (forward1[i + k][j] > 0.) {
              if (k < seq.length[index] - i - 1) {
                buff = backward1[i + k][j] * obs_product * occupancy->mass[k] /
                       forward1[i + k][j];
                occupancy_auxiliary[k] = buff * state_in[i][j];
                occupancy_entropy[j][k] += occupancy_auxiliary[k];

/*                if (occupancy->mass[k] > 0.) {
                  entropy2 -= occupancy_auxiliary[k] * log(occupancy->mass[k]);
                } */
              }

              else {
                buff = obs_product * (1. - occupancy->cumul[k - 1]);
                occupancy_auxiliary[k] = buff * state_in[i][j];
                if (occupancy->cumul[k - 1] < 1.) {
                  entropy2 -= occupancy_auxiliary[k] * log(1. - occupancy->cumul[k - 1]);
                }
              }

              auxiliary[j] += buff;
            }
          }

          sum = 0.;
          for (m = k - 1;m >= 1;m--) {
            sum += occupancy_auxiliary[m];
            if (backward[i + m][j] > 0.) {
              buff = sum / backward[i + m][j];
              if (buff > 0.) {
                state_entropy[i + m][j] += buff * (predicted_entropy[i][j] - log(buff));
              }
            }
          }
          break;
        }

        // case Markovian state

        case MARKOVIAN : {
          if (backward1[i + 1][j] > 0.) {
//          if (forward1[i + 1][j] > 0.) {
            auxiliary[j] = backward1[i + 1][j] / state_in[i][j];

/*            auxiliary[j] = backward1[i + 1][j] * observation[i + 1][j] /
                           (forward1[i + 1][j] * norm[i + 1]); */

            state_entropy[i + 1][j] = predicted_entropy[i][j];
          }
          break;
        }
        }
      }

      for (j = 0;j < nb_state;j++) {
        backward1[i][j] = 0.;

        for (k = 0;k < nb_state;k++) {
          buff = auxiliary[k] * transition[j][k] * forward1[i][j];
          backward1[i][j] += buff;
          transition_entropy[j][k] += buff;

/*          if (transition[j][k] > 0.) {
            entropy2 -= buff * log(transition[j][k]);
          } */
        }

        switch (sojourn_type[j]) {

        // case semi-Markovian state

        case SEMI_MARKOVIAN : {
          backward[i][j] = backward[i + 1][j] + backward1[i][j] - auxiliary[j] * state_in[i][j];
          if (backward[i][j] < 0.) {
            backward[i][j] = 0.;
          }
          if (backward[i][j] > 1.) {
            backward[i][j] = 1.;
          }
          break;
        }

        // case Markovian state

        case MARKOVIAN : {
          backward[i][j] = backward1[i][j];
          break;
        }
        }

        if (backward[i][j] > 0.) {
          for (k = 0;k < nb_output_process;k++) {
            if (categorical_process[k]) {
              if (categorical_process[k]->observation[j]->mass[*pioutput[k]] > 0.) {
                entropy2 -= backward[i][j] * log(categorical_process[k]->observation[j]->mass[*pioutput[k]]);
              }
            }

            else if (discrete_parametric_process[k]) {
              if (discrete_parametric_process[k]->observation[j]->mass[*pioutput[k]] > 0.) {
                entropy2 -= backward[i][j] * log(discrete_parametric_process[k]->observation[j]->mass[*pioutput[k]]);
              }
            }

            else {
              if (((continuous_parametric_process[k]->ident == GAMMA) ||
                  (continuous_parametric_process[k]->ident == ZERO_INFLATED_GAMMA)) && (seq.min_value[k + 1] < seq.min_interval[k + 1] / 2)) {
                switch (seq.type[k + 1]) {
                case INT_VALUE :
                  entropy2 -= backward[i][j] * log(continuous_parametric_process[k]->observation[j]->mass_computation(*pioutput[k] , *pioutput[k] + seq.min_interval[k + 1]));
                  break;
                case REAL_VALUE :
                  entropy2 -= backward[i][j] * log(continuous_parametric_process[k]->observation[j]->mass_computation(*proutput[k] , *proutput[k] + seq.min_interval[k + 1]));
                  break;
                }
              }

              else if (continuous_parametric_process[k]->ident == LINEAR_MODEL) {
                switch (seq.type[k + 1]) {
                case INT_VALUE :
                  residual = *pioutput[k] - (continuous_parametric_process[k]->observation[j]->intercept +
                              continuous_parametric_process[k]->observation[j]->slope *
                              (seq.index_param_type == IMPLICIT_TYPE ? i : seq.index_parameter[index][i]));
                  break;
                case REAL_VALUE :
                  residual = *proutput[k] - (continuous_parametric_process[k]->observation[j]->intercept +
                              continuous_parametric_process[k]->observation[j]->slope *
                              (seq.index_param_type == IMPLICIT_TYPE ? i : seq.index_parameter[index][i]));
                  break;
                }

                entropy2 -= backward[i][j] * log(continuous_parametric_process[k]->observation[j]->mass_computation(residual , residual));
              }

              else if (continuous_parametric_process[k]->ident == AUTOREGRESSIVE_MODEL) {
                if (i == 0) {
                  switch (seq.type[k + 1]) {
                  case INT_VALUE :
                    residual = *pioutput[k] - continuous_parametric_process[k]->observation[j]->location;
                    break;
                  case REAL_VALUE :
                    residual = *proutput[k] - continuous_parametric_process[k]->observation[j]->location;
                    break;
                  }
                }

                else {
                  switch (seq.type[k + 1]) {
                  case INT_VALUE :
                    residual = *pioutput[k] - (continuous_parametric_process[k]->observation[j]->location +
                                continuous_parametric_process[k]->observation[j]->autoregressive_coeff *
                                (*(pioutput[k] - 1) - continuous_parametric_process[k]->observation[j]->location));
                    break;
                  case REAL_VALUE :
                    residual = *proutput[k] - (continuous_parametric_process[k]->observation[j]->location +
                                continuous_parametric_process[k]->observation[j]->autoregressive_coeff *
                                (*(proutput[k] - 1) - continuous_parametric_process[k]->observation[j]->location));
                    break;
                  }
                }

                entropy2 -= backward[i][j] * log(continuous_parametric_process[k]->observation[j]->mass_computation(residual , residual));
              }

              else {
                switch (seq.type[k + 1]) {
                case INT_VALUE :
                  entropy2 -= backward[i][j] * log(continuous_parametric_process[k]->observation[j]->mass_computation(*pioutput[k] - seq.min_interval[k + 1] / 2 , *pioutput[k] + seq.min_interval[k + 1] / 2));
                  break;
                case REAL_VALUE :
                  entropy2 -= backward[i][j] * log(continuous_parametric_process[k]->observation[j]->mass_computation(*proutput[k] - seq.min_interval[k + 1] / 2 , *proutput[k] + seq.min_interval[k + 1] / 2));
                  break;
                }
              }
            }
          }
        }
      }

      switch (output) {

      case IN_STATE : {
        for (j = 0;j < nb_state;j++) {
          switch (sojourn_type[j]) {

          // case semi-Markovian state

          case SEMI_MARKOVIAN : {
            backward_output[i + 1][j] = auxiliary[j] * state_in[i][j];
            break;
          }

          // case Markovian state

          case MARKOVIAN : {
            backward_output[i + 1][j] = 0.;
            for (k = 0;k < nb_state;k++) {
              if (k != j) {
                backward_output[i + 1][j] += transition[k][j] * forward1[i][k];
              }
            }
            backward_output[i + 1][j] *= auxiliary[j];
            break;
          }
          }
        }
        break;
      }

      case OUT_STATE : {
        for (j = 0;j < nb_state;j++) {
          switch (sojourn_type[j]) {

          // case semi-Markovian state

          case SEMI_MARKOVIAN : {
            backward_output[i][j] = backward1[i][j];
            break;
          }

          // case Markovian state

          case MARKOVIAN : {
            backward_output[i][j] = 0.;
            for (k = 0;k < nb_state;k++) {
              if (k != j) {
                backward_output[i][j] += auxiliary[k] * transition[j][k];
              }
            }
            backward_output[i][j] *= forward1[i][j];
            break;
          }
          }
        }
        break;
      }
      }
    }

    if (output == IN_STATE) {
      for (i = 0;i < nb_state;i++) {
        backward_output[0][i] = backward[0][i];
      }
    }

    for (i = 0;i < nb_state;i++) {
      if (initial[i] > 0.) {
        entropy2 -= backward[0][i] * log(initial[i]);
      }
    }

    for (i = 0;i < nb_state;i++) {
      for (j = 0;j < nb_state;j++) {
        if (transition[i][j] > 0.) {
          entropy2 -= transition_entropy[i][j] * log(transition[i][j]);
        }
      }
    }

    for (i = 0;i < nb_state;i++) {
      if (sojourn_type[i] == SEMI_MARKOVIAN) {
        occupancy = state_process->sojourn_time[i];

        if (initial[i] > 0.) {
          obs_product = 1.;

#         ifdef DEBUG
          backward0[i] = 0.;
#         endif

          for (j = 1;j < MIN(seq.length[index] + 1 , occupancy->nb_value);j++) {
            obs_product *= observation[j - 1][i] / norm[j - 1];
            if (obs_product == 0.) {
              break;
            }

            occupancy_auxiliary[j] = 0.;

            if (backward1[j - 1][i] > 0.) {
//            if (forward1[j - 1][i] > 0.) {
              if (j < seq.length[index]) {
                switch (type) {

                case ORDINARY : {
                  occupancy_auxiliary[j] = backward1[j - 1][i] * obs_product * occupancy->mass[j] *
                                           initial[i] / forward1[j - 1][i];
                  occupancy_entropy[i][j] += occupancy_auxiliary[j];

/*                  if (occupancy->mass[j] > 0.) {
                    entropy2 -= occupancy_auxiliary[j] * log(occupancy->mass[j]);
                  } */
                  break;
                }

                case EQUILIBRIUM : {
                  occupancy_auxiliary[j] = backward1[j - 1][i] * obs_product * forward[i]->mass[j] *
                                           initial[i] / forward1[j - 1][i];
                  if (forward[i]->mass[j] > 0.) {
                    entropy2 -= occupancy_auxiliary[j] * log(forward[i]->mass[j]);
                  }
                  break;
                }
                }
              }

              else {
                switch (type) {

                case ORDINARY : {
                  occupancy_auxiliary[j] = obs_product * (1. - occupancy->cumul[j - 1]) * initial[i];
                  if (occupancy->cumul[j - 1] < 1.) {
                    entropy2 -= occupancy_auxiliary[j] * log(1. - occupancy->cumul[j - 1]);
                  }
                  break;
                }

                case EQUILIBRIUM : {
                  occupancy_auxiliary[j] = obs_product * (1. - forward[i]->cumul[j - 1]) * initial[i];
                  if (forward[i]->cumul[j - 1] < 1.) {
                    entropy2 -= occupancy_auxiliary[j] * log(1. - forward[i]->cumul[j - 1]);
                  }
                  break;
                }
                }
              }

#             ifdef DEBUG
              backward0[i] += occupancy_auxiliary[j];
#             endif

            }
          }

#         ifdef DEBUG
          cout << i << " " << backward[0][i] << " " << backward0[i] << endl;
#         endif

          sum = 0.;
          for (k = j - 1;k >= 1;k--) {
            sum += occupancy_auxiliary[k];
            if (backward[k - 1][i] > 0.) {
              buff = sum / backward[k - 1][i];
              if (buff > 0.) {
                state_entropy[k - 1][i] -= buff * log(buff);
              }
            }
          }
        }

        for (j = occupancy->offset;j < MIN(seq.length[index] , occupancy->nb_value);j++) {
          if (occupancy->mass[j] > 0.) {
            entropy2 -= occupancy_entropy[i][j] * log(occupancy->mass[j]);
          }
        }
      }
    }

    entropy2 += seq_likelihood;

#   ifdef MESSAGE
    if ((entropy2 < entropy1 - DOUBLE_ERROR) || (entropy2 > entropy1 + DOUBLE_ERROR)) {
      cout << "\nERROR: " << entropy1 << " " << entropy2 << endl;
    }
#   endif

    // restoration

    pstate = seq.int_sequence[index][0];

    for (i = 0;i < seq.length[index];i++) {
      backward_max = 0.;
      for (j = 0;j < nb_state;j++) {
        if (backward[i][j] > backward_max) {
          backward_max = backward[i][j];
          *pstate = j;
        }
      }

      pstate++;
    }

    seq.min_value[0] = 0;
    seq.max_value[0] = nb_state - 1;
    seq.build_marginal_frequency_distribution(0);

    state_seq_likelihood = SemiMarkov::likelihood_computation(seq , index);

    for (i = 0;i < seq.length[index];i++) {
      partial_entropy[i] = 0.;
      for (j = 0;j < nb_state;j++) {
        if (state_entropy[i][j] < 0.) {
          state_entropy[i][j] = 0.;
        }
        if (backward[i][j] > 0.) {
          partial_entropy[i] += backward[i][j] * (state_entropy[i][j] - log(backward[i][j]));
        }
      }
      if (partial_entropy[i] < 0.) {
        partial_entropy[i] = 0.;
      }
    }

    conditional_entropy[0] = partial_entropy[0];
    for (i = 1;i < seq.length[index];i++) {
      conditional_entropy[i] = partial_entropy[i] - partial_entropy[i - 1];
    }

    max_marginal_entropy = 0.;
    for (i = 0;i < seq.length[index];i++) {
      marginal_entropy[i] = 0.;
      for (j = 0;j < nb_state;j++) {
        if (backward[i][j] > 0.) {
          marginal_entropy[i] -= backward[i][j] * log(backward[i][j]);
        }
      }
      if (marginal_entropy[i] > max_marginal_entropy) {
        max_marginal_entropy = marginal_entropy[i];
      }
    }

    switch (format) {

    case ASCII : {
      switch (output) {
      case SSTATE :
        *os << "\n" << SEQ_label[SEQL_POSTERIOR_STATE_PROBABILITY] << "\n\n";
        break;
      case IN_STATE :
        *os << "\n" << SEQ_label[SEQL_POSTERIOR_IN_STATE_PROBABILITY] << "\n\n";
        break;
      case OUT_STATE :
        *os << "\n" << SEQ_label[SEQL_POSTERIOR_OUT_STATE_PROBABILITY] << "\n\n";
        break;
      }

//      seq.profile_ascii_print(*os , index , nb_state , backward_output ,
//                              STAT_label[STATL_STATE]);
      seq.profile_ascii_print(*os , index , nb_state , backward_output , conditional_entropy ,
                              marginal_entropy , partial_entropy);

      *os << "\n" << STAT_label[STATL_LIKELIHOOD] << ": " << seq_likelihood
          << "\n" << SEQ_label[SEQL_STATE_SEQUENCE_LIKELIHOOD] << ": " << state_seq_likelihood
          << "   (" << exp(state_seq_likelihood - seq_likelihood) << ")" << endl;
      break;
    }

    case SPREADSHEET : {
      switch (output) {
      case SSTATE :
        *os << "\n" << SEQ_label[SEQL_POSTERIOR_STATE_PROBABILITY] << "\n\n";
        break;
      case IN_STATE :
        *os << "\n" << SEQ_label[SEQL_POSTERIOR_IN_STATE_PROBABILITY] << "\n\n";
        break;
      case OUT_STATE :
        *os << "\n" << SEQ_label[SEQL_POSTERIOR_OUT_STATE_PROBABILITY] << "\n\n";
        break;
      }

//      seq.profile_spreadsheet_print(*os , index , nb_state , backward_output ,
//                                    STAT_label[STATL_STATE]);
      seq.profile_spreadsheet_print(*os , index , nb_state , backward_output , conditional_entropy ,
                                    marginal_entropy , partial_entropy);

      *os << "\n" << STAT_label[STATL_LIKELIHOOD] << "\t" << seq_likelihood
          << "\n" << SEQ_label[SEQL_STATE_SEQUENCE_LIKELIHOOD] << "\t" << state_seq_likelihood
          << "\t" << exp(state_seq_likelihood - seq_likelihood) << endl;
      break;
    }

    case GNUPLOT : {
//      seq.profile_plot_print(*os , index , nb_state , backward_output);
      seq.profile_plot_print(*os , index , nb_state , backward_output , conditional_entropy ,
                             marginal_entropy , partial_entropy);
      break;
    }

    case PLOT : {
      seq.profile_plotable_write((*plot_set)[1] , index , nb_state , backward_output);
      seq.entropy_profile_plotable_write((*plot_set)[2] , index , conditional_entropy , NULL ,
                                         marginal_entropy);
      seq.entropy_profile_plotable_write((*plot_set)[3] , index , partial_entropy);
      break;
    }
    }

    if (format != GNUPLOT) {
/*      double gini_index;

      gini_index = 0.;
      for (i = 0;i < seq.length[index];i++) {
        for (j = 0;j < nb_state;j++) {
          gini_index += backward[i][j] * (1. - backward[i][j]);
        }
      } */

      double entropy3 , nb_state_sequence;

      entropy3 = 0.;
      for (i = 0;i < seq.length[index];i++) {
        for (j = 0;j < nb_state;j++) {
          if (backward[i][j] > 0.) {
            entropy3 -= backward[i][j] * log(backward[i][j]);
          }
        }
      }

      // computation of the number of state sequences

      for (i = 0;i < nb_output_process;i++) {
        switch (seq.type[i + 1]) {
        case INT_VALUE :
          pioutput[i] = seq.int_sequence[index][i + 1];
          break;
        case REAL_VALUE :
          proutput[i] = seq.real_sequence[index][i + 1];
          break;
        }
      }

      // forward recurrence

      for (i = 0;i < seq.length[index];i++) {
        for (j = 0;j < nb_state;j++) {

          // computation of the indicator functions of the observation probabilities

          if (observation[i][j] > 0.) {
            observation[i][j] = 1.;
          }

          forward1[i][j] = 0.;

          switch (sojourn_type[j]) {

          // case semi-Markovian state

          case SEMI_MARKOVIAN : {
            occupancy = state_process->sojourn_time[j];

            if (i < seq.length[index] - 1) {
              for (k = 1;k <= MIN(i + 1 , occupancy->nb_value - 1);k++) {
                if (observation[i - k + 1][j] == 0.) {
                  break;
                }

                if (k < i + 1) {
                  if (occupancy->mass[k] > 0.) {
                    forward1[i][j] += state_in[i - k][j];
                  }
                }

                else {
                  if (initial[j] > 0.) {
                    switch (type) {

                    case ORDINARY : {
                      if (occupancy->mass[k] > 0.) {
                        forward1[i][j]++;
                      }
                      break;
                    }

                    case EQUILIBRIUM : {
                      if (forward[j]->mass[k] > 0.) {
                        forward1[i][j]++;
                      }
                      break;
                    }
                    }
                  }
                }
              }
            }

            else {
              for (k = 1;k <= MIN(i + 1 , occupancy->nb_value - 1);k++) {
                if (observation[i - k + 1][j] == 0.) {
                  break;
                }

                if (k < i + 1) {
                  if (1. - occupancy->cumul[k - 1] > 0.) {
                    forward1[i][j] += state_in[i - k][j];
                  }
                }

                else {
                  if (initial[j] > 0.) {
                    switch (type) {

                    case ORDINARY : {
                      if (1. - occupancy->cumul[k - 1] > 0.) {
                        forward1[i][j]++;
                      }
                      break;
                    }

                    case EQUILIBRIUM : {
                      if (1. - forward[j]->cumul[k - 1] > 0.) {
                        forward1[i][j]++;
                      }
                      break;
                    }
                    }
                  }
                }
              }
            }
            break;
          }

          // case Markovian state

          case MARKOVIAN : {
            if (observation[i][j] == 1.) {
              if (i == 0) {
                if (initial[j] > 0.) {
                  forward1[i][j] = 1.;
                }
              }
              else {
                forward1[i][j] = state_in[i - 1][j];
              }
            }
            break;
          }
          }
        }

        if (i < seq.length[index] - 1) {
          for (j = 0;j < nb_state;j++) {
            state_in[i][j] = 0.;
            for (k = 0;k < nb_state;k++) {
              if (transition[k][j] > 0.) {
                state_in[i][j] += forward1[i][k];
              }
            }
          }
        }

        for (j = 0;j < nb_output_process;j++) {
          switch (seq.type[j + 1]) {
          case INT_VALUE :
            pioutput[j]++;
            break;
          case REAL_VALUE :
            proutput[j]++;
            break;
          }
        }
      }

      nb_state_sequence = 0.;
      i = seq.length[index] - 1;
      for (j = 0;j < nb_state;j++) {
        nb_state_sequence += forward1[i][j];
      }

      switch (format) {
      case ASCII :
/*        *os << "\n" << SEQ_label[SEQL_GINI_INDEX] << ": " << gini_index << " ("
            << gini_index / seq.length[index] << ")   " << SEQ_label[SEQL_UPPER_BOUND] << ": "
            << seq.length[index] * (1. - 1. / nb_state) << " (" << 1. - 1. / nb_state
        *os << ")   " << SEQ_label[SEQL_UPPER_BOUND] << ": "
            << seq.length[index] * log((double)nb_state) << " (" << log((double)nb_state) */
        *os << "\n" << SEQ_label[SEQL_STATE_SEQUENCE_ENTROPY] << ": " << entropy1
            << " (" << entropy1 / seq.length[index] << ")   " << SEQ_label[SEQL_UPPER_BOUND] << ": "
            << log((double)nb_state_sequence) << " ("
            << log((double)nb_state_sequence) / seq.length[index]
            << ")\n" << SEQ_label[SEQL_MARGINAL_ENTROPY_SUM] << ": " << entropy3 << " ("
            << entropy3 / seq.length[index] << ")\n\n"
            << SEQ_label[SEQL_NB_STATE_SEQUENCE] << ": " << nb_state_sequence << endl;
        break;
      case SPREADSHEET :
/*        *os << "\n" << SEQ_label[SEQL_GINI_INDEX] << "\t" << gini_index << "\t"
            << gini_index / seq.length[index] << "\t" << SEQ_label[SEQL_UPPER_BOUND] << "\t"
            << seq.length[index] * (1. - 1. / nb_state) << "\t" << 1. - 1. / nb_state
        *os << "\t" << SEQ_label[SEQL_UPPER_BOUND] << "\t"
            << seq.length[index] * log((double)nb_state) << "\t" << log((double)nb_state) */
        *os << "\n" << SEQ_label[SEQL_STATE_SEQUENCE_ENTROPY] << "\t" << entropy1
            << "\t" << entropy1 / seq.length[index] << "\t" << SEQ_label[SEQL_UPPER_BOUND] << "\t"
            << log((double)nb_state_sequence) << "\t"
            << log((double)nb_state_sequence) / seq.length[index]
            << "\n" << SEQ_label[SEQL_MARGINAL_ENTROPY_SUM] << "\t" << entropy3 << "\t"
            << entropy3 / seq.length[index] << "\n\n"
            << SEQ_label[SEQL_NB_STATE_SEQUENCE] << "\t" << nb_state_sequence << endl;
        break;
      }

#     ifdef DEBUG
      int state;
      double min_nb_state_sequence , smoothed_proba , cumul_smoothed_proba ,
             max_smoothed_proba , **backward2;

      // backward recurrence

      min_nb_state_sequence = nb_state_sequence;

      backward2 = new double*[seq.length[index]];
      for (i = 0;i < seq.length[index];i++) {
        backward2[i] = new double[nb_state];
      }

      i = seq.length[index] - 1;
      for (j = 0;j < nb_state;j++) {
        backward2[i][j] = forward1[i][j];
        backward1[i][j] = 1.;
      }

      for (i = seq.length[index] - 2;i >= 0;i--) {
        for (j = 0;j < nb_state;j++) {
          auxiliary[j] = 0.;

          switch (sojourn_type[j]) {

          // case semi-Markovian state

          case SEMI_MARKOVIAN : {
            occupancy = state_process->sojourn_time[j];

            for (k = 1;k < MIN(seq.length[index] - i , occupancy->nb_value);k++) {
              if (observation[i + k][j] == 0.) {
                break;
              }

              if (k < seq.length[index] - i - 1) {
                if (occupancy->mass[k] > 0.) {
                  auxiliary[j] += backward1[i + k][j];
                }
              }
              else {
                if (1. - occupancy->cumul[k - 1] > 0.) {
                  auxiliary[j]++;
                }
              }
            }
            break;
          }

          // case Markovian state

          case MARKOVIAN : {
            if (observation[i + 1][j] == 1.) {
              auxiliary[j] = backward1[i + 1][j];
            }
            break;
          }
          }
        }

        for (j = 0;j < nb_state;j++) {
          backward1[i][j] = 0.;

          for (k = 0;k < nb_state;k++) {
            if (transition[j][k] > 0.) {
              backward1[i][j] += auxiliary[k];
            }
          }

          switch (sojourn_type[j]) {

          // case semi-Markovian state

          case SEMI_MARKOVIAN : {

#           ifdef DEBUG
            if ((i == 0) && (initial[j] > 0.)) {
              occupancy = state_process->sojourn_time[j];
              backward0[j] = 0.;

              for (k = 1;k < MIN(seq.length[index] + 1 , occupancy->nb_value);k++) {
                if (observation[k - 1][j] == 0.) {
                  break;
                }

                if (k < seq.length[index]) {
                  if (occupancy->mass[k] > 0.) {
                    backward0[j] += backward1[k - 1][j];
                  }
                }
                else {
                  if (1. - occupancy->cumul[k - 1] > 0.) {
                    backward0[j]++;
                  }
                }
              }
            }
#           endif

            backward2[i][j] = backward2[i + 1][j] + backward1[i][j] * forward1[i][j] -
                              auxiliary[j] * state_in[i][j];

#           ifdef DEBUG
            if ((i == 0) && (initial[j] > 0.)) {
              cout << j << " " << backward2[i][j] << " " << backward0[j] << endl;
            }
#           endif

            break;
          }

          // case Markovian state

          case MARKOVIAN : {
            backward2[i][j] = backward1[i][j] * forward1[i][j];
            break;
          }
          }
        }

        smoothed_proba = 1.1;
        cumul_smoothed_proba = 0.;
        nb_state_sequence = 0;

        for (j = 0;j < nb_state;j++) {
          max_smoothed_proba = 0.;
          for (k = 0;k < nb_state;k++) {
            if ((backward[i][k] > max_smoothed_proba) && (backward[i][k] < smoothed_proba)) {
              max_smoothed_proba = backward[i][k];
              state = k;
            }
          }
          cumul_smoothed_proba += max_smoothed_proba;
          nb_state_sequence += backward2[i][state];

          if (cumul_smoothed_proba < 1. - MIN_SMOOTHED_PROBABILITY) {
            smoothed_proba = max_smoothed_proba;
          }
          else {
            break;
          }
        }

        if (nb_state_sequence < min_nb_state_sequence) {
          min_nb_state_sequence = nb_state_sequence;
        }
      }

      cout << SEQ_label[SEQL_NB_STATE_SEQUENCE]
           << " (" << 1. - MIN_SMOOTHED_PROBABILITY << " beam)"
           << ": " << min_nb_state_sequence << endl;

      cout << "\n";
      for (i = 0;i < seq.length[index];i++) {
        obs_product = 0.;
        for (j = 0;j < nb_state;j++) {
          cout << backward2[i][j] << " (" << backward[i][j] << ")  ";
          obs_product += backward2[i][j];
        }
        cout << "| " << obs_product << endl;
      }

      for (i = 0;i < seq.length[index];i++) {
        delete [] backward2[i];
      }
      delete [] backward2;
#     endif

    }
  }

  for (i = 0;i < seq.length[index];i++) {
    delete [] observation[i];
  }
  delete [] observation;

  delete [] norm;
  delete [] state_norm;

  for (i = 0;i < seq.length[index];i++) {
    delete [] forward1[i];
  }
  delete [] forward1;

  for (i = 0;i < seq.length[index] - 1;i++) {
    delete [] state_in[i];
  }
  delete [] state_in;

  for (i = 0;i < seq.length[index];i++) {
    delete [] backward[i];
  }
  delete [] backward;

  for (i = 0;i < seq.length[index];i++) {
    delete [] backward1[i];
  }
  delete [] backward1;

  delete [] auxiliary;
  delete [] occupancy_auxiliary;

  if (output != SSTATE) {
    for (i = 0;i < seq.length[index];i++) {
      delete [] backward_output[i];
    }
    delete [] backward_output;
  }

  delete [] transition_predicted;
  delete [] occupancy_predicted;

  for (i = 0;i < seq.length[index];i++) {
    delete [] state_entropy[i];
  }
  delete [] state_entropy;

  for (i = 0;i < seq.length[index];i++) {
    delete [] predicted_entropy[i];
  }
  delete [] predicted_entropy;

  for (i = 0;i < nb_state;i++) {
    delete [] transition_entropy[i];
  }
  delete [] transition_entropy;

  for (i = 0;i < nb_state;i++) {
    delete [] occupancy_entropy[i];
  }
  delete [] occupancy_entropy;

  delete [] partial_entropy;
  delete [] conditional_entropy;
  delete [] marginal_entropy;

# ifdef DEBUG
  delete [] backward0;
# endif

  delete [] pioutput;
  delete [] proutput;

  return seq_likelihood;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Simulation of state sequences for an observed sequence using
 *         the forward-backward algorithm for sampling.
 *
 *  \param[in] seq               reference on a MarkovianSequences object,
 *  \param[in] index             sequence index,
 *  \param[in] os                stream,
 *  \param[in] format            file format (ASCII/SPREADSHEET),
 *  \param[in] nb_state_sequence number of state sequences.
 *
 *  \return                      log-likelihood for the observed sequence.
 */
/*--------------------------------------------------------------*/

double HiddenSemiMarkov::forward_backward_sampling(const MarkovianSequences &seq , int index ,
                                                   ostream &os , output_format format ,
                                                   int nb_state_sequence) const

{
  int i , j , k;
  int state_occupancy , *pstate , **pioutput;
  double seq_likelihood , state_seq_likelihood , obs_product , residual , **observation , *norm ,
         *state_norm , **forward1 , **state_in , *backward , *cumul_backward , **proutput;
  DiscreteParametric *occupancy;

# ifdef DEBUG
  int m;
  double sum;
# endif


  // initializations

  observation = new double*[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    observation[i] = new double[nb_state];
  }

  norm = new double[seq.length[index]];
  state_norm = new double[nb_state];

  forward1 = new double*[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    forward1[i] = new double[nb_state];
  }

  state_in = new double*[seq.length[index] - 1];
  for (i = 0;i < seq.length[index] - 1;i++) {
    state_in[i] = new double[nb_state];
  }

  backward = new double[seq.length[index] + 1];
  cumul_backward = new double[seq.length[index] + 1];

  pioutput = new int*[nb_output_process];
  proutput = new double*[nb_output_process];

  for (i = 0;i < nb_output_process;i++) {
    switch (seq.type[i + 1]) {
    case INT_VALUE :
      pioutput[i] = seq.int_sequence[index][i + 1];
      break;
    case REAL_VALUE :
      proutput[i] = seq.real_sequence[index][i + 1];
      break;
    }
  }

# ifdef DEBUG
  double **state_sequence_probability;


  state_sequence_probability = new double*[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    state_sequence_probability[i] = new double[nb_state];
    for (j = 0;j < nb_state;j++) {
      state_sequence_probability[i][j] = 0.;
    }
  }
# endif

  // forward recurrence

  seq_likelihood = 0.;
  for (i = 0;i < seq.length[index];i++) {
    norm[i] = 0.;

    for (j = 0;j < nb_state;j++) {

      // computation of the observation probabilities

      observation[i][j] = 1.;
      for (k = 0;k < nb_output_process;k++) {
        if (categorical_process[k]) {
          observation[i][j] *= categorical_process[k]->observation[j]->mass[*pioutput[k]];
        }

        else if (discrete_parametric_process[k]) {
          observation[i][j] *= discrete_parametric_process[k]->observation[j]->mass[*pioutput[k]];
        }

        else {
          if (((continuous_parametric_process[k]->ident == GAMMA) ||
               (continuous_parametric_process[k]->ident == ZERO_INFLATED_GAMMA)) && (seq.min_value[k + 1] < seq.min_interval[k + 1] / 2)) {
            switch (seq.type[k + 1]) {
            case INT_VALUE :
              observation[i][j] *= continuous_parametric_process[k]->observation[j]->mass_computation(*pioutput[k] , *pioutput[k] + seq.min_interval[k + 1]);
              break;
            case REAL_VALUE :
              observation[i][j] *= continuous_parametric_process[k]->observation[j]->mass_computation(*proutput[k] , *proutput[k] + seq.min_interval[k + 1]);
              break;
            }
          }

          else if (continuous_parametric_process[k]->ident == LINEAR_MODEL) {
            switch (seq.type[k + 1]) {
            case INT_VALUE :
              residual = *pioutput[k] - (continuous_parametric_process[k]->observation[j]->intercept +
                          continuous_parametric_process[k]->observation[j]->slope *
                          (seq.index_param_type == IMPLICIT_TYPE ? i : seq.index_parameter[index][i]));
              break;
            case REAL_VALUE :
              residual = *proutput[k] - (continuous_parametric_process[k]->observation[j]->intercept +
                          continuous_parametric_process[k]->observation[j]->slope *
                          (seq.index_param_type == IMPLICIT_TYPE ? i : seq.index_parameter[index][i]));
              break;
            }

            observation[i][j] *= continuous_parametric_process[k]->observation[j]->mass_computation(residual , residual);
          }

          else if (continuous_parametric_process[k]->ident == AUTOREGRESSIVE_MODEL) {
            if (i == 0) {
              switch (seq.type[k + 1]) {
              case INT_VALUE :
                residual = *pioutput[k] - continuous_parametric_process[k]->observation[j]->location;
                break;
              case REAL_VALUE :
                residual = *proutput[k] - continuous_parametric_process[k]->observation[j]->location;
                break;
              }
            }

            else {
              switch (seq.type[k + 1]) {
              case INT_VALUE :
                residual = *pioutput[k] - (continuous_parametric_process[k]->observation[j]->location +
                            continuous_parametric_process[k]->observation[j]->autoregressive_coeff *
                            (*(pioutput[k] - 1) - continuous_parametric_process[k]->observation[j]->location));
                break;
              case REAL_VALUE :
                residual = *proutput[k] - (continuous_parametric_process[k]->observation[j]->location +
                            continuous_parametric_process[k]->observation[j]->autoregressive_coeff *
                            (*(proutput[k] - 1) - continuous_parametric_process[k]->observation[j]->location));
                break;
              }
            }

            observation[i][j] *= continuous_parametric_process[k]->observation[j]->mass_computation(residual , residual);
          }

          else {
            switch (seq.type[k + 1]) {
            case INT_VALUE :
              observation[i][j] *= continuous_parametric_process[k]->observation[j]->mass_computation(*pioutput[k] - seq.min_interval[k + 1] / 2 , *pioutput[k] + seq.min_interval[k + 1] / 2);
              break;
            case REAL_VALUE :
              observation[i][j] *= continuous_parametric_process[k]->observation[j]->mass_computation(*proutput[k] - seq.min_interval[k + 1] / 2 , *proutput[k] + seq.min_interval[k + 1] / 2);
              break;
            }
          }
        }
      }

      switch (sojourn_type[j]) {

      // case semi-Markovian state 

      case SEMI_MARKOVIAN : {
        if (i == 0) {
          state_norm[j] = initial[j];
        }
        else {
          state_norm[j] += state_in[i - 1][j] - forward1[i - 1][j];
        }
        state_norm[j] *= observation[i][j];

        norm[i] += state_norm[j];
        break;
      }

      // case Markovian state

      case MARKOVIAN : {
        if (i == 0) {
          forward1[i][j] = initial[j];
        }
        else {
          forward1[i][j] = state_in[i - 1][j];
        }
        forward1[i][j] *= observation[i][j];

        norm[i] += forward1[i][j];
        break;
      }
      }
    }

    if (norm[i] > 0.) {
      for (j = 0;j < nb_state;j++) {
        switch (sojourn_type[j]) {
        case SEMI_MARKOVIAN :
          state_norm[j] /= norm[i];
          break;
        case MARKOVIAN :
          forward1[i][j] /= norm[i];
          break;
        }
      }

      seq_likelihood += log(norm[i]);
    }

    else {
      seq_likelihood = D_INF;
      break;
    }

    for (j = 0;j < nb_state;j++) {

      // case semi-Markovian state

      if (sojourn_type[j] == SEMI_MARKOVIAN) {
        occupancy = state_process->sojourn_time[j];
        obs_product = 1.;
        forward1[i][j] = 0.;

        if (i < seq.length[index] - 1) {
          for (k = 1;k <= MIN(i + 1 , occupancy->nb_value - 1);k++) {
            obs_product *= observation[i - k + 1][j] / norm[i - k + 1];
            if (obs_product == 0.) {
              break;
            }

            if (k < i + 1) {
              forward1[i][j] += obs_product * occupancy->mass[k] * state_in[i - k][j];
            }

            else {
              switch (type) {
              case ORDINARY :
                forward1[i][j] += obs_product * occupancy->mass[k] * initial[j];
                break;
              case EQUILIBRIUM :
                forward1[i][j] += obs_product * forward[j]->mass[k] * initial[j];
                break;
              }
            }
          }
        }

        else {
          for (k = 1;k <= MIN(i + 1 , occupancy->nb_value - 1);k++) {
            obs_product *= observation[i - k + 1][j] / norm[i - k + 1];
            if (obs_product == 0.) {
              break;
            }

            if (k < i + 1) {
              forward1[i][j] += obs_product * (1. - occupancy->cumul[k - 1]) * state_in[i - k][j];
            }

            else {
              switch (type) {
              case ORDINARY :
                forward1[i][j] += obs_product * (1. - occupancy->cumul[k - 1]) * initial[j];
                break;
              case EQUILIBRIUM :
                forward1[i][j] += obs_product * (1. - forward[j]->cumul[k - 1]) * initial[j];
                break;
              }
            }
          }
        }
      }
    }

    if (i < seq.length[index] - 1) {
      for (j = 0;j < nb_state;j++) {
        state_in[i][j] = 0.;
        for (k = 0;k < nb_state;k++) {
          state_in[i][j] += transition[k][j] * forward1[i][k];
        }
      }
    }

    for (j = 0;j < nb_output_process;j++) {
      switch (seq.type[j + 1]) {
      case INT_VALUE :
        pioutput[j]++;
        break;
      case REAL_VALUE :
        proutput[j]++;
        break;
      }
    }
  }

  if (seq_likelihood != D_INF) {

    // backward passes

#   ifdef MESSAGE
    cout << "\n";
#   endif

    for (i = 0;i < nb_state_sequence;i++) {
      j = seq.length[index] - 1;
      pstate = seq.int_sequence[index][0] + j;
      stat_tool::cumul_computation(nb_state , forward1[j] , cumul_backward);
      *pstate = cumul_method(nb_state , cumul_backward);

      do {

        // case semi-Markovian state

        if (sojourn_type[*pstate] == SEMI_MARKOVIAN) {
          occupancy = state_process->sojourn_time[*pstate];
          obs_product = 1.;

          if (j < seq.length[index] - 1) {
            for (k = 1;k <= MIN(j + 1 , occupancy->nb_value - 1);k++) {
              obs_product *= observation[j - k + 1][*pstate] / norm[j - k + 1];
              if (obs_product == 0.) {
                break;
              }

              if (k < j + 1) {
                backward[k] = obs_product * occupancy->mass[k] * state_in[j - k][*pstate] /
                              forward1[j][*pstate];
              }

              else {
                switch (type) {
                case ORDINARY :
                  backward[k] = obs_product * occupancy->mass[k] * initial[*pstate] /
                                forward1[j][*pstate];
                  break;
                case EQUILIBRIUM :
                  backward[k] = obs_product * forward[*pstate]->mass[k] * initial[*pstate] /
                                forward1[j][*pstate];
                  break;
                }
              }
            }
          }

          else {
            for (k = 1;k <= MIN(j + 1 , occupancy->nb_value - 1);k++) {
              obs_product *= observation[j - k + 1][*pstate] / norm[j - k + 1];
              if (obs_product == 0.) {
                break;
              }

              if (k < j + 1) {
                backward[k] = obs_product * (1. - occupancy->cumul[k - 1]) * state_in[j - k][*pstate] /
                              forward1[j][*pstate];
              }

              else {
                switch (type) {
                case ORDINARY :
                  backward[k] = obs_product * (1. - occupancy->cumul[k - 1]) * initial[*pstate] /
                                forward1[j][*pstate];
                  break;
                case EQUILIBRIUM :
                  backward[k] = obs_product * (1. - forward[*pstate]->cumul[k - 1]) * initial[*pstate] /
                                forward1[j][*pstate];
                  break;
                }
              }
            }
          }

          stat_tool::cumul_computation(k - 1 , backward + 1 , cumul_backward);
          state_occupancy = 1 + cumul_method(k - 1 , cumul_backward);

#         ifdef DEBUG
          sum = 0.;
          for (m = 1;m < k;m++) {
            sum += backward[m];
          }
          if ((sum < 1. - DOUBLE_ERROR) || (sum > 1. + DOUBLE_ERROR)) {
            cout << "\nERROR: " << j << " " << sum << endl;
          }
#         endif

          for (k = 1;k < state_occupancy;k++) {
            pstate--;
            *pstate = *(pstate + 1);
          }
          j -= (state_occupancy - 1);

          if (j == 0) {
            break;
          }
        }

        j--;
        for (k = 0;k < nb_state;k++) {
          backward[k] = transition[k][*pstate] * forward1[j][k] / state_in[j][*pstate];
        }
        stat_tool::cumul_computation(nb_state , backward , cumul_backward);
        *--pstate = cumul_method(nb_state , cumul_backward);

#       ifdef DEBUG
        sum = 0.;
        for (k = 0;k < nb_state;k++) {
          sum += backward[k];
        }
        if ((sum < 1. - DOUBLE_ERROR) || (sum > 1. + DOUBLE_ERROR)) {
          cout << "\nERROR: " << j << " " << sum << endl;
        }
#       endif

      }
      while (j > 0);

#     ifdef DEBUG
      pstate = seq.int_sequence[index][0];
      for (j = 0;j < seq.length[index];j++) {
        state_sequence_probability[j][*pstate++]++;
      }
#     endif

#     ifdef MESSAGE
      state_seq_likelihood = SemiMarkov::likelihood_computation(seq , index);

      pstate = seq.int_sequence[index][0];

      switch (format) {

      case ASCII : {
        for (j = 0;j < seq.length[index];j++) {
          os << *pstate++ << " ";
        }

        os << "  " << i + 1 << "  " << state_seq_likelihood
           << "   (" << exp(state_seq_likelihood - seq_likelihood) << ")" << endl;
        break;
      }

      case SPREADSHEET : {
        for (j = 0;j < seq.length[index];j++) {
          os << *pstate++ << "\t";
        }

        os << "\t" << i + 1 << "\t" << state_seq_likelihood
           << "\t" << exp(state_seq_likelihood - seq_likelihood) << endl;
        break;
      }
      }
#     endif

    }

#   ifdef DEBUG
    if (nb_state_sequence >= 1000) {
      for (i = 0;i < seq.length[index];i++) {
        for (j = 0;j < nb_state;j++) {
          state_sequence_probability[i][j] /= nb_state_sequence;
        }
      }

      pstate = seq.int_sequence[index][0];
      for (j = 0;j < seq.length[index];j++) {
        *pstate++ = I_DEFAULT;
      }

      os << "\n" << SEQ_label[SEQL_POSTERIOR_STATE_PROBABILITY] << "\n\n";
      seq.profile_ascii_print(os , index , nb_state , state_sequence_probability ,
                              STAT_label[STATL_STATE]);
    }
#   endif

  }

  for (i = 0;i < seq.length[index];i++) {
    delete [] observation[i];
  }
  delete [] observation;

  delete [] norm;
  delete [] state_norm;

  for (i = 0;i < seq.length[index];i++) {
    delete [] forward1[i];
  }
  delete [] forward1;

  for (i = 0;i < seq.length[index] - 1;i++) {
    delete [] state_in[i];
  }
  delete [] state_in;

  delete [] backward;
  delete [] cumul_backward;

  delete [] pioutput;
  delete [] proutput;

# ifdef DEBUG
  for (i = 0;i < seq.length[index];i++) {
    delete [] state_sequence_probability[i];
  }
  delete [] state_sequence_probability;
# endif

  return seq_likelihood;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the log-parameters of a hidden semi-Markov chain.
 */
/*--------------------------------------------------------------*/

void HiddenSemiMarkov::log_computation()

{
  int i , j;
  double *pcumul;
  DiscreteParametric *occupancy;


  Chain::log_computation();

  for (i = 0;i < nb_state;i++) {
    if (sojourn_type[i] == SEMI_MARKOVIAN) {
      occupancy = state_process->sojourn_time[i];

      if (occupancy->mass[occupancy->offset] > 0.) {
        stat_tool::log_computation(occupancy->nb_value , occupancy->mass , occupancy->mass);

        pcumul = occupancy->cumul;
        for (j = 0;j < occupancy->nb_value;j++) {
          *pcumul = 1. - *pcumul;
          pcumul++;
        }
        stat_tool::log_computation(occupancy->nb_value , occupancy->cumul , occupancy->cumul);

        if (type == EQUILIBRIUM) {
          stat_tool::log_computation(forward[i]->nb_value , forward[i]->mass , forward[i]->mass);

          pcumul = forward[i]->cumul;
          for (j = 0;j < forward[i]->nb_value;j++) {
            *pcumul = 1. - *pcumul;
            pcumul++;
          }
          stat_tool::log_computation(forward[i]->nb_value , forward[i]->cumul , forward[i]->cumul);
        }
      }
    }
  }

  for (i = 0;i < nb_output_process;i++) {
    if (categorical_process[i]) {
      for (j = 0;j < nb_state;j++) {
        categorical_process[i]->observation[j]->log_computation();
      }
    }

    else if (discrete_parametric_process[i]) {
      for (j = 0;j < nb_state;j++) {
        stat_tool::log_computation(discrete_parametric_process[i]->nb_value ,
                                   discrete_parametric_process[i]->observation[j]->mass ,
                                   discrete_parametric_process[i]->observation[j]->cumul);
      }
    }
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the most probable state sequences using the Viterbi algorithm.
 *
 *  \param[in] seq                   reference on a MarkovianSequences object,
 *  \param[in] posterior_probability pointer on the posterior probabilities of the most probable state sequences,
 *  \param[in] index                 sequence index.
 *
 *  \return                          log-likelihood for the most probable state sequences.
 */
/*--------------------------------------------------------------*/

double HiddenSemiMarkov::viterbi(const MarkovianSequences &seq ,
                                 double *posterior_probability , int index) const

{
  int i , j , k , m;
  int length , *pstate , **pioutput , **input_state , **optimal_state ,
      **optimal_occupancy;
  double likelihood = 0. , obs_product , buff , residual , forward_max , **observation ,
         *forward1 , **state_in , **proutput;
  DiscreteParametric *occupancy;


  // initializations

  length = (index == I_DEFAULT ? seq.max_length : seq.length[index]);

  observation = new double*[length];
  for (i = 0;i < length;i++) {
    observation[i] = new double[nb_state];
  }

  forward1 = new double[nb_state];

  state_in = new double*[length - 1];
  for (i = 0;i < length - 1;i++) {
    state_in[i] = new double[nb_state];
  }

  input_state = new int*[length - 1];
  for (i = 0;i < length - 1;i++) {
    input_state[i] = new int[nb_state];
  }

  optimal_state = new int*[length];
  for (i = 0;i < length;i++) {
    optimal_state[i] = new int[nb_state];
  }

  optimal_occupancy = new int*[length];
  for (i = 0;i < length;i++) {
    optimal_occupancy[i] = new int[nb_state];
  }

  pioutput = new int*[nb_output_process];
  proutput = new double*[nb_output_process];

  for (i = 0;i < seq.nb_sequence;i++) {
    if ((index == I_DEFAULT) || (index == i)) {
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

      // forward recurrence

      for (j = 0;j < seq.length[i];j++) {
        for (k = 0;k < nb_state;k++) {

          // computation of the observation probabilities

          observation[j][k] = 0.;
          for (m = 0;m < nb_output_process;m++) {
            if (categorical_process[m]) {
              buff = categorical_process[m]->observation[k]->cumul[*pioutput[m]];
            }

            else if (discrete_parametric_process[m]) {
              buff = discrete_parametric_process[m]->observation[k]->cumul[*pioutput[m]];
            }

            else {
              if (((continuous_parametric_process[m]->ident == GAMMA) ||
                   (continuous_parametric_process[m]->ident == ZERO_INFLATED_GAMMA)) && (seq.min_value[m + 1] < seq.min_interval[m + 1] / 2)) {
                switch (seq.type[m + 1]) {
                case INT_VALUE :
                  buff = continuous_parametric_process[m]->observation[k]->mass_computation(*pioutput[m] , *pioutput[m] + seq.min_interval[m + 1]);
                  break;
                case REAL_VALUE :
                  buff = continuous_parametric_process[m]->observation[k]->mass_computation(*proutput[m] , *proutput[m] + seq.min_interval[m + 1]);
                  break;
                }
              }

              else if (continuous_parametric_process[m]->ident == LINEAR_MODEL) {
                switch (seq.type[m + 1]) {
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

                buff = continuous_parametric_process[m]->observation[k]->mass_computation(residual , residual);
              }

              else if (continuous_parametric_process[m]->ident == AUTOREGRESSIVE_MODEL) {
                if (j == 0) {
                  switch (seq.type[m + 1]) {
                  case INT_VALUE :
                    residual = *pioutput[m] - continuous_parametric_process[m]->observation[k]->location;
                    break;
                  case REAL_VALUE :
                    residual = *proutput[m] - continuous_parametric_process[m]->observation[k]->location;
                    break;
                  }
                }

                else {
                  switch (seq.type[m + 1]) {
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

                buff = continuous_parametric_process[m]->observation[k]->mass_computation(residual , residual);
              }

              else {
                switch (seq.type[m + 1]) {
                case INT_VALUE :
                  buff = continuous_parametric_process[m]->observation[k]->mass_computation(*pioutput[m] - seq.min_interval[m + 1] / 2 , *pioutput[m] + seq.min_interval[m + 1] / 2);
                  break;
                case REAL_VALUE :
                  buff = continuous_parametric_process[m]->observation[k]->mass_computation(*proutput[m] - seq.min_interval[m + 1] / 2 , *proutput[m] + seq.min_interval[m + 1] / 2);
                  break;
                }
              }

              if (buff > 0.) {
                buff = log(buff);
              }
              else {
                buff = D_INF;
              }
            }

            if (buff == D_INF) {
              observation[j][k] = D_INF;
              break;
            }
            else {
              observation[j][k] += buff;
            }
          }

          switch (sojourn_type[k]) {

          // case semi-Markovian state

          case SEMI_MARKOVIAN : {
            occupancy = state_process->sojourn_time[k];
            obs_product = 0.;
            forward1[k] = D_INF;

            if (j < seq.length[i] - 1) {
              for (m = 1;m <= MIN(j + 1 , occupancy->nb_value - 1);m++) {
                if (observation[j - m + 1][k] == D_INF) {
                  break;
                }
                else {
                  obs_product += observation[j - m + 1][k];
                }

                if (m < j + 1) {
                  buff = obs_product + occupancy->mass[m] + state_in[j - m][k];
                }

                else {
                  switch (type) {
                  case ORDINARY :
                    buff = obs_product + occupancy->mass[m] + cumul_initial[k];
                    break;
                  case EQUILIBRIUM :
                    buff = obs_product + forward[k]->mass[m] + cumul_initial[k];
                    break;
                  }
                }

                if (buff > forward1[k]) {
                  forward1[k] = buff;
                  if (m < j + 1) {
                    optimal_state[j][k] = input_state[j - m][k];
                  }
                  optimal_occupancy[j][k] = m;
                }
              }
            }

            else {
              for (m = 1;m <= MIN(j + 1 , occupancy->nb_value - 1);m++) {
                if (observation[j - m + 1][k] == D_INF) {
                  break;
                }
                else {
                  obs_product += observation[j - m + 1][k];
                }

                if (m < j + 1) {
                  buff = obs_product + occupancy->cumul[m - 1] + state_in[j - m][k];
                }

                else {
                  switch (type) {
                  case ORDINARY :
                    buff = obs_product + occupancy->cumul[m - 1] + cumul_initial[k];
                    break;
                  case EQUILIBRIUM :
                    buff = obs_product + forward[k]->cumul[m - 1] + cumul_initial[k];
                    break;
                  }
                }

                if (buff > forward1[k]) {
                  forward1[k] = buff;
                  if (m < j + 1) {
                    optimal_state[j][k] = input_state[j - m][k];
                  }
                  optimal_occupancy[j][k] = m;
                }
              }
            }
            break;
          }

          // case Markovian state

          case MARKOVIAN : {
            if (j == 0) {
              forward1[k] = cumul_initial[k];
            }
            else {
              forward1[k] = state_in[j - 1][k];
              optimal_state[j][k] = input_state[j - 1][k];
            }
            optimal_occupancy[j][k] = 1;

            if (forward1[k] != D_INF) {
              if (observation[j][k] == D_INF) {
                forward1[k] = D_INF;
              }
              else {
                forward1[k] += observation[j][k];
              }
            }
            break;
          }
          }
        }

#       ifdef DEBUG
        cout << j << " : ";
        for (k = 0;k < nb_state;k++) {
          cout << forward1[k];
          if (forward1[k] != D_INF) {
            cout << " " << optimal_occupancy[j][k] << " " << optimal_state[j][k];
          }
          cout << " | ";
        }
        cout << endl;
#       endif

        if (j < seq.length[i] - 1) {
          for (k = 0;k < nb_state;k++) {
            state_in[j][k] = D_INF;
            for (m = 0;m < nb_state;m++) {
              buff = cumul_transition[m][k] + forward1[m];
              if (buff > state_in[j][k]) {
                state_in[j][k] = buff;
                input_state[j][k] = m;
              }
            }
          }
        }

        for (k = 0;k < nb_output_process;k++) {
          switch (seq.type[k + 1]) {
          case INT_VALUE :
            pioutput[k]++;
            break;
          case REAL_VALUE :
            proutput[k]++;
            break;
          }
        }
      }

      // extraction of the log-likelihood for the most probable state sequence

      pstate = seq.int_sequence[i][0] + seq.length[i] - 1;
      forward_max = D_INF;

      for (j = 0;j < nb_state;j++) {
        if (forward1[j] > forward_max) {
          forward_max = forward1[j];
          *pstate = j;
        }
      }

      if (forward_max != D_INF) {
        likelihood += forward_max;
        if (posterior_probability) {
          posterior_probability[i] = forward_max;
        }
      }

      else {
        likelihood = D_INF;
        if (posterior_probability) {
          posterior_probability[i] = 0.;
        }
        break;
      }

      // restoration of the most probable state sequence

      j = seq.length[i] - 1;

      do {
        for (k = 0;k < optimal_occupancy[j][*pstate] - 1;k++) {
          pstate--;
          *pstate = *(pstate + 1);
        }

        if (j >= optimal_occupancy[j][*pstate]) {
          pstate--;
          *pstate = optimal_state[j][*(pstate + 1)];
          j -= optimal_occupancy[j][*(pstate + 1)];
        }
        else {
          j -= optimal_occupancy[j][*pstate];
        }
      }
      while (j >= 0);
    }
  }

  for (i = 0;i < length;i++) {
    delete [] observation[i];
  }
  delete [] observation;

  delete [] forward1;

  for (i = 0;i < length - 1;i++) {
    delete [] state_in[i];
  }
  delete [] state_in;

  for (i = 0;i < length - 1;i++) {
    delete [] input_state[i];
  }
  delete [] input_state;

  for (i = 0;i < length;i++) {
    delete [] optimal_state[i];
  }
  delete [] optimal_state;

  for (i = 0;i < length;i++) {
    delete [] optimal_occupancy[i];
  }
  delete [] optimal_occupancy;

  delete [] pioutput;
  delete [] proutput;

  return likelihood;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the most probable state sequences using the Viterbi algorithm.
 *
 *  \param[in] seq reference on a SemiMarkovData object.
 */
/*--------------------------------------------------------------*/

void HiddenSemiMarkov::viterbi(SemiMarkovData &seq) const

{
  seq.posterior_probability = new double[seq.nb_sequence];
  seq.restoration_likelihood = viterbi(seq , seq.posterior_probability);
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the N most probable state sequences for
 *         an observed sequence using the generalized Viterbi algorithm.
 *
 *  \param[in] seq                reference on a MarkovianSequences object,
 *  \param[in] index              sequence index,
 *  \param[in] os                 stream,
 *  \param[in] seq_likelihood     log-likelihood for the observed sequence,
 *  \param[in] format             file format (ASCII/SPREADSHEET),
 *  \param[in] inb_state_sequence number of state sequences.
 *
 *  \return                       log-likelihood for the most probable state sequence.
 */
/*--------------------------------------------------------------*/

double HiddenSemiMarkov::generalized_viterbi(const MarkovianSequences &seq , int index ,
                                             ostream &os , double seq_likelihood ,
                                             output_format format , int inb_state_sequence) const

{
  bool **active_cell;
  int i , j , k , m;
  int nb_state_sequence , max_occupancy , brank , previous_rank , nb_cell , *rank ,
      *pstate , **pioutput , ***input_state , ***optimal_state , ***optimal_occupancy ,
      ***input_rank , ***optimal_rank;
  double buff , residual , forward_max , state_seq_likelihood , likelihood_cumul ,
         *obs_product , **observation , **forward1 , **proutput , ***state_in;
  DiscreteParametric *occupancy;


  // initializations

  observation = new double*[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    observation[i] = new double[nb_state];
  }

  obs_product = new double[seq.length[index] + 1];

  forward1 = new double*[nb_state];
  for (i = 0;i < nb_state;i++) {
    forward1[i] = new double[inb_state_sequence];
  }

  state_in = new double**[seq.length[index] - 1];
  for (i = 0;i < seq.length[index] - 1;i++) {
    state_in[i] = new double*[nb_state];
    for (j = 0;j < nb_state;j++) {
      state_in[i][j] = new double[inb_state_sequence];
    }
  }

  rank = new int[MAX(seq.length[index] + 1 , nb_state)];

  input_state = new int**[seq.length[index] - 1];
  for (i = 0;i < seq.length[index] - 1;i++) {
    input_state[i] = new int*[nb_state];
    for (j = 0;j < nb_state;j++) {
      input_state[i][j] = new int[inb_state_sequence];
    }
  }

  optimal_state = new int**[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    optimal_state[i] = new int*[nb_state];
    for (j = 0;j < nb_state;j++) {
      optimal_state[i][j] = new int[inb_state_sequence];
    }
  }

  optimal_occupancy = new int**[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    optimal_occupancy[i] = new int*[nb_state];
    for (j = 0;j < nb_state;j++) {
      optimal_occupancy[i][j] = new int[inb_state_sequence];
    }
  }

  input_rank = new int**[seq.length[index] - 1];
  for (i = 0;i < seq.length[index] - 1;i++) {
    input_rank[i] = new int*[nb_state];
    for (j = 0;j < nb_state;j++) {
      input_rank[i][j] = new int[inb_state_sequence];
    }
  }

  optimal_rank = new int**[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    optimal_rank[i] = new int*[nb_state];
    for (j = 0;j < nb_state;j++) {
      optimal_rank[i][j] = new int[inb_state_sequence];
    }
  }

  active_cell = new bool*[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    active_cell[i] = new bool[nb_state];
    for (j = 0;j < nb_state;j++) {
      active_cell[i][j] = false;
    }
  }

  pioutput = new int*[nb_output_process];
  proutput = new double*[nb_output_process];

  for (i = 0;i < nb_output_process;i++) {
    switch (seq.type[i + 1]) {
    case INT_VALUE :
      pioutput[i] = seq.int_sequence[index][i + 1];
      break;
    case REAL_VALUE :
      proutput[i] = seq.real_sequence[index][i + 1];
      break;
    }
  }
  nb_state_sequence = 1;

# ifdef DEBUG
  double entropy = 0. , **state_sequence_probability;


  state_sequence_probability = new double*[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    state_sequence_probability[i] = new double[nb_state];
    for (j = 0;j < nb_state;j++) {
//      state_sequence_probability[i][j] = 0.;
      state_sequence_probability[i][j] = D_INF;
    }
  }
# endif

  // forward recurrence

  for (i = 0;i < seq.length[index];i++) {
    for (j = 0;j < nb_state;j++) {

      // computation of the observation probabilities

      observation[i][j] = 0.;
      for (k = 0;k < nb_output_process;k++) {
        if (categorical_process[k]) {
          buff = categorical_process[k]->observation[j]->cumul[*pioutput[k]];
        }

        else if (discrete_parametric_process[k]) {
          buff = discrete_parametric_process[k]->observation[j]->cumul[*pioutput[k]];
        }

        else {
          if (((continuous_parametric_process[k]->ident == GAMMA) ||
               (continuous_parametric_process[k]->ident == ZERO_INFLATED_GAMMA)) && (seq.min_value[k + 1] < seq.min_interval[k + 1] / 2)) {
            switch (seq.type[k + 1]) {
            case INT_VALUE :
              buff = continuous_parametric_process[k]->observation[j]->mass_computation(*pioutput[k] , *pioutput[k] + seq.min_interval[k + 1]);
              break;
            case REAL_VALUE :
              buff = continuous_parametric_process[k]->observation[j]->mass_computation(*proutput[k] , *proutput[k] + seq.min_interval[k + 1]);
              break;
            }
          }

          else if (continuous_parametric_process[k]->ident == LINEAR_MODEL) {
            switch (seq.type[k + 1]) {
            case INT_VALUE :
              residual = *pioutput[k] - (continuous_parametric_process[k]->observation[j]->intercept +
                          continuous_parametric_process[k]->observation[j]->slope *
                          (seq.index_param_type == IMPLICIT_TYPE ? i : seq.index_parameter[index][i]));
              break;
            case REAL_VALUE :
              residual = *proutput[k] - (continuous_parametric_process[k]->observation[j]->intercept +
                          continuous_parametric_process[k]->observation[j]->slope *
                          (seq.index_param_type == IMPLICIT_TYPE ? i : seq.index_parameter[index][i]));
              break;
            }

            buff = continuous_parametric_process[k]->observation[j]->mass_computation(residual , residual);
          }

          else if (continuous_parametric_process[k]->ident == AUTOREGRESSIVE_MODEL) {
            if (i == 0) {
              switch (seq.type[k + 1]) {
              case INT_VALUE :
                residual = *pioutput[k] - continuous_parametric_process[k]->observation[j]->location;
                break;
              case REAL_VALUE :
                residual = *proutput[k] - continuous_parametric_process[k]->observation[j]->location;
                break;
              }
            }

            else {
              switch (seq.type[k + 1]) {
              case INT_VALUE :
                residual = *pioutput[k] - (continuous_parametric_process[k]->observation[j]->location +
                            continuous_parametric_process[k]->observation[j]->autoregressive_coeff *
                            (*(pioutput[k] - 1) - continuous_parametric_process[k]->observation[j]->location));
                break;
              case REAL_VALUE :
                residual = *proutput[k] - (continuous_parametric_process[k]->observation[j]->location +
                            continuous_parametric_process[k]->observation[j]->autoregressive_coeff *
                            (*(proutput[k] - 1) - continuous_parametric_process[k]->observation[j]->location));
                break;
              }
            }

            buff = continuous_parametric_process[k]->observation[j]->mass_computation(residual , residual);
          }

          else {
            switch (seq.type[k + 1]) {
            case INT_VALUE :
              buff = continuous_parametric_process[k]->observation[j]->mass_computation(*pioutput[k] - seq.min_interval[k + 1] / 2 , *pioutput[k] + seq.min_interval[k + 1] / 2);
              break;
            case REAL_VALUE :
              buff = continuous_parametric_process[k]->observation[j]->mass_computation(*proutput[k] - seq.min_interval[k + 1] / 2 , *proutput[k] + seq.min_interval[k + 1] / 2);
              break;
            }
          }

          if (buff > 0.) {
            buff = log(buff);
          }
          else {
            buff = D_INF;
          }
        }

        if (buff == D_INF) {
          observation[i][j] = D_INF;
          break;
        }
        else {
          observation[i][j] += buff;
        }
      }

      switch (sojourn_type[j]) {

      // case semi-Markovian state

      case SEMI_MARKOVIAN : {
        occupancy = state_process->sojourn_time[j];

        obs_product[0] = 0.;
        for (k = 1;k <= MIN(i + 1 , occupancy->nb_value - 1);k++) {
          if (observation[i - k + 1][j] == D_INF) {
            break;
          }
          else {
            obs_product[k] = obs_product[k - 1] + observation[i - k + 1][j];
          }
        }
        max_occupancy = k - 1;

        for (k = 1;k <= max_occupancy;k++) {
          rank[k] = 0;
        }

        for (k = 0;k < nb_state_sequence;k++) {
          forward1[j][k] = D_INF;

          if (i < seq.length[index] - 1) {
            for (m = 1;m <= max_occupancy;m++) {
              if (m < i + 1) {
                buff = obs_product[m] + occupancy->mass[m] + state_in[i - m][j][rank[m]];
              }

              else {
                if (rank[i + 1] == 0) {
                  switch (type) {
                  case ORDINARY :
                    buff = obs_product[m] + occupancy->mass[m] + cumul_initial[j];
                    break;
                  case EQUILIBRIUM :
                    buff = obs_product[m] + forward[j]->mass[m] + cumul_initial[j];
                    break;
                  }
                }

                else {
                  buff = D_INF;
                }
              }

              if (buff > forward1[j][k]) {
                forward1[j][k] = buff;
                if (m < i + 1) {
                  optimal_state[i][j][k] = input_state[i - m][j][rank[m]];
                  optimal_rank[i][j][k] = input_rank[i - m][j][rank[m]];
                }
                optimal_occupancy[i][j][k] = m;
              }
            }
          }

          else {
            for (m = 1;m <= max_occupancy;m++) {
              if (m < i + 1) {
                buff = obs_product[m] + occupancy->cumul[m - 1] + state_in[i - m][j][rank[m]];
              }

              else {
                if (rank[i + 1] == 0) {
                  switch (type) {
                  case ORDINARY :
                    buff = obs_product[m] + occupancy->cumul[m - 1] + cumul_initial[j];
                    break;
                  case EQUILIBRIUM :
                    buff = obs_product[m] + forward[j]->cumul[m - 1] + cumul_initial[j];
                    break;
                  }
                }

                else {
                  buff = D_INF;
                }
              }

              if (buff > forward1[j][k]) {
                forward1[j][k] = buff;
                if (m < i + 1) {
                  optimal_state[i][j][k] = input_state[i - m][j][rank[m]];
                  optimal_rank[i][j][k] = input_rank[i - m][j][rank[m]];
                }
                optimal_occupancy[i][j][k] = m;
              }
            }
          }

          if (forward1[j][k] != D_INF) {
            rank[optimal_occupancy[i][j][k]]++;
          }
        }
        break;
      }

      // case Markovian state

      case MARKOVIAN : {
        for (k = 0;k < nb_state_sequence;k++) {
          if (i == 0) {
            forward1[j][k] = cumul_initial[j];
          }
          else {
            forward1[j][k] = state_in[i - 1][j][k];
            optimal_state[i][j][k] = input_state[i - 1][j][k];
            optimal_rank[i][j][k] = input_rank[i - 1][j][k];
          }
          optimal_occupancy[i][j][k] = 1;

          if (forward1[j][k] != D_INF) {
            if (observation[i][j] == D_INF) {
              forward1[j][k] = D_INF;
            }
            else {
              forward1[j][k] += observation[i][j];
            }
          }
        }
        break;
      }
      }

      for (k = nb_state_sequence;k < inb_state_sequence;k++) {
        forward1[j][k] = D_INF;
      }
    }

#   ifdef DEBUG
    cout << i << " : ";
    for (j = 0;j < nb_state;j++) {
      cout << j << " :";
      for (k = 0;k < nb_state_sequence;k++) {
        cout << " " << forward1[j][k];
        if (forward1[j][k] != D_INF) {
          cout << " " << optimal_occupancy[i][j][k];
          if (optimal_occupancy[i][j][k] < i + 1) {
            cout << " " << optimal_state[i][j][k] << " " << optimal_rank[i][j][k];
          }
        }
        cout << " |";
      }
      cout << "| ";
    }
    cout << endl;
#   endif

    if (i < seq.length[index] - 1) {
      if (nb_state_sequence < inb_state_sequence) {
        if (nb_state_sequence * nb_state < inb_state_sequence) {
          nb_state_sequence *= nb_state;
        }
        else {
          nb_state_sequence = inb_state_sequence;
        }
      }

      for (j = 0;j < nb_state;j++) {
        for (k = 0;k < nb_state;k++) {
          rank[k] = 0;
        }

        for (k = 0;k < nb_state_sequence;k++) {
          state_in[i][j][k] = D_INF;
          for (m = 0;m < nb_state;m++) {
            buff = cumul_transition[m][j] + forward1[m][rank[m]];
            if (buff > state_in[i][j][k]) {
              state_in[i][j][k] = buff;
              input_state[i][j][k] = m;
              input_rank[i][j][k] = rank[m];
            }
          }

          if (state_in[i][j][k] != D_INF) {
            rank[input_state[i][j][k]]++;
          }
        }

        for (k = nb_state_sequence;k < inb_state_sequence;k++) {
          state_in[i][j][k] = D_INF;
        }
      }

#     ifdef DEBUG
      cout << i << " : ";
      for (j = 0;j < nb_state;j++) {
        cout << j << " :";
        for (k = 0;k < nb_state_sequence;k++) {
          cout << " " << state_in[i][j][k];
          if (state_in[i][j][k] != D_INF) {
            cout << " " << input_state[i][j][k] << " " << input_rank[i][j][k];
          }
          cout << " |";
        }
        cout << "| ";
      }
      cout << endl;
#     endif

    }

    for (j = 0;j < nb_output_process;j++) {
      switch (seq.type[j + 1]) {
      case INT_VALUE :
        pioutput[j]++;
        break;
      case REAL_VALUE :
        proutput[j]++;
        break;
      }
    }
  }

  // extraction of the log-likelihood for the most probable state sequence

  for (i = 0;i < nb_state;i++) {
    rank[i] = 0;
  }
  likelihood_cumul = 0.;

  for (i = 0;i < nb_state_sequence;i++) {
    pstate = seq.int_sequence[index][0] + seq.length[index] - 1;
    forward_max = D_INF;

    for (j = 0;j < nb_state;j++) {
      if (forward1[j][rank[j]] > forward_max) {
        forward_max = forward1[j][rank[j]];
        *pstate = j;
      }
    }

    if (i == 0) {
      state_seq_likelihood = forward_max;
    }

    if (forward_max == D_INF) {
      break;
    }

    // restoration of the most probable state sequence

    brank = rank[*pstate];
    rank[*pstate]++;
    j = seq.length[index] - 1;

#   ifdef DEBUG
    cout << "\n" << *pstate << " " << optimal_occupancy[j][*pstate][brank] << " " << brank << " | ";
#   endif

    do {
      for (k = 0;k < optimal_occupancy[j][*pstate][brank];k++) {
        active_cell[j - k][*pstate] = true;
      }

      for (k = 0;k < optimal_occupancy[j][*pstate][brank] - 1;k++) {
        pstate--;
        *pstate = *(pstate + 1);
      }

      if (j >= optimal_occupancy[j][*pstate][brank]) {
        pstate--;
        *pstate = optimal_state[j][*(pstate + 1)][brank];
        previous_rank = optimal_rank[j][*(pstate + 1)][brank];
        j -= optimal_occupancy[j][*(pstate + 1)][brank];
        brank = previous_rank;

#       ifdef DEBUG
        cout << *pstate << " " << optimal_occupancy[j][*pstate][brank] << " " << brank << " | ";
#       endif

      }
      else {
        j -= optimal_occupancy[j][*pstate][brank];
      }
    }
    while (j >= 0);

#   ifdef DEBUG
    cout << endl;
#   endif

    likelihood_cumul += exp(forward_max);

#   ifdef DEBUG
    pstate = seq.int_sequence[index][0];
    for (j = 0;j < seq.length[index];j++) {
/*      state_sequence_probability[j][*pstate++] += exp(forward_max - seq_likelihood); */

      if (forward_max > state_sequence_probability[j][*pstate]) {
        state_sequence_probability[j][*pstate] = forward_max;
      }
      pstate++;
    }
#   endif

    nb_cell = 0;
    for (j = 0;j < seq.length[index];j++) {
      for (k = 0;k < nb_state;k++) {
        if (active_cell[j][k]) {
          nb_cell++;
        }
      }
    }

#   ifdef MESSAGE
    if (i == 0) {
      os << "\n";
    }

    pstate = seq.int_sequence[index][0];

    switch (format) {

    case ASCII : {
      for (j = 0;j < seq.length[index];j++) {
        os << *pstate++ << " ";
      }

//      os << "  " << i + 1 << "  " << forward_max << "   (" << exp(forward_max - state_seq_likelihood)
      os << "  " << i + 1 << "  " << forward_max << "   (" << exp(forward_max - seq_likelihood)
         << "  " << likelihood_cumul / exp(seq_likelihood) << "  " << nb_cell << ")" << endl;

      if (nb_component == nb_state) {
        os << SEQ_label[SEQL_STATE_BEGIN] << ": ";

        pstate = seq.int_sequence[index][0] + 1;
        if (seq.index_parameter) {
          for (j = 1;j < seq.length[index];j++) {
            if (*pstate != *(pstate - 1)) {
              os << seq.index_parameter[index][j] << ", ";
            }
            pstate++;
          }
        }

        else {
          for (j = 1;j < seq.length[index];j++) {
            if (*pstate != *(pstate - 1)) {
              os << j << ", ";
            }
            pstate++;
          }
        }

        os << endl;
      }
      break;
    }

    case SPREADSHEET : {
      for (j = 0;j < seq.length[index];j++) {
        os << *pstate++ << "\t";
      }

//      os << "\t" << i + 1 << "\t" << forward_max << "\t" << exp(forward_max - state_seq_likelihood)
      os << "\t" << i + 1 << "\t" << forward_max << "\t" << exp(forward_max - seq_likelihood)
         << "\t" << likelihood_cumul / exp(seq_likelihood) << "\t" << nb_cell << endl;
      break;
    }
    }
#   endif

#   ifdef DEBUG
    entropy -= exp(forward_max - seq_likelihood) * forward_max;
#   endif

  }

# ifdef DEBUG
  os << "\n" << SEQ_label[SEQL_STATE_SEQUENCE_ENTROPY] << ": " << entropy + seq_likelihood << endl;

  if (likelihood_cumul / exp(seq_likelihood) > 0.8) {
    for (i = 0;i < seq.length[index];i++) {
      for (j = 0;j < nb_state;j++) {
        if (state_sequence_probability[i][j] != D_INF) {
          state_sequence_probability[i][j] = exp(state_sequence_probability[i][j] - seq_likelihood);
        }
        else {
          state_sequence_probability[i][j] = 0.;
        }
      }
    }

    pstate = seq.int_sequence[index][0];
    for (j = 0;j < seq.length[index];j++) {
      *pstate++ = I_DEFAULT;
    }

//    os << "\n" << SEQ_label[SEQL_POSTERIOR_STATE_PROBABILITY] << "\n\n";
    os << "\n" << SEQ_label[SEQL_MAX_POSTERIOR_STATE_PROBABILITY] << "\n\n";
    seq.profile_ascii_print(os , index , nb_state , state_sequence_probability ,
                            STAT_label[STATL_STATE]);
  }
# endif

  for (i = 0;i < seq.length[index];i++) {
    delete [] observation[i];
  }
  delete [] observation;

  delete [] obs_product;

  for (i = 0;i < nb_state;i++) {
    delete [] forward1[i];
  }
  delete [] forward1;

  for (i = 0;i < seq.length[index] - 1;i++) {
    for (j = 0;j < nb_state;j++) {
      delete [] state_in[i][j];
    }
    delete [] state_in[i];
  }
  delete [] state_in;

  delete [] rank;

  for (i = 0;i < seq.length[index] - 1;i++) {
    for (j = 0;j < nb_state;j++) {
      delete [] input_state[i][j];
    }
    delete [] input_state[i];
  }
  delete [] input_state;

  for (i = 0;i < seq.length[index];i++) {
    for (j = 0;j < nb_state;j++) {
      delete [] optimal_state[i][j];
    }
    delete [] optimal_state[i];
  }
  delete [] optimal_state;

  for (i = 0;i < seq.length[index];i++) {
    for (j = 0;j < nb_state;j++) {
      delete [] optimal_occupancy[i][j];
    }
    delete [] optimal_occupancy[i];
  }
  delete [] optimal_occupancy;

  for (i = 0;i < seq.length[index] - 1;i++) {
    for (j = 0;j < nb_state;j++) {
      delete [] input_rank[i][j];
    }
    delete [] input_rank[i];
  }
  delete [] input_rank;

  for (i = 0;i < seq.length[index];i++) {
    for (j = 0;j < nb_state;j++) {
      delete [] optimal_rank[i][j];
    }
    delete [] optimal_rank[i];
  }
  delete [] optimal_rank;

  for (i = 0;i < seq.length[index];i++) {
    delete [] active_cell[i];
  }
  delete [] active_cell;

  delete [] pioutput;
  delete [] proutput;

# ifdef DEBUG
  for (i = 0;i < seq.length[index];i++) {
    delete [] state_sequence_probability[i];
  }
  delete [] state_sequence_probability;
# endif

  return state_seq_likelihood;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of state profiles using the Viterbi forward-backward algorithm.
 *
 *  \param[in] seq            reference on a MarkovianSequences object,
 *  \param[in] index          sequence index,
 *  \param[in] os             stream,
 *  \param[in] plot           pointer on a MultiPlot object,
 *  \param[in] output         output type,
 *  \param[in] format         output format (ASCII/SPREADSHEET/GNUPLOT/PLOT),
 *  \param[in] seq_likelihood log-likelihood for the observed sequence.
 *
 *  \return                   log-likelihood for the most probable state sequence.
 */
/*--------------------------------------------------------------*/

double HiddenSemiMarkov::viterbi_forward_backward(const MarkovianSequences &seq , int index ,
                                                  ostream *os , MultiPlot *plot ,
                                                  state_profile output , output_format format ,
                                                  double seq_likelihood) const

{
  int i , j , k , m;
  int *pstate , **pioutput;
  double obs_product , buff , residual , state_seq_likelihood , backward_max , **observation ,
         **forward1 , **state_in , **backward , **backward1 , *auxiliary , *occupancy_auxiliary ,
         **backward_output , **proutput;
  DiscreteParametric *occupancy;


  // initializations

  observation = new double*[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    observation[i] = new double[nb_state];
  }

  forward1 = new double*[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    forward1[i] = new double[nb_state];
  }

  state_in = new double*[seq.length[index] - 1];
  for (i = 0;i < seq.length[index] - 1;i++) {
    state_in[i] = new double[nb_state];
  }

  backward = new double*[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    backward[i] = new double[nb_state];
  }

  backward1 = new double*[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    backward1[i] = new double[nb_state];
  }

  auxiliary = new double[nb_state];
  occupancy_auxiliary = new double[seq.length[index] + 1];

  if (output == SSTATE) {
    backward_output = backward;
  }
  else {
    backward_output = new double*[seq.length[index]];
    for (i = 0;i < seq.length[index];i++) {
      backward_output[i] = new double[nb_state];
    }
  }

  pioutput = new int*[nb_output_process];
  proutput = new double*[nb_output_process];

# ifdef MESSAGE
  int *state_sequence , **input_state , **optimal_state , **optimal_forward_occupancy;

  input_state = new int*[seq.length[index] - 1];
  for (i = 0;i < seq.length[index] - 1;i++) {
    input_state[i] = new int[nb_state];
  }

  optimal_state = new int*[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    optimal_state[i] = new int[nb_state];
  }

  optimal_forward_occupancy = new int*[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    optimal_forward_occupancy[i] = new int[nb_state];
  }

  state_sequence = new int[seq.length[index]];
# endif

  for (i = 0;i < nb_output_process;i++) {
    switch (seq.type[i + 1]) {
    case INT_VALUE :
      pioutput[i] = seq.int_sequence[index][i + 1];
      break;
    case REAL_VALUE :
      proutput[i] = seq.real_sequence[index][i + 1];
      break;
    }
  }

  // forward recurrence

  for (i = 0;i < seq.length[index];i++) {
    for (j = 0;j < nb_state;j++) {

      // computation of the observation probabilities

      observation[i][j] = 0.;
      for (k = 0;k < nb_output_process;k++) {
        if (categorical_process[k]) {
          buff = categorical_process[k]->observation[j]->cumul[*pioutput[k]];
        }

        else if (discrete_parametric_process[k]) {
          buff = discrete_parametric_process[k]->observation[j]->cumul[*pioutput[k]];
        }

        else {
          if (((continuous_parametric_process[k]->ident == GAMMA) ||
               (continuous_parametric_process[k]->ident == ZERO_INFLATED_GAMMA)) && (seq.min_value[k + 1] < seq.min_interval[k + 1] / 2)) {
            switch (seq.type[k + 1]) {
            case INT_VALUE :
              buff = continuous_parametric_process[k]->observation[j]->mass_computation(*pioutput[k] , *pioutput[k] + seq.min_interval[k + 1]);
              break;
            case REAL_VALUE :
              buff = continuous_parametric_process[k]->observation[j]->mass_computation(*proutput[k] , *proutput[k] + seq.min_interval[k + 1]);
              break;
            }
          }

          else if (continuous_parametric_process[k]->ident == LINEAR_MODEL) {
            switch (seq.type[k + 1]) {
            case INT_VALUE :
              residual = *pioutput[k] - (continuous_parametric_process[k]->observation[j]->intercept +
                          continuous_parametric_process[k]->observation[j]->slope *
                          (seq.index_param_type == IMPLICIT_TYPE ? i : seq.index_parameter[index][i]));
              break;
            case REAL_VALUE :
              residual = *proutput[k] - (continuous_parametric_process[k]->observation[j]->intercept +
                          continuous_parametric_process[k]->observation[j]->slope *
                          (seq.index_param_type == IMPLICIT_TYPE ? i : seq.index_parameter[index][i]));
              break;
            }

            buff = continuous_parametric_process[k]->observation[j]->mass_computation(residual , residual);
          }

          else if (continuous_parametric_process[k]->ident == AUTOREGRESSIVE_MODEL) {
            if (i == 0) {
              switch (seq.type[k + 1]) {
              case INT_VALUE :
                residual = *pioutput[k] - continuous_parametric_process[k]->observation[j]->location;
                break;
              case REAL_VALUE :
                residual = *proutput[k] - continuous_parametric_process[k]->observation[j]->location;
                break;
              }
            }

            else {
              switch (seq.type[k + 1]) {
              case INT_VALUE :
                residual = *pioutput[k] - (continuous_parametric_process[k]->observation[j]->location +
                            continuous_parametric_process[k]->observation[j]->autoregressive_coeff *
                            (*(pioutput[k] - 1) - continuous_parametric_process[k]->observation[j]->location));
                break;
              case REAL_VALUE :
                residual = *proutput[k] - (continuous_parametric_process[k]->observation[j]->location +
                            continuous_parametric_process[k]->observation[j]->autoregressive_coeff *
                            (*(proutput[k] - 1) - continuous_parametric_process[k]->observation[j]->location));
                break;
              }
            }

            buff = continuous_parametric_process[k]->observation[j]->mass_computation(residual , residual);
          }

          else {
            switch (seq.type[k + 1]) {
            case INT_VALUE :
              buff = continuous_parametric_process[k]->observation[j]->mass_computation(*pioutput[k] - seq.min_interval[k + 1] / 2 , *pioutput[k] + seq.min_interval[k + 1] / 2);
              break;
            case REAL_VALUE :
              buff = continuous_parametric_process[k]->observation[j]->mass_computation(*proutput[k] - seq.min_interval[k + 1] / 2 , *proutput[k] + seq.min_interval[k + 1] / 2);
              break;
            }
          }

          if (buff > 0.) {
            buff = log(buff);
          }
          else {
            buff = D_INF;
          }
        }

        if (buff == D_INF) {
          observation[i][j] = D_INF;
          break;
        }
        else {
          observation[i][j] += buff;
        }
      }

      switch (sojourn_type[j]) {

      // case semi-Markovian state

      case SEMI_MARKOVIAN : {
        occupancy = state_process->sojourn_time[j];
        obs_product = 0.;
        forward1[i][j] = D_INF;

        if (i < seq.length[index] - 1) {
          for (k = 1;k <= MIN(i + 1 , occupancy->nb_value - 1);k++) {
            if (observation[i - k + 1][j] == D_INF) {
              break;
            }
            else {
              obs_product += observation[i - k + 1][j];
            }

            if (k < i + 1) {
              buff = obs_product + occupancy->mass[k] + state_in[i - k][j];
            }

            else {
              switch (type) {
              case ORDINARY :
                buff = obs_product + occupancy->mass[k] + cumul_initial[j];
                break;
              case EQUILIBRIUM :
                buff = obs_product + forward[j]->mass[k] + cumul_initial[j];
                break;
              }
            }

            if (buff > forward1[i][j]) {
              forward1[i][j] = buff;

#             ifdef MESSAGE
              if (k < i + 1) {
                optimal_state[i][j] = input_state[i - k][j];
              }
              optimal_forward_occupancy[i][j] = k;
#             endif

            }
          }
        }

        else {
          for (k = 1;k <= MIN(i + 1 , occupancy->nb_value - 1);k++) {
            if (observation[i - k + 1][j] == D_INF) {
              break;
            }
            else {
              obs_product += observation[i - k + 1][j];
            }

            if (k < i + 1) {
              buff = obs_product + occupancy->cumul[k - 1] + state_in[i - k][j];
            }

            else {
              switch (type) {
              case ORDINARY :
                buff = obs_product + occupancy->cumul[k - 1] + cumul_initial[j];
                break;
              case EQUILIBRIUM :
                buff = obs_product + forward[j]->cumul[k - 1] + cumul_initial[j];
                break;
              }
            }

            if (buff > forward1[i][j]) {
              forward1[i][j] = buff;

#             ifdef MESSAGE
              if (k < i + 1) {
                optimal_state[i][j] = input_state[i - k][j];
              }
              optimal_forward_occupancy[i][j] = k;
#             endif

            }
          }
        }
        break;
      }

      // case Markovian state

      case MARKOVIAN : {
        if (i == 0) {
          forward1[i][j] = cumul_initial[j];
        }
        else {
          forward1[i][j] = state_in[i - 1][j];

#         ifdef MESSAGE
          optimal_state[i][j] = input_state[i - 1][j];
#         endif

        }

#       ifdef MESSAGE
        optimal_forward_occupancy[i][j] = 1;
#       endif

        if (forward1[i][j] != D_INF) {
          if (observation[i][j] == D_INF) {
            forward1[i][j] = D_INF;
          }
          else {
            forward1[i][j] += observation[i][j];
          }
        }
        break;
      }
      }
    }

    if (i < seq.length[index] - 1) {
      for (j = 0;j < nb_state;j++) {
        state_in[i][j] = D_INF;
        for (k = 0;k < nb_state;k++) {
          buff = cumul_transition[k][j] + forward1[i][k];
          if (buff > state_in[i][j]) {
            state_in[i][j] = buff;

#           ifdef MESSAGE
            input_state[i][j] = k;
#           endif

          }
        }
      }
    }

    for (j = 0;j < nb_output_process;j++) {
      switch (seq.type[j + 1]) {
      case INT_VALUE :
        pioutput[j]++;
        break;
      case REAL_VALUE :
        proutput[j]++;
        break;
      }
    }
  }

  // extraction of the log-likelihood for the most probable state sequence

# ifdef MESSAGE
  pstate = state_sequence + seq.length[index] - 1;
# endif

  state_seq_likelihood = D_INF;
  i = seq.length[index] - 1;
  for (j = 0;j < nb_state;j++) {
    if (forward1[i][j] > state_seq_likelihood) {
      state_seq_likelihood = forward1[i][j];

#     ifdef MESSAGE
      *pstate = j;
#     endif

    }
  }

  if (state_seq_likelihood != D_INF) {

#   ifdef MESSAGE
    i = seq.length[index] - 1;

    do {
      for (j = 0;j < optimal_forward_occupancy[i][*pstate] - 1;j++) {
        pstate--;
        *pstate = *(pstate + 1);
      }

      if (i >= optimal_forward_occupancy[i][*pstate]) {
        pstate--;
        *pstate = optimal_state[i][*(pstate + 1)];
        i -= optimal_forward_occupancy[i][*(pstate + 1)];
      }
      else {
        i -= optimal_forward_occupancy[i][*pstate];
      }
    }
    while (i >= 0);
#   endif

    // backward recurrence

    i = seq.length[index] - 1;
    for (j = 0;j < nb_state;j++) {
      backward1[i][j] = 0.;
      backward[i][j] = forward1[i][j];

      if (output == OUT_STATE) {
        backward_output[i][j] = backward[i][j];
      }
    }

    for (i = seq.length[index] - 2;i >= 0;i--) {
      for (j = 0;j < nb_state;j++) {
        switch (sojourn_type[j]) {

        // case semi-Markovian state

        case SEMI_MARKOVIAN : {
          occupancy = state_process->sojourn_time[j];
          obs_product = 0.;

          for (k = 1;k < MIN(seq.length[index] - i , occupancy->nb_value);k++) {
            if (observation[i + k][j] == D_INF) {
              break;
            }
            else {
              obs_product += observation[i + k][j];
            }

            if (k < seq.length[index] - i - 1) {
              occupancy_auxiliary[k] = backward1[i + k][j] + obs_product + occupancy->mass[k];
            }
            else {
              occupancy_auxiliary[k] = obs_product + occupancy->cumul[k - 1];
            }
          }

          auxiliary[j] = D_INF;
          for (m = k - 1;m >= 1;m--) {
            if (occupancy_auxiliary[m] > auxiliary[j]) {
              auxiliary[j] = occupancy_auxiliary[m];
            }

            // transformation of semi-Markovian log-likelihoods in Markovian log-likelihoods

            if ((auxiliary[j] != D_INF) && (state_in[i][j] != D_INF)) {
              buff = auxiliary[j] + state_in[i][j];
              if (buff > backward[i + m][j]) {
                backward[i + m][j] = buff;
              }
            }
          }
          break;
        }

        // case Markovian state

        case MARKOVIAN : {
          if ((backward1[i + 1][j] != D_INF) && (observation[i + 1][j] != D_INF)) {
            auxiliary[j] = backward1[i + 1][j] + observation[i + 1][j];
          }
          else {
            auxiliary[j] = D_INF;
          }
          break;
        }
        }
      }

      for (j = 0;j < nb_state;j++) {
        backward1[i][j] = D_INF;
        for (k = 0;k < nb_state;k++) {
          buff = auxiliary[k] + cumul_transition[j][k];
          if (buff > backward1[i][j]) {
            backward1[i][j] = buff;
          }
        }

        if ((backward1[i][j] != D_INF) && (forward1[i][j] != D_INF)) {
          backward[i][j] = backward1[i][j] + forward1[i][j];
        }
        else {
          backward[i][j] = D_INF;
        }
      }

      switch (output) {

      case IN_STATE : {
        for (j = 0;j < nb_state;j++) {
          switch (sojourn_type[j]) {

          // case semi-Markovian state

          case SEMI_MARKOVIAN : {
            if ((auxiliary[j] != D_INF) && (state_in[i][j] != D_INF)) {
              backward_output[i + 1][j] = auxiliary[j] + state_in[i][j];
            }
            else {
              backward_output[i + 1][j] = D_INF;
            }
            break;
          }

          // case Markovian state

          case MARKOVIAN : {
            backward_output[i + 1][j] = D_INF;

            if (auxiliary[j] != D_INF) {
              for (k = 0;k < nb_state;k++) {
                if (k != j) {
                  buff = cumul_transition[k][j] + forward1[i][k];
                  if (buff > backward_output[i + 1][j]) {
                    backward_output[i + 1][j] = buff;
                  }
                }
              }

              if (backward_output[i + 1][j] != D_INF) {
                backward_output[i + 1][j] += auxiliary[j];
              }
            }
            break;
          }
          }
        }
        break;
      }

      case OUT_STATE : {
        for (j = 0;j < nb_state;j++) {
          switch (sojourn_type[j]) {

          // case semi-Markovian state

          case SEMI_MARKOVIAN : {
            backward_output[i][j] = backward[i][j];
            break;
          }

          // case Markovian state

          case MARKOVIAN : {
            backward_output[i][j] = D_INF;

            if (forward1[i][j] != D_INF) {
              for (k = 0;k < nb_state;k++) {
                if (k != j) {
                  buff = auxiliary[k] + cumul_transition[j][k];
                  if (buff > backward_output[i][j]) {
                    backward_output[i][j] = buff;
                  }
                }
              }

              if (backward_output[i][j] != D_INF) {
                backward_output[i][j] += forward1[i][j];
              }
            }
            break;
          }
          }
        }
        break;
      }
      }
    }

    // particular case of staying in the initial state

    for (i = 0;i < nb_state;i++) {
      if ((sojourn_type[i] == SEMI_MARKOVIAN) && (cumul_initial[i] != D_INF)) {
        occupancy = state_process->sojourn_time[i];
        obs_product = 0.;

        for (j = 1;j < MIN(seq.length[index] + 1 , occupancy->nb_value);j++) {
          if (observation[j - 1][i] == D_INF) {
            break;
          }
          else {
            obs_product += observation[j - 1][i];
          }

          if (j < seq.length[index]) {
            switch (type) {
            case ORDINARY :
              occupancy_auxiliary[j] = backward1[j - 1][i] + obs_product + occupancy->mass[j];
              break;
            case EQUILIBRIUM :
              occupancy_auxiliary[j] = backward1[j - 1][i] + obs_product + forward[i]->mass[j];
              break;
            }
          }

          else {
            switch (type) {
            case ORDINARY :
              occupancy_auxiliary[j] = obs_product + occupancy->cumul[j - 1];
              break;
            case EQUILIBRIUM :
              occupancy_auxiliary[j] = obs_product + forward[i]->cumul[j - 1];
              break;
            }
          }
        }

        auxiliary[i] = D_INF;
        for (k = j - 1;k >= 1;k--) {
          if (occupancy_auxiliary[k] > auxiliary[i]) {
            auxiliary[i] = occupancy_auxiliary[k];
          }

          // transformation of semi-Markovian log-likelihoods in Markovian log-likelihoods

          if (auxiliary[i] != D_INF) {
            buff = auxiliary[i] + cumul_initial[i];
            if (buff > backward[k - 1][i]) {
              backward[k - 1][i] = buff;
            }
          }
        }
      }

      if (output == IN_STATE) {
        backward_output[0][i] = backward[0][i];
      }
    }

    // restoration of the most probable state sequence

    pstate = seq.int_sequence[index][0];

    for (i = 0;i < seq.length[index];i++) {
      backward_max = D_INF;
      for (j = 0;j < nb_state;j++) {
        if (backward[i][j] > backward_max) {
          backward_max = backward[i][j];
          *pstate = j;
        }
      }

#     ifdef MESSAGE
      if (*pstate != state_sequence[i]) {
        cout << "\nERROR: " << i << " | " << *pstate << " " << state_sequence[i] << endl;
      }
#     endif

      pstate++;
    }

    //  normalization

    for (i = 0;i < seq.length[index];i++) {
      for (j = 0;j < nb_state;j++) {
        if (backward_output[i][j] != D_INF) {
          backward_output[i][j] = exp(backward_output[i][j] - seq_likelihood);
//          backward_output[i][j] = exp(backward_output[i][j] - state_seq_likelihood);
        }
        else {
          backward_output[i][j] = 0.;
        }
      }
    }

    switch (format) {

    case ASCII : {
      switch (output) {
      case SSTATE :
        *os << "\n" << SEQ_label[SEQL_MAX_POSTERIOR_STATE_PROBABILITY] << "\n\n";
        break;
      case IN_STATE :
        *os << "\n" << SEQ_label[SEQL_MAX_POSTERIOR_IN_STATE_PROBABILITY] << "\n\n";
        break;
      case OUT_STATE :
        *os << "\n" << SEQ_label[SEQL_MAX_POSTERIOR_OUT_STATE_PROBABILITY] << "\n\n";
        break;
      }

      seq.profile_ascii_print(*os , index , nb_state , backward_output ,
                              STAT_label[STATL_STATE]);

      *os << "\n" << SEQ_label[SEQL_STATE_SEQUENCE_LIKELIHOOD] << ": " << state_seq_likelihood
          << "   (" << exp(state_seq_likelihood - seq_likelihood) << ")" << endl;
      break;
    }

    case SPREADSHEET : {
      switch (output) {
      case SSTATE :
        *os << "\n" << SEQ_label[SEQL_MAX_POSTERIOR_STATE_PROBABILITY] << "\n\n";
        break;
      case IN_STATE :
        *os << "\n" << SEQ_label[SEQL_MAX_POSTERIOR_IN_STATE_PROBABILITY] << "\n\n";
        break;
      case OUT_STATE :
        *os << "\n" << SEQ_label[SEQL_MAX_POSTERIOR_OUT_STATE_PROBABILITY] << "\n\n";
        break;
      }

      seq.profile_spreadsheet_print(*os , index , nb_state , backward_output ,
                                    STAT_label[STATL_STATE]);

      *os << "\n" << SEQ_label[SEQL_STATE_SEQUENCE_LIKELIHOOD] << "\t" << state_seq_likelihood
          << "\t" << exp(state_seq_likelihood - seq_likelihood) << endl;
      break;
    }

    case GNUPLOT : {
      seq.profile_plot_print(*os , index , nb_state , backward_output);
      break;
    }

    case PLOT : {
      seq.profile_plotable_write(*plot , index , nb_state , backward_output);
      break;
    }
    }

#   ifdef DEBUG
    if (format != GNUPLOT) {
      double ambiguity = 0.;

      pstate = seq.int_sequence[index][0];
//      if (output == SSTATE) {
      for (i = 0;i < seq.length[index];i++) {
        for (j = 0;j < nb_state;j++) {
          if (j != *pstate) {
            ambiguity += backward_output[i][j];
          }
        }
        pstate++;
      }
      ambiguity *= exp(seq_likelihood - state_seq_likelihood);
/*      }

      else {
        for (i = 0;i < seq.length[index];i++) {
          for (j = 0;j < nb_state;j++) {
            if ((backward[i][j] != D_INF) && (j != *pstate)) {
              ambiguity += exp(backward[i][j] - state_seq_likelihood);
            }
          }
          pstate++;
        }
      } */

      switch (format) {
      case ASCII :
        *os << "\n" << SEQ_label[SEQL_AMBIGUITY] << ": " << ambiguity
            << " (" << ambiguity / seq.length[index] << ")" << endl;
        break;
      case SPREADSHEET :
        *os << "\n" << SEQ_label[SEQL_AMBIGUITY] << "\t" << ambiguity
            << "\t" << ambiguity / seq.length[index] << "\t" << endl;
        break;
      }
    }
#   endif

  }

  for (i = 0;i < seq.length[index];i++) {
    delete [] observation[i];
  }
  delete [] observation;

  for (i = 0;i < seq.length[index];i++) {
    delete [] forward1[i];
  }
  delete [] forward1;

  for (i = 0;i < seq.length[index] - 1;i++) {
    delete [] state_in[i];
  }
  delete [] state_in;

  for (i = 0;i < seq.length[index];i++) {
    delete [] backward[i];
  }
  delete [] backward;

  for (i = 0;i < seq.length[index];i++) {
    delete [] backward1[i];
  }
  delete [] backward1;

  delete [] auxiliary;
  delete [] occupancy_auxiliary;

  if (output != SSTATE) {
    for (i = 0;i < seq.length[index];i++) {
      delete [] backward_output[i];
    }
    delete [] backward_output;
  }

  delete [] pioutput;
  delete [] proutput;

# ifdef MESSAGE
  for (i = 0;i < seq.length[index] - 1;i++) {
    delete [] input_state[i];
  }
  delete [] input_state;

  for (i = 0;i < seq.length[index];i++) {
    delete [] optimal_state[i];
  }
  delete [] optimal_state;

  for (i = 0;i < seq.length[index];i++) {
    delete [] optimal_forward_occupancy[i];
  }
  delete [] optimal_forward_occupancy;

  delete [] state_sequence;
# endif

  return state_seq_likelihood;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of state and entropy profiles using the forward-backward algorithm,
 *         of state profiles using the Viterbi forward-backward algorithm,
 *         computation of the N most probable state sequences using the generalized Viterbi algorithm or
 *         simulation of state sequences using the forward-backward algorithm for sampling and
 *         writing of the results.
 *
 *  \param[in] error             reference on a StatError object,
 *  \param[in] os                stream,
 *  \param[in] iseq              reference on a MarkovianSequences object,
 *  \param[in] identifier        sequence identifier,
 *  \param[in] output            output type,
 *  \param[in] format            format (ASCII/SPREADSHEET),
 *  \param[in] state_sequence    method for computing the state sequences (GENERALIZED_VITERBI/FORWARD_BACKWARD_SAMPLING),
 *  \param[in] nb_state_sequence number of state sequences.
 *
 *  \return                      error status.
 */
/*--------------------------------------------------------------*/

bool HiddenSemiMarkov::state_profile_write(StatError &error , ostream &os ,
                                           const MarkovianSequences &iseq , int identifier ,
                                           state_profile output , output_format format ,
                                           latent_structure_algorithm state_sequence ,
                                           int nb_state_sequence) const

{
  bool status = true;
  int i;
  int offset = I_DEFAULT , nb_value , index = I_DEFAULT;
  double seq_likelihood , max_marginal_entropy , entropy;
  HiddenSemiMarkov *hsmarkov1 , *hsmarkov2;
  SemiMarkovData *seq;


  error.init();

  for (i = 0;i < iseq.nb_variable;i++) {
    if ((iseq.type[i] != INT_VALUE) && (iseq.type[i] != REAL_VALUE) && (iseq.type[i] != STATE)) {
      status = false;
      ostringstream error_message , correction_message;
      error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                    << STAT_error[STATR_VARIABLE_TYPE];
      correction_message << STAT_variable_word[INT_VALUE] << " or "
                         << STAT_variable_word[REAL_VALUE];
      error.correction_update((error_message.str()).c_str() , (correction_message.str()).c_str());
    }
  }

  if (nb_output_process == iseq.nb_variable) {
    offset = 0;
  }
  else if ((iseq.type[0] == STATE) && (nb_output_process + 1 == iseq.nb_variable)) {
    offset = 1;
  }
  else {
    status = false;
    error.update(STAT_error[STATR_NB_OUTPUT_PROCESS]);
  }

  if (offset != I_DEFAULT) {
    for (i = 0;i < nb_output_process;i++) {
      if ((categorical_process[i]) || (discrete_parametric_process[i])) {
        if (iseq.type[i + offset] == REAL_VALUE) {
          status = false;
          ostringstream error_message;
          error_message << STAT_label[STATL_VARIABLE] << " " << i + offset + 1 << ": "
                        << STAT_error[STATR_VARIABLE_TYPE];
          error.correction_update((error_message.str()).c_str() , STAT_variable_word[INT_VALUE]);
        }

        else {
          if (iseq.min_value[i + offset] < 0) {
            status = false;
            ostringstream error_message;
            error_message << STAT_label[STATL_VARIABLE] << " " << i + offset + 1 << ": "
                          << STAT_error[STATR_POSITIVE_MIN_VALUE];
            error.update((error_message.str()).c_str());
          }

          if (!(iseq.marginal_distribution[i + offset])) {
            status = false;
            ostringstream error_message;
            error_message << STAT_label[STATL_VARIABLE] << " " << i + offset + 1 << ": "
                          << STAT_error[STATR_MARGINAL_FREQUENCY_DISTRIBUTION];
            error.update((error_message.str()).c_str());
          }

          else {
            if (categorical_process[i]) {
              nb_value = categorical_process[i]->nb_value;
            }
            else {
              nb_value = discrete_parametric_process[i]->nb_value;
            }

            if (nb_value < iseq.marginal_distribution[i + offset]->nb_value) {
              status = false;
              ostringstream error_message;
              error_message << STAT_label[STATL_OUTPUT_PROCESS] << " " << i + 1 << ": "
                            << STAT_error[STATR_NB_OUTPUT];
              error.update((error_message.str()).c_str());
            }
          }
        }
      }
    }
  }

  if (identifier != I_DEFAULT) {
    for (i = 0;i < iseq.nb_sequence;i++) {
      if (identifier == iseq.identifier[i]) {
        index = i;
        break;
      }
    }

    if (i == iseq.nb_sequence) {
      status = false;
      error.update(SEQ_error[SEQR_SEQUENCE_IDENTIFIER]);
    }
  }

  if (nb_state_sequence < 2) {
    status = false;
    error.update(SEQ_error[SEQR_NB_STATE_SEQUENCE]);
  }

  if (status) {
    if (nb_output_process == iseq.nb_variable) {
      seq = new SemiMarkovData(iseq);
    }
    else {
      seq = new SemiMarkovData(iseq , SEQUENCE_COPY , (type == EQUILIBRIUM ? true : false));
    }

    hsmarkov1 = new HiddenSemiMarkov(*this , false);

    hsmarkov2 = new HiddenSemiMarkov(*this , false);
    hsmarkov2->create_cumul();
    hsmarkov2->log_computation();

    for (i = 0;i < seq->nb_sequence;i++) {
      if ((index == I_DEFAULT) || (index == i)) {
        seq_likelihood = hsmarkov1->forward_backward(*seq , i , &os , NULL , output , format ,
                                                     max_marginal_entropy , entropy);

        if (seq_likelihood == D_INF) {
          status = false;

          if (index == I_DEFAULT) {
            ostringstream error_message;
            error_message << SEQ_label[SEQL_SEQUENCE] << " " << i + 1 << " "
                          << SEQ_error[SEQR_INCOMPATIBLE_MODEL];
            error.update((error_message.str()).c_str());
          }
          else {
            error.update(SEQ_error[SEQR_SEQUENCE_INCOMPATIBLE_MODEL]);
          }
        }

        else {
          hsmarkov2->viterbi_forward_backward(*seq , i , &os , NULL , output , format ,
                                              seq_likelihood);

          switch (state_sequence) {
          case GENERALIZED_VITERBI :
            hsmarkov2->generalized_viterbi(*seq , i , os , seq_likelihood , format ,
                                           nb_state_sequence);
            break;
          case FORWARD_BACKWARD_SAMPLING :
            hsmarkov1->forward_backward_sampling(*seq , i , os , format ,
                                                 nb_state_sequence);
            break;
          }
        }
      }
    }

    delete seq;

    delete hsmarkov1;
    delete hsmarkov2;
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of state and entropy profiles using the forward-backward algorithm,
 *         of state profiles using the Viterbi forward-backward algorithm,
 *         computation of the N most probable state sequences using the generalized Viterbi algorithm or
 *         simulation of state sequences using the forward-backward algorithm for sampling and
 *         displaying the results.
 *
 *  \param[in] error             reference on a StatError object,
 *  \param[in] iseq              reference on a MarkovianSequences object,
 *  \param[in] identifier        sequence identifier,
 *  \param[in] output            output type,
 *  \param[in] state_sequence    method for computing the state sequences (GENERALIZED_VITERBI/FORWARD_BACKWARD_SAMPLING),
 *  \param[in] nb_state_sequence number of state sequences.
 *
 *  \return                      error status.
 */
/*--------------------------------------------------------------*/

bool HiddenSemiMarkov::state_profile_ascii_write(StatError &error , const MarkovianSequences &iseq ,
                                                 int identifier , state_profile output ,
                                                 latent_structure_algorithm state_sequence ,
                                                 int nb_state_sequence) const

{
  return state_profile_write(error , cout , iseq , identifier ,
                             output , ASCII , state_sequence , nb_state_sequence);
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of state and entropy profiles using the forward-backward algorithm,
 *         of state profiles using the Viterbi forward-backward algorithm,
 *         computation of the N most probable state sequences using the generalized Viterbi algorithm or
 *         simulation of state sequences using the forward-backward algorithm for sampling and
 *         writing of the results in a file.
 *
 *  \param[in] error             reference on a StatError object,
 *  \param[in] path              file path,
 *  \param[in] iseq              reference on a MarkovianSequences object,
 *  \param[in] identifier        sequence identifier,
 *  \param[in] output            output type,
 *  \param[in] format            file format (ASCII/SPREADSHEET),
 *  \param[in] state_sequence    method for computing the state sequences (GENERALIZED_VITERBI/FORWARD_BACKWARD_SAMPLING),
 *  \param[in] nb_state_sequence number of state sequences.
 *
 *  \return                      error status.
 */
/*--------------------------------------------------------------*/

bool HiddenSemiMarkov::state_profile_write(StatError &error , const string path ,
                                           const MarkovianSequences &iseq , int identifier ,
                                           state_profile output , output_format format ,
                                           latent_structure_algorithm state_sequence ,
                                           int nb_state_sequence) const

{
  bool status = true;
  ofstream out_file(path.c_str());


  error.init();

  if (!out_file) {
    status = false;
    error.update(STAT_error[STATR_FILE_NAME]);
  }
  else {
    status = state_profile_write(error , out_file , iseq , identifier ,
                                 output , format , state_sequence , nb_state_sequence);
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of state and entropy profiles using the forward-backward algorithm,
 *         of state profiles using the Viterbi forward-backward algorithm,
 *         computation of the N most probable state sequences using the generalized Viterbi algorithm or
 *         simulation of state sequences using the forward-backward algorithm for sampling and
 *         displaying of the results.
 *
 *  \param[in] error             reference on a StatError object,
 *  \param[in] identifier        sequence identifier,
 *  \param[in] output            output type,
 *  \param[in] state_sequence    method for computing the state sequences (GENERALIZED_VITERBI/FORWARD_BACKWARD_SAMPLING),
 *  \param[in] nb_state_sequence number of state sequences.
 *
 *  \return                      error status.
 */
/*--------------------------------------------------------------*/

bool HiddenSemiMarkov::state_profile_ascii_write(StatError &error , int identifier , state_profile output , 
                                                 latent_structure_algorithm state_sequence ,
                                                 int nb_state_sequence) const

{
  bool status;


  error.init();

  if (!semi_markov_data) {
    status = false;
    error.update(STAT_error[STATR_NO_DATA]);
  }
  else {
    status = state_profile_write(error , cout , *semi_markov_data , identifier ,
                                 output , ASCII , state_sequence , nb_state_sequence);
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of state and entropy profiles using the forward-backward algorithm,
 *         of state profiles using the Viterbi forward-backward algorithm,
 *         computation of the N most probable state sequences using the generalized Viterbi algorithm or
 *         simulation of state sequences using the forward-backward algorithm for sampling and
 *         writing of the results in a file.
 *
 *  \param[in] error             reference on a StatError object,
 *  \param[in] path              file path,
 *  \param[in] identifier        sequence identifier,
 *  \param[in] output            output type,
 *  \param[in] format            file format (ASCII/SPREADSHEET),
 *  \param[in] state_sequence    method for computing the state sequences (GENERALIZED_VITERBI/FORWARD_BACKWARD_SAMPLING),
 *  \param[in] nb_state_sequence number of state sequences.
 *
 *  \return                      error status.
 */
/*--------------------------------------------------------------*/

bool HiddenSemiMarkov::state_profile_write(StatError &error , const string path , int identifier ,
                                           state_profile output , output_format format ,
                                           latent_structure_algorithm state_sequence ,
                                           int nb_state_sequence) const

{
  bool status = true;
  ofstream out_file(path.c_str());


  error.init();

  if (!out_file) {
    status = false;
    error.update(STAT_error[STATR_FILE_NAME]);
  }
  if (!semi_markov_data) {
    status = false;
    error.update(STAT_error[STATR_NO_DATA]);
  }

  if (status) {
    status = state_profile_write(error , out_file , *semi_markov_data , identifier ,
                                 output , format , state_sequence , nb_state_sequence);
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of state and entropy profiles using the forward-backward algorithm,
 *         of state profiles using the Viterbi forward-backward algorithm and
 *         plot of the results at the Gnuplot format.
 *
 *  \param[in] error      reference on a StatError object,
 *  \param[in] prefix     file prefix,
 *  \param[in] iseq       reference on a MarkovianSequences object,
 *  \param[in] identifier sequence identifier,
 *  \param[in] output     output type,
 *  \param[in] title      figure title.
 *
 *  \return               error status.
 */
/*--------------------------------------------------------------*/

bool HiddenSemiMarkov::state_profile_plot_write(StatError &error , const char *prefix ,
                                                const MarkovianSequences &iseq , int identifier ,
                                                state_profile output , const char *title) const

{
  bool status = true;
  int i , j;
  int offset = I_DEFAULT , nb_value , index;
  double seq_likelihood , max_marginal_entropy , entropy , state_seq_likelihood;
  HiddenSemiMarkov *hsmarkov;
  SemiMarkovData *seq;
  ostringstream data_file_name[2];
  ofstream *out_data_file;


  error.init();

  for (i = 0;i < iseq.nb_variable;i++) {
    if ((iseq.type[i] != INT_VALUE) && (iseq.type[i] != REAL_VALUE) && (iseq.type[i] != STATE)) {
      status = false;
      ostringstream error_message , correction_message;
      error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                    << STAT_error[STATR_VARIABLE_TYPE];
      correction_message << STAT_variable_word[INT_VALUE] << " or "
                         << STAT_variable_word[REAL_VALUE];
      error.correction_update((error_message.str()).c_str() , (correction_message.str()).c_str());
    }
  }

  if (nb_output_process == iseq.nb_variable) {
    offset = 0;
  }
  else if ((iseq.type[0] == STATE) && (nb_output_process + 1 == iseq.nb_variable)) {
    offset = 1;
  }
  else {
    status = false;
    error.update(STAT_error[STATR_NB_OUTPUT_PROCESS]);
  }

  if (offset != I_DEFAULT) {
    for (i = 0;i < nb_output_process;i++) {
      if ((categorical_process[i]) || (discrete_parametric_process[i])) {
        if (iseq.type[i + offset] == REAL_VALUE) {
          status = false;
          ostringstream error_message;
          error_message << STAT_label[STATL_VARIABLE] << " " << i + offset + 1 << ": "
                        << STAT_error[STATR_VARIABLE_TYPE];
          error.correction_update((error_message.str()).c_str() , STAT_variable_word[INT_VALUE]);
        }

        else {
          if (iseq.min_value[i + offset] < 0) {
            status = false;
            ostringstream error_message;
            error_message << STAT_label[STATL_VARIABLE] << " " << i + offset + 1 << ": "
                          << STAT_error[STATR_POSITIVE_MIN_VALUE];
            error.update((error_message.str()).c_str());
          }

          if (!(iseq.marginal_distribution[i + offset])) {
            status = false;
            ostringstream error_message;
            error_message << STAT_label[STATL_VARIABLE] << " " << i + offset + 1 << ": "
                          << STAT_error[STATR_MARGINAL_FREQUENCY_DISTRIBUTION];
            error.update((error_message.str()).c_str());
          }

          else {
            if (categorical_process[i]) {
              nb_value = categorical_process[i]->nb_value;
            }
            else {
              nb_value = discrete_parametric_process[i]->nb_value;
            }

            if (nb_value < iseq.marginal_distribution[i + offset]->nb_value) {
              status = false;
              ostringstream error_message;
              error_message << STAT_label[STATL_OUTPUT_PROCESS] << " " << i + 1 << ": "
                            << STAT_error[STATR_NB_OUTPUT];
              error.update((error_message.str()).c_str());
            }
          }
        }
      }
    }
  }

  for (i = 0;i < iseq.nb_sequence;i++) {
    if (identifier == iseq.identifier[i]) {
      index = i;
      break;
    }
  }

  if (i == iseq.nb_sequence) {
    status = false;
    error.update(SEQ_error[SEQR_SEQUENCE_IDENTIFIER]);
  }

  if (status) {

    // writing of the date files

    data_file_name[0] << prefix << 0 << ".dat";
    out_data_file = new ofstream((data_file_name[0].str()).c_str());

    if (!out_data_file) {
      status = false;
      error.update(STAT_error[STATR_FILE_PREFIX]);
    }

    else {
      if (iseq.type[0] != STATE) {
        seq = new SemiMarkovData(iseq);
      }
      else {
        seq = new SemiMarkovData(iseq , SEQUENCE_COPY , (type == EQUILIBRIUM ? true : false));
      }

      hsmarkov = new HiddenSemiMarkov(*this , false);

      seq_likelihood = hsmarkov->forward_backward(*seq , index , out_data_file , NULL , output ,
                                                  GNUPLOT , max_marginal_entropy , entropy);
      out_data_file->close();
      delete out_data_file;

      if (seq_likelihood == D_INF) {
        status = false;
        error.update(SEQ_error[SEQR_SEQUENCE_INCOMPATIBLE_MODEL]);
      }

      else {
        data_file_name[1] << prefix << 1 << ".dat";
        out_data_file = new ofstream((data_file_name[1].str()).c_str());

        hsmarkov->create_cumul();
        hsmarkov->log_computation();
        state_seq_likelihood = hsmarkov->viterbi_forward_backward(*seq , index , out_data_file , NULL ,
                                                                  output , GNUPLOT , seq_likelihood);
        out_data_file->close();
        delete out_data_file;

        // writing of the script files

        for (i = 0;i < 2;i++) {
          ostringstream file_name[2];

          switch (i) {
          case 0 :
            file_name[0] << prefix << ".plot";
            break;
          case 1 :
            file_name[0] << prefix << ".print";
            break;
          }

          ofstream out_file((file_name[0].str()).c_str());

          if (i == 1) {
            out_file << "set terminal postscript" << endl;
            file_name[1] << label(prefix) << ".ps";
            out_file << "set output \"" << file_name[1].str() << "\"\n\n";
          }

          out_file << "set border 15 lw 0\n" << "set tics out\n" << "set xtics nomirror\n"
                   << "set title \"";
          if (title) {
            out_file << title << " - ";
          }
          switch (output) {
          case SSTATE :
            out_file << SEQ_label[SEQL_MAX_POSTERIOR_STATE_PROBABILITY] << "\"\n\n";
            break;
          case IN_STATE :
            out_file << SEQ_label[SEQL_MAX_POSTERIOR_IN_STATE_PROBABILITY] << "\"\n\n";
            break;
          case OUT_STATE :
            out_file << SEQ_label[SEQL_MAX_POSTERIOR_OUT_STATE_PROBABILITY] << "\"\n\n";
            break;
          }

          if (seq->index_parameter) {
            if (seq->index_parameter[index][seq->length[index] - 1] - seq->index_parameter[index][0] < TIC_THRESHOLD) {
              out_file << "set xtics 0,1" << endl;
            }

            out_file << "plot [" << seq->index_parameter[index][0] << ":"
                     << seq->index_parameter[index][seq->length[index] - 1] << "] [0:"
                     << exp(state_seq_likelihood - seq_likelihood) << "] ";
            for (j = 0;j < nb_state;j++) {
              out_file << "\"" << label((data_file_name[1].str()).c_str()) << "\" using "
                       << 1 << " : " << j + 2 << " title \"" << STAT_label[STATL_STATE] << " "
                       << j << "\" with linespoints";
              if (j < nb_state - 1) {
                out_file << ",\\";
              }
              out_file << endl;
            }

            if (i == 0) {
              out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
            }
            out_file << endl;

            out_file << "set title \"";
            if (title) {
              out_file << title << " - ";
            }
            switch (output) {
            case SSTATE :
              out_file << SEQ_label[SEQL_POSTERIOR_STATE_PROBABILITY] << "\"\n\n";
              break;
            case IN_STATE :
              out_file << SEQ_label[SEQL_POSTERIOR_IN_STATE_PROBABILITY] << "\"\n\n";
              break;
            case OUT_STATE :
              out_file << SEQ_label[SEQL_POSTERIOR_OUT_STATE_PROBABILITY] << "\"\n\n";
              break;
            }

            out_file << "plot [" << seq->index_parameter[index][0] << ":"
                     << seq->index_parameter[index][seq->length[index] - 1] << "] [0:1] ";
            for (j = 0;j < nb_state;j++) {
              out_file << "\"" << label((data_file_name[0].str()).c_str()) << "\" using "
                       << 1 << " : " << j + 2 << " title \"" << STAT_label[STATL_STATE] << " "
                       << j << "\" with linespoints";
              if (j < nb_state - 1) {
                out_file << ",\\";
              }
              out_file << endl;
            }

            if (i == 0) {
              out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
            }
            out_file << endl;

            out_file << "set title";
            if (title) {
              out_file << " \"" << title << "\"";
            }
            out_file << "\n\n";

            out_file << "plot [" << seq->index_parameter[index][0] << ":"
                     << seq->index_parameter[index][seq->length[index] - 1] << "] [0:"
                     << max_marginal_entropy << "] "
                     << "\"" << label((data_file_name[0].str()).c_str()) << "\" using "
                     << 1 << " : " << nb_state + 2 << " title \"" << SEQ_label[SEQL_CONDITIONAL_ENTROPY]
                     << "\" with linespoints,\\" << endl;
            out_file << "\"" << label((data_file_name[0].str()).c_str()) << "\" using "
                     << 1 << " : " << nb_state + 3 << " title \"" << SEQ_label[SEQL_MARGINAL_ENTROPY]
                     << "\" with linespoints" << endl;

            if (i == 0) {
              out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
            }
            out_file << endl;

            out_file << "set title";
            if (title) {
              out_file << " \"" << title << "\"";
            }
            out_file << "\n\n";

            out_file << "plot [" << seq->index_parameter[index][0] << ":"
                     << seq->index_parameter[index][seq->length[index] - 1] << "] [0:" << entropy << "] "
                     << "\"" << label((data_file_name[0].str()).c_str()) << "\" using "
                     << 1 << " : " << nb_state + 4 << " title \""
                     << SEQ_label[SEQL_PARTIAL_STATE_SEQUENCE_ENTROPY] << "\" with linespoints" << endl;

            if (seq->index_parameter[index][seq->length[index] - 1] - seq->index_parameter[index][0] < TIC_THRESHOLD) {
              out_file << "set xtics autofreq" << endl;
            }
          }

          else {
            if (seq->length[index] - 1 < TIC_THRESHOLD) {
              out_file << "set xtics 0,1" << endl;
            }

            out_file << "plot [0:" << seq->length[index] - 1 << "] [0:"
                     << exp(state_seq_likelihood - seq_likelihood) << "] ";
            for (j = 0;j < nb_state;j++) {
              out_file << "\"" << label((data_file_name[1].str()).c_str()) << "\" using "
//                       << j + 1 << " title \"" << STAT_label[STATL_STATE] << " "
                       << 1 << " : " << j + 2 << " title \"" << STAT_label[STATL_STATE] << " "
                       << j << "\" with linespoints";
              if (j < nb_state - 1) {
                out_file << ",\\";
              }
              out_file << endl;
            }

            if (i == 0) {
              out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
            }
            out_file << endl;

            out_file << "set title \"";
            if (title) {
              out_file << title << " - ";
            }
            switch (output) {
            case SSTATE :
              out_file << SEQ_label[SEQL_POSTERIOR_STATE_PROBABILITY] << "\"\n\n";
              break;
            case IN_STATE :
              out_file << SEQ_label[SEQL_POSTERIOR_IN_STATE_PROBABILITY] << "\"\n\n";
              break;
            case OUT_STATE :
              out_file << SEQ_label[SEQL_POSTERIOR_OUT_STATE_PROBABILITY] << "\"\n\n";
              break;
            }

            out_file << "plot [0:" << seq->length[index] - 1 << "] [0:1] ";
            for (j = 0;j < nb_state;j++) {
              out_file << "\"" << label((data_file_name[0].str()).c_str()) << "\" using "
                       << j + 1 << " title \"" << STAT_label[STATL_STATE] << " "
                       << j << "\" with linespoints";
              if (j < nb_state - 1) {
                out_file << ",\\";
              }
              out_file << endl;
            }

            if (i == 0) {
              out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
            }
            out_file << endl;

            out_file << "set title";
            if (title) {
              out_file << " \"" << title << "\"";
            }
            out_file << "\n\n";

            out_file << "plot [0:" << seq->length[index] - 1 << "] [0:" << max_marginal_entropy << "] "
                     << "\"" << label((data_file_name[0].str()).c_str()) << "\" using "
                     << nb_state + 1 << " title \"" << SEQ_label[SEQL_CONDITIONAL_ENTROPY]
                     << "\" with linespoints,\\" << endl;
            out_file << "\"" << label((data_file_name[0].str()).c_str()) << "\" using "
                     << nb_state + 2 << " title \"" << SEQ_label[SEQL_MARGINAL_ENTROPY]
                     << "\" with linespoints" << endl;

            if (i == 0) {
              out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
            }
            out_file << endl;

            out_file << "set title";
            if (title) {
              out_file << " \"" << title << "\"";
            }
            out_file << "\n\n";

            out_file << "plot [0:" << seq->length[index] - 1 << "] [0:" << entropy << "] "
                     << "\"" << label((data_file_name[0].str()).c_str()) << "\" using "
                     << nb_state + 3 << " title \"" << SEQ_label[SEQL_PARTIAL_STATE_SEQUENCE_ENTROPY]
                     << "\" with linespoints" << endl;

            if (seq->length[index] - 1 < TIC_THRESHOLD) {
              out_file << "set xtics autofreq" << endl;
            }
          }

          if (i == 1) {
            out_file << "\nset terminal x11" << endl;
          }

          out_file << "\npause 0 \"" << STAT_label[STATL_END] << "\"" << endl;
        }
      }

      delete seq;
      delete hsmarkov;
    }
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of state and entropy profiles using the forward-backward algorithm,
 *         of state profiles using the Viterbi forward-backward algorithm and
 *         plot of the results at the Gnuplot format.
 *
 *  \param[in] error      reference on a StatError object,
 *  \param[in] prefix     file prefix,
 *  \param[in] identifier sequence identifier,
 *  \param[in] output     output type,
 *  \param[in] title      figure title.
 *
 *  \return               error status.
 */
/*--------------------------------------------------------------*/

bool HiddenSemiMarkov::state_profile_plot_write(StatError &error , const char *prefix , int identifier ,
                                                state_profile output , const char *title) const

{
  bool status;


  error.init();

  if (!semi_markov_data) {
    status = false;
    error.update(STAT_error[STATR_NO_DATA]);
  }
  else {
    status = state_profile_plot_write(error , prefix , *semi_markov_data , identifier ,
                                      output , title);
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of state and entropy profiles using the forward-backward algorithm,
 *         of state profiles using the Viterbi forward-backward algorithm and
 *         plot of the profiles.
 *
 *  \param[in] error      reference on a StatError object,
 *  \param[in] iseq       reference on a MarkovianSequences object,
 *  \param[in] identifier sequence identifier,
 *  \param[in] output     output type.
 *
 *  \return               MultiPlotSet object.
 */
/*--------------------------------------------------------------*/

MultiPlotSet* HiddenSemiMarkov::state_profile_plotable_write(StatError &error ,
                                                             const MarkovianSequences &iseq ,
                                                             int identifier , state_profile output) const

{
  bool status = true;
  int i;
  int offset = I_DEFAULT , nb_value , index;
  double seq_likelihood , max_marginal_entropy , entropy , state_seq_likelihood;
  HiddenSemiMarkov *hsmarkov;
  SemiMarkovData *seq;
  ostringstream legend;
  MultiPlotSet *plot_set;


  plot_set = NULL;
  error.init();

  for (i = 0;i < iseq.nb_variable;i++) {
    if ((iseq.type[i] != INT_VALUE) && (iseq.type[i] != REAL_VALUE) && (iseq.type[i] != STATE)) {
      status = false;
      ostringstream error_message , correction_message;
      error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                    << STAT_error[STATR_VARIABLE_TYPE];
      correction_message << STAT_variable_word[INT_VALUE] << " or "
                         << STAT_variable_word[REAL_VALUE];
      error.correction_update((error_message.str()).c_str() , (correction_message.str()).c_str());
    }
  }

  if (nb_output_process == iseq.nb_variable) {
    offset = 0;
  }
  else if ((iseq.type[0] == STATE) && (nb_output_process + 1 == iseq.nb_variable)) {
    offset = 1;
  }
  else {
    status = false;
    error.update(STAT_error[STATR_NB_OUTPUT_PROCESS]);
  }

  if (offset != I_DEFAULT) {
    for (i = 0;i < nb_output_process;i++) {
      if ((categorical_process[i]) || (discrete_parametric_process[i])) {
        if (iseq.type[i + offset] == REAL_VALUE) {
          status = false;
          ostringstream error_message;
          error_message << STAT_label[STATL_VARIABLE] << " " << i + offset + 1 << ": "
                        << STAT_error[STATR_VARIABLE_TYPE];
          error.correction_update((error_message.str()).c_str() , STAT_variable_word[INT_VALUE]);
        }

        else {
          if (iseq.min_value[i + offset] < 0) {
            status = false;
            ostringstream error_message;
            error_message << STAT_label[STATL_VARIABLE] << " " << i + offset + 1 << ": "
                          << STAT_error[STATR_POSITIVE_MIN_VALUE];
            error.update((error_message.str()).c_str());
          }

          if (!(iseq.marginal_distribution[i + offset])) {
            status = false;
            ostringstream error_message;
            error_message << STAT_label[STATL_VARIABLE] << " " << i + offset + 1 << ": "
                          << STAT_error[STATR_MARGINAL_FREQUENCY_DISTRIBUTION];
            error.update((error_message.str()).c_str());
          }

          else {
            if (categorical_process[i]) {
              nb_value = categorical_process[i]->nb_value;
            }
            else {
              nb_value = discrete_parametric_process[i]->nb_value;
            }

            if (nb_value < iseq.marginal_distribution[i + offset]->nb_value) {
              status = false;
              ostringstream error_message;
              error_message << STAT_label[STATL_OUTPUT_PROCESS] << " " << i + 1 << ": "
                            << STAT_error[STATR_NB_OUTPUT];
              error.update((error_message.str()).c_str());
            }
          }
        }
      }
    }
  }

  for (i = 0;i < iseq.nb_sequence;i++) {
    if (identifier == iseq.identifier[i]) {
      index = i;
      break;
    }
  }

  if (i == iseq.nb_sequence) {
    status = false;
    error.update(SEQ_error[SEQR_SEQUENCE_IDENTIFIER]);
  }

  if (status) {
    plot_set = new MultiPlotSet(4);

    MultiPlotSet &plot = *plot_set;

    plot.border = "15 lw 0";

    if (iseq.type[0] != STATE) {
      seq = new SemiMarkovData(iseq);
    }
    else {
      seq = new SemiMarkovData(iseq , SEQUENCE_COPY , (type == EQUILIBRIUM ? true : false));
    }

    hsmarkov = new HiddenSemiMarkov(*this , false);

    seq_likelihood = hsmarkov->forward_backward(*seq , index , NULL , plot_set , output ,
                                                PLOT , max_marginal_entropy , entropy);

    if (seq_likelihood == D_INF) {
      delete plot_set;
      plot_set = NULL;
      error.update(SEQ_error[SEQR_SEQUENCE_INCOMPATIBLE_MODEL]);
    }

    else {
      hsmarkov->create_cumul();
      hsmarkov->log_computation();
      state_seq_likelihood = hsmarkov->viterbi_forward_backward(*seq , index , NULL , &plot[0] ,
                                                                output , PLOT , seq_likelihood);

      // maximum posterior probabilities

      switch (output) {
      case SSTATE :
        plot[0].title = SEQ_label[SEQL_MAX_POSTERIOR_STATE_PROBABILITY];
        break;
      case IN_STATE :
        plot[0].title = SEQ_label[SEQL_MAX_POSTERIOR_IN_STATE_PROBABILITY];
        break;
      case OUT_STATE :
        plot[0].title = SEQ_label[SEQL_MAX_POSTERIOR_OUT_STATE_PROBABILITY];
        break;
      }

      if (seq->index_parameter) {
        plot[0].xrange = Range(seq->index_parameter[index][0] , seq->index_parameter[index][seq->length[index] - 1]);
        if (seq->index_parameter[index][seq->length[index] - 1] - seq->index_parameter[index][0] < TIC_THRESHOLD) {
          plot[0].xtics = 1;
        }
      }

      else {
        plot[0].xrange = Range(0 , seq->length[index] - 1);
        if (seq->length[index] - 1 < TIC_THRESHOLD) {
          plot[0].xtics = 1;
        }
      }

      plot[0].yrange = Range(0. , exp(state_seq_likelihood - seq_likelihood));

      for (i = 0;i < nb_state;i++) {
        legend.str("");
        legend << STAT_label[STATL_STATE] << " " << i;
        plot[0][i].legend = legend.str();

        plot[0][i].style = "linespoints";
      }

      // smoothed probabilities

      switch (output) {
      case SSTATE :
        plot[1].title = SEQ_label[SEQL_POSTERIOR_STATE_PROBABILITY];
        break;
      case IN_STATE :
        plot[1].title = SEQ_label[SEQL_POSTERIOR_IN_STATE_PROBABILITY];
        break;
      case OUT_STATE :
        plot[1].title = SEQ_label[SEQL_POSTERIOR_OUT_STATE_PROBABILITY];
        break;
      }

      if (seq->index_parameter) {
        plot[1].xrange = Range(seq->index_parameter[index][0] , seq->index_parameter[index][seq->length[index] - 1]);
        if (seq->index_parameter[index][seq->length[index] - 1] - seq->index_parameter[index][0] < TIC_THRESHOLD) {
          plot[1].xtics = 1;
        }
      }

      else {
        plot[1].xrange = Range(0 , seq->length[index] - 1);
        if (seq->length[index] - 1 < TIC_THRESHOLD) {
          plot[1].xtics = 1;
        }
      }

      plot[1].yrange = Range(0. , 1.);

      for (i = 0;i < nb_state;i++) {
        legend.str("");
        legend << STAT_label[STATL_STATE] << " " << i;
        plot[1][i].legend = legend.str();

        plot[1][i].style = "linespoints";
      }

      // conditional entropy profiles

      if (seq->index_parameter) {
        plot[2].xrange = Range(seq->index_parameter[index][0] , seq->index_parameter[index][seq->length[index] - 1]);
        if (seq->index_parameter[index][seq->length[index] - 1] - seq->index_parameter[index][0] < TIC_THRESHOLD) {
          plot[2].xtics = 1;
        }
      }

      else {
        plot[2].xrange = Range(0 , seq->length[index] - 1);
        if (seq->length[index] - 1 < TIC_THRESHOLD) {
          plot[2].xtics = 1;
        }
      }

      plot[2].yrange = Range(0. , max_marginal_entropy);

      plot[2][0].legend = SEQ_label[SEQL_CONDITIONAL_ENTROPY];
      plot[2][0].style = "linespoints";

      plot[2][1].legend = SEQ_label[SEQL_MARGINAL_ENTROPY];
      plot[2][1].style = "linespoints";

      // partial entropy profiles

      if (seq->index_parameter) {
        plot[3].xrange = Range(seq->index_parameter[index][0] , seq->index_parameter[index][seq->length[index] - 1]);
        if (seq->index_parameter[index][seq->length[index] - 1] - seq->index_parameter[index][0] < TIC_THRESHOLD) {
          plot[3].xtics = 1;
        }
      }

      else {
        plot[3].xrange = Range(0 , seq->length[index] - 1);
        if (seq->length[index] - 1 < TIC_THRESHOLD) {
          plot[3].xtics = 1;
        }
      }

      plot[3].yrange = Range(0. ,entropy);

      plot[3][0].legend = SEQ_label[SEQL_PARTIAL_STATE_SEQUENCE_ENTROPY];
      plot[3][0].style = "linespoints";
    }

    delete seq;
    delete hsmarkov;
  }

  return plot_set;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of state and entropy profiles using the forward-backward algorithm,
 *         of state profiles using the Viterbi forward-backward algorithm and
 *         plot of the results.
 *
 *  \param[in] error      reference on a StatError object,
 *  \param[in] identifier sequence identifier,
 *  \param[in] output     output type.
 *
 *  \return               MultiPlotSet object.
 */
/*--------------------------------------------------------------*/

MultiPlotSet* HiddenSemiMarkov::state_profile_plotable_write(StatError &error , int identifier ,
                                                             state_profile output) const

{
  MultiPlotSet *plot_set;


  error.init();

  if (!semi_markov_data) {
    plot_set = NULL;
    error.update(STAT_error[STATR_NO_DATA]);
  }
  else {
    plot_set = state_profile_plotable_write(error , *semi_markov_data , identifier , output);
  }

  return plot_set;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the most probable state sequences.
 *
 *  \param[in] error               reference on a StatError object,
 *  \param[in] display             flag for displaying the posterior state sequence probabilities,
 *  \param[in] iseq                reference on a MarkovianSequences object,
 *  \param[in] characteristic_flag flag on the computation of the characteristic distributions.
 *
 *  \return                        SemiMarkovData object.
 */
/*--------------------------------------------------------------*/

SemiMarkovData* HiddenSemiMarkov::state_sequence_computation(StatError &error , bool display ,
                                                             const MarkovianSequences &iseq ,
                                                             bool characteristic_flag) const

{
  bool status = true;
  int i;
  int nb_value;
  HiddenSemiMarkov *hsmarkov;
  SemiMarkovData *seq;


  seq = NULL;
  error.init();

  for (i = 0;i < iseq.nb_variable;i++) {
    if ((iseq.type[i] != INT_VALUE) && (iseq.type[i] != REAL_VALUE)) {
      status = false;
      ostringstream error_message , correction_message;
      error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                    << STAT_error[STATR_VARIABLE_TYPE];
      correction_message << STAT_variable_word[INT_VALUE] << " or "
                         << STAT_variable_word[REAL_VALUE];
      error.correction_update((error_message.str()).c_str() , (correction_message.str()).c_str());
    }
  }

  if (nb_output_process != iseq.nb_variable) {
    status = false;
    error.update(STAT_error[STATR_NB_OUTPUT_PROCESS]);
  }

  else {
    for (i = 0;i < nb_output_process;i++) {
      if ((categorical_process[i]) || (discrete_parametric_process[i])) {
        if (iseq.type[i] == REAL_VALUE) {
          status = false;
          ostringstream error_message;
          error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                        << STAT_error[STATR_VARIABLE_TYPE];
          error.correction_update((error_message.str()).c_str() , STAT_variable_word[INT_VALUE]);
        }

        else {
          if (iseq.min_value[i] < 0) {
            status = false;
            ostringstream error_message;
            error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                          << STAT_error[STATR_POSITIVE_MIN_VALUE];
            error.update((error_message.str()).c_str());
          }

          if (!(iseq.marginal_distribution[i])) {
            status = false;
            ostringstream error_message;
            error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                          << STAT_error[STATR_MARGINAL_FREQUENCY_DISTRIBUTION];
            error.update((error_message.str()).c_str());
          }

          else {
            if (categorical_process[i]) {
              nb_value = categorical_process[i]->nb_value;
            }
            else {
              nb_value = discrete_parametric_process[i]->nb_value;
            }

            if (nb_value < iseq.marginal_distribution[i]->nb_value) {
              status = false;
              ostringstream error_message;
              error_message << STAT_label[STATL_OUTPUT_PROCESS] << " " << i + 1 << ": "
                            << STAT_error[STATR_NB_OUTPUT];
              error.update((error_message.str()).c_str());
            }
          }
        }
      }
    }
  }

  if (status) {
    seq = new SemiMarkovData(iseq , ADD_STATE_VARIABLE , (type == EQUILIBRIUM ? true : false));

    seq->semi_markov = new SemiMarkov(*this , false);

    hsmarkov = new HiddenSemiMarkov(*this , false);

    hsmarkov->forward_backward(*seq);

    hsmarkov->create_cumul();
    hsmarkov->log_computation();
    hsmarkov->viterbi(*seq);

    delete hsmarkov;

    // extraction of the characteristics of the sequences and
    // computation of the characteristic distributions of the model

    if (seq->restoration_likelihood == D_INF) {
      delete seq;
      seq = NULL;
      error.update(SEQ_error[SEQR_STATE_SEQUENCE_COMPUTATION_FAILURE]);
    }

    else {
      seq->likelihood = likelihood_computation(iseq , seq->posterior_probability);

      if ((display) && (seq->nb_sequence <= POSTERIOR_PROBABILITY_NB_SEQUENCE)) {
        int j;
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

/*      seq->min_value_computation(0);
      seq->max_value_computation(0); */

      seq->min_value[0] = 0;
      seq->max_value[0] = nb_state - 1;
      seq->build_marginal_frequency_distribution(0);
      seq->build_characteristic(0 , true , (type == EQUILIBRIUM ? true : false));

      seq->build_transition_count(this);
      seq->build_observation_frequency_distribution(nb_state);

/*      if ((seq->max_value[0] < nb_state - 1) || (!(seq->characteristics[0]))) {
        delete seq;
        seq = NULL;
        error.update(SEQ_error[SEQR_STATES_NOT_REPRESENTED]);
      }

      else if (characteristic_flag) { */
      if (characteristic_flag) {
        seq->semi_markov->characteristic_computation(*seq , true);
      }
    }
  }

  return seq;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Comparison of hidden semi-Markov chains for a sample of sequences.
 *
 *  \param[in] error     reference on a StatError object,
 *  \param[in] display   flag for displaying the results of model comparison,
 *  \param[in] nb_model  number of hidden semi-Markov chains,
 *  \param[in] ihsmarkov pointer on HiddenSemiMarkov objects,
 *  \param[in] algorithm type of algorithm (FORWARD/VITERBI),
 *  \param[in] path      file path.
 *
 *  \return              error status.
 */
/*--------------------------------------------------------------*/

bool MarkovianSequences::comparison(StatError &error , bool display , int nb_model ,
                                    const HiddenSemiMarkov **ihsmarkov ,
                                    latent_structure_algorithm algorithm ,
                                    const string path) const

{
  bool status = true;
  int i , j;
  int nb_value;
  double **likelihood;
  HiddenSemiMarkov **hsmarkov;
  SemiMarkovData *seq;


  error.init();

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

  for (i = 0;i < nb_model;i++) {
    if (ihsmarkov[i]->nb_output_process != nb_variable) {
      status = false;
      ostringstream error_message;
      error_message << SEQ_label[SEQL_HIDDEN_SEMI_MARKOV_CHAIN] << " " << i + 1 << ": "
                    << STAT_error[STATR_NB_OUTPUT_PROCESS];
      error.update((error_message.str()).c_str());
    }

    else {
      for (j = 0;j < nb_variable;j++) {
        if ((ihsmarkov[i]->categorical_process[j]) || (ihsmarkov[i]->discrete_parametric_process[j])) {
          if (type[j] == REAL_VALUE) {
            status = false;
            ostringstream error_message;
            error_message << STAT_label[STATL_VARIABLE] << " " << j + 1 << ": "
                          << STAT_error[STATR_VARIABLE_TYPE];
            error.correction_update((error_message.str()).c_str() , STAT_variable_word[INT_VALUE]);
          }

          else {
            if (min_value[j] < 0) {
              status = false;
              ostringstream error_message;
              error_message << STAT_label[STATL_VARIABLE] << " " << j + 1 << ": "
                            << STAT_error[STATR_POSITIVE_MIN_VALUE];
              error.update((error_message.str()).c_str());
            }

            if (!marginal_distribution[j]) {
              status = false;
              ostringstream error_message;
              error_message << STAT_label[STATL_VARIABLE] << " " << j + 1 << ": "
                            << STAT_error[STATR_MARGINAL_FREQUENCY_DISTRIBUTION];
              error.update((error_message.str()).c_str());
            }

            else {
              if (ihsmarkov[i]->categorical_process[j]) {
                nb_value = ihsmarkov[i]->categorical_process[j]->nb_value;
              }
              else {
                nb_value = ihsmarkov[i]->discrete_parametric_process[j]->nb_value;
              }

              if (nb_value < marginal_distribution[j]->nb_value) {
                status = false;
                ostringstream error_message;
                error_message << SEQ_label[SEQL_HIDDEN_SEMI_MARKOV_CHAIN] << " " << i + 1 << ": "
                              << STAT_label[STATL_OUTPUT_PROCESS] << " " << j + 1 << ": "
                              << STAT_error[STATR_NB_OUTPUT];
                error.update((error_message.str()).c_str());
              }
            }
          }
        }
      }
    }
  }

  if (status) {
    likelihood = new double*[nb_sequence];
    for (i = 0;i < nb_sequence;i++) {
      likelihood[i] = new double[nb_model];
    }

    hsmarkov = new HiddenSemiMarkov*[nb_model];
    for (i = 0;i < nb_model;i++) {
      hsmarkov[i] = new HiddenSemiMarkov(*(ihsmarkov[i]) , false , false);
    }

    if (algorithm == VITERBI) {
      for (i = 0;i < nb_model;i++) {
        hsmarkov[i]->create_cumul();
        hsmarkov[i]->log_computation();
      }

      seq = new SemiMarkovData(*this);
    }

    // for each sequence, computation of the log-likelihood for the observed sequence (FORWARD) or
    // of the log-likelihood for the most probable state sequence (VITERBI) for each model

    for (i = 0;i < nb_sequence;i++) {
      for (j = 0;j < nb_model;j++) {
        switch (algorithm) {
        case FORWARD :
          likelihood[i][j] = hsmarkov[j]->likelihood_computation(*this , NULL , i);
          break;
        case VITERBI :
          likelihood[i][j] = hsmarkov[j]->viterbi(*seq , NULL , i);
          break;
        }
      }
    }

    if (display) {
      likelihood_write(cout , nb_model , likelihood ,
                       SEQ_label[SEQL_HIDDEN_SEMI_MARKOV_CHAIN] , true , algorithm);
    }
    if (!path.empty()) {
      status = likelihood_write(error , path , nb_model , likelihood ,
                                SEQ_label[SEQL_HIDDEN_SEMI_MARKOV_CHAIN] , algorithm);
    }

    for (i = 0;i < nb_sequence;i++) {
      delete [] likelihood[i];
    }
    delete [] likelihood;

    for (i = 0;i < nb_model;i++) {
      delete hsmarkov[i];
    }
    delete [] hsmarkov;

    if (algorithm == VITERBI) {
      delete seq;
    }
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Simulation using a hidden semi-Markov chain.
 *
 *  \param[in] error               reference on a StatError object,
 *  \param[in] length_distribution sequence length frequency distribution,
 *  \param[in] counting_flag       flag on the computation of the counting distributions,
 *  \param[in] divergence_flag     flag on the computation of the Kullback-Leibler divergence.
 *
 *  \return                        SemiMarkovData object.
 */
/*--------------------------------------------------------------*/

SemiMarkovData* HiddenSemiMarkov::simulation(StatError &error ,
                                             const FrequencyDistribution &length_distribution ,
                                             bool counting_flag , bool divergence_flag) const

{
  int i;
  MarkovianSequences *observed_seq;
  SemiMarkovData *seq;


  seq = SemiMarkov::simulation(error , length_distribution , counting_flag , divergence_flag);

  if ((seq) && (!divergence_flag)) {
    seq->posterior_probability = new double[seq->nb_sequence];
    for (i = 0;i < seq->nb_sequence;i++) {
      seq->posterior_probability[i] = SemiMarkov::likelihood_computation(*seq , i);
    }

    observed_seq = seq->remove_variable_1();
    seq->likelihood = likelihood_computation(*observed_seq , seq->posterior_probability);
    delete observed_seq;

    forward_backward(*seq);
  }

  return seq;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Simulation using a hidden semi-Markov chain.
 *
 *  \param[in] error         reference on a StatError object,
 *  \param[in] nb_sequence   number of sequences,
 *  \param[in] length        sequence length,
 *  \param[in] counting_flag flag on the computation of the counting distributions.
 *
 *  \return                  SemiMarkovData object.
 */
/*--------------------------------------------------------------*/

SemiMarkovData* HiddenSemiMarkov::simulation(StatError &error , int nb_sequence ,
                                             int length , bool counting_flag) const

{
  int i;
  MarkovianSequences *observed_seq;
  SemiMarkovData *seq;


  seq = SemiMarkov::simulation(error , nb_sequence , length , counting_flag);

  if (seq) {
    seq->posterior_probability = new double[seq->nb_sequence];
    for (i = 0;i < seq->nb_sequence;i++) {
      seq->posterior_probability[i] = SemiMarkov::likelihood_computation(*seq , i);
    }

    observed_seq = seq->remove_variable_1();
    seq->likelihood = likelihood_computation(*observed_seq , seq->posterior_probability);
    delete observed_seq;

    forward_backward(*seq);
  }

  return seq;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Simulation using a hidden semi-Markov chain.
 *
 *  \param[in] error         reference on a StatError object,
 *  \param[in] nb_sequence   number of sequences,
 *  \param[in] iseq          reference on a MarkovianSequences object,
 *  \param[in] counting_flag flag on the computation of the counting distributions.
 *
 *  \return                  SemiMarkovData object.
 */
/*--------------------------------------------------------------*/

SemiMarkovData* HiddenSemiMarkov::simulation(StatError &error , int nb_sequence ,
                                             const MarkovianSequences &iseq , bool counting_flag) const

{
  int i;
  MarkovianSequences *observed_seq;
  SemiMarkovData *seq;


  seq = SemiMarkov::simulation(error , nb_sequence , iseq , counting_flag);

  if (seq) {
    seq->posterior_probability = new double[seq->nb_sequence];
    for (i = 0;i < seq->nb_sequence;i++) {
      seq->posterior_probability[i] = SemiMarkov::likelihood_computation(*seq , i);
    }

    observed_seq = seq->remove_variable_1();
    seq->likelihood = likelihood_computation(*observed_seq , seq->posterior_probability);
    delete observed_seq;

    forward_backward(*seq);
  }

  return seq;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of Kullback-Leibler divergences between hidden semi-Markov chains.
 *
 *  \param[in] error               reference on a StatError object,
 *  \param[in] display             flag for displaying the matrix of pairwise distances between models,
 *  \param[in] nb_model            number of hidden semi-Markov chains,
 *  \param[in] ihsmarkov           pointer on HiddenSemiMarkov objects,
 *  \param[in] length_distribution sequence length frequency distribution,
 *  \param[in] path                file path.
 *
 *  \return                        DistanceMatrix object.
 */
/*--------------------------------------------------------------*/

DistanceMatrix* HiddenSemiMarkov::divergence_computation(StatError &error , bool display ,
                                                         int nb_model , const HiddenSemiMarkov **ihsmarkov ,
                                                         FrequencyDistribution **length_distribution ,
                                                         const string path) const

{
  bool status = true , lstatus;
  int i , j , k;
  int cumul_length , nb_failure;
  double **likelihood;
  long double divergence;
  const HiddenSemiMarkov **hsmarkov;
  MarkovianSequences *seq;
  SemiMarkovData *simul_seq;
  DistanceMatrix *dist_matrix;
  ofstream *out_file;


  dist_matrix = NULL;
  error.init();

  for (i = 0;i < nb_model - 1;i++) {
    if (ihsmarkov[i]->type != type) {
      status = false;
      ostringstream error_message;
      error_message << SEQ_label[SEQL_HIDDEN_SEMI_MARKOV_CHAIN] << " " << i + 2 << ": "
                    << SEQ_error[SEQR_MODEL_TYPE];
      error.update((error_message.str()).c_str());
    }

    if (ihsmarkov[i]->nb_output_process != nb_output_process) {
      status = false;
      ostringstream error_message;
      error_message << SEQ_label[SEQL_HIDDEN_SEMI_MARKOV_CHAIN] << " " << i + 2 << ": "
                    << STAT_error[STATR_NB_OUTPUT_PROCESS];
      error.update((error_message.str()).c_str());
    }

    else {
      for (j = 0;j < nb_output_process;j++) {
        if ((categorical_process[j]) && (ihsmarkov[i]->categorical_process[j]) &&
            (ihsmarkov[i]->categorical_process[j]->nb_value != categorical_process[j]->nb_value)) {
          status = false;
          ostringstream error_message;
          error_message << SEQ_label[SEQL_HIDDEN_SEMI_MARKOV_CHAIN] << " " << i + 2 << ": "
                        << STAT_label[STATL_OUTPUT_PROCESS] << " " << j << " "
                        << STAT_error[STATR_NB_OUTPUT];
          error.update((error_message.str()).c_str());
        }

        if (((continuous_parametric_process[j]) && (!(ihsmarkov[i]->continuous_parametric_process[j]))) ||
            ((!continuous_parametric_process[j]) && (ihsmarkov[i]->continuous_parametric_process[j]))) {
          status = false;
          ostringstream error_message;
          error_message << SEQ_label[SEQL_HIDDEN_MARKOV_CHAIN] << " " << i + 2 << ": "
                        << STAT_label[STATL_OUTPUT_PROCESS] << " " << j << " "
                        << SEQ_error[SEQR_OUTPUT_PROCESS_TYPE];
          error.update((error_message.str()).c_str());
        }
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

    hsmarkov = new const HiddenSemiMarkov*[nb_model];

    hsmarkov[0] = this;
    for (i = 1;i < nb_model;i++) {
      hsmarkov[i] = ihsmarkov[i - 1];
    }

    dist_matrix = new DistanceMatrix(nb_model , SEQ_label[SEQL_HIDDEN_SEMI_MARKOV_CHAIN]);

    for (i = 0;i < nb_model;i++) {

      // generation of a sample of sequences using a hidden semi-Markov chain

      simul_seq = hsmarkov[i]->simulation(error , *length_distribution[i] , false , true);
      seq = simul_seq->remove_variable_1();

      likelihood = new double*[seq->nb_sequence];
      for (j = 0;j < seq->nb_sequence;j++) {
        likelihood[j] = new double[nb_model];
      }

      for (j = 0;j < seq->nb_sequence;j++) {
        likelihood[j][i] = hsmarkov[i]->likelihood_computation(*seq , NULL , j);

        if ((display) && (likelihood[j][i] == D_INF)) {
          cout << "\nERROR - " << SEQ_error[SEQR_REFERENCE_MODEL] << ": " << i + 1 << endl;
        }
      }

      // computation of the log-likelihood of each hidden semi-Markov chain for the sample of sequences

      for (j = 0;j < nb_model;j++) {
        if (j != i) {
          divergence = 0.;
          cumul_length = 0;
          nb_failure = 0;

          for (k = 0;k < seq->nb_sequence;k++) {
            likelihood[k][j] = hsmarkov[j]->likelihood_computation(*seq , NULL , k);

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
        }
      }

      if (display) {
        cout << SEQ_label[SEQL_HIDDEN_SEMI_MARKOV_CHAIN] << " " << i + 1 << ": " << seq->nb_sequence << " "
             << SEQ_label[SEQL_SIMULATED] << " " << SEQ_label[seq->nb_sequence == 1 ? SEQL_SEQUENCE : SEQL_SEQUENCES] << endl;
        seq->likelihood_write(cout , nb_model , likelihood , SEQ_label[SEQL_HIDDEN_SEMI_MARKOV_CHAIN]);
      }
      if (out_file) {
        *out_file << SEQ_label[SEQL_HIDDEN_SEMI_MARKOV_CHAIN] << " " << i + 1 << ": " << seq->nb_sequence << " "
                  << SEQ_label[SEQL_SIMULATED] << " " << SEQ_label[seq->nb_sequence == 1 ? SEQL_SEQUENCE : SEQL_SEQUENCES] << endl;
        seq->likelihood_write(*out_file , nb_model , likelihood , SEQ_label[SEQL_HIDDEN_SEMI_MARKOV_CHAIN]);
      }

      for (j = 0;j < seq->nb_sequence;j++) {
        delete [] likelihood[j];
      }
      delete [] likelihood;

      delete seq;
      delete simul_seq;
    }

    if (out_file) {
      out_file->close();
      delete out_file;
    }

    delete hsmarkov;
  }

  return dist_matrix;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of Kullback-Leibler divergences between hidden semi-Markov chains.
 *
 *  \param[in] error       reference on a StatError object,
 *  \param[in] display     flag for displaying the matrix of pairwise distances between models,
 *  \param[in] nb_model    number of hidden semi-Markov chains,
 *  \param[in] hsmarkov    pointer on HiddenSemiMarkov objects,
 *  \param[in] nb_sequence number of generated sequences,
 *  \param[in] length      sequence length,
 *  \param[in] path        file path.
 *
 *  \return                DistanceMatrix object.
 */
/*--------------------------------------------------------------*/

DistanceMatrix* HiddenSemiMarkov::divergence_computation(StatError &error , bool display ,
                                                         int nb_model , const HiddenSemiMarkov **hsmarkov ,
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

    dist_matrix = divergence_computation(error , display , nb_model , hsmarkov , length_distribution , path);

    for (i = 0;i < nb_model;i++) {
      delete length_distribution[i];
    }
    delete [] length_distribution;
  }

  return dist_matrix;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of Kullback-Leibler divergences between hidden semi-Markov chains.
 *
 *  \param[in] error       reference on a StatError object,
 *  \param[in] display     flag for displaying the matrix of pairwise distances between models,
 *  \param[in] nb_model    number of hidden semi-Markov chains,
 *  \param[in] hsmarkov    pointer on HiddenSemiMarkov objects,
 *  \param[in] nb_sequence number of generated sequences,
 *  \param[in] seq         pointer on MarkovianSequences objects,
 *  \param[in] path        file path.
 *
 *  \return                DistanceMatrix object.
 */
/*--------------------------------------------------------------*/

DistanceMatrix* HiddenSemiMarkov::divergence_computation(StatError &error , bool display ,
                                                         int nb_model , const HiddenSemiMarkov **hsmarkov ,
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

    dist_matrix = divergence_computation(error , display , nb_model , hsmarkov , length_distribution , path);

    for (i = 0;i < nb_model;i++) {
      delete length_distribution[i];
    }
    delete [] length_distribution;
  }

  return dist_matrix;
}


};  // namespace sequence_analysis
