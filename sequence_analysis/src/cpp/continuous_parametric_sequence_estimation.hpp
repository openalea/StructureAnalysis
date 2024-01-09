/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       V-Plants: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2015 CIRAD/INRA/Inria Virtual Plants
 *
 *       File author(s): Yann Guedon (yann.guedon@cirad.fr)
 *
 *       $Source$
 *       $Id: continuous_distribution_sequence_estimation.hpp 12646 2012-08-03 08:12:47Z guedon $
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



#ifndef CONTINUOUS_PARAMETRIC_SEQUENCE_ESTIMATION_HPP
#define CONTINUOUS_PARAMETRIC_SEQUENCE_ESTIMATION_HPP


#include <math.h>

#include <boost/math/special_functions/digamma.hpp>

#include "stat_tool/stat_tools.h"
#include "stat_tool/curves.h"
#include "stat_tool/distribution.h"
#include "stat_tool/markovian.h"
#include "stat_tool/vectors.h"
#include "stat_tool/distance_matrix.h"
#include "stat_tool/stat_label.h"

#include "sequences.h"

using namespace std;
using namespace boost::math;
using namespace stat_tool;


namespace sequence_analysis {


// extern double von_mises_concentration_computation(double mean_direction);



/*--------------------------------------------------------------*
 *
 *  Estimation des parametres des lois gamma d'un processus d'observation.
 *
 *  arguments : comptage des etats, variable, pointeur sur un processus d'observation continu,
 *              iteration EM.
 *
 *--------------------------------------------------------------*/

template <typename Type>
void MarkovianSequences::gamma_estimation(Type ***state_sequence_count , int variable ,
                                          ContinuousParametricProcess *process , int iter) const

{
  register int i , j , k;
  double diff , log_geometric_mean , *zero_mass , *mean , *variance;
  Type *state_frequency;


  state_frequency = new Type[process->nb_state];
  zero_mass = new double[process->nb_state];
  mean = new double[process->nb_state];
  variance = new double[process->nb_state];

  for (i = 0;i < process->nb_state;i++) {
    zero_mass[i] = 0.;
    mean[i] = 0.;
    state_frequency[i] = 0;
  }

  switch (type[variable]) {

  case INT_VALUE : {
    for (i = 0;i < nb_sequence;i++) {
      for (j = 0;j < length[i];j++) {
        for (k = 0;k < process->nb_state;k++) {
          if (int_sequence[i][variable][j] == 0) {
            zero_mass[k] += state_sequence_count[i][j][k];
          }
          else {
            mean[k] += state_sequence_count[i][j][k] * int_sequence[i][variable][j];
          }
          state_frequency[k] += state_sequence_count[i][j][k];
        }
      }
    }
    break;
  }

  case REAL_VALUE : {
    for (i = 0;i < nb_sequence;i++) {
      for (j = 0;j < length[i];j++) {
        for (k = 0;k < process->nb_state;k++) {
          if (real_sequence[i][variable][j] == 0.) {
            zero_mass[k] += state_sequence_count[i][j][k];
          }
          else {
            mean[k] += state_sequence_count[i][j][k] * real_sequence[i][variable][j];
          }
          state_frequency[k] += state_sequence_count[i][j][k];
        }
      }
    }
    break;
  }
  }

  for (i = 0;i < process->nb_state;i++) {
    if (state_frequency[i] > 0) {
      mean[i] /= state_frequency[i];
    }
    variance[i] = 0.;
  }

  switch (type[variable]) {

  case INT_VALUE : {
    for (i = 0;i < nb_sequence;i++) {
      for (j = 0;j < length[i];j++) {
        for (k = 0;k < process->nb_state;k++) {
          diff = int_sequence[i][variable][j] - mean[k];
          variance[k] += state_sequence_count[i][j][k] * diff * diff;
        }
      }
    }
    break;
  }

  case REAL_VALUE : {
    for (i = 0;i < nb_sequence;i++) {
      for (j = 0;j < length[i];j++) {
        for (k = 0;k < process->nb_state;k++) {
          diff = real_sequence[i][variable][j] - mean[k];
          variance[k] += state_sequence_count[i][j][k] * diff * diff;
        }
      }
    }
    break;
  }
  }

  for (i = 0;i < process->nb_state;i++) {
    if (state_frequency[i] > 0) {
      variance[i] /= state_frequency[i];
//      variance[i] /= (state_frequency[i] - 1);
    }
  }

  for (i = 0;i < process->nb_state;i++) {

#   ifdef DEBUG
    if ((iter >= 5) && (state_frequency[i] > 0)) {
      cout << "\n" << STAT_word[STATW_STATE] << " " << i << " : "
           << zero_mass[i] << ", " << state_frequency[i] << " | "
           << zero_mass[i] / state_frequency[i] << endl;
    }
#   endif

    if (state_frequency[i] > 0) {
      if (zero_mass[i] / state_frequency[i] > GAMMA_ZERO_FREQUENCY_THRESHOLD) {
        process->observation[i]->shape = 0;
        process->observation[i]->scale = D_DEFAULT;
      }

      else {
        if (variance[i] > 0.) {
/*          if (sqrt(variance[i]) < mean[i] * GAMMA_VARIATION_COEFF_THRESHOLD) {
            variance[i] = mean[i] * mean[i] * GAMMA_VARIATION_COEFF_THRESHOLD * GAMMA_VARIATION_COEFF_THRESHOLD;
          } */

//          process->observation[i]->shape = mean[i] * mean[i] / variance[i];
//          process->observation[i]->scale = variance[i] / mean[i];

          // Hawang & Huang (2012), Ann. Inst. Statist. Math. 54(4), 840-847

          process->observation[i]->shape = mean[i] * mean[i] / variance[i] - 1. / (double)state_frequency[i];
          process->observation[i]->scale = mean[i] / process->observation[i]->shape;

#         ifdef DEBUG    // essai pour eviter les très petits parametres de forme
          if ((process->observation[i]->shape < 1.) && (mean[i] < 5.)) {
            process->observation[i]->shape = 1.;
            process->observation[i]->scale = mean[i] / process->observation[i]->shape;
          }
#         endif

          if ((process->observation[i]->shape >= GAMMA_SHAPE_PARAMETER_THRESHOLD) &&
              (state_frequency[i] < GAMMA_FREQUENCY_THRESHOLD)) {
            log_geometric_mean = 0.;
            state_frequency[i] = 0;

            switch (type[variable]) {

            case INT_VALUE : {
              for (j = 0;j < nb_sequence;j++) {
                for (k = 0;k < length[j];k++) {
                  if (int_sequence[j][variable][k] > 0) {
                    log_geometric_mean += state_sequence_count[j][k][i] * log(int_sequence[j][variable][k]);
                    state_frequency[i] += state_sequence_count[j][k][i];
                  }
                }
              }
              break;
            }

            case REAL_VALUE : {
              for (j = 0;j < nb_sequence;j++) {
                for (k = 0;k < length[j];k++) {
                  if (real_sequence[j][variable][k] > 0.) {
                    log_geometric_mean += state_sequence_count[j][k][i] * log(real_sequence[j][variable][k]);
                    state_frequency[i] += state_sequence_count[j][k][i];
                  }
                }
              }
              break;
            }
            }

            log_geometric_mean /= state_frequency[i];
/*            j = 0;   a revoir

#           ifdef DEBUG
            cout << "\n" << STAT_word[STATW_STATE] << " " << i << "   "
                 << STAT_word[STATW_SHAPE] << " : " << process->observation[i]->shape << "   "
                 << STAT_word[STATW_SCALE] << " : " << process->observation[i]->scale << endl;
#           endif

            do {
              process->observation[i]->scale = exp(log_geometric_mean - digamma(process->observation[i]->shape));
              process->observation[i]->shape = mean[i] / process->observation[i]->scale;
              j++;

#             ifdef DEBUG
              cout << STAT_word[STATW_SHAPE] << " : " << process->observation[i]->shape << "   "
                   << STAT_word[STATW_SCALE] << " : " << process->observation[i]->scale << endl;
#             endif

            }
            while (j < MIN(GAMMA_ITERATION_FACTOR * iter , GAMMA_MAX_NB_ITERATION)); */

            // approximations Johnson, Kotz & Balakrishnan, Continuous Univariate Distributions, vol. 1, 2nd ed., pp. 361-362

//            process->observation[i]->shape = mean[i] / (2 * (mean[i] - exp(log_geometric_mean))) - 1./12.;
            diff = log(mean[i]) - log_geometric_mean;
            process->observation[i]->shape = (1 + sqrt(1 + 4 * diff / 3)) / (4 * diff);
            process->observation[i]->scale = mean[i] / process->observation[i]->shape;
          }
        }

        else {
          process->observation[i]->shape = GAMMA_MIN_SHAPE_PARAMETER;
          process->observation[i]->scale = GAMMA_DEFAULT_SCALE_PARAMETER;
        }
      }
    }

    else {
      process->observation[i]->shape = D_DEFAULT;
      process->observation[i]->scale = D_DEFAULT;
    }
  }

  delete [] state_frequency;
  delete [] zero_mass;
  delete [] mean;
  delete [] variance;
}


/*--------------------------------------------------------------*
 *
 *  Estimation des parametres des lois zero-inflated gamma d'un processus d'observation.
 *
 *  arguments : comptage des etats, variable, pointeur sur un processus d'observation continu,
 *              iteration EM.
 *
 *--------------------------------------------------------------*/

template <typename Type>
void MarkovianSequences::zero_inflated_gamma_estimation(Type ***state_sequence_count , int variable ,
                                                        ContinuousParametricProcess *process , int iter) const

{
  register int i , j , k;
  double diff , log_geometric_mean , *zero_mass , *mean , *variance;
  Type *state_frequency;


  state_frequency = new Type[process->nb_state];
  zero_mass = new double[process->nb_state];
  mean = new double[process->nb_state];
  variance = new double[process->nb_state];

  for (i = 0;i < process->nb_state;i++) {
    zero_mass[i] = 0.;
    mean[i] = 0.;
    state_frequency[i] = 0;
  }

  switch (type[variable]) {

  case INT_VALUE : {
    for (i = 0;i < nb_sequence;i++) {
      for (j = 0;j < length[i];j++) {
        for (k = 0;k < process->nb_state;k++) {
          if (int_sequence[i][variable][j] == 0) {
            zero_mass[k] += state_sequence_count[i][j][k];
          }
          else {
            mean[k] += state_sequence_count[i][j][k] * int_sequence[i][variable][j];
            state_frequency[k] += state_sequence_count[i][j][k];
          }
        }
      }
    }
    break;
  }

  case REAL_VALUE : {
    for (i = 0;i < nb_sequence;i++) {
      for (j = 0;j < length[i];j++) {
        for (k = 0;k < process->nb_state;k++) {
          if (real_sequence[i][variable][j] == 0.) {
            zero_mass[k] += state_sequence_count[i][j][k];
          }
          else {
            mean[k] += state_sequence_count[i][j][k] * real_sequence[i][variable][j];
            state_frequency[k] += state_sequence_count[i][j][k];
          }
        }
      }
    }
    break;
  }
  }

  for (i = 0;i < process->nb_state;i++) {
    if (state_frequency[i] > 0) {
      mean[i] /= state_frequency[i];
    }
    variance[i] = 0.;
  }

  switch (type[variable]) {

  case INT_VALUE : {
    for (i = 0;i < nb_sequence;i++) {
      for (j = 0;j < length[i];j++) {
        for (k = 0;k < process->nb_state;k++) {
          if (int_sequence[i][variable][j] > 0) {
            diff = int_sequence[i][variable][j] - mean[k];
            variance[k] += state_sequence_count[i][j][k] * diff * diff;
          }
        }
      }
    }
    break;
  }

  case REAL_VALUE : {
    for (i = 0;i < nb_sequence;i++) {
      for (j = 0;j < length[i];j++) {
        for (k = 0;k < process->nb_state;k++) {
          if (real_sequence[i][variable][j] > 0.) {
            diff = real_sequence[i][variable][j] - mean[k];
            variance[k] += state_sequence_count[i][j][k] * diff * diff;
          }
        }
      }
    }
    break;
  }
  }

  for (i = 0;i < process->nb_state;i++) {
    if (state_frequency[i] > 0) {
      variance[i] /= state_frequency[i];
//      variance[i] /= (state_frequency[i] - 1);
    }
  }

  for (i = 0;i < process->nb_state;i++) {
    if (zero_mass[i] + state_frequency[i] > 0) {
      if (zero_mass[i] / (zero_mass[i] + state_frequency[i]) > GAMMA_ZERO_FREQUENCY_THRESHOLD) {
        process->observation[i]->zero_probability = 1.;
        process->observation[i]->shape = D_DEFAULT;
        process->observation[i]->scale = D_DEFAULT;
      }

      else {
        process->observation[i]->zero_probability = zero_mass[i] / (zero_mass[i] + state_frequency[i]);

        if (variance[i] > 0.) {
/*          if (sqrt(variance[i]) < mean[i] * GAMMA_VARIATION_COEFF_THRESHOLD) {
            variance[i] = mean[i] * mean[i] * GAMMA_VARIATION_COEFF_THRESHOLD * GAMMA_VARIATION_COEFF_THRESHOLD;
          } */

//          process->observation[i]->shape = mean[i] * mean[i] / variance[i];
//          process->observation[i]->scale = variance[i] / mean[i];

          // Hawang & Huang (2012), Ann. Inst. Statist. Math. 54(4), 840-847

          process->observation[i]->shape = mean[i] * mean[i] / variance[i] - 1. / (double)state_frequency[i];
          process->observation[i]->scale = mean[i] / process->observation[i]->shape;

#         ifdef DEBUG    // essai pour eviter la bimodalite
          if ((iter > 5) && (process->observation[i]->zero_probability > 0.5) &&
              (process->observation[i]->shape > 1.)) {
            process->observation[i]->shape = 1.;
            process->observation[i]->scale = mean[i] / process->observation[i]->shape;
          }
#         endif

          if ((process->observation[i]->shape >= GAMMA_SHAPE_PARAMETER_THRESHOLD) &&
              (state_frequency[i] < GAMMA_FREQUENCY_THRESHOLD)) {
            log_geometric_mean = 0.;

            switch (type[variable]) {

            case INT_VALUE : {
              for (j = 0;j < nb_sequence;j++) {
                for (k = 0;k < length[j];k++) {
                  if (int_sequence[j][variable][k] > 0) {
                    log_geometric_mean += state_sequence_count[j][k][i] * log(int_sequence[j][variable][k]);
                  }
                }
              }
              break;
            }

            case REAL_VALUE : {
              for (j = 0;j < nb_sequence;j++) {
                for (k = 0;k < length[j];k++) {
                  if (real_sequence[j][variable][k] > 0.) {
                    log_geometric_mean += state_sequence_count[j][k][i] * log(real_sequence[j][variable][k]);
                  }
                }
              }
              break;
            }
            }

            log_geometric_mean /= state_frequency[i];
/*            j = 0;   a revoir

#           ifdef DEBUG
            cout << "\n" << STAT_word[STATW_STATE] << " " << i << "   "
                 << STAT_word[STATW_SHAPE] << " : " << process->observation[i]->shape << "   "
                 << STAT_word[STATW_SCALE] << " : " << process->observation[i]->scale << endl;
#           endif

            do {
              process->observation[i]->scale = exp(log_geometric_mean - digamma(process->observation[i]->shape));
              process->observation[i]->shape = mean[i] / process->observation[i]->scale;
              j++;

#             ifdef DEBUG
              cout << STAT_word[STATW_SHAPE] << " : " << process->observation[i]->shape << "   "
                   << STAT_word[STATW_SCALE] << " : " << process->observation[i]->scale << endl;
#             endif

            }
            while (j < MIN(GAMMA_ITERATION_FACTOR * iter , GAMMA_MAX_NB_ITERATION)); */

            // approximations Johnson, Kotz & Balakrishnan, Continuous Univariate Distributions, vol. 1, 2nd ed., pp. 361-362

//            process->observation[i]->shape = mean[i] / (2 * (mean[i] - exp(log_geometric_mean))) - 1./12.;
            diff = log(mean[i]) - log_geometric_mean;
            process->observation[i]->shape = (1 + sqrt(1 + 4 * diff / 3)) / (4 * diff);
            process->observation[i]->scale = mean[i] / process->observation[i]->shape;
          }
        }

        else {
          process->observation[i]->shape = GAMMA_MIN_SHAPE_PARAMETER;
          process->observation[i]->scale = GAMMA_DEFAULT_SCALE_PARAMETER;
        }
      }
    }

    else {
      process->observation[i]->zero_probability = D_DEFAULT;
    }
  }

  delete [] state_frequency;
  delete [] zero_mass;
  delete [] mean;
  delete [] variance;
}


/*--------------------------------------------------------------*
 *
 *  Estimation des parametres des lois gaussiennes d'un processus d'observation.
 *
 *  arguments : comptage des etats, variable, pointeur sur un processus d'observation continu.
 *
 *--------------------------------------------------------------*/

template <typename Type>
void MarkovianSequences::gaussian_estimation(Type ***state_sequence_count , int variable ,
                                             ContinuousParametricProcess *process) const

{
  register int i , j , k;
  double diff , *mean , *variance;
  Type *state_frequency;


  state_frequency = new Type[process->nb_state];
  mean = new double[process->nb_state];

  for (i = 0;i < process->nb_state;i++) {
    mean[i] = 0.;
    state_frequency[i] = 0;
  }

  switch (type[variable]) {

  case INT_VALUE : {
    for (i = 0;i < nb_sequence;i++) {
      for (j = 0;j < length[i];j++) {
        for (k = 0;k < process->nb_state;k++) {
          mean[k] += state_sequence_count[i][j][k] * int_sequence[i][variable][j];
          state_frequency[k] += state_sequence_count[i][j][k];
        }
      }
    }
    break;
  }

  case REAL_VALUE : {
    for (i = 0;i < nb_sequence;i++) {
      for (j = 0;j < length[i];j++) {
        for (k = 0;k < process->nb_state;k++) {
          mean[k] += state_sequence_count[i][j][k] * real_sequence[i][variable][j];
          state_frequency[k] += state_sequence_count[i][j][k];
        }
      }
    }
    break;
  }
  }

  for (i = 0;i < process->nb_state;i++) {
    if (state_frequency[i] > 0) {
      mean[i] /= state_frequency[i];
      process->observation[i]->location = mean[i];
    }
    else {
      process->observation[i]->location = D_INF;
    }
  }

  switch (process->tied_dispersion) {

  case false : {
    variance = new double[process->nb_state];
    for (i = 0;i < process->nb_state;i++) {
      variance[i] = 0.;
    }

    switch (type[variable]) {

    case INT_VALUE : {
      for (i = 0;i < nb_sequence;i++) {
        for (j = 0;j < length[i];j++) {
          for (k = 0;k < process->nb_state;k++) {
            diff = int_sequence[i][variable][j] - mean[k];
            variance[k] += state_sequence_count[i][j][k] * diff * diff;
          }
        }
      }
      break;
    }

    case REAL_VALUE : {
      for (i = 0;i < nb_sequence;i++) {
        for (j = 0;j < length[i];j++) {
          for (k = 0;k < process->nb_state;k++) {
            diff = real_sequence[i][variable][j] - mean[k];
            variance[k] += state_sequence_count[i][j][k] * diff * diff;
          }
        }
      }
      break;
    }
    }

    for (i = 0;i < process->nb_state;i++) {
      if (state_frequency[i] > 1) {
//        variance[i] /= state_frequency[i];
        variance[i] /= (state_frequency[i] - 1);
        process->observation[i]->dispersion = sqrt(variance[i]);
        if (process->observation[i]->dispersion / process->observation[i]->location < GAUSSIAN_MIN_VARIATION_COEFF) {
          process->observation[i]->dispersion = process->observation[i]->location * GAUSSIAN_MIN_VARIATION_COEFF;
        }
      }
    }
    break;
  }

  case true : {
    for (i = 1;i < process->nb_state;i++) {
      state_frequency[0] += state_frequency[i];
    }

    variance = new double[1];
    variance[0] = 0.;

    switch (type[variable]) {

    case INT_VALUE : {
      for (i = 0;i < nb_sequence;i++) {
        for (j = 0;j < length[i];j++) {
          for (k = 0;k < process->nb_state;k++) {
            diff = int_sequence[i][variable][j] - mean[k];
            variance[0] += state_sequence_count[i][j][k] * diff * diff;
          }
        }
      }
      break;
    }

    case REAL_VALUE : {
      for (i = 0;i < nb_sequence;i++) {
        for (j = 0;j < length[i];j++) {
          for (k = 0;k < process->nb_state;k++) {
            diff = real_sequence[i][variable][j] - mean[k];
            variance[0] += state_sequence_count[i][j][k] * diff * diff;
          }
        }
      }
      break;
    }
    }

//    variance[0] /= state_frequency[0];
    variance[0] /= (state_frequency[0] - process->nb_state);

    process->observation[0]->dispersion = sqrt(variance[0]);
    for (i = 1;i < process->nb_state;i++) {
      process->observation[i]->dispersion = process->observation[0]->dispersion;
    }
    break;
  }
  }

  delete [] state_frequency;
  delete [] mean;
  delete [] variance;
}


/*--------------------------------------------------------------*
 *
 *  Estimation des parametres des lois de von Mises d'un processus d'observation.
 *
 *  arguments : comptage des etats, variable, pointeur sur un processus d'observation continu.
 *
 *--------------------------------------------------------------*/

template <typename Type>
void MarkovianSequences::von_mises_estimation(Type ***state_sequence_count , int variable ,
                                              ContinuousParametricProcess *process) const

{
  register int i , j , k;
  double buff , global_mean_direction , concentration , **mean_direction;
  Type *state_frequency;


  state_frequency = new Type[process->nb_state];
  mean_direction = new double*[process->nb_state];
  for (i = 0;i < process->nb_state;i++) {
    mean_direction[i] = new double[4];
  }

  for (i = 0;i < process->nb_state;i++) {
    mean_direction[i][0] = 0.;
    mean_direction[i][1] = 0.;

    state_frequency[i] = 0;
  }

  switch (type[variable]) {

  case INT_VALUE : {
    for (i = 0;i < nb_sequence;i++) {
      for (j = 0;j < length[i];j++) {
        for (k = 0;k < process->nb_state;k++) {
          mean_direction[k][0] += state_sequence_count[i][j][k] * cos(int_sequence[i][variable][j] * M_PI / 180);
          mean_direction[k][1] += state_sequence_count[i][j][k] * sin(int_sequence[i][variable][j] * M_PI / 180);
          state_frequency[k] += state_sequence_count[i][j][k];
        }
      }
    }
    break;
  }

  case REAL_VALUE : {
    for (i = 0;i < nb_sequence;i++) {
      for (j = 0;j < length[i];j++) {
        for (k = 0;k < process->nb_state;k++) {
          switch (process->unit) {
          case DEGREE :
            mean_direction[k][0] += state_sequence_count[i][j][k] * cos(real_sequence[i][variable][j] * M_PI / 180);
            mean_direction[k][1] += state_sequence_count[i][j][k] * sin(real_sequence[i][variable][j] * M_PI / 180);
            break;
          case RADIAN :
            mean_direction[k][0] += state_sequence_count[i][j][k] * cos(real_sequence[i][variable][j]);
            mean_direction[k][1] += state_sequence_count[i][j][k] * sin(real_sequence[i][variable][j]);
            break;
          }

          state_frequency[k] += state_sequence_count[i][j][k];
        }
      }
    }
    break;
  }
  }

  for (i = 0;i < process->nb_state;i++) {
    if (state_frequency[i] > 0) {
      mean_direction[i][0] /= state_frequency[i];
      mean_direction[i][1] /= state_frequency[i];

      mean_direction[i][2] = sqrt(mean_direction[i][0] * mean_direction[i][0] +
                                  mean_direction[i][1] * mean_direction[i][1]);

      if (mean_direction[i][2] > 0.) {
        mean_direction[i][3] = atan(mean_direction[i][1] / mean_direction[i][0]);

        if (mean_direction[i][0] < 0.) {
          mean_direction[i][3] += M_PI;
        }
        else if (mean_direction[i][1] < 0.) {
          mean_direction[i][3] += 2 * M_PI;
        }

        if (process->unit == DEGREE) {
          mean_direction[i][3] *= (180 / M_PI);
        }
      }

      process->observation[i]->location = mean_direction[i][3];
    }

    else {
      process->observation[i]->location = D_INF;
    }
  }

  switch (process->tied_dispersion) {

  case false : {
    for (i = 0;i < process->nb_state;i++) {
      if (state_frequency[i] > 0) {
        process->observation[i]->dispersion = von_mises_concentration_computation(mean_direction[i][2]);
      }
      else {
        process->observation[i]->dispersion = D_DEFAULT;
      }
    }
    break;
  }

  case true : {
    global_mean_direction = 0.;
    buff = 0.;

    for (i = 0;i < process->nb_state;i++) {
      global_mean_direction += state_frequency[i] * mean_direction[i][2];
      buff += state_frequency[i];
    }
    concentration = von_mises_concentration_computation(global_mean_direction / buff);

    for (i = 0;i < process->nb_state;i++) {
      process->observation[i]->dispersion = concentration;
    }
    break;
  }
  }

  for (i = 0;i < process->nb_state;i++) {
    delete [] mean_direction[i];
  }
  delete [] mean_direction;

  delete [] state_frequency;
}


/*--------------------------------------------------------------*
 *
 *  Estimation des parametres des modeles lineaires gaussiens de tendance
 *  d'un processus d'observation.
 *
 *  arguments : comptage des etats, variable, pointeur sur un processus d'observation continu.
 *
 *--------------------------------------------------------------*/

template <typename Type>
void MarkovianSequences::linear_model_estimation(Type ***state_sequence_count , int variable ,
                                                 ContinuousParametricProcess *process) const

{
  register int i , j , k;
  double diff , *mean , *index_parameter_mean , *variance , *index_parameter_variance , *covariance ,
         *residual_mean , *residual_square_sum , *residual_variance;
  Type *state_frequency;


  state_frequency = new Type[process->nb_state];
  mean = new double[process->nb_state];
  index_parameter_mean = new double[process->nb_state];

  for (i = 0;i < process->nb_state;i++) {
    mean[i] = 0.;
    index_parameter_mean[i] = 0.;
    state_frequency[i] = 0;
  }

  switch (type[variable]) {

  case INT_VALUE : {
    for (i = 0;i < nb_sequence;i++) {
      for (j = 0;j < length[i];j++) {
        for (k = 0;k < process->nb_state;k++) {
          mean[k] += state_sequence_count[i][j][k] * int_sequence[i][variable][j];
          state_frequency[k] += state_sequence_count[i][j][k];
        }
      }
    }
    break;
  }

  case REAL_VALUE : {
    for (i = 0;i < nb_sequence;i++) {
      for (j = 0;j < length[i];j++) {
        for (k = 0;k < process->nb_state;k++) {
          mean[k] += state_sequence_count[i][j][k] * real_sequence[i][variable][j];
          state_frequency[k] += state_sequence_count[i][j][k];
        }
      }
    }
    break;
  }
  }

  switch (index_parameter_type) {

  case IMPLICIT_TYPE : {
    for (i = 0;i < nb_sequence;i++) {
      for (j = 0;j < length[i];j++) {
        for (k = 0;k < process->nb_state;k++) {
          index_parameter_mean[k] += state_sequence_count[i][j][k] * j;
        }
      }
    }
    break;
  }

  case TIME : {
    for (i = 0;i < nb_sequence;i++) {
      for (j = 0;j < length[i];j++) {
        for (k = 0;k < process->nb_state;k++) {
          index_parameter_mean[k] += state_sequence_count[i][j][k] * index_parameter[i][j];
        }
      }
    }
    break;
  }
  }

  for (i = 0;i < process->nb_state;i++) {
    if (state_frequency[i] > 0) {
      mean[i] /= state_frequency[i];
      index_parameter_mean[i] /= state_frequency[i];
    }
  }

  variance = new double[process->nb_state];
  index_parameter_variance = new double[process->nb_state];
  covariance = new double[process->nb_state];
  for (i = 0;i < process->nb_state;i++) {
    variance[i] = 0.;
    index_parameter_variance[i] = 0.;
    covariance[i] = 0.;
  }

  switch (type[variable]) {

  case INT_VALUE : {
    for (i = 0;i < nb_sequence;i++) {
      for (j = 0;j < length[i];j++) {
        for (k = 0;k < process->nb_state;k++) {
          diff = int_sequence[i][variable][j] - mean[k];
          variance[k] += state_sequence_count[i][j][k] * diff * diff;
        }
      }
    }
    break;
  }

  case REAL_VALUE : {
    for (i = 0;i < nb_sequence;i++) {
      for (j = 0;j < length[i];j++) {
        for (k = 0;k < process->nb_state;k++) {
          diff = real_sequence[i][variable][j] - mean[k];
          variance[k] += state_sequence_count[i][j][k] * diff * diff;
        }
      }
    }
    break;
  }
  }

  switch (index_parameter_type) {

  case IMPLICIT_TYPE : {
    for (i = 0;i < nb_sequence;i++) {
      for (j = 0;j < length[i];j++) {
        for (k = 0;k < process->nb_state;k++) {
          diff = j - index_parameter_mean[k];
          index_parameter_variance[k] += state_sequence_count[i][j][k] * diff * diff;
        }
      }
    }
    break;
  }

  case TIME : {
    for (i = 0;i < nb_sequence;i++) {
      for (j = 0;j < length[i];j++) {
        for (k = 0;k < process->nb_state;k++) {
          diff = index_parameter[i][j] - index_parameter_mean[k];
          index_parameter_variance[k] += state_sequence_count[i][j][k] * diff * diff;
        }
      }
    }
    break;
  }
  }

  switch (type[variable]) {

  case INT_VALUE : {
    switch (index_parameter_type) {

    case IMPLICIT_TYPE : {
      for (i = 0;i < nb_sequence;i++) {
        for (j = 0;j < length[i];j++) {
          for (k = 0;k < process->nb_state;k++) {
            covariance[k] += state_sequence_count[i][j][k] * (int_sequence[i][variable][j] - mean[k]) *
                             (j - index_parameter_mean[k]);
          }
        }
      }
      break;
    }

    case TIME : {
      for (i = 0;i < nb_sequence;i++) {
        for (j = 0;j < length[i];j++) {
          for (k = 0;k < process->nb_state;k++) {
            covariance[k] += state_sequence_count[i][j][k] * (int_sequence[i][variable][j] - mean[k]) *
                             (index_parameter[i][j] - index_parameter_mean[k]);
          }
        }
      }
      break;
    }
    }
    break;
  }

  case REAL_VALUE : {
    switch (index_parameter_type) {

    case IMPLICIT_TYPE : {
      for (i = 0;i < nb_sequence;i++) {
        for (j = 0;j < length[i];j++) {
          for (k = 0;k < process->nb_state;k++) {
            covariance[k] += state_sequence_count[i][j][k] * (real_sequence[i][variable][j] - mean[k]) *
                             (j - index_parameter_mean[k]);
          }
        }
      }
      break;
    }

    case TIME : {
      for (i = 0;i < nb_sequence;i++) {
        for (j = 0;j < length[i];j++) {
          for (k = 0;k < process->nb_state;k++) {
            covariance[k] += state_sequence_count[i][j][k] * (real_sequence[i][variable][j] - mean[k]) *
                             (index_parameter[i][j] - index_parameter_mean[k]);
          }
        }
      }
      break;
    }
    }
    break;
  }
  }

  for (i = 0;i < process->nb_state;i++) {
    if (state_frequency[i] > 0) {
      process->observation[i]->slope = covariance[i] / index_parameter_variance[i];
      process->observation[i]->intercept = mean[i] - process->observation[i]->slope * index_parameter_mean[i];
      process->observation[i]->correlation = covariance[i] / sqrt(variance[i] * index_parameter_variance[i]);
    }
    else {
      process->observation[i]->slope = D_INF;
      process->observation[i]->intercept = D_INF;
    }
  }

  residual_mean = new double[process->nb_state];
  residual_square_sum = new double[process->nb_state];
  for (i = 0;i < process->nb_state;i++) {
    residual_mean[i] = 0.;
    residual_square_sum[i] = 0.;
  }

  switch (type[variable]) {

  case INT_VALUE : {
    switch (index_parameter_type) {

    case IMPLICIT_TYPE : {
      for (i = 0;i < nb_sequence;i++) {
        for (j = 0;j < length[i];j++) {
          for (k = 0;k < process->nb_state;k++) {
            diff = state_sequence_count[i][j][k] * (int_sequence[i][variable][j] -
                    (process->observation[k]->intercept + process->observation[k]->slope * j));
            residual_mean[k] += diff;
            residual_square_sum[k] += diff * diff;
          }
        }
      }
      break;
    }

    case TIME : {
      for (i = 0;i < nb_sequence;i++) {
        for (j = 0;j < length[i];j++) {
          for (k = 0;k < process->nb_state;k++) {
            diff = state_sequence_count[i][j][k] * (int_sequence[i][variable][j] -
                    (process->observation[k]->intercept + process->observation[k]->slope * index_parameter[i][j]));
            residual_mean[k] += diff;
            residual_square_sum[k] += diff * diff;
          }
        }
      }
      break;
    }
    }
    break;
  }

  case REAL_VALUE : {
    switch (index_parameter_type) {

    case IMPLICIT_TYPE : {
      for (i = 0;i < nb_sequence;i++) {
        for (j = 0;j < length[i];j++) {
          for (k = 0;k < process->nb_state;k++) {
            diff = state_sequence_count[i][j][k] * (real_sequence[i][variable][j] -
                    (process->observation[k]->intercept + process->observation[k]->slope * j));
            residual_mean[k] += diff;
            residual_square_sum[k] += diff * diff;
          }
        }
      }
      break;
    }

    case TIME : {
      for (i = 0;i < nb_sequence;i++) {
        for (j = 0;j < length[i];j++) {
          for (k = 0;k < process->nb_state;k++) {
            diff = state_sequence_count[i][j][k] * (real_sequence[i][variable][j] -
                    (process->observation[k]->intercept + process->observation[k]->slope * index_parameter[i][j]));
            residual_mean[k] += diff;
            residual_square_sum[k] += diff * diff;
          }
        }
      }
      break;
    }
    }
    break;
  }
  }

  for (i = 0;i < process->nb_state;i++) {
    if (state_frequency[i] > 0) {
      residual_mean[i] /= state_frequency[i];
    }

    if (state_frequency[i] > 2) {
      residual_square_sum[i] /= (state_frequency[i] - 2);
      process->observation[i]->slope_standard_deviation = sqrt(residual_square_sum[i] / index_parameter_variance[i]);
      process->observation[i]->sample_size = state_frequency[i] - 2;
    }
    else {
      process->observation[i]->slope_standard_deviation = 0;
      process->observation[i]->sample_size = 0;
    }
  }

  residual_variance = new double[process->nb_state];
  for (i = 0;i < process->nb_state;i++) {
    residual_variance[i] = 0.;
  }

  switch (type[variable]) {

  case INT_VALUE : {
    switch (index_parameter_type) {

    case IMPLICIT_TYPE : {
      for (i = 0;i < nb_sequence;i++) {
        for (j = 0;j < length[i];j++) {
          for (k = 0;k < process->nb_state;k++) {
            diff = int_sequence[i][variable][j] - (process->observation[k]->intercept +
                    process->observation[k]->slope * j) - residual_mean[k];
            residual_variance[k] += state_sequence_count[i][j][k] * diff * diff;
          }
        }
      }
      break;
    }

    case TIME : {
      for (i = 0;i < nb_sequence;i++) {
        for (j = 0;j < length[i];j++) {
          for (k = 0;k < process->nb_state;k++) {
            diff = int_sequence[i][variable][j] - (process->observation[k]->intercept +
                    process->observation[k]->slope * index_parameter[i][j]) - residual_mean[k];
            residual_variance[k] += state_sequence_count[i][j][k] * diff * diff;
          }
        }
      }
      break;
    }
    }
    break;
  }

  case REAL_VALUE : {
    switch (index_parameter_type) {

    case IMPLICIT_TYPE : {
      for (i = 0;i < nb_sequence;i++) {
        for (j = 0;j < length[i];j++) {
          for (k = 0;k < process->nb_state;k++) {
            diff = real_sequence[i][variable][j] - (process->observation[k]->intercept +
                    process->observation[k]->slope * j) - residual_mean[k];
            residual_variance[k] += state_sequence_count[i][j][k] * diff * diff;
          }
        }
      }
      break;
    }

    case TIME : {
      for (i = 0;i < nb_sequence;i++) {
        for (j = 0;j < length[i];j++) {
          for (k = 0;k < process->nb_state;k++) {
            diff = real_sequence[i][variable][j] - (process->observation[k]->intercept +
                    process->observation[k]->slope * index_parameter[i][j]) - residual_mean[k];
            residual_variance[k] += state_sequence_count[i][j][k] * diff * diff;
          }
        }
      }
      break;
    }
    }
    break;
  }
  }

  for (i = 0;i < process->nb_state;i++) {
    if (state_frequency[i] > 2) {
//      process->observation[i]->dispersion = sqrt(residual_variance[i] / state_frequency[i]);
      process->observation[i]->dispersion = sqrt(residual_variance[i] / (state_frequency[i] - 2));
      if (process->observation[i]->dispersion / mean[i] < GAUSSIAN_MIN_VARIATION_COEFF) {
        process->observation[i]->dispersion = mean[i] * GAUSSIAN_MIN_VARIATION_COEFF;
      }
    }

#   ifdef DEBUG
    cout << "\n" << STAT_label[STATL_VARIABLE] << " " << variable + 1 << "   "
         << STAT_word[STATW_STATE] << " " << i << "   "
         << STAT_word[STATW_INTERCEPT] << " : " << process->observation[i]->intercept << "   "
         << STAT_word[STATW_SLOPE] << " : " << process->observation[i]->slope << "   "
         << STAT_word[STATW_STANDARD_DEVIATION] << " : " << process->observation[i]->dispersion << endl;
#   endif

  }

  delete [] state_frequency;
  delete [] mean;
  delete [] index_parameter_mean;
  delete [] variance;
  delete [] index_parameter_variance;
  delete [] covariance;
  delete [] residual_mean;
  delete [] residual_square_sum;
  delete [] residual_variance;
}


};  // namespace sequence_analysis



#endif
