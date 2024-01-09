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
 *       $Id: continuous_distribution_estimation.hpp 12646 2012-08-03 08:12:47Z guedon $
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



#ifndef CONTINUOUS_PARAMETRIC_ESTIMATION_HPP
#define CONTINUOUS_PARAMETRIC_ESTIMATION_HPP


#include <math.h>

#include <boost/math/special_functions/digamma.hpp>

#include "stat_tool/stat_tools.h"
#include "stat_tool/curves.h"
#include "stat_tool/markovian.h"
#include "stat_tool/vectors.h"
#include "stat_tool/stat_label.h"

using namespace std;
using namespace boost::math;


namespace stat_tool {


extern double von_mises_concentration_computation(double mean_direction);



/*--------------------------------------------------------------*
 *
 *  Estimation des parametres des lois gamma d'un processus d'observation.
 *
 *  arguments : comptage des etats, variable, pointeur sur un processus d'observation continu,
 *              iteration EM.
 *
 *--------------------------------------------------------------*/

template <typename Type>
void Vectors::gamma_estimation(Type **component_vector_count , int variable ,
                               ContinuousParametricProcess *process , int iter) const

{
  register int i , j;
  double diff , log_geometric_mean , *zero_mass , *mean , *variance;
  Type *component_frequency;


  component_frequency = new Type[process->nb_state];
  zero_mass = new double[process->nb_state];
  mean = new double[process->nb_state];
  variance = new double[process->nb_state];

  for (i = 0;i < process->nb_state;i++) {
    zero_mass[i] = 0.;
    mean[i] = 0.;
    component_frequency[i] = 0;
  }

  switch (type[variable]) {

  case INT_VALUE : {
    for (i = 0;i < nb_vector;i++) {
      for (j = 0;j < process->nb_state;j++) {
        if (int_vector[i][variable] == 0) {
          zero_mass[j] += component_vector_count[i][j];
        }
        else {
          mean[j] += component_vector_count[i][j] * int_vector[i][variable];
        }
        component_frequency[j] += component_vector_count[i][j];
      }
    }
    break;
  }

  case REAL_VALUE : {
    for (i = 0;i < nb_vector;i++) {
      for (j = 0;j < process->nb_state;j++) {
        if (real_vector[i][variable] == 0.) {
          zero_mass[j] += component_vector_count[i][j];
        }
        else {
          mean[j] += component_vector_count[i][j] * real_vector[i][variable];
        }
        component_frequency[j] += component_vector_count[i][j];
      }
    }
    break;
  }
  }

  for (i = 0;i < process->nb_state;i++) {
    if (component_frequency[i] > 0) {
      mean[i] /= component_frequency[i];
    }
    variance[i] = 0.;
  }

  switch (type[variable]) {

  case INT_VALUE : {
    for (i = 0;i < nb_vector;i++) {
      for (j = 0;j < process->nb_state;j++) {
        diff = int_vector[i][variable] - mean[j];
        variance[j] += component_vector_count[i][j] * diff * diff;
      }
    }
    break;
  }

  case REAL_VALUE : {
    for (i = 0;i < nb_vector;i++) {
      for (j = 0;j < process->nb_state;j++) {
        diff = real_vector[i][variable] - mean[j];
        variance[j] += component_vector_count[i][j] * diff * diff;
      }
    }
    break;
  }
  }

  for (i = 0;i < process->nb_state;i++) {
    if (component_frequency[i] > 0) {
      variance[i] /= component_frequency[i];
//      variance[i] /= (component_frequency[i] - 1);
    }
  }

  for (i = 0;i < process->nb_state;i++) {

#   ifdef DEBUG
    if ((iter >= 5) && (component_frequency[i] > 0)) {
      cout << "\n" << STAT_word[STATW_COMPONENT] << " " << i << " : "
           << zero_mass[i] << ", " << component_frequency[i] << " | "
           << zero_mass[i] / component_frequency[i] << endl;
    }
#   endif

    if (component_frequency[i] > 0) {
      if (zero_mass[i] / component_frequency[i] > GAMMA_ZERO_FREQUENCY_THRESHOLD) {
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

          process->observation[i]->shape = mean[i] * mean[i] / variance[i] - 1. / (double)component_frequency[i];
          process->observation[i]->scale = mean[i] / process->observation[i]->shape;

#         ifdef DEBUG    // essai pour eviter les très petits parametres de forme
          if ((process->observation[i]->shape < 1.) && (mean[i] < 5.)) {
            process->observation[i]->shape = 1.;
            process->observation[i]->scale = mean[i] / process->observation[i]->shape;
          }
#         endif

          if ((process->observation[i]->shape >= GAMMA_SHAPE_PARAMETER_THRESHOLD) &&
              (component_frequency[i] < GAMMA_FREQUENCY_THRESHOLD)) {
            log_geometric_mean = 0.;
            component_frequency[i] = 0;

            switch (type[variable]) {

            case INT_VALUE : {
              for (j = 0;j < nb_vector;j++) {
                if (int_vector[j][variable] > 0) {
                  log_geometric_mean += component_vector_count[j][i] * log(int_vector[j][variable]);
                  component_frequency[i] += component_vector_count[j][i];
                }
              }
              break;
            }

            case REAL_VALUE : {
              for (j = 0;j < nb_vector;j++) {
                if (real_vector[j][variable] > 0.) {
                  log_geometric_mean += component_vector_count[j][i] * log(real_vector[j][variable]);
                  component_frequency[i] += component_vector_count[j][i];
                }
              }
              break;
            }
            }

            log_geometric_mean /= component_frequency[i];
/*            j = 0;   a revoir

#           ifdef DEBUG
            cout << "\n" << STAT_word[STATW_COMPONENT] << " " << i << "   "
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

  delete [] component_frequency;
  delete [] zero_mass;
  delete [] mean;
  delete [] variance;
}


/*--------------------------------------------------------------*
 *
 *  Estimation des parametres des lois gamma d'un processus d'observation :
 *  parametres de forme et eventuellement parametres d'echelle lies.
 *
 *  arguments : comptage des etats, variable, pointeur sur un processus d'observation continue,
 *              facteur pour la variance, iteration EM.
 *
 *--------------------------------------------------------------*/

template <typename Type>
void Vectors::tied_gamma_estimation(Type **component_vector_count , int variable ,
                                    ContinuousParametricProcess *process ,
                                    int variance_factor , int iter) const

{
  register int i , j;
  double diff , variance , log_geometric_mean , *mean , *factor;
  Type *component_frequency;


  component_frequency = new Type[process->nb_state];
  mean = new double[process->nb_state];
  factor = new double[process->nb_state];

  switch (process->tied_location) {

  case false : {
    for (i = 0;i < process->nb_state;i++) {
      mean[i] = 0.;
      component_frequency[i] = 0;
    }

    switch (type[variable]) {

    case INT_VALUE : {
      for (i = 0;i < nb_vector;i++) {
        for (j = 0;j < process->nb_state;j++) {
          mean[j] += component_vector_count[i][j] * int_vector[i][variable];
          component_frequency[j] += component_vector_count[i][j];
        }
      }
      break;
    }

    case REAL_VALUE : {
      for (i = 0;i < nb_vector;i++) {
        for (j = 0;j < process->nb_state;j++) {
          mean[j] += component_vector_count[i][j] * real_vector[i][variable];
          component_frequency[j] += component_vector_count[i][j];
        }
      }
      break;
    }
    }

    for (i = 0;i < process->nb_state;i++) {
      if (component_frequency[i] > 0) {
        mean[i] /= component_frequency[i];
      }
    }

    switch (variance_factor) {

    case CONVOLUTION_FACTOR : {
      factor[0] = 1.;
      for (i = 1;i < process->nb_state;i++) {
        if (component_frequency[i] > 0) {
          factor[i] = mean[i] / mean[0];
        }
        else {
          factor[i] = 1.;
        }
      }
      break;
    }

    case SCALING_FACTOR : {
      factor[0] = 1.;
      for (i = 1;i < process->nb_state;i++) {
        if (component_frequency[i] > 0) {
          factor[i] = (mean[i] * mean[i]) / (mean[0] * mean[0]);
        }
        else {
          factor[i] = 1.;
        }
      }
      break;
    }
    }

    for (i = 1;i < process->nb_state;i++) {
      component_frequency[0] += component_frequency[i];
    }
    break;
  }

  case true : {
    mean[0] = 0.;
    component_frequency[0] = 0.;

    factor[0] = 1.;
    for (i = 1;i < process->nb_state;i++) {
      factor[i] = factor[i - 1] * 2;
    }

    switch (type[variable]) {

    case INT_VALUE : {
      for (i = 0;i < nb_vector;i++) {
        for (j = 0;j < process->nb_state;j++) {
          mean[0] += component_vector_count[i][j] * int_vector[i][variable] / factor[j];
          component_frequency[0] += component_vector_count[i][j];
        }
      }
      break;
    }

    case REAL_VALUE : {
      for (i = 0;i < nb_vector;i++) {
        for (j = 0;j < process->nb_state;j++) {
          mean[0] += component_vector_count[i][j] * real_vector[i][variable] / factor[j];
          component_frequency[0] += component_vector_count[i][j];
        }
      }
      break;
    }
    }

    mean[0] /= component_frequency[0];
    for (i = 1;i < process->nb_state;i++) {
      mean[i] = mean[0] * factor[i];
    }

    if (variance_factor == SCALING_FACTOR) {
      for (i = 1;i < process->nb_state;i++) {
        factor[i] *= factor[i];
//        factor[i] = factor[i - 1] * 4;
      }
    }
    break;
  }
  }

  variance = 0.;

  switch (type[variable]) {

  case INT_VALUE : {
    for (i = 0;i < nb_vector;i++) {
      for (j = 0;j < process->nb_state;j++) {
        diff = int_vector[i][variable] - mean[j];
        variance += component_vector_count[i][j] * diff * diff / factor[j];
      }
    }
    break;
  }

  case REAL_VALUE : {
    for (i = 0;i < nb_vector;i++) {
      for (j = 0;j < process->nb_state;j++) {
        diff = real_vector[i][variable] - mean[j];
        variance += component_vector_count[i][j] * diff * diff / factor[j];
      }
    }
    break;
  }
  }

//  variance /= component_frequency[0];
  variance /= (component_frequency[0] - process->nb_state);

  switch (variance_factor) {

  case CONVOLUTION_FACTOR : {
    process->observation[0]->shape = mean[0] * mean[0] / variance;
    process->observation[0]->scale = variance / mean[0];

    for (i = 1;i < process->nb_state;i++) {
      process->observation[i]->shape = process->observation[0]->shape * factor[i];
      process->observation[i]->scale = process->observation[0]->scale;
    }
    break;
  }

  case SCALING_FACTOR : {
    switch (process->tied_location) {

    case false : {
//      process->observation[0]->shape = mean[0] * mean[0] / variance;
//      process->observation[0]->scale = variance / mean[0];

      // Hawang & Huang (2012), Ann. Inst. Statist. Math. 54(4), 840-847

      process->observation[0]->shape = mean[0] * mean[0] / variance - 1. / (double)component_frequency[0];
      process->observation[0]->scale = mean[0] / process->observation[0]->shape;

      for (i = 1;i < process->nb_state;i++) {
        process->observation[i]->shape = process->observation[0]->shape;
//        factor[i] = sqrt(factor[i]);
        factor[i] = mean[i] / mean[0];
        process->observation[i]->scale = process->observation[0]->scale * factor[i];
      }
      break;
    }

    case true : {
      process->observation[0]->shape = mean[0] * mean[0] / variance;
      process->observation[0]->scale = variance / mean[0];

      for (i = 1;i < process->nb_state;i++) {
        process->observation[i]->shape = process->observation[0]->shape;
//        factor[i] = sqrt(factor[i]);
        factor[i] = factor[i - 1] * 2;
        process->observation[i]->scale = process->observation[0]->scale * factor[i];
      }
      break;
    }
    }
    break;
  }
  }

  if ((process->observation[0]->shape >= GAMMA_SHAPE_PARAMETER_THRESHOLD) &&
      (component_frequency[0] < GAMMA_FREQUENCY_THRESHOLD)) {
    log_geometric_mean = 0.;
//    component_frequency[0] = 0.;
    for (i = 0;i < process->nb_state;i++) {
      component_frequency[i] = 0;
    }

    switch (type[variable]) {

    case INT_VALUE : {
      for (i = 0;i < nb_vector;i++) {
        for (j = 0;j < process->nb_state;j++) {
          if (int_vector[i][variable] > 0) {
//            log_geometric_mean += component_vector_count[i][j] * log(int_vector[i][variable] / factor[j]);
//            component_frequency[0] += component_vector_count[i][j];
            log_geometric_mean += component_vector_count[i][j] * log(int_vector[i][variable]);
            component_frequency[j] += component_vector_count[i][j];
          }
        }
      }
      break;
    }

    case REAL_VALUE : {
      for (i = 0;i < nb_vector;i++) {
        for (j = 0;j < process->nb_state;j++) {
          if (real_vector[i][variable] > 0.) {
//            log_geometric_mean += component_vector_count[i][j] * log(real_vector[i][variable] / factor[j]);
//            component_frequency[0] += component_vector_count[i][j];
            log_geometric_mean += component_vector_count[i][j] * log(real_vector[i][variable]);
            component_frequency[j] += component_vector_count[i][j];
          }
        }
      }
      break;
    }
    }

    for (i = 1;i < process->nb_state;i++) {
      log_geometric_mean -= component_frequency[i] * log(factor[i]);
      component_frequency[0] += component_frequency[i];
    }
    log_geometric_mean /= component_frequency[0];

#   ifdef MESSAGE
    cout << "\n" << STAT_word[STATW_SHAPE] << " : " << process->observation[0]->shape << "   "
         << STAT_word[STATW_SCALE] << " : " << process->observation[0]->scale << endl;
#   endif

    // approximations Johnson, Kotz & Balakrishnan, Continuous Univariate Distributions, vol. 1, 2nd ed., pp. 361-362

//    process->observation[0]->shape = mean[0] / (2 * (mean[0] - exp(log_geometric_mean))) - 1./12.;
    diff = log(mean[0]) - log_geometric_mean;
    process->observation[0]->shape = (1 + sqrt(1 + 4 * diff / 3)) / (4 * diff);

#   ifdef DEBUG
    cout << STAT_word[STATW_SHAPE] << " : " << process->observation[0]->shape << endl;
#   endif

    i = 0;

/*    do {
      process->observation[0]->scale = exp(log_geometric_mean - digamma(process->observation[0]->shape));
      process->observation[0]->shape = mean[0] / process->observation[0]->scale;
      i++;
    }
    while (i < MIN(GAMMA_ITERATION_FACTOR * iter , GAMMA_MAX_NB_ITERATION)); */

    do {
      process->observation[0]->shape = process->observation[0]->shape * (log(process->observation[0]->shape) -
                                        digamma(process->observation[0]->shape)) / diff;
      i++;

#     ifdef DEBUG
      cout << STAT_word[STATW_SHAPE] << " : " << process->observation[0]->shape << endl;
#     endif

    }
    while (i < MIN(GAMMA_ITERATION_FACTOR * iter , GAMMA_MAX_NB_ITERATION));

    process->observation[0]->scale = mean[0] / process->observation[0]->shape;

#   ifdef DEBUG
    cout << "\n" << STAT_word[STATW_SHAPE] << " : " << process->observation[0]->shape << "   "
         << STAT_word[STATW_SCALE] << " : " << process->observation[0]->scale << endl;
#   endif

    switch (variance_factor) {

    case CONVOLUTION_FACTOR : {
      for (i = 1;i < process->nb_state;i++) {
        process->observation[i]->shape = process->observation[0]->shape * factor[i];
        process->observation[i]->scale = process->observation[0]->scale;
      }
      break;
    }

    case SCALING_FACTOR : {
      for (i = 1;i < process->nb_state;i++) {
        process->observation[i]->shape = process->observation[0]->shape;
        process->observation[i]->scale = process->observation[0]->scale * factor[i];
      }
      break;
    }
    }
  }

  delete [] component_frequency;
  delete [] mean;
  delete [] factor;
}


/*--------------------------------------------------------------*
 *
 *  Estimation des parametres des lois gaussiennes d'un processus d'observation.
 *
 *  arguments : comptage des etats, variable, pointeur sur un processus d'observation continu.
 *
 *--------------------------------------------------------------*/

template <typename Type>
void Vectors::gaussian_estimation(Type **component_vector_count , int variable ,
                                  ContinuousParametricProcess *process) const

{
  register int i , j;
  double diff , *mean , *variance;
  Type *component_frequency;


  component_frequency = new Type[process->nb_state];
  mean = new double[process->nb_state];

  for (i = 0;i < process->nb_state;i++) {
    mean[i] = 0.;
    component_frequency[i] = 0;
  }

  switch (type[variable]) {

  case INT_VALUE : {
    for (i = 0;i < nb_vector;i++) {
      for (j = 0;j < process->nb_state;j++) {
        mean[j] += component_vector_count[i][j] * int_vector[i][variable];
        component_frequency[j] += component_vector_count[i][j];
      }
    }
    break;
  }

  case REAL_VALUE : {
    for (i = 0;i < nb_vector;i++) {
      for (j = 0;j < process->nb_state;j++) {
        mean[j] += component_vector_count[i][j] * real_vector[i][variable];
        component_frequency[j] += component_vector_count[i][j];
      }
    }
    break;
  }
  }

  for (i = 0;i < process->nb_state;i++) {
    if (component_frequency[i] > 0) {
      mean[i] /= component_frequency[i];
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
      for (i = 0;i < nb_vector;i++) {
        for (j = 0;j < process->nb_state;j++) {
          diff = int_vector[i][variable] - mean[j];
          variance[j] += component_vector_count[i][j] * diff * diff;
        }
      }
      break;
    }

    case REAL_VALUE : {
      for (i = 0;i < nb_vector;i++) {
        for (j = 0;j < process->nb_state;j++) {
          diff = real_vector[i][variable] - mean[j];
          variance[j] += component_vector_count[i][j] * diff * diff;
        }
      }
      break;
    }
    }

    for (i = 0;i < process->nb_state;i++) {
      if (component_frequency[i] > 1) {
//        variance[i] /= component_frequency[i];
        variance[i] /= (component_frequency[i] - 1);
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
      component_frequency[0] += component_frequency[i];
    }

    variance = new double[1];
    variance[0] = 0.;

    switch (type[variable]) {

    case INT_VALUE : {
      for (i = 0;i < nb_vector;i++) {
        for (j = 0;j < process->nb_state;j++) {
          diff = int_vector[i][variable] - mean[j];
          variance[0] += component_vector_count[i][j] * diff * diff;
        }
      }
      break;
    }

    case REAL_VALUE : {
      for (i = 0;i < nb_vector;i++) {
        for (j = 0;j < process->nb_state;j++) {
          diff = real_vector[i][variable] - mean[j];
          variance[0] += component_vector_count[i][j] * diff * diff;
        }
      }
      break;
    }
    }

//    variance[0] /= component_frequency[0];
    variance[0] /= (component_frequency[0] - process->nb_state);

    process->observation[0]->dispersion = sqrt(variance[0]);
    for (i = 1;i < process->nb_state;i++) {
      process->observation[i]->dispersion = process->observation[0]->dispersion;
    }
    break;
  }
  }

  delete [] component_frequency;
  delete [] mean;
  delete [] variance;
}


/*--------------------------------------------------------------*
 *
 *  Estimation des parametres des lois gaussiennes d'un processus d'observation :
 *  variance et eventuellement moyennes liees.
 *
 *  arguments : comptage des etats, variable, pointeur sur un processus d'observation continu,
 *              facteur pour la variance.
 *
 *--------------------------------------------------------------*/

template <typename Type>
void Vectors::tied_gaussian_estimation(Type **component_vector_count , int variable ,
                                       ContinuousParametricProcess *process ,
                                       int variance_factor) const

{
  register int i , j;
  double diff , variance , *mean , *factor;
  Type *component_frequency;


  component_frequency = new Type[process->nb_state];
  mean = new double[process->nb_state];
  factor = new double[process->nb_state];

  switch (process->tied_location) {

  case false : {
    for (i = 0;i < process->nb_state;i++) {
      mean[i] = 0.;
      component_frequency[i] = 0;
    }

    switch (type[variable]) {

    case INT_VALUE : {
      for (i = 0;i < nb_vector;i++) {
        for (j = 0;j < process->nb_state;j++) {
          mean[j] += component_vector_count[i][j] * int_vector[i][variable];
          component_frequency[j] += component_vector_count[i][j];
        }
      }
      break;
    }

    case REAL_VALUE : {
      for (i = 0;i < nb_vector;i++) {
        for (j = 0;j < process->nb_state;j++) {
          mean[j] += component_vector_count[i][j] * real_vector[i][variable];
          component_frequency[j] += component_vector_count[i][j];
        }
      }
      break;
    }
    }

    for (i = 0;i < process->nb_state;i++) {
      if (component_frequency[i] > 0) {
        mean[i] /= component_frequency[i];
        process->observation[i]->location = mean[i];
      }
      else {
        process->observation[i]->location = D_INF;
      }
    }

    switch (variance_factor) {

    case CONVOLUTION_FACTOR : {
      factor[0] = 1.;
      for (i = 1;i < process->nb_state;i++) {
        if (component_frequency[i] > 0) {
          factor[i] = mean[i] / mean[0];
        }
        else {
          factor[i] = 1.;
        }
      }
      break;
    }

    case SCALING_FACTOR : {
      factor[0] = 1.;
      for (i = 1;i < process->nb_state;i++) {
        if (component_frequency[i] > 0) {
          factor[i] = (mean[i] * mean[i]) / (mean[0] * mean[0]);
        }
        else {
          factor[i] = 1.;
        }
      }
      break;
    }
    }

    for (i = 1;i < process->nb_state;i++) {
      component_frequency[0] += component_frequency[i];
    }
    break;
  }

  case true : {
    mean[0] = 0.;
    component_frequency[0] = 0.;

    factor[0] = 1.;
    for (i = 1;i < process->nb_state;i++) {
      factor[i] = factor[i - 1] * 2;
    }

    switch (type[variable]) {

    case INT_VALUE : {
      for (i = 0;i < nb_vector;i++) {
        for (j = 0;j < process->nb_state;j++) {
          mean[0] += component_vector_count[i][j] * int_vector[i][variable] / factor[j];
          component_frequency[0] += component_vector_count[i][j];
        }
      }
      break;
    }

    case REAL_VALUE : {
      for (i = 0;i < nb_vector;i++) {
        for (j = 0;j < process->nb_state;j++) {
          mean[0] += component_vector_count[i][j] * real_vector[i][variable] / factor[j];
          component_frequency[0] += component_vector_count[i][j];
        }
      }
      break;
    }
    }

    mean[0] /= component_frequency[0];
    process->observation[0]->location = mean[0];
    for (i = 1;i < process->nb_state;i++) {
      mean[i] = mean[0] * factor[i];
      process->observation[i]->location = mean[i];
    }

    if (variance_factor == SCALING_FACTOR) {
      for (i = 1;i < process->nb_state;i++) {
//        factor[i] = factor[i - 1] * 4;
        factor[i] *= factor[i];
      }
    }
    break;
  }
  }

  variance = 0.;

  switch (type[variable]) {

  case INT_VALUE : {
    for (i = 0;i < nb_vector;i++) {
      for (j = 0;j < process->nb_state;j++) {
        diff = int_vector[i][variable] - mean[j];
        variance += component_vector_count[i][j] * diff * diff / factor[j];
      }
    }
    break;
  }

  case REAL_VALUE : {
    for (i = 0;i < nb_vector;i++) {
      for (j = 0;j < process->nb_state;j++) {
        diff = real_vector[i][variable] - mean[j];
        variance += component_vector_count[i][j] * diff * diff / factor[j];
      }
    }
    break;
  }
  }

//  variance /= component_frequency[0];
  variance /= (component_frequency[0] - process->nb_state);

  process->observation[0]->dispersion = sqrt(variance);
  for (i = 1;i < process->nb_state;i++) {
    process->observation[i]->dispersion = sqrt(variance * factor[i]);
  }

  delete [] component_frequency;
  delete [] mean;
  delete [] factor;
}


/*--------------------------------------------------------------*
 *
 *  Estimation des parametres des lois de von Mises d'un processus d'observation.
 *
 *  arguments : comptage des etats, variable, pointeur sur un processus d'observation continu.
 *
 *--------------------------------------------------------------*/

template <typename Type>
void Vectors::von_mises_estimation(Type **component_vector_count , int variable ,
                                   ContinuousParametricProcess *process) const

{
  register int i , j;
  double buff , global_mean_direction , concentration , **mean_direction;
  Type *component_frequency;


  component_frequency = new Type[process->nb_state];
  mean_direction = new double*[process->nb_state];
  for (i = 0;i < process->nb_state;i++) {
    mean_direction[i] = new double[4];
  }

  for (i = 0;i < process->nb_state;i++) {
    mean_direction[i][0] = 0.;
    mean_direction[i][1] = 0.;

    component_frequency[i] = 0;
  }

  switch (type[variable]) {

  case INT_VALUE : {
    for (i = 0;i < nb_vector;i++) {
      for (j = 0;j < process->nb_state;j++) {
        mean_direction[j][0] += component_vector_count[i][j] * cos(int_vector[i][variable] * M_PI / 180);
        mean_direction[j][1] += component_vector_count[i][j] * sin(int_vector[i][variable] * M_PI / 180);
        component_frequency[j] += component_vector_count[i][j];
      }
    }
    break;
  }

  case REAL_VALUE : {
    for (i = 0;i < nb_vector;i++) {
      for (j = 0;j < process->nb_state;j++) {
        switch (process->unit) {
        case DEGREE :
          mean_direction[j][0] += component_vector_count[i][j] * cos(real_vector[i][variable] * M_PI / 180);
          mean_direction[j][1] += component_vector_count[i][j] * sin(real_vector[i][variable] * M_PI / 180);
          break;
        case RADIAN :
          mean_direction[j][0] += component_vector_count[i][j] * cos(real_vector[i][variable]);
          mean_direction[j][1] += component_vector_count[i][j] * sin(real_vector[i][variable]);
          break;
        }

        component_frequency[j] += component_vector_count[i][j];
      }
    }
    break;
  }
  }

  for (i = 0;i < process->nb_state;i++) {
    if (component_frequency[i] > 0) {
      mean_direction[i][0] /= component_frequency[i];
      mean_direction[i][1] /= component_frequency[i];

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
      if (component_frequency[i] > 0) {
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
      global_mean_direction += component_frequency[i] * mean_direction[i][2];
      buff += component_frequency[i];
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

  delete [] component_frequency;
}


};  // namespace stat_tool



#endif
