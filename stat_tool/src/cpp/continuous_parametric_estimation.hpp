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

#include <boost/math/distributions/gamma.hpp>

#include "stat_tool/markovian.h"
#include "stat_tool/vectors.h"
#include "stat_tool/stat_label.h"

using namespace std;
using namespace boost::math;


namespace stat_tool {


extern double von_mises_concentration_computation(double mean_direction);



/*--------------------------------------------------------------*/
/**
 *  \brief Estimation of gamma observation distributions.
 *
 *  \param[in] component_vector_count component counts,
 *  \param[in] variable               variable index,
 *  \param[in] process                pointer on a ContinuousParametricProcess object,
 *  \param[in] iter                   EM iteration.
 */
/*--------------------------------------------------------------*/

template <typename Type>
void Vectors::gamma_estimation(Type **component_vector_count , int variable ,
                               ContinuousParametricProcess *process , int iter) const

{
  int i , j;
  double buff , diff , log_geometric_mean , *zero_mass , *bmean;
  long double *variance;
  Type *component_frequency;


  component_frequency = new Type[process->nb_state];
  zero_mass = new double[process->nb_state];
  bmean = new double[process->nb_state];
  variance = new long double[process->nb_state];

  for (i = 0;i < process->nb_state;i++) {
    zero_mass[i] = 0.;
    bmean[i] = 0.;
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
          bmean[j] += component_vector_count[i][j] * int_vector[i][variable];
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
          bmean[j] += component_vector_count[i][j] * real_vector[i][variable];
        }
        component_frequency[j] += component_vector_count[i][j];
      }
    }
    break;
  }
  }

  for (i = 0;i < process->nb_state;i++) {
    if (component_frequency[i] > 0) {
      bmean[i] /= component_frequency[i];
    }
    variance[i] = 0.;
  }

  switch (type[variable]) {

  case INT_VALUE : {
    for (i = 0;i < nb_vector;i++) {
      for (j = 0;j < process->nb_state;j++) {
        diff = int_vector[i][variable] - bmean[j];
        variance[j] += component_vector_count[i][j] * diff * diff;
      }
    }
    break;
  }

  case REAL_VALUE : {
    for (i = 0;i < nb_vector;i++) {
      for (j = 0;j < process->nb_state;j++) {
        diff = real_vector[i][variable] - bmean[j];
        variance[j] += component_vector_count[i][j] * diff * diff;
      }
    }
    break;
  }
  }

  for (i = 0;i < process->nb_state;i++) {
//    if (component_frequency[i] > 0) {
//      variance[i] /= component_frequency[i];
    if (component_frequency[i] > 1) {
      variance[i] /= (component_frequency[i] - 1);
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

//    if (component_frequency[i] > 0) {
    if (component_frequency[i] > 1) {
      if (zero_mass[i] / component_frequency[i] > GAMMA_ZERO_FREQUENCY_THRESHOLD) {
        process->observation[i]->shape = 0;
        process->observation[i]->scale = D_DEFAULT;
      }

      else {
        if (variance[i] > 0.) {
/*          if (sqrtl(variance[i]) < bmean[i] * GAMMA_VARIATION_COEFF_THRESHOLD) {
            variance[i] = bmean[i] * bmean[i] * GAMMA_VARIATION_COEFF_THRESHOLD * GAMMA_VARIATION_COEFF_THRESHOLD;
          }
          process->observation[i]->shape = bmean[i] * bmean[i] / variance[i];
          process->observation[i]->scale = variance[i] / bmean[i]; */

          // Hawang & Huang (2012), Ann. Inst. Statist. Math. 54(4), 840-847

          buff = bmean[i] * bmean[i] / variance[i];
          if (buff > GAMMA_INVERSE_SAMPLE_SIZE_FACTOR / (double)component_frequency[i]) {
            process->observation[i]->shape = buff - 1. / (double)component_frequency[i];
          }
          else {
            process->observation[i]->shape = buff;
          }
/*          if (process->observation[i]->shape < GAMMA_MIN_SHAPE_PARAMETER) {
            process->observation[i]->shape = GAMMA_MIN_SHAPE_PARAMETER;
          } */
          process->observation[i]->scale = bmean[i] / process->observation[i]->shape;

#         ifdef MESSAGE
          if (buff <= GAMMA_INVERSE_SAMPLE_SIZE_FACTOR / (double)component_frequency[i]) {
            cout << "\n" << STAT_word[STATW_COMPONENT] << " " << i << "   "
                 << STAT_word[STATW_SHAPE] << " : " << process->observation[i]->shape << "   "
                 << STAT_word[STATW_SCALE] << " : " << process->observation[i]->scale << " |  ";

            boost::math::gamma_distribution<double> dist(process->observation[i]->shape , process->observation[i]->scale);

            cout << pdf(dist , 0.) << " " << pdf(dist , CONTINUOUS_POSITIVE_INF_BOUND) << " |  "
                 << cdf(dist , 0.) << " " << cdf(dist , CONTINUOUS_POSITIVE_INF_BOUND) << endl;
          }
#         endif

          if ((process->observation[i]->shape >= GAMMA_SHAPE_PARAMETER_THRESHOLD) &&
              (component_frequency[i] < GAMMA_FREQUENCY_THRESHOLD)) {

#           ifdef DEBUG
            cout << STAT_word[STATW_COMPONENT] << " " << i << "   "
                 << STAT_word[STATW_SHAPE] << " : " << process->observation[i]->shape << "   "
                 << STAT_word[STATW_SCALE] << " : " << process->observation[i]->scale;
#           endif

            log_geometric_mean = 0.;

            switch (type[variable]) {

            case INT_VALUE : {
              for (j = 0;j < nb_vector;j++) {
                if (int_vector[j][variable] > 0) {
                  log_geometric_mean += component_vector_count[j][i] * log(int_vector[j][variable]);
                }
              }
              break;
            }

            case REAL_VALUE : {
              for (j = 0;j < nb_vector;j++) {
                if (real_vector[j][variable] > 0.) {
                  log_geometric_mean += component_vector_count[j][i] * log(real_vector[j][variable]);
                }
              }
              break;
            }
            }

            log_geometric_mean /= (component_frequency[i] - zero_mass[i]);
/*            j = 0;   to be reworked

            do {
              process->observation[i]->scale = exp(log_geometric_mean - digamma(process->observation[i]->shape));
              process->observation[i]->shape = bmean[i] / process->observation[i]->scale;
              j++;

#             ifdef DEBUG
              cout << STAT_word[STATW_SHAPE] << " : " << process->observation[i]->shape << "   "
                   << STAT_word[STATW_SCALE] << " : " << process->observation[i]->scale << endl;
#             endif

            }
            while (j < MIN(GAMMA_ITERATION_FACTOR * iter , GAMMA_MAX_NB_ITERATION)); */

            // approximations Johnson, Kotz & Balakrishnan, Continuous Univariate Distributions, vol. 1, 2nd ed., pp. 361-362

//            process->observation[i]->shape = bmean[i] / (2 * (bmean[i] - exp(log_geometric_mean))) - 1./12.;
            diff = log(bmean[i]) - log_geometric_mean;
            process->observation[i]->shape = (1 + sqrt(1 + 4 * diff / 3)) / (4 * diff);
            process->observation[i]->scale = bmean[i] / process->observation[i]->shape;

#           ifdef DEBUG
            cout << " | " << STAT_word[STATW_SHAPE] << " : " << process->observation[i]->shape << "   "
                 << STAT_word[STATW_SCALE] << " : " << process->observation[i]->scale << endl;
#           endif

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
  delete [] bmean;
  delete [] variance;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Estimation of gamma observation distributions:
 *         tied shape parameters and eventually tied scale parameters.
 *
 *  \param[in] component_vector_count component counts,
 *  \param[in] variable               variable index,
 *  \param[in] process                pointer on a ContinuousParametricProcess object,
 *  \param[in] variance_factor        tying rule,
 *  \param[in] iter                   EM iteration.
 */
/*--------------------------------------------------------------*/

template <typename Type>
void Vectors::tied_gamma_estimation(Type **component_vector_count , int variable ,
                                    ContinuousParametricProcess *process ,
                                    tying_rule variance_factor , int iter) const

{
  int i , j;
  double buff , diff , log_geometric_mean , *bmean , *factor;
  long double variance;
  Type *component_frequency;


  component_frequency = new Type[process->nb_state];
  bmean = new double[process->nb_state];
  factor = new double[process->nb_state];

  if (process->tied_location) {
    bmean[0] = 0.;
    component_frequency[0] = 0.;

    factor[0] = 1.;
    for (i = 1;i < process->nb_state;i++) {
      factor[i] = factor[i - 1] * 2;
    }

    switch (type[variable]) {

    case INT_VALUE : {
      for (i = 0;i < nb_vector;i++) {
        for (j = 0;j < process->nb_state;j++) {
          bmean[0] += component_vector_count[i][j] * int_vector[i][variable] / factor[j];
          component_frequency[0] += component_vector_count[i][j];
        }
      }
      break;
    }

    case REAL_VALUE : {
      for (i = 0;i < nb_vector;i++) {
        for (j = 0;j < process->nb_state;j++) {
          bmean[0] += component_vector_count[i][j] * real_vector[i][variable] / factor[j];
          component_frequency[0] += component_vector_count[i][j];
        }
      }
      break;
    }
    }

    bmean[0] /= component_frequency[0];
    for (i = 1;i < process->nb_state;i++) {
      bmean[i] = bmean[0] * factor[i];
    }

    if (variance_factor == SCALING_FACTOR) {
      for (i = 1;i < process->nb_state;i++) {
        factor[i] *= factor[i];
//        factor[i] = factor[i - 1] * 4;
      }
    }
  }

  else {
    for (i = 0;i < process->nb_state;i++) {
      bmean[i] = 0.;
      component_frequency[i] = 0;
    }

    switch (type[variable]) {

    case INT_VALUE : {
      for (i = 0;i < nb_vector;i++) {
        for (j = 0;j < process->nb_state;j++) {
          bmean[j] += component_vector_count[i][j] * int_vector[i][variable];
          component_frequency[j] += component_vector_count[i][j];
        }
      }
      break;
    }

    case REAL_VALUE : {
      for (i = 0;i < nb_vector;i++) {
        for (j = 0;j < process->nb_state;j++) {
          bmean[j] += component_vector_count[i][j] * real_vector[i][variable];
          component_frequency[j] += component_vector_count[i][j];
        }
      }
      break;
    }
    }

    for (i = 0;i < process->nb_state;i++) {
      if (component_frequency[i] > 0) {
        bmean[i] /= component_frequency[i];
      }
    }

    switch (variance_factor) {

    case CONVOLUTION_FACTOR : {
      factor[0] = 1.;
      for (i = 1;i < process->nb_state;i++) {
        if (component_frequency[i] > 0) {
          factor[i] = bmean[i] / bmean[0];
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
          factor[i] = (bmean[i] * bmean[i]) / (bmean[0] * bmean[0]);
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
  }

  variance = 0.;

  switch (type[variable]) {

  case INT_VALUE : {
    for (i = 0;i < nb_vector;i++) {
      for (j = 0;j < process->nb_state;j++) {
        diff = int_vector[i][variable] - bmean[j];
        variance += component_vector_count[i][j] * diff * diff / factor[j];
      }
    }
    break;
  }

  case REAL_VALUE : {
    for (i = 0;i < nb_vector;i++) {
      for (j = 0;j < process->nb_state;j++) {
        diff = real_vector[i][variable] - bmean[j];
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
    process->observation[0]->shape = bmean[0] * bmean[0] / variance;
    process->observation[0]->scale = variance / bmean[0];

    for (i = 1;i < process->nb_state;i++) {
      process->observation[i]->shape = process->observation[0]->shape * factor[i];
      process->observation[i]->scale = process->observation[0]->scale;
    }
    break;
  }

  case SCALING_FACTOR : {
    if (process->tied_location) {
      process->observation[0]->shape = bmean[0] * bmean[0] / variance;
      process->observation[0]->scale = variance / bmean[0];

      for (i = 1;i < process->nb_state;i++) {
        process->observation[i]->shape = process->observation[0]->shape;
//        factor[i] = sqrt(factor[i]);
        factor[i] = factor[i - 1] * 2;
        process->observation[i]->scale = process->observation[0]->scale * factor[i];
      }
    }

    else {
//      process->observation[0]->shape = bmean[0] * bmean[0] / variance;
//      process->observation[0]->scale = variance / bmean[0];

      // Hawang & Huang (2012), Ann. Inst. Statist. Math. 54(4), 840-847

      buff = bmean[0] * bmean[0] / variance;
      if (buff > GAMMA_INVERSE_SAMPLE_SIZE_FACTOR / (double)component_frequency[0]) {
        process->observation[0]->shape = buff - 1. / (double)component_frequency[0];
      }
      else {
        process->observation[0]->shape = buff;
      }
/*      if (process->observation[0]->shape < GAMMA_MIN_SHAPE_PARAMETER) {
        process->observation[0]->shape = GAMMA_MIN_SHAPE_PARAMETER;
      } */
      process->observation[0]->scale = bmean[0] / process->observation[0]->shape;

      for (i = 1;i < process->nb_state;i++) {
        process->observation[i]->shape = process->observation[0]->shape;
//        factor[i] = sqrt(factor[i]);
        factor[i] = bmean[i] / bmean[0];
        process->observation[i]->scale = process->observation[0]->scale * factor[i];
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

//    process->observation[0]->shape = bmean[0] / (2 * (bmean[0] - exp(log_geometric_mean))) - 1./12.;
    diff = log(bmean[0]) - log_geometric_mean;
    process->observation[0]->shape = (1 + sqrt(1 + 4 * diff / 3)) / (4 * diff);

#   ifdef DEBUG
    cout << STAT_word[STATW_SHAPE] << " : " << process->observation[0]->shape << endl;
#   endif

    i = 0;

/*    do {
      process->observation[0]->scale = exp(log_geometric_mean - digamma(process->observation[0]->shape));
      process->observation[0]->shape = bmean[0] / process->observation[0]->scale;
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

    process->observation[0]->scale = bmean[0] / process->observation[0]->shape;

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
  delete [] bmean;
  delete [] factor;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Estimation of inverse Gaussian observation distributions.
 *
 *  \param[in] component_vector_count component counts,
 *  \param[in] variable               variable index,
 *  \param[in] process                pointer on a ContinuousParametricProcess object.
 */
/*--------------------------------------------------------------*/

template <typename Type>
void Vectors::inverse_gaussian_estimation(Type **component_vector_count , int variable ,
                                          ContinuousParametricProcess *process) const

{
  int i , j;
  double *bmean , *inverse_scale;
  Type *component_frequency;


  component_frequency = new Type[process->nb_state];
  bmean = new double[process->nb_state];

  for (i = 0;i < process->nb_state;i++) {
    bmean[i] = 0.;
    component_frequency[i] = 0;
  }

  switch (type[variable]) {

  case INT_VALUE : {
    for (i = 0;i < nb_vector;i++) {
      for (j = 0;j < process->nb_state;j++) {
        bmean[j] += component_vector_count[i][j] * int_vector[i][variable];
        component_frequency[j] += component_vector_count[i][j];
      }
    }
    break;
  }

  case REAL_VALUE : {
    for (i = 0;i < nb_vector;i++) {
      for (j = 0;j < process->nb_state;j++) {
        bmean[j] += component_vector_count[i][j] * real_vector[i][variable];
        component_frequency[j] += component_vector_count[i][j];
      }
    }
    break;
  }
  }

  for (i = 0;i < process->nb_state;i++) {
    if (component_frequency[i] > 0) {
      bmean[i] /= component_frequency[i];
      process->observation[i]->location = bmean[i];
    }
    else {
      process->observation[i]->location = D_DEFAULT;
    }
  }

  inverse_scale = new double[process->nb_state];
  for (i = 0;i < process->nb_state;i++) {
    inverse_scale[i] = 0.;
  }

  switch (type[variable]) {

  case INT_VALUE : {
    for (i = 0;i < nb_vector;i++) {
      for (j = 0;j < process->nb_state;j++) {
        if ((bmean[j] > 0.) && (int_vector[i][variable] > 0)) {
          inverse_scale[j] += component_vector_count[i][j] * (1. / (double)int_vector[i][variable] - 1. / bmean[j]);
        }
      }
    }
    break;
  }

  case REAL_VALUE : {
    for (i = 0;i < nb_vector;i++) {
      for (j = 0;j < process->nb_state;j++) {
        if ((bmean[j] > 0.) && (real_vector[i][variable] > 0.)) {
          inverse_scale[j] += component_vector_count[i][j] * (1. / real_vector[i][variable] - 1. / bmean[j]);
        }
      }
    }
    break;
  }
  }

  for (i = 0;i < process->nb_state;i++) {
    if (inverse_scale[i] > 0.) {
      process->observation[i]->scale = component_frequency[i] / inverse_scale[i];
    }
    else {
      process->observation[i]->scale = D_DEFAULT;
    }
  }

  delete [] component_frequency;
  delete [] bmean;
  delete [] inverse_scale;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Estimation of parameters of inverse Gaussian observation distributions:
 *         tied variances and eventually tied means.
 *
 *  \param[in] component_vector_count component counts,
 *  \param[in] variable               variable index,
 *  \param[in] process                pointer on a ContinuousParametricProcess object,
 *  \param[in] variance_factor        tying rule for the variance.
 */
/*--------------------------------------------------------------*/

template <typename Type>
void Vectors::tied_inverse_gaussian_estimation(Type **component_vector_count , int variable ,
                                               ContinuousParametricProcess *process ,
                                               tying_rule variance_factor) const

{
  int i , j;
  double inverse_scale , *bmean , *factor;
  Type *component_frequency;

  double inverse_scale2;


  component_frequency = new Type[process->nb_state];
  bmean = new double[process->nb_state];
  factor = new double[process->nb_state];

  if (process->tied_location) {
    bmean[0] = 0.;
    component_frequency[0] = 0.;

    factor[0] = 1.;
    for (i = 1;i < process->nb_state;i++) {
      factor[i] = factor[i - 1] * 2;
    }

    switch (type[variable]) {

    case INT_VALUE : {
      for (i = 0;i < nb_vector;i++) {
        for (j = 0;j < process->nb_state;j++) {
          bmean[0] += component_vector_count[i][j] * int_vector[i][variable] / factor[j];
          component_frequency[0] += component_vector_count[i][j];
        }
      }
      break;
    }

    case REAL_VALUE : {
      for (i = 0;i < nb_vector;i++) {
        for (j = 0;j < process->nb_state;j++) {
          bmean[0] += component_vector_count[i][j] * real_vector[i][variable] / factor[j];
          component_frequency[0] += component_vector_count[i][j];
        }
      }
      break;
    }
    }

    bmean[0] /= component_frequency[0];
    process->observation[0]->location = bmean[0];
    for (i = 1;i < process->nb_state;i++) {
      bmean[i] = bmean[0] * factor[i];
      process->observation[i]->location = bmean[i];
    }

    if (variance_factor == CONVOLUTION_FACTOR) {
      for (i = 1;i < process->nb_state;i++) {
//        factor[i] = factor[i - 1] * 4;
        factor[i] *= factor[i];
      }
    }

    inverse_scale = 0.;

    switch (variance_factor) {

    case CONVOLUTION_FACTOR : {
      switch (type[variable]) {

      case INT_VALUE : {
        for (i = 0;i < nb_vector;i++) {
          for (j = 0;j < process->nb_state;j++) {
            if ((bmean[j] > 0.) && (int_vector[i][variable] > 0)) {
              inverse_scale += component_vector_count[i][j] * factor[j] * (1. / (double)int_vector[i][variable] - 2. / bmean[j] +
                                                                           int_vector[i][variable] / (bmean[j] * bmean[j]));
            }
          }
        }
        break;
      }

      case REAL_VALUE : {
        for (i = 0;i < nb_vector;i++) {
          for (j = 0;j < process->nb_state;j++) {
            if ((bmean[j] > 0.) && (real_vector[i][variable] > 0.)) {
              inverse_scale += component_vector_count[i][j] * factor[j] * (1. / real_vector[i][variable] - 2. / bmean[j] +
                                                                           real_vector[i][variable] / (bmean[j] * bmean[j]));
            }
          }
        }
        break;
      }
      }
      break;
    }

    case SCALING_FACTOR : {
      inverse_scale2 = 0.;

      switch (type[variable]) {

      case INT_VALUE : {
        for (i = 0;i < nb_vector;i++) {
          for (j = 0;j < process->nb_state;j++) {
            if ((bmean[j] > 0.) && (int_vector[i][variable] > 0)) {
              inverse_scale += component_vector_count[i][j] * factor[j] * (1. / (double)int_vector[i][variable] - 1. / bmean[j]);
              inverse_scale2 += component_vector_count[i][j] * factor[j] * (1. / (double)int_vector[i][variable] - 2. / bmean[j] +
                                                                            int_vector[i][variable] / (bmean[j] * bmean[j]));
            }
          }
        }
        break;
      }

      case REAL_VALUE : {
        for (i = 0;i < nb_vector;i++) {
          for (j = 0;j < process->nb_state;j++) {
            if ((bmean[j] > 0.) && (real_vector[i][variable] > 0.)) {
              inverse_scale += component_vector_count[i][j] * factor[j] * (1. / real_vector[i][variable] - 1. / bmean[j]);
              inverse_scale2 += component_vector_count[i][j] * factor[j] * (1. / real_vector[i][variable] - 2. / bmean[j] +
                                                                           real_vector[i][variable] / (bmean[j] * bmean[j]));
            }
          }
        }
        break;
      }
      }

#     ifdef MESSAGE
      if ((inverse_scale2 < inverse_scale - DOUBLE_ERROR) || (inverse_scale2 > inverse_scale + DOUBLE_ERROR)) {
        cout << "\nERROR: " << inverse_scale << " | " << inverse_scale2 << endl;
      }
#     endif

      break;
    }
    }
  }

  else {
    for (i = 0;i < process->nb_state;i++) {
      bmean[i] = 0.;
      component_frequency[i] = 0;
    }

    switch (type[variable]) {

    case INT_VALUE : {
      for (i = 0;i < nb_vector;i++) {
        for (j = 0;j < process->nb_state;j++) {
          bmean[j] += component_vector_count[i][j] * int_vector[i][variable];
          component_frequency[j] += component_vector_count[i][j];
        }
      }
      break;
    }

    case REAL_VALUE : {
      for (i = 0;i < nb_vector;i++) {
        for (j = 0;j < process->nb_state;j++) {
          bmean[j] += component_vector_count[i][j] * real_vector[i][variable];
          component_frequency[j] += component_vector_count[i][j];
        }
      }
      break;
    }
    }

    for (i = 0;i < process->nb_state;i++) {
      if (component_frequency[i] > 0) {
        bmean[i] /= component_frequency[i];
        process->observation[i]->location = bmean[i];
      }
      else {
        process->observation[i]->location = D_DEFAULT;
      }
    }

    switch (variance_factor) {

    case CONVOLUTION_FACTOR : {
      factor[0] = 1.;
      for (i = 1;i < process->nb_state;i++) {
        if (component_frequency[i] > 0) {
          factor[i] = (bmean[i] * bmean[i]) / (bmean[0] * bmean[0]);
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
          factor[i] = bmean[i] / bmean[0];
        }
        else {
          factor[i] = 1.;
        }
      }
      break;
    }
    }

    inverse_scale = 0.;
    inverse_scale2 = 0.;

    switch (type[variable]) {

    case INT_VALUE : {
      for (i = 0;i < nb_vector;i++) {
        for (j = 0;j < process->nb_state;j++) {
          if ((bmean[j] > 0.) && (int_vector[i][variable] > 0)) {
            inverse_scale += component_vector_count[i][j] * factor[j] * (1. / (double)int_vector[i][variable] - 1. / bmean[j]);
            inverse_scale2 += component_vector_count[i][j] * factor[j] * (1. / (double)int_vector[i][variable] - 2. / bmean[j] +
                                                                          int_vector[i][variable] / (bmean[j] * bmean[j]));
          }
        }
      }
      break;
    }

    case REAL_VALUE : {
      for (i = 0;i < nb_vector;i++) {
        for (j = 0;j < process->nb_state;j++) {
          if ((bmean[j] > 0.) && (real_vector[i][variable] > 0.)) {
            inverse_scale += component_vector_count[i][j] * factor[j] * (1. / real_vector[i][variable] - 1 / bmean[j]);
            inverse_scale2 += component_vector_count[i][j] * factor[j] * (1. / real_vector[i][variable] - 2. / bmean[j] +
                                                                          real_vector[i][variable] / (bmean[j] * bmean[j]));
          }
        }
      }
      break;
    }
    }

#   ifdef MESSAGE
/*    cout << "\n";
    for (i = 1;i < process->nb_state;i++) {
      cout << factor[i];
      if (variance_factor == CONVOLUTION_FACTOR) {
        cout << " (" << sqrt(factor[i]) << ")";
      }
      if (i < process->nb_state - 1) {
        cout << ", ";
      }
    }
    cout << endl; */

    if ((inverse_scale2 < inverse_scale - DOUBLE_ERROR) || (inverse_scale2 > inverse_scale + DOUBLE_ERROR)) {
      cout << "\nERROR: " << inverse_scale << " | " << inverse_scale2 << endl;
    }
#   endif

    for (i = 1;i < process->nb_state;i++) {
      component_frequency[0] += component_frequency[i];
    }
  }

  process->observation[0]->scale = component_frequency[0] / inverse_scale;
  for (i = 1;i < process->nb_state;i++) {
    process->observation[i]->scale = process->observation[0]->scale * factor[i];
  }

  delete [] component_frequency;
  delete [] bmean;
  delete [] factor;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Estimation of Gaussian observation distributions.
 *
 *  \param[in] component_vector_count component counts,
 *  \param[in] variable               variable index,
 *  \param[in] process                pointer on a ContinuousParametricProcess object.
 */
/*--------------------------------------------------------------*/

template <typename Type>
void Vectors::gaussian_estimation(Type **component_vector_count , int variable ,
                                  ContinuousParametricProcess *process) const

{
  int i , j;
  double diff , *bmean;
  long double *variance;
  Type *component_frequency;


  component_frequency = new Type[process->nb_state];
  bmean = new double[process->nb_state];

  for (i = 0;i < process->nb_state;i++) {
    bmean[i] = 0.;
    component_frequency[i] = 0;
  }

  switch (type[variable]) {

  case INT_VALUE : {
    for (i = 0;i < nb_vector;i++) {
      for (j = 0;j < process->nb_state;j++) {
        bmean[j] += component_vector_count[i][j] * int_vector[i][variable];
        component_frequency[j] += component_vector_count[i][j];
      }
    }
    break;
  }

  case REAL_VALUE : {
    for (i = 0;i < nb_vector;i++) {
      for (j = 0;j < process->nb_state;j++) {
        bmean[j] += component_vector_count[i][j] * real_vector[i][variable];
        component_frequency[j] += component_vector_count[i][j];
      }
    }
    break;
  }
  }

  for (i = 0;i < process->nb_state;i++) {
    if (component_frequency[i] > 0) {
      bmean[i] /= component_frequency[i];
      process->observation[i]->location = bmean[i];
    }
    else {
      process->observation[i]->location = D_INF;
    }
  }

  if (process->tied_dispersion) {
    for (i = 1;i < process->nb_state;i++) {
      component_frequency[0] += component_frequency[i];
    }

    variance = new long double[1];
    variance[0] = 0.;

    switch (type[variable]) {

    case INT_VALUE : {
      for (i = 0;i < nb_vector;i++) {
        for (j = 0;j < process->nb_state;j++) {
          diff = int_vector[i][variable] - bmean[j];
          variance[0] += component_vector_count[i][j] * diff * diff;
        }
      }
      break;
    }

    case REAL_VALUE : {
      for (i = 0;i < nb_vector;i++) {
        for (j = 0;j < process->nb_state;j++) {
          diff = real_vector[i][variable] - bmean[j];
          variance[0] += component_vector_count[i][j] * diff * diff;
        }
      }
      break;
    }
    }

//    variance[0] /= component_frequency[0];
    variance[0] /= (component_frequency[0] - process->nb_state);

    process->observation[0]->dispersion = sqrtl(variance[0]);
    for (i = 1;i < process->nb_state;i++) {
      process->observation[i]->dispersion = process->observation[0]->dispersion;
    }
  }

  else {
    variance = new long double[process->nb_state];
    for (i = 0;i < process->nb_state;i++) {
      variance[i] = 0.;
    }

    switch (type[variable]) {

    case INT_VALUE : {
      for (i = 0;i < nb_vector;i++) {
        for (j = 0;j < process->nb_state;j++) {
          diff = int_vector[i][variable] - bmean[j];
          variance[j] += component_vector_count[i][j] * diff * diff;
        }
      }
      break;
    }

    case REAL_VALUE : {
      for (i = 0;i < nb_vector;i++) {
        for (j = 0;j < process->nb_state;j++) {
          diff = real_vector[i][variable] - bmean[j];
          variance[j] += component_vector_count[i][j] * diff * diff;
        }
      }
      break;
    }
    }

    for (i = 0;i < process->nb_state;i++) {
//      if (component_frequency[i] > 0) {
//        variance[i] /= component_frequency[i];
      if (component_frequency[i] > 1) {
        variance[i] /= (component_frequency[i] - 1);
        process->observation[i]->dispersion = sqrtl(variance[i]);
        if ((process->observation[i]->location != 0.) &&
            (process->observation[i]->dispersion / process->observation[i]->location < GAUSSIAN_MIN_VARIATION_COEFF)) {
          process->observation[i]->dispersion = process->observation[i]->location * GAUSSIAN_MIN_VARIATION_COEFF;
        }
      }
    }
  }

  delete [] component_frequency;
  delete [] bmean;
  delete [] variance;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Estimation of Gaussian observation distributions with evenly spaced means. 
 *
 *  \param[in] component_vector_count component counts,
 *  \param[in] variable               variable index,
 *  \param[in] process                pointer on a ContinuousParametricProcess object.
 */
/*--------------------------------------------------------------*/

template <typename Type>
void Vectors::tied_gaussian_estimation(Type **component_vector_count , int variable ,
                                       ContinuousParametricProcess *process) const

{
  int i , j;
  double old_offset , diff , *bmean;
  long double *variance;
  Type *component_frequency;


  bmean = new double[2];
  component_frequency = new Type[process->nb_state];

  old_offset = process->offset;
  process->offset = 0.;
  bmean[0] = 0.;
  component_frequency[0] = 0;
  bmean[1] = 0.;
  component_frequency[1] = 0;

  switch (type[variable]) {

  case INT_VALUE : {
    for (i = 0;i < nb_vector;i++) {
      for (j = 0;j < process->nb_state - 1;j++) {
        process->offset += component_vector_count[i][j] * (int_vector[i][variable] - (process->observation[j]->location - old_offset));
        bmean[0] += component_vector_count[i][j] * (int_vector[i][variable] - old_offset) / (j + 1);
        component_frequency[0] += component_vector_count[i][j];
      }

      bmean[1] += component_vector_count[i][process->nb_state - 1] * int_vector[i][variable];
      component_frequency[1] += component_vector_count[i][process->nb_state - 1];
    }
    break;
  }

  case REAL_VALUE : {
    for (i = 0;i < nb_vector;i++) {
      for (j = 0;j < process->nb_state - 1;j++) {
        process->offset += component_vector_count[i][j] * (real_vector[i][variable] - (process->observation[j]->location - old_offset));
        bmean[0] += component_vector_count[i][j] * (real_vector[i][variable] - old_offset) / (j + 1);
        component_frequency[0] += component_vector_count[i][j];
      }

      bmean[1] += component_vector_count[i][process->nb_state - 1] * real_vector[i][variable];
      component_frequency[1] += component_vector_count[i][process->nb_state - 1];
    }
    break;
  }
  }

  process->offset /= component_frequency[0];

  bmean[0] /= component_frequency[0];
  for (i = 0;i < process->nb_state - 1;i++) {
    process->observation[i]->location = process->offset + (i + 1) * bmean[0];
  }
  process->observation[process->nb_state - 1]->location = bmean[1] / component_frequency[1];

  if (process->tied_dispersion) {
    variance = new long double[2];
    variance[0] = 0.;
    variance[1] = 0.;

    switch (type[variable]) {

    case INT_VALUE : {
      for (i = 0;i < nb_vector;i++) {
        for (j = 0;j < process->nb_state - 1;j++) {
          diff = int_vector[i][variable] - process->observation[j]->location;
          variance[0] += component_vector_count[i][j] * diff * diff;
        }

        diff = int_vector[i][variable] - process->observation[process->nb_state - 1]->location;
        variance[1] += component_vector_count[i][process->nb_state - 1] * diff * diff;
      }
      break;
    }

    case REAL_VALUE : {
      for (i = 0;i < nb_vector;i++) {
        for (j = 0;j < process->nb_state - 1;j++) {
          diff = real_vector[i][variable] - process->observation[j]->location;
          variance[0] += component_vector_count[i][j] * diff * diff;
        }

        diff = real_vector[i][variable] - process->observation[process->nb_state - 1]->location;
        variance[1] += component_vector_count[i][process->nb_state - 1] * diff * diff;
      }
      break;
    }
    }

//    variance[0] /= component_frequency[0];
    variance[0] /= (component_frequency[0] - (process->nb_state - 1));

    process->observation[0]->dispersion = sqrtl(variance[0]);
    for (i = 1;i < process->nb_state - 1;i++) {
      process->observation[i]->dispersion = process->observation[0]->dispersion;
    }

//    variance[1] /= component_frequency[1];
    variance[1] /= (component_frequency[1] - 1);
    if (sqrtl(variance[1]) > 2 * process->observation[0]->dispersion) {
      process->observation[process->nb_state - 1]->dispersion = sqrtl(variance[1]);
    }    

#   ifdef MESSAGE
    cout << "TEST: " << process->offset << " " << bmean[0] << " " << process->observation[0]->dispersion << endl;
#   endif

  }

  else {
    variance = new long double[process->nb_state];
    for (i = 0;i < process->nb_state;i++) {
      variance[i] = 0.;
      component_frequency[i] = 0;
    }

    switch (type[variable]) {

    case INT_VALUE : {
      for (i = 0;i < nb_vector;i++) {
        for (j = 0;j < process->nb_state;j++) {
          diff = int_vector[i][variable] - process->observation[j]->location;
          variance[j] += component_vector_count[i][j] * diff * diff;
          component_frequency[j] += component_vector_count[i][j];
        }
      }
      break;
    }

    case REAL_VALUE : {
      for (i = 0;i < nb_vector;i++) {
        for (j = 0;j < process->nb_state;j++) {
          diff = real_vector[i][variable] - process->observation[j]->location;
          variance[j] += component_vector_count[i][j] * diff * diff;
          component_frequency[j] += component_vector_count[i][j];
        }
      }
      break;
    }
    }

    for (i = 0;i < process->nb_state;i++) {
//      if (component_frequency[i] > 0) {
//        variance[i] /= component_frequency[i];
      if (component_frequency[i] > 1) {
        variance[i] /= (component_frequency[i] - 1);
        process->observation[i]->dispersion = sqrtl(variance[i]);
        if ((process->observation[i]->location != 0.) &&
            (process->observation[i]->dispersion / process->observation[i]->location < GAUSSIAN_MIN_VARIATION_COEFF)) {
          process->observation[i]->dispersion = process->observation[i]->location * GAUSSIAN_MIN_VARIATION_COEFF;
        }
      }
    }
  }

  delete [] component_frequency;
  delete [] bmean;
  delete [] variance;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Estimation of parameters of Gaussian observation distributions:
 *         tied variances and eventually tied means.
 *
 *  \param[in] component_vector_count component counts,
 *  \param[in] variable               variable index,
 *  \param[in] process                pointer on a ContinuousParametricProcess object,
 *  \param[in] variance_factor        tying rule for the variance.
 */
/*--------------------------------------------------------------*/

template <typename Type>
void Vectors::tied_gaussian_estimation(Type **component_vector_count , int variable ,
                                       ContinuousParametricProcess *process ,
                                       tying_rule variance_factor) const

{
  int i , j;
  double diff , *bmean , *factor;
  long double variance;
  Type *component_frequency;


  component_frequency = new Type[process->nb_state];
  bmean = new double[process->nb_state];
  factor = new double[process->nb_state];

  if (process->tied_location) {
    bmean[0] = 0.;
    component_frequency[0] = 0.;

    factor[0] = 1.;
    for (i = 1;i < process->nb_state;i++) {
      factor[i] = factor[i - 1] * 2;
    }

    switch (type[variable]) {

    case INT_VALUE : {
      for (i = 0;i < nb_vector;i++) {
        for (j = 0;j < process->nb_state;j++) {
          bmean[0] += component_vector_count[i][j] * int_vector[i][variable] / factor[j];
          component_frequency[0] += component_vector_count[i][j];
        }
      }
      break;
    }

    case REAL_VALUE : {
      for (i = 0;i < nb_vector;i++) {
        for (j = 0;j < process->nb_state;j++) {
          bmean[0] += component_vector_count[i][j] * real_vector[i][variable] / factor[j];
          component_frequency[0] += component_vector_count[i][j];
        }
      }
      break;
    }
    }

    bmean[0] /= component_frequency[0];
    process->observation[0]->location = bmean[0];
    for (i = 1;i < process->nb_state;i++) {
      bmean[i] = bmean[0] * factor[i];
      process->observation[i]->location = bmean[i];
    }

    if (variance_factor == SCALING_FACTOR) {
      for (i = 1;i < process->nb_state;i++) {
//        factor[i] = factor[i - 1] * 4;
        factor[i] *= factor[i];
      }
    }
  }

  else {
    for (i = 0;i < process->nb_state;i++) {
      bmean[i] = 0.;
      component_frequency[i] = 0;
    }

    switch (type[variable]) {

    case INT_VALUE : {
      for (i = 0;i < nb_vector;i++) {
        for (j = 0;j < process->nb_state;j++) {
          bmean[j] += component_vector_count[i][j] * int_vector[i][variable];
          component_frequency[j] += component_vector_count[i][j];
        }
      }
      break;
    }

    case REAL_VALUE : {
      for (i = 0;i < nb_vector;i++) {
        for (j = 0;j < process->nb_state;j++) {
          bmean[j] += component_vector_count[i][j] * real_vector[i][variable];
          component_frequency[j] += component_vector_count[i][j];
        }
      }
      break;
    }
    }

    for (i = 0;i < process->nb_state;i++) {
      if (component_frequency[i] > 0) {
        bmean[i] /= component_frequency[i];
        process->observation[i]->location = bmean[i];
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
          factor[i] = bmean[i] / bmean[0];
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
          factor[i] = (bmean[i] * bmean[i]) / (bmean[0] * bmean[0]);
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
  }

  variance = 0.;

  switch (type[variable]) {

  case INT_VALUE : {
    for (i = 0;i < nb_vector;i++) {
      for (j = 0;j < process->nb_state;j++) {
        diff = int_vector[i][variable] - bmean[j];
        variance += component_vector_count[i][j] * diff * diff / factor[j];
      }
    }
    break;
  }

  case REAL_VALUE : {
    for (i = 0;i < nb_vector;i++) {
      for (j = 0;j < process->nb_state;j++) {
        diff = real_vector[i][variable] - bmean[j];
        variance += component_vector_count[i][j] * diff * diff / factor[j];
      }
    }
    break;
  }
  }

//  variance /= component_frequency[0];
  variance /= (component_frequency[0] - process->nb_state);

  process->observation[0]->dispersion = sqrtl(variance);
  for (i = 1;i < process->nb_state;i++) {
    process->observation[i]->dispersion = sqrtl(variance * factor[i]);
  }

  delete [] component_frequency;
  delete [] bmean;
  delete [] factor;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Estimation of von Mises observation distributions.
 *
 *  \param[in] component_vector_count component counts,
 *  \param[in] variable               variable index,
 *  \param[in] process                pointer on a ContinuousParametricProcess object.
 */
/*--------------------------------------------------------------*/

template <typename Type>
void Vectors::von_mises_estimation(Type **component_vector_count , int variable ,
                                   ContinuousParametricProcess *process) const

{
  int i , j;
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

        process->observation[i]->location = mean_direction[i][3];
      }

      else {
        process->observation[i]->location = D_INF;
      }
    }

    else {
      process->observation[i]->location = D_INF;
    }
  }

  if (process->tied_dispersion) {
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
  }

  else {
    for (i = 0;i < process->nb_state;i++) {
      if (component_frequency[i] > 0) {
        process->observation[i]->dispersion = von_mises_concentration_computation(mean_direction[i][2]);
      }
      else {
        process->observation[i]->dispersion = D_DEFAULT;
      }
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
