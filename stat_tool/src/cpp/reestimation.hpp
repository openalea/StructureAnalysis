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



#ifndef REESTIMATION_HPP
#define REESTIMATION_HPP



#include <math.h>

#include <boost/math/special_functions/digamma.hpp>

#include "stat_label.h"

using namespace std;
using namespace boost::math;


namespace stat_tool {



/*--------------------------------------------------------------*/
/**
 *  \brief Construction of a Reestimation object.
 *
 *  \param[in] inb_value number of values from 0.
 */
/*--------------------------------------------------------------*/

template <typename Type>
void Reestimation<Type>::init(int inb_value)

{
  nb_element = 0;
  nb_value = inb_value;
  alloc_nb_value = nb_value;
  offset = 0;
  max = 0;
  mean = D_DEFAULT;
  variance = D_DEFAULT;

  if (nb_value == 0) {
    frequency = NULL;
  }

  else {
    int i;

    frequency = new Type[nb_value];

    for (i = 0;i < nb_value;i++) {
      frequency[i] = 0;
    }
  }
}


/*--------------------------------------------------------------*/
/**
 * \brief Copy of a Reestimation object.
 *
 * \param[in] histo reference on a Reestimation object.
 */
/*--------------------------------------------------------------*/

template <typename Type>
void Reestimation<Type>::copy(const Reestimation<Type> &histo)

{
  int i;


  nb_value = histo.nb_value;
  alloc_nb_value = nb_value;
  nb_element = histo.nb_element;
  offset = histo.offset;
  max = histo.max;
  mean = histo.mean;
  variance = histo.variance;

  frequency = new Type[nb_value];

  for (i = 0;i < nb_value;i++) {
    frequency[i] = histo.frequency[i];
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Constructor by copy of the Reestimation class.
 *
 *  \param[in] histo reference on a Reestimation object.
 */
/*--------------------------------------------------------------*/

template <typename Type>
Reestimation<Type>::Reestimation(const Reestimation<Type> &histo)

{
  copy(histo);
}


/*--------------------------------------------------------------*/
/**
 *  \brief Constructor by merging of the Reestimation class.
 *
 *  \param[in] nb_histo number of  Reestimation objects,
 *  \param[in] histo    pointer on the Reestimation objects.
 */
/*--------------------------------------------------------------*/

template <typename Type>
Reestimation<Type>::Reestimation(int nb_histo , const Reestimation<Type> **histo)

{
  int i , j;


  nb_value = 0;
  nb_element = 0;
  for (i = 0;i < nb_histo;i++) {
    if (histo[i]->nb_value > nb_value) {
      nb_value = histo[i]->nb_value;
    }
    nb_element += histo[i]->nb_element;
  }
  alloc_nb_value = nb_value;

  frequency = new Type[nb_value];

  for (i = 0;i < nb_value;i++) {
    frequency[i] = 0;
  }
  for (i = 0;i < nb_histo;i++) {
    for (j = histo[i]->offset;j < histo[i]->nb_value;j++) {
      frequency[j] += histo[i]->frequency[j];
    }
  }

  // computation of the characteristics of the frequency distributions

  offset_computation();
  max_computation();
  mean_computation();
  variance_computation();
}


/*--------------------------------------------------------------*/
/**
 *  \brief Destructor of the Reestimation class.
 */
/*--------------------------------------------------------------*/

template <typename Type>
Reestimation<Type>::~Reestimation()

{
  delete [] frequency;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Assignment operator of the Reestimation class.
 *
 *  \param[in] histo reference on a Reestimation object.
 *
 *  \return          Reestimation object.
 */
/*--------------------------------------------------------------*/

template <typename Type>
Reestimation<Type>& Reestimation<Type>::operator=(const Reestimation<Type> &histo)

{
  if (&histo != this) {
    delete [] frequency;

    copy(histo);
  }

  return *this;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of the characteristics of a Reestimation object.
 *
 *  \param[in,out] os           stream,
 *  \param[in]     shape        flag on the writing of the shape characteristics,
 *  \param[in]     comment_flag flag comments.
 */
/*--------------------------------------------------------------*/

template <typename Type>
ostream& Reestimation<Type>::ascii_characteristic_print(ostream &os , bool shape ,
                                                        bool comment_flag) const

{
  os << STAT_label[STATL_SAMPLE_SIZE] << ": " << nb_element << endl;

  if ((mean != D_DEFAULT) && (variance != D_DEFAULT)) {
    if (comment_flag) {
      os << "# ";
    }
    os << STAT_label[STATL_MEAN] << ": " << mean << "   "
       << STAT_label[STATL_MEDIAN] << ": " << quantile_computation() << "   "
       << STAT_label[STATL_MODE] << ": " << mode_computation() << endl;

    if (comment_flag) {
      os << "# ";
    }
    os << STAT_label[STATL_VARIANCE] << ": " << variance << "   "
       << STAT_label[STATL_STANDARD_DEVIATION] << ": " << sqrt(variance);
    if (variance > 0.) {
      os << "   " << STAT_label[STATL_LOWER_QUARTILE] << ": " << quantile_computation(0.25)
         << "   " << STAT_label[STATL_UPPER_QUARTILE] << ": " << quantile_computation(0.75);

#     ifdef DEBUG
      os << "   0.9 " << STAT_label[STATL_QUANTILE] << ": " << quantile_computation(0.9);
#     endif

    }
    os << endl;

    if ((shape) && (variance > 0.)) {
      if (comment_flag) {
        os << "# ";
      }
      os << STAT_label[STATL_SKEWNESS_COEFF] << ": " << skewness_computation() << "   "
         << STAT_label[STATL_KURTOSIS_COEFF] << ": " << kurtosis_computation() << endl;
    }
  }

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of the characteristics of a Reestimation object for a circular variable.
 *
 *  \param[in,out] os           stream,
 *  \param[in]     comment_flag flag comment.
 */
/*--------------------------------------------------------------*/

template <typename Type>
ostream& Reestimation<Type>::ascii_circular_characteristic_print(ostream &os ,
                                                                 bool comment_flag) const

{
  double mean_direction[4];


  os << STAT_label[STATL_SAMPLE_SIZE] << ": " << nb_element << endl;

  mean_direction_computation(mean_direction);

  if (comment_flag) {
    os << "# ";
  }
  os << STAT_label[STATL_MEAN_DIRECTION] << ": " << mean_direction[3];
  if (mean_direction[2] > 0.) {
    os << "   " << STAT_label[STATL_MEAN_RESULTANT_LENGTH] << ": " << mean_direction[2]
       << "   " << STAT_label[STATL_CIRCULAR_STANDARD_DEVIATION] << ": "
       << 180 * sqrt(-2 * log(mean_direction[2])) / M_PI;
  }
  os << endl;

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Display of a Reestimation object.
 *
 *  \param[in,out] os stream.
 */
/*--------------------------------------------------------------*/

template <typename Type>
ostream& Reestimation<Type>::print(ostream &os) const

{
  int i;


  os << endl;
  ascii_characteristic_print(os);

  os << "offset : " << offset;
  if (max > 0) {
    os << "   maximum : " << max;
  }
  os << endl;

  os << "frequencies (" << nb_value << " values) : ";
  for (i = 0;i < nb_value;i++) {
    os << frequency[i] << " ";
  }
  os << endl;

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the number of values from 0.
 */
/*--------------------------------------------------------------*/

template <typename Type>
void Reestimation<Type>::nb_value_computation()

{
  Type *pfrequency;


  pfrequency = frequency + alloc_nb_value;
  nb_value = alloc_nb_value;

  while ((*--pfrequency == 0) && (nb_value > 0)) {
    nb_value--;
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the number of values of null frequency from 0.
 */
/*--------------------------------------------------------------*/

template <typename Type>
void Reestimation<Type>::offset_computation()

{
  Type *pfrequency;


  pfrequency = frequency;
  offset = 0;

  while ((*pfrequency++ == 0) && (offset < nb_value - 1)) {
    offset++;
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the sample size.
 */
/*--------------------------------------------------------------*/

template <typename Type>
void Reestimation<Type>::nb_element_computation()

{
  int i;


  nb_element = 0;
  for (i = offset;i < nb_value;i++) {
    nb_element += frequency[i];
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Determination of the maximum frequency.
 */
/*--------------------------------------------------------------*/

template <typename Type>
void Reestimation<Type>::max_computation()

{
  int i;


  max = 0;
  for (i = offset;i < nb_value;i++) {
    if (frequency[i] > max) {
      max = frequency[i];
    }
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Determination of the mode.
 *
 *  \return mode.
 */
/*--------------------------------------------------------------*/

template <typename Type>
double Reestimation<Type>::mode_computation() const

{
  int i;
  int max_frequency;
  double mode;


  max_frequency = 0;
  for (i = offset;i < nb_value;i++) {
    if (frequency[i] > max_frequency) {
      max_frequency = frequency[i];
      mode = i;
    }
  }

  i = mode;
  while (frequency[i + 1] == frequency[i]) {
    i++;
  }
  if (i > mode) {
    mode = (i + mode) / 2.;
  }

  return mode;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Mean computation.
 */
/*--------------------------------------------------------------*/

template <typename Type>
void Reestimation<Type>::mean_computation()

{
  if (nb_element > 0) {
    int i;


    mean = 0.;
    for (i = offset;i < nb_value;i++) {
      mean += frequency[i] * i;
    }
    mean /= nb_element;
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of a quantile.
 *
 *  \param[in] icumul value of the cumulative distribution function.
 *
 *  \return           quantile.
 */
/*--------------------------------------------------------------*/

template <typename Type>
double Reestimation<Type>::quantile_computation(double icumul) const

{
  int i;
  double cumul , quantile;


  cumul = 0.;
  for (i = offset;i < nb_value;i++) {
    cumul += frequency[i] / (double)nb_element;
    if (cumul >= icumul) {
      quantile = i;
      if (cumul == icumul) {
        quantile += 0.5;
      }
      break;
    }
  }

  return quantile;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Variance computation.
 *
 *  \param[in] bias flag bias.
 */
/*--------------------------------------------------------------*/

template <typename Type>
void Reestimation<Type>::variance_computation(bool bias)

{
  if (mean != D_DEFAULT) {
    if (nb_element > 1) {
      int i;
      double diff;
      long double square_sum;


      square_sum = 0.;
      for (i = offset;i < nb_value;i++) {
        diff = i - mean;
        square_sum += frequency[i] * diff * diff;
      }

      if (bias) {
        variance = square_sum / nb_element;
      }
      else {
        variance = square_sum / (nb_element - 1);
      }
    }

    else {
      variance = 0.;
    }
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the mean absolute deviation.
 *
 *  \param[in] location location measure (e.g. mean or median).
 *
 *  \return             mean absolute deviation.
 */
/*--------------------------------------------------------------*/

template <typename Type>
double Reestimation<Type>::mean_absolute_deviation_computation(double location) const

{
  int i;
  double mean_absolute_deviation = D_DEFAULT;


  if (nb_element > 0) {
    mean_absolute_deviation = 0.;
    for (i = offset;i < nb_value;i++) {
      mean_absolute_deviation += frequency[i] * fabs(i - location);
    }
    mean_absolute_deviation /= nb_element;
  }

  return mean_absolute_deviation;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the log of the geometric mean.
 *
 *  \return log geometric mean.
 */
/*--------------------------------------------------------------*/

template <typename Type>
double Reestimation<Type>::log_geometric_mean_computation() const

{
  int i;
  double log_geometric_mean = D_DEFAULT;


  if (nb_element - frequency[0] > 0) {
    log_geometric_mean = 0.;
    for (i = MAX(offset , 2);i < nb_value;i++) {
      log_geometric_mean += frequency[i] * log(i);
    }
    log_geometric_mean /= (nb_element - frequency[0]);
  }

  return log_geometric_mean;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the coefficient of skewness.
 *
 *  \return coefficient of skewness.
 */
/*--------------------------------------------------------------*/

template <typename Type>
double Reestimation<Type>::skewness_computation() const

{
  int i;
  double skewness = D_INF , diff;
  long double cube_sum;


  if ((mean != D_DEFAULT) && (variance != D_DEFAULT)) {
    if ((nb_element > 2) && (variance > 0.)) {
      cube_sum = 0.;
      for (i = offset;i < nb_value;i++) {
        diff = i - mean;
        cube_sum += frequency[i] * diff * diff * diff;
      }
      skewness = cube_sum * nb_element /
                 ((nb_element - 1) * (double)(nb_element - 2) * pow(variance , 1.5));
    }

    else {
      skewness = 0.;
    }
  }

  return skewness;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the excess kurtosis:
 *         excess kurtosis = coefficient of kurtosis - 3.
 *
 *  \return excess kurtosis.
 */
/*--------------------------------------------------------------*/

template <typename Type>
double Reestimation<Type>::kurtosis_computation() const

{
  int i;
  double kurtosis = D_INF , diff;
  long double power_sum;


  if ((mean != D_DEFAULT) && (variance != D_DEFAULT)) {
    if (variance > 0.) {
      power_sum = 0.;
      for (i = offset;i < nb_value;i++) {
        diff = i - mean;
        power_sum += frequency[i] * diff * diff * diff * diff;
      }
      kurtosis = power_sum / ((nb_element - 1) * variance * variance) - 3.;
    }

    else {
      kurtosis = -2.;
    }
  }

  return kurtosis;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the mean direction for a circular variable.
 *
 *  \param[in] mean_direction pointer on the mean direction.
 */
/*--------------------------------------------------------------*/

template <typename Type>
void Reestimation<Type>::mean_direction_computation(double *mean_direction) const

{
  int i , j;


  mean_direction[0] = 0.;
  mean_direction[1] = 0.;

  for (i = offset;i < nb_value;i++) {
    mean_direction[0] += frequency[i] * cos(i * M_PI / 180);
    mean_direction[1] += frequency[i] * sin(i * M_PI / 180);
  }

  mean_direction[0] /= nb_element;
  mean_direction[1] /= nb_element;

  mean_direction[2] = sqrt(mean_direction[0] * mean_direction[0] +
                           mean_direction[1] * mean_direction[1]);

  if (mean_direction[2] > 0.) {
    mean_direction[3] = atan(mean_direction[1] / mean_direction[0]);

    if (mean_direction[0] < 0.) {
      mean_direction[3] += M_PI;
    }
    else if (mean_direction[1] < 0.) {
      mean_direction[3] += 2 * M_PI;
    }

    mean_direction[3] *= (180 / M_PI);
  }

  else {
    mean_direction[3] = D_DEFAULT;
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the information quantity.
 *
 *  \return information quantity.
 */
/*--------------------------------------------------------------*/

template <typename Type>
double Reestimation<Type>::information_computation() const

{
  int i;
  double information = D_INF;


  if (nb_element > 0) {
    information = 0.;
    for (i = offset;i < nb_value;i++) {
      if (frequency[i] > 0) {
        information += frequency[i] * log((double)frequency[i] / (double)nb_element);
      }
    }
  }

  return information;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the log-likelihood of a discrete distribution for a sample.
 *
 *  \param[in] dist reference on a Distribution object.
 *
 *  \return         log-likelihood.
 */
/*--------------------------------------------------------------*/

template <typename Type>
double Reestimation<Type>::likelihood_computation(const Distribution &dist) const

{
  int i;
  double likelihood = 0.;


  if (nb_element > 0) {
    if ((offset < dist.offset) || (nb_value > dist.nb_value)) {
      likelihood = D_INF;
    }

    else {
      for (i = offset;i < nb_value;i++) {
        if (frequency[i] > 0) {
          if (dist.mass[i] > 0.) {
            likelihood += frequency[i] * log(dist.mass[i]);
          }
          else {
            likelihood = D_INF;
            break;
          }
        }
      }

      if ((likelihood != D_INF) && (dist.complement > 0.)) {
        likelihood -= nb_element * log(1. - dist.complement);
      }
    }
  }

  return likelihood;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Estimation of a categorical distribution on the basis of a frequency distribution.
 *
 *  \param[in] dist pointer on a Distribution object.
 */
/*--------------------------------------------------------------*/

template <typename Type>
void Reestimation<Type>::distribution_estimation(Distribution *dist) const

{
  if (nb_element > 0) {
    int i;


    dist->offset = offset;
    dist->nb_value = nb_value;

    for (i = 0;i < offset;i++) {
      dist->mass[i] = 0.;
    }
    for (i = offset;i < nb_value;i++) {
      dist->mass[i] = (double)frequency[i] / (double)nb_element;
    }

    dist->cumul_computation();

    dist->max = (double)max / (double)nb_element;
    dist->mean = mean;
    dist->variance_computation();
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Estimation of a discrete distribution on the basis of a frequency distribution
 *         using a penalized likelihood estimator.
 *
 *  \param[in] dist     pointer on a Distribution object,
 *  \param[in] weight   penalty weight,
 *  \param[in] pen_type penalty type (first- or second-order difference or entropy),
 *  \param[in] penalty  penalty,
 *  \param[in] outside  management of side effects (zero outside the support or
 *                      continuation of the distribution).
 */
/*--------------------------------------------------------------*/

template <typename Type>
void Reestimation<Type>::penalized_likelihood_estimation(Distribution *dist , double weight ,
                                                         penalty_type pen_type , double *penalty ,
                                                         side_effect outside) const

{
  if (nb_element > 0) {
    int i;
    int iter;
    double ratio , inf_ratio , sup_ratio , norm , inf_norm , sup_norm;


    dist->penalty_computation(weight , pen_type , penalty , outside);

#   ifdef DEBUG
    {
      switch (pen_type) {
      case FIRST_DIFFERENCE :
        norm = dist->first_difference_norm_computation();
        break;
      case SECOND_DIFFERENCE :
        norm = dist->second_difference_norm_computation();
        break;
      case ENTROPY :
        norm = dist->information_computation();
        break;
      }

      norm = -norm * weight;
      for (i = dist->offset;i < dist->nb_value;i++) {
        norm += frequency[i] * log(dist->mass[i]);
      }

      cout << "\nPenalizedLikelihood: " << norm << endl;
    }
#   endif

    // computation of the normalization constant

    inf_norm = 0.;

    do {
      inf_norm += 0.05 * nb_element;

      for (i = offset;i < nb_value;i++) {
        if (inf_norm + penalty[i] <= 0.) {
          break;
        }
      }
    }
    while ((i < nb_value) && (inf_norm < nb_element / 2.));

    sup_norm = 2. * nb_element;

    if ((pen_type != ENTROPY) && (nb_element + weight > sup_norm)) {
      sup_norm = nb_element + weight;
    }

#   ifdef DEBUG
    cout << "Initialization : " << inf_norm << " " << sup_norm;
#   endif

    inf_ratio = 0.;
    sup_ratio = 0.;

    for (i = offset;i < nb_value;i++) {
      if (sup_norm + penalty[i] > 0.) {
        inf_ratio += frequency[i] / (sup_norm + penalty[i]);
      }
      if (inf_norm + penalty[i] > 0.) {
        sup_ratio += frequency[i] / (inf_norm + penalty[i]);
      }
    }

    iter = 0;
    do {
      iter++;

      ratio = 0.;
      norm = (inf_norm + sup_norm) / 2.;

      for (i = offset;i < nb_value;i++) {
        if (norm + penalty[i] > 0.) {
          ratio += frequency[i] / (norm + penalty[i]);
        }

        else {

#         ifdef MESSAGE
          cout << "\nRATIO ERROR" << "   " << STAT_label[STATL_ITERATION] << " " << iter
               << "   norm: " << norm << "   ratio: " << ratio << " | " << i << " ("
               << inf_ratio << ", " << sup_ratio << ")" << endl;
#         endif

          break;
        }
      }

      if (i < nb_value) {
        break;
      }

      else {
        if (ratio < 1.) {
          inf_ratio = ratio;
          sup_norm = norm;
        }
        else {
          sup_ratio = ratio;
          inf_norm = norm;
        }
      }
    }
    while ((sup_ratio - inf_ratio > BISECTION_RATIO_THRESHOLD) && (iter < BISECTION_NB_ITER));

    if ((sup_ratio - inf_ratio <= BISECTION_RATIO_THRESHOLD) || (iter >= BISECTION_NB_ITER)) {
      norm = (inf_norm + sup_norm) / 2.;

#     ifdef DEBUG
      cout << "   norm : " << norm << " | " << nb_element << "  (" << iter << ")  ";
#     endif

      // distribution reestimation

      for (i = 0;i < offset;i++) {
        dist->mass[i] = 0.;
      }
      for (i = offset;i < nb_value;i++) {
        dist->mass[i] = frequency[i] / (norm + penalty[i]);
      }
      for (i = nb_value;i < dist->nb_value;i++) {
        dist->mass[i] = 0.;
      }

      dist->offset_computation();
      dist->nb_value_computation();
      dist->cumul_computation();
      dist->max_computation();
      dist->mean_computation();
      dist->variance_computation();
    }
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Estimation of the parameters of a shifted binomial distribution on the basis of
 *         a frequency distribution.
 *
 *  \param[in] dist               pointer on a DiscreteParametric object,
 *  \param[in] min_inf_bound      minimum lower bound of the support,
 *  \param[in] min_inf_bound_flag flag on the distribution shift.
 *
 *  \return                       maximized log-likelihood.
 */
/*--------------------------------------------------------------*/

template <typename Type>
double Reestimation<Type>::binomial_estimation(DiscreteParametric *dist , int min_inf_bound ,
                                               bool min_inf_bound_flag) const

{
  int i , j;
  int max_inf_bound , inf_bound , sup_bound , sup_sup_bound , est_sup_bound ,
      min_sup_bound , max_sup_bound;
  double shift_mean , probability , likelihood , inf_bound_likelihood , max_likelihood = D_INF;


  if (mean - min_inf_bound >= variance) {

    // computation of the interval tested on the lower bound of the support

    if (!min_inf_bound_flag) {
      max_inf_bound = min_inf_bound;
    }
    else {
      max_inf_bound = MIN(offset , (int)(mean - variance));
    }

    if (variance == 0.) {
      if (mean <= max_inf_bound) {
        max_likelihood = 0.;
        dist->init((int)mean , (int)mean + 1 , D_DEFAULT , 0.);
      }
      else {
        if (mean < dist->alloc_nb_value) {
          max_likelihood = 0.;
          dist->init(min_inf_bound , (int)mean , D_DEFAULT , 1.);
        }
      }
    }

    else {
      sup_sup_bound = MIN(dist->alloc_nb_value ,
                          (int)(nb_value * SAMPLE_NB_VALUE_COEFF)) - 1;

      // estimation for each possible lower and upper bounds of the probability
      // on the basis of the sample mean

#     ifdef DEBUG
//      cout << "\nbound range: " << min_inf_bound
//           << " | " << max_inf_bound << " (" << sup_sup_bound << ")" << endl;

      int k = 0;
#     endif

      for (i = min_inf_bound;i <= max_inf_bound;i++) {
        shift_mean = mean - i;

        if (shift_mean > variance) {
          est_sup_bound = (int)round(i + shift_mean * shift_mean / (shift_mean - variance));

          min_sup_bound = MAX(nb_value - 1 , i + 1);
          if (min_sup_bound < MIN(sup_sup_bound , est_sup_bound - SUP_BOUND_MARGIN)) {
            min_sup_bound = MIN(sup_sup_bound , est_sup_bound - SUP_BOUND_MARGIN);
          }
          max_sup_bound = MIN(sup_sup_bound , est_sup_bound + SUP_BOUND_MARGIN);
          if (max_sup_bound < min_sup_bound) {
            max_sup_bound = min_sup_bound;
          }

#         ifdef DEBUG
//          cout << min_sup_bound << " " << max_sup_bound << " | " << est_sup_bound << endl;

          k += max_sup_bound - min_sup_bound + 1;
#         endif

          inf_bound_likelihood = D_INF;
          dist->inf_bound = i;

          for (j = min_sup_bound;j <= max_sup_bound;j++) {
            dist->sup_bound = j;
            dist->probability = (mean - i) / (j - i);

            dist->binomial_computation(1 , STANDARD);
            likelihood = dist->likelihood_computation(*this);

#           ifdef DEBUG
//            cout << i << " " << j << " " << dist->probability << " " << likelihood << endl;
#           endif

            if (likelihood > max_likelihood) {
              max_likelihood = likelihood;
              inf_bound = dist->inf_bound;
              sup_bound = dist->sup_bound;
              probability = dist->probability;
            }

            if (likelihood > inf_bound_likelihood) {
              inf_bound_likelihood = likelihood;
            }
          }

          if (inf_bound_likelihood < max_likelihood) {
            break;
          }
        }

        else {
          break;
        }
      }

      // update of the estimated parameters

      if (max_likelihood != D_INF) {
        dist->init(inf_bound , sup_bound , D_DEFAULT , probability);
      }

#     ifdef DEBUG
//      cout << "\nnumber of cases: "
//           << (max_inf_bound - min_inf_bound + 1) * (2 * SUP_BOUND_MARGIN + 1)
//           << " | number of computations: " << k << endl;
#     endif
    }
  }

  return max_likelihood;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Estimation of the parameters of a shifted Poisson distribution on the basis of
 *         a frequency distribution.
 *
 *  \param[in] dist               pointer on a DiscreteParametric object,
 *  \param[in] min_inf_bound      minimum lower bound of the support,
 *  \param[in] min_inf_bound_flag flag on the distribution shift,
 *  \param[in] cumul_threshold    threshold on the cumulative distribution function.
 *
 *  \return                       maximized log-likelihood.
 */
/*--------------------------------------------------------------*/

template <typename Type>
double Reestimation<Type>::poisson_estimation(DiscreteParametric *dist , int min_inf_bound ,
                                              bool min_inf_bound_flag , double cumul_threshold) const

{
  int i;
  int max_inf_bound , inf_bound;
  double diff , parameter , likelihood , max_likelihood = D_INF;


  // computation of the interval tested on the lower bound of the support

  if (!min_inf_bound_flag) {
    max_inf_bound = min_inf_bound;
  }

  else {
    diff = mean - variance;
    if (diff < 0.) {
      diff--;
    }

    min_inf_bound = MAX(min_inf_bound , (int)round(diff) - INF_BOUND_MARGIN);
    max_inf_bound = MIN(offset , (int)round(diff) + INF_BOUND_MARGIN);
  }

  if (((mean - max_inf_bound) * POISSON_RATIO < variance) &&
      ((mean - min_inf_bound) / POISSON_RATIO > variance)) {

    // estimation for each possible lower bound of the parameter on the basis of the sample mean

#   ifdef DEBUG
//    cout << "\ninf bound range: " << min_inf_bound << " | " << max_inf_bound << endl;
#   endif

    for (i = max_inf_bound;i >= min_inf_bound;i--) {
      dist->inf_bound = i;
      dist->parameter = mean - i;

      dist->poisson_computation(nb_value , cumul_threshold , STANDARD);
      likelihood = dist->likelihood_computation(*this);

      if (likelihood > max_likelihood) {
        max_likelihood = likelihood;
        inf_bound = dist->inf_bound;
        parameter = dist->parameter;
      }

      else {
        if (likelihood < max_likelihood) {
          break;
        }
      }
    }

    // update of the estimated parameters

    if (max_likelihood != D_INF) {
      dist->init(inf_bound , I_DEFAULT , parameter , D_DEFAULT);
    }

#   ifdef DEBUG
//    cout << "\nnumber of cases: " << (max_inf_bound - min_inf_bound + 1)
//         << " | number of computations: " << (max_inf_bound - i + 1) << endl;
#   endif
  }

  return max_likelihood;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Estimation of the parameters of a shifted negative binomial distribution
 *         on the basis of a frequency distribution.
 *
 *  \param[in] dist               pointer on a DiscreteParametric object,
 *  \param[in] min_inf_bound      minimum lower bound of the support,
 *  \param[in] min_inf_bound_flag flag on the distribution shift,
 *  \param[in] cumul_threshold    threshold on the cumulative distribution function.
 *
 *  \return                       maximized log-likelihood.
 */
/*--------------------------------------------------------------*/

template <typename Type>
double Reestimation<Type>::negative_binomial_estimation(DiscreteParametric *dist , int min_inf_bound ,
                                                        bool min_inf_bound_flag , double cumul_threshold) const

{
  int i;
  int max_inf_bound , inf_bound;
  double diff , shift_mean , parameter , probability , likelihood , max_likelihood = D_INF;


  // computation of the interval tested on the lower bound of the support

  if (!min_inf_bound_flag) {
    max_inf_bound = min_inf_bound;
  }

  else {
    diff = mean - variance;
    if (diff < 0.) {
      diff--;
    }

    min_inf_bound = MAX(min_inf_bound , (int)diff + 1);
    max_inf_bound = offset;
  }

  if (mean - max_inf_bound < variance) {

    // estimation for each possible lower bound of the shape parameter and
    // the probability on the basis of the sample mean and variance

#   ifdef DEBUG
    cout << "\ninf bound range: " << min_inf_bound << " | " << max_inf_bound << endl;
#   endif

    for (i = max_inf_bound;i >= min_inf_bound;i--) {
      dist->inf_bound = i;

      shift_mean = mean - i;
      dist->parameter = shift_mean * shift_mean / (variance - shift_mean);
      dist->probability = shift_mean / variance;

#     ifdef DEBUG
//      cout << i << " : " dist->parameter << " | " << dist->probability << endl;
#     endif

      dist->negative_binomial_computation(nb_value , cumul_threshold , STANDARD);
      likelihood = dist->likelihood_computation(*this);

      if (likelihood > max_likelihood) {
        max_likelihood = likelihood;
        inf_bound = dist->inf_bound;
        parameter = dist->parameter;
        probability = dist->probability;
      }

      else {
        if (likelihood < max_likelihood) {
          break;
        }
      }
    }

    // update of the estimated parameters

    if (max_likelihood != D_INF) {
      dist->init(inf_bound , I_DEFAULT , parameter , probability);
    }

#   ifdef DEBUG
//    cout << "\nnumber of cases : " << max_inf_bound - min_inf_bound + 1
//         << " | number of computations : " << (max_inf_bound - i + 1) << endl;
#   endif
  }

  return max_likelihood;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Estimation of the parameters of a shifted Poisson geometric distribution
 *         on the basis of a frequency distribution.
 *
 *  \param[in] dist               pointer on a DiscreteParametric object,
 *  \param[in] min_inf_bound      minimum lower bound of the support,
 *  \param[in] min_inf_bound_flag flag on the distribution shift,
 *  \param[in] cumul_threshold    threshold on the cumulative distribution function.
 *
 *  \return                       maximized log-likelihood.
 */
/*--------------------------------------------------------------*/

template <typename Type>
double Reestimation<Type>::poisson_geometric_estimation(DiscreteParametric *dist , int min_inf_bound ,
                                                        bool min_inf_bound_flag , double cumul_threshold) const

{
  int i;
  int max_inf_bound , inf_bound;
  double diff , buff , parameter , probability , likelihood , max_likelihood = D_INF;


  // computation of the interval tested on the lower bound of the support

  if (!min_inf_bound_flag) {
    max_inf_bound = min_inf_bound;
  }
  else {
    max_inf_bound = MIN(offset , (int)(mean * mean / (variance + mean)));
  }

  if (mean - max_inf_bound < variance) {

    // estimation for each possible lower bound of the shape parameter and
    // the probability on the basis of the sample mean and variance

#   ifdef DEBUG
    cout << "\ninf bound range: " << min_inf_bound << " | " << max_inf_bound << endl;
#   endif

    for (i = max_inf_bound;i >= min_inf_bound;i--) {
      dist->inf_bound = i;
      buff = mean * mean - (variance + mean) * i;
      dist->parameter = (buff + mean * sqrt(buff)) / (variance + mean);
      dist->probability = (i + dist->parameter) / mean;

#     ifdef DEBUG
//      cout << i << " : " dist->parameter << " | " << dist->probability << endl;
#     endif

      dist->poisson_geometric_computation(nb_value , cumul_threshold);
      likelihood = dist->likelihood_computation(*this);

      if (likelihood > max_likelihood) {
        max_likelihood = likelihood;
        inf_bound = dist->inf_bound;
        parameter = dist->parameter;
        probability = dist->probability;
      }

      else {
        if (likelihood < max_likelihood) {
          break;
        }
      }
    }

    // update of the estimated parameters

    if (max_likelihood != D_INF) {
      dist->init(inf_bound , I_DEFAULT , parameter , probability);
    }

#   ifdef DEBUG
//    cout << "\nnumber of cases: " << max_inf_bound - min_inf_bound + 1
//         << " | number of computations: " << (max_inf_bound - i + 1) << endl;
#   endif
  }

  return max_likelihood;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Estimation of the parameters of a discrete parametric distribution
 *         (binomial, Poisson, negative binomial, Poisson geometric) on the basis of
 *         a frequency distribution.
 *
 *  \param[in] dist               pointer on a DiscreteParametric object,
 *  \param[in] min_inf_bound      minimum lower bound of the support,
 *  \param[in] min_inf_bound_flag flag on the distribution shift,
 *  \param[in] cumul_threshold    threshold on the cumulative distribution function,
 *  \param[in] poisson_geometric  flag on the estimation of a Poisson geometric distribution.
 *
 *  \return                       maximized log-likelihood.
 */
/*--------------------------------------------------------------*/

template <typename Type>
double Reestimation<Type>::parametric_estimation(DiscreteParametric *dist , int min_inf_bound ,
                                                 bool min_inf_bound_flag , double cumul_threshold ,
                                                 bool poisson_geometric) const

{
  double likelihood;


  switch (dist->ident) {
  case BINOMIAL :
    likelihood = binomial_estimation(dist , min_inf_bound , min_inf_bound_flag);
    break;
  case POISSON :
    likelihood = poisson_estimation(dist , min_inf_bound , min_inf_bound_flag , cumul_threshold);
    break;
  case NEGATIVE_BINOMIAL :
    likelihood = negative_binomial_estimation(dist , min_inf_bound , min_inf_bound_flag , cumul_threshold);
    break;
  case POISSON_GEOMETRIC :
    likelihood = poisson_geometric_estimation(dist , min_inf_bound , min_inf_bound_flag , cumul_threshold);
    break;
  }

  return likelihood;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Estimation of the type and parameters of a discrete parametric distribution
 *         (binomial, Poisson, negative binomial, Poisson geometric) on the basis of
 *         a frequency distribution.
 *
 *  \param[in] dist               pointer on a DiscreteParametric object,
 *  \param[in] min_inf_bound      minimum lower bound of the support,
 *  \param[in] min_inf_bound_flag flag on the distribution shift,
 *  \param[in] cumul_threshold    threshold on the cumulative distribution function,
 *  \param[in] poisson_geometric  flag on the estimation of a Poisson geometric distribution.
 *
 *  \return                       maximized log-likelihood.
 */
/*--------------------------------------------------------------*/

template <typename Type>
double Reestimation<Type>::type_parametric_estimation(DiscreteParametric *dist , int min_inf_bound ,
                                                      bool min_inf_bound_flag , double cumul_threshold ,
                                                      bool poisson_geometric) const

{
  double likelihood , max_likelihood;
  DiscreteParametric *bdist;


  bdist = new DiscreteParametric(dist->alloc_nb_value);

  max_likelihood = binomial_estimation(dist , min_inf_bound , min_inf_bound_flag);
  if (max_likelihood != D_INF) {
    dist->ident = BINOMIAL;
  }

# ifdef DEBUG
//  max_likelihood = D_INF;   for the earthquake data (Durand & Guedon, 2016)
//  likelihood = poisson_estimation(bdist , 0 , false , cumul_threshold);
# endif

  likelihood = poisson_estimation(bdist , min_inf_bound , min_inf_bound_flag , cumul_threshold);
  if (likelihood > max_likelihood) {
    bdist->ident = POISSON;
    max_likelihood = likelihood;
    dist->equal_size_copy(*bdist);
    dist->copy(*bdist);
  }

  if ((min_inf_bound == 0) || (!poisson_geometric)) {  // for comparing negative binomial and poisson geometric distributions
  likelihood = negative_binomial_estimation(bdist , min_inf_bound , min_inf_bound_flag , cumul_threshold);
  if (likelihood > max_likelihood) {
    bdist->ident = NEGATIVE_BINOMIAL;
    max_likelihood = likelihood;
    dist->equal_size_copy(*bdist);
    dist->copy(*bdist);
  }
  }

  if ((min_inf_bound > 0) && (poisson_geometric)) {
    likelihood = poisson_geometric_estimation(bdist , min_inf_bound , min_inf_bound_flag , cumul_threshold);
    if (likelihood > max_likelihood) {
      bdist->ident = POISSON_GEOMETRIC;
      max_likelihood = likelihood;
      dist->equal_size_copy(*bdist);
      dist->copy(*bdist);
    }
  }

  delete bdist;

  return max_likelihood;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Estimation of the type and parameters of a discrete parametric distribution
 *         (binomial, Poisson, negative binomial) on the basis of a frequency distribution.
 *
 *  \param[in] min_inf_bound      minimum lower bound of the support,
 *  \param[in] min_inf_bound_flag flag on the distribution shift,
 *  \param[in] cumul_threshold    threshold on the cumulative distribution function.
 *
 *  \return                       discrete parametric distribution.
 */
/*--------------------------------------------------------------*/

template <typename Type>
DiscreteParametric* Reestimation<Type>::type_parametric_estimation(int min_inf_bound , bool min_inf_bound_flag ,
                                                                   double cumul_threshold) const

{
  double likelihood;
  DiscreteParametric *dist;


  // construction of a DiscreteParametric object

  dist = new DiscreteParametric((int)(nb_value * SAMPLE_NB_VALUE_COEFF));

  // parameter estimation

  likelihood = type_parametric_estimation(dist , min_inf_bound , min_inf_bound_flag , cumul_threshold);

  // computation of the probability mass function

  if (likelihood != D_INF) {
    dist->computation(nb_value , cumul_threshold);
  }
  else {
    delete dist;
    dist = NULL;
  }

  return dist;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Combination of the reestimation quantities of the basis time interval distribution and
 *         the length-biased distribution (equilibrium stochastic process).
 *
 *  \param[in] length_bias_reestim pointer on the reestimation quantities of the length-biased distribution,
 *  \param[in] imean               distribution mean.
 */
/*--------------------------------------------------------------*/

template <typename Type>
void Reestimation<Type>::equilibrium_process_combination(const Reestimation<Type> *length_bias_reestim ,
                                                         double imean)

{
  if (nb_element + length_bias_reestim->nb_element > 0) {
    int i;


    for (i = offset;i < nb_value;i++) {
      frequency[i] = (frequency[i] + length_bias_reestim->frequency[i]) *
                     (nb_element + length_bias_reestim->nb_element) /
                     (nb_element + length_bias_reestim->nb_element * i / imean);
    }

    nb_element_computation();
    max_computation();
    mean_computation();
    variance_computation();
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Estimation of a discrete distribution on the basis of the reestimation quantities of
 *         the basis time interval distribution and the length-biased distribution (equilibrium stochastic process).
 *
 *  \param[in] length_bias_reestim pointer on the reestimation quantities of the length-biased distribution,
 *  \param[in] dist                distribution,
 *  \param[in] imean               distribution mean.
 */
/*--------------------------------------------------------------*/

template <typename Type>
void Reestimation<Type>::equilibrium_process_estimation(const Reestimation<Type> *length_bias_reestim ,
                                                        Distribution *dist , double imean) const

{
  if (nb_element + length_bias_reestim->nb_element > 0) {
    int i;


    for (i = 0;i < offset;i++) {
      dist->mass[i] = 0.;
    }
    for (i = offset;i < nb_value;i++) {
      dist->mass[i] = (frequency[i] + length_bias_reestim->frequency[i]) /
                      (nb_element + length_bias_reestim->nb_element * i / imean);
    }
    for (i = nb_value;i < dist->nb_value;i++) {
      dist->mass[i] = 0.;
    }

    dist->offset_computation();
    dist->nb_value_computation();
    dist->cumul_computation();

    // renormalization of the distribution

    for (i = dist->offset;i < dist->nb_value;i++) {
      dist->mass[i] /= dist->cumul[nb_value - 1];
    }

    dist->cumul_computation();
    dist->max_computation();
    dist->mean_computation();
    dist->variance_computation();
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Estimation of a discrete distribution on the basis of the reestimation quantities of
 *         the basis time interval distribution and the length-biased distribution using
 *         a penalized likelihood estimator (equilibrium stochastic process).
 *
 *  \param[in] length_bias_reestim pointer on the reestimation quantities of the length-biased distribution,
 *  \param[in] dist                distribution,
 *  \param[in] imean               distribution mean,
 *  \param[in] weight              penalty weight ,
 *  \param[in] pen_type            penalty type (first- or second-order difference or entropy),
 *  \param[in] penalty             penalty,
 *  \param[in] outside             management of side effects (zero outside the support or
 *                                 continuation of the distribution).
 */
/*--------------------------------------------------------------*/

template <typename Type>
void Reestimation<Type>::penalized_likelihood_equilibrium_process_estimation(const Reestimation<Type> *length_bias_reestim ,
                                                                             Distribution *dist , double imean ,
                                                                             double weight , penalty_type pen_type ,
                                                                             double *penalty , side_effect outside) const

{
  if (nb_element + length_bias_reestim->nb_element > 0) {
    int i;
    int iter;
    double ratio , inf_ratio , sup_ratio , norm , inf_norm , sup_norm;


    dist->penalty_computation(weight , pen_type , penalty , outside);

    // computation of the normalization constant

    inf_norm = 0.;

    do {
      inf_norm += 0.05 * nb_element;

      for (i = offset;i < nb_value;i++) {
        if (inf_norm + length_bias_reestim->nb_element * i / imean + penalty[i] <= 0.) {
          break;
        }
      }
    }
    while ((i < nb_value) && (inf_norm < nb_element / 2.));

    sup_norm = 2. * nb_element;

    if ((pen_type != ENTROPY) && (nb_element + weight > sup_norm)) {
      sup_norm = nb_element + weight;
    }

    inf_ratio = 0.;
    sup_ratio = 0.;

    for (i = offset;i < nb_value;i++) {
      if (sup_norm + length_bias_reestim->nb_element * i / imean + penalty[i] > 0.) {
        inf_ratio += (frequency[i] + length_bias_reestim->frequency[i]) /
                     (sup_norm + length_bias_reestim->nb_element * i / imean + penalty[i]);
      }
      if (inf_norm + length_bias_reestim->nb_element * i / imean + penalty[i] > 0.) {
        sup_ratio += (frequency[i] + length_bias_reestim->frequency[i]) /
                     (inf_norm + length_bias_reestim->nb_element * i / imean + penalty[i]);
      }
    }

    iter = 0;
    do {
      iter++;

      ratio = 0.;
      norm = (inf_norm + sup_norm) / 2.;

      for (i = offset;i < nb_value;i++) {
        if (norm + length_bias_reestim->nb_element * i / imean + penalty[i] > 0.) {
          ratio += (frequency[i] + length_bias_reestim->frequency[i]) /
                   (norm + length_bias_reestim->nb_element * i / imean + penalty[i]);
        }

        else {

#         ifdef MESSAGE
          cout << "\nRATIO ERROR" << "   " << STAT_label[STATL_ITERATION] << " " << iter
               << "   norm: " << norm << "   ratio: " << ratio << " | " << i << " ("
               << inf_ratio << ", " << sup_ratio << ")" << endl;
#         endif

          break;
        }
      }

      if (i < nb_value) {
        break;
      }

      else {
        if (ratio < 1.) {
          inf_ratio = ratio;
          sup_norm = norm;
        }
        else {
          sup_ratio = ratio;
          inf_norm = norm;
        }
      }
    }
    while ((sup_ratio - inf_ratio > BISECTION_RATIO_THRESHOLD) && (iter < BISECTION_NB_ITER));

    if ((sup_ratio - inf_ratio <= BISECTION_RATIO_THRESHOLD) || (iter >= BISECTION_NB_ITER)) {
      norm = (inf_norm + sup_norm) / 2.;

      // distribution reestimation

      for (i = 0;i < offset;i++) {
        dist->mass[i] = 0.;
      }
      for (i = offset;i < nb_value;i++) {
        dist->mass[i] = (frequency[i] + length_bias_reestim->frequency[i]) /
                        (norm + length_bias_reestim->nb_element * i / imean + penalty[i]);
      }
      for (i = nb_value;i < dist->nb_value;i++) {
        dist->mass[i] = 0.;
      }

      dist->offset_computation();
      dist->nb_value_computation();
      dist->cumul_computation();
      dist->max_computation();
      dist->mean_computation();
      dist->variance_computation();
    }
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Estimation of a state occupancy distribution using the Kaplan-Meier estimator.
 *
 *  \param[in] final_run                   pointer on the right-censored sojourn times,
 *  \param[in] occupancy_reestim           pointer on the reestimation quantities,
 *  \param[in] occupancy_survivor          pointer on the survival function corresponding to the complete sojourn times,
 *  \param[in] censored_occupancy_survivor pointer on the survival function corresponding to the right-censored sojourn times,
 *  \param[in] characteristic_computation  flag for the computation of the characteristics of the reestimation quantities.
 */
/*--------------------------------------------------------------*/

template <typename Type>
void Reestimation<Type>::state_occupancy_estimation(const Reestimation<Type> *final_run ,
                                                    Reestimation<double> *occupancy_reestim ,
                                                    Type *occupancy_survivor ,
                                                    Type *censored_occupancy_survivor ,
                                                    bool characteristic_computation)

{
  int i;
  int max_nb_value;
  double hazard_rate , hazard_product;


  if (nb_value > 0) {
    occupancy_survivor[nb_value - 1] = frequency[nb_value - 1];
    for (i = nb_value - 2;i >= 1;i--) {
      occupancy_survivor[i] = occupancy_survivor[i + 1] + frequency[i];
    }
  }

  max_nb_value = MAX(nb_value + 1 , final_run->nb_value);
  for (i = max_nb_value - 1;i >= final_run->nb_value;i--) {
    censored_occupancy_survivor[i] = 0;
  }
  censored_occupancy_survivor[final_run->nb_value - 1] = final_run->frequency[final_run->nb_value - 1];
  for (i = final_run->nb_value - 2;i >= 2;i--) {
    censored_occupancy_survivor[i] = censored_occupancy_survivor[i + 1] + final_run->frequency[i];
  }

  hazard_product = 1.;
  for (i = 1;i < nb_value;i++) {
    if (occupancy_survivor[i] + censored_occupancy_survivor[i + 1] > 0) {
      hazard_rate = (double)frequency[i] / (double)(occupancy_survivor[i] +
                     censored_occupancy_survivor[i + 1]);
      occupancy_reestim->frequency[i] = hazard_rate * hazard_product * (nb_element +
                                        final_run->nb_element);
      hazard_product *= (1. - hazard_rate);
    }

    else {
      break;
    }
  }

  if (characteristic_computation) {
    occupancy_reestim->nb_value = i;
  }

  if (i == nb_value) {
    for (i = nb_value;i < final_run->nb_value - 1;i++) {
      occupancy_reestim->frequency[i] = 0.;
    }
    occupancy_reestim->frequency[final_run->nb_value - 1] = hazard_product *
                                                            (nb_element + final_run->nb_element);

    if ((characteristic_computation) &&
        (occupancy_reestim->frequency[final_run->nb_value - 1] > 0)) {
      occupancy_reestim->nb_value = final_run->nb_value;
    }
  }

  if (characteristic_computation) {
    occupancy_reestim->offset = offset;
    occupancy_reestim->nb_element_computation();
    occupancy_reestim->max_computation();
    occupancy_reestim->mean_computation();
    occupancy_reestim->variance_computation();
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Estimation of a gamma distribution on the basis of a frequency distribution.
 *
 *  \param[in] dist continuous distribution,
 *  \param[in] iter EM iteration.
 */
/*--------------------------------------------------------------*/

template <typename Type>
void Reestimation<Type>::gamma_estimation(ContinuousParametric *dist , int iter) const

{
//  int i;
  double buff , log_geometric_mean , diff;
//  double bvariance;


  if (frequency[0] / nb_element > GAMMA_ZERO_FREQUENCY_THRESHOLD) {
    dist->shape = 0;
    dist->scale = D_DEFAULT;
  }

  else {
    if (variance > 0.) {
/*      if (sqrt(variance) < mean * GAMMA_VARIATION_COEFF_THRESHOLD) {
        bvariance = mean * mean * GAMMA_VARIATION_COEFF_THRESHOLD * GAMMA_VARIATION_COEFF_THRESHOLD;
      }
      else {
        bvariance = variance;
      }

      dist->shape = mean * mean / bvariance;
      dist->scale = bvariance / mean; */

      // Hawang & Huang (2012), Ann. Inst. Statist. Math. 54(4), 840-847

      buff = mean * mean / variance;
      if (buff > GAMMA_INVERSE_SAMPLE_SIZE_FACTOR / (double)nb_element) {
        dist->shape = buff - 1. / (double)nb_element;
      }
      else {
        dist->shape = buff;
      }
/*      if (dist->shape < GAMMA_MIN_SHAPE_PARAMETER) {
        dist->shape = GAMMA_MIN_SHAPE_PARAMETER;
      } */
      dist->scale = mean / dist->shape;

      if ((dist->shape >= GAMMA_SHAPE_PARAMETER_THRESHOLD) &&
          (nb_element < GAMMA_FREQUENCY_THRESHOLD)) {
        log_geometric_mean = log_geometric_mean_computation();

/*        i = 0;  to be investigated

#       ifdef DEBUG
        cout << "\n" << STAT_word[STATW_SHAPE] << " : " << dist->shape << "   "
             << STAT_word[STATW_SCALE] << " : " << dist->scale << endl;
#       endif

        do {
          dist->scale = exp(log_geometric_mean - digamma(dist->shape));
          dist->shape = mean / dist->scale;
          i++;

#         ifdef DEBUG
          cout << STAT_word[STATW_SHAPE] << " : " << dist->shape  << "   "
               << STAT_word[STATW_SCALE] << " : " << dist->scale << endl;
#         endif

        }
        while (i < MIN(GAMMA_ITERATION_FACTOR * iter , GAMMA_MAX_NB_ITERATION)); */

        // approximations Johnson, Kotz & Balakrishnan, Continuous Univariate Distributions, vol. 1, 2nd ed., pp. 361-362

//        dist->shape = mean / (2 * (mean - exp(log_geometric_mean))) - 1./12.;
        diff = log(mean) - log_geometric_mean;
        dist->shape = (1 + sqrt(1 + 4 * diff / 3)) / (4 * diff);
        dist->scale = mean / dist->shape;
      }
    }

    else {
      dist->shape = GAMMA_MIN_SHAPE_PARAMETER;
      dist->scale = GAMMA_DEFAULT_SCALE_PARAMETER;
    }
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Estimation of a zero-inflated gamma distribution on the basis of
 *         a frequency distribution.
 *
 *  \param[in] dist continuous distribution,
 *  \param[in] iter EM iteration.
 */
/*--------------------------------------------------------------*/

template <typename Type>
void Reestimation<Type>::zero_inflated_gamma_estimation(ContinuousParametric *dist , int iter) const

{
  if (frequency[0] / nb_element > GAMMA_ZERO_FREQUENCY_THRESHOLD) {
    dist->zero_probability = 1.;
    dist->shape = D_DEFAULT;
    dist->scale = D_DEFAULT;
  }

  else {
    int i;
    double bmean , bvariance , buff , diff , log_geometric_mean;


    dist->zero_probability = frequency[0] / nb_element;

    bmean = 0.;
    for (i = MIN(offset , 1);i < nb_value;i++) {
      bmean += frequency[i] * i;
    }
    bmean /= (nb_element - frequency[0]);

    bvariance = 0.;

    if (nb_element - frequency[0] > 0) {
      for (i = MIN(offset , 1);i < nb_value;i++) {
        diff = i - bmean;
        bvariance += frequency[i] * diff * diff;
      }
      bvariance /= (nb_element - frequency[0]);
    }

    if (bvariance > 0.) {
/*      if (sqrt(bvariance) < bmean * GAMMA_VARIATION_COEFF_THRESHOLD) {
        bvariance = bmean * bmean * GAMMA_VARIATION_COEFF_THRESHOLD * GAMMA_VARIATION_COEFF_THRESHOLD;
      }
      dist->shape = bmean * bmean / bvariance;
      dist->scale = bvariance / bmean; */

      // Hawang & Huang (2012), Ann. Inst. Statist. Math. 54(4), 840-847

      buff = bmean * bmean / bvariance;
      if (buff > GAMMA_INVERSE_SAMPLE_SIZE_FACTOR / (double)(nb_element - frequency[0])) {
        dist->shape = buff - 1. / (double)(nb_element - frequency[0]);
      }
      else {
        dist->shape = buff;
      }
/*      if (dist->shape < GAMMA_MIN_SHAPE_PARAMETER) {
        dist->shape = GAMMA_MIN_SHAPE_PARAMETER;
      } */
      dist->scale = bmean / dist->shape;

      if ((dist->shape >= GAMMA_SHAPE_PARAMETER_THRESHOLD) &&
          (nb_element - frequency[0] < GAMMA_FREQUENCY_THRESHOLD)) {
        log_geometric_mean = log_geometric_mean_computation();
/*        i = 0;   a revoir

#       ifdef DEBUG
        cout << "\n" << STAT_word[STATW_SHAPE] << " : " << dist->shape << "   "
             << STAT_word[STATW_SCALE] << " : " << dist->scale << endl;
#       endif

        do {
          dist->scale = exp(log_geometric_mean - digamma(dist->shape));
          dist->shape = bmean / dist->scale;
          i++;

#         ifdef DEBUG
          cout << STAT_word[STATW_SHAPE] << " : " << dist->shape  << "   "
               << STAT_word[STATW_SCALE] << " : " << dist->scale << endl;
#         endif

        }
        while (i < MIN(GAMMA_ITERATION_FACTOR * iter , GAMMA_MAX_NB_ITERATION)); */

        // approximations Johnson, Kotz & Balakrishnan, Continuous Univariate Distributions, vol. 1, 2nd ed., pp. 361-362

//        dist->shape = bmean / (2 * (bmean - exp(log_geometric_mean))) - 1./12.;
        diff = log(bmean) - log_geometric_mean;
        dist->shape = (1 + sqrt(1 + 4 * diff / 3)) / (4 * diff);
        dist->scale = bmean / dist->shape;
      }
    }

    else {
      dist->shape = GAMMA_MIN_SHAPE_PARAMETER;
      dist->scale = GAMMA_DEFAULT_SCALE_PARAMETER;
    }
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Estimation of an inverse Gaussian distribution on the basis of
 *         a frequency distribution.
 *
 *  \param[in] dist continuous distribution.
 */
/*--------------------------------------------------------------*/

template <typename Type>
void Reestimation<Type>::inverse_gaussian_estimation(ContinuousParametric *dist) const

{
  int i;


  if (mean > 0.) {
    dist->location = mean;

    dist->scale = 0.;
    for (i = MIN(offset , 1);i < nb_value;i++) {
      dist->scale += frequency[i] * (1. / (double)i - 1. / mean);
    }
    dist->scale = nb_element / dist->scale;
  }

  else {
    dist->location = D_DEFAULT;
    dist->scale = D_DEFAULT;
  }
}


};  // namespace stat_tool



#endif
