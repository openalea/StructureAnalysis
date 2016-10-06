/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       V-Plants: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2016 CIRAD/INRA/Inria Virtual Plants
 *
 *       File author(s): Yann Guedon (yann.guedon@cirad.fr)
 *
 *       $Source$
 *       $Id: distribution_algorithms.cpp 18450 2015-07-29 09:42:43Z guedon $
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
#include <cstdlib>

#include <boost/math/distributions/binomial.hpp>
#include <boost/math/distributions/poisson.hpp>
#include <boost/math/distributions/negative_binomial.hpp>

#include "tool/config.h"

#include "distribution.h"
#include "stat_label.h"

using namespace std;
using namespace boost::math;


namespace stat_tool {



/*--------------------------------------------------------------*/
/**
 *  \brief Convolution of 2 distributions (the convolution can be put in one of the 2 distributions).
 *
 *  \param[in] dist1     reference on the 1st distribution,
 *  \param[in] dist2     reference on the 2nd distribution,
 *  \param[in] inb_value number of values of the convolution.
 */
/*--------------------------------------------------------------*/

void Distribution::convolution(Distribution &dist1 , Distribution &dist2 , int inb_value)

{
  register int i , j;
  int coffset , cnb_value , min , max;
  double sum , *pmass1 , *pmass2 , *pmass;


  cnb_value = MIN(dist1.nb_value + dist2.nb_value - 1 , alloc_nb_value);
  if ((inb_value != I_DEFAULT) && (inb_value < cnb_value)) {
    cnb_value = inb_value;
  }
  coffset = MIN(dist1.offset + dist2.offset , cnb_value - 1);

  pmass = mass + cnb_value;

  for (i = cnb_value - 1;i >= coffset;i--) {
    sum = 0.;
    min = MAX(dist1.offset , i - (dist2.nb_value - 1));
    max = MIN(dist1.nb_value - 1 , i - dist2.offset);

    if (max >= min) {
      pmass1 = dist1.mass + min;
      pmass2 = dist2.mass + i - min;

      for (j = min;j <= max;j++) {
        sum += *pmass1++ * *pmass2--;
      }
    }

    *--pmass = sum;
  }

  for (i = coffset - 1;i >= 0;i--) {
    *--pmass = 0.;
  }

  offset = coffset;
  nb_value = cnb_value;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the probability mass function of a binomial distribution.
 *
 *  \param[in] inb_value number of values,
 *  \param[in] mode      computation mode (STANDARD/RENEWAL).
 */
/*--------------------------------------------------------------*/

void DiscreteParametric::binomial_computation(int inb_value , distribution_computation mode)

{
  register int i;
  int set , subset;
  double failure = 1. - probability , success = probability , ratio , scale , term;


  switch (mode) {
  case STANDARD :
    offset = inf_bound;
    nb_value = sup_bound + 1;
    break;
  case RENEWAL :
    offset = MIN(inf_bound , inb_value - 1);
    nb_value = MIN(sup_bound + 1 , inb_value);
    break;
  }

  // null probability values before the lower bound of the support

  for (i = 0;i < MIN(nb_value , inf_bound);i++) {
    mass[i] = 0.;
  }

  if (nb_value > inf_bound) {
    set = sup_bound - inf_bound;

    if (probability <= B_PROBABILITY) {
      subset = 0;

      // case direct computation

      if ((sup_bound - inf_bound) / failure < B_THRESHOLD) {

        // computation of the lower bound probability

        term = 1.;
        for (i = 0;i < sup_bound - inf_bound;i++) {
          term *= failure;
        }
        mass[inf_bound] = term;

        // computation of the probabilities for the successive values (forward recurrence)

        ratio = success / failure;

        for (i = inf_bound + 1;i < nb_value;i++) {
          scale = (double)(set - subset) / (double)(subset + 1);
          subset++;
          term *= scale * ratio;
          mass[i] = term;
        }
      }

      // case computation in log

      else {

        // computation of the lower bound probability

        term = (sup_bound - inf_bound) * log(failure);
        mass[inf_bound] = exp(term);

        // computation of the probabilities for the successive values (forward recurrence)

        ratio = log(success / failure);

        for (i = inf_bound + 1;i < nb_value;i++) {
          scale = (double)(set - subset) / (double)(subset + 1);
          subset++;
          term += log(scale) + ratio;
          mass[i] = exp(term);
        }
      }
    }

    else {
      subset = set - 1;

      // case direct computation

      if ((sup_bound - inf_bound) / success < B_THRESHOLD) {

        // computation of the upper bound probability

        term = 1.;
        for (i = 0;i < sup_bound - inf_bound;i++) {
          term *= success;
        }
        if (sup_bound < nb_value) {
          mass[sup_bound] = term;
        }

        // computation of the probabilities for the successive values (backward recurrence)

        ratio = failure / success;

        for (i = sup_bound - 1;i >= inf_bound;i--) {
          scale = (double)(subset + 1) / (double)(set - subset);
          subset--;
          term *= scale * ratio;
          if (i < nb_value) {
            mass[i] = term;
          }
        }
      }

      // case computation in log

      else {

        // computation of the upper bound probability

        term = (sup_bound - inf_bound) * log(success);
        if (sup_bound < nb_value) {
          mass[sup_bound] = exp(term);
        }

        // computation of the probabilities for the successive values (backward recurrence)

        ratio = log(failure / success);

        for (i = sup_bound - 1;i >= inf_bound;i--) {
          scale = (double)(subset + 1) / (double)(set - subset);
          subset--;
          term += log(scale) + ratio;
          if (i < nb_value) {
            mass[i] = exp(term);
          }
        }
      }
    }
  }

  cumul_computation();

# ifdef DEBUG2
  if (mode == STANDARD) {
    binomial dist(sup_bound - inf_bound , probability);


    cout << "\nTEST binomial distribution" << endl;
    for (i = inf_bound;i <= sup_bound;i++) {
      cout << i << "  " << pdf(dist , i - inf_bound) << " | " << mass[i]
           << "   " << cdf(dist , i - inf_bound) << " | " << cumul[i] << endl;
    }
  }
# endif

}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the probability mass function of a Poisson distribution.
 *         The number of values is determined using a threshold on the cumulative
 *         distribution function or using a predefined bound.
 *
 *  \param[in] inb_value       number of values,
 *  \param[in] cumul_threshold threshold on the cumulative distribution function,
 *  \param[in] mode            computation mode (STANDARD/RENEWAL).
 */
/*--------------------------------------------------------------*/

void DiscreteParametric::poisson_computation(int inb_value , double cumul_threshold ,
                                             distribution_computation mode)

{
  register int i;
  double log_parameter , num , denom , *pmass , *pcumul;


  pmass = mass;
  pcumul = cumul;

  switch (mode) {

  // case complete computation

  case STANDARD : {

    // null probability values before the lower bound of the support

    for (i = 0;i < inf_bound;i++) {
      *pmass++ = 0.;
      *pcumul++ = 0.;
    }

    i = 1;

    // case direct computation

    if (parameter < P_THRESHOLD) {

      // computation of the lower bound probability

      *pmass = exp(-parameter);
      *pcumul = *pmass;

      // computation of the probabilities for the successive values (forward recurrence)

      while (((*pcumul < cumul_threshold) || (i + inf_bound < inb_value)) &&
             (i + inf_bound < alloc_nb_value)) {
        pmass++;
        pcumul++;
        *pmass = *(pmass - 1) * parameter / i;
        *pcumul = *(pcumul - 1) + *pmass;
        i++;
      }
    }

    // case computation in log

    else {

      // computation of the lower bound probability

      num = -parameter;
      *pmass = exp(num);
      *pcumul = *pmass;

      // computation of the probabilities for the successive values (forward recurrence)

      log_parameter = log(parameter);
      denom = 0.;

      while (((*pcumul < cumul_threshold) || (i + inf_bound < inb_value)) &&
             (i + inf_bound < alloc_nb_value)) {
        num += log_parameter;
        denom += log((double)i);
        *++pmass = exp(num - denom);
        pcumul++;
        *pcumul = *(pcumul - 1) + *pmass;
        i++;
      }
    }

    i += inf_bound;
    break;
  }

  // case incomplete computation (renewal process)

  case RENEWAL : {

    // null probability values before the lower bound of the support

    for (i = 0;i < MIN(inb_value , inf_bound);i++) {
      *pmass++ = 0.;
      *pcumul++ = 0.;
    }

    if (inb_value > inf_bound) {
      i = 1;

      // case direct computation

      if (parameter < P_THRESHOLD) {

        // computation of the lower bound probability

        *pmass = exp(-parameter);
        *pcumul = *pmass;

        // computation of the probabilities for the successive values (forward recurrence)

        while ((*pcumul < cumul_threshold) && (i + inf_bound < inb_value)) {
          pmass++;
          pcumul++;
          *pmass = *(pmass - 1) * parameter / i;
          *pcumul = *(pcumul - 1) + *pmass;
          i++;
        }
      }

      // case computation in log

      else {

        // computation of the lower bound probability

        num = -parameter;
        *pmass = exp(num);
        *pcumul = *pmass;

        // computation of the probabilities for the successive values (forward recurrence)

        log_parameter = log(parameter);
        denom = 0.;

        while ((*pcumul < cumul_threshold) && (i + inf_bound < inb_value)) {
          num += log_parameter;
          denom += log((double)i);
          *++pmass = exp(num - denom);
          pcumul++;
          *pcumul = *(pcumul - 1) + *pmass;
          i++;
        }
      }

      i += inf_bound;
    }
    break;
  }
  }

  offset = MIN(inf_bound , i - 1);
  nb_value = i;

# ifdef DEBUG2
  if (mode == STANDARD) {
    poisson dist(parameter);


    cout << "\nTEST Poisson distribution" << endl;
    for (i = inf_bound;i < nb_value;i++) {
      cout << i << "  " << pdf(dist , i - inf_bound) << " | " << mass[i]
           << "   " << cdf(dist , i - inf_bound) << " | " << cumul[i] << endl;
    }
  }
# endif

}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the probability mass function of a negative binomial distribution.
 *         The number of values is determined using a threshold on the cumulative
 *         distribution function or using a predefined bound.
 *
 *  \param[in] inb_value       number of values,
 *  \param[in] cumul_threshold threshold on the cumulative distribution function,
 *  \param[in] mode            computation mode (STANDARD/RENEWAL).
 */
/*--------------------------------------------------------------*/

void DiscreteParametric::negative_binomial_computation(int inb_value , double cumul_threshold ,
                                                       distribution_computation mode)

{
  register int i;
  double failure = 1. - probability , success = probability , log_failure ,
         set , subset , scale , term , *pmass , *pcumul;


  pmass = mass;
  pcumul = cumul;

  switch (mode) {

    // case complete computation

    case STANDARD : {

    // null probability values before the lower bound of the support

    for (i = 0;i < inf_bound;i++) {
      *pmass++ = 0.;
      *pcumul++ = 0.;
    }

    i++;
    subset = parameter - 1.;
    set = subset;

    // case direct computation

    if (sqrt(parameter) / success < NB_THRESHOLD) {

      // computation of the lower bound probability

      term = pow(success , parameter);
      *pmass = term;
      *pcumul = *pmass;

      // computation of the probabilities for the successive values (forward recurrence)

      while (((*pcumul < cumul_threshold) || (i < inb_value)) &&
             (i < alloc_nb_value)) {
        set++;
        scale = set / (set - subset);
        term *= scale * failure;
        *++pmass = term;
        pcumul++;
        *pcumul = *(pcumul - 1) + *pmass;
        i++;
      }
    }

    // case computation in log

    else {

      // computation of the lower bound probability

      term = parameter * log(success);
      *pmass = exp(term);
      *pcumul = *pmass;

      // computation of the probabilities for the successive values (forward recurrence)

      log_failure = log(failure);

      while (((*pcumul < cumul_threshold) || (i < inb_value)) &&
             (i < alloc_nb_value)) {
        set++;
        scale = set / (set - subset);
        term += log(scale) + log_failure;
        *++pmass = exp(term);
        pcumul++;
        *pcumul = *(pcumul - 1) + *pmass;
        i++;
      }
    }
    break;
  }

  // case incomplete computation (renewal process)

  case RENEWAL : {

    // null probability values before the lower bound of the support

    for (i = 0;i < MIN(inf_bound , inb_value);i++) {
      *pmass++ = 0.;
      *pcumul++ = 0.;
    }

    if (inb_value > inf_bound) {
      i++;
      subset = parameter - 1.;
      set = subset;

      // case direct computation

      if (sqrt(parameter) / success < NB_THRESHOLD) {

        // computation of the lower bound probability

        term = pow(success , parameter);
        *pmass = term;
        *pcumul = *pmass;

        // computation of the probabilities for the successive values (forward recurrence)

        while ((*pcumul < cumul_threshold) && (i < inb_value)) {
          set++;
          scale = set / (set - subset);
          term *= scale * failure;
          *++pmass = term;
          pcumul++;
          *pcumul = *(pcumul - 1) + *pmass;
          i++;
        }
      }

      // case computation in log

      else {

        // computation of the lower bound probability

        term = parameter * log(success);
        *pmass = exp(term);
        *pcumul = *pmass;

        // computation of the probabilities for the successive values (forward recurrence)

        log_failure = log(failure);

        while ((*pcumul < cumul_threshold) && (i < inb_value)) {
          set++;
          scale = set / (set - subset);
          term += log(scale) + log_failure;
          *++pmass = exp(term);
          pcumul++;
          *pcumul = *(pcumul - 1) + *pmass;
          i++;
        }
      }
    }
    break;
  }
  }

  offset = MIN(inf_bound , i - 1);
  nb_value = i;

# ifdef DEBUG2
  if (mode == STANDARD) {
    negative_binomial dist(parameter , probability);


    cout << "TEST negative binomial distribution" << endl;
    for (i = inf_bound;i < nb_value;i++) {
      cout << i << "  " << pdf(dist , i - inf_bound) << " | " << mass[i]
           << "   " << cdf(dist , i - inf_bound) << " | " << cumul[i] << endl;
    }
  }
# endif

}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the probability mass function of a discrete uniform distribution.
 */
/*--------------------------------------------------------------*/

void DiscreteParametric::uniform_computation()

{
  register int i;
  double proba;


  offset = inf_bound;
  nb_value = sup_bound + 1;

  for (i = 0;i < inf_bound;i++) {
    mass[i] = 0.;
  }

  proba = 1. / (double)(sup_bound - inf_bound + 1);
  for (i = inf_bound;i <= sup_bound;i++) {
    mass[i] = proba;
  }

  cumul_computation();
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the number of values of a parametric discrete distribution
 *         (binomial, Poisson, negative binomial, uniform).
 *
 *  \param[in] ident           distribution identifier,
 *  \param[in] inf_bound       lower bound of the support,
 *  \param[in] sup_bound       upper bound of the support (binomial or uniform distribution),
 *  \param[in] parameter       parameter (negative binomial distribution),
 *  \param[in] probability     probability (binomial or negative binomial distribution),
 *  \param[in] cumul_threshold threshold on the cumulative distribution function.
 *
 *  \return number of values
 */
/*--------------------------------------------------------------*/

int DiscreteParametric::nb_value_computation(discrete_parametric ident , int inf_bound , int sup_bound ,
                                             double parameter , double probability , double cumul_threshold)

{
  int nb_value = 0;


  if ((ident == BINOMIAL) || (ident == UNIFORM)) {
    nb_value = sup_bound + 1;
  }

  else {
    if ((ident == POISSON) || (ident == NEGATIVE_BINOMIAL)) {
      DiscreteParametric *dist;

      dist = new DiscreteParametric(ident , inf_bound , sup_bound , parameter ,
                                    probability , cumul_threshold);
      nb_value = dist->nb_value;
      delete dist;
    }
  }

  return nb_value;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the probability mass function of a parametric discrete distribution
 *         (binomial, Poisson, negative binomial, uniform).
 *
 *  \param[in] min_nb_value    minimum number of values,
 *  \param[in] cumul_threshold threshold on the cumulative distribution function.
 */
/*--------------------------------------------------------------*/

void DiscreteParametric::computation(int min_nb_value , double cumul_threshold)

{
  if (ident > 0) {
    switch (ident) {
    case BINOMIAL :
      binomial_computation(1 , STANDARD);
      break;
    case POISSON :
      poisson_computation(min_nb_value , cumul_threshold , STANDARD);
      break;
    case NEGATIVE_BINOMIAL :
      negative_binomial_computation(min_nb_value , cumul_threshold , STANDARD);
      break;
    case UNIFORM :
      uniform_computation();
      break;
    }

    max_computation();
    mean_computation();
    variance_computation();
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the forward recurrence or sojourn time distribution on the basis of
 *         the recurrence or sojourn time distribution.
 *
 *  \param[in] dist recurrence or sojourn time distribution.
 */
/*--------------------------------------------------------------*/

void Forward::computation(const DiscreteParametric &dist)

{
  register int i;
  double norm;


  offset = 1;
  nb_value = dist.nb_value;
  mass[0] = 0.;

  // computation of the normalization quantity

  if (ident == CATEGORICAL) {
    norm = dist.mean;
  }
  else {
    norm = parametric_mean_computation();
  }

  // computation of the probability mass function

  for (i = 1;i < nb_value;i++) {
    mass[i] = (1. - dist.cumul[i - 1]) / norm;
  }

  cumul_computation();

  max_computation();
  mean_computation();
  variance_computation();
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the log-likelihood of the survival function of a given distribution.
 *
 *  \param[in] histo reference on a FrequencyDistribution object.
 *
 *  \return          log-likelihood.
 */
/*--------------------------------------------------------------*/

double Distribution::survivor_likelihood_computation(const FrequencyDistribution &histo) const

{
  register int i;
  double likelihood = 0.;


  if (histo.nb_element > 0) {
    if ((histo.offset == 0) || (histo.nb_value > nb_value)) {
      likelihood = D_INF;
    }

    else {
      for (i = histo.offset;i < histo.nb_value;i++) {
        if (histo.frequency[i] > 0) {
          if (cumul[i - 1] < 1.) {
            likelihood += histo.frequency[i] * log(1. - cumul[i - 1]);
          }
          else {
            likelihood = D_INF;
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
 *  \brief Computation of the chi2 value for a discrete distribution fit (Chi2 goodness of fit test).
 *
 *  \param[in] histo reference on a FrequencyDistribution object.
 *
 *  \return          chi2 value.
 */
/*--------------------------------------------------------------*/

double Distribution::chi2_value_computation(const FrequencyDistribution &histo) const

{
  register int i;
  double value , var1 , var2;


  if ((histo.offset < offset) || (histo.nb_value > nb_value)) {
    value = -D_INF;
  }

  else {
    value = 0.;

    for (i = offset;i < histo.nb_value;i++) {
      if (mass[i] > 0.) {
        var1 = histo.nb_element * mass[i];
        if (cumul[nb_value - 1] < CUMUL_THRESHOLD) {
          var1 /= cumul[nb_value - 1];
        }

        var2 = histo.frequency[i] - var1;
        value += var2 * var2 / var1;
      }

      else {
        if (histo.frequency[i] > 0) {
          value = -D_INF;
          break;
        }
      }
    }

    if (value != -D_INF) {
      var1 = 0.;
      for (i = histo.nb_value;i < nb_value;i++) {
        var1 += histo.nb_element * mass[i];
      }
      if (cumul[nb_value - 1] < CUMUL_THRESHOLD) {
        var1 /= cumul[nb_value - 1];
      }
      value += var1;
    }
  }

  return value;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Grouping of successive values, computation of the degrees of freedom and
 *         the chi2 value (Chi2 goodness of fit test).
 *
 *  \param[in] histo reference on a FrequencyDistribution object,
 *  \param[in] test  reference on a Test object.
 */
/*--------------------------------------------------------------*/

void Distribution::chi2_degree_of_freedom(const FrequencyDistribution &histo , Test &test) const

{
  register int i , j;
  int *filter_frequency;
  double *filter_mass;
  Distribution *filter_dist;
  FrequencyDistribution *filter_histo;


  if ((histo.offset >= offset) && (histo.nb_value <= nb_value)) {

    // construction and initialization of a Distribution and a FrequencyDistribution object

    filter_dist = new Distribution(nb_value);
    filter_dist->offset = offset;
    filter_dist->nb_parameter = nb_parameter;

    filter_histo = new FrequencyDistribution(histo.nb_value);
    filter_histo->nb_element = histo.nb_element;

    // grouping of values

    filter_mass = filter_dist->mass + offset;
    filter_frequency = filter_histo->frequency + offset;
    *filter_mass = 0.;
    *filter_frequency = 0;
    j = offset + 1;

    for (i = offset;i < histo.nb_value - 1;i++) {
      *filter_mass += mass[i];
      *filter_frequency += histo.frequency[i];
      if ((*filter_mass * histo.nb_element / cumul[nb_value - 1] > CHI2_FREQUENCY) &&
          ((1. - cumul[i]) * histo.nb_element / cumul[nb_value - 1] > 1)) {
        *++filter_mass = 0.;
        *++filter_frequency = 0;
        j++;
      }
    }
    *filter_frequency += histo.frequency[histo.nb_value - 1];

    filter_histo->offset_computation();
    filter_histo->nb_value = j;

    for (i = histo.nb_value - 1;i < nb_value - 1;i++) {
      *filter_mass += mass[i];
      if ((*filter_mass * histo.nb_element / cumul[nb_value - 1] > CHI2_FREQUENCY) &&
          ((1. - cumul[i]) * histo.nb_element / cumul[nb_value - 1] > 1)) {
        *++filter_mass = 0.;
        j++;
      }
    }
    *filter_mass += mass[nb_value - 1];

    filter_dist->nb_value = j;
    filter_dist->cumul_computation();

    // update of the degrees of freedom and computation of the chi2 value

    test.df1 = filter_dist->nb_value - filter_dist->offset -
               filter_dist->nb_parameter - 1;
    if (test.df1 < 1) {
      test.df1 = 1;
    }

    test.value = filter_dist->chi2_value_computation(*filter_histo);

#   ifdef DEBUG
/*    cout << *filter_histo;
    cout << *filter_dist; */
#   endif

    delete filter_dist;
    delete filter_histo;
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Chi2 goodness of fit test for a discrete distribution.
 *
 *  \param[in] histo reference on a FrequencyDistribution object,
 *  \param[in] test  reference on a Test object.
 */
/*--------------------------------------------------------------*/

void Distribution::chi2_fit(const FrequencyDistribution &histo , Test &test) const

{
  if ((histo.offset >= offset) && (histo.nb_value <= nb_value)) {
    chi2_degree_of_freedom(histo , test);
    if ((test.df1 > 0) && (test.value > 0.)) {
      test.chi2_critical_probability_computation();
    }
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Truncation of a distribution.
 *
 *  \param[in] error      reference on a StatError object,
 *  \param[in] imax_value maximum value.
 *
 *  \return               DiscreteParametricModel object.
 */
/*--------------------------------------------------------------*/

DiscreteParametricModel* Distribution::truncate(StatError &error , int imax_value) const

{
  register int i;
  DiscreteParametricModel *dist;


  error.init();

  if (imax_value <= offset) {
    dist = NULL;
    error.update(STAT_error[STATR_MAX_VALUE]);
  }

  else {

    // construction of a DiscreteParametricModel object

    dist = new DiscreteParametricModel(MIN(imax_value + 1 , nb_value));

    for (i = 0;i < dist->nb_value - 1;i++) {
      dist->mass[i] = mass[i];
      dist->cumul[i] = cumul[i];
    }
    dist->mass[dist->nb_value - 1] = 1. - dist->cumul[dist->nb_value - 2];
    dist->cumul[dist->nb_value - 1] = 1.;

    dist->offset = offset;
    dist->max_computation();
    dist->mean_computation();
    dist->variance_computation();
  }

  return dist;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Fit of a distribution.
 *
 *  \param[in] error reference on a StatError object,
 *  \param[in] idist reference on a DiscreteParametric object.
 *
 *  \return          DiscreteParametricModel object.
 */
/*--------------------------------------------------------------*/

DiscreteParametricModel* FrequencyDistribution::fit(StatError &error ,
                                                    const DiscreteParametric &idist) const

{
  DiscreteParametricModel *dist;


  error.init();

  if ((offset < idist.offset) || (nb_value > idist.nb_value)) {
    dist = NULL;
    error.update(STAT_error[STATR_VALUE_RANGE]);
  }

  else {
    dist = new DiscreteParametricModel(idist , this);
  }

  return dist;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Estimation of a discrete parametric distribution (binomial, Poisson or negative binomial).
 *
 *  \param[in] ident           distribution identifier,
 *  \param[in] min_inf_bound   minimum lower bound,
 *  \param[in] flag            flag on the estimation of the lower bound,
 *  \param[in] cumul_threshold threshold on the cumulative distribution function.
 *
 *  \return                    DiscreteParametric object.
 */
/*--------------------------------------------------------------*/

DiscreteParametric* FrequencyDistribution::parametric_estimation(discrete_parametric ident , int min_inf_bound ,
                                                                 bool flag , double cumul_threshold) const

{
  double likelihood;
  DiscreteParametric *dist;


  // construction of a DiscreteParametric object

  dist = new DiscreteParametric((int)(nb_value * SAMPLE_NB_VALUE_COEFF) , ident);

  // parameter estimation

  likelihood = Reestimation<int>::parametric_estimation(dist , min_inf_bound ,
                                                        flag , cumul_threshold);

  // update of the estimated distribution

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
 *  \brief Estimation of a discrete parametric distribution (binomial, Poisson or negative binomial).
 *
 *  \param[in] error           reference on a StatError object,
 *  \param[in] ident           distribution identifier,
 *  \param[in] min_inf_bound   minimum lower bound,
 *  \param[in] flag            flag on the estimation of the lower bound,
 *  \param[in] cumul_threshold threshold on the cumulative distribution function.
 *
 *  \return                    DiscreteParametricModel object.
 */
/*--------------------------------------------------------------*/

DiscreteParametricModel* FrequencyDistribution::parametric_estimation(StatError &error , discrete_parametric ident ,
                                                                      int min_inf_bound , bool flag ,
                                                                      double cumul_threshold) const

{
  double likelihood;
  DiscreteParametricModel *dist;


  error.init();

  if ((min_inf_bound < 0) || (min_inf_bound > 1) || (min_inf_bound > offset)) {
    dist = NULL;
    error.update(STAT_error[STATR_MIN_INF_BOUND]);
  }

  else {

    // construction of a DiscreteParametricModel object

    dist = new DiscreteParametricModel((int)(nb_value * SAMPLE_NB_VALUE_COEFF) , ident);
    dist->frequency_distribution = new DiscreteDistributionData(*this);

    // parameter estimation

    likelihood = Reestimation<int>::parametric_estimation(dist , min_inf_bound ,
                                                          flag , cumul_threshold);

    if (likelihood != D_INF) {

      // update of the estimated distribution

      dist->computation(nb_value , cumul_threshold);

      // update of the number of free parameters

      dist->nb_parameter_update();
      if (!flag) {
        (dist->nb_parameter)--;
      }
    }

    else {
      delete dist;
      dist = NULL;
      error.update(STAT_error[STATR_ESTIMATION_FAILURE]);
    }
  }

  return dist;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Estimation of a discrete parametric distribution (binomial, Poisson or negative binomial).
 *
 *  \param[in] error           reference on a StatError object,
 *  \param[in] min_inf_bound   minimum lower bound,
 *  \param[in] flag            flag on the estimation of the lower bound,
 *  \param[in] cumul_threshold threshold on the cumulative distribution function.
 *
 *  \return                    DiscreteParametricModel object.
 */
/*--------------------------------------------------------------*/

DiscreteParametricModel* FrequencyDistribution::type_parametric_estimation(StatError &error ,
                                                                           int min_inf_bound , bool flag ,
                                                                           double cumul_threshold) const

{
  double likelihood;
  DiscreteParametricModel *dist;


  error.init();

  if ((min_inf_bound < 0) || (min_inf_bound > 1) || (min_inf_bound > offset)) {
    dist = NULL;
    error.update(STAT_error[STATR_MIN_INF_BOUND]);
  }

  else {

    // construction of a DiscreteParametricModel object

    dist = new DiscreteParametricModel((int)(nb_value * SAMPLE_NB_VALUE_COEFF));
    dist->frequency_distribution = new DiscreteDistributionData(*this);

    // parameter estimation

    likelihood = Reestimation<int>::type_parametric_estimation(dist , min_inf_bound ,
                                                               flag , cumul_threshold);

    if (likelihood != D_INF) {

      // update of the estimated distribution

      dist->computation(nb_value , cumul_threshold);

      // update of the number of free parameters

      dist->nb_parameter_update();
      if (!flag) {
        (dist->nb_parameter)--;
      }
    }

    else {
      delete dist;
      dist = NULL;
      error.update(STAT_error[STATR_ESTIMATION_FAILURE]);
    }
  }

  return dist;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the penalty terms in the framework of a penalized likelihood approach.
 *
 *  \param[in] weight   penalty weight,
 *  \param[in] pen_type penalty type (1st order, 2nd order difference or entropy),
 *  \param[in] penalty  penalties,
 *  \param[in] outside  management of side effects (zero outside the support or
 *                      continuation of the distribution).
 */
/*--------------------------------------------------------------*/

void Distribution::penalty_computation(double weight , penalty_type pen_type ,
                                       double *penalty , side_effect outside) const

{
  register int i;


  switch (pen_type) {

  case FIRST_DIFFERENCE : {
    switch (outside) {
    case ZERO :
      penalty[offset] = 2 * weight * (2 * mass[offset] - mass[offset + 1]);
      break;
    case CONTINUATION :
      penalty[offset] = 2 * weight * (mass[offset] - mass[offset + 1]);
      break;
    }

    for (i = offset + 1;i < nb_value - 1;i++) {
      penalty[i] = 2 * weight * (-mass[i - 1] + 2 * mass[i] - mass[i + 1]);
    }

    switch (outside) {
    case ZERO :
      penalty[nb_value - 1] = 2 * weight * (-mass[nb_value - 2] + 2 * mass[nb_value - 1]);
      break;
    case CONTINUATION :
      penalty[nb_value - 1] = 2 * weight * (-mass[nb_value - 2] + mass[nb_value - 1]);
      break;
    }
    break;
  }

  case SECOND_DIFFERENCE : {
    i = offset;

    switch (outside) {

    case ZERO : {
      penalty[offset] = 2 * weight * (6 * mass[i] - 4 * mass[i + 1] + mass[i + 2]);
      i++;
      penalty[offset + 1] = 2 * weight * (-4 * mass[i - 1] + 6 * mass[i] - 4 * mass[i + 1] +
                                          mass[i + 2]);
      break;
    }

    case CONTINUATION : {
      penalty[offset] = 2 * weight * (3 * mass[i] - 4 * mass[i + 1] + mass[i + 2]);
      i++;
      penalty[offset + 1] = 2 * weight * (-3 * mass[i - 1] + 6 * mass[i] - 4 * mass[i + 1] +
                                          mass[i + 2]);
      break;
    }
    }

    for (i = offset + 2;i < nb_value - 2;i++) {
      penalty[i] = 2 * weight * (mass[i - 2] - 4 * mass[i - 1] + 6 * mass[i] -
                                 4 * mass[i + 1] + mass[i + 2]);
    }

    i = nb_value - 2;

    switch (outside) {

    case ZERO : {
      penalty[nb_value - 2] = 2 * weight * (mass[i - 2] - 4 * mass[i - 1] + 6 * mass[i] -
                                            4 * mass[i + 1]);
      i++;
      penalty[nb_value - 1] = 2 * weight * (mass[i - 2] - 4 * mass[i - 1] + 6 * mass[i]);
      break;
    }

    case CONTINUATION : {
      penalty[nb_value - 2] = 2 * weight * (mass[i - 2] - 4 * mass[i - 1] + 6 * mass[i] -
                                            3 * mass[i + 1]);
      i++;
      penalty[nb_value - 1] = 2 * weight * (mass[i - 2] - 4 * mass[i - 1] + 3 * mass[i]);
      break;
    }
    }
    break;
  }

  case ENTROPY : {
    for (i = offset;i < nb_value;i++) {
      penalty[i] = weight * (log(mass[i]) + 1);
    }
    break;
  }
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Reestimation of the parameters of a discrete distribution
 *         (binomial, Poisson or negative binomial).
 *
 *  \param[in] reestim  reference on the reestimation quantities,
 *  \param[in] nb_estim number of reestimated parameters (negative binomial).
 */
/*--------------------------------------------------------------*/

void DiscreteParametric::reestimation(const Reestimation<double> *reestim , int nb_estim)

{
  switch (ident) {

  case BINOMIAL : {
    probability = (reestim->mean - inf_bound) / (sup_bound - inf_bound);
    break;
  }

  case POISSON : {
    parameter = reestim->mean - inf_bound;
    break;
  }

  case NEGATIVE_BINOMIAL : {
    switch (nb_estim) {

    case 1 : {
      if (reestim->mean - inf_bound + parameter > 0.) {
        probability = parameter / (reestim->mean - inf_bound + parameter);
      }
      break;
    }

    case 2 : {
/*      register int i;
      double previous_parameter = parameter , sum1 , sum2; */

      parameter = (reestim->mean - inf_bound) * probability / (1. - probability);

/*     sum1 = 0.;
      sum2 = 0.;
      for (i = inf_bound + 1;i < nb_value;i++) {
        sum2 += 1. / (i - inf_bound + previous_parameter - 1);
        sum1 += reestim->frequency[i] * sum2;
      }

      probability = exp(-sum1 / reestim->nb_element); */

#     ifdef DEBUG
//      cout << "<" << probability << "> ";
#     endif
      break;
    }
    }

    break;
  }
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Simulation using the cumulative distribution function.
 *         The median is used to determine the direction of the algorithm
 *         (forward from the lower bound or backward from the highest value).
 *
 *  \param[in] nb_value number of values,
 *  \param[in] cumul    pointer on the cumulative distribution function,
 *  \param[in] scale    scaling factor.
 *
 *  \return             generated value.
 */
/*--------------------------------------------------------------*/

int cumul_method(int nb_value , const double *cumul , double scale)

{
  register int i;
  double limit;


  limit = ((double)rand() / (RAND_MAX + 1.)) * scale;
//  limit = ((double)random() / (double)0x7fffffff) * scale;

  if ((limit < cumul[nb_value / 2])) {
    i = 0;
    while (*cumul++ <= limit) {
      i++;
    }
  }

  else {
    i = nb_value - 1;
    cumul += nb_value - 1;
    while (*--cumul > limit) {
      i--;
    }
  }

  return i;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Simulation using the cumulative distribution function.
 *
 *  \return generated value.
 */
/*--------------------------------------------------------------*/

int Distribution::simulation() const

{
  int value;


  value = offset + cumul_method(nb_value - offset , cumul + offset , 1. - complement);

  return value;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Simulation using the rejection sampling method.
 *         Principle: A point (x, Px) is drawn in the rectangle [xmin, xmax] x [0. ,Pmax].
 *         If the point is below the distribution, the x value is kept;
 *         If not, a new point is drawn.
 *
 *  \return generated value.
 */
/*--------------------------------------------------------------*/

int DiscreteParametric::simulation() const

{
  int range = nb_value - offset , value;
  double x , y;


  if ((ident == CATEGORICAL) || (range < MIN_RANGE) || (range * max > MAX_SURFACE)) {
    value = Distribution::simulation();
  }

  else {
    do {
      x = (double)rand() / (RAND_MAX + 1.);
      y = (double)rand() / (RAND_MAX + 1.);
//      x = (double)random() / (double)0x7fffffff;
//      y = (double)random() / (double)0x7fffffff;
      value = (int)(offset + range * x);
    }
    while (y * max > mass[value]);
  }

  return value;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Building of a frequency distribution by simulating a discrete parametric distribution.
 *
 *  \param[in] error      reference on a StatError object,
 *  \param[in] nb_element sample size.
 *
 *  \return               DiscreteDistributionData object.
 */
/*--------------------------------------------------------------*/

DiscreteDistributionData* DiscreteParametricModel::simulation(StatError &error ,
                                                              int nb_element) const

{
  register int i;
  DiscreteDistributionData *histo;


  error.init();

  if ((nb_element < 1) || (nb_element > DIST_NB_ELEMENT)) {
    histo = NULL;
    error.update(STAT_error[STATR_SAMPLE_SIZE]);
  }

  else {

    // construction of the frequency distribution

    histo = new DiscreteDistributionData(*this);
    histo->distribution = new DiscreteParametricModel(*this , false);

    // simulation

    for (i = 0;i < nb_element;i++) {
      (histo->frequency[DiscreteParametric::simulation()])++;
    }

    // extraction of the frequency distribution characteristics

    histo->nb_value_computation();
    histo->offset_computation();
    histo->nb_element = nb_element;
    histo->max_computation();
    histo->mean_computation();
    histo->variance_computation();
  }

  return histo;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the log-likelihood of a discrete-time renewal process for time interval data.
 *
 *  \param[in] forward_dist forward recurrence time distribution,
 *  \param[in] within       complete time interval frequency distribution,
 *  \param[in] backward     backward recurrence time frequency distribution,
 *  \param[in] forward      forward recurrence time frequency distribution,
 *  \param[in] no_event     observation period frequency distribution for the case of no event.
 *
 *  \return                 log-likelihood.
 */
/*--------------------------------------------------------------*/

double DiscreteParametric::renewal_likelihood_computation(const Forward &forward_dist ,
                                                          const FrequencyDistribution &within ,
                                                          const FrequencyDistribution &backward ,
                                                          const FrequencyDistribution &forward ,
                                                          const FrequencyDistribution *no_event) const

{
  double likelihood , buff;
  FrequencyDistribution *histo;


  likelihood = likelihood_computation(within);

  if (likelihood != D_INF) {
    histo = new FrequencyDistribution(backward , SHIFT , 1);
    buff = survivor_likelihood_computation(*histo);
    delete histo;

    if (buff != D_INF) {
      likelihood += buff;
      buff = forward_dist.likelihood_computation(forward);

      if (buff != D_INF) {
        likelihood += buff;

        if (no_event) {
          histo = new FrequencyDistribution(*no_event , SHIFT , 1);
          buff = forward_dist.survivor_likelihood_computation(*histo);
          delete histo;

          if (buff != D_INF) {
            likelihood += buff;
          }
          else {
            likelihood = D_INF;
          }
        }
      }

      else {
        likelihood = D_INF;
      }
    }

    else {
      likelihood = D_INF;
    }
  }

  return likelihood;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the reestimation quantities corresponding to the inter-event distribution
 *         (EM algorithm of an equilibrium renewal process estimated on the basis of time interval data).
 *
 *  \param[in] within              complete time interval frequency distribution,
 *  \param[in] backward            backward recurrence time frequency distribution,
 *  \param[in] forward             forward recurrence time frequency distribution,
 *  \param[in] no_event            observation period frequency distribution for the case of no event,
 *  \param[in] inter_event_reestim pointer on the reestimation quantities of the inter-event distribution
 *  \param[in] length_bias_reestim pointer on the reestimation quantities of the length-biased distribution,
 *  \param[in] iter                EM iteration.
 */
/*--------------------------------------------------------------*/

void DiscreteParametric::expectation_step(const FrequencyDistribution &within ,
                                          const FrequencyDistribution &backward ,
                                          const FrequencyDistribution &forward ,
                                          const FrequencyDistribution *no_event ,
                                          Reestimation<double> *inter_event_reestim ,
                                          Reestimation<double> *length_bias_reestim , int iter) const

{
  register int i , j;
  int reestim_offset , reestim_nb_value , *pfrequency;
  double sum , *ifrequency , *lfrequency , *pmass , *pcumul , *norm;


  // initializations

  ifrequency = inter_event_reestim->frequency;
  lfrequency = length_bias_reestim->frequency;
  for (i = 0;i < inter_event_reestim->alloc_nb_value;i++) {
    *ifrequency++ = 0.;
    *lfrequency++ = 0.;
  }

  // computation of the reestimation quantities of the inter-event distribution

  ifrequency = inter_event_reestim->frequency + within.offset;
  pfrequency = within.frequency + within.offset;
  for (i = within.offset;i < within.nb_value;i++) {
    *ifrequency++ += *pfrequency++;
  }

  pfrequency = backward.frequency + backward.offset;
  pcumul = cumul + backward.offset;
  ifrequency = inter_event_reestim->frequency + backward.offset + 1;
  pmass = mass + backward.offset + 1;
  sum = 0.;

  for (i = backward.offset;i < backward.nb_value;i++) {
    sum += *pfrequency++ / (1. - *pcumul++);
    *ifrequency++ += *pmass++ * sum;
  }
  for (i = backward.nb_value;i < nb_value - 1;i++) {
    *ifrequency++ += *pmass++ * sum;
  }

  // computation of the reestimation quantities of the length-biased distribution

  pfrequency = forward.frequency + forward.offset;
  pcumul = cumul + forward.offset - 1;
  lfrequency = length_bias_reestim->frequency + forward.offset;
  pmass = mass + forward.offset;
  sum = 0.;

  for (i = forward.offset;i < forward.nb_value;i++) {
    sum += *pfrequency++ / (1. - *pcumul++);
    *lfrequency++ += *pmass++ * sum;
  }
  for (i = forward.nb_value;i < nb_value;i++) {
    *lfrequency++ += *pmass++ * sum;
  }

  if (no_event) {
    norm = new double[no_event->nb_value];

    for (i = no_event->offset;i < no_event->nb_value;i++) {
      if (no_event->frequency[i] > 0) {
        pmass = mass + i + 1;
        norm[i] = 0.;
        for (j = i + 1;j < nb_value;j++) {
          norm[i] += (j - i) * *pmass++;
        }
      }
    }

    lfrequency = length_bias_reestim->frequency + no_event->offset + 1;
    pmass = mass + no_event->offset + 1;
    for (i = no_event->offset + 1;i < nb_value;i++) {
      pfrequency = no_event->frequency + no_event->offset;
      sum = 0.;
      for (j = no_event->offset;j < MIN(i , no_event->nb_value);j++) {
        if ((*pfrequency > 0) && (norm[j] > 0.)) {
          sum += *pfrequency * (i - j) / norm[j];
        }
        pfrequency++;
      }

      *lfrequency++ += *pmass++ * sum;
    }

    delete [] norm;
  }

  reestim_offset = 1;
  reestim_nb_value = inter_event_reestim->alloc_nb_value;

  ifrequency = inter_event_reestim->frequency + inter_event_reestim->alloc_nb_value;
  lfrequency = length_bias_reestim->frequency + inter_event_reestim->alloc_nb_value;
  while ((*--ifrequency == 0) && (*--lfrequency == 0) && (reestim_nb_value > 2)) {
    reestim_nb_value--;
  }
  inter_event_reestim->nb_value = reestim_nb_value;
  length_bias_reestim->nb_value = reestim_nb_value;

  ifrequency = inter_event_reestim->frequency + reestim_offset;
  lfrequency = length_bias_reestim->frequency + reestim_offset;
  while ((*ifrequency++ == 0) && (*lfrequency++ == 0) && (reestim_offset < reestim_nb_value - 1)) {
    reestim_offset++;
  }
  inter_event_reestim->offset = reestim_offset;
  length_bias_reestim->offset = reestim_offset;

  inter_event_reestim->nb_element_computation();
  length_bias_reestim->nb_element_computation();

# ifdef DEBUG
  if ((iter < 10) || ((iter < 100) && (iter % 10 == 0)) ||
      ((iter < 1000) && (iter % 100 == 0)) || (iter % 1000 == 0)) {
    inter_event_reestim->max_computation();
    inter_event_reestim->mean_computation();
    inter_event_reestim->variance_computation();

    length_bias_reestim->max_computation();
    length_bias_reestim->mean_computation();
    length_bias_reestim->variance_computation();

    cout << "\nreestimation quantities of the inter-event distribution:" << *inter_event_reestim << endl;
    cout << "\nreestimation quantities of the length-biased distribution:" << *length_bias_reestim << endl;
  }
# endif

}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the mean of a distribution by interval bissection.
 *
 *  \param[in] distribution_reestim pointer on the reestimation quantities of the distribution
 *  \param[in] length_bias_reestim  pointer on the reestimation quantities of the length-biased distribution.
 *
 *  \return                         distribution mean.
 */
/*--------------------------------------------------------------*/

double interval_bisection(Reestimation<double> *distribution_reestim ,
                          Reestimation<double> *length_bias_reestim)

{
  register int i;
  double ratio , inf_ratio , sup_ratio , mean , inf_mean , sup_mean , *dfrequency , *lfrequency;

# ifdef DEBUG
  int iter = 0;
# endif


  // initializations: computation of the first two values

  dfrequency = distribution_reestim->frequency + distribution_reestim->offset;
  lfrequency = length_bias_reestim->frequency + distribution_reestim->offset;
  inf_ratio = 0.;
  sup_ratio = 0.;
  inf_mean = distribution_reestim->offset;
  sup_mean = distribution_reestim->nb_value - 1;

  for (i = distribution_reestim->offset;i < distribution_reestim->nb_value;i++) {
    inf_ratio += (*dfrequency + *lfrequency) * i /
                 (distribution_reestim->nb_element * sup_mean + length_bias_reestim->nb_element * i);
    sup_ratio += (*dfrequency++ + *lfrequency++) * i /
                 (distribution_reestim->nb_element * inf_mean + length_bias_reestim->nb_element * i);
  }

  do {
    dfrequency = distribution_reestim->frequency + distribution_reestim->offset;
    lfrequency = length_bias_reestim->frequency + distribution_reestim->offset;
    ratio = 0.;
    mean = (inf_mean + sup_mean) / 2.;

    for (i = distribution_reestim->offset;i < distribution_reestim->nb_value;i++) {
      ratio += (*dfrequency++ + *lfrequency++) * i /
               (distribution_reestim->nb_element * mean + length_bias_reestim->nb_element * i);
    }

#   ifdef DEBUG
    cout << STAT_label[STATL_ITERATION] << " " << iter++ << ": " << mean << " " << ratio << endl;
#   endif

    if (ratio < 1.) {
      inf_ratio = ratio;
      sup_mean = mean;
    }
    else {
      sup_ratio = ratio;
      inf_mean = mean;
    }
  }
  while (sup_ratio - inf_ratio > BISECTION_RATIO_THRESHOLD);

  mean = (inf_mean + sup_mean) / 2.;

# ifdef DEBUG
  cout << STAT_label[STATL_MEAN] << ": " << mean << " " << ratio << endl;
# endif

  return mean;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Estimation of an equilibrium renewal process on the basis of time interval data
 *         using the EM algorithm.
 *
 *  \param[in] error             reference on a StatError object,
 *  \param[in] os                stream,
 *  \param[in] backward          backward recurrence time frequency distribution,
 *  \param[in] forward           forward recurrence time frequency distribution,
 *  \param[in] no_event          observation period frequency distribution for the case of no event,
 *  \param[in] iinter_event      reference on the initial inter-event distribution,
 *  \param[in] estimator         estimator type (likelihood or penalized likelihood),
 *  \param[in] nb_iter           number of iterations,
 *  \param[in] mean_estimator    method for the computation of the mean of the inter-event distribution,
 *  \param[in] weight            penalty weight,
 *  \param[in] pen_type          penalty type,
 *  \param[in] outside           management of side effects (zero outside the support or
 *                               continuation of the distribution),
 *  \param[in] iinter_event_mean mean of the inter-event distribution.
 *
 *  \return                      DiscreteParametricModel object.
 */
/*--------------------------------------------------------------*/

DiscreteParametricModel* FrequencyDistribution::estimation(StatError &error , ostream &os ,
                                                           const FrequencyDistribution &backward ,
                                                           const FrequencyDistribution &forward ,
                                                           const FrequencyDistribution *no_event ,
                                                           const DiscreteParametric &iinter_event ,
                                                           estimation_criterion estimator , int nb_iter ,
                                                           duration_distribution_mean_estimator mean_estimator ,
                                                           double weight , penalty_type pen_type , side_effect outside ,
                                                           double iinter_event_mean) const

{
  bool status = true;
  register int i;
  int inb_value , max_nb_value;
  double likelihood , previous_likelihood , inter_event_mean , *penalty;
  DiscreteParametricModel *inter_event;
  Forward *forward_dist;
  Reestimation<double> *inter_event_reestim , *length_bias_reestim;
  FrequencyDistribution *backward_forward;
  const FrequencyDistribution *phisto[2];


  inter_event = NULL;
  error.init();

  if (nb_element < NB_COMPLETE_INTERVAL) {
    status = false;
    error.update(STAT_error[STATR_NB_COMPLETE_INTERVAL_TOO_SMALL]);
  }

  if (offset == 0) {
    status = false;
    error.update(STAT_error[STATR_COMPLETE_MIN_VALUE]);
  }
  if (forward.offset == 0) {
    status = false;
    error.update(STAT_error[STATR_FORWARD_MIN_VALUE]);
  }
  if ((no_event) && (no_event->offset == 0)) {
    status = false;
    error.update(STAT_error[STATR_NO_EVENT_MIN_VALUE]);
  }

  if ((nb_iter != I_DEFAULT) && (nb_iter < 1)) {
    status = false;
    error.update(STAT_error[STATR_NB_ITERATION]);
  }

  if ((weight != D_DEFAULT) && (weight <= 0.)) {
    status = false;
    error.update(STAT_error[STATR_PENALTY_WEIGHT]);
  }

  if ((mean_estimator == ESTIMATED) && (iinter_event_mean == D_DEFAULT)) {
    status = false;
    error.update(STAT_error[STATR_MEAN_ESTIMATION]);
  }

  inb_value = nb_value;
  if (backward.nb_value + 1 > inb_value) {
    inb_value = backward.nb_value + 1;
  }
  if (forward.nb_value > inb_value) {
    inb_value = forward.nb_value;
  }

  if ((no_event) && (no_event->nb_value > inb_value)) {
    max_nb_value = no_event->nb_value;
  }
  else {
    max_nb_value = inb_value;
  }

  if ((iinter_event.offset > offset) || (iinter_event.nb_value < max_nb_value)) {
    status = false;
    error.update(STAT_error[STATR_INTER_EVENT_SUPPORT]);
  }

  if (status) {
    phisto[0] = new FrequencyDistribution(backward , SHIFT , 1);
    phisto[1] = &forward;
    backward_forward = new FrequencyDistribution(2 , phisto);
    delete phisto[0];

#   ifdef MESSAGE
    {
      int max_nb_element , width[2];
      long old_adjust;


      old_adjust = os.setf(ios::right , ios::adjustfield);

      width[0] = column_width(max_nb_value - 1);

      max_nb_element = nb_element;
      if (backward_forward->nb_element > max_nb_element) {
        max_nb_element = backward_forward->nb_element;
      }
      if ((no_event) && (no_event->nb_element > max_nb_element)) {
        max_nb_element = no_event->nb_element;
      }
      width[1] = column_width(max_nb_element) + ASCII_SPACE;

      os << "\n   | " << STAT_label[STATL_OBSERVATION_INTER_EVENT] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
         << " | " << STAT_label[STATL_BACKWARD] << "/" << STAT_label[STATL_FORWARD]
         << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
      if (no_event) {
        os << " | no-event " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
      }
      os << endl;

      for (i = 0;i < max_nb_value;i++) {
        os << setw(width[0]) << i;

        if (i < nb_value) {
          os << setw(width[1]) << frequency[i];
        }
        else {
          os << setw(width[1]) << " ";
        }

        if (i < backward_forward->nb_value) {
          os << setw(width[1]) << backward_forward->frequency[i];
        }
        else {
          os << setw(width[1]) << " ";
        }

        if (no_event) {
          if (i < no_event->nb_value) {
            os << setw(width[1]) << no_event->frequency[i];
          }
          else {
            os << setw(width[1]) << " ";
          }
        }

        os << "    |  ";
        if (i < backward.nb_value) {
          os << setw(width[1]) << backward.frequency[i];
        }
        else {
          os << setw(width[1]) << " ";
        }

        if (i < forward.nb_value) {
          os << setw(width[1]) << forward.frequency[i];
        }
        else {
          os << setw(width[1]) << " ";
        }

        os << endl;
      }
      os << endl;

      os << setw(width[0]) << " "
         << setw(width[1]) << nb_element
         << setw(width[1]) << backward_forward->nb_element;
      if (no_event) {
        os << setw(width[1]) << no_event->nb_element;
      }
      os << "    |  " << setw(width[1]) << backward.nb_element
         << setw(width[1]) << forward.nb_element << "\n" << endl;

      os.setf((FMTFLAGS)old_adjust , ios::adjustfield);
    }
#   endif

    // construction of the inter-event distribution

    inter_event = new DiscreteParametricModel(iinter_event , this);
    inter_event->init(CATEGORICAL , I_DEFAULT , I_DEFAULT , D_DEFAULT , D_DEFAULT);
    forward_dist = new Forward(*inter_event);

    if (estimator == PENALIZED_LIKELIHOOD) {
      penalty = new double[inter_event->nb_value];

      if (weight == D_DEFAULT) {
        if (pen_type != ENTROPY) {
          weight = RENEWAL_DIFFERENCE_WEIGHT;
        }
        else {
          weight = RENEWAL_ENTROPY_WEIGHT;
        }
      }

      if (no_event) {
        weight *= (nb_element + backward.nb_element + forward.nb_element + no_event->nb_element);
      }
      else {
        weight *= (nb_element + backward.nb_element + forward.nb_element);
      }
    }

    inter_event_reestim = new Reestimation<double>(inter_event->nb_value);
    length_bias_reestim = new Reestimation<double>(inter_event->nb_value);

    likelihood = D_INF;
    i = 0;

    do {
      i++;

      inter_event->expectation_step(*this , backward , forward , no_event ,
                                    inter_event_reestim , length_bias_reestim , i);

      switch (estimator) {

      case LIKELIHOOD : {
        switch (mean_estimator) {
        case ESTIMATED :
          inter_event_mean = iinter_event_mean;
          break;
        case COMPUTED :
          inter_event_mean = interval_bisection(inter_event_reestim , length_bias_reestim);
          break;
        case ONE_STEP_LATE :
          inter_event_mean = inter_event->mean;
          break;
        }

        inter_event_reestim->equilibrium_process_estimation(length_bias_reestim , inter_event ,
                                                            inter_event_mean);
        break;
      }

      case PENALIZED_LIKELIHOOD : {
        switch (mean_estimator) {
        case ESTIMATED :
          inter_event_mean = iinter_event_mean;
          break;
        case ONE_STEP_LATE :
          inter_event_mean = inter_event->mean;
          break;
        }

        inter_event_reestim->penalized_likelihood_equilibrium_process_estimation(length_bias_reestim ,
                                                                                 inter_event , inter_event_mean ,
                                                                                 weight , pen_type , penalty ,
                                                                                 outside);
        break;
      }
      }

      forward_dist->computation(*inter_event);
      previous_likelihood = likelihood;
      likelihood = inter_event->renewal_likelihood_computation(*forward_dist , *this , backward ,
                                                               forward , no_event);

#     ifdef MESSAGE
      if ((i < 10) || ((i < 100) && (i % 10 == 0)) || ((i < 1000) && (i % 100 == 0)) || (i % 1000 == 0)) {
        os << STAT_label[STATL_ITERATION] << " " << i << "   "
           << STAT_label[STATL_LIKELIHOOD] << ": " << likelihood << "   "
           << STAT_label[STATL_SMOOTHNESS] << ": " << inter_event->second_difference_norm_computation();
        if (estimator == PENALIZED_LIKELIHOOD) {
          os << "   cumul: " << inter_event->cumul[inter_event->nb_value - 1];
        }

        if ((no_event) && (no_event->offset + 1 == no_event->nb_value) && (backward_forward->nb_value > nb_value) &&
            ((forward.nb_element + no_event->nb_element) * (1. - inter_event->cumul[inb_value - 2]) > 0.)) {
          if (mean_estimator == ESTIMATED) {
            inter_event_mean = iinter_event_mean;
          }
          else {
            inter_event_mean = inter_event->mean;
          }

          os << "   smaller upper bound: "
             << inb_value - 1 + (no_event->nb_element * inter_event_mean) /
                                ((forward.nb_element + no_event->nb_element) * (1. - inter_event->cumul[inb_value - 2]));
        }

/*        if (backward_forward->nb_value > nb_value) {
          register int j;
          double term = forward.nb_element * (backward_forward->nb_value - 1) /
                        inter_event->mean + nb_element + backward.nb_element;
          if (no_event) {
            term += no_event->nb_element * (backward_forward->nb_value - 1) / inter_event->mean;
          }
          for (j = backward_forward->offset;j < backward_forward->nb_value;j++) {
            term -= backward_forward->frequency[j] / (1. - inter_event->cumul[j - 1]);
          }

          os << " |   " << term;
        } */

        os << endl;
      }
#     endif

    }
    while ((likelihood != D_INF) && (((nb_iter == I_DEFAULT) && (i < RENEWAL_NB_ITER) &&
             ((likelihood - previous_likelihood) / -likelihood > RENEWAL_LIKELIHOOD_DIFF)) ||
            ((nb_iter != I_DEFAULT) && (i < nb_iter))));

    if (likelihood != D_INF) {

#     ifdef MESSAGE
      os << "\n" << i << " " << STAT_label[STATL_ITERATIONS] << "   "
         << STAT_label[STATL_LIKELIHOOD] << ": " << likelihood << "   "
         << STAT_label[STATL_SMOOTHNESS] << ": " << inter_event->second_difference_norm_computation();
      if (estimator == PENALIZED_LIKELIHOOD) {
        os << "   cumul: " << inter_event->cumul[inter_event->nb_value - 1];
      }

      if ((no_event) && (no_event->offset + 1 == no_event->nb_value) && (backward_forward->nb_value > nb_value) &&
          ((forward.nb_element + no_event->nb_element) * (1. - inter_event->cumul[inb_value - 2]) > 0.)) {
        if (mean_estimator == ESTIMATED) {
          inter_event_mean = iinter_event_mean;
        }
        else {
          inter_event_mean = inter_event->mean;
        }

        os << "   smaller upper bound: "
           << inb_value - 1 + (no_event->nb_element * inter_event_mean) /
                              ((forward.nb_element + no_event->nb_element) * (1. - inter_event->cumul[inb_value - 2]));
      }
      os << endl;
#     endif

    }

    else {
      delete inter_event;
      inter_event = NULL;
      error.update(STAT_error[STATR_ESTIMATION_FAILURE]);
    }

    delete backward_forward;

    delete forward_dist;
    if (estimator == PENALIZED_LIKELIHOOD) {
      delete [] penalty;
    }

    delete inter_event_reestim;
    delete length_bias_reestim;
  }

  return inter_event;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Estimation of an equilibrium renewal process on the basis of time interval data
 *         using the EM algorithm.
 *
 *  \param[in] error          reference on a StatError object,
 *  \param[in] os             stream,
 *  \param[in] backward       backward recurrence time frequency distribution,
 *  \param[in] forward        forward recurrence time frequency distribution,
 *  \param[in] no_event       observation period frequency distribution for the case of no event,
 *  \param[in] estimator      estimator type (likelihood or penalized likelihood),
 *  \param[in] nb_iter        number of iterations,
 *  \param[in] mean_estimator method for the computation of the mean of the inter-event distribution,
 *  \param[in] weight         penalty weight,
 *  \param[in] pen_type       penalty type,
 *  \param[in] outside        management of side effects (zero outside the support or
 *                            continuation of the distribution).
 *
 *  \return                   DiscreteParametricModel object.
 */
/*--------------------------------------------------------------*/

DiscreteParametricModel* FrequencyDistribution::estimation(StatError &error , ostream &os ,
                                                           const FrequencyDistribution &backward ,
                                                           const FrequencyDistribution &forward ,
                                                           const FrequencyDistribution *no_event ,
                                                           estimation_criterion estimator , int nb_iter ,
                                                           duration_distribution_mean_estimator mean_estimator ,
                                                           double weight , penalty_type pen_type , side_effect outside) const

{
  register int i;
  int nb_histo , *pfrequency;
  double *pmass;
  DiscreteParametric *iinter_event;
  DiscreteParametricModel *inter_event;
  FrequencyDistribution *interval;
  const FrequencyDistribution *phisto[4];


  nb_histo = 3;
  phisto[0] = this;
  phisto[1] = new FrequencyDistribution(backward , SHIFT , 1);
  phisto[2] = &forward;
  if (no_event) {
    nb_histo++;
    phisto[3] = new FrequencyDistribution(*no_event , SHIFT , 1);
  }

  interval = new FrequencyDistribution(nb_histo , phisto);
  delete phisto[1];
  if (no_event) {
    delete phisto[3];
  }

  iinter_event = new DiscreteParametric((int)(interval->nb_value * MAX_VALUE_COEFF));

  iinter_event->offset = interval->offset;

  pmass = iinter_event->mass;
  for (i = 0;i < interval->offset;i++) {
    *pmass++ = 0.;
  }

  pfrequency = interval->frequency + interval->offset;
  for (i = interval->offset;i < interval->nb_value;i++) {
    *pmass++ = (double)*pfrequency++ / (double)(interval->nb_element + 1);
  }

  for (i = interval->nb_value;i < iinter_event->nb_value - 1;i++) {
    *pmass++ = 0.;
  }
  *pmass = 1. / (double)(interval->nb_element + 1);

  iinter_event->cumul_computation();

  iinter_event->max = (double)max / (double)(interval->nb_element + 1);
  iinter_event->mean_computation();
  iinter_event->variance_computation();

  delete interval;

# ifdef DEBUG
  iinter_event->ascii_print(cout);
# endif

  inter_event = estimation(error , os , backward , forward , no_event ,
                           *iinter_event , estimator , nb_iter , mean_estimator ,
                           weight , pen_type , outside);
  delete iinter_event;

  return inter_event;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the log-likelihood of a state occupancy distribution of
 *         an ordinary semi-Markov chain.
 *
 *  \param[in] sojourn_time frequency distribution of complete sojourn times,
 *  \param[in] final_run    frequency distribution of right-censored sojourn times.
 *
 *  \return                 log-likelihood of the estimated distribution.
 */
/*--------------------------------------------------------------*/

double DiscreteParametric::state_occupancy_likelihood_computation(const FrequencyDistribution &sojourn_time ,
                                                                  const FrequencyDistribution &final_run) const

{
  double likelihood , buff;


  likelihood = likelihood_computation(sojourn_time);

  if (likelihood != D_INF) {
    buff = survivor_likelihood_computation(final_run);

    if (buff != D_INF) {
      likelihood += buff;
    }
    else {
      likelihood = D_INF;
    }
  }

  return likelihood;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the log-likelihood of a state occupancy distribution of
 *         an equilibrium semi-Markov chain.
 *
 *  \param[in] forward      forward sojourn time distribution,
 *  \param[in] sojourn_time frequency distribution of complete sojourn times,
 *  \param[in] initial_run  frequency distribution of left-censored sojourn times,
 *  \param[in] final_run    frequency distribution of right-censored sojourn times,
 *  \param[in] single_run   frequency distribution of sequence lengths for the case of
 *                          a single visited state.
 *
 *  \return                 log-likelihood of the estimated distribution.
 */
/*--------------------------------------------------------------*/

double DiscreteParametric::state_occupancy_likelihood_computation(const Forward &forward ,
                                                                  const FrequencyDistribution &sojourn_time ,
                                                                  const FrequencyDistribution &initial_run ,
                                                                  const FrequencyDistribution &final_run ,
                                                                  const FrequencyDistribution &single_run) const

{
  double likelihood , buff;


  likelihood = likelihood_computation(sojourn_time);

  if (likelihood != D_INF) {
    buff = survivor_likelihood_computation(final_run);

    if (buff != D_INF) {
      likelihood += buff;
      buff = forward.likelihood_computation(initial_run);

      if (buff != D_INF) {
        likelihood += buff;
        buff = forward.survivor_likelihood_computation(single_run);

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
  }

  return likelihood;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the reestimation quantities of a state occupancy distribution
 *         (EM estimator of an ordinary semi-Markov chain).
 *
 *  \param[in] sojourn_time      frequency distribution of complete sojourn times,
 *  \param[in] final_run         frequency distribution of right-censored sojourn times,
 *  \param[in] occupancy_reestim pointer on the reestimation quantities of the state occupancy distribution,
 *  \param[in] iter              EM iteration.
 */
/*--------------------------------------------------------------*/

void DiscreteParametric::expectation_step(const FrequencyDistribution &sojourn_time ,
                                          const FrequencyDistribution &final_run ,
                                          Reestimation<double> *occupancy_reestim , int iter) const

{
  register int i;
  int *pfrequency;
  double sum , *ofrequency , *pmass , *pcumul;


  // initializations

  ofrequency = occupancy_reestim->frequency;
  for (i = 0;i < occupancy_reestim->alloc_nb_value;i++) {
    *ofrequency++ = 0.;
  }

  // computation of the reestimation quantities of the state occupancy distribution

  ofrequency = occupancy_reestim->frequency + sojourn_time.offset;
  pfrequency = sojourn_time.frequency + sojourn_time.offset;
  for (i = sojourn_time.offset;i < sojourn_time.nb_value;i++) {
    *ofrequency++ += *pfrequency++;
  }

  pfrequency = final_run.frequency + final_run.offset;
  pcumul = cumul + final_run.offset - 1;
  ofrequency = occupancy_reestim->frequency + final_run.offset;
  pmass = mass + final_run.offset;
  sum = 0.;

  for (i = final_run.offset;i < final_run.nb_value;i++) {
    sum += *pfrequency++ / (1. - *pcumul++);
    *ofrequency++ += *pmass++ * sum;
  }
  for (i = final_run.nb_value;i < nb_value;i++) {
    *ofrequency++ += *pmass++ * sum;
  }

  occupancy_reestim->nb_value_computation();
  occupancy_reestim->offset_computation();
  occupancy_reestim->nb_element_computation();

# ifdef DEBUG
  if ((iter < 10) || ((iter < 100) && (iter % 10 == 0)) ||
      ((iter < 1000) && (iter % 100 == 0)) || (iter % 1000 == 0)) {
    occupancy_reestim->max_computation();
    occupancy_reestim->mean_computation();
    occupancy_reestim->variance_computation();

    cout << "\nreestimation quantities of the state occupancy distribution:" << *occupancy_reestim << endl;
  }
# endif

}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the reestimation quantities of a state occupancy distribution
 *         (EM estimator of an equilibrium semi-Markov chain).
 *
 *  \param[in] sojourn_time        frequency distribution of complete sojourn times,
 *  \param[in] initial_run         frequency distribution of left-censored sojourn times,
 *  \param[in] final_run           frequency distribution of right-censored sojourn times,
 *  \param[in] single_run          frequency distribution of sequence lengths for the case of
 *                                 a single visited state,
 *  \param[in] occupancy_reestim   pointer on the reestimation quantities of the state occupancy distribution,
 *  \param[in] length_bias_reestim pointer on the reestimation quantities of the length-biased distribution,
 *  \param[in] iter                EM iteration,
 *  \param[in] combination         combination or not of the reestimation quantities,
 *  \param[in] mean_estimator      method for the computation of the mean of the state occupancy distribution.
 */
/*--------------------------------------------------------------*/

void DiscreteParametric::expectation_step(const FrequencyDistribution &sojourn_time ,
                                          const FrequencyDistribution &initial_run ,
                                          const FrequencyDistribution &final_run ,
                                          const FrequencyDistribution &single_run ,
                                          Reestimation<double> *occupancy_reestim ,
                                          Reestimation<double> *length_bias_reestim ,
                                          int iter , bool combination ,
                                          duration_distribution_mean_estimator mean_estimator) const

{
  register int i , j;
  int reestim_offset , reestim_nb_value , *pfrequency;
  double sum , occupancy_mean , *ofrequency , *lfrequency , *pmass , *pcumul , *norm;


  // initializations

  ofrequency = occupancy_reestim->frequency;
  lfrequency = length_bias_reestim->frequency;
  for (i = 0;i < occupancy_reestim->alloc_nb_value;i++) {
    *ofrequency++ = 0.;
    *lfrequency++ = 0.;
  }

  // computation of the reestimation quantities of the state occupancy distribution

  ofrequency = occupancy_reestim->frequency + sojourn_time.offset;
  pfrequency = sojourn_time.frequency + sojourn_time.offset;
  for (i = sojourn_time.offset;i < sojourn_time.nb_value;i++) {
    *ofrequency++ += *pfrequency++;
  }

  if (final_run.nb_element > 0) {
    pfrequency = final_run.frequency + final_run.offset;
    pcumul = cumul + final_run.offset - 1;
    ofrequency = occupancy_reestim->frequency + final_run.offset;
    pmass = mass + final_run.offset;
    sum = 0.;

    for (i = final_run.offset;i < final_run.nb_value;i++) {
      sum += *pfrequency++ / (1. - *pcumul++);
      *ofrequency++ += *pmass++ * sum;
    }
    for (i = final_run.nb_value;i < nb_value;i++) {
      *ofrequency++ += *pmass++ * sum;
    }
  }

  // computation of the reestimation quantities of the length-biased distribution

  if (initial_run.nb_element > 0) {
    pfrequency = initial_run.frequency + initial_run.offset;
    pcumul = cumul + initial_run.offset - 1;
    lfrequency = length_bias_reestim->frequency + initial_run.offset;
    pmass = mass + initial_run.offset;
    sum = 0.;

    for (i = initial_run.offset;i < initial_run.nb_value;i++) {
      sum += *pfrequency++ / (1. - *pcumul++);
      *lfrequency++ += *pmass++ * sum;
    }
    for (i = initial_run.nb_value;i < nb_value;i++) {
      *lfrequency++ += *pmass++ * sum;
    }
  }

  if (single_run.nb_element > 0) {
    norm = new double[single_run.nb_value];

    for (i = single_run.offset;i < single_run.nb_value;i++) {
      if (single_run.frequency[i] > 0) {
        pmass = mass + i;
        norm[i] = 0.;
        for (j = i;j < nb_value;j++) {
          norm[i] += (j + 1 - i) * *pmass++;
        }
      }
    }

    lfrequency = length_bias_reestim->frequency + single_run.offset;
    pmass = mass + single_run.offset;
    for (i = single_run.offset;i < nb_value;i++) {
      pfrequency = single_run.frequency + single_run.offset;
      sum = 0.;
      for (j = single_run.offset;j <= MIN(i , single_run.nb_value - 1);j++) {
        if ((*pfrequency > 0) && (norm[j] > 0.)) {
          sum += *pfrequency * (i + 1 - j) / norm[j];
               }
        pfrequency++;
      }

      *lfrequency++ += *pmass++ * sum;
    }

    delete [] norm;
  }

  reestim_offset = 1;
  reestim_nb_value = occupancy_reestim->alloc_nb_value;

  ofrequency = occupancy_reestim->frequency + occupancy_reestim->alloc_nb_value;
  lfrequency = length_bias_reestim->frequency + occupancy_reestim->alloc_nb_value;
  while ((*--ofrequency == 0) && (*--lfrequency == 0) && (reestim_nb_value > 2)) {
    reestim_nb_value--;
  }
  occupancy_reestim->nb_value = reestim_nb_value;
  length_bias_reestim->nb_value = reestim_nb_value;

  ofrequency = occupancy_reestim->frequency + reestim_offset;
  lfrequency = length_bias_reestim->frequency + reestim_offset;
  while ((*ofrequency++ == 0) && (*lfrequency++ == 0) && (reestim_offset < reestim_nb_value - 1)) {
    reestim_offset++;
  }
  occupancy_reestim->offset = reestim_offset;
  length_bias_reestim->offset = reestim_offset;

  occupancy_reestim->nb_element_computation();
  length_bias_reestim->nb_element_computation();

# ifdef DEBUG
  if ((iter < 10) || ((iter < 100) && (iter % 10 == 0)) ||
      ((iter < 1000) && (iter % 100 == 0)) || (iter % 1000 == 0)) {
    occupancy_reestim->max_computation();
    occupancy_reestim->mean_computation();
    occupancy_reestim->variance_computation();

    length_bias_reestim->max_computation();
    length_bias_reestim->mean_computation();
    length_bias_reestim->variance_computation();

    cout << "\nreestimation quantities of the state occupancy distribution:" << *occupancy_reestim << endl;
    cout << "\nreestimation quantities of the length-biased distribution:" << *length_bias_reestim << endl;
  }
# endif

  if (combination) {
    switch (mean_estimator) {
    case COMPUTED :
      occupancy_mean = interval_bisection(occupancy_reestim , length_bias_reestim);
      break;
    case ONE_STEP_LATE :
      occupancy_mean = mean;
      break;
    }

#   ifdef DEBUG
    cout << STAT_label[STATL_SOJOURN_TIME] << " " << STAT_label[STATL_MEAN] << ": " << occupancy_mean << endl;
#   endif

    occupancy_reestim->equilibrium_process_combination(length_bias_reestim , occupancy_mean);

#   ifdef DEBUG
    if ((iter < 10) || ((iter < 100) && (iter % 10 == 0)) ||
        ((iter < 1000) && (iter % 100 == 0)) || (iter % 1000 == 0)) {
      cout << "\nreestimation quantities of the state occupancy distribution:" << *occupancy_reestim << endl;
    }
#   endif

  }
}


};  // namespace stat_tool
