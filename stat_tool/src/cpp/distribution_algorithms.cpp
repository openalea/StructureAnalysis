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



/*--------------------------------------------------------------*
 *
 *  Convolution de 2 lois (le resultat de la convolution peut 
 *  etre mis dans l'une des 2 lois).
 *
 *  arguments : references sur les 2 lois, nombre de valeurs
 *              de la loi resultante.
 *
 *--------------------------------------------------------------*/

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


/*--------------------------------------------------------------*
 *
 *  Calcul d'une loi binomiale.
 *
 *  arguments : nombre de valeurs,
 *              mode de calcul ('s' : standard, 'r' : renouvellement).
 *
 *--------------------------------------------------------------*/

void DiscreteParametric::binomial_computation(int inb_value , char mode)

{
  register int i;
  int set , subset;
  double failure = 1. - probability , success = probability , ratio , scale , term;


  switch (mode) {
  case 's' :
    offset = inf_bound;
    nb_value = sup_bound + 1;
    break;
  case 'r' :
    offset = MIN(inf_bound , inb_value - 1);
    nb_value = MIN(sup_bound + 1 , inb_value);
    break;
  }

  // valeurs de probabilite nulle avant la borne inferieure

  for (i = 0;i < MIN(nb_value , inf_bound);i++) {
    mass[i] = 0.;
  }

  if (nb_value > inf_bound) {
    set = sup_bound - inf_bound;

    if (probability <= B_PROBABILITY) {
      subset = 0;

      // cas calcul direct

      if ((sup_bound - inf_bound) / failure < B_THRESHOLD) {

        // calcul de la probabilite de la borne inferieure

        term = 1.;
        for (i = 0;i < sup_bound - inf_bound;i++) {
          term *= failure;
        }
        mass[inf_bound] = term;

        // calcul des probabilites des valeurs suivantes

        ratio = success / failure;

        for (i = inf_bound + 1;i < nb_value;i++) {
          scale = (double)(set - subset) / (double)(subset + 1);
          subset++;
          term *= scale * ratio;
          mass[i] = term;
        }
      }

      // cas calcul en log

      else {

        // calcul de la probabilite de la borne inferieure

        term = (sup_bound - inf_bound) * log(failure);
        mass[inf_bound] = exp(term);

        // calcul des probabilites des valeurs suivantes

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

      // cas calcul direct

      if ((sup_bound - inf_bound) / success < B_THRESHOLD) {

        // calcul de la probabilite de la borne superieure

        term = 1.;
        for (i = 0;i < sup_bound - inf_bound;i++) {
          term *= success;
        }
        if (sup_bound < nb_value) {
          mass[sup_bound] = term;
        }

        // calcul des probabilites des valeurs precedentes

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

      // cas calcul en log

      else {

        // calcul de la probabilite de la borne superieure

        term = (sup_bound - inf_bound) * log(success);
        if (sup_bound < nb_value) {
          mass[sup_bound] = exp(term);
        }

        // calcul des probabilites des valeurs precedentes

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

# ifdef DEBUG
  if (mode == 's') {
    binomial dist(sup_bound - inf_bound , probability);


    cout << "\nTEST binomial distribution" << endl;
    for (i = inf_bound;i <= sup_bound;i++) {
      cout << i << "  " << pdf(dist , i - inf_bound) << " | " << mass[i]
           << "   " << cdf(dist , i - inf_bound) << " | " << cumul[i] << endl;
    }
  }
# endif

}


/*--------------------------------------------------------------*
 *
 *  Calcul d'une loi de Poisson, loi definie sur
 *  [borne inferieure , ... , infini]. On fixe le nombre de valeurs,
 *  soit a partir d'un seuil calcule sur la fonction de repartition,
 *  soit a partir d'une borne predefinie.
 *
 *  arguments : nombre de valeurs,
 *              seuil sur la fonction de repartition,
 *              mode de calcul ('s' : standard, 'r' : renouvellement).
 *
 *--------------------------------------------------------------*/

void DiscreteParametric::poisson_computation(int inb_value , double cumul_threshold ,
                                             char mode)

{
  register int i;
  double log_parameter , num , denom , *pmass , *pcumul;


  pmass = mass;
  pcumul = cumul;

  switch (mode) {

  // cas calcul complet

  case 's' : {

    // valeurs de probabilite nulle avant la borne inferieure

    for (i = 0;i < inf_bound;i++) {
      *pmass++ = 0.;
      *pcumul++ = 0.;
    }

    i = 1;

    // cas calcul direct

    if (parameter < P_THRESHOLD) {

      // calcul de la probabilite de la borne inferieure

      *pmass = exp(-parameter);
      *pcumul = *pmass;

      // calcul des probabilites des valeurs suivantes

      while (((*pcumul < cumul_threshold) || (i + inf_bound < inb_value)) &&
             (i + inf_bound < alloc_nb_value)) {
        pmass++;
        pcumul++;
        *pmass = *(pmass - 1) * parameter / i;
        *pcumul = *(pcumul - 1) + *pmass;
        i++;
      }
    }

    // cas calcul en log

    else {

      // calcul de la probabilite de la borne inferieure

      num = -parameter;
      *pmass = exp(num);
      *pcumul = *pmass;

      // calcul des probabilites des valeurs suivantes

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

  // cas calcul incomplet (renouvellement)

  case 'r' : {

    // valeurs de probabilite nulle avant la borne inferieure

    for (i = 0;i < MIN(inb_value , inf_bound);i++) {
      *pmass++ = 0.;
      *pcumul++ = 0.;
    }

    if (inb_value > inf_bound) {
      i = 1;

      // cas calcul direct

      if (parameter < P_THRESHOLD) {

        // calcul de la probabilite de la borne inferieure

        *pmass = exp(-parameter);
        *pcumul = *pmass;

        // calcul des probabilites des valeurs suivantes

        while ((*pcumul < cumul_threshold) && (i + inf_bound < inb_value)) {
          pmass++;
          pcumul++;
          *pmass = *(pmass - 1) * parameter / i;
          *pcumul = *(pcumul - 1) + *pmass;
          i++;
        }
      }

      // cas calcul en log

      else {

        // calcul de la probabilite de la borne inferieure

        num = -parameter;
        *pmass = exp(num);
        *pcumul = *pmass;

        // calcul des probabilites des valeurs suivantes

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

# ifdef DEBUG
  if (mode == 's') {
    poisson dist(parameter);


    cout << "\nTEST Poisson distribution" << endl;
    for (i = inf_bound;i < nb_value;i++) {
      cout << i << "  " << pdf(dist , i - inf_bound) << " | " << mass[i]
           << "   " << cdf(dist , i - inf_bound) << " | " << cumul[i] << endl;
    }
  }
# endif

}


/*--------------------------------------------------------------*
 *
 *  Calcul d'une loi binomiale negative, loi definie sur
 *  [borne inferieure , ... , infini]. On fixe le nombre de valeurs,
 *  soit a partir d'un seuil calcule sur la fonction de repartition,
 *  soit a partir d'une borne predefinie.
 *
 *  arguments : nombre de valeurs, seuil sur la fonction de repartition,
 *              mode de calcul ('s' : standard, 'r' : renouvellement).
 *
 *--------------------------------------------------------------*/

void DiscreteParametric::negative_binomial_computation(int inb_value , double cumul_threshold ,
                                                       char mode)

{
  register int i;
  double failure = 1. - probability , success = probability , log_failure ,
         set , subset , scale , term , *pmass , *pcumul;


  pmass = mass;
  pcumul = cumul;

  switch (mode) {

    // cas calcul complet

    case 's' : {

    // valeurs de probabilite nulle avant la borne inferieure

    for (i = 0;i < inf_bound;i++) {
      *pmass++ = 0.;
      *pcumul++ = 0.;
    }

    i++;
    subset = parameter - 1.;
    set = subset;

    // cas calcul direct

    if (sqrt(parameter) / success < NB_THRESHOLD) {

      // calcul de la probabilite de la borne inferieure

      term = pow(success , parameter);
      *pmass = term;
      *pcumul = *pmass;

      // calcul des probabilites des valeurs suivantes

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

    // cas calcul en log

    else {

      // calcul de la probabilite de la borne inferieure

      term = parameter * log(success);
      *pmass = exp(term);
      *pcumul = *pmass;

      // calcul des probabilites des valeurs suivantes

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

  // cas calcul incomplet (renouvellement)

  case 'r' : {

    // valeurs de probabilite nulle avant la borne inferieure

    for (i = 0;i < MIN(inf_bound , inb_value);i++) {
      *pmass++ = 0.;
      *pcumul++ = 0.;
    }

    if (inb_value > inf_bound) {
      i++;
      subset = parameter - 1.;
      set = subset;

      // cas calcul direct

      if (sqrt(parameter) / success < NB_THRESHOLD) {

        // calcul de la probabilite de la borne inferieure

        term = pow(success , parameter);
        *pmass = term;
        *pcumul = *pmass;

        // calcul des probabilites de valeurs suivantes

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

      // cas calcul en log

      else {

        // calcul de la probabilite de la borne inferieure

        term = parameter * log(success);
        *pmass = exp(term);
        *pcumul = *pmass;

        // calcul des probabilites des valeurs suivantes

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

# ifdef DEBUG
  if (mode == 's') {
    negative_binomial dist(parameter , probability);


    cout << "TEST negative binomial distribution" << endl;
    for (i = inf_bound;i < nb_value;i++) {
      cout << i << "  " << pdf(dist , i - inf_bound) << " | " << mass[i]
           << "   " << cdf(dist , i - inf_bound) << " | " << cumul[i] << endl;
    }
  }
# endif

}


/*--------------------------------------------------------------*
 *
 *  Calcul d'une loi uniforme.
 *
 *--------------------------------------------------------------*/

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


/*--------------------------------------------------------------*
 *
 *  Calcul du nombre de valeurs d'une loi discrete elementaire.
 *
 *  arguments : identificateur, bornes inferieure et superieure,
 *              parametre, probabilite,
 *              seuil sur la fonction de repartition.
 *
 *--------------------------------------------------------------*/

int nb_value_computation(int ident , int inf_bound , int sup_bound ,
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


/*--------------------------------------------------------------*
 *
 *  Calcul d'une loi discrete elementaire (binomiale, Poisson,
 *  binomiale negative, ou uniforme).
 *
 *  arguments : nombre minimum de valeurs,
 *              seuil sur la fonction de repartition.
 *
 *--------------------------------------------------------------*/

void DiscreteParametric::computation(int min_nb_value , double cumul_threshold)

{
  if (ident > 0) {
    switch (ident) {
    case BINOMIAL :
      binomial_computation(1 , 's');
      break;
    case POISSON :
      poisson_computation(min_nb_value , cumul_threshold , 's');
      break;
    case NEGATIVE_BINOMIAL :
      negative_binomial_computation(min_nb_value , cumul_threshold , 's');
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


/*--------------------------------------------------------------*
 *
 *  Calcul de la loi de l'intervalle de temps residuel a partir
 *  de la loi de temps de retour/sejour.
 *
 *  argument : loi de temps de retour/sejour.
 *
 *--------------------------------------------------------------*/

void Forward::computation(const DiscreteParametric &dist)

{
  register int i;
  double norm;


  offset = 1;
  nb_value = dist.nb_value;
  mass[0] = 0.;

  // calcul de la quantite de normalisation

  if (ident == CATEGORICAL) {
    norm = dist.mean;
  }
  else {
    norm = parametric_mean_computation();
  }

  // calcul des probabilites des valeurs

  for (i = 1;i < nb_value;i++) {
    mass[i] = (1. - dist.cumul[i - 1]) / norm;
  }

  cumul_computation();

  max_computation();
  mean_computation();
  variance_computation();
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la vraisemblance pour la fonction de survie d'une loi donnee.
 *
 *  argument : reference sur un objet FrequencyDistribution.
 *
 *--------------------------------------------------------------*/

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


/*--------------------------------------------------------------*
 *
 *  Calcul de la valeur du chi2 pour une loi.
 *
 *  argument : reference sur un objet FrequencyDistribution.
 *
 *--------------------------------------------------------------*/

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


/*--------------------------------------------------------------*
 *
 *  Regroupement des valeurs, calcul du nombre de degres de liberte et
 *  de la valeur du chi2.
 *
 *  arguments : reference sur un objet FrequencyDistribution et sur un objet Test.
 *
 *--------------------------------------------------------------*/

void Distribution::chi2_degree_of_freedom(const FrequencyDistribution &histo , Test &test) const

{
  register int i , j;
  int *filter_frequency;
  double *filter_mass;
  Distribution *filter_dist;
  FrequencyDistribution *filter_histo;


  if ((histo.offset >= offset) && (histo.nb_value <= nb_value)) {

    // creation et initialisation des objets Distribution et FrequencyDistribution

    filter_dist = new Distribution(nb_value);
    filter_dist->offset = offset;
    filter_dist->nb_parameter = nb_parameter;

    filter_histo = new FrequencyDistribution(histo.nb_value);
    filter_histo->nb_element = histo.nb_element;

    // regroupement des valeurs

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

    // mise a jour du nombre de degres de libertes et calcul de
    // la valeur du chi2

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


/*--------------------------------------------------------------*
 *
 *  Test d'ajustement d'une loi.
 *
 *  arguments : reference sur un objet FrequencyDistribution et sur un objet Test.
 *
 *--------------------------------------------------------------*/

void Distribution::chi2_fit(const FrequencyDistribution &histo , Test &test) const

{
  if ((histo.offset >= offset) && (histo.nb_value <= nb_value)) {
    chi2_degree_of_freedom(histo , test);
    if ((test.df1 > 0) && (test.value > 0.)) {
      test.chi2_critical_probability_computation();
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Troncature d'une loi.
 *
 *  arguments : reference sur un objet StatError, valeur maximum.
 *
 *--------------------------------------------------------------*/

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

    // creation d'un objet DiscreteParametricModel

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


/*--------------------------------------------------------------*
 *
 *  Ajustement d'une loi.
 *
 *  arguments : references sur un objet StatError et
 *              sur un objet DiscreteParametric.
 *
 *--------------------------------------------------------------*/

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


/*--------------------------------------------------------------*
 *
 *  Estimation des parametres d'une loi discrete elementaire
 *  (binomiale negative, binomiale, Poisson).
 *
 *  arguments : identificateur de la loi, borne inferieure minimum,
 *              flag sur la borne inferieure, seuil sur la fonction de repartition.
 *
 *--------------------------------------------------------------*/

DiscreteParametric* FrequencyDistribution::parametric_estimation(int ident , int min_inf_bound ,
                                                                 bool flag , double cumul_threshold) const

{
  double likelihood;
  DiscreteParametric *dist;


  // creation d'un objet DiscreteParametric

  dist = new DiscreteParametric((int)(nb_value * SAMPLE_NB_VALUE_COEFF) , ident);

  // estimation des parametres de la loi

  likelihood = Reestimation<int>::parametric_estimation(dist , min_inf_bound ,
                                                        flag , cumul_threshold);

  // mise a jour de la loi estimee

  if (likelihood != D_INF) {
    dist->computation(nb_value , cumul_threshold);
  }
  else {
    delete dist;
    dist = NULL;
  }

  return dist;
}


/*--------------------------------------------------------------*
 *
 *  Estimation des parametres d'une loi discrete elementaire
 *  (binomiale negative, binomiale, Poisson).
 *
 *  arguments : reference sur un objet StatError, identificateur de la loi,
 *              borne inferieure minimum, flag sur la borne inferieure,
 *              seuil sur la fonction de repartition.
 *
 *--------------------------------------------------------------*/

DiscreteParametricModel* FrequencyDistribution::parametric_estimation(StatError &error , int ident ,
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

    // creation d'un objet DiscreteParametricModel

    dist = new DiscreteParametricModel((int)(nb_value * SAMPLE_NB_VALUE_COEFF) , ident);
    dist->frequency_distribution = new DiscreteDistributionData(*this);

    // estimation des parametres de la loi

    likelihood = Reestimation<int>::parametric_estimation(dist , min_inf_bound ,
                                                          flag , cumul_threshold);

    if (likelihood != D_INF) {

      // mise a jour de la loi estimee

      dist->computation(nb_value , cumul_threshold);

      // mise a jour du nombre de parametres inconnus

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


/*--------------------------------------------------------------*
 *
 *  Estimation des parametres d'une loi discrete elementaire
 *  (binomiale negative, binomiale, Poisson).
 *
 *  arguments : reference sur un objet StatError, borne inferieure minimum,
 *              flag sur la borne inferieure, seuil sur la fonction de repartition.
 *
 *--------------------------------------------------------------*/

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

    // creation d'un objet DiscreteParametricModel

    dist = new DiscreteParametricModel((int)(nb_value * SAMPLE_NB_VALUE_COEFF));
    dist->frequency_distribution = new DiscreteDistributionData(*this);

    // estimation des parametres de la loi

    likelihood = Reestimation<int>::type_parametric_estimation(dist , min_inf_bound ,
                                                               flag , cumul_threshold);

    if (likelihood != D_INF) {

      // mise a jour de la loi estimee

      dist->computation(nb_value , cumul_threshold);

      // mise a jour du nombre de parametres inconnus

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


/*--------------------------------------------------------------*
 *
 *  Calcul des termes de penalisation dans le cadre d'une approche
 *  de penalisation de la vraisemblance.
 *
 *  arguments : poids de la penalisation, type de penalisation (difference 1ere,
 *              seconde ou entropy), type de gestion des effets de bord
 *              (zero a l'exterieur du support ou prolongation de la loi).
 *
 *--------------------------------------------------------------*/

void Distribution::penalty_computation(double weight , int type , double *penalty , int outside) const

{
  register int i;


  switch (type) {

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


/*--------------------------------------------------------------*
 *
 *  Reestimation des parametres d'une loi discrete elementaire
 *  (binomiale, Poisson, binomiale negative).
 *
 *  arguments : reference sur les quantites de reestimation,
 *              nombre de parametres reestimes (binomiale negative).
 *
 *--------------------------------------------------------------*/

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


/*--------------------------------------------------------------*
 *
 *  Simulation par une loi discrete en utilisant la fonction de repartition.
 *  On teste le terme median de la fonction de repartition pour savoir
 *  ou commencer la recherche.
 *
 *  arguments : nombre de valeurs, pointeur sur la fonction de repartition,
 *              facteur d'echelle.
 *
 *--------------------------------------------------------------*/

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


/*--------------------------------------------------------------*
 *
 *  Simulation par une loi discrete en utilisant la fonction de repartition
 *
 *--------------------------------------------------------------*/

int Distribution::simulation() const

{
  int value;


  value = offset + cumul_method(nb_value - offset , cumul + offset , 1. - complement);

  return value;
}


/*--------------------------------------------------------------*
 *
 *  Simulation par une loi discrete par la methode du rejet.
 *  Principe : on tire un point (x,Px) dans le rectangle
 *  [xmin,xmax] [0.,Pmax]. Si le point est sous la courbe, on
 *  garde la realisation correspondante (abscisse x), sinon,
 *  on retire un nouveau point.
 *
 *--------------------------------------------------------------*/

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


/*--------------------------------------------------------------*
 *
 *  Constitution d'un echantillon par simulation d'une loi parametrique.
 *
 *  arguments : reference sur un objet StatError, effectif.
 *
 *--------------------------------------------------------------*/

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

    // creation de la loi empirique

    histo = new DiscreteDistributionData(*this);
    histo->distribution = new DiscreteParametricModel(*this , false);

    // simulation

    for (i = 0;i < nb_element;i++) {
      (histo->frequency[DiscreteParametric::simulation()])++;
    }

    // extraction des caracteristiques de la loi empirique

    histo->nb_value_computation();
    histo->offset_computation();
    histo->nb_element = nb_element;
    histo->max_computation();
    histo->mean_computation();
    histo->variance_computation();
  }

  return histo;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la vraisemblance des donnees d'intervalles de temps.
 *
 *  arguments : loi de l'intervalle de temps residuel, lois empiriques des intervalles
 *              de temps complets, censures a gauche, a droite et
 *              de la longueur de la periode d'observation dans le cas 0 evenement.
 *
 *--------------------------------------------------------------*/

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
    histo = new FrequencyDistribution(backward , 's' , 1);
    buff = survivor_likelihood_computation(*histo);
    delete histo;

    if (buff != D_INF) {
      likelihood += buff;
      buff = forward_dist.likelihood_computation(forward);

      if (buff != D_INF) {
        likelihood += buff;

        if (no_event) {
          histo = new FrequencyDistribution(*no_event , 's' , 1);
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


/*--------------------------------------------------------------*
 *
 *  Calcul des quantites de reestimation correspondant a la loi inter-evenement
 *  (estimateur EM d'un processus de renouvellement en equilibre a partir
 *   de donnees d'intervalles de temps).
 *
 *  arguments : lois empiriques des intervalles de temps complets, censures a gauche,
 *              a droite et de la longueur de la periode d'observation dans le cas 0 evenement,
 *              pointeurs sur les quantites de reestimation de la loi inter-evenement et
 *              de la loi biaisee par la longueur.
 *
 *--------------------------------------------------------------*/

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


  // initialisations

  ifrequency = inter_event_reestim->frequency;
  lfrequency = length_bias_reestim->frequency;
  for (i = 0;i < inter_event_reestim->alloc_nb_value;i++) {
    *ifrequency++ = 0.;
    *lfrequency++ = 0.;
  }

  // calcul des quantites de reestimation de la loi inter-evenement

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

  // calcul des quantites de reestimation de la loi biaisee par la longueur

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

    cout << "\nquantites de reestimation loi inter_evenement :" << *inter_event_reestim << endl;
    cout << "\nquantites de reestimation loi biaisee par la longueur :" << *length_bias_reestim << endl;
  }
# endif

}


/*--------------------------------------------------------------*
 *
 *  Calcul de la moyenne d'une loi par bissection d'intervalle.
 *
 *  arguments : pointeurs sur les quantites de reestimation de la loi et
 *              de la loi biaisee par la longueur.
 *
 *--------------------------------------------------------------*/

double interval_bisection(Reestimation<double> *distribution_reestim ,
                          Reestimation<double> *length_bias_reestim)

{
  register int i;
  double ratio , inf_ratio , sup_ratio , mean , inf_mean , sup_mean , *dfrequency , *lfrequency;

# ifdef DEBUG
  int iter = 0;
# endif


  // initialisations : calculs des 2 premieres valeurs

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


/*--------------------------------------------------------------*
 *
 *  Estimation des parametres d'un processus de renouvellement en equilibre
 *  par l'algorithme EM a partir de donnees d'intervalles de temps.
 *
 *  arguments : reference sur un objet StatError, stream, lois empiriques des intervalles
 *              de temps censures a gauche, a droite et de la longueur
 *              de la periode d'observation dans le cas 0 evenement,
 *              reference sur la loi inter-evenement initiale, type d'estimateur
 *              (vraisemblance ou vraisemblance penalisee), nombre d'iterations,
 *              methode de calcul de la moyenne de la loi inter-evenement,
 *              poids de la penalisation, type de penalisation,
 *              type de gestion des effets de bord (zero a l'exterieur du support ou
 *              prolongation de la loi), moyenne de la loi inter-evenement.
 *
 *--------------------------------------------------------------*/

DiscreteParametricModel* FrequencyDistribution::estimation(StatError &error , ostream &os ,
                                                           const FrequencyDistribution &backward ,
                                                           const FrequencyDistribution &forward ,
                                                           const FrequencyDistribution *no_event ,
                                                           const DiscreteParametric &iinter_event ,
                                                           int estimator , int nb_iter ,
                                                           int mean_computation_method , double weight ,
                                                           int penalty_type , int outside ,
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

  if ((mean_computation_method == ESTIMATED) && (iinter_event_mean == D_DEFAULT)) {
    status = false;
    error.update(STAT_error[STATR_MEAN_COMPUTATION_METHOD]);
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
    phisto[0] = new FrequencyDistribution(backward , 's' , 1);
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

    // creation de la loi inter-evenement

    inter_event = new DiscreteParametricModel(iinter_event , this);
    inter_event->init(CATEGORICAL , I_DEFAULT , I_DEFAULT , D_DEFAULT , D_DEFAULT);
    forward_dist = new Forward(*inter_event);

    if (estimator == PENALIZED_LIKELIHOOD) {
      penalty = new double[inter_event->nb_value];

      if (weight == D_DEFAULT) {
        if (penalty_type != ENTROPY) {
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
        switch (mean_computation_method) {
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
        switch (mean_computation_method) {
        case ESTIMATED :
          inter_event_mean = iinter_event_mean;
          break;
        case ONE_STEP_LATE :
          inter_event_mean = inter_event->mean;
          break;
        }

        inter_event_reestim->penalized_likelihood_equilibrium_process_estimation(length_bias_reestim ,
                                                                                 inter_event , inter_event_mean ,
                                                                                 weight , penalty_type , penalty ,
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
          if (mean_computation_method == ESTIMATED) {
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
        if (mean_computation_method == ESTIMATED) {
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


/*--------------------------------------------------------------*
 *
 *  Estimation des parametres d'un processus de renouvellement en equilibre
 *  par l'algorithme EM a partir de donnees d'intervalles de temps.
 *
 *  arguments : reference sur un objet StatError, stream, lois empiriques des intervalles
 *              de temps censures a gauche, a droite et de la longueur de la periode
 *              d'observation dans le cas 0 evenement, type d'estimateur
 *              (vraisemblance ou vraisemblance penalisee), nombre d'iterations,
 *              methode de calcul de la moyenne de la loi inter-evenement,
 *              poids de la penalisation, type de penalisation,
 *              type de gestion des effets de bord (zero a l'exterieur du support ou
 *              prolongation de la loi).
 *
 *--------------------------------------------------------------*/

DiscreteParametricModel* FrequencyDistribution::estimation(StatError &error , ostream &os ,
                                                           const FrequencyDistribution &backward ,
                                                           const FrequencyDistribution &forward ,
                                                           const FrequencyDistribution *no_event ,
                                                           int estimator , int nb_iter ,
                                                           int mean_computation_method , double weight ,
                                                           int penalty_type , int outside) const

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
  phisto[1] = new FrequencyDistribution(backward , 's' , 1);
  phisto[2] = &forward;
  if (no_event) {
    nb_histo++;
    phisto[3] = new FrequencyDistribution(*no_event , 's' , 1);
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
                           *iinter_event , estimator , nb_iter , mean_computation_method ,
                           weight , penalty_type , outside);
  delete iinter_event;

  return inter_event;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la vraisemblance des temps de sejour dans un etat pour
 *  une semi-chaine de Markov ordinaire.
 *
 *  arguments : lois empiriques des temps de sejour complets et censures a droite.
 *
 *--------------------------------------------------------------*/

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


/*--------------------------------------------------------------*
 *
 *  Calcul de la vraisemblance des temps de sejour dans un etat pour
 *  une semi-chaine de Markov en equilibre.
 *
 *  arguments : loi de l'intervalle de temps residuel, lois empiriques des temps
 *              de sejour complets, censures a droite, a gauche et
 *              de la longueur des sequences dans le cas d'un seul etat visite.
 *
 *--------------------------------------------------------------*/

double DiscreteParametric::state_occupancy_likelihood_computation(const Forward &forward ,
                                                                  const FrequencyDistribution &sojourn_time ,
                                                                  const FrequencyDistribution &final_run ,
                                                                  const FrequencyDistribution &initial_run ,
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


/*--------------------------------------------------------------*
 *
 *  Calcul des quantites de reestimation correspondant a la loi d'occupation d'un etat
 *  (estimateur EM d'une semi-chaine de Markov ordinaire).
 *
 *  arguments : lois empiriques des temps de sejour complets et censures a droite,
 *              pointeurs sur les quantites de reestimation de la loi d'occupation de l'etat.
 *
 *--------------------------------------------------------------*/

void DiscreteParametric::expectation_step(const FrequencyDistribution &sojourn_time ,
                                          const FrequencyDistribution &final_run ,
                                          Reestimation<double> *occupancy_reestim , int iter) const

{
  register int i;
  int *pfrequency;
  double sum , *ofrequency , *pmass , *pcumul;


  // initialisations

  ofrequency = occupancy_reestim->frequency;
  for (i = 0;i < occupancy_reestim->alloc_nb_value;i++) {
    *ofrequency++ = 0.;
  }

  // calcul des quantites de reestimation de la loi d'occupation de l'etat

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

    cout << "\nquantites de reestimation loi d'occupation de l'etat :" << *occupancy_reestim << endl;
  }
# endif

}


/*--------------------------------------------------------------*
 *
 *  Calcul des quantites de reestimation correspondant a la loi d'occupation d'un etat
 *  (estimateur EM d'une semi-chaine de Markov en equilibre).
 *
 *  arguments : lois empiriques des temps de sejour complets, censures a droite, a gauche et
 *              de la longueur des sequences dans le cas d'un seul etat visite,
 *              pointeurs sur les quantites de reestimation de la loi d'occupation de l'etat et
 *              de la loi biaisee par la longueur, combinaison ou non des quantites de reestimation,
 *              methode de calcul de la moyenne de la loi d'occupation de l'etat.
 *
 *--------------------------------------------------------------*/

void DiscreteParametric::expectation_step(const FrequencyDistribution &sojourn_time ,
                                          const FrequencyDistribution &final_run ,
                                          const FrequencyDistribution &initial_run ,
                                          const FrequencyDistribution &single_run ,
                                          Reestimation<double> *occupancy_reestim ,
                                          Reestimation<double> *length_bias_reestim , int iter ,
                                          bool combination , int mean_computation_method) const

{
  register int i , j;
  int reestim_offset , reestim_nb_value , *pfrequency;
  double sum , occupancy_mean , *ofrequency , *lfrequency , *pmass , *pcumul , *norm;


  // initialisations

  ofrequency = occupancy_reestim->frequency;
  lfrequency = length_bias_reestim->frequency;
  for (i = 0;i < occupancy_reestim->alloc_nb_value;i++) {
    *ofrequency++ = 0.;
    *lfrequency++ = 0.;
  }

  // calcul des quantites de reestimation de la loi d'occupation de l'etat

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

  // calcul des quantites de reestimation de la loi biaisee par la longueur

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

    cout << "\nquantites de reestimation loi d'occupation de l'etat :" << *occupancy_reestim << endl;
    cout << "\nquantites de reestimation loi biaisee par la longueur :" << *length_bias_reestim << endl;
  }
# endif

  if (combination) {
    switch (mean_computation_method) {
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
      cout << "\nquantites de reestimation loi d'occupation de l'etat :" << *occupancy_reestim << endl;
    }
#   endif

  }
}


};  // namespace stat_tool
