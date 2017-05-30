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



#include "renewal.h"

using namespace std;
using namespace stat_tool;


namespace sequence_analysis {



/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the length-biased distribution from the inter-event distribution.
 *
 *  \param[in] inter_event inter-event distribution.
 */
/*--------------------------------------------------------------*/

void LengthBias::computation(const DiscreteParametric &inter_event)

{
  int i;
  double norm , *pmass , *imass;


  offset = inter_event.offset;
  nb_value = inter_event.nb_value;

  pmass = mass;
  for (i = 0;i < offset;i++) {
    *pmass++ = 0.;
  }
  imass = inter_event.mass + offset;

  // computation of the normalization quantity

  if (ident == CATEGORICAL) {
    norm = inter_event.mean;
  }
  else {
    norm = parametric_mean_computation();
  }

  // computation of the probability mass function

  for (i = offset;i < nb_value;i++) {
    *pmass++ = i * *imass++ / norm;
  }

  cumul_computation();

  max_computation();
  mean_computation();
  variance_computation();
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the backward recurrence time distribution from
 *         the inter-event distribution.
 *
 *  \param[in] inter_event inter-event distribution,
 *  \param[in] time        observation period distribution.
 */
/*--------------------------------------------------------------*/

void Backward::computation(const DiscreteParametric &inter_event , const Distribution &time)

{
  int i , j;
  double norm , sum , *pmass , *icumul , *scumul , *tmass , *tcumul;


  offset = 0;
  nb_value = MIN(inter_event.nb_value - 1 , time.nb_value);

  pmass = mass;
  icumul = inter_event.cumul;

  // computation of the normalization quantity

  if (ident == CATEGORICAL) {
    norm = inter_event.mean;
  }
  else {
    norm = parametric_mean_computation();
  }

  // computation of the probability mass function

  tmass = time.mass;
  tcumul = time.cumul;

  for (i = 0;i < nb_value;i++) {
    *pmass = (1. - *tcumul) * (1. - *icumul) / norm;

    if (*tmass > 0.) {
      scumul = icumul;
      sum = 0.;
      for (j = i;j < inter_event.nb_value - 1;j++) {
        sum += (1. - *scumul++);
      }

      *pmass += *tmass * sum / norm;
    }

    pmass++;
    icumul++;
    tmass++;
    tcumul++;
  }

  cumul_computation();

  max_computation();
  mean_computation();
  variance_computation();
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the number of events distribution from
 *         the inter-event distribution.
 *
 *  \param[in] inter_event reference on the inter-event distribution.
 */
/*--------------------------------------------------------------*/

void NbEvent::computation(DiscreteParametric &inter_event)

{
  int i , j;
  int time_nb_value;
  double bcumul = 1. , previous_cumul = 1. , *pmass , *pcumul;
  DiscreteParametric *nevent_time;
  Forward *forward;


  // computation of the number of values

  switch (type) {
  case ORDINARY :
    nb_value = time / inter_event.offset + 1;
    break;
  case EQUILIBRIUM :
    nb_value = (time - 1) / inter_event.offset + 2;
    break;
  }

  switch (type) {
  case ORDINARY :
    time_nb_value = time + 1;
    break;
  case EQUILIBRIUM :
    time_nb_value = time;
    break;
  }

  nevent_time = new DiscreteParametric(time + 1 , ident);

  pmass = mass - 1;
  pcumul = cumul - 1;

  for (i = 0;i < nb_value;i++) {
    if (i < nb_value - 1) {

      if (i == 0) {
        switch (type) {

        // ordinary renewal process: time to the 1st event distributed according to
        // the inter-event distribution

        case ORDINARY : {
          nevent_time->mass_copy(inter_event , time + 1);
          break;
        }

        // equilibrium renewal process: time to the 1st event distributed according to
        // the forward recurrence time distribution

        case EQUILIBRIUM : {
          forward = new Forward(inter_event);
          nevent_time->mass_copy(*forward , time + 1);
          break;
        }
        }

        nevent_time->cumul_computation();
      }

      else {
        switch (type) {
        case ORDINARY :
          j = i + 1;
          break;
        case EQUILIBRIUM :
          j = i;
          break;
        }

        // computation of the time to the (n+1)th events distribution

        if ((j == 1) && (ident != CATEGORICAL)) {
          nevent_time->mass_copy(inter_event , time + 1);
          nevent_time->cumul_computation();
        }

        else {
          switch (ident) {

          case CATEGORICAL : {
            nevent_time->convolution(inter_event , *nevent_time , time + 1);
            nevent_time->cumul_computation();
            break;
          }

          case BINOMIAL : {
            nevent_time->inf_bound = j * inter_event.inf_bound;
            nevent_time->sup_bound = j * inter_event.sup_bound;
            nevent_time->probability = inter_event.probability;
            nevent_time->binomial_computation(time_nb_value , RENEWAL);
            break;
          }

          case POISSON : {
            nevent_time->inf_bound = j * inter_event.inf_bound;
            nevent_time->parameter = j * inter_event.parameter;
            nevent_time->poisson_computation(time_nb_value , RENEWAL_THRESHOLD , RENEWAL);
            break;
          }

          case NEGATIVE_BINOMIAL : {
            nevent_time->inf_bound = j * inter_event.inf_bound;
            nevent_time->parameter = j * inter_event.parameter;
            nevent_time->probability = inter_event.probability;
            nevent_time->negative_binomial_computation(time_nb_value , RENEWAL_THRESHOLD , RENEWAL);
            break;
          }
          }
        }

        if ((type == EQUILIBRIUM) && (ident != CATEGORICAL)) {
          nevent_time->convolution(*forward , *nevent_time , time + 1);
          nevent_time->cumul_computation();
        }
      }

      if (time < nevent_time->nb_value) {
        bcumul = MIN(nevent_time->cumul[time] , 1.);
      }
      *++pmass = previous_cumul - bcumul;
      *++pcumul = 1. - bcumul;
      previous_cumul = bcumul;
    }

    else {
      *++pmass = previous_cumul;
      *++pcumul = 1.;
    }
  }

  delete nevent_time;
  if (type == EQUILIBRIUM) {
    delete forward;
  }

  offset_computation();
  max_computation();
  mean_computation();
  variance_computation();
}


/*--------------------------------------------------------------*/
/**
 *  \brief Fast computation (O(n)) in the case of an ordinary renewal process of
 *         the number of events distribution from a binomial inter-event distribution.
 */
/*--------------------------------------------------------------*/

void NbEvent::binomial_computation()

{
  int i , j;
  int main_set , main_subset , inf_bound_set , inf_bound_subset ,
      rapid_index , nb_term;
  double failure = 1. - probability , success = probability , k_success , main_term ,
         inf_bound_term , sum , scale , *pmass , *pcumul;


  nb_value = time / inf_bound + 1;

  pmass = mass;
  pcumul = cumul;

  // computation of the success probability at the power (upper bound - lower bound)

  k_success = 1.;
  for (j = 0;j < sup_bound - inf_bound;j++) {
    k_success *= success;
  }
  main_term = k_success;

  // values of null probability such that (observation period + 1) > (i + 1) * (upper bound)

  i = 0;
  while (time + 1 > (i + 1) * sup_bound) {
    main_term *= k_success;
    *pmass++ = 0.;
    *pcumul++ = 0.;
    i++;
  }

  // computation of the first non-null probability value (exhaustive computation)

  main_subset = (i + 1) * (sup_bound - inf_bound);
  main_set = main_subset;
  sum = main_term;

  nb_term = (i + 1) * sup_bound - time - 1;
  if ((i == nb_value - 1) && ((i + 1) * inf_bound > time + 1)) {
    nb_term = (i + 1) * (sup_bound - inf_bound);
  }

  for (j = 0;j < nb_term;j++) {
    scale = (double)main_subset / (double)(main_set - (main_subset - 1));
    main_subset--;
    main_term *= scale * failure / success;
    sum += main_term;
  }
  *pmass = sum;
  *pcumul = sum;
  i++;
  rapid_index = i;

  // fast computation of the probability masses for the following values

  while (i < nb_value) {

    // computation of the main terms

    // computation of the 1st term

    if (i > rapid_index) {
      if (inf_bound == 0) {
        main_set++;
        scale = (double)main_set / (double)(main_set - main_subset);
        main_term *= scale * failure;
      }
      else {
        main_subset = inf_bound_subset;
        main_set = inf_bound_set;
        main_term = inf_bound_term;
      }
    }
    if ((i == rapid_index) || (inf_bound > 0)) {
      scale = (double)main_subset / (double)(main_set - (main_subset - 1));
      main_subset--;
      main_term *= scale * failure;
    }
    *++pmass = main_term;

    // computation of the (j - lower bound - 1) following terms

    for (j = 1;j <= sup_bound - inf_bound - 1;j++) {
      main_set++;
      scale = (double)main_set / (double)(main_set - main_subset);
      main_term *= scale * failure;
      *pmass += main_term;
    }

    // computation of the terms corresponding to the lower bound

    // computation of the 1st term

    if (inf_bound > 0) {
      inf_bound_subset = main_subset;
      inf_bound_set = main_set + 1;
      scale = (double)inf_bound_set / (double)(inf_bound_set - inf_bound_subset);
      inf_bound_term = main_term * scale * failure / success;
      *pmass += inf_bound_term;
    }

    // computation of the (lower bound - 1) following terms

    if (inf_bound > 1) {
      nb_term = inf_bound - 1;
      if ((i == nb_value - 1) && ((i + 1) * inf_bound > time + 1)) {
        nb_term = time - i * inf_bound;
      }
      for (j = 0;j < nb_term;j++) {
        scale = (double)inf_bound_subset / (double)(inf_bound_set - (inf_bound_subset - 1));
        inf_bound_subset--;
        inf_bound_term *= scale * failure / success;
        *pmass += inf_bound_term;
      }
    }

    // update of the cumulative distribution function

    pcumul++;
    *pcumul = *(pcumul - 1) + *pmass;
    i++;
  }

  offset_computation();
  max_computation();
  mean_computation();
  variance_computation();
}


/*--------------------------------------------------------------*/
/**
 *  \brief Fast computation (O(n)) in the case of an ordinary renewal process of
 *         the number of events distribution from a negative binomial inter-event distribution.
 *         This computation requires an integer-valued negative binomial shape parameter.
 */
/*--------------------------------------------------------------*/

void NbEvent::negative_binomial_computation()

{
  int i , j;
  int main_set , main_subset , inf_bound_set , inf_bound_subset;
  double failure = 1. - probability , success = probability , main_term , inf_bound_term ,
         scale , *pmass;


  nb_value = time / inf_bound + 1;

  pmass = mass;

  // computation of the terms for number of events = 0

  main_term = 1.;
  for (i = 0;i < parameter;i++) {
    main_term *= success;
  }
  *pmass = main_term;

  main_subset = (int)parameter - 1;
  main_set = main_subset;
  for (i = inf_bound + 1;i <= time;i++) {
    main_set++;
    scale = (double)main_set / (double)(main_set - main_subset);
    main_term *= scale * failure;
    *pmass += main_term;
  }
  *pmass = 1. - *pmass;
  pmass++;

  // computation of probability masses for values from 1 to (observation period / lower bound)

  for (i = 1;i < nb_value;i++) {
    *pmass = 0.;

    // computation of the terms corresponding to the lower bound

    // computation of the 1st term

    if (inf_bound > 1) {
      if (i == 1) {
        inf_bound_subset = main_subset;
        inf_bound_set = main_set;
        inf_bound_term = main_term;
      }
      else {
        scale = (double)(main_set - main_subset) / (double)main_set;
        inf_bound_subset = main_subset;
        inf_bound_set = main_set - 1;
        inf_bound_term = main_term * scale * success / failure;
      }
      *pmass += inf_bound_term;

      // computation of the (lower bound - 2) following terms

      if (inf_bound > 2) {
        for (j = 0;j < inf_bound - 2;j++) {
          scale = (double)(inf_bound_set - inf_bound_subset) / (double)inf_bound_set;
          inf_bound_set--;
          inf_bound_term *= scale / failure;
          *pmass += inf_bound_term;
        }
      }
    }

    // computation of the main terms

    // computation of the 1st term

    if (i == 1) {
      main_subset++;
      main_set++;
      scale = (double)main_set / (double)main_subset;
      main_term *= scale * failure;
    }
    else {
      scale = (double)(main_set - main_subset) / (double)(main_subset + 1);
      main_subset++;
      main_term *= scale * success;
    }
    if (inf_bound > 1) {
      for (j = 0;j < inf_bound - 1;j++) {
        scale = (double)(main_set - main_subset) / (double)main_set;
        main_set--;
        main_term *= scale / failure;
      }
    }
    if (parameter == 1) {
      main_term /= failure;
    }
    *pmass += main_term;

    // computation of the (k - 1) following terms

    if (parameter > 1) {
      for (j = 1;j <= (int)parameter - 1;j++) {
        main_subset++;
        main_set++;
        scale = (double)main_set / (double)main_subset;
        main_term *= scale * success;
        if (j == parameter - 1) {
          main_term /= failure;
        }
        *pmass += main_term;
      }
    }

    pmass++;
  }

  offset_computation();
  cumul_computation();

  max_computation();
  mean_computation();
  variance_computation();
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation in the case of an ordinary renewal process of
 *         the number of events distribution from a parametric inter-event
 *         distribution (binomial, Poisson, negative binomial).
 *
 *  \param[in] inter_event reference on an inter-event distribution.
 */
/*--------------------------------------------------------------*/

void NbEvent::ordinary_computation(DiscreteParametric &inter_event)

{
  if (type == ORDINARY) {
    if ((ident == BINOMIAL) &&
        (time * (1. - probability) / (probability * sqrt((double)(MAX(1 , inf_bound)))) < RB_THRESHOLD)) {
      binomial_computation();
    }

    else if ((ident == NEGATIVE_BINOMIAL) && (parameter == (int)parameter) && (inf_bound > 0) &&
             (time - (time / inf_bound + 1) * inf_bound + 1 >= 0) &&
             (time * probability / ((1. - probability) * sqrt((double)(MAX(1 , inf_bound)))) < RNB_THRESHOLD)) {

#     ifdef DEBUG
      cout << "algebric calculation" << endl;
#     endif

      negative_binomial_computation();
    }

    else {
      computation(inter_event);
    }
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of non-event/event probabilities as a function of time.
 */
/*--------------------------------------------------------------*/

void Renewal::index_event_computation()

{
  int i , j;
  double no_event , event , *pindex_event;


  switch (type) {

  case ORDINARY : {
    index_event->offset = 0;

    index_event->point[0][0] = 0.;
    index_event->point[1][0] = 1.;
    break;
  }

  case EQUILIBRIUM : {
    index_event->offset = 1;

    index_event->point[0][0] = 0.;
    index_event->point[1][0] = 0.;
    break;
  }
  }

  for (i = 1;i < index_event->length;i++) {
    pindex_event = index_event->point[1] + i;
    no_event = 0.;
    event = 0.;

    for (j = 1;j <= MIN(i , inter_event->nb_value - 1);j++) {
      if (j < i) {
        no_event += (1. - inter_event->cumul[j]) * *--pindex_event;
        event += inter_event->mass[j] * *pindex_event;
      }

      else {
        switch (type) {
        case ORDINARY :
          no_event += 1. - inter_event->cumul[j];
          event += inter_event->mass[j];
          break;
        case EQUILIBRIUM :
          no_event += (1. - forward->cumul[j]);
          event += forward->mass[j];
          break;
        }
      }
    }

    index_event->point[0][i] = no_event;
    index_event->point[1][i] = event;
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the distributions of a renewal process from
 *         the inter-event distribution.
 *         time interval: length-biased distribution,
 *                        backward and forward recurrence time distributions,
 *        count: number of events distributions and resulting mixture,
 *        intensity: no-event/event probabilities as a function of time.
 *
 *  \param[in] inter_event_flag flag for the computation of the inter-event distribution,
 *  \param[in] itype            renewal process type (ORDINARY/EQUILIBRIUM),
 *  \param[in] dtime            pointer on the observation period distribution.
 */
/*--------------------------------------------------------------*/

void Renewal::computation(bool inter_event_flag , process_type itype , const Distribution *dtime)

{
  int i , j;
  int nb_value , time_nb_value;
  double sum , *tmass , *pmass , *cumul , *previous_cumul , *pcumul1 , *pcumul2;
  DiscreteParametric *pnevent_time , *power , *forward_power;


  if (itype == DEFAULT_TYPE) {
    itype = type;
  }

  // computation of the inter-event distribution, the length-biased distribution and
  // the forward recurrence time distribution

  if (inter_event_flag) {
    inter_event->computation(1 , RENEWAL_THRESHOLD);
    length_bias->computation(*inter_event);
    forward->computation(*inter_event);
  }

  if ((itype != type) || ((dtime) && (*dtime != *time))) {
    type_init(itype);

    tmass = time->mass + time->offset;
    for (i = time->offset;i < time->nb_value;i++) {
      if (*tmass++ > 0.) {
        delete nb_event[i];
      }
    }

    delete mixture;

    if ((dtime) && (*dtime != *time)) {
      delete time;

      for (i = 1;i <= nb_event_max;i++) {
        delete nevent_time[i];
      }
      delete [] nevent_time;

      delete [] nb_event;

      delete index_event;

      time = new Distribution(*dtime);

      switch (type) {
      case ORDINARY :
        nb_event_max = (time->nb_value - 1) / inter_event->offset;
        break;
      case EQUILIBRIUM :
        nb_event_max = (time->nb_value - 2) / inter_event->offset + 1;
        break;
      }

      nevent_time = new DiscreteParametric*[nb_event_max + 1];
      for (i = 0;i <= nb_event_max;i++) {
        nevent_time[i] = NULL;
      }

      nb_event = new NbEvent*[time->nb_value];

      index_event = new Curves(2 , time->nb_value);
    }

    for (i = 0;i < time->offset;i++) {
      nb_event[i] = NULL;
    }

    tmass = time->mass + time->offset;
    for (i = time->offset;i < time->nb_value;i++) {
      if (*tmass++ > 0.) {
        switch (type) {
        case ORDINARY :
          nb_value = i / inter_event->offset + 1;
          break;
        case EQUILIBRIUM :
          nb_value = (i - 1) / inter_event->offset + 2;
          break;
        }

        nb_event[i] = new NbEvent(type , i , nb_value , inter_event->ident);
      }

      else {
        nb_event[i] = NULL;
      }
    }

    init(inter_event->inf_bound , inter_event->sup_bound ,
         inter_event->parameter , inter_event->probability);

    mixture = new Distribution(nb_value);
  }

  // computation of the backward recurrence time distribution

  backward->computation(*inter_event , *time);

  // construction and initialization of the cumulative distribution functions

  cumul = new double[time->nb_value];
  previous_cumul = new double[time->nb_value];

  tmass = time->mass + time->offset;
  pcumul1 = cumul + time->offset;
  pcumul2 = previous_cumul + time->offset;

  for (i = time->offset;i < time->nb_value;i++) {
    if (*tmass++ > 0.) {
      *pcumul1 = 1.;
      *pcumul2 = 1.;
    }
    pcumul1++;
    pcumul2++;
  }

  // computation of the number of values of the number of events distributions and the resulting mixture

  tmass = time->mass + time->offset;

  for (i = time->offset;i < time->nb_value;i++) {
    if (*tmass++ > 0.) {
      switch (type) {
      case ORDINARY :
        nb_event[i]->nb_value = i / inter_event->offset + 1;
        break;
      case EQUILIBRIUM :
        nb_event[i]->nb_value = (i - 1) / inter_event->offset + 2;
        break;
      }
    }
  }

  mixture->nb_value = nb_event[time->nb_value - 1]->nb_value;

  switch (type) {
  case ORDINARY :
    time_nb_value = time->nb_value;
    break;
  case EQUILIBRIUM :
    time_nb_value = time->nb_value - 1;
    break;
  }

  if (type == EQUILIBRIUM) {
    pnevent_time = new DiscreteParametric(time->nb_value , inter_event->ident);
    power = pnevent_time;
  }

  // computation of the number of events distributions and the ressulting mixture

  pmass = mixture->mass;

  for (i = 0;i < mixture->nb_value;i++) {
    if (i < mixture->nb_value - 1) {
      j = i + 1;

      if (!nevent_time[j]) {
        nevent_time[j] = new DiscreteParametric(time->nb_value , inter_event->ident);
      }

      switch (type) {
      case ORDINARY :
        power = nevent_time[j];
        break;
      case EQUILIBRIUM :
        forward_power = nevent_time[j];
        break;
      }

      if (i == 0) {
        if (type == EQUILIBRIUM) {
          forward_power->mass_copy(*forward , time->nb_value);
          forward_power->cumul_computation();
        }

        power->mass_copy(*inter_event , time->nb_value);
        power->cumul_computation();
      }

      else {
        if (type == EQUILIBRIUM) {
          forward_power->convolution(*forward , *power , time->nb_value);
          forward_power->cumul_computation();
        }

        // computation of the time to the (n+1)th event distribution

        power->ident = inter_event->ident;

        switch (inter_event->ident) {

        case CATEGORICAL : {
          switch (type) {
          case ORDINARY :
            power->convolution(*inter_event , *nevent_time[j - 1] , time_nb_value);
            break;
          case EQUILIBRIUM :
            power->convolution(*inter_event , *pnevent_time , time_nb_value);
            break;
          }

          power->cumul_computation();
          break;
        }

        case BINOMIAL : {
          power->inf_bound = j * inter_event->inf_bound;
          power->sup_bound = j * inter_event->sup_bound;
          power->probability = inter_event->probability;
          power->binomial_computation(time_nb_value , RENEWAL);
          break;
        }

        case POISSON : {
          power->inf_bound = j * inter_event->inf_bound;
          power->parameter = j * inter_event->parameter;
          power->poisson_computation(time_nb_value , RENEWAL_THRESHOLD , RENEWAL);
          break;
        }

        case NEGATIVE_BINOMIAL : {
          power->inf_bound = j * inter_event->inf_bound;
          power->parameter = j * inter_event->parameter;
          power->probability = inter_event->probability;
          power->negative_binomial_computation(time_nb_value , RENEWAL_THRESHOLD , RENEWAL);
          break;
        }
        }
      }
    }

    // computation of the number of events distributions and the resulting mixture

    tmass = time->mass + time->offset;
    pcumul1 = cumul + time->offset;
    pcumul2 = previous_cumul + time->offset;
    sum = 0.;

    for (j = time->offset;j < time->nb_value;j++) {
      if (*tmass > 0.) {
        if (i == nb_event[j]->nb_value - 1) {
          *pcumul1 = 0.;
        }

        else {
          switch (type) {

          case ORDINARY : {
            if (j < power->nb_value) {
              *pcumul1 = MIN(power->cumul[j] , 1.);
            }
            break;
          }

          case EQUILIBRIUM : {
            if (j < forward_power->nb_value) {
              *pcumul1 = MIN(forward_power->cumul[j] , 1.);
            }
            break;
          }
          }
        }

        if (i < nb_event[j]->nb_value) {
          nb_event[j]->mass[i] = *pcumul2 - *pcumul1;
          nb_event[j]->cumul[i] = 1. - *pcumul1;
          *pcumul2 = *pcumul1;
          sum += *tmass * nb_event[j]->mass[i];
        }
      }

      tmass++;
      pcumul1++;
      pcumul2++;
    }

    *pmass++ = sum;
  }

  for (i = 1;i < mixture->nb_value;i++) {
    nevent_time[i]->max_computation();
  }

  tmass = time->mass + time->offset;

  for (i = time->offset;i < time->nb_value;i++) {
    if (*tmass++ > 0.) {
      nb_event[i]->offset_computation();
      nb_event[i]->max_computation();
      nb_event[i]->mean_computation();
      nb_event[i]->variance_computation();
    }
  }

  mixture->offset_computation();
  mixture->cumul_computation();

  mixture->max_computation();
  mixture->mean_computation();
  mixture->variance_computation();

  if (type == EQUILIBRIUM) {
    delete pnevent_time;
  }

  delete [] cumul;
  delete [] previous_cumul;

  index_event_computation();
}


};  // namespace sequence_analysis
