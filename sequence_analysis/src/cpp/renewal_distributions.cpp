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
 *       $Id: renewal_distributions.cpp 18066 2015-04-23 10:51:02Z guedon $
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



#include "stat_tool/stat_tools.h"
#include "stat_tool/curves.h"
#include "stat_tool/distribution.h"
#include "stat_tool/markovian.h"
#include "stat_tool/vectors.h"
#include "stat_tool/distance_matrix.h"

#include "renewal.h"

using namespace std;
using namespace stat_tool;


namespace sequence_analysis {



/*--------------------------------------------------------------*
 *
 *  Calcul de la loi biaisee par la longueur a partir de la loi inter-evenement.
 *
 *  argument : loi inter-evenement.
 *
 *--------------------------------------------------------------*/

void LengthBias::computation(const DiscreteParametric &inter_event)

{
  register int i;
  double norm , *pmass , *imass;


  offset = inter_event.offset;
  nb_value = inter_event.nb_value;

  pmass = mass;
  for (i = 0;i < offset;i++) {
    *pmass++ = 0.;
  }
  imass = inter_event.mass + offset;

  // calcul de la quantite de normalisation

  if (ident == CATEGORICAL) {
    norm = inter_event.mean;
  }
  else {
    norm = parametric_mean_computation();
  }

  // calcul des probabilites des valeurs

  for (i = offset;i < nb_value;i++) {
    *pmass++ = i * *imass++ / norm;
  }

  cumul_computation();

  max_computation();
  mean_computation();
  variance_computation();
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la loi de l'intervalle de temps apres le dernier evenement
 *  a partir de la loi inter-evenement.
 *
 *  arguments : loi inter-evenement, loi du temps d'observation.
 *
 *--------------------------------------------------------------*/

void Backward::computation(const DiscreteParametric &inter_event , const Distribution &time)

{
  register int i , j;
  double norm , sum , *pmass , *icumul , *scumul , *tmass , *tcumul;


  offset = 0;
  nb_value = MIN(inter_event.nb_value - 1 , time.nb_value);

  pmass = mass;
  icumul = inter_event.cumul;

  // calcul de la quantite de normalisation

  if (ident == CATEGORICAL) {
    norm = inter_event.mean;
  }
  else {
    norm = parametric_mean_computation();
  }

  // calcul des probabilites des valeurs

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


/*--------------------------------------------------------------*
 *
 *  Calcul de la loi du nombre d'evenements correspondant a
 *  une loi quelconque comme loi inter-evenement.
 *
 *  argument : reference sur la loi inter-evenement.
 *
 *--------------------------------------------------------------*/

void NbEvent::computation(DiscreteParametric &inter_event)

{
  register int i , j;
  int time_nb_value;
  double bcumul = 1. , previous_cumul = 1. , *pmass , *pcumul;
  DiscreteParametric *nevent_time;
  Forward *forward;


  // calcul du nombre de valeurs

  switch (type) {
  case 'o' :
    nb_value = time / inter_event.offset + 1;
    break;
  case 'e' :
    nb_value = (time - 1) / inter_event.offset + 2;
    break;
  }

  switch (type) {
  case 'o' :
    time_nb_value = time + 1;
    break;
  case 'e' :
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

        // processus de renouvellement ordinaire :
        // temps avant le 1er evenement distribue selon la loi inter-evenement

        case 'o' : {
          nevent_time->mass_copy(inter_event , time + 1);
          break;
        }

        // processus de renouvellement en equilibre :
        // temps avant le 1er evenement distribue selon la loi "forward"

        case 'e' : {
          forward = new Forward(inter_event);
          nevent_time->mass_copy(*forward , time + 1);
          break;
        }
        }

        nevent_time->cumul_computation();
      }

      else {
        switch (type) {
        case 'o' :
          j = i + 1;
          break;
        case 'e' :
          j = i;
          break;
        }

        // calcul de la loi du temps de (n+1) evenements

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
            nevent_time->binomial_computation(time_nb_value , 'r');
            break;
          }

          case POISSON : {
            nevent_time->inf_bound = j * inter_event.inf_bound;
            nevent_time->parameter = j * inter_event.parameter;
            nevent_time->poisson_computation(time_nb_value , RENEWAL_THRESHOLD , 'r');
            break;
          }

          case NEGATIVE_BINOMIAL : {
            nevent_time->inf_bound = j * inter_event.inf_bound;
            nevent_time->parameter = j * inter_event.parameter;
            nevent_time->probability = inter_event.probability;
            nevent_time->negative_binomial_computation(time_nb_value ,
                                                       RENEWAL_THRESHOLD , 'r');
            break;
          }
          }
        }

        if ((type == 'e') && (ident != CATEGORICAL)) {
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
  if (type == 'e') {
    delete forward;
  }

  offset_computation();
  max_computation();
  mean_computation();
  variance_computation();
}


/*--------------------------------------------------------------*
 *
 *  Calcul rapide (O(n)), dans le cas d'un processus de renouvellement
 *  ordinaire, de la loi du nombre d'evenements correspondant a
 *  une loi binomiale comme loi inter-evenement.
 *
 *--------------------------------------------------------------*/

void NbEvent::binomial_computation()

{
  register int i , j;
  int main_set , main_subset , inf_bound_set , inf_bound_subset ,
      rapid_index , nb_term;
  double failure = 1. - probability , success = probability , k_success , main_term ,
         inf_bound_term , sum , scale , *pmass , *pcumul;


  nb_value = time / inf_bound + 1;

  pmass = mass;
  pcumul = cumul;

  // calcul de la probabilite de succes a la puissance de
  // (borne superieure - borne inferieure)

  k_success = 1.;
  for (j = 0;j < sup_bound - inf_bound;j++) {
    k_success *= success;
  }
  main_term = k_success;

  // valeurs telles que (temps d'observation + 1) > (i + 1) * (borne superieure)
  // de probabilite nulle

  i = 0;
  while (time + 1 > (i + 1) * sup_bound) {
    main_term *= k_success;
    *pmass++ = 0.;
    *pcumul++ = 0.;
    i++;
  }

  // calcul de la premiere valeur de probabilite > 0. (calcul exhaustif)

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

  // calcul rapide des probabilites des valeurs suivantes

  while (i < nb_value) {

    // calcul des termes principaux

    // calcul du 1er terme

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

    // calcul des (j - borne inferieure - 1) termes suivants

    for (j = 1;j <= sup_bound - inf_bound - 1;j++) {
      main_set++;
      scale = (double)main_set / (double)(main_set - main_subset);
      main_term *= scale * failure;
      *pmass += main_term;
    }

    // calcul des termes correspondant a la borne inferieure

    // calcul du 1er terme

    if (inf_bound > 0) {
      inf_bound_subset = main_subset;
      inf_bound_set = main_set + 1;
      scale = (double)inf_bound_set / (double)(inf_bound_set - inf_bound_subset);
      inf_bound_term = main_term * scale * failure / success;
      *pmass += inf_bound_term;
    }

    // calcul des (borne inferieure - 1) termes suivants

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

    // mise a jour le la fonction de repartition

    pcumul++;
    *pcumul = *(pcumul - 1) + *pmass;
    i++;
  }

  offset_computation();
  max_computation();
  mean_computation();
  variance_computation();
}


/*--------------------------------------------------------------*
 *
 *  Calcul rapide (O(n)), dans le cas d'un processus de renouvellement
 *  ordinaire, de la loi du nombre  d'evenements correspondant
 *  a une loi binomiale negative comme loi inter-evenement. Ce calcul
 *  n'est valable que si le parametre est entier.
 *
 *--------------------------------------------------------------*/

void NbEvent::negative_binomial_computation()

{
  register int i , j;
  int main_set , main_subset , inf_bound_set , inf_bound_subset;
  double failure = 1. - probability , success = probability , main_term , inf_bound_term ,
         scale , *pmass;


  nb_value = time / inf_bound + 1;

  pmass = mass;

  // calcul des termes pour nombre d'evenements = 0

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

  // calcul des probabilites des valeurs de 1 a (temps d'observation / borne inferieure)

  for (i = 1;i < nb_value;i++) {
    *pmass = 0.;

    // calcul des termes correspondant a la borne inferieure

    // calcul du 1er terme

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

      // calcul des (borne inferieure - 2) termes suivants

      if (inf_bound > 2) {
        for (j = 0;j < inf_bound - 2;j++) {
          scale = (double)(inf_bound_set - inf_bound_subset) / (double)inf_bound_set;
          inf_bound_set--;
          inf_bound_term *= scale / failure;
          *pmass += inf_bound_term;
        }
      }
    }

    // calcul des termes principaux

    // calcul du 1er terme

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

    // calcul des (k-1) termes suivant

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


/*--------------------------------------------------------------*
 *
 *  Calcul dans le cas dans le cas d'un processus de renouvellement
 *  ordinaire de la loi du nombre d'evenements correspondant a une loi
 *  inter-evenement parametrique (binomiale negative, binomiale, Poisson).
 *
 *  argument : reference sur la loi inter-evenement.
 *
 *--------------------------------------------------------------*/

void NbEvent::ordinary_computation(DiscreteParametric &inter_event)

{
  if (type == 'o') {
    if ((ident == BINOMIAL) &&
        (time * (1. - probability) / (probability * sqrt((double)(MAX(1 , inf_bound)))) < RB_THRESHOLD)) {
      binomial_computation();
    }

    else if ((ident == NEGATIVE_BINOMIAL) && (parameter == (int)parameter) && (inf_bound > 0) &&
             (time - (time / inf_bound + 1) * inf_bound + 1 >= 0) &&
             (time * probability / ((1. - probability) * sqrt((double)(MAX(1 , inf_bound)))) < RNB_THRESHOLD)) {

#     ifdef DEBUG
      cout << "calcul algebrique" << endl;
#     endif

      negative_binomial_computation();
    }

    else {
      computation(inter_event);
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Calcul des probabilites de non-evenement/evenement fonction du temps.
 *
 *--------------------------------------------------------------*/

void Renewal::index_event_computation()

{
  register int i , j;
  double no_event , event , *pindex_event;


  switch (type) {

  case 'o' : {
    index_event->offset = 0;

    index_event->point[0][0] = 0.;
    index_event->point[1][0] = 1.;
    break;
  }

  case 'e' : {
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
        case 'o' :
          no_event += 1. - inter_event->cumul[j];
          event += inter_event->mass[j];
          break;
        case 'e' :
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


/*--------------------------------------------------------------*
 *
 *  Calcul des lois d'un processus de renouvellement a partir
 *  d'une loi quelconque comme loi inter-evenement.
 *  intervalle : loi inter-evenement, loi biaisee par la longueur,
 *               loi de l'intervalle de temps apres le dernier evenement,
 *               loi de l'intervalle de temps residuel,
 *  comptage : melange de lois du nombre d'evenements,
 *  intensite.
 *
 *  arguments : flag pour le calcul de la loi inter-evenement,
 *              type de processus ('o' : ordinaire, 'e' : en equilibre),
 *              pointeur sur la loi du temps d'observation.
 *
 *--------------------------------------------------------------*/

void Renewal::computation(bool inter_event_flag , char itype , const Distribution *dtime)

{
  register int i , j;
  int nb_value , time_nb_value;
  double sum , *tmass , *pmass , *cumul , *previous_cumul , *pcumul1 , *pcumul2;
  DiscreteParametric *pnevent_time , *power , *forward_power;


  if (itype == 'v') {
    itype = type;
  }

  // calcul de la loi inter-evenement, de la loi biaisee par la longueur et
  // de la loi de l'intervalle de temps residuel

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
      case 'o' :
        nb_event_max = (time->nb_value - 1) / inter_event->offset;
        break;
      case 'e' :
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
        case 'o' :
          nb_value = i / inter_event->offset + 1;
          break;
        case 'e' :
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

  // calcul de la loi du temps apres le dernier evenement

  backward->computation(*inter_event , *time);

  // creation et initialisation des variables "fonction de repartition"

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

  // calcul du nombre de valeurs des lois du nombre d'evenements et de la loi resultante

  tmass = time->mass + time->offset;

  for (i = time->offset;i < time->nb_value;i++) {
    if (*tmass++ > 0.) {
      switch (type) {
      case 'o' :
        nb_event[i]->nb_value = i / inter_event->offset + 1;
        break;
      case 'e' :
        nb_event[i]->nb_value = (i - 1) / inter_event->offset + 2;
        break;
      }
    }
  }

  mixture->nb_value = nb_event[time->nb_value - 1]->nb_value;

  switch (type) {
  case 'o' :
    time_nb_value = time->nb_value;
    break;
  case 'e' :
    time_nb_value = time->nb_value - 1;
    break;
  }

  if (type == 'e') {
    pnevent_time = new DiscreteParametric(time->nb_value , inter_event->ident);
    power = pnevent_time;
  }

  // calcul des lois du nombre d'evenements et du melange resultant


  pmass = mixture->mass;

  for (i = 0;i < mixture->nb_value;i++) {
    if (i < mixture->nb_value - 1) {
      j = i + 1;

      if (!nevent_time[j]) {
        nevent_time[j] = new DiscreteParametric(time->nb_value , inter_event->ident);
      }

      switch (type) {
      case 'o' :
        power = nevent_time[j];
        break;
      case 'e' :
        forward_power = nevent_time[j];
        break;
      }

      if (i == 0) {
        if (type == 'e') {
          forward_power->mass_copy(*forward , time->nb_value);
          forward_power->cumul_computation();
        }

        power->mass_copy(*inter_event , time->nb_value);
        power->cumul_computation();
      }

      else {
        if (type == 'e') {
          forward_power->convolution(*forward , *power , time->nb_value);
          forward_power->cumul_computation();
        }

        // calcul de la loi du temps de (n+1) evenements

        power->ident = inter_event->ident;

        switch (inter_event->ident) {

        case CATEGORICAL : {
          switch (type) {
          case 'o' :
            power->convolution(*inter_event , *nevent_time[j - 1] , time_nb_value);
            break;
          case 'e' :
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
          power->binomial_computation(time_nb_value , 'r');
          break;
        }

        case POISSON : {
          power->inf_bound = j * inter_event->inf_bound;
          power->parameter = j * inter_event->parameter;
          power->poisson_computation(time_nb_value , RENEWAL_THRESHOLD , 'r');
          break;
        }

        case NEGATIVE_BINOMIAL : {
          power->inf_bound = j * inter_event->inf_bound;
          power->parameter = j * inter_event->parameter;
          power->probability = inter_event->probability;
          power->negative_binomial_computation(time_nb_value ,
                                               RENEWAL_THRESHOLD , 'r');
          break;
        }
        }
      }
    }

    // calcul des lois du nombre d'evenements et de la loi resultante

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

          case 'o' : {
            if (j < power->nb_value) {
              *pcumul1 = MIN(power->cumul[j] , 1.);
            }
            break;
          }

          case 'e' : {
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

  if (type == 'e') {
    delete pnevent_time;
  }

  delete [] cumul;
  delete [] previous_cumul;

  index_event_computation();
}


};  // namespace sequence_analysis
