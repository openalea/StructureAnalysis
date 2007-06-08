/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       AMAPmod: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2002 UMR Cirad/Inra Modelisation des Plantes
 *
 *       File author(s): Y. Guedon (yann.guedon@cirad.fr)
 *
 *       $Source$
 *       $Id$
 *
 *       Forum for AMAPmod developers: amldevlp@cirad.fr
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
#include "stat_tool/regression.h"
#include "stat_tool/curves.h"
#include "stat_tool/markovian.h"
#include "sequences.h"
#include "markov.h"

#ifdef DEBUG
#include "sequence_label.h"
#endif

using namespace std;



/*--------------------------------------------------------------*
 *
 *  Calcul des probabilites de chaque etat en fonction de l'index
 *  pour une chaine de Markov.
 *
 *--------------------------------------------------------------*/

void Markov::index_state_distribution()

{
  register int i , j , k;
  int power[ORDER] , state_index[ORDER];
  double sum , *state_seq , *pstate_seq , *states , *pstates , **ptransition;
  Curves *index_state;


  index_state = process[0]->index_value;

  i = 1;
  for (j = 0;j < order;j++) {
    power[j] = i;
    i *= nb_state;
  }

  // initialisation des probabilites des sequences d'etats

  state_seq = new double[nb_row];
  pstate_seq = new double[nb_row];

  pstates = pstate_seq;
  i = 0;

  for (j = 0;j < nb_row;j++) {
    if (j == self_row[i]) {
      index_state->point[i][0] = initial[i];
      *pstates++ = initial[i++];
    }
    else {
      *pstates++ = 0.;
    }
  }

  // calcul des probabilites de chaque etat en fonction de l'index

  for (i = 1;i < index_state->length;i++) {

    // calcul des probabilites des sequences d'etats de
    // longueur = ordre de la chaine de Markov

    for (j = 0;j < order;j++) {
      state_index[j] = 0;
    }

    states = state_seq;

    for (j = 0;j < nb_row;j++) {
      ptransition = transition;
      pstates = pstate_seq;
      for (k = 0;k < order - 1;k++) {
        ptransition += state_index[k] * power[k + 1];
        pstates += state_index[k] * power[k + 1];
      }

      *states = 0.;
      for (k = 0;k < nb_state;k++) {
        *states += *(*ptransition + state_index[order - 1]) * *pstates++;
        ptransition++;
      }
      states++;

      // mise a jour des indices des etats

      for (k = 0;k < order;k++) {
        if (state_index[k] < nb_state - 1) {
          state_index[k]++;
          break;
        }
        else {
          state_index[k] = 0;
        }
      }
    }

    // mise a jour des probabilites des sequences d'etats et
    // calcul des probabilites des etats

    states = state_seq;
    pstates = pstate_seq;

    for (j = 0;j < nb_state;j++) {
      sum = 0.;
      for (k = 0;k < power[order - 1];k++) {
        sum += *states;
        *pstates++ = *states++;
      }
      index_state->point[j][i] = sum;
    }
  }

  delete [] state_seq;
  delete [] pstate_seq;
}


/*--------------------------------------------------------------*
 *
 *  Calcul des probabilites de chaque memoire pour une chaine de Markov
 *  en tenant compte de la distribution des longueurs des sequences.
 *
 *--------------------------------------------------------------*/

double* Markov::memory_computation() const

{
  register int i , j , k;
  int power[ORDER] , state_index[ORDER];
  double *memory , *state_seq , *pstate_seq , *states , *pstates , **ptransition;


  memory = new double[nb_row];
  for (i = 0;i < nb_row;i++) {
    memory[i] = 0.;
  }

  i = 1;
  for (j = 0;j < order;j++) {
    power[j] = i;
    i *= nb_state;
  }

  // initialisation des probabilites des sequences d'etats

  state_seq = new double[nb_row];
  pstate_seq = new double[nb_row];

  pstates = pstate_seq;
  i = 0;

  for (j = 0;j < nb_row;j++) {
    if (j == self_row[i]) {
      *pstates++ = initial[i++];
    }
    else {
      *pstates++ = 0.;
    }
  }

  // calcul des probabilites de chaque memoire en fonction de l'index

  for (i = 1;i < process[0]->length->nb_value - 2;i++) {

    // calcul des probabilites des sequences d'etats de
    // longueur = ordre de la chaine de Markov

    for (j = 0;j < order;j++) {
      state_index[j] = 0;
    }

    states = state_seq;

    for (j = 0;j < nb_row;j++) {
      ptransition = transition;
      pstates = pstate_seq;
      for (k = 0;k < order - 1;k++) {
        ptransition += state_index[k] * power[k + 1];
        pstates += state_index[k] * power[k + 1];
      }

      *states = 0.;
      for (k = 0;k < nb_state;k++) {
        *states += *(*ptransition + state_index[order - 1]) * *pstates++;
        ptransition++;
      }
      states++;

      // mise a jour des indices des etats

      for (k = 0;k < order;k++) {
        if (state_index[k] < nb_state - 1) {
          state_index[k]++;
          break;
        }
        else {
          state_index[k] = 0;
        }
      }
    }

    // mise a jour des probabilites des sequences d'etats et
    // accumulation des probabilites des memoires

    states = state_seq;
    pstates = pstate_seq;

    for (j = 0;j < nb_row;j++) {
      memory[j] += *states * (1. - process[0]->length->cumul[i]);
      *pstates++ = *states++;
    }
  }

  delete [] state_seq;
  delete [] pstate_seq;

  return memory;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la probabilite de ne pas observer un etat d'une chaine de Markov.
 *
 *  arguments : etat, seuil sur la somme des probabilites des sequences d'etats.
 *
 *--------------------------------------------------------------*/

void Markov::state_no_occurrence_probability(int state , double increment)

{
  register int i;

  for (i = 0;i < nb_state;i++) {
    if ((i != state) && (!accessibility[i][state])) {
      break;
    }
  }

  if (i < nb_state) {
    register int j , k;
    int power[ORDER] , state_index[ORDER];
    double state_sum , *state_seq , *pstate_seq , *states , *pstates , **ptransition ,
           &no_occurrence = process[0]->no_occurrence[state];


    i = 1;
    for (j = 0;j < order;j++) {
      power[j] = i;
      i *= nb_state;
    }

    // initialisation des probabilites des sequences d'etats

    state_seq = new double[nb_row];
    pstate_seq = new double[nb_row];

    states = state_seq;
    state_sum = 0.;
    no_occurrence = 0.;
    i = 0;

    for (j = 0;j < nb_row;j++) {
      if (j == self_row[i]) {
        if (i != state) {
          if (accessibility[i][state]) {
            *states = initial[i];
            state_sum += *states++;
          }
          else {
            *states++ = 0.;
            no_occurrence += initial[i];
          }
        }

        else {
          *states++ = 0.;
        }
        i++;
      }

      else {
        *states++ = 0.;
      }
    }

    i = 1;

    while ((state_sum > increment) || (i < (nb_state - 1) * order)) {

      // mise a jour des probabilites des sequences d'etats

      states = state_seq;
      pstates = pstate_seq;
      for (j = 0;j < nb_row;j++) {
        *pstates++ = *states++;
      }

      // calcul des probabilites des sequences d'etats et
      // mise a jour de la probabilite de ne pas observer l'etat

      for (j = 0;j < order;j++) {
        state_index[j] = 0;
      }

      states = state_seq;
      state_sum = 0.;

      for (j = 0;j < nb_row;j++) {
        *states = 0.;

        for (k = 1;k < MIN(i , order);k++) {
          if ((state_index[order - 1 - k] == state) ||
              (!accessibility[state_index[order - 1 - k]][state])) {
            break;
          }
        }

        if ((k == MIN(i , order)) && (state_index[order - 1] != state)) {
          ptransition = transition;
          pstates = pstate_seq;
          for (k = 0;k < order - 1;k++) {
            ptransition += state_index[k] * power[k + 1];
            pstates += state_index[k] * power[k + 1];
          }

          if (accessibility[state_index[order - 1]][state]) {
            for (k = 0;k < nb_state;k++) {
              if (k != state) {
                *states += *(*ptransition + state_index[order - 1]) * *pstates;
              }
              ptransition++;
              pstates++;
            }
          }

          else {
            for (k = 0;k < nb_state;k++) {
              if (k != state) {
                no_occurrence += *(*ptransition + state_index[order - 1]) * *pstates;
              }
              ptransition++;
              pstates++;
            }
          }
        }

        state_sum += *states++;

        // mise a jour des indices des etats

        for (k = 0;k < order;k++) {
          if (state_index[k] < nb_state - 1) {
            state_index[k]++;
            break;
          }
          else {
            state_index[k] = 0;
          }
        }
      }

      i++;
    }

    delete [] state_seq;
    delete [] pstate_seq;
  }
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la loi du temps avant la premiere occurrence d'un etat
 *  pour une chaine de Markov.
 *
 *  arguments : etat, nombre minimum de valeurs,
 *              seuil sur la fonction de repartition.
 *
 *--------------------------------------------------------------*/

void Markov::state_first_occurrence_distribution(int state , int min_nb_value ,
                                                 double cumul_threshold)

{
  register int i , j , k;
  int power[ORDER] , state_index[ORDER];
  double *state_seq , *pstate_seq , *states , *pstates , *pmass , *pcumul ,
         **ptransition;
  Distribution *first_occurrence;


  first_occurrence = process[0]->first_occurrence[state];
  first_occurrence->complement = process[0]->no_occurrence[state];

  pmass = first_occurrence->mass;
  pcumul = first_occurrence->cumul;

  i = 1;
  for (j = 0;j < order;j++) {
    power[j] = i;
    i *= nb_state;
  }

  // initialisation des probabilites des sequences d'etats

  state_seq = new double[nb_row];
  pstate_seq = new double[nb_row];

  states = state_seq;
  i = 0;

  for (j = 0;j < nb_row;j++) {
    if (j == self_row[i]) {
      if (i != state) {
        *states++ = initial[i++];
      }
      else {
        *states++ = 0.;
        *pmass = initial[i++];
      }
    }
    else {
      *states++ = 0.;
    }
  }
  *pcumul = *pmass;

  i = 1;

  while (((*pcumul < cumul_threshold - first_occurrence->complement) || (i < min_nb_value)) &&
         (i < first_occurrence->alloc_nb_value)) {

    // mise a jour des probabilites des sequences d'etats

    states = state_seq;
    pstates = pstate_seq;
    for (j = 0;j < nb_row;j++) {
      *pstates++ = *states++;
    }

    // calcul des probabilites des sequences d'etats et de la valeur courante

    for (j = 0;j < order;j++) {
      state_index[j] = 0;
    }

    states = state_seq;
    *++pmass = 0.;

    for (j = 0;j < nb_row;j++) {
      *states = 0.;

      for (k = 1;k < MIN(i , order);k++) {
        if (state_index[order - 1 - k] == state) {
          break;
        }
      }

      if (k == MIN(i , order)) {
        ptransition = transition;
        pstates = pstate_seq;
        for (k = 0;k < order - 1;k++) {
          ptransition += state_index[k] * power[k + 1];
          pstates += state_index[k] * power[k + 1];
        }

        if (state_index[order - 1] != state) {
          for (k = 0;k < nb_state;k++) {
            if (k != state) {
              *states += *(*ptransition + state_index[order - 1]) * *pstates;
            }
            ptransition++;
            pstates++;
          }
        }

        else {
          for (k = 0;k < nb_state;k++) {
            if (k != state) {
              *pmass += *(*ptransition + state_index[order - 1]) * *pstates;
            }
            ptransition++;
            pstates++;
          }
        }
      }

      states++;

      // mise a jour des indices des etats

      for (k = 0;k < order;k++) {
        if (state_index[k] < nb_state - 1) {
          state_index[k]++;
          break;
        }
        else {
          state_index[k] = 0;
        }
      }
    }

    // mise a jour de la fonction de repartition

    pcumul++;
    *pcumul = *(pcumul - 1) + *pmass;
    i++;
  }

  first_occurrence->nb_value = i;

# ifdef DEBUG
  if (first_occurrence->complement > 0.) {
    cout << "\n" << SEQ_label[SEQL_STATE_NO_OCCURRENCE] << " " << state << " : "
         << first_occurrence->complement << " | "
         << 1. - first_occurrence->cumul[first_occurrence->nb_value - 1] << endl;
  }
# endif

  first_occurrence->offset_computation();
  first_occurrence->max_computation();
  first_occurrence->mean_computation();
  first_occurrence->variance_computation();

  delete [] state_seq;
  delete [] pstate_seq;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la probabilite de quitter un etat sans possibilite
 *  d'y revenir pour une chaine de Markov.
 *
 *  arguments : loi des memoires, etat,
 *              seuil sur la somme des probabilites des sequences d'etats.
 *
 *--------------------------------------------------------------*/

void Markov::state_leave_probability(const double *memory , int state ,
                                     double increment)

{
  if (state_type[state] == 't') {
    register int i , j , k;
    int power[ORDER] , state_index[ORDER];
    double state_sum , *state_seq , *pstate_seq , *states , *pstates ,
           **ptransition , &leave = process[0]->leave[state];


    i = 1;
    for (j = 0;j < order;j++) {
      power[j] = i;
      i *= nb_state;
    }

    state_seq = new double[nb_row];
    pstate_seq = new double[nb_row];

    // initialisation des probabilites des sequences d'etats

    for (i = 0;i < order;i++) {
      state_index[i] = 0;
    }

    states = state_seq;
    state_sum = 0.;

    for (i = 0;i < nb_row;i++) {
      if (state_index[order - 1] == state) {
        if (i == self_row[state]) {
          *states = initial[state];
        }
        else {
          *states = 0.;
        }
        *states += memory[i];
        state_sum += *states++;
      }

      else {
        *states++ = 0.;
      }

      // mise a jour des indices des etats

      for (j = 0;j < order;j++) {
        if (state_index[j] < nb_state - 1) {
          state_index[j]++;
          break;
        }
        else {
          state_index[j] = 0;
        }
      }
    }

    states = state_seq;
    for (i = 0;i < nb_row;i++) {
      *states++ /= state_sum;
    }

    leave = 0.;
    i = 1;

    do {

      // mise a jour des probabilites des sequences d'etats

      pstates = pstate_seq;
      states = state_seq;
      for (j = 0;j < nb_row;j++) {
        *pstates++ = *states++;
      }

      // calcul des probabilites des sequences d'etats et mise a jour
      // de la probabilite de quitter l'etat

      for (j = 0;j < order;j++) {
        state_index[j] = 0;
      }

      states = state_seq;
      state_sum = 0.;

      for (j = 0;j < nb_row;j++) {
        *states = 0.;

        for (k = 1;k < MIN(i , order);k++) {
          if ((state_index[order - 1 - k] == state) || (!accessibility[state_index[order - 1 - k]][state])) {
            break;
          }
        }

        if ((((k == i) && (i < order) && (state_index[order - 1 - k] == state)) ||
            (k == order)) && (state_index[order - 1] != state)) {
          ptransition = transition;
          pstates = pstate_seq;
          for (k = 0;k < order - 1;k++) {
            ptransition += state_index[k] * power[k + 1];
            pstates += state_index[k] * power[k + 1];
          }

          if (accessibility[state_index[order - 1]][state]) {
            for (k = 0;k < nb_state;k++) {
              *states += *(*ptransition + state_index[order - 1]) * *pstates++;
              ptransition++;
            }
            state_sum += *states;
          }

          else {
            for (k = 0;k < nb_state;k++) {
              leave += *(*ptransition + state_index[order - 1]) * *pstates++;
              ptransition++;
            }
          }
        }

        states++;

        // mise a jour des indices des etats

        for (k = 0;k < order;k++) {
          if (state_index[k] < nb_state - 1) {
            state_index[k]++;
            break;
          }
          else {
            state_index[k] = 0;
          }
        }
      }

      i++;
    }
    while ((state_sum > increment) || (i < (nb_state - 1) * order));

    delete [] state_seq;
    delete [] pstate_seq;
  }
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la loi du temps de retour dans un etat
 *  pour une chaine de Markov.
 *
 *  arguments : loi des memoires, etat, nombre minimum de valeurs,
 *              seuil sur la fonction de repartition.
 *
 *--------------------------------------------------------------*/

void Markov::state_recurrence_time_distribution(const double *memory , int state ,
                                                int min_nb_value , double cumul_threshold)

{
  register int i , j , k;
  int power[ORDER] , state_index[ORDER];
  double sum , *state_seq , *pstate_seq , *states , *pstates ,
         *pmass , *pcumul , **ptransition;
  Distribution *recurrence_time;


  recurrence_time = process[0]->recurrence_time[state];
  recurrence_time->complement = process[0]->leave[state];

  pmass = recurrence_time->mass;
  pcumul = recurrence_time->cumul;
  *pmass = 0.;
  *pcumul = 0.;

  i = 1;
  for (j = 0;j < order;j++) {
    power[j] = i;
    i *= nb_state;
  }

  state_seq = new double[nb_row];
  pstate_seq = new double[nb_row];

  // initialisation des probabilites des sequences d'etats

  for (i = 0;i < order;i++) {
    state_index[i] = 0;
  }

  states = state_seq;
  sum = 0.;

  for (i = 0;i < nb_row;i++) {
    if (state_index[order - 1] == state) {
      if (i == self_row[state]) {
        *states = initial[state];
      }
      else {
        *states = 0.;
      }
      *states += memory[i];
      sum += *states++;
    }

    else {
      *states++ = 0.;
    }

    // mise a jour des indices des etats

    for (j = 0;j < order;j++) {
      if (state_index[j] < nb_state - 1) {
        state_index[j]++;
        break;
      }
      else {
        state_index[j] = 0;
      }
    }
  }

  states = state_seq;
  for (i = 0;i < nb_row;i++) {
    *states++ /= sum;
  }

  i = 1;

  do {

    // mise a jour des probabilites des sequences d'etats

    pstates = pstate_seq;
    states = state_seq;
    for (j = 0;j < nb_row;j++) {
      *pstates++ = *states++;
    }

    // calcul des probabilites des sequences d'etats et de la valeur courante

    for (j = 0;j < order;j++) {
      state_index[j] = 0;
    }

    states = state_seq;
    *++pmass = 0.;

    for (j = 0;j < nb_row;j++) {
      *states = 0.;

      for (k = 1;k < MIN(i , order);k++) {
        if (state_index[order - 1 - k] == state) {
          break;
        }
      }

      if (((k == i) && (i < order) && (state_index[order - 1 - k] == state)) ||
          (k == order)) {
        ptransition = transition;
        pstates = pstate_seq;
        for (k = 0;k < order - 1;k++) {
          ptransition += state_index[k] * power[k + 1];
          pstates += state_index[k] * power[k + 1];
        }

        if (state_index[order - 1] != state) {
          for (k = 0;k < nb_state;k++) {
            *states += *(*ptransition + state_index[order - 1]) * *pstates++;
            ptransition++;
          }
        }

        else {
          for (k = 0;k < nb_state;k++) {
            *pmass += *(*ptransition + state_index[order - 1]) * *pstates++;
            ptransition++;
          }
        }
      }

      states++;

      // mise a jour des indices des etats

      for (k = 0;k < order;k++) {
        if (state_index[k] < nb_state - 1) {
          state_index[k]++;
          break;
        }
        else {
          state_index[k] = 0;
        }
      }
    }

    // mise a jour de la fonction de repartition

    pcumul++;
    *pcumul = *(pcumul - 1) + *pmass;
    i++;
  }
  while (((*pcumul < cumul_threshold - recurrence_time->complement) || (i < min_nb_value)) &&
         (i < recurrence_time->alloc_nb_value));

  recurrence_time->nb_value = i;
  recurrence_time->nb_value_computation();

  if (recurrence_time->nb_value > 0) {
    recurrence_time->offset_computation();
    recurrence_time->max_computation();
    recurrence_time->mean_computation();
    recurrence_time->variance_computation();
  }

  else {
    delete process[0]->recurrence_time[state];
    process[0]->recurrence_time[state] = 0;
  }

  delete [] state_seq;
  delete [] pstate_seq;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la loi d'occupation d'un etat d'une chaine de Markov.
 *
 *  arguments : loi des memoires, etat, nombre minimum de valeurs,
 *              seuil sur la fonction de repartition.
 *
 *--------------------------------------------------------------*/

void Markov::state_sojourn_time_distribution(const double *memory , int state ,
                                             int min_nb_value , double cumul_threshold)

{
  register int i , j , k;
  int power[ORDER] , state_index[ORDER];
  double sum , *state_seq , *pstate_seq , *states , *pstates ,
         *pmass , *pcumul , **ptransition;
  Parametric *sojourn_time;


  sojourn_time = process[0]->sojourn_time[state];

  pmass = sojourn_time->mass;
  pcumul = sojourn_time->cumul;
  *pmass = 0.;
  *pcumul = 0.;

  i = 1;
  for (j = 0;j < order;j++) {
    power[j] = i;
    i *= nb_state;
  }

  state_seq = new double[nb_row];
  pstate_seq = new double[nb_row];

  // initialisation des probabilites des sequences d'etats

  for (i = 0;i < order;i++) {
    state_index[i] = 0;
  }

  states = state_seq;
  sum = 0.;

  for (i = 0;i < nb_row;i++) {
    if (state_index[order - 1] == state) {
      if (i == self_row[state]) {
        *states = initial[state];
      }
      else {
        *states = 0.;
      }
      if ((order == 1) || ((order > 1) && (state_index[order - 2] != state))) {
        *states += memory[i];
      }
      sum += *states++;
    }

    else {
      *states++ = 0.;
    }

    // mise a jour des indices des etats

    for (j = 0;j < order;j++) {
      if (state_index[j] < nb_state - 1) {
        state_index[j]++;
        break;
      }
      else {
        state_index[j] = 0;
      }
    }
  }

  states = state_seq;
  for (i = 0;i < nb_row;i++) {
    *states++ /= sum;
  }

  // temps < ordre de la chaine de Markov

  for (i = 1;i < order;i++) {

    // mise a jour des probabilites des sequences d'etats

    pstates = pstate_seq;
    states = state_seq;
    for (j = 0;j < nb_row;j++) {
      *pstates++ = *states++;
    }

    // calcul des probabilites des sequences d'etats et de la valeur courante

    for (j = 0;j < order;j++) {
      state_index[j] = 0;
    }

    states = state_seq;
    *++pmass = 0.;

    for (j = 0;j < nb_row;j++) {
      *states = 0.;

      for (k = 0;k <= i;k++) {
        if (state_index[order - 1 - k] != state) {
          break;
        }
      }

      if ((k == i + 1) && ((j == self_row[state]) || (i == order - 1) ||
           ((i < order - 1) && (state_index[order - 1 - k] != state)))) {
        ptransition = transition;
        pstates = pstate_seq;
        for (k = 0;k < order - 1;k++) {
          ptransition += state_index[k] * power[k + 1];
          pstates += state_index[k] * power[k + 1];
        }

        for (k = 0;k < nb_state;k++) {
          *states += *(*ptransition + state_index[order - 1]) * *pstates;
          *pmass += (1. - *(*ptransition + state_index[order - 1])) * *pstates++;
          ptransition++;
        }
      }

      states++;

      // mise a jour des indices des etats

      for (k = 0;k < order;k++) {
        if (state_index[k] < nb_state - 1) {
          state_index[k]++;
          break;
        }
        else {
          state_index[k] = 0;
        }
      }
    }

    // mise a jour de la fonction de repartition

    pcumul++;
    *pcumul = *(pcumul - 1) + *pmass;
  }

  *state_seq = *(state_seq + self_row[state]);

  // calcul des probabilites des valeurs de la traine geometrique

  ptransition = transition + self_row[state];

  while (((*pcumul < cumul_threshold) || (i < min_nb_value)) &&
         (i < sojourn_time->alloc_nb_value)) {
    *++pmass = *state_seq * (1. - *(*ptransition + state));
    *state_seq *= *(*ptransition + state);
    pcumul++;
    *pcumul = *(pcumul - 1) + *pmass;
    i++;
  }
  sojourn_time->nb_value = i;

  sojourn_time->offset_computation();
  sojourn_time->max_computation();
  sojourn_time->mean_computation();
  sojourn_time->variance_computation();

  delete [] state_seq;
  delete [] pstate_seq;
}


/*--------------------------------------------------------------*
 *
 *  Calcul des probabilites de chaque observation en fonction de l'index
 *  pour une chaine de Markov cachee.
 *
 *  argument : indice du processus d'observation.
 *
 *--------------------------------------------------------------*/

void Markov::index_output_distribution(int variable)

{
  register int i , j , k;
  Curves *index_state , *index_value;


  index_value = process[variable]->index_value;

  // calcul des probabilites des etats de la chaine de Markov
  // sous-jacente en fonction de l'index si necessaire

  if (!(process[0]->index_value)) {
    process[0]->index_value = new Curves(nb_state , index_value->length);
    index_state_distribution();
  }
  index_state = process[0]->index_value;

  // prise en compte des probabilites d'observation

  for (i = 0;i < index_value->length;i++) {
    for (j = 0;j < process[variable]->nb_value;j++) {
      index_value->point[j][i] = 0.;
      for (k = 0;k < nb_state;k++) {
        index_value->point[j][i] += process[variable]->observation[k]->mass[j] *
                                    index_state->point[k][i];
      }
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la probabilite de ne pas observer une valeur
 *  pour une chaine de Markov cachee.
 *
 *  arguments : indice du processus d'observation, observation,
 *              seuil sur la somme des probabilites des sequences d'etats.
 *
 *--------------------------------------------------------------*/

void Markov::output_no_occurrence_probability(int variable , int output , double increment)

{
  bool status = false , *output_accessibility;
  register int i , j , k;
  int power[ORDER] , state_index[ORDER];
  double state_sum , sum , *observation , *state_seq , *pstate_seq , *states , *pstates ,
         **ptransition , &no_occurrence = process[variable]->no_occurrence[output];


  observation = new double[nb_state];
  for (i = 0;i < nb_state;i++) {
    observation[i] = process[variable]->observation[i]->mass[output];
  }

  // calcul de l'accessibilite d'une observation a partir d'un etat donne

  output_accessibility = new bool[nb_state];

  for (i = 0;i < nb_state;i++) {
    output_accessibility[i] = false;

    for (j = 0;j < nb_state;j++) {
      if (j == i) {
        if (observation[j] > 0.) {
          output_accessibility[i] = true;
          break;
        }
      }

      else {
        if ((accessibility[i][j]) && (observation[j] > 0.)) {
          output_accessibility[i] = true;
          break;
        }
      }
    }

    if (!output_accessibility[i]) {
      status = true;
    }
  }

  if (status) {
    i = 1;
    for (j = 0;j < order;j++) {
      power[j] = i;
      i *= nb_state;
    }

    // initialisation des probabilites des sequences d'etats

    state_seq = new double[nb_row];
    pstate_seq = new double[nb_row];

    states = state_seq;
    state_sum = 0.;
    no_occurrence = 0.;
    i = 0;

    for (j = 0;j < nb_row;j++) {
      if (j == self_row[i]) {
        if (output_accessibility[i]) {
          *states = (1. - observation[i]) * initial[i];
          state_sum += *states++;
        }
        else {
          *states++ = 0.;
          no_occurrence += initial[i];
        }
        i++;
      }

      else {
        *states++ = 0.;
      }
    }

    i = 1;

    while ((state_sum > increment) || (i < nb_state * order)) {

      // mise a jour des probabilites des sequences d'etats

      states = state_seq;
      pstates = pstate_seq;
      for (j = 0;j < nb_row;j++) {
        *pstates++ = *states++;
      }

      // calcul des probabilites des sequences d'etats et
      // mise a jour de la probabilite de ne pas observer la valeur

      for (j = 0;j < order;j++) {
        state_index[j] = 0;
      }

      states = state_seq;
      state_sum = 0.;

      for (j = 0;j < nb_row;j++) {
        ptransition = transition;
        pstates = pstate_seq;
        for (k = 0;k < order - 1;k++) {
          ptransition += state_index[k] * power[k + 1];
          pstates += state_index[k] * power[k + 1];
        }

        sum = 0.;
        for (k = 0;k < nb_state;k++) {
          sum += *(*ptransition + state_index[order - 1]) * *pstates++;
          ptransition++;
        }

        if (output_accessibility[state_index[order - 1]]) {
          *states = (1. - observation[state_index[order - 1]]) * sum;
          state_sum += *states++;
        }
        else {
          *states++ = 0.;
          no_occurrence += sum;
        }

        // mise a jour des indices des etats

        for (k = 0;k < order;k++) {
          if (state_index[k] < nb_state - 1) {
            state_index[k]++;
            break;
          }
          else {
            state_index[k] = 0;
          }
        }
      }

      i++;
    }

    delete [] state_seq;
    delete [] pstate_seq;
  }

  delete [] observation;
  delete [] output_accessibility;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la loi du temps avant la premiere occurrence d'une observation
 *  pour une chaine de Markov cachee.
 *
 *  arguments : indice du processus d'observation, observation, nombre minimum de valeurs,
 *              seuil sur la fonction de repartition.
 *
 *--------------------------------------------------------------*/

void Markov::output_first_occurrence_distribution(int variable , int output ,
                                                  int min_nb_value , double cumul_threshold)

{
  register int i , j , k;
  int power[ORDER] , state_index[ORDER];
  double sum , *observation , *state_seq , *pstate_seq , *states , *pstates ,
         *pmass , *pcumul , **ptransition;
  Distribution *first_occurrence;


  first_occurrence = process[variable]->first_occurrence[output];
  first_occurrence->complement = process[variable]->no_occurrence[output];

  pmass = first_occurrence->mass;
  pcumul = first_occurrence->cumul;

  observation = new double[nb_state];
  for (i = 0;i < nb_state;i++) {
    observation[i] = process[variable]->observation[i]->mass[output];
  }

  i = 1;
  for (j = 0;j < order;j++) {
    power[j] = i;
    i *= nb_state;
  }

  // initialisation des probabilites des sequences d'etats

  state_seq = new double[nb_row];
  pstate_seq = new double[nb_row];

  states = state_seq;
  *pmass = 0.;
  i = 0;

  for (j = 0;j < nb_row;j++) {
    if (j == self_row[i]) {
      *states++ = (1. - observation[i]) * initial[i];
      *pmass += observation[i] * initial[i];
      i++;
    }
    else {
      *states++ = 0.;
    }
  }
  *pcumul = *pmass;

  i = 1;

  while (((*pcumul < cumul_threshold - first_occurrence->complement) || (i < min_nb_value)) &&
         (i < first_occurrence->alloc_nb_value)) {

    // mise a jour des probabilites des sequences d'etats

    states = state_seq;
    pstates = pstate_seq;
    for (j = 0;j < nb_row;j++) {
      *pstates++ = *states++;
    }

    // calcul des probabilites des sequences d'etats et de la valeur courante

    for (j = 0;j < order;j++) {
      state_index[j] = 0;
    }

    states = state_seq;
    *++pmass = 0.;

    for (j = 0;j < nb_row;j++) {
      ptransition = transition;
      pstates = pstate_seq;
      for (k = 0;k < order - 1;k++) {
        ptransition += state_index[k] * power[k + 1];
        pstates += state_index[k] * power[k + 1];
      }

      sum = 0.;
      for (k = 0;k < nb_state;k++) {
        sum += *(*ptransition + state_index[order - 1]) * *pstates++;
        ptransition++;
      }

      *states++ = (1. - observation[state_index[order - 1]]) * sum;
      *pmass += observation[state_index[order - 1]] * sum;

      // mise a jour des indices des etats

      for (k = 0;k < order;k++) {
        if (state_index[k] < nb_state - 1) {
          state_index[k]++;
          break;
        }
        else {
          state_index[k] = 0;
        }
      }
    }

    // mise a jour de la fonction de repartition

    pcumul++;
    *pcumul = *(pcumul - 1) + *pmass;
    i++;
  }

  first_occurrence->nb_value = i;

  first_occurrence->offset_computation();
  first_occurrence->max_computation();
  first_occurrence->mean_computation();
  first_occurrence->variance_computation();

  delete [] observation;
  delete [] state_seq;
  delete [] pstate_seq;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la probabilite de quitter une observation sans possibilite
 *  d'y revenir pour une chaine de Markov cachee.
 *
 *  arguments : loi des memoires, indice du processus d'observation, observation,
 *              seuil sur la somme des probabilites des sequences d'etats.
 *
 *--------------------------------------------------------------*/

void Markov::output_leave_probability(const double *memory , int variable ,
                                      int output , double increment)

{
  bool status = false , *output_accessibility;
  register int i , j , k;
  int power[ORDER] , state_index[ORDER];
  double state_sum , sum , *observation , *state_seq , *pstate_seq , *states , *pstates ,
         **ptransition , &leave = process[variable]->leave[output];


  observation = new double[nb_state];
  for (i = 0;i < nb_state;i++) {
    observation[i] = process[variable]->observation[i]->mass[output];
  }

  // calcul de l'accessibilite d'une observation a partir d'un etat donne

  output_accessibility = new bool[nb_state];

  for (i = 0;i < nb_state;i++) {
    output_accessibility[i] = false;

    for (j = 0;j < nb_state;j++) {
      if (j == i) {
        if (observation[j] > 0.) {
          output_accessibility[i] = true;
          break;
        }
      }

      else {
        if ((accessibility[i][j]) && (observation[j] > 0.)) {
          output_accessibility[i] = true;
          break;
        }
      }
    }

    if (!output_accessibility[i]) {
      status = true;
    }
  }

  if (status) {
    i = 1;
    for (j = 0;j < order;j++) {
      power[j] = i;
      i *= nb_state;
    }

    state_seq = new double[nb_row];
    pstate_seq = new double[nb_row];

    // initialisation des probabilites des sequences d'etats

    for (i = 0;i < order;i++) {
      state_index[i] = 0;
    }

    states = state_seq;
    state_sum = 0.;
    i = 0;

    for (j = 0;j < nb_row;j++) {
      if (j == self_row[i]) {
        *states = initial[i++];
      }
      else {
        *states = 0.;
      }
      *states += memory[j];
      *states *= observation[state_index[order - 1]];
      state_sum += *states++;

      // mise a jour des indices des etats

      for (k = 0;k < order;k++) {
        if (state_index[k] < nb_state - 1) {
          state_index[k]++;
          break;
        }
        else {
          state_index[k] = 0;
        }
      }
    }

    states = state_seq;
    for (i = 0;i < nb_row;i++) {
      *states++ /= state_sum;
    }

    leave = 0.;
    i = 1;

    do {

      // mise a jour des probabilites des sequences d'etats

      pstates = pstate_seq;
      states = state_seq;
      for (j = 0;j < nb_row;j++) {
        *pstates++ = *states++;
      }

      // calcul des probabilites des sequences d'etats et mise a jour
      // de la probabilite de quitter l'observation

      for (j = 0;j < order;j++) {
        state_index[j] = 0;
      }

      states = state_seq;
      state_sum = 0.;

      for (j = 0;j < nb_row;j++) {
        *states = 0.;

        if (observation[state_index[order - 1]] < 1.) {
          ptransition = transition;
          pstates = pstate_seq;
          for (k = 0;k < order - 1;k++) {
            ptransition += state_index[k] * power[k + 1];
            pstates += state_index[k] * power[k + 1];
          }

          sum = 0.;
          for (k = 0;k < nb_state;k++) {
            sum += *(*ptransition + state_index[order - 1]) * *pstates++;
            ptransition++;
          }

          if (output_accessibility[state_index[order - 1]]) {
            *states = (1. - observation[state_index[order - 1]]) * sum;
            state_sum += *states;
          }
          else {
            leave += sum;
          }
        }

        states++;

        // mise a jour des indices des etats

        for (k = 0;k < order;k++) {
          if (state_index[k] < nb_state - 1) {
            state_index[k]++;
            break;
          }
          else {
            state_index[k] = 0;
          }
        }
      }

      i++;
    }
    while ((state_sum > increment) || (i < nb_state * order));

    delete [] state_seq;
    delete [] pstate_seq;
  }

  delete [] observation;
  delete [] output_accessibility;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la loi du temps de retour dans une observation
 *  pour une chaine de Markov cachee.
 *
 *  arguments : loi des memoires, indice du processus d'observation,
 *              observation, nombre minimum de valeurs,
 *              seuil sur la fonction de repartition.
 *
 *--------------------------------------------------------------*/

void Markov::output_recurrence_time_distribution(const double *memory , int variable ,
                                                 int output , int min_nb_value , double cumul_threshold)

{
  register int i , j , k;
  int power[ORDER] , state_index[ORDER];
  double sum , *observation , *state_seq , *pstate_seq , *states , *pstates ,
         *pmass , *pcumul , **ptransition;
  Distribution *recurrence_time;


  recurrence_time = process[variable]->recurrence_time[output];
  recurrence_time->complement = process[variable]->leave[output];

  pmass = recurrence_time->mass;
  pcumul = recurrence_time->cumul;
  *pmass = 0.;
  *pcumul = 0.;

  observation = new double[nb_state];
  for (i = 0;i < nb_state;i++) {
    observation[i] = process[variable]->observation[i]->mass[output];
  }

  i = 1;
  for (j = 0;j < order;j++) {
    power[j] = i;
    i *= nb_state;
  }

  state_seq = new double[nb_row];
  pstate_seq = new double[nb_row];

  // initialisation des probabilites des sequences d'etats

  for (i = 0;i < order;i++) {
    state_index[i] = 0;
  }

  states = state_seq;
  sum = 0.;
  i = 0;

  for (j = 0;j < nb_row;j++) {
    if (j == self_row[i]) {
      *states = initial[i++];
    }
    else {
      *states = 0.;
    }
    *states += memory[j];
    *states *= observation[state_index[order - 1]];
    sum += *states++;

    // mise a jour des indices des etats

    for (k = 0;k < order;k++) {
      if (state_index[k] < nb_state - 1) {
        state_index[k]++;
        break;
      }
      else {
        state_index[k] = 0;
      }
    }
  }

  states = state_seq;
  for (i = 0;i < nb_row;i++) {
    *states++ /= sum;
  }

  i = 1;

  do {

    // mise a jour des probabilites des sequences d'etats

    pstates = pstate_seq;
    states = state_seq;
    for (j = 0;j < nb_row;j++) {
      *pstates++ = *states++;
    }

    // calcul des probabilites des sequences d'etats et de la valeur courante

    for (j = 0;j < order;j++) {
      state_index[j] = 0;
    }

    states = state_seq;
    *++pmass = 0.;

    for (j = 0;j < nb_row;j++) {
      ptransition = transition;
      pstates = pstate_seq;
      for (k = 0;k < order - 1;k++) {
        ptransition += state_index[k] * power[k + 1];
        pstates += state_index[k] * power[k + 1];
      }

      sum = 0.;
      for (k = 0;k < nb_state;k++) {
        sum += *(*ptransition + state_index[order - 1]) * *pstates++;
        ptransition++;
      }

      *states++ = (1. - observation[state_index[order - 1]]) * sum;
      *pmass += observation[state_index[order - 1]] * sum;

      // mise a jour des indices des etats

      for (k = 0;k < order;k++) {
        if (state_index[k] < nb_state - 1) {
          state_index[k]++;
          break;
        }
        else {
          state_index[k] = 0;
        }
      }
    }

    // mise a jour de la fonction de repartition

    pcumul++;
    *pcumul = *(pcumul - 1) + *pmass;
    i++;
  }
  while (((*pcumul < cumul_threshold - recurrence_time->complement) || (i < min_nb_value)) &&
         (i < recurrence_time->alloc_nb_value));

  recurrence_time->nb_value = i;
  recurrence_time->nb_value_computation();

  if (recurrence_time->nb_value > 0) {
    recurrence_time->offset_computation();
    recurrence_time->max_computation();
    recurrence_time->mean_computation();
    recurrence_time->variance_computation();
  }

  else {
    delete process[variable]->recurrence_time[output];
    process[variable]->recurrence_time[output] = 0;
  }

  delete [] observation;
  delete [] state_seq;
  delete [] pstate_seq;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la loi du temps de sejour dans une observation
 *  pour une chaine de Markov cachee.
 *
 *  arguments : loi des memoires, indice du processus d'observation,
 *              observation, nombre minimum de valeurs,
 *              seuil sur la fonction de repartition.
 *
 *--------------------------------------------------------------*/

void Markov::output_sojourn_time_distribution(const double *memory , int variable ,
                                              int output , int min_nb_value , double cumul_threshold)

{
  register int i , j , k;
  int power[ORDER] , state_index[ORDER];
  double sum , *observation , *state_seq , *pstate_seq , *states , *pstates ,
         *pmass , *pcumul , **ptransition , &absorption = process[variable]->absorption[output];
  Parametric *sojourn_time;


  sojourn_time = process[variable]->sojourn_time[output];

  pmass = sojourn_time->mass;
  pcumul = sojourn_time->cumul;
  *pmass = 0.;
  *pcumul = 0.;

  observation = new double[nb_state];
  for (i = 0;i < nb_state;i++) {
    observation[i] = process[variable]->observation[i]->mass[output];
  }

  i = 1;
  for (j = 0;j < order;j++) {
    power[j] = i;
    i *= nb_state;
  }

  state_seq = new double[nb_row];
  pstate_seq = new double[nb_row];

  // initialisation des probabilites des sequences d'etats

  states = state_seq;
  sum = 0.;

  if (order == 1) {
    for (i = 0;i < nb_state;i++) {
      *states = initial[i];
      for (j = 0;j < nb_state;j++) {
        *states += transition[j][i] * (1. - observation[j]) * (initial[j] + memory[j]);
      }
      *states *= observation[i];
      sum += *states++;
    }
  }

  else {
    for (i = 0;i < order;i++) {
      state_index[i] = 0;
    }

    i = 0;
    for (j = 0;j < nb_row;j++) {
      if (j == self_row[i]) {
        *states = initial[i++];
      }
      else {
        *states = 0.;
      }
      *states += (1. - observation[state_index[order - 2]]) * memory[j];
      *states *= observation[state_index[order - 1]];
      sum += *states++;

      // mise a jour des indices des etats

      for (k = 0;k < order;k++) {
        if (state_index[k] < nb_state - 1) {
          state_index[k]++;
          break;
        }
        else {
          state_index[k] = 0;
        }
      }
    }
  }

  states = state_seq;
  for (i = 0;i < nb_row;i++) {
    *states++ /= sum;
  }

  i = 1;

  do {

    // mise a jour des probabilites des sequences d'etats

    pstates = pstate_seq;
    states = state_seq;
    for (j = 0;j < nb_row;j++) {
      *pstates++ = *states++;
    }

    // calcul des probabilites des sequences d'etats et de la valeur courante

    for (j = 0;j < order;j++) {
      state_index[j] = 0;
    }

    states = state_seq;
    absorption = 0.;
    *++pmass = 0.;

    for (j = 0;j < nb_row;j++) {
      ptransition = transition;
      pstates = pstate_seq;
      for (k = 0;k < order - 1;k++) {
        ptransition += state_index[k] * power[k + 1];
        pstates += state_index[k] * power[k + 1];
      }

      sum = 0.;
      for (k = 0;k < nb_state;k++) {
        sum += *(*ptransition + state_index[order - 1]) * *pstates++;
        ptransition++;
      }

      if ((state_type[state_index[order - 1]] == 'a') &&
          (observation[state_index[order - 1]] == 1.)) {
        absorption += sum;
      }

      *states++ = observation[state_index[order - 1]] * sum;
      *pmass += (1. - observation[state_index[order - 1]]) * sum;

      // mise a jour des indices des etats

      for (k = 0;k < order;k++) {
        if (state_index[k] < nb_state - 1) {
          state_index[k]++;
          break;
        }
        else {
          state_index[k] = 0;
        }
      }
    }

    // mise a jour de la fonction de repartition

    pcumul++;
    *pcumul = *(pcumul - 1) + *pmass;
    i++;
  }
  while (((*pcumul < cumul_threshold - absorption) || (i < min_nb_value)) &&
         (i < sojourn_time->alloc_nb_value));

  if (*pcumul == 0.) {
    absorption = 1.;
    delete process[variable]->sojourn_time[output];
    process[variable]->sojourn_time[output] = 0;
  }

  else {
    sojourn_time->nb_value = i;
    sojourn_time->complement = absorption;

#   ifdef DEBUG
    if (absorption > 0.) {
      cout << "\n" << SEQ_label[SEQL_OUTPUT_ABSORPTION] << " " << output << " : "
           << absorption << " | " << 1. - sojourn_time->cumul[sojourn_time->nb_value - 1] << endl;
    }
#   endif

    sojourn_time->offset_computation();
    sojourn_time->max_computation();
    sojourn_time->mean_computation();
    sojourn_time->variance_computation();
  }

  delete [] observation;
  delete [] state_seq;
  delete [] pstate_seq;
}
