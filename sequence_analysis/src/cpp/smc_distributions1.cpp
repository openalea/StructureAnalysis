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
 *       $Id: smc_distributions1.cpp 18077 2015-04-23 10:57:21Z guedon $
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
#include "stat_tool/stat_label.h"

#include "sequences.h"
#include "semi_markov.h"

using namespace std;
using namespace stat_tool;


namespace sequence_analysis {



/*--------------------------------------------------------------*
 *
 *  Calcul des probabilites de chaque etat en fonction du temps
 *  pour une semi-chaine de Markov.
 *
 *--------------------------------------------------------------*/

void SemiMarkovChain::index_state_distribution()

{
  register int i , j , k;
  double sum , *state_out , **state_in;
  Curves *index_state;
  DiscreteParametric *occupancy;


  index_state = state_process->index_value;

  state_out = new double[nb_state];

  state_in = new double*[index_state->length - 1];
  for (i = 0;i < index_state->length - 1;i++) {
    state_in[i] = new double[nb_state];
  }

  for (i = 0;i < index_state->length;i++) {
    for (j = 0;j < nb_state;j++) {
      switch (state_subtype[j]) {

      // cas etat semi-markovien

      case SEMI_MARKOVIAN : {
        if (i == 0) {
          index_state->point[j][i] = initial[j];
        }
        else {
          index_state->point[j][i] = state_in[i - 1][j] - state_out[j] + index_state->point[j][i - 1];
        }

        if (i < index_state->length - 1) {
          occupancy = state_process->sojourn_time[j];
          state_out[j] = 0.;
//          istate = 0.;

          for (k = 1;k <= MIN(i + 1 , occupancy->nb_value - 1);k++) {
            if (k < i + 1) {
              state_out[j] += occupancy->mass[k] * state_in[i - k][j];
//              istate += (1. - occupancy->cumul[k - 1]) * state_in[i - k][j];
            }
            else {
              switch (type) {
              case 'o' :
                state_out[j] += occupancy->mass[k] * initial[j];
//                istate += (1. - occupancy->cumul[k - 1]) * initial[j];
                break;
              case 'e' :
                state_out[j] += forward[j]->mass[k] * initial[j];
//                istate += (1. - forward[j]->cumul[k - 1]) * initial[j];
                break;
              }
            }
          }

//          index_state->point[j][i] = istate;
        }
        break;
      }

      // cas etat markovien

      case MARKOVIAN : {
        if (i == 0) {
          index_state->point[j][i] = initial[j];
        }
        else {
          index_state->point[j][i] = state_in[i - 1][j];
        }

        if (i < index_state->length - 1) {
          state_out[j] = index_state->point[j][i];
        }
        break;
      }
      }
    }

    if (i < index_state->length - 1) {
      for (j = 0;j < nb_state;j++) {
        state_in[i][j] = 0.;
        for (k = 0;k < nb_state;k++) {
          state_in[i][j] += transition[k][j] * state_out[k];
        }
      }
    }

    // renormalisation pour tenir compte des seuils appliques sur
    // les fonctions de repartition des lois d'occupation des etats

    sum = 0.;
    for (j = 0;j < nb_state;j++) {
      sum += index_state->point[j][i];
    }

    if (sum < 1.) {
      for (j = 0;j < nb_state;j++) {
        index_state->point[j][i] /= sum;
      }
    }
  }

  delete [] state_out;

  for (i = 0;i < index_state->length - 1;i++) {
    delete [] state_in[i];
  }
  delete [] state_in;
}


/*--------------------------------------------------------------*
 *
 *  Calcul des probabilites de chaque memoire pour une semi-chaine de Markov.
 *
 *--------------------------------------------------------------*/

double* SemiMarkovChain::memory_computation() const

{
  register int i , j , k;
  double sum , *memory , *state_out , **state_in;
  DiscreteParametric *occupancy;


  memory = new double[nb_state];
  state_out = new double[nb_state];

  switch (type) {

  case 'o' : {
    state_in = new double*[state_process->length->nb_value - 3];
    for (i = 0;i < state_process->length->nb_value - 3;i++) {
      state_in[i] = new double[nb_state];
    }

    for (i = 0;i < nb_state;i++) {
      memory[i] = 0.;
    }

    for (i = 0;i < state_process->length->nb_value - 2;i++) {
      for (j = 0;j < nb_state;j++) {
        switch (state_subtype[j]) {

        // cas etat semi-markovien

        case SEMI_MARKOVIAN : {
          occupancy = state_process->sojourn_time[j];
          state_out[j] = 0.;

          for (k = 1;k <= MIN(i + 1 , occupancy->nb_value - 1);k++) {
            if (k < i + 1) {
              state_out[j] += occupancy->mass[k] * state_in[i - k][j];
            }
            else {
              state_out[j] += occupancy->mass[k] * initial[j];
            }
          }
          break;
        }

        // cas etat markovien

        case MARKOVIAN : {
          if (i == 0) {
            state_out[j] = initial[j];
          }
          else {
            state_out[j] = state_in[i - 1][j];
          }
          break;
        }
        }

        // accumulation des probabilites des memoires

        memory[j] += state_out[j] * (1. - state_process->length->cumul[i + 1]);
      }

      if (i < state_process->length->nb_value - 3) {
        for (j = 0;j < nb_state;j++) {
          state_in[i][j] = 0.;
          for (k = 0;k < nb_state;k++) {
            state_in[i][j] += transition[k][j] * state_out[k];
          }
        }
      }
    }

    for (i = 0;i < state_process->length->nb_value - 3;i++) {
      delete [] state_in[i];
    }
    delete [] state_in;

    break;
  }

  case 'e' : {
    state_in = new double*[STATIONARY_PROBABILITY_LENGTH];
    for (i = 0;i < STATIONARY_PROBABILITY_LENGTH;i++) {
      state_in[i] = new double[nb_state];
    }

    i = 0;

    do {
      if (i > 0) {
        sum = 0.;
      }

      for (j = 0;j < nb_state;j++) {
        if (i > 0) {
          sum += fabs(state_in[i - 1][j] - state_out[j]);
        }

        switch (state_subtype[j]) {

        // cas etat semi-markovien

        case SEMI_MARKOVIAN : {
          occupancy = state_process->sojourn_time[j];
          state_out[j] = 0.;

          for (k = 1;k <= MIN(i + 1 , occupancy->nb_value - 1);k++) {
            if (k < i + 1) {
              state_out[j] += occupancy->mass[k] * state_in[i - k][j];
            }
            else {
              state_out[j] += forward[j]->mass[k] * initial[j];
            }
          }
          break;
        }

        // cas etat markovien

        case MARKOVIAN : {
          if (i == 0) {
            state_out[j] = initial[j];
          }
          else {
            state_out[j] = state_in[i - 1][j];
          }
          break;
        }
        }
      }

      for (j = 0;j < nb_state;j++) {
        state_in[i][j] = 0.;
        for (k = 0;k < nb_state;k++) {
          state_in[i][j] += transition[k][j] * state_out[k];
        }
      }

#     ifdef DEBUG
//      if ((i > 0) && (i % 100 == 0)) {
        cout << i << "  ";
        for (j = 0;j < nb_state;j++) {
          cout << state[j] << " ";
        }
        cout << " | " << sum / nb_state << endl;
//      }
#     endif

      i++;
    }
    while (((i == 1) || (sum / nb_state > STATIONARY_PROBABILITY_THRESHOLD)) &&
           (i < STATIONARY_PROBABILITY_LENGTH));

#   ifdef DEBUG
    cout << "\n" << SEQ_label[SEQL_LENGTH] << ": "  << i << endl;
#   endif

    for (j = 0;j < nb_state;j++) {
      memory[j] = state_in[i - 1][j];
    }

    for (i = 0;i < STATIONARY_PROBABILITY_LENGTH;i++) {
      delete [] state_in[i];
    }
    delete [] state_in;

    break;
  }
  }

  delete [] state_out;

  return memory;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la probabilite de ne pas observer un etat
 *  d'une semi-chaine de Markov ordinaire.
 *
 *  arguments : etat, seuil sur la somme des probabilites de quitter un etat.
 *
 *--------------------------------------------------------------*/

void SemiMarkovChain::state_no_occurrence_probability(int state , double increment)

{
  register int i;

  for (i = 0;i < nb_state;i++) {
    if ((i != state) && (!accessibility[i][state])) {
      break;
    }
  }

  if (i < nb_state) {
    register int j , k;
    int min_time;
    double sum , *state_out , **state_in ,
           &no_occurrence = state_process->no_occurrence[state];
    DiscreteParametric *occupancy;


    state_out = new double[nb_state];

    state_in = new double*[LEAVE_LENGTH];
    state_in[0] = NULL;
    for (i = 1;i < LEAVE_LENGTH;i++) {
      state_in[i] = new double[nb_state];
    }

    no_occurrence = 0.;
    for (i = 0;i < nb_state;i++) {
      if ((i != state) && (!accessibility[i][state])) {
        no_occurrence += initial[i];
      }
    }

    sum = 0.;
    for (i = 0;i < nb_state;i++) {
      if (i != state) {
        switch (state_subtype[i]) {

        case SEMI_MARKOVIAN : {
          sum += state_process->sojourn_time[i]->mean;
          break;
        }

        case MARKOVIAN : {
          if (transition[i][i] < 1.) {
            sum += 1. / (1. - transition[i][i]);
          }
          break;
        }
        }
      }
    }
    min_time = (int)sum + 1;

    i = 1;

    do {

      // calcul des probabilites de quitter (semi-Markov) / d'etre dans (Markov) un etat et
      // mise a jour de la probabilite de ne pas observer l'etat selectionne

      sum = 0.;

      for (j = 0;j < nb_state;j++) {
        if ((j != state) && (accessibility[j][state])) {
          switch (state_subtype[j]) {

          // cas etat semi-markovien

          case SEMI_MARKOVIAN : {
            occupancy = state_process->sojourn_time[j];
            state_out[j] = 0.;

            for (k = 1;k <= MIN(i , occupancy->nb_value - 1);k++) {
              if (k < i) {
                state_out[j] += occupancy->mass[k] * state_in[i - k][j];
              }
              else {
                state_out[j] += occupancy->mass[k] * initial[j];
              }
            }
            break;
          }

          // cas etat markovien

          case MARKOVIAN : {
            if (i == 1) {
              state_out[j] = initial[j];
            }
            else {
              state_out[j] = state_in[i - 1][j];
            }
            break;
          }
          }

          if ((transition[j][j] == 0.) || (transition[j][j] == 1.)) {
            sum += state_out[j];
          }
          else {
            sum += state_out[j] * (1. - transition[j][j]);
          }

          for (k = 0;k < nb_state;k++) {
            if ((k != state) && (!accessibility[k][state])) {
              no_occurrence += transition[j][k] * state_out[j];
            }
          }
        }
      }

      for (j = 0;j < nb_state;j++) {
        if ((j != state) && (accessibility[j][state])) {
          state_in[i][j] = 0.;
          for (k = 0;k < nb_state;k++) {
            if ((k != state) && (accessibility[k][state])) {
              state_in[i][j] += transition[k][j] * state_out[k];
            }
          }
        }
      }

      i++;
    }
    while (((sum > increment) || (i <= min_time)) && (i < LEAVE_LENGTH));

    delete [] state_out;

    for (i = 1;i < LEAVE_LENGTH;i++) {
      delete [] state_in[i];
    }
    delete [] state_in;
  }
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la loi du temps avant la 1ere occurrence d'un etat
 *  pour une semi-chaine de Markov.
 *
 *  arguments : etat, nombre minimum de valeurs,
 *              seuil sur la fonction de repartition.
 *
 *--------------------------------------------------------------*/

void SemiMarkovChain::state_first_occurrence_distribution(int state , int min_nb_value ,
                                                          double cumul_threshold)

{
  register int i , j , k;
  double *state_out , **state_in , *pmass , *pcumul;
  DiscreteParametric *occupancy;
  Distribution *first_occurrence;


  first_occurrence = state_process->first_occurrence[state];
  first_occurrence->complement = state_process->no_occurrence[state];

  pmass = first_occurrence->mass;
  pcumul = first_occurrence->cumul;

  state_out = new double[nb_state];

  state_in = new double*[first_occurrence->alloc_nb_value];
  state_in[0] = NULL;
  for (i = 1;i < first_occurrence->alloc_nb_value;i++) {
    state_in[i] = new double[nb_state];
  }

  *pmass = initial[state];
  *pcumul = *pmass;

  i = 1;

  while (((*pcumul < cumul_threshold - first_occurrence->complement) || (i < min_nb_value)) &&
         (i < first_occurrence->alloc_nb_value)) {

    // calcul des probabilites de quitter (semi-Markov) / d'etre dans (Markov) un etat et
    // de la valeur courante

    *++pmass = 0.;

    for (j = 0;j < nb_state;j++) {
      if (j != state) {
        switch (state_subtype[j]) {

        // cas etat semi-markovien

        case SEMI_MARKOVIAN : {
          occupancy = state_process->sojourn_time[j];
          state_out[j] = 0.;

          for (k = 1;k <= MIN(i , occupancy->nb_value - 1);k++) {
            if (k < i) {
              state_out[j] += occupancy->mass[k] * state_in[i - k][j];
            }
            else {
              switch (type) {
              case 'o' :
                state_out[j] += occupancy->mass[k] * initial[j];
                break;
              case 'e' :
                state_out[j] += forward[j]->mass[k] * initial[j];
                break;
              }
            }
          }
          break;
        }

        // cas etat markovien

        case MARKOVIAN : {
          if (i == 1) {
            state_out[j] = initial[j];
          }
          else {
            state_out[j] = state_in[i - 1][j];
          }
          break;
        }
        }

        *pmass += transition[j][state] * state_out[j];
      }
    }

    for (j = 0;j < nb_state;j++) {
      if (j != state) {
        state_in[i][j] = 0.;
        for (k = 0;k < nb_state;k++) {
          if (k != state) {
            state_in[i][j] += transition[k][j] * state_out[k];
          }
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

  delete [] state_out;

  for (i = 1;i < first_occurrence->alloc_nb_value;i++) {
    delete [] state_in[i];
  }
  delete [] state_in;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la probabilite de quitter un etat sans possibilite
 *  d'y revenir pour une semi-chaine de Markov ordinaire.
 *
 *  arguments : etat, seuil sur la somme des probabilites de quitter un etat.
 *
 *--------------------------------------------------------------*/

void SemiMarkovChain::state_leave_probability(int state , double increment)

{
  if (state_type[state] == 't') {
    register int i , j , k;
    int min_time;
    double sum , *state_out , **state_in , &leave = state_process->leave[state];
    DiscreteParametric *occupancy;


    state_out = new double[nb_state];

    state_in = new double*[LEAVE_LENGTH];
    state_in[0] = NULL;
    state_in[1] = NULL;
    for (i = 2;i < LEAVE_LENGTH;i++) {
      state_in[i] = new double[nb_state];
    }

    leave = 0.;
    for (i = 0;i < nb_state;i++) {
      if ((i != state) && (!accessibility[i][state])) {
        leave += transition[state][i];
      }
    }

    sum = 0.;
    for (i = 0;i < nb_state;i++) {
      if (i != state) {
        switch (state_subtype[i]) {
        case SEMI_MARKOVIAN :
          sum += state_process->sojourn_time[i]->mean;
          break;
        case MARKOVIAN :
          sum += 1. / (1. - transition[i][i]);
          break;
        }
      }
    }
    min_time = (int)sum + 1;

    i = 2;

    do {

      // calcul des probabilites de quitter (semi-Markov) / d'etre dans (Markov) un etat et
      // mise a jour de la probabilite de quitter l'etat selectionne

      sum = 0.;

      for (j = 0;j < nb_state;j++) {
        if ((j != state) && (accessibility[j][state])) {
          switch (state_subtype[j]) {

          // cas etat semi-markovien

          case SEMI_MARKOVIAN : {
            occupancy = state_process->sojourn_time[j];
            state_out[j] = 0.;

            for (k = 1;k < MIN(i , occupancy->nb_value);k++) {
              if (k < i - 1) {
                state_out[j] += occupancy->mass[k] * state_in[i - k][j];
              }
              else {
                state_out[j] += occupancy->mass[k] * transition[state][j];
              }
            }
            break;
          }

          // cas etat markovien

          case MARKOVIAN : {
            if (i == 2) {
              state_out[j] = transition[state][j];
            }
            else {
              state_out[j] = state_in[i - 1][j];
            }
            break;
          }
          }

          switch (state_subtype[j]) {
          case SEMI_MARKOVIAN :
            sum += state_out[j];
            break;
          case MARKOVIAN :
            sum += state_out[j] * (1. - transition[j][j]);
            break;
          }

          for (k = 0;k < nb_state;k++) {
            if ((k != state) && (!accessibility[k][state])) {
              leave += transition[j][k] * state_out[j];
            }
          }
        }
      }

      if (transition[state][state] > 0.) {
        sum /= (1. - transition[state][state]);
      }

      for (j = 0;j < nb_state;j++) {
        if ((j != state) && (accessibility[j][state])) {
          state_in[i][j] = 0.;
          for (k = 0;k < nb_state;k++) {
            if ((k != state) && (accessibility[k][state])) {
              state_in[i][j] += transition[k][j] * state_out[k];
            }
          }
        }
      }

      i++;
    }
    while (((sum > increment) || (i <= min_time)) && (i < LEAVE_LENGTH));

    if (state_subtype[state] == SEMI_MARKOVIAN) {
      leave /= state_process->sojourn_time[state]->parametric_mean_computation();
    }

    delete [] state_out;

    for (i = 2;i < LEAVE_LENGTH;i++) {
      delete [] state_in[i];
    }
    delete [] state_in;
  }
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la loi du temps de retour dans un etat
 *  pour une semi-chaine de Markov.
 *
 *  arguments : etat, nombre minimum de valeurs,
 *              seuil sur la fonction de repartition.
 *
 *--------------------------------------------------------------*/

void SemiMarkovChain::state_recurrence_time_distribution(int state , int min_nb_value ,
                                                         double cumul_threshold)

{
  register int i , j , k;
  double occupancy_mean , *state_out , **state_in , *pmass , *pcumul;
  Distribution *recurrence_time;
  DiscreteParametric *occupancy;


  recurrence_time = state_process->recurrence_time[state];
  recurrence_time->complement = state_process->leave[state];

  pmass = recurrence_time->mass;
  pcumul = recurrence_time->cumul;
  *pmass = 0.;
  *pcumul = 0.;

  state_out = new double[nb_state];

  state_in = new double*[recurrence_time->alloc_nb_value];
  state_in[0] = NULL;
  state_in[1] = NULL;
  for (i = 2;i < recurrence_time->alloc_nb_value;i++) {
    state_in[i] = new double[nb_state];
  }

  // calcul de la probabilite de la valeur 1

  switch (state_subtype[state]) {
  case SEMI_MARKOVIAN :
    occupancy_mean = state_process->sojourn_time[state]->parametric_mean_computation();
    *++pmass = (occupancy_mean - 1.) / occupancy_mean;
    break;
  case MARKOVIAN :
    *++pmass = transition[state][state];
    break;
  }

  *++pcumul = *pmass;

  i = 2;

  while (((*pcumul < cumul_threshold - recurrence_time->complement) || (i < min_nb_value)) &&
         (i < recurrence_time->alloc_nb_value)) {

    // calcul des probabilites de quitter (semi-Markov) / d'etre dans (Markov) un etat et
    // de la valeur courante

    *++pmass = 0.;

    for (j = 0;j < nb_state;j++) {
      if (j != state) {
        switch (state_subtype[j]) {

        // cas etat semi-markovien

        case SEMI_MARKOVIAN : {
          occupancy = state_process->sojourn_time[j];
          state_out[j] = 0.;

          for (k = 1;k < MIN(i , occupancy->nb_value);k++) {
            if (k < i - 1) {
              state_out[j] += occupancy->mass[k] * state_in[i - k][j];
            }
            else {
              state_out[j] += occupancy->mass[k] * transition[state][j];
            }
          }
          break;
        }

        // cas etat markovien

        case MARKOVIAN : {
          if (i == 2) {
            state_out[j] = transition[state][j];
          }
          else {
            state_out[j] = state_in[i - 1][j];
          }
          break;
        }
        }

        *pmass += transition[j][state] * state_out[j];
      }
    }

    for (j = 0;j < nb_state;j++) {
      if (j != state) {
        state_in[i][j] = 0.;
        for (k = 0;k < nb_state;k++) {
          if (k != state) {
            state_in[i][j] += transition[k][j] * state_out[k];
          }
        }
      }
    }

    if (state_subtype[state] == SEMI_MARKOVIAN) {
      *pmass /= occupancy_mean;
    }
    pcumul++;
    *pcumul = *(pcumul - 1) + *pmass;
    i++;
  }

  recurrence_time->nb_value = i;
  recurrence_time->nb_value_computation();

  delete [] state_out;

  for (i = 2;i < recurrence_time->alloc_nb_value;i++) {
    delete [] state_in[i];
  }
  delete [] state_in;

  if (recurrence_time->nb_value > 0) {
    recurrence_time->offset_computation();
    recurrence_time->max_computation();
    recurrence_time->mean_computation();
    recurrence_time->variance_computation();
  }

  else {
    delete state_process->recurrence_time[state];
    state_process->recurrence_time[state] = NULL;
  }
}


/*--------------------------------------------------------------*
 *
 *  Calcul des probabilites de chaque observation en fonction de l'index
 *  pour une semi-chaine de Markov cachee.
 *
 *  argument : indice du processus d'observation.
 *
 *--------------------------------------------------------------*/

void SemiMarkov::index_output_distribution(int variable)

{
  register int i , j , k;
  Curves *index_state , *index_value;


  index_value = categorical_process[variable]->index_value;

  // calcul des probabilites des etats de la semi-chaine de Markov
  // sous-jacente en fonction de l'index si necessaire

  if (!(state_process->index_value)) {
    state_process->index_value = new Curves(nb_state , index_value->length);
    index_state_distribution();
  }
  index_state = state_process->index_value;

  // prise en compte des probabilites d'observation

  for (i = 0;i < index_value->length;i++) {
    for (j = 0;j < categorical_process[variable]->nb_value;j++) {
      index_value->point[j][i] = 0.;
      for (k = 0;k < nb_state;k++) {
        index_value->point[j][i] += categorical_process[variable]->observation[k]->mass[j] *
                                    index_state->point[k][i];
      }
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la probabilite de ne pas observer une valeur
 *  pour une semi-chaine de Markov cachee ordinaire.
 *
 *  arguments : indice du processus d'observation, observation,
 *              seuils sur la somme des probabilites de quitter un etat.
 *
 *--------------------------------------------------------------*/

void SemiMarkov::output_no_occurrence_probability(int variable , int output ,
                                                  double increment)

{
  bool status = false , *output_accessibility;
  register int i , j , k;
  int min_time;
  double sum , *state_out , **state_in , *observation , **obs_power ,
         &no_occurrence = categorical_process[variable]->no_occurrence[output];
  DiscreteParametric *occupancy;


  observation = new double[nb_state];
  for (i = 0;i < nb_state;i++) {
    observation[i] = categorical_process[variable]->observation[i]->mass[output];
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
    obs_power = new double*[nb_state];
    for (i = 0;i < nb_state;i++) {
      if (state_subtype[i] == SEMI_MARKOVIAN) {
        obs_power[i] = new double[LEAVE_LENGTH + 1];
        obs_power[i][0] = 1.;
      }
    }

    state_out = new double[nb_state];

    state_in = new double*[LEAVE_LENGTH];
    for (i = 0;i < LEAVE_LENGTH;i++) {
      state_in[i] = new double[nb_state];
    }

    no_occurrence = 0.;
    for (i = 0;i < nb_state;i++) {
      if (!output_accessibility[i]) {
        no_occurrence += initial[i];
      }
    }

    sum = 0.;
    for (i = 0;i < nb_state;i++) {
      switch (state_subtype[i]) {

      case SEMI_MARKOVIAN : {
        sum += state_process->sojourn_time[i]->mean;
        break;
      }

      case MARKOVIAN : {
        if (transition[i][i] < 1.) {
          sum += 1. / (1. - transition[i][i]);
        }
        break;
      }
      }
    }
    min_time = (int)sum + 1;

    i = 0;

    do {

      // calcul des probabilites de quitter (semi-Markov) / d'etre dans (Markov) un etat et
      // mise a jour de la probabilite de ne pas observer la valeur

      sum = 0.;

      for (j = 0;j < nb_state;j++) {
        if (output_accessibility[j]) {
          switch (state_subtype[j]) {

          // cas etat semi-markovien

          case SEMI_MARKOVIAN : {
            occupancy = state_process->sojourn_time[j];
            state_out[j] = 0.;

            // calcul des puissances des probabilites d'observation

            obs_power[j][i + 1] = obs_power[j][i] * (1. - observation[j]);

            for (k = 1;k <= MIN(i + 1 , occupancy->nb_value - 1);k++) {
              if (k < i + 1) {
                state_out[j] += obs_power[j][k] * occupancy->mass[k] * state_in[i - k][j];
              }
              else {
                state_out[j] += obs_power[j][k] * occupancy->mass[k] * initial[j];
              }
            }
            break;
          }

          // cas etat markovien

          case MARKOVIAN : {
            if (i == 0) {
              state_out[j] = (1. - observation[j]) * initial[j];
            }
            else {
              state_out[j] = (1. - observation[j]) * state_in[i - 1][j];
            }
            break;
          }
          }

          if ((transition[j][j] == 0.) || (transition[j][j] == 1.)) {
            sum += state_out[j];
          }
          else {
            sum += state_out[j] * (1. - transition[j][j]);
          }

          for (k = 0;k < nb_state;k++) {
            if (!output_accessibility[k]) {
              no_occurrence += transition[j][k] * state_out[j];
            }
          }
        }
      }

      for (j = 0;j < nb_state;j++) {
        if (output_accessibility[j]) {
          state_in[i][j] = 0.;
          for (k = 0;k < nb_state;k++) {
            if (output_accessibility[k]) {
              state_in[i][j] += transition[k][j] * state_out[k];
            }
          }
        }
      }

      i++;
    }
    while (((sum > increment) || (i < min_time)) && (i < LEAVE_LENGTH));

    for (i = 0;i < nb_state;i++) {
      if (state_subtype[i] == SEMI_MARKOVIAN) {
        delete [] obs_power[i];
      }
    }
    delete [] obs_power;

    delete [] state_out;

    for (i = 0;i < LEAVE_LENGTH;i++) {
      delete [] state_in[i];
    }
    delete [] state_in;
  }

  delete [] observation;
  delete [] output_accessibility;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la loi du temps avant la 1ere occurrence d'une observation
 *  pour une semi-chaine de Markov cachee.
 *
 *  arguments : indice du processus d'observation, observation,
 *              nombre minimum de valeurs, seuil sur la fonction de repartition.
 *
 *--------------------------------------------------------------*/

void SemiMarkov::output_first_occurrence_distribution(int variable , int output ,
                                                      int min_nb_value , double cumul_threshold)

{
  register int i , j , k;
  double sum , *state_out , **state_in , *observation , **obs_power , *pmass , *pcumul;
  DiscreteParametric *occupancy;
  Distribution *first_occurrence;


  first_occurrence = categorical_process[variable]->first_occurrence[output];
  first_occurrence->complement = categorical_process[variable]->no_occurrence[output];

  pmass = first_occurrence->mass - 1;
  pcumul = first_occurrence->cumul - 1;

  observation = new double[nb_state];
  for (i = 0;i < nb_state;i++) {
    observation[i] = categorical_process[variable]->observation[i]->mass[output];
  }

  obs_power = new double*[nb_state];
  for (i = 0;i < nb_state;i++) {
    if (state_subtype[i] == SEMI_MARKOVIAN) {
      obs_power[i] = new double[first_occurrence->alloc_nb_value + 1];
      obs_power[i][0] = 1.;
    }
  }

  state_out = new double[nb_state];

  state_in = new double*[first_occurrence->alloc_nb_value];
  for (i = 0;i < first_occurrence->alloc_nb_value;i++) {
    state_in[i] = new double[nb_state];
  }

  i = 0;

  do {

    // calcul des probabilites de quitter (semi-Markov) / d'etre dans (Markov) un etat et
    // de la valeur courante

    *++pmass = 0.;

    for (j = 0;j < nb_state;j++) {
      switch (state_subtype[j]) {

      // cas etat semi-markovien

      case SEMI_MARKOVIAN : {
        occupancy = state_process->sojourn_time[j];
        state_out[j] = 0.;
        sum = 0.;

        // calcul des puissances des probabilites d'observation

        obs_power[j][i + 1] = obs_power[j][i] * (1. - observation[j]);

        for (k = 1;k <= MIN(i + 1 , occupancy->nb_value - 1);k++) {
          if (k < i + 1) {
            state_out[j] += obs_power[j][k] * occupancy->mass[k] * state_in[i - k][j];
            sum += obs_power[j][k - 1] * (1. - occupancy->cumul[k - 1]) * state_in[i - k][j];
          }
          else {
            switch (type) {
            case 'o' :
              state_out[j] += obs_power[j][k] * occupancy->mass[k] * initial[j];
              sum += obs_power[j][k - 1] * (1. - occupancy->cumul[k - 1]) * initial[j];
              break;
            case 'e' :
              state_out[j] += obs_power[j][k] * forward[j]->mass[k] * initial[j];
              sum += obs_power[j][k - 1] * (1. - forward[j]->cumul[k - 1]) * initial[j];
              break;
            }
          }
        }
        break;
      }

      // cas etat markovien

      case MARKOVIAN : {
        if (i == 0) {
          state_out[j] = (1. - observation[j]) * initial[j];
          sum = initial[j];
        }
        else {
          state_out[j] = (1. - observation[j]) * state_in[i - 1][j];
          sum = state_in[i - 1][j];
        }
        break;
      }
      }

      *pmass += observation[j] * sum;
    }

    for (j = 0;j < nb_state;j++) {
      state_in[i][j] = 0.;
      for (k = 0;k < nb_state;k++) {
        state_in[i][j] += transition[k][j] * state_out[k];
      }
    }

    // mise a jour de la fonction de repartition

    pcumul++;
    if (i == 0) {
      *pcumul = *pmass;
    }
    else {
      *pcumul = *(pcumul - 1) + *pmass;
    }
    i++;
  }
  while (((*pcumul < cumul_threshold - first_occurrence->complement) || (i < min_nb_value)) &&
         (i < first_occurrence->alloc_nb_value));

  first_occurrence->nb_value = i;

  first_occurrence->offset_computation();
  first_occurrence->max_computation();
  first_occurrence->mean_computation();
  first_occurrence->variance_computation();

  delete [] observation;

  for (i = 0;i < nb_state;i++) {
    if (state_subtype[i] == SEMI_MARKOVIAN) {
      delete [] obs_power[i];
    }
  }
  delete [] obs_power;

  delete [] state_out;

  for (i = 0;i < first_occurrence->alloc_nb_value;i++) {
    delete [] state_in[i];
  }
  delete [] state_in;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la probabilite de quitter une observation sans possibilite
 *  d'y revenir pour une semi-chaine de Markov cachee ordinaire.
 *
 *  arguments : loi des memoires, indice du processus d'observation, observation,
 *              seuil sur la somme des probabilites de quitter un etat.
 *
 *--------------------------------------------------------------*/

void SemiMarkov::output_leave_probability(const double *memory , int variable ,
                                          int output , double increment)

{
  bool status = false , *output_accessibility;
  register int i , j , k;
  int min_time;
  double sum0 , sum1 , *observation , **obs_power , *input_proba , *state_out ,
         **state_in , &leave = categorical_process[variable]->leave[output];
  DiscreteParametric *occupancy;


  observation = new double[nb_state];
  for (i = 0;i < nb_state;i++) {
    observation[i] = categorical_process[variable]->observation[i]->mass[output];
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
    obs_power = new double*[nb_state];
    for (i = 0;i < nb_state;i++) {
      if (transition[i][i] < 1.) {
        obs_power[i] = new double[LEAVE_LENGTH];
        obs_power[i][0] = 1.;
      }
    }

    state_out = new double[nb_state];

    state_in = new double*[LEAVE_LENGTH];
    state_in[0] = NULL;
    for (i = 1;i < LEAVE_LENGTH;i++) {
      state_in[i] = new double[nb_state];
    }

    // calcul des probabilites d'entree et de sortie

    input_proba = new double[nb_state];
    sum0 = 0.;

    for (i = 0;i < nb_state;i++) {
      sum1 = 0.;
      for (j = 0;j < nb_state;j++) {
        if ((i != j) || (transition[j][j] == 1.)) {
           sum1 += transition[j][i] * memory[j];
        }
      }
      input_proba[i] = observation[i] * (initial[i] + sum1);

      // cas etat non-absorbant

      if (transition[i][i] < 1.) {
        switch (state_subtype[i]) {
        case SEMI_MARKOVIAN :
          sum0 += state_process->sojourn_time[i]->mean * input_proba[i];
          break;
        case MARKOVIAN :
          sum0 += input_proba[i] / (1. - transition[i][i]);
          break;
        }
      }

      // cas etat absorbant

      else {
        sum0 += input_proba[i];
      }
    }

    for (i = 0;i < nb_state;i++) {
      input_proba[i] /= sum0;
    }

    sum0 = 0.;
    for (i = 0;i < nb_state;i++) {
      switch (state_subtype[i]) {

      case SEMI_MARKOVIAN : {
        sum0 += state_process->sojourn_time[i]->mean;
        break;
      }

      case MARKOVIAN : {
        if (transition[i][i] < 1.) {
          sum0 += 1. / (1. - transition[i][i]);
        }
        break;
      }
      }
    }
    min_time = (int)sum0 + 1;

    leave = 0.;
    i = 1;

    do {

      // calcul des probabilites de quitter un etat et mise a jour
      // de la probabilite de quitter l'observation selectionnee

      sum0 = 0.;

      for (j = 0;j < nb_state;j++) {
        if (output_accessibility[j]) {
          state_out[j] = 0.;

          // cas etat non-absorbant

          if (transition[j][j] < 1.) {
            occupancy = state_process->sojourn_time[j];

            // calcul des puissances des probabilites d'observation

            obs_power[j][i] = obs_power[j][i - 1] * (1. - observation[j]);

            for (k = 1;k <= MIN(i , occupancy->nb_value - 1);k++) {
              if (k < i) {
                state_out[j] += obs_power[j][k] * occupancy->mass[k] * state_in[i - k][j];
              }
              else {
                state_out[j] += obs_power[j][k - 1] * (1. - occupancy->cumul[k - 1]) *
                                input_proba[j];
              }
            }

            sum0 += state_out[j];

            switch (state_subtype[j]) {

            case SEMI_MARKOVIAN : {
              for (k = 0;k < nb_state;k++) {
                if (!output_accessibility[k]) {
                  leave += transition[j][k] * state_out[j];
                }
              }
              break;
            }

            case MARKOVIAN : {
              for (k = 0;k < nb_state;k++) {
                if ((!output_accessibility[k]) && (k != j)) {
                  leave += transition[j][k] * state_out[j] / (1. - transition[j][j]);
                }
              }
              break;
            }
            }
          }
        }
      }

      for (j = 0;j < nb_state;j++) {
        if (output_accessibility[j]) {
          state_in[i][j] = 0.;
          for (k = 0;k < nb_state;k++) {
            if (output_accessibility[k]) {
              if ((transition[k][k] == 0.) || (transition[k][k] == 1.)) {
                state_in[i][j] += transition[k][j] * state_out[k];
              }
              else if (j != k) {
                state_in[i][j] += transition[k][j] * state_out[k] / (1. - transition[k][k]);
              }
            }
          }
        }
      }

      i++;
    }
    while (((sum0 > increment) || (i <= min_time)) && (i < LEAVE_LENGTH));

    for (i = 0;i < nb_state;i++) {
      if (transition[i][i] < 1.) {
        delete [] obs_power[i];
      }
    }
    delete [] obs_power;

    delete [] state_out;

    for (i = 1;i < LEAVE_LENGTH;i++) {
      delete [] state_in[i];
    }
    delete [] state_in;

    delete [] input_proba;
  }

  delete [] observation;
  delete [] output_accessibility;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la loi du temps de retour dans une observation
 *  pour une semi-chaine de Markov cachee.
 *
 *  arguments : loi des memoires, indice du processus d'observation,
 *              observation, nombre minimum de valeurs,
 *              seuil sur la fonction de repartition.
 *
 *--------------------------------------------------------------*/

void SemiMarkov::output_recurrence_time_distribution(const double *memory , int variable ,
                                                     int output , int min_nb_value ,
                                                     double cumul_threshold)

{
  register int i , j , k , m;
  double sum0 , sum1 , *observation , **obs_power , *input_proba , *output_proba ,
         *state_out , **state_in , *pmass , *pcumul;
  Distribution *recurrence_time;
  DiscreteParametric *occupancy;


  recurrence_time = categorical_process[variable]->recurrence_time[output];
  recurrence_time->complement = categorical_process[variable]->leave[output];

  pmass = recurrence_time->mass;
  pcumul = recurrence_time->cumul;
  *pmass = 0.;
  *pcumul = 0.;

  observation = new double[nb_state];
  for (i = 0;i < nb_state;i++) {
    observation[i] = categorical_process[variable]->observation[i]->mass[output];
  }

  obs_power = new double*[nb_state];
  for (i = 0;i < nb_state;i++) {
    if (transition[i][i] < 1.) {
      obs_power[i] = new double[recurrence_time->alloc_nb_value];
      obs_power[i][0] = 1.;
    }
  }

  state_out = new double[nb_state];

  state_in = new double*[recurrence_time->alloc_nb_value];
  state_in[0] = NULL;
  for (i = 1;i < recurrence_time->alloc_nb_value;i++) {
    state_in[i] = new double[nb_state];
  }

  // calcul des probabilites d'entree et de sortie

  input_proba = new double[nb_state];
  output_proba = new double[nb_state];
  sum0 = 0.;

  for (i = 0;i < nb_state;i++) {
    sum1 = 0.;
    for (j = 0;j < nb_state;j++) {
      if ((i != j) || (transition[j][j] == 1.)) {
        sum1 += transition[j][i] * memory[j];
      }
    }
    input_proba[i] = observation[i] * (initial[i] + sum1);

    // cas etat non-absorbant

    if (transition[i][i] < 1.) {
      switch (state_subtype[i]) {
      case SEMI_MARKOVIAN :
        sum0 += state_process->sojourn_time[i]->mean * input_proba[i];
        break;
      case MARKOVIAN :
        sum0 += input_proba[i] / (1. - transition[i][i]);
        break;
      }

      sum1 = 0.;

      switch (state_subtype[i]) {

      case SEMI_MARKOVIAN : {
        for (j = 0;j < nb_state;j++) {
          sum1 += observation[j] * transition[i][j];
        }
        break;
      }

      case MARKOVIAN : {
        for (j = 0;j < nb_state;j++) {
          if (j != i) {
            sum1 += observation[j] * transition[i][j] / (1. - transition[i][i]);
          }
        }
        break;
      }
      }

      output_proba[i] = sum1;
    }

    // cas etat absorbant

    else {
      sum0 += input_proba[i];
    }
  }

  for (i = 0;i < nb_state;i++) {
    input_proba[i] /= sum0;
  }

  i = 1;

  do {

    // calcul des probabilites de quitter un etat et de la valeur courante

    *++pmass = 0.;

    for (j = 0;j < nb_state;j++) {

      // cas etat non-absorbant

      if (transition[j][j] < 1.) {
        occupancy = state_process->sojourn_time[j];
        state_out[j] = 0.;
        sum0 = 0.;

        // calcul des puissances des probabilites d'observation

        obs_power[j][i] = obs_power[j][i - 1] * (1. - observation[j]);

        for (k = 1;k <= MIN(i , occupancy->nb_value - 1);k++) {
          if (k < i) {
            state_out[j] += obs_power[j][k] * occupancy->mass[k] * state_in[i - k][j];
            sum0 += obs_power[j][k] * (1. - occupancy->cumul[k]) * state_in[i - k][j];
          }

          else {
            state_out[j] += obs_power[j][k - 1] * (1. - occupancy->cumul[k - 1]) * input_proba[j];

            sum1 = 0.;
            for (m = k;m < occupancy->nb_value;m++) {
              sum1 += (1. - occupancy->cumul[m]);
            }
            sum0 += obs_power[j][k - 1] * sum1 * input_proba[j];
          }
        }

        *pmass += output_proba[j] * state_out[j] + observation[j] * sum0;
      }

      // cas etat absorbant

      else {
        if (i == 1) {
          state_out[j] = input_proba[j];
        }
        else {
          state_out[j] = (1. - observation[j]) * state_in[i - 1][j];
        }

        *pmass += observation[j] * state_out[j];
      }
    }

    for (j = 0;j < nb_state;j++) {
      state_in[i][j] = 0.;
      for (k = 0;k < nb_state;k++) {
        if ((transition[k][k] == 0.) || (transition[k][k] == 1.)) {
          state_in[i][j] += transition[k][j] * state_out[k];
        }
        else if (j != k) {
          state_in[i][j] += transition[k][j] * state_out[k] / (1. - transition[k][k]);
        }
      }
    }

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
    delete categorical_process[variable]->recurrence_time[output];
    categorical_process[variable]->recurrence_time[output] = NULL;
  }

  delete [] observation;

  for (i = 0;i < nb_state;i++) {
    if (transition[i][i] < 1.) {
      delete [] obs_power[i];
    }
  }
  delete [] obs_power;

  delete [] state_out;

  for (i = 1;i < recurrence_time->alloc_nb_value;i++) {
    delete [] state_in[i];
  }
  delete [] state_in;

  delete [] input_proba;
  delete [] output_proba;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la loi du temps de sejour dans une observation
 *  pour une semi-chaine de Markov cachee.
 *
 *  arguments : loi des memoires, indice du processus d'observation,
 *              observation, nombre minimum de valeurs,
 *              seuil sur la fonction de repartition.
 *
 *--------------------------------------------------------------*/

void SemiMarkov::output_sojourn_time_distribution(const double *memory , int variable ,
                                                  int output , int min_nb_value ,
                                                  double cumul_threshold)

{
  register int i , j , k , m;
  double sum0 , sum1 , *observation , **obs_power , **input_proba ,
         *output_proba , *state_out , **state_in , *pmass , *pcumul ,
         &absorption = categorical_process[variable]->absorption[output];
  DiscreteParametric *sojourn_time , *occupancy;


  sojourn_time = categorical_process[variable]->sojourn_time[output];

  pmass = sojourn_time->mass;
  pcumul = sojourn_time->cumul;
  *pmass = 0.;
  *pcumul = 0.;

  observation = new double[nb_state];
  for (i = 0;i < nb_state;i++) {
    observation[i] = categorical_process[variable]->observation[i]->mass[output];
  }

  obs_power = new double*[nb_state];
  for (i = 0;i < nb_state;i++) {
    if (transition[i][i] < 1.) {
      obs_power[i] = new double[sojourn_time->alloc_nb_value];
      obs_power[i][0] = 1.;
    }
  }

  state_out = new double[nb_state];

  state_in = new double*[sojourn_time->alloc_nb_value];
  state_in[0] = NULL;
  for (i = 1;i < sojourn_time->alloc_nb_value;i++) {
    state_in[i] = new double[nb_state];
  }

  // calcul des probabilites d'entree et de sortie

  input_proba = new double*[nb_state];
  output_proba = new double[nb_state];
  sum0 = 0.;

  for (i = 0;i < nb_state;i++) {
    input_proba[i] = new double[2];

    sum1 = 0.;
    for (j = 0;j < nb_state;j++) {
      if ((i != j) || (transition[j][j] == 1.)) {
        sum1 += transition[j][i] * (1. - observation[j]) * memory[j];
      }
    }
    input_proba[i][0] = observation[i] * (initial[i] + sum1);
    sum0 += input_proba[i][0];

    // cas etat non-absorbant

    if (transition[i][i] < 1.) {
      sum1 = 0.;
      for (j = 0;j < nb_state;j++) {
        if ((i != j) || (transition[j][j] == 1.)) {
          sum1 += transition[j][i] * memory[j];
        }
      }
      input_proba[i][1] = observation[i] * (1. - observation[i]) * (initial[i] + sum1);

      switch (state_subtype[i]) {
      case SEMI_MARKOVIAN :
        sum0 += (state_process->sojourn_time[i]->mean - 1) * input_proba[i][1];
        break;
      case MARKOVIAN :
        sum0 += (1. / (1. - transition[i][i]) - 1) * input_proba[i][1];
        break;
      }

      sum1 = 0.;

      switch (state_subtype[i]) {

      case SEMI_MARKOVIAN : {
        for (j = 0;j < nb_state;j++) {
          sum1 += (1. - observation[j]) * transition[i][j];
        }
        break;
      }

      case MARKOVIAN : {
        for (j = 0;j < nb_state;j++) {
          if (j != i) {
            sum1 += (1. - observation[j]) * transition[i][j] / (1. - transition[i][i]);
          }
        }
        break;
      }
      }

      output_proba[i] = sum1;
    }
  }

  for (i = 0;i < nb_state;i++) {
    input_proba[i][0] /= sum0;
    if (transition[i][i] < 1.) {
      input_proba[i][1] /= sum0;
    }
  }

  i = 1;

  do {

    // calcul des probabilites de quitter un etat et de la valeur courante

    absorption = 0.;
    *++pmass = 0.;

    for (j = 0;j < nb_state;j++) {
      state_out[j] = 0.;

      if (observation[j] > 0.) {

        // cas etat non-absorbant

        if (transition[j][j] < 1.) {
          occupancy = state_process->sojourn_time[j];
          sum0 = 0.;

          // calcul des puissances des probabilites d'observation

          obs_power[j][i] = obs_power[j][i - 1] * observation[j];

          for (k = 1;k <= MIN(i , occupancy->nb_value - 1);k++) {
            if (k < i) {
              state_out[j] += obs_power[j][k] * occupancy->mass[k] * state_in[i - k][j];
              sum0 += obs_power[j][k] * (1. - occupancy->cumul[k]) * state_in[i - k][j];
            }

            else {
              state_out[j] += obs_power[j][k - 1] * (occupancy->mass[k] * input_proba[j][0] +
                               (1. - occupancy->cumul[k]) * input_proba[j][1]);

              sum1 = 0.;
              for (m = k + 1;m < occupancy->nb_value;m++) {
                sum1 += (1. - occupancy->cumul[m]);
              }
              sum0 += obs_power[j][k - 1] * ((1. - occupancy->cumul[k]) * input_proba[j][0] +
                       sum1 * input_proba[j][1]);
            }
          }

          *pmass += output_proba[j] * state_out[j] + (1. - observation[j]) * sum0;
        }

        // cas etat absorbant

        else {
          if (i == 1) {
            state_out[j] = input_proba[j][0];
          }
          else {
            state_out[j] = observation[j] * state_in[i - 1][j];
          }

          *pmass += (1. - observation[j]) * state_out[j];
        }

        if ((transition[j][j] == 0.) || (transition[j][j] == 1.)) {
          for (k = 0;k < nb_state;k++) {
            if ((state_type[k] == 'a') && (observation[k] == 1.)) {
              absorption += transition[j][k] * state_out[j];
            }
          }
        }

        else {
          for (k = 0;k < nb_state;k++) {
            if ((state_type[k] == 'a') && (observation[k] == 1.) && (k != j)) {
              absorption += transition[j][k] * state_out[j] / (1. - transition[j][j]);
            }
          }
        }
      }
    }

    for (j = 0;j < nb_state;j++) {
      state_in[i][j] = 0.;
      for (k = 0;k < nb_state;k++) {
        if ((transition[k][k] == 0.) || (transition[k][k] == 1.)) {
          state_in[i][j] += transition[k][j] * state_out[k];
        }
        else if (j != k) {
          state_in[i][j] += transition[k][j] * state_out[k] / (1. - transition[k][k]);
        }
      }
    }

    pcumul++;
    *pcumul = *(pcumul - 1) + *pmass;
    i++;
  }
  while (((*pcumul < cumul_threshold - absorption) || (i < min_nb_value)) &&
         (i < sojourn_time->alloc_nb_value));

  if (*pcumul > 1.) {

#   ifdef MESSAGE
    cout << STAT_label[STATL_OUTPUT] << " " << output << ": CONVERGENCE ERROR" << endl;
#   endif

  }

  if (*pcumul == 0.) {
    absorption = 1.;
    delete categorical_process[variable]->sojourn_time[output];
    categorical_process[variable]->sojourn_time[output] = NULL;
  }

  else {
    sojourn_time->nb_value = i;
    sojourn_time->complement = absorption;

    sojourn_time->offset_computation();
    sojourn_time->max_computation();
    sojourn_time->mean_computation();
    sojourn_time->variance_computation();
  }

  delete [] observation;

  for (i = 0;i < nb_state;i++) {
    if (transition[i][i] < 1.) {
      delete [] obs_power[i];
    }
  }
  delete [] obs_power;

  delete [] state_out;

  for (i = 1;i < sojourn_time->alloc_nb_value;i++) {
    delete [] state_in[i];
  }
  delete [] state_in;

  for (i = 0;i < nb_state;i++) {
    delete [] input_proba[i];
  }
  delete [] input_proba;

  delete [] output_proba;
}


};  // namespace sequence_analysis
