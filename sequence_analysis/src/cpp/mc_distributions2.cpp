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

using namespace std;



/*--------------------------------------------------------------*
 *
 *  Calcul d'un melange de lois du nombre de series d'un etat ('r') ou
 *  du nombre d'occurrences d'un etat ('o') d'une chaine de Markov
 *  pour une distribution des longueurs de sequences donnee.
 *
 *  arguments : etat, type de forme.
 *
 *--------------------------------------------------------------*/

void Markov::state_nb_pattern_mixture(int state , char pattern)

{
  register int i , j , k , m;
  int max_length , nb_pattern , index_nb_pattern , increment , power[ORDER] ,
      state_index[ORDER];
  double sum , *state_seq , *pstate_seq , *states , *pstates , *pmass , *lmass ,
         **ptransition;
  Distribution *pdist;


  max_length = process[0]->length->nb_value - 1;

  switch (pattern) {
  case 'r' :
    pdist = process[0]->nb_run[state];
    nb_pattern = max_length / 2 + 2;
    break;
  case 'o' :
    pdist = process[0]->nb_occurrence[state];
    nb_pattern = max_length + 1;
    break;
  }

  pmass = pdist->mass;
  for (i = 0;i < pdist->nb_value;i++) {
    *pmass++ = 0.;
  }

  i = 1;
  for (j = 0;j < order;j++) {
    power[j] = i;
    i *= nb_state;
  }

  state_seq = new double[nb_row * nb_pattern];
  pstate_seq = new double[nb_row * nb_pattern];

  lmass = process[0]->length->mass;
  index_nb_pattern = 1;

  for (i = 0;i < max_length;i++) {

    // initialisation des probabilites des sequences d'etats
    // pour un nombre de formes donne de l'etat selectionne

    if (i == 0) {
      states = state_seq;
      j = 0;

      for (k = 0;k < nb_row;k++) {
        if (k == self_row[j]) {
          if (k / power[order - 1] == state) {
            *states++ = 0.;
            *states = initial[j++];
          }
          else {
            *states++ = initial[j++];
            *states = 0.;
          }
        }
        else {
          *states++ = 0.;
          *states = 0.;
        }
        states += nb_pattern - 1;
      }
    }

    else {

      // mise a jour des probabilites des sequences d'etats

      pstates = pstate_seq;
      states = state_seq;

      for (j = 0;j < nb_row;j++) {
        for (k = 0;k < index_nb_pattern;k++) {
          *pstates++ = *states;
          *states++ = 0.;
        }
        *states = 0.;
        pstates += nb_pattern - index_nb_pattern;
        states += nb_pattern - index_nb_pattern;
      }

      states = state_seq;
      for (j = 0;j < order;j++) {
        state_index[j] = 0;
      }

      for (j = 0;j < nb_row;j++) {

        // calcul des probabilites des sequences d'etats de
        // longueur = ordre de la chaine de Markov pour chaque nombre
        // de formes de l'etat selectionne

        ptransition = transition;
        pstates = pstate_seq;
        for (k = 0;k < order - 1;k++) {
          ptransition += state_index[k] * power[k + 1];
          pstates += state_index[k] * power[k + 1] * nb_pattern;
        }

        for (k = 0;k < nb_state;k++) {
          switch (pattern) {
          case 'r' :
            increment = (((state_index[order - 1] == state) &&
                          (((order == 1) && (k != state)) ||
                           ((order > 1) && (state_index[order - 2] != state)))) ? 1 : 0);
            break;
          case 'o' :
            increment = (state_index[order - 1] == state ? 1 : 0);
            break;
          }

          if (increment == 1) {
            states++;
          }
          for (m = 0;m < index_nb_pattern;m++) {
            *states++ += *(*ptransition + state_index[order - 1]) * *pstates++;
          }
          if (increment == 0) {
            states++;
          }
          states -= (index_nb_pattern + 1);
          ptransition++;
          pstates += nb_pattern - index_nb_pattern;
        }
        states += nb_pattern;

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

    if ((pattern == 'o') || (i % 2 == 0)) {
      index_nb_pattern++;
    }

    // mise a jour du melange de lois du nombre de formes de l'etat selectionne

    if (*++lmass > 0.) {
      pmass = pdist->mass;
      for (j = 0;j < index_nb_pattern;j++) {
        states = state_seq + j;
        sum = 0.;
        for (k = 0;k < nb_row;k++) {
          sum += *states;
          states += nb_pattern;
        }
        *pmass++ += *lmass * sum;
      }
    }
  }

  pdist->nb_value_computation();
  pdist->offset_computation();
  pdist->cumul_computation();

  pdist->max_computation();
  pdist->mean_computation();
  pdist->variance_computation();

  delete [] state_seq;
  delete [] pstate_seq;
}


/*--------------------------------------------------------------*
 *
 *  Calcul d'un melange de lois du nombre de series d'une observation
 *  emise par une chaine de Markov cachee pour une distribution
 *  des longueurs de sequences donnee.
 *
 *  arguments : indice du processus d'observation, observation.
 *
 *--------------------------------------------------------------*/

void Markov::output_nb_run_mixture(int variable , int output)

{
  register int i , j , k , m;
  int max_length , nb_pattern , index_nb_pattern , power[ORDER] , state_index[ORDER];
  double sum , *observation , *state_seq , *pstate_seq , *states , *pstates ,
         *pmass , *lmass , **ptransition;
  Distribution *nb_run;


  nb_run = process[variable]->nb_run[output];

  max_length = process[variable]->length->nb_value - 1;
  nb_pattern = max_length / 2 + 2;

  pmass = nb_run->mass;
  for (i = 0;i < nb_run->nb_value;i++) {
    *pmass++ = 0.;
  }

  observation = new double[nb_state];
  for (i = 0;i < nb_state;i++) {
    observation[i] = process[variable]->observation[i]->mass[output];
  }

  i = 1;
  for (j = 0;j < order;j++) {
    power[j] = i;
    i *= nb_state;
  }

  state_seq = new double[nb_row * nb_pattern * 2];
  pstate_seq = new double[nb_row * nb_pattern * 2];

  lmass = process[variable]->length->mass;
  index_nb_pattern = 1;

  for (i = 0;i < max_length;i++) {

    // initialisation des probabilites des sequences d'etats
    // pour un nombre de series donne de l'observation selectionnee

    if (i == 0) {
      states = state_seq;
      j = 0;

      for (k = 0;k < nb_row;k++) {
        if (k == self_row[j]) {
          *states++ = (1. - observation[j]) * initial[j];
          *states++ = 0.;
          *states++ = 0.;
          *states++ = observation[j] * initial[j];
          j++;
        }
        else {
          for (m = 0;m < 4;m++) {
            *states++ = 0.;
          }
        }
        states += (nb_pattern - 2) * 2;
      }
    }

    else {

      // mise a jour des probabilites des sequences d'etats

      pstates = pstate_seq;
      states = state_seq;

      for (j = 0;j < nb_row;j++) {
        for (k = 0;k < index_nb_pattern * 2;k++) {
          *pstates++ = *states;
          *states++ = 0.;
        }
        for (k = 0;k < 2;k++) {
          *pstates++ = 0.;
          *states++ = 0.;
        }
        pstates += (nb_pattern - (index_nb_pattern + 1)) * 2;
        states += (nb_pattern - (index_nb_pattern + 1)) * 2;
      }

      states = state_seq;
      for (j = 0;j < order;j++) {
        state_index[j] = 0;
      }

      for (j = 0;j < nb_row;j++) {

        // calcul des probabilites des sequences d'etats de
        // longueur = ordre de la chaine de Markov pour chaque nombre
        // de series de l'observation selectionnee

        ptransition = transition;
        pstates = pstate_seq;
        for (k = 0;k < order - 1;k++) {
          ptransition += state_index[k] * power[k + 1];
          pstates += state_index[k] * power[k + 1] * nb_pattern * 2;
        }

        for (k = 0;k < nb_state;k++) {
          for (m = 0;m <= index_nb_pattern;m++) {
            if (m < index_nb_pattern) {
              *states += (1. - observation[state_index[order - 1]]) *
                         *(*ptransition + state_index[order - 1]) *
                         (*pstates + *(pstates + 1));
            }
            states++;
            pstates -= 2;

            if (m > 0) {
              *states += observation[state_index[order - 1]] *
                         *(*ptransition + state_index[order - 1]) *
                         (*pstates + *(pstates + 3));
            }
            states++;
            pstates += 4;
          }
          states -= (index_nb_pattern + 1) * 2;
          ptransition++;
          pstates += (nb_pattern - (index_nb_pattern + 1)) * 2;
        }
        states += nb_pattern * 2;

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

    if (i % 2 == 0) {
      index_nb_pattern++;
    }

    // mise a jour du melange de lois du nombre de series de l'observation selectionnee

    if (*++lmass > 0.) {
      pmass = nb_run->mass;
      for (j = 0;j < index_nb_pattern;j++) {
        states = state_seq + j * 2;
        sum = 0.;
        for (k = 0;k < nb_row;k++) {
          sum += *states + *(states + 1);
          states += nb_pattern * 2;
        }
        *pmass++ += *lmass * sum;
      }
    }
  }

  nb_run->nb_value_computation();
  nb_run->offset_computation();
  nb_run->cumul_computation();

  nb_run->max_computation();
  nb_run->mean_computation();
  nb_run->variance_computation();

  delete [] observation;
  delete [] state_seq;
  delete [] pstate_seq;
}


/*--------------------------------------------------------------*
 *
 *  Calcul d'un melange de lois du nombre d'occurrences d'une observation
 *  emise par une chaine de Markov cachee pour une distribution
 *  des longueurs de sequences donnee.
 *
 *  arguments : indice du processus d'observation, observation.
 *
 *--------------------------------------------------------------*/

void Markov::output_nb_occurrence_mixture(int variable , int output)

{
  register int i , j , k , m;
  int max_length , nb_pattern , index_nb_pattern , power[ORDER] , state_index[ORDER];
  double sum , *observation , *state_seq , *pstate_seq , *states , *pstates ,
         *pmass , *lmass , **ptransition;
  Distribution *nb_occurrence;


  nb_occurrence = process[variable]->nb_occurrence[output];

  max_length = process[variable]->length->nb_value - 1;
  nb_pattern = max_length + 1;

  pmass = nb_occurrence->mass;
  for (i = 0;i < nb_occurrence->nb_value;i++) {
    *pmass++ = 0.;
  }

  observation = new double[nb_state];
  for (i = 0;i < nb_state;i++) {
    observation[i] = process[variable]->observation[i]->mass[output];
  }

  i = 1;
  for (j = 0;j < order;j++) {
    power[j] = i;
    i *= nb_state;
  }

  state_seq = new double[nb_row * nb_pattern];
  pstate_seq = new double[nb_row * nb_pattern];

  lmass = process[variable]->length->mass;
  index_nb_pattern = 1;

  for (i = 0;i < max_length;i++) {

    // initialisation des probabilites des sequences d'etats
    // pour un nombre d'occurrences donne de l'observation selectionnee

    if (i == 0) {
      states = state_seq;
      j = 0;

      for (k = 0;k < nb_row;k++) {
        if (k == self_row[j]) {
          *states++ = (1. - observation[j]) * initial[j];
          *states = observation[j] * initial[j];
          j++;
        }
        else {
          *states++ = 0.;
          *states = 0.;
        }
        states += nb_pattern - 1;
      }
    }

    else {

      // mise a jour des probabilites des sequences d'etats

      pstates = pstate_seq;
      states = state_seq;

      for (j = 0;j < nb_row;j++) {
        for (k = 0;k < index_nb_pattern;k++) {
          *pstates++ = *states;
          *states++ = 0.;
        }
        *states = 0.;
        pstates += nb_pattern - index_nb_pattern;
        states += nb_pattern - index_nb_pattern;
      }

      states = state_seq;
      for (j = 0;j < order;j++) {
        state_index[j] = 0;
      }

      for (j = 0;j < nb_row;j++) {

        // calcul des probabilites des sequences d'etats de
        // longueur = ordre de la chaine de Markov pour chaque nombre
        // d'occurrences de l'observation selectionnee

        ptransition = transition;
        pstates = pstate_seq;
        for (k = 0;k < order - 1;k++) {
          ptransition += state_index[k] * power[k + 1];
          pstates += state_index[k] * power[k + 1] * nb_pattern;
        }

        for (k = 0;k < nb_state;k++) {
          for (m = 0;m < index_nb_pattern;m++) {
            *states++ += (1. - observation[state_index[order - 1]]) *
                         *(*ptransition + state_index[order - 1]) * *pstates;
            *states += observation[state_index[order - 1]] *
                       *(*ptransition + state_index[order - 1]) * *pstates++;
          }
          states -= index_nb_pattern;
          ptransition++;
          pstates += nb_pattern - index_nb_pattern;
        }
        states += nb_pattern;

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

    index_nb_pattern++;

    // mise a jour du melange de lois du nombre d'occurrences
    // de l'observation selectionnee

    if (*++lmass > 0.) {
      pmass = nb_occurrence->mass;
      for (j = 0;j < index_nb_pattern;j++) {
        states = state_seq + j;
        sum = 0.;
        for (k = 0;k < nb_row;k++) {
          sum += *states;
          states += nb_pattern;
        }
        *pmass++ += *lmass * sum;
      }
    }
  }

  nb_occurrence->nb_value_computation();
  nb_occurrence->offset_computation();
  nb_occurrence->cumul_computation();

  nb_occurrence->max_computation();
  nb_occurrence->mean_computation();
  nb_occurrence->variance_computation();

  delete [] observation;
  delete [] state_seq;
  delete [] pstate_seq;
}


/*--------------------------------------------------------------*
 *
 *  Calcul des lois caracteristiques d'un objet Markov.
 *
 *  arguments : longueur des sequences, flag sur le calcul des lois de comptage,
 *              indice du processus d'observation.
 *
 *--------------------------------------------------------------*/

void Markov::characteristic_computation(int length , bool counting_flag , int variable)

{
  if (nb_component > 0) {
    register int i , j , k;
    double *memory;
    Parametric dlength(UNIFORM , length , length , D_DEFAULT , D_DEFAULT);


    memory = 0;

    // calcul des lois caracteristiques au niveau etat

    if (((variable == I_DEFAULT) || (variable == 0)) &&
        ((!(process[0]->length)) || (dlength != *(process[0]->length)))) {
      process[0]->create_characteristic(dlength , true , counting_flag);

      memory = memory_computation();

      index_state_distribution();

      for (i = 0;i < nb_state;i++) {
        state_no_occurrence_probability(i);
        state_first_occurrence_distribution(i);

        state_leave_probability(memory , i);
        state_recurrence_time_distribution(memory , i);

        if (state_type[i] != 'a') {
          state_sojourn_time_distribution(memory , i);
        }
        else {
          process[0]->absorption[i] = 1.;
          delete process[0]->sojourn_time[i];
          process[0]->sojourn_time[i] = 0;
        }

        if (counting_flag) {
          state_nb_pattern_mixture(i , 'r');
          state_nb_pattern_mixture(i , 'o');
        }
      }
    }

    // calcul des lois caracteristiques au niveau observation

    for (i = 1;i <= nb_output_process;i++) {
      if (((variable == I_DEFAULT) || (i == variable)) &&
          ((!(process[i]->length)) || (dlength != *(process[i]->length)))) {
        process[i]->create_characteristic(dlength , true , counting_flag);

        if (!memory) {
          memory = memory_computation();
        }

        index_output_distribution(i);

        for (j = 0;j < process[i]->nb_value;j++) {
          output_no_occurrence_probability(i , j);
          output_first_occurrence_distribution(i , j);

          output_leave_probability(memory , i , j);
          output_recurrence_time_distribution(memory , i , j);

          for (k = 0;k < nb_state;k++) {
            if ((process[i]->observation[k]->mass[j] > 0.) &&
                ((state_type[k] != 'a') || (process[i]->observation[k]->mass[j] < 1.))) {
              break;
            }
          }

          if (k < nb_state) {
            output_sojourn_time_distribution(memory , i , j);
          }
          else {
            process[i]->absorption[j] = 1.;
            delete process[i]->sojourn_time[j];
            process[i]->sojourn_time[j] = 0;
          }

          if (counting_flag) {
            output_nb_run_mixture(i , j);
            output_nb_occurrence_mixture(i , j);
          }
        }
      }
    }

    delete [] memory;
  }
}


/*--------------------------------------------------------------*
 *
 *  Calcul des lois caracteristiques d'un objet Markov.
 *
 *  arguments : reference sur un objet Markov_data,
 *              flag sur le calcul des lois de comptage,
 *              indice du processus d'observation,
 *              flag pour tenir compte des longueurs.
 *
 *--------------------------------------------------------------*/

void Markov::characteristic_computation(const Markov_data &seq , bool counting_flag ,
                                        int variable , bool length_flag)

{
  if (nb_component > 0) {
    register int i , j , k;
    int seq_variable;
    double *memory;
    Distribution dlength(*(seq.hlength));


    memory = 0;

    // calcul des lois caracteristiques au niveau etat

    if (((variable == I_DEFAULT) || (variable == 0)) && ((!length_flag) ||
         ((length_flag) && ((!(process[0]->length)) || (dlength != *(process[0]->length)))))) {
      process[0]->create_characteristic(dlength , true , counting_flag);

      memory = memory_computation();

      index_state_distribution();

      for (i = 0;i < nb_state;i++) {
        state_no_occurrence_probability(i);
        if (seq.type[0] == STATE) {
          state_first_occurrence_distribution(i , (seq.characteristics[0] ? seq.characteristics[0]->first_occurrence[i]->nb_value : 1));
        }
        else {
          state_first_occurrence_distribution(i);
        }

        state_leave_probability(memory , i);
        if (seq.type[0] == STATE) {
          state_recurrence_time_distribution(memory , i , (seq.characteristics[0] ? seq.characteristics[0]->recurrence_time[i]->nb_value : 1));
        }
        else {
          state_recurrence_time_distribution(memory , i);
        }

        if (state_type[i] != 'a') {
          if (seq.type[0] == STATE) {
            state_sojourn_time_distribution(memory , i , (seq.characteristics[0] ? seq.characteristics[0]->sojourn_time[i]->nb_value : 1));
          }
          else {
            state_sojourn_time_distribution(memory , i);
          }
        }

        else {
          process[0]->absorption[i] = 1.;
          delete process[0]->sojourn_time[i];
          process[0]->sojourn_time[i] = 0;
        }

        if (counting_flag) {
          state_nb_pattern_mixture(i , 'r');
          state_nb_pattern_mixture(i , 'o');
        }
      }
    }

    // calcul des lois caracteristiques au niveau observation

    for (i = 1;i <= nb_output_process;i++) {
      if (((variable == I_DEFAULT) || (i == variable)) && ((!length_flag) ||
           ((length_flag) && ((!(process[i]->length)) ||
             (dlength != *(process[i]->length)))))) {
        switch (seq.type[0]) {
        case INT_VALUE :
          seq_variable = i - 1;
          break;
        case STATE :
          seq_variable = i;
          break;
        }

        process[i]->create_characteristic(dlength , true , counting_flag);

        if (!memory) {
          memory = memory_computation();
        }

        index_output_distribution(i);

        for (j = 0;j < process[i]->nb_value;j++) {
          output_no_occurrence_probability(i , j);
          output_first_occurrence_distribution(i , j , (seq.characteristics[seq_variable] ? seq.characteristics[seq_variable]->first_occurrence[j]->nb_value : 1));

          output_leave_probability(memory , i , j);
          output_recurrence_time_distribution(memory , i , j , (seq.characteristics[seq_variable] ? seq.characteristics[seq_variable]->recurrence_time[j]->nb_value : 1));

          for (k = 0;k < nb_state;k++) {
            if ((process[i]->observation[k]->mass[j] > 0.) &&
                ((state_type[k] != 'a') || (process[i]->observation[k]->mass[j] < 1.))) {
              break;
            }
          }

          if (k < nb_state) {
            output_sojourn_time_distribution(memory , i , j , (seq.characteristics[seq_variable] ? seq.characteristics[seq_variable]->sojourn_time[j]->nb_value : 1));
          }
          else {
            process[i]->absorption[j] = 1.;
            delete process[i]->sojourn_time[j];
            process[i]->sojourn_time[j] = 0;
          }

          if (counting_flag) {
            output_nb_run_mixture(i , j);
            output_nb_occurrence_mixture(i , j);
          }
        }
      }
    }

    delete [] memory;
  }
}
