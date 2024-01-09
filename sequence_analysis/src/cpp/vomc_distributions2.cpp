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
 *       $Id: vomc_distributions2.cpp 18084 2015-04-23 11:00:58Z guedon $
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

#include "sequences.h"
#include "variable_order_markov.h"
#include "sequence_label.h"

using namespace std;
using namespace stat_tool;


namespace sequence_analysis {



/*--------------------------------------------------------------*
 *
 *  Calcul d'un melange de lois du nombre de series d'un etat ('r') ou
 *  du nombre d'occurrences d'un etat ('o') d'une chaine de Markov
 *  d'ordre variable pour une distribution des longueurs de sequences donnee.
 *
 *  arguments : etat, type de forme.
 *
 *--------------------------------------------------------------*/

void VariableOrderMarkovChain::state_nb_pattern_mixture(int istate , char pattern)

{
  register int i , j , k , m;
  int max_length , nb_pattern , index_nb_pattern , increment;
  double sum , *pmass , *lmass , **memory , **previous_memory;
  Distribution *pdist;


  max_length = state_process->length->nb_value - 1;

  switch (pattern) {
  case 'r' :
    pdist = state_process->nb_run[istate];
    nb_pattern = max_length / 2 + 2;
    break;
  case 'o' :
    pdist = state_process->nb_occurrence[istate];
    nb_pattern = max_length + 1;
    break;
  }

  pmass = pdist->mass;
  for (i = 0;i < pdist->nb_value;i++) {
    *pmass++ = 0.;
  }

  memory = new double*[nb_row];
  previous_memory = new double*[nb_row];
  for (i = 1;i < nb_row;i++) {
    memory[i] = new double[nb_pattern];
    previous_memory[i] = new double[nb_pattern];
  }

  lmass = state_process->length->mass;
  index_nb_pattern = 1;

  for (i = 0;i < max_length;i++) {

    // initialisation des probabilites des memoires
    // pour un nombre de formes donne de l'etat selectionne

    if (i == 0) {
      switch (type) {

      case 'o' : {
        for (j = 1;j < nb_row;j++) {
          if (order[j] == 1) {
            if (state[j][0] == istate) {
              memory[j][0] = 0.;
              memory[j][1] = initial[state[j][0]];
            }
            else {
              memory[j][0] = initial[state[j][0]];
              memory[j][1] = 0.;
            }
          }

          else {
            memory[j][0] = 0.;
            memory[j][1] = 0.;
          }
        }
        break;
      }

      case 'e' : {
        for (j = 1;j < nb_row;j++) {
          if (!child[j]) {
            if (state[j][0] == istate) {
              memory[j][0] = 0.;
              memory[j][1] = initial[j];
            }
            else {
              memory[j][0] = initial[j];
              memory[j][1] = 0.;
            }
          }

          else {
            memory[j][0] = 0.;
            memory[j][1] = 0.;
          }
        }
        break;
      }
      }
    }

    else {

      // mise a jour des probabilites des memoires

      for (j = 1;j < nb_row;j++) {
        for (k = 0;k < index_nb_pattern;k++) {
          previous_memory[j][k] = memory[j][k];
          memory[j][k] = 0.;
        }
        memory[j][index_nb_pattern] = 0.;
      }

      for (j = 1;j < nb_row;j++) {

        // calcul des probabilites des memoires pour chaque nombre de formes
        // de l'etat selectionne

        for (k = 0;k < nb_memory[j];k++) {
          switch (pattern) {
          case 'r' :
            increment = (((state[j][0] == istate) &&
                          (((order[j] == 1) && (state[previous[j][k]][0] != istate)) ||
                           ((order[j] > 1) && (state[j][1] != istate)))) ? 1 : 0);
            break;
          case 'o' :
            increment = (state[j][0] == istate ? 1 : 0);
            break;
          }

          for (m = 0;m < index_nb_pattern;m++) {
             memory[j][m + increment] += transition[previous[j][k]][state[j][0]] *
                                         previous_memory[previous[j][k]][m];
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
        sum = 0.;
        for (k = 1;k < nb_row;k++) {
          sum += memory[k][j];
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

  for (i = 1;i < nb_row;i++) {
    delete [] memory[i];
    delete [] previous_memory[i];
  }
  delete [] memory;
  delete [] previous_memory;
}


/*--------------------------------------------------------------*
 *
 *  Calcul d'un melange de lois du nombre de series d'une observation
 *  emise par une chaine de Markov d'ordre variable cachee
 *  pour une distribution des longueurs de sequences donnee.
 *
 *  arguments : indice du processus d'observation, observation.
 *
 *--------------------------------------------------------------*/

void VariableOrderMarkov::output_nb_run_mixture(int variable , int output)

{
  register int i , j , k , m;
  int max_length , nb_pattern , index_nb_pattern;
  double sum , *observation , *pmass , *lmass , **memory , **previous_memory;
  Distribution *nb_run;


  nb_run = categorical_process[variable]->nb_run[output];

  max_length = categorical_process[variable]->length->nb_value - 1;
  nb_pattern = max_length / 2 + 2;

  pmass = nb_run->mass;
  for (i = 0;i < nb_run->nb_value;i++) {
    *pmass++ = 0.;
  }

  observation = new double[nb_state];
  for (i = 0;i < nb_state;i++) {
    observation[i] = categorical_process[variable]->observation[i]->mass[output];
  }

  memory = new double*[nb_row];
  previous_memory = new double*[nb_row];
  for (i = 1;i < nb_row;i++) {
    memory[i] = new double[nb_pattern * 2];
    previous_memory[i] = new double[nb_pattern * 2];
  }

  lmass = categorical_process[variable]->length->mass;
  index_nb_pattern = 1;

  for (i = 0;i < max_length;i++) {

    // initialisation des probabilites des memoires
    // pour un nombre de series donne de l'observation selectionnee

    if (i == 0) {
      switch (type) {

      case 'o' : {
        for (j = 1;j < nb_row;j++) {
          if (order[j] == 1) {
            memory[j][0] = (1. - observation[state[j][0]]) * initial[state[j][0]];
            memory[j][1] = 0.;
            memory[j][2] = 0.;
            memory[j][3] = observation[state[j][0]] * initial[state[j][0]];
          }
          else {
            for (k = 0;k < 4;k++) {
              memory[j][k] = 0.;
            }
          }
        }
        break;
      }

      case 'e' : {
        for (j = 1;j < nb_row;j++) {
          if (!child[j]) {
            memory[j][0] = (1. - observation[state[j][0]]) * initial[j];
            memory[j][1] = 0.;
            memory[j][2] = 0.;
            memory[j][3] = observation[state[j][0]] * initial[j];
          }
          else {
            for (k = 0;k < 4;k++) {
              memory[j][k] = 0.;
            }
          }
        }
        break;
      }
      }
    }

    else {

      // mise a jour des probabilites des memoires

      for (j = 1;j < nb_row;j++) {
        for (k = 0;k < index_nb_pattern * 2;k++) {
          previous_memory[j][k] = memory[j][k];
          memory[j][k] = 0.;
        }
        for (k = index_nb_pattern * 2;k < (index_nb_pattern + 1) * 2;k++) {
          previous_memory[j][k] = 0.;
          memory[j][k] = 0.;
        }
      }

      for (j = 1;j < nb_row;j++) {

        // calcul des probabilites des memoires
        // pour chaque nombre de series de l'observation selectionnee

        for (k = 0;k < nb_memory[j];k++) {
          for (m = 0;m <= index_nb_pattern;m++) {
            if (m < index_nb_pattern) {
              memory[j][m * 2] += (1. - observation[state[j][0]]) * transition[previous[j][k]][state[j][0]] *
                                  (previous_memory[previous[j][k]][m * 2] + previous_memory[previous[j][k]][m * 2 + 1]);
            }
            if (m > 0) {
              memory[j][m * 2 + 1] += observation[state[j][0]] * transition[previous[j][k]][state[j][0]] *
                                      (previous_memory[previous[j][k]][(m - 1) * 2] + previous_memory[previous[j][k]][m * 2 + 1]);
            }
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
        sum = 0.;
        for (k = 1;k < nb_row;k++) {
          sum += memory[k][j * 2] + memory[k][j * 2 + 1];
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

  for (i = 1;i < nb_row;i++) {
    delete [] memory[i];
    delete [] previous_memory[i];
  }
  delete [] memory;
  delete [] previous_memory;
}


/*--------------------------------------------------------------*
 *
 *  Calcul d'un melange de lois du nombre d'occurrences d'une observation
 *  emise par une chaine de Markov d'ordre variable cachee
 *  pour une distribution des longueurs de sequences donnee.
 *
 *  arguments : indice du processus d'observation, observation.
 *
 *--------------------------------------------------------------*/

void VariableOrderMarkov::output_nb_occurrence_mixture(int variable , int output)

{
  register int i , j , k , m;
  int max_length , nb_pattern , index_nb_pattern;
  double sum , *observation , *pmass , *lmass , **memory , **previous_memory;
  Distribution *nb_occurrence;


  nb_occurrence = categorical_process[variable]->nb_occurrence[output];

  max_length = categorical_process[variable]->length->nb_value - 1;
  nb_pattern = max_length + 1;

  pmass = nb_occurrence->mass;
  for (i = 0;i < nb_occurrence->nb_value;i++) {
    *pmass++ = 0.;
  }

  observation = new double[nb_state];
  for (i = 0;i < nb_state;i++) {
    observation[i] = categorical_process[variable]->observation[i]->mass[output];
  }

  memory = new double*[nb_row];
  previous_memory = new double*[nb_row];
  for (i = 1;i < nb_row;i++) {
    memory[i] = new double[nb_pattern];
    previous_memory[i] = new double[nb_pattern];
  }

  lmass = categorical_process[variable]->length->mass;
  index_nb_pattern = 1;

  for (i = 0;i < max_length;i++) {

    // initialisation des probabilites des memoires
    // pour un nombre d'occurrences donne de l'observation selectionnee

    if (i == 0) {
      switch (type) {

      case 'o' : {
        for (j = 1;j < nb_row;j++) {
          if (order[j] == 1) {
            memory[j][0] = (1. - observation[state[j][0]]) * initial[state[j][0]];
            memory[j][1] = observation[state[j][0]] * initial[state[j][0]];
          }
          else {
            memory[j][0] = 0.;
            memory[j][1] = 0.;
          }
        }
        break;
      }

      case 'e' : {
        for (j = 1;j < nb_row;j++) {
          if (!child[j]) {
            memory[j][0] = (1. - observation[state[j][0]]) * initial[j];
            memory[j][1] = observation[state[j][0]] * initial[j];
          }
          else {
            memory[j][0] = 0.;
            memory[j][1] = 0.;
          }
        }
        break;
      }
      }
    }

    else {

      // mise a jour des probabilites des memoires

      for (j = 1;j < nb_row;j++) {
        for (k = 0;k < index_nb_pattern;k++) {
          previous_memory[j][k] = memory[j][k];
          memory[j][k] = 0.;
        }
        memory[j][index_nb_pattern] = 0.;
      }

      for (j = 0;j < nb_row;j++) {

        // calcul des probabilites des memoires
        // pour chaque nombre d'occurrences de l'observation selectionnee

        for (k = 0;k < nb_memory[j];k++) {
          for (m = 0;m < index_nb_pattern;m++) {
            memory[j][m] += (1. - observation[state[j][0]]) * transition[previous[j][k]][state[j][0]] *
                            previous_memory[previous[j][k]][m];
            memory[j][m + 1] += observation[state[j][0]] * transition[previous[j][k]][state[j][0]] *
                                previous_memory[previous[j][k]][m];
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
        sum = 0.;
        for (k = 1;k < nb_row;k++) {
          sum += memory[k][j];
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

  for (i = 1;i < nb_row;i++) {
    delete [] memory[i];
    delete [] previous_memory[i];
  }
  delete [] memory;
  delete [] previous_memory;
}


/*--------------------------------------------------------------*
 *
 *  Calcul des lois caracteristiques d'un objet VariableOrderMarkov.
 *
 *  arguments : longueur des sequences, flag sur le calcul des lois de comptage,
 *              indice du processus d'observation.
 *
 *--------------------------------------------------------------*/

void VariableOrderMarkov::characteristic_computation(int length , bool counting_flag ,
                                                     int variable)

{
  if (nb_component > 0) {
    bool computation[NB_OUTPUT_PROCESS + 1];
    register int i , j , k;
    double *memory;
    DiscreteParametric dlength(UNIFORM , length , length , D_DEFAULT , D_DEFAULT);


    memory = NULL;

    // calcul des lois de type intensite et intervalle au niveau etat

    if (((variable == I_DEFAULT) || (variable == 0)) &&
        ((!(state_process->length)) ||
         (dlength != *(state_process->length)))) {
      computation[0] = true;
      state_process->create_characteristic(dlength , true , counting_flag);

      switch (type) {

      case 'o' : {
        memory = memory_computation();
        break;
      }

      case 'e' : {
        memory = new double[nb_row];
        for (i = 1;i < nb_row;i++) {
          memory[i] = initial[i];
        }
        break;
      }
      }

      index_state_distribution();

      for (i = 0;i < nb_state;i++) {
        if (type == 'o') {
          state_no_occurrence_probability(i);
        }
        state_first_occurrence_distribution(i);

        if (type == 'o') {
          state_leave_probability(memory , i);
        }
        state_recurrence_time_distribution(memory , i);

        if (state_type[i] != 'a') {
          state_sojourn_time_distribution(memory , i);
        }
        else {
          state_process->absorption[i] = 1.;
          delete state_process->sojourn_time[i];
          state_process->sojourn_time[i] = NULL;
        }
      }

#     ifdef DEBUG
      if (type == 'e') {
        double sum = 0.;

        // calcul de la loi stationnaire dans le cas d'un processus en equilibre
        // avec renormalisation pour tenir compte des seuils appliques sur
        // les fonctions de repartition des lois de temps de retour dans les etats

        cout << "\n" << STAT_label[STATL_STATIONARY_PROBABILITIES] << endl;
        for (i = 1;i < nb_row;i++) {
          cout << initial[i] << " ";
        }
        cout << endl;

        for (i = 0;i < nb_state;i++) {
          sum += 1. / state_process->recurrence_time[i]->mean;
        }
        for (i = 0;i < nb_state;i++) {
          cout << 1. / (state_process->recurrence_time[i]->mean * sum) << " ";
        }
        cout << endl;
      }
#     endif

    }

    else {
      computation[0] = false;
    }

    // calcul des lois de type intensite et intervalle au niveau observation

    for (i = 0;i < nb_output_process;i++) {
      if ((categorical_process[i]) && ((variable == I_DEFAULT) || (i == variable)) &&
          ((!(categorical_process[i]->length)) ||
           (dlength != *(categorical_process[i]->length)))) {
        computation[i + 1] = true;
        categorical_process[i]->create_characteristic(dlength , true , counting_flag);

        index_output_distribution(i);

        if (!memory) {
          switch (type) {

          case 'o' : {
            memory = memory_computation();
            break;
          }

          case 'e' : {
            memory = new double[nb_row];
            for (j = 1;j < nb_row;j++) {
              memory[j] = initial[j];
            }
            break;
          }
          }
        }

        for (j = 0;j < categorical_process[i]->nb_value;j++) {
          if (type == 'o') {
            output_no_occurrence_probability(i , j);
          }
          if (categorical_process[i]->no_occurrence[j] < 1. - DOUBLE_ERROR) {
            output_first_occurrence_distribution(i , j);
          }
          else {
            delete categorical_process[i]->first_occurrence[j];
            categorical_process[i]->first_occurrence[j] = NULL;
            categorical_process[i]->leave[j] = 1.;
          }

          if ((type == 'o') && (categorical_process[i]->first_occurrence[j])) {
            output_leave_probability(memory , i , j);
          }
          if (categorical_process[i]->leave[j] < 1. - DOUBLE_ERROR) {
            output_recurrence_time_distribution(memory , i , j);
          }
          else {
            delete categorical_process[i]->recurrence_time[j];
            categorical_process[i]->recurrence_time[j] = NULL;
          }

          for (k = 0;k < nb_state;k++) {
            if ((categorical_process[i]->observation[k]->mass[j] > 0.) &&
                ((state_type[k] != 'a') || (categorical_process[i]->observation[k]->mass[j] < 1.))) {
              break;
            }
          }

          if (k < nb_state) {
            output_sojourn_time_distribution(memory , i , j);
          }
          else {
            categorical_process[i]->absorption[j] = 1.;
            delete categorical_process[i]->sojourn_time[j];
            categorical_process[i]->sojourn_time[j] = NULL;
          }
        }
      }

      else {
        computation[i + 1] = false;
      }
    }

    delete [] memory;

    if (counting_flag) {

      // calcul des lois de comptage au niveau etat

      if (computation[0]) {
        for (i = 0;i < nb_state;i++) {
          state_nb_pattern_mixture(i , 'r');
          state_nb_pattern_mixture(i , 'o');
        }
      }

      // calcul des lois de comptage au niveau observation

      for (i = 0;i < nb_output_process;i++) {
        if (computation[i + 1]) {
          for (j = 0;j < categorical_process[i]->nb_value;j++) {
            output_nb_run_mixture(i , j);
            output_nb_occurrence_mixture(i , j);
          }
        }
      }
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Calcul des lois caracteristiques d'un objet VariableOrderMarkov.
 *
 *  arguments : reference sur un objet VariableOrderMarkovData,
 *              flag sur le calcul des lois de comptage,
 *              indice du processus d'observation,
 *              flag pour tenir compte des longueurs.
 *
 *--------------------------------------------------------------*/

void VariableOrderMarkov::characteristic_computation(const VariableOrderMarkovData &seq ,
                                                     bool counting_flag , int variable ,
                                                     bool length_flag)

{
  if (nb_component > 0) {
    bool computation[NB_OUTPUT_PROCESS + 1];
    register int i , j , k;
    int seq_variable;
    double *memory;
    Distribution dlength(*(seq.length_distribution));


    memory = NULL;

    // calcul des lois caracteristiques au niveau etat

    if  (((variable == I_DEFAULT) || (variable == 0)) && ((!length_flag) ||
          ((length_flag) && ((!(state_process->length)) ||
            (dlength != *(state_process->length)))))) {
      computation[0] = true;
      state_process->create_characteristic(dlength , true , counting_flag);

      switch (type) {

      case 'o' : {
        memory = memory_computation();
        break;
      }

      case 'e' : {
        memory = new double[nb_row];
        for (i = 1;i < nb_row;i++) {
          memory[i] = initial[i];
        }
        break;
      }
      }

      index_state_distribution();

      for (i = 0;i < nb_state;i++) {
        if (type == 'o') {
          state_no_occurrence_probability(i);
        }
        if (seq.type[0] == STATE) {
          state_first_occurrence_distribution(i , ((seq.characteristics[0]) && (i < seq.marginal_distribution[0]->nb_value) && (seq.characteristics[0]->first_occurrence[i]->nb_element > 0) ? seq.characteristics[0]->first_occurrence[i]->nb_value : 1));
        }
        else {
          state_first_occurrence_distribution(i);
        }

        if (type == 'o') {
          state_leave_probability(memory , i);
        }
        if (seq.type[0] == STATE) {
          state_recurrence_time_distribution(memory , i , ((seq.characteristics[0]) && (i < seq.marginal_distribution[0]->nb_value) && (seq.characteristics[0]->recurrence_time[i]->nb_element > 0) ? seq.characteristics[0]->recurrence_time[i]->nb_value : 1));
        }
        else {
          state_recurrence_time_distribution(memory , i);
        }

        if (state_type[i] != 'a') {
          if (seq.type[0] == STATE) {
            state_sojourn_time_distribution(memory , i , ((seq.characteristics[0]) && (i < seq.marginal_distribution[0]->nb_value) && (seq.characteristics[0]->sojourn_time[i]->nb_element > 0) ? seq.characteristics[0]->sojourn_time[i]->nb_value : 1));
          }
          else {
            state_sojourn_time_distribution(memory , i);
          }
        }

        else {
          state_process->absorption[i] = 1.;
          delete state_process->sojourn_time[i];
          state_process->sojourn_time[i] = NULL;
        }
      }

#     ifdef DEBUG
      if (type == 'e') {
        double sum = 0.;

        // calcul de la loi stationnaire dans le cas d'un processus en equilibre
        // avec renormalisation pour tenir compte des seuils appliques sur
        // les fonctions de repartition des lois de temps de retour dans les etats

        cout << "\n" << STAT_label[STATL_STATIONARY_PROBABILITIES] << endl;
        for (i = 1;i < nb_row;i++) {
          cout << initial[i] << " ";
        }
        cout << endl;

        for (i = 0;i < nb_state;i++) {
          sum += 1. / state_process->recurrence_time[i]->mean;
        }
        for (i = 0;i < nb_state;i++) {
          cout << 1. / (state_process->recurrence_time[i]->mean * sum) << " ";
        }
        cout << endl;
      }
#     endif

    }

    else {
      computation[0] = false;
    }

    // calcul des lois de type intensite et intervalle au niveau observation

    for (i = 0;i < nb_output_process;i++) {
      if ((categorical_process[i]) && ((variable == I_DEFAULT) || (variable == 1)) &&
          ((!length_flag) || ((length_flag) && ((!(categorical_process[i]->length)) ||
           (dlength != *(categorical_process[i]->length)))))) {
        computation[i + 1] = true;
        categorical_process[i]->create_characteristic(dlength , true , counting_flag);

        switch (seq.type[0]) {
        case STATE :
          seq_variable = i + 1;
          break;
        default :
          seq_variable = i;
          break;
        }

        index_output_distribution(i);

        if (!memory) {
          switch (type) {

          case 'o' : {
            memory = memory_computation();
            break;
          }

          case 'e' : {
            memory = new double[nb_row];
            for (j = 1;j < nb_row;j++) {
              memory[j] = initial[j];
            }
            break;
          }
          }
        }

        for (j = 0;j < categorical_process[i]->nb_value;j++) {
          if (type == 'o') {
            output_no_occurrence_probability(i , j);
          }
          if (categorical_process[i]->no_occurrence[j] < 1. - DOUBLE_ERROR) {
            output_first_occurrence_distribution(i , j , ((seq.characteristics[seq_variable]) && (j < seq.characteristics[seq_variable]->nb_value) && (seq.characteristics[seq_variable]->first_occurrence[j]->nb_element > 0) ? seq.characteristics[seq_variable]->first_occurrence[j]->nb_value : 1));
          }
          else {
            delete categorical_process[i]->first_occurrence[j];
            categorical_process[i]->first_occurrence[j] = NULL;
            categorical_process[i]->leave[j] = 1.;
          }

          if ((type == 'o') && (categorical_process[i]->first_occurrence[j])) {
            output_leave_probability(memory , i , j);
          }
          if (categorical_process[i]->leave[j] < 1. - DOUBLE_ERROR) {
            output_recurrence_time_distribution(memory , i , j , ((seq.characteristics[seq_variable]) && (j < seq.characteristics[seq_variable]->nb_value) && (seq.characteristics[seq_variable]->recurrence_time[j]->nb_element > 0) ? seq.characteristics[seq_variable]->recurrence_time[j]->nb_value : 1));
          }
          else {
            delete categorical_process[i]->recurrence_time[j];
            categorical_process[i]->recurrence_time[j] = NULL;
          }

          for (k = 0;k < nb_state;k++) {
            if ((categorical_process[i]->observation[k]->mass[j] > 0.) &&
                ((state_type[k] != 'a') || (categorical_process[i]->observation[k]->mass[j] < 1.))) {
              break;
            }
          }

          if (k < nb_state) {
            output_sojourn_time_distribution(memory , i , j , ((seq.characteristics[seq_variable]) && (j < seq.characteristics[seq_variable]->nb_value) && (seq.characteristics[seq_variable]->sojourn_time[j]->nb_element > 0) ? seq.characteristics[seq_variable]->sojourn_time[j]->nb_value : 1));
          }
          else {
            categorical_process[i]->absorption[j] = 1.;
            delete categorical_process[i]->sojourn_time[j];
            categorical_process[i]->sojourn_time[j] = NULL;
          }
        }
      }

      else {
        computation[i + 1] = false;
      }
    }

    delete [] memory;

    if (counting_flag) {

      // calcul des lois de comptage au niveau etat

      if (computation[0]) {
        for (i = 0;i < nb_state;i++) {
          state_nb_pattern_mixture(i , 'r');
          state_nb_pattern_mixture(i , 'o');
        }
      }

      // calcul des lois de comptage au niveau observation

      for (i = 0;i < nb_output_process;i++) {
        if (computation[i + 1]) {
          for (j = 0;j < categorical_process[i]->nb_value;j++) {
            output_nb_run_mixture(i , j);
            output_nb_occurrence_mixture(i , j);
          }
        }
      }
    }
  }
}


};  // namespace sequence_analysis
