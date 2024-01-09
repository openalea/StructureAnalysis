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
 *       $Id: hvomc_algorithms1.cpp 18056 2015-04-23 09:47:19Z guedon $
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
#include <sstream>

#include "stat_tool/stat_tools.h"
#include "stat_tool/curves.h"
#include "stat_tool/distribution.h"
#include "stat_tool/markovian.h"
#include "stat_tool/vectors.h"
#include "stat_tool/distance_matrix.h"
#include "stat_tool/stat_label.h"

#include "stat_tool/distribution_reestimation.hpp"   // probleme compilateur C++ Windows

#include "sequences.h"
#include "variable_order_markov.h"
#include "hidden_variable_order_markov.h"
#include "sequence_label.h"

using namespace std;
using namespace stat_tool;


namespace sequence_analysis {



/*--------------------------------------------------------------*
 *
 *  Calcul de la vraisemblance de sequences pour une chaine de Markov
 *  d'ordre variable cachee par l'algorithme forward.
 *
 *  arguments : reference sur un objet MarkovianSequences, pointeur sur
 *              les probabilites a posteriori des sequences d'etats
 *              les plus probables, indice de la sequence.
 *
 *--------------------------------------------------------------*/

double HiddenVariableOrderMarkov::likelihood_computation(const MarkovianSequences &seq ,
                                                         double *posterior_probability , int index) const

{
  register int i , j , k , m;
  int nb_value , **pioutput;
  double likelihood = 0. , seq_likelihood , *forward , *auxiliary , norm , **proutput;


  // verification de la compatibilite entre le modele et les donnees

  if (nb_output_process == seq.nb_variable) {
    for (i = 0;i < nb_output_process;i++) {
      if ((categorical_process[i]) || (discrete_parametric_process[i])) {
        if (categorical_process[i]) {
          nb_value = categorical_process[i]->nb_value;
        }
        else  {
          nb_value = discrete_parametric_process[i]->nb_value;
        }

        if (nb_value < seq.marginal_distribution[i]->nb_value) {
          likelihood = D_INF;
          break;
        }
      }
    }
  }

  else {
    likelihood = D_INF;
  }

  if (likelihood != D_INF) {

    // initialisations

    forward = new double[nb_row];
    auxiliary = new double[nb_row];

    pioutput = new int*[seq.nb_variable];
    proutput = new double*[seq.nb_variable];

    for (i = 0;i < seq.nb_sequence;i++) {
      if ((index == I_DEFAULT) || (index == i)) {
        for (j = 0;j < seq.nb_variable;j++) {
          switch (seq.type[j]) {
          case INT_VALUE :
            pioutput[j] = seq.int_sequence[i][j];
            break;
          case REAL_VALUE :
            proutput[j] = seq.real_sequence[i][j];
            break;
          }
        }
        seq_likelihood = 0.;

        norm = 0.;

        switch (type) {

        case 'o' : {
          for (j = 1;j < nb_row;j++) {
            if (order[j] == 1) {
              forward[j] = initial[state[j][0]];

              for (k = 0;k < nb_output_process;k++) {
                if (categorical_process[k]) {
                  forward[j] *= categorical_process[k]->observation[state[j][0]]->mass[*pioutput[k]];
                }

                else if (discrete_parametric_process[k]) {
                  forward[j] *= discrete_parametric_process[k]->observation[state[j][0]]->mass[*pioutput[k]];
                }

                else {
                  if (((continuous_parametric_process[k]->ident == GAMMA) ||
                      (continuous_parametric_process[k]->ident == ZERO_INFLATED_GAMMA)) && (seq.min_value[k] < seq.min_interval[k] / 2)) {
                    switch (seq.type[k]) {
                    case INT_VALUE :
                      forward[j] *= continuous_parametric_process[k]->observation[state[j][0]]->mass_computation(*pioutput[k] , *pioutput[k] + seq.min_interval[k]);
                      break;
                    case REAL_VALUE :
                      forward[j] *= continuous_parametric_process[k]->observation[state[j][0]]->mass_computation(*proutput[k] , *proutput[k] + seq.min_interval[k]);
                      break;
                    }
                  }

                  else {
                    switch (seq.type[k]) {
                    case INT_VALUE :
                      forward[j] *= continuous_parametric_process[k]->observation[state[j][0]]->mass_computation(*pioutput[k] - seq.min_interval[k] / 2 , *pioutput[k] + seq.min_interval[k] / 2);
                      break;
                    case REAL_VALUE :
                      forward[j] *= continuous_parametric_process[k]->observation[state[j][0]]->mass_computation(*proutput[k] - seq.min_interval[k] / 2 , *proutput[k] + seq.min_interval[k] / 2);
                      break;
                    }
                  }
                }
              }

              norm += forward[j];
            }

            else {
              forward[j] = 0.;
            }
          }
          break;
        }

        case 'e' : {
          for (j = 1;j < nb_row;j++) {
            if (!child[j]) {
              forward[j] = initial[j];

              for (k = 0;k < nb_output_process;k++) {
                if (categorical_process[k]) {
                  forward[j] *= categorical_process[k]->observation[state[j][0]]->mass[*pioutput[k]];
                }

                else if (discrete_parametric_process[k]) {
                  forward[j] *= discrete_parametric_process[k]->observation[state[j][0]]->mass[*pioutput[k]];
                }

                else {
                  if (((continuous_parametric_process[k]->ident == GAMMA) ||
                      (continuous_parametric_process[k]->ident == ZERO_INFLATED_GAMMA)) && (seq.min_value[k] < seq.min_interval[k] / 2)) {
                    switch (seq.type[k]) {
                    case INT_VALUE :
                      forward[j] *= continuous_parametric_process[k]->observation[state[j][0]]->mass_computation(*pioutput[k] , *pioutput[k] + seq.min_interval[k]);
                      break;
                    case REAL_VALUE :
                      forward[j] *= continuous_parametric_process[k]->observation[state[j][0]]->mass_computation(*proutput[k] , *proutput[k] + seq.min_interval[k]);
                      break;
                    }
                  }

                  else {
                    switch (seq.type[k]) {
                    case INT_VALUE :
                      forward[j] *= continuous_parametric_process[k]->observation[state[j][0]]->mass_computation(*pioutput[k] - seq.min_interval[k] / 2 , *pioutput[k] + seq.min_interval[k] / 2);
                      break;
                    case REAL_VALUE :
                      forward[j] *= continuous_parametric_process[k]->observation[state[j][0]]->mass_computation(*proutput[k] - seq.min_interval[k] / 2 , *proutput[k] + seq.min_interval[k] / 2);
                      break;
                    }
                  }
                }
              }

              norm += forward[j];
            }

            else {
              forward[j] = 0.;
            }
          }
          break;
        }
        }

        if (norm > 0.) {
          for (j = 1;j < nb_row;j++) {
            forward[j] /= norm;
          }
          seq_likelihood += log(norm);
        }

        else {
          seq_likelihood = D_INF;
          break;
        }

        for (j = 1;j < seq.length[i];j++) {
          for (k = 0;k < seq.nb_variable;k++) {
            switch (seq.type[k]) {
            case INT_VALUE :
              pioutput[k]++;
              break;
            case REAL_VALUE :
              proutput[k]++;
              break;
            }
          }
          norm = 0.;

          for (k = 1;k < nb_row;k++) {
            auxiliary[k] = 0.;
            for (m = 0;m < nb_memory[k];m++) {
              auxiliary[k] += transition[previous[k][m]][state[k][0]] * forward[previous[k][m]];
            }

            for (m = 0;m < nb_output_process;m++) {
              if (categorical_process[m]) {
                auxiliary[k] *= categorical_process[m]->observation[state[k][0]]->mass[*pioutput[m]];
              }

              else if (discrete_parametric_process[m]) {
                auxiliary[k] *= discrete_parametric_process[m]->observation[state[k][0]]->mass[*pioutput[m]];
              }

              else {
                if (((continuous_parametric_process[m]->ident == GAMMA) ||
                    (continuous_parametric_process[m]->ident == ZERO_INFLATED_GAMMA)) && (seq.min_value[m] < seq.min_interval[m] / 2)) {
                  switch (seq.type[m]) {
                  case INT_VALUE :
                    auxiliary[k] *= continuous_parametric_process[m]->observation[state[k][0]]->mass_computation(*pioutput[m] , *pioutput[m] + seq.min_interval[m]);
                    break;
                  case REAL_VALUE :
                    auxiliary[k] *= continuous_parametric_process[m]->observation[state[k][0]]->mass_computation(*proutput[m] , *proutput[m] + seq.min_interval[m]);
                    break;
                  }
                }

                else {
                  switch (seq.type[m]) {
                  case INT_VALUE :
                    auxiliary[k] *= continuous_parametric_process[m]->observation[state[k][0]]->mass_computation(*pioutput[m] - seq.min_interval[m] / 2 , *pioutput[m] + seq.min_interval[m] / 2);
                    break;
                  case REAL_VALUE :
                    auxiliary[k] *= continuous_parametric_process[m]->observation[state[k][0]]->mass_computation(*proutput[m] - seq.min_interval[m] / 2 , *proutput[m] + seq.min_interval[m] / 2);
                    break;
                  }
                }
              }
            }

            norm += auxiliary[k];
          }

          if (norm > 0.) {
            for (k = 1;k < nb_row;k++) {
              forward[k] = auxiliary[k] / norm;
            }
            seq_likelihood += log(norm);
          }

          else {
            seq_likelihood = D_INF;
            break;
          }
        }

        if (seq_likelihood != D_INF) {
          likelihood += seq_likelihood;
          if (posterior_probability) {
            posterior_probability[i] = exp(posterior_probability[i] - seq_likelihood);
          }
        }

        else {
          likelihood = D_INF;
          break;
        }
      }
    }

    delete [] forward;
    delete [] auxiliary;

    delete [] pioutput;
    delete [] proutput;
  }

  return likelihood;
}


/*--------------------------------------------------------------*
 *
 *  Estimation des parametres d'une chaine de Markov d'ordre variable cachee
 *  a partir d'un echantillon de sequences par l'algorithme EM.
 *
 *  arguments : reference sur un objet StatError, stream, chaine de Markov cachee initiale,
 *              type d'estimation des probabilites de transition initiale (cas ordinaire),
 *              flag parametres de dispersion communs (processus d'observation continus),
 *              flags sur le calcul des lois de comptage et sur le calcul des sequences
 *              d'etats optimales, nombre d'iterations.
 *
 *--------------------------------------------------------------*/

HiddenVariableOrderMarkov* MarkovianSequences::hidden_variable_order_markov_estimation(StatError &error , ostream &os ,
                                                                                       const HiddenVariableOrderMarkov &ihmarkov ,
                                                                                       bool global_initial_transition ,
                                                                                       bool common_dispersion ,
                                                                                       bool counting_flag , bool state_sequence ,
                                                                                       int nb_iter) const

{
  bool status;
  register int i , j , k , m;
  int nb_terminal , max_nb_value , iter , **pioutput;
  double likelihood = D_INF , previous_likelihood , observation_likelihood , **forward , norm ,
         **predicted , buff , *backward , *auxiliary , ***state_sequence_count , diff , variance ,
         **mean_direction , global_mean_direction , concentration , **proutput;
  Distribution *weight;
  ChainReestimation<double> *chain_reestim;
  Reestimation<double> ***observation_reestim;
  FrequencyDistribution *hobservation;
  HiddenVariableOrderMarkov *hmarkov;
  VariableOrderMarkovData *seq;


  hmarkov = NULL;
  error.init();

  // test nombre de valeurs observees par variable

  status = false;
  for (i = 0;i < nb_variable;i++) {
    if (max_value[i] > min_value[i]) {
      status = true;
      break;
    }
  }

  if (!status) {
    error.update(STAT_error[STATR_VARIABLE_NB_VALUE]);
  }

  for (i = 0;i < nb_variable;i++) {
    if ((type[i] != INT_VALUE) && (type[i] != REAL_VALUE) && (type[i] != STATE)) {
      status = false;
      ostringstream error_message , correction_message;
      error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                    << STAT_error[STATR_VARIABLE_TYPE];
      correction_message << STAT_variable_word[INT_VALUE] << " or "
                         << STAT_variable_word[REAL_VALUE];
      error.correction_update((error_message.str()).c_str() , (correction_message.str()).c_str());
    }
  }

  if (ihmarkov.nb_output_process != nb_variable) {
    status = false;
    error.update(STAT_error[STATR_NB_OUTPUT_PROCESS]);
  }

  else {
    for (i = 0;i < nb_variable;i++) {
      if ((ihmarkov.categorical_process[i]) || (ihmarkov.discrete_parametric_process[i])) {
        if (type[i] == REAL_VALUE) {
          status = false;
          ostringstream error_message;
          error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                        << STAT_error[STATR_VARIABLE_TYPE];
          error.correction_update((error_message.str()).c_str() , STAT_variable_word[INT_VALUE]);
        }

        else {
          if (min_value[i] < 0) {
            status = false;
            ostringstream error_message;
            error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                          << STAT_error[STATR_POSITIVE_MIN_VALUE];
            error.update((error_message.str()).c_str());
          }

          if (!marginal_distribution[i]) {
            status = false;
            ostringstream error_message;
            error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                          << STAT_error[STATR_MARGINAL_FREQUENCY_DISTRIBUTION];
            error.update((error_message.str()).c_str());
          }

          else {
            if (((ihmarkov.categorical_process[i]) &&
                 (ihmarkov.categorical_process[i]->nb_value != marginal_distribution[i]->nb_value)) ||
                ((ihmarkov.discrete_parametric_process[i]) &&
                 (ihmarkov.discrete_parametric_process[i]->nb_value < marginal_distribution[i]->nb_value))) {
              status = false;
              ostringstream error_message;
              error_message << STAT_label[STATL_OUTPUT_PROCESS] << " " << i + 1 << ": "
                            << STAT_error[STATR_NB_OUTPUT];
              error.update((error_message.str()).c_str());
            }

            else if ((ihmarkov.categorical_process[i]) && (!characteristics[i])) {
              for (j = 0;j < marginal_distribution[i]->nb_value;j++) {
                if (marginal_distribution[i]->frequency[j] == 0) {
                  status = false;
                  ostringstream error_message;
                  error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                                << STAT_error[STATR_MISSING_VALUE] << " " << j;
                  error.update((error_message.str()).c_str());
                }
              }
            }
          }
        }
      }
    }
  }

  if ((nb_iter != I_DEFAULT) && (nb_iter < 1)) {
    status = false;
    error.update(STAT_error[STATR_NB_ITERATION]);
  }

  if (status) {

    // creation de la chaine de Markov cachee

    hmarkov = new HiddenVariableOrderMarkov(ihmarkov , false);

    if (hmarkov->type == 'e') {
      nb_terminal = (hmarkov->nb_row - 1) * (hmarkov->nb_state - 1) / hmarkov->nb_state + 1;

      for (i = 1;i < hmarkov->nb_row;i++) {
        if (!hmarkov->child[i]) {
          hmarkov->initial[i] = 1. / (double)nb_terminal;
        }
        else {
          hmarkov->initial[i] = 0.;
        }
      }
    }

    if (common_dispersion) {
      for (i = 0;i < hmarkov->nb_output_process;i++) {
        if (hmarkov->continuous_parametric_process[i]) {
          hmarkov->continuous_parametric_process[i]->tied_dispersion = true;
        }
      }
    }

#   ifdef DEBUG
    cout << *hmarkov;
#   endif

    // creation des structures de donnees de l'algorithme

    forward = new double*[max_length];
    for (i = 0;i < max_length;i++) {
      forward[i] = new double[hmarkov->nb_row];
    }

    predicted = new double*[max_length];
    for (i = 0;i < max_length;i++) {
      predicted[i] = new double[hmarkov->nb_row];
    }

    backward = new double[hmarkov->nb_row];

    auxiliary = new double[hmarkov->nb_row];

    chain_reestim = new ChainReestimation<double>((hmarkov->type == 'o' ?  'o' : 'v') ,
                                                  hmarkov->nb_state , hmarkov->nb_row);

    observation_reestim = new Reestimation<double>**[hmarkov->nb_output_process];
    for (i = 0;i < hmarkov->nb_output_process;i++) {
      if (marginal_distribution[i]) {
        observation_reestim[i] = new Reestimation<double>*[hmarkov->nb_state];
        for (j = 0;j < hmarkov->nb_state;j++) {
          observation_reestim[i][j] = new Reestimation<double>(marginal_distribution[i]->nb_value);
        }
      }

      else {
        observation_reestim[i] = NULL;
      }
    }

    max_nb_value = 0;
    for (i = 0;i < hmarkov->nb_output_process;i++) {
      if ((hmarkov->discrete_parametric_process[i]) &&
          (max_nb_value < marginal_distribution[i]->nb_value)) {
        max_nb_value = marginal_distribution[i]->nb_value;
      }
    }

    if (max_nb_value > 0) {
      hobservation = new FrequencyDistribution(max_nb_value);
    }
    else {
      hobservation = NULL;
    }

    for (i = 0;i < hmarkov->nb_output_process;i++) {
      if (!marginal_distribution[i]) {
        break;
      }
    }

    if (i < hmarkov->nb_output_process) {
      state_sequence_count = new double**[nb_sequence];
      for (i = 0;i < nb_sequence;i++) {
        state_sequence_count[i] = new double*[length[i]];
        for (j = 0;j < length[i];j++) {
          state_sequence_count[i][j] = new double[hmarkov->nb_state];
        }
      }
    }
    else {
      state_sequence_count = NULL;
    }

    for (i = 0;i < hmarkov->nb_output_process;i++) {
      if ((hmarkov->continuous_parametric_process[i]) &&
          (hmarkov->continuous_parametric_process[i]->ident == VON_MISES)) {
        break;
      }
    }

    if (i < hmarkov->nb_output_process) {
      mean_direction = new double*[hmarkov->nb_state];
      for (i = 0;i < hmarkov->nb_state;i++) {
        mean_direction[i] = new double[4];
      }
    }
    else {
      mean_direction = NULL;
    }

    pioutput = new int*[nb_variable];
    proutput = new double*[nb_variable];

    iter = 0;
    do {
      iter++;
      previous_likelihood = likelihood;
      likelihood = 0.;

      // initialisation des quantites de reestimation

      chain_reestim->init();

      for (i = 0;i < hmarkov->nb_output_process;i++) {
        if (observation_reestim[i]) {
          for (j = 0;j < hmarkov->nb_state;j++) {
            for (k = 0;k < marginal_distribution[i]->nb_value;k++) {
              observation_reestim[i][j]->frequency[k] = 0.;
            }
          }
        }
      }

      if (state_sequence_count) {
        for (i = 0;i < nb_sequence;i++) {
          for (j = 0;j < length[i];j++) {
            for (k = 0;k < hmarkov->nb_state;k++) {
              state_sequence_count[i][j][k] = 0.;
            }
          }
        }
      }

      for (i = 0;i < nb_sequence;i++) {

        // recurrence "forward"

        for (j = 0;j < nb_variable;j++) {
          switch (type[j]) {
          case INT_VALUE :
            pioutput[j] = int_sequence[i][j];
            break;
          case REAL_VALUE :
            proutput[j] = real_sequence[i][j];
            break;
          }
        }

        norm = 0.;

        switch (hmarkov->type) {

        case 'o' : {
          for (j = 1;j < hmarkov->nb_row;j++) {
            if (hmarkov->order[j] == 1) {
              forward[0][j] = hmarkov->initial[hmarkov->state[j][0]];

              for (k = 0;k < hmarkov->nb_output_process;k++) {
                if (hmarkov->categorical_process[k]) {
                  forward[0][j] *= hmarkov->categorical_process[k]->observation[hmarkov->state[j][0]]->mass[*pioutput[k]];
                }

                else if (hmarkov->discrete_parametric_process[k]) {
                  forward[0][j] *= hmarkov->discrete_parametric_process[k]->observation[hmarkov->state[j][0]]->mass[*pioutput[k]];
                }

                else {
                  if (((hmarkov->continuous_parametric_process[k]->ident == GAMMA) ||
                      (hmarkov->continuous_parametric_process[k]->ident == ZERO_INFLATED_GAMMA)) && (min_value[k] < min_interval[k] / 2)) {
                    switch (type[k]) {
                    case INT_VALUE :
                      forward[0][j] *= hmarkov->continuous_parametric_process[k]->observation[hmarkov->state[j][0]]->mass_computation(*pioutput[k] , *pioutput[k] + min_interval[k]);
                      break;
                    case REAL_VALUE :
                      forward[0][j] *= hmarkov->continuous_parametric_process[k]->observation[hmarkov->state[j][0]]->mass_computation(*proutput[k] , *proutput[k] + min_interval[k]);
                      break;
                    }
                  }

                  else {
                    switch (type[k]) {
                    case INT_VALUE :
                      forward[0][j] *= hmarkov->continuous_parametric_process[k]->observation[hmarkov->state[j][0]]->mass_computation(*pioutput[k] - min_interval[k] / 2 , *pioutput[k] + min_interval[k] / 2);
                      break;
                    case REAL_VALUE :
                      forward[0][j] *= hmarkov->continuous_parametric_process[k]->observation[hmarkov->state[j][0]]->mass_computation(*proutput[k] - min_interval[k] / 2 , *proutput[k] + min_interval[k] / 2);
                      break;
                    }
                  }
                }
              }

              norm += forward[0][j];
            }

            else {
              forward[0][j] = 0.;
            }
          }
          break;
        }

        case 'e' : {
          for (j = 1;j < hmarkov->nb_row;j++) {
            if (!(hmarkov->child[j])) {
              forward[0][j] = hmarkov->initial[j];

              for (k = 0;k < hmarkov->nb_output_process;k++) {
                if (hmarkov->categorical_process[k]) {
                  forward[0][j] *= hmarkov->categorical_process[k]->observation[hmarkov->state[j][0]]->mass[*pioutput[k]];
                }

                else if (hmarkov->discrete_parametric_process[k]) {
                  forward[0][j] *= hmarkov->discrete_parametric_process[k]->observation[hmarkov->state[j][0]]->mass[*pioutput[k]];
                }

                else {
                  if (((hmarkov->continuous_parametric_process[k]->ident == GAMMA) ||
                      (hmarkov->continuous_parametric_process[k]->ident == ZERO_INFLATED_GAMMA)) && (min_value[k] < min_interval[k] / 2)) {
                    switch (type[k]) {
                    case INT_VALUE :
                      forward[0][j] *= hmarkov->continuous_parametric_process[k]->observation[hmarkov->state[j][0]]->mass_computation(*pioutput[k] , *pioutput[k] + min_interval[k]);
                      break;
                    case REAL_VALUE :
                      forward[0][j] *= hmarkov->continuous_parametric_process[k]->observation[hmarkov->state[j][0]]->mass_computation(*proutput[k] , *proutput[k] + min_interval[k]);
                      break;
                    }
                  }

                  else {
                    switch (type[k]) {
                    case INT_VALUE :
                      forward[0][j] *= hmarkov->continuous_parametric_process[k]->observation[hmarkov->state[j][0]]->mass_computation(*pioutput[k] - min_interval[k] / 2 , *pioutput[k] + min_interval[k] / 2);
                      break;
                    case REAL_VALUE :
                      forward[0][j] *= hmarkov->continuous_parametric_process[k]->observation[hmarkov->state[j][0]]->mass_computation(*proutput[k] - min_interval[k] / 2 , *proutput[k] + min_interval[k] / 2);
                      break;
                    }
                  }
                }
              }

              norm += forward[0][j];
            }

            else {
              forward[0][j] = 0.;
            }
          }
          break;
        }
        }

        if (norm > 0.) {
          for (j = 1;j < hmarkov->nb_row;j++) {
            forward[0][j] /= norm;
          }

          likelihood += log(norm);
        }

        else {
          likelihood = D_INF;
          break;
        }

        for (j = 1;j < length[i];j++) {
          for (k = 0;k < nb_variable;k++) {
            switch (type[k]) {
            case INT_VALUE :
              pioutput[k]++;
              break;
            case REAL_VALUE :
              proutput[k]++;
              break;
            }
          }
          norm = 0.;

          for (k = 1;k < hmarkov->nb_row;k++) {
            forward[j][k] = 0.;
            for (m = 0;m < hmarkov->nb_memory[k];m++) {
              forward[j][k] += hmarkov->transition[hmarkov->previous[k][m]][hmarkov->state[k][0]] *
                               forward[j - 1][hmarkov->previous[k][m]];
            }
            predicted[j][k] = forward[j][k];

            for (m = 0;m < hmarkov->nb_output_process;m++) {
              if (hmarkov->categorical_process[m]) {
                forward[j][k] *= hmarkov->categorical_process[m]->observation[hmarkov->state[k][0]]->mass[*pioutput[m]];
              }

              else if (hmarkov->discrete_parametric_process[m]) {
                forward[j][k] *= hmarkov->discrete_parametric_process[m]->observation[hmarkov->state[k][0]]->mass[*pioutput[m]];
              }

              else {
                if (((hmarkov->continuous_parametric_process[m]->ident == GAMMA) ||
                    (hmarkov->continuous_parametric_process[m]->ident == ZERO_INFLATED_GAMMA)) && (min_value[m] < min_interval[m] / 2)) {
                  switch (type[m]) {
                  case INT_VALUE :
                    forward[j][k] *= hmarkov->continuous_parametric_process[m]->observation[hmarkov->state[k][0]]->mass_computation(*pioutput[m] , *pioutput[m] + min_interval[m]);
                    break;
                  case REAL_VALUE :
                    forward[j][k] *= hmarkov->continuous_parametric_process[m]->observation[hmarkov->state[k][0]]->mass_computation(*proutput[m] , *proutput[m] + min_interval[m]);
                    break;
                  }
                }

                else {
                  switch (type[m]) {
                  case INT_VALUE :
                    forward[j][k] *= hmarkov->continuous_parametric_process[m]->observation[hmarkov->state[k][0]]->mass_computation(*pioutput[m] - min_interval[m] / 2 , *pioutput[m] + min_interval[m] / 2);
                    break;
                  case REAL_VALUE :
                    forward[j][k] *= hmarkov->continuous_parametric_process[m]->observation[hmarkov->state[k][0]]->mass_computation(*proutput[m] - min_interval[m] / 2 , *proutput[m] + min_interval[m] / 2);
                    break;
                  }
                }
              }
            }

            norm += forward[j][k];
          }

          if (norm > 0.) {
            for (k = 1;k < hmarkov->nb_row;k++) {
              forward[j][k] /= norm;
            }
            likelihood += log(norm);
          }

          else {
            likelihood = D_INF;
            break;
          }
        }

        if (likelihood == D_INF) {
          break;
        }

        // recurrence "backward"

        j = length[i] - 1;
        for (k = 1;k < hmarkov->nb_row;k++) {
          backward[k] = forward[j][k];

          // accumulation des quantites de reestimation des lois d'observation

          for (m = 0;m < hmarkov->nb_output_process;m++) {
            if (observation_reestim[m]) {
              observation_reestim[m][hmarkov->state[k][0]]->frequency[*pioutput[m]] += backward[k];
            }
          }

          if (state_sequence_count) {
            state_sequence_count[i][j][hmarkov->state[k][0]] += backward[k];
          }
        }

        for (j = length[i] - 2;j >= 0;j--) {
          for (k = 0;k < nb_variable;k++) {
            if (type[k] == INT_VALUE) {
              pioutput[k]--;
            }
          }

          for (k = 1;k < hmarkov->nb_row;k++) {
            if (predicted[j + 1][k] > 0.) {
              auxiliary[k] = backward[k] / predicted[j + 1][k];
            }
            else {
              auxiliary[k] = 0.;
            }
          }

          for (k = 1;k < hmarkov->nb_row;k++) {
            backward[k] = 0.;

            if (hmarkov->next[k]) {
              for (m = 0;m < hmarkov->nb_state;m++) {
                buff = auxiliary[hmarkov->next[k][m]] * hmarkov->transition[k][m] * forward[j][k];
                backward[k] += buff;

                // accumulation des quantites de reestimation des probabilites de transition

                chain_reestim->transition[k][m] += buff;
              }

              // accumulation des quantites de reestimation des lois d'observation

              for (m = 0;m < hmarkov->nb_output_process;m++) {
                if (observation_reestim[m]) {
                  observation_reestim[m][hmarkov->state[k][0]]->frequency[*pioutput[m]] += backward[k];
                }
              }

              if (state_sequence_count) {
                state_sequence_count[i][j][hmarkov->state[k][0]] += backward[k];
              }
            }
          }

#         ifdef DEBUG
/*          cout << j << " : ";
          sum = 0.;
          for (k = 1;k < hmarkov->nb_row;k++) {
            sum += backward[k];
            cout << backward[k] << " ";
          }
          cout << "| " << sum << endl; */
#         endif

        }

        // accumulation des quantites de reestimation des probabilites initiales

        if (hmarkov->type == 'o') {
          for (j = 1;j < hmarkov->nb_row;j++) {
            if (hmarkov->order[j] == 1) {
              chain_reestim->initial[hmarkov->state[j][0]] += backward[j];
            }
          }
        }
      }

      if (likelihood != D_INF) {

        // reestimation des probabilites initiales

        if (hmarkov->type == 'o') {
          reestimation(hmarkov->nb_state , chain_reestim->initial ,
                       hmarkov->initial , MIN_PROBABILITY , false);
        }

        // reestimation des probabilites de transition

        for (i = hmarkov->nb_row - 1;i >= 1;i--) {
          if (hmarkov->memory_type[i] == COMPLETION) {
/*          if ((hmarkov->memory_type[i] == COMPLETION) || ((hmarkov->type == 'o') &&
                 (global_initial_transition) && (hmarkov->order[i] > 1))) { */
            for (j = 0;j < hmarkov->nb_state;j++) {
              chain_reestim->transition[hmarkov->parent[i]][j] += chain_reestim->transition[i][j];
            }
          }
        }

        for (i = 1;i < hmarkov->nb_row;i++) {
          if ((hmarkov->memory_type[i] == TERMINAL) || ((hmarkov->type == 'o') &&
               (hmarkov->memory_type[i] == NON_TERMINAL))) {
            reestimation(hmarkov->nb_state , chain_reestim->transition[i] ,
                         hmarkov->transition[i] , MIN_PROBABILITY , false);
          }
          else if (hmarkov->memory_type[i] == COMPLETION) {
            for (j = 0;j < hmarkov->nb_state;j++) {
              hmarkov->transition[i][j] = hmarkov->transition[hmarkov->parent[i]][j];
            }
          }
        }

        if (hmarkov->type == 'e') {
          hmarkov->initial_probability_computation();
        }

        // reestimation des lois d'observation

        for (i = 0;i < hmarkov->nb_output_process;i++) {
          if (hmarkov->categorical_process[i]) {
            for (j = 0;j < hmarkov->nb_state;j++) {
              reestimation(marginal_distribution[i]->nb_value , observation_reestim[i][j]->frequency ,
                           hmarkov->categorical_process[i]->observation[j]->mass ,
                           MIN_PROBABILITY , false);
            }
          }

          else if (observation_reestim[i]) {
            for (j = 0;j < hmarkov->nb_state;j++) {
              observation_reestim[i][j]->nb_value_computation();
              observation_reestim[i][j]->offset_computation();
              observation_reestim[i][j]->nb_element_computation();
              observation_reestim[i][j]->max_computation();
              if ((hmarkov->discrete_parametric_process[i]) ||
                  (hmarkov->continuous_parametric_process[i]->ident != ZERO_INFLATED_GAMMA)) {
                observation_reestim[i][j]->mean_computation();
                observation_reestim[i][j]->variance_computation(true);
//                observation_reestim[i][j]->variance_computation();
              }
            }

            if (hmarkov->discrete_parametric_process[i]) {
              for (j = 0;j < hmarkov->nb_state;j++) {
                hobservation->update(observation_reestim[i][j] ,
                                     MAX((int)(observation_reestim[i][j]->nb_element *
                                               MAX(sqrt(observation_reestim[i][j]->variance) , 1.) * OBSERVATION_COEFF) , MIN_NB_ELEMENT));
                observation_likelihood = hobservation->Reestimation<int>::type_parametric_estimation(hmarkov->discrete_parametric_process[i]->observation[j] ,
                                                                                                     0 , true , OBSERVATION_THRESHOLD);

                if (observation_likelihood != D_INF) {
                  hmarkov->discrete_parametric_process[i]->observation[j]->computation(marginal_distribution[i]->nb_value ,
                                                                                       OBSERVATION_THRESHOLD);

                  if (hmarkov->discrete_parametric_process[i]->observation[j]->ident == BINOMIAL) {
                    for (k = hmarkov->discrete_parametric_process[i]->observation[j]->nb_value;k < marginal_distribution[i]->nb_value;k++) {
                      hmarkov->discrete_parametric_process[i]->observation[j]->mass[k] = 0.;
                    }
                  }
                }
              }
            }

            else {
              switch (hmarkov->continuous_parametric_process[i]->ident) {

              case GAMMA : {
                for (j = 0;j < hmarkov->nb_state;j++) {
                  observation_reestim[i][j]->gamma_estimation(hmarkov->continuous_parametric_process[i]->observation[j] , iter);
                }
                break;
              }

              case ZERO_INFLATED_GAMMA : {
                for (j = 0;j < hmarkov->nb_state;j++) {
                  observation_reestim[i][j]->zero_inflated_gamma_estimation(hmarkov->continuous_parametric_process[i]->observation[j] , iter);
                }
                break;
              }

              case GAUSSIAN : {
                for (j = 0;j < hmarkov->nb_state;j++) {
                  hmarkov->continuous_parametric_process[i]->observation[j]->location = observation_reestim[i][j]->mean;
                }

                switch (common_dispersion) {

                case false : {
                  for (j = 0;j < hmarkov->nb_state;j++) {
                    hmarkov->continuous_parametric_process[i]->observation[j]->dispersion = sqrt(observation_reestim[i][j]->variance);
                    if (hmarkov->continuous_parametric_process[i]->observation[j]->dispersion /
                        hmarkov->continuous_parametric_process[i]->observation[j]->location < GAUSSIAN_MIN_VARIATION_COEFF) {
                      hmarkov->continuous_parametric_process[i]->observation[j]->dispersion = hmarkov->continuous_parametric_process[i]->observation[j]->location * GAUSSIAN_MIN_VARIATION_COEFF;
                    }
                  }
                  break;
                }

                case true : {
                  variance = 0.;
                  buff = 0.;

                  for (j = 0;j < hmarkov->nb_state;j++) {
                    for (k = observation_reestim[i][j]->offset;k < observation_reestim[i][j]->nb_value;k++) {
                      diff = k - observation_reestim[i][j]->mean;
                      variance += observation_reestim[i][j]->frequency[k] * diff * diff;
                    }

                    buff += observation_reestim[i][j]->nb_element;
                  }

                  variance /= buff;
//                  variance /= (buff - 1);

                  for (j = 0;j < hmarkov->nb_state;j++) {
                    hmarkov->continuous_parametric_process[i]->observation[j]->dispersion = sqrt(variance);
                  }
                  break;
                }
                }

                break;
              }

              case VON_MISES : {
                for (j = 0;j < hmarkov->nb_state;j++) {
                  observation_reestim[i][j]->mean_direction_computation(mean_direction[j]);
                  hmarkov->continuous_parametric_process[i]->observation[j]->location = mean_direction[j][3];
                }

                switch (common_dispersion) {

                case false : {
                  for (j = 0;j < hmarkov->nb_state;j++) {
                    hmarkov->continuous_parametric_process[i]->observation[j]->dispersion = von_mises_concentration_computation(mean_direction[j][2]);
                  }
                  break;
                }

                case true : {
                  global_mean_direction = 0.;
                  buff = 0.;

                  for (j = 0;j < hmarkov->nb_state;j++) {
                    global_mean_direction += observation_reestim[i][j]->nb_element * mean_direction[j][2];
                    buff += observation_reestim[i][j]->nb_element;
                  }
                  concentration = von_mises_concentration_computation(global_mean_direction / buff);

                  for (j = 0;j < hmarkov->nb_state;j++) {
                    hmarkov->continuous_parametric_process[i]->observation[j]->dispersion = concentration;
                  }
                  break;
                }
                }
                break;
              }
              }
            }
          }

          else {
            switch (hmarkov->continuous_parametric_process[i]->ident) {
            case GAMMA :
              gamma_estimation(state_sequence_count , i ,
                               hmarkov->continuous_parametric_process[i] , iter);
              break;
            case ZERO_INFLATED_GAMMA :
              zero_inflated_gamma_estimation(state_sequence_count , i ,
                                             hmarkov->continuous_parametric_process[i] , iter);
              break;
            case GAUSSIAN :
              gaussian_estimation(state_sequence_count , i ,
                                  hmarkov->continuous_parametric_process[i]);
              break;
            case VON_MISES :
              von_mises_estimation(state_sequence_count , i ,
                                   hmarkov->continuous_parametric_process[i]);
              break;
            }
          }
        }
      }

#     ifdef MESSAGE
      os << STAT_label[STATL_ITERATION] << " " << iter << "   "
         << SEQ_label[SEQL_OBSERVED_SEQUENCES_LIKELIHOOD] << ": " << likelihood << endl;
#     endif

#     ifdef DEBUG
      if (iter % 5 == 0) {
        cout << *hmarkov;
      }
#     endif

    }
    while ((likelihood != D_INF) && (((nb_iter == I_DEFAULT) && (iter < VARIABLE_ORDER_MARKOV_NB_ITER) &&
             ((likelihood - previous_likelihood) / -likelihood > VARIABLE_ORDER_MARKOV_LIKELIHOOD_DIFF)) ||
            ((nb_iter != I_DEFAULT) && (iter < nb_iter))));

    if (likelihood != D_INF) {

#     ifdef MESSAGE
      os << "\n" << iter << " " << STAT_label[STATL_ITERATIONS] << endl;
#     endif

      // reestimation des probabilites initiales

      if (hmarkov->type == 'o') {
        reestimation(hmarkov->nb_state , chain_reestim->initial ,
                     hmarkov->initial , MIN_PROBABILITY , true);
      }

      // reestimation des probabilites de transition

      if ((hmarkov->type == 'o') && (global_initial_transition)) {
        for (i = hmarkov->nb_row - 1;i >= 1;i--) {
          if ((hmarkov->memory_type[i] != COMPLETION) && (hmarkov->order[i] > 1)) {
            for (j = 0;j < hmarkov->nb_state;j++) {
              chain_reestim->transition[hmarkov->parent[i]][j] += chain_reestim->transition[i][j];
            }
          }
        }
      }

      for (i = 1;i < hmarkov->nb_row;i++) {
        if ((hmarkov->memory_type[i] == TERMINAL) || ((hmarkov->type == 'o') &&
             (hmarkov->memory_type[i] == NON_TERMINAL))) {
          reestimation(hmarkov->nb_state , chain_reestim->transition[i] ,
                       hmarkov->transition[i] , MIN_PROBABILITY , true);
        }
        else if (hmarkov->memory_type[i] == COMPLETION) {
          for (j = 0;j < hmarkov->nb_state;j++) {
            hmarkov->transition[i][j] = hmarkov->transition[hmarkov->parent[i]][j];
          }
        }
      }

      if (hmarkov->type == 'e') {
        hmarkov->initial_probability_computation();
      }

      // reestimation des lois d'observation categorielles

      for (i = 0;i < hmarkov->nb_output_process;i++) {
        if (hmarkov->categorical_process[i]) {
          for (j = 0;j < hmarkov->nb_state;j++) {
            reestimation(marginal_distribution[i]->nb_value , observation_reestim[i][j]->frequency ,
                         hmarkov->categorical_process[i]->observation[j]->mass ,
                         MIN_PROBABILITY , true);
          }
        }

        else if (hmarkov->discrete_parametric_process[i]) {
          hmarkov->discrete_parametric_process[i]->nb_value_computation();
        }
      }
    }

    // destruction des structures de donnees de l'algorithme

    for (i = 0;i < max_length;i++) {
      delete [] forward[i];
    }
    delete [] forward;

    for (i = 0;i < max_length;i++) {
      delete [] predicted[i];
    }
    delete [] predicted;

    delete [] backward;

    delete [] auxiliary;

    delete chain_reestim;

    for (i = 0;i < hmarkov->nb_output_process;i++) {
      if (observation_reestim[i]) {
        for (j = 0;j < hmarkov->nb_state;j++) {
          delete observation_reestim[i][j];
        }
        delete [] observation_reestim[i];
      }
    }
    delete [] observation_reestim;

    delete hobservation;

    if (state_sequence_count) {
      for (i = 0;i < nb_sequence;i++) {
        for (j = 0;j < length[i];j++) {
          delete [] state_sequence_count[i][j];
        }
        delete [] state_sequence_count[i];
      }
      delete [] state_sequence_count;
    }

    if (mean_direction) {
      for (i = 0;i < hmarkov->nb_state;i++) {
        delete [] mean_direction[i];
      }
      delete [] mean_direction;
    }

    delete [] pioutput;
    delete [] proutput;

    if (likelihood == D_INF) {
      delete hmarkov;
      hmarkov = NULL;
      error.update(STAT_error[STATR_ESTIMATION_FAILURE]);
    }

    else {
      if (state_sequence) {
        hmarkov->markov_data = new VariableOrderMarkovData(*this , 'a' , (hmarkov->type == 'e' ? true : false));
        seq = hmarkov->markov_data;

        for (i = 0;i < hmarkov->nb_output_process;i++) {
          if (((hmarkov->discrete_parametric_process[i]) || (hmarkov->continuous_parametric_process[i])) &&
              (seq->characteristics[i + 1])) {
            delete seq->characteristics[i + 1];
            seq->characteristics[i + 1] = NULL;
          }
        }

        hmarkov->forward_backward(*seq);

        hmarkov->create_cumul();
        hmarkov->log_computation();
        hmarkov->viterbi(*seq);
        hmarkov->remove_cumul();

        seq->min_value_computation(0);
        seq->max_value_computation(0);
        seq->build_marginal_frequency_distribution(0);
        seq->build_characteristic(0);

        seq->build_transition_count(*hmarkov);
        seq->build_observation_frequency_distribution(hmarkov->nb_state);
        seq->build_observation_histogram(hmarkov->nb_state);

        // calcul des melanges de lois d'observation (poids deduits de la restauration)

        weight = seq->weight_computation();

        for (i = 0;i < hmarkov->nb_output_process;i++) {
          if (hmarkov->categorical_process[i]) {
            hmarkov->categorical_process[i]->restoration_weight = new Distribution(*weight);
            hmarkov->categorical_process[i]->restoration_mixture = hmarkov->categorical_process[i]->mixture_computation(hmarkov->categorical_process[i]->restoration_weight);
          }

          else if (hmarkov->discrete_parametric_process[i]) {
            for (j = 0;j < hmarkov->nb_state;j++) {
              hmarkov->discrete_parametric_process[i]->observation[j]->cumul_computation();
            }

            hmarkov->discrete_parametric_process[i]->restoration_weight = new Distribution(*weight);
            hmarkov->discrete_parametric_process[i]->restoration_mixture = hmarkov->discrete_parametric_process[i]->mixture_computation(hmarkov->discrete_parametric_process[i]->restoration_weight);
          }

          else if (hmarkov->continuous_parametric_process[i]) {
            hmarkov->continuous_parametric_process[i]->restoration_weight = new Distribution(*weight);
          }
        }

        delete weight;

#       ifdef MESSAGE
        os << "\n" << SEQ_label[SEQL_STATE_SEQUENCES_LIKELIHOOD] << ": " << seq->restoration_likelihood;

        for (i = 0;i < nb_variable;i++) {
          if (type[i] == REAL_VALUE) {
            break;
          }
        }
        if (i == nb_variable) {
          os << " | " << hmarkov->VariableOrderMarkov::likelihood_computation(*seq);
        }
        os << endl;
#       endif

      }

      else {
        hmarkov->markov_data = new VariableOrderMarkovData(*this , 'c' , (hmarkov->type == 'e' ? true : false));
        seq = hmarkov->markov_data;
        if (seq->type[0] == STATE) {
          seq->state_variable_init(INT_VALUE);
        }

        for (i = 0;i < hmarkov->nb_output_process;i++) {
          if (((hmarkov->discrete_parametric_process[i]) || (hmarkov->continuous_parametric_process[i])) &&
              (seq->characteristics[i])) {
            delete seq->characteristics[i];
            seq->characteristics[i] = NULL;
          }
        }
      }

      for (i = 0;i < hmarkov->nb_output_process;i++) {
        if (hmarkov->categorical_process[i]) {
          for (j = 0;j < hmarkov->nb_state;j++) {
            hmarkov->categorical_process[i]->observation[j]->cumul_computation();

            hmarkov->categorical_process[i]->observation[j]->max_computation();
//            hmarkov->categorical_process[i]->observation[j]->mean_computation();
//            hmarkov->categorical_process[i]->observation[j]->variance_computation();
          }
        }
      }

      // calcul de la vraisemblance et des lois caracteristiques du modele

      seq->likelihood = hmarkov->likelihood_computation(*this , seq->posterior_probability);

#     ifdef DEBUG
//      cout << *hmarkov;
      cout << "iteration " << iter << "  "
           << SEQ_label[SEQL_OBSERVED_SEQUENCES_LIKELIHOOD] << ": " << seq->likelihood << endl;
#     endif

#     ifdef MESSAGE
      if  ((state_sequence) && (seq->nb_sequence <= POSTERIOR_PROBABILITY_NB_SEQUENCE)) {
        os << "\n" << SEQ_label[SEQL_POSTERIOR_STATE_SEQUENCE_PROBABILITY] << endl;
        for (i = 0;i < seq->nb_sequence;i++) {
          os << SEQ_label[SEQL_SEQUENCE] << " " << seq->identifier[i] << ": "
             << seq->posterior_probability[i] << endl;
        }
      }
#     endif

      hmarkov->component_computation();
      hmarkov->characteristic_computation(*seq , counting_flag , I_DEFAULT , false);

      // calcul des melanges de lois d'observation (poids theoriques)

      switch (hmarkov->type) {

      case 'o' : {
        weight = hmarkov->state_process->weight_computation();
        break;
      }

      case 'e' : {
        weight = new Distribution(hmarkov->nb_state);

        for (i = 0;i < hmarkov->nb_state;i++) {
          weight->mass[i] = 0.;
        }
        for (i = 1;i < hmarkov->nb_row;i++) {
          if ((hmarkov->memory_type[i] == TERMINAL) || (hmarkov->memory_type[i] == COMPLETION)) {
            weight->mass[hmarkov->state[i][0]] += hmarkov->initial[i];
          }
        }

        weight->cumul_computation();
        weight->max_computation();
        break;
      }
      }

      for (i = 0;i < hmarkov->nb_output_process;i++) {
        if (hmarkov->categorical_process[i]) {
          hmarkov->categorical_process[i]->weight = new Distribution(*weight);
          hmarkov->categorical_process[i]->mixture = hmarkov->categorical_process[i]->mixture_computation(hmarkov->categorical_process[i]->weight);
        }

        else if (hmarkov->discrete_parametric_process[i]) {
          hmarkov->discrete_parametric_process[i]->weight = new Distribution(*weight);
          hmarkov->discrete_parametric_process[i]->mixture = hmarkov->discrete_parametric_process[i]->mixture_computation(hmarkov->discrete_parametric_process[i]->weight);
        }

        else if (hmarkov->continuous_parametric_process[i]) {
          hmarkov->continuous_parametric_process[i]->weight = new Distribution(*weight);
        }
      }

      delete weight;
    }
  }

  return hmarkov;
}


/*--------------------------------------------------------------*
 *
 *  Estimation des parametres d'une chaine de Markov d'ordre variable cachee
 *  a partir d'un echantillon de sequences par l'algorithme MCEM.
 *
 *  arguments : reference sur un objet StatError, stream, chaine de Markov cachee initiale,
 *              type d'estimation des probabilites de transition initiale (cas ordinaire),
 *              flag parametres de dispersion communs (processus d'observation continus),
 *              parametres pour le nombre de sequences d'etats simulees, flags sur le calcul
 *              des lois de comptage et sur le calcul des sequences d'etats optimales,
 *              nombre d'iterations.
 *
 *--------------------------------------------------------------*/

HiddenVariableOrderMarkov* MarkovianSequences::hidden_variable_order_markov_stochastic_estimation(StatError &error , ostream &os ,
                                                                                                  const HiddenVariableOrderMarkov &ihmarkov ,
                                                                                                  bool global_initial_transition ,
                                                                                                  bool common_dispersion ,
                                                                                                  int min_nb_state_sequence ,
                                                                                                  int max_nb_state_sequence ,
                                                                                                  double parameter , bool counting_flag ,
                                                                                                  bool state_sequence , int nb_iter) const

{
  bool status;
  register int i , j , k , m;
  int nb_terminal , iter , nb_state_sequence , memory , *state_seq , *pstate ,
      ***state_sequence_count , nb_element , **pioutput;
  double likelihood = D_INF , previous_likelihood , observation_likelihood , **forward ,
         norm , **predicted , *backward , *cumul_backward , diff , variance ,
         **mean_direction , concentration , global_mean_direction , **proutput;
  Distribution *weight;
  ChainReestimation<double> *chain_reestim;
  Reestimation<double> ***observation_reestim;
  HiddenVariableOrderMarkov *hmarkov;
  VariableOrderMarkovData *seq;

# ifdef DEBUG
  double sum;
# endif


  hmarkov = NULL;
  error.init();

  // test nombre de valeurs observees par variable

  status = false;
  for (i = 0;i < nb_variable;i++) {
    if (max_value[i] > min_value[i]) {
      status = true;
      break;
    }
  }

  if (!status) {
    error.update(STAT_error[STATR_VARIABLE_NB_VALUE]);
  }

  for (i = 0;i < nb_variable;i++) {
    if ((type[i] != INT_VALUE) && (type[i] != REAL_VALUE) && (type[i] != STATE)) {
      status = false;
      ostringstream error_message , correction_message;
      error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                    << STAT_error[STATR_VARIABLE_TYPE];
      correction_message << STAT_variable_word[INT_VALUE] << " or "
                         << STAT_variable_word[STATE];
      error.correction_update((error_message.str()).c_str() , (correction_message.str()).c_str());
    }
  }

  if (ihmarkov.nb_output_process != nb_variable) {
    status = false;
    error.update(STAT_error[STATR_NB_OUTPUT_PROCESS]);
  }

  else {
    for (i = 0;i < nb_variable;i++) {
      if ((ihmarkov.categorical_process[i]) || (ihmarkov.discrete_parametric_process[i])) {
        if (type[i] == REAL_VALUE) {
          status = false;
          ostringstream error_message;
          error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                        << STAT_error[STATR_VARIABLE_TYPE];
          error.correction_update((error_message.str()).c_str() , STAT_variable_word[INT_VALUE]);
        }

        else {
          if (min_value[i] < 0) {
            status = false;
            ostringstream error_message;
            error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                          << STAT_error[STATR_POSITIVE_MIN_VALUE];
            error.update((error_message.str()).c_str());
          }

          if (!marginal_distribution[i]) {
            status = false;
            ostringstream error_message;
            error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                          << STAT_error[STATR_MARGINAL_FREQUENCY_DISTRIBUTION];
            error.update((error_message.str()).c_str());
          }

          else {
            if (((ihmarkov.categorical_process[i]) &&
                 (ihmarkov.categorical_process[i]->nb_value != marginal_distribution[i]->nb_value)) ||
                ((ihmarkov.discrete_parametric_process[i]) &&
                 (ihmarkov.discrete_parametric_process[i]->nb_value < marginal_distribution[i]->nb_value))) {
              status = false;
              ostringstream error_message;
              error_message << STAT_label[STATL_OUTPUT_PROCESS] << " " << i + 1 << ": "
                            << STAT_error[STATR_NB_OUTPUT];
              error.update((error_message.str()).c_str());
            }

            else if ((ihmarkov.categorical_process[i]) && (!characteristics[i])) {
              for (j = 0;j < marginal_distribution[i]->nb_value;j++) {
                if (marginal_distribution[i]->frequency[j] == 0) {
                  status = false;
                  ostringstream error_message;
                  error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                                << STAT_error[STATR_MISSING_VALUE] << " " << j;
                  error.update((error_message.str()).c_str());
                }
              }
            }
          }
        }
      }
    }
  }

  if ((min_nb_state_sequence < 1) || (min_nb_state_sequence > max_nb_state_sequence)) {
    status = false;
    error.update(SEQ_error[SEQR_MIN_NB_STATE_SEQUENCE]);
  }

  if ((nb_iter != I_DEFAULT) && (nb_iter < 1)) {
    status = false;
    error.update(STAT_error[STATR_NB_ITERATION]);
  }

  if (status) {

    // creation de la chaine de Markov cachee

    hmarkov = new HiddenVariableOrderMarkov(ihmarkov , false);

    if (hmarkov->type == 'e') {
      nb_terminal = (hmarkov->nb_row - 1) * (hmarkov->nb_state - 1) / hmarkov->nb_state + 1;

      for (i = 1;i < hmarkov->nb_row;i++) {
        if (!hmarkov->child[i]) {
          hmarkov->initial[i] = 1. / (double)nb_terminal;
        }
        else {
          hmarkov->initial[i] = 0.;
        }
      }
    }

    if (common_dispersion) {
      for (i = 0;i < hmarkov->nb_output_process;i++) {
        if (hmarkov->continuous_parametric_process[i]) {
          hmarkov->continuous_parametric_process[i]->tied_dispersion = true;
        }
      }
    }

#   ifdef DEBUG
    cout << *hmarkov;
#   endif

    // creation des structures de donnees de l'algorithme

    forward = new double*[max_length];
    for (i = 0;i < max_length;i++) {
      forward[i] = new double[hmarkov->nb_row];
    }

    predicted = new double*[max_length];
    for (i = 0;i < max_length;i++) {
      predicted[i] = new double[hmarkov->nb_row];
    }

    backward = new double[hmarkov->nb_row];
    cumul_backward = new double[hmarkov->nb_row];

    state_seq = new int[max_length];

    chain_reestim = new ChainReestimation<double>((hmarkov->type == 'o' ?  'o' : 'v') ,
                                                  hmarkov->nb_state , hmarkov->nb_row);

    observation_reestim = new Reestimation<double>**[hmarkov->nb_output_process];
    for (i = 0;i < hmarkov->nb_output_process;i++) {
      if (marginal_distribution[i]) {
        observation_reestim[i] = new Reestimation<double>*[hmarkov->nb_state];
        for (j = 0;j < hmarkov->nb_state;j++) {
          observation_reestim[i][j] = new Reestimation<double>(marginal_distribution[i]->nb_value);
        }
      }

      else {
        observation_reestim[i] = NULL;
      }
    }

    for (i = 0;i < hmarkov->nb_output_process;i++) {
      if (!marginal_distribution[i]) {
        break;
      }
    }

    if (i < hmarkov->nb_output_process) {
      state_sequence_count = new int**[nb_sequence];
      for (i = 0;i < nb_sequence;i++) {
        state_sequence_count[i] = new int*[length[i]];
        for (j = 0;j < length[i];j++) {
          state_sequence_count[i][j] = new int[hmarkov->nb_state];
        }
      }
    }
    else {
      state_sequence_count = NULL;
    }

    for (i = 0;i < hmarkov->nb_output_process;i++) {
      if ((hmarkov->continuous_parametric_process[i]) &&
          (hmarkov->continuous_parametric_process[i]->ident == VON_MISES)) {
        break;
      }
    }

    if (i < hmarkov->nb_output_process) {
      mean_direction = new double*[hmarkov->nb_state];
      for (i = 0;i < hmarkov->nb_state;i++) {
        mean_direction[i] = new double[4];
      }
    }
    else {
      mean_direction = NULL;
    }

    pioutput = new int*[nb_variable];
    proutput = new double*[nb_variable];

    iter = 0;
    do {
      previous_likelihood = likelihood;
      likelihood = 0.;

      // calcul du nombre de sequences d'etats simulees

      if (min_nb_state_sequence + (int)::round(parameter * iter) < max_nb_state_sequence) {
        nb_state_sequence = min_nb_state_sequence + (int)::round(parameter * iter);
      }
      else {
        nb_state_sequence = max_nb_state_sequence;
      }

/*      nb_state_sequence = max_nb_state_sequence - (int)round((max_nb_state_sequence - min_nb_state_sequence) *
                          exp(-parameter * iter)); */

      iter++;

      // initialisation des quantites de reestimation

      chain_reestim->init();

      for (i = 0;i < hmarkov->nb_output_process;i++) {
        if (observation_reestim[i]) {
          for (j = 0;j < hmarkov->nb_state;j++) {
            for (k = 0;k < marginal_distribution[i]->nb_value;k++) {
              observation_reestim[i][j]->frequency[k] = 0.;
            }
          }
        }
      }

      if (state_sequence_count) {
        for (i = 0;i < nb_sequence;i++) {
          for (j = 0;j < length[i];j++) {
            for (k = 0;k < hmarkov->nb_state;k++) {
              state_sequence_count[i][j][k] = 0;
            }
          }
        }
      }

      for (i = 0;i < nb_sequence;i++) {

        // recurrence "forward"

        for (j = 0;j < nb_variable;j++) {
          switch (type[j]) {
          case INT_VALUE :
            pioutput[j] = int_sequence[i][j];
            break;
          case REAL_VALUE :
            proutput[j] = real_sequence[i][j];
            break;
          }
        }

        norm = 0.;

        switch (hmarkov->type) {

        case 'o' : {
          for (j = 1;j < hmarkov->nb_row;j++) {
            if (hmarkov->order[j] == 1) {
              forward[0][j] = hmarkov->initial[hmarkov->state[j][0]];

              for (k = 0;k < hmarkov->nb_output_process;k++) {
                if (hmarkov->categorical_process[k]) {
                  forward[0][j] *= hmarkov->categorical_process[k]->observation[hmarkov->state[j][0]]->mass[*pioutput[k]];
                }

                else if (hmarkov->discrete_parametric_process[k]) {
                  forward[0][j] *= hmarkov->discrete_parametric_process[k]->observation[hmarkov->state[j][0]]->mass[*pioutput[k]];
                }

                else {
                  if (((hmarkov->continuous_parametric_process[k]->ident == GAMMA) ||
                      (hmarkov->continuous_parametric_process[k]->ident == ZERO_INFLATED_GAMMA)) && (min_value[k] < min_interval[k] / 2)) {
                    switch (type[k]) {
                    case INT_VALUE :
                      forward[0][j] *= hmarkov->continuous_parametric_process[k]->observation[hmarkov->state[j][0]]->mass_computation(*pioutput[k] , *pioutput[k] + min_interval[k]);
                      break;
                    case REAL_VALUE :
                      forward[0][j] *= hmarkov->continuous_parametric_process[k]->observation[hmarkov->state[j][0]]->mass_computation(*proutput[k] , *proutput[k] + min_interval[k]);
                      break;
                    }
                  }

                  else {
                    switch (type[k]) {
                    case INT_VALUE :
                      forward[0][j] *= hmarkov->continuous_parametric_process[k]->observation[hmarkov->state[j][0]]->mass_computation(*pioutput[k] - min_interval[k] / 2 , *pioutput[k] + min_interval[k] / 2);
                      break;
                    case REAL_VALUE :
                      forward[0][j] *= hmarkov->continuous_parametric_process[k]->observation[hmarkov->state[j][0]]->mass_computation(*proutput[k] - min_interval[k] / 2 , *proutput[k] + min_interval[k] / 2);
                      break;
                    }
                  }
                }
              }

              norm += forward[0][j];
            }

            else {
              forward[0][j] = 0.;
            }
          }
          break;
        }

        case 'e' : {
          for (j = 1;j < hmarkov->nb_row;j++) {
            if (!(hmarkov->child[j])) {
              forward[0][j] = hmarkov->initial[j];

              for (k = 0;k < hmarkov->nb_output_process;k++) {
                if (hmarkov->categorical_process[k]) {
                  forward[0][j] *= hmarkov->categorical_process[k]->observation[hmarkov->state[j][0]]->mass[*pioutput[k]];
                }

                else if (hmarkov->discrete_parametric_process[k]) {
                  forward[0][j] *= hmarkov->discrete_parametric_process[k]->observation[hmarkov->state[j][0]]->mass[*pioutput[k]];
                }

                else {
                  if (((hmarkov->continuous_parametric_process[k]->ident == GAMMA) ||
                      (hmarkov->continuous_parametric_process[k]->ident == ZERO_INFLATED_GAMMA)) && (min_value[k] < min_interval[k] / 2)) {
                    switch (type[k]) {
                    case INT_VALUE :
                      forward[0][j] *= hmarkov->continuous_parametric_process[k]->observation[hmarkov->state[j][0]]->mass_computation(*pioutput[k] , *pioutput[k] + min_interval[k]);
                      break;
                    case REAL_VALUE :
                      forward[0][j] *= hmarkov->continuous_parametric_process[k]->observation[hmarkov->state[j][0]]->mass_computation(*proutput[k] , *proutput[k] + min_interval[k]);
                      break;
                    }
                  }

                  else {
                    switch (type[k]) {
                    case INT_VALUE :
                      forward[0][j] *= hmarkov->continuous_parametric_process[k]->observation[hmarkov->state[j][0]]->mass_computation(*pioutput[k] - min_interval[k] / 2 , *pioutput[k] + min_interval[k] / 2);
                      break;
                    case REAL_VALUE :
                      forward[0][j] *= hmarkov->continuous_parametric_process[k]->observation[hmarkov->state[j][0]]->mass_computation(*proutput[k] - min_interval[k] / 2 , *proutput[k] + min_interval[k] / 2);
                      break;
                    }
                  }
                }
              }

              norm += forward[0][j];
            }

            else {
              forward[0][j] = 0.;
            }
          }
          break;
        }
        }

        if (norm > 0.) {
          for (j = 1;j < hmarkov->nb_row;j++) {
            forward[0][j] /= norm;
          }

          likelihood += log(norm);
        }

        else {
          likelihood = D_INF;
          break;
        }

        for (j = 1;j < length[i];j++) {
          for (k = 0;k < nb_variable;k++) {
            switch (type[k]) {
            case INT_VALUE :
              pioutput[k]++;
              break;
            case REAL_VALUE :
              proutput[k]++;
              break;
            }
          }
          norm = 0.;

          for (k = 1;k < hmarkov->nb_row;k++) {
            forward[j][k] = 0.;
            for (m = 0;m < hmarkov->nb_memory[k];m++) {
              forward[j][k] += hmarkov->transition[hmarkov->previous[k][m]][hmarkov->state[k][0]] *
                               forward[j - 1][hmarkov->previous[k][m]];
            }
            predicted[j][k] = forward[j][k];

            for (m = 0;m < hmarkov->nb_output_process;m++) {
              if (hmarkov->categorical_process[m]) {
                forward[j][k] *= hmarkov->categorical_process[m]->observation[hmarkov->state[k][0]]->mass[*pioutput[m]];
              }

              else if (hmarkov->discrete_parametric_process[m]) {
                forward[j][k] *= hmarkov->discrete_parametric_process[m]->observation[hmarkov->state[k][0]]->mass[*pioutput[m]];
              }

              else {
                if (((hmarkov->continuous_parametric_process[m]->ident == GAMMA) ||
                    (hmarkov->continuous_parametric_process[m]->ident == ZERO_INFLATED_GAMMA)) && (min_value[m] < min_interval[m] / 2)) {
                  switch (type[m]) {
                  case INT_VALUE :
                    forward[j][k] *= hmarkov->continuous_parametric_process[m]->observation[hmarkov->state[k][0]]->mass_computation(*pioutput[m] , *pioutput[m] + min_interval[m]);
                    break;
                  case REAL_VALUE :
                    forward[j][k] *= hmarkov->continuous_parametric_process[m]->observation[hmarkov->state[k][0]]->mass_computation(*proutput[m] , *proutput[m] + min_interval[m]);
                    break;
                  }
                }

                else {
                  switch (type[m]) {
                  case INT_VALUE :
                    forward[j][k] *= hmarkov->continuous_parametric_process[m]->observation[hmarkov->state[k][0]]->mass_computation(*pioutput[m] - min_interval[m] / 2 , *pioutput[m] + min_interval[m] / 2);
                    break;
                  case REAL_VALUE :
                    forward[j][k] *= hmarkov->continuous_parametric_process[m]->observation[hmarkov->state[k][0]]->mass_computation(*proutput[m] - min_interval[m] / 2 , *proutput[m] + min_interval[m] / 2);
                    break;
                  }
                }
              }
            }

            norm += forward[j][k];
          }

          if (norm > 0.) {
            for (k = 1;k < hmarkov->nb_row;k++) {
              forward[j][k] /= norm;
            }
            likelihood += log(norm);
          }

          else {
            likelihood = D_INF;
            break;
          }
        }

        if (likelihood == D_INF) {
          break;
        }

        // passes "backward"

        for (j = 0;j < nb_state_sequence;j++) {
          k = length[i] - 1;
          pstate = state_seq + k;
          for (m = 0;m < nb_variable;m++) {
            if (type[m] == INT_VALUE) {
              pioutput[m] = int_sequence[i][m] + k;
            }
          }

          cumul_computation(hmarkov->nb_row - 1 , forward[k] + 1 , cumul_backward);
          memory = 1 + cumul_method(hmarkov->nb_row - 1 , cumul_backward);
          *pstate = hmarkov->state[memory][0];

          // accumulation des quantites de reestimation des lois d'observation

          for (m = 0;m < hmarkov->nb_output_process;m++) {
            if (observation_reestim[m]) {
              (observation_reestim[m][*pstate]->frequency[*pioutput[m]])++;
            }
          }

          if (state_sequence_count) {
            (state_sequence_count[i][k][*pstate])++;
          }

          for (k = length[i] - 2;k >= 0;k--) {
            for (m = 0;m < hmarkov->nb_memory[memory];m++) {
              backward[m] = hmarkov->transition[hmarkov->previous[memory][m]][hmarkov->state[memory][0]] *
                            forward[k][hmarkov->previous[memory][m]] / predicted[k + 1][memory];
            }

#           ifdef DEBUG
            sum = 0.;
            for (m = 0;m < hmarkov->nb_memory[memory];m++) {
              sum += backward[m];
            }
            if ((sum < 1. - DOUBLE_ERROR) || (sum > 1. + DOUBLE_ERROR)) {
              cout << "\nERROR: " << k << " " << sum << endl;
            }
#           endif

            cumul_computation(hmarkov->nb_memory[memory] , backward , cumul_backward);
            memory = hmarkov->previous[memory][cumul_method(hmarkov->nb_memory[memory] , cumul_backward)];
            *--pstate = hmarkov->state[memory][0];

            // accumulation des quantites de reestimation des probabilites de transition et
            // des lois d'observation

            (chain_reestim->transition[memory][*(pstate + 1)])++;

            for (m = 0;m < hmarkov->nb_output_process;m++) {
              if (observation_reestim[m]) {
                (observation_reestim[m][*pstate]->frequency[*--pioutput[m]])++;
              }
            }

            if (state_sequence_count) {
              (state_sequence_count[i][k][*pstate])++;
            }
          }

          // accumulation des quantites de reestimation des probabilites initiales

          if (hmarkov->type == 'o') {
            (chain_reestim->initial[*pstate])++;
          }
        }
      }

      if (likelihood != D_INF) {

        // reestimation des probabilites initiales

        if (hmarkov->type == 'o') {
          reestimation(hmarkov->nb_state , chain_reestim->initial ,
                       hmarkov->initial , MIN_PROBABILITY , false);
        }

        // reestimation des probabilites de transition

        for (i = hmarkov->nb_row - 1;i >= 1;i--) {
          if (hmarkov->memory_type[i] == COMPLETION) {
/*          if ((hmarkov->memory_type[i] == COMPLETION) || ((hmarkov->type == 'o') &&
                 (global_initial_transition) && (hmarkov->order[i] > 1))) { */
            for (j = 0;j < hmarkov->nb_state;j++) {
              chain_reestim->transition[hmarkov->parent[i]][j] += chain_reestim->transition[i][j];
            }
          }
        }

        for (i = 1;i < hmarkov->nb_row;i++) {
          if ((hmarkov->memory_type[i] == TERMINAL) || ((hmarkov->type == 'o') &&
               (hmarkov->memory_type[i] == NON_TERMINAL))) {
            reestimation(hmarkov->nb_state , chain_reestim->transition[i] ,
                         hmarkov->transition[i] , MIN_PROBABILITY , false);
          }
          else if (hmarkov->memory_type[i] == COMPLETION) {
            for (j = 0;j < hmarkov->nb_state;j++) {
              hmarkov->transition[i][j] = hmarkov->transition[hmarkov->parent[i]][j];
            }
          }
        }

        if (hmarkov->type == 'e') {
          hmarkov->initial_probability_computation();
        }

        // reestimation des lois d'observation

        for (i = 0;i < hmarkov->nb_output_process;i++) {
          if (hmarkov->categorical_process[i]) {
            for (j = 0;j < hmarkov->nb_state;j++) {
              reestimation(marginal_distribution[i]->nb_value , observation_reestim[i][j]->frequency ,
                           hmarkov->categorical_process[i]->observation[j]->mass ,
                           MIN_PROBABILITY , false);
            }
          }

          else if (observation_reestim[i]) {
            for (j = 0;j < hmarkov->nb_state;j++) {
              observation_reestim[i][j]->nb_value_computation();
              observation_reestim[i][j]->offset_computation();
              observation_reestim[i][j]->nb_element_computation();
              observation_reestim[i][j]->max_computation();
              if ((hmarkov->discrete_parametric_process[i]) ||
                  (hmarkov->continuous_parametric_process[i]->ident != ZERO_INFLATED_GAMMA)) {
                observation_reestim[i][j]->mean_computation();
//                observation_reestim[i][j]->variance_computation();
                observation_reestim[i][j]->variance_computation(true);
              }
            }

            if (hmarkov->discrete_parametric_process[i]) {
              for (j = 0;j < hmarkov->nb_state;j++) {
                observation_likelihood = observation_reestim[i][j]->type_parametric_estimation(hmarkov->discrete_parametric_process[i]->observation[j] ,
                                                                                               0 , true , OBSERVATION_THRESHOLD);

                if (observation_likelihood != D_INF) {
                  hmarkov->discrete_parametric_process[i]->observation[j]->computation(marginal_distribution[i]->nb_value ,
                                                                                       OBSERVATION_THRESHOLD);

                  if (hmarkov->discrete_parametric_process[i]->observation[j]->ident == BINOMIAL) {
                    for (k = hmarkov->discrete_parametric_process[i]->observation[j]->nb_value;k < marginal_distribution[i]->nb_value;k++) {
                      hmarkov->discrete_parametric_process[i]->observation[j]->mass[k] = 0.;
                    }
                  }
                }
              }
            }

            else {
              switch (hmarkov->continuous_parametric_process[i]->ident) {

              case GAMMA : {
                for (j = 0;j < hmarkov->nb_state;j++) {
                  observation_reestim[i][j]->gamma_estimation(hmarkov->continuous_parametric_process[i]->observation[j] , iter);
                }
                break;
              }

              case ZERO_INFLATED_GAMMA : {
                for (j = 0;j < hmarkov->nb_state;j++) {
                  observation_reestim[i][j]->zero_inflated_gamma_estimation(hmarkov->continuous_parametric_process[i]->observation[j] , iter);
                }
                break;
              }

              case GAUSSIAN : {
                for (j = 0;j < hmarkov->nb_state;j++) {
                  hmarkov->continuous_parametric_process[i]->observation[j]->location = observation_reestim[i][j]->mean;
                }

                switch (common_dispersion) {

                case false : {
                  for (j = 0;j < hmarkov->nb_state;j++) {
                    hmarkov->continuous_parametric_process[i]->observation[j]->dispersion = sqrt(observation_reestim[i][j]->variance);
                    if (hmarkov->continuous_parametric_process[i]->observation[j]->dispersion /
                        hmarkov->continuous_parametric_process[i]->observation[j]->location < GAUSSIAN_MIN_VARIATION_COEFF) {
                      hmarkov->continuous_parametric_process[i]->observation[j]->dispersion = hmarkov->continuous_parametric_process[i]->observation[j]->location * GAUSSIAN_MIN_VARIATION_COEFF;
                    }
                  }
                  break;
                }

                case true : {
                  variance = 0.;
                  nb_element = 0;

                  for (j = 0;j < hmarkov->nb_state;j++) {
                    for (k = observation_reestim[i][j]->offset;k < observation_reestim[i][j]->nb_value;k++) {
                      diff = k - observation_reestim[i][j]->mean;
                      variance += observation_reestim[i][j]->frequency[k] * diff * diff;
                    }

                    nb_element += observation_reestim[i][j]->nb_element;
                  }

                  variance /= nb_element;
//                  variance /= (nb_element - 1);

                  for (j = 0;j < hmarkov->nb_state;j++) {
                    hmarkov->continuous_parametric_process[i]->observation[j]->dispersion = sqrt(variance);
                  }
                  break;
                }
                }
                break;
              }

              case VON_MISES : {
                for (j = 0;j < hmarkov->nb_state;j++) {
                  observation_reestim[i][j]->mean_direction_computation(mean_direction[j]);
                  hmarkov->continuous_parametric_process[i]->observation[j]->location = mean_direction[j][3];
                }

                switch (common_dispersion) {

                case false : {
                  for (j = 0;j < hmarkov->nb_state;j++) {
                    hmarkov->continuous_parametric_process[i]->observation[j]->dispersion = von_mises_concentration_computation(mean_direction[j][2]);
                  }
                  break;
                }

                case true : {
                  global_mean_direction = 0.;
                  nb_element = 0;

                  for (j = 0;j < hmarkov->nb_state;j++) {
                    global_mean_direction += observation_reestim[i][j]->nb_element * mean_direction[j][2];
                    nb_element += observation_reestim[i][j]->nb_element;
                  }
                  concentration = von_mises_concentration_computation(global_mean_direction / nb_element);

                  for (j = 0;j < hmarkov->nb_state;j++) {
                    hmarkov->continuous_parametric_process[i]->observation[j]->dispersion = concentration;
                  }
                  break;
                }
                }
                break;
              }
              }
            }
          }

          else {
            switch (hmarkov->continuous_parametric_process[i]->ident) {
            case GAMMA :
              gamma_estimation(state_sequence_count , i ,
                               hmarkov->continuous_parametric_process[i] , iter);
              break;
            case ZERO_INFLATED_GAMMA :
              zero_inflated_gamma_estimation(state_sequence_count , i ,
                                             hmarkov->continuous_parametric_process[i] , iter);
              break;
            case GAUSSIAN :
              gaussian_estimation(state_sequence_count , i ,
                                  hmarkov->continuous_parametric_process[i]);
              break;
            case VON_MISES :
              von_mises_estimation(state_sequence_count , i ,
                                   hmarkov->continuous_parametric_process[i]);
              break;
            }
          }
        }
      }

#     ifdef MESSAGE
      os << STAT_label[STATL_ITERATION] << " " << iter << "   "
         << SEQ_label[SEQL_OBSERVED_SEQUENCES_LIKELIHOOD] << ": " << likelihood
         << "   (" << nb_state_sequence << ")" << endl;
#     endif

#     ifdef DEBUG
      if (iter % 5 == 0) {
        cout << *hmarkov;
      }
#     endif

    }
    while ((likelihood != D_INF) && (((nb_iter == I_DEFAULT) && (iter < VARIABLE_ORDER_MARKOV_NB_ITER) &&
             ((likelihood - previous_likelihood) / -likelihood > VARIABLE_ORDER_MARKOV_LIKELIHOOD_DIFF)) ||
            ((nb_iter != I_DEFAULT) && (iter < nb_iter))));

    if (likelihood != D_INF) {

#     ifdef MESSAGE
      os << "\n" << iter << " " << STAT_label[STATL_ITERATIONS] << endl;
#     endif

      // reestimation des probabilites initiales

      if (hmarkov->type == 'o') {
        reestimation(hmarkov->nb_state , chain_reestim->initial ,
                     hmarkov->initial , MIN_PROBABILITY , true);
      }

      // reestimation des probabilites de transition

      if ((hmarkov->type == 'o') && (global_initial_transition)) {
        for (i = hmarkov->nb_row - 1;i >= 1;i--) {
          if ((hmarkov->memory_type[i] != COMPLETION) && (hmarkov->order[i] > 1)) {
            for (j = 0;j < hmarkov->nb_state;j++) {
              chain_reestim->transition[hmarkov->parent[i]][j] += chain_reestim->transition[i][j];
            }
          }
        }
      }

      for (i = 1;i < hmarkov->nb_row;i++) {
        if ((hmarkov->memory_type[i] == TERMINAL) || ((hmarkov->type == 'o') &&
             (hmarkov->memory_type[i] == NON_TERMINAL))) {
          reestimation(hmarkov->nb_state , chain_reestim->transition[i] ,
                       hmarkov->transition[i] , MIN_PROBABILITY , true);
        }
        else if (hmarkov->memory_type[i] == COMPLETION) {
          for (j = 0;j < hmarkov->nb_state;j++) {
            hmarkov->transition[i][j] = hmarkov->transition[hmarkov->parent[i]][j];
          }
        }
      }

      if (hmarkov->type == 'e') {
        hmarkov->initial_probability_computation();
      }

      // reestimation des lois d'observation categorielles

      for (i = 0;i < hmarkov->nb_output_process;i++) {
        if (hmarkov->categorical_process[i]) {
          for (j = 0;j < hmarkov->nb_state;j++) {
            reestimation(marginal_distribution[i]->nb_value , observation_reestim[i][j]->frequency ,
                         hmarkov->categorical_process[i]->observation[j]->mass ,
                         MIN_PROBABILITY , true);
          }
        }

        else if (hmarkov->discrete_parametric_process[i]) {
          hmarkov->discrete_parametric_process[i]->nb_value_computation();
        }
      }
    }

    // destruction des structures de donnees de l'algorithme

    for (i = 0;i < max_length;i++) {
      delete [] forward[i];
    }
    delete [] forward;

    for (i = 0;i < max_length;i++) {
      delete [] predicted[i];
    }
    delete [] predicted;

    delete [] backward;
    delete [] cumul_backward;

    delete [] state_seq;

    delete chain_reestim;

    for (i = 0;i < hmarkov->nb_output_process;i++) {
      if (observation_reestim[i]) {
        for (j = 0;j < hmarkov->nb_state;j++) {
          delete observation_reestim[i][j];
        }
        delete [] observation_reestim[i];
      }
    }
    delete [] observation_reestim;

    if (state_sequence_count) {
      for (i = 0;i < nb_sequence;i++) {
        for (j = 0;j < length[i];j++) {
          delete [] state_sequence_count[i][j];
        }
        delete [] state_sequence_count[i];
      }
      delete [] state_sequence_count;
    }

    if (mean_direction) {
      for (i = 0;i < hmarkov->nb_state;i++) {
        delete [] mean_direction[i];
      }
      delete [] mean_direction;
    }

    delete [] pioutput;
    delete [] proutput;

    if (likelihood == D_INF) {
      delete hmarkov;
      hmarkov = NULL;
      error.update(STAT_error[STATR_ESTIMATION_FAILURE]);
    }

    else {
      if  (state_sequence) {
        hmarkov->markov_data = new VariableOrderMarkovData(*this , 'a' , (hmarkov->type == 'e' ? true : false));
        seq = hmarkov->markov_data;

        for (i = 0;i < hmarkov->nb_output_process;i++) {
          if (((hmarkov->discrete_parametric_process[i]) || (hmarkov->continuous_parametric_process[i])) &&
              (seq->characteristics[i + 1])) {
            delete seq->characteristics[i + 1];
            seq->characteristics[i + 1] = NULL;
          }
        }

        hmarkov->forward_backward(*seq);

        hmarkov->create_cumul();
        hmarkov->log_computation();
        hmarkov->viterbi(*seq);
        hmarkov->remove_cumul();

        seq->min_value_computation(0);
        seq->max_value_computation(0);
        seq->build_marginal_frequency_distribution(0);
        seq->build_characteristic(0);

        seq->build_transition_count(*hmarkov);
        seq->build_observation_frequency_distribution(hmarkov->nb_state);
        seq->build_observation_histogram(hmarkov->nb_state);

        // calcul des melanges de lois d'observation (poids deduits de la restauration)

        weight = seq->weight_computation();

        for (i = 0;i < hmarkov->nb_output_process;i++) {
          if (hmarkov->categorical_process[i]) {
            hmarkov->categorical_process[i]->restoration_weight = new Distribution(*weight);
            hmarkov->categorical_process[i]->restoration_mixture = hmarkov->categorical_process[i]->mixture_computation(hmarkov->categorical_process[i]->restoration_weight);
          }

          else if (hmarkov->discrete_parametric_process[i]) {
            for (j = 0;j < hmarkov->nb_state;j++) {
              hmarkov->discrete_parametric_process[i]->observation[j]->cumul_computation();
            }

            hmarkov->discrete_parametric_process[i]->restoration_weight = new Distribution(*weight);
            hmarkov->discrete_parametric_process[i]->restoration_mixture = hmarkov->discrete_parametric_process[i]->mixture_computation(hmarkov->discrete_parametric_process[i]->restoration_weight);
          }

          else if (hmarkov->continuous_parametric_process[i]) {
            hmarkov->continuous_parametric_process[i]->restoration_weight = new Distribution(*weight);
          }
        }

        delete weight;

#       ifdef MESSAGE
        os << "\n" << SEQ_label[SEQL_STATE_SEQUENCES_LIKELIHOOD] << ": " << seq->restoration_likelihood;

        for (i = 0;i < nb_variable;i++) {
          if (type[i] == REAL_VALUE) {
            break;
          }
        }
        if (i == nb_variable) {
          os << " | " << hmarkov->VariableOrderMarkov::likelihood_computation(*seq);
        }
        os << endl;
#       endif

      }

      else {
        hmarkov->markov_data = new VariableOrderMarkovData(*this , 'c' , (hmarkov->type == 'e' ? true : false));
        seq = hmarkov->markov_data;
        if (seq->type[0] == STATE) {
          seq->state_variable_init(INT_VALUE);
        }

        for (i = 0;i < hmarkov->nb_output_process;i++) {
          if (((hmarkov->discrete_parametric_process[i]) || (hmarkov->continuous_parametric_process[i])) &&
              (seq->characteristics[i])) {
            delete seq->characteristics[i];
            seq->characteristics[i] = NULL;
          }
        }
      }

      for (i = 0;i < hmarkov->nb_output_process;i++) {
        if (hmarkov->categorical_process[i]) {
          for (j = 0;j < hmarkov->nb_state;j++) {
            hmarkov->categorical_process[i]->observation[j]->cumul_computation();

            hmarkov->categorical_process[i]->observation[j]->max_computation();
//            hmarkov->categorical_process[i]->observation[j]->mean_computation();
//            hmarkov->categorical_process[i]->observation[j]->variance_computation();
          }
        }
      }

      // calcul de la vraisemblance et des lois caracteristiques du modele

      seq->likelihood = hmarkov->likelihood_computation(*this , seq->posterior_probability);

#     ifdef DEBUG
//      cout << *hmarkov;
      cout << "iteration " << iter << "  "
           << SEQ_label[SEQL_OBSERVED_SEQUENCES_LIKELIHOOD] << ": " << seq->likelihood << endl;
#     endif

#     ifdef MESSAGE
      if  ((state_sequence) && (seq->nb_sequence <= POSTERIOR_PROBABILITY_NB_SEQUENCE)) {
        os << "\n" << SEQ_label[SEQL_POSTERIOR_STATE_SEQUENCE_PROBABILITY] << endl;
        for (i = 0;i < seq->nb_sequence;i++) {
          os << SEQ_label[SEQL_SEQUENCE] << " " << seq->identifier[i] << ": "
             << seq->posterior_probability[i] << endl;
        }
      }
#     endif

      hmarkov->component_computation();
      hmarkov->characteristic_computation(*seq , counting_flag , I_DEFAULT , false);

      // calcul des melanges de lois d'observation (poids theoriques)

      switch (hmarkov->type) {

      case 'o' : {
        weight = hmarkov->state_process->weight_computation();
        break;
      }

      case 'e' : {
        weight = new Distribution(hmarkov->nb_state);

        for (i = 0;i < hmarkov->nb_state;i++) {
          weight->mass[i] = 0.;
        }
        for (i = 1;i < hmarkov->nb_row;i++) {
          if ((hmarkov->memory_type[i] == TERMINAL) || (hmarkov->memory_type[i] == COMPLETION)) {
            weight->mass[hmarkov->state[i][0]] += hmarkov->initial[i];
          }
        }

        weight->cumul_computation();
        weight->max_computation();
        break;
      }
      }

      for (i = 0;i < hmarkov->nb_output_process;i++) {
        if (hmarkov->categorical_process[i]) {
          hmarkov->categorical_process[i]->weight = new Distribution(*weight);
          hmarkov->categorical_process[i]->mixture = hmarkov->categorical_process[i]->mixture_computation(hmarkov->categorical_process[i]->weight);
        }

        else if (hmarkov->discrete_parametric_process[i]) {
          hmarkov->discrete_parametric_process[i]->weight = new Distribution(*weight);
          hmarkov->discrete_parametric_process[i]->mixture = hmarkov->discrete_parametric_process[i]->mixture_computation(hmarkov->discrete_parametric_process[i]->weight);
        }

        else if (hmarkov->continuous_parametric_process[i]) {
          hmarkov->continuous_parametric_process[i]->weight = new Distribution(*weight);
        }
      }

      delete weight;
    }
  }

  return hmarkov;
}


};  // namespace sequence_analysis
