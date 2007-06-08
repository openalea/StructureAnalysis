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



#include <math.h>
#include <sstream>
#include "stat_tool/stat_tools.h"
#include "stat_tool/distance_matrix.h"
#include "stat_tool/regression.h"
#include "stat_tool/curves.h"
#include "stat_tool/markovian.h"
#include "stat_tool/stat_label.h"
#include "sequences.h"
#include "markov.h"
#include "hidden_markov.h"
#include "sequence_label.h"

#include "stat_tool/distribution_reestimation.h"   // probleme compilateur C++ Windows

using namespace std;


extern char* label(const char *file_name);



/*--------------------------------------------------------------*
 *
 *  Initialisation des parametres d'une chaine de Markov
 *
 *  arguments : flag sur la nature de la chaine de Markov,
 *              probabilite de rester dans un etat.
 *
 *--------------------------------------------------------------*/

void Hidden_markov::init(bool left_right , double self_transition)

{
  register int i , j;
  int power = nb_row / nb_state , state , state_index[ORDER];


  accessibility = new bool*[nb_state];
  for (i = 0;i < nb_state;i++) {
    accessibility[i] = new bool[nb_state];
  }

  state_type = new char[nb_state];

  switch (left_right) {

  // cas chaine de Markov ou toutes les transitions sont possibles

  case false : {
    nb_component = 1;
    component_nb_state = new int[nb_component];
    component_nb_state[0] = nb_state;
    component = new int*[nb_component];
    component[0] = new int[component_nb_state[0]];

    for (i = 0;i < nb_state;i++) {
      for (j = 0;j < nb_state;j++) {
        accessibility[i][j] = true;
      }

      component[0][i] = i;
      state_type[i] = 'r';
    }

    for (i = 0;i < nb_state;i++) {
      initial[i] = 1. / (double)nb_state;
    }

    for (i = 0;i < nb_row;i++) {
      state = i / power;
      for (j = 0;j < state;j++) {
        transition[i][j] = (1. - self_transition) / (nb_state - 1);
      }
      transition[i][state] = self_transition;
      for (j = state + 1;j < nb_state;j++) {
        transition[i][j] = (1. - self_transition) / (nb_state - 1);
      }
    }

    break;
  }

  // cas chaine de Markov gauche-droite

  case true : {
    nb_component = nb_state;
    component_nb_state = new int[nb_component];
    component = new int*[nb_component];

    for (i = 0;i < nb_state;i++) {
      for (j = 0;j <= i;j++) {
        accessibility[i][j] = false;
      }
      for (j = i + 1;j < nb_state;j++) {
        accessibility[i][j] = true;
      }

      component_nb_state[i] = 1;
      component[i] = new int[component_nb_state[i]];
      component[i][0] = i;

      if (i < nb_state - 1) {
        state_type[i] = 't';
      }
      else {
        state_type[i] = 'a';
      }
    }

    for (i = 0;i < nb_state - 1;i++) {
      initial[i] = 1. / (double)(nb_state - 1);
    }
    initial[nb_state - 1] = 0.;

    for (i = 0;i < order;i++) {
      state_index[i] = 0;
    }

    for (i = 0;i < nb_row;i++) {
      for (j = 0;j < state_index[order - 1];j++) {
        transition[i][j] = 0.;
      }

      for (j = 1;j < order;j++) {
        if (state_index[j] < state_index[j - 1]) {
          break;
        }
      }

      if ((j == order) && (state_index[order - 1] < nb_state - 1))  {
        transition[i][state_index[order - 1]] = self_transition;
        for (j = state_index[order - 1] + 1;j < nb_state;j++) {
          transition[i][j] = (1. - self_transition) / (nb_state - (state_index[order - 1] + 1));
        }
      }

      else {
        transition[i][state_index[order - 1]] = 1.;
        for (j = state_index[order - 1] + 1;j < nb_state;j++) {
           transition[i][j] = 0.;
        }
      }

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

    break;
  }
  }
}


/*--------------------------------------------------------------*
 *
 *  Initialisation des parametres d'une chaine de Markov en
 *  respectant les probabilites a 0.
 *
 *  argument : probabilite de rester dans un etat.
 *
 *--------------------------------------------------------------*/

void Hidden_markov::init(double self_transition)

{
  register int i , j;
  int nb_positive , power = nb_row / nb_state , state;
  double self;


  nb_positive = 0;
  for (i = 0;i < nb_state;i++) {
    if (initial[i] > 0.) {
      nb_positive++;
    }
  }

  for (i = 0;i < nb_state;i++) {
    if (initial[i] > 0.) {
      initial[i] = 1. / (double)nb_positive;
    }
  }

  for (i = 0;i < nb_row;i++) {
    state = i / power;
    if (transition[i][state] > 0.) {
      self = self_transition;
    }
    else {
      self = 0.;
    }

    nb_positive = 0;

    for (j = 0;j < state;j++) {
      if (transition[i][j] > 0.) {
        nb_positive++;
      }
    }
    for (j = state + 1;j < nb_state;j++) {
      if (transition[i][j] > 0.) {
        nb_positive++;
      }
    }

    for (j = 0;j < state;j++) {
      if (transition[i][j] > 0.) {
        transition[i][j] = (1. - self) / nb_positive;
      }
    }
    if (transition[i][j] < 1.) {
      transition[i][j] = self;
    }
    for (j = state + 1;j < nb_state;j++) {
      if (transition[i][j] > 0.) {
        transition[i][j] = (1. - self) / nb_positive;
      }
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la vraisemblance de sequences pour une chaine de Markov cachee
 *  par l'algorithme forward.
 *
 *  arguments : reference sur un objet Markovian_sequences, indice de la sequence.
 *
 *--------------------------------------------------------------*/

double Hidden_markov::likelihood_computation(const Markovian_sequences &seq , int index) const

{
  register int i , j , k , m;
  int **poutput , power[ORDER] , state_index[ORDER];
  double likelihood = 0. , **ptransition , *forward , *pforward , *auxiliary , norm;


  // verification de la compatibilite entre le modele et les donnees

  if (nb_output_process == seq.nb_variable) {
    for (i = 0;i < nb_output_process;i++) {
      if (process[i + 1]->nb_value < seq.marginal[i]->nb_value) {
        likelihood = D_INF;
        break;
      }
    }
  }

  else {
    likelihood = D_INF;
  }

  if (likelihood != D_INF) {

    // initialisations

    i = 1;
    for (j = 0;j < order;j++) {
      power[j] = i;
      i *= nb_state;
    }

    forward = new double[nb_row];
    auxiliary = new double[nb_row];

    poutput = new int*[seq.nb_variable];

    for (i = 0;i < seq.nb_sequence;i++) {
      if ((index == I_DEFAULT) || (index == i)) {
        for (j = 0;j < seq.nb_variable;j++) {
          poutput[j] = seq.sequence[i][j];
        }

        norm = 0.;
        j = 0;

        for (k = 0;k < nb_row;k++) {
          if (k == self_row[j]) {
            forward[k] = initial[j];
            for (m = 0;m < nb_output_process;m++) {
              forward[k] *= process[m + 1]->observation[j]->mass[*poutput[m]];
            }
            norm += forward[k];
            j++;
          }

          else {
            forward[k] = 0.;
          }
        }

        if (norm > 0.) {
          for (j = 0;j < nb_row;j++) {
            forward[j] /= norm;
          }
          likelihood += log(norm);
        }

        else {
          likelihood = D_INF;
          break;
        }

        for (j = 1;j < seq.length[i];j++) {
          for (k = 0;k < seq.nb_variable;k++) {
            poutput[k]++;
          }

          for (k = 0;k < order;k++) {
            state_index[k] = 0;
          }
          norm = 0.;

          for (k = 0;k < nb_row;k++) {
            ptransition = transition;
            pforward = forward;
            for (m = 0;m < order - 1;m++) {
              ptransition += state_index[m] * power[m + 1];
              pforward += state_index[m] * power[m + 1];
            }

            auxiliary[k] = 0.;
            for (m = 0;m < nb_state;m++) {
              auxiliary[k] += *(*ptransition + state_index[order - 1]) * *pforward++;
              ptransition++;
            }

            for (m = 0;m < nb_output_process;m++) {
              auxiliary[k] *= process[m + 1]->observation[state_index[order - 1]]->mass[*poutput[m]];
            }
            norm += auxiliary[k];

            // mise a jour des indices des etats

            for (m = 0;m < order;m++) {
              if (state_index[m] < nb_state - 1) {
                state_index[m]++;
                break;
              }
              else {
                state_index[m] = 0;
              }
            }
          }

          if (norm > 0.) {
            for (k = 0;k < nb_row;k++) {
              forward[k] = auxiliary[k] / norm;
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
      }
    }

    delete [] forward;
    delete [] auxiliary;

    delete [] poutput;
  }

  return likelihood;
}


/*--------------------------------------------------------------*
 *
 *  Estimation des parametres d'une chaine de Markov cachee a partir
 *  d'un echantillon de sequences par l'algorithme EM.
 *
 *  arguments : reference sur un objet Format_error, stream, chaine de Markov cachee initiale,
 *              flags sur le calcul des lois de comptage et sur le calcul des sequences
 *              d'etats optimales, nombre d'iterations, flag sur le calcul des lois caracteristiques.
 *
 *--------------------------------------------------------------*/

Hidden_markov* Markovian_sequences::hidden_markov_estimation(Format_error &error , ostream &os ,
                                                             const Hidden_markov &ihmarkov ,
                                                             bool counting_flag , bool state_sequence ,
                                                             int nb_iter , bool characteristic_flag) const

{
  bool status;
  register int i , j , k , m;
  int iter , **poutput , power[ORDER] , state_index[ORDER];
  double likelihood = D_INF , previous_likelihood , buff , **ptransition , **forward ,
         *pforward , norm , **predicted , *backward , *auxiliary , *pauxiliary ,
         **ptransition_reestim , ***observation_reestim , *reestim;
  Chain_reestimation<double> *chain_reestim;
  Hidden_markov *hmarkov;
  Markov_data *seq;


  hmarkov = 0;
  error.init();

  // test nombre de valeurs observees par variable

  status = false;
  for (i = 0;i < nb_variable;i++) {
    if (marginal[i]->nb_value > 1) {
      status = true;
      break;
    }
  }

  if (!status) {
    error.update(SEQ_error[SEQR_VARIABLE_NB_VALUE]);
  }

  for (i = 0;i < nb_variable;i++) {
    if (!characteristics[i]) {
      for (j = 0;j < marginal[i]->nb_value;j++) {
        if (marginal[i]->frequency[j] == 0) {
          status = false;
          ostringstream error_message;
          error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                        << STAT_error[STATR_MISSING_VALUE] << " " << j;
          error.update((error_message.str()).c_str());
        }
      }
    }
  }

  if (ihmarkov.nb_output_process != nb_variable) {
    status = false;
    error.update(SEQ_error[SEQR_NB_OUTPUT_PROCESS]);
  }

  else {
    for (i = 0;i < nb_variable;i++) {
      if (ihmarkov.process[i + 1]->nb_value != marginal[i]->nb_value) {
        status = false;
        ostringstream error_message;
        error_message << STAT_label[STATL_OUTPUT_PROCESS] << " " << i + 1 << ": "
                      << SEQ_error[SEQR_NB_OUTPUT];
        error.update((error_message.str()).c_str());
      }
    }
  }

  if ((nb_iter != I_DEFAULT) && (nb_iter < 1)) {
    status = false;
    error.update(STAT_error[STATR_NB_ITERATION]);
  }

  if (status) {

    // creation de la chaine de Markov cachee

    hmarkov = new Hidden_markov(ihmarkov , false , false);

#   ifdef DEBUG
    cout << *hmarkov;
#   endif

    // initialisations

    i = 1;
    for (j = 0;j < hmarkov->order;j++) {
      power[j] = i;
      i *= hmarkov->nb_state;
    }

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

    chain_reestim = new Chain_reestimation<double>(hmarkov->type , hmarkov->nb_state ,
                                                   (int)pow((double)hmarkov->nb_state , hmarkov->order));

    observation_reestim = new double**[hmarkov->nb_output_process];
    for (i = 0;i < hmarkov->nb_output_process;i++) {
      observation_reestim[i] = new double*[hmarkov->nb_state];
      for (j = 0;j < hmarkov->nb_state;j++) {
        observation_reestim[i][j] = new double[marginal[i]->nb_value];
      }
    }

    poutput = new int*[nb_variable];

    iter = 0;
    do {
      iter++;
      previous_likelihood = likelihood;
      likelihood = 0.;

      // initialisation des quantites de reestimation

      chain_reestim->init();

      for (i = 0;i < hmarkov->nb_output_process;i++) {
        for (j = 0;j < hmarkov->nb_state;j++) {
          reestim = observation_reestim[i][j];
          for (k = 0;k < marginal[i]->nb_value;k++) {
            *reestim++ = 0.;
          }
        }
      }

      for (i = 0;i < nb_sequence;i++) {
        for (j = 0;j < nb_variable;j++) {
          poutput[j] = sequence[i][j];
        }

        // recurrence "forward"

        norm = 0.;
        j = 0;

        for (k = 0;k < hmarkov->nb_row;k++) {
          if (k == hmarkov->self_row[j]) {
            forward[0][k] = hmarkov->initial[j];
            for (m = 0;m < hmarkov->nb_output_process;m++) {
              forward[0][k] *= hmarkov->process[m + 1]->observation[j]->mass[*poutput[m]];
            }
            norm += forward[0][k];
            j++;
          }

          else {
            forward[0][k] = 0.;
          }
        }

        if (norm > 0.) {
          for (j = 0;j < hmarkov->nb_row;j++) {
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
            poutput[k]++;
          }

          for (k = 0;k < hmarkov->order;k++) {
            state_index[k] = 0;
          }
          norm = 0.;

          for (k = 0;k < hmarkov->nb_row;k++) {
            ptransition = hmarkov->transition;
            pforward = forward[j - 1];
            for (m = 0;m < hmarkov->order - 1;m++) {
              ptransition += state_index[m] * power[m + 1];
              pforward += state_index[m] * power[m + 1];
            }

            forward[j][k] = 0.;
            for (m = 0;m < hmarkov->nb_state;m++) {
              forward[j][k] += *(*ptransition + state_index[hmarkov->order - 1]) * *pforward++;
              ptransition++;
            }
            predicted[j][k] = forward[j][k];

            for (m = 0;m < hmarkov->nb_output_process;m++) {
              forward[j][k] *= hmarkov->process[m + 1]->observation[state_index[hmarkov->order - 1]]->mass[*poutput[m]];
            }
            norm += forward[j][k];

            // mise a jour des indices des etats

            for (m = 0;m < hmarkov->order;m++) {
              if (state_index[m] < hmarkov->nb_state - 1) {
                state_index[m]++;
                break;
              }
              else {
                state_index[m] = 0;
              }
            }
          }

          if (norm > 0.) {
            for (k = 0;k < hmarkov->nb_row;k++) {
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
        for (k = 0;k < hmarkov->nb_row;k++) {
          backward[k] = forward[j][k];

          // accumulation des quantites de reestimation des probabilites d'observation

          for (m = 0;m < hmarkov->nb_output_process;m++) {
            observation_reestim[m][k / power[hmarkov->order - 1]][*poutput[m]] += backward[k];
          }
        }

        for (j = length[i] - 2;j >= 0;j--) {
          for (k = 0;k < hmarkov->nb_row;k++) {
            if (predicted[j + 1][k] > 0.) {
              auxiliary[k] = backward[k] / predicted[j + 1][k];
            }
            else {
              auxiliary[k] = 0.;
            }
          }

          for (k = 0;k < nb_variable;k++) {
            poutput[k]--;
          }

          for (k = 0;k < hmarkov->nb_row;k++) {
            ptransition = hmarkov->transition;
            ptransition_reestim = chain_reestim->transition;
            pauxiliary = auxiliary;

            ptransition += state_index[0] * power[0];
            ptransition_reestim += state_index[0] * power[0];
            for (m = 1;m < hmarkov->order;m++) {
              ptransition += state_index[m] * power[m];
              ptransition_reestim += state_index[m] * power[m];
              pauxiliary += state_index[m] * power[m - 1];
            }

            backward[k] = 0.;
            for (m = 0;m < hmarkov->nb_state;m++) {
              buff = *pauxiliary * *(*ptransition + m) * forward[j][k];
              pauxiliary += power[hmarkov->order - 1];
              backward[k] += buff;

              // accumulation des quantites de reestimation des probabilites de transition

              *(*ptransition_reestim + m) += buff;
            }

            // accumulation des quantites de reestimation des probabilites initiales et
            // des probabilites d'observation

            if ((j == 0) && (hmarkov->self_row[state_index[hmarkov->order - 1]] == k)) {
              chain_reestim->initial[state_index[hmarkov->order - 1]] += backward[k];
            }

            for (m = 0;m < hmarkov->nb_output_process;m++) {
              observation_reestim[m][state_index[hmarkov->order - 1]][*poutput[m]] += backward[k];
            }

            // mise a jour des indices des etats

            for (m = 0;m < hmarkov->order;m++) {
              if (state_index[m] < hmarkov->nb_state - 1) {
                state_index[m]++;
                break;
              }
              else {
                state_index[m] = 0;
              }
            }
          }

#         ifdef DEBUG
/*          cout << j << " : ";
          sum = 0.;
          for (k = 0;k < hmarkov->nb_row;k++) {
            sum += backward[k];
            cout << backward[k] << " ";
          }
          cout << "| " << sum << endl; */
#         endif

        }
      }

      if (likelihood != D_INF) {

        // reestimation des probabilites initiales

        reestimation(hmarkov->nb_state , chain_reestim->initial ,
                     hmarkov->initial , MIN_PROBABILITY , false);

        // reestimation des probabilites de transition

        for (i = 0;i < hmarkov->nb_row;i++) {
          reestimation(hmarkov->nb_state , chain_reestim->transition[i] ,
                       hmarkov->transition[i] , MIN_PROBABILITY , false);
        }

        // reestimation des probabilites d'observation

        for (i = 0;i < hmarkov->nb_output_process;i++) {
          for (j = 0;j < hmarkov->nb_state;j++) {
            reestimation(marginal[i]->nb_value , observation_reestim[i][j] ,
                         hmarkov->process[i + 1]->observation[j]->mass , MIN_PROBABILITY , false);
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
    while ((likelihood != D_INF) && (((nb_iter == I_DEFAULT) && (iter < MARKOV_NB_ITER) &&
             ((likelihood - previous_likelihood) / -likelihood > MARKOV_LIKELIHOOD_DIFF)) ||
            ((nb_iter != I_DEFAULT) && (iter < nb_iter))));

    if (likelihood != D_INF) {

#     ifdef MESSAGE
      os << "\n" << iter << " " << STAT_label[STATL_ITERATIONS] << endl;
#     endif

      // reestimation des probabilites initiales

      reestimation(hmarkov->nb_state , chain_reestim->initial ,
                   hmarkov->initial , MIN_PROBABILITY , true);

      // reestimation des probabilites de transition

      for (i = 0;i < hmarkov->nb_row;i++) {
        reestimation(hmarkov->nb_state , chain_reestim->transition[i] ,
                     hmarkov->transition[i] , MIN_PROBABILITY , true);
      }

      // reestimation des probabilites d'observation

      for (i = 0;i < hmarkov->nb_output_process;i++) {
        for (j = 0;j < hmarkov->nb_state;j++) {
          reestimation(marginal[i]->nb_value , observation_reestim[i][j] ,
                       hmarkov->process[i + 1]->observation[j]->mass , MIN_PROBABILITY , true);
        }
      }
    }

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
      for (j = 0;j < hmarkov->nb_state;j++) {
        delete [] observation_reestim[i][j];
      }
      delete [] observation_reestim[i];
    }
    delete [] observation_reestim;

    delete [] poutput;

    if (likelihood == D_INF) {
      delete hmarkov;
      hmarkov = 0;
      error.update(STAT_error[STATR_ESTIMATION_FAILURE]);
    }

    else {
      if (state_sequence) {
        hmarkov->markov_data = new Markov_data(*this , 0);
        seq = hmarkov->markov_data;
        seq->type[0] = STATE;

        hmarkov->create_cumul();
        hmarkov->log_computation();

        seq->likelihood = hmarkov->viterbi(*seq);

        hmarkov->remove_cumul();

        seq->max_value[0] = hmarkov->nb_state - 1;
        seq->build_marginal_histogram(0);
        seq->build_characteristic(0);

        seq->build_transition_count(hmarkov->order);
        seq->build_observation_histogram();

#       ifdef DEBUG
        cout << "\n" << SEQ_label[SEQL_STATE_SEQUENCES_LIKELIHOOD] << ": " << seq->likelihood
             << " | " << hmarkov->Markov::likelihood_computation(*seq) << endl;
#       endif

      }

      else {
        hmarkov->markov_data = new Markov_data(*this);
        seq = hmarkov->markov_data;
        seq->state_variable_init(INT_VALUE);
      }

      for (i = 1;i <= hmarkov->nb_output_process;i++) {
        for (j = 0;j < hmarkov->nb_state;j++) {
          hmarkov->process[i]->observation[j]->cumul_computation();

          hmarkov->process[i]->observation[j]->max_computation();
//          hmarkov->process[i]->observation[j]->mean_computation();
//          hmarkov->process[i]->observation[j]->variance_computation();
        }
      }

      // calcul de la vraisemblance et des lois caracteristiques du modele

      seq->hidden_likelihood = hmarkov->likelihood_computation(*this);

#     ifdef DEBUG
//      cout << *hmarkov;
      cout << "iteration " << iter << "  "
           << SEQ_label[SEQL_OBSERVED_SEQUENCES_LIKELIHOOD] << ": " << seq->hidden_likelihood << endl;
#     endif

      if (characteristic_flag) {
        hmarkov->component_computation();
        hmarkov->characteristic_computation(*seq , counting_flag , I_DEFAULT , false);
      }
    }
  }

  return hmarkov;
}


/*--------------------------------------------------------------*
 *
 *  Estimation des parametres d'une chaine de Markov cachee a partir
 *  d'un echantillon de sequences.
 *
 *  arguments : reference sur un objet Format_error, stream, nombre d'etats,
 *              flag sur la nature de la chaine de Markov, ordre, flags sur le calcul
 *              des lois de comptage et sur le calcul des sequences d'etats optimales,
 *              probabilite de rester dans un etat initiale, nombre d'iterations.
 *
 *--------------------------------------------------------------*/

Hidden_markov* Markovian_sequences::hidden_markov_estimation(Format_error &error , ostream &os ,
                                                             int nb_state , bool left_right , int order ,
                                                             bool counting_flag , bool state_sequence ,
                                                             double self_transition , int nb_iter) const

{
  bool status = true;
  register int i;
  int nb_value[SEQUENCE_NB_VARIABLE];
  Hidden_markov *ihmarkov , *hmarkov;


  hmarkov = 0;
  error.init();

  if ((nb_state < 2) || (nb_state > NB_STATE)) {
    status = false;
    error.update(SEQ_error[SEQR_NB_STATE]);
  }
  if ((order < 1) || (order > ORDER)) {
    status = false;
    error.update(SEQ_error[SEQR_ORDER]);
  }
  if ((self_transition != D_DEFAULT) && ((self_transition <= 0.) || (self_transition >= 1.))) {
    status = false;
    error.update(SEQ_error[SEQR_SELF_TRANSITION]);
  }

  if (status) {
    for (i = 0;i < nb_variable;i++) {
      nb_value[i] = marginal[i]->nb_value;
    }

    ihmarkov = new Hidden_markov(nb_state , order , nb_variable , nb_value);

    // initialisation des parametres de la chaine de Markov

    if (self_transition == D_DEFAULT) {
      self_transition = MAX(1. - 1. / hlength->mean , SELF_TRANSITION);
    }
    ihmarkov->init(left_right , self_transition);

    // initialisation des lois d'observation

    for (i = 0;i < ihmarkov->nb_output_process;i++) {
      ihmarkov->process[i + 1]->init();
    }

    hmarkov = hidden_markov_estimation(error , os , *ihmarkov , counting_flag ,
                                       state_sequence , nb_iter);

    delete ihmarkov;
  }

  return hmarkov;
}


/*--------------------------------------------------------------*
 *
 *  Estimation du nombre d'etats et des parametres d'une chaine
 *  de Markov cachee a partir d'un echantillon de sequences.
 *
 *  arguments : reference sur un objet Format_error, stream, chaine de Markov cachee initiale,
 *              etat servant de reference pour ajouter des etats, nombre maximum d'etats,
 *              type de penalisation (AIC(c)/BIC), flags sur le calcul des lois de comptage et
 *              sur le calcul des sequences d'etats optimales, probabilite de rester dans un etat initiale.
 *
 *--------------------------------------------------------------*/

Hidden_markov* Markovian_sequences::hidden_markov_nb_state_estimation(Format_error &error , ostream &os ,
                                                                      const Hidden_markov &ihmarkov ,
                                                                      int state , int max_nb_state ,
                                                                      int penalty_type , bool counting_flag ,
                                                                      bool state_sequence ,
                                                                      double self_transition) const

{
  bool status = true;
  register int i , j;
  int nb_parameter[NB_STATE + 1];
  double penalty , max_likelihood , likelihood[NB_STATE + 1] , penalized_likelihood[NB_STATE + 1];
  Hidden_markov *hmarkov , *initial_hmarkov , *phmarkov;


  hmarkov = 0;
  error.init();

  if ((state < 0) || (state > ihmarkov.nb_state - 1)) {
    status = false;
    error.update(SEQ_error[SEQR_STATE_INDEX]);
  }
  else {
    if (ihmarkov.state_type[state] != 'r') {
      status = false;
      error.update(SEQ_error[SEQR_STATE_TYPE]);
    }
  }
  if ((max_nb_state <= ihmarkov.nb_state) || (max_nb_state > NB_STATE)) {
    status = false;
    error.update(SEQ_error[SEQR_MAX_NB_STATE]);
  }
  if ((self_transition != D_DEFAULT) && ((self_transition <= 0.) || (self_transition >= 1.))) {
    status = false;
    error.update(SEQ_error[SEQR_SELF_TRANSITION]);
  }

  if (status) {

    // estimation des parametres d'une chaine de Markov cachee

    hmarkov = hidden_markov_estimation(error , os , ihmarkov , false , state_sequence ,
                                       I_DEFAULT , false);

    if (hmarkov) {

      // calcul du terme de compensation

      switch (penalty_type) {
      case AIC :
        penalty = 1.;
        break;
      case BIC :
        penalty = 0.5 * log((double)cumul_length);
        break;
      }

      nb_parameter[ihmarkov.nb_state] = hmarkov->nb_parameter_computation(MIN_PROBABILITY);
      likelihood[ihmarkov.nb_state] = hmarkov->markov_data->hidden_likelihood;

      if (penalty_type == AICc) {
        if (nb_parameter[ihmarkov.nb_state] < cumul_length - 1) {
          penalized_likelihood[ihmarkov.nb_state] = likelihood[ihmarkov.nb_state] - (double)(nb_parameter[ihmarkov.nb_state] * cumul_length)/
                                                    (double)(cumul_length - nb_parameter[ihmarkov.nb_state] - 1);
        }
        else {
          penalized_likelihood[ihmarkov.nb_state] = D_INF;
        }
      }

      else {
        penalized_likelihood[ihmarkov.nb_state] = likelihood[ihmarkov.nb_state] - nb_parameter[ihmarkov.nb_state] * penalty;
      }

      max_likelihood = penalized_likelihood[ihmarkov.nb_state];

      initial_hmarkov = new Hidden_markov(*hmarkov , state);

      if (self_transition == D_DEFAULT) {
        self_transition = MAX(1. - 1. / hlength->mean , SELF_TRANSITION);
      }

      for (i = ihmarkov.nb_state + 1;i <= max_nb_state;i++) {
        initial_hmarkov->init(self_transition);
        for (j = 0;j < initial_hmarkov->nb_output_process;j++) {
          initial_hmarkov->process[j + 1]->init();
        }

        phmarkov = hidden_markov_estimation(error , os , *initial_hmarkov , false , state_sequence ,
                                            I_DEFAULT , false);
        delete initial_hmarkov;

        if (phmarkov) {
          nb_parameter[i] = phmarkov->nb_parameter_computation(MIN_PROBABILITY);
          likelihood[i] = phmarkov->markov_data->hidden_likelihood;

          if (penalty_type == AICc) {
            if (nb_parameter[i] < cumul_length - 1) {
              penalized_likelihood[i] = likelihood[i] - (double)(nb_parameter[i] * cumul_length) /
                                        (double)(cumul_length - nb_parameter[i] - 1);
            }
            else {
              penalized_likelihood[i] = D_INF;
            }
          }

          else {
            penalized_likelihood[i] = likelihood[i] - nb_parameter[i] * penalty;
          }

          if (i < max_nb_state) {
            initial_hmarkov = new Hidden_markov(*phmarkov , state);
          }

          if (penalized_likelihood[i] > max_likelihood) {
            max_likelihood = penalized_likelihood[i];
            delete hmarkov;
            hmarkov = phmarkov;
          }
          else {
            delete phmarkov;
          }
        }

        else {
          likelihood[i] = D_INF;
          break;
        }
      }

#     ifdef MESSAGE
      {
        double norm = 0. , weight[NB_STATE + 1];

        for (i = ihmarkov.nb_state;i <= max_nb_state;i++) {
          if (likelihood[i] != D_INF) {
            weight[i] = exp(penalized_likelihood[i] - max_likelihood);
            norm += weight[i];
          }
          else {
            break;
          }
        }

        for (i = ihmarkov.nb_state;i <= max_nb_state;i++) {
          if (likelihood[i] != D_INF) {
            os << "\n" << i << " " << STAT_label[STATL_STATES]
               << "   2 * " << STAT_label[STATL_LIKELIHOOD] << ": " << 2 * likelihood[i] << "   "
               << nb_parameter[i] << " " << STAT_label[nb_parameter[i] == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS]
               << "   2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("
               << STAT_criterion_word[penalty_type] << "): " << 2 * penalized_likelihood[i] << "   "
               << STAT_label[STATL_WEIGHT] << ": " << weight[i] / norm << endl;
          }
          else {
            break;
          }
        }
      }
#     endif

      // calcul des lois caracteristiques du modele

      hmarkov->component_computation();
      hmarkov->characteristic_computation(*(hmarkov->markov_data) , counting_flag ,
                                          I_DEFAULT , false);
    }
  }

  return hmarkov;
}


/*--------------------------------------------------------------*
 *
 *  Estimation du nombre d'etats et des parametres d'une chaine
 *  de Markov cachee a partir d'un echantillon de sequences.
 *
 *  arguments : reference sur un objet Format_error, stream,
 *              nombres minimum et maximum d'etats, type de penalisation (AIC(c)/BIC), ordre,
 *              flags sur le calcul des lois de comptage et sur le calcul des sequences
 *              d'etats optimales, probabilite de rester dans un etat initiale.
 *
 *--------------------------------------------------------------*/

Hidden_markov* Markovian_sequences::hidden_markov_nb_state_estimation(Format_error &error , ostream &os ,
                                                                      int min_nb_state , int max_nb_state ,
                                                                      int penalty_type , int order ,
                                                                      bool counting_flag , bool state_sequence ,
                                                                      double self_transition) const

{
  bool status = true;
  register int i;
  int nb_value[SEQUENCE_NB_VARIABLE];
  Hidden_markov *ihmarkov , *hmarkov;


  hmarkov = 0;
  error.init();

  if ((min_nb_state < 2) || (min_nb_state >= max_nb_state)) {
    status = false;
    error.update(SEQ_error[SEQR_MIN_NB_STATE]);
  }
  if ((max_nb_state <= min_nb_state) || (max_nb_state > NB_STATE)) {
    status = false;
    error.update(SEQ_error[SEQR_MAX_NB_STATE]);
  }
  if ((order < 1) || (order > ORDER)) {
    status = false;
    error.update(SEQ_error[SEQR_ORDER]);
  }
  if ((self_transition != D_DEFAULT) && ((self_transition <= 0.) || (self_transition >= 1.))) {
    status = false;
    error.update(SEQ_error[SEQR_SELF_TRANSITION]);
  }

  if (status) {
    for (i = 0;i < nb_variable;i++) {
      nb_value[i] = marginal[i]->nb_value;
    }

    ihmarkov = new Hidden_markov(min_nb_state , order , nb_variable , nb_value);

    // initialisation des parametres de la chaine de Markov

    if (self_transition == D_DEFAULT) {
      self_transition = MAX(1. - 1. / hlength->mean , SELF_TRANSITION);
    }
    ihmarkov->init(false , self_transition);

    // initialisation des lois d'observation

    for (i = 0;i < ihmarkov->nb_output_process;i++) {
      ihmarkov->process[i + 1]->init();
    }

    hmarkov = hidden_markov_nb_state_estimation(error , os , *ihmarkov , 0 , max_nb_state ,
                                                penalty_type , counting_flag ,
                                                state_sequence , self_transition);

    delete ihmarkov;
  }

  return hmarkov;
}


/*--------------------------------------------------------------*
 *
 *  Calcul des probabilites a posteriori des etats par l'algorithme forward-backward.
 *
 *  arguments : reference sur un objet Semi_markov_data, indice de la sequence,
 *              stream, format de fichier ('a' : ASCII, 's' : Spreadsheet, 'g' : Gnuplot).
 *
 *--------------------------------------------------------------*/

double Hidden_markov::forward_backward(const Markov_data &seq , int index ,
                                       ostream &os , char format) const

{
  register int i , j , k;
  int *pstate , **poutput , power[ORDER] , state_index[ORDER];
  double seq_likelihood , state_seq_likelihood , **ptransition , **forward , *pforward , norm ,
         **predicted , **backward , *auxiliary , *pauxiliary , backward_max , **state_backward;


  // initialisations

  i = 1;
  for (j = 0;j < order;j++) {
    power[j] = i;
    i *= nb_state;
  }

  forward = new double*[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    forward[i] = new double[nb_row];
  }

  predicted = new double*[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    predicted[i] = new double[nb_row];
  }

  backward = new double*[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    backward[i] = new double[nb_row];
  }

  auxiliary = new double[nb_row];

  state_backward = new double*[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    state_backward[i] = new double[nb_state];
  }

  poutput = new int*[nb_output_process];

  for (i = 0;i < nb_output_process;i++) {
    poutput[i] = seq.sequence[index][i + 1];
  }

  // recurrence "forward"

  seq_likelihood = 0.;
  norm = 0.;
  i = 0;

  for (j = 0;j < nb_row;j++) {
    if (j == self_row[i]) {
      forward[0][j] = initial[i];
      for (k = 0;k < nb_output_process;k++) {
        forward[0][j] *= process[k + 1]->observation[i]->mass[*poutput[k]];
      }
      norm += forward[0][j];
      i++;
    }

    else {
      forward[0][j] = 0.;
    }
  }

  if (norm > 0.) {
    for (i = 0;i < nb_row;i++) {
      forward[0][i] /= norm;
    }
    seq_likelihood += log(norm);
  }

  else {
    seq_likelihood = D_INF;
  }

  if (seq_likelihood != D_INF) {
    for (i = 1;i < seq.length[index];i++) {
      for (j = 0;j < nb_output_process;j++) {
        poutput[j]++;
      }

      for (j = 0;j < order;j++) {
        state_index[j] = 0;
      }
      norm = 0.;

      for (j = 0;j < nb_row;j++) {
        ptransition = transition;
        pforward = forward[i - 1];
        for (k = 0;k < order - 1;k++) {
          ptransition += state_index[k] * power[k + 1];
          pforward += state_index[k] * power[k + 1];
        }

        forward[i][j] = 0.;
        for (k = 0;k < nb_state;k++) {
          forward[i][j] += *(*ptransition + state_index[order - 1]) * *pforward++;
          ptransition++;
        }
        predicted[i][j] = forward[i][j];

        for (k = 0;k < nb_output_process;k++) {
          forward[i][j] *= process[k + 1]->observation[state_index[order - 1]]->mass[*poutput[k]];
        }
        norm += forward[i][j];

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

      if (norm > 0.) {
        for (j = 0;j < nb_row;j++) {
          forward[i][j] /= norm;
        }
        seq_likelihood += log(norm);
      }

      else {
        seq_likelihood = D_INF;
        break;
      }
    }
  }

  // recurrence "backward"

  if (seq_likelihood != D_INF) {
    i = seq.length[index] - 1;
    for (j = 0;j < nb_row;j++) {
      backward[i][j] = forward[i][j];
    }

    for (i = seq.length[index] - 2;i >= 0;i--) {
      for (j = 0;j < nb_row;j++) {
        if (predicted[i + 1][j] > 0.) {
          auxiliary[j] = backward[i + 1][j] / predicted[i + 1][j];
        }
        else {
          auxiliary[j] = 0.;
        }
      }

      for (j = 0;j < nb_row;j++) {
        ptransition = transition;
        pauxiliary = auxiliary;

        ptransition += state_index[0] * power[0];
        for (k = 1;k < order;k++) {
          ptransition += state_index[k] * power[k];
          pauxiliary += state_index[k] * power[k - 1];
        }

        backward[i][j] = 0.;
        for (k = 0;k < nb_state;k++) {
          backward[i][j] += *pauxiliary * *(*ptransition + k);
          pauxiliary += power[order - 1];
        }
        backward[i][j] *= forward[i][j];

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

    // restauration

    pstate = seq.sequence[index][0];

    for (i = 0;i < seq.length[index];i++) {
      for (j = 0;j < nb_state;j++) {
        state_backward[i][j] = 0.;
      }
      k = 0;
      for (j = 0;j < nb_row;j++) {
        state_backward[i][k] += backward[i][j];
        if (k < nb_state - 1) {
          k++;
        }
        else {
          k = 0;
        }
      }

      backward_max = 0.;
      for (j = 0;j < nb_state;j++) {
        if (state_backward[i][j] > backward_max) {
          backward_max = state_backward[i][j];
          *pstate = j;
        }
      }

      pstate++;
    }

    state_seq_likelihood = Markov::likelihood_computation(seq , index);

    switch (format) {

    case 'a' : {
      seq.profile_ascii_print(os , index , nb_state , state_backward ,
                              STAT_label[STATL_STATE]);

      os << "\n" << STAT_label[STATL_LIKELIHOOD] << ": " << seq_likelihood
         << "\n" << SEQ_label[SEQL_STATE_SEQUENCE_LIKELIHOOD] << ": " << state_seq_likelihood
         << "\n" << endl;
      break;
    }

    case 's' : {
      seq.profile_spreadsheet_print(os , index , nb_state , state_backward ,
                                    STAT_label[STATL_STATE]);

      os << "\n" << STAT_label[STATL_LIKELIHOOD] << "\t" << seq_likelihood
         << "\n" << SEQ_label[SEQL_STATE_SEQUENCE_LIKELIHOOD] << "\t" << state_seq_likelihood
         << "\n" << endl;
      break;
    }

    case 'g' : {
      seq.profile_plot_print(os , index , nb_state , state_backward);
      break;
    }
    }
  }

  for (i = 0;i < seq.length[index];i++) {
    delete [] forward[i];
  }
  delete [] forward;

  for (i = 0;i < seq.length[index];i++) {
    delete [] predicted[i];
  }
  delete [] predicted;

  for (i = 0;i < seq.length[index];i++) {
    delete [] backward[i];
  }
  delete [] backward;

  delete [] auxiliary;

  for (i = 0;i < seq.length[index];i++) {
    delete [] state_backward[i];
  }
  delete [] state_backward;

  delete [] poutput;

  return seq_likelihood;
}


/*--------------------------------------------------------------*
 *
 *  Calcul des probabilites a posteriori des etats par l'algorithme forward-backward et
 *  ecriture du resultat.
 *
 *  arguments : reference sur un objet Format_error, stream, sequences,
 *              identificateur de la sequence, format ('a' : ASCII, 's' : Spreadsheet).
 *
 *--------------------------------------------------------------*/

bool Hidden_markov::state_profile_write(Format_error &error , ostream &os ,
                                        const Markov_data &iseq , int identifier ,
                                        char format) const

{
  bool status = true;
  register int i;
  int index = I_DEFAULT;
  Markov_data *seq;


  error.init();

  if (identifier != I_DEFAULT) {
    for (i = 0;i < iseq.nb_sequence;i++) {
      if (identifier == iseq.identifier[i]) {
        index = i;
        break;
      }
    }

    if (i == iseq.nb_sequence) {
      status = false;
      error.update(SEQ_error[SEQR_SEQUENCE_IDENTIFIER]);
    }
  }

  if (status) {
    if (iseq.type[0] == INT_VALUE) {
      seq = new Markov_data((Markovian_sequences&)iseq , 0);
      seq->type[0] = STATE;
    }
    else {
      seq = new Markov_data(iseq , false);
    }

    forward_backward(*seq , index , os , format);

    delete seq;
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Calcul des probabilites a posteriori des etats par l'algorithme forward-backward et
 *  ecriture du resultat.
 *
 *  arguments : reference sur un objet Format_error, stream,
 *              identificateur de la sequence.
 *
 *--------------------------------------------------------------*/

bool Hidden_markov::state_profile_ascii_write(Format_error &error , ostream &os ,
                                              int identifier) const

{
  bool status;


  error.init();

  if (!markov_data) {
    status = false;
    error.update(STAT_error[STATR_NO_DATA]);
  }
  else {
    status = state_profile_write(error , os , *markov_data , identifier , 'a');
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Calcul des probabilites a posteriori des etats par l'algorithme forward-backward et
 *  ecriture du resultat dans un fichier.
 *
 *  arguments : reference sur un objet Format_error, path, identificateur de la sequence,
 *              format de fichier ('a' : ASCII, 's' : Spreadsheet).
 *
 *--------------------------------------------------------------*/

bool Hidden_markov::state_profile_write(Format_error &error , const char *path ,
                                        int identifier , char format) const

{
  bool status = true;
  ofstream out_file(path);


  error.init();

  if (!out_file) {
    status = false;
    error.update(STAT_error[STATR_FILE_NAME]);
  }
  if (!markov_data) {
    status = false;
    error.update(STAT_error[STATR_NO_DATA]);
  }

  if (status) {
    status = state_profile_write(error , out_file , *markov_data , identifier , format);
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Calcul des probabilites a posteriori des etats par l'algorithme forward-backward et
 *  affichage du resultat au format Gnuplot.
 *
 *  arguments : reference sur un objet Format_error, prefixe des fichiers,
 *              sequences, identificateur de la sequence, titre des figures.
 *
 *--------------------------------------------------------------*/

bool Hidden_markov::state_profile_plot_write(Format_error &error , const char *prefix ,
                                             const Markov_data &iseq , int identifier ,
                                             const char *title) const

{
  bool status = true;
  register int i , j;
  int index;
  Markov_data *seq;
  ostringstream data_file_name;
  ofstream *data_out_file;


  error.init();

  for (i = 0;i < iseq.nb_sequence;i++) {
    if (identifier == iseq.identifier[i]) {
      index = i;
      break;
    }
  }

  if (i == iseq.nb_sequence) {
    status = false;
    error.update(SEQ_error[SEQR_SEQUENCE_IDENTIFIER]);
  }

  if (status) {

    // ecriture du fichier de donnees

    data_file_name << prefix << ".dat";
    data_out_file = new ofstream((data_file_name.str()).c_str());

    if (!data_out_file) {
      status = false;
      error.update(STAT_error[STATR_FILE_PREFIX]);
    }

    else {
      if (iseq.type[0] == INT_VALUE) {
        seq = new Markov_data((Markovian_sequences&)iseq , 0);
        seq->type[0] = STATE;
      }
      else {
        seq = new Markov_data(iseq , false);
      }

      forward_backward(*seq , index , *data_out_file , 'g');
      data_out_file->close();
      delete data_out_file;

      // ecriture du fichier de commandes et du fichier d'impression

      for (i = 0;i < 2;i++) {
        ostringstream file_name[2];

        switch (i) {
        case 0 :
          file_name[0] << prefix << ".plot";
          break;
        case 1 :
          file_name[0] << prefix << ".print";
          break;
        }

        ofstream out_file((file_name[0].str()).c_str());

        if (i == 1) {
          out_file << "set terminal postscript" << endl;
          file_name[1] << label(prefix) << ".ps";
          out_file << "set output \"" << file_name[1].str() << "\"\n\n";
        }

        out_file << "set border 15 lw 0\n" << "set tics out\n" << "set xtics nomirror\n"
                 << "set title";
        if (title) {
          out_file << " \"" << title << "\"";
        }
        out_file << "\n\n";

        if (seq->length[index] - 1 < TIC_THRESHOLD) {
          out_file << "set xtics 0,1" << endl;
        }

        out_file << "plot [0:" << seq->length[index] - 1 << "] [0:1] ";
        for (j = 0;j < nb_state;j++) {
          out_file << "\"" << label((data_file_name.str()).c_str()) << "\" using "
                   << j + 1 << " title \"" << STAT_label[STATL_STATE] << " "
                   << j << "\" with linespoints";
          if (j < nb_state - 1) {
            out_file << ",\\";
          }
          out_file << endl;
        }

        if (seq->length[index] - 1 < TIC_THRESHOLD) {
          out_file << "set xtics autofreq" << endl;
        }

        if (i == 1) {
          out_file << "\nset terminal x11" << endl;
        }

        out_file << "\npause 0 \"" << STAT_label[STATL_END] << "\"" << endl;
      }

      delete seq;
    }
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Calcul des probabilites a posteriori des etats par l'algorithme forward-backward et
 *  affichage du resultat au format Gnuplot.
 *
 *  arguments : reference sur un objet Format_error, prefixe des fichiers,
 *              identificateur de la sequence, titre des figures.
 *
 *--------------------------------------------------------------*/

bool Hidden_markov::state_profile_plot_write(Format_error &error , const char *prefix ,
                                             int identifier , const char *title) const

{
  bool status;


  error.init();

  if (!markov_data) {
    status = false;
    error.update(STAT_error[STATR_NO_DATA]);
  }
  else {
    status = state_profile_plot_write(error , prefix , *markov_data , identifier , title);
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Calcul des logarithmes des parametres d'une chaine de Markov cachee
 *
 *--------------------------------------------------------------*/

void Hidden_markov::log_computation()

{
  register int i , j;


  Chain::log_computation();

  for (i = 1;i <= nb_output_process;i++) {
    for (j = 0;j < nb_state;j++) {
      process[i]->observation[j]->log_computation();
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Calcul des sequences d'etats optimale par l'algorithme de Viterbi.
 *
 *  arguments : reference sur un objet Markov_data, indice de la sequence.
 *
 *--------------------------------------------------------------*/

double Hidden_markov::viterbi(const Markov_data &seq , int index) const

{
  register int i , j , k , m;
  int length , *pstate1 , *pstate2 , **poutput , **optimal_state , *poptimal ,
      power[ORDER] , state_index[ORDER];
  double likelihood = 0. , buff , forward_max , **ptransition ,
         *forward , *previous_forward , *pforward;


  // initialisations

  i = 1;
  for (j = 0;j < order;j++) {
    power[j] = i;
    i *= nb_state;
  }

  forward = new double[nb_row];
  previous_forward = new double[nb_row];

  length = (index == I_DEFAULT ? seq.max_length : seq.length[index]);

  optimal_state = new int*[length];
  for (i = 0;i < length;i++) {
    optimal_state[i] = new int[nb_row];
  }

  poutput = new int*[nb_output_process];

  for (i = 0;i < seq.nb_sequence;i++) {
    if ((index == I_DEFAULT) || (index == i)) {
      for (j = 0;j < nb_output_process;j++) {
        poutput[j] = seq.sequence[i][j + 1];
      }

      // recurrence "forward"

      j = 0;
      for (k = 0;k < nb_row;k++) {
        if (k == self_row[j]) {
          forward[k] = cumul_initial[j];
          if (forward[k] != D_INF) {
            for (m = 0;m < nb_output_process;m++) {
              buff = process[m + 1]->observation[j]->cumul[*poutput[m]];
              if (buff == D_INF) {
                forward[k] = D_INF;
                break;
              }
              else {
                forward[k] += buff;
              }
            }
          }
          j++;
        }

        else {
          forward[k] = D_INF;
        }
      }

#     ifdef DEBUG
/*      cout << "\n" << 0 << " : ";
      for (j = 0;j < nb_row;j++) {
        cout << forward[j] << " | ";
      }
      cout << endl; */
#     endif

      for (j = 1;j < seq.length[i];j++) {
        for (k = 0;k < nb_output_process;k++) {
          poutput[k]++;
        }

        for (k = 0;k < nb_row;k++) {
          previous_forward[k] = forward[k];
        }

        for (k = 0;k < order;k++) {
          state_index[k] = 0;
        }

        for (k = 0;k < nb_row;k++) {
          ptransition = cumul_transition;
          pforward = previous_forward;
          for (m = 0;m < order - 1;m++) {
            ptransition += state_index[m] * power[m + 1];
            pforward += state_index[m] * power[m + 1];
          }

          forward[k] = D_INF;
          for (m = 0;m < nb_state;m++) {
            buff = *(*ptransition + state_index[order - 1]) + *pforward++;
            ptransition++;
            if (buff > forward[k]) {
              forward[k] = buff;
              optimal_state[j][k] = m;
            }
          }

          if (forward[k] != D_INF) {
            for (m = 0;m < nb_output_process;m++) {
              buff = process[m + 1]->observation[state_index[order - 1]]->cumul[*poutput[m]];
              if (buff == D_INF) {
                forward[k] = D_INF;
                break;
              }
              else {
                forward[k] += buff;
              }
            }
          }

          // mise a jour des indices des etats

          for (m = 0;m < order;m++) {
            if (state_index[m] < nb_state - 1) {
              state_index[m]++;
              break;
            }
            else {
              state_index[m] = 0;
            }
          }
        }

#       ifdef DEBUG
/*        cout << j << " : ";
        for (k = 0;k < nb_row;k++) {
          cout << forward[k] << " " << optimal_state[j][k] << " | ";
        }
        cout << endl; */
#       endif

      }

      // extraction de la vraisemblance du chemin optimal

      pstate1 = seq.sequence[i][0] + seq.length[i] - 1;
      for (j = 0;j < order;j++) {
        state_index[j] = 0;
      }
      forward_max = D_INF;

      for (j = 0;j < nb_row;j++) {
        if (forward[j] > forward_max) {
          forward_max = forward[j];
          pstate2 = pstate1;
          for (k = 0;k < MIN(order , seq.length[i]);k++) {
            *pstate2-- = state_index[order - 1 - k];
          }
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

      if (forward_max != D_INF) {
        likelihood += forward_max;
      }
      else {
        likelihood = D_INF;
        break;
      }

#     ifdef DEBUG
      if (order == 1) {
        double forward_max2 , *next_forward;


        next_forward = new double[nb_state];

        j = seq.length[i] - 1;
        for (k = 0;k < nb_state;k++) {
          forward[k] = 0.;
          for (m = 0;m < nb_output_process;m++) {
            buff = process[m + 1]->observation[k]->cumul[*poutput[m]];
            if (buff == D_INF) {
              forward[k] = D_INF;
              break;
            }
            else {
              forward[k] += buff;
            }
          }
        }

        for (j = seq.length[i] - 2;j >= 0;j--) {
          for (k = 0;k < nb_output_process;k++) {
            poutput[k]--;
          }

          for (k = 0;k < nb_state;k++) {
            next_forward[k] = forward[k];
          }

          for (k = 0;k < nb_state;k++) {
            forward[k] = D_INF;

            for (m = 0;m < nb_state;m++) {
              buff = next_forward[m] + cumul_transition[k][m];
              if (buff > forward[k]) {
                forward[k] = buff;
              }
            }

            if (forward[k] != D_INF) {
              for (m = 0;m < nb_output_process;m++) {
                buff = process[m + 1]->observation[k]->cumul[*poutput[m]];
                if (buff == D_INF) {
                  forward[k] = D_INF;
                  break;
                }
                else {
                  forward[k] += buff;
                }
              }
            }
          }
        }

        forward_max2 = D_INF;
        for (j = 0;j < nb_state;j++) {
          if (forward[j] + cumul_initial[j] > forward_max2) {
            forward_max2 = forward[j] + cumul_initial[j];
          }
        }

        cout << "\nTest: " << forward_max << " | " << forward_max2 << endl;

        delete [] next_forward;
      }
#     endif

      // restauration

      pstate1 -= (order - 1);
      for (j = seq.length[i] - 1 - order;j >= 0;j--) {
        poptimal = optimal_state[j + order];
        pstate2 = pstate1;
        for (k = 0;k < order;k++) {
          poptimal += *pstate2++ * power[k];
        }
        *--pstate1 = *poptimal;
      }

#     ifdef DEBUG
/*      cout << "\n";
      for (j = seq.length[i] - 1;j >= 0;j--) {
         cout << seq.sequence[i][0][j] << " ";
      }
      cout << endl; */
#     endif

    }
  }

  delete [] forward;
  delete [] previous_forward;

  for (i = 0;i < length;i++) {
    delete [] optimal_state[i];
  }
  delete [] optimal_state;

  delete [] poutput;

  return likelihood;
}


/*--------------------------------------------------------------*
 *
 *  Estimation des parametres d'une chaine de Markov cachee a
 *  partir d'un echantillon de sequences par l'algorithme de Viterbi.
 *
 *  arguments : reference sur un objet Format_error, stream,
 *              chaine de Markov cachee initiale,
 *              flag sur le calcul des lois de comptage, nombre d'iterations.
 *
 *--------------------------------------------------------------*/

Hidden_markov* Markovian_sequences::hidden_markov_viterbi_estimation(Format_error &error , ostream &os ,
                                                                     const Hidden_markov &ihmarkov ,
                                                                     bool counting_flag , int nb_iter) const

{
  bool status;
  register int i , j;
  int iter;
  double likelihood = D_INF , previous_likelihood;
  Hidden_markov *hmarkov;
  Markov_data *seq;


  hmarkov = 0;
  error.init();

  // test nombre de valeurs observees par processus

  status = false;
  for (i = 0;i < nb_variable;i++) {
    if (marginal[i]->nb_value > 1) {
      status = true;
      break;
    }
  }

  if (!status) {
    error.update(SEQ_error[SEQR_VARIABLE_NB_VALUE]);
  }

  for (i = 0;i < nb_variable;i++) {
    if (!characteristics[i]) {
      for (j = 0;j < marginal[i]->nb_value;j++) {
        if (marginal[i]->frequency[j] == 0) {
          status = false;
          ostringstream error_message;
          error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                        << STAT_error[STATR_MISSING_VALUE] << " " << j;
          error.update((error_message.str()).c_str());
        }
      }
    }
  }

  if (status) {
    if (ihmarkov.nb_output_process != nb_variable) {
      status = false;
      error.update(SEQ_error[SEQR_NB_OUTPUT_PROCESS]);
    }

    else {
      for (i = 0;i < nb_variable;i++) {
        if (ihmarkov.process[i + 1]->nb_value != marginal[i]->nb_value) {
          status = false;
          ostringstream error_message;
          error_message << STAT_label[STATL_OUTPUT_PROCESS] << " " << i + 1 << ": "
                        << SEQ_error[SEQR_NB_OUTPUT];
          error.update((error_message.str()).c_str());
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

    hmarkov = new Hidden_markov(ihmarkov , false , false);

    hmarkov->create_cumul();
    hmarkov->log_computation();

    hmarkov->markov_data = new Markov_data(*this , 0);
    seq = hmarkov->markov_data;
    seq->type[0] = STATE;
    seq->max_value[0] = hmarkov->nb_state - 1;
    seq->marginal[0] = new Histogram(hmarkov->nb_state);

    seq->chain_data = new Chain_data(hmarkov->type , hmarkov->nb_state ,
                                     (int)pow((double)hmarkov->nb_state , hmarkov->order));
    seq->create_observation_histogram(hmarkov->nb_state);

    iter = 0;
    do {
      iter++;
      previous_likelihood = likelihood;

      likelihood = hmarkov->viterbi(*seq);

      if (likelihood != D_INF) {
        seq->marginal_histogram_computation(0);

        seq->transition_count_computation(*(seq->chain_data) , hmarkov->order);
        seq->observation_histogram_computation();

#       ifdef DEBUG
        cout << "iteration " << iter << "  "
             << SEQ_label[SEQL_STATE_SEQUENCES_LIKELIHOOD] << ": " << likelihood
             << " | " << hmarkov->Markov::likelihood_computation(*seq) << endl;
#       endif

        // reestimation des probabilites initiales

        reestimation(hmarkov->nb_state , seq->chain_data->initial ,
                     hmarkov->initial , MIN_PROBABILITY , false);

        // reestimation des probabilites de transition

        for (i = 0;i < hmarkov->nb_row;i++) {
          reestimation(hmarkov->nb_state , seq->chain_data->transition[i] ,
                       hmarkov->transition[i] , MIN_PROBABILITY , false);
        }

        // reestimation des probabilites d'observation

        for (i = 1;i <= hmarkov->nb_output_process;i++) {
          for (j = 0;j < hmarkov->nb_state;j++) {
            reestimation(seq->marginal[i]->nb_value , seq->observation[i][j]->frequency ,
                         hmarkov->process[i]->observation[j]->mass , MIN_PROBABILITY , false);
          }
        }

        hmarkov->log_computation();
      }

#     ifdef MESSAGE
      os << STAT_label[STATL_ITERATION] << " " << iter << "   "
         << SEQ_label[SEQL_STATE_SEQUENCES_LIKELIHOOD] << ": " << likelihood << endl;
#     endif

    }
    while ((likelihood != D_INF) && (((nb_iter == I_DEFAULT) && (iter < MARKOV_NB_ITER) &&
             ((likelihood - previous_likelihood) / -likelihood > MARKOV_LIKELIHOOD_DIFF)) ||
            ((nb_iter != I_DEFAULT) && (iter < nb_iter) && (likelihood > previous_likelihood))));

    if (likelihood == D_INF) {
      delete hmarkov;
      hmarkov = 0;
      error.update(STAT_error[STATR_ESTIMATION_FAILURE]);
    }

    else {

#     ifdef MESSAGE
      os << "\n" << iter << " " << STAT_label[STATL_ITERATIONS] << endl;
#     endif

      // reestimation des probabilites initiales

      reestimation(hmarkov->nb_state , seq->chain_data->initial ,
                   hmarkov->initial , MIN_PROBABILITY , true);

      // reestimation des probabilites de transition

      for (i = 0;i < hmarkov->nb_row;i++) {
        reestimation(hmarkov->nb_state , seq->chain_data->transition[i] ,
                     hmarkov->transition[i] , MIN_PROBABILITY , true);
      }

      // reestimation des probabilites d'observation

      for (i = 1;i <= hmarkov->nb_output_process;i++) {
        for (j = 0;j < hmarkov->nb_state;j++) {
          reestimation(seq->marginal[i]->nb_value , seq->observation[i][j]->frequency ,
                       hmarkov->process[i]->observation[j]->mass , MIN_PROBABILITY , true);
        }
      }

      hmarkov->log_computation();

      seq->likelihood = hmarkov->viterbi(*seq);

#     ifdef DEBUG
      cout << "\n" << SEQ_label[SEQL_STATE_SEQUENCES_LIKELIHOOD] << ": " << seq->likelihood << endl;
#     endif

      seq->marginal_histogram_computation(0);
      seq->build_characteristic(0);

      seq->transition_count_computation(*(seq->chain_data) , hmarkov->order);
      seq->observation_histogram_computation();

      hmarkov->remove_cumul();

      for (i = 1;i <= hmarkov->nb_output_process;i++) {
        for (j = 0;j < hmarkov->nb_state;j++) {
          hmarkov->process[i]->observation[j]->cumul_computation();

          hmarkov->process[i]->observation[j]->max_computation();
//          hmarkov->process[i]->observation[j]->mean_computation();
//          hmarkov->process[i]->observation[j]->variance_computation();
        }
      }

      // calcul de la vraisemblance et des lois caracteristiques du modele

      seq->hidden_likelihood = hmarkov->likelihood_computation(*this);

#     ifdef DEBUG
      cout << SEQ_label[SEQL_OBSERVED_SEQUENCES_LIKELIHOOD] << ": " << seq->hidden_likelihood << endl;
#     endif

      hmarkov->component_computation();
      hmarkov->characteristic_computation(*seq , counting_flag , I_DEFAULT , false);
    }
  }

  return hmarkov;
}


/*--------------------------------------------------------------*
 *
 *  Estimation des parametres d'une chaine de Markov cachee a
 *  partir d'un echantillon de sequences par l'algorithme de Viterbi.
 *
 *  arguments : reference sur un objet Format_error, stream, nombre d'etats,
 *              ordre, flag sur le calcul des lois de comptage,
 *              probabilite de rester dans un etat initiale, nombre d'iterations.
 *
 *--------------------------------------------------------------*/

Hidden_markov* Markovian_sequences::hidden_markov_viterbi_estimation(Format_error &error , ostream &os ,
                                                                     int nb_state , int order , bool counting_flag ,
                                                                     double self_transition , int nb_iter) const

{
  bool status = true;
  register int i;
  int nb_value[SEQUENCE_NB_VARIABLE];
  Hidden_markov *ihmarkov , *hmarkov;


  hmarkov = 0;
  error.init();

  if ((nb_state < 2) || (nb_state > NB_STATE)) {
    status = false;
    error.update(SEQ_error[SEQR_NB_STATE]);
  }
  if ((order < 1) || (order > ORDER)) {
    status = false;
    error.update(SEQ_error[SEQR_ORDER]);
  }
  if ((self_transition != D_DEFAULT) && ((self_transition <= 0.) || (self_transition >= 1.))) {
    status = false;
    error.update(SEQ_error[SEQR_SELF_TRANSITION]);
  }

  if (status) {
    for (i = 0;i < nb_variable;i++) {
      nb_value[i] = marginal[i]->nb_value;
    }

    ihmarkov = new Hidden_markov(nb_state , order , nb_variable , nb_value);

    // initialisation des parametres de la chaine de Markov

    if (self_transition == D_DEFAULT) {
      self_transition = MAX(1. - 1. / hlength->mean , SELF_TRANSITION);
    }
    ihmarkov->init(true , self_transition);

    // initialisation des lois d'observation

    for (i = 0;i < ihmarkov->nb_output_process;i++) {
      ihmarkov->process[i + 1]->init();
    }

    hmarkov = hidden_markov_viterbi_estimation(error , os , *ihmarkov , counting_flag , nb_iter);

    delete ihmarkov;
  }

  return hmarkov;
}


/*--------------------------------------------------------------*
 *
 *  Calcul des sequences d'etats optimales.
 *
 *  arguments : references sur un objet Format_error et sur un objet Markovian_sequences,
 *              flag sur le calcul des caracteristiques.
 *
 *--------------------------------------------------------------*/

Markov_data* Hidden_markov::state_sequence_computation(Format_error &error , const Markovian_sequences &iseq ,
                                                       bool characteristic_flag) const

{
  bool status = true;
  register int i;
  Hidden_markov *hmarkov;
  Markov_data *seq;


  seq = 0;
  error.init();

  if (nb_output_process != iseq.nb_variable) {
    status = false;
    error.update(SEQ_error[SEQR_NB_OUTPUT_PROCESS]);
  }

  else {
    for (i = 0;i < nb_output_process;i++) {
      if (process[i + 1]->nb_value < iseq.marginal[i]->nb_value) {
        status = false;
        ostringstream error_message;
        error_message << STAT_label[STATL_OUTPUT_PROCESS] << " " << i + 1 << ": "
                      << SEQ_error[SEQR_NB_OUTPUT];
        error.update((error_message.str()).c_str());
      }
    }
  }

  if (status) {
    seq = new Markov_data(iseq , 0);

    seq->type[0] = STATE;
    seq->markov = new Markov(*this , false , false);

    hmarkov = new Hidden_markov(*this , false , false);

    hmarkov->create_cumul();
    hmarkov->log_computation();

    seq->likelihood = hmarkov->viterbi(*seq);

    delete hmarkov;

    // extraction des caracteristiques des sequences et
    // calcul des lois caracteristiques du modele

    if (seq->likelihood == D_INF) {
      delete seq;
      seq = 0;
      error.update(SEQ_error[SEQR_STATE_SEQUENCE_COMPUTATION_FAILURE]);
    }

    else {
      seq->max_value[0] = nb_state - 1;
      seq->build_marginal_histogram(0);
      seq->build_characteristic(0);

      seq->build_transition_count(order);
      seq->build_observation_histogram();

      if (characteristic_flag) {
        seq->markov->characteristic_computation(*seq , true);
      }
    }
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Comparaison de differentes chaines de Markov cachees pour un ensemble
 *  de sequences.
 *
 *  arguments : reference sur un objet Format_error, stream, nombre de chaines
 *              de Markov cachees, pointeur sur les chaines de Markov cachees,
 *              type algorithme (forward ou Viterbi), path.
 *
 *--------------------------------------------------------------*/

bool Markovian_sequences::comparison(Format_error &error , ostream &os , int nb_model ,
                                     const Hidden_markov **ihmarkov , int algorithm ,
                                     const char *path) const

{
  bool status = true;
  register int i , j;
  double **likelihood;
  Hidden_markov **hmarkov;
  Markov_data *seq;


  error.init();

  for (i = 0;i < nb_variable;i++) {
    if (!characteristics[i]) {
      for (j = 0;j < marginal[i]->nb_value;j++) {
        if (marginal[i]->frequency[j] == 0) {
          status = false;
          ostringstream error_message;
          error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                        << STAT_error[STATR_MISSING_VALUE] << " " << j;
          error.update((error_message.str()).c_str());
        }
      }
    }
  }

  // test nombre de valeurs observees par processus

  for (i = 0;i < nb_model;i++) {
    if (ihmarkov[i]->nb_output_process != nb_variable) {
      status = false;
      ostringstream error_message;
      error_message << SEQ_label[SEQL_HIDDEN_MARKOV_CHAIN] << " " << i + 1 << ": "
                    << SEQ_error[SEQR_NB_OUTPUT_PROCESS];
      error.update((error_message.str()).c_str());
    }

    else {
      for (j = 0;j < nb_variable;j++) {
        if (ihmarkov[i]->process[j + 1]->nb_value < marginal[j]->nb_value) {
          status = false;
          ostringstream error_message;
          error_message << SEQ_label[SEQL_HIDDEN_MARKOV_CHAIN] << " " << i + 1 << ": "
                        << STAT_label[STATL_OUTPUT_PROCESS] << " " << j + 1 << ": "
                        << SEQ_error[SEQR_NB_OUTPUT];
          error.update((error_message.str()).c_str());
        }
      }
    }
  }

  if (status) {
    likelihood = new double*[nb_sequence];
    for (i = 0;i < nb_sequence;i++) {
      likelihood[i] = new double[nb_model];
    }

    if (algorithm == VITERBI) {
      hmarkov = new Hidden_markov*[nb_model];
      for (i = 0;i < nb_model;i++) {
        hmarkov[i] = new Hidden_markov(*(ihmarkov[i]) , false , false);
        hmarkov[i]->create_cumul();
        hmarkov[i]->log_computation();
      }

      seq = new Markov_data(*this , 0);
    }

    // pour chaque sequence, calcul de la vraisemblance (FORWARD) ou de la vraisemblance
    // du chemin optimal (VITERBI) pour chaque modele possible

    for (i = 0;i < nb_sequence;i++) {
      for (j = 0;j < nb_model;j++) {
        switch (algorithm) {
        case FORWARD :
          likelihood[i][j] = ihmarkov[j]->likelihood_computation(*this , i);
          break;
        case VITERBI :
          likelihood[i][j] = hmarkov[j]->viterbi(*seq , i);
          break;
        }
      }
    }

#   ifdef MESSAGE
    likelihood_write(os , nb_model , likelihood , SEQ_label[SEQL_HIDDEN_MARKOV_CHAIN] ,
                     true , algorithm);
#   endif

    if (path) {
      status = likelihood_write(error , path , nb_model , likelihood ,
                                SEQ_label[SEQL_HIDDEN_MARKOV_CHAIN] , algorithm);
    }

    for (i = 0;i < nb_sequence;i++) {
      delete [] likelihood[i];
    }
    delete [] likelihood;

    if (algorithm == VITERBI) {
      for (i = 0;i < nb_model;i++) {
        delete hmarkov[i];
      }
      delete [] hmarkov;

      delete seq;
    }
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Simulation par une chaine de Markov cachee.
 *
 *  arguments : reference sur un objet Format_error,
 *              histogramme des longueurs des sequences,
 *              flag calcul d'une divergence de Kullback-Leibler.
 *
 *--------------------------------------------------------------*/

Markov_data* Hidden_markov::simulation(Format_error &error , const Histogram &hlength ,
                                       bool divergence_flag) const

{
  Markovian_sequences *observ_seq;
  Markov_data *seq;


  seq = Markov::simulation(error , hlength , divergence_flag);

  if ((seq) && (!divergence_flag)) {
    observ_seq = seq->remove_variable_1();
    seq->hidden_likelihood = likelihood_computation(*observ_seq);
    delete observ_seq;
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Simulation par une chaine de Markov cachee.
 *
 *  arguments : reference sur un objet Format_error,
 *              nombre et longueur des sequences.
 *
 *--------------------------------------------------------------*/

Markov_data* Hidden_markov::simulation(Format_error &error , int nb_sequence , int length) const

{
  Markovian_sequences *observ_seq;
  Markov_data *seq;


  seq = Markov::simulation(error , nb_sequence , length);

  if (seq) {
    observ_seq = seq->remove_variable_1();
    seq->hidden_likelihood = likelihood_computation(*observ_seq);
    delete observ_seq;
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Simulation par une chaine de Markov cachee.
 *
 *  arguments : reference sur un objet Format_error, nombre de sequences,
 *              reference sur un objet Markovian_sequences.
 *
 *--------------------------------------------------------------*/

Markov_data* Hidden_markov::simulation(Format_error &error , int nb_sequence ,
                                       const Markovian_sequences &iseq) const

{
  Markovian_sequences *observ_seq;
  Markov_data *seq;


  seq = Markov::simulation(error , nb_sequence , iseq);

  if (seq) {
    observ_seq = seq->remove_variable_1();
    seq->hidden_likelihood = likelihood_computation(*observ_seq);
    delete observ_seq;
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Comparaison de chaines de Markov cachees par calcul de divergences
 *  de Kullback-Leibler.
 *
 *  arguments : reference sur un objet Format_error, stream, nombre de chaines
 *              de Markov cachees, pointeur sur les chaines de Markov cachees,
 *              histogramme des longueurs des sequences, path.
 *
 *--------------------------------------------------------------*/

Distance_matrix* Hidden_markov::divergence_computation(Format_error &error , ostream &os ,
                                                       int nb_model , const Hidden_markov **ihmarkov ,
                                                       Histogram **hlength , const char *path) const

{
  bool status = true , lstatus;
  register int i , j , k;
  int cumul_length;
  double ref_likelihood , target_likelihood , **likelihood;
  const Hidden_markov **hmarkov;
  Markovian_sequences *seq;
  Markov_data *simul_seq;
  Distance_matrix *dist_matrix;
  ofstream *out_file;


  dist_matrix = 0;
  error.init();

  for (i = 0;i < nb_model - 1;i++) {
    if (ihmarkov[i]->nb_output_process != nb_output_process) {
      status = false;
      ostringstream error_message;
      error_message << SEQ_label[SEQL_HIDDEN_MARKOV_CHAIN] << " " << i + 2 << ": "
                    << SEQ_error[SEQR_NB_OUTPUT_PROCESS];
      error.update((error_message.str()).c_str());
    }

    else {
      for (j = 1;j <= nb_output_process;j++) {
        if (ihmarkov[i]->process[j]->nb_value != process[j]->nb_value) {
          status = false;
          ostringstream error_message;
          error_message << SEQ_label[SEQL_HIDDEN_MARKOV_CHAIN] << " " << i + 2 << ": "
                        << STAT_label[STATL_OUTPUT_PROCESS] << " " << j << " "
                        << SEQ_error[SEQR_NB_OUTPUT];
          error.update((error_message.str()).c_str());
        }
      }
    }
  }

  for (i = 0;i < nb_model;i++) {
    lstatus = true;

    if ((hlength[i]->nb_element < 1) || (hlength[i]->nb_element > NB_SEQUENCE)) {
      lstatus = false;
      ostringstream error_message;
      error_message << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_HISTOGRAM] << " "
                    << i + 1 << ": "  << SEQ_error[SEQR_NB_SEQUENCE];
      error.update((error_message.str()).c_str());
    }
    if (hlength[i]->offset < 2) {
      lstatus = false;
      ostringstream error_message;
      error_message << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_HISTOGRAM] << " "
                    << i + 1 << ": "  << SEQ_error[SEQR_SHORT_SEQUENCE_LENGTH];
      error.update((error_message.str()).c_str());
    }
    if (hlength[i]->nb_value - 1 > MAX_LENGTH) {
      lstatus = false;
      ostringstream error_message;
      error_message << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_HISTOGRAM] << " "
                    << i + 1 << ": "  << SEQ_error[SEQR_LONG_SEQUENCE_LENGTH];
      error.update((error_message.str()).c_str());
    }

    if (!lstatus) {
      status = false;
    }

    else {
      cumul_length = 0;
      for (j = hlength[i]->offset;j < hlength[i]->nb_value;j++) {
        cumul_length += j * hlength[i]->frequency[j];
      }

      if (cumul_length > CUMUL_LENGTH) {
        status = false;
        ostringstream error_message;
        error_message << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_HISTOGRAM] << " "
                      << i + 1 << ": "  << SEQ_error[SEQR_CUMUL_SEQUENCE_LENGTH];
        error.update((error_message.str()).c_str());
      }
    }
  }

  if (status) {
    out_file = 0;

    if (path) {
      out_file = new ofstream(path);

      if (!out_file) {
        error.update(STAT_error[STATR_FILE_NAME]);

#       ifdef MESSAGE
        os << error;
#       endif

      }
    }

    hmarkov = new const Hidden_markov*[nb_model];

    hmarkov[0] = this;
    for (i = 1;i < nb_model;i++) {
      hmarkov[i] = ihmarkov[i - 1];
    }

    dist_matrix = new Distance_matrix(nb_model , SEQ_label[SEQL_HIDDEN_MARKOV_CHAIN]);

    for (i = 0;i < nb_model;i++) {

      // simulation d'un echantillon de sequences a partir d'une chaine de Markov cachee

      simul_seq = hmarkov[i]->simulation(error , *hlength[i] , true);
      seq = simul_seq->remove_variable_1();

      likelihood = new double*[seq->nb_sequence];
      for (j = 0;j < seq->nb_sequence;j++) {
        likelihood[j] = new double[nb_model];
      }

      ref_likelihood = 0.;
      for (j = 0;j < seq->nb_sequence;j++) {
        likelihood[j][i] = hmarkov[i]->likelihood_computation(*seq , j);
        ref_likelihood += likelihood[j][i];
      }

      // calcul des vraisemblances de l'echantillon pour chacune des chaines de Markov cachees

      for (j = 0;j < nb_model;j++) {
        if (j != i) {
          target_likelihood = 0.;
          for (k = 0;k < seq->nb_sequence;k++) {
            likelihood[k][j] = hmarkov[j]->likelihood_computation(*seq , k);
            if (target_likelihood != D_INF) {
              if (likelihood[k][j] != D_INF) {
                target_likelihood += likelihood[k][j];
              }
              else {
                target_likelihood = D_INF;
              }
            }
          }

          if (target_likelihood != D_INF) {
            dist_matrix->update(i + 1 , j + 1 , ref_likelihood - target_likelihood , seq->cumul_length);
          }
        }
      }

#     ifdef MESSAGE
      os << SEQ_label[SEQL_HIDDEN_MARKOV_CHAIN] << " " << i + 1 << ": " << seq->nb_sequence << " "
         << SEQ_label[SEQL_SIMULATED] << " " << SEQ_label[seq->nb_sequence == 1 ? SEQL_SEQUENCE : SEQL_SEQUENCES] << endl;
      seq->likelihood_write(os , nb_model , likelihood , SEQ_label[SEQL_HIDDEN_MARKOV_CHAIN]);
#     endif

      if (out_file) {
        *out_file << SEQ_label[SEQL_HIDDEN_MARKOV_CHAIN] << " " << i + 1 << ": " << seq->nb_sequence << " "
                  << SEQ_label[SEQL_SIMULATED] << " " << SEQ_label[seq->nb_sequence == 1 ? SEQL_SEQUENCE : SEQL_SEQUENCES] << endl;
        seq->likelihood_write(*out_file , nb_model , likelihood , SEQ_label[SEQL_HIDDEN_MARKOV_CHAIN]);
      }

      for (j = 0;j < seq->nb_sequence;j++) {
        delete [] likelihood[j];
      }
      delete [] likelihood;

      delete seq;
      delete simul_seq;
    }

    if (out_file) {
      out_file->close();
      delete out_file;
    }

    delete hmarkov;
  }

  return dist_matrix;
}


/*--------------------------------------------------------------*
 *
 *  Comparaison de chaines de Markov cachees par calcul de divergences
 *  de Kullback-Leibler.
 *
 *  arguments : reference sur un objet Format_error, stream, nombre de chaines
 *              de Markov cachees, pointeur sur les chaines de Markov cachees,
 *              nombre et longueur des sequences, path.
 *
 *--------------------------------------------------------------*/

Distance_matrix* Hidden_markov::divergence_computation(Format_error &error , ostream &os ,
                                                       int nb_model , const Hidden_markov **hmarkov ,
                                                       int nb_sequence , int length , const char *path) const

{
  bool status = true;
  register int i;
  Histogram **hlength;
  Distance_matrix *dist_matrix;


  dist_matrix = 0;
  error.init();

  if ((nb_sequence < 1) || (nb_sequence > NB_SEQUENCE)) {
    status = false;
    error.update(SEQ_error[SEQR_NB_SEQUENCE]);
  }
  if (length < 2) {
    status = false;
    error.update(SEQ_error[SEQR_SHORT_SEQUENCE_LENGTH]);
  }
  if (length > MAX_LENGTH) {
    status = false;
    error.update(SEQ_error[SEQR_LONG_SEQUENCE_LENGTH]);
  }

  if (status) {
    hlength = new Histogram*[nb_model];

    hlength[0] = new Histogram(length + 1);

    hlength[0]->nb_element = nb_sequence;
    hlength[0]->offset = length;
    hlength[0]->max = nb_sequence;
    hlength[0]->mean = length;
    hlength[0]->variance = 0.;
    hlength[0]->frequency[length] = nb_sequence;

    for (i = 1;i < nb_model;i++) {
      hlength[i] = new Histogram(*hlength[0]);
    }

    dist_matrix = divergence_computation(error , os , nb_model , hmarkov , hlength , path);

    for (i = 0;i < nb_model;i++) {
      delete hlength[i];
    }
    delete [] hlength;
  }

  return dist_matrix;
}


/*--------------------------------------------------------------*
 *
 *  Comparaison de chaines de Markov cachees par calcul de divergences
 *  de Kullback-Leibler.
 *
 *  arguments : reference sur un objet Format_error, stream, nombre de chaines
 *              de Markov cachees, pointeur sur les chaines de Markov cachees,
 *              pointeurs sur des objets Markovian_sequences, path.
 *
 *--------------------------------------------------------------*/

Distance_matrix* Hidden_markov::divergence_computation(Format_error &error , ostream &os ,
                                                       int nb_model , const Hidden_markov **hmarkov ,
                                                       int nb_sequence , const Markovian_sequences **seq ,
                                                       const char *path) const

{
  register int i;
  Histogram **hlength;
  Distance_matrix *dist_matrix;


  error.init();

  if ((nb_sequence < 1) || (nb_sequence > NB_SEQUENCE)) {
    dist_matrix = 0;
    error.update(SEQ_error[SEQR_NB_SEQUENCE]);
  }

  else {
    hlength = new Histogram*[nb_model];
    for (i = 0;i < nb_model;i++) {
      hlength[i] = seq[i]->hlength->frequency_scale(nb_sequence);
    }

    dist_matrix = divergence_computation(error , os , nb_model , hmarkov , hlength , path);

    for (i = 0;i < nb_model;i++) {
      delete hlength[i];
    }
    delete [] hlength;
  }

  return dist_matrix;
}
