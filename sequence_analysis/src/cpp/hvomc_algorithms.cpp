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
#include <iomanip>
#include "stat_tool/stat_tools.h"
#include "stat_tool/distance_matrix.h"
#include "stat_tool/curves.h"
#include "stat_tool/markovian.h"
#include "stat_tool/stat_label.h"
#include "sequences.h"
#include "variable_order_markov.h"
#include "hidden_variable_order_markov.h"
#include "sequence_label.h"
#include "tool/config.h"

#include "stat_tool/distribution_reestimation.h"   // probleme compilateur C++ Windows

using namespace std;


extern void cumul_computation(int nb_value , const double *pmass , double *pcumul);
extern int cumul_method(int nb_value , const double *cumul , double scale = 1.);
extern void log_computation(int nb_value , const double *pmass , double *plog);

extern int column_width(int value);
extern int column_width(int nb_value , const double *value , double scale = 1.);

extern char* label(const char *file_name);



/*--------------------------------------------------------------*
 *
 *  Calcul de la vraisemblance de sequences pour une chaine de Markov
 *  d'ordre variable cachee par l'algorithme forward.
 *
 *  arguments : reference sur un objet Markovian_sequences, pointeur sur
 *              les probabilites a posteriori des sequences d'etats
 *              les plus probables, indice de la sequence.
 *
 *--------------------------------------------------------------*/

double Hidden_variable_order_markov::likelihood_computation(const Markovian_sequences &seq ,
                                                            double *posterior_probability , int index) const

{
  register int i , j , k , m;
  int nb_value , **poutput;
  double likelihood = 0. , seq_likelihood , *forward , *auxiliary , norm;


  // verification de la compatibilite entre le modele et les donnees

  if (nb_output_process == seq.nb_variable) {
    for (i = 0;i < nb_output_process;i++) {
      if (nonparametric_process[i + 1]) {
        nb_value = nonparametric_process[i + 1]->nb_value;
      }
      else {
        nb_value = parametric_process[i + 1]->nb_value;
      }

      if (nb_value < seq.marginal[i]->nb_value) {
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

    forward = new double[nb_row];
    auxiliary = new double[nb_row];

    poutput = new int*[seq.nb_variable];

    for (i = 0;i < seq.nb_sequence;i++) {
      if ((index == I_DEFAULT) || (index == i)) {
        for (j = 0;j < seq.nb_variable;j++) {
          poutput[j] = seq.sequence[i][j];
        }
        seq_likelihood = 0.;

        norm = 0.;

        switch (type) {

        case 'o' : {
          for (j = 1;j < nb_row;j++) {
            if (order[j] == 1) {
              forward[j] = initial[state[j][0]];

              for (k = 0;k < nb_output_process;k++) {
                if (nonparametric_process[k + 1]) {
                  forward[j] *= nonparametric_process[k + 1]->observation[state[j][0]]->mass[*poutput[k]];
                }
                else {
                  forward[j] *= parametric_process[k + 1]->observation[state[j][0]]->mass[*poutput[k]];
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
                if (nonparametric_process[k + 1]) {
                  forward[j] *= nonparametric_process[k + 1]->observation[state[j][0]]->mass[*poutput[k]];
                }
                else {
                  forward[j] *= parametric_process[k + 1]->observation[state[j][0]]->mass[*poutput[k]];
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
            poutput[k]++;
          }
          norm = 0.;

          for (k = 1;k < nb_row;k++) {
            auxiliary[k] = 0.;
            for (m = 0;m < nb_memory[k];m++) {
              auxiliary[k] += transition[previous[k][m]][state[k][0]] * forward[previous[k][m]];
            }

            for (m = 0;m < nb_output_process;m++) {
              if (nonparametric_process[m + 1]) {
                auxiliary[k] *= nonparametric_process[m + 1]->observation[state[k][0]]->mass[*poutput[m]];
              }
              else {
                auxiliary[k] *= parametric_process[m + 1]->observation[state[k][0]]->mass[*poutput[m]];
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

    delete [] poutput;
  }

  return likelihood;
}


/*--------------------------------------------------------------*
 *
 *  Estimation des parametres d'une chaine de Markov d'ordre variable cachee
 *  a partir d'un echantillon de sequences par l'algorithme EM.
 *
 *  arguments : reference sur un objet Format_error, stream, chaine de Markov cachee initiale,
 *              type d'estimation des probabilites de transition initiale (cas ordinaire),
 *              flags sur le calcul des lois de comptage et sur le calcul des sequences
 *              d'etats optimales, nombre d'iterations.
 *
 *--------------------------------------------------------------*/

Hidden_variable_order_markov* Markovian_sequences::hidden_variable_order_markov_estimation(Format_error &error , ostream &os ,
                                                                                           const Hidden_variable_order_markov &ihmarkov ,
                                                                                           bool global_initial_transition ,
                                                                                           bool counting_flag , bool state_sequence ,
                                                                                           int nb_iter) const

{
  bool status;
  register int i , j , k , m;
  int nb_terminal , max_nb_value , iter , **poutput;
  double likelihood = D_INF , previous_likelihood , observation_likelihood , min_likelihood ,
         **forward , norm , **predicted , buff , *backward , *auxiliary , *reestim;
  Chain_reestimation<double> *chain_reestim;
  Reestimation<double> ***observation_reestim;
  Histogram *hobservation;
  Hidden_variable_order_markov *hmarkov;
  Variable_order_markov_data *seq;


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

  if (ihmarkov.nb_output_process != nb_variable) {
    status = false;
    error.update(SEQ_error[SEQR_NB_OUTPUT_PROCESS]);
  }

  else {
    for (i = 0;i < nb_variable;i++) {
      if (((ihmarkov.nonparametric_process[i + 1]) &&
           (ihmarkov.nonparametric_process[i + 1]->nb_value != marginal[i]->nb_value)) ||
          ((ihmarkov.parametric_process[i + 1]) &&
           (ihmarkov.parametric_process[i + 1]->nb_value < marginal[i]->nb_value))) {
        status = false;
        ostringstream error_message;
        error_message << STAT_label[STATL_OUTPUT_PROCESS] << " " << i + 1 << ": "
                      << SEQ_error[SEQR_NB_OUTPUT];
        error.update((error_message.str()).c_str());
      }

      else if ((ihmarkov.nonparametric_process[i + 1]) && (!characteristics[i])) {
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
  }

  if ((nb_iter != I_DEFAULT) && (nb_iter < 1)) {
    status = false;
    error.update(STAT_error[STATR_NB_ITERATION]);
  }

  if (status) {

    // creation de la chaine de Markov cachee

    hmarkov = new Hidden_variable_order_markov(ihmarkov , false);

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

#   ifdef DEBUG
    cout << *hmarkov;
#   endif

    // initialisations

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

    chain_reestim = new Chain_reestimation<double>((hmarkov->type == 'o' ?  'o' : 'v') ,
                                                   hmarkov->nb_state , hmarkov->nb_row);

    observation_reestim = new Reestimation<double>**[hmarkov->nb_output_process];
    for (i = 0;i < hmarkov->nb_output_process;i++) {
      observation_reestim[i] = new Reestimation<double>*[hmarkov->nb_state];
      for (j = 0;j < hmarkov->nb_state;j++) {
        observation_reestim[i][j] = new Reestimation<double>(marginal[i]->nb_value);
      }
    }

    max_nb_value = 0;
    for (i = 0;i < hmarkov->nb_output_process;i++) {
      if ((hmarkov->parametric_process[i + 1]) && (max_nb_value < marginal[i]->nb_value)) {
        max_nb_value = marginal[i]->nb_value;
      }
    }

    if (max_nb_value > 0) {
      hobservation = new Histogram(max_nb_value);
    }
    else {
      hobservation = 0;
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
          reestim = observation_reestim[i][j]->frequency;
          for (k = 0;k < marginal[i]->nb_value;k++) {
            *reestim++ = 0.;
          }
        }
      }

      for (i = 0;i < nb_sequence;i++) {

        // recurrence "forward"

        for (j = 0;j < nb_variable;j++) {
          poutput[j] = sequence[i][j];
        }
        norm = 0.;

        switch (hmarkov->type) {

        case 'o' : {
          for (j = 1;j < hmarkov->nb_row;j++) {
            if (hmarkov->order[j] == 1) {
              forward[0][j] = hmarkov->initial[hmarkov->state[j][0]];

              for (k = 0;k < hmarkov->nb_output_process;k++) {
                if (hmarkov->nonparametric_process[k + 1]) {
                  forward[0][j] *= hmarkov->nonparametric_process[k + 1]->observation[hmarkov->state[j][0]]->mass[*poutput[k]];
                }
                else {
                  forward[0][j] *= hmarkov->parametric_process[k + 1]->observation[hmarkov->state[j][0]]->mass[*poutput[k]];
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
                if (hmarkov->nonparametric_process[k + 1]) {
                  forward[0][j] *= hmarkov->nonparametric_process[k + 1]->observation[hmarkov->state[j][0]]->mass[*poutput[k]];
                }
                else {
                  forward[0][j] *= hmarkov->parametric_process[k + 1]->observation[hmarkov->state[j][0]]->mass[*poutput[k]];
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
            poutput[k]++;
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
              if (hmarkov->nonparametric_process[m + 1]) {
                forward[j][k] *= hmarkov->nonparametric_process[m + 1]->observation[hmarkov->state[k][0]]->mass[*poutput[m]];
              }
              else {
                forward[j][k] *= hmarkov->parametric_process[m + 1]->observation[hmarkov->state[k][0]]->mass[*poutput[m]];
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
            observation_reestim[m][hmarkov->state[k][0]]->frequency[*poutput[m]] += backward[k];
          }
        }

        for (j = length[i] - 2;j >= 0;j--) {
          for (k = 0;k < nb_variable;k++) {
            poutput[k]--;
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
                observation_reestim[m][hmarkov->state[k][0]]->frequency[*poutput[m]] += backward[k];
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
          if (hmarkov->nonparametric_process[i + 1]) {
            for (j = 0;j < hmarkov->nb_state;j++) {
              reestimation(marginal[i]->nb_value , observation_reestim[i][j]->frequency ,
                           hmarkov->nonparametric_process[i + 1]->observation[j]->mass ,
                           MIN_PROBABILITY , false);
            }
          }

          else {
            for (j = 0;j < hmarkov->nb_state;j++) {
              observation_reestim[i][j]->nb_value_computation();
              observation_reestim[i][j]->offset_computation();
              observation_reestim[i][j]->nb_element_computation();
              observation_reestim[i][j]->max_computation();
              observation_reestim[i][j]->mean_computation();
              observation_reestim[i][j]->variance_computation();

              hobservation->update(observation_reestim[i][j] ,
                                   MAX((int)(observation_reestim[i][j]->nb_element *
                                             MAX(sqrt(observation_reestim[i][j]->variance) , 1.) * OBSERVATION_COEFF) , MIN_NB_ELEMENT));
              observation_likelihood = hobservation->Reestimation<int>::type_parametric_estimation(hmarkov->parametric_process[i + 1]->observation[j] ,
                                                                                                   0 , true , OBSERVATION_THRESHOLD);

              if (observation_likelihood == D_INF) {
                min_likelihood = D_INF;
              }
              else {
                hmarkov->parametric_process[i + 1]->observation[j]->computation(marginal[i]->nb_value ,
                                                                                OBSERVATION_THRESHOLD);

                if (hmarkov->parametric_process[i + 1]->observation[j]->ident == BINOMIAL) {
                  for (k = hmarkov->parametric_process[i + 1]->observation[j]->nb_value;k < marginal[i]->nb_value;k++) {
                    hmarkov->parametric_process[i + 1]->observation[j]->mass[k] = 0.;
                  }
                }
              }
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

      // reestimation des lois d'observation non-parametriques

      for (i = 0;i < hmarkov->nb_output_process;i++) {
        if (hmarkov->nonparametric_process[i + 1]) {
          for (j = 0;j < hmarkov->nb_state;j++) {
            reestimation(marginal[i]->nb_value , observation_reestim[i][j]->frequency ,
                         hmarkov->nonparametric_process[i + 1]->observation[j]->mass ,
                         MIN_PROBABILITY , true);
          }
        }

        else {
          hmarkov->parametric_process[i + 1]->nb_value_computation();
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
      for (j = 0;j < hmarkov->nb_state;j++) {
        delete observation_reestim[i][j];
      }
      delete [] observation_reestim[i];
    }
    delete [] observation_reestim;

    delete hobservation;

    delete [] poutput;

    if (likelihood == D_INF) {
      delete hmarkov;
      hmarkov = 0;
      error.update(STAT_error[STATR_ESTIMATION_FAILURE]);
    }

    else {
      if (state_sequence) {
        hmarkov->markov_data = new Variable_order_markov_data(*this , 0 , (hmarkov->type == 'e' ? true : false));
        seq = hmarkov->markov_data;
        seq->type[0] = STATE;

        for (i = 1;i <= hmarkov->nb_output_process;i++) {
          if ((hmarkov->parametric_process[i]) && (seq->characteristics[i])) {
            delete seq->characteristics[i];
            seq->characteristics[i] = 0;
          }
        }

        hmarkov->create_cumul();
        hmarkov->log_computation();

        seq->posterior_probability = new double[seq->nb_sequence];
        seq->likelihood = hmarkov->viterbi(*seq , seq->posterior_probability);

        hmarkov->remove_cumul();

        seq->max_value[0] = hmarkov->nb_state - 1;
        seq->build_marginal_histogram(0);
        seq->build_characteristic(0);

        seq->build_transition_count(*hmarkov);
        seq->build_observation_histogram();

        // calcul des lois d'observation parametriques

        for (i = 1;i <= hmarkov->nb_output_process;i++) {
          if (hmarkov->parametric_process[i]) {
            for (j = 0;j < hmarkov->nb_state;j++) {
              hmarkov->parametric_process[i]->observation[j]->computation(seq->observation[i][j]->nb_value ,
                                                                          OBSERVATION_THRESHOLD);
            }
          }
        }

#       ifdef MESSAGE
        cout << "\n" << SEQ_label[SEQL_STATE_SEQUENCES_LIKELIHOOD] << ": " << seq->likelihood
             << " | " << hmarkov->Variable_order_markov::likelihood_computation(*seq) << endl;
#       endif

      }

      else {
        hmarkov->markov_data = new Variable_order_markov_data(*this , (hmarkov->type == 'e' ? true : false));
        seq = hmarkov->markov_data;
        seq->state_variable_init(INT_VALUE);

        for (i = 1;i <= hmarkov->nb_output_process;i++) {
          if ((hmarkov->parametric_process[i]) && (seq->characteristics[i - 1])) {
            delete seq->characteristics[i - 1];
            seq->characteristics[i - 1] = 0;
          }
        }
      }

      for (i = 1;i <= hmarkov->nb_output_process;i++) {
        if (hmarkov->nonparametric_process[i]) {
          for (j = 0;j < hmarkov->nb_state;j++) {
            hmarkov->nonparametric_process[i]->observation[j]->cumul_computation();

            hmarkov->nonparametric_process[i]->observation[j]->max_computation();
//            hmarkov->nonparametric_process[i]->observation[j]->mean_computation();
//            hmarkov->nonparametric_process[i]->observation[j]->variance_computation();
          }
        }
      }

      // calcul de la vraisemblance et des lois caracteristiques du modele

      seq->hidden_likelihood = hmarkov->likelihood_computation(*this , seq->posterior_probability);

#     ifdef DEBUG
//      cout << *hmarkov;
      cout << "iteration " << iter << "  "
           << SEQ_label[SEQL_OBSERVED_SEQUENCES_LIKELIHOOD] << ": " << seq->hidden_likelihood << endl;
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
    }
  }

  return hmarkov;
}


/*--------------------------------------------------------------*
 *
 *  Estimation des parametres d'une chaine de Markov d'ordre variable cachee
 *  a partir d'un echantillon de sequences par l'algorithme SEM/MCEM.
 *
 *  arguments : reference sur un objet Format_error, stream, chaine de Markov cachee initiale,
 *              type d'estimation des probabilites de transition initiale (cas ordinaire),
 *              parametres pour le nombre de sequences d'etats simulees, flags sur le calcul
 *              des lois de comptage et sur le calcul des sequences d'etats optimales,
 *              nombre d'iterations.
 *
 *--------------------------------------------------------------*/

Hidden_variable_order_markov* Markovian_sequences::hidden_variable_order_markov_stochastic_estimation(Format_error &error , ostream &os ,
                                                                                                      const Hidden_variable_order_markov &ihmarkov ,
                                                                                                      bool global_initial_transition ,
                                                                                                      int min_nb_state_sequence ,
                                                                                                      int max_nb_state_sequence ,
                                                                                                      double parameter , bool counting_flag ,
                                                                                                      bool state_sequence , int nb_iter) const

{
  bool status;
  register int i , j , k , m;
  int nb_terminal , iter , nb_state_sequence , memory , *state_seq , *pstate , **poutput;
  double likelihood = D_INF , previous_likelihood , observation_likelihood , min_likelihood ,
         **forward , norm , **predicted , *backward , *cumul_backward , *reestim;
  Chain_reestimation<double> *chain_reestim;
  Reestimation<double> ***observation_reestim;
  Hidden_variable_order_markov *hmarkov;
  Variable_order_markov_data *seq;

# ifdef DEBUG
  double sum;
# endif


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

  if (ihmarkov.nb_output_process != nb_variable) {
    status = false;
    error.update(SEQ_error[SEQR_NB_OUTPUT_PROCESS]);
  }

  else {
    for (i = 0;i < nb_variable;i++) {
      if (((ihmarkov.nonparametric_process[i + 1]) &&
           (ihmarkov.nonparametric_process[i + 1]->nb_value != marginal[i]->nb_value)) ||
          ((ihmarkov.parametric_process[i + 1]) &&
           (ihmarkov.parametric_process[i + 1]->nb_value < marginal[i]->nb_value))) {
        status = false;
        ostringstream error_message;
        error_message << STAT_label[STATL_OUTPUT_PROCESS] << " " << i + 1 << ": "
                      << SEQ_error[SEQR_NB_OUTPUT];
        error.update((error_message.str()).c_str());
      }

      else if ((ihmarkov.nonparametric_process[i + 1]) && (!characteristics[i])) {
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
  }

  if ((nb_iter != I_DEFAULT) && (nb_iter < 1)) {
    status = false;
    error.update(STAT_error[STATR_NB_ITERATION]);
  }

  if ((min_nb_state_sequence < 1) || (min_nb_state_sequence > max_nb_state_sequence)) {
    status = false;
    error.update(SEQ_error[SEQR_MIN_NB_STATE_SEQUENCE]);
  }

  if (status) {

    // creation de la chaine de Markov cachee

    hmarkov = new Hidden_variable_order_markov(ihmarkov , false);

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

#   ifdef DEBUG
    cout << *hmarkov;
#   endif

    // initialisations

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

    chain_reestim = new Chain_reestimation<double>((hmarkov->type == 'o' ?  'o' : 'v') ,
                                                   hmarkov->nb_state , hmarkov->nb_row);

    observation_reestim = new Reestimation<double>**[hmarkov->nb_output_process];
    for (i = 0;i < hmarkov->nb_output_process;i++) {
      observation_reestim[i] = new Reestimation<double>*[hmarkov->nb_state];
      for (j = 0;j < hmarkov->nb_state;j++) {
        observation_reestim[i][j] = new Reestimation<double>(marginal[i]->nb_value);
      }
    }

    poutput = new int*[nb_variable];

    iter = 0;
    do {
      previous_likelihood = likelihood;
      likelihood = 0.;

      // calcul du nombre de sequences d'etats simulees

      if (min_nb_state_sequence + (int)round(parameter * iter) < max_nb_state_sequence) {
        nb_state_sequence = min_nb_state_sequence + (int)round(parameter * iter);
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
        for (j = 0;j < hmarkov->nb_state;j++) {
          reestim = observation_reestim[i][j]->frequency;
          for (k = 0;k < marginal[i]->nb_value;k++) {
            *reestim++ = 0.;
          }
        }
      }

      for (i = 0;i < nb_sequence;i++) {

        // recurrence "forward"

        for (j = 0;j < nb_variable;j++) {
          poutput[j] = sequence[i][j];
        }
        norm = 0.;

        switch (hmarkov->type) {

        case 'o' : {
          for (j = 1;j < hmarkov->nb_row;j++) {
            if (hmarkov->order[j] == 1) {
              forward[0][j] = hmarkov->initial[hmarkov->state[j][0]];

              for (k = 0;k < hmarkov->nb_output_process;k++) {
                if (hmarkov->nonparametric_process[k + 1]) {
                  forward[0][j] *= hmarkov->nonparametric_process[k + 1]->observation[hmarkov->state[j][0]]->mass[*poutput[k]];
                }
                else {
                  forward[0][j] *= hmarkov->parametric_process[k + 1]->observation[hmarkov->state[j][0]]->mass[*poutput[k]];
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
                if (hmarkov->nonparametric_process[k + 1]) {
                  forward[0][j] *= hmarkov->nonparametric_process[k + 1]->observation[hmarkov->state[j][0]]->mass[*poutput[k]];
                }
                else {
                  forward[0][j] *= hmarkov->parametric_process[k + 1]->observation[hmarkov->state[j][0]]->mass[*poutput[k]];
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
            poutput[k]++;
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
              if (hmarkov->nonparametric_process[m + 1]) {
                forward[j][k] *= hmarkov->nonparametric_process[m + 1]->observation[hmarkov->state[k][0]]->mass[*poutput[m]];
              }
              else {
                forward[j][k] *= hmarkov->parametric_process[m + 1]->observation[hmarkov->state[k][0]]->mass[*poutput[m]];
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
            poutput[m] = sequence[i][m] + k;
          }

          cumul_computation(hmarkov->nb_row - 1 , forward[k] + 1 , cumul_backward);
          memory = 1 + cumul_method(hmarkov->nb_row - 1 , cumul_backward);
          *pstate = hmarkov->state[memory][0];

          // accumulation des quantites de reestimation des lois d'observation

          for (m = 0;m < hmarkov->nb_output_process;m++) {
            (observation_reestim[m][*pstate]->frequency[*poutput[m]])++;
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
              (observation_reestim[m][*pstate]->frequency[*--poutput[m]])++;
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
          if (hmarkov->nonparametric_process[i + 1]) {
            for (j = 0;j < hmarkov->nb_state;j++) {
              reestimation(marginal[i]->nb_value , observation_reestim[i][j]->frequency ,
                           hmarkov->nonparametric_process[i + 1]->observation[j]->mass ,
                           MIN_PROBABILITY , false);
            }
          }

          else {
            for (j = 0;j < hmarkov->nb_state;j++) {
              observation_reestim[i][j]->nb_value_computation();
              observation_reestim[i][j]->offset_computation();
              observation_reestim[i][j]->nb_element_computation();
              observation_reestim[i][j]->max_computation();
              observation_reestim[i][j]->mean_computation();
              observation_reestim[i][j]->variance_computation();

              observation_likelihood = observation_reestim[i][j]->type_parametric_estimation(hmarkov->parametric_process[i + 1]->observation[j] ,
                                                                                             0 , true , OBSERVATION_THRESHOLD);

              if (observation_likelihood == D_INF) {
                min_likelihood = D_INF;
              }
              else {
                hmarkov->parametric_process[i + 1]->observation[j]->computation(marginal[i]->nb_value ,
                                                                                OBSERVATION_THRESHOLD);

                if (hmarkov->parametric_process[i + 1]->observation[j]->ident == BINOMIAL) {
                  for (k = hmarkov->parametric_process[i + 1]->observation[j]->nb_value;k < marginal[i]->nb_value;k++) {
                    hmarkov->parametric_process[i + 1]->observation[j]->mass[k] = 0.;
                  }
                }
              }
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

      // reestimation des lois d'observation non-parametriques

      for (i = 0;i < hmarkov->nb_output_process;i++) {
        if (hmarkov->nonparametric_process[i + 1]) {
          for (j = 0;j < hmarkov->nb_state;j++) {
            reestimation(marginal[i]->nb_value , observation_reestim[i][j]->frequency ,
                         hmarkov->nonparametric_process[i + 1]->observation[j]->mass ,
                         MIN_PROBABILITY , true);
          }
        }

        else {
          hmarkov->parametric_process[i + 1]->nb_value_computation();
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
      for (j = 0;j < hmarkov->nb_state;j++) {
        delete observation_reestim[i][j];
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
      if  (state_sequence) {
        hmarkov->markov_data = new Variable_order_markov_data(*this , 0 , (hmarkov->type == 'e' ? true : false));
        seq = hmarkov->markov_data;
        seq->type[0] = STATE;

        for (i = 1;i <= hmarkov->nb_output_process;i++) {
          if ((hmarkov->parametric_process[i]) && (seq->characteristics[i])) {
            delete seq->characteristics[i];
            seq->characteristics[i] = 0;
          }
        }

        hmarkov->create_cumul();
        hmarkov->log_computation();

        seq->posterior_probability = new double[seq->nb_sequence];
        seq->likelihood = hmarkov->viterbi(*seq , seq->posterior_probability);

        hmarkov->remove_cumul();

        seq->max_value[0] = hmarkov->nb_state - 1;
        seq->build_marginal_histogram(0);
        seq->build_characteristic(0);

        seq->build_transition_count(*hmarkov);
        seq->build_observation_histogram();

        // calcul des lois d'observation parametriques

        for (i = 1;i <= hmarkov->nb_output_process;i++) {
          if (hmarkov->parametric_process[i]) {
            for (j = 0;j < hmarkov->nb_state;j++) {
              hmarkov->parametric_process[i]->observation[j]->computation(seq->observation[i][j]->nb_value ,
                                                                          OBSERVATION_THRESHOLD);
            }
          }
        }

#       ifdef MESSAGE
        cout << "\n" << SEQ_label[SEQL_STATE_SEQUENCES_LIKELIHOOD] << ": " << seq->likelihood
             << " | " << hmarkov->Variable_order_markov::likelihood_computation(*seq) << endl;
#       endif

      }

      else {
        hmarkov->markov_data = new Variable_order_markov_data(*this , (hmarkov->type == 'e' ? true : false));
        seq = hmarkov->markov_data;
        seq->state_variable_init(INT_VALUE);

        for (i = 1;i <= hmarkov->nb_output_process;i++) {
          if ((hmarkov->parametric_process[i]) && (seq->characteristics[i - 1])) {
            delete seq->characteristics[i - 1];
            seq->characteristics[i - 1] = 0;
          }
        }
      }

      for (i = 1;i <= hmarkov->nb_output_process;i++) {
        if (hmarkov->nonparametric_process[i]) {
          for (j = 0;j < hmarkov->nb_state;j++) {
            hmarkov->nonparametric_process[i]->observation[j]->cumul_computation();

            hmarkov->nonparametric_process[i]->observation[j]->max_computation();
//            hmarkov->nonparametric_process[i]->observation[j]->mean_computation();
//            hmarkov->nonparametric_process[i]->observation[j]->variance_computation();
          }
        }
      }

      // calcul de la vraisemblance et des lois caracteristiques du modele

      seq->hidden_likelihood = hmarkov->likelihood_computation(*this , seq->posterior_probability);

#     ifdef DEBUG
//      cout << *hmarkov;
      cout << "iteration " << iter << "  "
           << SEQ_label[SEQL_OBSERVED_SEQUENCES_LIKELIHOOD] << ": " << seq->hidden_likelihood << endl;
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
    }
  }

  return hmarkov;
}


/*--------------------------------------------------------------
 *
 *  Ecriture des profils d'etats et des profils d'entropie pour une sequence.
 *
 *  arguments : stream, indice de la sequence, nombre d'etats,
 *              pointeur sur les profils d'etats, 
 *              pointeurs sur les profils d'entropie.
 *
 *--------------------------------------------------------------*/

ostream& Sequences::profile_ascii_print(ostream &os , int index , int nb_state ,
                                        double **profiles , double *begin_conditional_entropy ,
                                        double *marginal_entropy , double *begin_partial_entropy ,
                                        double *end_conditional_entropy , double *end_partial_entropy) const

{
  register int i , j;
  int buff , *width;
  long old_adjust;


  old_adjust = os.flags(ios::adjustfield);

  // calcul des largeurs des colonnes

  width = new int[nb_variable + 8];

  for (i = 0;i <  nb_variable;i++) {
    width[i] = column_width(max_value[i]);
    if (i > 0) {
      width[i] += ASCII_SPACE;
    }
  }

  width[nb_variable] = column_width(max_length) + ASCII_SPACE;

  width[nb_variable + 1] = 0;
  for (i = 0;i < length[index];i++) {
    buff = column_width(nb_state , profiles[i]);
    if (buff > width[nb_variable + 1]) {
      width[nb_variable + 1] = buff;
    }
  }
  width[nb_variable + 1] += ASCII_SPACE;

  width[nb_variable + 2] = column_width(length[index] , begin_conditional_entropy) + ASCII_SPACE;
  width[nb_variable + 3] = column_width(length[index] , marginal_entropy) + ASCII_SPACE;
  width[nb_variable + 4] = column_width(length[index] , begin_partial_entropy) + ASCII_SPACE;

  width[nb_variable + 5] = column_width(nb_sequence);

  if ((end_conditional_entropy) && (end_partial_entropy)) {
    width[nb_variable + 6] = column_width(length[index] , end_conditional_entropy) + ASCII_SPACE;
    width[nb_variable + 7] = column_width(length[index] , end_partial_entropy) + ASCII_SPACE;
  }

  os << SEQ_label[SEQL_OPTIMAL] << " " << STAT_label[STATL_STATE];
  for (i = 1;i < nb_variable;i++) {
    os << " | " << STAT_label[STATL_VARIABLE] << " " << i;
  }
  os << " | " << SEQ_label[SEQL_INDEX];
  for (i = 0;i < nb_state;i++) {
    os << " | " << STAT_label[STATL_STATE] << " " << i;
  }
  os << " | " << SEQ_label[SEQL_CONDITIONAL_ENTROPY] << " | " << SEQ_label[SEQL_MARGINAL_ENTROPY]
     << " | " << SEQ_label[SEQL_PARTIAL_STATE_SEQUENCE_ENTROPY] << endl;

  for (i = 0;i < length[index];i++) {
    os.setf(ios::right , ios::adjustfield);
    for (j = 0;j < nb_variable;j++) {
      os << setw(width[j]) << sequence[index][j][i];
    }
    os << setw(width[nb_variable]) << i << "  ";

    os.setf(ios::left , ios::adjustfield);
    for (j = 0;j < nb_state;j++) {
      os << setw(width[nb_variable + 1]) << profiles[i][j];
    }

    if ((end_conditional_entropy) && (end_partial_entropy)) {
      os << setw(width[nb_variable + 2]) << begin_conditional_entropy[i];
      os << setw(width[nb_variable + 6]) << end_conditional_entropy[i];
      os << setw(width[nb_variable + 3]) << marginal_entropy[i];
      os << setw(width[nb_variable + 4]) << begin_partial_entropy[i];
      os << setw(width[nb_variable + 7]) << end_partial_entropy[i];
    }
    else {
      os << setw(width[nb_variable + 2]) << begin_conditional_entropy[i];
      os << setw(width[nb_variable + 3]) << marginal_entropy[i];
      os << setw(width[nb_variable + 4]) << begin_partial_entropy[i];
    }

    if (i == 0) {
      os.setf(ios::right , ios::adjustfield);
      os << setw(width[nb_variable + 5]) << identifier[index];
    }
    os << endl;
  }

  delete [] width;

  os.setf((FMTFLAGS)old_adjust , ios::adjustfield);

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture des profils d'etats et des profils d'entropie pour une sequence
 *  au format tableur.
 *
 *  arguments : stream, indice de la sequence, nombre d'etats,
 *              pointeur sur les profils d'etats, 
 *              pointeurs sur les profils d'entropie.
 *
 *--------------------------------------------------------------*/

ostream& Sequences::profile_spreadsheet_print(ostream &os , int index , int nb_state ,
                                              double **profiles , double *begin_conditional_entropy ,
                                              double *marginal_entropy , double *begin_partial_entropy ,
                                              double *end_conditional_entropy , double *end_partial_entropy) const

{
  register int i , j;


  os << SEQ_label[SEQL_OPTIMAL] << " " << STAT_label[STATL_STATE];
  for (i = 1;i < nb_variable;i++) {
    os << "\t" << STAT_label[STATL_VARIABLE] << " " << i;
  }
  os << "\t" << SEQ_label[SEQL_INDEX];
  for (i = 0;i < nb_state;i++) {
    os << "\t" << STAT_label[STATL_STATE] << " " << i;
  }
  os << "\t" << SEQ_label[SEQL_CONDITIONAL_ENTROPY] << "\t" << SEQ_label[SEQL_MARGINAL_ENTROPY]
     << "\t" << SEQ_label[SEQL_PARTIAL_STATE_SEQUENCE_ENTROPY] << endl;

  for (i = 0;i < length[index];i++) {
    for (j = 0;j < nb_variable;j++) {
      os << sequence[index][j][i] << "\t";
    }
    os << i;
    for (j = 0;j < nb_state;j++) {
      os << "\t" << profiles[i][j];
    }

    if ((end_conditional_entropy) && (end_partial_entropy)) {
      os << "\t" << begin_conditional_entropy[i] << "\t" << end_conditional_entropy[i] << "\t" << marginal_entropy[i]
         << "\t" << begin_partial_entropy[i] << "\t" << end_partial_entropy[i];
    }
    else {
      os << "\t" << begin_conditional_entropy[i] << "\t" << marginal_entropy[i]
         << "\t" << begin_partial_entropy[i];
    }

    if (i == 0) {
      os << "\t" << identifier[index];
    }
    os << endl;
  }

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture des profils d'etats et des profils d'entropie pour une sequence
 *  au format Gnuplot.
 *
 *  arguments : stream, indice de la sequence, nombre d'etats,
 *              pointeur sur les profils d'etat,
 *              pointeurs sur les profils d'entropies.
 *
 *--------------------------------------------------------------*/

ostream& Sequences::profile_plot_print(ostream &os , int index , int nb_state ,
                                       double **profiles , double *begin_conditional_entropy ,
                                       double *marginal_entropy , double *begin_partial_entropy ,
                                       double *end_conditional_entropy , double *end_partial_entropy) const

{
  register int i , j;


  for (i = 0;i < length[index];i++) {
    for (j = 0;j < nb_state;j++) {
      os << profiles[i][j] << " ";
    }

    if ((end_conditional_entropy) && (end_partial_entropy)) {
      os << begin_conditional_entropy[i] << " " << end_conditional_entropy[i] << " "
         << marginal_entropy[i] << " " << begin_partial_entropy[i] << " "
         << end_partial_entropy[i] << endl;
    }
    else {
      os << begin_conditional_entropy[i] << " " << marginal_entropy[i] << " "
         << begin_partial_entropy[i] << endl;
    }
  }

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Calcul des probabilites a posteriori des etats par l'algorithme forward-backward.
 *
 *  arguments : reference sur un objet Variable_order_markov_data, indice de la sequence,
 *              stream, format de fichier ('a' : ASCII, 's' : Spreadsheet, 'g' : Gnuplot),
 *              references sur l'entropie marginale maximum et l'entropie (pour la visualisation).
 *
 *--------------------------------------------------------------*/

double Hidden_variable_order_markov::forward_backward(const Variable_order_markov_data &seq ,
                                                      int index , ostream &os , char format ,
                                                      double &max_marginal_entropy , double &entropy1) const

{
  register int i , j , k;
  int *pstate , **poutput;
  double seq_likelihood , state_seq_likelihood , **forward , norm , **predicted , entropy2 ,
         buff , **backward , *auxiliary , backward_max , **state_backward , *transition_predicted ,
         **forward_state_entropy , **transition_entropy , *begin_partial_entropy , *begin_conditional_entropy ,
         *end_backward , *end_partial_entropy , *end_conditional_entropy , *marginal_entropy;

// double **backward_state_entropy;


  // initialisations

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

  transition_predicted = new double[nb_row];

  forward_state_entropy = new double*[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    forward_state_entropy[i] = new double[nb_row];
  }

  transition_entropy = new double*[nb_row];
  for (i = 1;i < nb_row;i++) {
    transition_entropy[i] = new double[nb_state];
  }

  begin_partial_entropy = new double[seq.length[index]];
  begin_conditional_entropy = new double[seq.length[index]];

/*  backward_state_entropy = new double*[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    backward_state_entropy[i] = new double[nb_row];
  } */

  end_backward = new double[nb_row];
  end_partial_entropy = new double[seq.length[index]];
  end_conditional_entropy = new double[seq.length[index]];

  marginal_entropy = new double[seq.length[index]];

  poutput = new int*[nb_output_process];

  for (i = 0;i < nb_output_process;i++) {
    poutput[i] = seq.sequence[index][i + 1];
  }

  // recurrence "forward"

  seq_likelihood = 0.;
  norm = 0.;

  switch (type) {

  case 'o' : {
    for (i = 1;i < nb_row;i++) {
      if (order[i] == 1) {
        forward[0][i] = initial[state[i][0]];

        for (j = 0;j < nb_output_process;j++) {
          if (nonparametric_process[j + 1]) {
            forward[0][i] *= nonparametric_process[j + 1]->observation[state[i][0]]->mass[*poutput[j]];
          }
          else {
            forward[0][i] *= parametric_process[j + 1]->observation[state[i][0]]->mass[*poutput[j]];
          }
        }

        norm += forward[0][i];
      }

      else {
        forward[0][i] = 0.;
      }
    }
    break;
  }

  case 'e' : {
    for (i = 1;i < nb_row;i++) {
      if (!child[i]) {
        forward[0][i] = initial[i];

        for (j = 0;j < nb_output_process;j++) {
          if (nonparametric_process[j + 1]) {
            forward[0][i] *= nonparametric_process[j + 1]->observation[state[i][0]]->mass[*poutput[j]];
          }
          else {
            forward[0][i] *= parametric_process[j + 1]->observation[state[i][0]]->mass[*poutput[j]];
          }
        }

        norm += forward[0][i];
      }

      else {
        forward[0][i] = 0.;
      }
    }
    break;
  }
  }

  if (norm > 0.) {
    for (i = 1;i < nb_row;i++) {
      forward[0][i] /= norm;
    }

    seq_likelihood += log(norm);
  }

  else {
    seq_likelihood = D_INF;
  }

  if (seq_likelihood != D_INF) {
    for (i = 1;i < nb_row;i++) {
      forward_state_entropy[0][i] = 0.;
    }

    for (i = 1;i < seq.length[index];i++) {
      for (j = 0;j < nb_output_process;j++) {
        poutput[j]++;
      }
      norm = 0.;

      for (j = 1;j < nb_row;j++) {
        forward[i][j] = 0.;
        for (k = 0;k < nb_memory[j];k++) {
          transition_predicted[k] = transition[previous[j][k]][state[j][0]] * forward[i - 1][previous[j][k]];
          forward[i][j] += transition_predicted[k];

//          forward[i][j] += transition[previous[j][k]][state[j][0]] * forward[i - 1][previous[j][k]];
        }
        predicted[i][j] = forward[i][j];

        forward_state_entropy[i][j] = 0.;
        if (predicted[i][j] > 0.) {
          for (k = 0;k < nb_memory[j];k++) {
            if (transition_predicted[k] > 0.) {
              buff = transition_predicted[k] / predicted[i][j];
              forward_state_entropy[i][j] += buff * (forward_state_entropy[i - 1][previous[j][k]] - log(buff));
            }
          }

          if (forward_state_entropy[i][j] < 0.) {
            forward_state_entropy[i][j] = 0.;
          }
        }

        for (k = 0;k < nb_output_process;k++) {
          if (nonparametric_process[k + 1]) {
            forward[i][j] *= nonparametric_process[k + 1]->observation[state[j][0]]->mass[*poutput[k]];
          }
          else {
            forward[i][j] *= parametric_process[k + 1]->observation[state[j][0]]->mass[*poutput[k]];
          }
        }

        norm += forward[i][j];
      }

      if (norm > 0.) {
        for (j = 1;j < nb_row;j++) {
          forward[i][j] /= norm;
        }
        seq_likelihood += log(norm);
      }

      else {
        seq_likelihood = D_INF;
        break;
      }
    }

    entropy1 = 0.;
    i = seq.length[index] - 1;
    for (j = 1;j < nb_row;j++) {
      if (forward[i][j] > 0.) {
        entropy1 += forward[i][j] * (forward_state_entropy[i][j] - log(forward[i][j]));
      }
    }
  }

  // recurrence "backward"

  if (seq_likelihood != D_INF) {
    entropy2 = 0.;

    for (i = 1;i < nb_row;i++) {
      for (j = 0;j < nb_state;j++) {
        transition_entropy[i][j] = 0.;
      }
    }

    i = seq.length[index] - 1;
    for (j = 1;j < nb_row;j++) {
      backward[i][j] = forward[i][j];

      if (backward[i][j] > 0.) {
        for (k = 0;k < nb_output_process;k++) {
          if (nonparametric_process[k + 1]) {
            if (nonparametric_process[k + 1]->observation[state[j][0]]->mass[*poutput[k]] > 0.) {
              entropy2 -= backward[i][j] * log(nonparametric_process[k + 1]->observation[state[j][0]]->mass[*poutput[k]]);
            }
          }
          else {
            if (parametric_process[k + 1]->observation[state[j][0]]->mass[*poutput[k]] > 0.) {
              entropy2 -= backward[i][j] * log(parametric_process[k + 1]->observation[state[j][0]]->mass[*poutput[k]]);
            }
          }
        }
      }

//      backward_state_entropy[i][j] = 0.;
    }

    for (i = seq.length[index] - 2;i >= 0;i--) {
      for (j = 0;j < nb_output_process;j++) {
        poutput[j]--;
      }

      for (j = 1;j < nb_row;j++) {
        if (predicted[i + 1][j] > 0.) {
          auxiliary[j] = backward[i + 1][j] / predicted[i + 1][j];
        }
        else {
          auxiliary[j] = 0.;
        }
      }

      for (j = 1;j < nb_row;j++) {
        backward[i][j] = 0.;
//        backward_state_entropy[i][j] = 0.;

        if (next[j]) {
//          norm = 0.;

          for (k = 0;k < nb_state;k++) {
/*            transition_predicted[k] = auxiliary[next[j][k]] * transition[j][k];
            norm += transition_predicted[k]; */

            buff = auxiliary[next[j][k]] * transition[j][k] * forward[i][j];
            backward[i][j] += buff;
            transition_entropy[j][k] += buff;

/*            if (transition[j][k] > 0.) {
              entropy2 -= buff * log(transition[j][k]);
            } */
          }

/*          if (norm > 0.) {
            for (k = 0;k < nb_state;k++) {
              if (transition_predicted[k] > 0.) {
                buff = transition_predicted[k] / norm;
                backward_state_entropy[i][j] += buff * (backward_state_entropy[i + 1][next[j][k]] - log(buff));
              }
            }

            if (backward_state_entropy[i][j] < 0.) {
              backward_state_entropy[i][j] = 0.;
            }
          } */

          if (backward[i][j] > 0.) {
            for (k = 0;k < nb_output_process;k++) {
              if (nonparametric_process[k + 1]) {
                if (nonparametric_process[k + 1]->observation[state[j][0]]->mass[*poutput[k]] > 0.) {
                  entropy2 -= backward[i][j] * log(nonparametric_process[k + 1]->observation[state[j][0]]->mass[*poutput[k]]);
                }
              }
              else {
                if (parametric_process[k + 1]->observation[state[j][0]]->mass[*poutput[k]] > 0.) {
                  entropy2 -= backward[i][j] * log(parametric_process[k + 1]->observation[state[j][0]]->mass[*poutput[k]]);
                }
              }
            }
          }
        }
      }
    }

    for (i = 1;i < nb_row;i++) {
      switch (type) {

      case 'o' : {
        if ((order[i] == 1) && (initial[state[i][0]] > 0.)) {
          entropy2 -= backward[0][i] * log(initial[state[i][0]]);
        }
        break;
      }

      case 'e' : {
        if ((!child[i]) && (initial[i] > 0.)) {
          entropy2 -= backward[0][i] * log(initial[i]);
        }
        break;
      }
      }
    }

    for (i = 1;i < nb_row;i++) {
      for (j = 0;j < nb_state;j++) {
        if (transition[i][j] > 0.) {
          entropy2 -= transition_entropy[i][j] * log(transition[i][j]);
        }
      }
    }

    entropy2 += seq_likelihood;

#   ifdef MESSAGE
    if ((entropy2 < entropy1 - DOUBLE_ERROR) || (entropy2 > entropy1 + DOUBLE_ERROR)) {
      cout << "\nERROR: " << entropy1 << " " << entropy2 << endl;
    }
#   endif

    for (i = 0;i < seq.length[index];i++) {
      begin_partial_entropy[i] = 0.;
      for (j = 1;j < nb_row;j++) {
        if (backward[i][j] > 0.) {
          begin_partial_entropy[i] += backward[i][j] * (forward_state_entropy[i][j] - log(backward[i][j]));
        }
      }
      if (begin_partial_entropy[i] < 0.) {
        begin_partial_entropy[i] = 0.;
      }
    }

    begin_conditional_entropy[0] = 0.;
    for (i = 1;i < nb_row;i++) {
      if (backward[0][i] > 0.) {
        begin_conditional_entropy[0] -= backward[0][i] * log(backward[0][i]);
      }
    }
    if (begin_conditional_entropy[0] < 0.) {
      begin_conditional_entropy[0] = 0.;
    }

    for (i = 1;i < seq.length[index];i++) {
      begin_conditional_entropy[i] = 0.;
      for (j = 1;j < nb_row;j++) {
        for (k = 0;k < nb_memory[j];k++) {
          if ((predicted[i][j] > 0.) && (backward[i - 1][previous[j][k]] > 0.)) {
            buff = backward[i][j] * transition[previous[j][k]][state[j][0]] *
                   forward[i - 1][previous[j][k]] / predicted[i][j];
            if (buff > 0.) {
              begin_conditional_entropy[i] -= buff * log(buff / backward[i - 1][previous[j][k]]);
            }
          }
        }
      }
      if (begin_conditional_entropy[i] < 0.) {
        begin_conditional_entropy[i] = 0.;
      }
    }

/*    for (i = 0;i < seq.length[index];i++) {
      end_partial_entropy[i] = begin_partial_entropy[seq.length[index] - 1];
      for (j = 1;j < nb_row;j++) {
        if (backward[i][j] > 0.) {
          end_partial_entropy[i - order[j] + 1] -= backward[i][j] * forward_state_entropy[i][j];
        }
      }
      if (end_partial_entropy[i] < 0.) {
        end_partial_entropy[i] = 0.;
      }
    }

    for (i = 0;i < seq.length[index];i++) {
      end_partial_entropy[i] = 0.;
    }

    for (i = 0;i < seq.length[index] - 1;i++) {
      for (j = 1;j < nb_row;j++) {
        if (backward[i][j] > 0.) {
          end_partial_entropy[i - order[j] + 1] += backward[i][j] * (backward_state_entropy[i][j] - log(backward[i][j]));
        }
      }
    }

    i = seq.length[index] - 1;
    for (j = 1;j < nb_row;j++) {
      if (end_backward[j] > 0.) {
        end_partial_entropy[i - order[j] + 1] -= end_backward[j] * log(end_backward[j]);
      }
    }

    for (i = 0;i < seq.length[index];i++) {
      if (end_partial_entropy[i] < 0.) {
        end_partial_entropy[i] = 0.;
      }
    } */

    for (i = 0;i < seq.length[index];i++) {
      end_conditional_entropy[i] = 0.;
    }

    for (i = 0;i < seq.length[index] - 1;i++) {
      for (j = 1;j < nb_row;j++) {
        if (next[j]) {
          for (k = 0;k < nb_state;k++) {
            if (predicted[i + 1][next[j][k]] > 0.) {
              buff = transition[j][k] * forward[i][j] / predicted[i + 1][next[j][k]];
              if (buff > 0.) {
                end_conditional_entropy[i - order[j] + 1] -= (backward[i + 1][next[j][k]] * buff) * log(buff);
              }
            }
          }
        }
      }
    }

    i = seq.length[index] - 1;
    end_backward[0] = 0.;
    for (j = 1;j < nb_row;j++) {
      end_backward[j] = backward[i][j];
    }
    for (j = nb_row - 1;j >= 1;j--) {
      end_backward[parent[j]] += end_backward[j];
    }

#   ifdef DEBUG
    cout << "\nTEST sum to 1: " << end_backward[0] << endl;
#   endif

    for (j = 1;j < nb_row;j++) {
      if ((end_backward[j] > 0.) && (end_backward[parent[j]] > 0.)) {
        end_conditional_entropy[i - order[j] + 1] -= end_backward[j] * log(end_backward[j] / end_backward[parent[j]]);
      }
    }

    for (i = 0;i < seq.length[index];i++) {
      if (end_conditional_entropy[i] < 0.) {
        end_conditional_entropy[i] = 0.;
      }
    }

#   ifdef MESSAGE
    buff = begin_conditional_entropy[0];
    if ((buff < begin_partial_entropy[0] - DOUBLE_ERROR) || (buff > begin_partial_entropy[0] + DOUBLE_ERROR)) {
      cout << "\nERROR: " << 0 << " | " << buff << " " << begin_partial_entropy[0] << endl;
    }
    for (i = 1;i < seq.length[index];i++) {
      buff += begin_conditional_entropy[i];
      if ((buff < begin_partial_entropy[i] - DOUBLE_ERROR) || (buff > begin_partial_entropy[i] + DOUBLE_ERROR)) {
        cout << "\nERROR: " << i << " | " << buff << " " << begin_partial_entropy[i] << endl;
      }
    }

/*    i = seq.length[index] - 1;
    buff = end_conditional_entropy[i];
    if ((buff < end_partial_entropy[i] - DOUBLE_ERROR) || (buff > end_partial_entropy[i] + DOUBLE_ERROR)) {
      cout << "\nERROR: " << i << " | " << buff << " " << end_partial_entropy[i] << endl;
    }
    for (i = seq.length[index] - 2;i >= 0;i--) {
      buff += end_conditional_entropy[i];
      if ((buff < end_partial_entropy[i] - DOUBLE_ERROR) || (buff > end_partial_entropy[i] + DOUBLE_ERROR)) {
        cout << "\nERROR: " << i << " | " << buff << " " << end_partial_entropy[i] << endl;
      }
    } */
#   endif

    // restauration

    pstate = seq.sequence[index][0];

    for (i = 0;i < seq.length[index];i++) {
      for (j = 0;j < nb_state;j++) {
        state_backward[i][j] = 0.;
      }
      for (j = 1;j < nb_row;j++) {
        state_backward[i][state[j][0]] += backward[i][j];
      }
    }

    for (i = 0;i < seq.length[index];i++) {
      backward_max = 0.;
      for (j = 0;j < nb_state;j++) {
        if (state_backward[i][j] > backward_max) {
          backward_max = state_backward[i][j];
          *pstate = j;
        }
      }

      pstate++;
    }

    state_seq_likelihood = Variable_order_markov::likelihood_computation(seq , index);

/*    begin_conditional_entropy[0] = begin_partial_entropy[0];
    for (i = 1;i < seq.length[index];i++) {
      begin_conditional_entropy[i] = begin_partial_entropy[i] - begin_partial_entropy[i - 1];
    } */

    begin_partial_entropy[0] = begin_conditional_entropy[0];
    for (i = 1;i < seq.length[index];i++) {
      begin_partial_entropy[i] = begin_partial_entropy[i - 1] + begin_conditional_entropy[i];
    }

    end_partial_entropy[seq.length[index] - 1] = end_conditional_entropy[seq.length[index] - 1];
    for (i = seq.length[index] - 2;i >= 0;i--) {
      end_partial_entropy[i] = end_partial_entropy[i + 1] + end_conditional_entropy[i];
    }

    max_marginal_entropy = 0.;
    for (i = 0;i < seq.length[index];i++) {
      marginal_entropy[i] = 0.;
      for (j = 0;j < nb_state;j++) {
        if (state_backward[i][j] > 0.) {
          marginal_entropy[i] -= state_backward[i][j] * log(state_backward[i][j]);
        }
      }
      if (marginal_entropy[i] > max_marginal_entropy) {
        max_marginal_entropy = marginal_entropy[i];
      }
      if (marginal_entropy[i] < 0.) {
        marginal_entropy[i] = 0.;
      }
    }

    switch (format) {

    case 'a' : {
      os << "\n" << SEQ_label[SEQL_POSTERIOR_STATE_PROBABILITY] << "\n\n";
//      seq.profile_ascii_print(os , index , nb_state , state_backward ,
//                              STAT_label[STATL_STATE]);
      seq.profile_ascii_print(os , index , nb_state , state_backward , begin_conditional_entropy ,
                              marginal_entropy , begin_partial_entropy , end_conditional_entropy ,
                              end_partial_entropy);

      os << "\n" << STAT_label[STATL_LIKELIHOOD] << ": " << seq_likelihood
         << "\n" << SEQ_label[SEQL_STATE_SEQUENCE_LIKELIHOOD] << ": " << state_seq_likelihood
         << "   (" << exp(state_seq_likelihood - seq_likelihood) << ")" << endl;
      break;
    }

    case 's' : {
      os << "\n" << SEQ_label[SEQL_POSTERIOR_STATE_PROBABILITY] << "\n\n";
//      seq.profile_spreadsheet_print(os , index , nb_state , state_backward ,
//                                    STAT_label[STATL_STATE]);
      seq.profile_spreadsheet_print(os , index , nb_state , state_backward , begin_conditional_entropy ,
                                    marginal_entropy , begin_partial_entropy , end_conditional_entropy ,
                                    end_partial_entropy);

      os << "\n" << STAT_label[STATL_LIKELIHOOD] << "\t" << seq_likelihood
         << "\n" << SEQ_label[SEQL_STATE_SEQUENCE_LIKELIHOOD] << "\t" << state_seq_likelihood
         << "\t" << exp(state_seq_likelihood - seq_likelihood) << endl;
      break;
    }

    case 'g' : {
//      seq.profile_plot_print(os , index , nb_state , state_backward);
      seq.profile_plot_print(os , index , nb_state , state_backward , begin_conditional_entropy ,
                             marginal_entropy , begin_partial_entropy , end_conditional_entropy ,
                             end_partial_entropy);
      break;
    }
    }

    if (format != 'g') {
/*      double gini_index;

      gini_index = 0.;
      for (i = 0;i < seq.length[index];i++) {
        for (j = 0;j < nb_state;j++) {
          gini_index += state_backward[i][j] * (1. - state_backward[i][j]);
        }
      } */

      double entropy3 , observation , nb_state_sequence;

      entropy3 = 0.;
      for (i = 0;i < seq.length[index];i++) {
        for (j = 0;j < nb_state;j++) {
          if (state_backward[i][j] > 0.) {
            entropy3 -= state_backward[i][j] * log(state_backward[i][j]);
          }
        }
      }

      // calcul du nombre de sequences d'etats possibles

      for (i = 0;i < nb_output_process;i++) {
        poutput[i] = seq.sequence[index][i + 1];
      }

      // recurrence "forward"

      switch (type) {

      case 'o' : {
        for (i = 1;i < nb_row;i++) {
          if (order[i] == 1) {
            forward[0][i] = initial[state[i][0]];

            for (j = 0;j < nb_output_process;j++) {
              if (nonparametric_process[j + 1]) {
                forward[0][i] *= nonparametric_process[j + 1]->observation[state[i][0]]->mass[*poutput[j]];
              }
              else {
                forward[0][i] *= parametric_process[j + 1]->observation[state[i][0]]->mass[*poutput[j]];
              }
            }

            if (forward[0][i] > 0.) {
              forward[0][i] = 1.;
            }
          }

          else {
            forward[0][i] = 0.;
          }
        }
        break;
      }

      case 'e' : {
        for (i = 1;i < nb_row;i++) {
          if (!child[i]) {
            forward[0][i] = initial[i];

            for (j = 0;j < nb_output_process;j++) {
              if (nonparametric_process[j + 1]) {
                forward[0][i] *= nonparametric_process[j + 1]->observation[state[i][0]]->mass[*poutput[j]];
              }
              else {
                forward[0][i] *= parametric_process[j + 1]->observation[state[i][0]]->mass[*poutput[j]];
              }
            }

            if (forward[0][i] > 0.) {
              forward[0][i] = 1.;
            }
          }

          else {
            forward[0][i] = 0.;
          }
        }
        break;
      }
      }

      for (i = 1;i < seq.length[index];i++) {
        for (j = 0;j < nb_output_process;j++) {
          poutput[j]++;
        }

        for (j = 1;j < nb_row;j++) {
          forward[i][j] = 0.;
          for (k = 0;k < nb_memory[j];k++) {
            if (transition[previous[j][k]][state[j][0]] > 0.) {
              forward[i][j] += forward[i - 1][previous[j][k]];
            }
          }

          observation = 1.;
          for (k = 0;k < nb_output_process;k++) {
            if (nonparametric_process[k + 1]) {
              observation *= nonparametric_process[k + 1]->observation[state[j][0]]->mass[*poutput[k]];
            }
            else {
              observation *= parametric_process[k + 1]->observation[state[j][0]]->mass[*poutput[k]];
            }
          }

          if (observation == 0.) {
            forward[i][j] = 0.;
          }
        }
      }

      nb_state_sequence = 0.;
      i = seq.length[index] - 1;
      for (j = 1;j < nb_row;j++) {
        nb_state_sequence += forward[i][j];
      }

      switch (format) {
      case 'a' :
/*        os << "\n" << SEQ_label[SEQL_GINI_INDEX] << ": " << gini_index << " ("
           << gini_index / seq.length[index] */
        os << "\n" << SEQ_label[SEQL_STATE_SEQUENCE_ENTROPY] << ": " << entropy1
           << " (" << entropy1 / seq.length[index] << ")   " << SEQ_label[SEQL_UPPER_BOUND] << ": "
           << log((double)nb_state_sequence) << " ("
           << log((double)nb_state_sequence) / seq.length[index]
           << ")\n" << SEQ_label[SEQL_MARGINAL_ENTROPY_SUM] << ": " << entropy3 << " ("
           << entropy3 / seq.length[index] << ")\n\n"
           << SEQ_label[SEQL_NB_STATE_SEQUENCE] << ": " << nb_state_sequence << endl;
        break;
      case 's' :
/*        os << "\n" << SEQ_label[SEQL_GINI_INDEX] << "\t" << gini_index << "\t"
           << gini_index / seq.length[index] */
        os << "\n" << SEQ_label[SEQL_STATE_SEQUENCE_ENTROPY] << "\t" << entropy1
           << "\t" << entropy1 / seq.length[index] << "\t" << SEQ_label[SEQL_UPPER_BOUND] << "\t"
           << log((double)nb_state_sequence) << "\t"
           << log((double)nb_state_sequence) / seq.length[index]
           << "\n" << SEQ_label[SEQL_MARGINAL_ENTROPY_SUM] << "\t" << entropy3 << "\t"
           << entropy3 / seq.length[index] << "\n\n"
           << SEQ_label[SEQL_NB_STATE_SEQUENCE] << "\t" << nb_state_sequence << endl;
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

  delete [] transition_predicted;

  for (i = 0;i < seq.length[index];i++) {
    delete [] forward_state_entropy[i];
  }
  delete [] forward_state_entropy;

  for (i = 1;i < nb_row;i++) {
    delete [] transition_entropy[i];
  }
  delete [] transition_entropy;

  delete [] begin_partial_entropy;
  delete [] begin_conditional_entropy;

/*  for (i = 0;i < seq.length[index];i++) {
    delete [] backward_state_entropy[i];
  }
  delete [] backward_state_entropy; */

  delete [] end_backward;
  delete [] end_partial_entropy;
  delete [] end_conditional_entropy;

  delete [] marginal_entropy;

  delete [] poutput;

  return (seq_likelihood);
}


/*--------------------------------------------------------------*
 *
 *  Simulation de L sequences d'etats correspondant a une sequence observee
 *  par l'algorithme forward-backward.
 *
 *  arguments : reference sur un objet Variable_order_markov_data, indice de la sequence,
 *              stream, format de fichier ('a' : ASCII, 's' : Spreadsheet),
 *              nombre de sequences d'etats.
 *
 *--------------------------------------------------------------*/

double Hidden_variable_order_markov::forward_backward_sampling(const Variable_order_markov_data &seq ,
                                                               int index , ostream &os ,
                                                               char format , int nb_state_sequence) const

{
  register int i , j , k;
  int memory , *pstate , **poutput;
  double seq_likelihood , state_seq_likelihood , **forward , norm , **predicted ,
         *backward , *cumul_backward;

# ifdef DEBUG
  double sum;
# endif


  // initialisations

  forward = new double*[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    forward[i] = new double[nb_row];
  }

  predicted = new double*[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    predicted[i] = new double[nb_row];
  }

  backward = new double[nb_row];
  cumul_backward = new double[nb_row];

  poutput = new int*[nb_output_process];

  for (i = 0;i < nb_output_process;i++) {
    poutput[i] = seq.sequence[index][i + 1];
  }

# ifdef DEBUG
  double **state_sequence_probability;


  state_sequence_probability = new double*[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    state_sequence_probability[i] = new double[nb_state];
    for (j = 0;j < nb_state;j++) {
      state_sequence_probability[i][j] = 0.;
    }
  }
# endif

  // recurrence "forward"

  seq_likelihood = 0.;
  norm = 0.;

  switch (type) {

  case 'o' : {
    for (i = 1;i < nb_row;i++) {
      if (order[i] == 1) {
        forward[0][i] = initial[state[i][0]];

        for (j = 0;j < nb_output_process;j++) {
          if (nonparametric_process[j + 1]) {
            forward[0][i] *= nonparametric_process[j + 1]->observation[state[i][0]]->mass[*poutput[j]];
          }
          else {
            forward[0][i] *= parametric_process[j + 1]->observation[state[i][0]]->mass[*poutput[j]];
          }
        }

        norm += forward[0][i];
      }

      else {
        forward[0][i] = 0.;
      }
    }
    break;
  }

  case 'e' : {
    for (i = 1;i < nb_row;i++) {
      if (!child[i]) {
        forward[0][i] = initial[i];

        for (j = 0;j < nb_output_process;j++) {
          if (nonparametric_process[j + 1]) {
            forward[0][i] *= nonparametric_process[j + 1]->observation[state[i][0]]->mass[*poutput[j]];
          }
          else {
            forward[0][i] *= parametric_process[j + 1]->observation[state[i][0]]->mass[*poutput[j]];
          }
        }

        norm += forward[0][i];
      }

      else {
        forward[0][i] = 0.;
      }
    }
    break;
  }
  }

  if (norm > 0.) {
    for (i = 1;i < nb_row;i++) {
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
      norm = 0.;

      for (j = 1;j < nb_row;j++) {
        forward[i][j] = 0.;
        for (k = 0;k < nb_memory[j];k++) {
          forward[i][j] += transition[previous[j][k]][state[j][0]] * forward[i - 1][previous[j][k]];
        }
        predicted[i][j] = forward[i][j];

        for (k = 0;k < nb_output_process;k++) {
          if (nonparametric_process[k + 1]) {
            forward[i][j] *= nonparametric_process[k + 1]->observation[state[j][0]]->mass[*poutput[k]];
          }
          else {
            forward[i][j] *= parametric_process[k + 1]->observation[state[j][0]]->mass[*poutput[k]];
          }
        }

        norm += forward[i][j];
      }

      if (norm > 0.) {
        for (j = 1;j < nb_row;j++) {
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

  if (seq_likelihood != D_INF) {

#   ifdef MESSAGE
    cout << "\n";
#   endif

    // passes "backward"

    for (i = 0;i < nb_state_sequence;i++) {
      j = seq.length[index] - 1;
      pstate = seq.sequence[index][0] + j;
      ::cumul_computation(nb_row - 1 , forward[j] + 1 , cumul_backward);
      memory = 1 + cumul_method(nb_row - 1 , cumul_backward);
      *pstate = state[memory][0];

      for (j = seq.length[index] - 2;j >= 0;j--) {
        for (k = 0;k < nb_memory[memory];k++) {
          backward[k] = transition[previous[memory][k]][state[memory][0]] *
                        forward[j][previous[memory][k]] / predicted[j + 1][memory];
        }

#       ifdef DEBUG
        sum = 0.;
        for (k = 0;k < nb_memory[memory];k++) {
          sum += backward[k];
        }
        if ((sum < 1. - DOUBLE_ERROR) || (sum > 1. + DOUBLE_ERROR)) {
          cout << "\nERROR: " << j << " " << sum << endl;
        }
#       endif

        ::cumul_computation(nb_memory[memory] , backward , cumul_backward);
        memory = previous[memory][cumul_method(nb_memory[memory] , cumul_backward)];
        *--pstate = state[memory][0];
      }

#     ifdef DEBUG
      pstate = seq.sequence[index][0];
      for (j = 0;j < seq.length[index];j++) {
        state_sequence_probability[j][*pstate++]++;
      }
#     endif

#     ifdef MESSAGE
      state_seq_likelihood = Variable_order_markov::likelihood_computation(seq , index);

      pstate = seq.sequence[index][0];

      switch (format) {

      case 'a' : {
        for (j = 0;j < seq.length[index];j++) {
          os << *pstate++ << " ";
        }

        os << "  " << i + 1 << "  " << state_seq_likelihood
           << "   (" << exp(state_seq_likelihood - seq_likelihood) << ")" << endl;
        break;
      }

      case 's' : {
        for (j = 0;j < seq.length[index];j++) {
          os << *pstate++ << "\t";
        }

        os << "\t" << i + 1 << "\t" << state_seq_likelihood
           << "\t" << exp(state_seq_likelihood - seq_likelihood) << endl;
        break;
      }
      }
#     endif

    }

#   ifdef DEBUG
    if (nb_state_sequence >= 1000) {
      for (i = 0;i < seq.length[index];i++) {
        for (j = 0;j < nb_state;j++) {
          state_sequence_probability[i][j] /= nb_state_sequence;
        }
      }

      pstate = seq.sequence[index][0];
      for (j = 0;j < seq.length[index];j++) {
        *pstate++ = I_DEFAULT;
      }

      os << "\n" << SEQ_label[SEQL_POSTERIOR_STATE_PROBABILITY] << "\n\n";
      seq.profile_ascii_print(os , index , nb_state , state_sequence_probability ,
                              STAT_label[STATL_STATE]);
    }
#   endif

  }

  for (i = 0;i < seq.length[index];i++) {
    delete [] forward[i];
  }
  delete [] forward;

  for (i = 0;i < seq.length[index];i++) {
    delete [] predicted[i];
  }
  delete [] predicted;

  delete [] backward;
  delete [] cumul_backward;

  delete [] poutput;

# ifdef DEBUG
  for (i = 0;i < seq.length[index];i++) {
    delete [] state_sequence_probability[i];
  }
  delete [] state_sequence_probability;
# endif

  return seq_likelihood;
}


/*--------------------------------------------------------------*
 *
 *  Calcul des probabilites a posteriori des etats par l'algorithme forward-backward,
 *  des vraisemblances des sequences d'etats optimales par l'algorithme de Viterbi
 *  forward-backward, des L sequences d'etats optimales par l'algorithme de Viterbi
 *  generalise ou l'algorithme forward-backward de simulation et ecriture des resultats.
 *
 *  arguments : reference sur un objet Format_error, stream, sequences,
 *              identificateur de la sequence, format ('a' : ASCII, 's' : Spreadsheet),
 *              methode de calcul des sequences d'etats (algorithme de Viterbi generalise ou
 *              algorithme forward-backward de simulation), nombre de sequences d'etats.
 *
 *--------------------------------------------------------------*/

bool Hidden_variable_order_markov::state_profile_write(Format_error &error , ostream &os ,
                                                       const Variable_order_markov_data &iseq ,
                                                       int identifier , char format ,
                                                       int state_sequence , int nb_state_sequence) const

{
  bool status = true;
  register int i;
  int index = I_DEFAULT;
  double seq_likelihood , max_marginal_entropy , entropy;
  Hidden_variable_order_markov *hmarkov;
  Variable_order_markov_data *seq;


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

  if (nb_state_sequence < 2) {
    status = false;
    error.update(SEQ_error[SEQR_NB_STATE_SEQUENCE]);
  }

  if (status) {
    if (iseq.type[0] == INT_VALUE) {
      seq = new Variable_order_markov_data((Markovian_sequences&)iseq , 0);
      seq->type[0] = STATE;
    }
    else {
      seq = new Variable_order_markov_data(iseq , false);
    }

    hmarkov = new Hidden_variable_order_markov(*this , false);
    hmarkov->create_cumul();
    hmarkov->log_computation();

    for (i = 0;i < seq->nb_sequence;i++) {
      if ((index == I_DEFAULT) || (index == i)) {
        seq_likelihood = forward_backward(*seq , i , os , format ,
                                          max_marginal_entropy , entropy);
        hmarkov->viterbi_forward_backward(*seq , i , os , format , seq_likelihood);

        switch (state_sequence) {
        case GENERALIZED_VITERBI :
          hmarkov->generalized_viterbi(*seq , i , os , seq_likelihood , format ,
                                       nb_state_sequence);
          break;
        case FORWARD_BACKWARD_SAMPLING :
          forward_backward_sampling(*seq , i , os , format , nb_state_sequence);
          break;
        }
      }
    }

    delete seq;
    delete hmarkov;
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Calcul des probabilites a posteriori des etats par l'algorithme forward-backward,
 *  des vraisemblances des sequences d'etats optimales par l'algorithme de Viterbi
 *  forward-backward, des L sequences d'etats optimales par l'algorithme de Viterbi
 *  generalise ou l'algorithme forward-backward de simulation et ecriture des resultats.
 *
 *  arguments : reference sur un objet Format_error, stream, identificateur de la sequence,
 *              methode de calcul des sequences d'etats (algorithme de Viterbi generalise ou
 *              algorithme forward-backward de simulation), nombre de sequences d'etats.
 *
 *--------------------------------------------------------------*/

bool Hidden_variable_order_markov::state_profile_ascii_write(Format_error &error , ostream &os ,
                                                             int identifier , int state_sequence ,
                                                             int nb_state_sequence) const

{
  bool status;


  error.init();

  if (!markov_data) {
    status = false;
    error.update(STAT_error[STATR_NO_DATA]);
  }
  else {
    status = state_profile_write(error , os , *markov_data , identifier , 'a' ,
                                 state_sequence , nb_state_sequence);
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Calcul des probabilites a posteriori des etats par l'algorithme forward-backward,
 *  des vraisemblances des sequences d'etats optimales par l'algorithme de Viterbi
 *  forward-backward, des L sequences d'etats optimales par l'algorithme de Viterbi
 *  generalise ou l'algorithme forward-backward de simulation et ecriture des resultats
 *  dans un fichier.
 *
 *  arguments : reference sur un objet Format_error, path, identificateur de la sequence,
 *              format de fichier ('a' : ASCII, 's' : Spreadsheet),
 *              methode de calcul des sequences d'etats (algorithme de Viterbi generalise ou
 *              algorithme forward-backward de simulation), nombre de sequences d'etats.
 *
 *--------------------------------------------------------------*/

bool Hidden_variable_order_markov::state_profile_write(Format_error &error , const char *path ,
                                                       int identifier , char format ,
                                                       int state_sequence , int nb_state_sequence) const

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
    status = state_profile_write(error , out_file , *markov_data , identifier ,
                                 format , state_sequence , nb_state_sequence);
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Calcul des probabilites a posteriori des etats par l'algorithme forward-backward,
 *  des vraisemblances des sequences d'etats optimales par l'algorithme de Viterbi
 *  forward-backward et affichage des resultats au format Gnuplot.
 *
 *  arguments : reference sur un objet Format_error, prefixe des fichiers,
 *              sequences, identificateur de la sequence, titre des figures.
 *
 *--------------------------------------------------------------*/

bool Hidden_variable_order_markov::state_profile_plot_write(Format_error &error , const char *prefix ,
                                                            const Variable_order_markov_data &iseq ,
                                                            int identifier , const char *title) const

{
  bool status = true;
  register int i , j;
  int index;
  double seq_likelihood , max_marginal_entropy , entropy , state_seq_likelihood;
  Hidden_variable_order_markov *hmarkov;
  Variable_order_markov_data *seq;
  ostringstream data_file_name[2];
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

    // ecriture des fichiers de donnees

    data_file_name[0] << prefix << 0 << ".dat";
    data_out_file = new ofstream((data_file_name[0].str()).c_str());

    if (!data_out_file) {
      status = false;
      error.update(STAT_error[STATR_FILE_PREFIX]);
    }

    else {
      if (iseq.type[0] == INT_VALUE) {
        seq = new Variable_order_markov_data((Markovian_sequences&)iseq , 0);
        seq->type[0] = STATE;
      }
      else {
        seq = new Variable_order_markov_data(iseq , false);
      }

      seq_likelihood = forward_backward(*seq , index , *data_out_file , 'g' ,
                                        max_marginal_entropy , entropy);
      data_out_file->close();
      delete data_out_file;

      data_file_name[1] << prefix << 1 << ".dat";
      data_out_file = new ofstream((data_file_name[1].str()).c_str());

      hmarkov = new Hidden_variable_order_markov(*this , false);

      hmarkov->create_cumul();
      hmarkov->log_computation();
      state_seq_likelihood = hmarkov->viterbi_forward_backward(*seq , index , *data_out_file ,
                                                               'g' , seq_likelihood);
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
                 << "set title \"";
        if (title) {
          out_file << title << " - ";
        }
        out_file << SEQ_label[SEQL_POSTERIOR_STATE_PROBABILITY] << "\"\n\n";

        if (seq->length[index] - 1 < TIC_THRESHOLD) {
          out_file << "set xtics 0,1" << endl;
        }

        out_file << "plot [0:" << seq->length[index] - 1 << "] [0:1] ";
        for (j = 0;j < nb_state;j++) {
          out_file << "\"" << label((data_file_name[0].str()).c_str()) << "\" using "
                   << j + 1 << " title \"" << STAT_label[STATL_STATE] << " "
                   << j << "\" with linespoints";
          if (j < nb_state - 1) {
            out_file << ",\\";
          }
          out_file << endl;
        }

        if (i == 0) {
          out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
        }
        out_file << endl;

        out_file << "set title \"";
        if (title) {
          out_file << title << " - ";
        }
        out_file << SEQ_label[SEQL_MAX_POSTERIOR_STATE_PROBABILITY] << "\"\n\n";

        out_file << "plot [0:" << seq->length[index] - 1 << "] [0:"
                 << exp(state_seq_likelihood - seq_likelihood) << "] ";
        for (j = 0;j < nb_state;j++) {
          out_file << "\"" << label((data_file_name[1].str()).c_str()) << "\" using "
                   << j + 1 << " title \"" << STAT_label[STATL_STATE] << " "
                   << j << "\" with linespoints";
          if (j < nb_state - 1) {
            out_file << ",\\";
          }
          out_file << endl;
        }

        if (i == 0) {
          out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
        }
        out_file << endl;

        out_file << "set title";
        if (title) {
          out_file << " \"" << title << "\"";
        }
        out_file << "\n\n";

        out_file << "plot [0:" << seq->length[index] - 1 << "] [0:" << max_marginal_entropy << "] "
                 << "\"" << label((data_file_name[0].str()).c_str()) << "\" using "
                 << nb_state + 1 << " title \"" << SEQ_label[SEQL_CONDITIONAL_ENTROPY]
                 << "\" with linespoints,\\" << endl;
        out_file << "\"" << label((data_file_name[0].str()).c_str()) << "\" using "
                 << nb_state + 2 << " title \"" << SEQ_label[SEQL_CONDITIONAL_ENTROPY]
                 << "\" with linespoints,\\" << endl;
        out_file << "\"" << label((data_file_name[0].str()).c_str()) << "\" using "
                 << nb_state + 3 << " title \"" << SEQ_label[SEQL_MARGINAL_ENTROPY]
                 << "\" with linespoints" << endl;

        if (i == 0) {
          out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
        }
        out_file << endl;

        out_file << "set title";
        if (title) {
          out_file << " \"" << title << "\"";
        }
        out_file << "\n\n";

        out_file << "plot [0:" << seq->length[index] - 1 << "] [0:" << entropy << "] "
                 << "\"" << label((data_file_name[0].str()).c_str()) << "\" using "
                 << nb_state + 4 << " title \"" << SEQ_label[SEQL_PARTIAL_STATE_SEQUENCE_ENTROPY]
                 << "\" with linespoints,\\" << endl;
        out_file << "\"" << label((data_file_name[0].str()).c_str()) << "\" using "
                 << nb_state + 5 << " title \"" << SEQ_label[SEQL_PARTIAL_STATE_SEQUENCE_ENTROPY]
                 << "\" with linespoints" << endl;

        if (seq->length[index] - 1 < TIC_THRESHOLD) {
          out_file << "set xtics autofreq" << endl;
        }

        if (i == 1) {
          out_file << "\nset terminal x11" << endl;
        }

        out_file << "\npause 0 \"" << STAT_label[STATL_END] << "\"" << endl;
      }

      delete seq;
      delete hmarkov;
    }
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Calcul des probabilites a posteriori des etats par l'algorithme forward-backward,
 *  des vraisemblances des sequences d'etats optimales par l'algorithme de Viterbi
 *  forward-backward et affichage des resultats au format Gnuplot.
 *
 *  arguments : reference sur un objet Format_error, prefixe des fichiers,
 *              identificateur de la sequence, titre des figures.
 *
 *--------------------------------------------------------------*/

bool Hidden_variable_order_markov::state_profile_plot_write(Format_error &error ,
                                                            const char *prefix , int identifier ,
                                                            const char *title) const

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
 *  Calcul des logarithmes des parametres d'une chaine de Markov
 *  d'ordre variable cachee.
 *
 *--------------------------------------------------------------*/

void Hidden_variable_order_markov::log_computation()

{
  register int i , j;


  Chain::log_computation();

  for (i = 1;i <= nb_output_process;i++) {
    if (nonparametric_process[i]) {
      for (j = 0;j < nb_state;j++) {
        nonparametric_process[i]->observation[j]->log_computation();
      }
    }

    else {
      for (j = 0;j < nb_state;j++) {
        ::log_computation(parametric_process[i]->nb_value , parametric_process[i]->observation[j]->mass ,
                          parametric_process[i]->observation[j]->cumul);
      }
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Calcul des sequences d'etats les plus probables par l'algorithme de Viterbi.
 *
 *  arguments : reference sur un objet Variable_order_markov_data,
 *              pointeur sur les probabilites a posteriori des sequences d'etats
 *              les plus probables, indice de la sequence.
 *
 *--------------------------------------------------------------*/

double Hidden_variable_order_markov::viterbi(const Variable_order_markov_data &seq ,
                                             double *posterior_probability , int index) const

{
  register int i , j , k , m;
  int length , memory , *pstate , **poutput , **optimal_memory;
  double likelihood = 0. , buff , forward_max , *forward , *previous_forward;


  // initialisations

  forward = new double[nb_row];
  previous_forward = new double[nb_row];

  length = (index == I_DEFAULT ? seq.max_length : seq.length[index]);

  optimal_memory = new int*[length];
  for (i = 0;i < length;i++) {
    optimal_memory[i] = new int[nb_row];
  }

  poutput = new int*[nb_output_process];

  for (i = 0;i < seq.nb_sequence;i++) {
    if ((index == I_DEFAULT) || (index == i)) {
      for (j = 0;j < nb_output_process;j++) {
        poutput[j] = seq.sequence[i][j + 1];
      }

      // recurrence "forward"

      switch (type) {

      case 'o' : {
        for (j = 1;j < nb_row;j++) {
          if (order[j] == 1) {
            forward[j] = cumul_initial[state[j][0]];

            if (forward[j] != D_INF) {
              for (k = 0;k < nb_output_process;k++) {
                if (nonparametric_process[k + 1]) {
                  buff = nonparametric_process[k + 1]->observation[state[j][0]]->cumul[*poutput[k]];
                }
                else {
                  buff = parametric_process[k + 1]->observation[state[j][0]]->cumul[*poutput[k]];
                }

                if (buff == D_INF) {
                  forward[j] = D_INF;
                  break;
                }
                else {
                  forward[j] += buff;
                }
              }
            }
          }

          else {
            forward[j] = D_INF;
          }
        }
        break;
      }

      case 'e' : {
        for (j = 1;j < nb_row;j++) {
          if (!child[j]) {
            forward[j] = cumul_initial[j];

            if (forward[j] != D_INF) {
              for (k = 0;k < nb_output_process;k++) {
                if (nonparametric_process[k + 1]) {
                  buff = nonparametric_process[k + 1]->observation[state[j][0]]->cumul[*poutput[k]];
                }
                else {
                  buff = parametric_process[k + 1]->observation[state[j][0]]->cumul[*poutput[k]];
                }

                if (buff == D_INF) {
                  forward[j] = D_INF;
                  break;
                }
                else {
                  forward[j] += buff;
                }
              }
            }
          }

          else {
            forward[j] = D_INF;
          }
        }
        break;
      }
      }

#     ifdef DEBUG
      cout << "\n" << 0 << " : ";
      for (j = 1;j < nb_row;j++) {
        cout << forward[j] << " | ";
      }
      cout << endl;
#     endif

      for (j = 1;j < seq.length[i];j++) {
        for (k = 0;k < nb_output_process;k++) {
          poutput[k]++;
        }

        for (k = 1;k < nb_row;k++) {
          previous_forward[k] = forward[k];
        }

        for (k = 1;k < nb_row;k++) {
          forward[k] = D_INF;
          for (m = 0;m < nb_memory[k];m++) {
            buff = cumul_transition[previous[k][m]][state[k][0]] + previous_forward[previous[k][m]];
            if (buff > forward[k]) {
              forward[k] = buff;
              optimal_memory[j][k] = previous[k][m];
            }
          }

          if (forward[k] != D_INF) {
            for (m = 0;m < nb_output_process;m++) {
              if (nonparametric_process[m + 1]) {
                buff = nonparametric_process[m + 1]->observation[state[k][0]]->cumul[*poutput[m]];
              }
              else {
                buff = parametric_process[m + 1]->observation[state[k][0]]->cumul[*poutput[m]];
              }

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

#       ifdef DEBUG
        cout << j << " : ";
        for (k = 1;k < nb_row;k++) {
          cout << forward[k] << " " << optimal_memory[j][k] << " | ";
        }
        cout << endl;
#       endif

      }

      // extraction de la vraisemblance du chemin optimal

      pstate = seq.sequence[i][0] + seq.length[i] - 1;
      forward_max = D_INF;

      for (j = 1;j < nb_row;j++) {
        if (forward[j] > forward_max) {
          forward_max = forward[j];
          memory = j;
        }
      }

      if (forward_max != D_INF) {
        likelihood += forward_max;
        *pstate = state[memory][0];
        if (posterior_probability) {
          posterior_probability[i] = forward_max;
        }
      }

      else {
        likelihood = D_INF;
        if (posterior_probability) {
          posterior_probability[i] = 0.;
        }
        break;
      }

      // restauration

      for (j = seq.length[i] - 1;j > 0;j--) {
        memory = optimal_memory[j][memory];
        *--pstate = state[memory][0];
      }

#     ifdef DEBUG
      cout << "\n";
      for (j = seq.length[i] - 1;j >= 0;j--) {
        cout << seq.sequence[i][0][j] << " ";
      }
      cout << endl;
#     endif

    }
  }

  delete [] forward;
  delete [] previous_forward;

  for (i = 0;i < length;i++) {
    delete [] optimal_memory[i];
  }
  delete [] optimal_memory;

  delete [] poutput;

  return likelihood;
}


/*--------------------------------------------------------------*
 *
 *  Calcul des L sequences d'etats optimales par l'algorithme de Viterbi generalise.
 *
 *  arguments : reference sur un objet Variable_order_markov_data, indice de la sequence,
 *              stream, vraisemblance des donnees, format de fichier
 *              ('a' : ASCII, 's' : Spreadsheet), nombre de sequences d'etats.
 *
 *--------------------------------------------------------------*/

double Hidden_variable_order_markov::generalized_viterbi(const Variable_order_markov_data &seq , int index ,
                                                         ostream &os , double seq_likelihood ,
                                                         char format , int inb_state_sequence) const

{
  bool **active_cell;
  register int i , j , k , m;
  int nb_state_sequence , memory , brank , previous_rank , nb_cell , *rank , *pstate ,
      **poutput , ***optimal_memory , ***optimal_rank;
  double buff , forward_max , state_seq_likelihood , likelihood_cumul , **forward ,
         **previous_forward;

  // initialisations

  forward = new double*[nb_row];
  forward[0] = 0;
  for (i = 1;i < nb_row;i++) {
    forward[i] = new double[inb_state_sequence];
  }

  previous_forward = new double*[nb_row];
  previous_forward[0] = 0;
  for (i = 1;i < nb_row;i++) {
    previous_forward[i] = new double[inb_state_sequence];
    for (j = 1;j < inb_state_sequence;j++) {
      previous_forward[i][j] = D_INF;
    }
  }

  rank = new int[nb_row];

  optimal_memory = new int**[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    optimal_memory[i] = new int*[nb_row];
    optimal_memory[i][0] = 0;
    for (j = 1;j < nb_row;j++) {
      optimal_memory[i][j] = new int[inb_state_sequence];
    }
  }

  optimal_rank = new int**[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    optimal_rank[i] = new int*[nb_row];
    optimal_rank[i][0] = 0;
    for (j = 1;j < nb_row;j++) {
      optimal_rank[i][j] = new int[inb_state_sequence];
    }
  }

  active_cell = new bool*[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    active_cell[i] = new bool[nb_state];
    for (j = 0;j < nb_state;j++) {
      active_cell[i][j] = false;
    }
  }

  poutput = new int*[nb_output_process];

  for (i = 0;i < nb_output_process;i++) {
    poutput[i] = seq.sequence[index][i + 1];
  }

# ifdef DEBUG
  double entropy = 0. , **state_sequence_probability;


  state_sequence_probability = new double*[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    state_sequence_probability[i] = new double[nb_state];
    for (j = 0;j < nb_state;j++) {
//      state_sequence_probability[i][j] = 0.;
      state_sequence_probability[i][j] = D_INF;
    }
  }
# endif

  // recurrence "forward"

  switch (type) {

  case 'o' : {
    for (i = 1;i < nb_row;i++) {
      if (order[i] == 1) {
        forward[i][0] = cumul_initial[state[i][0]];

        if (forward[i][0] != D_INF) {
          for (j = 0;j < nb_output_process;j++) {
            if (nonparametric_process[j + 1]) {
              buff = nonparametric_process[j + 1]->observation[state[i][0]]->cumul[*poutput[j]];
            }
            else {
              buff = parametric_process[j + 1]->observation[state[i][0]]->cumul[*poutput[j]];
            }

            if (buff == D_INF) {
              forward[i][0] = D_INF;
              break;
            }
            else {
              forward[i][0] += buff;
            }
          }
        }
      }

      else {
        forward[i][0] = D_INF;
      }
    }
    break;
  }

  case 'e' : {
    for (i = 1;i < nb_row;i++) {
      if (!child[i]) {
        forward[i][0] = cumul_initial[i];

        if (forward[i][0] != D_INF) {
          for (j = 0;j < nb_output_process;j++) {
            if (nonparametric_process[j + 1]) {
              buff = nonparametric_process[j + 1]->observation[state[i][0]]->cumul[*poutput[j]];
            }
            else {
              buff = parametric_process[j + 1]->observation[state[i][0]]->cumul[*poutput[j]];
            }

            if (buff == D_INF) {
              forward[i][0] = D_INF;
              break;
            }
            else {
              forward[i][0] += buff;
            }
          }
        }
      }

      else {
        forward[i][0] = D_INF;
      }
    }
    break;
  }
  }

  nb_state_sequence = 1;

  for (i = 1;i < seq.length[index];i++) {
    for (j = 0;j < nb_output_process;j++) {
      poutput[j]++;
    }

    for (j = 1;j < nb_row;j++) {
      for (k = 0;k < nb_state_sequence;k++) {
        previous_forward[j][k] = forward[j][k];
      }
    }

    if (nb_state_sequence < inb_state_sequence) {
      if (nb_state_sequence * nb_state < inb_state_sequence) {
        nb_state_sequence *= nb_state;
      }
      else {
        nb_state_sequence = inb_state_sequence;
      }
    }

    for (j = 1;j < nb_row;j++) {
      for (k = 1;k < nb_row;k++) {
        rank[k] = 0;
      }

      for (k = 0;k < nb_state_sequence;k++) {
        forward[j][k] = D_INF;
        for (m = 0;m < nb_memory[j];m++) {
          buff = cumul_transition[previous[j][m]][state[j][0]] +
                 previous_forward[previous[j][m]][rank[previous[j][m]]];
          if (buff > forward[j][k]) {
            forward[j][k] = buff;
            optimal_memory[i][j][k] = previous[j][m];
            optimal_rank[i][j][k] = rank[previous[j][m]];
          }
        }

        if (forward[j][k] != D_INF) {
          rank[optimal_memory[i][j][k]]++;

          for (m = 0;m < nb_output_process;m++) {
            if (nonparametric_process[m + 1]) {
              buff = nonparametric_process[m + 1]->observation[state[j][0]]->cumul[*poutput[m]];
            }
            else {
              buff = parametric_process[m + 1]->observation[state[j][0]]->cumul[*poutput[m]];
            }

            if (buff == D_INF) {
              forward[j][k] = D_INF;
              break;
            }
            else {
              forward[j][k] += buff;
            }
          }
        }
      }
    }
  }

  // extraction de la vraisemblance du chemin optimal

  for (i = 1;i < nb_row;i++) {
    rank[i] = 0;
  }
  likelihood_cumul = 0.;

  for (i = 0;i < nb_state_sequence;i++) {
    pstate = seq.sequence[index][0] + seq.length[index] - 1;
    forward_max = D_INF;

    for (j = 1;j < nb_row;j++) {
      if (forward[j][rank[j]] > forward_max) {
        forward_max = forward[j][rank[j]];
        memory = j;
      }
    }

    if (i == 0) {
      state_seq_likelihood = forward_max;
    }

    if (forward_max == D_INF) {
      break;
    }

    // restauration

    *pstate = state[memory][0];
    active_cell[seq.length[index] - 1][*pstate] = true;
    brank = rank[memory];
    rank[memory]++;

#   ifdef DEBUG
    cout << "\n" << *pstate << " " << brank << " | ";
#   endif

    for (j = seq.length[index] - 1;j > 0;j--) {
      previous_rank = optimal_rank[j][memory][brank];
      memory = optimal_memory[j][memory][brank];
      *--pstate = state[memory][0];
      active_cell[j - 1][*pstate] = true;
      brank = previous_rank;

#     ifdef DEBUG
      cout << *pstate << " " << brank << " | ";
#     endif
    }

#   ifdef DEBUG
    cout << endl;
#   endif

    likelihood_cumul += exp(forward_max);

#   ifdef DEBUG
    pstate = seq.sequence[index][0];
    for (j = 0;j < seq.length[index];j++) {
/*      state_sequence_probability[j][*pstate++] += exp(forward_max - seq_likelihood); */

      if (forward_max > state_sequence_probability[j][*pstate]) {
        state_sequence_probability[j][*pstate] = forward_max;
      }
      pstate++;
    }
#   endif

    nb_cell = 0;
    for (j = 0;j < seq.length[index];j++) {
      for (k = 0;k < nb_state;k++) {
        if (active_cell[j][k]) {
          nb_cell++;
        }
      }
    }

#   ifdef MESSAGE
    if (i == 0) {
      os << "\n";
    }

    pstate = seq.sequence[index][0];

    switch (format) {

    case 'a' : {
      for (j = 0;j < seq.length[index];j++) {
        os << *pstate++ << " ";
      }

//      os << "  " << i + 1 << "  " << forward_max << "   (" << exp(forward_max - state_seq_likelihood)
      os << "  " << i + 1 << "  " << forward_max << "   (" << exp(forward_max - seq_likelihood)
         << "  " << likelihood_cumul / exp(seq_likelihood) << "  " << nb_cell << ")" << endl;
      break;
    }

    case 's' : {
      for (j = 0;j < seq.length[index];j++) {
        os << *pstate++ << "\t";
      }

//      os << "\t" << i + 1 << "\t" << forward_max << "\t" << exp(forward_max - state_seq_likelihood)
      os << "\t" << i + 1 << "\t" << forward_max << "\t" << exp(forward_max - seq_likelihood)
         << "\t" << likelihood_cumul / exp(seq_likelihood) << "\t" << nb_cell << endl;
      break;
    }
    }
#   endif

#   ifdef DEBUG
    entropy -= exp(forward_max - seq_likelihood) * forward_max;
#   endif

  }

# ifdef DEBUG
  os << "\n" << SEQ_label[SEQL_STATE_SEQUENCE_ENTROPY] << ": " << entropy + seq_likelihood << endl;

  if (likelihood_cumul / exp(seq_likelihood) > 0.8) {
    for (i = 0;i < seq.length[index];i++) {
      for (j = 0;j < nb_state;j++) {
        if (state_sequence_probability[i][j] != D_INF) {
          state_sequence_probability[i][j] = exp(state_sequence_probability[i][j] - seq_likelihood);
        }
        else {
          state_sequence_probability[i][j] = 0.;
        }
      }
    }

    pstate = seq.sequence[index][0];
    for (j = 0;j < seq.length[index];j++) {
      *pstate++ = I_DEFAULT;
    }

//    os << "\n" << SEQ_label[SEQL_POSTERIOR_STATE_PROBABILITY] << "\n\n";
    os << "\n" << SEQ_label[SEQL_MAX_POSTERIOR_STATE_PROBABILITY] << "\n\n";
    seq.profile_ascii_print(os , index , nb_state , state_sequence_probability ,
                            STAT_label[STATL_STATE]);
  }
# endif

  for (i = 1;i < nb_row;i++) {
    delete [] forward[i];
  }
  delete [] forward;

  for (i = 1;i < nb_row;i++) {
    delete [] previous_forward[i];
  }
  delete [] previous_forward;

  delete [] rank;

  for (i = 0;i < seq.length[index];i++) {
    for (j = 1;j < nb_row;j++) {
      delete [] optimal_memory[i][j];
    }
    delete [] optimal_memory[i];
  }
  delete [] optimal_memory;

  for (i = 0;i < seq.length[index];i++) {
    for (j = 1;j < nb_row;j++) {
      delete [] optimal_rank[i][j];
    }
    delete [] optimal_rank[i];
  }
  delete [] optimal_rank;

  for (i = 0;i < seq.length[index];i++) {
    delete [] active_cell[i];
  }
  delete [] active_cell;

  delete [] poutput;

# ifdef DEBUG
  for (i = 0;i < seq.length[index];i++) {
    delete [] state_sequence_probability[i];
  }
  delete [] state_sequence_probability;
# endif

  return state_seq_likelihood;
}


/*--------------------------------------------------------------*
 *
 *  Calcul des vraisemblances des sequences d'etats optimales
 *  par l'algorithme de Viterbi forward-backward.
 *
 *  arguments : reference sur un objet Variable_order_markov_data, indice de la sequence,
 *              stream, format de fichier ('a' : ASCII, 's' : Spreadsheet, 'g' : Gnuplot),
 *              vraisemblance des donnees.
 *
 *--------------------------------------------------------------*/

double Hidden_variable_order_markov::viterbi_forward_backward(const Variable_order_markov_data &seq ,
                                                              int index , ostream &os , char format ,
                                                              double seq_likelihood) const

{
  register int i , j , k;
  int *pstate , **poutput;
  double buff , state_seq_likelihood , backward_max , **forward , **backward ,
         *auxiliary , **state_backward;


  // initialisations

  forward = new double*[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    forward[i] = new double[nb_row];
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

# ifdef MESSAGE
  int memory , *state_sequence , **optimal_memory;

  optimal_memory = new int*[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    optimal_memory[i] = new int[nb_row];
  }

  state_sequence = new int[seq.length[index]];
# endif

  for (i = 0;i < nb_output_process;i++) {
    poutput[i] = seq.sequence[index][i + 1];
  }

  // recurrence "forward"

  switch (type) {

  case 'o' : {
    for (i = 1;i < nb_row;i++) {
      if (order[i] == 1) {
        forward[0][i] = cumul_initial[state[i][0]];

        if (forward[0][i] != D_INF) {
          for (j = 0;j < nb_output_process;j++) {
            if (nonparametric_process[j + 1]) {
              buff = nonparametric_process[j + 1]->observation[state[i][0]]->cumul[*poutput[j]];
            }
            else {
              buff = parametric_process[j + 1]->observation[state[i][0]]->cumul[*poutput[j]];
            }

            if (buff == D_INF) {
              forward[0][i] = D_INF;
              break;
            }
            else {
              forward[0][i] += buff;
            }
          }
        }
      }

      else {
        forward[0][i] = D_INF;
      }
    }
    break;
  }

  case 'e' : {
    for (i = 1;i < nb_row;i++) {
      if (!child[i]) {
        forward[0][i] = cumul_initial[i];

        if (forward[0][i] != D_INF) {
          for (j = 0;j < nb_output_process;j++) {
            if (nonparametric_process[j + 1]) {
              buff = nonparametric_process[j + 1]->observation[state[i][0]]->cumul[*poutput[j]];
            }
            else {
              buff = parametric_process[j + 1]->observation[state[i][0]]->cumul[*poutput[j]];
            }

            if (buff == D_INF) {
              forward[0][i] = D_INF;
              break;
            }
            else {
              forward[0][i] += buff;
            }
          }
        }
      }

      else {
        forward[0][i] = D_INF;
      }
    }
    break;
  }
  }

  for (i = 1;i < seq.length[index];i++) {
    for (j = 0;j < nb_output_process;j++) {
      poutput[j]++;
    }

    for (j = 1;j < nb_row;j++) {
      forward[i][j] = D_INF;
      for (k = 0;k < nb_memory[j];k++) {
        buff = cumul_transition[previous[j][k]][state[j][0]] + forward[i - 1][previous[j][k]];
        if (buff > forward[i][j]) {
          forward[i][j] = buff;

#         ifdef MESSAGE
          optimal_memory[i][j] = previous[j][k];
#         endif
        }
      }

      if (forward[i][j] != D_INF) {
        for (k = 0;k < nb_output_process;k++) {
          if (nonparametric_process[k + 1]) {
            buff = nonparametric_process[k + 1]->observation[state[j][0]]->cumul[*poutput[k]];
          }
          else {
            buff = parametric_process[k + 1]->observation[state[j][0]]->cumul[*poutput[k]];
          }

          if (buff == D_INF) {
            forward[i][j] = D_INF;
            break;
          }
          else {
            forward[i][j] += buff;
          }
        }
      }
    }
  }

  // extraction de la vraisemblance du chemin optimal

# ifdef MESSAGE
  pstate = state_sequence + seq.length[index] - 1;
# endif

  state_seq_likelihood = D_INF;
  i = seq.length[index] - 1;
  for (j = 1;j < nb_row;j++) {
    if (forward[i][j] > state_seq_likelihood) {
      state_seq_likelihood = forward[i][j];

#     ifdef MESSAGE
      memory = j;
#     endif
    }
  }

  if (state_seq_likelihood != D_INF) {

#   ifdef MESSAGE
    *pstate = state[memory][0];
    for (i = seq.length[index] - 1;i > 0;i--) {
      memory = optimal_memory[i][memory];
      *--pstate = state[memory][0];
    }
#   endif

    // recurrence "backward"

    i = seq.length[index] - 1;
    for (j = 1;j < nb_row;j++) {
      backward[i][j] = 0.;
    }

    for (i = seq.length[index] - 2;i >= 0;i--) {
      for (j = 1;j < nb_row;j++) {
        auxiliary[j] = backward[i + 1][j];

        for (k = 0;k < nb_output_process;k++) {
          if (nonparametric_process[k + 1]) {
            auxiliary[j] += nonparametric_process[k + 1]->observation[state[j][0]]->cumul[*poutput[k]];
          }
          else {
            auxiliary[j] += parametric_process[k + 1]->observation[state[j][0]]->cumul[*poutput[k]];
          }
        }
      }

      for (j = 0;j < nb_output_process;j++) {
        poutput[j]--;
      }

      for (j = 1;j < nb_row;j++) {
        backward[i][j] = D_INF;
        if (next[j]) {
          for (k = 0;k < nb_state;k++) {
            buff = auxiliary[next[j][k]] + cumul_transition[j][k];
            if (buff > backward[i][j]) {
              backward[i][j] = buff;
            }
          }
        }
      }
    }

    // restauration

    pstate = seq.sequence[index][0];

    for (i = 0;i < seq.length[index];i++) {
      backward_max = D_INF;
      for (j = 1;j < nb_row;j++) {
        if (backward[i][j] != D_INF) {
          if (forward[i][j] != D_INF) {
            backward[i][j] += forward[i][j];
            if (backward[i][j] > backward_max) {
              backward_max = backward[i][j];
              *pstate = state[j][0];
            }
          }

          else {
            backward[i][j] = D_INF;
          }
        }
      }

#     ifdef MESSAGE
      if (*pstate != state_sequence[i]) {
        cout << "\nERROR: " << i << " | " << *pstate << " " << state_sequence[i] << endl;
      }
#     endif

      pstate++;
    }

    //  normalisation

    for (i = 0;i < seq.length[index];i++) {
      for (j = 0;j < nb_state;j++) {
        state_backward[i][j] = D_INF;
      }
      for (j = 1;j < nb_row;j++) {
        if (backward[i][j] > state_backward[i][state[j][0]]) {
          state_backward[i][state[j][0]] = backward[i][j];
        }
      }

      for (j = 0;j < nb_state;j++) {
        if (state_backward[i][j] != D_INF) {
          state_backward[i][j] = exp(state_backward[i][j] - seq_likelihood);
//          state_backward[i][j] = exp(state_backward[i][j] - state_seq_likelihood);
        }
        else {
          state_backward[i][j] = 0.;
        }
      }
    }

    switch (format) {

    case 'a' : {
      os << "\n" << SEQ_label[SEQL_MAX_POSTERIOR_STATE_PROBABILITY] << "\n\n";
      seq.profile_ascii_print(os , index , nb_state , state_backward ,
                              STAT_label[STATL_STATE]);

      os << "\n" << SEQ_label[SEQL_STATE_SEQUENCE_LIKELIHOOD] << ": " << state_seq_likelihood
         << "   (" << exp(state_seq_likelihood - seq_likelihood) << ")" << endl;
      break;
    }

    case 's' : {
      os << "\n" << SEQ_label[SEQL_MAX_POSTERIOR_STATE_PROBABILITY] << "\n\n";
      seq.profile_spreadsheet_print(os , index , nb_state , state_backward ,
                                    STAT_label[STATL_STATE]);

      os << "\n" << SEQ_label[SEQL_STATE_SEQUENCE_LIKELIHOOD] << "\t" << state_seq_likelihood
         << "\t" << exp(state_seq_likelihood - seq_likelihood) << endl;
      break;
    }

    case 'g' : {
      seq.profile_plot_print(os , index , nb_state , state_backward);
      break;
    }
    }

#   ifdef DEBUG
    if (format != 'g') {
      double ambiguity = 0.;

      pstate = seq.sequence[index][0];
      for (i = 0;i < seq.length[index];i++) {
        for (j = 0;j < nb_state;j++) {
          if (j != *pstate) {
            ambiguity += state_backward[i][j];
          }
        }
        pstate++;
      }
      ambiguity *= exp(seq_likelihood - state_seq_likelihood);

      switch (format) {
      case 'a' :
        os << "\n" << SEQ_label[SEQL_AMBIGUITY] << ": " << ambiguity
           << " (" << ambiguity / seq.length[index] << ")" << endl;
        break;
      case 's' :
        os << "\n" << SEQ_label[SEQL_AMBIGUITY] << "\t" << ambiguity
           << "\t" << ambiguity / seq.length[index] << "\t" << endl;
        break;
      }
    }
#   endif

  }

  for (i = 0;i < seq.length[index];i++) {
    delete [] forward[i];
  }
  delete [] forward;

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

# ifdef MESSAGE
  for (i = 0;i < seq.length[index];i++) {
    delete [] optimal_memory[i];
  }
  delete [] optimal_memory;

  delete [] state_sequence;
# endif

  return state_seq_likelihood;
}


/*--------------------------------------------------------------*
 *
 *  Calcul des sequences d'etats optimales.
 *
 *  arguments : references sur un objet Format_error et sur un objet Markovian_sequences,
 *              flag sur le calcul des caracteristiques.
 *
 *--------------------------------------------------------------*/

Variable_order_markov_data* Hidden_variable_order_markov::state_sequence_computation(Format_error &error ,
                                                                                     const Markovian_sequences &iseq ,
                                                                                     bool characteristic_flag) const

{
  bool status = true;
  register int i;
  int nb_value;
  Hidden_variable_order_markov *hmarkov;
  Variable_order_markov_data *seq;


  seq = 0;
  error.init();

  if (nb_output_process != iseq.nb_variable) {
    status = false;
    error.update(SEQ_error[SEQR_NB_OUTPUT_PROCESS]);
  }

  else {
    for (i = 0;i < nb_output_process;i++) {
      if (nonparametric_process[i + 1]) {
        nb_value = nonparametric_process[i + 1]->nb_value;
      }
      else {
        nb_value = parametric_process[i + 1]->nb_value;
      }

      if (nb_value < iseq.marginal[i]->nb_value) {
        status = false;
        ostringstream error_message;
        error_message << STAT_label[STATL_OUTPUT_PROCESS] << " " << i + 1 << ": "
                      << SEQ_error[SEQR_NB_OUTPUT];
        error.update((error_message.str()).c_str());
      }
    }
  }

  if (status) {
    seq = new Variable_order_markov_data(iseq , 0 , (type == 'e' ? true : false));

    seq->type[0] = STATE;
    seq->markov = new Variable_order_markov(*this , false);

    hmarkov = new Hidden_variable_order_markov(*this , false);

    hmarkov->create_cumul();
    hmarkov->log_computation();

    seq->likelihood = hmarkov->viterbi(*seq);

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

      seq->build_transition_count(*hmarkov);
      seq->build_observation_histogram();

      if (characteristic_flag) {
        seq->markov->characteristic_computation(*seq , true);
      }
    }

    delete hmarkov;
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Comparaison de differentes chaines de Markov d'ordre variable cachees
 *  pour un ensemble de sequences.
 *
 *  arguments : reference sur un objet Format_error, stream, nombre de chaines
 *              de Markov cachees, pointeur sur les chaines de Markov cachees,
 *              type algorithme (forward ou Viterbi), path.
 *
 *--------------------------------------------------------------*/

bool Markovian_sequences::comparison(Format_error &error , ostream &os , int nb_model ,
                                     const Hidden_variable_order_markov **ihmarkov ,
                                     int algorithm , const char *path) const

{
  bool status = true;
  register int i , j;
  int nb_value;
  double **likelihood;
  Hidden_variable_order_markov **hmarkov;
  Variable_order_markov_data *seq;


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
        if (ihmarkov[i]->nonparametric_process[j + 1]) {
          nb_value = ihmarkov[i]->nonparametric_process[j + 1]->nb_value;
        }
        else {
          nb_value = ihmarkov[i]->parametric_process[j + 1]->nb_value;
        }

        if (nb_value < marginal[j]->nb_value) {
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
      hmarkov = new Hidden_variable_order_markov*[nb_model];
      for (i = 0;i < nb_model;i++) {
        hmarkov[i] = new Hidden_variable_order_markov(*(ihmarkov[i]) , false);
        hmarkov[i]->create_cumul();
        hmarkov[i]->log_computation();
      }

      seq = new Variable_order_markov_data(*this , 0);
    }

    // pour chaque sequence, calcul de la vraisemblance (FORWARD) ou de la vraisemblance
    // du chemin optimal (VITERBI) pour chaque modele possible

    for (i = 0;i < nb_sequence;i++) {
      for (j = 0;j < nb_model;j++) {
        switch (algorithm) {
        case FORWARD :
          likelihood[i][j] = ihmarkov[j]->likelihood_computation(*this , 0 , i);
          break;
        case VITERBI :
          likelihood[i][j] = hmarkov[j]->viterbi(*seq , 0 , i);
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
 *  Simulation par une chaine de Markov d'ordre variable cachee.
 *
 *  arguments : reference sur un objet Format_error,
 *              histogramme des longueurs des sequences,
 *              flag sur le calcul des lois de comptage,
 *              flag calcul d'une divergence de Kullback-Leibler.
 *
 *--------------------------------------------------------------*/

Variable_order_markov_data* Hidden_variable_order_markov::simulation(Format_error &error ,
                                                                     const Histogram &hlength ,
                                                                     bool counting_flag ,
                                                                     bool divergence_flag) const

{
  Markovian_sequences *observ_seq;
  Variable_order_markov_data *seq;


  seq = Variable_order_markov::simulation(error , hlength , counting_flag , divergence_flag);

  if ((seq) && (!divergence_flag)) {
    observ_seq = seq->remove_variable_1();
    seq->hidden_likelihood = likelihood_computation(*observ_seq);
    delete observ_seq;
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Simulation par une chaine de Markov d'ordre variable cachee.
 *
 *  arguments : reference sur un objet Format_error,
 *              nombre et longueur des sequences.
 *
 *--------------------------------------------------------------*/

Variable_order_markov_data* Hidden_variable_order_markov::simulation(Format_error &error ,
                                                                     int nb_sequence , int length ,
                                                                     bool counting_flag) const

{
  Markovian_sequences *observ_seq;
  Variable_order_markov_data *seq;


  seq = Variable_order_markov::simulation(error , nb_sequence , length , counting_flag);

  if (seq) {
    observ_seq = seq->remove_variable_1();
    seq->hidden_likelihood = likelihood_computation(*observ_seq);
    delete observ_seq;
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Simulation par une chaine de Markov d'ordre variable cachee.
 *
 *  arguments : reference sur un objet Format_error, nombre de sequences,
 *              reference sur un objet Markovian_sequences.
 *
 *--------------------------------------------------------------*/

Variable_order_markov_data* Hidden_variable_order_markov::simulation(Format_error &error ,
                                                                     int nb_sequence ,
                                                                     const Markovian_sequences &iseq ,
                                                                     bool counting_flag) const

{
  Markovian_sequences *observ_seq;
  Variable_order_markov_data *seq;


  seq = Variable_order_markov::simulation(error , nb_sequence , iseq , counting_flag);

  if (seq) {
    observ_seq = seq->remove_variable_1();
    seq->hidden_likelihood = likelihood_computation(*observ_seq);
    delete observ_seq;
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Comparaison de chaines de Markov d'ordre variable cachees par calcul
 *  de divergences de Kullback-Leibler.
 *
 *  arguments : reference sur un objet Format_error, stream, nombre de chaines
 *              de Markov cachees, pointeur sur les chaines de Markov cachees,
 *              histogramme des longueurs des sequences, path.
 *
 *--------------------------------------------------------------*/

Distance_matrix* Hidden_variable_order_markov::divergence_computation(Format_error &error , ostream &os , int nb_model ,
                                                                      const Hidden_variable_order_markov **ihmarkov ,
                                                                      Histogram **hlength , const char *path) const

{
  bool status = true , lstatus;
  register int i , j , k;
  int cumul_length;
  double ref_likelihood , target_likelihood , **likelihood;
  const Hidden_variable_order_markov **hmarkov;
  Markovian_sequences *seq;
  Variable_order_markov_data *simul_seq;
  Distance_matrix *dist_matrix;
  ofstream *out_file;


  dist_matrix = 0;
  error.init();

  for (i = 0;i < nb_model - 1;i++) {
    if (ihmarkov[i]->type != type) {
      status = false;
      ostringstream error_message;
      error_message << SEQ_label[SEQL_HIDDEN_MARKOV_CHAIN] << " " << i + 2 << ": "
                    << SEQ_error[SEQR_MODEL_TYPE];
      error.update((error_message.str()).c_str());
    }

    if (ihmarkov[i]->nb_output_process != nb_output_process) {
      status = false;
      ostringstream error_message;
      error_message << SEQ_label[SEQL_HIDDEN_MARKOV_CHAIN] << " " << i + 2 << ": "
                    << SEQ_error[SEQR_NB_OUTPUT_PROCESS];
      error.update((error_message.str()).c_str());
    }

    else {
      for (j = 1;j <= nb_output_process;j++) {
        if ((ihmarkov[i]->nonparametric_process[j]) && (nonparametric_process[j]) &&
            (ihmarkov[i]->nonparametric_process[j]->nb_value != nonparametric_process[j]->nb_value)) {
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

    hmarkov = new const Hidden_variable_order_markov*[nb_model];

    hmarkov[0] = this;
    for (i = 1;i < nb_model;i++) {
      hmarkov[i] = ihmarkov[i - 1];
    }

    dist_matrix = new Distance_matrix(nb_model , SEQ_label[SEQL_HIDDEN_MARKOV_CHAIN]);

    for (i = 0;i < nb_model;i++) {

      // simulation d'un echantillon de sequences a partir d'une chaine de Markov cachee

      simul_seq = hmarkov[i]->simulation(error , *hlength[i] , false , true);
      seq = simul_seq->remove_variable_1();

      likelihood = new double*[seq->nb_sequence];
      for (j = 0;j < seq->nb_sequence;j++) {
        likelihood[j] = new double[nb_model];
      }

      ref_likelihood = 0.;
      for (j = 0;j < seq->nb_sequence;j++) {
        likelihood[j][i] = hmarkov[i]->likelihood_computation(*seq , 0 , j);
        ref_likelihood += likelihood[j][i];
      }

      // calcul des vraisemblances de l'echantillon pour chacune des chaines de Markov cachees

      for (j = 0;j < nb_model;j++) {
        if (j != i) {
          target_likelihood = 0.;
          for (k = 0;k < seq->nb_sequence;k++) {
            likelihood[k][j] = hmarkov[j]->likelihood_computation(*seq , 0 , k);
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
 *  Comparaison de chaines de Markov d'ordre variable cachees par calcul
 *  de divergences de Kullback-Leibler.
 *
 *  arguments : reference sur un objet Format_error, stream, nombre de chaines
 *              de Markov cachees, pointeur sur les chaines de Markov cachees,
 *              nombre et longueur des sequences, path.
 *
 *--------------------------------------------------------------*/

Distance_matrix* Hidden_variable_order_markov::divergence_computation(Format_error &error , ostream &os ,
                                                                      int nb_model , const Hidden_variable_order_markov **hmarkov ,
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
 *  Comparaison de chaines de Markov d'ordre variable cachees par calcul
 *  de divergences de Kullback-Leibler.
 *
 *  arguments : reference sur un objet Format_error, stream, nombre de chaines
 *              de Markov cachees, pointeur sur les chaines de Markov cachees,
 *              pointeurs sur des objets Markovian_sequences, path.
 *
 *--------------------------------------------------------------*/

Distance_matrix* Hidden_variable_order_markov::divergence_computation(Format_error &error , ostream &os ,
                                                                      int nb_model , const Hidden_variable_order_markov **hmarkov ,
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
