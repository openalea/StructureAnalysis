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
#include "stat_tool/curves.h"
#include "stat_tool/markovian.h"
#include "stat_tool/stat_label.h"
#include "sequences.h"
#include "semi_markov.h"
#include "hidden_semi_markov.h"
#include "sequence_label.h"

#include "stat_tool/distribution_reestimation.h"   // probleme compilateur C++ Windows

using namespace std;


/* template <typename Type>    - probleme compilateur C++ Windows
extern void reestimation(int nb_value , Type *reestim , double *pmass ,
                         double min_probability , bool null_probability); */

extern double interval_bisection(Reestimation<double> *distribution_reestim ,
                                 Reestimation<double> *length_bias_reestim);
extern void cumul_computation(int nb_value , const double *pmass , double *pcumul);
extern int cumul_method(int nb_value , const double *cumul , double scale = 1.);
extern void log_computation(int nb_value , const double *pmass , double *plog);

extern char* label(const char *file_name);



/*--------------------------------------------------------------*
 *
 *  Calcul de la vraisemblance de sequences pour une semi-chaine de Markov cachee
 *  par l'algorithme forward.
 *
 *  arguments : reference sur un objet Markovian_sequences, indice de la sequence.
 *
 *--------------------------------------------------------------*/

double Hidden_semi_markov::likelihood_computation(const Markovian_sequences &seq , int index) const

{
  register int i , j , k , m;
  int nb_value , length , **poutput;
  double likelihood = 0. , obs_product , **observation , *norm , *state_norm ,
         *forward1 , **state_in;
  Parametric *occupancy;


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

    length = (index == I_DEFAULT ? seq.max_length : seq.length[index]);

    observation = new double*[length];
    for (i = 0;i < length;i++) {
      observation[i] = new double[nb_state];
    }

    norm = new double[length];
    state_norm = new double[nb_state];
    forward1 = new double[nb_state];

    state_in = new double*[length - 1];
    for (i = 0;i < length - 1;i++) {
      state_in[i] = new double[nb_state];
    }

    poutput = new int*[seq.nb_variable];

    for (i = 0;i < seq.nb_sequence;i++) {
      if ((index == I_DEFAULT) || (index == i)) {
        for (j = 0;j < seq.nb_variable;j++) {
          poutput[j] = seq.sequence[i][j];
        }

        for (j = 0;j < seq.length[i];j++) {
          norm[j] = 0.;

          for (k = 0;k < nb_state;k++) {

            // calcul des probabilites d'observation

            observation[j][k] = 1.;
            for (m = 0;m < nb_output_process;m++) {
              if (nonparametric_process[m + 1]) {
                observation[j][k] *= nonparametric_process[m + 1]->observation[k]->mass[*poutput[m]];
              }
              else {
                observation[j][k] *= parametric_process[m + 1]->observation[k]->mass[*poutput[m]];
              }
            }

            switch (state_subtype[k]) {

            // cas etat semi-markovien

            case SEMI_MARKOVIAN : {
              if (j == 0) {
                state_norm[k] = initial[k];
              }
              else {
                state_norm[k] += state_in[j - 1][k] - forward1[k];
              }
              state_norm[k] *= observation[j][k];

              norm[j] += state_norm[k];
              break;
            }

            // cas etat markovien

            case MARKOVIAN : {
              if (j == 0) {
                forward1[k] = initial[k];
              }
              else {
                forward1[k] = state_in[j - 1][k];
              }
              forward1[k] *= observation[j][k];

              norm[j] += forward1[k];
              break;
            }
            }
          }

          if (norm[j] > 0.) {
            for (k = 0;k < nb_state;k++) {
              switch (state_subtype[k]) {
              case SEMI_MARKOVIAN :
                state_norm[k] /= norm[j];
                break;
              case MARKOVIAN :
                forward1[k] /= norm[j];
                break;
              }
            }

            likelihood += log(norm[j]);
          }

          else {
            likelihood = D_INF;
            break;
          }

          for (k = 0;k < nb_state;k++) {

            // cas etat semi-markovien

            if (state_subtype[k] == SEMI_MARKOVIAN) {
              occupancy = nonparametric_process[0]->sojourn_time[k];
              obs_product = 1.;
              forward1[k] = 0.;

              if (j < seq.length[i] - 1) {
                for (m = 1;m <= MIN(j + 1 , occupancy->nb_value - 1);m++) {
                  obs_product *= observation[j - m + 1][k] / norm[j - m + 1];
                  if (obs_product == 0.) {
                    break;
                  }

                  if (m < j + 1) {
                    forward1[k] += obs_product * occupancy->mass[m] * state_in[j - m][k];
                  }

                  else {
                    switch (type) {
                    case 'o' :
                      forward1[k] += obs_product * occupancy->mass[m] * initial[k];
                      break;
                    case 'e' :
                      forward1[k] += obs_product * forward[k]->mass[m] * initial[k];
                      break;
                    }
                  }
                }
              }

              else {
                for (m = 1;m <= MIN(j + 1 , occupancy->nb_value - 1);m++) {
                  obs_product *= observation[j - m + 1][k] / norm[j - m + 1];
                  if (obs_product == 0.) {
                    break;
                  }

                  if (m < j + 1) {
                    forward1[k] += obs_product * (1. - occupancy->cumul[m - 1]) * state_in[j - m][k];
                  }

                  else {
                    switch (type) {
                    case 'o' :
                      forward1[k] += obs_product * (1. - occupancy->cumul[m - 1]) * initial[k];
                      break;
                    case 'e' :
                      forward1[k] += obs_product * (1. - forward[k]->cumul[m - 1]) * initial[k];
                      break;
                    }
                  }
                }
              }
            }
          }

          if (j < seq.length[i] - 1) {
            for (k = 0;k < nb_state;k++) {
              state_in[j][k] = 0.;
              for (m = 0;m < nb_state;m++) {
                state_in[j][k] += transition[m][k] * forward1[m];
              }
            }
          }

          for (k = 0;k < seq.nb_variable;k++) {
            poutput[k]++;
          }
        }

        if (likelihood == D_INF) {
          break;
        }
      }
    }

    for (i = 0;i < length;i++) {
      delete [] observation[i];
    }
    delete [] observation;

    delete [] norm;
    delete [] state_norm;
    delete [] forward1;

    for (i = 0;i < length - 1;i++) {
      delete [] state_in[i];
    }
    delete [] state_in;

    delete [] poutput;
  }

  return likelihood;
}


/*--------------------------------------------------------------*
 *
 *  Estimation des parametres d'une semi-chaine de Markov cachee a partir d'un echantillon
 *  de sequences par l'algorithme EM.
 *
 *  arguments : reference sur un objet Format_error, stream, semi-chaine de Markov cachee initiale,
 *              type d'estimateur pour la reestimation des lois d'occupation des etats,
 *              flags sur le calcul des lois de comptage et sur le calcul des sequences
 *              d'etats optimales, nombre d'iterations, methode de calcul de la moyenne
 *              des lois d'occupation des etats (semi-chaine de Markov cachee en equilibre).
 *
 *--------------------------------------------------------------*/

Hidden_semi_markov* Markovian_sequences::hidden_semi_markov_estimation(Format_error &error , ostream &os ,
                                                                       const Hidden_semi_markov &ihsmarkov ,
                                                                       int estimator , bool counting_flag ,
                                                                       bool state_sequence , int nb_iter ,
                                                                       int mean_computation) const

{
  bool status;
  register int i , j , k , m , n;
  int max_nb_value , iter , nb_likelihood_decrease , offset , nb_value , *occupancy_nb_value ,
      *censored_occupancy_nb_value , **poutput;
  double likelihood = D_INF , previous_likelihood , occupancy_likelihood , observation_likelihood ,
         min_likelihood , obs_product , buff , sum , occupancy_mean , **observation , *norm ,
         *state_norm , **forward , **state_in , *backward , **backward1 , *auxiliary , *reestim ,
         *ofrequency , *lfrequency , *occupancy_survivor , *censored_occupancy_survivor;
  Chain_reestimation<double> *chain_reestim;
  Parametric *occupancy;
  Reestimation<double> **occupancy_reestim , **length_bias_reestim , **censored_occupancy_reestim ,
                       ***observation_reestim;
  Histogram *hoccupancy , *hobservation;
  Hidden_semi_markov *hsmarkov;
  Semi_markov_data *seq;

# ifdef DEBUG
  double test[NB_STATE][4];
# endif


  hsmarkov = 0;
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

  if (ihsmarkov.nb_output_process != nb_variable) {
    status = false;
    error.update(SEQ_error[SEQR_NB_OUTPUT_PROCESS]);
  }

  else {
    for (i = 0;i < nb_variable;i++) {
      if (((ihsmarkov.nonparametric_process[i + 1]) &&
           (ihsmarkov.nonparametric_process[i + 1]->nb_value != marginal[i]->nb_value)) ||
          ((ihsmarkov.parametric_process[i + 1]) &&
           (ihsmarkov.parametric_process[i + 1]->nb_value < marginal[i]->nb_value))) {
        status = false;
        ostringstream error_message;
        error_message << STAT_label[STATL_OUTPUT_PROCESS] << " " << i + 1 << ": "
                      << SEQ_error[SEQR_NB_OUTPUT];
        error.update((error_message.str()).c_str());
      }

      else if ((ihsmarkov.nonparametric_process[i + 1]) && (!characteristics[i])) {
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
    if (max_length > COUNTING_MAX_LENGTH) {
      counting_flag = false;
    }

    // creation de la semi-chaine de Markov cachee

    hsmarkov = new Hidden_semi_markov(ihsmarkov , false , (int)(max_length * SAMPLE_NB_VALUE_COEFF));

    if (hsmarkov->type == 'e') {
      for (i = 0;i < hsmarkov->nb_state;i++) {
        hsmarkov->initial[i] = 1. / (double)hsmarkov->nb_state;
      }
    }

#   ifdef DEBUG
    cout << *hsmarkov;
#   endif

    // creation des structures de donnees de l'algorithme

    observation = new double*[max_length];
    for (i = 0;i < max_length;i++) {
      observation[i] = new double[hsmarkov->nb_state];
    }

    norm = new double[max_length];
    state_norm = new double[hsmarkov->nb_state];

    forward = new double*[max_length];
    for (i = 0;i < max_length;i++) {
      forward[i] = new double[hsmarkov->nb_state];
    }

    state_in = new double*[max_length - 1];
    for (i = 0;i < max_length - 1;i++) {
      state_in[i] = new double[hsmarkov->nb_state];
    }

    backward = new double[hsmarkov->nb_state];

    backward1 = new double*[max_length];
    for (i = 0;i < max_length;i++) {
      backward1[i] = new double[hsmarkov->nb_state];
    }

    auxiliary = new double[hsmarkov->nb_state];

    chain_reestim = new Chain_reestimation<double>((hsmarkov->type == 'o' ?  'o' : 'v') ,
                                                   hsmarkov->nb_state , hsmarkov->nb_state);

    occupancy_nb_value = new int[hsmarkov->nb_state];
    occupancy_reestim = new Reestimation<double>*[hsmarkov->nb_state];
    if (hsmarkov->type == 'e') {
      length_bias_reestim = new Reestimation<double>*[hsmarkov->nb_state];
    }

    for (i = 0;i < hsmarkov->nb_state;i++) {
      switch (hsmarkov->state_subtype[i]) {

      case SEMI_MARKOVIAN : {
        if (estimator == COMPLETE_LIKELIHOOD) {
          occupancy_nb_value[i] = hsmarkov->nonparametric_process[0]->sojourn_time[i]->alloc_nb_value;
        }
        else {
          occupancy_nb_value[i] = MIN(hsmarkov->nonparametric_process[0]->sojourn_time[i]->alloc_nb_value ,
                                      max_length);
        }

        occupancy_reestim[i] = new Reestimation<double>(occupancy_nb_value[i]);
        if (hsmarkov->type == 'e') {
          length_bias_reestim[i] = new Reestimation<double>(occupancy_nb_value[i]);
        }
        break;
      }

      case MARKOVIAN : {
        occupancy_reestim[i] = 0;
        if (hsmarkov->type == 'e') {
          length_bias_reestim[i] = 0;
        }
        break;
      }
      }
    }

    max_nb_value = 0;
    for (i = 0;i < hsmarkov->nb_state;i++) {
      if ((hsmarkov->state_subtype[i] == SEMI_MARKOVIAN) && (occupancy_nb_value[i] > max_nb_value)) {
        max_nb_value = occupancy_nb_value[i];
      }
    }

    if (estimator == KAPLAN_MEIER) {
      censored_occupancy_nb_value = new int[hsmarkov->nb_state];
      censored_occupancy_reestim = new Reestimation<double>*[hsmarkov->nb_state];
      for (i = 0;i < hsmarkov->nb_state;i++) {
        switch (hsmarkov->state_subtype[i]) {
        case SEMI_MARKOVIAN :
          censored_occupancy_nb_value[i] = MIN(hsmarkov->nonparametric_process[0]->sojourn_time[i]->alloc_nb_value ,
                                               max_length + 1);
          censored_occupancy_reestim[i] = new Reestimation<double>(censored_occupancy_nb_value[i]);
          break;
        case MARKOVIAN :
          censored_occupancy_reestim[i] = 0;
          break;
        }
      }

      occupancy_survivor = new double[max_nb_value];
      censored_occupancy_survivor = new double[max_nb_value + 1];
    }

    hoccupancy = new Histogram(max_nb_value);

    observation_reestim = new Reestimation<double>**[hsmarkov->nb_output_process];
    for (i = 0;i < hsmarkov->nb_output_process;i++) {
      observation_reestim[i] = new Reestimation<double>*[hsmarkov->nb_state];
      for (j = 0;j < hsmarkov->nb_state;j++) {
        observation_reestim[i][j] = new Reestimation<double>(marginal[i]->nb_value);
      }
    }

    max_nb_value = 0;
    for (i = 0;i < hsmarkov->nb_output_process;i++) {
      if ((hsmarkov->parametric_process[i + 1]) && (max_nb_value < marginal[i]->nb_value)) {
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
    nb_likelihood_decrease = 0;

    do {
      iter++;
      previous_likelihood = likelihood;
      likelihood = 0.;

      // initialisation des quantites de reestimation

      chain_reestim->init();

      for (i = 0;i < hsmarkov->nb_state;i++) {
        if (hsmarkov->state_subtype[i] == SEMI_MARKOVIAN) {
          reestim = occupancy_reestim[i]->frequency;
          for (j = 0;j < occupancy_nb_value[i];j++) {
            *reestim++ = 0.;
          }

          if (hsmarkov->type == 'e') {
            reestim = length_bias_reestim[i]->frequency;
            for (j = 0;j < occupancy_nb_value[i];j++) {
              *reestim++ = 0.;
            }
          }

          if (estimator == KAPLAN_MEIER) {
            reestim = censored_occupancy_reestim[i]->frequency;
            for (j = 0;j < censored_occupancy_nb_value[i];j++) {
              *reestim++ = 0.;
            }
          }
        }
      }

      for (i = 0;i < hsmarkov->nb_output_process;i++) {
        for (j = 0;j < hsmarkov->nb_state;j++) {
          reestim = observation_reestim[i][j]->frequency;
          for (k = 0;k < marginal[i]->nb_value;k++) {
            *reestim++ = 0.;
          }
        }
      }

#     ifdef DEBUG
      for (i = 0;i < hsmarkov->nb_state;i++) {
        for (j = 0;j < 4;j++) {
          test[i][j] = 0.;
        }
      }
#     endif

      for (i = 0;i < nb_sequence;i++) {
        for (j = 0;j < nb_variable;j++) {
          poutput[j] = sequence[i][j];
        }

        // recurrence "forward"

        for (j = 0;j < length[i];j++) {
          norm[j] = 0.;

          for (k = 0;k < hsmarkov->nb_state;k++) {

            // calcul des probabilites d'observation

            observation[j][k] = 1.;
            for (m = 0;m < hsmarkov->nb_output_process;m++) {
              if (hsmarkov->nonparametric_process[m + 1]) {
                observation[j][k] *= hsmarkov->nonparametric_process[m + 1]->observation[k]->mass[*poutput[m]];
              }
              else {
                observation[j][k] *= hsmarkov->parametric_process[m + 1]->observation[k]->mass[*poutput[m]];
              }
            }

            switch (hsmarkov->state_subtype[k]) {

            // cas etat semi-markovien

            case SEMI_MARKOVIAN : {
              if (j == 0) {
                state_norm[k] = hsmarkov->initial[k];
              }
              else {
                state_norm[k] += state_in[j - 1][k] - forward[j - 1][k];
              }
              state_norm[k] *= observation[j][k];

              norm[j] += state_norm[k];
              break;
            }

            // cas etat markovien

            case MARKOVIAN : {
              if (j == 0) {
                forward[j][k] = hsmarkov->initial[k];
              }
              else {
                forward[j][k] = state_in[j - 1][k];
              }
              forward[j][k] *= observation[j][k];

              norm[j] += forward[j][k];
              break;
            }
            }
          }

          if (norm[j] > 0.) {
            for (k = 0;k < hsmarkov->nb_state;k++) {
              switch (hsmarkov->state_subtype[k]) {
              case SEMI_MARKOVIAN :
                state_norm[k] /= norm[j];
                break;
              case MARKOVIAN :
                forward[j][k] /= norm[j];
                break;
              }
            }

            likelihood += log(norm[j]);
          }

          else {
            likelihood = D_INF;
            break;
          }

          for (k = 0;k < hsmarkov->nb_state;k++) {

            // cas etat semi-markovien

            if (hsmarkov->state_subtype[k] == SEMI_MARKOVIAN) {
              occupancy = hsmarkov->nonparametric_process[0]->sojourn_time[k];
              obs_product = 1.;
              forward[j][k] = 0.;

              if (j < length[i] - 1) {
                for (m = 1;m <= MIN(j + 1 , occupancy->nb_value - 1);m++) {
                  obs_product *= observation[j - m + 1][k] / norm[j - m + 1];
                  if (obs_product == 0.) {
                    break;
                  }

                  if (m < j + 1) {
                    forward[j][k] += obs_product * occupancy->mass[m] * state_in[j - m][k];
                  }

                  else {
                    switch (hsmarkov->type) {
                    case 'o' :
                      forward[j][k] += obs_product * occupancy->mass[m] * hsmarkov->initial[k];
                      break;
                    case 'e' :
                      forward[j][k] += obs_product * hsmarkov->forward[k]->mass[m] * hsmarkov->initial[k];
                      break;
                    }
                  }
                }
              }

              else {
                for (m = 1;m <= MIN(j + 1 , occupancy->nb_value - 1);m++) {
                  obs_product *= observation[j - m + 1][k] / norm[j - m + 1];
                  if (obs_product == 0.) {
                    break;
                  }

                  if (m < j + 1) {
                    forward[j][k] += obs_product * (1. - occupancy->cumul[m - 1]) * state_in[j - m][k];
                  }

                  else {
                    switch (hsmarkov->type) {
                    case 'o' :
                      forward[j][k] += obs_product * (1. - occupancy->cumul[m - 1]) *
                                       hsmarkov->initial[k];
                      break;
                    case 'e' :
                      forward[j][k] += obs_product * (1. - hsmarkov->forward[k]->cumul[m - 1]) *
                                       hsmarkov->initial[k];
                      break;
                    }
                  }
                }
              }
            }
          }

          if (j < length[i] - 1) {
            for (k = 0;k < hsmarkov->nb_state;k++) {
              state_in[j][k] = 0.;
              for (m = 0;m < hsmarkov->nb_state;m++) {
                state_in[j][k] += hsmarkov->transition[m][k] * forward[j][m];
              }
            }
          }

          for (k = 0;k < nb_variable;k++) {
            poutput[k]++;
          }
        }

        if (likelihood == D_INF) {
          break;
        }

#       ifdef DEBUG
        for (j = 0;j < length[i];j++) {
          cout << j << " : ";
          for (k = 0;k < hsmarkov->nb_state;k++) {
            cout << forward[j][k] << " ";
          }
          cout << endl;
        }
        cout << endl;
#       endif

        // recurrence "backward"

        for (j = 0;j < nb_variable;j++) {
          poutput[j]--;
        }

        j = length[i] - 1;
        for (k = 0;k < hsmarkov->nb_state;k++) {
          backward[k] = forward[j][k];
          backward1[j][k] = backward[k];

          // accumulation des quantites de reestimation des lois d'observation

          for (m = 0;m < hsmarkov->nb_output_process;m++) {
            observation_reestim[m][k]->frequency[*poutput[m]] += backward[k];
          }
        }

        for (j = length[i] - 2;j >= 0;j--) {
          for (k = 0;k < nb_variable;k++) {
            poutput[k]--;
          }

          for (k = 0;k < hsmarkov->nb_state;k++) {
            auxiliary[k] = 0.;

            switch (hsmarkov->state_subtype[k]) {

            // cas etat semi-markovien

            case SEMI_MARKOVIAN : {
              occupancy = hsmarkov->nonparametric_process[0]->sojourn_time[k];
              obs_product = 1.;

              for (m = 1;m < MIN(length[i] - j , occupancy->nb_value);m++) {
                obs_product *= observation[j + m][k] / norm[j + m];
                if (obs_product == 0.) {
                  break;
                }

                if (backward1[j + m][k] > 0.) {
//                if (forward[j + m][k] > 0.) {
                  if (m < length[i] - j - 1) {
                    buff = backward1[j + m][k] * obs_product * occupancy->mass[m] /
                           forward[j + m][k];

                    // accumulation des quantites de reestimation des lois d'occupation des etats

                    occupancy_reestim[k]->frequency[m] += buff * state_in[j][k];
                  }

                  else {
                    buff = obs_product * (1. - occupancy->cumul[m - 1]);

                    // accumulation des quantites de reestimation des lois d'occupation des etats

                    switch (estimator) {

                    case COMPLETE_LIKELIHOOD : {
                      for (n = m;n < occupancy->nb_value;n++) {
                        occupancy_reestim[k]->frequency[n] += obs_product * occupancy->mass[n] *
                                                              state_in[j][k];
                      }
                      break;
                    }

                    case KAPLAN_MEIER : {
                      censored_occupancy_reestim[k]->frequency[m] += buff * state_in[j][k];
                      break;
                    }
                    }
                  }

                  auxiliary[k] += buff;
                }
              }
              break;
            }

            // cas etat markovien

            case MARKOVIAN : {
              if (backward1[j + 1][k] > 0.) {
//              if (forward[j + 1][k] > 0.) {
                auxiliary[k] = backward1[j + 1][k] / state_in[j][k];

/*                auxiliary[k] = backward1[j + 1][k] * observation[j + 1][k] /
                               (forward[j + 1][k] * norm[j + 1]); */

              }
              break;
            }
            }
          }

          for (k = 0;k < hsmarkov->nb_state;k++) {
            backward1[j][k] = 0.;

            for (m = 0;m < hsmarkov->nb_state;m++) {
              buff = auxiliary[m] * hsmarkov->transition[k][m] * forward[j][k];
              backward1[j][k] += buff;

              // accumulation des quantites de reestimation des probabilites de transition

              chain_reestim->transition[k][m] += buff;
            }

            switch (hsmarkov->state_subtype[k]) {

            // cas etat semi-markovien

            case SEMI_MARKOVIAN : {
              backward[k] = backward[k] + backward1[j][k] - auxiliary[k] * state_in[j][k];
              if (backward[k] < 0.) {
                backward[k] = 0.;
              }
              if (backward[k] > 1.) {
                backward[k] = 1.;
              }
              break;
            }

            // cas etat markovien

            case MARKOVIAN : {
              backward[k] = backward1[j][k];
              break;
            }
            }

            // accumulation des quantites de reestimation des lois d'observation

            for (m = 0;m < hsmarkov->nb_output_process;m++) {
              observation_reestim[m][k]->frequency[*poutput[m]] += backward[k];
            }
          }
        }

        // accumulation des quantites de reestimation des probabilites initiales

        if (hsmarkov->type == 'o') {
          for (j = 0;j < hsmarkov->nb_state;j++) {
            chain_reestim->initial[j] += backward[j];
          }
        }

        // accumulation des quantites de reestimation des lois d'occupation des etats initiaux

        if ((hsmarkov->type == 'o') || (estimator == COMPLETE_LIKELIHOOD)) {
          for (j = 0;j < hsmarkov->nb_state;j++) {
            if ((hsmarkov->state_subtype[j] == SEMI_MARKOVIAN) && (hsmarkov->initial[j] > 0.)) {
              occupancy = hsmarkov->nonparametric_process[0]->sojourn_time[j];
              obs_product = 1.;
              if (hsmarkov->type == 'e') {
                sum = 0.;
              }

              for (k = 1;k < MIN(length[i] + 1 , occupancy->nb_value);k++) {
                obs_product *= observation[k - 1][j] / norm[k - 1];
                if (obs_product == 0.) {
                  break;
                }

                if (backward1[k - 1][j] > 0.) {
//                if (forward[k - 1][j] > 0.) {
                  if (k < length[i]) {
                    switch (hsmarkov->type) {
                    case 'o' :
                      occupancy_reestim[j]->frequency[k] += backward1[k - 1][j] * obs_product *
                                                            occupancy->mass[k] * hsmarkov->initial[j] /
                                                            forward[k - 1][j];
                      break;
                    case 'e' :
                      sum += backward1[k - 1][j] * obs_product / forward[k - 1][j];
                      length_bias_reestim[j]->frequency[k] += sum * occupancy->mass[k] * hsmarkov->initial[j] /
                                                              occupancy->mean;
                      break;
                    }
                  }

                  else {
                    switch (estimator) {

                    case COMPLETE_LIKELIHOOD : {
                      for (m = k;m < occupancy->nb_value;m++) {
                        switch (hsmarkov->type) {
                        case 'o' :
                          occupancy_reestim[j]->frequency[m] += obs_product * occupancy->mass[m] *
                                                                hsmarkov->initial[j];
                          break;
                        case 'e' :
                          length_bias_reestim[j]->frequency[m] += (sum + obs_product * (m + 1 - k)) * occupancy->mass[m] *
                                                                  hsmarkov->initial[j] / occupancy->mean;
                          break;
                        }
                      }
                      break;
                    }

                    case KAPLAN_MEIER : {
                      censored_occupancy_reestim[j]->frequency[k] += obs_product *
                                                                     (1. - occupancy->cumul[k - 1]) *
                                                                     hsmarkov->initial[j];
                      break;
                    }
                    }
                  }
                }
              }
            }
          }
        }

#       ifdef DEBUG
        for (j = length[i] - 1;j >= 0;j--) {
          cout << j << " : ";
          double sum = 0.;
          for (k = 0;k < hsmarkov->nb_state;k++) {
            sum += backward[k];
            cout << backward[k];
            if ((hsmarkov->state_subtype[k] == SEMI_MARKOVIAN) && (j < length[i] - 1)){
              cout << " (" << backward1[j][k] << ") ";
            }
          }
          cout << "| " << sum << endl;

          for (k = 0;k < hsmarkov->nb_state;k++) {
            if (hsmarkov->state_subtype[k] == SEMI_MARKOVIAN) {
              if (j < length[i] - 1) {
                test[k][0] += backward1[j][k];
                test[k][1] += auxiliary[k] * state_in[j][k];
              }
              else {
                test[k][2] += backward[k];
              }
              if (j == 0) {
                test[k][3] += backward[k];
              }
            }
          }
        }
#       endif

      }

      if (likelihood != D_INF) {
        if (likelihood < previous_likelihood) {
          nb_likelihood_decrease++;
        }
        else {
          nb_likelihood_decrease = 0;
        }

        // reestimation des probabilites initiales

        if (hsmarkov->type == 'o') {
          reestimation(hsmarkov->nb_state , chain_reestim->initial ,
                       hsmarkov->initial , MIN_PROBABILITY , false);
        }

        // reestimation des probabilites de transition

        for (i = 0;i < hsmarkov->nb_state;i++) {
          reestimation(hsmarkov->nb_state , chain_reestim->transition[i] ,
                       hsmarkov->transition[i] , MIN_PROBABILITY , false);
        }

        // reestimation des lois d'occupation des etats

        min_likelihood = 0.;

        for (i = 0;i < hsmarkov->nb_state;i++) {
          if (hsmarkov->state_subtype[i] == SEMI_MARKOVIAN) {
            occupancy = hsmarkov->nonparametric_process[0]->sojourn_time[i];

            if (estimator == KAPLAN_MEIER) {
              occupancy_reestim[i]->nb_value_computation();
              occupancy_reestim[i]->offset_computation();
              occupancy_reestim[i]->nb_element_computation();

              censored_occupancy_reestim[i]->nb_value_computation();
              censored_occupancy_reestim[i]->offset_computation();
              censored_occupancy_reestim[i]->nb_element_computation();

              if (censored_occupancy_reestim[i]->nb_element > 0.) {

#               ifdef DEBUG
                cout << "\n" << STAT_label[STATL_STATE] << " " << i << " (" << test[i][2]
                     << " | " << censored_occupancy_reestim[i]->nb_element << ") - ";

                occupancy_reestim[i]->max_computation();
                occupancy_reestim[i]->mean_computation();
                occupancy_reestim[i]->variance_computation();

                occupancy_reestim[i]->ascii_characteristic_print(cout);
#               endif

                occupancy_reestim[i]->state_occupancy_estimation(censored_occupancy_reestim[i] ,
                                                                 occupancy_reestim[i] ,
                                                                 occupancy_survivor ,
                                                                 censored_occupancy_survivor , false);
              }
            }

#           ifdef DEBUG
            cout << STAT_label[STATL_STATE] << " " << i << " (";
#           endif

            if ((hsmarkov->type == 'o') || (estimator == PARTIAL_LIKELIHOOD)) {
              occupancy_reestim[i]->nb_value_computation();
              occupancy_reestim[i]->offset_computation();
              occupancy_reestim[i]->nb_element_computation();
              occupancy_reestim[i]->max_computation();
              occupancy_reestim[i]->mean_computation();
              occupancy_reestim[i]->variance_computation();

#             ifdef DEBUG
              if (hsmarkov->type == 'o') {
                switch (estimator) {
                case COMPLETE_LIKELIHOOD :
                  cout << test[i][0] + test[i][2] << " | " << test[i][1] + test[i][3];
                  break;
                case PARTIAL_LIKELIHOOD :
                  cout << test[i][0];
                  break;
                }
                cout << " | " << occupancy_reestim[i]->nb_element << ") - ";
                occupancy_reestim[i]->ascii_characteristic_print(cout);
              }
#             endif

            }

            else {
              offset = 1;
              nb_value = occupancy_nb_value[i];

              ofrequency = occupancy_reestim[i]->frequency + occupancy_nb_value[i];
              lfrequency = length_bias_reestim[i]->frequency + occupancy_nb_value[i];
              while ((*--ofrequency == 0) && (*--lfrequency == 0) && (nb_value > 2)) {
                nb_value--;
              }
              occupancy_reestim[i]->nb_value = nb_value;
              length_bias_reestim[i]->nb_value = nb_value;

              ofrequency = occupancy_reestim[i]->frequency + offset;
              lfrequency = length_bias_reestim[i]->frequency + offset;
              while ((*ofrequency++ == 0) && (*lfrequency++ == 0) && (offset < nb_value - 1)) {
                offset++;
              }
              occupancy_reestim[i]->offset = offset;
              length_bias_reestim[i]->offset = offset;

              occupancy_reestim[i]->nb_element_computation();
              length_bias_reestim[i]->nb_element_computation();

#             ifdef DEBUG
              occupancy_reestim[i]->max_computation();
              occupancy_reestim[i]->mean_computation();
              occupancy_reestim[i]->variance_computation();

              cout << test[i][1] << " | " << occupancy_reestim[i]->nb_element << ") - ";
              occupancy_reestim[i]->ascii_characteristic_print(cout);

              length_bias_reestim[i]->max_computation();
              length_bias_reestim[i]->mean_computation();
              length_bias_reestim[i]->variance_computation();

              cout << STAT_label[STATL_STATE] << " " << i << " (" << test[i][3] << " | "
                   << length_bias_reestim[i]->nb_element << ") - ";
              length_bias_reestim[i]->ascii_characteristic_print(cout);
#             endif

              switch (mean_computation) {
              case COMPUTED :
                occupancy_mean = interval_bisection(occupancy_reestim[i] , length_bias_reestim[i]);
                break;
              case ONE_STEP_LATE :
                occupancy_mean = occupancy->mean;
                break;
              }

#             ifdef DEBUG
              cout << SEQ_label[SEQL_SOJOURN_TIME] << " " << STAT_label[STATL_MEAN] << ": "
                   << occupancy_mean << endl;
#             endif

              occupancy_reestim[i]->equilibrium_process_combination(length_bias_reestim[i] , occupancy_mean);

#             ifdef DEBUG
              cout << test[i][0] + test[i][2] << " | " << test[i][1] + test[i][3] << " | "
                   << occupancy_reestim[i]->nb_element << ") - ";
              occupancy_reestim[i]->ascii_characteristic_print(cout);
#             endif
            }

            hoccupancy->update(occupancy_reestim[i] ,
                               MAX((int)(occupancy_reestim[i]->nb_element *
                                         MAX(sqrt(occupancy_reestim[i]->variance) , 1.) * OCCUPANCY_COEFF) , MIN_NB_ELEMENT));
            if (iter <= EXPLORATION_NB_ITER) {
              occupancy_likelihood = hoccupancy->Reestimation<int>::parametric_estimation(occupancy , 1 , true ,
                                                                                          OCCUPANCY_THRESHOLD);
            }
            else {
              occupancy_likelihood = hoccupancy->Reestimation<int>::type_parametric_estimation(occupancy , 1 , true ,
                                                                                               OCCUPANCY_THRESHOLD);
            }

            if (occupancy_likelihood == D_INF) {
              min_likelihood = D_INF;
            }
            else {
              occupancy->computation(hoccupancy->nb_value , OCCUPANCY_THRESHOLD);
              if (hsmarkov->type == 'e') {
                hsmarkov->forward[i]->copy(*occupancy);
                hsmarkov->forward[i]->computation(*occupancy);
              }
            }

#           ifdef DEBUG
            cout << STAT_word[STATW_STATE] << " " << i << " " << STAT_word[STATW_OCCUPANCY_DISTRIBUTION] << endl;
            occupancy->ascii_print(cout);
#           endif

          }
        }

        if (hsmarkov->type == 'e') {
          hsmarkov->initial_probability_computation();
        }

        // reestimation des lois d'observation

        for (i = 0;i < hsmarkov->nb_output_process;i++) {
          if (hsmarkov->nonparametric_process[i + 1]) {
            for (j = 0;j < hsmarkov->nb_state;j++) {
              reestimation(marginal[i]->nb_value , observation_reestim[i][j]->frequency ,
                           hsmarkov->nonparametric_process[i + 1]->observation[j]->mass ,
                           MIN_PROBABILITY , false);
            }
          }

          else {
            for (j = 0;j < hsmarkov->nb_state;j++) {
              observation_reestim[i][j]->nb_value_computation();
              observation_reestim[i][j]->offset_computation();
              observation_reestim[i][j]->nb_element_computation();
              observation_reestim[i][j]->max_computation();
              observation_reestim[i][j]->mean_computation();
              observation_reestim[i][j]->variance_computation();

              hobservation->update(observation_reestim[i][j] ,
                                   MAX((int)(observation_reestim[i][j]->nb_element *
                                             MAX(sqrt(observation_reestim[i][j]->variance) , 1.) * OBSERVATION_COEFF) , MIN_NB_ELEMENT));
              observation_likelihood = hobservation->Reestimation<int>::type_parametric_estimation(hsmarkov->parametric_process[i + 1]->observation[j] ,
                                                                                                   0 , true , OBSERVATION_THRESHOLD);

              if (observation_likelihood == D_INF) {
                min_likelihood = D_INF;
              }
              else {
                hsmarkov->parametric_process[i + 1]->observation[j]->computation(marginal[i]->nb_value ,
                                                                                 OBSERVATION_THRESHOLD);

                if (hsmarkov->parametric_process[i + 1]->observation[j]->ident == BINOMIAL) {
                  for (k = hsmarkov->parametric_process[i + 1]->observation[j]->nb_value;k < marginal[i]->nb_value;k++) {
                    hsmarkov->parametric_process[i + 1]->observation[j]->mass[k] = 0.;
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
        cout << *hsmarkov;
      }
#     endif

    }
    while ((likelihood != D_INF) && (((nb_iter == I_DEFAULT) && (iter < SEMI_MARKOV_NB_ITER) &&
             (((likelihood - previous_likelihood) / -likelihood > SEMI_MARKOV_LIKELIHOOD_DIFF) ||
              (min_likelihood == D_INF) || (nb_likelihood_decrease == 1))) ||
            ((nb_iter != I_DEFAULT) && (iter < nb_iter))));

    if (likelihood != D_INF) {

#     ifdef MESSAGE
      os << "\n" << iter << " " << STAT_label[STATL_ITERATIONS] << endl;
#     endif

      // reestimation des probabilites initiales

      if (hsmarkov->type == 'o') {
        reestimation(hsmarkov->nb_state , chain_reestim->initial ,
                     hsmarkov->initial , MIN_PROBABILITY , true);
      }

      // reestimation des probabilites de transition

      for (i = 0;i < hsmarkov->nb_state;i++) {
        reestimation(hsmarkov->nb_state , chain_reestim->transition[i] ,
                     hsmarkov->transition[i] , MIN_PROBABILITY , true);
      }

      if (hsmarkov->type == 'e') {
        hsmarkov->initial_probability_computation();
      }

      for (i = 0;i < hsmarkov->nb_state;i++) {
        if ((hsmarkov->state_subtype[i] == SEMI_MARKOVIAN) &&
            (hsmarkov->nonparametric_process[0]->sojourn_time[i]->mean == 1.)) {
          hsmarkov->state_subtype[i] = MARKOVIAN;
          delete hsmarkov->nonparametric_process[0]->sojourn_time[i];
          hsmarkov->nonparametric_process[0]->sojourn_time[i] = 0;
          delete hsmarkov->forward[i];
          hsmarkov->forward[i] = 0;
        }
      }

      // reestimation des lois d'observation non-parametriques

      for (i = 0;i < hsmarkov->nb_output_process;i++) {
        if (hsmarkov->nonparametric_process[i + 1]) {
          for (j = 0;j < hsmarkov->nb_state;j++) {
            reestimation(marginal[i]->nb_value , observation_reestim[i][j]->frequency ,
                         hsmarkov->nonparametric_process[i + 1]->observation[j]->mass ,
                         MIN_PROBABILITY , true);
          }
        }

        else {
          hsmarkov->parametric_process[i + 1]->nb_value_computation();
        }
      }
    }

    // destruction des structures de donnees de l'algorithme

    for (i = 0;i < max_length;i++) {
      delete [] observation[i];
    }
    delete [] observation;

    delete [] norm;
    delete [] state_norm;

    for (i = 0;i < max_length;i++) {
      delete [] forward[i];
    }
    delete [] forward;

    for (i = 0;i < max_length - 1;i++) {
      delete [] state_in[i];
    }
    delete [] state_in;

    delete [] backward;

    for (i = 0;i < max_length;i++) {
      delete [] backward1[i];
    }
    delete [] backward1;

    delete [] auxiliary;

    delete chain_reestim;

    delete [] occupancy_nb_value;

    for (i = 0;i < hsmarkov->nb_state;i++) {
      delete occupancy_reestim[i];
    }
    delete [] occupancy_reestim;

    if (hsmarkov->type == 'e') {
      for (i = 0;i < hsmarkov->nb_state;i++) {
        delete length_bias_reestim[i];
      }
      delete [] length_bias_reestim;
    }

    if (estimator == KAPLAN_MEIER) {
      delete [] censored_occupancy_nb_value;

      for (i = 0;i < hsmarkov->nb_state;i++) {
        delete censored_occupancy_reestim[i];
      }
      delete [] censored_occupancy_reestim;

      delete [] occupancy_survivor;
      delete [] censored_occupancy_survivor;
    }

    delete hoccupancy;

    for (i = 0;i < hsmarkov->nb_output_process;i++) {
      for (j = 0;j < hsmarkov->nb_state;j++) {
        delete observation_reestim[i][j];
      }
      delete [] observation_reestim[i];
    }
    delete [] observation_reestim;

    delete hobservation;

    delete [] poutput;

    if (likelihood == D_INF) {
      delete hsmarkov;
      hsmarkov = 0;
      error.update(STAT_error[STATR_ESTIMATION_FAILURE]);
    }

    else {
      if (state_sequence) {
        hsmarkov->semi_markov_data = new Semi_markov_data(*this , 0 , (hsmarkov->type == 'e' ? true : false));
        seq = hsmarkov->semi_markov_data;
        seq->type[0] = STATE;

        for (i = 1;i <= hsmarkov->nb_output_process;i++) {
          if ((hsmarkov->parametric_process[i]) && (seq->characteristics[i])) {
            delete seq->characteristics[i];
            seq->characteristics[i] = 0;
          }
        }

        hsmarkov->create_cumul();
        hsmarkov->log_computation();

        seq->likelihood = hsmarkov->viterbi(*seq);

        hsmarkov->remove_cumul();

        seq->max_value[0] = hsmarkov->nb_state - 1;
        seq->build_marginal_histogram(0);
        seq->build_characteristic(0 , true , (hsmarkov->type == 'e' ? true : false));

        seq->build_transition_count(hsmarkov);
        seq->build_observation_histogram();

        // calcul des lois d'occupation des etats

        for (i = 0;i < hsmarkov->nb_state;i++) {
          if (hsmarkov->state_subtype[i] == SEMI_MARKOVIAN) {
            hsmarkov->nonparametric_process[0]->sojourn_time[i]->computation((seq->characteristics[0] ? seq->characteristics[0]->sojourn_time[i]->nb_value : 1) ,
                                                                             OCCUPANCY_THRESHOLD);
            if (hsmarkov->state_type[i] == 'r') {
              if (hsmarkov->type == 'o') {
                hsmarkov->forward[i]->copy(*(hsmarkov->nonparametric_process[0]->sojourn_time[i]));
              }
              hsmarkov->forward[i]->computation(*(hsmarkov->nonparametric_process[0]->sojourn_time[i]));
            }
          }
        }

        // calcul des lois d'observation parametriques

        for (i = 1;i <= hsmarkov->nb_output_process;i++) {
          if (hsmarkov->parametric_process[i]) {
            for (j = 0;j < hsmarkov->nb_state;j++) {
              hsmarkov->parametric_process[i]->observation[j]->computation(seq->observation[i][j]->nb_value ,
                                                                           OBSERVATION_THRESHOLD);
            }
          }
        }

#       ifdef MESSAGE
        if (seq->characteristics[0]) {
          cout << "\n" << SEQ_label[SEQL_STATE_SEQUENCES_LIKELIHOOD] << ": " << seq->likelihood
               << " | " << hsmarkov->Semi_markov::likelihood_computation(*seq) << endl;
        }
#       endif

      }

      else {
        if (hsmarkov->type == 'o') {
          for (i = 0;i < hsmarkov->nb_state;i++) {
            if ((hsmarkov->state_subtype[i] == SEMI_MARKOVIAN) && (hsmarkov->state_type[i] == 'r')) {
              hsmarkov->forward[i]->copy(*(hsmarkov->nonparametric_process[0]->sojourn_time[i]));
              hsmarkov->forward[i]->computation(*(hsmarkov->nonparametric_process[0]->sojourn_time[i]));
            }
          }
        }

        hsmarkov->semi_markov_data = new Semi_markov_data(*this , (hsmarkov->type == 'e' ? true : false));
        seq = hsmarkov->semi_markov_data;
        seq->state_variable_init(INT_VALUE);

        for (i = 1;i <= hsmarkov->nb_output_process;i++) {
          if ((hsmarkov->parametric_process[i]) && (seq->characteristics[i - 1])) {
            delete seq->characteristics[i - 1];
            seq->characteristics[i - 1] = 0;
          }
        }
      }

      for (i = 1;i <= hsmarkov->nb_output_process;i++) {
        if (hsmarkov->nonparametric_process[i]) {
          for (j = 0;j < hsmarkov->nb_state;j++) {
            hsmarkov->nonparametric_process[i]->observation[j]->cumul_computation();

            hsmarkov->nonparametric_process[i]->observation[j]->max_computation();
//            hsmarkov->nonparametric_process[i]->observation[j]->mean_computation();
//            hsmarkov->nonparametric_process[i]->observation[j]->variance_computation();
          }
        }
      }

      // calcul de la vraisemblance et des lois caracteristiques du modele

      seq->hidden_likelihood = hsmarkov->likelihood_computation(*this);

      hsmarkov->component_computation();
      hsmarkov->characteristic_computation(*seq , counting_flag , I_DEFAULT , false);
    }
  }

  return hsmarkov;
}


/*--------------------------------------------------------------*
 *
 *  Estimation des parametres d'une semi-chaine de Markov cachee a partir
 *  d'un echantillon de sequences.
 *
 *  arguments : reference sur un objet Format_error, stream, type de processus
 *              ('o' : ordinaire, 'e' : en equilibre), nombre d'etats de la chaine de Markov,
 *              flag sur la nature de la chaine de Markov, type d'estimateur pour la reestimation
 *              des lois d'occupation des etats, flags sur le calcul des lois de comptage et
 *              sur le calcul des sequences d'etats optimales, temps moyen d'occupation d'un etat,
 *              nombre d'iterations, methode de calcul de la moyenne des lois d'occupation des etats
 *              (semi-chaine de Markov cachee en equilibre).
 *
 *--------------------------------------------------------------*/

Hidden_semi_markov* Markovian_sequences::hidden_semi_markov_estimation(Format_error &error , ostream &os ,
                                                                       char type , int nb_state , bool left_right ,
                                                                       int estimator , bool counting_flag ,
                                                                       bool state_sequence , double occupancy_mean ,
                                                                       int nb_iter , int mean_computation) const

{
  bool status = true;
  register int i;
  int nb_value[SEQUENCE_NB_VARIABLE];
  double proba;
  Hidden_semi_markov *ihsmarkov , *hsmarkov;


  hsmarkov = 0;
  error.init();

  if ((nb_state < 2) || (nb_state > NB_STATE)) {
    status = false;
    error.update(SEQ_error[SEQR_NB_STATE]);
  }
  if ((occupancy_mean != D_DEFAULT) && (occupancy_mean <= 1.)) {
    status = false;
    error.update(SEQ_error[SEQR_OCCUPANCY]);
  }

  if (status) {
    for (i = 0;i < nb_variable;i++) {
      nb_value[i] = marginal[i]->nb_value;
    }

    ihsmarkov = new Hidden_semi_markov(type , nb_state , nb_variable , nb_value);

    // initialisation des parametres de la chaine de Markov

    ihsmarkov->init(left_right , 0.);

    // initialisation des lois d'occupations des etats

    if (occupancy_mean == D_DEFAULT) {
      occupancy_mean = MAX(hlength->mean , OCCUPANCY_MEAN);
    }

    ihsmarkov->state_subtype = new int[nb_state];
    ihsmarkov->nonparametric_process[0]->absorption = new double[nb_state];
    ihsmarkov->nonparametric_process[0]->sojourn_time = new Parametric*[nb_state];
    ihsmarkov->forward = new Forward*[nb_state];

    for (i = 0;i < nb_state;i++) {
      if (ihsmarkov->state_type[i] != 'a') {
        ihsmarkov->state_subtype[i] = SEMI_MARKOVIAN;
        ihsmarkov->nonparametric_process[0]->absorption[i] = 0.;
        proba = 1. / occupancy_mean;
        ihsmarkov->nonparametric_process[0]->sojourn_time[i] = new Parametric(NEGATIVE_BINOMIAL , 1 , I_DEFAULT ,
                                                                              1. , proba , OCCUPANCY_THRESHOLD);

        if (ihsmarkov->state_type[i] == 'r') {
          ihsmarkov->forward[i] = new Forward(*(ihsmarkov->nonparametric_process[0]->sojourn_time[i]) ,
                                              ihsmarkov->nonparametric_process[0]->sojourn_time[i]->alloc_nb_value);
        }
        else {
          ihsmarkov->forward[i] = 0;
        }
      }

      else {
        ihsmarkov->state_subtype[i] = MARKOVIAN;
        ihsmarkov->nonparametric_process[0]->absorption[i] = 1.;
        ihsmarkov->nonparametric_process[0]->sojourn_time[i] = 0;
        ihsmarkov->forward[i] = 0;
      }
    }

    // initialisation des lois d'observation

    for (i = 0;i < ihsmarkov->nb_output_process;i++) {
      if (ihsmarkov->nonparametric_process[i + 1]) {
        ihsmarkov->nonparametric_process[i + 1]->init();
      }
      else {
        ihsmarkov->parametric_process[i + 1]->init();
      }
    }

    hsmarkov = hidden_semi_markov_estimation(error , os , *ihsmarkov , estimator , counting_flag ,
                                             state_sequence , nb_iter , mean_computation);
    delete ihsmarkov;
  }

  return hsmarkov;
}


/*--------------------------------------------------------------*
 *
 *  Estimation des parametres d'une semi-chaine de Markov cachee a partir d'un echantillon
 *  de sequences par l'algorithme SEM/MCEM.
 *
 *  arguments : reference sur un objet Format_error, stream, semi-chaine de Markov cachee initiale,
 *              parametres pour le nombre de sequences d'etats simulees, type d'estimateur
 *              pour la reestimation des lois d'occupation des etats, flags sur le calcul
 *              des lois de comptage et sur le calcul des sequences d'etats optimales,
 *              nombre d'iterations.
 *
 *--------------------------------------------------------------*/

Hidden_semi_markov* Markovian_sequences::hidden_semi_markov_stochastic_estimation(Format_error &error , ostream &os ,
                                                                                  const Hidden_semi_markov &ihsmarkov ,
                                                                                  int min_nb_state_sequence ,
                                                                                  int max_nb_state_sequence ,
                                                                                  double parameter , int estimator ,
                                                                                  bool counting_flag , bool state_sequence ,
                                                                                  int nb_iter) const

{
  bool status;
  register int i , j , k , m , n;
  int max_nb_value , iter , nb_state_sequence , state_occupancy , nb_likelihood_decrease ,
      *occupancy_nb_value , *state_seq , *pstate , **poutput;
  double likelihood = D_INF , previous_likelihood , occupancy_likelihood , observation_likelihood ,
         min_likelihood , obs_product , **observation , *norm , *state_norm , **forward ,
         **state_in , *backward , *cumul_backward , *reestim , *occupancy_survivor ,
         *censored_occupancy_survivor;
  Chain_reestimation<double> *chain_reestim;
  Parametric *occupancy;
  Reestimation<double> *bcomplete_run , *censored_run , **complete_run , **final_run , **initial_run ,
                       **single_run , ***observation_reestim;
  Hidden_semi_markov *hsmarkov;
  Semi_markov_data *seq;
  const Reestimation<double> *prun[3];

# ifdef DEBUG
  double sum;
# endif


  hsmarkov = 0;
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

  if (ihsmarkov.nb_output_process != nb_variable) {
    status = false;
    error.update(SEQ_error[SEQR_NB_OUTPUT_PROCESS]);
  }

  else {
    for (i = 0;i < nb_variable;i++) {
      if (((ihsmarkov.nonparametric_process[i + 1]) &&
           (ihsmarkov.nonparametric_process[i + 1]->nb_value != marginal[i]->nb_value)) ||
          ((ihsmarkov.parametric_process[i + 1]) &&
           (ihsmarkov.parametric_process[i + 1]->nb_value < marginal[i]->nb_value))) {
        status = false;
        ostringstream error_message;
        error_message << STAT_label[STATL_OUTPUT_PROCESS] << " " << i + 1 << ": "
                      << SEQ_error[SEQR_NB_OUTPUT];
        error.update((error_message.str()).c_str());
      }

      else if ((ihsmarkov.nonparametric_process[i + 1]) && (!characteristics[i])) {
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
    if (max_length > COUNTING_MAX_LENGTH) {
      counting_flag = false;
    }

    // creation de la semi-chaine de Markov cachee

    hsmarkov = new Hidden_semi_markov(ihsmarkov , false , (int)(max_length * SAMPLE_NB_VALUE_COEFF));

    if (hsmarkov->type == 'e') {
      for (i = 0;i < hsmarkov->nb_state;i++) {
        hsmarkov->initial[i] = 1. / (double)hsmarkov->nb_state;
      }
    }

#   ifdef DEBUG
    cout << *hsmarkov;
#   endif

    // creation des structures de donnees de l'algorithme

    observation = new double*[max_length];
    for (i = 0;i < max_length;i++) {
      observation[i] = new double[hsmarkov->nb_state];
    }

    norm = new double[max_length];
    state_norm = new double[hsmarkov->nb_state];

    forward = new double*[max_length];
    for (i = 0;i < max_length;i++) {
      forward[i] = new double[hsmarkov->nb_state];
    }

    state_in = new double*[max_length - 1];
    for (i = 0;i < max_length - 1;i++) {
      state_in[i] = new double[hsmarkov->nb_state];
    }

    backward = new double[max_length + 1];
    cumul_backward = new double[max_length + 1];

    state_seq = new int[max_length];

    chain_reestim = new Chain_reestimation<double>((hsmarkov->type == 'o' ?  'o' : 'v') ,
                                                   hsmarkov->nb_state , hsmarkov->nb_state);

    occupancy_nb_value = new int[hsmarkov->nb_state];
    complete_run = new Reestimation<double>*[hsmarkov->nb_state];
    final_run = new Reestimation<double>*[hsmarkov->nb_state];
    if (hsmarkov->type == 'e') {
      initial_run = new Reestimation<double>*[hsmarkov->nb_state];
      single_run = new Reestimation<double>*[hsmarkov->nb_state];
    }

    for (i = 0;i < hsmarkov->nb_state;i++) {
      switch (hsmarkov->state_subtype[i]) {

      case SEMI_MARKOVIAN : {
        occupancy_nb_value[i] = MIN(hsmarkov->nonparametric_process[0]->sojourn_time[i]->alloc_nb_value ,
                                    max_length + 1);

        complete_run[i] = new Reestimation<double>(occupancy_nb_value[i]);
        final_run[i] = new Reestimation<double>(occupancy_nb_value[i]);
        if (hsmarkov->type == 'e') {
          initial_run[i] = new Reestimation<double>(occupancy_nb_value[i]);
          single_run[i] = new Reestimation<double>(occupancy_nb_value[i]);
        }
        break;
      }

      case MARKOVIAN : {
        complete_run[i] = 0;
        final_run[i] = 0;
        if (hsmarkov->type == 'e') {
          initial_run[i] = 0;
          single_run[i] = 0;
        }
        break;
      }
      }
    }

    max_nb_value = 0;
    for (i = 0;i < hsmarkov->nb_state;i++) {
      if ((hsmarkov->state_subtype[i] == SEMI_MARKOVIAN) && (occupancy_nb_value[i] > max_nb_value)) {
        max_nb_value = occupancy_nb_value[i];
      }
    }

    if (estimator != PARTIAL_LIKELIHOOD) {
      occupancy_survivor = new double[max_nb_value];
      censored_occupancy_survivor = new double[max_nb_value + 1];
    }

    observation_reestim = new Reestimation<double>**[hsmarkov->nb_output_process];
    for (i = 0;i < hsmarkov->nb_output_process;i++) {
      observation_reestim[i] = new Reestimation<double>*[hsmarkov->nb_state];
      for (j = 0;j < hsmarkov->nb_state;j++) {
        observation_reestim[i][j] = new Reestimation<double>(marginal[i]->nb_value);
      }
    }

    poutput = new int*[nb_variable];

    iter = 0;
    nb_likelihood_decrease = 0;

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

      for (i = 0;i < hsmarkov->nb_state;i++) {
        if (hsmarkov->state_subtype[i] == SEMI_MARKOVIAN) {
          reestim = complete_run[i]->frequency;
          for (j = 0;j < occupancy_nb_value[i];j++) {
            *reestim++ = 0.;
          }

          reestim = final_run[i]->frequency;
          for (j = 0;j < occupancy_nb_value[i];j++) {
            *reestim++ = 0.;
          }

          if (hsmarkov->type == 'e') {
            reestim = initial_run[i]->frequency;
            for (j = 0;j < occupancy_nb_value[i];j++) {
              *reestim++ = 0.;
            }

            reestim = single_run[i]->frequency;
            for (j = 0;j < occupancy_nb_value[i];j++) {
              *reestim++ = 0.;
            }
          }
        }
      }

      for (i = 0;i < hsmarkov->nb_output_process;i++) {
        for (j = 0;j < hsmarkov->nb_state;j++) {
          reestim = observation_reestim[i][j]->frequency;
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

        for (j = 0;j < length[i];j++) {
          norm[j] = 0.;

          for (k = 0;k < hsmarkov->nb_state;k++) {

            // calcul des probabilites d'observation

            observation[j][k] = 1.;
            for (m = 0;m < hsmarkov->nb_output_process;m++) {
              if (hsmarkov->nonparametric_process[m + 1]) {
                observation[j][k] *= hsmarkov->nonparametric_process[m + 1]->observation[k]->mass[*poutput[m]];
              }
              else {
                observation[j][k] *= hsmarkov->parametric_process[m + 1]->observation[k]->mass[*poutput[m]];
              }
            }

            switch (hsmarkov->state_subtype[k]) {

            // cas etat semi-markovien

            case SEMI_MARKOVIAN : {
              if (j == 0) {
                state_norm[k] = hsmarkov->initial[k];
              }
              else {
                state_norm[k] += state_in[j - 1][k] - forward[j - 1][k];
              }
              state_norm[k] *= observation[j][k];

              norm[j] += state_norm[k];
              break;
            }

            // cas etat markovien

            case MARKOVIAN : {
              if (j == 0) {
                forward[j][k] = hsmarkov->initial[k];
              }
              else {
                forward[j][k] = state_in[j - 1][k];
              }
              forward[j][k] *= observation[j][k];

              norm[j] += forward[j][k];
              break;
            }
            }
          }

          if (norm[j] > 0.) {
            for (k = 0;k < hsmarkov->nb_state;k++) {
              switch (hsmarkov->state_subtype[k]) {
              case SEMI_MARKOVIAN :
                state_norm[k] /= norm[j];
                break;
              case MARKOVIAN :
                forward[j][k] /= norm[j];
                break;
              }
            }

            likelihood += log(norm[j]);
          }

          else {
            likelihood = D_INF;
            break;
          }

          for (k = 0;k < hsmarkov->nb_state;k++) {

            // cas etat semi-markovien

            if (hsmarkov->state_subtype[k] == SEMI_MARKOVIAN) {
              occupancy = hsmarkov->nonparametric_process[0]->sojourn_time[k];
              obs_product = 1.;
              forward[j][k] = 0.;

              if (j < length[i] - 1) {
                for (m = 1;m <= MIN(j + 1 , occupancy->nb_value - 1);m++) {
                  obs_product *= observation[j - m + 1][k] / norm[j - m + 1];
                  if (obs_product == 0.) {
                    break;
                  }

                  if (m < j + 1) {
                    forward[j][k] += obs_product * occupancy->mass[m] * state_in[j - m][k];
                  }

                  else {
                    switch (hsmarkov->type) {
                    case 'o' :
                      forward[j][k] += obs_product * occupancy->mass[m] * hsmarkov->initial[k];
                      break;
                    case 'e' :
                      forward[j][k] += obs_product * hsmarkov->forward[k]->mass[m] * hsmarkov->initial[k];
                      break;
                    }
                  }
                }
              }

              else {
                for (m = 1;m <= MIN(j + 1 , occupancy->nb_value - 1);m++) {
                  obs_product *= observation[j - m + 1][k] / norm[j - m + 1];
                  if (obs_product == 0.) {
                    break;
                  }

                  if (m < j + 1) {
                    forward[j][k] += obs_product * (1. - occupancy->cumul[m - 1]) * state_in[j - m][k];
                  }

                  else {
                    switch (hsmarkov->type) {
                    case 'o' :
                      forward[j][k] += obs_product * (1. - occupancy->cumul[m - 1]) *
                                       hsmarkov->initial[k];
                      break;
                    case 'e' :
                      forward[j][k] += obs_product * (1. - hsmarkov->forward[k]->cumul[m - 1]) *
                                       hsmarkov->initial[k];
                      break;
                    }
                  }
                }
              }
            }
          }

          if (j < length[i] - 1) {
            for (k = 0;k < hsmarkov->nb_state;k++) {
              state_in[j][k] = 0.;
              for (m = 0;m < hsmarkov->nb_state;m++) {
                state_in[j][k] += hsmarkov->transition[m][k] * forward[j][m];
              }
            }
          }

          for (k = 0;k < nb_variable;k++) {
            poutput[k]++;
          }
        }

        if (likelihood == D_INF) {
          break;
        }

#       ifdef DEBUG
        for (j = 0;j < length[i];j++) {
          cout << j << " : ";
          for (k = 0;k < hsmarkov->nb_state;k++) {
            cout << forward[j][k] << " ";
          }
          cout << endl;
        }
        cout << endl;
#       endif

        // passes "backward"

        for (j = 0;j < nb_state_sequence;j++) {
          k = length[i] - 1;
          pstate = state_seq + k;
          for (m = 0;m < nb_variable;m++) {
            poutput[m] = sequence[i][m] + k;
          }

          cumul_computation(hsmarkov->nb_state , forward[k] , cumul_backward);
          *pstate = cumul_method(hsmarkov->nb_state , cumul_backward);

          // accumulation des quantites de reestimation des lois d'observation

          for (m = 0;m < hsmarkov->nb_output_process;m++) {
            (observation_reestim[m][*pstate]->frequency[*poutput[m]])++;
          }

          do {

            // cas etat semi-markovien

            if (hsmarkov->state_subtype[*pstate] == SEMI_MARKOVIAN) {
              occupancy = hsmarkov->nonparametric_process[0]->sojourn_time[*pstate];
              obs_product = 1.;

              if (k < length[i] - 1) {
                for (m = 1;m <= MIN(k + 1 , occupancy->nb_value - 1);m++) {
                  obs_product *= observation[k - m + 1][*pstate] / norm[k - m + 1];
                  if (obs_product == 0.) {
                    break;
                  }

                  if (m < k + 1) {
                    backward[m] = obs_product * occupancy->mass[m] * state_in[k - m][*pstate] /
                                  forward[k][*pstate];
                  }

                  else {
                    switch (hsmarkov->type) {
                    case 'o' :
                      backward[m] = obs_product * occupancy->mass[m] * hsmarkov->initial[*pstate] /
                                    forward[k][*pstate];
                      break;
                    case 'e' :
                      backward[m] = obs_product * hsmarkov->forward[*pstate]->mass[m] * hsmarkov->initial[*pstate] /
                                    forward[k][*pstate];
                      break;
                    }
                  }
                }
              }

              else {
                for (m = 1;m <= MIN(k + 1 , occupancy->nb_value - 1);m++) {
                  obs_product *= observation[k - m + 1][*pstate] / norm[k - m + 1];
                  if (obs_product == 0.) {
                    break;
                  }

                  if (m < k + 1) {
                    backward[m] = obs_product * (1. - occupancy->cumul[m - 1]) * state_in[k - m][*pstate] /
                                  forward[k][*pstate];
                  }

                  else {
                    switch (hsmarkov->type) {
                    case 'o' :
                      backward[m] = obs_product * (1. - occupancy->cumul[m - 1]) *
                                    hsmarkov->initial[*pstate] / forward[k][*pstate];
                      break;
                    case 'e' :
                      backward[m] = obs_product * (1. - hsmarkov->forward[*pstate]->cumul[m - 1]) *
                                    hsmarkov->initial[*pstate] / forward[k][*pstate];
                      break;
                    }
                  }
                }
              }

              cumul_computation(m - 1 , backward + 1 , cumul_backward);
              state_occupancy = 1 + cumul_method(m - 1 , cumul_backward);

#             ifdef DEBUG
              sum = 0.;
              for (n = 1;n < m;n++) {
                sum += backward[n];
              }
              if ((sum < 1. - DOUBLE_ERROR) || (sum > 1. + DOUBLE_ERROR)) {
                cout << "\nERROR: " << k << " " << sum << endl;
              }
#             endif

              // accumulation des quantites de reestimation des lois d'occupation des etats

              if (k < length[i] - 1) {
                if (state_occupancy < k + 1) {
                  (complete_run[*pstate]->frequency[state_occupancy])++;
                }

                else {
                  switch (hsmarkov->type) {
                  case 'o' :
                    (complete_run[*pstate]->frequency[state_occupancy])++;
                    break;
                  case 'e' :
                    (initial_run[*pstate]->frequency[state_occupancy])++;
                    break;
                  }
                }
              }

              else {
                if (state_occupancy < k + 1) {
                  (final_run[*pstate]->frequency[state_occupancy])++;
                }

                else {
                  switch (hsmarkov->type) {
                  case 'o' :
                    (final_run[*pstate]->frequency[state_occupancy])++;
                    break;
                  case 'e' :
                    (single_run[*pstate]->frequency[state_occupancy])++;
                    break;
                  }
                }
              }

              for (m = 1;m < state_occupancy;m++) {
                pstate--;
                *pstate = *(pstate + 1);

                // accumulation des quantites de reestimation des lois d'observation

                for (n = 0;n < hsmarkov->nb_output_process;n++) {
                  (observation_reestim[n][*pstate]->frequency[*--poutput[n]])++;
                }
              }
              k -= (state_occupancy - 1);

              if (k == 0) {
                break;
              }
            }

            k--;
            for (m = 0;m < hsmarkov->nb_state;m++) {
              backward[m] = hsmarkov->transition[m][*pstate] * forward[k][m] / state_in[k][*pstate];
            }
            cumul_computation(hsmarkov->nb_state , backward , cumul_backward);
            *--pstate = cumul_method(hsmarkov->nb_state , cumul_backward);

            // accumulation des quantites de reestimation des probabilites de transition et
            // des lois d'observation

            (chain_reestim->transition[*pstate][*(pstate + 1)])++;

            for (m = 0;m < hsmarkov->nb_output_process;m++) {
              (observation_reestim[m][*pstate]->frequency[*--poutput[m]])++;
            }

#           ifdef DEBUG
            sum = 0.;
            for (m = 0;m < nb_state;m++) {
              sum += backward[m];
            }
            if ((sum < 1. - DOUBLE_ERROR) || (sum > 1. + DOUBLE_ERROR)) {
              cout << "\nERROR: " << k << " " << sum << endl;
            }
#           endif

          }
          while (k > 0);

          // accumulation des quantites de reestimation des probabilites initiales

          if (hsmarkov->type == 'o') {
            (chain_reestim->initial[*pstate])++;
          }
        }
      }

      if (likelihood != D_INF) {
        if (likelihood < previous_likelihood) {
          nb_likelihood_decrease++;
        }
        else {
          nb_likelihood_decrease = 0;
        }

        // reestimation des probabilites initiales

        if (hsmarkov->type == 'o') {
          reestimation(hsmarkov->nb_state , chain_reestim->initial ,
                       hsmarkov->initial , MIN_PROBABILITY , false);
        }

        // reestimation des probabilites de transition

        for (i = 0;i < hsmarkov->nb_state;i++) {
          reestimation(hsmarkov->nb_state , chain_reestim->transition[i] ,
                       hsmarkov->transition[i] , MIN_PROBABILITY , false);
        }

        // reestimation des lois d'occupation des etats

        min_likelihood = 0.;

        for (i = 0;i < hsmarkov->nb_state;i++) {
          if (hsmarkov->state_subtype[i] == SEMI_MARKOVIAN) {
            occupancy = hsmarkov->nonparametric_process[0]->sojourn_time[i];

            complete_run[i]->nb_value_computation();
            complete_run[i]->offset_computation();
            complete_run[i]->nb_element_computation();

#           ifdef DEBUG
            cout << "\n" << STAT_label[STATL_STATE] << " " << i << " ";

            complete_run[i]->max_computation();
            complete_run[i]->mean_computation();
            complete_run[i]->variance_computation();

            complete_run[i]->print(cout);
#           endif

            if ((iter > STOCHASTIC_EXPLORATION_NB_ITER) && (estimator == COMPLETE_LIKELIHOOD)) {
              final_run[i]->nb_value_computation();
              final_run[i]->offset_computation();
              final_run[i]->nb_element_computation();

              switch (hsmarkov->type) {
              case 'o' : {
                if (final_run[i]->nb_element > 0.) {
                  complete_run[i]->state_occupancy_estimation(final_run[i] , complete_run[i] ,
                                                              occupancy_survivor ,
                                                              censored_occupancy_survivor , false);
                }
                break;
              }

              case 'e' : {
                initial_run[i]->nb_value_computation();
                initial_run[i]->offset_computation();
                initial_run[i]->nb_element_computation();

                single_run[i]->nb_value_computation();
                single_run[i]->offset_computation();
                single_run[i]->nb_element_computation();

                prun[0] = complete_run[i];
                prun[1] = complete_run[i];
                bcomplete_run = new Reestimation<double>(2 , prun);

                prun[0] = initial_run[i];
                prun[1] = final_run[i];
                prun[2] = single_run[i];
                censored_run = new Reestimation<double>(3 , prun);

#               ifdef DEBUG
                censored_run->print(cout);
#               endif

                bcomplete_run->state_occupancy_estimation(censored_run , complete_run[i] ,
                                                          occupancy_survivor ,
                                                          censored_occupancy_survivor , false);
                delete bcomplete_run;
                delete censored_run;
                break;
              }
              }

              if ((hsmarkov->type == 'e') || (final_run[i]->nb_element > 0.)) {
                complete_run[i]->nb_value_computation();
                complete_run[i]->offset_computation();
                complete_run[i]->nb_element_computation();

#               ifdef DEBUG
                complete_run[i]->max_computation();
                complete_run[i]->mean_computation();
                complete_run[i]->variance_computation();

                complete_run[i]->print(cout);
#               endif

              }
            }

            complete_run[i]->max_computation();
            complete_run[i]->mean_computation();
            complete_run[i]->variance_computation();

            if (iter <= EXPLORATION_NB_ITER) {
              occupancy_likelihood = complete_run[i]->parametric_estimation(occupancy , 1 , true ,
                                                                            OCCUPANCY_THRESHOLD);
            }
            else {
              occupancy_likelihood = complete_run[i]->type_parametric_estimation(occupancy , 1 , true ,
                                                                                 OCCUPANCY_THRESHOLD);
            }

#           ifdef DEBUG
            if (i == 1) {
              occupancy->print(cout);
            }
#           endif

            if (occupancy_likelihood == D_INF) {
              min_likelihood = D_INF;
            }
            else {
              occupancy->computation(complete_run[i]->nb_value , OCCUPANCY_THRESHOLD);
              if (hsmarkov->type == 'e') {
                hsmarkov->forward[i]->copy(*occupancy);
                hsmarkov->forward[i]->computation(*occupancy);
              }
            }

#           ifdef DEBUG
            cout << STAT_word[STATW_STATE] << " " << i << " " << STAT_word[STATW_OCCUPANCY_DISTRIBUTION] << endl;
            occupancy->ascii_print(cout);
#           endif

          }
        }

        if (hsmarkov->type == 'e') {
          hsmarkov->initial_probability_computation();
        }

        // reestimation des lois d'observation

        for (i = 0;i < hsmarkov->nb_output_process;i++) {
          if (hsmarkov->nonparametric_process[i + 1]) {
            for (j = 0;j < hsmarkov->nb_state;j++) {
              reestimation(marginal[i]->nb_value , observation_reestim[i][j]->frequency ,
                           hsmarkov->nonparametric_process[i + 1]->observation[j]->mass ,
                           MIN_PROBABILITY , false);
            }
          }

          else {
            for (j = 0;j < hsmarkov->nb_state;j++) {
              observation_reestim[i][j]->nb_value_computation();
              observation_reestim[i][j]->offset_computation();
              observation_reestim[i][j]->nb_element_computation();
              observation_reestim[i][j]->max_computation();
              observation_reestim[i][j]->mean_computation();
              observation_reestim[i][j]->variance_computation();

              observation_likelihood = observation_reestim[i][j]->type_parametric_estimation(hsmarkov->parametric_process[i + 1]->observation[j] ,
                                                                                             0 , true , OBSERVATION_THRESHOLD);

              if (observation_likelihood == D_INF) {
                min_likelihood = D_INF;
              }
              else {
                hsmarkov->parametric_process[i + 1]->observation[j]->computation(marginal[i]->nb_value ,
                                                                                 OBSERVATION_THRESHOLD);

                if (hsmarkov->parametric_process[i + 1]->observation[j]->ident == BINOMIAL) {
                  for (k = hsmarkov->parametric_process[i + 1]->observation[j]->nb_value;k < marginal[i]->nb_value;k++) {
                    hsmarkov->parametric_process[i + 1]->observation[j]->mass[k] = 0.;
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
        cout << *hsmarkov;
      }
#     endif

    }
    while ((likelihood != D_INF) && ((iter <= STOCHASTIC_EXPLORATION_NB_ITER + 2) ||
            ((nb_iter == I_DEFAULT) && (iter < SEMI_MARKOV_NB_ITER) &&
             (((likelihood - previous_likelihood) / -likelihood > SEMI_MARKOV_LIKELIHOOD_DIFF) ||
              (min_likelihood == D_INF) || (nb_likelihood_decrease == 1))) ||
            ((nb_iter != I_DEFAULT) && (iter < nb_iter))));

    if (likelihood != D_INF) {

#     ifdef MESSAGE
      os << "\n" << iter << " " << STAT_label[STATL_ITERATIONS] << endl;

      if (hsmarkov->type == 'e') {
        for (i = 0;i < hsmarkov->nb_state;i++) {
          if (single_run[i]->nb_element > 0) {
            os << "\n" << SEQ_label[SEQL_BIASED] << " " << STAT_label[STATL_STATE] << " " << i
               << " " << SEQ_label[SEQL_OCCUPANCY] << " "  << STAT_label[STATL_DISTRIBUTION] << endl;
          }
        }
      }
#     endif

      // reestimation des probabilites initiales

      if (hsmarkov->type == 'o') {
        reestimation(hsmarkov->nb_state , chain_reestim->initial ,
                     hsmarkov->initial , MIN_PROBABILITY , true);
      }

      // reestimation des probabilites de transition

      for (i = 0;i < hsmarkov->nb_state;i++) {
        reestimation(hsmarkov->nb_state , chain_reestim->transition[i] ,
                     hsmarkov->transition[i] , MIN_PROBABILITY , true);
      }

      if (hsmarkov->type == 'e') {
        hsmarkov->initial_probability_computation();
      }

      for (i = 0;i < hsmarkov->nb_state;i++) {
        if ((hsmarkov->state_subtype[i] == SEMI_MARKOVIAN) &&
            (hsmarkov->nonparametric_process[0]->sojourn_time[i]->mean == 1.)) {
          hsmarkov->state_subtype[i] = MARKOVIAN;
          delete hsmarkov->nonparametric_process[0]->sojourn_time[i];
          hsmarkov->nonparametric_process[0]->sojourn_time[i] = 0;
          delete hsmarkov->forward[i];
          hsmarkov->forward[i] = 0;
        }
      }

      // reestimation des lois d'observation non-parametriques

      for (i = 0;i < hsmarkov->nb_output_process;i++) {
        if (hsmarkov->nonparametric_process[i + 1]) {
          for (j = 0;j < hsmarkov->nb_state;j++) {
            reestimation(marginal[i]->nb_value , observation_reestim[i][j]->frequency ,
                         hsmarkov->nonparametric_process[i + 1]->observation[j]->mass ,
                         MIN_PROBABILITY , true);
          }
        }

        else {
          hsmarkov->parametric_process[i + 1]->nb_value_computation();
        }
      }
    }

    // destruction des structures de donnees de l'algorithme

    for (i = 0;i < max_length;i++) {
      delete [] observation[i];
    }
    delete [] observation;

    delete [] norm;
    delete [] state_norm;

    for (i = 0;i < max_length;i++) {
      delete [] forward[i];
    }
    delete [] forward;

    for (i = 0;i < max_length - 1;i++) {
      delete [] state_in[i];
    }
    delete [] state_in;

    delete [] backward;
    delete [] cumul_backward;

    delete [] state_seq;

    delete chain_reestim;

    for (i = 0;i < hsmarkov->nb_state;i++) {
      delete complete_run[i];
    }
    delete [] complete_run;

    for (i = 0;i < hsmarkov->nb_state;i++) {
      delete final_run[i];
    }
    delete [] final_run;

    if (hsmarkov->type == 'e') {
      for (i = 0;i < hsmarkov->nb_state;i++) {
        delete initial_run[i];
      }
      delete [] initial_run;

      for (i = 0;i < hsmarkov->nb_state;i++) {
        delete single_run[i];
      }
      delete [] single_run;
    }

    delete [] occupancy_nb_value;

    if (estimator != PARTIAL_LIKELIHOOD) {
      delete [] occupancy_survivor;
      delete [] censored_occupancy_survivor;
    }

    for (i = 0;i < hsmarkov->nb_output_process;i++) {
      for (j = 0;j < hsmarkov->nb_state;j++) {
        delete observation_reestim[i][j];
      }
      delete [] observation_reestim[i];
    }
    delete [] observation_reestim;

    delete [] poutput;

    if (likelihood == D_INF) {
      delete hsmarkov;
      hsmarkov = 0;
      error.update(STAT_error[STATR_ESTIMATION_FAILURE]);
    }

    else {
      if (state_sequence) {
        hsmarkov->semi_markov_data = new Semi_markov_data(*this , 0 , (hsmarkov->type == 'e' ? true : false));
        seq = hsmarkov->semi_markov_data;
        seq->type[0] = STATE;

        for (i = 1;i <= hsmarkov->nb_output_process;i++) {
          if ((hsmarkov->parametric_process[i]) && (seq->characteristics[i])) {
            delete seq->characteristics[i];
            seq->characteristics[i] = 0;
          }
        }

        hsmarkov->create_cumul();
        hsmarkov->log_computation();

        seq->likelihood = hsmarkov->viterbi(*seq);

        hsmarkov->remove_cumul();

        seq->max_value[0] = hsmarkov->nb_state - 1;
        seq->build_marginal_histogram(0);
        seq->build_characteristic(0 , true , (hsmarkov->type == 'e' ? true : false));

        seq->build_transition_count(hsmarkov);
        seq->build_observation_histogram();

        // calcul des lois d'occupation des etats

        for (i = 0;i < hsmarkov->nb_state;i++) {
          if (hsmarkov->state_subtype[i] == SEMI_MARKOVIAN) {
            hsmarkov->nonparametric_process[0]->sojourn_time[i]->computation((seq->characteristics[0] ? seq->characteristics[0]->sojourn_time[i]->nb_value : 1) ,
                                                                             OCCUPANCY_THRESHOLD);
            if (hsmarkov->state_type[i] == 'r') {
              if (hsmarkov->type == 'o') {
                hsmarkov->forward[i]->copy(*(hsmarkov->nonparametric_process[0]->sojourn_time[i]));
              }
              hsmarkov->forward[i]->computation(*(hsmarkov->nonparametric_process[0]->sojourn_time[i]));
            }
          }
        }

        // calcul des lois d'observation parametriques

        for (i = 1;i <= hsmarkov->nb_output_process;i++) {
          if (hsmarkov->parametric_process[i]) {
            for (j = 0;j < hsmarkov->nb_state;j++) {
              hsmarkov->parametric_process[i]->observation[j]->computation(seq->observation[i][j]->nb_value ,
                                                                           OBSERVATION_THRESHOLD);
            }
          }
        }

#       ifdef MESSAGE
        if (seq->characteristics[0]) {
          cout << "\n" << SEQ_label[SEQL_STATE_SEQUENCES_LIKELIHOOD] << ": " << seq->likelihood
               << " | " << hsmarkov->Semi_markov::likelihood_computation(*seq) << endl;
        }
#       endif

      }

      else {
        if (hsmarkov->type == 'o') {
          for (i = 0;i < hsmarkov->nb_state;i++) {
            if ((hsmarkov->state_subtype[i] == SEMI_MARKOVIAN) && (hsmarkov->state_type[i] == 'r')) {
              hsmarkov->forward[i]->copy(*(hsmarkov->nonparametric_process[0]->sojourn_time[i]));
              hsmarkov->forward[i]->computation(*(hsmarkov->nonparametric_process[0]->sojourn_time[i]));
            }
          }
        }

        hsmarkov->semi_markov_data = new Semi_markov_data(*this , (hsmarkov->type == 'e' ? true : false));
        seq = hsmarkov->semi_markov_data;
        seq->state_variable_init(INT_VALUE);

        for (i = 1;i <= hsmarkov->nb_output_process;i++) {
          if ((hsmarkov->parametric_process[i]) && (seq->characteristics[i - 1])) {
            delete seq->characteristics[i - 1];
            seq->characteristics[i - 1] = 0;
          }
        }
      }

      for (i = 1;i <= hsmarkov->nb_output_process;i++) {
        if (hsmarkov->nonparametric_process[i]) {
          for (j = 0;j < hsmarkov->nb_state;j++) {
            hsmarkov->nonparametric_process[i]->observation[j]->cumul_computation();

            hsmarkov->nonparametric_process[i]->observation[j]->max_computation();
//            hsmarkov->nonparametric_process[i]->observation[j]->mean_computation();
//            hsmarkov->nonparametric_process[i]->observation[j]->variance_computation();
          }
        }
      }

      // calcul de la vraisemblance et des lois caracteristiques du modele

      seq->hidden_likelihood = hsmarkov->likelihood_computation(*this);

      hsmarkov->component_computation();
      hsmarkov->characteristic_computation(*seq , counting_flag , I_DEFAULT , false);
    }
  }

  return hsmarkov;
}


/*--------------------------------------------------------------*
 *
 *  Estimation des parametres d'une semi-chaine de Markov cachee a partir
 *  d'un echantillon de sequences.
 *
 *  arguments : reference sur un objet Format_error, stream, type de processus
 *              ('o' : ordinaire, 'e' : en equilibre), nombre d'etats de la chaine de Markov,
 *              flag sur la nature de la chaine de Markov, parametres pour nombre de sequences
 *              type d'estimateur pour la reestimation des lois d'occupation des etats,
 *              d'etats simulees, flags sur le calcul des lois de comptage et sur le calcul
 *              des sequences d'etats optimales, temps moyen d'occupation d'un etat,
 *              nombre d'iterations.
 *
 *--------------------------------------------------------------*/

Hidden_semi_markov* Markovian_sequences::hidden_semi_markov_stochastic_estimation(Format_error &error , ostream &os ,
                                                                                  char type , int nb_state , bool left_right ,
                                                                                  int min_nb_state_sequence ,
                                                                                  int max_nb_state_sequence ,
                                                                                  double parameter , int estimator ,
                                                                                  bool counting_flag , bool state_sequence ,
                                                                                  double occupancy_mean , int nb_iter) const

{
  bool status = true;
  register int i;
  int nb_value[SEQUENCE_NB_VARIABLE];
  double proba;
  Hidden_semi_markov *ihsmarkov , *hsmarkov;


  hsmarkov = 0;
  error.init();

  if ((nb_state < 2) || (nb_state > NB_STATE)) {
    status = false;
    error.update(SEQ_error[SEQR_NB_STATE]);
  }
  if ((occupancy_mean != D_DEFAULT) && (occupancy_mean <= 1.)) {
    status = false;
    error.update(SEQ_error[SEQR_OCCUPANCY]);
  }

  if (status) {
    for (i = 0;i < nb_variable;i++) {
      nb_value[i] = marginal[i]->nb_value;
    }

    ihsmarkov = new Hidden_semi_markov(type , nb_state , nb_variable , nb_value);

    // initialisation des parametres de la chaine de Markov

    ihsmarkov->init(left_right , 0.);

    // initialisation des lois d'occupations des etats

    if (occupancy_mean == D_DEFAULT) {
      occupancy_mean = MAX(hlength->mean , OCCUPANCY_MEAN);
    }

    ihsmarkov->state_subtype = new int[nb_state];
    ihsmarkov->nonparametric_process[0]->absorption = new double[nb_state];
    ihsmarkov->nonparametric_process[0]->sojourn_time = new Parametric*[nb_state];
    ihsmarkov->forward = new Forward*[nb_state];

    for (i = 0;i < nb_state;i++) {
      if (ihsmarkov->state_type[i] != 'a') {
        ihsmarkov->state_subtype[i] = SEMI_MARKOVIAN;
        ihsmarkov->nonparametric_process[0]->absorption[i] = 0.;
        proba = 1. / occupancy_mean;
        ihsmarkov->nonparametric_process[0]->sojourn_time[i] = new Parametric(NEGATIVE_BINOMIAL , 1 , I_DEFAULT ,
                                                                              1. , proba , OCCUPANCY_THRESHOLD);

        if (ihsmarkov->state_type[i] == 'r') {
          ihsmarkov->forward[i] = new Forward(*(ihsmarkov->nonparametric_process[0]->sojourn_time[i]) ,
                                              ihsmarkov->nonparametric_process[0]->sojourn_time[i]->alloc_nb_value);
        }
        else {
          ihsmarkov->forward[i] = 0;
        }
      }

      else {
        ihsmarkov->state_subtype[i] = MARKOVIAN;
        ihsmarkov->nonparametric_process[0]->absorption[i] = 1.;
        ihsmarkov->nonparametric_process[0]->sojourn_time[i] = 0;
        ihsmarkov->forward[i] = 0;
      }
    }

    // initialisation des lois d'observation

    for (i = 0;i < ihsmarkov->nb_output_process;i++) {
      if (ihsmarkov->nonparametric_process[i + 1]) {
        ihsmarkov->nonparametric_process[i + 1]->init();
      }
      else {
        ihsmarkov->parametric_process[i + 1]->init();
      }
    }

    hsmarkov = hidden_semi_markov_stochastic_estimation(error , os , *ihsmarkov , min_nb_state_sequence ,
                                                        max_nb_state_sequence , parameter , estimator ,
                                                        counting_flag , state_sequence , nb_iter);
    delete ihsmarkov;
  }

  return hsmarkov;
}


/*--------------------------------------------------------------*
 *
 *  Calcul des probabilites a posteriori des etats par l'algorithme forward-backward.
 *
 *  arguments : reference sur un objet Semi_markov_data, indice de la sequence,
 *              stream, type de sortie, format de fichier ('a' : ASCII,
 *              's' : Spreadsheet, 'g' : Gnuplot), references sur l'entropie marginale
 *              maximum et l'entropie (pour la visualisation).
 *
 *--------------------------------------------------------------*/

double Hidden_semi_markov::forward_backward(const Semi_markov_data &seq , int index ,
                                            ostream &os , int output , char format ,
                                            double &max_marginal_entropy , double &entropy1) const

{
  register int i , j , k , m;
  int *pstate , **poutput;
  double seq_likelihood , state_seq_likelihood , obs_product , entropy2 , buff , sum ,
         backward_max , **observation , *norm , *state_norm , **forward1 , **state_in ,
         **backward , **backward1 , *auxiliary , *occupancy_auxiliary , **backward_output ,
         *transition_predicted , *occupancy_predicted , **state_entropy , **predicted_entropy ,
         **transition_entropy , **occupancy_entropy , *partial_entropy , *conditional_entropy ,
         *marginal_entropy;
  Parametric *occupancy;


  // initialisations

  observation = new double*[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    observation[i] = new double[nb_state];
  }

  norm = new double[seq.length[index]];
  state_norm = new double[nb_state];

  forward1 = new double*[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    forward1[i] = new double[nb_state];
  }

  state_in = new double*[seq.length[index] - 1];
  for (i = 0;i < seq.length[index] - 1;i++) {
    state_in[i] = new double[nb_state];
  }

  backward = new double*[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    backward[i] = new double[nb_state];
  }

  backward1 = new double*[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    backward1[i] = new double[nb_state];
  }

  auxiliary = new double[nb_state];
  occupancy_auxiliary = new double[seq.length[index] + 1];

  if (output == SSTATE) {
    backward_output = backward;
  }
  else {
    backward_output = new double*[seq.length[index]];
    for (i = 0;i < seq.length[index];i++) {
      backward_output[i] = new double[nb_state];
    }
  }

  transition_predicted = new double[nb_state];
  occupancy_predicted = new double[seq.length[index] + 1];

  state_entropy = new double*[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    state_entropy[i] = new double[nb_state];
  }

  predicted_entropy = new double*[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    predicted_entropy[i] = new double[nb_state];
  }

  transition_entropy = new double*[nb_state];
  for (i = 0;i < nb_state;i++) {
    transition_entropy[i] = new double[nb_state];
  }

  occupancy_entropy = new double*[nb_state];
  for (i = 0;i < nb_state;i++) {
    switch (state_subtype[i]) {
    case SEMI_MARKOVIAN :
      occupancy = nonparametric_process[0]->sojourn_time[i];
      occupancy_entropy[i] = new double[MIN(seq.length[index] , occupancy->nb_value)];
      break;
    case MARKOVIAN :
      occupancy_entropy[i] = 0;
      break;
    }
  }

  partial_entropy = new double[seq.length[index]];
  conditional_entropy = new double[seq.length[index]];
  marginal_entropy = new double[seq.length[index]];

# ifdef DEBUG
  double *backward0;

  backward0 = new double[nb_state];
# endif

  poutput = new int*[nb_output_process];

  for (i = 0;i < nb_output_process;i++) {
    poutput[i] = seq.sequence[index][i + 1];
  }

  // recurrence "forward"

  seq_likelihood = 0.;
  for (i = 0;i < seq.length[index];i++) {
    norm[i] = 0.;

    for (j = 0;j < nb_state;j++) {

      // calcul des probabilites d'observation

      observation[i][j] = 1.;
      for (k = 0;k < nb_output_process;k++) {
        if (nonparametric_process[k + 1]) {
          observation[i][j] *= nonparametric_process[k + 1]->observation[j]->mass[*poutput[k]];
        }
        else {
          observation[i][j] *= parametric_process[k + 1]->observation[j]->mass[*poutput[k]];
        }
      }

      switch (state_subtype[j]) {

      // cas etat semi-markovien 

      case SEMI_MARKOVIAN : {
        if (i == 0) {
          state_norm[j] = initial[j];
        }
        else {
          state_norm[j] += state_in[i - 1][j] - forward1[i - 1][j];
        }
        state_norm[j] *= observation[i][j];

        norm[i] += state_norm[j];
        break;
      }

      // cas etat markovien

      case MARKOVIAN : {
        if (i == 0) {
          forward1[i][j] = initial[j];
          state_entropy[i][j] = 0.;
        }
        else {
          forward1[i][j] = state_in[i - 1][j];
          state_entropy[i][j] = predicted_entropy[i - 1][j];
        }
        forward1[i][j] *= observation[i][j];

        norm[i] += forward1[i][j];
        break;
      }
      }
    }

    if (norm[i] > 0.) {
      for (j = 0;j < nb_state;j++) {
        switch (state_subtype[j]) {
        case SEMI_MARKOVIAN :
          state_norm[j] /= norm[i];
          break;
        case MARKOVIAN :
          forward1[i][j] /= norm[i];
          break;
        }
      }

      seq_likelihood += log(norm[i]);
    }

    else {
      seq_likelihood = D_INF;
      break;
    }

    for (j = 0;j < nb_state;j++) {

      // cas etat semi-markovien

      if (state_subtype[j] == SEMI_MARKOVIAN) {
        occupancy = nonparametric_process[0]->sojourn_time[j];
        obs_product = 1.;
        forward1[i][j] = 0.;

        if (i < seq.length[index] - 1) {
          for (k = 1;k <= MIN(i + 1 , occupancy->nb_value - 1);k++) {
            obs_product *= observation[i - k + 1][j] / norm[i - k + 1];
            if (obs_product == 0.) {
              break;
            }

            if (k < i + 1) {
              occupancy_predicted[k] = obs_product * occupancy->mass[k] * state_in[i - k][j];
//              forward1[i][j] += obs_product * occupancy->mass[k] * state_in[i - k][j];
            }

            else {
              switch (type) {
              case 'o' :
                occupancy_predicted[k] = obs_product * occupancy->mass[k] * initial[j];
//                forward1[i][j] += obs_product * occupancy->mass[k] * initial[j];
                break;
              case 'e' :
                occupancy_predicted[k] = obs_product * forward[j]->mass[k] * initial[j];
//                forward1[i][j] += obs_product * forward[j]->mass[k] * initial[j];
                break;
              }
            }

            forward1[i][j] += occupancy_predicted[k];
          }
        }

        else {
          for (k = 1;k <= MIN(i + 1 , occupancy->nb_value - 1);k++) {
            obs_product *= observation[i - k + 1][j] / norm[i - k + 1];
            if (obs_product == 0.) {
              break;
            }

            if (k < i + 1) {
              occupancy_predicted[k] = obs_product * (1. - occupancy->cumul[k - 1]) * state_in[i - k][j];
//              forward1[i][j] += obs_product * (1. - occupancy->cumul[k - 1]) * state_in[i - k][j];
            }

            else {
              switch (type) {
              case 'o' :
                occupancy_predicted[k] = obs_product * (1. - occupancy->cumul[k - 1]) * initial[j];
//                forward1[i][j] += obs_product * (1. - occupancy->cumul[k - 1]) * initial[j];
                break;
              case 'e' :
                occupancy_predicted[k] = obs_product * (1. - forward[j]->cumul[k - 1]) * initial[j];
//                forward1[i][j] += obs_product * (1. - forward[j]->cumul[k - 1]) * initial[j];
                break;
              }
            }

            forward1[i][j] += occupancy_predicted[k];
          }
        }

        state_entropy[i][j] = 0.;

        if (forward1[i][j] > 0.) {
          for (m = 1;m < k;m++) {
            buff = occupancy_predicted[m] / forward1[i][j];
            if (buff > 0.) {
              if (m < i + 1) {
                state_entropy[i][j] += buff * (predicted_entropy[i - m][j] - log(buff));
              }
              else {
                state_entropy[i][j] -= buff * log(buff);
              }
            }
          }

          if (state_entropy[i][j] < 0.) {
            state_entropy[i][j] = 0.;
          }
        }
      }
    }

    if (i < seq.length[index] - 1) {
      for (j = 0;j < nb_state;j++) {
        state_in[i][j] = 0.;
        for (k = 0;k < nb_state;k++) {
          transition_predicted[k] = transition[k][j] * forward1[i][k];
          state_in[i][j] += transition_predicted[k];

//          state_in[i][j] += transition[k][j] * forward1[i][k];
        }

        predicted_entropy[i][j] = 0.;

        if (state_in[i][j] > 0.) {
          for (k = 0;k < nb_state;k++) {
            buff = transition_predicted[k] / state_in[i][j];
            if (buff > 0.) {
              predicted_entropy[i][j] += buff * (state_entropy[i][k] - log(buff));
            }
          }

          if (predicted_entropy[i][j] < 0.) {
            predicted_entropy[i][j] = 0.;
          }
        }
      }
    }

    for (j = 0;j < nb_output_process;j++) {
      poutput[j]++;
    }
  }

  if (seq_likelihood != D_INF) {
    entropy1 = 0.;
    i = seq.length[index] - 1;
    for (j = 0;j < nb_state;j++) {
      if (forward1[i][j] > 0.) {
        entropy1 += forward1[i][j] * (state_entropy[i][j] - log(forward1[i][j]));
      }
    }

    for (i = 0;i < nb_state;i++) {
      if (state_subtype[i] == SEMI_MARKOVIAN) {
        for (j = 0;j < seq.length[index];j++) {
          state_entropy[j][i] = 0.;
        }
      }
    }

    // recurrence "backward"

    for (i = 0;i < nb_output_process;i++) {
      poutput[i]--;
    }

    for (i = 0;i < nb_state;i++) {
      for (j = 0;j < nb_state;j++) {
        transition_entropy[i][j] = 0.;
      }
    }

    for (i = 0;i < nb_state;i++) {
      if (state_subtype[i] == SEMI_MARKOVIAN) {
        occupancy = nonparametric_process[0]->sojourn_time[i];
        for (j = occupancy->offset;j < MIN(seq.length[index] , occupancy->nb_value);j++) {
          occupancy_entropy[i][j] = 0.;
        }
      }
    }

    entropy2 = 0.;

    i = seq.length[index] - 1;
    for (j = 0;j < nb_state;j++) {
      backward[i][j] = forward1[i][j];
      backward1[i][j] = backward[i][j];

      if (output == OUT_STATE) {
        backward_output[i][j] = backward[i][j];
      }

      if (backward[i][j] > 0.) {
        for (k = 0;k < nb_output_process;k++) {
          if (nonparametric_process[k + 1]) {
            if (nonparametric_process[k + 1]->observation[j]->mass[*poutput[k]] > 0.) {
              entropy2 -= backward[i][j] * log(nonparametric_process[k + 1]->observation[j]->mass[*poutput[k]]);
            }
          }
          else {
            if (parametric_process[k + 1]->observation[j]->mass[*poutput[k]] > 0.) {
              entropy2 -= backward[i][j] * log(parametric_process[k + 1]->observation[j]->mass[*poutput[k]]);
            }
          }
        }
      }
    }

    for (i = seq.length[index] - 2;i >= 0;i--) {
      for (j = 0;j < nb_output_process;j++) {
        poutput[j]--;
      }

      for (j = 0;j < nb_state;j++) {
        auxiliary[j] = 0.;

        switch (state_subtype[j]) {

        // cas etat semi-markovien

        case SEMI_MARKOVIAN : {
          occupancy = nonparametric_process[0]->sojourn_time[j];
          obs_product = 1.;

          for (k = 1;k < MIN(seq.length[index] - i , occupancy->nb_value);k++) {
            obs_product *= observation[i + k][j] / norm[i + k];
            if (obs_product == 0.) {
              break;
            }

            occupancy_auxiliary[k] = 0.;

            if (backward1[i + k][j] > 0.) {
//            if (forward1[i + k][j] > 0.) {
              if (k < seq.length[index] - i - 1) {
                buff = backward1[i + k][j] * obs_product * occupancy->mass[k] /
                       forward1[i + k][j];
                occupancy_auxiliary[k] = buff * state_in[i][j];
                occupancy_entropy[j][k] += occupancy_auxiliary[k];

/*                if (occupancy->mass[k] > 0.) {
                  entropy2 -= occupancy_auxiliary[k] * log(occupancy->mass[k]);
                } */
              }

              else {
                buff = obs_product * (1. - occupancy->cumul[k - 1]);
                occupancy_auxiliary[k] = buff * state_in[i][j];
                if (occupancy->cumul[k - 1] < 1.) {
                  entropy2 -= occupancy_auxiliary[k] * log(1. - occupancy->cumul[k - 1]);
                }
              }

              auxiliary[j] += buff;
            }
          }

          sum = 0.;
          for (m = k - 1;m >= 1;m--) {
            sum += occupancy_auxiliary[m];
            if (backward[i + m][j] > 0.) {
              buff = sum / backward[i + m][j];
              if (buff > 0.) {
                state_entropy[i + m][j] += buff * (predicted_entropy[i][j] - log(buff));
              }
            }
          }
          break;
        }

        // cas etat markovien

        case MARKOVIAN : {
          if (backward1[i + 1][j] > 0.) {
//          if (forward1[i + 1][j] > 0.) {
            auxiliary[j] = backward1[i + 1][j] / state_in[i][j];

/*            auxiliary[j] = backward1[i + 1][j] * observation[i + 1][j] /
                           (forward1[i + 1][j] * norm[i + 1]); */

            state_entropy[i + 1][j] = predicted_entropy[i][j];
          }
          break;
        }
        }
      }

      for (j = 0;j < nb_state;j++) {
        backward1[i][j] = 0.;

        for (k = 0;k < nb_state;k++) {
          buff = auxiliary[k] * transition[j][k] * forward1[i][j];
          backward1[i][j] += buff;
          transition_entropy[j][k] += buff;

/*          if (transition[j][k] > 0.) {
            entropy2 -= buff * log(transition[j][k]);
          } */
        }

        switch (state_subtype[j]) {

        // cas etat semi-markovien

        case SEMI_MARKOVIAN : {
          backward[i][j] = backward[i + 1][j] + backward1[i][j] - auxiliary[j] * state_in[i][j];
          if (backward[i][j] < 0.) {
            backward[i][j] = 0.;
          }
          if (backward[i][j] > 1.) {
            backward[i][j] = 1.;
          }
          break;
        }

        // cas etat markovien

        case MARKOVIAN : {
          backward[i][j] = backward1[i][j];
          break;
        }
        }

        if (backward[i][j] > 0.) {
          for (k = 0;k < nb_output_process;k++) {
            if (nonparametric_process[k + 1]) {
              if (nonparametric_process[k + 1]->observation[j]->mass[*poutput[k]] > 0.) {
                entropy2 -= backward[i][j] * log(nonparametric_process[k + 1]->observation[j]->mass[*poutput[k]]);
              }
            }
            else {
              if (parametric_process[k + 1]->observation[j]->mass[*poutput[k]] > 0.) {
                entropy2 -= backward[i][j] * log(parametric_process[k + 1]->observation[j]->mass[*poutput[k]]);
              }
            }
          }
        }
      }

      switch (output) {

      case IN_STATE : {
        for (j = 0;j < nb_state;j++) {
          switch (state_subtype[j]) {

          // cas etat semi-markovien

          case SEMI_MARKOVIAN : {
            backward_output[i + 1][j] = auxiliary[j] * state_in[i][j];
            break;
          }

          // cas etat markovien

          case MARKOVIAN : {
            backward_output[i + 1][j] = 0.;
            for (k = 0;k < nb_state;k++) {
              if (k != j) {
                backward_output[i + 1][j] += transition[k][j] * forward1[i][k];
              }
            }
            backward_output[i + 1][j] *= auxiliary[j];
            break;
          }
          }
        }
        break;
      }

      case OUT_STATE : {
        for (j = 0;j < nb_state;j++) {
          switch (state_subtype[j]) {

          // cas etat semi-markovien

          case SEMI_MARKOVIAN : {
            backward_output[i][j] = backward1[i][j];
            break;
          }

          // cas etat markovien

          case MARKOVIAN : {
            backward_output[i][j] = 0.;
            for (k = 0;k < nb_state;k++) {
              if (k != j) {
                backward_output[i][j] += auxiliary[k] * transition[j][k];
              }
            }
            backward_output[i][j] *= forward1[i][j];
            break;
          }
          }
        }
        break;
      }
      }
    }

    if (output == IN_STATE) {
      for (i = 0;i < nb_state;i++) {
        backward_output[0][i] = backward[0][i];
      }
    }

    for (i = 0;i < nb_state;i++) {
      if (initial[i] > 0.) {
        entropy2 -= backward[0][i] * log(initial[i]);
      }
    }

    for (i = 0;i < nb_state;i++) {
      for (j = 0;j < nb_state;j++) {
        if (transition[i][j] > 0.) {
          entropy2 -= transition_entropy[i][j] * log(transition[i][j]);
        }
      }
    }

    for (i = 0;i < nb_state;i++) {
      if (state_subtype[i] == SEMI_MARKOVIAN) {
        occupancy = nonparametric_process[0]->sojourn_time[i];

        if (initial[i] > 0.) {
          obs_product = 1.;

#         ifdef DEBUG
          backward0[i] = 0.;
#         endif

          for (j = 1;j < MIN(seq.length[index] + 1 , occupancy->nb_value);j++) {
            obs_product *= observation[j - 1][i] / norm[j - 1];
            if (obs_product == 0.) {
              break;
            }

            occupancy_auxiliary[j] = 0.;

            if (backward1[j - 1][i] > 0.) {
//            if (forward1[j - 1][i] > 0.) {
              if (j < seq.length[index]) {
                switch (type) {

                case 'o' : {
                  occupancy_auxiliary[j] = backward1[j - 1][i] * obs_product * occupancy->mass[j] *
                                           initial[i] / forward1[j - 1][i];
                  occupancy_entropy[i][j] += occupancy_auxiliary[j];

/*                  if (occupancy->mass[j] > 0.) {
                    entropy2 -= occupancy_auxiliary[j] * log(occupancy->mass[j]);
                  } */
                  break;
                }

                case 'e' : {
                  occupancy_auxiliary[j] = backward1[j - 1][i] * obs_product * forward[i]->mass[j] *
                                           initial[i] / forward1[j - 1][i];
                  if (forward[i]->mass[j] > 0.) {
                    entropy2 -= occupancy_auxiliary[j] * log(forward[i]->mass[j]);
                  }
                  break;
                }
                }
              }

              else {
                switch (type) {

                case 'o' : {
                  occupancy_auxiliary[j] = obs_product * (1. - occupancy->cumul[j - 1]) * initial[i];
                  if (occupancy->cumul[j - 1] < 1.) {
                    entropy2 -= occupancy_auxiliary[j] * log(1. - occupancy->cumul[j - 1]);
                  }
                  break;
                }

                case 'e' : {
                  occupancy_auxiliary[j] = obs_product * (1. - forward[i]->cumul[j - 1]) * initial[i];
                  if (forward[i]->cumul[j - 1] < 1.) {
                    entropy2 -= occupancy_auxiliary[j] * log(1. - forward[i]->cumul[j - 1]);
                  }
                  break;
                }
                }
              }

#             ifdef DEBUG
              backward0[i] += occupancy_auxiliary[j];
#             endif

            }
          }

#         ifdef DEBUG
          cout << i << " " << backward[0][i] << " " << backward0[i] << endl;
#         endif

          sum = 0.;
          for (k = j - 1;k >= 1;k--) {
            sum += occupancy_auxiliary[k];
            if (backward[k - 1][i] > 0.) {
              buff = sum / backward[k - 1][i];
              if (buff > 0.) {
                state_entropy[k - 1][i] -= buff * log(buff);
              }
            }
          }
        }

        for (j = occupancy->offset;j < MIN(seq.length[index] , occupancy->nb_value);j++) {
          if (occupancy->mass[j] > 0.) {
            entropy2 -= occupancy_entropy[i][j] * log(occupancy->mass[j]);
          }
        }
      }
    }

    entropy2 += seq_likelihood;

#   ifdef MESSAGE
    if ((entropy2 < entropy1 - DOUBLE_ERROR) || (entropy2 > entropy1 + DOUBLE_ERROR)) {
      cout << "\nERROR: " << entropy1 << " " << entropy2 << endl;
    }
#   endif

    // restauration

    pstate = seq.sequence[index][0];

    for (i = 0;i < seq.length[index];i++) {
      backward_max = 0.;
      for (j = 0;j < nb_state;j++) {
        if (backward[i][j] > backward_max) {
          backward_max = backward[i][j];
          *pstate = j;
        }
      }

      pstate++;
    }

    state_seq_likelihood = Semi_markov::likelihood_computation(seq , index);

    for (i = 0;i < seq.length[index];i++) {
      partial_entropy[i] = 0.;
      for (j = 0;j < nb_state;j++) {
        if (state_entropy[i][j] < 0.) {
          state_entropy[i][j] = 0.;
        }
        if (backward[i][j] > 0.) {
          partial_entropy[i] += backward[i][j] * (state_entropy[i][j] - log(backward[i][j]));
        }
      }
      if (partial_entropy[i] < 0.) {
        partial_entropy[i] = 0.;
      }
    }

    conditional_entropy[0] = partial_entropy[0];
    for (i = 1;i < seq.length[index];i++) {
      conditional_entropy[i] = partial_entropy[i] - partial_entropy[i - 1];
    }

    max_marginal_entropy = 0.;
    for (i = 0;i < seq.length[index];i++) {
      marginal_entropy[i] = 0.;
      for (j = 0;j < nb_state;j++) {
        if (backward[i][j] > 0.) {
          marginal_entropy[i] -= backward[i][j] * log(backward[i][j]);
        }
      }
      if (marginal_entropy[i] > max_marginal_entropy) {
        max_marginal_entropy = marginal_entropy[i];
      }
    }

    switch (format) {

    case 'a' : {
      switch (output) {
      case SSTATE :
        os << "\n" << SEQ_label[SEQL_POSTERIOR_STATE_PROBABILITY] << "\n\n";
        break;
      case IN_STATE :
        os << "\n" << SEQ_label[SEQL_POSTERIOR_IN_STATE_PROBABILITY] << "\n\n";
        break;
      case OUT_STATE :
        os << "\n" << SEQ_label[SEQL_POSTERIOR_OUT_STATE_PROBABILITY] << "\n\n";
        break;
      }

//      seq.profile_ascii_print(os , index , nb_state , backward_output ,
//                              STAT_label[STATL_STATE]);
      seq.profile_ascii_print(os , index , nb_state , backward_output , conditional_entropy ,
                              marginal_entropy , partial_entropy);

      os << "\n" << STAT_label[STATL_LIKELIHOOD] << ": " << seq_likelihood
         << "\n" << SEQ_label[SEQL_STATE_SEQUENCE_LIKELIHOOD] << ": " << state_seq_likelihood
         << "   (" << exp(state_seq_likelihood - seq_likelihood) << ")" << endl;
      break;
    }

    case 's' : {
      switch (output) {
      case SSTATE :
        os << "\n" << SEQ_label[SEQL_POSTERIOR_STATE_PROBABILITY] << "\n\n";
        break;
      case IN_STATE :
        os << "\n" << SEQ_label[SEQL_POSTERIOR_IN_STATE_PROBABILITY] << "\n\n";
        break;
      case OUT_STATE :
        os << "\n" << SEQ_label[SEQL_POSTERIOR_OUT_STATE_PROBABILITY] << "\n\n";
        break;
      }

//      seq.profile_spreadsheet_print(os , index , nb_state , backward_output ,
//                                    STAT_label[STATL_STATE]);
      seq.profile_spreadsheet_print(os , index , nb_state , backward_output , conditional_entropy ,
                                    marginal_entropy , partial_entropy);

      os << "\n" << STAT_label[STATL_LIKELIHOOD] << "\t" << seq_likelihood
         << "\n" << SEQ_label[SEQL_STATE_SEQUENCE_LIKELIHOOD] << "\t" << state_seq_likelihood
         << "\t" << exp(state_seq_likelihood - seq_likelihood) << endl;
      break;
    }

    case 'g' : {
//      seq.profile_plot_print(os , index , nb_state , backward_output);
      seq.profile_plot_print(os , index , nb_state , backward_output , conditional_entropy ,
                             marginal_entropy , partial_entropy);
      break;
    }
    }

    if (format != 'g') {
/*      double gini_index;

      gini_index = 0.;
      for (i = 0;i < seq.length[index];i++) {
        for (j = 0;j < nb_state;j++) {
          gini_index += backward[i][j] * (1. - backward[i][j]);
        }
      } */

      double entropy3 , nb_state_sequence;

      entropy3 = 0.;
      for (i = 0;i < seq.length[index];i++) {
        for (j = 0;j < nb_state;j++) {
          if (backward[i][j] > 0.) {
            entropy3 -= backward[i][j] * log(backward[i][j]);
          }
        }
      }

      // calcul du nombre de sequences d'etats possibles

      for (i = 0;i < nb_output_process;i++) {
        poutput[i] = seq.sequence[index][i + 1];
      }

      // recurrence "forward"

      for (i = 0;i < seq.length[index];i++) {
        for (j = 0;j < nb_state;j++) {

          // calcul des probabilites d'observation

          observation[i][j] = 1.;
          for (k = 0;k < nb_output_process;k++) {
            if (nonparametric_process[k + 1]) {
              observation[i][j] *= nonparametric_process[k + 1]->observation[j]->mass[*poutput[k]];
            }
            else {
              observation[i][j] *= parametric_process[k + 1]->observation[j]->mass[*poutput[k]];
            }
          }

          if (observation[i][j] > 0.) {
            observation[i][j] = 1.;
          }

          forward1[i][j] = 0.;

          switch (state_subtype[j]) {

          // cas etat semi-markovien

          case SEMI_MARKOVIAN : {
            occupancy = nonparametric_process[0]->sojourn_time[j];

            if (i < seq.length[index] - 1) {
              for (k = 1;k <= MIN(i + 1 , occupancy->nb_value - 1);k++) {
                if (observation[i - k + 1][j] == 0.) {
                  break;
                }

                if (k < i + 1) {
                  if (occupancy->mass[k] > 0.) {
                    forward1[i][j] += state_in[i - k][j];
                  }
                }

                else {
                  if (initial[j] > 0.) {
                    switch (type) {

                    case 'o' : {
                      if (occupancy->mass[k] > 0.) {
                        forward1[i][j]++;
                      }
                      break;
                    }

                    case 'e' : {
                      if (forward[j]->mass[k] > 0.) {
                        forward1[i][j]++;
                      }
                      break;
                    }
                    }
                  }
                }
              }
            }

            else {
              for (k = 1;k <= MIN(i + 1 , occupancy->nb_value - 1);k++) {
                if (observation[i - k + 1][j] == 0.) {
                  break;
                }

                if (k < i + 1) {
                  if (1. - occupancy->cumul[k - 1] > 0.) {
                    forward1[i][j] += state_in[i - k][j];
                  }
                }

                else {
                  if (initial[j] > 0.) {
                    switch (type) {

                    case 'o' : {
                      if (1. - occupancy->cumul[k - 1] > 0.) {
                        forward1[i][j]++;
                      }
                      break;
                    }

                    case 'e' : {
                      if (1. - forward[j]->cumul[k - 1] > 0.) {
                        forward1[i][j]++;
                      }
                      break;
                    }
                    }
                  }
                }
              }
            }
            break;
          }

          // cas etat markovien

          case MARKOVIAN : {
            if (observation[i][j] == 1.) {
              if (i == 0) {
                if (initial[j] > 0.) {
                  forward1[i][j] = 1.;
                }
              }
              else {
                forward1[i][j] = state_in[i - 1][j];
              }
            }
            break;
          }
          }
        }

        if (i < seq.length[index] - 1) {
          for (j = 0;j < nb_state;j++) {
            state_in[i][j] = 0.;
            for (k = 0;k < nb_state;k++) {
              if (transition[k][j] > 0.) {
                state_in[i][j] += forward1[i][k];
              }
            }
          }
        }

        for (j = 0;j < nb_output_process;j++) {
          poutput[j]++;
        }
      }

      nb_state_sequence = 0.;
      i = seq.length[index] - 1;
      for (j = 0;j < nb_state;j++) {
        nb_state_sequence += forward1[i][j];
      }

      switch (format) {
      case 'a' :
/*        os << "\n" << SEQ_label[SEQL_GINI_INDEX] << ": " << gini_index << " ("
           << gini_index / seq.length[index] << ")   " << SEQ_label[SEQL_UPPER_BOUND] << ": "
           << seq.length[index] * (1. - 1. / nb_state) << " (" << 1. - 1. / nb_state
        os << ")   " << SEQ_label[SEQL_UPPER_BOUND] << ": "
           << seq.length[index] * log((double)nb_state) << " (" << log((double)nb_state) */
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
           << gini_index / seq.length[index] << "\t" << SEQ_label[SEQL_UPPER_BOUND] << "\t"
           << seq.length[index] * (1. - 1. / nb_state) << "\t" << 1. - 1. / nb_state
        os << "\t" << SEQ_label[SEQL_UPPER_BOUND] << "\t"
           << seq.length[index] * log((double)nb_state) << "\t" << log((double)nb_state) */
        os << "\n" << SEQ_label[SEQL_STATE_SEQUENCE_ENTROPY] << "\t" << entropy1
           << "\t" << entropy1 / seq.length[index] << "\t" << SEQ_label[SEQL_UPPER_BOUND] << "\t"
           << log((double)nb_state_sequence) << "\t"
           << log((double)nb_state_sequence) / seq.length[index]
           << "\n" << SEQ_label[SEQL_MARGINAL_ENTROPY_SUM] << "\t" << entropy3 << "\t"
           << entropy3 / seq.length[index] << "\n\n"
           << SEQ_label[SEQL_NB_STATE_SEQUENCE] << "\t" << nb_state_sequence << endl;
        break;
      }

#     ifdef DEBUG
      int state;
      double min_nb_state_sequence , smoothed_proba , cumul_smoothed_proba ,
             max_smoothed_proba , **backward2;

      // recurrence "backward"

      min_nb_state_sequence = nb_state_sequence;

      backward2 = new double*[seq.length[index]];
      for (i = 0;i < seq.length[index];i++) {
        backward2[i] = new double[nb_state];
      }

      i = seq.length[index] - 1;
      for (j = 0;j < nb_state;j++) {
        backward2[i][j] = forward1[i][j];
        backward1[i][j] = 1.;
      }

      for (i = seq.length[index] - 2;i >= 0;i--) {
        for (j = 0;j < nb_state;j++) {
          auxiliary[j] = 0.;

          switch (state_subtype[j]) {

          // cas etat semi-markovien

          case SEMI_MARKOVIAN : {
            occupancy = nonparametric_process[0]->sojourn_time[j];

            for (k = 1;k < MIN(seq.length[index] - i , occupancy->nb_value);k++) {
              if (observation[i + k][j] == 0.) {
                break;
              }

              if (k < seq.length[index] - i - 1) {
                if (occupancy->mass[k] > 0.) {
                  auxiliary[j] += backward1[i + k][j];
                }
              }
              else {
                if (1. - occupancy->cumul[k - 1] > 0.) {
                  auxiliary[j]++;
                }
              }
            }
            break;
          }

          // cas etat markovien

          case MARKOVIAN : {
            if (observation[i + 1][j] == 1.) {
              auxiliary[j] = backward1[i + 1][j];
            }
            break;
          }
          }
        }

        for (j = 0;j < nb_state;j++) {
          backward1[i][j] = 0.;

          for (k = 0;k < nb_state;k++) {
            if (transition[j][k] > 0.) {
              backward1[i][j] += auxiliary[k];
            }
          }

          switch (state_subtype[j]) {

          // cas etat semi-markovien

          case SEMI_MARKOVIAN : {

#           ifdef DEBUG
            if ((i == 0) && (initial[j] > 0.)) {
              occupancy = nonparametric_process[0]->sojourn_time[j];
              backward0[j] = 0.;

              for (k = 1;k < MIN(seq.length[index] + 1 , occupancy->nb_value);k++) {
                if (observation[k - 1][j] == 0.) {
                  break;
                }

                if (k < seq.length[index]) {
                  if (occupancy->mass[k] > 0.) {
                    backward0[j] += backward1[k - 1][j];
                  }
                }
                else {
                  if (1. - occupancy->cumul[k - 1] > 0.) {
                    backward0[j]++;
                  }
                }
              }
            }
#           endif

            backward2[i][j] = backward2[i + 1][j] + backward1[i][j] * forward1[i][j] -
                              auxiliary[j] * state_in[i][j];

#           ifdef DEBUG
            if ((i == 0) && (initial[j] > 0.)) {
              cout << j << " " << backward2[i][j] << " " << backward0[j] << endl;
            }
#           endif

            break;
          }

          // cas etat markovien

          case MARKOVIAN : {
            backward2[i][j] = backward1[i][j] * forward1[i][j];
            break;
          }
          }
        }

        smoothed_proba = 1.1;
        cumul_smoothed_proba = 0.;
        nb_state_sequence = 0;

        for (j = 0;j < nb_state;j++) {
          max_smoothed_proba = 0.;
          for (k = 0;k < nb_state;k++) {
            if ((backward[i][k] > max_smoothed_proba) && (backward[i][k] < smoothed_proba)) {
              max_smoothed_proba = backward[i][k];
              state = k;
            }
          }
          cumul_smoothed_proba += max_smoothed_proba;
          nb_state_sequence += backward2[i][state];

          if (cumul_smoothed_proba < 1. - MIN_SMOOTHED_PROBABILITY) {
            smoothed_proba = max_smoothed_proba;
          }
          else {
            break;
          }
        }

        if (nb_state_sequence < min_nb_state_sequence) {
          min_nb_state_sequence = nb_state_sequence;
        }
      }

      os << SEQ_label[SEQL_NB_STATE_SEQUENCE]
         << " (" << 1. - MIN_SMOOTHED_PROBABILITY << " beam)"
         << ": " << min_nb_state_sequence << endl;

      os << "\n";
      for (i = 0;i < seq.length[index];i++) {
        obs_product = 0.;
        for (j = 0;j < nb_state;j++) {
          os << backward2[i][j] << " (" << backward[i][j] << ")  ";
          obs_product += backward2[i][j];
        }
        os << "| " << obs_product << endl;
      }

      for (i = 0;i < seq.length[index];i++) {
        delete [] backward2[i];
      }
      delete [] backward2;
#     endif

    }
  }

  for (i = 0;i < seq.length[index];i++) {
    delete [] observation[i];
  }
  delete [] observation;

  delete [] norm;
  delete [] state_norm;

  for (i = 0;i < seq.length[index];i++) {
    delete [] forward1[i];
  }
  delete [] forward1;

  for (i = 0;i < seq.length[index] - 1;i++) {
    delete [] state_in[i];
  }
  delete [] state_in;

  for (i = 0;i < seq.length[index];i++) {
    delete [] backward[i];
  }
  delete [] backward;

  for (i = 0;i < seq.length[index];i++) {
    delete [] backward1[i];
  }
  delete [] backward1;

  delete [] auxiliary;
  delete [] occupancy_auxiliary;

  if (output != SSTATE) {
    for (i = 0;i < seq.length[index];i++) {
      delete [] backward_output[i];
    }
    delete [] backward_output;
  }

  delete [] transition_predicted;
  delete [] occupancy_predicted;

  for (i = 0;i < seq.length[index];i++) {
    delete [] state_entropy[i];
  }
  delete [] state_entropy;

  for (i = 0;i < seq.length[index];i++) {
    delete [] predicted_entropy[i];
  }
  delete [] predicted_entropy;

  for (i = 0;i < nb_state;i++) {
    delete [] transition_entropy[i];
  }
  delete [] transition_entropy;

  for (i = 0;i < nb_state;i++) {
    delete [] occupancy_entropy[i];
  }
  delete [] occupancy_entropy;

  delete [] partial_entropy;
  delete [] conditional_entropy;
  delete [] marginal_entropy;

# ifdef DEBUG
  delete [] backward0;
# endif

  delete [] poutput;

  return (seq_likelihood);
}


/*--------------------------------------------------------------*
 *
 *  Simulation de L sequences d'etats correspondant a une sequence observee
 *  par l'algorithme forward-backward.
 *
 *  arguments : reference sur un objet Semi_markov_data, indice de la sequence,
 *              stream, format de fichier ('a' : ASCII, 's' : Spreadsheet),
 *              nombre de sequences d'etats.
 *
 *--------------------------------------------------------------*/

double Hidden_semi_markov::forward_backward_sampling(const Semi_markov_data &seq , int index ,
                                                     ostream &os , char format ,
                                                     int nb_state_sequence) const

{
  register int i , j , k;
  int state_occupancy , *pstate , **poutput;
  double seq_likelihood , state_seq_likelihood , obs_product , **observation ,
         *norm , *state_norm , **forward1 , **state_in , *backward , *cumul_backward;
  Parametric *occupancy;

# ifdef DEBUG
  register int m;
  double sum;
# endif


  // initialisations

  observation = new double*[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    observation[i] = new double[nb_state];
  }

  norm = new double[seq.length[index]];
  state_norm = new double[nb_state];

  forward1 = new double*[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    forward1[i] = new double[nb_state];
  }

  state_in = new double*[seq.length[index] - 1];
  for (i = 0;i < seq.length[index] - 1;i++) {
    state_in[i] = new double[nb_state];
  }

  backward = new double[seq.length[index] + 1];
  cumul_backward = new double[seq.length[index] + 1];

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
  for (i = 0;i < seq.length[index];i++) {
    norm[i] = 0.;

    for (j = 0;j < nb_state;j++) {

      // calcul des probabilites d'observation

      observation[i][j] = 1.;
      for (k = 0;k < nb_output_process;k++) {
        if (nonparametric_process[k + 1]) {
          observation[i][j] *= nonparametric_process[k + 1]->observation[j]->mass[*poutput[k]];
        }
        else {
          observation[i][j] *= parametric_process[k + 1]->observation[j]->mass[*poutput[k]];
        }
      }

      switch (state_subtype[j]) {

      // cas etat semi-markovien 

      case SEMI_MARKOVIAN : {
        if (i == 0) {
          state_norm[j] = initial[j];
        }
        else {
          state_norm[j] += state_in[i - 1][j] - forward1[i - 1][j];
        }
        state_norm[j] *= observation[i][j];

        norm[i] += state_norm[j];
        break;
      }

      // cas etat markovien

      case MARKOVIAN : {
        if (i == 0) {
          forward1[i][j] = initial[j];
        }
        else {
          forward1[i][j] = state_in[i - 1][j];
        }
        forward1[i][j] *= observation[i][j];

        norm[i] += forward1[i][j];
        break;
      }
      }
    }

    if (norm[i] > 0.) {
      for (j = 0;j < nb_state;j++) {
        switch (state_subtype[j]) {
        case SEMI_MARKOVIAN :
          state_norm[j] /= norm[i];
          break;
        case MARKOVIAN :
          forward1[i][j] /= norm[i];
          break;
        }
      }

      seq_likelihood += log(norm[i]);
    }

    else {
      seq_likelihood = D_INF;
      break;
    }

    for (j = 0;j < nb_state;j++) {

      // cas etat semi-markovien

      if (state_subtype[j] == SEMI_MARKOVIAN) {
        occupancy = nonparametric_process[0]->sojourn_time[j];
        obs_product = 1.;
        forward1[i][j] = 0.;

        if (i < seq.length[index] - 1) {
          for (k = 1;k <= MIN(i + 1 , occupancy->nb_value - 1);k++) {
            obs_product *= observation[i - k + 1][j] / norm[i - k + 1];
            if (obs_product == 0.) {
              break;
            }

            if (k < i + 1) {
              forward1[i][j] += obs_product * occupancy->mass[k] * state_in[i - k][j];
            }

            else {
              switch (type) {
              case 'o' :
                forward1[i][j] += obs_product * occupancy->mass[k] * initial[j];
                break;
              case 'e' :
                forward1[i][j] += obs_product * forward[j]->mass[k] * initial[j];
                break;
              }
            }
          }
        }

        else {
          for (k = 1;k <= MIN(i + 1 , occupancy->nb_value - 1);k++) {
            obs_product *= observation[i - k + 1][j] / norm[i - k + 1];
            if (obs_product == 0.) {
              break;
            }

            if (k < i + 1) {
              forward1[i][j] += obs_product * (1. - occupancy->cumul[k - 1]) * state_in[i - k][j];
            }

            else {
              switch (type) {
              case 'o' :
                forward1[i][j] += obs_product * (1. - occupancy->cumul[k - 1]) * initial[j];
                break;
              case 'e' :
                forward1[i][j] += obs_product * (1. - forward[j]->cumul[k - 1]) * initial[j];
                break;
              }
            }
          }
        }
      }
    }

    if (i < seq.length[index] - 1) {
      for (j = 0;j < nb_state;j++) {
        state_in[i][j] = 0.;
        for (k = 0;k < nb_state;k++) {
          state_in[i][j] += transition[k][j] * forward1[i][k];
        }
      }
    }

    for (j = 0;j < nb_output_process;j++) {
      poutput[j]++;
    }
  }

  if (seq_likelihood != D_INF) {

    // passes "backward"

#   ifdef MESSAGE
    cout << "\n";
#   endif

    for (i = 0;i < nb_state_sequence;i++) {
      j = seq.length[index] - 1;
      pstate = seq.sequence[index][0] + j;
      ::cumul_computation(nb_state , forward1[j] , cumul_backward);
      *pstate = cumul_method(nb_state , cumul_backward);

      do {

        // cas etat semi-markovien

        if (state_subtype[*pstate] == SEMI_MARKOVIAN) {
          occupancy = nonparametric_process[0]->sojourn_time[*pstate];
          obs_product = 1.;

          if (j < seq.length[index] - 1) {
            for (k = 1;k <= MIN(j + 1 , occupancy->nb_value - 1);k++) {
              obs_product *= observation[j - k + 1][*pstate] / norm[j - k + 1];
              if (obs_product == 0.) {
                break;
              }

              if (k < j + 1) {
                backward[k] = obs_product * occupancy->mass[k] * state_in[j - k][*pstate] /
                              forward1[j][*pstate];
              }

              else {
                switch (type) {
                case 'o' :
                  backward[k] = obs_product * occupancy->mass[k] * initial[*pstate] /
                                forward1[j][*pstate];
                  break;
                case 'e' :
                  backward[k] = obs_product * forward[*pstate]->mass[k] * initial[*pstate] /
                                forward1[j][*pstate];
                  break;
                }
              }
            }
          }

          else {
            for (k = 1;k <= MIN(j + 1 , occupancy->nb_value - 1);k++) {
              obs_product *= observation[j - k + 1][*pstate] / norm[j - k + 1];
              if (obs_product == 0.) {
                break;
              }

              if (k < j + 1) {
                backward[k] = obs_product * (1. - occupancy->cumul[k - 1]) * state_in[j - k][*pstate] /
                              forward1[j][*pstate];
              }

              else {
                switch (type) {
                case 'o' :
                  backward[k] = obs_product * (1. - occupancy->cumul[k - 1]) * initial[*pstate] /
                                forward1[j][*pstate];
                  break;
                case 'e' :
                  backward[k] = obs_product * (1. - forward[*pstate]->cumul[k - 1]) * initial[*pstate] /
                                forward1[j][*pstate];
                  break;
                }
              }
            }
          }

          ::cumul_computation(k - 1 , backward + 1 , cumul_backward);
          state_occupancy = 1 + cumul_method(k - 1 , cumul_backward);

#         ifdef DEBUG
          sum = 0.;
          for (m = 1;m < k;m++) {
            sum += backward[m];
          }
          if ((sum < 1. - DOUBLE_ERROR) || (sum > 1. + DOUBLE_ERROR)) {
            cout << "\nERROR: " << j << " " << sum << endl;
          }
#         endif

          for (k = 1;k < state_occupancy;k++) {
            pstate--;
            *pstate = *(pstate + 1);
          }
          j -= (state_occupancy - 1);

          if (j == 0) {
            break;
          }
        }

        j--;
        for (k = 0;k < nb_state;k++) {
          backward[k] = transition[k][*pstate] * forward1[j][k] / state_in[j][*pstate];
        }
        ::cumul_computation(nb_state , backward , cumul_backward);
        *--pstate = cumul_method(nb_state , cumul_backward);

#       ifdef DEBUG
        sum = 0.;
        for (k = 0;k < nb_state;k++) {
          sum += backward[k];
        }
        if ((sum < 1. - DOUBLE_ERROR) || (sum > 1. + DOUBLE_ERROR)) {
          cout << "\nERROR: " << j << " " << sum << endl;
        }
#       endif

      }
      while (j > 0);

#     ifdef DEBUG
      pstate = seq.sequence[index][0];
      for (j = 0;j < seq.length[index];j++) {
        state_sequence_probability[j][*pstate++]++;
      }
#     endif

#     ifdef MESSAGE
      state_seq_likelihood = Semi_markov::likelihood_computation(seq , index);

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
    delete [] observation[i];
  }
  delete [] observation;

  delete [] norm;
  delete [] state_norm;

  for (i = 0;i < seq.length[index];i++) {
    delete [] forward1[i];
  }
  delete [] forward1;

  for (i = 0;i < seq.length[index] - 1;i++) {
    delete [] state_in[i];
  }
  delete [] state_in;

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
 *              identificateur de la sequence, type de sortie,
 *              format ('a' : ASCII, 's' : Spreadsheet),
 *              methode de calcul des sequences d'etats (algorithme de Viterbi generalise ou
 *              algorithme forward-backward de simulation), nombre de sequences d'etats.
 *
 *--------------------------------------------------------------*/

bool Hidden_semi_markov::state_profile_write(Format_error &error , ostream &os ,
                                             const Semi_markov_data &iseq , int identifier ,
                                             int output , char format , int state_sequence ,
                                             int nb_state_sequence) const

{
  bool status = true;
  register int i;
  int index = I_DEFAULT;
  double seq_likelihood , max_marginal_entropy , entropy;
  Hidden_semi_markov *hsmarkov1 , *hsmarkov2;
  Semi_markov_data *seq;


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
      seq = new Semi_markov_data((Markovian_sequences&)iseq , 0);
      seq->type[0] = STATE;
    }
    else {
      seq = new Semi_markov_data(iseq , false);
    }

    hsmarkov1 = new Hidden_semi_markov(*this , false , false);

    hsmarkov2 = new Hidden_semi_markov(*this , false , false);
    hsmarkov2->create_cumul();
    hsmarkov2->log_computation();

    for (i = 0;i < seq->nb_sequence;i++) {
      if ((index == I_DEFAULT) || (index == i)) {
        seq_likelihood = hsmarkov1->forward_backward(*seq , i , os , output , format ,
                                                     max_marginal_entropy , entropy);
        hsmarkov2->viterbi_forward_backward(*seq , i , os , output , format , seq_likelihood);

        switch (state_sequence) {
        case GENERALIZED_VITERBI :
          hsmarkov2->generalized_viterbi(*seq , i , os , seq_likelihood , format ,
                                         nb_state_sequence);
          break;
        case FORWARD_BACKWARD_SAMPLING :
          hsmarkov1->forward_backward_sampling(*seq , i , os , format ,
                                               nb_state_sequence);
          break;
        }
      }
    }

    delete seq;

    delete hsmarkov1;
    delete hsmarkov2;
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
 *  arguments : reference sur un objet Format_error, stream,
 *              identificateur de la sequence, type de sortie,
 *              methode de calcul des sequences d'etats (algorithme de Viterbi generalise ou
 *              algorithme forward-backward de simulation), nombre de sequences d'etats.
 *
 *--------------------------------------------------------------*/

bool Hidden_semi_markov::state_profile_ascii_write(Format_error &error , ostream &os ,
                                                   int identifier , int output ,
                                                   int state_sequence , int nb_state_sequence) const

{
  bool status;


  error.init();

  if (!semi_markov_data) {
    status = false;
    error.update(STAT_error[STATR_NO_DATA]);
  }
  else {
    status = state_profile_write(error , os , *semi_markov_data , identifier ,
                                 output , 'a' , state_sequence , nb_state_sequence);
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
 *              type de sortie, format de fichier ('a' : ASCII, 's' : Spreadsheet),
 *              methode de calcul des sequences d'etats (algorithme de Viterbi generalise ou
 *              algorithme forward-backward de simulation), nombre de sequences d'etats.
 *
 *--------------------------------------------------------------*/

bool Hidden_semi_markov::state_profile_write(Format_error &error , const char *path ,
                                             int identifier , int output , char format ,
                                             int state_sequence , int nb_state_sequence) const

{
  bool status = true;
  ofstream out_file(path);


  error.init();

  if (!out_file) {
    status = false;
    error.update(STAT_error[STATR_FILE_NAME]);
  }
  if (!semi_markov_data) {
    status = false;
    error.update(STAT_error[STATR_NO_DATA]);
  }

  if (status) {
    status = state_profile_write(error , out_file , *semi_markov_data , identifier ,
                                 output , format , state_sequence , nb_state_sequence);
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Calcul des probabilites a posteriori des etats par l'algorithme forward-backward,
 *  des vraisemblances des sequences d'etats optimales par l'algorithme de Viterbi
 *  forward-backward et affichage des resultats au format Gnuplot.
 *
 *  arguments : reference sur un objet Format_error, prefixe des fichiers, sequences,
 *              identificateur de la sequence, type de sortie, titre des figures.
 *
 *--------------------------------------------------------------*/

bool Hidden_semi_markov::state_profile_plot_write(Format_error &error , const char *prefix ,
                                                  const Semi_markov_data &iseq , int identifier ,
                                                  int output , const char *title) const

{
  bool status = true;
  register int i , j;
  int index;
  double seq_likelihood , max_marginal_entropy , entropy , state_seq_likelihood;
  Hidden_semi_markov *hsmarkov;
  Semi_markov_data *seq;
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
        seq = new Semi_markov_data((Markovian_sequences&)iseq , 0);
        seq->type[0] = STATE;
      }
      else {
        seq = new Semi_markov_data(iseq , false);
      }

      hsmarkov = new Hidden_semi_markov(*this , false , false);

      seq_likelihood = hsmarkov->forward_backward(*seq , index , *data_out_file , output , 'g' ,
                                                  max_marginal_entropy , entropy);
      data_out_file->close();
      delete data_out_file;

      data_file_name[1] << prefix << 1 << ".dat";
      data_out_file = new ofstream((data_file_name[1].str()).c_str());

      hsmarkov->create_cumul();
      hsmarkov->log_computation();
      state_seq_likelihood = hsmarkov->viterbi_forward_backward(*seq , index , *data_out_file ,
                                                                output , 'g' , seq_likelihood);
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
        switch (output) {
        case SSTATE :
          out_file << SEQ_label[SEQL_POSTERIOR_STATE_PROBABILITY] << "\"\n\n";
          break;
        case IN_STATE :
          out_file << SEQ_label[SEQL_POSTERIOR_IN_STATE_PROBABILITY] << "\"\n\n";
          break;
        case OUT_STATE :
          out_file << SEQ_label[SEQL_POSTERIOR_OUT_STATE_PROBABILITY] << "\"\n\n";
          break;
        }

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
        switch (output) {
        case SSTATE :
          out_file << SEQ_label[SEQL_MAX_POSTERIOR_STATE_PROBABILITY] << "\"\n\n";
          break;
        case IN_STATE :
          out_file << SEQ_label[SEQL_MAX_POSTERIOR_IN_STATE_PROBABILITY] << "\"\n\n";
          break;
        case OUT_STATE :
          out_file << SEQ_label[SEQL_MAX_POSTERIOR_OUT_STATE_PROBABILITY] << "\"\n\n";
          break;
        }

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
                 << nb_state + 2 << " title \"" << SEQ_label[SEQL_MARGINAL_ENTROPY]
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
                 << nb_state + 3 << " title \"" << SEQ_label[SEQL_PARTIAL_STATE_SEQUENCE_ENTROPY]
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
      delete hsmarkov;
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
 *              identificateur de la sequence, type de sortie, titre des figures.
 *
 *--------------------------------------------------------------*/

bool Hidden_semi_markov::state_profile_plot_write(Format_error &error , const char *prefix ,
                                                  int identifier , int output , const char *title) const

{
  bool status;


  error.init();

  if (!semi_markov_data) {
    status = false;
    error.update(STAT_error[STATR_NO_DATA]);
  }
  else {
    status = state_profile_plot_write(error , prefix , *semi_markov_data , identifier ,
                                      output , title);
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Calcul des logarithmes des parametres d'une semi-chaine de Markov cachee.
 *
 *--------------------------------------------------------------*/

void Hidden_semi_markov::log_computation()

{
  register int i , j;
  double *pcumul;
  Parametric *occupancy;


  Chain::log_computation();

  for (i = 0;i < nb_state;i++) {
    if (state_subtype[i] == SEMI_MARKOVIAN) {
      occupancy = nonparametric_process[0]->sojourn_time[i];

      if (occupancy->mass[occupancy->offset] > 0.) {
        ::log_computation(occupancy->nb_value , occupancy->mass , occupancy->mass);

        pcumul = occupancy->cumul;
        for (j = 0;j < occupancy->nb_value;j++) {
          *pcumul = 1. - *pcumul;
          pcumul++;
        }
        ::log_computation(occupancy->nb_value , occupancy->cumul , occupancy->cumul);

        if (type == 'e') {
          ::log_computation(forward[i]->nb_value , forward[i]->mass , forward[i]->mass);

          pcumul = forward[i]->cumul;
          for (j = 0;j < forward[i]->nb_value;j++) {
            *pcumul = 1. - *pcumul;
            pcumul++;
          }
          ::log_computation(forward[i]->nb_value , forward[i]->cumul , forward[i]->cumul);
        }
      }
    }
  }

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
 *  Calcul des sequences d'etats optimales par l'algorithme de Viterbi.
 *
 *  arguments : reference sur un objet Semi_markov_data, indice de la sequence.
 *
 *--------------------------------------------------------------*/

double Hidden_semi_markov::viterbi(const Semi_markov_data &seq , int index) const

{
  register int i , j , k , m;
  int length , *pstate , **poutput , **input_state , **optimal_state ,
      **optimal_occupancy;
  double likelihood = 0. , obs_product , buff , forward_max , **observation ,
         *forward1 , **state_in;
  Parametric *occupancy;


  // initialisations

  length = (index == I_DEFAULT ? seq.max_length : seq.length[index]);

  observation = new double*[length];
  for (i = 0;i < length;i++) {
    observation[i] = new double[nb_state];
  }

  forward1 = new double[nb_state];

  state_in = new double*[length - 1];
  for (i = 0;i < length - 1;i++) {
    state_in[i] = new double[nb_state];
  }

  input_state = new int*[length - 1];
  for (i = 0;i < length - 1;i++) {
    input_state[i] = new int[nb_state];
  }

  optimal_state = new int*[length];
  for (i = 0;i < length;i++) {
    optimal_state[i] = new int[nb_state];
  }

  optimal_occupancy = new int*[length];
  for (i = 0;i < length;i++) {
    optimal_occupancy[i] = new int[nb_state];
  }

  poutput = new int*[nb_output_process];

  for (i = 0;i < seq.nb_sequence;i++) {
    if ((index == I_DEFAULT) || (index == i)) {
      for (j = 0;j < nb_output_process;j++) {
        poutput[j] = seq.sequence[i][j + 1];
      }

      // recurrence "forward"

      for (j = 0;j < seq.length[i];j++) {
        for (k = 0;k < nb_state;k++) {

          // calcul des probabilites d'observation

          observation[j][k] = 0.;
          for (m = 0;m < nb_output_process;m++) {
            if (nonparametric_process[m + 1]) {
              buff = nonparametric_process[m + 1]->observation[k]->cumul[*poutput[m]];
            }
            else {
              buff = parametric_process[m + 1]->observation[k]->cumul[*poutput[m]];
            }

            if (buff == D_INF) {
              observation[j][k] = D_INF;
              break;
            }
            else {
              observation[j][k] += buff;
            }
          }

          switch (state_subtype[k]) {

          // cas etat semi-markovien

          case SEMI_MARKOVIAN : {
            occupancy = nonparametric_process[0]->sojourn_time[k];
            obs_product = 0.;
            forward1[k] = D_INF;

            if (j < seq.length[i] - 1) {
              for (m = 1;m <= MIN(j + 1 , occupancy->nb_value - 1);m++) {
                if (observation[j - m + 1][k] == D_INF) {
                  break;
                }
                else {
                  obs_product += observation[j - m + 1][k];
                }

                if (m < j + 1) {
                  buff = obs_product + occupancy->mass[m] + state_in[j - m][k];
                }

                else {
                  switch (type) {
                  case 'o' :
                    buff = obs_product + occupancy->mass[m] + cumul_initial[k];
                    break;
                  case 'e' :
                    buff = obs_product + forward[k]->mass[m] + cumul_initial[k];
                    break;
                  }
                }

                if (buff > forward1[k]) {
                  forward1[k] = buff;
                  if (m < j + 1) {
                    optimal_state[j][k] = input_state[j - m][k];
                  }
                  optimal_occupancy[j][k] = m;
                }
              }
            }

            else {
              for (m = 1;m <= MIN(j + 1 , occupancy->nb_value - 1);m++) {
                if (observation[j - m + 1][k] == D_INF) {
                  break;
                }
                else {
                  obs_product += observation[j - m + 1][k];
                }

                if (m < j + 1) {
                  buff = obs_product + occupancy->cumul[m - 1] + state_in[j - m][k];
                }

                else {
                  switch (type) {
                  case 'o' :
                    buff = obs_product + occupancy->cumul[m - 1] + cumul_initial[k];
                    break;
                  case 'e' :
                    buff = obs_product + forward[k]->cumul[m - 1] + cumul_initial[k];
                    break;
                  }
                }

                if (buff > forward1[k]) {
                  forward1[k] = buff;
                  if (m < j + 1) {
                    optimal_state[j][k] = input_state[j - m][k];
                  }
                  optimal_occupancy[j][k] = m;
                }
              }
            }
            break;
          }

          // cas etat markovien

          case MARKOVIAN : {
            if (j == 0) {
              forward1[k] = cumul_initial[k];
            }
            else {
              forward1[k] = state_in[j - 1][k];
              optimal_state[j][k] = input_state[j - 1][k];
            }
            optimal_occupancy[j][k] = 1;

            if (forward1[k] != D_INF) {
              if (observation[j][k] == D_INF) {
                forward1[k] = D_INF;
              }
              else {
                forward1[k] += observation[j][k];
              }
            }
            break;
          }
          }
        }

#       ifdef DEBUG
        cout << j << " : ";
        for (k = 0;k < nb_state;k++) {
          cout << forward1[k];
          if (forward1[k] != D_INF) {
            cout << " " << optimal_occupancy[j][k] << " " << optimal_state[j][k];
          }
          cout << " | ";
        }
        cout << endl;
#       endif

        if (j < seq.length[i] - 1) {
          for (k = 0;k < nb_state;k++) {
            state_in[j][k] = D_INF;
            for (m = 0;m < nb_state;m++) {
              buff = cumul_transition[m][k] + forward1[m];
              if (buff > state_in[j][k]) {
                state_in[j][k] = buff;
                input_state[j][k] = m;
              }
            }
          }
        }

        for (k = 0;k < nb_output_process;k++) {
          poutput[k]++;
        }
      }

      // extraction de la vraisemblance du chemin optimal

      pstate = seq.sequence[i][0] + seq.length[i] - 1;
      forward_max = D_INF;

      for (j = 0;j < nb_state;j++) {
        if (forward1[j] > forward_max) {
          forward_max = forward1[j];
          *pstate = j;
        }
      }

      if (forward_max != D_INF) {
        likelihood += forward_max;
      }
      else {
        likelihood = D_INF;
        break;
      }

      // restauration

      j = seq.length[i] - 1;

      do {
        for (k = 0;k < optimal_occupancy[j][*pstate] - 1;k++) {
          pstate--;
          *pstate = *(pstate + 1);
        }

        if (j >= optimal_occupancy[j][*pstate]) {
          pstate--;
          *pstate = optimal_state[j][*(pstate + 1)];
          j -= optimal_occupancy[j][*(pstate + 1)];
        }
        else {
          j -= optimal_occupancy[j][*pstate];
        }
      }
      while (j >= 0);
    }
  }

  for (i = 0;i < length;i++) {
    delete [] observation[i];
  }
  delete [] observation;

  delete [] forward1;

  for (i = 0;i < length - 1;i++) {
    delete [] state_in[i];
  }
  delete [] state_in;

  for (i = 0;i < length - 1;i++) {
    delete [] input_state[i];
  }
  delete [] input_state;

  for (i = 0;i < length;i++) {
    delete [] optimal_state[i];
  }
  delete [] optimal_state;

  for (i = 0;i < length;i++) {
    delete [] optimal_occupancy[i];
  }
  delete [] optimal_occupancy;

  delete [] poutput;

  return likelihood;
}


/*--------------------------------------------------------------*
 *
 *  Calcul des L sequences d'etats optimales par l'algorithme de Viterbi generalise.
 *
 *  arguments : reference sur un objet Semi_markov_data, indice de la sequence,
 *              stream, vraisemblance des donnees, format de fichier
 *              ('a' : ASCII, 's' : Spreadsheet), nombre de sequences d'etats.
 *
 *--------------------------------------------------------------*/

double Hidden_semi_markov::generalized_viterbi(const Semi_markov_data &seq , int index ,
                                               ostream &os , double seq_likelihood ,
                                               char format , int inb_state_sequence) const

{
  bool **active_cell;
  register int i , j , k , m;
  int nb_state_sequence , max_occupancy , brank , previous_rank , nb_cell , *rank ,
      *pstate , **poutput , ***input_state , ***optimal_state , ***optimal_occupancy ,
      ***input_rank , ***optimal_rank;
  double buff , forward_max , state_seq_likelihood , likelihood_cumul , *obs_product ,
         **observation , **forward1 , ***state_in;
  Parametric *occupancy;


  // initialisations

  observation = new double*[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    observation[i] = new double[nb_state];
  }

  obs_product = new double[seq.length[index] + 1];

  forward1 = new double*[nb_state];
  for (i = 0;i < nb_state;i++) {
    forward1[i] = new double[inb_state_sequence];
  }

  state_in = new double**[seq.length[index] - 1];
  for (i = 0;i < seq.length[index] - 1;i++) {
    state_in[i] = new double*[nb_state];
    for (j = 0;j < nb_state;j++) {
      state_in[i][j] = new double[inb_state_sequence];
    }
  }

  rank = new int[MAX(seq.length[index] + 1 , nb_state)];

  input_state = new int**[seq.length[index] - 1];
  for (i = 0;i < seq.length[index] - 1;i++) {
    input_state[i] = new int*[nb_state];
    for (j = 0;j < nb_state;j++) {
      input_state[i][j] = new int[inb_state_sequence];
    }
  }

  optimal_state = new int**[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    optimal_state[i] = new int*[nb_state];
    for (j = 0;j < nb_state;j++) {
      optimal_state[i][j] = new int[inb_state_sequence];
    }
  }

  optimal_occupancy = new int**[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    optimal_occupancy[i] = new int*[nb_state];
    for (j = 0;j < nb_state;j++) {
      optimal_occupancy[i][j] = new int[inb_state_sequence];
    }
  }

  input_rank = new int**[seq.length[index] - 1];
  for (i = 0;i < seq.length[index] - 1;i++) {
    input_rank[i] = new int*[nb_state];
    for (j = 0;j < nb_state;j++) {
      input_rank[i][j] = new int[inb_state_sequence];
    }
  }

  optimal_rank = new int**[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    optimal_rank[i] = new int*[nb_state];
    for (j = 0;j < nb_state;j++) {
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
  nb_state_sequence = 1;

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

  for (i = 0;i < seq.length[index];i++) {
    for (j = 0;j < nb_state;j++) {

      // calcul des probabilites d'observation

      observation[i][j] = 0.;
      for (k = 0;k < nb_output_process;k++) {
        if (nonparametric_process[k + 1]) {
          buff = nonparametric_process[k + 1]->observation[j]->cumul[*poutput[k]];
        }
        else {
          buff = parametric_process[k + 1]->observation[j]->cumul[*poutput[k]];
        }

        if (buff == D_INF) {
          observation[i][j] = D_INF;
          break;
        }
        else {
          observation[i][j] += buff;
        }
      }

      switch (state_subtype[j]) {

      // cas etat semi-markovien

      case SEMI_MARKOVIAN : {
        occupancy = nonparametric_process[0]->sojourn_time[j];

        obs_product[0] = 0.;
        for (k = 1;k <= MIN(i + 1 , occupancy->nb_value - 1);k++) {
          if (observation[i - k + 1][j] == D_INF) {
            break;
          }
          else {
            obs_product[k] = obs_product[k - 1] + observation[i - k + 1][j];
          }
        }
        max_occupancy = k - 1;

        for (k = 1;k <= max_occupancy;k++) {
          rank[k] = 0;
        }

        for (k = 0;k < nb_state_sequence;k++) {
          forward1[j][k] = D_INF;

          if (i < seq.length[index] - 1) {
            for (m = 1;m <= max_occupancy;m++) {
              if (m < i + 1) {
                buff = obs_product[m] + occupancy->mass[m] + state_in[i - m][j][rank[m]];
              }

              else {
                if (rank[i + 1] == 0) {
                  switch (type) {
                  case 'o' :
                    buff = obs_product[m] + occupancy->mass[m] + cumul_initial[j];
                    break;
                  case 'e' :
                    buff = obs_product[m] + forward[j]->mass[m] + cumul_initial[j];
                    break;
                  }
                }

                else {
                  buff = D_INF;
                }
              }

              if (buff > forward1[j][k]) {
                forward1[j][k] = buff;
                if (m < i + 1) {
                  optimal_state[i][j][k] = input_state[i - m][j][rank[m]];
                  optimal_rank[i][j][k] = input_rank[i - m][j][rank[m]];
                }
                optimal_occupancy[i][j][k] = m;
              }
            }
          }

          else {
            for (m = 1;m <= max_occupancy;m++) {
              if (m < i + 1) {
                buff = obs_product[m] + occupancy->cumul[m - 1] + state_in[i - m][j][rank[m]];
              }

              else {
                if (rank[i + 1] == 0) {
                  switch (type) {
                  case 'o' :
                    buff = obs_product[m] + occupancy->cumul[m - 1] + cumul_initial[j];
                    break;
                  case 'e' :
                    buff = obs_product[m] + forward[j]->cumul[m - 1] + cumul_initial[j];
                    break;
                  }
                }

                else {
                  buff = D_INF;
                }
              }

              if (buff > forward1[j][k]) {
                forward1[j][k] = buff;
                if (m < i + 1) {
                  optimal_state[i][j][k] = input_state[i - m][j][rank[m]];
                  optimal_rank[i][j][k] = input_rank[i - m][j][rank[m]];
                }
                optimal_occupancy[i][j][k] = m;
              }
            }
          }

          if (forward1[j][k] != D_INF) {
            rank[optimal_occupancy[i][j][k]]++;
          }
        }
        break;
      }

      // cas etat markovien

      case MARKOVIAN : {
        for (k = 0;k < nb_state_sequence;k++) {
          if (i == 0) {
            forward1[j][k] = cumul_initial[j];
          }
          else {
            forward1[j][k] = state_in[i - 1][j][k];
            optimal_state[i][j][k] = input_state[i - 1][j][k];
            optimal_rank[i][j][k] = input_rank[i - 1][j][k];
          }
          optimal_occupancy[i][j][k] = 1;

          if (forward1[j][k] != D_INF) {
            if (observation[i][j] == D_INF) {
              forward1[j][k] = D_INF;
            }
            else {
              forward1[j][k] += observation[i][j];
            }
          }
        }
        break;
      }
      }

      for (k = nb_state_sequence;k < inb_state_sequence;k++) {
        forward1[j][k] = D_INF;
      }
    }

#   ifdef DEBUG
    cout << i << " : ";
    for (j = 0;j < nb_state;j++) {
      cout << j << " :";
      for (k = 0;k < nb_state_sequence;k++) {
        cout << " " << forward1[j][k];
        if (forward1[j][k] != D_INF) {
          cout << " " << optimal_occupancy[i][j][k];
          if (optimal_occupancy[i][j][k] < i + 1) {
            cout << " " << optimal_state[i][j][k] << " " << optimal_rank[i][j][k];
          }
        }
        cout << " |";
      }
      cout << "| ";
    }
    cout << endl;
#   endif

    if (i < seq.length[index] - 1) {
      if (nb_state_sequence < inb_state_sequence) {
        if (nb_state_sequence * nb_state < inb_state_sequence) {
          nb_state_sequence *= nb_state;
        }
        else {
          nb_state_sequence = inb_state_sequence;
        }
      }

      for (j = 0;j < nb_state;j++) {
        for (k = 0;k < nb_state;k++) {
          rank[k] = 0;
        }

        for (k = 0;k < nb_state_sequence;k++) {
          state_in[i][j][k] = D_INF;
          for (m = 0;m < nb_state;m++) {
            buff = cumul_transition[m][j] + forward1[m][rank[m]];
            if (buff > state_in[i][j][k]) {
              state_in[i][j][k] = buff;
              input_state[i][j][k] = m;
              input_rank[i][j][k] = rank[m];
            }
          }

          if (state_in[i][j][k] != D_INF) {
            rank[input_state[i][j][k]]++;
          }
        }

        for (k = nb_state_sequence;k < inb_state_sequence;k++) {
          state_in[i][j][k] = D_INF;
        }
      }

#     ifdef DEBUG
      cout << i << " : ";
      for (j = 0;j < nb_state;j++) {
        cout << j << " :";
        for (k = 0;k < nb_state_sequence;k++) {
          cout << " " << state_in[i][j][k];
          if (state_in[i][j][k] != D_INF) {
            cout << " " << input_state[i][j][k] << " " << input_rank[i][j][k];
          }
          cout << " |";
        }
        cout << "| ";
      }
      cout << endl;
#     endif

    }

    for (j = 0;j < nb_output_process;j++) {
      poutput[j]++;
    }
  }

  // extraction de la vraisemblance du chemin optimal

  for (i = 0;i < nb_state;i++) {
    rank[i] = 0;
  }
  likelihood_cumul = 0.;

  for (i = 0;i < nb_state_sequence;i++) {
    pstate = seq.sequence[index][0] + seq.length[index] - 1;
    forward_max = D_INF;

    for (j = 0;j < nb_state;j++) {
      if (forward1[j][rank[j]] > forward_max) {
        forward_max = forward1[j][rank[j]];
        *pstate = j;
      }
    }

    if (i == 0) {
      state_seq_likelihood = forward_max;
    }

    if (forward_max == D_INF) {
      break;
    }

    // restauration

    brank = rank[*pstate];
    rank[*pstate]++;
    j = seq.length[index] - 1;

#   ifdef DEBUG
    cout << "\n" << *pstate << " " << optimal_occupancy[j][*pstate][brank] << " " << brank << " | ";
#   endif

    do {
      for (k = 0;k < optimal_occupancy[j][*pstate][brank];k++) {
        active_cell[j - k][*pstate] = true;
      }

      for (k = 0;k < optimal_occupancy[j][*pstate][brank] - 1;k++) {
        pstate--;
        *pstate = *(pstate + 1);
      }

      if (j >= optimal_occupancy[j][*pstate][brank]) {
        pstate--;
        *pstate = optimal_state[j][*(pstate + 1)][brank];
        previous_rank = optimal_rank[j][*(pstate + 1)][brank];
        j -= optimal_occupancy[j][*(pstate + 1)][brank];
        brank = previous_rank;

#       ifdef DEBUG
        cout << *pstate << " " << optimal_occupancy[j][*pstate][brank] << " " << brank << " | ";
#       endif

      }
      else {
        j -= optimal_occupancy[j][*pstate][brank];
      }
    }
    while (j >= 0);

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

  for (i = 0;i < seq.length[index];i++) {
    delete [] observation[i];
  }
  delete [] observation;

  delete [] obs_product;

  for (i = 0;i < nb_state;i++) {
    delete [] forward1[i];
  }
  delete [] forward1;

  for (i = 0;i < seq.length[index] - 1;i++) {
    for (j = 0;j < nb_state;j++) {
      delete [] state_in[i][j];
    }
    delete [] state_in[i];
  }
  delete [] state_in;

  delete [] rank;

  for (i = 0;i < seq.length[index] - 1;i++) {
    for (j = 0;j < nb_state;j++) {
      delete [] input_state[i][j];
    }
    delete [] input_state[i];
  }
  delete [] input_state;

  for (i = 0;i < seq.length[index];i++) {
    for (j = 0;j < nb_state;j++) {
      delete [] optimal_state[i][j];
    }
    delete [] optimal_state[i];
  }
  delete [] optimal_state;

  for (i = 0;i < seq.length[index];i++) {
    for (j = 0;j < nb_state;j++) {
      delete [] optimal_occupancy[i][j];
    }
    delete [] optimal_occupancy[i];
  }
  delete [] optimal_occupancy;

  for (i = 0;i < seq.length[index] - 1;i++) {
    for (j = 0;j < nb_state;j++) {
      delete [] input_rank[i][j];
    }
    delete [] input_rank[i];
  }
  delete [] input_rank;

  for (i = 0;i < seq.length[index];i++) {
    for (j = 0;j < nb_state;j++) {
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
 *  arguments : reference sur un objet Semi_markov_data, indice de la sequence,
 *              stream, type de sortie, format de fichier ('a' : ASCII,
 *              's' : Spreadsheet, 'g' : Gnuplot), vraisemblance des donnees.
 *
 *--------------------------------------------------------------*/

double Hidden_semi_markov::viterbi_forward_backward(const Semi_markov_data &seq ,
                                                    int index , ostream &os , int output ,
                                                    char format , double seq_likelihood) const

{
  register int i , j , k , m;
  int *pstate , **poutput;
  double obs_product , buff , state_seq_likelihood , backward_max , **observation ,
         **forward1 , **state_in , **backward , **backward1 , *auxiliary ,
         *occupancy_auxiliary , **backward_output;
  Parametric *occupancy;


  // initialisations

  observation = new double*[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    observation[i] = new double[nb_state];
  }

  forward1 = new double*[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    forward1[i] = new double[nb_state];
  }

  state_in = new double*[seq.length[index] - 1];
  for (i = 0;i < seq.length[index] - 1;i++) {
    state_in[i] = new double[nb_state];
  }

  backward = new double*[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    backward[i] = new double[nb_state];
  }

  backward1 = new double*[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    backward1[i] = new double[nb_state];
  }

  auxiliary = new double[nb_state];
  occupancy_auxiliary = new double[seq.length[index] + 1];

  if (output == SSTATE) {
    backward_output = backward;
  }
  else {
    backward_output = new double*[seq.length[index]];
    for (i = 0;i < seq.length[index];i++) {
      backward_output[i] = new double[nb_state];
    }
  }

  poutput = new int*[nb_output_process];

# ifdef MESSAGE
  int *state_sequence , **input_state , **optimal_state , **optimal_forward_occupancy;

  input_state = new int*[seq.length[index] - 1];
  for (i = 0;i < seq.length[index] - 1;i++) {
    input_state[i] = new int[nb_state];
  }

  optimal_state = new int*[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    optimal_state[i] = new int[nb_state];
  }

  optimal_forward_occupancy = new int*[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    optimal_forward_occupancy[i] = new int[nb_state];
  }

  state_sequence = new int[seq.length[index]];
# endif

  for (i = 0;i < nb_output_process;i++) {
    poutput[i] = seq.sequence[index][i + 1];
  }

  // recurrence "forward"

  for (i = 0;i < seq.length[index];i++) {
    for (j = 0;j < nb_state;j++) {

      // calcul des probabilites d'observation

      observation[i][j] = 0.;
      for (k = 0;k < nb_output_process;k++) {
        if (nonparametric_process[k + 1]) {
          buff = nonparametric_process[k + 1]->observation[j]->cumul[*poutput[k]];
        }
        else {
          buff = parametric_process[k + 1]->observation[j]->cumul[*poutput[k]];
        }

        if (buff == D_INF) {
          observation[i][j] = D_INF;
          break;
        }
        else {
          observation[i][j] += buff;
        }
      }

      switch (state_subtype[j]) {

      // cas etat semi-markovien

      case SEMI_MARKOVIAN : {
        occupancy = nonparametric_process[0]->sojourn_time[j];
        obs_product = 0.;
        forward1[i][j] = D_INF;

        if (i < seq.length[index] - 1) {
          for (k = 1;k <= MIN(i + 1 , occupancy->nb_value - 1);k++) {
            if (observation[i - k + 1][j] == D_INF) {
              break;
            }
            else {
              obs_product += observation[i - k + 1][j];
            }

            if (k < i + 1) {
              buff = obs_product + occupancy->mass[k] + state_in[i - k][j];
            }

            else {
              switch (type) {
              case 'o' :
                buff = obs_product + occupancy->mass[k] + cumul_initial[j];
                break;
              case 'e' :
                buff = obs_product + forward[j]->mass[k] + cumul_initial[j];
                break;
              }
            }

            if (buff > forward1[i][j]) {
              forward1[i][j] = buff;

#             ifdef MESSAGE
              if (k < i + 1) {
                optimal_state[i][j] = input_state[i - k][j];
              }
              optimal_forward_occupancy[i][j] = k;
#             endif

            }
          }
        }

        else {
          for (k = 1;k <= MIN(i + 1 , occupancy->nb_value - 1);k++) {
            if (observation[i - k + 1][j] == D_INF) {
              break;
            }
            else {
              obs_product += observation[i - k + 1][j];
            }

            if (k < i + 1) {
              buff = obs_product + occupancy->cumul[k - 1] + state_in[i - k][j];
            }

            else {
              switch (type) {
              case 'o' :
                buff = obs_product + occupancy->cumul[k - 1] + cumul_initial[j];
                break;
              case 'e' :
                buff = obs_product + forward[j]->cumul[k - 1] + cumul_initial[j];
                break;
              }
            }

            if (buff > forward1[i][j]) {
              forward1[i][j] = buff;

#             ifdef MESSAGE
              if (k < i + 1) {
                optimal_state[i][j] = input_state[i - k][j];
              }
              optimal_forward_occupancy[i][j] = k;
#             endif

            }
          }
        }
        break;
      }

      // cas etat markovien

      case MARKOVIAN : {
        if (i == 0) {
          forward1[i][j] = cumul_initial[j];
        }
        else {
          forward1[i][j] = state_in[i - 1][j];

#         ifdef MESSAGE
          optimal_state[i][j] = input_state[i - 1][j];
#         endif

        }

#       ifdef MESSAGE
        optimal_forward_occupancy[i][j] = 1;
#       endif

        if (forward1[i][j] != D_INF) {
          if (observation[i][j] == D_INF) {
            forward1[i][j] = D_INF;
          }
          else {
            forward1[i][j] += observation[i][j];
          }
        }
        break;
      }
      }
    }

    if (i < seq.length[index] - 1) {
      for (j = 0;j < nb_state;j++) {
        state_in[i][j] = D_INF;
        for (k = 0;k < nb_state;k++) {
          buff = cumul_transition[k][j] + forward1[i][k];
          if (buff > state_in[i][j]) {
            state_in[i][j] = buff;

#           ifdef MESSAGE
            input_state[i][j] = k;
#           endif

          }
        }
      }
    }

    for (j = 0;j < nb_output_process;j++) {
      poutput[j]++;
    }
  }

  // extraction de la vraisemblance du chemin optimal

# ifdef MESSAGE
  pstate = state_sequence + seq.length[index] - 1;
# endif

  state_seq_likelihood = D_INF;
  i = seq.length[index] - 1;
  for (j = 0;j < nb_state;j++) {
    if (forward1[i][j] > state_seq_likelihood) {
      state_seq_likelihood = forward1[i][j];

#     ifdef MESSAGE
      *pstate = j;
#     endif

    }
  }

  if (state_seq_likelihood != D_INF) {

#   ifdef MESSAGE
    i = seq.length[index] - 1;

    do {
      for (j = 0;j < optimal_forward_occupancy[i][*pstate] - 1;j++) {
        pstate--;
        *pstate = *(pstate + 1);
      }

      if (i >= optimal_forward_occupancy[i][*pstate]) {
        pstate--;
        *pstate = optimal_state[i][*(pstate + 1)];
        i -= optimal_forward_occupancy[i][*(pstate + 1)];
      }
      else {
        i -= optimal_forward_occupancy[i][*pstate];
      }
    }
    while (i >= 0);
#   endif

    // recurrence "backward"

    i = seq.length[index] - 1;
    for (j = 0;j < nb_state;j++) {
      backward1[i][j] = 0.;
      backward[i][j] = forward1[i][j];

      if (output == OUT_STATE) {
        backward_output[i][j] = backward[i][j];
      }
    }

    for (i = seq.length[index] - 2;i >= 0;i--) {
      for (j = 0;j < nb_state;j++) {
        switch (state_subtype[j]) {

        // cas etat semi-markovien

        case SEMI_MARKOVIAN : {
          occupancy = nonparametric_process[0]->sojourn_time[j];
          obs_product = 0.;

          for (k = 1;k < MIN(seq.length[index] - i , occupancy->nb_value);k++) {
            if (observation[i + k][j] == D_INF) {
              break;
            }
            else {
              obs_product += observation[i + k][j];
            }

            if (k < seq.length[index] - i - 1) {
              occupancy_auxiliary[k] = backward1[i + k][j] + obs_product + occupancy->mass[k];
            }
            else {
              occupancy_auxiliary[k] = obs_product + occupancy->cumul[k - 1];
            }
          }

          auxiliary[j] = D_INF;
          for (m = k - 1;m >= 1;m--) {
            if (occupancy_auxiliary[m] > auxiliary[j]) {
              auxiliary[j] = occupancy_auxiliary[m];
            }

            // transformation des vraisemblances "semi-markoviennes" en vraisemblances "markoviennes"

            if ((auxiliary[j] != D_INF) && (state_in[i][j] != D_INF)) {
              buff = auxiliary[j] + state_in[i][j];
              if (buff > backward[i + m][j]) {
                backward[i + m][j] = buff;
              }
            }
          }
          break;
        }

        // cas etat markovien

        case MARKOVIAN : {
          if ((backward1[i + 1][j] != D_INF) && (observation[i + 1][j] != D_INF)) {
            auxiliary[j] = backward1[i + 1][j] + observation[i + 1][j];
          }
          else {
            auxiliary[j] = D_INF;
          }
          break;
        }
        }
      }

      for (j = 0;j < nb_state;j++) {
        backward1[i][j] = D_INF;
        for (k = 0;k < nb_state;k++) {
          buff = auxiliary[k] + cumul_transition[j][k];
          if (buff > backward1[i][j]) {
            backward1[i][j] = buff;
          }
        }

        if ((backward1[i][j] != D_INF) && (forward1[i][j] != D_INF)) {
          backward[i][j] = backward1[i][j] + forward1[i][j];
        }
        else {
          backward[i][j] = D_INF;
        }
      }

      switch (output) {

      case IN_STATE : {
        for (j = 0;j < nb_state;j++) {
          switch (state_subtype[j]) {

          // cas etat semi-markovien

          case SEMI_MARKOVIAN : {
            if ((auxiliary[j] != D_INF) && (state_in[i][j] != D_INF)) {
              backward_output[i + 1][j] = auxiliary[j] + state_in[i][j];
            }
            else {
              backward_output[i + 1][j] = D_INF;
            }
            break;
          }

          // cas etat markovien

          case MARKOVIAN : {
            backward_output[i + 1][j] = D_INF;

            if (auxiliary[j] != D_INF) {
              for (k = 0;k < nb_state;k++) {
                if (k != j) {
                  buff = cumul_transition[k][j] + forward1[i][k];
                  if (buff > backward_output[i + 1][j]) {
                    backward_output[i + 1][j] = buff;
                  }
                }
              }

              if (backward_output[i + 1][j] != D_INF) {
                backward_output[i + 1][j] += auxiliary[j];
              }
            }
            break;
          }
          }
        }
        break;
      }

      case OUT_STATE : {
        for (j = 0;j < nb_state;j++) {
          switch (state_subtype[j]) {

          // cas etat semi-markovien

          case SEMI_MARKOVIAN : {
            backward_output[i][j] = backward[i][j];
            break;
          }

          // cas etat markovien

          case MARKOVIAN : {
            backward_output[i][j] = D_INF;

            if (forward1[i][j] != D_INF) {
              for (k = 0;k < nb_state;k++) {
                if (k != j) {
                  buff = auxiliary[k] + cumul_transition[j][k];
                  if (buff >  backward_output[i][j]) {
                    backward_output[i][j] = buff;
                  }
                }
              }

              if (backward_output[i][j] != D_INF) {
                backward_output[i][j] += forward1[i][j];
              }
            }
            break;
          }
          }
        }
        break;
      }
      }
    }

    // cas particulier de rester dans l'etat initial

    for (i = 0;i < nb_state;i++) {
      if ((state_subtype[i] == SEMI_MARKOVIAN) && (cumul_initial[i] != D_INF)) {
        occupancy = nonparametric_process[0]->sojourn_time[i];
        obs_product = 0.;

        for (j = 1;j < MIN(seq.length[index] + 1 , occupancy->nb_value);j++) {
          if (observation[j - 1][i] == D_INF) {
            break;
          }
          else {
            obs_product += observation[j - 1][i];
          }

          if (j < seq.length[index]) {
            switch (type) {
            case 'o' :
              occupancy_auxiliary[j] = backward1[j - 1][i] + obs_product + occupancy->mass[j];
              break;
            case 'e' :
              occupancy_auxiliary[j] = backward1[j - 1][i] + obs_product + forward[i]->mass[j];
              break;
            }
          }

          else {
            switch (type) {
            case 'o' :
              occupancy_auxiliary[j] = obs_product + occupancy->cumul[j - 1];
              break;
            case 'e' :
              occupancy_auxiliary[j] = obs_product + forward[i]->cumul[j - 1];
              break;
            }
          }
        }

        auxiliary[i] = D_INF;
        for (k = j - 1;k >= 1;k--) {
          if (occupancy_auxiliary[k] > auxiliary[i]) {
            auxiliary[i] = occupancy_auxiliary[k];
          }

          // transformation des vraisemblances "semi-markoviennes" en vraisemblances "markoviennes"

          if (auxiliary[i] != D_INF) {
            buff = auxiliary[i] + cumul_initial[i];
            if (buff > backward[k - 1][i]) {
              backward[k - 1][i] = buff;
            }
          }
        }
      }

      if (output == IN_STATE) {
        backward_output[0][i] = backward[0][i];
      }
    }

    // restauration

    pstate = seq.sequence[index][0];

    for (i = 0;i < seq.length[index];i++) {
      backward_max = D_INF;
      for (j = 0;j < nb_state;j++) {
        if (backward[i][j] > backward_max) {
          backward_max = backward[i][j];
          *pstate = j;
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
        if (backward_output[i][j] != D_INF) {
          backward_output[i][j] = exp(backward_output[i][j] - seq_likelihood);
//          backward_output[i][j] = exp(backward_output[i][j] - state_seq_likelihood);
        }
        else {
          backward_output[i][j] = 0.;
        }
      }
    }

    switch (format) {

    case 'a' : {
      switch (output) {
      case SSTATE :
        os << "\n" << SEQ_label[SEQL_MAX_POSTERIOR_STATE_PROBABILITY] << "\n\n";
        break;
      case IN_STATE :
        os << "\n" << SEQ_label[SEQL_MAX_POSTERIOR_IN_STATE_PROBABILITY] << "\n\n";
        break;
      case OUT_STATE :
        os << "\n" << SEQ_label[SEQL_MAX_POSTERIOR_OUT_STATE_PROBABILITY] << "\n\n";
        break;
      }

      seq.profile_ascii_print(os , index , nb_state , backward_output ,
                              STAT_label[STATL_STATE]);

      os << "\n" << SEQ_label[SEQL_STATE_SEQUENCE_LIKELIHOOD] << ": " << state_seq_likelihood
         << "   (" << exp(state_seq_likelihood - seq_likelihood) << ")" << endl;
      break;
    }

    case 's' : {
      switch (output) {
      case SSTATE :
        os << "\n" << SEQ_label[SEQL_MAX_POSTERIOR_STATE_PROBABILITY] << "\n\n";
        break;
      case IN_STATE :
        os << "\n" << SEQ_label[SEQL_MAX_POSTERIOR_IN_STATE_PROBABILITY] << "\n\n";
        break;
      case OUT_STATE :
        os << "\n" << SEQ_label[SEQL_MAX_POSTERIOR_OUT_STATE_PROBABILITY] << "\n\n";
        break;
      }

      seq.profile_spreadsheet_print(os , index , nb_state , backward_output ,
                                    STAT_label[STATL_STATE]);

      os << "\n" << SEQ_label[SEQL_STATE_SEQUENCE_LIKELIHOOD] << "\t" << state_seq_likelihood
         << "\t" << exp(state_seq_likelihood - seq_likelihood) << endl;
      break;
    }

    case 'g' : {
      seq.profile_plot_print(os , index , nb_state , backward_output);
      break;
    }
    }

#   ifdef DEBUG
    if (format != 'g') {
      double ambiguity = 0.;

      pstate = seq.sequence[index][0];
//      if (output == SSTATE) {
      for (i = 0;i < seq.length[index];i++) {
        for (j = 0;j < nb_state;j++) {
          if (j != *pstate) {
            ambiguity += backward_output[i][j];
          }
        }
        pstate++;
      }
      ambiguity *= exp(seq_likelihood - state_seq_likelihood);
/*      }

      else {
        for (i = 0;i < seq.length[index];i++) {
          for (j = 0;j < nb_state;j++) {
            if ((backward[i][j] != D_INF) && (j != *pstate)) {
              ambiguity += exp(backward[i][j] - state_seq_likelihood);
            }
          }
          pstate++;
        }
      } */

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
    delete [] observation[i];
  }
  delete [] observation;

  for (i = 0;i < seq.length[index];i++) {
    delete [] forward1[i];
  }
  delete [] forward1;

  for (i = 0;i < seq.length[index] - 1;i++) {
    delete [] state_in[i];
  }
  delete [] state_in;

  for (i = 0;i < seq.length[index];i++) {
    delete [] backward[i];
  }
  delete [] backward;

  for (i = 0;i < seq.length[index];i++) {
    delete [] backward1[i];
  }
  delete [] backward1;

  delete [] auxiliary;
  delete [] occupancy_auxiliary;

  if (output != SSTATE) {
    for (i = 0;i < seq.length[index];i++) {
      delete [] backward_output[i];
    }
    delete [] backward_output;
  }

  delete [] poutput;

# ifdef MESSAGE
  for (i = 0;i < seq.length[index] - 1;i++) {
    delete [] input_state[i];
  }
  delete [] input_state;

  for (i = 0;i < seq.length[index];i++) {
    delete [] optimal_state[i];
  }
  delete [] optimal_state;

  for (i = 0;i < seq.length[index];i++) {
    delete [] optimal_forward_occupancy[i];
  }
  delete [] optimal_forward_occupancy;

  delete [] state_sequence;
# endif

  return state_seq_likelihood;
}


/*--------------------------------------------------------------*
 *
 *  Estimation des parametres d'une semi-chaine de Markov ordinaire cachee
 *  a partir d'un echantillon de sequences par l'algorithme de Viterbi.
 *
 *  arguments : reference sur un objet Format_error, stream, semi-chaine de Markov cachee initiale,
 *              type d'estimateur pour la reestimation des lois d'occupation des etats,
 *              flag sur le calcul des lois de comptage, nombre d'iterations.
 *
 *--------------------------------------------------------------*/

Hidden_semi_markov* Markovian_sequences::hidden_semi_markov_viterbi_estimation(Format_error &error , ostream &os ,
                                                                               const Hidden_semi_markov &ihsmarkov ,
                                                                               int estimator , bool counting_flag ,
                                                                               int nb_iter) const

{
  bool status;
  register int i , j , k;
  int iter , nb_likelihood_decrease , max_nb_value , *occupancy_survivor ,
      *censored_occupancy_survivor;
  double likelihood = D_INF , previous_likelihood , occupancy_likelihood ,
         observation_likelihood , min_likelihood;
  Parametric *occupancy;
  Reestimation<double> *occupancy_reestim;
  Hidden_semi_markov *hsmarkov;
  Semi_markov_data *seq;


  hsmarkov = 0;
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

  else {
    if (ihsmarkov.type != 'o') {
      status = false;
      error.update(SEQ_error[SEQR_MODEL_TYPE]);
    }

    if (ihsmarkov.nb_output_process != nb_variable) {
      status = false;
      error.update(SEQ_error[SEQR_NB_OUTPUT_PROCESS]);
    }

    else {
      for (i = 0;i < nb_variable;i++) {
        if (((ihsmarkov.nonparametric_process[i + 1]) &&
             (ihsmarkov.nonparametric_process[i + 1]->nb_value != marginal[i]->nb_value)) ||
            ((ihsmarkov.parametric_process[i + 1]) &&
             (ihsmarkov.parametric_process[i + 1]->nb_value < marginal[i]->nb_value))) {
          status = false;
          ostringstream error_message;
          error_message << STAT_label[STATL_OUTPUT_PROCESS] << " " << i + 1 << ": "
                        << SEQ_error[SEQR_NB_OUTPUT];
          error.update((error_message.str()).c_str());
        }

        else if ((ihsmarkov.nonparametric_process[i + 1]) && (!characteristics[i])) {
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
  }

  if ((nb_iter != I_DEFAULT) && (nb_iter < 1)) {
    status = false;
    error.update(STAT_error[STATR_NB_ITERATION]);
  }

  if (status) {
    if (max_length > COUNTING_MAX_LENGTH) {
      counting_flag = false;
    }

    // creation de la semi-chaine de Markov cachee

    hsmarkov = new Hidden_semi_markov(ihsmarkov , false , (int)(max_length * SAMPLE_NB_VALUE_COEFF));

    hsmarkov->create_cumul();
    hsmarkov->log_computation();

    hsmarkov->semi_markov_data = new Semi_markov_data(*this , 0 , false);
    seq = hsmarkov->semi_markov_data;
    seq->type[0] = STATE;
    seq->max_value[0] = hsmarkov->nb_state - 1;
    seq->marginal[0] = new Histogram(hsmarkov->nb_state);

    seq->chain_data = new Chain_data(hsmarkov->type , hsmarkov->nb_state , hsmarkov->nb_state);
    seq->characteristics[0] = new Sequence_characteristics(hsmarkov->nb_state);
    seq->characteristics[0]->create_sojourn_time_histogram(max_length , false);
    seq->create_observation_histogram(hsmarkov->nb_state);

    for (i = 1;i <= hsmarkov->nb_output_process;i++) {
      if ((hsmarkov->parametric_process[i]) && (seq->characteristics[i])) {
        delete seq->characteristics[i];
        seq->characteristics[i] = 0;
      }
    }

    if (estimator == COMPLETE_LIKELIHOOD) {
      max_nb_value = 0;
      for (i = 0;i < hsmarkov->nb_state;i++) {
        if ((hsmarkov->state_subtype[i] == SEMI_MARKOVIAN) &&
            (hsmarkov->nonparametric_process[0]->sojourn_time[i]->alloc_nb_value > max_nb_value)) {
          max_nb_value = hsmarkov->nonparametric_process[0]->sojourn_time[i]->alloc_nb_value;
        }
      }

      occupancy_survivor = new int[MIN(max_nb_value , max_length)];
      censored_occupancy_survivor = new int[MIN(max_nb_value , max_length + 1)];
      occupancy_reestim = new Reestimation<double>(MIN(max_nb_value , max_length + 1));
    }

    iter = 0;
    do {
      iter++;
      previous_likelihood = likelihood;

      likelihood = hsmarkov->viterbi(*seq);

      if (likelihood != D_INF) {
        if (likelihood < previous_likelihood) {
          nb_likelihood_decrease++;
        }
        else {
          nb_likelihood_decrease = 0;
        }

        seq->marginal_histogram_computation(0);

        seq->transition_count_computation(*(seq->chain_data) , hsmarkov);
        seq->sojourn_time_histogram_computation(0);
        seq->observation_histogram_computation();

#       ifdef DEBUG
        for (i = 0;i < hsmarkov->nb_state;i++) {
          if (hsmarkov->state_subtype[i] == SEMI_MARKOVIAN) {
            hsmarkov->nonparametric_process[0]->sojourn_time[i]->computation(seq->sojourn_time[0][i]->nb_value ,
                                                                             OCCUPANCY_THRESHOLD);
          }
        }

        cout << likelihood << " | " << hsmarkov->Semi_markov::likelihood_computation(*seq) << endl;
#       endif

        // reestimation des probabilites initiales

        reestimation(hsmarkov->nb_state , seq->chain_data->initial ,
                     hsmarkov->initial , MIN_PROBABILITY , false);

        // reestimation des probabilites de transition

        for (i = 0;i < hsmarkov->nb_state;i++) {
          reestimation(hsmarkov->nb_state , seq->chain_data->transition[i] ,
                       hsmarkov->transition[i] , MIN_PROBABILITY , false);
        }

        // reestimation des lois d'occupation des etats

        min_likelihood = 0.;

        for (i = 0;i < hsmarkov->nb_state;i++) {
          if (hsmarkov->state_subtype[i] == SEMI_MARKOVIAN) {

#           ifdef DEBUG
            cout << STAT_label[STATL_STATE] << " " << i << " - ";
            seq->characteristics[0]->sojourn_time[i]->ascii_characteristic_print(cout);
#           endif

            occupancy = hsmarkov->nonparametric_process[0]->sojourn_time[i];

            if ((estimator == COMPLETE_LIKELIHOOD) &&
                (seq->characteristics[0]->final_run[i]->nb_element > 0)) {

#             ifdef DEBUG
              cout << STAT_label[STATL_STATE] << " " << i << " - ";
              seq->characteristics[0]->final_run[i]->ascii_characteristic_print(cout);
#             endif

              seq->characteristics[0]->sojourn_time[i]->state_occupancy_estimation(seq->characteristics[0]->final_run[i] ,
                                                                                   occupancy_reestim ,
                                                                                   occupancy_survivor ,
                                                                                   censored_occupancy_survivor);

#             ifdef DEBUG
              cout << STAT_label[STATL_STATE] << " " << i << " ";
              occupancy_reestim->ascii_characteristic_print(cout);
#             endif

              if (iter <= EXPLORATION_NB_ITER) {
                occupancy_likelihood = occupancy_reestim->parametric_estimation(occupancy , 1 , true ,
                                                                                OCCUPANCY_THRESHOLD);
              }
              else {
                occupancy_likelihood = occupancy_reestim->type_parametric_estimation(occupancy , 1 , true ,
                                                                                     OCCUPANCY_THRESHOLD);
              }
            }

            else {
              if (iter <= EXPLORATION_NB_ITER) {
                occupancy_likelihood = seq->characteristics[0]->sojourn_time[i]->Reestimation<int>::parametric_estimation(occupancy , 1 , true ,
                                                                                                                          OCCUPANCY_THRESHOLD);
              }
              else {
                occupancy_likelihood = seq->characteristics[0]->sojourn_time[i]->Reestimation<int>::type_parametric_estimation(occupancy , 1 , true ,
                                                                                                                               OCCUPANCY_THRESHOLD);
              }
            }

            if (occupancy_likelihood == D_INF) {
              min_likelihood = D_INF;
            }
            else {
              occupancy->computation(seq->characteristics[0]->sojourn_time[i]->nb_value , OCCUPANCY_THRESHOLD);
            }

#           ifdef DEBUG
            cout << STAT_word[STATW_STATE] << " " << i << " " << STAT_word[STATW_OCCUPANCY_DISTRIBUTION] << endl;
            occupancy->ascii_print(cout);
#           endif

          }
        }

        // reestimation des lois d'observation

        for (i = 1;i <= hsmarkov->nb_output_process;i++) {
          if (hsmarkov->nonparametric_process[i]) {
            for (j = 0;j < hsmarkov->nb_state;j++) {
              reestimation(seq->marginal[i]->nb_value , seq->observation[i][j]->frequency ,
                           hsmarkov->nonparametric_process[i]->observation[j]->mass ,
                           MIN_PROBABILITY , false);
            }
          }

          else {
            for (j = 0;j < hsmarkov->nb_state;j++) {
              observation_likelihood = seq->observation[i][j]->Reestimation<int>::type_parametric_estimation(hsmarkov->parametric_process[i]->observation[j] ,
                                                                                                             0 , true , OBSERVATION_THRESHOLD);

              if (observation_likelihood == D_INF) {
                min_likelihood = D_INF;
              }
              else {
                hsmarkov->parametric_process[i]->observation[j]->computation(seq->marginal[i]->nb_value ,
                                                                             OBSERVATION_THRESHOLD);

                if (hsmarkov->parametric_process[i]->observation[j]->ident == BINOMIAL) {
                  for (k = hsmarkov->parametric_process[i]->observation[j]->nb_value;k < seq->marginal[i]->nb_value;k++) {
                    hsmarkov->parametric_process[i]->observation[j]->mass[k] = 0.;
                  }
                }
              }
            }
          }
        }

        hsmarkov->log_computation();
      }

#     ifdef MESSAGE
      os << STAT_label[STATL_ITERATION] << " " << iter << "   "
         << SEQ_label[SEQL_STATE_SEQUENCES_LIKELIHOOD] << ": " << likelihood << endl;
#     endif

    }
    while ((likelihood != D_INF) && (((nb_iter == I_DEFAULT) && (iter < SEMI_MARKOV_NB_ITER) &&
             (((likelihood - previous_likelihood) / -likelihood > SEMI_MARKOV_LIKELIHOOD_DIFF) ||
              (min_likelihood == D_INF) || (nb_likelihood_decrease == 1))) ||
            ((nb_iter != I_DEFAULT) && (iter < nb_iter) && ((likelihood > previous_likelihood) ||
              (min_likelihood == D_INF) || (nb_likelihood_decrease == 1)))));

    if (estimator == COMPLETE_LIKELIHOOD) {
      delete [] occupancy_survivor;
      delete [] censored_occupancy_survivor;
      delete occupancy_reestim;
    }

    if (likelihood == D_INF) {
      delete hsmarkov;
      hsmarkov = 0;
      error.update(STAT_error[STATR_ESTIMATION_FAILURE]);
    }

    else {

#     ifdef MESSAGE
      os << "\n" << iter << " " << STAT_label[STATL_ITERATIONS] << endl;
#     endif

      // reestimation des probabilites initiales

      reestimation(hsmarkov->nb_state , seq->chain_data->initial ,
                   hsmarkov->initial , MIN_PROBABILITY , true);

      // reestimation des probabilites de transition

      for (i = 0;i < hsmarkov->nb_state;i++) {
        reestimation(hsmarkov->nb_state , seq->chain_data->transition[i] ,
                     hsmarkov->transition[i] , MIN_PROBABILITY , true);
      }

      for (i = 0;i < hsmarkov->nb_state;i++) {
        if ((hsmarkov->state_subtype[i] == SEMI_MARKOVIAN) &&
            (hsmarkov->nonparametric_process[0]->sojourn_time[i]->mean == 1.)) {
          hsmarkov->state_subtype[i] = MARKOVIAN;
          delete hsmarkov->nonparametric_process[0]->sojourn_time[i];
          hsmarkov->nonparametric_process[0]->sojourn_time[i] = 0;
          delete hsmarkov->forward[i];
          hsmarkov->forward[i] = 0;
        }
      }

      // reestimation des lois d'observation non-parametriques

      for (i = 1;i <= hsmarkov->nb_output_process;i++) {
        if (hsmarkov->nonparametric_process[i]) {
          for (j = 0;j < hsmarkov->nb_state;j++) {
            reestimation(seq->marginal[i]->nb_value , seq->observation[i][j]->frequency ,
                         hsmarkov->nonparametric_process[i]->observation[j]->mass ,
                         MIN_PROBABILITY , true);
          }
        }

        else {
          hsmarkov->parametric_process[i]->nb_value_computation();
        }
      }

      hsmarkov->Chain::log_computation();

      for (i = 1;i <= hsmarkov->nb_output_process;i++) {
        if (hsmarkov->nonparametric_process[i]) {
          for (j = 0;j < hsmarkov->nb_state;j++) {
            hsmarkov->nonparametric_process[i]->observation[j]->log_computation();
          }
        }
      }

      seq->likelihood = hsmarkov->viterbi(*seq);

      seq->marginal_histogram_computation(0);
      seq->build_characteristic(0 , false);

      seq->transition_count_computation(*(seq->chain_data) , hsmarkov);
      seq->observation_histogram_computation();

      hsmarkov->remove_cumul();

      // calcul des lois d'occupation des etats

      for (i = 0;i < hsmarkov->nb_state;i++) {
        if (hsmarkov->state_subtype[i] == SEMI_MARKOVIAN) {
          hsmarkov->nonparametric_process[0]->sojourn_time[i]->computation(seq->characteristics[0]->sojourn_time[i]->nb_value ,
                                                                           OCCUPANCY_THRESHOLD);
          if (hsmarkov->state_type[i] == 'r') {
            hsmarkov->forward[i]->copy(*(hsmarkov->nonparametric_process[0]->sojourn_time[i]));
            hsmarkov->forward[i]->computation(*(hsmarkov->nonparametric_process[0]->sojourn_time[i]));
          }
        }
      }

      for (i = 1;i <= hsmarkov->nb_output_process;i++) {
        if (hsmarkov->nonparametric_process[i]) {
          for (j = 0;j < hsmarkov->nb_state;j++) {
            hsmarkov->nonparametric_process[i]->observation[j]->cumul_computation();

            hsmarkov->nonparametric_process[i]->observation[j]->max_computation();
//            hsmarkov->nonparametric_process[i]->observation[j]->mean_computation();
//            hsmarkov->nonparametric_process[i]->observation[j]->variance_computation();
          }
        }

        // calcul des lois d'observation parametriques

        else {
          for (j = 0;j < hsmarkov->nb_state;j++) {
            hsmarkov->parametric_process[i]->observation[j]->computation(seq->observation[i][j]->nb_value ,
                                                                         OBSERVATION_THRESHOLD);
          }
        }
      }

#     ifdef MESSAGE
      cout << "\n" << SEQ_label[SEQL_STATE_SEQUENCES_LIKELIHOOD] << ": " << seq->likelihood
           << " | " << hsmarkov->Semi_markov::likelihood_computation(*seq) << endl;
#     endif

      // calcul de la vraisemblance et des lois caracteristiques du modele

      seq->hidden_likelihood = hsmarkov->likelihood_computation(*this);

#     ifdef DEBUG
      cout << SEQ_label[SEQL_OBSERVED_SEQUENCES_LIKELIHOOD] << ": " << seq->hidden_likelihood << endl;
#     endif

      hsmarkov->component_computation();
      hsmarkov->characteristic_computation(*seq , counting_flag , I_DEFAULT , false);
    }
  }

  return hsmarkov;
}


/*--------------------------------------------------------------*
 *
 *  Estimation des parametres d'une semi-chaine de Markov cachee a partir
 *  d'un echantillon de sequences par l'algorithme de Viterbi.
 *
 *  arguments : reference sur un objet Format_error, stream,
 *              nombre d'etats de la chaine de Markov, type d'estimateur
 *              pour la reestimation des lois d'occupation des etats,
 *              flag sur le calcul des lois de comptage,
 *              temps moyen d'occupation d'un etat, nombre d'iterations.
 *
 *--------------------------------------------------------------*/

Hidden_semi_markov* Markovian_sequences::hidden_semi_markov_viterbi_estimation(Format_error &error , ostream &os ,
                                                                               int nb_state , int estimator ,
                                                                               bool counting_flag , double occupancy_mean ,
                                                                               int nb_iter) const

{
  bool status = true;
  register int i;
  int nb_value[SEQUENCE_NB_VARIABLE];
  double proba;
  Hidden_semi_markov *ihsmarkov , *hsmarkov;


  hsmarkov = 0;
  error.init();

  if ((nb_state < 2) || (nb_state > NB_STATE)) {
    status = false;
    error.update(SEQ_error[SEQR_NB_STATE]);
  }
  if ((occupancy_mean != D_DEFAULT) && (occupancy_mean <= 1.)) {
    status = false;
    error.update(SEQ_error[SEQR_OCCUPANCY]);
  }

  if (status) {
    for (i = 0;i < nb_variable;i++) {
      nb_value[i] = marginal[i]->nb_value;
    }

    ihsmarkov = new Hidden_semi_markov('o' , nb_state , nb_variable , nb_value);

    // initialisation des parametres de la chaine de Markov

    ihsmarkov->init(true , 0.);

    // initialisation des lois d'occupations des etats

    if (occupancy_mean == D_DEFAULT) {
      occupancy_mean = MAX(hlength->mean , OCCUPANCY_MEAN);
    }

    ihsmarkov->nonparametric_process[0]->absorption = new double[nb_state];
    ihsmarkov->nonparametric_process[0]->sojourn_time = new Parametric*[nb_state];
    ihsmarkov->forward = new Forward*[nb_state];

    for (i = 0;i < nb_state;i++) {
      if (ihsmarkov->state_type[i] != 'a') {
        ihsmarkov->nonparametric_process[0]->absorption[i] = 0.;
        proba = 1. / occupancy_mean;
        ihsmarkov->nonparametric_process[0]->sojourn_time[i] = new Parametric(NEGATIVE_BINOMIAL , 1 , I_DEFAULT ,
                                                                              1. , proba , OCCUPANCY_THRESHOLD);

        if (ihsmarkov->state_type[i] == 'r') {
          ihsmarkov->forward[i] = new Forward(*(ihsmarkov->nonparametric_process[0]->sojourn_time[i]));
        }
        else {
          ihsmarkov->forward[i] = 0;
        }
      }

      else {
        ihsmarkov->nonparametric_process[0]->absorption[i] = 1.;
        ihsmarkov->nonparametric_process[0]->sojourn_time[i] = 0;
        ihsmarkov->forward[i] = 0;
      }
    }

    // initialisation des lois d'observation

    for (i = 0;i < ihsmarkov->nb_output_process;i++) {
      if (ihsmarkov->nonparametric_process[i + 1]) {
        ihsmarkov->nonparametric_process[i + 1]->init();
      }
      else {
        ihsmarkov->parametric_process[i + 1]->init();
      }
    }

    hsmarkov = hidden_semi_markov_viterbi_estimation(error , os , *ihsmarkov , estimator ,
                                                     counting_flag , nb_iter);
    delete ihsmarkov;
  }

  return hsmarkov;
}


/*--------------------------------------------------------------*
 *
 *  Calcul des sequences d'etats optimales.
 *
 *  arguments : references sur un objet Format_error et sur un objet Markovian_sequences,
 *              flag sur le calcul des lois caracteristiques.
 *
 *--------------------------------------------------------------*/

Semi_markov_data* Hidden_semi_markov::state_sequence_computation(Format_error &error ,
                                                                 const Markovian_sequences &iseq ,
                                                                 bool characteristic_flag) const

{
  bool status = true;
  register int i;
  int nb_value;
  Hidden_semi_markov *hsmarkov;
  Semi_markov_data *seq;


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
    seq = new Semi_markov_data(iseq , 0 , (type == 'e' ? true : false));

    seq->type[0] = STATE;
    seq->semi_markov = new Semi_markov(*this , false , false);

    hsmarkov = new Hidden_semi_markov(*this , false , false);

    hsmarkov->create_cumul();
    hsmarkov->log_computation();

    seq->likelihood = hsmarkov->viterbi(*seq);

    delete hsmarkov;

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
      seq->build_characteristic(0 , true , (type == 'e' ? true : false));

      seq->build_transition_count(this);
      seq->build_observation_histogram();

      if (characteristic_flag) {
        seq->semi_markov->characteristic_computation(*seq , true);
      }
    }
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Comparaison de differentes semi-chaines de Markov cachees pour un ensemble
 *  de sequences.
 *
 *  arguments : reference sur un objet Format_error, stream, nombre de semi-chaines
 *              de Markov cachees, pointeur sur les semi-chaines de Markov cachees,
 *              type d'algorithme (forward ou Viterbi), path.
 *
 *--------------------------------------------------------------*/

bool Markovian_sequences::comparison(Format_error &error , ostream &os , int nb_model ,
                                     const Hidden_semi_markov **ihsmarkov , int algorithm ,
                                     const char *path) const

{
  bool status = true;
  register int i , j;
  int nb_value;
  double **likelihood;
  Hidden_semi_markov **hsmarkov;
  Semi_markov_data *seq;


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
    if (ihsmarkov[i]->nb_output_process != nb_variable) {
      status = false;
      ostringstream error_message;
      error_message << SEQ_label[SEQL_HIDDEN_SEMI_MARKOV_CHAIN] << " " << i + 1 << ": "
                    << SEQ_error[SEQR_NB_OUTPUT_PROCESS];
      error.update((error_message.str()).c_str());
    }

    else {
      for (j = 0;j < nb_variable;j++) {
        if (ihsmarkov[i]->nonparametric_process[j + 1]) {
          nb_value = ihsmarkov[i]->nonparametric_process[j + 1]->nb_value;
        }
        else {
          nb_value = ihsmarkov[i]->parametric_process[j + 1]->nb_value;
        }

        if (nb_value < marginal[j]->nb_value) {
          status = false;
          ostringstream error_message;
          error_message << SEQ_label[SEQL_HIDDEN_SEMI_MARKOV_CHAIN] << " " << i + 1 << ": "
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

    hsmarkov = new Hidden_semi_markov*[nb_model];
    for (i = 0;i < nb_model;i++) {
      hsmarkov[i] = new Hidden_semi_markov(*(ihsmarkov[i]) , false , false);
    }

    if (algorithm == VITERBI) {
      for (i = 0;i < nb_model;i++) {
        hsmarkov[i]->create_cumul();
        hsmarkov[i]->log_computation();
      }

      seq = new Semi_markov_data(*this , 0);
    }

    // pour chaque sequence, calcul de la vraisemblance (FORWARD) ou de la vraisemblance
    // du chemin optimal (VITERBI) pour chaque modele possible

    for (i = 0;i < nb_sequence;i++) {
      for (j = 0;j < nb_model;j++) {
        switch (algorithm) {
        case FORWARD :
          likelihood[i][j] = hsmarkov[j]->likelihood_computation(*this , i);
          break;
        case VITERBI :
          likelihood[i][j] = hsmarkov[j]->viterbi(*seq , i);
          break;
        }
      }
    }

#   ifdef MESSAGE
    likelihood_write(os , nb_model , likelihood ,
                     SEQ_label[SEQL_HIDDEN_SEMI_MARKOV_CHAIN] , true , algorithm);
#   endif

    if (path) {
      status = likelihood_write(error , path , nb_model , likelihood ,
                                SEQ_label[SEQL_HIDDEN_SEMI_MARKOV_CHAIN] , algorithm);
    }

    for (i = 0;i < nb_sequence;i++) {
      delete [] likelihood[i];
    }
    delete [] likelihood;

    for (i = 0;i < nb_model;i++) {
      delete hsmarkov[i];
    }
    delete [] hsmarkov;

    if (algorithm == VITERBI) {
      delete seq;
    }
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Simulation par une semi-chaine de Markov cachee.
 *
 *  arguments : reference sur un objet Format_error,
 *              histogramme des longueurs des sequences,
 *              flag sur le calcul des lois de comptage,
 *              flag calcul d'une divergence de Kullback-Leibler.
 *
 *--------------------------------------------------------------*/

Semi_markov_data* Hidden_semi_markov::simulation(Format_error &error , const Histogram &hlength ,
                                                 bool counting_flag , bool divergence_flag) const

{
  Markovian_sequences *observ_seq;
  Semi_markov_data *seq;


  seq = Semi_markov::simulation(error , hlength , counting_flag , divergence_flag);

  if ((seq) && (!divergence_flag)) {
    observ_seq = seq->remove_variable_1();
    seq->hidden_likelihood = likelihood_computation(*observ_seq);
    delete observ_seq;
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Simulation par une semi-chaine de Markov cachee.
 *
 *  arguments : reference sur un objet Format_error,
 *              nombre et longueur des sequences,
 *              flag sur le calcul des lois de comptage.
 *
 *--------------------------------------------------------------*/

Semi_markov_data* Hidden_semi_markov::simulation(Format_error &error , int nb_sequence ,
                                                 int length , bool counting_flag) const

{
  Markovian_sequences *observ_seq;
  Semi_markov_data *seq;


  seq = Semi_markov::simulation(error , nb_sequence , length , counting_flag);

  if (seq) {
    observ_seq = seq->remove_variable_1();
    seq->hidden_likelihood = likelihood_computation(*observ_seq);
    delete observ_seq;
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Simulation par une semi-chaine de Markov cachee.
 *
 *  arguments : reference sur un objet Format_error, nombre de sequences,
 *              reference sur un objet Markovian_sequences,
 *              flag sur le calcul des lois de comptage.
 *
 *--------------------------------------------------------------*/

Semi_markov_data* Hidden_semi_markov::simulation(Format_error &error , int nb_sequence ,
                                                 const Markovian_sequences &iseq , bool counting_flag) const

{
  Markovian_sequences *observ_seq;
  Semi_markov_data *seq;


  seq = Semi_markov::simulation(error , nb_sequence , iseq , counting_flag);

  if (seq) {
    observ_seq = seq->remove_variable_1();
    seq->hidden_likelihood = likelihood_computation(*observ_seq);
    delete observ_seq;
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Comparaison de semi-chaines de Markov cachees par calcul de divergences
 *  de Kullback-Leibler.
 *
 *  arguments : reference sur un objet Format_error, stream, nombre de semi-chaines
 *              de Markov cachees, pointeur sur les semi-chaines de Markov cachees,
 *              histogramme des longueurs des sequences, path.
 *
 *--------------------------------------------------------------*/

Distance_matrix* Hidden_semi_markov::divergence_computation(Format_error &error , ostream &os ,
                                                            int nb_model , const Hidden_semi_markov **ihsmarkov ,
                                                            Histogram **hlength , const char *path) const

{
  bool status = true , lstatus;
  register int i , j , k;
  int cumul_length;
  double ref_likelihood , target_likelihood , **likelihood;
  const Hidden_semi_markov **hsmarkov;
  Markovian_sequences *seq;
  Semi_markov_data *simul_seq;
  Distance_matrix *dist_matrix;
  ofstream *out_file;


  dist_matrix = 0;
  error.init();

  for (i = 0;i < nb_model - 1;i++) {
    if (ihsmarkov[i]->type != type) {
      status = false;
      ostringstream error_message;
      error_message << SEQ_label[SEQL_HIDDEN_SEMI_MARKOV_CHAIN] << " " << i + 2 << ": "
                    << SEQ_error[SEQR_MODEL_TYPE];
      error.update((error_message.str()).c_str());
    }

    if (ihsmarkov[i]->nb_output_process != nb_output_process) {
      status = false;
      ostringstream error_message;
      error_message << SEQ_label[SEQL_HIDDEN_SEMI_MARKOV_CHAIN] << " " << i + 2 << ": "
                    << SEQ_error[SEQR_NB_OUTPUT_PROCESS];
      error.update((error_message.str()).c_str());
    }

    else {
      for (j = 1;j <= nb_output_process;j++) {
        if ((ihsmarkov[i]->nonparametric_process[j]) && (nonparametric_process[j]) &&
            (ihsmarkov[i]->nonparametric_process[j]->nb_value != nonparametric_process[j]->nb_value)) {
          status = false;
          ostringstream error_message;
          error_message << SEQ_label[SEQL_HIDDEN_SEMI_MARKOV_CHAIN] << " " << i + 2 << ": "
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

    hsmarkov = new const Hidden_semi_markov*[nb_model];

    hsmarkov[0] = this;
    for (i = 1;i < nb_model;i++) {
      hsmarkov[i] = ihsmarkov[i - 1];
    }

    dist_matrix = new Distance_matrix(nb_model , SEQ_label[SEQL_HIDDEN_SEMI_MARKOV_CHAIN]);

    for (i = 0;i < nb_model;i++) {

      // simulation d'un echantillon de sequences a partir d'une semi-chaine de Markov cachee

      simul_seq = hsmarkov[i]->simulation(error , *hlength[i] , false , true);
      seq = simul_seq->remove_variable_1();

      likelihood = new double*[seq->nb_sequence];
      for (j = 0;j < seq->nb_sequence;j++) {
        likelihood[j] = new double[nb_model];
      }

      ref_likelihood = 0.;
      for (j = 0;j < seq->nb_sequence;j++) {
        likelihood[j][i] = hsmarkov[i]->likelihood_computation(*seq , j);
        ref_likelihood += likelihood[j][i];
      }

      // calcul des vraisemblances de l'echantillon pour chacune des semi-chaines de Markov cachees

      for (j = 0;j < nb_model;j++) {
        if (j != i) {
          target_likelihood = 0.;
          for (k = 0;k < seq->nb_sequence;k++) {
            likelihood[k][j] = hsmarkov[j]->likelihood_computation(*seq , k);
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
      os << SEQ_label[SEQL_HIDDEN_SEMI_MARKOV_CHAIN] << " " << i + 1 << ": " << seq->nb_sequence << " "
         << SEQ_label[SEQL_SIMULATED] << " " << SEQ_label[seq->nb_sequence == 1 ? SEQL_SEQUENCE : SEQL_SEQUENCES] << endl;
      seq->likelihood_write(os , nb_model , likelihood , SEQ_label[SEQL_HIDDEN_SEMI_MARKOV_CHAIN]);
#     endif

      if (out_file) {
        *out_file << SEQ_label[SEQL_HIDDEN_SEMI_MARKOV_CHAIN] << " " << i + 1 << ": " << seq->nb_sequence << " "
                  << SEQ_label[SEQL_SIMULATED] << " " << SEQ_label[seq->nb_sequence == 1 ? SEQL_SEQUENCE : SEQL_SEQUENCES] << endl;
        seq->likelihood_write(*out_file , nb_model , likelihood , SEQ_label[SEQL_HIDDEN_SEMI_MARKOV_CHAIN]);
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

    delete hsmarkov;
  }

  return dist_matrix;
}


/*--------------------------------------------------------------*
 *
 *  Comparaison de semi-chaines de Markov cachees par calcul de divergences
 *  de Kullback-Leibler.
 *
 *  arguments : reference sur un objet Format_error, stream, nombre de semi-chaines
 *              de Markov cachees, pointeur sur les semi-chaines de Markov cachees,
 *              nombre et longueur des sequences, path.
 *
 *--------------------------------------------------------------*/

Distance_matrix* Hidden_semi_markov::divergence_computation(Format_error &error , ostream &os ,
                                                            int nb_model , const Hidden_semi_markov **hsmarkov ,
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

    dist_matrix = divergence_computation(error , os , nb_model , hsmarkov , hlength , path);

    for (i = 0;i < nb_model;i++) {
      delete hlength[i];
    }
    delete [] hlength;
  }

  return dist_matrix;
}


/*--------------------------------------------------------------*
 *
 *  Comparaison de semi-chaines de Markov cachees par calcul de divergences
 *  de Kullback-Leibler.
 *
 *  arguments : reference sur un objet Format_error, stream, nombre de semi-chaines
 *              de Markov cachees, pointeur sur les semi-chaines de Markov cachees,
 *              pointeurs sur des objets Markovian_sequences, path.
 *
 *--------------------------------------------------------------*/

Distance_matrix* Hidden_semi_markov::divergence_computation(Format_error &error , ostream &os ,
                                                            int nb_model , const Hidden_semi_markov **hsmarkov ,
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

    dist_matrix = divergence_computation(error , os , nb_model , hsmarkov , hlength , path);

    for (i = 0;i < nb_model;i++) {
      delete hlength[i];
    }
    delete [] hlength;
  }

  return dist_matrix;
}
