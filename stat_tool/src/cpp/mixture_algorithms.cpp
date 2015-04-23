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
 *       $Id: mixture_algorithms.cpp 15033 2013-09-30 14:54:49Z guedon $
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

#include "stat_tools.h"
#include "markovian.h"
#include "vectors.h"
#include "mixture.h"
#include "stat_label.h"

#include "stat_tool/distribution_reestimation.hpp"   // probleme compilateur C++ Windows

using namespace std;
using namespace boost::math;


namespace stat_tool {



/*--------------------------------------------------------------*
 *
 *  Calcul du parametre de concentration d'une loi de Von Mises
 *  a partir de la direction moyenne (d'apres Mardia & Jupp, 2000; pp. 85-86).
 *
 *  argument : direction moyenne.
 *
 *--------------------------------------------------------------*/

double von_mises_concentration_computation(double mean_direction)

{
  register int i;
  double concentration , power[6];


  if (mean_direction < 0.53) {
    power[1] = mean_direction;
    for (i = 2;i <= 5;i++) {
      power[i] = power[i - 1] * mean_direction;
    }

    concentration = 2 * power[1] + power[3] + 5 * power[5] / 6;
  }

  else if (mean_direction < 0.85) {
    concentration = -0.4 + 1.39 * mean_direction + 0.43 / (1. - mean_direction);
  }

  else {
    power[1] = 1. - mean_direction;
    for (i = 2;i <= 3;i++) {
      power[i] = power[i - 1] * (1. - mean_direction);
    }

    concentration = 1. / (2 * power[1] - power[2] - power[3]);
  }

# ifdef DEBUG
  cout << "mean direction: " << mean_direction << " | concentration : " << concentration << endl;
# endif

  return concentration;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la quantite d'information d'un objet MixtureData dans
 *  le cas de variables observees discretes.
 *
 *--------------------------------------------------------------*/

double MixtureData::classification_information_computation() const

{
  register int i , j;
  double information , buff;


  information = marginal_distribution[0]->information_computation();

  if (information != D_INF) {
    for (i = 1;i < nb_variable;i++) {
      if (observation_distribution[i]) {
        for (j = 0;j < marginal_distribution[0]->nb_value;j++) {
          if (marginal_distribution[0]->frequency[j] > 0) {
            buff = observation_distribution[i][j]->information_computation();

            if (buff != D_INF) {
              information += buff;
            }
            else {
              information = D_INF;
              break;
            }
          }
        }

        if (information == D_INF) {
          break;
        }
      }

      else {
        information = D_INF;
        break;
      }
    }
  }

  return information;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la vraisemblance de vecteurs affectes aux composantes
 *  pour un melange de lois multivaries.
 *
 *  argument : reference sur un objet MixtureData.
 *
 *--------------------------------------------------------------*/

double Mixture::classification_likelihood_computation(const MixtureData &vec) const

{
  register int i , j;
  int nb_value;
  double buff , likelihood = 0.;


  // verification de la compatibilite entre le modele et les donnees

  if (nb_output_process + 1 == vec.nb_variable) {
    if ((!(vec.marginal_distribution[0])) || (nb_component < vec.marginal_distribution[0]->nb_value)) {
      likelihood = D_INF;
    }

    for (i = 0;i < nb_output_process;i++) {
      if ((categorical_process[i]) || (discrete_parametric_process[i])) {
        if (categorical_process[i]) {
          nb_value = categorical_process[i]->nb_value;
        }
        else {
          nb_value = discrete_parametric_process[i]->nb_value;
        }

        if (nb_value < vec.marginal_distribution[i + 1]->nb_value) {
          likelihood = D_INF;
          break;
        }
      }

      else if (!(vec.marginal_distribution[i + 1])) {
        likelihood = D_INF;
        break;
      }
    }
  }

  else {
    likelihood = D_INF;
  }

  if (likelihood != D_INF) {
    likelihood = weight->likelihood_computation(*(vec.marginal_distribution[0]));

    if (likelihood != D_INF) {
      for (i = 0;i < nb_output_process;i++) {
        if (categorical_process[i]) {
          for (j = 0;j < nb_component;j++) {
            buff = categorical_process[i]->observation[j]->likelihood_computation(*(vec.observation_distribution[i + 1][j]));

            if (buff != D_INF) {
              likelihood += buff;
            }
            else {
              likelihood = D_INF;
              break;
            }
          }
        }

        else if (discrete_parametric_process[i]) {
          for (j = 0;j < nb_component;j++) {
            buff = discrete_parametric_process[i]->observation[j]->likelihood_computation(*(vec.observation_distribution[i + 1][j]));

            if (buff != D_INF) {
              likelihood += buff;
            }
            else {
              likelihood = D_INF;
              break;
            }
          }
        }

        else {
          for (j = 0;j < nb_component;j++) {
            buff = continuous_parametric_process[i]->observation[j]->likelihood_computation(*(vec.observation_distribution[i + 1][j]) ,
                                                                                            (int)vec.min_interval[i]);

            if (buff != D_INF) {
              likelihood += buff;
            }
            else {
              likelihood = D_INF;
              break;
            }
          }
        }

        if (likelihood == D_INF) {
          break;
        }
      }
    }
  }

  return likelihood;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la vraisemblance de vecteurs pour un melange de lois multivaries.
 *
 *  argument : reference sur un objet Vectors, indice du vecteur.
 *
 *--------------------------------------------------------------*/

double Mixture::likelihood_computation(const Vectors &vec , int index) const

{
  register int i , j , k;
  int nb_value;
  double likelihood = 0. , norm , component_proba;


  // verification de la compatibilite entre le modele et les donnees

  if (nb_output_process == vec.nb_variable) {
    for (i = 0;i < nb_output_process;i++) {
      if ((categorical_process[i]) || (discrete_parametric_process[i])) {
        if (categorical_process[i]) {
          nb_value = categorical_process[i]->nb_value;
        }
        else {
          nb_value = discrete_parametric_process[i]->nb_value;
        }

        if (nb_value < vec.marginal_distribution[i]->nb_value) {
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
    likelihood = 0.;

    for (i = 0;i < vec.nb_vector;i++) {
      if ((index == I_DEFAULT) || (index == i)) {
        norm = 0.;

        for (j = 0;j < nb_component;j++) {
          component_proba = 1.;

          for (k = 0;k < nb_output_process;k++) {
            if (categorical_process[k]) {
              component_proba *= categorical_process[k]->observation[j]->mass[vec.int_vector[i][k]];
            }

            else if (discrete_parametric_process[k]) {
              component_proba *= discrete_parametric_process[k]->observation[j]->mass[vec.int_vector[i][k]];
            }

            else {
              if ((continuous_parametric_process[k]->ident == GAMMA) && (vec.min_value[k] < vec.min_interval[k] / 2)) {
                switch (vec.type[k]) {
                case INT_VALUE :
                  component_proba *= continuous_parametric_process[k]->observation[j]->mass_computation(vec.int_vector[i][k] , vec.int_vector[i][k] + vec.min_interval[k]);
                  break;
                case REAL_VALUE :
                  component_proba *= continuous_parametric_process[k]->observation[j]->mass_computation(vec.real_vector[i][k] , vec.real_vector[i][k] + vec.min_interval[k]);
                  break;
                }
              }

              else {
                switch (vec.type[k]) {
                case INT_VALUE :
                  component_proba *= continuous_parametric_process[k]->observation[j]->mass_computation(vec.int_vector[i][k] - vec.min_interval[k] / 2 , vec.int_vector[i][k] + vec.min_interval[k] / 2);
                  break;
                case REAL_VALUE :
                  component_proba *= continuous_parametric_process[k]->observation[j]->mass_computation(vec.real_vector[i][k] - vec.min_interval[k] / 2 , vec.real_vector[i][k] + vec.min_interval[k] / 2);
                  break;
                }
              }
            }
          }

          component_proba *= weight->mass[j];
          norm += component_proba;
        }

        if (norm > 0.) {
          likelihood += log(norm);
        }
        else {
          likelihood = D_INF;
          break;
        }
      }
    }
  }

  return likelihood;
}


/*--------------------------------------------------------------*
 *
 *  Estimation des parametres d'un melange multivarie
 *  a partir d'un echantillon de vecteurs par l'algorithme EM.
 *
 *  arguments : reference sur un objet StatError, stream, melange initial,
 *              flags estimation des composantes,
 *              parametres de dispersion communs (processus d'observation continus),
 *              type de lien entre les variances des lois gaussiennes (CONVOLUTION_FACTOR / SCALING_FACTOR)
 *              flag sur le calcul des affectations optimales, nombre d'iterations.
 *
 *--------------------------------------------------------------*/

Mixture* Vectors::mixture_estimation(StatError &error , ostream &os , const Mixture &imixt ,
                                     bool known_component , bool common_dispersion ,
                                     int variance_factor , bool assignment , int nb_iter) const

{
  bool status;
  register int i , j , k;
  int max_nb_value , iter;
  double likelihood = D_INF , previous_likelihood , observation_likelihood , buff ,
         **component_vector_count , norm , *component_proba , diff , variance ,
         **mean_direction , global_mean_direction , concentration;
  Distribution *weight;
  Reestimation<double> *weight_reestim , ***observation_reestim;
  FrequencyDistribution *hobservation;
  Mixture *mixt;
  MixtureData *vec;


  mixt = NULL;
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

  if (imixt.nb_output_process != nb_variable) {
    status = false;
    error.update(STAT_error[STATR_NB_OUTPUT_PROCESS]);
  }

  else {
    for (i = 0;i < nb_variable;i++) {
      if ((imixt.categorical_process[i]) || (imixt.discrete_parametric_process[i])) {
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
            if (((imixt.categorical_process[i]) &&
                 (imixt.categorical_process[i]->nb_value != marginal_distribution[i]->nb_value)) ||
                ((imixt.discrete_parametric_process[i]) &&
                 (imixt.discrete_parametric_process[i]->nb_value < marginal_distribution[i]->nb_value))) {
              status = false;
              ostringstream error_message;
              error_message << STAT_label[STATL_OUTPUT_PROCESS] << " " << i + 1 << ": "
                            << STAT_error[STATR_NB_OUTPUT];
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

    // creation du melange multivarie

    mixt = new Mixture(imixt , false);

    if (common_dispersion) {
      for (i = 0;i < mixt->nb_output_process;i++) {
        if (mixt->continuous_parametric_process[i]) {
          mixt->continuous_parametric_process[i]->tied_dispersion = true;
        }
      }
    }

#   ifdef DEBUG
    cout << *mixt;
#   endif

    // creation des structures de donnees de l'algorithme
 
    component_proba = new double[mixt->nb_component];

    weight_reestim = new Reestimation<double>(mixt->nb_component);

    if (!known_component) {
      observation_reestim = new Reestimation<double>**[mixt->nb_output_process];
      for (i = 0;i < mixt->nb_output_process;i++) {
        if (marginal_distribution[i]) {
          observation_reestim[i] = new Reestimation<double>*[mixt->nb_component];
          for (j = 0;j < mixt->nb_component;j++) {
            observation_reestim[i][j] = new Reestimation<double>(marginal_distribution[i]->nb_value);
          }
        }

        else {
          observation_reestim[i] = NULL;
        }
      }

      max_nb_value = 0;
      for (i = 0;i < mixt->nb_output_process;i++) {
        if ((mixt->discrete_parametric_process[i]) &&
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

      for (i = 0;i < mixt->nb_output_process;i++) {
        if (!marginal_distribution[i]) {
          break;
        }
      }

      if (i < mixt->nb_output_process) {
        component_vector_count = new double*[nb_vector];
        for (i = 0;i < nb_vector;i++) {
          component_vector_count[i] = new double[mixt->nb_component];
        }
      }
      else {
        component_vector_count = NULL;
      }

      for (i = 0;i < mixt->nb_output_process;i++) {
        if ((mixt->continuous_parametric_process[i]) &&
            (mixt->continuous_parametric_process[i]->ident == VON_MISES)) {
          break;
        }
      }

      if (i < mixt->nb_output_process) {
        mean_direction = new double*[mixt->nb_component];
        for (i = 0;i < mixt->nb_component;i++) {
          mean_direction[i] = new double[4];
        }
      }
      else {
        mean_direction = NULL;
      }
    }

    iter = 0;
    do {
      iter++;
      previous_likelihood = likelihood;
      likelihood = 0.;

      // initialisation des quantites de reestimation

      for (i = 0;i < mixt->nb_component;i++) {
        weight_reestim->frequency[i] = 0.;
      }

      if (!known_component) {
        for (i = 0;i < mixt->nb_output_process;i++) {
          if (observation_reestim[i]) {
            for (j = 0;j < mixt->nb_component;j++) {
              for (k = 0;k < marginal_distribution[i]->nb_value;k++) {
                observation_reestim[i][j]->frequency[k] = 0.;
              }
            }
          }
        }

        if (component_vector_count) {
          for (i = 0;i < nb_vector;i++) {
            for (j = 0;j < mixt->nb_component;j++) {
              component_vector_count[i][j] = 0.;
            }
          }
        }
      }

      for (i = 0;i < nb_vector;i++) {
        norm = 0.;

        // calcul des probabilites des composantes

        for (j = 0;j < mixt->nb_component;j++) {
          component_proba[j] = 1.;

          for (k = 0;k < mixt->nb_output_process;k++) {
            if (mixt->categorical_process[k]) {
              component_proba[j] *= mixt->categorical_process[k]->observation[j]->mass[int_vector[i][k]];
            }

            else if (mixt->discrete_parametric_process[k]) {
              component_proba[j] *= mixt->discrete_parametric_process[k]->observation[j]->mass[int_vector[i][k]];
            }

            else {
              if ((mixt->continuous_parametric_process[k]->ident == GAMMA) && (min_value[k] < min_interval[k] / 2)) {
                switch (type[k]) {
                case INT_VALUE :
                  component_proba[j] *= mixt->continuous_parametric_process[k]->observation[j]->mass_computation(int_vector[i][k] , int_vector[i][k] + min_interval[k]);
                  break;
                case REAL_VALUE :
                  component_proba[j] *= mixt->continuous_parametric_process[k]->observation[j]->mass_computation(real_vector[i][k] , real_vector[i][k] + min_interval[k]);
                  break;
                }
              }

              else {
                switch (type[k]) {
                case INT_VALUE :
                  component_proba[j] *= mixt->continuous_parametric_process[k]->observation[j]->mass_computation(int_vector[i][k] - min_interval[k] / 2 , int_vector[i][k] + min_interval[k] / 2);
                  break;
                case REAL_VALUE :
                  component_proba[j] *= mixt->continuous_parametric_process[k]->observation[j]->mass_computation(real_vector[i][k] - min_interval[k] / 2 , real_vector[i][k] + min_interval[k] / 2);
                  break;
                }
              }
            }
          }

          component_proba[j] *= mixt->weight->mass[j];
          norm += component_proba[j];
        }

        if (norm > 0.) {
          likelihood += log(norm);

          for (j = 0;j < mixt->nb_component;j++) {
            component_proba[j] /= norm;

            // accumulation des quantites de reestimation des poids

            weight_reestim->frequency[j] += component_proba[j];

            // accumulation des quantites de reestimation des lois d'observation

            if (!known_component) {
              for (k = 0;k < mixt->nb_output_process;k++) {
                if (observation_reestim[k]) {
                  observation_reestim[k][j]->frequency[int_vector[i][k]] += component_proba[j];
                }
              }

              if (component_vector_count) {
                component_vector_count[i][j] += component_proba[j];
              }
            }
          }
        }

        else {
          likelihood = D_INF;
          break;
        }
      }

      if (likelihood != D_INF) {

        // reestimation des poids

        reestimation(mixt->nb_component , weight_reestim->frequency ,
                     mixt->weight->mass , MIN_PROBABILITY , false);

        // reestimation des lois d'observation

        if (!known_component) {
          for (i = 0;i < mixt->nb_output_process;i++) {
            if (mixt->categorical_process[i]) {
              for (j = 0;j < mixt->nb_component;j++) {
                reestimation(marginal_distribution[i]->nb_value , observation_reestim[i][j]->frequency ,
                             mixt->categorical_process[i]->observation[j]->mass ,
                             MIN_PROBABILITY , false);
              }
            }

            else if (observation_reestim[i]) {
              for (j = 0;j < mixt->nb_component;j++) {
                observation_reestim[i][j]->nb_value_computation();
                observation_reestim[i][j]->offset_computation();
                observation_reestim[i][j]->nb_element_computation();
                observation_reestim[i][j]->max_computation();
                observation_reestim[i][j]->mean_computation();
                observation_reestim[i][j]->variance_computation(true);
//                observation_reestim[i][j]->variance_computation();
              }

              if (mixt->discrete_parametric_process[i]) {
                for (j = 0;j < mixt->nb_component;j++) {
                  hobservation->update(observation_reestim[i][j] ,
                                       MAX((int)(observation_reestim[i][j]->nb_element *
                                                 MAX(sqrt(observation_reestim[i][j]->variance) , 1.) * OBSERVATION_COEFF) , MIN_NB_ELEMENT));
                  observation_likelihood = hobservation->Reestimation<int>::type_parametric_estimation(mixt->discrete_parametric_process[i]->observation[j] ,
                                                                                                       0 , true , OBSERVATION_THRESHOLD);

                  if (observation_likelihood != D_INF) {
                    mixt->discrete_parametric_process[i]->observation[j]->computation(marginal_distribution[i]->nb_value ,
                                                                                      OBSERVATION_THRESHOLD);

                    if (mixt->discrete_parametric_process[i]->observation[j]->ident == BINOMIAL) {
                      for (k = mixt->discrete_parametric_process[i]->observation[j]->nb_value;k < marginal_distribution[i]->nb_value;k++) {
                        mixt->discrete_parametric_process[i]->observation[j]->mass[k] = 0.;
                      }
                    }
                  }
                }
              }

              else {
                switch (mixt->continuous_parametric_process[i]->ident) {

                case GAMMA : {
                  for (j = 0;j < mixt->nb_component;j++) {
                    observation_reestim[i][j]->gamma_estimation(mixt->continuous_parametric_process[i]->observation[j] , iter);
                  }
                  break;
                }

                case GAUSSIAN : {
                  for (j = 0;j < mixt->nb_component;j++) {
                    mixt->continuous_parametric_process[i]->observation[j]->location = observation_reestim[i][j]->mean;
                  }

                  switch (common_dispersion) {

                  case false : {
                    for (j = 0;j < mixt->nb_component;j++) {
                      mixt->continuous_parametric_process[i]->observation[j]->dispersion = sqrt(observation_reestim[i][j]->variance);
                      if (mixt->continuous_parametric_process[i]->observation[j]->dispersion /
                          mixt->continuous_parametric_process[i]->observation[j]->location < GAUSSIAN_MIN_VARIATION_COEFF) {
                        mixt->continuous_parametric_process[i]->observation[j]->dispersion = mixt->continuous_parametric_process[i]->observation[j]->location * GAUSSIAN_MIN_VARIATION_COEFF;
                      }
                    }
                    break;
                  }

                  case true : {
                    variance = 0.;
                    buff = 0.;

                    for (j = 0;j < mixt->nb_component;j++) {
                      for (k = observation_reestim[i][j]->offset;k < observation_reestim[i][j]->nb_value;k++) {
                        diff = k - observation_reestim[i][j]->mean;
                        variance += observation_reestim[i][j]->frequency[k] * diff * diff;
                      }

                      buff += observation_reestim[i][j]->nb_element;
                    }

                    variance /= buff;
//                    variance /= (buff - 1);

                    for (j = 0;j < mixt->nb_component;j++) {
                      mixt->continuous_parametric_process[i]->observation[j]->dispersion = sqrt(variance);
                    }
                    break;
                  }
                  }

                  break;
                }

                case VON_MISES : {
                  for (j = 0;j < mixt->nb_component;j++) {
                    observation_reestim[i][j]->mean_direction_computation(mean_direction[j]);
                    mixt->continuous_parametric_process[i]->observation[j]->location = mean_direction[j][3];
                  }

                  switch (common_dispersion) {

                  case false : {
                    for (j = 0;j < mixt->nb_component;j++) {
                      mixt->continuous_parametric_process[i]->observation[j]->dispersion = von_mises_concentration_computation(mean_direction[j][2]);
                    }
                    break;
                  }

                  case true : {
                    global_mean_direction = 0.;
                    buff = 0.;

                    for (j = 0;j < mixt->nb_component;j++) {
                      global_mean_direction += observation_reestim[i][j]->nb_element * mean_direction[j][2];
                      buff += observation_reestim[i][j]->nb_element;
                    }
                    concentration = von_mises_concentration_computation(global_mean_direction / buff);

                    for (j = 0;j < mixt->nb_component;j++) {
                      mixt->continuous_parametric_process[i]->observation[j]->dispersion = concentration;
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
              if (variance_factor != INDEPENDENT) {
                switch (mixt->continuous_parametric_process[i]->ident) {
                case GAMMA :
                  tied_gamma_estimation(component_vector_count , i ,
                                        mixt->continuous_parametric_process[i] ,
                                        variance_factor , iter);
                  break;
                case GAUSSIAN :
                  tied_gaussian_estimation(component_vector_count , i ,
                                           mixt->continuous_parametric_process[i] ,
                                           variance_factor);
                  break;
                }
              }

              else {
                switch (mixt->continuous_parametric_process[i]->ident) {
                case GAMMA :
                  gamma_estimation(component_vector_count , i ,
                                   mixt->continuous_parametric_process[i] , iter);
                  break;
                case GAUSSIAN :
                  gaussian_estimation(component_vector_count , i ,
                                      mixt->continuous_parametric_process[i]);
                  break;
                case VON_MISES :
                  von_mises_estimation(component_vector_count , i ,
                                       mixt->continuous_parametric_process[i]);
                  break;
                }
              }
            }
          }
        }
      }

#     ifdef MESSAGE
      os << STAT_label[STATL_ITERATION] << " " << iter << "   "
         << STAT_label[STATL_LIKELIHOOD] << ": " << likelihood << endl;
#     endif

#     ifdef DEBUG
      if (iter % 5 == 0) {
        cout << *mixt;
      }
#     endif

    }
    while ((likelihood != D_INF) && (((nb_iter == I_DEFAULT) && (iter < MIXTURE_NB_ITER) &&
             ((likelihood - previous_likelihood) / -likelihood > MIXTURE_LIKELIHOOD_DIFF)) ||
            ((nb_iter != I_DEFAULT) && (iter < nb_iter))));

    if (likelihood != D_INF) {

#     ifdef MESSAGE
      os << "\n" << iter << " " << STAT_label[STATL_ITERATIONS] << endl;
#     endif

      reestimation(mixt->nb_component , weight_reestim->frequency ,
                   mixt->weight->mass , MIN_PROBABILITY , true);

      // reestimation des lois d'observation categorielles

      if (!known_component) {
        for (i = 0;i < mixt->nb_output_process;i++) {
          if (mixt->categorical_process[i]) {
            for (j = 0;j < mixt->nb_component;j++) {
              reestimation(marginal_distribution[i]->nb_value , observation_reestim[i][j]->frequency ,
                           mixt->categorical_process[i]->observation[j]->mass ,
                           MIN_PROBABILITY , true);
            }
          }

          else if (mixt->discrete_parametric_process[i]) {
            mixt->discrete_parametric_process[i]->nb_value_computation();
          }
        }
      }
    }

    // destruction des structures de donnees de l'algorithme

    delete [] component_proba;

    delete weight_reestim;

    if (!known_component) {
      for (i = 0;i < mixt->nb_output_process;i++) {
        if (observation_reestim[i]) {
          for (j = 0;j < mixt->nb_component;j++) {
            delete observation_reestim[i][j];
          }
          delete [] observation_reestim[i];
        }
      }
      delete [] observation_reestim;

      delete hobservation;

      if (component_vector_count) {
        for (i = 0;i < nb_vector;i++) {
          delete [] component_vector_count[i];
        }
        delete [] component_vector_count;
      }

      if (mean_direction) {
        for (i = 0;i < mixt->nb_component;i++) {
          delete [] mean_direction[i];
        }
        delete [] mean_direction;
      }
    }

    if (likelihood == D_INF) {
      delete mixt;
      mixt = NULL;
      error.update(STAT_error[STATR_ESTIMATION_FAILURE]);
    }

    else {
      if (assignment) {
        mixt->mixture_data = new MixtureData(*this , 'a');
        vec = mixt->mixture_data;

        mixt->individual_assignment(*vec , true);

#       ifdef MESSAGE
        os << "\n" << STAT_label[STATL_CLASSIFICATION_LIKELIHOOD] << ": " << vec->restoration_likelihood;

        for (i = 0;i < nb_variable;i++) {
          if (type[i] == REAL_VALUE) {
            break;
          }
        }
        if (i == nb_variable) {
          os << " | " << mixt->classification_likelihood_computation(*vec)
             << " (" << vec->classification_information_computation() << ")" << endl;
        }
        cout << endl;
#       endif

        // calcul des melanges de lois d'observation (poids deduits de la restauration)

        for (i = 0;i < mixt->nb_output_process;i++) {
          if (mixt->categorical_process[i]) {
            mixt->categorical_process[i]->restoration_weight = new Distribution(*(vec->marginal_distribution[0]));
            mixt->categorical_process[i]->restoration_mixture = mixt->categorical_process[i]->mixture_computation(mixt->categorical_process[i]->restoration_weight);
          }

          else if (mixt->discrete_parametric_process[i]) {
            mixt->discrete_parametric_process[i]->restoration_weight = new Distribution(*(vec->marginal_distribution[0]));
            mixt->discrete_parametric_process[i]->restoration_mixture = mixt->discrete_parametric_process[i]->mixture_computation(mixt->discrete_parametric_process[i]->restoration_weight);
          }

          else if (mixt->continuous_parametric_process[i]) {
            mixt->continuous_parametric_process[i]->restoration_weight = new Distribution(*(vec->marginal_distribution[0]));
          }
        }
      }

      else {
        mixt->mixture_data = new MixtureData(*this);
        vec = mixt->mixture_data;
        if (vec->type[0] == STATE) {
          vec->state_variable_init(INT_VALUE);
        }
      }

      mixt->weight->cumul_computation();
      mixt->weight->max_computation();
      mixt->weight->mean_computation();
      mixt->weight->variance_computation();

      if (!known_component) {
        for (i = 0;i < mixt->nb_output_process;i++) {
          if (mixt->categorical_process[i]) {
            for (j = 0;j < mixt->nb_component;j++) {
              mixt->categorical_process[i]->observation[j]->cumul_computation();

              mixt->categorical_process[i]->observation[j]->max_computation();
//              mixt->categorical_process[i]->observation[j]->mean_computation();
//              mixt->categorical_process[i]->observation[j]->variance_computation();
            }
          }
        }
      }

      // calcul de la vraisemblance du modele

      vec->likelihood = mixt->likelihood_computation(*this);

#     ifdef DEBUG
//      cout << *mixt;
      cout << "iteration " << iter << "  "
           << STAT_label[STATL_LIKELIHOOD] << ": " << vec->likelihood << endl;
#     endif

#     ifdef MESSAGE
      if  ((assignment) && (vec->nb_vector <= POSTERIOR_PROBABILITY_NB_VECTOR)) {
        os << "\n" << STAT_label[STATL_POSTERIOR_ASSIGNMENT_PROBABILITY] << endl;
        for (i = 0;i < vec->nb_vector;i++) {
          os << STAT_label[STATL_VECTOR] << " " << vec->identifier[i] << ": "
             << vec->posterior_probability[i] << endl;
        }
      }
#     endif

      // calcul des melanges de lois d'observation (poids estimes)

      for (i = 0;i < mixt->nb_output_process;i++) {
        if (mixt->categorical_process[i]) {
          mixt->categorical_process[i]->weight = new Distribution(*(mixt->weight));
          mixt->categorical_process[i]->mixture = mixt->categorical_process[i]->mixture_computation(mixt->categorical_process[i]->weight);
        }

        else if (mixt->discrete_parametric_process[i]) {
          mixt->discrete_parametric_process[i]->weight = new Distribution(*(mixt->weight));
          mixt->discrete_parametric_process[i]->mixture = mixt->discrete_parametric_process[i]->mixture_computation(mixt->discrete_parametric_process[i]->weight);
        }

        else if (mixt->continuous_parametric_process[i]) {
          mixt->continuous_parametric_process[i]->weight = new Distribution(*(mixt->weight));
        }
      }
    }
  }

  return mixt;
}


/*--------------------------------------------------------------*
 *
 *  Estimation des parametres d'un melange gamma ou gaussien univarie a parametres lies
 *  a partir d'un echantillon de vecteurs par l'algorithme EM.
 *
 *  arguments : reference sur un objet StatError, stream, nombre de composantes,
 *              identificateur des composantes (GAMMA / GAUSSIAN),
 *              moyenne et parametre de forme de la premiere composante (GAMMA) /
 *              moyenne et ecart-type de la premiere composante (GAUSSIAN), flag moyennes liees,
 *              type de lien entre les variances (CONVOLUTION_FACTOR / SCALING_FACTOR)
 *              flag sur le calcul des affectations optimales, nombre d'iterations.
 *
 *--------------------------------------------------------------*/

Mixture* Vectors::mixture_estimation(StatError &error , ostream &os , int nb_component ,
                                     int ident , double mean , double standard_deviation ,
                                     bool tied_mean , int variance_factor ,
                                     bool assignment , int nb_iter) const

{
  register int i , j;
  Mixture *imixt , *mixt;


  error.init();

  if ((nb_component < 2) || (nb_component > MIXTURE_NB_COMPONENT)) {
    mixt = NULL;
    error.update(STAT_error[STATR_NB_DISTRIBUTION]);
  }

  else {
    imixt = new Mixture(nb_component , ident , mean , standard_deviation , tied_mean , variance_factor);

    mixt = mixture_estimation(error , os , *imixt , false , false ,
                              variance_factor , assignment , nb_iter);
  }

  return mixt;
}


/*--------------------------------------------------------------*
 *
 *  Estimation des parametres d'un melange multivarie
 *  a partir d'un echantillon de vecteurs par l'algorithme MCEM.
 *
 *  arguments : reference sur un objet StatError, stream, melange initial,
 *              flags estimation des composantes,
 *              parametres de dispersion communs (processus d'observation continus),
 *              type de lien entre les variances des lois gaussiennes (CONVOLUTION_FACTOR / SCALING_FACTOR),
 *              parametres pour le nombre d'affectation des individus simulees,
 *              flag sur le calcul des affectations optimales, nombre d'iterations.
 *
 *--------------------------------------------------------------*/

Mixture* Vectors::mixture_stochastic_estimation(StatError &error , ostream &os , const Mixture &imixt ,
                                                bool known_component , bool common_dispersion ,
                                                int variance_factor , int min_nb_assignment ,
                                                int max_nb_assignment , double parameter ,
                                                bool assignment , int nb_iter) const

{
  bool status;
  register int i , j , k , m;
  int max_nb_value , iter , nb_assignment , **component_vector_count;
  double likelihood = D_INF , previous_likelihood , observation_likelihood , buff ,
         norm , *component_proba , *component_cumul , diff , variance ,
         **mean_direction , global_mean_direction , concentration;
  Distribution *weight;
  Reestimation<double> *weight_reestim , ***observation_reestim;
  FrequencyDistribution *hobservation;
  Mixture *mixt;
  MixtureData *vec;


  mixt = NULL;
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

  if (imixt.nb_output_process != nb_variable) {
    status = false;
    error.update(STAT_error[STATR_NB_OUTPUT_PROCESS]);
  }

  else {
    for (i = 0;i < nb_variable;i++) {
      if ((imixt.categorical_process[i]) || (imixt.discrete_parametric_process[i])) {
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
            if (((imixt.categorical_process[i]) &&
                 (imixt.categorical_process[i]->nb_value != marginal_distribution[i]->nb_value)) ||
                ((imixt.discrete_parametric_process[i]) &&
                 (imixt.discrete_parametric_process[i]->nb_value < marginal_distribution[i]->nb_value))) {
              status = false;
              ostringstream error_message;
              error_message << STAT_label[STATL_OUTPUT_PROCESS] << " " << i + 1 << ": "
                            << STAT_error[STATR_NB_OUTPUT];
              error.update((error_message.str()).c_str());
            }
          }
        }
      }
    }
  }

  if ((min_nb_assignment < 1) || (min_nb_assignment > max_nb_assignment)) {
    status = false;
    error.update(STAT_error[STATR_MIN_NB_ASSIGNMENT]);
  }

  if ((nb_iter != I_DEFAULT) && (nb_iter < 1)) {
    status = false;
    error.update(STAT_error[STATR_NB_ITERATION]);
  }

  if (status) {

    // creation du melange multivarie

    mixt = new Mixture(imixt , false);

    if (common_dispersion) {
      for (i = 0;i < mixt->nb_output_process;i++) {
        if (mixt->continuous_parametric_process[i]) {
          mixt->continuous_parametric_process[i]->tied_dispersion = true;
        }
      }
    }

#   ifdef DEBUG
    cout << *mixt;
#   endif

    // creation des structures de donnees de l'algorithme
 
    component_proba = new double[mixt->nb_component];
    component_cumul = new double[mixt->nb_component];

    weight_reestim = new Reestimation<double>(mixt->nb_component);

    if (!known_component) {
      observation_reestim = new Reestimation<double>**[mixt->nb_output_process];
      for (i = 0;i < mixt->nb_output_process;i++) {
        if (marginal_distribution[i]) {
          observation_reestim[i] = new Reestimation<double>*[mixt->nb_component];
          for (j = 0;j < mixt->nb_component;j++) {
            observation_reestim[i][j] = new Reestimation<double>(marginal_distribution[i]->nb_value);
          }
        }

        else {
          observation_reestim[i] = NULL;
        }
      }

      max_nb_value = 0;
      for (i = 0;i < mixt->nb_output_process;i++) {
        if ((mixt->discrete_parametric_process[i]) &&
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

      for (i = 0;i < mixt->nb_output_process;i++) {
        if (!marginal_distribution[i]) {
          break;
        }
      }

      if (i < mixt->nb_output_process) {
        component_vector_count = new int*[nb_vector];
        for (i = 0;i < nb_vector;i++) {
          component_vector_count[i] = new int[mixt->nb_component];
        }
      }
      else {
        component_vector_count = NULL;
      }

      for (i = 0;i < mixt->nb_output_process;i++) {
        if ((mixt->continuous_parametric_process[i]) &&
            (mixt->continuous_parametric_process[i]->ident == VON_MISES)) {
          break;
        }
      }

      if (i < mixt->nb_output_process) {
        mean_direction = new double*[mixt->nb_component];
        for (i = 0;i < mixt->nb_component;i++) {
          mean_direction[i] = new double[4];
        }
      }
      else {
        mean_direction = NULL;
      }
    }

    iter = 0;
    do {
      previous_likelihood = likelihood;
      likelihood = 0.;

      // calcul du nombre d'affectations des individus simulees

      if (min_nb_assignment + (int)::round(parameter * iter) < max_nb_assignment) {
        nb_assignment = min_nb_assignment + (int)::round(parameter * iter);
      }
      else {
        nb_assignment = max_nb_assignment;
      }

/*      nb_assignment = max_nb_assignment - (int)::round((max_nb_assignment - min_nb_assignment) *
                          exp(-parameter * iter)); */

      iter++;

      // initialisation des quantites de reestimation

      for (i = 0;i < mixt->nb_component;i++) {
        weight_reestim->frequency[i] = 0.;
      }

      if (!known_component) {
        for (i = 0;i < mixt->nb_output_process;i++) {
          if (observation_reestim[i]) {
            for (j = 0;j < mixt->nb_component;j++) {
              for (k = 0;k < marginal_distribution[i]->nb_value;k++) {
                observation_reestim[i][j]->frequency[k] = 0.;
              }
            }
          }
        }

        if (component_vector_count) {
          for (i = 0;i < nb_vector;i++) {
            for (j = 0;j < mixt->nb_component;j++) {
              component_vector_count[i][j] = 0;
            }
          }
        }
      }

      for (i = 0;i < nb_vector;i++) {
        norm = 0.;

        // calcul des probabilites des composantes

        for (j = 0;j < mixt->nb_component;j++) {
          component_proba[j] = 1.;

          for (k = 0;k < mixt->nb_output_process;k++) {
            if (mixt->categorical_process[k]) {
              component_proba[j] *= mixt->categorical_process[k]->observation[j]->mass[int_vector[i][k]];
            }

            else if (mixt->discrete_parametric_process[k]) {
              component_proba[j] *= mixt->discrete_parametric_process[k]->observation[j]->mass[int_vector[i][k]];
            }

            else {
              if ((mixt->continuous_parametric_process[k]->ident == GAMMA) && (min_value[k] < min_interval[k] / 2)) {
                switch (type[k]) {
                case INT_VALUE :
                  component_proba[j] *= mixt->continuous_parametric_process[k]->observation[j]->mass_computation(int_vector[i][k] , int_vector[i][k] + min_interval[k]);
                  break;
                case REAL_VALUE :
                  component_proba[j] *= mixt->continuous_parametric_process[k]->observation[j]->mass_computation(real_vector[i][k] , real_vector[i][k] + min_interval[k]);
                  break;
                }
              }

              else {
                switch (type[k]) {
                case INT_VALUE :
                  component_proba[j] *= mixt->continuous_parametric_process[k]->observation[j]->mass_computation(int_vector[i][k] - min_interval[k] / 2 , int_vector[i][k] + min_interval[k] / 2);
                  break;
                case REAL_VALUE :
                  component_proba[j] *= mixt->continuous_parametric_process[k]->observation[j]->mass_computation(real_vector[i][k] - min_interval[k] / 2 , real_vector[i][k] + min_interval[k] / 2);
                  break;
                }
              }
            }
          }

          component_proba[j] *= mixt->weight->mass[j];
          norm += component_proba[j];
        }

        if (norm > 0.) {
          likelihood += log(norm);

          for (j = 0;j < mixt->nb_component;j++) {
            component_proba[j] /= norm;
          }
          cumul_computation(mixt->nb_component , component_proba , component_cumul);

          for (j = 0;j < nb_assignment;j++) {
            k = cumul_method(mixt->nb_component , component_cumul);

            // accumulation des quantites de reestimation des poids

            (weight_reestim->frequency[k])++;

            // accumulation des quantites de reestimation des lois d'observation

            if (!known_component) {
              for (m = 0;m < mixt->nb_output_process;m++) {
                if (observation_reestim[m]) {
                  (observation_reestim[m][k]->frequency[int_vector[i][m]])++;
                }
              }

              if (component_vector_count) {
                (component_vector_count[i][k])++;
              }
            }
          }
        }

        else {
          likelihood = D_INF;
          break;
        }
      }

      if (likelihood != D_INF) {

        // reestimation des poids

        reestimation(mixt->nb_component , weight_reestim->frequency ,
                     mixt->weight->mass , MIN_PROBABILITY , false);

        // reestimation des lois d'observation

        if (!known_component) {
          for (i = 0;i < mixt->nb_output_process;i++) {
            if (mixt->categorical_process[i]) {
              for (j = 0;j < mixt->nb_component;j++) {
                reestimation(marginal_distribution[i]->nb_value , observation_reestim[i][j]->frequency ,
                             mixt->categorical_process[i]->observation[j]->mass ,
                             MIN_PROBABILITY , false);
              }
            }

            else if (observation_reestim[i]) {
              for (j = 0;j < mixt->nb_component;j++) {
                observation_reestim[i][j]->nb_value_computation();
                observation_reestim[i][j]->offset_computation();
                observation_reestim[i][j]->nb_element_computation();
                observation_reestim[i][j]->max_computation();
                observation_reestim[i][j]->mean_computation();
                observation_reestim[i][j]->variance_computation(true);
//                observation_reestim[i][j]->variance_computation();
              }

              if (mixt->discrete_parametric_process[i]) {
                for (j = 0;j < mixt->nb_component;j++) {
                  hobservation->update(observation_reestim[i][j] ,
                                       MAX((int)(observation_reestim[i][j]->nb_element *
                                                 MAX(sqrt(observation_reestim[i][j]->variance) , 1.) * OBSERVATION_COEFF) , MIN_NB_ELEMENT));
                  observation_likelihood = hobservation->Reestimation<int>::type_parametric_estimation(mixt->discrete_parametric_process[i]->observation[j] ,
                                                                                                       0 , true , OBSERVATION_THRESHOLD);

                  if (observation_likelihood != D_INF) {
                    mixt->discrete_parametric_process[i]->observation[j]->computation(marginal_distribution[i]->nb_value ,
                                                                                      OBSERVATION_THRESHOLD);

                    if (mixt->discrete_parametric_process[i]->observation[j]->ident == BINOMIAL) {
                      for (k = mixt->discrete_parametric_process[i]->observation[j]->nb_value;k < marginal_distribution[i]->nb_value;k++) {
                        mixt->discrete_parametric_process[i]->observation[j]->mass[k] = 0.;
                      }
                    }
                  }
                }
              }

              else {
                switch (mixt->continuous_parametric_process[i]->ident) {

                case GAMMA : {
                  for (j = 0;j < mixt->nb_component;j++) {
                    observation_reestim[i][j]->gamma_estimation(mixt->continuous_parametric_process[i]->observation[j] , iter);
                  }
                  break;
                }

                case GAUSSIAN : {
                  for (j = 0;j < mixt->nb_component;j++) {
                    mixt->continuous_parametric_process[i]->observation[j]->location = observation_reestim[i][j]->mean;
                  }

                  switch (common_dispersion) {

                  case false : {
                    for (j = 0;j < mixt->nb_component;j++) {
                      mixt->continuous_parametric_process[i]->observation[j]->dispersion = sqrt(observation_reestim[i][j]->variance);
                      if (mixt->continuous_parametric_process[i]->observation[j]->dispersion /
                          mixt->continuous_parametric_process[i]->observation[j]->location < GAUSSIAN_MIN_VARIATION_COEFF) {
                        mixt->continuous_parametric_process[i]->observation[j]->dispersion = mixt->continuous_parametric_process[i]->observation[j]->location * GAUSSIAN_MIN_VARIATION_COEFF;
                      }
                    }
                    break;
                  }

                  case true : {
                    variance = 0.;
                    buff = 0.;

                    for (j = 0;j < mixt->nb_component;j++) {
                      for (k = observation_reestim[i][j]->offset;k < observation_reestim[i][j]->nb_value;k++) {
                        diff = k - observation_reestim[i][j]->mean;
                        variance += observation_reestim[i][j]->frequency[k] * diff * diff;
                      }

                      buff += observation_reestim[i][j]->nb_element;
                    }

                    variance /= buff;
//                    variance /= (buff - 1);

                    for (j = 0;j < mixt->nb_component;j++) {
                      mixt->continuous_parametric_process[i]->observation[j]->dispersion = sqrt(variance);
                    }
                    break;
                  }
                  }

                  break;
                }

                case VON_MISES : {
                  for (j = 0;j < mixt->nb_component;j++) {
                    observation_reestim[i][j]->mean_direction_computation(mean_direction[j]);
                    mixt->continuous_parametric_process[i]->observation[j]->location = mean_direction[j][3];
                  }

                  switch (common_dispersion) {

                  case false : {
                    for (j = 0;j < mixt->nb_component;j++) {
                      mixt->continuous_parametric_process[i]->observation[j]->dispersion = von_mises_concentration_computation(mean_direction[j][2]);
                    }
                    break;
                  }

                  case true : {
                    global_mean_direction = 0.;
                    buff = 0.;

                    for (j = 0;j < mixt->nb_component;j++) {
                      global_mean_direction += observation_reestim[i][j]->nb_element * mean_direction[j][2];
                      buff += observation_reestim[i][j]->nb_element;
                    }
                    concentration = von_mises_concentration_computation(global_mean_direction / buff);

                    for (j = 0;j < mixt->nb_component;j++) {
                      mixt->continuous_parametric_process[i]->observation[j]->dispersion = concentration;
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
              if (variance_factor != INDEPENDENT) {
                switch (mixt->continuous_parametric_process[i]->ident) {
                case GAMMA :
                  tied_gamma_estimation(component_vector_count , i ,
                                        mixt->continuous_parametric_process[i] ,
                                        variance_factor , iter);
                  break;
                case GAUSSIAN :
                  tied_gaussian_estimation(component_vector_count , i ,
                                           mixt->continuous_parametric_process[i] ,
                                           variance_factor);
                  break;
                }
              }

              else {
                switch (mixt->continuous_parametric_process[i]->ident) {
                case GAMMA :
                  gamma_estimation(component_vector_count , i ,
                                   mixt->continuous_parametric_process[i] , iter);
                  break;
                case GAUSSIAN :
                  gaussian_estimation(component_vector_count , i ,
                                      mixt->continuous_parametric_process[i]);
                  break;
                case VON_MISES :
                  von_mises_estimation(component_vector_count , i ,
                                       mixt->continuous_parametric_process[i]);
                  break;
                }
              }
            }
          }
        }
      }

#     ifdef MESSAGE
      os << STAT_label[STATL_ITERATION] << " " << iter << "   "
         << STAT_label[STATL_LIKELIHOOD] << ": " << likelihood << endl;
#     endif

#     ifdef DEBUG
      if (iter % 5 == 0) {
        cout << *mixt;
      }
#     endif

    }
    while ((likelihood != D_INF) && (((nb_iter == I_DEFAULT) && (iter < MIXTURE_NB_ITER) &&
             ((likelihood - previous_likelihood) / -likelihood > MIXTURE_LIKELIHOOD_DIFF)) ||
            ((nb_iter != I_DEFAULT) && (iter < nb_iter))));

    if (likelihood != D_INF) {

#     ifdef MESSAGE
      os << "\n" << iter << " " << STAT_label[STATL_ITERATIONS] << endl;
#     endif

      reestimation(mixt->nb_component , weight_reestim->frequency ,
                   mixt->weight->mass , MIN_PROBABILITY , true);

      // reestimation des lois d'observation categorielles

      if (!known_component) {
        for (i = 0;i < mixt->nb_output_process;i++) {
          if (mixt->categorical_process[i]) {
            for (j = 0;j < mixt->nb_component;j++) {
              reestimation(marginal_distribution[i]->nb_value , observation_reestim[i][j]->frequency ,
                           mixt->categorical_process[i]->observation[j]->mass ,
                           MIN_PROBABILITY , true);
            }
          }

          else if (mixt->discrete_parametric_process[i]) {
            mixt->discrete_parametric_process[i]->nb_value_computation();
          }
        }
      }
    }

    // destruction des structures de donnees de l'algorithme

    delete [] component_proba;
    delete [] component_cumul;

    delete weight_reestim;

    if (!known_component) {
      for (i = 0;i < mixt->nb_output_process;i++) {
        if (observation_reestim[i]) {
          for (j = 0;j < mixt->nb_component;j++) {
            delete observation_reestim[i][j];
          }
          delete [] observation_reestim[i];
        }
      }
      delete [] observation_reestim;

      delete hobservation;

      if (component_vector_count) {
        for (i = 0;i < nb_vector;i++) {
          delete [] component_vector_count[i];
        }
        delete [] component_vector_count;
      }

      if (mean_direction) {
        for (i = 0;i < mixt->nb_component;i++) {
          delete [] mean_direction[i];
        }
        delete [] mean_direction;
      }
    }

    if (likelihood == D_INF) {
      delete mixt;
      mixt = NULL;
      error.update(STAT_error[STATR_ESTIMATION_FAILURE]);
    }

    else {
      if (assignment) {
        mixt->mixture_data = new MixtureData(*this , 'a');
        vec = mixt->mixture_data;

        mixt->individual_assignment(*vec , true);

#       ifdef MESSAGE
        os << "\n" << STAT_label[STATL_CLASSIFICATION_LIKELIHOOD] << ": " << vec->restoration_likelihood;

        for (i = 0;i < nb_variable;i++) {
          if (type[i] == REAL_VALUE) {
            break;
          }
        }
        if (i == nb_variable) {
          os << " | " << mixt->classification_likelihood_computation(*vec)
             << " (" << vec->classification_information_computation() << ")" << endl;
        }
        cout << endl;
#       endif

        // calcul des melanges de lois d'observation (poids deduits de la restauration)

        for (i = 0;i < mixt->nb_output_process;i++) {
          if (mixt->categorical_process[i]) {
            mixt->categorical_process[i]->restoration_weight = new Distribution(*(vec->marginal_distribution[0]));
            mixt->categorical_process[i]->restoration_mixture = mixt->categorical_process[i]->mixture_computation(mixt->categorical_process[i]->restoration_weight);
          }

          else if (mixt->discrete_parametric_process[i]) {
            mixt->discrete_parametric_process[i]->restoration_weight = new Distribution(*(vec->marginal_distribution[0]));
            mixt->discrete_parametric_process[i]->restoration_mixture = mixt->discrete_parametric_process[i]->mixture_computation(mixt->discrete_parametric_process[i]->restoration_weight);
          }

          else if (mixt->continuous_parametric_process[i]) {
            mixt->continuous_parametric_process[i]->restoration_weight = new Distribution(*(vec->marginal_distribution[0]));
          }
        }
      }

      else {
        mixt->mixture_data = new MixtureData(*this);
        vec = mixt->mixture_data;
        if (vec->type[0] == STATE) {
          vec->state_variable_init(INT_VALUE);
        }
      }

      mixt->weight->cumul_computation();
      mixt->weight->max_computation();
      mixt->weight->mean_computation();
      mixt->weight->variance_computation();

      if (!known_component) {
        for (i = 0;i < mixt->nb_output_process;i++) {
          if (mixt->categorical_process[i]) {
            for (j = 0;j < mixt->nb_component;j++) {
              mixt->categorical_process[i]->observation[j]->cumul_computation();

              mixt->categorical_process[i]->observation[j]->max_computation();
//              mixt->categorical_process[i]->observation[j]->mean_computation();
//              mixt->categorical_process[i]->observation[j]->variance_computation();
            }
          }
        }
      }

      // calcul de la vraisemblance du modele

      vec->likelihood = mixt->likelihood_computation(*this);

#     ifdef DEBUG
//      cout << *mixt;
      cout << "iteration " << iter << "  "
           << STAT_label[STATL_LIKELIHOOD] << ": " << vec->likelihood << endl;
#     endif

#     ifdef MESSAGE
      if  ((assignment) && (vec->nb_vector <= POSTERIOR_PROBABILITY_NB_VECTOR)) {
        os << "\n" << STAT_label[STATL_POSTERIOR_ASSIGNMENT_PROBABILITY] << endl;
        for (i = 0;i < vec->nb_vector;i++) {
          os << STAT_label[STATL_VECTOR] << " " << vec->identifier[i] << ": "
             << vec->posterior_probability[i] << endl;
        }
      }
#     endif

      // calcul des melanges de lois d'observation (poids estimes)

      for (i = 0;i < mixt->nb_output_process;i++) {
        if (mixt->categorical_process[i]) {
          mixt->categorical_process[i]->weight = new Distribution(*(mixt->weight));
          mixt->categorical_process[i]->mixture = mixt->categorical_process[i]->mixture_computation(mixt->categorical_process[i]->weight);
        }

        else if (mixt->discrete_parametric_process[i]) {
          mixt->discrete_parametric_process[i]->weight = new Distribution(*(mixt->weight));
          mixt->discrete_parametric_process[i]->mixture = mixt->discrete_parametric_process[i]->mixture_computation(mixt->discrete_parametric_process[i]->weight);
        }

        else if (mixt->continuous_parametric_process[i]) {
          mixt->continuous_parametric_process[i]->weight = new Distribution(*(mixt->weight));
        }
      }
    }
  }

  return mixt;
}


/*--------------------------------------------------------------*
 *
 *  Estimation des parametres d'un melange gamma ou gaussien univarie a parametres lies
 *  a partir d'un echantillon de vecteurs par l'algorithme MCEM.
 *
 *  arguments : reference sur un objet StatError, stream, nombre de composantes,
 *              identificateur des composantes (GAMMA / GAUSSIAN),
 *              moyenne et parametre de forme de la premiere composante (GAMMA) /
 *              moyenne et ecart-type de la premiere composante (GAUSSIAN), flag moyennes liees,
 *              type de lien entre les variances (CONVOLUTION_FACTOR / SCALING_FACTOR)
 *              parametres pour le nombre d'affectation des individus simulees,
 *              flag sur le calcul des affectations optimales, nombre d'iterations.
 *
 *--------------------------------------------------------------*/

Mixture* Vectors::mixture_stochastic_estimation(StatError &error , ostream &os , int nb_component ,
                                                int ident , double mean , double standard_deviation ,
                                                bool tied_mean , int variance_factor ,
                                                int min_nb_assignment , int max_nb_assignment ,
                                                double parameter , bool assignment , int nb_iter) const

{
  register int i , j;
  Mixture *imixt , *mixt;


  error.init();

  if ((nb_component < 2) || (nb_component > MIXTURE_NB_COMPONENT)) {
    mixt = NULL;
    error.update(STAT_error[STATR_NB_DISTRIBUTION]);
  }

  else {
    imixt = new Mixture(nb_component , ident , mean , standard_deviation , tied_mean , variance_factor);

    mixt = mixture_stochastic_estimation(error , os , *imixt , false , false , variance_factor ,
                                         min_nb_assignment , max_nb_assignment , parameter ,
                                         assignment , nb_iter);
  }

  return mixt;
}


/*--------------------------------------------------------------*
 *
 *  Affectation des individus aux composantes les plus probables
 *  pour un melange de lois multivariees et calcul des probabilites a posteriori et
 *  des entropies associees.
 *
 *  arguments : reference sur un objet MixtureData, flag affectation
 *              des individus aux composantes.
 *
 *--------------------------------------------------------------*/

void Mixture::individual_assignment(MixtureData &vec , bool assignment) const

{
  register int i , j , k;
  double norm , *component_proba;


  vec.posterior_probability = new double[vec.nb_vector];
  vec.entropy = new double[vec.nb_vector];
  vec.restoration_likelihood = 0.;

  component_proba = new double[nb_component];

  vec.sample_entropy = 0.;

  for (i = 0;i < vec.nb_vector;i++) {
    norm = 0.;

    for (j = 0;j < nb_component;j++) {
      component_proba[j] = 1.;

      for (k = 0;k < nb_output_process;k++) {
        if (categorical_process[k]) {
          component_proba[j] *= categorical_process[k]->observation[j]->mass[vec.int_vector[i][k + 1]];
        }

        else if (discrete_parametric_process[k]) {
          component_proba[j] *= discrete_parametric_process[k]->observation[j]->mass[vec.int_vector[i][k + 1]];
        }

        else {
          if ((continuous_parametric_process[k]->ident == GAMMA) && (vec.min_value[k + 1] < vec.min_interval[k + 1] / 2)) {
            switch (vec.type[k + 1]) {
            case INT_VALUE :
              component_proba[j] *= continuous_parametric_process[k]->observation[j]->mass_computation(vec.int_vector[i][k + 1] , vec.int_vector[i][k + 1] + vec.min_interval[k + 1]);
              break;
            case REAL_VALUE :
              component_proba[j] *= continuous_parametric_process[k]->observation[j]->mass_computation(vec.real_vector[i][k + 1] , vec.real_vector[i][k + 1] + vec.min_interval[k + 1]);
              break;
            }
          }

          else {
            switch (vec.type[k + 1]) {
            case INT_VALUE :
              component_proba[j] *= continuous_parametric_process[k]->observation[j]->mass_computation(vec.int_vector[i][k + 1] - vec.min_interval[k + 1] / 2 , vec.int_vector[i][k + 1] + vec.min_interval[k + 1] / 2);
              break;
            case REAL_VALUE :
              component_proba[j] *= continuous_parametric_process[k]->observation[j]->mass_computation(vec.real_vector[i][k + 1] - vec.min_interval[k + 1] / 2 , vec.real_vector[i][k + 1] + vec.min_interval[k + 1] / 2);
              break;
            }
          }
        }
      }

      component_proba[j] *= weight->mass[j];
      norm += component_proba[j];
    }
 
    switch (assignment) {

    case false : {
      vec.posterior_probability[i] = component_proba[vec.int_vector[i][0]];
      break;
    }

    case true : {
      vec.posterior_probability[i] = 0.;
      for (j = 0;j < nb_component;j++) {
        if (component_proba[j] > vec.posterior_probability[i]) {
          vec.posterior_probability[i] = component_proba[j];
          vec.int_vector[i][0] = j;
        }
      }
      break;
    }
    }

    if (vec.posterior_probability[i] > 0.) {
      vec.restoration_likelihood += log(vec.posterior_probability[i]);
    }

    if (norm > 0.) {
      vec.posterior_probability[i] /= norm;

      vec.entropy[i] = 0.;
      for (j = 0;j < nb_component;j++) {
        if (component_proba[j] > 0.) {
          component_proba[j] /= norm;
          vec.entropy[i] -= component_proba[j] * log(component_proba[j]);
        }
      }

      vec.sample_entropy += vec.entropy[i];
    }
  }

  if (assignment) {
    vec.min_value_computation(0);
    vec.max_value_computation(0);

    vec.build_marginal_frequency_distribution(0);

    vec.build_observation_frequency_distribution(nb_component);
    vec.build_observation_histogram(nb_component);
  }

  delete [] component_proba;
}


/*--------------------------------------------------------------*
 *
 *  Simulation par un melange de lois.
 *
 *  arguments : reference sur un objet StatError, effectif.
 *
 *--------------------------------------------------------------*/

MixtureData* Mixture::simulation(StatError &error , int nb_vector) const

{
  register int i , j , k;
  int *itype , *decimal_scale;
  double buff , min_location;
  Mixture *mixt;
  Vectors *observed_vec;
  MixtureData *vec;


  error.init();

  if ((nb_vector < 1) || (nb_vector > MIXTURE_NB_VECTOR)) {
    vec = NULL;
    error.update(STAT_error[STATR_SAMPLE_SIZE]);
  }

  else {

    // creation d'un objet MixtureData

    itype = new int[nb_output_process + 1];

    itype[0] = STATE;
    for (i = 0;i < nb_output_process;i++) {
      if (!continuous_parametric_process[i]) {
        itype[i + 1] = INT_VALUE;
      }
      else {
        itype[i + 1] = REAL_VALUE;
      }
    }

    vec = new MixtureData(nb_vector , nb_output_process + 1 , itype);
    delete [] itype;

    vec->mixture = new Mixture(*this , false);
    mixt = vec->mixture;

    decimal_scale = new int[mixt->nb_output_process];

    for (i = 0;i < mixt->nb_output_process;i++) {
      if (mixt->continuous_parametric_process[i]) {
        switch (mixt->continuous_parametric_process[i]->ident) {

        case GAMMA : {
          min_location = mixt->continuous_parametric_process[i]->observation[0]->location * mixt->continuous_parametric_process[i]->observation[0]->dispersion;
          for (j = 1;j < mixt->nb_component;j++) {
            buff = mixt->continuous_parametric_process[i]->observation[j]->location * mixt->continuous_parametric_process[i]->observation[0]->dispersion;
            if (buff < min_location) {
              min_location = buff;
            }
          }

          buff = (int)ceil(log(min_location) / log(10));
          if (buff < GAMMA_MAX_NB_DECIMAL) {
            decimal_scale[i] = pow(10 , (GAMMA_MAX_NB_DECIMAL - buff));
          }
          else {
           decimal_scale[i] = 1;
          }

#         ifdef MESSAGE
          cout << "\nScale: " << i + 1 << " " << decimal_scale[i] << endl;
#         endif

          break;
        }

        case GAUSSIAN : {
          min_location = fabs(mixt->continuous_parametric_process[i]->observation[0]->location);
          for (j = 1;j < mixt->nb_component;j++) {
            buff = fabs(mixt->continuous_parametric_process[i]->observation[j]->location);
            if (buff < min_location) {
              min_location = buff;
            }
          }

          buff = (int)ceil(log(min_location) / log(10));
          if (buff < GAUSSIAN_MAX_NB_DECIMAL) {
            decimal_scale[i] = pow(10 , (GAUSSIAN_MAX_NB_DECIMAL - buff));
          }
          else {
            decimal_scale[i] = 1;
          }

#         ifdef MESSAGE
          cout << "\nScale: " << i + 1 << " " << decimal_scale[i] << endl;
#         endif

          break;
        }

        case VON_MISES : {
          switch (mixt->continuous_parametric_process[i]->unit) {
          case DEGREE :
            decimal_scale[i] = DEGREE_DECIMAL_SCALE;
            break;
          case RADIAN :
            decimal_scale[i] = RADIAN_DECIMAL_SCALE;
            break;
          }

          for (j = 0;j < mixt->nb_component;j++) {
            mixt->continuous_parametric_process[i]->observation[j]->von_mises_cumul_computation();
          }
          break;
        }
        }
      }
    }

    for (i = 0;i < vec->nb_vector;i++) {

      // ponderation

      j = mixt->weight->simulation();
      vec->int_vector[i][0] = j;

      // composantes

      for (k = 0;k < mixt->nb_output_process;k++) {
        if (mixt->categorical_process[k]) {
          vec->int_vector[i][k + 1] = mixt->categorical_process[k]->observation[j]->simulation();
        }
        else if (mixt->discrete_parametric_process[k]) {
          vec->int_vector[i][k + 1] = mixt->discrete_parametric_process[k]->observation[j]->simulation();
        }
        else if (mixt->continuous_parametric_process[k]) {
          vec->real_vector[i][k + 1] = round(mixt->continuous_parametric_process[k]->observation[j]->simulation() * decimal_scale[k]) / decimal_scale[k];
        }
      }
    }

    delete [] decimal_scale;

    for (i = 0;i < mixt->nb_output_process;i++) {
      if ((mixt->continuous_parametric_process[i]) &&
          (mixt->continuous_parametric_process[i]->ident == VON_MISES)) {
        for (j = 0;j < mixt->nb_component;j++) {
          delete [] mixt->continuous_parametric_process[i]->observation[j]->cumul;
          mixt->continuous_parametric_process[i]->observation[j]->cumul = NULL;
        }
      }
    }

    // extraction des caracteristiques des vecteurs simules

    vec->min_value_computation(0);
    vec->max_value_computation(0);
    vec->build_marginal_frequency_distribution(0);

    for (i = 1;i < vec->nb_variable;i++) {
      vec->min_value_computation(i);
      vec->max_value_computation(i);

      vec->build_marginal_frequency_distribution(i);
      vec->min_interval_computation(i);
    }

    vec->build_observation_frequency_distribution(mixt->nb_component);
    vec->build_observation_histogram(mixt->nb_component);

    // calcul des vraisemblances, des probabilites a posteriori et des entropies

    mixt->individual_assignment(*vec , false);

    observed_vec = vec->remove_variable_1();
    vec->likelihood = mixt->likelihood_computation(*observed_vec);
    delete observed_vec;

    // calcul des melanges de lois d'observation (poids deduits de la restauration)

    for (i = 0;i < mixt->nb_output_process;i++) {
      if (mixt->categorical_process[i]) {
        delete mixt->categorical_process[i]->restoration_weight;
        delete mixt->categorical_process[i]->restoration_mixture;
        mixt->categorical_process[i]->restoration_weight = new Distribution(*(vec->marginal_distribution[0]));
        mixt->categorical_process[i]->restoration_mixture = mixt->categorical_process[i]->mixture_computation(mixt->categorical_process[i]->restoration_weight);
      }

      else if (mixt->discrete_parametric_process[i]) {
        delete mixt->discrete_parametric_process[i]->restoration_weight;
        delete mixt->discrete_parametric_process[i]->restoration_mixture;
        mixt->discrete_parametric_process[i]->restoration_weight = new Distribution(*(vec->marginal_distribution[0]));
        mixt->discrete_parametric_process[i]->restoration_mixture = mixt->discrete_parametric_process[i]->mixture_computation(mixt->discrete_parametric_process[i]->restoration_weight);
      }

      else if (mixt->continuous_parametric_process[i]) {
        delete mixt->continuous_parametric_process[i]->restoration_weight;
        mixt->continuous_parametric_process[i]->restoration_weight = new Distribution(*(vec->marginal_distribution[0]));
      }
    }
  }

  return vec;
}


};  // namespace stat_tool
