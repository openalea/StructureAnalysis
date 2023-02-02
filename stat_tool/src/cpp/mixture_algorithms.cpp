/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       V-Plants: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2017 CIRAD/INRA/Inria Virtual Plants
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

#include "mixture.h"
#include "stat_label.h"

#include "stat_tool/distribution_reestimation.hpp"   // problem compiler C++ Windows

using namespace std;
using namespace boost::math;


namespace stat_tool {



/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the concentration parameter of a von Mises distribution
 *         on the basis of the mean direction (Mardia & Jupp, 2000; pp. 85-86).
 *
 *  \param[in] mean_direction mean direction.
 *
 *  \return                   concentration parameter.
 */
/*--------------------------------------------------------------*/

double von_mises_concentration_computation(double mean_direction)

{
  int i;
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


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the information quantity of a MixtureData object
 *         in the case of discrete observed variables.
 *
 *  \return information quantity.
 */
/*--------------------------------------------------------------*/

double MixtureData::classification_information_computation() const

{
  int i , j;
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


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the log-likelihood of a mixture of multivariate
 *         distributions for individuals assigned to components.
 *
 *  \param[in] vec reference on a MixtureData object.
 *
 *  \return        log-likelihood.
 */
/*--------------------------------------------------------------*/

double Mixture::classification_likelihood_computation(const MixtureData &vec) const

{
  int i , j;
  int nb_value;
  double buff , likelihood = 0.;


  // checking of the compatibility between the model and the data

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


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the log-likelihood of a mixture of
           multivariate distributions for a sample of observed vectors.
 *
 *  \param[in] vec   reference on a Vectors object,
 *  \param[in] index vector index.
 *
 *  \return          log-likelihood.
 */
/*--------------------------------------------------------------*/

double Mixture::likelihood_computation(const Vectors &vec , int index) const

{
  int i , j , k;
  int nb_value;
  double likelihood = 0. , norm , component_proba;


  // checking of the compatibility between the model and the data

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
              if (((continuous_parametric_process[k]->ident == GAMMA) ||
                   (continuous_parametric_process[k]->ident == INVERSE_GAUSSIAN)) && (vec.min_value[k] < vec.min_interval[k] / 2)) {
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


/*--------------------------------------------------------------*/
/**
 *  \brief Estimation of a multivariate mixture of distributions using the EM algorithm.
 *
 *  \param[in] error             reference on a StatError object,
 *  \param[in] display           flag for displaying estimation intermediate results,
 *  \param[in] imixt             initial mixture,
 *  \param[in] known_component   flags component estimation,
 *  \param[in] common_dispersion common dispersion parameter (continuous observation processes),
 *  \param[in] variance_factor   type of relationship between variances (CONVOLUTION_FACTOR/SCALING_FACTOR)
 *  \param[in] assignment        flag on the computation of the optimal assignments,
 *  \param[in] nb_iter           number of iterations.
 *
 *  \return                      Mixture object.
 */
/*--------------------------------------------------------------*/

Mixture* Vectors::mixture_estimation(StatError &error , bool display , const Mixture &imixt ,
                                     bool known_component , bool common_dispersion ,
                                     tying_rule variance_factor , bool assignment , int nb_iter) const

{
  bool status;
  int i , j , k;
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

  // test number of values per variable

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

#     ifdef MESSAGE
      cout << STAT_label[STATL_VARIABLE] << " " << i + 1 << " - " << STAT_label[STATL_MIN_VALUE] << ": " << min_value[i]
           << "   minimum interval: " << min_interval[i] << endl;
#     endif

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

    // construction of a Mixture object

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

    // construction of the data structures of the algorithm
 
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

      // initialization of the reestimation quantities

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

        // computation of the component posterior probabilities

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
              if (((mixt->continuous_parametric_process[k]->ident == GAMMA) ||
                   (mixt->continuous_parametric_process[k]->ident == INVERSE_GAUSSIAN)) && (min_value[k] < min_interval[k] / 2)) {
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

            // accumulation of the reestimation quantities of the weights

            weight_reestim->frequency[j] += component_proba[j];

            // accumulation of the reestimation quantities of the observation distributions

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

        // reestimation of the weights

        reestimation(mixt->nb_component , weight_reestim->frequency ,
                     mixt->weight->mass , MIN_PROBABILITY , false);

        // reestimation of the observation distributions

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

                case INVERSE_GAUSSIAN : {
                  for (j = 0;j < mixt->nb_component;j++) {
                    observation_reestim[i][j]->inverse_gaussian_estimation(mixt->continuous_parametric_process[i]->observation[j]);
                  }
                  break;
                }

                case GAUSSIAN : {
                  for (j = 0;j < mixt->nb_component;j++) {
                    mixt->continuous_parametric_process[i]->observation[j]->location = observation_reestim[i][j]->mean;
                  }

                  if (common_dispersion) {
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
                  }

                  else {
                    for (j = 0;j < mixt->nb_component;j++) {
                      mixt->continuous_parametric_process[i]->observation[j]->dispersion = sqrt(observation_reestim[i][j]->variance);
                      if (mixt->continuous_parametric_process[i]->observation[j]->dispersion /
                          mixt->continuous_parametric_process[i]->observation[j]->location < GAUSSIAN_MIN_VARIATION_COEFF) {
                        mixt->continuous_parametric_process[i]->observation[j]->dispersion = mixt->continuous_parametric_process[i]->observation[j]->location * GAUSSIAN_MIN_VARIATION_COEFF;
                      }
                    }
                  }

                  break;
                }

                case VON_MISES : {
                  for (j = 0;j < mixt->nb_component;j++) {
                    observation_reestim[i][j]->mean_direction_computation(mean_direction[j]);
                    mixt->continuous_parametric_process[i]->observation[j]->location = mean_direction[j][3];
                  }

                  if (common_dispersion) {
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
                  }

                  else {
                    for (j = 0;j < mixt->nb_component;j++) {
                      mixt->continuous_parametric_process[i]->observation[j]->dispersion = von_mises_concentration_computation(mean_direction[j][2]);
                    }
                  }
                  break;
                }
                }
              }
            }

            else {
              if (mixt->continuous_parametric_process[i]->offset > 0.) {
                tied_gaussian_estimation(component_vector_count , i ,
                                         mixt->continuous_parametric_process[i]);
              }

              else if (variance_factor != INDEPENDENT) {
                switch (mixt->continuous_parametric_process[i]->ident) {
                case GAMMA :
                  tied_gamma_estimation(component_vector_count , i ,
                                        mixt->continuous_parametric_process[i] ,
                                        variance_factor , iter);
                  break;
                case INVERSE_GAUSSIAN :
                  tied_inverse_gaussian_estimation(component_vector_count , i ,
                                                   mixt->continuous_parametric_process[i] ,
                                                   variance_factor);
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
                case INVERSE_GAUSSIAN :
                  inverse_gaussian_estimation(component_vector_count , i ,
                                              mixt->continuous_parametric_process[i]);
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

      if (display) {
        cout << STAT_label[STATL_ITERATION] << " " << iter << "   "
             << STAT_label[STATL_LIKELIHOOD] << ": " << likelihood << endl;
      }

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
      if (display) {
        cout << "\n" << iter << " " << STAT_label[STATL_ITERATIONS] << endl;
      }

      reestimation(mixt->nb_component , weight_reestim->frequency ,
                   mixt->weight->mass , MIN_PROBABILITY , true);

      // reestimation of the categorical observation distributions

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

    // destruction of the data structures of the algorithm

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
        mixt->mixture_data = new MixtureData(*this , ADD_COMPONENT_VARIABLE);
        vec = mixt->mixture_data;

        mixt->individual_assignment(*vec , true);

        if (display) {
          cout << "\n" << STAT_label[STATL_CLASSIFICATION_LIKELIHOOD] << ": " << vec->restoration_likelihood;

          for (i = 0;i < nb_variable;i++) {
            if (type[i] == REAL_VALUE) {
              break;
            }
          }
          if (i == nb_variable) {
            cout << " | " << mixt->classification_likelihood_computation(*vec)
                 << " (" << vec->classification_information_computation() << ")" << endl;
          }
          cout << endl;
        }

        // computation of the mixtures of observation distributions (weights deduced from the restoration)

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

      // computation of the log-likelihood of the model

      vec->likelihood = mixt->likelihood_computation(*this);

#     ifdef DEBUG
//      cout << *mixt;
      cout << "iteration " << iter << "  "
           << STAT_label[STATL_LIKELIHOOD] << ": " << vec->likelihood << endl;
#     endif

      if ((display) && (assignment) && (vec->nb_vector <= POSTERIOR_PROBABILITY_NB_VECTOR)) {
        cout << "\n" << STAT_label[STATL_POSTERIOR_ASSIGNMENT_PROBABILITY] << endl;
        for (i = 0;i < vec->nb_vector;i++) {
          cout << STAT_label[STATL_VECTOR] << " " << vec->identifier[i] << ": "
               << vec->posterior_probability[i] << endl;
        }
      }

      // computation of the mixtures of observation distributions (estimated weights)

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

#     ifdef MESSAGE
      if ((nb_variable == 1) && (mixt->nb_component == 2) && (mixt->continuous_parametric_process[0])) {
        continuous_parametric_estimation(0 , mixt->continuous_parametric_process[0]->ident);
      }
#     endif

    }
  }

  return mixt;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Estimation of a mixture of univariate Gaussian distributions with
 *         evenly spaced means using the EM algorithm.
 *
 *  \param[in] error              reference on a StatError object,
 *  \param[in] display            flag for displaying estimation intermediate results,
 *  \param[in] nb_component       number of components,
 *  \param[in] offset             mixture offset,
 *  \param[in] mean               mean - offset of the 1st component,
 *  \param[in] standard_deviation standard deviation,
 *  \param[in] common_dispersion  common dispersion parameter,
 *  \param[in] assignment         flag on the computation of the optimal assignments,
 *  \param[in] nb_iter            number of iterations.
 *
 *  \return                       Mixture object.
 */
/*--------------------------------------------------------------*/

Mixture* Vectors::mixture_estimation(StatError &error , bool display , int nb_component ,
                                     double offset , double mean , double standard_deviation ,
                                     bool common_dispersion , bool assignment , int nb_iter) const

{
  int i , j;
  Mixture *imixt , *mixt;


  error.init();

  if ((nb_component < 2) || (nb_component > MIXTURE_NB_COMPONENT)) {
    mixt = NULL;
    error.update(STAT_error[STATR_NB_DISTRIBUTION]);
  }

  else {
    imixt = new Mixture(nb_component , offset , mean , standard_deviation , common_dispersion);

    mixt = mixture_estimation(error , display , *imixt , false , common_dispersion ,
                              INDEPENDENT , assignment , nb_iter);
  }

  return mixt;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Estimation of a mixture of univariate gamma, inverse Gaussian or
 *         Gaussian distributions with tied parameters using the EM algorithm.
 *
 *  \param[in] error              reference on a StatError object,
 *  \param[in] display            flag for displaying estimation intermediate results,
 *  \param[in] nb_component       number of components,
 *  \param[in] ident              component identifiers,
 *  \param[in] mean               mean of the 1st component,
 *  \param[in] standard_deviation shape parameter (gamma) / scale parameter (inverse Gaussian) /
 *                                standard deviation (Gaussian) of the 1st component,
 *  \param[in] tied_mean          flag tied means,
 *  \param[in] variance_factor    type of relationship between variances (CONVOLUTION_FACTOR/SCALING_FACTOR)
 *  \param[in] assignment         flag on the computation of the optimal assignments,
 *  \param[in] nb_iter            number of iterations.
 *
 *  \return                       Mixture object.
 */
/*--------------------------------------------------------------*/

Mixture* Vectors::mixture_estimation(StatError &error , bool display , int nb_component ,
                                     int ident , double mean , double standard_deviation ,
                                     bool tied_mean , tying_rule variance_factor ,
                                     bool assignment , int nb_iter) const

{
  int i , j;
  Mixture *imixt , *mixt;


  error.init();

  if ((nb_component < 2) || (nb_component > MIXTURE_NB_COMPONENT)) {
    mixt = NULL;
    error.update(STAT_error[STATR_NB_DISTRIBUTION]);
  }

  else {
    imixt = new Mixture(nb_component , ident , mean , standard_deviation , tied_mean , variance_factor);

#   ifdef MESSAGE
    cout << *mixt;
#   endif

    mixt = mixture_estimation(error , display , *imixt , false , false ,
                              variance_factor , assignment , nb_iter);
  }

  return mixt;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Estimation of a multivariate mixture of distributions using the MCEM algorithm.
 *
 *  \param[in] error             reference on a StatError object,
 *  \param[in] display           flag for displaying estimation intermediate results,
 *  \param[in] imixt             initial mixture,
 *  \param[in] known_component   flags component estimation,
 *  \param[in] common_dispersion common dispersion parameter (continuous observation processes),
 *  \param[in] variance_factor   type of relationship between variances (CONVOLUTION_FACTOR/SCALING_FACTOR),
 *  \param[in] min_nb_assignment minimum number of assignments of generated individuals,
 *  \param[in] max_nb_assignment maximum number of assignments of generated individuals,
 *  \param[in] parameter         parameter for the assignments of the generated individuals,
 *  \param[in] assignment        flag on the computation of the optimal assignments,
 *  \param[in] nb_iter           number of iterations.
 *
 *  \return                      Mixture object.
 */
/*--------------------------------------------------------------*/

Mixture* Vectors::mixture_stochastic_estimation(StatError &error , bool display , const Mixture &imixt ,
                                                bool known_component , bool common_dispersion ,
                                                tying_rule variance_factor , int min_nb_assignment ,
                                                int max_nb_assignment , double parameter ,
                                                bool assignment , int nb_iter) const

{
  bool status;
  int i , j , k , m;
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

  // test number of values per variable

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

    // construction of a Mixture object

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

    // construction of the data structures of the algorithm
 
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

      // computation of the number of assignments of generated individuals

      if (min_nb_assignment + (int)::round(parameter * iter) < max_nb_assignment) {
        nb_assignment = min_nb_assignment + (int)::round(parameter * iter);
      }
      else {
        nb_assignment = max_nb_assignment;
      }

/*      nb_assignment = max_nb_assignment - (int)::round((max_nb_assignment - min_nb_assignment) *
                          exp(-parameter * iter)); */

      iter++;

      // initialization of the reestimation quantities

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

        // computation of the component posterior probabilities

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
              if (((mixt->continuous_parametric_process[k]->ident == GAMMA) ||
                   (mixt->continuous_parametric_process[k]->ident == INVERSE_GAUSSIAN)) && (min_value[k] < min_interval[k] / 2)) {
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

            // accumulation of the reestimation quantities of the weights

            (weight_reestim->frequency[k])++;

            // accumulation of the reestimation quantities of the observation distributions

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

        // reestimation of the weights

        reestimation(mixt->nb_component , weight_reestim->frequency ,
                     mixt->weight->mass , MIN_PROBABILITY , false);

        // reestimation of the observation distributions

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

                case INVERSE_GAUSSIAN : {
                  for (j = 0;j < mixt->nb_component;j++) {
                    observation_reestim[i][j]->inverse_gaussian_estimation(mixt->continuous_parametric_process[i]->observation[j]);
                  }
                  break;
                }

                case GAUSSIAN : {
                  for (j = 0;j < mixt->nb_component;j++) {
                    mixt->continuous_parametric_process[i]->observation[j]->location = observation_reestim[i][j]->mean;
                  }

                  if (common_dispersion) {
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
                  }

                  else {
                    for (j = 0;j < mixt->nb_component;j++) {
                      mixt->continuous_parametric_process[i]->observation[j]->dispersion = sqrt(observation_reestim[i][j]->variance);
                      if (mixt->continuous_parametric_process[i]->observation[j]->dispersion /
                          mixt->continuous_parametric_process[i]->observation[j]->location < GAUSSIAN_MIN_VARIATION_COEFF) {
                        mixt->continuous_parametric_process[i]->observation[j]->dispersion = mixt->continuous_parametric_process[i]->observation[j]->location * GAUSSIAN_MIN_VARIATION_COEFF;
                      }
                    }
                  }

                  break;
                }

                case VON_MISES : {
                  for (j = 0;j < mixt->nb_component;j++) {
                    observation_reestim[i][j]->mean_direction_computation(mean_direction[j]);
                    mixt->continuous_parametric_process[i]->observation[j]->location = mean_direction[j][3];
                  }

                  if (common_dispersion) {
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
                  }

                  else {
                    for (j = 0;j < mixt->nb_component;j++) {
                      mixt->continuous_parametric_process[i]->observation[j]->dispersion = von_mises_concentration_computation(mean_direction[j][2]);
                    }
                  }
                  break;
                }
                }
              }
            }

            else {
              if (mixt->continuous_parametric_process[i]->offset > 0.) {
                tied_gaussian_estimation(component_vector_count , i ,
                                         mixt->continuous_parametric_process[i]);
              }

              else if (variance_factor != INDEPENDENT) {
                switch (mixt->continuous_parametric_process[i]->ident) {
                case GAMMA :
                  tied_gamma_estimation(component_vector_count , i ,
                                        mixt->continuous_parametric_process[i] ,
                                        variance_factor , iter);
                  break;
                case INVERSE_GAUSSIAN :
                  tied_inverse_gaussian_estimation(component_vector_count , i ,
                                                   mixt->continuous_parametric_process[i] ,
                                                   variance_factor);
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
                case INVERSE_GAUSSIAN :
                  inverse_gaussian_estimation(component_vector_count , i ,
                                              mixt->continuous_parametric_process[i]);
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

      if (display) {
        cout << STAT_label[STATL_ITERATION] << " " << iter << "   "
             << STAT_label[STATL_LIKELIHOOD] << ": " << likelihood << endl;
      }

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
      if (display) {
        cout << "\n" << iter << " " << STAT_label[STATL_ITERATIONS] << endl;
      }

      reestimation(mixt->nb_component , weight_reestim->frequency ,
                   mixt->weight->mass , MIN_PROBABILITY , true);

      // reestimation of the categorical observation distributions

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

    // destruction of the data structures of the algorithm

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
        mixt->mixture_data = new MixtureData(*this , ADD_COMPONENT_VARIABLE);
        vec = mixt->mixture_data;

        mixt->individual_assignment(*vec , true);

        if (display) {
          cout << "\n" << STAT_label[STATL_CLASSIFICATION_LIKELIHOOD] << ": " << vec->restoration_likelihood;

          for (i = 0;i < nb_variable;i++) {
            if (type[i] == REAL_VALUE) {
              break;
            }
          }
          if (i == nb_variable) {
            cout << " | " << mixt->classification_likelihood_computation(*vec)
                 << " (" << vec->classification_information_computation() << ")" << endl;
          }
          cout << endl;
        }

        // computation of the mixtures of observation distributions (weights deduced from the restoration)

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

      // computation of the log-likelihood of the model

      vec->likelihood = mixt->likelihood_computation(*this);

#     ifdef DEBUG
//      cout << *mixt;
      cout << "iteration " << iter << "  "
           << STAT_label[STATL_LIKELIHOOD] << ": " << vec->likelihood << endl;
#     endif

      if  ((display) && (assignment) && (vec->nb_vector <= POSTERIOR_PROBABILITY_NB_VECTOR)) {
        cout << "\n" << STAT_label[STATL_POSTERIOR_ASSIGNMENT_PROBABILITY] << endl;
        for (i = 0;i < vec->nb_vector;i++) {
          cout << STAT_label[STATL_VECTOR] << " " << vec->identifier[i] << ": "
               << vec->posterior_probability[i] << endl;
        }
      }

      // computation of the mixtures of observation distributions (estimated weights)

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


/*--------------------------------------------------------------*/
/**
 *  \brief Estimation of a mixture of univariate Gaussian distributions with
 *         evenly spaced means using the MCEM algorithm.
 *
 *  \param[in] error              reference on a StatError object,
 *  \param[in] display            flag for displaying estimation intermediate results,
 *  \param[in] nb_component       number of components,
 *  \param[in] offset             mixture offset,
 *  \param[in] mean               mean - offset of the 1st component,
 *  \param[in] standard_deviation standard deviation,
 *  \param[in] common_dispersion  common dispersion parameter,
 *  \param[in] min_nb_assignment  minimum number of assignments of generated individuals,
 *  \param[in] max_nb_assignment  maximum number of assignments of generated individuals,
 *  \param[in] parameter          parameter for the assignments of the generated individuals,
 *  \param[in] assignment         flag on the computation of the optimal assignments,
 *  \param[in] nb_iter            number of iterations.
 *
 *  \return                       Mixture object.
 */
/*--------------------------------------------------------------*/

Mixture* Vectors::mixture_stochastic_estimation(StatError &error , bool display , int nb_component ,
                                                double offset , double mean , double standard_deviation ,
                                                bool common_dispersion , int min_nb_assignment , int max_nb_assignment ,
                                                double parameter , bool assignment , int nb_iter) const

{
  int i , j;
  Mixture *imixt , *mixt;


  error.init();

  if ((nb_component < 2) || (nb_component > MIXTURE_NB_COMPONENT)) {
    mixt = NULL;
    error.update(STAT_error[STATR_NB_DISTRIBUTION]);
  }

  else {
    imixt = new Mixture(nb_component , offset , mean , standard_deviation , common_dispersion);

    mixt = mixture_stochastic_estimation(error , display , *imixt , false , common_dispersion , INDEPENDENT ,
                                         min_nb_assignment , max_nb_assignment , parameter ,
                                         assignment , nb_iter);
  }

  return mixt;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Estimation of a mixture of univariate gamma, inverse Gaussian or
 *         Gaussian distributions with tied parameters using the MCEM algorithm.
 *
 *  \param[in] error              reference on a StatError object,
 *  \param[in] display            flag for displaying estimation intermediate results,
 *  \param[in] nb_component       number of components,
 *  \param[in] ident              component identifiers,
 *  \param[in] mean               mean and of the 1st component,
 *  \param[in] standard_deviation shape parameter (gamma) / scale parameter (inverse Gaussian) /
 *                                standard deviation (Gaussian) of the 1st component,
 *  \param[in] tied_mean          flag tied means,
 *  \param[in] variance_factor    type of relationship between variances (CONVOLUTION_FACTOR/SCALING_FACTOR)
 *  \param[in] min_nb_assignment  minimum number of assignments of generated individuals,
 *  \param[in] max_nb_assignment  maximum number of assignments of generated individuals,
 *  \param[in] parameter          parameter for the assignments of the generated individuals,
 *  \param[in] assignment         flag on the computation of the optimal assignments,
 *  \param[in] nb_iter            number of iterations.
 *
 *  \return                       Mixture object.
 */
/*--------------------------------------------------------------*/

Mixture* Vectors::mixture_stochastic_estimation(StatError &error , bool display , int nb_component ,
                                                int ident , double mean , double standard_deviation ,
                                                bool tied_mean , tying_rule variance_factor ,
                                                int min_nb_assignment , int max_nb_assignment ,
                                                double parameter , bool assignment , int nb_iter) const

{
  int i , j;
  Mixture *imixt , *mixt;


  error.init();

  if ((nb_component < 2) || (nb_component > MIXTURE_NB_COMPONENT)) {
    mixt = NULL;
    error.update(STAT_error[STATR_NB_DISTRIBUTION]);
  }

  else {
    imixt = new Mixture(nb_component , ident , mean , standard_deviation , tied_mean , variance_factor);

    mixt = mixture_stochastic_estimation(error , display , *imixt , false , false , variance_factor ,
                                         min_nb_assignment , max_nb_assignment , parameter ,
                                         assignment , nb_iter);
  }

  return mixt;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Assignment of the individuals to the most probable component of
 *         a multivariate mixture of distributions and computation of
 *         the associated posterior probabilities and entropies
 *
 *  \param[in] vec        reference on a MixtureData object,
 *  \param[in] assignment flag on the assignment of individuals to components.
 */
/*--------------------------------------------------------------*/

void Mixture::individual_assignment(MixtureData &vec , bool assignment) const

{
  int i , j , k;
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
          if (((continuous_parametric_process[k]->ident == GAMMA) ||
               (continuous_parametric_process[k]->ident == INVERSE_GAUSSIAN)) && (vec.min_value[k + 1] < vec.min_interval[k + 1] / 2)) {
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
 
    if (assignment) {
      vec.posterior_probability[i] = 0.;
      for (j = 0;j < nb_component;j++) {
        if (component_proba[j] > vec.posterior_probability[i]) {
          vec.posterior_probability[i] = component_proba[j];
          vec.int_vector[i][0] = j;
        }
      }
    }

    else {
      vec.posterior_probability[i] = component_proba[vec.int_vector[i][0]];
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


/*--------------------------------------------------------------*/
/**
 *  \brief Simulation using a multivariate mixture of distributions.
 *
 *  \param[in] error     reference on a StatError object,
 *  \param[in] nb_vector sample size.
 *
 *  \return              MixtureData object.
 */
/*--------------------------------------------------------------*/

MixtureData* Mixture::simulation(StatError &error , int nb_vector) const

{
  int i , j , k;
  int *decimal_scale;
  variable_nature *itype;
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

    // construction of a MixtureData object

    itype = new variable_nature[nb_output_process + 1];

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

        case INVERSE_GAUSSIAN : {
          min_location = mixt->continuous_parametric_process[i]->observation[0]->location;
          for (j = 1;j < mixt->nb_component;j++) {
            buff = mixt->continuous_parametric_process[i]->observation[j]->location;
            if (buff < min_location) {
              min_location = buff;
            }
          }

          buff = (int)ceil(log(min_location) / log(10));
          if (buff < INVERSE_GAUSSIAN_MAX_NB_DECIMAL) {
            decimal_scale[i] = pow(10 , (INVERSE_GAUSSIAN_MAX_NB_DECIMAL - buff));
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

      // weight

      j = mixt->weight->simulation();
      vec->int_vector[i][0] = j;

      // component

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

    // extraction of the characteristics of the generated data

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

    // computation of the log-likelihoods, posterior probabilies and entropies

    mixt->individual_assignment(*vec , false);

    observed_vec = vec->remove_variable_1();
    vec->likelihood = mixt->likelihood_computation(*observed_vec);
    delete observed_vec;

    // computation of the mixtures of observation distributions (weights deduced from the restoration)

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
