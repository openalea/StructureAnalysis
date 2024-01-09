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
 *       $Id: convolution_algorithms.cpp 17986 2015-04-23 06:43:20Z guedon $
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
#include "convolution.h"
#include "stat_label.h"

using namespace std;


namespace stat_tool {



/*--------------------------------------------------------------*
 *
 *  Calcul d'un produit de convolution de lois.
 *
 *  arguments : nombre minimum de valeurs de chaque loi,
 *              seuil sur la fonction de repartition,
 *              flags pour calculer les lois elementaires.
 *
 *--------------------------------------------------------------*/

void Convolution::computation(int min_nb_value , double cumul_threshold , bool *dist_flag)

{
  register int i;


  // calcul des lois elementaires

  if (dist_flag) {
    for (i = 0;i < nb_distribution;i++) {
      if (dist_flag[i]) {
        distribution[i]->computation(min_nb_value , cumul_threshold);
      }
    }
  }

  // calcul de la loi resultante

  convolution(*distribution[0] , *distribution[1]);
  for (i = 2;i < nb_distribution;i++) {
    convolution(*this , *distribution[i]);
  }

  cumul_computation();

  max_computation();
  mean_computation();
  variance_computation();
}


/*--------------------------------------------------------------*
 *
 *  Calcul des quantites de reestimation (estimateur EM).
 *
 *  arguments : reference sur la loi empirique, pointeurs sur les produits
 *              de convolution "partiels" et sur les quantites de reestimation.
 *
 *--------------------------------------------------------------*/

void Convolution::expectation_step(const FrequencyDistribution &histo ,
                                   const Distribution **partial_convol ,
                                   Reestimation<double> **reestim) const

{
  register int i , j , k;
  int min , max;


  for (i = 0;i < nb_distribution;i++) {
    if (reestim[i]) {
      for (j = 0;j < reestim[i]->alloc_nb_value;j++) {
        reestim[i]->frequency[j] = 0.;
      }

      for (j = histo.offset;j < histo.nb_value;j++) {
        if ((histo.frequency[j] > 0) && (mass[j] > 0.)) {
          min = MAX(distribution[i]->offset , j - (partial_convol[i]->nb_value - 1));
          max = MIN(distribution[i]->nb_value - 1 , j - partial_convol[i]->offset);

          if (max >= min) {
            for (k = min;k <= max;k++) {
              reestim[i]->frequency[k] += histo.frequency[j] * distribution[i]->mass[k] *
                                          partial_convol[i]->mass[j - k] / mass[j];
            }
          }
        }
      }

      reestim[i]->nb_value_computation();
      reestim[i]->offset_computation();
      reestim[i]->nb_element_computation();
      reestim[i]->max_computation();
      reestim[i]->mean_computation();
      reestim[i]->variance_computation();

#     ifdef DEBUG
      cout << "\nquantites de reestimation loi " << i << " :" << *reestim[i] << endl;
#     endif

    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Estimation des parametres d'un produit de convolution de lois
 *  par l'algorithme EM.
 *
 *  arguments : reference sur un objet StatError, stream, references sur
 *              la loi connue et sur la loi inconnue, type d'estimateur (vraisemblance,
 *              vraisemblance penalisee ou estimation d'une loi parametrique),
 *              nombre d'iterations, poids de la penalisation, type de penalisation,
 *              type de gestion des effets de bord (zero a l'exterieur du support ou
 *              prolongation de la loi).
 *
 *--------------------------------------------------------------*/

Convolution* FrequencyDistribution::convolution_estimation(StatError &error , ostream &os ,
                                                           const DiscreteParametric &known_dist ,
                                                           const DiscreteParametric &unknown_dist ,
                                                           int estimator , int nb_iter , double weight ,
                                                           int penalty_type , int outside) const

{
  bool status = true , dist_flag[2];
  register int i;
  int nb_likelihood_decrease;
  double likelihood , previous_likelihood , hlikelihood , *penalty;
  const Distribution *partial_convol[2];
  Reestimation<double> *reestim[2];
  Convolution *convol;
  ConvolutionData *convol_histo;


  convol = NULL;
  error.init();

  if (known_dist.offset + unknown_dist.offset > offset) {
    status = false;
    error.update(STAT_error[STATR_KNOWN_DISTRIBUTION_MIN_VALUE]);
  }
  if (known_dist.mean + unknown_dist.offset >= mean) {
    status = false;
    error.update(STAT_error[STATR_KNOWN_DISTRIBUTION_MEAN]);
  }

  if ((nb_iter != I_DEFAULT) && (nb_iter < 1)) {
    status = false;
    error.update(STAT_error[STATR_NB_ITERATION]);
  }

  if ((weight != D_DEFAULT) && (weight <= 0.)) {
    status = false;
    error.update(STAT_error[STATR_PENALTY_WEIGHT]);
  }

  if (status) {

    // creation d'un objet Convolution

    convol = new Convolution(known_dist , unknown_dist);
    convol->convolution_data = new ConvolutionData(*this , convol->nb_distribution);
    convol_histo = convol->convolution_data;

    if (estimator == PENALIZED_LIKELIHOOD) {
      penalty = new double[convol->distribution[1]->alloc_nb_value];

      if (weight == D_DEFAULT) {
        if (penalty_type != ENTROPY) {
          weight = CONVOLUTION_DIFFERENCE_WEIGHT;
        }
        else {
          weight = CONVOLUTION_ENTROPY_WEIGHT;
        }
      }

      weight *= nb_element;
    }

    // initialisation des produits de convolution partiels et
    // creation des quantites de reestimation

    partial_convol[0] = convol->distribution[1];
    partial_convol[1] = convol->distribution[0];

    reestim[0] = NULL;
    reestim[1] = new Reestimation<double>(convol->distribution[1]->alloc_nb_value);

    convol->computation(nb_value);

    convol->distribution[1]->init(CATEGORICAL , I_DEFAULT , I_DEFAULT , D_DEFAULT , D_DEFAULT);

#   ifdef DEBUG
    os << " (" << convol->mean << " " << convol->variance << ")" << endl;
#   endif

    likelihood = D_INF;
    i = 0;

    do {
      i++;

      // calcul des quantites de reestimation

      convol->expectation_step(*this , partial_convol , reestim);

      if (estimator != PENALIZED_LIKELIHOOD) {
        reestim[1]->distribution_estimation(convol->distribution[1]);
      }
      else {
        reestim[1]->penalized_likelihood_estimation(convol->distribution[1] , weight ,
                                                    penalty_type , penalty , outside);
      }

      // calcul du produit de convolution estime et de la log-vraisemblance correspondante

      convol->computation(nb_value);
      previous_likelihood = likelihood;
      likelihood = convol->likelihood_computation(*this);

#     ifdef MESSAGE
      if ((i < 10) || ((i < 100) && (i % 10 == 0)) || ((i < 1000) && (i % 100 == 0)) || (i % 1000 == 0)) {
        os << STAT_label[STATL_ITERATION] << " " << i << "   "
           << STAT_label[STATL_LIKELIHOOD] << ": " << likelihood << "   "
           << STAT_label[STATL_SMOOTHNESS] << ": " << convol->distribution[1]->second_difference_norm_computation();
        if (estimator == PENALIZED_LIKELIHOOD) {
          os << "   cumul: " << convol->distribution[1]->cumul[convol->distribution[1]->nb_value - 1];
        }
        os << endl;
      }
#     endif

    }
    while ((likelihood != D_INF) && (((nb_iter == I_DEFAULT) && (i < CONVOLUTION_NB_ITER) && 
             ((likelihood - previous_likelihood) / -likelihood > CONVOLUTION_LIKELIHOOD_DIFF)) ||
            ((nb_iter != I_DEFAULT) && (i < nb_iter))));

    if (likelihood != D_INF) {

#     ifdef MESSAGE
      os << "\n" << i << " " << STAT_label[STATL_ITERATIONS] << "   "
         << STAT_label[STATL_LIKELIHOOD] << ": " << likelihood << "   "
         << STAT_label[STATL_SMOOTHNESS] << ": " << convol->distribution[1]->second_difference_norm_computation();
      if (estimator == PENALIZED_LIKELIHOOD) {
        os << "   cumul: " << convol->distribution[1]->cumul[convol->distribution[1]->nb_value - 1];
      }
      os << endl;
#     endif

      if (estimator == PARAMETRIC_REGULARIZATION) {
        dist_flag[0] = false;
        dist_flag[1] = true;

        likelihood = D_INF;
        nb_likelihood_decrease = 0;

        i = 0;
        do {
          i++;

          // calcul des quantites de reestimation

          convol->expectation_step(*this , partial_convol , reestim);
          convol_histo->frequency_distribution[1]->update(reestim[1] ,
                                                          (int)(reestim[1]->nb_element *
                                                          MAX(sqrt(reestim[1]->variance) , 1.) * CONVOLUTION_COEFF));
          hlikelihood = convol_histo->frequency_distribution[1]->Reestimation<int>::type_parametric_estimation(convol->distribution[1] ,
                                                                                                               MIN(unknown_dist.offset , 1) , true ,
                                                                                                               CONVOLUTION_THRESHOLD);

          if (hlikelihood == D_INF) {
            likelihood = D_INF;
          }

          // calcul du produit de convolution estime et de la log-vraisemblance correspondante

          else {
            convol->computation(nb_value , CONVOLUTION_THRESHOLD , dist_flag);
            previous_likelihood = likelihood;
            likelihood = convol->likelihood_computation(*this);

            if (likelihood < previous_likelihood) {
              nb_likelihood_decrease++;
            }
            else {
              nb_likelihood_decrease = 0;
            }

#           ifdef DEBUG
            if ((i < 10) || (i % 10 == 0)) {
              os << STAT_label[STATL_ITERATION] << " " << i << "   "
                 << STAT_label[STATL_LIKELIHOOD] << ": " << likelihood << "   "
                 << STAT_label[STATL_SMOOTHNESS] << ": " << convol->distribution[1]->second_difference_norm_computation() << endl;
            }
#           endif

          }
        }
        while ((likelihood != D_INF) && (i < CONVOLUTION_NB_ITER) &&
               (((likelihood - previous_likelihood) / -likelihood > CONVOLUTION_LIKELIHOOD_DIFF) ||
                (hlikelihood == D_INF) || (nb_likelihood_decrease == 1)));

#       ifdef MESSAGE
        if (likelihood != D_INF) {
          os << "\n" << i << " " << STAT_label[STATL_ITERATIONS] << "   "
             << STAT_label[STATL_LIKELIHOOD] << ": " << likelihood << "   "
             << STAT_label[STATL_SMOOTHNESS] << ": " << convol->distribution[1]->second_difference_norm_computation() << endl;
        }
#       endif

      }
    }

    if (likelihood != D_INF) {

      // mise a jour du nombre de parametres inconnus

      convol->distribution[1]->nb_parameter_update();
      convol->nb_parameter = convol->distribution[1]->nb_parameter;

      reestim[0] = new Reestimation<double>(convol->distribution[0]->nb_value);

      convol->expectation_step(*this , partial_convol , reestim);
      convol_histo->frequency_distribution[0]->update(reestim[0] , nb_element);
      convol_histo->frequency_distribution[1]->update(reestim[1] , nb_element);

      delete reestim[0];
    }

    else {
      delete convol;
      convol = NULL;
      error.update(STAT_error[STATR_ESTIMATION_FAILURE]);
    }

    if (estimator == PENALIZED_LIKELIHOOD) {
      delete [] penalty;
    }

    delete reestim[1];
  }

  return convol;
}


/*--------------------------------------------------------------*
 *
 *  Estimation des parametres d'un produit de convolution de lois
 *  par l'algorithme EM.
 *
 *  arguments : reference sur un objet StatError, stream, reference sur la loi connue,
 *              borne inferieure minimum, type d'estimateur (vraisemblance,
 *              vraisemblance penalisee ou estimation d'une loi parametrique),
 *              nombre d'iterations, poids de la penalisation, type de penalisation,
 *              type de gestion des effets de bord (zero a l'exterieur du support ou
 *              prolongation de la loi).
 *
 *--------------------------------------------------------------*/

Convolution* FrequencyDistribution::convolution_estimation(StatError &error , ostream &os ,
                                                           const DiscreteParametric &known_dist ,
                                                           int min_inf_bound , int estimator ,
                                                           int nb_iter , double weight ,
                                                           int penalty_type , int outside) const

{
  bool status = true;
  double proba;
  DiscreteParametric *unknown_dist;
  Convolution *convol;


  convol = NULL;
  error.init();

  if (((min_inf_bound != 0) && (min_inf_bound != 1)) ||
      (min_inf_bound + known_dist.offset > offset)) {
    status = false;
    error.update(STAT_error[STATR_MIN_INF_BOUND]);
  }
  if (known_dist.offset + min_inf_bound > offset) {
    status = false;
    error.update(STAT_error[STATR_KNOWN_DISTRIBUTION_MIN_VALUE]);
  }
  if (known_dist.mean + min_inf_bound >= mean) {
    status = false;
    error.update(STAT_error[STATR_KNOWN_DISTRIBUTION_MEAN]);
  }

  if (status) {
    proba = 1. / (mean - known_dist.mean - min_inf_bound + 1.);
    if (proba > 1. - CONVOLUTION_INIT_PROBABILITY) {
      proba = 1. - CONVOLUTION_INIT_PROBABILITY;
    }
    else if (proba < CONVOLUTION_INIT_PROBABILITY) {
      proba = CONVOLUTION_INIT_PROBABILITY;
    }

    unknown_dist = new DiscreteParametric(NEGATIVE_BINOMIAL , min_inf_bound , I_DEFAULT ,
                                          1. , proba , CONVOLUTION_THRESHOLD);

#   ifdef DEBUG
    unknown_dist->ascii_print(os);
#   endif

    convol = convolution_estimation(error , os , known_dist , *unknown_dist , estimator ,
                                    nb_iter , weight , penalty_type , outside);

    delete unknown_dist;
  }

  return convol;
}


/*--------------------------------------------------------------*
 *
 *  Simulation par un produit de convolution de lois.
 *
 *  arguments : reference sur un objet StatError, effectif.
 *
 *--------------------------------------------------------------*/

ConvolutionData* Convolution::simulation(StatError &error , int nb_element) const

{
  register int i , j;
  int value , sum;
  ConvolutionData *convol_histo;


  error.init();

  if ((nb_element < 1) || (nb_element > DIST_NB_ELEMENT)) {
    convol_histo = NULL;
    error.update(STAT_error[STATR_SAMPLE_SIZE]);
  }

  else {

    // creation d'un objet ConvolutionData

    convol_histo = new ConvolutionData(*this);
    convol_histo->convolution = new Convolution(*this , false);

    for (i = 0;i < nb_element;i++) {
      sum = 0;

      for (j = 0;j < nb_distribution;j++) {

        // loi elementaire

        value = distribution[j]->simulation();
        sum += value;

        (convol_histo->frequency_distribution[j]->frequency[value])++;
      }

      // loi resultante

      (convol_histo->frequency[sum])++;
    }

    // extraction des caracteristiques des lois empiriques

    convol_histo->nb_value_computation();
    convol_histo->offset_computation();
    convol_histo->nb_element = nb_element;
    convol_histo->max_computation();
    convol_histo->mean_computation();
    convol_histo->variance_computation();

    for (i = 0;i < convol_histo->nb_distribution;i++) {
      convol_histo->frequency_distribution[i]->nb_value_computation();
      convol_histo->frequency_distribution[i]->offset_computation();
      convol_histo->frequency_distribution[i]->nb_element = nb_element;
      convol_histo->frequency_distribution[i]->max_computation();
      convol_histo->frequency_distribution[i]->mean_computation();
      convol_histo->frequency_distribution[i]->variance_computation();
    }
  }

  return convol_histo;
}


};  // namespace stat_tool
