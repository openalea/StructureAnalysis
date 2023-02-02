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
 *       $Id$
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

#include "convolution.h"
#include "stat_label.h"

using namespace std;


namespace stat_tool {



/*--------------------------------------------------------------*/
/**
 *  \brief Computation of a convolution of discrete distributions.
 *
 *  \param[in] min_nb_value    upper bound of the elementary distribution supports,
 *  \param[in] cumul_threshold threshold on the cumulative distribution function,
 *  \param[in] dist_flag       flags for computing elementary distributions.
 */
/*--------------------------------------------------------------*/

void Convolution::computation(int min_nb_value , double cumul_threshold , bool *dist_flag)

{
  int i;


  // computation of the elementary distributions

  if (dist_flag) {
    for (i = 0;i < nb_distribution;i++) {
      if (dist_flag[i]) {
        distribution[i]->computation(min_nb_value , cumul_threshold);
      }
    }
  }

  // computation of the resulting convolution

  convolution(*distribution[0] , *distribution[1]);
  for (i = 2;i < nb_distribution;i++) {
    convolution(*this , *distribution[i]);
  }

  cumul_computation();

  max_computation();
  mean_computation();
  variance_computation();
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of reestimation quantities (E-step of EM algorithm).
 *
 *  \param[in] histo          reference on a FrequencyDistribution object,
 *  \param[in] partial_convol pointer on the partial convolutions, 
 *  \param[in] reestim        pointer on the reestimation quantities.
 */
/*--------------------------------------------------------------*/

void Convolution::expectation_step(const FrequencyDistribution &histo ,
                                   const Distribution **partial_convol ,
                                   Reestimation<double> **reestim) const

{
  int i , j , k;
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
      cout << "\nreestimation quantities distribution " << i << " :" << *reestim[i] << endl;
#     endif

    }
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Deconvolution of an elementary distribution using the EM algorithm.
 *
 *  \param[in] error        reference on a StatError object,
 *  \param[in] display      flag for displaying estimation intermediate results,
 *  \param[in] known_dist   reference on the known distribution,
 *  \param[in] unknown_dist reference on the unknown distribution,
 *  \param[in] estimator    estimator type (likelihood, penalized likelihood or 
 *                          estimation of a parametric distribution),
 *  \param[in] nb_iter      number of iterations,
 *  \param[in] weight       penalty weight,
 *  \param[in] pen_type     penalty type,
 *  \param[in] outside      management of side effects (zero outside the support or
 *                          continuation of the distribution).
 *
 *  \return                 Convolution object.
 */
/*--------------------------------------------------------------*/

Convolution* FrequencyDistribution::convolution_estimation(StatError &error , bool display ,
                                                           const DiscreteParametric &known_dist ,
                                                           const DiscreteParametric &unknown_dist ,
                                                           estimation_criterion estimator , int nb_iter ,
                                                           double weight , penalty_type pen_type ,
                                                           side_effect outside) const

{
  bool status = true , dist_flag[2];
  int i;
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

    // construction of a Convolution object

    convol = new Convolution(known_dist , unknown_dist);
    convol->convolution_data = new ConvolutionData(*this , convol->nb_distribution);
    convol_histo = convol->convolution_data;

    if (estimator == PENALIZED_LIKELIHOOD) {
      penalty = new double[convol->distribution[1]->alloc_nb_value];

      if (weight == D_DEFAULT) {
        if (pen_type != ENTROPY) {
          weight = CONVOLUTION_DIFFERENCE_WEIGHT;
        }
        else {
          weight = CONVOLUTION_ENTROPY_WEIGHT;
        }
      }

      weight *= nb_element;
    }

    // initialization of partial convolutions and construction of the reestimation quantities

    partial_convol[0] = convol->distribution[1];
    partial_convol[1] = convol->distribution[0];

    reestim[0] = NULL;
    reestim[1] = new Reestimation<double>(convol->distribution[1]->alloc_nb_value);

    convol->computation(nb_value);

    convol->distribution[1]->init(CATEGORICAL , I_DEFAULT , I_DEFAULT , D_DEFAULT , D_DEFAULT);

#   ifdef DEBUG
    cout << " (" << convol->mean << " " << convol->variance << ")" << endl;
#   endif

    likelihood = D_INF;
    i = 0;

    do {
      i++;

      // computation of the reestimation quantities

      convol->expectation_step(*this , partial_convol , reestim);

      if (estimator != PENALIZED_LIKELIHOOD) {
        reestim[1]->distribution_estimation(convol->distribution[1]);
      }
      else {
        reestim[1]->penalized_likelihood_estimation(convol->distribution[1] , weight ,
                                                    pen_type , penalty , outside);
      }

      // computation of the estimated convolution and the corresponding log-likelihood

      convol->computation(nb_value);
      previous_likelihood = likelihood;
      likelihood = convol->likelihood_computation(*this);

      // display of estimation results

      if ((display) && ((i < 10) || ((i < 100) && (i % 10 == 0)) || ((i < 1000) && (i % 100 == 0)) || (i % 1000 == 0))) {
        cout << STAT_label[STATL_ITERATION] << " " << i << "   "
             << STAT_label[STATL_LIKELIHOOD] << ": " << likelihood << "   "
             << STAT_label[STATL_SMOOTHNESS] << ": " << convol->distribution[1]->second_difference_norm_computation();
        if (estimator == PENALIZED_LIKELIHOOD) {
          cout << "   cumul: " << convol->distribution[1]->cumul[convol->distribution[1]->nb_value - 1];
        }
        cout << endl;
      }
    }
    while ((likelihood != D_INF) && (((nb_iter == I_DEFAULT) && (i < CONVOLUTION_NB_ITER) && 
             ((likelihood - previous_likelihood) / -likelihood > CONVOLUTION_LIKELIHOOD_DIFF)) ||
            ((nb_iter != I_DEFAULT) && (i < nb_iter))));

    if (likelihood != D_INF) {

      // display of estimation results

      if (display) {
        cout << "\n" << i << " " << STAT_label[STATL_ITERATIONS] << "   "
             << STAT_label[STATL_LIKELIHOOD] << ": " << likelihood << "   "
             << STAT_label[STATL_SMOOTHNESS] << ": " << convol->distribution[1]->second_difference_norm_computation();
        if (estimator == PENALIZED_LIKELIHOOD) {
          cout << "   cumul: " << convol->distribution[1]->cumul[convol->distribution[1]->nb_value - 1];
        }
        cout << endl;
      }

      if (estimator == PARAMETRIC_REGULARIZATION) {
        dist_flag[0] = false;
        dist_flag[1] = true;

        likelihood = D_INF;
        nb_likelihood_decrease = 0;

        i = 0;
        do {
          i++;

          // computation of the reestimation quantities

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

          // computation of the estimated convolution and the corresponding log-likelihood

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
            if ((display) && ((i < 10) || (i % 10 == 0))) {
              cout << STAT_label[STATL_ITERATION] << " " << i << "   "
                   << STAT_label[STATL_LIKELIHOOD] << ": " << likelihood << "   "
                   << STAT_label[STATL_SMOOTHNESS] << ": " << convol->distribution[1]->second_difference_norm_computation() << endl;
            }
#           endif

          }
        }
        while ((likelihood != D_INF) && (i < CONVOLUTION_NB_ITER) &&
               (((likelihood - previous_likelihood) / -likelihood > CONVOLUTION_LIKELIHOOD_DIFF) ||
                (hlikelihood == D_INF) || (nb_likelihood_decrease == 1)));

        if ((display) && (likelihood != D_INF)) {
          cout << "\n" << i << " " << STAT_label[STATL_ITERATIONS] << "   "
               << STAT_label[STATL_LIKELIHOOD] << ": " << likelihood << "   "
               << STAT_label[STATL_SMOOTHNESS] << ": " << convol->distribution[1]->second_difference_norm_computation() << endl;
        }
      }
    }

    if (likelihood != D_INF) {

      // update of the number of free parameters

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


/*--------------------------------------------------------------*/
/**
 *  \brief Deconvolution of an elementary distribution using the EM algorithm.
 *
 *  \param[in] error         reference on a StatError object,
 *  \param[in] display       flag for displaying estimation intermediate results,
 *  \param[in] known_dist    reference on the known distribution,
 *  \param[in] min_inf_bound minimum lower bound of the support,
 *  \param[in] estimator     estimator type (likelihood, penalized likelihood or 
 *                           estimation of a parametric distribution),
 *  \param[in] nb_iter       number of iterations,
 *  \param[in] weight        penalty weight,
 *  \param[in] pen_type      penalty type,
 *  \param[in] outside       management of side effects (zero outside the support or
 *                           continuation of the distribution).
 *
 *  \return                  Convolution object.
 */
/*--------------------------------------------------------------*/

Convolution* FrequencyDistribution::convolution_estimation(StatError &error , bool display ,
                                                           const DiscreteParametric &known_dist ,
                                                           int min_inf_bound , estimation_criterion estimator ,
                                                           int nb_iter , double weight ,
                                                           penalty_type pen_type , side_effect outside) const

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

    convol = convolution_estimation(error , display , known_dist , *unknown_dist , estimator ,
                                    nb_iter , weight , pen_type , outside);

    delete unknown_dist;
  }

  return convol;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Simulation using a convolution of discrete distributions.
 *
 *  \param[in] error       reference on a StatError object,
 *  \param[in] nb_element  sample size.
 *
 *  \return                ConvolutionData object.
 */
/*--------------------------------------------------------------*/

ConvolutionData* Convolution::simulation(StatError &error , int nb_element) const

{
  int i , j;
  int value , sum;
  ConvolutionData *convol_histo;


  error.init();

  if ((nb_element < 1) || (nb_element > DIST_NB_ELEMENT)) {
    convol_histo = NULL;
    error.update(STAT_error[STATR_SAMPLE_SIZE]);
  }

  else {

    // construction of a ConvolutionData object

    convol_histo = new ConvolutionData(*this);
    convol_histo->convolution = new Convolution(*this , false);

    for (i = 0;i < nb_element;i++) {
      sum = 0;

      for (j = 0;j < nb_distribution;j++) {

        // elementary distribution

        value = distribution[j]->simulation();
        sum += value;

        (convol_histo->frequency_distribution[j]->frequency[value])++;
      }

      // resulting convolution

      (convol_histo->frequency[sum])++;
    }

    // computation of frequency distribution characteristics

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
