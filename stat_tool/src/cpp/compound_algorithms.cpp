/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       StructureAnalysis: Identifying patterns in plant architecture and development
 *
 *       Copyright 1995-2019 CIRAD AGAP
 *
 *       File author(s): Yann Guedon (yann.guedon@cirad.fr)
 *
 *       $Source$
 *       $Id$
 *
 *       Forum for StructureAnalysis developers:
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



#include <cmath>

#include "compound.h"
#include "stat_label.h"

using namespace std;


namespace stat_tool {



/*--------------------------------------------------------------*/
/**
 *  \brief Computation of a compound distribution.
 *
 *  \param[in] min_nb_value    lower bound of the support,
 *  \param[in] cumul_threshold threshold on the cumulative distribution function,
 *  \param[in] sum_flag        flag for the computation of the sum distribution,
 *  \param[in] dist_flag       flag for the computation of the basis distribution.
 */
/*--------------------------------------------------------------*/

void Compound::computation(int min_nb_value , double cumul_threshold ,
                           bool sum_flag , bool dist_flag)

{
  int i , j;
  DiscreteParametric *power_dist;


  // computation of the sum distribution and the basis distribution

  if (sum_flag) {
    sum_distribution->computation(1 , CUMUL_THRESHOLD);
  }
  if (dist_flag) {
    distribution->computation(min_nb_value , cumul_threshold);
  }

  // computation of convolution powers of the basis distribution and
  // computation of the resulting compound distribution

  power_dist = new DiscreteParametric((sum_distribution->nb_value - 1) * (distribution->nb_value - 1) + 1 ,
                                      distribution->ident);

  for (i = 0;i < MIN(power_dist->nb_value , alloc_nb_value);i++) {
    mass[i] = 0.;
  }
  if (sum_distribution->offset == 0) {
    mass[0] += sum_distribution->mass[0];
  }

  if (((distribution->ident == CATEGORICAL) || (distribution->ident == UNIFORM)) &&
      (sum_distribution->offset > 1)) {
    power_dist->mass_copy(*distribution);
    for (i = 2;i < sum_distribution->offset;i++) {
      power_dist->convolution(*power_dist , *distribution);
    }
  }

  for (i = MAX(sum_distribution->offset , 1);i < sum_distribution->nb_value;i++) {
    if (i == 1) {
      power_dist->mass_copy(*distribution);
    }

    else {
      if ((distribution->ident == CATEGORICAL) || (distribution->ident == UNIFORM)) {
        power_dist->convolution(*power_dist , *distribution);
      }
 
      else {
        power_dist->inf_bound = i * distribution->inf_bound;
        power_dist->sup_bound = (distribution->sup_bound != I_DEFAULT ? i * distribution->sup_bound : I_DEFAULT);
        power_dist->parameter = (distribution->parameter != D_DEFAULT ? i * distribution->parameter : D_DEFAULT);
        power_dist->probability = distribution->probability;

        power_dist->computation(min_nb_value , cumul_threshold);
      }
    }

    for (j = power_dist->offset;j < MIN(power_dist->nb_value , alloc_nb_value);j++) {
      mass[j] += sum_distribution->mass[i] * power_dist->mass[j];
    }
  }

  offset_computation();
  nb_value = MIN(power_dist->nb_value , alloc_nb_value);

  cumul_computation();

  max_computation();
  mean_computation();
  variance_computation();

  delete power_dist;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of a compound distribution.
 *
 *  \param[in] power_dist      pointer on the convolution powers of the the basis distribution,
 *  \param[in] min_nb_value    lower bound of the support,
 *  \param[in] cumul_threshold threshold on the cumulative distribution function,
 *  \param[in] sum_flag        flag for the computation of the sum distribution,
 *  \param[in] dist_flag       flag for the computation of the basis distribution.
 */
/*--------------------------------------------------------------*/

void Compound::computation(DiscreteParametric **power_dist , int min_nb_value ,
                           double cumul_threshold , bool sum_flag , bool dist_flag)

{
  int i , j;
  int sum_nb_value , min;


  // computation of the sum distribution and the basis distribution

  sum_nb_value = sum_distribution->nb_value;
  if (sum_flag) {
    sum_distribution->computation(1 , CUMUL_THRESHOLD);
  }
  if (dist_flag) {
    distribution->computation(min_nb_value , cumul_threshold);
  }

  // computation of convolution powers of the basis distribution

  if (dist_flag) {
    min = MAX(sum_distribution->offset - 1 , 1);
  }
  else {
    min = sum_nb_value;
  }

  if ((distribution->ident == CATEGORICAL) || (distribution->ident == UNIFORM)) {
    power_dist[min]->ident = distribution->ident;

    if (dist_flag) {
      power_dist[min]->mass_copy(*distribution);
      for (i = 2;i <= min;i++) {
        power_dist[min]->convolution(*power_dist[min] , *distribution);
      }
    }

    else if (min < sum_distribution->nb_value) {
      power_dist[min]->convolution(*power_dist[min - 1] , *distribution);
    }

    for (i = min + 1;i < sum_distribution->nb_value;i++) {
      power_dist[i]->ident = distribution->ident;
      power_dist[i]->convolution(*power_dist[i - 1] , *distribution);
    }
  }

  else {
    if (min == 1) {
      power_dist[min]->mass_copy(*distribution);
      min++;
    }

    for (i = min;i < sum_distribution->nb_value;i++) {
      power_dist[i]->ident = distribution->ident;
      power_dist[i]->inf_bound = i * distribution->inf_bound;
      power_dist[i]->sup_bound = (distribution->sup_bound != I_DEFAULT ? i * distribution->sup_bound : I_DEFAULT);
      power_dist[i]->parameter = (distribution->parameter != D_DEFAULT ? i * distribution->parameter : D_DEFAULT);
      power_dist[i]->probability = distribution->probability;

      power_dist[i]->computation(min_nb_value , cumul_threshold);
    }
  }

  // computation of the resulting compound distribution

  if (sum_distribution->offset == 0) {
    offset = 0;
  }
  else {
    offset = power_dist[sum_distribution->offset]->offset;
  }
  nb_value = MIN(power_dist[sum_distribution->nb_value - 1]->nb_value , alloc_nb_value);

  for (i = 0;i < nb_value;i++) {
    mass[i] = 0.;
  }
  if (sum_distribution->offset == 0) {
    mass[0] += sum_distribution->mass[0];
  }

  for (i = MAX(sum_distribution->offset , 1);i < sum_distribution->nb_value;i++) {
    for (j = power_dist[i]->offset;j < power_dist[i]->nb_value;j++) {
      mass[j] += sum_distribution->mass[i] * power_dist[i]->mass[j];
    }
  }

  cumul_computation();

  max_computation();
  mean_computation();
  variance_computation();
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of reestimation quantities (E-step of the EM algorithm).
 *
 *  \param[in] histo       reference on a FrequencyDistribution object,
 *  \param[in] power_dist  pointer on the convolution powers of the basis distribution 
 *  \param[in] sum_reestim pointer on the reestimation quantities of the sum distribution
 *  \param[in] reestim     pointer on the reestimation quantities of the basis distribution.
 */
/*--------------------------------------------------------------*/

void Compound::expectation_step(const FrequencyDistribution &histo ,
                                DiscreteParametric **power_dist ,
                                Reestimation<double> *sum_reestim ,
                                Reestimation<double> *reestim) const

{
  int i , j , k;
  double num , denom , *term;


  term = new double[sum_distribution->nb_value];

  // initializations

  if (sum_reestim) {
    for (i = 0;i < sum_reestim->alloc_nb_value;i++) {
      sum_reestim->frequency[i] = 0.;
    }
  }

  if (reestim) {
    for (i = 0;i < reestim->alloc_nb_value;i++) {
      reestim->frequency[i] = 0.;
    }
  }

  for (i = histo.offset;i < histo.nb_value;i++) {
    if (histo.frequency[i] > 0) {

      // computation of the normalizing term

      denom = 0.;
      if (sum_distribution->offset == 0) {
        if (i == 0) {
          term[0] = sum_distribution->mass[0];
          denom += term[0];
        }
        else {
          term[0] = 0.;
        }
      }

      for (j = MAX(sum_distribution->offset , 1);j < sum_distribution->nb_value;j++) {
        if ((i >= power_dist[j]->offset) && (i < power_dist[j]->nb_value)) {
          term[j] = sum_distribution->mass[j] * power_dist[j]->mass[i];
          denom += term[j];
        }
        else {
          term[j] = 0.;
        }
      }

      if (denom > 0.) {

        // accumulation of reestimation quantities for the sum distribution

        if (sum_reestim) {
          for (j = sum_distribution->offset;j < sum_distribution->nb_value;j++) {
            sum_reestim->frequency[j] += histo.frequency[i] * term[j] / denom;
          }
        }

        // accumulation of reestimation quantities for the basis distribution

        if (reestim) {
          for (j = distribution->offset;j <= MIN(i , distribution->nb_value - 1);j++) {
            num = 0.;

            for (k = sum_distribution->nb_value - 1;k >= MAX(sum_distribution->offset , 1);k--) {
              if (k == 1) {
                if (i == j) {
                  num += sum_distribution->mass[k];
                }
              }

              else {
                if (i - j < power_dist[k - 1]->nb_value) {
                  if (i - j >= power_dist[k - 1]->offset) {
                    num += sum_distribution->mass[k] * k * power_dist[k - 1]->mass[i - j];
                  }
                }
                else {
                  break;
                }
              }
            }

            reestim->frequency[j] += histo.frequency[i] * distribution->mass[j] * num / denom;
          }
        }
      }
    }
  }

  delete [] term;

  if (sum_reestim) {
    sum_reestim->nb_value_computation();
    sum_reestim->offset_computation();
    sum_reestim->nb_element_computation();
    sum_reestim->max_computation();
    sum_reestim->mean_computation();
    sum_reestim->variance_computation();
  }

  if (reestim) {
    reestim->nb_value_computation();
    reestim->offset_computation();
    reestim->nb_element_computation();
    reestim->max_computation();
    reestim->mean_computation();
    reestim->variance_computation();
  }

# ifdef DEBUG
  if (sum_reestim) {
    cout << "\nsum distribution reestimation quantities:" << *sum_reestim << endl;
  }
  if (reestim) {
    cout << "\nbasis distribution reestimation quantities:" << *reestim << endl;
  }
# endif

}


/*--------------------------------------------------------------*/
/**
 *  \brief Estimation of a compound distribution using the EM algorithm.
 *
 *  \param[in] error     reference on a StatError object,
 *  \param[in] os        stream for displaying estimation intermediate results,
 *  \param[in] sum_dist  reference on the sum distribution,
 *  \param[in] dist      reference on the basis distribution,
 *  \param[in] type      unknown distribution type (SUM/ELEMENTARY),
 *  \param[in] estimator estimator type (likelihood, penalized likelihood or
 *                       estimation of a parametric distribution),
 *  \param[in] nb_iter   number of iterations,
 *  \param[in] weight    penalty weight,
 *  \param[in] pen_type  penalty type,
 *  \param[in] outside   management of side effects (zero outside the support or
 *                       continuation of the distribution).
 *
 *  \return              Compound object.
 */
/*--------------------------------------------------------------*/

Compound* FrequencyDistribution::compound_estimation(StatError &error , ostream *os ,
                                                     const DiscreteParametric &sum_dist ,
                                                     const DiscreteParametric &dist ,
                                                     compound_distribution type ,
                                                     estimation_criterion estimator ,
                                                     int nb_iter , double weight ,
                                                     penalty_type pen_type , side_effect outside) const

{
  bool status = true , sum_compute , dist_compute;
  int i;
  int sum_nb_value , nb_likelihood_decrease;
  double likelihood , previous_likelihood , hlikelihood , *penalty;
  DiscreteParametric **power_dist;
  Reestimation<double> *sum_reestim , *reestim;
  Compound *compound;
  CompoundData *compound_histo;


  compound = NULL;
  error.init();

  if ((nb_iter != I_DEFAULT) && (nb_iter < 1)) {
    status = false;
    error.update(STAT_error[STATR_NB_ITERATION]);
  }

  if ((weight != D_DEFAULT) && (weight <= 0.)) {
    status = false;
    error.update(STAT_error[STATR_PENALTY_WEIGHT]);
  }

  if (status) {

    // construction of a Compound object

    compound = new Compound(sum_dist , dist , type);
    compound->compound_data = new CompoundData(*this , *compound);
    compound_histo = compound->compound_data;

    if (estimator == PENALIZED_LIKELIHOOD) {
      switch (type) {
      case SUM :
        penalty = new double[compound->sum_distribution->alloc_nb_value];
        break;
      case ELEMENTARY :
        penalty = new double[compound->distribution->alloc_nb_value];
        break;
      }

      if (weight == D_DEFAULT) {
        if (pen_type != ENTROPY) {
          weight = COMPOUND_DIFFERENCE_WEIGHT;
        }
        else {
          weight = COMPOUND_ENTROPY_WEIGHT;
        }
      }

      switch (type) {
      case SUM :
        weight *= nb_element;
        break;
      case ELEMENTARY :
        weight *= compound->sum_distribution->mean * nb_element;
        break;
      }
    }

    switch (type) {
    case SUM :
      sum_compute = true;
      dist_compute = false;
      break;
    case ELEMENTARY :
      sum_compute = false;
      dist_compute = true;
      break;
    }

    // construction of the convolution powers of the basis distribution and
    // the reestimation quantities

    sum_nb_value = compound->sum_distribution->alloc_nb_value;

    power_dist = new DiscreteParametric*[sum_nb_value];
    for (i = MAX(sum_dist.offset - 1 , 1);i < sum_nb_value;i++) {
      power_dist[i] = new DiscreteParametric(i * (compound->distribution->alloc_nb_value - 1) + 1 ,
                                             compound->distribution->ident);
    }

    sum_reestim = new Reestimation<double>(sum_nb_value);
    reestim = new Reestimation<double>(compound->distribution->alloc_nb_value);

    compound->computation(power_dist , nb_value , COMPOUND_THRESHOLD ,
                          sum_compute , true);

#   ifdef DEBUG
    cout << " (" << compound->mean << " " << compound->variance << ")" << endl;
#   endif

    switch (type) {
    case SUM :
      compound->sum_distribution->init(CATEGORICAL , I_DEFAULT , I_DEFAULT , D_DEFAULT , D_DEFAULT);
      break;
    case ELEMENTARY :
      compound->distribution->init(CATEGORICAL , I_DEFAULT , I_DEFAULT , D_DEFAULT , D_DEFAULT);
      break;
    }

    likelihood = D_INF;
    i = 0;

    do {
      i++;

      // computation of the reestimation quantities

      switch (type) {

      case SUM : {
        compound->expectation_step(*this , power_dist , sum_reestim , NULL);

        if (estimator != PENALIZED_LIKELIHOOD) {
          sum_reestim->distribution_estimation(compound->sum_distribution);
        }
        else {
          sum_reestim->penalized_likelihood_estimation(compound->sum_distribution , weight ,
                                                       pen_type , penalty , outside);
        }
        break;
      }

      case ELEMENTARY : {
        compound->expectation_step(*this , power_dist , NULL , reestim);

        if (estimator != PENALIZED_LIKELIHOOD) {
          reestim->distribution_estimation(compound->distribution);
        }
        else {
          reestim->penalized_likelihood_estimation(compound->distribution , weight ,
                                                   pen_type , penalty , outside);
        }
        break;
      }
      }

      // computation of the estimated compound distribution and the corresponding log-likelihood

      compound->computation(power_dist , nb_value , COMPOUND_THRESHOLD ,
                            sum_compute , dist_compute);
      previous_likelihood = likelihood;
      likelihood = compound->likelihood_computation(*this);

      // display of estimation results

      if ((os) && ((i < 10) || ((i < 100) && (i % 10 == 0)) || ((i < 1000) && (i % 100 == 0)) || (i % 1000 == 0))) {
        *os << STAT_label[STATL_ITERATION] << " " << i << "   "
            << STAT_label[STATL_LIKELIHOOD] << ": " << likelihood << "   "
            << STAT_label[STATL_SMOOTHNESS] << ": ";

        switch (type) {

        case SUM : {
          *os << compound->sum_distribution->second_difference_norm_computation();
          if (estimator == PENALIZED_LIKELIHOOD) {
            *os << "   cumul: " << compound->sum_distribution->cumul[compound->sum_distribution->nb_value - 1];
          }
          break;
        }

        case ELEMENTARY : {
          *os << compound->distribution->second_difference_norm_computation();
          if (estimator == PENALIZED_LIKELIHOOD) {
            *os << "   cumul: " << compound->distribution->cumul[compound->distribution->nb_value - 1];
          }
          break;
        }
        }

        *os << endl;
      }
    }
    while ((likelihood != D_INF) && (((nb_iter == I_DEFAULT) && (i < COMPOUND_NB_ITER) && 
             ((likelihood - previous_likelihood) / -likelihood > COMPOUND_LIKELIHOOD_DIFF)) ||
            ((nb_iter != I_DEFAULT) && (i < nb_iter))));

    if (likelihood != D_INF) {
      if (os) {
        *os << "\n" << i << " " << STAT_label[STATL_ITERATIONS] << "   "
             << STAT_label[STATL_LIKELIHOOD] << ": " << likelihood << "   "
             << STAT_label[STATL_SMOOTHNESS] << ": ";

        switch (type) {

        case SUM : {
          *os << compound->sum_distribution->second_difference_norm_computation();
          if (estimator == PENALIZED_LIKELIHOOD) {
            *os << "   cumul: " << compound->sum_distribution->cumul[compound->sum_distribution->nb_value - 1];
          }
          break;
        }

        case ELEMENTARY : {
          *os << compound->distribution->second_difference_norm_computation();
          if (estimator == PENALIZED_LIKELIHOOD) {
            *os << "   cumul: " << compound->distribution->cumul[compound->distribution->nb_value - 1];
          }
          break;
        }
        }

        *os << endl;
      }

      if (estimator == PARAMETRIC_REGULARIZATION) {
        likelihood = D_INF;
        nb_likelihood_decrease = 0;

        i = 0;
        do {
          i++;

          // computation of the reestimation quantities

          switch (type) {

          case SUM : {
            compound->expectation_step(*this , power_dist , sum_reestim , NULL);
            compound_histo->sum_frequency_distribution->update(sum_reestim ,
                                                               (int)(sum_reestim->nb_element *
                                                               MAX(sqrt(sum_reestim->variance) , 1.) * COMPOUND_COEFF));
            hlikelihood = compound_histo->sum_frequency_distribution->Reestimation<int>::type_parametric_estimation(compound->sum_distribution ,
                                                                                                                    MIN(sum_dist.offset , 1) , true);
            break;
          }

          case ELEMENTARY : {
            compound->expectation_step(*this , power_dist , NULL , reestim);
            compound_histo->frequency_distribution->update(reestim ,
                                                           (int)(reestim->nb_element *
                                                           MAX(sqrt(reestim->variance) , 1.) * COMPOUND_COEFF));
            hlikelihood = compound_histo->frequency_distribution->Reestimation<int>::type_parametric_estimation(compound->distribution ,
                                                                                                                MIN(dist.offset , 1) , true ,
                                                                                                                COMPOUND_THRESHOLD);
            break;
          }
          }

          if (hlikelihood == D_INF) {
            likelihood = D_INF;
          }

          // computation of the estimated compound distribution and the corresponding log-likelihood

          else {
            compound->computation(power_dist , nb_value , COMPOUND_THRESHOLD ,
                                  sum_compute , dist_compute);
            previous_likelihood = likelihood;
            likelihood = compound->likelihood_computation(*this);

            if (likelihood < previous_likelihood) {
              nb_likelihood_decrease++;
            }
            else {
              nb_likelihood_decrease = 0;
            }

#           ifdef DEBUG
            if ((os) && ((i < 10) || (i % 10 == 0))) {
              *os << STAT_label[STATL_ITERATION] << " " << i << "   "
                  << STAT_label[STATL_LIKELIHOOD] << ": " << likelihood << "   "
                  << STAT_label[STATL_SMOOTHNESS] << ": ";

              switch (type) {
              case SUM :
                *os << compound->sum_distribution->second_difference_norm_computation() << endl;
                break;
              case ELEMENTARY :
                *os << compound->distribution->second_difference_norm_computation() << endl;
                break;
              }
            }
#           endif

          }
        }
        while ((likelihood != D_INF) && (i < COMPOUND_NB_ITER) &&
               (((likelihood - previous_likelihood) / -likelihood > COMPOUND_LIKELIHOOD_DIFF) ||
                (hlikelihood == D_INF) || (nb_likelihood_decrease == 1)));

        if ((os) && (likelihood != D_INF)) {
          *os << "\n" << i << " " << STAT_label[STATL_ITERATIONS] << "   "
              << STAT_label[STATL_LIKELIHOOD] << ": " << likelihood << "   "
              << STAT_label[STATL_SMOOTHNESS] << ": ";

          switch (type) {
          case SUM :
            *os << compound->sum_distribution->second_difference_norm_computation() << endl;
            break;
          case ELEMENTARY :
            *os << compound->distribution->second_difference_norm_computation() << endl;
            break;
          }
        }
      }
    }

    if (likelihood != D_INF) {

      // update of the number of free parameters

      switch (type) {
      case SUM :
        compound->sum_distribution->nb_parameter_update();
        compound->nb_parameter = compound->sum_distribution->nb_parameter;
        break;
      case ELEMENTARY :
        compound->distribution->nb_parameter_update();
        compound->nb_parameter = compound->distribution->nb_parameter;
        break;
      }

      compound->expectation_step(*this , power_dist , sum_reestim , reestim);
      compound_histo->sum_frequency_distribution->update(sum_reestim , nb_element);
      compound_histo->frequency_distribution->update(reestim , (int)(nb_element * compound->sum_distribution->mean));
    }

    else {
      delete compound;
      compound = NULL;
      error.update(STAT_error[STATR_ESTIMATION_FAILURE]);
    }

    if (estimator == PENALIZED_LIKELIHOOD) {
      delete [] penalty;
    }

    for (i = MAX(sum_dist.offset - 1 , 1);i < sum_nb_value;i++) {
      delete power_dist[i];
    }
    delete [] power_dist;

    delete sum_reestim;
    delete reestim;
  }

  return compound;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Estimation of a compound distribution using the EM algorithm.
 *
 *  \param[in] error         reference on a StatError object,
 *  \param[in] os            stream for displaying estimation intermediate results,
 *  \param[in] known_dist    reference on the known distribution,
 *  \param[in] type          unknown distribution type (SUM/ELEMENTARY),
 *  \param[in] min_inf_bound minimum lower bound of the support of the unknown distribution,
 *  \param[in] estimator     estimator type (likelihood, penalized likelihood or
 *                           estimation of a parametric distribution),
 *  \param[in] nb_iter       number of iterations,
 *  \param[in] weight        penalty weight ,
 *  \param[in] pen_type      penalty type,
 *  \param[in] outside       management of side effects (zero outside the support or
 *                           continuation of the distribution).
 *
 *  \return                  Compound object.
 */
/*--------------------------------------------------------------*/

Compound* FrequencyDistribution::compound_estimation(StatError &error , ostream *os ,
                                                     const DiscreteParametric &known_dist ,
                                                     compound_distribution type , int min_inf_bound ,
                                                     estimation_criterion estimator , int nb_iter ,
                                                     double weight , penalty_type pen_type ,
                                                     side_effect outside) const

{
  double proba;
  DiscreteParametric *unknown_dist;
  Compound *compound;


  error.init();

  if ((min_inf_bound != 0) && (min_inf_bound != 1)) {
    compound = NULL;
    error.update(STAT_error[STATR_MIN_INF_BOUND]);
  }

  else {
    proba = 1. / ((mean / known_dist.mean) - min_inf_bound + 1.);
    if (proba > 1. - COMPOUND_INIT_PROBABILITY) {
      proba = 1. - COMPOUND_INIT_PROBABILITY;
    }
    else if (proba < COMPOUND_INIT_PROBABILITY) {
      proba = COMPOUND_INIT_PROBABILITY;
    }

    switch (type) {
    case SUM :
      unknown_dist = new DiscreteParametric(NEGATIVE_BINOMIAL , min_inf_bound , I_DEFAULT ,
                                            1. , proba);
      break;
    case ELEMENTARY :
      unknown_dist = new DiscreteParametric(NEGATIVE_BINOMIAL , min_inf_bound , I_DEFAULT ,
                                            1. , proba , COMPOUND_THRESHOLD);
      break;
    }

#   ifdef DEBUG
    unknown_dist->ascii_print(cout);
#   endif

    switch (type) {
    case SUM :
      compound = compound_estimation(error , os , *unknown_dist , known_dist , type ,
                                     estimator , nb_iter , weight , pen_type , outside);
      break;
    case ELEMENTARY :
      compound = compound_estimation(error , os , known_dist , *unknown_dist , type ,
                                     estimator , nb_iter , weight , pen_type , outside);
      break;
    }

    delete unknown_dist;
  }

  return compound;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Simulation using a compound distribution.
 *
 *  \param[in] error      reference on a StatError object,
 *  \param[in] nb_element sample size.
 *
 *  \return               sample generated by a compound distribution.
 */
/*--------------------------------------------------------------*/

CompoundData* Compound::simulation(StatError &error , int nb_element) const

{
  int i , j;
  int nb_dist , sum , value;
  CompoundData *compound_histo;


  error.init();

  if ((nb_element < 1) || (nb_element > DIST_NB_ELEMENT)) {
    compound_histo = NULL;
    error.update(STAT_error[STATR_SAMPLE_SIZE]);
  }

  else {

    // construction of a CompoundData object

    compound_histo = new CompoundData(*this);
    compound_histo->compound = new Compound(*this , false);

    for (i = 0;i < nb_element;i++) {

      // sum distribution

      nb_dist = sum_distribution->simulation();
      (compound_histo->sum_frequency_distribution->frequency[nb_dist])++;

      sum = 0;
      for (j = 0;j < nb_dist;j++) {

        // basis distribution

        value = distribution->simulation();
        sum += value;
        (compound_histo->frequency_distribution->frequency[value])++;
      }

      // resulting compound distribution

      (compound_histo->frequency[sum])++;
    }

    // computation of frequency distribution characteristics

    compound_histo->nb_value_computation();
    compound_histo->offset_computation();
    compound_histo->nb_element = nb_element;
    compound_histo->max_computation();
    compound_histo->mean_computation();
    compound_histo->variance_computation();

    compound_histo->sum_frequency_distribution->nb_value_computation();
    compound_histo->sum_frequency_distribution->offset_computation();
    compound_histo->sum_frequency_distribution->nb_element_computation();
    compound_histo->sum_frequency_distribution->max_computation();
    compound_histo->sum_frequency_distribution->mean_computation();
    compound_histo->sum_frequency_distribution->variance_computation();

    compound_histo->frequency_distribution->nb_value_computation();
    compound_histo->frequency_distribution->offset_computation();
    compound_histo->frequency_distribution->nb_element_computation();
    compound_histo->frequency_distribution->max_computation();
    compound_histo->frequency_distribution->mean_computation();
    compound_histo->frequency_distribution->variance_computation();
  }

  return compound_histo;
}


};  // namespace stat_tool
