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

#include "stat_tool/stat_label.h"

#include "renewal.h"
#include "sequence_label.h"

using namespace std;
using namespace stat_tool;


namespace sequence_analysis {



/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the information quantity of a TimeEvents object.
 */
/*--------------------------------------------------------------*/

double TimeEvents::information_computation() const

{
  int i;
  double information , buff;


  information = htime->information_computation();

  if (information != D_INF) {
    for (i = htime->offset;i < htime->nb_value;i++) {
      if (htime->frequency[i] > 0) {
        buff = hnb_event[i]->information_computation();

        if (buff != D_INF) {
          information += buff;
        }
        else {
          information = D_INF;
          break;
        }
      }
    }
  }

  return information;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the log-likelihood of a mixture of number of events distributions
 *         for pairs {observation period, number of events}.
 *
 *  \param[in] timev reference on a TimeEvents object.
 */
/*--------------------------------------------------------------*/

double Renewal::likelihood_computation(const TimeEvents &timev) const

{
  int i;
  double likelihood , buff;


  likelihood = time->likelihood_computation(*(timev.htime));

  if (likelihood != D_INF) {
    for (i = timev.htime->offset;i < timev.htime->nb_value;i++) {
      if (timev.htime->frequency[i] > 0) {
        buff = nb_event[i]->likelihood_computation(*(timev.hnb_event[i]));

        if (buff != D_INF) {
          likelihood += buff;
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
 *  \brief Computation of the inter-event distribution reestimation quantities
 *         (EM estimator of an ordinary renewal process on the basis of count data).
 *
 *  \param[in] timev               reference on the pairs {observation period, number of events},
 *  \param[in] inter_event_reestim pointer on the reestimation quantities.
 */
/*--------------------------------------------------------------*/

void Renewal::expectation_step(const TimeEvents &timev ,
                               Reestimation<double> *inter_event_reestim) const

{
  int i , j;
  int min_time , max_time , *ptime , *pnb_event , *pfrequency;
  double num , denom , *ifrequency , *pmass;


  // initialization

  ifrequency = inter_event_reestim->frequency;
  for (i = 0;i < inter_event_reestim->alloc_nb_value;i++) {
    *ifrequency++ = 0.;
  }

  ptime = timev.time;
  pnb_event = timev.nb_event;
  pfrequency = timev.frequency;

  for (i = 0;i < timev.nb_class;i++) {

    // case no event

    if (*pnb_event == 0) {
      min_time = MAX(inter_event->offset , *ptime + 1);

      if (min_time < inter_event->nb_value) {
        denom = nb_event[*ptime]->mass[*pnb_event];

        if (denom > 0.) {
          ifrequency = inter_event_reestim->frequency + min_time;
          pmass = inter_event->mass + min_time;
          for (j = min_time;j < inter_event->nb_value;j++) {
            *ifrequency++ += *pfrequency * *pmass++ / denom;
          }
        }
      }
    }

    // case number of events > 0

    else {
      denom = 0.;
      if (*pnb_event < nb_event[*ptime]->nb_value) {
        denom = nb_event[*ptime]->mass[*pnb_event];
      }

      if (denom > 0.) {
        max_time = MIN(nevent_time[*pnb_event]->nb_value - 1 , *ptime);

        ifrequency = inter_event_reestim->frequency + inter_event->offset;
        pmass = inter_event->mass + inter_event->offset;

        for (j = inter_event->offset;j < inter_event->nb_value;j++) {
          num = 0.;

          // case number of events = 1: complet time interval

          if (*pnb_event == 1) {
            if ((*ptime - j >= 0) && (*ptime - j < nevent_time[*pnb_event]->nb_value)) {
              num = 1. - nevent_time[*pnb_event]->cumul[*ptime - j];
            }
          }

          // case number of events > 1: complete time intervals

          else {
            if (*ptime - j >= nevent_time[*pnb_event - 1]->offset) {
              if (*ptime - j < nevent_time[*pnb_event - 1]->nb_value) {
                num = *pnb_event * (nevent_time[*pnb_event - 1]->cumul[*ptime - j] -
                                    nevent_time[*pnb_event]->cumul[*ptime - j]);
              }
              else if (*ptime - j < nevent_time[*pnb_event]->nb_value) {
                num = *pnb_event * (1. - nevent_time[*pnb_event]->cumul[*ptime - j]);
              }
            }
          }

          // right-censored time interval

          if ((max_time > *ptime - j) && (max_time >= nevent_time[*pnb_event]->offset)) {
            num += nevent_time[*pnb_event]->cumul[max_time];
            if (*ptime - j >= nevent_time[*pnb_event]->offset) {
              num -= nevent_time[*pnb_event]->cumul[*ptime - j];
            }
          }

          *ifrequency++ += *pfrequency * *pmass++ * num / denom;
        }
      }
    }

    ptime++;
    pnb_event++;
    pfrequency++;
  }

  inter_event_reestim->nb_value_computation();
  inter_event_reestim->offset_computation();
  inter_event_reestim->nb_element_computation();
  inter_event_reestim->max_computation();
  inter_event_reestim->mean_computation();
  inter_event_reestim->variance_computation();

# ifdef DEBUG
  cout << "\n" << (timev.mixture->mean + 1) * timev.nb_element << " | "
       << " " << inter_event_reestim->nb_element << endl;

  cout << "\nthe reestimation quantities inter-event distribution:" << *inter_event_reestim << endl;
# endif

}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the inter-event distribution reestimation quantities
 *         (EM estimator of an equilibrium renewal process on the basis of count data).
 *
 *  \param[in] timev               reference on the pairs {observation period, number of events},
 *  \param[in] inter_event_reestim pointer on the inter-event distribution reestimation quantities,
 *  \param[in] length_bias_reestim pointer on the length-biased distribution reestimation quantities,
 *  \param[in] estimator           estimator type (complete or partial likelihood),
 *  \param[in] combination         combination or not of the reestimation quantities,
 *  \param[in] mean_estimator      method for the computation of the inter-event distribution mean (equilibrium renewel process).
 */
/*--------------------------------------------------------------*/

void Renewal::expectation_step(const TimeEvents &timev ,
                               Reestimation<double> *inter_event_reestim ,
                               Reestimation<double> *length_bias_reestim ,
                               censoring_estimator estimator , bool combination ,
                               duration_distribution_mean_estimator mean_estimator) const

{
  int i , j;
  int max_time , offset , nb_value , *ptime , *pnb_event , *pfrequency;
  double complete_num , censored_num , denom , inter_event_mean , *ifrequency ,
         *lfrequency , *pmass;


  // initializations

  ifrequency = inter_event_reestim->frequency;
  for (i = 0;i < inter_event_reestim->alloc_nb_value;i++) {
    *ifrequency++ = 0.;
  }

  if (estimator == COMPLETE_LIKELIHOOD) {
    lfrequency = length_bias_reestim->frequency;
    for (i = 0;i < length_bias_reestim->alloc_nb_value;i++) {
      *lfrequency++ = 0.;
    }
  }

  ptime = timev.time;
  pnb_event = timev.nb_event;
  pfrequency = timev.frequency;

  for (i = 0;i < timev.nb_class;i++) {

    // case no event

    if (*pnb_event == 0) {
      if ((estimator == COMPLETE_LIKELIHOOD) && (*ptime + 1 < inter_event->nb_value)) {
        denom = nb_event[*ptime]->mass[*pnb_event] * inter_event->mean;

        if (denom > 0.) {
          lfrequency = length_bias_reestim->frequency + *ptime + 1;
          pmass = inter_event->mass + *ptime + 1;
          for (j = *ptime + 1;j < inter_event->nb_value;j++) {
            *lfrequency++ += *pfrequency * (j - *ptime) * *pmass++ / denom;
          }
        }
      }
    }

    // case number of events > 0

    else {
      denom = 0.;
      if (*pnb_event < nb_event[*ptime]->nb_value) {
        denom = nb_event[*ptime]->mass[*pnb_event];
      }

      if (denom > 0.) {

        // left-censored time interval

/*        if (estimator == COMPLETE_LIKELIHOOD) {
          lfrequency = length_bias_reestim->frequency + inter_event->offset;
          pmass = inter_event->mass + inter_event->offset;
          num = 0.;

          for (j = 1;j < inter_event->nb_value;j++) {
            if (j <= *ptime) {

              // case 1 event

              if (*pnb_event == 1) {
                if (*ptime - j < aux_nevent_time[*pnb_event]->nb_value) {
                  num += 1. - aux_nevent_time[*pnb_event]->cumul[*ptime - j];
                }
              }

              // case number of events > 1

              else {
                if (*ptime - j >= aux_nevent_time[*pnb_event - 1]->offset) {
                  if (*ptime - j < aux_nevent_time[*pnb_event - 1]->nb_value) {
                    num += aux_nevent_time[*pnb_event - 1]->cumul[*ptime - j] -
                           aux_nevent_time[*pnb_event]->cumul[*ptime - j];
                  }
                  else if (*ptime - j < aux_nevent_time[*pnb_event]->nb_value) {
                    num += 1. - aux_nevent_time[*pnb_event]->cumul[*ptime - j];
                  }
                }
              }
            }

            if (j >= inter_event->offset) {
              *lfrequency++ += *pfrequency * *pmass++ * num / (denom * inter_event->mean);
            }
          }
        } */

        max_time = MIN(nevent_time[*pnb_event]->nb_value - 1 , *ptime);

        ifrequency = inter_event_reestim->frequency + inter_event->offset;
        if (estimator == COMPLETE_LIKELIHOOD) {
          lfrequency = length_bias_reestim->frequency + inter_event->offset;
        }
        pmass = inter_event->mass + inter_event->offset;

        for (j = inter_event->offset;j < inter_event->nb_value;j++) {
          complete_num = 0.;

          // case number of events > 1: complete time intervals

          if (*pnb_event > 1) {
            if (*ptime - j >= nevent_time[*pnb_event - 1]->offset) {
              if (*ptime - j < nevent_time[*pnb_event - 1]->nb_value) {
                complete_num = (*pnb_event - 1) * (nevent_time[*pnb_event - 1]->cumul[*ptime - j] -
                                                   nevent_time[*pnb_event]->cumul[*ptime - j]);
              }
              else if (*ptime - j < nevent_time[*pnb_event]->nb_value) {
                complete_num = (*pnb_event - 1) * (1. - nevent_time[*pnb_event]->cumul[*ptime - j]);
              }
            }
          }

          // left- and right-censored time intervals

          censored_num = 0.;
          if ((max_time > *ptime - j) && (max_time >= nevent_time[*pnb_event]->offset)) {
            censored_num += nevent_time[*pnb_event]->cumul[max_time];
            if (*ptime - j >= nevent_time[*pnb_event]->offset) {
              censored_num -= nevent_time[*pnb_event]->cumul[*ptime - j];
            }
          }

          *ifrequency++ += *pfrequency * *pmass * (complete_num + censored_num) / denom;
          if (estimator == COMPLETE_LIKELIHOOD) {
            *lfrequency++ += *pfrequency * *pmass * censored_num / denom;
          }
          pmass++;
        }
      }
    }

    ptime++;
    pnb_event++;
    pfrequency++;
  }

  switch (estimator) {

  case PARTIAL_LIKELIHOOD : {
    inter_event_reestim->nb_value_computation();
    inter_event_reestim->offset_computation();
    inter_event_reestim->nb_element_computation();
    break;
  }

  case COMPLETE_LIKELIHOOD : {
    offset = 1;
    nb_value = inter_event_reestim->alloc_nb_value;

    ifrequency = inter_event_reestim->frequency + inter_event_reestim->alloc_nb_value;
    lfrequency = length_bias_reestim->frequency + inter_event_reestim->alloc_nb_value;
    while ((*--ifrequency == 0) && (*--lfrequency == 0) && (nb_value > 2)) {
      nb_value--;
    }
    inter_event_reestim->nb_value = nb_value;
    length_bias_reestim->nb_value = nb_value;

    ifrequency = inter_event_reestim->frequency + offset;
    lfrequency = length_bias_reestim->frequency + offset;
    while ((*ifrequency++ == 0) && (*lfrequency++ == 0) && (offset < nb_value - 1)) {
      offset++;
    }
    inter_event_reestim->offset = offset;
    length_bias_reestim->offset = offset;

    inter_event_reestim->nb_element_computation();
    length_bias_reestim->nb_element_computation();
    break;

#   ifdef DEBUG
    inter_event_reestim->max_computation();
    inter_event_reestim->mean_computation();
    inter_event_reestim->variance_computation();

    length_bias_reestim->max_computation();
    length_bias_reestim->mean_computation();
    length_bias_reestim->variance_computation();

    cout << "\n" << timev.nb_element << " | "  << length_bias_reestim->nb_element << " || "
         << timev.mixture->mean * timev.nb_element << " | " << " " << inter_event_reestim->nb_element << endl;

    cout << "\nthe reestimation quantities inter-event distribution:" << *inter_event_reestim << endl;
    cout << "\nthe reestimation quantities length-biased distribution:" << *length_bias_reestim << endl;
#   endif

  }
  }

  if ((estimator == COMPLETE_LIKELIHOOD) && (combination)) {
    switch (mean_estimator) {
    case ESTIMATED :
      inter_event_mean = timev.htime->mean / timev.mixture->mean;
      break;
    case COMPUTED :
      inter_event_mean = interval_bisection(inter_event_reestim , length_bias_reestim);
      break;
    case ONE_STEP_LATE :
      inter_event_mean = inter_event->mean;
      break;
    }

#   ifdef DEBUG
    if (mean_estimator != ESTIMATED) {
      cout << SEQ_label[SEQL_INTER_EVENT] << " " << STAT_label[STATL_MEAN] << ": "
           << inter_event_mean << " (" << timev.htime->mean / timev.mixture->mean << ") | ";
    }
#   endif

    inter_event_reestim->equilibrium_process_combination(length_bias_reestim , inter_event_mean);
  }

  else {
    inter_event_reestim->max_computation();
    inter_event_reestim->mean_computation();
    inter_event_reestim->variance_computation();
  }

# ifdef DEBUG
  cout << "\nthe reestimation quantities inter-event distribution:" << *inter_event_reestim << endl;
# endif

}


/*--------------------------------------------------------------*/
/**
 *  \brief Estimation of a renewal process on the basis of pairs
 *         {observation period, number of events} using the EM algorithm.
 *
 *  \param[in] error                 reference on a StatError object,
 *  \param[in] display               flag for displaying estimation intermediate results,
 *  \param[in] type                  renewal process type (ORDINARY/EQUILIBRIUM),
 *  \param[in] iinter_event          reference on the initial inter-event distribution,
 *  \param[in] estimator             estimator type (maximum likelihood or penalized likelihood or
 *                                   estimation of a parametric distribution),
 *  \param[in] nb_iter               number of iterations,
 *  \param[in] equilibrium_estimator estimator type in the case of an equilibrium renewal process,
 *  \param[in] mean_estimator        method of computation of the inter-event distribution mean,
 *  \param[in] weight                penalty weight,
 *  \param[in] pen_type              penalty type,
 *  \param[in] outside               management of side effects (zero outside the support or
 *                                   continuation of the distribution).
 *
 *  \return                          Renewal object.
 */
/*--------------------------------------------------------------*/

Renewal* TimeEvents::estimation(StatError &error , bool display , process_type type ,
                                const DiscreteParametric &iinter_event , estimation_criterion estimator ,
                                int nb_iter , censoring_estimator equilibrium_estimator ,
                                duration_distribution_mean_estimator mean_estimator ,
                                double weight , penalty_type pen_type , side_effect outside) const

{
  bool status = true;
  int i;
  int nb_likelihood_decrease;
  double likelihood , previous_likelihood , information , hlikelihood , inter_event_mean , *penalty;
  DiscreteParametric *pinter_ev;
  Reestimation<double> *inter_event_reestim , *length_bias_reestim;
  FrequencyDistribution *hreestim;
  Renewal *renew;


  renew = NULL;
  error.init();

  if (mixture->nb_value <= 2) {
    status = false;
    error.update(SEQ_error[SEQR_MAX_NB_EVENT_TOO_SMALL]);
  }
  if (mixture->mean < MIN_NB_EVENT) {
    status = false;
    error.update(SEQ_error[SEQR_NB_EVENT_TOO_SMALL]);
  }

  if (min_inter_event_computation() < MIN_INTER_EVENT) {
    status = false;
    error.update(SEQ_error[SEQR_TIME_UNIT]);
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
    information = information_computation();

    // construction of a Renewal object

    renew = new Renewal(type , *htime , iinter_event);
    renew->renewal_data = new RenewalData(*this , type);

    pinter_ev = renew->inter_event;

    if (estimator == PENALIZED_LIKELIHOOD) {
      penalty = new double[pinter_ev->alloc_nb_value];

      if (weight == D_DEFAULT) {
        if (pen_type != ENTROPY) {
          weight = RENEWAL_DIFFERENCE_WEIGHT;
        }
        else {
          weight = RENEWAL_ENTROPY_WEIGHT;
        }
      }

      if (equilibrium_estimator == PARTIAL_LIKELIHOOD) {
        weight *= mixture->mean * nb_element;
      }
      else {
        weight *= (mixture->mean + 1) * nb_element;
      }
    }

    inter_event_reestim = new Reestimation<double>(pinter_ev->alloc_nb_value);

    if (type == EQUILIBRIUM) {
      if (equilibrium_estimator == COMPLETE_LIKELIHOOD) {
        length_bias_reestim = new Reestimation<double>(pinter_ev->alloc_nb_value);
      }
    }

    renew->computation();

    renew->init(CATEGORICAL , I_DEFAULT , I_DEFAULT , D_DEFAULT , D_DEFAULT);
    likelihood = D_INF;
    i = 0;

    do {
      i++;

      // computation of the reestimation quantities

      switch (type) {
      case ORDINARY :
        renew->expectation_step(*this , inter_event_reestim);
        break;
      case EQUILIBRIUM :
        renew->expectation_step(*this , inter_event_reestim ,
                                length_bias_reestim , equilibrium_estimator);
        break;
      }

      if (estimator != PENALIZED_LIKELIHOOD) {
        if ((type == ORDINARY) || (equilibrium_estimator == PARTIAL_LIKELIHOOD)) {
          inter_event_reestim->distribution_estimation(pinter_ev);
        }

        else {
          switch (mean_estimator) {
          case ESTIMATED :
            inter_event_mean = htime->mean / mixture->mean;
            break;
          case COMPUTED :
            inter_event_mean = interval_bisection(inter_event_reestim , length_bias_reestim);
            break;
          case ONE_STEP_LATE :
            inter_event_mean = pinter_ev->mean;
            break;
          }

          inter_event_reestim->equilibrium_process_estimation(length_bias_reestim , pinter_ev ,
                                                              inter_event_mean);
        }
      }

      else {
        if ((type == ORDINARY) || (equilibrium_estimator == PARTIAL_LIKELIHOOD)) {
          inter_event_reestim->penalized_likelihood_estimation(pinter_ev , weight ,
                                                               pen_type , penalty ,
                                                               outside);
        }

        else {
          switch (mean_estimator) {
          case ESTIMATED :
            inter_event_mean = htime->mean / mixture->mean;
            break;
          case ONE_STEP_LATE :
            inter_event_mean = pinter_ev->mean;
            break;
          }

          inter_event_reestim->penalized_likelihood_equilibrium_process_estimation(length_bias_reestim ,
                                                                                   pinter_ev , inter_event_mean ,
                                                                                   weight , pen_type ,
                                                                                   penalty , outside);
        }
      }

      // computation of the mixture of number of events distributions and the associated log-likelihood

      renew->computation();

      previous_likelihood = likelihood;
      likelihood = renew->likelihood_computation(*this);

      if ((display) && ((i < 10) || ((i < 100) && (i % 10 == 0)) || ((i < 1000) && (i % 100 == 0)) || (i % 1000 == 0))) {
        cout << STAT_label[STATL_ITERATION] << " " << i << "   "
             << STAT_label[STATL_LIKELIHOOD] << ": " << likelihood << "   "
             << STAT_label[STATL_DEVIANCE] << ": " << 2 * (information - likelihood) << "   "
             << STAT_label[STATL_SMOOTHNESS] << ": "  << pinter_ev->second_difference_norm_computation();
        if (estimator == PENALIZED_LIKELIHOOD) {
          cout << "   cumul: " << pinter_ev->cumul[pinter_ev->nb_value - 1];
        }
        cout << endl;
      }
    }
    while ((likelihood != D_INF) && (((nb_iter == I_DEFAULT) && (i < RENEWAL_NB_ITER) &&
             ((likelihood - previous_likelihood) / -likelihood > RENEWAL_LIKELIHOOD_DIFF)) ||
            ((nb_iter != I_DEFAULT) && (i < nb_iter))));

    if (likelihood != D_INF) {
      if (display) {
        cout << "\n" << i << " " << STAT_label[STATL_ITERATIONS] << "   "
             << STAT_label[STATL_LIKELIHOOD] << ": " << likelihood << "   "
             << STAT_label[STATL_DEVIANCE] << ": " << 2 * (information - likelihood) << "   "
             << STAT_label[STATL_SMOOTHNESS] << ": " << pinter_ev->second_difference_norm_computation();
        if (estimator == PENALIZED_LIKELIHOOD) {
          cout << "   cumul: " << pinter_ev->cumul[pinter_ev->nb_value - 1];
        }
        cout << endl;
      }

      if (estimator == PARAMETRIC_REGULARIZATION) {
        hreestim = new FrequencyDistribution(pinter_ev->alloc_nb_value);

        likelihood = D_INF;
        nb_likelihood_decrease = 0;

        i = 0;
        do {
          i++;

          // computation of the reestimation quantities

          switch (type) {
          case ORDINARY :
            renew->expectation_step(*this , inter_event_reestim);
            break;
          case EQUILIBRIUM :
            renew->expectation_step(*this , inter_event_reestim ,
                                    length_bias_reestim , equilibrium_estimator ,
                                    true , mean_estimator);
            break;
          }

          hreestim->update(inter_event_reestim , (int)(inter_event_reestim->nb_element *
                           MAX(sqrt(inter_event_reestim->variance) , 1.) * RENEWAL_COEFF));
          hlikelihood = hreestim->Reestimation<int>::type_parametric_estimation(pinter_ev , 1 , true ,
                                                                                RENEWAL_THRESHOLD);

          if (hlikelihood == D_INF) {
            likelihood = D_INF;
          }

          // computation of the mixture of number of events distributions and the associated log-likelihood

          else {
            renew->init(pinter_ev->ident , pinter_ev->inf_bound , pinter_ev->sup_bound ,
                        pinter_ev->parameter , pinter_ev->probability);
            renew->computation();

            previous_likelihood = likelihood;
            likelihood = renew->likelihood_computation(*this);

            if (likelihood < previous_likelihood) {
              nb_likelihood_decrease++;
            }
            else {
              nb_likelihood_decrease = 0;
            }

#           ifdef DEBUG
            if ((i < 10) || (i % 10 == 0)) {
              cout << STAT_label[STATL_ITERATION] << " " << i << "   "
                   << STAT_label[STATL_LIKELIHOOD] << ": " << likelihood << "   "
                   << STAT_label[STATL_SMOOTHNESS] << ": " << pinter_ev->second_difference_norm_computation() << endl;
            }
#           endif

          }
        }
        while ((likelihood != D_INF) && (i < RENEWAL_NB_ITER) &&
               (((likelihood - previous_likelihood) / -likelihood > RENEWAL_LIKELIHOOD_DIFF) ||
                (hlikelihood == D_INF) || (nb_likelihood_decrease == 1)));

        delete hreestim;

        if ((display) && (likelihood != D_INF)) {
          cout << "\n" << i << " " << STAT_label[STATL_ITERATIONS] << "   "
               << STAT_label[STATL_LIKELIHOOD] << ": " << likelihood << "   "
               << STAT_label[STATL_SMOOTHNESS] << ": " << pinter_ev->second_difference_norm_computation() << endl;
        }
      }
    }

    if (estimator == PENALIZED_LIKELIHOOD) {
      delete [] penalty;
    }

    delete inter_event_reestim;

    if (type == EQUILIBRIUM) {
      if (equilibrium_estimator == COMPLETE_LIKELIHOOD) {
        delete length_bias_reestim;
      }
    }

    if (likelihood != D_INF) {

      // update of the number of free parameters

      pinter_ev->nb_parameter_update();
      for (i = renew->time->offset;i < renew->time->nb_value;i++) {
        if (renew->time->mass[i] > 0.) {
          renew->nb_event[i]->nb_parameter = pinter_ev->nb_parameter;
        }
      }
      renew->mixture->nb_parameter = pinter_ev->nb_parameter;
    }

    else {
      delete renew;
      renew = NULL;
      error.update(STAT_error[STATR_ESTIMATION_FAILURE]);
    }
  }

  return renew;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Estimation a renewal process on the basis of pairs
 *         {observation period, number of events} using the EM algorithm.
 *
 *  \param[in] error                 reference on a StatError object,
 *  \param[in] display               flag for displaying estimation intermediate results,
 *  \param[in] type                  renewal process type (ORDINARY/EQUILIBRIUM),
 *  \param[in] estimator type        (maximum likelihood or penalized likelihood or
 *                                   estimation of a parametric distribution),
 *  \param[in] nb_iter               number of iterations,
 *  \param[in] equilibrium_estimator estimator type in the case of an equilibrium renewal process,
 *  \param[in] mean_estimator        method of computation of the inter-event distribution mean,
 *  \param[in] weight                penalty weight,
 *  \param[in] pen_type              penalty type,
 *  \param[in] outside               management of side effects (zero outside the support or
 *                                   continuation of the distribution).
 *
 *  \return                          Renewal object.
 */
/*--------------------------------------------------------------*/

Renewal* TimeEvents::estimation(StatError &error , bool display , process_type type ,
                                estimation_criterion estimator , int nb_iter ,
                                censoring_estimator equilibrium_estimator ,
                                duration_distribution_mean_estimator mean_estimator , double weight ,
                                penalty_type pen_type , side_effect outside) const

{
  double proba;
  DiscreteParametric *iinter_event;
  Renewal *renew;


  proba = mixture->mean / htime->mean;
  if (proba > 1. - RENEWAL_INIT_PROBABILITY) {
    proba = 1. - RENEWAL_INIT_PROBABILITY;
  }
  else if (proba < RENEWAL_INIT_PROBABILITY) {
    proba = RENEWAL_INIT_PROBABILITY;
  }

  iinter_event = new DiscreteParametric(NEGATIVE_BINOMIAL , 1 , I_DEFAULT , 1. ,
                                        proba , RENEWAL_THRESHOLD);

# ifdef DEBUG
  iinter_event->ascii_print(cout);
# endif

  renew = estimation(error , display , type , *iinter_event , estimator , nb_iter ,
                     equilibrium_estimator , mean_estimator , weight ,
                     pen_type , outside);
  delete iinter_event;

  return renew;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Estimation of an equilibrium renewal process on the basis of
 *         time interval data using the EM algorithm.
 *
 *  \param[in] error          reference on a StatError object,
 *  \param[in] display        flag for displaying estimation intermediate results,
 *  \param[in] iinter_event   reference on the initial inter-event distribution,
 *  \param[in] estimator      estimator type (maximum likelihood or penalized likelihood),
 *  \param[in] nb_iter        number of iterations,
 *  \param[in] mean_estimator method of computation of the inter-event distribution mean,
 *  \param[in] weight         penalty weight,
 *  \param[in] pen_type       penalty type,
 *  \param[in] outside        management of side effects (zero outside the support or
 *                            continuation of the distribution).
 *
 *  \return                   Renewal object.
 */
/*--------------------------------------------------------------*/

Renewal* RenewalData::estimation(StatError &error , bool display ,
                                 const DiscreteParametric &iinter_event ,
                                 estimation_criterion estimator , int nb_iter ,
                                 duration_distribution_mean_estimator mean_estimator , double weight ,
                                 penalty_type pen_type , side_effect outside) const

{
  int i , j;
  int *psequence;
  DiscreteParametricModel *inter_event;
  Renewal *renew;
  FrequencyDistribution *within_backward , *within_forward , *no_event;


  within_backward = new FrequencyDistribution(MIN(backward->nb_value , htime->nb_value - 1));
  within_forward = new FrequencyDistribution(MIN(forward->nb_value , htime->nb_value));
  no_event = new FrequencyDistribution(htime->nb_value);

  for (i = 0;i < nb_element;i++) {
    psequence = sequence[i] + length[i];
    for (j = 0;j < length[i];j++) {
      if (*--psequence == 1) {
        (within_backward->frequency[j])++;
        break;
      }
    }

    psequence = sequence[i];
    for (j = 0;j < length[i];j++) {
      if (*psequence++ == 1) {
        (within_forward->frequency[j + 1])++;
        break;
      }
    }

    if (j == length[i]) {
      (no_event->frequency[j])++;
    }
  }

  within_backward->nb_value_computation();
  within_backward->offset_computation();
  within_backward->nb_element_computation();

# ifdef DEBUG
  within_backward->max_computation();
  within_backward->mean_computation();
  within_backward->variance_computation();

  cout << *within_backward;
  cout << *backward;
# endif

  within_forward->nb_value_computation();
  within_forward->offset_computation();
  within_forward->nb_element_computation();

# ifdef DEBUG
  within_forward->max_computation();
  within_forward->mean_computation();
  within_forward->variance_computation();

  cout << *within_forward;
  cout << *forward;
# endif

  no_event->nb_value_computation();
  no_event->offset_computation();
  no_event->nb_element_computation();

# ifdef DEBUG
  no_event->max_computation();
  no_event->mean_computation();
  no_event->variance_computation();

  cout << *no_event;
# endif

  if (no_event->nb_element == 0) {
    delete no_event;
    no_event = NULL;
  }

  inter_event = within->estimation(error , display , *within_backward , *within_forward , no_event ,
                                   iinter_event , estimator , nb_iter , mean_estimator ,
                                   weight , pen_type , outside , htime->mean / mixture->mean);

  delete within_backward;
  delete within_forward;
  delete no_event;

  if (inter_event) {
    renew = new Renewal(*this , *((DiscreteParametric*)inter_event));
  }
  else {
    renew = NULL;
  }

  return renew;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Estimation of an equilibrium renewal process on the basis of
 *         time interval data using the EM algorithm.
 *
 *  \param[in] error          reference on a StatError object,
 *  \param[in] display        flag for displaying estimation intermediate results,
 *  \param[in] estimator      estimator type (maximum likelihood or penalized likelihood),
 *  \param[in] nb_iter        number of iterations,
 *  \param[in] mean_estimator method of computation of the inter-event distribution mean,
 *  \param[in] weight         penalty weight,
 *  \param[in] pen_type       penalty type,
 *  \param[in] outside        management of side effects (zero outside the support or
 *                            continuation of the distribution).
 *
 *  \return                   Renewal object.
 */
/*--------------------------------------------------------------*/

Renewal* RenewalData::estimation(StatError &error , bool display , estimation_criterion estimator ,
                                 int nb_iter , duration_distribution_mean_estimator mean_estimator ,
                                 double weight , penalty_type pen_type , side_effect outside) const

{
  double proba;
  DiscreteParametric *iinter_event;
  Renewal *renew;


  proba = mixture->mean / htime->mean;
  if (proba > 1. - RENEWAL_INIT_PROBABILITY) {
    proba = 1. - RENEWAL_INIT_PROBABILITY;
  }
  else if (proba < RENEWAL_INIT_PROBABILITY) {
    proba = RENEWAL_INIT_PROBABILITY;
  }

  iinter_event = new DiscreteParametric(NEGATIVE_BINOMIAL , 1 , I_DEFAULT , 1. ,
                                        proba , RENEWAL_THRESHOLD);

# ifdef DEBUG
  iinter_event->ascii_print(cout);
# endif

  renew = estimation(error , display , *iinter_event , estimator , nb_iter ,
                     mean_estimator , weight , pen_type , outside);
  delete iinter_event;

  return renew;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Simulation using a renewal process.
 *
 *  \param[in] error  reference on a StatError object,
 *  \param[in] itype  renewal process type (ORDINARY/EQUILIBRIUM),
 *  \param[in] ihtime observation frequency distribution.
 *
 *  \return           RenewalData object.
 */
/*--------------------------------------------------------------*/

RenewalData* Renewal::simulation(StatError &error , process_type itype ,
                                 const FrequencyDistribution &ihtime) const

{
  bool status = true , compute;
  int i , j , k , m;
  int offset , time_interval , cumul_time , *ptime , *pnb_event , *psequence;
  Distribution *dtime;
  Renewal *renew;
  RenewalData *timev;


  timev = NULL;
  error.init();

  if ((ihtime.nb_element < 1) || (ihtime.nb_element > RENEWAL_NB_ELEMENT)) {
    status = false;
    error.update(STAT_error[STATR_SAMPLE_SIZE]);
  }
  if (ihtime.offset < MAX(inter_event->offset , 2)) {
    status = false;
    error.update(SEQ_error[SEQR_SHORT_OBSERVATION_TIME]);
  }
  if (ihtime.nb_value - 1 > MAX_TIME) {
    status = false;
    error.update(SEQ_error[SEQR_LONG_OBSERVATION_TIME]);
  }

  if (status) {
    dtime = new Distribution(ihtime);

    if ((itype != type) || (*dtime != *time)) {
      compute = true;
    }
    else {
      compute = false;
    }

    // construction of a RenewalData object

    timev = new RenewalData(itype , *this);

    timev->renewal = new Renewal(*this , false);
    renew = timev->renewal;

    timev->length = new int[ihtime.nb_element];
    timev->sequence = new int*[ihtime.nb_element];

    switch (itype) {
    case ORDINARY :
      offset = 0;
      break;
    case EQUILIBRIUM :
      offset = 1;
      break;
    }

    // 1st to nth event

    ptime = new int[ihtime.nb_element];
    pnb_event = new int[ihtime.nb_element];

    i = 0;
    for (j = ihtime.offset;j < ihtime.nb_value;j++) {
      for (k = 0;k < ihtime.frequency[j];k++) {

        // time to the 1st event (equilibrium renewal process)

        if (itype == EQUILIBRIUM) {
          if (i == 0) {
            cumul_time = renew->forward->simulation();
          }
          else {
            cumul_time -= *(ptime - 1);
          }
          (timev->forward->frequency[cumul_time])++;

          time_interval = cumul_time;
        }

        // observation period

        *ptime = j;

        timev->length[i] = *ptime + 1 - offset;
        timev->sequence[i] = new int[timev->length[i]];

        // time to the 1st event (ordinary renewal process)

        if (itype == ORDINARY) {
          time_interval = renew->inter_event->simulation();
          cumul_time = time_interval;
          (timev->inter_event->frequency[time_interval])++;
          if (time_interval <= *ptime) {
            (timev->within->frequency[time_interval])++;
          }
          else {
            (timev->length_bias->frequency[time_interval])++;
          }
        }

        psequence = timev->sequence[i] - 1;
        if (itype == ORDINARY) {
          *++psequence = 1;
        }
        for (m = 1;m < MIN(time_interval , *ptime + 1);m++) {
          *++psequence = 0;
        }
        if (time_interval <= *ptime) {
          *++psequence = 1;
        }

        *pnb_event = 0;
        while (cumul_time <= *ptime) {
          (*pnb_event)++;
          time_interval = renew->inter_event->simulation();
          cumul_time += time_interval;

          for (m = cumul_time - time_interval + 1;m < MIN(cumul_time , *ptime + 1);m++) {
            *++psequence = 0;
          }

          (timev->inter_event->frequency[time_interval])++;
          if (cumul_time <= *ptime) {
            (timev->within->frequency[time_interval])++;
            *++psequence = 1;
          }
          else {
            (timev->length_bias->frequency[time_interval])++;
          }
        }

        (timev->backward->frequency[*ptime - (cumul_time - time_interval)])++;
        if (itype == ORDINARY) {
          (timev->forward->frequency[cumul_time - *ptime])++;
        }

        ptime++;
        pnb_event++;
        i++;
      }
    }
    ptime -= ihtime.nb_element;
    pnb_event -= ihtime.nb_element;

    // construction of the triplets {observation period, number of events, frequency} and of
    // the observation period frequency distribution and number of events frequency distributions

    timev->build(ihtime.nb_element , ptime , pnb_event);
    delete [] ptime;
    delete [] pnb_event;

    // extraction of the characteristics of the inter-event frequency distribution,
    // the frequency distribution of time intervals between events within the observation period,
    // the length-biased frequency distribution,
    // the backward and forward recurrence time frequency distributions,

    timev->inter_event->nb_value_computation();
    timev->inter_event->offset_computation();
    timev->inter_event->nb_element_computation();
    timev->inter_event->max_computation();
    timev->inter_event->mean_computation();
    timev->inter_event->variance_computation();

    timev->within->nb_value_computation();
    timev->within->offset_computation();
    timev->within->nb_element_computation();
    timev->within->max_computation();
    timev->within->mean_computation();
    timev->within->variance_computation();

    timev->length_bias->nb_value_computation();
    timev->length_bias->offset_computation();
    timev->length_bias->nb_element_computation();
    timev->length_bias->max_computation();
    timev->length_bias->mean_computation();
    timev->length_bias->variance_computation();

    timev->backward->nb_value_computation();
    timev->backward->offset_computation();
    timev->backward->nb_element_computation();
    timev->backward->max_computation();
    timev->backward->mean_computation();
    timev->backward->variance_computation();

    timev->forward->nb_value_computation();
    timev->forward->offset_computation();
    timev->forward->nb_element_computation();
    timev->forward->max_computation();
    timev->forward->mean_computation();
    timev->forward->variance_computation();

    // extraction of no-event/event probabilities as a function of time

    timev->build_index_event(offset);

    if (compute) {
      renew->computation(false , itype , dtime);
    }
    delete dtime;
  }

  return timev;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Simulation using a renewal process.
 *
 *  \param[in] error      reference on a StatError object,
 *  \param[in] itype      renewal process type (ORDINARY/EQUILIBRIUM),
 *  \param[in] nb_element sample size,
 *  \param[in] itime      observation period.
 *
 *  \return               RenewalData object.
 */
/*--------------------------------------------------------------*/

RenewalData* Renewal::simulation(StatError &error , process_type itype ,
                                 int nb_element , int itime) const

{
  bool status = true;
  RenewalData *timev;


  timev = NULL;
  error.init();

  if ((nb_element < 1) || (nb_element > RENEWAL_NB_ELEMENT)) {
    status = false;
    error.update(STAT_error[STATR_SAMPLE_SIZE]);
  }
  if (itime < MAX(inter_event->offset , 2)) {
    status = false;
    error.update(SEQ_error[SEQR_SHORT_OBSERVATION_TIME]);
  }
  if (itime > MAX_TIME) {
    status = false;
    error.update(SEQ_error[SEQR_LONG_OBSERVATION_TIME]);
  }

  if (status) {
    FrequencyDistribution htime(itime + 1);

    htime.nb_element = nb_element;
    htime.offset = itime;
    htime.max = nb_element;
    htime.mean = itime;
    htime.variance = 0.;
    htime.frequency[itime] = nb_element;

    timev = simulation(error , itype , htime);
  }

  return timev;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Simulation using a renewal process.
 *
 *  \param[in] error      reference on a StatError object,
 *  \param[in] itype      renewal process type (ORDINARY/EQUILIBRIUM),
 *  \param[in] nb_element sample size,
 *  \param[in] itimev     reference on a TimeEvents object.
 *
 *  \return               RenewalData object.
 */
/*--------------------------------------------------------------*/

RenewalData* Renewal::simulation(StatError &error , process_type itype ,
                                 int nb_element , const TimeEvents &itimev) const

{
  FrequencyDistribution *htime;
  RenewalData *timev;


  error.init();

  if ((nb_element < 1) || (nb_element > RENEWAL_NB_ELEMENT)) {
    timev = NULL;
    error.update(STAT_error[STATR_SAMPLE_SIZE]);
  }

  else {
    htime = itimev.htime->frequency_scale(nb_element);

    timev = simulation(error , itype , *htime);
    delete htime;
  }

  return timev;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Constructor of the RenewalIterator class.
 *
 *  \param[in] irenewal pointer on a Renewal object,
 *  \param[in] ilength  sequence length.
 */
/*--------------------------------------------------------------*/

RenewalIterator::RenewalIterator(Renewal *irenewal , int ilength)

{
  renewal = irenewal;
  (renewal->nb_iterator)++;

  interval = 0;
  counter = 0;

  length = ilength;
  sequence = new int[length];
}


/*--------------------------------------------------------------*/
/**
 *  \brief Copy of a RenewalIterator object.
 *
 *  \param[in] iter reference on a RenewalIterator object.
 */
/*--------------------------------------------------------------*/

void RenewalIterator::copy(const RenewalIterator &iter)

{
  int i;
  int *psequence , *isequence;


  renewal = iter.renewal;
  (renewal->nb_iterator)++;

  interval = iter.interval;
  counter = iter.counter;
  length = iter.length;

  sequence = new int[length];

  psequence = sequence;
  isequence = iter.sequence;
  for (i = 0;i < length;i++) {
    *psequence++ = *isequence++;
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Destructor of the RenewalIterator class.
 */
/*--------------------------------------------------------------*/

RenewalIterator::~RenewalIterator()

{
  (renewal->nb_iterator)--;
  delete [] sequence;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Assignment operator of the RenewalIterator class.
 *
 *  \param[in] iter reference on a RenewalIterator object.
 *
 *  \return         RenewalIterator object.
 */
/*--------------------------------------------------------------*/

RenewalIterator& RenewalIterator::operator=(const RenewalIterator &iter)

{
  if (&iter != this) {
    (renewal->nb_iterator)--;
    delete [] sequence;

    copy(iter);
  }

  return *this;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Simulation using a renewal process.
 *
 *  \param[in] ilength sequence length,
 *  \param[in] type    renewal process type (ORDINARY/EQUILIBRIUM).
 */
/*--------------------------------------------------------------*/

void RenewalIterator::simulation(int ilength , process_type type)

{
  int i;
  int offset , *psequence;


  switch (type) {
  case ORDINARY :
    offset = 1;
    break;
  default :
    offset = 0;
    break;
  }

  if (ilength + offset != length) {
    length = ilength + offset;
    delete [] sequence;
    sequence = new int[length];
  }

  psequence = sequence;

  switch (type) {
  case ORDINARY :
    interval = renewal->inter_event->simulation();
    *psequence++ = 1;
    counter = 1;
    break;
  case EQUILIBRIUM :
    interval = renewal->forward->simulation();
    counter = 1;
    break;
  }

  for (i = offset;i < length;i++) {
    if (counter < interval) {
      *psequence++ = 0;
      counter++;
    }

    else {
      interval = renewal->inter_event->simulation();
      *psequence++ = 1;
      counter = 1;
    }
  }
}


};  // namespace sequence_analysis
