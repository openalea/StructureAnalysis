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

#include "nonhomogeneous_markov.h"
#include "sequence_label.h"

using namespace std;
using namespace stat_tool;


namespace sequence_analysis {



/*--------------------------------------------------------------*/
/**
 *  \brief Update of the transition distribution of a state for a nonhomogeneous Markov chain.
 *
 *  \param[in] state       state,
 *  \param[in] index       sequence index,
 *  \param[in] index_chain reference on the transition probabilities.
 */
/*--------------------------------------------------------------*/

void NonhomogeneousMarkov::transition_update(int state , int index , Chain &index_chain) const


{
  int i;
  double scale , *pparam;


  pparam = self_transition[state]->parameter;

  // update of the self-transition probability

  switch (self_transition[state]->ident) {
  case LOGISTIC :
    index_chain.transition[state][state] = pparam[0] / (1. + pparam[1] * exp(-pparam[2] * index));
    break;
  case MONOMOLECULAR :
    index_chain.transition[state][state] = pparam[0] + pparam[1] * exp(-pparam[2] * index);
    break;
  }

  // update of the state change probabilities

  scale = (1. - index_chain.transition[state][state]) / (1. - transition[state][state]);
  for (i = 0;i < nb_state;i++) {
    if (i != state) {
      index_chain.transition[state][i] = scale * transition[state][i];
    }
  }

  if (index_chain.cumul_transition) {
    stat_tool::cumul_computation(nb_state , index_chain.transition[state] ,
                                 index_chain.cumul_transition[state]);
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the state probabilities as a function of
 *         the index parameter for a nonhomogeneous Markov chain.
 */
/*--------------------------------------------------------------*/

void NonhomogeneousMarkov::index_state_distribution()

{
  int i , j , k;
  Curves *index_state;
  Chain *index_chain;


  index_state = process->index_value;

  // initialization of the transition probability matrix

  index_chain = new Chain(*this);

  // initialization of the state probabilities

  for (i = 0;i < nb_state;i++) {
    index_state->point[i][0] = initial[i];
  }

  for (i = 1;i < index_state->length;i++) {

    // change in transition probabilities with the index parameter

    for (j = 0;j < nb_state;j++) {
      if (!homogeneity[j]) {
        transition_update(j , i - 1 , *index_chain);
      }
    }

    // computation of the state probabilities

    for (j = 0;j < nb_state;j++) {
      index_state->point[j][i] = 0.;
      for (k = 0;k < nb_state;k++) {
        index_state->point[j][i] += index_chain->transition[k][j] * index_state->point[k][i - 1];
      }
    }
  }

  delete index_chain;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the probability of not visiting a state
 *         for a nonhomogeneous Markov chain.
 *
 *  \param[in] state     state,
 *  \param[in] increment threshold on the sum of the state probabilities.
 */
/*--------------------------------------------------------------*/

void NonhomogeneousMarkov::state_no_occurrence_probability(int state , double increment)

{
  int i;

  for (i = 0;i < nb_state;i++) {
    if ((i != state) && (!accessibility[i][state])) {
      break;
    }
  }

  if (i < nb_state) {
    int j , k;
    double state_sum , *current_state , *previous_state ,
           &no_occurrence = process->no_occurrence[state];
    Chain *index_chain;


    // initialization of the transition probability matrix

    index_chain = new Chain(*this);

    // initialization of the state probabilities

    current_state = new double[nb_state];
    previous_state = new double[nb_state];

    state_sum = 0.;
    no_occurrence = 0.;

    for (i = 0;i < nb_state;i++) {
      if (i != state) {
        if (accessibility[i][state]) {
          current_state[i] = initial[i];
          state_sum += current_state[i];
        }
        else {
          current_state[i] = 0.;
          no_occurrence += initial[i];
        }
      }
    }

    i = 1;

    while ((state_sum > increment) || (i < nb_state - 1)) {

      // change in transition probabilities with the index parameter

      for (j = 0;j < nb_state;j++) {
        if (!homogeneity[j]) {
          transition_update(j , i - 1 , *index_chain);
        }
      }

      // update of the state probabilities

      for (j = 0;j < nb_state;j++) {
        previous_state[j] = current_state[j];
      }

      // computation of the state probabilities and update of
      // the probability of not visiting the selected state

      state_sum = 0.;

      for (j = 0;j < nb_state;j++) {
        if (j != state) {
          if (accessibility[j][state]) {
            current_state[j] = 0.;

            for (k = 0;k < nb_state;k++) {
              if (k != state) {
                current_state[j] += index_chain->transition[k][j] * previous_state[k];
              }
            }

            state_sum += current_state[j];
          }

          else {
            for (k = 0;k < nb_state;k++) {
              if (k != state) {
                no_occurrence += index_chain->transition[k][j] * previous_state[k];
              }
            }
          }
        }
      }

      i++;
    }

    delete index_chain;
    delete [] current_state;
    delete [] previous_state;
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the distribution of the time to the 1st occurrence of a state
 *         for a nonhomogeneous Markov chain.
 *
 *  \param[in] state           state,
 *  \param[in] min_nb_value    minimum number of values,
 *  \param[in] cumul_threshold threshold on the cumulative distribution function.
 */
/*--------------------------------------------------------------*/

void NonhomogeneousMarkov::state_first_occurrence_distribution(int state , int min_nb_value ,
                                                               double cumul_threshold)

{
  int i , j , k;
  double *current_state , *previous_state , *pmass , *pcumul;
  Chain *index_chain;
  Distribution *first_occurrence;


  first_occurrence = process->first_occurrence[state];
  first_occurrence->complement = process->no_occurrence[state];

  pmass = first_occurrence->mass;
  pcumul = first_occurrence->cumul;

  // initialization of the transition probability matrix

  index_chain = new Chain(*this);

  // initialization of the state probabilities

  current_state = new double[nb_state];
  previous_state = new double[nb_state];

  for (i = 0;i < nb_state;i++) {
    if (i != state) {
      current_state[i] = initial[i];
    }
    else {
      *pmass = initial[i];
    }
  }
  *pcumul = *pmass;

  i = 1;

  while (((*pcumul < cumul_threshold - first_occurrence->complement) || (i < min_nb_value)) &&
         (i < first_occurrence->alloc_nb_value)) {

    // change in transition probabilities with the index parameter

    for (j = 0;j < nb_state;j++) {
      if (!homogeneity[j]) {
        transition_update(j , i - 1 , *index_chain);
      }
    }

    // update of the state probabilities

    for (j = 0;j < nb_state;j++) {
      previous_state[j] = current_state[j];
    }

    // computation of the state probabilities and the current probabilty mass

    *++pmass = 0.;

    for (j = 0;j < nb_state;j++) {
      if (j != state) {
        current_state[j] = 0.;

        for (k = 0;k < nb_state;k++) {
          if (k != state) {
            current_state[j] += index_chain->transition[k][j] * previous_state[k];
          }
        }
      }

      else {
        for (k = 0;k < nb_state;k++) {
          if (k != state) {
            *pmass += index_chain->transition[k][j] * previous_state[k];
          }
        }
      }
    }

    // update of the cumulative distribution function

    pcumul++;
    *pcumul = *(pcumul - 1) + *pmass;
    i++;
  }

  first_occurrence->nb_value = i;

# ifdef DEBUG
  if (first_occurrence->complement > 0.) {
    cout << "\n" << SEQ_label[SEQL_STATE_NO_OCCURRENCE] << " " << state << " : "
         << first_occurrence->complement << " | "
         << 1. - first_occurrence->cumul[first_occurrence->nb_value - 1] << endl;
  }
# endif

  first_occurrence->offset_computation();
  first_occurrence->max_computation();
  first_occurrence->mean_computation();
  first_occurrence->variance_computation();

  delete index_chain;
  delete [] current_state;
  delete [] previous_state;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the mixture of the distributions of the number of runs (RUN) or
 *         occurrences (OCCURRENCE) of a state for a sequence length mixing distribution and
 *         a nonhomogeneous Markov chain.
 *
 *  \param[in] state   state,
 *  \param[in] pattern count pattern type.
 */
/*--------------------------------------------------------------*/

void NonhomogeneousMarkov::state_nb_pattern_mixture(int state , count_pattern pattern)

{
  int i , j , k , m;
  int max_length , nb_pattern , index_nb_pattern , increment;
  double sum , **current_state , **previous_state , *cstate , *pstate , *pmass , *lmass;
  Distribution *pdist;
  Chain *index_chain;


  max_length = process->length->nb_value - 1;

  switch (pattern) {
  case RUN :
    pdist = process->nb_run[state];
    nb_pattern = max_length / 2 + 2;
    break;
  case OCCURRENCE :
    pdist = process->nb_occurrence[state];
    nb_pattern = max_length + 1;
    break;
  }

  pmass = pdist->mass;
  for (i = 0;i < pdist->nb_value;i++) {
    *pmass++ = 0.;
  }

  // initialization of the transition probability matrix

  index_chain = new Chain(*this);

  current_state = new double*[nb_state];
  for (i = 0;i < nb_state;i++) {
    current_state[i] = new double[nb_pattern];
  }

  previous_state = new double*[nb_state];
  for (i = 0;i < nb_state;i++) {
    previous_state[i] = new double[nb_pattern];
  }

  lmass = process->length->mass;
  index_nb_pattern = 1;

  for (i = 0;i < max_length;i++) {

    // initialization of the state probabilities for a number of runs or occurrences of
    // the selected state

    if (i == 0) {
      for (j = 0;j < nb_state;j++) {
        if (j == state) {
          current_state[j][0] = 0.;
          current_state[j][1] = initial[j];
        }
        else {
          current_state[j][0] = initial[j];
          current_state[j][1] = 0.;
        }
      }
    }

    else {

      // change in transition probabilities with the index parameter

      for (j = 0;j < nb_state;j++) {
        if (!homogeneity[j]) {
          transition_update(j , i - 1 , *index_chain);
        }
      }

      // update of the state probabilities

      for (j = 0;j < nb_state;j++) {
        for (k = 0;k < index_nb_pattern;k++) {
          previous_state[j][k] = current_state[j][k];
          current_state[j][k] = 0.;
        }
        current_state[j][index_nb_pattern] = 0.;
      }

      for (j = 0;j < nb_state;j++) {

        // computation of the state probabilities for each number of runs or occurrences of
        // the selected state

        for (k = 0;k < nb_state;k++) {
          switch (pattern) {
          case RUN :
            increment = (((k != state) && (j == state)) ? 1 : 0);
            break;
          case OCCURRENCE :
            increment = (j == state ? 1 : 0);
            break;
          }

          cstate = current_state[j];
          pstate = previous_state[k];

          if (increment == 1) {
            cstate++;
          }
          for (m = 0;m < index_nb_pattern;m++) {
            *cstate++ += index_chain->transition[k][j] * *pstate++;
          }
        }
      }
    }

    if ((pattern == OCCURRENCE) || (i % 2 == 0)) {
      index_nb_pattern++;
    }

    // update of the mixture of the distributions of the number of runs or
    // occurrences of the selected state

    if (*++lmass > 0.) {
      pmass = pdist->mass;
      for (j = 0;j < index_nb_pattern;j++) {
        sum = 0.;
        for (k = 0;k < nb_state;k++) {
          sum += current_state[k][j];
        }
        *pmass++ += *lmass * sum;
      }
    }
  }

  pdist->nb_value_computation();
  pdist->offset_computation();
  pdist->cumul_computation();

  pdist->max_computation();
  pdist->mean_computation();
  pdist->variance_computation();

  delete index_chain;

  for (i = 0;i < nb_state;i++) {
    delete [] current_state[i];
  }
  delete [] current_state;

  for (i = 0;i < nb_state;i++) {
    delete [] previous_state[i];
  }
  delete [] previous_state;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the characteristic distributions of a NonhomogeneousMarkov object.
 *
 *  \param[in] length        sequence length,
 *  \param[in] counting_flag flag on the computation of the counting distributions.
 */
/*--------------------------------------------------------------*/

void NonhomogeneousMarkov::characteristic_computation(int length , bool counting_flag)

{
  if (nb_component > 0) {
    int i;
    DiscreteParametric dlength(UNIFORM , length , length , D_DEFAULT , D_DEFAULT);


    if ((!(process->length)) || (dlength != *(process->length))) {
      process->create_characteristic(dlength , homogeneity , counting_flag);

      index_state_distribution();

      for (i = 0;i < nb_state;i++) {
        state_no_occurrence_probability(i);
        state_first_occurrence_distribution(i);

        if (homogeneity[i]) {
          if (stype[i] != ABSORBING) {
            process->sojourn_time[i]->init(NEGATIVE_BINOMIAL , 1 , I_DEFAULT , 1. ,
                                           1. - transition[i][i]);
            process->sojourn_time[i]->computation(1 , OCCUPANCY_THRESHOLD);
            process->sojourn_time[i]->ident = CATEGORICAL;
          }

          else {
            process->absorption[i] = 1.;
            delete process->sojourn_time[i];
            process->sojourn_time[i] = NULL;
          }
        }

        if (counting_flag) {
          state_nb_pattern_mixture(i , RUN);
          state_nb_pattern_mixture(i , OCCURRENCE);
        }
      }
    }
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the characteristic distributions of a NonhomogeneousMarkov object.
 *
 *  \param[in] seq           reference on a NonhomogeneousMarkovData object,
 *  \param[in] counting_flag flag on the computation of the counting distributions,
 *  \param[in] length_flag   flag on the sequence length.
 */
/*--------------------------------------------------------------*/

void NonhomogeneousMarkov::characteristic_computation(const NonhomogeneousMarkovData &seq ,
                                                      bool counting_flag , bool length_flag)

{
  if (nb_component > 0) {
    int i;
    Distribution dlength(*(seq.length_distribution));


    if ((!length_flag) || ((length_flag) && ((!(process->length)) ||
          (dlength != *(process->length))))) {
      process->create_characteristic(dlength , homogeneity , counting_flag);

      index_state_distribution();

      for (i = 0;i < nb_state;i++) {
        state_no_occurrence_probability(i);
        state_first_occurrence_distribution(i , seq.characteristics[0]->first_occurrence[i]->nb_value);

        if (homogeneity[i]) {
          if (stype[i] != ABSORBING) {
            process->sojourn_time[i]->init(NEGATIVE_BINOMIAL , 1 , I_DEFAULT , 1. ,
                                           1. - transition[i][i]);
            process->sojourn_time[i]->computation(seq.characteristics[0]->sojourn_time[i]->nb_value ,
                                                  OCCUPANCY_THRESHOLD);
            process->sojourn_time[i]->ident = CATEGORICAL;
          }

          else {
            process->absorption[i] = 1.;
            delete process->sojourn_time[i];
            process->sojourn_time[i] = NULL;
          }
        }

        if (counting_flag) {
          state_nb_pattern_mixture(i , RUN);
          state_nb_pattern_mixture(i , OCCURRENCE);
        }
      }
    }
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the variation explained by the self-transition probability function.
 *
 *  \param[in] mean mean.
 *
 *  \return         regression square sum.
 */
/*--------------------------------------------------------------*/

double Function::regression_square_sum_computation(double mean) const

{
  int i;
  int *pfrequency;
  double regression_square_sum , diff , *ppoint;


  pfrequency = frequency;
  ppoint = point;
  regression_square_sum = 0.;

  for (i = 0;i <= max_value;i++) {
    if (*pfrequency > 0) {
      diff = *ppoint - mean;
      regression_square_sum += *pfrequency * diff * diff;
    }
    pfrequency++;
    ppoint++;
  }

  return regression_square_sum;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the residuals of a self-transition probability function.
 *
 *  \param[in] self_transition reference on the self-transition probability function for a state.
 */
/*--------------------------------------------------------------*/

void Function::residual_computation(const SelfTransition &self_transition)

{
  int i;
  int *pfrequency , *sfrequency;
  double *presidual , *ppoint , *spoint;


  pfrequency = frequency;
  sfrequency = self_transition.frequency;
  presidual = residual;
  ppoint = point;
  spoint = self_transition.point[0];

  for (i = 0;i <= max_value;i++) {
    *pfrequency = *sfrequency++;
    if (*pfrequency++ > 0) {
      *presidual++ = *spoint - *ppoint;
    }
    else {
      *presidual++ = -D_INF;
    }
    ppoint++;
    spoint++;
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the mean of the residuals of a self-transition probability function.
 *
 *  \return residual mean.
 */
/*--------------------------------------------------------------*/

double Function::residual_mean_computation() const

{
  int i;
  int nb_element , *pfrequency;
  double residual_mean , *presidual;


  pfrequency = frequency;
  presidual = residual;
  nb_element = 0;
  residual_mean = 0.;

  for (i = 0;i <= max_value;i++) {
    if (*pfrequency > 0) {
      nb_element += *pfrequency;
      residual_mean += *pfrequency * *presidual;
    }
    pfrequency++;
    presidual++;
  }
  residual_mean /= nb_element;

  return residual_mean;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the variance of the residuals of a self-transition probability function.
 *
 *  \param[in] residual_mean residual mean.
 *
 *  \return                  residual variance.
 */
/*--------------------------------------------------------------*/

double Function::residual_variance_computation(double residual_mean) const

{
  int i;
  int *pfrequency;
  double residual_variance = D_DEFAULT , diff , *presidual;


  if (residual_df > 0.) {
    pfrequency = frequency;
    presidual = residual;
    residual_variance = 0.;

    for (i = 0;i <= max_value;i++) {
      if (*pfrequency > 0) {
        diff = *presidual - residual_mean;
        residual_variance += *pfrequency * diff * diff;
      }
      pfrequency++;
      presidual++;
    }

    residual_variance /= residual_df;
  }

  return residual_variance;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the sum of squared residuals of a self-transition probability function.
 *
 *  \return residual square sum.
 */
/*--------------------------------------------------------------*/

double Function::residual_square_sum_computation() const

{
  int i;
  int *pfrequency;
  double residual_square_sum , *presidual;


  pfrequency = frequency;
  presidual = residual;
  residual_square_sum = 0.;

  for (i = 0;i <= max_value;i++) {
    if (*pfrequency > 0) {
      residual_square_sum += *pfrequency * *presidual * *presidual;
    }
    pfrequency++;
    presidual++;
  }

  return residual_square_sum;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Estimation of the parameters of the logistic function y = a / (1 + b * exp(-c * x)).
 */
/*--------------------------------------------------------------*/

Function* SelfTransition::logistic_regression() const

{
  int i;
  int iter , nb_element = nb_element_computation() , norm , init_nb_element , *pfrequency;
  double start_proba , denom , residual , residual_square_sum = -D_INF , previous_residual_square_sum ,
         correction[3] , *ppoint;
  Function *function;


  function = new Function(LOGISTIC , length);

  function->regression_df = function->nb_parameter;
  function->residual_df = nb_element - function->nb_parameter;

  // parameter initialization

  init_nb_element = (int)(START_RATIO *  nb_element);
  init_nb_element = MAX(init_nb_element , REGRESSION_NB_ELEMENT / 4);

  pfrequency = frequency;
  ppoint = point[0];
  start_proba = 0.;
  i = 0;

  do {
    start_proba += *pfrequency * *ppoint++;
    i += *pfrequency++;
  }
  while (i < init_nb_element);
  start_proba /= i;

  init_nb_element = (int)(END_RATIO *  nb_element);
  init_nb_element = MAX(init_nb_element , REGRESSION_NB_ELEMENT / 4);

  pfrequency = frequency + length;
  ppoint = point[0] + length;
  function->parameter[0] = 0.;
  i = 0;

  do {
    function->parameter[0] += *--pfrequency * *--ppoint;
    i += *pfrequency;
  }
  while (i < init_nb_element);
  function->parameter[0] /= i;
  function->parameter[1] = function->parameter[0] / start_proba - 1.;

  pfrequency = frequency + 1;
  ppoint = point[0] + 1;
  function->parameter[2] = 0.;
  norm = 0;
  for (i = 1;i < length;i++) {
    if ((*pfrequency > 0) && ((function->parameter[0] / *ppoint - 1.) / function->parameter[1] > 0.)) {
      function->parameter[2] -= *pfrequency * log((function->parameter[0] / *ppoint - 1.) / function->parameter[1]) / i;
      norm += *pfrequency;
    }
    pfrequency++;
    ppoint++;
  }
  function->parameter[2] /= norm;

# ifdef DEBUG
  cout << "\n";
  function->ascii_parameter_print(cout);
  cout << endl;
# endif

  // least-square iterations

  iter = 0;
  do {
    iter++;
    previous_residual_square_sum = residual_square_sum;

    pfrequency = frequency;
    ppoint = point[0];
    residual_square_sum = 0.;

    for (i = 0;i < function->nb_parameter;i++) {
      correction[i] = 0.;
    }
  
    for (i = 0;i < length;i++) {
      if (*pfrequency > 0) {
        denom = 1. + function->parameter[1] * exp(-function->parameter[2] * i);
        residual = *ppoint - function->parameter[0] / denom;
        residual_square_sum += *pfrequency * residual * residual;
        correction[0] += *pfrequency * residual / denom;
        correction[1] -= *pfrequency * residual * function->parameter[0] * exp(-function->parameter[2] * i) /
                         (denom * denom);
        if (i > 0) {
          correction[2] += *pfrequency * residual * function->parameter[0] * function->parameter[1] * i *
                           exp(-function->parameter[2] * i) / (denom * denom);
        }
      }
      pfrequency++;
      ppoint++;
    }
    residual_square_sum /= nb_element;

    function->parameter[0] += GRADIENT_DESCENT_COEFF * correction[0] / nb_element;
    function->parameter[1] += GRADIENT_DESCENT_COEFF * correction[1] / nb_element;
    function->parameter[2] += GRADIENT_DESCENT_COEFF * correction[2] / (nb_element - frequency[0]);

    // application of thresholds on parameters

    if (function->parameter[0] < MIN_PROBABILITY) {
      function->parameter[0] = MIN_PROBABILITY;
    }
    if (function->parameter[0] > 1. - MIN_PROBABILITY) {
      function->parameter[0] = 1. - MIN_PROBABILITY;
    }
    if (function->parameter[0] / (1. + function->parameter[1]) < MIN_PROBABILITY) {
      function->parameter[1] = function->parameter[0] / MIN_PROBABILITY - 1.;
    }
    if (function->parameter[0] / (1. + function->parameter[1]) > 1. - MIN_PROBABILITY) {
      function->parameter[1] = function->parameter[0] / (1. - MIN_PROBABILITY) - 1.;
    }

#   ifdef DEBUG
    if ((iter < 10) || (iter % 10 == 0)) {
      function->ascii_parameter_print(cout);
      cout << "\niteration " << iter << ", " << residual_square_sum << " | "
           << (previous_residual_square_sum - residual_square_sum) / residual_square_sum << endl;
    }
#   endif

  }
  while (((previous_residual_square_sum - residual_square_sum) / residual_square_sum > RESIDUAL_SQUARE_SUM_DIFF) &&
         (iter < REGRESSION_NB_ITER));

  // computation of the logistic function and the residuals

  function->computation();
  function->residual_computation(*this);

  return function;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Estimation of the parameters of the monomolecular function y = a + b * exp(-c * x).
 */
/*--------------------------------------------------------------*/

Function* SelfTransition::monomolecular_regression() const

{
  int i;
  int iter , nb_element = nb_element_computation() , norm , init_nb_element , *pfrequency;
  double start_proba , residual , residual_square_sum = -D_INF , previous_residual_square_sum ,
         correction[3] , *ppoint;
  Function *function;


  function = new Function(MONOMOLECULAR , length);

  function->regression_df = function->nb_parameter;
  function->residual_df = nb_element - function->nb_parameter;

  // parameter initialization

  init_nb_element = (int)(START_RATIO * nb_element);
  init_nb_element = MAX(init_nb_element , REGRESSION_NB_ELEMENT / 4);

  pfrequency = frequency;
  ppoint = point[0];
  start_proba = 0.;
  i = 0;

  do {
    start_proba += *pfrequency * *ppoint++;
    i += *pfrequency++;
  }
  while (i < init_nb_element);
  start_proba /= i;

  init_nb_element = (int)(END_RATIO * nb_element);
  init_nb_element = MAX(init_nb_element , REGRESSION_NB_ELEMENT / 4);

  pfrequency = frequency + length;
  ppoint = point[0] + length;
  function->parameter[0] = 0.;
  i = 0;

  do {
    function->parameter[0] += *--pfrequency * *--ppoint;
    i += *pfrequency;
  }
  while (i < init_nb_element);
  function->parameter[0] /= i;
  function->parameter[1] = start_proba - function->parameter[0];

  pfrequency = frequency + 1;
  ppoint = point[0] + 1;
  function->parameter[2] = 0.;
  norm = 0;
  for (i = 1;i < length;i++) {
    if ((*pfrequency > 0) && ((*ppoint - function->parameter[0]) / function->parameter[1] > 0.)) {
      function->parameter[2] -= *pfrequency * log((*ppoint - function->parameter[0]) / function->parameter[1]) / i;
      norm += *pfrequency;
    }
    pfrequency++;
    ppoint++;
  }
  function->parameter[2] /= norm;

# ifdef DEBUG
  cout << "\n";
  function->ascii_parameter_print(cout);
  cout << endl;
# endif

  // least-square iterations

  iter = 0;
  do {
    iter++;
    previous_residual_square_sum = residual_square_sum;

    pfrequency = frequency;
    ppoint = point[0];
    residual_square_sum = 0.;

    for (i = 0;i < function->nb_parameter;i++) {
      correction[i] = 0.;
    }
  
    for (i = 0;i < length;i++) {
      if (*pfrequency > 0) {
        residual = *ppoint - (function->parameter[0] + function->parameter[1] *
                   exp(-function->parameter[2] * i));
        residual_square_sum += *pfrequency * residual * residual;
        correction[0] += *pfrequency * residual;
        correction[1] += *pfrequency * residual * exp(-function->parameter[2] * i);
        if (i > 0) {
          correction[2] -= *pfrequency * residual * function->parameter[1] * i *
                           exp(-function->parameter[2] * i);
        }
      }
      pfrequency++;
      ppoint++;
    }
    residual_square_sum /= nb_element;

    function->parameter[0] += GRADIENT_DESCENT_COEFF * correction[0] / nb_element;
    function->parameter[1] += GRADIENT_DESCENT_COEFF * correction[1] / nb_element;
    function->parameter[2] += GRADIENT_DESCENT_COEFF * correction[2] / (nb_element - frequency[0]);

    // application of thresholds on parameters

    if (function->parameter[0] < MIN_PROBABILITY) {
      function->parameter[0] = MIN_PROBABILITY;
    }
    if (function->parameter[0] > 1. - MIN_PROBABILITY) {
      function->parameter[0] = 1. - MIN_PROBABILITY;
    }
    if (function->parameter[0] + function->parameter[1] < MIN_PROBABILITY) {
      function->parameter[1] = MIN_PROBABILITY - function->parameter[0];
    }
    if (function->parameter[0] + function->parameter[1] > 1. - MIN_PROBABILITY) {
      function->parameter[1] = 1. - MIN_PROBABILITY - function->parameter[0];
    }

#   ifdef DEBUG
    if ((iter < 10) || (iter % 10 == 0)) {
      function->ascii_parameter_print(cout);
      cout << "\niteration " << iter << ", " << residual_square_sum << " | "
           << (previous_residual_square_sum - residual_square_sum) / residual_square_sum << endl;
    }
#   endif

  }
  while (((previous_residual_square_sum - residual_square_sum) / residual_square_sum > RESIDUAL_SQUARE_SUM_DIFF) &&
         (iter < REGRESSION_NB_ITER));

  // computation of the monomolecular function and the residuals

  function->computation();
  function->residual_computation(*this);

  return function;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the log-likelihood of a nonhomogeneous Markov chain for sequences.
 *
 *  \param[in] seq   reference on a MarkovianSequences object,
 *  \param[in] index sequence index.
 *
 *  \return          log-likelihood.
 */
/*--------------------------------------------------------------*/

double NonhomogeneousMarkov::likelihood_computation(const MarkovianSequences &seq , int index) const

{
  int i , j , k;
  int *pstate;
  double likelihood = 0. , proba;
  Chain *index_chain;


  // checking of the compatibility of the model with the data

  if (seq.nb_variable == 1) {
    if ((seq.marginal_distribution[0]) &&
        (nb_state < seq.marginal_distribution[0]->nb_value)) {
      likelihood = D_INF;
    }
  }

  else {
    likelihood = D_INF;
  }

  if (likelihood != D_INF) {
    index_chain = new Chain(*this);

    for (i = 0;i < seq.nb_sequence;i++) {
      if ((index == I_DEFAULT) || (index == i)) {

        // initialization of the transition probability matrix

        if (i > 0) {
          index_chain->parameter_copy(*this);
        }

        pstate = seq.int_sequence[i][0];

        proba = initial[*pstate];
        if (proba > 0.) {
          likelihood += log(proba);
        }
        else {
          likelihood = D_INF;
          break;
        }

        for (j = 1;j < seq.length[i];j++) {

          // change in transition probabilities with the index parameter

          for (k = 0;k < nb_state;k++) {
            if (!homogeneity[k]) {
              transition_update(k , j - 1 , *index_chain);
            }
          }

          proba = index_chain->transition[*pstate][*(pstate + 1)];
          pstate++;

          if (proba > 0.) {
            likelihood += log(proba);
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

    delete index_chain;
  }

  return likelihood;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Construction of the initial state and transition counts.
 */
/*--------------------------------------------------------------*/

void NonhomogeneousMarkovData::build_transition_count()

{
  chain_data = new ChainData(ORDINARY , marginal_distribution[0]->nb_value ,
                             marginal_distribution[0]->nb_value);
  transition_count_computation(*chain_data);
}


/*--------------------------------------------------------------*/
/**
 *  \brief Estimation of a nonhomogeneous Markov chain.
 *
 *  \param[in] error         reference on a StatError object,
 *  \param[in] ident         identifiers of the self-transition probability functions,
 *  \param[in] counting_flag flag on the computation of the counting distributions.
 *
 *  \return                  NonhomogeneousMarkov object.
 */
/*--------------------------------------------------------------*/

NonhomogeneousMarkov* MarkovianSequences::nonhomogeneous_markov_estimation(StatError &error , parametric_function *ident ,
                                                                           bool counting_flag) const

{
  bool status = true;
  int i;
  NonhomogeneousMarkov *markov;
  NonhomogeneousMarkovData *seq;


  markov = NULL;
  error.init();

  if (nb_variable > 1) {
    status = false;
    error.correction_update(STAT_error[STATR_NB_VARIABLE] , 1);
  }
  if ((marginal_distribution[0]->nb_value < 2) ||
      (marginal_distribution[0]->nb_value > NB_STATE)) {
    status = false;
    error.update(SEQ_error[SEQR_NB_STATE]);
  }

  if (status) {
    markov = new NonhomogeneousMarkov(marginal_distribution[0]->nb_value , ident);
    markov->markov_data = new NonhomogeneousMarkovData(*this);

    seq = markov->markov_data;
    seq->state_variable_init();
    seq->build_transition_count();

    // estimation of the Markov chain parameters

    seq->chain_data->estimation(*markov);

    // estimation of the self-transition probability functions

    seq->self_transition_computation(markov->homogeneity);

    for (i = 0;i < markov->nb_state;i++) {
      if (!(markov->homogeneity[i])) {
 
#       ifdef DEBUG
        cout << *(seq->self_transition[i]);
#       endif

        if (seq->self_transition[i]->nb_element_computation() >= REGRESSION_NB_ELEMENT) {
          switch (ident[i]) {
          case LOGISTIC :
            markov->self_transition[i] = seq->self_transition[i]->logistic_regression();
            break;
          case MONOMOLECULAR :
            markov->self_transition[i] = seq->self_transition[i]->monomolecular_regression();
            break;
          }
        }

        else {
          markov->homogeneity[i] = true;
        }
      }
    }

    for (i = 0;i < markov->nb_state;i++) {
      if (!(markov->homogeneity[i])) {
        break;
      }
    }

    if (i == markov->nb_state) {
      delete [] markov->self_transition;
      markov->self_transition = NULL;
    }

    // computation of the log-likelihood and the characteristic distributions of the model

    seq->likelihood = markov->likelihood_computation(*seq);

    if (seq->likelihood == D_INF) {
      delete markov;
      markov = NULL;
      error.update(STAT_error[STATR_ESTIMATION_FAILURE]);
    }

    else {
      markov->component_computation();
      markov->characteristic_computation(*seq , counting_flag , false);
    }
  }

  return markov;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Simulation using a nonhomogeneous Markov chain.
 *
 *  \param[in] error               reference on a StatError object,
 *  \param[in] length_distribution sequence length frequency distribution,
 *  \param[in] counting_flag       flag on the computation of the counting distributions.
 *
 *  \return                        NonhomogeneousMarkovData object.
 */
/*--------------------------------------------------------------*/

NonhomogeneousMarkovData* NonhomogeneousMarkov::simulation(StatError &error ,
                                                           const FrequencyDistribution &length_distribution ,
                                                           bool counting_flag) const

{
  bool status = true;
  int i , j , k;
  int cumul_length , *pstate;
  Chain *index_chain;
  NonhomogeneousMarkov *markov;
  NonhomogeneousMarkovData *seq;


  seq = NULL;
  error.init();

  if ((length_distribution.nb_element < 1) || (length_distribution.nb_element > NB_SEQUENCE)) {
    status = false;
    error.update(SEQ_error[SEQR_NB_SEQUENCE]);
  }
  if (length_distribution.offset < 2) {
    status = false;
    error.update(SEQ_error[SEQR_SHORT_SEQUENCE_LENGTH]);
  }
  if (length_distribution.nb_value - 1 > MAX_LENGTH) {
    status = false;
    error.update(SEQ_error[SEQR_LONG_SEQUENCE_LENGTH]);
  }

  if (status) {
    cumul_length = 0;
    for (i = length_distribution.offset;i < length_distribution.nb_value;i++) {
      cumul_length += i * length_distribution.frequency[i];
    }

    if (cumul_length > CUMUL_LENGTH) {
      status = false;
      error.update(SEQ_error[SEQR_CUMUL_SEQUENCE_LENGTH]);
    }
  }

  if (status) {

    // initializations

    seq = new NonhomogeneousMarkovData(length_distribution);
    seq->type[0] = STATE;

    seq->markov = new NonhomogeneousMarkov(*this , false);

    markov = seq->markov;
    markov->create_cumul();
    markov->cumul_computation();

    for (i = 0;i < markov->nb_state;i++) {
      if (!(markov->homogeneity[i])) {
        if (markov->self_transition[i]->max_value < seq->max_length - 1) {
          delete [] markov->self_transition[i]->point;
          markov->self_transition[i]->max_value = seq->max_length - 1;
          markov->self_transition[i]->point = new double[markov->self_transition[i]->max_value + 1];
          markov->self_transition[i]->computation();
        }
      }
    }

    index_chain = new Chain(*markov);

    for (i = 0;i < seq->nb_sequence;i++) {

      // initialization of the transition probability matrix

      index_chain->parameter_copy(*this);

      pstate = seq->int_sequence[i][0];
      *pstate = cumul_method(markov->nb_state , markov->cumul_initial);
    
      for (j = 1;j < seq->length[i];j++) {

        // change in transition probabilities with the index parameter

        for (k = 0;k < markov->nb_state;k++) {
          if (!markov->homogeneity[k]) {
            transition_update(k , j - 1 , *index_chain);
          }
        }

        *(pstate + 1) = cumul_method(markov->nb_state , index_chain->cumul_transition[*pstate]);
        pstate++;
      }
    }

    markov->remove_cumul();
    delete index_chain;

    // computation of the characteristics of the generated sequences

    for (i = 0;i < seq->nb_variable;i++) {
      seq->max_value_computation(i);
      seq->build_marginal_frequency_distribution(i);
    }

    seq->self_transition_computation(markov->homogeneity);
    seq->build_transition_count();
    seq->build_characteristic();

    for (i = 0;i < markov->nb_state;i++) {
      if (!(markov->homogeneity[i])) {
        markov->self_transition[i]->regression_df = markov->self_transition[i]->nb_parameter;
        markov->self_transition[i]->residual_df = seq->self_transition[i]->nb_element_computation() -
                                                  markov->self_transition[i]->nb_parameter;

        delete [] markov->self_transition[i]->residual;
        delete [] markov->self_transition[i]->frequency;
        markov->self_transition[i]->residual = new double[markov->self_transition[i]->max_value + 1];
        markov->self_transition[i]->frequency = new int[markov->self_transition[i]->max_value + 1];

        markov->self_transition[i]->residual_computation(*(seq->self_transition[i]));
      }
    }

    markov->characteristic_computation(*seq , counting_flag);

    // computation of the log-likelihood of the model for the generated sequences

    seq->likelihood = markov->likelihood_computation(*seq);
  }

  return seq;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Simulation using a nonhomogeneous Markov chain.
 *
 *  \param[in] error         reference on a StatError object,
 *  \param[in] nb_sequence   number of sequences,
 *  \param[in] length        sequence length,
 *  \param[in] counting_flag flag on the computation of the counting distributions.
 *
 *  \return                  NonhomogeneousMarkovData object.
 */
/*--------------------------------------------------------------*/

NonhomogeneousMarkovData* NonhomogeneousMarkov::simulation(StatError &error , int nb_sequence ,
                                                           int length , bool counting_flag) const

{
  bool status = true;
  NonhomogeneousMarkovData *seq;


  seq = NULL;
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
    FrequencyDistribution length_distribution(length + 1);

    length_distribution.nb_element = nb_sequence;
    length_distribution.offset = length;
    length_distribution.max = nb_sequence;
    length_distribution.mean = length;
    length_distribution.variance = 0.;
    length_distribution.frequency[length] = nb_sequence;

    seq = simulation(error , length_distribution , counting_flag);
  }

  return seq;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Simulation using a nonhomogeneous Markov chain.
 *
 *  \param[in] error         reference on a StatError object,
 *  \param[in] nb_sequence   number of sequences,
 *  \param[in] iseq          reference on a MarkovianSequences object,
 *  \param[in] counting_flag flag on the computation of the counting distributions.
 *
 *  \return                  NonhomogeneousMarkovData object.
 */
/*--------------------------------------------------------------*/

NonhomogeneousMarkovData* NonhomogeneousMarkov::simulation(StatError &error , int nb_sequence ,
                                                           const MarkovianSequences &iseq ,
                                                           bool counting_flag) const

{
  FrequencyDistribution *length_distribution;
  NonhomogeneousMarkovData *seq;


  error.init();

  if ((nb_sequence < 1) || (nb_sequence > NB_SEQUENCE)) {
    seq = NULL;
    error.update(SEQ_error[SEQR_NB_SEQUENCE]);
  }

  else {
    length_distribution = iseq.length_distribution->frequency_scale(nb_sequence);

    seq = simulation(error , *length_distribution , counting_flag);
    delete length_distribution;
  }

  return seq;
}


};  // namespace sequence_analysis
