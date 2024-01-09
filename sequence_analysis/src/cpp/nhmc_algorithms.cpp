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
 *       $Id: nhmc_algorithms.cpp 18059 2015-04-23 10:47:57Z guedon $
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

#include "stat_tool/stat_tools.h"
#include "stat_tool/regression.h"
#include "stat_tool/curves.h"
#include "stat_tool/distribution.h"
#include "stat_tool/markovian.h"
#include "stat_tool/vectors.h"
#include "stat_tool/distance_matrix.h"
#include "stat_tool/stat_label.h"

#include "sequences.h"
#include "nonhomogeneous_markov.h"
#include "sequence_label.h"

using namespace std;
using namespace stat_tool;


namespace sequence_analysis {



/*--------------------------------------------------------------*
 *
 *  Mise a jour des probabilites de transition a partir d'un etat donne
 *  pour une chaine de Markov non-homogene.
 *
 *  arguments : etat, index, reference sur les probabilites de transition courante.
 *
 *--------------------------------------------------------------*/

void NonhomogeneousMarkov::transition_update(int state , int index , Chain &index_chain) const


{
  register int i;
  double scale , *pparam;


  pparam = self_transition[state]->parameter;

  // mise a jour de la probabilite de rester dans un etat

  switch (self_transition[state]->ident) {
  case STAT_MONOMOLECULAR :
    index_chain.transition[state][state] = pparam[0] + pparam[1] * exp(-pparam[2] * index);
    break;
  case STAT_LOGISTIC :
    index_chain.transition[state][state] = pparam[0] / (1. + pparam[1] * exp(-pparam[2] * index));
    break;
  }

  // mise a jour des probabilites de passage

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


/*--------------------------------------------------------------*
 *
 *  Calcul des probabilites de chaque etat en fonction de l'index
 *  pour une chaine de Markov non-homogene.
 *
 *--------------------------------------------------------------*/

void NonhomogeneousMarkov::index_state_distribution()

{
  register int i , j , k;
  Curves *index_state;
  Chain *index_chain;


  index_state = process->index_value;

  // initialisation de la matrice des probabilites de transition courante

  index_chain = new Chain(*this);

  // initialisation des probabilites des etats

  for (i = 0;i < nb_state;i++) {
    index_state->point[i][0] = initial[i];
  }

  // calcul des probabilites de chaque etat en fonction de l'index

  for (i = 1;i < index_state->length;i++) {

    // prise en compte de l'evolution des probabilites de transition

    for (j = 0;j < nb_state;j++) {
      if (!homogeneity[j]) {
        transition_update(j , i - 1 , *index_chain);
      }
    }

    // calcul des probabilites des etats

    for (j = 0;j < nb_state;j++) {
      index_state->point[j][i] = 0.;
      for (k = 0;k < nb_state;k++) {
        index_state->point[j][i] += index_chain->transition[k][j] * index_state->point[k][i - 1];
      }
    }
  }

  delete index_chain;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la probabilite de ne pas observer un etat
 *  d'une chaine de Markov non-homogene.
 *
 *  arguments : etat, seuil sur la somme des probabilites des etats.
 *
 *--------------------------------------------------------------*/

void NonhomogeneousMarkov::state_no_occurrence_probability(int state , double increment)

{
  register int i;

  for (i = 0;i < nb_state;i++) {
    if ((i != state) && (!accessibility[i][state])) {
      break;
    }
  }

  if (i < nb_state) {
    register int j , k;
    double state_sum , *current_state , *previous_state ,
           &no_occurrence = process->no_occurrence[state];
    Chain *index_chain;


    // initialisation de la matrice des probabilites de transition courante

    index_chain = new Chain(*this);

    // initialisation des probabilites des etats

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

      // prise en compte de l'evolution des probabilites de transition

      for (j = 0;j < nb_state;j++) {
        if (!homogeneity[j]) {
          transition_update(j , i - 1 , *index_chain);
        }
      }

      // mise a jour des probabilites des etats

      for (j = 0;j < nb_state;j++) {
        previous_state[j] = current_state[j];
      }

      // calcul des probabilites des etats et mise a jour
      // de la probabilite de ne pas observer l'etat

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


/*--------------------------------------------------------------*
 *
 *  Calcul de la loi du temps avant la premiere occurrence d'un etat
 *  pour une chaine de Markov non-homogene.
 *
 *  arguments : etat, nombre minimum de valeurs,
 *              seuil sur la fonction de repartition.
 *
 *--------------------------------------------------------------*/

void NonhomogeneousMarkov::state_first_occurrence_distribution(int state , int min_nb_value ,
                                                               double cumul_threshold)

{
  register int i , j , k;
  double *current_state , *previous_state , *pmass , *pcumul;
  Chain *index_chain;
  Distribution *first_occurrence;


  first_occurrence = process->first_occurrence[state];
  first_occurrence->complement = process->no_occurrence[state];

  pmass = first_occurrence->mass;
  pcumul = first_occurrence->cumul;

  // initialisation de la matrice des probabilites de transition courante

  index_chain = new Chain(*this);

  // initialisation des probabilites des etats

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

    // prise en compte de l'evolution des probabilites de transition

    for (j = 0;j < nb_state;j++) {
      if (!homogeneity[j]) {
        transition_update(j , i - 1 , *index_chain);
      }
    }

    // mise a jour des probabilites des etats

    for (j = 0;j < nb_state;j++) {
      previous_state[j] = current_state[j];
    }

    // calcul des probabilites des etats et de la valeur courante

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

    // mise a jour de la fonction de repartition

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


/*--------------------------------------------------------------*
 *
 *  Calcul d'un melange de lois du nombre de series d'un etat ('r')
 *  ou du nombre d'occurrences d'un etat ('o') d'une chaine de Markov non-homogene
 *  pour une distribution des longueurs de sequences donnee.
 *
 *  arguments : etat, type de forme.
 *
 *--------------------------------------------------------------*/

void NonhomogeneousMarkov::state_nb_pattern_mixture(int state , char pattern)

{
  register int i , j , k , m;
  int max_length , nb_pattern , index_nb_pattern , increment;
  double sum , **current_state , **previous_state , *cstate , *pstate , *pmass , *lmass;
  Distribution *pdist;
  Chain *index_chain;


  max_length = process->length->nb_value - 1;

  switch (pattern) {
  case 'r' :
    pdist = process->nb_run[state];
    nb_pattern = max_length / 2 + 2;
    break;
  case 'o' :
    pdist = process->nb_occurrence[state];
    nb_pattern = max_length + 1;
    break;
  }

  pmass = pdist->mass;
  for (i = 0;i < pdist->nb_value;i++) {
    *pmass++ = 0.;
  }

  // initialisation de la matrice des probabilites de transition courante

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

    // initialisation des probabilites des etats pour un nombre de formes donne
    // de l'etat selectionne

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

      // prise en compte de l'evolution des probabilites de transition

      for (j = 0;j < nb_state;j++) {
        if (!homogeneity[j]) {
          transition_update(j , i - 1 , *index_chain);
        }
      }

      // mise a jour des probabilites des etats

      for (j = 0;j < nb_state;j++) {
        for (k = 0;k < index_nb_pattern;k++) {
          previous_state[j][k] = current_state[j][k];
          current_state[j][k] = 0.;
        }
        current_state[j][index_nb_pattern] = 0.;
      }

      for (j = 0;j < nb_state;j++) {

        // calcul des probabilites des etats pour chaque nombre de formes
        // de l'etat selectionne

        for (k = 0;k < nb_state;k++) {
          switch (pattern) {
          case 'r' :
            increment = (((k != state) && (j == state)) ? 1 : 0);
            break;
          case 'o' :
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

    if ((pattern == 'o') || (i % 2 == 0)) {
      index_nb_pattern++;
    }

    // mise a jour du melange de lois du nombre de formes de l'etat selectionne

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


/*--------------------------------------------------------------*
 *
 *  Calcul des lois caracteristiques d'un objet NonhomogeneousMarkov.
 *
 *  arguments : longueur des sequences, flag sur le calcul des lois de comptage.
 *
 *--------------------------------------------------------------*/

void NonhomogeneousMarkov::characteristic_computation(int length , bool counting_flag)

{
  if (nb_component > 0) {
    register int i;
    DiscreteParametric dlength(UNIFORM , length , length , D_DEFAULT , D_DEFAULT);


    // calcul des lois caracteristiques au niveau etat

    if ((!(process->length)) || (dlength != *(process->length))) {
      process->create_characteristic(dlength , homogeneity , counting_flag);

      index_state_distribution();

      for (i = 0;i < nb_state;i++) {
        state_no_occurrence_probability(i);
        state_first_occurrence_distribution(i);

        if (homogeneity[i]) {
          if (state_type[i] != 'a') {
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
          state_nb_pattern_mixture(i , 'r');
          state_nb_pattern_mixture(i , 'o');
        }
      }
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Calcul des lois caracteristiques d'un objet NonhomogeneousMarkov.
 *
 *  arguments : reference sur un objet NonhomogeneousMarkovData,
 *              flag sur le calcul des lois de comptage,
 *              flag pour tenir compte des longueurs.
 *
 *--------------------------------------------------------------*/

void NonhomogeneousMarkov::characteristic_computation(const NonhomogeneousMarkovData &seq ,
                                                      bool counting_flag , bool length_flag)

{
  if (nb_component > 0) {
    register int i;
    Distribution dlength(*(seq.length_distribution));


    // calcul des lois caracteristiques au niveau etat

    if ((!length_flag) || ((length_flag) && ((!(process->length)) ||
          (dlength != *(process->length))))) {
      process->create_characteristic(dlength , homogeneity , counting_flag);

      index_state_distribution();

      for (i = 0;i < nb_state;i++) {
        state_no_occurrence_probability(i);
        state_first_occurrence_distribution(i , seq.characteristics[0]->first_occurrence[i]->nb_value);

        if (homogeneity[i]) {
          if (state_type[i] != 'a') {
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
          state_nb_pattern_mixture(i , 'r');
          state_nb_pattern_mixture(i , 'o');
        }
      }
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la variation expliquee par le modele.
 *
 *  argument : moyenne.
 *
 *--------------------------------------------------------------*/

double Function::regression_square_sum_computation(double mean) const

{
  register int i;
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


/*--------------------------------------------------------------*
 *
 *  Calcul des residus.
 *
 *  argument : reference sur les probabilites de rester dans un etat.
 *
 *--------------------------------------------------------------*/

void Function::residual_computation(const SelfTransition &self_transition)

{
  register int i;
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


/*--------------------------------------------------------------*
 *
 *  Calcul de la moyenne des residus.
 *
 *--------------------------------------------------------------*/

double Function::residual_mean_computation() const

{
  register int i;
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


/*--------------------------------------------------------------*
 *
 *  Calcul de la variance des residus.
 *
 *  argument : moyenne.
 *
 *--------------------------------------------------------------*/

double Function::residual_variance_computation(double residual_mean) const

{
  register int i;
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


/*--------------------------------------------------------------*
 *
 *  Calcul de la somme des carres des residus.
 *
 *--------------------------------------------------------------*/

double Function::residual_square_sum_computation() const

{
  register int i;
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


/*--------------------------------------------------------------*
 *
 *  Estimation des parametres de la fonction y = a + b * exp(-c * x)
 *  par regression.
 *
 *--------------------------------------------------------------*/

Function* SelfTransition::monomolecular_regression() const

{
  register int i;
  int iter , nb_element = nb_element_computation() , norm , init_nb_element , *pfrequency;
  double start_proba , residual , residual_square_sum = -D_INF , previous_residual_square_sum ,
         correction[3] , *ppoint;
  Function *function;


  function = new Function(STAT_MONOMOLECULAR , length);

  function->regression_df = function->nb_parameter;
  function->residual_df = nb_element - function->nb_parameter;

  // initialisations des parametres

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

  // iterations (gradient sur les moindres-carres)

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

    // applications de seuils sur les parametres

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

  // calcul de la fonction et des residus

  function->computation();
  function->residual_computation(*this);

  return function;
}


/*--------------------------------------------------------------*
 *
 *  Estimation des parametres de la fonction y = a / (1 + b * exp(-c * x))
 *  par regression.
 *
 *--------------------------------------------------------------*/

Function* SelfTransition::logistic_regression() const

{
  register int i;
  int iter , nb_element = nb_element_computation() , norm , init_nb_element , *pfrequency;
  double start_proba , denom , residual , residual_square_sum = -D_INF , previous_residual_square_sum ,
         correction[3] , *ppoint;
  Function *function;


  function = new Function(STAT_LOGISTIC , length);

  function->regression_df = function->nb_parameter;
  function->residual_df = nb_element - function->nb_parameter;

  // initialisations des parametres

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

  // iterations (gradient sur les moindres-carres)

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

    // applications de seuils sur les parametres

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

  // calcul de la fonction et des residus

  function->computation();
  function->residual_computation(*this);

  return function;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la vraisemblance de sequences pour une chaine de Markov non-homogene.
 *
 *  arguments : reference sur un objet MarkovianSequences, indice de la sequence.
 *
 *--------------------------------------------------------------*/

double NonhomogeneousMarkov::likelihood_computation(const MarkovianSequences &seq , int index) const

{
  register int i , j , k;
  int *pstate;
  double likelihood = 0. , proba;
  Chain *index_chain;


  // verification de la compatibilite entre le modele et les donnees

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

        // initialisation de la matrice des probabilites de transition courante

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

          // prise en compte de l'evolution des probabilites de transition

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


/*--------------------------------------------------------------*
 *
 *  Construction des comptages des etats initiaux et des transitions.
 *
 *--------------------------------------------------------------*/

void NonhomogeneousMarkovData::build_transition_count()

{
  chain_data = new ChainData('o' , marginal_distribution[0]->nb_value ,
                             marginal_distribution[0]->nb_value);
  transition_count_computation(*chain_data);
}


/*--------------------------------------------------------------*
 *
 *  Estimation des parametres d'une chaine de Markov non-homogene
 *  a partir d'un echantillon de sequences.
 *
 *  arguments : reference sur un objet StatError,
 *              identificateurs evolution des probabilites de rester dans un etat,
 *              flag sur le calcul des lois de comptage.
 *
 *--------------------------------------------------------------*/

NonhomogeneousMarkov* MarkovianSequences::nonhomogeneous_markov_estimation(StatError &error , int *ident ,
                                                                           bool counting_flag) const

{
  bool status = true;
  register int i;
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

    // estimation des parametres de la chaine de Markov

    seq->chain_data->estimation(*markov);

    // estimation des fonctions d'evolution des probabilites de rester
    // dans un etat

    seq->self_transition_computation(markov->homogeneity);

    for (i = 0;i < markov->nb_state;i++) {
      if (!(markov->homogeneity[i])) {
 
#       ifdef DEBUG
        cout << *(seq->self_transition[i]);
#       endif

        if (seq->self_transition[i]->nb_element_computation() >= REGRESSION_NB_ELEMENT) {
          switch (ident[i]) {
          case STAT_MONOMOLECULAR :
            markov->self_transition[i] = seq->self_transition[i]->monomolecular_regression();
            break;
          case STAT_LOGISTIC :
            markov->self_transition[i] = seq->self_transition[i]->logistic_regression();
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

    // calcul de la vraisemblance et des lois caracteristiques du modele

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


/*--------------------------------------------------------------*
 *
 *  Simulation par une chaine de Markov non-homogene.
 *
 *  arguments : reference sur un objet StatError,
 *              loi empirique des longueurs des sequences,
 *              flag sur le calcul des lois de comptage.
 *
 *--------------------------------------------------------------*/

NonhomogeneousMarkovData* NonhomogeneousMarkov::simulation(StatError &error ,
                                                           const FrequencyDistribution &length_distribution ,
                                                           bool counting_flag) const

{
  bool status = true;
  register int i , j , k;
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

    // initialisations

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

      // initialisation de la matrice des probabilites de transition courante

      index_chain->parameter_copy(*this);

      pstate = seq->int_sequence[i][0];
      *pstate = cumul_method(markov->nb_state , markov->cumul_initial);
    
      for (j = 1;j < seq->length[i];j++) {

        // prise en compte de l'evolution des probabilites de transition

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

    // extraction des caracteristiques des sequences simulees

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

    // calcul de la vraisemblance

    seq->likelihood = markov->likelihood_computation(*seq);
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Simulation par une chaine de Markov non-homogene.
 *
 *  arguments : reference sur un objet StatError,
 *              nombre et longueur des sequences,
 *              flag sur le calcul des lois de comptage.
 *
 *--------------------------------------------------------------*/

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


/*--------------------------------------------------------------*
 *
 *  Simulation par une chaine de Markov non-homogene.
 *
 *  arguments : reference sur un objet StatError, nombre de sequences,
 *              reference sur un objet MarkovianSequences,
 *              flag sur le calcul des lois de comptage.
 *
 *--------------------------------------------------------------*/

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
