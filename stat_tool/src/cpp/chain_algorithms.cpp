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
 *       $Id: chain_algorithms.cpp 17974 2015-04-23 06:33:49Z guedon $
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
#include <sstream>

#include "stat_tools.h"
#include "markovian.h"
#include "stat_label.h"

using namespace std;


namespace stat_tool {



/*--------------------------------------------------------------*
 *
 *  Initialisation des parametres d'une chaine de Markov
 *
 *  arguments : flag sur la nature de la chaine de Markov,
 *              probabilite de rester dans un etat.
 *
 *--------------------------------------------------------------*/

void Chain::init(bool left_right , double self_transition)

{
  register int i , j;


  accessibility = new bool*[nb_state];
  for (i = 0;i < nb_state;i++) {
    accessibility[i] = new bool[nb_state];
  }

  state_type = new char[nb_state];

  switch (left_right) {

  // cas chaine de Markov ou toutes les transitions sont possibles

  case false : {
    nb_component = 1;
    component_nb_state = new int[nb_component];
    component_nb_state[0] = nb_state;
    component = new int*[nb_component];
    component[0] = new int[component_nb_state[0]];

    for (i = 0;i < nb_state;i++) {
      for (j = 0;j < nb_state;j++) {
        accessibility[i][j] = true;
      }

      component[0][i] = i;
      state_type[i] = 'r';
    }

    for (i = 0;i < nb_state;i++) {
      initial[i] = 1. / (double)nb_state;
    }

    for (i = 0;i < nb_state;i++) {
      for (j = 0;j < i;j++) {
        transition[i][j] = (1. - self_transition) / (nb_state - 1);
      }
      transition[i][i] = self_transition;
      for (j = i + 1;j < nb_state;j++) {
        transition[i][j] = (1. - self_transition) / (nb_state - 1);
      }
    }

    break;
  }

  // cas chaine de Markov gauche-droite

  case true : {
    nb_component = nb_state;
    component_nb_state = new int[nb_component];
    component = new int*[nb_component];

    for (i = 0;i < nb_state;i++) {
      for (j = 0;j <= i;j++) {
        accessibility[i][j] = false;
      }
      for (j = i + 1;j < nb_state;j++) {
        accessibility[i][j] = true;
      }

      component_nb_state[i] = 1;
      component[i] = new int[component_nb_state[i]];
      component[i][0] = i;

      if (i < nb_state - 1) {
        state_type[i] = 't';
      }
      else {
        state_type[i] = 'a';
      }
    }

    for (i = 0;i < nb_state - 1;i++) {
      initial[i] = 1. / (double)(nb_state - 1);
    }
    initial[nb_state - 1] = 0.;

    for (i = 0;i < nb_state;i++) {
      for (j = 0;j < i;j++) {
        transition[i][j] = 0.;
      }

      if (i < nb_state - 1) {
        transition[i][i] = self_transition;
        for (j = i + 1;j < nb_state;j++) {
          transition[i][j] = (1. - self_transition) / (nb_state - (i + 1));
        }
      }

      else {
        transition[i][i] = 1.;
      }
    }

    break;
  }
  }
}


/*--------------------------------------------------------------*
 *
 *  Construction de la matrice des transitions possibles entre etats
 *  d'une chaine de Markov d'ordre 1.
 *
 *--------------------------------------------------------------*/

bool** Chain::logic_transition_computation() const

{
  bool **logic_transition;
  register int i , j;


  // calcul de la matrice des transitions possibles entre etats

  logic_transition = new bool*[nb_state];

  for (i = 0;i < nb_state;i++) {
    logic_transition[i] = new bool[nb_state];

    for (j = 0;j < nb_state;j++) {
      if (j == i) {
        logic_transition[i][j] = false;
      }
      else {
        logic_transition[i][j] = (transition[i][j] == 0. ? false : true);
      }
    }
  }

  return logic_transition;
}


/*--------------------------------------------------------------*
 *
 *  Recherche d'une composante connexe dans une chaine de Markov.
 *
 *  arguments : reference sur un objet StatError,
 *              matrice des transitions possibles entre etats.
 *
 *--------------------------------------------------------------*/

bool Chain::connex_component_research(StatError &error , bool **ilogic_transition) const

{
  bool status = true , **logic_transition;
  register int i , j;
  int state , test_state , nb_used_state , *state_transition , *used_transition ,
      *predecessor , *state_order;


  if (ilogic_transition) {
    logic_transition = ilogic_transition;
  }

  else {

    // construction de la matrice des transitions possibles entre etats

    logic_transition = logic_transition_computation();
  }

  // test etats initiaux tels que probabilite initiale > 0

  if (type == 'o') {
    for (i = 0;i < nb_state;i++) {
      for (j = 0;j < nb_state;j++) {
        if (logic_transition[j][i]) {
          break;
        }
      }

      if ((j == nb_state) && (initial[i] == 0.)) {
        status = false;
        ostringstream error_message;
        error_message << STAT_label[STATL_STATE] << " " << i << ": "
                      << STAT_error[STATR_NULL_INITIAL_PROBABILITY];
        error.update((error_message.str()).c_str());
      }
    }
  }

  // calcul de connexite (algorithme de Tarjan) : initialisations

  state_transition = new int[nb_state];

  for (i = 0;i < nb_state;i++) {
    state_transition[i] = 0;
    for (j = 0;j < nb_state;j++) {
      if ((logic_transition[i][j]) || (logic_transition[j][i])) {
        state_transition[i]++;
      }
    }
  }

  used_transition = new int[nb_state];
  predecessor = new int[nb_state];
  state_order = new int[nb_state];

  for (i = 0;i < nb_state;i++) {
    used_transition[i] = 0;
    predecessor[i] = I_DEFAULT;
    state_order[i] = I_DEFAULT;
  }

  state = 0;
  predecessor[state] = state;
  state_order[state] = 0;
  nb_used_state = 1;

  do {
    if (used_transition[state] == state_transition[state]) {
      state = predecessor[state];
    }

    else {
      used_transition[state]++;
      i = 0;

      // recherche du prochain etat

      for (j = 0;j < nb_state;j++) {
        if ((logic_transition[state][j]) || (logic_transition[j][state])) {
          i++;
          if (i == used_transition[state]) {
            break;
          }
        }
      }
      test_state = j;

      if (predecessor[test_state] == I_DEFAULT) {
        predecessor[test_state] = state;
        state = test_state;
        state_order[state] = nb_used_state++;
      }
    }
  }
  while ((state != 0) || (used_transition[0] < state_transition[0]));

  if (nb_used_state < nb_state) {
    status = false;
    error.update(STAT_parsing[STATP_CHAIN_STRUCTURE]);
  }

  if (!ilogic_transition) {
    for (i = 0;i < nb_state;i++) {
      delete [] logic_transition[i];
    }
    delete [] logic_transition;
  }

  delete [] state_transition;
  delete [] used_transition;
  delete [] predecessor;
  delete [] state_order;

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de l'accessibilite des etats d'une chaine de Markov.
 *
 *  argument : matrice des transitions possibles entre etats.
 *
 *--------------------------------------------------------------*/

void Chain::graph_accessibility_computation(bool **ilogic_transition)

{
  if (!accessibility) {
    bool **logic_transition;
    register int i , j , k;
    int state , test_state , nb_used_state , *state_transition , *used_transition ,
        *predecessor , *state_order;


    // creation de la matrice d'accessibilite

    accessibility = new bool*[nb_state];
    for (i = 0;i < nb_state;i++) {
      accessibility[i] = new bool[nb_state];
      for (j = 0;j < i;j++) {
        accessibility[i][j] = false;
      }
      accessibility[i][i] = true;
      for (j = i + 1;j < nb_state;j++) {
        accessibility[i][j] = false;
      }
    }

    if (ilogic_transition) {
      logic_transition = ilogic_transition;
    }

    else {

      // construction de la matrice des transitions possibles entre etats

      logic_transition = logic_transition_computation();
    }

    // calcul de connexite (algorithme de Tarjan): initialisations

    state_transition = new int[nb_state];

    for (i = 0;i < nb_state;i++) {
      state_transition[i] = 0;
      for (j = 0;j < nb_state;j++) {
        if (logic_transition[i][j]) {
          state_transition[i]++;
        }
      }
    }

    used_transition = new int[nb_state];
    predecessor = new int[nb_state];
    state_order = new int[nb_state];

    for (i = 0;i < nb_state;i++) {
      for (j = 0;j < nb_state;j++) {
        used_transition[j] = 0;
        predecessor[j] = I_DEFAULT;
        state_order[j] = I_DEFAULT;
      }

      state = i;
      predecessor[state] = state;
      state_order[state] = 0;
      nb_used_state = 1;

      do {
        if (used_transition[state] == state_transition[state]) {
          state = predecessor[state];
        }

        else {
          used_transition[state]++;
          j = 0;

          // recherche du prochain etat

          for (k = 0;k < nb_state;k++) {
            if (logic_transition[state][k]) {
              j++;
              if (j == used_transition[state]) {
                break;
              }
            }
          }
          test_state = k;

          if (predecessor[test_state] == I_DEFAULT) {
            predecessor[test_state] = state;
            state = test_state;
            state_order[state] = nb_used_state++;
            accessibility[i][state] = true;
          }
        }
      }
      while ((state != i) || (used_transition[i] < state_transition[i]));
    }

#   ifdef DEBUG
/*    cout << "\n\n";
    for (i = 0;i < nb_state;i++) {
      for (j = 0;j < nb_state;j++) {
        cout << accessibility[i][j] << " ";
      }
      cout << endl;
    }
    cout << endl; */
#   endif

    if (!ilogic_transition) {
      for (i = 0;i < nb_state;i++) {
        delete [] logic_transition[i];
      }
      delete [] logic_transition;
    }

    delete [] state_transition;
    delete [] used_transition;
    delete [] predecessor;
    delete [] state_order;
  }
}


/*--------------------------------------------------------------*
 *
 *  Calcul de l'accessibilite des etats d'une chaine de Markov d'ordre 1.
 *
 *--------------------------------------------------------------*/

void Chain::probability_accessibility_computation()

{
  if (!accessibility) {
    bool stop;
    register int i , j , k;
    int length;
    double sum , *state_seq , *pstate_seq , *states , *pstates ,
           **ptransition , **uniform_transition;


    // creation de la matrice d'accessibilite

    accessibility = new bool*[nb_state];
    for (i = 0;i < nb_state;i++) {
      accessibility[i] = new bool[nb_state];
      for (j = 0;j < i;j++) {
        accessibility[i][j] = false;
      }
      accessibility[i][i] = true;
      for (j = i + 1;j < nb_state;j++) {
        accessibility[i][j] = false;
      }
    }

    // calcul de la matrice des transitions possibles entre etats

    uniform_transition = new double*[nb_state];
    for (i = 0;i < nb_state;i++) {
      uniform_transition[i] = new double[nb_state];

      // cas etat non-absorbant

      if (transition[i][i] < 1.) {
        for (j = 0;j < i;j++) {
          uniform_transition[i][j] = 1. / (double)(nb_state - 1);
        }
        uniform_transition[i][i] = 0.;
        for (j = i + 1;j < nb_state;j++) {
          uniform_transition[i][j] = 1. / (double)(nb_state - 1);
        }
      }

      // cas etat absorbant

      else {
        for (j = 0;j < nb_state;j++) {
          uniform_transition[i][j] = transition[i][j];
        }
      }
    }

    // calcul accessibilite

    state_seq = new double[nb_state];
    pstate_seq = new double[nb_state];

    for (i = 0;i < nb_state;i++) {
      if (uniform_transition[i][i] == 0.) {
        pstates = pstate_seq;
        for (j = 0;j < nb_state;j++) {
          *pstates++ = 0.;
        }
        pstate_seq[i] = 1.;

        stop = false;
        length = 0;

        do {

          // calcul des probabilites des sequences d'etats

          states = state_seq;
          sum = 0.;

          for (j = 0;j < nb_state;j++) {
            pstates = pstate_seq;
            *states = 0.;
            for (k = 0;k < nb_state;k++) {
              *states += uniform_transition[k][j] * *pstates++;
            }

            if (*states > 0.) {
              accessibility[i][j] = true;
            }
            states++;
          }

          for (j = 0;j < nb_state;j++) {
            if (!accessibility[i][j]) {
              break;
            }
          }

          if (j == nb_state) {
            stop = true;
          }

          else {
            pstates = pstate_seq;
            states = state_seq;
            sum = 0.;

            for (j = 0;j < nb_state;j++) {
              sum += fabs(*states - *pstates);
              *pstates++ = *states++;
            }

            if (sum / nb_state < ACCESSIBILITY_THRESHOLD) {
              stop = true;
            }

            length++;
          }
        }
        while ((!stop) && (length < ACCESSIBILITY_LENGTH));
      }
    }

#   ifdef DEBUG
/*    cout << "\n\n";
    for (i = 0;i < nb_state;i++) {
      for (j = 0;j < nb_state;j++) {
        cout << accessibility[i][j] << " ";
      }
      cout << endl;
    }
    cout << endl; */
#   endif

    for (i = 0;i < nb_state;i++) {
      delete [] uniform_transition[i];
    }
    delete [] uniform_transition;

    delete [] state_seq;
    delete [] pstate_seq;
  }
}


/*--------------------------------------------------------------*
 *
 *  Extraction des classes d'une chaine de Markov a partir
 *  de l'accessibilite des etats.
 *
 *  argument : matrice des transitions possibles entre etats.
 *
 *--------------------------------------------------------------*/

void Chain::component_computation(bool **ilogic_transition)

{
  if (nb_component == 0) {
    bool *used_state;
    register int i , j , k;
    int nb_used_state = 0 , state , *bcomponent_nb_state , **bcomponent;


    // calcul de l'accessibilite des etats

    graph_accessibility_computation(ilogic_transition);

    // extraction des classes

    nb_component = 0;
    bcomponent_nb_state = new int[nb_state];

    bcomponent = new int*[nb_state];
    for (i = 0;i < nb_state;i++) {
      bcomponent[i] = new int[nb_state - i];
    }

    used_state = new bool[nb_state];
    for (i = 0;i < nb_state;i++) {
      used_state[i] = false;
    }

    while (nb_used_state < nb_state) {
      bcomponent_nb_state[nb_component] = 0;

      for (i = 0;i < nb_state;i++) {
        if (!used_state[i]) {
          state = i;
          used_state[i] = true;
          nb_used_state++;
          bcomponent[nb_component][bcomponent_nb_state[nb_component]++] = i;
          break;
        }
      }

      for (i = 0;i < nb_state;i++) {
        if ((!used_state[i]) && (accessibility[state][i]) &&
            (accessibility[i][state])) {
          used_state[i] = true;
          nb_used_state++;
          bcomponent[nb_component][bcomponent_nb_state[nb_component]++] = i;
        }
      }

      nb_component++;
    }

    // copie des classes

    component_nb_state = new int[nb_component];
    component = new int*[nb_component];

    for (i = 0;i < nb_component;i++) {
      component_nb_state[i] = bcomponent_nb_state[i];
      component[i] = new int[component_nb_state[i]];
      for (j = 0;j < component_nb_state[i];j++) {
        component[i][j] = bcomponent[i][j];
      }
    }

    delete [] bcomponent_nb_state;
    for (i = 0;i < nb_state;i++) {
      delete [] bcomponent[i];
    }
    delete [] bcomponent;
    delete [] used_state;

    // extraction des types des etats

    state_type = new char[nb_state];

    for (i = 0;i < nb_component;i++) {
      state = component[i][0];
      state_type[state] = 'r';

      for (j = 0;j < nb_component;j++) {
        if (j != i) {
          for (k = 0;k < component_nb_state[j];k++) {
            if (accessibility[state][component[j][k]]) {
              state_type[state] = 't';
              break;
            }
          }
          if (k < component_nb_state[j]) {
            break;
          }
        }
      }

      if ((state_type[state] == 'r') && (component_nb_state[i] == 1)) {
        state_type[state] = 'a';
      }

      for (j = 1;j < component_nb_state[i];j++) {
        state_type[component[i][j]] = state_type[state];
      }
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Application d'un seuil sur les parametres d'une chaine de Markov.
 *
 *  arguments : probabilite minimum, flag semi-Markov.
 *
 *--------------------------------------------------------------*/

void Chain::thresholding(double min_probability , bool semi_markov)

{
  bool stop;
  register int i , j;
  int nb_correction;
  double norm;


  if (min_probability > THRESHOLDING_FACTOR / (double)nb_state) {
    min_probability = THRESHOLDING_FACTOR / (double)nb_state;
  }

  if (type == 'o') {
    do {
      stop = true;
      nb_correction = 0;
      norm = 0.;

      for (i = 0;i < nb_state;i++) {
        if (initial[i] <= min_probability) {
          nb_correction++;
          initial[i] = min_probability;
        }
        else {
          norm += initial[i];
        }
      }

      if (nb_correction > 0) {
        for (i = 0;i < nb_state;i++) {
          if (initial[i] > min_probability) {
            initial[i] *= (1. - nb_correction * min_probability) / norm;
            if (initial[i] < min_probability) {
              stop = false;
            }
          }
        }
      }
    }
    while (!stop);
  }

  for (i = 0;i < nb_row;i++) {
    do {
      stop = true;
      nb_correction = 0;
      norm = 0.;

      for (j = 0;j < nb_state;j++) {
        if (((!semi_markov) || (j != i)) && (accessibility[i][j]) &&
            (transition[i][j] <= min_probability)) {
          nb_correction++;
          transition[i][j] = min_probability;
        }
        else {
          norm += transition[i][j];
        }
      }

      if (nb_correction > 0) {
        for (j = 0;j < nb_state;j++) {
          if (transition[i][j] > min_probability) {
            transition[i][j] *= (1. - nb_correction * min_probability) / norm;
            if (transition[i][j] < min_probability) {
              stop = false;
            }
          }
        }
      }
    }
    while (!stop);
  }

/*  if (accessibility) {
    for (i = 0;i < nb_state;i++) {
      delete [] accessibility[i];
    }
    delete [] accessibility;
  }
  accessibility = NULL;

  delete [] component_nb_state;

  if (component) {
    for (i = 0;i < nb_component;i++) {
      delete [] component[i];
    }
    delete [] component;
  }
  nb_component = 0;

  component_computation(); */
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la vraisemblance des etats initiaux et des transitions
 *  pour une chaine de Markov.
 *
 *  arguments : reference sur un objet ChainData, prise en compte ou non loi initiale.
 *
 *--------------------------------------------------------------*/

double Chain::likelihood_computation(const ChainData &chain_data , bool initial_flag) const

{
  register int i , j;
  double likelihood;


  if ((chain_data.nb_state != nb_state) || (chain_data.nb_row != nb_row)) {
    likelihood = D_INF;
  }

  else {
    likelihood = 0.;

    if (initial_flag) {
      for (i = 0;i < nb_state;i++) {
        if (chain_data.initial[i] > 0) {
          if (initial[i] > 0.) {
            likelihood += chain_data.initial[i] * log(initial[i]);
          }
          else {
            likelihood = D_INF;
            break;
          }
        }
      }
    }

    if (likelihood != D_INF) {
      for (i = 0;i < nb_row;i++) {
        for (j = 0;j < nb_state;j++) {
          if (chain_data.transition[i][j] > 0) {
            if (transition[i][j] > 0.) {
              likelihood += chain_data.transition[i][j] * log(transition[i][j]);
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
 *  Calcul du nombre de probabilites de transitions independantes.
 *
 *  argument : probabilite minimum.
 *
 *--------------------------------------------------------------*/

int Chain::nb_parameter_computation(double min_probability) const

{
  register int i , j;
  int nb_parameter = 0;


  for (i = 0;i < nb_row;i++) {
    for (j = 0;j < nb_state;j++) {
      if (transition[i][j] > min_probability) {
        nb_parameter++;
      }
    }
  }
  nb_parameter -= nb_row;

  return nb_parameter;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la valeur du chi2 pour une chaine de Markov.
 *
 *  argument : reference sur un objet ChainData.
 *
 *--------------------------------------------------------------*/

double Chain::chi2_value_computation(const ChainData &chain_data) const

{
  register int i , j;
  int sum;
  double value , var1 , var2;


  if ((chain_data.nb_state != nb_state) || (chain_data.nb_row != nb_row)) {
    value = -D_INF;
  }

  else {
    value = 0.;

    for (i = 0;i < nb_row;i++) {
      sum = 0;
      for (j = 0;j < nb_state;j++) {
        sum += chain_data.transition[i][j];
      }

      if (sum > 0) {
        for (j = 0;j < nb_state;j++) {
          if (chain_data.transition[i][j] > 0) {
            if (transition[i][j] > 0.) {
              var1 = sum * transition[i][j];
              var2 = chain_data.transition[i][j] - var1;
              value += var2 * var2 / var1;
            }
            else {
              value = -D_INF;
              break;
            }
          }
        }
      }
    }
  }

  return value;
}


/*--------------------------------------------------------------*
 *
 *  Mesure de l'ajustement d'une chaine de Markov a des sequences.
 *
 *  arguments : reference sur un objet ChainData et sur un objet Test.
 *
 *--------------------------------------------------------------*/

void Chain::chi2_fit(const ChainData &chain_data , Test &test) const

{
  if ((chain_data.nb_state == nb_state) || (chain_data.nb_row == nb_row)) {
    test.df1 = nb_parameter_computation();
    test.value = chi2_value_computation(chain_data);
    if ((test.df1 > 0) && (test.value > 0.)) {
      test.chi2_critical_probability_computation();
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Estimation des parametres d'une chaine de Markov d'ordre 1 a partir
 *  des etats initiaux et des transitions.
 *
 *  argument : reference sur un objet Chain.
 *
 *--------------------------------------------------------------*/

void ChainData::estimation(Chain &chain) const

{
  register int i , j;
  int sum;


  // estimation des probabilites initiales

  if (chain.type == 'o') {
    sum = 0;
    for (i = 0;i < nb_state;i++) {
      sum += initial[i];
    }

    for (i = 0;i < nb_state;i++) {
      chain.initial[i] = (double)initial[i] / (double)sum;
    }
  }

  // estimation des probabilites de transition

  for (i = 0;i < nb_state;i++) {
    sum = 0;
    for (j = 0;j < nb_state;j++) {
      sum += transition[i][j];
    }

    if (sum > 0) {
      for (j = 0;j < nb_state;j++) {
        chain.transition[i][j] = (double)transition[i][j] / (double)sum;
      }
    }

    else {
      for (j = 0;j < i;j++) {
        chain.transition[i][j] = 0.;
      }
      chain.transition[i][i] = 1.;
      for (j = i + 1;j < nb_state;j++) {
        chain.transition[i][j] = 0.;
      }
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Calcul du nombre de probabilites de transitions independantes.
 *
 *--------------------------------------------------------------*/

int ChainData::nb_parameter_computation() const

{
  register int i , j;
  int sum , nb_parameter = 0;


  for (i = 0;i < nb_row;i++) {
    sum = 0;
    for (j = 0;j < nb_state;j++) {
      if (transition[i][j] > 0) {
        sum++;
      }
    }

    if (sum > 1) {
      nb_parameter += sum - 1;
    }
  }

  return nb_parameter;
}


};  // namespace stat_tool
