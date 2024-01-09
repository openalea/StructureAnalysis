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
 *       $Id: smc_distributions2.cpp 18078 2015-04-23 10:57:34Z guedon $
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



#include "stat_tool/stat_tools.h"
#include "stat_tool/curves.h"
#include "stat_tool/distribution.h"
#include "stat_tool/markovian.h"
#include "stat_tool/vectors.h"
#include "stat_tool/distance_matrix.h"

#include "sequences.h"
#include "semi_markov.h"
#include "sequence_label.h"

using namespace std;
using namespace stat_tool;


namespace sequence_analysis {



/*--------------------------------------------------------------*
 *
 *  Calcul d'un melange de lois du nombre de series d'un etat ('r') ou
 *  du nombre d'occurrences d'un etat ('o') d'une semi-chaine de Markov
 *  pour une distribution des longueurs de sequences donnee.
 *
 *  arguments : etat, type de forme.
 *
 *--------------------------------------------------------------*/

void SemiMarkovChain::state_nb_pattern_mixture(int state , char pattern)

{
  register int i , j , k , m;
  int max_length , index_nb_pattern , previous_nb_pattern , increment;
  double sum , *pmass , *lmass , **state_out , *pstate_out , ***state_in;
  Distribution *pdist;
  DiscreteParametric *occupancy;


  switch (pattern) {
  case 'r' :
    pdist = state_process->nb_run[state];
    break;
  case 'o' :
    pdist = state_process->nb_occurrence[state];
    break;
  }

  pmass = pdist->mass;
  for (i = 0;i < pdist->nb_value;i++) {
    *pmass++ = 0.;
  }

  max_length = state_process->length->nb_value - 1;

  state_out = new double*[nb_state];
  for (i = 0;i < nb_state;i++) {
    state_out[i] = new double[pattern == 'o' ? max_length : (max_length + 1) / 2 + 1];
  }

  state_in = new double**[max_length - 1];
  index_nb_pattern = 1;

  for (i = 0;i < max_length - 1;i++) {
    state_in[i] = new double*[nb_state];
    for (j = 0;j < nb_state;j++) {
      state_in[i][j] = new double[index_nb_pattern + 1];
    }
    if ((pattern == 'o') || (i % 2 == 1)) {
      index_nb_pattern++;
    }
  }

  // calcul des probabilites de quitter (semi-Markov) / d'etre dans (Markov) un etat
  // fonction du nombre de formes de l'etat selectionne

  lmass = state_process->length->mass;
  index_nb_pattern = 1;

  for (i = 0;i < max_length;i++) {
    lmass++;

    for (j = 0;j < nb_state;j++) {

      // initialisation des probabilites de quitter un etat a l'instant i

      if (i < max_length - 1) {
        pstate_out = state_out[j];
        for (k = 0;k <= index_nb_pattern;k++) {
          *pstate_out++ = 0.;
        }
      }

      switch (state_subtype[j]) {

      // cas etat semi-markovien

      case SEMI_MARKOVIAN : {
        occupancy = state_process->sojourn_time[j];

        for (k = (*lmass > 0. ? 1 : occupancy->offset);k <= MIN(i + 1 , occupancy->nb_value - 1);k++) {
          switch (pattern) {
          case 'r' :
            increment = 1;
            break;
          case 'o' :
            increment = k;
            break;
          }

          if (i < max_length - 1) {
            pstate_out = state_out[j];
            if (j == state) {
              pstate_out += increment;
            }
          }
          if (*lmass > 0.) {
            pmass = pdist->mass;
            if (j == state) {
              pmass += increment;
            }
          }

          if (k < i + 1) {
            switch (pattern) {

            case 'r' : {
              if ((j == state) && (k == 1) && (i % 2 == 1)) {
                previous_nb_pattern = index_nb_pattern - 1;
              }
              else {
                previous_nb_pattern = (i - k) / 2 + 1;
              }
              break;
            }

            case 'o' : {
              previous_nb_pattern = i - k + 1;
              break;
            }
            }

            if (i < max_length - 1) {
              for (m = 0;m <= previous_nb_pattern;m++) {
                *pstate_out++ += occupancy->mass[k] * state_in[i - k][j][m];
              }
            }
            if (*lmass > 0.) {
              for (m = 0;m <= previous_nb_pattern;m++) {
                *pmass++ += *lmass * (1. - occupancy->cumul[k - 1]) * state_in[i - k][j][m];
              }
            }
          }

          else {
            if (i < max_length - 1) {
              switch (type) {
              case 'o' :
                *pstate_out += occupancy->mass[k] * initial[j];
                break;
              case 'e' :
                *pstate_out += forward[j]->mass[k] * initial[j];
                break;
              }
            }

            if (*lmass > 0.) {
              switch (type) {
              case 'o' :
                *pmass += *lmass * (1. - occupancy->cumul[k - 1]) * initial[j];
                break;
              case 'e' :
                *pmass += *lmass * (1. - forward[j]->cumul[k - 1]) * initial[j];
                break;
              }
            }
          }
        }
        break;
      }

      // cas etat markovien

      case MARKOVIAN : {
        if (i < max_length - 1) {
          pstate_out = state_out[j];
          if (j == state) {
            pstate_out++;
          }
        }

        if (*lmass > 0.) {
          pmass = pdist->mass;
          if (j == state) {
            pmass++;
          }
        }

        if (i == 0) {
          *pstate_out = initial[j];
          if (*lmass > 0.) {
            *pmass += *lmass * initial[j];
          }
        }

        else {
          switch (pattern) {

          case 'r' : {
            if ((j == state) && (i % 2 == 1)) {
              previous_nb_pattern = index_nb_pattern - 1;
            }
            else {
              previous_nb_pattern = (i - 1) / 2 + 1;
            }
            break;
          }
  
          case 'o' : {
            previous_nb_pattern = i;
            break;
          }
          }

          if (i < max_length - 1) {
            for (k = 0;k <= previous_nb_pattern;k++) {
              *pstate_out++ = state_in[i - 1][j][k];
            }
          }
          if (*lmass > 0.) {
            for (k = 0;k <= previous_nb_pattern;k++) {
              *pmass++ += *lmass * state_in[i - 1][j][k];
            }
          }
        }
        break;
      }
      }
    }

    if (i < max_length - 1) {
      for (j = 0;j < nb_state;j++) {
        for (k = 0;k <= index_nb_pattern;k++) {
          state_in[i][j][k] = 0.;
          for (m = 0;m < nb_state;m++) {
            if ((pattern == 'o') || (j != state) || (j != m)) {
              state_in[i][j][k] += transition[m][j] * state_out[m][k];
            }
            else if (k < index_nb_pattern) {
              state_in[i][j][k] += transition[m][j] * state_out[m][k + 1];
            }
          }
        }
      }
    }

    if ((pattern == 'o') || (i % 2 == 1)) {
      index_nb_pattern++;
    }
  }

  // renormalisation du melange de lois du nombre de formes de l'etat
  // selectionne pour tenir compte des seuils appliques sur les fonctions
  // de repartition des lois d'occupation des etats

  pmass = pdist->mass;
  sum = 0.;
  for (i = 0;i < pdist->nb_value;i++) {
    sum += *pmass++;
  }

  if (sum < 1.) {
    pmass = pdist->mass;
    for (i = 0;i < pdist->nb_value;i++) {
      *pmass++ /= sum;
    }
  }

  pdist->nb_value_computation();
  pdist->offset_computation();
  pdist->cumul_computation();

  pdist->max_computation();
  pdist->mean_computation();
  pdist->variance_computation();

  for (i = 0;i < nb_state;i++) {
    delete [] state_out[i];
  }
  delete [] state_out;

  for (i = 0;i < max_length - 1;i++) {
    for (j = 0;j < nb_state;j++) {
      delete [] state_in[i][j];
    }
    delete [] state_in[i];
  }
  delete [] state_in;
}


/*--------------------------------------------------------------*
 *
 *  Calcul d'un melange de lois du nombre de series d'une observation
 *  emise par une semi-chaine de Markov cachee pour une distribution
 *  des longueurs de sequences donnee.
 *
 *  arguments : indice du processus d'observation, observation.
 *
 *--------------------------------------------------------------*/

void SemiMarkov::output_nb_run_mixture(int variable , int output)

{
  register int i , j , k , m , n;
  int max_length , index_nb_pattern , min , max;
  double sum0 , sum1 , *pmass , *lmass , **state_out , ***state_in , ***state_nb_run;
  Distribution *nb_run;
  DiscreteParametric *occupancy;


  nb_run = categorical_process[variable]->nb_run[output];

  pmass = nb_run->mass;
  for (i = 0;i < nb_run->nb_value;i++) {
    *pmass++ = 0.;
  }

  max_length = categorical_process[variable]->length->nb_value - 1;

  state_out = new double*[nb_state * 2];
  for (i = 0;i < nb_state * 2;i++) {
    state_out[i] = new double[(max_length + 1) / 2 + 1];
  }

  state_in = new double**[max_length - 1];
  index_nb_pattern = 1;

  for (i = 0;i < max_length - 1;i++) {
    state_in[i] = new double*[nb_state * 2];
    for (j = 0;j < nb_state * 2;j++) {
      state_in[i][j] = new double[index_nb_pattern + 1];
    }
    if (i % 2 == 1) {
      index_nb_pattern++;
    }
  }

  // calcul des lois du nombre de series de l'observation selectionnee
  // pour les differents temps passes dans un etat en tenant compte
  // de l'observation emise avant d'entrer dans l'etat

  state_nb_run = new double**[nb_state * 2];

  for (i = 0;i < nb_state * 2;i++) {
    if (state_subtype[i / 2] == SEMI_MARKOVIAN) {
      occupancy = state_process->sojourn_time[i / 2];
      state_nb_run[i] = new double*[MIN(max_length + 1 , occupancy->nb_value) * 2];

      state_nb_run[i][0] = NULL;
      state_nb_run[i][1] = NULL;
      index_nb_pattern = 1;
      for (j = 1;j < MIN(max_length + 1 , occupancy->nb_value);j++) {
        state_nb_run[i][j * 2] = new double[index_nb_pattern + 1];
        state_nb_run[i][j * 2 + 1] = new double[index_nb_pattern + 1];
        if (j % 2 == 1) {
          index_nb_pattern++;
        }
      }
    }
  }

  for (i = 0;i < nb_state * 2;i++) {
    if (state_subtype[i / 2] == SEMI_MARKOVIAN) {
      switch (i % 2) {

      case 0 : {
        state_nb_run[i][2][0] = 1. - categorical_process[variable]->observation[i / 2]->mass[output];
        state_nb_run[i][2][1] = 0.;
        state_nb_run[i][3][0] = 0.;
        state_nb_run[i][3][1] = categorical_process[variable]->observation[i / 2]->mass[output];
        break;
      }

      case 1 : {
        state_nb_run[i][2][0] = 1. - categorical_process[variable]->observation[i / 2]->mass[output];
        state_nb_run[i][2][1] = 0.;
        state_nb_run[i][3][0] = categorical_process[variable]->observation[i / 2]->mass[output];
        state_nb_run[i][3][1] = 0.;
        break;
      }
      }

      occupancy = state_process->sojourn_time[i / 2];
      index_nb_pattern = 1;

      for (j = 2;j < MIN(max_length + 1 , occupancy->nb_value);j++) {
        for (k = 0;k <= index_nb_pattern;k++) {
          state_nb_run[i][j * 2][k] = (1. - categorical_process[variable]->observation[i / 2]->mass[output]) *
                                      (state_nb_run[i][j * 2 - 2][k] + state_nb_run[i][j * 2 - 1][k]);

          sum0 = state_nb_run[i][j * 2 - 1][k];
          if (k > 0) {
            sum0 += state_nb_run[i][j * 2 - 2][k - 1];
          }
          state_nb_run[i][j * 2 + 1][k] = categorical_process[variable]->observation[i / 2]->mass[output] * sum0;
        }

        if (j % 2 == 0) {
          index_nb_pattern++;
          state_nb_run[i][j * 2][index_nb_pattern] = 0.;
          state_nb_run[i][j * 2 + 1][index_nb_pattern] = 0.;
        }
      }
    }
  }

  // calcul des probabilites de quitter (semi-Markov) / d'etre dans (Markov) un etat
  // fonction du nombre de series de l'observation selectionnee

  lmass = categorical_process[variable]->length->mass;
  index_nb_pattern = 1;

  for (i = 0;i < max_length;i++) {
    lmass++;

    // initialisation des probabilites de quitter un etat a l'instant i

    for (j = 0;j < nb_state * 2;j++) {
      for (k = 0;k <= index_nb_pattern;k++) {
        state_out[j][k] = 0.;
      }
    }

    for (j = 0;j < nb_state;j++) {
      switch (state_subtype[j]) {

      // cas etat semi-markovien

      case SEMI_MARKOVIAN : {
        occupancy = state_process->sojourn_time[j];

        for (k = (*lmass > 0. ? 1 : occupancy->offset);k <= MIN(i + 1 , occupancy->nb_value - 1);k++) {
          if (*lmass > 0.) {
            pmass = nb_run->mass;
          }

          for (m = 0;m <= index_nb_pattern;m++) {
            if (k < i + 1) {
              min = MAX(m - ((i - k) / 2 + 1) , 0);
              max = MIN((k % 2 == 0 ? k / 2 : k / 2 + 1) , m);

              if (max >= min) {
                sum0 = 0.;
                sum1 = 0.;
                for (n = min;n <= max;n++) {
                  sum0 += state_nb_run[j * 2][k * 2][n] * state_in[i - k][j * 2][m - n] +
                          state_nb_run[j * 2 + 1][k * 2][n] * state_in[i - k][j * 2 + 1][m - n];
                  sum1 += state_nb_run[j * 2][k * 2 + 1][n] * state_in[i - k][j * 2][m - n] +
                          state_nb_run[j * 2 + 1][k * 2 + 1][n] * state_in[i - k][j * 2 + 1][m - n];
                }

                if (i < max_length - 1) {
                  state_out[j * 2][m] += occupancy->mass[k] * sum0;
                  state_out[j * 2 + 1][m] += occupancy->mass[k] * sum1;
                }
                if (*lmass > 0.) {
                  *pmass += *lmass * (1. - occupancy->cumul[k - 1]) * (sum0 + sum1);
                }
              }
            }

            else {
              sum0 = state_nb_run[j * 2][k * 2][m] * initial[j];
              sum1 = state_nb_run[j * 2][k * 2 + 1][m] * initial[j];

              if (i < max_length - 1) {
                switch (type) {
                case 'o' :
                  state_out[j * 2][m] += occupancy->mass[k] * sum0;
                  state_out[j * 2 + 1][m] += occupancy->mass[k] * sum1;
                  break;
                case 'e' :
                  state_out[j * 2][m] += forward[j]->mass[k] * sum0;
                  state_out[j * 2 + 1][m] += forward[j]->mass[k] * sum1;
                  break;
                }
              }

              if (*lmass > 0.) {
                switch (type) {
                case 'o' :
                  *pmass += *lmass * (1. - occupancy->cumul[k - 1]) * (sum0 + sum1);
                  break;
                case 'e' :
                  *pmass += *lmass * (1. - forward[j]->cumul[k - 1]) * (sum0 + sum1);
                  break;
                }
              }
            }

            if (*lmass > 0.) {
              pmass++;
            }
          }
        }
        break;
      }

      // cas etat markovien

      case MARKOVIAN : {
        if (*lmass > 0.) {
          pmass = nb_run->mass;
        }

        if (i == 0) {
          state_out[j * 2][0] = (1. - categorical_process[variable]->observation[j]->mass[output]) * initial[j];
          state_out[j * 2 + 1][1] = categorical_process[variable]->observation[j]->mass[output] * initial[j];

           if (*lmass > 0.) {
            *pmass++ += *lmass * state_out[j * 2][0];
            *pmass += *lmass * state_out[j * 2 + 1][1];
          }
        }

        else {
          for (k = 0;k <= index_nb_pattern;k++) {
            sum0 = 0.;
            if ((k < index_nb_pattern) || (i % 2 == 1)) {
              state_out[j * 2][k] = (1. - categorical_process[variable]->observation[j]->mass[output]) *
                                    (state_in[i - 1][j * 2][k] + state_in[i - 1][j * 2 + 1][k]);
              sum0 += state_in[i - 1][j * 2 + 1][k];
            }
            if (k > 0) {
              sum0 += state_in[i - 1][j * 2][k - 1];
            }
            state_out[j * 2 + 1][k] = categorical_process[variable]->observation[j]->mass[output] * sum0;

            if (*lmass > 0.) {
              *pmass++ += *lmass * (state_out[j * 2][k] + state_out[j * 2 + 1][k]);
            }
          }
        }
        break;
      }
      }
    }

    if (i < max_length - 1) {
      for (j = 0;j < nb_state;j++) {
        for (k = 0;k <= index_nb_pattern;k++) {
          state_in[i][j * 2][k] = 0.;
          state_in[i][j * 2 + 1][k] = 0.;
          for (m = 0;m < nb_state;m++) {
            state_in[i][j * 2][k] += transition[m][j] * state_out[m * 2][k];
            state_in[i][j * 2 + 1][k] += transition[m][j] * state_out[m * 2 + 1][k];
          }
        }
      }
    }

    if (i % 2 == 1) {
      index_nb_pattern++;
    }
  }

  // renormalisation du melange de lois du nombre de series de l'etat
  // selectionne pour tenir compte des seuils appliques sur les fonctions
  // de repartition des lois d'occupation des etats

  pmass = nb_run->mass;
  sum0 = 0.;
  for (i = 0;i < nb_run->nb_value;i++) {
    sum0 += *pmass++;
  }

  if (sum0 < 1.) {
    pmass = nb_run->mass;
    for (i = 0;i < nb_run->nb_value;i++) {
      *pmass++ /= sum0;
    }
  }

  nb_run->nb_value_computation();
  nb_run->offset_computation();
  nb_run->cumul_computation();

  nb_run->max_computation();
  nb_run->mean_computation();
  nb_run->variance_computation();

  for (i = 0;i < nb_state * 2;i++) {
    if (state_subtype[i / 2] == SEMI_MARKOVIAN) {
      occupancy = state_process->sojourn_time[i / 2];
      for (j = 1;j < MIN(max_length + 1 , occupancy->nb_value);j++) {
        delete [] state_nb_run[i][j * 2];
        delete [] state_nb_run[i][j * 2 + 1];
      }
      delete [] state_nb_run[i];
    }
  }
  delete [] state_nb_run;

  for (i = 0;i < nb_state * 2;i++) {
    delete [] state_out[i];
  }
  delete [] state_out;

  for (i = 0;i < max_length - 1;i++) {
    for (j = 0;j < nb_state * 2;j++) {
      delete [] state_in[i][j];
    }
    delete [] state_in[i];
  }
  delete [] state_in;
}


/*--------------------------------------------------------------*
 *
 *  Calcul d'un melange de lois du nombre d'occurrences d'une observation
 *  emise par une semi-chaine de Markov cachee pour une distribution
 *  des longueurs de sequences donnee.
 *
 *  arguments : indice du processus d'observation, observation.
 *
 *--------------------------------------------------------------*/

void SemiMarkov::output_nb_occurrence_mixture(int variable , int output)

{
  register int i , j , k , m , n;
  int max_length , min , max;
  double sum , *pmass , *omass , *lmass , **state_out , ***state_in;
  Distribution *nb_occurrence;
  DiscreteParametric *occupancy , ***observation;


  nb_occurrence = categorical_process[variable]->nb_occurrence[output];

  pmass = nb_occurrence->mass;
  for (i = 0;i < nb_occurrence->nb_value;i++) {
    *pmass++ = 0.;
  }

  max_length = categorical_process[variable]->length->nb_value - 1;

  state_out = new double*[nb_state];
  for (i = 0;i < nb_state;i++) {
    state_out[i] = new double[max_length + 1];
  }

  state_in = new double**[max_length - 1];
  for (i = 0;i < max_length - 1;i++) {
    state_in[i] = new double*[nb_state];
    for (j = 0;j < nb_state;j++) {
      state_in[i][j] = new double[i + 2];
    }
  }

  // calcul des lois du nombre d'occurrences de l'observation selectionnee
  // pour les differents temps passes dans un etat donne

  observation = new DiscreteParametric**[nb_state];
  for (i = 0;i < nb_state;i++) {
    if (state_subtype[i] == SEMI_MARKOVIAN) {
      occupancy = state_process->sojourn_time[i];
      observation[i] = new DiscreteParametric*[MIN(max_length + 1 , occupancy->nb_value)];

      observation[i][0] = NULL;
      for (j = 1;j < MIN(max_length + 1 , occupancy->nb_value);j++) {
        observation[i][j] = new DiscreteParametric(BINOMIAL , 0 , j , D_DEFAULT ,
                                                   categorical_process[variable]->observation[i]->mass[output]);
      }
    }
  }

  // calcul des probabilites de quitter (semi-Markov) / d'etre dans (Markov) un etat
  // fonction du nombre d'occurrences de l'observation selectionnee

  lmass = categorical_process[variable]->length->mass;

  for (i = 0;i < max_length;i++) {
    lmass++;

    for (j = 0;j < nb_state;j++) {

      // initialisation des probabilites de quitter un etat a l'instant i

      for (k = 0;k <= i + 1;k++) {
        state_out[j][k] = 0.;
      }

      switch (state_subtype[j]) {

      // cas etat semi-markovien

      case SEMI_MARKOVIAN : {
        occupancy = state_process->sojourn_time[j];

        for (k = (*lmass > 0. ? 1 : occupancy->offset);k <= MIN(i + 1 , occupancy->nb_value - 1);k++) {
          if (*lmass > 0.) {
            pmass = nb_occurrence->mass;
          }

          for (m = 0;m <= i + 1;m++) {
            if (k < i + 1) {
              min = MAX(m - (i - k + 1) , 0);
              max = MIN(k , m);

              if (max >= min) {
                omass = observation[j][k]->mass + min;
                sum = 0.;
                for (n = min;n <= max;n++) {
                  sum += *omass++ * state_in[i - k][j][m - n];
                }

                if (i < max_length - 1) {
                  state_out[j][m] += occupancy->mass[k] * sum;
                }
                if (*lmass > 0.) {
                  *pmass += *lmass * (1. - occupancy->cumul[k - 1]) * sum;
                }
              }
            }

            else {
              sum = observation[j][k]->mass[m] * initial[j];

              if (i < max_length - 1) {
                switch (type) {
                case 'o' :
                  state_out[j][m] += occupancy->mass[k] * sum;
                  break;
                case 'e' :
                  state_out[j][m] += forward[j]->mass[k] * sum;
                  break;
                }
              }

              if (*lmass > 0.) {
                switch (type) {
                case 'o' :
                  *pmass += *lmass * (1. - occupancy->cumul[k - 1]) * sum;
                  break;
                case 'e' :
                  *pmass += *lmass * (1. - forward[j]->cumul[k - 1]) * sum;
                  break;
                }
              }
            }

            if (*lmass > 0.) {
              pmass++;
            }
          }
        }
        break;
      }

      // cas etat markovien

      case MARKOVIAN : {
        if (*lmass > 0.) {
          pmass = nb_occurrence->mass;
        }

        if (i == 0) {
          state_out[j][0] = (1. - categorical_process[variable]->observation[j]->mass[output]) * initial[j];
          state_out[j][1] = categorical_process[variable]->observation[j]->mass[output] * initial[j];

           if (*lmass > 0.) {
            *pmass++ += *lmass * state_out[j][0];
            *pmass += *lmass * state_out[j][1];
          }
        }

        else {
          for (k = 0;k <= i + 1;k++) {
            if (k < i + 1) {
              state_out[j][k] += (1. - categorical_process[variable]->observation[j]->mass[output]) *
                                 state_in[i - 1][j][k];
            }
            if (k > 0) {
              state_out[j][k] += categorical_process[variable]->observation[j]->mass[output] *
                                 state_in[i - 1][j][k - 1];
            }

            if (*lmass > 0.) {
              *pmass++ += *lmass * state_out[j][k];
            }
          }
        }
        break;
      }
      }
    }

    if (i < max_length - 1) {
      for (j = 0;j < nb_state;j++) {
        for (k = 0;k <= i + 1;k++) {
          state_in[i][j][k] = 0.;
          for (m = 0;m < nb_state;m++) {
            state_in[i][j][k] += transition[m][j] * state_out[m][k];
          }
        }
      }
    }
  }

  // renormalisation du melange de lois du nombre d'occurrences de l'etat
  // selectionne pour tenir compte des seuils appliques sur les fonctions
  // de repartition des lois d'occupation des etats

  pmass = nb_occurrence->mass;
  sum = 0.;
  for (i = 0;i < nb_occurrence->nb_value;i++) {
    sum += *pmass++;
  }

  if (sum < 1.) {
    pmass = nb_occurrence->mass;
    for (i = 0;i < nb_occurrence->nb_value;i++) {
      *pmass++ /= sum;
    }
  }

  nb_occurrence->nb_value_computation();
  nb_occurrence->offset_computation();
  nb_occurrence->cumul_computation();

  nb_occurrence->max_computation();
  nb_occurrence->mean_computation();
  nb_occurrence->variance_computation();

  for (i = 0;i < nb_state;i++) {
    if (state_subtype[i] == SEMI_MARKOVIAN) {
      for (j = 1;j < MIN(max_length + 1 , state_process->sojourn_time[i]->nb_value);j++) {
        delete observation[i][j];
      }
      delete [] observation[i];
    }
  }
  delete [] observation;

  for (i = 0;i < nb_state;i++) {
    delete [] state_out[i];
  }
  delete [] state_out;

  for (i = 0;i < max_length - 1;i++) {
    for (j = 0;j < nb_state;j++) {
      delete [] state_in[i][j];
    }
    delete [] state_in[i];
  }
  delete [] state_in;
}


/*--------------------------------------------------------------*
 *
 *  Calcul des lois caracteristiques d'un objet SemiMarkov.
 *
 *  arguments : longueur des sequences, flag sur le calcul des lois de comptage,
 *              indice du processus d'observation.
 *
 *--------------------------------------------------------------*/

void SemiMarkov::characteristic_computation(int length , bool counting_flag , int variable)

{
  if (nb_component > 0) {
    bool computation[NB_OUTPUT_PROCESS + 1];
    register int i , j , k;
    double *memory;
    DiscreteParametric dlength(UNIFORM , length , length , D_DEFAULT , D_DEFAULT);


    memory = NULL;

    // calcul des lois de type intensite et intervalle au niveau etat

    if (((variable == I_DEFAULT) || (variable == 0)) &&
        ((!(state_process->length)) ||
         (dlength != *(state_process->length)))) {
      computation[0] = true;
      state_process->create_characteristic(dlength , false , counting_flag);

      index_state_distribution();

      for (i = 0;i < nb_state;i++) {
        if (type == 'o') {
          state_no_occurrence_probability(i);
        }
        state_first_occurrence_distribution(i);

        if (type == 'o') {
          state_leave_probability(i);
        }
        if (state_process->leave[i] < 1. - DOUBLE_ERROR) {
          state_recurrence_time_distribution(i);
        }
        else {
          delete state_process->recurrence_time[i];
          state_process->recurrence_time[i] = NULL;
        }

        if ((state_subtype[i] == MARKOVIAN) && (transition[i][i] < 1.)) {
          if (transition[i][i] > 0.) {
            state_process->sojourn_time[i] = new DiscreteParametric(NEGATIVE_BINOMIAL , 1 ,
                                                                    I_DEFAULT , 1. , 1. - transition[i][i] ,
                                                                    OCCUPANCY_THRESHOLD);
            state_process->sojourn_time[i]->parameter = D_DEFAULT;
            state_process->sojourn_time[i]->probability = D_DEFAULT;
          }

          else {
            state_process->sojourn_time[i] = new DiscreteParametric(UNIFORM , 1 , 1 ,
                                                                    D_DEFAULT , D_DEFAULT);
            state_process->sojourn_time[i]->sup_bound = I_DEFAULT;
          }

          state_process->sojourn_time[i]->ident = CATEGORICAL;
          state_process->sojourn_time[i]->inf_bound = I_DEFAULT;
        }
      }

#     ifdef MESSAGE
      if (type == 'e') {
        double sum = 0.;

        // calcul de la loi stationnaire dans le cas d'un processus en equilibre
        // avec renormalisation pour tenir compte des seuils appliques sur
        // les fonctions de repartition des lois de temps de retour dans les etats

        for (i = 0;i < nb_state;i++) {
          sum += 1. / state_process->recurrence_time[i]->mean;
        }

        cout << "\n" << STAT_label[STATL_STATIONARY_PROBABILITIES] << endl;
        for (i = 0;i < nb_state;i++) {
          cout << initial[i] << " | "
               << 1. / (state_process->recurrence_time[i]->mean * sum) << endl;
        }
      }
#     endif

    }

    else {
      computation[0] = false;
    }

    // calcul des lois de type intensite et intervalle au niveau observation

    for (i = 0;i < nb_output_process;i++) {
      if ((categorical_process[i]) && ((variable == I_DEFAULT) || (i == variable)) &&
          ((!(categorical_process[i]->length)) ||
           (dlength != *(categorical_process[i]->length)))) {
        computation[i + 1] = true;
        categorical_process[i]->create_characteristic(dlength , true , counting_flag);

        index_output_distribution(i);

        if (!memory) {
          memory = memory_computation();
        }

        for (j = 0;j < categorical_process[i]->nb_value;j++) {
          if (type == 'o') {
            output_no_occurrence_probability(i , j);
          }
          if (categorical_process[i]->no_occurrence[j] < 1. - DOUBLE_ERROR) {
            output_first_occurrence_distribution(i , j);
          }
          else {
            delete categorical_process[i]->first_occurrence[j];
            categorical_process[i]->first_occurrence[j] = NULL;
            categorical_process[i]->leave[j] = 1.;
          }

          if ((type == 'o') && (categorical_process[i]->first_occurrence[j])) {
            output_leave_probability(memory , i , j);
          }
          if (categorical_process[i]->leave[j] < 1. - DOUBLE_ERROR) {
            output_recurrence_time_distribution(memory , i , j);
          }
          else {
            delete categorical_process[i]->recurrence_time[j];
            categorical_process[i]->recurrence_time[j] = NULL;
          }

          for (k = 0;k < nb_state;k++) {
            if ((categorical_process[i]->observation[k]->mass[j] > 0.) &&
                ((state_type[k] != 'a') || (categorical_process[i]->observation[k]->mass[j] < 1.))) {
              break;
            }
          }

          if (k < nb_state) {
            output_sojourn_time_distribution(memory , i , j);
          }
          else {
            categorical_process[i]->absorption[j] = 1.;
            delete categorical_process[i]->sojourn_time[j];
            categorical_process[i]->sojourn_time[j] = NULL;
          }
        }
      }

      else {
        computation[i + 1] = false;
      }
    }

    delete [] memory;

    if (counting_flag) {

      // calcul des lois de comptage au niveau etat

      if (computation[0]) {
        for (i = 0;i < nb_state;i++) {
          state_nb_pattern_mixture(i , 'r');
          state_nb_pattern_mixture(i , 'o');
        }
      }

      // calcul des lois de comptage au niveau observation

      for (i = 0;i < nb_output_process;i++) {
        if (computation[i + 1]) {
          for (j = 0;j < categorical_process[i]->nb_value;j++) {
            output_nb_run_mixture(i , j);
            output_nb_occurrence_mixture(i , j);
          }
        }
      }
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Calcul des lois caracteristiques d'un objet SemiMarkov.
 *
 *  arguments : reference sur un objet SemiMarkovData,
 *              flag sur le calcul des lois de comptage,
 *              indice du processus d'observation,
 *              flag pour tenir compte des longueurs.
 *
 *--------------------------------------------------------------*/

void SemiMarkov::characteristic_computation(const SemiMarkovData &seq , bool counting_flag ,
                                            int variable , bool length_flag)

{
  if (nb_component > 0) {
    bool computation[NB_OUTPUT_PROCESS + 1];
    register int i , j , k;
    int seq_variable;
    double *memory;
    Distribution dlength(*(seq.length_distribution));


    memory = NULL;

    // calcul des lois de type intensite et intervalle au niveau etat

    if (((variable == I_DEFAULT) || (variable == 0)) && ((!length_flag) ||
         ((length_flag) && ((!(state_process->length)) ||
           (dlength != *(state_process->length)))))) {
      computation[0] = true;
      state_process->create_characteristic(dlength , false , counting_flag);

      index_state_distribution();

      for (i = 0;i < nb_state;i++) {
        if (type == 'o') {
          state_no_occurrence_probability(i);
        }
        if (seq.type[0] == STATE) {
          state_first_occurrence_distribution(i , ((seq.characteristics[0]) && (i < seq.marginal_distribution[0]->nb_value) && (seq.characteristics[0]->first_occurrence[i]) && (seq.characteristics[0]->first_occurrence[i]->nb_element > 0) ? seq.characteristics[0]->first_occurrence[i]->nb_value : 1));
        }
        else {
          state_first_occurrence_distribution(i);
        }

        if (type == 'o') {
          state_leave_probability(i);
        }
        if (state_process->leave[i] < 1. - DOUBLE_ERROR) {
          if (seq.type[0] == STATE) {
            state_recurrence_time_distribution(i , ((seq.characteristics[0]) && (i < seq.marginal_distribution[0]->nb_value) && (seq.characteristics[0]->recurrence_time[i]->nb_element > 0) ? seq.characteristics[0]->recurrence_time[i]->nb_value : 1));
          }
          else {
            state_recurrence_time_distribution(i);
          }
        }
        else {
          delete state_process->recurrence_time[i];
          state_process->recurrence_time[i] = NULL;
        }

        if ((state_subtype[i] == MARKOVIAN) && (transition[i][i] < 1.)) {
          if (transition[i][i] > 0.) {
            state_process->sojourn_time[i] = new DiscreteParametric(NEGATIVE_BINOMIAL , 1 ,
                                                                    I_DEFAULT , 1. , 1. - transition[i][i] ,
                                                                    OCCUPANCY_THRESHOLD);

            if ((seq.type[0] == STATE) && (seq.characteristics[0]) && (i < seq.marginal_distribution[0]->nb_value) &&
                (seq.characteristics[0]->sojourn_time[i]->nb_value > state_process->sojourn_time[i]->nb_value)) {
              state_process->sojourn_time[i]->computation(seq.characteristics[0]->sojourn_time[i]->nb_value , OCCUPANCY_THRESHOLD);
            }
            state_process->sojourn_time[i]->parameter = D_DEFAULT;
            state_process->sojourn_time[i]->probability = D_DEFAULT;
          }

          else {
            state_process->sojourn_time[i] = new DiscreteParametric(UNIFORM , 1 , 1 ,
                                                                    D_DEFAULT , D_DEFAULT);
            state_process->sojourn_time[i]->sup_bound = I_DEFAULT;
          }

          state_process->sojourn_time[i]->ident = CATEGORICAL;
          state_process->sojourn_time[i]->inf_bound = I_DEFAULT;
        }
      }

#     ifdef MESSAGE
      if (type == 'e') {
        double sum = 0.;

        // calcul de la loi stationnaire dans le cas d'un processus en equilibre
        // avec renormalisation pour tenir compte des seuils appliques sur
        // les fonctions de repartition des lois de temps de retour dans les etats

        for (i = 0;i < nb_state;i++) {
          sum += 1. / state_process->recurrence_time[i]->mean;
        }

        cout << "\n" << STAT_label[STATL_STATIONARY_PROBABILITIES] << endl;
        for (i = 0;i < nb_state;i++) {
          cout << initial[i] << " | "
               << 1. / (state_process->recurrence_time[i]->mean * sum) << endl;
        }
      }
#     endif

    }

    else {
      computation[0] = false;
    }

    // calcul des lois de type intensite et intervalle au niveau observation

    for (i = 0;i < nb_output_process;i++) {
      if ((categorical_process[i]) && ((variable == I_DEFAULT) || (i == variable)) &&
          ((!length_flag) || ((length_flag) && ((!(categorical_process[i]->length)) ||
             (dlength != *(categorical_process[i]->length)))))) {
        computation[i + 1] = true;
        categorical_process[i]->create_characteristic(dlength , true , counting_flag);

        switch (seq.type[0]) {
        case STATE :
          seq_variable = i + 1;
          break;
        default :
          seq_variable = i;
          break;
        }

        index_output_distribution(i);

        if (!memory) {
          memory = memory_computation();
        }

        for (j = 0;j < categorical_process[i]->nb_value;j++) {
          if (type == 'o') {
            output_no_occurrence_probability(i , j);
          }
          if (categorical_process[i]->no_occurrence[j] < 1. - DOUBLE_ERROR) {
            output_first_occurrence_distribution(i , j , ((seq.characteristics[seq_variable]) && (j < seq.characteristics[seq_variable]->nb_value) && (seq.characteristics[seq_variable]->first_occurrence[j]->nb_element > 0) ? seq.characteristics[seq_variable]->first_occurrence[j]->nb_value : 1));
          }
          else {
            delete categorical_process[i]->first_occurrence[j];
            categorical_process[i]->first_occurrence[j] = NULL;
            categorical_process[i]->leave[j] = 1.;
          }

          if ((type == 'o') && (categorical_process[i]->first_occurrence[j])) {
            output_leave_probability(memory , i , j);
          }
          if (categorical_process[i]->leave[j] < 1. - DOUBLE_ERROR) {
            output_recurrence_time_distribution(memory , i , j , ((seq.characteristics[seq_variable]) && (j < seq.characteristics[seq_variable]->nb_value) && (seq.characteristics[seq_variable]->recurrence_time[j]->nb_element > 0) ? seq.characteristics[seq_variable]->recurrence_time[j]->nb_value : 1));
          }
          else {
            delete categorical_process[i]->recurrence_time[j];
            categorical_process[i]->recurrence_time[j] = NULL;
          }

          for (k = 0;k < nb_state;k++) {
            if ((categorical_process[i]->observation[k]->mass[j] > 0.) &&
                ((state_type[k] != 'a') || (categorical_process[i]->observation[k]->mass[j] < 1.))) {
              break;
            }
          }

          if (k < nb_state) {
            output_sojourn_time_distribution(memory , i , j , ((seq.characteristics[seq_variable]) && (j < seq.characteristics[seq_variable]->nb_value) && (seq.characteristics[seq_variable]->sojourn_time[j]->nb_element > 0) ? seq.characteristics[seq_variable]->sojourn_time[j]->nb_value : 1));
          }
          else {
            categorical_process[i]->absorption[j] = 1.;
            delete categorical_process[i]->sojourn_time[j];
            categorical_process[i]->sojourn_time[j] = NULL;
          }
        }
      }

      else {
        computation[i + 1] = false;
      }
    }

    delete [] memory;

    if (counting_flag) {

      // calcul des lois de comptage au niveau etat

      if (computation[0]) {
        for (i = 0;i < nb_state;i++) {
          state_nb_pattern_mixture(i , 'r');
          state_nb_pattern_mixture(i , 'o');
        }
      }

      // calcul des lois de comptage au niveau observation

      for (i = 0;i < nb_output_process;i++) {
        if (computation[i + 1]) {
          for (j = 0;j < categorical_process[i]->nb_value;j++) {
            output_nb_run_mixture(i , j);
            output_nb_occurrence_mixture(i , j);
          }
        }
      }
    }
  }
}


};  // namespace sequence_analysis
