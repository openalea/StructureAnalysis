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
 *       $Id: hvomc_algorithms2.cpp 18057 2015-04-23 09:47:33Z guedon $
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
#include <iomanip>

#include "tool/config.h"

#include "stat_tool/stat_tools.h"
#include "stat_tool/curves.h"
#include "stat_tool/distribution.h"
#include "stat_tool/markovian.h"
#include "stat_tool/vectors.h"
#include "stat_tool/distance_matrix.h"
#include "stat_tool/stat_label.h"

#include "sequences.h"
#include "variable_order_markov.h"
#include "hidden_variable_order_markov.h"
#include "sequence_label.h"

using namespace std;
using namespace stat_tool;


namespace sequence_analysis {



/*--------------------------------------------------------------*
 *
 *  Calcul des entropies des sequences d'etats par l'algorithme forward-backward.
 *
 *  argument : reference sur un objet VariableOrderMarkovData.
 *
 *--------------------------------------------------------------*/

void HiddenVariableOrderMarkov::forward_backward(VariableOrderMarkovData &seq) const

{
  register int i , j , k , m;
  int **pioutput;
  double seq_likelihood , observation , **forward , norm , **predicted , buff ,
         *transition_predicted , **forward_state_entropy , **proutput;

# ifdef MESSAGE
  double entropy , **backward , *auxiliary , **transition_entropy;
# endif


  // initialisations

  seq.entropy = new double[seq.nb_sequence];
  seq.nb_state_sequence = new double[seq.nb_sequence];

  forward = new double*[seq.max_length];
  for (i = 0;i < seq.max_length;i++) {
    forward[i] = new double[nb_row];
  }

  predicted = new double*[seq.max_length];
  for (i = 0;i < seq.max_length;i++) {
    predicted[i] = new double[nb_row];
  }

  transition_predicted = new double[nb_row];

  forward_state_entropy = new double*[seq.max_length];
  for (i = 0;i < seq.max_length;i++) {
    forward_state_entropy[i] = new double[nb_row];
  }

# ifdef MESSAGE
  backward = new double*[seq.max_length];
  for (i = 0;i < seq.max_length;i++) {
    backward[i] = new double[nb_row];
  }

  auxiliary = new double[nb_row];

  transition_entropy = new double*[nb_row];
  for (i = 1;i < nb_row;i++) {
    transition_entropy[i] = new double[nb_state];
  }
# endif

  pioutput = new int*[nb_output_process];
  proutput = new double*[nb_output_process];

  seq.sample_entropy = 0.;

  for (i = 0;i < seq.nb_sequence;i++) {
    for (j = 0;j < nb_output_process;j++) {
      switch (seq.type[j + 1]) {
      case INT_VALUE :
        pioutput[j] = seq.int_sequence[i][j + 1];
        break;
      case REAL_VALUE :
        proutput[j] = seq.real_sequence[i][j + 1];
        break;
      }
    }

    // recurrence "forward"

    seq_likelihood = 0.;
    norm = 0.;

    switch (type) {

    case 'o' : {
      for (j = 1;j < nb_row;j++) {
        if (order[j] == 1) {
          forward[0][j] = initial[state[j][0]];

          for (k = 0;k < nb_output_process;k++) {
            if (categorical_process[k]) {
              forward[0][j] *= categorical_process[k]->observation[state[j][0]]->mass[*pioutput[k]];
            }

            else if (discrete_parametric_process[k]) {
              forward[0][j] *= discrete_parametric_process[k]->observation[state[j][0]]->mass[*pioutput[k]];
            }

            else {
              if (((continuous_parametric_process[k]->ident == GAMMA) ||
                  (continuous_parametric_process[k]->ident == ZERO_INFLATED_GAMMA)) && (seq.min_value[k + 1] < seq.min_interval[k + 1] / 2)) {
                switch (seq.type[k + 1]) {
                case INT_VALUE :
                  forward[0][j] *= continuous_parametric_process[k]->observation[state[j][0]]->mass_computation(*pioutput[k] , *pioutput[k] + seq.min_interval[k + 1]);
                  break;
                case REAL_VALUE :
                  forward[0][j] *= continuous_parametric_process[k]->observation[state[j][0]]->mass_computation(*proutput[k] , *proutput[k] + seq.min_interval[k + 1]);
                  break;
                }
              }

              else {
                switch (seq.type[k + 1]) {
                case INT_VALUE :
                  forward[0][j] *= continuous_parametric_process[k]->observation[state[j][0]]->mass_computation(*pioutput[k] - seq.min_interval[k + 1] / 2 , *pioutput[k] + seq.min_interval[k + 1] / 2);
                  break;
                case REAL_VALUE :
                  forward[0][j] *= continuous_parametric_process[k]->observation[state[j][0]]->mass_computation(*proutput[k] - seq.min_interval[k + 1] / 2 , *proutput[k] + seq.min_interval[k + 1] / 2);
                  break;
                }
              }
            }
          }

          norm += forward[0][j];
        }

        else {
          forward[0][j] = 0.;
        }
      }
      break;
    }

    case 'e' : {
      for (j = 1;j < nb_row;j++) {
        if (!child[j]) {
          forward[0][j] = initial[j];

          for (k = 0;k < nb_output_process;k++) {
            if (categorical_process[k]) {
              forward[0][j] *= categorical_process[k]->observation[state[j][0]]->mass[*pioutput[k]];
            }

            else if (discrete_parametric_process[k]) {
              forward[0][j] *= discrete_parametric_process[k]->observation[state[j][0]]->mass[*pioutput[k]];
            }

            else {
              if (((continuous_parametric_process[k]->ident == GAMMA) ||
                  (continuous_parametric_process[k]->ident == ZERO_INFLATED_GAMMA)) && (seq.min_value[k + 1] < seq.min_interval[k + 1] / 2)) {
                switch (seq.type[k + 1]) {
                case INT_VALUE :
                  forward[0][j] *= continuous_parametric_process[k]->observation[state[j][0]]->mass_computation(*pioutput[k] , *pioutput[k] + seq.min_interval[k + 1]);
                  break;
                case REAL_VALUE :
                  forward[0][j] *= continuous_parametric_process[k]->observation[state[j][0]]->mass_computation(*proutput[k] , *proutput[k] + seq.min_interval[k + 1]);
                  break;
                }
              }

              else {
                switch (seq.type[k + 1]) {
                case INT_VALUE :
                  forward[0][j] *= continuous_parametric_process[k]->observation[state[j][0]]->mass_computation(*pioutput[k] - seq.min_interval[k + 1] / 2 , *pioutput[k] + seq.min_interval[k + 1] / 2);
                  break;
                case REAL_VALUE :
                  forward[0][j] *= continuous_parametric_process[k]->observation[state[j][0]]->mass_computation(*proutput[k] - seq.min_interval[k + 1] / 2 , *proutput[k] + seq.min_interval[k + 1] / 2);
                  break;
                }
              }
            }
          }

          norm += forward[0][j];
        }

        else {
          forward[0][j] = 0.;
        }
      }
      break;
    }
    }

    if (norm > 0.) {
      for (j = 1;j < nb_row;j++) {
        forward[0][j] /= norm;
      }

      seq_likelihood += log(norm);
    }

    else {
      seq_likelihood = D_INF;
    }

    if (seq_likelihood != D_INF) {
      for (j = 1;j < nb_row;j++) {
        forward_state_entropy[0][j] = 0.;
      }

      for (j = 1;j < seq.length[i];j++) {
        for (k = 0;k < nb_output_process;k++) {
          switch (seq.type[k + 1]) {
          case INT_VALUE :
            pioutput[k]++;
            break;
          case REAL_VALUE :
            proutput[k]++;
            break;
          }
        }
        norm = 0.;

        for (k = 1;k < nb_row;k++) {
          forward[j][k] = 0.;
          for (m = 0;m < nb_memory[k];m++) {
            transition_predicted[m] = transition[previous[k][m]][state[k][0]] * forward[j - 1][previous[k][m]];
            forward[j][k] += transition_predicted[m];

//            forward[j][k] += transition[previous[k][m]][state[k][0]] * forward[j - 1][previous[k][m]];
          }
          predicted[j][k] = forward[j][k];

          forward_state_entropy[j][k] = 0.;
          if (predicted[j][k] > 0.) {
            for (m = 0;m < nb_memory[k];m++) {
              if (transition_predicted[m] > 0.) {
                buff = transition_predicted[m] / predicted[j][k];
                forward_state_entropy[j][k] += buff * (forward_state_entropy[j - 1][previous[k][m]] - log(buff));
              }
            }

            if (forward_state_entropy[j][k] < 0.) {
              forward_state_entropy[j][k] = 0.;
            }
          }

          for (m = 0;m < nb_output_process;m++) {
            if (categorical_process[m]) {
              forward[j][k] *= categorical_process[m]->observation[state[k][0]]->mass[*pioutput[m]];
            }

            else if (discrete_parametric_process[m]) {
              forward[j][k] *= discrete_parametric_process[m]->observation[state[k][0]]->mass[*pioutput[m]];
            }

            else {
              if (((continuous_parametric_process[m]->ident == GAMMA) ||
                  (continuous_parametric_process[m]->ident == ZERO_INFLATED_GAMMA)) && (seq.min_value[m + 1] < seq.min_interval[m + 1] / 2)) {
                switch (seq.type[m + 1]) {
                case INT_VALUE :
                  forward[j][k] *= continuous_parametric_process[m]->observation[state[k][0]]->mass_computation(*pioutput[m] , *pioutput[m] + seq.min_interval[m + 1]);
                  break;
                case REAL_VALUE :
                  forward[j][k] *= continuous_parametric_process[m]->observation[state[k][0]]->mass_computation(*proutput[m] , *proutput[m] + seq.min_interval[m + 1]);
                  break;
                }
              }

              else {
                switch (seq.type[m + 1]) {
                case INT_VALUE :
                  forward[j][k] *= continuous_parametric_process[m]->observation[state[k][0]]->mass_computation(*pioutput[m] - seq.min_interval[m + 1] / 2 , *pioutput[m] + seq.min_interval[m + 1] / 2);
                  break;
                case REAL_VALUE :
                  forward[j][k] *= continuous_parametric_process[m]->observation[state[k][0]]->mass_computation(*proutput[m] - seq.min_interval[m + 1] / 2 , *proutput[m] + seq.min_interval[m + 1] / 2);
                  break;
                }
              }
            }
          }

          norm += forward[j][k];
        }

        if (norm > 0.) {
          for (k = 1;k < nb_row;k++) {
            forward[j][k] /= norm;
          }
          seq_likelihood += log(norm);
        }

        else {
          seq_likelihood = D_INF;
          break;
        }
      }

      seq.entropy[i] = 0.;
      j = seq.length[i] - 1;
      for (k = 1;k < nb_row;k++) {
        if (forward[j][k] > 0.) {
          seq.entropy[i] += forward[j][k] * (forward_state_entropy[j][k] - log(forward[j][k]));
        }
      }
      seq.sample_entropy += seq.entropy[i];

#     ifdef DEBUG
      cout << "\n";
      for (j = 0;j < seq.length[i];j++) {
        cout << j << " |";
        for (k = 1;k < nb_row;k++) {
          cout << " " << forward_state_entropy[j][k];
        }
        cout << endl;
      }
#     endif

    }

    // recurrence "backward"

    if (seq_likelihood != D_INF) {

#     ifdef MESSAGE
      entropy = 0.;

      for (j = 1;j < nb_row;j++) {
        for (k = 0;k < nb_state;k++) {
          transition_entropy[j][k] = 0.;
        }
      }

      j = seq.length[i] - 1;
      for (k = 1;k < nb_row;k++) {
        backward[j][k] = forward[j][k];

        if (backward[j][k] > 0.) {
          for (m = 0;m < nb_output_process;m++) {
            if (categorical_process[m]) {
              if (categorical_process[m]->observation[state[k][0]]->mass[*pioutput[m]] > 0.) {
                entropy -= backward[j][k] * log(categorical_process[m]->observation[state[k][0]]->mass[*pioutput[m]]);
              }
            }

            else if (discrete_parametric_process[m]) {
              if (discrete_parametric_process[m]->observation[state[k][0]]->mass[*pioutput[m]] > 0.) {
                entropy -= backward[j][k] * log(discrete_parametric_process[m]->observation[state[k][0]]->mass[*pioutput[m]]);
              }
            }

            else {
              if (((continuous_parametric_process[m]->ident == GAMMA) ||
                  (continuous_parametric_process[m]->ident == ZERO_INFLATED_GAMMA)) && (seq.min_value[m + 1] < seq.min_interval[m + 1] / 2)) {
                switch (seq.type[m + 1]) {
                case INT_VALUE :
                  entropy -= backward[j][k] * log(continuous_parametric_process[m]->observation[state[k][0]]->mass_computation(*pioutput[m] , *pioutput[m] + seq.min_interval[m + 1]));
                  break;
                case REAL_VALUE :
                  entropy -= backward[j][k] * log(continuous_parametric_process[m]->observation[state[k][0]]->mass_computation(*proutput[m] , *proutput[m] + seq.min_interval[m + 1]));
                  break;
                }
              }

              else {
                switch (seq.type[m + 1]) {
                case INT_VALUE :
                  entropy -= backward[j][k] * log(continuous_parametric_process[m]->observation[state[k][0]]->mass_computation(*pioutput[m] - seq.min_interval[m + 1] / 2 , *pioutput[m] + seq.min_interval[m + 1] / 2));
                  break;
                case REAL_VALUE :
                  entropy -= backward[j][k] * log(continuous_parametric_process[m]->observation[state[k][0]]->mass_computation(*proutput[m] - seq.min_interval[m + 1] / 2 , *proutput[m] + seq.min_interval[m + 1] / 2));
                  break;
                }
              }
            }
          }
        }
      }

      for (j = seq.length[i] - 2;j >= 0;j--) {
        for (k = 0;k < nb_output_process;k++) {
          switch (seq.type[k + 1]) {
          case INT_VALUE :
            pioutput[k]--;
            break;
          case REAL_VALUE :
            proutput[k]--;
            break;
          }
        }

        for (k = 1;k < nb_row;k++) {
          if (predicted[j + 1][k] > 0.) {
            auxiliary[k] = backward[j + 1][k] / predicted[j + 1][k];
          }
          else {
            auxiliary[k] = 0.;
          }
        }

        for (k = 1;k < nb_row;k++) {
          backward[j][k] = 0.;

          if (next[k]) {
            for (m = 0;m < nb_state;m++) {
              buff = auxiliary[next[k][m]] * transition[k][m] * forward[j][k];
              backward[j][k] += buff;
              transition_entropy[k][m] += buff;
            }

            if (backward[j][k] > 0.) {
              for (m = 0;m < nb_output_process;m++) {
                if (categorical_process[m]) {
                  if (categorical_process[m]->observation[state[k][0]]->mass[*pioutput[m]] > 0.) {
                    entropy -= backward[j][k] * log(categorical_process[m]->observation[state[k][0]]->mass[*pioutput[m]]);
                  }
                }

                else if (discrete_parametric_process[m]) {
                  if (discrete_parametric_process[m]->observation[state[k][0]]->mass[*pioutput[m]] > 0.) {
                    entropy -= backward[j][k] * log(discrete_parametric_process[m]->observation[state[k][0]]->mass[*pioutput[m]]);
                  }
                }

                else {
                  if (((continuous_parametric_process[m]->ident == GAMMA) ||
                      (continuous_parametric_process[m]->ident == ZERO_INFLATED_GAMMA)) && (seq.min_value[m + 1] < seq.min_interval[m + 1] / 2)) {
                    switch (seq.type[m + 1]) {
                    case INT_VALUE :
                      entropy -= backward[j][k] * log(continuous_parametric_process[m]->observation[state[k][0]]->mass_computation(*pioutput[m] , *pioutput[m] + seq.min_interval[m + 1]));
                      break;
                    case REAL_VALUE :
                      entropy -= backward[j][k] * log(continuous_parametric_process[m]->observation[state[k][0]]->mass_computation(*proutput[m] , *proutput[m] + seq.min_interval[m + 1]));
                      break;
                    }
                  }

                  else {
                    switch (seq.type[m + 1]) {
                    case INT_VALUE :
                      entropy -= backward[j][k] * log(continuous_parametric_process[m]->observation[state[k][0]]->mass_computation(*pioutput[m] - seq.min_interval[m + 1] / 2 , *pioutput[m] + seq.min_interval[m + 1] / 2));
                      break;
                    case REAL_VALUE :
                      entropy -= backward[j][k] * log(continuous_parametric_process[m]->observation[state[k][0]]->mass_computation(*proutput[m] - seq.min_interval[m + 1] / 2 , *proutput[m] + seq.min_interval[m + 1] / 2));
                      break;
                    }
                  }
                }
              }
            }
          }
        }
      }

      for (j = 1;j < nb_row;j++) {
        switch (type) {

        case 'o' : {
          if ((order[j] == 1) && (initial[state[j][0]] > 0.)) {
            entropy -= backward[0][j] * log(initial[state[j][0]]);
          }
          break;
        }

        case 'e' : {
          if ((!child[j]) && (initial[j] > 0.)) {
            entropy -= backward[0][j] * log(initial[j]);
          }
          break;
        }
        }
      }

      for (j = 1;j < nb_row;j++) {
        for (k = 0;k < nb_state;k++) {
          if (transition[j][k] > 0.) {
            entropy -= transition_entropy[j][k] * log(transition[j][k]);
          }
        }
      }

      entropy += seq_likelihood;

      if ((entropy < seq.entropy[i] - DOUBLE_ERROR) || (entropy > seq.entropy[i] + DOUBLE_ERROR)) {
        cout << "\nERROR: " << i << " " << seq.entropy[i] << " " << entropy << endl;
      }
#     endif

      // calcul du nombre de sequences d'etats possibles

      for (j = 0;j < nb_output_process;j++) {
        switch (seq.type[j + 1]) {
        case INT_VALUE :
          pioutput[j] = seq.int_sequence[i][j + 1];
          break;
        case REAL_VALUE :
          proutput[j] = seq.real_sequence[i][j + 1];
          break;
        }
      }

      // recurrence "forward"

      switch (type) {

      case 'o' : {
        for (j = 1;j < nb_row;j++) {
          if (order[j] == 1) {
            forward[0][j] = initial[state[j][0]];

            for (k = 0;k < nb_output_process;k++) {
              if (categorical_process[k]) {
                forward[0][j] *= categorical_process[k]->observation[state[j][0]]->mass[*pioutput[k]];
              }

              else if (discrete_parametric_process[k]) {
                forward[0][j] *= discrete_parametric_process[k]->observation[state[j][0]]->mass[*pioutput[k]];
              }

              else {
                if (((continuous_parametric_process[k]->ident == GAMMA) ||
                    (continuous_parametric_process[k]->ident == ZERO_INFLATED_GAMMA)) && (seq.min_value[k + 1] < seq.min_interval[k + 1] / 2)) {
                  switch (seq.type[k + 1]) {
                  case INT_VALUE :
                    forward[0][j] *= continuous_parametric_process[k]->observation[state[j][0]]->mass_computation(*pioutput[k] , *pioutput[k] + seq.min_interval[k + 1]);
                    break;
                  case REAL_VALUE :
                    forward[0][j] *= continuous_parametric_process[k]->observation[state[j][0]]->mass_computation(*proutput[k] , *proutput[k] + seq.min_interval[k + 1]);
                    break;
                  }
                }

                else {
                  switch (seq.type[k + 1]) {
                  case INT_VALUE :
                    forward[0][j] *= continuous_parametric_process[k]->observation[state[j][0]]->mass_computation(*pioutput[k] - seq.min_interval[k + 1] / 2 , *pioutput[k] + seq.min_interval[k + 1] / 2);
                    break;
                  case REAL_VALUE :
                    forward[0][j] *= continuous_parametric_process[k]->observation[state[j][0]]->mass_computation(*proutput[k] - seq.min_interval[k + 1] / 2 , *proutput[k] + seq.min_interval[k + 1] / 2);
                    break;
                  }
                }
              }
            }

            if (forward[0][j] > 0.) {
              forward[0][j] = 1.;
            }
          }

          else {
            forward[0][j] = 0.;
          }
        }
        break;
      }

      case 'e' : {
        for (j = 1;j < nb_row;j++) {
          if (!child[j]) {
            forward[0][j] = initial[j];

            for (k = 0;k < nb_output_process;k++) {
              if (categorical_process[k]) {
                forward[0][j] *= categorical_process[k]->observation[state[j][0]]->mass[*pioutput[k]];
              }

              else if (discrete_parametric_process[k]) {
                forward[0][j] *= discrete_parametric_process[k]->observation[state[j][0]]->mass[*pioutput[k]];
              }

              else {
                if (((continuous_parametric_process[k]->ident == GAMMA) ||
                    (continuous_parametric_process[k]->ident == ZERO_INFLATED_GAMMA)) && (seq.min_value[k + 1] < seq.min_interval[k + 1] / 2)) {
                  switch (seq.type[k + 1]) {
                  case INT_VALUE :
                    forward[0][j] *= continuous_parametric_process[k]->observation[state[j][0]]->mass_computation(*pioutput[k] , *pioutput[k] + seq.min_interval[k + 1]);
                    break;
                  case REAL_VALUE :
                    forward[0][j] *= continuous_parametric_process[k]->observation[state[j][0]]->mass_computation(*proutput[k] , *proutput[k] + seq.min_interval[k + 1]);
                    break;
                  }
                }

                else {
                  switch (seq.type[k + 1]) {
                  case INT_VALUE :
                    forward[0][j] *= continuous_parametric_process[k]->observation[state[j][0]]->mass_computation(*pioutput[k] - seq.min_interval[k + 1] / 2 , *pioutput[k] + seq.min_interval[k + 1] / 2);
                    break;
                  case REAL_VALUE :
                    forward[0][j] *= continuous_parametric_process[k]->observation[state[j][0]]->mass_computation(*proutput[k] - seq.min_interval[k + 1] / 2 , *proutput[k] + seq.min_interval[k + 1] / 2);
                    break;
                  }
                }
              }
            }

            if (forward[0][j] > 0.) {
              forward[0][j] = 1.;
            }
          }

          else {
            forward[0][j] = 0.;
          }
        }
        break;
      }
      }

      for (j = 1;j < seq.length[i];j++) {
        for (k = 0;k < nb_output_process;k++) {
          switch (seq.type[k + 1]) {
          case INT_VALUE :
            pioutput[k]++;
            break;
          case REAL_VALUE :
            proutput[k]++;
            break;
          }
        }

        for (k = 1;k < nb_row;k++) {
          forward[j][k] = 0.;
          for (m = 0;m < nb_memory[k];m++) {
            if (transition[previous[k][m]][state[k][0]] > 0.) {
              forward[j][k] += forward[j - 1][previous[k][m]];
            }
          }

          observation = 1.;
          for (m = 0;m < nb_output_process;m++) {
            if (categorical_process[m]) {
              observation *= categorical_process[m]->observation[state[k][0]]->mass[*pioutput[m]];
            }

            else if (discrete_parametric_process[m]) {
              observation *= discrete_parametric_process[m]->observation[state[k][0]]->mass[*pioutput[m]];
            }

            else {
              if (((continuous_parametric_process[m]->ident == GAMMA) ||
                  (continuous_parametric_process[m]->ident == ZERO_INFLATED_GAMMA)) && (seq.min_value[m + 1] < seq.min_interval[m + 1] / 2)) {
                switch (seq.type[m + 1]) {
                case INT_VALUE :
                  observation *= continuous_parametric_process[m]->observation[state[k][0]]->mass_computation(*pioutput[m] , *pioutput[m] + seq.min_interval[m + 1]);
                  break;
                case REAL_VALUE :
                  observation *= continuous_parametric_process[m]->observation[state[k][0]]->mass_computation(*proutput[m] , *proutput[m] + seq.min_interval[m + 1]);
                  break;
                }
              }

              else {
                switch (seq.type[m + 1]) {
                case INT_VALUE :
                  observation *= continuous_parametric_process[m]->observation[state[k][0]]->mass_computation(*pioutput[m] - seq.min_interval[m + 1] / 2 , *pioutput[m] + seq.min_interval[m + 1] / 2);
                  break;
                case REAL_VALUE :
                  observation *= continuous_parametric_process[m]->observation[state[k][0]]->mass_computation(*proutput[m] - seq.min_interval[m + 1] / 2 , *proutput[m] + seq.min_interval[m + 1] / 2);
                  break;
                }
              }
            }
          }

          if (observation == 0.) {
            forward[j][k] = 0.;
          }
        }
      }

      seq.nb_state_sequence[i] = 0.;
      j = seq.length[i] - 1;
      for (k = 1;k < nb_row;k++) {
        seq.nb_state_sequence[i] += forward[j][k];
      }
    }
  }

  for (i = 0;i < seq.max_length;i++) {
    delete [] forward[i];
  }
  delete [] forward;

  for (i = 0;i < seq.max_length;i++) {
    delete [] predicted[i];
  }
  delete [] predicted;

  delete [] transition_predicted;

  for (i = 0;i < seq.max_length;i++) {
    delete [] forward_state_entropy[i];
  }
  delete [] forward_state_entropy;

# ifdef MESSAGE
  for (i = 0;i < seq.max_length;i++) {
    delete [] backward[i];
  }
  delete [] backward;

  delete [] auxiliary;

  for (i = 1;i < nb_row;i++) {
    delete [] transition_entropy[i];
  }
  delete [] transition_entropy;
# endif

  delete [] pioutput;
  delete [] proutput;
}


/*--------------------------------------------------------------
 *
 *  Ecriture des profils d'etats et d'entropie.
 *
 *  arguments : stream, indice de la sequence, nombre d'etats,
 *              pointeur sur les profils d'etats, 
 *              pointeurs sur les profils d'entropie.
 *
 *--------------------------------------------------------------*/

ostream& Sequences::profile_ascii_print(ostream &os , int index , int nb_state ,
                                        double **profiles , double *begin_conditional_entropy ,
                                        double *marginal_entropy , double *begin_partial_entropy ,
                                        double *end_conditional_entropy , double *end_partial_entropy) const

{
  register int i , j;
  int buff , *width;
  long old_adjust;


  old_adjust = os.flags(ios::adjustfield);

  // calcul des largeurs des colonnes

  width = new int[nb_variable + 8];

  for (i = 0;i < nb_variable;i++) {
    if (type[i] != REAL_VALUE) {
      width[i] = column_width((int)max_value[i]);
    }
    else {
      width[i] = column_width(length[index] , real_sequence[index][i]);
    }

    if (i > 0) {
      width[i] += ASCII_SPACE;
    }
  }

  if (index_parameter) {
    width[nb_variable] = column_width(index_parameter_distribution->nb_value - 1) + ASCII_SPACE;
  }
  else {
    width[nb_variable] = column_width(max_length) + ASCII_SPACE;
  }

  width[nb_variable + 1] = 0;
  for (i = 0;i < length[index];i++) {
    buff = column_width(nb_state , profiles[i]);
    if (buff > width[nb_variable + 1]) {
      width[nb_variable + 1] = buff;
    }
  }
  width[nb_variable + 1] += ASCII_SPACE;

  width[nb_variable + 2] = column_width(length[index] , begin_conditional_entropy) + ASCII_SPACE;
  width[nb_variable + 3] = column_width(length[index] , marginal_entropy) + ASCII_SPACE;
  width[nb_variable + 4] = column_width(length[index] , begin_partial_entropy) + ASCII_SPACE;

  width[nb_variable + 5] = column_width(nb_sequence);

  if ((end_conditional_entropy) && (end_partial_entropy)) {
    width[nb_variable + 6] = column_width(length[index] , end_conditional_entropy) + ASCII_SPACE;
    width[nb_variable + 7] = column_width(length[index] , end_partial_entropy) + ASCII_SPACE;
  }

  os << SEQ_label[SEQL_OPTIMAL] << " " << STAT_label[STATL_STATE];
  for (i = 1;i < nb_variable;i++) {
    os << " | " << STAT_label[STATL_VARIABLE] << " " << i;
  }
  if (index_parameter_type == TIME) {
    os << " | " << SEQ_label[SEQL_TIME];
  }
  else {
    os << " | " << SEQ_label[SEQL_INDEX];
  }
  for (i = 0;i < nb_state;i++) {
    os << " | " << STAT_label[STATL_STATE] << " " << i;
  }
  os << " | " << SEQ_label[SEQL_CONDITIONAL_ENTROPY] << " | " << SEQ_label[SEQL_MARGINAL_ENTROPY]
     << " | " << SEQ_label[SEQL_PARTIAL_STATE_SEQUENCE_ENTROPY] << endl;

  for (i = 0;i < length[index];i++) {
    os.setf(ios::right , ios::adjustfield);
    for (j = 0;j < nb_variable;j++) {
      if (type[j] != REAL_VALUE) {
        os << setw(width[j]) << int_sequence[index][j][i];
      }
      else {
        os << setw(width[j]) << real_sequence[index][j][i];
      }
    }
    os << setw(width[nb_variable]) << (index_parameter ? index_parameter[index][i] : i) << "  ";

    os.setf(ios::left , ios::adjustfield);
    for (j = 0;j < nb_state;j++) {
      os << setw(width[nb_variable + 1]) << profiles[i][j];
    }

    if ((end_conditional_entropy) && (end_partial_entropy)) {
      os << setw(width[nb_variable + 2]) << begin_conditional_entropy[i];
      os << setw(width[nb_variable + 6]) << end_conditional_entropy[i];
      os << setw(width[nb_variable + 3]) << marginal_entropy[i];
      os << setw(width[nb_variable + 4]) << begin_partial_entropy[i];
      os << setw(width[nb_variable + 7]) << end_partial_entropy[i];
    }
    else {
      os << setw(width[nb_variable + 2]) << begin_conditional_entropy[i];
      os << setw(width[nb_variable + 3]) << marginal_entropy[i];
      os << setw(width[nb_variable + 4]) << begin_partial_entropy[i];
    }

    if (i == 0) {
      os.setf(ios::right , ios::adjustfield);
      os << setw(width[nb_variable + 5]) << identifier[index];
    }
    os << endl;
  }

  delete [] width;

  os.setf((FMTFLAGS)old_adjust , ios::adjustfield);

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture des profils d'etats et d'entropie au format tableur.
 *
 *  arguments : stream, indice de la sequence, nombre d'etats,
 *              pointeur sur les profils d'etats, 
 *              pointeurs sur les profils d'entropie.
 *
 *--------------------------------------------------------------*/

ostream& Sequences::profile_spreadsheet_print(ostream &os , int index , int nb_state ,
                                              double **profiles , double *begin_conditional_entropy ,
                                              double *marginal_entropy , double *begin_partial_entropy ,
                                              double *end_conditional_entropy , double *end_partial_entropy) const

{
  register int i , j;


  os << SEQ_label[SEQL_OPTIMAL] << " " << STAT_label[STATL_STATE];
  for (i = 1;i < nb_variable;i++) {
    os << "\t" << STAT_label[STATL_VARIABLE] << " " << i;
  }
  if (index_parameter_type == TIME) {
    os << "\t" << SEQ_label[SEQL_TIME];
  }
  else {
    os << "\t" << SEQ_label[SEQL_INDEX];
  }
  for (i = 0;i < nb_state;i++) {
    os << "\t" << STAT_label[STATL_STATE] << " " << i;
  }
  os << "\t" << SEQ_label[SEQL_CONDITIONAL_ENTROPY] << "\t" << SEQ_label[SEQL_MARGINAL_ENTROPY]
     << "\t" << SEQ_label[SEQL_PARTIAL_STATE_SEQUENCE_ENTROPY] << endl;

  for (i = 0;i < length[index];i++) {
    for (j = 0;j < nb_variable;j++) {
      if (type[j] != REAL_VALUE) {
        os << int_sequence[index][j][i] << "\t";
      }
      else {
        os << real_sequence[index][j][i] << "\t";
      }
    }
    os << (index_parameter ? index_parameter[index][i] : i);

    for (j = 0;j < nb_state;j++) {
      os << "\t" << profiles[i][j];
    }

    if ((end_conditional_entropy) && (end_partial_entropy)) {
      os << "\t" << begin_conditional_entropy[i] << "\t" << end_conditional_entropy[i] << "\t" << marginal_entropy[i]
         << "\t" << begin_partial_entropy[i] << "\t" << end_partial_entropy[i];
    }
    else {
      os << "\t" << begin_conditional_entropy[i] << "\t" << marginal_entropy[i]
         << "\t" << begin_partial_entropy[i];
    }

    if (i == 0) {
      os << "\t" << identifier[index];
    }
    os << endl;
  }

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture des profils d'etats et d'entropie au format Gnuplot.
 *
 *  arguments : stream, indice de la sequence, nombre d'etats,
 *              pointeur sur les profils d'etat,
 *              pointeurs sur les profils d'entropies.
 *
 *--------------------------------------------------------------*/

ostream& Sequences::profile_plot_print(ostream &os , int index , int nb_state ,
                                       double **profiles , double *begin_conditional_entropy ,
                                       double *marginal_entropy , double *begin_partial_entropy ,
                                       double *end_conditional_entropy , double *end_partial_entropy) const

{
  register int i , j;


  for (i = 0;i < length[index];i++) {
    if (index_parameter) {
      os << index_parameter[index][i] << " ";
    }

    for (j = 0;j < nb_state;j++) {
      os << profiles[i][j] << " ";
    }

    if ((end_conditional_entropy) && (end_partial_entropy)) {
      os << begin_conditional_entropy[i] << " " << end_conditional_entropy[i] << " "
         << marginal_entropy[i] << " " << begin_partial_entropy[i] << " "
         << end_partial_entropy[i] << endl;
    }
    else {
      os << begin_conditional_entropy[i] << " " << marginal_entropy[i] << " "
         << begin_partial_entropy[i] << endl;
    }
  }

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture des profils d'etats au format "plotable".
 *
 *  arguments : reference sur un objet MultiPlot, indice de la sequence,
 *              nombre d'etats, pointeur sur les profils d'etat.
 *
 *--------------------------------------------------------------*/

void Sequences::profile_plotable_write(MultiPlot &plot , int index , int nb_state ,
                                       double **profiles) const

{
  register int i , j;


  plot.resize(nb_state);

  if (index_parameter) {
    for (i = 0;i < length[index];i++) {
      for (j = 0;j < nb_state;j++) {
        plot[j].add_point(index_parameter[index][i] , profiles[i][j]);
      }
    }
  }

  else {
    for (i = 0;i < length[index];i++) {
      for (j = 0;j < nb_state;j++) {
        plot[j].add_point(i , profiles[i][j]);
      }
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Ecriture des profils d'entropie au format "plotable".
 *
 *  arguments : reference sur un objet MultiPlot, indice de la sequence,
 *              pointeurs sur les profils d'entropies.
 *
 *--------------------------------------------------------------*/

void Sequences::entropy_profile_plotable_write(MultiPlot &plot , int index ,
                                               double *begin_entropy , double *end_entropy ,
                                               double *marginal_entropy) const

{
  register int i , j;
  int nb_plot;


  nb_plot = 1;
  if (end_entropy) {
    nb_plot++;
  }
  if (marginal_entropy) {
    nb_plot++;
  }
  plot.resize(nb_plot);

  if (index_parameter) {
    for (i = 0;i < length[index];i++) {
      plot[0].add_point(index_parameter[index][i] , begin_entropy[i]);
    }

    if (end_entropy) {
      for (i = 0;i < length[index];i++) {
        plot[1].add_point(index_parameter[index][i] , end_entropy[i]);
      }
    }

    if (marginal_entropy) {
      i = (end_entropy ? 2 : 1);
      for (j = 0;j < length[index];j++) {
        plot[i].add_point(index_parameter[index][j] , marginal_entropy[j]);
      }
    }
  }

  else {
    for (i = 0;i < length[index];i++) {
      plot[0].add_point(i , begin_entropy[i]);
    }

    if (end_entropy) {
      for (i = 0;i < length[index];i++) {
        plot[1].add_point(i , end_entropy[i]);
      }
    }

    if (marginal_entropy) {
      i = (end_entropy ? 2 : 1);
      for (j = 0;j < length[index];j++) {
        plot[i].add_point(j , marginal_entropy[j]);
      }
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Calcul des profils d'etats et d'entropie par l'algorithme forward-backward.
 *
 *  arguments : reference sur un objet MarkovianSequences, indice de la sequence,
 *              stream, pointeur sur un objet MultiPlotSet, format de sortie ('a' : ASCII,
 *              's' : Spreadsheet, 'g' : Gnuplot, 'p' : plotable), references sur
 *              l'entropie marginale maximum et l'entropie (pour la visualisation).
 *
 *--------------------------------------------------------------*/

double HiddenVariableOrderMarkov::forward_backward(const MarkovianSequences &seq , int index ,
                                                   ostream *os , MultiPlotSet *plot_set ,
                                                   char format , double &max_marginal_entropy ,
                                                   double &entropy1) const

{
  register int i , j , k;
  int *pstate , **pioutput;
  double seq_likelihood , state_seq_likelihood , **forward , norm , **predicted ,
         entropy2 , buff , **backward , *auxiliary , backward_max , **state_backward ,
         *transition_predicted , **forward_state_entropy , **transition_entropy ,
         *begin_partial_entropy , *begin_conditional_entropy , *end_backward ,
         *end_partial_entropy , *end_conditional_entropy , *marginal_entropy , **proutput;
//  double **backward_state_entropy;


  // initialisations

  forward = new double*[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    forward[i] = new double[nb_row];
  }

  predicted = new double*[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    predicted[i] = new double[nb_row];
  }

  backward = new double*[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    backward[i] = new double[nb_row];
  }

  auxiliary = new double[nb_row];

  state_backward = new double*[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    state_backward[i] = new double[nb_state];
  }

  transition_predicted = new double[nb_row];

  forward_state_entropy = new double*[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    forward_state_entropy[i] = new double[nb_row];
  }

  transition_entropy = new double*[nb_row];
  for (i = 1;i < nb_row;i++) {
    transition_entropy[i] = new double[nb_state];
  }

  begin_partial_entropy = new double[seq.length[index]];
  begin_conditional_entropy = new double[seq.length[index]];

/*  backward_state_entropy = new double*[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    backward_state_entropy[i] = new double[nb_row];
  } */

  end_backward = new double[nb_row];
  end_partial_entropy = new double[seq.length[index]];
  end_conditional_entropy = new double[seq.length[index]];

  marginal_entropy = new double[seq.length[index]];

  pioutput = new int*[nb_output_process];
  proutput = new double*[nb_output_process];

  for (i = 0;i < nb_output_process;i++) {
    switch (seq.type[i + 1]) {
    case INT_VALUE :
      pioutput[i] = seq.int_sequence[index][i + 1];
      break;
    case REAL_VALUE :
      proutput[i] = seq.real_sequence[index][i + 1];
      break;
    }
  }

  // recurrence "forward"

  seq_likelihood = 0.;
  norm = 0.;

  switch (type) {

  case 'o' : {
    for (i = 1;i < nb_row;i++) {
      if (order[i] == 1) {
        forward[0][i] = initial[state[i][0]];

        for (j = 0;j < nb_output_process;j++) {
          if (categorical_process[j]) {
            forward[0][i] *= categorical_process[j]->observation[state[i][0]]->mass[*pioutput[j]];
          }

          else if (discrete_parametric_process[j]) {
            forward[0][i] *= discrete_parametric_process[j]->observation[state[i][0]]->mass[*pioutput[j]];
          }

          else {
            if (((continuous_parametric_process[j]->ident == GAMMA) ||
                (continuous_parametric_process[j]->ident == ZERO_INFLATED_GAMMA)) && (seq.min_value[j + 1] < seq.min_interval[j + 1] / 2)) {
              switch (seq.type[j + 1]) {
              case INT_VALUE :
                forward[0][i] *= continuous_parametric_process[j]->observation[state[i][0]]->mass_computation(*pioutput[j] , *pioutput[j] + seq.min_interval[j + 1]);
                break;
              case REAL_VALUE :
                forward[0][i] *= continuous_parametric_process[j]->observation[state[i][0]]->mass_computation(*proutput[j] , *proutput[j] + seq.min_interval[j + 1]);
                break;
              }
            }

            else {
              switch (seq.type[j + 1]) {
              case INT_VALUE :
                forward[0][i] *= continuous_parametric_process[j]->observation[state[i][0]]->mass_computation(*pioutput[j] - seq.min_interval[j + 1] / 2 , *pioutput[j] + seq.min_interval[j + 1] / 2);
                break;
              case REAL_VALUE :
                forward[0][i] *= continuous_parametric_process[j]->observation[state[i][0]]->mass_computation(*proutput[j] - seq.min_interval[j + 1] / 2 , *proutput[j] + seq.min_interval[j + 1] / 2);
                break;
              }
            }
          }
        }

        norm += forward[0][i];
      }

      else {
        forward[0][i] = 0.;
      }
    }
    break;
  }

  case 'e' : {
    for (i = 1;i < nb_row;i++) {
      if (!child[i]) {
        forward[0][i] = initial[i];

        for (j = 0;j < nb_output_process;j++) {
          if (categorical_process[j]) {
            forward[0][i] *= categorical_process[j]->observation[state[i][0]]->mass[*pioutput[j]];
          }

          else if (discrete_parametric_process[j]) {
            forward[0][i] *= discrete_parametric_process[j]->observation[state[i][0]]->mass[*pioutput[j]];
          }

          else {
            if (((continuous_parametric_process[j]->ident == GAMMA) ||
                (continuous_parametric_process[j]->ident == ZERO_INFLATED_GAMMA)) && (seq.min_value[j + 1] < seq.min_interval[j + 1] / 2)) {
              switch (seq.type[j + 1]) {
              case INT_VALUE :
                forward[0][i] *= continuous_parametric_process[j]->observation[state[i][0]]->mass_computation(*pioutput[j] , *pioutput[j] + seq.min_interval[j + 1]);
                break;
              case REAL_VALUE :
                forward[0][i] *= continuous_parametric_process[j]->observation[state[i][0]]->mass_computation(*proutput[j] , *proutput[j] + seq.min_interval[j + 1]);
                break;
              }
            }

            else {
              switch (seq.type[j + 1]) {
              case INT_VALUE :
                forward[0][i] *= continuous_parametric_process[j]->observation[state[i][0]]->mass_computation(*pioutput[j] - seq.min_interval[j + 1] / 2 , *pioutput[j] + seq.min_interval[j + 1] / 2);
                break;
              case REAL_VALUE :
                forward[0][i] *= continuous_parametric_process[j]->observation[state[i][0]]->mass_computation(*proutput[j] - seq.min_interval[j + 1] / 2 , *proutput[j] + seq.min_interval[j + 1] / 2);
                break;
              }
            }
          }
        }

        norm += forward[0][i];
      }

      else {
        forward[0][i] = 0.;
      }
    }
    break;
  }
  }

  if (norm > 0.) {
    for (i = 1;i < nb_row;i++) {
      forward[0][i] /= norm;
    }

    seq_likelihood += log(norm);
  }

  else {
    seq_likelihood = D_INF;
  }

  if (seq_likelihood != D_INF) {
    for (i = 1;i < nb_row;i++) {
      forward_state_entropy[0][i] = 0.;
    }

    for (i = 1;i < seq.length[index];i++) {
      for (j = 0;j < nb_output_process;j++) {
        switch (seq.type[j + 1]) {
        case INT_VALUE :
          pioutput[j]++;
          break;
        case REAL_VALUE :
          proutput[j]++;
          break;
        }
      }
      norm = 0.;

      for (j = 1;j < nb_row;j++) {
        forward[i][j] = 0.;
        for (k = 0;k < nb_memory[j];k++) {
          transition_predicted[k] = transition[previous[j][k]][state[j][0]] * forward[i - 1][previous[j][k]];
          forward[i][j] += transition_predicted[k];

//          forward[i][j] += transition[previous[j][k]][state[j][0]] * forward[i - 1][previous[j][k]];
        }
        predicted[i][j] = forward[i][j];

        forward_state_entropy[i][j] = 0.;
        if (predicted[i][j] > 0.) {
          for (k = 0;k < nb_memory[j];k++) {
            if (transition_predicted[k] > 0.) {
              buff = transition_predicted[k] / predicted[i][j];
              forward_state_entropy[i][j] += buff * (forward_state_entropy[i - 1][previous[j][k]] - log(buff));
            }
          }

          if (forward_state_entropy[i][j] < 0.) {
            forward_state_entropy[i][j] = 0.;
          }
        }

        for (k = 0;k < nb_output_process;k++) {
          if (categorical_process[k]) {
            forward[i][j] *= categorical_process[k]->observation[state[j][0]]->mass[*pioutput[k]];
          }

          else if (discrete_parametric_process[k]) {
            forward[i][j] *= discrete_parametric_process[k]->observation[state[j][0]]->mass[*pioutput[k]];
          }

          else {
            if (((continuous_parametric_process[k]->ident == GAMMA) ||
                (continuous_parametric_process[k]->ident == ZERO_INFLATED_GAMMA)) && (seq.min_value[k + 1] < seq.min_interval[k + 1] / 2)) {
              switch (seq.type[k + 1]) {
              case INT_VALUE :
                forward[i][j] *= continuous_parametric_process[k]->observation[state[j][0]]->mass_computation(*pioutput[k] , *pioutput[k] + seq.min_interval[k + 1]);
                break;
              case REAL_VALUE :
                forward[i][j] *= continuous_parametric_process[k]->observation[state[j][0]]->mass_computation(*proutput[k] , *proutput[k] + seq.min_interval[k + 1]);
                break;
              }
            }

            else {
              switch (seq.type[k + 1]) {
              case INT_VALUE :
                forward[i][j] *= continuous_parametric_process[k]->observation[state[j][0]]->mass_computation(*pioutput[k] - seq.min_interval[k + 1] / 2 , *pioutput[k] + seq.min_interval[k + 1] / 2);
                break;
              case REAL_VALUE :
                forward[i][j] *= continuous_parametric_process[k]->observation[state[j][0]]->mass_computation(*proutput[k] - seq.min_interval[k + 1] / 2 , *proutput[k] + seq.min_interval[k + 1] / 2);
                break;
              }
            }
          }
        }

        norm += forward[i][j];
      }

      if (norm > 0.) {
        for (j = 1;j < nb_row;j++) {
          forward[i][j] /= norm;
        }
        seq_likelihood += log(norm);
      }

      else {
        seq_likelihood = D_INF;
        break;
      }
    }

    entropy1 = 0.;
    i = seq.length[index] - 1;
    for (j = 1;j < nb_row;j++) {
      if (forward[i][j] > 0.) {
        entropy1 += forward[i][j] * (forward_state_entropy[i][j] - log(forward[i][j]));
      }
    }

#   ifdef DEBUG
    cout << "\n";
    for (i = 0;i < seq.length[index];i++) {
      cout << i << " |";
      for (j = 1;j < nb_row;j++) {
        cout << " " << forward_state_entropy[i][j];
      }
      cout << endl;
    }
#   endif

  }

  // recurrence "backward"

  if (seq_likelihood != D_INF) {
    entropy2 = 0.;

    for (i = 1;i < nb_row;i++) {
      for (j = 0;j < nb_state;j++) {
        transition_entropy[i][j] = 0.;
      }
    }

    i = seq.length[index] - 1;
    for (j = 1;j < nb_row;j++) {
      backward[i][j] = forward[i][j];

      if (backward[i][j] > 0.) {
        for (k = 0;k < nb_output_process;k++) {
          if (categorical_process[k]) {
            if (categorical_process[k]->observation[state[j][0]]->mass[*pioutput[k]] > 0.) {
              entropy2 -= backward[i][j] * log(categorical_process[k]->observation[state[j][0]]->mass[*pioutput[k]]);
            }
          }

          else if (discrete_parametric_process[k]) {
            if (discrete_parametric_process[k]->observation[state[j][0]]->mass[*pioutput[k]] > 0.) {
              entropy2 -= backward[i][j] * log(discrete_parametric_process[k]->observation[state[j][0]]->mass[*pioutput[k]]);
            }
          }

          else {
            if (((continuous_parametric_process[k]->ident == GAMMA) ||
                (continuous_parametric_process[k]->ident == ZERO_INFLATED_GAMMA)) && (seq.min_value[k + 1] < seq.min_interval[k + 1] / 2)) {
              switch (seq.type[k + 1]) {
              case INT_VALUE :
                entropy2 -= backward[i][j] * log(continuous_parametric_process[k]->observation[state[j][0]]->mass_computation(*pioutput[k] , *pioutput[k] + seq.min_interval[k + 1]));
                break;
              case REAL_VALUE :
                entropy2 -= backward[i][j] * log(continuous_parametric_process[k]->observation[state[j][0]]->mass_computation(*proutput[k] , *proutput[k] + seq.min_interval[k + 1]));
                break;
              }
            }

            else {
              switch (seq.type[k + 1]) {
              case INT_VALUE :
                entropy2 -= backward[i][j] * log(continuous_parametric_process[k]->observation[state[j][0]]->mass_computation(*pioutput[k] - seq.min_interval[k + 1] / 2 , *pioutput[k] + seq.min_interval[k + 1] / 2));
                break;
              case REAL_VALUE :
                entropy2 -= backward[i][j] * log(continuous_parametric_process[k]->observation[state[j][0]]->mass_computation(*proutput[k] - seq.min_interval[k + 1] / 2 , *proutput[k] + seq.min_interval[k + 1] / 2));
                break;
              }
            }
          }
        }
      }

//      backward_state_entropy[i][j] = 0.;
    }

    for (i = seq.length[index] - 2;i >= 0;i--) {
      for (j = 0;j < nb_output_process;j++) {
        switch (seq.type[j + 1]) {
        case INT_VALUE :
          pioutput[j]--;
          break;
        case REAL_VALUE :
          proutput[j]--;
          break;
        }
      }

      for (j = 1;j < nb_row;j++) {
        if (predicted[i + 1][j] > 0.) {
          auxiliary[j] = backward[i + 1][j] / predicted[i + 1][j];
        }
        else {
          auxiliary[j] = 0.;
        }
      }

      for (j = 1;j < nb_row;j++) {
        backward[i][j] = 0.;
//        backward_state_entropy[i][j] = 0.;

        if (next[j]) {
//          norm = 0.;

          for (k = 0;k < nb_state;k++) {
/*            transition_predicted[k] = auxiliary[next[j][k]] * transition[j][k];
            norm += transition_predicted[k]; */

            buff = auxiliary[next[j][k]] * transition[j][k] * forward[i][j];
            backward[i][j] += buff;
            transition_entropy[j][k] += buff;

/*            if (transition[j][k] > 0.) {
              entropy2 -= buff * log(transition[j][k]);
            } */
          }

/*          if (norm > 0.) {
            for (k = 0;k < nb_state;k++) {
              if (transition_predicted[k] > 0.) {
                buff = transition_predicted[k] / norm;
                backward_state_entropy[i][j] += buff * (backward_state_entropy[i + 1][next[j][k]] - log(buff));
              }
            }

            if (backward_state_entropy[i][j] < 0.) {
              backward_state_entropy[i][j] = 0.;
            }
          } */

          if (backward[i][j] > 0.) {
            for (k = 0;k < nb_output_process;k++) {
              if (categorical_process[k]) {
                if (categorical_process[k]->observation[state[j][0]]->mass[*pioutput[k]] > 0.) {
                  entropy2 -= backward[i][j] * log(categorical_process[k]->observation[state[j][0]]->mass[*pioutput[k]]);
                }
              }

              else if (discrete_parametric_process[k]) {
                if (discrete_parametric_process[k]->observation[state[j][0]]->mass[*pioutput[k]] > 0.) {
                  entropy2 -= backward[i][j] * log(discrete_parametric_process[k]->observation[state[j][0]]->mass[*pioutput[k]]);
                }
              }

              else {
                if (((continuous_parametric_process[k]->ident == GAMMA) ||
                    (continuous_parametric_process[k]->ident == ZERO_INFLATED_GAMMA)) && (seq.min_value[k + 1] < seq.min_interval[k + 1] / 2)) {
                  switch (seq.type[k + 1]) {
                  case INT_VALUE :
                    entropy2 -= backward[i][j] * log(continuous_parametric_process[k]->observation[state[j][0]]->mass_computation(*pioutput[k] , *pioutput[k] + seq.min_interval[k + 1]));
                    break;
                  case REAL_VALUE :
                    entropy2 -= backward[i][j] * log(continuous_parametric_process[k]->observation[state[j][0]]->mass_computation(*proutput[k] , *proutput[k] + seq.min_interval[k + 1]));
                    break;
                  }
                }

                else {
                  switch (seq.type[k + 1]) {
                  case INT_VALUE :
                    entropy2 -= backward[i][j] * log(continuous_parametric_process[k]->observation[state[j][0]]->mass_computation(*pioutput[k] - seq.min_interval[k + 1] / 2 , *pioutput[k] + seq.min_interval[k + 1] / 2));
                    break;
                  case REAL_VALUE :
                    entropy2 -= backward[i][j] * log(continuous_parametric_process[k]->observation[state[j][0]]->mass_computation(*proutput[k] - seq.min_interval[k + 1] / 2 , *proutput[k] + seq.min_interval[k + 1] / 2));
                    break;
                  }
                }
              }
            }
          }
        }
      }
    }

    for (i = 1;i < nb_row;i++) {
      switch (type) {

      case 'o' : {
        if ((order[i] == 1) && (initial[state[i][0]] > 0.)) {
          entropy2 -= backward[0][i] * log(initial[state[i][0]]);
        }
        break;
      }

      case 'e' : {
        if ((!child[i]) && (initial[i] > 0.)) {
          entropy2 -= backward[0][i] * log(initial[i]);
        }
        break;
      }
      }
    }

    for (i = 1;i < nb_row;i++) {
      for (j = 0;j < nb_state;j++) {
        if (transition[i][j] > 0.) {
          entropy2 -= transition_entropy[i][j] * log(transition[i][j]);
        }
      }
    }

    entropy2 += seq_likelihood;

#   ifdef MESSAGE
    if ((entropy2 < entropy1 - DOUBLE_ERROR) || (entropy2 > entropy1 + DOUBLE_ERROR)) {
      cout << "\nERROR: " << entropy1 << " " << entropy2 << endl;
    }
#   endif

    for (i = 0;i < seq.length[index];i++) {
      begin_partial_entropy[i] = 0.;
      for (j = 1;j < nb_row;j++) {
        if (backward[i][j] > 0.) {
          begin_partial_entropy[i] += backward[i][j] * (forward_state_entropy[i][j] - log(backward[i][j]));
        }
      }
      if (begin_partial_entropy[i] < 0.) {
        begin_partial_entropy[i] = 0.;
      }
    }

    begin_conditional_entropy[0] = 0.;
    for (i = 1;i < nb_row;i++) {
      if (backward[0][i] > 0.) {
        begin_conditional_entropy[0] -= backward[0][i] * log(backward[0][i]);
      }
    }
    if (begin_conditional_entropy[0] < 0.) {
      begin_conditional_entropy[0] = 0.;
    }

    for (i = 1;i < seq.length[index];i++) {
      begin_conditional_entropy[i] = 0.;
      for (j = 1;j < nb_row;j++) {
        for (k = 0;k < nb_memory[j];k++) {
          if ((predicted[i][j] > 0.) && (backward[i - 1][previous[j][k]] > 0.)) {
            buff = backward[i][j] * transition[previous[j][k]][state[j][0]] *
                   forward[i - 1][previous[j][k]] / predicted[i][j];
            if (buff > 0.) {
              begin_conditional_entropy[i] -= buff * log(buff / backward[i - 1][previous[j][k]]);
            }
          }
        }
      }
      if (begin_conditional_entropy[i] < 0.) {
        begin_conditional_entropy[i] = 0.;
      }
    }

/*    for (i = 0;i < seq.length[index];i++) {
      end_partial_entropy[i] = begin_partial_entropy[seq.length[index] - 1];
      for (j = 1;j < nb_row;j++) {
        if (backward[i][j] > 0.) {
          end_partial_entropy[i - order[j] + 1] -= backward[i][j] * forward_state_entropy[i][j];
        }
      }
      if (end_partial_entropy[i] < 0.) {
        end_partial_entropy[i] = 0.;
      }
    }

    for (i = 0;i < seq.length[index];i++) {
      end_partial_entropy[i] = 0.;
    }

    for (i = 0;i < seq.length[index] - 1;i++) {
      for (j = 1;j < nb_row;j++) {
        if (backward[i][j] > 0.) {
          end_partial_entropy[i - order[j] + 1] += backward[i][j] * (backward_state_entropy[i][j] - log(backward[i][j]));
        }
      }
    }

    i = seq.length[index] - 1;
    for (j = 1;j < nb_row;j++) {
      if (end_backward[j] > 0.) {
        end_partial_entropy[i - order[j] + 1] -= end_backward[j] * log(end_backward[j]);
      }
    }

    for (i = 0;i < seq.length[index];i++) {
      if (end_partial_entropy[i] < 0.) {
        end_partial_entropy[i] = 0.;
      }
    } */

    for (i = 0;i < seq.length[index];i++) {
      end_conditional_entropy[i] = 0.;
    }

    for (i = 0;i < seq.length[index] - 1;i++) {
      for (j = 1;j < nb_row;j++) {
        if ((next[j]) && (i - order[j] + 1 >= 0)) {
          for (k = 0;k < nb_state;k++) {
            if (predicted[i + 1][next[j][k]] > 0.) {
              buff = transition[j][k] * forward[i][j] / predicted[i + 1][next[j][k]];
              if (buff > 0.) {
                end_conditional_entropy[i - order[j] + 1] -= (backward[i + 1][next[j][k]] * buff) * log(buff);
              }
            }
          }
        }
      }
    }

    i = seq.length[index] - 1;
    end_backward[0] = 0.;
    for (j = 1;j < nb_row;j++) {
      end_backward[j] = backward[i][j];
    }
    for (j = nb_row - 1;j >= 1;j--) {
      end_backward[parent[j]] += end_backward[j];
    }

#   ifdef DEBUG
    cout << "\nTEST sum to 1: " << end_backward[0] << endl;
#   endif

    for (j = 1;j < nb_row;j++) {
      if ((i - order[j] + 1 >= 0) && (end_backward[j] > 0.) && (end_backward[parent[j]] > 0.)) {
        end_conditional_entropy[i - order[j] + 1] -= end_backward[j] * log(end_backward[j] / end_backward[parent[j]]);
      }
    }

    for (i = 0;i < seq.length[index];i++) {
      if (end_conditional_entropy[i] < 0.) {
        end_conditional_entropy[i] = 0.;
      }
    }

#   ifdef MESSAGE
    buff = begin_conditional_entropy[0];
    if ((buff < begin_partial_entropy[0] - DOUBLE_ERROR) || (buff > begin_partial_entropy[0] + DOUBLE_ERROR)) {
      cout << "\nERROR: " << 0 << " | " << buff << " " << begin_partial_entropy[0] << endl;
    }
    for (i = 1;i < seq.length[index];i++) {
      buff += begin_conditional_entropy[i];
      if ((buff < begin_partial_entropy[i] - DOUBLE_ERROR) || (buff > begin_partial_entropy[i] + DOUBLE_ERROR)) {
        cout << "\nERROR: " << i << " | " << buff << " " << begin_partial_entropy[i] << endl;
      }
    }

/*    i = seq.length[index] - 1;
    buff = end_conditional_entropy[i];
    if ((buff < end_partial_entropy[i] - DOUBLE_ERROR) || (buff > end_partial_entropy[i] + DOUBLE_ERROR)) {
      cout << "\nERROR: " << i << " | " << buff << " " << end_partial_entropy[i] << endl;
    }
    for (i = seq.length[index] - 2;i >= 0;i--) {
      buff += end_conditional_entropy[i];
      if ((buff < end_partial_entropy[i] - DOUBLE_ERROR) || (buff > end_partial_entropy[i] + DOUBLE_ERROR)) {
        cout << "\nERROR: " << i << " | " << buff << " " << end_partial_entropy[i] << endl;
      }
    } */
#   endif

    // restauration

    pstate = seq.int_sequence[index][0];

    for (i = 0;i < seq.length[index];i++) {
      for (j = 0;j < nb_state;j++) {
        state_backward[i][j] = 0.;
      }
      for (j = 1;j < nb_row;j++) {
        state_backward[i][state[j][0]] += backward[i][j];
      }
    }

    for (i = 0;i < seq.length[index];i++) {
      backward_max = 0.;
      for (j = 0;j < nb_state;j++) {
        if (state_backward[i][j] > backward_max) {
          backward_max = state_backward[i][j];
          *pstate = j;
        }
      }

      pstate++;
    }

    state_seq_likelihood = VariableOrderMarkov::likelihood_computation(seq , index);

/*    begin_conditional_entropy[0] = begin_partial_entropy[0];
    for (i = 1;i < seq.length[index];i++) {
      begin_conditional_entropy[i] = begin_partial_entropy[i] - begin_partial_entropy[i - 1];
    } */

    begin_partial_entropy[0] = begin_conditional_entropy[0];
    for (i = 1;i < seq.length[index];i++) {
      begin_partial_entropy[i] = begin_partial_entropy[i - 1] + begin_conditional_entropy[i];
    }

    end_partial_entropy[seq.length[index] - 1] = end_conditional_entropy[seq.length[index] - 1];
    for (i = seq.length[index] - 2;i >= 0;i--) {
      end_partial_entropy[i] = end_partial_entropy[i + 1] + end_conditional_entropy[i];
    }

    max_marginal_entropy = 0.;
    for (i = 0;i < seq.length[index];i++) {
      marginal_entropy[i] = 0.;
/*      for (j = 0;j < nb_state;j++) {
        if (state_backward[i][j] > 0.) {
          marginal_entropy[i] -= state_backward[i][j] * log(state_backward[i][j]);
        }
      } */
      for (j = 1;j < nb_row;j++) {
        if (backward[i][j] > 0.) {
          marginal_entropy[i] -= backward[i][j] * log(backward[i][j]);
        }
      }
      if (marginal_entropy[i] > max_marginal_entropy) {
        max_marginal_entropy = marginal_entropy[i];
      }
      if (marginal_entropy[i] < 0.) {
        marginal_entropy[i] = 0.;
      }
    }

    switch (format) {

    case 'a' : {
      *os << "\n" << SEQ_label[SEQL_POSTERIOR_STATE_PROBABILITY] << "\n\n";
//      seq.profile_ascii_print(*os , index , nb_state , state_backward ,
//                              STAT_label[STATL_STATE]);
      seq.profile_ascii_print(*os , index , nb_state , state_backward , begin_conditional_entropy ,
                              marginal_entropy , begin_partial_entropy , end_conditional_entropy ,
                              end_partial_entropy);

      *os << "\n" << STAT_label[STATL_LIKELIHOOD] << ": " << seq_likelihood
          << "\n" << SEQ_label[SEQL_STATE_SEQUENCE_LIKELIHOOD] << ": " << state_seq_likelihood
          << "   (" << exp(state_seq_likelihood - seq_likelihood) << ")" << endl;
      break;
    }

    case 's' : {
      *os << "\n" << SEQ_label[SEQL_POSTERIOR_STATE_PROBABILITY] << "\n\n";
//      seq.profile_spreadsheet_print(*os , index , nb_state , state_backward ,
//                                    STAT_label[STATL_STATE]);
      seq.profile_spreadsheet_print(*os , index , nb_state , state_backward , begin_conditional_entropy ,
                                    marginal_entropy , begin_partial_entropy , end_conditional_entropy ,
                                    end_partial_entropy);

      *os << "\n" << STAT_label[STATL_LIKELIHOOD] << "\t" << seq_likelihood
          << "\n" << SEQ_label[SEQL_STATE_SEQUENCE_LIKELIHOOD] << "\t" << state_seq_likelihood
          << "\t" << exp(state_seq_likelihood - seq_likelihood) << endl;
      break;
    }

    case 'g' : {
//      seq.profile_plot_print(*os , index , nb_state , state_backward);
      seq.profile_plot_print(*os , index , nb_state , state_backward , begin_conditional_entropy ,
                             marginal_entropy , begin_partial_entropy , end_conditional_entropy ,
                             end_partial_entropy);
      break;
    }

    case 'p' : {
      seq.profile_plotable_write((*plot_set)[1] , index , nb_state , state_backward);
      seq.entropy_profile_plotable_write((*plot_set)[2] , index , begin_conditional_entropy ,
                                         end_conditional_entropy , marginal_entropy);
      seq.entropy_profile_plotable_write((*plot_set)[3] , index , begin_partial_entropy ,
                                         end_partial_entropy);
      break;
    }
    }

    if (format != 'g') {
/*      double gini_index;

      gini_index = 0.;
      for (i = 0;i < seq.length[index];i++) {
        for (j = 0;j < nb_state;j++) {
          gini_index += state_backward[i][j] * (1. - state_backward[i][j]);
        }
      } */

      double entropy3 , observation , nb_state_sequence;

      entropy3 = 0.;
      for (i = 0;i < seq.length[index];i++) {
        for (j = 0;j < nb_state;j++) {
          if (state_backward[i][j] > 0.) {
            entropy3 -= state_backward[i][j] * log(state_backward[i][j]);
          }
        }
      }

      // calcul du nombre de sequences d'etats possibles

      for (i = 0;i < nb_output_process;i++) {
        switch (seq.type[i + 1]) {
        case INT_VALUE :
          pioutput[i] = seq.int_sequence[index][i + 1];
          break;
        case REAL_VALUE :
          proutput[i] = seq.real_sequence[index][i + 1];
          break;
        }
      }

      // recurrence "forward"

      switch (type) {

      case 'o' : {
        for (i = 1;i < nb_row;i++) {
          if (order[i] == 1) {
            forward[0][i] = initial[state[i][0]];

            for (j = 0;j < nb_output_process;j++) {
              if (categorical_process[j]) {
                forward[0][i] *= categorical_process[j]->observation[state[i][0]]->mass[*pioutput[j]];
              }

              else if (discrete_parametric_process[j]) {
                forward[0][i] *= discrete_parametric_process[j]->observation[state[i][0]]->mass[*pioutput[j]];
              }

              else {
                if (((continuous_parametric_process[j]->ident == GAMMA) ||
                    (continuous_parametric_process[j]->ident == ZERO_INFLATED_GAMMA)) && (seq.min_value[j + 1] < seq.min_interval[j + 1] / 2)) {
                  switch (seq.type[j + 1]) {
                  case INT_VALUE :
                    forward[0][i] *= continuous_parametric_process[j]->observation[state[i][0]]->mass_computation(*pioutput[j] , *pioutput[j] + seq.min_interval[j + 1]);
                    break;
                  case REAL_VALUE :
                    forward[0][i] *= continuous_parametric_process[j]->observation[state[i][0]]->mass_computation(*proutput[j] , *proutput[j] + seq.min_interval[j + 1]);
                    break;
                  }
                }

                else {
                  switch (seq.type[j + 1]) {
                  case INT_VALUE :
                    forward[0][i] *= continuous_parametric_process[j]->observation[state[i][0]]->mass_computation(*pioutput[j] - seq.min_interval[j + 1] / 2 , *pioutput[j] + seq.min_interval[j + 1] / 2);
                    break;
                  case REAL_VALUE :
                    forward[0][i] *= continuous_parametric_process[j]->observation[state[i][0]]->mass_computation(*proutput[j] - seq.min_interval[j + 1] / 2 , *proutput[j] + seq.min_interval[j + 1] / 2);
                    break;
                  }
                }
              }
            }

            if (forward[0][i] > 0.) {
              forward[0][i] = 1.;
            }
          }

          else {
            forward[0][i] = 0.;
          }
        }
        break;
      }

      case 'e' : {
        for (i = 1;i < nb_row;i++) {
          if (!child[i]) {
            forward[0][i] = initial[i];

            for (j = 0;j < nb_output_process;j++) {
              if (categorical_process[j]) {
                forward[0][i] *= categorical_process[j]->observation[state[i][0]]->mass[*pioutput[j]];
              }

              else if (discrete_parametric_process[j]) {
                forward[0][i] *= discrete_parametric_process[j]->observation[state[i][0]]->mass[*pioutput[j]];
              }

              else {
                if (((continuous_parametric_process[j]->ident == GAMMA) ||
                    (continuous_parametric_process[j]->ident == ZERO_INFLATED_GAMMA)) && (seq.min_value[j + 1] < seq.min_interval[j + 1] / 2)) {
                  switch (seq.type[j + 1]) {
                  case INT_VALUE :
                    forward[0][i] *= continuous_parametric_process[j]->observation[state[i][0]]->mass_computation(*pioutput[j] , *pioutput[j] + seq.min_interval[j + 1]);
                    break;
                  case REAL_VALUE :
                    forward[0][i] *= continuous_parametric_process[j]->observation[state[i][0]]->mass_computation(*proutput[j] , *proutput[j] + seq.min_interval[j + 1]);
                    break;
                  }
                }

                else {
                  switch (seq.type[j + 1]) {
                  case INT_VALUE :
                    forward[0][i] *= continuous_parametric_process[j]->observation[state[i][0]]->mass_computation(*pioutput[j] - seq.min_interval[j + 1] / 2 , *pioutput[j] + seq.min_interval[j + 1] / 2);
                    break;
                  case REAL_VALUE :
                    forward[0][i] *= continuous_parametric_process[j]->observation[state[i][0]]->mass_computation(*proutput[j] - seq.min_interval[j + 1] / 2 , *proutput[j] + seq.min_interval[j + 1] / 2);
                    break;
                  }
                }
              }
            }

            if (forward[0][i] > 0.) {
              forward[0][i] = 1.;
            }
          }

          else {
            forward[0][i] = 0.;
          }
        }
        break;
      }
      }

      for (i = 1;i < seq.length[index];i++) {
        for (j = 0;j < nb_output_process;j++) {
          switch (seq.type[j + 1]) {
          case INT_VALUE :
            pioutput[j]++;
            break;
          case REAL_VALUE :
            proutput[j]++;
            break;
          }
        }

        for (j = 1;j < nb_row;j++) {
          forward[i][j] = 0.;
          for (k = 0;k < nb_memory[j];k++) {
            if (transition[previous[j][k]][state[j][0]] > 0.) {
              forward[i][j] += forward[i - 1][previous[j][k]];
            }
          }

          observation = 1.;
          for (k = 0;k < nb_output_process;k++) {
            if (categorical_process[k]) {
              observation *= categorical_process[k]->observation[state[j][0]]->mass[*pioutput[k]];
            }

            else if (discrete_parametric_process[k]) {
              observation *= discrete_parametric_process[k]->observation[state[j][0]]->mass[*pioutput[k]];
            }

            else {
              if (((continuous_parametric_process[k]->ident == GAMMA) ||
                  (continuous_parametric_process[k]->ident == ZERO_INFLATED_GAMMA)) && (seq.min_value[k + 1] < seq.min_interval[k + 1] / 2)) {
                switch (seq.type[k + 1]) {
                case INT_VALUE :
                  observation *= continuous_parametric_process[k]->observation[state[j][0]]->mass_computation(*pioutput[k] , *pioutput[k] + seq.min_interval[k + 1]);
                  break;
                case REAL_VALUE :
                  observation *= continuous_parametric_process[k]->observation[state[j][0]]->mass_computation(*proutput[k] , *proutput[k] + seq.min_interval[k + 1]);
                  break;
                }
              }

              else {
                switch (seq.type[k + 1]) {
                case INT_VALUE :
                  observation *= continuous_parametric_process[k]->observation[state[j][0]]->mass_computation(*pioutput[k] - seq.min_interval[k + 1] / 2 , *pioutput[k] + seq.min_interval[k + 1] / 2);
                  break;
                case REAL_VALUE :
                  observation *= continuous_parametric_process[k]->observation[state[j][0]]->mass_computation(*proutput[k] - seq.min_interval[k + 1] / 2 , *proutput[k] + seq.min_interval[k + 1] / 2);
                  break;
                }
              }
            }
          }

          if (observation == 0.) {
            forward[i][j] = 0.;
          }
        }
      }

      nb_state_sequence = 0.;
      i = seq.length[index] - 1;
      for (j = 1;j < nb_row;j++) {
        nb_state_sequence += forward[i][j];
      }

      switch (format) {
      case 'a' :
/*        *os << "\n" << SEQ_label[SEQL_GINI_INDEX] << ": " << gini_index << " ("
            << gini_index / seq.length[index] */
        *os << "\n" << SEQ_label[SEQL_STATE_SEQUENCE_ENTROPY] << ": " << entropy1
            << " (" << entropy1 / seq.length[index] << ")   " << SEQ_label[SEQL_UPPER_BOUND] << ": "
            << log((double)nb_state_sequence) << " ("
            << log((double)nb_state_sequence) / seq.length[index]
            << ")\n" << SEQ_label[SEQL_MARGINAL_ENTROPY_SUM] << ": " << entropy3 << " ("
            << entropy3 / seq.length[index] << ")\n\n"
            << SEQ_label[SEQL_NB_STATE_SEQUENCE] << ": " << nb_state_sequence << endl;
        break;
      case 's' :
/*        *os << "\n" << SEQ_label[SEQL_GINI_INDEX] << "\t" << gini_index << "\t"
            << gini_index / seq.length[index] */
        *os << "\n" << SEQ_label[SEQL_STATE_SEQUENCE_ENTROPY] << "\t" << entropy1
            << "\t" << entropy1 / seq.length[index] << "\t" << SEQ_label[SEQL_UPPER_BOUND] << "\t"
            << log((double)nb_state_sequence) << "\t"
            << log((double)nb_state_sequence) / seq.length[index]
            << "\n" << SEQ_label[SEQL_MARGINAL_ENTROPY_SUM] << "\t" << entropy3 << "\t"
            << entropy3 / seq.length[index] << "\n\n"
            << SEQ_label[SEQL_NB_STATE_SEQUENCE] << "\t" << nb_state_sequence << endl;
        break;
      }
    }
  }

  for (i = 0;i < seq.length[index];i++) {
    delete [] forward[i];
  }
  delete [] forward;

  for (i = 0;i < seq.length[index];i++) {
    delete [] predicted[i];
  }
  delete [] predicted;

  for (i = 0;i < seq.length[index];i++) {
    delete [] backward[i];
  }
  delete [] backward;

  delete [] auxiliary;

  for (i = 0;i < seq.length[index];i++) {
    delete [] state_backward[i];
  }
  delete [] state_backward;

  delete [] transition_predicted;

  for (i = 0;i < seq.length[index];i++) {
    delete [] forward_state_entropy[i];
  }
  delete [] forward_state_entropy;

  for (i = 1;i < nb_row;i++) {
    delete [] transition_entropy[i];
  }
  delete [] transition_entropy;

  delete [] begin_partial_entropy;
  delete [] begin_conditional_entropy;

/*  for (i = 0;i < seq.length[index];i++) {
    delete [] backward_state_entropy[i];
  }
  delete [] backward_state_entropy; */

  delete [] end_backward;
  delete [] end_partial_entropy;
  delete [] end_conditional_entropy;

  delete [] marginal_entropy;

  delete [] pioutput;
  delete [] proutput;

  return seq_likelihood;
}


/*--------------------------------------------------------------*
 *
 *  Simulation de sequences d'etats correspondant a une sequence observee
 *  par l'algorithme forward-backward de simulation.
 *
 *  arguments : reference sur un objet MarkovianSequences, indice de la sequence,
 *              stream, format de fichier ('a' : ASCII, 's' : Spreadsheet),
 *              nombre de sequences d'etats.
 *
 *--------------------------------------------------------------*/

double HiddenVariableOrderMarkov::forward_backward_sampling(const MarkovianSequences &seq ,
                                                            int index , ostream &os , char format ,
                                                            int nb_state_sequence) const

{
  register int i , j , k;
  int memory , *pstate , **pioutput;
  double seq_likelihood , state_seq_likelihood , **forward , norm , **predicted ,
         *backward , *cumul_backward , **proutput;

# ifdef DEBUG
  double sum;
# endif


  // initialisations

  forward = new double*[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    forward[i] = new double[nb_row];
  }

  predicted = new double*[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    predicted[i] = new double[nb_row];
  }

  backward = new double[nb_row];
  cumul_backward = new double[nb_row];

  pioutput = new int*[nb_output_process];
  proutput = new double*[nb_output_process];

  for (i = 0;i < nb_output_process;i++) {
    switch (seq.type[i + 1]) {
    case INT_VALUE :
      pioutput[i] = seq.int_sequence[index][i + 1];
      break;
    case REAL_VALUE :
      proutput[i] = seq.real_sequence[index][i + 1];
      break;
    }
  }

# ifdef DEBUG
  double **state_sequence_probability;


  state_sequence_probability = new double*[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    state_sequence_probability[i] = new double[nb_state];
    for (j = 0;j < nb_state;j++) {
      state_sequence_probability[i][j] = 0.;
    }
  }
# endif

  // recurrence "forward"

  seq_likelihood = 0.;
  norm = 0.;

  switch (type) {

  case 'o' : {
    for (i = 1;i < nb_row;i++) {
      if (order[i] == 1) {
        forward[0][i] = initial[state[i][0]];

        for (j = 0;j < nb_output_process;j++) {
          if (categorical_process[j]) {
            forward[0][i] *= categorical_process[j]->observation[state[i][0]]->mass[*pioutput[j]];
          }

          else if (discrete_parametric_process[j]) {
            forward[0][i] *= discrete_parametric_process[j]->observation[state[i][0]]->mass[*pioutput[j]];
          }

          else {
            if (((continuous_parametric_process[j]->ident == GAMMA) ||
                (continuous_parametric_process[j]->ident == ZERO_INFLATED_GAMMA)) && (seq.min_value[j + 1] < seq.min_interval[j + 1] / 2)) {
              switch (seq.type[j + 1]) {
              case INT_VALUE :
                forward[0][i] *= continuous_parametric_process[j]->observation[state[i][0]]->mass_computation(*pioutput[j] , *pioutput[j] + seq.min_interval[j + 1]);
                break;
              case REAL_VALUE :
                forward[0][i] *= continuous_parametric_process[j]->observation[state[i][0]]->mass_computation(*proutput[j] , *proutput[j] + seq.min_interval[j + 1]);
                break;
              }
            }

            else {
              switch (seq.type[j + 1]) {
              case INT_VALUE :
                forward[0][i] *= continuous_parametric_process[j]->observation[state[i][0]]->mass_computation(*pioutput[j] - seq.min_interval[j + 1] / 2 , *pioutput[j] + seq.min_interval[j + 1] / 2);
                break;
              case REAL_VALUE :
                forward[0][i] *= continuous_parametric_process[j]->observation[state[i][0]]->mass_computation(*proutput[j] - seq.min_interval[j + 1] / 2 , *proutput[j] + seq.min_interval[j + 1] / 2);
                break;
              }
            }
          }
        }

        norm += forward[0][i];
      }

      else {
        forward[0][i] = 0.;
      }
    }
    break;
  }

  case 'e' : {
    for (i = 1;i < nb_row;i++) {
      if (!child[i]) {
        forward[0][i] = initial[i];

        for (j = 0;j < nb_output_process;j++) {
          if (categorical_process[j]) {
            forward[0][i] *= categorical_process[j]->observation[state[i][0]]->mass[*pioutput[j]];
          }

          else if (discrete_parametric_process[j]) {
            forward[0][i] *= discrete_parametric_process[j]->observation[state[i][0]]->mass[*pioutput[j]];
          }

          else {
            if (((continuous_parametric_process[j]->ident == GAMMA) ||
                (continuous_parametric_process[j]->ident == ZERO_INFLATED_GAMMA)) && (seq.min_value[j + 1] < seq.min_interval[j + 1] / 2)) {
              switch (seq.type[j + 1]) {
              case INT_VALUE :
                forward[0][i] *= continuous_parametric_process[j]->observation[state[i][0]]->mass_computation(*pioutput[j] , *pioutput[j] + seq.min_interval[j + 1]);
                break;
              case REAL_VALUE :
                forward[0][i] *= continuous_parametric_process[j]->observation[state[i][0]]->mass_computation(*proutput[j] , *proutput[j] + seq.min_interval[j + 1]);
                break;
              }
            }

            else {
              switch (seq.type[j + 1]) {
              case INT_VALUE :
                forward[0][i] *= continuous_parametric_process[j]->observation[state[i][0]]->mass_computation(*pioutput[j] - seq.min_interval[j + 1] / 2 , *pioutput[j] + seq.min_interval[j + 1] / 2);
                break;
              case REAL_VALUE :
                forward[0][i] *= continuous_parametric_process[j]->observation[state[i][0]]->mass_computation(*proutput[j] - seq.min_interval[j + 1] / 2 , *proutput[j] + seq.min_interval[j + 1] / 2);
                break;
              }
            }
          }
        }

        norm += forward[0][i];
      }

      else {
        forward[0][i] = 0.;
      }
    }
    break;
  }
  }

  if (norm > 0.) {
    for (i = 1;i < nb_row;i++) {
      forward[0][i] /= norm;
    }

    seq_likelihood += log(norm);
  }

  else {
    seq_likelihood = D_INF;
  }

  if (seq_likelihood != D_INF) {
    for (i = 1;i < seq.length[index];i++) {
      for (j = 0;j < nb_output_process;j++) {
        switch (seq.type[j + 1]) {
        case INT_VALUE :
          pioutput[j]++;
          break;
        case REAL_VALUE :
          proutput[j]++;
          break;
        }
      }
      norm = 0.;

      for (j = 1;j < nb_row;j++) {
        forward[i][j] = 0.;
        for (k = 0;k < nb_memory[j];k++) {
          forward[i][j] += transition[previous[j][k]][state[j][0]] * forward[i - 1][previous[j][k]];
        }
        predicted[i][j] = forward[i][j];

        for (k = 0;k < nb_output_process;k++) {
          if (categorical_process[k]) {
            forward[i][j] *= categorical_process[k]->observation[state[j][0]]->mass[*pioutput[k]];
          }

          else if (discrete_parametric_process[k]) {
            forward[i][j] *= discrete_parametric_process[k]->observation[state[j][0]]->mass[*pioutput[k]];
          }

          else {
            if (((continuous_parametric_process[k]->ident == GAMMA) ||
                (continuous_parametric_process[k]->ident == ZERO_INFLATED_GAMMA)) && (seq.min_value[k + 1] < seq.min_interval[k + 1] / 2)) {
              switch (seq.type[k + 1]) {
              case INT_VALUE :
                forward[i][j] *= continuous_parametric_process[k]->observation[state[j][0]]->mass_computation(*pioutput[k] , *pioutput[k] + seq.min_interval[k + 1]);
                break;
              case REAL_VALUE :
                forward[i][j] *= continuous_parametric_process[k]->observation[state[j][0]]->mass_computation(*proutput[k] , *proutput[k] + seq.min_interval[k + 1]);
                break;
              }
            }

            else {
              switch (seq.type[k + 1]) {
              case INT_VALUE :
                forward[i][j] *= continuous_parametric_process[k]->observation[state[j][0]]->mass_computation(*pioutput[k] - seq.min_interval[k + 1] / 2 , *pioutput[k] + seq.min_interval[k + 1] / 2);
                break;
              case REAL_VALUE :
                forward[i][j] *= continuous_parametric_process[k]->observation[state[j][0]]->mass_computation(*proutput[k] - seq.min_interval[k + 1] / 2 , *proutput[k] + seq.min_interval[k + 1] / 2);
                break;
              }
            }
          }
        }

        norm += forward[i][j];
      }

      if (norm > 0.) {
        for (j = 1;j < nb_row;j++) {
          forward[i][j] /= norm;
        }
        seq_likelihood += log(norm);
      }

      else {
        seq_likelihood = D_INF;
        break;
      }
    }
  }

  if (seq_likelihood != D_INF) {

#   ifdef MESSAGE
    cout << "\n";
#   endif

    // passes "backward"

    for (i = 0;i < nb_state_sequence;i++) {
      j = seq.length[index] - 1;
      pstate = seq.int_sequence[index][0] + j;
      stat_tool::cumul_computation(nb_row - 1 , forward[j] + 1 , cumul_backward);
      memory = 1 + cumul_method(nb_row - 1 , cumul_backward);
      *pstate = state[memory][0];

      for (j = seq.length[index] - 2;j >= 0;j--) {
        for (k = 0;k < nb_memory[memory];k++) {
          backward[k] = transition[previous[memory][k]][state[memory][0]] *
                        forward[j][previous[memory][k]] / predicted[j + 1][memory];
        }

#       ifdef DEBUG
        sum = 0.;
        for (k = 0;k < nb_memory[memory];k++) {
          sum += backward[k];
        }
        if ((sum < 1. - DOUBLE_ERROR) || (sum > 1. + DOUBLE_ERROR)) {
          cout << "\nERROR: " << j << " " << sum << endl;
        }
#       endif

        stat_tool::cumul_computation(nb_memory[memory] , backward , cumul_backward);
        memory = previous[memory][cumul_method(nb_memory[memory] , cumul_backward)];
        *--pstate = state[memory][0];
      }

#     ifdef DEBUG
      pstate = seq.int_sequence[index][0];
      for (j = 0;j < seq.length[index];j++) {
        state_sequence_probability[j][*pstate++]++;
      }
#     endif

#     ifdef MESSAGE
      state_seq_likelihood = VariableOrderMarkov::likelihood_computation(seq , index);

      pstate = seq.int_sequence[index][0];

      switch (format) {

      case 'a' : {
        for (j = 0;j < seq.length[index];j++) {
          os << *pstate++ << " ";
        }

        os << "  " << i + 1 << "  " << state_seq_likelihood
           << "   (" << exp(state_seq_likelihood - seq_likelihood) << ")" << endl;
        break;
      }

      case 's' : {
        for (j = 0;j < seq.length[index];j++) {
          os << *pstate++ << "\t";
        }

        os << "\t" << i + 1 << "\t" << state_seq_likelihood
           << "\t" << exp(state_seq_likelihood - seq_likelihood) << endl;
        break;
      }
      }
#     endif

    }

#   ifdef DEBUG
    if (nb_state_sequence >= 1000) {
      for (i = 0;i < seq.length[index];i++) {
        for (j = 0;j < nb_state;j++) {
          state_sequence_probability[i][j] /= nb_state_sequence;
        }
      }

      pstate = seq.int_sequence[index][0];
      for (j = 0;j < seq.length[index];j++) {
        *pstate++ = I_DEFAULT;
      }

      os << "\n" << SEQ_label[SEQL_POSTERIOR_STATE_PROBABILITY] << "\n\n";
      seq.profile_ascii_print(os , index , nb_state , state_sequence_probability ,
                              STAT_label[STATL_STATE]);
    }
#   endif

  }

  for (i = 0;i < seq.length[index];i++) {
    delete [] forward[i];
  }
  delete [] forward;

  for (i = 0;i < seq.length[index];i++) {
    delete [] predicted[i];
  }
  delete [] predicted;

  delete [] backward;
  delete [] cumul_backward;

  delete [] pioutput;
  delete [] proutput;

# ifdef DEBUG
  for (i = 0;i < seq.length[index];i++) {
    delete [] state_sequence_probability[i];
  }
  delete [] state_sequence_probability;
# endif

  return seq_likelihood;
}


/*--------------------------------------------------------------*
 *
 *  Calcul des logarithmes des parametres d'une chaine de Markov
 *  d'ordre variable cachee.
 *
 *--------------------------------------------------------------*/

void HiddenVariableOrderMarkov::log_computation()

{
  register int i , j;


  Chain::log_computation();

  for (i = 0;i < nb_output_process;i++) {
    if (categorical_process[i]) {
      for (j = 0;j < nb_state;j++) {
        categorical_process[i]->observation[j]->log_computation();
      }
    }

    else if (discrete_parametric_process[i]) {
      for (j = 0;j < nb_state;j++) {
        stat_tool::log_computation(discrete_parametric_process[i]->nb_value ,
                                   discrete_parametric_process[i]->observation[j]->mass ,
                                   discrete_parametric_process[i]->observation[j]->cumul);
      }
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Calcul des sequences d'etats les plus probables par l'algorithme de Viterbi.
 *
 *  arguments : reference sur un objet MarkovianSequences,
 *              pointeur sur les probabilites a posteriori des sequences d'etats
 *              les plus probables, indice de la sequence.
 *
 *--------------------------------------------------------------*/

double HiddenVariableOrderMarkov::viterbi(const MarkovianSequences &seq ,
                                          double *posterior_probability , int index) const

{
  register int i , j , k , m;
  int length , memory , *pstate , **pioutput , **optimal_memory;
  double likelihood = 0. , buff , forward_max , *forward , *previous_forward , **proutput;


  // initialisations

  forward = new double[nb_row];
  previous_forward = new double[nb_row];

  length = (index == I_DEFAULT ? seq.max_length : seq.length[index]);

  optimal_memory = new int*[length];
  for (i = 0;i < length;i++) {
    optimal_memory[i] = new int[nb_row];
  }

  pioutput = new int*[nb_output_process];
  proutput = new double*[nb_output_process];

  for (i = 0;i < seq.nb_sequence;i++) {
    if ((index == I_DEFAULT) || (index == i)) {
      for (j = 0;j < nb_output_process;j++) {
        switch (seq.type[j + 1]) {
        case INT_VALUE :
          pioutput[j] = seq.int_sequence[i][j + 1];
          break;
        case REAL_VALUE :
          proutput[j] = seq.real_sequence[i][j + 1];
          break;
        }
      }

      // recurrence "forward"

      switch (type) {

      case 'o' : {
        for (j = 1;j < nb_row;j++) {
          if (order[j] == 1) {
            forward[j] = cumul_initial[state[j][0]];

            if (forward[j] != D_INF) {
              for (k = 0;k < nb_output_process;k++) {
                if (categorical_process[k]) {
                  buff = categorical_process[k]->observation[state[j][0]]->cumul[*pioutput[k]];
                }

                else if (discrete_parametric_process[k]) {
                  buff = discrete_parametric_process[k]->observation[state[j][0]]->cumul[*pioutput[k]];
                }

                else {
                  if (((continuous_parametric_process[k]->ident == GAMMA) ||
                      (continuous_parametric_process[k]->ident == ZERO_INFLATED_GAMMA)) && (seq.min_value[k + 1] < seq.min_interval[k + 1] / 2)) {
                    switch (seq.type[k + 1]) {
                    case INT_VALUE :
                      buff = continuous_parametric_process[k]->observation[state[j][0]]->mass_computation(*pioutput[k] , *pioutput[k] + seq.min_interval[k + 1]);
                      break;
                    case REAL_VALUE :
                      buff = continuous_parametric_process[k]->observation[state[j][0]]->mass_computation(*proutput[k] , *proutput[k] + seq.min_interval[k + 1]);
                      break;
                    }
                  }

                  else {
                    switch (seq.type[k + 1]) {
                    case INT_VALUE :
                      buff = continuous_parametric_process[k]->observation[state[j][0]]->mass_computation(*pioutput[k] - seq.min_interval[k + 1] / 2 , *pioutput[k] + seq.min_interval[k + 1] / 2);
                      break;
                    case REAL_VALUE :
                      buff = continuous_parametric_process[k]->observation[state[j][0]]->mass_computation(*proutput[k] - seq.min_interval[k + 1] / 2 , *proutput[k] + seq.min_interval[k + 1] / 2);
                      break;
                    }
                  }

                  if (buff > 0.) {
                    buff = log(buff);
                  }
                  else {
                    buff = D_INF;
                  }
                }

                if (buff == D_INF) {
                  forward[j] = D_INF;
                  break;
                }
                else {
                  forward[j] += buff;
                }
              }
            }
          }

          else {
            forward[j] = D_INF;
          }
        }
        break;
      }

      case 'e' : {
        for (j = 1;j < nb_row;j++) {
          if (!child[j]) {
            forward[j] = cumul_initial[j];

            if (forward[j] != D_INF) {
              for (k = 0;k < nb_output_process;k++) {
                if (categorical_process[k]) {
                  buff = categorical_process[k]->observation[state[j][0]]->cumul[*pioutput[k]];
                }

                else if (discrete_parametric_process[k]) {
                  buff = discrete_parametric_process[k]->observation[state[j][0]]->cumul[*pioutput[k]];
                }

                else {
                  if (((continuous_parametric_process[k]->ident == GAMMA) ||
                      (continuous_parametric_process[k]->ident == ZERO_INFLATED_GAMMA)) && (seq.min_value[k + 1] < seq.min_interval[k + 1] / 2)) {
                    switch (seq.type[k + 1]) {
                    case INT_VALUE :
                      buff = continuous_parametric_process[k]->observation[state[j][0]]->mass_computation(*pioutput[k] , *pioutput[k] + seq.min_interval[k + 1]);
                      break;
                    case REAL_VALUE :
                      buff = continuous_parametric_process[k]->observation[state[j][0]]->mass_computation(*proutput[k] , *proutput[k] + seq.min_interval[k + 1]);
                      break;
                    }
                  }

                  else {
                    switch (seq.type[k + 1]) {
                    case INT_VALUE :
                      buff = continuous_parametric_process[k]->observation[state[j][0]]->mass_computation(*pioutput[k] - seq.min_interval[k + 1] / 2 , *pioutput[k] + seq.min_interval[k + 1] / 2);
                      break;
                    case REAL_VALUE :
                      buff = continuous_parametric_process[k]->observation[state[j][0]]->mass_computation(*proutput[k] - seq.min_interval[k + 1] / 2 , *proutput[k] + seq.min_interval[k + 1] / 2);
                      break;
                    }
                  }

                  if (buff > 0.) {
                    buff = log(buff);
                  }
                  else {
                    buff = D_INF;
                  }
                }

                if (buff == D_INF) {
                  forward[j] = D_INF;
                  break;
                }
                else {
                  forward[j] += buff;
                }
              }
            }
          }

          else {
            forward[j] = D_INF;
          }
        }
        break;
      }
      }

#     ifdef DEBUG
      cout << "\n" << 0 << " : ";
      for (j = 1;j < nb_row;j++) {
        cout << forward[j] << " | ";
      }
      cout << endl;
#     endif

      for (j = 1;j < seq.length[i];j++) {
        for (k = 0;k < nb_output_process;k++) {
          switch (seq.type[k + 1]) {
          case INT_VALUE :
            pioutput[k]++;
            break;
          case REAL_VALUE :
            proutput[k]++;
            break;
          }
        }

        for (k = 1;k < nb_row;k++) {
          previous_forward[k] = forward[k];
        }

        for (k = 1;k < nb_row;k++) {
          forward[k] = D_INF;
          for (m = 0;m < nb_memory[k];m++) {
            buff = cumul_transition[previous[k][m]][state[k][0]] + previous_forward[previous[k][m]];
            if (buff > forward[k]) {
              forward[k] = buff;
              optimal_memory[j][k] = previous[k][m];
            }
          }

          if (forward[k] != D_INF) {
            for (m = 0;m < nb_output_process;m++) {
              if (categorical_process[m]) {
                buff = categorical_process[m]->observation[state[k][0]]->cumul[*pioutput[m]];
              }

              else if (discrete_parametric_process[m]) {
                buff = discrete_parametric_process[m]->observation[state[k][0]]->cumul[*pioutput[m]];
              }

              else {
                if (((continuous_parametric_process[m]->ident == GAMMA) ||
                    (continuous_parametric_process[m]->ident == ZERO_INFLATED_GAMMA)) && (seq.min_value[m + 1] < seq.min_interval[m + 1] / 2)) {
                  switch (seq.type[m + 1]) {
                  case INT_VALUE :
                    buff = continuous_parametric_process[m]->observation[state[k][0]]->mass_computation(*pioutput[m] , *pioutput[m] + seq.min_interval[m + 1]);
                    break;
                  case REAL_VALUE :
                    buff = continuous_parametric_process[m]->observation[state[k][0]]->mass_computation(*proutput[m] , *proutput[m] + seq.min_interval[m + 1]);
                    break;
                  }
                }

                else {
                  switch (seq.type[m + 1]) {
                  case INT_VALUE :
                    buff = continuous_parametric_process[m]->observation[state[k][0]]->mass_computation(*pioutput[m] - seq.min_interval[m + 1] / 2 , *pioutput[m] + seq.min_interval[m + 1] / 2);
                    break;
                  case REAL_VALUE :
                    buff = continuous_parametric_process[m]->observation[state[k][0]]->mass_computation(*proutput[m] - seq.min_interval[m + 1] / 2 , *proutput[m] + seq.min_interval[m + 1] / 2);
                    break;
                  }
                }

                if (buff > 0.) {
                  buff = log(buff);
                }
                else {
                  buff = D_INF;
                }
              }

              if (buff == D_INF) {
                forward[k] = D_INF;
                break;
              }
              else {
                forward[k] += buff;
              }
            }
          }
        }

#       ifdef DEBUG
        cout << j << " : ";
        for (k = 1;k < nb_row;k++) {
          cout << forward[k] << " " << optimal_memory[j][k] << " | ";
        }
        cout << endl;
#       endif

      }

      // extraction de la vraisemblance du chemin optimal

      pstate = seq.int_sequence[i][0] + seq.length[i] - 1;
      forward_max = D_INF;

      for (j = 1;j < nb_row;j++) {
        if (forward[j] > forward_max) {
          forward_max = forward[j];
          memory = j;
        }
      }

      if (forward_max != D_INF) {
        likelihood += forward_max;
        *pstate = state[memory][0];
        if (posterior_probability) {
          posterior_probability[i] = forward_max;
        }
      }

      else {
        likelihood = D_INF;
        if (posterior_probability) {
          posterior_probability[i] = 0.;
        }
        break;
      }

      // restauration

      for (j = seq.length[i] - 1;j > 0;j--) {
        memory = optimal_memory[j][memory];
        *--pstate = state[memory][0];
      }

#     ifdef DEBUG
      cout << "\n";
      for (j = seq.length[i] - 1;j >= 0;j--) {
        cout << seq.int_sequence[i][0][j] << " ";
      }
      cout << endl;
#     endif

    }
  }

  delete [] forward;
  delete [] previous_forward;

  for (i = 0;i < length;i++) {
    delete [] optimal_memory[i];
  }
  delete [] optimal_memory;

  delete [] pioutput;
  delete [] proutput;

  return likelihood;
}


/*--------------------------------------------------------------*
 *
 *  Calcul des sequences d'etats les plus probables par l'algorithme de Viterbi.
 *
 *  argument : reference sur un objet VariableOrderMarkovData.
 *
 *--------------------------------------------------------------*/

void HiddenVariableOrderMarkov::viterbi(VariableOrderMarkovData &seq) const

{
  seq.posterior_probability = new double[seq.nb_sequence];
  seq.restoration_likelihood = viterbi(seq , seq.posterior_probability);
}


/*--------------------------------------------------------------*
 *
 *  Calcul des N sequences d'etats les plus probables correspondant a
 *  une sequence observee par l'algorithme de Viterbi generalise.
 *
 *  arguments : reference sur un objet MarkovianSequences, indice de la sequence,
 *              stream, vraisemblance des donnees, format de fichier
 *              ('a' : ASCII, 's' : Spreadsheet), nombre de sequences d'etats.
 *
 *--------------------------------------------------------------*/

double HiddenVariableOrderMarkov::generalized_viterbi(const MarkovianSequences &seq , int index ,
                                                      ostream &os , double seq_likelihood ,
                                                      char format , int inb_state_sequence) const

{
  bool **active_cell;
  register int i , j , k , m;
  int nb_state_sequence , memory , brank , previous_rank , nb_cell , *rank , *pstate ,
      **pioutput , ***optimal_memory , ***optimal_rank;
  double buff , observation , forward_max , state_seq_likelihood , likelihood_cumul ,
         **forward , **previous_forward , **proutput;


  // initialisations

  forward = new double*[nb_row];
  forward[0] = NULL;
  for (i = 1;i < nb_row;i++) {
    forward[i] = new double[inb_state_sequence];
  }

  previous_forward = new double*[nb_row];
  previous_forward[0] = NULL;
  for (i = 1;i < nb_row;i++) {
    previous_forward[i] = new double[inb_state_sequence];
    for (j = 1;j < inb_state_sequence;j++) {
      previous_forward[i][j] = D_INF;
    }
  }

  rank = new int[nb_row];

  optimal_memory = new int**[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    optimal_memory[i] = new int*[nb_row];
    optimal_memory[i][0] = NULL;
    for (j = 1;j < nb_row;j++) {
      optimal_memory[i][j] = new int[inb_state_sequence];
    }
  }

  optimal_rank = new int**[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    optimal_rank[i] = new int*[nb_row];
    optimal_rank[i][0] = NULL;
    for (j = 1;j < nb_row;j++) {
      optimal_rank[i][j] = new int[inb_state_sequence];
    }
  }

  active_cell = new bool*[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    active_cell[i] = new bool[nb_state];
    for (j = 0;j < nb_state;j++) {
      active_cell[i][j] = false;
    }
  }

  pioutput = new int*[nb_output_process];
  proutput = new double*[nb_output_process];

  for (i = 0;i < nb_output_process;i++) {
    switch (seq.type[i + 1]) {
    case INT_VALUE :
      pioutput[i] = seq.int_sequence[index][i + 1];
      break;
    case REAL_VALUE :
      proutput[i] = seq.real_sequence[index][i + 1];
      break;
    }
  }

# ifdef DEBUG
  double entropy = 0. , **state_sequence_probability;


  state_sequence_probability = new double*[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    state_sequence_probability[i] = new double[nb_state];
    for (j = 0;j < nb_state;j++) {
//      state_sequence_probability[i][j] = 0.;
      state_sequence_probability[i][j] = D_INF;
    }
  }
# endif

  // recurrence "forward"

  switch (type) {

  case 'o' : {
    for (i = 1;i < nb_row;i++) {
      if (order[i] == 1) {
        forward[i][0] = cumul_initial[state[i][0]];

        if (forward[i][0] != D_INF) {
          for (j = 0;j < nb_output_process;j++) {
            if (categorical_process[j]) {
              buff = categorical_process[j]->observation[state[i][0]]->cumul[*pioutput[j]];
            }

            else if (discrete_parametric_process[j]) {
              buff = discrete_parametric_process[j]->observation[state[i][0]]->cumul[*pioutput[j]];
            }

            else {
              if (((continuous_parametric_process[j]->ident == GAMMA) ||
                  (continuous_parametric_process[j]->ident == ZERO_INFLATED_GAMMA)) && (seq.min_value[j + 1] < seq.min_interval[j + 1] / 2)) {
                switch (seq.type[j + 1]) {
                case INT_VALUE :
                  buff = continuous_parametric_process[j]->observation[state[i][0]]->mass_computation(*pioutput[j] , *pioutput[j] + seq.min_interval[j + 1]);
                  break;
                case REAL_VALUE :
                   buff = continuous_parametric_process[j]->observation[state[i][0]]->mass_computation(*proutput[j] , *proutput[j] + seq.min_interval[j + 1]);
                  break;
                }
              }

              else {
                switch (seq.type[j + 1]) {
                case INT_VALUE :
                  buff = continuous_parametric_process[j]->observation[state[i][0]]->mass_computation(*pioutput[j] - seq.min_interval[j + 1] / 2 , *pioutput[j] + seq.min_interval[j + 1] / 2);
                  break;
                case REAL_VALUE :
                   buff = continuous_parametric_process[j]->observation[state[i][0]]->mass_computation(*proutput[j] - seq.min_interval[j + 1] / 2 , *proutput[j] + seq.min_interval[j + 1] / 2);
                  break;
                }
              }

              if (buff > 0.) {
                buff = log(buff);
              }
              else {
                buff = D_INF;
              }
            }

            if (buff == D_INF) {
              forward[i][0] = D_INF;
              break;
            }
            else {
              forward[i][0] += buff;
            }
          }
        }
      }

      else {
        forward[i][0] = D_INF;
      }
    }
    break;
  }

  case 'e' : {
    for (i = 1;i < nb_row;i++) {
      if (!child[i]) {
        forward[i][0] = cumul_initial[i];

        if (forward[i][0] != D_INF) {
          for (j = 0;j < nb_output_process;j++) {
            if (categorical_process[j]) {
              buff = categorical_process[j]->observation[state[i][0]]->cumul[*pioutput[j]];
            }

            else if (discrete_parametric_process[j]) {
              buff = discrete_parametric_process[j]->observation[state[i][0]]->cumul[*pioutput[j]];
            }

            else {
              if (((continuous_parametric_process[j]->ident == GAMMA) ||
                  (continuous_parametric_process[j]->ident == ZERO_INFLATED_GAMMA)) && (seq.min_value[j + 1] < seq.min_interval[j + 1] / 2)) {
                switch (seq.type[j + 1]) {
                case INT_VALUE :
                  buff = continuous_parametric_process[j]->observation[state[i][0]]->mass_computation(*pioutput[j] , *pioutput[j] + seq.min_interval[j + 1]);
                  break;
                case REAL_VALUE :
                   buff = continuous_parametric_process[j]->observation[state[i][0]]->mass_computation(*proutput[j] , *proutput[j] + seq.min_interval[j + 1]);
                  break;
                }
              }

              else {
                switch (seq.type[j + 1]) {
                case INT_VALUE :
                  buff = continuous_parametric_process[j]->observation[state[i][0]]->mass_computation(*pioutput[j] - seq.min_interval[j + 1] / 2 , *pioutput[j] + seq.min_interval[j + 1] / 2);
                  break;
                case REAL_VALUE :
                   buff = continuous_parametric_process[j]->observation[state[i][0]]->mass_computation(*proutput[j] - seq.min_interval[j + 1] / 2 , *proutput[j] + seq.min_interval[j + 1] / 2);
                  break;
                }
              }

              if (buff > 0.) {
                buff = log(buff);
              }
              else {
                buff = D_INF;
              }
            }

            if (buff == D_INF) {
              forward[i][0] = D_INF;
              break;
            }
            else {
              forward[i][0] += buff;
            }
          }
        }
      }

      else {
        forward[i][0] = D_INF;
      }
    }
    break;
  }
  }

  nb_state_sequence = 1;

  for (i = 1;i < seq.length[index];i++) {
    for (j = 0;j < nb_output_process;j++) {
      switch (seq.type[j + 1]) {
      case INT_VALUE :
        pioutput[j]++;
        break;
      case REAL_VALUE :
        proutput[j]++;
        break;
      }
    }

    for (j = 1;j < nb_row;j++) {
      for (k = 0;k < nb_state_sequence;k++) {
        previous_forward[j][k] = forward[j][k];
      }
    }

    if (nb_state_sequence < inb_state_sequence) {
      if (nb_state_sequence * nb_state < inb_state_sequence) {
        nb_state_sequence *= nb_state;
      }
      else {
        nb_state_sequence = inb_state_sequence;
      }
    }

    for (j = 1;j < nb_row;j++) {
      observation = 0.;

      for (k = 0;k < nb_output_process;k++) {
        if (categorical_process[k]) {
          buff = categorical_process[k]->observation[state[j][0]]->cumul[*pioutput[k]];
        }

        else if (discrete_parametric_process[k]) {
          buff = discrete_parametric_process[k]->observation[state[j][0]]->cumul[*pioutput[k]];
        }

        else {
          if (((continuous_parametric_process[k]->ident == GAMMA) ||
              (continuous_parametric_process[k]->ident == ZERO_INFLATED_GAMMA)) && (seq.min_value[k + 1] < seq.min_interval[k + 1] / 2)) {
            switch (seq.type[k + 1]) {
            case INT_VALUE :
              buff = continuous_parametric_process[k]->observation[state[j][0]]->mass_computation(*pioutput[k] , *pioutput[k] + seq.min_interval[k + 1]);
              break;
            case REAL_VALUE :
              buff = continuous_parametric_process[k]->observation[state[j][0]]->mass_computation(*proutput[k] , *proutput[k] + seq.min_interval[k + 1]);
              break;
            }
          }

          else {
            switch (seq.type[k + 1]) {
            case INT_VALUE :
              buff = continuous_parametric_process[k]->observation[state[j][0]]->mass_computation(*pioutput[k] - seq.min_interval[k + 1] / 2 , *pioutput[k] + seq.min_interval[k + 1] / 2);
              break;
            case REAL_VALUE :
              buff = continuous_parametric_process[k]->observation[state[j][0]]->mass_computation(*proutput[k] - seq.min_interval[k + 1] / 2 , *proutput[k] + seq.min_interval[k + 1] / 2);
              break;
            }
          }

          if (buff > 0.) {
            buff = log(buff);
          }
          else {
            buff = D_INF;
          }
        }

        if (buff == D_INF) {
          observation = D_INF;
          break;
        }
        else {
          observation += buff;
        }
      }

      for (k = 1;k < nb_row;k++) {
        rank[k] = 0;
      }

      for (k = 0;k < nb_state_sequence;k++) {
        forward[j][k] = D_INF;
        for (m = 0;m < nb_memory[j];m++) {
          buff = cumul_transition[previous[j][m]][state[j][0]] +
                 previous_forward[previous[j][m]][rank[previous[j][m]]];
          if (buff > forward[j][k]) {
            forward[j][k] = buff;
            optimal_memory[i][j][k] = previous[j][m];
            optimal_rank[i][j][k] = rank[previous[j][m]];
          }
        }

        if (forward[j][k] != D_INF) {
          rank[optimal_memory[i][j][k]]++;

/*          for (m = 0;m < nb_output_process;m++) {
            if (categorical_process[m]) {
              buff = categorical_process[m]->observation[state[j][0]]->cumul[*pioutput[m]];
            }

            else if (discrete_parametric_process[m]) {
              buff = discrete_parametric_process[m]->observation[state[j][0]]->cumul[*pioutput[m]];
            }

            else {
              if (((continuous_parametric_process[m]->ident == GAMMA) ||
                  (continuous_parametric_process[m]->ident == ZERO_INFLATED_GAMMA)) && (seq.min_value[m + 1] < seq.min_interval[m + 1] / 2)) {
                switch (seq.type[m + 1]) {
                case INT_VALUE :
                  buff = continuous_parametric_process[m]->observation[state[j][0]]->mass_computation(*pioutput[m] , *pioutput[m] + seq.min_interval[m + 1]);
                  break;
                case REAL_VALUE :
                  buff = continuous_parametric_process[m]->observation[state[j][0]]->mass_computation(*proutput[m] , *proutput[m] + seq.min_interval[m + 1]);
                  break;
                }
              }

              else {
                switch (seq.type[m + 1]) {
                case INT_VALUE :
                  buff = continuous_parametric_process[m]->observation[state[j][0]]->mass_computation(*pioutput[m] - seq.min_interval[m + 1] / 2 , *pioutput[m] + seq.min_interval[m + 1] / 2);
                  break;
                case REAL_VALUE :
                  buff = continuous_parametric_process[m]->observation[state[j][0]]->mass_computation(*proutput[m] - seq.min_interval[m + 1] / 2 , *proutput[m] + seq.min_interval[m + 1] / 2);
                  break;
                }
              }

              if (buff > 0.) {
                buff = log(buff);
              }
              else {
                buff = D_INF;
              }
            }

            if (buff == D_INF) {
              forward[j][k] = D_INF;
              break;
            }
            else {
              forward[j][k] += buff;
            }
          } */

          if (observation == D_INF) {
            forward[j][k] = D_INF;
            break;
          }
          else {
            forward[j][k] += observation;
          }
        }
      }
    }
  }

  // extraction de la vraisemblance du chemin optimal

  for (i = 1;i < nb_row;i++) {
    rank[i] = 0;
  }
  likelihood_cumul = 0.;

  for (i = 0;i < nb_state_sequence;i++) {
    pstate = seq.int_sequence[index][0] + seq.length[index] - 1;
    forward_max = D_INF;

    for (j = 1;j < nb_row;j++) {
      if (forward[j][rank[j]] > forward_max) {
        forward_max = forward[j][rank[j]];
        memory = j;
      }
    }

    if (i == 0) {
      state_seq_likelihood = forward_max;
    }

    if (forward_max == D_INF) {
      break;
    }

    // restauration

    *pstate = state[memory][0];
    active_cell[seq.length[index] - 1][*pstate] = true;
    brank = rank[memory];
    rank[memory]++;

#   ifdef DEBUG
    cout << "\n" << *pstate << " " << brank << " | ";
#   endif

    for (j = seq.length[index] - 1;j > 0;j--) {
      previous_rank = optimal_rank[j][memory][brank];
      memory = optimal_memory[j][memory][brank];
      *--pstate = state[memory][0];
      active_cell[j - 1][*pstate] = true;
      brank = previous_rank;

#     ifdef DEBUG
      cout << *pstate << " " << brank << " | ";
#     endif
    }

#   ifdef DEBUG
    cout << endl;
#   endif

    likelihood_cumul += exp(forward_max);

#   ifdef DEBUG
    pstate = seq.int_sequence[index][0];
    for (j = 0;j < seq.length[index];j++) {
/*      state_sequence_probability[j][*pstate++] += exp(forward_max - seq_likelihood); */

      if (forward_max > state_sequence_probability[j][*pstate]) {
        state_sequence_probability[j][*pstate] = forward_max;
      }
      pstate++;
    }
#   endif

    nb_cell = 0;
    for (j = 0;j < seq.length[index];j++) {
      for (k = 0;k < nb_state;k++) {
        if (active_cell[j][k]) {
          nb_cell++;
        }
      }
    }

#   ifdef MESSAGE
    if (i == 0) {
      os << "\n";
    }

    pstate = seq.int_sequence[index][0];

    switch (format) {

    case 'a' : {
      for (j = 0;j < seq.length[index];j++) {
        os << *pstate++ << " ";
      }

//      os << "  " << i + 1 << "  " << forward_max << "   (" << exp(forward_max - state_seq_likelihood)
      os << "  " << i + 1 << "  " << forward_max << "   (" << exp(forward_max - seq_likelihood)
         << "  " << likelihood_cumul / exp(seq_likelihood) << "  " << nb_cell << ")" << endl;
      break;
    }

    case 's' : {
      for (j = 0;j < seq.length[index];j++) {
        os << *pstate++ << "\t";
      }

//      os << "\t" << i + 1 << "\t" << forward_max << "\t" << exp(forward_max - state_seq_likelihood)
      os << "\t" << i + 1 << "\t" << forward_max << "\t" << exp(forward_max - seq_likelihood)
         << "\t" << likelihood_cumul / exp(seq_likelihood) << "\t" << nb_cell << endl;
      break;
    }
    }
#   endif

#   ifdef DEBUG
    entropy -= exp(forward_max - seq_likelihood) * forward_max;
#   endif

  }

# ifdef DEBUG
  os << "\n" << SEQ_label[SEQL_STATE_SEQUENCE_ENTROPY] << ": " << entropy + seq_likelihood << endl;

  if (likelihood_cumul / exp(seq_likelihood) > 0.8) {
    for (i = 0;i < seq.length[index];i++) {
      for (j = 0;j < nb_state;j++) {
        if (state_sequence_probability[i][j] != D_INF) {
          state_sequence_probability[i][j] = exp(state_sequence_probability[i][j] - seq_likelihood);
        }
        else {
          state_sequence_probability[i][j] = 0.;
        }
      }
    }

    pstate = seq.int_sequence[index][0];
    for (j = 0;j < seq.length[index];j++) {
      *pstate++ = I_DEFAULT;
    }

//    os << "\n" << SEQ_label[SEQL_POSTERIOR_STATE_PROBABILITY] << "\n\n";
    os << "\n" << SEQ_label[SEQL_MAX_POSTERIOR_STATE_PROBABILITY] << "\n\n";
    seq.profile_ascii_print(os , index , nb_state , state_sequence_probability ,
                            STAT_label[STATL_STATE]);
  }
# endif

  for (i = 1;i < nb_row;i++) {
    delete [] forward[i];
  }
  delete [] forward;

  for (i = 1;i < nb_row;i++) {
    delete [] previous_forward[i];
  }
  delete [] previous_forward;

  delete [] rank;

  for (i = 0;i < seq.length[index];i++) {
    for (j = 1;j < nb_row;j++) {
      delete [] optimal_memory[i][j];
    }
    delete [] optimal_memory[i];
  }
  delete [] optimal_memory;

  for (i = 0;i < seq.length[index];i++) {
    for (j = 1;j < nb_row;j++) {
      delete [] optimal_rank[i][j];
    }
    delete [] optimal_rank[i];
  }
  delete [] optimal_rank;

  for (i = 0;i < seq.length[index];i++) {
    delete [] active_cell[i];
  }
  delete [] active_cell;

  delete [] pioutput;
  delete [] proutput;

# ifdef DEBUG
  for (i = 0;i < seq.length[index];i++) {
    delete [] state_sequence_probability[i];
  }
  delete [] state_sequence_probability;
# endif

  return state_seq_likelihood;
}


/*--------------------------------------------------------------*
 *
 *  Calcul des profils d'etats par l'algorithme de Viterbi forward-backward.
 *
 *  arguments : reference sur un objet MarkovianSequences, indice de la sequence,
 *              stream, pointeur sur un objet MultiPlot, format de sortie ('a' : ASCII,
 *              's' : Spreadsheet, 'g' : Gnuplot, 'p' : plotable), vraisemblance des donnees.
 *
 *--------------------------------------------------------------*/

double HiddenVariableOrderMarkov::viterbi_forward_backward(const MarkovianSequences &seq ,
                                                           int index , ostream *os , MultiPlot *plot ,
                                                           char format , double seq_likelihood) const

{
  register int i , j , k;
  int *pstate , **pioutput;
  double buff , state_seq_likelihood , backward_max , **forward , **backward , *auxiliary ,
         **state_backward , **proutput;


  // initialisations

  forward = new double*[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    forward[i] = new double[nb_row];
  }

  backward = new double*[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    backward[i] = new double[nb_row];
  }

  auxiliary = new double[nb_row];

  state_backward = new double*[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    state_backward[i] = new double[nb_state];
  }

  pioutput = new int*[nb_output_process];
  proutput = new double*[nb_output_process];

# ifdef MESSAGE
  int memory , *state_sequence , **optimal_memory;

  optimal_memory = new int*[seq.length[index]];
  for (i = 0;i < seq.length[index];i++) {
    optimal_memory[i] = new int[nb_row];
  }

  state_sequence = new int[seq.length[index]];
# endif

  for (i = 0;i < nb_output_process;i++) {
    switch (seq.type[i + 1]) {
    case INT_VALUE :
      pioutput[i] = seq.int_sequence[index][i + 1];
      break;
    case REAL_VALUE :
      proutput[i] = seq.real_sequence[index][i + 1];
      break;
    }
  }

  // recurrence "forward"

  switch (type) {

  case 'o' : {
    for (i = 1;i < nb_row;i++) {
      if (order[i] == 1) {
        forward[0][i] = cumul_initial[state[i][0]];

        if (forward[0][i] != D_INF) {
          for (j = 0;j < nb_output_process;j++) {
            if (categorical_process[j]) {
              buff = categorical_process[j]->observation[state[i][0]]->cumul[*pioutput[j]];
            }
            else if (discrete_parametric_process[j]) {
              buff = discrete_parametric_process[j]->observation[state[i][0]]->cumul[*pioutput[j]];
            }

            else {
              if (((continuous_parametric_process[j]->ident == GAMMA) ||
                  (continuous_parametric_process[j]->ident == ZERO_INFLATED_GAMMA)) && (seq.min_value[j + 1] < seq.min_interval[j + 1] / 2)) {
                switch (seq.type[j + 1]) {
                case INT_VALUE :
                  buff = continuous_parametric_process[j]->observation[state[i][0]]->mass_computation(*pioutput[j] , *pioutput[j] + seq.min_interval[j + 1]);
                  break;
                case REAL_VALUE :
                  buff = continuous_parametric_process[j]->observation[state[i][0]]->mass_computation(*proutput[j] , *proutput[j] + seq.min_interval[j + 1]);
                  break;
                }
              }

              else {
                switch (seq.type[j + 1]) {
                case INT_VALUE :
                  buff = continuous_parametric_process[j]->observation[state[i][0]]->mass_computation(*pioutput[j] - seq.min_interval[j + 1] / 2 , *pioutput[j] + seq.min_interval[j + 1] / 2);
                  break;
                case REAL_VALUE :
                  buff = continuous_parametric_process[j]->observation[state[i][0]]->mass_computation(*proutput[j] - seq.min_interval[j + 1] / 2 , *proutput[j] + seq.min_interval[j + 1] / 2);
                  break;
                }
              }

              if (buff > 0.) {
                buff = log(buff);
              }
              else {
                buff = D_INF;
              }
            }

            if (buff == D_INF) {
              forward[0][i] = D_INF;
              break;
            }
            else {
              forward[0][i] += buff;
            }
          }
        }
      }

      else {
        forward[0][i] = D_INF;
      }
    }
    break;
  }

  case 'e' : {
    for (i = 1;i < nb_row;i++) {
      if (!child[i]) {
        forward[0][i] = cumul_initial[i];

        if (forward[0][i] != D_INF) {
          for (j = 0;j < nb_output_process;j++) {
            if (categorical_process[j]) {
              buff = categorical_process[j]->observation[state[i][0]]->cumul[*pioutput[j]];
            }

            else if (discrete_parametric_process[j]) {
              buff = discrete_parametric_process[j]->observation[state[i][0]]->cumul[*pioutput[j]];
            }

            else {
              if (((continuous_parametric_process[j]->ident == GAMMA) ||
                  (continuous_parametric_process[j]->ident == ZERO_INFLATED_GAMMA)) && (seq.min_value[j + 1] < seq.min_interval[j + 1] / 2)) {
                switch (seq.type[j + 1]) {
                case INT_VALUE :
                  buff = continuous_parametric_process[j]->observation[state[i][0]]->mass_computation(*pioutput[j] , *pioutput[j] + seq.min_interval[j + 1]);
                  break;
                case REAL_VALUE :
                  buff = continuous_parametric_process[j]->observation[state[i][0]]->mass_computation(*proutput[j] , *proutput[j] + seq.min_interval[j + 1]);
                  break;
                }
              }

              else {
                switch (seq.type[j + 1]) {
                case INT_VALUE :
                  buff = continuous_parametric_process[j]->observation[state[i][0]]->mass_computation(*pioutput[j] - seq.min_interval[j + 1] / 2 , *pioutput[j] + seq.min_interval[j + 1] / 2);
                  break;
                case REAL_VALUE :
                  buff = continuous_parametric_process[j]->observation[state[i][0]]->mass_computation(*proutput[j] - seq.min_interval[j + 1] / 2 , *proutput[j] + seq.min_interval[j + 1] / 2);
                  break;
                }
              }

              if (buff > 0.) {
                buff = log(buff);
              }
              else {
                buff = D_INF;
              }
            }

            if (buff == D_INF) {
              forward[0][i] = D_INF;
              break;
            }
            else {
              forward[0][i] += buff;
            }
          }
        }
      }

      else {
        forward[0][i] = D_INF;
      }
    }
    break;
  }
  }

  for (i = 1;i < seq.length[index];i++) {
    for (j = 0;j < nb_output_process;j++) {
      switch (seq.type[j + 1]) {
      case INT_VALUE :
        pioutput[j]++;
        break;
      case REAL_VALUE :
        proutput[j]++;
        break;
      }
    }

    for (j = 1;j < nb_row;j++) {
      forward[i][j] = D_INF;
      for (k = 0;k < nb_memory[j];k++) {
        buff = cumul_transition[previous[j][k]][state[j][0]] + forward[i - 1][previous[j][k]];
        if (buff > forward[i][j]) {
          forward[i][j] = buff;

#         ifdef MESSAGE
          optimal_memory[i][j] = previous[j][k];
#         endif
        }
      }

      if (forward[i][j] != D_INF) {
        for (k = 0;k < nb_output_process;k++) {
          if (categorical_process[k]) {
            buff = categorical_process[k]->observation[state[j][0]]->cumul[*pioutput[k]];
          }

          else if (discrete_parametric_process[k]) {
            buff = discrete_parametric_process[k]->observation[state[j][0]]->cumul[*pioutput[k]];
          }

          else {
            if (((continuous_parametric_process[k]->ident == GAMMA) ||
                (continuous_parametric_process[k]->ident == ZERO_INFLATED_GAMMA)) && (seq.min_value[k + 1] < seq.min_interval[k + 1] / 2)) {
              switch (seq.type[k + 1]) {
              case INT_VALUE :
                buff = continuous_parametric_process[k]->observation[state[j][0]]->mass_computation(*pioutput[k] , *pioutput[k] + seq.min_interval[k + 1]);
                break;
              case REAL_VALUE :
                buff = continuous_parametric_process[k]->observation[state[j][0]]->mass_computation(*proutput[k] , *proutput[k] + seq.min_interval[k + 1]);
                break;
              }
            }

            else {
              switch (seq.type[k + 1]) {
              case INT_VALUE :
                buff = continuous_parametric_process[k]->observation[state[j][0]]->mass_computation(*pioutput[k] - seq.min_interval[k + 1] / 2 , *pioutput[k] + seq.min_interval[k + 1] / 2);
                break;
              case REAL_VALUE :
                buff = continuous_parametric_process[k]->observation[state[j][0]]->mass_computation(*proutput[k] - seq.min_interval[k + 1] / 2 , *proutput[k] + seq.min_interval[k + 1] / 2);
                break;
              }
            }

            if (buff > 0.) {
              buff = log(buff);
            }
            else {
              buff = D_INF;
            }
          }

          if (buff == D_INF) {
            forward[i][j] = D_INF;
            break;
          }
          else {
            forward[i][j] += buff;
          }
        }
      }
    }
  }

  // extraction de la vraisemblance du chemin optimal

# ifdef MESSAGE
  pstate = state_sequence + seq.length[index] - 1;
# endif

  state_seq_likelihood = D_INF;
  i = seq.length[index] - 1;
  for (j = 1;j < nb_row;j++) {
    if (forward[i][j] > state_seq_likelihood) {
      state_seq_likelihood = forward[i][j];

#     ifdef MESSAGE
      memory = j;
#     endif
    }
  }

  if (state_seq_likelihood != D_INF) {

#   ifdef MESSAGE
    *pstate = state[memory][0];
    for (i = seq.length[index] - 1;i > 0;i--) {
      memory = optimal_memory[i][memory];
      *--pstate = state[memory][0];
    }
#   endif

    // recurrence "backward"

    i = seq.length[index] - 1;
    for (j = 1;j < nb_row;j++) {
      backward[i][j] = 0.;
    }

    for (i = seq.length[index] - 2;i >= 0;i--) {
      for (j = 1;j < nb_row;j++) {
        auxiliary[j] = backward[i + 1][j];

        for (k = 0;k < nb_output_process;k++) {
          if (categorical_process[k]) {
            auxiliary[j] += categorical_process[k]->observation[state[j][0]]->cumul[*pioutput[k]];
          }

          else if (discrete_parametric_process[k]) {
            auxiliary[j] += discrete_parametric_process[k]->observation[state[j][0]]->cumul[*pioutput[k]];
          }

          else {
            if (((continuous_parametric_process[k]->ident == GAMMA) ||
                (continuous_parametric_process[k]->ident == ZERO_INFLATED_GAMMA)) && (seq.min_value[k + 1] < seq.min_interval[k + 1] / 2)) {
              switch (seq.type[k + 1]) {
              case INT_VALUE :
                 buff = continuous_parametric_process[k]->observation[state[j][0]]->mass_computation(*pioutput[k] , *pioutput[k] + seq.min_interval[k + 1]);
                break;
              case REAL_VALUE :
                buff  = continuous_parametric_process[k]->observation[state[j][0]]->mass_computation(*proutput[k] , *proutput[k] + seq.min_interval[k + 1]);
                break;
              }
            }

            else {
              switch (seq.type[k + 1]) {
              case INT_VALUE :
                 buff = continuous_parametric_process[k]->observation[state[j][0]]->mass_computation(*pioutput[k] - seq.min_interval[k + 1] / 2 , *pioutput[k] + seq.min_interval[k + 1] / 2);
                break;
              case REAL_VALUE :
                buff  = continuous_parametric_process[k]->observation[state[j][0]]->mass_computation(*proutput[k] - seq.min_interval[k + 1] / 2 , *proutput[k] + seq.min_interval[k + 1] / 2);
                break;
              }
            }

            if (buff > 0.) {
              auxiliary[j] += log(buff);
            }
            else {
              auxiliary[j] = D_INF;
            }
          }
        }
      }

      for (j = 0;j < nb_output_process;j++) {
        switch (seq.type[j + 1]) {
        case INT_VALUE :
          pioutput[j]--;
          break;
        case REAL_VALUE :
          proutput[j]--;
          break;
        }
      }

      for (j = 1;j < nb_row;j++) {
        backward[i][j] = D_INF;
        if (next[j]) {
          for (k = 0;k < nb_state;k++) {
            buff = auxiliary[next[j][k]] + cumul_transition[j][k];
            if (buff > backward[i][j]) {
              backward[i][j] = buff;
            }
          }
        }
      }
    }

    // restauration

    pstate = seq.int_sequence[index][0];

    for (i = 0;i < seq.length[index];i++) {
      backward_max = D_INF;
      for (j = 1;j < nb_row;j++) {
        if (backward[i][j] != D_INF) {
          if (forward[i][j] != D_INF) {
            backward[i][j] += forward[i][j];
            if (backward[i][j] > backward_max) {
              backward_max = backward[i][j];
              *pstate = state[j][0];
            }
          }

          else {
            backward[i][j] = D_INF;
          }
        }
      }

#     ifdef MESSAGE
      if (*pstate != state_sequence[i]) {
        cout << "\nERROR: " << i << " | " << *pstate << " " << state_sequence[i] << endl;
      }
#     endif

      pstate++;
    }

    //  normalisation

    for (i = 0;i < seq.length[index];i++) {
      for (j = 0;j < nb_state;j++) {
        state_backward[i][j] = D_INF;
      }
      for (j = 1;j < nb_row;j++) {
        if (backward[i][j] > state_backward[i][state[j][0]]) {
          state_backward[i][state[j][0]] = backward[i][j];
        }
      }

      for (j = 0;j < nb_state;j++) {
        if (state_backward[i][j] != D_INF) {
          state_backward[i][j] = exp(state_backward[i][j] - seq_likelihood);
//          state_backward[i][j] = exp(state_backward[i][j] - state_seq_likelihood);
        }
        else {
          state_backward[i][j] = 0.;
        }
      }
    }

    switch (format) {

    case 'a' : {
      *os << "\n" << SEQ_label[SEQL_MAX_POSTERIOR_STATE_PROBABILITY] << "\n\n";
      seq.profile_ascii_print(*os , index , nb_state , state_backward ,
                              STAT_label[STATL_STATE]);

      *os << "\n" << SEQ_label[SEQL_STATE_SEQUENCE_LIKELIHOOD] << ": " << state_seq_likelihood
          << "   (" << exp(state_seq_likelihood - seq_likelihood) << ")" << endl;
      break;
    }

    case 's' : {
      *os << "\n" << SEQ_label[SEQL_MAX_POSTERIOR_STATE_PROBABILITY] << "\n\n";
      seq.profile_spreadsheet_print(*os , index , nb_state , state_backward ,
                                    STAT_label[STATL_STATE]);

      *os << "\n" << SEQ_label[SEQL_STATE_SEQUENCE_LIKELIHOOD] << "\t" << state_seq_likelihood
          << "\t" << exp(state_seq_likelihood - seq_likelihood) << endl;
      break;
    }

    case 'g' : {
      seq.profile_plot_print(*os , index , nb_state , state_backward);
      break;
    }

    case 'p' : {
      seq.profile_plotable_write(*plot , index , nb_state , state_backward);
      break;
    }
    }

#   ifdef DEBUG
    if (format != 'g') {
      double ambiguity = 0.;

      pstate = seq.int_sequence[index][0];
      for (i = 0;i < seq.length[index];i++) {
        for (j = 0;j < nb_state;j++) {
          if (j != *pstate) {
            ambiguity += state_backward[i][j];
          }
        }
        pstate++;
      }
      ambiguity *= exp(seq_likelihood - state_seq_likelihood);

      switch (format) {
      case 'a' :
        *os << "\n" << SEQ_label[SEQL_AMBIGUITY] << ": " << ambiguity
            << " (" << ambiguity / seq.length[index] << ")" << endl;
        break;
      case 's' :
        *os << "\n" << SEQ_label[SEQL_AMBIGUITY] << "\t" << ambiguity
            << "\t" << ambiguity / seq.length[index] << "\t" << endl;
        break;
      }
    }
#   endif

  }

  for (i = 0;i < seq.length[index];i++) {
    delete [] forward[i];
  }
  delete [] forward;

  for (i = 0;i < seq.length[index];i++) {
    delete [] backward[i];
  }
  delete [] backward;

  delete [] auxiliary;

  for (i = 0;i < seq.length[index];i++) {
    delete [] state_backward[i];
  }
  delete [] state_backward;

  delete [] pioutput;
  delete [] proutput;

# ifdef MESSAGE
  for (i = 0;i < seq.length[index];i++) {
    delete [] optimal_memory[i];
  }
  delete [] optimal_memory;

  delete [] state_sequence;
# endif

  return state_seq_likelihood;
}


/*--------------------------------------------------------------*
 *
 *  Calcul des profils d'etats et d'entropie par l'algorithme forward-backward,
 *  des profils d'etats par l'algorithme de Viterbi forward-backward,
 *  calcul des N sequences d'etats les plus probables par l'algorithme de Viterbi generalise ou
 *  simulation de sequences d'etats par l'algorithme forward-backward de simulation et
 *  ecriture des resultats.
 *
 *  arguments : reference sur un objet StatError, stream, sequences,
 *              identificateur de la sequence, format ('a' : ASCII, 's' : Spreadsheet),
 *              methode de calcul des sequences d'etats (algorithme de Viterbi generalise ou
 *              algorithme forward-backward de simulation), nombre de sequences d'etats.
 *
 *--------------------------------------------------------------*/

bool HiddenVariableOrderMarkov::state_profile_write(StatError &error , ostream &os ,
                                                    const MarkovianSequences &iseq ,
                                                    int identifier , char format ,
                                                    int state_sequence , int nb_state_sequence) const

{
  bool status = true;
  register int i;
  int offset = I_DEFAULT , nb_value , index = I_DEFAULT;
  double seq_likelihood , max_marginal_entropy , entropy;
  HiddenVariableOrderMarkov *hmarkov;
  VariableOrderMarkovData *seq;


  error.init();

  for (i = 0;i < iseq.nb_variable;i++) {
    if ((iseq.type[i] != INT_VALUE) && (iseq.type[i] != REAL_VALUE) && (iseq.type[i] != STATE)) {
      status = false;
      ostringstream error_message , correction_message;
      error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                    << STAT_error[STATR_VARIABLE_TYPE];
      correction_message << STAT_variable_word[INT_VALUE] << " or "
                         << STAT_variable_word[REAL_VALUE];
      error.correction_update((error_message.str()).c_str() , (correction_message.str()).c_str());
    }
  }

  if (nb_output_process == iseq.nb_variable) {
    offset = 0;
  }
  else if ((iseq.type[0] == STATE) && (nb_output_process + 1 == iseq.nb_variable)) {
    offset = 1;
  }
  else {
    status = false;
    error.update(STAT_error[STATR_NB_OUTPUT_PROCESS]);
  }

  if (offset != I_DEFAULT) {
    for (i = 0;i < nb_output_process;i++) {
      if ((categorical_process[i]) || (discrete_parametric_process[i])) {
        if (iseq.type[i + offset] == REAL_VALUE) {
          status = false;
          ostringstream error_message;
          error_message << STAT_label[STATL_VARIABLE] << " " << i + offset + 1 << ": "
                        << STAT_error[STATR_VARIABLE_TYPE];
          error.correction_update((error_message.str()).c_str() , STAT_variable_word[INT_VALUE]);
        }

        else {
          if (iseq.min_value[i + offset] < 0) {
            status = false;
            ostringstream error_message;
            error_message << STAT_label[STATL_VARIABLE] << " " << i + offset + 1 << ": "
                          << STAT_error[STATR_POSITIVE_MIN_VALUE];
            error.update((error_message.str()).c_str());
          }

          if (!(iseq.marginal_distribution[i + offset])) {
            status = false;
            ostringstream error_message;
            error_message << STAT_label[STATL_VARIABLE] << " " << i + offset + 1 << ": "
                          << STAT_error[STATR_MARGINAL_FREQUENCY_DISTRIBUTION];
            error.update((error_message.str()).c_str());
          }

          else {
            if (categorical_process[i]) {
              nb_value = categorical_process[i]->nb_value;
            }
            else {
              nb_value = discrete_parametric_process[i]->nb_value;
            }

            if (nb_value < iseq.marginal_distribution[i + offset]->nb_value) {
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

  if (identifier != I_DEFAULT) {
    for (i = 0;i < iseq.nb_sequence;i++) {
      if (identifier == iseq.identifier[i]) {
        index = i;
        break;
      }
    }

    if (i == iseq.nb_sequence) {
      status = false;
      error.update(SEQ_error[SEQR_SEQUENCE_IDENTIFIER]);
    }
  }

  if (nb_state_sequence < 2) {
    status = false;
    error.update(SEQ_error[SEQR_NB_STATE_SEQUENCE]);
  }

  if (status) {
    if (nb_output_process == iseq.nb_variable) {
      seq = new VariableOrderMarkovData(iseq);
    }
    else {
      seq = new VariableOrderMarkovData(iseq , 'c' , (type == 'e' ? true : false));
    }

    hmarkov = new HiddenVariableOrderMarkov(*this , false);
    hmarkov->create_cumul();
    hmarkov->log_computation();

    for (i = 0;i < seq->nb_sequence;i++) {
      if ((index == I_DEFAULT) || (index == i)) {
        seq_likelihood = forward_backward(*seq , i , &os , NULL , format ,
                                          max_marginal_entropy , entropy);

        if (seq_likelihood == D_INF) {
          status = false;

          if (index == I_DEFAULT) {
            ostringstream error_message;
            error_message << SEQ_label[SEQL_SEQUENCE] << " " << i + 1 << " "
                          << SEQ_error[SEQR_INCOMPATIBLE_MODEL];
            error.update((error_message.str()).c_str());
          }
          else {
            error.update(SEQ_error[SEQR_SEQUENCE_INCOMPATIBLE_MODEL]);
          }
        }

        else {
          hmarkov->viterbi_forward_backward(*seq , i , &os , NULL , format ,
                                            seq_likelihood);

          switch (state_sequence) {
          case GENERALIZED_VITERBI :
            hmarkov->generalized_viterbi(*seq , i , os , seq_likelihood , format ,
                                         nb_state_sequence);
            break;
          case FORWARD_BACKWARD_SAMPLING :
            forward_backward_sampling(*seq , i , os , format , nb_state_sequence);
            break;
          }
        }
      }
    }

    delete seq;
    delete hmarkov;
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Calcul des profils d'etats et d'entropie par l'algorithme forward-backward,
 *  des profils d'etats par l'algorithme de Viterbi forward-backward,
 *  calcul des N sequences d'etats les plus probables par l'algorithme de Viterbi generalise ou
 *  simulation de sequences d'etats par l'algorithme forward-backward de simulation et
 *  ecriture des resultats dans un fichier.
 *
 *  arguments : reference sur un objet StatError, path, sequences,
 *              identificateur de la sequence, format de fichier ('a' : ASCII, 's' : Spreadsheet),
 *              methode de calcul des sequences d'etats (algorithme de Viterbi generalise ou
 *              algorithme forward-backward de simulation), nombre de sequences d'etats.
 *
 *--------------------------------------------------------------*/

bool HiddenVariableOrderMarkov::state_profile_write(StatError &error , const char *path ,
                                                    const MarkovianSequences &iseq ,
                                                    int identifier , char format ,
                                                    int state_sequence , int nb_state_sequence) const

{
  bool status = true;
  ofstream out_file(path);


  error.init();

  if (!out_file) {
    status = false;
    error.update(STAT_error[STATR_FILE_NAME]);
  }
  else {
    status = state_profile_write(error , out_file , iseq , identifier ,
                                 format , state_sequence , nb_state_sequence);
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Calcul des profils d'etats et d'entropie par l'algorithme forward-backward,
 *  des profils d'etats par l'algorithme de Viterbi forward-backward,
 *  calcul des N sequences d'etats les plus probables par l'algorithme de Viterbi generalise ou
 *  simulation de sequences d'etats par l'algorithme forward-backward de simulation et
 *  ecriture des resultats.
 *
 *  arguments : reference sur un objet StatError, stream, identificateur de la sequence,
 *              methode de calcul des sequences d'etats (algorithme de Viterbi generalise ou
 *              algorithme forward-backward de simulation), nombre de sequences d'etats.
 *
 *--------------------------------------------------------------*/

bool HiddenVariableOrderMarkov::state_profile_ascii_write(StatError &error , ostream &os ,
                                                          int identifier , int state_sequence ,
                                                          int nb_state_sequence) const

{
  bool status;


  error.init();

  if (!markov_data) {
    status = false;
    error.update(STAT_error[STATR_NO_DATA]);
  }
  else {
    status = state_profile_write(error , os , *markov_data , identifier , 'a' ,
                                 state_sequence , nb_state_sequence);
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Calcul des profils d'etats et d'entropie par l'algorithme forward-backward,
 *  des profils d'etats par l'algorithme de Viterbi forward-backward,
 *  calcul des N sequences d'etats les plus probables par l'algorithme de Viterbi generalise ou
 *  simulation de sequences d'etats par l'algorithme forward-backward de simulation et
 *  ecriture des resultats dans un fichier.
 *
 *  arguments : reference sur un objet StatError, path, identificateur de la sequence,
 *              format de fichier ('a' : ASCII, 's' : Spreadsheet),
 *              methode de calcul des sequences d'etats (algorithme de Viterbi generalise ou
 *              algorithme forward-backward de simulation), nombre de sequences d'etats.
 *
 *--------------------------------------------------------------*/

bool HiddenVariableOrderMarkov::state_profile_write(StatError &error , const char *path ,
                                                    int identifier , char format ,
                                                    int state_sequence , int nb_state_sequence) const

{
  bool status = true;
  ofstream out_file(path);


  error.init();

  if (!out_file) {
    status = false;
    error.update(STAT_error[STATR_FILE_NAME]);
  }
  if (!markov_data) {
    status = false;
    error.update(STAT_error[STATR_NO_DATA]);
  }

  if (status) {
    status = state_profile_write(error , out_file , *markov_data , identifier ,
                                 format , state_sequence , nb_state_sequence);
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Calcul des profils d'etats et d'entropie par l'algorithme forward-backward,
 *  des profils d'etats par l'algorithme de Viterbi forward-backward et
 *  affichage des resultats au format Gnuplot.
 *
 *  arguments : reference sur un objet StatError, prefixe des fichiers,
 *              sequences, identificateur de la sequence, titre des figures.
 *
 *--------------------------------------------------------------*/

bool HiddenVariableOrderMarkov::state_profile_plot_write(StatError &error , const char *prefix ,
                                                         const MarkovianSequences &iseq ,
                                                         int identifier , const char *title) const

{
  bool status = true;
  register int i , j;
  int offset = I_DEFAULT , nb_value , index;
  double seq_likelihood , max_marginal_entropy , entropy , state_seq_likelihood;
  HiddenVariableOrderMarkov *hmarkov;
  VariableOrderMarkovData *seq;
  ostringstream data_file_name[2];
  ofstream *out_data_file;


  error.init();

  for (i = 0;i < iseq.nb_variable;i++) {
    if ((iseq.type[i] != INT_VALUE) && (iseq.type[i] != REAL_VALUE) && (iseq.type[i] != STATE)) {
      status = false;
      ostringstream error_message , correction_message;
      error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                    << STAT_error[STATR_VARIABLE_TYPE];
      correction_message << STAT_variable_word[INT_VALUE] << " or "
                         << STAT_variable_word[REAL_VALUE];
      error.correction_update((error_message.str()).c_str() , (correction_message.str()).c_str());
    }
  }

  if (nb_output_process == iseq.nb_variable) {
    offset = 0;
  }
  else if ((iseq.type[0] == STATE) && (nb_output_process + 1 == iseq.nb_variable)) {
    offset = 1;
  }
  else {
    status = false;
    error.update(STAT_error[STATR_NB_OUTPUT_PROCESS]);
  }

  if (offset != I_DEFAULT) {
    for (i = 0;i < nb_output_process;i++) {
      if ((categorical_process[i]) || (discrete_parametric_process[i])) {
        if (iseq.type[i + offset] == REAL_VALUE) {
          status = false;
          ostringstream error_message;
          error_message << STAT_label[STATL_VARIABLE] << " " << i + offset + 1 << ": "
                        << STAT_error[STATR_VARIABLE_TYPE];
          error.correction_update((error_message.str()).c_str() , STAT_variable_word[INT_VALUE]);
        }

        else {
          if (iseq.min_value[i + offset] < 0) {
            status = false;
            ostringstream error_message;
            error_message << STAT_label[STATL_VARIABLE] << " " << i + offset + 1 << ": "
                          << STAT_error[STATR_POSITIVE_MIN_VALUE];
            error.update((error_message.str()).c_str());
          }

          if (!(iseq.marginal_distribution[i + offset])) {
            status = false;
            ostringstream error_message;
            error_message << STAT_label[STATL_VARIABLE] << " " << i + offset + 1 << ": "
                          << STAT_error[STATR_MARGINAL_FREQUENCY_DISTRIBUTION];
            error.update((error_message.str()).c_str());
          }

          else {
            if (categorical_process[i]) {
              nb_value = categorical_process[i]->nb_value;
            }
            else {
              nb_value = discrete_parametric_process[i]->nb_value;
            }

            if (nb_value < iseq.marginal_distribution[i + offset]->nb_value) {
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

  for (i = 0;i < iseq.nb_sequence;i++) {
    if (identifier == iseq.identifier[i]) {
      index = i;
      break;
    }
  }

  if (i == iseq.nb_sequence) {
    status = false;
    error.update(SEQ_error[SEQR_SEQUENCE_IDENTIFIER]);
  }

  if (status) {

    // ecriture des fichiers de donnees

    data_file_name[0] << prefix << 0 << ".dat";
    out_data_file = new ofstream((data_file_name[0].str()).c_str());

    if (!out_data_file) {
      status = false;
      error.update(STAT_error[STATR_FILE_PREFIX]);
    }

    else {
      if (iseq.type[0] != STATE) {
        seq = new VariableOrderMarkovData(iseq);
      }
      else {
        seq = new VariableOrderMarkovData(iseq , 'c' , (type == 'e' ? true : false));
      }

      seq_likelihood = forward_backward(*seq , index , out_data_file , NULL , 'g' ,
                                        max_marginal_entropy , entropy);
      out_data_file->close();
      delete out_data_file;

      if (seq_likelihood == D_INF) {
        status = false;
        error.update(SEQ_error[SEQR_SEQUENCE_INCOMPATIBLE_MODEL]);
      }

      else {
        data_file_name[1] << prefix << 1 << ".dat";
        out_data_file = new ofstream((data_file_name[1].str()).c_str());

        hmarkov = new HiddenVariableOrderMarkov(*this , false);

        hmarkov->create_cumul();
        hmarkov->log_computation();
        state_seq_likelihood = hmarkov->viterbi_forward_backward(*seq , index , out_data_file , NULL ,
                                                                 'g' , seq_likelihood);
        out_data_file->close();
        delete out_data_file;

        // ecriture du fichier de commandes et du fichier d'impression

        for (i = 0;i < 2;i++) {
          ostringstream file_name[2];

          switch (i) {
          case 0 :
            file_name[0] << prefix << ".plot";
            break;
          case 1 :
            file_name[0] << prefix << ".print";
            break;
          }

          ofstream out_file((file_name[0].str()).c_str());

          if (i == 1) {
            out_file << "set terminal postscript" << endl;
            file_name[1] << label(prefix) << ".ps";
            out_file << "set output \"" << file_name[1].str() << "\"\n\n";
          }

          out_file << "set border 15 lw 0\n" << "set tics out\n" << "set xtics nomirror\n"
                   << "set title \"";
          if (title) {
            out_file << title << " - ";
          }
          out_file << SEQ_label[SEQL_MAX_POSTERIOR_STATE_PROBABILITY] << "\"\n\n";

          if (seq->index_parameter) {
            if (seq->index_parameter[index][seq->length[index] - 1] - seq->index_parameter[index][0] < TIC_THRESHOLD) {
              out_file << "set xtics 0,1" << endl;
            }

            out_file << "plot [" << seq->index_parameter[index][0] << ":"
                     << seq->index_parameter[index][seq->length[index] - 1] << "] [0:"
                     << exp(state_seq_likelihood - seq_likelihood) << "] ";
            for (j = 0;j < nb_state;j++) {
              out_file << "\"" << label((data_file_name[1].str()).c_str()) << "\" using "
                       << 1 << " : " << j + 2 << " title \"" << STAT_label[STATL_STATE] << " "
                       << j << "\" with linespoints";
              if (j < nb_state - 1) {
                out_file << ",\\";
              }
              out_file << endl;
            }

            if (i == 0) {
              out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
            }
            out_file << endl;

            out_file << "set title \"";
            if (title) {
              out_file << title << " - ";
            }
            out_file << SEQ_label[SEQL_POSTERIOR_STATE_PROBABILITY] << "\"\n\n";

            out_file << "plot [" << seq->index_parameter[index][0] << ":"
                     << seq->index_parameter[index][seq->length[index] - 1] << "] [0:1] ";
            for (j = 0;j < nb_state;j++) {
              out_file << "\"" << label((data_file_name[0].str()).c_str()) << "\" using "
                       << 1 << " : " << j + 2 << " title \"" << STAT_label[STATL_STATE] << " "
                       << j << "\" with linespoints";
              if (j < nb_state - 1) {
                out_file << ",\\";
              }
              out_file << endl;
            }

            if (i == 0) {
              out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
            }
            out_file << endl;

            out_file << "set title";
            if (title) {
              out_file << " \"" << title << "\"";
            }
            out_file << "\n\n";

            out_file << "plot [" << seq->index_parameter[index][0] << ":"
                     << seq->index_parameter[index][seq->length[index] - 1] << "] [0:"
                     << max_marginal_entropy << "] "
                     << "\"" << label((data_file_name[0].str()).c_str()) << "\" using "
                     << 1 << " : " << nb_state + 2 << " title \"" << SEQ_label[SEQL_CONDITIONAL_ENTROPY]
                     << "\" with linespoints,\\" << endl;
            out_file << "\"" << label((data_file_name[0].str()).c_str()) << "\" using "
                     << 1 << " : " << nb_state + 3 << " title \"" << SEQ_label[SEQL_CONDITIONAL_ENTROPY]
                     << "\" with linespoints,\\" << endl;
            out_file << "\"" << label((data_file_name[0].str()).c_str()) << "\" using "
                     << 1 << " : " << nb_state + 4 << " title \"" << SEQ_label[SEQL_MARGINAL_ENTROPY]
                     << "\" with linespoints" << endl;

            if (i == 0) {
              out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
            }
            out_file << endl;

            out_file << "set title";
            if (title) {
              out_file << " \"" << title << "\"";
            }
            out_file << "\n\n";

            out_file << "plot [" << seq->index_parameter[index][0] << ":"
                     << seq->index_parameter[index][seq->length[index] - 1] << "] [0:" << entropy << "] "
                     << "\"" << label((data_file_name[0].str()).c_str()) << "\" using "
                     << 1 << " : " << nb_state + 5 << " title \""
                     << SEQ_label[SEQL_PARTIAL_STATE_SEQUENCE_ENTROPY] << "\" with linespoints,\\" << endl;
            out_file << "\"" << label((data_file_name[0].str()).c_str()) << "\" using "
                     << 1 << " : " << nb_state + 6 << " title \""
                     << SEQ_label[SEQL_PARTIAL_STATE_SEQUENCE_ENTROPY] << "\" with linespoints" << endl;

            if (seq->index_parameter[index][seq->length[index] - 1] - seq->index_parameter[index][0] < TIC_THRESHOLD) {
              out_file << "set xtics autofreq" << endl;
            }
          }

          else {
            if (seq->length[index] - 1 < TIC_THRESHOLD) {
              out_file << "set xtics 0,1" << endl;
            }

            out_file << "plot [0:" << seq->length[index] - 1 << "] [0:"
                     << exp(state_seq_likelihood - seq_likelihood) << "] ";
            for (j = 0;j < nb_state;j++) {
              out_file << "\"" << label((data_file_name[1].str()).c_str()) << "\" using "
                       << j + 1 << " title \"" << STAT_label[STATL_STATE] << " "
                       << j << "\" with linespoints";
              if (j < nb_state - 1) {
                out_file << ",\\";
              }
              out_file << endl;
            }

            if (i == 0) {
              out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
            }
            out_file << endl;

            out_file << "set title \"";
            if (title) {
              out_file << title << " - ";
            }
            out_file << SEQ_label[SEQL_POSTERIOR_STATE_PROBABILITY] << "\"\n\n";

            out_file << "plot [0:" << seq->length[index] - 1 << "] [0:1] ";
            for (j = 0;j < nb_state;j++) {
              out_file << "\"" << label((data_file_name[0].str()).c_str()) << "\" using "
                       << j + 1 << " title \"" << STAT_label[STATL_STATE] << " "
                       << j << "\" with linespoints";
              if (j < nb_state - 1) {
                out_file << ",\\";
              }
              out_file << endl;
            }

            if (i == 0) {
              out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
            }
            out_file << endl;

            out_file << "set title";
            if (title) {
              out_file << " \"" << title << "\"";
            }
            out_file << "\n\n";

            out_file << "plot [0:" << seq->length[index] - 1 << "] [0:" << max_marginal_entropy << "] "
                     << "\"" << label((data_file_name[0].str()).c_str()) << "\" using "
                     << nb_state + 1 << " title \"" << SEQ_label[SEQL_CONDITIONAL_ENTROPY]
                     << "\" with linespoints,\\" << endl;
            out_file << "\"" << label((data_file_name[0].str()).c_str()) << "\" using "
                     << nb_state + 2 << " title \"" << SEQ_label[SEQL_CONDITIONAL_ENTROPY]
                     << "\" with linespoints,\\" << endl;
            out_file << "\"" << label((data_file_name[0].str()).c_str()) << "\" using "
                     << nb_state + 3 << " title \"" << SEQ_label[SEQL_MARGINAL_ENTROPY]
                     << "\" with linespoints" << endl;

            if (i == 0) {
              out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
            }
            out_file << endl;

            out_file << "set title";
            if (title) {
              out_file << " \"" << title << "\"";
            }
            out_file << "\n\n";

            out_file << "plot [0:" << seq->length[index] - 1 << "] [0:" << entropy << "] "
                     << "\"" << label((data_file_name[0].str()).c_str()) << "\" using "
                     << nb_state + 4 << " title \"" << SEQ_label[SEQL_PARTIAL_STATE_SEQUENCE_ENTROPY]
                     << "\" with linespoints,\\" << endl;
            out_file << "\"" << label((data_file_name[0].str()).c_str()) << "\" using "
                     << nb_state + 5 << " title \"" << SEQ_label[SEQL_PARTIAL_STATE_SEQUENCE_ENTROPY]
                     << "\" with linespoints" << endl;

            if (seq->length[index] - 1 < TIC_THRESHOLD) {
              out_file << "set xtics autofreq" << endl;
            }
          }

          if (i == 1) {
            out_file << "\nset terminal x11" << endl;
          }

          out_file << "\npause 0 \"" << STAT_label[STATL_END] << "\"" << endl;
        }

        delete hmarkov;
      }

      delete seq;
    }
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Calcul des profils d'etats et d'entropie par l'algorithme forward-backward,
 *  des profils d'etats par l'algorithme de Viterbi forward-backward et
 *  affichage des resultats au format Gnuplot.
 *
 *  arguments : reference sur un objet StatError, prefixe des fichiers,
 *              identificateur de la sequence, titre des figures.
 *
 *--------------------------------------------------------------*/

bool HiddenVariableOrderMarkov::state_profile_plot_write(StatError &error ,
                                                         const char *prefix , int identifier ,
                                                         const char *title) const

{
  bool status;


  error.init();

  if (!markov_data) {
    status = false;
    error.update(STAT_error[STATR_NO_DATA]);
  }
  else {
    status = state_profile_plot_write(error , prefix , *markov_data , identifier , title);
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Calcul des profils d'etats et d'entropie par l'algorithme forward-backward,
 *  des profils d'etats par l'algorithme de Viterbi forward-backward et
 *  sortie graphique des resultats.
 *
 *  arguments : reference sur un objet StatError, sequences,
 *              identificateur de la sequence.
 *
 *--------------------------------------------------------------*/

MultiPlotSet* HiddenVariableOrderMarkov::state_profile_plotable_write(StatError &error ,
                                                                      const MarkovianSequences &iseq ,
                                                                      int identifier) const

{
  bool status = true;
  register int i;
  int offset = I_DEFAULT , nb_value , index;
  double seq_likelihood , max_marginal_entropy , entropy , state_seq_likelihood;
  HiddenVariableOrderMarkov *hmarkov;
  VariableOrderMarkovData *seq;
  ostringstream legend;
  MultiPlotSet *plot_set;


  plot_set = NULL;
  error.init();

  for (i = 0;i < iseq.nb_variable;i++) {
    if ((iseq.type[i] != INT_VALUE) && (iseq.type[i] != REAL_VALUE) && (iseq.type[i] != STATE)) {
      status = false;
      ostringstream error_message , correction_message;
      error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                    << STAT_error[STATR_VARIABLE_TYPE];
      correction_message << STAT_variable_word[INT_VALUE] << " or "
                         << STAT_variable_word[REAL_VALUE];
      error.correction_update((error_message.str()).c_str() , (correction_message.str()).c_str());
    }
  }

  if (nb_output_process == iseq.nb_variable) {
    offset = 0;
  }
  else if ((iseq.type[0] == STATE) && (nb_output_process + 1 == iseq.nb_variable)) {
    offset = 1;
  }
  else {
    status = false;
    error.update(STAT_error[STATR_NB_OUTPUT_PROCESS]);
  }

  if (offset != I_DEFAULT) {
    for (i = 0;i < nb_output_process;i++) {
      if ((categorical_process[i]) || (discrete_parametric_process[i])) {
        if (iseq.type[i + offset] == REAL_VALUE) {
          status = false;
          ostringstream error_message;
          error_message << STAT_label[STATL_VARIABLE] << " " << i + offset + 1 << ": "
                        << STAT_error[STATR_VARIABLE_TYPE];
          error.correction_update((error_message.str()).c_str() , STAT_variable_word[INT_VALUE]);
        }

        else {
          if (iseq.min_value[i + offset] < 0) {
            status = false;
            ostringstream error_message;
            error_message << STAT_label[STATL_VARIABLE] << " " << i + offset + 1 << ": "
                          << STAT_error[STATR_POSITIVE_MIN_VALUE];
            error.update((error_message.str()).c_str());
          }

          if (!(iseq.marginal_distribution[i + offset])) {
            status = false;
            ostringstream error_message;
            error_message << STAT_label[STATL_VARIABLE] << " " << i + offset + 1 << ": "
                          << STAT_error[STATR_MARGINAL_FREQUENCY_DISTRIBUTION];
            error.update((error_message.str()).c_str());
          }

          else {
            if (categorical_process[i]) {
              nb_value = categorical_process[i]->nb_value;
            }
            else {
              nb_value = discrete_parametric_process[i]->nb_value;
            }

            if (nb_value < iseq.marginal_distribution[i + offset]->nb_value) {
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

  for (i = 0;i < iseq.nb_sequence;i++) {
    if (identifier == iseq.identifier[i]) {
      index = i;
      break;
    }
  }

  if (i == iseq.nb_sequence) {
    status = false;
    error.update(SEQ_error[SEQR_SEQUENCE_IDENTIFIER]);
  }

  if (status) {
    plot_set = new MultiPlotSet(4);

    MultiPlotSet &plot = *plot_set;

    plot.border = "15 lw 0";

    if (iseq.type[0] != STATE) {
      seq = new VariableOrderMarkovData(iseq);
    }
    else {
      seq = new VariableOrderMarkovData(iseq , 'c' , (type == 'e' ? true : false));
    }

    seq_likelihood = forward_backward(*seq , index , NULL , plot_set , 'p' ,
                                      max_marginal_entropy , entropy);

    if (seq_likelihood == D_INF) {
      delete plot_set;
      plot_set = NULL;
      error.update(SEQ_error[SEQR_SEQUENCE_INCOMPATIBLE_MODEL]);
    }

    else {
      hmarkov = new HiddenVariableOrderMarkov(*this , false);

      hmarkov->create_cumul();
      hmarkov->log_computation();
      state_seq_likelihood = hmarkov->viterbi_forward_backward(*seq , index , NULL , &plot[0] ,
                                                               'p' , seq_likelihood);

      // 1ere vue : probabilitees maximum

      plot[0].title = SEQ_label[SEQL_MAX_POSTERIOR_STATE_PROBABILITY];

      if (seq->index_parameter) {
        plot[0].xrange = Range(seq->index_parameter[index][0] , seq->index_parameter[index][seq->length[index] - 1]);
        if (seq->index_parameter[index][seq->length[index] - 1] - seq->index_parameter[index][0] < TIC_THRESHOLD) {
          plot[0].xtics = 1;
        }
      }

      else {
        plot[0].xrange = Range(0 , seq->length[index] - 1);
        if (seq->length[index] - 1 < TIC_THRESHOLD) {
          plot[0].xtics = 1;
        }
      }

      plot[0].yrange = Range(0. , exp(state_seq_likelihood - seq_likelihood));

      for (i = 0;i < nb_state;i++) {
        legend.str("");
        legend << STAT_label[STATL_STATE] << " " << i;
        plot[0][i].legend = legend.str();

        plot[0][i].style = "linespoints";
      }

      // 2eme vue : probabilitees lissees

      plot[1].title = SEQ_label[SEQL_POSTERIOR_STATE_PROBABILITY];

      if (seq->index_parameter) {
        plot[1].xrange = Range(seq->index_parameter[index][0] , seq->index_parameter[index][seq->length[index] - 1]);
        if (seq->index_parameter[index][seq->length[index] - 1] - seq->index_parameter[index][0] < TIC_THRESHOLD) {
          plot[1].xtics = 1;
        }
      }

      else {
        plot[1].xrange = Range(0 , seq->length[index] - 1);
        if (seq->length[index] - 1 < TIC_THRESHOLD) {
          plot[1].xtics = 1;
        }
      }

      plot[1].yrange = Range(0. , 1.);

      for (i = 0;i < nb_state;i++) {
        legend.str("");
        legend << STAT_label[STATL_STATE] << " " << i;
        plot[1][i].legend = legend.str();

        plot[1][i].style = "linespoints";
      }

      // 3eme vue : profils d'entropies conditionnelles

      if (seq->index_parameter) {
        plot[2].xrange = Range(seq->index_parameter[index][0] , seq->index_parameter[index][seq->length[index] - 1]);
        if (seq->index_parameter[index][seq->length[index] - 1] - seq->index_parameter[index][0] < TIC_THRESHOLD) {
          plot[2].xtics = 1;
        }
      }

      else {
        plot[2].xrange = Range(0 , seq->length[index] - 1);
        if (seq->length[index] - 1 < TIC_THRESHOLD) {
          plot[2].xtics = 1;
        }
      }

      plot[2].yrange = Range(0. , max_marginal_entropy);

      plot[2][0].legend = SEQ_label[SEQL_CONDITIONAL_ENTROPY];
      plot[2][0].style = "linespoints";

      plot[2][1].legend = SEQ_label[SEQL_CONDITIONAL_ENTROPY];
      plot[2][1].style = "linespoints";

      plot[2][2].legend = SEQ_label[SEQL_MARGINAL_ENTROPY];
      plot[2][2].style = "linespoints";

      // 4eme vue : profils d'entropies partielles

      if (seq->index_parameter) {
        plot[3].xrange = Range(seq->index_parameter[index][0] , seq->index_parameter[index][seq->length[index] - 1]);
        if (seq->index_parameter[index][seq->length[index] - 1] - seq->index_parameter[index][0] < TIC_THRESHOLD) {
          plot[3].xtics = 1;
        }
      }

      else {
        plot[3].xrange = Range(0 , seq->length[index] - 1);
        if (seq->length[index] - 1 < TIC_THRESHOLD) {
          plot[3].xtics = 1;
        }
      }

      plot[3].yrange = Range(0. ,entropy);

      plot[3][0].legend = SEQ_label[SEQL_PARTIAL_STATE_SEQUENCE_ENTROPY];
      plot[3][0].style = "linespoints";

      plot[3][1].legend = SEQ_label[SEQL_PARTIAL_STATE_SEQUENCE_ENTROPY];
      plot[3][1].style = "linespoints";

      delete hmarkov;
    }

    delete seq;
  }

  return plot_set;
}


/*--------------------------------------------------------------*
 *
 *  Calcul des profils d'etats et d'entropie par l'algorithme forward-backward,
 *  des profils d'etats par l'algorithme de Viterbi forward-backward et
 *  sortie graphique des resultats.
 *
 *  arguments : reference sur un objet StatError,
 *              identificateur de la sequence.
 *
 *--------------------------------------------------------------*/

MultiPlotSet* HiddenVariableOrderMarkov::state_profile_plotable_write(StatError &error ,
                                                                      int identifier) const

{
  MultiPlotSet *plot_set;


  error.init();

  if (!markov_data) {
    plot_set = NULL;
    error.update(STAT_error[STATR_NO_DATA]);
  }
  else {
    plot_set = state_profile_plotable_write(error , *markov_data , identifier);
  }

  return plot_set;
}


/*--------------------------------------------------------------*
 *
 *  Calcul des sequences d'etats les plus probables.
 *
 *  arguments : references sur un objet StatError et sur un objet MarkovianSequences,
 *              flag sur le calcul des caracteristiques.
 *
 *--------------------------------------------------------------*/

VariableOrderMarkovData* HiddenVariableOrderMarkov::state_sequence_computation(StatError &error ,
                                                                               const MarkovianSequences &iseq ,
                                                                               bool characteristic_flag) const

{
  bool status = true;
  register int i;
  int nb_value;
  HiddenVariableOrderMarkov *hmarkov;
  VariableOrderMarkovData *seq;


  seq = NULL;
  error.init();

  for (i = 0;i < iseq.nb_variable;i++) {
    if ((iseq.type[i] != INT_VALUE) && (iseq.type[i] != REAL_VALUE)) {
      status = false;
      ostringstream error_message , correction_message;
      error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                    << STAT_error[STATR_VARIABLE_TYPE];
      correction_message << STAT_variable_word[INT_VALUE] << " or "
                         << STAT_variable_word[REAL_VALUE];
      error.correction_update((error_message.str()).c_str() , (correction_message.str()).c_str());
    }
  }

  if (nb_output_process != iseq.nb_variable) {
    status = false;
    error.update(STAT_error[STATR_NB_OUTPUT_PROCESS]);
  }

  else {
    for (i = 0;i < nb_output_process;i++) {
      if ((categorical_process[i]) || (discrete_parametric_process[i])) {
        if (iseq.type[i] == REAL_VALUE) {
          status = false;
          ostringstream error_message;
          error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                        << STAT_error[STATR_VARIABLE_TYPE];
          error.correction_update((error_message.str()).c_str() , STAT_variable_word[INT_VALUE]);
        }

        else {
          if (iseq.min_value[i] < 0) {
            status = false;
            ostringstream error_message;
            error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                          << STAT_error[STATR_POSITIVE_MIN_VALUE];
            error.update((error_message.str()).c_str());
          }

          if (!(iseq.marginal_distribution[i])) {
            status = false;
            ostringstream error_message;
            error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                          << STAT_error[STATR_MARGINAL_FREQUENCY_DISTRIBUTION];
            error.update((error_message.str()).c_str());
          }

          else {
            if (categorical_process[i]) {
              nb_value = categorical_process[i]->nb_value;
            }
            else {
              nb_value = discrete_parametric_process[i]->nb_value;
            }

            if (nb_value < iseq.marginal_distribution[i]->nb_value) {
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

  if (status) {
    seq = new VariableOrderMarkovData(iseq , 'a' , (type == 'e' ? true : false));

    seq->markov = new VariableOrderMarkov(*this , false);

    hmarkov = new HiddenVariableOrderMarkov(*this , false);

    hmarkov->forward_backward(*seq);

    hmarkov->create_cumul();
    hmarkov->log_computation();
    hmarkov->viterbi(*seq);

    // extraction des caracteristiques des sequences et
    // calcul des lois caracteristiques du modele

    if (seq->restoration_likelihood == D_INF) {
      delete seq;
      seq = NULL;
      error.update(SEQ_error[SEQR_STATE_SEQUENCE_COMPUTATION_FAILURE]);
    }

    else {
      seq->likelihood = likelihood_computation(iseq , seq->posterior_probability);

/*      seq->min_value_computation(0);
      seq->max_value_computation(0); */

      seq->min_value[0] = 0;
      seq->max_value[0] = nb_state - 1;
      seq->build_marginal_frequency_distribution(0);
      seq->build_characteristic(0);

      seq->build_transition_count(*hmarkov);
      seq->build_observation_frequency_distribution(nb_state);

/*      if ((seq->max_value[0] < nb_state - 1) || (!(seq->characteristics[0]))) {
        delete seq;
        seq = NULL;
        error.update(SEQ_error[SEQR_STATES_NOT_REPRESENTED]);
      }

      else if (characteristic_flag) { */
      if (characteristic_flag) {
        seq->markov->characteristic_computation(*seq , true);
      }
    }

    delete hmarkov;
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Comparaison de differentes chaines de Markov d'ordre variable cachees
 *  pour un ensemble de sequences.
 *
 *  arguments : reference sur un objet StatError, stream, nombre de chaines
 *              de Markov cachees, pointeur sur les chaines de Markov cachees,
 *              type algorithme (forward ou Viterbi), path.
 *
 *--------------------------------------------------------------*/

bool MarkovianSequences::comparison(StatError &error , ostream &os , int nb_model ,
                                    const HiddenVariableOrderMarkov **ihmarkov ,
                                    int algorithm , const char *path) const

{
  bool status = true;
  register int i , j;
  int nb_value;
  double **likelihood;
  HiddenVariableOrderMarkov **hmarkov;
  VariableOrderMarkovData *seq;


  error.init();

  for (i = 0;i < nb_variable;i++) {
    if ((type[i] != INT_VALUE) && (type[i] != REAL_VALUE) && (type[i] != STATE)) {
      status = false;
      ostringstream error_message , correction_message;
      error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                    << STAT_error[STATR_VARIABLE_TYPE];
      correction_message << STAT_variable_word[INT_VALUE] << " or "
                         << STAT_variable_word[REAL_VALUE];
      error.correction_update((error_message.str()).c_str() , (correction_message.str()).c_str());
    }
  }

  for (i = 0;i < nb_model;i++) {
    if (ihmarkov[i]->nb_output_process != nb_variable) {
      status = false;
      ostringstream error_message;
      error_message << SEQ_label[SEQL_HIDDEN_MARKOV_CHAIN] << " " << i + 1 << ": "
                    << STAT_error[STATR_NB_OUTPUT_PROCESS];
      error.update((error_message.str()).c_str());
    }

    else {
      for (j = 0;j < nb_variable;j++) {
        if ((ihmarkov[i]->categorical_process[j]) || (ihmarkov[i]->discrete_parametric_process[j])) {
          if (type[j] == REAL_VALUE) {
            status = false;
            ostringstream error_message;
            error_message << STAT_label[STATL_VARIABLE] << " " << j + 1 << ": "
                          << STAT_error[STATR_VARIABLE_TYPE];
            error.correction_update((error_message.str()).c_str() , STAT_variable_word[INT_VALUE]);
          }

          else {
            if (min_value[j] < 0) {
              status = false;
              ostringstream error_message;
              error_message << STAT_label[STATL_VARIABLE] << " " << j + 1 << ": "
                            << STAT_error[STATR_POSITIVE_MIN_VALUE];
              error.update((error_message.str()).c_str());
            }

            if (!marginal_distribution[j]) {
              status = false;
              ostringstream error_message;
              error_message << STAT_label[STATL_VARIABLE] << " " << j + 1 << ": "
                            << STAT_error[STATR_MARGINAL_FREQUENCY_DISTRIBUTION];
              error.update((error_message.str()).c_str());
            }

            else {
              if (ihmarkov[i]->categorical_process[j]) {
                nb_value = ihmarkov[i]->categorical_process[j]->nb_value;
              }
              else {
                nb_value = ihmarkov[i]->discrete_parametric_process[j]->nb_value;
              }

              if (nb_value < marginal_distribution[j]->nb_value) {
                status = false;
                ostringstream error_message;
                error_message << SEQ_label[SEQL_HIDDEN_MARKOV_CHAIN] << " " << i + 1 << ": "
                              << STAT_label[STATL_OUTPUT_PROCESS] << " " << j + 1 << ": "
                              << STAT_error[STATR_NB_OUTPUT];
                error.update((error_message.str()).c_str());
              }
            }
          }
        }
      }
    }
  }

  if (status) {
    likelihood = new double*[nb_sequence];
    for (i = 0;i < nb_sequence;i++) {
      likelihood[i] = new double[nb_model];
    }

    if (algorithm == VITERBI) {
      hmarkov = new HiddenVariableOrderMarkov*[nb_model];
      for (i = 0;i < nb_model;i++) {
        hmarkov[i] = new HiddenVariableOrderMarkov(*(ihmarkov[i]) , false);
        hmarkov[i]->create_cumul();
        hmarkov[i]->log_computation();
      }

      seq = new VariableOrderMarkovData(*this);
    }

    // pour chaque sequence, calcul de la vraisemblance (FORWARD) ou de la vraisemblance
    // de la sequence d'etats la plus probable (VITERBI) pour chaque modele possible

    for (i = 0;i < nb_sequence;i++) {
      for (j = 0;j < nb_model;j++) {
        switch (algorithm) {
        case FORWARD :
          likelihood[i][j] = ihmarkov[j]->likelihood_computation(*this , NULL , i);
          break;
        case VITERBI :
          likelihood[i][j] = hmarkov[j]->viterbi(*seq , NULL , i);
          break;
        }
      }
    }

#   ifdef MESSAGE
    likelihood_write(os , nb_model , likelihood , SEQ_label[SEQL_HIDDEN_MARKOV_CHAIN] ,
                     true , algorithm);
#   endif

    if (path) {
      status = likelihood_write(error , path , nb_model , likelihood ,
                                SEQ_label[SEQL_HIDDEN_MARKOV_CHAIN] , algorithm);
    }

    for (i = 0;i < nb_sequence;i++) {
      delete [] likelihood[i];
    }
    delete [] likelihood;

    if (algorithm == VITERBI) {
      for (i = 0;i < nb_model;i++) {
        delete hmarkov[i];
      }
      delete [] hmarkov;

      delete seq;
    }
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Simulation par une chaine de Markov d'ordre variable cachee.
 *
 *  arguments : reference sur un objet StatError,
 *              loi empirique des longueurs des sequences,
 *              flag sur le calcul des lois de comptage,
 *              flag calcul d'une divergence de Kullback-Leibler.
 *
 *--------------------------------------------------------------*/

VariableOrderMarkovData* HiddenVariableOrderMarkov::simulation(StatError &error ,
                                                               const FrequencyDistribution &length_distribution ,
                                                               bool counting_flag ,
                                                               bool divergence_flag) const

{
  register int i;
  MarkovianSequences *observed_seq;
  VariableOrderMarkovData *seq;


  seq = VariableOrderMarkov::simulation(error , length_distribution , counting_flag , divergence_flag);

  if ((seq) && (!divergence_flag)) {
    seq->posterior_probability = new double[seq->nb_sequence];
    for (i = 0;i < seq->nb_sequence;i++) {
      seq->posterior_probability[i] = VariableOrderMarkov::likelihood_computation(*seq , i);
    }

/*    for (i = 0;i < nb_output_process;i++) {
      if (continuous_parametric_process[i]) {
        seq->restoration_likelihood = VariableOrderMarkov::likelihood_computation(*seq , I_DEFAULT);
        break;
      }
    } */

    observed_seq = seq->remove_variable_1();
    seq->likelihood = likelihood_computation(*observed_seq , seq->posterior_probability);
    delete observed_seq;

    forward_backward(*seq);
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Simulation par une chaine de Markov d'ordre variable cachee.
 *
 *  arguments : reference sur un objet StatError,
 *              nombre et longueur des sequences.
 *
 *--------------------------------------------------------------*/

VariableOrderMarkovData* HiddenVariableOrderMarkov::simulation(StatError &error ,
                                                               int nb_sequence , int length ,
                                                               bool counting_flag) const

{
  register int i;
  MarkovianSequences *observed_seq;
  VariableOrderMarkovData *seq;


  seq = VariableOrderMarkov::simulation(error , nb_sequence , length , counting_flag);

  if (seq) {
    seq->posterior_probability = new double[seq->nb_sequence];
    for (i = 0;i < seq->nb_sequence;i++) {
      seq->posterior_probability[i] = VariableOrderMarkov::likelihood_computation(*seq , i);
    }

/*    for (i = 0;i < nb_output_process;i++) {
      if (continuous_parametric_process[i]) {
        seq->restoration_likelihood = VariableOrderMarkov::likelihood_computation(*seq , I_DEFAULT);
        break;
      }
    } */

    observed_seq = seq->remove_variable_1();
    seq->likelihood = likelihood_computation(*observed_seq , seq->posterior_probability);
    delete observed_seq;

    forward_backward(*seq);
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Simulation par une chaine de Markov d'ordre variable cachee.
 *
 *  arguments : reference sur un objet StatError, nombre de sequences,
 *              reference sur un objet MarkovianSequences.
 *
 *--------------------------------------------------------------*/

VariableOrderMarkovData* HiddenVariableOrderMarkov::simulation(StatError &error ,
                                                               int nb_sequence ,
                                                               const MarkovianSequences &iseq ,
                                                               bool counting_flag) const

{
  register int i;
  MarkovianSequences *observed_seq;
  VariableOrderMarkovData *seq;


  seq = VariableOrderMarkov::simulation(error , nb_sequence , iseq , counting_flag);

  if (seq) {
    seq->posterior_probability = new double[seq->nb_sequence];
    for (i = 0;i < seq->nb_sequence;i++) {
      seq->posterior_probability[i] = VariableOrderMarkov::likelihood_computation(*seq , i);
    }

/*    for (i = 0;i < nb_output_process;i++) {
      if (continuous_parametric_process[i]) {
        seq->restoration_likelihood = VariableOrderMarkov::likelihood_computation(*seq , I_DEFAULT);
        break;
      }
    } */

    observed_seq = seq->remove_variable_1();
    seq->likelihood = likelihood_computation(*observed_seq , seq->posterior_probability);
    delete observed_seq;

    forward_backward(*seq);
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Comparaison de chaines de Markov d'ordre variable cachees par calcul
 *  de divergences de Kullback-Leibler.
 *
 *  arguments : reference sur un objet StatError, stream, nombre de chaines
 *              de Markov cachees, pointeur sur les chaines de Markov cachees,
 *              loi empirique des longueurs des sequences, path.
 *
 *--------------------------------------------------------------*/

DistanceMatrix* HiddenVariableOrderMarkov::divergence_computation(StatError &error , ostream &os , int nb_model ,
                                                                  const HiddenVariableOrderMarkov **ihmarkov ,
                                                                  FrequencyDistribution **length_distribution ,
                                                                  const char *path) const

{
  bool status = true , lstatus;
  register int i , j , k;
  int cumul_length , nb_failure;
  double **likelihood;
  long double divergence;
  const HiddenVariableOrderMarkov **hmarkov;
  MarkovianSequences *seq;
  VariableOrderMarkovData *simul_seq;
  DistanceMatrix *dist_matrix;
  ofstream *out_file;


  dist_matrix = NULL;
  error.init();

  for (i = 0;i < nb_model - 1;i++) {
    if (ihmarkov[i]->type != type) {
      status = false;
      ostringstream error_message;
      error_message << SEQ_label[SEQL_HIDDEN_MARKOV_CHAIN] << " " << i + 2 << ": "
                    << SEQ_error[SEQR_MODEL_TYPE];
      error.update((error_message.str()).c_str());
    }

    if (ihmarkov[i]->nb_output_process != nb_output_process) {
      status = false;
      ostringstream error_message;
      error_message << SEQ_label[SEQL_HIDDEN_MARKOV_CHAIN] << " " << i + 2 << ": "
                    << STAT_error[STATR_NB_OUTPUT_PROCESS];
      error.update((error_message.str()).c_str());
    }

    else {
      for (j = 0;j < nb_output_process;j++) {
        if ((categorical_process[j]) && (ihmarkov[i]->categorical_process[j]) &&
            (ihmarkov[i]->categorical_process[j]->nb_value != categorical_process[j]->nb_value)) {
          status = false;
          ostringstream error_message;
          error_message << SEQ_label[SEQL_HIDDEN_MARKOV_CHAIN] << " " << i + 2 << ": "
                        << STAT_label[STATL_OUTPUT_PROCESS] << " " << j << " "
                        << STAT_error[STATR_NB_OUTPUT];
          error.update((error_message.str()).c_str());
        }

        if (((continuous_parametric_process[j]) && (!(ihmarkov[i]->continuous_parametric_process[j]))) ||
            ((!continuous_parametric_process[j]) && (ihmarkov[i]->continuous_parametric_process[j]))) {
          status = false;
          ostringstream error_message;
          error_message << SEQ_label[SEQL_HIDDEN_MARKOV_CHAIN] << " " << i + 2 << ": "
                        << STAT_label[STATL_OUTPUT_PROCESS] << " " << j << " "
                        << SEQ_error[SEQR_OUTPUT_PROCESS_TYPE];
          error.update((error_message.str()).c_str());
        }
      }
    }
  }

  for (i = 0;i < nb_model;i++) {
    lstatus = true;

    if ((length_distribution[i]->nb_element < 1) || (length_distribution[i]->nb_element > NB_SEQUENCE)) {
      lstatus = false;
      ostringstream error_message;
      error_message << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " "
                    << i + 1 << ": "  << SEQ_error[SEQR_NB_SEQUENCE];
      error.update((error_message.str()).c_str());
    }
    if (length_distribution[i]->offset < 2) {
      lstatus = false;
      ostringstream error_message;
      error_message << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " "
                    << i + 1 << ": "  << SEQ_error[SEQR_SHORT_SEQUENCE_LENGTH];
      error.update((error_message.str()).c_str());
    }
    if (length_distribution[i]->nb_value - 1 > MAX_LENGTH) {
      lstatus = false;
      ostringstream error_message;
      error_message << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " "
                    << i + 1 << ": "  << SEQ_error[SEQR_LONG_SEQUENCE_LENGTH];
      error.update((error_message.str()).c_str());
    }

    if (!lstatus) {
      status = false;
    }

    else {
      cumul_length = 0;
      for (j = length_distribution[i]->offset;j < length_distribution[i]->nb_value;j++) {
        cumul_length += j * length_distribution[i]->frequency[j];
      }

      if (cumul_length > CUMUL_LENGTH) {
        status = false;
        ostringstream error_message;
        error_message << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " "
                      << i + 1 << ": "  << SEQ_error[SEQR_CUMUL_SEQUENCE_LENGTH];
        error.update((error_message.str()).c_str());
      }
    }
  }

  if (status) {
    out_file = NULL;

    if (path) {
      out_file = new ofstream(path);

      if (!out_file) {
        error.update(STAT_error[STATR_FILE_NAME]);

#       ifdef MESSAGE
        os << error;
#       endif

      }
    }

    hmarkov = new const HiddenVariableOrderMarkov*[nb_model];

    hmarkov[0] = this;
    for (i = 1;i < nb_model;i++) {
      hmarkov[i] = ihmarkov[i - 1];
    }

    dist_matrix = new DistanceMatrix(nb_model , SEQ_label[SEQL_HIDDEN_MARKOV_CHAIN]);

    for (i = 0;i < nb_model;i++) {

      // simulation d'un echantillon de sequences a partir d'une chaine de Markov cachee

      simul_seq = hmarkov[i]->simulation(error , *length_distribution[i] , false , true);
      seq = simul_seq->remove_variable_1();

      likelihood = new double*[seq->nb_sequence];
      for (j = 0;j < seq->nb_sequence;j++) {
        likelihood[j] = new double[nb_model];
      }

      for (j = 0;j < seq->nb_sequence;j++) {
        likelihood[j][i] = hmarkov[i]->likelihood_computation(*seq , NULL , j);

#       ifdef MESSAGE
        if (likelihood[j][i] == D_INF) {
          os << "\nERROR - " << SEQ_error[SEQR_REFERENCE_MODEL] << ": " << i + 1 << endl;
        }
#       endif

      }

      // calcul des vraisemblances de l'echantillon pour chacune des chaines de Markov cachees

      for (j = 0;j < nb_model;j++) {
        if (j != i) {
          divergence = 0.;
          cumul_length = 0;
          nb_failure = 0;

          for (k = 0;k < seq->nb_sequence;k++) {
            likelihood[k][j] = hmarkov[j]->likelihood_computation(*seq , NULL , k);

//            if (divergence != -D_INF) {
              if (likelihood[k][j] != D_INF) {
                divergence += likelihood[k][i] - likelihood[k][j];
                cumul_length += seq->length[k];
              }
              else {
                nb_failure++;
//                divergence = -D_INF;
              }
//            }
          }

#         ifdef MESSAGE
          if (nb_failure > 0) {
            os << "\nWARNING - " << SEQ_error[SEQR_REFERENCE_MODEL] << ": " << i + 1 << ", "
               << SEQ_error[SEQR_TARGET_MODEL] << ": " << j + 1 << " - "
               << SEQ_error[SEQR_DIVERGENCE_NB_FAILURE] << ": " << nb_failure << endl;
          }
#         endif

//          if (divergence != -D_INF) {
            dist_matrix->update(i + 1 , j + 1 , divergence , cumul_length);
//          }
        }
      }

#     ifdef MESSAGE
      os << SEQ_label[SEQL_HIDDEN_MARKOV_CHAIN] << " " << i + 1 << ": " << seq->nb_sequence << " "
         << SEQ_label[SEQL_SIMULATED] << " " << SEQ_label[seq->nb_sequence == 1 ? SEQL_SEQUENCE : SEQL_SEQUENCES] << endl;
      seq->likelihood_write(os , nb_model , likelihood , SEQ_label[SEQL_HIDDEN_MARKOV_CHAIN]);
#     endif

      if (out_file) {
        *out_file << SEQ_label[SEQL_HIDDEN_MARKOV_CHAIN] << " " << i + 1 << ": " << seq->nb_sequence << " "
                  << SEQ_label[SEQL_SIMULATED] << " " << SEQ_label[seq->nb_sequence == 1 ? SEQL_SEQUENCE : SEQL_SEQUENCES] << endl;
        seq->likelihood_write(*out_file , nb_model , likelihood , SEQ_label[SEQL_HIDDEN_MARKOV_CHAIN]);
      }

      for (j = 0;j < seq->nb_sequence;j++) {
        delete [] likelihood[j];
      }
      delete [] likelihood;

      delete seq;
      delete simul_seq;
    }

    if (out_file) {
      out_file->close();
      delete out_file;
    }

    delete hmarkov;
  }

  return dist_matrix;
}


/*--------------------------------------------------------------*
 *
 *  Comparaison de chaines de Markov d'ordre variable cachees par calcul
 *  de divergences de Kullback-Leibler.
 *
 *  arguments : reference sur un objet StatError, stream, nombre de chaines
 *              de Markov cachees, pointeur sur les chaines de Markov cachees,
 *              nombre et longueur des sequences, path.
 *
 *--------------------------------------------------------------*/

DistanceMatrix* HiddenVariableOrderMarkov::divergence_computation(StatError &error , ostream &os ,
                                                                  int nb_model , const HiddenVariableOrderMarkov **hmarkov ,
                                                                  int nb_sequence , int length , const char *path) const

{
  bool status = true;
  register int i;
  FrequencyDistribution **length_distribution;
  DistanceMatrix *dist_matrix;


  dist_matrix = NULL;
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
    length_distribution = new FrequencyDistribution*[nb_model];

    length_distribution[0] = new FrequencyDistribution(length + 1);

    length_distribution[0]->nb_element = nb_sequence;
    length_distribution[0]->offset = length;
    length_distribution[0]->max = nb_sequence;
    length_distribution[0]->mean = length;
    length_distribution[0]->variance = 0.;
    length_distribution[0]->frequency[length] = nb_sequence;

    for (i = 1;i < nb_model;i++) {
      length_distribution[i] = new FrequencyDistribution(*length_distribution[0]);
    }

    dist_matrix = divergence_computation(error , os , nb_model , hmarkov , length_distribution , path);

    for (i = 0;i < nb_model;i++) {
      delete length_distribution[i];
    }
    delete [] length_distribution;
  }

  return dist_matrix;
}


/*--------------------------------------------------------------*
 *
 *  Comparaison de chaines de Markov d'ordre variable cachees par calcul
 *  de divergences de Kullback-Leibler.
 *
 *  arguments : reference sur un objet StatError, stream, nombre de chaines
 *              de Markov cachees, pointeur sur les chaines de Markov cachees,
 *              pointeurs sur des objets MarkovianSequences, path.
 *
 *--------------------------------------------------------------*/

DistanceMatrix* HiddenVariableOrderMarkov::divergence_computation(StatError &error , ostream &os ,
                                                                  int nb_model , const HiddenVariableOrderMarkov **hmarkov ,
                                                                  int nb_sequence , const MarkovianSequences **seq ,
                                                                  const char *path) const

{
  register int i;
  FrequencyDistribution **length_distribution;
  DistanceMatrix *dist_matrix;


  error.init();

  if ((nb_sequence < 1) || (nb_sequence > NB_SEQUENCE)) {
    dist_matrix = NULL;
    error.update(SEQ_error[SEQR_NB_SEQUENCE]);
  }

  else {
    length_distribution = new FrequencyDistribution*[nb_model];
    for (i = 0;i < nb_model;i++) {
      length_distribution[i] = seq[i]->length_distribution->frequency_scale(nb_sequence);
    }

    dist_matrix = divergence_computation(error , os , nb_model , hmarkov , length_distribution , path);

    for (i = 0;i < nb_model;i++) {
      delete length_distribution[i];
    }
    delete [] length_distribution;
  }

  return dist_matrix;
}


};  // namespace sequence_analysis
