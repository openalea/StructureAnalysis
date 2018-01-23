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
 *       $Id: change_points2.cpp 18669 2015-11-09 12:08:08Z guedon $
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

#include "sequences.h"
#include "sequence_label.h"

using namespace std;
using namespace stat_tool;


namespace sequence_analysis {


extern double log_factorial(int value);
extern double log_binomial_coefficient(int inf_bound , double parameter , int value);



/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the contrast functions within a forward recursion.
 *
 *  \param[in] time                time instant,
 *  \param[in] index               sequence index,
 *  \param[in] model_type          segment model types,
 *  \param[in] common_contrast     flag contrast functions common to the individuals,
 *  \param[in] factorial           log factorial for Poisson models,
 *  \param[in] shape_parameter     negative binomial shape parameters,
 *  \param[in] binomial_coeff      log binomial coefficients for negative binomial models,
 *  \param[in] seq_mean            sequence means for Gaussian change in the variance models or
 *                                 stationary piecewise autoregressive models,
 *  \param[in] seq_index_parameter index parameters,
 *  \param[in] hyperparam          hyperparameters for Bayesian models,
 *  \param[in] rank                ranks for ordinal variables,
 *  \param[in] contrast            contrast functions,
 *  \param[in] nb_segment          number of segments for bounding time loops.
 */
/*--------------------------------------------------------------*/

void Sequences::forward_contrast(int time , int index , segment_model *model_type , bool common_contrast ,
                                 double ***factorial , double *shape_parameter , double ***binomial_coeff ,
                                 double **seq_mean , int *seq_index_parameter , double **hyperparam ,
                                 double **rank , long double *contrast , int nb_segment) const

{
  int i , j , k , m;
  int max_nb_value , count , *frequency , *inf_bound_parameter;
  double sum , factorial_sum , proba , binomial_coeff_sum , diff , index_parameter_sum ,
         index_parameter_diff , shifted_diff , range_diff , mean , buff;
  long double index_parameter_square_sum , square_sum , mix_square_sum , shifted_square_sum ,
              autocovariance , prior_contrast , square_sum_term[3] , **residual;


  // initializations

  max_nb_value = 0;
  inf_bound_parameter = new int[nb_variable];
  residual = NULL;

  for (i = 1;i < nb_variable;i++) {
    if ((model_type[i - 1] == CATEGORICAL_CHANGE) && (marginal_distribution[i]->nb_value > max_nb_value)) {
      max_nb_value = marginal_distribution[i]->nb_value;
    }

    if ((model_type[i - 1] == NEGATIVE_BINOMIAL_0_CHANGE) || (model_type[i - 1] == NEGATIVE_BINOMIAL_1_CHANGE)) {
      switch (model_type[i - 1]) {
      case NEGATIVE_BINOMIAL_0_CHANGE :
        inf_bound_parameter[i - 1] = 0;
        break;
      case NEGATIVE_BINOMIAL_1_CHANGE :
        inf_bound_parameter[i - 1] = 1;
        break;
      }
    }

    if (((i == 1) && ((model_type[0] == MEAN_CHANGE) || (model_type[0] == INTERCEPT_SLOPE_CHANGE))) ||
        (((model_type[i - 1] == GAUSSIAN_CHANGE) || (model_type[i - 1] == VARIANCE_CHANGE) ||
          (model_type[i - 1] == ORDINAL_GAUSSIAN_CHANGE) || (model_type[i - 1] == LINEAR_MODEL_CHANGE) ||
          (model_type[i - 1] == AUTOREGRESSIVE_MODEL_CHANGE) ||
          (model_type[i - 1] == STATIONARY_AUTOREGRESSIVE_MODEL_CHANGE)) && (!residual))) {
      residual = new long double*[MAX(nb_sequence , 2)];
      if ((index != I_DEFAULT) || (!common_contrast)) {
        for (j = 0;j < nb_sequence;j++) {
          if ((index == I_DEFAULT) || (index == j)) {
            residual[j] = new long double[time + 1];
          }
          else {
            residual[j] = NULL;
          }
        }
      }
      else {
        residual[0] = new long double[time + 1];
      }
    }
  }

  if (max_nb_value > 0) {
    frequency = new int[max_nb_value];
  }
  else {
    frequency = NULL;
  }

  // computation of segment contrast functions (log-likelihoods or sum of squared deviations)

  for (i = nb_segment;i <= time;i++) {
    contrast[i] = 0.;
  }

  for (i = 1;i < nb_variable;i++) {
    if (model_type[i - 1] == CATEGORICAL_CHANGE) {
      if ((index != I_DEFAULT) || (!common_contrast)) {
        for (j = 0;j < nb_sequence;j++) {
          if ((index == I_DEFAULT) || (index == j)) {
            for (k = 0;k < marginal_distribution[i]->nb_value;k++) {
              frequency[k] = 0;
            }

            frequency[int_sequence[j][i][time]]++;

#           ifdef MESSAGE
            sum = 0.;
#           endif

            for (k = time - 1;k >= nb_segment;k--) {

#             ifdef MESSAGE
              sum += (time - k) * log((double)(time - k) / (double)(time - k + 1)) +
                     log((double)(frequency[int_sequence[j][i][k]] + 1) / (double)(time - k + 1));
              if (frequency[int_sequence[j][i][k]] > 0) {
                sum -= frequency[int_sequence[j][i][k]] *
                       log((double)frequency[int_sequence[j][i][k]] / (double)(frequency[int_sequence[j][i][k]] + 1));
              }
/*              frequency[int_sequence[j][i][k]]++;

              if (contrast[k] != D_INF) {
                contrast[k] += sum;
              } */
#             endif

              frequency[int_sequence[j][i][k]]++;
              if (contrast[k] != D_INF) {
                buff = 0.;
                for (m = 0;m < marginal_distribution[i]->nb_value;m++) {
                  if (frequency[m] > 0) {
//                    contrast[k] += frequency[m] * log((double)frequency[m] / (double)(time - k + 1));
                    buff += frequency[m] * log((double)frequency[m] / (double)(time - k + 1));
                  }
                }
                contrast[k] += buff;

#               ifdef MESSAGE
                if ((buff < sum - DOUBLE_ERROR) || (buff > sum + DOUBLE_ERROR)) {
                  cout << "\nERROR: " << k << " " << time << " " << j << " | " << buff << " " << sum << endl;
                }
#               endif

              }
            }
          }
        }
      }

      else {
        for (j = 0;j < marginal_distribution[i]->nb_value;j++) {
          frequency[j] = 0;
        }

        for (j = time;j >= nb_segment;j--) {
          for (k = 0;k < nb_sequence;k++) {
            frequency[int_sequence[k][i][j]]++;
          }

          if (contrast[j] != D_INF) {
            for (k = 0;k < marginal_distribution[i]->nb_value;k++) {
              if (frequency[k] > 0) {
                contrast[j] += frequency[k] * log((double)frequency[k] / (double)(nb_sequence * (time - j + 1)));
              }
            }
          }
        }
      }
    }

    else if (model_type[i - 1] == POISSON_CHANGE) {
/*      for (j = 0;j < nb_sequence;j++) {
        if ((index == I_DEFAULT) || (index == j)) {
          factorial[i][j][time] = log_factorial(int_sequence[j][i][time]);
        }
      } */

      if ((index != I_DEFAULT) || (!common_contrast)) {
        for (j = 0;j < nb_sequence;j++) {
          if ((index == I_DEFAULT) || (index == j)) {
            sum = 0.;
            factorial_sum = 0.;

            for (k = time;k >= nb_segment;k--) {
              sum += int_sequence[j][i][k];
              factorial_sum += factorial[i][j][k];
              if ((contrast[k] != D_INF) && (sum > 0.)) {
                contrast[k] += sum * (log(sum / (time - k + 1)) - 1) - factorial_sum;
              }
            }
          }
        }
      }

      else {
        sum = 0.;
        factorial_sum = 0.;

        for (j = time;j >= nb_segment;j--) {
          for (k = 0;k < nb_sequence;k++) {
            sum += int_sequence[k][i][j];
            factorial_sum += factorial[i][k][j];
          }

          if ((contrast[j] != D_INF) && (sum > 0.)) {
            contrast[j] += sum * (log(sum / (nb_sequence * (time - j + 1))) - 1) - factorial_sum;
          }
        }
      }
    }

    else if ((model_type[i - 1] == NEGATIVE_BINOMIAL_0_CHANGE) || (model_type[i - 1] == NEGATIVE_BINOMIAL_1_CHANGE)) {
/*      for (j = 0;j < nb_sequence;j++) {
        if ((index == I_DEFAULT) || (index == j)) {
          binomial_coeff[i][j][time] = log_binomial_coefficient(inf_bound_parameter[i - 1] , shape_parameter[i - 1] ,
                                                                int_sequence[j][i][time]);
        }
      } */

      if ((index != I_DEFAULT) || (!common_contrast)) {
        for (j = 0;j < nb_sequence;j++) {
          if ((index == I_DEFAULT) || (index == j)) {
            sum = 0.;
            binomial_coeff_sum = 0.;

            for (k = time;k >= nb_segment;k--) {
              sum += int_sequence[j][i][k];
              binomial_coeff_sum += binomial_coeff[i][j][k];

              if (contrast[k] != D_INF) {
                if (sum > inf_bound_parameter[i - 1] * (time - k + 1)) {
                  proba = shape_parameter[i - 1] * (time - k + 1) /
                          ((shape_parameter[i - 1] - inf_bound_parameter[i - 1]) * (time - k + 1) + sum);
                  contrast[k] += binomial_coeff_sum + shape_parameter[i - 1] * (time - k + 1) * log(proba) +
                                 (sum - inf_bound_parameter[i - 1] * (time - k + 1)) * log(1. - proba);
                }
                else {
                  contrast[k] = D_INF;
                }
              }
            }
          }
        }
      }

      else {
        sum = 0.;
        binomial_coeff_sum = 0.;

        for (j = time;j >= nb_segment;j--) {
          for (k = 0;k < nb_sequence;k++) {
            sum += int_sequence[k][i][j];
            binomial_coeff_sum += binomial_coeff[i][k][j];
          }

          if (contrast[j] != D_INF) {
            if (sum > inf_bound_parameter[i - 1] * nb_sequence * (time - j + 1)) {
              proba = shape_parameter[i - 1] * nb_sequence * (time - j + 1) /
                      ((shape_parameter[i - 1] - inf_bound_parameter[i - 1]) * nb_sequence * (time - j + 1) + sum);
              contrast[j] += binomial_coeff_sum + shape_parameter[i - 1] * nb_sequence * (time - j + 1) * log(proba) +
                             (sum - inf_bound_parameter[i - 1] * nb_sequence * (time - j + 1)) * log(1. - proba);
            }
            else {
              contrast[j] = D_INF;
            }
          }
        }
      }
    }

    else if ((model_type[i - 1] == GAUSSIAN_CHANGE) || (model_type[0] == MEAN_CHANGE)) {
      if ((index != I_DEFAULT) || (!common_contrast)) {
        if (type[i] != REAL_VALUE) {
          for (j = 0;j < nb_sequence;j++) {
            if ((index == I_DEFAULT) || (index == j)) {
              square_sum = 0.;
              sum = int_sequence[j][i][time];
              residual[j][time] = 0.;

              for (k = time - 1;k >= nb_segment;k--) {
                diff = int_sequence[j][i][k] - sum / (time - k);
                square_sum += ((double)(time - k) / (double)(time - k + 1)) * diff * diff;
                sum += int_sequence[j][i][k];
                residual[j][k] = square_sum;
              }
            }
          }
        }

        else {
          for (j = 0;j < nb_sequence;j++) {
            if ((index == I_DEFAULT) || (index == j)) {
              square_sum = 0.;
              sum = real_sequence[j][i][time];
              residual[j][time] = 0.;

              for (k = time - 1;k >= nb_segment;k--) {
                diff = real_sequence[j][i][k] - sum / (time - k);
                square_sum += ((double)(time - k) / (double)(time - k + 1)) * diff * diff;
                sum += real_sequence[j][i][k];
                residual[j][k] = square_sum;
              }
            }
          }
        }

#       ifdef DEBUG
        for (j = 0;j < nb_sequence;j++) {
          if ((index == I_DEFAULT) || (index == j)) {
            cout << time << " | ";
            for (k = time;k >= nb_segment;k--) {
              cout << residual[j][k] << " ";
            }
            cout << endl;
          }
        }
#       endif

      }

      else {
        square_sum = 0.;
        sum = 0.;
        count = 0;

        if (type[i] != REAL_VALUE) {
          for (j = time;j >= nb_segment;j--) {
            for (k = 0;k < nb_sequence;k++) {
              if (count > 0) {
                diff = int_sequence[k][i][j] - sum / count;
                square_sum += ((double)count / (double)(count + 1)) * diff * diff;
              }
              count++;
              sum += int_sequence[k][i][j];
            }
            residual[0][j] = square_sum;
          }
        }

        else {
          for (j = time;j >= nb_segment;j--) {
            for (k = 0;k < nb_sequence;k++) {
              if (count > 0) {
                diff = real_sequence[k][i][j] - sum / count;
                square_sum += ((double)count / (double)(count + 1)) * diff * diff;
              }
              count++;
              sum += real_sequence[k][i][j];
            }
            residual[0][j] = square_sum;
          }
        }

#       ifdef MESSAGE

        // alternative implementation

        square_sum = 0.;
        sum = 0.;

        residual[1] = new long double[length[0]];

        if (type[i] != REAL_VALUE) {
          for (j = time;j >= nb_segment;j--) {
            for (k = 0;k < nb_sequence;k++) {
              sum += int_sequence[k][i][j];
              square_sum += int_sequence[k][i][j] * int_sequence[k][i][j];
            }
            residual[1][j] = square_sum - sum * sum / (nb_sequence * (time - j + 1));
          }
        }

        else {
          for (j = time;j >= nb_segment;j--) {
            for (k = 0;k < nb_sequence;k++) {
              sum += real_sequence[k][i][j];
              square_sum += real_sequence[k][i][j] * real_sequence[k][i][j];
            }
            residual[1][j] = square_sum - sum * sum / (nb_sequence * (time - j + 1));
          }
        }

        for (j = time;j >= nb_segment;j--) {
          if ((residual[1][j] < residual[0][j] - DOUBLE_ERROR) || (residual[1][j] > residual[0][j] + DOUBLE_ERROR)) {
            cout << "\nERROR: " << time << " " << j << " | " << residual[1][j] << " " << residual[0][j] << endl;
          }
        }

        delete [] residual[1];

#       endif
      }
    }

    else if (model_type[i - 1] == VARIANCE_CHANGE) {
      if ((index != I_DEFAULT) || (!common_contrast)) {
        if (type[i] != REAL_VALUE) {
          for (j = 0;j < nb_sequence;j++) {
            if ((index == I_DEFAULT) || (index == j)) {
              square_sum = 0.;
              for (k = time;k >= nb_segment;k--) {
                diff = int_sequence[j][i][k] - seq_mean[i][j];
                square_sum += diff * diff;
                residual[j][k] = square_sum;
              }
            }
          }
        }

        else {
          for (j = 0;j < nb_sequence;j++) {
            if ((index == I_DEFAULT) || (index == j)) {
              square_sum = 0.;
              for (k = time;k >= nb_segment;k--) {
                diff = real_sequence[j][i][k] - seq_mean[i][j];
                square_sum += diff * diff;
                residual[j][k] = square_sum;
              }
            }
          }
        }
      }

      else {
        square_sum = 0.;

        if (type[i] != REAL_VALUE) {
          for (j = time;j >= nb_segment;j--) {
            for (k = 0;k < nb_sequence;k++) {
              diff = int_sequence[k][i][j] - seq_mean[i][0];
              square_sum += diff * diff;
            }
            residual[0][j] = square_sum;
          }
        }

        else {
          for (j = time;j >= nb_segment;j--) {
            for (k = 0;k < nb_sequence;k++) {
              diff = real_sequence[k][i][j] - seq_mean[i][0];
              square_sum += diff * diff;
            }
            residual[0][j] = square_sum;
          }
        }
      }
    }

    else if (model_type[i - 1] == ORDINAL_GAUSSIAN_CHANGE) {
      if ((index != I_DEFAULT) || (!common_contrast)) {
        for (j = 0;j < nb_sequence;j++) {
          if ((index == I_DEFAULT) || (index == j)) {
            square_sum = 0.;
            sum = rank[i][int_sequence[j][i][time]];
            residual[j][time] = 0.;

            for (k = time - 1;k >= nb_segment;k--) {
              diff = rank[i][int_sequence[j][i][k]] - sum / (time - k);
              square_sum += ((double)(time - k) / (double)(time - k + 1)) * diff * diff;
              sum += rank[i][int_sequence[j][i][k]];
              residual[j][k] = square_sum;

              if (residual[j][k] == 0.) {
                residual[j][k] = (time - k + 1) * MIN_RANK_SQUARE_SUM;
              }
            }
          }
        }
      }

      else {
        square_sum = 0.;
        sum = 0.;
        count = 0;

        for (j = time;j >= nb_segment;j--) {
          for (k = 0;k < nb_sequence;k++) {
            if (count > 0) {
              diff = rank[i][int_sequence[k][i][j]] - sum / count;
              square_sum += ((double)count / (double)(count + 1)) * diff * diff;
            }
            count++;
            sum += rank[i][int_sequence[k][i][j]];
          }
          residual[0][j] = square_sum;

          if (residual[0][j] == 0.) {
            residual[0][j] = count * MIN_RANK_SQUARE_SUM;
          }
        }
      }
    }

    else if ((model_type[i - 1] == LINEAR_MODEL_CHANGE) || (model_type[0] == INTERCEPT_SLOPE_CHANGE)) {
      if ((index != I_DEFAULT) || (!common_contrast)) {
        if (type[i] != REAL_VALUE) {
          for (j = 0;j < nb_sequence;j++) {
            if ((index == I_DEFAULT) || (index == j)) {
              index_parameter_square_sum = 0.;
              square_sum = 0.;
              mix_square_sum = 0.;
              index_parameter_sum = seq_index_parameter[time];
              sum = int_sequence[j][i][time];
              residual[j][time] = 0.;

              for (k = time - 1;k >= nb_segment;k--) {
                index_parameter_diff = seq_index_parameter[k] - index_parameter_sum / (time - k);
                index_parameter_square_sum += ((double)(time - k) / (double)(time - k + 1)) *
                                               index_parameter_diff * index_parameter_diff;
                diff = int_sequence[j][i][k] - sum / (time - k);
                square_sum += ((double)(time - k) / (double)(time - k + 1)) * diff * diff;
                mix_square_sum += ((double)(time - k) / (double)(time - k + 1)) * index_parameter_diff * diff;
                index_parameter_sum += seq_index_parameter[k];
                sum += int_sequence[j][i][k];

                if ((k < time - 1) && (index_parameter_square_sum > 0.)) {
                  residual[j][k] = square_sum - mix_square_sum * mix_square_sum / index_parameter_square_sum;
                }
                else {
                  residual[j][k] = 0.;
                }
              }
            }
          }
        }

        else {
          for (j = 0;j < nb_sequence;j++) {
            if ((index == I_DEFAULT) || (index == j)) {
              index_parameter_square_sum = 0.;
              square_sum = 0.;
              mix_square_sum = 0.;
              index_parameter_sum = seq_index_parameter[time];
              sum = real_sequence[j][i][time];
              residual[j][time] = 0.;

              for (k = time - 1;k >= nb_segment;k--) {
                index_parameter_diff = seq_index_parameter[k] - index_parameter_sum / (time - k);
                index_parameter_square_sum += ((double)(time - k) / (double)(time - k + 1)) *
                                              index_parameter_diff * index_parameter_diff;
                diff = real_sequence[j][i][k] - sum / (time - k);
                square_sum += ((double)(time - k) / (double)(time - k + 1)) * diff * diff;
                mix_square_sum += ((double)(time - k) / (double)(time - k + 1)) * index_parameter_diff * diff;
                index_parameter_sum += seq_index_parameter[k];
                sum += real_sequence[j][i][k];

                if ((k < time - 1) && (index_parameter_square_sum > 0.)) {
                  residual[j][k] = square_sum - mix_square_sum * mix_square_sum / index_parameter_square_sum;
                }
                else {
                  residual[j][k] = 0.;
                }
              }
            }
          }
        }
      }

      else {
        index_parameter_square_sum = 0.;
        index_parameter_sum = nb_sequence * seq_index_parameter[time];
        square_sum = 0.;
        mix_square_sum = 0.;
        count = 1;
        residual[0][time] = 0.;

        if (type[i] != REAL_VALUE) {
          sum = int_sequence[0][i][time];
          for (j = 1;j < nb_sequence;j++) {
            diff = int_sequence[j][i][time] - sum / count;
            square_sum += ((double)count / (double)(count + 1)) * diff * diff;
            count++;
            sum += int_sequence[j][i][time];
          }

          for (j = time - 1;j >= nb_segment;j--) {
            for (k = 0;k < nb_sequence;k++) {
              index_parameter_diff = seq_index_parameter[j] - index_parameter_sum / count;
              index_parameter_square_sum += ((double)count / (double)(count + 1)) *
                                            index_parameter_diff * index_parameter_diff;
              diff = int_sequence[k][i][j] - sum / count;
              square_sum += ((double)count / (double)(count + 1)) * diff * diff;
              mix_square_sum += ((double)count / (double)(count + 1)) * index_parameter_diff * diff;
              count++;
              index_parameter_sum += seq_index_parameter[j];
              sum += int_sequence[k][i][j];
            }

            if (index_parameter_square_sum > 0.) {
              residual[0][j] = square_sum - mix_square_sum * mix_square_sum / index_parameter_square_sum;
            }
            else {
              residual[0][j] = 0.;
            }
          }
        }

        else {
          sum = real_sequence[0][i][time];
          for (j = 1;j < nb_sequence;j++) {
            diff = real_sequence[j][i][time] - sum / count;
            square_sum += ((double)count / (double)(count + 1)) * diff * diff;
            count++;
            sum += real_sequence[j][i][time];
          }

          for (j = time - 1;j >= nb_segment;j--) {
            for (k = 0;k < nb_sequence;k++) {
              index_parameter_diff = seq_index_parameter[j] - index_parameter_sum / count;
              index_parameter_square_sum += ((double)count / (double)(count + 1)) *
                                            index_parameter_diff * index_parameter_diff;
              diff = real_sequence[k][i][j] - sum / count;
              square_sum += ((double)count / (double)(count + 1)) * diff * diff;
              mix_square_sum += ((double)count / (double)(count + 1)) * index_parameter_diff * diff;
              count++;
              index_parameter_sum += seq_index_parameter[j];
              sum += real_sequence[k][i][j];
            }

            if (index_parameter_square_sum > 0.) {
              residual[0][j] = square_sum - mix_square_sum * mix_square_sum / index_parameter_square_sum;
            }
            else {
              residual[0][j] = 0.;
            }
          }
        }
      }
    }

    else if (model_type[i - 1] == AUTOREGRESSIVE_MODEL_CHANGE) {
      if ((index != I_DEFAULT) || (!common_contrast)) {
        if (type[i] != REAL_VALUE) {
          for (j = 0;j < nb_sequence;j++) {
            if ((index == I_DEFAULT) || (index == j)) {
              sum = int_sequence[j][i][time];

#             ifdef DEBUG
              if (time == 10) {
                cout << "\n";
              }
#             endif

              if (time - 1 >= nb_segment) {
                diff = int_sequence[j][i][time - 1] - int_sequence[j][i][time];
                square_sum = diff * diff / 4.;
                shifted_square_sum = square_sum;
                autocovariance = -square_sum;
                sum += int_sequence[j][i][time - 1];
                residual[j][time - 1] = 0.;

#               ifdef DEBUG
                if (time == 10) {
                  cout << time - 1 << "   " << square_sum << " " << shifted_square_sum << " " << autocovariance << endl;
                }
#               endif
              }

              for (k = time - 2;k >= nb_segment;k--) {
                diff = int_sequence[j][i][k + 1] - sum / (time - k);
                shifted_diff = int_sequence[j][i][k] - sum / (time - k);
                square_sum += diff * diff +
                              ((double)(time - k) / ((double)(time - k + 1) * (time - k + 1))) * shifted_diff * shifted_diff;
                shifted_square_sum += (1. + (double)(time - k) / ((double)(time - k + 1) * (time - k + 1))) * shifted_diff * shifted_diff -
                                      (2. / (double)(time - k + 1)) * shifted_diff * (int_sequence[j][i][k] - int_sequence[j][i][time]);
                autocovariance += diff * shifted_diff +
                                  ((double)(time - k) / ((double)(time - k + 1) * (time - k + 1))) * shifted_diff * shifted_diff -
                                  (1. / (double)(time - k + 1)) * shifted_diff * (int_sequence[j][i][k] - int_sequence[j][i][time]);
                sum += int_sequence[j][i][k];

                residual[j][k] = square_sum;
                if (shifted_square_sum > 0.) {
                  residual[j][k] -= autocovariance * autocovariance / shifted_square_sum;
                }

#               ifdef DEBUG
                if (time == 10) {
                  cout << k << "   " << square_sum << " " << shifted_square_sum << " " << autocovariance << endl;
                }
#               endif

              }

#             ifdef DEBUG
              if (time == 10) {
                cout << "\n";
              }

              sum = int_sequence[j][i][time];

              for (k = time - 1;k >= nb_segment;k--) {
                sum += int_sequence[j][i][k];
                mean = sum / (time - k + 1);

                square_sum = 0.;
                shifted_square_sum = 0.;
                autocovariance = 0.;
                for (m = k + 1;m <= time;m++) {
                  diff = int_sequence[j][i][m] - mean;
                  shifted_diff = int_sequence[j][i][m - 1] - mean;
                  square_sum += diff * diff;
                  shifted_square_sum += shifted_diff * shifted_diff;
                  autocovariance += diff * shifted_diff;
                }

                buff = square_sum;
                if (shifted_square_sum > 0.) {
                  buff -= autocovariance * autocovariance / shifted_square_sum;
                }

                if (time == 10) {
                  cout << k << "   " << square_sum << " " << shifted_square_sum << " " << autocovariance << "   "
                       << residual[j][k] << " | " << buff << endl; 
                }
              }
#             endif

            }
          }
        }

        else {
          for (j = 0;j < nb_sequence;j++) {
            if ((index == I_DEFAULT) || (index == j)) {
              sum = real_sequence[j][i][time];

              if (time - 1 >= nb_segment) {
                diff = real_sequence[j][i][time - 1] - real_sequence[j][i][time];
                square_sum = diff * diff / 4.;
                shifted_square_sum = square_sum;
                autocovariance = -square_sum;
                sum += real_sequence[j][i][time - 1];
                residual[j][time - 1] = 0.;
              }

              for (k = time - 2;k >= nb_segment;k--) {
                diff = real_sequence[j][i][k + 1] - sum / (time - k);
                shifted_diff = real_sequence[j][i][k] - sum / (time - k);
                square_sum += diff * diff +
                              ((double)(time - k) / ((double)(time - k + 1) * (time - k + 1))) * shifted_diff * shifted_diff;
                shifted_square_sum += (1. + (double)(time - k) / ((double)(time - k + 1) * (time - k + 1))) * shifted_diff * shifted_diff -
                                      (2. / (double)(time - k + 1)) * shifted_diff * (real_sequence[j][i][k] - real_sequence[j][i][time]);
                autocovariance += diff * shifted_diff +
                                  ((double)(time - k) / ((double)(time - k + 1) * (time - k + 1))) * shifted_diff * shifted_diff -
                                  (1. / (double)(time - k + 1)) * shifted_diff * (real_sequence[j][i][k] - real_sequence[j][i][time]);
                sum += real_sequence[j][i][k];

                residual[j][k] = square_sum;
                if (shifted_square_sum > 0.) {
                  residual[j][k] -= autocovariance * autocovariance / shifted_square_sum;
                }
              }
            }
          }
        }
      }

      else {
        if (type[i] != REAL_VALUE) {
          sum = 0.;
          for (j = 0;j < nb_sequence;j++) {
            sum += int_sequence[j][i][time];
          }

          if (time - 1 >= nb_segment) {
            for (j = 0;j < nb_sequence;j++) {
              sum += int_sequence[j][i][time - 1];
            }
            mean = sum / (nb_sequence * 2);

            square_sum = 0.;
            shifted_square_sum = 0.;
            autocovariance = 0.;
            for (j = 0;j < nb_sequence;j++) {
              diff = int_sequence[j][i][time] - mean;
              shifted_diff = int_sequence[j][i][time - 1] - mean;
              square_sum += diff * diff;
              shifted_square_sum += shifted_diff * shifted_diff;
              autocovariance += diff * shifted_diff;
            }

            residual[0][time - 1] = square_sum;
            if (shifted_square_sum > 0.) {
              residual[0][time - 1] -= autocovariance * autocovariance / shifted_square_sum;
            }
          }

          for (j = time - 2;j >= nb_segment;j--) {
            mean = sum / (nb_sequence * (time - j));
            square_sum_term[0] = 0.;
            square_sum_term[1] = 0.;
            square_sum_term[2] = 0.;
            shifted_diff = 0.;
            range_diff = 0.;

            for (k = 0;k < nb_sequence;k++) {
              sum += int_sequence[k][i][j];
              diff = int_sequence[k][i][j + 1] - mean;
              square_sum_term[0] += diff * diff;
              buff = int_sequence[k][i][j] - mean;
              shifted_diff += buff;
              square_sum_term[1] += buff * buff;
              square_sum_term[2] += diff * buff;
              range_diff += int_sequence[k][i][j] - int_sequence[k][i][time];
            }

            square_sum += square_sum_term[0] +
                          ((double)(time - j) / ((double)nb_sequence * (time - j + 1) * (time - j + 1))) * shifted_diff * shifted_diff;
            shifted_square_sum += square_sum_term[1] +
                                  ((double)(time - j) / ((double)nb_sequence * (time - j + 1) * (time - j + 1))) * shifted_diff * shifted_diff -
                                  (2. / ((double)nb_sequence * (time - j + 1))) * shifted_diff * range_diff;
            autocovariance += square_sum_term[2] +
                              ((double)(time - j) / ((double)nb_sequence * (time - j + 1) * (time - j + 1))) * shifted_diff * shifted_diff -
                              (1. / ((double)nb_sequence * (time - j + 1))) * shifted_diff * range_diff;

            residual[0][j] = square_sum;
            if (shifted_square_sum > 0.) {
              residual[0][j] -= autocovariance * autocovariance / shifted_square_sum;
            }
          }
        }

        else {
          sum = 0.;
          for (j = 0;j < nb_sequence;j++) {
            sum += real_sequence[j][i][time];
          }

          if (time - 1 >= nb_segment) {
            for (j = 0;j < nb_sequence;j++) {
              sum += real_sequence[j][i][time - 1];
            }
            mean = sum / (nb_sequence * 2);

            square_sum = 0.;
            shifted_square_sum = 0.;
            autocovariance = 0.;
            for (j = 0;j < nb_sequence;j++) {
              diff = real_sequence[j][i][time] - mean;
              shifted_diff = real_sequence[j][i][time - 1] - mean;
              square_sum += diff * diff;
              shifted_square_sum += shifted_diff * shifted_diff;
              autocovariance += diff * shifted_diff;
            }

            residual[0][time - 1] = square_sum;
            if (shifted_square_sum > 0.) {
              residual[0][time - 1] -= autocovariance * autocovariance / shifted_square_sum;
            }
          }

          for (j = time - 2;j >= nb_segment;j--) {
            mean = sum / (nb_sequence * (time - j));
            square_sum_term[0] = 0.;
            square_sum_term[1] = 0.;
            square_sum_term[2] = 0.;
            shifted_diff = 0.;
            range_diff = 0.;

            for (k = 0;k < nb_sequence;k++) {
              sum += real_sequence[k][i][j];
              diff = real_sequence[k][i][j + 1] - mean;
              square_sum_term[0] += diff * diff;
              buff = real_sequence[k][i][j] - mean;
              shifted_diff += buff;
              square_sum_term[1] += buff * buff;
              square_sum_term[2] += diff * buff;
              range_diff += real_sequence[k][i][j] - real_sequence[k][i][time];
            }

            square_sum += square_sum_term[0] +
                          ((double)(time - j) / ((double)nb_sequence * (time - j + 1) * (time - j + 1))) * shifted_diff * shifted_diff;
            shifted_square_sum += square_sum_term[1] +
                                  ((double)(time - j) / ((double)nb_sequence * (time - j + 1) * (time - j + 1))) * shifted_diff * shifted_diff -
                                  (2. / ((double)nb_sequence * (time - j + 1))) * shifted_diff * range_diff;
            autocovariance += square_sum_term[2] +
                              ((double)(time - j) / ((double)nb_sequence * (time - j + 1) * (time - j + 1))) * shifted_diff * shifted_diff -
                              (1. / ((double)nb_sequence * (time - j + 1))) * shifted_diff * range_diff;

            residual[0][j] = square_sum;
            if (shifted_square_sum > 0.) {
              residual[0][j] -= autocovariance * autocovariance / shifted_square_sum;
            }
          }
        }
      }
    }

    else if (model_type[i - 1] == STATIONARY_AUTOREGRESSIVE_MODEL_CHANGE) {
      if ((index != I_DEFAULT) || (!common_contrast)) {
        if (type[i] != REAL_VALUE) {
          for (j = 0;j < nb_sequence;j++) {
            if ((index == I_DEFAULT) || (index == j)) {
              square_sum = 0.;
              shifted_square_sum = 0.;
              autocovariance = 0.;

              for (k = time - 1;k >= nb_segment;k--) {
                diff = int_sequence[j][i][k + 1] - seq_mean[i][j];
                shifted_diff = int_sequence[j][i][k] - seq_mean[i][j];
                square_sum += diff * diff;
                shifted_square_sum += shifted_diff * shifted_diff;
                autocovariance += diff * shifted_diff;

                residual[j][k] = square_sum;
                if (shifted_square_sum > 0.) {
                  residual[j][k] -= autocovariance * autocovariance / shifted_square_sum;
                }
              }
            }
          }
        }

        else {
          for (j = 0;j < nb_sequence;j++) {
            if ((index == I_DEFAULT) || (index == j)) {
              square_sum = 0.;
              shifted_square_sum = 0.;
              autocovariance = 0.;

              for (k = time - 1;k >= nb_segment;k--) {
                diff = real_sequence[j][i][k + 1] - seq_mean[i][j];
                shifted_diff = real_sequence[j][i][k] - seq_mean[i][j];
                square_sum += diff * diff;
                shifted_square_sum += shifted_diff * shifted_diff;
                autocovariance += diff * shifted_diff;

                residual[j][k] = square_sum;
                if (shifted_square_sum > 0.) {
                  residual[j][k] -= autocovariance * autocovariance / shifted_square_sum;
                }
              }
            }
          }
        }
      }

      else {
        square_sum = 0.;
        shifted_square_sum = 0.;
        autocovariance = 0.;

        if (type[i] != REAL_VALUE) {
          for (j = time - 1;j >= nb_segment;j--) {
            for (k = 0;k < nb_sequence;k++) {
              diff = int_sequence[k][i][j + 1] - seq_mean[i][0];
              shifted_diff = int_sequence[k][i][j] - seq_mean[i][0];
              square_sum += diff * diff;
              shifted_square_sum += shifted_diff * shifted_diff;
              autocovariance += diff * shifted_diff;
            }

            residual[0][j] = square_sum;
            if (shifted_square_sum > 0.) {
              residual[0][j] -= autocovariance * autocovariance / shifted_square_sum;
            }
          }
        }

        else {
          for (j = time - 1;j >= nb_segment;j--) {
            for (k = 0;k < nb_sequence;k++) {
              diff = real_sequence[k][i][j + 1] - seq_mean[i][0];
              shifted_diff = real_sequence[k][i][j] - seq_mean[i][0];
              square_sum += diff * diff;
              shifted_square_sum += shifted_diff * shifted_diff;
              autocovariance += diff * shifted_diff;
            }

            residual[0][j] = square_sum;
            if (shifted_square_sum > 0.) {
              residual[0][j] -= autocovariance * autocovariance / shifted_square_sum;
            }
          }
        }
      }
    }

    else if (model_type[i - 1] == BAYESIAN_POISSON_CHANGE) {
      prior_contrast = -lgamma(hyperparam[i][0]) + hyperparam[i][0] * log(hyperparam[i][1]);
      factorial[i][index][time] = log_factorial(int_sequence[index][i][time]);

      sum = 0.;
      factorial_sum = 0.;
      for (j = time;j >= nb_segment;j--) {
        sum += int_sequence[index][i][j];
        factorial_sum += factorial[i][index][j];
        if (contrast[j] != D_INF) {
          contrast[j] += prior_contrast - factorial_sum + lgamma(hyperparam[i][0] + sum) -
                         (hyperparam[i][0] + sum) * log(hyperparam[i][1] + time - j + 1);
        }
      }
    }

    else if (model_type[i - 1] == BAYESIAN_GAUSSIAN_CHANGE) {
      prior_contrast = log(hyperparam[i][1]) / 2 - lgamma(hyperparam[i][2] / 2) +
                       hyperparam[i][2] * log(hyperparam[i][3] / 2) / 2;

      if (type[i] != REAL_VALUE) {
        square_sum = 0.;
        sum = int_sequence[index][i][time];
        if (contrast[time] != D_INF) {
          diff = hyperparam[i][0] - sum;
          contrast[time] += prior_contrast - log(2 * M_PI) / 2 -
                            log(hyperparam[i][1] + 1) / 2 + lgamma((hyperparam[i][2] + 1) / 2) -
                            (hyperparam[i][2] + 1) *
                            log((hyperparam[i][3] + hyperparam[i][1] *
                                 diff * diff / (hyperparam[i][1] + 1)) / 2) / 2;
        }

        for (j = time - 1;j >= nb_segment;j--) {
          diff = int_sequence[index][i][j] - sum / (time - j);
          square_sum += ((double)(time - j) / (double)(time - j + 1)) * diff * diff;
          sum += int_sequence[index][i][j];
          if (contrast[j] != D_INF) {
            diff = hyperparam[i][0] - sum / (time - j + 1);
            contrast[j] += prior_contrast - (time - j + 1) * log(2 * M_PI) / 2 -
                           log(hyperparam[i][1] + time - j + 1) / 2 +
                           lgamma((hyperparam[i][2] + time - j + 1) / 2) -
                           (hyperparam[i][2] + time - j + 1) *
                           logl((hyperparam[i][3] + square_sum + hyperparam[i][1] * (time - j + 1) *
                                 diff * diff / (hyperparam[i][1] + time - j + 1)) / 2) / 2;
          }
        }
      }

      else {
        square_sum = 0.;
        sum = real_sequence[index][i][time];
        if (contrast[time] != D_INF) {
          diff = hyperparam[i][0] - sum;
          contrast[time] += prior_contrast - log(2 * M_PI) / 2 -
                            log(hyperparam[i][1] + 1) / 2 + lgamma((hyperparam[i][2] + 1) / 2) -
                            (hyperparam[i][2] + 1) *
                            log((hyperparam[i][3] + hyperparam[i][1] *
                                 diff * diff / (hyperparam[i][1] + 1)) / 2) / 2;
        }

        for (j = time - 1;j >= nb_segment;j--) {
          diff = real_sequence[index][i][j] - sum / (time - j);
          square_sum += ((double)(time - j) / (double)(time - j + 1)) * diff * diff;
          sum += real_sequence[index][i][j];
          if (contrast[j] != D_INF) {
            diff = hyperparam[i][0] - sum / (time - j + 1);
            contrast[j] += prior_contrast - (time - j + 1) * log(2 * M_PI) / 2 -
                           log(hyperparam[i][1] + time - j + 1) / 2 +
                           lgamma((hyperparam[i][2] + time - j + 1) / 2) -
                           (hyperparam[i][2] + time - j + 1) *
                           logl((hyperparam[i][3] + square_sum + hyperparam[i][1] * (time - j + 1) *
                                 diff * diff / (hyperparam[i][1] + time - j + 1)) / 2) / 2;
          }
        }
      }
    }

    if ((model_type[0] == MEAN_CHANGE) || (model_type[0] == INTERCEPT_SLOPE_CHANGE)) {
      if ((index != I_DEFAULT) || (!common_contrast)) {
        for (j = 0;j < nb_sequence;j++) {
          if ((index == I_DEFAULT) || (index == j)) {
            for (k = time - 1;k >= nb_segment;k--) {
              contrast[k] -= residual[j][k];
            }
          }
        }
      }

      else {
        for (j = time;j >= nb_segment;j--) {
          contrast[j] -= residual[0][j];
        }
      }
    }

    else if ((model_type[i - 1] == GAUSSIAN_CHANGE) || (model_type[i - 1] == VARIANCE_CHANGE) ||
             (model_type[i - 1] == ORDINAL_GAUSSIAN_CHANGE) || (model_type[i - 1] == LINEAR_MODEL_CHANGE)) {
      if ((index != I_DEFAULT) || (!common_contrast)) {
        for (j = 0;j < nb_sequence;j++) {
          if ((index == I_DEFAULT) || (index == j)) {
            for (k = time;k >= nb_segment;k--) {
              if (contrast[k] != D_INF) {
//                if (residual[j][k] > 0.) {
                if (residual[j][k] > (time - k + 1) * ROUNDOFF_ERROR) {
                  contrast[k] -= ((double)(time - k + 1) / 2.) * (logl(residual[j][k] /
                                   (time - k + 1)) + log(2 * M_PI) + 1);
/*                  contrast[k] -= ((double)(time - k + 1) / 2.) * (logl(residual[j][k] /
                                   (time - k)) + log(2 * M_PI)) + (double)(time - k) / 2.; */
                }
                else {
                  contrast[k] = D_INF;
                }
              }
            }
          }
        }
      }

      else {
        for (j = time;j >= nb_segment;j--) {
          if (contrast[j] != D_INF) {
//            if (residual[0][j] > 0.) {
            if (residual[0][j] > nb_sequence * (time - j + 1) * ROUNDOFF_ERROR) {
              contrast[j] -= ((double)(nb_sequence * (time - j + 1)) / 2.) * (logl(residual[0][j] /
                               (nb_sequence * (time - j + 1))) + log(2 * M_PI) + 1);
            }
            else {
              contrast[j] = D_INF;
            }
          }
        }
      }
    }

    else if ((model_type[i - 1] == AUTOREGRESSIVE_MODEL_CHANGE) ||
             (model_type[i - 1] == STATIONARY_AUTOREGRESSIVE_MODEL_CHANGE)) {
      if ((index != I_DEFAULT) || (!common_contrast)) {
        for (j = 0;j < nb_sequence;j++) {
          if ((index == I_DEFAULT) || (index == j)) {
            for (k = time - 1;k >= nb_segment;k--) {
              if (contrast[k] != D_INF) {
//                if (residual[j][k] > 0.) {
                if (residual[j][k] > (time - k) * ROUNDOFF_ERROR) {
                  contrast[k] -= ((double)(time - k) / 2.) * (logl(residual[j][k] /
                                   (time - k)) + log(2 * M_PI) + 1);
                }
                else {
                  contrast[k] = D_INF;
                }
              }
            }
          }
        }
      }

      else {
        for (j = time - 1;j >= nb_segment;j--) {
          if (contrast[j] != D_INF) {
//            if (residual[0][j] > 0.) {
            if (residual[0][j] > nb_sequence * (time - j) * ROUNDOFF_ERROR) {
              contrast[j] -= ((double)(nb_sequence * (time - j)) / 2.) * (logl(residual[0][j] /
                               (nb_sequence * (time - j))) + log(2 * M_PI) + 1);
            }
            else {
              contrast[j] = D_INF;
            }
          }
        }
      }
    }
  }

# ifdef DEBUG
  for (i = time - 1;i >= nb_segment;i--) {
    cout << contrast[i] << "  ";
  }
  cout << endl;
# endif

  delete [] frequency;
  delete [] inf_bound_parameter;

  if (residual) {
    if ((index != I_DEFAULT) || (!common_contrast)) {
      for (i = 0;i < nb_sequence;i++) {
        delete [] residual[i];
      }
    }
    else {
      delete [] residual[0];
    }
    delete [] residual;
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the contrast functions within a backward recursion.
 *
 *  \param[in] time                time instant,
 *  \param[in] index               sequence index,
 *  \param[in] model_type          segment model types,
 *  \param[in] common_contrast     flag contrast functions common to the individuals,
 *  \param[in] factorial           log factorials for Poisson models,
 *  \param[in] shape_parameter     negative binomial shape parameters,
 *  \param[in] binomial_coeff      log binomial coefficients for negative binomial models,
 *  \param[in] seq_mean            sequence means for Gaussian change in the variance models or
 *                                 stationary piecewise autoregressive models,
 *  \param[in] seq_index_parameter index parameters,
 *  \param[in] hyperparam          hyperparameters for Bayesian models,
 *  \param[in] rank                ranks for ordinal variables,
 *  \param[in] contrast            contrast functions.
 */
/*--------------------------------------------------------------*/

void Sequences::backward_contrast(int time , int index , segment_model *model_type , bool common_contrast ,
                                  double ***factorial , double *shape_parameter , double ***binomial_coeff ,
                                  double **seq_mean , int *seq_index_parameter , double **hyperparam ,
                                  double **rank , long double *contrast) const

{
  int i , j , k , m;
  int max_nb_value , count , *frequency , *inf_bound_parameter;
  double sum , factorial_sum , proba , binomial_coeff_sum , diff , index_parameter_sum ,
         index_parameter_diff , shifted_diff , range_diff , mean , buff;
  long double index_parameter_square_sum , square_sum , mix_square_sum , shifted_square_sum ,
              autocovariance , prior_contrast , square_sum_term[3] , **residual;


  // initializations

  max_nb_value = 0;
  inf_bound_parameter = new int[nb_variable];
  residual = NULL;

  for (i = 1;i < nb_variable;i++) {
    if ((model_type[i - 1] == CATEGORICAL_CHANGE) && (marginal_distribution[i]->nb_value > max_nb_value)) {
      max_nb_value = marginal_distribution[i]->nb_value;
    }

    if ((model_type[i - 1] == NEGATIVE_BINOMIAL_0_CHANGE) || (model_type[i - 1] == NEGATIVE_BINOMIAL_1_CHANGE)) {
      switch (model_type[i - 1]) {
      case NEGATIVE_BINOMIAL_0_CHANGE :
        inf_bound_parameter[i - 1] = 0;
        break;
      case NEGATIVE_BINOMIAL_1_CHANGE :
        inf_bound_parameter[i - 1] = 1;
        break;
      }
    }

    if (((i == 1) && ((model_type[0] == MEAN_CHANGE) || (model_type[0] == INTERCEPT_SLOPE_CHANGE))) ||
        (((model_type[i - 1] == GAUSSIAN_CHANGE) || (model_type[i - 1] == VARIANCE_CHANGE) ||
          (model_type[i - 1] == ORDINAL_GAUSSIAN_CHANGE) || (model_type[i - 1] == LINEAR_MODEL_CHANGE) ||
          (model_type[i - 1] == STATIONARY_AUTOREGRESSIVE_MODEL_CHANGE) ||
          (model_type[i - 1] == AUTOREGRESSIVE_MODEL_CHANGE)) && (!residual))) {
      residual = new long double*[MAX(nb_sequence , 2)];
      if ((index != I_DEFAULT) || (!common_contrast)) {
        for (j = 0;j < nb_sequence;j++) {
          if ((index == I_DEFAULT) || (index == j)) {
            residual[j] = new long double[length[j]];
          }
          else {
            residual[j] = NULL;
          }
        }
      }
      else {
        residual[0] = new long double[length[0]];
      }
    }
  }

  if (max_nb_value > 0) {
    frequency = new int[max_nb_value];
  }
  else {
    frequency = NULL;
  }

  // computation of segment contrast functions (log-likelihoods or sum of squared deviations)

  for (i = time;i < length[index == I_DEFAULT ? 0 : index];i++) {
    contrast[i] = 0.;
  }

  for (i = 1;i < nb_variable;i++) {
    if (model_type[i - 1] == CATEGORICAL_CHANGE) {
      if ((index != I_DEFAULT) || (!common_contrast)) {
        for (j = 0;j < nb_sequence;j++) {
          if ((index == I_DEFAULT) || (index == j)) {
            for (k = 0;k < marginal_distribution[i]->nb_value;k++) {
              frequency[k] = 0;
            }

            frequency[int_sequence[j][i][time]]++;
//            sum = 0.;

            for (k = time + 1;k < length[j];k++) {
              frequency[int_sequence[j][i][k]]++;
              if (contrast[k] != D_INF) {
                for (m = 0;m < marginal_distribution[i]->nb_value;m++) {
                  if (frequency[m] > 0) {
                    contrast[k] += frequency[m] * log((double)frequency[m] / (double)(k - time + 1));
                  }
                }
              }

/*              sum += (k - time) * log((double)(k - time) / (double)(k - time + 1)) +
                     log((double)(frequency[int_sequence[j][i][k]] + 1) / (double)(k - time + 1));
              if (frequency[int_sequence[j][i][k]] > 0) {
                sum -= frequency[int_sequence[j][i][k]] *
                       log((double)frequency[int_sequence[j][i][k]] / (double)(frequency[int_sequence[j][i][k]] + 1));
              }
              frequency[int_sequence[j][i][k]]++;

              if (contrast[k] != D_INF) {
                contrast[k] += sum;
              } */
            }
          }
        }
      }

      else {
        for (j = 0;j < marginal_distribution[i]->nb_value;j++) {
          frequency[j] = 0;
        }

        for (j = time;j < length[0];j++) {
          for (k = 0;k < nb_sequence;k++) {
            frequency[int_sequence[k][i][j]]++;
          }

          if (contrast[j] != D_INF) {
            for (k = 0;k < marginal_distribution[i]->nb_value;k++) {
              if (frequency[k] > 0) {
                contrast[j] += frequency[k] * log((double)frequency[k] / (double)(nb_sequence * (j - time + 1)));
              }
            }
          }
        }
      }
    }

    else if (model_type[i - 1] == POISSON_CHANGE) {
      if ((index != I_DEFAULT) || (!common_contrast)) {
        for (j = 0;j < nb_sequence;j++) {
          if ((index == I_DEFAULT) || (index == j)) {
            sum = 0.;
            factorial_sum = 0.;

            for (k = time;k < length[j];k++) {
              sum += int_sequence[j][i][k];
              factorial_sum += factorial[i][j][k];
              if ((contrast[k] != D_INF) && (sum > 0.)) {
                contrast[k] += sum * (log(sum / (k - time + 1)) - 1) - factorial_sum;
              }
            }
          }
        }
      }

      else {
        sum = 0.;
        factorial_sum = 0.;

        for (j = time;j < length[0];j++) {
          for (k = 0;k < nb_sequence;k++) {
            sum += int_sequence[k][i][j];
            factorial_sum += factorial[i][k][j];
          }

          if ((contrast[j] != D_INF) && (sum > 0.)) {
            contrast[j] += sum * (log(sum / (nb_sequence * (j - time + 1))) - 1) - factorial_sum;
          }
        }
      }
    }

    else if ((model_type[i - 1] == NEGATIVE_BINOMIAL_0_CHANGE) || (model_type[i - 1] == NEGATIVE_BINOMIAL_1_CHANGE)) {
      if ((index != I_DEFAULT) || (!common_contrast)) {
        for (j = 0;j < nb_sequence;j++) {
          if ((index == I_DEFAULT) || (index == j)) {
            sum = 0.;
            binomial_coeff_sum = 0.;

            for (k = time;k < length[j];k++) {
              sum += int_sequence[j][i][k];
              binomial_coeff_sum += binomial_coeff[i][j][k];

              if (contrast[k] != D_INF) {
                if (sum > inf_bound_parameter[i - 1] * (k - time + 1)) {
                  proba = shape_parameter[i - 1] * (k - time + 1) /
                          ((shape_parameter[i - 1] - inf_bound_parameter[i - 1]) * (k - time + 1) + sum);
                  contrast[k] += binomial_coeff_sum + shape_parameter[i - 1] * (k - time + 1) * log(proba) +
                                 (sum - inf_bound_parameter[i - 1] * (k - time + 1)) * log(1. - proba);
                }
                else {
                  contrast[k] = D_INF;
                }
              }
            }
          }
        }
      }

      else {
        sum = 0.;
        binomial_coeff_sum = 0.;

        for (j = time;j < length[0];j++) {
          for (k = 0;k < nb_sequence;k++) {
            sum += int_sequence[k][i][j];
            binomial_coeff_sum += binomial_coeff[i][k][j];
          }

          if (contrast[j] != D_INF) {
            if (sum > inf_bound_parameter[i - 1] * nb_sequence * (j - time + 1)) {
              proba = shape_parameter[i - 1] * nb_sequence * (j - time + 1) /
                      ((shape_parameter[i - 1] - inf_bound_parameter[i - 1]) * nb_sequence * (j - time + 1) + sum);
              contrast[j] += binomial_coeff_sum + shape_parameter[i - 1] * nb_sequence * (j - time + 1) * log(proba) +
                             (sum - inf_bound_parameter[i - 1] * nb_sequence * (j - time + 1)) * log(1. - proba);
            }
            else {
              contrast[j] = D_INF;
            }
          }
        }
      }
    }

    else if ((model_type[i - 1] == GAUSSIAN_CHANGE) || (model_type[0] == MEAN_CHANGE)) {
      if ((index != I_DEFAULT) || (!common_contrast)) {
        if (type[i] != REAL_VALUE) {
          for (j = 0;j < nb_sequence;j++) {
            if ((index == I_DEFAULT) || (index == j)) {
              square_sum = 0.;
              sum = int_sequence[j][i][time];
              residual[j][time] = 0.;

              for (k = time + 1;k < length[j];k++) {
                diff = int_sequence[j][i][k] - sum / (k - time);
                square_sum += ((double)(k - time) / (double)(k - time + 1)) * diff * diff;
                sum += int_sequence[j][i][k];
                residual[j][k] = square_sum;
              }
            }
          }
        }

        else {
          for (j = 0;j < nb_sequence;j++) {
            if ((index == I_DEFAULT) || (index == j)) {
              square_sum = 0.;
              sum = real_sequence[j][i][time];
              residual[j][time] = 0.;

              for (k = time + 1;k < length[j];k++) {
                diff = real_sequence[j][i][k] - sum / (k - time);
                square_sum += ((double)(k - time) / (double)(k - time + 1)) * diff * diff;
                sum += real_sequence[j][i][k];
                residual[j][k] = square_sum;
              }
            }
          }
        }
      }

      else {
        square_sum = 0.;
        sum = 0.;
        count = 0;

        if (type[i] != REAL_VALUE) {
          for (j = time;j < length[0];j++) {
            for (k = 0;k < nb_sequence;k++) {
              if (count > 0) {
                diff = int_sequence[k][i][j] - sum / count;
                square_sum += ((double)count / (double)(count + 1)) * diff * diff;
              }
              count++;
              sum += int_sequence[k][i][j];
            }
            residual[0][j] = square_sum;
          }
        }

        else {
          for (j = time;j < length[0];j++) {
            for (k = 0;k < nb_sequence;k++) {
              if (count > 0) {
                diff = real_sequence[k][i][j] - sum / count;
                square_sum += ((double)count / (double)(count + 1)) * diff * diff;
              }
              count++;
              sum += real_sequence[k][i][j];
            }
            residual[0][j] = square_sum;
          }
        }
      }
    }

    else if (model_type[i - 1] == VARIANCE_CHANGE) {
      if ((index != I_DEFAULT) || (!common_contrast)) {
        if (type[i] != REAL_VALUE) {
          for (j = 0;j < nb_sequence;j++) {
            if ((index == I_DEFAULT) || (index == j)) {
              square_sum = 0.;
              for (k = time;k < length[j];k++) {
                diff = int_sequence[j][i][k] - seq_mean[i][j];
                square_sum += diff * diff;
                residual[j][k] = square_sum;
              }
            }
          }
        }

        else {
          for (j = 0;j < nb_sequence;j++) {
            if ((index == I_DEFAULT) || (index == j)) {
              square_sum = 0.;
              for (k = time;k < length[j];k++) {
                diff = real_sequence[j][i][k] - seq_mean[i][j];
                square_sum += diff * diff;
                residual[j][k] = square_sum;
              }
            }
          }
        }
      }

      else {
        square_sum = 0.;

        if (type[i] != REAL_VALUE) {
          for (j = time;j < length[0];j++) {
            for (k = 0;k < nb_sequence;k++) {
              diff = int_sequence[k][i][j] - seq_mean[i][0];
              square_sum += diff * diff;
            }
            residual[0][j] = square_sum;
          }
        }

        else {
          for (j = time;j < length[0];j++) {
            for (k = 0;k < nb_sequence;k++) {
              diff = real_sequence[k][i][j] - seq_mean[i][0];
              square_sum += diff * diff;
            }
            residual[0][j] = square_sum;
          }
        }
      }
    }

    else if (model_type[i - 1] == ORDINAL_GAUSSIAN_CHANGE) {
      if ((index != I_DEFAULT) || (!common_contrast)) {
        for (j = 0;j < nb_sequence;j++) {
          if ((index == I_DEFAULT) || (index == j)) {
            square_sum = 0.;
            sum = rank[i][int_sequence[j][i][time]];
            residual[j][time] = 0.;

            for (k = time + 1;k < length[j];k++) {
              diff = rank[i][int_sequence[j][i][k]] - sum / (k - time);
              square_sum += ((double)(k - time) / (double)(k - time + 1)) * diff * diff;
              sum += rank[i][int_sequence[j][i][k]];
              residual[j][k] = square_sum;

              if (residual[j][k] == 0.) {
                residual[j][k] = (k - time + 1) * MIN_RANK_SQUARE_SUM;
              }
            }
          }
        }
      }

      else {
        square_sum = 0.;
        sum = 0.;
        count = 0;

        for (j = time;j < length[0];j++) {
          for (k = 0;k < nb_sequence;k++) {
            if (count > 0) {
              diff = rank[i][int_sequence[k][i][j]] - sum / count;
              square_sum += ((double)count / (double)(count + 1)) * diff * diff;
            }
            count++;
            sum += rank[i][int_sequence[k][i][j]];
          }
          residual[0][j] = square_sum;

          if (residual[0][j] == 0.) {
            residual[0][j] = count * MIN_RANK_SQUARE_SUM;
          }
        }
      }
    }

    else if ((model_type[i - 1] == LINEAR_MODEL_CHANGE) || (model_type[0] == INTERCEPT_SLOPE_CHANGE)) {
      if ((index != I_DEFAULT) || (!common_contrast)) {
        if (type[i] != REAL_VALUE) {
          for (j = 0;j < nb_sequence;j++) {
            if ((index == I_DEFAULT) || (index == j)) {
              index_parameter_square_sum = 0.;
              square_sum = 0.;
              mix_square_sum = 0.;
              index_parameter_sum = seq_index_parameter[time];
              sum = int_sequence[j][i][time];
              residual[j][time] = 0.;

              for (k = time + 1;k < length[j];k++) {
                index_parameter_diff = seq_index_parameter[k] - index_parameter_sum / (k - time);
                index_parameter_square_sum += ((double)(k - time) / (double)(k - time + 1)) *
                                              index_parameter_diff * index_parameter_diff;
                diff = int_sequence[j][i][k] - sum / (k - time);
                square_sum += ((double)(k - time) / (double)(k - time + 1)) * diff * diff;
                mix_square_sum += ((double)(k - time) / (double)(k - time + 1)) * index_parameter_diff * diff;
                index_parameter_sum += seq_index_parameter[k];
                sum += int_sequence[j][i][k];

                if ((k > time + 1) && (index_parameter_square_sum > 0.)) {
                  residual[j][k] = square_sum - mix_square_sum * mix_square_sum / index_parameter_square_sum;
                }
                else {
                  residual[j][k] = 0.;
                }
              }
            }
          }
        }

        else {
          for (j = 0;j < nb_sequence;j++) {
            if ((index == I_DEFAULT) || (index == j)) {
              index_parameter_square_sum = 0.;
              square_sum = 0.;
              mix_square_sum = 0.;
              index_parameter_sum = seq_index_parameter[time];
              sum = real_sequence[j][i][time];
              residual[j][time] = 0.;

              for (k = time + 1;k < length[j];k++) {
                index_parameter_diff = seq_index_parameter[k] - index_parameter_sum / (k - time);
                index_parameter_square_sum += ((double)(k - time) / (double)(k - time + 1)) *
                                              index_parameter_diff * index_parameter_diff;
                diff = real_sequence[j][i][k] - sum / (k - time);
                square_sum += ((double)(k - time) / (double)(k - time + 1)) * diff * diff;
                mix_square_sum += ((double)(k - time) / (double)(k - time + 1)) * index_parameter_diff * diff;
                index_parameter_sum += seq_index_parameter[k];
                sum += real_sequence[j][i][k];

                if ((k > time + 1) && (index_parameter_square_sum > 0.)) {
                  residual[j][k] = square_sum - mix_square_sum * mix_square_sum / index_parameter_square_sum;
                }
                else {
                  residual[j][k] = 0.;
                }
              }
            }
          }
        }
      }

      else {
        index_parameter_square_sum = 0.;
        index_parameter_sum = nb_sequence * seq_index_parameter[time];
        square_sum = 0.;
        mix_square_sum = 0.;
        count = 1;
        residual[0][time] = 0.;

        if (type[i] != REAL_VALUE) {
          sum = int_sequence[0][i][time];
          for (j = 1;j < nb_sequence;j++) {
            diff = int_sequence[j][i][time] - sum / count;
            square_sum += ((double)count / (double)(count + 1)) * diff * diff;
            count++;
            sum += int_sequence[j][i][time];
          }

          for (j = time + 1;j < length[0];j++) {
            for (k = 0;k < nb_sequence;k++) {
              index_parameter_diff = seq_index_parameter[j] - index_parameter_sum / count;
              index_parameter_square_sum += ((double)count / (double)(count + 1)) *
                                            index_parameter_diff * index_parameter_diff;
              diff = int_sequence[k][i][j] - sum / count;
              square_sum += ((double)count / (double)(count + 1)) * diff * diff;
              mix_square_sum += ((double)count / (double)(count + 1)) * index_parameter_diff * diff;
              count++;
              index_parameter_sum += seq_index_parameter[j];
              sum += int_sequence[k][i][j];
            }

            if (index_parameter_square_sum > 0.) {
              residual[0][j] = square_sum - mix_square_sum * mix_square_sum / index_parameter_square_sum;
            }
            else {
              residual[0][j] = 0.;
            }
          }
        }

        else {
          sum = real_sequence[0][i][time];
          for (j = 1;j < nb_sequence;j++) {
            diff = real_sequence[j][i][time] - sum / count;
            square_sum += ((double)count / (double)(count + 1)) * diff * diff;
            count++;
            sum += real_sequence[j][i][time];
          }

          for (j = time + 1;j < length[0];j++) {
            for (k = 0;k < nb_sequence;k++) {
              index_parameter_diff = seq_index_parameter[j] - index_parameter_sum / count;
              index_parameter_square_sum += ((double)count / (double)(count + 1)) *
                                            index_parameter_diff * index_parameter_diff;
              diff = real_sequence[k][i][j] - sum / count;
              square_sum += ((double)count / (double)(count + 1)) * diff * diff;
              mix_square_sum += ((double)count / (double)(count + 1)) * index_parameter_diff * diff;
              count++;
              index_parameter_sum += seq_index_parameter[j];
              sum += real_sequence[k][i][j];
            }

            if (index_parameter_square_sum > 0.) {
              residual[0][j] = square_sum - mix_square_sum * mix_square_sum / index_parameter_square_sum;
            }
            else {
              residual[0][j] = 0.;
            }
          }
        }
      }
    }

    else if (model_type[i - 1] == AUTOREGRESSIVE_MODEL_CHANGE) {
      if ((index != I_DEFAULT) || (!common_contrast)) {
        if (type[i] != REAL_VALUE) {
          for (j = 0;j < nb_sequence;j++) {
            if ((index == I_DEFAULT) || (index == j)) {
              sum = int_sequence[j][i][time];

              if (time + 1 < length[j]) {
                diff = int_sequence[j][i][time + 1] - int_sequence[j][i][time];
                square_sum = diff * diff / 4.;
                shifted_square_sum = square_sum;
                autocovariance = -square_sum;
                sum += int_sequence[j][i][time + 1];
                residual[j][time + 1] = 0.;
              }

              for (k = time + 2;k < length[j];k++) {
                diff = int_sequence[j][i][k] - sum / (k - time);
                shifted_diff = int_sequence[j][i][k - 1] - sum / (k - time);
                square_sum += (1. + (double)(k - time) / ((double)(k - time + 1) * (k - time + 1))) * diff * diff -
                              (2. / (double)(k - time + 1)) * diff * (int_sequence[j][i][k] - int_sequence[j][i][time]);
                shifted_square_sum += shifted_diff * shifted_diff +
                                      ((double)(k - time) / ((double)(k - time + 1) * (k - time + 1))) * diff * diff;
                autocovariance += diff * shifted_diff +
                                  ((double)(k - time) / ((double)(k - time + 1) * (k - time + 1))) * diff * diff -
                                  (1. / (double)(k - time + 1)) * diff * (int_sequence[j][i][k] - int_sequence[j][i][time]);
                sum += int_sequence[j][i][k];

                residual[j][k] = square_sum;
                if (shifted_square_sum > 0.) {
                  residual[j][k] -= autocovariance * autocovariance / shifted_square_sum;
                }
              }
            }
          }
        }

        else {
          for (j = 0;j < nb_sequence;j++) {
            if ((index == I_DEFAULT) || (index == j)) {
              sum = real_sequence[j][i][time];

              if (time + 1 < length[j]) {
                diff = real_sequence[j][i][time + 1] - real_sequence[j][i][time];
                square_sum = diff * diff / 4.;
                shifted_square_sum = square_sum;
                autocovariance = -square_sum;
                sum += real_sequence[j][i][time + 1];
                residual[j][time + 1] = 0.;
              }

              for (k = time + 2;k < length[j];k++) {
                diff = real_sequence[j][i][k] - sum / (k - time);
                shifted_diff = real_sequence[j][i][k - 1] - sum / (k - time);
                square_sum += (1. + (double)(k - time) / ((double)(k - time + 1) * (k - time + 1))) * diff * diff -
                              (2. / (double)(k - time + 1)) * diff * (real_sequence[j][i][k] - real_sequence[j][i][time]);
                shifted_square_sum += shifted_diff * shifted_diff +
                                      ((double)(k - time) / ((double)(k - time + 1) * (k - time + 1))) * diff * diff;
                autocovariance += diff * shifted_diff +
                                  ((double)(k - time) / ((double)(k - time + 1) * (k - time + 1))) * diff * diff -
                                  (1. / (double)(k - time + 1)) * diff * (real_sequence[j][i][k] - real_sequence[j][i][time]);
                sum += real_sequence[j][i][k];

                residual[j][k] = square_sum;
                if (shifted_square_sum > 0.) {
                  residual[j][k] -= autocovariance * autocovariance / shifted_square_sum;
                }
              }
            }
          }
        }
      }

      else {
        if (type[i] != REAL_VALUE) {
          sum = 0.;
          for (j = 0;j < nb_sequence;j++) {
            sum += int_sequence[j][i][time];
          }

          if (time + 1 < length[0]) {
            for (j = 0;j < nb_sequence;j++) {
              sum += int_sequence[j][i][time + 1];
            }
            mean = sum / (nb_sequence * 2);

            square_sum = 0.;
            shifted_square_sum = 0.;
            autocovariance = 0.;
            for (j = 0;j < nb_sequence;j++) {
              diff = int_sequence[j][i][time + 1] - mean;
              shifted_diff = int_sequence[j][i][time] - mean;
              square_sum += diff * diff;
              shifted_square_sum += shifted_diff * shifted_diff;
              autocovariance += diff * shifted_diff;
            }

            residual[0][time + 1] = square_sum;
            if (shifted_square_sum > 0.) {
              residual[0][time + 1] -= autocovariance * autocovariance / shifted_square_sum;
            }
          }

          for (j = time + 2;j < length[0];j++) {
            mean = sum / (nb_sequence * (j - time));
            square_sum_term[0] = 0.;
            square_sum_term[1] = 0.;
            square_sum_term[2] = 0.;
            diff = 0.;
            range_diff = 0.;

            for (k = 0;k < nb_sequence;k++) {
              sum += int_sequence[k][i][j];
              buff = int_sequence[k][i][j] - mean;
              diff += buff;
              square_sum_term[0] += buff * buff;
              shifted_diff = int_sequence[k][i][j - 1] - mean;
              square_sum_term[1] += shifted_diff * shifted_diff;
              square_sum_term[2] += buff * shifted_diff;
              range_diff += int_sequence[k][i][j] - int_sequence[k][i][time];
            }

            square_sum += square_sum_term[0] +
                          ((double)(j - time) / ((double)nb_sequence * (j - time + 1) * (j - time + 1))) * diff * diff -
                          (2. / ((double)nb_sequence * (j - time + 1))) * diff * range_diff;
            shifted_square_sum += square_sum_term[1] +
                                  ((double)(j - time) / ((double)nb_sequence * (j - time + 1) * (j - time + 1))) * diff * diff;
            autocovariance += square_sum_term[2] +
                              ((double)(j - time) / ((double)nb_sequence * (j - time + 1) * (j - time + 1))) * diff * diff -
                              (1. / ((double)nb_sequence * (j - time + 1))) * diff * range_diff;

            residual[0][j] = square_sum;
            if (shifted_square_sum > 0.) {
              residual[0][j] -= autocovariance * autocovariance / shifted_square_sum;
            }
          }
        }

        else {
          sum = 0.;
          for (j = 0;j < nb_sequence;j++) {
            sum += real_sequence[j][i][time];
          }

          if (time + 1 < length[0]) {
            for (j = 0;j < nb_sequence;j++) {
              sum += real_sequence[j][i][time + 1];
            }
            mean = sum / (nb_sequence * 2);

            square_sum = 0.;
            shifted_square_sum = 0.;
            autocovariance = 0.;
            for (j = 0;j < nb_sequence;j++) {
              diff = real_sequence[j][i][time + 1] - mean;
              shifted_diff = real_sequence[j][i][time] - mean;
              square_sum += diff * diff;
              shifted_square_sum += shifted_diff * shifted_diff;
              autocovariance += diff * shifted_diff;
            }

            residual[0][time + 1] = square_sum ;
            if (shifted_square_sum > 0.) {
              residual[0][time + 1] -= autocovariance * autocovariance / shifted_square_sum;
            }
          }

          for (j = time + 2;j < length[0];j++) {
            mean = sum / (nb_sequence * (j - time));
            square_sum_term[0] = 0.;
            square_sum_term[1] = 0.;
            square_sum_term[2] = 0.;
            diff = 0.;
            range_diff = 0.;

            for (k = 0;k < nb_sequence;k++) {
              sum += real_sequence[k][i][j];
              buff = real_sequence[k][i][j] - mean;
              diff += buff;
              square_sum_term[0] += buff * buff;
              shifted_diff = real_sequence[k][i][j - 1] - mean;
              square_sum_term[1] += shifted_diff * shifted_diff;
              square_sum_term[2] += buff * shifted_diff;
              range_diff += real_sequence[k][i][j] - real_sequence[k][i][time];
            }

            square_sum += square_sum_term[0] +
                          ((double)(j - time) / ((double)nb_sequence * (j - time + 1) * (j - time + 1))) * diff * diff -
                          (2. / ((double)nb_sequence * (j - time + 1))) * diff * range_diff;
            shifted_square_sum += square_sum_term[1] +
                                  ((double)(j - time) / ((double)nb_sequence * (j - time + 1) * (j - time + 1))) * diff * diff;
            autocovariance += square_sum_term[2] +
                              ((double)(j - time) / ((double)nb_sequence * (j - time + 1) * (j - time + 1))) * diff * diff -
                              (1. / ((double)nb_sequence * (j - time + 1))) * diff * range_diff;

            residual[0][j] = square_sum;
            if (shifted_square_sum > 0.) {
              residual[0][j] -= autocovariance * autocovariance / shifted_square_sum;
            }
          }
        }
      }
    }

    else if (model_type[i - 1] == STATIONARY_AUTOREGRESSIVE_MODEL_CHANGE) {
      if ((index != I_DEFAULT) || (!common_contrast)) {
        if (type[i] != REAL_VALUE) {
          for (j = 0;j < nb_sequence;j++) {
            if ((index == I_DEFAULT) || (index == j)) {
              square_sum = 0.;
              shifted_square_sum = 0.;
              autocovariance = 0.;

              for (k = time + 1;k < length[j];k++) {
                diff = int_sequence[j][i][k] - seq_mean[i][j];
                shifted_diff = int_sequence[j][i][k - 1] - seq_mean[i][j];
                square_sum += diff * diff;
                shifted_square_sum += shifted_diff * shifted_diff;
                autocovariance += diff * shifted_diff;

                residual[j][k] = square_sum;
                if (shifted_square_sum > 0.) {
                  residual[j][k] -= autocovariance * autocovariance / shifted_square_sum;
                }
              }
            }
          }
        }

        else {
          for (j = 0;j < nb_sequence;j++) {
            if ((index == I_DEFAULT) || (index == j)) {
              square_sum = 0.;
              shifted_square_sum = 0.;
              autocovariance = 0.;

              for (k = time + 1;k < length[j];k++) {
                diff = real_sequence[j][i][k] - seq_mean[i][j];
                shifted_diff = real_sequence[j][i][k - 1] - seq_mean[i][j];
                square_sum += diff * diff;
                shifted_square_sum += shifted_diff * shifted_diff;
                autocovariance += diff * shifted_diff;

                residual[j][k] = square_sum;
                if (shifted_square_sum > 0.) {
                  residual[j][k] -= autocovariance * autocovariance / shifted_square_sum;
                }
              }
            }
          }
        }
      }

      else {
        square_sum = 0.;
        shifted_square_sum = 0.;
        autocovariance = 0.;

        if (type[i] != REAL_VALUE) {
          for (j = time + 1;j < length[0];j++) {
            for (k = 0;k < nb_sequence;k++) {
              diff = int_sequence[k][i][j] - seq_mean[i][0];
              shifted_diff = int_sequence[k][i][j - 1] - seq_mean[i][0];
              square_sum += diff * diff;
              shifted_square_sum += shifted_diff * shifted_diff;
              autocovariance += diff * shifted_diff;
            }

            residual[0][j] = square_sum;
            if (shifted_square_sum > 0.) {
              residual[0][j] -= autocovariance * autocovariance / shifted_square_sum;
            }
          }
        }

        else {
          for (j = time + 1;j < length[0];j++) {
            for (k = 0;k < nb_sequence;k++) {
              diff = real_sequence[k][i][j] - seq_mean[i][0];
              shifted_diff = real_sequence[k][i][j - 1] - seq_mean[i][0];
              square_sum += diff * diff;
              shifted_square_sum += shifted_diff * shifted_diff;
              autocovariance += diff * shifted_diff;
            }

            residual[0][j] = square_sum;
            if (shifted_square_sum > 0.) {
              residual[0][j] -= autocovariance * autocovariance / shifted_square_sum;
            }
          }
        }
      }
    }

    else if (model_type[i - 1] == BAYESIAN_POISSON_CHANGE) {
      prior_contrast = -lgamma(hyperparam[i][0]) + hyperparam[i][0] * log(hyperparam[i][1]);

      sum = 0.;
      factorial_sum = 0.;
      for (j = time;j < length[index];j++) {
        sum += int_sequence[index][i][j];
        factorial_sum += factorial[i][index][j];
        if (contrast[j] != D_INF) {
          contrast[j] += prior_contrast - factorial_sum + lgamma(hyperparam[i][0] + sum) -
                         (hyperparam[i][0] + sum) * log(hyperparam[i][1] + j - time + 1);
        }
      }
    }

    else if (model_type[i - 1] == BAYESIAN_GAUSSIAN_CHANGE) {
      prior_contrast = log(hyperparam[i][1]) / 2 - lgamma(hyperparam[i][2] / 2) +
                       hyperparam[i][2] * log(hyperparam[i][3] / 2) / 2;

      if (type[i] != REAL_VALUE) {
        square_sum = 0.;
        sum = int_sequence[index][i][time];
        if (contrast[time] != D_INF) {
          diff = hyperparam[i][0] - sum;
          contrast[time] += prior_contrast - log(2 * M_PI) / 2 -
                            log(hyperparam[i][1] + 1) / 2 + lgamma((hyperparam[i][2] + 1) / 2) -
                            (hyperparam[i][2] + 1) *
                            log((hyperparam[i][3] + hyperparam[i][1] *
                                 diff * diff / (hyperparam[i][1] + 1)) / 2) / 2;
        }

        for (j = time + 1;j < length[index];j++) {
          diff = int_sequence[index][i][j] - sum / (j - time);
          square_sum += ((double)(j - time) / (double)(j - time + 1)) * diff * diff;
          sum += int_sequence[index][i][j];
          if (contrast[j] != D_INF) {
            diff = hyperparam[i][0] - sum / (j - time + 1);
            contrast[j] += prior_contrast - (j - time + 1) * log(2 * M_PI) / 2 -
                           log(hyperparam[i][1] + j - time + 1) / 2 +
                           lgamma((hyperparam[i][2] + j - time + 1) / 2) -
                           (hyperparam[i][2] + j - time + 1) *
                           logl((hyperparam[i][3] + square_sum + hyperparam[i][1] * (j - time + 1) *
                                 diff * diff / (hyperparam[i][1] + j - time + 1)) / 2) / 2;
          }
        }
      }

      else {
        square_sum = 0.;
        sum = real_sequence[index][i][time];
        if (contrast[time] != D_INF) {
          diff = hyperparam[i][0] - sum;
          contrast[time] += prior_contrast - log(2 * M_PI) / 2 -
                            log(hyperparam[i][1] + 1) / 2 + lgamma((hyperparam[i][2] + 1) / 2) -
                            (hyperparam[i][2] + 1) *
                            log((hyperparam[i][3] + hyperparam[i][1] *
                                 diff * diff / (hyperparam[i][1] + 1)) / 2) / 2;
        }

        for (j = time + 1;j < length[index];j++) {
          diff = real_sequence[index][i][j] - sum / (j - time);
          square_sum += ((double)(j - time) / (double)(j - time + 1)) * diff * diff;
          sum += real_sequence[index][i][j];
          if (contrast[j] != D_INF) {
            diff = hyperparam[i][0] - sum / (j - time + 1);
            contrast[j] += prior_contrast - (j - time + 1) * log(2 * M_PI) / 2 -
                           log(hyperparam[i][1] + j - time + 1) / 2 +
                           lgamma((hyperparam[i][2] + j - time + 1) / 2) -
                           (hyperparam[i][2] + j - time + 1) *
                           logl((hyperparam[i][3] + square_sum + hyperparam[i][1] * (j - time + 1) *
                                 diff * diff / (hyperparam[i][1] + j - time + 1)) / 2) / 2;
          }
        }
      }
    }

    if ((model_type[0] == MEAN_CHANGE) || (model_type[0] == INTERCEPT_SLOPE_CHANGE)) {
      if ((index != I_DEFAULT) || (!common_contrast)) {
        for (j = 0;j < nb_sequence;j++) {
          if ((index == I_DEFAULT) || (index == j)) {
            for (k = time + 1;k < length[j];k++) {
              contrast[k] -= residual[j][k];
            }
          }
        }
      }

      else {
        for (j = time;j < length[0];j++) {
          contrast[j] -= residual[0][j];
        }
      }
    }

    else if ((model_type[i - 1] == GAUSSIAN_CHANGE) || (model_type[i - 1] == VARIANCE_CHANGE) ||
             (model_type[i - 1] == ORDINAL_GAUSSIAN_CHANGE) || (model_type[i - 1] == LINEAR_MODEL_CHANGE)) {
      if ((index != I_DEFAULT) || (!common_contrast)) {
        for (j = 0;j < nb_sequence;j++) {
          if ((index == I_DEFAULT) || (index == j)) {
            for (k = time;k < length[j];k++) {
              if (contrast[k] != D_INF) {
//                if (residual[j][k] > 0.) {
                if (residual[j][k] > (k - time + 1) * ROUNDOFF_ERROR) {
                  contrast[k] -= ((double)(k - time + 1) / 2.) * (logl(residual[j][k] /
                                   (k - time + 1)) + log(2 * M_PI) + 1);
/*                  contrast[k] -= ((double)(k - time + 1) / 2.) * (logl(residual[j][k] /
                                   (k - time)) + log(2 * M_PI)) + (double)(k - time) / 2.; */
                }
                else {
                  contrast[k] = D_INF;
                }
              }
            }
          }
        }
      }

       else {
        for (j = time;j < length[0];j++) {
          if (contrast[j] != D_INF) {
//            if (residual[0][j] > 0.) {
            if (residual[0][j] > nb_sequence * (j - time + 1) * ROUNDOFF_ERROR) {
              contrast[j] -= ((double)(nb_sequence * (j - time + 1)) / 2.) * (logl(residual[0][j] /
                               (nb_sequence * (j - time + 1))) + log(2 * M_PI) + 1);
            }
            else {
              contrast[j] = D_INF;
            }
          }
        }
      }
    }

    else if ((model_type[i - 1] == AUTOREGRESSIVE_MODEL_CHANGE) ||
             (model_type[i - 1] == STATIONARY_AUTOREGRESSIVE_MODEL_CHANGE)) {
      if ((index != I_DEFAULT) || (!common_contrast)) {
        for (j = 0;j < nb_sequence;j++) {
          if ((index == I_DEFAULT) || (index == j)) {
            for (k = time + 1;k < length[j];k++) {
              if (contrast[k] != D_INF) {
//                if (residual[j][k] > 0.) {
                if (residual[j][k] > (k - time) * ROUNDOFF_ERROR) {
                  contrast[k] -= ((double)(k - time) / 2.) * (logl(residual[j][k] /
                                   (k - time)) + log(2 * M_PI) + 1);
                }
                else {
                  contrast[k] = D_INF;
                }
              }
            }
          }
        }
      }

      else {
        for (j = time + 1;j < length[0];j++) {
          if (contrast[j] != D_INF) {
//            if (residual[0][j] > 0.) {
            if (residual[0][j] > nb_sequence * (j - time) * ROUNDOFF_ERROR) {
              contrast[j] -= ((double)(nb_sequence * (j - time)) / 2.) * (logl(residual[0][j] /
                               (nb_sequence * (j - time))) + log(2 * M_PI) + 1);
            }
            else {
              contrast[j] = D_INF;
            }
          }
        }
      }
    }
  }

# ifdef DEBUG
  for (i = time;i < length[index == I_DEFAULT ? 0 : index];i++) {
    cout << contrast[i] << "  ";
  }
  cout << endl;
# endif

  delete [] frequency;
  delete [] inf_bound_parameter;

  if (residual) {
    if ((index != I_DEFAULT) || (!common_contrast)) {
      for (i = 0;i < nb_sequence;i++) {
        delete [] residual[i];
      }
    }
    else {
      delete [] residual[0];
    }
    delete [] residual;
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Optimal segmentation of a single sequence or a sample of sequences.
 *
 *  \param[in] index                    sequence index,
 *  \param[in] nb_segment               number of segments,
 *  \param[in] model_type               segment model types,
 *  \param[in] common_contrast          flag contrast functions common to the individuals,
 *  \param[in] shape_parameter          negative binomial shape parameters,
 *  \param[in] rank                     ranks (for ordinal variables),
 *  \param[in] isegmentation_likelihood pointer on the segmentation log-likelihoods,
 *  \param[in] nb_parameter             pointer on the number of free parameters of models,
 *  \param[in] segment_penalty          pointer on the penalties related to segment lengths (for mBIC).
 *
 *  \return                             log-likelihood of the optimal segmentation.
 */
/*--------------------------------------------------------------*/

double Sequences::segmentation(int index , int nb_segment , segment_model *model_type ,
                               bool common_contrast , double *shape_parameter , double **rank ,
                               double *isegmentation_likelihood , int *nb_parameter ,
                               double *segment_penalty)

{
  bool *used_output;
  int i , j , k , m , n , p;
  int max_nb_value , seq_length , count , *inf_bound_parameter , *seq_index_parameter ,
      *psegment , **optimal_length;
  double buff , segmentation_likelihood , **seq_mean , **hyperparam , **forward ,
         ***factorial , ***binomial_coeff;
  long double *contrast; 


  max_nb_value = 0;
  factorial = new double**[nb_variable];
  inf_bound_parameter = new int[nb_variable];
  binomial_coeff = new double**[nb_variable];
  seq_mean = new double*[nb_variable];
  seq_index_parameter = NULL;
  hyperparam = new double*[nb_variable];

  for (i = 1;i < nb_variable;i++) {
    if ((model_type[i - 1] == CATEGORICAL_CHANGE) && (marginal_distribution[i]->nb_value > max_nb_value)) {
      max_nb_value = marginal_distribution[i]->nb_value;
    }

    // computation of log of factorials for Poisson models

    if ((model_type[i - 1] == POISSON_CHANGE) || (model_type[i - 1] == BAYESIAN_POISSON_CHANGE)) {
      factorial[i] = new double*[nb_sequence];
      for (j = 0;j < nb_sequence;j++) {
        if ((index == I_DEFAULT) || (index == j)) {
          factorial[i][j] = new double[length[j]];
          for (k = 0;k < length[j];k++) {
            factorial[i][j][k] = log_factorial(int_sequence[j][i][k]);
          }
        }
        else {
          factorial[i][j] = NULL;
        }
      }
    }

    else {
      factorial[i] = NULL;
    }

    // computation of log of binomial coefficients for negative binomial models

    if ((model_type[i - 1] == NEGATIVE_BINOMIAL_0_CHANGE) || (model_type[i - 1] == NEGATIVE_BINOMIAL_1_CHANGE)) {
      switch (model_type[i - 1]) {
      case NEGATIVE_BINOMIAL_0_CHANGE :
        inf_bound_parameter[i - 1] = 0;
        break;
      case NEGATIVE_BINOMIAL_1_CHANGE :
        inf_bound_parameter[i - 1] = 1;
        break;
      }

      binomial_coeff[i] = new double*[nb_sequence];
      for (j = 0;j < nb_sequence;j++) {
        if ((index == I_DEFAULT) || (index == j)) {
          binomial_coeff[i][j] = new double[length[j]];
          for (k = 0;k < length[j];k++) {
            binomial_coeff[i][j][k] = log_binomial_coefficient(inf_bound_parameter[i - 1] , shape_parameter[i - 1] ,
                                                               int_sequence[j][i][k]);
          }
        }
        else {
          binomial_coeff[i][j] = NULL;
        }
      }
    }

    else {
      binomial_coeff[i] = NULL;
    }

    // computation of sequence means for Gaussian change in the variance models or
    // stationary piecewise autoregressive models

    if ((model_type[i - 1] == VARIANCE_CHANGE) || (model_type[i - 1] == STATIONARY_AUTOREGRESSIVE_MODEL_CHANGE)) {
      if ((index != I_DEFAULT) || (!common_contrast)) {
        seq_mean[i] = new double[nb_sequence];

        if (type[i] != REAL_VALUE) {
          for (j = 0;j < nb_sequence;j++) {
            if ((index == I_DEFAULT) || (index == j)) {
              seq_mean[i][j] = 0.;
              for (k = 0;k < length[j];k++) {
                seq_mean[i][j] += int_sequence[j][i][k];
              }
              seq_mean[i][j] /= length[j];
            }
          }
        }

        else {
          for (j = 0;j < nb_sequence;j++) {
            if ((index == I_DEFAULT) || (index == j)) {
              seq_mean[i][j] = 0.;
              for (k = 0;k < length[j];k++) {
                seq_mean[i][j] += real_sequence[j][i][k];
              }
              seq_mean[i][j] /= length[j];
            }
          }
        }
      }

      else {
        seq_mean[i] = new double[1];
        seq_mean[i][0] = 0.;

        if (type[i] != REAL_VALUE) {
          for (j = 0;j < length[0];j++) {
            for (k = 0;k < nb_sequence;k++) {
              seq_mean[i][0] += int_sequence[k][i][j];
            }
          }
        }

        else {
          for (j = 0;j < length[0];j++) {
            for (k = 0;k < nb_sequence;k++) {
              seq_mean[i][0] += real_sequence[k][i][j];
            }
          }
        }

        seq_mean[i][0] /= (nb_sequence * length[0]);
      }
    }

    else {
      seq_mean[i] = NULL;
    }

    if (((i == 1) && (model_type[0] == INTERCEPT_SLOPE_CHANGE)) ||
        ((model_type[i - 1] == LINEAR_MODEL_CHANGE) && (!seq_index_parameter))) {
      if (index_param_type == IMPLICIT_TYPE) {
        seq_index_parameter = new int[seq_length];
        for (j = 0;j < seq_length;j++) {
          seq_index_parameter[j] = j;
        }
      }
      else {
        seq_index_parameter = index_parameter[index == I_DEFAULT ? 0 : index];
      }
    }

    // computation of hyperparameters for Bayesian Poisson and Gaussian models

    if (model_type[i - 1] == BAYESIAN_POISSON_CHANGE) {
      hyperparam[i] = new double[2];
      gamma_hyperparameter_computation(index , i , hyperparam[i]);
    }
    else if (model_type[i - 1] == BAYESIAN_GAUSSIAN_CHANGE) {
      hyperparam[i] = new double[4];
      gaussian_gamma_hyperparameter_computation(index , i , hyperparam[i]);
    }
    else {
      hyperparam[i] = NULL;
    }
  }

  seq_length = length[index == I_DEFAULT ? 0 : index];
  contrast = new long double[seq_length];

  forward = new double*[seq_length];
  for (i = 0;i < seq_length;i++) {
    forward[i] = new double[nb_segment];
  }

  optimal_length = new int*[seq_length];
  for (i = 0;i < seq_length;i++) {
    optimal_length[i] = new int[nb_segment];
  }

  if ((nb_parameter) && (max_nb_value > 0)) {
    used_output = new bool[max_nb_value];
  }
  else {
    used_output = NULL;
  }

  // forward recurrence

  for (i = 0;i < seq_length;i++) {

    // computation of segment contrast functions (log-likelihoods or sum of squared deviations)

    forward_contrast(i , index , model_type , common_contrast , factorial ,
                     shape_parameter , binomial_coeff , seq_mean , seq_index_parameter ,
                     hyperparam , rank , contrast);

    for (j = 0;j < MIN((i < seq_length - 1 ? nb_segment - 1 : nb_segment) , i + 1);j++) {
//    for (j = MAX(0 , nb_segment + i - seq_length);j < MIN((i < seq_length - 1 ? nb_segment - 1 : nb_segment) , i + 1);j++) {
      if (j == 0) {
        forward[i][j] = contrast[0];
        if (forward[i][j] != D_INF) {
          optimal_length[i][j] = i + 1;
        }
      }

      else {
        forward[i][j] = D_INF;
        for (k = i;k >= j;k--) {
          if ((contrast[k] != D_INF) && (forward[k - 1][j - 1] != D_INF)) {
            buff = contrast[k] + forward[k - 1][j - 1];
            if (buff > forward[i][j]) {
              forward[i][j] = buff;
              optimal_length[i][j] = i - k + 1;
            }
          }
        }
      }
    }
  }

  if ((model_type[0] != MEAN_CHANGE) && (model_type[0] != INTERCEPT_SLOPE_CHANGE)) {
    if (isegmentation_likelihood) {
      for (i = 0;i < nb_segment;i++) {
        isegmentation_likelihood[i] = forward[seq_length - 1][i];
      }
    }

    segmentation_likelihood = forward[seq_length - 1][nb_segment - 1];
  }

  else {
    count = (index == I_DEFAULT ? nb_sequence : 1);

    if (isegmentation_likelihood) {
      for (i = 0;i < nb_segment;i++) {
        if (forward[seq_length - 1][i] < 0.) {
          isegmentation_likelihood[i] = -((double)(count * seq_length) / 2.) *
                                         (log(-forward[seq_length - 1][i] /
                                           (count * seq_length)) + log(2 * M_PI) + 1);
/*          isegmentation_likelihood[i] = -(((double)(count * seq_length) / 2.) *
                                          (log(-forward[seq_length - 1][i] /
                                            (count * (seq_length - nb_segment))) + log(2 * M_PI)) +
                                          (double)(count * (seq_length - nb_segment)) / 2.); */
        }
        else {
          isegmentation_likelihood[i] = D_INF;
        }
      }
    }

    if (forward[seq_length - 1][nb_segment - 1] < 0.) {
      segmentation_likelihood = -((double)(count * seq_length) / 2.) *
                                 (log(-forward[seq_length - 1][nb_segment - 1] /
                                   (count * seq_length)) + log(2 * M_PI) + 1);
/*      segmentation_likelihood = -(((double)(count * seq_length) / 2.) *
                                  (log(-forward[seq_length - 1][nb_segment - 1] /
                                    (count * (seq_length - nb_segment))) + log(2 * M_PI)) +
                                  (double)(count * (seq_length - nb_segment)) / 2.); */
    }
    else {
      segmentation_likelihood = D_INF;
    }
  }

  // computation of the penalty term related to the change-point distribution (modified BIC)

  if (segment_penalty) {

#   ifdef DEBUG
    int cumul_segment_length;
    cout << "\n";
#   endif

    for (i = 0;i < nb_segment;i++) {
      segment_penalty[i] = 0.;
      j = seq_length - 1;

#     ifdef DEBUG
      cumul_segment_length = 0;
#     endif

      for (k = i;k >= 0;k--) {

#       ifdef DEBUG
        cout << optimal_length[j][k] << " ";
        cumul_segment_length += optimal_length[j][k];
#       endif

        segment_penalty[i] += log((double)optimal_length[j][k]);
        j -= optimal_length[j][k];
      }

#     ifdef DEBUG
      cout << "| " << segment_penalty[i] << endl;
      if (cumul_segment_length != seq_length) {
        cout << "\nERROR: " << i << "   " << cumul_segment_length << " | " << seq_length << endl;
      }
#     endif

    }
  }

  // computation of the number of free parameters

  if (nb_parameter) {
    for (i = 0;i < nb_segment;i++) {
//      nb_parameter[i] = 0;
      nb_parameter[i] = i;

      if (model_type[0] == MEAN_CHANGE) {
        if ((index != I_DEFAULT) || (common_contrast)) {
          nb_parameter[i] += i + 2;
        }
        else {
          nb_parameter[i] += nb_sequence * (i + 1) + 1;
        }
      }

      else if (model_type[0] == INTERCEPT_SLOPE_CHANGE) {
        if ((index != I_DEFAULT) || (common_contrast)) {
          nb_parameter[i] += (i + 1) * 2 + 1;
        }
        else {
          nb_parameter[i] += nb_sequence * (i + 1) * 2 + 1;
        }
      }

      else {
        for (j = 1;j < nb_variable;j++) {
          if (model_type[j - 1] == CATEGORICAL_CHANGE) {
            if ((index != I_DEFAULT) || (!common_contrast)) {
              for (k = 0;k < nb_sequence;k++) {
                if ((index == I_DEFAULT) || (index == k)) {
                  m = length[k] - 1;

                  for (n = i;n >= 0;n--) {
                    for (p = 0;p < marginal_distribution[j]->nb_value;p++) {
                      used_output[p] = false;
                    }
                    nb_parameter[i]--;

                    for (p = m;p > m - optimal_length[m][n];p--) {
                      if (!used_output[int_sequence[k][j][p]]) {
                        nb_parameter[i]++;
                        used_output[int_sequence[k][j][p]] = true;
                      }
                    }

                    m -= optimal_length[m][n];
                  }
                }
              }
            }

            else {
              k = length[0] - 1;

              for (m = i;m >= 0;m--) {
                for (n = 0;n < marginal_distribution[j]->nb_value;n++) {
                  used_output[n] = false;
                }
                nb_parameter[i]--;

                for (n = k;n > k - optimal_length[k][m];n--) {
                  for (p = 0;p < nb_sequence;p++) {
                    if (!used_output[int_sequence[p][j][n]]) {
                      nb_parameter[i]++;
                      used_output[int_sequence[p][j][n]] = true;
                    }
                  }
                }

                k -= optimal_length[k][m];
              }
            }
          }

          else if ((model_type[j - 1] == POISSON_CHANGE) || (model_type[j - 1] == NEGATIVE_BINOMIAL_0_CHANGE) ||
                   (model_type[j - 1] == NEGATIVE_BINOMIAL_1_CHANGE) || (model_type[j - 1] == BAYESIAN_POISSON_CHANGE)) {
            if ((index != I_DEFAULT) || (common_contrast)) {
              nb_parameter[i] += i + 1;
            }
            else {
              nb_parameter[i] += nb_sequence * (i + 1);
            }
          }

          else if ((model_type[j - 1] == GAUSSIAN_CHANGE) || (model_type[j - 1] == ORDINAL_GAUSSIAN_CHANGE) ||
                   (model_type[j - 1] == BAYESIAN_GAUSSIAN_CHANGE)) {
            if ((index != I_DEFAULT) || (common_contrast)) {
              nb_parameter[i] += 2 * (i + 1);
            }
            else {
              nb_parameter[i] += nb_sequence * 2 * (i + 1);
            }
          }

          else if (model_type[j - 1] == VARIANCE_CHANGE) {
            if ((index != I_DEFAULT) || (common_contrast)) {
              nb_parameter[i] += i + 2;
            }
            else {
              nb_parameter[i] += nb_sequence * (i + 2);
            }
          }

          else if ((model_type[j - 1] == LINEAR_MODEL_CHANGE) || (model_type[j - 1] == AUTOREGRESSIVE_MODEL_CHANGE)) {
            if ((index != I_DEFAULT) || (common_contrast)) {
              nb_parameter[i] += 3 * (i + 1);
            }
            else {
              nb_parameter[i] += nb_sequence * 3 * (i + 1);
            }
          }

          else if (model_type[j - 1] == STATIONARY_AUTOREGRESSIVE_MODEL_CHANGE) {
            if ((index != I_DEFAULT) || (common_contrast)) {
              nb_parameter[i] += 2 * (i + 1) + 1;
            }
            else {
              nb_parameter[i] += nb_sequence * (2 * (i + 1) + 1);
            }
          }
        }
      }
    }
  }

  // restoration

  i = seq_length - 1;
  psegment = int_sequence[index == I_DEFAULT ? 0 : index][0] + i;

  for (j = nb_segment - 1;j >= 0;j--) {
//    for (k = 0;k < optimal_length[i][j];m++) {
    for (k = i;k > i - optimal_length[i][j];k--) {
      *psegment-- = j;
    }
    i -= optimal_length[i][j];
  }

  if (index == I_DEFAULT) {
    for (i = 1;i < nb_sequence;i++) {
      for (j = 0;j < length[0];j++) {
        int_sequence[i][0][j] = int_sequence[0][0][j];
      }
    }
  }

  min_value[0] = 0;
  max_value[0] = nb_segment - 1;
  delete marginal_distribution[0];
  build_marginal_frequency_distribution(0);

  for (i = 1;i < nb_variable;i++) {
    if ((model_type[i - 1] == POISSON_CHANGE) || (model_type[i - 1] == BAYESIAN_POISSON_CHANGE)) {
      for (j = 0;j < nb_sequence;j++) {
        delete [] factorial[i][j];
      }
      delete [] factorial[i];
    }

    if ((model_type[i - 1] == NEGATIVE_BINOMIAL_0_CHANGE) || (model_type[i - 1] == NEGATIVE_BINOMIAL_1_CHANGE)) {
      for (j = 0;j < nb_sequence;j++) {
        delete [] binomial_coeff[i][j];
      }
      delete [] binomial_coeff[i];
    }

    delete [] seq_mean[i];
    delete [] hyperparam[i];
  }
  delete [] factorial;
  delete [] inf_bound_parameter;
  delete [] binomial_coeff;
  delete [] seq_mean;
  delete [] hyperparam;

  if (index_param_type == IMPLICIT_TYPE) {
    delete [] seq_index_parameter;
  }

  delete [] contrast;

  for (i = 0;i < seq_length;i++) {
    delete [] forward[i];
  }
  delete [] forward;

  for (i = 0;i < seq_length;i++) {
    delete [] optimal_length[i];
  }
  delete [] optimal_length;

  delete [] used_output;

  return segmentation_likelihood;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Optimal segmentation of a single sequence or a sample of sequences.
 *
 *  \param[in] error           reference on a StatError object,
 *  \param[in] display         flag for displaying the segmentation,
 *  \param[in] iidentifier     sequence identifier,
 *  \param[in] nb_segment      number of segments,
 *  \param[in] model_type      segment model types,
 *  \param[in] common_contrast flag contrast functions common to the individuals,
 *  \param[in] shape_parameter negative binomial shape parameters,
 *  \param[in] output          output (sequence or residuals).
 *  \param[in] continuity      flag continuous piecewise linear function.
 *
 *  \return                    Sequences object.
 */
/*--------------------------------------------------------------*/

Sequences* Sequences::segmentation(StatError &error , bool display , int iidentifier ,
                                   int nb_segment , segment_model *model_type ,
                                   bool common_contrast , double *shape_parameter ,
                                   sequence_type output , bool continuity) const

{
  bool status = true;
  int i , j;
  int index , nb_parameter;
  double segmentation_likelihood , segment_penalty , penalized_likelihood , **rank;
  FrequencyDistribution *marginal;
  Sequences *seq , *iseq , *oseq;


  oseq = NULL;
  error.init();

/*  if (((index_param_type == TIME) && (index_interval->variance > 0.)) ||
      (index_param_type == POSITION)) {
    status = false;
    error.update(SEQ_error[SEQR_INDEX_PARAMETER_TYPE]);
  }
  if (index_param_type == POSITION) {
    status = false;
    error.correction_update(SEQ_error[SEQR_INDEX_PARAMETER_TYPE] , SEQ_index_parameter_word[TIME]);
  } */

  for (i = 0;i < nb_variable;i++) {
    if ((model_type[i] == CATEGORICAL_CHANGE) || (model_type[i] == POISSON_CHANGE) ||
        (model_type[i] == NEGATIVE_BINOMIAL_0_CHANGE) || (model_type[i] == NEGATIVE_BINOMIAL_1_CHANGE) ||
        (model_type[i] == ORDINAL_GAUSSIAN_CHANGE) || (model_type[i] == BAYESIAN_POISSON_CHANGE)) {
      if ((type[i] != INT_VALUE) && (type[i] != STATE)) {
        status = false;
        ostringstream error_message , correction_message;
        error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                      << STAT_error[STATR_VARIABLE_TYPE];
        correction_message << STAT_variable_word[INT_VALUE] << " or "
                           << STAT_variable_word[STATE];
        error.correction_update((error_message.str()).c_str() , (correction_message.str()).c_str());
      }

      else {
        if (((model_type[i] != NEGATIVE_BINOMIAL_1_CHANGE) && (min_value[i] < 0)) ||
            ((model_type[i] == NEGATIVE_BINOMIAL_1_CHANGE) && (min_value[i] < 1))) {
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

        else if (model_type[i] == CATEGORICAL_CHANGE) {
          if ((marginal_distribution[i]->nb_value < 2) ||
              (marginal_distribution[i]->nb_value > NB_OUTPUT)) {
            status = false;
            ostringstream error_message;
            error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                          << STAT_error[STATR_NB_VALUE];
            error.update((error_message.str()).c_str());
          }

          else {
            for (j = 0;j < marginal_distribution[i]->nb_value;j++) {
              if (marginal_distribution[i]->frequency[j] == 0) {
                status = false;
                ostringstream error_message;
                error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                              << STAT_error[STATR_MISSING_VALUE] << " " << j;
                error.update((error_message.str()).c_str());
              }
            }
          }
        }
      }

      if (((model_type[i] == CATEGORICAL_CHANGE) || (model_type[i] == ORDINAL_GAUSSIAN_CHANGE)) &&
          ((output == SUBTRACTION_RESIDUAL) || (output == ABSOLUTE_RESIDUAL) || (output == DIVISION_RESIDUAL))) {
        status = false;
        error.update(SEQ_error[SEQR_FORBIDDEN_OUTPUT]);
      }
    }

    else if ((type[i] != INT_VALUE) && (type[i] != STATE) && (type[i] != REAL_VALUE)) {
      status = false;
      ostringstream error_message , correction_message;
      error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                    << STAT_error[STATR_VARIABLE_TYPE];
      correction_message << STAT_variable_word[INT_VALUE] << " or "
                         << STAT_variable_word[STATE] << " or "
                         << STAT_variable_word[REAL_VALUE];
      error.correction_update((error_message.str()).c_str() , (correction_message.str()).c_str());
    }

    else if (((model_type[i] == AUTOREGRESSIVE_MODEL_CHANGE) || (model_type[i] == STATIONARY_AUTOREGRESSIVE_MODEL_CHANGE)) &&
             (index_param_type != IMPLICIT_TYPE) && (index_interval->variance > 0.)) {
      status = false;
      error.update(SEQ_error[SEQR_INDEX_PARAMETER_TYPE]);
    }

    if (((model_type[i] == CATEGORICAL_CHANGE) || (model_type[i] == ORDINAL_GAUSSIAN_CHANGE) ||
         (model_type[i] == AUTOREGRESSIVE_MODEL_CHANGE) || (model_type[i] == STATIONARY_AUTOREGRESSIVE_MODEL_CHANGE)) &&
        (output == SEQUENCE_SAMPLE)) {
      status = false;
      error.update(SEQ_error[SEQR_FORBIDDEN_OUTPUT]);
    }
  }

  if (iidentifier != I_DEFAULT) {
    for (i = 0;i < nb_sequence;i++) {
      if (iidentifier == identifier[i]) {
        index = i;
        break;
      }
    }

    if (i == nb_sequence) {
      status = false;
      error.update(SEQ_error[SEQR_SEQUENCE_IDENTIFIER]);
    }
  }

  else {
    index = I_DEFAULT;
    if (length_distribution->variance > 0.) {
      status = false;
      error.update(SEQ_error[SEQR_VARIABLE_SEQUENCE_LENGTH]);
    }
  }

  if (((index != I_DEFAULT) || (!common_contrast)) && (output == SEQUENCE_SAMPLE)) {
    status = false;
    error.update(SEQ_error[SEQR_FORBIDDEN_OUTPUT]);
  }

  if ((status) && ((nb_segment < 1) || (nb_segment > length[index == I_DEFAULT ? 0 : index] / 2))) {
    status = false;
    error.update(SEQ_error[SEQR_NB_SEGMENT]);
  }

  if (status) {
    if (index != I_DEFAULT) {
      iseq = new Sequences(*this , 1 , &index);
      seq = new Sequences(*iseq , ADD_STATE_VARIABLE);
      delete iseq;
    }
    else {
      seq = new Sequences(*this , ADD_STATE_VARIABLE);
    }

    // rank computation for ordinal variables

    rank = new double*[seq->nb_variable];

    for (i = 1;i < seq->nb_variable;i++) {
      if (model_type[i - 1] == ORDINAL_GAUSSIAN_CHANGE) {
        rank[i] = seq->marginal_distribution[i]->rank_computation();
      }
      else {
        rank[i] = NULL;
      }
    }

    segmentation_likelihood = seq->segmentation((index == I_DEFAULT ? index : 0) , nb_segment , model_type ,
                                                common_contrast , shape_parameter , rank);

    for (i = 1;i < seq->nb_variable;i++) {
      delete [] rank[i];
    }
    delete [] rank;

    if (segmentation_likelihood != D_INF) {
      if (display) {
        segment_penalty = 0.;
        i = 0;
        for (j = 1;j < seq->length[0];j++) {
          if (seq->int_sequence[0][0][j] != seq->int_sequence[0][0][j - 1]) {
            segment_penalty += log((double)(j - i));
            i = j;
          }
        }
        segment_penalty += log((double)(seq->length[0] - i));

        nb_parameter = seq->nb_parameter_computation((index == I_DEFAULT ? index : 0) , nb_segment , model_type ,
                                                     common_contrast);

        penalized_likelihood = 2 * segmentation_likelihood - nb_parameter *
                               log((double)((seq->nb_variable - 1) * seq->length[0])) - segment_penalty;

        cout << "\n" << nb_segment << " " << (nb_segment == 1 ? SEQ_label[SEQL_SEGMENT] : SEQ_label[SEQL_SEGMENTS])
             << "   2 * " << STAT_label[STATL_LIKELIHOOD] << ": " << 2 * segmentation_likelihood << "   "
             << nb_parameter << " " << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS]
             << "   2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " (Modified "  << STAT_criterion_word[BIC] << "): "
             << penalized_likelihood << endl;
      }

      oseq = seq->segmentation_output(nb_segment , model_type , common_contrast , display , output ,
                                      NULL , continuity);

      if ((output == SEQUENCE) || (output == ABSOLUTE_RESIDUAL)) {
        delete seq;
      }
    }

    else {
      delete seq;
      oseq = NULL;
      error.update(SEQ_error[SEQR_SEGMENTATION_FAILURE]);
    }
  }

  return oseq;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Optimal segmentation of a single sequence or a sample of sequences.
 *
 *  \param[in] error           reference on a StatError object,
 *  \param[in] display         flag for displaying the segmentation,
 *  \param[in] iidentifier     sequence identifier,
 *  \param[in] nb_segment      number of segments,
 *  \param[in] model_type      segment model types,
 *  \param[in] common_contrast flag contrast functions common to the individuals,
 *  \param[in] shape_parameter negative binomial shape parameters,
 *  \param[in] output          output (sequence or residuals).
 *  \param[in] continuity      flag continuous piecewise linear function.
 *
 *  \return                    Sequences object.
 */
/*--------------------------------------------------------------*/

Sequences* Sequences::segmentation(StatError &error , bool display , int iidentifier ,
                                   int nb_segment , vector<segment_model> model_type ,
                                   bool common_contrast , vector<double> shape_parameter ,
                                   sequence_type output , bool continuity) const

{
  return segmentation(error , display , iidentifier , nb_segment , model_type.data() ,
                      common_contrast , shape_parameter.data() , output , continuity);
}


};  // namespace sequence_analysis
