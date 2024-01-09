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
 *       $Id: change_points3.cpp 11914 2012-03-26 06:29:13Z guedon $
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

#include <boost/math/special_functions/gamma.hpp>

#include "tool/config.h"

#include "stat_tool/stat_tools.h"
#include "stat_tool/curves.h"
#include "stat_tool/distribution.h"
#include "stat_tool/markovian.h"
#include "stat_tool/vectors.h"
#include "stat_tool/distance_matrix.h"

#include "sequences.h"
#include "sequence_label.h"

using namespace std;
using namespace boost::math;
using namespace stat_tool;


namespace sequence_analysis {


#if defined (SYSTEM_IS__CYGWIN)
#define expl exp
#endif



/*--------------------------------------------------------------*
 *
 *  Calcul des N segmentations les plus probables d'une sequence.
 *
 *  arguments : indice de la sequence, nombre de segments, types des modeles,
 *              rangs (variables ordinales), stream, format de fichier ('a' : ASCII,
 *              's' : Spreadsheet), nombre de segmentation, vraisemblance des donnees.
 *
 *--------------------------------------------------------------*/

double Sequences::N_segmentation(int index , int nb_segment , int *model_type ,
                                 double **irank , ostream &os , char format ,
                                 int inb_segmentation , double likelihood) const

{
  bool **active_cell;
  register int i , j , k , m , n;
  int max_nb_value , brank , previous_rank , nb_cell , *frequency , *seq_index_parameter = NULL ,
      *rank , *psegment , ***optimal_length , ***optimal_rank;
  double sum , factorial_sum , proba , diff , index_parameter_sum , index_parameter_diff , buff ,
         segmentation_likelihood , response_mean , index_parameter_mean , index_parameter_variance ,
         covariance , residual_mean , *nb_segmentation , *sequence_mean , **hyperparam , **factorial ,
         **nb_segmentation_forward , **mean , **variance , **intercept , **slope , ***forward;
  long double square_sum , index_parameter_square_sum , mix_square_sum , prior_contrast , *residual ,
              *contrast , likelihood_cumul , posterior_probability_cumul;

# ifdef MESSAGE
  int *change_point;
  double sum2;
  long double norm;
# endif


  // initialisations

  max_nb_value = 0;
  factorial = new double*[nb_variable];
  hyperparam = new double*[nb_variable];

  for (i = 1;i < nb_variable;i++) {
    if (((model_type[i - 1] == CATEGORICAL_CHANGE) || (model_type[0] == MULTIVARIATE_CATEGORICAL_CHANGE)) &&
        (marginal_distribution[i]->nb_value > max_nb_value)) {
      max_nb_value = marginal_distribution[i]->nb_value;
    }

    if ((model_type[i - 1] == POISSON_CHANGE) || (model_type[0] == MULTIVARIATE_POISSON_CHANGE) ||
        (model_type[i - 1] == BAYESIAN_POISSON_CHANGE)) {
      factorial[i] = new double[length[index]];
    }
    else {
      factorial[i] = NULL;
    }

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

  if (max_nb_value > 0) {
    frequency = new int[max_nb_value];
  }
  else {
    frequency = NULL;
  }

  sequence_mean = new double[nb_variable];
  residual = new long double[length[index]];

  contrast = new long double[length[index]];

  nb_segmentation_forward = new double*[length[index]];
  for (i = 0;i < length[index];i++) {
    nb_segmentation_forward[i] = new double[nb_segment];
  }

  forward = new double**[length[index]];
  for (i = 0;i < length[index];i++) {
    forward[i] = new double*[nb_segment];
    for (j = 0;j < nb_segment;j++) {
      forward[i][j] = new double[inb_segmentation];
    }
  }

  nb_segmentation = new double[nb_segment];
  rank = new int[length[index] + 1];

  optimal_length = new int**[length[index]];
  for (i = 0;i < length[index];i++) {
    optimal_length[i] = new int*[nb_segment];
    for (j = 0;j < nb_segment;j++) {
      optimal_length[i][j] = new int[inb_segmentation];
    }
  }

  optimal_rank = new int**[length[index]];
  for (i = 0;i < length[index];i++) {
    optimal_rank[i] = new int*[nb_segment];
    for (j = 0;j < nb_segment;j++) {
      optimal_rank[i][j] = new int[inb_segmentation];
    }
  }

  mean = new double*[nb_variable];
  variance = new double*[nb_variable];
  intercept = new double*[nb_variable];
  slope = new double*[nb_variable];

  for (i = 1;i < nb_variable;i++) {
    if ((model_type[i - 1] == POISSON_CHANGE) || (model_type[0] == MULTIVARIATE_POISSON_CHANGE) ||
        (model_type[i - 1] == GAUSSIAN_CHANGE) || (model_type[i - 1] == VARIANCE_CHANGE) ||
        (model_type[0] == MEAN_CHANGE) || (model_type[0] == MEAN_VARIANCE_CHANGE) ||
        (model_type[i - 1] == BAYESIAN_POISSON_CHANGE) || (model_type[i - 1] == BAYESIAN_GAUSSIAN_CHANGE)) {
      mean[i] = new double[nb_segment];
      variance[i] = new double[nb_segment];
    }
    else if ((model_type[i - 1] == GEOMETRIC_0_CHANGE) || (model_type[i - 1] == GEOMETRIC_1_CHANGE) ||
             (model_type[0] == MULTIVARIATE_GEOMETRIC_0_CHANGE)) {
      mean[i] = new double[nb_segment];
      variance[i] = NULL;
    }
    else if ((model_type[i - 1] == LINEAR_MODEL_CHANGE) || (model_type[0] == INTERCEPT_SLOPE_CHANGE)) {
      mean[i] = NULL;
      variance[i] = new double[nb_segment];
    }
    else {
      mean[i] = NULL;
      variance[i] = NULL;
    }

    if ((model_type[i - 1] == LINEAR_MODEL_CHANGE) || (model_type[0] == INTERCEPT_SLOPE_CHANGE)) {
      intercept[i] = new double[nb_segment];
      slope[i] = new double[nb_segment];
    }
    else {
      intercept[i] = NULL;
      slope[i] = NULL;
    }
  }

  active_cell = new bool*[length[index]];
  for (i = 0;i < length[index];i++) {
    active_cell[i] = new bool[nb_segment];
    for (j = 0;j < nb_segment;j++) {
      active_cell[i][j] = false;
    }
  }

# ifdef MESSAGE
  double **mean_square_diff;

  mean_square_diff = new double*[nb_variable];
  for (i = 1;i < nb_variable;i++) {
    if ((model_type[i - 1] == GAUSSIAN_CHANGE) || (model_type[0] == MEAN_CHANGE) ||
        (model_type[0] == MEAN_VARIANCE_CHANGE)) {
      mean_square_diff[i] = new double[length[index]];
    }
    else {
      mean_square_diff[i] = NULL;
    }
  }
# endif

# ifdef DEBUG
  double **segment_probability;

  segment_probability = new double*[length[index]];
  for (i = 0;i < length[index];i++) {
    segment_probability[i] = new double[nb_segment];
    for (j = 0;j < nb_segment;j++) {
      segment_probability[i][j] = D_INF;
    }
  }
# endif

# ifdef MESSAGE
  double **smoothed_probability;
  long double approximated_likelihood = 0.;

  if (inb_segmentation >= 1000) {
    smoothed_probability = new double*[length[index]];
    for (i = 0;i < length[index];i++) {
      smoothed_probability[i] = new double[nb_segment];
      for (j = 0;j < nb_segment;j++) {
        smoothed_probability[i][j] = 0.;
      }
    }
  }
# endif

  for (i = 1;i < nb_variable;i++) {
    if (model_type[i - 1] == VARIANCE_CHANGE) {
      sequence_mean[i] = 0.;
      if (type[i] != REAL_VALUE) {
        for (j = 0;j < length[index];j++) {
          sequence_mean[i] += int_sequence[index][i][j];
        }
      }
      else {
        for (j = 0;j < length[index];j++) {
          sequence_mean[i] += real_sequence[index][i][j];
        }
      }
      sequence_mean[i] /= length[index];
    }

    if (((model_type[i - 1] == LINEAR_MODEL_CHANGE) && (!seq_index_parameter)) ||
        ((i == 1) && (model_type[0] == INTERCEPT_SLOPE_CHANGE))) {
      if (index_parameter_type == IMPLICIT_TYPE) {
        seq_index_parameter = new int[length[index]];
        for (j = 0;j < length[index];j++) {
          seq_index_parameter[j] = j;
        }
      }
      else {
        seq_index_parameter = index_parameter[index];
      }
    }
  }

# ifdef DEBUG
  for (i = 0;i < nb_segment;i++) {
    nb_segmentation[i] = 1;
  }
# endif

  // recurrence "forward"

  for (i = 0;i < length[index];i++) {

    // calcul des log-vraisemblances des segments

    if (model_type[0] == MULTIVARIATE_CATEGORICAL_CHANGE) {
      for (j = 0;j < max_nb_value;j++) {
        frequency[j] = 0;
      }

      for (j = i;j >= 0;j--) {
        for (k = 1;k < nb_variable;k++) {
          frequency[int_sequence[index][k][j]]++;
        }

        contrast[j] = 0.;
        for (k = 0;k < max_nb_value;k++) {
          if (frequency[k] > 0) {
            contrast[j] += frequency[k] * log((double)frequency[k] / (double)((nb_variable - 1) * (i - j + 1)));
          }
        }
      }
    }

    else if (model_type[0] == MULTIVARIATE_POISSON_CHANGE) {
      for (j = 1;j < nb_variable;j++) {
        factorial[j][i] = 0.;
        for (k = 2;k <= int_sequence[index][j][i];k++) {
          factorial[j][i] += log((double)k);
        }
      }

      sum = 0.;
      factorial_sum = 0.;
      for (j = i;j >= 0;j--) {
        for (k = 1;k < nb_variable;k++) {
          sum += int_sequence[index][k][j];
          factorial_sum += factorial[k][j];
        }
        if (sum > 0.) {
          contrast[j] = sum * (log(sum / ((nb_variable - 1) * (i - j + 1))) - 1) - factorial_sum;
        }
        else {
          contrast[j] = 0.;
        }
      }
    }

    else if (model_type[0] == MULTIVARIATE_GEOMETRIC_0_CHANGE) {
      sum = 0.;
      for (j = i;j >= 0;j--) {
        for (k = 1;k < nb_variable;k++) {
          sum += int_sequence[index][k][j];
        }
        if (sum > 0.) {
          proba = (nb_variable - 1) * (i - j + 1) / ((nb_variable - 1) * (i - j + 1) + sum);
          contrast[j] = (nb_variable - 1) * (i - j + 1) * log(proba) + sum * log(1. - proba);
        }
        else {
          contrast[j] = 0.;
        }
      }
    }

    else {
      for (j = 0;j <= i;j++) {
        contrast[j] = 0.;
      }

      for (j = 1;j < nb_variable;j++) {
        if (model_type[j - 1] == CATEGORICAL_CHANGE) {
          for (k = 0;k < marginal_distribution[j]->nb_value;k++) {
            frequency[k] = 0;
          }
          sum = 0.;

          frequency[int_sequence[index][j][i]]++;
          for (k = i - 1;k >= 0;k--) {
            sum += (i - k) * log((double)(i - k) / (double)(i - k + 1)) +
                   log((double)(frequency[int_sequence[index][j][k]] + 1) / (double)(i - k + 1));
            if (frequency[int_sequence[index][j][k]] > 0) {
              sum -= frequency[int_sequence[index][j][k]] *
                     log((double)frequency[int_sequence[index][j][k]] / (double)(frequency[int_sequence[index][j][k]] + 1));
            }
            frequency[int_sequence[index][j][k]]++;

            if (contrast[k] != D_INF) {
              contrast[k] += sum;
            }

#           ifdef MESSAGE
            sum2 = 0.;
            for (m = 0;m < marginal_distribution[j]->nb_value;m++) {
              if (frequency[m] > 0) {
                sum2 += frequency[m] * log((double)frequency[m] / (double)(i - k + 1));
              }
            }

            if ((sum2 < sum - DOUBLE_ERROR) || (sum2 > sum + DOUBLE_ERROR)) {
              cout << "\nERROR: " << i << " " << k << " | " << sum2 << " " << sum << endl;
            }
#           endif

          }
        }

        else if (model_type[j - 1] == POISSON_CHANGE) {
          factorial[j][i] = 0.;
          for (k = 2;k <= int_sequence[index][j][i];k++) {
            factorial[j][i] += log((double)k);
          }

          sum = 0.;
          factorial_sum = 0.;
          for (k = i;k >= 0;k--) {
            sum += int_sequence[index][j][k];
            factorial_sum += factorial[j][k];
            if ((contrast[k] != D_INF) && (sum > 0.)) {
              contrast[k] += sum * (log(sum / (i - k + 1)) - 1) - factorial_sum;
            }
          }
        }

        else if (model_type[j - 1] == GEOMETRIC_0_CHANGE) {
          sum = 0.;
          for (k = i;k >= 0;k--) {
            sum += int_sequence[index][j][k];
            if ((contrast[k] != D_INF) && (sum > 0.)) {
              proba = (i - k + 1) / (i - k + 1 + sum);
              contrast[k] += (i - k + 1) * log(proba) + sum * log(1. - proba);
            }
          }
        }

        else if (model_type[j - 1] == GEOMETRIC_1_CHANGE) {
          sum = 0.;
          for (k = i;k >= 0;k--) {
            sum += int_sequence[index][j][k];
            if ((contrast[k] != D_INF) && (sum > i - k + 1)) {
              proba = (i - k + 1) / sum;
              contrast[k] += (i - k + 1) * log(proba) + (sum - (i - k + 1)) * log(1. - proba);
            }
          }
        }

        else if (model_type[j - 1] == BAYESIAN_POISSON_CHANGE) {
          prior_contrast = -lgamma(hyperparam[j][0]) + hyperparam[j][0] * log(hyperparam[j][1]);

          factorial[j][i] = 0.;
          for (k = 2;k <= int_sequence[index][j][i];k++) {
            factorial[j][i] += log((double)k);
          }

          sum = 0.;
          factorial_sum = 0.;
          for (k = i;k >= 0;k--) {
            sum += int_sequence[index][j][k];
            factorial_sum += factorial[j][k];
            if (contrast[k] != D_INF) {
              contrast[k] += prior_contrast - factorial_sum + lgamma(hyperparam[j][0] + sum) -
                             (hyperparam[j][0] + sum) * log(hyperparam[j][1] + i - k + 1);
            }
          }
        }

        else if (model_type[j - 1] == BAYESIAN_GAUSSIAN_CHANGE) {
          prior_contrast = log(hyperparam[j][1]) / 2 - lgamma(hyperparam[j][2] / 2) +
                           hyperparam[j][2] * log(hyperparam[j][3] / 2) / 2;

          if (type[j] != REAL_VALUE) {
            square_sum = 0.;
            sum = int_sequence[index][j][i];
            if (contrast[i] != D_INF) {
              diff = hyperparam[j][0] - sum;
              contrast[i] += prior_contrast - log(2 * M_PI) / 2 -
                             log(hyperparam[j][1] + 1) / 2 + lgamma((hyperparam[j][2] + 1) / 2) -
                             (hyperparam[j][2] + 1) *
                             log((hyperparam[j][3] + hyperparam[j][1] *
                                  diff * diff / (hyperparam[j][1] + 1)) / 2) / 2;
            }

            for (k = i - 1;k >= 0;k--) {
              diff = int_sequence[index][j][k] - sum / (i - k);
              square_sum += ((double)(i - k) / (double)(i - k + 1)) * diff * diff;
              sum += int_sequence[index][j][k];
              if (contrast[k] != D_INF) {
                diff = hyperparam[j][0] - sum / (i - k + 1);
                contrast[k] += prior_contrast - (i - k + 1) * log(2 * M_PI) / 2 -
                               log(hyperparam[j][1] + i - k + 1) / 2 +
                               lgamma((hyperparam[j][2] + i - k + 1) / 2) -
                               (hyperparam[j][2] + i - k + 1) *
                               logl((hyperparam[j][3] + square_sum + hyperparam[j][1] * (i - k + 1) *
                                     diff * diff / (hyperparam[j][1] + i - k + 1)) / 2) / 2;
              }
            }
          }

          else {
            square_sum = 0.;
            sum = real_sequence[index][j][i];
            if (contrast[i] != D_INF) {
              diff = hyperparam[j][0] - sum;
              contrast[i] += prior_contrast - log(2 * M_PI) / 2 -
                             log(hyperparam[j][1] + 1) / 2 + lgamma((hyperparam[j][2] + 1) / 2) -
                             (hyperparam[j][2] + 1) *
                             log((hyperparam[j][3] + hyperparam[j][1] *
                                  diff * diff / (hyperparam[j][1] + 1)) / 2) / 2;
            }

            for (k = i - 1;k >= 0;k--) {
              diff = real_sequence[index][j][k] - sum / (i - k);
              square_sum += ((double)(i - k) / (double)(i - k + 1)) * diff * diff;
              sum += real_sequence[index][j][k];
              if (contrast[k] != D_INF) {
                diff = hyperparam[j][0] - sum / (i - k + 1);
                contrast[k] += prior_contrast - (i - k + 1) * log(2 * M_PI) / 2 -
                               log(hyperparam[j][1] + i - k + 1) / 2 +
                               lgamma((hyperparam[j][2] + i - k + 1) / 2) -
                               (hyperparam[j][2] + i - k + 1) *
                               logl((hyperparam[j][3] + square_sum + hyperparam[j][1] * (i - k + 1) *
                                     diff * diff / (hyperparam[j][1] + i - k + 1)) / 2) / 2;
              }
            }
          }
        }

        else {
          if (model_type[j - 1] == VARIANCE_CHANGE) {
            square_sum = 0.;

            if (type[j] != REAL_VALUE) {
              for (k = i;k >= 0;k--) {
                diff = int_sequence[index][j][k] - sequence_mean[j];
                square_sum += diff * diff;
                residual[k] = square_sum;
              }
            }

            else {
              for (k = i;k >= 0;k--) {
                diff = real_sequence[index][j][k] - sequence_mean[j];
                square_sum += diff * diff;
                residual[k] = square_sum;
              }
            }
          }

          else if (model_type[j - 1] == ORDINAL_GAUSSIAN_CHANGE) {
            square_sum = 0.;
            sum = irank[j][int_sequence[index][j][i]];
            residual[i] = 0.;

            for (k = i - 1;k >= 0;k--) {
              diff = irank[j][int_sequence[index][j][k]] - sum / (i - k);
              square_sum += ((double)(i - k) / (double)(i - k + 1)) * diff * diff;
              sum += irank[j][int_sequence[index][j][k]];
              residual[k] = square_sum;
            }
          }

          else if ((model_type[j - 1] == LINEAR_MODEL_CHANGE) || (model_type[0] == INTERCEPT_SLOPE_CHANGE)) {
            if (type[j] != REAL_VALUE) {
              square_sum = 0.;
              index_parameter_square_sum = 0.;
              mix_square_sum = 0.;
              sum = int_sequence[index][j][i];
              index_parameter_sum = seq_index_parameter[i];
              residual[i] = 0.;

              for (k = i - 1;k >= 0;k--) {
                diff = int_sequence[index][j][k] - sum / (i - k);
                square_sum += ((double)(i - k) / (double)(i - k + 1)) * diff * diff;
                index_parameter_diff = seq_index_parameter[k] - index_parameter_sum / (i - k);
                index_parameter_square_sum += ((double)(i - k) / (double)(i - k + 1)) *
                                              index_parameter_diff * index_parameter_diff;
                mix_square_sum += ((double)(i - k) / (double)(i - k + 1)) * diff * index_parameter_diff;
                sum += int_sequence[index][j][k];
                index_parameter_sum += seq_index_parameter[k];

                if ((k < i - 1) && (index_parameter_square_sum > 0.)) {
                  residual[k] = square_sum - mix_square_sum * mix_square_sum / index_parameter_square_sum;
                }
                else {
                  residual[k] = 0.;
                }
              }
            }

            else {
              square_sum = 0.;
              index_parameter_square_sum = 0.;
              mix_square_sum = 0.;
              sum = real_sequence[index][j][i];
              index_parameter_sum = seq_index_parameter[i];
              residual[i] = 0.;

              for (k = i - 1;k >= 0;k--) {
                diff = real_sequence[index][j][k] - sum / (i - k);
                square_sum += ((double)(i - k) / (double)(i - k + 1)) * diff * diff;
                index_parameter_diff = seq_index_parameter[k] - index_parameter_sum / (i - k);
                index_parameter_square_sum += ((double)(i - k) / (double)(i - k + 1)) *
                                              index_parameter_diff * index_parameter_diff;
                mix_square_sum += ((double)(i - k) / (double)(i - k + 1)) * diff *  index_parameter_diff;
                sum += real_sequence[index][j][k];
                index_parameter_sum += seq_index_parameter[k];

                if ((k < i - 1) && (index_parameter_square_sum > 0.)) {
                  residual[k] = square_sum - mix_square_sum * mix_square_sum / index_parameter_square_sum;
                }
                else {
                  residual[k] = 0.;
                }
              }
            }
          }

          else {
            if (type[j] != REAL_VALUE) {
              square_sum = 0.;
              sum = int_sequence[index][j][i];
              residual[i] = 0.;

              for (k = i - 1;k >= 0;k--) {
                diff = int_sequence[index][j][k] - sum / (i - k);
                square_sum += ((double)(i - k) / (double)(i - k + 1)) * diff * diff;
                sum += int_sequence[index][j][k];
                residual[k] = square_sum;
              }

#             ifdef MESSAGE
//              cout << "\n";

              square_sum = int_sequence[index][j][i] * int_sequence[index][j][i];
              sum = int_sequence[index][j][i];

              for (k = i - 1;k >= 0;k--) {
                square_sum += int_sequence[index][j][k] * int_sequence[index][j][k];
                sum += int_sequence[index][j][k];

                if ((square_sum - sum * sum / (i - k + 1) < residual[k] - DOUBLE_ERROR) ||
                    (square_sum - sum * sum / (i - k + 1) > residual[k] + DOUBLE_ERROR)) {
                  cout << "\nERROR: " << i << " " << k << " | " << square_sum - sum * sum / (i - k + 1)
                       << " " << residual[k] << endl;
                }
              }

/*              mean_square_diff[j][i] = 0.;
              sum = 0.;

              for (k = i - 1;k >= 0;k--) {
                diff = int_sequence[index][j][i] - int_sequence[index][j][k];
                sum += diff * diff;
                mean_square_diff[j][k] += sum;

                if ((mean_square_diff[j][k] / (i - k + 1) < residual[k] - DOUBLE_ERROR) ||
                    (mean_square_diff[j][k] / (i - k + 1) > residual[k] + DOUBLE_ERROR)) {
                  cout << "\nERROR: " << i << " " << k << " | " << mean_square_diff[j][k] / (i - k + 1)
                       << " " << residual[k] << endl;
                }
              } */
#             endif

            }

            else {
              square_sum = 0.;
              sum = real_sequence[index][j][i];
              residual[i] = 0.;

              for (k = i - 1;k >= 0;k--) {
                diff = real_sequence[index][j][k] - sum / (i - k);
                square_sum += ((double)(i - k) / (double)(i - k + 1)) * diff * diff;
                sum += real_sequence[index][j][k];
                residual[k] = square_sum;
              }
            }
          }

          if ((model_type[0] == MEAN_CHANGE) || (model_type[0] == INTERCEPT_SLOPE_CHANGE)) {
            for (k = i - 1;k >= 0;k--) {
              contrast[k] -= residual[k];
            }
          }

          else if (model_type[0] == MEAN_VARIANCE_CHANGE) {
            for (k = i - 1;k >= 0;k--) {
              contrast[k] += residual[k];
            }
          }

          else {
            for (k = i;k >= 0;k--) {
//              if ((contrast[k] != D_INF) && (residual[k] > 0.)) {
              if ((contrast[k] != D_INF) && (residual[k] > sqrt((double)(i - k + 1)) * ROUNDOFF_ERROR)) {
                contrast[k] -= ((double)(i - k + 1) / 2.) * (logl(residual[k] /
                                 (i - k + 1)) + log(2 * M_PI) + 1);
/*                contrast[k] -= ((double)(i - k + 1) / 2.) * (logl(residual[k] /
                                 (i - k)) + log(2 * M_PI)) + (double)(i - k) / 2.; */
              }
              else {
                contrast[k] = D_INF;
              }
            }
          }
        }
      }

      if (model_type[0] == MEAN_VARIANCE_CHANGE) {
        contrast[i] = D_INF;
        for (j = i - 1;j >= 0;j--) {
//          if (contrast[j] > 0.) {
          if (contrast[j] > sqrt((double)((nb_variable - 1) * (i - j + 1))) * ROUNDOFF_ERROR) {
            contrast[j] = -((double)((nb_variable - 1) * (i - j + 1)) / 2.) * (logl(contrast[j] /
                             ((nb_variable - 1) * (i - j + 1))) + log(2 * M_PI) + 1);
/*            contrast[j] = -((double)((nb_variable - 1) * (i - j + 1)) / 2.) * (logl(contrast[j] /
                             ((nb_variable - 1) * (i - j))) + log(2 * M_PI)) +
                           (double)((nb_variable - 1) * (i - j)) / 2.; */
          }
          else {
            contrast[j] = D_INF;
          }
        }
      }
    }

#   ifdef DEBUG
    for (j = i - 1;j >= 0;j--) {
      cout << contrast[j] << "  ";
    }
    cout << endl;
#   endif

    // calcul du nombre de segmentations

    for (j = 0;j < nb_segment;j++) {
      nb_segmentation_forward[i][j] = 0;
    }

    for (j = MAX(0 , nb_segment + i - length[index]);j < MIN((i < length[index] - 1 ? nb_segment - 1 : nb_segment) , i + 1);j++) {
      if (j == 0) {
        if (contrast[0] != D_INF) {
          nb_segmentation_forward[i][j]++;
        }
      }

      else {
        for (k = i;k >= j;k--) {
          if (contrast[k] != D_INF) {
            nb_segmentation_forward[i][j] += nb_segmentation_forward[k - 1][j - 1];
          }
        }
      }
    }

#   ifdef DEBUG
    nb_segmentation[0] = 1;
    for (j = 1;j < MIN((i < length[index] - 1 ? nb_segment - 1 : nb_segment) , i + 1);j++) {
      nb_segmentation[j] = nb_segmentation[j - 1] * (i - j + 1) / j;
    }

    if (i < inb_segmentation) {
      cout << i << ": ";
      for (j = 0;j < MIN((i < length[index] - 1 ? nb_segment - 1 : nb_segment) , i + 1);j++) {
        cout << nb_segmentation_forward[i][j] << " " << nb_segmentation[j] << " | ";
      }
      cout << endl;
    }
#   endif

    for (j = 0;j < MIN((i < length[index] - 1 ? nb_segment - 1 : nb_segment) , i + 1);j++) {
      nb_segmentation[j] = nb_segmentation_forward[i][j];
      if (nb_segmentation[j] > inb_segmentation) {
        nb_segmentation[j] = inb_segmentation;
      }
    }

    for (j = MAX(0 , nb_segment + i - length[index]);j < MIN((i < length[index] - 1 ? nb_segment - 1 : nb_segment) , i + 1);j++) {
      if (j == 0) {
        forward[i][j][0] = contrast[0];
        if (forward[i][j][0] != D_INF) {
          optimal_length[i][j][0] = i + 1;
        }
      }

      else {
/*        if (j < nb_segment - 1) {
          if ((i > j) && (nb_segmentation[j] < inb_segmentation)) {
            nb_segmentation[j] = nb_segmentation[j] * i / (i - j);
            if (nb_segmentation[j] > inb_segmentation) {
              nb_segmentation[j] = inb_segmentation;
            }
          }
        }

        else {
          nb_segmentation[j] = nb_segmentation[j - 1] * (i - j + 1) / j;
          if (nb_segmentation[j] > inb_segmentation) {
            nb_segmentation[j] = inb_segmentation;
          }
        } */

#       ifdef DEBUG
        cout << "TEST: " << i << " " << j << ": " << nb_segmentation[j] << endl;
#       endif

        for (k = i;k >= j;k--) {
          rank[k] = 0;
        }

        for (k = 0;k < nb_segmentation[j];k++) {
          forward[i][j][k] = D_INF;
          for (m = i;m >= j;m--) {
            if ((contrast[m] != D_INF) && (forward[m - 1][j - 1][rank[m]] != D_INF)) {
              buff = contrast[m] + forward[m - 1][j - 1][rank[m]];
              if (buff > forward[i][j][k]) {
                forward[i][j][k] = buff;
                optimal_length[i][j][k] = i - m + 1;
                optimal_rank[i][j][k] = rank[m];
              }
            }
          }

          if (forward[i][j][k] != D_INF) {
            rank[i - optimal_length[i][j][k] + 1]++;
          }

#         ifdef DEBUG
          else {
            cout << "\ntest utile" << endl;
          }
#         endif

        }
      }

      for (k = (int)nb_segmentation[j];k < inb_segmentation;k++) {
        forward[i][j][k] = D_INF;
      }
    }

#   ifdef DEBUG
    cout << i << " : ";
    for (j = MAX(0 , nb_segment + i - length[index]);j < MIN((i < length[index] - 1 ? nb_segment - 1 : nb_segment) , i + 1);j++) {
      cout << j << " :";
      for (k = 0;k < nb_segmentation[j];k++) {
        cout << " " << forward[i][j][k];
        if (forward[i][j][k] != D_INF) {
          cout << " " << optimal_length[i][j][k];
        }
        cout << " |";
      }
      cout << "| ";
    }
    cout << endl;
#   endif

  }

# ifdef MESSAGE
  buff = 1.;
  for (i = 1;i < nb_segment;i++) {
    buff *= (double)(length[index] - i) / (double)i;
//    buff = buff * (length[index] - i) / i;
  }

  os.precision(10);

  os << "\n" << SEQ_label[SEQL_NB_SEGMENTATION] << ": "
     << nb_segmentation_forward[length[index] - 1][nb_segment - 1] << " (" << buff << ")" << endl;

  os.precision(6);

  if (((model_type[0] == MEAN_CHANGE) || (model_type[0] == INTERCEPT_SLOPE_CHANGE)) &&
      (format == 's') && (nb_segment == 2) && (inb_segmentation >= length[index] - 1)) {
    change_point = new int[length[index]];
  }
# endif

  // restauration

  likelihood_cumul = 0.;
  posterior_probability_cumul = 0.;

  for (i = 0;i < nb_segmentation[nb_segment - 1];i++) {
    if (forward[length[index] - 1][nb_segment - 1][i] == D_INF) {
      break;
    }

#   ifdef DEBUG
    cout << "\n";
#   endif

    j = length[index] - 1;
    psegment = int_sequence[index][0] + j;
    brank = i;

    for (k = nb_segment - 1;k >= 0;k--) {
      for (m = j;m > j - optimal_length[j][k][brank];m--) {
        active_cell[m][k] = true;
        *psegment-- = k;
      }

      for (m = 1;m < nb_variable;m++) {
        if ((model_type[m - 1] == POISSON_CHANGE) || (model_type[0] == MULTIVARIATE_POISSON_CHANGE) ||
            (model_type[m - 1] == GEOMETRIC_0_CHANGE) || (model_type[m - 1] == GEOMETRIC_1_CHANGE) ||
            (model_type[0] == MULTIVARIATE_GEOMETRIC_0_CHANGE) || (model_type[m - 1] == GAUSSIAN_CHANGE) ||
            (model_type[m - 1] == VARIANCE_CHANGE) || (model_type[0] == MEAN_CHANGE) ||
            (model_type[0] == MEAN_VARIANCE_CHANGE) || (model_type[m - 1] == BAYESIAN_POISSON_CHANGE) ||
            (model_type[m - 1] == BAYESIAN_GAUSSIAN_CHANGE)) {
          mean[m][k] = 0.;
          if (type[m] != REAL_VALUE) {
            for (n = j;n > j - optimal_length[j][k][brank];n--) {
              mean[m][k] += int_sequence[index][m][n];
            }
          }
          else {
            for (n = j;n > j - optimal_length[j][k][brank];n--) {
              mean[m][k] += real_sequence[index][m][n];
            }
          }
          mean[m][k] /= optimal_length[j][k][brank];

          if ((model_type[m - 1] != GEOMETRIC_0_CHANGE) && (model_type[m - 1] != GEOMETRIC_1_CHANGE) &&
              (model_type[0] != MULTIVARIATE_GEOMETRIC_0_CHANGE)) {
            variance[m][k] = 0.;
            if (optimal_length[j][k][brank] > 1) {
              if (type[m] != REAL_VALUE) {
                for (n = j;n > j - optimal_length[j][k][brank];n--) {
                  diff = int_sequence[index][m][n] - mean[m][k];
                  variance[m][k] += diff * diff;
                }
              }
              else {
                for (n = j;n > j - optimal_length[j][k][brank];n--) {
                  diff = real_sequence[index][m][n] - mean[m][k];
                  variance[m][k] += diff * diff;
                }
              }

//              variance[m][k] /= optimal_length[j][k][brank];
              variance[m][k] /= (optimal_length[j][k][brank] - 1);
            }
          }
        }

        if ((model_type[m - 1] == LINEAR_MODEL_CHANGE) || (model_type[0] == INTERCEPT_SLOPE_CHANGE)) {
          response_mean = 0.;
          if (type[m] != REAL_VALUE) {
            for (n = j;n > j - optimal_length[j][k][brank];n--) {
              response_mean += int_sequence[index][m][n];
            }
          }
          else {
            for (n = j;n > j - optimal_length[j][k][brank];n--) {
              response_mean += real_sequence[index][m][n];
            }
          }
          response_mean /= optimal_length[j][k][brank];

          index_parameter_mean = 0.;
          for (n = j;n > j - optimal_length[j][k][brank];n--) {
            index_parameter_mean += seq_index_parameter[n];
          }
          index_parameter_mean /= optimal_length[j][k][brank];

          index_parameter_variance = 0.;
          for (n = j;n > j - optimal_length[j][k][brank];n--) {
            diff = seq_index_parameter[n] - index_parameter_mean;
            index_parameter_variance += diff * diff;                
          }

          covariance = 0.;
          if (type[m] != REAL_VALUE) {
            for (n = j;n > j - optimal_length[j][k][brank];n--) {
              covariance += (int_sequence[index][m][n] - response_mean) * (seq_index_parameter[n] - index_parameter_mean);
            }
          }
          else {
            for (n = j;n > j - optimal_length[j][k][brank];n--) {
              covariance += (real_sequence[index][m][n] - response_mean) * (seq_index_parameter[n] - index_parameter_mean);
            }
          }

          slope[m][k] = covariance / index_parameter_variance;
          intercept[m][k] = response_mean - slope[m][k] * index_parameter_mean;

          residual_mean = 0.;
          if (type[m] != REAL_VALUE) {
            for (n = j;n > j - optimal_length[j][k][brank];n--) {
              residual_mean += int_sequence[index][m][n] - (intercept[m][k] + slope[m][k] * seq_index_parameter[n]);
            }
          }
          else {
            for (n = j;n > j - optimal_length[j][k][brank];n--) {
              residual_mean += real_sequence[index][m][n] - (intercept[m][k] + slope[m][k] * seq_index_parameter[n]);
            }
          }
          residual_mean /= optimal_length[j][k][brank];

          variance[m][k] = 0.;
          if (optimal_length[j][k][brank] > 2) {
            if (type[m] != REAL_VALUE) {
              for (n = j;n > j - optimal_length[j][k][brank];n--) {
                diff = int_sequence[index][m][n] - (intercept[m][k] + slope[m][k] * seq_index_parameter[n]) - residual_mean;
                variance[m][k] += diff * diff;
              }
            }
            else {
              for (n = j;n > j - optimal_length[j][k][brank];n--) {
                diff = real_sequence[index][m][n] - (intercept[m][k] + slope[m][k] * seq_index_parameter[n]) - residual_mean;
                variance[m][k] += diff * diff;
              }
            }
            variance[m][k] /= (optimal_length[j][k][brank] - 2);
          }
        }
      }

#     ifdef DEBUG
      cout << k << " " << optimal_length[j][k][brank] << " " << brank << " | ";
#     endif

      if (k > 0) {
        previous_rank = optimal_rank[j][k][brank];
        j -= optimal_length[j][k][brank];
        brank = previous_rank;
      }
    }

#   ifdef DEBUG
    cout << endl;
#   endif

    if ((model_type[0] == MEAN_CHANGE) || (model_type[0] == INTERCEPT_SLOPE_CHANGE)) {
      if (forward[length[index] - 1][nb_segment - 1][i] < 0.) {
        forward[length[index] - 1][nb_segment - 1][i] = -((double)((nb_variable - 1) * length[index]) / 2.) *
                                                         (log(-forward[length[index] - 1][nb_segment - 1][i] /
                                                           ((nb_variable - 1) * length[index])) + log(2 * M_PI) + 1);
/*        forward[length[index] - 1][nb_segment - 1][i] = -((double)((nb_variable - 1) * length[index]) / 2.) *
                                                         (log(-forward[length[index] - 1][nb_segment - 1][i] /
                                                           ((nb_variable - 1) * (length[index] - nb_segment))) + log(2 * M_PI)) -
                                                         (double)((nb_variable - 1) * (length[index] - nb_segment)) / 2.; */
      }
      else {
        forward[length[index] - 1][nb_segment - 1][i] = D_INF;
      }
    }

    if (i == 0) {
      segmentation_likelihood = forward[length[index] - 1][nb_segment - 1][i];
    }

    if (forward[length[index] - 1][nb_segment - 1][i] != D_INF) {
      likelihood_cumul += exp(forward[length[index] - 1][nb_segment - 1][i]);
      if (likelihood != D_INF) {
        posterior_probability_cumul += exp(forward[length[index] - 1][nb_segment - 1][i] - likelihood);
      }
    }

#   ifdef DEBUG
    psegment = int_sequence[index][0];
    for (j = 0;j < length[index];j++) {
//      if (((i == 0) || (*psegment != *(psegment - 1))) &&
//          (forward[length[index] - 1][nb_segment - 1][i] > segment_probability[j][*psegment])) {
      if (forward[length[index] - 1][nb_segment - 1][i] > segment_probability[j][*psegment]) {
        segment_probability[j][*psegment] = forward[length[index] - 1][nb_segment - 1][i];
      }
      psegment++;
    }
#   endif

#   ifdef MESSAGE
    if (inb_segmentation >= 1000) {

      // approximation des probabilites lissees

      buff = exp(forward[length[index] - 1][nb_segment - 1][i]);
      approximated_likelihood += buff;
      psegment = int_sequence[index][0];
      for (j = 0;j < length[index];j++) {
        smoothed_probability[j][*psegment++] += buff;
      }
    }
#   endif

    nb_cell = 0;
    for (j = 0;j < length[index];j++) {
      for (k = 0;k < nb_segment;k++) {
        if (active_cell[j][k]) {
          nb_cell++;
        }
      }
    }

#   ifdef MESSAGE
    if (i == 0) {
      os << "\n";
    }

    switch (format) {

    case 'a' : {
      if (inb_segmentation <= 200) {
      psegment = int_sequence[index][0];
      for (j = 0;j < length[index];j++) {
        os << *psegment++ << " ";
      }

      os << "  " << i + 1 << "  " << forward[length[index] - 1][nb_segment - 1][i] << "   (";
      if (likelihood != D_INF) {
        os << exp(forward[length[index] - 1][nb_segment - 1][i] - likelihood) << "  ";
        if (boost::math::isnan(likelihood_cumul)) {
          os << likelihood_cumul / exp(likelihood);
        }
        else {
          os << posterior_probability_cumul;
        }
      }
      else {
        os << exp(forward[length[index] - 1][nb_segment - 1][i] - segmentation_likelihood);
      }
      os << "  " << nb_cell << ")" << endl;

      os << (nb_segment == 2 ? SEQ_label[SEQL_CHANGE_POINT] : SEQ_label[SEQL_CHANGE_POINTS]) << ": ";

      psegment = int_sequence[index][0] + 1;
      if (index_parameter) {
        for (j = 1;j < length[index];j++) {
          if (*psegment != *(psegment - 1)) {
            os << index_parameter[index][j] << ", ";
          }
          psegment++;
        }
      }

      else {
        for (j = 1;j < length[index];j++) {
          if (*psegment != *(psegment - 1)) {
            os << j << ", ";
          }
          psegment++;
        }
      }
      os << endl;

      for (j = 1;j < nb_variable;j++) {
        if ((model_type[j - 1] == POISSON_CHANGE) || (model_type[0] == MULTIVARIATE_POISSON_CHANGE) ||
            (model_type[j - 1] == BAYESIAN_POISSON_CHANGE)) {
          if (nb_variable > 2) {
            os << STAT_label[STATL_VARIABLE] << " " << j << "   ";
          }
          os << SEQ_label[SEQL_SEGMENT] << " " << STAT_label[STATL_MEAN] << ", "
             << STAT_label[STATL_VARIANCE] << ": ";
          for (k = 0;k < nb_segment;k++) {
            os << mean[j][k] << " " << variance[j][k];
            if (k < nb_segment - 1) {
              os << " | ";
            }
          }
          os << endl;
        }

        else if ((model_type[j - 1] == GEOMETRIC_0_CHANGE) || (model_type[j - 1] == GEOMETRIC_1_CHANGE) ||
                 (model_type[0] == MULTIVARIATE_GEOMETRIC_0_CHANGE)) {
          if (nb_variable > 2) {
            os << STAT_label[STATL_VARIABLE] << " " << j << "   ";
          }
          os << SEQ_label[SEQL_SEGMENT] << " " << STAT_label[STATL_MEAN] << ": ";
          for (k = 0;k < nb_segment;k++) {
            os << mean[j][k];
            if (k < nb_segment - 1) {
              os << " ";
            }
          }
          os << endl;
        }

        else if ((model_type[j - 1] == GAUSSIAN_CHANGE) || (model_type[j - 1] == VARIANCE_CHANGE) ||
                 (model_type[0] == MEAN_CHANGE) || (model_type[0] == MEAN_VARIANCE_CHANGE) ||
                 (model_type[j - 1] == BAYESIAN_GAUSSIAN_CHANGE)) {
          if (nb_variable > 2) {
            os << STAT_label[STATL_VARIABLE] << " " << j << "   ";
          }
          os << SEQ_label[SEQL_SEGMENT] << " " << STAT_label[STATL_MEAN] << ", "
             << STAT_label[STATL_STANDARD_DEVIATION] << ": ";
          for (k = 0;k < nb_segment;k++) {
            os << mean[j][k] << " " << sqrt(variance[j][k]);
            if (k < nb_segment - 1) {
              os << " | ";
            }
          }
          os << endl;
        }

        else if ((model_type[j - 1] == LINEAR_MODEL_CHANGE) || (model_type[0] == INTERCEPT_SLOPE_CHANGE)) {
          if (nb_variable > 2) {
            os << STAT_label[STATL_VARIABLE] << " " << j << "   ";
          }
          os << SEQ_label[SEQL_SEGMENT] << " " << STAT_label[STATL_INTERCEPT] << ", "
             << STAT_label[STATL_SLOPE] << ", " << STAT_label[STATL_RESIDUAL] << " "
             << STAT_label[STATL_STANDARD_DEVIATION] << ": ";
          for (k = 0;k < nb_segment;k++) {
            os << intercept[j][k] << " " << slope[j][k] << " " << sqrt(variance[j][k]);
            if (k < nb_segment - 1) {
              os << " | ";
            }
          }
          os << endl;
        }
      }
      }
      break;
    }

    case 's' : {
      psegment = int_sequence[index][0];
      for (j = 0;j < length[index];j++) {
        os << *psegment++ << "\t";
      }

      os << "\t" << i + 1 << "\t" << forward[length[index] - 1][nb_segment - 1][i] << "\t";
      if (likelihood != D_INF) {
        os << exp(forward[length[index] - 1][nb_segment - 1][i] - likelihood) << "\t";
        if (boost::math::isnan(likelihood_cumul)) {
          os << posterior_probability_cumul;
        }
        else {
          os << likelihood_cumul / exp(likelihood);
        }
      }
      else {
        os << exp(forward[length[index] - 1][nb_segment - 1][i] - segmentation_likelihood);
      }
      os << "\t" << nb_cell << endl;

      for (j = 1;j < nb_variable;j++) {
        if ((model_type[j - 1] == POISSON_CHANGE) || (model_type[0] == MULTIVARIATE_POISSON_CHANGE) ||
            (model_type[j - 1] == BAYESIAN_POISSON_CHANGE)) {
          if (nb_variable > 2) {
            os << STAT_label[STATL_VARIABLE] << "\t" << j << "\t";
          }
          os << SEQ_label[SEQL_SEGMENT] << "\t" << STAT_label[STATL_MEAN] << "\t"
             << STAT_label[STATL_VARIANCE];
          for (k = 0;k < nb_segment;k++) {
            os << "\t" << mean[j][k] << "\t" << variance[j][k];
          }
          os << endl;
        }

        else if ((model_type[j - 1] == GEOMETRIC_0_CHANGE) || (model_type[j - 1] == GEOMETRIC_1_CHANGE) ||
                 (model_type[0] == MULTIVARIATE_GEOMETRIC_0_CHANGE)) {
          if (nb_variable > 2) {
            os << STAT_label[STATL_VARIABLE] << "\t" << j << "\t";
          }
          os << SEQ_label[SEQL_SEGMENT] << "\t" << STAT_label[STATL_MEAN];
          for (k = 0;k < nb_segment;k++) {
            os << "\t" << mean[j][k];
          }
          os << endl;
        }

        else if ((model_type[j - 1] == GAUSSIAN_CHANGE) || (model_type[j - 1] == VARIANCE_CHANGE) ||
                 (model_type[0] == MEAN_CHANGE) || (model_type[0] == MEAN_VARIANCE_CHANGE) ||
                 (model_type[j - 1] == BAYESIAN_GAUSSIAN_CHANGE)) {
          if (nb_variable > 2) {
            os << STAT_label[STATL_VARIABLE] << "\t" << j << "\t";
          }
          os << SEQ_label[SEQL_SEGMENT] << "\t" << STAT_label[STATL_MEAN] << "\t"
             << STAT_label[STATL_STANDARD_DEVIATION];
          for (k = 0;k < nb_segment;k++) {
            os << "\t" << mean[j][k] << "\t" << sqrt(variance[j][k]);
          }
          os << endl;
        }

        else if ((model_type[j - 1] == LINEAR_MODEL_CHANGE) || (model_type[0] == INTERCEPT_SLOPE_CHANGE)) {
          if (nb_variable > 2) {
            os << STAT_label[STATL_VARIABLE] << "\t" << j << "\t";
          }
          os << SEQ_label[SEQL_SEGMENT] << "\t" << STAT_label[STATL_INTERCEPT] << "\t"
             << STAT_label[STATL_SLOPE] << "\t" << STAT_label[STATL_RESIDUAL] << " "
             << STAT_label[STATL_STANDARD_DEVIATION];
          for (k = 0;k < nb_segment;k++) {
            os << "\t" << intercept[j][k] << "\t" << slope[j][k] << "\t" << sqrt(variance[j][k]);
          }
          os << endl;
        }
      }

      if (((model_type[0] == MEAN_CHANGE) || (model_type[0] == INTERCEPT_SLOPE_CHANGE)) &&
          (nb_segment == 2) && (inb_segmentation >= length[index] - 1)) {
        psegment = int_sequence[index][0] + 1;
        for (j = 1;j < length[index];j++) {
          if (*psegment != *(psegment - 1)) {
            change_point[i] = j;
            break;
          }
          psegment++;
        }
      }
      break;
    }
    }
#   endif

  }

# ifdef MESSAGE
  if (((model_type[0] == MEAN_CHANGE) || (model_type[0] == INTERCEPT_SLOPE_CHANGE)) &&
      (format == 's') && (nb_segment == 2) && (inb_segmentation >= length[index] - 1)) {
    norm = 0.;
    for (i = 0;i < length[index] - 1;i++) {
      norm += exp(forward[length[index] - 1][nb_segment - 1][i]);
    }

    os << "\n" << SEQ_label[SEQL_POSTERIOR_CHANGE_POINT_PROBABILITY] << "\n\n";

    if (index_parameter) {
      for (i = 1;i < length[index];i++) {
        for (j = 0;j < length[index] - 1;j++) {
          if (change_point[j] == i) {
            os << index_parameter[index][i] << "\t"
               << exp(forward[length[index] - 1][nb_segment - 1][j]) / norm << endl;
            break;
          }
        }
      }
    }

    else {
      for (i = 1;i < length[index];i++) {
        for (j = 0;j < length[index] - 1;j++) {
          if (change_point[j] == i) {
            os << i << "\t" << exp(forward[length[index] - 1][nb_segment - 1][j]) / norm << endl;
            break;
          }
        }
      }
    }

    delete [] change_point;
  }

  if (((model_type[0] == MEAN_CHANGE) || (model_type[0] == INTERCEPT_SLOPE_CHANGE)) &&
       (nb_segment > 2) && (inb_segmentation >= nb_segmentation_forward[length[index] - 1][nb_segment - 1])) {
    norm = 0.;
    for (i = 0;i < nb_segmentation_forward[length[index] - 1][nb_segment - 1];i++) {
      norm += exp(forward[length[index] - 1][nb_segment - 1][i]);
    }

    os << SEQ_label[SEQL_POSTERIOR_PROBABILITY] << ": " << exp(forward[length[index] - 1][nb_segment - 1][0]) / norm << endl;
  }
# endif

# ifdef DEBUG
  if (((likelihood != D_INF) && (likelihood_cumul / exp(likelihood) > 0.8)) ||
      (segmentation_likelihood != D_INF)) {
    if (likelihood != D_INF) {
      for (i = 0;i < length[index];i++) {
        for (j = 0;j < nb_segment;j++) {
          if (segment_probability[i][j] != D_INF) {
            segment_probability[i][j] = exp(segment_probability[i][j] - likelihood);
          }
          else {
            segment_probability[i][j] = 0.;
          }
        }
      }

      os << "\n" << SEQ_label[SEQL_MAX_POSTERIOR_SEGMENT_PROBABILITY] << "\n\n";
    }

    else {
      for (i = 0;i < length[index];i++) {
        for (j = 0;j < nb_segment;j++) {
          if (segment_probability[i][j] != D_INF) {
            segment_probability[i][j] = exp(segment_probability[i][j] - segmentation_likelihood);
          }
          else {
            segment_probability[i][j] = 0.;
          }
        }
      }

      os << "\n" << SEQ_label[SEQL_MAX_SEGMENT_LIKELIHOOD] << "\n\n";
    }

    psegment = int_sequence[index][0];
    for (j = 0;j < length[index];j++) {
      *psegment++ = I_DEFAULT;
    }

    profile_ascii_print(os , index , nb_segment , segment_probability ,
                        SEQ_label[SEQL_SEGMENT]);
  }
# endif

# ifdef MESSAGE
  if (inb_segmentation >= 1000) {
    double previous_cumul[3] ,
           cdf[10] = {0.5 , 0.75 , 0.9 , 0.95 , 0.975 , 0.99 , 0.995, 0.9975 , 0.999 , 1};
    long double divergence;
    Distribution *segmentation;


    for (i = 0;i < length[index];i++) {
      for (j = 0;j < nb_segment;j++) {
        smoothed_probability[i][j] /= approximated_likelihood;
      }
    }

    psegment = int_sequence[index][0];
    for (i = 0;i < length[index];i++) {
      *psegment++ = I_DEFAULT;
    }

    os << "\n" << SEQ_label[SEQL_POSTERIOR_SEGMENT_PROBABILITY] << "\n\n";

    profile_ascii_print(os , index , nb_segment , smoothed_probability , SEQ_label[SEQL_SEGMENT]);

    // approximation de la divergence de Kullback-Leibler de la loi uniforme par rapport
    // a la loi des segmentations

    segmentation = new Distribution(inb_segmentation);
    likelihood_cumul = 0.;
    divergence = 0.;

    if (likelihood != D_INF) {
      for (i = 0;i < inb_segmentation;i++) {
        segmentation->mass[i] = exp(forward[length[index] - 1][nb_segment - 1][i] - likelihood);
        likelihood_cumul += exp(forward[length[index] - 1][nb_segment - 1][i] - likelihood);
        segmentation->cumul[i] = likelihood_cumul;

        divergence += exp(forward[length[index] - 1][nb_segment - 1][i] - likelihood) *
                      (forward[length[index] - 1][nb_segment - 1][i] - likelihood +
                       log(nb_segmentation_forward[length[index] - 1][nb_segment - 1]));
      }

      segmentation->complement = 1. - likelihood_cumul;
    }

    else {
      for (i = 0;i < inb_segmentation;i++) {
        segmentation->mass[i] = exp(forward[length[index] - 1][nb_segment - 1][i]) / approximated_likelihood;

        divergence += exp(forward[length[index] - 1][nb_segment - 1][i]) / approximated_likelihood *
                      (forward[length[index] - 1][nb_segment - 1][i] - log(approximated_likelihood) +
                       log(nb_segmentation_forward[length[index] - 1][nb_segment - 1]));
      }

      segmentation->cumul_computation();
    }

    segmentation->max_computation();
    segmentation->mean_computation();
    segmentation->variance_computation();

    os << "\n" << SEQ_label[SEQL_SEGMENTATION_DIVERGENCE] << ": " << divergence << endl;

    os << "\n";
    segmentation->ascii_characteristic_print(os , true);
    if (likelihood != D_INF) {
      os << STAT_label[STATL_CUMULATIVE] << " " << SEQ_label[SEQL_POSTERIOR_PROBABILITY] << ": "
         << likelihood_cumul << " (" << inb_segmentation << " "
         << SEQ_label[SEQL_SEGMENTATIONS]<< ")" << endl;
    }

    delete segmentation;

    os << "\n";
    likelihood_cumul = 0.;
    i = 0;
    if (likelihood != D_INF) {
      for (j = 0;j < inb_segmentation;j++) {
        previous_cumul[0] = likelihood_cumul;
        likelihood_cumul += exp(forward[length[index] - 1][nb_segment - 1][j]);
        if (likelihood_cumul / exp(likelihood) > cdf[i]) {
          os << j << " " << previous_cumul[0] / exp(likelihood) << " "
             << likelihood_cumul / exp(likelihood) << " (";
//          os << j + 1 << " " << likelihood_cumul / exp(likelihood) << " (";
          if (i == 0) {
            os << (j + 1) / nb_segmentation_forward[length[index] - 1][nb_segment - 1] << ")" << endl;
          }
          else {
            os << likelihood_cumul / exp(likelihood) - previous_cumul[1] << " "
               << (j + 1) / nb_segmentation_forward[length[index] - 1][nb_segment - 1] - previous_cumul[2] << ")" << endl;
          }

          if (cdf[i] == 1) {
            break;
          }
          previous_cumul[1] = likelihood_cumul / exp(likelihood);
          previous_cumul[2] = (j + 1) / nb_segmentation_forward[length[index] - 1][nb_segment - 1];
          i++;
        }
      }
    }

    else {
      for (j = 0;j < inb_segmentation;j++) {
        previous_cumul[0] = likelihood_cumul;
        likelihood_cumul += exp(forward[length[index] - 1][nb_segment - 1][j]);
        if (likelihood_cumul / approximated_likelihood > cdf[i]) {
          os << j << " " << previous_cumul[0] / approximated_likelihood << " "
             << likelihood_cumul / approximated_likelihood << " (";
//          os << j + 1 << " " << likelihood_cumul / approximated_likelihood << " (";
          if (i == 0) {
            os << (j + 1) / nb_segmentation_forward[length[index] - 1][nb_segment - 1] << ")" << endl;
          }
          else {
            os << likelihood_cumul / approximated_likelihood - previous_cumul[1] << " "
               << (j + 1) / nb_segmentation_forward[length[index] - 1][nb_segment - 1] - previous_cumul[2] << ")" << endl;
          }

          if (cdf[i] == 1) {
            break;
          }
          previous_cumul[1] = likelihood_cumul / approximated_likelihood;
          previous_cumul[2] = (j + 1) / nb_segmentation_forward[length[index] - 1][nb_segment - 1];
          i++;
        }
      }
    }

/*    ofstream out_file("Spreadsheet/segmentation_probability.xld");
  
    likelihood_cumul = 0.;
    if (likelihood != D_INF) {
      for (i = 0;i < inb_segmentation;i++) {
        likelihood_cumul += exp(forward[length[index] - 1][nb_segment - 1][i]);
        out_file << i + 1 << "\t" << exp(forward[length[index] - 1][nb_segment - 1][i] - likelihood) << "\t"
                 << likelihood_cumul / exp(likelihood) << "\t"
                 << 1. / nb_segmentation_forward[length[index] - 1][nb_segment - 1] << endl;
      }
    }

    else {
      for (i = 0;i < inb_segmentation;i++) {
        likelihood_cumul += exp(forward[length[index] - 1][nb_segment - 1][i]);
        out_file << i + 1 << "\t" << exp(forward[length[index] - 1][nb_segment - 1][i]) / approximated_likelihood << "\t"
                 << likelihood_cumul / approximated_likelihood << "\t"
                 << 1. / nb_segmentation_forward[length[index] - 1][nb_segment - 1] << endl;
      }
    } */
  }
# endif

  delete [] frequency;

  for (i = 1;i < nb_variable;i++) {
    delete [] factorial[i];
  }
  delete [] factorial;

  for (i = 1;i < nb_variable;i++) {
    delete [] hyperparam[i];
  }
  delete [] hyperparam;

  delete [] sequence_mean;
  delete [] residual;

  if (index_parameter_type == IMPLICIT_TYPE) {
    delete [] seq_index_parameter;
  }

  delete [] contrast;

  for (i = 0;i < length[index];i++) {
    delete [] nb_segmentation_forward[i];
  }
  delete [] nb_segmentation_forward;

  for (i = 0;i < length[index];i++) {
    for (j = 0;j < nb_segment;j++) {
      delete [] forward[i][j];
    }
    delete [] forward[i];
  }
  delete [] forward;

  delete [] nb_segmentation;
  delete [] rank;

  for (i = 0;i < length[index];i++) {
    for (j = 0;j < nb_segment;j++) {
      delete [] optimal_length[i][j];
    }
    delete [] optimal_length[i];
  }
  delete [] optimal_length;

  for (i = 0;i < length[index];i++) {
    for (j = 0;j < nb_segment;j++) {
      delete [] optimal_rank[i][j];
    }
    delete [] optimal_rank[i];
  }
  delete [] optimal_rank;

  for (i = 1;i < nb_variable;i++) {
    delete [] mean[i];
    delete [] variance[i];
    delete [] intercept[i];
    delete [] slope[i];
  }
  delete [] mean;
  delete [] variance;
  delete [] intercept;
  delete [] slope;

  for (i = 0;i < length[index];i++) {
    delete [] active_cell[i];
  }
  delete [] active_cell;

# ifdef MESSAGE
  for (i = 1;i < nb_variable;i++) {
    delete [] mean_square_diff[i];
  }
  delete [] mean_square_diff;
# endif

# ifdef DEBUG
  for (i = 0;i < length[index];i++) {
    delete [] segment_probability[i];
  }
  delete [] segment_probability;
# endif

# ifdef MESSAGE
  if (inb_segmentation >= 1000) {
    for (i = 0;i < length[index];i++) {
      delete [] smoothed_probability[i];
    }
    delete [] smoothed_probability;
  }
# endif

  return segmentation_likelihood;
}


/*--------------------------------------------------------------*
 *
 *  Calcul par maximisation des profils de segments/ruptures d'une sequence.
 *
 *  arguments : indice de la sequence, nombre de segments, types des modeles,
 *              rangs (variables ordinales), stream, pointeur sur un objet MultiPlotSet,
 *              type de sortie, format de sortie ('a' : ASCII, 's' : Spreadsheet,
 *              'g' : Gnuplot, 'p' : plotable), vraisemblance des donnees.
 *
 *--------------------------------------------------------------*/

double Sequences::forward_backward_dynamic_programming(int index , int nb_segment ,
                                                       int *model_type , double **rank ,
                                                       ostream *os , MultiPlotSet *plot_set ,
                                                       int output , char format ,
                                                       double likelihood) const

{
  register int i , j , k , m;
  int max_nb_value , *frequency , *seq_index_parameter = NULL , *psegment , **optimal_length;
  double sum , factorial_sum , proba , diff , index_parameter_sum , index_parameter_diff , buff ,
         segmentation_likelihood , backward_max , response_mean , index_parameter_mean ,
         index_parameter_variance , covariance , intercept , slope , *sequence_mean , **hyperparam ,
         **factorial , **forward , **backward , **backward_output , **piecewise_function;
  long double square_sum , index_parameter_square_sum , mix_square_sum , prior_contrast ,
              *residual , *contrast;


  max_nb_value = 0;
  factorial = new double*[nb_variable];
  hyperparam = new double*[nb_variable];

  for (i = 1;i < nb_variable;i++) {
    if (((model_type[i - 1] == CATEGORICAL_CHANGE) || (model_type[0] == MULTIVARIATE_CATEGORICAL_CHANGE)) &&
        (marginal_distribution[i]->nb_value > max_nb_value)) {
      max_nb_value = marginal_distribution[i]->nb_value;
    }

    if ((model_type[i - 1] == POISSON_CHANGE) || (model_type[0] == MULTIVARIATE_POISSON_CHANGE) ||
        (model_type[i - 1] == BAYESIAN_POISSON_CHANGE)) {
      factorial[i] = new double[length[index]];
    }
    else {
      factorial[i] = NULL;
    }

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

  if (max_nb_value > 0) {
    frequency = new int[max_nb_value];
  }
  else {
    frequency = NULL;
  }

  sequence_mean = new double[nb_variable];
  residual = new long double[length[index]];

  contrast = new long double[length[index]];

  forward = new double*[length[index]];
  for (i = 0;i < length[index];i++) {
    forward[i] = new double[nb_segment];
  }

  optimal_length = new int*[length[index]];
  for (i = 0;i < length[index];i++) {
    optimal_length[i] = new int[nb_segment];
  }

  backward = new double*[length[index]];
  for (i = 0;i < length[index];i++) {
    backward[i] = new double[nb_segment];
  }

  backward_output = new double*[length[index]];
  for (i = 0;i < length[index];i++) {
    backward_output[i] = new double[nb_segment];
  }

  piecewise_function = new double*[nb_variable];
  for (i = 1;i < nb_variable;i++) {
    if ((model_type[i - 1] == POISSON_CHANGE) || (model_type[0] == MULTIVARIATE_POISSON_CHANGE) ||
        (model_type[i - 1] == GEOMETRIC_0_CHANGE) || (model_type[i - 1] == GEOMETRIC_1_CHANGE) ||
        (model_type[0] == MULTIVARIATE_GEOMETRIC_0_CHANGE) || (model_type[i - 1] == GAUSSIAN_CHANGE) ||
        (model_type[i - 1] == VARIANCE_CHANGE) || (model_type[i - 1] == LINEAR_MODEL_CHANGE) ||
        (model_type[0] == MEAN_CHANGE) || (model_type[0] == MEAN_VARIANCE_CHANGE) ||
        (model_type[0] == INTERCEPT_SLOPE_CHANGE) || (model_type[i - 1] == BAYESIAN_POISSON_CHANGE) ||
        (model_type[i - 1] == BAYESIAN_GAUSSIAN_CHANGE)) {
      piecewise_function[i] = new double[length[index]];
    }
    else {
      piecewise_function[i] = NULL;
    }
  }

  for (i = 1;i < nb_variable;i++) {
    if (model_type[i - 1] == VARIANCE_CHANGE) {
      sequence_mean[i] = 0.;
      if (type[i] != REAL_VALUE) {
        for (j = 0;j < length[index];j++) {
          sequence_mean[i] += int_sequence[index][i][j];
        }
      }
      else {
        for (j = 0;j < length[index];j++) {
          sequence_mean[i] += real_sequence[index][i][j];
        }
      }
      sequence_mean[i] /= length[index];
    }

    if (((model_type[i - 1] == LINEAR_MODEL_CHANGE) && (!seq_index_parameter)) ||
        ((i == 1) && (model_type[0] == INTERCEPT_SLOPE_CHANGE))) {
      if (index_parameter_type == IMPLICIT_TYPE) {
        seq_index_parameter = new int[length[index]];
        for (j = 0;j < length[index];j++) {
          seq_index_parameter[j] = j;
        }
      }
      else {
        seq_index_parameter = index_parameter[index];
      }
    }
  }

  // recurrence "forward"

  for (i = 0;i < length[index];i++) {

    // calcul des log-vraisemblances des segments

    if (model_type[0] == MULTIVARIATE_CATEGORICAL_CHANGE) {
      for (j = 0;j < max_nb_value;j++) {
        frequency[j] = 0;
      }

      for (j = i;j >= 0;j--) {
        for (k = 1;k < nb_variable;k++) {
          frequency[int_sequence[index][k][j]]++;
        }

        contrast[j] = 0.;
        for (k = 0;k < max_nb_value;k++) {
          if (frequency[k] > 0) {
            contrast[j] += frequency[k] * log((double)frequency[k] / (double)((nb_variable - 1) * (i - j + 1)));
          }
        }
      }
    }

    else if (model_type[0] == MULTIVARIATE_POISSON_CHANGE) {
      for (j = 1;j < nb_variable;j++) {
        factorial[j][i] = 0.;
        for (k = 2;k <= int_sequence[index][j][i];k++) {
          factorial[j][i] += log((double)k);
        }
      }

      sum = 0.;
      factorial_sum = 0.;
      for (j = i;j >= 0;j--) {
        for (k = 1;k < nb_variable;k++) {
          sum += int_sequence[index][k][j];
          factorial_sum += factorial[k][j];
        }
        if (sum > 0.) {
          contrast[j] = sum * (log(sum / ((nb_variable - 1) * (i - j + 1))) - 1) - factorial_sum;
        }
        else {
          contrast[j] = 0.;
        }
      }
    }

    else if (model_type[0] == MULTIVARIATE_GEOMETRIC_0_CHANGE) {
      sum = 0.;
      for (j = i;j >= 0;j--) {
        for (k = 1;k < nb_variable;k++) {
          sum += int_sequence[index][k][j];
        }
        if (sum > 0.) {
          proba = (nb_variable - 1) * (i - j + 1) / ((nb_variable - 1) * (i - j + 1) + sum);
          contrast[j] = (nb_variable - 1) * (i - j + 1) * log(proba) + sum * log(1. - proba);
        }
        else {
          contrast[j] = 0.;
        }
      }
    }

    else {
      for (j = 0;j <= i;j++) {
        contrast[j] = 0.;
      }

      for (j = 1;j < nb_variable;j++) {
        if (model_type[j - 1] == CATEGORICAL_CHANGE) {
          for (k = 0;k < marginal_distribution[j]->nb_value;k++) {
            frequency[k] = 0;
          }
          sum = 0.;

          frequency[int_sequence[index][j][i]]++;
          for (k = i - 1;k >= 0;k--) {
            sum += (i - k) * log((double)(i - k) / (double)(i - k + 1)) +
                   log((double)(frequency[int_sequence[index][j][k]] + 1) / (double)(i - k + 1));
            if (frequency[int_sequence[index][j][k]] > 0) {
              sum -= frequency[int_sequence[index][j][k]] *
                     log((double)frequency[int_sequence[index][j][k]] / (double)(frequency[int_sequence[index][j][k]] + 1));
            }
            frequency[int_sequence[index][j][k]]++;

            if (contrast[k] != D_INF) {
              contrast[k] += sum;
            }

/*            frequency[int_sequence[index][j][k]]++;
            if (contrast[k] != D_INF) {
              for (m = 0;m < marginal_distribution[j]->nb_value;m++) {
                if (frequency[m] > 0) {
                  contrast[k] += frequency[m] * log((double)frequency[m] / (double)(i - k + 1));
                }
              }
            } */
          }
        }

        else if (model_type[j - 1] == POISSON_CHANGE) {
          factorial[j][i] = 0.;
          for (k = 2;k <= int_sequence[index][j][i];k++) {
            factorial[j][i] += log((double)k);
          }

          sum = 0.;
          factorial_sum = 0.;
          for (k = i;k >= 0;k--) {
            sum += int_sequence[index][j][k];
            factorial_sum += factorial[j][k];
            if ((contrast[k] != D_INF) && (sum > 0.)) {
              contrast[k] += sum * (log(sum / (i - k + 1)) - 1) - factorial_sum;
            }
          }
        }

        else if (model_type[j - 1] == GEOMETRIC_0_CHANGE) {
          sum = 0.;
          for (k = i;k >= 0;k--) {
            sum += int_sequence[index][j][k];
            if ((contrast[k] != D_INF) && (sum > 0.)) {
              proba = (i - k + 1) / (i - k + 1 + sum);
              contrast[k] += (i - k + 1) * log(proba) + sum * log(1. - proba);
            }
          }
        }

        else if (model_type[j - 1] == GEOMETRIC_1_CHANGE) {
          sum = 0.;
          for (k = i;k >= 0;k--) {
            sum += int_sequence[index][j][k];
            if ((contrast[k] != D_INF) && (sum > i - k + 1)) {
              proba = (i - k + 1) / sum;
              contrast[k] += (i - k + 1) * log(proba) + (sum - (i - k + 1)) * log(1. - proba);
            }
          }
        }

        else if (model_type[j - 1] == BAYESIAN_POISSON_CHANGE) {
          prior_contrast = -lgamma(hyperparam[j][0]) + hyperparam[j][0] * log(hyperparam[j][1]);

          factorial[j][i] = 0.;
          for (k = 2;k <= int_sequence[index][j][i];k++) {
            factorial[j][i] += log((double)k);
          }

          sum = 0.;
          factorial_sum = 0.;
          for (k = i;k >= 0;k--) {
            sum += int_sequence[index][j][k];
            factorial_sum += factorial[j][k];
            if (contrast[k] != D_INF) {
              contrast[k] += prior_contrast - factorial_sum + lgamma(hyperparam[j][0] + sum) -
                             (hyperparam[j][0] + sum) * log(hyperparam[j][1] + i - k + 1);
            }
          }
        }

        else if (model_type[j - 1] == BAYESIAN_GAUSSIAN_CHANGE) {
          prior_contrast = log(hyperparam[j][1]) / 2 - lgamma(hyperparam[j][2] / 2) +
                           hyperparam[j][2] * log(hyperparam[j][3] / 2) / 2;

          if (type[j] != REAL_VALUE) {
            square_sum = 0.;
            sum = int_sequence[index][j][i];
            if (contrast[i] != D_INF) {
              diff = hyperparam[j][0] - sum;
              contrast[i] += prior_contrast - log(2 * M_PI) / 2 -
                             log(hyperparam[j][1] + 1) / 2 + lgamma((hyperparam[j][2] + 1) / 2) -
                             (hyperparam[j][2] + 1) *
                             log((hyperparam[j][3] + hyperparam[j][1] *
                                  diff * diff / (hyperparam[j][1] + 1)) / 2) / 2;
            }

            for (k = i - 1;k >= 0;k--) {
              diff = int_sequence[index][j][k] - sum / (i - k);
              square_sum += ((double)(i - k) / (double)(i - k + 1)) * diff * diff;
              sum += int_sequence[index][j][k];
              if (contrast[k] != D_INF) {
                diff = hyperparam[j][0] - sum / (i - k + 1);
                contrast[k] += prior_contrast - (i - k + 1) * log(2 * M_PI) / 2 -
                               log(hyperparam[j][1] + i - k + 1) / 2 +
                               lgamma((hyperparam[j][2] + i - k + 1) / 2) -
                               (hyperparam[j][2] + i - k + 1) *
                               logl((hyperparam[j][3] + square_sum + hyperparam[j][1] * (i - k + 1) *
                                     diff * diff / (hyperparam[j][1] + i - k + 1)) / 2) / 2;
              }
            }
          }

          else {
            square_sum = 0.;
            sum = real_sequence[index][j][i];
            if (contrast[i] != D_INF) {
              diff = hyperparam[j][0] - sum;
              contrast[i] += prior_contrast - log(2 * M_PI) / 2 -
                             log(hyperparam[j][1] + 1) / 2 + lgamma((hyperparam[j][2] + 1) / 2) -
                             (hyperparam[j][2] + 1) *
                             log((hyperparam[j][3] + hyperparam[j][1] *
                                  diff * diff / (hyperparam[j][1] + 1)) / 2) / 2;
            }

            for (k = i - 1;k >= 0;k--) {
              diff = real_sequence[index][j][k] - sum / (i - k);
              square_sum += ((double)(i - k) / (double)(i - k + 1)) * diff * diff;
              sum += real_sequence[index][j][k];
              if (contrast[k] != D_INF) {
                diff = hyperparam[j][0] - sum / (i - k + 1);
                contrast[k] += prior_contrast - (i - k + 1) * log(2 * M_PI) / 2 -
                               log(hyperparam[j][1] + i - k + 1) / 2 +
                               lgamma((hyperparam[j][2] + i - k + 1) / 2) -
                               (hyperparam[j][2] + i - k + 1) *
                               logl((hyperparam[j][3] + square_sum + hyperparam[j][1] * (i - k + 1) *
                                     diff * diff / (hyperparam[j][1] + i - k + 1)) / 2) / 2;
              }
            }
          }
        }

        else {
          if (model_type[j - 1] == VARIANCE_CHANGE) {
            square_sum = 0.;

            if (type[j] != REAL_VALUE) {
              for (k = i;k >= 0;k--) {
                diff = int_sequence[index][j][k] - sequence_mean[j];
                square_sum += diff * diff;
                residual[k] = square_sum;
              }
            }

            else {
              for (k = i;k >= 0;k--) {
                diff = real_sequence[index][j][k] - sequence_mean[j];
                square_sum += diff * diff;
                residual[k] = square_sum;
              }
            }
          }

          else if (model_type[j - 1] == ORDINAL_GAUSSIAN_CHANGE) {
            square_sum = 0.;
            sum = rank[j][int_sequence[index][j][i]];
            residual[i] = 0.;

            for (k = i - 1;k >= 0;k--) {
              diff = rank[j][int_sequence[index][j][k]] - sum / (i - k);
              square_sum += ((double)(i - k) / (double)(i - k + 1)) * diff * diff;
              sum += rank[j][int_sequence[index][j][k]];
              residual[k] = square_sum;
            }
          }

          else if ((model_type[j - 1] == LINEAR_MODEL_CHANGE) || (model_type[0] == INTERCEPT_SLOPE_CHANGE)) {
            if (type[j] != REAL_VALUE) {
              square_sum = 0.;
              index_parameter_square_sum = 0.;
              mix_square_sum = 0.;
              sum = int_sequence[index][j][i];
              index_parameter_sum = seq_index_parameter[i];
              residual[i] = 0.;

              for (k = i - 1;k >= 0;k--) {
                diff = int_sequence[index][j][k] - sum / (i - k);
                square_sum += ((double)(i - k) / (double)(i - k + 1)) * diff * diff;
                index_parameter_diff = seq_index_parameter[k] - index_parameter_sum / (i - k);
                index_parameter_square_sum += ((double)(i - k) / (double)(i - k + 1)) *
                                              index_parameter_diff * index_parameter_diff;
                mix_square_sum += ((double)(i - k) / (double)(i - k + 1)) * diff * index_parameter_diff;
                sum += int_sequence[index][j][k];
                index_parameter_sum += seq_index_parameter[k];

                if ((k < i - 1) && (index_parameter_square_sum > 0.)) {
                  residual[k] = square_sum - mix_square_sum * mix_square_sum / index_parameter_square_sum;
                }
                else {
                  residual[k] = 0.;
                }
              }
            }

            else {
              square_sum = 0.;
              index_parameter_square_sum = 0.;
              mix_square_sum = 0.;
              sum = real_sequence[index][j][i];
              index_parameter_sum = seq_index_parameter[i];
              residual[i] = 0.;

              for (k = i - 1;k >= 0;k--) {
                diff = real_sequence[index][j][k] - sum / (i - k);
                square_sum += ((double)(i - k) / (double)(i - k + 1)) * diff * diff;
                index_parameter_diff = seq_index_parameter[k] - index_parameter_sum / (i - k);
                index_parameter_square_sum += ((double)(i - k) / (double)(i - k + 1)) *
                                              index_parameter_diff * index_parameter_diff;
                mix_square_sum += ((double)(i - k) / (double)(i - k + 1)) * diff *  index_parameter_diff;
                sum += real_sequence[index][j][k];
                index_parameter_sum += seq_index_parameter[k];

                if ((k < i - 1) && (index_parameter_square_sum > 0.)) {
                  residual[k] = square_sum - mix_square_sum * mix_square_sum / index_parameter_square_sum;
                }
                else {
                  residual[k] = 0.;
                }
              }
            }
          }

          else {
            if (type[j] != REAL_VALUE) {
              square_sum = 0.;
              sum = int_sequence[index][j][i];
              residual[i] = 0.;

              for (k = i - 1;k >= 0;k--) {
                diff = int_sequence[index][j][k] - sum / (i - k);
                square_sum += ((double)(i - k) / (double)(i - k + 1)) * diff * diff;
                sum += int_sequence[index][j][k];
                residual[k] = square_sum;
              }
            }

            else {
              square_sum = 0.;
              sum = real_sequence[index][j][i];
              residual[i] = 0.;

              for (k = i - 1;k >= 0;k--) {
                diff = real_sequence[index][j][k] - sum / (i - k);
                square_sum += ((double)(i - k) / (double)(i - k + 1)) * diff * diff;
                sum += real_sequence[index][j][k];
                residual[k] = square_sum;
              }
            }
          }

          if ((model_type[0] == MEAN_CHANGE) || (model_type[0] == INTERCEPT_SLOPE_CHANGE)) {
            for (k = i - 1;k >= 0;k--) {
              contrast[k] -= residual[k];
            }
          }

          else if (model_type[0] == MEAN_VARIANCE_CHANGE) {
            for (k = i - 1;k >= 0;k--) {
              contrast[k] += residual[k];
            }
          }

          else {
            for (k = i;k >= 0;k--) {
//              if ((contrast[k] != D_INF) && (residual[k] > 0.)) {
              if ((contrast[k] != D_INF) && (residual[k] > sqrt((double)(i - k + 1)) * ROUNDOFF_ERROR)) {
                contrast[k] -= ((double)(i - k + 1) / 2.) * (logl(residual[k] /
                                 (i - k + 1)) + log(2 * M_PI) + 1);
/*                contrast[k] -= ((double)(i - k + 1) / 2.) * (logl(residual[k] /
                                 (i - k)) + log(2 * M_PI)) + (double)(i - k) / 2.; */
              }
              else {
                contrast[k] = D_INF;
              }
            }
          }
        }
      }

      if (model_type[0] == MEAN_VARIANCE_CHANGE) {
        contrast[i] = D_INF;
        for (j = i - 1;j >= 0;j--) {
//          if (contrast[j] > 0.) {
          if (contrast[j] > sqrt((double)((nb_variable - 1) * (i - j + 1))) * ROUNDOFF_ERROR) {
            contrast[j] = -((double)((nb_variable - 1) * (i - j + 1)) / 2.) * (logl(contrast[j] /
                             ((nb_variable - 1) * (i - j + 1))) + log(2 * M_PI) + 1);
/*            contrast[j] = -((double)((nb_variable - 1) * (i - j + 1)) / 2.) * (logl(contrast[j] /
                             ((nb_variable - 1) * (i - j))) + log(2 * M_PI)) +
                           (double)((nb_variable - 1) * (i - j)) / 2.; */
          }
          else {
            contrast[j] = D_INF;
          }
        }
      }
    }

    for (j = 0;j < nb_segment;j++) {
      forward[i][j] = D_INF;
    }

    for (j = MAX(0 , nb_segment + i - length[index]);j < MIN((i < length[index] - 1 ? nb_segment - 1 : nb_segment) , i + 1);j++) {
      if (j == 0) {
        forward[i][j] = contrast[0];
        if (forward[i][j] != D_INF) {
          optimal_length[i][j] = i + 1;
        }
      }

      else {
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

  if (forward[length[index] - 1][nb_segment - 1] == D_INF) {
    segmentation_likelihood = D_INF;
  }

  else {

    // restauration

    i = length[index] - 1;
    psegment = int_sequence[index][0] + i;

    for (j = nb_segment - 1;j >= 0;j--) {
      for (k = i;k > i - optimal_length[i][j];k--) {
        *psegment-- = j;
      }
      i -= optimal_length[i][j];
    }

    // recurrence "backward"

    for (i = length[index] - 1;i >= 0;i--) {

      // calcul des log-vraisemblances des segments

      if (model_type[0] == MULTIVARIATE_CATEGORICAL_CHANGE) {
        for (j = 0;j < max_nb_value;j++) {
          frequency[j] = 0;
        }

        for (j = i;j < length[index];j++) {
          for (k = 1;k < nb_variable;k++) {
            frequency[int_sequence[index][k][j]]++;
          }

          contrast[j] = 0.;
          for (k = 0;k < max_nb_value;k++) {
            if (frequency[k] > 0) {
              contrast[j] += frequency[k] * log((double)frequency[k] / (double)((nb_variable - 1) * (j - i + 1)));
            }
          }
        }
      }

      else if (model_type[0] == MULTIVARIATE_POISSON_CHANGE) {
        sum = 0.;
        factorial_sum = 0.;
        for (j = i;j < length[index];j++) {
          for (k = 1;k < nb_variable;k++) {
            sum += int_sequence[index][k][j];
            factorial_sum += factorial[k][j];
          }
          if (sum > 0.) {
            contrast[j] = sum * (log(sum / ((nb_variable - 1) * (j - i + 1))) - 1) - factorial_sum;
          }
          else {
            contrast[j] = 0.;
          }
        }
      }

      else if (model_type[0] == MULTIVARIATE_GEOMETRIC_0_CHANGE) {
        sum = 0.;
        for (j = i;j < length[index];j++) {
          for (k = 1;k < nb_variable;k++) {
            sum += int_sequence[index][k][j];
          }
          if (sum > 0.) {
            proba = (nb_variable - 1) * (j - i + 1) / ((nb_variable - 1) * (j - i + 1) + sum);
            contrast[j] = (nb_variable - 1) * (j - i + 1) * log(proba) + sum * log(1. - proba);
          }
          else {
            contrast[j] = 0.;
          }
        }
      }

      else {
        for (j = i;j < length[index];j++) {
          contrast[j] = 0.;
        }

        for (j = 1;j < nb_variable;j++) {
          if (model_type[j - 1] == CATEGORICAL_CHANGE) {
            for (k = 0;k < marginal_distribution[j]->nb_value;k++) {
              frequency[k] = 0;
            }
            sum = 0.;

            frequency[int_sequence[index][j][i]]++;
            for (k = i + 1;k < length[index];k++) {
              sum += (k - i) * log((double)(k - i) / (double)(k - i + 1)) +
                     log((double)(frequency[int_sequence[index][j][k]] + 1) / (double)(k - i + 1));
              if (frequency[int_sequence[index][j][k]] > 0) {
                sum -= frequency[int_sequence[index][j][k]] *
                       log((double)frequency[int_sequence[index][j][k]] / (double)(frequency[int_sequence[index][j][k]] + 1));
              }
              frequency[int_sequence[index][j][k]]++;

              if (contrast[k] != D_INF) {
                contrast[k] += sum;
              }

/*              frequency[int_sequence[index][j][k]]++;
              if (contrast[k] != D_INF) {
                for (m = 0;m < marginal_distribution[j]->nb_value;m++) {
                  if (frequency[m] > 0) {
                    contrast[k] += frequency[m] * log((double)frequency[m] / (double)(k - i + 1));
                  }
                }
              } */
            }
          }

          else if (model_type[j - 1] == POISSON_CHANGE) {
            sum = 0.;
            factorial_sum = 0.;
            for (k = i;k < length[index];k++) {
              sum += int_sequence[index][j][k];
              factorial_sum += factorial[j][k];
              if ((contrast[k] != D_INF) && (sum > 0.)) {
                contrast[k] += sum * (log(sum / (k - i + 1)) - 1) - factorial_sum;
              }
            }
          }

          else if (model_type[j - 1] == GEOMETRIC_0_CHANGE) {
            sum = 0.;
            for (k = i;k < length[index];k++) {
              sum += int_sequence[index][j][k];
              if ((contrast[k] != D_INF) && (sum > 0.)) {
                proba = (k - i + 1) / (k - i + 1 + sum);
                contrast[k] += (k - i + 1) * log(proba) + sum * log(1. - proba);
              }
            }
          }

          else if (model_type[j - 1] == GEOMETRIC_1_CHANGE) {
            sum = 0.;
            for (k = i;k < length[index];k++) {
              sum += int_sequence[index][j][k];
              if ((contrast[k] != D_INF) && (sum > k - i + 1)) {
                proba = (k - i + 1) / sum;
                contrast[k] += (k - i + 1) * log(proba) + (sum - (k - i + 1)) * log(1. - proba);
              }
            }
          }

          else if (model_type[j - 1] == BAYESIAN_POISSON_CHANGE) {
            prior_contrast = -lgamma(hyperparam[j][0]) + hyperparam[j][0] * log(hyperparam[j][1]);

            sum = 0.;
            factorial_sum = 0.;
            for (k = i;k < length[index];k++) {
              sum += int_sequence[index][j][k];
              factorial_sum += factorial[j][k];
              if (contrast[k] != D_INF) {
                contrast[k] += prior_contrast - factorial_sum + lgamma(hyperparam[j][0] + sum) -
                               (hyperparam[j][0] + sum) * log(hyperparam[j][1] + k - i + 1);
              }
            }
          }

          else if (model_type[j - 1] == BAYESIAN_GAUSSIAN_CHANGE) {
            prior_contrast = log(hyperparam[j][1]) / 2 - lgamma(hyperparam[j][2] / 2) +
                             hyperparam[j][2] * log(hyperparam[j][3] / 2) / 2;

            if (type[j] != REAL_VALUE) {
              square_sum = 0.;
              sum = int_sequence[index][j][i];
              if (contrast[i] != D_INF) {
                diff = hyperparam[j][0] - sum;
                contrast[i] += prior_contrast - log(2 * M_PI) / 2 -
                               log(hyperparam[j][1] + 1) / 2 + lgamma((hyperparam[j][2] + 1) / 2) -
                               (hyperparam[j][2] + 1) *
                               log((hyperparam[j][3] + hyperparam[j][1] *
                                    diff * diff / (hyperparam[j][1] + 1)) / 2) / 2;
              }

              for (k = i + 1;k < length[index];k++) {
                diff = int_sequence[index][j][k] - sum / (k - i);
                square_sum += ((double)(k - i) / (double)(k - i + 1)) * diff * diff;
                sum += int_sequence[index][j][k];
                if (contrast[k] != D_INF) {
                  diff = hyperparam[j][0] - sum / (k - i + 1);
                  contrast[k] += prior_contrast - (k - i + 1) * log(2 * M_PI) / 2 -
                                 log(hyperparam[j][1] + k - i + 1) / 2 +
                                 lgamma((hyperparam[j][2] + k - i + 1) / 2) -
                                 (hyperparam[j][2] + k - i + 1) *
                                 logl((hyperparam[j][3] + square_sum + hyperparam[j][1] * (k - i + 1) *
                                       diff * diff / (hyperparam[j][1] + k - i + 1)) / 2) / 2;
                }
              }
            }

            else {
              square_sum = 0.;
              sum = real_sequence[index][j][i];
              if (contrast[i] != D_INF) {
                diff = hyperparam[j][0] - sum;
                contrast[i] += prior_contrast - log(2 * M_PI) / 2 -
                               log(hyperparam[j][1] + 1) / 2 + lgamma((hyperparam[j][2] + 1) / 2) -
                               (hyperparam[j][2] + 1) *
                               log((hyperparam[j][3] + hyperparam[j][1] *
                                    diff * diff / (hyperparam[j][1] + 1)) / 2) / 2;
              }

              for (k = i + 1;k < length[index];k++) {
                diff = real_sequence[index][j][k] - sum / (k - i);
                square_sum += ((double)(k - i) / (double)(k - i + 1)) * diff * diff;
                sum += real_sequence[index][j][k];
                if (contrast[k] != D_INF) {
                  diff = hyperparam[j][0] - sum / (k - i + 1);
                  contrast[k] += prior_contrast - (k - i + 1) * log(2 * M_PI) / 2 -
                                 log(hyperparam[j][1] + k - i + 1) / 2 +
                                 lgamma((hyperparam[j][2] + k - i + 1) / 2) -
                                 (hyperparam[j][2] + k - i + 1) *
                                 logl((hyperparam[j][3] + square_sum + hyperparam[j][1] * (k - i + 1) *
                                       diff * diff / (hyperparam[j][1] + k - i + 1)) / 2) / 2;
                }
              }
            }
          }

          else {
            if (model_type[j - 1] == VARIANCE_CHANGE) {
              square_sum = 0.;

              if (type[j] != REAL_VALUE) {
                for (k = i;k < length[index];k++) {
                  diff = int_sequence[index][j][k] - sequence_mean[j];
                  square_sum += diff * diff;
                  residual[k] = square_sum;
                }
              }

              else {
                for (k = i;k < length[index];k++) {
                  diff = real_sequence[index][j][k] - sequence_mean[j];
                  square_sum += diff * diff;
                  residual[k] = square_sum;
                }
              }
            }

            else if (model_type[j - 1] == ORDINAL_GAUSSIAN_CHANGE) {
              square_sum = 0.;
              sum = rank[j][int_sequence[index][j][i]];
              residual[i] = 0.;

              for (k = i + 1;k < length[index];k++) {
                diff = rank[j][int_sequence[index][j][k]] - sum / (k - i);
                square_sum += ((double)(k - i) / (double)(k - i + 1)) * diff * diff;
                sum += rank[j][int_sequence[index][j][k]];
                residual[k] = square_sum;
              }
            }

            else if ((model_type[j - 1] == LINEAR_MODEL_CHANGE) || (model_type[0] == INTERCEPT_SLOPE_CHANGE)) {
              if (type[j] != REAL_VALUE) {
                square_sum = 0.;
                index_parameter_square_sum = 0.;
                mix_square_sum = 0.;
                sum = int_sequence[index][j][i];
                index_parameter_sum = seq_index_parameter[i];
                residual[i] = 0.;

                for (k = i + 1;k < length[index];k++) {
                  diff = int_sequence[index][j][k] - sum / (k - i);
                  square_sum += ((double)(k - i) / (double)(k - i + 1)) * diff * diff;
                  index_parameter_diff = seq_index_parameter[k] - index_parameter_sum / (k - i);
                  index_parameter_square_sum += ((double)(k - i) / (double)(k - i + 1)) *
                                                index_parameter_diff * index_parameter_diff;
                  mix_square_sum += ((double)(k - i) / (double)(k - i + 1)) * diff * index_parameter_diff;
                  sum += int_sequence[index][j][k];
                  index_parameter_sum += seq_index_parameter[k];

                  if ((k > i + 1) && (index_parameter_square_sum > 0.)) {
                    residual[k] = square_sum - mix_square_sum * mix_square_sum / index_parameter_square_sum;
                  }
                  else {
                    residual[k] = 0.;
                  }
                }
              }

              else {
                square_sum = 0.;
                index_parameter_square_sum = 0.;
                mix_square_sum = 0.;
                sum = real_sequence[index][j][i];
                index_parameter_sum = seq_index_parameter[i];
                residual[i] = 0.;

                for (k = i + 1;k < length[index];k++) {
                  diff = real_sequence[index][j][k] - sum / (k - i);
                  square_sum += ((double)(k - i) / (double)(k - i + 1)) * diff * diff;
                  index_parameter_diff = seq_index_parameter[k] - index_parameter_sum / (k - i);
                  index_parameter_square_sum += ((double)(k - i) / (double)(k - i + 1)) *
                                                index_parameter_diff * index_parameter_diff;
                  mix_square_sum += ((double)(k - i) / (double)(k - i + 1)) * diff * index_parameter_diff;
                  sum += real_sequence[index][j][k];
                  index_parameter_sum += seq_index_parameter[k];

                  if ((k > i + 1) && (index_parameter_square_sum > 0.)) {
                    residual[k] = square_sum - mix_square_sum * mix_square_sum / index_parameter_square_sum;
                  }
                  else {
                    residual[k] = 0.;
                  }
                }
              }
            }

            else {
              if (type[j] != REAL_VALUE) {
                square_sum = 0.;
                sum = int_sequence[index][j][i];
                residual[i] = 0.;

                for (k = i + 1;k < length[index];k++) {
                  diff = int_sequence[index][j][k] - sum / (k - i);
                  square_sum += ((double)(k - i) / (double)(k - i + 1)) * diff * diff;
                  sum += int_sequence[index][j][k];
                  residual[k] = square_sum;
                }
              }

              else {
                square_sum = 0.;
                sum = real_sequence[index][j][i];
                residual[i] = 0.;

                for (k = i + 1;k < length[index];k++) {
                  diff = real_sequence[index][j][k] - sum / (k - i);
                  square_sum += ((double)(k - i) / (double)(k - i + 1)) * diff * diff;
                  sum += real_sequence[index][j][k];
                  residual[k] = square_sum;
                }
              }
            }

            if ((model_type[0] == MEAN_CHANGE) || (model_type[0] == INTERCEPT_SLOPE_CHANGE)) {
              for (k = i + 1;k < length[index];k++) {
                contrast[k] -= residual[k];
              }
            }

            else if (model_type[0] == MEAN_VARIANCE_CHANGE) {
              for (k = i + 1;k < length[index];k++) {
                contrast[k] += residual[k];
              }
            }

            else {
              for (k = i;k < length[index];k++) {
//                if ((contrast[k] != D_INF) && (residual[k] > 0.)) {
                if ((contrast[k] != D_INF) && (residual[k] > sqrt((double)(k - i + 1)) * ROUNDOFF_ERROR)) {
                  contrast[k] -= ((double)(k - i + 1) / 2.) * (logl(residual[k] /
                                   (k - i + 1)) + log(2 * M_PI) + 1);
/*                  contrast[k] -= ((double)(k - i + 1) / 2.) * (logl(residual[k] /
                                   (k - i)) + log(2 * M_PI)) + (double)(k - i) / 2.; */
                }
                else {
                  contrast[k] = D_INF;
                }
              }
            }
          }
        }

        if (model_type[0] == MEAN_VARIANCE_CHANGE) {
          contrast[i] = D_INF;
          for (j = i + 1;j < length[index];j++) {
//            if (contrast[j] > 0.) {
            if (contrast[j] > sqrt((double)((nb_variable - 1) * (j - i + 1))) * ROUNDOFF_ERROR) {
              contrast[j] = -((double)((nb_variable - 1) * (j - i + 1)) / 2.) * (logl(contrast[j] /
                               ((nb_variable - 1) * (j - i + 1))) + log(2 * M_PI) + 1);
/*              contrast[j] = -((double)((nb_variable - 1) * (j - i + 1)) / 2.) * (logl(contrast[j] /
                               ((nb_variable - 1) * (j - i))) + log(2 * M_PI)) +
                             (double)((nb_variable - 1) * (j - i)) / 2.; */
            }
            else {
              contrast[j] = D_INF;
            }
          }
        }
      }

      for (j = 0;j < nb_segment;j++) {
        backward_output[i][j] = D_INF;
      }

      for (j = MAX((i == 0 ? 0 : 1) , nb_segment + i - length[index]);j < MIN(nb_segment , i + 1);j++) {
        if (j < nb_segment - 1) {
          backward[i][j] = D_INF;
          for (k = length[index] + j - nb_segment;k >= i;k--) {
            if ((contrast[k] != D_INF) && (backward[k + 1][j + 1] != D_INF)) {
              buff = contrast[k] + backward[k + 1][j + 1];
              if (buff > backward[i][j]) {
                backward[i][j] = buff;
              }
            }

            if ((output == SEGMENT) && (k > i) && (backward[i][j] != D_INF)) {
              if (i == 0) {
                if (backward[i][j] > backward_output[k][j]) {
                  backward_output[k][j] = backward[i][j];
                }
              }
              else if (forward[i - 1][j - 1] != D_INF) {
                buff = forward[i - 1][j - 1] + backward[i][j];
                if (buff > backward_output[k][j]) {
                  backward_output[k][j] = buff;
                }
              }
            }
          }
        }

        else {
          backward[i][j] = contrast[length[index] - 1];

          if ((output == SEGMENT) && (forward[i - 1][j - 1] != D_INF) &&
              (backward[i][j] != D_INF)) {
            buff = forward[i - 1][j - 1] + backward[i][j];
            for (k = length[index] - 1;k > i;k--) {
              if (buff > backward_output[k][j]) {
                backward_output[k][j] = buff;
              }
            }
          }
        }

        if (backward[i][j] != D_INF) {
          if (i == 0) {
            backward_output[i][j] = backward[i][j];
          }
          else if (forward[i - 1][j - 1] != D_INF) {
            backward_output[i][j] = forward[i - 1][j - 1] + backward[i][j];
          }
        }
      }
    }

#   ifdef DEBUG
    cout << "\n";
    for (i = 1;i < length[index];i++) {
      cout << i;
      for (j = 0;j < nb_segment;j++) {
        if (j == 0) {
          cout << " | " << backward[i][j];
        }
        else {
          cout << " | " << ((forward[i - 1][j - 1] != D_INF) && (backward[i][j] != D_INF) ? forward[i - 1][j - 1] + backward[i][j] : D_INF);
        }
        cout << " " << backward_output[i][j];
      }
      cout << endl;
    }
    cout << endl;
#   endif

    // restauration

#   ifdef MESSAGE
    if (output == SEGMENT) {
      int optimal_segment;

      psegment = int_sequence[index][0];

      for (i = 0;i < length[index];i++) {
        backward_max = D_INF;
        for (j = 0;j < nb_segment;j++) {
          if (backward_output[i][j] > backward_max) {
            backward_max = backward_output[i][j];
            optimal_segment = j;
          }
        }

        if (optimal_segment != *psegment) {
          cout << "\nERROR: " << i << " | " << *psegment << " " << optimal_segment << endl;
        }

        psegment++;
      }
    }
#   endif

    for (i = 1;i < nb_variable;i++) {
      if ((model_type[i - 1] == POISSON_CHANGE) || (model_type[0] == MULTIVARIATE_POISSON_CHANGE) ||
          (model_type[i - 1] == GEOMETRIC_0_CHANGE) || (model_type[i - 1] == GEOMETRIC_1_CHANGE) ||
          (model_type[0] == MULTIVARIATE_GEOMETRIC_0_CHANGE) || (model_type[i - 1] == GAUSSIAN_CHANGE) ||
          (model_type[i - 1] == VARIANCE_CHANGE) || (model_type[0] == MEAN_CHANGE) ||
          (model_type[0] == MEAN_VARIANCE_CHANGE) || (model_type[i - 1] == BAYESIAN_POISSON_CHANGE) ||
          (model_type[i - 1] == BAYESIAN_GAUSSIAN_CHANGE)) {
        psegment = int_sequence[index][0] + 1;

        if (type[i] != REAL_VALUE) {
          piecewise_function[i][0] = int_sequence[index][i][0];
          j = 0;

          for (k = 1;k < length[index];k++) {
            if (*psegment != *(psegment - 1)) {
              piecewise_function[i][j] /= (k - j);
              for (m = j + 1;m < k;m++) {
                piecewise_function[i][m] = piecewise_function[i][j];
              }
              j = k;
              piecewise_function[i][j] = int_sequence[index][i][k];
            }
            else {
              piecewise_function[i][j] += int_sequence[index][i][k];
            }
            psegment++;
          }
        }

        else {
          piecewise_function[i][0] = real_sequence[index][i][0];
          j = 0;

          for (k = 1;k < length[index];k++) {
            if (*psegment != *(psegment - 1)) {
              piecewise_function[i][j] /= (k - j);
              for (m = j + 1;m < k;m++) {
                piecewise_function[i][m] = piecewise_function[i][j];
              }
              j = k;
              piecewise_function[i][j] = real_sequence[index][i][k];
            }
            else {
              piecewise_function[i][j] += real_sequence[index][i][k];
            }
            psegment++;
          }
        }

        piecewise_function[i][j] /= (length[index] - j);
        for (k = j + 1;k < length[index];k++) {
          piecewise_function[i][k] = piecewise_function[i][j];
        }
      }

      else if ((model_type[i - 1] == LINEAR_MODEL_CHANGE) || (model_type[0] == INTERCEPT_SLOPE_CHANGE)) {
        psegment = int_sequence[index][0] + 1;

        if (type[i] != REAL_VALUE) {
          response_mean = int_sequence[index][i][0];
          index_parameter_mean = seq_index_parameter[0];                
          j = 0;

          for (k = 1;k < length[index];k++) {
            if (*psegment != *(psegment - 1)) {
              response_mean /= (k - j);
              index_parameter_mean /= (k - j);

              index_parameter_variance = 0.;
              covariance = 0.;
              for (m = j;m < k;m++) {
                diff = seq_index_parameter[m] - index_parameter_mean;
                index_parameter_variance += diff * diff;                
                covariance += (int_sequence[index][i][m] - response_mean) * (seq_index_parameter[m] - index_parameter_mean);
              }

              slope = covariance / index_parameter_variance;
              intercept = response_mean - slope * index_parameter_mean;

              for (m = j;m < k;m++) {
                piecewise_function[i][m] = intercept + slope * seq_index_parameter[m];
              }

              j = k;
              response_mean = int_sequence[index][i][k];
              index_parameter_mean = seq_index_parameter[k];                
            }

            else {
              response_mean += int_sequence[index][i][k];
              index_parameter_mean += seq_index_parameter[k];                
            }
            psegment++;
          }

          response_mean /= (length[index] - j);
          index_parameter_mean /= (length[index] - j);

          index_parameter_variance = 0.;
          covariance = 0.;
          for (k = j;k < length[index];k++) {
            diff = seq_index_parameter[k] - index_parameter_mean;
            index_parameter_variance += diff * diff;                
            covariance += (int_sequence[index][i][k] - response_mean) * (seq_index_parameter[k] - index_parameter_mean);
          }

          slope = covariance / index_parameter_variance;
          intercept = response_mean - slope * index_parameter_mean;

          for (k = j;k < length[index];k++) {
            piecewise_function[i][k] = intercept + slope * seq_index_parameter[k];
          }
        }

        else {
          response_mean = real_sequence[index][i][0];
          index_parameter_mean = seq_index_parameter[0];                
          j = 0;

          for (k = 1;k < length[index];k++) {
            if (*psegment != *(psegment - 1)) {
              response_mean /= (k - j);
              index_parameter_mean /= (k - j);

              index_parameter_variance = 0.;
              covariance = 0.;
              for (m = j;m < k;m++) {
                diff = seq_index_parameter[m] - index_parameter_mean;
                index_parameter_variance += diff * diff;                
                covariance += (real_sequence[index][i][m] - response_mean) * (seq_index_parameter[m] - index_parameter_mean);
              }

              slope = covariance / index_parameter_variance;
              intercept = response_mean - slope * index_parameter_mean;

              for (m = j;m < k;m++) {
                piecewise_function[i][m] = intercept + slope * seq_index_parameter[m];
              }

              j = k;
              response_mean = real_sequence[index][i][k];
              index_parameter_mean = seq_index_parameter[k];                
            }

            else {
              response_mean += real_sequence[index][i][k];
              index_parameter_mean += seq_index_parameter[k];                
            }
            psegment++;
          }

          response_mean /= (length[index] - j);
          index_parameter_mean /= (length[index] - j);

          index_parameter_variance = 0.;
          covariance = 0.;
          for (k = j;k < length[index];k++) {
            diff = seq_index_parameter[k] - index_parameter_mean;
            index_parameter_variance += diff * diff;                
            covariance += (real_sequence[index][i][k] - response_mean) * (seq_index_parameter[k] - index_parameter_mean);
          }

          slope = covariance / index_parameter_variance;
          intercept = response_mean - slope * index_parameter_mean;

          for (k = j;k < length[index];k++) {
            piecewise_function[i][k] = intercept + slope * seq_index_parameter[k];
          }
        }
      }
    }

#   ifdef MESSAGE
    if ((backward[0][0] < forward[length[index] - 1][nb_segment - 1] - DOUBLE_ERROR) ||
        (backward[0][0] > forward[length[index] - 1][nb_segment - 1] + DOUBLE_ERROR)) {
      cout << "\nERROR: " << backward[0][0] << " | " << forward[length[index] - 1][nb_segment - 1] << endl;
    }
/*    if ((backward_output[0][0] < backward[0][0] - DOUBLE_ERROR) ||
        (backward_output[0][0] > backward[0][0] + DOUBLE_ERROR)) {
      cout << "\nERROR: " << backward_output[0][0] << " | " << backward[0][0] << endl;
    } */
#   endif

    if ((model_type[0] != MEAN_CHANGE) && (model_type[0] != INTERCEPT_SLOPE_CHANGE)) {
      segmentation_likelihood = forward[length[index] - 1][nb_segment - 1];

      if (likelihood != D_INF) {
        for (i = 0;i < length[index];i++) {
          for (j = 0;j < nb_segment;j++) {
            if (backward_output[i][j] != D_INF) {
              backward_output[i][j] = exp(backward_output[i][j] - likelihood);
            }
            else {
              backward_output[i][j] = 0.;
            }
          }
        }
      }

      else if (segmentation_likelihood != D_INF) {
        for (i = 0;i < length[index];i++) {
          for (j = 0;j < nb_segment;j++) {
            if (backward_output[i][j] != D_INF) {
              backward_output[i][j] = exp(backward_output[i][j] - segmentation_likelihood);
            }
            else {
              backward_output[i][j] = 0.;
            }
          }
        }
      }
    }

    else {
      if (forward[length[index] - 1][nb_segment - 1] < 0.) {
        segmentation_likelihood = -((double)((nb_variable - 1) * length[index]) / 2.) *
                                   (log(-forward[length[index] - 1][nb_segment - 1] /
                                     ((nb_variable - 1) * length[index])) + log(2 * M_PI) + 1);
/*        segmentation_likelihood = -((double)((nb_variable - 1) * length[index]) / 2.) *
                                   (log(-forward[length[index] - 1][nb_segment - 1] /
                                     ((nb_variable - 1) * (length[index] - nb_segment))) + log(2 * M_PI)) -
                                   (double)((nb_variable - 1) * (length[index] - nb_segment)) / 2.; */

        for (i = 0;i < length[index];i++) {
          for (j = 0;j < nb_segment;j++) {
            if (backward_output[i][j] < 0.) {
              backward_output[i][j] = pow(backward_output[i][j] / forward[length[index] - 1][nb_segment - 1] ,
                                          -((double)((nb_variable - 1) * length[index]) / 2.));
/*              backward_output[i][j] = exp(-((double)((nb_variable - 1) * length[index]) / 2.) *
                                          log(backward_output[i][j] / forward[length[index] - 1][nb_segment - 1])); */
            }
            else {
              backward_output[i][j] = 0.;
            }
          }
        }
      }

      else {
        segmentation_likelihood = D_INF;
      }
    }

    if (segmentation_likelihood == D_INF) {
      for (i = 0;i < length[index];i++) {
        for (j = 0;j < nb_segment;j++) {
          backward_output[i][j] = 0.;
        }
      }
    }

    switch (format) {

    case 'a' : {
      if (likelihood != D_INF) {
        switch (output) {
        case CHANGE_POINT :
          *os << "\n" << SEQ_label[SEQL_MAX_POSTERIOR_CHANGE_POINT_PROBABILITY] << "\n\n";
          break;
        case SEGMENT :
          *os << "\n" << SEQ_label[SEQL_MAX_POSTERIOR_SEGMENT_PROBABILITY] << "\n\n";
          break;
        }
      }

      else {
        switch (output) {
        case CHANGE_POINT :
          *os << "\n" << SEQ_label[SEQL_MAX_CHANGE_POINT_LIKELIHOOD] << "\n\n";
          break;
        case SEGMENT :
          *os << "\n" << SEQ_label[SEQL_MAX_SEGMENT_LIKELIHOOD] << "\n\n";
          break;
        }
      }

      profile_ascii_print(*os , index , nb_segment , backward_output ,
                          (output == CHANGE_POINT ? SEQ_label[SEQL_CHANGE_POINT] : SEQ_label[SEQL_SEGMENT]) ,
                          piecewise_function);

      *os << "\n" << SEQ_label[SEQL_SEGMENTATION_LIKELIHOOD] << ": " << segmentation_likelihood;
      if (likelihood != D_INF) {
        *os << "   (" << exp(segmentation_likelihood - likelihood) << ")";
      }
      *os << endl;
      break;
    }

    case 's' : {
      if (likelihood != D_INF) {
        switch (output) {
        case CHANGE_POINT :
          *os << "\n" << SEQ_label[SEQL_MAX_POSTERIOR_CHANGE_POINT_PROBABILITY] << "\n\n";
          break;
        case SEGMENT :
          *os << "\n" << SEQ_label[SEQL_MAX_POSTERIOR_SEGMENT_PROBABILITY] << "\n\n";
          break;
        }
      }

      else {
        switch (output) {
        case CHANGE_POINT :
          *os << "\n" << SEQ_label[SEQL_MAX_CHANGE_POINT_LIKELIHOOD] << "\n\n";
          break;
        case SEGMENT :
          *os << "\n" << SEQ_label[SEQL_MAX_SEGMENT_LIKELIHOOD] << "\n\n";
          break;
        }
      }

      profile_spreadsheet_print(*os , index , nb_segment , backward_output ,
                                (output == CHANGE_POINT ? SEQ_label[SEQL_CHANGE_POINT] : SEQ_label[SEQL_SEGMENT]) ,
                                piecewise_function);

      *os << "\n" << SEQ_label[SEQL_SEGMENTATION_LIKELIHOOD] << "\t" << segmentation_likelihood;
      if (likelihood != D_INF) {
        *os << "\t" << exp(segmentation_likelihood - likelihood);
      }
      *os << endl;
      break;
    }

    case 'g' : {
      profile_plot_print(*os , index , nb_segment , backward_output , piecewise_function);
      break;
    }

    case 'p' : {
      MultiPlotSet &plot = *plot_set;

      i = 0;
      for (j = 1;j < nb_variable;j++) {
        if ((model_type[j - 1] == POISSON_CHANGE) || (model_type[0] == MULTIVARIATE_POISSON_CHANGE) ||
            (model_type[j - 1] == GEOMETRIC_0_CHANGE) || (model_type[j - 1] == GEOMETRIC_1_CHANGE) ||
            (model_type[0] == MULTIVARIATE_GEOMETRIC_0_CHANGE) || (model_type[j - 1] == GAUSSIAN_CHANGE) ||
            (model_type[j - 1] == VARIANCE_CHANGE) || (model_type[i - 1] == LINEAR_MODEL_CHANGE) ||
            (model_type[0] == MEAN_CHANGE) || (model_type[0] == MEAN_VARIANCE_CHANGE) ||
            (model_type[0] == INTERCEPT_SLOPE_CHANGE) || (model_type[j - 1] == BAYESIAN_POISSON_CHANGE) ||
            (model_type[j - 1] == BAYESIAN_GAUSSIAN_CHANGE)) {
          plot[i].resize(2);

          if (index_parameter) {
            if (type[j] != REAL_VALUE) {
              for (k = 0;k < length[index];k++) {
                plot[i][0].add_point(index_parameter[index][k] , int_sequence[index][j][k]);
              }
            }
            else {
              for (k = 0;k < length[index];k++) {
                plot[i][0].add_point(index_parameter[index][k] , real_sequence[index][j][k]);
              }
            }

            for (k = 0;k < length[index];k++) {
              plot[i][1].add_point(index_parameter[index][k] , piecewise_function[j][k]);
            }
          }

          else {
            if (type[j] != REAL_VALUE) {
              for (k = 0;k < length[index];k++) {
                plot[i][0].add_point(k , int_sequence[index][j][k]);
              }
            }
            else {
              for (k = 0;k < length[index];k++) {
                plot[i][0].add_point(k , real_sequence[index][j][k]);
              }
            }

            for (k = 0;k < length[index];k++) {
              plot[i][1].add_point(k , piecewise_function[j][k]);
            }
          }

          i++;
        }
      }

      profile_plotable_write(plot[i] , index , nb_segment , backward_output);
      break;
    }
    }

#   ifdef MESSAGE
    if (format != 'g') {
      double ambiguity = 0.;

      psegment = int_sequence[index][0];
      for (i = 0;i < length[index];i++) {
        for (j = 0;j < nb_segment;j++) {
          if (j != *psegment) {
            ambiguity += backward_output[i][j];
          }
        }
        psegment++;
      }

      if (likelihood != D_INF) {
        ambiguity *= exp(likelihood - segmentation_likelihood);
      }

      switch (format) {
      case 'a' :
        *os << "\n" << SEQ_label[SEQL_AMBIGUITY] << ": " << ambiguity
            << " (" << ambiguity / length[index] << ")" << endl;
        break;
      case 's' :
        *os << "\n" << SEQ_label[SEQL_AMBIGUITY] << "\t" << ambiguity
            << "\t" << ambiguity / length[index] << "\t" << endl;
        break;
      }
    }
#   endif

  }

  delete [] frequency;

  for (i = 1;i < nb_variable;i++) {
    delete [] factorial[i];
  }
  delete [] factorial;

  for (i = 1;i < nb_variable;i++) {
    delete [] hyperparam[i];
  }
  delete [] hyperparam;

  delete [] sequence_mean;
  delete [] residual;

  if (index_parameter_type == IMPLICIT_TYPE) {
    delete [] seq_index_parameter;
  }

  delete [] contrast;

  for (i = 0;i < length[index];i++) {
    delete [] forward[i];
  }
  delete [] forward;

  for (i = 0;i < length[index];i++) {
    delete [] optimal_length[i];
  }
  delete [] optimal_length;

  for (i = 0;i < length[index];i++) {
    delete [] backward[i];
  }
  delete [] backward;

  for (i = 0;i < length[index];i++) {
    delete [] backward_output[i];
  }
  delete [] backward_output;

  for (i = 1;i < nb_variable;i++) {
    delete [] piecewise_function[i];
  }
  delete [] piecewise_function;

  return segmentation_likelihood;
}


/*--------------------------------------------------------------*
 *
 *  Calcul des N segmentations les plus probables, des profils de segments/ruptures et
 *  des profils d'entropies pour une sequence.
 *
 *  arguments : reference sur un objet StatError, stream,
 *              identificateur de la sequence, nombre de segments, types des modeles,
 *              type de sortie, format ('a' : ASCII, 's' : Spreadsheet),
 *              methode de calcul des segmentations (algorithme de programmation dynamique ou
 *              algorithme forward-backward de simulation), nombre de segmentations.
 *
 *--------------------------------------------------------------*/

bool Sequences::segment_profile_write(StatError &error , ostream &os , int iidentifier ,
                                      int nb_segment , int *model_type , int output ,
                                      char format , int segmentation , int nb_segmentation) const

{
  bool status = true;
  register int i , j;
  int index = I_DEFAULT;
  double likelihood = D_INF , segmentation_likelihood , **rank;
  const FrequencyDistribution **pmarginal;
  FrequencyDistribution *marginal;
  Sequences *seq;


  error.init();

/*  if (((index_parameter_type == TIME) && (index_interval->variance > 0.)) ||
      (index_parameter_type == POSITION)) {
    status = false;
    error.update(SEQ_error[SEQR_INDEX_PARAMETER_TYPE]);
  }
  if (index_parameter_type == POSITION) {
    status = false;
    error.correction_update(SEQ_error[SEQR_INDEX_PARAMETER_TYPE] , SEQ_index_parameter_word[TIME]);
  } */

  for (i = 0;i < nb_variable;i++) {
    if ((model_type[i] == CATEGORICAL_CHANGE) || (model_type[0] == MULTIVARIATE_CATEGORICAL_CHANGE) ||
        (model_type[i] == POISSON_CHANGE) || (model_type[0] == MULTIVARIATE_POISSON_CHANGE) ||
        (model_type[i] == GEOMETRIC_0_CHANGE) || (model_type[i] == GEOMETRIC_1_CHANGE) ||
        (model_type[0] == MULTIVARIATE_GEOMETRIC_0_CHANGE) || (model_type[i] == ORDINAL_GAUSSIAN_CHANGE) ||
        (model_type[i] == BAYESIAN_POISSON_CHANGE)) {
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
        if (((model_type[i] != GEOMETRIC_1_CHANGE) && (min_value[i] < 0)) ||
            ((model_type[i] == GEOMETRIC_1_CHANGE) && (min_value[i] < 1))) {
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
  }

  if (model_type[0] == MULTIVARIATE_CATEGORICAL_CHANGE) {
    pmarginal = new const FrequencyDistribution*[nb_variable];

    for (i = 0;i < nb_variable;i++) {
      pmarginal[i] = marginal_distribution[i];
    }
    marginal = new FrequencyDistribution(nb_variable , pmarginal);

    if ((marginal->nb_value < 2) || (marginal->nb_value > NB_OUTPUT)) {
      status = false;
      error.update(STAT_error[STATR_NB_VALUE]);
    }

    else {
      for (i = 0;i < marginal->nb_value;i++) {
        if (marginal->frequency[i] == 0) {
          status = false;
          ostringstream error_message;
          error_message << STAT_error[STATR_MISSING_VALUE] << " " << i;
          error.update((error_message.str()).c_str());
        }
      }
    }

    delete [] pmarginal;
    delete marginal;
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

  if ((nb_segment < 2) || (nb_segment > (index == I_DEFAULT ? length_distribution->offset : length[index]) / 2)) {
    status = false;
    error.update(SEQ_error[SEQR_NB_SEGMENT]);
  }

  if (nb_segmentation < 2) {
    status = false;
    error.update(SEQ_error[SEQR_NB_SEGMENTATION]);
  }

  if (status) {
    seq = new Sequences(*this , 'a');

    // calcul des rangs pour les variables ordinales

    rank = new double*[seq->nb_variable];

    for (i = 1;i < seq->nb_variable;i++) {
      if (model_type[i - 1] == ORDINAL_GAUSSIAN_CHANGE) {
        rank[i] = seq->marginal_distribution[i]->rank_computation();
      }
      else {
        rank[i] = NULL;
      }
    }

    for (i = 0;i < seq->nb_sequence;i++) {
      if ((index == I_DEFAULT) || (index == i)) {
        if ((model_type[0] != MEAN_CHANGE) && (model_type[0] != INTERCEPT_SLOPE_CHANGE)) {
          likelihood = seq->forward_backward(i , nb_segment , model_type , rank , &os , NULL ,
                                             output , format);
        }
        segmentation_likelihood = seq->forward_backward_dynamic_programming(i , nb_segment , model_type ,
                                                                            rank , &os , NULL , output , format ,
                                                                            likelihood);
        if (segmentation_likelihood == D_INF) {
          status = false;
          error.update(SEQ_error[SEQR_SEGMENTATION_FAILURE]);
        }

        else if ((format == 'a') || (length[index == I_DEFAULT ? 0 : index] <= 256)) {
          switch (segmentation) {
          case FORWARD_DYNAMIC_PROGRAMMING :
            seq->N_segmentation(i , nb_segment , model_type , rank , os , format ,
                                nb_segmentation , likelihood);
            break;
          case FORWARD_BACKWARD_SAMPLING :
            seq->forward_backward_sampling(i , nb_segment , model_type , rank , os ,
                                           format , nb_segmentation);
            break;
          }
        }
      }
    }

    delete seq;

    for (i = 1;i < seq->nb_variable;i++) {
      delete [] rank[i];
    }
    delete [] rank;
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Calcul des N segmentations les plus probables, des profils de segments/ruptures et
 *  des profils d'entropies pour une sequence et ecriture des resultats dans un fichier.
 *
 *  arguments : reference sur un objet StatError, path, identificateur de la sequence,
 *              nombre de segments, types des modeles, type de sortie, format de fichier
 *              ('a' : ASCII, 's' : Spreadsheet), methode de calcul des segmentations
 *              (algorithme de programmation dynamique ou algorithme forward-backward
 *               de simulation), nombre de segmentations.
 *
 *--------------------------------------------------------------*/

bool Sequences::segment_profile_write(StatError &error , const char *path ,
                                      int iidentifier , int nb_segment , int *model_type ,
                                      int output , char format , int segmentation ,
                                      int nb_segmentation) const

{
  bool status = true;
  ofstream out_file(path);


  error.init();

  if (!out_file) {
    status = false;
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {
    status = segment_profile_write(error , out_file , iidentifier , nb_segment , model_type ,
                                   output , format , segmentation , nb_segmentation);
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Calcul des profils de segments/ruptures et des profils d'entropies
 *  pour une sequence et affichage des resultats au format Gnuplot.
 *
 *  arguments : reference sur un objet StatError, prefixe des fichiers,
 *              identificateur de la sequence, nombre de segments, types des modeles,
 *              type de sortie, titre des figures.
 *
 *--------------------------------------------------------------*/

bool Sequences::segment_profile_plot_write(StatError &error , const char *prefix ,
                                           int iidentifier , int nb_segment , int *model_type ,
                                           int output , const char *title) const

{
  bool status = true;
  register int i , j , k;
  int index;
  double likelihood = D_INF , segmentation_likelihood , **rank;
  const FrequencyDistribution **pmarginal;
  FrequencyDistribution *marginal;
  Sequences *seq;
  ostringstream data_file_name[2];
  ofstream *out_data_file;


  error.init();

/*  if (((index_parameter_type == TIME) && (index_interval->variance > 0.)) ||
      (index_parameter_type == POSITION)) {
    status = false;
    error.update(SEQ_error[SEQR_INDEX_PARAMETER_TYPE]);
  }
  if (index_parameter_type == POSITION) {
    status = false;
    error.correction_update(SEQ_error[SEQR_INDEX_PARAMETER_TYPE] , SEQ_index_parameter_word[TIME]);
  } */

  for (i = 0;i < nb_variable;i++) {
    if ((model_type[i] == CATEGORICAL_CHANGE) || (model_type[0] == MULTIVARIATE_CATEGORICAL_CHANGE) ||
        (model_type[i] == POISSON_CHANGE) || (model_type[0] == MULTIVARIATE_POISSON_CHANGE) ||
        (model_type[i] == GEOMETRIC_0_CHANGE) || (model_type[i] == GEOMETRIC_1_CHANGE) ||
        (model_type[0] == MULTIVARIATE_GEOMETRIC_0_CHANGE) || (model_type[i] == ORDINAL_GAUSSIAN_CHANGE) ||
        (model_type[i] == BAYESIAN_POISSON_CHANGE)) {
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
        if (((model_type[i] != GEOMETRIC_1_CHANGE) && (min_value[i] < 0)) ||
            ((model_type[i] == GEOMETRIC_1_CHANGE) && (min_value[i] < 1))) {
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
  }

  if (model_type[0] == MULTIVARIATE_CATEGORICAL_CHANGE) {
    pmarginal = new const FrequencyDistribution*[nb_variable];

    for (i = 0;i < nb_variable;i++) {
      pmarginal[i] = marginal_distribution[i];
    }
    marginal = new FrequencyDistribution(nb_variable , pmarginal);

    if ((marginal->nb_value < 2) || (marginal->nb_value > NB_OUTPUT)) {
      status = false;
      error.update(STAT_error[STATR_NB_VALUE]);
    }

    else {
      for (i = 0;i < marginal->nb_value;i++) {
        if (marginal->frequency[i] == 0) {
          status = false;
          ostringstream error_message;
          error_message << STAT_error[STATR_MISSING_VALUE] << " " << i;
          error.update((error_message.str()).c_str());
        }
      }
    }

    delete [] pmarginal;
    delete marginal;
  }

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

  else if ((nb_segment < 2) || (nb_segment > length[index] / 2)) {
    status = false;
    error.update(SEQ_error[SEQR_NB_SEGMENT]);
  }

  if (status) {

    // ecriture des fichiers de donnees

    i = (((model_type[0] == MEAN_CHANGE) || (model_type[0] == INTERCEPT_SLOPE_CHANGE)) ? 0 : 1);
    data_file_name[i] << prefix << i << ".dat";
    out_data_file = new ofstream((data_file_name[i].str()).c_str());

    if (!out_data_file) {
      status = false;
      error.update(STAT_error[STATR_FILE_PREFIX]);
    }

    else {
      seq = new Sequences(*this , 'a');

      // calcul des rangs pour les variables ordinales

      rank = new double*[seq->nb_variable];

      for (i = 1;i < seq->nb_variable;i++) {
        if (model_type[i - 1] == ORDINAL_GAUSSIAN_CHANGE) {
          rank[i] = seq->marginal_distribution[i]->rank_computation();
        }
        else {
          rank[i] = NULL;
        }
      }

      if ((model_type[0] != MEAN_CHANGE) && (model_type[0] != INTERCEPT_SLOPE_CHANGE)) {
        likelihood = seq->forward_backward(index , nb_segment , model_type , rank ,
                                           out_data_file , NULL , output , 'g');
        out_data_file->close();
        delete out_data_file;

        data_file_name[0] << prefix << 0 << ".dat";
        out_data_file = new ofstream((data_file_name[0].str()).c_str());
      }

#     ifdef DEBUG
      likelihood = D_INF;
#     endif

      segmentation_likelihood = seq->forward_backward_dynamic_programming(index , nb_segment , model_type ,
                                                                          rank , out_data_file , NULL ,
                                                                          output , 'g' , likelihood);
      out_data_file->close();
      delete out_data_file;

      if (segmentation_likelihood == D_INF) {
        status = false;
        error.update(SEQ_error[SEQR_SEGMENTATION_FAILURE]);
      }

      else {

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

          out_file << "set border 15 lw 0\n" << "set tics out\n" << "set xtics nomirror\n";

          if (index_parameter) {
            if (seq->index_parameter[index][seq->length[index] - 1] - seq->index_parameter[index][0] < TIC_THRESHOLD) {
              out_file << "set xtics 0,1" << endl;
            }

            j = 2;
            for (k = 1;k < seq->nb_variable;k++) {
              if ((model_type[k - 1] == POISSON_CHANGE) || (model_type[0] == MULTIVARIATE_POISSON_CHANGE) ||
                  (model_type[k - 1] == GEOMETRIC_0_CHANGE) || (model_type[k - 1] == GEOMETRIC_1_CHANGE) ||
                  (model_type[0] == MULTIVARIATE_GEOMETRIC_0_CHANGE) || (model_type[k - 1] == GAUSSIAN_CHANGE) ||
                  (model_type[k - 1] == VARIANCE_CHANGE) || (model_type[k - 1] == LINEAR_MODEL_CHANGE) ||
                  (model_type[0] == MEAN_CHANGE) || (model_type[0] == MEAN_VARIANCE_CHANGE) ||
                  (model_type[0] == INTERCEPT_SLOPE_CHANGE) || (model_type[k - 1] == BAYESIAN_POISSON_CHANGE) ||
                  (model_type[k - 1] == BAYESIAN_GAUSSIAN_CHANGE)) {
                out_file << "set title \"";
                if (title) {
                  out_file << title;
                  if (seq->nb_variable > 2) {
                    out_file << " - ";
                  }
                }

                if (seq->nb_variable > 2) {
                  out_file << STAT_label[STATL_VARIABLE] << " " << k;
                }
                out_file << "\n\n";

                out_file << "plot [" << seq->index_parameter[index][0] << ":"
                         << seq->index_parameter[index][seq->length[index] - 1] << "] ["
                         << MIN(seq->min_value[k] , 0) << ":"
                         << MAX(seq->max_value[k] , seq->min_value[k] + 1) << "] "
                         << "\"" << label((data_file_name[0].str()).c_str()) << "\" using " << 1 << " : " << j++
                         << " title \"" << SEQ_label[SEQL_SEQUENCE] << " " << iidentifier
                         << "\" with " << (index_interval->variance > index_interval->mean ? "points" : "linespoints") << ",\\" << endl;
                out_file << "\"" << label((data_file_name[0].str()).c_str()) << "\" using " << 1 << " : " << j++
                         << " title \"" << SEQ_label[SEQL_PIECEWISE_LINEAR_FUNCTION] << "\" with lines" << endl;

                if (i == 0) {
                  out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
                }
                out_file << endl;
              }
            }

            out_file << "set title \"";
            if (title) {
              out_file << title << " - ";
            }

            if (likelihood != D_INF) {
              switch (output) {
              case CHANGE_POINT :
                out_file << SEQ_label[SEQL_MAX_POSTERIOR_CHANGE_POINT_PROBABILITY] << "\"\n\n";
                break;
              case SEGMENT :
                out_file << SEQ_label[SEQL_MAX_POSTERIOR_SEGMENT_PROBABILITY] << "\"\n\n";
                break;
              }
            }

            else {
              switch (output) {
              case CHANGE_POINT :
                out_file << SEQ_label[SEQL_MAX_CHANGE_POINT_LIKELIHOOD] << "\"\n\n";
                break;
              case SEGMENT :
                out_file << SEQ_label[SEQL_MAX_SEGMENT_LIKELIHOOD] << "\"\n\n";
                break;
              }
            }

            out_file << "plot [" << seq->index_parameter[index][0] << ":"
                     << seq->index_parameter[index][seq->length[index] - 1];
            if (likelihood != D_INF) {
              out_file << "] [0:"  << exp(segmentation_likelihood - likelihood) << "] ";
            }
            else {
              out_file << "] [0:1] ";
            }
            for (k = 0;k < nb_segment;k++) {
              out_file << "\"" << label((data_file_name[0].str()).c_str()) << "\" using "
                       << 1 << " : " << j++ << " title \""
                       << (output == CHANGE_POINT ? SEQ_label[SEQL_CHANGE_POINT] : SEQ_label[SEQL_SEGMENT])
                       << " " << k << "\" with linespoints";
              if (k < nb_segment - 1) {
                out_file << ",\\";
              }
              out_file << endl;
            }

            if (likelihood != D_INF) {
              if (i == 0) {
                out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
              }
              out_file << endl;

              out_file << "set title \"";
              if (title) {
                out_file << title << " - ";
              }
              switch (output) {
              case CHANGE_POINT :
                out_file << SEQ_label[SEQL_POSTERIOR_CHANGE_POINT_PROBABILITY] << "\"\n\n";
                break;
              case SEGMENT :
                out_file << SEQ_label[SEQL_POSTERIOR_SEGMENT_PROBABILITY] << "\"\n\n";
                break;
              }

              out_file << "plot [" << seq->index_parameter[index][0] << ":"
                       << seq->index_parameter[index][seq->length[index] - 1] << "] [0:1] ";
              for (k = 0;k < nb_segment;k++) {
                out_file << "\"" << label((data_file_name[1].str()).c_str()) << "\" using "
                         << 1 << " : " << k + 2 << " title \""
                         << (output == CHANGE_POINT ? SEQ_label[SEQL_CHANGE_POINT] : SEQ_label[SEQL_SEGMENT])
                         << " " << k << "\" with linespoints";
                if (k < nb_segment - 1) {
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
              out_file << SEQ_label[SEQL_POSTERIOR_CHANGE_POINT_PROBABILITY] << "\"\n\n";

              out_file << "plot [" << seq->index_parameter[index][0] << ":"
                       << seq->index_parameter[index][seq->length[index] - 1] << "] [0:1] ";
              for (k = MAX(1 , nb_segment - 3);k < nb_segment;k++) {
                out_file << "\"" << label((data_file_name[1].str()).c_str()) << "\" using "
                         << 1 << " : " << nb_segment + k + 1 << " title \"" << k + 1 << " "
                         << SEQ_label[SEQL_SEGMENTS] << "\" with linespoints";
                if (k < nb_segment - 1) {
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
              out_file << SEQ_label[SEQL_BEGIN_CONDITIONAL_ENTROPY] << "\"\n\n";

              out_file << "plot [" << seq->index_parameter[index][0] << ":"
                       << seq->index_parameter[index][seq->length[index] - 1]
                       << "] [0:" << log(2.) << "] ";
              for (k = MAX(1 , nb_segment - 3);k < nb_segment;k++) {
                out_file << "\"" << label((data_file_name[1].str()).c_str()) << "\" using "
                         << 1 << " : " << 2 * nb_segment + k << " title \"" << k + 1 << " "
                         << SEQ_label[SEQL_SEGMENTS] << "\" with linespoints";
                if (k < nb_segment - 1) {
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
              out_file << SEQ_label[SEQL_END_CONDITIONAL_ENTROPY] << "\"\n\n";

              out_file << "plot [" << seq->index_parameter[index][0] << ":"
                       << seq->index_parameter[index][seq->length[index] - 1]
                       << "] [0:" << log(2.) << "] ";
              for (k = MAX(1 , nb_segment - 3);k < nb_segment;k++) {
                out_file << "\"" << label((data_file_name[1].str()).c_str()) << "\" using "
                         << 1 << " : " << 3 * nb_segment + k - 1 << " title \"" << k + 1 << " "
                         << SEQ_label[SEQL_SEGMENTS] << "\" with linespoints";
                if (k < nb_segment - 1) {
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
                       << seq->index_parameter[index][seq->length[index] - 1]
                       << "] [0:" << log(2.) << "] "
                       << "\"" << label((data_file_name[1].str()).c_str()) << "\" using "
                       << 1 << " : " << 3 * nb_segment - 1 << " title \""
                       << SEQ_label[SEQL_BEGIN_CONDITIONAL_ENTROPY] << "\" with linespoints,\\" << endl;
              out_file << "\"" << label((data_file_name[1].str()).c_str()) << "\" using "
                       << 1 << " : " << 4 * nb_segment - 2 << " title \""
                       << SEQ_label[SEQL_END_CONDITIONAL_ENTROPY] << "\" with linespoints,\\" << endl;
              out_file << "\"" << label((data_file_name[1].str()).c_str()) << "\" using "
                       << 1 << " : " << 4 * nb_segment - 1 << " title \""
                       << SEQ_label[SEQL_CHANGE_POINT_ENTROPY] << "\" with linespoints" << endl;
            }

            if (seq->index_parameter[index][seq->length[index] - 1] - seq->index_parameter[index][0] < TIC_THRESHOLD) {
              out_file << "set xtics autofreq" << endl;
            }
          }

          else {
            if (seq->length[index] - 1 < TIC_THRESHOLD) {
              out_file << "set xtics 0,1" << endl;
            }

            j = 1;
            for (k = 1;k < seq->nb_variable;k++) {
              if ((model_type[k - 1] == POISSON_CHANGE) || (model_type[0] == MULTIVARIATE_POISSON_CHANGE) ||
                  (model_type[k - 1] == GEOMETRIC_0_CHANGE) || (model_type[k - 1] == GEOMETRIC_1_CHANGE) ||
                  (model_type[0] == MULTIVARIATE_GEOMETRIC_0_CHANGE) || (model_type[k - 1] == GAUSSIAN_CHANGE) ||
                  (model_type[k - 1] == VARIANCE_CHANGE) || (model_type[k - 1] == LINEAR_MODEL_CHANGE) ||
                  (model_type[0] == MEAN_CHANGE) || (model_type[0] == MEAN_VARIANCE_CHANGE) ||
                  (model_type[0] == INTERCEPT_SLOPE_CHANGE) || (model_type[k - 1] == BAYESIAN_POISSON_CHANGE) ||
                  (model_type[k - 1] == BAYESIAN_GAUSSIAN_CHANGE)) {
                out_file << "set title \"";
                if (title) {
                  out_file << title;
                  if (seq->nb_variable > 2) {
                    out_file << " - ";
                  }
                }

                if (seq->nb_variable > 2) {
                  out_file << STAT_label[STATL_VARIABLE] << " " << k;
                }
                out_file << "\n\n";

                out_file << "plot [0:" << seq->length[index] - 1 << "] ["
                         << MIN(seq->min_value[k] , 0) << ":"
                         << MAX(seq->max_value[k] , seq->min_value[k] + 1) << "] "
                         << "\"" << label((data_file_name[0].str()).c_str()) << "\" using " << j++
                         << " title \"" << SEQ_label[SEQL_SEQUENCE] << " " << iidentifier
                         << "\" with linespoints,\\" << endl;
                out_file << "\"" << label((data_file_name[0].str()).c_str()) << "\" using " << j++
                         << " title \"" << SEQ_label[SEQL_PIECEWISE_LINEAR_FUNCTION] << "\" with lines" << endl;

                if (i == 0) {
                  out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
                }
                out_file << endl;
              }
            }

            out_file << "set title \"";
            if (title) {
              out_file << title << " - ";
            }

            if (likelihood != D_INF) {
              switch (output) {
              case CHANGE_POINT :
                out_file << SEQ_label[SEQL_MAX_POSTERIOR_CHANGE_POINT_PROBABILITY] << "\"\n\n";
                break;
              case SEGMENT :
                out_file << SEQ_label[SEQL_MAX_POSTERIOR_SEGMENT_PROBABILITY] << "\"\n\n";
                break;
              }
            }

            else {
              switch (output) {
              case CHANGE_POINT :
                out_file << SEQ_label[SEQL_MAX_CHANGE_POINT_LIKELIHOOD] << "\"\n\n";
                break;
              case SEGMENT :
                out_file << SEQ_label[SEQL_MAX_SEGMENT_LIKELIHOOD] << "\"\n\n";
                break;
              }
            }

            out_file << "plot [0:" << seq->length[index] - 1;
            if (likelihood != D_INF) {
              out_file << "] [0:"  << exp(segmentation_likelihood - likelihood) << "] ";
            }
            else {
              out_file << "] [0:1] ";
            }
            for (k = 0;k < nb_segment;k++) {
              out_file << "\"" << label((data_file_name[0].str()).c_str()) << "\" using "
                       << j++ << " title \""
                       << (output == CHANGE_POINT ? SEQ_label[SEQL_CHANGE_POINT] : SEQ_label[SEQL_SEGMENT])
                       << " " << k << "\" with linespoints";
              if (k < nb_segment - 1) {
                out_file << ",\\";
              }
              out_file << endl;
            }

            if (likelihood != D_INF) {
              if (i == 0) {
                out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
              }
              out_file << endl;

              out_file << "set title \"";
              if (title) {
                out_file << title << " - ";
              }
              switch (output) {
              case CHANGE_POINT :
                out_file << SEQ_label[SEQL_POSTERIOR_CHANGE_POINT_PROBABILITY] << "\"\n\n";
                break;
              case SEGMENT :
                out_file << SEQ_label[SEQL_POSTERIOR_SEGMENT_PROBABILITY] << "\"\n\n";
                break;
              }

              out_file << "plot [0:" << seq->length[index] - 1 << "] [0:1] ";
              for (k = 0;k < nb_segment;k++) {
                out_file << "\"" << label((data_file_name[1].str()).c_str()) << "\" using "
                         << k + 1 << " title \""
                         << (output == CHANGE_POINT ? SEQ_label[SEQL_CHANGE_POINT] : SEQ_label[SEQL_SEGMENT])
                         << " " << k << "\" with linespoints";
                if (k < nb_segment - 1) {
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
              out_file << SEQ_label[SEQL_POSTERIOR_CHANGE_POINT_PROBABILITY] << "\"\n\n";

              out_file << "plot [0:" << seq->length[index] - 1 << "] [0:1] ";
              for (k = MAX(1 , nb_segment - 3);k < nb_segment;k++) {
                out_file << "\"" << label((data_file_name[1].str()).c_str()) << "\" using "
                         << nb_segment + k << " title \"" << k + 1 << " "
                         << SEQ_label[SEQL_SEGMENTS] << "\" with linespoints";
                if (k < nb_segment - 1) {
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
              out_file << SEQ_label[SEQL_BEGIN_CONDITIONAL_ENTROPY] << "\"\n\n";

              out_file << "plot [0:" << seq->length[index] - 1
                       << "] [0:" << log(2.) << "] ";
              for (k = MAX(1 , nb_segment - 3);k < nb_segment;k++) {
                out_file << "\"" << label((data_file_name[1].str()).c_str()) << "\" using "
                         << 2 * nb_segment + k - 1 << " title \"" << k + 1 << " "
                         << SEQ_label[SEQL_SEGMENTS] << "\" with linespoints";
                if (k < nb_segment - 1) {
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
              out_file << SEQ_label[SEQL_END_CONDITIONAL_ENTROPY] << "\"\n\n";

              out_file << "plot [0:" << seq->length[index] - 1
                       << "] [0:" << log(2.) << "] ";
              for (k = MAX(1 , nb_segment - 3);k < nb_segment;k++) {
                out_file << "\"" << label((data_file_name[1].str()).c_str()) << "\" using "
                         << 3 * nb_segment + k - 2 << " title \"" << k + 1 << " "
                         << SEQ_label[SEQL_SEGMENTS] << "\" with linespoints";
                if (k < nb_segment - 1) {
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

              out_file << "plot [0:" << seq->length[index] - 1 << "] [0:" << log(2.) << "] "
                       << "\"" << label((data_file_name[1].str()).c_str()) << "\" using "
                       << 3 * nb_segment - 2 << " title \""
                       << SEQ_label[SEQL_BEGIN_CONDITIONAL_ENTROPY] << "\" with linespoints,\\" << endl;
              out_file << "\"" << label((data_file_name[1].str()).c_str()) << "\" using "
                       << 4 * nb_segment - 3 << " title \""
                       << SEQ_label[SEQL_END_CONDITIONAL_ENTROPY] << "\" with linespoints,\\" << endl;
              out_file << "\"" << label((data_file_name[1].str()).c_str()) << "\" using "
                       << 4 * nb_segment - 2 << " title \""
                       << SEQ_label[SEQL_CHANGE_POINT_ENTROPY] << "\" with linespoints" << endl;
            }

            if (seq->length[index] - 1 < TIC_THRESHOLD) {
              out_file << "set xtics autofreq" << endl;
            }
          }

          if (i == 1) {
            out_file << "\nset terminal x11" << endl;
          }

          out_file << "\npause 0 \"" << STAT_label[STATL_END] << "\"" << endl;
        }
      }

      delete seq;

      for (i = 1;i < seq->nb_variable;i++) {
        delete [] rank[i];
      }
      delete [] rank;
    }
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Calcul des profils de segments/ruptures et des profils d'entropies
 *  pour une sequence et sortie graphique des resultats.
 *
 *  arguments : reference sur un objet StatError, identificateur de la sequence,
 *              nombre de segments, types des modeles, type de sortie.
 *
 *--------------------------------------------------------------*/

MultiPlotSet* Sequences::segment_profile_plotable_write(StatError &error , int iidentifier ,
                                                        int nb_segment , int *model_type ,
                                                        int output) const

{
  bool status = true;
  register int i , j , k;
  int index , nb_plot_set;
  double likelihood = D_INF , segmentation_likelihood , **rank;
  const FrequencyDistribution **pmarginal;
  FrequencyDistribution *marginal;
  Sequences *seq;
  ostringstream title , legend;
  MultiPlotSet *plot_set;


  plot_set = NULL;
  error.init();

/*  if (((index_parameter_type == TIME) && (index_interval->variance > 0.)) ||
      (index_parameter_type == POSITION)) {
    status = false;
    error.update(SEQ_error[SEQR_INDEX_PARAMETER_TYPE]);
  }
  if (index_parameter_type == POSITION) {
    status = false;
    error.correction_update(SEQ_error[SEQR_INDEX_PARAMETER_TYPE] , SEQ_index_parameter_word[TIME]);
  } */

  for (i = 0;i < nb_variable;i++) {
    if ((model_type[i] == CATEGORICAL_CHANGE) || (model_type[0] == MULTIVARIATE_CATEGORICAL_CHANGE) ||
        (model_type[i] == POISSON_CHANGE) || (model_type[0] == MULTIVARIATE_POISSON_CHANGE) ||
        (model_type[i] == GEOMETRIC_0_CHANGE) || (model_type[i] == GEOMETRIC_1_CHANGE) ||
        (model_type[0] == MULTIVARIATE_GEOMETRIC_0_CHANGE) || (model_type[i] == ORDINAL_GAUSSIAN_CHANGE) ||
        (model_type[i] == BAYESIAN_POISSON_CHANGE)) {
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
        if (((model_type[i] != GEOMETRIC_1_CHANGE) && (min_value[i] < 0)) ||
            ((model_type[i] == GEOMETRIC_1_CHANGE) && (min_value[i] < 1))) {
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
  }

  if (model_type[0] == MULTIVARIATE_CATEGORICAL_CHANGE) {
    pmarginal = new const FrequencyDistribution*[nb_variable];

    for (i = 0;i < nb_variable;i++) {
      pmarginal[i] = marginal_distribution[i];
    }
    marginal = new FrequencyDistribution(nb_variable , pmarginal);

    if ((marginal->nb_value < 2) || (marginal->nb_value > NB_OUTPUT)) {
      status = false;
      error.update(STAT_error[STATR_NB_VALUE]);
    }

    else {
      for (i = 0;i < marginal->nb_value;i++) {
        if (marginal->frequency[i] == 0) {
          status = false;
          ostringstream error_message;
          error_message << STAT_error[STATR_MISSING_VALUE] << " " << i;
          error.update((error_message.str()).c_str());
        }
      }
    }

    delete [] pmarginal;
    delete marginal;
  }

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

  else if ((nb_segment < 2) || (nb_segment > length[index] / 2)) {
    status = false;
    error.update(SEQ_error[SEQR_NB_SEGMENT]);
  }

  if (status) {
    seq = new Sequences(*this , 'a');

    // calcul du nombre de vues

    nb_plot_set = 1;
    for (i = 1;i < seq->nb_variable;i++) {
      if ((model_type[i - 1] == POISSON_CHANGE) || (model_type[0] == MULTIVARIATE_POISSON_CHANGE) ||
          (model_type[i - 1] == GEOMETRIC_0_CHANGE) || (model_type[i - 1] == GEOMETRIC_1_CHANGE) ||
          (model_type[0] == MULTIVARIATE_GEOMETRIC_0_CHANGE) || (model_type[i - 1] == GAUSSIAN_CHANGE) ||
          (model_type[i - 1] == VARIANCE_CHANGE) || (model_type[i - 1] == LINEAR_MODEL_CHANGE) ||
          (model_type[0] == MEAN_CHANGE) || (model_type[0] == MEAN_VARIANCE_CHANGE) ||
          (model_type[0] == INTERCEPT_SLOPE_CHANGE) || (model_type[i - 1] == BAYESIAN_POISSON_CHANGE) ||
          (model_type[i - 1] == BAYESIAN_GAUSSIAN_CHANGE)) {
        nb_plot_set++;
      }
    }
    if ((model_type[0] != MEAN_CHANGE) && (model_type[0] != INTERCEPT_SLOPE_CHANGE)) {
      nb_plot_set += 5;
    }

    plot_set = new MultiPlotSet(nb_plot_set);

    MultiPlotSet &plot = *plot_set;

    plot.border = "15 lw 0";

    // calcul des rangs pour les variables ordinales

    rank = new double*[seq->nb_variable];

    for (i = 1;i < seq->nb_variable;i++) {
      if (model_type[i - 1] == ORDINAL_GAUSSIAN_CHANGE) {
        rank[i] = seq->marginal_distribution[i]->rank_computation();
      }
      else {
        rank[i] = NULL;
      }
    }

    if ((model_type[0] != MEAN_CHANGE) && (model_type[0] != INTERCEPT_SLOPE_CHANGE)) {
      likelihood = seq->forward_backward(index , nb_segment , model_type , rank ,
                                         NULL , plot_set , output , 'p');
    }

#   ifdef DEBUG
    likelihood = D_INF;
#   endif

    segmentation_likelihood = seq->forward_backward_dynamic_programming(index , nb_segment , model_type ,
                                                                        rank , NULL , plot_set ,
                                                                        output , 'p' , likelihood);

    if (segmentation_likelihood == D_INF) {
      delete plot_set;
      plot_set = NULL;
      error.update(SEQ_error[SEQR_SEGMENTATION_FAILURE]);
    }

    else {
      i = 0;

      // vues : sequence et fonction lineaire par morceaux

      for (j = 1;j < seq->nb_variable;j++) {
        if ((model_type[j - 1] == POISSON_CHANGE) || (model_type[0] == MULTIVARIATE_POISSON_CHANGE) ||
            (model_type[j - 1] == GEOMETRIC_0_CHANGE) || (model_type[j - 1] == GEOMETRIC_1_CHANGE) ||
            (model_type[0] == MULTIVARIATE_GEOMETRIC_0_CHANGE) || (model_type[j - 1] == GAUSSIAN_CHANGE) ||
            (model_type[j - 1] == VARIANCE_CHANGE) || (model_type[j - 1] == LINEAR_MODEL_CHANGE) ||
            (model_type[0] == MEAN_CHANGE) || (model_type[0] == MEAN_VARIANCE_CHANGE) ||
            (model_type[0] == INTERCEPT_SLOPE_CHANGE) || (model_type[j - 1] == BAYESIAN_POISSON_CHANGE) ||
            (model_type[j - 1] == BAYESIAN_GAUSSIAN_CHANGE)) {
          if (seq->nb_variable > 2) {
            title.str("");
            title << STAT_label[STATL_VARIABLE] << " " << j;
            plot[i].title = title.str();
          }

          if (seq->index_parameter) {
            plot[i].xrange = Range(seq->index_parameter[index][0] , seq->index_parameter[index][seq->length[index] - 1]);
            if (seq->index_parameter[index][seq->length[index] - 1] - seq->index_parameter[index][0] < TIC_THRESHOLD) {
              plot[i].xtics = 1;
            }
          }

          else {
            plot[i].xrange = Range(0 , seq->length[index] - 1);
            if (seq->length[index] - 1 < TIC_THRESHOLD) {
              plot[i].xtics = 1;
            }
          }

          plot[i].yrange = Range(MIN(seq->min_value[j] , 0) , MAX(seq->max_value[j] , seq->min_value[j] + 1));

          legend.str("");
          legend << SEQ_label[SEQL_SEQUENCE] << " " << iidentifier;
          plot[i][0].legend = legend.str();

          plot[i][0].style = "linespoints";

          plot[i][1].legend = SEQ_label[SEQL_PIECEWISE_LINEAR_FUNCTION];
          plot[i][1].style = "lines";
          i++;
        }
      }

      // vue : probabilitees maximum

      if (likelihood != D_INF) {
        switch (output) {
        case CHANGE_POINT :
          plot[i].title = SEQ_label[SEQL_MAX_POSTERIOR_CHANGE_POINT_PROBABILITY];
          break;
        case SEGMENT :
          plot[i].title = SEQ_label[SEQL_MAX_POSTERIOR_SEGMENT_PROBABILITY];
          break;
        }

        plot[i].yrange = Range(0. , exp(segmentation_likelihood - likelihood));
      }

      else {
        switch (output) {
        case CHANGE_POINT :
          plot[i].title = SEQ_label[SEQL_MAX_CHANGE_POINT_LIKELIHOOD];
          break;
        case SEGMENT :
          plot[i].title = SEQ_label[SEQL_MAX_SEGMENT_LIKELIHOOD];
          break;
        }

        plot[i].yrange = Range(0. , 1.);
      }

      if (seq->index_parameter) {
        plot[i].xrange = Range(seq->index_parameter[index][0] , seq->index_parameter[index][seq->length[index] - 1]);
        if (seq->index_parameter[index][seq->length[index] - 1] - seq->index_parameter[index][0] < TIC_THRESHOLD) {
          plot[i].xtics = 1;
        }
      }

      else {
        plot[i].xrange = Range(0 , seq->length[index] - 1);
        if (seq->length[index] - 1 < TIC_THRESHOLD) {
          plot[i].xtics = 1;
        }
      }

      for (j = 0;j < nb_segment;j++) {
        legend.str("");
        legend << (output == CHANGE_POINT ? SEQ_label[SEQL_CHANGE_POINT] : SEQ_label[SEQL_SEGMENT])
               << " " << j;
        plot[i][j].legend = legend.str();

        plot[i][j].style = "linespoints";
      }
      i++;

      // vue : probabilitees lissees

      if (likelihood != D_INF) {
        switch (output) {
        case CHANGE_POINT :
          plot[i].title = SEQ_label[SEQL_POSTERIOR_CHANGE_POINT_PROBABILITY];
          break;
        case SEGMENT :
          plot[i].title = SEQ_label[SEQL_POSTERIOR_SEGMENT_PROBABILITY];
          break;
        }

        if (seq->index_parameter) {
          plot[i].xrange = Range(seq->index_parameter[index][0] , seq->index_parameter[index][seq->length[index] - 1]);
          if (seq->index_parameter[index][seq->length[index] - 1] - seq->index_parameter[index][0] < TIC_THRESHOLD) {
            plot[i].xtics = 1;
          }
        }

        else {
          plot[i].xrange = Range(0 , seq->length[index] - 1);
          if (seq->length[index] - 1 < TIC_THRESHOLD) {
            plot[i].xtics = 1;
          }
        }

        plot[i].yrange = Range(0. , 1.);

        for (j = 0;j < nb_segment;j++) {
          legend.str("");
          legend << (output == CHANGE_POINT ? SEQ_label[SEQL_CHANGE_POINT] : SEQ_label[SEQL_SEGMENT])
                 << " " << j;
          plot[i][j].legend = legend.str();

          plot[i][j].style = "linespoints";
        }
        i++;

        // vue : profils de ruptures

        plot[i].title = SEQ_label[SEQL_POSTERIOR_CHANGE_POINT_PROBABILITY];

        if (seq->index_parameter) {
          plot[i].xrange = Range(seq->index_parameter[index][0] , seq->index_parameter[index][seq->length[index] - 1]);
          if (seq->index_parameter[index][seq->length[index] - 1] - seq->index_parameter[index][0] < TIC_THRESHOLD) {
            plot[i].xtics = 1;
          }
        }

        else {
          plot[i].xrange = Range(0 , seq->length[index] - 1);
          if (seq->length[index] - 1 < TIC_THRESHOLD) {
            plot[i].xtics = 1;
          }
        }

        plot[i].yrange = Range(0. , 1.);

        j = 0;
        for (k = MAX(1 , nb_segment - 3);k < nb_segment;k++) {
          legend.str("");
          legend << k + 1 << " " << SEQ_label[SEQL_SEGMENTS];
          plot[i][j].legend = legend.str();

          plot[i][j].style = "linespoints";
          j++;
        }
        i++;

        // vue : profils d'entropies conditonnees par le passe

        plot[i].title = SEQ_label[SEQL_BEGIN_CONDITIONAL_ENTROPY];

        if (seq->index_parameter) {
          plot[i].xrange = Range(seq->index_parameter[index][0] , seq->index_parameter[index][seq->length[index] - 1]);
          if (seq->index_parameter[index][seq->length[index] - 1] - seq->index_parameter[index][0] < TIC_THRESHOLD) {
            plot[i].xtics = 1;
          }
        }

        else {
          plot[i].xrange = Range(0 , seq->length[index] - 1);
          if (seq->length[index] - 1 < TIC_THRESHOLD) {
            plot[i].xtics = 1;
          }
        }

        plot[i].yrange = Range(0. , log(2.));

        j = 0;
        for (k = MAX(1 , nb_segment - 3);k < nb_segment;k++) {
          legend.str("");
          legend << k + 1 << " " << SEQ_label[SEQL_SEGMENTS];
          plot[i][j].legend = legend.str();

          plot[i][j].style = "linespoints";
          j++;
        }
        i++;

        // vue : profils d'entropies par le futur

        plot[i].title = SEQ_label[SEQL_END_CONDITIONAL_ENTROPY];

        if (seq->index_parameter) {
          plot[i].xrange = Range(seq->index_parameter[index][0] , seq->index_parameter[index][seq->length[index] - 1]);
          if (seq->index_parameter[index][seq->length[index] - 1] - seq->index_parameter[index][0] < TIC_THRESHOLD) {
            plot[i].xtics = 1;
          }
        }

        else {
          plot[i].xrange = Range(0 , seq->length[index] - 1);
          if (seq->length[index] - 1 < TIC_THRESHOLD) {
            plot[i].xtics = 1;
          }
        }

        plot[i].yrange = Range(0. , log(2.));

        j = 0;
        for (k = MAX(1 , nb_segment - 3);k < nb_segment;k++) {
          legend.str("");
          legend << k + 1 << " " << SEQ_label[SEQL_SEGMENTS];
          plot[i][j].legend = legend.str();

          plot[i][j].style = "linespoints";
          j++;
        }
        i++;

        // vue : profils d'entropies

        if (seq->index_parameter) {
          plot[i].xrange = Range(seq->index_parameter[index][0] , seq->index_parameter[index][seq->length[index] - 1]);
          if (seq->index_parameter[index][seq->length[index] - 1] - seq->index_parameter[index][0] < TIC_THRESHOLD) {
            plot[i].xtics = 1;
          }
        }

        else {
          plot[i].xrange = Range(0 , seq->length[index] - 1);
          if (seq->length[index] - 1 < TIC_THRESHOLD) {
            plot[i].xtics = 1;
          }
        }

        plot[i].yrange = Range(0. , log(2.));

        plot[i][0].legend = SEQ_label[SEQL_BEGIN_CONDITIONAL_ENTROPY];
        plot[i][0].style = "linespoints";

        plot[i][1].legend = SEQ_label[SEQL_END_CONDITIONAL_ENTROPY];
        plot[i][1].style = "linespoints";

        plot[i][2].legend = SEQ_label[SEQL_CHANGE_POINT_ENTROPY];
        plot[i][2].style = "linespoints";
      }
    }

    delete seq;

    for (i = 1;i < seq->nb_variable;i++) {
      delete [] rank[i];
    }
    delete [] rank;
  }

  return plot_set;
}


/*--------------------------------------------------------------*
 *
 *  Segmentation hierarchique d'une sequence.
 *
 *  arguments : reference sur un objet StatError, stream,
 *              identificateur de la sequence, nombre de segments maximum,
 *              types des modeles.
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::hierarchical_segmentation(StatError &error , ostream &os , int iidentifier ,
                                                int max_nb_segment , int *model_type) const

{
  bool status = true;
  register int i , j , k;
  int index , max_nb_value , nb_segment , segment_index , split_change_point ,
      begin_change_point , end_change_point , merge , *frequency , *pisequence , *psegment ,
      *nb_parameter , **change_point , ***begin_frequency , ***end_frequency;
  double sum , factorial_sum , diff , buff , max_likelihood , *mean , *prsequence ,
         *likelihood , *penalty , *penalized_likelihood , **rank , **factorial ,
         **begin_sum , **end_sum , **begin_factorial_sum , **end_factorial_sum;
  long double square_sum , merge_contrast , *residual , *begin_contrast , *end_contrast ,
              **begin_square_sum , **end_square_sum;
  Sequences *seq , *iseq;


  seq = NULL;
  error.init();

/*  if (((index_parameter_type == TIME) && (index_interval->variance > 0.)) ||
      (index_parameter_type == POSITION)) {
    status = false;
    error.update(SEQ_error[SEQR_INDEX_PARAMETER_TYPE]);
  } */
  if (index_parameter_type == POSITION) {
    status = false;
    error.correction_update(SEQ_error[SEQR_INDEX_PARAMETER_TYPE] , SEQ_index_parameter_word[TIME]);
  }

  if ((model_type[0] == MEAN_CHANGE) || (model_type[0] == MEAN_VARIANCE_CHANGE) ||
      (model_type[0] == INTERCEPT_SLOPE_CHANGE)) {
    status = false;
    error.update(SEQ_error[SEQR_CHANGE_POINT_MODEL]);
  }

  for (i = 0;i < nb_variable;i++) {
    if ((model_type[i] == CATEGORICAL_CHANGE) || (model_type[i] == ORDINAL_GAUSSIAN_CHANGE) ||
        (model_type[i] == POISSON_CHANGE)) {
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
        if (min_value[i] < 0) {
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
  }

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

  else if ((max_nb_segment <= 2) || (max_nb_segment > length[index] / 2)) {
    status = false;
    error.update(SEQ_error[SEQR_MAX_NB_SEGMENT]);
  }

  if (status) {
    iseq = new Sequences(*this , 1 , &index);
    seq = new Sequences(*iseq , 'a');
    delete iseq;

    // calcul des rangs pour les variables ordinales

    rank = new double*[seq->nb_variable];

    for (i = 1;i < seq->nb_variable;i++) {
      if (model_type[i - 1] == ORDINAL_GAUSSIAN_CHANGE) {
        rank[i] = marginal_distribution[i - 1]->rank_computation();
      }
      else {
        rank[i] = NULL;
      }
    }

    max_nb_value = 0;
    factorial = new double*[seq->nb_variable];

    begin_frequency = new int**[seq->nb_variable];
    end_frequency = new int**[seq->nb_variable];
    begin_sum = new double*[seq->nb_variable];
    end_sum = new double*[seq->nb_variable];
    begin_factorial_sum = new double*[seq->nb_variable];
    end_factorial_sum = new double*[seq->nb_variable];
    begin_square_sum = new long double*[seq->nb_variable];
    end_square_sum = new long double*[seq->nb_variable];

    for (i = 1;i < seq->nb_variable;i++) {
      if ((model_type[i - 1] == CATEGORICAL_CHANGE) &&
          (seq->marginal_distribution[i]->nb_value > max_nb_value)) {
        max_nb_value = seq->marginal_distribution[i]->nb_value;
      }

      if (model_type[i - 1] == CATEGORICAL_CHANGE) {
        begin_frequency[i] = new int*[seq->length[0]];
        end_frequency[i] = new int*[seq->length[0]];
        for (j = 0;j < seq->length[0];j++) {
          begin_frequency[i][j] = new int[seq->marginal_distribution[i]->nb_value];
          end_frequency[i][j] = new int[seq->marginal_distribution[i]->nb_value];
        }

        begin_sum[i] = NULL;
        end_sum[i] = NULL;
      }

      else {
        begin_sum[i] = new double[seq->length[0]];
        end_sum[i] = new double[seq->length[0]];

        begin_frequency[i] = NULL;
        end_frequency[i] = NULL;
      }

      if (model_type[i - 1] == POISSON_CHANGE) {
        begin_factorial_sum[i] = new double[seq->length[0]];
        end_factorial_sum[i] = new double[seq->length[0]];

        factorial[i] = new double[seq->length[0]];
      }

      else {
        begin_factorial_sum[i] = NULL;
        end_factorial_sum[i] = NULL;

        factorial[i] = NULL;
      }

      if ((model_type[i - 1] == ORDINAL_GAUSSIAN_CHANGE) || (model_type[i - 1] == GAUSSIAN_CHANGE) ||
          (model_type[i - 1] == VARIANCE_CHANGE)){
        begin_square_sum[i] = new long double[seq->length[0]];
        end_square_sum[i] = new long double[seq->length[0]];
      }
      else {
        begin_square_sum[i] = NULL;
        end_square_sum[i] = NULL;
      }
    }

    if (max_nb_value > 0) {
      frequency = new int[max_nb_value];
    }
    else {
      frequency = NULL;
    }

    mean = new double[seq->nb_variable];
    residual = new long double[seq->length[0]];

    begin_contrast = new long double[seq->length[0]];
    end_contrast = new long double[seq->length[0]];

    change_point = new int*[max_nb_segment + 1];
    for (i = 1;i <= max_nb_segment;i++) {
      change_point[i] = new int[i + 1];
      change_point[i][0] = 0;
      change_point[i][i] = seq->length[0];
    }

    likelihood = new double[max_nb_segment + 1];
    nb_parameter = new int[max_nb_segment + 1];
    penalty = new double[max_nb_segment + 1];
    penalized_likelihood = new double[max_nb_segment + 1];

    for (i = 1;i < seq->nb_variable;i++) {
      if (model_type[i - 1] == VARIANCE_CHANGE) {
        mean[i] = 0.;

        if (seq->type[i] != REAL_VALUE) {
          pisequence = seq->int_sequence[0][i];
          for (j = 0;j < seq->length[0];j++) {
            mean[i] += *pisequence++;
          }
        }

        else {
          prsequence = seq->real_sequence[0][i];
          for (j = 0;j < seq->length[0];j++) {
            mean[i] += *prsequence++;
          }
        }

        mean[i] /= seq->length[0];
      }
    }

    // calcul des log-vraisemblances des segments

    for (i = 0;i < seq->length[0];i++) {
      begin_contrast[i] = 0.;
    }

    for (i = 1;i < seq->nb_variable;i++) {
      if (model_type[i - 1] == CATEGORICAL_CHANGE) {
        for (j = 0;j < seq->marginal_distribution[i]->nb_value;j++) {
          frequency[j] = 0;
        }
        sum = 0.;

        pisequence = seq->int_sequence[0][i];
        frequency[*pisequence++]++;

        for (j = 0;j < seq->marginal_distribution[i]->nb_value;j++) {
          begin_frequency[i][0][j] = frequency[j];
        }

        for (j = 1;j < seq->length[0];j++) {
          sum += j * log((double)j / (double)(j + 1)) +
                 log((double)(frequency[*pisequence] + 1) / (double)(j + 1));
          if (frequency[*pisequence] > 0) {
            sum -= frequency[*pisequence] *
                   log((double)frequency[*pisequence] / (double)(frequency[*pisequence] + 1));
          }
          frequency[*pisequence++]++;

          for (k = 0;k < seq->marginal_distribution[i]->nb_value;k++) {
            begin_frequency[i][j][k] = frequency[k];
          }

          if (begin_contrast[j] != D_INF) {
            begin_contrast[j] += sum;
          }

/*          frequency[*pisequence++]++;

          for (k = 0;k < seq->marginal_distribution[i]->nb_value;k++) {
            begin_frequency[i][j][k] = frequency[k];
          }

          if (begin_contrast[j] != D_INF) {
            for (k = 0;k < seq->marginal_distribution[i]->nb_value;k++) {
              if (frequency[k] > 0) {
                begin_contrast[j] += frequency[k] * log((double)frequency[k] / (double)(j + 1));
              }
            }
          } */
        }
      }

      else if (model_type[i - 1] == POISSON_CHANGE) {
        sum = 0.;
        factorial_sum = 0.;

        pisequence = seq->int_sequence[0][i];
        for (j = 0;j < seq->length[0];j++) {
          factorial[i][j] = 0.;
          for (k = 2;k <= seq->int_sequence[0][i][j];k++) {
            factorial[i][j] += log((double)k);
          }

          sum += *pisequence++;
          factorial_sum += factorial[i][j];

          begin_sum[i][j] = sum;
          begin_factorial_sum[i][j] = factorial_sum;

          if ((begin_contrast[j] != D_INF) && (sum > 0.)) {
            begin_contrast[j] += sum * (log(sum / (j + 1)) - 1) - factorial_sum;
          }
        }
      }

      else {
        switch (model_type[i - 1]) {

        case VARIANCE_CHANGE : {
          square_sum = 0.;

          if (seq->type[i] != REAL_VALUE) {
            pisequence = seq->int_sequence[0][i];
            for (j = 0;j < seq->length[0];j++) {
              diff = *pisequence++ - mean[i];
              square_sum += diff * diff;

              begin_square_sum[i][j] = square_sum;

              residual[j] = square_sum;
            }
          }

          else {
            prsequence = seq->real_sequence[0][i];
            for (j = 0;j < seq->length[0];j++) {
              diff = *prsequence++ - mean[i];
              square_sum += diff * diff;

              begin_square_sum[i][j] = square_sum;

              residual[j] = square_sum;
            }
          }
          break;
        }

        case ORDINAL_GAUSSIAN_CHANGE : {
          pisequence = seq->int_sequence[0][i];
          square_sum = 0.;
          sum = rank[i][*pisequence++];

          begin_square_sum[i][0] = square_sum;
          begin_sum[i][0] = sum;

          residual[0] = 0.;

          for (j = 1;j < seq->length[0];j++) {
            diff = rank[i][*pisequence] - sum / j;
            square_sum += ((double)j / (double)(j + 1)) * diff * diff;
            sum += rank[i][*pisequence++];

            begin_square_sum[i][j] = square_sum;
            begin_sum[i][j] = sum;

            residual[j] = square_sum;
          }
          break;
        }

        case GAUSSIAN_CHANGE : {
          if (seq->type[i] != REAL_VALUE) {
            pisequence = seq->int_sequence[0][i];
            square_sum = 0.;
            sum = *pisequence++;

            begin_square_sum[i][0] = square_sum;
            begin_sum[i][0] = sum;

            residual[0] = 0.;

            for (j = 1;j < seq->length[0];j++) {
              diff = *pisequence - sum / j;
              square_sum += ((double)j / (double)(j + 1)) * diff * diff;
              sum += *pisequence++;

              begin_square_sum[i][j] = square_sum;
              begin_sum[i][j] = sum;

              residual[j] = square_sum;
            }
          }

          else {
            prsequence = seq->real_sequence[0][i];
            square_sum = 0.;
            sum = *prsequence++;

            begin_square_sum[i][0] = square_sum;
            begin_sum[i][0] = sum;

            residual[0] = 0.;

            for (j = 1;j < seq->length[0];j++) {
              diff = *prsequence - sum / j;
              square_sum += ((double)j / (double)(j + 1)) * diff * diff;
              sum += *prsequence++;

              begin_square_sum[i][j] = square_sum;
              begin_sum[i][j] = sum;

              residual[j] = square_sum;
            }
          }
          break;
        }
        }

        for (j = 0;j < seq->length[0];j++) {
//          if ((begin_contrast[j] != D_INF) && (residual[j] > 0.)) {
          if ((begin_contrast[j] != D_INF) && (residual[j] > sqrt((double)(j + 1)) * ROUNDOFF_ERROR)) {
            begin_contrast[j] -= ((double)(j + 1) / 2.) * (logl(residual[j] /
                                   (j + 1)) + log(2 * M_PI) + 1);
/*            begin_contrast[j] -= ((double)(j + 1) / 2.) * (logl(residual[j] /
                                   j) + log(2 * M_PI)) + (double)j / 2.; */
          }
          else {
            begin_contrast[j] = D_INF;
          }
        }
      }
    }

    for (i = seq->length[0] - 1;i >= 0;i--) {
      end_contrast[i] = 0.;
    }

    for (i = 1;i < seq->nb_variable;i++) {
      if (model_type[i - 1] == CATEGORICAL_CHANGE) {
/*        for (j = 0;j < seq->marginal_distribution[i]->nb_value;j++) {
          frequency[j] = 0;
        }
        sum = 0.;

        pisequence = seq->int_sequence[0][i] + seq->length[0] - 1;
        frequency[*pisequence--]++;

        for (j = 0;j < seq->marginal_distribution[i]->nb_value;j++) {
          end_frequency[i][seq->length[0] - 1][j] = frequency[j];
        }

        for (j = seq->length[0] - 2;j >= 0;j--) {
          sum += (seq->length[0] - j - 1) *
                 log((double)(seq->length[0] - j - 1) / (double)(seq->length[0] - j)) +
                 log((double)(frequency[*pisequence] + 1) / (double)(seq->length[0] - j));
          if (frequency[*pisequence] > 0) {
            sum -= frequency[*pisequence] *
                   log((double)frequency[*pisequence] / (double)(frequency[*pisequence] + 1));
          }
          frequency[*pisequence--]++;

          for (k = 0;k < seq->marginal_distribution[i]->nb_value;k++) {
            end_frequency[i][j][k] = frequency[k];
          }

          if (end_contrast[j] != D_INF) {
            end_contrast[j] += sum;
          }

          frequency[*pisequence--]++;

          for (k = 0;k < seq->marginal_distribution[i]->nb_value;k++) {
            end_frequency[i][j][k] = frequency[k];
          }

          if (end_contrast[j] != D_INF) {
            for (k = 0;k < seq->marginal_distribution[i]->nb_value;k++) {
              if (frequency[k] > 0) {
                end_contrast[j] += frequency[k] * log((double)frequency[k] / (double)(seq->length[0] - j));
              }
            }
          }
        } */

        for (j = seq->length[0] - 1;j > 0;j--) {
          for (k = 0;k < seq->marginal_distribution[i]->nb_value;k++) {
            end_frequency[i][j][k] = begin_frequency[i][seq->length[0] - 1][k] - begin_frequency[i][j - 1][k];
          }

          if (end_contrast[j] != D_INF) {
            for (k = 0;k < seq->marginal_distribution[i]->nb_value;k++) {
              if (end_frequency[i][j][k] > 0) {
                end_contrast[j] += end_frequency[i][j][k] * log((double)end_frequency[i][j][k] /
                                   (double)(seq->length[0] - j));
              }
            }
          }
        }

        for (k = 0;k < seq->marginal_distribution[i]->nb_value;k++) {
          end_frequency[i][0][k] = begin_frequency[i][seq->length[0] - 1][k];
        }

        if (end_contrast[0] != D_INF) {
          for (k = 0;k < seq->marginal_distribution[i]->nb_value;k++) {
            if (end_frequency[i][0][k] > 0) {
              end_contrast[0] += end_frequency[i][0][k] * log((double)end_frequency[i][0][k] /
                                 (double)seq->length[0]);
            }
          }
        }
      }

      else if (model_type[i - 1] == POISSON_CHANGE) {
/*        sum = 0.;
        factorial_sum = 0.;

        pisequence = seq->int_sequence[0][i] + seq->length[0] - 1;
        for (j = seq->length[0] - 1;j >= 0;j--) {
          sum += *pisequence--;
          factorial_sum += factorial[i][j];

          end_sum[i][j] = sum;
          end_factorial_sum[i][j] = factorial_sum;

          if ((end_contrast[j] != D_INF) && (sum > 0.)) {
            end_contrast[j] += sum * (log(sum / (seq->length[0] - j)) - 1) - factorial_sum;
          }
        } */

        for (j = seq->length[0] - 1;j > 0;j--) {
          end_sum[i][j] = begin_sum[i][seq->length[0] - 1] - begin_sum[i][j - 1];
          end_factorial_sum[i][j] = begin_factorial_sum[i][seq->length[0] - 1] - begin_factorial_sum[i][j - 1];

          if ((end_contrast[j] != D_INF) && (end_sum[i][j] > 0.)) {
            end_contrast[j] += end_sum[i][j] * (log(end_sum[i][j] / (seq->length[0] - j)) - 1) -
                               end_factorial_sum[i][j];
          }
        }

        end_sum[i][0] = begin_sum[i][seq->length[0] - 1];
        end_factorial_sum[i][0] = begin_factorial_sum[i][seq->length[0] - 1];

        if ((end_contrast[0] != D_INF) && (end_sum[i][0] > 0.)) {
          end_contrast[0] += end_sum[i][0] * (log(end_sum[i][0] / seq->length[0]) - 1) -
                             end_factorial_sum[i][0];
        }
      }

      else {
/*        switch (model_type[i - 1]) {

        case VARIANCE_CHANGE : {
          square_sum = 0.;

          if (seq->type[i] != REAL_VALUE) {
            pisequence = seq->int_sequence[0][i] + seq->length[0] - 1;
            for (j = seq->length[0] - 1;j >= 0;j--) {
              diff = *pisequence-- - mean[i];
              square_sum += diff * diff;

              end_square_sum[i][j] = square_sum;

              residual[j] = square_sum;
            }
          }

          else {
            prsequence = seq->real_sequence[0][i] + seq->length[0] - 1;
            for (j = seq->length[0] - 1;j >= 0;j--) {
              diff = *prsequence-- - mean[i];
              square_sum += diff * diff;

              end_square_sum[i][j] = square_sum;

              residual[j] = square_sum;
            }
          }
          break;
        }

        case ORDINAL_GAUSSIAN_CHANGE : {
          pisequence = seq->int_sequence[0][i] + seq->length[0] - 1;
          square_sum = 0.;
          sum = rank[i][*pisequence--];

          end_square_sum[i][seq->length[0] - 1] = square_sum;
          end_sum[i][seq->length[0] - 1] = sum;

          residual[seq->length[0] - 1] = 0.;

          for (j = seq->length[0] - 2;j >= 0;j--) {
            diff = rank[i][*pisequence] - sum / (seq->length[0] - j - 1);
            square_sum += ((double)(seq->length[0] - j - 1) / (double)(seq->length[0] - j)) * diff * diff;
            sum += rank[i][*pisequence--];

            end_square_sum[i][j] = square_sum;
            end_sum[i][j] = sum;

            residual[j] = square_sum;
          }
          break;
        }

        case GAUSSIAN_CHANGE : {
          if (seq->type[i] != REAL_VALUE) {
            pisequence = seq->int_sequence[0][i] + seq->length[0] - 1;
            square_sum = 0.
            sum = *pisequence--;

            end_square_sum[i][seq->length[0] - 1] = square_sum;
            end_sum[i][seq->length[0] - 1] = sum;

            residual[seq->length[0] - 1] = 0.;

            for (j = seq->length[0] - 2;j >= 0;j--) {
              diff = *pisequence - sum / (seq->length[0] - j - 1);
              square_sum += ((double)(seq->length[0] - j - 1) / (double)(seq->length[0] - j)) * diff * diff;
              sum += *pisequence--;

              end_square_sum[i][j] = square_sum;
              end_sum[i][j] = sum;

              residual[j] = square_sum;
            }
          }

          else {
            prsequence = seq->real_sequence[0][i] + seq->length[0] - 1;
            square_sum = 0.;
            sum = *prsequence--;

            end_square_sum[i][seq->length[0] - 1] = square_sum;
            end_sum[i][seq->length[0] - 1] = sum;

            residual[seq->length[0] - 1] = 0.;

            for (j = seq->length[0] - 2;j >= 0;j--) {
              diff = *prsequence - sum / (seq->length[0] - j - 1);
              square_sum += ((double)(seq->length[0] - j - 1) / (double)(seq->length[0] - j)) * diff * diff;
              sum += *prsequence--;

              end_square_sum[i][j] = square_sum;
              end_sum[i][j] = sum;

              residual[j] = square_sum;
            }
          }
          break;
        }
        } */

        if (model_type[i - 1] == VARIANCE_CHANGE) {
          for (j = seq->length[0] - 1;j > 0;j--) {
            end_square_sum[i][j] = begin_square_sum[i][seq->length[0] - 1] - begin_square_sum[i][j - 1];
            residual[j] = end_square_sum[i][j];
          }

          end_square_sum[i][0] = begin_square_sum[i][seq->length[0] - 1];
          residual[0] = end_square_sum[i][0];
        }

        else {
          for (j = seq->length[0] - 1;j > 0;j--) {
            end_sum[i][j] = begin_sum[i][seq->length[0] - 1] - begin_sum[i][j - 1];
            diff = begin_sum[i][j - 1] / j  - end_sum[i][j] / (seq->length[0] - j);
            end_square_sum[i][j] = begin_square_sum[i][seq->length[0] - 1] - begin_square_sum[i][j - 1] -
                                   ((double)(j * (seq->length[0] - j)) / (double)seq->length[0]) * diff * diff;

            residual[j] = end_square_sum[i][j];
          }

          end_square_sum[i][0] = begin_square_sum[i][seq->length[0] - 1];
          end_sum[i][0] = begin_sum[i][seq->length[0] - 1];

          residual[0] = end_square_sum[i][0];
        }

        for (j = seq->length[0] - 1;j >= 0;j--) {
//          if ((end_contrast[j] != D_INF) && (residual[j] > 0.)) {
          if ((end_contrast[j] != D_INF) && (residual[j] > sqrt((double)(seq->length[0] - j)) * ROUNDOFF_ERROR)) {
            end_contrast[j] -= ((double)(seq->length[0] - j) / 2.) * (logl(residual[j] /
                                 (seq->length[0] - j)) + log(2 * M_PI) + 1);
/*            end_contrast[j] -= ((double)(seq->length[0] - j) / 2.) * (logl(residual[j] /
                                 (seq->length[0] - j - 1)) + log(2 * M_PI)) +
                               (double)(seq->length[0] - j - 1) / 2.; */
          }
          else {
            end_contrast[j] = D_INF;
          }
        }
      }
    }

    nb_segment = 1;
    likelihood[nb_segment] = begin_contrast[seq->length[0] - 1];
//    likelihood[nb_segment] = end_contrast[0];

    // calcul de la  vraisemblance penalisee au sens du BIC modifie (Zhang & Siegmund, 2007)

    psegment = seq->int_sequence[0][0];
    for (i = 0;i < seq->length[0];i++) {
      *psegment++ = 0;
    }

    nb_parameter[nb_segment] = seq->nb_parameter_computation(0 , nb_segment , model_type);

    penalty[nb_segment] = log((double)seq->length[0]);

    penalized_likelihood[nb_segment] = 2 * likelihood[nb_segment] - nb_parameter[nb_segment] *
                                       log((double)((seq->nb_variable - 1) * seq->length[0])) - penalty[nb_segment];

    // segmentation optimale en 2 segments

    likelihood[nb_segment + 1] = likelihood[nb_segment];
    for (i = 0;i < seq->length[0] - 1;i++) {
      if ((begin_contrast[i] != D_INF) && (end_contrast[i + 1] != D_INF)) {
        buff = begin_contrast[i] + end_contrast[i + 1];
        if (buff > likelihood[nb_segment + 1]) {
          likelihood[nb_segment + 1] = buff;
          split_change_point = i + 1;
        }
      }
    }

    if (likelihood[nb_segment + 1] > likelihood[nb_segment]) {
      segment_index = 1;
      begin_change_point = split_change_point;
      end_change_point = split_change_point;

      nb_segment++;
      change_point[nb_segment][segment_index] = split_change_point;

      // calcul de la  vraisemblance penalisee au sens du BIC modifie

      psegment = seq->int_sequence[0][0];
      penalty[nb_segment] = 0.;
      for (i = 0;i < nb_segment;i++) {
        for (j = change_point[nb_segment][i];j < change_point[nb_segment][i + 1];j++) {
          *psegment++ = i;
        }

        penalty[nb_segment] += log((double)(change_point[nb_segment][i + 1] -
                                            change_point[nb_segment][i]));
      }

      nb_parameter[nb_segment] = seq->nb_parameter_computation(0 , nb_segment , model_type);

      penalized_likelihood[nb_segment] = 2 * likelihood[nb_segment] - nb_parameter[nb_segment] *
                                         log((double)((seq->nb_variable - 1) * seq->length[0])) - penalty[nb_segment];

#     ifdef MESSAGE
      os << "\n" << nb_segment << " " << SEQ_label[SEQL_SEGMENTS] << ":";
      os << " " << likelihood[nb_segment] << " " << penalty[nb_segment]
         << " " << penalized_likelihood[nb_segment] << " ||";

      if (seq->index_parameter) {
        for (i = 0;i < nb_segment;i++) {
          os << " " << seq->index_parameter[0][change_point[nb_segment][i]];
        }
        os << endl;
      }

      else {
        for (i = 0;i <= nb_segment;i++) {
          os << " " << change_point[nb_segment][i];
        }
        os << endl;
      }
#     endif

    }

//    while ((penalized_likelihood[nb_segment] > penalized_likelihood[nb_segment - 1]) &&
    while (nb_segment < max_nb_segment) {

      // calcul des log-vraisemblances des segments modifies

      for (i = begin_change_point;i < change_point[nb_segment][segment_index + 1];i++) {
        begin_contrast[i] = 0.;
      }

      for (i = 1;i < seq->nb_variable;i++) {
        if (model_type[i - 1] == CATEGORICAL_CHANGE) {
          if (begin_change_point < split_change_point) {
            for (j = 0;j < seq->marginal_distribution[i]->nb_value;j++) {
              frequency[j] = begin_frequency[i][begin_change_point - 1][j];
            }
            sum = 0.;

            pisequence = seq->int_sequence[0][i] + begin_change_point;
            for (j = begin_change_point;j < split_change_point;j++) {
              sum += (j - change_point[nb_segment][segment_index - 1]) *
                     log((double)(j - change_point[nb_segment][segment_index - 1]) /
                         (double)(j - change_point[nb_segment][segment_index - 1] + 1)) +
                     log((double)(frequency[*pisequence] + 1) /
                         (double)(j - change_point[nb_segment][segment_index - 1] + 1));
              if (frequency[*pisequence] > 0) {
                sum -= frequency[*pisequence] *
                       log((double)frequency[*pisequence] / (double)(frequency[*pisequence] + 1));
              }
              frequency[*pisequence++]++;

              for (k = 0;k < seq->marginal_distribution[i]->nb_value;k++) {
                begin_frequency[i][j][k] = frequency[k];
              }

              if (begin_contrast[j] != D_INF) {
                begin_contrast[j] += sum;
              }

/*              frequency[*pisequence++]++;

              for (k = 0;k < seq->marginal_distribution[i]->nb_value;k++) {
                begin_frequency[i][j][k] = frequency[k];
              }

              if (begin_contrast[j] != D_INF) {
                for (k = 0;k < seq->marginal_distribution[i]->nb_value;k++) {
                  if (frequency[k] > 0) {
                    begin_contrast[j] += frequency[k] * log((double)frequency[k] /
                                                            (double)(j - change_point[nb_segment][segment_index - 1] + 1));
                  }
                }
              } */
            }
          }

          else {
            pisequence = seq->int_sequence[0][i] + split_change_point;
          }

          for (j = 0;j < seq->marginal_distribution[i]->nb_value;j++) {
            frequency[j] = 0;
          }
          sum = 0.;

          frequency[*pisequence++]++;

          for (j = 0;j < seq->marginal_distribution[i]->nb_value;j++) {
            begin_frequency[i][split_change_point][j] = frequency[j];
          }

          for (j = split_change_point + 1;j < change_point[nb_segment][segment_index + 1];j++) {
            sum += (j - split_change_point) *
                   log((double)(j - split_change_point) / (double)(j - split_change_point + 1)) +
                   log((double)(frequency[*pisequence] + 1) / (double)(j - split_change_point + 1));
            if (frequency[*pisequence] > 0) {
              sum -= frequency[*pisequence] *
                     log((double)frequency[*pisequence] / (double)(frequency[*pisequence] + 1));
            }
            frequency[*pisequence++]++;

            for (k = 0;k < seq->marginal_distribution[i]->nb_value;k++) {
              begin_frequency[i][j][k] = frequency[k];
            }

            if (begin_contrast[j] != D_INF) {
              begin_contrast[j] += sum;
            }

/*            frequency[*pisequence++]++;

            for (k = 0;k < seq->marginal_distribution[i]->nb_value;k++) {
              begin_frequency[i][j][k] = frequency[k];
            }

            if (begin_contrast[j] != D_INF) {
              for (k = 0;k < seq->marginal_distribution[i]->nb_value;k++) {
                if (frequency[k] > 0) {
                  begin_contrast[j] += frequency[k] * log((double)frequency[k] /
                                                          (double)(j - split_change_point + 1));
                }
              }
            } */
          }
        }

        else if (model_type[i - 1] == POISSON_CHANGE) {
          if (begin_change_point < split_change_point) {
            sum = begin_sum[i][begin_change_point - 1];
            factorial_sum = begin_factorial_sum[i][begin_change_point - 1];

            pisequence = seq->int_sequence[0][i] + begin_change_point;
            for (j = begin_change_point;j < split_change_point;j++) {
              sum += *pisequence++;
              factorial_sum += factorial[i][j];

              begin_sum[i][j] = sum;
              begin_factorial_sum[i][j] = factorial_sum;

              if ((begin_contrast[j] != D_INF) && (sum > 0.)) {
                begin_contrast[j] += sum * (log(sum / (j - change_point[nb_segment][segment_index - 1] + 1)) - 1) -
                                     factorial_sum;
              }
            }
          }

          else {
            pisequence = seq->int_sequence[0][i] + split_change_point;
          }

          sum = 0.;
          factorial_sum = 0.;

          for (j = split_change_point;j < change_point[nb_segment][segment_index + 1];j++) {
            sum += *pisequence++;
            factorial_sum += factorial[i][j];

            begin_sum[i][j] = sum;
            begin_factorial_sum[i][j] = factorial_sum;

            if ((begin_contrast[j] != D_INF) && (sum > 0.)) {
              begin_contrast[j] += sum * (log(sum / (j - split_change_point + 1)) - 1) -
                                   factorial_sum;
            }
          }
        }

        else {
          switch (model_type[i - 1]) {

          case VARIANCE_CHANGE : {
            if (begin_change_point < split_change_point) {
              square_sum = begin_square_sum[i][begin_change_point - 1];

              if (seq->type[i] != REAL_VALUE) {
                pisequence = seq->int_sequence[0][i] + begin_change_point;
                for (j = begin_change_point;j < split_change_point;j++) {
                  diff = *pisequence++ - mean[i];
                  square_sum += diff * diff;

                  begin_square_sum[i][j] = square_sum;

                  residual[j] = square_sum;
                }
              }

              else {
                prsequence = seq->real_sequence[0][i] + begin_change_point;
                for (j = begin_change_point;j < split_change_point;j++) {
                  diff = *prsequence++ - mean[i];
                  square_sum += diff * diff;

                  begin_square_sum[i][j] = square_sum;

                  residual[j] = square_sum;
                }
              }
            }

            else {
              if (seq->type[i] != REAL_VALUE) {
                pisequence = seq->int_sequence[0][i] + split_change_point;
              }
              else {
                prsequence = seq->real_sequence[0][i] + split_change_point;
              }
            }

            square_sum = 0.;

            if (seq->type[i] != REAL_VALUE) {
              for (j = split_change_point;j < change_point[nb_segment][segment_index + 1];j++) {
                diff = *pisequence++ - mean[i];
                square_sum += diff * diff;

                begin_square_sum[i][j] = square_sum;

                residual[j] = square_sum;
              }
            }

            else {
              for (j = split_change_point;j < change_point[nb_segment][segment_index + 1];j++) {
                diff = *prsequence++ - mean[i];
                square_sum += diff * diff;

                begin_square_sum[i][j] = square_sum;

                residual[j] = square_sum;
              }
            }
            break;
          }

          case ORDINAL_GAUSSIAN_CHANGE : {
            if (begin_change_point < split_change_point) {
              square_sum = begin_square_sum[i][begin_change_point - 1];
              sum = begin_sum[i][begin_change_point - 1];

              pisequence = seq->int_sequence[0][i] + begin_change_point;
              for (j = begin_change_point;j < split_change_point;j++) {
                diff = rank[i][*pisequence] - sum / (j - change_point[nb_segment][segment_index - 1]);
                square_sum += ((double)(j - change_point[nb_segment][segment_index - 1]) /
                               (double)(j - change_point[nb_segment][segment_index - 1] + 1)) * diff * diff;
                sum += rank[i][*pisequence++];
 
                begin_square_sum[i][j] = square_sum;
                begin_sum[i][j] = sum;

                residual[j] = square_sum;
              }
            }

            else {
              pisequence = seq->int_sequence[0][i] + split_change_point;
            }

            square_sum = 0.;
            sum = rank[i][*pisequence++];

            begin_square_sum[i][split_change_point] = square_sum;
            begin_sum[i][split_change_point] = sum;

            residual[split_change_point] = 0.;

            for (j = split_change_point + 1;j < change_point[nb_segment][segment_index + 1];j++) {
              diff = rank[i][*pisequence] - sum / (j - split_change_point);
              square_sum += ((double)(j - split_change_point) /
                             (double)(j - split_change_point + 1)) * diff * diff;
              sum += rank[i][*pisequence++];
 
              begin_square_sum[i][j] = square_sum;
              begin_sum[i][j] = sum;

              residual[j] = square_sum;
            }
            break;
          }

          case GAUSSIAN_CHANGE : {
            if (begin_change_point < split_change_point) {
              square_sum = begin_square_sum[i][begin_change_point - 1];
              sum = begin_sum[i][begin_change_point - 1];

              pisequence = seq->int_sequence[0][i] + begin_change_point;
              if (seq->type[i] != REAL_VALUE) {
                for (j = begin_change_point;j < split_change_point;j++) {
                  diff = *pisequence - sum / (j - change_point[nb_segment][segment_index - 1]);
                  square_sum += ((double)(j - change_point[nb_segment][segment_index - 1]) /
                                 (double)(j - change_point[nb_segment][segment_index - 1] + 1)) * diff * diff;
                  sum += *pisequence++;

                  begin_square_sum[i][j] = square_sum;
                  begin_sum[i][j] = sum;

                  residual[j] = square_sum;
                }
              }

              else {
                prsequence = seq->real_sequence[0][i] + begin_change_point;
                for (j = begin_change_point;j < split_change_point;j++) {
                  diff = *prsequence - sum / (j - change_point[nb_segment][segment_index - 1]);
                  square_sum += ((double)(j - change_point[nb_segment][segment_index - 1]) /
                                 (double)(j - change_point[nb_segment][segment_index - 1] + 1)) * diff * diff;
                  sum += *prsequence++;

                  begin_square_sum[i][j] = square_sum;
                  begin_sum[i][j] = sum;

                  residual[j] = square_sum;
                }
              }
            }

            else {
              if (seq->type[i] != REAL_VALUE) {
                pisequence = seq->int_sequence[0][i] + split_change_point;
              }
              else {
                prsequence = seq->real_sequence[0][i] + split_change_point;
              }
            }

            if (seq->type[i] != REAL_VALUE) {
              square_sum = 0.;
              sum = *pisequence++;

              begin_square_sum[i][split_change_point] = square_sum;
              begin_sum[i][split_change_point] = sum;

              residual[split_change_point] = 0.;

              for (j = split_change_point + 1;j < change_point[nb_segment][segment_index + 1];j++) {
                diff = *pisequence - sum / (j - split_change_point);
                square_sum += ((double)(j - split_change_point) /
                               (double)(j - split_change_point + 1)) * diff * diff;
                sum += *pisequence++;

                begin_square_sum[i][j] = square_sum;
                begin_sum[i][j] = sum;

                residual[j] = square_sum;
              }
            }

            else {
              square_sum = 0.;
              sum = *prsequence++;

              begin_square_sum[i][split_change_point] = square_sum;
              begin_sum[i][split_change_point] = sum;

              residual[split_change_point] = 0.;

              for (j = split_change_point + 1;j < change_point[nb_segment][segment_index + 1];j++) {
                diff = *prsequence - sum / (j - split_change_point);
                square_sum += ((double)(j - split_change_point) /
                               (double)(j - split_change_point + 1)) * diff * diff;
                sum += *prsequence++;

                begin_square_sum[i][j] = square_sum;
                begin_sum[i][j] = sum;

                residual[j] = square_sum;
              }
            }
            break;
          }
          }

          if (begin_change_point < split_change_point) {
            for (j = begin_change_point;j < split_change_point;j++) {
//              if ((begin_contrast[j] != D_INF) && (residual[j] > 0.)) {
              if ((begin_contrast[j] != D_INF) && (residual[j] > sqrt((double)(j - change_point[nb_segment][segment_index - 1] + 1)) * ROUNDOFF_ERROR)) {
                begin_contrast[j] -= ((double)(j - change_point[nb_segment][segment_index - 1] + 1) / 2.) * (logl(residual[j] /
                                       (j - change_point[nb_segment][segment_index - 1] + 1)) + log(2 * M_PI) + 1);
/*                begin_contrast[j] -= ((double)(j - change_point[nb_segment][segment_index - 1] + 1) / 2.) * (logl(residual[j] /
                                       (j - change_point[nb_segment][segment_index - 1])) + log(2 * M_PI)) +
                                     (double)(j - change_point[nb_segment][segment_index - 1]) / 2.; */
              }
              else {
                begin_contrast[j] = D_INF;
              }
            }
          }

          for (j = split_change_point;j < change_point[nb_segment][segment_index + 1];j++) {
//            if ((begin_contrast[j] != D_INF) && (residual[j] > 0.)) {
            if ((begin_contrast[j] != D_INF) && (residual[j] > sqrt((double)(j - split_change_point + 1)) * ROUNDOFF_ERROR)) {
              begin_contrast[j] -= ((double)(j - split_change_point + 1) / 2.) * (logl(residual[j] /
                                     (j - split_change_point + 1)) + log(2 * M_PI) + 1);
/*              begin_contrast[j] -= ((double)(j - split_change_point + 1) / 2.) * (logl(residual[j] /
                                     (j - split_change_point)) + log(2 * M_PI)) +
                                   (double)(j - split_change_point) / 2.; */
            }
            else {
              begin_contrast[j] = D_INF;
            }
          }
        }
      }

      for (i = end_change_point - 1;i >= change_point[nb_segment][segment_index - 1];i--) {
        end_contrast[i] = 0.;
      }

      for (i = 1;i < seq->nb_variable;i++) {
        if (model_type[i - 1] == CATEGORICAL_CHANGE) {
/*          if (end_change_point > split_change_point) {
            for (j = 0;j < seq->marginal_distribution[i]->nb_value;j++) {
              frequency[j] = end_frequency[i][end_change_point][j];
            }
            sum = 0.;

            pisequence = seq->int_sequence[0][i] + end_change_point - 1;
            for (j = end_change_point - 1;j >= split_change_point;j--) {
              sum += (change_point[nb_segment][segment_index + 1] - j - 1) *
                     log((double)(change_point[nb_segment][segment_index + 1] - j - 1) /
                         (double)(change_point[nb_segment][segment_index + 1] - j)) +
                     log((double)(frequency[*pisequence] + 1) /
                         (double)(change_point[nb_segment][segment_index + 1] - j));
              if (frequency[*pisequence] > 0) {
                sum -= frequency[*pisequence] *
                       log((double)frequency[*pisequence] / (double)(frequency[*pisequence] + 1));
              }
              frequency[*pisequence--]++;

              for (k = 0;k < seq->marginal_distribution[i]->nb_value;k++) {
                end_frequency[i][j][k] = frequency[k];
              }

              if (end_contrast[j] != D_INF) {
                end_contrast[j] += sum;
              }

              frequency[*pisequence--]++;

              for (k = 0;k < seq->marginal_distribution[i]->nb_value;k++) {
                end_frequency[i][j][k] = frequency[k];
              }

              if (end_contrast[j] != D_INF) {
                for (k = 0;k < seq->marginal_distribution[i]->nb_value;k++) {
                  if (frequency[k] > 0) {
                    end_contrast[j] += frequency[k] * log((double)frequency[k] /
                                                          (double)(change_point[nb_segment][segment_index + 1] - j));
                  }
                }
              }
            }
          }

          else {
            pisequence = seq->int_sequence[0][i] + split_change_point - 1;
          }

          for (j = 0;j < seq->marginal_distribution[i]->nb_value;j++) {
            frequency[j] = 0;
          }
          sum = 0.;

          frequency[*pisequence--]++;

          for (j = 0;j < seq->marginal_distribution[i]->nb_value;j++) {
            end_frequency[i][split_change_point - 1][j] = frequency[j];
          }

          for (j = split_change_point - 2;j >= change_point[nb_segment][segment_index - 1];j--) {
            sum += (split_change_point - j - 1) *
                   log((double)(split_change_point - j - 1) / (double)(split_change_point - j)) +
                   log((double)(frequency[*pisequence] + 1) / (double)(split_change_point - j));
            if (frequency[*pisequence] > 0) {
              sum -= frequency[*pisequence] *
                     log((double)frequency[*pisequence] / (double)(frequency[*pisequence] + 1));
            }
            frequency[*pisequence--]++;

            for (k = 0;k < seq->marginal_distribution[i]->nb_value;k++) {
              end_frequency[i][j][k] = frequency[k];
            }

            if (end_contrast[j] != D_INF) {
              end_contrast[j] += sum;
            }

            frequency[*pisequence--]++;

            for (k = 0;k < seq->marginal_distribution[i]->nb_value;k++) {
              end_frequency[i][j][k] = frequency[k];
            }

            if (end_contrast[j] != D_INF) {
              for (k = 0;k < seq->marginal_distribution[i]->nb_value;k++) {
                if (frequency[k] > 0) {
                  end_contrast[j] += frequency[k] * log((double)frequency[k] /
                                                        (double)(split_change_point - j));
                }
              }
            }
          } */

          if (end_change_point > split_change_point) {
            for (j = end_change_point - 1;j > split_change_point;j--) {
              for (k = 0;k < seq->marginal_distribution[i]->nb_value;k++) {
                end_frequency[i][j][k] = begin_frequency[i][change_point[nb_segment][segment_index + 1] - 1][k] -
                                         begin_frequency[i][j - 1][k];
              }

              if (end_contrast[j] != D_INF) {
                for (k = 0;k < seq->marginal_distribution[i]->nb_value;k++) {
                  if (end_frequency[i][j][k] > 0) {
                    end_contrast[j] += end_frequency[i][j][k] * log((double)end_frequency[i][j][k] /
                                                                    (double)(change_point[nb_segment][segment_index + 1] - j));
                  }
                }
              }
            }

            for (k = 0;k < seq->marginal_distribution[i]->nb_value;k++) {
              end_frequency[i][j][k] = begin_frequency[i][change_point[nb_segment][segment_index + 1] - 1][k];
            }

            if (end_contrast[j] != D_INF) {
              for (k = 0;k < seq->marginal_distribution[i]->nb_value;k++) {
                if (end_frequency[i][j][k] > 0) {
                  end_contrast[j] += end_frequency[i][j][k] * log((double)end_frequency[i][j][k] /
                                                                  (double)(change_point[nb_segment][segment_index + 1] - j));
                }
              }
            }
          }

          for (j = split_change_point - 1;j > change_point[nb_segment][segment_index - 1];j--) {
            for (k = 0;k < seq->marginal_distribution[i]->nb_value;k++) {
              end_frequency[i][j][k] = begin_frequency[i][split_change_point - 1][k] - begin_frequency[i][j - 1][k];
            }

            if (end_contrast[j] != D_INF) {
              for (k = 0;k < seq->marginal_distribution[i]->nb_value;k++) {
                if (end_frequency[i][j][k] > 0) {
                  end_contrast[j] += end_frequency[i][j][k] * log((double)end_frequency[i][j][k] /
                                                                  (double)(split_change_point - j));
                }
              }
            }
          }

          for (k = 0;k < seq->marginal_distribution[i]->nb_value;k++) {
            end_frequency[i][j][k] = begin_frequency[i][split_change_point - 1][k];
          }

          if (end_contrast[j] != D_INF) {
            for (k = 0;k < seq->marginal_distribution[i]->nb_value;k++) {
              if (end_frequency[i][j][k] > 0) {
                end_contrast[j] += end_frequency[i][j][k] * log((double)end_frequency[i][j][k] /
                                                                (double)(split_change_point - j));
              }
            }
          }
        }

        else if (model_type[i - 1] == POISSON_CHANGE) {
/*          if (end_change_point > split_change_point) {
            sum = end_sum[i][end_change_point];
            factorial_sum = end_factorial_sum[i][end_change_point];

            pisequence = seq->int_sequence[0][i] + end_change_point - 1;
            for (j = end_change_point - 1;j >= split_change_point;j--) {
              sum += *pisequence--;
              factorial_sum += factorial[i][j];

              end_sum[i][j] = sum;
              end_factorial_sum[i][j] = factorial_sum;

              if ((end_contrast[j] != D_INF) && (sum > 0.)) {
                end_contrast[j] += sum * (log(sum / (change_point[nb_segment][segment_index + 1] - j)) - 1) -
                                   factorial_sum;
              }
            }
          }

          else {
            pisequence = seq->int_sequence[0][i] + split_change_point - 1;
          }

          sum = 0.;
          factorial_sum = 0.;

          for (j = split_change_point - 1;j >= change_point[nb_segment][segment_index - 1];j--) {
            sum += *pisequence--;
            factorial_sum += factorial[i][j];

            end_sum[i][j] = sum;
            end_factorial_sum[i][j] = factorial_sum;

            if ((end_contrast[j] != D_INF) && (sum > 0.)) {
              end_contrast[j] += sum * (log(sum / (split_change_point - j)) - 1) -
                                 factorial_sum;
            }
          } */

          if (end_change_point > split_change_point) {
            for (j = end_change_point - 1;j > split_change_point;j--) {
              end_sum[i][j] = begin_sum[i][change_point[nb_segment][segment_index + 1] - 1] -
                              begin_sum[i][j - 1];
              end_factorial_sum[i][j] = begin_factorial_sum[i][change_point[nb_segment][segment_index + 1] - 1] -
                                        begin_factorial_sum[i][j - 1];

              if ((end_contrast[j] != D_INF) && (end_sum[i][j] > 0.)) {
                end_contrast[j] += end_sum[i][j] * (log(end_sum[i][j] / (change_point[nb_segment][segment_index + 1] - j)) - 1) -
                                   end_factorial_sum[i][j];
              }
            }

            end_sum[i][j] = begin_sum[i][change_point[nb_segment][segment_index + 1] - 1];
            end_factorial_sum[i][j] = begin_factorial_sum[i][change_point[nb_segment][segment_index + 1] - 1];

            if ((end_contrast[j] != D_INF) && (end_sum[i][j] > 0.)) {
              end_contrast[j] += end_sum[i][j] * (log(end_sum[i][j] / (change_point[nb_segment][segment_index + 1] - j)) - 1) -
                                 end_factorial_sum[i][j];
            }
          }

          for (j = split_change_point - 1;j > change_point[nb_segment][segment_index - 1];j--) {
            end_sum[i][j] = begin_sum[i][split_change_point - 1] - begin_sum[i][j - 1];
            end_factorial_sum[i][j] = begin_factorial_sum[i][split_change_point - 1] -
                                      begin_factorial_sum[i][j - 1];

            if ((end_contrast[j] != D_INF) && (end_sum[i][j] > 0.)) {
              end_contrast[j] += end_sum[i][j] * (log(end_sum[i][j] / (split_change_point - j)) - 1) -
                                 end_factorial_sum[i][j];
            }
          }

          end_sum[i][j] = begin_sum[i][split_change_point - 1];
          end_factorial_sum[i][j] = begin_factorial_sum[i][split_change_point - 1];

          if ((end_contrast[j] != D_INF) && (end_sum[i][j] > 0.)) {
            end_contrast[j] += end_sum[i][j] * (log(end_sum[i][j] / (split_change_point - j)) - 1) -
                               end_factorial_sum[i][j];
          }
        }

        else {
/*          switch (model_type[i - 1]) {

          case VARIANCE_CHANGE : {
            if (end_change_point > split_change_point) {
              square_sum = end_square_sum[i][end_change_point];

              if (seq->type[i] != REAL_VALUE) {
                pisequence = seq->int_sequence[0][i] + end_change_point - 1;
                for (j = end_change_point - 1;j >= split_change_point;j--) {
                  diff = *pisequence-- - mean[i];
                  square_sum += diff * diff;

                  end_square_sum[i][j] = square_sum;

                  residual[j] = square_sum;
                }
              }

              else {
                prsequence = seq->real_sequence[0][i] + end_change_point - 1;
                for (j = end_change_point - 1;j >= split_change_point;j--) {
                  diff = *prsequence-- - mean[i];
                  square_sum += diff * diff;

                  end_square_sum[i][j] = square_sum;

                  residual[j] = square_sum;
                }
              }
            }

            else {
              if (seq->type[i] != REAL_VALUE) {
                pisequence = seq->int_sequence[0][i] + split_change_point - 1;
              }
              else {
                prsequence = seq->real_sequence[0][i] + split_change_point - 1;
              }
            }

            square_sum = 0.;

            if (seq->type[i] != REAL_VALUE) {
              for (j = split_change_point - 1;j >= change_point[nb_segment][segment_index - 1];j--) {
                diff = *pisequence-- - mean[i];
                square_sum += diff * diff;

                end_square_sum[i][j] = square_sum;

                residual[j] = square_sum;
              }
            }

            else {
              for (j = split_change_point - 1;j >= change_point[nb_segment][segment_index - 1];j--) {
                diff = *prsequence-- - mean[i];
                square_sum += diff * diff;

                end_square_sum[i][j] = square_sum;

                residual[j] = square_sum;
              }
            }
            break;
          }

          case ORDINAL_GAUSSIAN_CHANGE : {
            if (end_change_point > split_change_point) {
              square_sum = end_square_sum[i][end_change_point];
              sum = end_sum[i][end_change_point];

              pisequence = seq->int_sequence[0][i] + end_change_point - 1;
              for (j = end_change_point - 1;j >= split_change_point;j--) {
                diff = rank[i][*pisequence] - sum / (change_point[nb_segment][segment_index + 1] - j - 1);
                square_sum += ((double)(change_point[nb_segment][segment_index + 1] - j - 1) /
                               (double)(change_point[nb_segment][segment_index + 1] - j)) * diff * diff;
                sum += rank[i][*pisequence--];

                end_square_sum[i][j] = square_sum;
                end_sum[i][j] = sum;

                residual[j] = square_sum;
              }
            }

            else {
              pisequence = seq->int_sequence[0][i] + split_change_point - 1;
            }

            square_sum = 0.;
            sum = rank[i][*pisequence--];

            end_square_sum[i][split_change_point - 1] = square_sum;
            end_sum[i][split_change_point - 1] = sum;

            residual[split_change_point - 1] = 0.;

            for (j = split_change_point - 2;j >= change_point[nb_segment][segment_index - 1];j--) {
              diff = rank[i][*pisequence] - sum / (split_change_point - j - 1);
              square_sum += ((double)(split_change_point - j - 1) /
                             (double)(split_change_point - j)) * diff * diff;
              sum += rank[i][*pisequence--];

              end_square_sum[i][j] = square_sum;
              end_sum[i][j] = sum;

              residual[j] = square_sum;
            }
            break;
          }

          case GAUSSIAN_CHANGE : {
            if (end_change_point > split_change_point) {
              square_sum = end_square_sum[i][end_change_point];
              sum = end_sum[i][end_change_point];

              if (seq->type[i] != REAL_VALUE) {
                pisequence = seq->int_sequence[0][i] + end_change_point - 1;
                for (j = end_change_point - 1;j >= split_change_point;j--) {
                  diff = *pisequence - sum / (change_point[nb_segment][segment_index + 1] - j - 1);
                  square_sum += ((double)(change_point[nb_segment][segment_index + 1] - j - 1) /
                                 (double)(change_point[nb_segment][segment_index + 1] - j)) * diff * diff;
                  sum += *pisequence--;

                  end_square_sum[i][j] = square_sum;
                  end_sum[i][j] = sum;

                  residual[j] = square_sum;
                }
              }

              else {
                prsequence = seq->real_sequence[0][i] + end_change_point - 1;
                for (j = end_change_point - 1;j >= split_change_point;j--) {
                  diff = *prsequence - sum / (change_point[nb_segment][segment_index + 1] - j - 1);
                  square_sum += ((double)(change_point[nb_segment][segment_index + 1] - j - 1) /
                                 (double)(change_point[nb_segment][segment_index + 1] - j)) * diff * diff;
                  sum += *prsequence--;

                  end_square_sum[i][j] = square_sum;
                  end_sum[i][j] = sum;

                  residual[j] = square_sum;
                }
              }
            }

            else {
              if (seq->type[i] != REAL_VALUE) {
                pisequence = seq->int_sequence[0][i] + split_change_point - 1;
              }
              else {
                prsequence = seq->real_sequence[0][i] + split_change_point - 1;
              }
            }

            if (seq->type[i] != REAL_VALUE) {
              square_sum = 0.;
              sum = *pisequence--;

              end_square_sum[i][split_change_point - 1] = square_sum;
              end_sum[i][split_change_point - 1] = sum;

              residual[split_change_point - 1] = 0.;

              for (j = split_change_point - 2;j >= change_point[nb_segment][segment_index - 1];j--) {
                diff = *pisequence - sum / (split_change_point - j - 1);
                square_sum += ((double)(split_change_point - j - 1) /
                               (double)(split_change_point - j)) * diff * diff;
                sum += *pisequence--;

                end_square_sum[i][j] = square_sum;
                end_sum[i][j] = sum;

                residual[j] = square_sum;
              }
            }

            else {
              square_sum = 0.;
              sum = *prsequence--;

              end_square_sum[i][split_change_point - 1] = square_sum;
              end_sum[i][split_change_point - 1] = sum;

              residual[split_change_point - 1] = 0.;

              for (j = split_change_point - 2;j >= change_point[nb_segment][segment_index - 1];j--) {
                diff = *prsequence - sum / (split_change_point - j - 1);
                square_sum += ((double)(split_change_point - j - 1) /
                               (double)(split_change_point - j)) * diff * diff;
                sum += *prsequence--;

                end_square_sum[i][j] = square_sum;
                end_sum[i][j] = sum;

                residual[j] = square_sum;
              }
            }
            break;
          }
          } */

          if (model_type[i - 1] == VARIANCE_CHANGE) {
            if (end_change_point > split_change_point) {
              for (j = end_change_point - 1;j > split_change_point;j--) {
                end_square_sum[i][j] = begin_square_sum[i][change_point[nb_segment][segment_index + 1] - 1] -
                                       begin_square_sum[i][j - 1];
                residual[j] = end_square_sum[i][j];
              }

              end_square_sum[i][j] = begin_square_sum[i][change_point[nb_segment][segment_index + 1] - 1];
              residual[j] = end_square_sum[i][j];
            }

            for (j = split_change_point - 1;j > change_point[nb_segment][segment_index - 1];j--) {
              end_square_sum[i][j] = begin_square_sum[i][split_change_point - 1] -
                                     begin_square_sum[i][j - 1];
              residual[j] = end_square_sum[i][j];
            }

            end_square_sum[i][j] = begin_square_sum[i][split_change_point - 1];
            residual[j] = end_square_sum[i][j];
          }

          else {
            if (end_change_point > split_change_point) {
              for (j = end_change_point - 1;j > split_change_point;j--) {
                end_sum[i][j] = begin_sum[i][change_point[nb_segment][segment_index + 1] - 1] -
                                begin_sum[i][j - 1];
                diff = begin_sum[i][j - 1] / (j - split_change_point) -
                       end_sum[i][j] / (change_point[nb_segment][segment_index + 1] - j);
                end_square_sum[i][j] = begin_square_sum[i][change_point[nb_segment][segment_index + 1] - 1] -
                                       begin_square_sum[i][j - 1] - ((double)((j - split_change_point) *
                                         (change_point[nb_segment][segment_index + 1] - j)) /
                                        (double)(change_point[nb_segment][segment_index + 1] - split_change_point)) * diff * diff;

                residual[j] = end_square_sum[i][j];
              }

              end_square_sum[i][j] = begin_square_sum[i][change_point[nb_segment][segment_index + 1] - 1];
              end_sum[i][j] = begin_sum[i][change_point[nb_segment][segment_index + 1] - 1];

              residual[j] = end_square_sum[i][j];
            }

            for (j = split_change_point - 1;j > change_point[nb_segment][segment_index - 1];j--) {
              end_sum[i][j] = begin_sum[i][split_change_point - 1] - begin_sum[i][j - 1];
              diff = begin_sum[i][j - 1] / (j - change_point[nb_segment][segment_index - 1]) -
                     end_sum[i][j] / (split_change_point - j);
              end_square_sum[i][j] = begin_square_sum[i][split_change_point - 1] - begin_square_sum[i][j - 1] -
                                     ((double)((j - change_point[nb_segment][segment_index - 1]) *
                                       (split_change_point - j)) /
                                      (double)(split_change_point - change_point[nb_segment][segment_index - 1])) * diff * diff;

              residual[j] = end_square_sum[i][j];
            }

            end_square_sum[i][j] = begin_square_sum[i][split_change_point - 1];
            end_sum[i][j] = begin_sum[i][split_change_point - 1];

            residual[j] = end_square_sum[i][j];
          }

          if (end_change_point > split_change_point) {
            for (j = end_change_point - 1;j >= split_change_point;j--) {
//              if ((end_contrast[j] != D_INF) && (residual[j] > 0.)) {
              if ((end_contrast[j] != D_INF) && (residual[j] > sqrt((double)(change_point[nb_segment][segment_index + 1] - j)) * ROUNDOFF_ERROR)) {
                end_contrast[j] -= ((double)(change_point[nb_segment][segment_index + 1] - j) / 2.) * (logl(residual[j] /
                                     (change_point[nb_segment][segment_index + 1] - j)) + log(2 * M_PI) + 1);
/*                end_contrast[j] -= ((double)(change_point[nb_segment][segment_index + 1] - j) / 2.) * (logl(residual[j] /
                                     (change_point[nb_segment][segment_index + 1] - j - 1)) + log(2 * M_PI)) +
                                   (double)(change_point[nb_segment][segment_index + 1] - j - 1) / 2.; */
              }
              else {
                end_contrast[j] = D_INF;
              }
            }
          }

          for (j = split_change_point - 1;j >= change_point[nb_segment][segment_index - 1];j--) {
//            if ((end_contrast[j] != D_INF) && (residual[j] > 0.)) {
            if ((end_contrast[j] != D_INF) && (residual[j] > sqrt((double)(split_change_point - j)) * ROUNDOFF_ERROR)) {
              end_contrast[j] -= ((double)(split_change_point - j) / 2.) * (logl(residual[j] /
                                   (split_change_point - j)) + log(2 * M_PI) + 1);
/*              end_contrast[j] -= ((double)(split_change_point - j) / 2.) * (logl(residual[j] /
                                   (split_change_point - j - 1)) + log(2 * M_PI)) +
                                 (double)(split_change_point - j - 1) / 2.; */
            }
            else {
              end_contrast[j] = D_INF;
            }
          }
        }
      }

#     ifdef MESSAGE
      if (begin_change_point < split_change_point) {
        if (index_parameter) {
          cout << "\nBegin merge: " << seq->index_parameter[0][change_point[nb_segment][segment_index - 1]]
               << " " << seq->index_parameter[0][begin_change_point] << " " << seq->index_parameter[0][split_change_point]
               << " | " << begin_contrast[begin_change_point - 1] << " " << begin_contrast[begin_change_point] << "\n" << endl;
        }

        else {
          cout << "\nBegin merge: " << change_point[nb_segment][segment_index - 1] << " "
               << begin_change_point << " " << split_change_point
               << " | " << begin_contrast[begin_change_point - 1] << " " << begin_contrast[begin_change_point] << "\n" << endl;
        }
      }

      if (end_change_point > split_change_point) {
        if ((index_parameter) && (segment_index < nb_segment - 1)) {
          cout << "\nEnd merge: " << seq->index_parameter[0][split_change_point] << " " << seq->index_parameter[0][end_change_point]
               << " " << seq->index_parameter[0][change_point[nb_segment][segment_index + 1]]
               << " | " << end_contrast[end_change_point] << " " << end_contrast[end_change_point - 1] << "\n" << endl;
        }

        else {
          cout << "\nEnd merge: " << split_change_point << " " << end_change_point
               << " " << change_point[nb_segment][segment_index + 1]
               << " | " << end_contrast[end_change_point] << " " << end_contrast[end_change_point - 1] << "\n" << endl;
        }
      }
#     endif

      // division d'un segment

      likelihood[nb_segment + 1] = likelihood[nb_segment];
      for (i = 0;i < nb_segment;i++) {
        for (j = change_point[nb_segment][i];j < change_point[nb_segment][i + 1] - 1;j++) {
          if ((begin_contrast[j] != D_INF) && (end_contrast[j + 1] != D_INF)) {
//            buff = likelihood[nb_segment] - begin_contrast[change_point[nb_segment][i + 1] - 1] +
            buff = likelihood[nb_segment] - end_contrast[change_point[nb_segment][i]] +
                   begin_contrast[j] + end_contrast[j + 1];
            if (buff > likelihood[nb_segment + 1]) {
              likelihood[nb_segment + 1] = buff;
              split_change_point = j + 1;
              segment_index = i + 1;
            }
          }
        }
      }

      // fusion de deux segments

      if (likelihood[nb_segment + 1] > likelihood[nb_segment]) {
        merge = 0;

        if (segment_index > 1) {
          merge_contrast = 0.;

          for (i = 1;i < seq->nb_variable;i++) {
            if (model_type[i - 1] == CATEGORICAL_CHANGE) {
              for (j = 0;j < seq->marginal_distribution[i]->nb_value;j++) {
                frequency[j] = begin_frequency[i][change_point[nb_segment][segment_index - 1] - 1][j] +
                               begin_frequency[i][split_change_point - 1][j];
              }

              if (merge_contrast != D_INF) {
                for (j = 0;j < seq->marginal_distribution[i]->nb_value;j++) {
                  if (frequency[j] > 0) {
                    merge_contrast += frequency[j] * log((double)frequency[j] /
                                                         (double)(split_change_point - change_point[nb_segment][segment_index - 2]));
                  }
                }
              }
            }

            else if (model_type[i - 1] == POISSON_CHANGE) {
              sum = begin_sum[i][change_point[nb_segment][segment_index - 1] - 1] +
                    begin_sum[i][split_change_point - 1];
              factorial_sum = begin_factorial_sum[i][change_point[nb_segment][segment_index - 1] - 1] +
                              begin_factorial_sum[i][split_change_point - 1];

              if ((merge_contrast != D_INF) && (sum > 0.)) {
                merge_contrast += sum * (log(sum / (split_change_point - change_point[nb_segment][segment_index - 2])) - 1) -
                                  factorial_sum;
              }
            }

            else {
              if (model_type[i - 1] == VARIANCE_CHANGE) {
                buff = begin_square_sum[i][change_point[nb_segment][segment_index - 1] - 1] +
                       begin_square_sum[i][split_change_point - 1];
              }

              else {
                diff = begin_sum[i][change_point[nb_segment][segment_index - 1] - 1] /
                       (change_point[nb_segment][segment_index - 1] - change_point[nb_segment][segment_index - 2]) -
                       begin_sum[i][split_change_point - 1] /
                       (split_change_point - change_point[nb_segment][segment_index - 1]);
                buff = begin_square_sum[i][change_point[nb_segment][segment_index - 1] - 1] +
                       begin_square_sum[i][split_change_point - 1] +
                       ((double)((change_point[nb_segment][segment_index - 1] - change_point[nb_segment][segment_index - 2]) *
                         (split_change_point - change_point[nb_segment][segment_index - 1])) /
                        (double)(split_change_point - change_point[nb_segment][segment_index - 2])) * diff * diff;
              }

              if ((merge_contrast != D_INF) && (buff > 0.)) {
                merge_contrast -= ((double)(split_change_point - change_point[nb_segment][segment_index - 2]) / 2.) * (log(buff /
                                    (split_change_point - change_point[nb_segment][segment_index - 2])) + log(2 * M_PI) + 1);
/*                merge_contrast -= ((double)(split_change_point - change_point[nb_segment][segment_index - 2]) / 2.) * (log(buff /
                                    (split_change_point - change_point[nb_segment][segment_index - 2] - 1)) + log(2 * M_PI)) +
                                  (double)(split_change_point - change_point[nb_segment][segment_index - 2] - 1) / 2.; */
              }
              else {
                merge_contrast = D_INF;
              }
            }
          }

          buff = likelihood[nb_segment + 1] - begin_contrast[change_point[nb_segment][segment_index - 1] - 1] -
                 begin_contrast[split_change_point - 1] + merge_contrast;
          if (buff > likelihood[nb_segment]) {
            likelihood[nb_segment] = buff;
            merge = -1;
          }
        }

        if (segment_index < nb_segment) {
          merge_contrast = 0.;

          for (i = 1;i < seq->nb_variable;i++) {
            if (model_type[i - 1] == CATEGORICAL_CHANGE) {
              for (j = 0;j < seq->marginal_distribution[i]->nb_value;j++) {
                frequency[j] = end_frequency[i][split_change_point][j] +
                               end_frequency[i][change_point[nb_segment][segment_index]][j];
              }

              if (merge_contrast != D_INF) {
                for (j = 0;j < seq->marginal_distribution[i]->nb_value;j++) {
                  if (frequency[j] > 0) {
                    merge_contrast += frequency[j] * log((double)frequency[j] /
                                                         (double)(change_point[nb_segment][segment_index + 1] - split_change_point));
                  }
                }
              }
            }

            else if (model_type[i - 1] == POISSON_CHANGE) {
              sum = end_sum[i][split_change_point] +
                    end_sum[i][change_point[nb_segment][segment_index]];
              factorial_sum = end_factorial_sum[i][split_change_point] +
                              end_factorial_sum[i][change_point[nb_segment][segment_index]];

              if ((merge_contrast != D_INF) && (sum > 0.)) {
                merge_contrast += sum * (log(sum / (change_point[nb_segment][segment_index + 1] - split_change_point)) - 1) -
                                  factorial_sum;
              }
            }

            else {
              if (model_type[i - 1] == VARIANCE_CHANGE) {
                buff = end_square_sum[i][split_change_point] +
                       end_square_sum[i][change_point[nb_segment][segment_index]];
              }

              else {
                diff = end_sum[i][split_change_point] /
                      (change_point[nb_segment][segment_index] - split_change_point) -
                      end_sum[i][change_point[nb_segment][segment_index]] /
                      (change_point[nb_segment][segment_index + 1] - change_point[nb_segment][segment_index]);
                buff = end_square_sum[i][split_change_point] +
                       end_square_sum[i][change_point[nb_segment][segment_index]] +
                       ((double)((change_point[nb_segment][segment_index] - split_change_point) *
                         (change_point[nb_segment][segment_index + 1] - change_point[nb_segment][segment_index])) /
                        (double)(change_point[nb_segment][segment_index + 1] - split_change_point)) * diff * diff;
              }

              if ((merge_contrast != D_INF) && (buff > 0.)) {
                merge_contrast -= ((double)(change_point[nb_segment][segment_index + 1] - split_change_point) / 2.) * (log(buff /
                                    (change_point[nb_segment][segment_index + 1] - split_change_point)) + log(2 * M_PI) + 1);
/*                merge_contrast -= ((double)(change_point[nb_segment][segment_index + 1] - split_change_point) / 2.) * (log(buff /
                                    (change_point[nb_segment][segment_index + 1] - split_change_point - 1)) + log(2 * M_PI)) +
                                  (double)(change_point[nb_segment][segment_index + 1] - split_change_point - 1) / 2.; */
              }
              else {
                merge_contrast = D_INF;
              }
            }
          }

          buff = likelihood[nb_segment + 1] - end_contrast[split_change_point] -
                 end_contrast[change_point[nb_segment][segment_index]] + merge_contrast;
          if (buff > likelihood[nb_segment]) {
            likelihood[nb_segment] = buff;
            begin_change_point = split_change_point;
            end_change_point = change_point[nb_segment][segment_index];
            change_point[nb_segment][segment_index] = split_change_point;
            merge = 1;
          }
        }

        if (merge == -1) {
          begin_change_point = change_point[nb_segment][segment_index - 1];
          end_change_point = split_change_point;
          change_point[nb_segment][segment_index - 1] = split_change_point;
          segment_index--;
        }

        else if (merge == 0) {
          begin_change_point = split_change_point;
          end_change_point = split_change_point;

          nb_segment++;
          for (i = nb_segment - 1;i > segment_index;i--) {
            change_point[nb_segment][i] = change_point[nb_segment - 1][i - 1];
          }
          change_point[nb_segment][segment_index] = split_change_point;
          for (i = segment_index - 1;i > 0;i--) {
            change_point[nb_segment][i] = change_point[nb_segment - 1][i];
          }
        }

        // calcul de la  vraisemblance penalisee au sens du BIC modifie

        psegment = seq->int_sequence[0][0];
        penalty[nb_segment] = 0.;
        for (i = 0;i < nb_segment;i++) {
          for (j = change_point[nb_segment][i];j < change_point[nb_segment][i + 1];j++) {
            *psegment++ = i;
          }

          penalty[nb_segment] += log((double)(change_point[nb_segment][i + 1] -
                                              change_point[nb_segment][i]));
        }

        nb_parameter[nb_segment] = seq->nb_parameter_computation(0 , nb_segment , model_type);

        penalized_likelihood[nb_segment] = 2 * likelihood[nb_segment] - nb_parameter[nb_segment] *
                                           log((double)((seq->nb_variable - 1) * seq->length[0])) - penalty[nb_segment];

#       ifdef MESSAGE
        os << nb_segment << " " << SEQ_label[SEQL_SEGMENTS] << ":"
           << " " << likelihood[nb_segment] << " " << penalty[nb_segment]
           << " " << penalized_likelihood[nb_segment] << " ||";
        os << " " << segment_index << " " << split_change_point
           << " | " << begin_change_point << " " << change_point[nb_segment][segment_index + 1]
           << " | " << end_change_point << " " << change_point[nb_segment][segment_index - 1] << " ||";

        if (seq->index_parameter) {
          for (i = 0;i < nb_segment;i++) {
            os << " " << seq->index_parameter[0][change_point[nb_segment][i]];
          }
          os << endl;
        }

        else {
          for (i = 0;i <= nb_segment;i++) {
            os << " " << change_point[nb_segment][i];
          }
          os << endl;
        }
#       endif

      }

      else {
        break;
      }
    }

    max_likelihood = D_INF;
    for (i = 1;i <= max_nb_segment;i++) {
      if (penalized_likelihood[i] > max_likelihood) {
        max_likelihood = penalized_likelihood[i];
        nb_segment = i;
      }
    }

    psegment = seq->int_sequence[0][0];
    for (i = 0;i < nb_segment;i++) {
      for (j = change_point[nb_segment][i];j < change_point[nb_segment][i + 1];j++) {
        *psegment++ = i;
      }
    }

#   ifdef MESSAGE
    {
      int width[6];
      long old_adjust;
      double norm , *weight;


      weight = new double[max_nb_segment + 1];

      old_adjust = os.flags(ios::adjustfield);

      norm = 0.;
      for (i = 1;i <= max_nb_segment;i++) {
        weight[i] = exp((penalized_likelihood[i] - max_likelihood) / 2);
        norm += weight[i];
      }
      for (i = 1;i <= max_nb_segment;i++) {
        weight[i] /= norm;
      }

      for (i = 1;i <= max_nb_segment;i++) {
        penalty[i] += nb_parameter[i] * log((double)seq->length[0]);
      }

      width[0] = column_width(max_nb_segment) + ASCII_SPACE;
      width[1] = column_width(max_nb_segment , likelihood + 1 , 2.) + ASCII_SPACE;
      width[2] = column_width(nb_parameter[max_nb_segment]) + ASCII_SPACE;
      width[3] = column_width(max_nb_segment , penalty + 1) + ASCII_SPACE;
      width[4] = column_width(max_nb_segment , penalized_likelihood + 1) + ASCII_SPACE;
      width[5] = column_width(max_nb_segment , weight + 1) + ASCII_SPACE;

      os << "\n" << SEQ_label[SEQL_NB_SEGMENT] << " | 2 * " << STAT_label[STATL_LIKELIHOOD]
         << " | " << STAT_label[STATL_FREE_PARAMETERS] << " | " << SEQ_label[SEQL_PENALTY]
         << " | Modified " << STAT_criterion_word[BIC] << " - "  << STAT_label[STATL_WEIGHT] << endl;

      os.setf(ios::left , ios::adjustfield);

      for (i = 1;i <= max_nb_segment;i++) {
        os << setw(width[0]) << i
           << setw(width[1]) << 2 * likelihood[i]
           << setw(width[2]) << nb_parameter[i]
           << setw(width[3]) << penalty[i]
           << setw(width[4]) << penalized_likelihood[i]
           << setw(width[5]) << weight[i] << endl;
      }
      os << endl;

      os.setf((FMTFLAGS)old_adjust , ios::adjustfield);

/*      for (i = 1;i <= max_nb_segment;i++) {
        os << "\n" << i << " " << (i == 1 ? SEQ_label[SEQL_SEGMENT] : SEQ_label[SEQL_SEGMENTS])
           << "   2 * " << STAT_label[STATL_LIKELIHOOD] << ": " << 2 * segmentation_likelihood[i] << "   "
           << nb_parameter[i] << " " << STAT_label[nb_parameter[i] == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS]
           << "   2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " (Modified "  << STAT_criterion_word[BIC] << "): "
           << penalized_likelihood[1][i] << "   " << STAT_label[STATL_WEIGHT] << ": " << weight[1][i] << endl;
      } */

      delete [] weight;
    }
#   endif

    for (i = 1;i < seq->nb_variable;i++) {
      delete [] rank[i];
    }
    delete [] rank;

    delete [] frequency;

    for (i = 1;i < seq->nb_variable;i++) {
      if (model_type[i - 1] == CATEGORICAL_CHANGE) {
        for (j = 0;j < seq->length[0];j++) {
          delete [] begin_frequency[i][j];
          delete [] end_frequency[i][j];
        }
      }
      delete [] begin_frequency[i];
      delete [] end_frequency[i];

      delete [] begin_sum[i];
      delete [] end_sum[i];

      delete [] begin_factorial_sum[i];
      delete [] end_factorial_sum[i];

      delete [] begin_square_sum[i];
      delete [] end_square_sum[i];

      delete [] factorial[i];
    }
    delete [] begin_frequency;
    delete [] end_frequency;

    delete [] begin_sum;
    delete [] end_sum;

    delete [] begin_factorial_sum;
    delete [] end_factorial_sum;

    delete [] begin_square_sum;
    delete [] end_square_sum;

    delete [] factorial;

    delete [] mean;
    delete [] residual;

    delete [] begin_contrast;
    delete [] end_contrast;

    for (i = 1;i <= max_nb_segment;i++) {
      delete [] change_point[i];
    }
    delete [] change_point;

    delete [] likelihood;
    delete [] nb_parameter;
    delete [] penalty;
    delete [] penalized_likelihood;
  }

  return seq;
}


};  // namespace sequence_analysis
