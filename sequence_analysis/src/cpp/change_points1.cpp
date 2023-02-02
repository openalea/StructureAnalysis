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
 *       $Id: change_points1.cpp 18669 2015-11-09 12:08:08Z guedon $
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
#include <boost/math/special_functions/binomial.hpp>
#include <boost/math/distributions/normal.hpp>

#include "sequences.h"
#include "sequence_label.h"

using namespace std;
using namespace boost::math;
using namespace stat_tool;


namespace sequence_analysis {



/*--------------------------------------------------------------*/
/**
 *  \brief Determination of the width of a column of reals.
 *
 *  \param[in] nb_value number of values,
 *  \param[in] value    pointer on real values.
 *
 *  \return             column width.
 */
/*--------------------------------------------------------------*/

int column_width(int nb_value , const long double *value)

{
  int i;
  int width , max_width = 0;


  for (i = 0;i < nb_value;i++) {
    ostringstream ostring;
    ostring << value[i];
    width = (ostring.str()).size();
    if (width > max_width) {
      max_width = width;
    }
  }

  return max_width;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the log of the factorial of a value.
 *
 *  \param[in] value value.
 *
 *  \return          log of the factorial of a value.
 */
/*--------------------------------------------------------------*/

double log_factorial(int value)

{
  int i;
  double log_factorial;


  log_factorial = 0.;
  for (i = 2;i <= value;i++) {
    log_factorial += log((double)i);
  }

  return log_factorial;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the log of a binomial coefficient for
 *         a negative binomial distribution.
 *
 *  \param[in] inf_bound inf bound to the support,
 *  \param[in] parameter shape parameter,
 *  \param[in] value     value.
 *
 *  \return              log of a binomial coefficient.
 */
/*--------------------------------------------------------------*/

double log_binomial_coefficient(int inf_bound , double parameter , int value)

{
  int i;
  double set , subset , log_coeff;


  subset = parameter - 1.;
  set = subset;
  log_coeff = 0.;

  for (i = inf_bound;i < value;i++) {
    set++;
    log_coeff += log(set / (set - subset));
  }

# ifdef MESSAGE
  if (parameter == (int)parameter) {
    double ilog_coeff = log(binomial_coefficient<double>(value - inf_bound + parameter - 1 , parameter - 1));

    if ((log_coeff < ilog_coeff - DOUBLE_ERROR) || (log_coeff > ilog_coeff + DOUBLE_ERROR)) {
      cout << "TEST binomial coeff: " << log_coeff << " " << ilog_coeff << endl;
    }
  }
# endif

  return (log_coeff);
}


/*--------------------------------------------------------------*/
/**
 *  \brief Empirical determination of the hyperparameters of a gamma prior
 *         distribution for a Poisson distribution.
 *
 *  \param[in] index      sequence index,
 *  \param[in] variable   variable index,
 *  \param[in] hyperparam pointer on the hyperparameters.
 */
/*--------------------------------------------------------------*/

void Sequences::gamma_hyperparameter_computation(int index , int variable ,
                                                 double *hyperparam) const

{
  int i;
  double mean , diff , variance;


  if (length[index] > 1) {
    mean = 0.;
    for (i = 0;i < length[index];i++) {
      mean += int_sequence[index][variable][i];
    }
    mean /= length[index];

    variance = 0.;
    for (i = 0;i < length[i];i++) {
      diff = int_sequence[index][variable][i] - mean;
      variance += diff * diff;
    }
    variance /= (length[index] - 1);

    hyperparam[0] = mean * mean / (variance * PRIOR_VARIANCE_FACTOR);
    hyperparam[1] = mean / (variance * PRIOR_VARIANCE_FACTOR);

#   ifdef MESSAGE
    hyperparam[0] = 1.;
    hyperparam[1] = 200. / 365.;
#   endif

  }

  else {
    hyperparam[0] = D_DEFAULT;
    hyperparam[1] = D_DEFAULT;
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Empirical determination of the hyperparameters of a Gaussian-gamma prior
 *         distribution for a Gaussian distribution.
 *
 *  \param[in] index      sequence index,
 *  \param[in] variable   variable index,
 *  \param[in] hyperparam pointer on the hyperparameters.
 */
/*--------------------------------------------------------------*/

void Sequences::gaussian_gamma_hyperparameter_computation(int index , int variable ,
                                                          double *hyperparam) const

{
  int i;
  int magnitude;
  double mean , diff , dispersion , round_factor;


  if (length[index] > 1) {
    mean = 0.;

    if (type[variable] != REAL_VALUE) {
      for (i = 0;i < length[index];i++) {
        mean += int_sequence[index][variable][i];
      }
    }
    else {
      for (i = 0;i < length[index];i++) {
        mean += real_sequence[index][variable][i];
      }
    }
    mean /= length[index];

    dispersion = 0.;

    if (type[variable] != REAL_VALUE) {
      for (i = 1;i < length[index];i++) {
        diff = int_sequence[index][variable][i] - int_sequence[index][variable][i - 1];
        dispersion += diff * diff;
      }
    }
    else {
      for (i = 1;i < length[index];i++) {
        diff = real_sequence[index][variable][i] - real_sequence[index][variable][i - 1];
        dispersion += diff * diff;
      }
    }
    dispersion /= (2 * (length[index] - 1));

    hyperparam[0] = mean;
    hyperparam[1] = PRIOR_SAMPLE_SIZE;
    hyperparam[2] = PRIOR_DEGREES_OF_FREEDOM;
    hyperparam[3] = dispersion / PRIOR_DISPERSION_FACTOR;

    magnitude = (int)(log10(hyperparam[0])) + 1;
    if (magnitude > PRIOR_PRECISION) {
      round_factor = pow(10.0 , magnitude - PRIOR_PRECISION);

#     ifdef DEBUG
      cout << "\nTEST 0: " << magnitude << " " << round_factor << endl;
#     endif

      hyperparam[0] = round_factor * ::round(hyperparam[0] / round_factor);
    }

    magnitude = (int)(log10(hyperparam[3])) + 1;
    if (magnitude > PRIOR_PRECISION) {
      round_factor = pow(10.0 , magnitude - PRIOR_PRECISION);

#     ifdef DEBUG
      cout << "\nTEST 3: " << magnitude << " " << round_factor << endl;
#     endif

      hyperparam[3] = round_factor * ::round(hyperparam[3] / round_factor);
    }

//    hyperparam[3] /= 10;
  }

  else {
    hyperparam[0] = D_DEFAULT;
    hyperparam[1] = D_DEFAULT;
    hyperparam[2] = D_DEFAULT;
    hyperparam[3] = D_DEFAULT;
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the number of free parameters.
 *
 *  \param[in] index           sequence index,
 *  \param[in] nb_segment      number of segments,
 *  \param[in] model_type      segment model types,
 *  \param[in] common_contrast flag contrast functions common to the individuals.
 *
 *  \return                    number of free parameters.
 */
/*--------------------------------------------------------------*/

int Sequences::nb_parameter_computation(int index , int nb_segment , segment_model *model_type ,
                                        bool common_contrast) const

{
  bool *used_output;
  int i , j , k , m;
  int nb_parameter , max_nb_value;


//  nb_parameter = 0;
  nb_parameter = nb_segment - 1;

  if (model_type[0] == MEAN_CHANGE) {
    if ((index != I_DEFAULT) || (common_contrast)) {
      nb_parameter += nb_segment + 1;
    }
    else {
      nb_parameter += nb_sequence * nb_segment + 1;
    }
  }

  else if (model_type[0] == INTERCEPT_SLOPE_CHANGE) {
    if ((index != I_DEFAULT) || (common_contrast)) {
      nb_parameter += nb_segment * 2 + 1;
    }
    else {
      nb_parameter += nb_sequence * nb_segment * 2 + 1;
    }
  }

  else {
    max_nb_value = 0;
    for (i = 1;i < nb_variable;i++) {
      if ((model_type[i - 1] == CATEGORICAL_CHANGE) && (marginal_distribution[i]->nb_value > max_nb_value)) {
        max_nb_value = marginal_distribution[i]->nb_value;
      }
    }

    if (max_nb_value > 0) {
      used_output = new bool[max_nb_value];
    }
    else {
      used_output = NULL;
    }

    for (i = 1;i < nb_variable;i++) {
      if (model_type[i - 1] == CATEGORICAL_CHANGE) {
        if ((index != I_DEFAULT) || (!common_contrast)) {
          for (j = 0;j < nb_sequence;j++) {
            if ((index == I_DEFAULT) || (index == j)) {
              for (k = 0;k < length[j];k++) {
                if ((k == 0) || ((k > 0) && (int_sequence[j][0][k] != int_sequence[j][0][k - 1]))) {
                  for (m = 0;m < marginal_distribution[i]->nb_value;m++) {
                    used_output[m] = false;
                  }
                  nb_parameter--;
                }

                if (!used_output[int_sequence[j][i][k]]) {
                  nb_parameter++;
                  used_output[int_sequence[j][i][k]] = true;
                }
              }
            }
          }
        }

        else {
          for (j = 0;j < length[0];j++) {
            if ((j == 0) || ((j > 0) && (int_sequence[0][0][j] != int_sequence[0][0][j - 1]))) {
              for (k = 0;k < marginal_distribution[i]->nb_value;k++) {
                used_output[k] = false;
              }
              nb_parameter--;
            }

            for (k = 0;k < nb_sequence;k++) {
              if (!used_output[int_sequence[k][i][j]]) {
                nb_parameter++;
                used_output[int_sequence[k][i][j]] = true;
              }
            }

          }
        }
      }

      else if ((model_type[i - 1] == POISSON_CHANGE) || (model_type[i - 1] == NEGATIVE_BINOMIAL_0_CHANGE) ||
               (model_type[i - 1] == NEGATIVE_BINOMIAL_1_CHANGE) || (model_type[i - 1] == BAYESIAN_POISSON_CHANGE)) {
        if ((index != I_DEFAULT) || (common_contrast)) {
          nb_parameter += nb_segment;
        }
        else {
          nb_parameter += nb_sequence * nb_segment;
        }
      }

      else if ((model_type[i - 1] == GAUSSIAN_CHANGE) || (model_type[i - 1] == ORDINAL_GAUSSIAN_CHANGE) ||
               (model_type[i - 1] == BAYESIAN_GAUSSIAN_CHANGE)) {
        if ((index != I_DEFAULT) || (common_contrast)) {
          nb_parameter += 2 * nb_segment;
        }
        else {
          nb_parameter += nb_sequence * 2 * nb_segment;
        }
      }

      else if (model_type[i - 1] == VARIANCE_CHANGE) {
        if ((index != I_DEFAULT) || (common_contrast)) {
          nb_parameter += nb_segment + 1;
        }
        else {
          nb_parameter += nb_sequence * (nb_segment + 1);
        }
      }

      else if ((model_type[i - 1] == LINEAR_MODEL_CHANGE) || (model_type[i - 1] == AUTOREGRESSIVE_MODEL_CHANGE)) {
        if ((index != I_DEFAULT) || (common_contrast)) {
          nb_parameter += 3 * nb_segment;
        }
        else {
          nb_parameter += nb_sequence * 3 * nb_segment;
        }
      }

      else if (model_type[i - 1] == STATIONARY_AUTOREGRESSIVE_MODEL_CHANGE) {
        if ((index != I_DEFAULT) || (common_contrast)) {
          nb_parameter += 2 * nb_segment + 1;
        }
        else {
          nb_parameter += nb_sequence * (2 * nb_segment + 1);
        }
      }
    }

    delete [] used_output;
  }

  return nb_parameter;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the log-likelihood in the case of a single segment.
 *
 *  \param[in] index           sequence index,
 *  \param[in] model_type      segment model types,
 *  \param[in] common_contrast flag contrast functions common to the individuals,
 *  \param[in] shape_parameter negative binomial shape parameters,
 *  \param[in] rank            ranks (for ordinal variables).
 *
 *  \return                    log-likelihood of the single-segment model.
 */
/*--------------------------------------------------------------*/

double Sequences::one_segment_likelihood(int index , segment_model *model_type , bool common_contrast ,
                                         double *shape_parameter , double **rank) const

{
  int i , j , k;
  int max_nb_value , seq_length , count , *frequency , *inf_bound_parameter , *seq_index_parameter;
  double sum , factorial_sum , binomial_coeff_sum , proba , mean , diff , index_parameter_mean ,
         index_parameter_diff , index_parameter_sum , shifted_diff , likelihood;
  long double index_parameter_square_sum , square_sum , mix_square_sum , shifted_square_sum ,
              autocovariance , *residual;


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
      residual = new long double[nb_sequence];
    }
  }

  if (max_nb_value > 0) {
    frequency = new int[max_nb_value];
  }
  else {
    frequency = NULL;
  }

  seq_length = length[index == I_DEFAULT ? 0 : index];
  seq_index_parameter = NULL;

  for (i = 1;i < nb_variable;i++) {
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
  }

  if ((model_type[0] != MEAN_CHANGE) && (model_type[0] != INTERCEPT_SLOPE_CHANGE)) {
    likelihood = 0.;
  }

  for (i = 1;i < nb_variable;i++) {
    if (model_type[i - 1] == CATEGORICAL_CHANGE) {
      if ((index != I_DEFAULT) || (!common_contrast)) {
        for (j = 0;j < nb_sequence;j++) {
          if ((index == I_DEFAULT) || (index == j)) {
            for (k = 0;k < marginal_distribution[i]->nb_value;k++) {
              frequency[k] = 0;
            }

            for (k = 0;k < length[j];k++) {
              frequency[int_sequence[j][i][k]]++;
            }

            for (k = 0;k < marginal_distribution[i]->nb_value;k++) {
              if (frequency[k] > 0) {
                likelihood += frequency[k] * log((double)frequency[k] / (double)length[j]);
              }
            }
          }
        }
      }

      else {
        for (j = 0;j < marginal_distribution[i]->nb_value;j++) {
          frequency[j] = 0;
        }

        for (j = 0;j < length[0];j++) {
          for (k = 0;k < nb_sequence;k++) {
            frequency[int_sequence[k][i][j]]++;
          }
        }

        for (j = 0;j < marginal_distribution[i]->nb_value;j++) {
          if (frequency[j] > 0) {
            likelihood += frequency[j] * log((double)frequency[j] / (double)(nb_sequence * length[0]));
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
            for (k = 0;k < length[j];k++) {
              sum += int_sequence[j][i][k];
              factorial_sum += log_factorial(int_sequence[j][i][k]);
            }

            if (sum > 0.) {
              likelihood += sum * (log(sum / length[j]) - 1) - factorial_sum;
            }
          }
        }
      }

      else {
        sum = 0.;
        factorial_sum = 0.;
        for (j = 0;j < length[0];j++) {
          for (k = 0;k < nb_sequence;k++) {
            sum += int_sequence[k][i][j];
            factorial_sum += log_factorial(int_sequence[k][i][j]);
          }
        }

        if (sum > 0.) {
          likelihood += sum * (log(sum / (nb_sequence * length[0])) - 1) - factorial_sum;
        }
      }
    }

    else if ((model_type[i - 1] == NEGATIVE_BINOMIAL_0_CHANGE) || (model_type[i - 1] == NEGATIVE_BINOMIAL_1_CHANGE)) {
      if ((index != I_DEFAULT) || (!common_contrast)) {
        for (j = 0;j < nb_sequence;j++) {
          if ((index == I_DEFAULT) || (index == j)) {
            sum = 0.;
            binomial_coeff_sum = 0.;
            for (k = 0;k < length[j];k++) {
              sum += int_sequence[j][i][k];
              binomial_coeff_sum += log_binomial_coefficient(inf_bound_parameter[i - 1] , shape_parameter[i - 1] ,
                                                             int_sequence[j][i][k]);
            }

            if (sum > inf_bound_parameter[i - 1] * length[j]) {
              proba = shape_parameter[i - 1] * length[j] /
                      ((shape_parameter[i - 1] - inf_bound_parameter[i - 1]) * length[j] + sum);
              likelihood += binomial_coeff_sum + shape_parameter[i - 1] * length[j] * log(proba) +
                            (sum - inf_bound_parameter[i - 1] * length[j]) * log(1. - proba);
            }
            else {
              likelihood = D_INF;
              break;
            }
          }
        }
      }

      else {
        sum = 0.;
        binomial_coeff_sum = 0.;
        for (j = 0;j < length[0];j++) {
          for (k = 0;k < nb_sequence;k++) {
            sum += int_sequence[k][i][j];
            binomial_coeff_sum += log_binomial_coefficient(inf_bound_parameter[i - 1] , shape_parameter[i - 1] ,
                                                           int_sequence[k][i][j]);
          }
        }

        if (sum > inf_bound_parameter[i - 1] * nb_sequence * length[0]) {
          proba = shape_parameter[i - 1] * nb_sequence * length[0] /
                  ((shape_parameter[i - 1] - inf_bound_parameter[i - 1]) * nb_sequence * length[0] + sum);
          likelihood += binomial_coeff_sum + shape_parameter[i - 1] * nb_sequence * length[0] * log(proba) +
                        (sum - inf_bound_parameter[i - 1] * nb_sequence * length[0]) * log(1. - proba);
        }
        else {
          likelihood = D_INF;
          break;
        }
      }
    }

    else if ((model_type[i - 1] == GAUSSIAN_CHANGE) || (model_type[0] == MEAN_CHANGE) ||
             (model_type[i - 1] == VARIANCE_CHANGE)) {
      if ((index != I_DEFAULT) || (!common_contrast)) {
        if (type[i] != REAL_VALUE) {
          for (j = 0;j < nb_sequence;j++) {
            if ((index == I_DEFAULT) || (index == j)) {
/*              residual[j] = 0.;
              sum = int_sequence[j][i][0];

              for (k = 1;k < length[j];k++) {
                diff = int_sequence[j][i][k] - sum / k;
                residual[j] += ((double)k / (double)(k + 1)) * diff * diff;
                sum += int_sequence[j][i][k];
              } */

              mean = 0.;
              for (k = 0;k < length[j];k++) {
                mean += int_sequence[j][i][k];
              }
              mean /= length[j];

              residual[j] = 0.;
              for (k = 0;k < length[j];k++) {
                diff = int_sequence[j][i][k] - mean;
                residual[j] += diff * diff;
              }
            }
          }
        }

        else {
          for (j = 0;j < nb_sequence;j++) {
            if ((index == I_DEFAULT) || (index == j)) {
/*              residual[j] = 0.;
              sum = real_sequence[j][i][0];

              for (k = 1;k < length[j];k++) {
                diff = real_sequence[j][i][k] - sum / k;
                residual[j] += ((double)k / (double)(k + 1)) * diff * diff;
                sum += real_sequence[j][i][k];
              } */

              mean = 0.;
              for (k = 0;k < length[j];k++) {
                mean += real_sequence[j][i][k];
              }
              mean /= length[j];

              residual[j] = 0.;
              for (k = 0;k < length[j];k++) {
                diff = real_sequence[j][i][k] - mean;
                residual[j] += diff * diff;
              }
            }
          }
        }
      }

      else {
        if (type[i] != REAL_VALUE) {
/*          residual[0] = 0.;
          sum = 0.;
          count = 0;

          for (j = 0;j < length[0];j++) {
            for (k = 0;k < nb_sequence;k++) {
              if (count > 0) {
                diff = int_sequence[k][i][j] - sum / count;
                residual[0] += ((double)count / (double)(count + 1)) * diff * diff;
              }
              count++;
              sum += int_sequence[k][i][j];
            }
          } */

          mean = 0.;
          for (j = 0;j < length[0];j++) {
            for (k = 0;k < nb_sequence;k++) {
              mean += int_sequence[k][i][j];
            }
          }
          mean /= nb_sequence * length[0];

          residual[0] = 0.;
          for (j = 0;j < length[0];j++) {
            for (k = 0;k < nb_sequence;k++) {
              diff = int_sequence[k][i][j] - mean;
              residual[0] += diff * diff;
            }
          }
        }

        else {
/*          residual[0] = 0.;
          sum = 0.;
          count = 0;

          for (j = 0;j < length[0];j++) {
            for (k = 0;k < nb_sequence;k++) {
              if (count > 0) {
                diff = real_sequence[k][i][j] - sum / count;
                residual[0] += ((double)count / (double)(count + 1)) * diff * diff;
              }
              count++;
              sum += real_sequence[k][i][j];
            }
          } */

          mean = 0.;
          for (j = 0;j < length[0];j++) {
            for (k = 0;k < nb_sequence;k++) {
              mean += real_sequence[k][i][j];
            }
          }
          mean /= nb_sequence * length[0];

          residual[0] = 0.;
          for (j = 0;j < length[0];j++) {
            for (k = 0;k < nb_sequence;k++) {
              diff = real_sequence[k][i][j] - mean;
              residual[0] += diff * diff;
            }
          }
        }
      }
    }

    else if (model_type[i - 1] == ORDINAL_GAUSSIAN_CHANGE) {
      if ((index != I_DEFAULT) || (!common_contrast)) {
        for (j = 0;j < nb_sequence;j++) {
          if ((index == I_DEFAULT) || (index == j)) {
/*            residual[j] = 0.;
            sum = rank[i][int_sequence[j][i][0]];

            for (k = 1;k < length[j];k++) {
              diff = rank[i][int_sequence[j][i][k]] - sum / k;
              residual[j] += ((double)k / (double)(k + 1)) * diff * diff;
              sum += rank[i][int_sequence[j][i][k]];
            } */

            mean = 0.;
            for (k = 0;k < length[j];k++) {
              mean += rank[i][int_sequence[j][i][k]];
            }
            mean /= length[j];

            residual[j] = 0.;
            for (k = 0;k < length[j];k++) {
              diff = rank[i][int_sequence[j][i][k]] - mean;
              residual[j] += diff * diff;
            }
          }
        }
      }

      else {
/*        residual[0] = 0.;
        sum = 0.;
        count = 0;

        for (j = 0;j < length[0];j++) {
          for (k = 0;k < nb_sequence;k++) {
            if (count > 0) {
              diff = rank[i][int_sequence[k][i][j]] - sum / count;
              residual[0] += ((double)count / (double)(count + 1)) * diff * diff;
            }
            count++;
            sum += rank[i][int_sequence[k][i][j]];
          }
        } */

        mean = 0.;
        for (j = 0;j < length[0];j++) {
          for (k = 0;k < nb_sequence;k++) {
            mean += rank[i][int_sequence[k][i][j]];
          }
        }
        mean /= nb_sequence * length[0];

        residual[0] = 0.;
        for (j = 0;j < length[0];j++) {
          for (k = 0;k < nb_sequence;k++) {
            diff = rank[i][int_sequence[k][i][j]] - mean;
            residual[0] += diff * diff;
          }
        }
      }
    }

    else if ((model_type[i - 1] == LINEAR_MODEL_CHANGE) || (model_type[0] == INTERCEPT_SLOPE_CHANGE)) {
      if ((index != I_DEFAULT) || (!common_contrast)) {
        if (type[i] != REAL_VALUE) {
          for (j = 0;j < nb_sequence;j++) {
            if ((index == I_DEFAULT) || (index == j)) {
/*              index_parameter_square_sum = 0.;
              square_sum = 0.;
              mix_square_sum = 0.;
              index_parameter_sum = seq_index_parameter[0];
              sum = int_sequence[j][i][0];

              for (k = 1;k < length[j];k++) {
                index_parameter_diff = seq_index_parameter[k] - index_parameter_sum / k;
                index_parameter_square_sum += ((double)k / (double)(k + 1)) *
                                               index_parameter_diff * index_parameter_diff;
                diff = int_sequence[j][i][k] - sum / k;
                square_sum += ((double)k / (double)(k + 1)) * diff * diff;
                mix_square_sum += ((double)k / (double)(k + 1)) * index_parameter_diff * diff;
                index_parameter_sum += seq_index_parameter[k];
                sum += int_sequence[j][i][k];
              } */

              index_parameter_mean = 0.;
              mean = 0.;
              for (k = 0;k < length[j];k++) {
                index_parameter_mean += seq_index_parameter[k];
                mean += int_sequence[j][i][k];
              }
              index_parameter_mean /= length[j];
              mean /= length[j];

              index_parameter_square_sum = 0.;
              square_sum = 0.;
              mix_square_sum = 0.;
              for (k = 0;k < length[j];k++) {
                index_parameter_diff = seq_index_parameter[k] - index_parameter_mean;
                diff = int_sequence[j][i][k] - mean;
                index_parameter_square_sum += index_parameter_diff * index_parameter_diff;
                square_sum += diff * diff;
                mix_square_sum += index_parameter_diff * diff;
              }

              if (index_parameter_square_sum > 0.) {
                residual[j] = square_sum - mix_square_sum * mix_square_sum / index_parameter_square_sum;
              }
              else {
                residual[j] = 0.;
              }
            }
          }
        }

        else {
          for (j = 0;j < nb_sequence;j++) {
            if ((index == I_DEFAULT) || (index == j)) {
/*              index_parameter_square_sum = 0.;
              square_sum = 0.;
              mix_square_sum = 0.;
              index_parameter_sum = seq_index_parameter[0];
              sum = real_sequence[j][i][0];

              for (k = 1;k < length[j];k++) {
                index_parameter_diff = seq_index_parameter[k] - index_parameter_sum / k;
                index_parameter_square_sum += ((double)k / (double)(k + 1)) *
                                               index_parameter_diff * index_parameter_diff;
                diff = real_sequence[j][i][k] - sum / k;
                square_sum += ((double)k / (double)(k + 1)) * diff * diff;
                mix_square_sum += ((double)k / (double)(k + 1)) * index_parameter_diff * diff;
                index_parameter_sum += seq_index_parameter[k];
                sum += real_sequence[j][i][k];
              } */

              index_parameter_mean = 0.;
              mean = 0.;
              for (k = 0;k < length[j];k++) {
                index_parameter_mean += seq_index_parameter[k];
                mean += real_sequence[j][i][k];
              }
              index_parameter_mean /= length[j];
              mean /= length[j];

              index_parameter_square_sum = 0.;
              square_sum = 0.;
              mix_square_sum = 0.;
              for (k = 0;k < length[j];k++) {
                index_parameter_diff = seq_index_parameter[k] - index_parameter_mean;
                diff = real_sequence[j][i][k] - mean;
                index_parameter_square_sum += index_parameter_diff * index_parameter_diff;
                square_sum += diff * diff;
                mix_square_sum += index_parameter_diff * diff;
              }

              if (index_parameter_square_sum > 0.) {
                residual[j] += square_sum - mix_square_sum * mix_square_sum / index_parameter_square_sum;
              }
              else {
                residual[j] = 0.;
              }
            }
          }
        }
      }

      else {
        if (type[i] != REAL_VALUE) {
/*          index_parameter_square_sum = 0.;
          square_sum = 0.;
          mix_square_sum = 0.;
          count = 1;

          index_parameter_sum = nb_sequence * seq_index_parameter[0];
          sum = int_sequence[0][i][0];
          for (j = 1;j < nb_sequence;j++) {
            diff = int_sequence[j][i][0] - sum / count;
            square_sum += ((double)count / (double)(count + 1)) * diff * diff;
            count++;
            sum += int_sequence[j][i][0];
          }

          for (j = 1;j < length[0];j++) {
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
          } */

          index_parameter_mean = 0.;
          mean = 0.;
          for (j = 0;j < length[0];j++) {
            index_parameter_mean += seq_index_parameter[j];
            for (k = 0;k < nb_sequence;k++) {
              mean += int_sequence[k][i][j];
            }
          }
          index_parameter_mean /= length[0];
          mean /= nb_sequence * length[0];

          index_parameter_square_sum = 0.;
          square_sum = 0.;
          mix_square_sum = 0.;
          for (j = 0;j < length[0];j++) {
            index_parameter_diff = seq_index_parameter[j] - index_parameter_mean;
            index_parameter_square_sum += index_parameter_diff * index_parameter_diff;
            for (k = 0;k < nb_sequence;k++) {
              diff = int_sequence[k][i][j] - mean;
              square_sum += diff * diff;
              mix_square_sum += index_parameter_diff * diff;
            }
          }
        }

        else {
/*          index_parameter_square_sum = 0.;
          square_sum = 0.;
          mix_square_sum = 0.;
          count = 1;

          index_parameter_sum = nb_sequence * seq_index_parameter[0];
          sum = real_sequence[0][i][0];
          for (j = 1;j < nb_sequence;j++) {
            diff = real_sequence[j][i][0] - sum / count;
            square_sum += ((double)count / (double)(count + 1)) * diff * diff;
            count++;
            sum += real_sequence[j][i][0];
          }

          for (j = 1;j < length[0];j++) {
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
          } */

          index_parameter_mean = 0.;
          mean = 0.;
          for (j = 0;j < length[0];j++) {
            index_parameter_mean += seq_index_parameter[j];
            for (k = 0;k < nb_sequence;k++) {
              mean += real_sequence[k][i][j];
            }
          }
          index_parameter_mean /= length[0];
          mean /= nb_sequence * length[0];

          index_parameter_square_sum = 0.;
          square_sum = 0.;
          mix_square_sum = 0.;
          for (j = 0;j < length[0];j++) {
            index_parameter_diff = seq_index_parameter[j] - index_parameter_mean;
            index_parameter_square_sum += index_parameter_diff * index_parameter_diff;
            for (k = 0;k < nb_sequence;k++) {
              diff = real_sequence[k][i][j] - mean;
              square_sum += diff * diff;
              mix_square_sum += index_parameter_diff * diff;
            }
          }
        }

        index_parameter_square_sum *= nb_sequence;
        if (index_parameter_square_sum > 0.) {
          residual[0] = square_sum - mix_square_sum * mix_square_sum / index_parameter_square_sum;
        }
        else {
          residual[0] = 0.;
        }
      }
    }

    else if ((model_type[i - 1] == AUTOREGRESSIVE_MODEL_CHANGE) ||
             (model_type[i - 1] == STATIONARY_AUTOREGRESSIVE_MODEL_CHANGE)) {
      if ((index != I_DEFAULT) || (!common_contrast)) {
        if (type[i] != REAL_VALUE) {
          for (j = 0;j < nb_sequence;j++) {
            if ((index == I_DEFAULT) || (index == j)) {
              mean = 0.;
              for (k = 0;k < length[j];k++) {
                mean += int_sequence[j][i][k];
              }
              mean /= length[j];

              square_sum = 0.;
              shifted_square_sum = 0.;
              autocovariance = 0.;
              for (k = 1;k < length[j];k++) {
                diff = int_sequence[j][i][k] - mean;
                shifted_diff = int_sequence[j][i][k - 1] - mean;
                square_sum += diff * diff;
                shifted_square_sum += shifted_diff * shifted_diff;
                autocovariance += diff * shifted_diff;
              }

              residual[j] = square_sum;
              if (shifted_square_sum > 0.) {
                residual[j] -= autocovariance * autocovariance / shifted_square_sum;
              }
            }
          }
        }

        else {
          for (j = 0;j < nb_sequence;j++) {
            if ((index == I_DEFAULT) || (index == j)) {
              mean = 0.;
              for (k = 0;k < length[j];k++) {
                mean += real_sequence[j][i][k];
              }
              mean /= length[j];

              square_sum = 0.;
              shifted_square_sum = 0.;
              autocovariance = 0.;
              for (k = 1;k < length[j];k++) {
                diff = real_sequence[j][i][k] - mean;
                shifted_diff = real_sequence[j][i][k - 1] - mean;
                square_sum += diff * diff;
                shifted_square_sum += shifted_diff * shifted_diff;
                autocovariance += diff * shifted_diff;
              }

              residual[j] = square_sum;
              if (shifted_square_sum > 0.) {
                residual[j] -= autocovariance * autocovariance / shifted_square_sum;
              }
            }
          }
        }

        for (j = 0;j < nb_sequence;j++) {
          if ((index == I_DEFAULT) || (index == j)) {
//            if (residual[j] > 0.) {
            if (residual[j] > (length[j] - 1) * ROUNDOFF_ERROR) {
              likelihood -= ((double)(length[j] - 1) / 2.) * (logl(residual[j] / (length[j] - 1)) +
                             log(2 * M_PI) + 1);
            }
            else {
              likelihood = D_INF;
              break;
            }
          }
        }
      }

      else {
        if (type[i] != REAL_VALUE) {
          mean = 0.;
          for (j = 0;j < length[0];j++) {
            for (k = 0;k < nb_sequence;k++) {
              mean += int_sequence[k][i][j];
            }
          }
          mean /= nb_sequence * length[0];

          square_sum = 0.;
          shifted_square_sum = 0.;
          autocovariance = 0.;
          for (j = 1;j < length[0];j++) {
            for (k = 0;k < nb_sequence;k++) {
              diff = int_sequence[k][i][j] - mean;
              shifted_diff = int_sequence[k][i][j - 1] - mean;
              square_sum += diff * diff;
              shifted_square_sum += shifted_diff * shifted_diff;
              autocovariance += diff * shifted_diff;
            }
          }
        }

        else {
          mean = 0.;
          for (j = 0;j < length[0];j++) {
            for (k = 0;k < nb_sequence;k++) {
              mean += real_sequence[k][i][j];
            }
          }
          mean /= nb_sequence * length[0];

          square_sum = 0.;
          shifted_square_sum = 0.;
          autocovariance = 0.;
          for (j = 1;j < length[0];j++) {
            for (k = 0;k < nb_sequence;k++) {
              diff = real_sequence[k][i][j] - mean;
              shifted_diff = real_sequence[k][i][j - 1] - mean;
              square_sum += diff * diff;
              shifted_square_sum += shifted_diff * shifted_diff;
              autocovariance += diff * shifted_diff;
            }
          }
        }

        residual[0] = square_sum;
        if (shifted_square_sum > 0.) {
          residual[0] -= autocovariance * autocovariance / shifted_square_sum;
        }

//        if (residual[0] > 0.) {
        if (residual[0] > nb_sequence * (length[0] - 1) * ROUNDOFF_ERROR) {
          likelihood -= ((double)(nb_sequence * (length[0] - 1)) / 2.) *
                        (logl(residual[0] / (nb_sequence * (length[0] - 1))) + log(2 * M_PI) + 1);
        }
        else {
          likelihood = D_INF;
          break;
        }
      }
    }

    if ((model_type[i - 1] == GAUSSIAN_CHANGE) || (model_type[i - 1] == VARIANCE_CHANGE) ||
        (model_type[i - 1] == ORDINAL_GAUSSIAN_CHANGE) || (model_type[i - 1] == LINEAR_MODEL_CHANGE)) {
      if ((index != I_DEFAULT) || (!common_contrast)) {
        for (j = 0;j < nb_sequence;j++) {
          if ((index == I_DEFAULT) || (index == j)) {
//            if (residual[j] > 0.) {
            if (residual[j] > length[j] * ROUNDOFF_ERROR) {
              likelihood -= ((double)length[j] / 2.) * (logl(residual[j] / length[j]) +
                             log(2 * M_PI) + 1);
/*              likelihood -= ((double)length[j] / 2.) * (logl(residual[j] / (length[j] - 1)) +
                             log(2 * M_PI)) - (double)(length[j] - 1) / 2.; */
            }
            else {
              likelihood = D_INF;
              break;
            }
          }
        }
      }

      else {
//        if (residual[0] > 0.) {
        if (residual[0] > nb_sequence * length[0] * ROUNDOFF_ERROR) {
          likelihood -= ((double)(nb_sequence * length[0]) / 2.) *
                        (logl(residual[0] / (nb_sequence * length[0])) + log(2 * M_PI) + 1);
/*          likelihood -= ((double)(nb_sequence * length[0]) / 2.) *
                        (logl(residual[0] / (nb_sequence * length[0] - 1)) +
                         log(2 * M_PI)) - (double)(nb_sequence * length[0] - 1) / 2.; */
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

  if ((model_type[0] == MEAN_CHANGE) || (model_type[0] == INTERCEPT_SLOPE_CHANGE)) {
    if (index != I_DEFAULT) {
//      if (residual[index] > 0.) {
      if (residual[index] > length[index] * ROUNDOFF_ERROR) {
        likelihood -= ((double)length[index] / 2.) * (logl(residual[index] / length[index]) +
                       log(2 * M_PI) + 1);
/*        likelihood -= ((double)length[index] / 2.) * (logl(residual[index] / (length[index] - 1)) +
                       log(2 * M_PI)) - (double)(length[index] - 1) / 2.; */
      }
      else {
        likelihood = D_INF;
      }
    }

    else {
      if (!common_contrast) {
        for (i = 1;i < nb_sequence;i++) {
          residual[0] += residual[i];
        }
      }

//      if (residual[0] > 0.) {
      if (residual[0] > nb_sequence * length[0] * ROUNDOFF_ERROR) {
        likelihood = -((double)(nb_sequence * length[0]) / 2.) *
                      (logl(residual[0] / (nb_sequence * length[0])) + log(2 * M_PI) + 1);
/*        likelihood = -((double)(nb_sequence * length[0]) / 2.) *
                      (logl(residual[0] / (nb_sequence * (length[0] - 1))) +
                       log(2 * M_PI)) - (double)(nb_sequence * (length[0] - 1)) / 2.; */
      }
      else {
        likelihood = D_INF;
      }
    }
  }

  for (i = 0;i < nb_sequence;i++) {
    if ((index == I_DEFAULT) || (index == i)) {
      for (j = 0;j < length[i];j++) {
        int_sequence[i][0][j] = 0;
      }
    }
  }

  delete [] frequency;
  delete [] residual;

  if (index_param_type == IMPLICIT_TYPE) {
    delete [] seq_index_parameter;
  }

  return likelihood;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of piecewise linear functions.
 *
 *  \param[in] index                     sequence index,
 *  \param[in] variable                  variable index,
 *  \param[in] nb_segment                number of segments,
 *  \param[in] model_type                segment model type,
 *  \param[in] common_contrast           flag contrast functions common to the individuals,
 *  \param[in] change_point              change points,
 *  \param[in] seq_index_parameter       index parameters,
 *  \param[in] piecewise_function        piecewise linear functions,
 *  \param[in] imean                     segment means,
 *  \param[in] variance                  segment variances or residual variances,
 *  \param[in] global_variance           global variance or residual variance,
 *  \param[in] iintercept                segment intercepts,
 *  \param[in] islope                    segment slopes,
 *  \param[in] autoregressive_coeff      segment autoregressive coefficient,
 *  \param[in] correlation               segment correlation coefficients (for linear models),
 *  \param[in] slope_standard_deviation  segment slope standard deviations (for linear models),
 *  \param[in] iindex_parameter_mean     segment index parameter mean (for linear models),
 *  \param[in] iindex_parameter_variance segment index parameter variance (for linear models),
 *  \param[in] determination_coeff       coefficient of determination (for autoregressive models).
 *
 *  \return                              log-likelihood of the piecewise linear function.
 */
/*--------------------------------------------------------------*/

double Sequences::piecewise_linear_function(int index , int variable , int nb_segment , segment_model model_type ,
                                            bool common_contrast , int *change_point , int *seq_index_parameter ,
                                            double **piecewise_function , double **imean , double **variance ,
                                            double *global_variance , double **iintercept , double **islope ,
                                            double **iautoregressive_coeff , double **correlation ,
                                            double **slope_standard_deviation , double **iindex_parameter_mean ,
                                            long double **iindex_parameter_variance , double **determination_coeff) const

{
  int i , j , k;
  double likelihood , mean , diff , diff_sum , index_parameter_mean , response_mean , shifted_diff ,
    slope , intercept , autoregressive_coeff , *individual_mean , *rank;
  long double square_sum , global_square_sum , index_parameter_variance , response_variance , covariance ,
              shifted_square_sum , autocovariance , residual_square_sum , mean_squared_error_1;

# ifdef MESSAGE
  long double mean_squared_error;
# endif


  if ((model_type == POISSON_CHANGE) || (model_type == NEGATIVE_BINOMIAL_0_CHANGE) ||
      (model_type == NEGATIVE_BINOMIAL_1_CHANGE) || (model_type == GAUSSIAN_CHANGE) ||
      (model_type == MEAN_CHANGE) || (model_type == VARIANCE_CHANGE) ||
      (model_type == BAYESIAN_POISSON_CHANGE) || (model_type == BAYESIAN_GAUSSIAN_CHANGE)) {
    if (((model_type == GAUSSIAN_CHANGE) || (model_type == VARIANCE_CHANGE)) && ((variance) || (global_variance))) {
      likelihood = 0.;
    }
    else {
      likelihood = D_INF;
    }

    if (global_variance) {
      global_square_sum = 0.;
    }

    if ((index != I_DEFAULT) || (!common_contrast)) {
      for (i = 0;i < nb_sequence;i++) {
        if ((index == I_DEFAULT) || (index == i)) {
          if (model_type == VARIANCE_CHANGE) {
            mean = 0.;

            if (type[variable] != REAL_VALUE) {
              for (j = 0;j < length[i];j++) {
                mean += int_sequence[i][variable][j];
              }
            }
            else {
              for (j = 0;j < length[i];j++) {
                mean += real_sequence[i][variable][j];
              }
            }
            mean /= length[i];
          }

          for (j = 0;j < nb_segment;j++) {
            if (model_type != VARIANCE_CHANGE) {
              mean = 0.;

              if (type[variable] != REAL_VALUE) {
                for (k = change_point[j];k < change_point[j + 1];k++) {
                  mean += int_sequence[i][variable][k];
                }
              }
              else {
                for (k = change_point[j];k < change_point[j + 1];k++) {
                  mean += real_sequence[i][variable][k];
                }
              }
              mean /= (change_point[j + 1] - change_point[j]);
            }

            if (imean) {
              imean[i][j] = mean;
            }
            if (piecewise_function) {
              for (k = change_point[j];k < change_point[j + 1];k++) {
                piecewise_function[i][k] = mean;
              }
            }

            if ((variance) || (global_variance)) {
              square_sum = 0.;

              if (change_point[j + 1] > change_point[j] + 1) {
                if (type[variable] != REAL_VALUE) {
                  for (k = change_point[j];k < change_point[j + 1];k++) {
                    diff = int_sequence[i][variable][k] - mean;
                    square_sum += diff * diff;
                  }
                }
                else {
                  for (k = change_point[j];k < change_point[j + 1];k++) {
                    diff = real_sequence[i][variable][k] - mean;
                    square_sum += diff * diff;
                  }
                }

                if (global_variance) {
                  global_square_sum += square_sum;
                }
                variance[i][j] = square_sum / ((change_point[j + 1] - change_point[j]) - 1);

                if (((model_type == GAUSSIAN_CHANGE) || (model_type == VARIANCE_CHANGE)) && (likelihood != D_INF)) {
                  if (square_sum > (change_point[j + 1] - change_point[j]) * ROUNDOFF_ERROR) {
                    likelihood -= ((double)(change_point[j + 1] - change_point[j]) / 2.) * (log(square_sum /
                                    (change_point[j + 1] - change_point[j])) + log(2 * M_PI) + 1);
                  }
                  else {
                    likelihood = D_INF;
                  }
                }
              }
            }
          }
        }
      }
    }

    else {
      individual_mean = new double[nb_sequence];

      // rank variance decomposition

      if (((model_type == POISSON_CHANGE) || (model_type == NEGATIVE_BINOMIAL_0_CHANGE) ||
           (model_type == NEGATIVE_BINOMIAL_1_CHANGE) || (model_type == BAYESIAN_POISSON_CHANGE)) && (variance)) {
        rank = marginal_distribution[variable]->rank_computation();

        for (i = 0;i < nb_segment;i++) {
          mean = 0.;
          for (j = 0;j < nb_sequence;j++) {
            individual_mean[j] = 0.;
            for (k = change_point[i];k < change_point[i + 1];k++) {
              individual_mean[j] += rank[int_sequence[j][variable][k]];
            }
            mean += individual_mean[j];
            individual_mean[j] /= (change_point[i + 1] - change_point[i]);
          }
          mean /= (nb_sequence * (change_point[i + 1] - change_point[i]));

          square_sum = 0.;
          for (j = 0;j < nb_sequence;j++) {
            for (k = change_point[i];k < change_point[i + 1];k++) {
              diff = rank[int_sequence[j][variable][k]] - mean;
              square_sum += diff * diff;
            }
          }
          variance[1][i] = square_sum / (nb_sequence * (change_point[i + 1] - change_point[i]));

          square_sum = 0.;
          for (j = 0;j < nb_sequence;j++) {
            diff = individual_mean[j] - mean;
            square_sum += diff * diff;
          }
          variance[2][i] = square_sum / nb_sequence;

          square_sum = 0.;
          for (j = 0;j < nb_sequence;j++) {
            for (k = change_point[i];k < change_point[i + 1];k++) {
              diff = rank[int_sequence[j][variable][k]] - individual_mean[j];
              square_sum += diff * diff;
            }
          }
          variance[3][i] = square_sum / (nb_sequence * (change_point[i + 1] - change_point[i]));
        }

        delete [] rank;
      }

      if (model_type == VARIANCE_CHANGE) {
        mean = 0.;

        if (type[variable] != REAL_VALUE) {
          for (i = 0;i < nb_sequence;i++) {
            individual_mean[i] = 0.;
            for (j = 0;j < length[0];j++) {
              individual_mean[i] += int_sequence[i][variable][j];
            }
            mean += individual_mean[i];
            individual_mean[i] /= length[0];
          }
        }
        else {
          for (i = 0;i < nb_sequence;i++) {
            individual_mean[i] = 0.;
            for (j = 0;j < length[0];j++) {
              individual_mean[i] += real_sequence[i][variable][j];
            }
            mean += individual_mean[i];
            individual_mean[i] /= length[0];
          }
        }
        mean /= (nb_sequence * length[0]);
      }

      for (i = 0;i < nb_segment;i++) {
        if (model_type != VARIANCE_CHANGE) {
          mean = 0.;

          if (type[variable] != REAL_VALUE) {
            for (j = 0;j < nb_sequence;j++) {
              individual_mean[j] = 0.;
              for (k = change_point[i];k < change_point[i + 1];k++) {
                individual_mean[j] += int_sequence[j][variable][k];
              }
              mean += individual_mean[j];
              individual_mean[j] /= (change_point[i + 1] - change_point[i]);
            }
          }
          else {
            for (j = 0;j < nb_sequence;j++) {
              individual_mean[j] = 0.;
              for (k = change_point[i];k < change_point[i + 1];k++) {
                individual_mean[j] += real_sequence[j][variable][k];
              }
              mean += individual_mean[j];
              individual_mean[j] /= (change_point[i + 1] - change_point[i]);
            }
          }
          mean /= (nb_sequence * (change_point[i + 1] - change_point[i]));
        }

        if (imean) {
          imean[0][i] = mean;
        }
        if (piecewise_function) {
          for (j = 0;j < nb_sequence;j++) {
            for (k = change_point[i];k < change_point[i + 1];k++) {
              piecewise_function[j][k] = mean;
            }
          }
        }

        if ((variance) || (global_variance)) {
          square_sum = 0.;

          if (type[variable] != REAL_VALUE) {
            for (j = 0;j < nb_sequence;j++) {
              for (k = change_point[i];k < change_point[i + 1];k++) {
                diff = int_sequence[j][variable][k] - mean;
                square_sum += diff * diff;
              }
            }
          }
          else {
            for (j = 0;j < nb_sequence;j++) {
              for (k = change_point[i];k < change_point[i + 1];k++) {
                diff = real_sequence[j][variable][k] - mean;
                square_sum += diff * diff;
              }
            }
          }

          if (global_variance) {
            global_square_sum += square_sum;
          }
          variance[0][i] = square_sum / (nb_sequence * (change_point[i + 1] - change_point[i]) - 1);

          if (((model_type == GAUSSIAN_CHANGE) || (model_type == VARIANCE_CHANGE)) && (likelihood != D_INF)) {
            if (square_sum > nb_sequence * (change_point[i + 1] - change_point[i]) * ROUNDOFF_ERROR) {
              likelihood -= ((double)(nb_sequence * (change_point[i + 1] - change_point[i])) / 2.) * (log(square_sum /
                              (nb_sequence * (change_point[i + 1] - change_point[i]))) + log(2 * M_PI) + 1);
            }
            else {
              likelihood = D_INF;
            }
          }

          // variance decomposition

          if ((model_type == GAUSSIAN_CHANGE) || (model_type == MEAN_CHANGE) ||
              (model_type == BAYESIAN_GAUSSIAN_CHANGE)) {
            variance[1][i] = square_sum / (nb_sequence * (change_point[i + 1] - change_point[i]));

            square_sum = 0.;
            for (j = 0;j < nb_sequence;j++) {
              diff = individual_mean[j] - mean;
              square_sum += diff * diff;
            }
            variance[2][i] = square_sum / nb_sequence;

            square_sum = 0.;

            if (type[variable] != REAL_VALUE) {
              for (j = 0;j < nb_sequence;j++) {
                for (k = change_point[i];k < change_point[i + 1];k++) {
                  diff = int_sequence[j][variable][k] - individual_mean[j];
                  square_sum += diff * diff;
                }
              }
            }
            else {
              for (j = 0;j < nb_sequence;j++) {
                for (k = change_point[i];k < change_point[i + 1];k++) {
                  diff = real_sequence[j][variable][k] - individual_mean[j];
                  square_sum += diff * diff;
                }
              }
            }

            variance[3][i] = square_sum / (nb_sequence * (change_point[i + 1] - change_point[i]));
          }
        }
      }

      delete [] individual_mean;
    }

    if (global_variance) {
      if (model_type == MEAN_CHANGE) {
        if (index != I_DEFAULT) {
          global_variance[variable] = global_square_sum / (length[index] - nb_segment);
        }
        else {
          global_variance[variable] = global_square_sum / (nb_sequence * length[0] - nb_segment);
        }

        if (index != I_DEFAULT) {
          if (global_square_sum > length[index] * ROUNDOFF_ERROR) {
            likelihood = -((double)length[index] / 2.) * (log(global_square_sum /
                            length[index]) + log(2 * M_PI) + 1);
          }
          else {
            likelihood = D_INF;
          }
        }

        else {
          if (global_square_sum > nb_sequence * length[0] * ROUNDOFF_ERROR) {
            likelihood = -((double)(nb_sequence * length[0]) / 2.) * (log(global_square_sum /
                            (nb_sequence * length[0])) + log(2 * M_PI) + 1);
          }
          else {
            likelihood = D_INF;
          }
        }
      }

      // computation of mean squared error

      else {
        if (index != I_DEFAULT) {
          global_variance[variable] = global_square_sum / length[index];
        }
        else {
          global_variance[variable] = global_square_sum / (nb_sequence * length[0]);
        }
      }
    }
  }

  else if ((model_type == LINEAR_MODEL_CHANGE) || (model_type == INTERCEPT_SLOPE_CHANGE)) {
    if ((model_type == LINEAR_MODEL_CHANGE) && ((variance) || (global_variance))) {
      likelihood = 0.;
    }
    else {
      likelihood = D_INF;
    }

    if (global_variance) {
      global_square_sum = 0.;
    }

    if ((index != I_DEFAULT) || (!common_contrast)) {
      for (i = 0;i < nb_sequence;i++) {
        if ((index == I_DEFAULT) || (index == i)) {
          for (j = 0;j < nb_segment;j++) {
            index_parameter_mean = 0.;
            for (k = change_point[j];k < change_point[j + 1];k++) {
              index_parameter_mean += seq_index_parameter[k];                
            }
            index_parameter_mean /= (change_point[j + 1] - change_point[j]);
            if (iindex_parameter_mean) {
              iindex_parameter_mean[i][j] = index_parameter_mean;
            }

            index_parameter_variance = 0.;
            for (k = change_point[j];k < change_point[j + 1];k++) {
              diff = seq_index_parameter[k] - index_parameter_mean;
              index_parameter_variance += diff * diff;                
            }
            if (iindex_parameter_variance) {
              iindex_parameter_variance[i][j] = index_parameter_variance;
            }

            response_mean = 0.;
            response_variance = 0.;
            covariance = 0.;

            if (type[variable] != REAL_VALUE) {
              for (k = change_point[j];k < change_point[j + 1];k++) {
                response_mean += int_sequence[i][variable][k];
              }
              response_mean /= (change_point[j + 1] - change_point[j]);

              for (k = change_point[j];k < change_point[j + 1];k++) {
                diff = int_sequence[i][variable][k] - response_mean;
                response_variance += diff * diff;
                covariance += (seq_index_parameter[k] - index_parameter_mean) * diff;
              }
            }

            else {
              for (k = change_point[j];k < change_point[j + 1];k++) {
                response_mean += real_sequence[i][variable][k];
              }
              response_mean /= (change_point[j + 1] - change_point[j]);

              for (k = change_point[j];k < change_point[j + 1];k++) {
                diff = real_sequence[i][variable][k] - response_mean;
                response_variance += diff * diff;
                covariance += (seq_index_parameter[k] - index_parameter_mean) * diff;
              }
            }

            slope = covariance / index_parameter_variance;
            intercept = response_mean - slope * index_parameter_mean;

            if ((islope) && (iintercept)) {
              iintercept[i][j] = intercept;
              islope[i][j] = slope;
            }
            if (correlation) {
              correlation[i][j] = covariance / sqrt(response_variance * index_parameter_variance);
            }

            if (piecewise_function) {
              for (k = change_point[j];k < change_point[j + 1];k++) {
                piecewise_function[i][k] = intercept + slope * seq_index_parameter[k];
              }
            }

            if ((variance) || (global_variance)) {
              if (change_point[j + 1] > change_point[j] + 2) {
                square_sum = 0.;

                if (type[variable] != REAL_VALUE) {
                  for (k = change_point[j];k < change_point[j + 1];k++) {
                    diff = int_sequence[i][variable][k] - (intercept + slope * seq_index_parameter[k]);
                    square_sum += diff * diff;
                  }
                }

                else {
                  for (k = change_point[j];k < change_point[j + 1];k++) {
                    diff = real_sequence[i][variable][k] - (intercept + slope * seq_index_parameter[k]);
                    square_sum += diff * diff;
                  }
                }

                if ((global_variance) || (model_type == INTERCEPT_SLOPE_CHANGE)) {
                  global_square_sum += square_sum;
                }

                variance[i][j] = square_sum / (change_point[j + 1] - change_point[j] - 2);

                if ((model_type == LINEAR_MODEL_CHANGE) && (likelihood != D_INF)) {
                  if (square_sum > (change_point[j + 1] - change_point[j]) * ROUNDOFF_ERROR) {
                    likelihood -= ((double)(change_point[j + 1] - change_point[j]) / 2.) * (log(square_sum /
                                    (change_point[j + 1] - change_point[j])) + log(2 * M_PI) + 1);
                  }
                  else {
                    likelihood = D_INF;
                  }
                }

                if (slope_standard_deviation) {
                  square_sum /= (change_point[j + 1] - change_point[j] - 2);
                  slope_standard_deviation[i][j] = sqrt(square_sum / index_parameter_variance);
                }
              }

              else {
                variance[i][j] = 0.;
              }
            }
          }
        }
      }
    }

    else {
      for (i = 0;i < nb_segment;i++) {
        index_parameter_mean = 0.;
        for (j = change_point[i];j < change_point[i + 1];j++) {
          index_parameter_mean += seq_index_parameter[j];                
        }
        index_parameter_mean /= (change_point[i + 1] - change_point[i]);
        if (iindex_parameter_mean) {
          iindex_parameter_mean[0][i] = index_parameter_mean;
        }

        index_parameter_variance = 0.;
        for (j = change_point[i];j < change_point[i + 1];j++) {
          diff = seq_index_parameter[j] - index_parameter_mean;
          index_parameter_variance += diff * diff;                
        }
        index_parameter_variance *= nb_sequence;
        if (iindex_parameter_variance) {
          iindex_parameter_variance[0][i] = index_parameter_variance;
        }

        response_mean = 0.;
        response_variance = 0.;
        covariance = 0.;

        if (type[variable] != REAL_VALUE) {
          for (j = change_point[i];j < change_point[i + 1];j++) {
            for (k = 0;k < nb_sequence;k++) {
              response_mean += int_sequence[k][variable][j];
            }
          }
          response_mean /= (nb_sequence * (change_point[i + 1] - change_point[i]));

          for (j = change_point[i];j < change_point[i + 1];j++) {
            diff_sum = 0.;
            for (k = 0;k < nb_sequence;k++) {
              diff = int_sequence[k][variable][j] - response_mean;
              response_variance += diff * diff;
              diff_sum += diff;
            }
            covariance += (seq_index_parameter[j] - index_parameter_mean) * diff_sum;
          }
        }

        else {
          for (j = change_point[i];j < change_point[i + 1];j++) {
            for (k = 0;k < nb_sequence;k++) {
              response_mean += real_sequence[k][variable][j];
            }
          }
          response_mean /= (nb_sequence * (change_point[i + 1] - change_point[i]));

          for (j = change_point[i];j < change_point[i + 1];j++) {
            diff_sum = 0.;
            for (k = 0;k < nb_sequence;k++) {
              diff = real_sequence[k][variable][j] - response_mean;
              response_variance += diff * diff;
              diff_sum += diff;
            }
            covariance += (seq_index_parameter[j] - index_parameter_mean) * diff_sum;
          }
        }

        slope = covariance / index_parameter_variance;
        intercept = response_mean - slope * index_parameter_mean;

        if ((islope) && (iintercept)) {
          iintercept[0][i] = intercept;
          islope[0][i] = slope;
        }
        if (correlation) {
          correlation[0][i] = covariance / sqrt(response_variance * index_parameter_variance);
        }

        if (piecewise_function) {
          for (j = 0;j < nb_sequence;j++) {
            for (k = change_point[i];k < change_point[i + 1];k++) {
              piecewise_function[j][k] = intercept + slope * seq_index_parameter[k];
            }
          }
        }

        if ((variance) || (global_variance)) {
          if (nb_sequence * (change_point[i + 1] - change_point[i]) > 2) {
            square_sum = 0.;

            if (type[variable] != REAL_VALUE) {
              for (j = change_point[i];j < change_point[i + 1];j++) {
                for (k = 0;k < nb_sequence;k++) {
                  diff = int_sequence[k][variable][j] - (intercept + slope * seq_index_parameter[j]);
                  square_sum += diff * diff;
                }
              }
            }

            else {
              for (j = change_point[i];j < change_point[i + 1];j++) {
                for (k = 0;k < nb_sequence;k++) {
                  diff = real_sequence[k][variable][j] - (intercept + slope * seq_index_parameter[j]);
                  square_sum += diff * diff;
                }
              }
            }

            if ((global_variance) || (model_type == INTERCEPT_SLOPE_CHANGE)) {
              global_square_sum += square_sum;
            }

            variance[0][i] = square_sum / (nb_sequence * (change_point[i + 1] - change_point[i]) - 2);

            if ((model_type == LINEAR_MODEL_CHANGE) && (likelihood != D_INF)) {
              if (square_sum > nb_sequence * (change_point[i + 1] - change_point[i]) * ROUNDOFF_ERROR) {
                likelihood -= ((double)(nb_sequence * (change_point[i + 1] - change_point[i])) / 2.) * (log(square_sum /
                                (nb_sequence * (change_point[i + 1] - change_point[i]))) + log(2 * M_PI) + 1);
              }
              else {
                likelihood = D_INF;
              }
            }

            if (slope_standard_deviation) {
              square_sum /= (nb_sequence * (change_point[i + 1] - change_point[i]) - 2);
              slope_standard_deviation[0][i] = sqrt(square_sum / index_parameter_variance);
            }
          }

          else {
            variance[0][i] = 0.;
          }
        }
      }
    }

    if (global_variance) {
      if (model_type == INTERCEPT_SLOPE_CHANGE) {
        if (index != I_DEFAULT) {
          global_variance[variable] = global_square_sum / (length[index] - 2 * nb_segment);
        }
        else {
          global_variance[variable] = global_square_sum / (nb_sequence * length[0] - 2 * nb_segment);
        }

        if (index != I_DEFAULT) {
          if (global_square_sum > length[index] * ROUNDOFF_ERROR) {
            likelihood = -((double)length[index] / 2.) * (log(global_square_sum /
                            length[index]) + log(2 * M_PI) + 1);
          }
          else {
            likelihood = D_INF;
          }
        }

        else {
          if (global_square_sum > nb_sequence * length[0] * ROUNDOFF_ERROR) {
            likelihood = -((double)(nb_sequence * length[0]) / 2.) * (log(global_square_sum /
                            (nb_sequence * length[0])) + log(2 * M_PI) + 1);
          }
          else {
            likelihood = D_INF;
          }
        }
      }

      // computation of mean squared error

      else {
        if (index != I_DEFAULT) {
          global_variance[variable] = global_square_sum / length[index];
        }
        else {
          global_variance[variable] = global_square_sum / (nb_sequence * length[0]);
        }
      }
    }
  }

  else if ((model_type == AUTOREGRESSIVE_MODEL_CHANGE) || (model_type == STATIONARY_AUTOREGRESSIVE_MODEL_CHANGE)) {
    if (variance) {
      likelihood = 0.;
    }
    else {
      likelihood = D_INF;
    }

    if (global_variance) {
      mean_squared_error_1 = 0.;
      global_square_sum = 0.;

#     ifdef MESSAGE
      mean_squared_error = 0.;
#     endif

    }

    if ((index != I_DEFAULT) || (!common_contrast)) {
      for (i = 0;i < nb_sequence;i++) {
        if ((index == I_DEFAULT) || (index == i)) {
          if (model_type == STATIONARY_AUTOREGRESSIVE_MODEL_CHANGE) {
            mean = 0.;

            if (type[variable] != REAL_VALUE) {
              for (j = 0;j < length[i];j++) {
                mean += int_sequence[i][variable][j];
              }
            }
            else {
              for (j = 0;j < length[i];j++) {
                mean += real_sequence[i][variable][j];
              }
            }
            mean /= length[i];
          }

          for (j = 0;j < nb_segment;j++) {
            if (model_type == AUTOREGRESSIVE_MODEL_CHANGE) {
              mean = 0.;

              if (type[variable] != REAL_VALUE) {
                for (k = change_point[j];k < change_point[j + 1];k++) {
                  mean += int_sequence[i][variable][k];
                }
              }
              else {
                for (k = change_point[j];k < change_point[j + 1];k++) {
                  mean += real_sequence[i][variable][k];
                }
              }
              mean /= (change_point[j + 1] - change_point[j]);
            }

            if (imean) {
              imean[i][j] = mean;
            }

            square_sum = 0.;
            shifted_square_sum = 0.;
            autocovariance = 0.;
            if (type[variable] != REAL_VALUE) {
              for (k = change_point[j] + 1;k < change_point[j + 1];k++) {
                diff = int_sequence[i][variable][k] - mean;
                shifted_diff = int_sequence[i][variable][k - 1] - mean;
                square_sum += diff * diff;
                shifted_square_sum += shifted_diff * shifted_diff;
                autocovariance += diff * shifted_diff;
              }
            }
            else {
              for (k = change_point[j] + 1;k < change_point[j + 1];k++) {
                diff = real_sequence[i][variable][k] - mean;
                shifted_diff = real_sequence[i][variable][k - 1] - mean;
                square_sum += diff * diff;
                shifted_square_sum += shifted_diff * shifted_diff;
                autocovariance += diff * shifted_diff;
              }
            }

            if (shifted_square_sum > 0.) {
              autoregressive_coeff = autocovariance / shifted_square_sum;
              if (autoregressive_coeff < -1.) {
                autoregressive_coeff = -1.;
              }
              else if (autoregressive_coeff > 1.) {
                autoregressive_coeff = 1.;
              }

              if (iautoregressive_coeff) {
                iautoregressive_coeff[i][j] = autoregressive_coeff;
              }

              if (piecewise_function) {
                piecewise_function[i][change_point[j]] = mean;

                if (type[variable] != REAL_VALUE) {
                  for (k = change_point[j] + 1;k < change_point[j + 1];k++) {
                    piecewise_function[i][k] = mean + autoregressive_coeff * (int_sequence[i][variable][k - 1] - mean);
                  }
                }
                else {
                  for (k = change_point[j] + 1;k < change_point[j + 1];k++) {
                    piecewise_function[i][k] = mean + autoregressive_coeff * (real_sequence[i][variable][k - 1] - mean);
                  }
                }
              }

              if (global_variance) {
                if (type[variable] != REAL_VALUE) {
                  diff = int_sequence[i][variable][change_point[j]] - mean;
                  mean_squared_error_1 += diff * diff;

#                 ifdef MESSAGE
                  mean_squared_error += diff * diff;
                  for (k = change_point[j] + 1;k < change_point[j + 1];k++) {
                    diff = int_sequence[i][variable][k] - (mean + autoregressive_coeff * (int_sequence[i][variable][k - 1] - mean));
                    mean_squared_error += diff * diff;
                  }
#                 endif

                }

                else {
                  diff = real_sequence[i][variable][change_point[j]] - mean;
                  mean_squared_error_1 += diff * diff;

#                 ifdef MESSAGE
                  mean_squared_error += diff * diff;
                  for (k = change_point[j] + 1;k < change_point[j + 1];k++) {
                    diff = real_sequence[i][variable][k] - (mean + autoregressive_coeff * (real_sequence[i][variable][k - 1] - mean));
                    mean_squared_error += diff * diff;
                  }
#                 endif

                }
              }

              if ((variance) || (global_variance)) {
                residual_square_sum = square_sum - autocovariance * autocovariance / shifted_square_sum;
                if (global_variance) {
                  global_square_sum += residual_square_sum;
                }
//                variance[i][j] = residual_square_sum / (change_point[j + 1] - change_point[j] - 3);
                variance[i][j] = residual_square_sum / (change_point[j + 1] - change_point[j] - 2);

                if (determination_coeff) {
                  determination_coeff[i][j] = 1.;
                  if (square_sum > 0.) {
                    determination_coeff[i][j] -= residual_square_sum / square_sum;
                  }
                }
 
                if (likelihood != D_INF) {
                  if (residual_square_sum > (change_point[j + 1] - change_point[j] - 1) * ROUNDOFF_ERROR) {
                    likelihood -= ((double)(change_point[j + 1] - change_point[j] - 1) / 2.) *(log(residual_square_sum /
                                    (change_point[j + 1] - change_point[j] - 1)) + log(2 * M_PI) + 1);
                  }
                  else {
                    likelihood = D_INF;
                  }
                }
              }
            }

            else {
              if (iautoregressive_coeff) {
                iautoregressive_coeff[i][j] = 0.;
              }

              if (piecewise_function) {
                for (k = change_point[j];k < change_point[j + 1];k++) {
                  piecewise_function[i][k] = mean;
                }
              }

              if (variance) {
                variance[i][j] = D_DEFAULT;
                if (determination_coeff) {
                  determination_coeff[i][j] = 0.;
                }
                likelihood = D_INF;
              }
            }
          }
        }
      }
    }

    else {
      if (model_type == STATIONARY_AUTOREGRESSIVE_MODEL_CHANGE) {
        mean = 0.;

        if (type[variable] != REAL_VALUE) {
          for (i = 0;i < nb_sequence;i++) {
            for (j = 0;j < length[0];j++) {
              mean += int_sequence[i][variable][j];
            }
          }
        }
        else {
          for (i = 0;i < nb_sequence;i++) {
            for (j = 0;j < length[0];j++) {
              mean += real_sequence[i][variable][j];
            }
          }
        }
        mean /= (nb_sequence * length[0]);
      }

      for (i = 0;i < nb_segment;i++) {
        if (model_type == AUTOREGRESSIVE_MODEL_CHANGE) {
          mean = 0.;

          if (type[variable] != REAL_VALUE) {
            for (j = 0;j < nb_sequence;j++) {
              for (k = change_point[i];k < change_point[i + 1];k++) {
                mean += int_sequence[j][variable][k];
              }
            }
          }
          else {
            for (j = 0;j < nb_sequence;j++) {
              for (k = change_point[i];k < change_point[i + 1];k++) {
                mean += real_sequence[j][variable][k];
              }
            }
          }
          mean /= (nb_sequence * (change_point[i + 1] - change_point[i]));
        }

        if (imean) {
          imean[0][i] = mean;
        }

        square_sum = 0.;
        shifted_square_sum = 0.;
        autocovariance = 0.;
        if (type[variable] != REAL_VALUE) {
          for (j = 0;j < nb_sequence;j++) {
            for (k = change_point[i] + 1;k < change_point[i + 1];k++) {
              diff = int_sequence[j][variable][k] - mean;
              shifted_diff = int_sequence[j][variable][k - 1] - mean;
              square_sum += diff * diff;
              shifted_square_sum += shifted_diff * shifted_diff;
              autocovariance += diff * shifted_diff;
            }
          }
        }
        else {
          for (j = 0;j < nb_sequence;j++) {
            for (k = change_point[i] + 1;k < change_point[i + 1];k++) {
              diff = real_sequence[j][variable][k] - mean;
              shifted_diff = real_sequence[j][variable][k - 1] - mean;
              square_sum += diff * diff;
              shifted_square_sum += shifted_diff * shifted_diff;
              autocovariance += diff * shifted_diff;
            }
          }
        }

        if (shifted_square_sum > 0.) {
          autoregressive_coeff = autocovariance / shifted_square_sum;
          if (autoregressive_coeff < -1.) {
            autoregressive_coeff = -1.;
          }
          else if (autoregressive_coeff > 1.) {
            autoregressive_coeff = 1.;
          }

          if (iautoregressive_coeff) {
            iautoregressive_coeff[0][i] = autoregressive_coeff;
          }

          if (piecewise_function) {
            for (j = 0;j < nb_sequence;j++) {
              piecewise_function[j][change_point[i]] = mean;

              if (type[variable] != REAL_VALUE) {
                for (k = change_point[i] + 1;k < change_point[i + 1];k++) {
                  piecewise_function[j][k] = mean + autoregressive_coeff * (int_sequence[j][variable][k - 1] - mean);
                }
              }
              else {
                for (k = change_point[i] + 1;k < change_point[i + 1];k++) {
                  piecewise_function[j][k] = mean + autoregressive_coeff * (real_sequence[j][variable][k - 1] - mean);
                }
              }
            }
          }

          if (global_variance) {
            for (j = 0;j < nb_sequence;j++) {
              if (type[variable] != REAL_VALUE) {
                diff = int_sequence[j][variable][change_point[i]] - mean;
                mean_squared_error_1 += diff * diff;

#               ifdef MESSAGE
                mean_squared_error += diff * diff;
                for (k = change_point[i] + 1;k < change_point[i + 1];k++) {
                  diff = int_sequence[j][variable][k] - (mean + autoregressive_coeff * (int_sequence[j][variable][k - 1] - mean));
                  mean_squared_error += diff * diff;
                }
#               endif

              }
              else {
                diff = real_sequence[j][variable][change_point[i]] - mean;
                mean_squared_error_1 += diff * diff;

#               ifdef MESSAGE
                mean_squared_error += diff * diff;
                for (k = change_point[i] + 1;k < change_point[i + 1];k++) {
                  diff= real_sequence[j][variable][k] - (mean + autoregressive_coeff * (real_sequence[j][variable][k - 1] - mean));
                  mean_squared_error += diff * diff;
                }
#               endif

              }
            }
          }

          if ((variance) || (global_variance)) {
            residual_square_sum = square_sum - autocovariance * autocovariance / shifted_square_sum;
            if (global_variance) {
              global_square_sum += residual_square_sum;
            }

//            variance[0][i] = residual_square_sum / (nb_sequence * (change_point[i + 1] - change_point[i] - 1) - 2);
            variance[0][i] = residual_square_sum / (nb_sequence * (change_point[i + 1] - change_point[i] - 1) - 1);

            if (determination_coeff) {
              determination_coeff[0][i] = 1.;
              if (square_sum > 0.) {
                determination_coeff[0][i] -= residual_square_sum / square_sum;
              }
            }

            if (likelihood != D_INF) {
              if (residual_square_sum > nb_sequence * (change_point[i + 1] - change_point[i] - 1) * ROUNDOFF_ERROR) {
                likelihood -= ((double)(nb_sequence * (change_point[i + 1] - change_point[i] - 1)) / 2.) * (log(residual_square_sum /
                                (nb_sequence * (change_point[i + 1] - change_point[i] - 1))) + log(2 * M_PI) + 1);
              }
              else {
                likelihood = D_INF;
              }
            }
          }
        }

        else {
          if (iautoregressive_coeff) {
            iautoregressive_coeff[0][i] = 0.;
          }

          if (piecewise_function) {
            for (j = 0;j < nb_sequence;j++) {
              for (k = change_point[i];k < change_point[i + 1];k++) {
                piecewise_function[j][k] = mean;
              }
            }
          }

          if (variance) {
            variance[0][i] = D_DEFAULT;
            if (determination_coeff) {
              determination_coeff[0][i] = 0.;
            }
            likelihood = D_INF;
          }
        }
      }
    }

    // computation of mean squared error

    if (global_variance) {
      if (index != I_DEFAULT) {
//        global_variance[variable] = global_square_sum / (length[index] - 3 * nb_segment);
//        global_variance[variable] = global_square_sum / (length[index] - 2 * nb_segment);
        global_variance[variable] = (mean_squared_error_1 + global_square_sum) / length[index];
       
      }

      else {
//        global_variance[variable] = global_square_sum / (nb_sequence * (length[0] - nb_segment) - 2 * nb_segment);
//        global_variance[variable] = global_square_sum / (nb_sequence * (length[0] - nb_segment) - nb_segment);
        global_variance[variable] = (mean_squared_error_1 + global_square_sum) / (nb_sequence * length[0]);
      }

#     ifdef MESSAGE
      if ((model_type == AUTOREGRESSIVE_MODEL_CHANGE) || (model_type == STATIONARY_AUTOREGRESSIVE_MODEL_CHANGE)) {
        if (index != I_DEFAULT) {
          mean_squared_error /= length[index];
        }
        else {
          mean_squared_error /= (nb_sequence * length[0]);
        }

        if ((global_variance[variable] < mean_squared_error - DOUBLE_ERROR) ||
            (global_variance[variable] > mean_squared_error + DOUBLE_ERROR)) {
          cout << "\nERROR " << SEQ_label[SEQL_ROOT_MEAN_SQUARE_ERROR] << ": " << sqrt(global_variance[variable]) << " | "
               << sqrtl(mean_squared_error) << endl;
        }
      }
#     endif

    }
  }

  return likelihood;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of piecewise linear functions.
 *
 *  \param[in,out] os                       stream,
 *  \param[in]     index                    sequence index,
 *  \param[in]     variable                 variable index,
 *  \param[in]     nb_segment               number of segments,
 *  \param[in]     model_type               segment model type,
 *  \param[in]     common_contrast          flag contrast functions common to the individuals,
 *  \param[in]     change_point             change points,
 *  \param[in]     seq_index_parameter      index parameters,
 *  \param[in]     mean                     segment means,
 *  \param[in]     variance                 segment variances or residual variances,
 *  \param[in]     intercept                segment intercepts,
 *  \param[in]     slope                    segment slopes,
 *  \param[in]     autoregressive_coeff     segment autoregressive coefficient,
 *  \param[in]     correlation              segment correlation coefficients (for linear models),
 *  \param[in]     slope_standard_deviation segment slope standard deviations (for linear models).
 *  \param[in]     index_parameter_mean     segment index parameter mean (for linear models),
 *  \param[in]     index_parameter_variance segment index parameter variance (for linear models),
 *  \param[in]     determination_coeff      coefficient of determination (for autoregressive models).
 */
/*--------------------------------------------------------------*/

ostream& Sequences::piecewise_linear_function_ascii_print(ostream &os , int index , int variable , int nb_segment ,
                                                          segment_model model_type , bool common_contrast , int *change_point ,
                                                          int *seq_index_parameter , double **mean , double **variance ,
                                                          double **intercept , double **slope , double **autoregressive_coeff ,
                                                          double **correlation , double **slope_standard_deviation ,
                                                          double **index_parameter_mean , long double **index_parameter_variance ,
                                                          double **determination_coeff) const

{
  int i , j;
  double diff , buff;
  Test *test;


  if ((model_type == POISSON_CHANGE) || (model_type == BAYESIAN_POISSON_CHANGE)) {
    if (nb_variable > 2) {
      os << STAT_label[STATL_VARIABLE] << " " << variable << "   ";
    }
    os << SEQ_label[SEQL_SEGMENT] << " " << STAT_label[STATL_MEAN] << ", "
       << STAT_label[STATL_VARIANCE] << ": ";

    if ((index != I_DEFAULT) || (!common_contrast)) {
      for (i = 0;i < nb_sequence;i++) {
        if ((index == I_DEFAULT) || (index == i)) {
          for (j = 0;j < nb_segment;j++) {
            os << mean[i][j] << " " << variance[i][j];
            if (j < nb_segment - 1) {
              os << " | ";
            }
            else if ((index == I_DEFAULT) && (i < nb_sequence - 1)) {
              os << endl;
            }
          }
        }
      }
    }

    else {
      for (i = 0;i < nb_segment;i++) {
        os << mean[0][i] << " " << variance[0][i];
        if (variance[1][i] > 0.) {
          os << " (" << 100 * variance[2][i] / variance[1][i] << "%, " << 100 * variance[3][i] / variance[1][i] << "%)";
        }
        if (i < nb_segment - 1) {
          os << " | ";
        }
      }
    }
    os << endl;
  }

  else if ((model_type == NEGATIVE_BINOMIAL_0_CHANGE) || (model_type == NEGATIVE_BINOMIAL_1_CHANGE) ||
           (model_type == GAUSSIAN_CHANGE) || (model_type == MEAN_CHANGE) ||
           (model_type == VARIANCE_CHANGE) || (model_type == BAYESIAN_GAUSSIAN_CHANGE)) {
    if (nb_variable > 2) {
      os << STAT_label[STATL_VARIABLE] << " " << variable << "   ";
    }

    if ((index != I_DEFAULT) || (!common_contrast)) {
      os << SEQ_label[SEQL_SEGMENT] << " " << STAT_label[STATL_MEAN] << ", "
         << STAT_label[STATL_STANDARD_DEVIATION] << ": ";

      for (i = 0;i < nb_sequence;i++) {
        if ((index == I_DEFAULT) || (index == i)) {
          for (j = 0;j < nb_segment;j++) {
            os << mean[i][j] << " " << sqrt(variance[i][j]);
            if (j < nb_segment - 1) {
              os << " | ";
            }
            else if ((index == I_DEFAULT) && (i < nb_sequence - 1)) {
              os << endl;
            }
          }
        }
      }
    }

    else {
      os << SEQ_label[SEQL_SEGMENT] << " " << STAT_label[STATL_MEAN] << ", "
         << STAT_label[STATL_STANDARD_DEVIATION];
      if (model_type != VARIANCE_CHANGE) {
        os << ", " << STAT_label[STATL_VARIANCE];
      }
      os << ": ";

      for (i = 0;i < nb_segment;i++) {
        os << mean[0][i] << " " << sqrt(variance[0][i]);
        if ((model_type != VARIANCE_CHANGE) && (variance[1][i] > 0.)) {
          os << " " << variance[0][i] << " (" << 100 * variance[2][i] / variance[1][i] << "%, "
             << 100 * variance[3][i] / variance[1][i] << "%)";
        }
        if (i < nb_segment - 1) {
          os << " | ";
        }
      }
    }
    os << endl;
  }

  else if ((model_type == LINEAR_MODEL_CHANGE) || (model_type == INTERCEPT_SLOPE_CHANGE)) {
    if (nb_variable > 2) {
      os << STAT_label[STATL_VARIABLE] << " " << variable << "   ";
    }

    if ((index == I_DEFAULT) && (!common_contrast)) {
      os << SEQ_label[SEQL_SEGMENT] << " " << STAT_label[STATL_INTERCEPT] << ", "
         << STAT_label[STATL_SLOPE] << ", " << STAT_label[STATL_RESIDUAL] << " "
         << STAT_label[STATL_STANDARD_DEVIATION] << ": ";
      for (i = 0;i < nb_sequence;i++) {
        for (j = 0;j < nb_segment;j++) {
          os << intercept[i][j] << " " << slope[i][j] << " " << sqrt(variance[i][j]);
          if (j < nb_segment - 1) {
            os << " | ";
          }
        }
        os << endl;
      }
    }

    else {
      if ((correlation) && (slope_standard_deviation)) {
        os << SEQ_label[SEQL_SEGMENT] << " " << STAT_label[STATL_INTERCEPT] << ", "
           << STAT_label[STATL_SLOPE] << " (" << SEQ_label[SEQL_CONFIDENCE_INTERVAL] << "), "
           << STAT_label[STATL_CORRELATION_COEFF] << " (" << STAT_label[STATL_LIMIT_CORRELATION_COEFF] << "), "
           << STAT_label[STATL_RESIDUAL] << " " << STAT_label[STATL_STANDARD_DEVIATION] << ", "
           << SEQ_label[SEQL_CHANGE_POINT_AMPLITUDE] << " (" << SEQ_label[SEQL_CONFIDENCE_INTERVALS] << ")" << endl;

        if (index != I_DEFAULT) {
          test = new Test(STUDENT , false , change_point[1] - change_point[0] - 2 , I_DEFAULT , D_DEFAULT);
          test->critical_probability = ref_critical_probability[0];
          test->t_value_computation();

          for (i = 0;i < nb_segment;i++) {
            os << intercept[index][i] << ", " << slope[index][i];
            if (slope_standard_deviation[index][i] > 0.) {
              os << " (" << slope[index][i] - test->value * slope_standard_deviation[index][i] << ", "
                 << slope[index][i] + test->value * slope_standard_deviation[index][i] << ")";
//               << " (slope_standard_deviation: " << slope_standard_deviation[index][i] << ")";
            }
            os << ", " << correlation[index][i] << " (-/+"
               << test->value / sqrt(test->value * test->value + change_point[i + 1] - change_point[i] - 2)
               << "), " << sqrt(variance[index][i]);

            if (i < nb_segment - 1) {
              os << ", " << intercept[index][i + 1] + slope[index][i + 1] * seq_index_parameter[change_point[i + 1]] -
                           (intercept[index][i] + slope[index][i] * seq_index_parameter[change_point[i + 1]]) << " (";

              diff = seq_index_parameter[change_point[i + 1]] - index_parameter_mean[index][i];
              buff = test->value * sqrt(variance[index][i] * (1. / (double)(change_point[i + 1] - change_point[i]) +
                     diff * diff / index_parameter_variance[index][i]));

//              os << test->value << ", ";
              os << MAX(intercept[index][i] + slope[index][i] * seq_index_parameter[change_point[i + 1]] - buff , 0) << ", "
                 << intercept[index][i] + slope[index][i] * seq_index_parameter[change_point[i + 1]] + buff << " | ";

              delete test;

              test = new Test(STUDENT , false , change_point[i + 2] - change_point[i + 1] - 2 , I_DEFAULT , D_DEFAULT);
              test->critical_probability = ref_critical_probability[0];
              test->t_value_computation();

              diff = seq_index_parameter[change_point[i + 1]] - index_parameter_mean[index][i + 1];
              buff = test->value * sqrt(variance[index][i + 1] * (1. / (double)(change_point[i + 2] - change_point[i + 1]) +
                     diff * diff / index_parameter_variance[index][i + 1]));
              os << MAX(intercept[index][i + 1] + slope[index][i + 1] * seq_index_parameter[change_point[i + 1]] - buff , 0) << ", "
                 << intercept[index][i + 1] + slope[index][i + 1] * seq_index_parameter[change_point[i + 1]] + buff << ")";
            }
            os << endl;
          }

          delete test;

          os << SEQ_label[SEQL_PIECEWISE_LINEAR_FUNCTION] << ": ";
          for (i = 0;i < nb_segment;i++) {
            os << intercept[index][i] + slope[index][i] * seq_index_parameter[change_point[i]] << " -> ";
            if (i < nb_segment - 1) {
              os << intercept[index][i] + slope[index][i] * seq_index_parameter[change_point[i + 1]] << " | ";
            }
            else {
              os << intercept[index][i] + slope[index][i] * seq_index_parameter[change_point[i + 1] - 1] << endl;
            }
          }
        }

        else if (common_contrast) {
          test = new Test(STUDENT , false , nb_sequence * (change_point[1] - change_point[0]) - 2 , I_DEFAULT , D_DEFAULT);
          test->critical_probability = ref_critical_probability[0];
          test->t_value_computation();

          for (i = 0;i < nb_segment;i++) {
            os << intercept[0][i] << ", " << slope[0][i];
            if (slope_standard_deviation[0][i] > 0.) {
              os << " (" << slope[0][i] - test->value * slope_standard_deviation[0][i] << ", "
                 << slope[0][i] + test->value * slope_standard_deviation[0][i] << ")";
//               << " (slope_standard_deviation: " << slope_standard_deviation[0][i] << ")";
            }
            os << ", " << correlation[0][i] << " (-/+"
               << test->value / sqrt(test->value * test->value + nb_sequence * (change_point[i + 1] - change_point[i]) - 2)
               << "), " << sqrt(variance[0][i]);

            if (i < nb_segment - 1) {
              os << ", " << intercept[0][i + 1] + slope[0][i + 1] * seq_index_parameter[change_point[i + 1]] -
                           (intercept[0][i] + slope[0][i] * seq_index_parameter[change_point[i + 1]]) << " (";

              diff = seq_index_parameter[change_point[i + 1]] - index_parameter_mean[0][i];
              buff = test->value * sqrt(variance[0][i] * (1. / (double)(nb_sequence * (change_point[i + 1] - change_point[i])) +
                     diff * diff / index_parameter_variance[0][i]));

//              os << test->value << ", ";
              os << MAX(intercept[0][i] + slope[0][i] * seq_index_parameter[change_point[i + 1]] - buff , 0) << ", "
                 << intercept[0][i] + slope[0][i] * seq_index_parameter[change_point[i + 1]] + buff << " | ";

              delete test;

              test = new Test(STUDENT , false , nb_sequence * (change_point[i + 2] - change_point[i + 1]) - 2 , I_DEFAULT , D_DEFAULT);
              test->critical_probability = ref_critical_probability[0];
              test->t_value_computation();

              diff = seq_index_parameter[change_point[i + 1]] - index_parameter_mean[0][i + 1];
              buff = test->value * sqrt(variance[0][i + 1] * (1. / (double)(nb_sequence * (change_point[i + 2] - change_point[i + 1])) +
                     diff * diff / index_parameter_variance[0][i + 1]));
              os << MAX(intercept[0][i + 1] + slope[0][i + 1] * seq_index_parameter[change_point[i + 1]] - buff , 0) << ", "
                 << intercept[0][i + 1] + slope[0][i + 1] * seq_index_parameter[change_point[i + 1]] + buff << ")";
            }
            os << endl;
          }

          delete test;

          os << SEQ_label[SEQL_PIECEWISE_LINEAR_FUNCTION] << ": ";
          for (i = 0;i < nb_segment;i++) {
            os << intercept[0][i] + slope[0][i] * seq_index_parameter[change_point[i]] << " -> ";
            if (i < nb_segment - 1) {
              os << intercept[0][i] + slope[0][i] * seq_index_parameter[change_point[i + 1]] << " | ";
            }
            else {
              os << intercept[0][i] + slope[0][i] * seq_index_parameter[change_point[i + 1] - 1] << endl;
            }
          }
        }
      }

      else {
        os << SEQ_label[SEQL_PIECEWISE_LINEAR_FUNCTION] << ", " 
           << SEQ_label[SEQL_SEGMENT] << " " << STAT_label[STATL_INTERCEPT] << ", "
           << STAT_label[STATL_SLOPE] << ", " << STAT_label[STATL_RESIDUAL] << " "
           << STAT_label[STATL_STANDARD_DEVIATION] << ": ";

        if (index != I_DEFAULT) {
          for (i = 0;i < nb_segment;i++) {
            os << intercept[index][i] + slope[index][i] * seq_index_parameter[change_point[i]] << " -> ";
            if (i < nb_segment - 1) {
              os << intercept[index][i] + slope[index][i] * seq_index_parameter[change_point[i + 1]] << " | ";
            }
            else {
              os << intercept[index][i] + slope[index][i] * seq_index_parameter[change_point[i + 1] - 1];
            }
          }
          os << " || ";
          for (i = 0;i < nb_segment;i++) {
            os << intercept[index][i] << " " << slope[index][i] << " " << sqrt(variance[index][i]);
            if (i < nb_segment - 1) {
              os << " | ";
            }
          }
        }

        else if (common_contrast) {
          for (i = 0;i < nb_segment;i++) {
            os << intercept[0][i] + slope[0][i] * seq_index_parameter[change_point[i]] << " -> ";
            if (i < nb_segment - 1) {
              os << intercept[0][i] + slope[0][i] * seq_index_parameter[change_point[i + 1]] << " | ";
            }
            else {
              os << intercept[0][i] + slope[0][i] * seq_index_parameter[change_point[i + 1] - 1];
            }
          }
          os << " || ";
          for (i = 0;i < nb_segment;i++) {
            os << intercept[0][i] << " " << slope[0][i] << " " << sqrt(variance[0][i]);
            if (i < nb_segment - 1) {
              os << " | ";
            }
          }
        }
        os << endl;
      }
    }
  }

  else if ((model_type == AUTOREGRESSIVE_MODEL_CHANGE) || (model_type == STATIONARY_AUTOREGRESSIVE_MODEL_CHANGE)) {
    normal dist;
    double standard_normal_value = quantile(complement(dist , 0.025)) , standard_error;

    if (nb_variable > 2) {
      os << STAT_label[STATL_VARIABLE] << " " << variable << "   ";
    }
    os << SEQ_label[SEQL_SEGMENT] << " " << STAT_label[STATL_MEAN] << ", "
       << SEQ_label[SEQL_AUTOREGRESSIVE_COEFF] << " ("  << SEQ_label[SEQL_CONFIDENCE_INTERVAL] << " | " 
       << STAT_label[STATL_NULL_AUTOREGRESSIVE_COEFF_95_CONFIDENCE_LIMIT] << "), "
       << STAT_label[STATL_RESIDUAL] << " " << STAT_label[STATL_STANDARD_DEVIATION];
    if (determination_coeff) {
      os << ", " << STAT_label[STATL_STANDARD_DEVIATION];
      os << ", " << STAT_label[STATL_DETERMINATION_COEFF];
    }
    os << endl;

    if ((index != I_DEFAULT) || (!common_contrast)) {
      for (i = 0;i < nb_sequence;i++) {
        if ((index == I_DEFAULT) || (index == i)) {
/*          if (index == I_DEFAULT) {
            os << endl;
          }
          else {
            os << ": ";
          } */

          for (j = 0;j < nb_segment;j++) {
            standard_error = standard_normal_value * sqrt((1. - autoregressive_coeff[i][j] * autoregressive_coeff[i][j]) /
                                                          (change_point[j + 1] - change_point[j]));
            os << mean[i][j] << " " << autoregressive_coeff[i][j] << " ("
               << MAX(autoregressive_coeff[i][j] - standard_error , -1.) << ", " << MIN(autoregressive_coeff[i][j] + standard_error , 1.)
               << " | -/+" << standard_normal_value / sqrt((double)(change_point[j + 1] - change_point[j])) << ") "
               << sqrt(variance[i][j]);
            if (determination_coeff) {
              os << " " << sqrt(variance[i][j] / (1. - determination_coeff[i][j]));
              os << " " << determination_coeff[i][j];
            }
            os << endl;
/*            if (j < nb_segment - 1) {
              os << " || ";
            } */
          }
        }
      }
    }

    else {
//      os << ": ";
      for (i = 0;i < nb_segment;i++) {
        standard_error = standard_normal_value * sqrt((1. - autoregressive_coeff[0][i] * autoregressive_coeff[0][i]) /
                                                      (nb_sequence * (change_point[i + 1] - change_point[i])));
        os << mean[0][i] << " " << autoregressive_coeff[0][i] << " ("
           << MAX(autoregressive_coeff[0][i] - standard_error , -1.) << ", " << MIN(autoregressive_coeff[0][i] + standard_error , 1.)
           << " | -/+" << standard_normal_value / sqrt((double)nb_sequence * (change_point[i + 1] - change_point[i])) << ") "
           << sqrt(variance[0][i]);
        if (determination_coeff) {
          os << " " << sqrt(variance[0][i] / (1. - determination_coeff[0][i]));
          os << "  " << determination_coeff[0][i];
        }
        os << endl;
/*        if (i < nb_segment - 1) {
          os << " || ";
        } */
      }
    }
//    os << endl;
  }

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of piecewise linear functions at the spreadsheet format.
 *
 *  \param[in,out] os                       stream,
 *  \param[in]     index                    sequence index,
 *  \param[in]     variable                 variable index,
 *  \param[in]     nb_segment               number of segments,
 *  \param[in]     model_type               segment model type,
 *  \param[in]     common_contrast          flag contrast functions common to the individuals,
 *  \param[in]     change_point             change points,
 *  \param[in]     seq_index_parameter      index parameters,
 *  \param[in]     mean                     segment means,
 *  \param[in]     variance                 segment variances or residual variances,
 *  \param[in]     intercept                segment intercepts,
 *  \param[in]     slope                    segment slopes,
 *  \param[in]     autoregressive_coeff     segment autoregressive coefficient,
 *  \param[in]     correlation              segment correlation coefficients (for linear models),
 *  \param[in]     slope_standard_deviation segment slope standard deviations (for linear models).
 *  \param[in]     index_parameter_mean     segment index parameter mean (for linear models),
 *  \param[in]     index_parameter_variance segment index parameter variance (for linear models),
 *  \param[in]     determination_coeff      coefficient of determination (for autoregressive models).
 */
/*--------------------------------------------------------------*/

ostream& Sequences::piecewise_linear_function_spreadsheet_print(ostream &os , int index , int variable , int nb_segment ,
                                                                segment_model model_type , bool common_contrast , int *change_point ,
                                                                int *seq_index_parameter , double **mean , double **variance ,
                                                                double **intercept , double **slope , double **autoregressive_coeff ,
                                                                double **correlation , double **slope_standard_deviation ,
                                                                double **index_parameter_mean , long double **index_parameter_variance ,
                                                                double **determination_coeff) const

{
  int i , j;
  double diff , buff;
  Test *test;


  if ((model_type == POISSON_CHANGE) || (model_type == BAYESIAN_POISSON_CHANGE)) {
    if (nb_variable > 2) {
      os << STAT_label[STATL_VARIABLE] << "\t" << variable << "\t";
    }
    os << SEQ_label[SEQL_SEGMENT] << "\t" << STAT_label[STATL_MEAN] << "\t"
       << STAT_label[STATL_VARIANCE];

    if ((index != I_DEFAULT) || (!common_contrast)) {
      for (i = 0;i < nb_sequence;i++) {
        if ((index == I_DEFAULT) || (index == i)) {
          for (j = 0;j < nb_segment;j++) {
            os << "\t" << mean[i][j] << "\t" << variance[i][j] << "\t";
          }
          if ((index == I_DEFAULT) && (i < nb_sequence - 1)) {
            os << "\t";
          }
        }
      }
    }

    else {
      for (i = 0;i < nb_segment;i++) {
        os << "\t" << mean[0][i] << "\t" << variance[0][i];
        if (variance[1][i] > 0.) {
          os << "\t" << 100 * variance[2][i] / variance[1][i] << "%\t" << 100 * variance[3][i] / variance[1][i] << "%";
        }
        if (i < nb_segment - 1) {
          os << "\t\t";
        }
      }
    }
    os << endl;
  }

  else if ((model_type == NEGATIVE_BINOMIAL_0_CHANGE) || (model_type == NEGATIVE_BINOMIAL_1_CHANGE) ||
           (model_type == GAUSSIAN_CHANGE) || (model_type == MEAN_CHANGE) ||
           (model_type == VARIANCE_CHANGE) || (model_type == BAYESIAN_GAUSSIAN_CHANGE)) {
    if (nb_variable > 2) {
      os << STAT_label[STATL_VARIABLE] << "\t" << variable << "\t";
    }

    if ((index != I_DEFAULT) || (!common_contrast)) {
      os << SEQ_label[SEQL_SEGMENT] << "\t" << STAT_label[STATL_MEAN] << "\t"
         << STAT_label[STATL_STANDARD_DEVIATION];

      for (i = 0;i < nb_sequence;i++) {
        if ((index == I_DEFAULT) || (index == i)) {
          for (j = 0;j < nb_segment;j++) {
            os << "\t" << mean[i][j] << "\t" << sqrt(variance[i][j]) << "\t";
          }
          if ((index == I_DEFAULT) && (i < nb_sequence - 1)) {
            os << "\t";
          }
        }
      }
    }

    else {
      os << SEQ_label[SEQL_SEGMENT] << "\t" << STAT_label[STATL_MEAN] << "\t"
         << STAT_label[STATL_STANDARD_DEVIATION];
      if (model_type != VARIANCE_CHANGE) {
        os << "\t" << STAT_label[STATL_VARIANCE];
      }

      for (i = 0;i < nb_segment;i++) {
        os << "\t" << mean[0][i] << "\t" << sqrt(variance[0][i]);
        if ((model_type != VARIANCE_CHANGE) && (variance[1][i] > 0.)) {
          os << "\t" << variance[0][i] << "\t" << 100 * variance[2][i] / variance[1][i] << "%\t"
             << 100 * variance[3][i] / variance[1][i] << "%";
        }
        if (i < nb_segment - 1) {
          os << "\t\t";
        }
      }
    }
    os << endl;
  }

  else if ((model_type == LINEAR_MODEL_CHANGE) || (model_type == INTERCEPT_SLOPE_CHANGE)) {
    if (nb_variable > 2) {
      os << STAT_label[STATL_VARIABLE] << "\t" << variable << "\t";
    }

    if ((index == I_DEFAULT) && (!common_contrast)) {
      os << SEQ_label[SEQL_PIECEWISE_LINEAR_FUNCTION] << "\t"
         << SEQ_label[SEQL_SEGMENT] << " " << STAT_label[STATL_INTERCEPT] << "\t"
         << SEQ_label[SEQL_SEGMENT] << " " << STAT_label[STATL_SLOPE] << "\t"
         << SEQ_label[SEQL_SEGMENT] << " " << STAT_label[STATL_RESIDUAL] << " "
         << STAT_label[STATL_STANDARD_DEVIATION] << endl;
      for (i = 0;i < nb_sequence;i++) {
        for (j = 0;j < nb_segment;j++) {
          os << intercept[i][j] + slope[i][j] * seq_index_parameter[change_point[j]] << "->";
          if (j < nb_segment - 1) {
            os << intercept[i][j] + slope[i][j] * seq_index_parameter[change_point[j + 1]] << "\t";
          }
          else {
            os << intercept[i][j] + slope[i][j] * seq_index_parameter[change_point[j + 1] - 1];
          }
        }
        for (j = 0;j < nb_segment;j++) {
          os << "\t\t" << intercept[i][j] << "\t" << slope[i][j] << "\t" << sqrt(variance[i][j]);
        }
        os << endl;
      }
    }

    else {
      if ((correlation) && (slope_standard_deviation)) {
        os << SEQ_label[SEQL_SEGMENT] << " " << STAT_label[STATL_INTERCEPT] << "\t"
           << SEQ_label[SEQL_SEGMENT] << " " << STAT_label[STATL_SLOPE] << "\t"
           << SEQ_label[SEQL_CONFIDENCE_INTERVAL] << "\t\t"
           << SEQ_label[SEQL_SEGMENT] << " " << STAT_label[STATL_CORRELATION_COEFF] << "\t"
           << STAT_label[STATL_LIMIT_CORRELATION_COEFF] << "\t" << SEQ_label[SEQL_SEGMENT] << " "
           << STAT_label[STATL_RESIDUAL] << " " << STAT_label[STATL_STANDARD_DEVIATION] << "\t"
           << SEQ_label[SEQL_CHANGE_POINT_AMPLITUDE] << "\t" << SEQ_label[SEQL_CONFIDENCE_INTERVALS] << endl;

        if (index != I_DEFAULT) {
          test = new Test(STUDENT , false , change_point[1] - change_point[0] - 2 , I_DEFAULT , D_DEFAULT);
          test->critical_probability = ref_critical_probability[0];
          test->t_value_computation();

          for (i = 0;i < nb_segment;i++) {
            os << intercept[index][i] << "\t" << slope[index][i];
            if (slope_standard_deviation[index][i] > 0.) {
              os << "\t" << slope[index][i] - test->value * slope_standard_deviation[index][i]
                 << "\t" << slope[index][i] + test->value * slope_standard_deviation[index][i];
            }
            os << "\t" << correlation[index][i] << "\t-/+"
               << test->value / sqrt(test->value * test->value + change_point[i + 1] - change_point[i] - 2)
               << "\t" << sqrt(variance[index][i]);

            if (i < nb_segment - 1) {
              os << "\t" << intercept[index][i + 1] + slope[index][i + 1] * seq_index_parameter[change_point[i + 1]] -
                           (intercept[index][i] + slope[index][i] * seq_index_parameter[change_point[i + 1]]);

              diff = seq_index_parameter[change_point[i + 1]] - index_parameter_mean[index][i];
              buff = test->value * sqrt(variance[index][i] * (1. / (double)(change_point[i + 1] - change_point[i]) +
                     diff * diff / index_parameter_variance[index][i]));

//              os << "\t" << test->value;
              os << "\t" << MAX(intercept[index][i] + slope[index][i] * seq_index_parameter[change_point[i + 1]] - buff , 0)
                 << "\t" << intercept[index][i] + slope[index][i] * seq_index_parameter[change_point[i + 1]] + buff;

              delete test;

              test = new Test(STUDENT , false , change_point[i + 2] - change_point[i + 1] - 2 , I_DEFAULT , D_DEFAULT);
              test->critical_probability = ref_critical_probability[0];
              test->t_value_computation();

              diff = seq_index_parameter[change_point[i + 1]] - index_parameter_mean[index][i + 1];
              buff = test->value * sqrt(variance[index][i + 1] * (1. / (double)(change_point[i + 2] - change_point[i + 1]) +
                     diff * diff / index_parameter_variance[index][i + 1]));
              os << "\t" << MAX(intercept[index][i + 1] + slope[index][i + 1] * seq_index_parameter[change_point[i + 1]] - buff , 0)
                 << "\t" << intercept[index][i + 1] + slope[index][i + 1] * seq_index_parameter[change_point[i + 1]] + buff;
            }
            os << endl;
          }

          delete test;

          os << SEQ_label[SEQL_PIECEWISE_LINEAR_FUNCTION] << "\t";
          for (i = 0;i < nb_segment;i++) {
            os << intercept[index][i] + slope[index][i] * seq_index_parameter[change_point[i]] << " -> ";
            if (i < nb_segment - 1) {
              os << intercept[index][i] + slope[index][i] * seq_index_parameter[change_point[i + 1]] << "\t";
            }
            else {
              os << intercept[index][i] + slope[index][i] * seq_index_parameter[change_point[i + 1] - 1] << endl;
            }
          }
        }

        else if (common_contrast) {
          test = new Test(STUDENT , false , nb_sequence * (change_point[1] - change_point[0]) - 2 , I_DEFAULT , D_DEFAULT);
          test->critical_probability = ref_critical_probability[0];
          test->t_value_computation();

          for (i = 0;i < nb_segment;i++) {
            os << intercept[0][i] << "\t" << slope[0][i];
            if (slope_standard_deviation[0][i] > 0.) {
              os << "\t" << slope[0][i] - test->value * slope_standard_deviation[0][i]
                 << "\t" << slope[0][i] + test->value * slope_standard_deviation[0][i];
            }
            os << "\t" << correlation[0][i] << "\t-/+"
               << test->value / sqrt(test->value * test->value + nb_sequence * (change_point[i + 1] - change_point[i]) - 2)
               << "\t" << sqrt(variance[0][i]);

            if (i < nb_segment - 1) {
              os << "\t" << intercept[0][i + 1] + slope[0][i + 1] * seq_index_parameter[change_point[i + 1]] -
                           (intercept[0][i] + slope[0][i] * seq_index_parameter[change_point[i + 1]]);

              diff = seq_index_parameter[change_point[i + 1]] - index_parameter_mean[0][i];
              buff = test->value * sqrt(variance[0][i] * (1. / (double)(nb_sequence * (change_point[i + 1] - change_point[i])) +
                     diff * diff / index_parameter_variance[0][i]));

//              os << "\t" << test->value;
              os << "\t" << MAX(intercept[0][i] + slope[0][i] * seq_index_parameter[change_point[i + 1]] - buff , 0)
                 << "\t" << intercept[0][i] + slope[0][i] * seq_index_parameter[change_point[i + 1]] + buff;

              delete test;

              test = new Test(STUDENT , false , nb_sequence * (change_point[i + 2] - change_point[i + 1]) - 2 , I_DEFAULT , D_DEFAULT);
              test->critical_probability = ref_critical_probability[0];
              test->t_value_computation();

              diff = seq_index_parameter[change_point[i + 1]] - index_parameter_mean[0][i + 1];
              buff = test->value * sqrt(variance[0][i + 1] * (1. / (double)(nb_sequence * (change_point[i + 2] - change_point[i + 1])) +
                     diff * diff / index_parameter_variance[0][i + 1]));
              os << "\t" << MAX(intercept[0][i + 1] + slope[0][i + 1] * seq_index_parameter[change_point[i + 1]] - buff , 0)
                 << "\t" << intercept[0][i + 1] + slope[0][i + 1] * seq_index_parameter[change_point[i + 1]] + buff;
            }
            os << endl;
          }

          delete test;

          os << SEQ_label[SEQL_PIECEWISE_LINEAR_FUNCTION] << "\t";
          for (i = 0;i < nb_segment;i++) {
            os << intercept[0][i] + slope[0][i] * seq_index_parameter[change_point[i]] << "->";
            if (i < nb_segment - 1) {
              os << intercept[0][i] + slope[0][i] * seq_index_parameter[change_point[i + 1]] << "\t";
            }
            else {
              os << intercept[0][i] + slope[0][i] * seq_index_parameter[change_point[i + 1] - 1] << endl;
            }
          }
        }
      }

      else {
        os << SEQ_label[SEQL_PIECEWISE_LINEAR_FUNCTION] << "\t" 
           << SEQ_label[SEQL_SEGMENT] << " " << STAT_label[STATL_INTERCEPT] << "\t"
           << SEQ_label[SEQL_SEGMENT] << " " << STAT_label[STATL_SLOPE] << "\t"
           << SEQ_label[SEQL_SEGMENT] << " " << STAT_label[STATL_RESIDUAL] << " "
           << STAT_label[STATL_STANDARD_DEVIATION] << "\t";

        if (index != I_DEFAULT) {
          for (i = 0;i < nb_segment;i++) {
            os << intercept[index][i] + slope[index][i] * seq_index_parameter[change_point[i]] << " -> ";
            if (i < nb_segment - 1) {
              os << intercept[index][i] + slope[index][i] * seq_index_parameter[change_point[i + 1]] << "\t";
            }
            else {
              os << intercept[index][i] + slope[index][i] * seq_index_parameter[change_point[i + 1] - 1];
            }
          }
          for (i = 0;i < nb_segment;i++) {
            os << "\t\t" << intercept[index][i] << "\t" << slope[index][i] << "\t" << sqrt(variance[index][i]);
          }
        }

        else if (common_contrast) {
          for (i = 0;i < nb_segment;i++) {
            os << intercept[0][i] + slope[0][i] * seq_index_parameter[change_point[i]] << " -> ";
            if (i < nb_segment - 1) {
              os << intercept[0][i] + slope[0][i] * seq_index_parameter[change_point[i + 1]] << "\t";
            }
            else {
              os << intercept[0][i] + slope[0][i] * seq_index_parameter[change_point[i + 1] - 1];
            }
          }
          for (i = 0;i < nb_segment;i++) {
            os << "\t\t" << intercept[0][i] << "\t" << slope[0][i] << "\t" << sqrt(variance[0][i]);
          }
        }
        os << endl;
      }
    }
  }

  else if ((model_type == AUTOREGRESSIVE_MODEL_CHANGE) || (model_type == STATIONARY_AUTOREGRESSIVE_MODEL_CHANGE)) {
    normal dist;
    double standard_normal_value = quantile(complement(dist , 0.025)) , standard_error;

    if (nb_variable > 2) {
      os << STAT_label[STATL_VARIABLE] << "\t" << variable << "\t";
    }
    os << SEQ_label[SEQL_SEGMENT] << " " << STAT_label[STATL_MEAN] << "\t"
       << SEQ_label[SEQL_SEGMENT] << " " << SEQ_label[SEQL_AUTOREGRESSIVE_COEFF] << "\t"
       << SEQ_label[SEQL_CONFIDENCE_INTERVAL] << "\t"
       << STAT_label[STATL_NULL_AUTOREGRESSIVE_COEFF_95_CONFIDENCE_LIMIT] << "\t"
       << SEQ_label[SEQL_SEGMENT] << " " << STAT_label[STATL_RESIDUAL] << " "
       << STAT_label[STATL_STANDARD_DEVIATION];
    if (determination_coeff) {
      os << "\t" << SEQ_label[SEQL_SEGMENT] << " " << STAT_label[STATL_STANDARD_DEVIATION];
      os << "\t" << SEQ_label[SEQL_SEGMENT] << " " << STAT_label[STATL_DETERMINATION_COEFF];
    }
    os << endl;

    if ((index != I_DEFAULT) || (!common_contrast)) {
      for (i = 0;i < nb_sequence;i++) {
        if ((index == I_DEFAULT) || (index == i)) {
/*          if (index == I_DEFAULT) {
            os << endl;
          }
          else {
            os << "\t";
          } */

          for (j = 0;j < nb_segment;j++) {
            standard_error = standard_normal_value * sqrt((1. - autoregressive_coeff[i][j] * autoregressive_coeff[i][j]) /
                                                          (change_point[j + 1] - change_point[j]));
            os << mean[i][j] << "\t"  << autoregressive_coeff[i][j] << "\t"
               << MAX(autoregressive_coeff[i][j] - standard_error , -1.) << "\t" << MIN(autoregressive_coeff[i][j] + standard_error , 1.)
               << "\t-/+" << standard_normal_value / sqrt((double)(change_point[j + 1] - change_point[j])) << "\t"
               << sqrt(variance[i][j]);
            if (determination_coeff) {
              os << "\t" << sqrt(variance[i][j] / (1. - determination_coeff[i][j]));
              os << "\t" << determination_coeff[i][j];
            }
            os << endl;
/*            if (j < nb_segment - 1) {
              os << "\t\t";
            } */
          }
        }
      }
    }

    else {
//      os << "\t";
      for (i = 0;i < nb_segment;i++) {
        standard_error = standard_normal_value * sqrt((1. - autoregressive_coeff[0][i] * autoregressive_coeff[0][i]) /
                                                      (nb_sequence * (change_point[i + 1] - change_point[i])));
        os << mean[0][i] << "\t" << autoregressive_coeff[0][i] << "\t"
           << MAX(autoregressive_coeff[0][i] - standard_error , -1.) << "\t" << MIN(autoregressive_coeff[0][i] + standard_error , 1.)
           << "\t-/+" << standard_normal_value / sqrt((double)nb_sequence * (change_point[i + 1] - change_point[i])) << "\t"
           << sqrt(variance[0][i]);
        if (determination_coeff) {
          os << "\t" << sqrt(variance[0][i] / (1. - determination_coeff[0][i]));
          os << "\t" << determination_coeff[0][i];
        }
        os << endl;
/*        if (i < nb_segment - 1) {
          os << "\t\t";
        } */
      }
    }
//    os << endl;
  }

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of continuous piecewise linear functions.
 *
 *  \param[in,out] os                  stream,
 *  \param[in]     index               sequence index,
 *  \param[in]     variable            variable index,
 *  \param[in]     nb_segment          number of segments,
 *  \param[in]     model_type          segment model type,
 *  \param[in]     common_contrast     flag contrast functions common to the individuals,
 *  \param[in]     change_point        change points,
 *  \param[in]     seq_index_parameter index parameters,
 *  \param[in]     intercept           segment intercepts,
 *  \param[in]     slope               segment  slopes,
 *  \param[in]     corrected_intercept corrected segment intercepts
 *  \param[in]     corrected_slope     corrected segment slopes.
 */
/*--------------------------------------------------------------*/

double Sequences::continuous_piecewise_linear_function(ostream &os , int index , int variable , int nb_segment ,
                                                       segment_model model_type , bool common_contrast ,
                                                       int *change_point , int *seq_index_parameter ,
                                                       double *intercept , double *slope ,
                                                       double *corrected_intercept , double *corrected_slope) const

{
  int i , j , k;
  double likelihood , diff , residual_mean , global_variance , *predicted_value , *variance;
  long double square_sum , global_square_sum , residual_square_sum , global_residual_square_sum;


  predicted_value = new double[nb_segment + 1];

  predicted_value[0] = intercept[0] + slope[0] * seq_index_parameter[0];
  for (i = 1;i < nb_segment;i++) {
    predicted_value[i] = (fabs(slope[i - 1]) * (intercept[i - 1] + slope[i - 1] * seq_index_parameter[change_point[i]]) +
                          fabs(slope[i]) * (intercept[i] + slope[i] * seq_index_parameter[change_point[i]])) /
                         (fabs(slope[i - 1]) + fabs(slope[i]));
  }
  predicted_value[nb_segment] = intercept[nb_segment - 1] + slope[nb_segment - 1] *
                                seq_index_parameter[length[index == I_DEFAULT ? 0 : index] - 1];

  for (i = 0;i < nb_segment;i++) {
    corrected_slope[i] = (predicted_value[i + 1] - predicted_value[i]) /
                         (seq_index_parameter[change_point[i + 1]] - seq_index_parameter[change_point[i]]);
    corrected_intercept[i] = predicted_value[i] - corrected_slope[i] * seq_index_parameter[change_point[i]];
  }

  variance = new double[nb_segment];

  likelihood = 0.;

  if (model_type == INTERCEPT_SLOPE_CHANGE) {
    global_square_sum = 0.;
    global_residual_square_sum = 0.;
  }

  if (index != I_DEFAULT) {
    for (i = 0;i < nb_segment;i++) {
      residual_mean = 0.;
      residual_square_sum = 0.;
      square_sum = 0.;
  
      if (change_point[i + 1] - change_point[i] > 2) {
        if (type[variable] != REAL_VALUE) {
          for (j = change_point[i];j < change_point[i + 1];j++) {
            diff = int_sequence[index][variable][j] - (corrected_intercept[i] + corrected_slope[i] * seq_index_parameter[j]);
            residual_mean += diff;
            square_sum += diff * diff;
          }
          residual_mean /= (change_point[i + 1] - change_point[i]);

          for (j = change_point[i];j < change_point[i + 1];j++) {
            diff = int_sequence[index][variable][j] - (corrected_intercept[i] + corrected_slope[i] * seq_index_parameter[j]) - residual_mean;
            residual_square_sum += diff * diff;
          }
        }

        else {
          for (j = change_point[i];j < change_point[i + 1];j++) {
            diff = real_sequence[index][variable][j] - (corrected_intercept[i] + corrected_slope[i] * seq_index_parameter[j]);
            residual_mean += diff;
            square_sum += diff * diff;
          }
          residual_mean /= (change_point[i + 1] - change_point[i]);

          for (j = change_point[i];j < change_point[i + 1];j++) {
            diff = real_sequence[index][variable][j] - (corrected_intercept[i] + corrected_slope[i] * seq_index_parameter[j]) - residual_mean;
            residual_square_sum += diff * diff;
          }
        }

        if (model_type == INTERCEPT_SLOPE_CHANGE) {
          global_square_sum += square_sum;
          global_residual_square_sum += residual_square_sum;
        }

#       ifdef DEBUG
        cout << "\nTEST " << STAT_label[STATL_RESIDUAL] << " " << STAT_label[STATL_STANDARD_DEVIATION] << ": " 
             << sqrt(square_sum / (change_point[i + 1] - change_point[i] - 2)) << " "
             << sqrt(residual_square_sum / (change_point[i + 1] - change_point[i] - 2)) << " | " << residual_mean << endl;
#       endif

        variance[i] = residual_square_sum / (change_point[i + 1] - change_point[i] - 2);

        if ((model_type == LINEAR_MODEL_CHANGE) && (likelihood != D_INF)) {
          if (square_sum > (change_point[i + 1] - change_point[i]) * ROUNDOFF_ERROR) {
            likelihood -= ((double)(change_point[i + 1] - change_point[i]) / 2.) * (log(square_sum /
                            (change_point[i + 1] - change_point[i])) + log(2 * M_PI) + 1);
          }
          else {
            likelihood = D_INF;
          }
        }
      }

      else {
        variance[i] = 0.;
      }
    }

    if (model_type == INTERCEPT_SLOPE_CHANGE) {
      if (global_square_sum > length[index] * ROUNDOFF_ERROR) {
        likelihood -= ((double)length[index] / 2.) * (log(global_square_sum / length[index]) +
                       log(2 * M_PI) + 1);

#       ifdef DEBUG
        cout << "\nTEST " << STAT_label[STATL_RESIDUAL] << " " << STAT_label[STATL_STANDARD_DEVIATION] << ": " 
             << sqrt(global_square_sum / (length[index] - 2 * nb_segment)) << " "
             << sqrt(global_residual_square_sum / (length[index] - 2 * nb_segment)) << endl;
#       endif

      }
      else {
        likelihood = D_INF;
      }
    }
  }

  else if (common_contrast) {
    for (i = 0;i < nb_segment;i++) {
      residual_mean = 0.;
      residual_square_sum = 0.;
      square_sum = 0.;

      if (change_point[i + 1] - change_point[i] > 2) {
        if (type[variable] != REAL_VALUE) {
          for (j = change_point[i];j < change_point[i + 1];j++) {
            for (k = 0;k < nb_sequence;k++) {
              diff = int_sequence[k][variable][j] - (corrected_intercept[i] + corrected_slope[i] * seq_index_parameter[j]);
              residual_mean += diff;
              square_sum += diff * diff;
            }
          }
          residual_mean /= (change_point[i + 1] - change_point[i]);

          for (j = change_point[i];j < change_point[i + 1];j++) {
            for (k = 0;k < nb_sequence;k++) {
              diff = int_sequence[k][variable][j] - (corrected_intercept[i] + corrected_slope[i] * seq_index_parameter[j]) - residual_mean;
              residual_square_sum += diff * diff;
            }
          }
        }

        else {
          for (j = change_point[i];j < change_point[i + 1];j++) {
            for (k = 0;k < nb_sequence;k++) {
              diff = real_sequence[k][variable][j] - (corrected_intercept[i] + corrected_slope[i] * seq_index_parameter[j]);
              residual_mean += diff;
              square_sum += diff * diff;
            }
          }
          residual_mean /= (change_point[i + 1] - change_point[i]);

          for (j = change_point[i];j < change_point[i + 1];j++) {
            for (k = 0;k < nb_sequence;k++) {
              diff = real_sequence[k][variable][j] - (corrected_intercept[i] + corrected_slope[i] * seq_index_parameter[j]) - residual_mean;
              residual_square_sum += diff * diff;
            }
          }
        }

        if (model_type == INTERCEPT_SLOPE_CHANGE) {
          global_square_sum += square_sum;
          global_residual_square_sum += residual_square_sum;
        }

        variance[i] = residual_square_sum / (nb_sequence * (change_point[i + 1] - change_point[i]) - 2);

#       ifdef DEBUG
        cout << "\nTEST " << STAT_label[STATL_RESIDUAL] << " " << STAT_label[STATL_STANDARD_DEVIATION] << ": " 
             << sqrt(square_sum / (nb_sequence * (change_point[i + 1] - change_point[i]) - 2)) << " "
             << sqrt(residual_square_sum / (nb_sequence * (change_point[i + 1] - change_point[i]) - 2)) << " | " << residual_mean << endl;
#       endif

        if ((model_type == LINEAR_MODEL_CHANGE) && (likelihood != D_INF)) {
          if (square_sum > nb_sequence * (change_point[i + 1] - change_point[i]) * ROUNDOFF_ERROR) {
            likelihood -= ((double)(nb_sequence * (change_point[i + 1] - change_point[i])) / 2.) * (log(square_sum /
                            (nb_sequence * (change_point[i + 1] - change_point[i]))) + log(2 * M_PI) + 1);
          }
          else {
            likelihood = D_INF;
          }
        }
      }

      else {
        variance[i] = 0.;
      }
    }

    if (model_type == INTERCEPT_SLOPE_CHANGE) {
      if (global_square_sum > nb_sequence * length[0] * ROUNDOFF_ERROR) {
        likelihood -= ((double)(nb_sequence * length[0]) / 2.) * (log(global_square_sum /
                        (nb_sequence * length[0])) + log(2 * M_PI) + 1);

#       ifdef DEBUG
        cout << "\nTEST " << STAT_label[STATL_RESIDUAL] << " " << STAT_label[STATL_STANDARD_DEVIATION] << ": " 
             << sqrt(global_square_sum / (nb_sequence * length[0] - 2 * nb_segment)) << " "
             << sqrt(global_residual_square_sum / (nb_sequence * length[0] - 2 * nb_segment)) << endl;
#       endif

      }
      else {
        likelihood = D_INF;
      }
    }
  }

  os << "\n" << SEQ_label[SEQL_SEGMENT] << " " << STAT_label[STATL_INTERCEPT] << ", "
     << STAT_label[STATL_SLOPE];
  if (model_type == LINEAR_MODEL_CHANGE) {
    os << ", " << STAT_label[STATL_RESIDUAL] << " " << STAT_label[STATL_STANDARD_DEVIATION];
  }
  os << ": ";
  for (i = 0;i < nb_segment;i++) {
    os << corrected_intercept[i] << " " << corrected_slope[i];
    if (model_type == LINEAR_MODEL_CHANGE) {
      os  << " " << sqrt(variance[i]);
    }
    if (i < nb_segment - 1) {
      os << " | ";
    }
  }

  if (model_type == INTERCEPT_SLOPE_CHANGE) {
    os << "   " << STAT_label[STATL_RESIDUAL] << " " << STAT_label[STATL_STANDARD_DEVIATION] << ": " << sqrt(global_variance);
  }
  os << endl;

  os << SEQ_label[SEQL_PIECEWISE_LINEAR_FUNCTION] << ": ";
  for (i = 0;i < nb_segment;i++) {
    os << corrected_intercept[i] + corrected_slope[i] * seq_index_parameter[change_point[i]] << " -> ";
    if (i < nb_segment - 1) {
      os << corrected_intercept[i] + corrected_slope[i] * seq_index_parameter[change_point[i + 1]] << " | ";
    }
    else {
      os << corrected_intercept[i] + corrected_slope[i] * seq_index_parameter[change_point[i + 1] - 1] << endl;
    }
  }

  delete [] predicted_value;
  delete [] variance;

  return likelihood;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Output of a segmentation of a single sequence or a sample of sequences.
 *
 *  \param[in] nb_segment      number of segments,
 *  \param[in] model_type      segment model types,
 *  \param[in] common_contrast flag contrast functions common to the individuals,
 *  \param[in] display         flag for displaying the segmentation,
 *  \param[in] output          output type (sequence or residuals),
 *  \param[in] ichange_point   change points,
 *  \param[in] continuity      flag continuous piecewise linear function.
 *
 *  \return                    Sequences object.
 */
/*--------------------------------------------------------------*/

Sequences* Sequences::segmentation_output(int nb_segment , segment_model *model_type ,
                                          bool common_contrast , bool display , sequence_type output ,
                                          int *ichange_point , bool continuity)

{
  bool *piecewise_function_flag;
  int i , j , k , m , n;
  int inb_variable , min_identifier , max_identifier , *iidentifier , *ilength , *change_point ,
       *seq_index_parameter = NULL;
  variable_nature *itype;
  double likelihood , corrected_likelihood , diff , buff , change_point_amplitude , mean_absolute_deviation ,
         *global_variance , ***mean , ***variance , ***index_parameter_mean , ***intercept , ***slope ,
         ***correlation , ***slope_standard_deviation , ***corrected_intercept , ***corrected_slope ,
         ***autoregressive_coeff , ***determination_coeff;
  long double ***index_parameter_variance;
  Test *test;
  Sequences *seq;


  if (ichange_point) {
    change_point = ichange_point;
  }

  else {
    change_point = new int[nb_segment + 1];

    change_point[0] = 0;
    i = 1;
    for (j = 1;j < length[0];j++) {
      if (int_sequence[0][0][j] != int_sequence[0][0][j - 1]) {
        change_point[i++] = j;
      }
    }
    change_point[i] = length[0];
  }

  mean = new double**[nb_variable];
  variance = new double**[nb_variable];
  intercept = new double**[nb_variable];
  slope = new double**[nb_variable];
  correlation = new double**[nb_variable];
  slope_standard_deviation = new double**[nb_variable];
  index_parameter_mean = new double**[nb_variable];
  index_parameter_variance = new long double**[nb_variable];
  autoregressive_coeff = new double**[nb_variable];
  determination_coeff = new double**[nb_variable];

  global_variance = NULL;

  if (continuity) {
    corrected_intercept = new double**[nb_variable];
    corrected_slope = new double**[nb_variable];
  }

  if (index_param_type == IMPLICIT_TYPE) {
    seq_index_parameter = new int[length[0]];
    for (j = 0;j < length[0];j++) {
      seq_index_parameter[j] = j;
    }
  }
  else {
    seq_index_parameter = index_parameter[0];
  }

  for (i = 1;i < nb_variable;i++) {
    mean[i] = NULL;
    intercept[i] = NULL;
    slope[i] = NULL;
    correlation[i] = NULL;
    slope_standard_deviation[i] = NULL;

    if (continuity) {
      corrected_intercept[i] = NULL;
      corrected_slope[i] = NULL;
    }

    if ((model_type[i - 1] == POISSON_CHANGE) || (model_type[i - 1] == NEGATIVE_BINOMIAL_0_CHANGE) ||
        (model_type[i - 1] == NEGATIVE_BINOMIAL_1_CHANGE) || (model_type[i - 1] == GAUSSIAN_CHANGE) ||
        (model_type[0] == MEAN_CHANGE) || (model_type[i - 1] == VARIANCE_CHANGE) ||
        (model_type[i - 1] == BAYESIAN_POISSON_CHANGE) || (model_type[i - 1] == BAYESIAN_GAUSSIAN_CHANGE)) {
      if (common_contrast) {
        mean[i] = new double*[1];
        mean[i][0] = new double[nb_segment];
        variance[i] = new double*[4];
        for (j = 0;j < 4;j++) {
          variance[i][j] = new double[nb_segment];
        }
      }

      else {
        mean[i] = new double*[nb_sequence];
        variance[i] = new double*[nb_sequence];

        for (j = 0;j < nb_sequence;j++) {
          mean[i][j] = new double[nb_segment];
          variance[i][j] = new double[nb_segment];
        }
      }
    }

    else if ((model_type[i - 1] == LINEAR_MODEL_CHANGE) || (model_type[0] == INTERCEPT_SLOPE_CHANGE)) {
      if (common_contrast) {
        intercept[i] = new double*[1];
        intercept[i][0] = new double[nb_segment];
        slope[i] = new double*[1];
        slope[i][0] = new double[nb_segment];
        variance[i] = new double*[1];
        variance[i][0] = new double[nb_segment];
        correlation[i] = new double*[1];
        correlation[i][0] = new double[nb_segment];
        slope_standard_deviation[i] = new double*[1];
        slope_standard_deviation[i][0] = new double[nb_segment];
        index_parameter_mean[i] = new double*[1];
        index_parameter_mean[i][0] = new double[nb_segment];
        index_parameter_variance[i] = new long double*[1];
        index_parameter_variance[i][0] = new long double[nb_segment];

        if (continuity) {
          corrected_intercept[i] = new double*[1];
          corrected_intercept[i][0] = new double[nb_segment];
          corrected_slope[i] = new double*[1];
          corrected_slope[i][0] = new double[nb_segment];
        }
      }

      else {
        intercept[i] = new double*[nb_sequence];
        slope[i] = new double*[nb_sequence];
        variance[i] = new double*[nb_sequence];
        correlation[i] = new double*[nb_sequence];
        slope_standard_deviation[i] = new double*[nb_sequence];
        index_parameter_mean[i] = new double*[nb_sequence];
        index_parameter_variance[i] = new long double*[nb_sequence];

        for (j = 0;j < nb_sequence;j++) {
          intercept[i][j] = new double[nb_segment];
          slope[i][j] = new double[nb_segment];
          variance[i][j] = new double[nb_segment];
          correlation[i][j] = new double[nb_segment];
          slope_standard_deviation[i][j] = new double[nb_segment];
          index_parameter_mean[i][j] = new double[nb_segment];
          index_parameter_variance[i][j] = new long double[nb_segment];
        }
      }
    }

    else if ((model_type[i - 1] == AUTOREGRESSIVE_MODEL_CHANGE) || (model_type[i - 1] == STATIONARY_AUTOREGRESSIVE_MODEL_CHANGE)) {
      if (common_contrast) {
        mean[i] = new double*[1];
        mean[i][0] = new double[nb_segment];
        autoregressive_coeff[i] = new double*[1];
        autoregressive_coeff[i][0] = new double[nb_segment];
        variance[i] = new double*[1];
        variance[i][0] = new double[nb_segment];
        determination_coeff[i] = new double*[1];
        determination_coeff[i][0] = new double[nb_segment];
      }

      else {
        mean[i] = new double*[nb_sequence];
        autoregressive_coeff[i] = new double*[nb_sequence];
        variance[i] = new double*[nb_sequence];
        determination_coeff[i] = new double*[nb_sequence];

        for (j = 0;j < nb_sequence;j++) {
          mean[i][j] = new double[nb_segment];
          autoregressive_coeff[i][j] = new double[nb_segment];
          variance[i][j] = new double[nb_segment];
          determination_coeff[i][j] = new double[nb_segment];
        }
      }
    }

    if ((((i == 1) && ((model_type[0] == MEAN_CHANGE) || (model_type[0] == INTERCEPT_SLOPE_CHANGE))) ||
         (model_type[i - 1] == GAUSSIAN_CHANGE) || (model_type[i - 1] == LINEAR_MODEL_CHANGE) ||
         (model_type[i - 1] == AUTOREGRESSIVE_MODEL_CHANGE) || (model_type[i - 1] == STATIONARY_AUTOREGRESSIVE_MODEL_CHANGE) ||
         (model_type[i - 1] == BAYESIAN_GAUSSIAN_CHANGE)) && (!global_variance)) {
      global_variance = new double[nb_variable];
    }
  }

  if (output == SEQUENCE) {
    piecewise_function_flag = new bool[nb_variable];

    piecewise_function_flag[0] = false;
    for (i = 1;i < nb_variable;i++) {
      if ((model_type[i - 1] == POISSON_CHANGE) || (model_type[i - 1] == NEGATIVE_BINOMIAL_0_CHANGE) ||
          (model_type[i - 1] == NEGATIVE_BINOMIAL_1_CHANGE) || (model_type[i - 1] == GAUSSIAN_CHANGE) ||
          (model_type[0] == MEAN_CHANGE) || (model_type[i - 1] == LINEAR_MODEL_CHANGE) ||
          (model_type[0] == INTERCEPT_SLOPE_CHANGE) || (model_type[i - 1] == AUTOREGRESSIVE_MODEL_CHANGE) ||
          (model_type[i - 1] == STATIONARY_AUTOREGRESSIVE_MODEL_CHANGE) ||
          (model_type[i - 1] == BAYESIAN_POISSON_CHANGE) || (model_type[i - 1] == BAYESIAN_GAUSSIAN_CHANGE)) {
        piecewise_function_flag[i] = true;
      }
      else {
        piecewise_function_flag[i] = false;
      }
    }

    seq = new Sequences(*this , piecewise_function_flag);
  }

  else if (output == SEQUENCE_SAMPLE) {
    iidentifier = new int[nb_sequence + 2];

    min_identifier = identifier[0];
    for (i = 0;i < nb_sequence;i++) {
      if (identifier[i] < min_identifier) {
        min_identifier = identifier[i];
      }

      iidentifier[i + 1] = identifier[i];
    }

    iidentifier[0] = min_identifier - 1;

    max_identifier = identifier[nb_sequence - 1];
    for (i = 0;i < nb_sequence - 1;i++) {
      if (identifier[i] > max_identifier) {
        max_identifier = identifier[i];
      }
    }

    iidentifier[nb_sequence + 1] = max_identifier + 1;

    ilength = new int[nb_sequence + 2];
    for (i = 0;i < nb_sequence + 2;i++) {
      ilength[i] = length[0];
    }

    itype = new variable_nature[nb_variable - 1];
    for (i = 0;i < nb_variable - 1;i++) {
      itype[i] = REAL_VALUE;
    }

    seq = new Sequences(nb_sequence + 2 , iidentifier , ilength , NULL ,
                        index_param_type , nb_variable - 1 , itype);
    delete [] iidentifier;
    delete [] ilength;
    delete [] itype;
  }

  else if (output == ABSOLUTE_RESIDUAL) {
    inb_variable = 1 + 2 * (nb_variable - 1);
    itype = new variable_nature[inb_variable];

    itype[0] = type[0];
    i = 1;
    for (j = 1;j < nb_variable;j++) {
      itype[i++] = REAL_VALUE;
      itype[i++] = AUXILIARY;
    }

    seq = new Sequences(nb_sequence , identifier , length , vertex_identifier ,
                        index_param_type , inb_variable , itype);
    delete [] itype;
  }

  for (i = 1;i < nb_variable;i++) {
    if ((model_type[i - 1] == POISSON_CHANGE) || (model_type[i - 1] == NEGATIVE_BINOMIAL_0_CHANGE) ||
        (model_type[i - 1] == NEGATIVE_BINOMIAL_1_CHANGE) || (model_type[i - 1] == GAUSSIAN_CHANGE) ||
        (model_type[0] == MEAN_CHANGE) || (model_type[i - 1] == VARIANCE_CHANGE) ||
        (model_type[i - 1] == LINEAR_MODEL_CHANGE) || (model_type[0] == INTERCEPT_SLOPE_CHANGE) ||
        (model_type[i - 1] == AUTOREGRESSIVE_MODEL_CHANGE) || (model_type[i - 1] == STATIONARY_AUTOREGRESSIVE_MODEL_CHANGE) ||
        (model_type[i - 1] == BAYESIAN_POISSON_CHANGE) || (model_type[i - 1] == BAYESIAN_GAUSSIAN_CHANGE)) {
      likelihood = piecewise_linear_function((nb_sequence == 1 ? 0 : I_DEFAULT) , i , nb_segment ,
                                              model_type[i - 1] , common_contrast , change_point ,
                                              seq_index_parameter , NULL , mean[i] , variance[i] ,
                                              global_variance , intercept[i] , slope[i] ,
                                              autoregressive_coeff[i] , correlation[i] ,
                                              slope_standard_deviation[i] , index_parameter_mean[i] ,
                                              index_parameter_variance[i] , determination_coeff[i]);

      if ((display) && (((i == 1) && ((model_type[0] == MEAN_CHANGE) || (model_type[0] == INTERCEPT_SLOPE_CHANGE))) ||
           (model_type[i - 1] == GAUSSIAN_CHANGE) || (model_type[i - 1] == VARIANCE_CHANGE) ||
           (model_type[i - 1] == LINEAR_MODEL_CHANGE) || (model_type[i - 1] == AUTOREGRESSIVE_MODEL_CHANGE) ||
           (model_type[i - 1] == STATIONARY_AUTOREGRESSIVE_MODEL_CHANGE))) {
        cout << "\n2 * " << STAT_label[STATL_LIKELIHOOD] << ": " << 2 * likelihood << endl;
      }
    }
  }

  if (display) {
    if (!ichange_point) {
      cout << (nb_segment == 2 ? SEQ_label[SEQL_CHANGE_POINT] : SEQ_label[SEQL_CHANGE_POINTS]) << ": ";

      for (i = 1;i < nb_segment;i++) {
        cout << seq_index_parameter[change_point[i]];
        if (i < nb_segment - 1) {
          cout << ", ";
        }
      }
    }

    if ((index_interval) && (index_interval->variance > 0.)) {
      if (!ichange_point) {
        cout << "   ";
      }
      cout << SEQ_label[SEQL_SEGMENT_SAMPLE_SIZE] << ": ";
      for (i = 0;i < nb_segment;i++) {
        cout << nb_sequence * (change_point[i + 1] - change_point[i]);
        if (i < nb_segment - 1) {
          cout << ", ";
        }
      }
      cout << endl;
    }

    else if (!ichange_point) {
      cout << endl;
    }

    if (nb_variable > 2) {
      cout << "\n";
    }

    for (i = 1;i < nb_variable;i++) {
      piecewise_linear_function_ascii_print(cout , (nb_sequence == 1 ? 0 : I_DEFAULT) , i , nb_segment , model_type[i - 1] ,
                                            common_contrast , change_point , seq_index_parameter ,
                                            mean[i] , variance[i] , intercept[i] , slope[i] ,
                                            autoregressive_coeff[i] , correlation[i] , slope_standard_deviation[i] ,
                                            index_parameter_mean[i] , index_parameter_variance[i] , determination_coeff[i]);

      if ((model_type[i - 1] == GAUSSIAN_CHANGE) || (model_type[0] == MEAN_CHANGE) ||
          (model_type[i - 1] == AUTOREGRESSIVE_MODEL_CHANGE) || (model_type[i - 1] == BAYESIAN_GAUSSIAN_CHANGE)) {
        if (nb_segment > 1) {
          change_point_amplitude = 0.;

          if ((nb_sequence == 1) || (common_contrast)) {
            for (j = 1;j < nb_segment;j++) {
              change_point_amplitude += fabs(mean[i][0][j] - mean[i][0][j - 1]);
            }
            change_point_amplitude /= (nb_segment - 1);
          }

          else {
            for (j = 0;j < nb_sequence;j++) {
              for (k = 1;k < nb_segment;k++) {
                change_point_amplitude += fabs(mean[i][j][k] - mean[i][j][k - 1]);
              }
            }
            change_point_amplitude /= (nb_sequence * (nb_segment - 1));
          }

          cout << STAT_label[STATL_MEAN] << " " << SEQ_label[SEQL_CHANGE_POINT_AMPLITUDE] << ": "
               << change_point_amplitude << "   ";
        }

        cout << SEQ_label[SEQL_ROOT_MEAN_SQUARE_ERROR] << ": " << sqrt(global_variance[i]);
        if (nb_segment > 1) {
          cout << "   " << STAT_label[STATL_RATIO] << ": "
               << change_point_amplitude / sqrt(global_variance[i]);
        }
        cout << endl;
      }

      else if ((model_type[i - 1] == LINEAR_MODEL_CHANGE) || (model_type[i - 1] == STATIONARY_AUTOREGRESSIVE_MODEL_CHANGE)) {
        cout << SEQ_label[SEQL_ROOT_MEAN_SQUARE_ERROR] << ": " << sqrt(global_variance[i]) << endl;
      }

      else if (model_type[0] == MEAN_CHANGE) {
        cout << SEQ_label[SEQL_GLOBAL_STANDARD_DEVIATION] << ": "  << sqrt(global_variance[i]) << endl;
      }
      else if (model_type[0] == INTERCEPT_SLOPE_CHANGE) {
        cout << SEQ_label[SEQL_GLOBAL_RESIDUAL_STANDARD_DEVIATION] << ": "  << sqrt(global_variance[i]) << endl;
      }

      if (continuity) {
        corrected_likelihood = continuous_piecewise_linear_function(cout , (nb_sequence == 1 ? 0 : I_DEFAULT) , i ,
                                                                    nb_segment , model_type[i - 1] , common_contrast ,
                                                                    change_point , seq_index_parameter , intercept[i][0] ,
                                                                    slope[i][0] , corrected_intercept[i][0] , corrected_slope[i][0]);

        cout << "2 * " << STAT_label[STATL_LIKELIHOOD] << ": "
              << 2 * corrected_likelihood << " | " << 2 * likelihood << endl;
      }
    }
  }

  switch (output) {

  case SEQUENCE : {
    if (common_contrast) {
      for (i = 0;i < nb_sequence;i++) {
        j = 1;
        for (k = 1;k < nb_variable;k++) {
          j++;
          if (piecewise_function_flag[k]) {
            if ((model_type[k - 1] == LINEAR_MODEL_CHANGE) || (model_type[k - 1] == INTERCEPT_SLOPE_CHANGE)) {
              if (continuity) {
                for (m = 0;m < nb_segment;m++) {
                  for (n = change_point[m];n < change_point[m + 1];n++) {
                    seq->real_sequence[i][j][n] = corrected_intercept[k][0][m] + corrected_slope[k][0][m] * seq_index_parameter[n];
                  }
                }
              }

              else{
                for (m = 0;m < nb_segment;m++) {
                  for (n = change_point[m];n < change_point[m + 1];n++) {
                    seq->real_sequence[i][j][n] = intercept[k][0][m] + slope[k][0][m] * seq_index_parameter[n];
                  }
                }
              }
            }

            else if ((model_type[k - 1] == AUTOREGRESSIVE_MODEL_CHANGE) || (model_type[k - 1] == STATIONARY_AUTOREGRESSIVE_MODEL_CHANGE)) {
              if (type[k] != REAL_VALUE) {
                for (m = 0;m < nb_segment;m++) {
                  seq->real_sequence[i][j][change_point[m]] = mean[k][0][m];
                  for (n = change_point[m] + 1;n < change_point[m + 1];n++) {
                    seq->real_sequence[i][j][n] = mean[k][0][m] + autoregressive_coeff[k][0][m] * (int_sequence[i][k][n - 1] - mean[k][0][m]);
                  }
                }
              }

              else {
                for (m = 0;m < nb_segment;m++) {
                  seq->real_sequence[i][j][change_point[m]] = mean[k][0][m];
                  for (n = change_point[m] + 1;n < change_point[m + 1];n++) {
                    seq->real_sequence[i][j][n] = mean[k][0][m] + autoregressive_coeff[k][0][m] * (real_sequence[i][k][n - 1] - mean[k][0][m]);
                  }
                }
              }
            }

            else {
              for (m = 0;m < nb_segment;m++) {
                for (n = change_point[m];n < change_point[m + 1];n++) {
                  seq->real_sequence[i][j][n] = mean[k][0][m];
                }
              }
            }

            j++;
          }
        }
      }
    }

    else {
      for (i = 0;i < nb_sequence;i++) {
        j = 1;
        for (k = 1;k < nb_variable;k++) {
          j++;
          if (piecewise_function_flag[k]) {
            if ((model_type[k - 1] == LINEAR_MODEL_CHANGE) || (model_type[k - 1] == INTERCEPT_SLOPE_CHANGE)) {
              if (continuity) {
                for (m = 0;m < nb_segment;m++) {
                  for (n = change_point[m];n < change_point[m + 1];n++) {
                    seq->real_sequence[i][j][n] = corrected_intercept[k][i][m] + corrected_slope[k][i][m] * seq_index_parameter[n];
                  }
                }
              }

              else {
                for (m = 0;m < nb_segment;m++) {
                  for (n = change_point[m];n < change_point[m + 1];n++) {
                    seq->real_sequence[i][j][n] = intercept[k][i][m] + slope[k][i][m] * seq_index_parameter[n];
                  }
                }
              }
            }

            else if ((model_type[k - 1] == AUTOREGRESSIVE_MODEL_CHANGE) || (model_type[k - 1] == STATIONARY_AUTOREGRESSIVE_MODEL_CHANGE)) {
              if (type[k] != REAL_VALUE) {
                for (m = 0;m < nb_segment;m++) {
                  seq->real_sequence[i][j][change_point[m]] = mean[k][i][m];
                  for (n = change_point[m] + 1;n < change_point[m + 1];n++) {
                    seq->real_sequence[i][j][n] = mean[k][i][m] + autoregressive_coeff[k][i][m] * (int_sequence[i][k][n - 1] - mean[k][i][m]);
                  }
                }
              }

              else {
                for (m = 0;m < nb_segment;m++) {
                  seq->real_sequence[i][j][change_point[m]] = mean[k][i][m];
                  for (n = change_point[m] + 1;n < change_point[m + 1];n++) {
                    seq->real_sequence[i][j][n] = mean[k][i][m] + autoregressive_coeff[k][i][m] * (real_sequence[i][k][n - 1] - mean[k][i][m]);
                  }
                }
              }
            }

            else {
              for (m = 0;m < nb_segment;m++) {
                for (n = change_point[m];n < change_point[m + 1];n++) {
                  seq->real_sequence[i][j][n] = mean[k][i][m];
                }
              }
            }
            j++;
          }
        }
      }
    }
    break;
  }

  case SEQUENCE_SAMPLE : {

    // copy of index parameters

    if (index_parameter) {
      for (i = 0;i < seq->nb_sequence;i++) {
        for (j = 0;j < (index_param_type == POSITION ? length[0] + 1 : length[0]);j++) {
          seq->index_parameter[i][j] = index_parameter[0][j];
        }
      }

      seq->build_index_parameter_frequency_distribution();
      seq->index_interval_computation();
    }

    // copy of sequences

    for (i = 0;i < nb_sequence;i++) {
      for (j = 1;j < nb_variable;j++) {
        if (type[j] != REAL_VALUE) {
          for (k = 0;k < length[i];k++) {
            seq->real_sequence[i + 1][j - 1][k] = int_sequence[i][j][k];
          }
        }

        else {
          for (k = 0;k < length[i];k++) {
            seq->real_sequence[i + 1][j - 1][k] = real_sequence[i][j][k];
          }
        }
      }
    }

    for (i = 1;i < nb_variable;i++) {
      if ((model_type[i - 1] == LINEAR_MODEL_CHANGE) || (model_type[i - 1] == INTERCEPT_SLOPE_CHANGE)) {
        if (continuity) {
          for (j = 0;j < nb_segment;j++) {
            for (k = change_point[j];k < change_point[j + 1];k++) {
              seq->real_sequence[0][i - 1][k] = corrected_intercept[i][0][j] + corrected_slope[i][0][j] * seq_index_parameter[k];
            }
          }
        }

        else{
          for (j = 0;j < nb_segment;j++) {
            for (k = change_point[j];k < change_point[j + 1];k++) {
              seq->real_sequence[0][i - 1][k] = intercept[i][0][j] + slope[i][0][j] * seq_index_parameter[k];
            }
          }
        }
      }

      else {
        for (j = 0;j < nb_segment;j++) {
          for (k = change_point[j];k < change_point[j + 1];k++) {
            seq->real_sequence[0][i - 1][k] = mean[i][0][j];
          }
        }
      }

      if ((model_type[i - 1] == MEAN_CHANGE) || (model_type[i - 1] == INTERCEPT_SLOPE_CHANGE)) {
        buff = sqrt(global_variance[0]);
        for (j = 0;j < length[i];j++) {
          seq->real_sequence[nb_sequence + 1][i - 1][j] = buff;
        }
      }

      else {
        for (j = 0;j < nb_segment;j++) {
          buff = sqrt(variance[i][0][j]);
          for (k = change_point[j];k < change_point[j + 1];k++) {
            seq->real_sequence[nb_sequence + 1][i - 1][k] = buff;
          }
        }
      }
    }

    for (i = 0;i < seq->nb_variable;i++) {
      seq->min_value_computation(i);
      seq->max_value_computation(i);

      seq->build_marginal_histogram(i);
    }
    break;
  }

  // residual computation

  case SUBTRACTION_RESIDUAL : {
    if (common_contrast) {
      for (i = 0;i < nb_sequence;i++) {
        for (j = 1;j < nb_variable;j++) {
          if (type[j] != REAL_VALUE) {
            real_sequence[i][j] = new double[length[i]];

            if ((model_type[j - 1] == LINEAR_MODEL_CHANGE) || (model_type[j - 1] == INTERCEPT_SLOPE_CHANGE)) {
              for (k = 0;k < nb_segment;k++) {
                for (m = change_point[k];m < change_point[k + 1];m++) {
                  real_sequence[i][j][m] = int_sequence[i][j][m] - (intercept[j][0][k] + slope[j][0][k] * seq_index_parameter[m]);
                }
              }
            }

            else if ((model_type[j - 1] == AUTOREGRESSIVE_MODEL_CHANGE) || (model_type[j - 1] == STATIONARY_AUTOREGRESSIVE_MODEL_CHANGE)) {
              for (k = 0;k < nb_segment;k++) {
                real_sequence[i][j][change_point[k]] = int_sequence[i][j][change_point[k]] - mean[j][0][k];
                for (m = change_point[k] + 1;m < change_point[k + 1];m++) {
                  real_sequence[i][j][m] = int_sequence[i][j][m] - (mean[j][0][k] + autoregressive_coeff[j][0][k] *
                                           (int_sequence[i][j][m - 1] - mean[j][0][k]));
                }
              }
            }

            else {
              for (k = 0;k < nb_segment;k++) {
                for (m = change_point[k];m < change_point[k + 1];m++) {
                  real_sequence[i][j][m] = int_sequence[i][j][m] - mean[j][0][k];
                }
              }
            }

            delete [] int_sequence[i][j];
            int_sequence[i][j] = NULL;
          }

          else {
            if ((model_type[j - 1] == LINEAR_MODEL_CHANGE) || (model_type[j - 1] == INTERCEPT_SLOPE_CHANGE)) {
              for (k = 0;k < nb_segment;k++) {
                for (m = change_point[k];m < change_point[k + 1];m++) {
                  real_sequence[i][j][m] -= (intercept[j][0][k] + slope[j][0][k] * seq_index_parameter[m]);
                }
              }
            }

            else if ((model_type[j - 1] == AUTOREGRESSIVE_MODEL_CHANGE) || (model_type[j - 1] == STATIONARY_AUTOREGRESSIVE_MODEL_CHANGE)) {
              for (k = 0;k < nb_segment;k++) {
                for (m = change_point[k + 1] - 1;m > change_point[k];m--) {
                  real_sequence[i][j][m] -= (mean[j][0][k] + autoregressive_coeff[j][0][k] *
                                             (real_sequence[i][j][m - 1] - mean[j][0][k]));
                }
                real_sequence[i][j][change_point[k]] -= mean[j][0][k];
              }
            }

            else {
              for (k = 0;k < nb_segment;k++) {
                for (m = change_point[k];m < change_point[k + 1];m++) {
                  real_sequence[i][j][m] -= mean[j][0][k];
                }
              }
            }
          }
        }
      }
    }

    else {
      for (i = 0;i < nb_sequence;i++) {
        for (j = 1;j < nb_variable;j++) {
          if (type[j] != REAL_VALUE) {
            real_sequence[i][j] = new double[length[i]];

            if ((model_type[j - 1] == LINEAR_MODEL_CHANGE) || (model_type[j - 1] == INTERCEPT_SLOPE_CHANGE)) {
              for (k = 0;k < nb_segment;k++) {
                for (m = change_point[k];m < change_point[k + 1];m++) {
                  real_sequence[i][j][m] = int_sequence[i][j][m] - (intercept[j][i][k] + slope[j][i][k] * seq_index_parameter[m]);
                }
              }
            }

            else if ((model_type[j - 1] == AUTOREGRESSIVE_MODEL_CHANGE) || (model_type[j - 1] == STATIONARY_AUTOREGRESSIVE_MODEL_CHANGE)) {
              for (k = 0;k < nb_segment;k++) {
                real_sequence[i][j][change_point[k]] = int_sequence[i][j][change_point[k]] - mean[j][i][k];
                for (m = change_point[k] + 1;m < change_point[k + 1];m++) {
                  real_sequence[i][j][m] = int_sequence[i][j][m] - (mean[j][i][k] + autoregressive_coeff[j][i][k] *
                                            (int_sequence[i][j][m - 1] - mean[j][i][k]));
                }
              }
            }

            else {
              for (k = 0;k < nb_segment;k++) {
                for (m = change_point[k];m < change_point[k + 1];m++) {
                  real_sequence[i][j][m] = int_sequence[i][j][m] - mean[j][i][k];
                }
              }
            }

            delete [] int_sequence[i][j];
            int_sequence[i][j] = NULL;
          }

          else {
            if ((model_type[j - 1] == LINEAR_MODEL_CHANGE) || (model_type[j - 1] == INTERCEPT_SLOPE_CHANGE)) {
              for (k = 0;k < nb_segment;k++) {
                for (m = change_point[k];m < change_point[k + 1];m++) {
                  real_sequence[i][j][m] -= (intercept[j][i][k] + slope[j][i][k] * seq_index_parameter[m]);
                }
              }
            }

            else if ((model_type[j - 1] == AUTOREGRESSIVE_MODEL_CHANGE) || (model_type[j - 1] == STATIONARY_AUTOREGRESSIVE_MODEL_CHANGE)) {
              for (k = 0;k < nb_segment;k++) {
                for (m = change_point[k + 1] - 1;m > change_point[k];m--) {
                  real_sequence[i][j][m] -= (mean[j][i][k] + autoregressive_coeff[j][i][k] *
                                             (real_sequence[i][j][m - 1] - mean[j][i][k]));
                }
                real_sequence[i][j][change_point[k]] -= mean[j][i][k];
              }
            }

            else {
              for (k = 0;k < nb_segment;k++) {
                for (m = change_point[k];m < change_point[k + 1];m++) {
                  real_sequence[i][j][m] -= mean[j][i][k];
                }
              }
            }
          }
        }
      }
    }
    break;
  }

  case ABSOLUTE_RESIDUAL : {
    if (index_parameter) {
      for (i = 0;i < nb_sequence;i++) {
        for (j = 0;j < (index_param_type == POSITION ? length[i] + 1 : length[i]);j++) {
          seq->index_parameter[i][j] = index_parameter[i][j];
        }
      }
    }

    if (index_parameter_distribution) {
      seq->index_parameter_distribution = new FrequencyDistribution(*index_parameter_distribution);
    }
    if (index_interval) {
      seq->index_interval = new FrequencyDistribution(*index_interval);
    }

    seq->min_value[0] = min_value[0];
    seq->max_value[0] = max_value[0];
    seq->marginal_distribution[0] = new FrequencyDistribution(*marginal_distribution[0]);

    for (i = 0;i < nb_sequence;i++) {
      for (j = 0;j < length[i];j++) {
        seq->int_sequence[i][0][j] = int_sequence[i][0][j];
      }
    }

    if (common_contrast) {
      for (i = 0;i < nb_sequence;i++) {
        j = 1;
        for (k = 1;k < nb_variable;k++) {
          if (type[k] != REAL_VALUE) {
            if ((model_type[k - 1] == LINEAR_MODEL_CHANGE) || (model_type[k - 1] == INTERCEPT_SLOPE_CHANGE)) {
              for (m = 0;m < nb_segment;m++) {
                for (n = change_point[m];n < change_point[m + 1];n++) {
                  seq->real_sequence[i][j][n] = fabs(int_sequence[i][k][n] - (intercept[k][0][m] + slope[k][0][m] * seq_index_parameter[n]));
                }
              }
            }

            else if ((model_type[k - 1] == AUTOREGRESSIVE_MODEL_CHANGE) || (model_type[k - 1] == STATIONARY_AUTOREGRESSIVE_MODEL_CHANGE)) {
              for (m = 0;m < nb_segment;m++) {
                seq->real_sequence[i][j][change_point[m]] = fabs(int_sequence[i][k][change_point[m]] - mean[k][0][m]);
                for (n = change_point[m] + 1;n < change_point[m + 1];n++) {
                  seq->real_sequence[i][j][n] = fabs(int_sequence[i][k][n] - (mean[k][0][m] + autoregressive_coeff[k][0][m] *
                                                     (int_sequence[i][k][n - 1] - mean[k][0][m])));
                }
              }
            }

            else {
              for (m = 0;m < nb_segment;m++) {
                for (n = change_point[m];n < change_point[m + 1];n++) {
                  seq->real_sequence[i][j][n] = fabs(int_sequence[i][k][n] - mean[k][0][m]);
                }
              }
            }
          }

          else {
            if ((model_type[k - 1] == LINEAR_MODEL_CHANGE) || (model_type[k - 1] == INTERCEPT_SLOPE_CHANGE)) {
              for (m = 0;m < nb_segment;m++) {
                for (n = change_point[m];n < change_point[m + 1];n++) {
                  seq->real_sequence[i][j][n] = fabs(real_sequence[i][k][n] - (intercept[k][0][m] + slope[k][0][m] * seq_index_parameter[n]));
                }
              }
            }

            else if ((model_type[k - 1] == AUTOREGRESSIVE_MODEL_CHANGE) || (model_type[k - 1] == STATIONARY_AUTOREGRESSIVE_MODEL_CHANGE)) {
              for (m = 0;m < nb_segment;m++) {
                seq->real_sequence[i][j][change_point[m]] = fabs(real_sequence[i][k][change_point[m]] - mean[k][0][m]);
                for (n = change_point[m] + 1;n < change_point[m + 1];n++) {
                  seq->real_sequence[i][j][n] = fabs(real_sequence[i][k][n] - (mean[k][0][m] + autoregressive_coeff[k][0][m] *
                                                     (real_sequence[i][k][n - 1] - mean[k][0][m])));
                }
              }
            }

            else {
              for (m = 0;m < nb_segment;m++) {
                for (n = change_point[m];n < change_point[m + 1];n++) {
                  seq->real_sequence[i][j][n] = fabs(real_sequence[i][k][n] - mean[k][0][m]);
                }
              }
            }
          }

          j += 2;
        }
      }

      i = 1;
      for (j = 1;j < nb_variable;j++) {
        for (k = 0;k < nb_segment;k++) {
          mean_absolute_deviation = 0.;

          if ((model_type[j - 1] == AUTOREGRESSIVE_MODEL_CHANGE) || (model_type[j - 1] == STATIONARY_AUTOREGRESSIVE_MODEL_CHANGE)) {
            for (m = change_point[k] + 1;m < change_point[k + 1];m++) {
              for (n = 0;n < nb_sequence;n++) {
                mean_absolute_deviation += seq->real_sequence[n][i][m];
              }
            }
//            mean_absolute_deviation /= (nb_sequence * (change_point[k + 1] - change_point[k] - 1) - 2);
            mean_absolute_deviation /= (nb_sequence * (change_point[k + 1] - change_point[k] - 1) - 1);
          }

          else {
            for (m = change_point[k];m < change_point[k + 1];m++) {
              for (n = 0;n < nb_sequence;n++) {
                mean_absolute_deviation += seq->real_sequence[n][i][m];
              }
            }
            if ((model_type[j - 1] == LINEAR_MODEL_CHANGE) || (model_type[j - 1] == INTERCEPT_SLOPE_CHANGE)) {
              mean_absolute_deviation /= (nb_sequence * (change_point[k + 1] - change_point[k]) - 2);
            }
            else {
              mean_absolute_deviation /= (nb_sequence * (change_point[k + 1] - change_point[k]) - 1);
            }
          }

          for (m = change_point[k];m < change_point[k + 1];m++) {
            for (n = 0;n < nb_sequence;n++) {
              seq->real_sequence[n][i + 1][m] = mean_absolute_deviation;
            }
          }
        }

        i += 2;
      }
    }

    else {
      for (i = 0;i < nb_sequence;i++) {
        j = 1;
        for (k = 1;k < nb_variable;k++) {
          if (type[k] != REAL_VALUE) {
            if ((model_type[k - 1] == LINEAR_MODEL_CHANGE) || (model_type[k - 1] == INTERCEPT_SLOPE_CHANGE)) {
              for (m = 0;m < nb_segment;m++) {
                for (n = change_point[m];n < change_point[m + 1];n++) {
                  seq->real_sequence[i][j][n] = fabs(int_sequence[i][k][n] - (intercept[k][i][m] + slope[k][i][m] * seq_index_parameter[n]));
                }
              }
            }

            else if ((model_type[k - 1] == AUTOREGRESSIVE_MODEL_CHANGE) || (model_type[k - 1] == STATIONARY_AUTOREGRESSIVE_MODEL_CHANGE)) {
              for (m = 0;m < nb_segment;m++) {
                seq->real_sequence[i][j][change_point[m]] = fabs(int_sequence[i][k][change_point[m]] - mean[k][i][m]);
                for (n = change_point[m] + 1;n < change_point[m + 1];n++) {
                  seq->real_sequence[i][j][n] = fabs(int_sequence[i][k][n] - (mean[k][i][m] + autoregressive_coeff[k][i][m] *
                                                     (int_sequence[i][k][n - 1] - mean[k][i][m])));
                }
              }
            }

            else {
              for (m = 0;m < nb_segment;m++) {
                for (n = change_point[m];n < change_point[m + 1];n++) {
                  seq->real_sequence[i][j][n] = fabs(int_sequence[i][k][n] - mean[k][i][m]);
                }
              }
            }
          }

          else {
            if ((model_type[k - 1] == LINEAR_MODEL_CHANGE) || (model_type[k - 1] == INTERCEPT_SLOPE_CHANGE)) {
              for (m = 0;m < nb_segment;m++) {
                for (n = change_point[m];n < change_point[m + 1];n++) {
                  seq->real_sequence[i][j][n] = fabs(real_sequence[i][k][n] - (intercept[k][i][m] + slope[k][i][m] * seq_index_parameter[n]));
                }
              }
            }

            else if ((model_type[k - 1] == AUTOREGRESSIVE_MODEL_CHANGE) || (model_type[k - 1] == STATIONARY_AUTOREGRESSIVE_MODEL_CHANGE)) {
              for (m = 0;m < nb_segment;m++) {
                seq->real_sequence[i][j][change_point[m]] = fabs(real_sequence[i][k][change_point[m]] - mean[k][i][m]);
                for (n = change_point[m] + 1;n < change_point[m + 1];n++) {
                  seq->real_sequence[i][j][n] = fabs(real_sequence[i][k][n] - (mean[k][i][m] + autoregressive_coeff[k][i][m] *
                                                     (real_sequence[i][k][n - 1] - mean[k][i][m])));
                }
              }
            }

            else {
              for (m = 0;m < nb_segment;m++) {
                for (n = change_point[m];n < change_point[m + 1];n++) {
                  seq->real_sequence[i][j][n] = fabs(real_sequence[i][k][n] - mean[k][i][m]);
                }
              }
            }
          }

          for (m = 0;m < nb_segment;m++) {
            mean_absolute_deviation = 0.;

            if ((model_type[k - 1] == AUTOREGRESSIVE_MODEL_CHANGE) || (model_type[k - 1] == STATIONARY_AUTOREGRESSIVE_MODEL_CHANGE)) {
              for (n = change_point[m] + 1;n < change_point[m + 1];n++) {
                mean_absolute_deviation += seq->real_sequence[i][j][n];
              }
//              mean_absolute_deviation /= (change_point[m + 1] - change_point[m] - 3);
              mean_absolute_deviation /= (change_point[m + 1] - change_point[m] - 2);
            }

            else {
              for (n = change_point[m];n < change_point[m + 1];n++) {
                mean_absolute_deviation += seq->real_sequence[i][j][n];
              }
              if ((model_type[k - 1] == LINEAR_MODEL_CHANGE) || (model_type[k - 1] == INTERCEPT_SLOPE_CHANGE)) {
                mean_absolute_deviation /= (change_point[m + 1] - change_point[m] - 2);
              }
              else {
                mean_absolute_deviation /= (change_point[m + 1] - change_point[m] - 1);
              }
            }

            for (n = change_point[m];n < change_point[m + 1];n++) {
              seq->real_sequence[i][j + 1][n] = mean_absolute_deviation;
            }
          }

          j += 2;
        }
      }
    }
    break;
  }

  case DIVISION_RESIDUAL : {
    if (common_contrast) {
      for (i = 0;i < nb_sequence;i++) {
        for (j = 1;j < nb_variable;j++) {
          if (type[j] != REAL_VALUE) {
            real_sequence[i][j] = new double[length[i]];

            if ((model_type[j - 1] == LINEAR_MODEL_CHANGE) || (model_type[j - 1] == INTERCEPT_SLOPE_CHANGE)) {
              for (k = 0;k < nb_segment;k++) {
                for (m = change_point[k];m < change_point[k + 1];m++) {
                  if (intercept[j][0][k] + slope[j][0][k] * seq_index_parameter[m] != 0.) {
                    real_sequence[i][j][m] = int_sequence[i][j][m] / (intercept[j][0][k] + slope[j][0][k] * seq_index_parameter[m]);
                  }
                  else {
                    real_sequence[i][j][m] = int_sequence[i][j][m];
                  }
                }
              }
            }

            else if ((model_type[j - 1] == AUTOREGRESSIVE_MODEL_CHANGE) || (model_type[j - 1] == STATIONARY_AUTOREGRESSIVE_MODEL_CHANGE)) {
              for (k = 0;k < nb_segment;k++) {
                if (mean[j][0][k] != 0.) {
                  real_sequence[i][j][change_point[k]] = int_sequence[i][j][change_point[k]] / mean[j][0][k];
                }
                else {
                  real_sequence[i][j][change_point[k]] = int_sequence[i][j][change_point[k]];
                }
                for (m = change_point[k] + 1;m < change_point[k + 1];m++) {
                  if (mean[j][0][k] + autoregressive_coeff[j][0][k] * (int_sequence[i][j][m - 1] - mean[j][0][k]) != 0.) {
                    real_sequence[i][j][m] = int_sequence[i][j][m] / (mean[j][0][k] + autoregressive_coeff[j][0][k] *
                                             (int_sequence[i][j][m - 1] - mean[j][0][k]));
                  }
                  else {
                    real_sequence[i][j][m] = int_sequence[i][j][m];
                  }
                }
              }
            }

            else {
              for (k = 0;k < nb_segment;k++) {
                if (mean[j][0][k] != 0.) {
                  for (m = change_point[k];m < change_point[k + 1];m++) {
                    real_sequence[i][j][m] = int_sequence[i][j][m] / mean[j][0][k];
                  }
                }
                else {
                  for (m = change_point[k];m < change_point[k + 1];m++) {
                    real_sequence[i][j][m] = int_sequence[i][j][m];
                  }
                }
              }
            }

            delete [] int_sequence[i][j];
            int_sequence[i][j] = NULL;
          }

          else {
            if ((model_type[j - 1] == LINEAR_MODEL_CHANGE) || (model_type[j - 1] == INTERCEPT_SLOPE_CHANGE)) {
              for (k = 0;k < nb_segment;k++) {
                for (m = change_point[k];m < change_point[k + 1];m++) {
                  if (intercept[j][0][k] + slope[j][0][k] * seq_index_parameter[m] != 0.) {
                    real_sequence[i][j][m] /= (intercept[j][0][k] + slope[j][0][k] * seq_index_parameter[m]);
                  }
                }
              }
            }

            else if ((model_type[j - 1] == AUTOREGRESSIVE_MODEL_CHANGE) || (model_type[j - 1] == STATIONARY_AUTOREGRESSIVE_MODEL_CHANGE)) {
              for (k = 0;k < nb_segment;k++) {
                for (m = change_point[k + 1] - 1;m > change_point[k];m--) {
                  if (mean[j][0][k] + autoregressive_coeff[j][0][k] * (real_sequence[i][j][m - 1] - mean[j][0][k]) != 0.) {
                    real_sequence[i][j][m] /= (mean[j][0][k] + autoregressive_coeff[j][0][k] * (real_sequence[i][j][m - 1] - mean[j][0][k]));
                  }
                }
                if (mean[j][0][k] != 0.) {
                  real_sequence[i][j][change_point[k]] /= mean[j][0][k];
                }

              }
            }

            else {
              for (k = 0;k < nb_segment;k++) {
                if (mean[j][0][k] != 0.) {
                  for (m = change_point[k];m < change_point[k + 1];m++) {
                    real_sequence[i][j][m] /= mean[j][0][k];
                  }
                }
              }
            }
          }
        }
      }
    }

    else {
      for (i = 0;i < nb_sequence;i++) {
        for (j = 1;j < nb_variable;j++) {
          if (type[j] != REAL_VALUE) {
            real_sequence[i][j] = new double[length[i]];

            if ((model_type[j - 1] == LINEAR_MODEL_CHANGE) || (model_type[j - 1] == INTERCEPT_SLOPE_CHANGE)) {
              for (k = 0;k < nb_segment;k++) {
                for (m = change_point[k];m < change_point[k + 1];m++) {
                  if (intercept[j][i][k] + slope[j][i][k] * seq_index_parameter[m] != 0.) {
                    real_sequence[i][j][m] = int_sequence[i][j][m] / (intercept[j][i][k] + slope[j][i][k] * seq_index_parameter[m]);
                  }
                  else {
                    real_sequence[i][j][m] = int_sequence[i][j][m];
                  }
                }
              }
            }

            else if ((model_type[j - 1] == AUTOREGRESSIVE_MODEL_CHANGE) || (model_type[j - 1] == STATIONARY_AUTOREGRESSIVE_MODEL_CHANGE)) {
              for (k = 0;k < nb_segment;k++) {
                if (mean[j][i][k] != 0.) {
                  real_sequence[i][j][change_point[k]] = int_sequence[i][j][change_point[k]] / mean[j][i][k];
                }
                else {
                  real_sequence[i][j][change_point[k]] = int_sequence[i][j][change_point[k]];
                }
                for (m = change_point[k] + 1;m < change_point[k + 1];m++) {
                  if (mean[j][i][k] + autoregressive_coeff[j][i][k] * (int_sequence[i][j][m - 1] - mean[j][i][k]) != 0.) {
                    real_sequence[i][j][m] = int_sequence[i][j][m] / (mean[j][i][k] + autoregressive_coeff[j][i][k] *
                                             (int_sequence[i][j][m - 1] - mean[j][i][k]));
                  }
                  else {
                    real_sequence[i][j][m] = int_sequence[i][j][m];
                  }
                }
              }
            }

            else {
              for (k = 0;k < nb_segment;k++) {
                if (mean[j][i][k] != 0.) {
                  for (m = change_point[k];m < change_point[k + 1];m++) {
                    real_sequence[i][j][m] = int_sequence[i][j][m] / mean[j][i][k];
                  }
                }
                else {
                  for (m = change_point[k];m < change_point[k + 1];m++) {
                    real_sequence[i][j][m] = int_sequence[i][j][m];
                  }
                }
              }
            }

            delete [] int_sequence[i][j];
            int_sequence[i][j] = NULL;
          }

          else {
            if ((model_type[j - 1] == LINEAR_MODEL_CHANGE) || (model_type[j - 1] == INTERCEPT_SLOPE_CHANGE)) {
              for (k = 0;k < nb_segment;k++) {
                for (m = change_point[k];m < change_point[k + 1];m++) {
                  if (intercept[j][i][k] + slope[j][i][k] * seq_index_parameter[m] != 0.) {
                    real_sequence[i][j][m] /= (intercept[j][i][k] + slope[j][i][k] * seq_index_parameter[m]);
                  }
                }
              }
            }

            else if ((model_type[j - 1] == AUTOREGRESSIVE_MODEL_CHANGE) || (model_type[j - 1] == STATIONARY_AUTOREGRESSIVE_MODEL_CHANGE)) {
              for (k = 0;k < nb_segment;k++) {
                for (m = change_point[k + 1] - 1;m > change_point[k];m--) {
                  if (mean[j][i][k] + autoregressive_coeff[j][i][k] * (real_sequence[i][j][m - 1] - mean[j][i][k]) != 0.) {
                    real_sequence[i][j][m] /= (mean[j][i][k] + autoregressive_coeff[j][i][k] *
                                               (real_sequence[i][j][m - 1] - mean[j][i][k]));
                  }
                }
              }
              if (mean[j][i][change_point[k]] != 0.) {
                real_sequence[i][j][change_point[k]] /= mean[j][i][k];
              }
            }

            else {
              for (k = 0;k < nb_segment;k++) {
                if (mean[j][i][k] != 0.) {
                  for (m = change_point[k];m < change_point[k + 1];m++) {
                    real_sequence[i][j][m] /= mean[j][i][k];
                  }
                }
              }
            }
          }
        }
      }
    }
    break;
  }
  }

  if (!ichange_point) {
    delete [] change_point;
  }

  for (i = 1;i < nb_variable;i++) {
    if ((model_type[i - 1] == POISSON_CHANGE) || (model_type[i - 1] == NEGATIVE_BINOMIAL_0_CHANGE) ||
        (model_type[i - 1] == NEGATIVE_BINOMIAL_1_CHANGE) || (model_type[i - 1] == GAUSSIAN_CHANGE) ||
        (model_type[0] == MEAN_CHANGE) || (model_type[i - 1] == VARIANCE_CHANGE) ||
        (model_type[i - 1] == BAYESIAN_POISSON_CHANGE) || (model_type[i - 1] == BAYESIAN_GAUSSIAN_CHANGE)) {
      if (common_contrast) {
        delete [] mean[i][0];
        for (j = 0;j < 4;j++) {
          delete [] variance[i][j];
        }
      }
      else{
        for (j = 0;j < nb_sequence;j++) {
          delete [] mean[i][j];
          delete [] variance[i][j];
        }
      }

      delete [] mean[i];
      delete [] variance[i];
    }

    else if ((model_type[i - 1] == LINEAR_MODEL_CHANGE) || (model_type[0] == INTERCEPT_SLOPE_CHANGE)) {
     if (common_contrast) {
        delete [] intercept[i][0];
        delete [] slope[i][0];
        delete [] variance[i][0];
        delete [] correlation[i][0];
        delete [] slope_standard_deviation[i][0];
        delete [] index_parameter_mean[i][0];
        delete [] index_parameter_variance[i][0];

        if (continuity) {
          delete [] corrected_intercept[i][0];
          delete [] corrected_intercept[i];
          delete [] corrected_slope[i][0];
          delete [] corrected_slope[i];
        }
      }

      else {
        for (j = 0;j < nb_sequence;j++) {
          delete [] intercept[i][j];
          delete [] slope[i][j];
          delete [] variance[i][j];
          delete [] correlation[i][j];
          delete [] slope_standard_deviation[i][j];
          delete [] index_parameter_mean[i][j];
          delete [] index_parameter_variance[i][j];
        }
      }

      delete [] intercept[i];
      delete [] slope[i];
      delete [] variance[i];
      delete [] correlation[i];
      delete [] slope_standard_deviation[i];
      delete [] index_parameter_mean[i];
      delete [] index_parameter_variance[i];
    }

    else if ((model_type[i - 1] == AUTOREGRESSIVE_MODEL_CHANGE) || (model_type[i - 1] == STATIONARY_AUTOREGRESSIVE_MODEL_CHANGE)) {
      if (common_contrast) {
        delete [] mean[i][0];
        delete [] autoregressive_coeff[i][0];
        delete [] variance[i][0];
        delete [] determination_coeff[i][0];
      }

      else {
        for (j = 0;j < nb_sequence;j++) {
          delete [] mean[i][j];
          delete [] autoregressive_coeff[i][j];
          delete [] variance[i][j];
          delete [] determination_coeff[i][j];
        }
      }

      delete [] mean[i];
      delete [] autoregressive_coeff[i];
      delete [] variance[i];
      delete [] determination_coeff[i];
    }
  }

  delete [] mean;
  delete [] variance;
  delete [] global_variance;
  delete [] intercept;
  delete [] slope;
  delete [] correlation;
  delete [] slope_standard_deviation;
  delete [] index_parameter_mean;
  delete [] index_parameter_variance;
  delete [] autoregressive_coeff;
  delete [] determination_coeff;

  if (continuity) {
    delete [] corrected_intercept;
    delete [] corrected_slope;
  }

  if (index_param_type == IMPLICIT_TYPE) {
    delete [] seq_index_parameter;
  }

  if (output == SEQUENCE) {
    delete [] piecewise_function_flag;
  }

  if ((output == SUBTRACTION_RESIDUAL) || (output == DIVISION_RESIDUAL)) {
    for (i = 1;i < nb_variable;i++) {
      type[i] = REAL_VALUE;
    }

    seq = this;
  }

  if (output == SEQUENCE) {
    for (i = 1;i < seq->nb_variable;i++) {
      if (seq->type[i] == AUXILIARY) {
        seq->min_value_computation(i);
        seq->max_value_computation(i);
      }
    }
  }

  else if (output == ABSOLUTE_RESIDUAL) {
    for (i = 1;i < seq->nb_variable;i++) {
      seq->min_value_computation(i);
      seq->max_value_computation(i);

      if (seq->type[i] != AUXILIARY) {
        seq->build_marginal_histogram(i);
      }
    }
  }

  else {
    for (i = 1;i < seq->nb_variable;i++) {
      seq->min_value_computation(i);
      seq->max_value_computation(i);

      if (seq->marginal_distribution[i]) {
        delete seq->marginal_distribution[i];
        seq->marginal_distribution[i] = NULL;
      }

      seq->build_marginal_histogram(i);
    }
  }

  return seq;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Segmentation of a single sequence or a sample of sequences.
 *
 *  \param[in] error           reference on a StatError object,
 *  \param[in] display         flag for displaying the segmentation,
 *  \param[in] iidentifier     sequence identifier,
 *  \param[in] nb_segment      number of segments,
 *  \param[in] ichange_point   change points,
 *  \param[in] model_type      segment model types,
 *  \param[in] common_contrast flag contrast functions common to the individuals,
 *  \param[in] shape_parameter negative binomial shape parameters,
 *  \param[in] output          output (sequence or residuals),
 *  \param[in] continuity      flag continuous piecewise linear function.
 *
 *  \return                    Sequences object.
 */
/*--------------------------------------------------------------*/

Sequences* Sequences::segmentation(StatError &error , bool display , int iidentifier ,
                                   int nb_segment , int *ichange_point , segment_model *model_type ,
                                   bool common_contrast , double *shape_parameter ,
                                   sequence_type output , bool continuity) const

{
  bool status = true;
  int i , j , k , m;
  int index , segmentation_index , seq_length , count , max_nb_value , nb_parameter ,
      *change_point = NULL , *inf_bound_parameter , *frequency , *seq_index_parameter;
  double sum , factorial_sum , binomial_coeff_sum , proba , mean , diff , index_parameter_mean ,
         index_parameter_diff , index_parameter_sum , shifted_diff , segmentation_likelihood ,
         segment_penalty , penalized_likelihood , **rank , **seq_mean;
  long double index_parameter_square_sum , square_sum , mix_square_sum , shifted_square_sum ,
              autocovariance , **residual;
  FrequencyDistribution *marginal;
  Sequences *iseq , *seq , *oseq;


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
        (model_type[i] == ORDINAL_GAUSSIAN_CHANGE)) {
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
        if (((model_type[i] != NEGATIVE_BINOMIAL_0_CHANGE) && (min_value[i] < 0)) ||
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

  if (status) {
    segmentation_index = (index == I_DEFAULT ? 0 : index);
    seq_length = length[segmentation_index];

    if ((nb_segment < 1) || (nb_segment > seq_length / 2)) {
      status = false;
      error.update(SEQ_error[SEQR_NB_SEGMENT]);
    }

    else {
      change_point = new int[nb_segment + 1];

      if (index_parameter) {
        change_point[0] = index_parameter[segmentation_index][0];
        change_point[nb_segment] = index_parameter[segmentation_index][seq_length - 1] + 1;
      }
      else {
        change_point[0] = 0;
        change_point[nb_segment] = seq_length;
      }

      for (i = 1;i < nb_segment;i++) {
        change_point[i] = ichange_point[i - 1];
      }

      for (i = 1;i < nb_segment - 1;i++) {
        if (change_point[i] >= change_point[i + 1]) {
          status = false;
          error.update(SEQ_error[SEQR_CHANGE_POINT]);
        }
      }

      if (index_parameter) {
        change_point[0] = 0;
        i = 1;
        for (j = 1;j < seq_length;j++) {
          if (index_parameter[segmentation_index][j] == change_point[i]) {
            change_point[i++] = j;
          }
        }

        if (i < nb_segment) {
          status = false;
          error.update(SEQ_error[SEQR_CHANGE_POINT]);
        }
        else {
          change_point[nb_segment] = seq_length;
        }
      }
    }
  }

  if (status) {
    max_nb_value = 0;
    inf_bound_parameter = new int[nb_variable];
    seq_mean = new double*[nb_variable];
    seq_index_parameter = NULL;
    rank = new double*[nb_variable];
    residual = NULL;

    for (i = 0;i < nb_variable;i++) {
      if ((model_type[i] == CATEGORICAL_CHANGE) && (marginal_distribution[i]->nb_value > max_nb_value)) {
        max_nb_value = marginal_distribution[i]->nb_value;
      }

      if ((model_type[i] == NEGATIVE_BINOMIAL_0_CHANGE) || (model_type[i] == NEGATIVE_BINOMIAL_1_CHANGE)) {
        switch (model_type[i]) {
        case NEGATIVE_BINOMIAL_0_CHANGE :
          inf_bound_parameter[i] = 0;
          break;
        case NEGATIVE_BINOMIAL_1_CHANGE :
          inf_bound_parameter[i] = 1;
          break;
        }
      }

      // computation of sequence means for Gaussian change in the variance model or
      // stationary piecewise autoregressive models

      if ((model_type[i] == VARIANCE_CHANGE) || (model_type[i] == STATIONARY_AUTOREGRESSIVE_MODEL_CHANGE)) {
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

      // rank computation for ordinal variables

      if (model_type[i] == ORDINAL_GAUSSIAN_CHANGE) {
        rank[i] = marginal_distribution[i]->rank_computation();
      }
      else {
        rank[i] = NULL;
      }

      if (((i == 0) && (model_type[0] == INTERCEPT_SLOPE_CHANGE)) ||
          ((model_type[i] == LINEAR_MODEL_CHANGE) && (!seq_index_parameter))) {
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

      if (((i == 0) && ((model_type[0] == MEAN_CHANGE) || (model_type[0] == INTERCEPT_SLOPE_CHANGE))) ||
          (((model_type[i] == GAUSSIAN_CHANGE) || (model_type[i] == VARIANCE_CHANGE) ||
            (model_type[i] == ORDINAL_GAUSSIAN_CHANGE) || (model_type[i] == LINEAR_MODEL_CHANGE) ||
            (model_type[i] == AUTOREGRESSIVE_MODEL_CHANGE) ||
            (model_type[i] == STATIONARY_AUTOREGRESSIVE_MODEL_CHANGE)) && (!residual))) {
        residual = new long double*[nb_sequence];
        if ((index != I_DEFAULT) || (!common_contrast)) {
          for (j = 0;j < nb_sequence;j++) {
            if ((index == I_DEFAULT) || (index == j)) {
              residual[j] = new long double[nb_segment];
            }
            else {
              residual[j] = NULL;
            }
          }
        }
        else {
          residual[0] = new long double[nb_segment];
        }
      }
    }

    if (max_nb_value > 0) {
      frequency = new int[max_nb_value];
    }
    else {
      frequency = NULL;
    }

    if ((model_type[0] != MEAN_CHANGE) && (model_type[0] != INTERCEPT_SLOPE_CHANGE)) {
      segmentation_likelihood = 0.;
    }

    for (i = 0;i < nb_variable;i++) {

      // computation of segment log-likelihoods

      if (model_type[i] == CATEGORICAL_CHANGE) {
        if ((index != I_DEFAULT) || (!common_contrast)) {
          for (j = 0;j < nb_sequence;j++) {
            if ((index == I_DEFAULT) || (index == j)) {
              for (k = 0;k < nb_segment;k++) {
                for (m = 0;m < marginal_distribution[i]->nb_value;m++) {
                  frequency[m] = 0;
                }

                for (m = change_point[k];m < change_point[k + 1];m++) {
                  frequency[int_sequence[j][i][m]]++;
                }

                for (m = 0;m < marginal_distribution[i]->nb_value;m++) {
                  if (frequency[m] > 0) {
                    segmentation_likelihood += frequency[m] * log((double)frequency[m] /
                                                (double)(change_point[k + 1] - change_point[k]));
                  }
                }
              }
            }
          }
        }

        else {
          for (j = 0;j < nb_segment;j++) {
            for (k = 0;k < marginal_distribution[i]->nb_value;k++) {
              frequency[k] = 0;
            }

            for (k = change_point[j];k < change_point[j + 1];k++) {
              for (m = 0;m < nb_sequence;m++) {
                frequency[int_sequence[m][i][k]]++;
              }
            }

            for (k = 0;k < marginal_distribution[i]->nb_value;k++) {
              if (frequency[k] > 0) {
                segmentation_likelihood += frequency[k] * log((double)frequency[k] /
                                            (double)(nb_sequence * (change_point[j + 1] - change_point[j])));
              }
            }
          }
        }
      }

      else if (model_type[i] == POISSON_CHANGE) {
        if ((index != I_DEFAULT) || (!common_contrast)) {
          for (j = 0;j < nb_sequence;j++) {
            if ((index == I_DEFAULT) || (index == j)) {
              for (k = 0;k < nb_segment;k++) {
                sum = 0.;
                factorial_sum = 0.;
                for (m = change_point[k];m < change_point[k + 1];m++) {
                  sum += int_sequence[j][i][m];
                  factorial_sum += log_factorial(int_sequence[j][i][m]);
                }

                if (sum > 0.) {
                  segmentation_likelihood += sum * (log(sum / (change_point[k + 1] - change_point[k])) - 1) -
                                             factorial_sum;
                }
              }
            }
          }
        }

        else {
          for (j = 0;j < nb_segment;j++) {
            sum = 0.;
            factorial_sum = 0.;
            for (k = change_point[j];k < change_point[j + 1];k++) {
              for (m = 0;m < nb_sequence;m++) {
                sum += int_sequence[m][i][k];
                factorial_sum += log_factorial(int_sequence[m][i][k]);
              }
            }

            if (sum > 0.) {
              segmentation_likelihood += sum * (log(sum / (nb_sequence * (change_point[j + 1] - change_point[j]))) - 1) -
                                         factorial_sum;
            }
          }
        }
      }

      else if ((model_type[i] == NEGATIVE_BINOMIAL_0_CHANGE) || (model_type[i] == NEGATIVE_BINOMIAL_1_CHANGE)) {
        if ((index != I_DEFAULT) || (!common_contrast)) {
          for (j = 0;j < nb_sequence;j++) {
            if ((index == I_DEFAULT) || (index == j)) {
              for (k = 0;k < nb_segment;k++) {
                sum = 0.;
                binomial_coeff_sum = 0.;
                for (m = change_point[k];m < change_point[k + 1];m++) {
                  sum += int_sequence[j][i][m];
                  binomial_coeff_sum += log_binomial_coefficient(inf_bound_parameter[i] , shape_parameter[i] ,
                                                                 int_sequence[j][i][m]);
                }

                if (sum > inf_bound_parameter[i] * (change_point[k + 1] - change_point[k])) {
                  proba = shape_parameter[i] * (change_point[k + 1] - change_point[k]) /
                          ((shape_parameter[i] - inf_bound_parameter[i]) * (change_point[k + 1] - change_point[k]) + sum);
                  segmentation_likelihood += binomial_coeff_sum + shape_parameter[i] * (change_point[k + 1] - change_point[k]) * log(proba) +
                                             (sum - inf_bound_parameter[i] *  (change_point[k + 1] - change_point[k])) * log(1. - proba);
                }
                else {
                  segmentation_likelihood = D_INF;
                  break;
                }
              }
            }

            if (segmentation_likelihood == D_INF) {
              break;
            }
          }
        }

        else {
          for (j = 0;j < nb_segment;j++) {
            sum = 0.;
            binomial_coeff_sum = 0.;
            for (k = change_point[j];k < change_point[j + 1];k++) {
              for (m = 0;m < nb_sequence;m++) {
                sum += int_sequence[m][i][k];
                binomial_coeff_sum += log_binomial_coefficient(inf_bound_parameter[i] , shape_parameter[i] ,
                                                               int_sequence[m][i][k]);
              }
            }

            if (sum > inf_bound_parameter[i] * nb_sequence * (change_point[j + 1] - change_point[j])) {
              proba = shape_parameter[i] * nb_sequence * (change_point[j + 1] - change_point[j]) /
                      ((shape_parameter[i] - inf_bound_parameter[i]) * nb_sequence * (change_point[j + 1] - change_point[j]) + sum);
              segmentation_likelihood += binomial_coeff_sum + shape_parameter[i] * nb_sequence * (change_point[j + 1] - change_point[j]) * log(proba) +
                                         (sum - inf_bound_parameter[i] * nb_sequence * (change_point[j + 1] - change_point[j])) * log(1. - proba);
            }
            else {
              segmentation_likelihood = D_INF;
              break;
            }
          }
        }
      }

      else if ((model_type[i] == GAUSSIAN_CHANGE)|| (model_type[0] == MEAN_CHANGE)) {
        if ((index != I_DEFAULT) || (!common_contrast)) {
          if (type[i] != REAL_VALUE) {
            for (j = 0;j < nb_sequence;j++) {
              if ((index == I_DEFAULT) || (index == j)) {
                for (k = 0;k < nb_segment;k++) {
/*                  residual[j][k] = 0.;
                  sum = int_sequence[j][i][change_point[k]];

                  for (m = change_point[k] + 1;m < change_point[k + 1];m++) {
                    diff = int_sequence[j][i][m] - sum / (m - change_point[k]);
                    residual[j][k] += ((double)(m - change_point[k]) / (double)(m - change_point[k] + 1)) *
                                      diff * diff;
                    sum += int_sequence[j][i][m];
                  } */

                  mean = 0.;
                  for (m = change_point[k];m < change_point[k + 1];m++) {
                    mean += int_sequence[j][i][m];
                  }
                  mean /= (change_point[k + 1] - change_point[k]);

                  residual[j][k] = 0.;
                  for (m = change_point[k];m < change_point[k + 1];m++) {
                    diff = int_sequence[j][i][m] - mean;
                    residual[j][k] += diff * diff;
                  }
                }
              }
            }
          }

          else {
            for (j = 0;j < nb_sequence;j++) {
              if ((index == I_DEFAULT) || (index == j)) {
                for (k = 0;k < nb_segment;k++) {
/*                  residual[j][k] = 0.;
                  sum = real_sequence[j][i][change_point[k]];

                  for (m = change_point[k] + 1;m < change_point[k + 1];m++) {
                    diff = real_sequence[j][i][m] - sum / (m - change_point[k]);
                    residual[j][k] += ((double)(m - change_point[k]) / (double)(m - change_point[k] + 1)) *
                                      diff * diff;
                    sum += real_sequence[j][i][m];
                  } */

                  mean = 0.;
                  for (m = change_point[k];m < change_point[k + 1];m++) {
                    mean += real_sequence[j][i][m];
                  }
                  mean /= (change_point[k + 1] - change_point[k]);

                  residual[j][k] = 0.;
                  for (m = change_point[k];m < change_point[k + 1];m++) {
                    diff = real_sequence[j][i][m] - mean;
                    residual[j][k] += diff * diff;
                  }
                }
              }
            }
          }
        }

        else {
          if (type[i] != REAL_VALUE) {
            for (j = 0;j < nb_segment;j++) {
/*              residual[0][j] = 0.;
              sum = 0.;
              count = 0;

              for (k = change_point[j];k < change_point[j + 1];k++) {
                for (m = 0;m < nb_sequence;m++) {
                  if (count > 0) {
                    diff = int_sequence[m][i][k] - sum / count;
                    residual[0][j] += ((double)count / (double)(count + 1)) * diff * diff;
                  }
                  count++;
                  sum += int_sequence[m][i][k];
                }
              } */

              mean = 0.;
              for (k = change_point[j];k < change_point[j + 1];k++) {
                for (m = 0;m < nb_sequence;m++) {
                  mean += int_sequence[m][i][k];
                }
              }
              mean /= nb_sequence * (change_point[j + 1] -  change_point[j]);

              residual[0][j] = 0.;
              for (k = change_point[j];k < change_point[j + 1];k++) {
                for (m = 0;m < nb_sequence;m++) {
                  diff = int_sequence[m][i][k] - mean;
                  residual[0][j] += diff * diff;
                }
              }
            }
          }

          else {
            for (j = 0;j < nb_segment;j++) {
/*              residual[0][j] = 0.;
              sum = 0.;
              count = 0;

              for (k = change_point[j];k < change_point[j + 1];k++) {
                for (m = 0;m < nb_sequence;m++) {
                  if (count > 0) {
                    diff = real_sequence[m][i][k] - sum / count;
                    residual[0][j] += ((double)count / (double)(count + 1)) * diff * diff;
                  }
                  count++;
                  sum += real_sequence[m][i][k];
                }
              } */

              mean = 0.;
              for (k = change_point[j];k < change_point[j + 1];k++) {
                for (m = 0;m < nb_sequence;m++) {
                  mean += real_sequence[m][i][k];
                }
              }
              mean /= nb_sequence * (change_point[j + 1] -  change_point[j]);

              residual[0][j] = 0.;
              for (k = change_point[j];k < change_point[j + 1];k++) {
                for (m = 0;m < nb_sequence;m++) {
                   diff = real_sequence[m][i][k] - mean;
                   residual[0][j] += diff * diff;
                }
              }
            }
          }
        }
      }

      else if (model_type[i] == VARIANCE_CHANGE) {
        if ((index != I_DEFAULT) || (!common_contrast)) {
          if (type[i] != REAL_VALUE) {
            for (j = 0;j < nb_sequence;j++) {
              if ((index == I_DEFAULT) || (index == j)) {
                for (k = 0;k < nb_segment;k++) {
                  residual[j][k] = 0.;

                  for (m = change_point[k];m < change_point[k + 1];m++) {
                    diff = int_sequence[j][i][m] - seq_mean[i][j];
                    residual[j][k] += diff * diff;
                  }
                }
              }
            }
          }

          else {
            for (j = 0;j < nb_sequence;j++) {
              if ((index == I_DEFAULT) || (index == j)) {
                for (k = 0;k < nb_segment;k++) {
                  residual[j][k] = 0.;

                  for (m = change_point[k];m < change_point[k + 1];m++) {
                    diff = real_sequence[j][i][m] - seq_mean[i][j];
                    residual[j][k] += diff * diff;
                  }
                }
              }
            }
          }
        }

        else {
          if (type[i] != REAL_VALUE) {
            for (j = 0;j < nb_segment;j++) {
              residual[0][j] = 0.;

              for (k = change_point[j];k < change_point[j + 1];k++) {
                for (m = 0;m < nb_sequence;m++) {
                  diff = int_sequence[m][i][k] - seq_mean[i][0];
                  residual[0][j] += diff * diff;
                }
              }
            }
          }

          else {
            for (j = 0;j < nb_segment;j++) {
              residual[0][j] = 0.;

              for (k = change_point[j];k < change_point[j + 1];k++) {
                for (m = 0;m < nb_sequence;m++) {
                  diff = real_sequence[m][i][k] - seq_mean[i][0];
                  residual[0][j] += diff * diff;
                }
              }
            }
          }
        }
      }

      else if (model_type[i] == ORDINAL_GAUSSIAN_CHANGE) {
        if ((index != I_DEFAULT) || (!common_contrast)) {
          for (j = 0;j < nb_sequence;j++) {
            if ((index == I_DEFAULT) || (index == j)) {
              for (k = 0;k < nb_segment;k++) {
/*                residual[j][k] = 0.;
                sum = rank[i][int_sequence[j][i][change_point[k]]];

                for (m = change_point[k] + 1;m < change_point[k + 1];m++) {
                  diff = rank[i][int_sequence[j][i][m]] - sum / (m - change_point[k]);
                  residual[j][k] += ((double)(m - change_point[k]) / (double)(m - change_point[k] + 1)) *
                                    diff * diff;
                  sum += rank[i][int_sequence[j][i][m]];
                } */

                mean = 0.;
                for (m = change_point[k];m < change_point[k + 1];m++) {
                  mean += rank[i][int_sequence[j][i][m]];
                }
                mean /= (change_point[k + 1] - change_point[k]);

                residual[j][k] = 0.;
                for (m = change_point[k];m < change_point[k + 1];m++) {
                  diff = rank[i][int_sequence[j][i][m]] - mean;
                  residual[j][k] += diff * diff;
                }

                if (residual[j][k] == 0.) {
                  residual[j][k] = (change_point[k + 1] - change_point[k]) * MIN_RANK_SQUARE_SUM;
                }
              }
            }
          }
        }

        else {
          for (j = 0;j < nb_segment;j++) {
/*            residual[0][j] = 0.;
            sum = 0.;
            count = 0;

            for (k = change_point[j];k < change_point[j + 1];k++) {
              for (m = 0;m < nb_sequence;m++) {
                if (count > 0) {
                  diff = rank[i][int_sequence[m][i][k]] - sum / count;
                  residual[0][j] += ((double)count / (double)(count + 1)) * diff * diff;
                }
                count++;
                sum += rank[i][int_sequence[m][i][k]];
              }
            } */

            mean = 0.;
            for (k = change_point[j];k < change_point[j + 1];k++) {
              for (m = 0;m < nb_sequence;m++) {
                mean += rank[i][int_sequence[m][i][k]];
              }
            }
            mean /= nb_sequence * (change_point[j + 1] -  change_point[j]);

            residual[0][j] = 0.;
            for (k = change_point[j];k < change_point[j + 1];k++) {
              for (m = 0;m < nb_sequence;m++) {
                diff = rank[i][int_sequence[m][i][k]] - mean;
                residual[0][j] += diff * diff;
              }
            }

            if (residual[0][j] == 0.) {
              residual[0][j] = nb_sequence * (change_point[k + 1] - change_point[k]) * MIN_RANK_SQUARE_SUM;
            }
          }
        }
      }

      else if ((model_type[i] == LINEAR_MODEL_CHANGE) || (model_type[0] == INTERCEPT_SLOPE_CHANGE)) {
        if ((index != I_DEFAULT) || (!common_contrast)) {
          if (type[i] != REAL_VALUE) {
            for (j = 0;j < nb_sequence;j++) {
              if ((index == I_DEFAULT) || (index == j)) {
                for (k = 0;k < nb_segment;k++) {
/*                  index_parameter_square_sum = 0.;
                  square_sum = 0.;
                  mix_square_sum = 0.;
                  index_parameter_sum = seq_index_parameter[change_point[k]];
                  sum = int_sequence[j][i][change_point[k]];

                  for (m = change_point[k] + 1;m < change_point[k + 1];m++) {
                    index_parameter_diff = seq_index_parameter[m] - index_parameter_sum / (m - change_point[k]);
                    index_parameter_square_sum += ((double)(m - change_point[k]) / (double)(m - change_point[k] + 1)) *
                                                  index_parameter_diff * index_parameter_diff;
                    diff = int_sequence[j][i][m] - sum / (m - change_point[k]);
                    square_sum += ((double)(m - change_point[k]) / (double)(m - change_point[k] + 1)) *
                                     diff * diff;
                    mix_square_sum += ((double)(m - change_point[k]) / (double)(m - change_point[k] + 1)) *
                                      index_parameter_diff * diff;
                    index_parameter_sum += seq_index_parameter[m];
                    sum += int_sequence[j][i][m];
                  } */

                  index_parameter_mean = 0.;
                  mean = 0.;
                  for (m = change_point[k];m < change_point[k + 1];m++) {
                    index_parameter_mean += seq_index_parameter[m];
                    mean += int_sequence[j][i][m];
                  }
                  index_parameter_mean /= (change_point[k + 1] - change_point[k]);
                  mean /= (change_point[k + 1] - change_point[k]);

                  index_parameter_square_sum = 0.;
                  square_sum = 0.;
                  mix_square_sum = 0.;
                  for (m = change_point[k];m < change_point[k + 1];m++) {
                    index_parameter_diff = seq_index_parameter[m] - index_parameter_mean;
                    diff = int_sequence[j][i][m] - mean;
                    index_parameter_square_sum += index_parameter_diff * index_parameter_diff;
                    square_sum += diff * diff;
                    mix_square_sum += index_parameter_diff * diff;
                  }

                  if ((change_point[k + 1] - change_point[k] > 2) && (index_parameter_square_sum > 0.)) {
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
                for (k = 0;k < nb_segment;k++) {
/*                  index_parameter_square_sum = 0.;
                  square_sum = 0.;
                  mix_square_sum = 0.;
                  index_parameter_sum = seq_index_parameter[change_point[k]];
                  sum = real_sequence[j][i][change_point[k]];

                  for (m = change_point[k] + 1;m < change_point[k + 1];m++) {
                    index_parameter_diff = seq_index_parameter[m] - index_parameter_sum / (m - change_point[k]);
                    index_parameter_square_sum += ((double)(m - change_point[k]) / (double)(m - change_point[k] + 1)) *
                                                  index_parameter_diff * index_parameter_diff;
                    diff = real_sequence[j][i][m] - sum / (m - change_point[k]);
                    square_sum += ((double)(m - change_point[k]) / (double)(m - change_point[k] + 1)) *
                                     diff * diff;
                    mix_square_sum += ((double)(m - change_point[k]) / (double)(m - change_point[k] + 1)) *
                                      index_parameter_diff * diff;
                    index_parameter_sum += seq_index_parameter[m];
                    sum += real_sequence[j][i][m];
                  } */

                  index_parameter_mean = 0.;
                  mean = 0.;
                  for (m = change_point[k];m < change_point[k + 1];m++) {
                    index_parameter_mean += seq_index_parameter[m];
                    mean += real_sequence[j][i][m];
                  }
                  index_parameter_mean /= (change_point[k + 1] - change_point[k]);
                  mean /= (change_point[k + 1] - change_point[k]);

                  index_parameter_square_sum = 0.;
                  square_sum = 0.;
                  mix_square_sum = 0.;
                  for (m = change_point[k];m < change_point[k + 1];m++) {
                    index_parameter_diff = seq_index_parameter[m] - index_parameter_mean;
                    diff = real_sequence[j][i][m] - mean;
                    index_parameter_square_sum += index_parameter_diff * index_parameter_diff;
                    square_sum += diff * diff;
                    mix_square_sum += index_parameter_diff * diff;
                  }

                  if ((change_point[k + 1] - change_point[k] > 2) && (index_parameter_square_sum > 0.)) {
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
          if (type[i] != REAL_VALUE) {
            for (j = 0;j < nb_segment;j++) {
/*              index_parameter_square_sum = 0.;
              square_sum = 0.;
              mix_square_sum = 0.;
              count = 1;

              index_parameter_sum = nb_sequence * seq_index_parameter[change_point[j]];
              sum = int_sequence[0][i][change_point[j]];
              for (k = 1;k < nb_sequence;k++) {
                diff = int_sequence[k][i][change_point[j]] - sum / count;
                square_sum += ((double)count / (double)(count + 1)) * diff * diff;
                count++;
                sum += int_sequence[k][i][change_point[j]];
              }

              for (k = change_point[j] + 1;k < change_point[j + 1];k++) {
                for (m = 0;m < nb_sequence;m++) {
                  index_parameter_diff = seq_index_parameter[k] - index_parameter_sum / count;
                  index_parameter_square_sum += ((double)count / (double)(count + 1)) *
                                                index_parameter_diff * index_parameter_diff;
                  diff = int_sequence[m][i][k] - sum / count;
                  square_sum += ((double)count / (double)(count + 1)) * diff * diff;
                  mix_square_sum += ((double)count / (double)(count + 1)) * index_parameter_diff * diff;
                  count++;
                  index_parameter_sum += seq_index_parameter[k];
                  sum += int_sequence[m][i][k];
                }
              } */

              index_parameter_mean = 0.;
              mean = 0.;
              for (k = change_point[j];k < change_point[j + 1];k++) {
                index_parameter_mean += seq_index_parameter[k];
                for (m = 0;m < nb_sequence;m++) {
                  mean += int_sequence[m][i][k];
                }
              }
              index_parameter_mean /= (change_point[j + 1] - change_point[j]);
              mean /= nb_sequence * (change_point[j + 1] - change_point[j]);

              index_parameter_square_sum = 0.;
              square_sum = 0.;
              mix_square_sum = 0.;
              for (k = change_point[j];k < change_point[j + 1];k++) {
                index_parameter_diff = seq_index_parameter[k] - index_parameter_mean;
                index_parameter_square_sum += index_parameter_diff * index_parameter_diff;
                for (m = 0;m < nb_sequence;m++) {
                  diff = int_sequence[m][i][k] - mean;
                  square_sum += diff * diff;
                  mix_square_sum += index_parameter_diff * diff;
                }
              }
              index_parameter_square_sum *= nb_sequence;

              if (index_parameter_square_sum > 0.) {
                residual[0][j] = square_sum - mix_square_sum * mix_square_sum / index_parameter_square_sum;
              }
              else {
                residual[0][j] = 0.;
              }
            }
          }

          else {
            for (j = 0;j < nb_segment;j++) {
/*              index_parameter_square_sum = 0.;
              square_sum = 0.;
              mix_square_sum = 0.;
              count = 1;

              index_parameter_sum = nb_sequence * seq_index_parameter[change_point[j]];
              sum = real_sequence[0][i][change_point[j]];
              for (k = 1;k < nb_sequence;k++) {
                diff = real_sequence[k][i][change_point[j]] - sum / count;
                square_sum += ((double)count / (double)(count + 1)) * diff * diff;
                count++;
                sum += real_sequence[k][i][change_point[j]];
              }

              for (k = change_point[j] + 1;k < change_point[j + 1];k++) {
                for (m = 0;m < nb_sequence;m++) {
                  index_parameter_diff = seq_index_parameter[k] - index_parameter_sum / count;
                  index_parameter_square_sum += ((double)count / (double)(count + 1)) *
                                                index_parameter_diff * index_parameter_diff;
                  diff = real_sequence[m][i][k] - sum / count;
                  square_sum += ((double)count / (double)(count + 1)) * diff * diff;
                  mix_square_sum += ((double)count / (double)(count + 1)) * index_parameter_diff * diff;
                  count++;
                  index_parameter_sum += seq_index_parameter[k];
                  sum += real_sequence[m][i][k];
                }
              } */

              index_parameter_mean = 0.;
              mean = 0.;
              for (k = change_point[j];k < change_point[j + 1];k++) {
                index_parameter_mean += seq_index_parameter[k];
                for (m = 0;m < nb_sequence;m++) {
                  mean += real_sequence[m][i][k];
                }
              }
              index_parameter_mean /= (change_point[j + 1] - change_point[j]);
              mean /= nb_sequence * (change_point[j + 1] - change_point[j]);

              index_parameter_square_sum = 0.;
              square_sum = 0.;
              mix_square_sum = 0.;
              for (k = change_point[j];k < change_point[j + 1];k++) {
                index_parameter_diff = seq_index_parameter[k] - index_parameter_mean;
                index_parameter_square_sum += index_parameter_diff * index_parameter_diff;
                for (m = 0;m < nb_sequence;m++) {
                  diff = real_sequence[m][i][k] - mean;
                  square_sum += diff * diff;
                  mix_square_sum += index_parameter_diff * diff;
                }
              }
              index_parameter_square_sum *= nb_sequence;

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

      else if (model_type[i] == AUTOREGRESSIVE_MODEL_CHANGE) {
        if ((index != I_DEFAULT) || (!common_contrast)) {
          if (type[i] != REAL_VALUE) {
            for (j = 0;j < nb_sequence;j++) {
              if ((index == I_DEFAULT) || (index == j)) {
                for (k = 0;k < nb_segment;k++) {
                  mean = 0.;
                  for (m = change_point[k];m < change_point[k + 1];m++) {
                    mean += int_sequence[j][i][m];
                  }
                  mean /= (change_point[k + 1] - change_point[k]);

                  square_sum = 0.;
                  shifted_square_sum = 0.;
                  autocovariance = 0.;
                  for (m = change_point[k] + 1;m < change_point[k + 1];m++) {
                    diff = int_sequence[j][i][m] - mean;
                    shifted_diff = int_sequence[j][i][m - 1] - mean;
                    square_sum += diff * diff;
                    shifted_square_sum += shifted_diff * shifted_diff;
                    autocovariance += diff * shifted_diff;
                  }

                  if (change_point[k + 1] - change_point[k] > 2) {
                    residual[j][k] = square_sum; 
                    if (shifted_square_sum > 0.) {
                      residual[j][k] -= autocovariance * autocovariance / shifted_square_sum; 
                    }
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
                for (k = 0;k < nb_segment;k++) {
                  mean = 0.;
                  for (m = change_point[k];m < change_point[k + 1];m++) {
                    mean += real_sequence[j][i][m];
                  }
                  mean /= (change_point[k + 1] - change_point[k]);

                  square_sum = 0.;
                  shifted_square_sum = 0.;
                  autocovariance = 0.;
                  for (m = change_point[k] + 1;m < change_point[k + 1];m++) {
                    diff = real_sequence[j][i][m] - mean;
                    shifted_diff = real_sequence[j][i][m - 1] - mean;
                    square_sum += diff * diff;
                    shifted_square_sum += shifted_diff * shifted_diff;
                    autocovariance += diff * shifted_diff;
                  }

                  if (change_point[k + 1] - change_point[k] > 2) {
                    residual[j][k] = square_sum; 
                    if (shifted_square_sum > 0.) {
                      residual[j][k] -= autocovariance * autocovariance / shifted_square_sum; 
                    }
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
          if (type[i] != REAL_VALUE) {
            for (j = 0;j < nb_segment;j++) {
              mean = 0.;
              for (k = change_point[j];k < change_point[j + 1];k++) {
                for (m = 0;m < nb_sequence;m++) {
                  mean += int_sequence[m][i][k];
                }
              }
              mean /= nb_sequence * (change_point[j + 1] - change_point[j]);

              square_sum = 0.;
              shifted_square_sum = 0.;
              autocovariance = 0.;
              for (k = change_point[j] + 1;k < change_point[j + 1];k++) {
                for (m = 0;m < nb_sequence;m++) {
                  diff = int_sequence[m][i][k] - mean;
                  shifted_diff = int_sequence[m][i][k - 1] - mean;
                  square_sum += diff * diff;
                  shifted_square_sum += shifted_diff * shifted_diff;
                  autocovariance += diff * shifted_diff;
                }
              }

              residual[0][j] = square_sum;
              if (shifted_square_sum > 0.) {
                residual[0][j] -= autocovariance * autocovariance / shifted_square_sum;
              }
            }
          }

          else {
            for (j = 0;j < nb_segment;j++) {
              mean = 0.;
              for (k = change_point[j];k < change_point[j + 1];k++) {
                for (m = 0;m < nb_sequence;m++) {
                  mean += real_sequence[m][i][k];
                }
              }
              mean /= nb_sequence * (change_point[j + 1] - change_point[j]);

              square_sum = 0.;
              shifted_square_sum = 0.;
              autocovariance = 0.;
              for (k = change_point[j] + 1;k < change_point[j + 1];k++) {
                for (m = 0;m < nb_sequence;m++) {
                  diff = real_sequence[m][i][k] - mean;
                  shifted_diff = real_sequence[m][i][k - 1] - mean;
                  square_sum += diff * diff;
                  shifted_square_sum += shifted_diff * shifted_diff;
                  autocovariance += diff * shifted_diff;
                }
              }

              residual[0][j] = square_sum;
              if (shifted_square_sum > 0.) {
                residual[0][j] -= autocovariance * autocovariance / shifted_square_sum;
              }
            }
          }
        }
      }


      else if (model_type[i] == STATIONARY_AUTOREGRESSIVE_MODEL_CHANGE) {
        if ((index != I_DEFAULT) || (!common_contrast)) {
          if (type[i] != REAL_VALUE) {
            for (j = 0;j < nb_sequence;j++) {
              if ((index == I_DEFAULT) || (index == j)) {
                for (k = 0;k < nb_segment;k++) {
                  square_sum = 0.;
                  shifted_square_sum = 0.;
                  autocovariance = 0.;
                  for (m = change_point[k] + 1;m < change_point[k + 1];m++) {
                    diff = int_sequence[j][i][m] - seq_mean[i][j];
                    shifted_diff = int_sequence[j][i][m - 1] - seq_mean[i][j];
                    square_sum += diff * diff;
                    shifted_square_sum += shifted_diff * shifted_diff;
                    autocovariance += diff * shifted_diff;
                  }

                  if (change_point[k + 1] - change_point[k] > 2) {
                    residual[j][k] = square_sum; 
                    if (shifted_square_sum > 0.) {
                      residual[j][k] -= autocovariance * autocovariance / shifted_square_sum; 
                    }
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
                for (k = 0;k < nb_segment;k++) {
                  square_sum = 0.;
                  shifted_square_sum = 0.;
                  autocovariance = 0.;
                  for (m = change_point[k] + 1;m < change_point[k + 1];m++) {
                    diff = real_sequence[j][i][m] - seq_mean[i][j];
                    shifted_diff = real_sequence[j][i][m - 1] - seq_mean[i][j];
                    square_sum += diff * diff;
                    shifted_square_sum += shifted_diff * shifted_diff;
                    autocovariance += diff * shifted_diff;
                  }

                  if (change_point[k + 1] - change_point[k] > 2) {
                    residual[j][k] = square_sum; 
                    if (shifted_square_sum > 0.) {
                      residual[j][k] -= autocovariance * autocovariance / shifted_square_sum; 
                    }
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
          if (type[i] != REAL_VALUE) {
            for (j = 0;j < nb_segment;j++) {
              square_sum = 0.;
              shifted_square_sum = 0.;
              autocovariance = 0.;
              for (k = change_point[j] + 1;k < change_point[j + 1];k++) {
                for (m = 0;m < nb_sequence;m++) {
                  diff = int_sequence[m][i][k] - seq_mean[i][0];
                  shifted_diff = int_sequence[m][i][k - 1] - seq_mean[i][0];
                  square_sum += diff * diff;
                  shifted_square_sum += shifted_diff * shifted_diff;
                  autocovariance += diff * shifted_diff;
                }
              }

              residual[0][j] = square_sum;
              if (shifted_square_sum > 0.) {
                residual[0][j] -= autocovariance * autocovariance / shifted_square_sum;
              }
            }
          }

          else {
            for (j = 0;j < nb_segment;j++) {
              square_sum = 0.;
              shifted_square_sum = 0.;
              autocovariance = 0.;
              for (k = change_point[j] + 1;k < change_point[j + 1];k++) {
                for (m = 0;m < nb_sequence;m++) {
                  diff = real_sequence[m][i][k] - seq_mean[i][0];
                  shifted_diff = real_sequence[m][i][k - 1] - seq_mean[i][0];
                  square_sum += diff * diff;
                  shifted_square_sum += shifted_diff * shifted_diff;
                  autocovariance += diff * shifted_diff;
                }
              }

              residual[0][j] = square_sum;
              if (shifted_square_sum > 0.) {
                residual[0][j] -= autocovariance * autocovariance / shifted_square_sum;
              }
            }
          }
        }
      }

      if ((model_type[i] == GAUSSIAN_CHANGE) || (model_type[i] == VARIANCE_CHANGE) ||
          (model_type[i] == ORDINAL_GAUSSIAN_CHANGE) || (model_type[i] == LINEAR_MODEL_CHANGE)) {
        if ((index != I_DEFAULT) || (!common_contrast)) {
          for (j = 0;j < nb_sequence;j++) {
            if ((index == I_DEFAULT) || (index == j)) {
              for (k = 0;k < nb_segment;k++) {
//                if (residual[j][k] > 0.) {
                if (residual[j][k] > (change_point[k + 1] - change_point[k]) * ROUNDOFF_ERROR) {
                  segmentation_likelihood -= ((double)(change_point[k + 1] - change_point[k]) / 2.) * (logl(residual[j][k] /
                                               (change_point[k + 1] - change_point[k])) + log(2 * M_PI) + 1);
/*                  segmentation_likelihood -= ((double)(change_point[k + 1] - change_point[k]) / 2.) * (logl(residual[j][k] /
                                               (change_point[k + 1] - change_point[k])) + log(2 * M_PI)) +
                                             (double)(change_point[k + 1] - change_point[k]) / 2.; */
                }
                else {
                  segmentation_likelihood = D_INF;
                  break;
                }
              }
            }

            if (segmentation_likelihood == D_INF) {
              break;
            }
          }
        }

        else {
          for (j = 0;j < nb_segment;j++) {
//            if (residual[0][j] > 0.) {
            if (residual[0][j] > nb_sequence * (change_point[j + 1] - change_point[j]) * ROUNDOFF_ERROR) {
              segmentation_likelihood -= ((double)(nb_sequence * (change_point[j + 1] - change_point[j])) / 2.) * (logl(residual[0][j] /
                                           (nb_sequence * (change_point[j + 1] - change_point[j]))) + log(2 * M_PI) + 1);
/*              segmentation_likelihood -= ((double)(nb_sequence * (change_point[j + 1] - change_point[j])) / 2.) * (logl(residual[0][j] /
                                           (nb_sequence * (change_point[j + 1] - change_point[j]))) + log(2 * M_PI)) +
                                         (double)(nb_sequence * (change_point[j + 1] - change_point[j])) / 2.; */
            }
            else {
              segmentation_likelihood = D_INF;
              break;
            }
          }
        }
      }

      else if ((model_type[i] == AUTOREGRESSIVE_MODEL_CHANGE) || (model_type[i] == STATIONARY_AUTOREGRESSIVE_MODEL_CHANGE)) {
        if ((index != I_DEFAULT) || (!common_contrast)) {
          for (j = 0;j < nb_sequence;j++) {
            if ((index == I_DEFAULT) || (index == j)) {
              for (k = 0;k < nb_segment;k++) {
//                if (residual[j][k] > 0.) {
                if (residual[j][k] > (change_point[k + 1] - change_point[k] - 1) * ROUNDOFF_ERROR) {
                  segmentation_likelihood -= ((double)(change_point[k + 1] - change_point[k] - 1) / 2.) * (logl(residual[j][k] /
                                               (change_point[k + 1] - change_point[k] - 1)) + log(2 * M_PI) + 1);
                }
                else {
                  segmentation_likelihood = D_INF;
                  break;
                }
              }
            }

            if (segmentation_likelihood == D_INF) {
              break;
            }
          }
        }

        else {
          for (j = 0;j < nb_segment;j++) {
//            if (residual[0][j] > 0.) {
            if (residual[0][j] > nb_sequence * (change_point[j + 1] - change_point[j] - 1) * ROUNDOFF_ERROR) {
              segmentation_likelihood -= ((double)(nb_sequence * (change_point[j + 1] - change_point[j] - 1)) / 2.) * (logl(residual[0][j] /
                                           (nb_sequence * (change_point[j + 1] - change_point[j] - 1))) + log(2 * M_PI) + 1);
            }
            else {
              segmentation_likelihood = D_INF;
              break;
            }
          }
        }
      }

      if (segmentation_likelihood == D_INF) {
        break;
      }
    }

    if ((model_type[0] == MEAN_CHANGE) || (model_type[0] == INTERCEPT_SLOPE_CHANGE)) {
      if (index != I_DEFAULT) {
        segmentation_likelihood = 0.;

        for (i = 1;i < nb_segment;i++) {
          residual[index][0] += residual[index][i];
        }

//        if (residual[index][0] > 0.) {
        if (residual[index][0] > length[index] * ROUNDOFF_ERROR) {
          segmentation_likelihood = -((double)length[index] / 2.) * (logl(residual[index][0] / length[index]) +
                                      log(2 * M_PI) + 1);
/*          segmentation_likelihood = -((double)length[index] / 2.) * (logl(residual[index][0] / (length[index] - nb_segment)) +
                                      log(2 * M_PI)) - (double)(length[index] - nb_segment) / 2.; */
        }
        else {
          segmentation_likelihood = D_INF;
        }
      }

      else {
        if (!common_contrast) {
          for (i = 0;i < nb_segment;i++) {
            for (j = 1;j < nb_sequence;j++) {
              residual[0][i] += residual[j][i];
            }
          }
        }

        for (i = 1;i < nb_segment;i++) {
          residual[0][0] += residual[0][i];
        }

//        if (residual[0][0] > 0.) {
        if (residual[0][0] > nb_sequence * length[0] * ROUNDOFF_ERROR) {
          segmentation_likelihood = -((double)(nb_sequence * length[0]) / 2.) *
                                     (logl(residual[0][0] / (nb_sequence * length[0])) + log(2 * M_PI) + 1);
/*          segmentation_likelihood = -((double)(nb_sequence * length[0]) / 2.) *
                                     (logl(residual[0][0] / (nb_sequence * (length[0] - nb_segment))) + log(2 * M_PI)) -
                                     (double)(nb_sequence * (length[0] - nb_segment)) / 2.; */
        }
        else {
          segmentation_likelihood = D_INF;
        }
      }
    }

    if (segmentation_likelihood != D_INF) {
      if (index != I_DEFAULT) {
        iseq = new Sequences(*this , 1 , &index);
        seq = new Sequences(*iseq , ADD_STATE_VARIABLE);
        delete iseq;
      }
      else {
        seq = new Sequences(*this , ADD_STATE_VARIABLE);
      }

      for (i = 0;i < seq->nb_sequence;i++) {
        for (j = 0;j < nb_segment;j++) {
          for (k = change_point[j];k < change_point[j + 1];k++) {
            seq->int_sequence[i][0][k] = j;
          }
        }
      }

      seq->min_value[0] = 0;
      seq->max_value[0] = nb_segment - 1;

      seq->build_marginal_frequency_distribution(0);

      if (display) {
        segment_penalty = 0.;
        for (i = 0;i < nb_segment;i++) {
          segment_penalty += log((double)(change_point[i + 1] - change_point[i]));
        }

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
                                      change_point , continuity);

      if ((output == SEQUENCE) || (output == ABSOLUTE_RESIDUAL)) {
        delete seq;
      }
    }

    else {
      oseq = NULL;
      error.update(SEQ_error[SEQR_SEGMENTATION_FAILURE]);
    }

    for (i = 0;i < nb_variable;i++) {
      delete [] seq_mean[i];
      delete [] rank[i];
    }
    delete [] seq_mean;
    delete [] rank;

    delete [] frequency;

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

    if (index_param_type == IMPLICIT_TYPE) {
      delete [] seq_index_parameter;
    }
  }

  delete [] change_point;

  return oseq;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Segmentation of a single sequence or a sample of sequences.
 *
 *  \param[in] error           reference on a StatError object,
 *  \param[in] display         flag for displaying the segmentation,
 *  \param[in] iidentifier     sequence identifier,
 *  \param[in] nb_segment      number of segments,
 *  \param[in] ichange_point   change points,
 *  \param[in] model_type      segment model types,
 *  \param[in] common_contrast flag contrast functions common to the individuals,
 *  \param[in] shape_parameter negative binomial shape parameters,
 *  \param[in] output          output (sequence or residuals),
 *  \param[in] continuity      flag continuous piecewise linear function.
 *
 *  \return                    Sequences object.
 */
/*--------------------------------------------------------------*/

Sequences* Sequences::segmentation(StatError &error , bool display , int iidentifier ,
                                   int nb_segment , vector<int> ichange_point , vector<segment_model> model_type ,
                                   bool common_contrast , vector<double> shape_parameter ,
                                   sequence_type output , bool continuity) const

{
  return segmentation(error , display , iidentifier , nb_segment , ichange_point.data() , model_type.data() ,
                      common_contrast , shape_parameter.data() , output , continuity);
}


};  // namespace sequence_analysis
