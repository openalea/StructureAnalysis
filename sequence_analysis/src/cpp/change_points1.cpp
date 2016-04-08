/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       V-Plants: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2016 CIRAD/INRA/Inria Virtual Plants
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

#include "tool/config.h"

#include "sequences.h"
#include "sequence_label.h"

using namespace std;
using namespace boost::math;
using namespace stat_tool;


namespace sequence_analysis {



/*--------------------------------------------------------------*
 *
 *  Determination of the width of a column of reals.
 *
 *  arguments: number of values, pointer on real values.
 *
 *--------------------------------------------------------------*/

int column_width(int nb_value , const long double *value)

{
  register int i;
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


/*--------------------------------------------------------------*
 *
 *  Computation of the log of the factorial of a value.
 *
 *  argument: value.
 *
 *--------------------------------------------------------------*/

double log_factorial(int value)

{
  register int i;
  double log_factorial;


  log_factorial = 0.;
  for (i = 2;i <= value;i++) {
    log_factorial += log((double)i);
  }

  return log_factorial;
}


/*--------------------------------------------------------------*
 *
 *  Computation of the log of a binomial coefficient for a negative binomial distribution.
 *
 *  arguments: inf bound, shape parameter, value.
 *
 *--------------------------------------------------------------*/

double log_binomial_coefficient(int inf_bound , double parameter , int value)

{
  register int i;
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


/*--------------------------------------------------------------*
 *
 *  Empirical determination of the hyperparameters of a gamma prior
 *  distribution for a Poisson distribution.
 *
 *  arguments: sequence index, variable index, pointers on the hyperparameters.
 *
 *--------------------------------------------------------------*/

void Sequences::gamma_hyperparameter_computation(int index , int variable ,
                                                 double *hyperparam) const

{
  register int i;
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


/*--------------------------------------------------------------*
 *
 *  Empirical determination of the hyperparameters of a Gaussian-gamma prior
 *  distribution for a Gaussian distribution.
 *
 *  arguments: sequence index, variable index, pointers on the hyperparameters.
 *
 *--------------------------------------------------------------*/

void Sequences::gaussian_gamma_hyperparameter_computation(int index , int variable ,
                                                          double *hyperparam) const

{
  register int i;
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


/*--------------------------------------------------------------*
 *
 *  Computation of the number of free parameters.
 *
 *  arguments: sequence index, number of segments, segment model types,
 *             flag contrast functions common to the individuals.
 *
 *--------------------------------------------------------------*/

int Sequences::nb_parameter_computation(int index , int nb_segment , segment_model *model_type ,
                                        bool common_contrast) const

{
  bool *used_output;
  register int i , j , k , m;
  int nb_parameter , max_nb_value , *psegment;


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
              psegment = int_sequence[j][0];

              for (k = 0;k < length[j];k++) {
                if ((k == 0) || ((k > 0) && (*psegment != *(psegment - 1)))) {
                  for (m = 0;m < marginal_distribution[i]->nb_value;m++) {
                    used_output[m] = false;
                  }
                  nb_parameter--;
                }
                psegment++;

                if (!used_output[int_sequence[j][i][k]]) {
                  nb_parameter++;
                  used_output[int_sequence[j][i][k]] = true;
                }
              }
            }
          }
        }

        else {
          psegment = int_sequence[0][0];

          for (j = 0;j < length[0];j++) {
            if ((j == 0) || ((j > 0) && (*psegment != *(psegment - 1)))) {
              for (k = 0;k < marginal_distribution[i]->nb_value;k++) {
                used_output[k] = false;
              }
              nb_parameter--;
            }
            psegment++;

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
          nb_parameter += 2 * nb_sequence * nb_segment;
        }
      }

      else if (model_type[i - 1] == LINEAR_MODEL_CHANGE) {
        if ((index != I_DEFAULT) || (common_contrast)) {
          nb_parameter += 3 * nb_segment;
        }
        else {
          nb_parameter += 3 * nb_sequence * nb_segment;
        }
      }

      else if (model_type[i - 1] == VARIANCE_CHANGE) {
        nb_parameter += nb_segment + 1;
      }
    }

    delete [] used_output;
  }

  return nb_parameter;
}


/*--------------------------------------------------------------*
 *
 *  Computation of the log-likelihood in the case of a single segment.
 *
 *  arguments: sequence index, segment model types,
 *             flag contrast functions common to the individuals,
 *             negative binomial distribution parameters,
 *             ranks (ordinal variables).
 *
 *--------------------------------------------------------------*/

double Sequences::one_segment_likelihood(int index , segment_model *model_type , bool common_contrast ,
                                         double *shape_parameter , double **rank) const

{
  register int i , j , k;
  int max_nb_value , seq_length , count , *frequency , *inf_bound_parameter , *seq_index_parameter;
  double sum , factorial_sum , binomial_coeff_sum , proba , diff , index_parameter_sum ,
         index_parameter_diff , likelihood;
  long double index_parameter_square_sum , square_sum , mix_square_sum , *residual;


  max_nb_value = 0;
  inf_bound_parameter = new int[nb_variable];

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
        (((model_type[i - 1] == GAUSSIAN_CHANGE) || (model_type[i - 1] == LINEAR_MODEL_CHANGE) ||
          (model_type[i - 1] == VARIANCE_CHANGE) || (model_type[i - 1] == ORDINAL_GAUSSIAN_CHANGE)) &&
         (!square_sum))) {
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

    else if (model_type[j - 1] == POISSON_CHANGE) {
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
              sum = int_sequence[j][i][0];
              residual[j] = 0.;

              for (k = 1;k < length[j];k++) {
                diff = int_sequence[j][i][k] - sum / k;
                residual[j] += ((double)k / (double)(k + 1)) * diff * diff;
                sum += int_sequence[j][i][k];
              }
            }
          }
        }

        else {
          for (j = 0;j < nb_sequence;j++) {
            if ((index == I_DEFAULT) || (index == j)) {
              sum = real_sequence[j][i][0];
              residual[j] = 0.;

              for (k = 1;k < length[j];k++) {
                diff = real_sequence[j][i][k] - sum / k;
                residual[j] += ((double)k / (double)(k + 1)) * diff * diff;
                sum += real_sequence[j][i][k];
              }
            }
          }
        }
      }

      else {
        residual[0] = 0.;
        sum = 0.;
        count = 0;

        if (type[i] != REAL_VALUE) {
          for (j = 0;j < length[0];j++) {
            for (k = 0;k < nb_sequence;k++) {
              if (count > 0) {
                diff = int_sequence[k][i][j] - sum / count;
                residual[0] += ((double)count / (double)(count + 1)) * diff * diff;
              }
              count++;
              sum += int_sequence[k][i][j];
            }
          }
        }

        else {
          for (j = 0;j < length[0];j++) {
            for (k = 0;k < nb_sequence;k++) {
              if (count > 0) {
                diff = real_sequence[k][i][j] - sum / count;
                residual[0] += ((double)count / (double)(count + 1)) * diff * diff;
              }
              count++;
              sum += real_sequence[k][i][j];
            }
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
              index_parameter_square_sum = 0.;
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
        index_parameter_square_sum = 0.;
        index_parameter_sum = nb_sequence * seq_index_parameter[0];
        square_sum = 0.;
        mix_square_sum = 0.;
        count = 1;

        if (type[i] != REAL_VALUE) {
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
          }
        }

        else {
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
          }
        }

        if (index_parameter_square_sum > 0.) {
          residual[0] = square_sum - mix_square_sum * mix_square_sum / index_parameter_square_sum;
        }
        else {
          residual[0] = 0.;
        }
      }
    }

    else if (model_type[i - 1] == ORDINAL_GAUSSIAN_CHANGE) {
      residual[index] = 0.;
      sum = rank[i][int_sequence[index][i][0]];

      for (j = 1;j < length[index];j++) {
        diff = rank[i][int_sequence[index][i][j]] - sum / j;
        residual[index] += ((double)j / (double)(j + 1)) * diff * diff;
        sum += rank[i][int_sequence[index][i][j]];
      }
    }

    if ((model_type[i - 1] == GAUSSIAN_CHANGE) || (model_type[i - 1] == LINEAR_MODEL_CHANGE) ||
        (model_type[i - 1] == VARIANCE_CHANGE) || (model_type[i - 1] == ORDINAL_GAUSSIAN_CHANGE)) {
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


/*--------------------------------------------------------------*
 *
 *  Computation of piecewise linear functions.
 *
 *  arguments: sequence index, variable, number of segments, segment model type,
 *             flag contrast functions common to the individuals, change points,
 *             index parameters, piecewise linear functions, means,
 *             variances or residual variances, global variance or residual variance,
 *             slopes, intercepts, correlation coefficients, slope standard deviations,
 *             index parameter mean and variance.
 *
 *--------------------------------------------------------------*/

double Sequences::piecewise_linear_function(int index , int variable , int nb_segment , segment_model model_type ,
                                            bool common_contrast , int *change_point , int *seq_index_parameter ,
                                            double **piecewise_function , double **imean , double **variance ,
                                            double *global_variance , double **iintercept , double **islope ,
                                            double **correlation , double **slope_standard_deviation ,
                                            double **iindex_parameter_mean , long double **iindex_parameter_variance) const

{
  register int i , j , k;
  double likelihood , mean , diff , diff_sum , index_parameter_mean , response_mean , slope , intercept;
  long double square_sum , global_square_sum , index_parameter_variance , response_variance , covariance;


  if ((model_type == POISSON_CHANGE) || (model_type == NEGATIVE_BINOMIAL_0_CHANGE) ||
      (model_type == NEGATIVE_BINOMIAL_1_CHANGE) || (model_type == GAUSSIAN_CHANGE) ||
      (model_type == MEAN_CHANGE) || (model_type == VARIANCE_CHANGE) ||
      (model_type == BAYESIAN_POISSON_CHANGE) || (model_type == BAYESIAN_GAUSSIAN_CHANGE)) {
    if (((model_type == GAUSSIAN_CHANGE) || (model_type == MEAN_CHANGE)) &&
        ((variance) || (global_variance))) {
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

                if ((model_type == GAUSSIAN_CHANGE) && (likelihood != D_INF)) {
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
      for (i = 0;i < nb_segment;i++) {
        mean = 0.;

        if (type[variable] != REAL_VALUE) {
          for (j = change_point[i];j < change_point[i + 1];j++) {
            for (k = 0;k < nb_sequence;k++) {
              mean += int_sequence[k][variable][j];
            }
          }
        }
        else {
          for (j = change_point[i];j < change_point[i + 1];j++) {
            for (k = 0;k < nb_sequence;k++) {
              mean += real_sequence[k][variable][j];
            }
          }
        }
        mean /= (nb_sequence * (change_point[i + 1] - change_point[i]));

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
            for (j = change_point[i];j < change_point[i + 1];j++) {
              for (k = 0;k < nb_sequence;k++) {
                diff = int_sequence[k][variable][j] - mean;
                square_sum += diff * diff;
              }
            }
          }
          else {
            for (j = change_point[i];j < change_point[i + 1];j++) {
              for (k = 0;k < nb_sequence;k++) {
                diff = real_sequence[k][variable][j] - mean;
                square_sum += diff * diff;
              }
            }
          }

          if (global_variance) {
            global_square_sum += square_sum;
          }
          variance[0][i] = square_sum / (nb_sequence * (change_point[i + 1] - change_point[i]) - 1);

          if ((model_type == GAUSSIAN_CHANGE) && (likelihood != D_INF)) {
            if (square_sum > nb_sequence * (change_point[i + 1] - change_point[i]) * ROUNDOFF_ERROR) {
              likelihood -= ((double)(nb_sequence * (change_point[i + 1] - change_point[i])) / 2.) * (log(square_sum /
                              (nb_sequence * (change_point[i + 1] - change_point[i]))) + log(2 * M_PI) + 1);
            }
            else {
              likelihood = D_INF;
            }
          }
        }
      }
    }

    if (global_variance) {
      if (index != I_DEFAULT) {
        global_variance[variable] = global_square_sum / (length[index] - nb_segment);
      }
      else {
        global_variance[variable] = global_square_sum / (nb_sequence * length[0] - nb_segment);
      }

      if (model_type == MEAN_CHANGE) {
        if (index != I_DEFAULT) {
          likelihood = -((double)length[index] / 2.) * (log(global_square_sum /
                          length[index]) + log(2 * M_PI) + 1);
        }
        else {
          likelihood = -((double)(nb_sequence * length[0]) / 2.) * (log(global_square_sum /
                          (nb_sequence * length[0])) + log(2 * M_PI) + 1);
        }
      }
    }
  }

  if ((model_type == LINEAR_MODEL_CHANGE) || (model_type == INTERCEPT_SLOPE_CHANGE)) {
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
      if (index != I_DEFAULT) {
        global_variance[variable] = global_square_sum / (length[index] - 2 * nb_segment);
      }
      else {
        global_variance[variable] = global_square_sum / (nb_sequence * length[0] - 2 * nb_segment);
      }

      if (model_type == INTERCEPT_SLOPE_CHANGE) {
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
    }
  }

  return likelihood;
}


/*--------------------------------------------------------------*
 *
 *  Writing of piecewise linear functions.
 *
 *  arguments: stream, sequence index, variable, number of segments, segment model type,
 *             flag contrast functions common to the individuals, change points,
 *             index parameters, means, variances or residual variances,
 *             slopes, intercepts, correlation coefficients, slope standard deviations.
 *
 *--------------------------------------------------------------*/

ostream& Sequences::piecewise_linear_function_ascii_print(ostream &os , int index , int variable , int nb_segment ,
                                                          segment_model model_type , bool common_contrast , int *change_point ,
                                                          int *seq_index_parameter , double **mean , double **variance ,
                                                          double **intercept , double **slope , double **correlation ,
                                                          double **slope_standard_deviation , double **index_parameter_mean ,
                                                          long double **index_parameter_variance) const

{
  register int i , j;
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
    os << SEQ_label[SEQL_SEGMENT] << " " << STAT_label[STATL_MEAN] << ", "
       << STAT_label[STATL_STANDARD_DEVIATION] << ": ";

    if ((index != I_DEFAULT) || (!common_contrast)) {
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
      for (i = 0;i < nb_segment;i++) {
        os << mean[0][i] << " " << sqrt(variance[0][i]);
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
           << STAT_label[STATL_SLOPE] << ", " << STAT_label[STATL_CORRELATION_COEFF] << " ("
           << STAT_label[STATL_LIMIT_CORRELATION_COEFF] << "), "
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
                 << slope[index][i] + test->value * slope_standard_deviation[index][i] << "), ";
//               << "slope_standard_deviation: " << slope_standard_deviation[index][i] << " "
            }
            os << correlation[index][i] << " (-/+"
               << test->value / sqrt(test->value * test->value + change_point[i + 1] - change_point[i] - 2) << "), ";
            os << sqrt(variance[index][i]);

            if (i < nb_segment - 1) {
              os << ", " <<  intercept[index][i + 1] + slope[index][i + 1] * seq_index_parameter[change_point[i + 1]] -
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
                 << slope[0][i] + test->value * slope_standard_deviation[0][i] << "), ";
//               << "slope_standard_deviation: " << slope_standard_deviation[0][i] << " "
            }
            os << correlation[0][i] << " (-/+"
               << test->value / sqrt(test->value * test->value + nb_sequence * (change_point[i + 1] - change_point[i]) - 2) << "), ";
            os << sqrt(variance[0][i]);

            if (i < nb_segment - 1) {
              os << ", " <<  intercept[0][i + 1] + slope[0][i + 1] * seq_index_parameter[change_point[i + 1]] -
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

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Writing of piecewise linear functions at the spreadsheet format.
 *
 *  arguments: stream, sequence index, variable, number of segments, segment model type,
 *             flag contrast functions common to the individuals, change points,
 *             index parameters, means, means, variances or residual variances,
 *             global variances or residual variances, slopes, intercepts.
 *
 *--------------------------------------------------------------*/

ostream& Sequences::piecewise_linear_function_spreadsheet_print(ostream &os , int index , int variable , int nb_segment ,
                                                                segment_model model_type , bool common_contrast ,
                                                                int *change_point , int *seq_index_parameter ,
                                                                double **mean , double **variance ,
                                                                double **intercept , double **slope) const

{
  register int i , j;


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
        os << "\t" << mean[0][i] << "\t" << variance[0][i] << "\t";
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
    os << SEQ_label[SEQL_SEGMENT] << "\t" << STAT_label[STATL_MEAN] << "\t"
       << STAT_label[STATL_STANDARD_DEVIATION];

    if ((index != I_DEFAULT) || (!common_contrast)) {
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
      for (i = 0;i < nb_segment;i++) {
        os << "\t" << mean[0][i] << "\t" << sqrt(variance[0][i]) << "\t";
      }
    }
    os << endl;
  }

  else if ((model_type == LINEAR_MODEL_CHANGE) || (model_type == INTERCEPT_SLOPE_CHANGE)) {
    if (nb_variable > 2) {
      os << STAT_label[STATL_VARIABLE] << "\t" << variable << "\t";
    }
    os << SEQ_label[SEQL_PIECEWISE_LINEAR_FUNCTION] << "\t"
       << SEQ_label[SEQL_SEGMENT] << "\t" << STAT_label[STATL_INTERCEPT] << "\t"
       << STAT_label[STATL_SLOPE] << "\t" << STAT_label[STATL_RESIDUAL] << " "
       << STAT_label[STATL_STANDARD_DEVIATION];

    if ((index != I_DEFAULT) || (!common_contrast)) {
      for (i = 0;i < nb_sequence;i++) {
        if ((index == I_DEFAULT) || (index == i)) {
          for (j = 0;j < nb_segment;j++) {
            os << "\t" << intercept[i][j] + slope[i][j] * seq_index_parameter[change_point[j]] << "->";
            if (i < nb_segment - 1) {
              os << intercept[i][j] + slope[i][j] * seq_index_parameter[change_point[j + 1]];
            }
            else {
              os << intercept[i][j] + slope[i][j] * seq_index_parameter[change_point[j + 1] - 1];
            }
          }
          for (j = 0;j < nb_segment;j++) {
            os << "\t\t" << intercept[i][j] << "\t" << slope[i][j] << "\t" << sqrt(variance[i][j]);
          }
          if ((index == I_DEFAULT) && (i < nb_sequence - 1)) {
            os << "\t\t";
          }
        }
      }
    }

    else {
      for (i = 0;i < nb_segment;i++) {
        os << "\t" << intercept[0][i] + slope[0][i] * seq_index_parameter[change_point[i]] << "->";
        if (i < nb_segment - 1) {
          os << intercept[0][i] + slope[0][i] * seq_index_parameter[change_point[i + 1]];
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


/*--------------------------------------------------------------*
 *
 *  Computation of continuous piecewise linear functions.
 *
 *  arguments: stream, sequence index, variable, number of segments, segment model type,
 *             flag contrast functions common to the individuals, change points,
 *             index parameters, piecewise linear functions,
 *             means, slopes, intercepts, corrected intercepts and slopes.
 *
 *--------------------------------------------------------------*/

double Sequences::continuous_piecewise_linear_function(ostream &os , int index , int variable , int nb_segment ,
                                                       segment_model model_type , bool common_contrast ,
                                                       int *change_point , int *seq_index_parameter ,
                                                       double *intercept , double *slope ,
                                                       double *corrected_intercept , double *corrected_slope) const

{
  register int i , j , k;
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

# ifdef MESSAGE
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
# endif

  delete [] predicted_value;
  delete [] variance;

  return likelihood;
}


/*--------------------------------------------------------------*
 *
 *  Output of a segmentation of a single sequence or a sample of sequences.
 *
 *  arguments: sequence index, number of segments, segment model types,
 *             flag contrast functions common to the individuals, stream,
 *             output (sequence or residuals), change points,
 *             flag continuous piecewise linear function.
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::segmentation_output(int nb_segment , segment_model *model_type ,
                                          bool common_contrast , ostream &os , sequence_type output ,
                                          int *ichange_point , bool continuity)

{
  bool *piecewise_function_flag;
  register int i , j , k , m , n;
  int *change_point , *psegment , *seq_index_parameter = NULL;
  double likelihood , corrected_likelihood , diff , buff , change_point_amplitude , *global_variance ,
         ***mean , ***variance , ***index_parameter_mean , ***intercept , ***slope , ***correlation ,
         ***slope_standard_deviation , ***corrected_intercept , ***corrected_slope;
  long double ***index_parameter_variance;
  Test *test;
  Sequences *seq;


  if (ichange_point) {
    change_point = ichange_point;
  }

  else {
    change_point = new int[nb_segment + 1];
    change_point[0] = 0;

    psegment = int_sequence[0][0] + 1;
    i = 1;
    for (j = 1;j < length[0];j++) {
      if (*psegment != *(psegment - 1)) {
        change_point[i++] = j;
      }
      psegment++;
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
      switch (common_contrast) {

      case false : {
        mean[i] = new double*[nb_sequence];
        variance[i] = new double*[nb_sequence];

        for (j = 0;j < nb_sequence;j++) {
          mean[i][j] = new double[nb_segment];
          variance[i][j] = new double[nb_segment];
        }
        break;
      }

      case true : {
        mean[i] = new double*[1];
        mean[i][0] = new double[nb_segment];
        variance[i] = new double*[1];
        variance[i][0] = new double[nb_segment];
        break;
      }
      }
    }

    else if ((model_type[i - 1] == LINEAR_MODEL_CHANGE) || (model_type[0] == INTERCEPT_SLOPE_CHANGE)) {
      switch (common_contrast) {

      case false : {
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
        break;
      }

      case true : {
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
        break;
      }
      }
    }

    if ((((i == 1) && ((model_type[0] == MEAN_CHANGE) || (model_type[0] == INTERCEPT_SLOPE_CHANGE))) ||
         (model_type[i - 1] == GAUSSIAN_CHANGE) || (model_type[i - 1] == LINEAR_MODEL_CHANGE) ||
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
          (model_type[0] == INTERCEPT_SLOPE_CHANGE) || (model_type[i - 1] == BAYESIAN_POISSON_CHANGE) ||
          (model_type[i - 1] == BAYESIAN_GAUSSIAN_CHANGE)) {
        piecewise_function_flag[i] = true;
      }
      else {
        piecewise_function_flag[i] = false;
      }
    }

    seq = new Sequences(*this , piecewise_function_flag);
  }

  for (i = 1;i < nb_variable;i++) {
    if ((model_type[i - 1] == POISSON_CHANGE) || (model_type[i - 1] == NEGATIVE_BINOMIAL_0_CHANGE) ||
        (model_type[i - 1] == NEGATIVE_BINOMIAL_1_CHANGE) || (model_type[i - 1] == GAUSSIAN_CHANGE) ||
        (model_type[0] == MEAN_CHANGE) || (model_type[i - 1] == VARIANCE_CHANGE) ||
        (model_type[i - 1] == LINEAR_MODEL_CHANGE) || (model_type[0] == INTERCEPT_SLOPE_CHANGE) ||
        (model_type[i - 1] == BAYESIAN_POISSON_CHANGE) || (model_type[i - 1] == BAYESIAN_GAUSSIAN_CHANGE)) {
      likelihood = piecewise_linear_function((nb_sequence == 1 ? 0 : I_DEFAULT) , i , nb_segment ,
                                              model_type[i - 1] , common_contrast , change_point ,
                                              seq_index_parameter , NULL , mean[i] , variance[i] ,
                                              global_variance , intercept[i] , slope[i] ,
                                              correlation[i] , slope_standard_deviation[i] ,
                                              index_parameter_mean[i] , index_parameter_variance[i]);
#     ifdef MESSAGE
      if (((i == 1) && ((model_type[0] == MEAN_CHANGE) || (model_type[0] == INTERCEPT_SLOPE_CHANGE))) ||
          (model_type[i - 1] == GAUSSIAN_CHANGE) || (model_type[i - 1] == LINEAR_MODEL_CHANGE)) {
        os << "2 * " << STAT_label[STATL_LIKELIHOOD] << ": " << 2 * likelihood << endl;
      }
#     endif

    }
  }

# ifdef MESSAGE
  if (!ichange_point) {
    os << (nb_segment == 2 ? SEQ_label[SEQL_CHANGE_POINT] : SEQ_label[SEQL_CHANGE_POINTS]) << ": ";

    for (i = 1;i < nb_segment;i++) {
      os << seq_index_parameter[change_point[i]];
      if (i < nb_segment - 1) {
        os << ", ";
      }
    }
  }

  if ((index_interval) && (index_interval->variance > 0.)) {
    if (!ichange_point) {
      os << "   ";
    }
    os << SEQ_label[SEQL_SEGMENT_SAMPLE_SIZE] << ": ";
    for (i = 0;i < nb_segment;i++) {
      os << nb_sequence * (change_point[i + 1] - change_point[i]);
      if (i < nb_segment - 1) {
        os << ", ";
      }
    }
    os << endl;
  }

  else if (!ichange_point) {
    os << endl;
  }

  if (nb_variable > 2) {
    os << "\n";
  }

  for (i = 1;i < nb_variable;i++) {
    piecewise_linear_function_ascii_print(os , (nb_sequence == 1 ? 0 : I_DEFAULT) , i , nb_segment , model_type[i - 1] ,
                                          common_contrast , change_point , seq_index_parameter ,
                                          mean[i] , variance[i] , intercept[i] , slope[i] ,
                                          correlation[i] , slope_standard_deviation[i] ,
                                          index_parameter_mean[i] , index_parameter_variance[i]);

    if ((model_type[i - 1] == GAUSSIAN_CHANGE) || (model_type[0] == MEAN_CHANGE) ||
        (model_type[0] == VARIANCE_CHANGE) || (model_type[i - 1] == BAYESIAN_GAUSSIAN_CHANGE)) {
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

        os << STAT_label[STATL_MEAN] << " " << SEQ_label[SEQL_CHANGE_POINT_AMPLITUDE] << ": "
           << change_point_amplitude << "   ";
      }

      if ((model_type[i - 1] == GAUSSIAN_CHANGE) || (model_type[0] == MEAN_CHANGE) ||
          (model_type[i - 1] == BAYESIAN_GAUSSIAN_CHANGE)) {
        os << SEQ_label[SEQL_GLOBAL_STANDARD_DEVIATION] << ": " << sqrt(global_variance[i]);
        if (nb_segment > 1) {
          os << "   " << STAT_label[STATL_RATIO] << ": "
             << change_point_amplitude / sqrt(global_variance[i]);
        }
      }
      os << endl;
    }

    if (model_type[0] == MEAN_CHANGE) {
      os << SEQ_label[SEQL_GLOBAL_STANDARD_DEVIATION] << ": "  << sqrt(global_variance[i]) << endl;
    }
    else if (model_type[0] == INTERCEPT_SLOPE_CHANGE) {
      os << SEQ_label[SEQL_GLOBAL_RESIDUAL_STANDARD_DEVIATION] << ": "  << sqrt(global_variance[i]) << endl;
    }

    if (continuity) {
      corrected_likelihood = continuous_piecewise_linear_function(os , (nb_sequence == 1 ? 0 : I_DEFAULT) , i ,
                                                                  nb_segment , model_type[i - 1] , common_contrast ,
                                                                  change_point , seq_index_parameter , intercept[i][0] ,
                                                                  slope[i][0] , corrected_intercept[i][0] , corrected_slope[i][0]);

      os << "2 * " << STAT_label[STATL_LIKELIHOOD] << ": "
          << 2 * corrected_likelihood << " | " << 2 * likelihood << endl;
    }
  }

# endif

  switch (output) {

  case SEQUENCE : {
    switch (common_contrast) {

    case false : {
      for (i = 0;i < nb_sequence;i++) {
        j = 1;
        for (k = 1;k < nb_variable;k++) {
          j++;
          if (piecewise_function_flag[k]) {
            if ((model_type[k - 1] != LINEAR_MODEL_CHANGE) && (model_type[k - 1] != INTERCEPT_SLOPE_CHANGE)) {
              for (m = 0;m < nb_segment;m++) {
                for (n = change_point[m];n < change_point[m + 1];n++) {
                  seq->real_sequence[i][j][n] = mean[k][i][m];
                }
              }
            }

            else {
              switch (continuity) {

              case false : {
                for (m = 0;m < nb_segment;m++) {
                  for (n = change_point[m];n < change_point[m + 1];n++) {
                    seq->real_sequence[i][j][n] = intercept[k][i][m] + slope[k][i][m] * seq_index_parameter[n];
                  }
                }
                break;
              }

              case true : {
                for (m = 0;m < nb_segment;m++) {
                  for (n = change_point[m];n < change_point[m + 1];n++) {
                    seq->real_sequence[i][j][n] = corrected_intercept[k][i][m] + corrected_slope[k][i][m] * seq_index_parameter[n];
                  }
                }
                break;
              }
              }
            }
            j++;
          }
        }
      }
      break;
    }

    case true : {
      for (i = 0;i < nb_sequence;i++) {
        j = 1;
        for (k = 1;k < nb_variable;k++) {
          j++;
          if (piecewise_function_flag[k]) {
            if ((model_type[k - 1] != LINEAR_MODEL_CHANGE) && (model_type[k - 1] != INTERCEPT_SLOPE_CHANGE)) {
              for (m = 0;m < nb_segment;m++) {
                for (n = change_point[m];n < change_point[m + 1];n++) {
                  seq->real_sequence[i][j][n] = mean[k][0][m];
                }
              }
            }

            else {
              switch (continuity) {

              case false : {
                for (m = 0;m < nb_segment;m++) {
                  for (n = change_point[m];n < change_point[m + 1];n++) {
                    seq->real_sequence[i][j][n] = intercept[k][0][m] + slope[k][0][m] * seq_index_parameter[n];
                  }
                }
                break;
              }

              case true : {
                for (m = 0;m < nb_segment;m++) {
                  for (n = change_point[m];n < change_point[m + 1];n++) {
                    seq->real_sequence[i][j][n] = corrected_intercept[k][0][m] + corrected_slope[k][0][m] * seq_index_parameter[n];
                  }
                }
                break;
              }
              }
            }
            j++;
          }
        }
      }
      break;
    }
    }
    break;
  }

  // residual computation

  case SUBTRACTION_RESIDUAL : {
    switch (common_contrast) {

    case false : {
      for (i = 0;i < nb_sequence;i++) {
        for (j = 1;j < nb_variable;j++) {
          if (type[j] != REAL_VALUE) {
            real_sequence[i][j] = new double[length[i]];

            if ((model_type[j - 1] != LINEAR_MODEL_CHANGE) && (model_type[j - 1] != INTERCEPT_SLOPE_CHANGE)) {
              for (k = 0;k < nb_segment;k++) {
                for (m = change_point[k];m < change_point[k + 1];m++) {
                  real_sequence[i][j][m] = int_sequence[i][j][m] - mean[j][i][k];
                }
              }
            }

            else {
              for (k = 0;k < nb_segment;k++) {
                for (m = change_point[k];m < change_point[k + 1];m++) {
                  real_sequence[i][j][m] = int_sequence[i][j][m] - (intercept[j][i][k] + slope[j][i][k] * seq_index_parameter[m]);
                }
              }
            }

            delete [] int_sequence[i][j];
            int_sequence[i][j] = NULL;
          }

          else {
            if ((model_type[j - 1] != LINEAR_MODEL_CHANGE) && (model_type[j - 1] != INTERCEPT_SLOPE_CHANGE)) {
              for (k = 0;k < nb_segment;k++) {
                for (m = change_point[k];m < change_point[k + 1];m++) {
                  real_sequence[i][j][m] -= mean[j][i][k];
                }
              }
            }

            else {
              for (k = 0;k < nb_segment;k++) {
                for (m = change_point[k];m < change_point[k + 1];m++) {
                  real_sequence[i][j][m] -= (intercept[j][i][k] + slope[j][i][k] * seq_index_parameter[m]);
                }
              }
            }
          }
        }
      }
      break;
    }

    case true : {
      for (i = 0;i < nb_sequence;i++) {
        for (j = 1;j < nb_variable;j++) {
          if (type[j] != REAL_VALUE) {
            real_sequence[i][j] = new double[length[i]];

            if ((model_type[j - 1] != LINEAR_MODEL_CHANGE) && (model_type[j - 1] != INTERCEPT_SLOPE_CHANGE)) {
              for (k = 0;k < nb_segment;k++) {
                for (m = change_point[k];m < change_point[k + 1];m++) {
                  real_sequence[i][j][m] = int_sequence[i][j][m] - mean[j][0][k];
                }
              }
            }

            else {
              for (k = 0;k < nb_segment;k++) {
                for (m = change_point[k];m < change_point[k + 1];m++) {
                  real_sequence[i][j][m] = int_sequence[i][j][m] - (intercept[j][0][k] + slope[j][0][k] * seq_index_parameter[m]);
                }
              }
            }

            delete [] int_sequence[i][j];
            int_sequence[i][j] = NULL;
          }

          else {
            if ((model_type[j - 1] != LINEAR_MODEL_CHANGE) && (model_type[j - 1] != INTERCEPT_SLOPE_CHANGE)) {
              for (k = 0;k < nb_segment;k++) {
                for (m = change_point[k];m < change_point[k + 1];m++) {
                  real_sequence[i][j][m] -= mean[j][0][k];
                }
              }
            }

            else {
              for (k = 0;k < nb_segment;k++) {
                for (m = change_point[k];m < change_point[k + 1];m++) {
                  real_sequence[i][j][m] -= (intercept[j][0][k] + slope[j][0][k] * seq_index_parameter[m]);
                }
              }
            }
          }
        }
      }
      break;
    }
    }
    break;
  }

  case DIVISION_RESIDUAL : {
    switch (common_contrast) {

    case false : {
      for (i = 0;i < nb_sequence;i++) {
        for (j = 1;j < nb_variable;j++) {
          if (type[j] != REAL_VALUE) {
            real_sequence[i][j] = new double[length[i]];

            if ((model_type[j - 1] != LINEAR_MODEL_CHANGE) && (model_type[j - 1] != INTERCEPT_SLOPE_CHANGE)) {
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

            else {
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

            delete [] int_sequence[i][j];
            int_sequence[i][j] = NULL;
          }

          else {
            if ((model_type[j - 1] != LINEAR_MODEL_CHANGE) && (model_type[j - 1] != INTERCEPT_SLOPE_CHANGE)) {
              for (k = 0;k < nb_segment;k++) {
                if (mean[j][i][k] != 0.) {
                  for (m = change_point[k];m < change_point[k + 1];m++) {
                    real_sequence[i][j][m] /= mean[j][i][k];
                  }
                }
              }
            }

            else {
              for (k = 0;k < nb_segment;k++) {
                for (m = change_point[k];m < change_point[k + 1];m++) {
                  if (intercept[j][i][k] + slope[j][i][k] * seq_index_parameter[m] != 0.) {
                    real_sequence[i][j][m] /= (intercept[j][i][k] + slope[j][i][k] * seq_index_parameter[m]);
                  }
                }
              }
            }
          }
        }
      }
      break;
    }

    case true : {
      for (i = 0;i < nb_sequence;i++) {
        for (j = 1;j < nb_variable;j++) {
          if (type[j] != REAL_VALUE) {
            real_sequence[i][j] = new double[length[i]];

            if ((model_type[j - 1] != LINEAR_MODEL_CHANGE) && (model_type[j - 1] != INTERCEPT_SLOPE_CHANGE)) {
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

            else {
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

            delete [] int_sequence[i][j];
            int_sequence[i][j] = NULL;
          }

          else {
            if ((model_type[j - 1] != LINEAR_MODEL_CHANGE) && (model_type[j - 1] != INTERCEPT_SLOPE_CHANGE)) {
              for (k = 0;k < nb_segment;k++) {
                if (mean[j][0][k] != 0.) {
                  for (m = change_point[k];m < change_point[k + 1];m++) {
                    real_sequence[i][j][m] /= mean[j][0][k];
                  }
                }
              }
            }

            else {
              for (k = 0;k < nb_segment;k++) {
                for (m = change_point[k];m < change_point[k + 1];m++) {
                  if (intercept[j][0][k] + slope[j][0][k] * seq_index_parameter[m] != 0.) {
                    real_sequence[i][j][m] /= (intercept[j][0][k] + slope[j][0][k] * seq_index_parameter[m]);
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
      switch (common_contrast) {

      case false : {
        for (j = 0;j < nb_sequence;j++) {
          delete [] mean[i][j];
          delete [] variance[i][j];
        }
        break;
      }

      case true : {
        delete [] mean[i][0];
        delete [] variance[i][0];
        break;
      }
      }

      delete [] mean[i];
      delete [] variance[i];
    }

    else if ((model_type[i - 1] == LINEAR_MODEL_CHANGE) || (model_type[0] == INTERCEPT_SLOPE_CHANGE)) {
     switch (common_contrast) {

      case false : {
        for (j = 0;j < nb_sequence;j++) {
          delete [] intercept[i][j];
          delete [] slope[i][j];
          delete [] variance[i][j];
          delete [] correlation[i][j];
          delete [] slope_standard_deviation[i][j];
          delete [] index_parameter_mean[i][j];
          delete [] index_parameter_variance[i][j];
        }
        break;
      }

      case true : {
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
        break;
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

  else {
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


/*--------------------------------------------------------------*
 *
 *  Segmentation of a single sequence or a sample of sequences.
 *
 *  arguments: reference on a StatError object, stream, sequence identifier,
 *             number of segments, change points, segment model types,
 *             flag contrast functions common to the individuals,
 *             negative binomial distribution parameters,
 *             output (sequence or residuals).
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::segmentation(StatError &error , ostream &os , int iidentifier ,
                                   int nb_segment , int *ichange_point , segment_model *model_type ,
                                   bool common_contrast , double *shape_parameter ,
                                   sequence_type output , bool continuity) const

{
  bool status = true;
  register int i , j , k , m;
  int index , segmentation_index , seq_length , count , max_nb_value , nb_parameter ,
      *change_point = NULL , *inf_bound_parameter , *frequency , *seq_index_parameter;
  double sum , factorial_sum , binomial_coeff_sum , proba , diff , index_parameter_sum ,
         index_parameter_diff , segmentation_likelihood , segment_penalty , penalized_likelihood ,
         *seq_mean , **rank;
  long double index_parameter_square_sum , square_sum , mix_square_sum , **residual;
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
          ((output == SUBTRACTION_RESIDUAL) || (output == DIVISION_RESIDUAL))) {
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

    // rank computation for ordinal variables

    rank = new double*[nb_variable];

    for (i = 0;i < nb_variable;i++) {
      if (model_type[i] == ORDINAL_GAUSSIAN_CHANGE) {
        rank[i] = marginal_distribution[i]->rank_computation();
      }
      else {
        rank[i] = NULL;
      }
    }

    max_nb_value = 0;
    inf_bound_parameter = new int[nb_variable];
    residual = NULL;

    seq_index_parameter = NULL;

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

      if (((i == 0) && ((model_type[0] == MEAN_CHANGE) || (model_type[0] == INTERCEPT_SLOPE_CHANGE))) ||
          (((model_type[i] == GAUSSIAN_CHANGE) || (model_type[i] == LINEAR_MODEL_CHANGE) ||
            (model_type[i] == VARIANCE_CHANGE) || (model_type[i] == ORDINAL_GAUSSIAN_CHANGE)) && (!residual))) {
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

    seq_mean = new double[nb_variable];

    for (i = 0;i < nb_variable;i++) {
      if (model_type[i] == VARIANCE_CHANGE) {
        seq_mean[i] = 0.;
        if (type[i] != REAL_VALUE) {
          for (j = 0;j < length[index];j++) {
            seq_mean[i] += int_sequence[index][i][j];
          }
        }
        else {
          for (j = 0;j < length[index];j++) {
            seq_mean[i] += real_sequence[index][i][j];
          }
        }
        seq_mean[i] /= length[index];
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

        else  {
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
                  residual[j][k] = 0.;
                  sum = int_sequence[j][i][change_point[k]];

                  for (m = change_point[k] + 1;m < change_point[k + 1];m++) {
                    diff = int_sequence[j][i][m] - sum / (m - change_point[k]);
                    residual[j][k] += ((double)(m - change_point[k]) / (double)(m - change_point[k] + 1)) *
                                      diff * diff;
                    sum += int_sequence[j][i][m];
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
                  sum = real_sequence[j][i][change_point[k]];

                  for (m = change_point[k] + 1;m < change_point[k + 1];m++) {
                    diff = real_sequence[j][i][m] - sum / (m - change_point[k]);
                    residual[j][k] += ((double)(m - change_point[k]) / (double)(m - change_point[k] + 1)) *
                                      diff * diff;
                    sum += real_sequence[j][i][m];
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
              }
            }
          }

          else {
            for (j = 0;j < nb_segment;j++) {
              residual[0][j] = 0.;
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
              }
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
                  index_parameter_square_sum = 0.;
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
                  index_parameter_square_sum = 0.;
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
              index_parameter_square_sum = 0.;
              index_parameter_sum = nb_sequence * seq_index_parameter[change_point[j]];
              square_sum = 0.;
              mix_square_sum = 0.;
              count = 1;
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
            for (j = 0;j < nb_segment;j++) {
              index_parameter_square_sum = 0.;
              index_parameter_sum = nb_sequence * seq_index_parameter[change_point[j]];
              square_sum = 0.;
              mix_square_sum = 0.;
              count = 1;
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

      else if (model_type[i] == VARIANCE_CHANGE) {
        for (j = 0;j < nb_segment;j++) {
          residual[index][j] = 0.;

          if (type[i] != REAL_VALUE) {
            for (k = change_point[j];k < change_point[j + 1];k++) {
              diff = int_sequence[index][i][k] - seq_mean[i];
              residual[index][j] += diff * diff;
            }
          }

          else {
            for (k = change_point[j];k < change_point[j + 1];k++) {
              diff = real_sequence[index][i][k] - seq_mean[i];
              residual[index][j] += diff * diff;
            }
          }
        }
      }

      else if (model_type[i] == ORDINAL_GAUSSIAN_CHANGE) {
        for (j = 0;j < nb_segment;j++) {
          residual[index][j] = 0.;
          sum = rank[i][int_sequence[index][i][change_point[j]]];

          for (k = change_point[j] + 1;k < change_point[j + 1];k++) {
            diff = rank[i][int_sequence[index][i][k]] - sum / (k - change_point[j]);
            residual[index][j] += ((double)(k - change_point[j]) / (double)(k - change_point[j] + 1)) *
                                  diff * diff;
            sum += rank[i][int_sequence[index][i][k]];
          }
        }
      }

      if ((model_type[i] == GAUSSIAN_CHANGE) || (model_type[i] == LINEAR_MODEL_CHANGE) ||
          (model_type[i] == VARIANCE_CHANGE) || (model_type[i] == ORDINAL_GAUSSIAN_CHANGE)) {
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

#   ifdef MESSAGE
    if (segmentation_likelihood != D_INF) {
      segment_penalty = 0.;
      for (i = 0;i < nb_segment;i++) {
        segment_penalty += log((double)(change_point[i + 1] - change_point[i]));
      }

      nb_parameter = seq->nb_parameter_computation((index == I_DEFAULT ? index : 0) , nb_segment , model_type ,
                                                   common_contrast);

      penalized_likelihood = 2 * segmentation_likelihood - nb_parameter *
                             log((double)((seq->nb_variable - 1) * seq->length[0])) - segment_penalty;

      os << "\n" << nb_segment << " " << (nb_segment == 1 ? SEQ_label[SEQL_SEGMENT] : SEQ_label[SEQL_SEGMENTS])
         << "   2 * " << STAT_label[STATL_LIKELIHOOD] << ": " << 2 * segmentation_likelihood << "   "
         << nb_parameter << " " << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS]
         << "   2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " (Modified "  << STAT_criterion_word[BIC] << "): "
         << penalized_likelihood << endl;
    }

    else {
      os << "\n" << nb_segment << " " << (nb_segment == 1 ? SEQ_label[SEQL_SEGMENT] : SEQ_label[SEQL_SEGMENTS])
         << "   " << STAT_label[STATL_LIKELIHOOD] << ": " << segmentation_likelihood << endl;
    }
#   endif

    oseq = seq->segmentation_output(nb_segment , model_type , common_contrast , os , output ,
                                    change_point , continuity);

    if (output == SEQUENCE) {
      delete seq;
    }

    for (i = 0;i < nb_variable;i++) {
      delete [] rank[i];
    }
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

    delete [] seq_mean;
    if (index_param_type == IMPLICIT_TYPE) {
      delete [] seq_index_parameter;
    }
  }

  delete [] change_point;

  return oseq;
}


/*--------------------------------------------------------------*
 *
 *  Computation of the contrast functions within a forward recursion.
 *
 *  arguments: time instant, sequence index, segment model types,
 *             flag contrast functions common to the individuals,
 *             log factorials for Poisson models, negative binomial distribution parameters,
 *             log binomial coefficients for negative binomial models,
 *             index parameters, sequence means for Gaussian change in the variance models.
 *             hyperparameters for Bayesian models, ranks for ordinal variables, contrasts,
 *             number of segments for bounding time loops.
 *
 *--------------------------------------------------------------*/

void Sequences::forward_contrast(int time , int index , segment_model *model_type , bool common_contrast ,
                                 double ***factorial , double *shape_parameter , double ***binomial_coeff ,
                                 int *seq_index_parameter , double *seq_mean , double **hyperparam ,
                                 double **rank , long double *contrast , int nb_segment) const

{
  register int i , j , k , m;
  int max_nb_value , count , *frequency , *inf_bound_parameter;
  double sum , factorial_sum , proba , binomial_coeff_sum , diff , index_parameter_sum ,
         index_parameter_diff , buff;
  long double index_parameter_square_sum , square_sum , mix_square_sum , prior_contrast , **residual;


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
        (((model_type[i - 1] == GAUSSIAN_CHANGE) || (model_type[i - 1] == LINEAR_MODEL_CHANGE) ||
          (model_type[i - 1] == VARIANCE_CHANGE) || (model_type[i - 1] == ORDINAL_GAUSSIAN_CHANGE)) && (!residual))) {
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

    else if (model_type[i - 1] == VARIANCE_CHANGE) {
      square_sum = 0.;

      if (type[i] != REAL_VALUE) {
        for (j = time;j >= nb_segment;j--) {
          diff = int_sequence[index][i][j] - seq_mean[i];
          square_sum += diff * diff;
          residual[index][j] = square_sum;
        }
      }

      else {
        for (j = time;j >= nb_segment;j--) {
          diff = real_sequence[index][i][j] - seq_mean[i];
          square_sum += diff * diff;
          residual[index][j] = square_sum;
        }
      }
    }

    else if (model_type[i - 1] == ORDINAL_GAUSSIAN_CHANGE) {
      square_sum = 0.;
      sum = rank[i][int_sequence[index][i][time]];
      residual[index][time] = 0.;

      for (j = time - 1;j >= nb_segment;j--) {
        diff = rank[i][int_sequence[index][i][j]] - sum / (time - j);
        square_sum += ((double)(time - j) / (double)(time - j + 1)) * diff * diff;
        sum += rank[i][int_sequence[index][i][j]];
        residual[index][j] = square_sum;
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

    else if ((model_type[i - 1] == GAUSSIAN_CHANGE) || (model_type[i - 1] == LINEAR_MODEL_CHANGE) ||
             (model_type[i - 1] == VARIANCE_CHANGE) || (model_type[i - 1] == ORDINAL_GAUSSIAN_CHANGE)) {
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


/*--------------------------------------------------------------*
 *
 *  Computation of the contrast functions within a backward recursion.
 *
 *  arguments: time instant, sequence index, segment model types,
 *             flag contrast functions common to the individuals,
 *             log factorials for Poisson models, negative binomial distribution parameters,
 *             log binomial coefficients for negative binomial models,
 *             index parameters, sequence means for Gaussian change in the variance models.
 *             hyperparameters for Bayesian models, ranks for ordinal variables, contrasts.
 *
 *--------------------------------------------------------------*/

void Sequences::backward_contrast(int time , int index , segment_model *model_type , bool common_contrast ,
                                  double ***factorial , double *shape_parameter , double ***binomial_coeff ,
                                  int *seq_index_parameter , double *seq_mean , double **hyperparam ,
                                  double **rank , long double *contrast) const

{
  register int i , j , k , m;
  int max_nb_value , count , *frequency , *inf_bound_parameter;
  double sum , factorial_sum , proba , binomial_coeff_sum , diff , index_parameter_sum ,
         index_parameter_diff , buff;
  long double index_parameter_square_sum , square_sum , mix_square_sum , prior_contrast , **residual;


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
        (((model_type[i - 1] == GAUSSIAN_CHANGE) || (model_type[i - 1] == LINEAR_MODEL_CHANGE) ||
          (model_type[i - 1] == VARIANCE_CHANGE) || (model_type[i - 1] == ORDINAL_GAUSSIAN_CHANGE)) && (!residual))) {
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

    else if (model_type[i - 1] == VARIANCE_CHANGE) {
      square_sum = 0.;

      if (type[i] != REAL_VALUE) {
        for (j = time;j < length[index];j++) {
          diff = int_sequence[index][i][j] - seq_mean[i];
          square_sum += diff * diff;
          residual[index][j] = square_sum;
        }
      }

      else {
        for (j = time;j < length[index];j++) {
          diff = real_sequence[index][i][j] - seq_mean[i];
          square_sum += diff * diff;
          residual[index][j] = square_sum;
        }
      }
    }

    else if (model_type[i - 1] == ORDINAL_GAUSSIAN_CHANGE) {
      square_sum = 0.;
      sum = rank[i][int_sequence[index][i][time]];
      residual[index][time] = 0.;

      for (j = time + 1;j < length[index];j++) {
        diff = rank[i][int_sequence[index][i][j]] - sum / (j - time);
        square_sum += ((double)(j - time) / (double)(j - time + 1)) * diff * diff;
        sum += rank[i][int_sequence[index][i][j]];
        residual[index][j] = square_sum;
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

    else if ((model_type[i - 1] == GAUSSIAN_CHANGE) || (model_type[i - 1] == LINEAR_MODEL_CHANGE) ||
             (model_type[i - 1] == VARIANCE_CHANGE) || (model_type[i - 1] == ORDINAL_GAUSSIAN_CHANGE)) {
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


/*--------------------------------------------------------------*
 *
 *  Optimal segmentation of a single sequence or a sample of sequences.
 *
 *  arguments: sequence index, number of segments, segment model types,
 *             flag contrast functions common to the individuals,
 *             negative binomial distribution parameters,
 *             ranks (ordinal variables), pointers on the segmentation likelihoods,
 *             on the number of free parameters of models and on the penalties
 *             related to segment lengths (for mBIC).
 *
 *--------------------------------------------------------------*/

double Sequences::segmentation(int index , int nb_segment , segment_model *model_type ,
                               bool common_contrast , double *shape_parameter , double **rank ,
                               double *isegmentation_likelihood , int *nb_parameter ,
                               double *segment_penalty)

{
  bool *used_output;
  register int i , j , k , m , n , p;
  int max_nb_value , seq_length , count , *inf_bound_parameter , *seq_index_parameter ,
      *psegment , **optimal_length;
  double buff , segmentation_likelihood , *seq_mean , **hyperparam , **forward ,
         ***factorial , ***binomial_coeff;
  long double *contrast; 


  max_nb_value = 0;
  factorial = new double**[nb_variable];
  inf_bound_parameter = new int[nb_variable];
  binomial_coeff = new double**[nb_variable];

  hyperparam = new double*[nb_variable];
  seq_index_parameter = NULL;

  for (i = 1;i < nb_variable;i++) {
    if ((model_type[i - 1] == CATEGORICAL_CHANGE) && (marginal_distribution[i]->nb_value > max_nb_value)) {
      max_nb_value = marginal_distribution[i]->nb_value;
    }

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

  seq_mean = new double[nb_variable];

  for (i = 1;i < nb_variable;i++) {
    if (model_type[i - 1] == VARIANCE_CHANGE) {
      seq_mean[i] = 0.;
      if (type[i] != REAL_VALUE) {
        for (j = 0;j < length[index];j++) {
          seq_mean[i] += int_sequence[index][i][j];
        }
      }
      else {
        for (j = 0;j < length[index];j++) {
          seq_mean[i] += real_sequence[index][i][j];
        }
      }
      seq_mean[i] /= length[index];
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
  }

  // forward recurrence

  for (i = 0;i < seq_length;i++) {

    // computation of segment contrast functions (log-likelihoods or sum of squared deviations)

    forward_contrast(i , index , model_type , common_contrast , factorial ,
                     shape_parameter , binomial_coeff , seq_index_parameter ,
                     seq_mean , hyperparam , rank , contrast);

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
              nb_parameter[i] += (i + 1) * nb_sequence;
            }
          }

          else if ((model_type[j - 1] == GAUSSIAN_CHANGE) || (model_type[j - 1] == ORDINAL_GAUSSIAN_CHANGE) ||
                   (model_type[j - 1] == BAYESIAN_GAUSSIAN_CHANGE)) {
            if ((index != I_DEFAULT) || (common_contrast)) {
              nb_parameter[i] += 2 * (i + 1);
            }
            else {
              nb_parameter[i] += 2 * (i + 1) * nb_sequence;
            }
          }

          else if (model_type[j - 1] == LINEAR_MODEL_CHANGE) {
            if ((index != I_DEFAULT) || (common_contrast)) {
              nb_parameter[i] += 3 * (i + 1);
            }
            else {
              nb_parameter[i] += 3 * (i + 1) * nb_sequence;
            }
          }

          else if (model_type[j - 1] == VARIANCE_CHANGE) {
            nb_parameter[i] += i + 2;
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
  }
  delete [] factorial;
  delete [] inf_bound_parameter;
  delete [] binomial_coeff;

  for (i = 1;i < nb_variable;i++) {
    delete [] hyperparam[i];
  }
  delete [] hyperparam;

  delete [] seq_mean;
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


/*--------------------------------------------------------------*
 *
 *  Optimal segmentation of a single sequence or a sample of sequences.
 *
 *  arguments: reference on a StatError object, stream, sequence identifier,
 *             number of segments, segment model types,
 *             flag contrast functions common to the individuals,
 *             negative binomial distribution parameters,
 *             output (sequence or residuals).
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::segmentation(StatError &error , ostream &os , int iidentifier ,
                                   int nb_segment , segment_model *model_type ,
                                   bool common_contrast , double *shape_parameter ,
                                   sequence_type output , bool continuity) const

{
  bool status = true;
  register int i , j;
  int index , nb_parameter , *psegment;
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
          ((output == SUBTRACTION_RESIDUAL) || (output == DIVISION_RESIDUAL))) {
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

#     ifdef MESSAGE
      psegment = seq->int_sequence[0][0] + 1;
      segment_penalty = 0.;
      i = 0;
      for (j = 1;j < seq->length[0];j++) {
        if (*psegment != *(psegment - 1)) {
          segment_penalty += log((double)(j - i));
          i = j;
        }
        psegment++;
      }
      segment_penalty += log((double)(seq->length[0] - i));

      nb_parameter = seq->nb_parameter_computation((index == I_DEFAULT ? index : 0) , nb_segment , model_type ,
                                                   common_contrast);

      penalized_likelihood = 2 * segmentation_likelihood - nb_parameter *
                             log((double)((seq->nb_variable - 1) * seq->length[0])) - segment_penalty;

      os << "\n" << nb_segment << " " << (nb_segment == 1 ? SEQ_label[SEQL_SEGMENT] : SEQ_label[SEQL_SEGMENTS])
         << "   2 * " << STAT_label[STATL_LIKELIHOOD] << ": " << 2 * segmentation_likelihood << "   "
         << nb_parameter << " " << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS]
         << "   2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " (Modified "  << STAT_criterion_word[BIC] << "): "
         << penalized_likelihood << endl;
#     endif

      oseq = seq->segmentation_output(nb_segment , model_type , common_contrast , os , output ,
                                      NULL , continuity);

      if (output == SEQUENCE) {
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


};  // namespace sequence_analysis
