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
 *       $Id: change_points1.cpp 18045 2015-04-23 09:33:56Z guedon $
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



/*--------------------------------------------------------------*
 *
 *  Calcul de la largeur d'une colonne de reels.
 *
 *  arguments : nombre de valeurs, pointeur sur des valeurs reelles.
 *
 *--------------------------------------------------------------*/

int column_width(int nb_value , const long double *value)

{
  register int i;
  int width , max_width = 0;


  for (i = 0;i < nb_value;i++) {
    ostringstream ostring;
    ostring << *value++;
    width = (ostring.str()).size();
    if (width > max_width) {
      max_width = width;
    }
  }

  return max_width;
}


/*--------------------------------------------------------------*
 *
 *  Calcul empirique des hyperparametres de la loi gamma a priori
 *  sur le parametre de la loi de Poisson.
 *
 *  arguments : indice de la sequence, indice de la variable,
 *              pointeurs sur les hyper-parametres.
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
 *  Calcul empirique des hyperparametres de la loi gaussiennne-gamma a priori
 *  sur les parametres de la loi gaussienne.
 *
 *  arguments : indice de la sequence, indice de la variable,
 *              pointeurs sur les hyper-parametres.
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
 *  Calcul du nombre de parametres independants.
 *
 *  arguments : indice de la sequence, nombre de segments, types des modeles.
 *
 *--------------------------------------------------------------*/

int Sequences::nb_parameter_computation(int index , int nb_segment , int *model_type) const

{
  bool *used_output;
  register int i , j , k;
  int nb_parameter , max_nb_value , *psegment;


  max_nb_value = 0;
  for (i = 1;i < nb_variable;i++) {
    if (((model_type[i - 1] == CATEGORICAL_CHANGE) || (model_type[0] == MULTIVARIATE_CATEGORICAL_CHANGE)) &&
        (marginal_distribution[i]->nb_value > max_nb_value)) {
      max_nb_value = marginal_distribution[i]->nb_value;
    }
  }

  if (max_nb_value > 0) {
    used_output = new bool[max_nb_value];
  }
  else {
    used_output = NULL;
  }

//  nb_parameter = 0;
  nb_parameter = nb_segment - 1;

  if (model_type[0] == MULTIVARIATE_CATEGORICAL_CHANGE) {
    psegment = int_sequence[index][0] + 1;

    for (i = 0;i < max_nb_value;i++) {
      used_output[i] = false;
    }

    used_output[int_sequence[index][1][0]] = true;
    for (i = 2;i < nb_variable;i++) {
      if (!used_output[int_sequence[index][i][0]]) {
        nb_parameter++;
        used_output[int_sequence[index][i][0]] = true;
      }
    }

    for (i = 1;i < length[index];i++) {
      if (*psegment != *(psegment - 1)) {
        for (j = 0;j < max_nb_value;j++) {
          used_output[j] = false;
        }

        used_output[int_sequence[index][1][i]] = true;
        for (j = 2;j < nb_variable;j++) {
          if (!used_output[int_sequence[index][j][i]]) {
            nb_parameter++;
            used_output[int_sequence[index][j][i]] = true;
          }
        }
      }

      else {
        for (j = 1;j < nb_variable;j++) {
          if (!used_output[int_sequence[index][j][i]]) {
            nb_parameter++;
            used_output[int_sequence[index][j][i]] = true;
          }
        }
      }

      psegment++;
    }
  }

  else if ((model_type[0] == MULTIVARIATE_POISSON_CHANGE) || (model_type[0] == MULTIVARIATE_GEOMETRIC_0_CHANGE)) {
    nb_parameter += nb_segment;
  }
  else if (model_type[0] == MEAN_CHANGE) {
    nb_parameter += (nb_variable - 1) * nb_segment + 1;
  }
  else if (model_type[0] == MEAN_VARIANCE_CHANGE) {
    nb_parameter += nb_variable * nb_segment;
  }
  else if (model_type[0] == INTERCEPT_SLOPE_CHANGE) {
    nb_parameter += (nb_variable - 1) * nb_segment * 2 + 1;
  }

  else {
    for (i = 1;i < nb_variable;i++) {
      if (model_type[i - 1] == CATEGORICAL_CHANGE) {
        psegment = int_sequence[index][0] + 1;

        for (j = 0;j < marginal_distribution[i]->nb_value;j++) {
          used_output[j] = false;
        }
        used_output[int_sequence[index][i][0]] = true;

        for (j = 1;j < length[index];j++) {
          if (*psegment != *(psegment - 1)) {
            for (k = 0;k < marginal_distribution[i]->nb_value;k++) {
              used_output[k] = false;
            }
            used_output[int_sequence[index][i][j]] = true;
          }

          else if (!used_output[int_sequence[index][i][j]]) {
            nb_parameter++;
            used_output[int_sequence[index][i][j]] = true;
          }

          psegment++;
        }
      }

      else if ((model_type[i - 1] == POISSON_CHANGE) || (model_type[i - 1] == GEOMETRIC_0_CHANGE) ||
               (model_type[i - 1] == GEOMETRIC_1_CHANGE) || (model_type[i - 1] == BAYESIAN_POISSON_CHANGE)) {
        nb_parameter += nb_segment;
      }

      else if (model_type[i - 1] == VARIANCE_CHANGE) {
        nb_parameter += nb_segment + 1;
      }

      else if (model_type[i - 1] == LINEAR_MODEL_CHANGE) {
        nb_parameter += 3 * nb_segment;
      }

      else {
        nb_parameter += 2 * nb_segment;
      }
    }
  }

  delete [] used_output;

  return nb_parameter;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la vraisemblance pour un segment.
 *
 *  arguments : indice de la sequence, types des modeles, rangs (variables ordinales).
 *
 *--------------------------------------------------------------*/

double Sequences::one_segment_likelihood(int index , int *model_type , double **rank) const

{
  register int i , j , k;
  int max_nb_value , *frequency , *seq_index_parameter = NULL , *psegment;
  double sum , factorial_sum , proba , diff , index_parameter_sum , index_parameter_diff , likelihood;
  long double square_sum , index_parameter_square_sum , mix_square_sum , residual;


  max_nb_value = 0;
  for (i = 1;i < nb_variable;i++) {
    if (((model_type[i - 1] == CATEGORICAL_CHANGE) || (model_type[0] == MULTIVARIATE_CATEGORICAL_CHANGE)) &&
        (marginal_distribution[i]->nb_value > max_nb_value)) {
      max_nb_value = marginal_distribution[i]->nb_value;
    }
  }

  if (max_nb_value > 0) {
    frequency = new int[max_nb_value];
  }
  else {
    frequency = NULL;
  }

  for (i = 1;i < nb_variable;i++) {
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

  if ((model_type[0] == MEAN_CHANGE) || (model_type[0] == MEAN_VARIANCE_CHANGE) ||
      (model_type[0] == INTERCEPT_SLOPE_CHANGE)) {
    residual = 0.;
  }
  else {
    likelihood = 0.;
  }

  if (model_type[0] == MULTIVARIATE_CATEGORICAL_CHANGE) {
    for (i = 0;i < max_nb_value;i++) {
      frequency[i] = 0;
    }

    for (i = 0;i < length[index];i++) {
      for (j = 1;j < nb_variable;j++) {
        frequency[int_sequence[index][j][i]]++;
      }
    }

    for (i = 0;i < max_nb_value;i++) {
      if (frequency[i] > 0) {
        likelihood += frequency[i] * log((double)frequency[i] / (double)((nb_variable - 1) * length[index]));
      }
    }
  }

  else if (model_type[0] == MULTIVARIATE_POISSON_CHANGE) {
    sum = 0.;
    factorial_sum = 0.;
    for (i = 0;i < length[index];i++) {
      for (j = 1;j < nb_variable;j++) {
        sum += int_sequence[index][j][i];
        for (k = 2;k <= int_sequence[index][j][i];k++) {
          factorial_sum += log((double)k);
        }
      }
    }

    if (sum > 0.) {
      likelihood = sum * (log(sum / ((nb_variable - 1) * length[index])) - 1) - factorial_sum;
    }
    else {
      likelihood = 0.;
    }
  }

  else if (model_type[0] == MULTIVARIATE_GEOMETRIC_0_CHANGE) {
    sum = 0.;
    for (i = 0;i < length[index];i++) {
      for (j = 1;j < nb_variable;j++) {
        sum += int_sequence[index][j][i];
      }
    }

    if (sum > 0.) {
      proba = (nb_variable - 1) * length[index] / ((nb_variable - 1) * length[index] + sum);
      likelihood = (nb_variable - 1) * length[index] * log(proba) + sum * log(1. - proba);
    }
    else {
      likelihood = 0.;
    }
  }

  else {
    for (i = 1;i < nb_variable;i++) {
      if (model_type[i - 1] == CATEGORICAL_CHANGE) {
        for (j = 0;j < marginal_distribution[i]->nb_value;j++) {
          frequency[j] = 0;
        }
        for (j = 0;j < length[index];j++) {
          frequency[int_sequence[index][i][j]]++;
        }

        for (j = 0;j < marginal_distribution[i]->nb_value;j++) {
          if (frequency[j] > 0) {
            likelihood += frequency[j] * log((double)frequency[j] / (double)(length[index]));
          }
        }
      }

      else if (model_type[i - 1] == POISSON_CHANGE) {
        sum = 0.;
        factorial_sum = 0.;
        for (j = 0;j < length[index];j++) {
          sum += int_sequence[index][i][j];
          for (k = 2;k <= int_sequence[index][i][j];k++) {
            factorial_sum += log((double)k);
          }
        }

        if (sum > 0.) {
          likelihood += sum * (log(sum / (length[index])) - 1) - factorial_sum;
        }
      }

      else if (model_type[i - 1] == GEOMETRIC_0_CHANGE) {
        sum = 0.;
        for (j = 0;j < length[index];j++) {
          sum += int_sequence[index][i][j];
        }

        if (sum > 0.) {
          proba = length[index] / (length[index] + sum);
          likelihood += length[index] * log(proba) + sum * log(1. - proba);
        }
      }

      else if (model_type[i - 1] == GEOMETRIC_1_CHANGE) {
        sum = 0.;
        for (j = 0;j < length[index];j++) {
          sum += int_sequence[index][i][j];
        }

        if (sum > 0.) {
          proba = length[index] / sum;
          likelihood += length[index] * log(proba) + (sum - length[index]) * log(1. - proba);
        }
      }

      else {
        if ((model_type[i - 1] != MEAN_CHANGE) && (model_type[i - 1] != MEAN_VARIANCE_CHANGE) &&
            (model_type[i - 1] != INTERCEPT_SLOPE_CHANGE)) {
          residual = 0.;
        }

        if (model_type[i - 1] == ORDINAL_GAUSSIAN_CHANGE) {
          sum = rank[i][int_sequence[index][i][0]];
          for (j = 1;j < length[index];j++) {
            diff = rank[i][int_sequence[index][i][j]] - sum / j;
            residual += ((double)j / (double)(j + 1)) * diff * diff;
            sum += rank[i][int_sequence[index][i][j]];
          }
        }

        else if ((model_type[i - 1] == LINEAR_MODEL_CHANGE) || (model_type[0] == INTERCEPT_SLOPE_CHANGE)) {
          if (type[i] != REAL_VALUE) {
            square_sum = 0.;
            index_parameter_square_sum = 0.;
            mix_square_sum = 0.;
            sum = int_sequence[index][i][0];
            index_parameter_sum = seq_index_parameter[0];

            for (j = 1;j < length[index];j++) {
              diff = int_sequence[index][i][j] - sum / j;
              square_sum += ((double)j / (double)(j + 1)) * diff * diff;
              index_parameter_diff = seq_index_parameter[j] - index_parameter_sum / j;
              index_parameter_square_sum += ((double)j / (double)(j + 1)) *
                                            index_parameter_diff * index_parameter_diff;
              mix_square_sum += ((double)j / (double)(j + 1)) * diff * index_parameter_diff;
              sum += int_sequence[index][i][j];
              index_parameter_sum += seq_index_parameter[j];
            }

            residual = square_sum - mix_square_sum * mix_square_sum / index_parameter_square_sum;
          }

          else {
            square_sum = 0.;
            index_parameter_square_sum = 0.;
            mix_square_sum = 0.;
            sum = real_sequence[index][i][0];
            index_parameter_sum = seq_index_parameter[0];

            for (j = 1;j < length[index];j++) {
              diff = real_sequence[index][i][j] - sum / j;
              square_sum += ((double)j / (double)(j + 1)) * diff * diff;
              index_parameter_diff = seq_index_parameter[j] - index_parameter_sum / j;
              index_parameter_square_sum += ((double)j / (double)(j + 1)) *
                                            index_parameter_diff * index_parameter_diff;
              mix_square_sum += ((double)j / (double)(j + 1)) * diff * index_parameter_diff;
              sum += real_sequence[index][i][j];
              index_parameter_sum += seq_index_parameter[j];
            }

            residual = square_sum - mix_square_sum * mix_square_sum / index_parameter_square_sum;
          }
        }

        else {
          if (type[i] != REAL_VALUE) {
            sum = int_sequence[index][i][0];
            for (j = 1;j < length[index];j++) {
              diff = int_sequence[index][i][j] - sum / j;
              residual += ((double)j / (double)(j + 1)) * diff * diff;
              sum += int_sequence[index][i][j];
            }

#           ifdef MESSAGE
            if ((model_type[i - 1] != MEAN_CHANGE) && (model_type[i - 1] != MEAN_VARIANCE_CHANGE) &&
                (model_type[i - 1] != INTERCEPT_SLOPE_CHANGE)) {
              double mean , diff;
              long double residual2;

              mean = 0.;
              for (j = 0;j < length[index];j++) {
                mean += int_sequence[index][i][j];
              }
              mean /= length[index];

              residual2 = 0.;
              for (j = 0;j < length[index];j++) {
                diff = int_sequence[index][i][j] - mean;
                residual2 += diff * diff;
              }

              if ((residual < residual2 - DOUBLE_ERROR) || (residual > residual2 + DOUBLE_ERROR)) {
                cout << "\nERROR: " << residual << " " << residual2 << endl;
              }
            }
#           endif

          }

          else {
            sum = real_sequence[index][i][0];
            for (j = 1;j < length[index];j++) {
              diff = real_sequence[index][i][j] - sum / j;
              residual += ((double)j / (double)(j + 1)) * diff * diff;
              sum += real_sequence[index][i][j];
            }
          }
        }

        if ((model_type[i - 1] != MEAN_CHANGE) && (model_type[i - 1] != MEAN_VARIANCE_CHANGE) &&
            (model_type[i - 1] != INTERCEPT_SLOPE_CHANGE)) {
//          if (residual > 0.) {
          if (residual > sqrt((double)length[index]) * ROUNDOFF_ERROR) {
            likelihood -= ((double)length[index] / 2.) * (logl(residual / length[index]) +
                           log(2 * M_PI) + 1);
/*            likelihood -= ((double)length[index] / 2.) * (logl(residual / (length[index] - 1)) +
                           log(2 * M_PI)) - (double)(length[index] - 1) / 2.; */
          }
          else {
            likelihood = D_INF;
            break;
          }
        }
      }
    }
  }

  if ((model_type[0] == MEAN_CHANGE) || (model_type[0] == MEAN_VARIANCE_CHANGE) ||
      (model_type[0] == INTERCEPT_SLOPE_CHANGE)) {
//    if (residual > 0.) {
    if (residual > sqrt((double)((nb_variable - 1) * length[index])) * ROUNDOFF_ERROR) {
      likelihood = -((double)((nb_variable - 1) * length[index]) / 2.) *
                    (logl(residual / ((nb_variable - 1) * length[index])) + log(2 * M_PI) + 1);
/*      likelihood = -((double)((nb_variable - 1) * length[index]) / 2.) *
                    (logl(residual / ((nb_variable - 1) * (length[index] - 1))) +
                     log(2 * M_PI)) - (double)((nb_variable - 1) * (length[index] - 1)) / 2.; */
    }
    else {
      likelihood = D_INF;
    }
  }

  for (i = 0;i < length[index];i++) {
    int_sequence[index][0][i] = 0;
  }

  delete [] frequency;

  if (index_parameter_type == IMPLICIT_TYPE) {
    delete [] seq_index_parameter;
  }

  return likelihood;
}


/*--------------------------------------------------------------*
 *
 *  Sortie des segmentations de sequences.
 *
 *  arguments : nombres de segments, types des modeles, stream,
 *              sortie (sequence ou residus), cas 1 sequence : ruptures.
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::segmentation_output(int *nb_segment , int *model_type , ostream &os ,
                                          int output , int *ichange_point)

{
  bool *piecewise_function;
  register int i , j , k , m , n;
  int max_nb_segment , *change_point , *psegment , *seq_index_parameter = NULL;
  double diff , response_mean , index_parameter_mean , response_variance , index_parameter_variance ,
         covariance , residual_mean , residual_square_sum , *change_point_amplitude , *global_variance ,
         *segment_variance , **mean , **variance , **intercept , **slope , **slope_standard_deviation ,
         **correlation;
  Test *test;
  Sequences *seq;


  max_nb_segment = nb_segment[0];
  for (i = 1;i < nb_sequence;i++) {
    if (nb_segment[i] > max_nb_segment) {
      max_nb_segment = nb_segment[i];
    }
  }

  if (ichange_point) {
    change_point = ichange_point;
  }
  else {
    change_point = new int[max_nb_segment + 1];
    change_point[0] = 0;
  }

  mean = new double*[nb_variable];
  for (i = 1;i < nb_variable;i++) {
    if ((model_type[i - 1] == POISSON_CHANGE) || (model_type[0] == MULTIVARIATE_POISSON_CHANGE) ||
        (model_type[i - 1] == GEOMETRIC_0_CHANGE) || (model_type[i - 1] == GEOMETRIC_1_CHANGE) ||
        (model_type[0] == MULTIVARIATE_GEOMETRIC_0_CHANGE) || (model_type[i - 1] == GAUSSIAN_CHANGE) ||
        (model_type[i - 1] == VARIANCE_CHANGE) || (model_type[0] == MEAN_CHANGE) ||
        (model_type[0] == MEAN_VARIANCE_CHANGE) || (model_type[i - 1] == BAYESIAN_POISSON_CHANGE) ||
        (model_type[i - 1] == BAYESIAN_GAUSSIAN_CHANGE)) {
      mean[i] = new double[max_nb_segment];
    }
    else {
      mean[i] = NULL;
    }
  }

  if (nb_sequence == 1) {
    variance = new double*[nb_variable];
    intercept = new double*[nb_variable];
    slope = new double*[nb_variable];
    slope_standard_deviation = new double*[nb_variable];
    correlation = new double*[nb_variable];

    for (i = 1;i < nb_variable;i++) {
      if ((model_type[i - 1] == POISSON_CHANGE) || (model_type[0] == MULTIVARIATE_POISSON_CHANGE) ||
          (model_type[i - 1] == GAUSSIAN_CHANGE) || (model_type[i - 1] == VARIANCE_CHANGE) ||
          (model_type[i - 1] == LINEAR_MODEL_CHANGE) || (model_type[0] == MEAN_CHANGE) ||
          (model_type[0] == MEAN_VARIANCE_CHANGE) || (model_type[0] == INTERCEPT_SLOPE_CHANGE) ||
          (model_type[i - 1] == BAYESIAN_POISSON_CHANGE) || (model_type[i - 1] == BAYESIAN_GAUSSIAN_CHANGE)) {
        variance[i] = new double[nb_segment[0]];
      }
      else {
        variance[i] = NULL;
      }

      if ((model_type[i - 1] == LINEAR_MODEL_CHANGE) || (model_type[0] == INTERCEPT_SLOPE_CHANGE)) {
        if (!seq_index_parameter) {
          if (index_parameter_type == IMPLICIT_TYPE) {
            seq_index_parameter = new int[length[0]];
            for (j = 0;j < length[0];j++) {
              seq_index_parameter[j] = j;
            }
          }
          else {
            seq_index_parameter = index_parameter[0];
          }
        }

        intercept[i] = new double[nb_segment[0]];
        slope[i] = new double[nb_segment[0]];
        slope_standard_deviation[i] = new double[nb_segment[0]];
        correlation[i] = new double[nb_segment[0]];
      }

      else {
        intercept[i] = NULL;
        slope[i] = NULL;
        slope_standard_deviation[i] = NULL;
        correlation[i] = NULL;
      }
    }

    change_point_amplitude = new double[nb_variable];
    global_variance = new double[nb_variable];

    if ((model_type[0] == MEAN_CHANGE) || (model_type[0] == INTERCEPT_SLOPE_CHANGE)) {
      global_variance[0] = 0.;
    }

    if (model_type[0] == MEAN_VARIANCE_CHANGE) {
      segment_variance = new double[nb_segment[0]];
      for (i = 0;i < nb_segment[0];i++) {
        segment_variance[i] = 0.;
      }
    }
    else {
      segment_variance = NULL;
    }
  }

  if (output == SEQUENCE) {
    piecewise_function = new bool[nb_variable];

    piecewise_function[0] = false;
    for (i = 1;i < nb_variable;i++) {
      if ((model_type[i - 1] == POISSON_CHANGE) || (model_type[0] == MULTIVARIATE_POISSON_CHANGE) ||
          (model_type[i - 1] == GEOMETRIC_0_CHANGE) || (model_type[i - 1] == GEOMETRIC_1_CHANGE) ||
          (model_type[0] == MULTIVARIATE_GEOMETRIC_0_CHANGE) || (model_type[i - 1] == GAUSSIAN_CHANGE) ||
          (model_type[i - 1] == LINEAR_MODEL_CHANGE) || (model_type[0] == MEAN_CHANGE) ||
          (model_type[0] == MEAN_VARIANCE_CHANGE) || (model_type[0] == INTERCEPT_SLOPE_CHANGE) ||
          (model_type[i - 1] == BAYESIAN_POISSON_CHANGE) || (model_type[i - 1] == BAYESIAN_GAUSSIAN_CHANGE)) {
        piecewise_function[i] = true;
      }
      else {
        piecewise_function[i] = false;
      }
    }

    seq = new Sequences(*this , piecewise_function);
  }

  for (i = 0;i < nb_sequence;i++) {
    if (!ichange_point) {
      psegment = int_sequence[i][0] + 1;
      j = 1;
      for (k = 1;k < length[i];k++) {
        if (*psegment != *(psegment - 1)) {
          change_point[j++] = k;
        }
        psegment++;
      }
      change_point[j] = length[i];
    }

    for (j = 1;j < nb_variable;j++) {
      if ((model_type[j - 1] == POISSON_CHANGE) || (model_type[0] == MULTIVARIATE_POISSON_CHANGE) ||
          (model_type[j - 1] == GEOMETRIC_0_CHANGE) || (model_type[j - 1] == GEOMETRIC_1_CHANGE) ||
          (model_type[0] == MULTIVARIATE_GEOMETRIC_0_CHANGE) || (model_type[j - 1] == GAUSSIAN_CHANGE) ||
          (model_type[j - 1] == VARIANCE_CHANGE) || (model_type[0] == MEAN_CHANGE) ||
          (model_type[0] == MEAN_VARIANCE_CHANGE) || (model_type[j - 1] == BAYESIAN_POISSON_CHANGE) ||
          (model_type[j - 1] == BAYESIAN_GAUSSIAN_CHANGE)) {
        if (type[j] != REAL_VALUE) {
          for (k = 0;k < nb_segment[i];k++) {
            mean[j][k] = 0.;
            for (m = change_point[k];m < change_point[k + 1];m++) {
              mean[j][k] += int_sequence[i][j][m];
            }
            mean[j][k] /= (change_point[k + 1] - change_point[k]);
          }
        }

        else {
          for (k = 0;k < nb_segment[i];k++) {
            mean[j][k] = 0.;
            for (m = change_point[k];m < change_point[k + 1];m++) {
              mean[j][k] += real_sequence[i][j][m];
            }
            mean[j][k] /= (change_point[k + 1] - change_point[k]);
          }
        }
      }

      if (nb_sequence == 1) {
        if ((model_type[j - 1] == GAUSSIAN_CHANGE) || (model_type[j - 1] == VARIANCE_CHANGE) ||
            (model_type[0] == MEAN_CHANGE) || (model_type[0] == MEAN_VARIANCE_CHANGE) ||
            (model_type[j - 1] == BAYESIAN_GAUSSIAN_CHANGE)) {
          if (nb_segment[i] > 1) {
            change_point_amplitude[j] = 0.;
            for (k = 1;k < nb_segment[i];k++) {
              change_point_amplitude[j] += fabs(mean[j][k] - mean[j][k - 1]);
            }
            change_point_amplitude[j] /= (nb_segment[i] - 1);
          }
        }

        if ((model_type[j - 1] == GAUSSIAN_CHANGE) || (model_type[j - 1] == BAYESIAN_GAUSSIAN_CHANGE)) {
          global_variance[j] = 0.;
        }

        if ((model_type[j - 1] == POISSON_CHANGE) || (model_type[0] == MULTIVARIATE_POISSON_CHANGE) ||
            (model_type[j - 1] == GAUSSIAN_CHANGE) || (model_type[j - 1] == VARIANCE_CHANGE) ||
            (model_type[0] == MEAN_CHANGE) || (model_type[0] == MEAN_VARIANCE_CHANGE) ||
            (model_type[j - 1] == BAYESIAN_POISSON_CHANGE) || (model_type[j - 1] == BAYESIAN_GAUSSIAN_CHANGE)) {
          if (type[j] != REAL_VALUE) {
            for (k = 0;k < nb_segment[i];k++) {
              variance[j][k] = 0.;
              if (change_point[k + 1] - change_point[k] > 1) {
                for (m = change_point[k];m < change_point[k + 1];m++) {
                  diff = int_sequence[i][j][m] - mean[j][k];
                  variance[j][k] += diff * diff;
                }

                if ((model_type[j - 1] == GAUSSIAN_CHANGE) || (model_type[j - 1] == BAYESIAN_GAUSSIAN_CHANGE)) {
                  global_variance[j] += variance[j][k];
                }
                else if (model_type[0] == MEAN_CHANGE) {
                  global_variance[0] += variance[j][k];
                }
                else if (model_type[0] == MEAN_VARIANCE_CHANGE) {
                  segment_variance[k] += variance[j][k];
                }

//                variance[j][k] /= (change_point[k + 1] - change_point[k]);
                variance[j][k] /= (change_point[k + 1] - change_point[k] - 1);
              }
            }
          }

          else {
            for (k = 0;k < nb_segment[i];k++) {
              variance[j][k] = 0.;
              if (change_point[k + 1] - change_point[k] > 1) {
                for (m = change_point[k];m < change_point[k + 1];m++) {
                  diff = real_sequence[i][j][m] - mean[j][k];
                  variance[j][k] += diff * diff;
                }

                if ((model_type[j - 1] == GAUSSIAN_CHANGE) || (model_type[j - 1] == BAYESIAN_GAUSSIAN_CHANGE)) {
                  global_variance[j] += variance[j][k];
                }
                else if (model_type[0] == MEAN_CHANGE) {
                  global_variance[0] += variance[j][k];
                }
                else if (model_type[0] == MEAN_VARIANCE_CHANGE) {
                  segment_variance[k] += variance[j][k];
                }

//                variance[j][k] /= (change_point[k + 1] - change_point[k]);
                variance[j][k] /= (change_point[k + 1] - change_point[k] - 1);
              }
            }
          }
        }

        if ((model_type[j - 1] == LINEAR_MODEL_CHANGE) || (model_type[0] == INTERCEPT_SLOPE_CHANGE)) {
          for (k = 0;k < nb_segment[i];k++) {
            response_mean = 0.;
            if (type[j] != REAL_VALUE) {
              for (m = change_point[k];m < change_point[k + 1];m++) {
                response_mean += int_sequence[i][j][m];
              }
            }
            else {
              for (m = change_point[k];m < change_point[k + 1];m++) {
                response_mean += real_sequence[i][j][m];
              }
            }
            response_mean /= (change_point[k + 1] - change_point[k]);

            index_parameter_mean = 0.;
            for (m = change_point[k];m < change_point[k + 1];m++) {
              index_parameter_mean += seq_index_parameter[m];                
            }
            index_parameter_mean /= (change_point[k + 1] - change_point[k]);

            response_variance = 0.;
            if (type[j] != REAL_VALUE) {
              for (m = change_point[k];m < change_point[k + 1];m++) {
                diff = int_sequence[i][j][m] - response_mean;
                response_variance += diff * diff;
              }
            }
            else {
              for (m = change_point[k];m < change_point[k + 1];m++) {
                diff = real_sequence[i][j][m] - response_mean;
                response_variance += diff * diff;
              }
            }

            index_parameter_variance = 0.;
            for (m = change_point[k];m < change_point[k + 1];m++) {
              diff = seq_index_parameter[m] - index_parameter_mean;
              index_parameter_variance += diff * diff;                
            }

            covariance = 0.;
            if (type[j] != REAL_VALUE) {
              for (m = change_point[k];m < change_point[k + 1];m++) {
                covariance += (int_sequence[i][j][m] - response_mean) * (seq_index_parameter[m] - index_parameter_mean);
              }
            }
            else {
              for (m = change_point[k];m < change_point[k + 1];m++) {
                covariance += (real_sequence[i][j][m] - response_mean) * (seq_index_parameter[m] - index_parameter_mean);
              }
            }

            slope[j][k] = covariance / index_parameter_variance;
            intercept[j][k] = response_mean - slope[j][k] * index_parameter_mean;
            correlation[j][k] = covariance / sqrt(response_variance * index_parameter_variance);

            residual_mean = 0.;
            residual_square_sum = 0.;
            if (type[j] != REAL_VALUE) {
              for (m = change_point[k];m < change_point[k + 1];m++) {
                diff = int_sequence[i][j][m] - (intercept[j][k] + slope[j][k] * seq_index_parameter[m]);
                residual_mean += diff;
                residual_square_sum += diff * diff;
              }
            }
            else {
              for (m = change_point[k];m < change_point[k + 1];m++) {
                diff = real_sequence[i][j][m] - (intercept[j][k] + slope[j][k] * seq_index_parameter[m]);
                residual_mean += diff;
                residual_square_sum += diff * diff;
              }
            }
            residual_mean /= (change_point[k + 1] - change_point[k]);
            if (change_point[k + 1] - change_point[k] > 2) {
              residual_square_sum /= (change_point[k + 1] - change_point[k] - 2);
            }
            else {
              residual_square_sum = 0.;
            }

            slope_standard_deviation[j][k] = sqrt(residual_square_sum / index_parameter_variance);

            variance[j][k] = 0.;
            if (change_point[k + 1] - change_point[k] > 2) {
              if (type[j] != REAL_VALUE) {
                for (m = change_point[k];m < change_point[k + 1];m++) {
                  diff = int_sequence[i][j][m] - (intercept[j][k] + slope[j][k] * seq_index_parameter[m]) - residual_mean;
                  variance[j][k] += diff * diff;
                }
              }
              else {
                for (m = change_point[k];m < change_point[k + 1];m++) {
                  diff = real_sequence[i][j][m] - (intercept[j][k] + slope[j][k] * seq_index_parameter[m]) - residual_mean;
                  variance[j][k] += diff * diff;
                }
              }

              if (model_type[0] == INTERCEPT_SLOPE_CHANGE) {
                global_variance[0] += variance[j][k];
              }

              variance[j][k] /= (change_point[k + 1] - change_point[k] - 2);
            }
          }
        }
      }
    }

#   ifdef MESSAGE
    if (nb_sequence == 1) {
      if (!ichange_point) {
        os << (nb_segment[i] == 2 ? SEQ_label[SEQL_CHANGE_POINT] : SEQ_label[SEQL_CHANGE_POINTS]) << ": ";

        if (index_parameter) {
          for (j = 1;j < nb_segment[i];j++) {
            os << index_parameter[i][change_point[j]];
            if (j < nb_segment[i] - 1) {
              os << ", ";
            }
          }
        }

        else {
          for (j = 1;j < nb_segment[i];j++) {
            os << change_point[j];
            if (j < nb_segment[i] - 1) {
              os << ", ";
            }
          }
        }
      }

      if ((index_interval) && (index_interval->variance > 0.)) {
        if (!ichange_point) {
          os << "   ";
        }
        os << SEQ_label[SEQL_SEGMENT_SAMPLE_SIZE] << ": ";
        for (j = 0;j < nb_segment[i];j++) {
          os << change_point[j + 1] - change_point[j];
          if (j < nb_segment[i] - 1) {
            os << ", ";
          }
        }
        os << endl;
      }

      else if (!ichange_point) {
        os << endl;
      }

      if ((model_type[0] == MEAN_CHANGE) || (model_type[0] == INTERCEPT_SLOPE_CHANGE)) {
        os << SEQ_label[SEQL_GLOBAL_STANDARD_DEVIATION] << ": ";
        if (model_type[0] == MEAN_CHANGE) {
          os << sqrt(global_variance[0] / ((nb_variable - 1) * (length[i] - nb_segment[i]))) << endl;
        }
        else if (model_type[0] == INTERCEPT_SLOPE_CHANGE) {
          os << sqrt(global_variance[0] / ((nb_variable - 1) * (length[i] - 2 * nb_segment[i]))) << endl;
        }
      }

      else if (model_type[0] == MEAN_VARIANCE_CHANGE) {
        os << SEQ_label[SEQL_SEGMENT] << " " << STAT_label[STATL_STANDARD_DEVIATION] << ": ";
        for (j = 0;j < nb_segment[i];j++) {
          if (change_point[j + 1] - change_point[j] > 1) {
            os << sqrt(segment_variance[j] / ((nb_variable - 1) *
                        (change_point[j + 1] - change_point[j] - 1))) << ", ";
          }
        }
        os << endl;
      }

      if (nb_variable > 2) {
        os << "\n";
      }
      for (j = 1;j < nb_variable;j++) {
        if ((model_type[j - 1] == POISSON_CHANGE) || (model_type[0] == MULTIVARIATE_POISSON_CHANGE) ||
            (model_type[j - 1] == BAYESIAN_POISSON_CHANGE)) {
          if (nb_variable > 2) {
            os << STAT_label[STATL_VARIABLE] << " " << j << "   ";
          }
          os << SEQ_label[SEQL_SEGMENT] << " " << STAT_label[STATL_MEAN] << ", "
             << STAT_label[STATL_VARIANCE] << ": ";
          for (k = 0;k < nb_segment[i];k++) {
            os << mean[j][k] << " " << variance[j][k];
            if (k < nb_segment[i] - 1) {
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
          for (k = 0;k < nb_segment[i];k++) {
            os << mean[j][k] << " ";
            if (k < nb_segment[i] - 1) {
              os << " | ";
            }
          }
          os << endl;
        }

        else if ((model_type[j - 1] == GAUSSIAN_CHANGE) || (model_type[0] == VARIANCE_CHANGE) ||
                 (model_type[0] == MEAN_CHANGE) || (model_type[0] == MEAN_VARIANCE_CHANGE) ||
                 (model_type[j - 1] == BAYESIAN_GAUSSIAN_CHANGE)) {
          if (nb_variable > 2) {
            os << STAT_label[STATL_VARIABLE] << " " << j << "   ";
          }
          os << SEQ_label[SEQL_SEGMENT] << " " << STAT_label[STATL_MEAN] << ", "
             << STAT_label[STATL_STANDARD_DEVIATION] << ": ";
          for (k = 0;k < nb_segment[i];k++) {
            os << mean[j][k] << " " << sqrt(variance[j][k]);
            if (k < nb_segment[i] - 1) {
              os << " | ";
            }
          }
          os << endl;

          if (nb_segment[i] > 1) {
            os << SEQ_label[SEQL_MEAN_CHANGE_POINT_AMPLITUDE] << ": " << change_point_amplitude[j] << "   ";
          }
          if ((model_type[j - 1] == GAUSSIAN_CHANGE) || (model_type[j - 1] == BAYESIAN_GAUSSIAN_CHANGE)) {
            os << SEQ_label[SEQL_GLOBAL_STANDARD_DEVIATION] << ": " << sqrt(global_variance[j] / (length[i] - nb_segment[i]));
            if (nb_segment[i] > 1) {
              os << "   " << STAT_label[STATL_RATIO] << ": "
                 << change_point_amplitude[j] / sqrt(global_variance[j] / (length[i] - nb_segment[i]));
            }
          }
          os << endl;
        }

        else if ((model_type[j - 1] == LINEAR_MODEL_CHANGE) || (model_type[0] == INTERCEPT_SLOPE_CHANGE)) {
          if (nb_variable > 2) {
            os << STAT_label[STATL_VARIABLE] << " " << j << "   ";
          }
          os << SEQ_label[SEQL_SEGMENT] << " " << STAT_label[STATL_INTERCEPT] << ", "
             << STAT_label[STATL_SLOPE] << ", " << STAT_label[STATL_CORRELATION_COEFF] << " ("
             << STAT_label[STATL_LIMIT_CORRELATION_COEFF] << "), "
             << STAT_label[STATL_RESIDUAL] << " " << STAT_label[STATL_STANDARD_DEVIATION] << ": ";
          for (k = 0;k < nb_segment[i];k++) {
            test = new Test(STUDENT , false , change_point[k + 1] - change_point[k] - 2 , I_DEFAULT , D_DEFAULT);
            test->critical_probability = ref_critical_probability[0];
            test->t_value_computation();

            os << intercept[j][k] << " " << slope[j][k];
            if (slope_standard_deviation[j][k] > 0.) {
              os << " (" << slope[j][k] - test->value * slope_standard_deviation[j][k] << ", "
                 << slope[j][k] + test->value * slope_standard_deviation[j][k] << ") ";
//               << "slope_standard_deviation: " << slope_standard_deviation[j][k] << " "
            }
            os << correlation[j][k] << " (-/+"
               << test->value / sqrt(test->value * test->value + change_point[k + 1] - change_point[k] - 2) << ") ";
            os << sqrt(variance[j][k]);

            delete test;
            if (k < nb_segment[i] - 1) {
              os << " | ";
            }
          }
          os << endl;
        }
      }
    }
#   endif

    switch (output) {

    case SEQUENCE : {
      j = 1;
      for (k = 1;k < nb_variable;k++) {
        j++;
        if (piecewise_function[k]) {
          if ((model_type[k - 1] != LINEAR_MODEL_CHANGE) && (model_type[k - 1] != INTERCEPT_SLOPE_CHANGE)) {
            for (m = 0;m < nb_segment[i];m++) {
              for (n = change_point[m];n < change_point[m + 1];n++) {
                seq->real_sequence[i][j][n] = mean[k][m];
              }
            }
          }

          else {
            for (m = 0;m < nb_segment[i];m++) {
              for (n = change_point[m];n < change_point[m + 1];n++) {
                seq->real_sequence[i][j][n] = intercept[k][m] + slope[k][m] * seq_index_parameter[n];
              }
            }
          }
          j++;
        }
      }
      break;
    }

    // calcul des residus

    case SUBTRACTION_RESIDUAL : {
      for (j = 1;j < nb_variable;j++) {
        if (type[j] != REAL_VALUE) {
          real_sequence[i][j] = new double[length[i]];

          if ((model_type[j - 1] != LINEAR_MODEL_CHANGE) && (model_type[j - 1] != INTERCEPT_SLOPE_CHANGE)) {
            for (k = 0;k < nb_segment[i];k++) {
              for (m = change_point[k];m < change_point[k + 1];m++) {
                real_sequence[i][j][m] = int_sequence[i][j][m] - mean[j][k];
              }
            }
          }

          else {
            for (k = 0;k < nb_segment[i];k++) {
              for (m = change_point[k];m < change_point[k + 1];m++) {
                real_sequence[i][j][m] = int_sequence[i][j][m] - (intercept[j][k] + slope[j][k] * seq_index_parameter[m]);
              }
            }
          }

          delete [] int_sequence[i][j];
          int_sequence[i][j] = NULL;
        }

        else {
          if ((model_type[j - 1] != LINEAR_MODEL_CHANGE) && (model_type[j - 1] != INTERCEPT_SLOPE_CHANGE)) {
            for (k = 0;k < nb_segment[i];k++) {
              for (m = change_point[k];m < change_point[k + 1];m++) {
                real_sequence[i][j][m] -= mean[j][k];
              }
            }
          }

          else {
            for (k = 0;k < nb_segment[i];k++) {
              for (m = change_point[k];m < change_point[k + 1];m++) {
                real_sequence[i][j][m] -= (intercept[j][k] + slope[j][k] * seq_index_parameter[m]);
              }
            }
          }
        }
      }
      break;
    }

    case DIVISION_RESIDUAL : {
      for (j = 1;j < nb_variable;j++) {
        if (type[j] != REAL_VALUE) {
          real_sequence[i][j] = new double[length[i]];

          if ((model_type[j - 1] != LINEAR_MODEL_CHANGE) && (model_type[j - 1] != INTERCEPT_SLOPE_CHANGE)) {
            for (k = 0;k < nb_segment[i];k++) {
              if (mean[j][k] != 0.) {
                for (m = change_point[k];m < change_point[k + 1];m++) {
                  real_sequence[i][j][m] = int_sequence[i][j][m] / mean[j][k];
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
            for (k = 0;k < nb_segment[i];k++) {
              for (m = change_point[k];m < change_point[k + 1];m++) {
                if (intercept[j][k] + slope[j][k] * seq_index_parameter[m] != 0.) {
                  real_sequence[i][j][m] = int_sequence[i][j][m] / (intercept[j][k] + slope[j][k] * seq_index_parameter[m]);
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
            for (k = 0;k < nb_segment[i];k++) {
              if (mean[j][k] != 0.) {
                for (m = change_point[k];m < change_point[k + 1];m++) {
                  real_sequence[i][j][m] /= mean[j][k];
                }
              }
            }
          }

          else {
            for (k = 0;k < nb_segment[i];k++) {
              for (m = change_point[k];m < change_point[k + 1];m++) {
                if (intercept[j][k] + slope[j][k] * seq_index_parameter[m] != 0.) {
                  real_sequence[i][j][m] /= (intercept[j][k] + slope[j][k] * seq_index_parameter[m]);
                }
              }
            }
          }
        }
      }
      break;
    }
    }
  }

  if (!ichange_point) {
    delete [] change_point;
  }

  for (i = 1;i < nb_variable;i++) {
    delete [] mean[i];
  }
  delete [] mean;

  if (nb_sequence == 1) {
    for (i = 1;i < nb_variable;i++) {
      delete [] variance[i];
      delete [] intercept[i];
      delete [] slope[i];
      delete [] slope_standard_deviation[i];
      delete [] correlation[i];
    }
    delete [] variance;
    delete [] intercept;
    delete [] slope;
    delete [] slope_standard_deviation;
    delete [] correlation;

    if (index_parameter_type == IMPLICIT_TYPE) {
      delete [] seq_index_parameter;
    }

    delete [] change_point_amplitude;
    delete [] global_variance;

    delete [] segment_variance;
  }

  if (output == SEQUENCE) {
    delete [] piecewise_function;
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
 *  Segmentation d'une sequence.
 *
 *  arguments : reference sur un objet StatError, stream, identificateur de la sequence,
 *              nombre de segments, ruptures, types des modeles,
 *              sortie (sequence ou residus).
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::segmentation(StatError &error , ostream &os , int iidentifier ,
                                   int nb_segment , int *ichange_point , int *model_type ,
                                   int output) const

{
  bool status = true;
  register int i , j , k , m;
  int index , max_nb_value , nb_parameter , *change_point = NULL , *frequency ,
      *seq_index_parameter = NULL , *psegment , inb_segment[1];
  double sum , factorial_sum , proba , diff , index_parameter_sum , index_parameter_diff ,
         segmentation_likelihood , segment_penalty , penalized_likelihood , *mean ,
         *segment_variance , **rank;
  long double square_sum , index_parameter_square_sum , mix_square_sum , *residual;
  const FrequencyDistribution **pmarginal;
  FrequencyDistribution *marginal;
  Sequences *iseq , *seq , *oseq;


  oseq = NULL;
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
        (model_type[0] == MULTIVARIATE_GEOMETRIC_0_CHANGE) || (model_type[i] == ORDINAL_GAUSSIAN_CHANGE)) {
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
        if (((model_type[i] != GEOMETRIC_0_CHANGE) && (min_value[i] < 0)) ||
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

  else {
    if ((nb_segment < 1) || (nb_segment > length[index] / 2)) {
      status = false;
      error.update(SEQ_error[SEQR_NB_SEGMENT]);
    }

    else {
      change_point = new int[nb_segment + 1];

      if (index_parameter) {
        change_point[0] = index_parameter[index][0];
        change_point[nb_segment] = index_parameter[index][length[index] - 1] + 1;
      }
      else {
        change_point[0] = 0;
        change_point[nb_segment] = length[index];
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
        for (j = 1;j < length[index];j++) {
          if (index_parameter[index][j] == change_point[i]) {
            change_point[i++] = j;
          }
        }

        if (i < nb_segment - 1) {
          status = false;
          error.update(SEQ_error[SEQR_CHANGE_POINT]);
        }
        else {
          change_point[nb_segment] = length[index];
        }
      }
    }
  }

  if (status) {

    // calcul des rangs pour les variables ordinales

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

    for (i = 0;i < nb_variable;i++) {
      if (((model_type[i] == CATEGORICAL_CHANGE) || (model_type[0] == MULTIVARIATE_CATEGORICAL_CHANGE)) &&
          (marginal_distribution[i]->nb_value > max_nb_value)) {
        max_nb_value = marginal_distribution[i]->nb_value;
      }
    }

    if (max_nb_value > 0) {
      frequency = new int[max_nb_value];
    }
    else {
      frequency = NULL;
    }

    mean = new double[nb_variable];
    residual = new long double[nb_segment];

    if (model_type[0] == MEAN_VARIANCE_CHANGE) {
      segment_variance = new double[nb_segment];
    }
    else {
      segment_variance = NULL;
    }

    for (i = 0;i < nb_variable;i++) {
      if (model_type[i] == VARIANCE_CHANGE) {
        mean[i] = 0.;
        if (type[i] != REAL_VALUE) {
          for (j = 0;j < length[index];j++) {
            mean[i] += int_sequence[index][i][j];
          }
        }
        else {
          for (j = 0;j < length[index];j++) {
            mean[i] += real_sequence[index][i][j];
          }
        }
        mean[i] /= length[index];
      }

      if (((model_type[i] == LINEAR_MODEL_CHANGE) && (!seq_index_parameter)) ||
          ((i == 0) && (model_type[0] == INTERCEPT_SLOPE_CHANGE))) {
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

    if ((model_type[0] == MEAN_CHANGE) || (model_type[0] == INTERCEPT_SLOPE_CHANGE)) {
      for (i = 0;i < nb_segment;i++) {
        residual[i] = 0.;
      }
    }
    else {
      segmentation_likelihood = 0.;
    }

    if (model_type[0] == MULTIVARIATE_CATEGORICAL_CHANGE) {
      for (i = 0;i < nb_segment;i++) {
        for (j = 0;j < max_nb_value;j++) {
          frequency[j] = 0;
        }

        for (j = change_point[i];j < change_point[i + 1];j++) {
          for (k = 0;k < nb_variable;k++) {
            frequency[int_sequence[index][k][j]]++;
          }
        }

        for (j = 0;j < max_nb_value;j++) {
          if (frequency[j] > 0) {
            segmentation_likelihood += frequency[j] * log((double)frequency[j] /
                                        (double)(nb_variable * (change_point[i + 1] - change_point[i])));
          }
        }
      }
    }

    else if (model_type[0] == MULTIVARIATE_POISSON_CHANGE) {
      for (i = 0;i < nb_segment;i++) {
        sum = 0.;
        factorial_sum = 0.;
        for (j = change_point[i];j < change_point[i + 1];j++) {
          for (k = 0;k < nb_variable;k++) {
            sum += int_sequence[index][k][j];
            for (m = 2;m <= int_sequence[index][k][j];m++) {
              factorial_sum += log((double)m);
            }
          }
        }

        if (sum > 0.) {
          segmentation_likelihood += sum * (log(sum / (nb_variable * (change_point[i + 1] - change_point[i]))) - 1) -
                                     factorial_sum;
        }
      }
    }

    else if (model_type[0] == MULTIVARIATE_GEOMETRIC_0_CHANGE) {
      for (i = 0;i < nb_segment;i++) {
        sum = 0.;
        for (j = change_point[i];j < change_point[i + 1];j++) {
          for (k = 0;k < nb_variable;k++) {
            sum += int_sequence[index][k][j];
          }
        }

        if (sum > 0.) {
          proba = nb_variable * (change_point[i + 1] - change_point[i]) /
                  (nb_variable * (change_point[i + 1] - change_point[i]) + sum);
          segmentation_likelihood += nb_variable * (change_point[i + 1] - change_point[i]) * log(proba) +
                                     sum * log(1. - proba);
        }
      }
    }

    else {
      if (model_type[0] == MEAN_VARIANCE_CHANGE) {
        for (i = 0;i < nb_segment;i++) {
          residual[i] = 0.;
        }
      }

      for (i = 0;i < nb_variable;i++) {

        // calcul des log-vraisemblances des segments

        if (model_type[i] == CATEGORICAL_CHANGE) {
          for (j = 0;j < nb_segment;j++) {
            for (k = 0;k < marginal_distribution[i]->nb_value;k++) {
              frequency[k] = 0;
            }
            for (k = change_point[j];k < change_point[j + 1];k++) {
              frequency[int_sequence[index][i][k]]++;
            }

            for (k = 0;k < marginal_distribution[i]->nb_value;k++) {
              if (frequency[k] > 0) {
                segmentation_likelihood += frequency[k] * log((double)frequency[k] /
                                            (double)(change_point[j + 1] - change_point[j]));
              }
            }
          }
        }

        else if (model_type[i] == POISSON_CHANGE) {
          for (j = 0;j < nb_segment;j++) {
            sum = 0.;
            factorial_sum = 0.;
            for (k = change_point[j];k < change_point[j + 1];k++) {
              sum += int_sequence[index][i][k];
              for (m = 2;m <= int_sequence[index][i][k];m++) {
                factorial_sum += log((double)m);
              }
            }

            if (sum > 0.) {
              segmentation_likelihood += sum * (log(sum / (change_point[j + 1] - change_point[j])) - 1) -
                                         factorial_sum;
            }
          }
        }

        else if (model_type[i] == GEOMETRIC_0_CHANGE) {
          for (j = 0;j < nb_segment;j++) {
            sum = 0.;
            for (k = change_point[j];k < change_point[j + 1];k++) {
              sum += int_sequence[index][i][k];
            }

            if (sum > 0.) {
              proba = (change_point[j + 1] - change_point[j]) / (change_point[j + 1] - change_point[j] + sum);
              segmentation_likelihood += (change_point[j + 1] - change_point[j]) * log(proba) + sum * log(1. - proba);
            }
          }
        }

        else if (model_type[i] == GEOMETRIC_1_CHANGE) {
          for (j = 0;j < nb_segment;j++) {
            sum = 0.;
            for (k = change_point[j];k < change_point[j + 1];k++) {
              sum += int_sequence[index][i][k];
            }

            if (sum > 0.) {
              proba = (change_point[j + 1] - change_point[j]) / sum;
              segmentation_likelihood += (change_point[j + 1] - change_point[j]) * log(proba) +
                                         (sum - (change_point[j + 1] - change_point[j])) * log(1. - proba);
            }
          }
        }

        else {
          for (j = 0;j < nb_segment;j++) {
            if (model_type[i] == VARIANCE_CHANGE) {
              residual[j] = 0.;

              if (type[i] != REAL_VALUE) {
                for (k = change_point[j];k < change_point[j + 1];k++) {
                  diff = int_sequence[index][i][k] - mean[i];
                  residual[j] += diff * diff;
                }
              }

              else {
                for (k = change_point[j];k < change_point[j + 1];k++) {
                  diff = real_sequence[index][i][k] - mean[i];
                  residual[j] += diff * diff;
                }
              }
            }

            else if (model_type[i] == ORDINAL_GAUSSIAN_CHANGE) {
              residual[j] = 0.;
              sum = rank[i][int_sequence[index][i][change_point[j]]];

              for (k = change_point[j] + 1;k < change_point[j + 1];k++) {
                diff = rank[i][int_sequence[index][i][k]] - sum / (k - change_point[j]);
                residual[j] += ((double)(k - change_point[j]) / (double)(k - change_point[j] + 1)) *
                               diff * diff;
                sum += rank[i][int_sequence[index][i][k]];
              }
            }

            else if ((model_type[i] == LINEAR_MODEL_CHANGE) || (model_type[0] == INTERCEPT_SLOPE_CHANGE)) {
              if (type[i] != REAL_VALUE) {
                square_sum = 0.;
                index_parameter_square_sum = 0.;
                mix_square_sum = 0.;
                sum = int_sequence[index][i][change_point[j]];
                index_parameter_sum = seq_index_parameter[change_point[j]];

                for (k = change_point[j] + 1;k < change_point[j + 1];k++) {
                  diff = int_sequence[index][i][k] - sum / (k - change_point[j]);
                  square_sum += ((double)(k - change_point[j]) / (double)(k - change_point[j] + 1)) *
                                diff * diff;
                  index_parameter_diff = seq_index_parameter[k] - index_parameter_sum / (k - change_point[j]);
                  index_parameter_square_sum += ((double)(k - change_point[j]) / (double)(k - change_point[j] + 1)) *
                                                index_parameter_diff * index_parameter_diff;
                  mix_square_sum += ((double)(k - change_point[j]) / (double)(k - change_point[j] + 1)) *
                                    diff * index_parameter_diff;
                  sum += int_sequence[index][i][k];
                  index_parameter_sum += seq_index_parameter[k];
                }

                if ((change_point[j + 1] - change_point[j] > 2) && (index_parameter_square_sum > 0.)) {
                  residual[j] = square_sum - mix_square_sum * mix_square_sum / index_parameter_square_sum;
                }
                else {
                  residual[j] = 0.;
                }
              }

              else {
                square_sum = 0.;
                index_parameter_square_sum = 0.;
                mix_square_sum = 0.;
                sum = real_sequence[index][i][change_point[j]];
                index_parameter_sum = seq_index_parameter[change_point[j]];

                for (k = change_point[j] + 1;k < change_point[j + 1];k++) {
                  diff = real_sequence[index][i][k] - sum / (k - change_point[j]);
                  square_sum += ((double)(k - change_point[j]) / (double)(k - change_point[j] + 1)) *
                                 diff * diff;
                  index_parameter_diff = seq_index_parameter[k] - index_parameter_sum / (k - change_point[j]);
                  index_parameter_square_sum += ((double)(k - change_point[j]) / (double)(k - change_point[j] + 1)) *
                                                index_parameter_diff * index_parameter_diff;
                  mix_square_sum += ((double)(k - change_point[j]) / (double)(k - change_point[j] + 1)) *
                                    diff * index_parameter_diff;
                  sum += real_sequence[index][i][k];
                  index_parameter_sum += seq_index_parameter[k];
                }

                if ((change_point[j + 1] - change_point[j] > 2) && (index_parameter_square_sum > 0.)) {
                  residual[j] = square_sum - mix_square_sum * mix_square_sum / index_parameter_square_sum;
                }
                else {
                  residual[j] = 0.;
                }
              }
            }

            else {
              if (model_type[i] == GAUSSIAN_CHANGE) {
                residual[j] = 0.;
              }

              if (type[i] != REAL_VALUE) {
                sum = int_sequence[index][i][change_point[j]];
                for (k = change_point[j] + 1;k < change_point[j + 1];k++) {
                  diff = int_sequence[index][i][k] - sum / (k - change_point[j]);
                  residual[j] += ((double)(k - change_point[j]) / (double)(k - change_point[j] + 1)) *
                                 diff * diff;
                  sum += int_sequence[index][i][k];
                }
              }

              else {
                sum = real_sequence[index][i][change_point[j]];
                for (k = change_point[j] + 1;k < change_point[j + 1];k++) {
                  diff = real_sequence[index][i][k] - sum / (k - change_point[j]);
                  residual[j] += ((double)(k - change_point[j]) / (double)(k - change_point[j] + 1)) *
                                 diff * diff;
                  sum += real_sequence[index][i][k];
                }
              }
            }

            if ((model_type[i] != MEAN_CHANGE) && (model_type[i] != MEAN_VARIANCE_CHANGE) &&
                (model_type[i] != INTERCEPT_SLOPE_CHANGE)) {
//              if (residual[j] > 0.) {
              if (residual[j] > sqrt((double)(change_point[j + 1] - change_point[j])) * ROUNDOFF_ERROR) {
                segmentation_likelihood -= ((double)(change_point[j + 1] - change_point[j]) / 2.) * (logl(residual[j] /
                                             (change_point[j + 1] - change_point[j])) + log(2 * M_PI) + 1);
/*                segmentation_likelihood -= ((double)(change_point[j + 1] - change_point[j]) / 2.) * (logl(residual[j] /
                                             (change_point[j + 1] - change_point[j])) + log(2 * M_PI)) +
                                           (double)(change_point[j + 1] - change_point[j]) / 2.; */
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
        for (i = 1;i < nb_segment;i++) {
          residual[0] += residual[i];
        }
//         if (residual[0] > 0.) {
        if (residual[0] > sqrt((double)(nb_variable * length[index])) * ROUNDOFF_ERROR) {
          segmentation_likelihood = -((double)(nb_variable * length[index]) / 2.) *
                                     (logl(residual[0] / (nb_variable * length[index])) + log(2 * M_PI) + 1);
/*          segmentation_likelihood = -((double)(nb_variable * length[index]) / 2.) *
                                     (logl(residual[0] / (nb_variable * (length[index] - nb_segment))) + log(2 * M_PI)) -
                                     (double)(nb_variable * (length[index] - nb_segment)) / 2.; */
        }
        else {
          segmentation_likelihood = D_INF;
        }
      }

      else if (model_type[0] == MEAN_VARIANCE_CHANGE) {
        for (i = 0;i < nb_segment;i++) {
//          if (residual[i] > 0.) {
          if (residual[i] > sqrt((double)(nb_variable * (change_point[i + 1] - change_point[i]))) * ROUNDOFF_ERROR) {
            segmentation_likelihood -= ((double)(nb_variable * (change_point[i + 1] - change_point[i])) / 2.) *
                                       (logl(residual[i] / (nb_variable * (change_point[i + 1] - change_point[i]))) +
                                        log(2 * M_PI) + 1);
/*            segmentation_likelihood -= ((double)(nb_variable * (change_point[i + 1] - change_point[i])) / 2.) *
                                       (logl(residual[i] / (nb_variable * (change_point[i + 1] - change_point[i]))) +
                                        log(2 * M_PI)) + (double)(nb_variable * (change_point[i + 1] - change_point[i])) / 2.; */
            segment_variance[i] = residual[i] / (nb_variable * (change_point[i + 1] - change_point[i]));
          }
          else {
            segmentation_likelihood = D_INF;
            break;
          }
        }
      }
    }

    iseq = new Sequences(*this , 1 , &index);
    seq = new Sequences(*iseq , 'a');
    delete iseq;

    psegment = seq->int_sequence[0][0];
    for (i = 0;i < nb_segment;i++) {
      for (j = change_point[i];j < change_point[i + 1];j++) {
        *psegment++ = i;
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

      nb_parameter = seq->nb_parameter_computation(0 , nb_segment , model_type);

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

    inb_segment[0] = nb_segment;
    oseq = seq->segmentation_output(inb_segment , model_type , os , output , change_point);

    if (output == SEQUENCE) {
      delete seq;
    }

    for (i = 0;i < nb_variable;i++) {
      delete [] rank[i];
    }
    delete [] rank;

    delete [] frequency;

    delete [] mean;
    delete [] residual;

    delete [] segment_variance;

    if (index_parameter_type == IMPLICIT_TYPE) {
      delete [] seq_index_parameter;
    }
  }

  delete [] change_point;

  return oseq;
}


/*--------------------------------------------------------------*
 *
 *  Segmentation optimale de sequences.
 *
 *  arguments : nombres de segments, rang (variable ordinale), types des modeles,
 *              cas 1 sequence : pointeurs sur les vraisemblances des segmentations,
 *              sur les nombre de parametres des modeles et sur les penalites
 *              liees aux longueurs de segments.
 *
 *--------------------------------------------------------------*/

double Sequences::segmentation(int *nb_segment , int *model_type , double **rank ,
                               double *isegmentation_likelihood , int *nb_parameter ,
                               double *segment_penalty)

{
  bool *used_output;
  register int i , j , k , m , n , r;
  int max_nb_segment , max_nb_value , *frequency , *seq_index_parameter = NULL , *psegment ,
      **optimal_length;
  double sum , factorial_sum , proba , diff , index_parameter_sum , index_parameter_diff , buff ,
         segmentation_likelihood , *mean , **hyperparam , **factorial , **forward;
  long double square_sum , index_parameter_square_sum , mix_square_sum , prior_contrast ,
       *residual , *contrast;


  max_nb_segment = nb_segment[0];
  for (i = 1;i < nb_sequence;i++) {
    if (nb_segment[i] > max_nb_segment) {
      max_nb_segment = nb_segment[i];
    }
  }

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
      factorial[i] = new double[max_length];
    }
    else {
      factorial[i] = NULL;
    }

    if (model_type[i - 1] == BAYESIAN_POISSON_CHANGE) {
      hyperparam[i] = new double[2];
    }
    else if (model_type[i - 1] == BAYESIAN_GAUSSIAN_CHANGE) {
      hyperparam[i] = new double[4];
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

  mean = new double[nb_variable];
  residual = new long double[max_length];

  contrast = new long double[max_length];

  forward = new double*[max_length];
  for (i = 0;i < max_length;i++) {
    forward[i] = new double[max_nb_segment];
  }

  optimal_length = new int*[max_length];
  for (i = 0;i < max_length;i++) {
    optimal_length[i] = new int[max_nb_segment];
  }

  if ((nb_parameter) && (max_nb_value > 0)) {
    used_output = new bool[max_nb_value];
  }
  else {
    used_output = NULL;
  }

  segmentation_likelihood = 0.;

  for (i = 0;i < nb_sequence;i++) {
    for (j = 1;j < nb_variable;j++) {
      if (model_type[j - 1] == VARIANCE_CHANGE) {
        mean[j] = 0.;
        if (type[j] != REAL_VALUE) {
          for (k = 0;k < length[i];k++) {
            mean[j] += int_sequence[i][j][k];
          }
        }
        else {
          for (k = 0;k < length[i];k++) {
            mean[j] += real_sequence[i][j][k];
          }
        }
        mean[j] /= length[i];
      }

      if (((model_type[j - 1] == LINEAR_MODEL_CHANGE) && (!seq_index_parameter)) ||
          ((j == 1) && (model_type[0] == INTERCEPT_SLOPE_CHANGE))) {
        if (index_parameter_type == IMPLICIT_TYPE) {
          seq_index_parameter = new int[max_length];
          for (k = 0;k < max_length;k++) {
            seq_index_parameter[k] = k;
          }
        }
        else {
          seq_index_parameter = index_parameter[i];
        }
      }

      if (model_type[j - 1] == BAYESIAN_POISSON_CHANGE) {
        gamma_hyperparameter_computation(i , j , hyperparam[j]);
      }
      if (model_type[j - 1] == BAYESIAN_GAUSSIAN_CHANGE) {
        gaussian_gamma_hyperparameter_computation(i , j , hyperparam[j]);
      }
    }

    // recurrence "forward"

    for (j = 0;j < length[i];j++) {

      // calcul des log-vraisemblances des segments

      if (model_type[0] == MULTIVARIATE_CATEGORICAL_CHANGE) {
        for (k = 0;k < max_nb_value;k++) {
          frequency[k] = 0;
        }

        for (k = j;k >= 0;k--) {
          for (m = 1;m < nb_variable;m++) {
            frequency[int_sequence[i][m][k]]++;
          }

          contrast[k] = 0.;
          for (m = 0;m < max_nb_value;m++) {
            if (frequency[m] > 0) {
              contrast[k] += frequency[m] * log((double)frequency[m] / (double)((nb_variable - 1) * (j - k + 1)));
            }
          }
        }
      }

      else  if (model_type[0] == MULTIVARIATE_POISSON_CHANGE) {
        for (k = 1;k < nb_variable;k++) {
          factorial[k][j] = 0.;
          for (m = 2;m <= int_sequence[i][k][j];m++) {
            factorial[k][j] += log((double)m);
          }
        }

        sum = 0.;
        factorial_sum = 0.;
        for (k = j;k >= 0;k--) {
          for (m = 1;m < nb_variable;m++) {
            sum += int_sequence[i][m][k];
            factorial_sum += factorial[m][k];
          }
          if (sum > 0.) {
            contrast[k] = sum * (log(sum / ((nb_variable - 1) * (j - k + 1))) - 1) - factorial_sum;
          }
          else {
            contrast[k] = 0.;
          }
        }
      }

      else if (model_type[0] == MULTIVARIATE_GEOMETRIC_0_CHANGE) {
        sum = 0.;
        for (k = j;k >= 0;k--) {
          for (m = 1;m < nb_variable;m++) {
            sum += int_sequence[i][m][k];
          }
          if (sum > 0.) {
            proba = (nb_variable - 1) * (j - k + 1) / ((nb_variable - 1) * (j - k + 1) + sum);
            contrast[k] = (nb_variable - 1) * (j - k + 1) * log(proba) + sum * log(1. - proba);
          }
          else {
            contrast[k] = 0.;
          }
        }
      }

      else {
        for (k = 0;k <= j;k++) {
          contrast[k] = 0.;
        }

        for (k = 1;k < nb_variable;k++) {
          if (model_type[k - 1] == CATEGORICAL_CHANGE) {
            for (m = 0;m < marginal_distribution[k]->nb_value;m++) {
              frequency[m] = 0;
            }
            sum = 0.;

            frequency[int_sequence[i][k][j]]++;
            for (m = j - 1;m >= 0;m--) {
              sum += (j - m) * log((double)(j - m) / (double)(j - m + 1)) +
                     log((double)(frequency[int_sequence[i][k][m]] + 1) / (double)(j - m + 1));
              if (frequency[int_sequence[i][k][m]] > 0) {
                sum -= frequency[int_sequence[i][k][m]] *
                       log((double)frequency[int_sequence[i][k][m]] / (double)(frequency[int_sequence[i][k][m]] + 1));
              }
              frequency[int_sequence[i][k][m]]++;

              if (contrast[m] != D_INF) {
                contrast[m] += sum;
              }

/*              frequency[int_sequence[i][k][m]]++;
              if (contrast[m] != D_INF) {
                for (n = 0;n < marginal_distribution[k]->nb_value;n++) {
                  if (frequency[n] > 0) {
                    contrast[m] += frequency[n] * log((double)frequency[n] / (double)(j - m + 1));
                  }
                }
              } */
            }
          }

          else if (model_type[k - 1] == POISSON_CHANGE) {
            factorial[k][j] = 0.;
            for (m = 2;m <= int_sequence[i][k][j];m++) {
              factorial[k][j] += log((double)m);
            }

            sum = 0.;
            factorial_sum = 0.;
            for (m = j;m >= 0;m--) {
              sum += int_sequence[i][k][m];
              factorial_sum += factorial[k][m];
              if ((contrast[m] != D_INF) && (sum > 0.)) {
                contrast[m] += sum * (log(sum / (j - m + 1)) - 1) - factorial_sum;
              }
            }
          }

          else if (model_type[k - 1] == GEOMETRIC_0_CHANGE) {
            sum = 0.;
            for (m = j;m >= 0;m--) {
              sum += int_sequence[i][k][m];
              if ((contrast[m] != D_INF) && (sum > 0.)) {
                proba = (j - m + 1) / (j - m + 1 + sum);
                contrast[m] += (j - m + 1) * log(proba) + sum * log(1. - proba);
              }
            }
          }

          else if (model_type[k - 1] == GEOMETRIC_1_CHANGE) {
            sum = 0.;
            for (m = j;m >= 0;m--) {
              sum += int_sequence[i][k][m];
              if ((contrast[m] != D_INF) && (sum > j - m + 1)) {
                proba = (j - m + 1) / sum;
                contrast[m] += (j - m + 1) * log(proba) + (sum - (j - m + 1)) * log(1. - proba);
              }
            }
          }

          else if (model_type[k - 1] == BAYESIAN_POISSON_CHANGE) {
            prior_contrast = -lgamma(hyperparam[k][0]) + hyperparam[k][0] * log(hyperparam[k][1]);

            factorial[k][j] = 0.;
            for (m = 2;m <= int_sequence[i][k][m];m++) {
              factorial[k][j] += log((double)m);
            }

            sum = 0.;
            factorial_sum = 0.;
            for (m = j;m >= 0;m--) {
              sum += int_sequence[i][k][m];
              factorial_sum += factorial[k][m];
              if (contrast[m] != D_INF) {
                contrast[m] += prior_contrast - factorial_sum + lgamma(hyperparam[k][0] + sum) -
                               (hyperparam[k][0] + sum) * log(hyperparam[k][1] + j - m + 1);
              }
            }
          }

          else if (model_type[k - 1] == BAYESIAN_GAUSSIAN_CHANGE) {
            prior_contrast = log(hyperparam[k][1]) / 2 - lgamma(hyperparam[k][2] / 2) +
                             hyperparam[k][2] * log(hyperparam[k][3] / 2) / 2;

            if (type[k] != REAL_VALUE) {
              square_sum = 0.;
              sum = int_sequence[i][k][j];
              if (contrast[j] != D_INF) {
                diff = hyperparam[k][0] - sum;
                contrast[j] += prior_contrast - log(2 * M_PI) / 2 -
                               log(hyperparam[k][1] + 1) / 2 + lgamma((hyperparam[k][2] + 1) / 2) -
                               (hyperparam[k][2] + 1) *
                               log((hyperparam[k][3] + hyperparam[k][1] *
                                    diff * diff / (hyperparam[k][1] + 1)) / 2) / 2;
              }

              for (m = j - 1;m >= 0;m--) {
                diff = int_sequence[i][k][m] - sum / (j - m);
                square_sum += ((double)(j - m) / (double)(j - m + 1)) * diff * diff;
                sum += int_sequence[i][k][m];
                if (contrast[m] != D_INF) {
                  diff = hyperparam[k][0] - sum / (j - m + 1);
                  contrast[m] += prior_contrast - (j - m + 1) * log(2 * M_PI) / 2 -
                                 log(hyperparam[k][1] + j - m + 1) / 2 +
                                 lgamma((hyperparam[k][2] + j - m + 1) / 2) -
                                 (hyperparam[k][2] + j - m + 1) *
                                 logl((hyperparam[k][3] + square_sum + hyperparam[k][1] * (j - m + 1) *
                                       diff * diff / (hyperparam[k][1] + j - m + 1)) / 2) / 2;
                }
              }
            }

            else {
              square_sum = 0.;
              sum = real_sequence[i][k][j];
              if (contrast[j] != D_INF) {
                diff = hyperparam[k][0] - sum;
                contrast[j] += prior_contrast - log(2 * M_PI) / 2 -
                               log(hyperparam[k][1] + 1) / 2 + lgamma((hyperparam[k][2] + 1) / 2) -
                               (hyperparam[k][2] + 1) *
                               log((hyperparam[k][3] + hyperparam[k][1] *
                                    diff * diff / (hyperparam[k][1] + 1)) / 2) / 2;
              }

              for (m = j - 1;m >= 0;m--) {
                diff = real_sequence[i][k][m] - sum / (j - m);
                square_sum += ((double)(j - m) / (double)(j - m + 1)) * diff * diff;
                sum += real_sequence[i][k][m];
                if (contrast[m] != D_INF) {
                  diff = hyperparam[k][0] - sum / (j - m + 1);
                  contrast[m] += prior_contrast - (j - m + 1) * log(2 * M_PI) / 2 -
                                 log(hyperparam[k][1] + j - m + 1) / 2 +
                                 lgamma((hyperparam[k][2] + j - m + 1) / 2) -
                                 (hyperparam[k][2] + j - m + 1) *
                                 logl((hyperparam[k][3] + square_sum + hyperparam[k][1] * (j - m + 1) *
                                       diff * diff / (hyperparam[k][1] + j - m + 1)) / 2) / 2;
                }
              }
            }
          }

          else {
            if (model_type[k - 1] == VARIANCE_CHANGE) {
              square_sum = 0.;

              if (type[k] != REAL_VALUE) {
                for (m = j;m >= 0;m--) {
                  diff = int_sequence[i][k][m] - mean[k];
                  square_sum += diff * diff;
                  residual[m] = square_sum;
                }
              }
              else {
                for (m = j;m >= 0;m--) {
                  diff = real_sequence[i][k][m] - mean[k];
                  square_sum += diff * diff;
                  residual[m] = square_sum;
                }
              }
            }

            else if (model_type[k - 1] == ORDINAL_GAUSSIAN_CHANGE) {
              square_sum = 0.;
              sum = rank[k][int_sequence[i][k][j]];
              residual[j] = 0.;

              for (m = j - 1;m >= 0;m--) {
                diff = rank[k][int_sequence[i][k][m]] - sum / (j - m);
                square_sum += ((double)(j - m) / (double)(j - m + 1)) * diff * diff;
                sum += rank[k][int_sequence[i][k][m]];
                residual[m] = square_sum;
              }
            }

            else if ((model_type[k - 1] == LINEAR_MODEL_CHANGE) || (model_type[0] == INTERCEPT_SLOPE_CHANGE)) {
              if (type[k] != REAL_VALUE) {
                square_sum = 0.;
                index_parameter_square_sum = 0.;
                mix_square_sum = 0.;
                sum = int_sequence[i][k][j];
                index_parameter_sum = seq_index_parameter[j];
                residual[j] = 0.;

                for (m = j - 1;m >= 0;m--) {
                  diff = int_sequence[i][k][m] - sum / (j - m);
                  square_sum += ((double)(j - m) / (double)(j - m + 1)) * diff * diff;
                  index_parameter_diff = seq_index_parameter[m] - index_parameter_sum / (j - m);
                  index_parameter_square_sum += ((double)(j - m) / (double)(j - m + 1)) *
                                                index_parameter_diff * index_parameter_diff;
                  mix_square_sum += ((double)(j - m) / (double)(j - m + 1)) * diff *  index_parameter_diff;
                  sum += int_sequence[i][k][m];
                  index_parameter_sum += seq_index_parameter[m];

                  if ((m < j - 1) && (index_parameter_square_sum > 0.)) {
                    residual[m] = square_sum - mix_square_sum * mix_square_sum / index_parameter_square_sum;
                  }
                  else {
                    residual[m] = 0.;
                  }
                }
              }

              else {
                square_sum = 0.;
                index_parameter_square_sum = 0.;
                mix_square_sum = 0.;
                sum = real_sequence[i][k][j];
                index_parameter_sum = seq_index_parameter[j];
                residual[j] = 0.;

                for (m = j - 1;m >= 0;m--) {
                  diff = real_sequence[i][k][m] - sum / (j - m);
                  square_sum += ((double)(j - m) / (double)(j - m + 1)) * diff * diff;
                  index_parameter_diff = seq_index_parameter[m] - index_parameter_sum / (j - m);
                  index_parameter_square_sum += ((double)(j - m) / (double)(j - m + 1)) *
                                                index_parameter_diff * index_parameter_diff;
                  mix_square_sum += ((double)(j - m) / (double)(j - m + 1)) * diff * index_parameter_diff;
                  sum += real_sequence[i][k][m];
                  index_parameter_sum += seq_index_parameter[m];

                  if ((m < j - 1) && (index_parameter_square_sum > 0.)) {
                    residual[m] = square_sum - mix_square_sum * mix_square_sum / index_parameter_square_sum;
                  }
                  else {
                    residual[m] = 0.;
                  }
                }
              }
            }

            else {
              if (type[k] != REAL_VALUE) {
                square_sum = 0.;
                sum = int_sequence[i][k][j];
                residual[j] = 0.;

                for (m = j - 1;m >= 0;m--) {
                  diff = int_sequence[i][k][m] - sum / (j - m);
                  square_sum += ((double)(j - m) / (double)(j - m + 1)) * diff * diff;
                  sum += int_sequence[i][k][m];
                  residual[m] = square_sum;
                }
              }

              else {
                square_sum = 0.;
                sum = real_sequence[i][k][j];
                residual[j] = 0.;

                for (m = j - 1;m >= 0;m--) {
                  diff = real_sequence[i][k][m] - sum / (j - m);
                  square_sum += ((double)(j - m) / (double)(j - m + 1)) * diff * diff;
                  sum += real_sequence[i][k][m];
                  residual[m] = square_sum;
                }
              }
            }

            if ((model_type[0] == MEAN_CHANGE) || (model_type[0] == INTERCEPT_SLOPE_CHANGE)) {
              for (m = j - 1;m >= 0;m--) {
                contrast[m] -= residual[m];
              }
            }

            else if (model_type[0] == MEAN_VARIANCE_CHANGE) {
              for (m = j - 1;m >= 0;m--) {
                contrast[m] += residual[m];
              }
            }

            else {
              for (m = j;m >= 0;m--) {
//                if ((contrast[m] != D_INF) && (residual[m] > 0.)) {
                if ((contrast[m] != D_INF) && (residual[m] > sqrt((double)(j - m + 1)) * ROUNDOFF_ERROR)) {
                  contrast[m] -= ((double)(j - m + 1) / 2.) * (logl(residual[m] /
                                   (j - m + 1)) + log(2 * M_PI) + 1);
/*                  contrast[m] -= ((double)(j - m + 1) / 2.) * (logl(residual[m] /
                                   (j - m)) + log(2 * M_PI)) + (double)(j - m) / 2.; */
                }
                else {
                  contrast[m] = D_INF;
                }
              }
            }
          }
        }

        if (model_type[0] == MEAN_VARIANCE_CHANGE) {
          contrast[j] = D_INF;
          for (k = j - 1;k >= 0;k--) {
//             if (contrast[k] > 0.) {
            if (contrast[k] > sqrt((double)((nb_variable - 1) * (j - k + 1))) * ROUNDOFF_ERROR) {
              contrast[k] = -((double)((nb_variable - 1) * (j - k + 1)) / 2.) * (logl(contrast[k] /
                               ((nb_variable - 1) * (j - k + 1))) + log(2 * M_PI) + 1);
/*              contrast[k] = -((double)((nb_variable - 1) * (j - k + 1)) / 2.) * (logl(contrast[k] /
                               ((nb_variable - 1) * (j - k))) + log(2 * M_PI)) +
                             (double)((nb_variable - 1) * (j - k)) / 2.; */
            }
            else {
              contrast[k] = D_INF;
            }
          }
        }
      }

      for (k = 0;k < MIN((j < length[i] - 1 ? nb_segment[i] - 1 : nb_segment[i]) , j + 1);k++) {
//      for (k = MAX(0 , nb_segment[i] + j - length[i]);k < MIN((j < length[i] - 1 ? nb_segment[i] - 1 : nb_segment[i]) , j + 1);k++) {
        if (k == 0) {
          forward[j][k] = contrast[0];
          if (forward[j][k] != D_INF) {
            optimal_length[j][k] = j + 1;
          }
        }

        else {
          forward[j][k] = D_INF;
          for (m = j;m >= k;m--) {
            if ((contrast[m] != D_INF) && (forward[m - 1][k - 1] != D_INF)) {
              buff = contrast[m] + forward[m - 1][k - 1];
              if (buff > forward[j][k]) {
                forward[j][k] = buff;
                optimal_length[j][k] = j - m + 1;
              }
            }
          }
        }
      }
    }

    if ((model_type[0] != MEAN_CHANGE) && (model_type[0] != INTERCEPT_SLOPE_CHANGE)) {
      if (isegmentation_likelihood) {
        for (j = 0;j < nb_segment[i];j++) {
          isegmentation_likelihood[j] = forward[length[i] - 1][j];
        }
      }

      if (forward[length[i] - 1][nb_segment[i] - 1] != D_INF) {
        segmentation_likelihood += forward[length[i] - 1][nb_segment[i] - 1];
      }
      else {
        segmentation_likelihood = D_INF;
        break;
      }
    }

    else {
      if (isegmentation_likelihood) {
        for (j = 0;j < nb_segment[i];j++) {
          if (forward[length[i] - 1][j] < 0.) {
            isegmentation_likelihood[j] = -((double)((nb_variable - 1) * length[i]) / 2.) *
                                           (log(-forward[length[i] - 1][j] /
                                             ((nb_variable - 1) * length[i])) + log(2 * M_PI) + 1);
/*            isegmentation_likelihood[j] = -(((double)((nb_variable - 1) * length[i]) / 2.) *
                                            (log(-forward[length[i] - 1][j] /
                                              ((nb_variable - 1) * (length[i] - nb_segment[i]))) + log(2 * M_PI)) +
                                            (double)((nb_variable - 1) * (length[i] - nb_segment[i])) / 2.); */
          }
          else {
            isegmentation_likelihood[j] = D_INF;
          }
        }
      }

      if (forward[length[i] - 1][nb_segment[i] - 1] < 0.) {
        segmentation_likelihood -= ((double)((nb_variable - 1) * length[i]) / 2.) *
                                   (log(-forward[length[i] - 1][nb_segment[i] - 1] /
                                     ((nb_variable - 1) * length[i])) + log(2 * M_PI) + 1);
/*        segmentation_likelihood -= (((double)((nb_variable - 1) * length[i]) / 2.) *
                                    (log(-forward[length[i] - 1][nb_segment[i] - 1] /
                                      ((nb_variable - 1) * (length[i] - nb_segment[i]))) + log(2 * M_PI)) +
                                    (double)((nb_variable - 1) * (length[i] - nb_segment[i])) / 2.); */
      }
      else {
        segmentation_likelihood = D_INF;
        break;
      }
    }

    // calcul du terme de penalite lie a la repartition des ruptures (BIC modifie)

    if (segment_penalty) {

#     ifdef DEBUG
      int cumul_segment_length;
      cout << "\n";
#     endif

      for (j = 0;j < nb_segment[i];j++) {
        segment_penalty[j] = 0.;
        k = length[i] - 1;

#       ifdef DEBUG
        cumul_segment_length = 0;
#       endif

        for (m = j;m >= 0;m--) {

#         ifdef DEBUG
          cout << optimal_length[k][m] << " ";
          cumul_segment_length += optimal_length[k][m];
#         endif

          segment_penalty[j] += log((double)optimal_length[k][m]);
          k -= optimal_length[k][m];
        }

#       ifdef DEBUG
        cout << "| " << segment_penalty[j] << endl;
        if (cumul_segment_length != length[i]) {
          cout << "\nERROR: " << j << "   " << cumul_segment_length << " | " << length[i] << endl;
        }
#       endif

      }
    }

    // calcul du nombre de parametres independants

    if (nb_parameter) {
      for (j = 0;j < nb_segment[i];j++) {
//        nb_parameter[j] = 0;
        nb_parameter[j] = j;

        if (model_type[0] == MULTIVARIATE_CATEGORICAL_CHANGE) {
          k = length[i] - 1;

          for (m = j;m >= 0;m--) {
            for (n = 0;n < max_nb_value;n++) {
              used_output[n] = false;
            }

            used_output[int_sequence[i][1][k]] = true;
            for (n = 2;n < nb_variable;n++) {
              if (!used_output[int_sequence[i][n][k]]) {
                nb_parameter[j]++;
                used_output[int_sequence[i][n][k]] = true;
              }
            }

            for (n = k - 1;n > k - optimal_length[k][m];n--) {
              for (r = 1;r < nb_variable;r++) {
                if (!used_output[int_sequence[i][r][n]]) {
                  nb_parameter[j]++;
                  used_output[int_sequence[i][r][n]] = true;
                }
              }
            }

            k -= optimal_length[k][m];
          }
        }

        else if ((model_type[0] == MULTIVARIATE_POISSON_CHANGE) || (model_type[0] == MULTIVARIATE_GEOMETRIC_0_CHANGE)) {
          nb_parameter[j] += j + 1;
        }
        else if (model_type[0] == MEAN_CHANGE) {
          nb_parameter[j] += (nb_variable - 1) * (j + 1) + 1;
        }
        else if (model_type[0] == MEAN_VARIANCE_CHANGE) {
          nb_parameter[j] += nb_variable * (j + 1);
        }
        else if (model_type[0] == INTERCEPT_SLOPE_CHANGE) {
          nb_parameter[j] += (nb_variable - 1) * (j + 1) * 2 + 1;
        }

        else {
          for (k = 1;k < nb_variable;k++) {
            if (model_type[k - 1] == CATEGORICAL_CHANGE) {
              m = length[i] - 1;

              for (n = j;n >= 0;n--) {
                for (r = 0;r < marginal_distribution[k]->nb_value;r++) {
                  used_output[r] = false;
                }
                used_output[int_sequence[i][k][m]] = true;

                for (r = m - 1;r > m - optimal_length[m][n];r--) {
                  if (!used_output[int_sequence[i][k][r]]) {
                    nb_parameter[j]++;
                    used_output[int_sequence[i][k][r]] = true;
                  }
                }
                m -= optimal_length[m][n];
              }
            }

            else if ((model_type[k - 1] == POISSON_CHANGE) || (model_type[k - 1] == GEOMETRIC_0_CHANGE) ||
                     (model_type[k - 1] == GEOMETRIC_1_CHANGE) || (model_type[k - 1] == BAYESIAN_POISSON_CHANGE)) {
              nb_parameter[j] += j + 1;
            }

            else if (model_type[k - 1] == VARIANCE_CHANGE) {
              nb_parameter[j] += j + 2;
            }

            else if (model_type[k - 1] == LINEAR_MODEL_CHANGE) {
              nb_parameter[j] += 3 * (j + 1);
            }

            else {
              nb_parameter[j] += 2 * (j + 1);
            }
          }
        }
      }
    }

    // restauration

    j = length[i] - 1;
    psegment = int_sequence[i][0] + j;

    for (k = nb_segment[i] - 1;k >= 0;k--) {
//      for (m = 0;m < optimal_length[j][k];m++) {
      for (m = j;m > j - optimal_length[j][k];m--) {
        *psegment-- = k;
      }
      j -= optimal_length[j][k];
    }
  }

  min_value[0] = 0;
  max_value[0] = max_nb_segment - 1;
  delete marginal_distribution[0];
  build_marginal_frequency_distribution(0);

  delete [] frequency;

  for (i = 1;i < nb_variable;i++) {
    delete [] factorial[i];
  }
  delete [] factorial;

  if (index_parameter_type == IMPLICIT_TYPE) {
    delete [] seq_index_parameter;
  }

  for (i = 1;i < nb_variable;i++) {
    delete [] hyperparam[i];
  }
  delete [] hyperparam;

  delete [] mean;
  delete [] residual;

  delete [] contrast;

  for (i = 0;i < max_length;i++) {
    delete [] forward[i];
  }
  delete [] forward;

  for (i = 0;i < max_length;i++) {
    delete [] optimal_length[i];
  }
  delete [] optimal_length;

  delete [] used_output;

  return segmentation_likelihood;
}


/*--------------------------------------------------------------*
 *
 *  Segmentation optimale de sequences.
 *
 *  arguments : reference sur un objet StatError, stream, nombres de segments,
 *              types des modeles, identificateur de la sequence,
 *              sortie (sequence ou residus).
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::segmentation(StatError &error , ostream &os , int *nb_segment ,
                                   int *model_type , int iidentifier , int output) const

{
  bool status = true;
  register int i , j;
  int index = I_DEFAULT , nb_parameter , *psegment;
  double segmentation_likelihood , segment_penalty , penalized_likelihood , **rank;
  const FrequencyDistribution **pmarginal;
  FrequencyDistribution *marginal;
  Sequences *seq , *iseq , *oseq;


  oseq = NULL;
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
    if ((model_type[i] == CATEGORICAL_CHANGE) || (model_type[i] == POISSON_CHANGE) ||
        (model_type[i] == GEOMETRIC_0_CHANGE) || (model_type[i] == GEOMETRIC_1_CHANGE) ||
        (model_type[i] == MULTIVARIATE_GEOMETRIC_0_CHANGE) || (model_type[i] == ORDINAL_GAUSSIAN_CHANGE) ||
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

  for (i = 0;i < nb_sequence;i++) {
    if ((index == I_DEFAULT) || (index == i)) {
      if ((nb_segment[i] < 1) || (nb_segment[i] > length[i] / 2)) {
        status = false;
        ostringstream error_message;
        error_message << SEQ_label[SEQL_SEQUENCE] << " " << identifier[i] << ": "
                      << SEQ_error[SEQR_NB_SEGMENT];
        error.update((error_message.str()).c_str());
      }
    }
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

    if (index != I_DEFAULT) {
      iseq = new Sequences(*seq , 1 , &index);
      delete seq;
      seq = iseq;
    }

    segmentation_likelihood = seq->segmentation(nb_segment , model_type , rank);

    for (i = 1;i < seq->nb_variable;i++) {
      delete [] rank[i];
    }
    delete [] rank;

    if (segmentation_likelihood != D_INF) {

#     ifdef MESSAGE
      if (seq->nb_sequence == 1) {
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

        nb_parameter = seq->nb_parameter_computation(0 , nb_segment[0] , model_type);

        penalized_likelihood = 2 * segmentation_likelihood - nb_parameter *
                               log((double)((seq->nb_variable - 1) * seq->length[0])) - segment_penalty;

        os << "\n" << nb_segment[0] << " " << (nb_segment[0] == 1 ? SEQ_label[SEQL_SEGMENT] : SEQ_label[SEQL_SEGMENTS])
           << "   2 * " << STAT_label[STATL_LIKELIHOOD] << ": " << 2 * segmentation_likelihood << "   "
           << nb_parameter << " " << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS]
           << "   2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " (Modified "  << STAT_criterion_word[BIC] << "): "
           << penalized_likelihood << endl;
      }
#     endif

      oseq = seq->segmentation_output(nb_segment , model_type , os , output);

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


/*--------------------------------------------------------------*
 *
 *  Calcul des l'entropies des segmentations et des ruptures.
 *
 *  arguments : indice de la sequence, nombre de segments, types des modeles,
 *              rangs (variables ordinales), pointeurs sur
 *              les vraisemblances de toutes les segmentations possibles,
 *              les entropies des segmentations, les entropies des ruptures
 *              en tenant compte ou non de rangs des ruptures,
 *              les entropies sous une hypothese de segmentations equiprobables et
 *              les entropies marginales.
 *
 *--------------------------------------------------------------*/

double Sequences::forward_backward(int index , int nb_segment , int *model_type ,
                                   double **rank , double *likelihood ,
                                   long double *segmentation_entropy ,
                                   long double *first_order_entropy ,
                                   long double *change_point_entropy , double *uniform_entropy ,
                                   long double *marginal_entropy) const

{
  register int i , j , k , m;
  int max_nb_value , *frequency , *seq_index_parameter = NULL;
  double sum , factorial_sum , proba , diff , index_parameter_sum , index_parameter_diff , buff ,
         rlikelihood , *mean , **hyperparam , **factorial , **nb_segmentation_forward ,
         **nb_segmentation_backward , **smoothed;
  long double square_sum , index_parameter_square_sum , mix_square_sum , segment_norm , sequence_norm ,
              lbuff , prior_contrast , *residual , *contrast , *normalized_contrast , *norm ,
              *backward_norm , *entropy_smoothed , *segment_predicted , **forward ,
              *forward_norm , **backward , **change_point , **forward_predicted_entropy ,
              **backward_predicted_entropy;


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

#     ifdef MESSAGE
      cout << "\nGamma hyperparameters: " << hyperparam[i][0] << " " << hyperparam[i][1] << endl;
#     endif

    }

    else if (model_type[i - 1] == BAYESIAN_GAUSSIAN_CHANGE) {
      hyperparam[i] = new double[4];
      gaussian_gamma_hyperparameter_computation(index , i , hyperparam[i]);

#     ifdef MESSAGE
      cout << "\nGaussian gamma hyperparameters: " << hyperparam[i][0] << " " << hyperparam[i][1]
           << " " << hyperparam[i][2] << " " << hyperparam[i][3] << endl;
#     endif

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

  mean = new double[nb_variable];
  residual = new long double[length[index]];

  contrast = new long double[length[index]];
  normalized_contrast = new long double[length[index]];

  nb_segmentation_forward = new double*[length[index]];
  for (i = 0;i < length[index];i++) {
    nb_segmentation_forward[i] = new double[nb_segment];
  }

  forward = new long double*[length[index]];
  for (i = 0;i < length[index];i++) {
    forward[i] = new long double[nb_segment];
  }

  segment_predicted = new long double[length[index]];

  forward_predicted_entropy = new long double*[length[index]];
  for (i = 0;i < length[index];i++) {
    forward_predicted_entropy[i] = new long double[nb_segment];
  }

  norm = new long double[length[index]];
  forward_norm = new long double[length[index]];

  nb_segmentation_backward = new double*[length[index]];
  for (i = 0;i < length[index];i++) {
    nb_segmentation_backward[i] = new double[nb_segment];
  }

  backward = new long double*[length[index]];
  for (i = 0;i < length[index];i++) {
    backward[i] = new long double[nb_segment];
  }

  backward_predicted_entropy = new long double*[length[index]];
  for (i = 0;i < length[index];i++) {
    backward_predicted_entropy[i] = new long double[nb_segment];
  }

  backward_norm = new long double[length[index]];

  smoothed = new double*[length[index]];
  for (i = 0;i < length[index];i++) {
    smoothed[i] = new double[nb_segment];
  }

  change_point = new long double*[nb_segment];
  for (i = 1;i < nb_segment;i++) {
    change_point[i] = new long double[length[index]];
  }

  entropy_smoothed = new long double[nb_segment];

  for (i = 1;i < nb_variable;i++) {
    if (model_type[i - 1] == VARIANCE_CHANGE) {
      mean[i] = 0.;
      if (type[i] != REAL_VALUE) {
        for (j = 0;j < length[index];j++) {
          mean[i] += int_sequence[index][i][j];
        }
      }
      else {
        for (j = 0;j < length[index];j++) {
          mean[i] += real_sequence[index][i][j];
        }
      }
      mean[i] /= length[index];
    }

    if ((model_type[i - 1] == LINEAR_MODEL_CHANGE) && (!seq_index_parameter)) {
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
                diff = int_sequence[index][j][k] - mean[j];
                square_sum += diff * diff;
                residual[k] = square_sum;
              }
            }

            else {
              for (k = i;k >= 0;k--) {
                diff = real_sequence[index][j][k] - mean[j];
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

          else if (model_type[j - 1] == LINEAR_MODEL_CHANGE) {
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

          if (model_type[0] == MEAN_VARIANCE_CHANGE) {
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

    // calcul du nombre de segmentations

    for (j = 0;j < nb_segment;j++) {
      nb_segmentation_forward[i][j] = 0;
    }

    for (j = 0;j < MIN((i < length[index] - 1 ? nb_segment - 1 : nb_segment) , i + 1);j++) {
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

    // recurrence et calcul des entropies predites

    if (contrast[i] != D_INF) {
      contrast[i] = expl(contrast[i]);
    }
    else {
      contrast[i] = 0.;
    }

    segment_norm = 0.;
    for (j = i - 1;j >= 0;j--) {
      segment_norm += norm[j];
      if (contrast[j] != D_INF) {
        contrast[j] = expl(contrast[j] - segment_norm);
      }
      else {
        contrast[j] = 0.;
      }
    }

    for (j = 0;j < nb_segment;j++) {
      forward[i][j] = 0.;
      forward_predicted_entropy[i][j] = 0.;
    }
    norm[i] = 0.;

    for (j = 0;j < MIN((i < length[index] - 1 ? nb_segment - 1 : nb_segment) , i + 1);j++) {
      if (j == 0) {
        forward[i][j] = contrast[0];
      }

      else {
        for (k = i;k >= j;k--) {
//          forward[i][j] += contrast[k] * forward[k - 1][j - 1];
          segment_predicted[k] = contrast[k] * forward[k - 1][j - 1];
          forward[i][j] += segment_predicted[k];
        }

        if (forward[i][j] > 0.) {
          for (k = i;k >= j;k--) {
            lbuff = segment_predicted[k] / forward[i][j];
            if (lbuff > 0.) {
              forward_predicted_entropy[i][j] += lbuff * (forward_predicted_entropy[k - 1][j - 1] - logl(lbuff));
            }
          }
        }
      }

      norm[i] += forward[i][j];
    }

    if (norm[i] > 0.) {
      for (j = 0;j < MIN((i < length[index] - 1 ? nb_segment - 1 : nb_segment) , i + 1);j++) {
        forward[i][j] /= norm[i];
      }

      norm[i] = logl(norm[i]);
    }

    forward_norm[i] = segment_norm + norm[i];
  }

  // calcul des entropies correspondant a une hypothèse de segmentations equiprobables
  // pour les differents nombres de segments possibles

  if (uniform_entropy) {
    for (i = 1;i < nb_segment;i++) {
      uniform_entropy[i] = log(nb_segmentation_forward[length[index] - 1][i]);
    }
  }

# ifdef DEBUG
  cout << "\n";
  buff = 1.;
  for (i = 1;i < nb_segment;i++) {
    buff *= (double)(length[index] - i) / (double)i;
    cout << i + 1 << " " << SEQ_label[SEQL_SEGMENTS] << ": "
         << nb_segmentation_forward[length[index] - 1][i] << " (" << buff << ") | "
         << log(nb_segmentation_forward[length[index] - 1][i]) << endl;
  }
# endif

  // extraction des log-vraisemblances de la sequence observee
  // pour les differents nombres de segments possibles

  for (i = 0;i < nb_segment;i++) {
    if (forward[length[index] - 1][i] > 0.) {
      likelihood[i] = logl(forward[length[index] - 1][i]) + forward_norm[length[index] - 1];
    }
    else {
      likelihood[i] = D_INF;
    }
  }

  rlikelihood = likelihood[nb_segment - 1];

  if (rlikelihood != D_INF) {
    for (i = 1;i < nb_segment;i++) {
      segmentation_entropy[i] = likelihood[i];
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
                  diff = int_sequence[index][j][k] - mean[j];
                  square_sum += diff * diff;
                  residual[k] = square_sum;
                }
              }

              else {
                for (k = i;k < length[index];k++) {
                  diff = real_sequence[index][j][k] - mean[j];
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

            else if (model_type[j - 1] == LINEAR_MODEL_CHANGE) {
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

            if (model_type[0] == MEAN_VARIANCE_CHANGE) {
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

      // calcul du nombre de segmentations

      for (j = 0;j < nb_segment;j++) {
        nb_segmentation_backward[i][j] = 0;
      }

      for (j = MAX((i == 0 ? 0 : 1) , nb_segment + i - length[index]);j < nb_segment;j++) {
        if (j < nb_segment - 1) {
          for (k = i;k <= length[index] + j - nb_segment;k++) {
            if (contrast[k] != D_INF) {
              nb_segmentation_backward[i][j] += nb_segmentation_backward[k + 1][j + 1];
            }
          }
        }

        else {
          if (contrast[length[index] - 1] != D_INF) {
            nb_segmentation_backward[i][j]++;
          }
        }
      }

      // recurrence et calcul des entropies predites

      if (contrast[i] != D_INF) {
        normalized_contrast[i] = expl(contrast[i]);
      }
      else {
        normalized_contrast[i] = 0.;
      }

      segment_norm = 0.;
      for (j = i + 1;j < length[index];j++) {
        segment_norm += norm[j];
        if (contrast[j] != D_INF) {
          normalized_contrast[j] = expl(contrast[j] - segment_norm);
        }
        else {
          normalized_contrast[j] = 0.;
        }
      }

      for (j = 0;j < nb_segment;j++) {
        backward[i][j] = 0.;
        backward_predicted_entropy[i][j] = 0.;
        smoothed[i][j] = 0.;
      }
      norm[i] = 0.;

      for (j = MAX((i == 0 ? 0 : 1) , nb_segment + i - length[index]);j < nb_segment;j++) {
        if (j < nb_segment - 1) {
          for (k = i;k <= length[index] + j - nb_segment;k++) {
//            backward[i][j] += normalized_contrast[k] * backward[k + 1][j + 1];
            segment_predicted[k] = normalized_contrast[k] * backward[k + 1][j + 1];
            backward[i][j] += segment_predicted[k];
          }

          if (backward[i][j] > 0.) {
            for (k = i;k <= length[index] + j - nb_segment;k++) {
              lbuff = segment_predicted[k] / backward[i][j];
              if (lbuff > 0.) {
                backward_predicted_entropy[i][j] += lbuff * (backward_predicted_entropy[k + 1][j + 1] - logl(lbuff));
              }
            }
          }
        }

        else {
          backward[i][j] = normalized_contrast[length[index] - 1];
        }

        norm[i] += backward[i][j];
      }

      if (norm[i] > 0.) {
        for (j = MAX((i == 0 ? 0 : 1) , nb_segment + i - length[index]);j < nb_segment;j++) {
          backward[i][j] /= norm[i];
        }

        norm[i] = logl(norm[i]);
      }

      backward_norm[i] = segment_norm + norm[i];

      // extraction des probabilitees lissees

      if (i < length[index] - 1) {
        for (j = MAX(0 , nb_segment + i - length[index]);j < MIN(nb_segment , i + 1);j++) {
          smoothed[i][j] = smoothed[i + 1][j];
          if (j > 0) {
            smoothed[i][j] -= forward[i][j - 1] * backward[i + 1][j] * sequence_norm;
          }
          if (j < nb_segment - 1) {
            smoothed[i][j] += forward[i][j] * backward[i + 1][j + 1] * sequence_norm;
          }

          if (smoothed[i][j] < 0.) {
            smoothed[i][j] = 0.;
          }
          if (smoothed[i][j] > 1.) {
            smoothed[i][j] = 1.;
          }
        }
      }

      else {
        smoothed[i][nb_segment - 1] = 1.;
      }

      if (i == 0) {
        sequence_norm = expl(backward_norm[i] - rlikelihood);
      }
      else {
        sequence_norm = expl(forward_norm[i - 1] + backward_norm[i] - rlikelihood);
      }

      // calcul des probabilites a posteriori des ruptures
      // pour les differents nombres de segments possibles

      if (i == 0) {

#       ifdef MESSAGE
        lbuff = backward[i][0] * sequence_norm;
        if ((lbuff < 1. - DOUBLE_ERROR) || (lbuff > 1. + DOUBLE_ERROR)) {
          cout << "\nERROR: " << lbuff << " | " << 1 << endl;
        }
#       endif

        for (j = 1;j < nb_segment;j++) {
          change_point[j][i] = 1.;
        }
      }

      else {
        for (j = 1;j < nb_segment;j++) {
          change_point[j][i] = 0.;
          for (k = MAX(1 , j + 1 + i - length[index]);k <= MIN(j , i);k++) {
            change_point[j][i] += forward[i - 1][k - 1] * backward[i][k + nb_segment - j - 1];
          }
          change_point[j][i] *= expl(forward_norm[i - 1] + backward_norm[i] - likelihood[j]);
        }
      }

      // calcul de l'entropie des segmentations

      segment_norm = 0.;
      for (j = i;j < length[index];j++) {
        segment_norm += norm[j];
        if (contrast[j] != D_INF) {
          normalized_contrast[j] = expl(contrast[j] - segment_norm);
        }
        else {
          normalized_contrast[j] = 0.;
        }
      }

      // calcul des entropies des segmentations pour les differents nombres de segments possibles

      if (i == 0) {
        for (j = i;j < length[index] - 1;j++) {
          if (contrast[j] != D_INF) {
            lbuff = normalized_contrast[j] * contrast[j];
            for (k = 1;k < MIN(nb_segment , length[index] - j);k++) {
              segmentation_entropy[k] -= backward[j + 1][nb_segment - k] *
                                         expl(backward_norm[i] - likelihood[k]) * lbuff;
//              segmentation_entropy[k] -= normalized_contrast[j] * backward[j + 1][nb_segment - k] *
//                                         expl(backward_norm[i] - likelihood[k]) * contrast[j];
            }
          }
        }
      }

      else {
        for (j = 1;j < nb_segment;j++) {
          if (j < nb_segment - 1) {
            for (k = i;k < length[index];k++) {
              if (contrast[k] != D_INF) {
                lbuff = forward[i - 1][j - 1] * normalized_contrast[k] * contrast[k];
                for (m = j + 1;m < MIN(nb_segment , j + length[index] - k);m++) {
                  segmentation_entropy[m] -= backward[k + 1][j + nb_segment - m] *
                                             expl(forward_norm[i - 1] + backward_norm[i] - likelihood[m]) * lbuff;
//                  segmentation_entropy[m] -= forward[i - 1][j - 1] * normalized_contrast[k] * backward[k + 1][j + nb_segment - m] *
//                                             expl(forward_norm[i - 1] + backward_norm[i] - likelihood[m]) * contrast[k];
                }
              }
            }
          }

          else {
            if (contrast[length[index] - 1] != D_INF) {
              lbuff = normalized_contrast[length[index] - 1] * contrast[length[index] - 1];
              for (k = 1;k < MIN(nb_segment , i + 1);k++) {
                segmentation_entropy[k] -= forward[i - 1][k - 1] *
                                           expl(forward_norm[i - 1] + backward_norm[i] - likelihood[k]) * lbuff;
//                segmentation_entropy[k] -= forward[i - 1][k - 1] * normalized_contrast[length[index] - 1] *
//                                           expl(forward_norm[i - 1] + backward_norm[i] - likelihood[k]) * contrast[length[index] - 1];
              }
            }
          }
        }
      }
    }

#   ifdef DEBUG
    cout << "\n" << SEQ_label[SEQL_SEGMENTATION_ENTROPY] << endl;
    for (i = 1;i < nb_segment;i++) {
      cout << i + 1 << " " << SEQ_label[SEQL_SEGMENTS] << ": "
           << forward_predicted_entropy[length[index] - 1][i] << ", "
           << backward_predicted_entropy[0][nb_segment - 1 - i] << ", "
           << segmentation_entropy[i] << endl;
    }
#   endif

#   ifdef MESSAGE
    for (i = 1;i < nb_segment;i++) {
      if (nb_segmentation_backward[0][nb_segment - 1 - i] != nb_segmentation_forward[length[index] - 1][i]) {
        cout << "\nERROR: " << i << "  " << nb_segmentation_forward[length[index] - 1][i]
             << " | " << nb_segmentation_backward[0][nb_segment - 1 - i] << endl;
      }
    }
#   endif

#   ifdef MESSAGE
    for (i = 0;i < length[index] - 1;i++) {
      sum = 0.;
      for (j = 0;j < nb_segment;j++) {
        sum += smoothed[i][j];
      }
      if ((sum < 1. - DOUBLE_ERROR) || (sum > 1. + DOUBLE_ERROR)) {
        cout << "\nERROR: " << i << " | " << sum << endl;
      }
    }

    for (i = 1;i < nb_segment;i++) {
      sum = 0.;
      for (j = 0;j < length[index];j++) {
        sum += change_point[i][j];
      }
      if ((sum < i + 1 - DOUBLE_ERROR) || (sum > i + 1 + DOUBLE_ERROR)) {
        cout << "\nERROR: " << sum << " | " << i + 1 << endl;
      }
    }
#   endif

    // calcul des entropies des ruptures ordonnees et des entropies marginales
    // pour les differents nombres de segments possibles

    for (i = 1;i < nb_segment;i++) {
      for (j = 0;j < i;j++) {
        entropy_smoothed[j] = 0.;
      }
      entropy_smoothed[i] = 1.;

      first_order_entropy[i] = 0.;
      marginal_entropy[i] = 0.;

      for (j = length[index] - 2;j >= 0;j--) {
        sequence_norm = expl(forward_norm[j] + backward_norm[j + 1] - likelihood[i]);

/*        for (k = MIN(i , j + 1) + 1;k <= i;k++) {
          entropy_smoothed[k] = 0.;
        } */

//        for (k = 0;k <= i;k++) {
        for (k = MAX(0 , i + 1 + j - length[index]);k <= MIN(i , j + 1);k++) {
          if (k > 0) {
//            entropy_smoothed[k] -= forward[j][k - 1] * backward[j + 1][k + nb_segment - i - 1] * sequence_norm;
            lbuff = forward[j][k - 1] * backward[j + 1][k + nb_segment - i - 1] * sequence_norm;
            entropy_smoothed[k] -= lbuff;
            if ((lbuff > 0.) && (lbuff < 1.)) {
              first_order_entropy[i] -= lbuff * logl(lbuff);
            }
          }
          if ((entropy_smoothed[k] > 0.) && (entropy_smoothed[k] < 1.)) {
            first_order_entropy[i] -= entropy_smoothed[k] * logl(entropy_smoothed[k]);
          }

          if (k < i) {
            entropy_smoothed[k] += forward[j][k] * backward[j + 1][k + nb_segment - i] * sequence_norm;
/*            lbuff = forward[j][k] * backward[j + 1][k + nb_segment - i] * sequence_norm;
            entropy_smoothed[k] += lbuff;
            if ((lbuff > 0.) && (lbuff < 1.)) {
              first_order_entropy[i] -= lbuff * logl(lbuff);
            } */
          }

          if (entropy_smoothed[k] < 0.) {
            entropy_smoothed[k] = 0.;
          }
          if (entropy_smoothed[k] > 1.) {
            entropy_smoothed[k] = 1.;
          }

          if (entropy_smoothed[k] > 0.) {
            first_order_entropy[i] += entropy_smoothed[k] * logl(entropy_smoothed[k]);
            marginal_entropy[i] -= entropy_smoothed[k] * logl(entropy_smoothed[k]);
          }
        }

#       ifdef MESSAGE
        sum = 0.;
        for (k = 0;k <= i;k++) {
          sum += entropy_smoothed[k];
        }
        if ((sum < 1. - DOUBLE_ERROR) || (sum > 1. + DOUBLE_ERROR)) {
          cout << "\nERROR: " << i + 1 << " " << j << " | " << sum << endl;
        }
#       endif

      }
    }

    // calcul des entropies des ruptures pour les differents nombres de segments possibles

    for (i = 1;i < nb_segment;i++) {
      change_point_entropy[i] = 0.;
      for (j = 1;j < length[index];j++) {
        if ((change_point[i][j] > 0.) && (change_point[i][j] < 1.)) {
          change_point_entropy[i] -= change_point[i][j] * logl(change_point[i][j]) +
                                     (1 - change_point[i][j]) * logl(1 - change_point[i][j]);
        }
      }
    }
  }

# ifdef DEBUG
  double **segment_length;
  Distribution *prior;


  // calcul des lois a priori des longueurs de segments sous une hypothese de loi a priori
  // uniforme sur les segmentations possibles

  segment_length = new double*[nb_segment];
  for (i = 0;i < nb_segment;i++) {
    segment_length[i] = new double[length[index] - nb_segment + 2];
    for (j = 0;j <= length[index] - nb_segment + 1;j++) {
      segment_length[i][j] = 0.;
    }
  }

  // recurrence "forward"

  for (i = 0;i < length[index];i++) {
    for (j = 0;j < nb_segment;j++) {
      nb_segmentation_forward[i][j] = 0;
    }

    for (j = MAX(0 , nb_segment + i - length[index]);j < MIN((i < length[index] - 1 ? nb_segment - 1 : nb_segment) , i + 1);j++) {
      if (j == 0) {
        nb_segmentation_forward[i][j]++;
      }

      else {
        for (k = i;k >= j;k--) {
          nb_segmentation_forward[i][j] += nb_segmentation_forward[k - 1][j - 1];
        }
      }
    }
  }

  // recurrence "backward"

  for (i = length[index] - 1;i > 0;i--) {
    for (j = 0;j < nb_segment;j++) {
      nb_segmentation_backward[i][j] = 0;
    }

    for (j = MAX(1 , nb_segment + i - length[index]);j < MIN(nb_segment , i + 1);j++) {
      if (j < nb_segment - 1) {
        for (k = i;k <= length[index] + j - nb_segment;k++) {
          nb_segmentation_backward[i][j] += nb_segmentation_backward[k + 1][j + 1];
          segment_length[j][k - i + 1] += nb_segmentation_forward[i - 1][j - 1] *
                                          nb_segmentation_backward[k + 1][j + 1];
        }
      }

      else {
        nb_segmentation_backward[i][j]++;
        segment_length[j][length[index] - i] += nb_segmentation_forward[i - 1][j - 1];
      }
    }
  }

  for (i = 0;i <= length[index] - nb_segment;i++) {
    segment_length[0][i + 1] += nb_segmentation_backward[i + 1][1];
  }

  for (i = 0;i < nb_segment;i++) {
    cout << "\n" << SEQ_label[SEQL_SEGMENT] << " " << i << ":";

    sum = 0.;
    for (j = 1;j <= length[index] - nb_segment + 1;j++) {
      sum += segment_length[i][j];
    }

    for (j = 1;j <= length[index] - nb_segment + 1;j++) {
      cout << " " << segment_length[i][j] / sum;
    }
    cout << endl;
  }

  prior = new Distribution(length[index] - nb_segment + 2);
  prior->mass[0] = 0.;

  buff = 1.;
  for (i = 1;i < nb_segment - 1;i++) {
    buff *= (double)(length[index] - i - 1) / (double)i;
  }
  sum = buff * (double)(length[index] - 1) / (double)(nb_segment - 1);

//  cout << "\nAll " << SEQ_label[SEQL_SEGMENTS] << ":";
  for (i = 1;i <= length[index] - nb_segment + 1;i++) {
//    cout << " " << buff / sum;
    prior->mass[i] = buff / sum;
    buff *= (double)(length[index] - i - nb_segment + 1) /
            (double)(length[index] - i - 1);
  }
//  cout << endl;

  prior->offset = 1.;
  prior->cumul_computation();
  prior->max = prior->mass[1];
  prior->mean_computation();
  prior->variance_computation();

  cout << "\n";
  prior->ascii_characteristic_print(cout , true);
//  prior->spreadsheet_characteristic_print(cout , true);
  prior->ascii_print(cout , false , true , false);

  for (i = 0;i < nb_segment;i++) {
    delete [] segment_length[i];
  }
  delete [] segment_length;

  delete prior;
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

  delete [] mean;
  delete [] residual;

  if (index_parameter_type == IMPLICIT_TYPE) {
    delete [] seq_index_parameter;
  }

  delete [] contrast;
  delete [] normalized_contrast;

  for (i = 0;i < length[index];i++) {
    delete [] nb_segmentation_forward[i];
  }
  delete [] nb_segmentation_forward;

  for (i = 0;i < length[index];i++) {
    delete [] forward[i];
  }
  delete [] forward;

  delete [] segment_predicted;

  for (i = 0;i < length[index];i++) {
    delete [] forward_predicted_entropy[i];
  }
  delete [] forward_predicted_entropy;

  delete [] norm;
  delete [] forward_norm;

  for (i = 0;i < length[index];i++) {
    delete [] nb_segmentation_backward[i];
  }
  delete [] nb_segmentation_backward;

  for (i = 0;i < length[index];i++) {
    delete [] backward[i];
  }
  delete [] backward;

  for (i = 0;i < length[index];i++) {
    delete [] backward_predicted_entropy[i];
  }
  delete [] backward_predicted_entropy;

  delete [] backward_norm;

  for (i = 0;i < length[index];i++) {
    delete [] smoothed[i];
  }
  delete [] smoothed;

  for (i = 1;i < nb_segment;i++) {
    delete [] change_point[i];
  }
  delete [] change_point;

  delete [] entropy_smoothed;

  return rlikelihood;
}


/*--------------------------------------------------------------*
 *
 *  Heuristique de la pente : estimation de la pente des log-vraisemblances
 *  fonction du nombre de segments.
 *
 *  arguments : nombres de segments minimum et maximum, type de forme de la penalite,
 *              longueur de la sequence, pointeur sur les formes de la penalite et
 *              sur les vraisemblances, references sur l'intercept,
 *              la pente, l'ecart type associe et l'ecart type residuel.
 *
 *--------------------------------------------------------------*/

void log_likelihood_slope(int min_nb_segment , int max_nb_segment , int penalty_shape_type ,
                          int length , double *penalty_shape , double *likelihood ,
                          double &intercept , double &slope , double &slope_standard_deviation ,
                          double &residual_standard_deviation)

{
  register int i;
  int nb_model;
  double diff , likelihood_mean , penalty_shape_mean , penalty_shape_variance , covariance ,
         residual_mean , residual_square_sum;


  nb_model = max_nb_segment - min_nb_segment + 1;

  likelihood_mean = 0.;
  for (i = min_nb_segment;i <= max_nb_segment;i++) {
    likelihood_mean += likelihood[i];
  }
  likelihood_mean /= nb_model;

  switch (penalty_shape_type) {

  case 0 : {
    penalty_shape_mean = (double)(min_nb_segment + max_nb_segment) / 2.;
    penalty_shape_variance = (double)(nb_model * (nb_model * nb_model - 1)) / 12.;

#   ifdef DEBUG
    double buff = 0.;
    for (i = min_nb_segment;i <= max_nb_segment;i++) {
      diff = i - penalty_shape_mean;
      buff += diff * diff;
    }

    cout << "TEST: " << penalty_shape_variance << " | " << buff << endl;
#   endif

    break;
  }

  default : {
    penalty_shape_mean = 0.;
    for (i = min_nb_segment;i <= max_nb_segment;i++) {
      penalty_shape_mean += penalty_shape[i];
    }
    penalty_shape_mean /= nb_model;

    penalty_shape_variance = 0.;
    for (i = min_nb_segment;i <= max_nb_segment;i++) {
      diff = penalty_shape[i] - penalty_shape_mean;
      penalty_shape_variance += diff * diff;
    }
    break;
  }
  }

  covariance = 0.;
  for (i = min_nb_segment;i <= max_nb_segment;i++) {
    covariance += (likelihood[i] - likelihood_mean) * (penalty_shape[i] - penalty_shape_mean);
  }

  slope = covariance / penalty_shape_variance;
  intercept = likelihood_mean - slope * penalty_shape_mean;

  residual_mean = 0.;
  residual_square_sum = 0.;
  for (i = min_nb_segment;i <= max_nb_segment;i++) {
    diff = likelihood[i] - (intercept + slope * penalty_shape[i]);
    residual_mean += diff;
    residual_square_sum += diff * diff;
  }
  residual_mean /= nb_model;

  if (nb_model > 2) {
    residual_square_sum /= (nb_model - 2);
    slope_standard_deviation = sqrt(residual_square_sum / penalty_shape_variance);
  }

  residual_standard_deviation = 0.;
  for (i = min_nb_segment;i <= max_nb_segment;i++) {
    diff = likelihood[i] - (intercept + slope * penalty_shape[i]) - residual_mean;
    residual_standard_deviation += diff * diff;
  }
  if (nb_model > 2) {
    residual_standard_deviation = sqrt(residual_standard_deviation / (nb_model - 2));
  }
}


/*--------------------------------------------------------------*
 *
 *  Segmentation optimale d'une sequence.
 *
 *  arguments : reference sur un objet StatError, stream, identificateur de la sequence,
 *              nombre de segments maximum, types des modeles, critere de selection
 *              du nombre de segments, nombre de segments minimum et
 *              forme de la penalite pour l'heuristique de la pente,
 *              sortie (sequence, entropies ou divergences).
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::segmentation(StatError &error , ostream &os , int iidentifier ,
                                   int max_nb_segment , int *model_type , int criterion ,
                                   int min_nb_segment , int penalty_shape_type , int output) const

{
  bool status = true , bayesian;
  register int i , j;
  int index , nb_segment , inb_sequence , *nb_parameter , inb_segment[1] , ilength[4] , itype[1];
  double buff , segmentation_intercept , segmentation_slope , segmentation_slope_standard_deviation ,
         segmentation_residual_standard_deviation , intercept , slope , slope_standard_deviation ,
         residual_standard_deviation , scaling_factor , max_likelihood[4] , *segmentation_likelihood ,
         *segment_penalty , *penalty_shape , **penalized_likelihood , *likelihood , *uniform_entropy ,
         *segmentation_divergence , **rank;
  long double *segmentation_entropy , *first_order_entropy , *change_point_entropy , *marginal_entropy;
  const FrequencyDistribution **pmarginal;
  FrequencyDistribution *marginal;
  Sequences *iseq , *seq , *oseq;


  oseq = NULL;
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

  if ((model_type[0] == MEAN_CHANGE) || (model_type[0] == INTERCEPT_SLOPE_CHANGE)) {
    if ((output != SEQUENCE) && (output != LOG_LIKELIHOOD_SLOPE)) {
      status = false;
      ostringstream correction_message;
      correction_message << SEQ_label[SEQL_SEQUENCE] << " or "
                         << STAT_criterion_word[LOG_LIKELIHOOD_SLOPE];
      error.correction_update(SEQ_error[SEQR_FORBIDDEN_OUTPUT] , (correction_message.str()).c_str());
    }

    if (criterion == LIKELIHOOD_SLOPE) {
      criterion = SEGMENTATION_LIKELIHOOD_SLOPE;
    }
    else if (criterion == ICL) {
      criterion = mBIC;
    }

/*    if ((criterion == LIKELIHOOD_SLOPE) || (criterion == ICL)) {
      status = false;
      ostringstream correction_message;
      correction_message << STAT_criterion_word[SEGMENTATION_LIKELIHOOD_SLOPE] << " or "
                         << STAT_criterion_word[mBIC];
      error.correction_update(SEQ_error[SEQR_FORBIDDEN_CRITERION] , (correction_message.str()).c_str());
    } */
  }

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

  else {
    if ((max_nb_segment < 2) || (max_nb_segment > length[index] / 2)) {
      status = false;
      error.update(SEQ_error[SEQR_MAX_NB_SEGMENT]);
    }

    if (min_nb_segment == 0) {
      min_nb_segment = max_nb_segment / 2;

      if (criterion == LIKELIHOOD_SLOPE) {
        criterion = ICL;
      }
      else if (criterion == SEGMENTATION_LIKELIHOOD_SLOPE) {
        criterion = mBIC;
      }
    }

    else if (min_nb_segment < 2) {
      status = false;
      error.update(SEQ_error[SEQR_MIN_NB_SEGMENT]);
    }
  }

  if (status) {
    if (max_nb_segment - min_nb_segment < SLOPE_NB_SEGMENT_RANGE) {
      if (criterion == LIKELIHOOD_SLOPE) {
        criterion = ICL;
      }
      else if (criterion == SEGMENTATION_LIKELIHOOD_SLOPE) {
        criterion = mBIC;
      }

      if (output == LOG_LIKELIHOOD_SLOPE) {
        output = SEQUENCE;
      }
    }

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

    segmentation_likelihood = new double[max_nb_segment + 2];
    nb_parameter = new int[max_nb_segment + 2];

    if ((model_type[0] == BAYESIAN_POISSON_CHANGE) || (model_type[0] == BAYESIAN_GAUSSIAN_CHANGE)) {
      bayesian = true;
//      nb_segment = 2;

      penalized_likelihood = new double*[1];
      penalized_likelihood[0] = new double[max_nb_segment + 1];
    }

    else {
      bayesian = false;

      segment_penalty = new double[max_nb_segment + 2];

      penalized_likelihood = new double*[4];
      for (i = 0;i < 4;i++) {
        penalized_likelihood[i] = new double[max_nb_segment + 1];
      }
    }

    if ((model_type[0] != MEAN_CHANGE) && (model_type[0] != INTERCEPT_SLOPE_CHANGE)) {
      likelihood = new double[max_nb_segment + 1];
      segmentation_entropy = new long double[max_nb_segment + 1];
      first_order_entropy = new long double[max_nb_segment + 1];
      change_point_entropy = new long double[max_nb_segment + 1];
      uniform_entropy = new double[max_nb_segment + 1];
      segmentation_divergence = new double[max_nb_segment + 1];
      marginal_entropy = new long double[max_nb_segment + 1];
    }

    inb_segment[0] = max_nb_segment + 1;

    if (max_nb_segment - min_nb_segment >= SLOPE_NB_SEGMENT_RANGE) {
      penalty_shape = new double[max_nb_segment + 1];

      switch (penalty_shape_type) {

      case 0 : {
        for (i = 1;i <= max_nb_segment;i++) {
          penalty_shape[i] = i;
        }
        break;
      }

      case 1 : {
        for (i = 1;i <= max_nb_segment;i++) {
          penalty_shape[i] = i - (double)(i * (i - 1)) / (double)(2 * seq->length[0]);
        }
        break;
      }

      case 2 : {
        buff = 1.;
        for (i = 1;i <= max_nb_segment;i++) {
          penalty_shape[i] = log(buff);
          buff *= (double)seq->length[0] / (double)i;
        }
        break;
      }

      case 3 : {
        for (i = 1;i <= max_nb_segment;i++) {
          penalty_shape[i] = i - (double)(i * (i - 1)) / (double)seq->length[0];
        }
        break;
      }

      case 4 : {
        buff = 1.;
        for (i = 1;i <= max_nb_segment;i++) {
          penalty_shape[i] = log(buff);
          buff *= (double)(seq->length[0] - i) / (double)i;
        }
        break;
      }
      }

#     ifdef DEBUG
      if (penalty_shape_type != 0) {
        os << "\nPenalty shape: ";
        for (i = 1;i <= max_nb_segment;i++) {
          os << " " << penalty_shape[i];
        }
        os << endl;
      }
#     endif

    }

    seq->segmentation(inb_segment , model_type , rank , segmentation_likelihood + 1 ,
                      nb_parameter + 1 , (bayesian ? NULL : segment_penalty + 1));

    if (max_nb_segment - min_nb_segment >= SLOPE_NB_SEGMENT_RANGE) {
      log_likelihood_slope(min_nb_segment , max_nb_segment , penalty_shape_type , seq->length[0] ,
                           penalty_shape , segmentation_likelihood , segmentation_intercept , segmentation_slope ,
                           segmentation_slope_standard_deviation , segmentation_residual_standard_deviation);
    }

    if ((model_type[0] != MEAN_CHANGE) && (model_type[0] != INTERCEPT_SLOPE_CHANGE)) {
      seq->forward_backward(0 , max_nb_segment , model_type , rank ,
                            likelihood + 1 , segmentation_entropy + 1 ,
                            first_order_entropy + 1 , change_point_entropy + 1 ,
                            uniform_entropy + 1 , marginal_entropy + 1);

      segmentation_divergence[1] = 0;
      for (i = 2;i <= max_nb_segment;i++) {
        segmentation_divergence[i] = uniform_entropy[i] - segmentation_entropy[i];
      }

      if (max_nb_segment - min_nb_segment >= SLOPE_NB_SEGMENT_RANGE) {
        log_likelihood_slope(min_nb_segment , max_nb_segment , penalty_shape_type , seq->length[0] ,
                             penalty_shape , likelihood , intercept , slope ,
                             slope_standard_deviation , residual_standard_deviation);

#       ifdef MESSAGE
        if (penalty_shape_type != 0) {
          os << "\nTEST, " << STAT_criterion_word[LIKELIHOOD_SLOPE] << ": " << slope << " | "
             << (likelihood[max_nb_segment] - likelihood[max_nb_segment - 1]) /
                (penalty_shape[max_nb_segment] - penalty_shape[max_nb_segment - 1]) << endl;
        }
#       endif

      }
    }

    switch (bayesian) {

    case false : {

      // calcul des vraisemblances penalisees au sens du critere ICL, de l'heuristique de pente et
      // du BIC modifie (Zhang & Siegmund, 2007)

//      segmentation_likelihood[1] = seq->one_segment_likelihood(0 , model_type , rank);
//      nb_parameter[1] = seq->nb_parameter_computation(0 , 1 , model_type);

      if ((model_type[0] != MEAN_CHANGE) && (model_type[0] != INTERCEPT_SLOPE_CHANGE) &&
          (likelihood[1] != D_INF)) {
        penalized_likelihood[0][1] = 2 * likelihood[1] - nb_parameter[1] *
                                     log((double)((seq->nb_variable - 1) * seq->length[0]));
        max_likelihood[0] = penalized_likelihood[0][1];

        if (max_nb_segment - min_nb_segment >= SLOPE_NB_SEGMENT_RANGE) {
          penalized_likelihood[1][1] = 2 * (likelihood[1] - 2 * penalty_shape[1] * slope);
          max_likelihood[1] = penalized_likelihood[1][1];
        }

        nb_segment = 1;
      }

      if (segmentation_likelihood[1] != D_INF) {
        penalized_likelihood[2][1] = 2 * segmentation_likelihood[1] - nb_parameter[1] *
                                     log((double)((seq->nb_variable - 1) * seq->length[0])) - segment_penalty[1];
        max_likelihood[2] = penalized_likelihood[2][1];

        if (max_nb_segment - min_nb_segment >= SLOPE_NB_SEGMENT_RANGE) {
          penalized_likelihood[3][1] = 2 * (segmentation_likelihood[1] - 2 * penalty_shape[1] * segmentation_slope);
          max_likelihood[3] = penalized_likelihood[3][1];
        }

        if ((model_type[0] == MEAN_CHANGE) || (model_type[0] == INTERCEPT_SLOPE_CHANGE)) {
          nb_segment = 1;
        }
      }

      if (((model_type[0] != MEAN_CHANGE) && (model_type[0] != INTERCEPT_SLOPE_CHANGE) &&
           (likelihood[1] == D_INF)) || (segmentation_likelihood[1] == D_INF)) {
        max_nb_segment = 0;
        nb_segment = 0;
      }

//      inb_segment[0] = 2;
//      segmentation_likelihood[2] = seq->segmentation(inb_segment , model_type , rank);
//      nb_parameter[2] = seq->nb_parameter_computation(0 , 2 , model_type);

      for (i = 2;i <= max_nb_segment;i++) {
//        inb_segment[0] = i + 1;
//        segmentation_likelihood[i + 1] = seq->segmentation(inb_segment , model_type , rank);
//        nb_parameter[i + 1] = seq->nb_parameter_computation(0 , i + 1 , model_type);

        if ((model_type[0] != MEAN_CHANGE) && (model_type[0] != INTERCEPT_SLOPE_CHANGE) &&
            (likelihood[i] != D_INF)) {
          penalized_likelihood[0][i] = 2 * (likelihood[i] - segmentation_entropy[i]) - nb_parameter[i] *
                                       log((double)((seq->nb_variable - 1) * seq->length[0]));
          if (penalized_likelihood[0][i] > max_likelihood[0]) {
            max_likelihood[0] = penalized_likelihood[0][i];
            if (criterion == ICL) {
              nb_segment = i;
            }
          }

          if (max_nb_segment - min_nb_segment >= SLOPE_NB_SEGMENT_RANGE) {
            penalized_likelihood[1][i] = 2 * (likelihood[i] - 2 * penalty_shape[i] * slope);
            if (penalized_likelihood[1][i] > max_likelihood[1]) {
              max_likelihood[1] = penalized_likelihood[1][i];
              if (criterion == LIKELIHOOD_SLOPE) {
                nb_segment = i;
              }
            }
          }
        }

        if (segmentation_likelihood[i] != D_INF) {
          penalized_likelihood[2][i] = 2 * segmentation_likelihood[i] - nb_parameter[i] *
                                       log((double)((seq->nb_variable - 1) * seq->length[0])) - segment_penalty[i];
          if (penalized_likelihood[2][i] > max_likelihood[2]) {
            max_likelihood[2] = penalized_likelihood[2][i];
            if (criterion == mBIC) {
              nb_segment = i;
            }
          }

          if (max_nb_segment - min_nb_segment >= SLOPE_NB_SEGMENT_RANGE) {
            penalized_likelihood[3][i] = 2 * (segmentation_likelihood[i] -
                                          2 * penalty_shape[i] * segmentation_slope);
            if (penalized_likelihood[3][i] > max_likelihood[3]) {
              max_likelihood[3] = penalized_likelihood[3][i];
              if (criterion == SEGMENTATION_LIKELIHOOD_SLOPE) {
                nb_segment = i;
              }
            }
          }
        }

        if (((model_type[0] != MEAN_CHANGE) && (model_type[0] != INTERCEPT_SLOPE_CHANGE) &&
             (likelihood[i] == D_INF)) || (segmentation_likelihood[i] == D_INF)) {
          max_nb_segment = i - 1;
          break;
        }
      }
      break;
    }

    case true : {
      if (likelihood[1] != D_INF) {
        penalized_likelihood[0][1] = 2 * likelihood[1];
        max_likelihood[0] = penalized_likelihood[0][1];
        nb_segment = 1;
      }

      else {
        max_nb_segment = 0;
        nb_segment = 0;
      }

      for (i = 2;i <= max_nb_segment;i++) {
        if (likelihood[i] != D_INF) {
//          penalized_likelihood[0][i] = 2 * (likelihood[i] - segmentation_entropy[i]);
          penalized_likelihood[0][i] = 2 * (likelihood[i] - segmentation_entropy[i] -
                                            uniform_entropy[i]);

          if (penalized_likelihood[0][i] > max_likelihood[0]) {
            max_likelihood[0] = penalized_likelihood[0][i];
            nb_segment = i;
          }
        }

        else {
          max_nb_segment = i - 1;
          break;
        }
      }
      break;
    }
    }

    if (nb_segment > 0) {

#     ifdef MESSAGE
      int width[20];
      long old_adjust;
      double norm , *posterior_probability , **weight;
      Test *test;


      if ((model_type[0] != MEAN_CHANGE) && (model_type[0] != INTERCEPT_SLOPE_CHANGE)) {
        posterior_probability = new double[max_nb_segment + 1];

        likelihood[1] = segmentation_likelihood[1];
        posterior_probability[1] = 1.;
        for (i = 2;i <= max_nb_segment;i++) {
          posterior_probability[i] = exp(segmentation_likelihood[i] - likelihood[i]);
        }
      }

      switch (bayesian) {

      case false : {
        weight = new double*[4];
        for (i = 0;i < 4;i++) {
          weight[i] = new double[max_nb_segment + 1];
        }
        break;
      }

      case true : {
        weight = new double*[1];
        weight[0] = new double[max_nb_segment + 1];
        break;
      }
      }

      if ((model_type[0] != MEAN_CHANGE) && (model_type[0] != INTERCEPT_SLOPE_CHANGE)) {
        norm = 0.;
        for (i = 1;i <= max_nb_segment;i++) {
          weight[0][i] = exp((penalized_likelihood[0][i] - max_likelihood[0]) / 2);
          norm += weight[0][i];
        }
        for (i = 1;i <= max_nb_segment;i++) {
          weight[0][i] /= norm;
        }
      }

      if (!bayesian) {
        if ((model_type[0] != MEAN_CHANGE) && (model_type[0] != INTERCEPT_SLOPE_CHANGE) &&
            (max_nb_segment - min_nb_segment >= SLOPE_NB_SEGMENT_RANGE)) {
          norm = 0.;
          for (i = 1;i <= max_nb_segment;i++) {
            weight[1][i] = exp((penalized_likelihood[1][i] - max_likelihood[1]) / 2);
            norm += weight[1][i];
          }
          for (i = 1;i <= max_nb_segment;i++) {
            weight[1][i] /= norm;
          }

          test = new Test(STUDENT , false , max_nb_segment - min_nb_segment - 1 , I_DEFAULT , D_DEFAULT);
          test->critical_probability = ref_critical_probability[0];
          test->t_value_computation();

          os << STAT_criterion_word[LIKELIHOOD_SLOPE] << ": " << slope << " ("
             << slope - test->value * slope_standard_deviation << ", "
             << slope + test->value * slope_standard_deviation << ") | "
             << STAT_label[STATL_RESIDUAL] << " " << STAT_label[STATL_STANDARD_DEVIATION] << ": "
             << residual_standard_deviation << endl;

          delete test;
        }

        norm = 0.;
        for (i = 1;i <= max_nb_segment;i++) {
          weight[2][i] = exp((penalized_likelihood[2][i] - max_likelihood[2]) / 2);
          norm += weight[2][i];
        }
        for (i = 1;i <= max_nb_segment;i++) {
          weight[2][i] /= norm;
        }

        if (max_nb_segment - min_nb_segment >= SLOPE_NB_SEGMENT_RANGE) {
          norm = 0.;
          for (i = 1;i <= max_nb_segment;i++) {
            weight[3][i] = exp((penalized_likelihood[3][i] - max_likelihood[3]) / 2);
            norm += weight[3][i];
          }
          for (i = 1;i <= max_nb_segment;i++) {
            weight[3][i] /= norm;
          }

          test = new Test(STUDENT , false , max_nb_segment - min_nb_segment - 1 , I_DEFAULT , D_DEFAULT);
          test->critical_probability = ref_critical_probability[0];
          test->t_value_computation();

          os << STAT_criterion_word[SEGMENTATION_LIKELIHOOD_SLOPE] << ": " << segmentation_slope << " ("
             << segmentation_slope - test->value * segmentation_slope_standard_deviation << ", "
             << segmentation_slope + test->value * segmentation_slope_standard_deviation << ") | "
             << STAT_label[STATL_RESIDUAL] << " " << STAT_label[STATL_STANDARD_DEVIATION] << ": "
             << segmentation_residual_standard_deviation << endl;

          delete test;
        }
      }

      old_adjust = os.flags(ios::adjustfield);

      if ((model_type[0] != MEAN_CHANGE) && (model_type[0] != INTERCEPT_SLOPE_CHANGE)) {
        width[0] = stat_tool::column_width(max_nb_segment) + ASCII_SPACE;
        width[1] = stat_tool::column_width(max_nb_segment , segmentation_likelihood + 1 , 2.) + ASCII_SPACE;
        width[2] = stat_tool::column_width(max_nb_segment , likelihood + 1 , 2.) + ASCII_SPACE;
        width[3] = stat_tool::column_width(max_nb_segment , posterior_probability + 1) + ASCII_SPACE;
        width[4] = column_width(max_nb_segment - 1 , segmentation_entropy + 2) + ASCII_SPACE;
        width[5] = column_width(max_nb_segment - 1 , first_order_entropy + 2) + ASCII_SPACE;
        width[6] = column_width(max_nb_segment - 1 , change_point_entropy + 2) + ASCII_SPACE;
        width[7] = stat_tool::column_width(max_nb_segment - 1 , uniform_entropy + 2) + ASCII_SPACE;
        width[8] = stat_tool::column_width(max_nb_segment - 1 , segmentation_divergence + 2) + ASCII_SPACE;
//        width[9] = column_width(max_nb_segment - 1 , marginal_entropy + 2) + ASCII_SPACE
        width[10] = stat_tool::column_width(nb_parameter[max_nb_segment]) + ASCII_SPACE;
        width[11] = stat_tool::column_width(max_nb_segment , penalized_likelihood[0] + 1) + ASCII_SPACE;
        width[12] = stat_tool::column_width(max_nb_segment , weight[0] + 1) + ASCII_SPACE;

        if (!bayesian) {
          if (max_nb_segment - min_nb_segment >= SLOPE_NB_SEGMENT_RANGE) {
            width[13] = stat_tool::column_width(max_nb_segment , penalized_likelihood[1] + 1) + ASCII_SPACE;
            width[14] = stat_tool::column_width(max_nb_segment , weight[1] + 1) + ASCII_SPACE;
          }
          width[15] = stat_tool::column_width(max_nb_segment , penalized_likelihood[2] + 1) + ASCII_SPACE;
          width[16] = stat_tool::column_width(max_nb_segment , weight[2] + 1) + ASCII_SPACE;
          if (max_nb_segment - min_nb_segment >= SLOPE_NB_SEGMENT_RANGE) {
            width[17] = stat_tool::column_width(max_nb_segment , penalized_likelihood[3] + 1) + ASCII_SPACE;
            width[18] = stat_tool::column_width(max_nb_segment , weight[3] + 1) + ASCII_SPACE;
          }
        }

        os << "\n" << SEQ_label[SEQL_NB_SEGMENT] << " | 2 * " << STAT_label[STATL_LIKELIHOOD]
           << " | 2 * " << SEQ_label[SEQL_POSSIBLE_SEGMENTATION_LIKELIHOOD]
           << " | " << SEQ_label[SEQL_POSTERIOR_PROBABILITY]
           << " | " << SEQ_label[SEQL_SEGMENTATION_ENTROPY]
           << " | " << SEQ_label[SEQL_FIRST_ORDER_ENTROPY]
           << " | " << SEQ_label[SEQL_CHANGE_POINT_ENTROPY]
           << " | " << SEQ_label[SEQL_UNIFORM_ENTROPY]
           << " | " << SEQ_label[SEQL_SEGMENTATION_DIVERGENCE]
//           << " | " << SEQ_label[SEQL_MARGINAL_ENTROPY]
           << " | " << STAT_label[STATL_FREE_PARAMETERS]
           << " | "  << STAT_criterion_word[ICL] << " - "  << STAT_label[STATL_WEIGHT];

        if (!bayesian) {
          if (max_nb_segment - min_nb_segment >= SLOPE_NB_SEGMENT_RANGE) {
            os << " | " << STAT_criterion_word[LIKELIHOOD_SLOPE] << " - "  << STAT_label[STATL_WEIGHT];
          }
          os << " | " << STAT_criterion_word[mBIC] << " - "  << STAT_label[STATL_WEIGHT];
          if (max_nb_segment - min_nb_segment >= SLOPE_NB_SEGMENT_RANGE) {
            os << " | " << STAT_criterion_word[SEGMENTATION_LIKELIHOOD_SLOPE] << " - "  << STAT_label[STATL_WEIGHT];
          }
        }
        os << endl;

        os.setf(ios::left , ios::adjustfield);

        os << setw(width[0]) << 1
           << setw(width[1]) << 2 * segmentation_likelihood[1]
           << setw(width[2]) << 2 * likelihood[1]
           << setw(width[3]) << posterior_probability[1]
           << setw(width[4]) << " "
           << setw(width[5]) << " "
           << setw(width[6]) << " "
           << setw(width[7]) << " "
           << setw(width[8]) << segmentation_divergence[1]
//           << setw(width[9]) << " "
           << setw(width[10]) << nb_parameter[1]
           << setw(width[11]) << penalized_likelihood[0][1]
           << setw(width[12]) << weight[0][1];

        if (!bayesian) {
          if (max_nb_segment - min_nb_segment >= SLOPE_NB_SEGMENT_RANGE) {
            os << setw(width[13]) << penalized_likelihood[1][1]
               << setw(width[14]) << weight[1][1];
          }
          os << setw(width[15]) << penalized_likelihood[2][1]
             << setw(width[16]) << weight[2][1];
          if (max_nb_segment - min_nb_segment >= SLOPE_NB_SEGMENT_RANGE) {
            os << setw(width[17]) << penalized_likelihood[3][1]
               << setw(width[18]) << weight[3][1];
          }
        }
        os << endl;

        for (i = 2;i <= max_nb_segment;i++) {
          os << setw(width[0]) << i
             << setw(width[1]) << 2 * segmentation_likelihood[i]
             << setw(width[2]) << 2 * likelihood[i]
             << setw(width[3]) << posterior_probability[i]
             << setw(width[4]) << segmentation_entropy[i]
             << setw(width[5]) << first_order_entropy[i]
             << setw(width[6]) << change_point_entropy[i]
             << setw(width[7]) << uniform_entropy[i]
             << setw(width[8]) << segmentation_divergence[i]
//             << setw(width[9]) << marginal_entropy[i]
             << setw(width[10]) << nb_parameter[i]
             << setw(width[11]) << penalized_likelihood[0][i]
             << setw(width[12]) << weight[0][i];

          if (!bayesian) {
            if (max_nb_segment - min_nb_segment >= SLOPE_NB_SEGMENT_RANGE) {
              os << setw(width[13]) << penalized_likelihood[1][i]
                 << setw(width[14]) << weight[1][i];
            }
            os << setw(width[15]) << penalized_likelihood[2][i]
               << setw(width[16]) << weight[2][i];
            if (max_nb_segment - min_nb_segment >= SLOPE_NB_SEGMENT_RANGE) {
              os << setw(width[17]) << penalized_likelihood[3][i]
                 << setw(width[18]) << weight[3][i];
            }
          }
          os << endl;
        }
        os << endl;
      }

      else {
        width[0] = stat_tool::column_width(max_nb_segment) + ASCII_SPACE;
        width[1] = stat_tool::column_width(max_nb_segment , segmentation_likelihood + 1 , 2.) + ASCII_SPACE;
        width[10] = stat_tool::column_width(nb_parameter[max_nb_segment]) + ASCII_SPACE;
        width[15] = stat_tool::column_width(max_nb_segment , penalized_likelihood[2] + 1) + ASCII_SPACE;
        width[16] = stat_tool::column_width(max_nb_segment , weight[2] + 1) + ASCII_SPACE;
        if (max_nb_segment - min_nb_segment >= SLOPE_NB_SEGMENT_RANGE) {
          width[17] = stat_tool::column_width(max_nb_segment , penalized_likelihood[3] + 1) + ASCII_SPACE;
          width[18] = stat_tool::column_width(max_nb_segment , weight[3] + 1) + ASCII_SPACE;
        }

        os << "\n" << SEQ_label[SEQL_NB_SEGMENT] << " | 2 * " << STAT_label[STATL_LIKELIHOOD]
           << " | " << STAT_label[STATL_FREE_PARAMETERS]
           << " | " << STAT_criterion_word[mBIC] << " - "  << STAT_label[STATL_WEIGHT];
        if (max_nb_segment - min_nb_segment >= SLOPE_NB_SEGMENT_RANGE) {
          os << " | " << STAT_criterion_word[SEGMENTATION_LIKELIHOOD_SLOPE] << " - "  << STAT_label[STATL_WEIGHT];
        }
        os << endl;

        os.setf(ios::left , ios::adjustfield);

        os << setw(width[0]) << 1
           << setw(width[1]) << 2 * segmentation_likelihood[1]
           << setw(width[10]) << nb_parameter[1]
           << setw(width[15]) << penalized_likelihood[2][1]
           << setw(width[16]) << weight[2][1];
        if (max_nb_segment - min_nb_segment >= SLOPE_NB_SEGMENT_RANGE) {
          os << setw(width[17]) << penalized_likelihood[3][1]
             << setw(width[18]) << weight[3][1];
        }
        os << endl;

        for (i = 2;i <= max_nb_segment;i++) {
          os << setw(width[0]) << i
             << setw(width[1]) << 2 * segmentation_likelihood[i]
             << setw(width[10]) << nb_parameter[i]
             << setw(width[15]) << penalized_likelihood[2][i]
             << setw(width[16]) << weight[2][i];
          if (max_nb_segment - min_nb_segment >= SLOPE_NB_SEGMENT_RANGE) {
            os << setw(width[17]) << penalized_likelihood[3][i]
               << setw(width[18]) << weight[3][i];
          }
          os << endl;
        }
        os << endl;
      }

      os.setf((FMTFLAGS)old_adjust , ios::adjustfield);

/*      if (!bayesian) {
        penalized_likelihood[0][1] = 2 * likelihood[1] - nb_parameter[1] *
                                     log((double)seq->length[0]);
//        penalized_likelihood[0][1] += nb_parameter[1] *
//                                      log((double)((seq->nb_variable - 2) * seq->length[0]));
        for (i = 2;i <= max_nb_segment;i++) {
          penalized_likelihood[0][i] = 2 * (likelihood[i] - segmentation_entropy[i]) -
                                       nb_parameter[i] * log((double)seq->length[0]);
//          penalized_likelihood[0][i] += nb_parameter[i] *
//                                        log((double)((seq->nb_variable - 2) * seq->length[0]));
        }

        max_likelihood[0] = penalized_likelihood[0][1];
        for (i = 2;i <= max_nb_segment;i++) {
          if (penalized_likelihood[0][i] > max_likelihood[0]) {
            max_likelihood[0] = penalized_likelihood[0][i];
          }
        }

        norm = 0.;
        for (i = 1;i <= max_nb_segment;i++) {
          weight[0][i] = exp((penalized_likelihood[0][i] - max_likelihood[0]) / 2);
          norm += weight[0][i];
        }
        for (i = 1;i <= max_nb_segment;i++) {
          weight[0][i] /= norm;
        }

        for (i = 1;i <= max_nb_segment;i++) {
          os << "\n" << i << " " << (i == 1 ? SEQ_label[SEQL_SEGMENT] : SEQ_label[SEQL_SEGMENTS])
             << "   2 * " << STAT_label[STATL_LIKELIHOOD] << ": " << 2 * segmentation_likelihood[i] << "   "
             << nb_parameter[i] << " " << STAT_label[nb_parameter[i] == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS]
             << "   2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << STAT_criterion_word[ICL] << "): "
             << penalized_likelihood[0][i] << "   " << STAT_label[STATL_WEIGHT] << ": " << weight[0][i] << endl;
        }
        os << endl;
      } */

      if ((model_type[0] != MEAN_CHANGE) && (model_type[0] != INTERCEPT_SLOPE_CHANGE)) {
        delete [] posterior_probability;
      }

      switch (bayesian) {

      case false : {
        for (i = 0;i < 4;i++) {
          delete [] weight[i];
        }
        break;
      }

      case true : {
        delete [] weight[0];
        break;
      }
      }

      delete [] weight;
#     endif

      inb_segment[0] = nb_segment;

      if (nb_segment == 1) {
        seq->one_segment_likelihood(0 , model_type , rank);
        seq->min_value[0] = 0;
        seq->max_value[0] = 0;
        seq->build_marginal_frequency_distribution(0);
      }

      else {
        seq->segmentation(inb_segment , model_type , rank);
      }

      switch (output) {

      case SEQUENCE : {
        oseq = seq->segmentation_output(inb_segment , model_type , os);
        break;
      }

      case SEGMENTATION_ENTROPY : {
        for (i = 0;i < 4;i++) {
          ilength[i] = max_nb_segment - 1;
        }
        itype[0] = REAL_VALUE;

        oseq = new Sequences(4 , NULL , ilength , NULL , TIME , 1 , itype);

        for (i = 2;i <= max_nb_segment;i++) {
          oseq->index_parameter[0][i - 2] = i;
          oseq->real_sequence[0][0][i - 2] = segmentation_entropy[i];
          oseq->index_parameter[1][i - 2] = i;
          oseq->real_sequence[1][0][i - 2] = first_order_entropy[i];
          oseq->index_parameter[2][i - 2] = i;
          oseq->real_sequence[2][0][i - 2] = change_point_entropy[i];
          oseq->index_parameter[3][i - 2] = i;
          oseq->real_sequence[3][0][i - 2] = uniform_entropy[i];
        }

        oseq->build_index_parameter_frequency_distribution();
        oseq->index_interval_computation();

        oseq->min_value_computation(0);
        oseq->max_value_computation(0);

        oseq->build_marginal_histogram(0);
        break;
      }

      case SEGMENTATION_DIVERGENCE : {
        ilength[0] = max_nb_segment;
        itype[0] = REAL_VALUE;

        oseq = new Sequences(1 , NULL , ilength , NULL , TIME , 1 , itype);

        for (i = 1;i <= max_nb_segment;i++) {
          oseq->index_parameter[0][i - 1] = i;
          oseq->real_sequence[0][0][i - 1] = segmentation_divergence[i];
        }

        oseq->build_index_parameter_frequency_distribution();
        oseq->index_interval_computation();

        oseq->min_value_computation(0);
        oseq->max_value_computation(0);

        oseq->build_marginal_histogram(0);
        break;
      }

      case LOG_LIKELIHOOD_SLOPE : {

#       ifdef DEBUG
        os << "\n";
        for (i = 1;i <= max_nb_segment;i++) {
          os << i << "\t" << penalty_shape[i];
          if ((model_type[0] != MEAN_CHANGE) && (model_type[0] != INTERCEPT_SLOPE_CHANGE)) {
            os << "\t" << likelihood[i] << "\t" << intercept + slope * penalty_shape[i];
          }
          os << "\t" << segmentation_likelihood[i] << "\t" << segmentation_intercept + segmentation_slope * penalty_shape[i] << endl;
        }
#       endif

        if ((model_type[0] != MEAN_CHANGE) && (model_type[0] != INTERCEPT_SLOPE_CHANGE)) {
          inb_sequence = 2;
        }
        else {
          inb_sequence = 1;
        }

        for (i = 0;i < inb_sequence;i++) {
          ilength[i] = max_nb_segment;
        }
        itype[0] = REAL_VALUE;
        itype[1] = AUXILIARY;

        oseq = new Sequences(inb_sequence , NULL , ilength , NULL , TIME , 2 , itype);

        switch (penalty_shape_type) {
        case 0 :
          scaling_factor = 1;
          break;
        case 1 :
          scaling_factor = seq->length[0];
          break;
        case 2 :
          scaling_factor = PENALTY_SHAPE_SCALING_FACTOR;
          break;
        case 3 :
          scaling_factor = seq->length[0];
          break;
        case 4 :
          scaling_factor = PENALTY_SHAPE_SCALING_FACTOR;
          break;
        }

        for (i = 1;i <= max_nb_segment;i++) {
          oseq->index_parameter[0][i - 1] = (int)::round(penalty_shape[i] * scaling_factor);
          oseq->real_sequence[0][0][i - 1] = segmentation_likelihood[i];
          oseq->real_sequence[0][1][i - 1] = segmentation_intercept + segmentation_slope * penalty_shape[i];
          if ((model_type[0] != MEAN_CHANGE) && (model_type[0] != INTERCEPT_SLOPE_CHANGE)) {
            oseq->index_parameter[1][i - 1] = (int)::round(penalty_shape[i] * scaling_factor);
            oseq->real_sequence[1][0][i - 1] = likelihood[i];
            oseq->real_sequence[1][1][i - 1] = intercept + slope * penalty_shape[i];
          }
        }

        oseq->build_index_parameter_frequency_distribution();
        oseq->index_interval_computation();

        oseq->min_value_computation(0);
        oseq->max_value_computation(0);
        oseq->min_value_computation(1);
        oseq->max_value_computation(1);

        oseq->build_marginal_histogram(0);
        break;
      }
      }

#     ifdef DEBUG
      if ((model_type[0] != MEAN_CHANGE) && (model_type[0] != MEAN_VARIANCE_CHANGE) &&
          (model_type[0] != INTERCEPT_SLOPE_CHANGE)) {
        hierarchical_segmentation(error , os , iidentifier , max_nb_segment , model_type);
      }
#     endif

    }

    else {
      oseq = NULL;
      error.update(SEQ_error[SEQR_SEGMENTATION_FAILURE]);
    }

    if (max_nb_segment - min_nb_segment >= SLOPE_NB_SEGMENT_RANGE) {
      delete [] penalty_shape;
    }

    for (i = 1;i < seq->nb_variable;i++) {
      delete [] rank[i];
    }
    delete [] rank;

    delete seq;

    delete [] segmentation_likelihood;
    delete [] nb_parameter;

    switch (bayesian) {

    case false : {
      delete [] segment_penalty;

      for (i = 0;i < 4;i++) {
        delete [] penalized_likelihood[i];
      }
      break;
    }

    case true : {
      delete [] penalized_likelihood[0];
      break;
    }
    }

    delete [] penalized_likelihood;

    if ((model_type[0] != MEAN_CHANGE) && (model_type[0] != INTERCEPT_SLOPE_CHANGE)) {
      delete [] likelihood;
      delete [] segmentation_entropy;
      delete [] first_order_entropy;
      delete [] change_point_entropy;
      delete [] uniform_entropy;
      delete [] segmentation_divergence;
      delete [] marginal_entropy;
    }
  }

  return oseq;
}


};  // namespace sequence_analysis
