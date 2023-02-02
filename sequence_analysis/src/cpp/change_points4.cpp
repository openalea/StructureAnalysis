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
 *       $Id: change_points4.cpp 11914 2012-03-26 06:29:13Z guedon $
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

#include <string>
#include <sstream>
#include <iomanip>

#include "sequences.h"
#include "sequence_label.h"

using namespace std;
using namespace stat_tool;


namespace sequence_analysis {


extern double log_factorial(int value);
extern double log_binomial_coefficient(int inf_bound , double parameter , int value);


#if defined (SYSTEM_IS__CYGWIN)
#define expl exp
#endif



/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the N most probable segmentations of a single sequence or
 *         a sample of sequences.
 *
 *  \param[in] index            sequence index,
 *  \param[in] nb_segment       number of segments,
 *  \param[in] model_type       segment model types,
 *  \param[in] common_contrast  flag contrast functions common to the individuals,
 *  \param[in] shape_parameter  negative binomial shape parameters,
 *  \param[in] irank            ranks (for ordinal variables),
 *  \param[in] os               stream,
 *  \param[in] format           file format (ASCII/SPREADSHEET),
 *  \param[in] inb_segmentation number of segmentations,
 *  \param[in] likelihood       log-likelihood of the multiple change-point model.
 *
 *  \return                     log-likelihood of the optimal segmentation.
 */
/*--------------------------------------------------------------*/

double Sequences::N_segmentation(int index , int nb_segment , segment_model *model_type ,
                                 bool common_contrast , double *shape_parameter ,
                                 double **irank , ostream &os , output_format format ,
                                 int inb_segmentation , double likelihood) const

{
  bool **active_cell;
  int i , j , k , m;
  int seq_length , brank , previous_rank , count , nb_cell , *inf_bound_parameter ,
      *seq_index_parameter , *rank , *change_point , *psegment , ***optimal_length , ***optimal_rank;
  double buff , segmentation_likelihood , *nb_segmentation , **hyperparam , **seq_mean ,
         **nb_segmentation_forward , ***factorial , ***binomial_coeff , ***forward , ***mean ,
         ***variance , ***intercept , ***slope , ***autoregressive_coeff;
  long double *contrast , likelihood_cumul , posterior_probability_cumul;

# ifdef MESSAGE
  int *first_change_point;
  long double norm;
# endif


  // initializations

  factorial = new double**[nb_variable];
  inf_bound_parameter = new int[nb_variable];
  binomial_coeff = new double**[nb_variable];
  seq_mean = new double*[nb_variable];
  hyperparam = new double*[nb_variable];

  for (i = 1;i < nb_variable;i++) {

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

  if (index_param_type == IMPLICIT_TYPE) {
    seq_index_parameter = new int[seq_length];
    for (j = 0;j < seq_length;j++) {
      seq_index_parameter[j] = j;
    }
  }
  else {
    seq_index_parameter = index_parameter[index == I_DEFAULT ? 0 : index];
  }

  contrast = new long double[seq_length];

  nb_segmentation_forward = new double*[seq_length];
  for (i = 0;i < seq_length;i++) {
    nb_segmentation_forward[i] = new double[nb_segment];
  }

  forward = new double**[seq_length];
  for (i = 0;i < seq_length;i++) {
    forward[i] = new double*[nb_segment];
    for (j = 0;j < nb_segment;j++) {
      forward[i][j] = new double[inb_segmentation];
    }
  }

  nb_segmentation = new double[nb_segment];
  rank = new int[seq_length + 1];

  optimal_length = new int**[seq_length];
  for (i = 0;i < seq_length;i++) {
    optimal_length[i] = new int*[nb_segment];
    for (j = 0;j < nb_segment;j++) {
      optimal_length[i][j] = new int[inb_segmentation];
    }
  }

  optimal_rank = new int**[seq_length];
  for (i = 0;i < seq_length;i++) {
    optimal_rank[i] = new int*[nb_segment];
    for (j = 0;j < nb_segment;j++) {
      optimal_rank[i][j] = new int[inb_segmentation];
    }
  }

  change_point = new int[nb_segment + 1];

  mean = new double**[nb_variable];
  variance = new double**[nb_variable];
  intercept = new double**[nb_variable];
  slope = new double**[nb_variable];
  autoregressive_coeff = new double**[nb_variable];

  for (i = 1;i < nb_variable;i++) {
    if ((model_type[i - 1] == POISSON_CHANGE) || (model_type[i - 1] == NEGATIVE_BINOMIAL_0_CHANGE) ||
        (model_type[i - 1] == NEGATIVE_BINOMIAL_1_CHANGE) || (model_type[i - 1] == GAUSSIAN_CHANGE) ||
        (model_type[0] == MEAN_CHANGE) || (model_type[i - 1] == VARIANCE_CHANGE) ||
        (model_type[i - 1] == BAYESIAN_POISSON_CHANGE) || (model_type[i - 1] == BAYESIAN_GAUSSIAN_CHANGE)) {
      if ((index != I_DEFAULT) || (!common_contrast)) {
        mean[i] = new double*[nb_sequence];
        variance[i] = new double*[nb_sequence];

        for (j = 0;j < nb_sequence;j++) {
          if ((index == I_DEFAULT) || (index == j)) {
            mean[i][j] = new double[nb_segment];
            variance[i][j] = new double[nb_segment];
          }
          else {
            mean[i][j] = NULL;
            variance[i][j] = NULL;
          }
        }
      }

      else {
        mean[i] = new double*[1];
        mean[i][0] = new double[nb_segment];
        variance[i] = new double*[1];
        variance[i][0] = new double[nb_segment];
      }
    }

    else if ((model_type[i - 1] == LINEAR_MODEL_CHANGE) || (model_type[0] == INTERCEPT_SLOPE_CHANGE)) {
      if ((index != I_DEFAULT) || (!common_contrast)) {
        intercept[i] = new double*[nb_sequence];
        slope[i] = new double*[nb_sequence];
        variance[i] = new double*[nb_sequence];

        for (j = 0;j < nb_sequence;j++) {
          if ((index == I_DEFAULT) || (index == j)) {
            intercept[i][j] = new double[nb_segment];
            slope[i][j] = new double[nb_segment];
            variance[i][j] = new double[nb_segment];
          }
          else {
            intercept[i][j] = NULL;
            slope[i][j] = NULL;
            variance[i][j] = NULL;
          }
        }
      }

      else {
        intercept[i] = new double*[1];
        intercept[i][0] = new double[nb_segment];
        slope[i] = new double*[1];
        slope[i][0] = new double[nb_segment];
        variance[i] = new double*[1];
        variance[i][0] = new double[nb_segment];
      }
    }

    else if ((model_type[i - 1] == AUTOREGRESSIVE_MODEL_CHANGE) || (model_type[i - 1] == STATIONARY_AUTOREGRESSIVE_MODEL_CHANGE)) {
      if ((index != I_DEFAULT) || (!common_contrast)) {
        mean[i] = new double*[nb_sequence];
        autoregressive_coeff[i] = new double*[nb_sequence];
        variance[i] = new double*[nb_sequence];

        for (j = 0;j < nb_sequence;j++) {
          if ((index == I_DEFAULT) || (index == j)) {
            mean[i][j] = new double[nb_segment];
            autoregressive_coeff[i][j] = new double[nb_segment];
            variance[i][j] = new double[nb_segment];
          }
          else {
            mean[i][j] = NULL;
            autoregressive_coeff[i][j] = NULL;
            variance[i][j] = NULL;
          }
        }
      }

      else {
        mean[i] = new double*[1];
        mean[i][0] = new double[nb_segment];
        autoregressive_coeff[i] = new double*[1];
        autoregressive_coeff[i][0] = new double[nb_segment];
        variance[i] = new double*[1];
        variance[i][0] = new double[nb_segment];
      }
    }
  }

  active_cell = new bool*[seq_length];
  for (i = 0;i < seq_length;i++) {
    active_cell[i] = new bool[nb_segment];
    for (j = 0;j < nb_segment;j++) {
      active_cell[i][j] = false;
    }
  }

# ifdef DEBUG
  double **segment_probability;

  segment_probability = new double*[seq_length];
  for (i = 0;i < seq_length;i++) {
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
    smoothed_probability = new double*[seq_length];
    for (i = 0;i < seq_length;i++) {
      smoothed_probability[i] = new double[nb_segment];
      for (j = 0;j < nb_segment;j++) {
        smoothed_probability[i][j] = 0.;
      }
    }
  }
# endif

# ifdef DEBUG
  for (i = 0;i < nb_segment;i++) {
    nb_segmentation[i] = 1;
  }
# endif

  // forward recurrence

  for (i = 0;i < seq_length;i++) {

    // computation of segment contrast functions (log-likelihoods or sum of squared deviations)

    forward_contrast(i , index , model_type , common_contrast , factorial ,
                     shape_parameter , binomial_coeff , seq_mean , seq_index_parameter ,
                     hyperparam , irank , contrast);

#   ifdef DEBUG
    for (j = i - 1;j >= 0;j--) {
      cout << contrast[j] << "  ";
    }
    cout << endl;
#   endif

    // computation of the number of segmentations

    for (j = 0;j < nb_segment;j++) {
      nb_segmentation_forward[i][j] = 0;
    }

    for (j = MAX(0 , nb_segment + i - seq_length);j < MIN((i < seq_length - 1 ? nb_segment - 1 : nb_segment) , i + 1);j++) {
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
    for (j = 1;j < MIN((i < seq_length - 1 ? nb_segment - 1 : nb_segment) , i + 1);j++) {
      nb_segmentation[j] = nb_segmentation[j - 1] * (i - j + 1) / j;
    }

    if (i < inb_segmentation) {
      cout << i << ": ";
      for (j = 0;j < MIN((i < seq_length - 1 ? nb_segment - 1 : nb_segment) , i + 1);j++) {
        cout << nb_segmentation_forward[i][j] << " " << nb_segmentation[j] << " | ";
      }
      cout << endl;
    }
#   endif

    for (j = 0;j < MIN((i < seq_length - 1 ? nb_segment - 1 : nb_segment) , i + 1);j++) {
      nb_segmentation[j] = nb_segmentation_forward[i][j];
      if (nb_segmentation[j] > inb_segmentation) {
        nb_segmentation[j] = inb_segmentation;
      }
    }

    for (j = MAX(0 , nb_segment + i - seq_length);j < MIN((i < seq_length - 1 ? nb_segment - 1 : nb_segment) , i + 1);j++) {
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
            cout << "\nuseful test" << endl;
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
    for (j = MAX(0 , nb_segment + i - seq_length);j < MIN((i < seq_length - 1 ? nb_segment - 1 : nb_segment) , i + 1);j++) {
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
  streamsize nb_digits;

  buff = 1.;
  for (i = 1;i < nb_segment;i++) {
    buff *= (double)(seq_length - i) / (double)i;
//    buff = buff * (seq_length - i) / i;
  }

  nb_digits = os.precision(10);

  os << "\n" << SEQ_label[SEQL_NB_SEGMENTATION] << ": "
     << nb_segmentation_forward[seq_length - 1][nb_segment - 1] << " (" << buff << ")" << endl;

  os.precision(nb_digits);

  if (((model_type[0] == MEAN_CHANGE) || (model_type[0] == INTERCEPT_SLOPE_CHANGE)) &&
      (format == SPREADSHEET) && (nb_segment == 2) && (inb_segmentation >= seq_length - 1)) {
    first_change_point = new int[seq_length];
  }
# endif

  // restoration

  likelihood_cumul = 0.;
  posterior_probability_cumul = 0.;

  for (i = 0;i < nb_segmentation[nb_segment - 1];i++) {
    if (forward[seq_length - 1][nb_segment - 1][i] == D_INF) {
      break;
    }

#   ifdef DEBUG
    cout << "\n";
#   endif

    j = seq_length - 1;
    change_point[nb_segment] = seq_length;
    psegment = int_sequence[index == I_DEFAULT ? 0 : index][0] + j;
    brank = i;

    for (k = nb_segment - 1;k >= 0;k--) {
      for (m = j;m > j - optimal_length[j][k][brank];m--) {
        active_cell[m][k] = true;
        *psegment-- = k;
      }

#     ifdef DEBUG
      cout << k << " " << optimal_length[j][k][brank] << " " << brank << " | ";
#     endif

      if (k > 0) {
        previous_rank = optimal_rank[j][k][brank];
      }
      j -= optimal_length[j][k][brank];
      change_point[k] = j + 1;
      if (k > 0) {
        brank = previous_rank;
      }
    }

#   ifdef DEBUG
    cout << endl;
#   endif

    for (j = 1;j < nb_variable;j++) {
      piecewise_linear_function(index , j , nb_segment , model_type[j - 1] , common_contrast ,
                                change_point , seq_index_parameter , NULL , mean[j] , variance[j] ,
                                NULL , intercept[j] , slope[j] , autoregressive_coeff[j]);
    }

    if ((model_type[0] == MEAN_CHANGE) || (model_type[0] == INTERCEPT_SLOPE_CHANGE)) {
      if (forward[seq_length - 1][nb_segment - 1][i] < 0.) {
        count = (index == I_DEFAULT ? nb_sequence : 1);

        forward[seq_length - 1][nb_segment - 1][i] = -((double)(count * seq_length) / 2.) *
                                                      (log(-forward[seq_length - 1][nb_segment - 1][i] /
                                                        (count * seq_length)) + log(2 * M_PI) + 1);
/*        forward[seq_length - 1][nb_segment - 1][i] = -((double)(count * seq_length) / 2.) *
                                                      (log(-forward[seq_length - 1][nb_segment - 1][i] /
                                                        (count * (seq_length - nb_segment))) + log(2 * M_PI)) -
                                                       (double)(count * (seq_length - nb_segment)) / 2.; */
      }
      else {
        forward[seq_length - 1][nb_segment - 1][i] = D_INF;
      }
    }

    if (i == 0) {
      segmentation_likelihood = forward[seq_length - 1][nb_segment - 1][i];
    }

    if (forward[seq_length - 1][nb_segment - 1][i] != D_INF) {
      likelihood_cumul += exp(forward[seq_length - 1][nb_segment - 1][i]);
      if (likelihood != D_INF) {
        posterior_probability_cumul += exp(forward[seq_length - 1][nb_segment - 1][i] - likelihood);
      }
    }

#   ifdef DEBUG
    psegment = int_sequence[index == I_DEFAULT ? 0 : index][0];
    for (j = 0;j < seq_length;j++) {
//      if (((i == 0) || (*psegment != *(psegment - 1))) &&
//          (forward[seq_length - 1][nb_segment - 1][i] > segment_probability[j][*psegment])) {
      if (forward[seq_length - 1][nb_segment - 1][i] > segment_probability[j][*psegment]) {
        segment_probability[j][*psegment] = forward[seq_length - 1][nb_segment - 1][i];
      }
      psegment++;
    }
#   endif

#   ifdef MESSAGE
    if (inb_segmentation >= 1000) {

      // approximation of smoothed probabilities

      buff = exp(forward[seq_length - 1][nb_segment - 1][i]);
      approximated_likelihood += buff;
      psegment = int_sequence[index == I_DEFAULT ? 0 : index][0];
      for (j = 0;j < seq_length;j++) {
        smoothed_probability[j][*psegment++] += buff;
      }
    }
#   endif

    nb_cell = 0;
    for (j = 0;j < seq_length;j++) {
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

    case ASCII : {
      if (inb_segmentation <= 200) {
      psegment = int_sequence[index == I_DEFAULT ? 0 : index][0];
      for (j = 0;j < seq_length;j++) {
        os << *psegment++ << " ";
      }

      os << "  " << i + 1 << "  " << forward[seq_length - 1][nb_segment - 1][i] << "   (";
      if (likelihood != D_INF) {
        os << exp(forward[seq_length - 1][nb_segment - 1][i] - likelihood) << "  ";
        if (boost::math::isnan(likelihood_cumul)) {
          os << likelihood_cumul / exp(likelihood);
        }
        else {
          os << posterior_probability_cumul;
        }
      }
      else {
        os << exp(forward[seq_length - 1][nb_segment - 1][i] - segmentation_likelihood);
      }
      os << "  " << nb_cell << ")" << endl;

      os << (nb_segment == 2 ? SEQ_label[SEQL_CHANGE_POINT] : SEQ_label[SEQL_CHANGE_POINTS]) << ": ";

      for (j = 1;j < nb_segment;j++) {
        os << seq_index_parameter[change_point[j]];
        if (j < nb_segment - 1) {
          os << ", ";
        }
      }
      os << endl;

/*      psegment = int_sequence[index == I_DEFAULT ? 0 : index][0] + 1;
      for (j = 1;j < seq_length;j++) {
        if (*psegment != *(psegment - 1)) {
          os << seq_index_parameter[j] << ", ";
        }
        psegment++;
      }
      os << endl; */

      for (j = 1;j < nb_variable;j++) {
        piecewise_linear_function_ascii_print(os , index , j , nb_segment , model_type[j - 1] ,
                                              common_contrast , change_point , seq_index_parameter ,
                                              mean[j] , variance[j] , intercept[j] , slope[j] ,
                                              autoregressive_coeff[j]);
      }
      }
      break;
    }

    case SPREADSHEET : {
      psegment = int_sequence[index == I_DEFAULT ? 0 : index][0];
      for (j = 0;j < seq_length;j++) {
        os << *psegment++ << "\t";
      }

      os << "\t" << i + 1 << "\t" << forward[seq_length - 1][nb_segment - 1][i] << "\t";
      if (likelihood != D_INF) {
        os << exp(forward[seq_length - 1][nb_segment - 1][i] - likelihood) << "\t";
        if (boost::math::isnan(likelihood_cumul)) {
          os << likelihood_cumul / exp(likelihood);
        }
        else {
          os << posterior_probability_cumul;
        }
      }
      else {
        os << exp(forward[seq_length - 1][nb_segment - 1][i] - segmentation_likelihood);
      }
      os << "\t" << nb_cell << endl;

      for (j = 1;j < nb_variable;j++) {
        piecewise_linear_function_spreadsheet_print(os , index , j , nb_segment , model_type[j - 1] ,
                                                    common_contrast , change_point , seq_index_parameter ,
                                                    mean[j] , variance[j] , intercept[j] , slope[j] ,
                                                    autoregressive_coeff[j]);
      }

      if (((model_type[0] == MEAN_CHANGE) || (model_type[0] == INTERCEPT_SLOPE_CHANGE)) &&
          (nb_segment == 2) && (inb_segmentation >= seq_length - 1)) {
       first_change_point[i] = change_point[1];
      }
      break;
    }
    }

#   endif

  }

# ifdef MESSAGE
  if (((model_type[0] == MEAN_CHANGE) || (model_type[0] == INTERCEPT_SLOPE_CHANGE)) &&
      (format == SPREADSHEET) && (nb_segment == 2) && (inb_segmentation >= seq_length - 1)) {
    norm = 0.;
    for (i = 0;i < seq_length - 1;i++) {
      norm += exp(forward[seq_length - 1][nb_segment - 1][i]);
    }

    os << "\n" << SEQ_label[SEQL_POSTERIOR_CHANGE_POINT_PROBABILITY] << "\n\n";

    for (i = 1;i < seq_length;i++) {
      for (j = 0;j < seq_length - 1;j++) {
        if (first_change_point[j] == i) {
          os << seq_index_parameter[i] << "\t"
             << exp(forward[seq_length - 1][nb_segment - 1][j]) / norm << endl;
          break;
        }
      }
    }

    delete [] first_change_point;
  }

  if (((model_type[0] == MEAN_CHANGE) || (model_type[0] == INTERCEPT_SLOPE_CHANGE)) &&
       (nb_segment > 2) && (inb_segmentation >= nb_segmentation_forward[seq_length - 1][nb_segment - 1])) {
    norm = 0.;
    for (i = 0;i < nb_segmentation_forward[seq_length - 1][nb_segment - 1];i++) {
      norm += exp(forward[seq_length - 1][nb_segment - 1][i]);
    }

    os << SEQ_label[SEQL_POSTERIOR_PROBABILITY] << ": " << exp(forward[seq_length - 1][nb_segment - 1][0]) / norm << endl;
  }
# endif

# ifdef DEBUG
  if (((likelihood != D_INF) && (likelihood_cumul / exp(likelihood) > 0.8)) ||
      (segmentation_likelihood != D_INF)) {
    if (likelihood != D_INF) {
      for (i = 0;i < seq_length;i++) {
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
      for (i = 0;i < seq_length;i++) {
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

    psegment = int_sequence[index == I_DEFAULT ? 0 : index][0];
    for (j = 0;j < seq_length;j++) {
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


    for (i = 0;i < seq_length;i++) {
      for (j = 0;j < nb_segment;j++) {
        smoothed_probability[i][j] /= approximated_likelihood;
      }
    }

    psegment = int_sequence[index == I_DEFAULT ? 0 : index][0];
    for (i = 0;i < seq_length;i++) {
      *psegment++ = I_DEFAULT;
    }

    os << "\n" << SEQ_label[SEQL_POSTERIOR_SEGMENT_PROBABILITY] << "\n\n";

    profile_ascii_print(os , index , nb_segment , smoothed_probability , SEQ_label[SEQL_SEGMENT]);

    // approximation of the Kullback-Leibler divergence of the uniform distribution from the
    // segmentation distribution

    segmentation = new Distribution(inb_segmentation);
    likelihood_cumul = 0.;
    divergence = 0.;

    if (likelihood != D_INF) {
      for (i = 0;i < inb_segmentation;i++) {
        segmentation->mass[i] = exp(forward[seq_length - 1][nb_segment - 1][i] - likelihood);
        likelihood_cumul += exp(forward[seq_length - 1][nb_segment - 1][i] - likelihood);
        segmentation->cumul[i] = likelihood_cumul;

        divergence += exp(forward[seq_length - 1][nb_segment - 1][i] - likelihood) *
                      (forward[seq_length - 1][nb_segment - 1][i] - likelihood +
                       log(nb_segmentation_forward[seq_length - 1][nb_segment - 1]));
      }

      segmentation->complement = 1. - likelihood_cumul;
    }

    else {
      for (i = 0;i < inb_segmentation;i++) {
        segmentation->mass[i] = exp(forward[seq_length - 1][nb_segment - 1][i]) / approximated_likelihood;

        divergence += exp(forward[seq_length - 1][nb_segment - 1][i]) / approximated_likelihood *
                      (forward[seq_length - 1][nb_segment - 1][i] - log(approximated_likelihood) +
                       log(nb_segmentation_forward[seq_length - 1][nb_segment - 1]));
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
        likelihood_cumul += exp(forward[seq_length - 1][nb_segment - 1][j]);
        if (likelihood_cumul / exp(likelihood) > cdf[i]) {
          os << j << " " << previous_cumul[0] / exp(likelihood) << " "
             << likelihood_cumul / exp(likelihood) << " (";
//          os << j + 1 << " " << likelihood_cumul / exp(likelihood) << " (";
          if (i == 0) {
            os << (j + 1) / nb_segmentation_forward[seq_length - 1][nb_segment - 1] << ")" << endl;
          }
          else {
            os << likelihood_cumul / exp(likelihood) - previous_cumul[1] << " "
               << (j + 1) / nb_segmentation_forward[seq_length - 1][nb_segment - 1] - previous_cumul[2] << ")" << endl;
          }

          if (cdf[i] == 1) {
            break;
          }
          previous_cumul[1] = likelihood_cumul / exp(likelihood);
          previous_cumul[2] = (j + 1) / nb_segmentation_forward[seq_length - 1][nb_segment - 1];
          i++;
        }
      }
    }

    else {
      for (j = 0;j < inb_segmentation;j++) {
        previous_cumul[0] = likelihood_cumul;
        likelihood_cumul += exp(forward[seq_length - 1][nb_segment - 1][j]);
        if (likelihood_cumul / approximated_likelihood > cdf[i]) {
          os << j << " " << previous_cumul[0] / approximated_likelihood << " "
             << likelihood_cumul / approximated_likelihood << " (";
//          os << j + 1 << " " << likelihood_cumul / approximated_likelihood << " (";
          if (i == 0) {
            os << (j + 1) / nb_segmentation_forward[seq_length - 1][nb_segment - 1] << ")" << endl;
          }
          else {
            os << likelihood_cumul / approximated_likelihood - previous_cumul[1] << " "
               << (j + 1) / nb_segmentation_forward[seq_length - 1][nb_segment - 1] - previous_cumul[2] << ")" << endl;
          }

          if (cdf[i] == 1) {
            break;
          }
          previous_cumul[1] = likelihood_cumul / approximated_likelihood;
          previous_cumul[2] = (j + 1) / nb_segmentation_forward[seq_length - 1][nb_segment - 1];
          i++;
        }
      }
    }

/*    ofstream out_file("Spreadsheet/segmentation_probability.xld");
  
    likelihood_cumul = 0.;
    if (likelihood != D_INF) {
      for (i = 0;i < inb_segmentation;i++) {
        likelihood_cumul += exp(forward[seq_length - 1][nb_segment - 1][i]);
        out_file << i + 1 << "\t" << exp(forward[seq_length - 1][nb_segment - 1][i] - likelihood) << "\t"
                 << likelihood_cumul / exp(likelihood) << "\t"
                 << 1. / nb_segmentation_forward[seq_length - 1][nb_segment - 1] << endl;
      }
    }

    else {
      for (i = 0;i < inb_segmentation;i++) {
        likelihood_cumul += exp(forward[seq_length - 1][nb_segment - 1][i]);
        out_file << i + 1 << "\t" << exp(forward[seq_length - 1][nb_segment - 1][i]) / approximated_likelihood << "\t"
                 << likelihood_cumul / approximated_likelihood << "\t"
                 << 1. / nb_segmentation_forward[seq_length - 1][nb_segment - 1] << endl;
      }
    } */
  }
# endif

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
    delete [] nb_segmentation_forward[i];
  }
  delete [] nb_segmentation_forward;

  for (i = 0;i < seq_length;i++) {
    for (j = 0;j < nb_segment;j++) {
      delete [] forward[i][j];
    }
    delete [] forward[i];
  }
  delete [] forward;

  delete [] nb_segmentation;
  delete [] rank;

  for (i = 0;i < seq_length;i++) {
    for (j = 0;j < nb_segment;j++) {
      delete [] optimal_length[i][j];
    }
    delete [] optimal_length[i];
  }
  delete [] optimal_length;

  for (i = 0;i < seq_length;i++) {
    for (j = 0;j < nb_segment;j++) {
      delete [] optimal_rank[i][j];
    }
    delete [] optimal_rank[i];
  }
  delete [] optimal_rank;

  delete [] change_point;

  for (i = 1;i < nb_variable;i++) {
    if ((model_type[i - 1] == POISSON_CHANGE) || (model_type[i - 1] == NEGATIVE_BINOMIAL_0_CHANGE) ||
        (model_type[i - 1] == NEGATIVE_BINOMIAL_1_CHANGE) || (model_type[i - 1] == GAUSSIAN_CHANGE) ||
        (model_type[0] == MEAN_CHANGE) || (model_type[i - 1] == VARIANCE_CHANGE) ||
        (model_type[i - 1] == BAYESIAN_POISSON_CHANGE) || (model_type[i - 1] == BAYESIAN_GAUSSIAN_CHANGE)) {
      if ((index != I_DEFAULT) || (!common_contrast)) {
        for (j = 0;j < nb_sequence;j++) {
          if ((index == I_DEFAULT) || (index == j)) {
            delete [] mean[i][j];
            delete [] variance[i][j];
          }
        }
      }

      else {
        delete [] mean[i][0];
        delete [] variance[i][0];
      }

      delete [] mean[i];
      delete [] variance[i];
    }

    else if ((model_type[i - 1] == LINEAR_MODEL_CHANGE) || (model_type[0] == INTERCEPT_SLOPE_CHANGE)) {
      if ((index != I_DEFAULT) || (!common_contrast)) {
        for (j = 0;j < nb_sequence;j++) {
          if ((index == I_DEFAULT) || (index == j)) {
            delete [] intercept[i][j];
            delete [] slope[i][j];
            delete [] variance[i][j];
          }
        }
      }

      else {
        delete [] intercept[i][0];
        delete [] slope[i][0];
        delete [] variance[i][0];
      }

      delete [] intercept[i];
      delete [] slope[i];
      delete [] variance[i];
    }

    else if ((model_type[i - 1] == AUTOREGRESSIVE_MODEL_CHANGE) || (model_type[i - 1] == STATIONARY_AUTOREGRESSIVE_MODEL_CHANGE)) {
      if ((index != I_DEFAULT) || (!common_contrast)) {
        for (j = 0;j < nb_sequence;j++) {
          if ((index == I_DEFAULT) || (index == j)) {
            delete [] mean[i][j];
            delete [] autoregressive_coeff[i][j];
            delete [] variance[i][j];
          }
        }
      }

      else {
        delete [] mean[i][0];
        delete [] autoregressive_coeff[i][0];
        delete [] variance[i][0];
      }

      delete [] mean[i];
      delete [] autoregressive_coeff[i];
      delete [] variance[i];
    }
  }

  delete [] mean;
  delete [] variance;
  delete [] intercept;
  delete [] slope;
  delete [] autoregressive_coeff;

  for (i = 0;i < seq_length;i++) {
    delete [] active_cell[i];
  }
  delete [] active_cell;

# ifdef DEBUG
  for (i = 0;i < seq_length;i++) {
    delete [] segment_probability[i];
  }
  delete [] segment_probability;
# endif

# ifdef MESSAGE
  if (inb_segmentation >= 1000) {
    for (i = 0;i < seq_length;i++) {
      delete [] smoothed_probability[i];
    }
    delete [] smoothed_probability;
  }
# endif

  return segmentation_likelihood;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation by maximization of segment or change-point profiles for
 *         a single sequence or a sample of sequences.
 *
 *  \param[in] index           sequence index,
 *  \param[in] nb_segment      number of segments,
 *  \param[in] model_type      segment model types,
 *  \param[in] common_contrast flag contrast functions common to the individuals,
 *  \param[in] shape_parameter negative binomial shape parameters,
 *  \param[in] rank            ranks (for ordinal variables),
 *  \param[in] os              stream,
 *  \param[in] plot_set        pointer on a MultiPlotSet object,
 *  \param[in] output          output type,
 *  \param[in] format          output format (ASCII/SPREADSHEET/GNUPLOT/PLOT),
 *  \param[in] likelihood      log-likelihood of the multiple change-point model.
 *
 *  \return                    log-likelihood of the optimal segmentation.
 */
/*--------------------------------------------------------------*/

double Sequences::forward_backward_dynamic_programming(int index , int nb_segment ,
                                                       segment_model *model_type , bool common_contrast ,
                                                       double *shape_parameter , double **rank ,
                                                       ostream *os , MultiPlotSet *plot_set ,
                                                       change_point_profile output , output_format format ,
                                                       double likelihood) const

{
  int i , j , k , m;
  int seq_length , count , *inf_bound_parameter , *seq_index_parameter , *change_point , *psegment ,
       **optimal_length;
  double buff , segmentation_likelihood , backward_max , **seq_mean , **hyperparam , **forward ,
         **backward , **backward_output , **output_piecewise_function , ***piecewise_function ,
         ***factorial , ***binomial_coeff;
  long double *contrast;


  factorial = new double**[nb_variable];
  inf_bound_parameter = new int[nb_variable];
  binomial_coeff = new double**[nb_variable];
  seq_mean = new double*[nb_variable];
  seq_index_parameter = NULL;
  hyperparam = new double*[nb_variable];

  for (i = 1;i < nb_variable;i++) {

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

  backward = new double*[seq_length];
  for (i = 0;i < seq_length;i++) {
    backward[i] = new double[nb_segment];
  }

  backward_output = new double*[seq_length];
  for (i = 0;i < seq_length;i++) {
    backward_output[i] = new double[nb_segment];
  }

  piecewise_function = new double**[nb_variable];
  for (i = 1;i < nb_variable;i++) {
    if ((model_type[i - 1] == POISSON_CHANGE) || (model_type[i - 1] == NEGATIVE_BINOMIAL_0_CHANGE) ||
        (model_type[i - 1] == NEGATIVE_BINOMIAL_1_CHANGE) || (model_type[i - 1] == GAUSSIAN_CHANGE) ||
        (model_type[0] == MEAN_CHANGE) || (model_type[i - 1] == VARIANCE_CHANGE) ||
        (model_type[i - 1] == LINEAR_MODEL_CHANGE) || (model_type[0] == INTERCEPT_SLOPE_CHANGE) ||
        (model_type[i - 1] == AUTOREGRESSIVE_MODEL_CHANGE) || (model_type[i - 1] == STATIONARY_AUTOREGRESSIVE_MODEL_CHANGE) ||
        (model_type[i - 1] == BAYESIAN_POISSON_CHANGE) || (model_type[i - 1] == BAYESIAN_GAUSSIAN_CHANGE)) {
      piecewise_function[i] = new double*[nb_sequence];
      for (j = 0;j < nb_sequence;j++) {
        if ((index == I_DEFAULT) || (index == j)) {
          piecewise_function[i][j] = new double[length[j]];
        }
        else {
          piecewise_function[i][j] = NULL;
        }
      }
    }

    else {
      piecewise_function[i] = NULL;
    }
  }

  // forward recurrence

  for (i = 0;i < seq_length;i++) {

    // computation of segment contrast functions (log-likelihoods or sum of squared deviations)

    forward_contrast(i , index , model_type , common_contrast , factorial ,
                     shape_parameter , binomial_coeff , seq_mean , seq_index_parameter ,
                     hyperparam , rank , contrast);

    for (j = 0;j < nb_segment;j++) {
      forward[i][j] = D_INF;
    }

    for (j = MAX(0 , nb_segment + i - seq_length);j < MIN((i < seq_length - 1 ? nb_segment - 1 : nb_segment) , i + 1);j++) {
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

  if (forward[seq_length - 1][nb_segment - 1] == D_INF) {
    segmentation_likelihood = D_INF;
  }

  else {

    // restoration

    change_point = new int[nb_segment + 1];
    i = seq_length - 1;
    change_point[nb_segment] = seq_length;
    psegment = int_sequence[index == I_DEFAULT ? 0 : index][0] + i;

    for (j = nb_segment - 1;j >= 0;j--) {
      for (k = i;k > i - optimal_length[i][j];k--) {
        *psegment-- = j;
      }
      i -= optimal_length[i][j];
      change_point[j] = i + 1;
    }

    if (index == I_DEFAULT) {
      for (i = 1;i < nb_sequence;i++) {
        for (j = 0;j < length[0];j++) {
          int_sequence[i][0][j] = int_sequence[0][0][j];
        }
      }
    }

    // backward recurrence

    for (i = seq_length - 1;i >= 0;i--) {

      // computation of segment contrast functions (log-likelihoods or sum of squared deviations)

      backward_contrast(i , index , model_type , common_contrast , factorial ,
                        shape_parameter , binomial_coeff , seq_mean , seq_index_parameter ,
                        hyperparam , rank , contrast);

      for (j = 0;j < nb_segment;j++) {
        backward_output[i][j] = D_INF;
      }

      for (j = MAX((i == 0 ? 0 : 1) , nb_segment + i - seq_length);j < MIN(nb_segment , i + 1);j++) {
        if (j < nb_segment - 1) {
          backward[i][j] = D_INF;
          for (k = seq_length + j - nb_segment;k >= i;k--) {
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
          backward[i][j] = contrast[seq_length - 1];

          if ((output == SEGMENT) && (forward[i - 1][j - 1] != D_INF) &&
              (backward[i][j] != D_INF)) {
            buff = forward[i - 1][j - 1] + backward[i][j];
            for (k = seq_length - 1;k > i;k--) {
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
    for (i = 1;i < seq_length;i++) {
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

    // restoration

#   ifdef MESSAGE
    if (output == SEGMENT) {
      int optimal_segment;

      psegment = int_sequence[index == I_DEFAULT ? 0 : index][0];

      for (i = 0;i < seq_length;i++) {
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
      piecewise_linear_function(index , i , nb_segment , model_type[i - 1] , common_contrast ,
                                change_point , seq_index_parameter , piecewise_function[i]);
    }

#   ifdef MESSAGE
    if ((backward[0][0] < forward[seq_length - 1][nb_segment - 1] - DOUBLE_ERROR) ||
        (backward[0][0] > forward[seq_length - 1][nb_segment - 1] + DOUBLE_ERROR)) {
      cout << "\nERROR: " << backward[0][0] << " | " << forward[seq_length - 1][nb_segment - 1] << endl;
    }
/*    if ((backward_output[0][0] < backward[0][0] - DOUBLE_ERROR) ||
        (backward_output[0][0] > backward[0][0] + DOUBLE_ERROR)) {
      cout << "\nERROR: " << backward_output[0][0] << " | " << backward[0][0] << endl;
    } */
#   endif

    if ((model_type[0] != MEAN_CHANGE) && (model_type[0] != INTERCEPT_SLOPE_CHANGE)) {
      segmentation_likelihood = forward[seq_length - 1][nb_segment - 1];

      if (likelihood != D_INF) {
        for (i = 0;i < seq_length;i++) {
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
        for (i = 0;i < seq_length;i++) {
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
      if (forward[seq_length - 1][nb_segment - 1] < 0.) {
        count = (index == I_DEFAULT ? nb_sequence : 1);

        segmentation_likelihood = -((double)(count * seq_length) / 2.) *
                                   (log(-forward[seq_length - 1][nb_segment - 1] /
                                     (count * seq_length)) + log(2 * M_PI) + 1);
/*        segmentation_likelihood = -((double)(count * seq_length) / 2.) *
                                   (log(-forward[seq_length - 1][nb_segment - 1] /
                                     (count * (seq_length - nb_segment))) + log(2 * M_PI)) -
                                   (double)(count * (seq_length - nb_segment)) / 2.; */

        for (i = 0;i < seq_length;i++) {
          for (j = 0;j < nb_segment;j++) {
            if (backward_output[i][j] < 0.) {
              backward_output[i][j] = pow(backward_output[i][j] / forward[seq_length - 1][nb_segment - 1] ,
                                          -((double)(count * seq_length) / 2.));
/*              backward_output[i][j] = exp(-((double)(count * seq_length) / 2.) *
                                          log(backward_output[i][j] / forward[seq_length - 1][nb_segment - 1])); */
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
      for (i = 0;i < seq_length;i++) {
        for (j = 0;j < nb_segment;j++) {
          backward_output[i][j] = 0.;
        }
      }
    }

    switch (format) {

    case ASCII : {
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

      if ((index != I_DEFAULT) || (common_contrast)) {
        output_piecewise_function = new double*[nb_variable];

        for (i = 1;i < nb_variable;i++) {
          if (piecewise_function[i]) {
            output_piecewise_function[i] = piecewise_function[i][index == I_DEFAULT ? 0 : index];
          }
          else {
            output_piecewise_function[i] = NULL;
          }
        }
      }

      else {
        output_piecewise_function = NULL;
      }

      profile_ascii_print(*os , index , nb_segment , backward_output ,
                          (output == CHANGE_POINT ? SEQ_label[SEQL_CHANGE_POINT] : SEQ_label[SEQL_SEGMENT]) ,
                          output_piecewise_function);
      delete [] output_piecewise_function;

      *os << "\n" << SEQ_label[SEQL_SEGMENTATION_LIKELIHOOD] << ": " << segmentation_likelihood;
      if (likelihood != D_INF) {
        *os << "   (" << exp(segmentation_likelihood - likelihood) << ")";
      }
      *os << endl;
      break;
    }

    case SPREADSHEET : {
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
                                common_contrast , piecewise_function);

      *os << "\n" << SEQ_label[SEQL_SEGMENTATION_LIKELIHOOD] << "\t" << segmentation_likelihood;
      if (likelihood != D_INF) {
        *os << "\t" << exp(segmentation_likelihood - likelihood);
      }
      *os << endl;
      break;
    }

    case GNUPLOT : {
      profile_plot_print(*os , index , nb_segment , backward_output , common_contrast , piecewise_function);
      break;
    }

    case PLOT : {
      MultiPlotSet &plot = *plot_set;

      i = 0;
      for (j = 1;j < nb_variable;j++) {
        if ((piecewise_function) && (piecewise_function[j])) {
          plot[i].resize(2);

          if ((index != I_DEFAULT) || (!common_contrast)) {
            if (type[j] != REAL_VALUE) {
              for (k = 0;k < nb_sequence;k++) {
                if ((index == I_DEFAULT) || (index == k)) {
                  for (m = 0;m < length[k];m++) {
                    plot[i][0].add_point(seq_index_parameter[m] , int_sequence[k][j][m]);
                    plot[i][1].add_point(seq_index_parameter[m] , piecewise_function[j][k][m]);
                  }
                }
              }
            }
            else {
              for (k = 0;k < nb_sequence;k++) {
                if ((index == I_DEFAULT) || (index == k)) {
                  for (m = 0;m < length[k];m++) {
                    plot[i][0].add_point(seq_index_parameter[m] , real_sequence[k][j][m]);
                    plot[i][1].add_point(seq_index_parameter[m] , piecewise_function[j][k][m]);
                  }
                }
              }
            }
          }

          else {
            if (type[j] != REAL_VALUE) {
              for (k = 0;k < nb_sequence;k++) {
                for (m = 0;m < length[k];m++) {
                  plot[i][0].add_point(seq_index_parameter[m] , int_sequence[k][j][m]);
                }
              }
            }
            else {
              for (k = 0;k < nb_sequence;k++) {
                for (m = 0;m < length[k];m++) {
                  plot[i][0].add_point(seq_index_parameter[m] , real_sequence[k][j][m]);
                }
              }
            }
            for (k = 0;k < length[0];k++) {
              plot[i][1].add_point(seq_index_parameter[k] , piecewise_function[j][0][k]);
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
    if (format != GNUPLOT) {
      double ambiguity = 0.;

      psegment = int_sequence[index == I_DEFAULT ? 0 : index][0];
      for (i = 0;i < seq_length;i++) {
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
      case ASCII :
        *os << "\n" << SEQ_label[SEQL_AMBIGUITY] << ": " << ambiguity
            << " (" << ambiguity / seq_length << ")" << endl;
        break;
      case SPREADSHEET :
        *os << "\n" << SEQ_label[SEQL_AMBIGUITY] << "\t" << ambiguity
            << "\t" << ambiguity / seq_length << "\t" << endl;
        break;
      }
    }
#   endif

    delete [] change_point;
  }

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

  for (i = 0;i < seq_length;i++) {
    delete [] backward[i];
  }
  delete [] backward;

  for (i = 0;i < seq_length;i++) {
    delete [] backward_output[i];
  }
  delete [] backward_output;

  for (i = 1;i < nb_variable;i++) {
    if ((model_type[i - 1] == POISSON_CHANGE) || (model_type[i - 1] == NEGATIVE_BINOMIAL_0_CHANGE) ||
        (model_type[i - 1] == NEGATIVE_BINOMIAL_1_CHANGE) || (model_type[i - 1] == GAUSSIAN_CHANGE) ||
        (model_type[0] == MEAN_CHANGE) || (model_type[i - 1] == VARIANCE_CHANGE) ||
        (model_type[i - 1] == LINEAR_MODEL_CHANGE) || (model_type[0] == INTERCEPT_SLOPE_CHANGE) ||
        (model_type[i - 1] == AUTOREGRESSIVE_MODEL_CHANGE) || (model_type[i - 1] == STATIONARY_AUTOREGRESSIVE_MODEL_CHANGE) ||
        (model_type[i - 1] == BAYESIAN_POISSON_CHANGE) || (model_type[i - 1] == BAYESIAN_GAUSSIAN_CHANGE)) {
      for (j = 0;j < nb_sequence;j++) {
        delete [] piecewise_function[i][j];
      }
      delete [] piecewise_function[i];
    }
  }
  delete [] piecewise_function;

  return segmentation_likelihood;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the N most probable segmentations, of segment/change-point profiles and
 *         entropy profiles for a single sequence or a sample of sequences.
 *
 *  \param[in] error           reference on a StatError object,
 *  \param[in] os              stream,
 *  \param[in] iidentifier     sequence identifier,
 *  \param[in] nb_segment      number of segments,
 *  \param[in] model_type      segment model types,
 *  \param[in] common_contrast flag contrast functions common to the individuals,
 *  \param[in] shape_parameter negative binomial shape parameters,
 *  \param[in] output          output type,
 *  \param[in] format          output format (ASCII/SPREADSHEET),
 *  \param[in] segmentation    method for computing segmentations (FORWARD_DYNAMIC_PROGRAMMING/ FORWARD_BACKWARD_SAMPLING),
 *  \param[in] nb_segmentation number of segmentations.
 *
 *  \return                    error status.
 */
/*--------------------------------------------------------------*/

bool Sequences::segment_profile_write(StatError &error , ostream &os , int iidentifier ,
                                      int nb_segment , segment_model *model_type ,
                                      bool common_contrast , double *shape_parameter ,
                                      change_point_profile output , output_format format ,
                                      latent_structure_algorithm segmentation , int nb_segmentation) const

{
  bool status = true;
  int i , j;
  int index;
  double segment_length_max , likelihood = D_INF , segmentation_likelihood , **rank;
  Sequences *seq;


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
    if ((nb_segment < 2) || (nb_segment > length[index == I_DEFAULT ? 0 : index] / 2)) {
      status = false;
      error.update(SEQ_error[SEQR_NB_SEGMENT]);
    }

    if (nb_segmentation < 2) {
      status = false;
      error.update(SEQ_error[SEQR_NB_SEGMENTATION]);
    }
  }

  if (status) {
    seq = new Sequences(*this , ADD_STATE_VARIABLE);

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

    if ((model_type[0] != MEAN_CHANGE) && (model_type[0] != INTERCEPT_SLOPE_CHANGE)) {
      likelihood = seq->forward_backward(index , nb_segment , model_type , common_contrast ,
                                         shape_parameter , rank , &os , NULL ,
                                         segment_length_max , output , format);
    }
    segmentation_likelihood = seq->forward_backward_dynamic_programming(index , nb_segment , model_type ,
                                                                        common_contrast , shape_parameter ,
                                                                        rank , &os , NULL , output , format ,
                                                                        likelihood);
    if (segmentation_likelihood == D_INF) {
      status = false;
      error.update(SEQ_error[SEQR_SEGMENTATION_FAILURE]);
    }

    else if ((format == ASCII) || (length[index == I_DEFAULT ? 0 : index] <= 400)) {
      switch (segmentation) {
      case FORWARD_DYNAMIC_PROGRAMMING :
        seq->N_segmentation(index , nb_segment , model_type , common_contrast , shape_parameter ,
                            rank , os , format , nb_segmentation , likelihood);
        break;
      case FORWARD_BACKWARD_SAMPLING :
        seq->forward_backward_sampling(index , nb_segment , model_type , common_contrast ,
                                       shape_parameter , rank , os , format , nb_segmentation);
        break;
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


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the N most probable segmentations, of segment/change-point profiles and
 *         entropy profiles for a single sequence or a sample of sequences and displaying the results.
 *
 *  \param[in] error           reference on a StatError object,
 *  \param[in] iidentifier     sequence identifier,
 *  \param[in] nb_segment      number of segments,
 *  \param[in] model_type      segment model types,
 *  \param[in] common_contrast flag contrast functions common to the individuals,
 *  \param[in] shape_parameter negative binomial shape parameters,
 *  \param[in] output          output type,
 *  \param[in] segmentation    method for computing segmentations (FORWARD_DYNAMIC_PROGRAMMING/ FORWARD_BACKWARD_SAMPLING),
 *  \param[in] nb_segmentation number of segmentations.
 *
 *  \return                    error status.
 */
/*--------------------------------------------------------------*/

bool Sequences::segment_profile_ascii_write(StatError &error , int iidentifier ,
                                            int nb_segment , vector<segment_model> model_type ,
                                            bool common_contrast , vector<double> shape_parameter ,
                                            change_point_profile output ,
                                            latent_structure_algorithm segmentation , int nb_segmentation) const

{
  return segment_profile_write(error , cout , iidentifier , nb_segment , model_type.data() ,
                               common_contrast , shape_parameter.data() , output ,
                               ASCII , segmentation , nb_segmentation);
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the N most probable segmentations, of segment/change-point profiles and
 *         entropy profiles for a single sequence or a sample of sequences and
 *         writing of the results in a file.
 *
 *  \param[in] error           reference on a StatError object,
 *  \param[in] path            file path,
 *  \param[in] iidentifier     sequence identifier,
 *  \param[in] nb_segment      number of segments,
 *  \param[in] model_type      segment model types,
 *  \param[in] common_contrast flag contrast functions common to the individuals,
 *  \param[in] shape_parameter negative binomial shape parameters,
 *  \param[in] output          output type,
 *  \param[in] format          file format (ASCII/SPREADSHEET),
 *  \param[in] segmentation    method for computing segmentations (FORWARD_DYNAMIC_PROGRAMMING/FORWARD_BACKWARD_SAMPLING),
 *  \param[in] nb_segmentation number of segmentations.
 *
 *  \return                    error status.
 */
/*--------------------------------------------------------------*/

bool Sequences::segment_profile_write(StatError &error , const string path , int iidentifier ,
                                      int nb_segment , vector<segment_model> model_type ,
                                      bool common_contrast , vector<double> shape_parameter ,
                                      change_point_profile output , output_format format ,
                                      latent_structure_algorithm segmentation , int nb_segmentation) const

{
  bool status = true;
  ofstream out_file(path.c_str());


  error.init();

  if (!out_file) {
    status = false;
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {
    status = segment_profile_write(error , out_file , iidentifier , nb_segment , model_type.data() ,
                                   common_contrast , shape_parameter.data() , output , format ,
                                   segmentation , nb_segmentation);
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of segment/change-point profiles and of entropy profiles for
 *         a single sequence or a sample of sequences and plot of the profiles 
 *         at the Gnuplot format.
 *
 *  \param[in] error           reference on a StatError object,
 *  \param[in] prefix          file prefix,
 *  \param[in] iidentifier     sequence identifier,
 *  \param[in] nb_segment      number of segments,
 *  \param[in] model_type      segment model types,
 *  \param[in] common_contrast flag contrast functions common to the individuals,
 *  \param[in] shape_parameter negative binomial shape parameters,
 *  \param[in] output          output type,
 *  \param[in] title           figure title.
 *
 *  \return                    error status.
 */
/*--------------------------------------------------------------*/

bool Sequences::segment_profile_plot_write(StatError &error , const char *prefix , int iidentifier ,
                                           int nb_segment , segment_model *model_type ,
                                           bool common_contrast , double *shape_parameter ,
                                           change_point_profile output , const char *title) const

{
  bool status = true;
  int i , j , k , m;
  int index , seq_length , *seq_index_parameter;
  double segment_length_max , likelihood = D_INF , segmentation_likelihood , **rank;
  Sequences *seq;
  ostringstream data_file_name[2];
  ofstream *out_data_file;


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

  if ((status) && ((nb_segment < 2) || (nb_segment > length[index == I_DEFAULT ? 0 : index] / 2))) {
    status = false;
    error.update(SEQ_error[SEQR_NB_SEGMENT]);
  }

  if (status) {

    // writing of the data files

    i = (((model_type[0] == MEAN_CHANGE) || (model_type[0] == INTERCEPT_SLOPE_CHANGE)) ? 0 : 1);
    data_file_name[i] << prefix << i << ".dat";
    out_data_file = new ofstream((data_file_name[i].str()).c_str());

    if (!out_data_file) {
      status = false;
      error.update(STAT_error[STATR_FILE_PREFIX]);
    }

    else {
      seq = new Sequences(*this , ADD_STATE_VARIABLE);

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

      if ((model_type[0] != MEAN_CHANGE) && (model_type[0] != INTERCEPT_SLOPE_CHANGE)) {
        likelihood = seq->forward_backward(index , nb_segment , model_type , common_contrast ,
                                           shape_parameter , rank , out_data_file , NULL ,
                                           segment_length_max , output , GNUPLOT);
        out_data_file->close();
        delete out_data_file;

        data_file_name[0] << prefix << 0 << ".dat";
        out_data_file = new ofstream((data_file_name[0].str()).c_str());
      }

#     ifdef DEBUG
      likelihood = D_INF;
#     endif

      segmentation_likelihood = seq->forward_backward_dynamic_programming(index , nb_segment , model_type ,
                                                                          common_contrast , shape_parameter ,
                                                                          rank , out_data_file , NULL ,
                                                                          output , GNUPLOT , likelihood);
      out_data_file->close();
      delete out_data_file;

      if (segmentation_likelihood == D_INF) {
        status = false;
        error.update(SEQ_error[SEQR_SEGMENTATION_FAILURE]);
      }

      else {
        seq_length = seq->length[index == I_DEFAULT ? 0 : index];

        if (index_param_type == IMPLICIT_TYPE) {
          seq_index_parameter = new int[seq_length];
          for (j = 0;j < seq_length;j++) {
            seq_index_parameter[j] = j;
          }
        }
        else {
          seq_index_parameter = seq->index_parameter[index == I_DEFAULT ? 0 : index];
        }

        // writing of the script files

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

//          if (index_parameter) {
          if (seq_index_parameter[seq_length - 1] - seq_index_parameter[0] < TIC_THRESHOLD) {
            out_file << "set xtics 0,1" << endl;
          }

          j = 2;
          for (k = 1;k < seq->nb_variable;k++) {
            if ((model_type[k - 1] == POISSON_CHANGE) || (model_type[k - 1] == NEGATIVE_BINOMIAL_0_CHANGE) ||
                (model_type[k - 1] == NEGATIVE_BINOMIAL_1_CHANGE) || (model_type[k - 1] == GAUSSIAN_CHANGE) ||
                (model_type[0] == MEAN_CHANGE) || (model_type[k - 1] == VARIANCE_CHANGE) ||
                (model_type[k - 1] == LINEAR_MODEL_CHANGE) || (model_type[0] == INTERCEPT_SLOPE_CHANGE) ||
                (model_type[k - 1] == AUTOREGRESSIVE_MODEL_CHANGE) || (model_type[k - 1] == STATIONARY_AUTOREGRESSIVE_MODEL_CHANGE) ||
                (model_type[k - 1] == BAYESIAN_POISSON_CHANGE) || (model_type[k - 1] == BAYESIAN_GAUSSIAN_CHANGE)) {
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

              out_file << "plot [" << seq_index_parameter[0] << ":"
                       << seq_index_parameter[seq_length - 1] << "] ["
                       << MIN(seq->min_value[k] , 0) << ":"
                       << MAX(seq->max_value[k] , seq->min_value[k] + 1) << "] ";

              for (m = 0;m < nb_sequence;m++) {
                if ((index == I_DEFAULT) || (index == m)) {
                  out_file << "\"" << label((data_file_name[0].str()).c_str()) << "\" using " << 1 << " : " << j++;

                  if (index == m) {
                    out_file << " title \"" << SEQ_label[SEQL_SEQUENCE] << " " << iidentifier << "\"";
                    if ((index_interval) && (index_interval->variance > index_interval->mean)) {
                      out_file << " with points,\\" << endl;
                    }
                    else {
                      out_file << " with linespoints,\\" << endl; 
                    }
                    out_file << "\"" << label((data_file_name[0].str()).c_str()) << "\" using " << 1 << " : " << j++
                             << " title \"" << SEQ_label[SEQL_PIECEWISE_LINEAR_FUNCTION] << "\" with lines" << endl;
                  }

                  else {
                    out_file << " notitle with points,\\" << endl;
                    if (!common_contrast) {
                      out_file << "\"" << label((data_file_name[0].str()).c_str()) << "\" using " << 1 << " : " << j++
                               << " notitle with lines";
                      if (m < nb_sequence - 1) {
                        out_file << ",\\";
                      }
                      out_file << endl;
                    }
                  }
                }
              }

              if ((index == I_DEFAULT) && (common_contrast)) {
                out_file << "\"" << label((data_file_name[0].str()).c_str()) << "\" using " << 1 << " : " << j++
                         << " notitle with lines" << endl;
              }

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

          out_file << "plot [" << seq_index_parameter[0] << ":"
                   << seq_index_parameter[seq_length - 1];
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

            out_file << "plot [" << seq_index_parameter[0] << ":"
                     << seq_index_parameter[seq_length - 1] << "] [0:1] ";
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

            out_file << "plot [" << seq_index_parameter[0] << ":"
                     << seq_index_parameter[seq_length - 1] << "] [0:1] ";
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
            out_file << SEQ_label[SEQL_SEGMENT_LENGTH] << " " << STAT_label[STATL_DISTRIBUTIONS] << "\"\n\n";

            out_file << "plot [" << 0 << ":" << seq_length - nb_segment + 1 << "] [0:" << segment_length_max << "] ";
            for (k = 0;k < nb_segment;k++) {
              out_file << "\"" << label((data_file_name[1].str()).c_str()) << "\" using "
                       << 2 * nb_segment + k + 1 << " title \""
                       << SEQ_label[SEQL_SEGMENT] << " " << k << "\" with linespoints" << ",\\" << endl;
            }
            out_file << "\"" << label((data_file_name[1].str()).c_str()) << "\" using "
                     << 3 * nb_segment + 1 << " title \""
                     << SEQ_label[SEQL_PRIOR_SEGMENT_LENGTH] << "\" with linespoints" << endl;

            if (i == 0) {
              out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
            }
            out_file << endl;

            out_file << "set title \"";
            if (title) {
              out_file << title << " - ";
            }
            out_file << SEQ_label[SEQL_BEGIN_CONDITIONAL_ENTROPY] << "\"\n\n";

            out_file << "plot [" << seq_index_parameter[0] << ":"
                     << seq_index_parameter[seq_length - 1]
                     << "] [0:" << log(2.) << "] ";
            for (k = MAX(1 , nb_segment - 3);k < nb_segment;k++) {
              out_file << "\"" << label((data_file_name[1].str()).c_str()) << "\" using "
                       << 1 << " : " << 3 * nb_segment + k + 1 << " title \"" << k + 1 << " "
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

            out_file << "plot [" << seq_index_parameter[0] << ":"
                     << seq_index_parameter[seq_length - 1]
                     << "] [0:" << log(2.) << "] ";
            for (k = MAX(1 , nb_segment - 3);k < nb_segment;k++) {
              out_file << "\"" << label((data_file_name[1].str()).c_str()) << "\" using "
                       << 1 << " : " << 4 * nb_segment + k << " title \"" << k + 1 << " "
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

            out_file << "plot [" << seq_index_parameter[0] << ":"
                     << seq_index_parameter[seq_length - 1]
                     << "] [0:" << log(2.) << "] "
                     << "\"" << label((data_file_name[1].str()).c_str()) << "\" using "
                     << 1 << " : " << 4 * nb_segment << " title \""
                     << SEQ_label[SEQL_BEGIN_CONDITIONAL_ENTROPY] << "\" with linespoints,\\" << endl;
            out_file << "\"" << label((data_file_name[1].str()).c_str()) << "\" using "
                     << 1 << " : " << 5 * nb_segment - 1 << " title \""
                     << SEQ_label[SEQL_END_CONDITIONAL_ENTROPY] << "\" with linespoints,\\" << endl;
            out_file << "\"" << label((data_file_name[1].str()).c_str()) << "\" using "
                     << 1 << " : " << 5 * nb_segment << " title \""
                     << SEQ_label[SEQL_CHANGE_POINT_ENTROPY] << "\" with linespoints" << endl;
          }

          if (seq_index_parameter[seq_length - 1] - seq_index_parameter[0] < TIC_THRESHOLD) {
            out_file << "set xtics autofreq" << endl;
          }

          out_file << "\npause 0 \"" << STAT_label[STATL_END] << "\"" << endl;
        }

        if (index_param_type == IMPLICIT_TYPE) {
          delete [] seq_index_parameter;
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


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of segment/change-point profiles and of entropy profiles for
 *         a single sequence or a sample of sequences and plot of the profiles.
 *
 *  \param[in] error           reference on a StatError object,
 *  \param[in] iidentifier     sequence identifier,
 *  \param[in] nb_segment      number of segments,
 *  \param[in] model_type      segment model types,
 *  \param[in] common_contrast flag contrast functions common to the individuals,
 *  \param[in] shape_parameter negative binomial shape parameters,
 *  \param[in] output          output type.
 *
 *  \return                    MultiPlotSet object.
 */
/*--------------------------------------------------------------*/

MultiPlotSet* Sequences::segment_profile_plotable_write(StatError &error , int iidentifier ,
                                                        int nb_segment , segment_model *model_type ,
                                                        bool common_contrast , double *shape_parameter ,
                                                        change_point_profile output) const

{
  bool status = true;
  int i , j , k;
  int index , nb_plot_set , segmentation_index , seq_length;
  double segment_length_max , likelihood = D_INF , segmentation_likelihood , **rank;
  Sequences *seq;
  ostringstream title , legend;
  MultiPlotSet *plot_set;


  plot_set = NULL;
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

  if ((status)  && ((nb_segment < 2) || (nb_segment > length[index == I_DEFAULT ? 0 : index] / 2))) {
    status = false;
    error.update(SEQ_error[SEQR_NB_SEGMENT]);
  }

  if (status) {
    seq = new Sequences(*this , ADD_STATE_VARIABLE);

    // computation of the number of plots

    nb_plot_set = 1;
    for (i = 1;i < seq->nb_variable;i++) {
      if ((model_type[i - 1] == POISSON_CHANGE) || (model_type[i - 1] == NEGATIVE_BINOMIAL_0_CHANGE) ||
          (model_type[i - 1] == NEGATIVE_BINOMIAL_1_CHANGE) || (model_type[i - 1] == GAUSSIAN_CHANGE) ||
          (model_type[0] == MEAN_CHANGE) || (model_type[i - 1] == VARIANCE_CHANGE) ||
          (model_type[i - 1] == LINEAR_MODEL_CHANGE) || (model_type[0] == INTERCEPT_SLOPE_CHANGE) ||
          (model_type[i - 1] == AUTOREGRESSIVE_MODEL_CHANGE) || (model_type[i - 1] == STATIONARY_AUTOREGRESSIVE_MODEL_CHANGE) ||
          (model_type[i - 1] == BAYESIAN_POISSON_CHANGE) || (model_type[i - 1] == BAYESIAN_GAUSSIAN_CHANGE)) {
        nb_plot_set++;
      }
    }
    if ((model_type[0] != MEAN_CHANGE) && (model_type[0] != INTERCEPT_SLOPE_CHANGE)) {
      nb_plot_set += 6;
    }

    plot_set = new MultiPlotSet(nb_plot_set);

    MultiPlotSet &plot = *plot_set;

    plot.border = "15 lw 0";

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

    if ((model_type[0] != MEAN_CHANGE) && (model_type[0] != INTERCEPT_SLOPE_CHANGE)) {
      likelihood = seq->forward_backward(index , nb_segment , model_type ,common_contrast ,
                                         shape_parameter , rank , NULL , plot_set ,
                                         segment_length_max , output , PLOT);
    }

#   ifdef DEBUG
    likelihood = D_INF;
#   endif

    segmentation_likelihood = seq->forward_backward_dynamic_programming(index , nb_segment , model_type ,
                                                                        common_contrast , shape_parameter ,
                                                                        rank , NULL , plot_set ,
                                                                        output , PLOT , likelihood);

    if (segmentation_likelihood == D_INF) {
      delete plot_set;
      plot_set = NULL;
      error.update(SEQ_error[SEQR_SEGMENTATION_FAILURE]);
    }

    else {
      segmentation_index = (index == I_DEFAULT ? 0 : index);
      seq_length = seq->length[segmentation_index];

      i = 0;

      // scatterplot of sequences and piecewise linear fonction

      for (j = 1;j < seq->nb_variable;j++) {
        if ((model_type[j - 1] == POISSON_CHANGE) || (model_type[j - 1] == NEGATIVE_BINOMIAL_0_CHANGE) ||
            (model_type[j - 1] == NEGATIVE_BINOMIAL_1_CHANGE) || (model_type[j - 1] == GAUSSIAN_CHANGE) ||
            (model_type[0] == MEAN_CHANGE) || (model_type[j - 1] == VARIANCE_CHANGE) ||
            (model_type[j - 1] == LINEAR_MODEL_CHANGE) || (model_type[0] == INTERCEPT_SLOPE_CHANGE) ||
            (model_type[j - 1] == AUTOREGRESSIVE_MODEL_CHANGE) || (model_type[j - 1] == STATIONARY_AUTOREGRESSIVE_MODEL_CHANGE) ||
            (model_type[j - 1] == BAYESIAN_POISSON_CHANGE) || (model_type[j - 1] == BAYESIAN_GAUSSIAN_CHANGE)) {
          if (seq->nb_variable > 2) {
            title.str("");
            title << STAT_label[STATL_VARIABLE] << " " << j;
            plot[i].title = title.str();
          }

          if (seq->index_parameter) {
            plot[i].xrange = Range(seq->index_parameter[index][0] , seq->index_parameter[index][seq_length - 1]);
            if (seq->index_parameter[index][seq_length - 1] - seq->index_parameter[index][0] < TIC_THRESHOLD) {
              plot[i].xtics = 1;
            }
          }

          else {
            plot[i].xrange = Range(0 , seq_length - 1);
            if (seq_length - 1 < TIC_THRESHOLD) {
              plot[i].xtics = 1;
            }
          }

          plot[i].yrange = Range(MIN(seq->min_value[j] , 0) , MAX(seq->max_value[j] , seq->min_value[j] + 1));

          legend.str("");
          legend << SEQ_label[SEQL_SEQUENCE] << " " << iidentifier;
          plot[i][0].legend = legend.str();

          if ((index == I_DEFAULT) || ((index_interval) && (index_interval->variance > index_interval->mean))) {
            plot[i][0].style = "points";
          }
          else {
            plot[i][0].style = "linespoints";
          }

          plot[i][1].legend = SEQ_label[SEQL_PIECEWISE_LINEAR_FUNCTION];
          plot[i][1].style = "lines";
          i++;
        }
      }

      // maximum posterior probabilities

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
        plot[i].xrange = Range(seq->index_parameter[index][0] , seq->index_parameter[index][seq_length - 1]);
        if (seq->index_parameter[index][seq_length - 1] - seq->index_parameter[index][0] < TIC_THRESHOLD) {
          plot[i].xtics = 1;
        }
      }

      else {
        plot[i].xrange = Range(0 , seq_length - 1);
        if (seq_length - 1 < TIC_THRESHOLD) {
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

      // smoothed probabilities

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
          plot[i].xrange = Range(seq->index_parameter[index][0] , seq->index_parameter[index][seq_length - 1]);
          if (seq->index_parameter[index][seq_length - 1] - seq->index_parameter[index][0] < TIC_THRESHOLD) {
            plot[i].xtics = 1;
          }
        }

        else {
          plot[i].xrange = Range(0 , seq_length - 1);
          if (seq_length - 1 < TIC_THRESHOLD) {
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

        // change-point profiles

        plot[i].title = SEQ_label[SEQL_POSTERIOR_CHANGE_POINT_PROBABILITY];

        if (seq->index_parameter) {
          plot[i].xrange = Range(seq->index_parameter[index][0] , seq->index_parameter[index][seq_length - 1]);
          if (seq->index_parameter[index][seq_length - 1] - seq->index_parameter[index][0] < TIC_THRESHOLD) {
            plot[i].xtics = 1;
          }
        }

        else {
          plot[i].xrange = Range(0 , seq_length - 1);
          if (seq_length - 1 < TIC_THRESHOLD) {
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

        // segment length distributions

        title.str("");
        title << SEQ_label[SEQL_SEGMENT_LENGTH] << " " << STAT_label[STATL_DISTRIBUTIONS];
        plot[i].title = title.str();

        for (j = 0;j < nb_segment;j++) {
          legend  << SEQ_label[SEQL_SEGMENT] << " " << i;
          plot[i][j].legend = legend.str();

          plot[i][j].style = "linespoints";
        }

        legend  << SEQ_label[SEQL_PRIOR_SEGMENT_LENGTH];
        plot[i][nb_segment].legend = legend.str();

        plot[i][nb_segment].style = "linespoints";
        i++;

        // profiles of entropies conditional on the past

        plot[i].title = SEQ_label[SEQL_BEGIN_CONDITIONAL_ENTROPY];

        if (seq->index_parameter) {
          plot[i].xrange = Range(seq->index_parameter[index][0] , seq->index_parameter[index][seq_length - 1]);
          if (seq->index_parameter[index][seq_length - 1] - seq->index_parameter[index][0] < TIC_THRESHOLD) {
            plot[i].xtics = 1;
          }
        }

        else {
          plot[i].xrange = Range(0 , seq_length - 1);
          if (seq_length - 1 < TIC_THRESHOLD) {
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

        // profiles of entropies conditional on the future

        plot[i].title = SEQ_label[SEQL_END_CONDITIONAL_ENTROPY];

        if (seq->index_parameter) {
          plot[i].xrange = Range(seq->index_parameter[index][0] , seq->index_parameter[index][seq_length - 1]);
          if (seq->index_parameter[index][seq_length - 1] - seq->index_parameter[index][0] < TIC_THRESHOLD) {
            plot[i].xtics = 1;
          }
        }

        else {
          plot[i].xrange = Range(0 , seq_length - 1);
          if (seq_length - 1 < TIC_THRESHOLD) {
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

        // entropy profiles

        if (seq->index_parameter) {
          plot[i].xrange = Range(seq->index_parameter[index][0] , seq->index_parameter[index][seq_length - 1]);
          if (seq->index_parameter[index][seq_length - 1] - seq->index_parameter[index][0] < TIC_THRESHOLD) {
            plot[i].xtics = 1;
          }
        }

        else {
          plot[i].xrange = Range(0 , seq_length - 1);
          if (seq_length - 1 < TIC_THRESHOLD) {
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


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of segment/change-point profiles and of entropy profiles for
 *         a single sequence or a sample of sequences and plot of the profiles.
 *
 *  \param[in] error           reference on a StatError object,
 *  \param[in] iidentifier     sequence identifier,
 *  \param[in] nb_segment      number of segments,
 *  \param[in] model_type      segment model types,
 *  \param[in] common_contrast flag contrast functions common to the individuals,
 *  \param[in] shape_parameter negative binomial shape parameters,
 *  \param[in] output          output type.
 *
 *  \return                    MultiPlotSet object.
 */
/*--------------------------------------------------------------*/

MultiPlotSet* Sequences::segment_profile_plotable_write(StatError &error , int iidentifier ,
                                                        int nb_segment , vector<segment_model> model_type ,
                                                        bool common_contrast , vector<double> shape_parameter ,
                                                        change_point_profile output) const

{
  return segment_profile_plotable_write(error , iidentifier , nb_segment , model_type.data() ,
                                        common_contrast , shape_parameter.data() , output);
}


};  // namespace sequence_analysis
