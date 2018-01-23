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
 *       $Id$
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

#include "sequences.h"
#include "sequence_label.h"

using namespace std;
using namespace stat_tool;


namespace sequence_analysis {


extern double log_factorial(int value);
extern double log_binomial_coefficient(int inf_bound , double parameter , int value);
extern int column_width(int nb_value , const long double *value);


#if defined (SYSTEM_IS__CYGWIN)
#define expl exp
#endif



/*--------------------------------------------------------------*/
/**
 *  \brief Computation of segmentation and change-point entropies.
 *
 *  \param[in] index                sequence index,
 *  \param[in] nb_segment           number of segments,
 *  \param[in] model_type           segment model types,
 *  \param[in] common_contrast      flag contrast functions common to the individuals,
 *  \param[in] shape_parameter      negative binomial shape parameters,
 *  \param[in] rank                 ranks (for ordinal variables),
 *  \param[in] likelihood           pointer on the log-likelihoods of all the possibles segmentations,
 *  \param[in] segmentation_entropy pointer on the segmentation entropies,
 *  \param[in] first_order_entropy  pointer on the entropies assuming first-order dependencies,
 *  \param[in] change_point_entropy pointer on the change-point entropies considering or not change-point ranks,
 *  \param[in] uniform_entropy      pointer on the entropies corresponding to a uniform distribution assumption,
 *  \param[in] marginal_entropy     pointer on the marginal entropies.
 *
 *  \return                         log-likelihood of the multiple change-point model.
 */
/*--------------------------------------------------------------*/

double Sequences::forward_backward(int index , int nb_segment , segment_model *model_type ,
                                   bool common_contrast , double *shape_parameter ,
                                   double **rank , double *likelihood ,
                                   long double *segmentation_entropy ,
                                   long double *first_order_entropy ,
                                   long double *change_point_entropy , double *uniform_entropy ,
                                   long double *marginal_entropy) const

{
  int i , j , k , m;
  int seq_length , *inf_bound_parameter , *seq_index_parameter;
  double sum , buff , rlikelihood , **seq_mean , **hyperparam , **nb_segmentation_forward ,
         **nb_segmentation_backward , **smoothed , ***factorial , ***binomial_coeff;
  long double segment_norm , sequence_norm , lbuff , *contrast , *normalized_contrast , *norm ,
              *backward_norm , *entropy_smoothed , *segment_predicted , **forward ,
              *forward_norm , **backward , **change_point , **forward_predicted_entropy ,
              **backward_predicted_entropy;


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

    if ((model_type[i - 1] == LINEAR_MODEL_CHANGE) && (!seq_index_parameter)) {
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

  seq_length = length[index == I_DEFAULT ? 0 : index];
  contrast = new long double[seq_length];
  normalized_contrast = new long double[seq_length];

  nb_segmentation_forward = new double*[seq_length];
  for (i = 0;i < seq_length;i++) {
    nb_segmentation_forward[i] = new double[nb_segment];
  }

  forward = new long double*[seq_length];
  for (i = 0;i < seq_length;i++) {
    forward[i] = new long double[nb_segment];
  }

  segment_predicted = new long double[seq_length];

  forward_predicted_entropy = new long double*[seq_length];
  for (i = 0;i < seq_length;i++) {
    forward_predicted_entropy[i] = new long double[nb_segment];
  }

  norm = new long double[seq_length];
  forward_norm = new long double[seq_length];

  nb_segmentation_backward = new double*[seq_length];
  for (i = 0;i < seq_length;i++) {
    nb_segmentation_backward[i] = new double[nb_segment];
  }

  backward = new long double*[seq_length];
  for (i = 0;i < seq_length;i++) {
    backward[i] = new long double[nb_segment];
  }

  backward_predicted_entropy = new long double*[seq_length];
  for (i = 0;i < seq_length;i++) {
    backward_predicted_entropy[i] = new long double[nb_segment];
  }

  backward_norm = new long double[seq_length];

  smoothed = new double*[seq_length];
  for (i = 0;i < seq_length;i++) {
    smoothed[i] = new double[nb_segment];
  }

  change_point = new long double*[nb_segment];
  for (i = 1;i < nb_segment;i++) {
    change_point[i] = new long double[seq_length];
  }

  entropy_smoothed = new long double[nb_segment];

  // forward recurrence

  for (i = 0;i < seq_length;i++) {

    // computation of segment log-likelihoods

    forward_contrast(i , index , model_type , common_contrast , factorial ,
                     shape_parameter , binomial_coeff , seq_mean , seq_index_parameter ,
                     hyperparam , rank , contrast);

    // computation of the number of segmentations

    for (j = 0;j < nb_segment;j++) {
      nb_segmentation_forward[i][j] = 0;
    }

    for (j = 0;j < MIN((i < seq_length - 1 ? nb_segment - 1 : nb_segment) , i + 1);j++) {
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

    // recurrence and computation of predicted entropies

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

    for (j = 0;j < MIN((i < seq_length - 1 ? nb_segment - 1 : nb_segment) , i + 1);j++) {
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
      for (j = 0;j < MIN((i < seq_length - 1 ? nb_segment - 1 : nb_segment) , i + 1);j++) {
        forward[i][j] /= norm[i];
      }

      norm[i] = logl(norm[i]);
    }

    forward_norm[i] = segment_norm + norm[i];
  }

  // computation of the entropies corresponding to a uniform distribution assumption for the possible segmentations
  // for the different numbers of segments

  if (uniform_entropy) {
    for (i = 1;i < nb_segment;i++) {
      uniform_entropy[i] = log(nb_segmentation_forward[seq_length - 1][i]);
    }
  }

# ifdef DEBUG
  cout << "\n";
  buff = 1.;
  for (i = 1;i < nb_segment;i++) {
    buff *= (double)(seq_length - i) / (double)i;
    cout << i + 1 << " " << SEQ_label[SEQL_SEGMENTS] << ": "
         << nb_segmentation_forward[seq_length - 1][i] << " (" << buff << ") | "
         << log(nb_segmentation_forward[seq_length - 1][i]) << endl;
  }
# endif

  // extraction of the log-likelihoods of the observed sequence for the different numbers of segments

  for (i = 0;i < nb_segment;i++) {
    if (forward[seq_length - 1][i] > 0.) {
      likelihood[i] = logl(forward[seq_length - 1][i]) + forward_norm[seq_length - 1];
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

    // backward recurrence

    for (i = seq_length - 1;i >= 0;i--) {

      // computation of segment log-likelihoods

      backward_contrast(i , index , model_type , common_contrast , factorial ,
                        shape_parameter , binomial_coeff , seq_mean , seq_index_parameter ,
                        hyperparam , rank , contrast);

      // computation of the number of possible segmentations

      for (j = 0;j < nb_segment;j++) {
        nb_segmentation_backward[i][j] = 0;
      }

      for (j = MAX((i == 0 ? 0 : 1) , nb_segment + i - seq_length);j < nb_segment;j++) {
        if (j < nb_segment - 1) {
          for (k = i;k <= seq_length + j - nb_segment;k++) {
            if (contrast[k] != D_INF) {
              nb_segmentation_backward[i][j] += nb_segmentation_backward[k + 1][j + 1];
            }
          }
        }

        else {
          if (contrast[seq_length - 1] != D_INF) {
            nb_segmentation_backward[i][j]++;
          }
        }
      }

      // recurrence and computation of predicted entropies

      if (contrast[i] != D_INF) {
        normalized_contrast[i] = expl(contrast[i]);
      }
      else {
        normalized_contrast[i] = 0.;
      }

      segment_norm = 0.;
      for (j = i + 1;j < seq_length;j++) {
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

      for (j = MAX((i == 0 ? 0 : 1) , nb_segment + i - seq_length);j < nb_segment;j++) {
        if (j < nb_segment - 1) {
          for (k = i;k <= seq_length + j - nb_segment;k++) {
//            backward[i][j] += normalized_contrast[k] * backward[k + 1][j + 1];
            segment_predicted[k] = normalized_contrast[k] * backward[k + 1][j + 1];
            backward[i][j] += segment_predicted[k];
          }

          if (backward[i][j] > 0.) {
            for (k = i;k <= seq_length + j - nb_segment;k++) {
              lbuff = segment_predicted[k] / backward[i][j];
              if (lbuff > 0.) {
                backward_predicted_entropy[i][j] += lbuff * (backward_predicted_entropy[k + 1][j + 1] - logl(lbuff));
              }
            }
          }
        }

        else {
          backward[i][j] = normalized_contrast[seq_length - 1];
        }

        norm[i] += backward[i][j];
      }

      if (norm[i] > 0.) {
        for (j = MAX((i == 0 ? 0 : 1) , nb_segment + i - seq_length);j < nb_segment;j++) {
          backward[i][j] /= norm[i];
        }

        norm[i] = logl(norm[i]);
      }

      backward_norm[i] = segment_norm + norm[i];

      // extraction of the smoothed probabilities

      if (i < seq_length - 1) {
        for (j = MAX(0 , nb_segment + i - seq_length);j < MIN(nb_segment , i + 1);j++) {
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

      // computation of posterior change-point probabilities for the different numbers of segments

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
          for (k = MAX(1 , j + 1 + i - seq_length);k <= MIN(j , i);k++) {
            change_point[j][i] += forward[i - 1][k - 1] * backward[i][k + nb_segment - j - 1];
          }
          change_point[j][i] *= expl(forward_norm[i - 1] + backward_norm[i] - likelihood[j]);
        }
      }

      // computation of the segmentation entropy

      segment_norm = 0.;
      for (j = i;j < seq_length;j++) {
        segment_norm += norm[j];
        if (contrast[j] != D_INF) {
          normalized_contrast[j] = expl(contrast[j] - segment_norm);
        }
        else {
          normalized_contrast[j] = 0.;
        }
      }

      // computation of the segmentation entropies for the different numbers of segments

      if (i == 0) {
        for (j = i;j < seq_length - 1;j++) {
          if (contrast[j] != D_INF) {
            lbuff = normalized_contrast[j] * contrast[j];
            for (k = 1;k < MIN(nb_segment , seq_length - j);k++) {
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
            for (k = i;k < seq_length;k++) {
              if (contrast[k] != D_INF) {
                lbuff = forward[i - 1][j - 1] * normalized_contrast[k] * contrast[k];
                for (m = j + 1;m < MIN(nb_segment , j + seq_length - k);m++) {
                  segmentation_entropy[m] -= backward[k + 1][j + nb_segment - m] *
                                             expl(forward_norm[i - 1] + backward_norm[i] - likelihood[m]) * lbuff;
//                  segmentation_entropy[m] -= forward[i - 1][j - 1] * normalized_contrast[k] * backward[k + 1][j + nb_segment - m] *
//                                             expl(forward_norm[i - 1] + backward_norm[i] - likelihood[m]) * contrast[k];
                }
              }
            }
          }

          else {
            if (contrast[seq_length - 1] != D_INF) {
              lbuff = normalized_contrast[seq_length - 1] * contrast[seq_length - 1];
              for (k = 1;k < MIN(nb_segment , i + 1);k++) {
                segmentation_entropy[k] -= forward[i - 1][k - 1] *
                                           expl(forward_norm[i - 1] + backward_norm[i] - likelihood[k]) * lbuff;
//                segmentation_entropy[k] -= forward[i - 1][k - 1] * normalized_contrast[seq_length - 1] *
//                                           expl(forward_norm[i - 1] + backward_norm[i] - likelihood[k]) * contrast[seq_length - 1];
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
           << forward_predicted_entropy[seq_length - 1][i] << ", "
           << backward_predicted_entropy[0][nb_segment - 1 - i] << ", "
           << segmentation_entropy[i] << endl;
    }
#   endif

#   ifdef MESSAGE
    for (i = 1;i < nb_segment;i++) {
      if (nb_segmentation_backward[0][nb_segment - 1 - i] != nb_segmentation_forward[seq_length - 1][i]) {
        cout << "\nERROR: " << i << "  " << nb_segmentation_forward[seq_length - 1][i]
             << " | " << nb_segmentation_backward[0][nb_segment - 1 - i] << endl;
      }
    }
#   endif

#   ifdef MESSAGE
    for (i = 0;i < seq_length - 1;i++) {
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
      for (j = 0;j < seq_length;j++) {
        sum += change_point[i][j];
      }
      if ((sum < i + 1 - DOUBLE_ERROR) || (sum > i + 1 + DOUBLE_ERROR)) {
        cout << "\nERROR: " << sum << " | " << i + 1 << endl;
      }
    }
#   endif

    // computation of the ordered change-point entropies and the marginal entropies for
    // the different numbers of segments

    for (i = 1;i < nb_segment;i++) {
      for (j = 0;j < i;j++) {
        entropy_smoothed[j] = 0.;
      }
      entropy_smoothed[i] = 1.;

      first_order_entropy[i] = 0.;
      marginal_entropy[i] = 0.;

      for (j = seq_length - 2;j >= 0;j--) {
        sequence_norm = expl(forward_norm[j] + backward_norm[j + 1] - likelihood[i]);

/*        for (k = MIN(i , j + 1) + 1;k <= i;k++) {
          entropy_smoothed[k] = 0.;
        } */

//        for (k = 0;k <= i;k++) {
        for (k = MAX(0 , i + 1 + j - seq_length);k <= MIN(i , j + 1);k++) {
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

    // computation of change-point entropies for the different numbers of segments

    for (i = 1;i < nb_segment;i++) {
      change_point_entropy[i] = 0.;
      for (j = 1;j < seq_length;j++) {
        if ((change_point[i][j] > 0.) && (change_point[i][j] < 1.)) {
          change_point_entropy[i] -= change_point[i][j] * logl(change_point[i][j]) +
                                     (1 - change_point[i][j]) * logl(1 - change_point[i][j]);
        }
      }
    }
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
  delete [] normalized_contrast;

  for (i = 0;i < seq_length;i++) {
    delete [] nb_segmentation_forward[i];
  }
  delete [] nb_segmentation_forward;

  for (i = 0;i < seq_length;i++) {
    delete [] forward[i];
  }
  delete [] forward;

  delete [] segment_predicted;

  for (i = 0;i < seq_length;i++) {
    delete [] forward_predicted_entropy[i];
  }
  delete [] forward_predicted_entropy;

  delete [] norm;
  delete [] forward_norm;

  for (i = 0;i < seq_length;i++) {
    delete [] nb_segmentation_backward[i];
  }
  delete [] nb_segmentation_backward;

  for (i = 0;i < seq_length;i++) {
    delete [] backward[i];
  }
  delete [] backward;

  for (i = 0;i < seq_length;i++) {
    delete [] backward_predicted_entropy[i];
  }
  delete [] backward_predicted_entropy;
  seq_length = length[index == I_DEFAULT ? 0 : index];

  delete [] backward_norm;

  for (i = 0;i < seq_length;i++) {
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


/*--------------------------------------------------------------*/
/**
 *  \brief Slope heuristic: date-driven slope estimation (log-likelihood as
 *         a function of the penalty shape).
 *
 *  \param[in] min_nb_segment              minimum number of segments,
 *  \param[in] max_nb_segment              maximum number of segments,
 *  \param[in] penalty_shape_type          penalty shape type,
 *  \param[in] penalty_shape               pointer on the penalty shapes,
 *  \param[in] likelihood                  pointer on the log-likelihoods,
 *  \param[in] intercept                   references on the intercept,
 *  \param[in] slope                       references on the slope,
 *  \param[in] slope_standard_deviation    slope standard deviation,
 *  \param[in] residual_standard_deviation residual standard deviation.
 */
/*--------------------------------------------------------------*/

void log_likelihood_slope(int min_nb_segment , int max_nb_segment , int penalty_shape_type ,
                          double *penalty_shape , double *likelihood ,
                          double &intercept , double &slope , double &slope_standard_deviation ,
                          double &residual_standard_deviation)

{
  int i;
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
    penalty_shape_mean = (double)(min_nb_segment + max_nb_segment - 2) / 2.;
    penalty_shape_variance = (double)(nb_model * (nb_model * nb_model - 1)) / 12.;

#   ifdef DEBUG
    double buff = 0.;
    for (i = min_nb_segment;i <= max_nb_segment;i++) {
      diff = i - 1 - penalty_shape_mean;
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


/*--------------------------------------------------------------*/
/**
 *  \brief Slope heuristic: dimension jump method.
 *
 *  \param[in] display        flag for displaying the outputs of the dimension jump method,
 *  \param[in] max_nb_segment maximum number of segments,
 *  \param[in] slope_step     step on the slope,
 *  \param[in] penalty_shape  pointer on the penalty shapes,
 *  \param[in] likelihood     pointer on the log-likelihoods.
 *
 *  \return                   optimal slope.
 */
/*--------------------------------------------------------------*/

double dimension_jump(bool display , int max_nb_segment , double slope_step ,
                      double *penalty_shape , double *likelihood)

{
  int i , j;
  int max_diff_nb_segment , nb_segment , previous_nb_segment , nb_step , *step_nb_segment;
  double slope , optimal_slope , max_likelihood , penalized_likelihood , *begin_slope , *end_slope;


  step_nb_segment = new int[max_nb_segment + 1];
  begin_slope = new double[max_nb_segment + 1];
  end_slope = new double[max_nb_segment + 1];

  max_diff_nb_segment = 0;
  slope = slope_step;
  i = 0;

  do {
    max_likelihood = D_INF;
    for (j = 1;j <= max_nb_segment;j++) {
      penalized_likelihood = 2 * (likelihood[j] - 2 * penalty_shape[j] * slope);
      if (penalized_likelihood > max_likelihood) {
        max_likelihood = penalized_likelihood;
        nb_segment = j;
      }
    }

    if (slope == slope_step) {
      step_nb_segment[i] = nb_segment;
      begin_slope[i] = slope;
    }

    else {
      if (previous_nb_segment > nb_segment) {
        end_slope[i] = slope - slope_step;
        i++;
        step_nb_segment[i] = nb_segment;
        begin_slope[i] = slope;

        if (previous_nb_segment - nb_segment > max_diff_nb_segment) {
          max_diff_nb_segment = previous_nb_segment - nb_segment;
          optimal_slope = slope;
        }
      }
    }

    previous_nb_segment = nb_segment;
    slope += slope_step;
  }
  while ((nb_segment > 1) && (slope <= MAX_SLOPE));

  end_slope[i] = ceil(slope);
  nb_step = i + 1;

  if (max_diff_nb_segment < MIN_DIMENSION_JUMP) {
    optimal_slope = D_DEFAULT;
  }
  else {
    optimal_slope *= 2;
  }

  // display of the outputs of the dimension jump method

  if (display) {
    int width[2];
    ios_base::fmtflags format_flags;

    format_flags = cout.setf(ios::left , ios::adjustfield);

    width[0] = stat_tool::column_width(max_nb_segment);
    width[1] = stat_tool::column_width(1 , begin_slope + nb_step - 1);

    cout << "\n" << SEQ_label[SEQL_PIECEWISE_STEP_FUNCTION] << endl;
    for (i = 0;i < nb_step;i++) {
      cout << setw(width[0]) << step_nb_segment[i] << "   "
           << setw(width[1]) << begin_slope[i] << " -> "
           << setw(width[1]) << end_slope[i] << endl;
    }

    cout << "\n" << SEQ_label[SEQL_DIMENSION_JUMP] << ": " << max_diff_nb_segment;
    if (optimal_slope > 0.) {
      cout << "   " << SEQ_label[SEQL_OPTIMAL_SLOPE] << ": " << optimal_slope;
    }
    cout << endl;

    cout.setf(format_flags , ios::adjustfield);
  }

  delete [] step_nb_segment;
  delete [] begin_slope;
  delete [] end_slope;

  return optimal_slope;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Optimal segmentation of a single sequence or a sample of sequences.
 *
 *  \param[in] error              reference on a StatError object,
 *  \param[in] display            flag for displaying the results of multiple change-point inference,
 *  \param[in] iidentifier        sequence identifier,
 *  \param[in] max_nb_segment     maximum number of segments,
 *  \param[in] model_type         segment model types,
 *  \param[in] common_contrast    flag contrast functions common to the individuals,
 *  \param[in] shape_parameter    negative binomial shape parameters,
 *  \param[in] criterion          criterion for the selection of the number of segments,
 *  \param[in] min_nb_segment     minimum  number of segments,
 *  \param[in] penalty_shape_type penalty shape for the slope heuristic,
 *  \param[in] output             output (sequence, entropies or Kullback-Leibler divergences).
 *
 *  \return                       Sequences object.
 */
/*--------------------------------------------------------------*/

Sequences* Sequences::segmentation(StatError &error , bool display , int iidentifier ,
                                   int max_nb_segment , segment_model *model_type ,
                                   bool common_contrast , double *shape_parameter ,
                                   model_selection_criterion criterion , int min_nb_segment ,
                                   int penalty_shape_type , sequence_type output) const

{
  bool status = true , bayesian;
  int i , j;
  int index , nb_segment , inb_sequence , *nb_parameter , ilength[4];
  variable_nature itype[1];
  double buff , segmentation_intercept , segmentation_slope , segmentation_slope_standard_deviation ,
         segmentation_residual_standard_deviation , segmentation_dimension_jump_slope , intercept ,
         slope , slope_standard_deviation , residual_standard_deviation , dimension_jump_slope ,
         scaling_factor , max_likelihood[6] , *segmentation_likelihood , *segment_penalty , *penalty_shape ,
         **penalized_likelihood , *likelihood , *uniform_entropy , *segmentation_divergence , **rank;
  long double *segmentation_entropy , *first_order_entropy , *change_point_entropy , *marginal_entropy;
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
    else if (criterion == DIMENSION_JUMP) {
      criterion = SEGMENTATION_DIMENSION_JUMP;
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
    if ((max_nb_segment < 2) || (max_nb_segment > length[index == I_DEFAULT ? 0 : index] / 2)) {
      status = false;
      error.update(SEQ_error[SEQR_MAX_NB_SEGMENT]);
    }

    if (min_nb_segment == 0) {
      min_nb_segment = max_nb_segment / 2;

      if (criterion == LIKELIHOOD_SLOPE) {
        criterion = DIMENSION_JUMP;
      }
      else if (criterion == SEGMENTATION_LIKELIHOOD_SLOPE) {
        criterion = SEGMENTATION_DIMENSION_JUMP;
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
        criterion = DIMENSION_JUMP;
      }
      else if (criterion == SEGMENTATION_LIKELIHOOD_SLOPE) {
        criterion = SEGMENTATION_DIMENSION_JUMP;
      }

      if (output == LOG_LIKELIHOOD_SLOPE) {
        output = SEQUENCE;
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

    // rank computation for ordinal variables

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

      penalized_likelihood = new double*[3];
      penalized_likelihood[2] = new double[max_nb_segment + 1];
    }

    else {
      bayesian = false;

      segment_penalty = new double[max_nb_segment + 2];

      penalized_likelihood = new double*[6];
      for (i = 0;i < 6;i++) {
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

    if (max_nb_segment >= DIMENSION_JUMP_NB_SEGMENT) {
      penalty_shape = new double[max_nb_segment + 1];

      switch (penalty_shape_type) {

      case 0 : {
        for (i = 1;i <= max_nb_segment;i++) {
          penalty_shape[i] = i - 1;
        }
        break;
      }

      case 1 : {
        buff = 1.;
        for (i = 1;i <= max_nb_segment;i++) {
          penalty_shape[i] = i - 1 + log(buff);
          buff *= (double)(seq->length[0] - i) / (double)i;
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
      }

#     ifdef MESSAGE
      double buff1;

      cout << "\nPenalty shapes" << endl;
      buff = 1.;
      buff1 = 1.;
      for (i = 1;i <= max_nb_segment;i++) {
        cout << i - 1 << "  " << i - 1 + log(buff) << "  " << log(buff1) << " | " << i - 1 + log(buff) - log(buff1) << endl;
        buff *= (double)(seq->length[0] - i) / (double)i;
        buff1 *= (double)seq->length[0] / (double)i;
      }
      cout << endl;
#     endif

    }

    seq->segmentation((index == I_DEFAULT ? index : 0) , max_nb_segment , model_type ,
                      common_contrast , shape_parameter , rank ,
                      segmentation_likelihood + 1 , nb_parameter + 1 ,
                      (bayesian ? NULL : segment_penalty + 1));

    if (max_nb_segment - min_nb_segment >= SLOPE_NB_SEGMENT_RANGE) {
      log_likelihood_slope(min_nb_segment , max_nb_segment , penalty_shape_type , penalty_shape ,
                           segmentation_likelihood , segmentation_intercept , segmentation_slope ,
                           segmentation_slope_standard_deviation , segmentation_residual_standard_deviation);
    }

    if (max_nb_segment >= DIMENSION_JUMP_NB_SEGMENT) {
      segmentation_dimension_jump_slope = dimension_jump(display , max_nb_segment , SLOPE_STEP , penalty_shape ,
                                                         segmentation_likelihood);
      if ((segmentation_dimension_jump_slope < 0.) && (criterion == SEGMENTATION_DIMENSION_JUMP)) {
        criterion = mBIC;
      }
    }

    if ((model_type[0] != MEAN_CHANGE) && (model_type[0] != INTERCEPT_SLOPE_CHANGE)) {
      seq->forward_backward((index == I_DEFAULT ? index : 0) , max_nb_segment , model_type ,
                            common_contrast , shape_parameter , rank ,
                            likelihood + 1 , segmentation_entropy + 1 ,
                            first_order_entropy + 1 , change_point_entropy + 1 ,
                            uniform_entropy + 1 , marginal_entropy + 1);

      segmentation_divergence[1] = 0;
      for (i = 2;i <= max_nb_segment;i++) {
        segmentation_divergence[i] = uniform_entropy[i] - segmentation_entropy[i];
      }

      if (max_nb_segment - min_nb_segment >= SLOPE_NB_SEGMENT_RANGE) {
        log_likelihood_slope(min_nb_segment , max_nb_segment , penalty_shape_type ,
                             penalty_shape , likelihood , intercept , slope ,
                             slope_standard_deviation , residual_standard_deviation);

#       ifdef MESSAGE
        if (penalty_shape_type != 0) {
          cout << "\nTEST, " << STAT_criterion_word[LIKELIHOOD_SLOPE] << ": " << slope << " | "
             << (likelihood[max_nb_segment] - likelihood[max_nb_segment - 1]) /
                (penalty_shape[max_nb_segment] - penalty_shape[max_nb_segment - 1]) << endl;
        }
#       endif

      }

      if (max_nb_segment >= DIMENSION_JUMP_NB_SEGMENT) {
        dimension_jump_slope = dimension_jump(display , max_nb_segment , SLOPE_STEP , penalty_shape , likelihood);
        if ((dimension_jump_slope < 0.) && (criterion == DIMENSION_JUMP)) {
          criterion = ICL;
        }
      }
    }

    if (bayesian) {
      if (likelihood[1] != D_INF) {
        penalized_likelihood[2][1] = 2 * likelihood[1];
        max_likelihood[2] = penalized_likelihood[2][1];
        nb_segment = 1;
      }

      else {
        max_nb_segment = 0;
        nb_segment = 0;
      }

      for (i = 2;i <= max_nb_segment;i++) {
        if (likelihood[i] != D_INF) {
//          penalized_likelihood[2][i] = 2 * (likelihood[i] - segmentation_entropy[i]);
          penalized_likelihood[2][i] = 2 * (likelihood[i] - segmentation_entropy[i] -
                                            uniform_entropy[i]);

          if (penalized_likelihood[2][i] > max_likelihood[2]) {
            max_likelihood[2] = penalized_likelihood[2][i];
            nb_segment = i;
          }
        }

        else {
          max_nb_segment = i - 1;
          break;
        }
      }
    }

    else {

      // computation of penalized likelihoods corresponding to the slope heuristic (data-driven slope estimation and
      // dimension jump), the ICL criterion and the modified BIC (Zhang & Siegmund, 2007)

/*      segmentation_likelihood[1] = seq->one_segment_likelihood((index == I_DEFAULT ? index : 0) , model_type ,
                                                               common_contrast , shape_parameter , rank);
      nb_parameter[1] = seq->nb_parameter_computation((index == I_DEFAULT ? index : 0) , 1 , model_type ,
                                                      common_contrast); */

      if ((model_type[0] != MEAN_CHANGE) && (model_type[0] != INTERCEPT_SLOPE_CHANGE) &&
          (likelihood[1] != D_INF)) {
        if (max_nb_segment - min_nb_segment >= SLOPE_NB_SEGMENT_RANGE) {
//          penalized_likelihood[0][1] = 2 * (likelihood[1] - 1.5 * penalty_shape[1] * slope);
          penalized_likelihood[0][1] = 2 * (likelihood[1] - 2 * penalty_shape[1] * slope);
          max_likelihood[0] = penalized_likelihood[0][1];
        }

        if ((max_nb_segment >= DIMENSION_JUMP_NB_SEGMENT) && (dimension_jump_slope > 0.)) {
          penalized_likelihood[1][1] = 2 * (likelihood[1] - 2 * penalty_shape[1] * dimension_jump_slope);
          max_likelihood[1] = penalized_likelihood[1][1];
        }

        penalized_likelihood[2][1] = 2 * likelihood[1] - nb_parameter[1] *
                                     log((double)(seq->nb_sequence * seq->length[0]));
        max_likelihood[2] = penalized_likelihood[2][1];

        nb_segment = 1;
      }

      if (segmentation_likelihood[1] != D_INF) {
        if (max_nb_segment - min_nb_segment >= SLOPE_NB_SEGMENT_RANGE) {
//          penalized_likelihood[3][1] = 2 * (segmentation_likelihood[1] - 1.5 * penalty_shape[1] * segmentation_slope);
          penalized_likelihood[3][1] = 2 * (segmentation_likelihood[1] - 2 * penalty_shape[1] * segmentation_slope);
          max_likelihood[3] = penalized_likelihood[3][1];
        }

        if ((max_nb_segment >= DIMENSION_JUMP_NB_SEGMENT) && (segmentation_dimension_jump_slope > 0.)) {
          penalized_likelihood[4][1] = 2 * (segmentation_likelihood[1] - 2 * penalty_shape[1] * segmentation_dimension_jump_slope);
          max_likelihood[4] = penalized_likelihood[4][1];
        }

        penalized_likelihood[5][1] = 2 * segmentation_likelihood[1] - nb_parameter[1] *
                                     log((double)(seq->nb_sequence * seq->length[0])) - segment_penalty[1];
        max_likelihood[5] = penalized_likelihood[5][1];

        if ((model_type[0] == MEAN_CHANGE) || (model_type[0] == INTERCEPT_SLOPE_CHANGE)) {
          nb_segment = 1;
        }
      }

      if (((model_type[0] != MEAN_CHANGE) && (model_type[0] != INTERCEPT_SLOPE_CHANGE) &&
           (likelihood[1] == D_INF)) || (segmentation_likelihood[1] == D_INF)) {
        max_nb_segment = 0;
        nb_segment = 0;
      }

/*      segmentation_likelihood[2] = seq->segmentation((index == I_DEFAULT ? index : 0) , 2 , model_type ,
                                                     common_contrast , shape_parameter , rank);
      nb_parameter[2] = seq->nb_parameter_computation((index == I_DEFAULT ? index : 0) , 2 , model_type ,
                                                      common_contrast); */

      for (i = 2;i <= max_nb_segment;i++) {
/*        segmentation_likelihood[i + 1] = seq->segmentation((index == I_DEFAULT ? index : 0) , i + 1 , model_type ,
                                                           common_contrast , shape_parameter , rank);
        nb_parameter[i + 1] = seq->nb_parameter_computation((index == I_DEFAULT ? index : 0) , i + 1 , model_type ,
                                                            common_contrast); */

        if ((model_type[0] != MEAN_CHANGE) && (model_type[0] != INTERCEPT_SLOPE_CHANGE) &&
            (likelihood[i] != D_INF)) {
          if (max_nb_segment - min_nb_segment >= SLOPE_NB_SEGMENT_RANGE) {
//            penalized_likelihood[0][i] = 2 * (likelihood[i] - 1.5 * penalty_shape[i] * slope);
            penalized_likelihood[0][i] = 2 * (likelihood[i] - 2 * penalty_shape[i] * slope);
            if (penalized_likelihood[0][i] > max_likelihood[0]) {
              max_likelihood[0] = penalized_likelihood[0][i];
              if (criterion == LIKELIHOOD_SLOPE) {
                nb_segment = i;
              }
            }
          }

          if ((max_nb_segment >= DIMENSION_JUMP_NB_SEGMENT) && (dimension_jump_slope > 0.)) {
            penalized_likelihood[1][i] = 2 * (likelihood[i] - 2 * penalty_shape[i] * dimension_jump_slope);
            if (penalized_likelihood[1][i] > max_likelihood[1]) {
              max_likelihood[1] = penalized_likelihood[1][i];
              if (criterion == DIMENSION_JUMP) {
                nb_segment = i;
              }
            }
          }

          penalized_likelihood[2][i] = 2 * (likelihood[i] - segmentation_entropy[i]) - nb_parameter[i] *
                                       log((double)(seq->nb_sequence * seq->length[0]));
          if (penalized_likelihood[2][i] > max_likelihood[2]) {
            max_likelihood[2] = penalized_likelihood[2][i];
            if (criterion == ICL) {
              nb_segment = i;
            }
          }
        }

        if (segmentation_likelihood[i] != D_INF) {
          if (max_nb_segment - min_nb_segment >= SLOPE_NB_SEGMENT_RANGE) {
//            penalized_likelihood[3][i] = 2 * (segmentation_likelihood[i] -
//                                          1.5 * penalty_shape[i] * segmentation_slope);
            penalized_likelihood[3][i] = 2 * (segmentation_likelihood[i] -
                                          2 * penalty_shape[i] * segmentation_slope);
            if (penalized_likelihood[3][i] > max_likelihood[3]) {
              max_likelihood[3] = penalized_likelihood[3][i];
              if (criterion == SEGMENTATION_LIKELIHOOD_SLOPE) {
                nb_segment = i;
              }
            }
          }

          if ((max_nb_segment >= DIMENSION_JUMP_NB_SEGMENT) && (segmentation_dimension_jump_slope > 0.)) {
            penalized_likelihood[4][i] = 2 * (segmentation_likelihood[i] -
                                          2 * penalty_shape[i] * segmentation_dimension_jump_slope);
            if (penalized_likelihood[4][i] > max_likelihood[4]) {
              max_likelihood[4] = penalized_likelihood[4][i];
              if (criterion == SEGMENTATION_DIMENSION_JUMP) {
                nb_segment = i;
              }
            }
          }

          penalized_likelihood[5][i] = 2 * segmentation_likelihood[i] - nb_parameter[i] *
                                       log((double)(seq->nb_sequence * seq->length[0])) - segment_penalty[i];
          if (penalized_likelihood[5][i] > max_likelihood[5]) {
            max_likelihood[5] = penalized_likelihood[5][i];
            if (criterion == mBIC) {
              nb_segment = i;
            }
          }
        }

        if (((model_type[0] != MEAN_CHANGE) && (model_type[0] != INTERCEPT_SLOPE_CHANGE) &&
             (likelihood[i] == D_INF)) || (segmentation_likelihood[i] == D_INF)) {
          max_nb_segment = i - 1;
          break;
        }
      }
    }

    if (nb_segment > 0) {
      if (display) {
        int width[23];
        ios_base::fmtflags format_flags;
        double norm , *posterior_probability , **weight;
        Test *test;


        format_flags = cout.setf(ios::left , ios::adjustfield);

        if ((model_type[0] != MEAN_CHANGE) && (model_type[0] != INTERCEPT_SLOPE_CHANGE)) {
          posterior_probability = new double[max_nb_segment + 1];

          likelihood[1] = segmentation_likelihood[1];
          posterior_probability[1] = 1.;
          for (i = 2;i <= max_nb_segment;i++) {
            posterior_probability[i] = exp(segmentation_likelihood[i] - likelihood[i]);
          }
        }

        if (bayesian) {
          weight = new double*[3];
          weight[2] = new double[max_nb_segment + 1];
        }
        else {
          weight = new double*[6];
          for (i = 0;i < 6;i++) {
            weight[i] = new double[max_nb_segment + 1];
          }
        }

        if ((model_type[0] != MEAN_CHANGE) && (model_type[0] != INTERCEPT_SLOPE_CHANGE)) {
          norm = 0.;
          for (i = 1;i <= max_nb_segment;i++) {
            weight[2][i] = exp((penalized_likelihood[2][i] - max_likelihood[2]) / 2);
            norm += weight[2][i];
          }
          for (i = 1;i <= max_nb_segment;i++) {
            weight[2][i] /= norm;
          }
        }

        if (!bayesian) {
          if ((model_type[0] != MEAN_CHANGE) && (model_type[0] != INTERCEPT_SLOPE_CHANGE)) {
            if (max_nb_segment - min_nb_segment >= SLOPE_NB_SEGMENT_RANGE) {
              norm = 0.;
              for (i = 1;i <= max_nb_segment;i++) {
                weight[0][i] = exp((penalized_likelihood[0][i] - max_likelihood[0]) / 2);
                norm += weight[0][i];
              }
              for (i = 1;i <= max_nb_segment;i++) {
                weight[0][i] /= norm;
              }

              test = new Test(STUDENT , false , max_nb_segment - min_nb_segment - 1 , I_DEFAULT , D_DEFAULT);
              test->critical_probability = ref_critical_probability[0];
              test->t_value_computation();

              cout << STAT_criterion_word[LIKELIHOOD_SLOPE] << ": " << slope << " ("
                   << slope - test->value * slope_standard_deviation << ", "
                   << slope + test->value * slope_standard_deviation << ") | "
                   << STAT_label[STATL_RESIDUAL] << " " << STAT_label[STATL_STANDARD_DEVIATION] << ": "
                   << residual_standard_deviation << endl;

              delete test;
            }

            if ((max_nb_segment >= DIMENSION_JUMP_NB_SEGMENT) && (dimension_jump_slope > 0.)) {
              norm = 0.;
              for (i = 1;i <= max_nb_segment;i++) {
                weight[1][i] = exp((penalized_likelihood[1][i] - max_likelihood[1]) / 2);
                norm += weight[1][i];
              }
              for (i = 1;i <= max_nb_segment;i++) {
                weight[1][i] /= norm;
              }
            }
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

            cout << STAT_criterion_word[SEGMENTATION_LIKELIHOOD_SLOPE] << ": " << segmentation_slope << " ("
                 << segmentation_slope - test->value * segmentation_slope_standard_deviation << ", "
                 << segmentation_slope + test->value * segmentation_slope_standard_deviation << ") | "
                 << STAT_label[STATL_RESIDUAL] << " " << STAT_label[STATL_STANDARD_DEVIATION] << ": "
                 << segmentation_residual_standard_deviation << endl;

            delete test;
          }

          if ((max_nb_segment >= DIMENSION_JUMP_NB_SEGMENT) && (segmentation_dimension_jump_slope > 0.)) {
            norm = 0.;
            for (i = 1;i <= max_nb_segment;i++) {
              weight[4][i] = exp((penalized_likelihood[4][i] - max_likelihood[4]) / 2);
              norm += weight[4][i];
            }
            for (i = 1;i <= max_nb_segment;i++) {
              weight[4][i] /= norm;
            }
          }

          norm = 0.;
          for (i = 1;i <= max_nb_segment;i++) {
            weight[5][i] = exp((penalized_likelihood[5][i] - max_likelihood[5]) / 2);
            norm += weight[5][i];
          }
          for (i = 1;i <= max_nb_segment;i++) {
            weight[5][i] /= norm;
          }
        }

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
//          width[9] = column_width(max_nb_segment - 1 , marginal_entropy + 2) + ASCII_SPACE
          width[10] = stat_tool::column_width(nb_parameter[max_nb_segment]) + ASCII_SPACE;

          if (!bayesian) {
            if (max_nb_segment - min_nb_segment >= SLOPE_NB_SEGMENT_RANGE) {
              width[11] = stat_tool::column_width(max_nb_segment , penalized_likelihood[0] + 1) + ASCII_SPACE;
              width[12] = stat_tool::column_width(max_nb_segment , weight[0] + 1) + ASCII_SPACE;
            }
            if ((max_nb_segment >= DIMENSION_JUMP_NB_SEGMENT) && (dimension_jump_slope > 0.)) {
              width[13] = stat_tool::column_width(max_nb_segment , penalized_likelihood[1] + 1) + ASCII_SPACE;
              width[14] = stat_tool::column_width(max_nb_segment , weight[1] + 1) + ASCII_SPACE;
            }
          }

          width[15] = stat_tool::column_width(max_nb_segment , penalized_likelihood[2] + 1) + ASCII_SPACE;
          width[16] = stat_tool::column_width(max_nb_segment , weight[2] + 1) + ASCII_SPACE;

          if (!bayesian) {
            if (max_nb_segment - min_nb_segment >= SLOPE_NB_SEGMENT_RANGE) {
              width[17] = stat_tool::column_width(max_nb_segment , penalized_likelihood[3] + 1) + ASCII_SPACE;
              width[18] = stat_tool::column_width(max_nb_segment , weight[3] + 1) + ASCII_SPACE;
            }
            if ((max_nb_segment >= DIMENSION_JUMP_NB_SEGMENT) && (segmentation_dimension_jump_slope > 0.)) {
              width[19] = stat_tool::column_width(max_nb_segment , penalized_likelihood[4] + 1) + ASCII_SPACE;
              width[20] = stat_tool::column_width(max_nb_segment , weight[4] + 1) + ASCII_SPACE;
            }
            width[21] = stat_tool::column_width(max_nb_segment , penalized_likelihood[5] + 1) + ASCII_SPACE;
            width[22] = stat_tool::column_width(max_nb_segment , weight[5] + 1) + ASCII_SPACE;
          }

          cout << "\n" << SEQ_label[SEQL_NB_SEGMENT] << " | 2 * " << STAT_label[STATL_LIKELIHOOD]
               << " | 2 * " << SEQ_label[SEQL_POSSIBLE_SEGMENTATION_LIKELIHOOD]
               << " | " << SEQ_label[SEQL_POSTERIOR_PROBABILITY]
               << " | " << SEQ_label[SEQL_SEGMENTATION_ENTROPY]
               << " | " << SEQ_label[SEQL_FIRST_ORDER_ENTROPY]
               << " | " << SEQ_label[SEQL_CHANGE_POINT_ENTROPY]
               << " | " << SEQ_label[SEQL_UNIFORM_ENTROPY]
               << " | " << SEQ_label[SEQL_SEGMENTATION_DIVERGENCE] << endl;
//               << " | " << SEQ_label[SEQL_MARGINAL_ENTROPY]

          cout << setw(width[0]) << 1
               << setw(width[1]) << 2 * segmentation_likelihood[1]
               << setw(width[2]) << 2 * likelihood[1]
               << setw(width[3]) << posterior_probability[1]
               << setw(width[4]) << " "
               << setw(width[5]) << " "
               << setw(width[6]) << " "
               << setw(width[7]) << " "
               << setw(width[8]) << segmentation_divergence[1] << endl;
//               << setw(width[9]) << " "

          for (i = 2;i <= max_nb_segment;i++) {
            cout << setw(width[0]) << i
                 << setw(width[1]) << 2 * segmentation_likelihood[i]
                 << setw(width[2]) << 2 * likelihood[i]
                 << setw(width[3]) << posterior_probability[i]
                 << setw(width[4]) << segmentation_entropy[i]
                 << setw(width[5]) << first_order_entropy[i]
                 << setw(width[6]) << change_point_entropy[i]
                 << setw(width[7]) << uniform_entropy[i]
                 << setw(width[8]) << segmentation_divergence[i] << endl;
//                 << setw(width[9]) << marginal_entropy[i]
          }

          cout << "\n" << SEQ_label[SEQL_NB_SEGMENT] << " | 2 * " << STAT_label[STATL_LIKELIHOOD]
               << " | 2 * " << SEQ_label[SEQL_POSSIBLE_SEGMENTATION_LIKELIHOOD]
               << " | " << SEQ_label[SEQL_POSTERIOR_PROBABILITY]
               << " | " << SEQ_label[SEQL_SEGMENTATION_DIVERGENCE]
               << " | " << STAT_label[STATL_FREE_PARAMETERS];

          if (!bayesian) {
            if (max_nb_segment - min_nb_segment >= SLOPE_NB_SEGMENT_RANGE) {
              cout << " | " << STAT_criterion_word[LIKELIHOOD_SLOPE] << " - "  << STAT_label[STATL_WEIGHT];
            }
            if ((max_nb_segment >= DIMENSION_JUMP_NB_SEGMENT) && (dimension_jump_slope > 0.)) {
              cout << " | " << STAT_criterion_word[DIMENSION_JUMP] << " - "  << STAT_label[STATL_WEIGHT];
            }
          }

          cout << " | "  << STAT_criterion_word[ICL] << " - "  << STAT_label[STATL_WEIGHT];

          if (!bayesian) {
            if (max_nb_segment - min_nb_segment >= SLOPE_NB_SEGMENT_RANGE) {
              cout << " | " << STAT_criterion_word[SEGMENTATION_LIKELIHOOD_SLOPE] << " - "  << STAT_label[STATL_WEIGHT];
            }
            if ((max_nb_segment >= DIMENSION_JUMP_NB_SEGMENT) && (segmentation_dimension_jump_slope > 0.)) {
              cout << " | " << STAT_criterion_word[SEGMENTATION_DIMENSION_JUMP] << " - "  << STAT_label[STATL_WEIGHT];
            }
            cout << " | " << STAT_criterion_word[mBIC] << " - "  << STAT_label[STATL_WEIGHT];
          }
          cout << endl;

          for (i = 1;i <= max_nb_segment;i++) {
            cout << setw(width[0]) << i
                 << setw(width[1]) << 2 * segmentation_likelihood[i]
                 << setw(width[2]) << 2 * likelihood[i]
                 << setw(width[3]) << posterior_probability[i]
                 << setw(width[8]) << segmentation_divergence[i]
                 << setw(width[10]) << nb_parameter[i];

            if (!bayesian) {
              if (max_nb_segment - min_nb_segment >= SLOPE_NB_SEGMENT_RANGE) {
                cout << setw(width[11]) << penalized_likelihood[0][i]
                     << setw(width[12]) << weight[0][i];
              }
              if ((max_nb_segment >= DIMENSION_JUMP_NB_SEGMENT) && (dimension_jump_slope > 0.)) {
                cout << setw(width[13]) << penalized_likelihood[1][i]
                     << setw(width[14]) << weight[1][i];
              }
            }

            cout << setw(width[15]) << penalized_likelihood[2][i]
                 << setw(width[16]) << weight[2][i];

            if (!bayesian) {
              if (max_nb_segment - min_nb_segment >= SLOPE_NB_SEGMENT_RANGE) {
                cout << setw(width[17]) << penalized_likelihood[3][i]
                     << setw(width[18]) << weight[3][i];
              }
              if ((max_nb_segment >= DIMENSION_JUMP_NB_SEGMENT) && (segmentation_dimension_jump_slope > 0.)) {
                cout << setw(width[19]) << penalized_likelihood[4][i]
                     << setw(width[20]) << weight[4][i];
              }
              cout << setw(width[21]) << penalized_likelihood[5][i]
                   << setw(width[22]) << weight[5][i];
            }
            cout << endl;
          }
        }

        else {
          width[0] = stat_tool::column_width(max_nb_segment) + ASCII_SPACE;
          width[1] = stat_tool::column_width(max_nb_segment , segmentation_likelihood + 1 , 2.) + ASCII_SPACE;
          width[10] = stat_tool::column_width(nb_parameter[max_nb_segment]) + ASCII_SPACE;

          if (max_nb_segment - min_nb_segment >= SLOPE_NB_SEGMENT_RANGE) {
            width[17] = stat_tool::column_width(max_nb_segment , penalized_likelihood[3] + 1) + ASCII_SPACE;
            width[18] = stat_tool::column_width(max_nb_segment , weight[3] + 1) + ASCII_SPACE;
          }
          if ((max_nb_segment >= DIMENSION_JUMP_NB_SEGMENT) && (segmentation_dimension_jump_slope > 0.)) {
            width[19] = stat_tool::column_width(max_nb_segment , penalized_likelihood[4] + 1) + ASCII_SPACE;
            width[20] = stat_tool::column_width(max_nb_segment , weight[4] + 1) + ASCII_SPACE;
          }
          width[21] = stat_tool::column_width(max_nb_segment , penalized_likelihood[5] + 1) + ASCII_SPACE;
          width[22] = stat_tool::column_width(max_nb_segment , weight[5] + 1) + ASCII_SPACE;

          cout << "\n" << SEQ_label[SEQL_NB_SEGMENT] << " | 2 * " << STAT_label[STATL_LIKELIHOOD]
               << " | " << STAT_label[STATL_FREE_PARAMETERS];
          if (max_nb_segment - min_nb_segment >= SLOPE_NB_SEGMENT_RANGE) {
            cout << " | " << STAT_criterion_word[SEGMENTATION_LIKELIHOOD_SLOPE] << " - "  << STAT_label[STATL_WEIGHT];
          }
          if ((max_nb_segment >= DIMENSION_JUMP_NB_SEGMENT) && (segmentation_dimension_jump_slope > 0.)) {
            cout << " | " << STAT_criterion_word[SEGMENTATION_DIMENSION_JUMP] << " - "  << STAT_label[STATL_WEIGHT];
          }
          cout << " | " << STAT_criterion_word[mBIC] << " - "  << STAT_label[STATL_WEIGHT] << endl;

          for (i = 1;i <= max_nb_segment;i++) {
            cout << setw(width[0]) << i
                 << setw(width[1]) << 2 * segmentation_likelihood[i]
                 << setw(width[10]) << nb_parameter[i];
            if (max_nb_segment - min_nb_segment >= SLOPE_NB_SEGMENT_RANGE) {
              cout << setw(width[17]) << penalized_likelihood[3][i]
                   << setw(width[18]) << weight[3][i];
            }
            if ((max_nb_segment >= DIMENSION_JUMP_NB_SEGMENT) && (segmentation_dimension_jump_slope > 0.)) {
              cout << setw(width[19]) << penalized_likelihood[4][i]
                   << setw(width[20]) << weight[4][i];
            }
            cout << setw(width[21]) << penalized_likelihood[5][i]
                 << setw(width[22]) << weight[5][i] << endl;
          }
        }

        if ((model_type[0] != MEAN_CHANGE) && (model_type[0] != INTERCEPT_SLOPE_CHANGE)) {
          delete [] posterior_probability;
        }

        if (bayesian) {
          delete [] weight[2];
        }
        else {
          for (i = 0;i < 6;i++) {
            delete [] weight[i];
          }
        }

        delete [] weight;

        cout.setf(format_flags , ios::adjustfield);
      }

      if (nb_segment == 1) {
        seq->one_segment_likelihood((index == I_DEFAULT ? index : 0) , model_type , common_contrast ,
                                    shape_parameter , rank);
        seq->min_value[0] = 0;
        seq->max_value[0] = 0;
        seq->build_marginal_frequency_distribution(0);
      }

      else {
        seq->segmentation((index == I_DEFAULT ? index : 0) , nb_segment , model_type , common_contrast ,
                          shape_parameter , rank);
      }

      switch (output) {

      case SEQUENCE : {
        oseq = seq->segmentation_output(nb_segment , model_type , common_contrast , display);
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
        cout << "\n";
        for (i = 1;i <= max_nb_segment;i++) {
          cout << i << "\t" << penalty_shape[i];
          if ((model_type[0] != MEAN_CHANGE) && (model_type[0] != INTERCEPT_SLOPE_CHANGE)) {
            cout << "\t" << likelihood[i] << "\t" << intercept + slope * penalty_shape[i];
          }
          else {
            cout << "\t" << segmentation_likelihood[i] << "\t" << segmentation_intercept + segmentation_slope * penalty_shape[i];
          }
          cout << endl;
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
          scaling_factor = PENALTY_SHAPE_SCALING_FACTOR;
          break;
        case 2 :
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
/*      if ((model_type[0] != MEAN_CHANGE) && (model_type[0] != INTERCEPT_SLOPE_CHANGE)) {
        hierarchical_segmentation(error , cout , iidentifier , max_nb_segment , model_type);
      } */
#     endif

    }

    else {
      oseq = NULL;
      error.update(SEQ_error[SEQR_SEGMENTATION_FAILURE]);
    }

    if (max_nb_segment >= DIMENSION_JUMP_NB_SEGMENT) {
      delete [] penalty_shape;
    }

    for (i = 1;i < seq->nb_variable;i++) {
      delete [] rank[i];
    }
    delete [] rank;

    delete seq;

    delete [] segmentation_likelihood;
    delete [] nb_parameter;

    if (bayesian) {
      delete [] penalized_likelihood[2];
    }

    else {
      delete [] segment_penalty;

      for (i = 0;i < 6;i++) {
        delete [] penalized_likelihood[i];
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


/*--------------------------------------------------------------*/
/**
 *  \brief Optimal segmentation of a single sequence or a sample of sequences.
 *
 *  \param[in] error              reference on a StatError object,
 *  \param[in] display            flag for displaying the results of multiple change-point inference,
 *  \param[in] iidentifier        sequence identifier,
 *  \param[in] max_nb_segment     maximum number of segments,
 *  \param[in] model_type         segment model types,
 *  \param[in] common_contrast    flag contrast functions common to the individuals,
 *  \param[in] shape_parameter    negative binomial shape parameters,
 *  \param[in] criterion          criterion for the selection of the number of segments,
 *  \param[in] min_nb_segment     minimum  number of segments,
 *  \param[in] penalty_shape_type penalty shape for the slope heuristic,
 *  \param[in] output             output (sequence, entropies or Kullback-Leibler divergences).
 *
 *  \return                       Sequences object.
 */
/*--------------------------------------------------------------*/

Sequences* Sequences::segmentation(StatError &error , bool display , int iidentifier ,
                                   int max_nb_segment , vector<segment_model> model_type ,
                                   bool common_contrast , vector<double> shape_parameter ,
                                   model_selection_criterion criterion , int min_nb_segment ,
                                   int penalty_shape_type , sequence_type output) const

{
  return segmentation(error , display , iidentifier , max_nb_segment , model_type.data() ,
                      common_contrast , shape_parameter.data() , criterion , min_nb_segment ,
                      penalty_shape_type , output);
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of segment/state, change-point and entropy profiles for
 *         a single sequence or a sample of sequences (in the case of
 *         multiple change-point models).
 *
 *  \param[in,out] os                       stream,
 *  \param[in]     index                    sequence index,
 *  \param[in]     nb_segment               number of segments/states,
 *  \param[in]     profiles                 pointer on the segment/state profiles,
 *  \param[in]     label                    profile type label,
 *  \param[in]     piecewise_function       pointer on the piecewise linear functions,
 *  \param[in]     change_point             pointer on the change-point profiles,
 *  \param[in]     segment_length           pointer on the segment length distributions,
 *  \param[in]     prior_segment_length     pointer on the prior segment length distribution assuming a
 *                                          uniform prior on the possible segmentations,
 *  \param[in]     begin_conditonal_entropy pointer on the profiles of entropies conditional on the past,
 *  \param[in]     end_conditional_entropy  pointer on the profiles of entropies conditional on the future,
 *  \param[in]     change_point_entropy     pointer on the change-point entropy profiles.
 */
/*--------------------------------------------------------------*/

ostream& Sequences::profile_ascii_print(ostream &os , int index , int nb_segment ,
                                        double **profiles , const char *label ,
                                        double **piecewise_function , long double **change_point ,
                                        Distribution **segment_length ,
                                        Distribution *prior_segment_length ,
                                        long double **begin_conditonal_entropy ,
                                        long double **end_conditional_entropy ,
                                        long double **change_point_entropy) const

{
  int i , j , k;
  int seq_length , start , buff , max_nb_value , *seq_index_parameter , *width;
  ios_base::fmtflags format_flags;


  format_flags = os.flags(ios::adjustfield);

  seq_length = length[index == I_DEFAULT ? 0 : index];

  if (index_param_type == IMPLICIT_TYPE) {
    seq_index_parameter = new int[seq_length];
    for (i = 0;i < seq_length;i++) {
      seq_index_parameter[i] = i;
    }
  }
  else {
    seq_index_parameter = index_parameter[index == I_DEFAULT ? 0 : index];
  }

  // computation of the column width

  width = new int[2 * nb_variable + 6];

  start = 0;
  if (change_point) {
    start++;
  }

  for (i = start;i < nb_variable;i++) {
    if (type[i] != REAL_VALUE) {
      width[i] = stat_tool::column_width((int)min_value[i] , (int)max_value[i]);
    }

    else {
      if ((index == I_DEFAULT) && (nb_sequence * (nb_variable - 1) <= SEQUENCE_MAX_NB_COLUMN)) {
        for (j = 0;j < nb_sequence;j++) {
          buff = stat_tool::column_width(length[j] , real_sequence[j][i]);
          if (buff > width[i]) {
            width[i] = buff;
          }
        }
      }

      else if (index != I_DEFAULT) {
        width[i] = stat_tool::column_width(length[index] , real_sequence[index][i]);
      }
    }

    if ((i > start) || (index == I_DEFAULT)) {
      width[i] += ASCII_SPACE;
    }
  }

  if (index_parameter) {
    width[nb_variable] = stat_tool::column_width(index_parameter_distribution->nb_value - 1) + ASCII_SPACE;
  }
  else {
    width[nb_variable] = stat_tool::column_width(seq_length) + ASCII_SPACE;
  }

  width[nb_variable + 1] = 0;
  for (i = 0;i < seq_length;i++) {
    buff = stat_tool::column_width(nb_segment , profiles[i]);
    if (buff > width[nb_variable + 1]) {
      width[nb_variable + 1] = buff;
    }
  }
  width[nb_variable + 1] += ASCII_SPACE;

  width[nb_variable + 2] = stat_tool::column_width(nb_sequence);

  if (piecewise_function) {
    for (i = 1;i < nb_variable;i++) {
      if (piecewise_function[i]) {
        width[nb_variable + 2 + i] = stat_tool::column_width(seq_length , piecewise_function[i]) + ASCII_SPACE;
      }
    }
  }

  if (change_point) {
    width[2 * nb_variable + 2] = 0;
    for (i = 1;i < nb_segment;i++) {
      buff = column_width(seq_length , change_point[i]);
      if (buff > width[2 * nb_variable + 2]) {
        width[2 * nb_variable + 2] = buff;
      }
    }
    width[2 * nb_variable + 2] += ASCII_SPACE;
  }

  if (segment_length) {
    width[2 * nb_variable + 3] = stat_tool::column_width(segment_length[0]->nb_value , segment_length[0]->mass);
    for (i = 1;i < nb_segment;i++) {
      buff = stat_tool::column_width(segment_length[i]->nb_value , segment_length[i]->mass);
      if (buff > width[2 * nb_variable + 3]) {
        width[2 * nb_variable + 3] = buff;
      }
    }
    if (prior_segment_length) {
      buff = stat_tool::column_width(prior_segment_length->nb_value , prior_segment_length->mass);
      if (buff > width[2 * nb_variable + 3]) {
        width[2 * nb_variable + 3] = buff;
      }
    }
    width[2 * nb_variable + 3] += ASCII_SPACE;

    width[2 * nb_variable + 4] = stat_tool::column_width(segment_length[0]->nb_value , segment_length[0]->cumul);
    for (i = 1;i < nb_segment;i++) {
      buff = stat_tool::column_width(segment_length[i]->nb_value , segment_length[i]->cumul);
      if (buff > width[2 * nb_variable + 4]) {
        width[2 * nb_variable + 4] = buff;
      }
    }
    if (prior_segment_length) {
      buff = stat_tool::column_width(prior_segment_length->nb_value , prior_segment_length->cumul);
      if (buff > width[2 * nb_variable + 4]) {
        width[2 * nb_variable + 4] = buff;
      }
    }
    width[2 * nb_variable + 4] += ASCII_SPACE;
  }

  if ((begin_conditonal_entropy) && (end_conditional_entropy) &&
      (change_point_entropy)) {
    width[2 * nb_variable + 5] = 0;

    for (i = 1;i < nb_segment;i++) {
      buff = column_width(seq_length , begin_conditonal_entropy[i]);
      if (buff > width[2 * nb_variable + 5]) {
        width[2 * nb_variable + 5] = buff;
      }
    }

    for (i = 1;i < nb_segment;i++) {
      buff = column_width(seq_length , end_conditional_entropy[i]);
      if (buff > width[2 * nb_variable + 5]) {
        width[2 * nb_variable + 5] = buff;
      }
    }

    for (i = 1;i < nb_segment;i++) {
      buff = column_width(seq_length , change_point_entropy[i]);
      if (buff > width[2 * nb_variable + 5]) {
        width[2 * nb_variable + 5] = buff;
      }
    }

    width[2 * nb_variable + 5] += ASCII_SPACE;
  }

  if (!change_point) {
    os << SEQ_label[SEQL_OPTIMAL] << " " << label << " | ";
  }
  for (i = 1;i < nb_variable;i++) {
    if ((index == I_DEFAULT) && (nb_sequence * (nb_variable - 1) <= SEQUENCE_MAX_NB_COLUMN)) {
      for (j = 0;j < nb_sequence;j++) {
        os << STAT_label[STATL_VARIABLE] << " " << i << " | ";
      }
    }
    else if (index != I_DEFAULT) {
      os << STAT_label[STATL_VARIABLE] << " " << i << " | ";
    }

    if ((piecewise_function) && (piecewise_function[i])) {
      os << SEQ_label[SEQL_PIECEWISE_LINEAR_FUNCTION] << " " << i << " | ";
    }
  }

  switch (index_param_type) {
  case TIME :
    os << SEQ_label[SEQL_TIME];
    break;
  case POSITION :
    os << SEQ_label[SEQL_POSITION];
    break;
  default :
    os << SEQ_label[SEQL_INDEX];
    break;
  }

  for (i = 0;i < nb_segment;i++) {
    os << " | " << label << " " << i;
  }
  if (change_point) {
    os << "   ";
    for (i = 1;i < nb_segment;i++) {
      os << " | " << i + 1 << " " << SEQ_label[SEQL_SEGMENTS];
    }
  }
  os << endl;

  for (i = 0;i < seq_length;i++) {
    os.setf(ios::right , ios::adjustfield);

    if (!change_point) {
      os << setw(width[0]) << int_sequence[index == I_DEFAULT ? 0 : index][0][i];
    }

    for (j = 1;j < nb_variable;j++) {
      if ((index == I_DEFAULT) && (nb_sequence * (nb_variable - 1) <= SEQUENCE_MAX_NB_COLUMN)) {
        if (type[j] != REAL_VALUE) {
          for (k = 0;k < nb_sequence;k++) {
            os << setw(width[j]) << int_sequence[k][j][i];
          }
        }
        else {
          for (k = 0;k < nb_sequence;k++) {
            os << setw(width[j]) << real_sequence[k][j][i];
          }
        }
      }

      else if (index != I_DEFAULT) {
        if (type[j] != REAL_VALUE) {
          os << setw(width[j]) << int_sequence[index][j][i];
        }
        else {
          os << setw(width[j]) << real_sequence[index][j][i];
        }
      }

      if ((piecewise_function) && (piecewise_function[j])) {
        os << setw(width[nb_variable + 2 + j]) << piecewise_function[j][i];
      }
    }

    os << setw(width[nb_variable]) << seq_index_parameter[i] << "  ";

    os.setf(ios::left , ios::adjustfield);
    for (j = 0;j < nb_segment;j++) {
      os << setw(width[nb_variable + 1]) << profiles[i][j];
    }

    if (change_point) {
      os << "   ";
      for (j = 1;j < nb_segment;j++) {
        os << setw(width[2 * nb_variable + 2]) << change_point[j][i];
      }
    }

    if (i == 0) {
      os.setf(ios::right , ios::adjustfield);
      if (index != I_DEFAULT) {
        os << setw(width[nb_variable + 2]) << identifier[index];
      }
    }
    os << endl;
  }

  if (segment_length) {
    if (prior_segment_length) {
      max_nb_value = prior_segment_length->nb_value;
    }

    else {
      max_nb_value = segment_length[0]->nb_value;
      for (i = 1;i < nb_segment;i++) {
        if (segment_length[i]->nb_value > max_nb_value) {
          max_nb_value = segment_length[i]->nb_value;
        }
      }
    }

    os << "\n\n" << SEQ_label[SEQL_SEGMENT_LENGTH] << " " << STAT_label[STATL_DISTRIBUTIONS] << endl;

    for (i = 0;i < nb_segment;i++) {
      os << "\n" << SEQ_label[SEQL_SEGMENT] << " " << i << " " << SEQ_label[SEQL_LENGTH]
         << " " << STAT_label[STATL_DISTRIBUTION] << endl;
      segment_length[i]->ascii_characteristic_print(os , true);
    }
    if (prior_segment_length) {
      os << "\n" << SEQ_label[SEQL_PRIOR_SEGMENT_LENGTH] << " " << STAT_label[STATL_DISTRIBUTION] << endl;
      prior_segment_length->ascii_characteristic_print(os , true);
    }

    os << "\n  ";
    for (i = 0;i < nb_segment;i++) {
      os << " | " << SEQ_label[SEQL_SEGMENT] << " " << i << " " << SEQ_label[SEQL_LENGTH]
         << " " << STAT_label[STATL_DISTRIBUTION];
    }
    if (prior_segment_length) {
      os << " | " << SEQ_label[SEQL_PRIOR_SEGMENT_LENGTH] << " " << STAT_label[STATL_DISTRIBUTION];
    }
    for (i = 0;i < nb_segment;i++) {
      os << " | " << SEQ_label[SEQL_SEGMENT] << " " << i << " " << STAT_label[STATL_CUMULATIVE]
         << " " << STAT_label[STATL_DISTRIBUTION] << " " << STAT_label[STATL_FUNCTION];
    }
    if (prior_segment_length) {
      os << " | " << SEQ_label[SEQL_PRIOR_SEGMENT_LENGTH] << " " << i << " " << STAT_label[STATL_CUMULATIVE]
         << " " << STAT_label[STATL_DISTRIBUTION] << " " << STAT_label[STATL_FUNCTION];
    }
    os << endl;

    for (i = 0;i < max_nb_value;i++) {
      os << setw(width[nb_variable]) << i;

      for (j = 0;j < nb_segment;j++) {
        if (i < segment_length[j]->nb_value) {
          os << setw(width[2 * nb_variable + 3]) << segment_length[j]->mass[i];
        }
        else {
          os << setw(width[2 * nb_variable + 3]) << " ";
        }
      }
      if (prior_segment_length) {
        os << setw(width[2 * nb_variable + 3]) << prior_segment_length->mass[i];
      }
      for (j = 0;j < nb_segment;j++) {
        if (i < segment_length[j]->nb_value) {
          os << setw(width[2 * nb_variable + 4]) << segment_length[j]->cumul[i];
        }
        else {
          os << setw(width[2 * nb_variable + 4]) << " ";
        }
      }
      if (prior_segment_length) {
        os << setw(width[2 * nb_variable + 4]) << prior_segment_length->cumul[i];
      }
      os << endl;
    }
    os << endl;
  }

  if ((begin_conditonal_entropy) && (end_conditional_entropy) &&
      (change_point_entropy)) {
    os << "\n" << SEQ_label[SEQL_BEGIN_CONDITIONAL_ENTROPY] << ", "
       << SEQ_label[SEQL_END_CONDITIONAL_ENTROPY] << ", "
       << SEQ_label[SEQL_CHANGE_POINT_ENTROPY] << endl;

    os << "\n";
    switch (index_param_type) {
    case TIME :
      os << SEQ_label[SEQL_TIME];
      break;
    case POSITION :
      os << SEQ_label[SEQL_POSITION];
      break;
    default :
      os << SEQ_label[SEQL_INDEX];
      break;
    }

    for (i = 1;i < nb_segment;i++) {
      os << " | " << i + 1 << " " << SEQ_label[SEQL_SEGMENTS];
    }
    os << "   ";
    for (i = 1;i < nb_segment;i++) {
      os << " | " << i + 1 << " " << SEQ_label[SEQL_SEGMENTS];
    }
    os << "   ";
    for (i = 1;i < nb_segment;i++) {
      os << " | " << i + 1 << " " << SEQ_label[SEQL_SEGMENTS];
    }
    os << endl;

    buff = width[nb_variable] - ASCII_SPACE;

    for (i = 0;i < seq_length;i++) {
      os.setf(ios::right , ios::adjustfield);
      os << setw(buff) << seq_index_parameter[i] << "  ";

      os.setf(ios::left , ios::adjustfield);
      for (j = 1;j < nb_segment;j++) {
        os << setw(width[2 * nb_variable + 5]) << begin_conditonal_entropy[j][i];
      }
      os << "   ";
      for (j = 1;j < nb_segment;j++) {
        os << setw(width[2 * nb_variable + 5]) << end_conditional_entropy[j][i];
      }
      os << "   ";
      for (j = 1;j < nb_segment;j++) {
        os << setw(width[2 * nb_variable + 5]) << change_point_entropy[j][i];
      }

      if (i == 0) {
        os.setf(ios::right , ios::adjustfield);
        if (index != I_DEFAULT) {
          os << setw(width[nb_variable + 2]) << identifier[index];
        }
      }
      os << endl;
    }
  }

  if (index_param_type == IMPLICIT_TYPE) {
    delete [] seq_index_parameter;
  }

  delete [] width;

  os.setf(format_flags , ios::adjustfield);

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of segment/state, change-point and entropy profiles for
 *         a single sequence or a sample of sequences (in the case of
 *         multiple change-point models) at the spreadsheet format.
 *
 *  \param[in,out] os                       stream,
 *  \param[in]     index                    sequence index,
 *  \param[in]     nb_segment               number of segments/states,
 *  \param[in]     profiles                 pointer on the segment/state profiles,
 *  \param[in]     label                    profile type label,
 *  \param[in]     common_contrast          flag contrast functions common to the individuals,
 *  \param[in]     piecewise_function       pointer on the piecewise linear functions,
 *  \param[in]     change_point             pointer on the change-point profiles,
 *  \param[in]     segment_length           pointer on the segment length distributions,
 *  \param[in]     prior_segment_length     pointer on the prior segment length distribution assuming a
 *                                          uniform prior on the possible segmentations,
 *  \param[in]     begin_conditonal_entropy pointer on the profiles of entropies conditional on the past,
 *  \param[in]     end_conditional_entropy  pointer on the profiles of entropies conditional on the future,
 *  \param[in]     change_point_entropy     pointer on the change-point entropy profiles.
 */
/*--------------------------------------------------------------*/

ostream& Sequences::profile_spreadsheet_print(ostream &os , int index , int nb_segment ,
                                              double **profiles , const char *label ,
                                              bool common_contrast , double ***piecewise_function ,
                                              long double **change_point ,
                                              Distribution **segment_length ,
                                              Distribution *prior_segment_length ,
                                              long double **begin_conditonal_entropy ,
                                              long double **end_conditional_entropy ,
                                              long double **change_point_entropy) const

{
  int i , j , k;
  int seq_length , max_nb_value , *seq_index_parameter;


  seq_length = length[index == I_DEFAULT ? 0 : index];

  if (index_param_type == IMPLICIT_TYPE) {
    seq_index_parameter = new int[seq_length];
    for (i = 0;i < seq_length;i++) {
      seq_index_parameter[i] = i;
    }
  }
  else {
    seq_index_parameter = index_parameter[index == I_DEFAULT ? 0 : index];
  }

  if (!change_point) {
    os << SEQ_label[SEQL_OPTIMAL] << " " << label << "\t";
  }

  for (i = 1;i < nb_variable;i++) {
//    if ((index == I_DEFAULT) && (nb_sequence * (nb_variable - 1) <= SEQUENCE_MAX_NB_COLUMN)) {
    if (index == I_DEFAULT) {
      if (common_contrast) {
        for (j = 0;j < nb_sequence;j++) {
          os << STAT_label[STATL_VARIABLE] << " " << i << "\t";
        }
        if ((piecewise_function) && (piecewise_function[i])) {
          os << SEQ_label[SEQL_PIECEWISE_LINEAR_FUNCTION] << " " << i << "\t";
        }
      }

      else {
        for (j = 0;j < nb_sequence;j++) {
          os << STAT_label[STATL_VARIABLE] << " " << i << "\t";
          if ((piecewise_function) && (piecewise_function[i])) {
            os << SEQ_label[SEQL_PIECEWISE_LINEAR_FUNCTION] << " " << i << "\t";
          }
        }
      }
    }

//    else if (index != I_DEFAULT) {
    else {
      os << STAT_label[STATL_VARIABLE] << " " << i << "\t";
      if ((piecewise_function) && (piecewise_function[i])) {
        os << SEQ_label[SEQL_PIECEWISE_LINEAR_FUNCTION] << " " << i << "\t";
      }
    }
  }

  switch (index_param_type) {
  case TIME :
    os << SEQ_label[SEQL_TIME];
    break;
  case POSITION :
    os << SEQ_label[SEQL_POSITION];
    break;
  default :
    os << SEQ_label[SEQL_INDEX];
    break;
  }

  for (i = 0;i < nb_segment;i++) {
    os << "\t" << label << " " << i;
  }
  if (change_point) {
    os << "\t";
    for (i = 1;i < nb_segment;i++) {
      os << "\t" << i + 1 << " " << SEQ_label[SEQL_SEGMENTS];
    }
  }
  os << endl;

  for (i = 0;i < seq_length;i++) {
    if (!change_point) {
      os << int_sequence[index == I_DEFAULT ? 0 : index][0][i] << "\t";
    }

    for (j = 1;j < nb_variable;j++) {
//      if ((index == I_DEFAULT) && (nb_sequence * (nb_variable - 1) <= SEQUENCE_MAX_NB_COLUMN)) {
      if (index == I_DEFAULT) {
        if (common_contrast) {
          if (type[j] != REAL_VALUE) {
            for (k = 0;k < nb_sequence;k++) {
              os << int_sequence[k][j][i] << "\t";
            }
          }
          else {
            for (k = 0;k < nb_sequence;k++) {
              os << real_sequence[k][j][i] << "\t";
            }
          }
          if ((piecewise_function) && (piecewise_function[j])) {
            os << piecewise_function[j][0][i] << "\t";
          }
        }

        else {
          if (type[j] != REAL_VALUE) {
            for (k = 0;k < nb_sequence;k++) {
              os << int_sequence[k][j][i] << "\t";
              if ((piecewise_function) && (piecewise_function[j])) {
                os << piecewise_function[j][k][i] << "\t";
              }
            }
          }
          else {
            for (k = 0;k < nb_sequence;k++) {
              os << real_sequence[k][j][i] << "\t";
              if ((piecewise_function) && (piecewise_function[j])) {
                os << piecewise_function[j][k][i] << "\t";
              }
            }
          }
        }
      }

//      else if (index != I_DEFAULT) {
      else {
        if (type[j] != REAL_VALUE) {
          os << int_sequence[index][j][i] << "\t";
        }
        else {
          os << real_sequence[index][j][i] << "\t";
        }
        if ((piecewise_function) && (piecewise_function[j])) {
          os << piecewise_function[j][index][i] << "\t";
        }
      }
    }

    os << seq_index_parameter[i];
    for (j = 0;j < nb_segment;j++) {
      os << "\t" << profiles[i][j];
    }

    if (change_point) {
      os << "\t";
      for (j = 1;j < nb_segment;j++) {
        os << "\t" << change_point[j][i];
      }
    }

    if ((index != I_DEFAULT) && (i == 0)) {
      os << "\t" << identifier[index];
    }
    os << endl;
  }

  if (segment_length) {
    if (prior_segment_length) {
      max_nb_value = prior_segment_length->nb_value;
    }

    else {
      max_nb_value = segment_length[0]->nb_value;
      for (i = 1;i < nb_segment;i++) {
        if (segment_length[i]->nb_value > max_nb_value) {
          max_nb_value = segment_length[i]->nb_value;
        }
      }
    }

    os << "\n\n" << SEQ_label[SEQL_SEGMENT_LENGTH] << " " << STAT_label[STATL_DISTRIBUTIONS] << endl;

    for (i = 0;i < nb_segment;i++) {
      os << "\n" << SEQ_label[SEQL_SEGMENT] << " " << i << " " << SEQ_label[SEQL_LENGTH]
         << " " << STAT_label[STATL_DISTRIBUTION] << endl;
      segment_length[i]->spreadsheet_characteristic_print(os , true);
    }
    if (prior_segment_length) {
      os << "\n" << SEQ_label[SEQL_PRIOR_SEGMENT_LENGTH] << " " << STAT_label[STATL_DISTRIBUTION] << endl;
      prior_segment_length->spreadsheet_characteristic_print(os , true);
    }

    os << "\n";
    for (i = 0;i < nb_segment;i++) {
      os << "\t" << SEQ_label[SEQL_SEGMENT] << " " << i << " " << SEQ_label[SEQL_LENGTH]
         << " " << STAT_label[STATL_DISTRIBUTION];
    }
    if (prior_segment_length) {
      os << "\t" << SEQ_label[SEQL_PRIOR_SEGMENT_LENGTH] << " " << STAT_label[STATL_DISTRIBUTION];
    }
    for (i = 0;i < nb_segment;i++) {
      os << "\t" << SEQ_label[SEQL_SEGMENT] << " " << i << " " << STAT_label[STATL_CUMULATIVE]
         << " " << STAT_label[STATL_DISTRIBUTION] << " " << STAT_label[STATL_FUNCTION];
    }
    if (prior_segment_length) {
      os << "\t" << SEQ_label[SEQL_PRIOR_SEGMENT_LENGTH] << " " << i << " " << STAT_label[STATL_CUMULATIVE]
         << " " << STAT_label[STATL_DISTRIBUTION] << " " << STAT_label[STATL_FUNCTION];
    }
    os << endl;

    for (i = 0;i < max_nb_value;i++) {
      os << i;

      for (j = 0;j < nb_segment;j++) {
        os << "\t";
        if (i < segment_length[j]->nb_value) {
          os << segment_length[j]->mass[i];
        }
      }
      if (prior_segment_length) {
        os << "\t" << prior_segment_length->mass[i];
      }
      for (j = 0;j < nb_segment;j++) {
        os << "\t";
        if (i < segment_length[j]->nb_value) {
          os << segment_length[j]->cumul[i];
        }
      }
      if (prior_segment_length) {
        os << "\t" << prior_segment_length->cumul[i];
      }
      os << endl;
    }
    os << endl;
  }

  if ((begin_conditonal_entropy) && (end_conditional_entropy) &&
      (change_point_entropy)) {
    os << "\n" << SEQ_label[SEQL_BEGIN_CONDITIONAL_ENTROPY];
    for (i = 1;i < nb_segment;i++) {
      os << "\t";
    }
    os << SEQ_label[SEQL_END_CONDITIONAL_ENTROPY];
    for (i = 1;i < nb_segment;i++) {
      os << "\t";
    }
    os << SEQ_label[SEQL_CHANGE_POINT_ENTROPY] << endl;

    os << "\n";
    switch (index_param_type) {
      os << SEQ_label[SEQL_TIME];
      break;
    case POSITION :
      os << SEQ_label[SEQL_POSITION];
      break;
    default :
      os << SEQ_label[SEQL_INDEX];
      break;
    }

    for (i = 1;i < nb_segment;i++) {
      os << "\t" << i + 1 << " " << SEQ_label[SEQL_SEGMENTS];
    }
    os << "\t";
    for (i = 1;i < nb_segment;i++) {
      os << "\t" << i + 1 << " " << SEQ_label[SEQL_SEGMENTS];
    }
    os << "\t";
    for (i = 1;i < nb_segment;i++) {
      os << "\t" << i + 1 << " " << SEQ_label[SEQL_SEGMENTS];
    }
    os << endl;

    for (i = 0;i < seq_length;i++) {
      os << seq_index_parameter[i];

      for (j = 1;j < nb_segment;j++) {
        os << "\t" << begin_conditonal_entropy[j][i];
      }
      os << "\t";
      for (j = 1;j < nb_segment;j++) {
        os << "\t" << end_conditional_entropy[j][i];
      }
      os << "\t";
      for (j = 1;j < nb_segment;j++) {
        os << "\t" << change_point_entropy[j][i];
      }

      if ((index != I_DEFAULT) && (i == 0)) {
        os << "\t" << identifier[index];
      }
      os << endl;
    }
  }

  if ((index != I_DEFAULT) && (index_param_type == IMPLICIT_TYPE)) {
    delete [] seq_index_parameter;
  }

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of segment/state, change-point and entropy profiles for
 *         a single sequence or a sample of sequences (in the case of
 *         multiple change-point models) at the Gnuplot format.
 *
 *  \param[in,out] os                       stream,
 *  \param[in]     index                    sequence index,
 *  \param[in]     nb_segment               number of segments/states,
 *  \param[in]     profiles                 pointer on the segment/state profiles,
 *  \param[in]     common_contrast          flag contrast functions common to the individuals,
 *  \param[in]     piecewise_function       pointer on the piecewise linear functions,
 *  \param[in]     change_point             pointer on the change-point profiles,
 *  \param[in]     segment_length           pointer on the segment length distributions,
 *  \param[in]     prior_segment_length     pointer on the prior segment length distribution assuming a
 *                                          uniform prior on the possible segmentations,
 *  \param[in]     begin_conditonal_entropy pointer on the profiles of entropies conditional on the past,
 *  \param[in]     end_conditional_entropy  pointer on the profiles of entropies conditional on the future,
 *  \param[in]     change_point_entropy     pointer on the change-point entropy profiles.
 */
/*--------------------------------------------------------------*/

ostream& Sequences::profile_plot_print(ostream &os , int index , int nb_segment ,
                                       double **profiles , bool common_contrast ,
                                       double ***piecewise_function , long double **change_point ,
                                       Distribution **segment_length ,
                                       Distribution *prior_segment_length ,
                                       long double **begin_conditonal_entropy ,
                                       long double **end_conditional_entropy ,
                                       long double **change_point_entropy) const

{
  int i , j , k;
  int seq_length , *seq_index_parameter;


  seq_length = length[index == I_DEFAULT ? 0 : index];

  if (index_param_type == IMPLICIT_TYPE) {
    seq_index_parameter = new int[seq_length];
    for (i = 0;i < seq_length;i++) {
      seq_index_parameter[i] = i;
    }
  }
  else {
    seq_index_parameter = index_parameter[index == I_DEFAULT ? 0 : index];
  }

  for (i = 0;i < seq_length;i++) {
    os << seq_index_parameter[i] << " ";

    for (j = 1;j < nb_variable;j++) {
      if ((piecewise_function) && (piecewise_function[j])) {
        if ((index != I_DEFAULT) || (!common_contrast)) {
          if (type[j] != REAL_VALUE) {
            for (k = 0;k < nb_sequence;k++) {
              if ((index == I_DEFAULT) || (index == k)) {
                os << int_sequence[k][j][i] << " " << piecewise_function[j][k][i] << " ";
              }
            }
          }
          else {
            for (k = 0;k < nb_sequence;k++) {
              if ((index == I_DEFAULT) || (index == k)) {
                os << real_sequence[k][j][i] << " " << piecewise_function[j][k][i] << " ";
              }
            }
          }
        }

        else {
          if (type[j] != REAL_VALUE) {
            for (k = 0;k < nb_sequence;k++) {
              os << int_sequence[k][j][i] << " ";
            }
          }
          else {
            for (k = 0;k < nb_sequence;k++) {
              os << real_sequence[k][j][i] << " ";
            }
          }

          os << piecewise_function[j][0][i] << " ";
        }
      }
    }

    for (j = 0;j < nb_segment;j++) {
      os << profiles[i][j] << " ";
    }

    if (change_point) {
      for (j = 1;j < nb_segment;j++) {
        os << change_point[j][i] << " ";
      }
    }

    if (segment_length) {
      for (j = 0;j < nb_segment;j++) {
        if (i < segment_length[j]->nb_value) {
          os << segment_length[j]->mass[i] << " ";
        }
        else {
          os << 0. << " ";
        }
      }
      if (prior_segment_length) {
        if (i < prior_segment_length->nb_value) {
          os << prior_segment_length->mass[i] << " ";
        }
        else {
          os << 0. << " ";
        }
      }

/*      for (j = 0;j < nb_segment;j++) {
        if (i < segment_length[j]->nb_value) {
          os << segment_length[j]->cumul[i] << " ";
        }
        else {
          os << 1. << " ";
        }
      }
      if (prior_segment_length) {
        if (i < prior_segment_length->nb_value) {
          os << prior_segment_length->cumul[i] << " ";
        }
        else {
          os << 1. << " ";
        }
      } */
    }

    if ((begin_conditonal_entropy) && (end_conditional_entropy) &&
        (change_point_entropy)) {
      for (j = 1;j < nb_segment;j++) {
        os << begin_conditonal_entropy[j][i] << " ";
      }
      for (j = 1;j < nb_segment;j++) {
        os << end_conditional_entropy[j][i] << " ";
      }
/*      for (j = 1;j < nb_segment;j++) {
        os << change_point_entropy[j][i] << " ";
      } */
      os << change_point_entropy[nb_segment - 1][i];
    }
    os << endl;
  }

  if (index_param_type == IMPLICIT_TYPE) {
    delete [] seq_index_parameter;
  }

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of change-point profiles for a single sequence or
 *         a sample of sequences (in the case of multiple change-point models)
 *         at the "plotable" format.
 *
 *  \param[in] plot         reference on a MultiPlot object,
 *  \param[in] index        sequence index,
 *  \param[in] nb_segment   number of segments,
 *  \param[in] change_point pointer on the change-point profiles.
 */
/*--------------------------------------------------------------*/

void Sequences::change_point_profile_plotable_write(MultiPlot &plot , int index , int nb_segment ,
                                                    long double **change_point) const

{
  int i , j , k;
  int seq_length , *seq_index_parameter;


  seq_length = length[index == I_DEFAULT ? 0 : index];

  if (index_param_type == IMPLICIT_TYPE) {
    seq_index_parameter = new int[seq_length];
    for (i = 0;i < seq_length;i++) {
      seq_index_parameter[i] = i;
    }
  }
  else {
    seq_index_parameter = index_parameter[index == I_DEFAULT ? 0 : index];
  }

  plot.resize(MAX(nb_segment - 1 , 3));
  i = 0;

  for (j = MAX(1 , nb_segment - 3);j < nb_segment;j++) {
    for (k = 0;k < seq_length;k++) {
      plot[i].add_point(seq_index_parameter[k] , change_point[j][k]);
    }
    i++;
  }

  if (index_param_type == IMPLICIT_TYPE) {
    delete [] seq_index_parameter;
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of segment length distributions including the prior segment length
 *         distribution corresponding to a uniform prior on the possible segmentations 
 *         at the "plotable" format.
 *
 *  \param[in] plot                 reference on a MultiPlot object,
 *  \param[in] nb_segment           number of segments,
 *  \param[in] segment_length_max   maximum probability of segment length distribution,
 *  \param[in] segment_length       pointer on the segment length distributions,
 *  \param[in] prior_segment_length pointer on the prior segment length distribution.
 */
/*--------------------------------------------------------------*/

void Sequences::segment_length_distribution_plotable_write(MultiPlot &plot , int nb_segment ,
                                                           double segment_length_max ,
                                                           Distribution **segment_length ,
                                                           Distribution *prior_segment_length) const

{
  int i;
  int max_nb_value;


  if (prior_segment_length) {
    max_nb_value = prior_segment_length->nb_value;
 
    plot.resize(nb_segment + 1);
  }

  else {
    max_nb_value = segment_length[0]->nb_value;
    for (i = 1;i < nb_segment;i++) {
      if (segment_length[i]->nb_value > max_nb_value) {
        max_nb_value = segment_length[i]->nb_value;
      }
    }

    plot.resize(nb_segment);
  }

  plot.xrange = Range(0 , max_nb_value - 1);
  if (max_nb_value - 1 < TIC_THRESHOLD) {
    plot.xtics = 1;
  }

  plot.yrange = Range(0 , segment_length_max * YSCALE);

  for (i = 0;i < nb_segment;i++) {
    segment_length[i]->plotable_mass_write(plot[i]);
  }
  if (prior_segment_length) {
    prior_segment_length->plotable_mass_write(plot[nb_segment]);
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of entropy profiles for a single sequence or
 *         a sample of sequences (in the case of multiple change-point models)
 *         at the "plotable" format.
 *
 *  \param[in] plot                     reference on a MultiPlot object,
 *  \param[in] index                    sequence index,
 *  \param[in] begin_conditonal_entropy pointer on the profiles of entropies conditional on the past,
 *  \param[in] end_conditional_entropy  pointer on the profiles of entropies conditional on the future,
 *  \param[in] change_point_entropy     pointer on the change-point entropy profiles.
 */
/*--------------------------------------------------------------*/

void Sequences::entropy_profile_plotable_write(MultiPlot &plot , int index ,
                                               long double *begin_conditional_entropy ,
                                               long double *end_conditional_entropy ,
                                               long double *change_point_entropy) const

{
  int i;
  int seq_length , *seq_index_parameter;


  seq_length = length[index == I_DEFAULT ? 0 : index];

  if (index_param_type == IMPLICIT_TYPE) {
    seq_index_parameter = new int[seq_length];
    for (i = 0;i < seq_length;i++) {
      seq_index_parameter[i] = i;
    }
  }
  else {
    seq_index_parameter = index_parameter[index == I_DEFAULT ? 0 : index];
  }

  plot.resize(3);

  for (i = 0;i < seq_length;i++) {
    plot[0].add_point(seq_index_parameter[i] , begin_conditional_entropy[i]);
    plot[1].add_point(seq_index_parameter[i] , end_conditional_entropy[i]);
    plot[2].add_point(seq_index_parameter[i] , change_point_entropy[i]);
  }

  if (index_param_type == IMPLICIT_TYPE) {
    delete [] seq_index_parameter;
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation by summation of segments/change-point profiles for
 *         a single sequence or a sample of sequences and of entropy profiles.
 *
 *  \param[in] index              sequence index,
 *  \param[in] nb_segment         number of segments,
 *  \param[in] model_type         segment model types,
 *  \param[in] common_contrast    flag contrast functions common to the individuals,
 *  \param[in] shape_parameter    negative binomial shape parameters,
 *  \param[in] rank               ranks (for ordinal variables),
 *  \param[in] os                 stream,
 *  \param[in] plot_set           pointer on a MultiPlotSet object,
 *  \param[in] segment_length_max maximum probability of segment length distribution,
 *  \param[in] output             output type,
 *  \param[in] format             output format (ASCII/SPREADSHEET/GNUPLOT/PLOT).
 *
 *  \return                       log-likelihood of the multiple change-point model.
 */
/*--------------------------------------------------------------*/

double Sequences::forward_backward(int index , int nb_segment , segment_model *model_type ,
                                   bool common_contrast , double *shape_parameter ,
                                   double **rank , ostream *os , MultiPlotSet *plot_set ,
                                   double &segment_length_max , change_point_profile output ,
                                   output_format format) const

{
  int i , j , k , m;
  int seq_length , *inf_bound_parameter , *seq_index_parameter;
  double sum , buff , rlikelihood , *likelihood , **seq_mean , **hyperparam , **backward_output ,
         ***factorial , ***binomial_coeff , ***smoothed;
  long double segment_norm , sequence_norm , lbuff , lsum , segmentation_entropy , first_order_entropy ,
              change_point_entropy_sum , marginal_entropy , *contrast , *normalized_contrast ,
              *norm , *forward_norm , *backward_norm , *entropy_smoothed , *segment_predicted ,
              **forward , **backward , **change_point , **forward_predicted_entropy ,
              **backward_predicted_entropy , **forward_partial_entropy , **backward_partial_entropy ,
              **change_point_entropy , ***state_entropy;
  Distribution **segment_length;
  DiscreteParametric *prior_segment_length;

# ifdef DEBUG
  long double *entropy_norm;
# endif


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

    if (((model_type[i - 1] == LINEAR_MODEL_CHANGE) || (output == CHANGE_POINT)) && (!seq_index_parameter)) {
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

  seq_length = length[index == I_DEFAULT ? 0 : index];
  contrast = new long double[seq_length];
  normalized_contrast = new long double[seq_length];

  forward = new long double*[seq_length];
  for (i = 0;i < seq_length;i++) {
    forward[i] = new long double[nb_segment];
  }

  segment_predicted = new long double[seq_length];

  forward_predicted_entropy = new long double*[seq_length];
  for (i = 0;i < seq_length;i++) {
    forward_predicted_entropy[i] = new long double[nb_segment];
  }

  norm = new long double[seq_length];
  forward_norm = new long double[seq_length];
//  entropy_norm = new long double[seq_length];

  likelihood = new double[nb_segment];

  backward = new long double*[seq_length];
  for (i = 0;i < seq_length;i++) {
    backward[i] = new long double[nb_segment];
  }

  backward_predicted_entropy = new long double*[seq_length];
  for (i = 0;i < seq_length;i++) {
    backward_predicted_entropy[i] = new long double[nb_segment];
  }

  backward_norm = new long double[seq_length];

  smoothed = new double**[nb_segment];
  for (i = 1;i < nb_segment;i++) {
    smoothed[i] = new double*[seq_length];
    for (j = 0;j < seq_length;j++) {
      smoothed[i][j] = new double[nb_segment];
    }
  }

  backward_output = new double*[seq_length];
  for (i = 0;i < seq_length;i++) {
    backward_output[i] = new double[nb_segment];
  }

  segment_length = new Distribution*[nb_segment];
  for (i = 0;i < nb_segment;i++) {
    segment_length[i] = new Distribution(seq_length);
    for (j = 0;j < seq_length;j++) {
      segment_length[i]->mass[j] = 0.;
    }
  }

  change_point = new long double*[nb_segment];
  for (i = 1;i < nb_segment;i++) {
    change_point[i] = new long double[seq_length];
  }

  entropy_smoothed = new long double[nb_segment];

  state_entropy = new long double**[nb_segment];
  for (i = 1;i < nb_segment;i++) {
    state_entropy[i] = new long double*[seq_length];
    for (j = 0;j < seq_length;j++) {
      state_entropy[i][j] = new long double[nb_segment];
    }
  }

  forward_partial_entropy = new long double*[nb_segment];
  for (i = 1;i < nb_segment;i++) {
    forward_partial_entropy[i] = new long double[seq_length];
  }

  backward_partial_entropy = new long double*[nb_segment];
  for (i = 1;i < nb_segment;i++) {
    backward_partial_entropy[i] = new long double[seq_length];
  }

  change_point_entropy = new long double*[nb_segment];
  for (i = 1;i < nb_segment;i++) {
    change_point_entropy[i] = new long double[seq_length];
  }

  // forward recurrence

  for (i = 0;i < seq_length;i++) {

    // computation of segment log-likelihoods

    forward_contrast(i , index , model_type , common_contrast , factorial ,
                     shape_parameter , binomial_coeff , seq_mean , seq_index_parameter ,
                     hyperparam , rank , contrast);

    // recurrence and computation of predicted entropies

    if (contrast[i] != D_INF) {
      contrast[i] = expl(contrast[i]);
    }
    else {
      contrast[i] = 0.;
    }

    segment_norm = 0.;
    for (j = i - 1;j >= 0;j--) {
      segment_norm += norm[j];

#     ifdef DEBUG
      if (i == seq_length - 1) {
        cout << j << ": " << contrast[j] << " " << segment_norm << " | ";
      }
#     endif

      if (contrast[j] != D_INF) {
        contrast[j] = expl(contrast[j] - segment_norm);
      }
      else {
        contrast[j] = 0.;
      }

#     ifdef DEBUG
      if (i == seq_length - 1) {
        cout << contrast[j];
        if (j > 0) {
          cout << " " << forward[j - 1][nb_segment - 2] << " | "
               << contrast[j] * forward[j - 1][nb_segment - 2];
        }
        cout << endl;
      }
#     endif

    }

#   ifdef DEBUG
    if (i == seq_length - 1) {
      cout << endl;
    }
#   endif

    for (j = 0;j < nb_segment;j++) {
      forward[i][j] = 0.;
      forward_predicted_entropy[i][j] = 0.;
    }
    norm[i] = 0.;

//    for (j = MAX(0 , nb_segment + i - seq_length);j < MIN((i < seq_length - 1 ? nb_segment - 1 : nb_segment) , i + 1);j++) {
    for (j = 0;j < MIN((i < seq_length - 1 ? nb_segment - 1 : nb_segment) , i + 1);j++) {
      if (j == 0) {
        forward[i][j] = contrast[0];
      }

      else {
        for (k = i;k >= j;k--) {
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
//      for (j = MAX(0 , nb_segment + i - seq_length);j < MIN((i < seq_length - 1 ? nb_segment - 1 : nb_segment) , i + 1);j++) {
      for (j = 0;j < MIN((i < seq_length - 1 ? nb_segment - 1 : nb_segment) , i + 1);j++) {
        forward[i][j] /= norm[i];
      }

      norm[i] = logl(norm[i]);
    }
//    entropy_norm[i] = norm[i];

    forward_norm[i] = segment_norm + norm[i];
  }

# ifdef DEBUG
  cout << "\n";
  for (i = 0;i < seq_length;i++) {
    cout << i << " |";
    lsum = 0.;
    for (j = 0;j < nb_segment;j++) {
      lsum += forward[i][j];
      cout << " " << forward[i][j];
    }
    cout << " | " << lsum << ", " << expl(norm[i]) << endl;
  }
# endif

# ifdef DEBUG
  cout << "\n";
  for (i = 0;i < seq_length;i++) {
    cout << i << " |";
    for (j = 0;j < nb_segment;j++) {
      cout << " " << forward_predicted_entropy[i][j];
    }
    cout << endl;
  }
# endif

  // extraction of the log-likelihoods of the observed sequence for the different numbers of segments

  for (i = 0;i < nb_segment;i++) {
    if (forward[seq_length - 1][i] > 0.) {
      likelihood[i] = logl(forward[seq_length - 1][i]) + forward_norm[seq_length - 1];
    }
    else {
      likelihood[i] = D_INF;
    }
  }

  rlikelihood = likelihood[nb_segment - 1];

  if (rlikelihood != D_INF) {

#   ifdef MESSAGE
    segmentation_entropy = rlikelihood;
#   endif

    for (i = 1;i < nb_segment;i++) {
      for (j = 0;j < seq_length;j++) {
        for (k = 0;k < nb_segment;k++) {
          state_entropy[i][j][k] = 0.;
        }
      }
    }

    // backward recurrence

    for (i = seq_length - 1;i >= 0;i--) {

      // computation of segment log-likelihoods

      backward_contrast(i , index , model_type , common_contrast , factorial ,
                        shape_parameter , binomial_coeff , seq_mean , seq_index_parameter ,
                        hyperparam , rank , contrast);

      // recurrence and computation of predicted entropies

      if (contrast[i] != D_INF) {
        normalized_contrast[i] = expl(contrast[i]);
      }
      else {
        normalized_contrast[i] = 0.;
      }

      segment_norm = 0.;
      for (j = i + 1;j < seq_length;j++) {
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
        backward_output[i][j] = 0.;

        for (k = 1;k < nb_segment;k++) {
          smoothed[k][i][j] = 0.;
        }
      }
      norm[i] = 0.;

//      for (j = MAX((i == 0 ? 0 : 1) , nb_segment + i - seq_length);j < MIN(nb_segment , i + 1);j++) {
      for (j = MAX((i == 0 ? 0 : 1) , nb_segment + i - seq_length);j < nb_segment;j++) {
        if (j < nb_segment - 1) {
          for (k = i;k <= seq_length + j - nb_segment;k++) {
            segment_predicted[k] = normalized_contrast[k] * backward[k + 1][j + 1];
            backward[i][j] += segment_predicted[k];
          }

          if (backward[i][j] > 0.) {
            for (k = i;k <= seq_length + j - nb_segment;k++) {
              lbuff = segment_predicted[k] / backward[i][j];
              if (lbuff > 0.) {
                backward_predicted_entropy[i][j] += lbuff * (backward_predicted_entropy[k + 1][j + 1] - logl(lbuff));
              }
            }
          }
        }

        else {
          backward[i][j] = normalized_contrast[seq_length - 1];
        }

        norm[i] += backward[i][j];
      }

      if (norm[i] > 0.) {
//        for (j = MAX((i == 0 ? 0 : 1) , nb_segment + i - seq_length);j < MIN(nb_segment , i + 1);j++) {
        for (j = MAX((i == 0 ? 0 : 1) , nb_segment + i - seq_length);j < nb_segment;j++) {
          backward[i][j] /= norm[i];
        }

        norm[i] = logl(norm[i]);
      }

      backward_norm[i] = segment_norm + norm[i];

      // extraction of the smoothed probabilities for the different numbers of segments

      if (i < seq_length - 1) {
        for (j = 1;j < nb_segment;j++) {
          sequence_norm = expl(forward_norm[i] + backward_norm[i + 1] - likelihood[j]);

          for (k = MAX(0 , j + i - seq_length - 1);k <= MIN(j , i);k++) {
            smoothed[j][i][k] = smoothed[j][i + 1][k];
            if (k > 0) {
              smoothed[j][i][k] -= forward[i][k - 1] * backward[i + 1][k + nb_segment - j - 1] *
                                   sequence_norm;
            }
            if (k < j) {
              smoothed[j][i][k] += forward[i][k] * backward[i + 1][k + nb_segment - j] *
                                   sequence_norm;
            }

            if (smoothed[j][i][k] < 0.) {
              smoothed[j][i][k] = 0.;
            }
            if (smoothed[j][i][k] > 1.) {
              smoothed[j][i][k] = 1.;
            }
          }
        }
      }

      else {
        for (j = 1;j < nb_segment;j++) {
          smoothed[j][i][j] = 1.;
        }
      }

      if (i == 0) {
        sequence_norm = expl(backward_norm[i] - rlikelihood);
      }
      else {
        sequence_norm = expl(forward_norm[i - 1] + backward_norm[i] - rlikelihood);

#       ifdef DEBUG
        cout << i << ": " << forward_norm[i - 1] << " " << backward_norm[i] << " | "
             << forward_norm[i - 1] + backward_norm[i] - rlikelihood << " " << sequence_norm << endl;
#       endif

      }

      if (output == SEGMENT) {
        for (j = 0;j < nb_segment;j++) {
          backward_output[i][j] = smoothed[nb_segment - 1][i][j];
        }
      }

      // computation of posterior change-point probabilities for the different numbers of segments

      if (i == 0) {

#       ifdef MESSAGE
        lbuff = backward[i][0] * sequence_norm;
        if ((lbuff < 1. - DOUBLE_ERROR) || (lbuff > 1. + DOUBLE_ERROR)) {
          cout << "\nERROR: " << lbuff << " | " << 1 << endl;
        }
#       endif

        if (output == CHANGE_POINT) {
          backward_output[i][0] = 1.;
        }
        for (j = 1;j < nb_segment;j++) {
          change_point[j][i] = 1.;
        }
      }

      else {
        change_point[nb_segment - 1][i] = 0.;
        for (j = MAX(1 , nb_segment + i - seq_length);j < MIN(nb_segment , i + 1);j++) {
          if (output == CHANGE_POINT) {
            backward_output[i][j] = forward[i - 1][j - 1] * backward[i][j] * sequence_norm;
          }
          change_point[nb_segment - 1][i] += forward[i - 1][j - 1] * backward[i][j];
        }
        change_point[nb_segment - 1][i] *= sequence_norm;

        for (j = 1;j < nb_segment - 1;j++) {
          change_point[j][i] = 0.;
          for (k = MAX(1 , j + 1 + i - seq_length);k <= MIN(j , i);k++) {
            change_point[j][i] += forward[i - 1][k - 1] * backward[i][k + nb_segment - j - 1];
          }
          change_point[j][i] *= expl(forward_norm[i - 1] + backward_norm[i] - likelihood[j]);
        }
      }

      segment_norm = 0.;
      for (j = i;j < seq_length;j++) {
        segment_norm += norm[j];
        if (contrast[j] != D_INF) {
          normalized_contrast[j] = expl(contrast[j] - segment_norm);
        }
        else {
          normalized_contrast[j] = 0.;
        }
      }

      // extraction of segment length distributions

      if (i == 0) {
        for (j = i;j <= seq_length - nb_segment;j++) {
          segment_length[0]->mass[j + 1] = normalized_contrast[j] * backward[j + 1][1] * sequence_norm;
        }
      }

      else {
        for (j = MAX(1 , nb_segment + i - seq_length);j < MIN(nb_segment , i + 1);j++) {
          if (j < nb_segment - 1) {
            if (backward[i][j] > 0.) {
              for (k = i;k <= seq_length + j - nb_segment;k++) {
                segment_length[j]->mass[k - i + 1] += forward[i - 1][j - 1] * normalized_contrast[k] *
                                                      backward[k + 1][j + 1] * sequence_norm;

              }
            }
          }

          else {
//            segment_length[j]->mass[seq_length - i] = forward[i - 1][j - 1] * backward[i][j] * sequence_norm;
            segment_length[j]->mass[seq_length - i] = forward[i - 1][j - 1] * normalized_contrast[seq_length - 1] *
                                                      sequence_norm;
          }
        }
      }

      // computation of partial entropies for the different numbers of segments

      if (i > 0) {
        for (j = 1;j < nb_segment;j++) {
          sequence_norm = expl(forward_norm[i - 1] + backward_norm[i] - likelihood[j]);

          for (k = MAX(1 , j + 1 + i - seq_length);k <= MIN(j , i);k++) {
            if (k < j) {
              lsum = 0.;
              for (m = seq_length + k - nb_segment;m >= i;m--) {
                lsum += normalized_contrast[m] * backward[m + 1][k + nb_segment - j];
                if (smoothed[j][m][k] > 0.) {
                  lbuff = forward[i - 1][k - 1] * lsum * sequence_norm / smoothed[j][m][k];
                  if (lbuff > 0.) {
                    state_entropy[j][m][k] += lbuff * (forward_predicted_entropy[i - 1][k - 1] - logl(lbuff));
                  }
                }
              }
            }
 
            else {
              lsum = forward[i - 1][k - 1] * normalized_contrast[seq_length - 1] *
                     sequence_norm;
              for (m = seq_length - 1;m >= i;m--) {
                if (smoothed[j][m][k] > 0.) {
                  lbuff = lsum / smoothed[j][m][k];
                  if (lbuff > 0.) {
                    state_entropy[j][m][k] += lbuff * (forward_predicted_entropy[i - 1][k - 1] - logl(lbuff));
                  }
                }
              }
            }
          }
        }
      }

      // computation of the segmentation entropy

#     ifdef MESSAGE
      if (i == 0) {
        for (j = i;j <= seq_length - nb_segment;j++) {
          if (contrast[j] != D_INF) {
            segmentation_entropy -= normalized_contrast[j] * backward[j + 1][1] *
                                    sequence_norm * contrast[j];
          }
        }
      }

      else {
        for (j = MAX(1 , nb_segment + i - seq_length);j < MIN(nb_segment , i + 1);j++) {
          if (j < nb_segment - 1) {
            for (k = i;k <= seq_length + j - nb_segment;k++) {
              if (contrast[k] != D_INF) {
                segmentation_entropy -= forward[i - 1][j - 1] * normalized_contrast[k] * backward[k + 1][j + 1] *
                                        sequence_norm * contrast[k];
              }
            }
          }

          else {
            if (contrast[seq_length - 1] != D_INF) {
              segmentation_entropy -= forward[i - 1][j - 1] * normalized_contrast[seq_length - 1] *
                                      sequence_norm * contrast[seq_length - 1];
            }
          }
        }
      }
#     endif

    }

//    segmentation_entropy = forward_predicted_entropy[seq_length - 1][nb_segment - 1];
//    segmentation_entropy = backward_predicted_entropy[0][0];

#   ifdef DEBUG
    cout << "\n";
//    for (i = seq_length - 1;i >= 0;i--) {
    for (i = 0;i < seq_length;i++) {
      cout << i << " |";
      lsum = 0.;
      for (j = 0;j < nb_segment;j++) {
        lsum += backward[i][j];
        cout << " " << backward[i][j];
      }
      cout << " | " << lsum << ", " << expl(norm[i]) << endl;
    }
#   endif

#   ifdef DEBUG
    cout << "\n";
//    for (i = seq_length - 1;i >= 0;i--) {
    for (i = 0;i < seq_length;i++) {
      cout << i << " |";
      for (j = 0;j < nb_segment;j++) {
        cout << " " << backward_predicted_entropy[i][j];
      }
      cout << endl;
    }
#   endif

#   ifdef DEBUG
    for (i = 1;i < nb_segment;i++) {
      cout << "\n" << i + 1 << " " << SEQ_label[SEQL_SEGMENTS] << endl;
      for (j = 0;j < seq_length;j++) {
        cout << j << " |";
        for (k = 0;k < nb_segment;k++) {
          cout << " " << state_entropy[i][j][k];
        }
        cout << endl;
      }
    }
#   endif

#   ifdef MESSAGE
    for (i = 1;i < nb_segment;i++) {
      for (j = 0;j < seq_length - 1;j++) {
        sum = 0.;
        for (k = 0;k < nb_segment;k++) {
          sum += smoothed[i][j][k];
        }
        if ((sum < 1. - DOUBLE_ERROR) || (sum > 1. + DOUBLE_ERROR)) {
          cout << "\nERROR: " << i << " | " << sum << endl;
        }
      }
    }

    for (i = 1;i < nb_segment;i++) {
      sum = 0.;
      for (j = 0;j < seq_length;j++) {
        sum += change_point[i][j];
      }
      if ((sum < i + 1 - DOUBLE_ERROR) || (sum > i + 1 + DOUBLE_ERROR)) {
        cout << "\nERROR: " << sum << " | " << i + 1 << endl;
      }
    }
#   endif

    for (i = 0;i < nb_segment;i++) {
      segment_length[i]->nb_value_computation();
      segment_length[i]->offset_computation();
      segment_length[i]->cumul_computation();
      segment_length[i]->max_computation();
      segment_length[i]->mean_computation();
      segment_length[i]->variance_computation();

#     ifdef DEBUG
      cout << "\n" << SEQ_label[SEQL_SEGMENT] << " " << i << " " << SEQ_label[SEQL_LENGTH] << " "
           << STAT_label[STATL_DISTRIBUTION] << "   ";
      segment_length[i]->ascii_characteristic_print(cout);
      cout << "\n";
      segment_length[i]->ascii_print(cout , false , true , false);
#     endif

    }

    prior_segment_length = new DiscreteParametric(PRIOR_SEGMENT_LENGTH , nb_segment , seq_length , D_DEFAULT , D_DEFAULT);

    segment_length_max = prior_segment_length->max;
    for (i = 0;i < nb_segment;i++) {
      if (segment_length[i]->max > segment_length_max) {
        segment_length_max = segment_length[i]->max;
      }
    }

    // extraction of partial entropies for the different numbers of segments

    for (i = 1;i < nb_segment;i++) {
      for (j = 0;j < seq_length;j++) {
        forward_partial_entropy[i][j] = 0.;
        for (k = 0;k < nb_segment;k++) {
          if (state_entropy[i][j][k] < 0.) {
            state_entropy[i][j][k] = 0.;
          }
          if (smoothed[i][j][k] > 0.) {
            forward_partial_entropy[i][j] += smoothed[i][j][k] * (state_entropy[i][j][k] -
                                              log(smoothed[i][j][k]));
          }
        }
        if (forward_partial_entropy[i][j] < 0.) {
          forward_partial_entropy[i][j] = 0.;
        }
      }
    }

    // computation of ordered change-point entropy and of marginal entropy

    for (i = 0;i < nb_segment - 1;i++) {
      entropy_smoothed[i] = 0.;
    }
    entropy_smoothed[nb_segment - 1] = 1.;

    first_order_entropy = 0.;
    marginal_entropy = 0.;

    for (i = seq_length - 2;i >= 0;i--) {
      sequence_norm = expl(forward_norm[i] + backward_norm[i + 1] - rlikelihood);

/*      for (j = MIN(nb_segment - 1 , i + 1) + 1;j < nb_segment;j++) {
        entropy_smoothed[j] = 0.;
      } */

//      for (j = 0;j < nb_segment;j++) {
      for (j = MAX(0 , nb_segment + i - seq_length);j <= MIN(nb_segment - 1 , i + 1);j++) {
        if (j > 0) {
//          entropy_smoothed[j] -= forward[i][j - 1] * backward[i + 1][j] * sequence_norm;
          lbuff = forward[i][j - 1] * backward[i + 1][j] * sequence_norm;
          entropy_smoothed[j] -= lbuff;
          if ((lbuff > 0.) && (lbuff < 1.)) {
            first_order_entropy -= lbuff * logl(lbuff);
          }
        }
        if ((entropy_smoothed[j] > 0.) && (entropy_smoothed[j] < 1.)) {
          first_order_entropy -= entropy_smoothed[j] * logl(entropy_smoothed[j]);
        }

        if (j < nb_segment - 1) {
          entropy_smoothed[j] += forward[i][j] * backward[i + 1][j + 1] * sequence_norm;
/*          lbuff = forward[i][j] * backward[i + 1][j + 1] * sequence_norm;
          entropy_smoothed[j] += lbuff;
          if ((lbuff > 0.) && (lbuff < 1.)) {
            first_order_entropy -= lbuff * logl(lbuff);
          } */
        }

        if (entropy_smoothed[j] < 0.) {
          entropy_smoothed[j] = 0.;
        }
        if (entropy_smoothed[j] > 1.) {
          entropy_smoothed[j] = 1.;
        }

        if (entropy_smoothed[j] > 0.) {
          first_order_entropy += entropy_smoothed[j] * logl(entropy_smoothed[j]);
          marginal_entropy -= entropy_smoothed[j] * logl(entropy_smoothed[j]);
        }
      }

#     ifdef MESSAGE
      sum = 0.;
      for (j = 0;j < nb_segment;j++) {
        sum += entropy_smoothed[j];
      }
      if ((sum < 1. - DOUBLE_ERROR) || (sum > 1. + DOUBLE_ERROR)) {
        cout << "\nERROR: " << nb_segment << " " << i << " | " << sum << endl;
      }
#     endif

    }

    // computation of change-point entropy profile and change-point entropy

    change_point_entropy[nb_segment - 1][0] = 0.;
    change_point_entropy_sum = 0.;
    for (i = 1;i < seq_length;i++) {
      if ((change_point[nb_segment - 1][i] > 0.) && (change_point[nb_segment - 1][i] < 1.)) {

        change_point_entropy[nb_segment - 1][i] = -change_point[nb_segment - 1][i] * logl(change_point[nb_segment - 1][i]) -
                                                  (1 - change_point[nb_segment - 1][i]) * logl(1 - change_point[nb_segment - 1][i]);
        change_point_entropy_sum += change_point_entropy[nb_segment - 1][i];
      }
      else {
        change_point_entropy[nb_segment - 1][i] = 0.;
      }
    }

    // computation of change-point entropy profiles for the different numbers of segments

    for (i = 1;i < nb_segment - 1;i++) {
      change_point_entropy[i][0] = 0.;
      for (j = 1;j < seq_length;j++) {
        if ((change_point[i][j] > 0.) && (change_point[i][j] < 1.)) {
          change_point_entropy[i][j] = -change_point[i][j] * logl(change_point[i][j]) -
                                       (1 - change_point[i][j]) * logl(1 - change_point[i][j]);
        }
        else {
          change_point_entropy[i][j] = 0.;
        }
      }
    }

    // supplementary forward recurrence for computing partial entropies

    for (i = 1;i < nb_segment;i++) {
      for (j = 0;j < seq_length;j++) {
        for (k = 0;k < nb_segment;k++) {
          state_entropy[i][j][k] = 0.;
        }
      }
    }

    for (i = 0;i < seq_length;i++) {

      // computation of segment log-likelihoods

      forward_contrast(i , index , model_type , common_contrast , factorial ,
                       shape_parameter , binomial_coeff , seq_mean , seq_index_parameter ,
                       hyperparam , rank , contrast);

      // recurrence

      if (contrast[i] != D_INF) {
        normalized_contrast[i] = expl(contrast[i]);
      }
      else {
        normalized_contrast[i] = 0.;
      }

      segment_norm = 0.;
      for (j = i - 1;j >= 0;j--) {
        segment_norm += norm[j];
        if (contrast[j] != D_INF) {
          normalized_contrast[j] = expl(contrast[j] - segment_norm);
        }
        else {
          normalized_contrast[j] = 0.;
        }
      }

      for (j = 0;j < nb_segment;j++) {
        forward[i][j] = 0.;
      }
      norm[i] = 0.;

//      for (j = MAX(0 , nb_segment + i - seq_length);j < MIN((i < seq_length - 1 ? nb_segment - 1 : nb_segment) , i + 1);j++) {
      for (j = 0;j < MIN((i < seq_length - 1 ? nb_segment - 1 : nb_segment) , i + 1);j++) {
        if (j == 0) {
          forward[i][j] = normalized_contrast[0];
        }
        else {
          for (k = i;k >= j;k--) {
            forward[i][j] += normalized_contrast[k] * forward[k - 1][j - 1];
          }
        }

        norm[i] += forward[i][j];
      }

      if (norm[i] > 0.) {
//        for (j = MAX(0 , nb_segment + i - seq_length);j < MIN((i < seq_length - 1 ? nb_segment - 1 : nb_segment) , i + 1);j++) {
        for (j = 0;j < MIN((i < seq_length - 1 ? nb_segment - 1 : nb_segment) , i + 1);j++) {
          forward[i][j] /= norm[i];
        }

        norm[i] = logl(norm[i]);
      }

      // computation of partial entropies

      segment_norm = 0.;
      for (j = i;j >= 0;j--) {
//        segment_norm += entropy_norm[j];
        segment_norm += norm[j];
        if (contrast[j] != D_INF) {
          normalized_contrast[j] = expl(contrast[j] - segment_norm);
        }
        else {
          normalized_contrast[j] = 0.;
        }
      }

      if (i < seq_length - 1) {
        for (j = 1;j < nb_segment;j++) {
          sequence_norm = expl(forward_norm[i] + backward_norm[i + 1] - likelihood[j]);

          for (k = MAX(nb_segment - 1 - j , nb_segment + i - seq_length);k <= MIN(nb_segment - 1 , i + nb_segment - 1 - j);k++) {
            if (k == nb_segment - 1 - j) {
              lsum = normalized_contrast[0] * backward[i + 1][k + 1] * sequence_norm;
              for (m = 0;m <= i;m++) {
                if (smoothed[j][m][0] > 0.) {
                  lbuff = lsum / smoothed[j][m][0];
                  if (lbuff > 0.) {
                    state_entropy[j][m][0] += lbuff * (backward_predicted_entropy[i + 1][k + 1] - logl(lbuff));
                  }
                }
              }
            }

            else {
              lsum = 0.;
              for (m = k + j - nb_segment + 1;m <= i;m++) {
                lsum += forward[m - 1][k + j - nb_segment] * normalized_contrast[m];
                if (smoothed[j][m][k + j - nb_segment + 1] > 0.) {
                  lbuff = lsum * backward[i + 1][k + 1] * sequence_norm /
                          smoothed[j][m][k + j - nb_segment + 1];
                  if (lbuff > 0.) {
                    state_entropy[j][m][k + j - nb_segment + 1] += lbuff * (backward_predicted_entropy[i + 1][k + 1] - logl(lbuff));
                  }
                }
              }
            }
          }
        }
      }
    }

#   ifdef DEBUG
    for (i = 1;i < nb_segment;i++) {
      cout << "\n" << i + 1 << " " << SEQ_label[SEQL_SEGMENTS] << endl;
      for (j = 0;j < seq_length;j++) {
        cout << j << " |";
        for (k = 0;k < nb_segment;k++) {
          cout << " " << state_entropy[i][j][k];
        }
        cout << endl;
      }
    }
#   endif

    // extraction of partial entropies for the different numbers of segments

    for (i = 1;i < nb_segment;i++) {
      for (j = 0;j < seq_length - 1;j++) {
        backward_partial_entropy[i][j + 1] = 0.;
        for (k = 0;k < nb_segment;k++) {
          if (state_entropy[i][j][k] < 0.) {
            state_entropy[i][j][k] = 0.;
          }
          if (smoothed[i][j][k] > 0.) {
            backward_partial_entropy[i][j + 1] += smoothed[i][j][k] * (state_entropy[i][j][k] -
                                                   log(smoothed[i][j][k]));
          }
        }
        if (backward_partial_entropy[i][j + 1] < 0.) {
          backward_partial_entropy[i][j + 1] = 0.;
        }
      }
    }

#   ifdef DEBUG
    cout << "\n" << SEQ_label[SEQL_SEGMENTATION_ENTROPY] << endl;
    for (i = 1;i < nb_segment - 1;i++) {
      cout << i + 1 << " " << SEQ_label[SEQL_SEGMENTS] << ": "
           << forward_predicted_entropy[seq_length - 1][i] << ", "
           << backward_predicted_entropy[0][nb_segment - 1 - i] << ", "
           << forward_partial_entropy[i][seq_length - 1] << ", "
           << backward_partial_entropy[i][1] << endl;
    }
    cout << nb_segment << " " << SEQ_label[SEQL_SEGMENTS] << ": "
         << forward_predicted_entropy[seq_length - 1][nb_segment - 1] << ", "
         << backward_predicted_entropy[0][0] << ", "
         << forward_partial_entropy[nb_segment - 1][seq_length - 1] << ", "
         << backward_partial_entropy[nb_segment - 1][1]
         << " | " << segmentation_entropy << endl;
#   endif

#   ifdef DEBUG
    for (i = 1;i < nb_segment;i++) {
      cout << "\n";
      for (j = 0;j < seq_length;j++) {
        cout << forward_partial_entropy[i][j] << " ";
      }
      if (i == nb_segment - 1) {
        cout << " | " << segmentation_entropy;
      }
      cout << endl;
    }

    for (i = 1;i < nb_segment;i++) {
      cout << "\n";
      if (i == nb_segment - 1) {
        cout << segmentation_entropy << " | ";
      }
      for (j = 1;j < seq_length;j++) {
        cout << backward_partial_entropy[i][j] << " ";
      }
      cout << endl;
    }
#   endif

    for (i = 1;i < nb_segment;i++) {
      for (j = seq_length - 1;j >= 1;j--) {
        forward_partial_entropy[i][j] -= forward_partial_entropy[i][j - 1];
        if (forward_partial_entropy[i][j] < 0.) {
          forward_partial_entropy[i][j] = 0.;
        }
      }
    }

    for (i = 1;i < nb_segment;i++) {
      backward_partial_entropy[i][0] = 0.;
      for (j = 1;j < seq_length - 1;j++) {
        backward_partial_entropy[i][j] -= backward_partial_entropy[i][j + 1];
        if (backward_partial_entropy[i][j] < 0.) {
          backward_partial_entropy[i][j] = 0.;
        }
      }
    }

#   ifdef MESSAGE
    for (i = 1;i < nb_segment;i++) {
      for (j = 0;j < seq_length;j++) {
        if (forward_partial_entropy[i][j] > change_point_entropy[i][j]) {
          cout << "\n" << SEQ_label[SEQL_BEGIN_CONDITIONAL_ENTROPY] << " ERROR: "
               << forward_partial_entropy[i][j] << " " << change_point_entropy[i][j]
               << " | " << i << ", " << j + 1 << endl;
        }

        if (backward_partial_entropy[i][j] > change_point_entropy[i][j]) {
          cout << "\n" << SEQ_label[SEQL_END_CONDITIONAL_ENTROPY] << " ERROR: "
               << backward_partial_entropy[i][j] << " " << change_point_entropy[i][j]
               << " | " << i << ", " << j + 1 << endl;
        }
      }
    }
#   endif

    if ((os) || (plot_set)) {
      switch (format) {

      case ASCII : {
        switch (output) {
        case CHANGE_POINT :
          *os << "\n" << SEQ_label[SEQL_POSTERIOR_CHANGE_POINT_PROBABILITY] << "\n\n";
          break;
        case SEGMENT :
          *os << "\n" << SEQ_label[SEQL_POSTERIOR_SEGMENT_PROBABILITY] << "\n\n";
          break;
        }

        profile_ascii_print(*os , index , nb_segment , backward_output ,
                            (output == CHANGE_POINT ? SEQ_label[SEQL_CHANGE_POINT] : SEQ_label[SEQL_SEGMENT]) ,
                            NULL , change_point , segment_length , prior_segment_length ,
                            forward_partial_entropy , backward_partial_entropy , change_point_entropy);

        *os << "\n" << SEQ_label[SEQL_POSSIBLE_SEGMENTATION_LIKELIHOOD] << ": " << rlikelihood << endl;
        *os << "\n" << SEQ_label[SEQL_SEGMENTATION_ENTROPY] << ": " << segmentation_entropy
            << "\n" << SEQ_label[SEQL_FIRST_ORDER_ENTROPY] << ": " << first_order_entropy
            << "\n" << SEQ_label[SEQL_CHANGE_POINT_ENTROPY] << ": " << change_point_entropy_sum
//            << " (" << change_point_entropy_sum / nb_segment << ")";
            << "\n" << SEQ_label[SEQL_MARGINAL_ENTROPY] << ": " << marginal_entropy << endl;

        // extraction of change-point uncertainty intervals

        if (output == CHANGE_POINT) {
          *os << "\n" << SEQ_label[SEQL_CHANGE_POINT_UNCERTAINTY_INTERVALS] << endl;
          for (i = 1;i < nb_segment;i++) {
            *os << SEQ_label[SEQL_CHANGE_POINT] << " " << i << " (";

            sum = 0.;
            j = 0;
            while (sum < CHANGE_POINT_UNCERTAINTY_PROBABILITY / 2) {
              j++;
              sum += backward_output[j][i];
            }
            *os << seq_index_parameter[j] << ", ";

#           ifdef MESSAGE
            while (sum <= 1. - CHANGE_POINT_UNCERTAINTY_PROBABILITY / 2) {
              j++;
              sum += backward_output[j][i];
            }
            *os << seq_index_parameter[j] << " | ";
#           endif

            sum = 0.;
            j = seq_length;
            while (sum < CHANGE_POINT_UNCERTAINTY_PROBABILITY / 2) {
              j--;
              sum += backward_output[j][i];
            }
            *os << seq_index_parameter[j] << ")" << endl;
          }
        }

        break;
      }

      case SPREADSHEET : {
        switch (output) {
        case CHANGE_POINT :
          *os << "\n" << SEQ_label[SEQL_POSTERIOR_CHANGE_POINT_PROBABILITY] << "\n\n";
          break;
        case SEGMENT :
          *os << "\n" << SEQ_label[SEQL_POSTERIOR_SEGMENT_PROBABILITY] << "\n\n";
          break;
        }

        profile_spreadsheet_print(*os , index , nb_segment , backward_output ,
                                  (output == CHANGE_POINT ? SEQ_label[SEQL_CHANGE_POINT] : SEQ_label[SEQL_SEGMENT]) ,
                                  common_contrast , NULL , change_point , segment_length , prior_segment_length ,
                                  forward_partial_entropy , backward_partial_entropy , change_point_entropy);

        *os << "\n" << SEQ_label[SEQL_POSSIBLE_SEGMENTATION_LIKELIHOOD] << "\t" << rlikelihood << endl;
        *os << "\n" << SEQ_label[SEQL_SEGMENTATION_ENTROPY] << "\t" << segmentation_entropy
            << "\n" << SEQ_label[SEQL_FIRST_ORDER_ENTROPY] << "\t" << first_order_entropy
            << "\n" << SEQ_label[SEQL_CHANGE_POINT_ENTROPY] << "\t" << change_point_entropy_sum
//            << "\t" << change_point_entropy_sum / nb_segment;
            << "\n" << SEQ_label[SEQL_MARGINAL_ENTROPY] << "\t" << marginal_entropy << endl;
        break;
      }

      case GNUPLOT : {
        profile_plot_print(*os , index , nb_segment , backward_output ,
                           common_contrast , NULL , change_point ,
                           segment_length , prior_segment_length , forward_partial_entropy ,
                           backward_partial_entropy , change_point_entropy);
        break;
      }

      case PLOT : {
        MultiPlotSet &plot = *plot_set;

        i = 1;
        for (j = 1;j < nb_variable;j++) {
          if ((model_type[j - 1] == POISSON_CHANGE) ||
              (model_type[j - 1] == NEGATIVE_BINOMIAL_0_CHANGE) || (model_type[j - 1] == NEGATIVE_BINOMIAL_1_CHANGE) ||
              (model_type[j - 1] == GAUSSIAN_CHANGE) ||
              (model_type[j - 1] == VARIANCE_CHANGE) || (model_type[j - 1] == BAYESIAN_POISSON_CHANGE) ||
              (model_type[j - 1] == BAYESIAN_GAUSSIAN_CHANGE)) {
            i++;
          }
        }

        profile_plotable_write(plot[i] , index , nb_segment , backward_output);
        i++;
        change_point_profile_plotable_write(plot[i] , index , nb_segment , change_point);
        i++;
        segment_length_distribution_plotable_write(plot[i] , nb_segment , segment_length_max ,
                                                   segment_length , prior_segment_length);
        i++;
        change_point_profile_plotable_write(plot[i] , index , nb_segment , forward_partial_entropy);
        i++;
        change_point_profile_plotable_write(plot[i] , index , nb_segment , backward_partial_entropy);
        i++;
        entropy_profile_plotable_write(plot[i] , index , forward_partial_entropy[nb_segment - 1] ,
                                       backward_partial_entropy[nb_segment - 1] ,
                                       change_point_entropy[nb_segment - 1]);
        break;
      }
      }
    }
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
  delete [] normalized_contrast;

  for (i = 0;i < seq_length;i++) {
    delete [] forward[i];
  }
  delete [] forward;

  delete [] segment_predicted;

  for (i = 0;i < seq_length;i++) {
    delete [] forward_predicted_entropy[i];
  }
  delete [] forward_predicted_entropy;

  delete [] norm;
  delete [] forward_norm;
//  delete [] entropy_norm;

  delete [] likelihood;

  for (i = 0;i < seq_length;i++) {
    delete [] backward[i];
  }
  delete [] backward;

  for (i = 0;i < seq_length;i++) {
    delete [] backward_predicted_entropy[i];
  }
  delete [] backward_predicted_entropy;

  delete [] backward_norm;

  for (i = 1;i < nb_segment;i++) {
    for (j = 0;j < seq_length;j++) {
      delete [] smoothed[i][j];
    }
    delete [] smoothed[i];
  }
  delete [] smoothed;

  for (i = 0;i < seq_length;i++) {
    delete [] backward_output[i];
  }
  delete [] backward_output;

  for (i = 0;i < nb_segment;i++) {
    delete segment_length[i];
  }
  delete [] segment_length;

  delete prior_segment_length;

  for (i = 1;i < nb_segment;i++) {
    delete [] change_point[i];
  }
  delete [] change_point;

  delete [] entropy_smoothed;

  for (i = 1;i < nb_segment;i++) {
    for (j = 0;j < seq_length;j++) {
      delete [] state_entropy[i][j];
    }
    delete [] state_entropy[i];
  }
  delete [] state_entropy;

  for (i = 1;i < nb_segment;i++) {
    delete [] forward_partial_entropy[i];
  }
  delete [] forward_partial_entropy;

  for (i = 1;i < nb_segment;i++) {
    delete [] backward_partial_entropy[i];
  }
  delete [] backward_partial_entropy;

  for (i = 1;i < nb_segment;i++) {
    delete [] change_point_entropy[i];
  }
  delete [] change_point_entropy;

  return rlikelihood;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Simulation of segmentations of a single sequence  or a sample of sequences.
 *
 *  \param[in] index           sequence index,
 *  \param[in] nb_segment      number of segments,
 *  \param[in] model_type      segment model types,
 *  \param[in] common_contrast flag contrast functions common to the individuals,
 *  \param[in] shape_parameter negative binomial shape parameters,
 *  \param[in] rank            ranks (for ordinal variables),
 *  \param[in] os              stream,
 *  \param[in] format          file format (ASCII/SPREADSHEET),
 *  \param[in] nb_segmentation number of segmentations.
 *
 *  \return                    log-likelihood of the multiple change-point model.
 */
/*--------------------------------------------------------------*/

double Sequences::forward_backward_sampling(int index , int nb_segment , segment_model *model_type ,
                                            bool common_contrast , double *shape_parameter ,
                                            double **rank , ostream &os , output_format format ,
                                            int nb_segmentation) const

{
  int i , j , k , m;
  int seq_length , segment_length , *inf_bound_parameter , *seq_index_parameter , *change_point ,
      *psegment;
  double sum , likelihood , segmentation_likelihood , *backward , *cumul_backward , **seq_mean ,
         **hyperparam , ***factorial , ***binomial_coeff , ***mean , ***variance ,
         ***intercept , ***slope , ***autoregressive_coeff;
  long double segment_norm , sequence_norm , *contrast , *norm , **forward;


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

  forward = new long double*[seq_length];
  for (i = 0;i < seq_length;i++) {
    forward[i] = new long double[nb_segment];
  }

  norm = new long double[seq_length];

  backward = new double[seq_length];
  cumul_backward = new double[seq_length];

  change_point = new int[nb_segment + 1];

  mean = new double**[nb_variable];
  variance = new double**[nb_variable];
  intercept = new double**[nb_variable];
  slope = new double**[nb_variable];
  autoregressive_coeff = new double**[nb_variable];

  for (i = 1;i < nb_variable;i++) {
    if ((model_type[i - 1] == POISSON_CHANGE) || (model_type[i - 1] == NEGATIVE_BINOMIAL_0_CHANGE) ||
        (model_type[i - 1] == NEGATIVE_BINOMIAL_1_CHANGE) || (model_type[i - 1] == GAUSSIAN_CHANGE) ||
        (model_type[i - 1] == VARIANCE_CHANGE) || (model_type[i - 1] == BAYESIAN_POISSON_CHANGE) ||
        (model_type[i - 1] == BAYESIAN_GAUSSIAN_CHANGE)) {
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

    else if (model_type[i - 1] == LINEAR_MODEL_CHANGE) {
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

# ifdef DEBUG
  double **segment_probability;

  segment_probability = new double*[seq_length];
  for (i = 0;i < seq_length;i++) {
    segment_probability[i] = new double[nb_segment];
    for (j = 0;j < nb_segment;j++) {
      segment_probability[i][j] = 0.;
    }
  }
# endif

  // forward recurrence

  for (i = 0;i < seq_length;i++) {

    // computation of segment log-likelihoods

    forward_contrast(i , index , model_type , common_contrast , factorial ,
                     shape_parameter , binomial_coeff , seq_mean , seq_index_parameter ,
                     hyperparam , rank , contrast);

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
    }
    norm[i] = 0.;

    for (j = MAX(0 , nb_segment + i - seq_length);j < MIN((i < seq_length - 1 ? nb_segment - 1 : nb_segment) , i + 1);j++) {
      if (j == 0) {
        forward[i][j] = contrast[0];
      }
      else {
        for (k = i;k >= j;k--) {
          forward[i][j] += contrast[k] * forward[k - 1][j - 1];
        }
      }

      norm[i] += forward[i][j];
    }

    if (norm[i] > 0.) {
      for (j = MAX(0 , nb_segment + i - seq_length);j < MIN((i < seq_length - 1 ? nb_segment - 1 : nb_segment) , i + 1);j++) {
        forward[i][j] /= norm[i];
      }

      norm[i] = logl(norm[i]);
    }

#   ifdef DEBUG
    cout << i << " |";
    for (j = 0;j < nb_segment;j++) {
      cout << " " << forward[i][j];
    }
    cout << " | " << expl(norm[i]) << endl;
#   endif

  }

  sequence_norm = segment_norm + norm[seq_length - 1];

  if (forward[seq_length - 1][nb_segment - 1] > 0.) {
    likelihood = logl(forward[seq_length - 1][nb_segment - 1]) + sequence_norm;
  }
  else {
    likelihood = D_INF;
  }

  if (likelihood != D_INF) {

    // backward pass

#   ifdef MESSAGE
    cout << "\n";
#   endif

    for (i = 0;i < nb_segmentation;i++) {
      j = seq_length - 1;
      change_point[nb_segment] = seq_length;
      psegment = int_sequence[index == I_DEFAULT ? 0 : index][0] + j;
      segmentation_likelihood = sequence_norm;

      for (k = nb_segment - 1;k >= 0;k--) {

        // computation of segment log-likelihoods

        forward_contrast(j , index , model_type , common_contrast , factorial ,
                         shape_parameter , binomial_coeff , seq_mean , seq_index_parameter ,
                         hyperparam , rank , contrast , k);

        segment_norm = 0.;
        for (m = j;m >= k;m--) {
          segment_norm += norm[m];
          if (contrast[m] != D_INF) {
            contrast[m] = expl(contrast[m] - segment_norm);
          }
          else {
            contrast[m] = 0.;
          }
        }

        if (k > 0) {
          for (m = j;m >= k;m--) {
            backward[m] = contrast[m] * forward[m - 1][k - 1] / forward[j][k];
          }
          stat_tool::cumul_computation(j - k , backward + k , cumul_backward);
          segment_length = j - (k + cumul_method(j - k , cumul_backward)) + 1;

#         ifdef MESSAGE
          sum = 0.;
          for (m = j;m >= k;m--) {
            sum += backward[m];
          }
          if ((sum < 1. - DOUBLE_ERROR) || (sum > 1. + DOUBLE_ERROR)) {
            cout << "\nERROR: " << j << " " << sum << endl;
          }
#         endif

        }

        else {
          segment_length = j + 1;
        }

        segmentation_likelihood += logl(contrast[j - segment_length + 1]);

        for (m = j;m > j - segment_length;m--) {
          *psegment-- = k;
        }
        j -= segment_length;
        change_point[k] = j + 1;
      }

      for (j = 1;j < nb_variable;j++) {
        piecewise_linear_function(index , j , nb_segment , model_type[j - 1] , common_contrast ,
                                  change_point , seq_index_parameter , NULL , mean[j] , variance[j] ,
                                  NULL , intercept[j] , slope[j] , autoregressive_coeff[j]);
      }

#     ifdef DEBUG

      // approximation of smoothed probabilities

      psegment = int_sequence[index == I_DEFAULT ? 0 : index][0];
      for (j = 0;j < seq_length;j++) {
        segment_probability[j][*psegment++]++;
      }
#     endif

#     ifdef MESSAGE
      if (i == 0) {
        os << "\n";
      }

      switch (format) {

      case ASCII : {
        psegment = int_sequence[index == I_DEFAULT ? 0 : index][0];
        for (j = 0;j < seq_length;j++) {
          os << *psegment++ << " ";
        }

        os << "  " << i + 1 << "  " << segmentation_likelihood << "   ("
           << exp(segmentation_likelihood - likelihood) << ")" << endl;

        os << (nb_segment == 2 ? SEQ_label[SEQL_CHANGE_POINT] : SEQ_label[SEQL_CHANGE_POINTS]) << ": ";

        for (j = 1;j < nb_segment;j++) {
          os << seq_index_parameter[change_point[j]];
          if (j < nb_segment - 1) {
            os << ", ";
          }
        }
        os << endl;

        for (j = 1;j < nb_variable;j++) {
          piecewise_linear_function_ascii_print(os , index , j , nb_segment , model_type[j - 1] ,
                                                common_contrast , change_point , seq_index_parameter ,
                                                mean[j] , variance[j] , intercept[j] , slope[j] ,
                                                autoregressive_coeff[j]);
        }
        break;
      }

      case SPREADSHEET : {
        psegment = int_sequence[index == I_DEFAULT ? 0 : index][0];
        for (j = 0;j < seq_length;j++) {
          os << *psegment++ << "\t";
        }

        os << "\t" << i + 1 << "\t" << segmentation_likelihood  << "\t"
           << exp(segmentation_likelihood - likelihood) << endl;

        for (j = 1;j < nb_variable;j++) {
          piecewise_linear_function_spreadsheet_print(os , index , j , nb_segment , model_type[j - 1] ,
                                                      common_contrast , change_point , seq_index_parameter ,
                                                      mean[j] , variance[j] , intercept[j] , slope[j] ,
                                                      autoregressive_coeff[j]);
        }
        break;
      }
      }
#     endif

    }

#   ifdef DEBUG
    if (nb_segmentation >= 1000) {
      for (i = 0;i < seq_length;i++) {
        for (j = 0;j < nb_segment;j++) {
          segment_probability[i][j] /= nb_segmentation;
        }
      }

      psegment = int_sequence[index == I_DEFAULT ? 0 : index][0];
      for (i = 0;i < seq_length;i++) {
        *psegment++ = I_DEFAULT;
      }

      os << "\n" << SEQ_label[SEQL_POSTERIOR_SEGMENT_PROBABILITY] << "\n\n";

      profile_ascii_print(os , index , nb_segment , segment_probability ,
                          SEQ_label[SEQL_SEGMENT]);
    }
#   endif

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

  delete [] norm;

  delete [] backward;
  delete [] cumul_backward;

  delete [] change_point;

  for (i = 1;i < nb_variable;i++) {
    if ((model_type[i - 1] == POISSON_CHANGE) || (model_type[i - 1] == NEGATIVE_BINOMIAL_0_CHANGE) ||
        (model_type[i - 1] == NEGATIVE_BINOMIAL_1_CHANGE) || (model_type[i - 1] == GAUSSIAN_CHANGE) ||
        (model_type[i - 1] == VARIANCE_CHANGE) || (model_type[i - 1] == BAYESIAN_POISSON_CHANGE) ||
        (model_type[i - 1] == BAYESIAN_GAUSSIAN_CHANGE)) {
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

    else if (model_type[i - 1] == LINEAR_MODEL_CHANGE) {
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

# ifdef DEBUG
  for (i = 0;i < seq_length;i++) {
    delete [] segment_probability[i];
  }
  delete [] segment_probability;
# endif

  return likelihood;
}


};  // namespace sequence_analysis
