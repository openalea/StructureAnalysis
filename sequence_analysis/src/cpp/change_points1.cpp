/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       V-Plants: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2010 CIRAD/INRIA Virtual Plants
 *
 *       File author(s): Y. Guedon (yann.guedon@cirad.fr)
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

#include "tool/config.h"

#include "stat_tool/stat_tools.h"
#include "stat_tool/curves.h"
#include "stat_tool/markovian.h"

#include "sequences.h"
#include "sequence_label.h"

using namespace std;


extern int column_width(int value);
extern int column_width(int nb_value , const double *value , double scale = 1.);
extern int column_width(int nb_value , const long double *value);



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
  int nb_parameter , *psegment , *pisequence;


  used_output = new bool[NB_OUTPUT];

  nb_parameter = nb_segment - 1;
//  nb_parameter = 0;

  for (i = 1;i < nb_variable;i++) {
    if (model_type[i - 1] == MULTINOMIAL_CHANGE) {
      psegment = int_sequence[index][0] + 1;
      pisequence = int_sequence[index][i];

      for (j = 0;j < marginal_distribution[i]->nb_value;j++) {
        used_output[j] = false;
      }
      used_output[*pisequence++] = true;

      for (j = 1;j < length[index];j++) {
        if (*psegment != *(psegment - 1)) {
          for (k = 0;k < marginal_distribution[i]->nb_value;k++) {
            used_output[k] = false;
          }
          used_output[*pisequence] = true;
        }

        else if (!used_output[*pisequence]) {
          nb_parameter++;
          used_output[*pisequence] = true;
        }

        psegment++;
        pisequence++;
      }
    }

    else if ((model_type[i - 1] == POISSON_CHANGE) ||
             (model_type[i - 1] == MEAN_CHANGE) ||
             (model_type[i - 1] == MEAN_VARIANCE_CHANGE)) {
      nb_parameter += nb_segment;
    }

    else if (model_type[i - 1] == VARIANCE_CHANGE) {
      nb_parameter += nb_segment + 1;
    }

    else {
      nb_parameter += 2 * nb_segment;
    }
  }

  if (model_type[0] == MEAN_CHANGE) {
    nb_parameter++;
  }
  else if (model_type[0] == MEAN_VARIANCE_CHANGE) {
    nb_parameter += nb_segment;
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
  int max_nb_value , *frequency , *pisequence , *psegment;
  double sum , factorial_sum , diff , likelihood , *prsequence;
  long double residual;


  max_nb_value = 0;
  for (i = 1;i < nb_variable;i++) {
    if (model_type[i - 1] == MULTINOMIAL_CHANGE) {
      if (marginal_distribution[i]->nb_value > max_nb_value) {
        max_nb_value = marginal_distribution[i]->nb_value;
      }
    }
  }

  if (max_nb_value > 0) {
    frequency = new int[max_nb_value];
  }
  else {
    frequency = NULL;
  }

  if ((model_type[0] == MEAN_CHANGE) || (model_type[0] == MEAN_VARIANCE_CHANGE)) {
    residual = 0.;
  }
  else {
    likelihood = 0.;
  }

  for (i = 1;i < nb_variable;i++) {
    if (model_type[i - 1] == MULTINOMIAL_CHANGE) {
      for (j = 0;j < marginal_distribution[i]->nb_value;j++) {
        frequency[j] = 0;
      }

      pisequence = int_sequence[index][i];
      for (j = 0;j < length[index];j++) {
        frequency[*pisequence++]++;
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

      pisequence = int_sequence[index][i];
      for (j = 0;j < length[index];j++) {
        sum += *pisequence++;
        for (k = 2;k <= int_sequence[index][i][j];k++) {
          factorial_sum += log((double)k);
        }
      }

      if (sum > 0.) {
        likelihood += sum * (log(sum / (length[index])) - 1) - factorial_sum;
      }
    }

    else {
      if ((model_type[i - 1] != MEAN_CHANGE) &&
          (model_type[i - 1] != MEAN_VARIANCE_CHANGE)) {
        residual = 0.;
      }

      if (model_type[i - 1] == ORDINAL_GAUSSIAN_CHANGE) {
        pisequence = int_sequence[index][i];
        sum = rank[i][*pisequence++];

        for (j = 1;j < length[index];j++) {
          diff = rank[i][*pisequence] - sum / j;
          residual += ((double)j / (double)(j + 1)) * diff * diff;
          sum += rank[i][*pisequence++];
        }
      }

      else {
        if (type[i] != REAL_VALUE) {
          pisequence = int_sequence[index][i];
          sum = *pisequence++;

          for (j = 1;j < length[index];j++) {
            diff = *pisequence - sum / j;
            residual += ((double)j / (double)(j + 1)) * diff * diff;
            sum += *pisequence++;
          }

#         ifdef MESSAGE
          if ((model_type[i - 1] != MEAN_CHANGE) &&
              (model_type[i - 1] != MEAN_VARIANCE_CHANGE)) {
            double mean , diff;
            long double residual2;


            mean = 0.;
            pisequence = int_sequence[index][i];
            for (j = 0;j < length[index];j++) {
              mean += *pisequence++;
            }
            mean /= length[index];

            residual2 = 0.;
            pisequence = int_sequence[index][i];
            for (j = 0;j < length[index];j++) {
              diff = *pisequence++ - mean;
              residual2 += diff * diff;
            }

            if ((residual < residual2 - DOUBLE_ERROR) || (residual > residual2 + DOUBLE_ERROR)) {
              cout << "\nERROR: " << residual << " " << residual2 << endl;
            }
          }
#         endif

        }

        else {
          prsequence = real_sequence[index][i];
          sum = *prsequence++;

          for (j = 1;j < length[index];j++) {
            diff = *prsequence - sum / j;
            residual += ((double)j / (double)(j + 1)) * diff * diff;
            sum += *prsequence++;
          }
        }
      }

      if ((model_type[i - 1] != MEAN_CHANGE) &&
          (model_type[i - 1] != MEAN_VARIANCE_CHANGE)) {
//        if (residual > 0.) {
        if (residual > sqrt((double)length[index]) * ROUNDOFF_ERROR) {
          likelihood -= ((double)length[index] / 2.) * (logl(residual / length[index]) +
                         log(2 * M_PI) + 1);
/*          likelihood -= ((double)length[index] / 2.) * (logl(residual / (length[index] - 1)) +
                         log(2 * M_PI)) - (double)(length[index] - 1) / 2.; */
        }
        else {
          likelihood = D_INF;
          break;
        }
      }
    }
  }

  if ((model_type[0] == MEAN_CHANGE) || (model_type[0] == MEAN_VARIANCE_CHANGE)) {
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

  psegment = int_sequence[index][0];
  for (i = 0;i < length[index];i++) {
    *psegment++ = 0;
  }

  delete [] frequency;

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
  bool *segment_mean;
  register int i , j , k , m , n;
  int max_nb_segment , *change_point , *psegment , *pisequence;
  double diff , global_variance , *prsequence , *segment_variance , **mean , **variance;
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
    if ((model_type[i - 1] == POISSON_CHANGE) || (model_type[i - 1] == GAUSSIAN_CHANGE) ||
        (model_type[i - 1] == MEAN_CHANGE) || (model_type[i - 1] == VARIANCE_CHANGE) ||
        (model_type[i - 1] == MEAN_VARIANCE_CHANGE)) {
      mean[i] = new double[max_nb_segment];
    }
    else {
      mean[i] = NULL;
    }
  }

  if (nb_sequence == 1) {
    variance = new double*[nb_variable];
    for (i = 1;i < nb_variable;i++) {
      if ((model_type[i - 1] == POISSON_CHANGE) || (model_type[i - 1] == GAUSSIAN_CHANGE) ||
          (model_type[i - 1] == MEAN_CHANGE) || (model_type[i - 1] == VARIANCE_CHANGE) ||
          (model_type[i - 1] == MEAN_VARIANCE_CHANGE)) {
        variance[i] = new double[nb_segment[0]];
      }
      else {
        variance[i] = NULL;
      }
    }

    if (model_type[0] == MEAN_CHANGE) {
      global_variance = 0.;
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
    segment_mean = new bool[nb_variable];

    segment_mean[0] = false;
    for (i = 1;i < nb_variable;i++) {
      if ((model_type[i - 1] == POISSON_CHANGE) || (model_type[i - 1] == GAUSSIAN_CHANGE) ||
          (model_type[i - 1] == MEAN_CHANGE) || (model_type[i - 1] == MEAN_VARIANCE_CHANGE)) {
        segment_mean[i] = true;
      }
      else {
        segment_mean[i] = false;
      }
    }

    seq = new Sequences(*this , segment_mean);
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
      if ((model_type[j - 1] == POISSON_CHANGE) || (model_type[j - 1] == GAUSSIAN_CHANGE) ||
          (model_type[j - 1] == MEAN_CHANGE) || (model_type[j - 1] == VARIANCE_CHANGE) ||
          (model_type[j - 1] == MEAN_VARIANCE_CHANGE)) {
        if (type[j] != REAL_VALUE) {
          pisequence = int_sequence[i][j];
          for (k = 0;k < nb_segment[i];k++) {
            mean[j][k] = 0.;
            for (m = change_point[k];m <  change_point[k + 1];m++) {
              mean[j][k] += *pisequence++;
            }
            mean[j][k] /= (change_point[k + 1] - change_point[k]);
          }
        }

        else {
          prsequence = real_sequence[i][j];
          for (k = 0;k < nb_segment[i];k++) {
            mean[j][k] = 0.;
            for (m = change_point[k];m <  change_point[k + 1];m++) {
              mean[j][k] += *prsequence++;
            }
            mean[j][k] /= (change_point[k + 1] - change_point[k]);
          }
        }

        if (nb_sequence == 1) {
          if (type[j] != REAL_VALUE) {
            pisequence = int_sequence[i][j];
            for (k = 0;k < nb_segment[i];k++) {
              variance[j][k] = 0.;
              for (m = change_point[k];m <  change_point[k + 1];m++) {
                diff = *pisequence++ - mean[j][k];
                variance[j][k] += diff * diff;
              }

              if (model_type[j - 1] == MEAN_CHANGE) {
                global_variance += variance[j][k];
              }
              if (model_type[j - 1] == MEAN_VARIANCE_CHANGE) {
                segment_variance[k] += variance[j][k];
              }

              variance[j][k] /= (change_point[k + 1] - change_point[k]);
//              variance[j][k] /= (change_point[k + 1] - change_point[k] - 1);
            }
          }

          else {
            prsequence = real_sequence[i][j];
            for (k = 0;k < nb_segment[i];k++) {
              variance[j][k] = 0.;
              for (m = change_point[k];m <  change_point[k + 1];m++) {
                diff = *prsequence++ - mean[j][k];
                variance[j][k] += diff * diff;
              }

              if (model_type[j - 1] == MEAN_CHANGE) {
                global_variance += variance[j][k];
              }
              if (model_type[j - 1] == MEAN_VARIANCE_CHANGE) {
                segment_variance[k] += variance[j][k];
              }

              variance[j][k] /= (change_point[k + 1] - change_point[k]);
//              variance[j][k] /= (change_point[k + 1] - change_point[k] - 1);
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
        os << endl;

        if (model_type[0] == MEAN_CHANGE) {
          os << SEQ_label[SEQL_GLOBAL_STANDARD_DEVIATION] << ": "
             << sqrt(global_variance / ((nb_variable - 1) * length[i])) << endl;
        }

        else if (model_type[0] == MEAN_VARIANCE_CHANGE) {
          os << SEQ_label[SEQL_SEGMENT_STANDARD_DEVIATION] << ": ";
          for (j = 0;j < nb_segment[i];j++) {
            os << sqrt(segment_variance[j] / ((nb_variable - 1) *
                        (change_point[j + 1] - change_point[j]))) << ", ";
          }
          os << endl;
        }
      }

      if (nb_variable > 2) {
        os << "\n";
      }
      for (j = 1;j < nb_variable;j++) {
        if ((model_type[j - 1] == GAUSSIAN_CHANGE) || (model_type[j - 1] == MEAN_CHANGE) ||
            (model_type[j - 1] == VARIANCE_CHANGE) || (model_type[j - 1] == MEAN_VARIANCE_CHANGE)) {
          if (nb_variable > 2) {
            os << STAT_label[STATL_VARIABLE] << " " << j << "   ";
          }
          os << SEQ_label[SEQL_SEGMENT_MEAN] << "  " << SEQ_label[SEQL_SEGMENT_STANDARD_DEVIATION] << ": ";
          for (k = 0;k < nb_segment[i];k++) {
            os << mean[j][k] << " " << sqrt(variance[j][k]) << " | ";
          }
          os << endl;
        }

        else if (model_type[j - 1] == POISSON_CHANGE) {
          if (nb_variable > 2) {
            os << STAT_label[STATL_VARIABLE] << " " << j << "   ";
          }
          os << SEQ_label[SEQL_SEGMENT_MEAN] << "  " << SEQ_label[SEQL_SEGMENT_VARIANCE] << ": ";
          for (k = 0;k < nb_segment[i];k++) {
            os << mean[j][k] << " " << variance[j][k] << " | ";
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
        if (segment_mean[k]) {
          prsequence = seq->real_sequence[i][j++];
          for (m = 0;m < nb_segment[i];m++) {
            for (n = change_point[m];n < change_point[m + 1];n++) {
              *prsequence++ = mean[k][m];
            }
          }
        }
      }
      break;
    }

    // calcul des residus

    case SUBTRACTION_RESIDUAL : {
      for (j = 1;j < nb_variable;j++) {
        if (type[j] != REAL_VALUE) {
          real_sequence[i][j] = new double[length[i]];

          prsequence = real_sequence[i][j];
          pisequence = int_sequence[i][j];
          for (k = 0;k < nb_segment[i];k++) {
            for (m = change_point[k];m < change_point[k + 1];m++) {
              *prsequence++ = *pisequence++ - mean[j][k];
            }
          }

          delete [] int_sequence[i][j];
          int_sequence[i][j] = NULL;
        }

        else {
          prsequence = real_sequence[i][j];
          for (k = 0;k < nb_segment[i];k++) {
            for (m = change_point[k];m < change_point[k + 1];m++) {
              *prsequence++ -= mean[j][k];
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

          prsequence = real_sequence[i][j];
          pisequence = int_sequence[i][j];
          for (k = 0;k < nb_segment[i];k++) {
            if (mean[j][k] != 0.) {
              for (m = change_point[k];m < change_point[k + 1];m++) {
                *prsequence++ = *pisequence++ / mean[j][k];
              }
            }
            else {
              for (m = change_point[k];m < change_point[k + 1];m++) {
                *prsequence++ = *pisequence++;
              }
            }
          }

          delete [] int_sequence[i][j];
          int_sequence[i][j] = NULL;
        }

        else {
          prsequence = real_sequence[i][j];
          for (k = 0;k < nb_segment[i];k++) {
            if (mean[j][k] != 0.) {
              for (m = change_point[k];m < change_point[k + 1];m++) {
                *prsequence++ /= mean[j][k];
              }
            }
            else {
              prsequence += (change_point[k + 1] - change_point[k]);
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
    }
    delete [] variance;

    delete [] segment_variance;
  }

  if (output == SEQUENCE) {
    delete [] segment_mean;
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
  int index , max_nb_value , nb_parameter , *change_point = NULL , *pindex_param ,
      *frequency , *pisequence , *psegment , inb_segment[1];
  double sum , factorial_sum , diff , segmentation_likelihood , segment_penalty ,
         penalized_likelihood , *mean , *prsequence , *segment_variance , **rank;
  long double residual;
  Sequences *iseq , *seq , *oseq;


  oseq = NULL;
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

  for (i = 0;i < nb_variable;i++) {
    if ((model_type[i] == MULTINOMIAL_CHANGE) || (model_type[i] == ORDINAL_GAUSSIAN_CHANGE) ||
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
                        << SEQ_error[SEQR_POSITIVE_MIN_VALUE];
          error.update((error_message.str()).c_str());
        }

        if (!marginal_distribution[i]) {
          status = false;
          ostringstream error_message;
          error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                        << STAT_error[STATR_MARGINAL_FREQUENCY_DISTRIBUTION];
          error.update((error_message.str()).c_str());
        }

        else if (model_type[i] == MULTINOMIAL_CHANGE) {
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

      if (((model_type[i] == MULTINOMIAL_CHANGE) || (model_type[i] == ORDINAL_GAUSSIAN_CHANGE)) &&
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

        pindex_param = index_parameter[index] + 1;
        i = 1;
        for (j = 1;j < length[index];j++) {
          if (*pindex_param++ == change_point[i]) {
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
      if ((model_type[i] == MULTINOMIAL_CHANGE) &&
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
          pisequence = int_sequence[index][i];
          for (j = 0;j < length[index];j++) {
            mean[i] += *pisequence++;
          }
        }

        else {
          prsequence = real_sequence[index][i];
          for (j = 0;j < length[index];j++) {
            mean[i] += *prsequence++;
          }
        }

        mean[i] /= length[index];
      }
    }

    if (model_type[0] == MEAN_CHANGE) {
      residual = 0.;
    }
    else {
      segmentation_likelihood = 0.;
    }

    for (i = 0;i < nb_segment;i++) {
      if (model_type[0] == MEAN_VARIANCE_CHANGE) {
        residual = 0.;
      }

      for (j = 0;j < nb_variable;j++) {

        // calcul des log-vraisemblances des segments

        if (model_type[j] == MULTINOMIAL_CHANGE) {

          for (k = 0;k < marginal_distribution[j]->nb_value;k++) {
            frequency[k] = 0;
          }

          pisequence = int_sequence[index][j] + change_point[i];
          for (k = change_point[i];k < change_point[i + 1];k++) {
            frequency[*pisequence++]++;
          }

          for (k = 0;k < marginal_distribution[j]->nb_value;k++) {
            if (frequency[k] > 0) {
              segmentation_likelihood += frequency[k] * log((double)frequency[k] /
                                          (double)(change_point[i + 1] - change_point[i]));
            }
          }
        }

        else if (model_type[j] == POISSON_CHANGE) {
          sum = 0.;
          factorial_sum = 0.;

          pisequence = int_sequence[index][j] + change_point[i];
          for (k = change_point[i];k < change_point[i + 1];k++) {
            sum += *pisequence++;
            for (m = 2;m <= int_sequence[index][j][k];m++) {
              factorial_sum += log((double)m);
            }
          }

          if (sum > 0.) {
            segmentation_likelihood += sum * (log(sum / (change_point[i + 1] - change_point[i])) - 1) -
                                       factorial_sum;
          }
          else {
            segmentation_likelihood = D_INF;
            break;
          }
        }

        else {
          if (model_type[j] == VARIANCE_CHANGE) {
            residual = 0.;

            if (type[j] != REAL_VALUE) {
              pisequence = int_sequence[index][j] + change_point[i];
              for (k = change_point[i];k < change_point[i + 1];k++) {
                diff = *pisequence++ - mean[j];
                residual += diff * diff;
              }
            }

            else {
              prsequence = real_sequence[index][j] + change_point[i];
              for (k = change_point[i];k < change_point[i + 1];k++) {
                diff = *prsequence++ - mean[j];
                residual += diff * diff;
              }
            }
          }

          else if (model_type[j] == ORDINAL_GAUSSIAN_CHANGE) {
            pisequence = int_sequence[index][j] + change_point[i];
            residual = 0.;
            sum = rank[j][*pisequence++];

            for (k = change_point[i] + 1;k < change_point[i + 1];k++) {
              diff = rank[j][*pisequence] - sum / (k - change_point[i]);
              residual += ((double)(k - change_point[i]) / (double)(k - change_point[i] + 1)) *
                          diff * diff;
              sum += rank[j][*pisequence++];
            }
          }

          else {
            if (model_type[j] == GAUSSIAN_CHANGE) {
              residual = 0.;
            }

            if (type[j] != REAL_VALUE) {
              pisequence = int_sequence[index][j] + change_point[i];
              sum = *pisequence++;

              for (k = change_point[i] + 1;k < change_point[i + 1];k++) {
                diff = *pisequence - sum / (k - change_point[i]);
                residual += ((double)(k - change_point[i]) / (double)(k - change_point[i] + 1)) *
                            diff * diff;
                sum += *pisequence++;
              }
            }

            else {
              prsequence = real_sequence[index][j] + change_point[i];
              sum = *prsequence++;

              for (k = change_point[i] + 1;k < change_point[i + 1];k++) {
                diff = *prsequence - sum / (k - change_point[i]);
                residual += ((double)(k - change_point[i]) / (double)(k - change_point[i] + 1)) *
                            diff * diff;
                sum += *prsequence++;
              }
            }
          }

          if ((model_type[j] != MEAN_CHANGE) && (model_type[j] != MEAN_VARIANCE_CHANGE)) {
//            if (residual > 0.) {
            if (residual > sqrt((double)(change_point[i + 1] - change_point[i])) * ROUNDOFF_ERROR) {
              segmentation_likelihood -= ((double)(change_point[i + 1] - change_point[i]) / 2.) * (logl(residual /
                                           (change_point[i + 1] - change_point[i])) + log(2 * M_PI) + 1);
/*              segmentation_likelihood -= ((double)(change_point[i + 1] - change_point[i]) / 2.) * (logl(residual /
                                           (change_point[i + 1] - change_point[i])) + log(2 * M_PI)) + (double)(change_point[i + 1] - change_point[i]) / 2.; */
            }
            else {
              segmentation_likelihood = D_INF;
              break;
            }
          }
        }
      }

      if (model_type[0] == MEAN_VARIANCE_CHANGE) {
//        if (residual > 0.) {
        if (residual > sqrt((double)(nb_variable * (change_point[i + 1] - change_point[i]))) * ROUNDOFF_ERROR) {
          segmentation_likelihood -= ((double)(nb_variable * (change_point[i + 1] - change_point[i])) / 2.) *
                                     (logl(residual / (nb_variable * (change_point[i + 1] - change_point[i]))) +
                                      log(2 * M_PI) + 1);
/*          segmentation_likelihood -= ((double)(nb_variable * (change_point[i + 1] - change_point[i])) / 2.) *
                                     (logl(residual / (nb_variable * (change_point[i + 1] - change_point[i]))) +
                                      log(2 * M_PI)) + (double)(nb_variable * (change_point[i + 1] - change_point[i])) / 2.; */
          segment_variance[i] = residual / (nb_variable * (change_point[i + 1] - change_point[i]));
        }
        else {
          segmentation_likelihood = D_INF;
          break;
        }
      }

      if (segmentation_likelihood == D_INF) {
        break;
      }
    }

    if (model_type[0] == MEAN_CHANGE) {
//      if (residual > 0.) {
      if (residual > sqrt((double)(nb_variable * length[index])) * ROUNDOFF_ERROR) {
        segmentation_likelihood = -((double)(nb_variable * length[index]) / 2.) *
                                   (logl(residual / (nb_variable * length[index])) + log(2 * M_PI) + 1);
/*        segmentation_likelihood = -((double)(nb_variable * length[index]) / 2.) *
                                   (logl(residual / (nb_variable * (length[index] - nb_segment))) + log(2 * M_PI)) -
                                   (double)(nb_variable * (length[index] - nb_segment)) / 2.; */
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

      if (model_type[0] == MEAN_CHANGE) {
        os << SEQ_label[SEQL_GLOBAL_STANDARD_DEVIATION] << ": "
           << sqrt(residual / (nb_variable * length[index])) << endl;
      }
      else if (model_type[0] == MEAN_VARIANCE_CHANGE) {
        os << SEQ_label[SEQL_SEGMENT_STANDARD_DEVIATION] << ": ";
        for (i = 0;i < nb_segment;i++) {
          os << sqrt(segment_variance[i]) << ", ";
        }
        os << endl;
      }
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
    delete [] segment_variance;
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
  int max_nb_segment , max_nb_value , *frequency , *pisequence , *psegment , **optimal_length;
  double sum , factorial_sum , diff , buff , segmentation_likelihood , *mean ,
         *prsequence , **factorial , **forward;
  long double sum_square , *residual , *contrast;


  max_nb_segment = nb_segment[0];
  for (i = 1;i < nb_sequence;i++) {
    if (nb_segment[i] > max_nb_segment) {
      max_nb_segment = nb_segment[i];
    }
  }

  max_nb_value = 0;
  factorial = new double*[nb_variable];

  for (i = 1;i < nb_variable;i++) {
    if ((model_type[i - 1] == MULTINOMIAL_CHANGE) &&
        (marginal_distribution[i]->nb_value > max_nb_value)) {
      max_nb_value = marginal_distribution[i]->nb_value;
    }

    if (model_type[i - 1] == POISSON_CHANGE) {
      factorial[i] = new double[max_length];
    }
    else {
      factorial[i] = NULL;
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

  if (nb_parameter) {
    used_output = new bool[NB_OUTPUT];
  }

  segmentation_likelihood = 0.;

  for (i = 0;i < nb_sequence;i++) {
    for (j = 1;j < nb_variable;j++) {
      if (model_type[j - 1] == VARIANCE_CHANGE) {
        mean[j] = 0.;

        if (type[j] != REAL_VALUE) {
          pisequence = int_sequence[i][j];
          for (k = 0;k < length[i];k++) {
            mean[j] += *pisequence++;
          }
        }

        else {
          prsequence = real_sequence[i][j];
          for (k = 0;k < length[i];k++) {
            mean[j] += *prsequence++;
          }
        }

        mean[j] /= length[i];
      }
    }

    // recurrence "forward"

    for (j = 0;j < length[i];j++) {

      // calcul des log-vraisemblances des segments

      for (k = 0;k <= j;k++) {
        contrast[k] = 0.;
      }

      for (k = 1;k < nb_variable;k++) {
        if (model_type[k - 1] == MULTINOMIAL_CHANGE) {
          for (m = 0;m < marginal_distribution[k]->nb_value;m++) {
            frequency[m] = 0;
          }
          sum = 0.;

          pisequence = int_sequence[i][k] + j;
          frequency[*pisequence--]++;
          for (m = j - 1;m >= 0;m--) {
            sum += (j - m) * log((double)(j - m) / (double)(j - m + 1)) +
                   log((double)(frequency[*pisequence] + 1) / (double)(j - m + 1));
            if (frequency[*pisequence] > 0) {
              sum -= frequency[*pisequence] *
                     log((double)frequency[*pisequence] / (double)(frequency[*pisequence] + 1));
            }
            frequency[*pisequence--]++;

            if (contrast[m] != D_INF) {
              contrast[m] += sum;
            }

/*            frequency[*pisequence--]++;
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

          pisequence = int_sequence[i][k] + j;
          for (m = j;m >= 0;m--) {
            sum += *pisequence--;
            factorial_sum += factorial[k][m];
            if ((contrast[m] != D_INF) && (sum > 0.)) {
              contrast[m] += sum * (log(sum / (j - m + 1)) - 1) - factorial_sum;
            }
          }
        }

        else {
          if (model_type[k - 1] == VARIANCE_CHANGE) {
            sum_square = 0.;

            if (type[k] != REAL_VALUE) {
              pisequence = int_sequence[i][k] + j;
              for (m = j;m >= 0;m--) {
                diff = *pisequence-- - mean[k];
                sum_square += diff * diff;
                residual[m] = sum_square;
              }
            }
            else {
              prsequence = real_sequence[i][k] + j;
              for (m = j;m >= 0;m--) {
                diff = *prsequence-- - mean[k];
                sum_square += diff * diff;
                residual[m] = sum_square;
              }
            }
          }

          else if (model_type[k - 1] == ORDINAL_GAUSSIAN_CHANGE) {
            pisequence = int_sequence[i][k] + j;
            sum_square = 0.;
            sum = rank[k][*pisequence--];
            residual[j] = 0.;

            for (m = j - 1;m >= 0;m--) {
              diff = rank[k][*pisequence] - sum / (j - m);
              sum_square += ((double)(j - m) / (double)(j - m + 1)) * diff * diff;
              sum += rank[k][*pisequence--];
              residual[m] = sum_square;
            }
          }

          else {
            if (type[k] != REAL_VALUE) {
              pisequence = int_sequence[i][k] + j;
              sum_square = 0.;
              sum = *pisequence--;
              residual[j] = 0.;

              for (m = j - 1;m >= 0;m--) {
                diff = *pisequence - sum / (j - m);
                sum_square += ((double)(j - m) / (double)(j - m + 1)) * diff * diff;
                sum += *pisequence--;
                residual[m] = sum_square;
              }
            }

            else {
              prsequence = real_sequence[i][k] + j;
              sum_square = 0.;
              sum = *prsequence--;
              residual[j] = 0.;

              for (m = j - 1;m >= 0;m--) {
                diff = *prsequence - sum / (j - m);
                sum_square += ((double)(j - m) / (double)(j - m + 1)) * diff * diff;
                sum += *prsequence--;
                residual[m] = sum_square;
              }
            }
          }

          if (model_type[k - 1] == MEAN_CHANGE) {
            for (m = j - 1;m >= 0;m--) {
              contrast[m] -= residual[m];
            }
          }

          else if (model_type[k - 1] == MEAN_VARIANCE_CHANGE) {
            for (m = j - 1;m >= 0;m--) {
              contrast[m] += residual[m];
            }
          }

          else {
            for (m = j;m >= 0;m--) {
//              if ((contrast[m] != D_INF) && (residual[m] > 0.)) {
              if ((contrast[m] != D_INF) && (residual[m] > sqrt((double)(j - m + 1)) * ROUNDOFF_ERROR)) {
                contrast[m] -= ((double)(j - m + 1) / 2.) * (logl(residual[m] /
                                 (j - m + 1)) + log(2 * M_PI) + 1);
/*                contrast[m] -= ((double)(j - m + 1) / 2.) * (logl(residual[m] /
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
//           if (contrast[k] > 0.) {
          if (contrast[k] > sqrt((double)((nb_variable - 1) * (j - k + 1))) * ROUNDOFF_ERROR) {
            contrast[k] = -((double)((nb_variable - 1) * (j - k + 1)) / 2.) * (logl(contrast[k] /
                             ((nb_variable - 1) * (j - k + 1))) + log(2 * M_PI) + 1);
/*            contrast[k] = -((double)((nb_variable - 1) * (j - k + 1)) / 2.) * (logl(contrast[k] /
                             ((nb_variable - 1) * (j - k))) + log(2 * M_PI)) +
                           (double)((nb_variable - 1) * (j - k)) / 2.; */
          }
          else {
            contrast[k] = D_INF;
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

    if (model_type[0] != MEAN_CHANGE) {
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
        nb_parameter[j] = j;

        for (k = 1;k < nb_variable;k++) {
          if (model_type[k - 1] == MULTINOMIAL_CHANGE) {
            m = length[i] - 1;
            pisequence = int_sequence[i][k] + m;

            for (n = j;n >= 0;n--) {
              for (r = 0;r < marginal_distribution[k]->nb_value;r++) {
                used_output[r] = false;
              }

              used_output[*pisequence--] = true;
//              for (r = 0;r < optimal_length[m][n] - 1;r++) {
              for (r = m - 1;r > m - optimal_length[m][n];r--) {
                if (!used_output[*pisequence]) {
                  nb_parameter[j]++;
                  used_output[*pisequence] = true;
                }
                *pisequence--;
              }
              m -= optimal_length[m][n];
            }
          }

          else if ((model_type[k - 1] == POISSON_CHANGE) ||
                   (model_type[k - 1] == MEAN_CHANGE) ||
                   (model_type[k - 1] == MEAN_VARIANCE_CHANGE)) {
            nb_parameter[j] += j + 1;
          }

          else if (model_type[k - 1] == VARIANCE_CHANGE) {
            nb_parameter[j] += j + 2;
          }

          else {
            nb_parameter[j] += 2 * (j + 1);
          }
        }

        if (model_type[0] == MEAN_CHANGE) {
          nb_parameter[j]++;
        }
        else if (model_type[0] == MEAN_VARIANCE_CHANGE) {
          nb_parameter[j] += j + 1;
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

  if (nb_parameter) {
    delete [] used_output;
  }

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
  Sequences *seq , *iseq , *oseq;


  oseq = NULL;
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

  for (i = 0;i < nb_variable;i++) {
    if ((model_type[i] == MULTINOMIAL_CHANGE) || (model_type[i] == ORDINAL_GAUSSIAN_CHANGE) ||
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
                        << SEQ_error[SEQR_POSITIVE_MIN_VALUE];
          error.update((error_message.str()).c_str());
        }

        if (!marginal_distribution[i]) {
          status = false;
          ostringstream error_message;
          error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                        << STAT_error[STATR_MARGINAL_FREQUENCY_DISTRIBUTION];
          error.update((error_message.str()).c_str());
        }

        else if (model_type[i] == MULTINOMIAL_CHANGE) {
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

      if (((model_type[i] == MULTINOMIAL_CHANGE) || (model_type[i] == ORDINAL_GAUSSIAN_CHANGE)) &&
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

  for (i = 0;i < nb_sequence;i++) {
    if ((index == I_DEFAULT) || (index == i)) {
      if ((nb_segment[i] < 2) || (nb_segment[i] > length[i] / 2)) {
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
 *  Segmentation optimale d'une sequence.
 *
 *  arguments : reference sur un objet StatError, stream,
 *              identificateur de la sequence, nombre de segments maximum,
 *              types des modeles, sortie (sequence, entropies ou divergences).
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::segmentation(StatError &error , ostream &os , int iidentifier ,
                                   int max_nb_segment , int *model_type , int output) const

{
  bool status = true;
  register int i , j;
  int index , nb_segment , *nb_parameter , inb_segment[1] , ilength[4] , itype[1];
  double max_likelihood[2] , *segmentation_likelihood , *segment_penalty , *slope_penalty ,
         **penalized_likelihood , *likelihood , *uniform_entropy ,
         *segmentation_divergence , **rank;
  long double *segmentation_entropy , *ranked_change_point_entropy ,
              *change_point_entropy , *marginal_entropy;
  Sequences *iseq , *seq , *oseq;


  oseq = NULL;
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

  if ((model_type[0] == MEAN_CHANGE) && (output != SEQUENCE)) {
    status = false;
    error.correction_update(SEQ_error[SEQR_FORBIDDEN_OUTPUT] , SEQ_label[SEQL_SEQUENCE]);
  }

  for (i = 0;i < nb_variable;i++) {
    if ((model_type[i] == MULTINOMIAL_CHANGE) || (model_type[i] == ORDINAL_GAUSSIAN_CHANGE) ||
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
                        << SEQ_error[SEQR_POSITIVE_MIN_VALUE];
          error.update((error_message.str()).c_str());
        }

        if (!marginal_distribution[i]) {
          status = false;
          ostringstream error_message;
          error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                        << STAT_error[STATR_MARGINAL_FREQUENCY_DISTRIBUTION];
          error.update((error_message.str()).c_str());
        }

        else if (model_type[i] == MULTINOMIAL_CHANGE) {
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

  else if ((max_nb_segment < 2) || (max_nb_segment > length[index] / 2)) {
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

    segmentation_likelihood = new double[max_nb_segment + 2];
    nb_parameter = new int[max_nb_segment + 2];
    segment_penalty = new double[max_nb_segment + 2];
    slope_penalty = new double[max_nb_segment + 1];

    penalized_likelihood = new double*[3];
    for (i = 0;i < 3;i++) {
      penalized_likelihood[i] = new double[max_nb_segment + 1];
    }

    if (model_type[0] != MEAN_CHANGE) {
      likelihood = new double[max_nb_segment + 1];
      segmentation_entropy = new long double[max_nb_segment + 1];
      ranked_change_point_entropy = new long double[max_nb_segment + 1];
      change_point_entropy = new long double[max_nb_segment + 1];
      uniform_entropy = new double[max_nb_segment + 1];
      segmentation_divergence = new double[max_nb_segment + 1];
      marginal_entropy = new long double[max_nb_segment + 1];
    }

    inb_segment[0] = max_nb_segment + 1;

    seq->segmentation(inb_segment , model_type , rank , segmentation_likelihood + 1 ,
                      nb_parameter + 1 , segment_penalty + 1);
    if (model_type[0] != MEAN_CHANGE) {
      seq->forward_backward(0 , max_nb_segment , model_type , rank , NULL , NULL , SEGMENT , 'a' ,
                            likelihood + 1 , segmentation_entropy + 1 ,
                            ranked_change_point_entropy + 1 , change_point_entropy + 1 ,
                            uniform_entropy + 1 , marginal_entropy + 1);

      segmentation_divergence[1] = 0;
      for (i = 2;i <= max_nb_segment;i++) {
        segmentation_divergence[i] = uniform_entropy[i] - segmentation_entropy[i];
      }
    }

    // calcul des vraisemblances penalisees au sens du BIC, du BIC modifie (Zhang & Siegmund, 2007) et
    // du critere de Lavielle (2005)

//    segmentation_likelihood[1] = seq->one_segment_likelihood(0 , model_type , rank);
//    nb_parameter[1] = seq->nb_parameter_computation(0 , 1 , model_type);

    if (segmentation_likelihood[1] != D_INF) {
      penalized_likelihood[0][1] = 2 * segmentation_likelihood[1] - nb_parameter[1] *
                                   log((double)((seq->nb_variable - 1) * seq->length[0]));
      max_likelihood[0] = penalized_likelihood[0][1];

      penalized_likelihood[1][1] = penalized_likelihood[0][1] - segment_penalty[1];
      max_likelihood[1] = penalized_likelihood[1][1];
      nb_segment = 1;
    }

    else {
      max_nb_segment = 0;
      nb_segment = 0;
    }

//    inb_segment[0] = 2;
//    segmentation_likelihood[2] = seq->segmentation(inb_segment , model_type , rank);
//    nb_parameter[2] = seq->nb_parameter_computation(0 , 2 , model_type);

    if ((segmentation_likelihood[1] != D_INF) && (segmentation_likelihood[2] != D_INF)) {
      slope_penalty[1] = (segmentation_likelihood[2] - segmentation_likelihood[1]) /
                         MAX(nb_parameter[2] - nb_parameter[1] - 1 , 1);
    }
    else {
      slope_penalty[1] = D_INF;
    }

    if ((segmentation_likelihood[1] != D_INF) && (slope_penalty[1] != D_INF)) {
      penalized_likelihood[2][1] = 2 * (segmentation_likelihood[1] - nb_parameter[1] * slope_penalty[1]);
    }
    else {
      penalized_likelihood[2][1] = D_INF;
    }

    for (i = 2;i <= max_nb_segment;i++) {
//      inb_segment[0] = i + 1;
//      segmentation_likelihood[i + 1] = seq->segmentation(inb_segment , model_type , rank);
//      nb_parameter[i + 1] = seq->nb_parameter_computation(0 , i + 1 , model_type);

      if (segmentation_likelihood[i] != D_INF) {
        penalized_likelihood[0][i] = 2 * segmentation_likelihood[i] - nb_parameter[i] *
                                     log((double)((seq->nb_variable - 1) * seq->length[0]));
        if (penalized_likelihood[0][i] > max_likelihood[0]) {
          max_likelihood[0] = penalized_likelihood[0][i];
//          nb_segment = i;
        }

        penalized_likelihood[1][i] = penalized_likelihood[0][i] - segment_penalty[i];
        if (penalized_likelihood[1][i] > max_likelihood[1]) {
          max_likelihood[1] = penalized_likelihood[1][i];
          nb_segment = i;
        }
      }

      else {
        max_nb_segment = i - 1;
        break;
      }

      if ((segmentation_likelihood[i] != D_INF) && (segmentation_likelihood[i + 1] != D_INF)) {
        slope_penalty[i] = (segmentation_likelihood[i + 1] - segmentation_likelihood[i]) /
                           MAX(nb_parameter[i + 1] - nb_parameter[i] - 1 , 1);
      }
      else {
        slope_penalty[i] = D_INF;
      }

      if ((segmentation_likelihood[i] != D_INF) && (slope_penalty[i] != D_INF)) {
        penalized_likelihood[2][i] = 2 * (segmentation_likelihood[i] - (nb_parameter[i] - nb_segment - 1) *
                                     slope_penalty[i]);
      }
      else {
        penalized_likelihood[2][i] = D_INF;
      }
    }

    if (nb_segment > 0) {

#     ifdef MESSAGE
      int width[16];
      long old_adjust;
      double norm , *posterior_probability , **weight , *normalized_likelihood , *curvature;


      if (model_type[0] != MEAN_CHANGE) {
        posterior_probability = new double[max_nb_segment + 1];

        likelihood[1] = segmentation_likelihood[1];
        posterior_probability[1] = 1.;
        for (i = 2;i <= max_nb_segment;i++) {
          posterior_probability[i] = exp(segmentation_likelihood[i] - likelihood[i]);
        }
      }

      weight = new double*[2];
      for (i = 0;i < 2;i++) {
        weight[i] = new double[max_nb_segment + 1];
      }

      normalized_likelihood = new double[max_nb_segment + 2];
      curvature = new double[max_nb_segment + 1];

      norm = 0.;
      for (i = 1;i <= max_nb_segment;i++) {
        weight[0][i] = exp((penalized_likelihood[0][i] - max_likelihood[0]) / 2);
        norm += weight[0][i];
      }
      for (i = 1;i <= max_nb_segment;i++) {
        weight[0][i] /= norm;
      }

      norm = 0.;
      for (i = 1;i <= max_nb_segment;i++) {
        weight[1][i] = exp((penalized_likelihood[1][i] - max_likelihood[1]) / 2);
        norm += weight[1][i];
      }
      for (i = 1;i <= max_nb_segment;i++) {
        weight[1][i] /= norm;
      }

      for (i = 1;i <= max_nb_segment + 1;i++) {
        normalized_likelihood[i] = (segmentation_likelihood[max_nb_segment + 1] - segmentation_likelihood[i]) * max_nb_segment /
                                   (segmentation_likelihood[max_nb_segment + 1] - segmentation_likelihood[1]) + 1;
      }
      for (i = 2;i <= max_nb_segment;i++) {
        curvature[i] = normalized_likelihood[i + 1] - 2 * normalized_likelihood[i] + normalized_likelihood[i - 1];
      }

      old_adjust = os.flags(ios::adjustfield);

      if (model_type[0] != MEAN_CHANGE) {
        width[0] = column_width(max_nb_segment) + ASCII_SPACE;
        width[1] = column_width(max_nb_segment , segmentation_likelihood + 1 , 2.) + ASCII_SPACE;
        width[2] = column_width(max_nb_segment , likelihood + 1 , 2.) + ASCII_SPACE;
        width[3] = column_width(max_nb_segment , posterior_probability + 1) + ASCII_SPACE;
        width[4] = column_width(max_nb_segment - 1 , segmentation_entropy + 2) + ASCII_SPACE;
        width[5] = column_width(max_nb_segment - 1 , ranked_change_point_entropy + 2) + ASCII_SPACE;
        width[6] = column_width(max_nb_segment - 1 , change_point_entropy + 2) + ASCII_SPACE;
        width[7] = column_width(max_nb_segment - 1 , uniform_entropy + 2) + ASCII_SPACE;
        width[8] = column_width(max_nb_segment - 1 , segmentation_divergence + 2) + ASCII_SPACE;
        width[9] = column_width(max_nb_segment - 1 , marginal_entropy + 2) + ASCII_SPACE;
        width[10] = column_width(nb_parameter[max_nb_segment]) + ASCII_SPACE;
//        width[11] = column_width(max_nb_segment , penalized_likelihood[0] + 1) + ASCII_SPACE;
//        width[12] = column_width(max_nb_segment , weight[0] + 1) + ASCII_SPACE;
        width[13] = column_width(max_nb_segment , penalized_likelihood[1] + 1) + ASCII_SPACE;
        width[14] = column_width(max_nb_segment , weight[1] + 1) + ASCII_SPACE;
        width[15] = column_width(max_nb_segment - 1 , curvature + 2) + ASCII_SPACE;

        os << "\n" << SEQ_label[SEQL_NB_SEGMENT] << " | 2 * " << STAT_label[STATL_LIKELIHOOD]
           << " | 2 * " << SEQ_label[SEQL_POSSIBLE_SEGMENTATION_LIKELIHOOD]
           << " | " << SEQ_label[SEQL_POSTERIOR_PROBABILITY]
           << " | " << SEQ_label[SEQL_SEGMENTATION_ENTROPY]
           << " | " << SEQ_label[SEQL_RANKED_CHANGE_POINT_ENTROPY]
           << " | " << SEQ_label[SEQL_CHANGE_POINT_ENTROPY]
           << " | " << SEQ_label[SEQL_UNIFORM_ENTROPY]
           << " | " << SEQ_label[SEQL_SEGMENTATION_DIVERGENCE]
           << " | " << SEQ_label[SEQL_MARGINAL_ENTROPY]
           << " | " << STAT_label[STATL_FREE_PARAMETERS]
//           << " | "  << STAT_criterion_word[BIC] << " - "  << STAT_label[STATL_WEIGHT]
           << " | Modified " << STAT_criterion_word[BIC] << " - "  << STAT_label[STATL_WEIGHT]
           << " | Lavielle criterion (standardized curvature)" << endl;

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
           << setw(width[9]) << " "
           << setw(width[10]) << nb_parameter[1]
//           << setw(width[11]) << penalized_likelihood[0][1]
//           << setw(width[12]) << weight[0][1]
           << setw(width[13]) << penalized_likelihood[1][1]
           << setw(width[14]) << weight[1][1]
           << setw(width[15]) << " " << endl;

        for (i = 2;i <= max_nb_segment;i++) {
          os << setw(width[0]) << i
             << setw(width[1]) << 2 * segmentation_likelihood[i]
             << setw(width[2]) << 2 * likelihood[i]
             << setw(width[3]) << posterior_probability[i]
             << setw(width[4]) << segmentation_entropy[i]
             << setw(width[5]) << ranked_change_point_entropy[i]
             << setw(width[6]) << change_point_entropy[i]
             << setw(width[7]) << uniform_entropy[i]
             << setw(width[8]) << segmentation_divergence[i]
             << setw(width[9]) << marginal_entropy[i]
             << setw(width[10]) << nb_parameter[i]
//             << setw(width[11]) << penalized_likelihood[0][i]
//             << setw(width[12]) << weight[0][i]
             << setw(width[13]) << penalized_likelihood[1][i]
             << setw(width[14]) << weight[1][i]
             << setw(width[15]) << curvature[i] << endl;
        }
        os << endl;
      }

      else {
        width[0] = column_width(max_nb_segment) + ASCII_SPACE;
        width[1] = column_width(max_nb_segment , segmentation_likelihood + 1 , 2.) + ASCII_SPACE;
        width[10] = column_width(nb_parameter[max_nb_segment]) + ASCII_SPACE;
//        width[11] = column_width(max_nb_segment , penalized_likelihood[0] + 1) + ASCII_SPACE;
//        width[12] = column_width(max_nb_segment , weight[0] + 1) + ASCII_SPACE;
        width[13] = column_width(max_nb_segment , penalized_likelihood[1] + 1) + ASCII_SPACE;
        width[14] = column_width(max_nb_segment , weight[1] + 1) + ASCII_SPACE;
        width[15] = column_width(max_nb_segment - 1 , curvature + 2) + ASCII_SPACE;

        os << "\n" << SEQ_label[SEQL_NB_SEGMENT] << " | 2 * " << STAT_label[STATL_LIKELIHOOD]
           << " | " << STAT_label[STATL_FREE_PARAMETERS]
//           << " | "  << STAT_criterion_word[BIC] << " - "  << STAT_label[STATL_WEIGHT]
           << " | Modified " << STAT_criterion_word[BIC] << " - "  << STAT_label[STATL_WEIGHT]
           << " | Lavielle criterion (standardized curvature)" << endl;

        os.setf(ios::left , ios::adjustfield);

        os << setw(width[0]) << 1
           << setw(width[1]) << 2 * segmentation_likelihood[1]
           << setw(width[10]) << nb_parameter[1]
//           << setw(width[11]) << penalized_likelihood[0][1]
//           << setw(width[12]) << weight[0][1]
           << setw(width[13]) << penalized_likelihood[1][1]
           << setw(width[14]) << weight[1][1]
           << setw(width[15]) << " " << endl;

        for (i = 2;i <= max_nb_segment;i++) {
          os << setw(width[0]) << i
             << setw(width[1]) << 2 * segmentation_likelihood[i]
             << setw(width[10]) << nb_parameter[i]
//             << setw(width[11]) << penalized_likelihood[0][i]
//             << setw(width[12]) << weight[0][i]
             << setw(width[13]) << penalized_likelihood[1][i]
             << setw(width[14]) << weight[1][i]
             << setw(width[15]) << curvature[i] << endl;
        }
        os << endl;
      }

      os.setf((FMTFLAGS)old_adjust , ios::adjustfield);

/*      for (i = 1;i <= max_nb_segment;i++) {
        penalized_likelihood[1][i] = 2 * segmentation_likelihood[i] - nb_parameter[i] *
                                     log((double)seq->length[0]) - segment_penalty[i];
//        penalized_likelihood[1][i] += nb_parameter[i] *
//                                      log((double)((seq->nb_variable - 2) * seq->length[0]));
      }

      max_likelihood[1] = penalized_likelihood[1][1];
      for (i = 2;i <= max_nb_segment;i++) {
        if (penalized_likelihood[1][i] > max_likelihood[1]) {
          max_likelihood[1] = penalized_likelihood[1][i];
        }
      }

      norm = 0.;
      for (i = 1;i <= max_nb_segment;i++) {
        weight[1][i] = exp((penalized_likelihood[1][i] - max_likelihood[1]) / 2);
        norm += weight[1][i];
      }
      for (i = 1;i <= max_nb_segment;i++) {
        weight[1][i] /= norm;
      } */

/*      for (i = 1;i <= max_nb_segment;i++) {
        os << "\n" << i << " " << (i == 1 ? SEQ_label[SEQL_SEGMENT] : SEQ_label[SEQL_SEGMENTS])
           << "   2 * " << STAT_label[STATL_LIKELIHOOD] << ": " << 2 * segmentation_likelihood[i] << "   "
           << nb_parameter[i] << " " << STAT_label[nb_parameter[i] == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS]
           << "   2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " (Modified "  << STAT_criterion_word[BIC] << "): "
           << penalized_likelihood[1][i] << "   " << STAT_label[STATL_WEIGHT] << ": " << weight[1][i] << endl;
      }

      for (i = 1;i <= max_nb_segment;i++) {
        os << "\n" << i << " " << (i == 1 ? SEQ_label[SEQL_SEGMENT] : SEQ_label[SEQL_SEGMENTS])
           << "   2 * " << STAT_label[STATL_LIKELIHOOD] << ": " << 2 * segmentation_likelihood[i] << "   "
           << nb_parameter[i] << " " << STAT_label[nb_parameter[i] == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS]
           << "   " << SEQ_label[SEQL_PENALTY] << ": " << slope_penalty[i] << "\n"
           << "2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " (" << "Lavielle criterion" << "): "
           << penalized_likelihood[2][i] << "   " << normalized_likelihood[i];
        if (i > 1) {
          os << "   " << "curvature" << ": " << curvature[i];
//           << segmentation_likelihood[i - 1] - 2 * segmentation_likelihood[i] + segmentation_likelihood[i + 1];
        }
        os << endl;
      }

      os << endl; */

      if (model_type[0] != MEAN_CHANGE) {
        delete [] posterior_probability;
      }

      for (i = 0;i < 2;i++) {
        delete [] weight[i];
      }
      delete [] weight;

      delete [] normalized_likelihood;
      delete [] curvature;
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
          oseq->real_sequence[1][0][i - 2] = ranked_change_point_entropy[i];
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
      }

#     ifdef DEBUG
      if ((model_type[0] != MEAN_CHANGE) && (model_type[0] != MEAN_VARIANCE_CHANGE)) {
        hierarchical_segmentation(error , os , iidentifier , max_nb_segment , model_type);
      }
#     endif

    }

    else {
      oseq = NULL;
      error.update(SEQ_error[SEQR_SEGMENTATION_FAILURE]);
    }

    for (i = 1;i < seq->nb_variable;i++) {
      delete [] rank[i];
    }
    delete [] rank;

    delete seq;

    delete [] segmentation_likelihood;
    delete [] nb_parameter;
    delete [] segment_penalty;
    delete [] slope_penalty;

    for (i = 0;i < 3;i++) {
      delete [] penalized_likelihood[i];
    }
    delete [] penalized_likelihood;

    if (model_type[0] != MEAN_CHANGE) {
      delete [] likelihood;
      delete [] segmentation_entropy;
      delete [] ranked_change_point_entropy;
      delete [] change_point_entropy;
      delete [] uniform_entropy;
      delete [] segmentation_divergence;
      delete [] marginal_entropy;
    }
  }

  return oseq;
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
  long double sum_square , merge_contrast , *residual , *begin_contrast , *end_contrast ,
              **begin_sum_square , **end_sum_square;
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

  if ((model_type[0] == MEAN_CHANGE) || (model_type[0] == MEAN_VARIANCE_CHANGE)) {
    status = false;
    error.update(SEQ_error[SEQR_CHANGE_POINT_MODEL]);
  }

  for (i = 0;i < nb_variable;i++) {
    if ((model_type[i] == MULTINOMIAL_CHANGE) || (model_type[i] == ORDINAL_GAUSSIAN_CHANGE) ||
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
                        << SEQ_error[SEQR_POSITIVE_MIN_VALUE];
          error.update((error_message.str()).c_str());
        }

        if (!marginal_distribution[i]) {
          status = false;
          ostringstream error_message;
          error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                        << STAT_error[STATR_MARGINAL_FREQUENCY_DISTRIBUTION];
          error.update((error_message.str()).c_str());
        }

        else if (model_type[i] == MULTINOMIAL_CHANGE) {
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
    begin_sum_square = new long double*[seq->nb_variable];
    end_sum_square = new long double*[seq->nb_variable];

    for (i = 1;i < seq->nb_variable;i++) {
      if ((model_type[i - 1] == MULTINOMIAL_CHANGE) &&
          (seq->marginal_distribution[i]->nb_value > max_nb_value)) {
        max_nb_value = seq->marginal_distribution[i]->nb_value;
      }

      if (model_type[i - 1] == MULTINOMIAL_CHANGE) {
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
        begin_sum_square[i] = new long double[seq->length[0]];
        end_sum_square[i] = new long double[seq->length[0]];
      }
      else {
        begin_sum_square[i] = NULL;
        end_sum_square[i] = NULL;
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
      if (model_type[i - 1] == MULTINOMIAL_CHANGE) {
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
          sum_square = 0.;

          if (seq->type[i] != REAL_VALUE) {
            pisequence = seq->int_sequence[0][i];
            for (j = 0;j < seq->length[0];j++) {
              diff = *pisequence++ - mean[i];
              sum_square += diff * diff;

              begin_sum_square[i][j] = sum_square;

              residual[j] = sum_square;
            }
          }

          else {
            prsequence = seq->real_sequence[0][i];
            for (j = 0;j < seq->length[0];j++) {
              diff = *prsequence++ - mean[i];
              sum_square += diff * diff;

              begin_sum_square[i][j] = sum_square;

              residual[j] = sum_square;
            }
          }
          break;
        }

        case ORDINAL_GAUSSIAN_CHANGE : {
          pisequence = seq->int_sequence[0][i];
          sum_square = 0.;
          sum = rank[i][*pisequence++];

          begin_sum_square[i][0] = sum_square;
          begin_sum[i][0] = sum;

          residual[0] = 0.;

          for (j = 1;j < seq->length[0];j++) {
            diff = rank[i][*pisequence] - sum / j;
            sum_square += ((double)j / (double)(j + 1)) * diff * diff;
            sum += rank[i][*pisequence++];

            begin_sum_square[i][j] = sum_square;
            begin_sum[i][j] = sum;

            residual[j] = sum_square;
          }
          break;
        }

        case GAUSSIAN_CHANGE : {
          if (seq->type[i] != REAL_VALUE) {
            pisequence = seq->int_sequence[0][i];
            sum_square = 0.;
            sum = *pisequence++;

            begin_sum_square[i][0] = sum_square;
            begin_sum[i][0] = sum;

            residual[0] = 0.;

            for (j = 1;j < seq->length[0];j++) {
              diff = *pisequence - sum / j;
              sum_square += ((double)j / (double)(j + 1)) * diff * diff;
              sum += *pisequence++;

              begin_sum_square[i][j] = sum_square;
              begin_sum[i][j] = sum;

              residual[j] = sum_square;
            }
          }

          else {
            prsequence = seq->real_sequence[0][i];
            sum_square = 0.;
            sum = *prsequence++;

            begin_sum_square[i][0] = sum_square;
            begin_sum[i][0] = sum;

            residual[0] = 0.;

            for (j = 1;j < seq->length[0];j++) {
              diff = *prsequence - sum / j;
              sum_square += ((double)j / (double)(j + 1)) * diff * diff;
              sum += *prsequence++;

              begin_sum_square[i][j] = sum_square;
              begin_sum[i][j] = sum;

              residual[j] = sum_square;
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
      if (model_type[i - 1] == MULTINOMIAL_CHANGE) {
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
          sum_square = 0.;

          if (seq->type[i] != REAL_VALUE) {
            pisequence = seq->int_sequence[0][i] + seq->length[0] - 1;
            for (j = seq->length[0] - 1;j >= 0;j--) {
              diff = *pisequence-- - mean[i];
              sum_square += diff * diff;

              end_sum_square[i][j] = sum_square;

              residual[j] = sum_square;
            }
          }

          else {
            prsequence = seq->real_sequence[0][i] + seq->length[0] - 1;
            for (j = seq->length[0] - 1;j >= 0;j--) {
              diff = *prsequence-- - mean[i];
              sum_square += diff * diff;

              end_sum_square[i][j] = sum_square;

              residual[j] = sum_square;
            }
          }
          break;
        }

        case ORDINAL_GAUSSIAN_CHANGE : {
          pisequence = seq->int_sequence[0][i] + seq->length[0] - 1;
          sum_square = 0.;
          sum = rank[i][*pisequence--];

          end_sum_square[i][seq->length[0] - 1] = sum_square;
          end_sum[i][seq->length[0] - 1] = sum;

          residual[seq->length[0] - 1] = 0.;

          for (j = seq->length[0] - 2;j >= 0;j--) {
            diff = rank[i][*pisequence] - sum / (seq->length[0] - j - 1);
            sum_square += ((double)(seq->length[0] - j - 1) / (double)(seq->length[0] - j)) * diff * diff;
            sum += rank[i][*pisequence--];

            end_sum_square[i][j] = sum_square;
            end_sum[i][j] = sum;

            residual[j] = sum_square;
          }
          break;
        }

        case GAUSSIAN_CHANGE : {
          if (seq->type[i] != REAL_VALUE) {
            pisequence = seq->int_sequence[0][i] + seq->length[0] - 1;
            sum_square = 0.
            sum = *pisequence--;

            end_sum_square[i][seq->length[0] - 1] = sum_square;
            end_sum[i][seq->length[0] - 1] = sum;

            residual[seq->length[0] - 1] = 0.;

            for (j = seq->length[0] - 2;j >= 0;j--) {
              diff = *pisequence - sum / (seq->length[0] - j - 1);
              sum_square += ((double)(seq->length[0] - j - 1) / (double)(seq->length[0] - j)) * diff * diff;
              sum += *pisequence--;

              end_sum_square[i][j] = sum_square;
              end_sum[i][j] = sum;

              residual[j] = sum_square;
            }
          }

          else {
            prsequence = seq->real_sequence[0][i] + seq->length[0] - 1;
            sum_square = 0.;
            sum = *prsequence--;

            end_sum_square[i][seq->length[0] - 1] = sum_square;
            end_sum[i][seq->length[0] - 1] = sum;

            residual[seq->length[0] - 1] = 0.;

            for (j = seq->length[0] - 2;j >= 0;j--) {
              diff = *prsequence - sum / (seq->length[0] - j - 1);
              sum_square += ((double)(seq->length[0] - j - 1) / (double)(seq->length[0] - j)) * diff * diff;
              sum += *prsequence--;

              end_sum_square[i][j] = sum_square;
              end_sum[i][j] = sum;

              residual[j] = sum_square;
            }
          }
          break;
        }
        } */

        if (model_type[i - 1] == VARIANCE_CHANGE) {
          for (j = seq->length[0] - 1;j > 0;j--) {
            end_sum_square[i][j] = begin_sum_square[i][seq->length[0] - 1] - begin_sum_square[i][j - 1];
            residual[j] = end_sum_square[i][j];
          }

          end_sum_square[i][0] = begin_sum_square[i][seq->length[0] - 1];
          residual[0] = end_sum_square[i][0];
        }

        else {
          for (j = seq->length[0] - 1;j > 0;j--) {
            end_sum[i][j] = begin_sum[i][seq->length[0] - 1] - begin_sum[i][j - 1];
            diff = begin_sum[i][j - 1] / j  - end_sum[i][j] / (seq->length[0] - j);
            end_sum_square[i][j] = begin_sum_square[i][seq->length[0] - 1] - begin_sum_square[i][j - 1] -
                                   ((double)(j * (seq->length[0] - j)) / (double)seq->length[0]) * diff * diff;

            residual[j] = end_sum_square[i][j];
          }

          end_sum_square[i][0] = begin_sum_square[i][seq->length[0] - 1];
          end_sum[i][0] = begin_sum[i][seq->length[0] - 1];

          residual[0] = end_sum_square[i][0];
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
      for (i = 0;i <= nb_segment;i++) {
        os << " " << change_point[nb_segment][i];
      }
      os << endl;
#     endif

    }

//    while ((penalized_likelihood[nb_segment] > penalized_likelihood[nb_segment - 1]) &&
    while (nb_segment < max_nb_segment) {

      // calcul des log-vraisemblances des segments modifies

      for (i = begin_change_point;i < change_point[nb_segment][segment_index + 1];i++) {
        begin_contrast[i] = 0.;
      }

      for (i = 1;i < seq->nb_variable;i++) {
        if (model_type[i - 1] == MULTINOMIAL_CHANGE) {
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
              sum_square = begin_sum_square[i][begin_change_point - 1];

              if (seq->type[i] != REAL_VALUE) {
                pisequence = seq->int_sequence[0][i] + begin_change_point;
                for (j = begin_change_point;j < split_change_point;j++) {
                  diff = *pisequence++ - mean[i];
                  sum_square += diff * diff;

                  begin_sum_square[i][j] = sum_square;

                  residual[j] = sum_square;
                }
              }

              else {
                prsequence = seq->real_sequence[0][i] + begin_change_point;
                for (j = begin_change_point;j < split_change_point;j++) {
                  diff = *prsequence++ - mean[i];
                  sum_square += diff * diff;

                  begin_sum_square[i][j] = sum_square;

                  residual[j] = sum_square;
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

            sum_square = 0.;

            if (seq->type[i] != REAL_VALUE) {
              for (j = split_change_point;j < change_point[nb_segment][segment_index + 1];j++) {
                diff = *pisequence++ - mean[i];
                sum_square += diff * diff;

                begin_sum_square[i][j] = sum_square;

                residual[j] = sum_square;
              }
            }

            else {
              for (j = split_change_point;j < change_point[nb_segment][segment_index + 1];j++) {
                diff = *prsequence++ - mean[i];
                sum_square += diff * diff;

                begin_sum_square[i][j] = sum_square;

                residual[j] = sum_square;
              }
            }
            break;
          }

          case ORDINAL_GAUSSIAN_CHANGE : {
            if (begin_change_point < split_change_point) {
              sum_square = begin_sum_square[i][begin_change_point - 1];
              sum = begin_sum[i][begin_change_point - 1];

              pisequence = seq->int_sequence[0][i] + begin_change_point;
              for (j = begin_change_point;j < split_change_point;j++) {
                diff = rank[i][*pisequence] - sum / (j - change_point[nb_segment][segment_index - 1]);
                sum_square += ((double)(j - change_point[nb_segment][segment_index - 1]) /
                               (double)(j - change_point[nb_segment][segment_index - 1] + 1)) * diff * diff;
                sum += rank[i][*pisequence++];
 
                begin_sum_square[i][j] = sum_square;
                begin_sum[i][j] = sum;

                residual[j] = sum_square;
              }
            }

            else {
              pisequence = seq->int_sequence[0][i] + split_change_point;
            }

            sum_square = 0.;
            sum = rank[i][*pisequence++];

            begin_sum_square[i][split_change_point] = sum_square;
            begin_sum[i][split_change_point] = sum;

            residual[split_change_point] = 0.;

            for (j = split_change_point + 1;j < change_point[nb_segment][segment_index + 1];j++) {
              diff = rank[i][*pisequence] - sum / (j - split_change_point);
              sum_square += ((double)(j - split_change_point) /
                             (double)(j - split_change_point + 1)) * diff * diff;
              sum += rank[i][*pisequence++];
 
              begin_sum_square[i][j] = sum_square;
              begin_sum[i][j] = sum;

              residual[j] = sum_square;
            }
            break;
          }

          case GAUSSIAN_CHANGE : {
            if (begin_change_point < split_change_point) {
              sum_square = begin_sum_square[i][begin_change_point - 1];
              sum = begin_sum[i][begin_change_point - 1];

              pisequence = seq->int_sequence[0][i] + begin_change_point;
              if (seq->type[i] != REAL_VALUE) {
                for (j = begin_change_point;j < split_change_point;j++) {
                  diff = *pisequence - sum / (j - change_point[nb_segment][segment_index - 1]);
                  sum_square += ((double)(j - change_point[nb_segment][segment_index - 1]) /
                                 (double)(j - change_point[nb_segment][segment_index - 1] + 1)) * diff * diff;
                  sum += *pisequence++;

                  begin_sum_square[i][j] = sum_square;
                  begin_sum[i][j] = sum;

                  residual[j] = sum_square;
                }
              }

              else {
                prsequence = seq->real_sequence[0][i] + begin_change_point;
                for (j = begin_change_point;j < split_change_point;j++) {
                  diff = *prsequence - sum / (j - change_point[nb_segment][segment_index - 1]);
                  sum_square += ((double)(j - change_point[nb_segment][segment_index - 1]) /
                                 (double)(j - change_point[nb_segment][segment_index - 1] + 1)) * diff * diff;
                  sum += *prsequence++;

                  begin_sum_square[i][j] = sum_square;
                  begin_sum[i][j] = sum;

                  residual[j] = sum_square;
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
              sum_square = 0.;
              sum = *pisequence++;

              begin_sum_square[i][split_change_point] = sum_square;
              begin_sum[i][split_change_point] = sum;

              residual[split_change_point] = 0.;

              for (j = split_change_point + 1;j < change_point[nb_segment][segment_index + 1];j++) {
                diff = *pisequence - sum / (j - split_change_point);
                sum_square += ((double)(j - split_change_point) /
                               (double)(j - split_change_point + 1)) * diff * diff;
                sum += *pisequence++;

                begin_sum_square[i][j] = sum_square;
                begin_sum[i][j] = sum;

                residual[j] = sum_square;
              }
            }

            else {
              sum_square = 0.;
              sum = *prsequence++;

              begin_sum_square[i][split_change_point] = sum_square;
              begin_sum[i][split_change_point] = sum;

              residual[split_change_point] = 0.;

              for (j = split_change_point + 1;j < change_point[nb_segment][segment_index + 1];j++) {
                diff = *prsequence - sum / (j - split_change_point);
                sum_square += ((double)(j - split_change_point) /
                               (double)(j - split_change_point + 1)) * diff * diff;
                sum += *prsequence++;

                begin_sum_square[i][j] = sum_square;
                begin_sum[i][j] = sum;

                residual[j] = sum_square;
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
        if (model_type[i - 1] == MULTINOMIAL_CHANGE) {
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
              sum_square = end_sum_square[i][end_change_point];

              if (seq->type[i] != REAL_VALUE) {
                pisequence = seq->int_sequence[0][i] + end_change_point - 1;
                for (j = end_change_point - 1;j >= split_change_point;j--) {
                  diff = *pisequence-- - mean[i];
                  sum_square += diff * diff;

                  end_sum_square[i][j] = sum_square;

                  residual[j] = sum_square;
                }
              }

              else {
                prsequence = seq->real_sequence[0][i] + end_change_point - 1;
                for (j = end_change_point - 1;j >= split_change_point;j--) {
                  diff = *prsequence-- - mean[i];
                  sum_square += diff * diff;

                  end_sum_square[i][j] = sum_square;

                  residual[j] = sum_square;
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

            sum_square = 0.;

            if (seq->type[i] != REAL_VALUE) {
              for (j = split_change_point - 1;j >= change_point[nb_segment][segment_index - 1];j--) {
                diff = *pisequence-- - mean[i];
                sum_square += diff * diff;

                end_sum_square[i][j] = sum_square;

                residual[j] = sum_square;
              }
            }

            else {
              for (j = split_change_point - 1;j >= change_point[nb_segment][segment_index - 1];j--) {
                diff = *prsequence-- - mean[i];
                sum_square += diff * diff;

                end_sum_square[i][j] = sum_square;

                residual[j] = sum_square;
              }
            }
            break;
          }

          case ORDINAL_GAUSSIAN_CHANGE : {
            if (end_change_point > split_change_point) {
              sum_square = end_sum_square[i][end_change_point];
              sum = end_sum[i][end_change_point];

              pisequence = seq->int_sequence[0][i] + end_change_point - 1;
              for (j = end_change_point - 1;j >= split_change_point;j--) {
                diff = rank[i][*pisequence] - sum / (change_point[nb_segment][segment_index + 1] - j - 1);
                sum_square += ((double)(change_point[nb_segment][segment_index + 1] - j - 1) /
                               (double)(change_point[nb_segment][segment_index + 1] - j)) * diff * diff;
                sum += rank[i][*pisequence--];

                end_sum_square[i][j] = sum_square;
                end_sum[i][j] = sum;

                residual[j] = sum_square;
              }
            }

            else {
              pisequence = seq->int_sequence[0][i] + split_change_point - 1;
            }

            sum_square = 0.;
            sum = rank[i][*pisequence--];

            end_sum_square[i][split_change_point - 1] = sum_square;
            end_sum[i][split_change_point - 1] = sum;

            residual[split_change_point - 1] = 0.;

            for (j = split_change_point - 2;j >= change_point[nb_segment][segment_index - 1];j--) {
              diff = rank[i][*pisequence] - sum / (split_change_point - j - 1);
              sum_square += ((double)(split_change_point - j - 1) /
                             (double)(split_change_point - j)) * diff * diff;
              sum += rank[i][*pisequence--];

              end_sum_square[i][j] = sum_square;
              end_sum[i][j] = sum;

              residual[j] = sum_square;
            }
            break;
          }

          case GAUSSIAN_CHANGE : {
            if (end_change_point > split_change_point) {
              sum_square = end_sum_square[i][end_change_point];
              sum = end_sum[i][end_change_point];

              if (seq->type[i] != REAL_VALUE) {
                pisequence = seq->int_sequence[0][i] + end_change_point - 1;
                for (j = end_change_point - 1;j >= split_change_point;j--) {
                  diff = *pisequence - sum / (change_point[nb_segment][segment_index + 1] - j - 1);
                  sum_square += ((double)(change_point[nb_segment][segment_index + 1] - j - 1) /
                                 (double)(change_point[nb_segment][segment_index + 1] - j)) * diff * diff;
                  sum += *pisequence--;

                  end_sum_square[i][j] = sum_square;
                  end_sum[i][j] = sum;

                  residual[j] = sum_square;
                }
              }

              else {
                prsequence = seq->real_sequence[0][i] + end_change_point - 1;
                for (j = end_change_point - 1;j >= split_change_point;j--) {
                  diff = *prsequence - sum / (change_point[nb_segment][segment_index + 1] - j - 1);
                  sum_square += ((double)(change_point[nb_segment][segment_index + 1] - j - 1) /
                                 (double)(change_point[nb_segment][segment_index + 1] - j)) * diff * diff;
                  sum += *prsequence--;

                  end_sum_square[i][j] = sum_square;
                  end_sum[i][j] = sum;

                  residual[j] = sum_square;
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
              sum_square = 0.;
              sum = *pisequence--;

              end_sum_square[i][split_change_point - 1] = sum_square;
              end_sum[i][split_change_point - 1] = sum;

              residual[split_change_point - 1] = 0.;

              for (j = split_change_point - 2;j >= change_point[nb_segment][segment_index - 1];j--) {
                diff = *pisequence - sum / (split_change_point - j - 1);
                sum_square += ((double)(split_change_point - j - 1) /
                               (double)(split_change_point - j)) * diff * diff;
                sum += *pisequence--;

                end_sum_square[i][j] = sum_square;
                end_sum[i][j] = sum;

                residual[j] = sum_square;
              }
            }

            else {
              sum_square = 0.;
              sum = *prsequence--;

              end_sum_square[i][split_change_point - 1] = sum_square;
              end_sum[i][split_change_point - 1] = sum;

              residual[split_change_point - 1] = 0.;

              for (j = split_change_point - 2;j >= change_point[nb_segment][segment_index - 1];j--) {
                diff = *prsequence - sum / (split_change_point - j - 1);
                sum_square += ((double)(split_change_point - j - 1) /
                               (double)(split_change_point - j)) * diff * diff;
                sum += *prsequence--;

                end_sum_square[i][j] = sum_square;
                end_sum[i][j] = sum;

                residual[j] = sum_square;
              }
            }
            break;
          }
          } */

          if (model_type[i - 1] == VARIANCE_CHANGE) {
            if (end_change_point > split_change_point) {
              for (j = end_change_point - 1;j > split_change_point;j--) {
                end_sum_square[i][j] = begin_sum_square[i][change_point[nb_segment][segment_index + 1] - 1] -
                                       begin_sum_square[i][j - 1];
                residual[j] = end_sum_square[i][j];
              }

              end_sum_square[i][j] = begin_sum_square[i][change_point[nb_segment][segment_index + 1] - 1];
              residual[j] = end_sum_square[i][j];
            }

            for (j = split_change_point - 1;j > change_point[nb_segment][segment_index - 1];j--) {
              end_sum_square[i][j] = begin_sum_square[i][split_change_point - 1] -
                                     begin_sum_square[i][j - 1];
              residual[j] = end_sum_square[i][j];
            }

            end_sum_square[i][j] = begin_sum_square[i][split_change_point - 1];
            residual[j] = end_sum_square[i][j];
          }

          else {
            if (end_change_point > split_change_point) {
              for (j = end_change_point - 1;j > split_change_point;j--) {
                end_sum[i][j] = begin_sum[i][change_point[nb_segment][segment_index + 1] - 1] -
                                begin_sum[i][j - 1];
                diff = begin_sum[i][j - 1] / (j - split_change_point) -
                       end_sum[i][j] / (change_point[nb_segment][segment_index + 1] - j);
                end_sum_square[i][j] = begin_sum_square[i][change_point[nb_segment][segment_index + 1] - 1] -
                                       begin_sum_square[i][j - 1] - ((double)((j - split_change_point) *
                                         (change_point[nb_segment][segment_index + 1] - j)) /
                                        (double)(change_point[nb_segment][segment_index + 1] - split_change_point)) * diff * diff;

                residual[j] = end_sum_square[i][j];
              }

              end_sum_square[i][j] = begin_sum_square[i][change_point[nb_segment][segment_index + 1] - 1];
              end_sum[i][j] = begin_sum[i][change_point[nb_segment][segment_index + 1] - 1];

              residual[j] = end_sum_square[i][j];
            }

            for (j = split_change_point - 1;j > change_point[nb_segment][segment_index - 1];j--) {
              end_sum[i][j] = begin_sum[i][split_change_point - 1] - begin_sum[i][j - 1];
              diff = begin_sum[i][j - 1] / (j - change_point[nb_segment][segment_index - 1]) -
                     end_sum[i][j] / (split_change_point - j);
              end_sum_square[i][j] = begin_sum_square[i][split_change_point - 1] - begin_sum_square[i][j - 1] -
                                     ((double)((j - change_point[nb_segment][segment_index - 1]) *
                                       (split_change_point - j)) /
                                      (double)(split_change_point - change_point[nb_segment][segment_index - 1])) * diff * diff;

              residual[j] = end_sum_square[i][j];
            }

            end_sum_square[i][j] = begin_sum_square[i][split_change_point - 1];
            end_sum[i][j] = begin_sum[i][split_change_point - 1];

            residual[j] = end_sum_square[i][j];
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
        cout << "\nBegin merge: " << change_point[nb_segment][segment_index - 1] << " " << begin_change_point << " " << split_change_point
             << " | " << begin_contrast[begin_change_point - 1] << " " << begin_contrast[begin_change_point] << "\n" << endl;
      }
      if (end_change_point > split_change_point) {
        cout << "\nEnd merge: " << split_change_point << " " << end_change_point << " " << change_point[nb_segment][segment_index + 1]
             << " | " << end_contrast[end_change_point] << " " << end_contrast[end_change_point - 1] << "\n" << endl;
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
            if (model_type[i - 1] == MULTINOMIAL_CHANGE) {
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
                buff = begin_sum_square[i][change_point[nb_segment][segment_index - 1] - 1] +
                       begin_sum_square[i][split_change_point - 1];
              }

              else {
                diff = begin_sum[i][change_point[nb_segment][segment_index - 1] - 1] /
                       (change_point[nb_segment][segment_index - 1] - change_point[nb_segment][segment_index - 2]) -
                       begin_sum[i][split_change_point - 1] /
                       (split_change_point - change_point[nb_segment][segment_index - 1]);
                buff = begin_sum_square[i][change_point[nb_segment][segment_index - 1] - 1] +
                       begin_sum_square[i][split_change_point - 1] +
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
            if (model_type[i - 1] == MULTINOMIAL_CHANGE) {
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
                buff = end_sum_square[i][split_change_point] +
                       end_sum_square[i][change_point[nb_segment][segment_index]];
              }

              else {
                diff = end_sum[i][split_change_point] /
                      (change_point[nb_segment][segment_index] - split_change_point) -
                      end_sum[i][change_point[nb_segment][segment_index]] /
                      (change_point[nb_segment][segment_index + 1] - change_point[nb_segment][segment_index]);
                buff = end_sum_square[i][split_change_point] +
                       end_sum_square[i][change_point[nb_segment][segment_index]] +
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
        os << nb_segment << " " << SEQ_label[SEQL_SEGMENTS] << ":";
        os << " " << likelihood[nb_segment] << " " << penalty[nb_segment]
           << " " << penalized_likelihood[nb_segment] << " ||";
        os << " " << segment_index << " " << split_change_point
           << " | " << begin_change_point << " " << change_point[nb_segment][segment_index + 1]
           << " | " << end_change_point << " " << change_point[nb_segment][segment_index - 1] << " ||";
        for (i = 0;i <= nb_segment;i++) {
          os << " " << change_point[nb_segment][i];
        }
        os << endl;
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
      if (model_type[i - 1] == MULTINOMIAL_CHANGE) {
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

      delete [] begin_sum_square[i];
      delete [] end_sum_square[i];

      delete [] factorial[i];
    }
    delete [] begin_frequency;
    delete [] end_frequency;

    delete [] begin_sum;
    delete [] end_sum;

    delete [] begin_factorial_sum;
    delete [] end_factorial_sum;

    delete [] begin_sum_square;
    delete [] end_sum_square;

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
