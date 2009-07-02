/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       AMAPmod: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2002 UMR Cirad/Inra Modelisation des Plantes
 *
 *       File author(s): Y. Guedon (yann.guedon@cirad.fr)
 *
 *       $Source$
 *       $Id$
 *
 *       Forum for AMAPmod developers: amldevlp@cirad.fr
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
#include "stat_tool/stat_tools.h"
#include "stat_tool/vectors.h"
#include "stat_tool/curves.h"
#include "stat_tool/markovian.h"
#include "stat_tool/stat_label.h"
#include "sequences.h"
#include "sequence_label.h"
#include "tool/config.h"

using namespace std;


extern int column_width(int value);
extern int column_width(int min_value , int max_value);
extern int column_width(int nb_value , const double *value , double scale = 1.);

extern void cumul_computation(int nb_value , const double *pmass , double *pcumul);
extern int cumul_method(int nb_value , const double *cumul , double scale = 1.);

extern char* label(const char *file_name);


#if defined (SYSTEM_IS__CYGWIN)
#define expl exp
#endif



/*--------------------------------------------------------------*/
/**
 *  Ecriture des profils de segments/d'etats pour une sequence.
 *
 *  arguments : stream, indice de la sequence, nombre de segments/d'etats,
 *              pointeur sur les profils de segments/d'etats, label,
 *              pointeur sur les moyennes et sur les profils de ruptures.
 */
/*--------------------------------------------------------------*/

ostream& Sequences::profile_ascii_print(ostream &os , int index , int nb_segment ,
                                        double **profiles , const char *label ,
                                        double **mean , double **change_point) const

{
  register int i , j;
  int start , buff , *width;
  long old_adjust;


  old_adjust = os.flags(ios::adjustfield);

  // calcul des largeurs des colonnes

  width = new int[2 * nb_variable + 4];

  start = 0;
  if (change_point) {
    start++;
  }
  for (i = start;i < nb_variable;i++) {
    if (type[i] != REAL_VALUE) {
      width[i] = column_width((int)min_value[i] , (int)max_value[i]);
    }
    else {
      width[i] = column_width(length[index] , real_sequence[index][i]);
    }
    if (i > start) {
      width[i] += ASCII_SPACE;
    }
  }

  if (index_parameter) {
    width[nb_variable] = column_width(hindex_parameter->nb_value - 1) + ASCII_SPACE;
  }
  else {
    width[nb_variable] = column_width(max_length) + ASCII_SPACE;
  }

  width[nb_variable + 1] = 0;
  for (i = 0;i < length[index];i++) {
    buff = column_width(nb_segment , profiles[i]);
    if (buff > width[nb_variable + 1]) {
      width[nb_variable + 1] = buff;
    }
  }
  width[nb_variable + 1] += ASCII_SPACE;

  width[nb_variable + 2] = column_width(nb_sequence);

  if (mean) {
    for (i = 1;i < nb_variable;i++) {
      if (mean[i]) {
        width[nb_variable + 2 + i] = column_width(length[index] , mean[i]) + ASCII_SPACE;
      }
    }
  }

  if (change_point) {
    width[2 * nb_variable + 2] = 0;
    for (i = 1;i < nb_segment;i++) {
      buff = column_width(length[index] , change_point[i]);
      if (buff > width[2 * nb_variable + 2]) {
        width[2 * nb_variable + 2] = buff;
      }
    }
    width[2 * nb_variable + 2] += ASCII_SPACE;
  }

  else {
    os << SEQ_label[SEQL_OPTIMAL] << " " << label << " | ";
  }
  for (i = 1;i < nb_variable;i++) {
    if ((mean) && (mean[i])) {
      os << SEQ_label[SEQL_SEGMENT_MEAN] << " " << i << " | ";
    }
    os << STAT_label[STATL_VARIABLE] << " " << i << " | ";
  }
  if (index_parameter_type == TIME) {
    os << SEQ_label[SEQL_TIME];
  }
  else {
    os << SEQ_label[SEQL_INDEX];
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

  for (i = 0;i < length[index];i++) {
    os.setf(ios::right , ios::adjustfield);
    if (!change_point) {
      os << setw(width[0]) << int_sequence[index][0][i];
    }

    for (j = 1;j < nb_variable;j++) {
      if ((mean) && (mean[j])) {
        os << setw(width[nb_variable + 2 + j]) << mean[j][i];
      }
      if (type[j] != REAL_VALUE) {
        os << setw(width[j]) << int_sequence[index][j][i];
      }
      else {
        os << setw(width[j]) << real_sequence[index][j][i];
      }
    }
    os << setw(width[nb_variable]) << (index_parameter ? index_parameter[index][i] : i) << "  ";

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
      os << setw(width[nb_variable + 2]) << identifier[index];
    }
    os << endl;
  }

  delete [] width;

  os.setf((FMTFLAGS)old_adjust , ios::adjustfield);

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture des profils de segments/d'etats pour une sequence
 *  au format tableur.
 *
 *  arguments : stream, indice de la sequence, nombre de segments/d'etats,
 *              pointeur sur les profils de segments/d'etats, label,
 *              pointeur sur les moyennes et sur les profils de ruptures.
 *
 *--------------------------------------------------------------*/

ostream& Sequences::profile_spreadsheet_print(ostream &os , int index , int nb_segment ,
                                              double **profiles , const char *label ,
                                              double **mean , double **change_point) const

{
  register int i , j;


  if (!change_point) {
    os << SEQ_label[SEQL_OPTIMAL] << " " << label << "\t";
  }
  for (i = 1;i < nb_variable;i++) {
    if ((mean) && (mean[i])) {
      os << SEQ_label[SEQL_SEGMENT_MEAN] << " " << i << "\t";
    }
    os << STAT_label[STATL_VARIABLE] << " " << i << "\t";
  }
  if (index_parameter_type == TIME) {
    os << SEQ_label[SEQL_TIME];
  }
  else {
    os << SEQ_label[SEQL_INDEX];
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

  for (i = 0;i < length[index];i++) {
    if (!change_point) {
      os << int_sequence[index][0][i] << "\t";
    }
    for (j = 1;j < nb_variable;j++) {
      if ((mean) && (mean[j])) {
        os << mean[j][i] << "\t";
      }
      if (type[j] != REAL_VALUE) {
        os << int_sequence[index][j][i] << "\t";
      }
      else {
        os << real_sequence[index][j][i] << "\t";
      }
    }
    os << (index_parameter ? index_parameter[index][i] : i);
    for (j = 0;j < nb_segment;j++) {
      os << "\t" << profiles[i][j];
    }

    if (change_point) {
      os << "\t";
      for (j = 1;j < nb_segment;j++) {
        os << "\t" << change_point[j][i];
      }
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
 *  Ecriture des profils de segments/d'etats pour une sequence
 *  au format Gnuplot.
 *
 *  arguments : stream, indice de la sequence, nombre de segments/d'etats,
 *              pointeur sur les profils de segments/d'etats, pointeur sur les moyennes et
 *              sur les profils de rupture.
 *
 *--------------------------------------------------------------*/

ostream& Sequences::profile_plot_print(ostream &os , int index , int nb_segment ,
                                       double **profiles , double **mean ,
                                       double **change_point) const

{
  register int i , j;


  for (i = 0;i < length[index];i++) {
    if (index_parameter) {
      os << index_parameter[index][i] << " ";
    }
    for (j = 1;j < nb_variable;j++) {
      if ((mean) && (mean[j])) {
        if (type[j] != REAL_VALUE) {
          os << int_sequence[index][j][i];
        }
        else {
          os << real_sequence[index][j][i];
        }
        os << " " << mean[j][i] << " ";
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
    os << endl;
  }

  return os;
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
  int nb_parameter , *psegment , *pisequence;


  used_output = new bool[NB_OUTPUT];

  nb_parameter = nb_segment - 1;
//  nb_parameter = 0;

  for (i = 1;i < nb_variable;i++) {
    if (model_type[i - 1] == MULTINOMIAL_CHANGE) {
      psegment = int_sequence[index][0] + 1;
      pisequence = int_sequence[index][i];

      for (j = 0;j < marginal[i]->nb_value;j++) {
        used_output[j] = false;
      }
      used_output[*pisequence++] = true;

      for (j = 1;j < length[index];j++) {
        if (*psegment != *(psegment - 1)) {
          for (k = 0;k < marginal[i]->nb_value;k++) {
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
  double sum , factorial_sum , diff , residual , likelihood , *prsequence;


  max_nb_value = 0;
  for (i = 1;i < nb_variable;i++) {
    if (model_type[i - 1] == MULTINOMIAL_CHANGE) {
      if (marginal[i]->nb_value > max_nb_value) {
        max_nb_value = marginal[i]->nb_value;
      }
    }
  }

  if (max_nb_value > 0) {
    frequency = new int[max_nb_value];
  }
  else {
    frequency = 0;
  }

  if ((model_type[0] == MEAN_CHANGE) || (model_type[0] == MEAN_VARIANCE_CHANGE)) {
    residual = 0.;
  }
  else {
    likelihood = 0.;
  }

  for (i = 1;i < nb_variable;i++) {
    if (model_type[i - 1] == MULTINOMIAL_CHANGE) {
      for (j = 0;j < marginal[i]->nb_value;j++) {
        frequency[j] = 0;
      }

      pisequence = int_sequence[index][i];
      for (j = 0;j < length[index];j++) {
        frequency[*pisequence++]++;
      }

      for (j = 0;j < marginal[i]->nb_value;j++) {
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
            double mean , diff , residual2;


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
          likelihood -= ((double)length[index] / 2.) * (log(residual / length[index]) +
                         log(2 * M_PI) + 1);
/*          likelihood -= ((double)length[index] / 2.) * (log(residual / (length[index] - 1)) +
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
                    (log(residual / ((nb_variable - 1) * length[index])) + log(2 * M_PI) + 1);
/*      likelihood = -((double)((nb_variable - 1) * length[index]) / 2.) *
                    (log(residual / ((nb_variable - 1) * (length[index] - 1))) +
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
                                          int output , int* ichange_point)

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
      mean[i] = 0;
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
        variance[i] = 0;
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
      segment_variance = 0;
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
          int_sequence[i][j] = 0;
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
          int_sequence[i][j] = 0;
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

      if (seq->marginal[i]) {
        delete seq->marginal[i];
        seq->marginal[i] = 0;
      }
    }
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Segmentation d'une sequence.
 *
 *  arguments : reference sur un objet Format_error, stream, identificateur de la sequence,
 *              nombre de segments, ruptures, types des modeles,
 *              sortie (sequence ou residus).
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::segmentation(Format_error &error , ostream &os , int iidentifier ,
                                   int nb_segment , int* ichange_point , int *model_type ,
                                   int output) const

{
  bool status = true;
  register int i , j , k , m;
  int index , max_nb_value , nb_parameter , *change_point = 0 , *pindex_param ,
      *frequency , *pisequence , *psegment , inb_segment[1];
  double sum , factorial_sum , diff , residual , segmentation_likelihood ,
         segment_penalty , penalized_likelihood , *mean , *prsequence , *segment_variance , **rank;
  Sequences *iseq , *seq , *oseq;


  oseq = 0;
  error.init();

  if (((index_parameter_type == TIME) && (index_interval->variance > 0.)) ||
      (index_parameter_type == POSITION)) {
    status = false;
    error.update(SEQ_error[SEQR_INDEX_PARAMETER_TYPE]);
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

        if (!marginal[i]) {
          status = false;
          ostringstream error_message;
          error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                        << STAT_error[STATR_MARGINAL_HISTOGRAM];
          error.update((error_message.str()).c_str());
        }

        else if (model_type[i] == MULTINOMIAL_CHANGE) {
          if ((marginal[i]->nb_value < 2) || (marginal[i]->nb_value > NB_OUTPUT)) {
            status = false;
            ostringstream error_message;
            error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                          << STAT_error[STATR_NB_VALUE];
            error.update((error_message.str()).c_str());
          }

          else {
            for (j = 0;j < marginal[i]->nb_value;j++) {
              if (marginal[i]->frequency[j] == 0) {
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
        rank[i] = marginal[i]->rank_computation();
      }
      else {
        rank[i] = 0;
      }
    }

    max_nb_value = 0;

    for (i = 0;i < nb_variable;i++) {
      if ((model_type[i] == MULTINOMIAL_CHANGE) && (marginal[i]->nb_value > max_nb_value)) {
        max_nb_value = marginal[i]->nb_value;
      }
    }

    if (max_nb_value > 0) {
      frequency = new int[max_nb_value];
    }
    else {
      frequency = 0;
    }

    mean = new double[nb_variable];

    if (model_type[0] == MEAN_VARIANCE_CHANGE) {
      segment_variance = new double[nb_segment];
    }
    else {
      segment_variance = 0;
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

          for (k = 0;k < marginal[j]->nb_value;k++) {
            frequency[k] = 0;
          }

          pisequence = int_sequence[index][j] + change_point[i];
          for (k = change_point[i];k < change_point[i + 1];k++) {
            frequency[*pisequence++]++;
          }

          for (k = 0;k < marginal[j]->nb_value;k++) {
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
              segmentation_likelihood -= ((double)(change_point[i + 1] - change_point[i]) / 2.) * (log(residual /
                                           (change_point[i + 1] - change_point[i])) + log(2 * M_PI) + 1);
/*              segmentation_likelihood -= ((double)(change_point[i + 1] - change_point[i]) / 2.) * (log(residual /
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
                                     (log(residual / (nb_variable * (change_point[i + 1] - change_point[i]))) +
                                      log(2 * M_PI) + 1);
/*          segmentation_likelihood -= ((double)(nb_variable * (change_point[i + 1] - change_point[i])) / 2.) *
                                     (log(residual / (nb_variable * (change_point[i + 1] - change_point[i]))) +
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
                                    (log(residual / (nb_variable * length[index])) + log(2 * M_PI) + 1);
/*        segmentation_likelihood = -((double)(nb_variable * length[index]) / 2.) *
                                    (log(residual / (nb_variable * (length[index] - nb_segment))) + log(2 * M_PI)) -
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
    seq->build_marginal_histogram(0);

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
  double sum , factorial_sum , diff , sum_square , buff , segmentation_likelihood ,
         *mean , *residual , *prsequence , *contrast , **factorial , **forward;


  max_nb_segment = nb_segment[0];
  for (i = 1;i < nb_sequence;i++) {
    if (nb_segment[i] > max_nb_segment) {
      max_nb_segment = nb_segment[i];
    }
  }

  max_nb_value = 0;
  factorial = new double*[nb_variable];

  for (i = 1;i < nb_variable;i++) {
    if ((model_type[i - 1] == MULTINOMIAL_CHANGE) && (marginal[i]->nb_value > max_nb_value)) {
      max_nb_value = marginal[i]->nb_value;
    }

    if (model_type[i - 1] == POISSON_CHANGE) {
      factorial[i] = new double[max_length];
    }
    else {
      factorial[i] = 0;
    }
  }

  if (max_nb_value > 0) {
    frequency = new int[max_nb_value];
  }
  else {
    frequency = 0;
  }

  mean = new double[nb_variable];
  residual = new double[max_length];

  contrast = new double[max_length];

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
          for (m = 0;m < marginal[k]->nb_value;m++) {
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
              for (n = 0;n < marginal[k]->nb_value;n++) {
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
                contrast[m] -= ((double)(j - m + 1) / 2.) * (log(residual[m] /
                                 (j - m + 1)) + log(2 * M_PI) + 1);
/*                contrast[m] -= ((double)(j - m + 1) / 2.) * (log(residual[m] /
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
            contrast[k] = -((double)((nb_variable - 1) * (j - k + 1)) / 2.) * (log(contrast[k] /
                             ((nb_variable - 1) * (j - k + 1))) + log(2 * M_PI) + 1);
/*            contrast[k] = -((double)((nb_variable - 1) * (j - k + 1)) / 2.) * (log(contrast[k] /
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
              for (r = 0;r < marginal[k]->nb_value;r++) {
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
  delete marginal[0];
  build_marginal_histogram(0);

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
 *  arguments : reference sur un objet Format_error, stream, nombres de segments,
 *              types des modeles, identificateur de la sequence,
 *              sortie (sequence ou residus).
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::segmentation(Format_error &error , ostream &os , int *nb_segment ,
                                   int *model_type , int iidentifier , int output) const

{
  bool status = true;
  register int i , j;
  int index = I_DEFAULT , nb_parameter , *psegment;
  double segmentation_likelihood , segment_penalty , penalized_likelihood , **rank;
  Sequences *seq , *iseq , *oseq;


  oseq = 0;
  error.init();

  if (((index_parameter_type == TIME) && (index_interval->variance > 0.)) ||
      (index_parameter_type == POSITION)) {
    status = false;
    error.update(SEQ_error[SEQR_INDEX_PARAMETER_TYPE]);
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

        if (!marginal[i]) {
          status = false;
          ostringstream error_message;
          error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                        << STAT_error[STATR_MARGINAL_HISTOGRAM];
          error.update((error_message.str()).c_str());
        }

        else if (model_type[i] == MULTINOMIAL_CHANGE) {
          if ((marginal[i]->nb_value < 2) || (marginal[i]->nb_value > NB_OUTPUT)) {
            status = false;
            ostringstream error_message;
            error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                          << STAT_error[STATR_NB_VALUE];
            error.update((error_message.str()).c_str());
          }

          else {
            for (j = 0;j < marginal[i]->nb_value;j++) {
              if (marginal[i]->frequency[j] == 0) {
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
        rank[i] = seq->marginal[i]->rank_computation();
      }
      else {
        rank[i] = 0;
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
      oseq = 0;
      error.update(SEQ_error[SEQR_SEGMENTATION_FAILURE]);
    }
  }

  return oseq;
}


/*--------------------------------------------------------------*
 *
 *  Segmentation optimale d'une sequence.
 *
 *  arguments : reference sur un objet Format_error, stream,
 *              identificateur de la sequence, nombre de segments maximum,
 *              types des modeles.
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::segmentation(Format_error &error , ostream &os , int iidentifier ,
                                   int max_nb_segment , int *model_type) const

{
  bool status = true;
  register int i , j;
  int index , nb_segment , *nb_parameter , inb_segment[1];
  double max_likelihood[2] , *segmentation_likelihood , *segment_penalty , *slope_penalty ,
         **penalized_likelihood , *likelihood , *change_point_entropy , *segment_entropy , **rank;
  Sequences *iseq , *seq , *oseq;

// double *normalized_change_point_entropy;


  oseq = 0;
  error.init();

  if (((index_parameter_type == TIME) && (index_interval->variance > 0.)) ||
      (index_parameter_type == POSITION)) {
    status = false;
    error.update(SEQ_error[SEQR_INDEX_PARAMETER_TYPE]);
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

        if (!marginal[i]) {
          status = false;
          ostringstream error_message;
          error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                        << STAT_error[STATR_MARGINAL_HISTOGRAM];
          error.update((error_message.str()).c_str());
        }

        else if (model_type[i] == MULTINOMIAL_CHANGE) {
          if ((marginal[i]->nb_value < 2) || (marginal[i]->nb_value > NB_OUTPUT)) {
            status = false;
            ostringstream error_message;
            error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                          << STAT_error[STATR_NB_VALUE];
            error.update((error_message.str()).c_str());
          }

          else {
            for (j = 0;j < marginal[i]->nb_value;j++) {
              if (marginal[i]->frequency[j] == 0) {
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
        rank[i] = marginal[i - 1]->rank_computation();
      }
      else {
        rank[i] = 0;
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
      change_point_entropy = new double[max_nb_segment + 1];
//      normalized_change_point_entropy = new double[max_nb_segment + 1];
      segment_entropy = new double[max_nb_segment + 1];
    }

    inb_segment[0] = max_nb_segment + 1;

    seq->segmentation(inb_segment , model_type , rank , segmentation_likelihood + 1 ,
                      nb_parameter + 1 , segment_penalty + 1);
    if (model_type[0] != MEAN_CHANGE) {
      seq->forward_backward(0 , max_nb_segment , model_type , rank , 0 , SEGMENT , 'a' ,
                            likelihood + 1 , change_point_entropy + 1 , segment_entropy + 1);
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

/*      if (model_type[0] != MEAN_CHANGE) {
        normalized_change_point_entropy[i] = change_point_entropy[i] / i;
      } */
    }

    if (nb_segment > 0) {

#     ifdef MESSAGE
      int width[13];
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
        width[4] = column_width(max_nb_segment - 1 , change_point_entropy + 2) + ASCII_SPACE;
//        width[5] = column_width(max_nb_segment - 1 , normalized_change_point_entropy + 2) + ASCII_SPACE;
        width[6] = column_width(max_nb_segment - 1 , segment_entropy + 2) + ASCII_SPACE;
        width[7] = column_width(nb_parameter[max_nb_segment]) + ASCII_SPACE;
        width[8] = column_width(max_nb_segment , penalized_likelihood[0] + 1) + ASCII_SPACE;
        width[9] = column_width(max_nb_segment , weight[0] + 1) + ASCII_SPACE;
        width[10] = column_width(max_nb_segment , penalized_likelihood[1] + 1) + ASCII_SPACE;
        width[11] = column_width(max_nb_segment , weight[1] + 1) + ASCII_SPACE;
        width[12] = column_width(max_nb_segment - 1 , curvature + 2) + ASCII_SPACE;

        os << "\n" << SEQ_label[SEQL_NB_SEGMENT] << " |  2 * " << STAT_label[STATL_LIKELIHOOD]
           << " |  2 * " << SEQ_label[SEQL_POSSIBLE_SEGMENTATION_LIKELIHOOD]
           << " | " << SEQ_label[SEQL_POSTERIOR_PROBABILITY] << " | " << SEQ_label[SEQL_CHANGE_POINT_ENTROPY]
//           << " | " << SEQ_label[SEQL_NORMALIZED_CHANGE_POINT_ENTROPY]
           << " | " << SEQ_label[SEQL_SEGMENT_ENTROPY] << " | " << STAT_label[STATL_FREE_PARAMETERS]
           << " | "  << STAT_criterion_word[BIC] << " - "  << STAT_label[STATL_WEIGHT]
           << " | Modified " << STAT_criterion_word[BIC] << " - "  << STAT_label[STATL_WEIGHT]
           << " | Lavielle criterion (standardized curvature)" << endl;

        os.setf(ios::left , ios::adjustfield);

        os << setw(width[0]) << 1
           << setw(width[1]) << 2 * segmentation_likelihood[1]
           << setw(width[2]) << 2 * likelihood[1]
           << setw(width[3]) << posterior_probability[1]
           << setw(width[4]) << " "
           << setw(width[6]) << " "
           << setw(width[7]) << nb_parameter[1]
           << setw(width[8]) << penalized_likelihood[0][1]
           << setw(width[9]) << weight[0][1]
           << setw(width[10]) << penalized_likelihood[1][1]
           << setw(width[11]) << weight[1][1]
           << setw(width[12]) << " " << endl;

        for (i = 2;i <= max_nb_segment;i++) {
          os << setw(width[0]) << i
             << setw(width[1]) << 2 * segmentation_likelihood[i]
             << setw(width[2]) << 2 * likelihood[i]
             << setw(width[3]) << posterior_probability[i]
             << setw(width[4]) << change_point_entropy[i]
//             << setw(width[5]) << normalized_change_point_entropy[i]
             << setw(width[6]) << segment_entropy[i]
             << setw(width[7]) << nb_parameter[i]
             << setw(width[8]) << penalized_likelihood[0][i]
             << setw(width[9]) << weight[0][i]
             << setw(width[10]) << penalized_likelihood[1][i]
             << setw(width[11]) << weight[1][i]
             << setw(width[12]) << curvature[i] << endl;
        }
        os << endl;
      }

      else {
        width[0] = column_width(max_nb_segment) + ASCII_SPACE;
        width[1] = column_width(max_nb_segment , segmentation_likelihood + 1 , 2.) + ASCII_SPACE;
        width[7] = column_width(nb_parameter[max_nb_segment]) + ASCII_SPACE;
        width[8] = column_width(max_nb_segment , penalized_likelihood[0] + 1) + ASCII_SPACE;
        width[9] = column_width(max_nb_segment , weight[0] + 1) + ASCII_SPACE;
        width[10] = column_width(max_nb_segment , penalized_likelihood[1] + 1) + ASCII_SPACE;
        width[11] = column_width(max_nb_segment , weight[1] + 1) + ASCII_SPACE;
        width[12] = column_width(max_nb_segment - 1 , curvature + 2) + ASCII_SPACE;

        os << "\n" << SEQ_label[SEQL_NB_SEGMENT] << " |  2 * " << STAT_label[STATL_LIKELIHOOD]
           << " | " << STAT_label[STATL_FREE_PARAMETERS]
           << " | "  << STAT_criterion_word[BIC] << " - "  << STAT_label[STATL_WEIGHT]
           << " | Modified " << STAT_criterion_word[BIC] << " - "  << STAT_label[STATL_WEIGHT]
           << " | Lavielle criterion (standardized curvature)" << endl;

        os.setf(ios::left , ios::adjustfield);

        os << setw(width[0]) << 1
           << setw(width[1]) << 2 * segmentation_likelihood[1]
           << setw(width[7]) << nb_parameter[1]
           << setw(width[8]) << penalized_likelihood[0][1]
           << setw(width[9]) << weight[0][1]
           << setw(width[10]) << penalized_likelihood[1][1]
           << setw(width[11]) << weight[1][1]
           << setw(width[12]) << " " << endl;

        for (i = 2;i <= max_nb_segment;i++) {
          os << setw(width[0]) << i
             << setw(width[1]) << 2 * segmentation_likelihood[i]
             << setw(width[7]) << nb_parameter[i]
             << setw(width[8]) << penalized_likelihood[0][i]
             << setw(width[9]) << weight[0][i]
             << setw(width[10]) << penalized_likelihood[1][i]
             << setw(width[11]) << weight[1][i]
             << setw(width[12]) << curvature[i] << endl;
        }
        os << endl;
      }

      os.setf((FMTFLAGS)old_adjust , ios::adjustfield);

#     ifdef DEBUG
      for (i = 1;i <= max_nb_segment;i++) {
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
      }
#     endif

      for (i = 1;i <= max_nb_segment;i++) {
        os << "\n" << i << " " << (i == 1 ? SEQ_label[SEQL_SEGMENT] : SEQ_label[SEQL_SEGMENTS])
           << "   2 * " << STAT_label[STATL_LIKELIHOOD] << ": " << 2 * segmentation_likelihood[i] << "   "
           << nb_parameter[i] << " " << STAT_label[nb_parameter[i] == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS]
           << "   2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " (Modified "  << STAT_criterion_word[BIC] << "): "
           << penalized_likelihood[1][i] << "   " << STAT_label[STATL_WEIGHT] << ": " << weight[1][i] << endl;
      }

/*      for (i = 1;i <= max_nb_segment;i++) {
        os << "\n" << i << " " << (i == 1 ? SEQ_label[SEQL_SEGMENT] : SEQ_label[SEQL_SEGMENTS])
           << "   2 * " << STAT_label[STATL_LIKELIHOOD] << ": " << 2 * segmentation_likelihood[i] << "   "
           << nb_parameter[i] << " " << STAT_label[nb_parameter[i] == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS]
           << "   2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[BIC] << "): "
           << penalized_likelihood[0][i] << "   " << STAT_label[STATL_WEIGHT] << ": " << weight[0][i] << endl;
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
      } */

      os << endl;

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
        seq->build_marginal_histogram(0);
      }

      else {
        seq->segmentation(inb_segment , model_type , rank);
      }

      oseq = seq->segmentation_output(inb_segment , model_type , os);

#     ifdef MESSAGE
      if ((model_type[0] != MEAN_CHANGE) && (model_type[0] != MEAN_VARIANCE_CHANGE)) {
        hierarchical_segmentation(error , os , iidentifier , max_nb_segment , model_type);
      }
#     endif

    }

    else {
      oseq = 0;
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
      delete [] change_point_entropy;
//      delete [] normalized_change_point_entropy;
      delete [] segment_entropy;
    }
  }

  return oseq;
}


/*--------------------------------------------------------------*
 *
 *  Calcul par sommation des profils de segmentation d'une sequence.
 *
 *  arguments : indice de la sequence, nombre de segments, types des modeles,
 *              rangs (variables ordinales), stream, type de sortie, format
 *              de fichier ('a' : ASCII, 's' : Spreadsheet, 'g' : Gnuplot),
 *              pointeurs sur les vraisemblances de toutes les segmentations possibles,
 *              les entropies des ruptures et les entropies marginale.
 *
 *--------------------------------------------------------------*/

double Sequences::forward_backward(int index , int nb_segment , int *model_type ,
                                   double **rank , ostream *os , int output , char format ,
                                   double *ilikelihood , double *ichange_point_entropy ,
                                   double *isegment_entropy) const

{
  register int i , j , k;
  int max_nb_value , *frequency , *pisequence;
  double sum , factorial_sum , diff , sum_square , buff , segment_norm , sequence_norm ,
         rlikelihood , change_point_entropy , segment_entropy , *mean , *residual , *prsequence ,
         *contrast , *likelihood , *norm , *forward_norm , *backward_norm , *entropy_backward ,
         **factorial , **forward , **backward , **backward_output , **change_point;
//  long double **forward , **backward;


  max_nb_value = 0;
  factorial = new double*[nb_variable];

  for (i = 1;i < nb_variable;i++) {
    if ((model_type[i - 1] == MULTINOMIAL_CHANGE) && (marginal[i]->nb_value > max_nb_value)) {
      max_nb_value = marginal[i]->nb_value;
    }

    if (model_type[i - 1] == POISSON_CHANGE) {
      factorial[i] = new double[length[index]];
    }
    else {
      factorial[i] = 0;
    }
  }

  if (max_nb_value > 0) {
    frequency = new int[max_nb_value];
  }
  else {
    frequency = 0;
  }

  mean = new double[nb_variable];
  residual = new double[length[index]];

  contrast = new double[length[index]];

  forward = new double*[length[index]];
  for (i = 0;i < length[index];i++) {
    forward[i] = new double[nb_segment];
  }

  norm = new double[length[index]];
  forward_norm = new double[length[index]];

  if (ilikelihood) {
    likelihood = ilikelihood;
  }
  else {
    likelihood = new double[nb_segment];
  }

  backward = new double*[length[index]];
  for (i = 0;i < length[index];i++) {
    backward[i] = new double[nb_segment];
  }

  backward_norm = new double[length[index]];

/*  forward = new long double*[length[index]];
  for (i = 0;i < length[index];i++) {
    forward[i] = new long double[nb_segment];
  }

  backward = new long double*[length[index]];
  for (i = 0;i < length[index];i++) {
    backward[i] = new long double[nb_segment];
  } */

  backward_output = new double*[length[index]];
  for (i = 0;i < length[index];i++) {
    backward_output[i] = new double[nb_segment];
  }

  change_point = new double*[nb_segment];
  for (i = 1;i < nb_segment;i++) {
    change_point[i] = new double[length[index]];
  }

  if (isegment_entropy) {
    entropy_backward = new double[nb_segment];
  }

  for (i = 1;i < nb_variable;i++) {
    if (model_type[i - 1] == VARIANCE_CHANGE) {
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

  // recurrence "forward"

  for (i = 0;i < length[index];i++) {

    // calcul des log-vraisemblances des segments

    for (j = 0;j <= i;j++) {
      contrast[j] = 0.;
    }

    for (j = 1;j < nb_variable;j++) {
      if (model_type[j - 1] == MULTINOMIAL_CHANGE) {
        for (k = 0;k < marginal[j]->nb_value;k++) {
          frequency[k] = 0;
        }
        sum = 0.;

        pisequence = int_sequence[index][j] + i;
        frequency[*pisequence--]++;
        for (k = i - 1;k >= 0;k--) {
          sum += (i - k) * log((double)(i - k) / (double)(i - k + 1)) +
                 log((double)(frequency[*pisequence] + 1) / (double)(i - k + 1));
          if (frequency[*pisequence] > 0) {
            sum -= frequency[*pisequence] *
                   log((double)frequency[*pisequence] / (double)(frequency[*pisequence] + 1));
          }
          frequency[*pisequence--]++;

          if (contrast[k] != D_INF) {
            contrast[k] += sum;
          }

/*          frequency[*pisequence--]++;
          if (contrast[k] != D_INF) {
            for (m = 0;m < marginal[j]->nb_value;m++) {
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

        pisequence = int_sequence[index][j] + i;
        for (k = i;k >= 0;k--) {
          sum += *pisequence--;
          factorial_sum += factorial[j][k];
          if ((contrast[k] != D_INF) && (sum > 0.)) {
            contrast[k] += sum * (log(sum / (i - k + 1)) - 1) - factorial_sum;
          }
        }
      }

      else {
        if (model_type[j - 1] == VARIANCE_CHANGE) {
          sum_square = 0.;

          if (type[j] != REAL_VALUE) {
            pisequence = int_sequence[index][j] + i;
            for (k = i;k >= 0;k--) {
              diff = *pisequence-- - mean[j];
              sum_square += diff * diff;
              residual[k] = sum_square;
            }
          }

          else {
            prsequence = real_sequence[index][j] + i;
            for (k = i;k >= 0;k--) {
              diff = *prsequence-- - mean[j];
              sum_square += diff * diff;
              residual[k] = sum_square;
            }
          }
        }

        else if (model_type[j - 1] == ORDINAL_GAUSSIAN_CHANGE) {
          pisequence = int_sequence[index][j] + i;
          sum_square = 0.;
          sum = rank[j][*pisequence--];
          residual[i] = 0.;

          for (k = i - 1;k >= 0;k--) {
            diff = rank[j][*pisequence] - sum / (i - k);
            sum_square += ((double)(i - k) / (double)(i - k + 1)) * diff * diff;
            sum += rank[j][*pisequence--];
            residual[k] = sum_square;
          }
        }

        else {
          if (type[j] != REAL_VALUE) {
            pisequence = int_sequence[index][j] + i;
            sum_square = 0.;
            sum = *pisequence--;
            residual[i] = 0.;

            for (k = i - 1;k >= 0;k--) {
              diff = *pisequence - sum / (i - k);
              sum_square += ((double)(i - k) / (double)(i - k + 1)) * diff * diff;
              sum += *pisequence--;
              residual[k] = sum_square;
            }
          }

          else {
            prsequence = real_sequence[index][j] + i;
            sum_square = 0.;
            sum = *prsequence--;
            residual[i] = 0.;

            for (k = i - 1;k >= 0;k--) {
              diff = *prsequence - sum / (i - k);
              sum_square += ((double)(i - k) / (double)(i - k + 1)) * diff * diff;
              sum += *prsequence--;
              residual[k] = sum_square;
            }
          }
        }

        if (model_type[j - 1] == MEAN_VARIANCE_CHANGE) {
          for (k = i - 1;k >= 0;k--) {
            contrast[k] += residual[k];
          }
        }

        else {
          for (k = i;k >= 0;k--) {
//            if ((contrast[k] != D_INF) && (residual[k] > 0.)) {
            if ((contrast[k] != D_INF) && (residual[k] > sqrt((double)(i - k + 1)) * ROUNDOFF_ERROR)) {
              contrast[k] -= ((double)(i - k + 1) / 2.) * (log(residual[k] /
                               (i - k + 1)) + log(2 * M_PI) + 1);
/*              contrast[k] -= ((double)(i - k + 1) / 2.) * (log(residual[k] /
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
//        if (contrast[j] > 0.) {
        if (contrast[j] > sqrt((double)((nb_variable - 1) * (i - j + 1))) * ROUNDOFF_ERROR) {
          contrast[j] = -((double)((nb_variable - 1) * (i - j + 1)) / 2.) * (log(contrast[j] /
                           ((nb_variable - 1) * (i - j + 1))) + log(2 * M_PI) + 1);
/*          contrast[j] = -((double)((nb_variable - 1) * (i - j + 1)) / 2.) * (log(contrast[j] /
                           ((nb_variable - 1) * (i - j))) + log(2 * M_PI)) +
                         (double)((nb_variable - 1) * (i - j)) / 2.; */
        }
        else {
          contrast[j] = D_INF;
        }
      }
    }

/*    for (j = 0;j < nb_segment;j++) {
      forward[i][j] = D_INF;
    }

    for (j = MAX(0 , nb_segment + i - length[index]);j < MIN((i < length[index] - 1 ? nb_segment - 1 : nb_segment) , i + 1);j++) {
      if (j == 0) {
        forward[i][j] = contrast[0];
      }

      else {
        forward[i][j] = 0.;
        for (k = i;k >= j;k--) {
          if ((contrast[k] != D_INF) && (forward[k - 1][j - 1] != D_INF)) {
            forward[i][j] += expl(contrast[k] + forward[k - 1][j - 1]);
          }
        }

        if (forward[i][j] > 0.) {
          forward[i][j] = logl(forward[i][j]);
        }
        else {
          forward[i][j] = D_INF;
        }
      }
    } */

    if (contrast[i] != D_INF) {
      contrast[i] = exp(contrast[i]);
    }
    else {
      contrast[i] = 0.;
    }

    segment_norm = 0.;
    for (j = i - 1;j >= 0;j--) {
      segment_norm += norm[j];

#     ifdef DEBUG
      if (i == length[index] - 1) {
        cout << j << ": " << contrast[j] << " " << segment_norm << " | ";
      }
#     endif

      if (contrast[j] != D_INF) {
        contrast[j] = exp(contrast[j] - segment_norm);
      }
      else {
        contrast[j] = 0.;
      }

#     ifdef DEBUG
      if (i == length[index] - 1) {
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
    if (i == length[index] - 1) {
      cout << endl;
    }
#   endif

    for (j = 0;j < nb_segment;j++) {
      forward[i][j] = 0.;
    }
    norm[i] = 0.;

    for (j = 0;j < MIN((i < length[index] - 1 ? nb_segment - 1 : nb_segment) , i + 1);j++) {
//    for (j = MAX(0 , nb_segment + i - length[index]);j < MIN((i < length[index] - 1 ? nb_segment - 1 : nb_segment) , i + 1);j++) {
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
      for (j = 0;j < MIN((i < length[index] - 1 ? nb_segment - 1 : nb_segment) , i + 1);j++) {
//      for (j = MAX(0 , nb_segment + i - length[index]);j < MIN((i < length[index] - 1 ? nb_segment - 1 : nb_segment) , i + 1);j++) {
        forward[i][j] /= norm[i];
      }

      norm[i] = log(norm[i]);
    }

    forward_norm[i] = segment_norm + norm[i];

#   ifdef DEBUG
    cout << i << " |";
    for (j = 0;j < nb_segment;j++) {
      cout << " " << forward[i][j];
    }
    cout << " | " << exp(norm[i]) << endl;
#   endif

  }

//  rlikelihood = forward[length[index] - 1][nb_segment - 1];

  for (i = 0;i < nb_segment;i++) {
    if (forward[length[index] - 1][i] > 0.) {
      likelihood[i] = log(forward[length[index] - 1][i]) + forward_norm[length[index] - 1];
    }
    else {
      likelihood[i] = D_INF;
    }
  }

  rlikelihood = likelihood[nb_segment - 1];

  if (rlikelihood != D_INF) {

    // recurrence "backward"

    for (i = length[index] - 1;i >= 0;i--) {

      // calcul des log-vraisemblances des segments

      for (j = i;j < length[index];j++) {
        contrast[j] = 0.;
      }

      for (j = 1;j < nb_variable;j++) {
        if (model_type[j - 1] == MULTINOMIAL_CHANGE) {
          for (k = 0;k < marginal[j]->nb_value;k++) {
            frequency[k] = 0;
          }
          sum = 0.;

          pisequence = int_sequence[index][j] + i;
          frequency[*pisequence++]++;
          for (k = i + 1;k < length[index];k++) {
            sum += (k - i) * log((double)(k - i) / (double)(k - i + 1)) +
                   log((double)(frequency[*pisequence] + 1) / (double)(k - i + 1));
            if (frequency[*pisequence] > 0) {
              sum -= frequency[*pisequence] *
                     log((double)frequency[*pisequence] / (double)(frequency[*pisequence] + 1));
            }
            frequency[*pisequence++]++;

            if (contrast[k] != D_INF) {
              contrast[k] += sum;
            }

/*            frequency[*pisequence++]++;
            if (contrast[k] != D_INF) {
              for (m = 0;m < marginal[j]->nb_value;m++) {
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

          pisequence = int_sequence[index][j] + i;
          for (k = i;k < length[index];k++) {
            sum += *pisequence++;
            factorial_sum += factorial[j][k];
            if ((contrast[k] != D_INF) && (sum > 0.)) {
              contrast[k] += sum * (log(sum / (k - i + 1)) - 1) - factorial_sum;
            }
          }
        }

        else {
          if (model_type[j - 1] == VARIANCE_CHANGE) {
            sum_square = 0.;

            if (type[j] != REAL_VALUE) {
              pisequence = int_sequence[index][j] + i;
              for (k = i;k < length[index];k++) {
                diff = *pisequence++ - mean[j];
                sum_square += diff * diff;
                residual[k] = sum_square;
              }
            }

            else {
              prsequence = real_sequence[index][j] + i;
              for (k = i;k < length[index];k++) {
                diff = *prsequence++ - mean[j];
                sum_square += diff * diff;
                residual[k] = sum_square;
              }
            }
          }

          else if (model_type[j - 1] == ORDINAL_GAUSSIAN_CHANGE) {
            pisequence = int_sequence[index][j] + i;
            sum_square = 0.;
            sum = rank[j][*pisequence++];
            residual[i] = 0.;

            for (k = i + 1;k < length[index];k++) {
              diff = rank[j][*pisequence] - sum / (k - i);
              sum_square += ((double)(k - i) / (double)(k - i + 1)) * diff * diff;
              sum += rank[j][*pisequence++];
              residual[k] = sum_square;
            }
          }

          else {
            if (type[j] != REAL_VALUE) {
              pisequence = int_sequence[index][j] + i;
              sum_square = 0.;
              sum = *pisequence++;
              residual[i] = 0.;

              for (k = i + 1;k < length[index];k++) {
                diff = *pisequence - sum / (k - i);
                sum_square += ((double)(k - i) / (double)(k - i + 1)) * diff * diff;
                sum += *pisequence++;
                residual[k] = sum_square;
              }
            }

            else {
              prsequence = real_sequence[index][j] + i;
              sum_square = 0.;
              sum = *prsequence++;
              residual[i] = 0.;

              for (k = i + 1;k < length[index];k++) {
                diff = *prsequence - sum / (k - i);
                sum_square += ((double)(k - i) / (double)(k - i + 1)) * diff * diff;
                sum += *prsequence++;
                residual[k] = sum_square;
              }
            }
          }
 
          if (model_type[j - 1] == MEAN_VARIANCE_CHANGE) {
            for (k = i + 1;k < length[index];k++) {
              contrast[k] += residual[k];
            }
          }

          else {
            for (k = i;k < length[index];k++) {
//              if ((contrast[k] != D_INF) && (residual[k] > 0.)) {
              if ((contrast[k] != D_INF) && (residual[k] > sqrt((double)(k - i + 1)) * ROUNDOFF_ERROR)) {
                contrast[k] -= ((double)(k - i + 1) / 2.) * (log(residual[k] /
                                 (k - i + 1)) + log(2 * M_PI) + 1);
/*                contrast[k] -= ((double)(k - i + 1) / 2.) * (log(residual[k] /
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
//          if (contrast[j] > 0.) {
          if (contrast[j] > sqrt((double)((nb_variable - 1) * (j - i + 1))) * ROUNDOFF_ERROR) {
            contrast[j] = -((double)((nb_variable - 1) * (j - i + 1)) / 2.) * (log(contrast[j] /
                             ((nb_variable - 1) * (j - i + 1))) + log(2 * M_PI) + 1);
/*            contrast[j] = -((double)((nb_variable - 1) * (j - i + 1)) / 2.) * (log(contrast[j] /
                             ((nb_variable - 1) * (j - i))) + log(2 * M_PI)) +
                           (double)((nb_variable - 1) * (j - i)) / 2.; */
          }
          else {
            contrast[j] = D_INF;
          }
        }
      }

/*      for (j = 0;j < nb_segment;j++) {
        backward[i][j] = D_INF;
        backward_output[i][j] = 0.;
      }

      for (j = MAX((i == 0 ? 0 : 1) , nb_segment + i - length[index]);j < MIN(nb_segment , i + 1);j++) {
        if (j < nb_segment - 1) {
          backward[i][j] = 0.;
          for (k = i;k <= length[index] + j - nb_segment;k++) {
            if ((contrast[k] != D_INF) && (backward[k + 1][j + 1] != D_INF)) {
              backward[i][j] += expl(contrast[k] + backward[k + 1][j + 1]);
            }
          }

          if (backward[i][j] > 0.) {
            backward[i][j] = logl(backward[i][j]);
          }
          else {
            backward[i][j] = D_INF;
          }
        }

        else {
          backward[i][j] = contrast[length[index] - 1];
        }

        if (backward[i][j] != D_INF) {
          if ((i == 0) && (j == 0)) {
            buff = 1.;
          }
          else if ((j > 0) && (forward[i - 1][j - 1] != D_INF)) {
            buff = exp(forward[i - 1][j - 1] + backward[i][j] - rlikelihood);
          }

          if (output == CHANGE_POINT) {
            backward_output[i][j] = buff;
          }
          change_point[nb_segment - 1][i] += buff;
        }
      }

      if (output == SEGMENT) {
        if (i < length[index] - 1) {
          for (j = MAX(0 , nb_segment + i - length[index]);j < MIN(nb_segment , i + 1);j++) {
            backward_output[i][j] = backward_output[i + 1][j];
            if ((j < nb_segment - 1) && (forward[i][j] != D_INF) && (backward[i + 1][j + 1] != D_INF)) {
              backward_output[i][j] += exp(forward[i][j] + backward[i + 1][j + 1] - rlikelihood);
            }
            if ((j > 0) && (forward[i][j - 1] != D_INF) && (backward[i + 1][j] != D_INF)) {
              backward_output[i][j] -= exp(forward[i][j - 1] + backward[i + 1][j] - rlikelihood);
            }

            if (backward_output[i][j] < 0.) {
              backward_output[i][j] = 0.;
            }
            if (backward_output[i][j] > 1.) {
              backward_output[i][j] = 1.;
            }
          }
        }

        else {
          backward_output[i][nb_segment - 1] = 1.;
        }
      } */

      if (contrast[i] != D_INF) {
        contrast[i] = exp(contrast[i]);
      }
      else {
        contrast[i] = 0.;
      }

      segment_norm = 0.;
      for (j = i + 1;j < length[index];j++) {
        segment_norm += norm[j];
        if (contrast[j] != D_INF) {
          contrast[j] = exp(contrast[j] - segment_norm);
        }
        else {
          contrast[j] = 0.;
        }
      }

      for (j = 0;j < nb_segment;j++) {
        backward[i][j] = 0.;
        backward_output[i][j] = 0.;
      }
      norm[i] = 0.;

      for (j = MAX((i == 0 ? 0 : 1) , nb_segment + i - length[index]);j < nb_segment;j++) {
//      for (j = MAX((i == 0 ? 0 : 1) , nb_segment + i - length[index]);j < MIN(nb_segment , i + 1);j++) {
        if (j < nb_segment - 1) {
          for (k = i;k <= length[index] + j - nb_segment;k++) {
            backward[i][j] += contrast[k] * backward[k + 1][j + 1];
          }
        }
        else {
          backward[i][j] = contrast[length[index] - 1];
        }

        norm[i] += backward[i][j];
      }

      if (norm[i] > 0.) {
        for (j = MAX((i == 0 ? 0 : 1) , nb_segment + i - length[index]);j < nb_segment;j++) {
//        for (j = MAX((i == 0 ? 0 : 1) , nb_segment + i - length[index]);j < MIN(nb_segment , i + 1);j++) {
          backward[i][j] /= norm[i];
        }

        norm[i] = log(norm[i]);
      }

      backward_norm[i] = segment_norm + norm[i];

      if (output == SEGMENT) {
        if (i < length[index] - 1) {
          for (j = MAX(0 , nb_segment + i - length[index]);j < MIN(nb_segment , i + 1);j++) {
            backward_output[i][j] = backward_output[i + 1][j];
            if (j < nb_segment - 1) {
              backward_output[i][j] += forward[i][j] * backward[i + 1][j + 1] * sequence_norm;
            }
            if (j > 0) {
              backward_output[i][j] -= forward[i][j - 1] * backward[i + 1][j] * sequence_norm;
            }

            if (backward_output[i][j] < 0.) {
              backward_output[i][j] = 0.;
            }
            if (backward_output[i][j] > 1.) {
              backward_output[i][j] = 1.;
            }
          }
        }

        else {
          backward_output[i][nb_segment - 1] = 1.;
        }
      }

      if (i == 0) {

#       ifdef MESSAGE
        buff = backward[i][0] * exp(backward_norm[i] - rlikelihood);
        if ((buff < 1. - DOUBLE_ERROR) || (buff > 1. + DOUBLE_ERROR)) {
          cout << "\nERROR: " << buff << " | " << 1 << endl;
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
        for (j = 1;j < nb_segment - 1;j++) {
          change_point[j][i] = 0.;
          for (k = MAX(1 , j + 1 + i - length[index]);k <= MIN(j , i);k++) {
            change_point[j][i] += forward[i - 1][k - 1] * backward[i][k + nb_segment - j - 1];
          }
          change_point[j][i] *= exp(forward_norm[i - 1] + backward_norm[i] - likelihood[j]);
        }

        sequence_norm = exp(forward_norm[i - 1] + backward_norm[i] - rlikelihood);

#       ifdef DEBUG
        cout << i << ": " << forward_norm[i - 1] << " " << backward_norm[i] << " | "
             << forward_norm[i - 1] + backward_norm[i] - rlikelihood << " " << sequence_norm << endl;
#       endif

        change_point[nb_segment - 1][i] = 0.;
        for (j = MAX(1 , nb_segment + i - length[index]);j < MIN(nb_segment , i + 1);j++) {
          if (output == CHANGE_POINT) {
            backward_output[i][j] = forward[i - 1][j - 1] * backward[i][j] * sequence_norm;
          }
          change_point[nb_segment - 1][i] += forward[i - 1][j - 1] * backward[i][j];
        }
        change_point[nb_segment - 1][i] *= sequence_norm;
      }
    }

#   ifdef MESSAGE
    if (output == SEGMENT) {
      for (i = 0;i < length[index] - 1;i++) {
        sum = 0.;
        for (j = 0;j < nb_segment;j++) {
          sum += backward_output[i][j];
        }
        if ((sum < 1. - DOUBLE_ERROR) || (sum > 1. + DOUBLE_ERROR)) {
          cout << "\nERROR: " << i << " | " << sum << endl;
        }
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

    change_point_entropy = 0.;
    for (i = 1;i < length[index];i++) {
      if ((change_point[nb_segment - 1][i] > 0.) && (change_point[nb_segment - 1][i] < 1.)) {
        change_point_entropy -= change_point[nb_segment - 1][i] * log(change_point[nb_segment - 1][i]) +
                                (1 - change_point[nb_segment - 1][i]) * log(1 - change_point[nb_segment - 1][i]);
      }
    }

    segment_entropy = 0.;
    for (i = 1;i < length[index] - 1;i++) {
      for (j = 0;j < nb_segment;j++) {
        if (backward_output[i][j] > 0.) {
          segment_entropy -= backward_output[i][j] * log(backward_output[i][j]);
        }
      }
    }

    if (ichange_point_entropy) {
      for (i = 1;i < nb_segment;i++) {
        ichange_point_entropy[i] = 0.;
        for (j = 1;j < length[index];j++) {
          if ((change_point[i][j] > 0.) && (change_point[i][j] < 1.)) {
            ichange_point_entropy[i] -= change_point[i][j] * log(change_point[i][j]) +
                                        (1 - change_point[i][j]) * log(1 - change_point[i][j]);
          }
        }
      }

      ichange_point_entropy[nb_segment - 1] = change_point_entropy;
    }

    if (isegment_entropy) {
      for (i = 1;i < nb_segment - 1;i++) {
        isegment_entropy[i] = 0.;
        for (j = 0;j < i;j++) {
          entropy_backward[j] = 0.;
        }
        entropy_backward[i] = 1.;

        for (j = length[index] - 2;j > 0;j--) {
          for (k = MAX(0 , i + 1 + j - length[index]);k <= MIN(i , j);k++) {
            sequence_norm = exp(forward_norm[j] + backward_norm[j + 1] - likelihood[i]);
            if (k < i) {
              entropy_backward[k] += forward[j][k] * backward[j + 1][k + nb_segment - i] * sequence_norm;
            }
            if (k > 0) {
              entropy_backward[k] -= forward[j][k - 1] * backward[j + 1][k + nb_segment - i - 1] * sequence_norm;
            }

            if (entropy_backward[k] < 0.) {
              entropy_backward[k] = 0.;
            }
            if (entropy_backward[k] > 1.) {
              entropy_backward[k] = 1.;
            }

            if (entropy_backward[k] > 0.) {
              isegment_entropy[i] -= entropy_backward[k] * log(entropy_backward[k]);
            }
          }

#         ifdef MESSAGE
          sum = 0.;
          for (k = MAX(0 , i + 1 + j - length[index]);k <= MIN(i , j);k++) {
            sum += entropy_backward[k];
          }
          if ((sum < 1. - DOUBLE_ERROR) || (sum > 1. + DOUBLE_ERROR)) {
            cout << "\nERROR: " << i + 1 << " " << j << " | " << sum << endl;
          }
#         endif

        }
      }

      isegment_entropy[nb_segment - 1] = segment_entropy;
    }

    if (os) {
      switch (format) {

      case 'a' : {
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
                            0 , change_point);

        *os << "\n" << SEQ_label[SEQL_POSSIBLE_SEGMENTATION_LIKELIHOOD] << ": " << rlikelihood << endl;
        *os << "\n" << SEQ_label[SEQL_CHANGE_POINT_ENTROPY] << ": " << change_point_entropy << " ("
            << change_point_entropy / nb_segment << ")";
        if (output == SEGMENT) {
          *os << "\n" << SEQ_label[SEQL_SEGMENT_ENTROPY] << ": " << segment_entropy;
        }
        *os << endl;
        break;
      }

      case 's' : {
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
                                  0 , change_point);

        *os << "\n" << SEQ_label[SEQL_POSSIBLE_SEGMENTATION_LIKELIHOOD] << "\t" << rlikelihood << endl;
        *os << "\n" << SEQ_label[SEQL_CHANGE_POINT_ENTROPY] << "\t" << change_point_entropy << "\t"
            << change_point_entropy / nb_segment;
        if (output == SEGMENT) {
          *os << "\n" << SEQ_label[SEQL_SEGMENT_ENTROPY] << "\t" << segment_entropy;
        }
        *os << endl;
        break;
      }

      case 'g' : {
        profile_plot_print(*os , index , nb_segment , backward_output , 0 , change_point);
        break;
      }
      }
    }
  }

  delete [] frequency;

  for (i = 1;i < nb_variable;i++) {
    delete [] factorial[i];
  }
  delete [] factorial;

  delete [] mean;
  delete [] residual;

  delete [] contrast;

  for (i = 0;i < length[index];i++) {
    delete [] forward[i];
  }
  delete [] forward;

  delete [] norm;
  delete [] forward_norm;

  if (!ilikelihood) {
    delete [] likelihood;
  }

  for (i = 0;i < length[index];i++) {
    delete [] backward[i];
  }
  delete [] backward;

  delete [] backward_norm;

  for (i = 0;i < length[index];i++) {
    delete [] backward_output[i];
  }
  delete [] backward_output;

  for (i = 1;i < nb_segment;i++) {
    delete [] change_point[i];
  }
  delete [] change_point;

  if (isegment_entropy) {
    delete [] entropy_backward;
  }

  return rlikelihood;
}


/*--------------------------------------------------------------*
 *
 *  Simulation des segmentations d'une sequence.
 *
 *  arguments : indice de la sequence, nombre de segments, types des modeles,
 *              rangs (variables ordinales), stream, format de fichier ('a' : ASCII,
 *              's' : Spreadsheet), nombre de segmentation.
 *
 *--------------------------------------------------------------*/

double Sequences::forward_backward_sampling(int index , int nb_segment , int *model_type ,
                                            double **rank , ostream &os , char format ,
                                            int nb_segmentation) const

{
  register int i , j , k , m , n;
  int max_nb_value , segment_length , *frequency , *pisequence , *psegment , **sisequence;
  double sum , factorial_sum , diff , sum_square , segment_norm , sequence_norm ,
         likelihood , segmentation_likelihood , *sequence_mean , *residual , *prsequence ,
         *contrast , *norm , **factorial , **forward , *backward , *cumul_backward ,
         **srsequence , **mean , **variance;


  max_nb_value = 0;
  factorial = new double*[nb_variable];

  for (i = 1;i < nb_variable;i++) {
    if ((model_type[i - 1] == MULTINOMIAL_CHANGE) && (marginal[i]->nb_value > max_nb_value)) {
      max_nb_value = marginal[i]->nb_value;
    }

    if (model_type[i - 1] == POISSON_CHANGE) {
      factorial[i] = new double[length[index]];
    }
    else {
      factorial[i] = 0;
    }
  }

  if (max_nb_value > 0) {
    frequency = new int[max_nb_value];
  }
  else {
    frequency = 0;
  }

  sequence_mean = new double[nb_variable];
  residual = new double[length[index]];

  contrast = new double[length[index]];

  forward = new double*[length[index]];
  for (i = 0;i < length[index];i++) {
    forward[i] = new double[nb_segment];
  }

  norm = new double[length[index]];

  backward = new double[length[index]];
  cumul_backward = new double[length[index]];

  mean = new double*[nb_variable];
  variance = new double*[nb_variable];
  for (i = 1;i < nb_variable;i++) {
    if ((model_type[i - 1] == POISSON_CHANGE) || (model_type[i - 1] == GAUSSIAN_CHANGE) ||
        (model_type[i - 1] == VARIANCE_CHANGE) || (model_type[i - 1] == MEAN_VARIANCE_CHANGE)) {
      mean[i] = new double[nb_segment];
      variance[i] = new double[nb_segment];
    }
    else {
      mean[i] = 0;
      variance[i] = 0;
    }
  }

  sisequence = new int*[nb_variable];
  srsequence = new double*[nb_variable];

# ifdef DEBUG
  double **segmentation_probability;


  segmentation_probability = new double*[length[index]];
  for (i = 0;i < length[index];i++) {
    segmentation_probability[i] = new double[nb_segment];
    for (j = 0;j < nb_segment;j++) {
      segmentation_probability[i][j] = 0.;
    }
  }
# endif

  for (i = 1;i < nb_variable;i++) {
    if (model_type[i - 1] == VARIANCE_CHANGE) {
      sequence_mean[i] = 0.;

      if (type[i] != REAL_VALUE) {
        pisequence = int_sequence[index][i];
        for (j = 0;j < length[index];j++) {
          sequence_mean[i] += *pisequence++;
        }
      }

      else {
        prsequence = real_sequence[index][i];
        for (j = 0;j < length[index];j++) {
          sequence_mean[i] += *prsequence++;
        }
      }

      sequence_mean[i] /= length[index];
    }
  }

  // recurrence "forward"

  for (i = 0;i < length[index];i++) {

    // calcul des log-vraisemblances des segments

    for (j = 0;j <= i;j++) {
      contrast[j] = 0.;
    }

    for (j = 1;j < nb_variable;j++) {
      if (model_type[j - 1] == MULTINOMIAL_CHANGE) {
        for (k = 0;k < marginal[j]->nb_value;k++) {
          frequency[k] = 0;
        }
        sum = 0.;

        pisequence = int_sequence[index][j] + i;
        frequency[*pisequence--]++;
        for (k = i - 1;k >= 0;k--) {
          sum += (i - k) * log((double)(i - k) / (double)(i - k + 1)) +
                 log((double)(frequency[*pisequence] + 1) / (double)(i - k + 1));
          if (frequency[*pisequence] > 0) {
            sum -= frequency[*pisequence] *
                   log((double)frequency[*pisequence] / (double)(frequency[*pisequence] + 1));
          }
          frequency[*pisequence--]++;

          if (contrast[k] != D_INF) {
            contrast[k] += sum;
          }

/*          frequency[*pisequence--]++;
          if (contrast[k] != D_INF) {
            for (m = 0;m < marginal[j]->nb_value;m++) {
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

        pisequence = int_sequence[index][j] + i;
        for (k = i;k >= 0;k--) {
          sum += *pisequence--;
          factorial_sum += factorial[j][k];
          if ((contrast[k] != D_INF) && (sum > 0.)) {
            contrast[k] += sum * (log(sum / (i - k + 1)) - 1) - factorial_sum;
          }
        }
      }

      else {
        if (model_type[j - 1] == VARIANCE_CHANGE) {
          sum_square = 0.;

          if (type[j] != REAL_VALUE) {
            pisequence = int_sequence[index][j] + i;
            for (k = i;k >= 0;k--) {
              diff = *pisequence-- - sequence_mean[j];
              sum_square += diff * diff;
              residual[k] = sum_square;
            }
          }

          else {
            prsequence = real_sequence[index][j] + i;
            for (k = i;k >= 0;k--) {
              diff = *prsequence-- - sequence_mean[j];
              sum_square += diff * diff;
              residual[k] = sum_square;
            }
          }
        }

        else if (model_type[j - 1] == ORDINAL_GAUSSIAN_CHANGE) {
          pisequence = int_sequence[index][j] + i;
          sum_square = 0.;
          sum = rank[j][*pisequence--];
          residual[i] = 0.;

          for (k = i - 1;k >= 0;k--) {
            diff = rank[j][*pisequence] - sum / (i - k);
            sum_square += ((double)(i - k) / (double)(i - k + 1)) * diff * diff;
            sum += rank[j][*pisequence--];
            residual[k] = sum_square;
          }
        }

        else {
          if (type[j] != REAL_VALUE) {
            pisequence = int_sequence[index][j] + i;
            sum_square = 0.;
            sum = *pisequence--;
            residual[i] = 0.;

            for (k = i - 1;k >= 0;k--) {
              diff = *pisequence - sum / (i - k);
              sum_square += ((double)(i - k) / (double)(i - k + 1)) * diff * diff;
              sum += *pisequence--;
              residual[k] = sum_square;
            }
          }

          else {
            prsequence = real_sequence[index][j] + i;
            sum_square = 0.;
            sum = *prsequence--;
            residual[i] = 0.;

            for (k = i - 1;k >= 0;k--) {
              diff = *prsequence - sum / (i - k);
              sum_square += ((double)(i - k) / (double)(i - k + 1)) * diff * diff;
              sum += *prsequence--;
              residual[k] = sum_square;
            }
          }
        }

        if (model_type[j - 1] == MEAN_VARIANCE_CHANGE) {
          for (k = i - 1;k >= 0;k--) {
            contrast[k] += residual[k];
          }
        }

        else {
          for (k = i;k >= 0;k--) {
//            if ((contrast[k] != D_INF) && (residual[k] > 0.)) {
            if ((contrast[k] != D_INF) && (residual[k] > sqrt((double)(i - k + 1)) * ROUNDOFF_ERROR)) {
              contrast[k] -= ((double)(i - k + 1) / 2.) * (log(residual[k] /
                              (i - k + 1)) + log(2 * M_PI) + 1);
/*              contrast[k] -= ((double)(i - k + 1) / 2.) * (log(residual[k] /
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
//        if (contrast[j] > 0.) {
        if (contrast[j] > sqrt((double)((nb_variable - 1) * (i - j + 1))) * ROUNDOFF_ERROR) {
          contrast[j] = -((double)((nb_variable - 1) * (i - j + 1)) / 2.) * (log(contrast[j] /
                           ((nb_variable - 1) * (i - j + 1))) + log(2 * M_PI) + 1);
/*          contrast[j] = -((double)((nb_variable - 1) * (i - j + 1)) / 2.) * (log(contrast[j] /
                           ((nb_variable - 1) * (i - j))) + log(2 * M_PI)) +
                         (double)((nb_variable - 1) * (i - j)) / 2.; */
        }
        else {
          contrast[j] = D_INF;
        }
      }
    }

    if (contrast[i] != D_INF) {
      contrast[i] = exp(contrast[i]);
    }
    else {
      contrast[i] = 0.;
    }

    segment_norm = 0.;
    for (j = i - 1;j >= 0;j--) {
      segment_norm += norm[j];
      if (contrast[j] != D_INF) {
        contrast[j] = exp(contrast[j] - segment_norm);
      }
      else {
        contrast[j] = 0.;
      }
    }

    for (j = 0;j < nb_segment;j++) {
      forward[i][j] = 0.;
    }
    norm[i] = 0.;

    for (j = MAX(0 , nb_segment + i - length[index]);j < MIN((i < length[index] - 1 ? nb_segment - 1 : nb_segment) , i + 1);j++) {
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
      for (j = MAX(0 , nb_segment + i - length[index]);j < MIN((i < length[index] - 1 ? nb_segment - 1 : nb_segment) , i + 1);j++) {
        forward[i][j] /= norm[i];
      }

      norm[i] = log(norm[i]);
    }

#   ifdef DEBUG
    cout << i << " |";
    for (j = 0;j < nb_segment;j++) {
      cout << " " << forward[i][j];
    }
    cout << " | " << exp(norm[i]) << endl;
#   endif

  }

  sequence_norm = segment_norm + norm[length[index] - 1];

  if (forward[length[index] - 1][nb_segment - 1] > 0.) {
    likelihood = log(forward[length[index] - 1][nb_segment - 1]) + sequence_norm;
  }
  else {
    likelihood = D_INF;
  }

  if (likelihood != D_INF) {

    // passes "backward"

#   ifdef MESSAGE
    cout << "\n";
#   endif

    for (i = 0;i < nb_segmentation;i++) {
      j = length[index] - 1;
      psegment = int_sequence[index][0] + j;
      for (k = 1;k < nb_variable;k++) {
        if (type[k] != REAL_VALUE) {
          sisequence[k] = int_sequence[index][k] + j;
        }
        else {
          srsequence[k] = real_sequence[index][k] + j;
        }
      }
      segmentation_likelihood = sequence_norm;

      for (k = nb_segment - 1;k >= 0;k--) {

        // calcul des log-vraisemblances des segments

        for (m = k;m <= j;m++) {
          contrast[m] = 0.;
        }

        for (m = 1;m < nb_variable;m++) {
          if (model_type[m - 1] == MULTINOMIAL_CHANGE) {
            for (n = 0;n < marginal[m]->nb_value;n++) {
              frequency[n] = 0;
            }
            sum = 0.;

            pisequence = int_sequence[index][m] + j;
            frequency[*pisequence--]++;
            for (n = j - 1;n >= k;n--) {
              sum += (j - n) * log((double)(j - n) / (double)(j - n + 1)) +
                     log((double)(frequency[*pisequence] + 1) / (double)(j - n + 1));
              if (frequency[*pisequence] > 0) {
                sum -= frequency[*pisequence] *
                       log((double)frequency[*pisequence] / (double)(frequency[*pisequence] + 1));
              }
              frequency[*pisequence--]++;

              if (contrast[n] != D_INF) {
                contrast[n] += sum;
              }

/*              frequency[*pisequence--]++;
              if (contrast[n] != D_INF) {
                for (r = 0;r < marginal[m]->nb_value;r++) {
                  if (frequency[r] > 0) {
                    contrast[n] += frequency[r] * log((double)frequency[r] / (double)(j - n + 1));
                  }
                }
              } */
            }
          }

          else if (model_type[m - 1] == POISSON_CHANGE) {
            sum = 0.;
            factorial_sum = 0.;

            pisequence = int_sequence[index][m] + j;
            for (n = j;n >= k;n--) {
              sum += *pisequence--;
              factorial_sum += factorial[m][n];
              if ((contrast[n] != D_INF) && (sum > 0.)) {
                contrast[n] += sum * (log(sum / (j - n + 1)) - 1) - factorial_sum;
              }
            }
          }

          else {
            if (model_type[m - 1] == VARIANCE_CHANGE) {
              sum_square = 0.;

              if (type[m] != REAL_VALUE) {
                pisequence = int_sequence[index][m] + j;
                for (n = j;n >= k;n--) {
                  diff = *pisequence-- - sequence_mean[m];
                  sum_square += diff * diff;
                  residual[n] = sum_square;
                }
              }

              else {
                prsequence = real_sequence[index][m] + j;
                for (n = j;n >= k;n--) {
                  diff = *prsequence-- - sequence_mean[m];
                  sum_square += diff * diff;
                  residual[n] = sum_square;
                }
              }
            }

            else if (model_type[m - 1] == ORDINAL_GAUSSIAN_CHANGE) {
              pisequence = int_sequence[index][m] + j;
              sum_square = 0.;
              sum = rank[m][*pisequence--];
              residual[j] = 0.;

              for (n = j - 1;n >= k;n--) {
                diff = rank[m][*pisequence] - sum / (j - n);
                sum_square += ((double)(j - n) / (double)(j - n + 1)) * diff * diff;
                sum += rank[m][*pisequence--];
                residual[n] = sum_square;
              }
            }

            else {
              if (type[m] != REAL_VALUE) {
                pisequence = int_sequence[index][m] + j;
                sum_square = 0.;
                sum = *pisequence--;
                residual[j] = 0.;

                for (n = j - 1;n >= k;n--) {
                  diff = *pisequence - sum / (j - n);
                  sum_square += ((double)(j - n) / (double)(j - n + 1)) * diff * diff;
                  sum += *pisequence--;
                  residual[n] = sum_square;
                }
              }

              else {
                prsequence = real_sequence[index][m] + j;
                sum_square = 0.;
                sum = *prsequence--;
                residual[j] = 0.;

                for (n = j - 1;n >= k;n--) {
                  diff = *prsequence - sum / (j - n);
                  sum_square += ((double)(j - n) / (double)(j - n + 1)) * diff * diff;
                  sum += *prsequence--;
                  residual[n] = sum_square;
                }
              }
            }

            if (model_type[m - 1] == MEAN_VARIANCE_CHANGE) {
              for (n = j - 1;n >= k;n--) {
                contrast[n] += residual[n];
              }
            }

            else {
              for (n = j;n >= k;n--) {
//                if ((contrast[n] != D_INF) && (residual[n] > 0.)) {
                if ((contrast[n] != D_INF) && (residual[n] > sqrt((double)(j - n + 1)) * ROUNDOFF_ERROR)) {
                  contrast[n] -= ((double)(j - n + 1) / 2.) * (log(residual[n] /
                                   (j - n + 1)) + log(2 * M_PI) + 1);
/*                  contrast[n] -= ((double)(j - n + 1) / 2.) * (log(residual[n] /
                                   (j - n)) + log(2 * M_PI)) + (double)(j - n) / 2.; */
                }
                else {
                  contrast[n] = D_INF;
                }
              }
            }
          }
        }

        if (model_type[0] == MEAN_VARIANCE_CHANGE) {
          contrast[j] = D_INF;
          for (m = j - 1;m >= k;m--) {
//            if (contrast[m] > 0.)) {
            if (contrast[m] > sqrt((double)((nb_variable - 1) * (j - m + 1))) * ROUNDOFF_ERROR) {
              contrast[m] = -((double)((nb_variable - 1) * (j - m + 1)) / 2.) * (log(contrast[m] /
                               ((nb_variable - 1) * (j - m + 1))) + log(2 * M_PI) + 1);
/*              contrast[m] = -((double)((nb_variable - 1) * (j - m + 1)) / 2.) * (log(contrast[m] /
                               ((nb_variable - 1) * (j - m))) + log(2 * M_PI)) +
                             (double)((nb_variable - 1) * (j - m)) / 2.; */
            }
            else {
              contrast[m] = D_INF;
            }
          }
        }

        segment_norm = 0.;
        for (m = j;m >= k;m--) {
          segment_norm += norm[m];
          if (contrast[m] != D_INF) {
            contrast[m] = exp(contrast[m] - segment_norm);
          }
          else {
            contrast[m] = 0.;
          }
        }

        if (k > 0) {
          for (m = j;m >= k;m--) {
            backward[m] = contrast[m] * forward[m - 1][k - 1] / forward[j][k];
          }
          ::cumul_computation(j - k , backward + k , cumul_backward);
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

        segmentation_likelihood += log(contrast[j - segment_length + 1]);

        for (m = j;m > j - segment_length;m--) {
          *psegment-- = k;
        }

        for (m = 1;m < nb_variable;m++) {
          if ((model_type[m - 1] == POISSON_CHANGE) || (model_type[m - 1] == GAUSSIAN_CHANGE) ||
              (model_type[m - 1] == VARIANCE_CHANGE) || (model_type[m - 1] == MEAN_VARIANCE_CHANGE)) {
            mean[m][k] = 0.;
            variance[m][k] = 0.;

            if (type[m] != REAL_VALUE) {
              for (n = j;n > j - segment_length;n--) {
                mean[m][k] += *sisequence[m]--;
              }
              mean[m][k] /= segment_length;

              if (segment_length > 1) {
                sisequence[m] += segment_length;
                for (n = j;n > j - segment_length;n--) {
                  diff = *sisequence[m]-- - mean[m][k];
                  variance[m][k] += diff * diff;
                }
                variance[m][k] /= segment_length;
//                variance[m][k] /= (segment_length - 1);
              }
            }

            else {
              for (n = j;n > j - segment_length;n--) {
                mean[m][k] += *srsequence[m]--;
              }
              mean[m][k] /= segment_length;

              if (segment_length > 1) {
                srsequence[m] += segment_length;
                for (n = j;n > j - segment_length;n--) {
                  diff = *srsequence[m]-- - mean[m][k];
                  variance[m][k] += diff * diff;
                }
                variance[m][k] /= segment_length;
//                variance[m][k] /= (segment_length - 1);
              }
            }
          }
        }

        j -= segment_length;
      }

#     ifdef DEBUG
      psegment = int_sequence[index][0];
      for (j = 0;j < length[index];j++) {
        segmentation_probability[j][*psegment++]++;
      }
#     endif

#     ifdef MESSAGE
      if (i == 0) {
        os << "\n";
      }

      switch (format) {

      case 'a' : {
        psegment = int_sequence[index][0];
        for (j = 0;j < length[index];j++) {
          os << *psegment++ << " ";
        }

        os << "  #";
        os << "  " << i + 1 << "  " << segmentation_likelihood << "   ("
           << exp(segmentation_likelihood - likelihood) << ")" << endl;

        os << "# ";
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
          if (model_type[j - 1] == POISSON_CHANGE) {
            os << "# ";
            if (nb_variable > 2) {
              os << STAT_label[STATL_VARIABLE] << " " << j << "   ";
            }
            os << SEQ_label[SEQL_SEGMENT_MEAN] << "  " << SEQ_label[SEQL_SEGMENT_VARIANCE] << ": ";
            for (k = 0;k < nb_segment;k++) {
              os << mean[j][k] << " " << variance[j][k] << " | ";
            }
            os << endl;
          }

          else if ((model_type[j - 1] == GAUSSIAN_CHANGE) || (model_type[j - 1] == VARIANCE_CHANGE) ||
                   (model_type[j - 1] == MEAN_VARIANCE_CHANGE)) {
            os << "# ";
            if (nb_variable > 2) {
              os << STAT_label[STATL_VARIABLE] << " " << j << "   ";
            }
            os << SEQ_label[SEQL_SEGMENT_MEAN] << "  " << SEQ_label[SEQL_SEGMENT_STANDARD_DEVIATION] << ": ";
            for (k = 0;k < nb_segment;k++) {
              os << mean[j][k] << " " << sqrt(variance[j][k]) << " | ";
            }
            os << endl;
          }
        }
        break;
      }

      case 's' : {
        psegment = int_sequence[index][0];
        for (j = 0;j < length[index];j++) {
          os << *psegment++ << "\t";
        }

        os << "\t" << i + 1 << "\t" << segmentation_likelihood  << "\t"
           << exp(segmentation_likelihood - likelihood) << endl;

        for (j = 1;j < nb_variable;j++) {
          if (model_type[j - 1] == POISSON_CHANGE) {
            if (nb_variable > 2) {
              os << STAT_label[STATL_VARIABLE] << "\t" << j << "\t";
            }
            os << SEQ_label[SEQL_SEGMENT_MEAN] << "\t" << SEQ_label[SEQL_SEGMENT_VARIANCE] << "\t";
            for (k = 0;k < nb_segment;k++) {
              os << mean[j][k] << "\t" << variance[j][k] << "\t";
            }
            os << endl;
          }

          else if ((model_type[j - 1] == GAUSSIAN_CHANGE) || (model_type[j - 1] == VARIANCE_CHANGE) ||
                   (model_type[j - 1] == MEAN_VARIANCE_CHANGE)) {
            if (nb_variable > 2) {
              os << STAT_label[STATL_VARIABLE] << "\t" << j << "\t";
            }
            os << SEQ_label[SEQL_SEGMENT_MEAN] << "\t" << SEQ_label[SEQL_SEGMENT_STANDARD_DEVIATION] << "\t";
            for (k = 0;k < nb_segment;k++) {
              os << mean[j][k] << "\t" << sqrt(variance[j][k]) << "\t";
            }
            os << endl;
          }
        }
        break;
      }
      }
#     endif

    }

#   ifdef DEBUG
    if (nb_segmentation >= 1000) {
      for (i = 0;i < length[index];i++) {
        for (j = 0;j < nb_segment;j++) {
          segmentation_probability[i][j] /= nb_segmentation;
        }
      }

      os << "\n" << SEQ_label[SEQL_POSTERIOR_SEGMENT_PROBABILITY] << "\n\n";

      psegment = int_sequence[index][0];
      for (j = 0;j < length[index];j++) {
        *psegment++ = I_DEFAULT;
      }

      profile_ascii_print(os , index , nb_segment , segmentation_probability ,
                          SEQ_label[SEQL_SEGMENT]);
    }
#   endif

  }

  delete [] frequency;

  for (i = 1;i < nb_variable;i++) {
    delete [] factorial[i];
  }
  delete [] factorial;

  delete [] sequence_mean;
  delete [] residual;

  delete [] contrast;

  for (i = 0;i < length[index];i++) {
    delete [] forward[i];
  }
  delete [] forward;

  delete [] norm;

  delete [] backward;
  delete [] cumul_backward;

  for (i = 1;i < nb_variable;i++) {
    delete [] mean[i];
    delete [] variance[i];
  }
  delete [] mean;
  delete [] variance;

  delete [] sisequence;
  delete [] srsequence;

# ifdef DEBUG
  for (i = 0;i < length[index];i++) {
    delete [] segmentation_probability[i];
  }
  delete [] segmentation_probability;
# endif

  return likelihood;
}


/*--------------------------------------------------------------*
 *
 *  Calcul des L segmentations optimales d'une sequence.
 *
 *  arguments : indice de la sequence, nombre de segments, types des modeles,
 *              rangs (variables ordinales), stream, format de fichier ('a' : ASCII,
 *              's' : Spreadsheet), nombre de segmentation, vraisemblance des donnees.
 *
 *--------------------------------------------------------------*/

double Sequences::L_segmentation(int index , int nb_segment , int *model_type ,
                                 double **irank , ostream &os , char format ,
                                 int inb_segmentation , double likelihood) const

{
  bool **active_cell;
  register int i , j , k , m , n;
  int max_nb_value , brank , previous_rank , nb_cell , *frequency , *pisequence , *rank ,
      *psegment , **sisequence , ***optimal_length , ***optimal_rank;
  double sum , factorial_sum , diff , sum_square , buff , segmentation_likelihood ,
         *nb_segmentation , *sequence_mean , *residual , *prsequence , *contrast ,
         **factorial , **srsequence , **mean , **variance , ***forward;
  long double likelihood_cumul;


# ifdef MESSAGE
  double max_nb_segmentation , sum2;

  max_nb_segmentation = 1.;
  for (i = 1;i < nb_segment;i++) {
    max_nb_segmentation *= (double)(length[index] - i) / (double)i;
//    max_nb_segmentation = max_nb_segmentation * (length[index] - i) / i;
  }
  os << "\n" << SEQ_label[SEQL_NB_SEGMENTATION] << ": " << max_nb_segmentation << endl;
# endif

  // initialisations

  max_nb_value = 0;
  factorial = new double*[nb_variable];

  for (i = 1;i < nb_variable;i++) {
    if ((model_type[i - 1] == MULTINOMIAL_CHANGE) && (marginal[i]->nb_value > max_nb_value)) {
      max_nb_value = marginal[i]->nb_value;
    }

    if (model_type[i - 1] == POISSON_CHANGE) {
      factorial[i] = new double[length[index]];
    }
    else {
      factorial[i] = 0;
    }
  }

  if (max_nb_value > 0) {
    frequency = new int[max_nb_value];
  }
  else {
    frequency = 0;
  }

  sequence_mean = new double[nb_variable];
  residual = new double[length[index]];

  contrast = new double[length[index]];

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
  for (i = 1;i < nb_variable;i++) {
    if ((model_type[i - 1] == POISSON_CHANGE) || (model_type[i - 1] == GAUSSIAN_CHANGE) ||
        (model_type[i - 1] == MEAN_CHANGE) || (model_type[i - 1] == VARIANCE_CHANGE) ||
        (model_type[i - 1] == MEAN_VARIANCE_CHANGE)) {
      mean[i] = new double[nb_segment];
      variance[i] = new double[nb_segment];
    }
    else {
      mean[i] = 0;
      variance[i] = 0;
    }
  }

  active_cell = new bool*[length[index]];
  for (i = 0;i < length[index];i++) {
    active_cell[i] = new bool[nb_segment];
    for (j = 0;j < nb_segment;j++) {
      active_cell[i][j] = false;
    }
  }

  sisequence = new int*[nb_variable];
  srsequence = new double*[nb_variable];

# ifdef MESSAGE
  double **mean_square_diff;


  mean_square_diff = new double*[nb_variable];
  for (i = 1;i < nb_variable;i++) {
    if ((model_type[i - 1] == GAUSSIAN_CHANGE) || (model_type[i - 1] == MEAN_CHANGE) ||
        (model_type[i - 1] == MEAN_VARIANCE_CHANGE)) {
      mean_square_diff[i] = new double[length[index]];
    }
    else {
      mean_square_diff[i] = 0;
    }
  }
# endif

# ifdef DEBUG
  double **segmentation_probability;


  segmentation_probability = new double*[length[index]];
  for (i = 0;i < length[index];i++) {
    segmentation_probability[i] = new double[nb_segment];
    for (j = 0;j < nb_segment;j++) {
//      segmentation_probability[i][j] = 0.;
      segmentation_probability[i][j] = D_INF;
    }
  }
# endif

  for (i = 1;i < nb_variable;i++) {
    if (model_type[i - 1] == VARIANCE_CHANGE) {
      sequence_mean[i] = 0.;

      if (type[i] != REAL_VALUE) {
        pisequence = int_sequence[index][i];
        for (j = 0;j < length[index];j++) {
          sequence_mean[i] += *pisequence++;
        }
      }

      else {
        prsequence = real_sequence[index][i];
        for (j = 0;j < length[index];j++) {
          sequence_mean[i] += *prsequence++;
        }
      }

      sequence_mean[i] /= length[index];
    }
  }

/*  for (i = 0;i < nb_segment;i++) {
    nb_segmentation[i] = 1;
  } */

  // recurrence "forward"

  for (i = 0;i < length[index];i++) {

    // calcul des log-vraisemblances des segments

    for (j = 0;j <= i;j++) {
      contrast[j] = 0.;
    }

    for (j = 1;j < nb_variable;j++) {
      if (model_type[j - 1] == MULTINOMIAL_CHANGE) {
        for (k = 0;k < marginal[j]->nb_value;k++) {
          frequency[k] = 0;
        }
        sum = 0.;

        pisequence = int_sequence[index][j] + i;
        frequency[*pisequence--]++;
        for (k = i - 1;k >= 0;k--) {
          sum += (i - k) * log((double)(i - k) / (double)(i - k + 1)) +
                 log((double)(frequency[*pisequence] + 1) / (double)(i - k + 1));
          if (frequency[*pisequence] > 0) {
            sum -= frequency[*pisequence] *
                   log((double)frequency[*pisequence] / (double)(frequency[*pisequence] + 1));
          }
          frequency[*pisequence--]++;

          if (contrast[k] != D_INF) {
            contrast[k] += sum;
          }

#         ifdef MESSAGE
          sum2 = 0.;
          for (m = 0;m < marginal[j]->nb_value;m++) {
            if (frequency[m] > 0) {
              sum2 += frequency[m] * log((double)frequency[m] / (double)(i - k + 1));
            }
          }

          if ((sum2 < sum - DOUBLE_ERROR) || (sum2 > sum + DOUBLE_ERROR)) {
            cout << "\nERROR: " << i << " " << k << " | " << sum2 << " " << sum << endl;
          }
#         endif

        }
      }

      else if (model_type[j - 1] == POISSON_CHANGE) {
        factorial[j][i] = 0.;
        for (k = 2;k <= int_sequence[index][j][i];k++) {
          factorial[j][i] += log((double)k);
        }

        sum = 0.;
        factorial_sum = 0.;

        pisequence = int_sequence[index][j] + i;
        for (k = i;k >= 0;k--) {
          sum += *pisequence--;
          factorial_sum += factorial[j][k];
          if ((contrast[k] != D_INF) && (sum > 0.)) {
            contrast[k] += sum * (log(sum / (i - k + 1)) - 1) - factorial_sum;
          }
        }
      }

      else {
        if (model_type[j - 1] == VARIANCE_CHANGE) {
          sum_square = 0.;

          if (type[j] != REAL_VALUE) {
            pisequence = int_sequence[index][j] + i;
            for (k = i;k >= 0;k--) {
              diff = *pisequence-- - sequence_mean[j];
              sum_square += diff * diff;
              residual[k] = sum_square;
            }
          }

          else {
            prsequence = real_sequence[index][j] + i;
            for (k = i;k >= 0;k--) {
              diff = *prsequence-- - sequence_mean[j];
              sum_square += diff * diff;
              residual[k] = sum_square;
            }
          }
        }

        else if (model_type[j - 1] == ORDINAL_GAUSSIAN_CHANGE) {
          pisequence = int_sequence[index][j] + i;
          sum_square = 0.;
          sum = irank[j][*pisequence--];
          residual[i] = 0.;

          for (k = i - 1;k >= 0;k--) {
            diff = irank[j][*pisequence] - sum / (i - k);
            sum_square += ((double)(i - k) / (double)(i - k + 1)) * diff * diff;
            sum += irank[j][*pisequence--];
            residual[k] = sum_square;
          }
        }

        else {
          if (type[j] != REAL_VALUE) {
            pisequence = int_sequence[index][j] + i;
            sum_square = 0.;
            sum = *pisequence--;
            residual[i] = 0.;

            for (k = i - 1;k >= 0;k--) {
              diff = *pisequence - sum / (i - k);
              sum_square += ((double)(i - k) / (double)(i - k + 1)) * diff * diff;
              sum += *pisequence--;
              residual[k] = sum_square;
            }

#           ifdef MESSAGE
//            cout << "\n";

            pisequence = int_sequence[index][j] + i;
            sum_square = *pisequence * *pisequence;
            sum = *pisequence--;

            for (k = i - 1;k >= 0;k--) {
              sum_square += *pisequence * *pisequence;
              sum += *pisequence--;

              if ((sum_square - sum * sum / (i - k + 1) < residual[k] - DOUBLE_ERROR) ||
                  (sum_square - sum * sum / (i - k + 1) > residual[k] + DOUBLE_ERROR)) {
                cout << "\nERROR: " << i << " " << k << " | " << sum_square - sum * sum / (i - k + 1)
                     << " " << residual[k] << endl;
              }
            }

/*            mean_square_diff[j][i] = 0.;
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
#           endif

          }

          else {
            prsequence = real_sequence[index][j] + i;
            sum_square = 0.;
            sum = *prsequence--;
            residual[i] = 0.;

            for (k = i - 1;k >= 0;k--) {
              diff = *prsequence - sum / (i - k);
              sum_square += ((double)(i - k) / (double)(i - k + 1)) * diff * diff;
              sum += *prsequence--;
              residual[k] = sum_square;
            }
          }
        }

        if (model_type[j - 1] == MEAN_CHANGE) {
          for (k = i - 1;k >= 0;k--) {
            contrast[k] -= residual[k];
          }
        }

        else if (model_type[j - 1] == MEAN_VARIANCE_CHANGE) {
          for (k = i - 1;k >= 0;k--) {
            contrast[k] += residual[k];
          }
        }

        else {
          for (k = i;k >= 0;k--) {
//            if ((contrast[k] != D_INF) && (residual[k] > 0.)) {
            if ((contrast[k] != D_INF) && (residual[k] > sqrt((double)(i - k + 1)) * ROUNDOFF_ERROR)) {
              contrast[k] -= ((double)(i - k + 1) / 2.) * (log(residual[k] /
                               (i - k + 1)) + log(2 * M_PI) + 1);
/*              contrast[k] -= ((double)(i - k + 1) / 2.) * (log(residual[k] /
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
//        if (contrast[j] > 0.) {
        if (contrast[j] > sqrt((double)((nb_variable - 1) * (i - j + 1))) * ROUNDOFF_ERROR) {
          contrast[j] = -((double)((nb_variable - 1) * (i - j + 1)) / 2.) * (log(contrast[j] /
                           ((nb_variable - 1) * (i - j + 1))) + log(2 * M_PI) + 1);
/*          contrast[j] = -((double)((nb_variable - 1) * (i - j + 1)) / 2.) * (log(contrast[j] /
                           ((nb_variable - 1) * (i - j))) + log(2 * M_PI)) +
                         (double)((nb_variable - 1) * (i - j)) / 2.; */
        }
        else {
          contrast[j] = D_INF;
        }
      }
    }

#   ifdef DEBUG
    for (j = i - 1;j >= 0;j--) {
      cout << contrast[j] << "  ";
    }
    cout << endl;
#   endif

    nb_segmentation[0] = 1;
    for (j = 1;j < MIN((i < length[index] - 1 ? nb_segment - 1 : nb_segment) , i + 1);j++) {
      nb_segmentation[j] = nb_segmentation[j - 1] * (i - j + 1) / j;
    }
    for (j = 1;j < MIN((i < length[index] - 1 ? nb_segment - 1 : nb_segment) , i + 1);j++) {
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

  // restauration

  likelihood_cumul = 0.;

  for (i = 0;i < nb_segmentation[nb_segment - 1];i++) {
    if (forward[length[index] - 1][nb_segment - 1][i] == D_INF) {
      break;
    }

#   ifdef DEBUG
    cout << "\n";
#   endif

    j = length[index] - 1;
    psegment = int_sequence[index][0] + j;
    for (k = 1;k < nb_variable;k++) {
      if (type[k] != REAL_VALUE) {
        sisequence[k] = int_sequence[index][k] + j;
      }
      else {
        srsequence[k] = real_sequence[index][k] + j;
      }
    }
    brank = i;

    for (k = nb_segment - 1;k >= 0;k--) {
      for (m = j;m > j - optimal_length[j][k][brank];m--) {
        active_cell[m][k] = true;
        *psegment-- = k;
      }

      for (m = 1;m < nb_variable;m++) {
        if ((model_type[m - 1] == POISSON_CHANGE) || (model_type[m - 1] == GAUSSIAN_CHANGE) ||
            (model_type[m - 1] == MEAN_CHANGE) || (model_type[m - 1] == VARIANCE_CHANGE) ||
            (model_type[m - 1] == MEAN_VARIANCE_CHANGE)) {
          mean[m][k] = 0.;
          variance[m][k] = 0.;

          if (type[m] != REAL_VALUE) {
            for (n = j;n > j - optimal_length[j][k][brank];n--) {
              mean[m][k] += *sisequence[m]--;
            }
            mean[m][k] /= optimal_length[j][k][brank];

            if (optimal_length[j][k][brank] > 1) {
              sisequence[m] += optimal_length[j][k][brank];
              for (n = j;n > j - optimal_length[j][k][brank];n--) {
                diff = *sisequence[m]-- - mean[m][k];
                variance[m][k] += diff * diff;
              }
              variance[m][k] /= optimal_length[j][k][brank];
//              variance[m][k] /= (optimal_length[j][k][brank] - 1);
            }
          }

          else {
            for (n = j;n > j - optimal_length[j][k][brank];n--) {
              mean[m][k] += *srsequence[m]--;
            }
            mean[m][k] /= optimal_length[j][k][brank];

            if (optimal_length[j][k][brank] > 1) {
              srsequence[m] += optimal_length[j][k][brank];
              for (n = j;n > j - optimal_length[j][k][brank];n--) {
                diff = *srsequence[m]-- - mean[m][k];
                variance[m][k] += diff * diff;
              }
              variance[m][k] /= optimal_length[j][k][brank];
//              variance[m][k] /= (optimal_length[j][k][brank] - 1);
            }
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

    if (model_type[0] == MEAN_CHANGE) {
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
      likelihood_cumul += expl(forward[length[index] - 1][nb_segment - 1][i]);
    }

#   ifdef DEBUG
    psegment = int_sequence[index][0];
    for (j = 0;j < length[index];j++) {
/*      segmentation_probability[j][*psegment++] += exp(forward[length[index] - 1][nb_segment - 1][i] - likelihood);

      if (((i == 0) || (*psegment != *(psegment - 1))) &&
          (forward[length[index] - 1][nb_segment - 1][i] > segmentation_probability[j][*psegment])) { */
      if (forward[length[index] - 1][nb_segment - 1][i] > segmentation_probability[j][*psegment]) {
        segmentation_probability[j][*psegment] = forward[length[index] - 1][nb_segment - 1][i];
      }
      psegment++;
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
      psegment = int_sequence[index][0];
      for (j = 0;j < length[index];j++) {
        os << *psegment++ << " ";
      }

      os << "  " << i + 1 << "  " << forward[length[index] - 1][nb_segment - 1][i] << "   (";
      if (likelihood != D_INF) {
        os << exp(forward[length[index] - 1][nb_segment - 1][i] - likelihood)
           << "  " << likelihood_cumul / expl(likelihood);
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
        if (model_type[j - 1] == POISSON_CHANGE) {
          if (nb_variable > 2) {
            os << STAT_label[STATL_VARIABLE] << " " << j << "   ";
          }
          os << SEQ_label[SEQL_SEGMENT_MEAN] << "  " << SEQ_label[SEQL_SEGMENT_VARIANCE] << ": ";
          for (k = 0;k < nb_segment;k++) {
            os << mean[j][k] << " " << variance[j][k] << " | ";
          }
          os << endl;
        }

        else if ((model_type[j - 1] == GAUSSIAN_CHANGE) || (model_type[j - 1] == MEAN_CHANGE) ||
                 (model_type[j - 1] == VARIANCE_CHANGE) || (model_type[j - 1] == MEAN_VARIANCE_CHANGE)) {
          if (nb_variable > 2) {
            os << STAT_label[STATL_VARIABLE] << " " << j << "   ";
          }
          os << SEQ_label[SEQL_SEGMENT_MEAN] << "  " << SEQ_label[SEQL_SEGMENT_STANDARD_DEVIATION] << ": ";
          for (k = 0;k < nb_segment;k++) {
            os << mean[j][k] << " " << sqrt(variance[j][k]) << " | ";
          }
          os << endl;
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
        os << exp(forward[length[index] - 1][nb_segment - 1][i] - likelihood)
           << "\t" << likelihood_cumul / expl(likelihood);
      }
      else {
        os << exp(forward[length[index] - 1][nb_segment - 1][i] - segmentation_likelihood);
      }
      os << "\t" << nb_cell << endl;

      for (j = 1;j < nb_variable;j++) {
        if (model_type[j - 1] == POISSON_CHANGE) {
          if (nb_variable > 2) {
            os << STAT_label[STATL_VARIABLE] << "\t" << j << "\t";
          }
          os << SEQ_label[SEQL_SEGMENT_MEAN] << "\t" << SEQ_label[SEQL_SEGMENT_VARIANCE] << "\t";
          for (k = 0;k < nb_segment;k++) {
            os << mean[j][k] << "\t" << variance[j][k] << "\t";
          }
          os << endl;
        }

        else if ((model_type[j - 1] == GAUSSIAN_CHANGE) || (model_type[j - 1] == MEAN_CHANGE) ||
                 (model_type[j - 1] == VARIANCE_CHANGE) || (model_type[j - 1] == MEAN_VARIANCE_CHANGE)) {
          if (nb_variable > 2) {
            os << STAT_label[STATL_VARIABLE] << "\t" << j << "\t";
          }
          os << SEQ_label[SEQL_SEGMENT_MEAN] << "\t" << SEQ_label[SEQL_SEGMENT_STANDARD_DEVIATION] << "\t";
          for (k = 0;k < nb_segment;k++) {
            os << mean[j][k] << "\t" << sqrt(variance[j][k]) << "\t";
          }
          os << endl;
        }
      }
      break;
    }
    }
#   endif

  }

# ifdef DEBUG
  if (((likelihood != D_INF) && (likelihood_cumul / expl(likelihood) > 0.8)) ||
      (segmentation_likelihood != D_INF)) {
    if (likelihood != D_INF) {
      for (i = 0;i < length[index];i++) {
        for (j = 0;j < nb_segment;j++) {
          if (segmentation_probability[i][j] != D_INF) {
            segmentation_probability[i][j] = exp(segmentation_probability[i][j] - likelihood);
          }
          else {
            segmentation_probability[i][j] = 0.;
          }
        }
      }

//      os << "\n" << SEQ_label[SEQL_POSTERIOR_SEGMENT_PROBABILITY] << "\n\n";
      os << "\n" << SEQ_label[SEQL_MAX_POSTERIOR_SEGMENT_PROBABILITY] << "\n\n";
    }

    else {
      for (i = 0;i < length[index];i++) {
        for (j = 0;j < nb_segment;j++) {
          if (segmentation_probability[i][j] != D_INF) {
            segmentation_probability[i][j] = exp(segmentation_probability[i][j] - segmentation_likelihood);
          }
          else {
            segmentation_probability[i][j] = 0.;
          }
        }
      }

      os << "\n" << SEQ_label[SEQL_MAX_SEGMENT_LIKELIHOOD] << "\n\n";
    }

    psegment = int_sequence[index][0];
    for (j = 0;j < length[index];j++) {
      *psegment++ = I_DEFAULT;
    }

    profile_ascii_print(os , index , nb_segment , segmentation_probability ,
                        SEQ_label[SEQL_SEGMENT]);
  }
# endif

  delete [] frequency;

  for (i = 1;i < nb_variable;i++) {
    delete [] factorial[i];
  }
  delete [] factorial;

  delete [] sequence_mean;
  delete [] residual;

  delete [] contrast;

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
  }
  delete [] mean;
  delete [] variance;

  for (i = 0;i < length[index];i++) {
    delete [] active_cell[i];
  }
  delete [] active_cell;

  delete [] sisequence;
  delete [] srsequence;

# ifdef MESSAGE
  for (i = 1;i < nb_variable;i++) {
    delete [] mean_square_diff[i];
  }
  delete [] mean_square_diff;
# endif

# ifdef DEBUG
  for (i = 0;i < length[index];i++) {
    delete [] segmentation_probability[i];
  }
  delete [] segmentation_probability;
# endif

  return segmentation_likelihood;
}


/*--------------------------------------------------------------*
 *
 *  Calcul par maximisation des profils de segmentation d'une sequence.
 *
 *  arguments : indice de la sequence, nombre de segments, types des modeles,
 *              rangs (variables ordinales), stream, type de sortie, format
 *              de fichier ('a' : ASCII, 's' : Spreadsheet, 'g' : Gnuplot),
 *              vraisemblance des donnees.
 *
 *--------------------------------------------------------------*/

double Sequences::forward_backward_dynamic_programming(int index , int nb_segment , int *model_type ,
                                                       double **rank , ostream &os , int output ,
                                                       char format , double likelihood) const

{
  register int i , j , k , m;
  int max_nb_value , *frequency , *pisequence , *psegment , **optimal_length;
  double sum , factorial_sum , diff , sum_square , buff , segmentation_likelihood ,
         backward_max , *sequence_mean , *residual , *prsequence , *contrast ,
         **factorial , **forward , **backward , **backward_output , **mean;


  max_nb_value = 0;
  factorial = new double*[nb_variable];

  for (i = 1;i < nb_variable;i++) {
    if ((model_type[i - 1] == MULTINOMIAL_CHANGE) && (marginal[i]->nb_value > max_nb_value)) {
      max_nb_value = marginal[i]->nb_value;
    }

    if (model_type[i - 1] == POISSON_CHANGE) {
      factorial[i] = new double[length[index]];
    }
    else {
      factorial[i] = 0;
    }
  }

  if (max_nb_value > 0) {
    frequency = new int[max_nb_value];
  }
  else {
    frequency = 0;
  }

  sequence_mean = new double[nb_variable];
  residual = new double[length[index]];

  contrast = new double[length[index]];

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

  mean = new double*[nb_variable];
  for (i = 1;i < nb_variable;i++) {
    if ((model_type[i - 1] == POISSON_CHANGE) || (model_type[i - 1] == GAUSSIAN_CHANGE) ||
        (model_type[i - 1] == MEAN_CHANGE) || (model_type[i - 1] == VARIANCE_CHANGE) ||
        (model_type[i - 1] == MEAN_VARIANCE_CHANGE)) {
      mean[i] = new double[length[index]];
    }
    else {
      mean[i] = 0;
    }
  }

  for (i = 1;i < nb_variable;i++) {
    if (model_type[i - 1] == VARIANCE_CHANGE) {
      sequence_mean[i] = 0.;

      if (type[i] != REAL_VALUE) {
        pisequence = int_sequence[index][i];
        for (j = 0;j < length[index];j++) {
          sequence_mean[i] += *pisequence++;
        }
      }

      else {
        prsequence = real_sequence[index][i];
        for (j = 0;j < length[index];j++) {
          sequence_mean[i] += *prsequence++;
        }
      }

      sequence_mean[i] /= length[index];
    }
  }

  // recurrence "forward"

  for (i = 0;i < length[index];i++) {

    // calcul des log-vraisemblances des segments

    for (j = 0;j <= i;j++) {
      contrast[j] = 0.;
    }

    for (j = 1;j < nb_variable;j++) {
      if (model_type[j - 1] == MULTINOMIAL_CHANGE) {
        for (k = 0;k < marginal[j]->nb_value;k++) {
          frequency[k] = 0;
        }
        sum = 0.;

        pisequence = int_sequence[index][j] + i;
        frequency[*pisequence--]++;
        for (k = i - 1;k >= 0;k--) {
          sum += (i - k) * log((double)(i - k) / (double)(i - k + 1)) +
                 log((double)(frequency[*pisequence] + 1) / (double)(i - k + 1));
          if (frequency[*pisequence] > 0) {
            sum -= frequency[*pisequence] *
                   log((double)frequency[*pisequence] / (double)(frequency[*pisequence] + 1));
          }
          frequency[*pisequence--]++;

          if (contrast[k] != D_INF) {
            contrast[k] += sum;
          }

/*          frequency[*pisequence--]++;
          if (contrast[k] != D_INF) {
            for (m = 0;m < marginal[j]->nb_value;m++) {
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

        pisequence = int_sequence[index][j] + i;
        for (k = i;k >= 0;k--) {
          sum += *pisequence--;
          factorial_sum += factorial[j][k];
          if ((contrast[k] != D_INF) && (sum > 0.)) {
            contrast[k] += sum * (log(sum / (i - k + 1)) - 1) - factorial_sum;
          }
        }
      }

      else {
        if (model_type[j - 1] == VARIANCE_CHANGE) {
          sum_square = 0.;

          if (type[j] != REAL_VALUE) {
            pisequence = int_sequence[index][j] + i;
            for (k = i;k >= 0;k--) {
              diff = *pisequence-- - sequence_mean[j];
              sum_square += diff * diff;
              residual[k] = sum_square;
            }
          }

          else {
            prsequence = real_sequence[index][j] + i;
            for (k = i;k >= 0;k--) {
              diff = *prsequence-- - sequence_mean[j];
              sum_square += diff * diff;
              residual[k] = sum_square;
            }
          }
        }

        else if (model_type[j - 1] == ORDINAL_GAUSSIAN_CHANGE) {
          pisequence = int_sequence[index][j] + i;
          sum_square = 0.;
          sum = rank[j][*pisequence--];
          residual[i] = 0.;

          for (k = i - 1;k >= 0;k--) {
            diff = rank[j][*pisequence] - sum / (i - k);
            sum_square += ((double)(i - k) / (double)(i - k + 1)) * diff * diff;
            sum += rank[j][*pisequence--];
            residual[k] = sum_square;
          }
        }

        else {
          if (type[j] != REAL_VALUE) {
            pisequence = int_sequence[index][j] + i;
            sum_square = 0.;
            sum = *pisequence--;
            residual[i] = 0.;

            for (k = i - 1;k >= 0;k--) {
              diff = *pisequence - sum / (i - k);
              sum_square += ((double)(i - k) / (double)(i - k + 1)) * diff * diff;
              sum += *pisequence--;
              residual[k] = sum_square;
            }
          }

          else {
            prsequence = real_sequence[index][j] + i;
            sum_square = 0.;
            sum = *prsequence--;
            residual[i] = 0.;

            for (k = i - 1;k >= 0;k--) {
              diff = *prsequence - sum / (i - k);
              sum_square += ((double)(i - k) / (double)(i - k + 1)) * diff * diff;
              sum += *prsequence--;
              residual[k] = sum_square;
            }
          }
        }

        if (model_type[j - 1] == MEAN_CHANGE) {
          for (k = i - 1;k >= 0;k--) {
            contrast[k] -= residual[k];
          }
        }

        else if (model_type[j - 1] == MEAN_VARIANCE_CHANGE) {
          for (k = i - 1;k >= 0;k--) {
            contrast[k] += residual[k];
          }
        }

        else {
          for (k = i;k >= 0;k--) {
//            if ((contrast[k] != D_INF) && (residual[k] > 0.)) {
            if ((contrast[k] != D_INF) && (residual[k] > sqrt((double)(i - k + 1)) * ROUNDOFF_ERROR)) {
              contrast[k] -= ((double)(i - k + 1) / 2.) * (log(residual[k] /
                               (i - k + 1)) + log(2 * M_PI) + 1);
/*              contrast[k] -= ((double)(i - k + 1) / 2.) * (log(residual[k] /
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
//        if (contrast[j] > 0.) {
        if (contrast[j] > sqrt((double)((nb_variable - 1) * (i - j + 1))) * ROUNDOFF_ERROR) {
          contrast[j] = -((double)((nb_variable - 1) * (i - j + 1)) / 2.) * (log(contrast[j] /
                           ((nb_variable - 1) * (i - j + 1))) + log(2 * M_PI) + 1);
/*          contrast[j] = -((double)((nb_variable - 1) * (i - j + 1)) / 2.) * (log(contrast[j] /
                           ((nb_variable - 1) * (i - j))) + log(2 * M_PI)) +
                         (double)((nb_variable - 1) * (i - j)) / 2.; */
        }
        else {
          contrast[j] = D_INF;
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

      for (j = i;j < length[index];j++) {
        contrast[j] = 0.;
      }

      for (j = 1;j < nb_variable;j++) {
        if (model_type[j - 1] == MULTINOMIAL_CHANGE) {
          for (k = 0;k < marginal[j]->nb_value;k++) {
            frequency[k] = 0;
          }
          sum = 0.;

          pisequence = int_sequence[index][j] + i;
          frequency[*pisequence++]++;
          for (k = i + 1;k < length[index];k++) {
            sum += (k - i) * log((double)(k - i) / (double)(k - i + 1)) +
                   log((double)(frequency[*pisequence] + 1) / (double)(k - i + 1));
            if (frequency[*pisequence] > 0) {
              sum -= frequency[*pisequence] *
                     log((double)frequency[*pisequence] / (double)(frequency[*pisequence] + 1));
            }
            frequency[*pisequence++]++;

            if (contrast[k] != D_INF) {
              contrast[k] += sum;
            }

/*            frequency[*pisequence++]++;
            if (contrast[k] != D_INF) {
              for (m = 0;m < marginal[j]->nb_value;m++) {
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

          pisequence = int_sequence[index][j] + i;
          for (k = i;k < length[index];k++) {
            sum += *pisequence++;
            factorial_sum += factorial[j][k];
            if ((contrast[k] != D_INF) && (sum > 0.)) {
              contrast[k] += sum * (log(sum / (k - i + 1)) - 1) - factorial_sum;
            }
          }
        }

        else {
          if (model_type[j - 1] == VARIANCE_CHANGE) {
            sum_square = 0.;

            if (type[j] != REAL_VALUE) {
              pisequence = int_sequence[index][j] + i;
              for (k = i;k < length[index];k++) {
                diff = *pisequence++ - sequence_mean[j];
                sum_square += diff * diff;
                residual[k] = sum_square;
              }
            }

            else {
              prsequence = real_sequence[index][j] + i;
              for (k = i;k < length[index];k++) {
                diff = *prsequence++ - sequence_mean[j];
                sum_square += diff * diff;
                residual[k] = sum_square;
              }
            }
          }

          else if (model_type[j - 1] == ORDINAL_GAUSSIAN_CHANGE) {
            pisequence = int_sequence[index][j] + i;
            sum_square = 0.;
            sum = rank[j][*pisequence++];
            residual[i] = 0.;

            for (k = i + 1;k < length[index];k++) {
              diff = rank[j][*pisequence] - sum / (k - i);
              sum_square += ((double)(k - i) / (double)(k - i + 1)) * diff * diff;
              sum += rank[j][*pisequence++];
              residual[k] = sum_square;
            }
          }

          else {
            if (type[j] != REAL_VALUE) {
              pisequence = int_sequence[index][j] + i;
              sum_square = 0.;
              sum = *pisequence++;
              residual[i] = 0.;

              for (k = i + 1;k < length[index];k++) {
                diff = *pisequence - sum / (k - i);
                sum_square += ((double)(k - i) / (double)(k - i + 1)) * diff * diff;
                sum += *pisequence++;
                residual[k] = sum_square;
              }
            }

            else {
              prsequence = real_sequence[index][j] + i;
              sum_square = 0.;
              sum = *prsequence++;
              residual[i] = 0.;

              for (k = i + 1;k < length[index];k++) {
                diff = *prsequence - sum / (k - i);
                sum_square += ((double)(k - i) / (double)(k - i + 1)) * diff * diff;
                sum += *prsequence++;
                residual[k] = sum_square;
              }
            }
          }

          if (model_type[j - 1] == MEAN_CHANGE) {
            for (k = i + 1;k < length[index];k++) {
              contrast[k] -= residual[k];
            }
          }

          else if (model_type[j - 1] == MEAN_VARIANCE_CHANGE) {
            for (k = i + 1;k < length[index];k++) {
              contrast[k] += residual[k];
            }
          }

          else {
            for (k = i;k < length[index];k++) {
//              if ((contrast[k] != D_INF) && (residual[k] > 0.)) {
              if ((contrast[k] != D_INF) && (residual[k] > sqrt((double)(k - i + 1)) * ROUNDOFF_ERROR)) {
                contrast[k] -= ((double)(k - i + 1) / 2.) * (log(residual[k] /
                                 (k - i + 1)) + log(2 * M_PI) + 1);
/*                contrast[k] -= ((double)(k - i + 1) / 2.) * (log(residual[k] /
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
//          if (contrast[j] > 0.) {
          if (contrast[j] > sqrt((double)((nb_variable - 1) * (j - i + 1))) * ROUNDOFF_ERROR) {
            contrast[j] = -((double)((nb_variable - 1) * (j - i + 1)) / 2.) * (log(contrast[j] /
                             ((nb_variable - 1) * (j - i + 1))) + log(2 * M_PI) + 1);
/*            contrast[j] = -((double)((nb_variable - 1) * (j - i + 1)) / 2.) * (log(contrast[j] /
                             ((nb_variable - 1) * (j - i))) + log(2 * M_PI)) +
                           (double)((nb_variable - 1) * (j - i)) / 2.; */
          }
          else {
            contrast[j] = D_INF;
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
      if ((model_type[i - 1] == POISSON_CHANGE) || (model_type[i - 1] == GAUSSIAN_CHANGE) ||
          (model_type[i - 1] == MEAN_CHANGE) || (model_type[i - 1] == VARIANCE_CHANGE) ||
          (model_type[i - 1] == MEAN_VARIANCE_CHANGE)) {
        psegment = int_sequence[index][0] + 1;

        if (type[i] != REAL_VALUE) {
          pisequence = int_sequence[index][i];
          mean[i][0] = *pisequence++;
          j = 0;

          for (k = 1;k < length[index];k++) {
            if (*psegment != *(psegment - 1)) {
              mean[i][j] /= (k - j);
              for (m = j + 1;m < k;m++) {
                mean[i][m] = mean[i][j];
              }
              j = k;
              mean[i][j] = *pisequence++;
            }
            else {
              mean[i][j] += *pisequence++;
            }
            psegment++;
          }
        }

        else {
          prsequence = real_sequence[index][i];
          mean[i][0] = *prsequence++;
          j = 0;

          for (k = 1;k < length[index];k++) {
            if (*psegment != *(psegment - 1)) {
              mean[i][j] /= (k - j);
              for (m = j + 1;m < k;m++) {
                mean[i][m] = mean[i][j];
              }
              j = k;
              mean[i][j] = *prsequence++;
            }
            else {
              mean[i][j] += *prsequence++;
            }
            psegment++;
          }
        }

        mean[i][j] /= (length[index] - j);
        for (k = j + 1;k < length[index];k++) {
          mean[i][k] = mean[i][j];
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

    if (model_type[0] != MEAN_CHANGE) {
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
          os << "\n" << SEQ_label[SEQL_MAX_POSTERIOR_CHANGE_POINT_PROBABILITY] << "\n\n";
          break;
        case SEGMENT :
          os << "\n" << SEQ_label[SEQL_MAX_POSTERIOR_SEGMENT_PROBABILITY] << "\n\n";
          break;
        }
      }

      else {
        switch (output) {
        case CHANGE_POINT :
          os << "\n" << SEQ_label[SEQL_MAX_CHANGE_POINT_LIKELIHOOD] << "\n\n";
          break;
        case SEGMENT :
          os << "\n" << SEQ_label[SEQL_MAX_SEGMENT_LIKELIHOOD] << "\n\n";
          break;
        }
      }

      profile_ascii_print(os , index , nb_segment , backward_output ,
                          (output == CHANGE_POINT ? SEQ_label[SEQL_CHANGE_POINT] : SEQ_label[SEQL_SEGMENT]) ,
                          mean);

      os << "\n" << SEQ_label[SEQL_SEGMENTATION_LIKELIHOOD] << ": " << segmentation_likelihood;
      if (likelihood != D_INF) {
        os << "   (" << exp(segmentation_likelihood - likelihood) << ")";
      }
      os << endl;
      break;
    }

    case 's' : {
      if (likelihood != D_INF) {
        switch (output) {
        case CHANGE_POINT :
          os << "\n" << SEQ_label[SEQL_MAX_POSTERIOR_CHANGE_POINT_PROBABILITY] << "\n\n";
          break;
        case SEGMENT :
          os << "\n" << SEQ_label[SEQL_MAX_POSTERIOR_SEGMENT_PROBABILITY] << "\n\n";
          break;
        }
      }

      else {
        switch (output) {
        case CHANGE_POINT :
          os << "\n" << SEQ_label[SEQL_MAX_CHANGE_POINT_LIKELIHOOD] << "\n\n";
          break;
        case SEGMENT :
          os << "\n" << SEQ_label[SEQL_MAX_SEGMENT_LIKELIHOOD] << "\n\n";
          break;
        }
      }

      profile_spreadsheet_print(os , index , nb_segment , backward_output ,
                                (output == CHANGE_POINT ? SEQ_label[SEQL_CHANGE_POINT] : SEQ_label[SEQL_SEGMENT]) ,
                                mean);

      os << "\n" << SEQ_label[SEQL_SEGMENTATION_LIKELIHOOD] << "\t" << segmentation_likelihood;
      if (likelihood != D_INF) {
        os << "\t" << exp(segmentation_likelihood - likelihood);
      }
      os << endl;
      break;
    }

    case 'g' : {
      profile_plot_print(os , index , nb_segment , backward_output , mean);
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
        os << "\n" << SEQ_label[SEQL_AMBIGUITY] << ": " << ambiguity
           << " (" << ambiguity / length[index] << ")" << endl;
        break;
      case 's' :
        os << "\n" << SEQ_label[SEQL_AMBIGUITY] << "\t" << ambiguity
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

  delete [] sequence_mean;
  delete [] residual;

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
    delete [] mean[i];
  }
  delete [] mean;

  return segmentation_likelihood;
}


/*--------------------------------------------------------------*
 *
 *  Calcul des L segmentations optimales et des profils de segmentation d'une sequence.
 *
 *  arguments : reference sur un objet Format_error, stream,
 *              identificateur de la sequence, nombre de segments, types des modeles,
 *              type de sortie, format ('a' : ASCII, 's' : Spreadsheet),
 *              methode de calcul des segmentations (algorithme de programmation dynamique ou
 *              algorithme forward-backward de simulation), nombre de segmentations.
 *
 *--------------------------------------------------------------*/

bool Sequences::segment_profile_write(Format_error &error , ostream &os , int iidentifier ,
                                      int nb_segment , int *model_type , int output ,
                                      char format , int segmentation , int nb_segmentation) const

{
  bool status = true;
  register int i , j;
  int index = I_DEFAULT;
  double likelihood = D_INF , segmentation_likelihood , **rank;
  Sequences *seq;


  error.init();

  if (((index_parameter_type == TIME) && (index_interval->variance > 0.)) ||
      (index_parameter_type == POSITION)) {
    status = false;
    error.update(SEQ_error[SEQR_INDEX_PARAMETER_TYPE]);
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

        if (!marginal[i]) {
          status = false;
          ostringstream error_message;
          error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                        << STAT_error[STATR_MARGINAL_HISTOGRAM];
          error.update((error_message.str()).c_str());
        }

        else if (model_type[i] == MULTINOMIAL_CHANGE) {
          if ((marginal[i]->nb_value < 2) || (marginal[i]->nb_value > NB_OUTPUT)) {
            status = false;
            ostringstream error_message;
            error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                          << STAT_error[STATR_NB_VALUE];
            error.update((error_message.str()).c_str());
          }

          else {
            for (j = 0;j < marginal[i]->nb_value;j++) {
              if (marginal[i]->frequency[j] == 0) {
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

  if ((nb_segment < 2) || (nb_segment > (index == I_DEFAULT ? hlength->offset : length[index]) / 2)) {
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
        rank[i] = seq->marginal[i]->rank_computation();
      }
      else {
        rank[i] = 0;
      }
    }

    for (i = 0;i < seq->nb_sequence;i++) {
      if ((index == I_DEFAULT) || (index == i)) {
        if (model_type[0] != MEAN_CHANGE) {
          likelihood = seq->forward_backward(i , nb_segment , model_type , rank , &os ,
                                             output , format);
        }
        segmentation_likelihood = seq->forward_backward_dynamic_programming(i , nb_segment , model_type ,
                                                                            rank , os , output , format ,
                                                                            likelihood);
        if (segmentation_likelihood == D_INF) {
          status = false;
          error.update(SEQ_error[SEQR_SEGMENTATION_FAILURE]);
        }

        else {
          switch (segmentation) {
          case FORWARD_DYNAMIC_PROGRAMMING :
            seq->L_segmentation(i , nb_segment , model_type , rank , os , format ,
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
 *  Calcul des L segmentations optimales et des profils de segmentation d'une sequence et
 *  ecriture des resultats dans un fichier.
 *
 *  arguments : reference sur un objet Format_error, path, identificateur de la sequence,
 *              nombre de segments, types des modeles, type de sortie, format de fichier
 *              ('a' : ASCII, 's' : Spreadsheet), methode de calcul des segmentations
 *              (algorithme de programmation dynamique ou algorithme forward-backward
 *               de simulation), nombre de segmentations.
 *
 *--------------------------------------------------------------*/

bool Sequences::segment_profile_write(Format_error &error , const char *path ,
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
 *  Calcul des profils de segmentation d'une sequence et affichage
 *  des resultats au format Gnuplot.
 *
 *  arguments : reference sur un objet Format_error, prefixe des fichiers,
 *              identificateur de la sequence, nombre de segments, types des modeles,
 *              type de sortie, titre des figures.
 *
 *--------------------------------------------------------------*/

bool Sequences::segment_profile_plot_write(Format_error &error , const char *prefix ,
                                           int iidentifier , int nb_segment , int *model_type ,
                                           int output , const char *title) const

{
  bool status = true;
  register int i , j , k;
  int index;
  double likelihood = D_INF , segmentation_likelihood , **rank;
  Sequences *seq;
  ostringstream data_file_name[2];
  ofstream *out_data_file;


  error.init();

  if (((index_parameter_type == TIME) && (index_interval->variance > 0.)) ||
      (index_parameter_type == POSITION)) {
    status = false;
    error.update(SEQ_error[SEQR_INDEX_PARAMETER_TYPE]);
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

        if (!marginal[i]) {
          status = false;
          ostringstream error_message;
          error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                        << STAT_error[STATR_MARGINAL_HISTOGRAM];
          error.update((error_message.str()).c_str());
        }

        else if (model_type[i] == MULTINOMIAL_CHANGE) {
          if ((marginal[i]->nb_value < 2) || (marginal[i]->nb_value > NB_OUTPUT)) {
            status = false;
            ostringstream error_message;
            error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                          << STAT_error[STATR_NB_VALUE];
            error.update((error_message.str()).c_str());
          }

          else {
            for (j = 0;j < marginal[i]->nb_value;j++) {
              if (marginal[i]->frequency[j] == 0) {
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

  else if ((nb_segment < 2) || (nb_segment > length[index] / 2)) {
    status = false;
    error.update(SEQ_error[SEQR_NB_SEGMENT]);
  }

  if (status) {

    // ecriture des fichiers de donnees

    i = (model_type[0] == MEAN_CHANGE ? 0 : 1);
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
          rank[i] = seq->marginal[i]->rank_computation();
        }
        else {
          rank[i] = 0;
        }
      }

      if (model_type[0] != MEAN_CHANGE) {
        likelihood = seq->forward_backward(index , nb_segment , model_type , rank ,
                                           out_data_file , output , 'g');
        out_data_file->close();
        delete out_data_file;

        data_file_name[0] << prefix << 0 << ".dat";
        out_data_file = new ofstream((data_file_name[0].str()).c_str());
      }

#     ifdef DEBUG
      likelihood = D_INF;
#     endif

      segmentation_likelihood = seq->forward_backward_dynamic_programming(index , nb_segment , model_type ,
                                                                          rank , *out_data_file ,
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

          for (j = 1;j < seq->nb_variable;j++) {
            if ((model_type[j - 1] == POISSON_CHANGE) || (model_type[j - 1] == GAUSSIAN_CHANGE) ||
                (model_type[j - 1] == MEAN_CHANGE) || (model_type[j - 1] == VARIANCE_CHANGE) ||
                (model_type[j - 1] == MEAN_VARIANCE_CHANGE)) {
              out_file << "set title";
              if (title) {
                out_file << " \"" << title << "\"";
              }
              out_file << "\n\n";
              break;
            }
          }

          if (index_parameter) {
            if (seq->index_parameter[index][seq->length[index] - 1] - seq->index_parameter[index][0] < TIC_THRESHOLD) {
              out_file << "set xtics 0,1" << endl;
            }

            j = 2;
            for (k = 1;k < seq->nb_variable;k++) {
              if ((model_type[k - 1] == POISSON_CHANGE) || (model_type[k - 1] == GAUSSIAN_CHANGE) ||
                  (model_type[k - 1] == MEAN_CHANGE) || (model_type[k - 1] == VARIANCE_CHANGE) ||
                  (model_type[k - 1] == MEAN_VARIANCE_CHANGE)) {
                out_file << "plot [" << seq->index_parameter[index][0] << ":"
                         << seq->index_parameter[index][seq->length[index] - 1] << "] ["
                         << MIN(seq->min_value[k] , 0) << ":"
                         << MAX(seq->max_value[k] , seq->min_value[k] + 1) << "] "
                         << "\"" << label((data_file_name[0].str()).c_str()) << "\" using " << 1 << " : " << j++
                         << " title \"" << SEQ_label[SEQL_SEQUENCE] << "\" with linespoints" << ",\\" << endl;
                out_file << "\"" << label((data_file_name[0].str()).c_str()) << "\" using " << 1 << " : " << j++
                         << " title \"" << SEQ_label[SEQL_SEGMENT_MEAN] << "\" with linespoints" << endl;

                if (i == 0) {
                  out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
                }
                out_file << endl;
              }
            }

            if (likelihood != D_INF) {
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
                out_file << " \"" << title << " - ";
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
              if ((model_type[k - 1] == POISSON_CHANGE) || (model_type[k - 1] == GAUSSIAN_CHANGE) ||
                  (model_type[k - 1] == MEAN_CHANGE) || (model_type[k - 1] == VARIANCE_CHANGE) ||
                  (model_type[k - 1] == MEAN_VARIANCE_CHANGE)) {
                out_file << "plot [0:" << seq->length[index] - 1 << "] ["
                         << MIN(seq->min_value[k] , 0) << ":"
                         << MAX(seq->max_value[k] , seq->min_value[k] + 1) << "] "
                         << "\"" << label((data_file_name[0].str()).c_str()) << "\" using " << j++
                         << " title \"" << iidentifier << "\" with linespoints,\\" << endl;
                out_file << "\"" << label((data_file_name[0].str()).c_str()) << "\" using " << j++
                         << " title \"" << SEQ_label[SEQL_SEGMENT_MEAN] << "\" with linespoints" << endl;

                if (i == 0) {
                  out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
                }
                out_file << endl;
              }
            }

            if (likelihood != D_INF) {
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
                out_file << " \"" << title << " - ";
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
 *  Segmentation hierarchique d'une sequence.
 *
 *  arguments : reference sur un objet Format_error, stream,
 *              identificateur de la sequence, nombre de segments maximum,
 *              types des modeles.
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::hierarchical_segmentation(Format_error &error , ostream &os , int iidentifier ,
                                                int max_nb_segment , int *model_type) const

{
  bool status = true;
  register int i , j , k;
  int index , max_nb_value , nb_segment , segment_index , split_change_point ,
      begin_change_point , end_change_point , merge , *frequency , *pisequence , *psegment ,
      *nb_parameter , **change_point , ***begin_frequency , ***end_frequency;
  double sum , factorial_sum , diff , sum_square , buff , merge_contrast , max_likelihood ,
         *mean , *residual , *prsequence , *begin_contrast , *end_contrast , *likelihood ,
         *penalty , *penalized_likelihood , **rank , **factorial , **begin_sum , **end_sum ,
         **begin_factorial_sum , **end_factorial_sum , **begin_sum_square , **end_sum_square;
  Sequences *seq , *iseq;


  seq = 0;
  error.init();

  if (((index_parameter_type == TIME) && (index_interval->variance > 0.)) ||
      (index_parameter_type == POSITION)) {
    status = false;
    error.update(SEQ_error[SEQR_INDEX_PARAMETER_TYPE]);
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

        if (!marginal[i]) {
          status = false;
          ostringstream error_message;
          error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                        << STAT_error[STATR_MARGINAL_HISTOGRAM];
          error.update((error_message.str()).c_str());
        }

        else if (model_type[i] == MULTINOMIAL_CHANGE) {
          if ((marginal[i]->nb_value < 2) || (marginal[i]->nb_value > NB_OUTPUT)) {
            status = false;
            ostringstream error_message;
            error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                          << STAT_error[STATR_NB_VALUE];
            error.update((error_message.str()).c_str());
          }

          else {
            for (j = 0;j < marginal[i]->nb_value;j++) {
              if (marginal[i]->frequency[j] == 0) {
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
        rank[i] = marginal[i - 1]->rank_computation();
      }
      else {
        rank[i] = 0;
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
    begin_sum_square = new double*[seq->nb_variable];
    end_sum_square = new double*[seq->nb_variable];

    for (i = 1;i < seq->nb_variable;i++) {
      if ((model_type[i - 1] == MULTINOMIAL_CHANGE) && (seq->marginal[i]->nb_value > max_nb_value)) {
        max_nb_value = seq->marginal[i]->nb_value;
      }

      if (model_type[i - 1] == MULTINOMIAL_CHANGE) {
        begin_frequency[i] = new int*[seq->length[0]];
        end_frequency[i] = new int*[seq->length[0]];
        for (j = 0;j < seq->length[0];j++) {
          begin_frequency[i][j] = new int[seq->marginal[i]->nb_value];
          end_frequency[i][j] = new int[seq->marginal[i]->nb_value];
        }

        begin_sum[i] = 0;
        end_sum[i] = 0;
      }

      else {
        begin_sum[i] = new double[seq->length[0]];
        end_sum[i] = new double[seq->length[0]];

        begin_frequency[i] = 0;
        end_frequency[i] = 0;
      }

      if (model_type[i - 1] == POISSON_CHANGE) {
        begin_factorial_sum[i] = new double[seq->length[0]];
        end_factorial_sum[i] = new double[seq->length[0]];

        factorial[i] = new double[seq->length[0]];
      }

      else {
        begin_factorial_sum[i] = 0;
        end_factorial_sum[i] = 0;

        factorial[i] = 0;
      }

      if ((model_type[i - 1] == ORDINAL_GAUSSIAN_CHANGE) || (model_type[i - 1] == GAUSSIAN_CHANGE) ||
          (model_type[i - 1] == VARIANCE_CHANGE)){
        begin_sum_square[i] = new double[seq->length[0]];
        end_sum_square[i] = new double[seq->length[0]];
      }
      else {
        begin_sum_square[i] = 0;
        end_sum_square[i] = 0;
      }
    }

    if (max_nb_value > 0) {
      frequency = new int[max_nb_value];
    }
    else {
      frequency = 0;
    }

    mean = new double[seq->nb_variable];
    residual = new double[seq->length[0]];

    begin_contrast = new double[seq->length[0]];
    end_contrast = new double[seq->length[0]];

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
        for (j = 0;j < seq->marginal[i]->nb_value;j++) {
          frequency[j] = 0;
        }
        sum = 0.;

        pisequence = seq->int_sequence[0][i];
        frequency[*pisequence++]++;

        for (j = 0;j < seq->marginal[i]->nb_value;j++) {
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

          for (k = 0;k < seq->marginal[i]->nb_value;k++) {
            begin_frequency[i][j][k] = frequency[k];
          }

          if (begin_contrast[j] != D_INF) {
            begin_contrast[j] += sum;
          }

/*          frequency[*pisequence++]++;

          for (k = 0;k < seq->marginal[i]->nb_value;k++) {
            begin_frequency[i][j][k] = frequency[k];
          }

          if (begin_contrast[j] != D_INF) {
            for (k = 0;k < seq->marginal[i]->nb_value;k++) {
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
            begin_contrast[j] -= ((double)(j + 1) / 2.) * (log(residual[j] /
                                   (j + 1)) + log(2 * M_PI) + 1);
/*            begin_contrast[j] -= ((double)(j + 1) / 2.) * (log(residual[j] /
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
/*        for (j = 0;j < seq->marginal[i]->nb_value;j++) {
          frequency[j] = 0;
        }
        sum = 0.;

        pisequence = seq->int_sequence[0][i] + seq->length[0] - 1;
        frequency[*pisequence--]++;

        for (j = 0;j < seq->marginal[i]->nb_value;j++) {
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

          for (k = 0;k < seq->marginal[i]->nb_value;k++) {
            end_frequency[i][j][k] = frequency[k];
          }

          if (end_contrast[j] != D_INF) {
            end_contrast[j] += sum;
          }

          frequency[*pisequence--]++;

          for (k = 0;k < seq->marginal[i]->nb_value;k++) {
            end_frequency[i][j][k] = frequency[k];
          }

          if (end_contrast[j] != D_INF) {
            for (k = 0;k < seq->marginal[i]->nb_value;k++) {
              if (frequency[k] > 0) {
                end_contrast[j] += frequency[k] * log((double)frequency[k] / (double)(seq->length[0] - j));
              }
            }
          }
        } */

        for (j = seq->length[0] - 1;j > 0;j--) {
          for (k = 0;k < seq->marginal[i]->nb_value;k++) {
            end_frequency[i][j][k] = begin_frequency[i][seq->length[0] - 1][k] - begin_frequency[i][j - 1][k];
          }

          if (end_contrast[j] != D_INF) {
            for (k = 0;k < seq->marginal[i]->nb_value;k++) {
              if (end_frequency[i][j][k] > 0) {
                end_contrast[j] += end_frequency[i][j][k] * log((double)end_frequency[i][j][k] /
                                   (double)(seq->length[0] - j));
              }
            }
          }
        }

        for (k = 0;k < seq->marginal[i]->nb_value;k++) {
          end_frequency[i][0][k] = begin_frequency[i][seq->length[0] - 1][k];
        }

        if (end_contrast[0] != D_INF) {
          for (k = 0;k < seq->marginal[i]->nb_value;k++) {
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
            end_contrast[j] -= ((double)(seq->length[0] - j) / 2.) * (log(residual[j] /
                                 (seq->length[0] - j)) + log(2 * M_PI) + 1);
/*            end_contrast[j] -= ((double)(seq->length[0] - j) / 2.) * (log(residual[j] /
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
            for (j = 0;j < seq->marginal[i]->nb_value;j++) {
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

              for (k = 0;k < seq->marginal[i]->nb_value;k++) {
                begin_frequency[i][j][k] = frequency[k];
              }

              if (begin_contrast[j] != D_INF) {
                begin_contrast[j] += sum;
              }

/*              frequency[*pisequence++]++;

              for (k = 0;k < seq->marginal[i]->nb_value;k++) {
                begin_frequency[i][j][k] = frequency[k];
              }

              if (begin_contrast[j] != D_INF) {
                for (k = 0;k < seq->marginal[i]->nb_value;k++) {
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

          for (j = 0;j < seq->marginal[i]->nb_value;j++) {
            frequency[j] = 0;
          }
          sum = 0.;

          frequency[*pisequence++]++;

          for (j = 0;j < seq->marginal[i]->nb_value;j++) {
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

            for (k = 0;k < seq->marginal[i]->nb_value;k++) {
              begin_frequency[i][j][k] = frequency[k];
            }

            if (begin_contrast[j] != D_INF) {
              begin_contrast[j] += sum;
            }

/*            frequency[*pisequence++]++;

            for (k = 0;k < seq->marginal[i]->nb_value;k++) {
              begin_frequency[i][j][k] = frequency[k];
            }

            if (begin_contrast[j] != D_INF) {
              for (k = 0;k < seq->marginal[i]->nb_value;k++) {
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
                begin_contrast[j] -= ((double)(j - change_point[nb_segment][segment_index - 1] + 1) / 2.) * (log(residual[j] /
                                       (j - change_point[nb_segment][segment_index - 1] + 1)) + log(2 * M_PI) + 1);
/*                begin_contrast[j] -= ((double)(j - change_point[nb_segment][segment_index - 1] + 1) / 2.) * (log(residual[j] /
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
              begin_contrast[j] -= ((double)(j - split_change_point + 1) / 2.) * (log(residual[j] /
                                     (j - split_change_point + 1)) + log(2 * M_PI) + 1);
/*              begin_contrast[j] -= ((double)(j - split_change_point + 1) / 2.) * (log(residual[j] /
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
            for (j = 0;j < seq->marginal[i]->nb_value;j++) {
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

              for (k = 0;k < seq->marginal[i]->nb_value;k++) {
                end_frequency[i][j][k] = frequency[k];
              }

              if (end_contrast[j] != D_INF) {
                end_contrast[j] += sum;
              }

              frequency[*pisequence--]++;

              for (k = 0;k < seq->marginal[i]->nb_value;k++) {
                end_frequency[i][j][k] = frequency[k];
              }

              if (end_contrast[j] != D_INF) {
                for (k = 0;k < seq->marginal[i]->nb_value;k++) {
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

          for (j = 0;j < seq->marginal[i]->nb_value;j++) {
            frequency[j] = 0;
          }
          sum = 0.;

          frequency[*pisequence--]++;

          for (j = 0;j < seq->marginal[i]->nb_value;j++) {
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

            for (k = 0;k < seq->marginal[i]->nb_value;k++) {
              end_frequency[i][j][k] = frequency[k];
            }

            if (end_contrast[j] != D_INF) {
              end_contrast[j] += sum;
            }

            frequency[*pisequence--]++;

            for (k = 0;k < seq->marginal[i]->nb_value;k++) {
              end_frequency[i][j][k] = frequency[k];
            }

            if (end_contrast[j] != D_INF) {
              for (k = 0;k < seq->marginal[i]->nb_value;k++) {
                if (frequency[k] > 0) {
                  end_contrast[j] += frequency[k] * log((double)frequency[k] /
                                                        (double)(split_change_point - j));
                }
              }
            }
          } */

          if (end_change_point > split_change_point) {
            for (j = end_change_point - 1;j > split_change_point;j--) {
              for (k = 0;k < seq->marginal[i]->nb_value;k++) {
                end_frequency[i][j][k] = begin_frequency[i][change_point[nb_segment][segment_index + 1] - 1][k] -
                                         begin_frequency[i][j - 1][k];
              }

              if (end_contrast[j] != D_INF) {
                for (k = 0;k < seq->marginal[i]->nb_value;k++) {
                  if (end_frequency[i][j][k] > 0) {
                    end_contrast[j] += end_frequency[i][j][k] * log((double)end_frequency[i][j][k] /
                                                                    (double)(change_point[nb_segment][segment_index + 1] - j));
                  }
                }
              }
            }

            for (k = 0;k < seq->marginal[i]->nb_value;k++) {
              end_frequency[i][j][k] = begin_frequency[i][change_point[nb_segment][segment_index + 1] - 1][k];
            }

            if (end_contrast[j] != D_INF) {
              for (k = 0;k < seq->marginal[i]->nb_value;k++) {
                if (end_frequency[i][j][k] > 0) {
                  end_contrast[j] += end_frequency[i][j][k] * log((double)end_frequency[i][j][k] /
                                                                  (double)(change_point[nb_segment][segment_index + 1] - j));
                }
              }
            }
          }

          for (j = split_change_point - 1;j > change_point[nb_segment][segment_index - 1];j--) {
            for (k = 0;k < seq->marginal[i]->nb_value;k++) {
              end_frequency[i][j][k] = begin_frequency[i][split_change_point - 1][k] - begin_frequency[i][j - 1][k];
            }

            if (end_contrast[j] != D_INF) {
              for (k = 0;k < seq->marginal[i]->nb_value;k++) {
                if (end_frequency[i][j][k] > 0) {
                  end_contrast[j] += end_frequency[i][j][k] * log((double)end_frequency[i][j][k] /
                                                                  (double)(split_change_point - j));
                }
              }
            }
          }

          for (k = 0;k < seq->marginal[i]->nb_value;k++) {
            end_frequency[i][j][k] = begin_frequency[i][split_change_point - 1][k];
          }

          if (end_contrast[j] != D_INF) {
            for (k = 0;k < seq->marginal[i]->nb_value;k++) {
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
                end_contrast[j] -= ((double)(change_point[nb_segment][segment_index + 1] - j) / 2.) * (log(residual[j] /
                                     (change_point[nb_segment][segment_index + 1] - j)) + log(2 * M_PI) + 1);
/*                end_contrast[j] -= ((double)(change_point[nb_segment][segment_index + 1] - j) / 2.) * (log(residual[j] /
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
              end_contrast[j] -= ((double)(split_change_point - j) / 2.) * (log(residual[j] /
                                   (split_change_point - j)) + log(2 * M_PI) + 1);
/*              end_contrast[j] -= ((double)(split_change_point - j) / 2.) * (log(residual[j] /
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
              for (j = 0;j < seq->marginal[i]->nb_value;j++) {
                frequency[j] = begin_frequency[i][change_point[nb_segment][segment_index - 1] - 1][j] +
                               begin_frequency[i][split_change_point - 1][j];
              }

              if (merge_contrast != D_INF) {
                for (j = 0;j < seq->marginal[i]->nb_value;j++) {
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
              for (j = 0;j < seq->marginal[i]->nb_value;j++) {
                frequency[j] = end_frequency[i][split_change_point][j] +
                               end_frequency[i][change_point[nb_segment][segment_index]][j];
              }

              if (merge_contrast != D_INF) {
                for (j = 0;j < seq->marginal[i]->nb_value;j++) {
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

      os << "\n" << SEQ_label[SEQL_NB_SEGMENT] << " |  2 * " << STAT_label[STATL_LIKELIHOOD]
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


/*--------------------------------------------------------------*
 *
 *  Segmentation optimale d'une sequence.
 *
 *  arguments : reference sur un objet Format_error, identificateur de la sequence,
 *              nombre de segments, reference sur un objet Vector_distance,
 *              stream, type de sortie.
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::segmentation(Format_error &error , int iidentifier ,
                                   int nb_segment , const Vector_distance &ivector_dist ,
                                   ostream &os , int output) const

{
  bool status = true;
  register int i , j , k , m;
  int index , *psegment , *pisequence , **optimal_length;
  double sum , diff , buff , backward_min , exponent , *mean_square_diff , *contrast ,
         *prsequence , **rank , **forward , **backward , **backward_output , **mean;
  Vector_distance *vector_dist;
  Sequences *seq , *oseq;


  oseq = 0;
  error.init();

  if (((index_parameter_type == TIME) && (index_interval->variance > 0.)) ||
      (index_parameter_type == POSITION)) {
    status = false;
    error.update(SEQ_error[SEQR_INDEX_PARAMETER_TYPE]);
  }

  if (ivector_dist.nb_variable != nb_variable) {
    status = false;
    error.update(STAT_error[STATR_NB_VARIABLE]);
  }

  else {
    for (i = 0;i < nb_variable;i++) {
      if (ivector_dist.variable_type[i] != NUMERIC) {
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

          if (!marginal[i]) {
            status = false;
            ostringstream error_message;
            error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                          << STAT_error[STATR_MARGINAL_HISTOGRAM];
            error.update((error_message.str()).c_str());
          }

          if ((ivector_dist.variable_type[i] == SYMBOLIC) && ((max_value[i] >= NB_SYMBOL) ||
               ((ivector_dist.symbol_distance[i]) && (ivector_dist.nb_value[i] != max_value[i] + 1)))) {
            status = false;
             ostringstream error_message;
            error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                          << STAT_error[STATR_NB_SYMBOL];
            error.update((error_message.str()).c_str());
          }

          if ((ivector_dist.variable_type[i] == CIRCULAR) &&
              (max_value[i] - min_value[i] >= ivector_dist.period[i])) {
            status = false;
            ostringstream error_message;
            error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                          << STAT_error[STATR_NB_VALUE_PERIOD];
            error.update((error_message.str()).c_str());
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

  if ((nb_segment < 2) || (nb_segment > length[index] / 2)) {
    status = false;
    error.update(SEQ_error[SEQR_NB_SEGMENT]);
  }

  if (status) {
    seq = new Sequences(*this , 'a');

    vector_dist = new Vector_distance(ivector_dist);

    // calcul des rangs pour les variables ordinales

    rank = new double*[seq->nb_variable];

    for (i = 1;i < seq->nb_variable;i++) {
      if (vector_dist->variable_type[i - 1] == ORDINAL) {
        rank[i] = seq->marginal[i]->rank_computation();
      }
      else {
        rank[i] = 0;
      }

      // calcul des dispersions pour la standardisation

      if (seq->marginal[i]) {
        vector_dist->dispersion_computation(i - 1 , seq->marginal[i] , rank[i]);
      }

      else {
        switch (vector_dist->distance_type) {
        case ABSOLUTE_VALUE :
          vector_dist->dispersion[i - 1] = seq->mean_absolute_difference_computation(i);
          break;
        case QUADRATIC :
          vector_dist->dispersion[i - 1] = 2 * seq->variance_computation(i , seq->mean_computation(i));
          break;
        }

        if (vector_dist->dispersion[i - 1] == 0.) {
          vector_dist->dispersion[i - 1] = 1.;
        }
      }
    }

    mean_square_diff = new double[seq->length[index]];
    contrast = new double[seq->length[index]];

    forward = new double*[seq->length[index]];
    for (i = 0;i < seq->length[index];i++) {
      forward[i] = new double[nb_segment];
    }

    optimal_length = new int*[seq->length[index]];
    for (i = 0;i < seq->length[index];i++) {
      optimal_length[i] = new int[nb_segment];
    }

    backward = new double*[seq->length[index]];
    for (i = 0;i < seq->length[index];i++) {
      backward[i] = new double[nb_segment];
    }

    backward_output = new double*[seq->length[index]];
    for (i = 0;i < seq->length[index];i++) {
      backward_output[i] = new double[nb_segment];
    }

    mean = new double*[seq->nb_variable];
    for (i = 1;i < seq->nb_variable;i++) {
      if (vector_dist->variable_type[i - 1] == NUMERIC) {
        mean[i] = new double[seq->length[index]];
      }
      else {
        mean[i] = 0;
      }
    }

    // recurrence "forward"

    for (i = 0;i < seq->length[index];i++) {

      // calcul des residus

      mean_square_diff[i] = 0.;
      sum = 0.;

      for (j = i - 1;j >= 0;j--) {
        for (k = 1;k < seq->nb_variable;k++) {
          switch (vector_dist->variable_type[k - 1]) {

          case SYMBOLIC : {
            if (!vector_dist->symbol_distance[k - 1]) {
              diff = (seq->int_sequence[index][k][i] == seq->int_sequence[index][k][j] ? 0. : 1.);
            }
            else {
              diff = vector_dist->symbol_distance[k - 1][seq->int_sequence[index][k][i]][seq->int_sequence[index][k][j]];
            }
            break;
          }

          case ORDINAL : {
            diff = rank[k][seq->int_sequence[index][k][i]] - rank[k][seq->int_sequence[index][k][j]];
            break;
          }

          case NUMERIC : {
            if (seq->type[k] != REAL_VALUE) {
              diff = seq->int_sequence[index][k][i] - seq->int_sequence[index][k][j];
            }
            else {
              diff = seq->real_sequence[index][k][i] - seq->real_sequence[index][k][j];
            }
            break;
          }

          case CIRCULAR : {
            if (seq->int_sequence[index][k][i] <= seq->int_sequence[index][k][j]) {
              diff = MIN(seq->int_sequence[index][k][j] - seq->int_sequence[index][k][i] ,
                         seq->int_sequence[index][k][i] + vector_dist->period[k - 1] -
                         seq->int_sequence[index][k][j]);
            }
            else {
              diff = MIN(seq->int_sequence[index][k][i] - seq->int_sequence[index][k][j] ,
                         seq->int_sequence[index][k][j] + vector_dist->period[k - 1] -
                         seq->int_sequence[index][k][i]);
            }
            break;
          }
          }

          switch (vector_dist->distance_type) {
          case ABSOLUTE_VALUE :
            sum += vector_dist->weight[k - 1] * fabs(diff) / vector_dist->dispersion[k - 1];
//                   (seq->nb_variable == 2 ? 1. : vector_dist->dispersion[k - 1]);
            break;
          case QUADRATIC :
            sum += vector_dist->weight[k - 1] * diff * diff / vector_dist->dispersion[k - 1];
//                   (seq->nb_variable == 2 ? 1. : vector_dist->dispersion[k - 1]);
            break;
          }
        }

        mean_square_diff[j] += sum;
      }
 
/*      contrast[i] = 0.;
      for (j = i - 1;j >= 0;j--) {
        contrast[j] = mean_square_diff[j] / (i - j + 1);
      } */

      contrast[i] = -D_INF;
      for (j = i - 1;j >= 0;j--) {
        if (mean_square_diff[j] > 0.) {
          switch (vector_dist->distance_type) {
          case ABSOLUTE_VALUE :
            contrast[j] = (i - j + 1) * log(sqrt(2.) * mean_square_diff[j] /
                           ((i - j + 1) * (i - j + 1)));
            break;
          case QUADRATIC :
            contrast[j] = ((double)(i - j + 1) / 2.) * log(mean_square_diff[j] /
                           ((i - j + 1) * (i - j + 1)));
            break;
          }
        }
        else {
          contrast[j] = -D_INF;
        }
      }

#     ifdef DEBUG
      for (j = i - 1;j >= 0;j--) {
        cout << contrast[j] << "  ";
      }
      cout << endl;
#     endif

      for (j = 0;j < nb_segment;j++) {
        forward[i][j] = -D_INF;
      }

      for (j = MAX(0 , nb_segment + i - seq->length[index]);j < MIN((i < seq->length[index] - 1 ? nb_segment - 1 : nb_segment) , i + 1);j++) {
        if (j == 0) {
          forward[i][j] = contrast[0];
          if (forward[i][j] != -D_INF) {
            optimal_length[i][j] = i + 1;
          }
        }

        else {
          for (k = i;k >= j;k--) {
            if ((contrast[k] != -D_INF) && (forward[k - 1][j - 1] != -D_INF)) {
              buff = contrast[k] + forward[k - 1][j - 1];
              if (buff < forward[i][j]) {
                forward[i][j] = buff;
                optimal_length[i][j] = i - k + 1;
              }
            }
          }
        }
      }
    }

    if (forward[seq->length[index] - 1][nb_segment - 1] != -D_INF) {

      // restauration

      i = seq->length[index] - 1;
      psegment = seq->int_sequence[index][0] + i;

      for (j = nb_segment - 1;j >= 0;j--) {
        for (k = i;k > i - optimal_length[i][j];k--) {
          *psegment-- = j;
        }
        i -= optimal_length[i][j];
      }

#     ifdef DEBUG
      switch (vector_dist->distance_type) {
      case ABSOLUTE_VALUE :
        os << "\nsegment standardized mean absolute difference";
        break;
      case QUADRATIC :
        os << "\nsegment standardized mean square difference";
        break;
      }
      os << ": " << forward[seq->length[index] - 1][nb_segment - 1] << endl;
#     endif

      // recurrence "backward"

      for (i = seq->length[index] - 1;i >= 0;i--) {

        // calcul des residus

        mean_square_diff[i] = 0.;
        sum = 0.;

        for (j = i + 1;j < seq->length[index];j++) {
          for (k = 1;k < seq->nb_variable;k++) {
            switch (vector_dist->variable_type[k - 1]) {

            case SYMBOLIC : {
              if (!vector_dist->symbol_distance[k - 1]) {
                diff = (seq->int_sequence[index][k][i] == seq->int_sequence[index][k][j] ? 0. : 1.);
              }
              else {
                diff = vector_dist->symbol_distance[k - 1][seq->int_sequence[index][k][i]][seq->int_sequence[index][k][j]];
              }
              break;
            }

            case ORDINAL : {
              diff = rank[k][seq->int_sequence[index][k][i]] - rank[k][seq->int_sequence[index][k][j]];
              break;
            }

            case NUMERIC : {
              if (seq->type[k] != REAL_VALUE) {
                diff = seq->int_sequence[index][k][i] - seq->int_sequence[index][k][j];
              }
              else {
                diff = seq->real_sequence[index][k][i] - seq->real_sequence[index][k][j];
              }
              break;
            }

            case CIRCULAR : {
              if (seq->int_sequence[index][k][i] <= seq->int_sequence[index][k][j]) {
                diff = MIN(seq->int_sequence[index][k][j] - seq->int_sequence[index][k][i] ,
                           seq->int_sequence[index][k][i] + vector_dist->period[k - 1] -
                           seq->int_sequence[index][k][j]);
              }
              else {
                diff = MIN(seq->int_sequence[index][k][i] - seq->int_sequence[index][k][j] ,
                           seq->int_sequence[index][k][j] + vector_dist->period[k - 1] -
                           seq->int_sequence[index][k][i]);
              }
              break;
            }
            }

            switch (vector_dist->distance_type) {
            case ABSOLUTE_VALUE :
              sum += vector_dist->weight[k - 1] * fabs(diff) / vector_dist->dispersion[k - 1];
//                     (seq->nb_variable == 2 ? 1. : vector_dist->dispersion[k - 1]);
              break;
            case QUADRATIC :
              sum += vector_dist->weight[k - 1] * diff * diff / vector_dist->dispersion[k - 1];
//                     (seq->nb_variable == 2 ? 1. : vector_dist->dispersion[k - 1]);
              break;
            }
          }

          mean_square_diff[j] += sum;
        }
 
/*        contrast[i] = 0.;
        for (j = i + 1;j < seq->length[index];j++) {
          contrast[j] = mean_square_diff[j] / (j - i + 1);
        } */

        contrast[i] = -D_INF;
        for (j = i + 1;j < seq->length[index];j++) {
          if (mean_square_diff[j] > 0.) {
            switch (vector_dist->distance_type) {
            case ABSOLUTE_VALUE :
              contrast[j] = (j - i + 1) * log(sqrt(2.) * mean_square_diff[j] /
                             ((j - i + 1) * (j - i + 1)));
              break;
            case QUADRATIC :
              contrast[j] = ((double)(j - i + 1) / 2.) * log(mean_square_diff[j] /
                             ((j - i + 1) * (j - i + 1)));
              break;
            }
          }
          else {
            contrast[j] = -D_INF;
          }
        }

        for (j = 0;j < nb_segment;j++) {
          backward_output[i][j] = -D_INF;
        }

        for (j = MAX((i == 0 ? 0 : 1) , nb_segment + i - seq->length[index]);j < MIN(nb_segment , i + 1);j++) {
          if (j < nb_segment - 1) {
            backward[i][j] = -D_INF;
            for (k = seq->length[index] + j - nb_segment;k >= i;k--) {
              if ((contrast[k] != -D_INF) && (backward[k + 1][j + 1] != -D_INF)) {
                buff = contrast[k] + backward[k + 1][j + 1];
                if (buff < backward[i][j]) {
                  backward[i][j] = buff;
                }
              }

              if ((output == SEGMENT) && (k > i) && (backward[i][j] != -D_INF)) {
                if (i == 0) {
                  if (backward[i][j] < backward_output[k][j]) {
                    backward_output[k][j] = backward[i][j];
                  }
                }
                else if (forward[i - 1][j - 1] != -D_INF) {
                  buff = forward[i - 1][j - 1] + backward[i][j];
                  if (buff < backward_output[k][j]) {
                    backward_output[k][j] = buff;
                  }
                }
              }
            }
          }

          else {
            backward[i][j] = contrast[seq->length[index] - 1];

            if ((output == SEGMENT) && (forward[i - 1][j - 1] != D_INF) &&
                (backward[i][j] != -D_INF)) {
              buff = forward[i - 1][j - 1] + backward[i][j];
              for (k = seq->length[index] - 1;k > i;k--) {
                if (buff < backward_output[k][j]) {
                  backward_output[k][j] = buff;
                }
              }
            }
          }

          if (backward[i][j] != -D_INF) {
            if (i == 0) {
              backward_output[i][j] = backward[i][j];
            }
            else if (forward[i - 1][j - 1] != -D_INF) {
              backward_output[i][j] = forward[i - 1][j - 1] + backward[i][j];
            }
          }
        }
      }

#     ifdef DEBUG
      cout << "\n";
      for (i = 1;i < seq->length[index];i++) {
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
#     endif

      // restauration

#     ifdef MESSAGE
      if (output == SEGMENT) {
        int optimal_segment;

        psegment = seq->int_sequence[index][0];

        for (i = 0;i < seq->length[index];i++) {
          backward_min = -D_INF;
          for (j = 0;j < nb_segment;j++) {
            if (backward_output[i][j] < backward_min) {
              backward_min = backward_output[i][j];
              optimal_segment = j;
            }
          }

          if (optimal_segment != *psegment) {
            cout << "\nERROR: " << i << " | " << *psegment << " " << optimal_segment << endl;
          }

          psegment++;
        }
      }
#     endif

      for (i = 1;i < seq->nb_variable;i++) {
        if (vector_dist->variable_type[i - 1] == NUMERIC) {
          psegment = seq->int_sequence[index][0] + 1;

          if (seq->type[i] != REAL_VALUE) {
            pisequence = seq->int_sequence[index][i];
            mean[i][0] = *pisequence++;
            j = 0;

            for (k = 1;k < seq->length[index];k++) {
              if (*psegment != *(psegment - 1)) {
                mean[i][j] /= (k - j);
                for (m = j + 1;m < k;m++) {
                  mean[i][m] = mean[i][j];
                }
                j = k;
                mean[i][j] = *pisequence++;
              }
              else {
                mean[i][j] += *pisequence++;
              }
              psegment++;
            }
          }

          else {
            prsequence = seq->real_sequence[index][i];
            mean[i][0] = *pisequence++;
            j = 0;

            for (k = 1;k < seq->length[index];k++) {
              if (*psegment != *(psegment - 1)) {
                mean[i][j] /= (k - j);
                for (m = j + 1;m < k;m++) {
                  mean[i][m] = mean[i][j];
                }
                j = k;
                mean[i][j] = *prsequence++;
              }
              else {
                mean[i][j] += *prsequence++;
              }
              psegment++;
            }
          }

          mean[i][j] /= (seq->length[index] - j);
          for (k = j + 1;k < seq->length[index];k++) {
            mean[i][k] = mean[i][j];
          }
        }
      }

#     ifdef MESSAGE
      if ((backward[0][0] < forward[seq->length[index] - 1][nb_segment - 1] - DOUBLE_ERROR) ||
          (backward[0][0] > forward[seq->length[index] - 1][nb_segment - 1] + DOUBLE_ERROR)) {
        cout << "\nERROR: " << backward[0][0] << " | " << forward[seq->length[index] - 1][nb_segment - 1] << endl;
      }
      if ((backward_output[0][0] < backward[0][0] - DOUBLE_ERROR) ||
          (backward_output[0][0] > backward[0][0] + DOUBLE_ERROR)) {
        cout << "\nERROR: " << backward_output[0][0] << " | " << backward[0][0] << endl;
      }
#     endif

/*      if (forward[seq->length[index] - 1][nb_segment - 1] > 0.) {
        switch (vector_dist->distance_type) {
        case ABSOLUTE_VALUE :
          exponent = -seq->length[index];
          break;
        case QUADRATIC :
          exponent = -((double)seq->length[index] / 2.);
          break;
        }

        for (i = 0;i < seq->length[index];i++) {
          for (j = 0;j < nb_segment;j++) {
            if (backward_output[i][j] > 0.) {
              backward_output[i][j] = pow(backward_output[i][j] / forward[seq->length[index] - 1][nb_segment - 1] ,
                                          exponent);
              backward_output[i][j] = exp(exponent * log(backward_output[i][j] /
                                             forward[seq->length[index] - 1][nb_segment - 1]));
            }
            else {
              backward_output[i][j] = 0.;
            }
          }
        } */

      if (forward[seq->length[index] - 1][nb_segment - 1] != -D_INF) {
        for (i = 0;i < seq->length[index];i++) {
          for (j = 0;j < nb_segment;j++) {
            if (backward_output[i][j] != -D_INF) {
              backward_output[i][j] = exp(-backward_output[i][j] + forward[seq->length[index] - 1][nb_segment - 1]);
            }
            else {
              backward_output[i][j] = 0.;
            }
          }
        }

        switch (output) {
        case CHANGE_POINT :
          os << "\n" << SEQ_label[SEQL_MAX_CHANGE_POINT_LIKELIHOOD] << "\n\n";
          break;
        case SEGMENT :
          os << "\n" << SEQ_label[SEQL_MAX_SEGMENT_LIKELIHOOD] << "\n\n";
          break;
        }

        seq->profile_ascii_print(os , index , nb_segment , backward_output ,
                                 (output == CHANGE_POINT ? SEQ_label[SEQL_CHANGE_POINT] : SEQ_label[SEQL_SEGMENT]) ,
                                 mean);
      }

      oseq = new Sequences(*seq , 1 , &index);

      oseq->min_value[0] = 0;
      oseq->max_value[0] = nb_segment - 1;
      oseq->build_marginal_histogram(0);
    }

    delete vector_dist;

    for (i = 1;i < seq->nb_variable;i++) {
      delete [] rank[i];
    }
    delete [] rank;

    delete [] mean_square_diff;
    delete [] contrast;

    for (i = 0;i < seq->length[index];i++) {
      delete [] forward[i];
    }
    delete [] forward;

    for (i = 0;i < seq->length[index];i++) {
      delete [] backward[i];
    }
    delete [] backward;

    for (i = 0;i < seq->length[index];i++) {
      delete [] backward_output[i];
    }
    delete [] backward_output;

    for (i = 0;i < seq->length[index];i++) {
      delete [] optimal_length[i];
    }
    delete [] optimal_length;

    for (i = 1;i < seq->nb_variable;i++) {
      delete [] mean[i];
    }
    delete [] mean;

    delete seq;
  }

  return oseq;
}
