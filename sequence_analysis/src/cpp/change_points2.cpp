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
extern int column_width(int min_value , int max_value);
extern int column_width(int nb_value , const double *value , double scale = 1.);

extern void cumul_computation(int nb_value , const double *pmass , double *pcumul);
extern int cumul_method(int nb_value , const double *cumul , double scale = 1.);

extern char* label(const char *file_name);


#if defined (SYSTEM_IS__CYGWIN)
#define expl exp
#endif



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


/*--------------------------------------------------------------*/
/**
 *  Ecriture des profils de segments/d'etats, des profils de ruptures et
 *  des profils d'entropies pour une sequence.
 *
 *  arguments : stream, indice de la sequence, nombre de segments/d'etats,
 *              pointeur sur les profils de segments/d'etats, label,
 *              pointeur sur les moyennes, les profils de ruptures et
 *              les profils d'entropies.
 */
/*--------------------------------------------------------------*/

ostream& Sequences::profile_ascii_print(ostream &os , int index , int nb_segment ,
                                        double **profiles , const char *label ,
                                        double **mean , long double **change_point ,
                                        long double **begin_conditonal_entropy ,
                                        long double **end_conditional_entropy ,
                                        long double **change_point_entropy) const

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

  if ((begin_conditonal_entropy) && (end_conditional_entropy) &&
      (change_point_entropy)) {
    width[2 * nb_variable + 3] = 0;

    for (i = 1;i < nb_segment;i++) {
      buff = column_width(length[index] , begin_conditonal_entropy[i]);
      if (buff > width[2 * nb_variable + 3]) {
        width[2 * nb_variable + 3] = buff;
      }
    }

    for (i = 1;i < nb_segment;i++) {
      buff = column_width(length[index] , end_conditional_entropy[i]);
      if (buff > width[2 * nb_variable + 3]) {
        width[2 * nb_variable + 3] = buff;
      }
    }

    for (i = 1;i < nb_segment;i++) {
      buff = column_width(length[index] , change_point_entropy[i]);
      if (buff > width[2 * nb_variable + 3]) {
        width[2 * nb_variable + 3] = buff;
      }
    }

    width[2 * nb_variable + 3] += ASCII_SPACE;
  }

  if (!change_point) {
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

  if ((begin_conditonal_entropy) && (end_conditional_entropy) &&
      (change_point_entropy)) {
    os << "\n" << SEQ_label[SEQL_BEGIN_CONDITIONAL_ENTROPY] << ", "
       << SEQ_label[SEQL_END_CONDITIONAL_ENTROPY] << ", "
       << SEQ_label[SEQL_CHANGE_POINT_ENTROPY] << endl;

    os << "\n";
    if (index_parameter_type == TIME) {
      os << SEQ_label[SEQL_TIME];
    }
    else {
      os << SEQ_label[SEQL_INDEX];
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

    for (i = 0;i < length[index];i++) {
      os.setf(ios::right , ios::adjustfield);
      os << setw(buff) << (index_parameter ? index_parameter[index][i] : i) << "  ";

      os.setf(ios::left , ios::adjustfield);
      for (j = 1;j < nb_segment;j++) {
        os << setw(width[2 * nb_variable + 3]) << begin_conditonal_entropy[j][i];
      }
      os << "   ";
      for (j = 1;j < nb_segment;j++) {
        os << setw(width[2 * nb_variable + 3]) << end_conditional_entropy[j][i];
      }
      os << "   ";
      for (j = 1;j < nb_segment;j++) {
        os << setw(width[2 * nb_variable + 3]) << change_point_entropy[j][i];
      }

      if (i == 0) {
        os.setf(ios::right , ios::adjustfield);
        os << setw(width[nb_variable + 2]) << identifier[index];
      }
      os << endl;
    }
  }

  delete [] width;

  os.setf((FMTFLAGS)old_adjust , ios::adjustfield);

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture des profils de segments/d'etats, des profils de ruptures et
 *  des profils d'entropies pour une sequence au format tableur.
 *
 *  arguments : stream, indice de la sequence, nombre de segments/d'etats,
 *              pointeur sur les profils de segments/d'etats, label,
 *              pointeur sur les moyennes, les profils de ruptures et
 *              les profils d'entropies.
 *
 *--------------------------------------------------------------*/

ostream& Sequences::profile_spreadsheet_print(ostream &os , int index , int nb_segment ,
                                              double **profiles , const char *label ,
                                              double **mean , long double **change_point ,
                                              long double **begin_conditonal_entropy ,
                                              long double **end_conditional_entropy ,
                                              long double **change_point_entropy) const

{
  register int i , j , k;


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

    // pour les donnees manquantes

    if ((index_parameter) && (index_interval->variance > 0.) && (i < length[index] - 1)) {
      for (k = index_parameter[index][i] + index_interval->offset;k < index_parameter[index][i + 1];k += index_interval->offset) {
        if (!change_point) {
          os << int_sequence[index][0][i] << "\t";
        }
        for (j = 1;j < nb_variable;j++) {
          if ((mean) && (mean[j])) {
            os << mean[j][i] << "\t";
          }
          os << "\t";
        }

        os << k;
        for (j = 0;j < nb_segment;j++) {
          os << "\t" << profiles[i][j];
        }

        if (change_point) {
          os << "\t";
          for (j = 1;j < nb_segment;j++) {
            os << "\t" << change_point[j][i];
          }
        }
        os << endl;
      }
    }
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
    if (index_parameter_type == TIME) {
      os << SEQ_label[SEQL_TIME];
    }
    else {
      os << SEQ_label[SEQL_INDEX];
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

    for (i = 0;i < length[index];i++) {
      os << (index_parameter ? index_parameter[index][i] : i);

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

      if (i == 0) {
        os << "\t" << identifier[index];
      }
      os << endl;

      // pour les donnees manquantes

      if ((index_parameter) && (index_interval->variance > 0.) && (i < length[index] - 1)) {
        for (k = index_parameter[index][i] + index_interval->offset;k < index_parameter[index][i + 1];k += index_interval->offset) {
          os << k;

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
          os << endl;
        }
      }
    }
  }

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture des profils de segments/d'etats, des profils de ruptures et
 *  des profils d'entropies pour une sequence au format Gnuplot.
 *
 *  arguments : stream, indice de la sequence, nombre de segments/d'etats,
 *              pointeur sur les profils de segments/d'etats, pointeur sur les moyennes,
 *              les profils de rupture et les profils d'entropies.
 *
 *--------------------------------------------------------------*/

ostream& Sequences::profile_plot_print(ostream &os , int index , int nb_segment ,
                                       double **profiles , double **mean ,
                                       long double **change_point ,
                                       long double **begin_conditonal_entropy ,
                                       long double **end_conditional_entropy ,
                                       long double **change_point_entropy) const

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

    if ((change_point) && (begin_conditonal_entropy) &&
        (end_conditional_entropy) && (change_point_entropy)) {
      for (j = 1;j < nb_segment;j++) {
        os << change_point[j][i] << " ";
      }
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

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture des profils de ruptures ou d'entropie pour une sequence au format "plotable".
 *
 *  arguments : reference sur un objet MultiPlot, indice de la sequence,
 *              nombre de segments, pointeur sur les profils de rupture.
 *
 *--------------------------------------------------------------*/

void Sequences::change_point_profile_plotable_write(MultiPlot &plot , int index , int nb_segment ,
                                                    long double **change_point) const

{
  register int i , j , k;


  plot.resize(MAX(nb_segment - 1 , 3));
  i = 0;

  if (index_parameter) {
    for (j = MAX(1 , nb_segment - 3);j < nb_segment;j++) {
      for (k = 0;k < length[index];k++) {
        plot[i].add_point(index_parameter[index][k] , change_point[j][k]);
      }
      i++;
    }
  }

  else {
    for (j =  MAX(1 , nb_segment - 3);j < nb_segment;j++) {
      for (k = 0;k < length[index];k++) {
        plot[i].add_point(k , change_point[j][k]);
      }
      i++;
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Ecriture des profils d'entropies au format "plotable".
 *
 *  arguments : reference sur un objet MultiPlot, indice de la sequence,
 *              pointeurs sur les profils d'entropies.
 *
 *--------------------------------------------------------------*/

void Sequences::entropy_profile_plotable_write(MultiPlot &plot , int index ,
                                               long double *begin_conditional_entropy ,
                                               long double *end_conditional_entropy ,
                                               long double *change_point_entropy) const

{
  register int i;


  plot.resize(3);

  if (index_parameter) {
    for (i = 0;i < length[index];i++) {
      plot[0].add_point(index_parameter[index][i] , begin_conditional_entropy[i]);
      plot[1].add_point(index_parameter[index][i] , end_conditional_entropy[i]);
      plot[2].add_point(index_parameter[index][i] , change_point_entropy[i]);
    }
  }

  else {
    for (i = 0;i < length[index];i++) {
      plot[0].add_point(i , begin_conditional_entropy[i]);
      plot[1].add_point(i , end_conditional_entropy[i]);
      plot[2].add_point(i , change_point_entropy[i]);
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Calcul par sommation des profils de segments/ruptures d'une sequence et
 *  des profils d'entropie.
 *
 *  arguments : indice de la sequence, nombre de segments, types des modeles,
 *              rangs (variables ordinales), stream, pointeur sur un objet MultiPlotSet,
 *              type de sortie, format de sortie ('a' : ASCII, 's' : Spreadsheet,
 *              'g' : Gnuplot, 'p' : plotable).
 *
 *--------------------------------------------------------------*/

double Sequences::forward_backward(int index , int nb_segment , int *model_type ,
                                   double **rank , ostream *os , MultiPlotSet *plot_set ,
                                   int output , char format) const

{
  register int i , j , k , m;
  int max_nb_value , *frequency , *pisequence;
  double sum , factorial_sum , diff , buff , rlikelihood , *mean , *prsequence ,
         *likelihood , **factorial , **backward_output , ***smoothed;
  long double sum_square , segment_norm , sequence_norm , lbuff , lsum ,
              segmentation_entropy , first_order_entropy , change_point_entropy_sum ,
              marginal_entropy , *residual , *contrast , *normalized_contrast , *norm ,
              *forward_norm , *backward_norm , *entropy_smoothed , *segment_predicted ,
              **forward , **backward , **change_point , **forward_predicted_entropy ,
              **backward_predicted_entropy , **forward_partial_entropy ,
              **backward_partial_entropy , **change_point_entropy , ***state_entropy;

# ifdef DEBUG
  long double *entropy_norm;
# endif


  max_nb_value = 0;
  factorial = new double*[nb_variable];

  for (i = 1;i < nb_variable;i++) {
    if ((model_type[i - 1] == MULTINOMIAL_CHANGE) &&
        (marginal_distribution[i]->nb_value > max_nb_value)) {
      max_nb_value = marginal_distribution[i]->nb_value;
    }

    if (model_type[i - 1] == POISSON_CHANGE) {
      factorial[i] = new double[length[index]];
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
  residual = new long double[length[index]];

  contrast = new long double[length[index]];
  normalized_contrast = new long double[length[index]];

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
//  entropy_norm = new long double[length[index]];

  likelihood = new double[nb_segment];

  backward = new long double*[length[index]];
  for (i = 0;i < length[index];i++) {
    backward[i] = new long double[nb_segment];
  }

  backward_predicted_entropy = new long double*[length[index]];
  for (i = 0;i < length[index];i++) {
    backward_predicted_entropy[i] = new long double[nb_segment];
  }

  backward_norm = new long double[length[index]];

  smoothed = new double**[nb_segment];
  for (i = 1;i < nb_segment;i++) {
    smoothed[i] = new double*[length[index]];
    for (j = 0;j < length[index];j++) {
      smoothed[i][j] = new double[nb_segment];
    }
  }

  backward_output = new double*[length[index]];
  for (i = 0;i < length[index];i++) {
    backward_output[i] = new double[nb_segment];
  }

  change_point = new long double*[nb_segment];
  for (i = 1;i < nb_segment;i++) {
    change_point[i] = new long double[length[index]];
  }

  entropy_smoothed = new long double[nb_segment];

  state_entropy = new long double**[nb_segment];
  for (i = 1;i < nb_segment;i++) {
    state_entropy[i] = new long double*[length[index]];
    for (j = 0;j < length[index];j++) {
      state_entropy[i][j] = new long double[nb_segment];
    }
  }

  forward_partial_entropy = new long double*[nb_segment];
  for (i = 1;i < nb_segment;i++) {
    forward_partial_entropy[i] = new long double[length[index]];
  }

  backward_partial_entropy = new long double*[nb_segment];
  for (i = 1;i < nb_segment;i++) {
    backward_partial_entropy[i] = new long double[length[index]];
  }

  change_point_entropy = new long double*[nb_segment];
  for (i = 1;i < nb_segment;i++) {
    change_point_entropy[i] = new long double[length[index]];
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
        for (k = 0;k < marginal_distribution[j]->nb_value;k++) {
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
              contrast[k] -= ((double)(i - k + 1) / 2.) * (logl(residual[k] /
                               (i - k + 1)) + log(2 * M_PI) + 1);
/*              contrast[k] -= ((double)(i - k + 1) / 2.) * (logl(residual[k] /
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
          contrast[j] = -((double)((nb_variable - 1) * (i - j + 1)) / 2.) * (logl(contrast[j] /
                           ((nb_variable - 1) * (i - j + 1))) + log(2 * M_PI) + 1);
/*          contrast[j] = -((double)((nb_variable - 1) * (i - j + 1)) / 2.) * (logl(contrast[j] /
                           ((nb_variable - 1) * (i - j))) + log(2 * M_PI)) +
                         (double)((nb_variable - 1) * (i - j)) / 2.; */
        }
        else {
          contrast[j] = D_INF;
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

#     ifdef DEBUG
      if (i == length[index] - 1) {
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
      forward_predicted_entropy[i][j] = 0.;
    }
    norm[i] = 0.;

//    for (j = MAX(0 , nb_segment + i - length[index]);j < MIN((i < length[index] - 1 ? nb_segment - 1 : nb_segment) , i + 1);j++) {
    for (j = 0;j < MIN((i < length[index] - 1 ? nb_segment - 1 : nb_segment) , i + 1);j++) {
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
//      for (j = MAX(0 , nb_segment + i - length[index]);j < MIN((i < length[index] - 1 ? nb_segment - 1 : nb_segment) , i + 1);j++) {
      for (j = 0;j < MIN((i < length[index] - 1 ? nb_segment - 1 : nb_segment) , i + 1);j++) {
        forward[i][j] /= norm[i];
      }

      norm[i] = logl(norm[i]);
    }
//    entropy_norm[i] = norm[i];

    forward_norm[i] = segment_norm + norm[i];
  }

# ifdef DEBUG
  cout << "\n";
  for (i = 0;i < length[index];i++) {
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
  for (i = 0;i < length[index];i++) {
    cout << i << " |";
    for (j = 0;j < nb_segment;j++) {
      cout << " " << forward_predicted_entropy[i][j];
    }
    cout << endl;
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

#   ifdef MESSAGE
    segmentation_entropy = rlikelihood;
#   endif

    for (i = 1;i < nb_segment;i++) {
      for (j = 0;j < length[index];j++) {
        for (k = 0;k < nb_segment;k++) {
          state_entropy[i][j][k] = 0.;
        }
      }
    }

    // recurrence "backward"

    for (i = length[index] - 1;i >= 0;i--) {

      // calcul des log-vraisemblances des segments

      for (j = i;j < length[index];j++) {
        contrast[j] = 0.;
      }

      for (j = 1;j < nb_variable;j++) {
        if (model_type[j - 1] == MULTINOMIAL_CHANGE) {
          for (k = 0;k < marginal_distribution[j]->nb_value;k++) {
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
                contrast[k] -= ((double)(k - i + 1) / 2.) * (logl(residual[k] /
                                 (k - i + 1)) + log(2 * M_PI) + 1);
/*                contrast[k] -= ((double)(k - i + 1) / 2.) * (logl(residual[k] /
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
            contrast[j] = -((double)((nb_variable - 1) * (j - i + 1)) / 2.) * (logl(contrast[j] /
                             ((nb_variable - 1) * (j - i + 1))) + log(2 * M_PI) + 1);
/*            contrast[j] = -((double)((nb_variable - 1) * (j - i + 1)) / 2.) * (logl(contrast[j] /
                             ((nb_variable - 1) * (j - i))) + log(2 * M_PI)) +
                           (double)((nb_variable - 1) * (j - i)) / 2.; */
          }
          else {
            contrast[j] = D_INF;
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
        backward_output[i][j] = 0.;

        for (k = 1;k < nb_segment;k++) {
          smoothed[k][i][j] = 0.;
        }
      }
      norm[i] = 0.;

//      for (j = MAX((i == 0 ? 0 : 1) , nb_segment + i - length[index]);j < MIN(nb_segment , i + 1);j++) {
      for (j = MAX((i == 0 ? 0 : 1) , nb_segment + i - length[index]);j < nb_segment;j++) {
        if (j < nb_segment - 1) {
          for (k = i;k <= length[index] + j - nb_segment;k++) {
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
//        for (j = MAX((i == 0 ? 0 : 1) , nb_segment + i - length[index]);j < MIN(nb_segment , i + 1);j++) {
        for (j = MAX((i == 0 ? 0 : 1) , nb_segment + i - length[index]);j < nb_segment;j++) {
          backward[i][j] /= norm[i];
        }

        norm[i] = logl(norm[i]);
      }

      backward_norm[i] = segment_norm + norm[i];

      // extraction des probabilitees lissees pour les differents nombres de segments possibles

      if (i < length[index] - 1) {
        for (j = 1;j < nb_segment;j++) {
          sequence_norm = expl(forward_norm[i] + backward_norm[i + 1] - likelihood[j]);

          for (k = MAX(0 , j + i - length[index] - 1);k <= MIN(j , i);k++) {
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

      // calcul des probabilites a posteriori des ruptures
      // pour les differents nombres de segments possibles

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
        for (j = MAX(1 , nb_segment + i - length[index]);j < MIN(nb_segment , i + 1);j++) {
          if (output == CHANGE_POINT) {
            backward_output[i][j] = forward[i - 1][j - 1] * backward[i][j] * sequence_norm;
          }
          change_point[nb_segment - 1][i] += forward[i - 1][j - 1] * backward[i][j];
        }
        change_point[nb_segment - 1][i] *= sequence_norm;

        for (j = 1;j < nb_segment - 1;j++) {
          change_point[j][i] = 0.;
          for (k = MAX(1 , j + 1 + i - length[index]);k <= MIN(j , i);k++) {
            change_point[j][i] += forward[i - 1][k - 1] * backward[i][k + nb_segment - j - 1];
          }
          change_point[j][i] *= expl(forward_norm[i - 1] + backward_norm[i] - likelihood[j]);
        }
      }

      // calcul des entropies partielles pour les differents nombres de segments possibles

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

      if (i > 0) {
        for (j = 1;j < nb_segment;j++) {
          sequence_norm = expl(forward_norm[i - 1] + backward_norm[i] - likelihood[j]);

          for (k = MAX(1 , j + i - length[index] - 1);k <= MIN(j , i);k++) {
            if (k < j) {
              lsum = 0.;
              for (m = length[index] + k - nb_segment;m >= i;m--) {
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
              lsum = forward[i - 1][k - 1] * normalized_contrast[length[index] - 1] *
                     sequence_norm;
              for (m = length[index] - 1;m >= i;m--) {
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

      // calcul de l'entropie des segmentations

#     ifdef MESSAGE
      if (i == 0) {
        for (j = i;j <= length[index] - nb_segment;j++) {
          if (contrast[j] != D_INF) {
            segmentation_entropy -= normalized_contrast[j] * backward[j + 1][1] *
                                    sequence_norm * contrast[j];
          }
        }
      }

      else {
        for (j = MAX(1 , nb_segment + i - length[index]);j < MIN(nb_segment , i + 1);j++) {
          if (j < nb_segment - 1) {
            for (k = i;k <= length[index] + j - nb_segment;k++) {
              if (contrast[k] != D_INF) {
                segmentation_entropy -= forward[i - 1][j - 1] * normalized_contrast[k] * backward[k + 1][j + 1] *
                                        sequence_norm * contrast[k];
              }
            }
          }

          else {
            if (contrast[length[index] - 1] != D_INF) {
              segmentation_entropy -= forward[i - 1][j - 1] * normalized_contrast[length[index] - 1] *
                                      sequence_norm * contrast[length[index] - 1];
            }
          }
        }
      }
#     endif

    }

//    segmentation_entropy = forward_predicted_entropy[length[index] - 1][nb_segment - 1];
//    segmentation_entropy = backward_predicted_entropy[0][0];

#   ifdef DEBUG
    cout << "\n";
//    for (i = length[index] - 1;i >= 0;i--) {
    for (i = 0;i < length[index];i++) {
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
//    for (i = length[index] - 1;i >= 0;i--) {
    for (i = 0;i < length[index];i++) {
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
      for (j = 0;j < length[index];j++) {
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
      for (j = 0;j < length[index] - 1;j++) {
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
      for (j = 0;j < length[index];j++) {
        sum += change_point[i][j];
      }
      if ((sum < i + 1 - DOUBLE_ERROR) || (sum > i + 1 + DOUBLE_ERROR)) {
        cout << "\nERROR: " << sum << " | " << i + 1 << endl;
      }
    }
#   endif

    // extraction des entropies partielles pour les differents nombres de segments possibles

    for (i = 1;i < nb_segment;i++) {
      for (j = 0;j < length[index];j++) {
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

    // calcul de l'entropie des ruptures ordonnees et de l'entropie marginale

    for (i = 0;i < nb_segment - 1;i++) {
      entropy_smoothed[i] = 0.;
    }
    entropy_smoothed[nb_segment - 1] = 1.;

    first_order_entropy = 0.;
    marginal_entropy = 0.;

    for (i = length[index] - 2;i >= 0;i--) {
      sequence_norm = expl(forward_norm[i] + backward_norm[i + 1] - rlikelihood);

/*      for (j = MIN(nb_segment - 1 , i + 1) + 1;j < nb_segment;j++) {
        entropy_smoothed[j] = 0.;
      } */

//      for (j = 0;j < nb_segment;j++) {
      for (j = MAX(0 , nb_segment + i - length[index]);j <= MIN(nb_segment - 1 , i + 1);j++) {
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

    // calcul du profil d'entropies des ruptures et de l'entropie des ruptures

    change_point_entropy[nb_segment - 1][0] = 0.;
    change_point_entropy_sum = 0.;
    for (i = 1;i < length[index];i++) {
      if ((change_point[nb_segment - 1][i] > 0.) && (change_point[nb_segment - 1][i] < 1.)) {

        change_point_entropy[nb_segment - 1][i] = -change_point[nb_segment - 1][i] * logl(change_point[nb_segment - 1][i]) -
                                                  (1 - change_point[nb_segment - 1][i]) * logl(1 - change_point[nb_segment - 1][i]);
         change_point_entropy_sum += change_point_entropy[nb_segment - 1][i];
      }
      else {
        change_point_entropy[nb_segment - 1][i] = 0.;
      }
    }

    // calcul des profils d'entropies des ruptures pour les different nombre de segments possibles

    for (i = 1;i < nb_segment - 1;i++) {
      change_point_entropy[i][0] = 0.;
      for (j = 1;j < length[index];j++) {
        if ((change_point[i][j] > 0.) && (change_point[i][j] < 1.)) {
          change_point_entropy[i][j] = -change_point[i][j] * logl(change_point[i][j]) -
                                       (1 - change_point[i][j]) * logl(1 - change_point[i][j]);
        }
        else {
          change_point_entropy[i][j] = 0.;
        }
      }
    }

    // recurrence "forward" supplementaire pour le calcul des entropies partielles

    for (i = 1;i < nb_segment;i++) {
      for (j = 0;j < length[index];j++) {
        for (k = 0;k < nb_segment;k++) {
          state_entropy[i][j][k] = 0.;
        }
      }
    }

    for (i = 0;i < length[index];i++) {

      // calcul des log-vraisemblances des segments

      for (j = 0;j <= i;j++) {
        contrast[j] = 0.;
      }

      for (j = 1;j < nb_variable;j++) {
        if (model_type[j - 1] == MULTINOMIAL_CHANGE) {
          for (k = 0;k < marginal_distribution[j]->nb_value;k++) {
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

/*            frequency[*pisequence--]++;
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

//      for (j = MAX(0 , nb_segment + i - length[index]);j < MIN((i < length[index] - 1 ? nb_segment - 1 : nb_segment) , i + 1);j++) {
      for (j = 0;j < MIN((i < length[index] - 1 ? nb_segment - 1 : nb_segment) , i + 1);j++) {
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
//        for (j = MAX(0 , nb_segment + i - length[index]);j < MIN((i < length[index] - 1 ? nb_segment - 1 : nb_segment) , i + 1);j++) {
        for (j = 0;j < MIN((i < length[index] - 1 ? nb_segment - 1 : nb_segment) , i + 1);j++) {
          forward[i][j] /= norm[i];
        }

        norm[i] = logl(norm[i]);
      }

      // calcul des entropies partielles

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

      if (i < length[index] - 1) {
        for (j = 1;j < nb_segment;j++) {
          sequence_norm = expl(forward_norm[i] + backward_norm[i + 1] - likelihood[j]);

          for (k = MAX(nb_segment - 1 - j , nb_segment + i - length[index]);k <= MIN(nb_segment - 1 , i+ nb_segment - 1 - j);k++) {
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
      for (j = 0;j < length[index];j++) {
        cout << j << " |";
        for (k = 0;k < nb_segment;k++) {
          cout << " " << state_entropy[i][j][k];
        }
        cout << endl;
      }
    }
#   endif

    // extraction des entropies partielles pour les differents nombres de segments possibles

    for (i = 1;i < nb_segment;i++) {
      for (j = 0;j < length[index] - 1;j++) {
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
           << forward_predicted_entropy[length[index] - 1][i] << ", "
           << backward_predicted_entropy[0][nb_segment - 1 - i] << ", "
           << forward_partial_entropy[i][length[index] - 1] << ", "
           << backward_partial_entropy[i][1] << endl;
    }
    cout << nb_segment << " " << SEQ_label[SEQL_SEGMENTS] << ": "
         << forward_predicted_entropy[length[index] - 1][nb_segment - 1] << ", "
         << backward_predicted_entropy[0][0] << ", "
         << forward_partial_entropy[nb_segment - 1][length[index] - 1] << ", "
         << backward_partial_entropy[nb_segment - 1][1]
         << " | " << segmentation_entropy << endl;
#   endif

#   ifdef DEBUG
    for (i = 1;i < nb_segment;i++) {
      cout << "\n";
      for (j = 0;j < length[index];j++) {
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
      for (j = 1;j < length[index];j++) {
        cout << backward_partial_entropy[i][j] << " ";
      }
      cout << endl;
    }
#   endif

    for (i = 1;i < nb_segment;i++) {
      for (j = length[index] - 1;j >= 1;j--) {
        forward_partial_entropy[i][j] -= forward_partial_entropy[i][j - 1];
        if (forward_partial_entropy[i][j] < 0.) {
          forward_partial_entropy[i][j] = 0.;
        }
      }
    }

    for (i = 1;i < nb_segment;i++) {
      backward_partial_entropy[i][0] = 0.;
      for (j = 1;j < length[index] - 1;j++) {
        backward_partial_entropy[i][j] -= backward_partial_entropy[i][j + 1];
        if (backward_partial_entropy[i][j] < 0.) {
          backward_partial_entropy[i][j] = 0.;
        }
      }
    }

#   ifdef MESSAGE
    for (i = 1;i < nb_segment;i++) {
      for (j = 0;j < length[index];j++) {
        if (forward_partial_entropy[i][j] > change_point_entropy[i][j]) {
          cout << "\n" << SEQ_label[SEQL_BEGIN_CONDITIONAL_ENTROPY] << " ERROR: "
               << forward_partial_entropy[i][j] << " " <<  change_point_entropy[i][j]
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
                            0 , change_point , forward_partial_entropy ,
                            backward_partial_entropy , change_point_entropy);

        *os << "\n" << SEQ_label[SEQL_POSSIBLE_SEGMENTATION_LIKELIHOOD] << ": " << rlikelihood << endl;
        *os << "\n" << SEQ_label[SEQL_SEGMENTATION_ENTROPY] << ": " << segmentation_entropy
            << "\n" << SEQ_label[SEQL_FIRST_ORDER_ENTROPY] << ": " << first_order_entropy
            << "\n" << SEQ_label[SEQL_CHANGE_POINT_ENTROPY] << ": " << change_point_entropy_sum
//            << " (" << change_point_entropy_sum / nb_segment << ")";
            << "\n" << SEQ_label[SEQL_MARGINAL_ENTROPY] << ": " << marginal_entropy << endl;
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
                                  0 , change_point , forward_partial_entropy ,
                                  backward_partial_entropy , change_point_entropy);

        *os << "\n" << SEQ_label[SEQL_POSSIBLE_SEGMENTATION_LIKELIHOOD] << "\t" << rlikelihood << endl;
        *os << "\n" << SEQ_label[SEQL_SEGMENTATION_ENTROPY] << "\t" << segmentation_entropy
            << "\n" << SEQ_label[SEQL_FIRST_ORDER_ENTROPY] << "\t" << first_order_entropy
            << "\n" << SEQ_label[SEQL_CHANGE_POINT_ENTROPY] << "\t" << change_point_entropy_sum
//            << "\t" << change_point_entropy_sum / nb_segment;
            << "\n" << SEQ_label[SEQL_MARGINAL_ENTROPY] << "\t" << marginal_entropy << endl;
        break;
      }

      case 'g' : {
        profile_plot_print(*os , index , nb_segment , backward_output , 0 , change_point ,
                           forward_partial_entropy , backward_partial_entropy ,
                           change_point_entropy);
        break;
      }

      case 'p' : {
        MultiPlotSet &plot = *plot_set;

        i = 1;
        for (j = 1;j < nb_variable;j++) {
          if ((model_type[j - 1] == POISSON_CHANGE) || (model_type[j - 1] == GAUSSIAN_CHANGE) ||
              (model_type[j - 1] == MEAN_CHANGE) || (model_type[j - 1] == VARIANCE_CHANGE) ||
              (model_type[j - 1] == MEAN_VARIANCE_CHANGE)) {
            i++;
          }
        }

        profile_plotable_write(plot[i] , index , nb_segment , backward_output);
        i++;
        change_point_profile_plotable_write(plot[i] , index , nb_segment , change_point);
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

# ifdef DEBUG
  Distribution *prior;


  // calcul des lois a priori des longueurs de segments sous une hypothese de loi a priori
  // uniforme sur les segmentations possibles

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

  delete prior;
# endif

  delete [] frequency;

  for (i = 1;i < nb_variable;i++) {
    delete [] factorial[i];
  }
  delete [] factorial;

  delete [] mean;
  delete [] residual;

  delete [] contrast;
  delete [] normalized_contrast;

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
//  delete [] entropy_norm;

  delete [] likelihood;

  for (i = 0;i < length[index];i++) {
    delete [] backward[i];
  }
  delete [] backward;

  for (i = 0;i < length[index];i++) {
    delete [] backward_predicted_entropy[i];
  }
  delete [] backward_predicted_entropy;

  delete [] backward_norm;

  for (i = 1;i < nb_segment;i++) {
    for (j = 0;j < length[index];j++) {
      delete [] smoothed[i][j];
    }
    delete [] smoothed[i];
  }
  delete [] smoothed;

  for (i = 0;i < length[index];i++) {
    delete [] backward_output[i];
  }
  delete [] backward_output;

  for (i = 1;i < nb_segment;i++) {
    delete [] change_point[i];
  }
  delete [] change_point;

  delete [] entropy_smoothed;

  for (i = 1;i < nb_segment;i++) {
    for (j = 0;j < length[index];j++) {
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
  double sum , factorial_sum , diff , likelihood , segmentation_likelihood , *sequence_mean ,
         *prsequence , *backward , *cumul_backward , **factorial , **srsequence ,
         **mean , **variance;
  long double sum_square , segment_norm , sequence_norm , *residual , *contrast ,
              *norm , **forward;


  max_nb_value = 0;
  factorial = new double*[nb_variable];

  for (i = 1;i < nb_variable;i++) {
    if ((model_type[i - 1] == MULTINOMIAL_CHANGE) &&
        (marginal_distribution[i]->nb_value > max_nb_value)) {
      max_nb_value = marginal_distribution[i]->nb_value;
    }

    if (model_type[i - 1] == POISSON_CHANGE) {
      factorial[i] = new double[length[index]];
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

  sequence_mean = new double[nb_variable];
  residual = new long double[length[index]];

  contrast = new long double[length[index]];

  forward = new long double*[length[index]];
  for (i = 0;i < length[index];i++) {
    forward[i] = new long double[nb_segment];
  }

  norm = new long double[length[index]];

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
      mean[i] = NULL;
      variance[i] = NULL;
    }
  }

  sisequence = new int*[nb_variable];
  srsequence = new double*[nb_variable];

# ifdef DEBUG
  double **segment_probability;

  segment_probability = new double*[length[index]];
  for (i = 0;i < length[index];i++) {
    segment_probability[i] = new double[nb_segment];
    for (j = 0;j < nb_segment;j++) {
      segment_probability[i][j] = 0.;
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
        for (k = 0;k < marginal_distribution[j]->nb_value;k++) {
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
              contrast[k] -= ((double)(i - k + 1) / 2.) * (logl(residual[k] /
                              (i - k + 1)) + log(2 * M_PI) + 1);
/*              contrast[k] -= ((double)(i - k + 1) / 2.) * (logl(residual[k] /
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
          contrast[j] = -((double)((nb_variable - 1) * (i - j + 1)) / 2.) * (logl(contrast[j] /
                           ((nb_variable - 1) * (i - j + 1))) + log(2 * M_PI) + 1);
/*          contrast[j] = -((double)((nb_variable - 1) * (i - j + 1)) / 2.) * (logl(contrast[j] /
                           ((nb_variable - 1) * (i - j))) + log(2 * M_PI)) +
                         (double)((nb_variable - 1) * (i - j)) / 2.; */
        }
        else {
          contrast[j] = D_INF;
        }
      }
    }

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

  sequence_norm = segment_norm + norm[length[index] - 1];

  if (forward[length[index] - 1][nb_segment - 1] > 0.) {
    likelihood = logl(forward[length[index] - 1][nb_segment - 1]) + sequence_norm;
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
            for (n = 0;n < marginal_distribution[m]->nb_value;n++) {
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
                for (r = 0;r < marginal_distribution[m]->nb_value;r++) {
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
                  contrast[n] -= ((double)(j - n + 1) / 2.) * (logl(residual[n] /
                                   (j - n + 1)) + log(2 * M_PI) + 1);
/*                  contrast[n] -= ((double)(j - n + 1) / 2.) * (logl(residual[n] /
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
              contrast[m] = -((double)((nb_variable - 1) * (j - m + 1)) / 2.) * (logl(contrast[m] /
                               ((nb_variable - 1) * (j - m + 1))) + log(2 * M_PI) + 1);
/*              contrast[m] = -((double)((nb_variable - 1) * (j - m + 1)) / 2.) * (logl(contrast[m] /
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

        segmentation_likelihood += logl(contrast[j - segment_length + 1]);

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

      // approximation des probabilites lissees

      psegment = int_sequence[index][0];
      for (j = 0;j < length[index];j++) {
        segment_probability[j][*psegment++]++;
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
          segment_probability[i][j] /= nb_segmentation;
        }
      }

      psegment = int_sequence[index][0];
      for (i = 0;i < length[index];i++) {
        *psegment++ = I_DEFAULT;
      }

      os << "\n" << SEQ_label[SEQL_POSTERIOR_SEGMENT_PROBABILITY] << "\n\n";

      profile_ascii_print(os , index , nb_segment , segment_probability ,
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
    delete [] segment_probability[i];
  }
  delete [] segment_probability;
# endif

  return likelihood;
}


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
  int max_nb_value , brank , previous_rank , nb_cell , *frequency , *pisequence , *rank ,
      *psegment , **sisequence , ***optimal_length , ***optimal_rank;
  double sum , factorial_sum , diff , buff , segmentation_likelihood , *nb_segmentation ,
         *sequence_mean , *prsequence , **factorial , **nb_segmentation_forward ,
         **srsequence , **mean , **variance , ***forward;
  long double sum_square , *residual , *contrast , likelihood_cumul;

# ifdef MESSAGE
  double sum2;
# endif


  // initialisations

  max_nb_value = 0;
  factorial = new double*[nb_variable];

  for (i = 1;i < nb_variable;i++) {
    if ((model_type[i - 1] == MULTINOMIAL_CHANGE) &&
        (marginal_distribution[i]->nb_value > max_nb_value)) {
      max_nb_value = marginal_distribution[i]->nb_value;
    }

    if (model_type[i - 1] == POISSON_CHANGE) {
      factorial[i] = new double[length[index]];
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
  for (i = 1;i < nb_variable;i++) {
    if ((model_type[i - 1] == POISSON_CHANGE) || (model_type[i - 1] == GAUSSIAN_CHANGE) ||
        (model_type[i - 1] == MEAN_CHANGE) || (model_type[i - 1] == VARIANCE_CHANGE) ||
        (model_type[i - 1] == MEAN_VARIANCE_CHANGE)) {
      mean[i] = new double[nb_segment];
      variance[i] = new double[nb_segment];
    }
    else {
      mean[i] = NULL;
      variance[i] = NULL;
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

# ifdef DEBUG
  for (i = 0;i < nb_segment;i++) {
    nb_segmentation[i] = 1;
  }
# endif

  // recurrence "forward"

  for (i = 0;i < length[index];i++) {

    // calcul des log-vraisemblances des segments

    for (j = 0;j <= i;j++) {
      contrast[j] = 0.;
    }

    for (j = 1;j < nb_variable;j++) {
      if (model_type[j - 1] == MULTINOMIAL_CHANGE) {
        for (k = 0;k < marginal_distribution[j]->nb_value;k++) {
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
          for (m = 0;m < marginal_distribution[j]->nb_value;m++) {
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
              contrast[k] -= ((double)(i - k + 1) / 2.) * (logl(residual[k] /
                               (i - k + 1)) + log(2 * M_PI) + 1);
/*              contrast[k] -= ((double)(i - k + 1) / 2.) * (logl(residual[k] /
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
          contrast[j] = -((double)((nb_variable - 1) * (i - j + 1)) / 2.) * (logl(contrast[j] /
                           ((nb_variable - 1) * (i - j + 1))) + log(2 * M_PI) + 1);
/*          contrast[j] = -((double)((nb_variable - 1) * (i - j + 1)) / 2.) * (logl(contrast[j] /
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
# endif

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
      likelihood_cumul += exp(forward[length[index] - 1][nb_segment - 1][i]);
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
        os << exp(forward[length[index] - 1][nb_segment - 1][i] - likelihood)
           << "  " << likelihood_cumul / exp(likelihood);
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
           << "\t" << likelihood_cumul / exp(likelihood);
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

  delete [] sequence_mean;
  delete [] residual;

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
  int max_nb_value , *frequency , *pisequence , *psegment , **optimal_length;
  double sum , factorial_sum , diff , buff , segmentation_likelihood , backward_max ,
         *sequence_mean , *prsequence , **factorial , **forward , **backward ,
         **backward_output , **mean;
  long double sum_square , *residual , *contrast;


  max_nb_value = 0;
  factorial = new double*[nb_variable];

  for (i = 1;i < nb_variable;i++) {
    if ((model_type[i - 1] == MULTINOMIAL_CHANGE) &&
        (marginal_distribution[i]->nb_value > max_nb_value)) {
      max_nb_value = marginal_distribution[i]->nb_value;
    }

    if (model_type[i - 1] == POISSON_CHANGE) {
      factorial[i] = new double[length[index]];
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

  mean = new double*[nb_variable];
  for (i = 1;i < nb_variable;i++) {
    if ((model_type[i - 1] == POISSON_CHANGE) || (model_type[i - 1] == GAUSSIAN_CHANGE) ||
        (model_type[i - 1] == MEAN_CHANGE) || (model_type[i - 1] == VARIANCE_CHANGE) ||
        (model_type[i - 1] == MEAN_VARIANCE_CHANGE)) {
      mean[i] = new double[length[index]];
    }
    else {
      mean[i] = NULL;
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
        for (k = 0;k < marginal_distribution[j]->nb_value;k++) {
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
              contrast[k] -= ((double)(i - k + 1) / 2.) * (logl(residual[k] /
                               (i - k + 1)) + log(2 * M_PI) + 1);
/*              contrast[k] -= ((double)(i - k + 1) / 2.) * (logl(residual[k] /
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
          contrast[j] = -((double)((nb_variable - 1) * (i - j + 1)) / 2.) * (logl(contrast[j] /
                           ((nb_variable - 1) * (i - j + 1))) + log(2 * M_PI) + 1);
/*          contrast[j] = -((double)((nb_variable - 1) * (i - j + 1)) / 2.) * (logl(contrast[j] /
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
          for (k = 0;k < marginal_distribution[j]->nb_value;k++) {
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
                contrast[k] -= ((double)(k - i + 1) / 2.) * (logl(residual[k] /
                                 (k - i + 1)) + log(2 * M_PI) + 1);
/*                contrast[k] -= ((double)(k - i + 1) / 2.) * (logl(residual[k] /
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
            contrast[j] = -((double)((nb_variable - 1) * (j - i + 1)) / 2.) * (logl(contrast[j] /
                             ((nb_variable - 1) * (j - i + 1))) + log(2 * M_PI) + 1);
/*            contrast[j] = -((double)((nb_variable - 1) * (j - i + 1)) / 2.) * (logl(contrast[j] /
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
                          mean);

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
                                mean);

      *os << "\n" << SEQ_label[SEQL_SEGMENTATION_LIKELIHOOD] << "\t" << segmentation_likelihood;
      if (likelihood != D_INF) {
        *os << "\t" << exp(segmentation_likelihood - likelihood);
      }
      *os << endl;
      break;
    }

    case 'g' : {
      profile_plot_print(*os , index , nb_segment , backward_output , mean);
      break;
    }

    case 'p' : {
      MultiPlotSet &plot = *plot_set;

      i = 0;
      for (j = 1;j < nb_variable;j++) {
        if ((model_type[j - 1] == POISSON_CHANGE) || (model_type[j - 1] == GAUSSIAN_CHANGE) ||
            (model_type[j - 1] == MEAN_CHANGE) || (model_type[j - 1] == VARIANCE_CHANGE) ||
            (model_type[j - 1] == MEAN_VARIANCE_CHANGE)) {
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
              plot[i][1].add_point(index_parameter[index][k] , mean[j][k]);
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
              plot[i][1].add_point(k , mean[j][k]);
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
  Sequences *seq;


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
        rank[i] = seq->marginal_distribution[i]->rank_computation();
      }
      else {
        rank[i] = NULL;
      }
    }

    for (i = 0;i < seq->nb_sequence;i++) {
      if ((index == I_DEFAULT) || (index == i)) {
        if (model_type[0] != MEAN_CHANGE) {
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
  Sequences *seq;
  ostringstream data_file_name[2];
  ofstream *out_data_file;


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
          rank[i] = seq->marginal_distribution[i]->rank_computation();
        }
        else {
          rank[i] = NULL;
        }
      }

      if (model_type[0] != MEAN_CHANGE) {
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
              if ((model_type[k - 1] == POISSON_CHANGE) || (model_type[k - 1] == GAUSSIAN_CHANGE) ||
                  (model_type[k - 1] == MEAN_CHANGE) || (model_type[k - 1] == VARIANCE_CHANGE) ||
                  (model_type[k - 1] == MEAN_VARIANCE_CHANGE)) {
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
                         << "\" with linespoints" << ",\\" << endl;
                out_file << "\"" << label((data_file_name[0].str()).c_str()) << "\" using " << 1 << " : " << j++
                         << " title \"" << SEQ_label[SEQL_SEGMENT_MEAN] << "\" with lines" << endl;

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
              if ((model_type[k - 1] == POISSON_CHANGE) || (model_type[k - 1] == GAUSSIAN_CHANGE) ||
                  (model_type[k - 1] == MEAN_CHANGE) || (model_type[k - 1] == VARIANCE_CHANGE) ||
                  (model_type[k - 1] == MEAN_VARIANCE_CHANGE)) {
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
                         << " title \"" << SEQ_label[SEQL_SEGMENT_MEAN] << "\" with lines" << endl;

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
  Sequences *seq;
  ostringstream title , legend;
  MultiPlotSet *plot_set;


  plot_set = NULL;
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
    seq = new Sequences(*this , 'a');

    // calcul du nombre de vues

    nb_plot_set = 1;
    for (i = 1;i < seq->nb_variable;i++) {
      if ((model_type[i - 1] == POISSON_CHANGE) || (model_type[i - 1] == GAUSSIAN_CHANGE) ||
          (model_type[i - 1] == MEAN_CHANGE) || (model_type[i - 1] == VARIANCE_CHANGE) ||
          (model_type[i - 1] == MEAN_VARIANCE_CHANGE)) {
        nb_plot_set++;
      }
    }
    if (model_type[0] != MEAN_CHANGE) {
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

    if (model_type[0] != MEAN_CHANGE) {
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

      // vues : sequence et fonction en escalier

      for (j = 1;j < seq->nb_variable;j++) {
        if ((model_type[j - 1] == POISSON_CHANGE) || (model_type[j - 1] == GAUSSIAN_CHANGE) ||
            (model_type[j - 1] == MEAN_CHANGE) || (model_type[j - 1] == VARIANCE_CHANGE) ||
            (model_type[j - 1] == MEAN_VARIANCE_CHANGE)) {
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

          plot[i][1].legend = SEQ_label[SEQL_SEGMENT_MEAN];
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
