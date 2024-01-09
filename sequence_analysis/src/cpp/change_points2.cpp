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
 *       $Id: change_points2.cpp 18046 2015-04-23 09:34:21Z guedon $
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


extern int column_width(int nb_value , const long double *value);


#if defined (SYSTEM_IS__CYGWIN)
#define expl exp
#endif



/*--------------------------------------------------------------*/
/**
 *  Ecriture des profils de segments/d'etats, des profils de ruptures et
 *  des profils d'entropies pour une sequence.
 *
 *  arguments : stream, indice de la sequence, nombre de segments/d'etats,
 *              pointeur sur les profils de segments/d'etats, label,
 *              pointeur sur les fonctions lineares par morceau,
 *              les profils de ruptures et les profils d'entropies.
 */
/*--------------------------------------------------------------*/

ostream& Sequences::profile_ascii_print(ostream &os , int index , int nb_segment ,
                                        double **profiles , const char *label ,
                                        double **piecewise_function , long double **change_point ,
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
      width[i] = stat_tool::column_width((int)min_value[i] , (int)max_value[i]);
    }
    else {
      width[i] = stat_tool::column_width(length[index] , real_sequence[index][i]);
    }
    if (i > start) {
      width[i] += ASCII_SPACE;
    }
  }

  if (index_parameter) {
    width[nb_variable] = stat_tool::column_width(index_parameter_distribution->nb_value - 1) + ASCII_SPACE;
  }
  else {
    width[nb_variable] = stat_tool::column_width(max_length) + ASCII_SPACE;
  }

  width[nb_variable + 1] = 0;
  for (i = 0;i < length[index];i++) {
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
        width[nb_variable + 2 + i] = stat_tool::column_width(length[index] , piecewise_function[i]) + ASCII_SPACE;
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
    if ((piecewise_function) && (piecewise_function[i])) {
      os << SEQ_label[SEQL_PIECEWISE_LINEAR_FUNCTION] << " " << i << " | ";
    }
    os << STAT_label[STATL_VARIABLE] << " " << i << " | ";
  }

  switch (index_parameter_type) {
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

  for (i = 0;i < length[index];i++) {
    os.setf(ios::right , ios::adjustfield);
    if (!change_point) {
      os << setw(width[0]) << int_sequence[index][0][i];
    }

    for (j = 1;j < nb_variable;j++) {
      if ((piecewise_function) && (piecewise_function[j])) {
        os << setw(width[nb_variable + 2 + j]) << piecewise_function[j][i];
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
    switch (index_parameter_type) {
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
 *              pointeur sur les fonctions lineares par morceaux,
 *              les profils de ruptures et les profils d'entropies.
 *
 *--------------------------------------------------------------*/

ostream& Sequences::profile_spreadsheet_print(ostream &os , int index , int nb_segment ,
                                              double **profiles , const char *label ,
                                              double **piecewise_function , long double **change_point ,
                                              long double **begin_conditonal_entropy ,
                                              long double **end_conditional_entropy ,
                                              long double **change_point_entropy) const

{
  register int i , j , k;


  if (!change_point) {
    os << SEQ_label[SEQL_OPTIMAL] << " " << label << "\t";
  }
  for (i = 1;i < nb_variable;i++) {
    if ((piecewise_function) && (piecewise_function[i])) {
      os << SEQ_label[SEQL_PIECEWISE_LINEAR_FUNCTION] << " " << i << "\t";
    }
    os << STAT_label[STATL_VARIABLE] << " " << i << "\t";
  }

  switch (index_parameter_type) {
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

  for (i = 0;i < length[index];i++) {
    if (!change_point) {
      os << int_sequence[index][0][i] << "\t";
    }
    for (j = 1;j < nb_variable;j++) {
      if ((piecewise_function) && (piecewise_function[j])) {
        os << piecewise_function[j][i] << "\t";
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
    switch (index_parameter_type) {
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
 *              pointeur sur les profils de segments/d'etats,
 *              pointeur sur les fonctions lineares par morceaux,
 *              les profils de rupture et les profils d'entropies.
 *
 *--------------------------------------------------------------*/

ostream& Sequences::profile_plot_print(ostream &os , int index , int nb_segment ,
                                       double **profiles , double **piecewise_function ,
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
      if ((piecewise_function) && (piecewise_function[j])) {
        if (type[j] != REAL_VALUE) {
          os << int_sequence[index][j][i];
        }
        else {
          os << real_sequence[index][j][i];
        }
        os << " " << piecewise_function[j][i] << " ";
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
  int max_nb_value , *frequency , *seq_index_parameter = NULL;
  double sum , factorial_sum , proba , diff , index_parameter_sum , index_parameter_diff , buff ,
         rlikelihood , *mean , *likelihood , **hyperparam , **factorial , **backward_output , ***smoothed;
  long double square_sum , index_parameter_square_sum , mix_square_sum , segment_norm , sequence_norm ,
              lbuff , lsum , prior_contrast , segmentation_entropy , first_order_entropy ,
              change_point_entropy_sum , marginal_entropy , *residual , *contrast , *normalized_contrast ,
              *norm , *forward_norm , *backward_norm , *entropy_smoothed , *segment_predicted ,
              **forward , **backward , **change_point , **forward_predicted_entropy ,
              **backward_predicted_entropy , **forward_partial_entropy ,
              **backward_partial_entropy , **change_point_entropy , ***state_entropy;

# ifdef DEBUG
  long double *entropy_norm;
# endif


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

          for (k = MAX(1 , j + 1 + i - length[index]);k <= MIN(j , i);k++) {
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

/*              frequency[int_sequence[index][j][k]]++;
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
//                if ((contrast[k] != D_INF) && (residual[k] > 0.)) {
                if ((contrast[k] != D_INF) && (residual[k] > sqrt((double)(i - k + 1)) * ROUNDOFF_ERROR)) {
                  contrast[k] -= ((double)(i - k + 1) / 2.) * (logl(residual[k] /
                                   (i - k + 1)) + log(2 * M_PI) + 1);
/*                  contrast[k] -= ((double)(i - k + 1) / 2.) * (logl(residual[k] /
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
//            if (contrast[j] > 0.) {
            if (contrast[j] > sqrt((double)((nb_variable - 1) * (i - j + 1))) * ROUNDOFF_ERROR) {
              contrast[j] = -((double)((nb_variable - 1) * (i - j + 1)) / 2.) * (logl(contrast[j] /
                               ((nb_variable - 1) * (i - j + 1))) + log(2 * M_PI) + 1);
/*              contrast[j] = -((double)((nb_variable - 1) * (i - j + 1)) / 2.) * (logl(contrast[j] /
                               ((nb_variable - 1) * (i - j))) + log(2 * M_PI)) +
                             (double)((nb_variable - 1) * (i - j)) / 2.; */
            }
            else {
              contrast[j] = D_INF;
            }
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

          for (k = MAX(nb_segment - 1 - j , nb_segment + i - length[index]);k <= MIN(nb_segment - 1 , i + nb_segment - 1 - j);k++) {
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
                            NULL , change_point , forward_partial_entropy ,
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
                                  NULL , change_point , forward_partial_entropy ,
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
        profile_plot_print(*os , index , nb_segment , backward_output ,
                           NULL , change_point , forward_partial_entropy ,
                           backward_partial_entropy , change_point_entropy);
        break;
      }

      case 'p' : {
        MultiPlotSet &plot = *plot_set;

        i = 1;
        for (j = 1;j < nb_variable;j++) {
          if ((model_type[j - 1] == POISSON_CHANGE) || (model_type[0] == MULTIVARIATE_POISSON_CHANGE) ||
              (model_type[j - 1] == GEOMETRIC_0_CHANGE) || (model_type[j - 1] == GEOMETRIC_1_CHANGE) ||
              (model_type[0] == MULTIVARIATE_GEOMETRIC_0_CHANGE) || (model_type[j - 1] == GAUSSIAN_CHANGE) ||
              (model_type[j - 1] == VARIANCE_CHANGE) || (model_type[0] == MEAN_CHANGE) ||
              (model_type[0] == MEAN_VARIANCE_CHANGE)|| (model_type[j - 1] == BAYESIAN_POISSON_CHANGE) ||
              (model_type[j - 1] == BAYESIAN_GAUSSIAN_CHANGE)) {
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
  int max_nb_value , segment_length , *frequency , *seq_index_parameter = NULL , *psegment;
  double sum , factorial_sum , proba , diff , index_parameter_sum , index_parameter_diff ,
         likelihood , segmentation_likelihood , response_mean , index_parameter_mean ,
         index_parameter_variance , covariance , residual_mean , *sequence_mean , *backward ,
         *cumul_backward , **hyperparam , **factorial , **mean , **variance , **intercept , **slope;
  long double square_sum , index_parameter_square_sum , mix_square_sum , segment_norm , sequence_norm ,
              prior_contrast , *residual , *contrast , *norm , **forward;


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

  forward = new long double*[length[index]];
  for (i = 0;i < length[index];i++) {
    forward[i] = new long double[nb_segment];
  }

  norm = new long double[length[index]];

  backward = new double[length[index]];
  cumul_backward = new double[length[index]];

  mean = new double*[nb_variable];
  variance = new double*[nb_variable];
  intercept = new double*[nb_variable];
  slope = new double*[nb_variable];

  for (i = 1;i < nb_variable;i++) {
    if ((model_type[i - 1] == POISSON_CHANGE) || (model_type[0] == MULTIVARIATE_POISSON_CHANGE) ||
        (model_type[i - 1] == GAUSSIAN_CHANGE) || (model_type[i - 1] == VARIANCE_CHANGE) ||
        (model_type[0] == MEAN_VARIANCE_CHANGE) || (model_type[i - 1] == BAYESIAN_POISSON_CHANGE) ||
        (model_type[i - 1] == BAYESIAN_GAUSSIAN_CHANGE)) {
      mean[i] = new double[nb_segment];
      variance[i] = new double[nb_segment];
    }
    else if ((model_type[i - 1] == GEOMETRIC_0_CHANGE) || (model_type[i - 1] == GEOMETRIC_1_CHANGE) ||
             (model_type[0] == MULTIVARIATE_GEOMETRIC_0_CHANGE)) {
      mean[i] = new double[nb_segment];
      variance[i] = NULL;
    }
    else if (model_type[i - 1] == LINEAR_MODEL_CHANGE) {
      mean[i] = NULL;
      variance[i] = new double[nb_segment];
    }
    else {
      mean[i] = NULL;
      variance[i] = NULL;
    }

    if (model_type[i - 1] == LINEAR_MODEL_CHANGE) {
      intercept[i] = new double[nb_segment];
      slope[i] = new double[nb_segment];
    }
    else {
      intercept[i] = NULL;
      slope[i] = NULL;
    }
  }

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
      segmentation_likelihood = sequence_norm;

      for (k = nb_segment - 1;k >= 0;k--) {

        // calcul des log-vraisemblances des segments

        if (model_type[0] == MULTIVARIATE_CATEGORICAL_CHANGE) {
          for (m = 0;m < max_nb_value;m++) {
            frequency[m] = 0;
          }

          for (m = j;m >= k;m--) {
            for (n = 1;n < nb_variable;n++) {
              frequency[int_sequence[index][n][m]]++;
            }

            contrast[m] = 0.;
            for (n = 0;n < max_nb_value;n++) {
              if (frequency[n] > 0) {
                contrast[m] += frequency[n] * log((double)frequency[n] / (double)((nb_variable - 1) * (j - m + 1)));
              }
            }
          }
        }

        else if (model_type[0] == MULTIVARIATE_POISSON_CHANGE) {
          sum = 0.;
          factorial_sum = 0.;
          for (m = j;m >= k;m--) {
            for (n = 1;n < nb_variable;n++) {
              sum += int_sequence[index][n][m];
              factorial_sum += factorial[n][m];
            }
            if (sum > 0.) {
              contrast[m] = sum * (log(sum / ((nb_variable - 1) * (j - m + 1))) - 1) - factorial_sum;
            }
            else {
              contrast[m] = 0.;
            }
          }
        }

        else if (model_type[0] == MULTIVARIATE_GEOMETRIC_0_CHANGE) {
          sum = 0.;
          for (m = j;m >= k;m--) {
            for (n = 1;n < nb_variable;n++) {
              sum += int_sequence[index][n][m];
            }
            if (sum > 0.) {
              proba = (nb_variable - 1) * (j - m + 1) / ((nb_variable - 1) * (j - m + 1) + sum);
              contrast[m] = (nb_variable - 1) * (j - m + 1) * log(proba) + sum * log(1. - proba);
            }
            else {
              contrast[m] = 0.;
            }
          }
        }

        else {
          for (m = k;m <= j;m++) {
            contrast[m] = 0.;
          }

          for (m = 1;m < nb_variable;m++) {
            if (model_type[m - 1] == CATEGORICAL_CHANGE) {
              for (n = 0;n < marginal_distribution[m]->nb_value;n++) {
                frequency[n] = 0;
              }
              sum = 0.;

              frequency[int_sequence[index][m][j]]++;
              for (n = j - 1;n >= k;n--) {
                sum += (j - n) * log((double)(j - n) / (double)(j - n + 1)) +
                       log((double)(frequency[int_sequence[index][m][n]] + 1) / (double)(j - n + 1));
                if (frequency[int_sequence[index][m][n]] > 0) {
                  sum -= frequency[int_sequence[index][m][n]] *
                         log((double)frequency[int_sequence[index][m][n]] / (double)(frequency[int_sequence[index][m][n]] + 1));
                }
                frequency[int_sequence[index][m][n]]++;

                if (contrast[n] != D_INF) {
                  contrast[n] += sum;
                }

/*                frequency[int_sequence[index][m][n]]++;
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
              for (n = j;n >= k;n--) {
                sum += int_sequence[index][m][n];
                factorial_sum += factorial[m][n];
                if ((contrast[n] != D_INF) && (sum > 0.)) {
                  contrast[n] += sum * (log(sum / (j - n + 1)) - 1) - factorial_sum;
                }
              }
            }

            else if (model_type[m - 1] == GEOMETRIC_0_CHANGE) {
              sum = 0.;
              for (n = j;n >= k;n--) {
                sum += int_sequence[index][m][n];
                if ((contrast[n] != D_INF) && (sum > 0.)) {
                  proba = (j - n + 1) / (j - n + 1 + sum);
                  contrast[n] += (j - n + 1) * log(proba) + sum * log(1. - proba);
                }
              }
            }

            else if (model_type[m - 1] == GEOMETRIC_1_CHANGE) {
              sum = 0.;
              for (n = j;n >= k;n--) {
                sum += int_sequence[index][m][n];
                if ((contrast[n] != D_INF) && (sum > j - n + 1)) {
                  proba = (j - n + 1) / sum;
                  contrast[n] += (j - n + 1) * log(proba) + (sum - (j - n + 1)) * log(1. - proba);
                }
              }
            }

            else if (model_type[m - 1] == BAYESIAN_POISSON_CHANGE) {
              prior_contrast = -lgamma(hyperparam[m][0]) + hyperparam[m][0] * log(hyperparam[m][1]);

              sum = 0.;
              factorial_sum = 0.;
              for (n = j;n >= k;n--) {
                sum += int_sequence[index][m][n]++;
                factorial_sum += factorial[m][n];
                if (contrast[n] != D_INF) {
                  contrast[n] += prior_contrast - factorial_sum + lgamma(hyperparam[m][0] + sum) -
                                 (hyperparam[m][0] + sum) * log(hyperparam[m][1] + j - n + 1);
                }
              }
            }

            else if (model_type[m - 1] == BAYESIAN_GAUSSIAN_CHANGE) {
              prior_contrast = log(hyperparam[m][1]) / 2 - lgamma(hyperparam[m][2] / 2) +
                               hyperparam[m][2] * log(hyperparam[m][3] / 2) / 2;

              if (type[m] != REAL_VALUE) {
                square_sum = 0.;
                sum = int_sequence[index][m][j];
                if (contrast[j] != D_INF) {
                  diff = hyperparam[m][0] - sum;
                  contrast[j] += prior_contrast - log(2 * M_PI) / 2 -
                                 log(hyperparam[m][1] + 1) / 2 + lgamma((hyperparam[m][2] + 1) / 2) -
                                 (hyperparam[m][2] + 1) *
                                 log((hyperparam[m][3] + hyperparam[m][1] *
                                      diff * diff / (hyperparam[m][1] + 1)) / 2) / 2;
                }

                for (n = j - 1;n >= k;n--) {
                  diff = int_sequence[index][m][n] - sum / (j - n);
                  square_sum += ((double)(j - n) / (double)(j - n + 1)) * diff * diff;
                  sum += int_sequence[index][m][n];
                  if (contrast[n] != D_INF) {
                    diff = hyperparam[m][0] - sum / (j - n + 1);
                    contrast[n] += prior_contrast - (j - n + 1) * log(2 * M_PI) / 2 -
                                   log(hyperparam[m][1] + j - n + 1) / 2 +
                                   lgamma((hyperparam[m][2] + j - n + 1) / 2) -
                                   (hyperparam[m][2] + j - n + 1) *
                                   logl((hyperparam[m][3] + square_sum + hyperparam[m][1] * (j - n + 1) *
                                         diff * diff / (hyperparam[m][1] + j - n + 1)) / 2) / 2;
                  }
                }
              }

              else {
                square_sum = 0.;
                sum = real_sequence[index][m][j];
                if (contrast[j] != D_INF) {
                  diff = hyperparam[m][0] - sum;
                  contrast[j] += prior_contrast - log(2 * M_PI) / 2 -
                                 log(hyperparam[m][1] + 1) / 2 + lgamma((hyperparam[m][2] + 1) / 2) -
                                 (hyperparam[m][2] + 1) *
                                 log((hyperparam[m][3] + hyperparam[m][1] *
                                      diff * diff / (hyperparam[m][1] + 1)) / 2) / 2;
                }

                for (n = j - 1;n >= k;n--) {
                  diff = real_sequence[index][m][n] - sum / (j - n);
                  square_sum += ((double)(j - n) / (double)(j - n + 1)) * diff * diff;
                  sum += real_sequence[index][m][n];
                  if (contrast[n] != D_INF) {
                    diff = hyperparam[m][0] - sum / (j - n + 1);
                    contrast[n] += prior_contrast - (j - n + 1) * log(2 * M_PI) / 2 -
                                   log(hyperparam[m][1] + j - n + 1) / 2 +
                                   lgamma((hyperparam[m][2] + j - n + 1) / 2) -
                                   (hyperparam[m][2] + j - n + 1) *
                                   logl((hyperparam[m][3] + square_sum + hyperparam[m][1] * (j - n + 1) *
                                         diff * diff / (hyperparam[m][1] + j - n + 1)) / 2) / 2;
                  }
                }
              }
            }

            else {
              if (model_type[m - 1] == VARIANCE_CHANGE) {
                square_sum = 0.;

                if (type[m] != REAL_VALUE) {
                  for (n = j;n >= k;n--) {
                    diff = int_sequence[index][m][n] - sequence_mean[m];
                    square_sum += diff * diff;
                    residual[n] = square_sum;
                  }
                }

                else {
                  for (n = j;n >= k;n--) {
                    diff = real_sequence[index][m][n] - sequence_mean[m];
                    square_sum += diff * diff;
                    residual[n] = square_sum;
                  }
                }
              }

              else if (model_type[m - 1] == ORDINAL_GAUSSIAN_CHANGE) {
                square_sum = 0.;
                sum = rank[m][int_sequence[index][m][j]];
                residual[j] = 0.;

                for (n = j - 1;n >= k;n--) {
                  diff = rank[m][int_sequence[index][m][n]] - sum / (j - n);
                  square_sum += ((double)(j - n) / (double)(j - n + 1)) * diff * diff;
                  sum += rank[m][int_sequence[index][m][n]];
                  residual[n] = square_sum;
                }
              }

              else if (model_type[m - 1] == LINEAR_MODEL_CHANGE) {
                if (type[m] != REAL_VALUE) {
                  square_sum = 0.;
                  index_parameter_square_sum = 0.;
                  mix_square_sum = 0.;
                  sum = int_sequence[index][m][j];
                  index_parameter_sum = seq_index_parameter[j];
                  residual[j] = 0.;

                  for (n = j - 1;n >= k;n--) {
                    diff = int_sequence[index][m][n] - sum / (j - n);
                    square_sum += ((double)(j - n) / (double)(j - n + 1)) * diff * diff;
                    index_parameter_diff = seq_index_parameter[n] - index_parameter_sum / (j - n);
                    index_parameter_square_sum += ((double)(j - n) / (double)(j - n + 1)) *
                                                  index_parameter_diff * index_parameter_diff;
                    mix_square_sum += ((double)(j - n) / (double)(j - n + 1)) * diff * index_parameter_diff;
                    sum += int_sequence[index][m][n];
                    index_parameter_sum += seq_index_parameter[n];

                    if ((n < j - 1) && (index_parameter_square_sum > 0.)) {
                      residual[n] = square_sum - mix_square_sum * mix_square_sum / index_parameter_square_sum;
                    }
                    else {
                      residual[n] = 0.;
                    }
                  }
                }

                else {
                  square_sum = 0.;
                  index_parameter_square_sum = 0.;
                  mix_square_sum = 0.;
                  sum = real_sequence[index][m][j];
                  index_parameter_sum = seq_index_parameter[j];
                  residual[j] = 0.;

                  for (n = j - 1;n >= k;n--) {
                    diff = real_sequence[index][m][n] - sum / (j - n);
                    square_sum += ((double)(j - n) / (double)(j - n + 1)) * diff * diff;
                    index_parameter_diff = seq_index_parameter[n] - index_parameter_sum / (j - n);
                    index_parameter_square_sum += ((double)(j - n) / (double)(j - n + 1)) *
                                                  index_parameter_diff * index_parameter_diff;
                    mix_square_sum += ((double)(j - n) / (double)(j - n + 1)) * diff * index_parameter_diff;
                    sum += real_sequence[index][m][n];
                    index_parameter_sum += seq_index_parameter[n];

                    if ((n < j - 1) && (index_parameter_square_sum > 0.)) {
                      residual[n] = square_sum - mix_square_sum * mix_square_sum / index_parameter_square_sum;
                    }
                    else {
                      residual[n] = 0.;
                    }
                  }
                }
              }

              else {
                if (type[m] != REAL_VALUE) {
                  square_sum = 0.;
                  sum = int_sequence[index][m][j];
                  residual[j] = 0.;

                  for (n = j - 1;n >= k;n--) {
                    diff = int_sequence[index][m][n] - sum / (j - n);
                    square_sum += ((double)(j - n) / (double)(j - n + 1)) * diff * diff;
                    sum += int_sequence[index][m][n];
                    residual[n] = square_sum;
                  }
                }

                else {
                  square_sum = 0.;
                  sum = real_sequence[index][m][j];
                  residual[j] = 0.;

                  for (n = j - 1;n >= k;n--) {
                    diff = real_sequence[index][m][n] - sum / (j - n);
                    square_sum += ((double)(j - n) / (double)(j - n + 1)) * diff * diff;
                    sum += real_sequence[index][m][n];
                    residual[n] = square_sum;
                  }
                }
              }

              if (model_type[0] == MEAN_VARIANCE_CHANGE) {
                for (n = j - 1;n >= k;n--) {
                  contrast[n] += residual[n];
                }
              }

              else {
                for (n = j;n >= k;n--) {
//                  if ((contrast[n] != D_INF) && (residual[n] > 0.)) {
                  if ((contrast[n] != D_INF) && (residual[n] > sqrt((double)(j - n + 1)) * ROUNDOFF_ERROR)) {
                    contrast[n] -= ((double)(j - n + 1) / 2.) * (logl(residual[n] /
                                     (j - n + 1)) + log(2 * M_PI) + 1);
/*                    contrast[n] -= ((double)(j - n + 1) / 2.) * (logl(residual[n] /
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
//              if (contrast[m] > 0.)) {
              if (contrast[m] > sqrt((double)((nb_variable - 1) * (j - m + 1))) * ROUNDOFF_ERROR) {
                contrast[m] = -((double)((nb_variable - 1) * (j - m + 1)) / 2.) * (logl(contrast[m] /
                                 ((nb_variable - 1) * (j - m + 1))) + log(2 * M_PI) + 1);
/*                contrast[m] = -((double)((nb_variable - 1) * (j - m + 1)) / 2.) * (logl(contrast[m] /
                                 ((nb_variable - 1) * (j - m))) + log(2 * M_PI)) +
                               (double)((nb_variable - 1) * (j - m)) / 2.; */
              }
              else {
                contrast[m] = D_INF;
              }
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

        for (m = 1;m < nb_variable;m++) {
          if ((model_type[m - 1] == POISSON_CHANGE) || (model_type[0] == MULTIVARIATE_POISSON_CHANGE) ||
              (model_type[m - 1] == GEOMETRIC_0_CHANGE) || (model_type[m - 1] == GEOMETRIC_1_CHANGE) ||
              (model_type[0] == MULTIVARIATE_GEOMETRIC_0_CHANGE) || (model_type[m - 1] == GAUSSIAN_CHANGE) ||
              (model_type[m - 1] == VARIANCE_CHANGE) || (model_type[0] == MEAN_VARIANCE_CHANGE) ||
              (model_type[m - 1] == BAYESIAN_POISSON_CHANGE) || (model_type[m - 1] == BAYESIAN_GAUSSIAN_CHANGE)) {
            mean[m][k] = 0.;
            if (type[m] != REAL_VALUE) {
              for (n = j;n > j - segment_length;n--) {
                mean[m][k] += int_sequence[index][m][n];
              }
            }
            else {
              for (n = j;n > j - segment_length;n--) {
                mean[m][k] += real_sequence[index][m][n];
              }
            }
            mean[m][k] /= segment_length;

            if ((model_type[m - 1] != GEOMETRIC_0_CHANGE) && (model_type[m - 1] != GEOMETRIC_1_CHANGE) &&
                (model_type[0] != MULTIVARIATE_GEOMETRIC_0_CHANGE)) {
              variance[m][k] = 0.;
              if (segment_length > 1) {
                if (type[m] != REAL_VALUE) {
                  for (n = j;n > j - segment_length;n--) {
                    diff = int_sequence[index][m][n] - mean[m][k];
                    variance[m][k] += diff * diff;
                  }
                }
                else {
                  for (n = j;n > j - segment_length;n--) {
                    diff = real_sequence[index][m][n] - mean[m][k];
                    variance[m][k] += diff * diff;
                  }
                }

//                variance[m][k] /= segment_length;
                variance[m][k] /= (segment_length - 1);
              }
            }
          }

          if ( model_type[m - 1] == LINEAR_MODEL_CHANGE) {
            response_mean = 0.;
            if (type[m] != REAL_VALUE) {
              for (n = j;n > j - segment_length;n--) {
                response_mean += int_sequence[index][m][n];
              }
            }
            else {
              for (n = j;n > j - segment_length;n--) {
                response_mean += real_sequence[index][m][n];
              }
            }
            response_mean /= segment_length;

            index_parameter_mean = 0.;
            for (n = j;n > j - segment_length;n--) {
              index_parameter_mean += seq_index_parameter[n];
            }
            index_parameter_mean /= segment_length;

            index_parameter_variance = 0.;
            for (n = j;n > j - segment_length;n--) {
              diff = seq_index_parameter[n] - index_parameter_mean;
              index_parameter_variance += diff * diff;                
            }

            covariance = 0.;
            if (type[m] != REAL_VALUE) {
              for (n = j;n > j - segment_length;n--) {
                covariance += (int_sequence[index][m][n] - response_mean) * (seq_index_parameter[n] - index_parameter_mean);
              }
            }
            else {
              for (n = j;n > j - segment_length;n--) {
                covariance += (real_sequence[index][m][n] - response_mean) * (seq_index_parameter[n] - index_parameter_mean);
              }
            }

            slope[m][k] = covariance / index_parameter_variance;
            intercept[m][k] = response_mean - slope[m][k] * index_parameter_mean;

            residual_mean = 0.;
            if (type[m] != REAL_VALUE) {
              for (n = j;n > j - segment_length;n--) {
                residual_mean += int_sequence[index][m][n] - (intercept[m][k] + slope[m][k] * seq_index_parameter[n]);
              }
            }
            else {
              for (n = j;n > j - segment_length;n--) {
                residual_mean += real_sequence[index][m][n] - (intercept[m][k] + slope[m][k] * seq_index_parameter[n]);
              }
            }
            residual_mean /= segment_length;

            variance[m][k] = 0.;
            if (segment_length > 2) {
              if (type[m] != REAL_VALUE) {
                for (n = j;n > j - segment_length;n--) {
                  diff = int_sequence[index][m][n] - (intercept[m][k] + slope[m][k] * seq_index_parameter[n]) - residual_mean;
                  variance[m][k] += diff * diff;
                }
              }
              else {
                for (n = j;n > j - segment_length;n--) {
                  diff = real_sequence[index][m][n] - (intercept[m][k] + slope[m][k] * seq_index_parameter[n]) - residual_mean;
                  variance[m][k] += diff * diff;
                }
              }
              variance[m][k] /= (segment_length - 2);
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
          if ((model_type[j - 1] == POISSON_CHANGE) || (model_type[0] == MULTIVARIATE_POISSON_CHANGE) ||
              (model_type[j - 1] == BAYESIAN_POISSON_CHANGE)) {
            os << "# ";
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
            os << "# ";
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
                   (model_type[0] == MEAN_VARIANCE_CHANGE)|| (model_type[j - 1] == BAYESIAN_GAUSSIAN_CHANGE)) {
            os << "# ";
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

          else if (model_type[j - 1] == LINEAR_MODEL_CHANGE) {
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
                   (model_type[0] == MEAN_VARIANCE_CHANGE) || (model_type[j - 1] == BAYESIAN_GAUSSIAN_CHANGE)) {
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

          else if (model_type[j - 1] == LINEAR_MODEL_CHANGE) {
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

  delete [] norm;

  delete [] backward;
  delete [] cumul_backward;

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

# ifdef DEBUG
  for (i = 0;i < length[index];i++) {
    delete [] segment_probability[i];
  }
  delete [] segment_probability;
# endif

  return likelihood;
}


};  // namespace sequence_analysis
