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
 *       $Id: vector_algorithms.cpp 18022 2015-04-23 07:07:33Z guedon $
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

#include "stat_tools.h"
#include "markovian.h"
#include "vectors.h"
#include "distance_matrix.h"
#include "stat_label.h"

using namespace std;


namespace stat_tool {



/*--------------------------------------------------------------*
 *
 *  Calcul de la matrice des coefficients de correlation des rangs de Spearman.
 *
 *--------------------------------------------------------------*/

double** Vectors::spearman_rank_correlation_computation() const

{
  register int i , j , k;
  int *pfrequency;
  double main_term , rank_mean , rank_diff , *correction , **correlation = NULL , **rank;


  if (nb_vector > 2) {
    correlation = new double*[nb_variable];
    for (i = 0;i < nb_variable;i++) {
      correlation[i] = new double[nb_variable];
    }

    // calcul du terme principal et des termes de correction pour les ex-aequo

    main_term = nb_vector * ((double)nb_vector * (double)nb_vector - 1);

    correction = new double[nb_variable];
    for (i = 0;i < nb_variable;i++) {
      pfrequency = marginal_distribution[i]->frequency + marginal_distribution[i]->offset;
      correction[i] = 0.;
      for (j = marginal_distribution[i]->offset;j < marginal_distribution[i]->nb_value;j++) {
        if (*pfrequency > 1) {
          correction[i] += *pfrequency * ((double)*pfrequency * (double)*pfrequency - 1);
        }
        pfrequency++;
      }
    }

    // calcul des rangs

    rank_mean = (double)(nb_vector + 1) / 2.;

    rank = new double*[nb_variable];

    for (i = 0;i < nb_variable;i++) {
      rank[i] = marginal_distribution[i]->rank_computation();
    }

    for (i = 0;i < nb_variable;i++) {
      correlation[i][i] = 1.;
      for (j = i + 1;j < nb_variable;j++) {

        // calcul des differences de rangs centrees

        rank_diff = 0.;
        for (k = 0;k < nb_vector;k++) {
          rank_diff += (rank[i][int_vector[k][i]] - rank_mean) * (rank[j][int_vector[k][j]] - rank_mean);
        }

        correlation[i][j] = 12. * rank_diff /
                            sqrt((main_term - correction[i]) * (main_term - correction[j]));

        correlation[j][i] = correlation[i][j];
      }
    }

    delete [] correction;

    for (i = 0;i < nb_variable;i++) {
      delete [] rank[i];
    }
    delete [] rank;
  }

  return correlation;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la matrice des coefficients de correlation des rangs de Kendall.
 *
 *--------------------------------------------------------------*/

double** Vectors::kendall_rank_correlation_computation() const

{
  register int i , j , k , m , n , p;
  int diff , *current_frequency , *pfrequency , *cumul_frequency , *index , **frequency;
  double nb_pair , sum , rank_diff_sign , *correction , **correlation = NULL;


  if (nb_vector > 2) {
    correlation = new double*[nb_variable];
    for (i = 0;i < nb_variable;i++) {
      correlation[i] = new double[nb_variable];
    }

    // calcul du terme principal et des termes de correction pour les ex-aequo

    nb_pair = nb_vector * ((double)nb_vector - 1);

    correction = new double[nb_variable];
    for (i = 0;i < nb_variable;i++) {
      pfrequency = marginal_distribution[i]->frequency + marginal_distribution[i]->offset;
      correction[i] = 0.;
      for (j = marginal_distribution[i]->offset;j < marginal_distribution[i]->nb_value;j++) {
        if (*pfrequency > 1) {
          correction[i] += *pfrequency * ((double)*pfrequency - 1);
        }
        pfrequency++;
      }
    }

    index = new int[nb_vector];

    for (i = 0;i < nb_variable;i++) {
      correlation[i][i] = 1.;

      // calcul d'un ordre sur les vecteurs a partir des valeurs de la 1ere variable

      cumul_frequency = new int[marginal_distribution[i]->nb_value];
      current_frequency = new int[marginal_distribution[i]->nb_value];

      cumul_frequency[marginal_distribution[i]->offset] = 0;
      for (j = marginal_distribution[i]->offset + 1;j < marginal_distribution[i]->nb_value;j++) {
        cumul_frequency[j] = cumul_frequency[j - 1] + marginal_distribution[i]->frequency[j - 1];
      }

      for (j = marginal_distribution[i]->offset;j < marginal_distribution[i]->nb_value;j++) {
        current_frequency[j] = 0;
      }

      for (j = 0;j < nb_vector;j++) {
        index[cumul_frequency[int_vector[j][i]] + current_frequency[int_vector[j][i]]] = j;
        current_frequency[int_vector[j][i]]++;
      }

      for (j = i + 1;j < nb_variable;j++) {
        rank_diff_sign = 0.;

        if ((marginal_distribution[i]->nb_value - marginal_distribution[i]->offset) *
            (marginal_distribution[j]->nb_value - marginal_distribution[j]->offset) < nb_vector) {

          // calcul des frequences correspondant a la loi jointe des 2 variables

          frequency = joint_frequency_computation(i , j);

          // calcul des paires concordantes et discordantes

          for (k = marginal_distribution[i]->offset;k < marginal_distribution[i]->nb_value;k++) {
            for (m = marginal_distribution[j]->offset;m < marginal_distribution[j]->nb_value;m++) {
              if (frequency[k][m] > 0) {
                sum = 0.;
                for (n = k + 1;n < marginal_distribution[i]->nb_value;n++) {
                  for (p = marginal_distribution[j]->offset;p < m;p++) {
                    sum -= frequency[n][p];
                  }
                  for (p = m + 1;p < marginal_distribution[j]->nb_value;p++) {
                    sum += frequency[n][p];
                  }
                }
              }

              rank_diff_sign += frequency[k][m] * sum;
            }
          }

          for (k = 0;k < marginal_distribution[i]->nb_value;k++) {
            delete [] frequency[k];
          }
          delete [] frequency;
        }

        else {

          // calcul des signes des differences de rangs

          for (k = 0;k < nb_vector;k++) {
            for (m = cumul_frequency[int_vector[index[k]][i]];m < nb_vector;m++) {
              diff = int_vector[index[m]][j] - int_vector[index[k]][j];
              if (diff > 0) {
                rank_diff_sign++;
              }
              else if (diff < 0) {
                rank_diff_sign--;
              }
            }
          }
        }

        correlation[i][j] = 2. * rank_diff_sign /
                            sqrt((nb_pair - correction[i]) * (nb_pair - correction[j]));

        correlation[j][i] = correlation[i][j];
      }

      delete [] cumul_frequency;
      delete [] current_frequency;
    }

    delete [] correction;
    delete [] index;
  }

  return correlation;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'une matrice des coefficients de correlation des rangs.
 *
 *  arguments : stream, type de coefficient (SPEARMAN/KENDALL),
 *              pointeur sur les coefficients.
 *
 *--------------------------------------------------------------*/

ostream& Vectors::rank_correlation_ascii_write(ostream &os , int correlation_type ,
                                               double **correlation) const

{
  register int i , j;
  int buff , width[2];
  Test *test;
  long old_adjust;


  old_adjust = os.setf(ios::right , ios::adjustfield);

  width[0] = column_width(nb_variable);

  // calcul des largeurs des colonnes

  width[1] = 0;
  for (i = 0;i < nb_variable;i++) {
    buff = column_width(nb_variable , correlation[i]);
    if (buff > width[1]) {
      width[1] = buff;
    }
  }
  width[1] += ASCII_SPACE;

  // ecriture de la matrice des coefficients de correlation

  switch (correlation_type) {
  case SPEARMAN :
    os << STAT_label[STATL_SPEARMAN_RANK_CORRELATION_MATRIX] << endl;
    break;
  case KENDALL :
    os << STAT_label[STATL_KENDALL_RANK_CORRELATION_MATRIX] << endl;
    break;
  }

  os << "\n" << setw(width[0] + width[1]) << 1;
  for (i = 1;i < nb_variable;i++) {
    os << setw(width[1]) << i + 1;
  }
  for (i = 0;i < nb_variable;i++) {
    os << "\n" << setw(width[0]) << i + 1;
    for (j = 0;j < nb_variable;j++) {
      os << setw(width[1]) << correlation[i][j];
    }
  }
  os << endl;

  // test du caractere significatif des coefficients de correlation

  if (nb_vector > 2) {
    switch (correlation_type) {
    case SPEARMAN :
      test = new Test(STUDENT , false , nb_vector - 2 , I_DEFAULT , D_DEFAULT);
      break;
    case KENDALL :
      test = new Test(STANDARD_NORMAL , false , I_DEFAULT , I_DEFAULT , D_DEFAULT);
      break;
    }

    for (i = 0;i < NB_CRITICAL_PROBABILITY;i++) {
      test->critical_probability = ref_critical_probability[i];

      switch (test->ident) {
      case STANDARD_NORMAL :
        test->standard_normal_value_computation();
        break;
      case STUDENT :
        test->t_value_computation();
        break;
      }

      os << "\n" << STAT_label[STATL_REFERENCE] << " ";
      switch (test->ident) {
        case STANDARD_NORMAL :
        os << STAT_label[STATL_STANDARD_NORMAL_VALUE];
        break;
      case STUDENT :
        os << STAT_label[STATL_T_VALUE];
        break;
      }
      os << ": " << test->value << "   " << STAT_label[STATL_REFERENCE] << " "
         << STAT_label[STATL_CRITICAL_PROBABILITY] << ": " << test->critical_probability << endl;

      switch (correlation_type) {
      case SPEARMAN :
        os << STAT_label[STATL_SPEARMAN_LIMIT_RANK_CORRELATION_COEFF] << ": "
           << test->value / sqrt(test->value * test->value + nb_vector - 2) << endl;
        break;
      case KENDALL :
        os << STAT_label[STATL_KENDALL_LIMIT_RANK_CORRELATION_COEFF] << ": "
           << test->value * sqrt((2 * (2 * (double)nb_vector + 5)) / (9 * (double)nb_vector * (double)(nb_vector - 1))) << endl;
        break;
      }
    }

    delete test;
  }

  os.setf((FMTFLAGS)old_adjust , ios::adjustfield);

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'une matrice des coefficients de correlation des rangs dans un fichier.
 *
 *  arguments : reference sur un objet StatError, path, type de coefficient
 *              (SPEARMAN/KENDALL), pointeur sur les coefficients.
 *
 *--------------------------------------------------------------*/

bool Vectors::rank_correlation_ascii_write(StatError &error , const char *path ,
                                           int correlation_type , double **correlation) const

{
  bool status;
  ofstream out_file(path);


  error.init();

  if (!out_file) {
    status = false;
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {
    status = true;
    rank_correlation_ascii_write(out_file , correlation_type , correlation);
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la matrice des coefficients de correlation des rangs
 *  de Spearman ou de Kendall.
 *
 *  arguments : reference sur un objet StatError, stream,
 *              type de coefficient (SPEARMAN / KENDALL), path.
 *
 *--------------------------------------------------------------*/

bool Vectors::rank_correlation_computation(StatError &error , ostream &os ,
                                           int correlation_type , const char *path) const

{
  bool status = true;
  register int i;
  double **correlation;


  error.init();

  if (nb_vector <= 2) {
    status = false;
    ostringstream correction_message;
    correction_message << STAT_error[STATR_GREATER_THAN] << " " << 2;
    error.correction_update(STAT_error[STATR_NB_VECTOR] , (correction_message.str()).c_str());
  }

  for (i = 0;i < nb_variable;i++) {
    if (type[i] != INT_VALUE) {
      status = false;
      ostringstream error_message;
      error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                    << STAT_error[STATR_VARIABLE_TYPE];
      error.correction_update((error_message.str()).c_str() , STAT_variable_word[INT_VALUE]);
    }

    else if (!marginal_distribution[i]) {
      status = false;
      ostringstream error_message;
      error_message << STAT_error[STATR_RANK_CORRELATION_COMPUTATION] << ": "
                    << STAT_label[STATL_VARIABLE] << " " << i + 1 << " "
                    << STAT_error[STATR_SHIFTED_SCALED];
      error.update((error_message.str()).c_str());
    }
  }

  if (status) {
    switch (correlation_type) {
    case SPEARMAN :
      correlation = spearman_rank_correlation_computation();
      break;
    case KENDALL :
      correlation = kendall_rank_correlation_computation();
      break;
    }

#   ifdef MESSAGE
    rank_correlation_ascii_write(os , correlation_type , correlation);
#   endif

    if (path) {
      status = rank_correlation_ascii_write(error , path , correlation_type , correlation);
    }

    for (i = 0;i < nb_variable;i++) {
      delete [] correlation[i];
    }
    delete [] correlation;
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Calcul du coefficient de correlation des rangs de Spearman entre 2 variables.
 *
 *--------------------------------------------------------------*/

double Vectors::spearman_rank_single_correlation_computation() const

{
  register int i , j;
  int *pfrequency;
  double correlation , main_term , rank_mean , correction[2] , *rank[2];


  // calcul du terme principal et des termes de correction pour les ex-aequo

  main_term = nb_vector * ((double)nb_vector * (double)nb_vector - 1);

  for (i = 0;i < 2;i++) {
    pfrequency = marginal_distribution[i]->frequency + marginal_distribution[i]->offset;
    correction[i] = 0.;
    for (j = marginal_distribution[i]->offset;j < marginal_distribution[i]->nb_value;j++) {
      if (*pfrequency > 1) {
        correction[i] += *pfrequency * ((double)*pfrequency * (double)*pfrequency - 1);
      }
      pfrequency++;
    }
  }

  // calcul des rangs

  rank_mean = (double)(nb_vector + 1) / 2.;

  for (i = 0;i < 2;i++) {
    rank[i] = marginal_distribution[i]->rank_computation();
  }

  // calcul des differences de rangs centrees

  correlation = 0.;
  for (i = 0;i < nb_vector;i++) {
    correlation += (rank[0][int_vector[i][0]] - rank_mean) *
                   (rank[1][int_vector[i][1]] - rank_mean);
  }

  correlation = 12. * correlation / sqrt((main_term - correction[0]) * (main_term - correction[1]));

  for (i = 0;i < 2;i++) {
    delete [] rank[i];
  }

  return correlation;
}


/*--------------------------------------------------------------*
 *
 *  Calcul du coefficient de correlation des rangs de Kendall entre 2 variables.
 *
 *--------------------------------------------------------------*/

double Vectors::kendall_rank_single_correlation_computation() const

{
  register int i , j , k , m;
  int diff , *current_frequency , *pfrequency , *cumul_frequency , *index ,
      **frequency;
  double sum , correlation , nb_pair , correction[2];


  // calcul du terme principal et des termes de correction pour les ex-aequo

  nb_pair = nb_vector * ((double)nb_vector - 1);

  for (i = 0;i < 2;i++) {
    pfrequency = marginal_distribution[i]->frequency + marginal_distribution[i]->offset;
    correction[i] = 0.;
    for (j = marginal_distribution[i]->offset;j < marginal_distribution[i]->nb_value;j++) {
      if (*pfrequency > 1) {
        correction[i] += *pfrequency * ((double)*pfrequency - 1);
      }
      pfrequency++;
    }
  }

  correlation = 0.;

  if ((marginal_distribution[0]->nb_value - marginal_distribution[0]->offset) *
      (marginal_distribution[1]->nb_value - marginal_distribution[1]->offset) < nb_vector) {

    // calcul des frequences correspondant a la loi jointe des 2 variables

    frequency = joint_frequency_computation(0 , 1);

    // calcul des paires concordantes et discordantes

    for (i = marginal_distribution[0]->offset;i < marginal_distribution[0]->nb_value;i++) {
      for (j = marginal_distribution[1]->offset;j < marginal_distribution[1]->nb_value;j++) {
        if (frequency[i][j] > 0) {
          sum = 0.;
          for (k = i + 1;k < marginal_distribution[0]->nb_value;k++) {
            for (m = marginal_distribution[1]->offset;m < j;m++) {
              sum -= frequency[k][m];
            }
            for (m = j + 1;m < marginal_distribution[1]->nb_value;m++) {
              sum += frequency[k][m];
            }
          }
        }

        correlation += frequency[i][j] * sum;
      }
    }

    for (i = 0;i < marginal_distribution[0]->nb_value;i++) {
      delete [] frequency[i];
    }
    delete [] frequency;
  }

  else {

    // calcul d'un ordre sur les vecteurs a partir des valeurs de la 1ere variable

    cumul_frequency = new int[marginal_distribution[0]->nb_value];
    current_frequency = new int[marginal_distribution[0]->nb_value];
    index = new int[nb_vector];

    cumul_frequency[marginal_distribution[0]->offset] = 0;
    for (i = marginal_distribution[0]->offset + 1;i < marginal_distribution[0]->nb_value;i++) {
      cumul_frequency[i] = cumul_frequency[i - 1] + marginal_distribution[0]->frequency[i - 1];
    }

    for (i = marginal_distribution[0]->offset;i < marginal_distribution[0]->nb_value;i++) {
      current_frequency[i] = 0;
    }

    for (i = 0;i < nb_vector;i++) {
      index[cumul_frequency[int_vector[i][0]] + current_frequency[int_vector[i][0]]] = i;
      current_frequency[int_vector[i][0]]++;
    }

    // calcul des signes des differences de rangs

    for (i = 0;i < nb_vector;i++) {
      for (j = cumul_frequency[int_vector[index[i]][0]];j < nb_vector;j++) {
        diff = int_vector[index[j]][1] - int_vector[index[i]][1];
        if (diff > 0) {
          correlation++;
        }
        else if (diff < 0) {
          correlation--;
        }
      }
    }

    delete [] cumul_frequency;
    delete [] current_frequency;
    delete [] index;
  }

  correlation = 2. * correlation / sqrt((nb_pair - correction[0]) * (nb_pair - correction[1]));

  return correlation;
}


/*--------------------------------------------------------------*
 *
 *  Comparaison de vecteurs.
 *
 *  arguments : references sur un objet StatError et sur un objet VectorDistance,
 *              flag standardisation (uniquement variables de meme type).
 *
 *--------------------------------------------------------------*/

DistanceMatrix* Vectors::comparison(StatError &error , const VectorDistance &ivector_dist ,
                                    bool standardization) const
{
  bool status = true;
  register int i , j , k;
  double distance , ldistance , **rank;
  FrequencyDistribution *merged_marginal;
  VectorDistance *vector_dist;
  DistanceMatrix *dist_matrix;


  dist_matrix = NULL;
  error.init();

  if ((nb_vector < 2) || (nb_vector > DISTANCE_NB_VECTOR)) {
    status = false;
    error.update(STAT_error[STATR_NB_VECTOR]);
  }

  if (ivector_dist.nb_variable != nb_variable) {
    status = false;
    error.update(STAT_error[STATR_NB_VARIABLE]);
  }

  else {
    for (i = 0;i < nb_variable;i++) {
      if (ivector_dist.variable_type[i] != NUMERIC) {
        if (type[i] != INT_VALUE) {
          status = false;
          ostringstream error_message;
          error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                        << STAT_error[STATR_VARIABLE_TYPE];
          error.correction_update((error_message.str()).c_str() , STAT_variable_word[INT_VALUE]);
        }

        else if (!marginal_distribution[i]) {
          status = false;
          ostringstream error_message;
          error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                        << STAT_error[STATR_MARGINAL_FREQUENCY_DISTRIBUTION];
          error.update((error_message.str()).c_str());
        }

        if ((ivector_dist.variable_type[i] == SYMBOLIC) &&
            ((min_value[i] < 0) || (max_value[i] >= NB_SYMBOL) ||
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

    if (!standardization) {
      for (i = 1;i < nb_variable;i++) {
        if (ivector_dist.variable_type[i] != ivector_dist.variable_type[0]) {
          status = false;
          ostringstream error_message;
          error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                        << STAT_error[STATR_VARIABLE_TYPE];
          error.correction_update((error_message.str()).c_str() , STAT_variable_type_word[0]);
        }
      }
    }
  }

  if (status) {
    vector_dist = new VectorDistance(ivector_dist);

    // calcul des rangs pour les variables ordinales

    rank = new double*[nb_variable];

    switch (standardization) {

    case false : {
      if (vector_dist->variable_type[0] == ORDINAL) {
        merged_marginal = new FrequencyDistribution(nb_variable , (const FrequencyDistribution**)marginal_distribution);
        rank[0] = merged_marginal->rank_computation();

#       ifdef MESSAGE
        cout << "\nRanks:";
        for (i = merged_marginal->offset;i < merged_marginal->nb_value;i++) {
          cout << " " << rank[0][i];
        }
        cout << endl;
#       endif

        for (i = 1;i < nb_variable;i++) {
          rank[i] = new double[merged_marginal->nb_value];
          for (j = merged_marginal->offset;j < merged_marginal->nb_value;j++) {
            rank[i][j] = rank[0][j];
          }
        }

        delete merged_marginal;
      }

      else {
        for (i = 0;i < nb_variable;i++) {
          rank[i] = NULL;
        }
      }

#     ifdef DEBUG
      for (i = 0;i < nb_variable;i++) {
        vector_dist->weight[i] = 1.;
      }
#     endif

      break;
    }

    case true : {
      for (i = 0;i < nb_variable;i++) {
        if (vector_dist->variable_type[i] == ORDINAL) {
          rank[i] = marginal_distribution[i]->rank_computation();
        }
        else {
          rank[i] = NULL;
        }
      }
      break;
    }
    }

    // calcul des dispersions pour la standardisation

    if (standardization) {
      for (i = 0;i < nb_variable;i++) {
        if (marginal_distribution[i]) {
          vector_dist->dispersion_computation(i , marginal_distribution[i] , rank[i]);
        }

        else {
          switch (vector_dist->distance_type) {
          case ABSOLUTE_VALUE :
            vector_dist->dispersion[i] = mean_absolute_difference_computation(i);
            break;
          case QUADRATIC :
            vector_dist->dispersion[i] = 2 * covariance[i][i];
            break;
          }

          if (vector_dist->dispersion[i] == 0.) {
            vector_dist->dispersion[i] = 1.;
          }
        }
      }
    }

#   ifdef DEBUG
    double *variable_distance;

    variable_distance = new double[nb_variable];
    for (i = 0;i < nb_variable;i++) {
      variable_distance[i] = 0.;
    }

    cout << *vector_dist;
    if (vector_dist->distance_type == ABSOLUTE_VALUE) {
      for (i = 0;i < nb_variable;i++) {
        if (vector_dist->variable_type[i] == NUMERIC) {
          cout << "\n" << STAT_label[STATL_VARIABLE] << " " << i << "   mean absolute difference: "
               << mean_absolute_difference_computation(i) << endl;
        }
      }
    }
#   endif

    dist_matrix = new DistanceMatrix(nb_vector , STAT_label[STATL_VECTOR] , identifier);

    for (i = 0;i < nb_vector;i++) {
      for (j = i + 1;j < nb_vector;j++) {
        distance = 0.;

        for (k = 0;k < vector_dist->nb_variable;k++) {
          switch (vector_dist->variable_type[k]) {

          case SYMBOLIC : {
            if (!(vector_dist->symbol_distance[k])) {
              if (int_vector[i][k] == int_vector[j][k]) {
                ldistance = 0.;
              }
              else {
                ldistance = 1.;
              }
            }

            else {
              ldistance = vector_dist->symbol_distance[k][int_vector[i][k]][int_vector[j][k]];
            }
            break;
          }

          case ORDINAL : {
            ldistance = rank[k][int_vector[i][k]] - rank[k][int_vector[j][k]];
            break;
          }

          case NUMERIC : {
            if (type[k] == INT_VALUE) {
              ldistance = int_vector[i][k] - int_vector[j][k];
            }
            else {
              ldistance = real_vector[i][k] - real_vector[j][k];
            }
            break;
          }

          case CIRCULAR : {
            if (int_vector[i][k] <= int_vector[j][k]) {
              ldistance = MIN(int_vector[j][k] - int_vector[i][k] , int_vector[i][k] + vector_dist->period[k] - int_vector[j][k]);
            }
            else {
              ldistance = MIN(int_vector[i][k] - int_vector[j][k] , int_vector[j][k] + vector_dist->period[k] - int_vector[i][k]);
            }
            break;
          }
          }

          switch (standardization) {

          case false : {

#           ifdef DEBUG
            switch (vector_dist->distance_type) {
            case ABSOLUTE_VALUE :
              variable_distance[k] += fabs(ldistance);
              break;
            case QUADRATIC :
              variable_distance[k] += ldistance * ldistance;
              break;
            }
#           endif

            switch (vector_dist->distance_type) {
            case ABSOLUTE_VALUE :
              distance += vector_dist->weight[k] * fabs(ldistance);
              break;
            case QUADRATIC :
              distance += vector_dist->weight[k] * ldistance * ldistance;
              break;
            }
            break;
          }

          case true : {

#           ifdef DEBUG
            switch (vector_dist->distance_type) {
            case ABSOLUTE_VALUE :
              variable_distance[k] += fabs(ldistance) / vector_dist->dispersion[k];
              break;
            case QUADRATIC :
              variable_distance[k] += ldistance * ldistance / vector_dist->dispersion[k];
              break;
            }
#           endif

            switch (vector_dist->distance_type) {
            case ABSOLUTE_VALUE :
              distance += vector_dist->weight[k] * fabs(ldistance) / vector_dist->dispersion[k];
              break;
            case QUADRATIC :
              distance += vector_dist->weight[k] * ldistance * ldistance / vector_dist->dispersion[k];
              break;
            }
            break;
          }
          }
        }

        if (vector_dist->distance_type == QUADRATIC) {
          distance = sqrt(distance);
        }

        dist_matrix->update(identifier[i] , identifier[j] , distance , 1);
        dist_matrix->update(identifier[j] , identifier[i] , distance , 1);
      }
    }

#   ifdef DEBUG
    cout << "\nAverage distances" << endl;
    for (i = 0;i < nb_variable;i++) {
      cout << STAT_label[STATL_VARIABLE] << " " << i << ": "
           << 2. * variable_distance[i] / (nb_vector * (nb_vector - 1)) << endl;
    }
    delete [] variable_distance;
#   endif

    delete vector_dist;

    for (i = 0;i < nb_variable;i++) {
      delete [] rank[i];
    }
    delete [] rank;
  }

  return dist_matrix;
}


/*--------------------------------------------------------------*
 *
 *  Calcul des frequences correspondant a la loi jointe de 2 variables.
 *
 *  arguments : indices des 2 variables.
 *
 *--------------------------------------------------------------*/

int** Vectors::joint_frequency_computation(int variable1 , int variable2) const

{
  register int i , j;
  int **frequency = NULL;


  if ((marginal_distribution[variable1]) && (marginal_distribution[variable2])) {
    frequency = new int*[marginal_distribution[variable1]->nb_value];
    for (i = 0;i < marginal_distribution[variable1]->nb_value;i++) {
      frequency[i] = new int[marginal_distribution[variable2]->nb_value];
      for (j = 0;j < marginal_distribution[variable2]->nb_value;j++) {
        frequency[i][j] = 0;
      }
    }

    for (i = 0;i < nb_vector;i++) {
      frequency[int_vector[i][variable1]][int_vector[i][variable2]]++;
    }
  }

  return frequency;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un tableau de contingence.
 *
 *  arguments : stream, indices des 2 variables, pointeurs sur les tableaux de
 *              contingence, d'ecarts et de contributions au chi2,
 *              reference sur le resultat du test du chi2, flag fichier.
 *
 *--------------------------------------------------------------*/

ostream& Vectors::contingency_table_ascii_write(ostream &os , int variable1 , int variable2 ,
                                                int **frequency , double **deviation ,
                                                double **chi2_contribution , const Test &test ,
                                                bool file_flag) const

{
  register int i , j;
  int buff , width[2];
  long old_adjust;


  if ((file_flag) || ((marginal_distribution[variable1]->nb_value - marginal_distribution[variable1]->offset <= DISPLAY_CONTINGENCY_NB_VALUE) &&
       (marginal_distribution[variable2]->nb_value - marginal_distribution[variable2]->offset <= DISPLAY_CONTINGENCY_NB_VALUE))) {
    old_adjust = os.setf(ios::right , ios::adjustfield);

    // calcul des largeurs des colonnes

    width[0] = column_width(MAX(marginal_distribution[variable1]->nb_value , marginal_distribution[variable2]->nb_value) - 1);
    width[1] = column_width(nb_vector) + ASCII_SPACE;

    // ecriture du tableau de contingence

    os << STAT_label[STATL_CONTINGENCY_TABLE] << endl;

    os << "\n" << setw(width[0] + width[1]) << marginal_distribution[variable2]->offset;
    for (i = marginal_distribution[variable2]->offset + 1;i < marginal_distribution[variable2]->nb_value;i++) {
      os << setw(width[1]) << i;
    }
    for (i = marginal_distribution[variable1]->offset;i < marginal_distribution[variable1]->nb_value;i++) {
      os << "\n" << setw(width[0]) << i;
      for (j = marginal_distribution[variable2]->offset;j < marginal_distribution[variable2]->nb_value;j++) {
        os << setw(width[1]) << frequency[i][j];
      }
      os << setw(width[1]) << marginal_distribution[variable1]->frequency[i];
    }
    os << "\n" << setw(width[0] + width[1]) << marginal_distribution[variable2]->frequency[marginal_distribution[variable2]->offset];
    for (i = marginal_distribution[variable2]->offset + 1;i < marginal_distribution[variable2]->nb_value;i++) {
      os << setw(width[1]) << marginal_distribution[variable2]->frequency[i];
    }
    os << setw(width[1]) << nb_vector << endl;

    // calcul des largeurs des colonnes

    width[1] = 0;
    for (i = marginal_distribution[variable1]->offset;i < marginal_distribution[variable1]->nb_value;i++) {
      buff = column_width(marginal_distribution[variable2]->nb_value - marginal_distribution[variable2]->offset ,
                          deviation[i] + marginal_distribution[variable2]->offset);
      if (buff > width[1]) {
        width[1] = buff;
      }
    }
    width[1] += ASCII_SPACE;

    // ecriture du tableau des ecarts

    os << "\n" << STAT_label[STATL_DEVIATION_TABLE] << endl;

    os << "\n" << setw(width[0] + width[1]) << marginal_distribution[variable2]->offset;
    for (i = marginal_distribution[variable2]->offset + 1;i < marginal_distribution[variable2]->nb_value;i++) {
      os << setw(width[1]) << i;
    }
    for (i = marginal_distribution[variable1]->offset;i < marginal_distribution[variable1]->nb_value;i++) {
      os << "\n" << setw(width[0]) << i;
      for (j = marginal_distribution[variable2]->offset;j < marginal_distribution[variable2]->nb_value;j++) {
        os << setw(width[1]) << deviation[i][j];
      }
    }
    os << endl;

    // calcul des largeurs des colonnes

    width[1] = 0;
    for (i = marginal_distribution[variable1]->offset;i < marginal_distribution[variable1]->nb_value;i++) {
      buff = column_width(marginal_distribution[variable2]->nb_value - marginal_distribution[variable2]->offset ,
                          chi2_contribution[i] + marginal_distribution[variable2]->offset);
      if (buff > width[1]) {
        width[1] = buff;
      }
    }
    width[1] += ASCII_SPACE;

    // ecriture du tableau de contributions au chi2

    if (test.value > 0.) {
      os << "\n" << STAT_label[STATL_CHI2_CONTRBUTION_TABLE] << endl;

      os << "\n" << setw(width[0] + width[1]) << marginal_distribution[variable2]->offset;
      for (i = marginal_distribution[variable2]->offset + 1;i < marginal_distribution[variable2]->nb_value;i++) {
        os << setw(width[1]) << i;
      }
      for (i = marginal_distribution[variable1]->offset;i < marginal_distribution[variable1]->nb_value;i++) {
        os << "\n" << setw(width[0]) << i;
        for (j = marginal_distribution[variable2]->offset;j < marginal_distribution[variable2]->nb_value;j++) {
          os << setw(width[1]) << chi2_contribution[i][j];
        }
      }
      os << endl;
    }

    os.setf((FMTFLAGS)old_adjust , ios::adjustfield);
  }

  os << "\n";
  test.ascii_print(os);

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un tableau de contingence dans un fichier.
 *
 *  arguments : reference sur un objet StatError, path, indices des 2 variables,
 *              pointeurs sur les tableaux de contingence, d'ecarts et
 *              de contributions au chi2, reference sur le resultat du test du chi2.
 *
 *--------------------------------------------------------------*/

bool Vectors::contingency_table_ascii_write(StatError &error , const char *path ,
                                            int variable1 , int variable2 , int **frequency ,
                                            double **deviation , double **chi2_contribution ,
                                            const Test &test) const

{
  bool status;
  ofstream out_file(path);


  error.init();

  if (!out_file) {
    status = false;
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {
    status = true;
    contingency_table_ascii_write(out_file , variable1 , variable2 , frequency ,
                                  deviation , chi2_contribution , test , true);
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un tableau de contingence dans un fichier au format tableur.
 *
 *  arguments : reference sur un objet StatError, path, indices des 2 variables,
 *              pointeurs sur les tableaux de contingence, d'ecarts et
 *              de contributions au chi2, reference sur le resultat du test du chi2.
 *
 *--------------------------------------------------------------*/

bool Vectors::contingency_table_spreadsheet_write(StatError &error , const char *path ,
                                                  int variable1 , int variable2 , int **frequency ,
                                                  double **deviation , double **chi2_contribution ,
                                                  const Test &test) const

{
  bool status;
  register int i , j;
  ofstream out_file(path);


  error.init();

  if (!out_file) {
    status = false;
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {
    status = true;

    // ecriture du tableau de contingence

    out_file << STAT_label[STATL_CONTINGENCY_TABLE] << endl;

    out_file << "\n";
    for (i = marginal_distribution[variable2]->offset;i < marginal_distribution[variable2]->nb_value;i++) {
      out_file << "\t" << i;
    }
    for (i = marginal_distribution[variable1]->offset;i < marginal_distribution[variable1]->nb_value;i++) {
      out_file << "\n" << i;
      for (j = marginal_distribution[variable2]->offset;j < marginal_distribution[variable2]->nb_value;j++) {
        out_file << "\t" << frequency[i][j];
      }
      out_file << "\t" << marginal_distribution[variable1]->frequency[i];
    }
    out_file << "\n";
    for (i = marginal_distribution[variable2]->offset;i < marginal_distribution[variable2]->nb_value;i++) {
      out_file << "\t" << marginal_distribution[variable2]->frequency[i];
    }
    out_file << "\t" << nb_vector << endl;

    // ecriture du tableau des ecarts

    out_file << "\n" << STAT_label[STATL_DEVIATION_TABLE] << endl;

    out_file << "\n";
    for (i = marginal_distribution[variable2]->offset;i < marginal_distribution[variable2]->nb_value;i++) {
      out_file << "\t" << i;
    }
    for (i = marginal_distribution[variable1]->offset;i < marginal_distribution[variable1]->nb_value;i++) {
      out_file << "\n" << i;
      for (j = marginal_distribution[variable2]->offset;j < marginal_distribution[variable2]->nb_value;j++) {
        out_file << "\t" << deviation[i][j];
      }
    }
    out_file << endl;

    // ecriture du tableau des contributions au chi2

    if (test.value > 0.) {
      out_file << "\n" << STAT_label[STATL_CHI2_CONTRBUTION_TABLE] << endl;

      out_file << "\n";
      for (i = marginal_distribution[variable2]->offset;i < marginal_distribution[variable2]->nb_value;i++) {
        out_file << "\t" << i;
      }
      for (i = marginal_distribution[variable1]->offset;i < marginal_distribution[variable1]->nb_value;i++) {
        out_file << "\n" << i;
        for (j = marginal_distribution[variable2]->offset;j < marginal_distribution[variable2]->nb_value;j++) {
          out_file << "\t" << chi2_contribution[i][j];
        }
      }
      out_file << endl;
    }

    out_file << "\n";
    test.spreadsheet_print(out_file);
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Tableau de contingence entre 2 variables.
 *
 *  arguments : reference sur un objet StatError, stream, indices des 2 variables,
 *              path, format de fichier ('a' : ASCII, 's' : Spreadsheet).
 *
 *--------------------------------------------------------------*/

bool Vectors::contingency_table(StatError &error , ostream &os , int variable1 ,
                                int variable2 , const char *path , char format) const

{
  bool status = true;
  register int i , j;
  int df , **frequency;
  double value , var , **deviation , **chi2_contribution;
  Test *test;


  error.init();

  if (variable1 == variable2) {
    status = false;
    error.update(STAT_error[STATR_VARIABLE_INDICES]);
  }

  if ((variable1 < 1) || (variable1 > nb_variable)) {
    status = false;
    ostringstream error_message;
    error_message << STAT_label[STATL_VARIABLE] << " " << variable1 << ": "
                  << STAT_error[STATR_VARIABLE_INDEX];
    error.update((error_message.str()).c_str());
  }

  else {
    variable1--;

    if (type[variable1] != INT_VALUE) {
      status = false;
      ostringstream error_message;
      error_message << STAT_label[STATL_VARIABLE] << " " << variable1 + 1 << ": "
                    << STAT_error[STATR_VARIABLE_TYPE];
      error.correction_update((error_message.str()).c_str() , STAT_variable_word[INT_VALUE]);
    }

    else if ((!marginal_distribution[variable1]) ||
             (marginal_distribution[variable1]->nb_value - marginal_distribution[variable1]->offset > CONTINGENCY_NB_VALUE)) {
      status = false;
      ostringstream error_message;
      error_message << STAT_label[STATL_VARIABLE] << " " << variable1 + 1 << ": "
                    << STAT_error[STATR_NB_VALUE];
      error.update((error_message.str()).c_str());
    }
  }

  if ((variable2 < 1) || (variable2 > nb_variable)) {
    status = false;
    ostringstream error_message;
    error_message << STAT_label[STATL_VARIABLE] << " " << variable2 << ": "
                  << STAT_error[STATR_VARIABLE_INDEX];
    error.update((error_message.str()).c_str());
  }

  else {
    variable2--;

    if (type[variable2] != INT_VALUE) {
      status = false;
      ostringstream error_message;
      error_message << STAT_label[STATL_VARIABLE] << " " << variable2 + 1 << ": "
                    << STAT_error[STATR_VARIABLE_TYPE];
      error.correction_update((error_message.str()).c_str() , STAT_variable_word[INT_VALUE]);
    }

    else if ((!marginal_distribution[variable2]) ||
             (marginal_distribution[variable2]->nb_value - marginal_distribution[variable2]->offset > CONTINGENCY_NB_VALUE)) {
      status = false;
      ostringstream error_message;
      error_message << STAT_label[STATL_VARIABLE] << " " << variable2 + 1 << ": "
                    << STAT_error[STATR_NB_VALUE];
      error.update((error_message.str()).c_str());
    }
  }

  if (status) {

    // calcul des frequences correspondant a la loi jointe des 2 variables

    frequency = joint_frequency_computation(variable1 , variable2);

    deviation = new double*[marginal_distribution[variable1]->nb_value];
    for (i = 0;i < marginal_distribution[variable1]->nb_value;i++) {
      deviation[i] = new double[marginal_distribution[variable2]->nb_value];
      for (j = 0;j < marginal_distribution[variable2]->nb_value;j++) {
        deviation[i][j] = 0.;
      }
    }

    chi2_contribution = new double*[marginal_distribution[variable1]->nb_value];
    for (i = 0;i < marginal_distribution[variable1]->nb_value;i++) {
      chi2_contribution[i] = new double[marginal_distribution[variable2]->nb_value];
      for (j = 0;j < marginal_distribution[variable2]->nb_value;j++) {
        chi2_contribution[i][j] = 0.;
      }
    }

    // calcul de la valeur du chi2

    df = 1;
    value = 0.;

    for (i = marginal_distribution[variable1]->offset;i < marginal_distribution[variable1]->nb_value;i++) {
      if (marginal_distribution[variable1]->frequency[i] > 0) {
        for (j = marginal_distribution[variable2]->offset;j < marginal_distribution[variable2]->nb_value;j++) {
          if (marginal_distribution[variable2]->frequency[j] > 0) {
            var = (double)(marginal_distribution[variable1]->frequency[i] * marginal_distribution[variable2]->frequency[j]) /
                  (double)nb_vector;

            df++;
            deviation[i][j] = frequency[i][j] - var;
            chi2_contribution[i][j] = deviation[i][j] * deviation[i][j] / var;
            value += chi2_contribution[i][j];
          }
        }

        df--;
      }
    }

    for (i = marginal_distribution[variable2]->offset;i < marginal_distribution[variable2]->nb_value;i++) {
      if (marginal_distribution[variable2]->frequency[i] > 0) {
        df--;
      }
    }

    if (value > 0.) {
      for (i = marginal_distribution[variable1]->offset;i < marginal_distribution[variable1]->nb_value;i++) {
        if (marginal_distribution[variable1]->frequency[i] > 0) {
          for (j = marginal_distribution[variable2]->offset;j < marginal_distribution[variable2]->nb_value;j++) {
            if (marginal_distribution[variable2]->frequency[j] > 0) {
              chi2_contribution[i][j] /= value;
            }
          }
        }
      }
    }

    test = new Test(CHI2 , true , df , I_DEFAULT , value);
    test->chi2_critical_probability_computation();

#   ifdef MESSAGE
    contingency_table_ascii_write(os , variable1 , variable2 , frequency ,
                                  deviation , chi2_contribution , *test);
#   endif

    if (path) {
      switch (format) {
      case 'a' :
        status = contingency_table_ascii_write(error , path , variable1 , variable2 , frequency ,
                                               deviation , chi2_contribution , *test);
        break;
      case 's' :
        status = contingency_table_spreadsheet_write(error , path , variable1 , variable2 , frequency ,
                                                     deviation , chi2_contribution , *test);
        break;
      }
    }

    for (i = 0;i < marginal_distribution[variable1]->nb_value;i++) {
      delete [] frequency[i];
      delete [] deviation[i];
      delete [] chi2_contribution[i];
    }
    delete [] frequency;
    delete [] deviation;
    delete [] chi2_contribution;

    delete test;
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture des resultats d'une analyse de variance a un facteur.
 *
 *  arguments : stream, type de la variable reponse (ORDINAL/NUMERIC),
 *              pointeurs sur les sous-echantillons pour chaque niveau possible,
 *              flag niveau de detail.
 *
 *--------------------------------------------------------------*/

ostream& Vectors::variance_analysis_ascii_write(ostream &os , int type , const Vectors **value_vec ,
                                                bool exhaustive) const

{
  register int i , j;
  int buff , width[4] , df[3];
  double diff , square_sum[3] , mean_square[3] , *value_mean , *variance , *standard_deviation ,
         *mean_absolute_deviation , *concentration_coeff , *skewness_coeff , *kurtosis_coeff , **cumul;
  Test *test;
  long old_adjust;


  old_adjust = os.setf(ios::right , ios::adjustfield);

  // lois conditionnelles

  value_mean = new double[marginal_distribution[0]->nb_value];
  variance = new double[marginal_distribution[0]->nb_value];
  standard_deviation = new double[marginal_distribution[0]->nb_value];
  mean_absolute_deviation = new double[marginal_distribution[0]->nb_value];
  concentration_coeff = new double[marginal_distribution[0]->nb_value];
  skewness_coeff = new double[marginal_distribution[0]->nb_value];
  kurtosis_coeff = new double[marginal_distribution[0]->nb_value];

  for (i = marginal_distribution[0]->offset;i < marginal_distribution[0]->nb_value;i++) {
    if (marginal_distribution[0]->frequency[i] > 0) {
      value_mean[i] = value_vec[i]->mean[1];
      variance[i] = value_vec[i]->covariance[1][1];
      standard_deviation[i] = sqrt(value_vec[i]->covariance[1][1]);
      mean_absolute_deviation[i] = value_vec[i]->mean_absolute_deviation_computation(1);
      if (value_vec[i]->nb_vector > 1) {
        concentration_coeff[i] = (value_vec[i]->mean_absolute_difference_computation(1) / (2 * value_vec[i]->mean[1])) *
                                 ((double)(value_vec[i]->nb_vector - 1) / (double)value_vec[i]->nb_vector);
      }
      else {
        concentration_coeff[i] = 1.;
      }
      skewness_coeff[i] = value_vec[i]->skewness_computation(1);
      kurtosis_coeff[i] = value_vec[i]->kurtosis_computation(1);
    }

    else {
      value_mean[i] = 0.;
      variance[i] = 0.;
      standard_deviation[i] = 0.;
      mean_absolute_deviation[i] = 0.;
      concentration_coeff[i] = 0.;
      skewness_coeff[i] = 0.;
      kurtosis_coeff[i] = 0.;
    }
  }

  width[0] = column_width(marginal_distribution[0]->nb_value - 1);
  buff = column_width(marginal_distribution[0]->max);
  if (buff > width[0]) {
    width[0] = buff;
  }
  buff = column_width(marginal_distribution[0]->nb_value - marginal_distribution[0]->offset , value_mean + marginal_distribution[0]->offset);
  if (buff > width[0]) {
    width[0] = buff;
  }
  buff = column_width(marginal_distribution[0]->nb_value - marginal_distribution[0]->offset , variance + marginal_distribution[0]->offset);
  if (buff > width[0]) {
    width[0] = buff;
  }
  buff = column_width(marginal_distribution[0]->nb_value - marginal_distribution[0]->offset , standard_deviation + marginal_distribution[0]->offset);
  if (buff > width[0]) {
    width[0] = buff;
  }
  buff = column_width(marginal_distribution[0]->nb_value - marginal_distribution[0]->offset , mean_absolute_deviation + marginal_distribution[0]->offset);
  if (buff > width[0]) {
    width[0] = buff;
  }
  buff = column_width(marginal_distribution[0]->nb_value - marginal_distribution[0]->offset , concentration_coeff + marginal_distribution[0]->offset);
  if (buff > width[0]) {
    width[0] = buff;
  }
  buff = column_width(marginal_distribution[0]->nb_value - marginal_distribution[0]->offset , skewness_coeff + marginal_distribution[0]->offset);
  if (buff > width[0]) {
    width[0] = buff;
  }
  buff = column_width(marginal_distribution[0]->nb_value - marginal_distribution[0]->offset , kurtosis_coeff + marginal_distribution[0]->offset);
  if (buff > width[0]) {
    width[0] = buff;
  }
  width[0] += ASCII_SPACE;

  os << STAT_label[STATL_VALUE] << "                       ";
  for (i = marginal_distribution[0]->offset;i < marginal_distribution[0]->nb_value;i++) {
    if (marginal_distribution[0]->frequency[i] > 0) {
      os << setw(width[0]) << i;
    }
  }
  os << endl;

  os << STAT_label[STATL_SAMPLE_SIZE] << "                 ";
  for (i = marginal_distribution[0]->offset;i < marginal_distribution[0]->nb_value;i++) {
    if (marginal_distribution[0]->frequency[i] > 0) {
      os << setw(width[0]) << marginal_distribution[0]->frequency[i];
    }
  }
  os << endl;

  os << STAT_label[STATL_MEAN] << "                        ";
  for (i = marginal_distribution[0]->offset;i < marginal_distribution[0]->nb_value;i++) {
    if (marginal_distribution[0]->frequency[i] > 0) {
      os << setw(width[0]) << value_mean[i];
    }
  }
  os << endl;

  os << STAT_label[STATL_VARIANCE] << "                    ";
  for (i = marginal_distribution[0]->offset;i < marginal_distribution[0]->nb_value;i++) {
    if (marginal_distribution[0]->frequency[i] > 0) {
      os << setw(width[0]) << variance[i];
    }
  }
  os << endl;

  os << STAT_label[STATL_STANDARD_DEVIATION] << "          ";
  for (i = marginal_distribution[0]->offset;i < marginal_distribution[0]->nb_value;i++) {
    if (marginal_distribution[0]->frequency[i] > 0) {
      os << setw(width[0]) << standard_deviation[i];
    }
  }
  os << endl;

  os << STAT_label[STATL_MEAN_ABSOLUTE_DEVIATION] << "     ";
  for (i = marginal_distribution[0]->offset;i < marginal_distribution[0]->nb_value;i++) {
    if (marginal_distribution[0]->frequency[i] > 0) {
      os << setw(width[0]) << mean_absolute_deviation[i];
    }
  }
  os << endl;

  os << STAT_label[STATL_CONCENTRATION_COEFF];
  for (i = marginal_distribution[0]->offset;i < marginal_distribution[0]->nb_value;i++) {
    if (marginal_distribution[0]->frequency[i] > 0) {
      os << setw(width[0]) << concentration_coeff[i];
    }
  }
  os << endl;

  os << STAT_label[STATL_SKEWNESS_COEFF] << "     ";
  for (i = marginal_distribution[0]->offset;i < marginal_distribution[0]->nb_value;i++) {
    if (marginal_distribution[0]->frequency[i] > 0) {
      os << setw(width[0]) << skewness_coeff[i];
    }
  }
  os << endl;

  os << STAT_label[STATL_KURTOSIS_COEFF] << "     ";
  for (i = marginal_distribution[0]->offset;i < marginal_distribution[0]->nb_value;i++) {
    if (marginal_distribution[0]->frequency[i] > 0) {
      os << setw(width[0]) << kurtosis_coeff[i];
    }
  }
  os << endl;

  delete [] value_mean;
  delete [] variance;
  delete [] standard_deviation;
  delete [] mean_absolute_deviation;
  delete [] concentration_coeff;
  delete [] skewness_coeff;
  delete [] kurtosis_coeff;

  if ((marginal_distribution[1]) && ((exhaustive) || (marginal_distribution[1]->nb_value <= DISPLAY_CONDITIONAL_NB_VALUE))) {
    cumul = new double*[marginal_distribution[0]->nb_value];
    for (i = marginal_distribution[0]->offset;i < marginal_distribution[0]->nb_value;i++) {
      if (marginal_distribution[0]->frequency[i] > 0) {
        cumul[i] = value_vec[i]->marginal_distribution[1]->cumul_computation();
      }
    }

    width[0] = column_width(marginal_distribution[1]->nb_value - 1);

    width[1] = 0;
    for (i = marginal_distribution[0]->offset;i < marginal_distribution[0]->nb_value;i++) {
      if (marginal_distribution[0]->frequency[i] > 0) {
        buff = column_width(value_vec[i]->marginal_distribution[1]->max);
        if (buff > width[1]) {
          width[1] = buff;
        }
      }
    }
    width[1] += ASCII_SPACE;

    width[2] = 0;
    for (i = marginal_distribution[0]->offset;i < marginal_distribution[0]->nb_value;i++) {
      if (marginal_distribution[0]->frequency[i] > 0) {
        buff = column_width(value_vec[i]->marginal_distribution[1]->nb_value - value_vec[i]->marginal_distribution[1]->offset ,
                            cumul[i] + value_vec[i]->marginal_distribution[1]->offset);
        if (buff > width[2]) {
          width[2] = buff;
        }
      }
    }
    width[2] += ASCII_SPACE;

    // ecriture des lois empiriques et des fonctions de repartition

    os << "\n  ";
    for (i = marginal_distribution[0]->offset;i < marginal_distribution[0]->nb_value;i++) {
      if (marginal_distribution[0]->frequency[i] > 0) {
        os << " | " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " " << i;
      }
    }
    for (i = marginal_distribution[0]->offset;i < marginal_distribution[0]->nb_value;i++) {
      if (marginal_distribution[0]->frequency[i] > 0) {
        os << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_DISTRIBUTION] << " "
           << i << " " << STAT_label[STATL_FUNCTION];
      }
    }
    os << endl;

    for (i = 0;i < marginal_distribution[1]->nb_value;i++) {
      os << setw(width[0]) << i;
      for (j = marginal_distribution[0]->offset;j < marginal_distribution[0]->nb_value;j++) {
        if (marginal_distribution[0]->frequency[j] > 0) {
          if (i < value_vec[j]->marginal_distribution[1]->nb_value) {
            os << setw(width[1]) << value_vec[j]->marginal_distribution[1]->frequency[i];
          }
          else {
            os << setw(width[1]) << " ";
          }
        }
      }
      for (j = marginal_distribution[0]->offset;j < marginal_distribution[0]->nb_value;j++) {
        if (marginal_distribution[0]->frequency[j] > 0) {
          if (i < value_vec[j]->marginal_distribution[1]->nb_value) {
            os << setw(width[2]) << cumul[j][i];
          }
          else {
            os << setw(width[2]) << " ";
          }
        }
      }
      os << endl;
    }

    for (i = marginal_distribution[0]->offset;i < marginal_distribution[0]->nb_value;i++) {
      if (marginal_distribution[0]->frequency[i] > 0) {
        delete [] cumul[i];
      }
    }
    delete [] cumul;
  }

  switch (type) {

  case ORDINAL : {
    int nb_histo;
    const FrequencyDistribution **value_marginal;


    value_marginal = new const FrequencyDistribution*[marginal_distribution[0]->nb_value];
    nb_histo = 0;
    for (i = marginal_distribution[0]->offset;i < marginal_distribution[0]->nb_value;i++) {
      if (marginal_distribution[0]->frequency[i] > 0) {
        value_marginal[nb_histo++] = value_vec[i]->marginal_distribution[1];
      }
    }

    test = value_marginal[0]->kruskal_wallis_test(nb_histo - 1 , value_marginal + 1);

    os << "\n" << STAT_label[STATL_KRUSKAL_WALLIS_TEST] << endl;
    test->ascii_print(os);

    delete [] value_marginal;
    break;
  }

  case NUMERIC : {

    // tableau d'analyse de variance

    square_sum[0] = 0.;
    square_sum[1] = 0.;
    df[0] = -1;

    for (i = marginal_distribution[0]->offset;i < marginal_distribution[0]->nb_value;i++) {
      if (marginal_distribution[0]->frequency[i] > 0) {
        diff = value_vec[i]->mean[1] - mean[1];
        square_sum[0] += diff * diff * marginal_distribution[0]->frequency[i];
        square_sum[1] += value_vec[i]->covariance[1][1] * (marginal_distribution[0]->frequency[i] - 1);
        df[0]++;
      }
    }

    square_sum[2] = covariance[1][1] * (nb_vector - 1);

    df[1] = nb_vector - 1 - df[0];
    df[2] = nb_vector - 1;

    for (i = 0;i < 3;i++) {
      mean_square[i] = square_sum[i] / df[i];
    }

#   ifdef DEBUG
    os << "\ntest: " << square_sum[0] + square_sum[1] << " | " << square_sum[2] << endl;
#   endif

    width[0] = column_width(nb_vector - 1) + ASCII_SPACE;
    width[1] = column_width(3 , square_sum) + ASCII_SPACE;
    width[2] = column_width(3 , mean_square) + ASCII_SPACE;

    os << "\n" << STAT_label[STATL_VARIATION_SOURCE] << " | " << STAT_label[STATL_FREEDOM_DEGREES]
       << " | " << STAT_label[STATL_SQUARE_SUM] << " | " << STAT_label[STATL_MEAN_SQUARE] << endl;
    for (i = 0;i < 3;i++) {
      switch (i) {
      case 0 :
        os << STAT_label[STATL_BETWEEN_SAMPLES];
        break;
      case 1 :
        os << STAT_label[STATL_WITHIN_SAMPLES] << " ";
        break;
      case 2 :
        os << STAT_label[STATL_TOTAL] << "          ";
        break;
      }

      os << setw(width[0]) << df[i];
      os << setw(width[1]) << square_sum[i];
      os << setw(width[2]) << mean_square[i] << endl;
    }
    os << endl;

    test = new Test(FISHER , true , df[0] , df[1] , mean_square[0] / mean_square[1]);
    test->F_critical_probability_computation();

    test->ascii_print(os , false , (df[0] == 1 ? false : true));
  }
  }

  delete test;

  os.setf((FMTFLAGS)old_adjust , ios::adjustfield);

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture des resultats d'une analyse de variance a un facteur dans un fichier.
 *
 *  arguments : reference sur un objet StatError, path,
 *              type de la variable reponse (ORDINAL/NUMERIC), pointeurs sur les sous-echantillons
 *              pour chaque niveau possible, flag niveau de detail.
 *
 *--------------------------------------------------------------*/

bool Vectors::variance_analysis_ascii_write(StatError &error , const char *path , int response_type ,
                                            const Vectors **value_vec , bool exhaustive) const

{
  bool status;
  ofstream out_file(path);


  error.init();

  if (!out_file) {
    status = false;
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {
    status = true;
    variance_analysis_ascii_write(out_file , response_type , value_vec , exhaustive);
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture des resultats d'une analyse de variance a un facteur
 *  dans un fichier au format tableur.
 *
 *  arguments : reference sur un objet StatError, path,
 *              type de la variable reponse (ORDINAL/NUMERIC),
 *              pointeurs sur les sous-echantillons pour chaque niveau possible.
 *
 *--------------------------------------------------------------*/

bool Vectors::variance_analysis_spreadsheet_write(StatError &error , const char *path ,
                                                  int response_type , const Vectors **value_vec) const

{
  bool status;
  register int i , j;
  int df[3];
  double diff , square_sum[3] , mean_square[3] , **cumul;
  Test *test;
  ofstream out_file(path);


  error.init();

  if (!out_file) {
    status = false;
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {
    status = true;

    // lois conditionnelles

    out_file << STAT_label[STATL_VALUE];
    for (i = marginal_distribution[0]->offset;i < marginal_distribution[0]->nb_value;i++) {
      if (marginal_distribution[0]->frequency[i] > 0) {
        out_file << "\t"  << i;
      }
    }
    out_file << endl;

    out_file << STAT_label[STATL_SAMPLE_SIZE];
    for (i = marginal_distribution[0]->offset;i < marginal_distribution[0]->nb_value;i++) {
      if (marginal_distribution[0]->frequency[i] > 0) {
        out_file << "\t" << marginal_distribution[0]->frequency[i];
      }
    }
    out_file << endl;

    out_file << STAT_label[STATL_MEAN];
    for (i = marginal_distribution[0]->offset;i < marginal_distribution[0]->nb_value;i++) {
      if (marginal_distribution[0]->frequency[i] > 0) {
        out_file << "\t" << value_vec[i]->mean[1];
      }
    }
    out_file << endl;

    out_file << STAT_label[STATL_VARIANCE];
    for (i = marginal_distribution[0]->offset;i < marginal_distribution[0]->nb_value;i++) {
      if (marginal_distribution[0]->frequency[i] > 0) {
        out_file << "\t" << value_vec[i]->covariance[1][1];
      }
    }
    out_file << endl;

    out_file << STAT_label[STATL_STANDARD_DEVIATION];
    for (i = marginal_distribution[0]->offset;i < marginal_distribution[0]->nb_value;i++) {
      if (marginal_distribution[0]->frequency[i] > 0) {
        out_file << "\t" << sqrt(value_vec[i]->covariance[1][1]);
      }
    }
    out_file << endl;

    out_file << STAT_label[STATL_MEAN_ABSOLUTE_DEVIATION];
    for (i = marginal_distribution[0]->offset;i < marginal_distribution[0]->nb_value;i++) {
      if (marginal_distribution[0]->frequency[i] > 0) {
        out_file << "\t" << value_vec[i]->mean_absolute_deviation_computation(1);
      }
    }
    out_file << endl;

    out_file << STAT_label[STATL_CONCENTRATION_COEFF];
    for (i = marginal_distribution[0]->offset;i < marginal_distribution[0]->nb_value;i++) {
      if (marginal_distribution[0]->frequency[i] > 0) {
        out_file << "\t" << (value_vec[i]->nb_vector > 1 ? (value_vec[i]->mean_absolute_difference_computation(1) /
                              (2 * value_vec[i]->mean[1])) * ((double)(value_vec[i]->nb_vector - 1) / (double)value_vec[i]->nb_vector) : 1.);
      }
    }
    out_file << endl;

    out_file << STAT_label[STATL_SKEWNESS_COEFF];
    for (i = marginal_distribution[0]->offset;i < marginal_distribution[0]->nb_value;i++) {
      if (marginal_distribution[0]->frequency[i] > 0) {
        out_file << "\t" << value_vec[i]->skewness_computation(1);
      }
    }
    out_file << endl;

    out_file << STAT_label[STATL_KURTOSIS_COEFF];
    for (i = marginal_distribution[0]->offset;i < marginal_distribution[0]->nb_value;i++) {
      if (marginal_distribution[0]->frequency[i] > 0) {
        out_file << "\t" << value_vec[i]->kurtosis_computation(1);
      }
    }
    out_file << endl;

    if (marginal_distribution[1]) {
      cumul = new double*[marginal_distribution[0]->nb_value];
      for (i = marginal_distribution[0]->offset;i < marginal_distribution[0]->nb_value;i++) {
        if (marginal_distribution[0]->frequency[i] > 0) {
          cumul[i] = value_vec[i]->marginal_distribution[1]->cumul_computation();
        }
      }

      out_file << "\n";
      for (i = marginal_distribution[0]->offset;i < marginal_distribution[0]->nb_value;i++) {
        if (marginal_distribution[0]->frequency[i] > 0) {
          out_file << "\t" << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " " << i;
        }
      }
      for (i = marginal_distribution[0]->offset;i < marginal_distribution[0]->nb_value;i++) {
        if (marginal_distribution[0]->frequency[i] > 0) {
          out_file << "\t" << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_DISTRIBUTION] << " "
                   << i << " " << STAT_label[STATL_FUNCTION];
        }
      }
      out_file << endl;

      for (i = 0;i < marginal_distribution[1]->nb_value;i++) {
        out_file << i;
        for (j = marginal_distribution[0]->offset;j < marginal_distribution[0]->nb_value;j++) {
          if (marginal_distribution[0]->frequency[j] > 0) {
            out_file << "\t";
            if (i < value_vec[j]->marginal_distribution[1]->nb_value) {
              out_file << value_vec[j]->marginal_distribution[1]->frequency[i];
            }
          }
        }
        for (j = marginal_distribution[0]->offset;j < marginal_distribution[0]->nb_value;j++) {
          if (marginal_distribution[0]->frequency[j] > 0) {
            out_file << "\t";
            if (i < value_vec[j]->marginal_distribution[1]->nb_value) {
              out_file << cumul[j][i];
            }
          }
        }
        out_file << endl;
      }

      for (i = marginal_distribution[0]->offset;i < marginal_distribution[0]->nb_value;i++) {
        if (marginal_distribution[0]->frequency[i] > 0) {
          delete [] cumul[i];
        }
      }
      delete [] cumul;
    }

    switch (response_type) {

    case ORDINAL : {
      int nb_histo;
      const FrequencyDistribution **value_marginal;


      value_marginal = new const FrequencyDistribution*[marginal_distribution[0]->nb_value];
      nb_histo = 0;
      for (i = marginal_distribution[0]->offset;i < marginal_distribution[0]->nb_value;i++) {
        if (marginal_distribution[0]->frequency[i] > 0) {
          value_marginal[nb_histo++] = value_vec[i]->marginal_distribution[1];
        }
      }

      test = value_marginal[0]->kruskal_wallis_test(nb_histo - 1 , value_marginal + 1);

      out_file << "\n" << STAT_label[STATL_KRUSKAL_WALLIS_TEST] << endl;
      test->spreadsheet_print(out_file);

      delete [] value_marginal;
      break;
    }

    case NUMERIC : {

      // tableau d'analyse de variance

      square_sum[0] = 0.;
      square_sum[1] = 0.;
      df[0] = -1;

      for (i = marginal_distribution[0]->offset;i < marginal_distribution[0]->nb_value;i++) {
        if (marginal_distribution[0]->frequency[i] > 0) {
          diff = value_vec[i]->mean[1] - mean[1];
          square_sum[0] += diff * diff * marginal_distribution[0]->frequency[i];
          square_sum[1] += value_vec[i]->covariance[1][1] * (marginal_distribution[0]->frequency[i] - 1);
          df[0]++;
        }
      }

      square_sum[2] = covariance[1][1] * (nb_vector - 1);

      df[1] = nb_vector - 1 - df[0];
      df[2] = nb_vector - 1;

      for (i = 0;i < 3;i++) {
        mean_square[i] = square_sum[i] / df[i];
      }

      out_file << "\n" << STAT_label[STATL_VARIATION_SOURCE] << "\t" << STAT_label[STATL_FREEDOM_DEGREES]
               << "\t" << STAT_label[STATL_SQUARE_SUM] << "\t" << STAT_label[STATL_MEAN_SQUARE] << endl;
      for (i = 0;i < 3;i++) {
        switch (i) {
        case 0 :
          out_file << STAT_label[STATL_BETWEEN_SAMPLES];
          break;
        case 1 :
          out_file << STAT_label[STATL_WITHIN_SAMPLES];
          break;
        case 2 :
          out_file << STAT_label[STATL_TOTAL];
          break;
        }

        out_file << "\t" << df[i] << "\t" << square_sum[i] << "\t" << mean_square[i] << endl;
      }
      out_file << endl;

      test = new Test(FISHER , true , df[0] , df[1] , mean_square[0] / mean_square[1]);
      test->F_critical_probability_computation();

      test->spreadsheet_print(out_file , (df[0] == 1 ? false : true));
      break;
    }
    }

    delete test;
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Analyse de variance a un facteur.
 *
 *  arguments : reference sur un objet StatError, stream, indices des 2 variables,
 *              type de la variable reponse (ORDINAL/NUMERIC), path,
 *              format de fichier ('a' : ASCII, 's' : Spreadsheet).
 *
 *--------------------------------------------------------------*/

bool Vectors::variance_analysis(StatError &error , ostream &os , int class_variable ,
                                int response_variable , int response_type ,
                                const char *path , char format) const

{
  bool status = true;
  register int i;
  int *value_nb_vector , **index;
  Vectors *vec;
  const Vectors **value_vec;


  error.init();

  if (class_variable == response_variable) {
    status = false;
    error.update(STAT_error[STATR_VARIABLE_INDICES]);
  }

  if ((class_variable < 1) || (class_variable > nb_variable)) {
    status = false;
    ostringstream error_message;
    error_message << STAT_label[STATL_VARIABLE] << " " << class_variable << ": "
                  << STAT_error[STATR_VARIABLE_INDEX];
    error.update((error_message.str()).c_str());
  }

  else {
    class_variable--;

    if (type[class_variable] != INT_VALUE) {
      status = false;
      ostringstream error_message;
      error_message << STAT_label[STATL_VARIABLE] << " " << class_variable + 1 << ": "
                    << STAT_error[STATR_VARIABLE_TYPE];
      error.correction_update((error_message.str()).c_str() , STAT_variable_word[INT_VALUE]);
    }

    else if ((!marginal_distribution[class_variable]) ||
             (marginal_distribution[class_variable]->nb_value - marginal_distribution[class_variable]->offset > VARIANCE_ANALYSIS_NB_VALUE)) {
      status = false;
      ostringstream error_message;
      error_message << STAT_label[STATL_VARIABLE] << " " << class_variable + 1 << ": "
                    << STAT_error[STATR_NB_VALUE];
      error.update((error_message.str()).c_str());
    }
  }

  if ((response_variable < 1) || (response_variable > nb_variable)) {
    status = false;
    ostringstream error_message;
    error_message << STAT_label[STATL_VARIABLE] << " " << response_variable << ": "
                  << STAT_error[STATR_VARIABLE_INDEX];
    error.update((error_message.str()).c_str());
  }

  else {
    response_variable--;

    if (response_type == ORDINAL) {
      if (type[response_variable] != INT_VALUE) {
        status = false;
        ostringstream error_message;
        error_message << STAT_label[STATL_VARIABLE] << " " << response_variable + 1 << ": "
                      << STAT_error[STATR_VARIABLE_TYPE];
        error.correction_update((error_message.str()).c_str() , STAT_variable_word[INT_VALUE]);
      }

      else if (!marginal_distribution[response_variable]) {
        ostringstream error_message;
        error_message << STAT_label[STATL_VARIABLE] << " " << response_variable + 1 << " "
                      << STAT_error[STATR_SHIFTED_SCALED];
        error.update((error_message.str()).c_str());
      }
    }
  }

  if (status) {

    // extraction des sous-echantillons pour chaque niveau possible

    vec = select_variable(class_variable , response_variable);

    value_nb_vector = new int[vec->marginal_distribution[0]->nb_value];
    index = new int*[vec->marginal_distribution[0]->nb_value];
    for (i = vec->marginal_distribution[0]->offset;i < vec->marginal_distribution[0]->nb_value;i++) {
      if (vec->marginal_distribution[0]->frequency[i] > 0) {
        value_nb_vector[i] = 0;
        index[i] = new int[vec->marginal_distribution[0]->frequency[i]];
      }
    }

    for (i = 0;i < nb_vector;i++) {
      index[vec->int_vector[i][0]][(value_nb_vector[vec->int_vector[i][0]])++] = i;
    }

    value_vec = new const Vectors*[vec->marginal_distribution[0]->nb_value];
    for (i = vec->marginal_distribution[0]->offset;i < vec->marginal_distribution[0]->nb_value;i++) {
      if (vec->marginal_distribution[0]->frequency[i] > 0) {
        value_vec[i] = new Vectors(*vec , value_nb_vector[i] , index[i]);
      }
    }
    delete [] value_nb_vector;

#   ifdef MESSAGE
    vec->variance_analysis_ascii_write(os , response_type , value_vec, false);
#   endif

    if (path) {
      switch (format) {
      case 'a' :
        status = vec->variance_analysis_ascii_write(error , path , response_type ,
                                                    value_vec , true);
        break;
      case 's' :
        status = vec->variance_analysis_spreadsheet_write(error , path , response_type ,
                                                          value_vec);
        break;
      }
    }

    for (i = vec->marginal_distribution[0]->offset;i < vec->marginal_distribution[0]->nb_value;i++) {
      if (vec->marginal_distribution[0]->frequency[i] > 0) {
        delete [] index[i];
        delete value_vec[i];
      }
    }
    delete [] index;
    delete [] value_vec;

    delete vec;
  }

  return status;
}


};  // namespace stat_tool
