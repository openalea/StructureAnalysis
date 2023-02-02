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

#include <string>
#include <sstream>
#include <iomanip>

#include "vectors.h"
#include "stat_label.h"

using namespace std;


namespace stat_tool {



/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the Spearman rank correlation coefficient matrix.
 *
 *  \return correlation matrix.
 */
/*--------------------------------------------------------------*/

double** Vectors::spearman_rank_correlation_computation() const

{
  int i , j , k;
  int *pfrequency;
  double main_term , rank_mean , rank_diff , *correction , **correlation = NULL , **rank;


  if (nb_vector > 2) {
    correlation = new double*[nb_variable];
    for (i = 0;i < nb_variable;i++) {
      correlation[i] = new double[nb_variable];
    }

    // computation of the main term and the correction terms for ties

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

    // rank computation

    rank_mean = (double)(nb_vector + 1) / 2.;

    rank = new double*[nb_variable];

    for (i = 0;i < nb_variable;i++) {
      rank[i] = marginal_distribution[i]->rank_computation();
    }

    for (i = 0;i < nb_variable;i++) {
      correlation[i][i] = 1.;
      for (j = i + 1;j < nb_variable;j++) {

        // computation of the differences between centered ranks

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


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the Kendall rank correlation coefficient matrix.
 *
 *  \return correlation matrix.
 */
/*--------------------------------------------------------------*/

double** Vectors::kendall_rank_correlation_computation() const

{
  int i , j , k , m , n , p;
  int diff , *current_frequency , *pfrequency , *cumul_frequency , *index , **frequency;
  double nb_pair , sum , rank_diff_sign , *correction , **correlation = NULL;


  if (nb_vector > 2) {
    correlation = new double*[nb_variable];
    for (i = 0;i < nb_variable;i++) {
      correlation[i] = new double[nb_variable];
    }

    // computation of the main term and the correction terms for ties

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

      // computation of an order on the individuals on the basis of the values taken by the 1st variable

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

          // computation of the frequencies corresponding to the joint distribution of the 2 variables

          frequency = joint_frequency_computation(i , j);

          // computation of concordant and discordant pairs

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

          // computation of the signs of the rank differences

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


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a rank correlation coefficient matrix.
 *
 *  \param[in,out] os          stream,
 *  \param[in]     correl_type rank correlation coefficient type (SPEARMAN/KENDALL),
 *  \param[in]     correlation pointer on the rank correlation coefficients.
 */
/*--------------------------------------------------------------*/

ostream& Vectors::rank_correlation_ascii_write(ostream &os , correlation_type correl_type ,
                                               double **correlation) const

{
  int i , j;
  int buff , width[2];
  Test *test;
  ios_base::fmtflags format_flags;


  format_flags = os.setf(ios::right , ios::adjustfield);

  width[0] = column_width(nb_variable);

  // computation of the column width

  width[1] = 0;
  for (i = 0;i < nb_variable;i++) {
    buff = column_width(nb_variable , correlation[i]);
    if (buff > width[1]) {
      width[1] = buff;
    }
  }
  width[1] += ASCII_SPACE;

  // writing of the rank correlation coefficient matrix

  switch (correl_type) {
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

  // test of the significant character of the rank correlation coefficients

  if (nb_vector > 2) {
    switch (correl_type) {
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

      switch (correl_type) {
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

  os.setf(format_flags , ios::adjustfield);

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a rank correlation coefficient matrix in a file.
 *
 *  \param[in] error       reference on a StatError object,
 *  \param[in] path        file path,
 *  \param[in] correl_type rank correlation coefficient type (SPEARMAN/KENDALL),
 *  \param[in] correlation pointer on the rank correlation coefficients.
 *
 *  \return                error status.
 */
/*--------------------------------------------------------------*/

bool Vectors::rank_correlation_ascii_write(StatError &error , const string path ,
                                           correlation_type correl_type , double **correlation) const

{
  bool status;
  ofstream out_file(path.c_str());


  error.init();

  if (!out_file) {
    status = false;
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {
    status = true;
    rank_correlation_ascii_write(out_file , correl_type , correlation);
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of a rank correlation coefficient matrix
 *         (either in the Spearman or in the Kendall sense).
 *
 *  \param[in] error       reference on a StatError object,
 *  \param[in] display     flag for displaying the rank correlation coefficient matrix,
 *  \param[in] correl_type rank correlation coefficient type (SPEARMAN/KENDALL),
 *  \param[in] path        file path.
 *
 *  \return                error status.
 */
/*--------------------------------------------------------------*/

bool Vectors::rank_correlation_computation(StatError &error , bool display ,
                                           correlation_type correl_type , const string path) const

{
  bool status = true;
  int i;
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
    switch (correl_type) {
    case SPEARMAN :
      correlation = spearman_rank_correlation_computation();
      break;
    case KENDALL :
      correlation = kendall_rank_correlation_computation();
      break;
    }

    if (display) {
      rank_correlation_ascii_write(cout , correl_type , correlation);
    }
    if (!path.empty()) {
      status = rank_correlation_ascii_write(error , path , correl_type , correlation);
    }

    for (i = 0;i < nb_variable;i++) {
      delete [] correlation[i];
    }
    delete [] correlation;
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the Spearman rank correlation coefficient between 2 variables.
 *
 *  \return correlation coefficient.
 */
/*--------------------------------------------------------------*/

double Vectors::spearman_rank_single_correlation_computation() const

{
  int i , j;
  int *pfrequency;
  double correlation , main_term , rank_mean , correction[2] , *rank[2];


  // computation of the main term and the correction terms for ties

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

  // rank computation

  rank_mean = (double)(nb_vector + 1) / 2.;

  for (i = 0;i < 2;i++) {
    rank[i] = marginal_distribution[i]->rank_computation();
  }

  // computation of the differences between centered ranks

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


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the Kendall rank correlation coefficient between 2 variables.
 *
 *  \return correlation coefficient.
 */
/*--------------------------------------------------------------*/

double Vectors::kendall_rank_single_correlation_computation() const

{
  int i , j , k , m;
  int diff , *current_frequency , *pfrequency , *cumul_frequency , *index ,
      **frequency;
  double sum , correlation , nb_pair , correction[2];


  // computation of the main term and the correction terms for ties

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

    // computation of the frequencies corresponding to the joint distribution of the 2 variables

    frequency = joint_frequency_computation(0 , 1);

    // computation of concordant and discordant pairs

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

    // computation of an order on the individuals on the basis of the values taken by the 1st variable

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

    // computation of the signs of the rank differences

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


/*--------------------------------------------------------------*/
/**
 *  \brief Comparison of vectors (computation of the matrix of pairwise distances between vectors).
 *
 *  \param[in] error           reference on a StatError object,
 *  \param[in] ivector_dist    reference on a VectorDistance object,
 *  \param[in] standardization flag standardization (only for variables of the same type).
 *
 *  \return                    DistanceMatrix object.
 */
/*--------------------------------------------------------------*/

DistanceMatrix* Vectors::comparison(StatError &error , const VectorDistance &ivector_dist ,
                                    bool standardization) const
{
  bool status = true;
  int i , j , k;
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
      if (ivector_dist.var_type[i] != NUMERIC) {
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

        if ((ivector_dist.var_type[i] == NOMINAL) &&
            ((min_value[i] < 0) || (max_value[i] >= NB_CATEGORY) ||
             ((ivector_dist.category_distance[i]) && (ivector_dist.nb_value[i] != max_value[i] + 1)))) {
          status = false;
          ostringstream error_message;
          error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                        << STAT_error[STATR_NB_CATEGORY];
          error.update((error_message.str()).c_str());
        }

        if ((ivector_dist.var_type[i] == CIRCULAR) &&
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
        if (ivector_dist.var_type[i] != ivector_dist.var_type[0]) {
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

    // rank computation for ordinal variables

    rank = new double*[nb_variable];

    if (standardization) {
      for (i = 0;i < nb_variable;i++) {
        if (vector_dist->var_type[i] == ORDINAL) {
          rank[i] = marginal_distribution[i]->rank_computation();
        }
        else {
          rank[i] = NULL;
        }
      }
    }

    else {
      if (vector_dist->var_type[0] == ORDINAL) {
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

    }

    // computation of dispersion measures for standardization of variables

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
        if (vector_dist->var_type[i] == NUMERIC) {
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
          switch (vector_dist->var_type[k]) {

          case NOMINAL : {
            if (!(vector_dist->category_distance[k])) {
              if (int_vector[i][k] == int_vector[j][k]) {
                ldistance = 0.;
              }
              else {
                ldistance = 1.;
              }
            }

            else {
              ldistance = vector_dist->category_distance[k][int_vector[i][k]][int_vector[j][k]];
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

          if (standardization) {

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
          }

          else {

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


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the frequencies corresponding to the joint distribution of
 *         2 variables.
 *
 *  \param[in] variable1 variable 1 index,
 *  \param[in] variable2 variable 2 index.
 *
 *  \return              joint frequency distribution.
 */
/*--------------------------------------------------------------*/

int** Vectors::joint_frequency_computation(int variable1 , int variable2) const

{
  int i , j;
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


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a contingency table.
 *
 *  \param[in,out] os                stream,
 *  \param[in]     variable1         variable 1 index,
 *  \param[in]     variable2         variable 2 index,
 *  \param[in]     frequency         pointer on the contingency table,
 *  \param[in]     deviation         pointer on the deviations to the independence assumption,
 *  \param[in]     chi2_contribution pointer on the contributions to the chi2 value,
 *  \param[in]     test              reference on a Test object,
 *  \param[in]     file_flag         file flag.
 */
/*--------------------------------------------------------------*/

ostream& Vectors::contingency_table_ascii_write(ostream &os , int variable1 , int variable2 ,
                                                int **frequency , double **deviation ,
                                                double **chi2_contribution , const Test &test ,
                                                bool file_flag) const

{
  int i , j;
  int buff , width[2];
  ios_base::fmtflags format_flags;


  if ((file_flag) || ((marginal_distribution[variable1]->nb_value - marginal_distribution[variable1]->offset <= DISPLAY_CONTINGENCY_NB_VALUE) &&
       (marginal_distribution[variable2]->nb_value - marginal_distribution[variable2]->offset <= DISPLAY_CONTINGENCY_NB_VALUE))) {
    format_flags = os.setf(ios::right , ios::adjustfield);

    // computation of the column widths

    width[0] = column_width(MAX(marginal_distribution[variable1]->nb_value , marginal_distribution[variable2]->nb_value) - 1);
    width[1] = column_width(nb_vector) + ASCII_SPACE;

    // writing of the contingency table

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

    // computation of the column width

    width[1] = 0;
    for (i = marginal_distribution[variable1]->offset;i < marginal_distribution[variable1]->nb_value;i++) {
      buff = column_width(marginal_distribution[variable2]->nb_value - marginal_distribution[variable2]->offset ,
                          deviation[i] + marginal_distribution[variable2]->offset);
      if (buff > width[1]) {
        width[1] = buff;
      }
    }
    width[1] += ASCII_SPACE;

    // writing of the deviation table

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

    // computation of the column width

    width[1] = 0;
    for (i = marginal_distribution[variable1]->offset;i < marginal_distribution[variable1]->nb_value;i++) {
      buff = column_width(marginal_distribution[variable2]->nb_value - marginal_distribution[variable2]->offset ,
                          chi2_contribution[i] + marginal_distribution[variable2]->offset);
      if (buff > width[1]) {
        width[1] = buff;
      }
    }
    width[1] += ASCII_SPACE;

    // writing of the table of contributions to the chi2 value

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

    os.setf(format_flags , ios::adjustfield);
  }

  os << "\n";
  test.ascii_print(os);

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a contingency table in a file.
 *
 *  \param[in] error             reference on a StatError object,
 *  \param[in] path              file path,
 *  \param[in] variable1         variable 1  index,
 *  \param[in] variable2         variable 2  index,
 *  \param[in] frequency         pointer on the contingency table,
 *  \param[in] deviation         pointer on the deviations to the independence assumption,
 *  \param[in] chi2_contribution pointer on the contributions to the chi2 value,
 *  \param[in] test              reference on a Test object.
 *
 *  \return                      error status.
 */
/*--------------------------------------------------------------*/

bool Vectors::contingency_table_ascii_write(StatError &error , const string path ,
                                            int variable1 , int variable2 , int **frequency ,
                                            double **deviation , double **chi2_contribution ,
                                            const Test &test) const

{
  bool status;
  ofstream out_file(path.c_str());


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


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a contingency table in a file at the spreadsheet format.
 *
 *  \param[in] error             reference on a StatError object,
 *  \param[in] path              file path,
 *  \param[in] variable1         variable 1 index,
 *  \param[in] variable2         variable 2 index,
 *  \param[in] frequency         pointer on the contingency table,
 *  \param[in] deviation         pointer on the deviations to the independence assumption,
 *  \param[in] chi2_contribution pointer on the contributions to the chi2 value,
 *  \param[in] test              reference on the result of a chi2 test.
 *
 *  \return                      error status.
 */
/*--------------------------------------------------------------*/

bool Vectors::contingency_table_spreadsheet_write(StatError &error , const string path ,
                                                  int variable1 , int variable2 , int **frequency ,
                                                  double **deviation , double **chi2_contribution ,
                                                  const Test &test) const

{
  bool status;
  int i , j;
  ofstream out_file(path.c_str());


  error.init();

  if (!out_file) {
    status = false;
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {
    status = true;

    // writing of the contingency table

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

    // writing of the deviation table

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

    // writing of the table of contributions to the chi2 value

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


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of a contingency table for 2 categorical variables.
 *
 *  \param[in] error     reference on a StatError object,
 *  \param[in] display   flag for displaying the contingency table,
 *  \param[in] variable1 variable 1 index,
 *  \param[in] variable2 variable 2 index,
 *  \param[in] path      file path,
 *  \param[in] format    file format (ASCII/SPREADSHEET).
 *
 *  \return              error status.
 */
/*--------------------------------------------------------------*/

bool Vectors::contingency_table(StatError &error , bool display , int variable1 ,
                                int variable2 , const string path , output_format format) const

{
  bool status = true;
  int i , j;
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

    // computation of the frequencies corresponding to the joint distribution of the 2 variables

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

    // computation of the chi2 value

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

    if (display) {
      contingency_table_ascii_write(cout , variable1 , variable2 , frequency ,
                                    deviation , chi2_contribution , *test);
    }

    if (!path.empty()) {
      switch (format) {
      case ASCII :
        status = contingency_table_ascii_write(error , path , variable1 , variable2 , frequency ,
                                               deviation , chi2_contribution , *test);
        break;
      case SPREADSHEET :
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


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of the results of a one-way analysis of variance.
 *
 *  \param[in,out] os         stream,
 *  \param[in]     type       response variable type (ORDINAL/NUMERIC),
 *  \param[in]     value_vec  pointer on the sub-samples for each level of the explanatory variable,
 *  \param[in]     exhaustive flag detail level.
 */
/*--------------------------------------------------------------*/

ostream& Vectors::variance_analysis_ascii_write(ostream &os , int type , const Vectors **value_vec ,
                                                bool exhaustive) const

{
  int i , j;
  int buff , width[4] , df[3];
  double diff , square_sum[3] , mean_square[3] , *value_mean , *variance , *standard_deviation ,
         *mean_absolute_deviation , *concentration_coeff , *skewness_coeff , *kurtosis_coeff , **cumul;
  Test *test;
  ios_base::fmtflags format_flags;


  format_flags = os.setf(ios::right , ios::adjustfield);

  // conditional distributions

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
      mean_absolute_deviation[i] = value_vec[i]->mean_absolute_deviation_computation(1 , value_vec[i]->mean[1]);
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

    // writing of the frequency distributions and the cumulative distribution functions

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

    // table of variance analysis

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

  os.setf(format_flags , ios::adjustfield);

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of the results of a one-way analysis of variance in a file.
 *
 *  \param[in] error         reference on a StatError object,
 *  \param[in] path          file path,
 *  \param[in] response_type response variable type (ORDINAL/NUMERIC),
 *  \param[in] value_vec     pointer on the sub-samples for each level of the explanatory variable,
 *  \param[in] exhaustive    flag detail level.
 *
 *  \return                  error status.
 */
/*--------------------------------------------------------------*/

bool Vectors::variance_analysis_ascii_write(StatError &error , const string path , int response_type ,
                                            const Vectors **value_vec , bool exhaustive) const

{
  bool status;
  ofstream out_file(path.c_str());


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


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of the results of a one-way variance analysis in a file
 *         at the spreadsheet format.
 *
 *  \param[in] error         reference on a StatError object,
 *  \param[in] path          file path,
 *  \param[in] response_type response variable type (ORDINAL/NUMERIC),
 *  \param[in] value_vec     pointer on the sub-samples for each level of the explanatory variable.
 *
 *  \return                  error status.
 */
/*--------------------------------------------------------------*/

bool Vectors::variance_analysis_spreadsheet_write(StatError &error , const string path ,
                                                  int response_type , const Vectors **value_vec) const

{
  bool status;
  int i , j;
  int df[3];
  double diff , square_sum[3] , mean_square[3] , **cumul;
  Test *test;
  ofstream out_file(path.c_str());


  error.init();

  if (!out_file) {
    status = false;
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {
    status = true;

    // conditional distributions

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
        out_file << "\t" << value_vec[i]->mean_absolute_deviation_computation(1 , value_vec[i]->mean[1]);
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

      // table of variance analysis

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


/*--------------------------------------------------------------*/
/**
 *  \brief One-way analysis of variance.
 *
 *  \param[in] error             reference on a StatError object,
 *  \param[in] display           flag for displaying the ANOVA results,
 *  \param[in] class_variable    explanatory variable index,
 *  \param[in] response_variable response variable index,
 *  \param[in] response_type     response variable type (ORDINAL/NUMERIC),
 *  \param[in] path              file path,
 *  \param[in] format            file format (ASCII/SPREADSHEET).
 *
 *  \return                      error status.
 */
/*--------------------------------------------------------------*/

bool Vectors::variance_analysis(StatError &error , bool display , int class_variable ,
                                int response_variable , int response_type ,
                                const string path , output_format format) const

{
  bool status = true;
  int i;
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

    // extraction of the sub-samples for each level of the explanatory variable

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

    if (display) {
      vec->variance_analysis_ascii_write(cout , response_type , value_vec , false);
    }

    if (!path.empty()) {
      switch (format) {
      case ASCII :
        status = vec->variance_analysis_ascii_write(error , path , response_type ,
                                                    value_vec , true);
        break;
      case SPREADSHEET :
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


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of sup norm distance between two empirical continuous distributions.
 *
 *  \param[in] error   reference on a StatError object,
 *  \param[in] display flag for displaying the sup norm distance,
 *  \param[in] ivec    reference on a Vector object.
 *
 *  \return           error status.
 */
/*--------------------------------------------------------------*/

bool Vectors::sup_norm_distance(StatError &error , bool display , const Vectors &ivec) const

{
  bool status = true , **selected_value;
  int i , j;
  int nb_crossing , previous_sign , sign , int_value , sample_size[2] , **rank;
  double min , max , sup_norm , diff_cumul , previous_diff_cumul , overlap , real_value ,
         sup_value , *cumul[2];


  if (nb_vector < SUP_NORM_DISTANCE_NB_VECTOR) {
    status = false;
    error.update(STAT_error[STATR_NB_VECTOR]);
  }
  if (nb_variable != 1) {
    status = false;
    error.correction_update(STAT_error[STATR_NB_VARIABLE] , 1);
  }

  if (ivec.nb_vector < SUP_NORM_DISTANCE_NB_VECTOR) {
    status = false;
    error.update(STAT_error[STATR_NB_VECTOR]);
  }
  if (ivec.nb_variable != 1) {
    status = false;
    error.correction_update(STAT_error[STATR_NB_VARIABLE] , 1);
  }

  if (type[0] != ivec.type[0]) {
    status = false;
    error.update(STAT_error[STATR_VARIABLE_TYPE]);
  }

  if (status) {
    selected_value = new bool*[2];
    selected_value[0] = new bool[nb_vector];
    selected_value[1] = new bool[ivec.nb_vector];

    for (i = 0;i < nb_vector;i++) {
      selected_value[0][i] = false;
    }
    for (i = 0;i < ivec.nb_vector;i++) {
      selected_value[1][i] = false;
    }

    rank = new int*[nb_vector + ivec.nb_vector];
    for (i = 0;i < nb_vector + ivec.nb_vector;i++) {
      rank[i] = new int[2];
    }

    switch (type[0]) {

    case INT_VALUE : {
      if ((marginal_distribution[0]) && (ivec.marginal_distribution[0])) {
        cumul[0] = marginal_distribution[0]->cumul_computation();
        cumul[1] = ivec.marginal_distribution[0]->cumul_computation();

        sup_norm = 0.;
        if (marginal_distribution[0]->offset < ivec.marginal_distribution[0]->offset) {
          max = cumul[0][ivec.marginal_distribution[0]->offset - 1];
          sup_value = ivec.marginal_distribution[0]->offset - 1;
          sign = 1;
        }
        else if (marginal_distribution[0]->offset > ivec.marginal_distribution[0]->offset) {
          max = cumul[1][marginal_distribution[0]->offset - 1];
          sup_value = marginal_distribution[0]->offset - 1;
          sign = -1;
        }
        else {
          max = 0.;
          sign = 0;
        }
        nb_crossing = 0;
        overlap = 0.;

#       ifdef DEBUG
        cout << "\n";
#       endif

        for (i = MAX(marginal_distribution[0]->offset , ivec.marginal_distribution[0]->offset);i < MIN(marginal_distribution[0]->nb_value , ivec.marginal_distribution[0]->nb_value);i++) {
          diff_cumul = cumul[0][i] - cumul[1][i];
          if (fabs(diff_cumul) > max) {
            max = fabs(diff_cumul);
            sup_value = i;
          }

          // test of non-crossing of empirical cumulative distribution functions

          previous_sign = sign;
          if (diff_cumul < 0.) {
            sign = -1;
          }
          else if (diff_cumul == 0.) {
            sign = 0;
          }
          else {
            sign = 1;
          }

          if ((previous_sign != 0) && (sign != 0) && (sign != previous_sign)) {
            if (display) {
              cout << "\n" << STAT_label[STATL_SUP_NORM_DISTANCE] << ": " << max
                   << "   sup norm value: " << sup_value << "   crossing value: " << i;
            }

            nb_crossing++;
            sup_norm += max;
            max = 0.;
          }

          overlap += MIN((double)marginal_distribution[0]->frequency[i] / (double)nb_vector ,
                         (double)ivec.marginal_distribution[0]->frequency[i] / (double)ivec.nb_vector);

#         ifdef DEBUG
          cout << i << "  ";
          if ((double)marginal_distribution[0]->frequency[i] / (double)nb_vector < (double)ivec.marginal_distribution[0]->frequency[i] / (double)ivec.nb_vector)) {
            cout << marginal_distribution[0]->frequency[i];
          }
          else {
            cout << ivec.marginal_distribution[0]->frequency[i];
          }
          cout << "   " << MIN((double)marginal_distribution[0]->frequency[i] / (double)nb_vector ,
                               (double)ivec.marginal_distribution[0]->frequency[i] / (double)ivec.nb_vector) << "  " << overlap << endl;
#         endif

        }
        sup_norm += max;

#       ifdef DEBUG
        cout << endl;
#       endif

        delete [] cumul[0];
        delete [] cumul[1];

        if (display) {
          cout << "\n" << STAT_label[STATL_SUP_NORM_DISTANCE] << ": " << sup_norm << "   "
               << STAT_label[STATL_OVERLAP] << ": " << overlap
               << "   number of crossings of cdfs: " << nb_crossing << "   sup norm value: " << sup_value
               << "   (" << STAT_label[STATL_MARGINAL] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << ")";
        }
      }

      // computation of ranks on the basis of the two samples

      i = 0;

      do {
        min = MAX(max_value[0], ivec.max_value[0]) + 1.;

        for (j = 0;j < nb_vector;j++) {
          if ((!selected_value[0][j]) && (int_vector[j][0] < min)) {
            min = int_vector[j][0];
            rank[i][0] = 0;
            rank[i][1] = j;
          }
        }
        for (j = 0;j < ivec.nb_vector;j++) {
          if ((!selected_value[1][j]) && (ivec.int_vector[j][0] < min)) {
            min = ivec.int_vector[j][0];
            rank[i][0] = 1;
            rank[i][1] = j;
          }
        }

        selected_value[rank[i][0]][rank[i][1]] = true;
        i++;

        for (j = 0;j < nb_vector;j++) {
          if ((!selected_value[0][j]) && (int_vector[j][0] == min)) {
            selected_value[0][j] = true;
            rank[i][0] = 0;
            rank[i][1] = j;
            i++;
          }
        }
        for (j = 0;j < ivec.nb_vector;j++) {
          if ((!selected_value[1][j]) && (ivec.int_vector[j][0] == min)) {
            selected_value[1][j] = true;
            rank[i][0] = 1;
            rank[i][1] = j;
            i++;
          }
        }
      }
      while (i < nb_vector + ivec.nb_vector);

      // computation of the sup-norm distance

      sample_size[0] = 0;
      sample_size[1] = 0;
      max = 0.;
      sup_norm = 0.;
      diff_cumul = 0.;
      sign = 0;
      nb_crossing = 0;

#     ifdef DEBUG
      cout << "\n";
#     endif

      i = 0;

      do {
        switch (rank[i][0]) {
        case 0 :
          int_value = int_vector[rank[i][1]][0];
          break;
        case 1 :
          int_value = ivec.int_vector[rank[i][1]][0];
          break;
        }
        sample_size[rank[i][0]]++;
        i++;

        while ((i < nb_vector + ivec.nb_vector) && (rank[i][0] == 0) &&
               (int_vector[rank[i][1]][0] == int_value)) {
          sample_size[0]++;
          i++;
        }
        while ((i < nb_vector + ivec.nb_vector) && (rank[i][0] == 1) &&
               (ivec.int_vector[rank[i][1]][0] == int_value)) {
          sample_size[1]++;
          i++;
        }

#       ifdef DEBUG
        previous_diff_cumul = diff_cumul;
#       endif

        diff_cumul = (double)sample_size[0] / (double)nb_vector - (double)sample_size[1] / (double)ivec.nb_vector;
        if (fabs(diff_cumul) > max) {
          max = fabs(diff_cumul);
          sup_value = int_value;
        }

#       ifdef DEBUG
        cout << i << "  " << int_value << "   "
             << (double)sample_size[0] / (double)nb_vector << "  " << (double)sample_size[1] / (double)ivec.nb_vector << endl;
#       endif

        // test of non-crossing of empirical cumulative distribution functions

#       ifdef DEBUG
        cout << i << "  " << diff_cumul << "   " << diff_cumul - previous_diff_cumul << endl;
#       endif

        previous_sign = sign;
        if (diff_cumul < 0.) {
          sign = -1;
        }
        else if (diff_cumul == 0.) {
          sign = 0;
        }
        else {
          sign = 1;
        }

        if ((previous_sign != 0) && (sign != 0) && (sign != previous_sign)) {
          if (display) {
            cout << "\n" << STAT_label[STATL_SUP_NORM_DISTANCE] << ": " << max
                 << "   sup norm value: " << sup_value << "   crossing value: " << int_value;
          }

          nb_crossing++;
          sup_norm += max;
          max = 0.;
        }
      }
      while (i < nb_vector + ivec.nb_vector);

      sup_norm += max;

#     ifdef DEBUG
      cout << endl;
#     endif

      break;
    }

    case REAL_VALUE : {

      // computation of ranks on the basis of the two samples

      i = 0;

      do {
        min = MAX(max_value[0], ivec.max_value[0]) + 1.;

        for (j = 0;j < nb_vector;j++) {
          if ((!selected_value[0][j]) && (real_vector[j][0] < min)) {
            min = real_vector[j][0];
            rank[i][0] = 0;
            rank[i][1] = j;
          }
        }
        for (j = 0;j < ivec.nb_vector;j++) {
          if ((!selected_value[1][j]) && (ivec.real_vector[j][0] < min)) {
            min = ivec.real_vector[j][0];
            rank[i][0] = 1;
            rank[i][1] = j;
          }
        }

        selected_value[rank[i][0]][rank[i][1]] = true;
        i++;

        for (j = 0;j < nb_vector;j++) {
          if ((!selected_value[0][j]) && (real_vector[j][0] == min)) {
            selected_value[0][j] = true;
            rank[i][0] = 0;
            rank[i][1] = j;
            i++;
          }
        }
        for (j = 0;j < ivec.nb_vector;j++) {
          if ((!selected_value[1][j]) && (ivec.real_vector[j][0] == min)) {
            selected_value[1][j] = true;
            rank[i][0] = 1;
            rank[i][1] = j;
            i++;
          }
        }
      }
      while (i < nb_vector + ivec.nb_vector);

      // computation of the sup-norm distance

      sample_size[0] = 0;
      sample_size[1] = 0;
      max = 0.;
      sup_norm = 0.;
      sign = 0;
      nb_crossing = 0;

#     ifdef DEBUG
      cout << "\n";
#     endif

      i = 0;

      do {
        switch (rank[i][0]) {
        case 0 :
          real_value = real_vector[rank[i][1]][0];
          break;
        case 1 :
          real_value = ivec.real_vector[rank[i][1]][0];
          break;
        }
        sample_size[rank[i][0]]++;
        i++;

        while ((i < nb_vector + ivec.nb_vector) && (rank[i][0] == 0) &&
               (real_vector[rank[i][1]][0] == real_value)) {
          sample_size[0]++;
          i++;
        }
        while ((i < nb_vector + ivec.nb_vector) && (rank[i][0] == 1) &&
               (ivec.real_vector[rank[i][1]][0] == real_value)) {
          sample_size[1]++;
          i++;
        }

        diff_cumul = (double)sample_size[0] / (double)nb_vector - (double)sample_size[1] / (double)ivec.nb_vector;
        if (fabs(diff_cumul) > max) {
          max = fabs(diff_cumul);
          sup_value = real_value;
        }

#       ifdef DEBUG
        cout << i << "  " << real_value << "   "
             << (double)sample_size[0] / (double)nb_vector << "  " << (double)sample_size[1] / (double)ivec.nb_vector << endl;
#       endif

        // test of non-crossing of empirical cumulative distribution functions

        previous_sign = sign;
        if (diff_cumul < 0.) {
          sign = -1;
        }
        else if (diff_cumul == 0.) {
          sign = 0;
        }
        else {
          sign = 1;
        }

        if ((previous_sign != 0) && (sign != 0) && (sign != previous_sign)) {
          if (display) {
            cout << "\n" << STAT_label[STATL_SUP_NORM_DISTANCE] << ": " << max
                 << "   sup norm value: " << sup_value << "   crossing value: " << real_value;
          }

          nb_crossing++;
          sup_norm += max;
          max = 0.;
        }
      }
      while (i < nb_vector + ivec.nb_vector);

      sup_norm += max;

#     ifdef DEBUG
      cout << endl;
#     endif

      break;
    }
    }

    if (display) {
      cout << "\n" << STAT_label[STATL_SUP_NORM_DISTANCE] << ": " << sup_norm << "   "
           << STAT_label[STATL_OVERLAP] << ": " << 1. - sup_norm
           << "   number of crossings of cdfs: " << nb_crossing << "   sup norm value: " << sup_value << endl;
    }

    delete [] selected_value[0];
    delete [] selected_value[1];
    delete [] selected_value;

    for (i = 0;i < nb_vector + ivec.nb_vector;i++) {
      delete [] rank[i];
    }
    delete [] rank;
  }

  return status;
}


};  // namespace stat_tool
