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
#include <iostream>

#include "stat_tool/stat_label.h"

#include "sequences.h"
#include "sequence_label.h"

using namespace std;
using namespace stat_tool;


namespace sequence_analysis {



/*--------------------------------------------------------------*/
/**
 *  \brief Writing of the alignment of 2 sequences.
 *
 *  \param[in,out] os              stream,
 *  \param[in]     width           column width,
 *  \param[in]     ref_index       reference sequence index,
 *  \param[in]     test_index      test sequence index,
 *  \param[in]     alignment       reference on the alignment,
 *  \param[in]     alignment_index alignment index.
 */
/*--------------------------------------------------------------*/

ostream& Sequences::alignment_ascii_print(ostream &os , int width , int ref_index , int test_index ,
                                          const Sequences &alignment , int alignment_index) const

{
  int i , j , k , m , n;
  int ref_rank , test_rank , alignment_rank;
  ios_base::fmtflags format_flags;


  format_flags = os.setf(ios::right , ios::adjustfield);

  os << "\n" << SEQ_label[SEQL_SEQUENCE] << " " << identifier[test_index]
     << " (" << SEQ_label[SEQL_LENGTH] << " " << length[test_index] << ") " << SEQ_label[SEQL_ALIGNED_ON]
     << " " << SEQ_label[SEQL_SEQUENCE] << " " << identifier[ref_index]
     << " (" << SEQ_label[SEQL_LENGTH] << " " << length[ref_index] << ")" << endl;

  ref_rank = 0;
  test_rank = 0;
  alignment_rank = 0;

  os << "\n";
  i = 0;
  for (j = 0;j < alignment.length[alignment_index];j++) {
    if ((alignment.int_sequence[alignment_index][0][j] != DELETION) &&
        (alignment.int_sequence[alignment_index][0][j] != BEGIN_END_DELETION)) {
      if (index_parameter) {
        os << setw(width) << index_parameter[test_index][i++];
      }

      else {
        if (type[0] != REAL_VALUE) {
          os << setw(width) << int_sequence[test_index][0][i++];
        }
        else {
          os << setw(width) << real_sequence[test_index][0][i++];
        }
      }
    }

    else if (alignment.int_sequence[alignment_index][0][j] == DELETION) {
      os << setw(width) << "-";
    }

    else {
      os << setw(width) << " ";
    }
    os << " ";

    if (((j - alignment_rank) * (width + 1) > LINE_NB_CHARACTER) ||
        (j == alignment.length[alignment_index] - 1)) {
      if (j < alignment.length[alignment_index] - 1) {
        os << "\\";
      }
      os << endl;

      // test sequence

      for (k = (index_parameter ? 0 : 1);k < nb_variable;k++) {
        m = test_rank;
        for (n = alignment_rank;n <= j;n++) {
          if ((alignment.int_sequence[alignment_index][0][n] != DELETION) &&
              (alignment.int_sequence[alignment_index][0][n] != BEGIN_END_DELETION)) {
            if (type[k] != REAL_VALUE) {
              os << setw(width) << int_sequence[test_index][k][m++];
            }
            else {
              os << setw(width) << real_sequence[test_index][k][m++];
            }
          }

          else if (alignment.int_sequence[alignment_index][0][n] == DELETION) {
            os << setw(width) << "-";
          }

          else {
            os << setw(width) << " ";
          }
          os << " ";
        }

        if (j < alignment.length[alignment_index] - 1) {
          os << "\\";
        }
        os << endl;
      }
      os << endl;

      if (j < alignment.length[alignment_index] - 1) {
        test_rank = m;
      }

      // edit operations

      for (k = alignment_rank;k <= j;k++) {
        switch (alignment.int_sequence[alignment_index][0][k]) {
        case DELETION :
          os << setw(width) << "d";
          break;
        case INSERTION :
          os << setw(width) << "i";
          break;
        case MATCH :
          os << setw(width) << "|";
          break;
        case SUBSTITUTION :
          os << setw(width) << "s";
          break;
        case TRANSPOSITION :
          os << setw(width) << "t";
          break;
        default :
          os << setw(width) << " ";
          break;
        }
        os << " ";
      }

      if (j < alignment.length[alignment_index] - 1) {
        os << "\\";
      }
      os << endl;

      // reference sequence

      os << "\n";
      if (index_parameter) {
        k = ref_rank;
        for (m = alignment_rank;m <= j;m++) {
          if ((alignment.int_sequence[alignment_index][0][m] != INSERTION) &&
              (alignment.int_sequence[alignment_index][0][m] != BEGIN_END_INSERTION)) {
            os << setw(width) << index_parameter[ref_index][k++];
          }
          else if (alignment.int_sequence[alignment_index][0][m] == INSERTION) {
            os << setw(width) << "-";
          }
          else {
            os << setw(width) << " ";
          }
          os << " ";
        }

        if (j < alignment.length[alignment_index] - 1) {
          os << "\\";
        }
        os << endl;
      }

      for (k = 0;k < nb_variable;k++) {
        m = ref_rank;
        for (n = alignment_rank;n <= j;n++) {
          if ((alignment.int_sequence[alignment_index][0][n] != INSERTION) &&
              (alignment.int_sequence[alignment_index][0][n] != BEGIN_END_INSERTION)) {
            if (type[k] != REAL_VALUE) {
              os << setw(width) << int_sequence[ref_index][k][m++];
            }
            else {
              os << setw(width) << real_sequence[ref_index][k][m++];
            }
          }

          else if (alignment.int_sequence[alignment_index][0][n] == INSERTION) {
            os << setw(width) << "-";
          }

          else {
            os << setw(width) << " ";
          }
          os << " ";
        }

        if (j < alignment.length[alignment_index] - 1) {
          os << "\\";
        }
        os << endl;
      }

      if (j < alignment.length[alignment_index] - 1) {
        alignment_rank = j + 1;
        ref_rank = m;
      }
      os << endl;
    }
  }

  os.setf(format_flags , ios::adjustfield);

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of the alignment of 2 sequences at the spreadsheet format.
 *
 *  \param[in,out] os              stream,
 *  \param[in]     ref_index       reference sequence index,
 *  \param[in]     test_index      test sequence index,
 *  \param[in]     alignment       reference on the alignment,
 *  \param[in]     alignment_index alignment index.
 */
/*--------------------------------------------------------------*/

ostream& Sequences::alignment_spreadsheet_print(ostream &os , int ref_index , int test_index ,
                                                const Sequences &alignment , int alignment_index) const

{
  int i , j , k;


  os << "\n" << SEQ_label[SEQL_SEQUENCE] << " " << identifier[test_index]
     << "\t" << SEQ_label[SEQL_LENGTH] << " " << length[test_index] << "\t" << SEQ_label[SEQL_ALIGNED_ON]
     << "\t" << SEQ_label[SEQL_SEQUENCE] << " " << identifier[ref_index]
     << "\t" << SEQ_label[SEQL_LENGTH] << " " << length[ref_index] << endl;

  // test sequence

  os << "\n";
  if (index_parameter) {
    i = 0;
    for (j = 0;j < alignment.length[alignment_index];j++) {
      if ((alignment.int_sequence[alignment_index][0][j] != DELETION) &&
          (alignment.int_sequence[alignment_index][0][j] != BEGIN_END_DELETION)) {
        os << index_parameter[test_index][i++];
      }
      else if (alignment.int_sequence[alignment_index][0][j] == DELETION) {
        os << "-";
      }
      os << "\t";
    }
    os << endl;
  }

  for (i = 0;i < nb_variable;i++) {
    j = 0;
    for (k = 0;k < alignment.length[alignment_index];k++) {
      if ((alignment.int_sequence[alignment_index][0][k] != DELETION) &&
          (alignment.int_sequence[alignment_index][0][k] != BEGIN_END_DELETION)) {
        if (type[i] != REAL_VALUE) {
          os << int_sequence[test_index][i][j++];
        }
        else {
          os << real_sequence[test_index][i][j++];
        }
      }

      else if (alignment.int_sequence[alignment_index][0][k] == DELETION) {
        os << "-";
      }
      os << "\t";
    }
    os << endl;
  }
  os << endl;

  // edit operations

  for (i = 0;i < alignment.length[alignment_index];i++) {
    switch (alignment.int_sequence[alignment_index][0][i]) {
    case DELETION :
      os << "d";
      break;
    case INSERTION :
      os << "i";
      break;
    case MATCH :
      os << "|";
      break;
    case SUBSTITUTION :
      os << "s";
      break;
    case TRANSPOSITION :
      os << "t";
      break;
    default :
      os << " ";
      break;
    }
    os << "\t";
  }
  os << endl;

  // reference sequence

  os << "\n";
  if (index_parameter) {
    i = 0;
    for (j = 0;j < alignment.length[alignment_index];j++) {
      if ((alignment.int_sequence[alignment_index][0][j] != INSERTION) &&
          (alignment.int_sequence[alignment_index][0][j] != BEGIN_END_INSERTION)) {
        os << index_parameter[ref_index][i++];
      }
      else if (alignment.int_sequence[alignment_index][0][j] == INSERTION) {
        os << "-";
      }
      os << "\t";
    }
    os << endl;
  }

  for (i = 0;i < nb_variable;i++) {
    j = 0;
    for (k = 0;k < alignment.length[alignment_index];k++) {
      if ((alignment.int_sequence[alignment_index][0][k] != INSERTION) &&
          (alignment.int_sequence[alignment_index][0][k] != BEGIN_END_INSERTION)) {
        if (type[i] != REAL_VALUE) {
          os << int_sequence[ref_index][i][j++];
        }
        else {
          os << real_sequence[ref_index][i][j++];
        }
      }

      else if (alignment.int_sequence[alignment_index][0][k] == INSERTION) {
        os << "-";
      }
      os << "\t";
    }
    os << endl;
  }
  os << endl;

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of a fixed insertion/deletion cost.
 *
 *  \param[in] vector_dist           reference on a VectorDistance object,
 *  \param[in] rank                  ranks (for ordinal variables),
 *  \param[in] max_category_distance maximum distances between categories.
 *
 *  \return                          insertion/deletion cost.
 */
/*--------------------------------------------------------------*/

double Sequences::indel_distance_computation(const VectorDistance &vector_dist ,
                                             double **rank , double **max_category_distance) const

{
  int i , j;
  double ldistance , distance = 0.;


  for (i = 0;i < vector_dist.get_nb_variable();i++) {
    switch (vector_dist.get_var_type(i)) {

    case NOMINAL : {
      if (!max_category_distance[i]) {
        ldistance = 1.;
      }

      else {
        ldistance = 0.;
        for (j = (int)min_value[i];j <= (int)max_value[i];j++) {
          if (max_category_distance[i][j] > ldistance) {
            ldistance = max_category_distance[i][j];
          }
        }
      }
      break;
    }

    case ORDINAL : {
      ldistance = rank[i][(int)max_value[i]] - rank[i][(int)min_value[i]];
      break;
    }

    case NUMERIC : {
      ldistance = max_value[i] - min_value[i];
      break;
    }

    case CIRCULAR : {
      ldistance = MIN(max_value[i] - min_value[i] , vector_dist.get_period(i) / 2.);
      break;
    }
    }

    switch (vector_dist.get_distance_type()) {
    case ABSOLUTE_VALUE :
      distance += vector_dist.get_weight(i) * fabs(ldistance) / vector_dist.get_dispersion(i);
      break;
    case QUADRATIC :
      distance += vector_dist.get_weight(i) * ldistance * ldistance / vector_dist.get_dispersion(i);
      break;
    }
  }

  if (vector_dist.get_distance_type() == QUADRATIC) {
    distance = sqrt(distance);
  }

  return distance;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the insertion/deletion cost of a vector.
 *
 *  \param[in] vector_dist           reference on a VectorDistance object,
 *  \param[in] index                 sequence index,
 *  \param[in] position              position in the sequence,
 *  \param[in] rank                  ranks (for ordinal variables),
 *  \param[in] max_category_distance maximum distances between categories.
 *
 *  \return                          insertion/deletion cost.
 */
/*--------------------------------------------------------------*/

double Sequences::indel_distance_computation(const VectorDistance &vector_dist ,
                                             int index , int position , double **rank ,
                                             double **max_category_distance) const

{
  int i;
  double ldistance , distance = 0.;


  for (i = 0;i < vector_dist.get_nb_variable();i++) {
    switch (vector_dist.get_var_type(i)) {

    case NOMINAL : {
      if (!max_category_distance[i]) {
        ldistance = 1.;
      }
      else {
        ldistance = max_category_distance[i][int_sequence[index][i][position]];
      }
      break;
    }

    case ORDINAL : {
      ldistance = MAX(rank[i][int_sequence[index][i][position]] - rank[i][(int)min_value[i]] ,
                      rank[i][(int)max_value[i]] - rank[i][int_sequence[index][i][position]]);
      break;
    }

    case NUMERIC : {
      if (type[i] != REAL_VALUE) {
        ldistance = MAX(int_sequence[index][i][position] - (int)min_value[i] ,
                        (int)max_value[i] - int_sequence[index][i][position]);
      }
      else {
        ldistance = MAX(real_sequence[index][i][position] - min_value[i] ,
                        max_value[i] - real_sequence[index][i][position]);
      }
      break;
    }

    case CIRCULAR : {
      ldistance = MAX(int_sequence[index][i][position] - (int)min_value[i] ,
                      (int)max_value[i] - int_sequence[index][i][position]);
      if (ldistance > vector_dist.get_period(i) / 2.) {
        ldistance = vector_dist.get_period(i) / 2.;
      }
      break;
    }
    }

    switch (vector_dist.get_distance_type()) {
    case ABSOLUTE_VALUE :
      distance += vector_dist.get_weight(i) * fabs(ldistance) / vector_dist.get_dispersion(i);
      break;
    case QUADRATIC :
      distance += vector_dist.get_weight(i) * ldistance * ldistance / vector_dist.get_dispersion(i);
      break;
    }
  }

  if (vector_dist.get_distance_type() == QUADRATIC) {
    distance = sqrt(distance);
  }

  return distance;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the cost of substitution of one vector by another.
 *
 *  \param[in] vector_dist   reference on a VectorDistance object,
 *  \param[in] ref_index     reference sequence index,
 *  \param[in] test_index    test sequence index,
 *  \param[in] ref_position  position in the reference sequence,
 *  \param[in] test_position position in the test sequence,
 *  \param[in] rank          ranks (for ordinal variables),
 *  \param[in] test_seq      pointer on the test sequences (multiple alignment).
 *
 *  \return                  substitution cost.
 */
/*--------------------------------------------------------------*/

double Sequences::substitution_distance_computation(const VectorDistance &vector_dist ,
                                                    int ref_index , int test_index ,
                                                    int ref_position , int test_position ,
                                                    double **rank , const Sequences *test_seq) const

{
  int i;
  double ldistance , distance = 0.;


  if (!test_seq) {
    test_seq = this;
  }

  for (i = 0;i < vector_dist.get_nb_variable();i++) {
    switch (vector_dist.get_var_type(i)) {

    case NOMINAL : {
      if (!vector_dist.get_category_distance(i)) {
        ldistance = (int_sequence[ref_index][i][ref_position] == test_seq->int_sequence[test_index][i][test_position] ? 0. : 1.);
      }
      else {
        ldistance = vector_dist.get_category_distance(i , int_sequence[ref_index][i][ref_position] ,
                                                      test_seq->int_sequence[test_index][i][test_position]);
      }
      break;
    }

    case ORDINAL : {
      ldistance = rank[i][int_sequence[ref_index][i][ref_position]] - rank[i][test_seq->int_sequence[test_index][i][test_position]];
      break;
    }

    case NUMERIC : {
      if (type[i] != REAL_VALUE) {
        ldistance = int_sequence[ref_index][i][ref_position] - test_seq->int_sequence[test_index][i][test_position];
      }
      else {
        ldistance = real_sequence[ref_index][i][ref_position] - test_seq->real_sequence[test_index][i][test_position];
      }
      break;
    }

    case CIRCULAR : {
      if (int_sequence[ref_index][i][ref_position] <= test_seq->int_sequence[test_index][i][test_position]) {
        ldistance = MIN(test_seq->int_sequence[test_index][i][test_position] - int_sequence[ref_index][i][ref_position] ,
                        int_sequence[ref_index][i][ref_position] + vector_dist.get_period(i) -
                        test_seq->int_sequence[test_index][i][test_position]);
      }
      else {
        ldistance = MIN(int_sequence[ref_index][i][ref_position] - test_seq->int_sequence[test_index][i][test_position] ,
                        test_seq->int_sequence[test_index][i][test_position] + vector_dist.get_period(i) -
                        int_sequence[ref_index][i][ref_position]);
      }
      break;
    }
    }

    switch (vector_dist.get_distance_type()) {
    case ABSOLUTE_VALUE :
      distance += vector_dist.get_weight(i) * fabs(ldistance) / vector_dist.get_dispersion(i);
      break;
    case QUADRATIC :
      distance += vector_dist.get_weight(i) * ldistance * ldistance / vector_dist.get_dispersion(i);
      break;
    }
  }

  if (vector_dist.get_distance_type() == QUADRATIC) {
    distance = sqrt(distance);
  }

  return distance;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Alignment of sequences.
 *
 *  \param[in] error                reference on a StatError object,
 *  \param[in] display              flag for displaying the alignments,
 *  \param[in] ivector_dist         reference on a VectorDistance object,
 *  \param[in] ref_identifier       reference sequence identifier,
 *  \param[in] test_identifier      test sequence identifier,
 *  \param[in] begin_free           flag begin-free alignment,
 *  \param[in] end_free             flag end-free alignment,
 *  \param[in] indel_cost           insertion/deletion costs adaptative or fixed,
 *  \param[in] indel_factor         factor for deducing the insertion/deletion costs,
 *  \param[in] transposition_flag   flag transposition,
 *  \param[in] transposition_factor factor for deducing the transposition costs,
 *  \param[in] result_path          result file path,
 *  \param[in] result_format        result file format (ASCII/SPREADSHEET),
 *  \param[in] alignment_path       alignment file path.
 *
 *  \return                         DistanceMatrix object.
 */
/*--------------------------------------------------------------*/

DistanceMatrix* Sequences::alignment(StatError &error , bool display , const VectorDistance &ivector_dist ,
                                     int ref_identifier , int test_identifier , bool begin_free ,
                                     bool end_free , insertion_deletion_cost indel_cost , double indel_factor ,
                                     bool transposition_flag , double transposition_factor ,
                                     const string result_path , output_format result_format ,
                                     const string alignment_path) const

{
  bool status = true , half_matrix;
  int i , j , k , m;
  int nb_alignment , ilength , alignment_index , var , width , ref_position , pref_position ,
      test_position , ptest_position , gap_length , max_gap_length , nb_deletion , nb_insertion ,
      nb_match , nb_substitution , nb_transposition , nb_begin_end , offset , *palignment ,
      *calignment , *category , **path_length , ***back_pointers;
  double buff , deletion_distance , insertion_distance , substitution_distance ,
         transposition_distance , max_transposition_cost , **rank , **max_category_distance ,
         **local_indel_distance , **local_substitution_distance , **cumul_distance;
  VectorDistance *vector_dist;
  DistanceMatrix *dist_matrix;
  Sequences *alignment;
  ofstream *out_file;


  dist_matrix = NULL;
  error.init();

  if (nb_sequence < 2) {
    status = false;
    error.update(SEQ_error[SEQR_NB_SEQUENCE]);
  }

  if (((index_param_type == TIME) && (index_interval->variance > 0.)) ||
      (index_param_type == POSITION)) {
    status = false;
    error.update(SEQ_error[SEQR_INDEX_PARAMETER_TYPE]);
  }

  if (ivector_dist.get_nb_variable() != nb_variable) {
    status = false;
    error.update(STAT_error[STATR_NB_VARIABLE]);
  }

  else {
    for (i = 0;i < nb_variable;i++) {
      if (ivector_dist.get_var_type(i) != NUMERIC) {
        if ((type[i] != INT_VALUE) && (type[i] != STATE)) {
          status = false;
          ostringstream error_message , correction_message;
          error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                        << STAT_error[STATR_VARIABLE_TYPE];
          correction_message << STAT_variable_word[INT_VALUE] << " or "
                             << STAT_variable_word[STATE];
          error.correction_update((error_message.str()).c_str() , (correction_message.str()).c_str());
        }

        else if (!marginal_distribution[i]) {
          status = false;
          ostringstream error_message;
          error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                        << STAT_error[STATR_MARGINAL_FREQUENCY_DISTRIBUTION];
          error.update((error_message.str()).c_str());
        }

        if ((ivector_dist.get_var_type(i) == NOMINAL) &&
            ((min_value[i] < 0) || (max_value[i] >= NB_CATEGORY) ||
             ((ivector_dist.get_category_distance(i)) && (ivector_dist.get_nb_value(i) != max_value[i] + 1)))) {
          status = false;
          ostringstream error_message;
          error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                        << STAT_error[STATR_NB_CATEGORY];
          error.update((error_message.str()).c_str());
        }

        if ((ivector_dist.get_var_type(i) == CIRCULAR) &&
            (max_value[i] - min_value[i] >= ivector_dist.get_period(i))) {
          status = false;
          ostringstream error_message;
          error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                        << STAT_error[STATR_NB_VALUE_PERIOD];
          error.update((error_message.str()).c_str());
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

  if (ref_identifier != I_DEFAULT) {
    for (i = 0;i < nb_sequence;i++) {
      if (ref_identifier == identifier[i]) {
        break;
      }
    }

    if (i == nb_sequence) {
      status = false;
      error.update(SEQ_error[SEQR_REF_SEQUENCE_IDENTIFIER]);
    }
  }

  if (test_identifier != I_DEFAULT) {
    for (i = 0;i < nb_sequence;i++) {
      if (test_identifier == identifier[i]) {
        break;
      }
    }

    if (i == nb_sequence) {
      status = false;
      error.update(SEQ_error[SEQR_TEST_SEQUENCE_IDENTIFIER]);
    }
  }

  if ((ref_identifier != I_DEFAULT) && (test_identifier != I_DEFAULT) && (ref_identifier == test_identifier)) {
    status = false;
    error.correction_update(SEQ_error[SEQR_SEQUENCE_IDENTIFIERS] , "different");
  }

  nb_alignment = 1;
  if (ref_identifier == I_DEFAULT) {
    nb_alignment *= (nb_sequence - 1);
  }
  if (test_identifier == I_DEFAULT) {
    nb_alignment *= (ref_identifier != I_DEFAULT ? nb_sequence - 1 : nb_sequence);
  }

  if (nb_alignment > NB_ALIGNMENT) {
    status = false;
    error.update(SEQ_error[SEQR_NB_ALIGNMENT]);
  }

  if (indel_factor <= 0.5) {
    status = false;
    error.update(SEQ_error[SEQR_INDEL_FACTOR]);
  }

  if ((transposition_factor < 0.) || (transposition_factor >= 2.)) {
    status = false;
    error.update(SEQ_error[SEQR_TRANSPOSITION_FACTOR]);
  }

  if (status) {
    if ((ref_identifier == I_DEFAULT) && (test_identifier == I_DEFAULT) && ((result_path.empty()) ||
         (nb_alignment > FILE_NB_ALIGNMENT)) && (alignment_path.empty())) {
      half_matrix = true;
    }
    else {
      half_matrix = false;
    }

    vector_dist = new VectorDistance(ivector_dist);

    // computation of the maximum substitution distance for nominal variables and
    // the ranks for ordinal variables

    rank = new double*[nb_variable];
    max_category_distance = new double*[nb_variable];

    for (i = 0;i < nb_variable;i++) {
      if ((vector_dist->get_var_type(i) == NOMINAL) && (vector_dist->get_category_distance(i))) {
        max_category_distance[i] = vector_dist->max_category_distance_computation(i);
      }
      else {
        max_category_distance[i] = 0;
      }

      if (vector_dist->get_var_type(i) == ORDINAL) {
        rank[i] = marginal_distribution[i]->rank_computation();
      }
      else {
        rank[i] = NULL;
      }

      // computation of dispersion measures for the standardization of variables

      if (marginal_distribution[i]) {
        vector_dist->dispersion_computation(i , marginal_distribution[i] , rank[i]);
      }

      else {
        switch (vector_dist->get_distance_type()) {
        case ABSOLUTE_VALUE :
          vector_dist->dispersion_update(i , mean_absolute_difference_computation(i));
          break;
        case QUADRATIC :
          vector_dist->dispersion_update(i , 2 * variance_computation(i , mean_computation(i)));
          break;
        }

        if (vector_dist->get_dispersion(i) == 0.) {
          vector_dist->dispersion_update(i , 1.);
        }
      }
    }

#   ifdef DEBUG
    cout << *vector_dist;
    if (vector_dist->get_distance_type() == ABSOLUTE_VALUE) {
      for (i = 0;i < nb_variable;i++) {
        if (vector_dist->get_var_type(i) == NUMERIC) {
          cout << "\n" << STAT_label[STATL_VARIABLE] << " " << i << "   mean absolute difference: "
               << mean_absolute_difference_computation(i) << endl;
        }
      }
    }
#   endif

    if (index_parameter) {
      width = column_width(index_parameter_distribution->nb_value - 1);
    }
    else {
      width = 0;
    }

    for (i = 0;i < nb_variable;i++) {
      if (type[i] != REAL_VALUE) {
        var = column_width((int)min_value[i] , (int)max_value[i]);
        if (var > width) {
          width = var;
        }
      }

      else {
        for (j = 0;j < nb_sequence;j++) {
          if ((ref_identifier == I_DEFAULT) || (ref_identifier == identifier[j]) ||
              (test_identifier == I_DEFAULT) || (test_identifier == identifier[j])) {
            var = column_width(length[j] , real_sequence[j][i]);
            if (var > width) {
              width = var;
            }
          }
        }
      }
    }

    out_file = NULL;

    if (!result_path.empty()) {
      out_file = new ofstream(result_path.c_str());

      if (!out_file) {
        error.update(STAT_error[STATR_FILE_NAME]);
        if (display) {
          cout << error;
        }
      }
    }

    // construction of the result data structures

    dist_matrix = new DistanceMatrix(nb_sequence , ref_identifier , test_identifier ,
                                     SEQ_label[SEQL_SEQUENCE] , identifier ,
                                     true , transposition_flag);

    if (!alignment_path.empty()) {
      alignment = new Sequences(nb_alignment , 1);
    }
    else {
      ilength = max_length + max_length;
      alignment = new Sequences(1 , NULL , &ilength , 1 , false);
    }

    // construction of the algorithm data structures - computation of the insertion/deletion costs

    if (indel_cost == FIXED) {
      buff = indel_distance_computation(*vector_dist , rank , max_category_distance) * indel_factor;
    }

    local_indel_distance = new double*[nb_sequence];
    for (i = 0;i < nb_sequence;i++) {
      if ((ref_identifier == I_DEFAULT) || (ref_identifier == identifier[i]) ||
          (test_identifier == I_DEFAULT) || (test_identifier == identifier[i])) {
        local_indel_distance[i] = new double[length[i] + 1];

        switch (indel_cost) {

        case ADAPTATIVE : {
          for (j = 1;j <= length[i];j++) {
            local_indel_distance[i][j] = indel_distance_computation(*vector_dist , i , j - 1 , rank , max_category_distance) * indel_factor;
          }

#         ifdef DEBUG
/*          cout << "\ninsertion/deletion cost of the vectors of the sequence " << i << ": ";
          for (j = 1;j <= length[i];j++) {
            cout << local_indel_distance[i][j] << " ";
          }
          cout << endl; */
#         endif

          break;
        }

        case FIXED : {
          for (j = 1;j <= length[i];j++) {
            local_indel_distance[i][j] = buff;
          }
          break;
        }
        }
      }

      else {
        local_indel_distance[i] = NULL;
      }
    }

    local_substitution_distance = new double*[max_length + 1];
    local_substitution_distance[0] = NULL;
    for (i = 1;i <= max_length;i++) {
      local_substitution_distance[i] = new double[max_length + 1];
    }

    cumul_distance = new double*[max_length + 1];
    for (i = 0;i <= max_length;i++) {
      cumul_distance[i] = new double[max_length + 1];
    }

    path_length = new int*[max_length + 1];
    for (i = 0;i <= max_length;i++) {
      path_length[i] = new int[max_length + 1];
    }

    back_pointers = new int**[max_length + 1];
    for (i = 0;i <= max_length;i++) {
      back_pointers[i] = new int*[max_length + 1];
      for (j = 0;j <= max_length;j++) {
        back_pointers[i][j] = new int[2];
      }
    }

    // alignment of sequences

#   ifdef DEBUG
    int nb_local_substitution_distance = 0;
    double mean_local_substitution_distance = 0.;
#   endif

    alignment_index = 0;

    for (i = 0;i < nb_sequence;i++) {
      if ((ref_identifier == I_DEFAULT) || (ref_identifier == identifier[i])) {
        for (j = (half_matrix ? i + 1 : 0);j < nb_sequence;j++) {
          if (((test_identifier == I_DEFAULT) || (test_identifier == identifier[j])) && (j != i)) {

            // initialization of the cumulative distances and the corresponding alignment lengths

            cumul_distance[0][0] = 0.;
            path_length[0][0] = 0;

            for (k = 1;k <= length[i];k++) {

              // deletion

              if (begin_free) {
                cumul_distance[k][0] = cumul_distance[k - 1][0];
              }
              else {
                cumul_distance[k][0] = cumul_distance[k - 1][0] + local_indel_distance[i][k];
              }

              path_length[k][0] = k;
              back_pointers[k][0][0] = k - 1;
              back_pointers[k][0][1] = 0;
            }

            for (k = 1;k <= length[j];k++) {

              // insertion

              if (begin_free) {
                cumul_distance[0][k] = cumul_distance[0][k - 1];
              }
              else {
                cumul_distance[0][k] = cumul_distance[0][k - 1] + local_indel_distance[j][k];
              }

              path_length[0][k] = k;
              back_pointers[0][k][0] = 0;
              back_pointers[0][k][1] = k - 1;
            }

            // computation of the cumulative distances and the corresponding alignment lengths

            for (k = 1;k <= length[i];k++) {
              for (m = 1;m <= length[j];m++) {

                // computation of the distance of substitution of one vector by another

                local_substitution_distance[k][m] = substitution_distance_computation(*vector_dist , i , j , k - 1 , m - 1 , rank);

#               ifdef DEBUG
                nb_local_substitution_distance++;
                mean_local_substitution_distance += local_substitution_distance[k][m];
#               endif

                // match/substitution

                cumul_distance[k][m] = cumul_distance[k - 1][m - 1] + local_substitution_distance[k][m];
                path_length[k][m] = path_length[k - 1][m - 1] + 1;
                back_pointers[k][m][0] = k - 1;
                back_pointers[k][m][1] = m - 1;

                // deletion

                if ((m < length[j]) || (!end_free)) {
                  buff = cumul_distance[k - 1][m] + local_indel_distance[i][k];
                }
                else {
                  buff = cumul_distance[k - 1][m];
                }

                if (buff < cumul_distance[k][m]) {
                  cumul_distance[k][m] = buff;
                  path_length[k][m] = path_length[k - 1][m] + 1;
                  back_pointers[k][m][0] = k - 1;
                  back_pointers[k][m][1] = m;
                }

                // insertion

                if ((k < length[i]) || (!end_free)) {
                  buff = cumul_distance[k][m - 1] + local_indel_distance[j][m];
                }
                else {
                  buff = cumul_distance[k][m - 1];
                }

                if (buff < cumul_distance[k][m]) {
                  cumul_distance[k][m] = buff;
                  path_length[k][m] = path_length[k][m - 1] + 1;
                  back_pointers[k][m][0] = k;
                  back_pointers[k][m][1] = m - 1;
                }

                // transposition

                if ((transposition_flag) && (k > 1) && (m > 1) &&
                    (local_substitution_distance[k][m] > 0.) &&
                    (local_substitution_distance[k - 1][m] == 0.) &&
                    (local_substitution_distance[k][m - 1] == 0.)) {
                  max_transposition_cost = local_substitution_distance[k][m];
                  if (local_indel_distance[i][k - 1] < max_transposition_cost) {
                    max_transposition_cost = local_indel_distance[i][k - 1];
                  }
                  if (local_indel_distance[i][k] < max_transposition_cost) {
                    max_transposition_cost = local_indel_distance[i][k];
                  }
                  buff = cumul_distance[k - 2][m - 2] + max_transposition_cost * transposition_factor;

                  if (buff < cumul_distance[k][m]) {
                    cumul_distance[k][m] = buff;
                    path_length[k][m] = path_length[k - 2][m - 2] + 2;
                    back_pointers[k][m][0] = k - 2;
                    back_pointers[k][m][1] = m - 2;
                  }
                }
              }
            }

            // end free (alternative implementation)

/*            if (end_free) {
              buff = cumul_distance[length[i]][length[j]];

              k = length[i];
              for (m = 1;m < length[i];m++) {
                if (cumul_distance[m][length[j]] < buff) {
                  buff = cumul_distance[m][length[j]];
                  k = m;
                }
              }
              for (m = k + 1;m <= length[i];m++) {
                cumul_distance[m][length[j]] = buff;
                path_length[m][length[j]] = path_length[m - 1][length[j]] + 1;
                back_pointers[m][length[j]][0] = m - 1;
                back_pointers[m][length[j]][1] = length[j];
              }

              k = length[j];
              for (m = 1;m < length[j];m++) {
                if (cumul_distance[length[i]][m] < buff) {
                  buff = cumul_distance[length[i]][m];
                  k = m;
                }
              }
              for (m = k + 1;m <= length[j];m++) {
                cumul_distance[length[i]][m] = buff;
                path_length[length[i]][m] = path_length[length[i]][m - 1] + 1;
                back_pointers[length[i]][m][0] = length[i];
                back_pointers[length[i]][m][1] = m - 1;
              }
            } */

#           ifdef DEBUG
/*            cout << "\n";
            for (k = length[i];k >= 0;k--) {
              for (m = 0;m <= length[j];m++) {
                cout << cumul_distance[k][m] << "  ";
              }
              cout << endl;
            }
            cout << endl; */
#           endif

            alignment->length[alignment_index] = path_length[length[i]][length[j]];
            if (!alignment_path.empty()) {
              alignment->int_sequence[alignment_index][0] = new int[alignment->length[alignment_index]];
            }

            // backtracking

            deletion_distance = 0.;
            insertion_distance = 0.;
            substitution_distance = 0.;
            transposition_distance = 0.;

            palignment = alignment->int_sequence[alignment_index][0] + alignment->length[alignment_index];
            k = path_length[length[i]][length[j]];
            pref_position = length[i];
            ptest_position = length[j];

#           ifdef DEBUG
//            cout << pref_position << " " << ptest_position << endl;
#           endif

            do {
              ref_position = pref_position;
              test_position = ptest_position;
              pref_position = back_pointers[ref_position][test_position][0];
              ptest_position = back_pointers[ref_position][test_position][1];

#             ifdef DEBUG
//              cout << pref_position << " " << ptest_position << endl;
#             endif

              if (test_position == ptest_position) {
                if (((test_position > 0) || (!begin_free)) && ((test_position < length[j]) || (!end_free))) {
                  *--palignment = DELETION;
                  deletion_distance += local_indel_distance[i][ref_position];
                }
                else {
                  *--palignment = BEGIN_END_DELETION;
                }
                k--;
              }

              else if (ref_position == pref_position) {
                if (((ref_position > 0) || (!begin_free)) && ((ref_position < length[i]) || (!end_free))) {
                  *--palignment = INSERTION;
                  insertion_distance += local_indel_distance[j][test_position];
                }
                else {
                  *--palignment = BEGIN_END_INSERTION;
                }
                k--;
              }

              else if ((ref_position == pref_position + 1) && (test_position == ptest_position + 1)) {
                if (local_substitution_distance[ref_position][test_position] == 0.) {
                  *--palignment = MATCH;
                }
                else {
                  *--palignment = SUBSTITUTION;
                  substitution_distance += local_substitution_distance[ref_position][test_position];
                }
                k--;
              }

              else if ((ref_position == pref_position + 2) && (test_position == ptest_position + 2)) {
                *--palignment = TRANSPOSITION;
                *--palignment = TRANSPOSITION;
                transposition_distance += local_substitution_distance[ref_position][test_position] *
                                          transposition_factor;
                k -= 2;
              }
            }
            while (k > 0);

            // search for the maximum number of successive insertions/deletions

            palignment = alignment->int_sequence[alignment_index][0];
            max_gap_length = 0;
            gap_length = 0;

            if ((*palignment == DELETION) || (*palignment == INSERTION)) {
              gap_length++;
            }

            for (k = 1;k < alignment->length[alignment_index];k++) {
              if (*(palignment + 1) != *palignment) {
                if (((*palignment == DELETION) || (*palignment == INSERTION)) &&
                    (gap_length > max_gap_length)) {
                  max_gap_length = gap_length;
                }
                gap_length = 0;
              }

              palignment++;
              if ((*palignment == DELETION) || (*palignment == INSERTION)) {
                gap_length++;
              }
            }

            if (((*palignment == DELETION) || (*palignment == INSERTION)) &&
                (gap_length > max_gap_length)) {
              max_gap_length = gap_length;
            }

            // update of the numbers of deletions, insertions, matchs, substitutions and transpositions

            palignment = alignment->int_sequence[alignment_index][0];

            nb_deletion = 0;
            nb_insertion = 0;
            nb_match = 0;
            nb_substitution = 0;
            nb_transposition = 0;
            nb_begin_end = 0;

            for (k = 0;k < alignment->length[alignment_index];k++) {
              switch (*palignment++) {
              case DELETION :
                nb_deletion++;
                break;
              case INSERTION :
                nb_insertion++;
                break;
              case MATCH :
                nb_match++;
                break;
              case SUBSTITUTION :
                nb_substitution++;
                break;
              case TRANSPOSITION :
                nb_transposition++;
                break;
              default :
                nb_begin_end++;
                break;
              }
            }

            dist_matrix->update(identifier[i] , identifier[j] , cumul_distance[length[i]][length[j]] ,
                                length[i] + length[j] - nb_begin_end , deletion_distance , nb_deletion ,
                                insertion_distance , nb_insertion , nb_match ,
                                substitution_distance , nb_substitution ,
                                transposition_distance , nb_transposition);

            if (half_matrix) {
              dist_matrix->update(identifier[j] , identifier[i] , cumul_distance[length[i]][length[j]] ,
                                  length[i] + length[j] - nb_begin_end , insertion_distance , nb_insertion ,
                                  deletion_distance , nb_deletion , nb_match ,
                                  substitution_distance , nb_substitution ,
                                  transposition_distance , nb_transposition);
            }

#           ifdef DEBUG
            {
              double sum = deletion_distance + insertion_distance + substitution_distance +
                           transposition_distance;
              if ((sum < cumul_distance[length[i]][length[j]] - DOUBLE_ERROR) ||
                  (sum > cumul_distance[length[i]][length[j]] + DOUBLE_ERROR)) {
                cout << "\nERROR: " << SEQ_label[SEQL_SEQUENCE] << " " << j + 1 << " aligned on "
                     << SEQ_label[SEQL_SEQUENCE] << " " << i + 1 << ": " << cumul_distance[length[i]][length[j]]
                     << " | " << sum << endl;
              }
            }
#           endif

            // writing of the alignment

            if ((display) && (nb_alignment <= DISPLAY_NB_ALIGNMENT)) {
              alignment_ascii_print(cout , width , i , j , *alignment , alignment_index);

              if (length[i] + length[j] - nb_begin_end > 0) {
                cout << STAT_label[STATL_DISTANCE] << " (" << SEQ_label[SEQL_ALIGNMENT_LENGTH]
                     << "): " << cumul_distance[length[i]][length[j]] / (length[i] + length[j] - nb_begin_end)
                     << " (" << path_length[length[i]][length[j]] - nb_begin_end << ") = "
                     << deletion_distance / (length[i] + length[j] - nb_begin_end) << " (" << nb_deletion << " d) + "
                     << insertion_distance / (length[i] + length[j] - nb_begin_end) << " (" << nb_insertion << " i) + 0 ("
                     << nb_match << " m) + " << substitution_distance / (length[i] + length[j] - nb_begin_end)
                     << " (" << nb_substitution << " s)";
                if (transposition_flag) {
                  cout << " + " << transposition_distance / (length[i] + length[j] - nb_begin_end)
                       << " (" << nb_transposition << " t)";
                }
                cout << endl;

                cout << SEQ_label[SEQL_MAX_GAP_LENGTH] << ": " << max_gap_length << endl;
              }
            }

            if ((out_file) && (nb_alignment <= FILE_NB_ALIGNMENT)) {
              switch (result_format) {

              case ASCII : {
                alignment_ascii_print(*out_file , width , i , j , *alignment , alignment_index);

                if (length[i] + length[j] - nb_begin_end > 0) {
                  *out_file << STAT_label[STATL_DISTANCE] << " (" << SEQ_label[SEQL_ALIGNMENT_LENGTH]
                            << "): " << cumul_distance[length[i]][length[j]] / (length[i] + length[j] - nb_begin_end)
                            << " (" << path_length[length[i]][length[j]] - nb_begin_end << ") = "
                            << deletion_distance / (length[i] + length[j] - nb_begin_end) << " (" << nb_deletion << " d) + "
                            << insertion_distance / (length[i] + length[j] - nb_begin_end) << " (" << nb_insertion << " i) + 0 ("
                            << nb_match << " m) + " << substitution_distance / (length[i] + length[j] - nb_begin_end)
                            << " (" << nb_substitution << " s)";
                  if (transposition_flag) {
                    *out_file << " + " << transposition_distance / (length[i] + length[j] - nb_begin_end)
                              << " (" << nb_transposition << " t)";
                  }
                  *out_file << endl;

                  *out_file << SEQ_label[SEQL_MAX_GAP_LENGTH] << ": " << max_gap_length << endl;
                }
                break;
              }

              case SPREADSHEET : {
                alignment_spreadsheet_print(*out_file , i , j , *alignment , alignment_index);

                if (length[i] + length[j] - nb_begin_end > 0) {
                  *out_file << STAT_label[STATL_DISTANCE] << " (" << SEQ_label[SEQL_ALIGNMENT_LENGTH]
                            << ")\t" << cumul_distance[length[i]][length[j]] / (length[i] + length[j] - nb_begin_end)
                            << " (" << path_length[length[i]][length[j]] - nb_begin_end << ")\t"
                            << deletion_distance / (length[i] + length[j] - nb_begin_end) << " (" << nb_deletion << " d)\t"
                            << insertion_distance / (length[i] + length[j] - nb_begin_end) << " (" << nb_insertion << " i)\t0 ("
                            << nb_match << " m)\t" << substitution_distance / (length[i] + length[j] - nb_begin_end)
                            << " (" << nb_substitution << " s)";
                  if (transposition_flag) {
                    *out_file << "\t" << transposition_distance / (length[i] + length[j] - nb_begin_end)
                              << " (" << nb_transposition << " t)";
                  }
                  *out_file << endl;

                  *out_file << SEQ_label[SEQL_MAX_GAP_LENGTH] << "\t" << max_gap_length << endl;
                }
                break;
              }
              }
            }

            if (!alignment_path.empty()) {
              alignment_index++;
            }
          }
        }
      }
    }

#   ifdef DEBUG
    cout << "\nlocal distance: "
         << mean_local_substitution_distance / nb_local_substitution_distance << endl;
#   endif

    if ((ref_identifier == I_DEFAULT) || (test_identifier == I_DEFAULT)) {

      // writing of distances, alignment lengths, and numbers of deletions, insertions,
      // matchs, substitutions and transpositions

      if ((display) && (nb_alignment <= DISPLAY_NB_ALIGNMENT)) {
        cout << "\n";
        dist_matrix->ascii_write(cout);
      }

      if (out_file) {
        *out_file << "\n";

        switch (result_format) {
        case ASCII :
          dist_matrix->ascii_write(*out_file);
          break;
        case SPREADSHEET :
          dist_matrix->spreadsheet_write(*out_file);
          break;
        }
      }
    }

    if (!alignment_path.empty()) {

      // grouping of insertions/deletions and removing of insertions/deletions corresponding to
      // begin of begin-free alignment or end of end-free alignment

      category = new int[(transposition_flag ? TRANSPOSITION : SUBSTITUTION) + 1];
      category[DELETION] = 0;
      category[INSERTION] = 0;
      category[MATCH] = 1;
      category[SUBSTITUTION] = 2;
      if (transposition_flag) {
        category[TRANSPOSITION] = 3;
      }

      if (display) {
        cout << "\n" << SEQ_label[SEQL_ALIGNMENT_CODING] << "\n" << STAT_label[STATL_INDEL] << ": " << 0
             << "\n" << STAT_label[STATL_MATCH] << ": " << 1 << "\n" << STAT_label[STATL_SUBSTITUTION] << ": " << 2;
        if (transposition_flag) {
          cout << "\n" << STAT_label[STATL_TRANSPOSITION] << ": " << 3;
        }
        cout << endl;
      }

      for (i = 0;i < alignment->nb_sequence;i++) {
        calignment = alignment->int_sequence[i][0];
        offset = 0;
        while ((*calignment == BEGIN_END_DELETION) || (*calignment == BEGIN_END_INSERTION)) {
          offset++;
          calignment++;
        }
        palignment = alignment->int_sequence[i][0];
        for (j = offset;j < alignment->length[i];j++) {
          if ((*calignment == BEGIN_END_DELETION) || (*calignment == BEGIN_END_INSERTION)) {
            break;
          }
          *palignment++ = category[*calignment++];
        }
        alignment->length[i] = j - offset;
      }

      delete [] category;

      alignment->min_value_computation(0);
      alignment->max_value_computation(0);
      alignment->build_marginal_frequency_distribution(0);

      alignment->max_length_computation();
      alignment->cumul_length_computation();
      alignment->build_length_frequency_distribution();

      // writing of alignment sequences

      status = alignment->ascii_data_write(error , alignment_path);
      if ((!status) && (display)) {
        cout << error;
      }
    }

    if (out_file) {
      out_file->close();
      delete out_file;
    }

    delete vector_dist;

    for (i = 0;i < nb_variable;i++) {
      delete [] rank[i];
      delete [] max_category_distance[i];
    }
    delete [] rank;
    delete [] max_category_distance;

    delete alignment;

    for (i = 0;i < nb_sequence;i++) {
      if ((ref_identifier == I_DEFAULT) || (ref_identifier == identifier[i]) ||
          (test_identifier == I_DEFAULT) || (test_identifier == identifier[i])) {
        delete [] local_indel_distance[i];
      }
    }
    delete [] local_indel_distance;

    for (i = 1;i <= max_length;i++) {
      delete [] local_substitution_distance[i];
    }
    delete [] local_substitution_distance;

    for (i = 0;i <= max_length;i++) {
      delete [] cumul_distance[i];
    }
    delete [] cumul_distance;

    for (i = 0;i <= max_length;i++) {
      delete [] path_length[i];
    }
    delete [] path_length;

    for (i = 0;i <= max_length;i++) {
      for (j = 0;j <= max_length;j++) {
        delete [] back_pointers[i][j];
      }
      delete [] back_pointers[i];
    }
    delete [] back_pointers;
  }

  return dist_matrix;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of a fixed substitution cost between 2 different vectors.
 *
 *  \param[in] ref_index             reference sequence index,
 *  \param[in] test_index            test sequence index,
 *  \param[in] ref_position          position in the reference sequence,
 *  \param[in] test_position         position in the test sequence,
 *  \param[in] substitution_distance substitution cost.
 *
 *  \return                          substitution cost.
 */
/*--------------------------------------------------------------*/

double Sequences::substitution_distance_computation(int ref_index , int test_index , int ref_position ,
                                                    int test_position , double substitution_distance) const

{
  int i;
  double distance = 0.;


  for (i = 0;i < nb_variable;i++) {
    if (int_sequence[ref_index][i][ref_position] != int_sequence[test_index][i][test_position]) {
      distance = substitution_distance;
      break;
    }
  }

  return distance;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Alignment of sequences.
 *
 *  \param[in] error           reference on a StatError object, 
 *  \param[in] display         flag for displaying the alignments,
 *  \param[in] ref_identifier  reference sequence identifier,
 *  \param[in] test_identifier test sequence identifier,
 *  \param[in] begin_free      flag begin-free alignment,
 *  \param[in] end_free        flag end-free alignment,
 *  \param[in] result_path     result file path,
 *  \param[in] result_format   result file format (ASCII/SPREADSHEET),
 *  \param[in] alignment_path  alignment file path.
 *
 *  \return                    DistanceMatrix object.
 */
/*--------------------------------------------------------------*/

DistanceMatrix* Sequences::alignment(StatError &error , bool display , int ref_identifier ,
                                     int test_identifier , bool begin_free , bool end_free ,
                                     const string result_path , output_format result_format ,
                                     const string alignment_path) const

{
  bool status = true , half_matrix;
  int i , j , k , m;
  int nb_alignment , ilength , alignment_index , var , width , ref_position ,
      pref_position , test_position , ptest_position , gap_length , max_gap_length ,
      nb_deletion , nb_insertion , nb_match , nb_begin_end , offset , *palignment ,
      *calignment , *category , **path_length , ***back_pointers;
  double buff , deletion_distance , insertion_distance , substitution_distance ,
         **cumul_distance;
  DistanceMatrix *dist_matrix;
  Sequences *alignment;
  ofstream *out_file;


  dist_matrix = NULL;
  error.init();

  if (nb_sequence < 2) {
    status = false;
    error.update(SEQ_error[SEQR_NB_SEQUENCE]);
  }

  if (((index_param_type == TIME) && (index_interval->variance > 0.)) ||
      (index_param_type == POSITION)) {
    status = false;
    error.update(SEQ_error[SEQR_INDEX_PARAMETER_TYPE]);
  }

  for (i = 0;i < nb_variable;i++) {
    if ((type[i] != INT_VALUE) && (type[i] != STATE)) {
      status = false;
      ostringstream error_message , correction_message;
      error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                    << STAT_error[STATR_VARIABLE_TYPE];
      correction_message << STAT_variable_word[INT_VALUE] << " or "
                         << STAT_variable_word[STATE];
      error.correction_update((error_message.str()).c_str() , (correction_message.str()).c_str());
    }
  }

  if (ref_identifier != I_DEFAULT) {
    for (i = 0;i < nb_sequence;i++) {
      if (ref_identifier == identifier[i]) {
        break;
      }
    }

    if (i == nb_sequence) {
      status = false;
      error.update(SEQ_error[SEQR_REF_SEQUENCE_IDENTIFIER]);
    }
  }

  if (test_identifier != I_DEFAULT) {
    for (i = 0;i < nb_sequence;i++) {
      if (test_identifier == identifier[i]) {
        break;
      }
    }

    if (i == nb_sequence) {
      status = false;
      error.update(SEQ_error[SEQR_TEST_SEQUENCE_IDENTIFIER]);
    }
  }

  if ((ref_identifier != I_DEFAULT) && (test_identifier != I_DEFAULT) && (ref_identifier == test_identifier)) {
    status = false;
    error.correction_update(SEQ_error[SEQR_SEQUENCE_IDENTIFIERS] , "different");
  }

  nb_alignment = 1;
  if (ref_identifier == I_DEFAULT) {
    nb_alignment *= (nb_sequence - 1);
  }
  if (test_identifier == I_DEFAULT) {
    nb_alignment *= (ref_identifier != I_DEFAULT ? nb_sequence - 1 : nb_sequence);
  }

  if (nb_alignment > NB_ALIGNMENT) {
    status = false;
    error.update(SEQ_error[SEQR_NB_ALIGNMENT]);
  }

  if (status) {
    if ((ref_identifier == I_DEFAULT) && (test_identifier == I_DEFAULT) && ((result_path.empty()) ||
         (nb_alignment > FILE_NB_ALIGNMENT)) && (alignment_path.empty())) {
      half_matrix = true;
    }
    else {
      half_matrix = false;
    }

    substitution_distance = INDEL_DISTANCE * 2.1;

    if (index_parameter) {
      width = column_width(index_parameter_distribution->nb_value - 1);
    }
    else {
      width = 0;
    }
    for (i = 0;i < nb_variable;i++) {
      var = column_width((int)min_value[i] , (int)max_value[i]);
      if (var > width) {
        width = var;
      }
    }

    out_file = NULL;

    if (!result_path.empty()) {
      out_file = new ofstream(result_path.c_str());

      if (!out_file) {
        error.update(STAT_error[STATR_FILE_NAME]);
        if (display) {
          cout << error;
        }
      }
    }

    // construction of the result data structures

    dist_matrix = new DistanceMatrix(nb_sequence , ref_identifier , test_identifier ,
                                     SEQ_label[SEQL_SEQUENCE] , identifier , false);

    if (!alignment_path.empty()) {
      alignment = new Sequences(nb_alignment , 1);
    }
    else {
      ilength = max_length + max_length;
      alignment = new Sequences(1 , NULL , &ilength , 1 , false);
    }

    // construction of the algorithm data structures

    cumul_distance = new double*[max_length + 1];
    for (i = 0;i <= max_length;i++) {
      cumul_distance[i] = new double[max_length + 1];
    }

    path_length = new int*[max_length + 1];
    for (i = 0;i <= max_length;i++) {
      path_length[i] = new int[max_length + 1];
    }

    back_pointers = new int**[max_length + 1];
    for (i = 0;i <= max_length;i++) {
      back_pointers[i] = new int*[max_length + 1];
      for (j = 0;j <= max_length;j++) {
        back_pointers[i][j] = new int[2];
      }
    }

    // alignment of sequences

    alignment_index = 0;

    for (i = 0;i < nb_sequence;i++) {
      if ((ref_identifier == I_DEFAULT) || (ref_identifier == identifier[i])) {
        for (j = (half_matrix ? i + 1 : 0);j < nb_sequence;j++) {
          if (((test_identifier == I_DEFAULT) || (test_identifier == identifier[j])) && (j != i)) {

            // initialization of the cumulative distances and the corresponding alignment lengths

            cumul_distance[0][0] = 0.;
            path_length[0][0] = 0;

            for (k = 1;k <= length[i];k++) {

              // deletion

              if (begin_free) {
                cumul_distance[k][0] = cumul_distance[k - 1][0];
              }
              else {
                cumul_distance[k][0] = cumul_distance[k - 1][0] + INDEL_DISTANCE;
              }

              path_length[k][0] = k;
              back_pointers[k][0][0] = k - 1;
              back_pointers[k][0][1] = 0;
            }

            for (k = 1;k <= length[j];k++) {

              // insertion

              if (begin_free) {
                cumul_distance[0][k] = cumul_distance[0][k - 1];
              }
              else {
                cumul_distance[0][k] = cumul_distance[0][k - 1] + INDEL_DISTANCE;
              }

              path_length[0][k] = k;
              back_pointers[0][k][0] = 0;
              back_pointers[0][k][1] = k - 1;
            }

            // computation of the cumulative distances and the corresponding alignment lengths

            for (k = 1;k <= length[i];k++) {
              for (m = 1;m <= length[j];m++) {

                // match/substitution

                cumul_distance[k][m] = cumul_distance[k - 1][m - 1] +
                                       substitution_distance_computation(i , j , k - 1 , m - 1 , substitution_distance);
                path_length[k][m] = path_length[k - 1][m - 1] + 1;
                back_pointers[k][m][0] = k - 1;
                back_pointers[k][m][1] = m - 1;

                // deletion

                if ((m < length[j]) || (!end_free)) {
                  buff = cumul_distance[k - 1][m] + INDEL_DISTANCE;
                }
                else {
                  buff = cumul_distance[k - 1][m];
                }

                if (buff < cumul_distance[k][m]) {
                  cumul_distance[k][m] = buff;
                  path_length[k][m] = path_length[k - 1][m] + 1;
                  back_pointers[k][m][0] = k - 1;
                  back_pointers[k][m][1] = m;
                }

                // insertion

                if ((k < length[i]) || (!end_free)) {
                  buff = cumul_distance[k][m - 1] + INDEL_DISTANCE;
                }
                else {
                  buff = cumul_distance[k][m - 1];
                }

                if (buff < cumul_distance[k][m]) {
                  cumul_distance[k][m] = buff;
                  path_length[k][m] = path_length[k][m - 1] + 1;
                  back_pointers[k][m][0] = k;
                  back_pointers[k][m][1] = m - 1;
                }
              }
            }

#           ifdef DEBUG
/*            cout << "\n";
            for (k = length[i];k >= 0;k--) {
              for (m = 0;m <= length[j];m++) {
                cout << cumul_distance[k][m] << "  ";
              }
              cout << endl;
            }
            cout << endl; */
#           endif

            alignment->length[alignment_index] = path_length[length[i]][length[j]];
            if (!alignment_path.empty()) {
              alignment->int_sequence[alignment_index][0] = new int[alignment->length[alignment_index]];
            }

            // backtracking

            deletion_distance = 0.;
            insertion_distance = 0.;

            palignment = alignment->int_sequence[alignment_index][0] + alignment->length[alignment_index];
            pref_position = length[i];
            ptest_position = length[j];

#           ifdef DEBUG
//            cout << pref_position << " " << ptest_position << endl;
#           endif

            for (k = path_length[length[i]][length[j]];k > 0;k--) {
              ref_position = pref_position;
              test_position = ptest_position;
              pref_position = back_pointers[ref_position][test_position][0];
              ptest_position = back_pointers[ref_position][test_position][1];

#             ifdef DEBUG
//              cout << pref_position << " " << ptest_position << endl;
#             endif

              if (test_position == ptest_position) {
                if (((test_position > 0) || (!begin_free)) && ((test_position < length[j]) || (!end_free))) {
                  *--palignment = DELETION;
                  deletion_distance += INDEL_DISTANCE;
                }
                else {
                  *--palignment = BEGIN_END_DELETION;
                }
              }

              else if (ref_position == pref_position) {
                if (((ref_position > 0) || (!begin_free)) && ((ref_position < length[i]) || (!end_free))) {
                  *--palignment = INSERTION;
                  insertion_distance += INDEL_DISTANCE;
                }
                else {
                  *--palignment = BEGIN_END_INSERTION;
                }
              }

              else if ((ref_position == pref_position + 1) && (test_position == ptest_position + 1)) {
                *--palignment = MATCH;
              }
            }

            // search for the maximum number of successive insertions/deletions

            palignment = alignment->int_sequence[alignment_index][0];
            max_gap_length = 0;
            gap_length = 0;

            if ((*palignment == DELETION) || (*palignment == INSERTION)) {
              gap_length++;
            }

            for (k = 1;k < alignment->length[alignment_index];k++) {
              if (*(palignment + 1) != *palignment) {
                if (((*palignment == DELETION) || (*palignment == INSERTION)) &&
                    (gap_length > max_gap_length)) {
                  max_gap_length = gap_length;
                }
                gap_length = 0;
              }

              palignment++;
              if (((*palignment == DELETION) || (*palignment == INSERTION)) &&
                  ((gap_length == 0) || (*palignment == *(palignment - 1)))) {
                gap_length++;
              }
            }

            if (((*palignment == DELETION) || (*palignment == INSERTION)) &&
                (gap_length > max_gap_length)) {
              max_gap_length = gap_length;
            }

            // update of the numbers of deletions, insertions and matchs

            palignment = alignment->int_sequence[alignment_index][0];

            nb_deletion = 0;
            nb_insertion = 0;
            nb_match = 0;
            nb_begin_end = 0;

            for (k = 0;k < alignment->length[alignment_index];k++) {
              switch (*palignment++) {
              case DELETION :
                nb_deletion++;
                break;
              case INSERTION :
                nb_insertion++;
                break;
              case MATCH :
                nb_match++;
                break;
              default :
                nb_begin_end++;
                break;
              }
            }

            dist_matrix->update(identifier[i] , identifier[j] , cumul_distance[length[i]][length[j]] ,
                                length[i] + length[j] - nb_begin_end , deletion_distance , nb_deletion ,
                                insertion_distance , nb_insertion , nb_match);

            if (half_matrix) {
              dist_matrix->update(identifier[j] , identifier[i] , cumul_distance[length[i]][length[j]] ,
                                  length[i] + length[j] - nb_begin_end , insertion_distance , nb_insertion ,
                                  deletion_distance , nb_deletion , nb_match);
            }

#           ifdef DEBUG
            {
              double sum = deletion_distance + insertion_distance;
              if ((sum < cumul_distance[length[i]][length[j]] - DOUBLE_ERROR) ||
                  (sum > cumul_distance[length[i]][length[j]] + DOUBLE_ERROR)) {
                cout << "\nERROR: " << SEQ_label[SEQL_SEQUENCE] << " " << j + 1 << " aligned on "
                     << SEQ_label[SEQL_SEQUENCE] << " " << i + 1 << ": " << cumul_distance[length[i]][length[j]]
                     << " | " << sum << endl;
              }
            }
#           endif

            // writing of the alignment

            if ((display) && (nb_alignment <= DISPLAY_NB_ALIGNMENT)) {
              alignment_ascii_print(cout , width , i , j , *alignment , alignment_index);

              if (length[i] + length[j] - nb_begin_end > 0) {
                cout << STAT_label[STATL_DISTANCE] << " (" << SEQ_label[SEQL_ALIGNMENT_LENGTH]
                     << "): " << cumul_distance[length[i]][length[j]] / (length[i] + length[j] - nb_begin_end)
                     << " (" << path_length[length[i]][length[j]] - nb_begin_end << ") = "
                     << deletion_distance / (length[i] + length[j] - nb_begin_end) << " (" << nb_deletion << " d) + "
                     << insertion_distance / (length[i] + length[j] - nb_begin_end) << " (" << nb_insertion << " i) + 0 ("
                     << nb_match << " m)" << endl;

                cout << SEQ_label[SEQL_MAX_GAP_LENGTH] << ": " << max_gap_length << endl;
              }
            }

            if ((out_file) && (nb_alignment <= FILE_NB_ALIGNMENT)) {
              switch (result_format) {

              case ASCII : {
                alignment_ascii_print(*out_file , width , i , j , *alignment , alignment_index);

                if (length[i] + length[j] - nb_begin_end > 0) {
                  *out_file << STAT_label[STATL_DISTANCE] << " (" << SEQ_label[SEQL_ALIGNMENT_LENGTH]
                            << "): " << cumul_distance[length[i]][length[j]] / (length[i] + length[j] - nb_begin_end)
                            << " (" << path_length[length[i]][length[j]] - nb_begin_end << ") = "
                            << deletion_distance / (length[i] + length[j] - nb_begin_end) << " (" << nb_deletion << " d) + "
                            << insertion_distance / (length[i] + length[j] - nb_begin_end) << " (" << nb_insertion << " i) + 0 ("
                            << nb_match << " m)" << endl;

                  *out_file << SEQ_label[SEQL_MAX_GAP_LENGTH] << ": " << max_gap_length << endl;
                }
                break;
              }

              case SPREADSHEET : {
                alignment_spreadsheet_print(*out_file , i , j , *alignment , alignment_index);

                if (length[i] + length[j] - nb_begin_end > 0) {
                  *out_file << STAT_label[STATL_DISTANCE] << " (" << SEQ_label[SEQL_ALIGNMENT_LENGTH]
                            << ")\t" << cumul_distance[length[i]][length[j]] / (length[i] + length[j] - nb_begin_end)
                            << " (" << path_length[length[i]][length[j]] - nb_begin_end << ")\t"
                            << deletion_distance / (length[i] + length[j] - nb_begin_end) << " (" << nb_deletion << " d)\t"
                            << insertion_distance / (length[i] + length[j] - nb_begin_end) << " (" << nb_insertion << " i)\t0 ("
                            << nb_match << " m)" << endl;

                  *out_file << SEQ_label[SEQL_MAX_GAP_LENGTH] << "\t" << max_gap_length << endl;
                }
                break;
              }
              }
            }

            if (!alignment_path.empty()) {
              alignment_index++;
            }
          }
        }
      }
    }

    if ((ref_identifier == I_DEFAULT) || (test_identifier == I_DEFAULT)) {

      // writing of distances, alignment lengths, and numbers of deletions, insertions and matchs

      if ((display) && (nb_alignment <= DISPLAY_NB_ALIGNMENT)) {
        cout << "\n";
        dist_matrix->ascii_write(cout);
      }

      if (out_file) {
        *out_file << "\n";

        switch (result_format) {
        case ASCII :
          dist_matrix->ascii_write(*out_file);
          break;
        case SPREADSHEET :
          dist_matrix->spreadsheet_write(*out_file);
          break;
        }
      }
    }

    if (!alignment_path.empty()) {

      // grouping of insertions/deletions and removing of insertions/deletions corresponding
      // to begin of begin-free alignment or end of end-free alignment

      category = new int[MATCH + 1];
      category[DELETION] = 0;
      category[INSERTION] = 0;
      category[MATCH] = 1;

      if (display) {
        cout << "\n" << SEQ_label[SEQL_ALIGNMENT_CODING] << "\n" << STAT_label[STATL_INDEL] << ": " << 0
             << "\n" << STAT_label[STATL_MATCH] << ": " << 1 << endl;
      }

      for (i = 0;i < alignment->nb_sequence;i++) {
        calignment = alignment->int_sequence[i][0];
        offset = 0;
        while ((*calignment == BEGIN_END_DELETION) || (*calignment == BEGIN_END_INSERTION)) {
          offset++;
          calignment++;
        }
        palignment = alignment->int_sequence[i][0];
        for (j = offset;j < alignment->length[i];j++) {
          if ((*calignment == BEGIN_END_DELETION) || (*calignment == BEGIN_END_INSERTION)) {
            break;
          }
          *palignment++ = category[*calignment++];
        }
        alignment->length[i] = j - offset;
      }

      delete [] category;

      alignment->min_value_computation(0);
      alignment->max_value_computation(0);
      alignment->build_marginal_frequency_distribution(0);

      alignment->max_length_computation();
      alignment->cumul_length_computation();
      alignment->build_length_frequency_distribution();

      // writing of alignment sequences

      status = alignment->ascii_data_write(error , alignment_path);
      if ((!status) && (display)) {
        cout << error;
      }
    }

    if (out_file) {
      out_file->close();
      delete out_file;
    }

    delete alignment;

    for (i = 0;i <= max_length;i++) {
      delete [] cumul_distance[i];
    }
    delete [] cumul_distance;

    for (i = 0;i <= max_length;i++) {
      delete [] path_length[i];
    }
    delete [] path_length;

    for (i = 0;i <= max_length;i++) {
      for (j = 0;j <= max_length;j++) {
        delete [] back_pointers[i][j];
      }
      delete [] back_pointers[i];
    }
    delete [] back_pointers;
  }

  return dist_matrix;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a multiple alignment of sequences.
 *
 *  \param[in,out] os stream.
 */
/*--------------------------------------------------------------*/

ostream& Sequences::multiple_alignment_ascii_print(ostream &os) const

{
  int i , j , k , m;
  int var , width , rank;
  ios_base::fmtflags format_flags;


  format_flags = os.setf(ios::right , ios::adjustfield);

  width = 0;
  for (i = 0;i < nb_variable - 1;i++) {
    var = column_width((int)min_value[i] , (int)max_value[i]);
    if (var > width) {
      width = var;
    }
  }

  rank = 0;
  for (i = 0;i < max_length;i++) {
    if ((int_sequence[0][nb_variable - 1][i] != GAP) &&
        (int_sequence[0][nb_variable - 1][i] != BEGIN_END_GAP)) {
      os << setw(width) << int_sequence[0][0][i];
    }
    else if (int_sequence[0][nb_variable - 1][i] == GAP) {
      os << setw(width) << "-";
    }
    else {
      os << setw(width) << " ";
    }
    os << " ";

    if (((i - rank) * (width + 1) > LINE_NB_CHARACTER) || (i == max_length - 1)) {
      if (i < max_length - 1) {
        os << "\\";
      }
      else {
        os << "   (" << identifier[0] << ")";
      }
      os << endl;

      for (j = 1;j < nb_variable - 1;j++) {
        for (k = rank;k <= i;k++) {
          if ((int_sequence[0][nb_variable - 1][k] != GAP) &&
              (int_sequence[0][nb_variable - 1][k] != BEGIN_END_GAP)) {
            os << setw(width) << int_sequence[0][j][k];
          }
          else if (int_sequence[0][nb_variable - 1][k] == GAP) {
            os << setw(width) << "-";
          }
          else {
            os << setw(width) << " ";
          }
          os << " ";
        }

        if (i < max_length - 1) {
          os << "\\";
        }
        os << endl;
      }
      os << endl;

      for (j = 1;j < nb_sequence;j++) {
        for (k = 0;k < nb_variable - 1;k++) {
          for (m = rank;m <= i;m++) {
            if ((int_sequence[j][nb_variable - 1][m] != GAP) &&
                (int_sequence[j][nb_variable - 1][m] != BEGIN_END_GAP)) {
              os << setw(width) << int_sequence[j][k][m];
            }
            else if (int_sequence[j][nb_variable - 1][m] == GAP) {
              os << setw(width) << "-";
            }
            else {
              os << setw(width) << " ";
            }
            os << " ";
          }

          if (i < max_length - 1) {
            os << "\\";
          }
          else if (k == 0) {
            os << "   (" << identifier[j] << ")";
          }
          os << endl;
        }
        os << endl;
      }
      os << endl;

      // writing of the consensus sequence

      for (j = 0;j < nb_variable - 1;j++) {
        for (k = rank;k <= i;k++) {
          for (m = 0;m < nb_sequence;m++) {
            if ((int_sequence[m][nb_variable - 1][k] == GAP) ||
                (int_sequence[m][nb_variable - 1][k] == BEGIN_END_GAP) ||
                (int_sequence[m][j][k] != int_sequence[0][j][k])) {
              break;
            }
          }

          if (m == nb_sequence) {
            os << setw(width) << int_sequence[0][j][k];
          }
          else {
            os << setw(width) << ".";
          }
          os << " ";
        }

        if (i < max_length - 1) {
          os << "\\";
        }
        else if (j == 0) {
          os << "   " << SEQ_label[SEQL_CONSENSUS];
        }
        os << endl;
      }
      os << "\n" << endl;

      if (i < max_length - 1) {
        rank = i + 1;
      }
    }
  }

  os.setf(format_flags , ios::adjustfield);

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a multiple alignment of sequences in a file.
 *
 *  \param[in] error reference on a StatError object,
 *  \param[in] path  file path.
 *
 *  \return          error status.
 */
/*--------------------------------------------------------------*/

bool Sequences::multiple_alignment_ascii_print(StatError &error , const string path) const

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
    multiple_alignment_ascii_print(out_file);
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Multiple alignment of sequences.
 *
 *  \param[in] test_seq              reference on the test sequences
 *  \param[in] vector_dist           reference on a VectorDistance object,
 *  \param[in] rank                  ranks (for ordinal variables),
 *  \param[in] max_category_distance maximum distances between categories,
 *  \param[in] begin_free            flag begin-free alignment,
 *  \param[in] end_free              flag end-free alignment,
 *  \param[in] indel_cost            insertion/deletion costs adaptative or fixed,
 *  \param[in] indel_factor          factor for deducing the insertion/deletion costs.
 *
 *  \return                          Sequences object.
 */
/*--------------------------------------------------------------*/

Sequences* Sequences::multiple_alignment(const Sequences &test_seq , const VectorDistance &vector_dist ,
                                         double **rank , double **max_category_distance , bool begin_free ,
                                         bool end_free , insertion_deletion_cost indel_cost , double indel_factor) const

{
  int i , j , k , m;
  int ref_position , pref_position , test_position , ptest_position , *alignment , *palignment ,
      *ilength , **path_length , ***back_pointers;
  double buff , sum , **ref_local_indel_distance , **test_local_indel_distance , **cumul_distance;
  Sequences *seq;


  // construction of the algorithm data structures - computation of the insertion/deletion costs

  if (indel_cost == FIXED) {
    buff = indel_distance_computation(vector_dist , rank , max_category_distance) * indel_factor;
  }

  ref_local_indel_distance = new double*[nb_sequence];
  for (i = 0;i < nb_sequence;i++) {
    ref_local_indel_distance[i] = new double[max_length + 1];

    switch (indel_cost) {

    case ADAPTATIVE : {
      for (j = 1;j <= max_length;j++) {
        if (int_sequence[i][nb_variable - 1][j - 1] == DATA) {
          ref_local_indel_distance[i][j] = indel_distance_computation(vector_dist , i , j - 1 , rank , max_category_distance) * indel_factor;
        }
        else {
          ref_local_indel_distance[i][j] = 0.;
        }
      }
      break;
    }

    case FIXED : {
      for (j = 1;j <= max_length;j++) {
        if (int_sequence[i][nb_variable - 1][j - 1] == DATA) {
          ref_local_indel_distance[i][j] = buff;
        }
        else {
          ref_local_indel_distance[i][j] = 0.;
        }
      }
      break;
    }
    }
  }

  test_local_indel_distance = new double*[test_seq.nb_sequence];
  for (i = 0;i < test_seq.nb_sequence;i++) {
    test_local_indel_distance[i] = new double[test_seq.max_length + 1];

    // computation of the insertion/deletion costs

    switch (indel_cost) {

    case ADAPTATIVE : {
      for (j = 1;j <= test_seq.max_length;j++) {
        if (test_seq.int_sequence[i][nb_variable - 1][j - 1] == DATA) {
          test_local_indel_distance[i][j] = test_seq.indel_distance_computation(vector_dist , i , j - 1 , rank , max_category_distance) * indel_factor;
        }
        else {
          test_local_indel_distance[i][j] = 0.;
        }
      }
      break;
    }

    case FIXED : {
      for (j = 1;j <= test_seq.max_length;j++) {
        if (test_seq.int_sequence[i][nb_variable - 1][j - 1] == DATA) {
          test_local_indel_distance[i][j] = buff;
        }
        else {
          test_local_indel_distance[i][j] = 0.;
        }
      }
      break;
    }
    }
  }

  cumul_distance = new double*[max_length + 1];
  for (i = 0;i <= max_length;i++) {
    cumul_distance[i] = new double[test_seq.max_length + 1];
  }

  path_length = new int*[max_length + 1];
  for (i = 0;i <= max_length;i++) {
    path_length[i] = new int[test_seq.max_length + 1];
  }

  back_pointers = new int**[max_length + 1];
  for (i = 0;i <= max_length;i++) {
    back_pointers[i] = new int*[test_seq.max_length + 1];
    for (j = 0;j <= test_seq.max_length;j++) {
      back_pointers[i][j] = new int[2];
    }
  }

  alignment = new int[max_length + test_seq.max_length];

  // initialization of the cumulative distances and the corresponding alignment lengths

  cumul_distance[0][0] = 0.;
  path_length[0][0] = 0;

  for (i = 1;i <= max_length;i++) {

    // deletion

    if (begin_free) {
      cumul_distance[i][0] = cumul_distance[i - 1][0];
    }

    else {
      sum = 0.;
      for (j = 0;j < nb_sequence;j++) {
        if (int_sequence[j][nb_variable - 1][i - 1] == DATA) {
          sum += ref_local_indel_distance[j][i];
        }
      }
      cumul_distance[i][0] = cumul_distance[i - 1][0] + sum / nb_sequence;
    }

    path_length[i][0] = i;
    back_pointers[i][0][0] = i - 1;
    back_pointers[i][0][1] = 0;
  }

  for (i = 1;i <= test_seq.max_length;i++) {

    // insertion

    if (begin_free) {
      cumul_distance[0][i] = cumul_distance[0][i - 1];
    }

    else {
      sum = 0.;
      for (j = 0;j < test_seq.nb_sequence;j++) {
        if (test_seq.int_sequence[j][nb_variable - 1][i - 1] == DATA) {
          sum += test_local_indel_distance[j][i];
        }
      }
      cumul_distance[0][i] = cumul_distance[0][i - 1] + sum / test_seq.nb_sequence;
    }

    path_length[0][i] = i;
    back_pointers[0][i][0] = 0;
    back_pointers[0][i][1] = i - 1;
  }

  // computation of the cumulative distances and the corresponding alignment lengths

  for (i = 1;i <= max_length;i++) {
    for (j = 1;j <= test_seq.max_length;j++) {

      sum = 0.;
      for (k = 0;k < nb_sequence;k++) {
        for (m = 0;m < test_seq.nb_sequence;m++) {

          // computation of the distance of substitution of one vector by another

          if ((int_sequence[k][nb_variable - 1][i - 1] == DATA) && (test_seq.int_sequence[m][nb_variable - 1][j - 1] == GAP)) {
            sum += ref_local_indel_distance[k][i];
          }
          else if ((int_sequence[k][nb_variable - 1][i - 1] == GAP) && (test_seq.int_sequence[m][nb_variable - 1][j - 1] == DATA)) {
            sum += test_local_indel_distance[m][j];
          }
          else if ((int_sequence[k][nb_variable - 1][i - 1] == DATA) && (test_seq.int_sequence[m][nb_variable - 1][j - 1] == DATA)) {
            sum += substitution_distance_computation(vector_dist , k , m , i - 1 , j - 1 , rank , &test_seq);
          }
        }
      }

      // match/substitution

      cumul_distance[i][j] = cumul_distance[i - 1][j - 1] + sum / (nb_sequence * test_seq.nb_sequence);
      path_length[i][j] = path_length[i - 1][j - 1] + 1;
      back_pointers[i][j][0] = i - 1;
      back_pointers[i][j][1] = j - 1;

      // deletion

      if ((j < test_seq.max_length) || (!end_free)) {
        sum = 0.;
        for (k = 0;k < nb_sequence;k++) {
          if (int_sequence[k][nb_variable - 1][i - 1] == DATA) {
            sum += ref_local_indel_distance[k][i];
          }
        }
        buff = cumul_distance[i - 1][j] + sum / nb_sequence;
      }

      else {
        buff = cumul_distance[i - 1][j];
      }

      if (buff < cumul_distance[i][j]) {
        cumul_distance[i][j] = buff;
        path_length[i][j] = path_length[i - 1][j] + 1;
        back_pointers[i][j][0] = i - 1;
        back_pointers[i][j][1] = j;
      }

      // insertion

      if ((i < max_length) || (!end_free)) {
        sum = 0.;
        for (k = 0;k < test_seq.nb_sequence;k++) {
          if (test_seq.int_sequence[k][nb_variable - 1][j - 1] == DATA) {
            sum += test_local_indel_distance[k][j];
          }
        }
        buff = cumul_distance[i][j - 1] + sum / test_seq.nb_sequence;
      }

      else {
        buff = cumul_distance[i][j - 1];
      }

      if (buff < cumul_distance[i][j]) {
        cumul_distance[i][j] = buff;
        path_length[i][j] = path_length[i][j - 1] + 1;
        back_pointers[i][j][0] = i;
        back_pointers[i][j][1] = j - 1;
      }
    }
  }

# ifdef DEBUG
  cout << "\nMultiple alignment distance: " << cumul_distance[max_length][test_seq.max_length] << endl;
# endif

  // backtracking

  palignment = alignment + path_length[max_length][test_seq.max_length];
  pref_position = max_length;
  ptest_position = test_seq.max_length;

# ifdef DEBUG
//   cout << pref_position << " " << ptest_position << endl;
# endif

  for (i = path_length[max_length][test_seq.max_length];i > 0;i--) {
    ref_position = pref_position;
    test_position = ptest_position;
    pref_position = back_pointers[ref_position][test_position][0];
    ptest_position = back_pointers[ref_position][test_position][1];

#   ifdef DEBUG
//    cout << pref_position << " " << ptest_position << endl;
#   endif

    if (test_position == ptest_position) {
      if (((test_position > 0) || (!begin_free)) && ((test_position < test_seq.max_length) || (!end_free))) {
        *--palignment = DELETION;
      }
      else {
        *--palignment = BEGIN_END_DELETION;
      }
    }

    else if (ref_position == pref_position) {
      if (((ref_position > 0) || (!begin_free)) && ((ref_position < max_length) || (!end_free))) {
        *--palignment = INSERTION;
      }
      else {
        *--palignment = BEGIN_END_INSERTION;
      }
    }

    else if ((ref_position == pref_position + 1) && (test_position == ptest_position + 1)) {
      *--palignment = SUBSTITUTION;
    }
  }

  // construction of the group of sequences

  ilength = new int[nb_sequence + test_seq.nb_sequence];
  for (i = 0;i < nb_sequence + test_seq.nb_sequence;i++) {
    ilength[i] = path_length[max_length][test_seq.max_length];
  }
  seq = new Sequences(nb_sequence + test_seq.nb_sequence , NULL , ilength ,
                      nb_variable , type);
  delete [] ilength;

  for (i = 0;i < nb_variable;i++) {
    seq->min_value[i] = min_value[i];
    seq->max_value[i] = max_value[i];
  }

  for (i = 0;i < nb_sequence;i++) {
    seq->identifier[i] = identifier[i];
  }
  for (i = 0;i < test_seq.nb_sequence;i++) {
    seq->identifier[nb_sequence + i] = test_seq.identifier[i];
  }

  palignment = alignment;
  ref_position = 0;
  test_position = 0;

  for (i = 0;i < path_length[max_length][test_seq.max_length];i++) {
    if ((*palignment != INSERTION) && (*palignment != BEGIN_END_INSERTION)) {
      for (j = 0;j < nb_sequence;j++) {
        for (k = 0;k < nb_variable;k++) {
          seq->int_sequence[j][k][i] = int_sequence[j][k][ref_position];
        }
      }
      ref_position++;
    }

    else {
      for (j = 0;j < nb_sequence;j++) {
        for (k = 0;k < nb_variable - 1;k++) {
          seq->int_sequence[j][k][i] = (int)max_value[k] + 1;
        }

        switch (*palignment) {

        case INSERTION : {
          if (((i > 0) && (int_sequence[j][nb_variable - 1][ref_position - 1] == BEGIN_END_GAP)) ||
              ((i < path_length[max_length][test_seq.max_length] - 1) && (int_sequence[j][nb_variable - 1][ref_position + 1] == BEGIN_END_GAP))) {
            seq->int_sequence[j][nb_variable - 1][i] = BEGIN_END_GAP;
          }
          else {
            seq->int_sequence[j][nb_variable - 1][i] = GAP;
          }
          break;
        }

        case BEGIN_END_INSERTION : {
          seq->int_sequence[j][nb_variable - 1][i] = BEGIN_END_GAP;
          break;
        }
        }
      }
    }

    if ((*palignment != DELETION) && (*palignment != BEGIN_END_DELETION)) {
      for (j = 0;j < test_seq.nb_sequence;j++) {
        for (k = 0;k < nb_variable;k++) {
          seq->int_sequence[j + nb_sequence][k][i] = test_seq.int_sequence[j][k][test_position];
        }
      }
      test_position++;
    }

    else {
      for (j = 0;j < test_seq.nb_sequence;j++) {
        for (k = 0;k < nb_variable - 1;k++) {
          seq->int_sequence[j + nb_sequence][k][i] = (int)max_value[k] + 1;
        }

        switch (*palignment) {

        case DELETION : {
          if (((i > 0) && (test_seq.int_sequence[j][nb_variable - 1][test_position - 1] == BEGIN_END_GAP)) ||
              ((i < path_length[max_length][test_seq.max_length] - 1) && (test_seq.int_sequence[j][nb_variable - 1][test_position + 1] == BEGIN_END_GAP))) {
            seq->int_sequence[j + nb_sequence][nb_variable - 1][i] = BEGIN_END_GAP;
          }
          else {
            seq->int_sequence[j + nb_sequence][nb_variable - 1][i] = GAP;
          }
          break;
        }

        case BEGIN_END_DELETION : {
          seq->int_sequence[j + nb_sequence][nb_variable - 1][i] = BEGIN_END_GAP;
          break;
        }
        }
      }
    }

    palignment++;
  }

  for (i = 0;i < nb_sequence;i++) {
    delete [] ref_local_indel_distance[i];
  }
  delete [] ref_local_indel_distance;

  for (i = 0;i < test_seq.nb_sequence;i++) {
    delete [] test_local_indel_distance[i];
  }
  delete [] test_local_indel_distance;

  for (i = 0;i <= max_length;i++) {
    delete [] cumul_distance[i];
  }
  delete [] cumul_distance;

  for (i = 0;i <= max_length;i++) {
    delete [] path_length[i];
  }
  delete [] path_length;

  for (i = 0;i <= max_length;i++) {
    for (j = 0;j <= test_seq.max_length;j++) {
      delete [] back_pointers[i][j];
    }
    delete [] back_pointers[i];
  }
  delete [] back_pointers;

  delete [] alignment;

  return seq;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Multiple alignment of sequences.
 *
 *  \param[in] error        reference on a StatError object,
 *  \param[in] display      flag for displaying the alignment,
 *  \param[in] ivector_dist reference on a VectorDistance object,
 *  \param[in] begin_free   flag begin-free alignment,
 *  \param[in] end_free     flag end-free alignment,
 *  \param[in] indel_cost   insertion/deletion costs adaptative or fixed,
 *  \param[in] indel_factor factor for deducing the insertion/deletion costs,
 *  \param[in] strategy     type of algorithm for the dendrogram construction,
 *  \param[in] path         file path.
 *
 *  \return                 Sequences object.
 */
/*--------------------------------------------------------------*/

Sequences* Sequences::multiple_alignment(StatError &error , bool display ,
                                         const VectorDistance &ivector_dist ,
                                         bool begin_free , bool end_free ,
                                         insertion_deletion_cost indel_cost , double indel_factor ,
                                         hierarchical_strategy strategy ,
                                         const string path) const

{
  bool status = true;
  int i , j , k;
  int *itype , *psequence , *csequence , *variable;
  double **rank , **max_category_distance;
  VectorDistance *vector_dist;
  DistanceMatrix *dist_matrix;
  Dendrogram *dendrogram;
  Sequences *seq , **clustered_seq;


  seq = NULL;
  error.init();

  if (index_parameter) {
    status = false;
    error.update(SEQ_error[SEQR_INDEX_PARAMETER_TYPE]);
  }

  for (i = 0;i < nb_variable;i++) {
    if ((type[i] != INT_VALUE) && (type[i] != STATE)) {
      status = false;
      ostringstream error_message , correction_message;
      error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                    << STAT_error[STATR_VARIABLE_TYPE];
      correction_message << STAT_variable_word[INT_VALUE] << " or "
                         << STAT_variable_word[STATE];
      error.correction_update((error_message.str()).c_str() , (correction_message.str()).c_str());
    }
  }

  if (status) {

    // pairwise alignment of sequences

    dist_matrix = alignment(error , NULL , ivector_dist , I_DEFAULT , I_DEFAULT ,
                            begin_free , end_free , indel_cost , indel_factor);

    if (dist_matrix) {

      // construction of a dendrogram on the basis of the matrix of pairwise alignment of sequences

      if (strategy != DIVISIVE) {
        dendrogram = dist_matrix->agglomerative_hierarchical_clustering(strategy);
      }
      else {
        dendrogram = dist_matrix->divisive_hierarchical_clustering();
      }

      if (display) {
        cout << *dendrogram << "\n";
      }

      vector_dist = new VectorDistance(ivector_dist);

      // computation of the maximum substitution distance for nominal variables and
      // the ranks for ordinal variables

      rank = new double*[nb_variable];
      max_category_distance = new double*[nb_variable];

      for (i = 0;i < nb_variable;i++) {
        if ((vector_dist->get_var_type(i) == NOMINAL) && (vector_dist->get_category_distance(i))) {
          max_category_distance[i] = vector_dist->max_category_distance_computation(i);
        }
        else {
          max_category_distance[i] = 0;
        }

        if (vector_dist->get_var_type(i) == ORDINAL) {
          rank[i] = marginal_distribution[i]->rank_computation();
        }
        else {
          rank[i] = NULL;
        }

        // computation of the dispersion measures for standardization

        if (marginal_distribution[i]) {
          vector_dist->dispersion_computation(i , marginal_distribution[i] , rank[i]);
        }

        else {
          switch (vector_dist->get_distance_type()) {
          case ABSOLUTE_VALUE :
            vector_dist->dispersion_update(i , mean_absolute_difference_computation(i));
            break;
          case QUADRATIC :
            vector_dist->dispersion_update(i ,  2 * variance_computation(i , mean_computation(i)));
            break;
          }

          if (vector_dist->get_dispersion(i) == 0.) {
            vector_dist->dispersion_update(i , 1.);
          }
        }
      }

      // construction of the initial groups

      itype = new int[nb_variable + 1];
      for (i = 0;i < nb_variable;i++) {
        itype[i] = type[i];
      }
      itype[nb_variable] = INT_VALUE;

      clustered_seq = new Sequences*[2 * nb_sequence - 1];
      for (i = 0;i < nb_sequence;i++) {
        clustered_seq[i] = new Sequences(1 , &identifier[i] , &length[i] ,
                                         nb_variable + 1 , itype);

        for (j = 0;j < nb_variable;j++) {
          clustered_seq[i]->min_value[j] = min_value[j];
          clustered_seq[i]->max_value[j] = max_value[j];

          psequence = clustered_seq[i]->int_sequence[0][j];
          csequence = int_sequence[i][j];
          for (k = 0;k < length[i];k++) {
            *psequence++ = *csequence++;
          }
        }

        psequence = clustered_seq[i]->int_sequence[0][nb_variable];
        for (j = 0;j < length[i];j++) {
          *psequence++ = DATA;
        }
      }

      // multiple alignment of sequences

      for (i = nb_sequence;i < 2 * nb_sequence - 1;i++) {
        clustered_seq[i] = clustered_seq[dendrogram->get_child(i , 0)]->multiple_alignment(*(clustered_seq[dendrogram->get_child(i , 1)]) ,
                                                                                           *vector_dist , rank ,
                                                                                           max_category_distance , begin_free ,
                                                                                           end_free , indel_cost , indel_factor);

#       ifdef DEBUG
        if (i < 2 * nb_sequence - 2) {
          clustered_seq[i]->multiple_alignment_ascii_print(cout);
        }
#       endif

      }

      // writing of the multiple alignment

      if (display) {
        clustered_seq[2 * nb_sequence - 2]->multiple_alignment_ascii_print(cout);
      }

      if (!path.empty()) {
        status = clustered_seq[2 * nb_sequence - 2]->multiple_alignment_ascii_print(error , path);
        if ((!status) && (display)) {
          cout << error;
        }
      }

      seq = new Sequences(nb_sequence , clustered_seq[2 * nb_sequence - 2]->identifier ,
                          clustered_seq[2 * nb_sequence - 2]->length , nb_variable , itype);

      variable = new int[nb_variable];
      for (i = 0;i < nb_variable;i++) {
        variable[i] = i;
      }
      seq->select_variable(*(clustered_seq[2 * nb_sequence - 2]) , variable);

      for (i = 0;i < nb_variable;i++) {
        (seq->max_value[i])++;
        seq->build_marginal_frequency_distribution(i);
      }

      delete [] variable;

      for (i = 0;i < 2 * nb_sequence - 1;i++) {
        delete clustered_seq[i];
      }
      delete [] clustered_seq;

      delete [] itype;

      delete dendrogram;

      delete vector_dist;

      for (i = 0;i < nb_variable;i++) {
        delete [] rank[i];
        delete [] max_category_distance[i];
      }
      delete [] rank;
      delete [] max_category_distance;
    }

    delete dist_matrix;
  }

  return seq;
}


};  // namespace sequence_analysis
