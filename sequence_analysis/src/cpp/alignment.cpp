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
 *       $Id: alignment.cpp 18042 2015-04-23 09:31:47Z guedon $
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
#include <iostream>

#include "stat_tool/stat_tools.h"
#include "stat_tool/curves.h"
#include "stat_tool/distribution.h"
#include "stat_tool/markovian.h"
#include "stat_tool/vectors.h"
#include "stat_tool/distance_matrix.h"
#include "stat_tool/stat_label.h"

#include "sequences.h"
#include "sequence_label.h"
#include "tool/config.h"

using namespace std;
using namespace stat_tool;


namespace sequence_analysis {



/*--------------------------------------------------------------*
 *
 *  Ecriture de l'alignement de 2 sequences.
 *
 *  arguments : stream, largeur des colonnes, indices des 2 sequences,
 *              reference sur l'alignement, indice de l'alignement.
 *
 *--------------------------------------------------------------*/

ostream& Sequences::alignment_ascii_print(ostream &os , int width , int ref_index , int test_index ,
                                          const Sequences &alignment , int alignment_index) const

{
  register int i , j , k , m , n;
  int ref_rank , test_rank , alignment_rank;
  long old_adjust;


  old_adjust = os.setf(ios::right , ios::adjustfield);

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

      // sequence de test

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

      // operations

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

      // sequence de reference

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

  os.setf((FMTFLAGS)old_adjust , ios::adjustfield);
  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture de l'alignement de 2 sequences au format tableur.
 *
 *  arguments : stream, indices des 2 sequences,
 *              reference sur l'alignement, indice de l'alignement.
 *
 *--------------------------------------------------------------*/

ostream& Sequences::alignment_spreadsheet_print(ostream &os , int ref_index , int test_index ,
                                                const Sequences &alignment , int alignment_index) const

{
  register int i , j , k;


  os << "\n" << SEQ_label[SEQL_SEQUENCE] << " " << identifier[test_index]
     << "\t" << SEQ_label[SEQL_LENGTH] << " " << length[test_index] << "\t" << SEQ_label[SEQL_ALIGNED_ON]
     << "\t" << SEQ_label[SEQL_SEQUENCE] << " " << identifier[ref_index]
     << "\t" << SEQ_label[SEQL_LENGTH] << " " << length[ref_index] << endl;

  // sequence de test

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

  // operations

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

  // sequence de reference

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


/*--------------------------------------------------------------*
 *
 *  Calcul du cout d'elision/insertion fixe.
 *
 *  arguments : reference sur un objet VectorDistance, rangs,
 *              distances maximum aux symboles.
 *
 *--------------------------------------------------------------*/

double Sequences::indel_distance_computation(const VectorDistance &vector_dist ,
                                             double **rank , double **max_symbol_distance) const

{
  register int i , j;
  double ldistance , distance = 0.;


  for (i = 0;i < vector_dist.get_nb_variable();i++) {
    switch (vector_dist.get_variable_type(i)) {

    case SYMBOLIC : {
      if (!max_symbol_distance[i]) {
        ldistance = 1.;
      }

      else {
        ldistance = 0.;
        for (j = (int)min_value[i];j <= (int)max_value[i];j++) {
          if (max_symbol_distance[i][j] > ldistance) {
            ldistance = max_symbol_distance[i][j];
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


/*--------------------------------------------------------------*
 *
 *  Calcul du cout d'elision/insertion d'un vecteur.
 *
 *  arguments : reference sur un objet VectorDistance, indice de la sequence,
 *              position dans la sequence, rangs, distances maximum aux symboles.
 *
 *--------------------------------------------------------------*/

double Sequences::indel_distance_computation(const VectorDistance &vector_dist ,
                                             int index , int position , double **rank ,
                                             double **max_symbol_distance) const

{
  register int i;
  double ldistance , distance = 0.;


  for (i = 0;i < vector_dist.get_nb_variable();i++) {
    switch (vector_dist.get_variable_type(i)) {

    case SYMBOLIC : {
      if (!max_symbol_distance[i]) {
        ldistance = 1.;
      }
      else {
        ldistance = max_symbol_distance[i][int_sequence[index][i][position]];
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


/*--------------------------------------------------------------*
 *
 *  Calcul du cout de substitution d'un vecteur par un autre vecteur.
 *
 *  arguments : reference sur un objet VectorDistance, indices des 2 sequences,
 *              positions dans les 2 sequence, rangs, pointeur sur
 *              les sequences de test (alignement multiple).
 *
 *--------------------------------------------------------------*/

double Sequences::substitution_distance_computation(const VectorDistance &vector_dist ,
                                                    int ref_index , int test_index ,
                                                    int ref_position , int test_position ,
                                                    double **rank , const Sequences *test_seq) const

{
  register int i;
  double ldistance , distance = 0.;


  if (!test_seq) {
    test_seq = this;
  }

  for (i = 0;i < vector_dist.get_nb_variable();i++) {
    switch (vector_dist.get_variable_type(i)) {

    case SYMBOLIC : {
      if (!vector_dist.get_symbol_distance(i)) {
        ldistance = (int_sequence[ref_index][i][ref_position] == test_seq->int_sequence[test_index][i][test_position] ? 0. : 1.);
      }
      else {
        ldistance = vector_dist.get_symbol_distance(i , int_sequence[ref_index][i][ref_position] ,
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


/*--------------------------------------------------------------*
 *
 *  Alignement de sequences.
 *
 *  arguments : reference sur un objet StatError, stream, reference sur un objet VectorDistance,
 *              identificateurs de 2 sequences, flags debut/fin libres,
 *              couts d'elision/insertion adaptatifs ou fixes,
 *              facteur pour deduire les couts d'elision/insertion,
 *              flag transposition, facteur pour deduire les couts de transposition,
 *              path resultat, format du fichier resultat ('a' : ASCII, 's' : Spreadsheet),
 *              path alignement, format du fichier alignement ('a' : ASCII, 'b' : Binary).
 *
 *--------------------------------------------------------------*/

DistanceMatrix* Sequences::alignment(StatError &error , ostream *os , const VectorDistance &ivector_dist ,
                                     int ref_identifier , int test_identifier , bool begin_free ,
                                     bool end_free , int indel_cost , double indel_factor ,
                                     bool transposition_flag , double transposition_factor ,
                                     const char *result_path , char result_format ,
                                     const char *alignment_path , char alignment_format) const

{
  bool status = true , half_matrix;
  register int i , j , k , m;
  int nb_alignment , ilength , alignment_index , var , width , ref_position , pref_position ,
      test_position , ptest_position , gap_length , max_gap_length , nb_deletion , nb_insertion ,
      nb_match , nb_substitution , nb_transposition , nb_begin_end , offset , *palignment ,
      *calignment , *symbol , **path_length , ***back_pointers;
  double buff , deletion_distance , insertion_distance , substitution_distance ,
         transposition_distance , max_transposition_cost , **rank , **max_symbol_distance ,
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

  if (((index_parameter_type == TIME) && (index_interval->variance > 0.)) ||
      (index_parameter_type == POSITION)) {
    status = false;
    error.update(SEQ_error[SEQR_INDEX_PARAMETER_TYPE]);
  }

  if (ivector_dist.get_nb_variable() != nb_variable) {
    status = false;
    error.update(STAT_error[STATR_NB_VARIABLE]);
  }

  else {
    for (i = 0;i < nb_variable;i++) {
      if (ivector_dist.get_variable_type(i) != NUMERIC) {
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

        if ((ivector_dist.get_variable_type(i) == SYMBOLIC) &&
            ((min_value[i] < 0) || (max_value[i] >= NB_SYMBOL) ||
             ((ivector_dist.get_symbol_distance(i)) && (ivector_dist.get_nb_value(i) != max_value[i] + 1)))) {
          status = false;
          ostringstream error_message;
          error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                        << STAT_error[STATR_NB_SYMBOL];
          error.update((error_message.str()).c_str());
        }

        if ((ivector_dist.get_variable_type(i) == CIRCULAR) &&
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
    if ((ref_identifier == I_DEFAULT) && (test_identifier == I_DEFAULT) && ((!result_path) ||
         (nb_alignment > FILE_NB_ALIGNMENT)) && (!alignment_path)) {
      half_matrix = true;
    }
    else {
      half_matrix = false;
    }

    vector_dist = new VectorDistance(ivector_dist);

    // calcul des distance maximum de substitution pour les variables symboliques et
    // des rangs pour les variables ordinales

    rank = new double*[nb_variable];
    max_symbol_distance = new double*[nb_variable];

    for (i = 0;i < nb_variable;i++) {
      if ((vector_dist->get_variable_type(i) == SYMBOLIC) && (vector_dist->get_symbol_distance(i))) {
        max_symbol_distance[i] = vector_dist->max_symbol_distance_computation(i);
      }
      else {
        max_symbol_distance[i] = 0;
      }

      if (vector_dist->get_variable_type(i) == ORDINAL) {
        rank[i] = marginal_distribution[i]->rank_computation();
      }
      else {
        rank[i] = NULL;
      }

      // calcul des dispersions pour la standardisation

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
    if (vector_dist->distance_type == ABSOLUTE_VALUE) {
      for (i = 0;i < nb_variable;i++) {
        if (vector_dist->variable_type[i] == NUMERIC) {
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

    if (result_path) {
      out_file = new ofstream(result_path);

      if (!out_file) {
        error.update(STAT_error[STATR_FILE_NAME]);

#       ifdef MESSAGE
        if (os) {
          *os << error;
        }
#       endif

      }
    }

    // creation des structures de donnees resultat

    dist_matrix = new DistanceMatrix(nb_sequence , ref_identifier , test_identifier ,
                                     SEQ_label[SEQL_SEQUENCE] , identifier ,
                                     true , transposition_flag);

    if (alignment_path) {
      alignment = new Sequences(nb_alignment , 1);
    }
    else {
      ilength = max_length + max_length;
      alignment = new Sequences(1 , NULL , &ilength , 1 , false);
    }

    // creation des structures de donnees de l'algorithme - calcul des couts d'elision/insertion

    if (indel_cost == FIXED) {
      buff = indel_distance_computation(*vector_dist , rank , max_symbol_distance) * indel_factor;
    }

    local_indel_distance = new double*[nb_sequence];
    for (i = 0;i < nb_sequence;i++) {
      if ((ref_identifier == I_DEFAULT) || (ref_identifier == identifier[i]) ||
          (test_identifier == I_DEFAULT) || (test_identifier == identifier[i])) {
        local_indel_distance[i] = new double[length[i] + 1];

        switch (indel_cost) {

        case ADAPTATIVE : {
          for (j = 1;j <= length[i];j++) {
            local_indel_distance[i][j] = indel_distance_computation(*vector_dist , i , j - 1 , rank , max_symbol_distance) * indel_factor;
          }

#         ifdef DEBUG
/*          cout << "\ncout d'elision/insertion des vecteurs de la sequence " << i << ": ";
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

    // alignement des sequences

#   ifdef DEBUG
    int nb_local_substitution_distance = 0;
    double mean_local_substitution_distance = 0.;
#   endif

    alignment_index = 0;

    for (i = 0;i < nb_sequence;i++) {
      if ((ref_identifier == I_DEFAULT) || (ref_identifier == identifier[i])) {
        for (j = (half_matrix ? i + 1 : 0);j < nb_sequence;j++) {
          if (((test_identifier == I_DEFAULT) || (test_identifier == identifier[j])) && (j != i)) {

            // initialisation des distances cumulees et de la longueur des chemins correspondants

            cumul_distance[0][0] = 0.;
            path_length[0][0] = 0;

            for (k = 1;k <= length[i];k++) {

              // elision

              switch (begin_free) {
              case false :
                cumul_distance[k][0] = cumul_distance[k - 1][0] + local_indel_distance[i][k];
                break;
              case true :
                cumul_distance[k][0] = cumul_distance[k - 1][0];
                break;
              }

              path_length[k][0] = k;
              back_pointers[k][0][0] = k - 1;
              back_pointers[k][0][1] = 0;
            }

            for (k = 1;k <= length[j];k++) {

              // insertion

              switch (begin_free) {
              case false :
                cumul_distance[0][k] = cumul_distance[0][k - 1] + local_indel_distance[j][k];
                break;
              case true :
                cumul_distance[0][k] = cumul_distance[0][k - 1];
                break;
              }

              path_length[0][k] = k;
              back_pointers[0][k][0] = 0;
              back_pointers[0][k][1] = k - 1;
            }

            // calcul des distances cumulees et de la longueur des chemins correspondants

            for (k = 1;k <= length[i];k++) {
              for (m = 1;m <= length[j];m++) {

                // calcul distance locale

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

                // elision

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

            // fin libre (autre implementation)

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
            if (alignment_path) {
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

            // recherche du nombre d'elisions/insertions successives maximum

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

            // mise a jour des nombre d'elisions, d'insertions, de matchs,
            // de substitutions et de transpositions

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

            // ecriture de l'alignement

#           ifdef MESSAGE
            if ((os) && (nb_alignment <= DISPLAY_NB_ALIGNMENT)) {
              alignment_ascii_print(*os , width , i , j , *alignment , alignment_index);

              if (length[i] + length[j] - nb_begin_end > 0) {
                *os << STAT_label[STATL_DISTANCE] << " (" << SEQ_label[SEQL_ALIGNMENT_LENGTH]
                    << "): " << cumul_distance[length[i]][length[j]] / (length[i] + length[j] - nb_begin_end)
                    << " (" << path_length[length[i]][length[j]] - nb_begin_end << ") = "
                    << deletion_distance / (length[i] + length[j] - nb_begin_end) << " (" << nb_deletion << " d) + "
                    << insertion_distance / (length[i] + length[j] - nb_begin_end) << " (" << nb_insertion << " i) + 0 ("
                    << nb_match << " m) + " << substitution_distance / (length[i] + length[j] - nb_begin_end)
                    << " (" << nb_substitution << " s)";
                if (transposition_flag) {
                  *os << " + " << transposition_distance / (length[i] + length[j] - nb_begin_end)
                      << " (" << nb_transposition << " t)";
                }
                *os << endl;

                *os << SEQ_label[SEQL_MAX_GAP_LENGTH] << ": " << max_gap_length << endl;
              }
            }
#           endif

            if ((out_file) && (nb_alignment <= FILE_NB_ALIGNMENT)) {
              switch (result_format) {

              case 'a' : {
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

              case 's' : {
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

            if (alignment_path) {
              alignment_index++;
            }
          }
        }
      }
    }

#   ifdef DEBUG
    cout << "\ndistance locale : "
         << mean_local_substitution_distance / nb_local_substitution_distance << endl;
#   endif

    if ((ref_identifier == I_DEFAULT) || (test_identifier == I_DEFAULT)) {

      // ecriture des distances, des longueurs des alignements et des nombres d'elisions,
      // d'insertions, de matchs, de substitutions et de transpositions

#     ifdef MESSAGE
      if ((os) && (nb_alignment <= DISPLAY_NB_ALIGNMENT)) {
        *os << "\n";
        dist_matrix->ascii_write(*os);
      }
#     endif

      if (out_file) {
        *out_file << "\n";

        switch (result_format) {
        case 'a' :
          dist_matrix->ascii_write(*out_file);
          break;
        case 's' :
          dist_matrix->spreadsheet_write(*out_file);
          break;
        }
      }
    }

    if (alignment_path) {

      // regroupement elision/insertion et supression des elisions/insertions correspondant
      // aux debuts/fins d'alignement libres

      symbol = new int[(transposition_flag ? TRANSPOSITION : SUBSTITUTION) + 1];
      symbol[DELETION] = 0;
      symbol[INSERTION] = 0;
      symbol[MATCH] = 1;
      symbol[SUBSTITUTION] = 2;
      if (transposition_flag) {
        symbol[TRANSPOSITION] = 3;
      }

#     ifdef MESSAGE
      if (os) {
        *os << "\n" << SEQ_label[SEQL_ALIGNMENT_CODING] << "\n" << STAT_label[STATL_INDEL] << ": " << 0
            << "\n" << STAT_label[STATL_MATCH] << ": " << 1 << "\n" << STAT_label[STATL_SUBSTITUTION] << ": " << 2;
        if (transposition_flag) {
          *os << "\n" << STAT_label[STATL_TRANSPOSITION] << ": " << 3;
        }
        *os << endl;
      }
#     endif

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
          *palignment++ = symbol[*calignment++];
        }
        alignment->length[i] = j - offset;
      }

      delete [] symbol;

      alignment->min_value_computation(0);
      alignment->max_value_computation(0);
      alignment->build_marginal_frequency_distribution(0);

      alignment->max_length_computation();
      alignment->cumul_length_computation();
      alignment->build_length_frequency_distribution();

      // ecriture des sequences d'alignement

      switch (alignment_format) {
      case 'a' :
        status = alignment->ascii_data_write(error , alignment_path);
        break;
      case 'b' :
//        status = alignment->binary_write(error , alignment_path);
        break;
      }

#     ifdef MESSAGE
      if ((!status) && (os)) {
        *os << error;
      }
#     endif

    }

    if (out_file) {
      out_file->close();
      delete out_file;
    }

    delete vector_dist;

    for (i = 0;i < nb_variable;i++) {
      delete [] rank[i];
      delete [] max_symbol_distance[i];
    }
    delete [] rank;
    delete [] max_symbol_distance;

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


/*--------------------------------------------------------------*
 *
 *  Calcul de la distance locale entre deux vecteurs.
 *
 *  arguments : indices des 2 sequences, rangs dans les 2 sequences,
 *              distance de substitution.
 *
 *--------------------------------------------------------------*/

double Sequences::substitution_distance_computation(int ref_index , int test_index , int ref_position ,
                                                    int test_position , double substitution_distance) const

{
  register int i;
  double distance = 0.;


  for (i = 0;i < nb_variable;i++) {
    if (int_sequence[ref_index][i][ref_position] != int_sequence[test_index][i][test_position]) {
      distance = substitution_distance;
      break;
    }
  }

  return distance;
}


/*--------------------------------------------------------------*
 *
 *  Alignement de sequences.
 *
 *  arguments : reference sur un objet StatError, stream, identificateurs de 2 sequences,
 *              flags debut/fin libres, path resultat, format du fichier resultat ('a' : ASCII, 's' : Spreadsheet),
 *              path alignement, format du fichier alignement ('a' : ASCII, 'b' : Binary),
 *
 *--------------------------------------------------------------*/

DistanceMatrix* Sequences::alignment(StatError &error , ostream *os , int ref_identifier ,
                                     int test_identifier , bool begin_free , bool end_free ,
                                     const char *result_path , char result_format ,
                                     const char *alignment_path , char alignment_format) const

{
  bool status = true , half_matrix;
  register int i , j , k , m;
  int nb_alignment , ilength , alignment_index , var , width , ref_position ,
      pref_position , test_position , ptest_position , gap_length , max_gap_length ,
      nb_deletion , nb_insertion , nb_match , nb_begin_end , offset , *palignment ,
      *calignment , *symbol , **path_length , ***back_pointers;
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

  if (((index_parameter_type == TIME) && (index_interval->variance > 0.)) ||
      (index_parameter_type == POSITION)) {
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
    if ((ref_identifier == I_DEFAULT) && (test_identifier == I_DEFAULT) && ((!result_path) ||
         (nb_alignment > FILE_NB_ALIGNMENT)) && (!alignment_path)) {
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

    if (result_path) {
      out_file = new ofstream(result_path);

      if (!out_file) {
        error.update(STAT_error[STATR_FILE_NAME]);

#       ifdef MESSAGE
        if (os) {
          *os << error;
        }
#       endif

      }
    }

    // creation des structures de donnees resultat

    dist_matrix = new DistanceMatrix(nb_sequence , ref_identifier , test_identifier ,
                                     SEQ_label[SEQL_SEQUENCE] , identifier , false);

    if (alignment_path) {
      alignment = new Sequences(nb_alignment , 1);
    }
    else {
      ilength = max_length + max_length;
      alignment = new Sequences(1 , NULL , &ilength , 1 , false);
    }

    // creation des structures de donnees de l'algorithme

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

    // alignement des sequences

    alignment_index = 0;

    for (i = 0;i < nb_sequence;i++) {
      if ((ref_identifier == I_DEFAULT) || (ref_identifier == identifier[i])) {
        for (j = (half_matrix ? i + 1 : 0);j < nb_sequence;j++) {
          if (((test_identifier == I_DEFAULT) || (test_identifier == identifier[j])) && (j != i)) {

            // initialisation des distances cumulees et de la longueur des chemins correspondants

            cumul_distance[0][0] = 0.;
            path_length[0][0] = 0;

            for (k = 1;k <= length[i];k++) {

              // elision

              switch (begin_free) {
              case false :
                cumul_distance[k][0] = cumul_distance[k - 1][0] + INDEL_DISTANCE;
                break;
              case true :
                cumul_distance[k][0] = cumul_distance[k - 1][0];
                break;
              }

              path_length[k][0] = k;
              back_pointers[k][0][0] = k - 1;
              back_pointers[k][0][1] = 0;
            }

            for (k = 1;k <= length[j];k++) {

              // insertion

              switch (begin_free) {
              case false :
                cumul_distance[0][k] = cumul_distance[0][k - 1] + INDEL_DISTANCE;
                break;
              case true :
                cumul_distance[0][k] = cumul_distance[0][k - 1];
                break;
              }

              path_length[0][k] = k;
              back_pointers[0][k][0] = 0;
              back_pointers[0][k][1] = k - 1;
            }

            // calcul des distances cumulees et de la longueur des chemins correspondants

            for (k = 1;k <= length[i];k++) {
              for (m = 1;m <= length[j];m++) {

                // match/substitution

                cumul_distance[k][m] = cumul_distance[k - 1][m - 1] +
                                       substitution_distance_computation(i , j , k - 1 , m - 1 , substitution_distance);
                path_length[k][m] = path_length[k - 1][m - 1] + 1;
                back_pointers[k][m][0] = k - 1;
                back_pointers[k][m][1] = m - 1;

                // elision

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
            if (alignment_path) {
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

            // recherche du nombre d'elisions/insertions successives maximum

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

            // mise a jour des nombre d'elisions, d'insertions et de matchs

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

            // ecriture de l'alignement

#           ifdef MESSAGE
            if ((os) && (nb_alignment <= DISPLAY_NB_ALIGNMENT)) {
              alignment_ascii_print(*os , width , i , j , *alignment , alignment_index);

              if (length[i] + length[j] - nb_begin_end > 0) {
                *os << STAT_label[STATL_DISTANCE] << " (" << SEQ_label[SEQL_ALIGNMENT_LENGTH]
                    << "): " << cumul_distance[length[i]][length[j]] / (length[i] + length[j] - nb_begin_end)
                    << " (" << path_length[length[i]][length[j]] - nb_begin_end << ") = "
                    << deletion_distance / (length[i] + length[j] - nb_begin_end) << " (" << nb_deletion << " d) + "
                    << insertion_distance / (length[i] + length[j] - nb_begin_end) << " (" << nb_insertion << " i) + 0 ("
                    << nb_match << " m)" << endl;

                *os << SEQ_label[SEQL_MAX_GAP_LENGTH] << ": " << max_gap_length << endl;
              }
            }
#           endif

            if ((out_file) && (nb_alignment <= FILE_NB_ALIGNMENT)) {
              switch (result_format) {

              case 'a' : {
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

              case 's' : {
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

            if (alignment_path) {
              alignment_index++;
            }
          }
        }
      }
    }

    if ((ref_identifier == I_DEFAULT) || (test_identifier == I_DEFAULT)) {

      // ecriture des distances, des longueurs des alignements et
      // des nombres d'elisions, d'insertions et de matchs

#     ifdef MESSAGE
      if ((os) && (nb_alignment <= DISPLAY_NB_ALIGNMENT)) {
        *os << "\n";
        dist_matrix->ascii_write(*os);
      }
#     endif

      if (out_file) {
        *out_file << "\n";

        switch (result_format) {
        case 'a' :
          dist_matrix->ascii_write(*out_file);
          break;
        case 's' :
          dist_matrix->spreadsheet_write(*out_file);
          break;
        }
      }
    }

    if (alignment_path) {

      // regroupement elision/insertion et supression des elisions/insertions correspondant
      // aux debuts/fins d'alignement libres

      symbol = new int[MATCH + 1];
      symbol[DELETION] = 0;
      symbol[INSERTION] = 0;
      symbol[MATCH] = 1;

#     ifdef MESSAGE
      if (os) {
        *os << "\n" << SEQ_label[SEQL_ALIGNMENT_CODING] << "\n" << STAT_label[STATL_INDEL] << ": " << 0
            << "\n" << STAT_label[STATL_MATCH] << ": " << 1 << endl;
      }
#     endif

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
          *palignment++ = symbol[*calignment++];
        }
        alignment->length[i] = j - offset;
      }

      delete [] symbol;

      alignment->min_value_computation(0);
      alignment->max_value_computation(0);
      alignment->build_marginal_frequency_distribution(0);

      alignment->max_length_computation();
      alignment->cumul_length_computation();
      alignment->build_length_frequency_distribution();

      // ecriture des sequences d'alignement

      switch (alignment_format) {
      case 'a' :
        status = alignment->ascii_data_write(error , alignment_path);
        break;
      case 'b' :
//        status = alignment->binary_write(error , alignment_path);
        break;
      }

#     ifdef MESSAGE
      if ((!status) && (os)) {
        *os << error;
      }
#     endif

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


/*--------------------------------------------------------------*
 *
 *  Ecriture de l'alignement multiple de sequences.
 *
 *  argument : stream.
 *
 *--------------------------------------------------------------*/

ostream& Sequences::multiple_alignment_ascii_print(ostream &os) const

{
  register int i , j , k , m;
  int var , width , rank;
  long old_adjust;


  old_adjust = os.setf(ios::right , ios::adjustfield);

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

      // ecriture de la sequence consensus

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

  os.setf((FMTFLAGS)old_adjust , ios::adjustfield);
  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture de l'alignement multiple de sequences dans un fichier.
 *
 *  arguments : reference sur un objet StatError, path.
 *
 *--------------------------------------------------------------*/

bool Sequences::multiple_alignment_ascii_print(StatError &error , const char *path) const

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
    multiple_alignment_ascii_print(out_file);
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Alignement multiple de sequences.
 *
 *  arguments : references sur les sequences de test et sur un objet VectorDistance,
 *              rangs, distances maximum aux symboles, flags debut/fin libres,
 *              couts d'elision/insertion adaptatifs ou fixes,
 *              facteur pour deduire les couts d'elision/insertion.
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::multiple_alignment(const Sequences &test_seq , const VectorDistance &vector_dist ,
                                         double **rank , double **max_symbol_distance , bool begin_free ,
                                         bool end_free , int indel_cost , double indel_factor) const

{
  register int i , j , k , m;
  int ref_position , pref_position , test_position , ptest_position , *alignment , *palignment ,
      *ilength , **path_length , ***back_pointers;
  double buff , sum , **ref_local_indel_distance , **test_local_indel_distance , **cumul_distance;
  Sequences *seq;


  // creation des structures de donnees de l'algorithme - calcul des couts d'elision/insertion

  if (indel_cost == FIXED) {
    buff = indel_distance_computation(vector_dist , rank , max_symbol_distance) * indel_factor;
  }

  ref_local_indel_distance = new double*[nb_sequence];
  for (i = 0;i < nb_sequence;i++) {
    ref_local_indel_distance[i] = new double[max_length + 1];

    switch (indel_cost) {

    case ADAPTATIVE : {
      for (j = 1;j <= max_length;j++) {
        if (int_sequence[i][nb_variable - 1][j - 1] == DATA) {
          ref_local_indel_distance[i][j] = indel_distance_computation(vector_dist , i , j - 1 , rank , max_symbol_distance) * indel_factor;
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

    // calcul des couts d'elision/insertion

    switch (indel_cost) {

    case ADAPTATIVE : {
      for (j = 1;j <= test_seq.max_length;j++) {
        if (test_seq.int_sequence[i][nb_variable - 1][j - 1] == DATA) {
          test_local_indel_distance[i][j] = test_seq.indel_distance_computation(vector_dist , i , j - 1 , rank , max_symbol_distance) * indel_factor;
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

  // initialisation des distances cumulees et de la longueur des chemins correspondants

  cumul_distance[0][0] = 0.;
  path_length[0][0] = 0;

  for (i = 1;i <= max_length;i++) {

    // elision

    switch (begin_free) {

    case false : {
      sum = 0.;
      for (j = 0;j < nb_sequence;j++) {
        if (int_sequence[j][nb_variable - 1][i - 1] == DATA) {
          sum += ref_local_indel_distance[j][i];
        }
      }
      cumul_distance[i][0] = cumul_distance[i - 1][0] + sum / nb_sequence;
      break;
    }

    case true : {
      cumul_distance[i][0] = cumul_distance[i - 1][0];
      break;
    }
    }

    path_length[i][0] = i;
    back_pointers[i][0][0] = i - 1;
    back_pointers[i][0][1] = 0;
  }

  for (i = 1;i <= test_seq.max_length;i++) {

    // insertion

    switch (begin_free) {

    case false : {
      sum = 0.;
      for (j = 0;j < test_seq.nb_sequence;j++) {
        if (test_seq.int_sequence[j][nb_variable - 1][i - 1] == DATA) {
          sum += test_local_indel_distance[j][i];
        }
      }
      cumul_distance[0][i] = cumul_distance[0][i - 1] + sum / test_seq.nb_sequence;
      break;
    }

    case true : {
      cumul_distance[0][i] = cumul_distance[0][i - 1];
      break;
    }
    }

    path_length[0][i] = i;
    back_pointers[0][i][0] = 0;
    back_pointers[0][i][1] = i - 1;
  }

  // calcul des distances cumulees et de la longueur des chemins correspondants

  for (i = 1;i <= max_length;i++) {
    for (j = 1;j <= test_seq.max_length;j++) {

      sum = 0.;
      for (k = 0;k < nb_sequence;k++) {
        for (m = 0;m < test_seq.nb_sequence;m++) {

          // calcul distance locale

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

      // elision

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

  // construction du groupe de sequences

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


/*--------------------------------------------------------------*
 *
 *  Alignement multiple de sequences.
 *
 *  arguments : reference sur un objet StatError, stream, reference sur un objet VectorDistance,
 *              flags debut/fin libres, couts d'elision/insertion adaptatifs ou fixes,
 *              facteur pour deduire les couts d'elision/insertion, type d'algorihme pour
 *              la construction du dendrogramme, path.
 *
 *--------------------------------------------------------------*/

Sequences* Sequences::multiple_alignment(StatError &error , ostream &os ,
                                         const VectorDistance &ivector_dist ,
                                         bool begin_free , bool end_free , int indel_cost ,
                                         double indel_factor , int algorithm , const char *path) const

{
  bool status = true;
  register int i , j , k;
  int *itype , *psequence , *csequence , *variable;
  double **rank , **max_symbol_distance;
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

    // alignement des sequences 2 a 2

    dist_matrix = alignment(error , NULL , ivector_dist , I_DEFAULT , I_DEFAULT ,
                            begin_free , end_free , indel_cost , indel_factor);

    if (dist_matrix) {

      // construction d'un dendrogramme a partir de la matrice des distances entre sequences

      if (algorithm != DIVISIVE) {
        dendrogram = dist_matrix->agglomerative_hierarchical_clustering(algorithm);
      }
      else {
        dendrogram = dist_matrix->divisive_hierarchical_clustering();
      }

#     ifdef MESSAGE
      cout << *dendrogram << "\n";
#     endif

      vector_dist = new VectorDistance(ivector_dist);

      // calcul des distance maximum de substitution pour les variables symboliques et
      // des rangs pour les variables ordinales

      rank = new double*[nb_variable];
      max_symbol_distance = new double*[nb_variable];

      for (i = 0;i < nb_variable;i++) {
        if ((vector_dist->get_variable_type(i) == SYMBOLIC) && (vector_dist->get_symbol_distance(i))) {
          max_symbol_distance[i] = vector_dist->max_symbol_distance_computation(i);
        }
        else {
          max_symbol_distance[i] = 0;
        }

        if (vector_dist->get_variable_type(i) == ORDINAL) {
          rank[i] = marginal_distribution[i]->rank_computation();
        }
        else {
          rank[i] = NULL;
        }

        // calcul des dispersions pour la standardisation

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

      // construction des groupes initiaux

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

      // alignement multiple des sequences

      for (i = nb_sequence;i < 2 * nb_sequence - 1;i++) {
        clustered_seq[i] = clustered_seq[dendrogram->get_child(i , 0)]->multiple_alignment(*(clustered_seq[dendrogram->get_child(i , 1)]) ,
                                                                                           *vector_dist , rank ,
                                                                                           max_symbol_distance , begin_free ,
                                                                                           end_free , indel_cost , indel_factor);

#       ifdef DEBUG
        if (i < 2 * nb_sequence - 2) {
          clustered_seq[i]->multiple_alignment_ascii_print(os);
        }
#       endif

      }

      // ecriture de l'alignement multiple

#     ifdef MESSAGE
      clustered_seq[2 * nb_sequence - 2]->multiple_alignment_ascii_print(os);
#     endif

      if (path) {
        status = clustered_seq[2 * nb_sequence - 2]->multiple_alignment_ascii_print(error , path);

#       ifdef MESSAGE
        if (!status) {
          os << error;
        }
#       endif

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
        delete [] max_symbol_distance[i];
      }
      delete [] rank;
      delete [] max_symbol_distance;
    }

    delete dist_matrix;
  }

  return seq;
}


};  // namespace sequence_analysis
