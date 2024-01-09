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
 *       $Id: sequences3.cpp 11060 2011-09-02 16:28:11Z guedon $
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

#include "tool/rw_tokenizer.h"
#include "tool/rw_cstring.h"
#include "tool/rw_locale.h"
#include "tool/config.h"

#include "stat_tool/stat_tools.h"
#include "stat_tool/curves.h"
#include "stat_tool/distribution.h"
#include "stat_tool/markovian.h"
#include "stat_tool/vectors.h"
#include "stat_tool/distance_matrix.h"
#include "stat_tool/stat_label.h"

#include "sequences.h"
#include "sequence_label.h"

using namespace std;
using namespace stat_tool;


namespace sequence_analysis {



/*--------------------------------------------------------------*
 *
 *  Construction d'un objet Sequences a partir d'un fichier.
 *
 *  arguments : reference sur un objet StatError, path, flag format.
 *
 *--------------------------------------------------------------*/

Sequences* sequences_ascii_read(StatError &error , const char *path , bool old_format)

{
  RWLocaleSnapshot locale("en");
  RWCString buffer , token;
  size_t position;
  bool status , lstatus;
  register int i , j , k , m;
  int line , read_line , offset , initial_nb_line , max_length , nb_variable = 0 ,
      index_parameter_type = IMPLICIT_TYPE , vector_size , nb_sequence , index , line_continue ,
      *type , *length;
  long int_value;
  double real_value;
  Sequences *seq;
  ifstream in_file(path);


  seq = NULL;
  error.init();

  if (!in_file) {
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {
    status = true;
    line = 0;
    type = NULL;
    length = NULL;

    // 1ere passe : analyse de la ligne definissant le parametre d'index et
    // de la ligne definissant le nombre de variables

    read_line = 0;

    while (buffer.readLine(in_file , false)) {
      line++;

      position = buffer.first('#');
      if (position != RW_NPOS) {
        buffer.remove(position);
      }
      i = 0;

      RWCTokenizer next(buffer);

      while (!((token = next()).isNull())) {
        switch (i) {

        case 0 : {

          // test mot cle INDEX_PARAMETER

          if ((!old_format) && (read_line == 0) && (token == SEQ_word[SEQW_INDEX_PARAMETER])) {
            index_parameter_type--;
          }

          // test nombre de variables

          else {
            lstatus = locale.stringToNum(token , &int_value);
            if (lstatus) {
              if ((int_value < 1) || (int_value > SEQUENCE_NB_VARIABLE)) {
                lstatus = false;
              }
              else {
                nb_variable = int_value;
              }
            }

            if (!lstatus) {
              status = false;
              error.update(STAT_parsing[STATP_NB_VARIABLE] , line , i + 1);
            }
          }
          break;
        }

        case 1 : {

          // test separateur

          if ((!old_format) && (read_line == 0) && (index_parameter_type != IMPLICIT_TYPE)) {
            if (token != ":") {
              status = false;
              error.update(STAT_parsing[STATP_SEPARATOR] , line , i + 1);
            }
          }

          // test mot cle VARIABLE(S)

          else if (token != STAT_word[nb_variable == 1 ? STATW_VARIABLE : STATW_VARIABLES]) {
            status = false;
            error.correction_update(STAT_parsing[STATP_KEY_WORD] ,
                                    STAT_word[nb_variable == 1 ? STATW_VARIABLE : STATW_VARIABLES] , line , i + 1);
          }
          break;
        }

        // test mot cle correspondant au type du parametre d'index

        case 2 : {
          if ((!old_format) && (read_line == 0) && (index_parameter_type != IMPLICIT_TYPE)) {
            for (j = TIME;j <= POSITION_INTERVAL;j++) {
              if (token == SEQ_index_parameter_word[j]) {
                index_parameter_type = j;
                break;
              }
            }

            if (j == POSITION_INTERVAL + 1) {
              status = false;
              error.update(STAT_parsing[STATP_KEY_WORD] , line , i + 1);
            }
          }
          break;
        }
        }

        i++;
      }

      if (i > 0) {
        if (((!old_format) && (read_line == 0) && (index_parameter_type != IMPLICIT_TYPE) && (i != 3)) ||
            (((old_format) || (read_line == 1) || (index_parameter_type == IMPLICIT_TYPE)) && (i != 2))) {
          status = false;
          error.update(STAT_parsing[STATP_FORMAT] , line);
        }

        read_line++;
//        if (((!old_format) && (index_parameter_type != IMPLICIT_TYPE) && (read_line == 2)) ||
//            (((old_format) || (index_parameter_type == IMPLICIT_TYPE)) && (read_line == 1)) {
        if ((((old_format) || (index_parameter_type == IMPLICIT_TYPE)) && (read_line == 1)) ||
            (read_line == 2)) {
          break;
        }
      }
    }

    if (read_line < 1) {
      status = false;
      error.update(STAT_parsing[STATP_FORMAT]);
    }

    // analyse des lignes definissant le type de chaque variable

    if (status) {
      type = new int[nb_variable];
      for (i = 0;i < nb_variable;i++) {
        type[i] = I_DEFAULT;
      }

      read_line = 0;
      offset = (old_format ? 1 : 0);

      while (buffer.readLine(in_file , false)) {
        line++;

        position = buffer.first('#');
        if (position != RW_NPOS) {
          buffer.remove(position);
        }
        i = 0;

        RWCTokenizer next(buffer);

        while (!((token = next()).isNull())) {
          switch (i) {

          // test mot cle VARIABLE

          case 0 : {
            if (token != STAT_word[STATW_VARIABLE]) {
              status = false;
              error.correction_update(STAT_parsing[STATP_KEY_WORD] , STAT_word[STATW_VARIABLE] , line , i + 1);
            }
            break;
          }

          // test index de la variable

          case 1 : {
            lstatus = locale.stringToNum(token , &int_value);
            if ((lstatus) && (int_value != read_line + 1)) {
              lstatus = false;
            }

            if (!lstatus) {
              status = false;
              error.correction_update(STAT_parsing[STATP_VARIABLE_INDEX] , read_line + 1 , line , i + 1);
            }
            break;
          }

          // test separateur

          case 2 : {
            if (token != ":") {
              status = false;
              error.update(STAT_parsing[STATP_SEPARATOR] , line , i + 1);
            }
            break;
          }

          // test mot cle correspondant au type de la variable

          case 3 : {
            if ((old_format) && (read_line == 0)) {
              for (j = TIME;j <= POSITION_INTERVAL;j++) {
                if (token == SEQ_index_parameter_word[j]) {
                  index_parameter_type = j;
                  break;
                }
              }

              if (j == POSITION_INTERVAL + 1) {
//                for (j = INT_VALUE;j <= STATE;j++) {
                for (j = INT_VALUE;j <= OLD_INT_VALUE;j++) {
                  if (token == STAT_variable_word[j]) {
//                    if (j == STATE) {
                    if ((j == STATE) || (j == OLD_INT_VALUE)) {
                      j = INT_VALUE;
                    }
                    type[read_line] = j;
                    break;
                  }
                }

//                if (j == STATE + 1) {
                if (j == OLD_INT_VALUE + 1) {
                  status = false;
                  error.update(STAT_parsing[STATP_KEY_WORD] , line , i + 1);
                }
              }
            }

            else {
//              for (j = INT_VALUE;j <= NB_INTERNODE;j++) {
              for (j = INT_VALUE;j <= OLD_INT_VALUE;j++) {
                if (token == STAT_variable_word[j]) {
//                  if ((j == NB_INTERNODE) && ((read_line != offset) || ((read_line == offset) &&
                  if ((j == OLD_INT_VALUE) && ((read_line != offset) || ((read_line == offset) &&
                        (index_parameter_type != POSITION) && (index_parameter_type != POSITION_INTERVAL)))) {
                    status = false;
                    error.update(STAT_parsing[STATP_VARIABLE_TYPE] , line , i + 1);
                  }

                  else {
//                    if (j == STATE) {
                    if ((j == STATE) || (j == OLD_INT_VALUE)) {
                      j = INT_VALUE;
                    }

                    if ((old_format) && (index_parameter_type != IMPLICIT_TYPE)) {
                      type[read_line - 1] = j;
                    }
                    else {
                      type[read_line] = j;
                    }
                  }
                  break;
                }
              }
            }

//            if (j == NB_INTERNODE + 1) {
            if (j == OLD_INT_VALUE + 1) {
              status = false;
              error.update(STAT_parsing[STATP_KEY_WORD] , line , i + 1);
            }
            break;
          }
          }

          i++;
        }

        if (i > 0) {
          if (i != 4) {
            status = false;
            error.update(STAT_parsing[STATP_FORMAT] , line);
          }

          read_line++;
          if (read_line == nb_variable) {
            break;
          }
        }
      }

      if (read_line < nb_variable) {
        status = false;
        error.update(STAT_parsing[STATP_FORMAT]);
      }

      else {
//        if ((((index_parameter_type == TIME) || (index_parameter_type == TIME_INTERVAL)) && (read_line < offset + 1)) ||
//            (((index_parameter_type == POSITION) || (index_parameter_type == POSITION_INTERVAL)) &&
//             ((read_line < offset + 1) || ((read_line > offset + 1) && (type[0] == NB_INTERNODE))))) {
        if (((index_parameter_type == TIME) || (index_parameter_type == TIME_INTERVAL) ||
             (index_parameter_type == POSITION) || (index_parameter_type == POSITION_INTERVAL)) &&
             (read_line < offset + 1)) {
          status = false;
          error.update(STAT_parsing[STATP_VARIABLE_TYPE]);
        }
      }

      initial_nb_line = line;
    }

    if (status) {
      vector_size = nb_variable;

      if (index_parameter_type != IMPLICIT_TYPE) {
        if (old_format) {
          nb_variable--;
        }
        else {
          vector_size++;
        }
      }

      nb_sequence = 0;
      lstatus = true;

      while (buffer.readLine(in_file , true)) {
        position = buffer.first('#');
        if (position != RW_NPOS) {
          buffer.remove(position);
        }

        if (!(buffer.isNull())) {
          if (buffer.last('\\') == RW_NPOS) {
            nb_sequence++;
            lstatus = true;
          }
          else {
            lstatus = false;
          }
        }
      }

      if ((nb_sequence == 0) || (!lstatus)) {
        status = false;
        error.update(STAT_parsing[STATP_FORMAT]);
      }

#     ifdef DEBUG
      cout << "\nnumber of sequences : " << nb_sequence << endl;
#     endif
    }

    // 2eme passe : analyse du format des sequences

    if (status) {
//      in_file.close();
//      in_file.open(path , ios::in);

      in_file.clear();
      in_file.seekg(0 , ios::beg);

      offset = (index_parameter_type == IMPLICIT_TYPE ? 0 : 1);

      length = new int[nb_sequence];
      for (i = 0;i < nb_sequence;i++) {
        if (vector_size == 1) {
          length[i] = 0;
        }
        else {
          length[i] = 1;
        }
      }

      line = 0;

      do {
        buffer.readLine(in_file , false);
        line++;

#       ifdef DEBUG
        cout << line << "  " << buffer << endl;
#       endif

      }
      while (line < initial_nb_line);

      max_length = 0;

      switch (index_parameter_type) {
      case TIME :
        index = -1;
        break;
      case POSITION :
        index = 0;
        break;
      }

      i = 0;

      while (buffer.readLine(in_file , false)) {
        line++;

#       ifdef DEBUG
        cout << line << "  " << buffer << endl;
#       endif

        position = buffer.first('#');
        if (position != RW_NPOS) {
          buffer.remove(position);
        }

        j = 0;
        k = 0;
        line_continue = false;

        RWCTokenizer next(buffer);

        while (!((token = next()).isNull())) {
          if (line_continue) {
            status = false;
            error.update(STAT_parsing[STATP_FORMAT] , line , j);
            break;
          }

          if ((vector_size > 1) && (j % (vector_size + 1) == vector_size)) {
            if (token == "\\") {
              line_continue = true;
              length[i]++;
            }

            else {
              if (token != "|") {
                status = false;
                error.update(STAT_parsing[STATP_SEPARATOR] , line , j + 1);
              }
              else {
                k = 0;
                length[i]++;
              }
            }
          }

          else {
            if ((vector_size == 1) && (token == "\\")) {
              line_continue = true;
            }

            else {
              if (((index_parameter_type != IMPLICIT_TYPE) && (k == 0)) ||
                  (type[k - offset] != REAL_VALUE)) {
                lstatus = locale.stringToNum(token , &int_value);
                if ((lstatus) && (((k == 0) && (((index_parameter_type == TIME) || (index_parameter_type == POSITION) ||
                        (index_parameter_type == POSITION_INTERVAL)) && (int_value < 0)) ||
                      ((index_parameter_type == TIME_INTERVAL) && (int_value <= 0))))) {
//                     ((k == 1) && (type[k - 1] == NB_INTERNODE) && (int_value < 0)))) {
                  lstatus = false;
                }

                else {
                  lstatus = locale.stringToNum(token , &real_value);
                }
              }

              if (!lstatus) {
                status = false;
                error.update(STAT_parsing[STATP_DATA_TYPE] , line , j + 1);
              }

              else if (k == 0) {
                switch (index_parameter_type) {

                case TIME : {
                  if (int_value <= index) {
                    status = false;
                    error.update(SEQ_parsing[SEQP_TIME_INDEX_ORDER] , line , j + 1);
                  }
                  else {
                    index = int_value;
                  }
                  break;
                }

                case POSITION : {
                  if (int_value < index) {
                    status = false;
                    error.update(SEQ_parsing[SEQP_POSITION_ORDER] , line , j + 1);
                  }
                  else {
                    index = int_value;
                  }
                  break;
                }
                }
              }

              if (vector_size == 1) {
                length[i]++;
              }
              else {
                k++;
              }
            }
          }

          j++;
        }

        if (j > 0) {
          if (vector_size > 1) {
            if (((line_continue) || ((index_parameter_type != POSITION) && (index_parameter_type != POSITION_INTERVAL))) &&
                (k != vector_size)) {
              status = false;
              error.update(STAT_parsing[STATP_FORMAT] , line , j);
            }
          }

          if (!line_continue) {
            if ((index_parameter_type == POSITION) || (index_parameter_type == POSITION_INTERVAL)) {
              if (k != 1) {
                status = false;
                error.update(STAT_parsing[STATP_FORMAT] , line , j);
              }

              length[i]--;
              if (length[i] == 0) {
                status = false;
                error.update(STAT_parsing[STATP_FORMAT] , line , j);
              }
            }

            switch (index_parameter_type) {
            case TIME :
              index = -1;
              break;
            case POSITION :
              index = 0;
              break;
            }

            if (length[i] > max_length) {
              max_length = length[i];
            }

            if (i < nb_sequence - 1) {
              i++;
            }
          }
        }
      }

      if (max_length <= 1) {
        status = false;
        error.update(SEQ_parsing[SEQP_MAX_SEQUENCE_LENGTH]);
      }

#     ifdef DEBUG
      for (i = 0;i < nb_sequence;i++) {
        cout << length[i] << " ";
      }
      cout << endl;
#     endif
    }

    // 3eme passe : copie des sequences

    if (status) {
//      in_file.close();
//      in_file.open(path , ios::in);

      in_file.clear();
      in_file.seekg(0 , ios::beg);

      seq = new Sequences(nb_sequence , NULL , length , NULL ,
                          index_parameter_type , nb_variable , type);

      line = 0;

      do {
        buffer.readLine(in_file , false);
        line++;
      }
      while (line < initial_nb_line);

      i = 0;
      j = 0;

      while (buffer.readLine(in_file , false)) {
        position = buffer.first('#');
        if (position != RW_NPOS) {
          buffer.remove(position);
        }

        k = 0;
        m = 0;
        line_continue = false;

        RWCTokenizer next(buffer);

        while (!((token = next()).isNull())) {
          if ((vector_size > 1) && (m % (vector_size + 1) == vector_size)) {
            if (token == "\\") {
              line_continue = true;
            }
            k = 0;
            j++;
          }

          else {
            if ((vector_size == 1) && (token == "\\")) {
              line_continue = true;
            }

            else {
              if (((index_parameter_type != IMPLICIT_TYPE) && (k == 0)) ||
                  (type[k - offset] != REAL_VALUE)) {
                locale.stringToNum(token , &int_value);

                if ((index_parameter_type != IMPLICIT_TYPE) && (k == 0)) {
                  seq->index_parameter[i][j] = int_value;
                }
                else {
                  seq->int_sequence[i][k - offset][j] = int_value;
                }
              }

              else {
                locale.stringToNum(token , seq->real_sequence[i][k - offset] + j);
              }

              if (vector_size == 1) {
                j++;
              }
              else {
                k++;
              }
            }
          }

          m++;
        }

        if ((m > 0) && (!line_continue)) {
          i++;
          j = 0;
        }
      }

      if ((seq->index_parameter_type == TIME_INTERVAL) || (seq->index_parameter_type == POSITION_INTERVAL)) {
        seq->index_parameter_computation();
      }

      if (seq->index_parameter) {
        seq->build_index_parameter_frequency_distribution();
      }
//      if ((seq->index_parameter_type == TIME) || ((seq->index_parameter_type == POSITION) &&
//          (seq->type[0] != NB_INTERNODE))) {
      if ((seq->index_parameter_type == TIME) || (seq->index_parameter_type == POSITION)) {
        seq->index_interval_computation();
      }

      for (i = 0;i < nb_variable;i++) {
        seq->min_value_computation(i);
        seq->max_value_computation(i);

        seq->build_marginal_frequency_distribution(i);
      }
    }

    delete [] type;
    delete [] length;
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture sur une ligne d'un objet Sequences.
 *
 *  argument : stream.
 *
 *--------------------------------------------------------------*/

ostream& Sequences::line_write(ostream &os) const

{
  os << nb_sequence << " " << SEQ_label[nb_sequence == 1 ? SEQL_SEQUENCE : SEQL_SEQUENCES] << "   "
     << nb_variable << " " << STAT_word[nb_variable == 1 ? STATW_VARIABLE : STATW_VARIABLES] << "   "
     << SEQ_label[nb_sequence == 1 ? SEQL_LENGTH : SEQL_CUMUL_LENGTH] << ": " << cumul_length;

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Sequences.
 *
 *  arguments : stream, flag niveau de detail, flag commentaire.
 *
 *--------------------------------------------------------------*/

ostream& Sequences::ascii_write(ostream &os , bool exhaustive , bool comment_flag) const

{
  register int i;
  double mean , variance;


  if (index_parameter) {
    os << SEQ_word[SEQW_INDEX_PARAMETER] << " : "
       << SEQ_index_parameter_word[index_parameter_type];
  }

  if (index_parameter_distribution) {
    os << "   ";
    if (comment_flag) {
      os << "# ";
    }
    os << "(" << SEQ_label[SEQL_MIN_INDEX_PARAMETER] << ": " << index_parameter_distribution->offset << ", "
       << SEQ_label[SEQL_MAX_INDEX_PARAMETER] << ": " << index_parameter_distribution->nb_value - 1 << ")" << endl;

    os << "\n";
    if (comment_flag) {
      os << "# ";
    }
    os << (index_parameter_type == TIME ? SEQ_label[SEQL_TIME] : SEQ_label[SEQL_POSITION])
       << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";
    index_parameter_distribution->ascii_characteristic_print(os , false , comment_flag);

    if (exhaustive) {
      os << "\n";
      if (comment_flag) {
        os << "# ";
      }
      os << "   | " << (index_parameter_type == TIME ? SEQ_label[SEQL_TIME] : SEQ_label[SEQL_POSITION])
         << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << endl;
      index_parameter_distribution->ascii_print(os , comment_flag);
    }
  }

  else {
    os << endl;
  }

  if (index_interval) {
    os << "\n";
    if (comment_flag) {
      os << "# ";
    }
    os << (index_parameter_type == TIME ? SEQ_label[SEQL_TIME_INTERVAL] : SEQ_label[SEQL_POSITION_INTERVAL])
       << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";
    index_interval->ascii_characteristic_print(os , false , comment_flag);

    if (exhaustive) {
      os << "\n";
      if (comment_flag) {
        os << "# ";
      }
      os << "   | " << (index_parameter_type == TIME ? SEQ_label[SEQL_TIME_INTERVAL] : SEQ_label[SEQL_POSITION_INTERVAL])
         << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << endl;
      index_interval->ascii_print(os , comment_flag);
    }
  }

  if (index_parameter) {
    os << "\n";
  }
  os << nb_variable << " " << STAT_word[nb_variable == 1 ? STATW_VARIABLE : STATW_VARIABLES] << endl;

  for (i = 0;i < nb_variable;i++) {
    os << "\n" << STAT_word[STATW_VARIABLE] << " " << i + 1 << " : "
       << STAT_variable_word[type[i]];

    if (type[i] != AUXILIARY) {
      os << "   ";
      if (comment_flag) {
        os << "# ";
      }
      os << "(" << STAT_label[STATL_MIN_VALUE] << ": " << min_value[i] << ", "
         << STAT_label[STATL_MAX_VALUE] << ": " << max_value[i] << ")" << endl;

      if (marginal_distribution[i]) {
        os << "\n";
        if (comment_flag) {
          os << "# ";
        }
        os << STAT_label[STATL_MARGINAL] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";

        marginal_distribution[i]->ascii_characteristic_print(os , exhaustive , comment_flag);

        if ((marginal_distribution[i]->nb_value <= ASCII_NB_VALUE) || (exhaustive)) {
          os << "\n";
          if (comment_flag) {
            os << "# ";
          }
          os << "   | " << STAT_label[STATL_FREQUENCY] << endl;
          marginal_distribution[i]->ascii_print(os , comment_flag);
        }
      }

      else {
        os << "\n";
        if (comment_flag) {
          os << "# ";
        }
        os << STAT_label[STATL_SAMPLE_SIZE] << ": " << cumul_length << endl;

        mean = mean_computation(i);
        variance = variance_computation(i , mean);

        if (comment_flag) {
          os << "# ";
        }
        os << STAT_label[STATL_MEAN] << ": " << mean << "   "
           << STAT_label[STATL_VARIANCE] << ": " << variance << "   "
           << STAT_label[STATL_STANDARD_DEVIATION] << ": " << sqrt(variance) << endl;

        if ((variance > 0.) && (exhaustive)) {
          if (comment_flag) {
            os << "# ";
          }
          os << STAT_label[STATL_SKEWNESS_COEFF] << ": " << skewness_computation(i , mean , variance) << "   "
             << STAT_label[STATL_KURTOSIS_COEFF] << ": " << kurtosis_computation(i , mean , variance) << endl;
        }

        if ((marginal_histogram[i]) && (exhaustive)) {
          os << "\n";
          if (comment_flag) {
            os << "# ";
          }
          os << STAT_label[STATL_MARGINAL] << " " << STAT_label[STATL_HISTOGRAM] << endl;

          os << "\n";
          if (comment_flag) {
            os << "# ";
          }
          os << STAT_label[STATL_VALUE] << " | " << STAT_label[STATL_FREQUENCY] << endl;
          marginal_histogram[i]->ascii_print(os , comment_flag);
        }
      }
    }

    else {

#     ifdef MESSAGE
      mean = mean_computation(i);
      variance = variance_computation(i , mean);

      os << "\n";
      if (comment_flag) {
        os << "# ";
      }
      os << STAT_label[STATL_MEAN] << ": " << mean << "   "
         << STAT_label[STATL_VARIANCE] << ": " << variance << "   "
         << STAT_label[STATL_STANDARD_DEVIATION] << ": " << sqrt(variance) << endl;
#     endif

//      os << endl;
    }
  }

  os << "\n";
  if (comment_flag) {
    os << "# ";
  }
  os << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";
  length_distribution->ascii_characteristic_print(os , false , comment_flag);

  if (exhaustive) {
    os << "\n";
    if (comment_flag) {
      os << "# ";
    }
    os << "   | " << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << endl;
    length_distribution->ascii_print(os , comment_flag);
  }

  os << "\n";
  if (comment_flag) {
    os << "# ";
  }
  os << SEQ_label[SEQL_CUMUL_LENGTH] << ": " << cumul_length << endl;

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Sequences.
 *
 *  arguments : stream, flag niveau de detail.
 *
 *--------------------------------------------------------------*/

ostream& Sequences::ascii_write(ostream &os , bool exhaustive) const

{
  return ascii_write(os , exhaustive , false);
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Sequences dans un fichier.
 *
 *  arguments : reference sur un objet StatError, path,
 *              flag niveau de detail.
 *
 *--------------------------------------------------------------*/

bool Sequences::ascii_write(StatError &error , const char *path ,
                            bool exhaustive) const

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
    ascii_write(out_file , exhaustive , false);
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture des sequences.
 *
 *  arguments : stream, format ('c' : column / 'l' : line / 'a' : array / 'p' : posterior),
 *              flag commentaire, probabilites a posteriori des sequences d'etats
 *              les plus probables, entropies des sequences d'etats, nombres de sequences d'etats
 *              (modeles markoviens caches), nombre de caracteres par ligne.
 *
 *--------------------------------------------------------------*/

ostream& Sequences::ascii_print(ostream &os , char format , bool comment_flag ,
                                double *posterior_probability , double *entropy ,
                                double *nb_state_sequence , int line_nb_character) const

{
  register int i , j , k , m;


  switch (format) {

  case 'c' : {
    for (i = 0;i < nb_sequence;i++) {
      os << "\n";

#     ifdef DEBUG
      for (j = 0;j < length[i];j++) {
        if (index_parameter) {
          os << index_parameter[i][j] << " ";
        }
        for (k = 0;k < nb_variable;k++) {
          if ((type[k] != REAL_VALUE) && (type[k] != AUXILIARY)) {
            os << int_sequence[i][k][j] << " ";
          }
          else {
            os << real_sequence[i][k][j] << " ";
          }
        }

        if (j < length[i] - 1) {
          if ((os.rdbuf())->out_waiting() > line_nb_character) {
            os << "\\" << endl;
          }

          else {
            if ((index_parameter) || (nb_variable > 1)) {
              os << "| ";
            }
          }
        }
      }

#     else
      ostringstream sos;

      for (j = 0;j < length[i];j++) {
        if (index_parameter) {
          sos << index_parameter[i][j] << " ";
        }
        for (k = 0;k < nb_variable;k++) {
          if ((type[k] != REAL_VALUE) && (type[k] != AUXILIARY)) {
            sos << int_sequence[i][k][j] << " ";
          }
          else {
            sos << real_sequence[i][k][j] << " ";
          }
        }

        if (j < length[i] - 1) {
          if (sos.str().size() > line_nb_character) {
            os << sos.str() << "\\" << endl;
            sos.str("");
          }

          else {
            if ((index_parameter) || (nb_variable > 1)) {
              sos << "| ";
            }
          }
        }
      }

      os << sos.str();
#     endif

      if (index_parameter_type == POSITION) {
        os << "| " << index_parameter[i][length[i]];
      }

      os << "   ";
      if (comment_flag) {
        os << "# ";
      }
      os << "(" << identifier[i] << ")" << endl;

      if ((posterior_probability) && (entropy) && (nb_state_sequence)) {
        if (comment_flag) {
          os << "# ";
        }
        os << SEQ_label[SEQL_POSTERIOR_STATE_SEQUENCE_PROBABILITY]
           << ": " << posterior_probability[i] << "   "
           << SEQ_label[SEQL_STATE_SEQUENCE_ENTROPY] << ": " << entropy[i] << "   "
           << SEQ_label[SEQL_NB_STATE_SEQUENCE] << ": " << nb_state_sequence[i] << endl;
        if (comment_flag) {
          os << "# ";
        }
//        os << SEQ_label[SEQL_POSTERIOR_STATE_SEQUENCE_PROBABILITY_LOG_RATIO]
//           << ": " << log(nb_state_sequence[i]) + log(posterior_probability[i]) << "   "
        os << SEQ_label[SEQL_STATE_SEQUENCE_DIVERGENCE] << ": "
           << log(nb_state_sequence[i]) - entropy[i] << endl;
      }
    }
    break;
  }

  case 'v' : {
    for (i = 0;i < nb_sequence;i++) {
      os << "\n";

      for (j = 0;j < length[i];j++) {
        if (index_parameter) {
          os << index_parameter[i][j] << " ";
        }
        for (k = 0;k < nb_variable;k++) {
          if ((type[k] != REAL_VALUE) && (type[k] != AUXILIARY)) {
            os << int_sequence[i][k][j] << " ";
          }
          else {
            os << real_sequence[i][k][j] << " ";
          }
        }

        if (j < length[i] - 1) {
          os << "\\" << endl;
        }
      }

      if (index_parameter_type == POSITION) {
        os << "| " << index_parameter[i][length[i]];
      }

      os << "   ";
      if (comment_flag) {
        os << "# ";
      }
      os << "(" << identifier[i] << ")" << endl;

      if ((posterior_probability) && (entropy) && (nb_state_sequence)) {
        if (comment_flag) {
          os << "# ";
        }
        os << SEQ_label[SEQL_POSTERIOR_STATE_SEQUENCE_PROBABILITY]
           << ": " << posterior_probability[i] << "   "
           << SEQ_label[SEQL_STATE_SEQUENCE_ENTROPY] << ": " << entropy[i] << "   "
           << SEQ_label[SEQL_NB_STATE_SEQUENCE] << ": " << nb_state_sequence[i] << endl;
        if (comment_flag) {
          os << "# ";
        }
//        os << SEQ_label[SEQL_POSTERIOR_STATE_SEQUENCE_PROBABILITY_LOG_RATIO]
//           << ": " << log(nb_state_sequence[i]) + log(posterior_probability[i]) << "   "
        os << SEQ_label[SEQL_STATE_SEQUENCE_DIVERGENCE] << ": "
           << log(nb_state_sequence[i]) - entropy[i] << endl;
      }
    }
    break;
  }

  case 'l' : {
    int buff , start , width;
    long old_adjust;


    old_adjust = os.setf(ios::right , ios::adjustfield);

    if (index_parameter) {
      width = column_width(index_parameter_distribution->nb_value - 1);
    }
    else {
      width = 0;
    }

    for (i = 0;i < nb_variable;i++) {
      if ((type[i] != REAL_VALUE) && (type[i] != AUXILIARY)) {
        buff = column_width((int)min_value[i] , (int)max_value[i]);
        if (buff > width) {
          width = buff;
        }
      }

      else {
        for (j = 0;j < nb_sequence;j++) {
          buff = column_width(length[j] , real_sequence[j][i]);
          if (buff > width) {
            width = buff;
          }
        }
      }
    }

    for (i = 0;i < nb_sequence;i++) {
      os << "\n";
      start = 0;

      for (j = 0;j < length[i];j++) {
        os << setw(j == start ? width : width + 1);
        if (index_parameter) {
          os << index_parameter[i][j];
        }
        else if (type[0] != REAL_VALUE) {
          os << int_sequence[i][0][j];
        }
        else {
          os << real_sequence[i][0][j];
        }

        if (j < length[i] - 1) {
          if ((j - start) * (width + 1) > 1000) {
//          if ((j - start) * (width + 1) > line_nb_character) {
            os << " \\" << endl;

            for (k = (index_parameter ? 0 : 1);k < nb_variable;k++) {
              if ((type[k] != REAL_VALUE) && (type[k] != AUXILIARY)) {
                os << setw(width) << int_sequence[i][k][start];
                for (m = start + 1;m <= j;m++) {
                  os << setw(width + 1) << int_sequence[i][k][m];
                }
              }

              else {
                os << setw(width) << real_sequence[i][k][start];
                for (m = start + 1;m <= j;m++) {
                  os << setw(width + 1) << real_sequence[i][k][m];
                }
              }

              os << " \\" << endl;
            }
            start = j + 1;
          }
        }

        else {
          if (index_parameter_type == POSITION) {
            os << setw(width + 1) << index_parameter[i][length[i]];
          }

          for (k = (index_parameter ? 0 : 1);k < nb_variable;k++) {
            if ((type[k] != REAL_VALUE) && (type[k] != AUXILIARY)) {
              os << endl;
              os << setw(width) << int_sequence[i][k][start];
              for (m = start + 1;m <= j;m++) {
                os << setw(width + 1) << int_sequence[i][k][m];
              }
            }

            else {
              os << endl;
              os << setw(width) << real_sequence[i][k][start];
              for (m = start + 1;m <= j;m++) {
                os << setw(width + 1) << real_sequence[i][k][m];
              }
            }
          }
        }
      }

      os << "   ";
      if (comment_flag) {
        os << "# ";
      }
      os << "(" << identifier[i] << ")" << endl;

      if ((posterior_probability) && (entropy) && (nb_state_sequence)) {
        if (comment_flag) {
          os << "# ";
        }
        os << SEQ_label[SEQL_POSTERIOR_STATE_SEQUENCE_PROBABILITY]
           << ": " << posterior_probability[i] << "   "
           << SEQ_label[SEQL_STATE_SEQUENCE_ENTROPY] << ": " << entropy[i] << "   "
           << SEQ_label[SEQL_NB_STATE_SEQUENCE] << ": " << nb_state_sequence[i] << endl;
        if (comment_flag) {
          os << "# ";
        }
//        os << SEQ_label[SEQL_POSTERIOR_STATE_SEQUENCE_PROBABILITY_LOG_RATIO]
//           << ": " << log(nb_state_sequence[i]) + log(posterior_probability[i]) << "   "
        os << SEQ_label[SEQL_STATE_SEQUENCE_DIVERGENCE] << ": "
           << log(nb_state_sequence[i]) - entropy[i] << endl;
      }
    }

    os.setf((FMTFLAGS)old_adjust , ios::adjustfield);
    break;
  }

  case 'a' : {
    os << "[";
    for (i = 0;i < nb_sequence;i++) {

#     ifdef DEBUG
      os << "[";
      for (j = 0;j < length[i];j++) {
        if ((!index_parameter) && (nb_variable == 1)) {
          if (type[0] != REAL_VALUE) {
            os << int_sequence[i][0][j];
          }
          else {
            os << real_sequence[i][0][j];
          }
        }

        else {
          os << "[";
          if (index_parameter) {
            os << index_parameter[i][j] << ",";
          }

          for (k = 0;k < nb_variable - 1;k++) {
            if ((type[k] != REAL_VALUE) && (type[k] != AUXILIARY)) {
              os << int_sequence[i][k][j] << ",";
            }
            else {
              os << real_sequence[i][k][j] << ",";
            }
          }

          if ((type[nb_variable - 1] != REAL_VALUE) && (type[nb_variable - 1] != AUXILIARY)) {
            os << int_sequence[i][nb_variable - 1][j] << "]";
          }
          else {
            os << real_sequence[i][nb_variable - 1][j] << "]";
          }
        }

        if (j < length[i] - 1) {
          os << ",";
          if ((os.rdbuf())->out_waiting() > line_nb_character) {
            os << "\\" << endl;
            os << "  ";
          }
        }
      }

#     else
      ostringstream sos;

      sos << "[";
      for (j = 0;j < length[i];j++) {
        if ((!index_parameter) && (nb_variable == 1)) {
          if (type[0] == REAL_VALUE) {
            sos << int_sequence[i][0][j];
          }
          else {
            sos << real_sequence[i][0][j];
          }
        }

        else {
          sos << "[";
          if (index_parameter) {
            sos << index_parameter[i][j] << ",";
          }

          for (k = 0;k < nb_variable - 1;k++) {
            if ((type[k] != REAL_VALUE) && (type[k] != AUXILIARY)) {
              sos << int_sequence[i][k][j] << ",";
            }
            else {
              sos << real_sequence[i][k][j] << ",";
            }
          }

          if ((type[nb_variable - 1] != REAL_VALUE) && (type[nb_variable - 1] != AUXILIARY)) {
            sos << int_sequence[i][nb_variable - 1][j] << "]";
          }
          else {
            sos << real_sequence[i][nb_variable - 1][j] << "]";
          }
        }

        if (j < length[i] - 1) {
          sos << ",";
          if (sos.str().size() > line_nb_character) {
            os << sos.str() << "\\" << endl;
            os << "  ";
            sos.str("");
          }
        }
      }

      os << sos.str();
#     endif

      if (index_parameter_type == POSITION) {
        os << ",[" << index_parameter[i][length[i]];
        for (j = 1;j < nb_variable;j++) {
          os << "," << I_DEFAULT;
        }
        os << "]";
      }

      os << "]";
      if (i < nb_sequence - 1) {
        os << ",\\" << endl;
        os << " ";
      }
    }
    os << "]" << endl;
    break;
  }

  case 'p' : {
    if ((posterior_probability) && (entropy) && (nb_state_sequence)) {
      bool *selected_sequence;
      int index , width[6];
      long old_adjust;
      double max , *divergence;


      divergence = new double[nb_sequence];
      for (i = 0;i < nb_sequence;i++) {
        divergence[i] = log(nb_state_sequence[i]) - entropy[i];
      }

      width[0] = column_width(nb_sequence) + ASCII_SPACE;
      width[1] = column_width(nb_sequence , posterior_probability) + ASCII_SPACE;
      width[2] = column_width(nb_sequence , entropy) + ASCII_SPACE;
      width[3] = column_width(nb_sequence , divergence) + ASCII_SPACE;
      width[4] = column_width(nb_sequence , nb_state_sequence) + ASCII_SPACE;
      width[5] = column_width(max_length) + ASCII_SPACE;

      selected_sequence = new bool[nb_sequence];
      for (i = 0;i < nb_sequence;i++) {
        selected_sequence[i] = false;
      }

      old_adjust = os.setf(ios::left , ios::adjustfield);

      os << "\n6 " << STAT_word[STATW_VARIABLES] << endl;

      os << "\n" << STAT_word[STATW_VARIABLE] << " 1 : " << STAT_variable_word[INT_VALUE] << endl;

      os << STAT_word[STATW_VARIABLE] << " 2 : " << STAT_variable_word[REAL_VALUE] << "   ";
      if (comment_flag) {
        os << "# ";
      }
      os << SEQ_label[SEQL_POSTERIOR_STATE_SEQUENCE_PROBABILITY] << endl;

      os << STAT_word[STATW_VARIABLE] << " 3 : " << STAT_variable_word[REAL_VALUE] << "   ";
      if (comment_flag) {
        os << "# ";
      }
      os << SEQ_label[SEQL_STATE_SEQUENCE_ENTROPY] << endl;

      os << STAT_word[STATW_VARIABLE] << " 4 : " << STAT_variable_word[REAL_VALUE] << "   ";
      if (comment_flag) {
        os << "# ";
      }
      os << SEQ_label[SEQL_STATE_SEQUENCE_DIVERGENCE] << endl;

      os << STAT_word[STATW_VARIABLE] << " 5 : " << STAT_variable_word[REAL_VALUE] << "   ";
      if (comment_flag) {
        os << "# ";
      }
      os << SEQ_label[SEQL_NB_STATE_SEQUENCE] << endl;

      os << STAT_word[STATW_VARIABLE] << " 6 : " << STAT_variable_word[INT_VALUE] << "    ";
      if (comment_flag) {
        os << "# ";
      }
      os << SEQ_label[SEQL_SEQUENCE_LENGTH] << endl;

      os << "\n";
      for (i = 0;i < nb_sequence;i++) {

#       ifdef DEBUG
        for (j = 0;j < length[i];j++) {  // pour UCs Fuji/Braeburn
          if (int_sequence[i][0][j] <= 1) {
            break;
          }
        }

        if (j < length[i]) {
          os << setw(width[1]) << posterior_probability[i]
             << setw(width[2]) << entropy[i]
             << setw(width[3]) << divergence[i]
             << setw(width[4]) << nb_state_sequence[i]
             << setw(width[5]) << length[i];
          if (comment_flag) {
            os << "# ";
          }
          os << "(" << identifier[i] << ")" << endl;
        }
#       endif

        max = 0.;
        for (j = 0;j < nb_sequence;j++) {
/*          if ((!selected_sequence[j]) && (entropy[j] > max)) {
            max = entropy[j]; */
          if ((!selected_sequence[j]) && (posterior_probability[j] > max)) {
            max = posterior_probability[j];
            index = j;
          }
        }

        selected_sequence[index] = true;

        os << setw(width[0]) << i + 1
           << setw(width[1]) << posterior_probability[index]
           << setw(width[2]) << entropy[index]
           << setw(width[3]) << divergence[index]
           << setw(width[4]) << nb_state_sequence[index]
           << setw(width[5]) << length[index];
        if (comment_flag) {
          os << "# ";
        }
        os << "(" << identifier[index] << ")" << endl;
      }

      delete [] divergence;
      delete [] selected_sequence;

      os.setf((FMTFLAGS)old_adjust , ios::adjustfield);
    }
    break;
  }
  }

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Sequences.
 *
 *  arguments : stream, format (ligne/colonne), flag niveau de detail.
 *
 *--------------------------------------------------------------*/

ostream& Sequences::ascii_data_write(ostream &os , char format , bool exhaustive) const

{
  ascii_write(os , exhaustive , false);
  ascii_print(os , format , false);

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Sequences dans un fichier.
 *
 *  arguments : reference sur un objet StatError, path,
 *              format (ligne/colonne), flag niveau de detail.
 *
 *--------------------------------------------------------------*/

bool Sequences::ascii_data_write(StatError &error , const char *path ,
                                 char format , bool exhaustive) const

{
  bool status = false;
  ofstream out_file(path);


  error.init();

  if (!out_file) {
    status = false;
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {
    status = true;
    if (format != 'a') {
      ascii_write(out_file , exhaustive , true);
    }
    ascii_print(out_file , format , true);
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Sequences dans un fichier au format tableur.
 *
 *  arguments : reference sur un objet StatError, path.
 *
 *--------------------------------------------------------------*/

bool Sequences::spreadsheet_write(StatError &error , const char *path) const

{
  bool status;
  register int i;
  double mean , variance;
  ofstream out_file(path);


  error.init();

  if (!out_file) {
    status = false;
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {
    status = true;

    if (index_parameter) {
      out_file << SEQ_word[SEQW_INDEX_PARAMETER] << "\t"
               << SEQ_index_parameter_word[index_parameter_type];
    }

    if (index_parameter_distribution) {
      out_file << "\t\t" << SEQ_label[SEQL_MIN_INDEX_PARAMETER] << "\t" << index_parameter_distribution->offset
               << "\t\t" << SEQ_label[SEQL_MAX_INDEX_PARAMETER] << "\t" << index_parameter_distribution->nb_value - 1 << endl;

      out_file << "\n" << (index_parameter_type == TIME ? SEQ_label[SEQL_TIME] : SEQ_label[SEQL_POSITION])
               << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\t";
      index_parameter_distribution->spreadsheet_characteristic_print(out_file);

      out_file << "\n\t" << (index_parameter_type == TIME ? SEQ_label[SEQL_TIME] : SEQ_label[SEQL_POSITION])
               << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << endl;
      index_parameter_distribution->spreadsheet_print(out_file);
    }

    else {
      out_file << endl;
    }

    if (index_interval) {
      out_file << "\n" << (index_parameter_type == TIME ? SEQ_label[SEQL_TIME_INTERVAL] : SEQ_label[SEQL_POSITION_INTERVAL])
               << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\t";
      index_interval->spreadsheet_characteristic_print(out_file);

      out_file << "\n\t" << (index_parameter_type == TIME ? SEQ_label[SEQL_TIME_INTERVAL] : SEQ_label[SEQL_POSITION_INTERVAL])
               << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << endl;
      index_interval->spreadsheet_print(out_file);
    }

    if (index_parameter) {
      out_file << "\n";
    }
    out_file << nb_variable << "\t" << STAT_word[nb_variable == 1 ? STATW_VARIABLE : STATW_VARIABLES] << endl;

    for (i = 0;i < nb_variable;i++) {
      out_file << "\n" << STAT_word[STATW_VARIABLE] << "\t" << i + 1 << "\t"
               << STAT_variable_word[type[i]];

      if (type[i] != AUXILIARY) {
        out_file << "\t\t" << STAT_label[STATL_MIN_VALUE] << "\t" << min_value[i]
                 << "\t\t" << STAT_label[STATL_MAX_VALUE] << "\t" << max_value[i] << endl;

        if (marginal_distribution[i]) {
          out_file << "\n" << STAT_label[STATL_MARGINAL] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\t";
          marginal_distribution[i]->spreadsheet_characteristic_print(out_file);

          out_file << "\n\t" << STAT_label[STATL_FREQUENCY] << endl;
          marginal_distribution[i]->spreadsheet_print(out_file);
        }

        else {
          out_file << "\n" << STAT_label[STATL_SAMPLE_SIZE] << "\t" << cumul_length << endl;

          mean = mean_computation(i);
          variance = variance_computation(i , mean);

          out_file << STAT_label[STATL_MEAN] << "\t" << mean << "\t\t"
                   << STAT_label[STATL_VARIANCE] << "\t" << variance << "\t\t"
                   << STAT_label[STATL_STANDARD_DEVIATION] << "\t" << sqrt(variance) << endl;

          if (variance > 0.) {
            out_file << STAT_label[STATL_SKEWNESS_COEFF] << "\t" << skewness_computation(i , mean , variance) << "\t\t"
                     << STAT_label[STATL_KURTOSIS_COEFF] << "\t" << kurtosis_computation(i , mean , variance) << endl;
          }

          if (marginal_histogram[i]) {
            out_file << "\n" << STAT_label[STATL_MARGINAL] << " " << STAT_label[STATL_HISTOGRAM] << endl;
            out_file << "\n" << STAT_label[STATL_VALUE] << "\t" << STAT_label[STATL_FREQUENCY] << endl;
            marginal_histogram[i]->spreadsheet_print(out_file);
          }
        }
      }

      else {
        out_file << endl;
      }
    }

    out_file << "\n" << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\t";
    length_distribution->spreadsheet_characteristic_print(out_file);

    out_file << "\n\t" << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << endl;
    length_distribution->spreadsheet_print(out_file);

    out_file << "\n" << SEQ_label[SEQL_CUMUL_LENGTH] << "\t" << cumul_length << endl;
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Sortie Gnuplot d'un objet Sequences.
 *
 *  arguments : reference sur un objet StatError, prefixe des fichiers,
 *              titre des figures.
 *
 *--------------------------------------------------------------*/

bool Sequences::plot_write(StatError &error , const char *prefix ,
                           const char *title) const

{
  bool status;
  register int i , j;
  int nb_histo;
  const FrequencyDistribution *phisto[2];
  ostringstream *data_file_name;


  error.init();

  // ecriture des fichier de donnees

  data_file_name = new ostringstream[nb_variable + 1];
  data_file_name[0] << prefix << 0 << ".dat";

  nb_histo = 0;

  if (index_parameter_distribution) {
    phisto[nb_histo++] = index_parameter_distribution;
  }
  if (index_interval) {
    phisto[nb_histo++] = index_interval;
  }

  status = length_distribution->plot_print((data_file_name[0].str()).c_str() , nb_histo , phisto);

  if (!status) {
    error.update(STAT_error[STATR_FILE_PREFIX]);
  }

  else {
    for (i = 0;i < nb_variable;i++) {
      if (marginal_distribution[i]) {
        data_file_name[i + 1] << prefix << i + 1 << ".dat";
        marginal_distribution[i]->plot_print((data_file_name[i + 1].str()).c_str());
      }
      else if (marginal_histogram[i]) {
        data_file_name[i + 1] << prefix << i + 1 << ".dat";
        marginal_histogram[i]->plot_print((data_file_name[i + 1].str()).c_str());
      }
    }

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

      out_file << "set border 15 lw 0\n" << "set tics out\n" << "set xtics nomirror\n"
               << "set title";
      if (title) {
        out_file << " \"" << title << "\"";
      }
      out_file << "\n\n";

      j = 2;

      if (index_parameter_distribution) {
        if (index_parameter_distribution->nb_value - 1 < TIC_THRESHOLD) {
          out_file << "set xtics 0,1" << endl;
        }
        if ((int)(index_parameter_distribution->max * YSCALE) + 1 < TIC_THRESHOLD) {
          out_file << "set ytics 0,1" << endl;
        }

        out_file << "plot [" << index_parameter_distribution->offset << ":"
                 << index_parameter_distribution->nb_value - 1 << "] [0:"
                 << (int)(index_parameter_distribution->max * YSCALE) + 1 << "] \""
                 << label((data_file_name[0].str()).c_str()) << "\" using " << j++ << " title \""
                 << (index_parameter_type == TIME ? SEQ_label[SEQL_TIME] : SEQ_label[SEQL_POSITION])
                 << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\" with impulses" << endl;

        if (index_parameter_distribution->nb_value - 1 < TIC_THRESHOLD) {
          out_file << "set xtics autofreq" << endl;
        }
        if ((int)(index_parameter_distribution->max * YSCALE) + 1 < TIC_THRESHOLD) {
          out_file << "set ytics autofreq" << endl;
        }

        if (i == 0) {
          out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
        }
        out_file << endl;
      }

      if (index_interval) {
        if (index_interval->nb_value - 1 < TIC_THRESHOLD) {
          out_file << "set xtics 0,1" << endl;
        }
        if ((int)(index_interval->max * YSCALE) + 1 < TIC_THRESHOLD) {
          out_file << "set ytics 0,1" << endl;
        }

        out_file << "plot [0:" << index_interval->nb_value - 1 << "] [0:"
                 << (int)(index_interval->max * YSCALE) + 1 << "] \""
                 << label((data_file_name[0].str()).c_str()) << "\" using " << j++ << " title \""
                 << (index_parameter_type == TIME ? SEQ_label[SEQL_TIME_INTERVAL] : SEQ_label[SEQL_POSITION_INTERVAL])
                 << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\" with impulses" << endl;

        if (index_interval->nb_value - 1 < TIC_THRESHOLD) {
          out_file << "set xtics autofreq" << endl;
        }
        if ((int)(index_interval->max * YSCALE) + 1 < TIC_THRESHOLD) {
          out_file << "set ytics autofreq" << endl;
        }

        if (i == 0) {
          out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
        }
        out_file << endl;
      }

      for (j = 0;j < nb_variable;j++) {
        if (marginal_distribution[j]) {
          if (marginal_distribution[j]->nb_value - 1 < TIC_THRESHOLD) {
            out_file << "set xtics 0,1" << endl;
          }
          if ((int)(marginal_distribution[j]->max * YSCALE) + 1 < TIC_THRESHOLD) {
            out_file << "set ytics 0,1" << endl;
          }

          out_file << "plot [0:" << MAX(marginal_distribution[j]->nb_value - 1 , 1) << "] [0:"
                   << (int)(marginal_distribution[j]->max * YSCALE) + 1 << "] \""
                   << label((data_file_name[j + 1].str()).c_str()) << "\" using 1 title \"";
          if (nb_variable > 1) {
            out_file << STAT_label[STATL_VARIABLE] << " " << j + 1 << " - ";
          }
          out_file << STAT_label[STATL_MARGINAL] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
                   << "\" with impulses" << endl;

          if (marginal_distribution[j]->nb_value - 1 < TIC_THRESHOLD) {
            out_file << "set xtics autofreq" << endl;
          }
          if ((int)(marginal_distribution[j]->max * YSCALE) + 1 < TIC_THRESHOLD) {
            out_file << "set ytics autofreq" << endl;
          }

          if (i == 0) {
            out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
          }
          out_file << endl;
        }

        else if (marginal_histogram[j]) {
          if ((int)(marginal_histogram[j]->max * YSCALE) + 1 < TIC_THRESHOLD) {
            out_file << "set ytics 0,1" << endl;
          }

          out_file << "plot [" << marginal_histogram[j]->min_value - marginal_histogram[j]->step << ":"
                   << marginal_histogram[j]->max_value + marginal_histogram[j]->step << "] [0:"
                   << (int)(marginal_histogram[j]->max * YSCALE) + 1 << "] \""
                   << label((data_file_name[j + 1].str()).c_str()) << "\" using 1:2 title \""
                   << STAT_label[STATL_VARIABLE] << " " << j + 1 << " - "
                   << STAT_label[STATL_MARGINAL] << " " << STAT_label[STATL_HISTOGRAM]
                   << "\" with histeps" << endl;

          if ((int)(marginal_histogram[j]->max * YSCALE) + 1 < TIC_THRESHOLD) {
            out_file << "set ytics autofreq" << endl;
          }

          if (i == 0) {
            out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
          }
          out_file << endl;
        }
      }

      if (length_distribution->nb_value - 1 < TIC_THRESHOLD) {
        out_file << "set xtics 0,1" << endl;
      }
      if ((int)(length_distribution->max * YSCALE) + 1 < TIC_THRESHOLD) {
        out_file << "set ytics 0,1" << endl;
      }

      out_file << "plot [0:" << length_distribution->nb_value - 1 << "] [0:"
               << (int)(length_distribution->max * YSCALE) + 1 << "] \""
               << label((data_file_name[0].str()).c_str()) << "\" using 1 title \""
               << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
               << "\" with impulses" << endl;

      if (length_distribution->nb_value - 1 < TIC_THRESHOLD) {
        out_file << "set xtics autofreq" << endl;
      }
      if ((int)(length_distribution->max * YSCALE) + 1 < TIC_THRESHOLD) {
        out_file << "set ytics autofreq" << endl;
      }

      if (i == 1) {
        out_file << "\nset terminal x11" << endl;
      }

      out_file << "\npause 0 \"" << STAT_label[STATL_END] << "\"" << endl;
    }
  }

  delete [] data_file_name;

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Sortie graphique d'un objet Sequences.
 *
 *--------------------------------------------------------------*/

MultiPlotSet* Sequences::get_plotable() const

{
  register int i , j;
  int nb_plot_set;
  ostringstream legend;
  MultiPlotSet *plot_set;


  nb_plot_set = 1;
  if (index_parameter_distribution) {
    nb_plot_set++;
  }
  if (index_interval) {
    nb_plot_set++;
  }
  for (i = 0;i < nb_variable;i++) {
    if ((marginal_distribution[i]) || (marginal_histogram[i])) {
      nb_plot_set++;
    }
  }

  plot_set = new MultiPlotSet(nb_plot_set);

  MultiPlotSet &plot = *plot_set;

  plot.border = "15 lw 0";

  i = 0;

  if (index_parameter_distribution) {

    // vue : loi empirique des parametres d'index

    plot[i].xrange = Range(index_parameter_distribution->offset , index_parameter_distribution->nb_value - 1);
    plot[i].yrange = Range(0 , ceil(index_parameter_distribution->max * YSCALE));

    if (index_parameter_distribution->nb_value - 1 < TIC_THRESHOLD) {
      plot[i].xtics = 1;
    }
    if (ceil(index_parameter_distribution->max * YSCALE) < TIC_THRESHOLD) {
      plot[i].ytics = 1;
    }

    plot[i].resize(1);

    legend.str("");
    legend << (index_parameter_type == TIME ? SEQ_label[SEQL_TIME] : SEQ_label[SEQL_POSITION])
           << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
    plot[i][0].legend = legend.str();

    plot[i][0].style = "impulses";

    index_parameter_distribution->plotable_frequency_write(plot[i][0]);
    i++;
  }

  if (index_interval) {

    // vue : loi empirique des intervalles entre parametres d'index successifs

    plot[i].xrange = Range(0 , index_interval->nb_value - 1);
    plot[i].yrange = Range(0 , ceil(index_interval->max * YSCALE));

    if (index_interval->nb_value - 1 < TIC_THRESHOLD) {
      plot[i].xtics = 1;
    }
    if (ceil(index_interval->max * YSCALE) < TIC_THRESHOLD) {
      plot[i].ytics = 1;
    }

    plot[i].resize(1);

    legend.str("");
    legend << (index_parameter_type == TIME ? SEQ_label[SEQL_TIME_INTERVAL] : SEQ_label[SEQL_POSITION_INTERVAL])
           << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
    plot[i][0].legend = legend.str();

    plot[i][0].style = "impulses";

    index_interval->plotable_frequency_write(plot[i][0]);
    i++;
  }

  for (j = 0;j < nb_variable;j++) {
    if (marginal_distribution[j]) {

      // vue : loi marginale empirique

      plot[i].xrange = Range(0 , MAX(marginal_distribution[j]->nb_value - 1 , 1));
      plot[i].yrange = Range(0 , ceil(marginal_distribution[j]->max * YSCALE));

      if (marginal_distribution[j]->nb_value - 1 < TIC_THRESHOLD) {
        plot[i].xtics = 1;
      }
      if (ceil(marginal_distribution[j]->max * YSCALE) < TIC_THRESHOLD) {
        plot[i].ytics = 1;
      }

      plot[i].resize(1);

      legend.str("");
      if (nb_variable > 1) {
        legend << STAT_label[STATL_VARIABLE] << " " << j + 1 << " - ";
      }
      legend << STAT_label[STATL_MARGINAL] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
      plot[i][0].legend = legend.str();

      plot[i][0].style = "impulses";

      marginal_distribution[j]->plotable_frequency_write(plot[i][0]);
      i++;
    }

    else if (marginal_histogram[j]) {

      // vue : histogramme marginal

      plot[i].xrange = Range(marginal_histogram[j]->min_value - marginal_histogram[j]->step ,
                             marginal_histogram[j]->max_value + marginal_histogram[j]->step);
      plot[i].yrange = Range(0 , ceil(marginal_histogram[j]->max * YSCALE));

      if (ceil(marginal_histogram[j]->max * YSCALE) < TIC_THRESHOLD) {
        plot[i].ytics = 1;
      }

      plot[i].resize(1);

      legend.str("");
      legend << STAT_label[STATL_VARIABLE] << " " << j + 1 << " "
             << STAT_label[STATL_MARGINAL] << " " << STAT_label[STATL_HISTOGRAM];
      plot[i][0].legend = legend.str();

      plot[i][0].style = "histeps";

      marginal_histogram[j]->plotable_write(plot[i][0]);
      i++;
    }
  }

  // vue : loi empirique des longueurs des sequences

  plot[i].xrange = Range(0 , length_distribution->nb_value - 1);
  plot[i].yrange = Range(0 , ceil(length_distribution->max * YSCALE));

  if (length_distribution->nb_value - 1 < TIC_THRESHOLD) {
    plot[i].xtics = 1;
  }
  if (ceil(length_distribution->max * YSCALE) < TIC_THRESHOLD) {
    plot[i].ytics = 1;
  }

  plot[i].resize(1);

  legend.str("");
  legend << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
  plot[i][0].legend = legend.str();

  plot[i][0].style = "impulses";

  length_distribution->plotable_frequency_write(plot[i][0]);

  return plot_set;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture des sequences au format Gnuplot.
 *
 *  arguments : path, longueur.
 *
 *--------------------------------------------------------------*/

bool Sequences::plot_print(const char *path , int ilength) const

{
  bool status = false;
  register int i , j , k;
  int length_nb_sequence , *index , *plength;
  ofstream out_file(path);


  if (out_file) {
    status = true;

    index = new int[nb_sequence];

    length_nb_sequence = 0;
    plength = length;
    for (i = 0;i < nb_sequence;i++) {
      if (*plength++ == ilength) {
        index[length_nb_sequence++] = i;
      }
    }

    for (i = 0;i < ilength;i++) {
      for (j = 0;j < length_nb_sequence;j++) {
        if (index_parameter) {
          out_file << index_parameter[index[j]][i] << " ";
        }
        for (k = 0;k < nb_variable;k++) {
          if ((type[k] != REAL_VALUE) && (type[k] != AUXILIARY)) {
            out_file << int_sequence[index[j]][k][i] << " ";
          }
          else {
            out_file << real_sequence[index[j]][k][i] << " ";
          }
        }
      }
      out_file << endl;
    }

    delete [] index;
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Sortie Gnuplot d'un objet Sequences.
 *
 *  arguments : reference sur un objet StatError, prefixe des fichiers,
 *              titre des figures.
 *
 *--------------------------------------------------------------*/

bool Sequences::plot_data_write(StatError &error , const char *prefix ,
                                const char *title) const

{
  bool status;
  register int i , j , k;
  int min_index_parameter , max_index_parameter , *pfrequency , *length_nb_sequence;
  double min , max;
  ostringstream *data_file_name;


  error.init();

  if (nb_sequence > PLOT_NB_SEQUENCE) {
    status = false;
    error.update(SEQ_error[SEQR_NB_SEQUENCE]);
  }

  else {

    // ecriture des fichiers de donnees

    data_file_name = new ostringstream[length_distribution->nb_value];

    data_file_name[length_distribution->offset] << prefix << length_distribution->offset << ".dat";
    status = plot_print((data_file_name[length_distribution->offset].str()).c_str() , length_distribution->offset);

    if (!status) {
      error.update(STAT_error[STATR_FILE_PREFIX]);
    }

    else {
      pfrequency = length_distribution->frequency + length_distribution->offset + 1;
      for (i = length_distribution->offset + 1;i < length_distribution->nb_value;i++) {
        if (*pfrequency++ > 0) {
          data_file_name[i] << prefix << i << ".dat";
          plot_print((data_file_name[i].str()).c_str() , i);
        }
      }

      length_nb_sequence = new int[length_distribution->nb_value];

      if (index_parameter) {
        min_index_parameter = index_parameter_distribution->offset;
        max_index_parameter = max_index_parameter_computation(true);
      }

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

        for (j = 0;j < nb_variable;j++) {
          for (k = 0;k < length_distribution->nb_value;k++) {
            length_nb_sequence[k] = 0;
          }

          out_file << "set title \"";
          if (title) {
            out_file << title;
            if (nb_variable > 1) {
              out_file << " - ";
            }
          }

          if (nb_variable > 1) {
            out_file << STAT_label[STATL_VARIABLE] << " " << j + 1;
          }
          out_file << "\"\n\n";

          min = min_value[j];
          max = max_value[j];
          if ((j + 1 < nb_variable) && (type[j + 1] == AUXILIARY)) {
            if (min_value[j + 1] < min) {
              min = min_value[j + 1];
            }
            if (max_value[j + 1] > max) {
              max = max_value[j + 1];
            }
          }

          if (index_parameter) {
            if (max_index_parameter - min_index_parameter < TIC_THRESHOLD) {
              out_file << "set xtics 0,1" << endl;
            }
            if (max - min < TIC_THRESHOLD) {
              out_file << "set ytics " << MIN(min , 0) << ",1" << endl;
            }

            out_file << "plot [" << min_index_parameter << ":" << max_index_parameter << "] ["
                     << MIN(min , 0) << ":" << MAX(max , min + 1) << "] ";
            for (k = 0;k < nb_sequence;k++) {
              out_file << "\"" << label((data_file_name[length[k]].str()).c_str()) << "\" using "
                       << length_nb_sequence[length[k]] * (nb_variable + 1) + 1 << " : "
                       << length_nb_sequence[length[k]] * (nb_variable + 1) + j + 2;
              if (nb_sequence <= PLOT_LEGEND_NB_SEQUENCE) {
                out_file << " title \"" << identifier[k] << "\" with linespoints";
              }
              else {
                out_file << " notitle with linespoints";
              }

              if ((j + 1 < nb_variable) && (type[j + 1] == AUXILIARY)) {
                out_file << ",\\" << endl;
                out_file << "\"" << label((data_file_name[length[k]].str()).c_str()) << "\" using "
                         << length_nb_sequence[length[k]] * (nb_variable + 1) + 1 << " : "
                         << length_nb_sequence[length[k]] * (nb_variable + 1) + j + 3;
                if (nb_sequence <= PLOT_LEGEND_NB_SEQUENCE) {
                  out_file << " title \"" << identifier[k] << "\" with lines";
                }
                else {
                  out_file << " notitle with lines";
                }
              }

              if (k < nb_sequence - 1) {
                out_file << ",\\";
              }
              out_file << endl;
              length_nb_sequence[length[k]]++;
            }

            if (max_index_parameter - min_index_parameter < TIC_THRESHOLD) {
              out_file << "set xtics autofreq" << endl;
            }
            if (max_value[j] - min_value[j] < TIC_THRESHOLD) {
              out_file << "set ytics autofreq" << endl;
            }
          }

          else {
            if (max_length - 1 < TIC_THRESHOLD) {
              out_file << "set xtics 0,1" << endl;
            }
            if (max_value[j] - min_value[j] < TIC_THRESHOLD) {
              out_file << "set ytics " << MIN(min_value[j] , 0) << ",1" << endl;
            }

            out_file << "plot [0:" << max_length - 1 << "] [" << MIN(min , 0)
                     << ":" << MAX(max , min + 1) << "] ";
            for (k = 0;k < nb_sequence;k++) {
              out_file << "\"" << label((data_file_name[length[k]].str()).c_str()) << "\" using "
                       << length_nb_sequence[length[k]] * nb_variable + j + 1;
              if (nb_sequence <= PLOT_LEGEND_NB_SEQUENCE) {
                out_file << " title \"" << identifier[k] << "\" with linespoints";
              }
              else {
                out_file << " notitle with linespoints";
              }

              if ((j + 1 < nb_variable) && (type[j + 1] == AUXILIARY)) {
                out_file << ",\\" << endl;
                out_file << "\"" << label((data_file_name[length[k]].str()).c_str()) << "\" using "
                         << length_nb_sequence[length[k]] * nb_variable + j + 2;
                if (nb_sequence <= PLOT_LEGEND_NB_SEQUENCE) {
                  out_file << " title \"" << identifier[k] << "\" with lines";
                }
                else {
                  out_file << " notitle with lines";
                }
              }

              if (k < nb_sequence - 1) {
                out_file << ",\\";
              }
              out_file << endl;
              length_nb_sequence[length[k]]++;
            }

            if (max_length - 1 < TIC_THRESHOLD) {
              out_file << "set xtics autofreq" << endl;
            }
            if (max_value[j] - min_value[j] - 1 < TIC_THRESHOLD) {
              out_file << "set ytics autofreq" << endl;
            }
          }

          if ((j + 1 < nb_variable) && (type[j + 1] == AUXILIARY)) {
            j++;
          }

          if ((i == 0) && (((type[nb_variable - 1] != AUXILIARY) && (j < nb_variable - 1)) ||
               ((type[nb_variable - 1] == AUXILIARY) && (j < nb_variable - 2)))) {
            out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
          }
          out_file << endl;
        }

        if (i == 1) {
          out_file << "\nset terminal x11" << endl;
        }

        out_file << "\npause 0 \"" << STAT_label[STATL_END] << "\"" << endl;
      }

      delete [] length_nb_sequence;
    }

    delete [] data_file_name;
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Sortie graphique d'un objet Sequences.
 *
 *  argument : reference sur un objet StatError.
 *
 *--------------------------------------------------------------*/

MultiPlotSet* Sequences::get_plotable_data(StatError &error) const

{
  register int i , j , k , m , n;
  int nb_plot_set , min_index_parameter , max_index_parameter;
  double min , max;
  ostringstream title , legend;
  MultiPlotSet *plot_set;


  error.init();

  if (nb_sequence > PLOT_NB_SEQUENCE) {
    plot_set = NULL;
    error.update(SEQ_error[SEQR_NB_SEQUENCE]);
  }

  else {
    nb_plot_set = 0;
    for (i = 0;i < nb_variable;i++) {
      if (type[i] != AUXILIARY) {
        nb_plot_set++;
      }
    }

    plot_set = new MultiPlotSet(nb_plot_set);

    MultiPlotSet &plot = *plot_set;

    plot.border = "15 lw 0";

    if (index_parameter) {
      min_index_parameter = index_parameter_distribution->offset;
      max_index_parameter = max_index_parameter_computation(true);
    }

    i = 0;
    for (j = 0;j < nb_variable;j++) {
      if (type[j] != AUXILIARY) {
        if (nb_variable > 1) {
          title.str("");
          title << STAT_label[STATL_VARIABLE] << " " << i + 1;
          plot[i].title = title.str();
        }

        min = min_value[j];
        max = max_value[j];
        if ((j + 1 < nb_variable) && (type[j + 1] == AUXILIARY)) {
          if (min_value[j + 1] < min) {
            min = min_value[j + 1];
          }
          if (max_value[j + 1] > max) {
            max = max_value[j + 1];
          }
        }

        plot[i].yrange = Range(MIN(min , 0) , MAX(max , min + 1));
        if (max - min < TIC_THRESHOLD) {
          plot[i].ytics = 1;
        }

        if ((j + 1 < nb_variable) && (type[j + 1] == AUXILIARY)) {
          plot[i].resize(nb_sequence * 2);
        }
        else {
          plot[i].resize(nb_sequence);
        }

        if (index_parameter) {
          plot[i].xrange = Range(min_index_parameter , max_index_parameter);
          if (max_index_parameter - min_index_parameter < TIC_THRESHOLD) {
            plot[i].xtics = 1;
          }

          k = 0;
          for (m = 0;m < nb_sequence;m++) {
            plot[i][k].style = "linespoints";

            if (nb_sequence <= PLOT_LEGEND_NB_SEQUENCE) {
              legend.str("");
              legend << identifier[m];
              plot[i][k].legend = legend.str();
            }

            if (type[j] != REAL_VALUE) {
              for (n = 0;n < length[m];n++) {
                plot[i][k].add_point(index_parameter[m][n] , int_sequence[m][j][n]);
              }
            }

            else {
              for (n = 0;n < length[m];n++) {
                plot[i][k].add_point(index_parameter[m][n] , real_sequence[m][j][n]);
              }
            }
            k++;

            if ((j + 1 < nb_variable) && (type[j + 1] == AUXILIARY)) {
              plot[i][k].style = "lines";

              if (nb_sequence <= PLOT_LEGEND_NB_SEQUENCE) {
                legend.str("");
                legend << identifier[m];
                plot[i][k].legend = legend.str();
              }

              for (n = 0;n < length[m];n++) {
                plot[i][k].add_point(index_parameter[m][n] , real_sequence[m][j + 1][n]);
              }

              k++;
            }
          }
        }

        else {
          plot[i].xrange = Range(0 , max_length - 1);
          if (max_length - 1 < TIC_THRESHOLD) {
            plot[i].xtics = 1;
          }

          k = 0;
          for (m = 0;m < nb_sequence;m++) {
            plot[i][k].style = "linespoints";

            if (nb_sequence <= PLOT_LEGEND_NB_SEQUENCE) {
              legend.str("");
              legend << identifier[m];
              plot[i][k].legend = legend.str();
            }

            if (type[j] != REAL_VALUE) {
              for (n = 0;n < length[m];n++) {
                plot[i][k].add_point(n , int_sequence[m][j][n]);
              }
            }

            else {
              for (n = 0;n < length[m];n++) {
                plot[i][k].add_point(n , real_sequence[m][j][n]);
              }
            }
            k++;

            if ((j + 1 < nb_variable) && (type[j + 1] == AUXILIARY)) {
              plot[i][k].style = "lines";

              if (nb_sequence <= PLOT_LEGEND_NB_SEQUENCE) {
                legend.str("");
                legend << identifier[m];
                plot[i][k].legend = legend.str();
              }

              for (n = 0;n < length[m];n++) {
                plot[i][k].add_point(n , real_sequence[m][j + 1][n]);
              }

              k++;
            }
          }
        }

        i++;
      }
    }
  }

  return plot_set;
}


};  // namespace sequence_analysis
