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
#include "tool/rw_tokenizer.h"
#include "tool/rw_cstring.h"
#include "tool/rw_locale.h"

#include "stat_tool/stat_tools.h"
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
extern char* label(const char *file_name);



/*--------------------------------------------------------------*
 *
 *  Construction d'un objet Sequences a partir d'un fichier.
 *
 *  arguments : reference sur un objet Format_error, path, flag format.
 *
 *--------------------------------------------------------------*/

Sequences* sequences_ascii_read(Format_error &error , const char *path , bool old_format)

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
              for (j = INT_VALUE;j <= NB_INTERNODE;j++) {
                if (token == STAT_variable_word[j]) {
                  if ((j == NB_INTERNODE) && ((read_line != offset) || ((read_line == offset) &&
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

            if (j == NB_INTERNODE + 1) {
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
        if ((((index_parameter_type == TIME) || (index_parameter_type == TIME_INTERVAL)) && (read_line < offset + 1)) ||
            (((index_parameter_type == POSITION) || (index_parameter_type == POSITION_INTERVAL)) &&
             ((read_line < offset + 1) || ((read_line > offset + 1) && (type[0] == NB_INTERNODE))))) {
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
                      ((index_parameter_type == TIME_INTERVAL) && (int_value <= 0)))  ||
                     ((k == 1) && (type[k - 1] == NB_INTERNODE) && (int_value < 0)))) {
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

      seq = new Sequences(nb_sequence , 0 , length , index_parameter_type ,
                          nb_variable , type);

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
        seq->build_index_parameter_histogram();
      }
      if ((seq->index_parameter_type == TIME) || ((seq->index_parameter_type == POSITION) &&
          (seq->type[0] != NB_INTERNODE))) {
        seq->index_interval_computation();
      }

      for (i = 0;i < nb_variable;i++) {
        seq->min_value_computation(i);
        seq->max_value_computation(i);
        seq->build_marginal_histogram(i);
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

  if (hindex_parameter) {
    os << "   ";
    if (comment_flag) {
      os << "# ";
    }
    os << "(" << SEQ_label[SEQL_MIN_INDEX_PARAMETER] << ": " << hindex_parameter->offset << ", "
       << SEQ_label[SEQL_MAX_INDEX_PARAMETER] << ": " << hindex_parameter->nb_value - 1 << ")" << endl;

    os << "\n";
    if (comment_flag) {
      os << "# ";
    }
    os << (index_parameter_type == TIME ? SEQ_label[SEQL_TIME] : SEQ_label[SEQL_POSITION])
       << " " << STAT_label[STATL_HISTOGRAM] << " - ";
    hindex_parameter->ascii_characteristic_print(os , false , comment_flag);

    if (exhaustive) {
      os << "\n";
      if (comment_flag) {
        os << "# ";
      }
      os << "   | " << (index_parameter_type == TIME ? SEQ_label[SEQL_TIME] : SEQ_label[SEQL_POSITION])
         << " " << STAT_label[STATL_HISTOGRAM] << endl;
      hindex_parameter->ascii_print(os , comment_flag);
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
       << " " << STAT_label[STATL_HISTOGRAM] << " - ";
    index_interval->ascii_characteristic_print(os , false , comment_flag);

    if (exhaustive) {
      os << "\n";
      if (comment_flag) {
        os << "# ";
      }
      os << "   | " << (index_parameter_type == TIME ? SEQ_label[SEQL_TIME_INTERVAL] : SEQ_label[SEQL_POSITION_INTERVAL])
         << " " << STAT_label[STATL_HISTOGRAM] << endl;
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

      if (marginal[i]) {
        os << "\n";
        if (comment_flag) {
          os << "# ";
        }
        os << STAT_label[STATL_MARGINAL] << " " << STAT_label[STATL_HISTOGRAM] << " - ";

        marginal[i]->ascii_characteristic_print(os , exhaustive , comment_flag);

        if ((marginal[i]->nb_value <= ASCII_NB_VALUE) || (exhaustive)) {
          os << "\n";
          if (comment_flag) {
            os << "# ";
          }
          os << "   | " << STAT_label[STATL_MARGINAL] << " " << STAT_label[STATL_HISTOGRAM] << endl;
          marginal[i]->ascii_print(os , comment_flag);
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

        if ((exhaustive) && (variance > 0.)) {
          if (comment_flag) {
            os << "# ";
          }
          os << STAT_label[STATL_SKEWNESS_COEFF] << ": " << skewness_computation(i , mean , variance) << "   "
             << STAT_label[STATL_KURTOSIS_COEFF] << ": " << kurtosis_computation(i , mean , variance) << endl;
        }
      }
    }

    else {
      os << endl;
    }
  }

  os << "\n";
  if (comment_flag) {
    os << "# ";
  }
  os << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_HISTOGRAM] << " - ";
  hlength->ascii_characteristic_print(os , false , comment_flag);

  if (exhaustive) {
    os << "\n";
    if (comment_flag) {
      os << "# ";
    }
    os << "   | " << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_HISTOGRAM] << endl;
    hlength->ascii_print(os , comment_flag);
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
 *  arguments : reference sur un objet Format_error, path,
 *              flag niveau de detail.
 *
 *--------------------------------------------------------------*/

bool Sequences::ascii_write(Format_error &error , const char *path ,
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
 *  arguments : stream, format ('c' : column / 'l' : line / 'a' : array),
 *              flag commentaire, probabilites a posteriori des sequences
 *              d'etats les plus probables (modeles markoviens caches),
 *              nombre de caracteres par ligne.
 *
 *--------------------------------------------------------------*/

ostream& Sequences::ascii_print(ostream &os , char format , bool comment_flag ,
                                double *posterior_probability , int line_nb_character) const

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
      std::ostringstream sos;

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

      if (posterior_probability) {
        if (comment_flag) {
          os << "# ";
        }
        os << SEQ_label[SEQL_POSTERIOR_STATE_SEQUENCE_PROBABILITY]
           << ": " << posterior_probability[i] << endl;
      }
    }
    break;
  }

  case 'l' : {
    int buff , start , width;
    long old_adjust;


    old_adjust = os.setf(ios::right , ios::adjustfield);

    if (index_parameter) {
      width = column_width(hindex_parameter->nb_value - 1);
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

      if (posterior_probability) {
        if (comment_flag) {
          os << "# ";
        }
        os << SEQ_label[SEQL_POSTERIOR_STATE_SEQUENCE_PROBABILITY]
           << ": " << posterior_probability[i] << endl;
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
      std::ostringstream sos;

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
 *  arguments : reference sur un objet Format_error, path,
 *              format (ligne/colonne), flag niveau de detail.
 *
 *--------------------------------------------------------------*/

bool Sequences::ascii_data_write(Format_error &error , const char *path ,
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
 *  arguments : reference sur un objet Format_error, path.
 *
 *--------------------------------------------------------------*/

bool Sequences::spreadsheet_write(Format_error &error , const char *path) const

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

    if (hindex_parameter) {
      out_file << "\t\t" << SEQ_label[SEQL_MIN_INDEX_PARAMETER] << "\t" << hindex_parameter->offset
               << "\t\t" << SEQ_label[SEQL_MAX_INDEX_PARAMETER] << "\t" << hindex_parameter->nb_value - 1 << endl;

      out_file << "\n" << (index_parameter_type == TIME ? SEQ_label[SEQL_TIME] : SEQ_label[SEQL_POSITION])
               << " " << STAT_label[STATL_HISTOGRAM] << "\t";
      hindex_parameter->spreadsheet_characteristic_print(out_file);

      out_file << "\n\t" << (index_parameter_type == TIME ? SEQ_label[SEQL_TIME] : SEQ_label[SEQL_POSITION])
               << " " << STAT_label[STATL_HISTOGRAM] << endl;
      hindex_parameter->spreadsheet_print(out_file);
    }

    else {
      out_file << endl;
    }

    if (index_interval) {
      out_file << "\n" << (index_parameter_type == TIME ? SEQ_label[SEQL_TIME_INTERVAL] : SEQ_label[SEQL_POSITION_INTERVAL])
               << " " << STAT_label[STATL_HISTOGRAM] << "\t";
      index_interval->spreadsheet_characteristic_print(out_file);

      out_file << "\n\t" << (index_parameter_type == TIME ? SEQ_label[SEQL_TIME_INTERVAL] : SEQ_label[SEQL_POSITION_INTERVAL])
               << " " << STAT_label[STATL_HISTOGRAM] << endl;
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

        if (marginal[i]) {
          out_file << "\n" << STAT_label[STATL_MARGINAL] << " " << STAT_label[STATL_HISTOGRAM] << "\t";
          marginal[i]->spreadsheet_characteristic_print(out_file);

          out_file << "\n\t" << STAT_label[STATL_MARGINAL] << " " << STAT_label[STATL_HISTOGRAM] << endl;
          marginal[i]->spreadsheet_print(out_file);
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
        }
      }

      else {
        out_file << endl;
      }
    }

    out_file << "\n" << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_HISTOGRAM] << "\t";
    hlength->spreadsheet_characteristic_print(out_file);

    out_file << "\n\t" << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_HISTOGRAM] << endl;
    hlength->spreadsheet_print(out_file);

    out_file << "\n" << SEQ_label[SEQL_CUMUL_LENGTH] << "\t" << cumul_length << endl;
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Sortie Gnuplot d'un objet Sequences.
 *
 *  arguments : reference sur un objet Format_error, prefixe des fichiers,
 *              titre des figures.
 *
 *--------------------------------------------------------------*/

bool Sequences::plot_write(Format_error &error , const char *prefix ,
                           const char *title) const

{
  bool status;
  register int i , j , k;
  int nb_histo;
  const Histogram **phisto;
  ostringstream data_file_name;


  error.init();

  // ecriture du fichier de donnees

  data_file_name << prefix << ".dat";

  phisto = new const Histogram*[nb_variable + 2];

  nb_histo = 0;

  if (hindex_parameter) {
    phisto[nb_histo++] = hindex_parameter;
  }
  if (index_interval) {
    phisto[nb_histo++] = index_interval;
  }

  for (i = 0;i < nb_variable;i++) {
    if (marginal[i]) {
      phisto[nb_histo++] = marginal[i];
    }
  }

  status = hlength->plot_print((data_file_name.str()).c_str() , nb_histo , phisto);

  delete [] phisto;

  if (!status) {
    error.update(STAT_error[STATR_FILE_PREFIX]);
  }

  // ecriture du fichier de commandes et du fichier d'impression

  else {
    for (i = 0;i < 2;i++) {
      j = 2;

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

      if (hindex_parameter) {
        if (hindex_parameter->nb_value - 1 < TIC_THRESHOLD) {
          out_file << "set xtics 0,1" << endl;
        }
        if ((int)(hindex_parameter->max * YSCALE) + 1 < TIC_THRESHOLD) {
          out_file << "set ytics 0,1" << endl;
        }

        out_file << "plot [" << hindex_parameter->offset << ":"
                 << hindex_parameter->nb_value - 1 << "] [0:"
                 << (int)(hindex_parameter->max * YSCALE) + 1 << "] \""
                 << label((data_file_name.str()).c_str()) << "\" using " << j++ << " title \""
                 << (index_parameter_type == TIME ? SEQ_label[SEQL_TIME] : SEQ_label[SEQL_POSITION])
                 << " " << STAT_label[STATL_HISTOGRAM] << "\" with impulses" << endl;

        if (hindex_parameter->nb_value - 1 < TIC_THRESHOLD) {
          out_file << "set xtics autofreq" << endl;
        }
        if ((int)(hindex_parameter->max * YSCALE) + 1 < TIC_THRESHOLD) {
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
                 << label((data_file_name.str()).c_str()) << "\" using " << j++ << " title \""
                 << (index_parameter_type == TIME ? SEQ_label[SEQL_TIME_INTERVAL] : SEQ_label[SEQL_POSITION_INTERVAL])
                 << " " << STAT_label[STATL_HISTOGRAM] << "\" with impulses" << endl;

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

      for (k = 0;k < nb_variable;k++) {
        if (marginal[k]) {
          if (marginal[k]->nb_value - 1 < TIC_THRESHOLD) {
            out_file << "set xtics 0,1" << endl;
          }
          if ((int)(marginal[k]->max * YSCALE) + 1 < TIC_THRESHOLD) {
            out_file << "set ytics 0,1" << endl;
          }

          out_file << "plot [0:" << MAX(marginal[k]->nb_value - 1 , 1) << "] [0:"
                   << (int)(marginal[k]->max * YSCALE) + 1 << "] \""
                   << label((data_file_name.str()).c_str()) << "\" using " << j++ << " title \"";
          if (nb_variable > 1) {
            out_file << STAT_label[STATL_VARIABLE] << " " << k + 1 << " - ";
          }
          out_file << STAT_label[STATL_MARGINAL] << " " << STAT_label[STATL_HISTOGRAM]
                   << "\" with impulses" << endl;

          if (marginal[k]->nb_value - 1 < TIC_THRESHOLD) {
            out_file << "set xtics autofreq" << endl;
          }
          if ((int)(marginal[k]->max * YSCALE) + 1 < TIC_THRESHOLD) {
            out_file << "set ytics autofreq" << endl;
          }

          if (i == 0) {
            out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
          }
          out_file << endl;
        }
      }

      if (hlength->nb_value - 1 < TIC_THRESHOLD) {
        out_file << "set xtics 0,1" << endl;
      }
      if ((int)(hlength->max * YSCALE) + 1 < TIC_THRESHOLD) {
        out_file << "set ytics 0,1" << endl;
      }

      out_file << "plot [0:" << hlength->nb_value - 1 << "] [0:"
               << (int)(hlength->max * YSCALE) + 1 << "] \""
               << label((data_file_name.str()).c_str()) << "\" using 1 title \""
               << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_HISTOGRAM]
               << "\" with impulses" << endl;

      if (hlength->nb_value - 1 < TIC_THRESHOLD) {
        out_file << "set xtics autofreq" << endl;
      }
      if ((int)(hlength->max * YSCALE) + 1 < TIC_THRESHOLD) {
        out_file << "set ytics autofreq" << endl;
      }

      if (i == 1) {
        out_file << "\nset terminal x11" << endl;
      }

      out_file << "\npause 0 \"" << STAT_label[STATL_END] << "\"" << endl;
    }
  }

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
  if (hindex_parameter) {
    nb_plot_set++;
  }
  if (index_interval) {
    nb_plot_set++;
  }

  for (i = 0;i < nb_variable;i++) {
    if (marginal[i]) {
      nb_plot_set++;
    }
  }

  plot_set = new MultiPlotSet(nb_plot_set);

  MultiPlotSet &plot = *plot_set;

  plot.border = "15 lw 0";

  i = 0;

  if (hindex_parameter) {

    // vue : histogramme des parametres d'index

    plot[i].xrange = Range(hindex_parameter->offset , hindex_parameter->nb_value - 1);
    plot[i].yrange = Range(0 , ceil(hindex_parameter->max * YSCALE));

    if (hindex_parameter->nb_value - 1 < TIC_THRESHOLD) {
      plot[i].xtics = 1;
    }
    if (ceil(hindex_parameter->max * YSCALE) < TIC_THRESHOLD) {
      plot[i].ytics = 1;
    }

    plot[i].resize(1);

    legend.str("");
    legend << (index_parameter_type == TIME ? SEQ_label[SEQL_TIME] : SEQ_label[SEQL_POSITION])
           << " " << STAT_label[STATL_HISTOGRAM];
    plot[i][0].legend = legend.str();

    plot[i][0].style = "impulses";

    hindex_parameter->plotable_frequency_write(plot[i][0]);
    i++;
  }

  if (index_interval) {

    // vue : histogramme des intervalles entre parametres d'index successifs

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
           << " " << STAT_label[STATL_HISTOGRAM];
    plot[i][0].legend = legend.str();

    plot[i][0].style = "impulses";

    index_interval->plotable_frequency_write(plot[i][0]);
    i++;
  }

  for (j = 0;j < nb_variable;j++) {
    if (marginal[j]) {

      // vue : loi marginale empirique

      plot[i].xrange = Range(0 , MAX(marginal[j]->nb_value - 1 , 1));
      plot[i].yrange = Range(0 , ceil(marginal[j]->max * YSCALE));

      if (marginal[j]->nb_value - 1 < TIC_THRESHOLD) {
        plot[i].xtics = 1;
      }
      if (ceil(marginal[j]->max * YSCALE) < TIC_THRESHOLD) {
        plot[i].ytics = 1;
      }

      plot[i].resize(1);

      legend.str("");
      if (nb_variable > 1) {
        legend << STAT_label[STATL_VARIABLE] << " " << j + 1 << " - ";
      }
      legend << STAT_label[STATL_MARGINAL] << " " << STAT_label[STATL_HISTOGRAM];
      plot[i][0].legend = legend.str();

      plot[i][0].style = "impulses";

      marginal[j]->plotable_frequency_write(plot[i][0]);
      i++;
    }
  }

  // vue : histogramme des longueurs des sequences

  plot[i].xrange = Range(0 , hlength->nb_value - 1);
  plot[i].yrange = Range(0 , ceil(hlength->max * YSCALE));

  if (hlength->nb_value - 1 < TIC_THRESHOLD) {
    plot[i].xtics = 1;
  }
  if (ceil(hlength->max * YSCALE) < TIC_THRESHOLD) {
    plot[i].ytics = 1;
  }

  plot[i].resize(1);

  legend.str("");
  legend << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_HISTOGRAM];
  plot[i][0].legend = legend.str();

  plot[i][0].style = "impulses";

  hlength->plotable_frequency_write(plot[i][0]);

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
 *  arguments : reference sur un objet Format_error, prefixe des fichiers,
 *              titre des figures.
 *
 *--------------------------------------------------------------*/

bool Sequences::plot_data_write(Format_error &error , const char *prefix ,
                                const char *title) const

{
  bool status;
  register int i , j , k;
  int min_index_parameter , max_index_parameter , *pfrequency , *length_nb_sequence;
  ostringstream *data_file_name;


  error.init();

  if (nb_sequence > PLOT_NB_SEQUENCE) {
    status = false;
    error.update(SEQ_error[SEQR_NB_SEQUENCE]);
  }

  else {

    // ecriture des fichiers de donnees

    data_file_name = new ostringstream[hlength->nb_value];

    data_file_name[hlength->offset] << prefix << hlength->offset << ".dat";
    status = plot_print((data_file_name[hlength->offset].str()).c_str() , hlength->offset);

    if (!status) {
      error.update(STAT_error[STATR_FILE_PREFIX]);
    }

    else {
      pfrequency = hlength->frequency + hlength->offset + 1;
      for (i = hlength->offset + 1;i < hlength->nb_value;i++) {
        if (*pfrequency++ > 0) {
          data_file_name[i] << prefix << i << ".dat";
          plot_print((data_file_name[i].str()).c_str() , i);
        }
      }

      length_nb_sequence = new int[hlength->nb_value];

      if (index_parameter) {
        min_index_parameter = hindex_parameter->offset;
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
          for (k = 0;k < hlength->nb_value;k++) {
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

          if (index_parameter) {
            if (max_index_parameter - min_index_parameter < TIC_THRESHOLD) {
              out_file << "set xtics 0,1" << endl;
            }
            if (max_value[j] - min_value[j] < TIC_THRESHOLD) {
              out_file << "set ytics " << MIN(min_value[j] , 0) << ",1" << endl;
            }

            out_file << "plot [" << min_index_parameter << ":" << max_index_parameter << "] ["
                     << MIN(min_value[j] , 0) << ":" << MAX(max_value[j] , min_value[j] + 1) << "] ";
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
                  out_file << " title \"" << identifier[k] << "\" with linespoints";
                }
                else {
                  out_file << " notitle with linespoints";
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

            out_file << "plot [0:" << max_length - 1 << "] [" << MIN(min_value[j] , 0)
                     << ":" << MAX(max_value[j] , min_value[j] + 1) << "] ";
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
                  out_file << " title \"" << identifier[k] << "\" with linespoints";
                }
                else {
                  out_file << " notitle with linespoints";
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
 *  argument : reference sur un objet Format_error.
 *
 *--------------------------------------------------------------*/

MultiPlotSet* Sequences::get_plotable_data(Format_error &error) const

{
  register int i , j , k , m , n;
  int nb_plot_set , min_index_parameter , max_index_parameter;
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
      min_index_parameter = hindex_parameter->offset;
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

        plot[i].yrange = Range(MIN(min_value[j] , 0) , MAX(max_value[j] , min_value[j] + 1));
        if (max_value[j] - min_value[j] < TIC_THRESHOLD) {
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
