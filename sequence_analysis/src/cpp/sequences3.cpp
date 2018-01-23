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

#include <string>
#include <sstream>
#include <iomanip>

#include <boost/tokenizer.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/classification.hpp>

#include "stat_tool/stat_label.h"

#include "stat_tool/quantile_computation.hpp"

#include "sequences.h"
#include "sequence_label.h"

using namespace std;
using namespace boost;
using namespace stat_tool;


namespace sequence_analysis {



/*--------------------------------------------------------------*/
/**
 *  \brief Construction of a Sequences object from a file.
 *
 *  \param[in] error      reference on a StatError object,
 *  \param[in] path       file path,
 *  \param[in] old_format flag format.
 *
 *  \return               Sequences object.
 */
/*--------------------------------------------------------------*/

Sequences* Sequences::ascii_read(StatError &error , const string path , bool old_format)

{
  string buffer , trimmed_buffer;
  size_t position;
  typedef tokenizer<char_separator<char>> tokenizer;
  char_separator<char> separator(" \t");
  bool status , lstatus;
  int i , j , k , m;
  int line , read_line , offset , initial_nb_line , max_length , nb_variable = 0 ,
      vector_size , nb_sequence , index , int_value , line_continue , *length;
  variable_nature *type;
  index_parameter_type index_param_type = IMPLICIT_TYPE;
  double real_value;
  Sequences *seq;
  ifstream in_file(path.c_str());


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

    // 1st pass: analysis of the optional line defining the index parameter and
    // of the mandatory line defining the number of variables

    read_line = 0;

    while (getline(in_file , buffer)) {
      line++;

      position = buffer.find('#');
      if (position != string::npos) {
        buffer.erase(position);
      }
      i = 0;

      tokenizer tok_buffer(buffer , separator);

      for (tokenizer::iterator token = tok_buffer.begin();token != tok_buffer.end();token++) {
        switch (i) {

        case 0 : {

          // test INDEX_PARAMETER keyword

          if ((!old_format) && (read_line == 0) && (*token == SEQ_word[SEQW_INDEX_PARAMETER])) {
            index_param_type = TIME;
          }

          // test number of variables

          else {
            lstatus = true;

/*            try {
              int_value = stoi(*token);   in C++ 11
            }
            catch(invalid_argument &arg) {
              lstatus = false;
            } */
            int_value = atoi(token->c_str());

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

          // test separator

          if ((!old_format) && (read_line == 0) && (index_param_type != IMPLICIT_TYPE)) {
            if (*token != ":") {
              status = false;
              error.update(STAT_parsing[STATP_SEPARATOR] , line , i + 1);
            }
          }

          // test VARIABLE(S) keyword

          else if (*token != STAT_word[nb_variable == 1 ? STATW_VARIABLE : STATW_VARIABLES]) {
            status = false;
            error.correction_update(STAT_parsing[STATP_KEYWORD] ,
                                    STAT_word[nb_variable == 1 ? STATW_VARIABLE : STATW_VARIABLES] , line , i + 1);
          }
          break;
        }

        // test keyword defining the index parameter type

        case 2 : {
          if ((!old_format) && (read_line == 0) && (index_param_type != IMPLICIT_TYPE)) {
            for (j = TIME;j <= POSITION_INTERVAL;j++) {
              if (*token == SEQ_index_parameter_word[j]) {
                index_param_type = (index_parameter_type)j;
                break;
              }
            }

            if (j == POSITION_INTERVAL + 1) {
              status = false;
              error.update(STAT_parsing[STATP_KEYWORD] , line , i + 1);
            }
          }
          break;
        }
        }

        i++;
      }

      if (i > 0) {
        if (((!old_format) && (read_line == 0) && (index_param_type != IMPLICIT_TYPE) && (i != 3)) ||
            (((old_format) || (read_line == 1) || (index_param_type == IMPLICIT_TYPE)) && (i != 2))) {
          status = false;
          error.update(STAT_parsing[STATP_FORMAT] , line);
        }

        read_line++;
//        if (((!old_format) && (index_param_type != IMPLICIT_TYPE) && (read_line == 2)) ||
//            (((old_format) || (index_param_type == IMPLICIT_TYPE)) && (read_line == 1)) {
        if ((((old_format) || (index_param_type == IMPLICIT_TYPE)) && (read_line == 1)) ||
            (read_line == 2)) {
          break;
        }
      }
    }

    if (read_line < 1) {
      status = false;
      error.update(STAT_parsing[STATP_FORMAT]);
    }

    // analysis of the lines defining the variable types

    if (status) {
      type = new variable_nature[nb_variable];
      for (i = 0;i < nb_variable;i++) {
        type[i] = AUXILIARY;
      }

      read_line = 0;
      offset = (old_format ? 1 : 0);

      while (getline(in_file , buffer)) {
        line++;

        position = buffer.find('#');
        if (position != string::npos) {
          buffer.erase(position);
        }
        i = 0;

        tokenizer tok_buffer(buffer , separator);

        for (tokenizer::iterator token = tok_buffer.begin();token != tok_buffer.end();token++) {
          switch (i) {

          // test VARIABLE keyword

          case 0 : {
            if (*token != STAT_word[STATW_VARIABLE]) {
              status = false;
              error.correction_update(STAT_parsing[STATP_KEYWORD] , STAT_word[STATW_VARIABLE] , line , i + 1);
            }
            break;
          }

          // test variable index

          case 1 : {
            lstatus = true;

/*            try {
              int_value = stoi(*token);   in C++ 11
            }
            catch(invalid_argument &arg) {
              lstatus = false;
            } */
            int_value = atoi(token->c_str());

            if ((lstatus) && (int_value != read_line + 1)) {
              lstatus = false;
            }

            if (!lstatus) {
              status = false;
              error.correction_update(STAT_parsing[STATP_VARIABLE_INDEX] , read_line + 1 , line , i + 1);
            }
            break;
          }

          // test separator

          case 2 : {
            if (*token != ":") {
              status = false;
              error.update(STAT_parsing[STATP_SEPARATOR] , line , i + 1);
            }
            break;
          }

          // test keyword defining the variable type

          case 3 : {
            if ((old_format) && (read_line == 0)) {
              for (j = TIME;j <= POSITION_INTERVAL;j++) {
                if (*token == SEQ_index_parameter_word[j]) {
                  index_param_type = (index_parameter_type)j;
                  break;
                }
              }

              if (j == POSITION_INTERVAL + 1) {
//                for (j = INT_VALUE;j <= STATE;j++) {
                for (j = INT_VALUE;j <= OLD_INT_VALUE;j++) {
                  if (*token == STAT_variable_word[j]) {
//                    if (j == STATE) {
                    if ((j == STATE) || (j == OLD_INT_VALUE)) {
                      j = INT_VALUE;
                    }
                    type[read_line] = (variable_nature)j;
                    break;
                  }
                }

//                if (j == STATE + 1) {
                if (j == OLD_INT_VALUE + 1) {
                  status = false;
                  error.update(STAT_parsing[STATP_KEYWORD] , line , i + 1);
                }
              }
            }

            else {
//              for (j = INT_VALUE;j <= NB_INTERNODE;j++) {
              for (j = INT_VALUE;j <= OLD_INT_VALUE;j++) {
                if (*token == STAT_variable_word[j]) {
//                  if ((j == NB_INTERNODE) && ((read_line != offset) || ((read_line == offset) &&
                  if ((j == OLD_INT_VALUE) && ((read_line != offset) || ((read_line == offset) &&
                        (index_param_type != POSITION) && (index_param_type != POSITION_INTERVAL)))) {
                    status = false;
                    error.update(STAT_parsing[STATP_VARIABLE_TYPE] , line , i + 1);
                  }

                  else {
//                    if (j == STATE) {
                    if ((j == STATE) || (j == OLD_INT_VALUE)) {
                      j = INT_VALUE;
                    }

                    if ((old_format) && (index_param_type != IMPLICIT_TYPE)) {
                      type[read_line - 1] = (variable_nature)j;
                    }
                    else {
                      type[read_line] = (variable_nature)j;
                    }
                  }
                  break;
                }
              }
            }

//            if (j == NB_INTERNODE + 1) {
            if (j == OLD_INT_VALUE + 1) {
              status = false;
              error.update(STAT_parsing[STATP_KEYWORD] , line , i + 1);
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
//        if ((((index_param_type == TIME) || (index_param_type == TIME_INTERVAL)) && (read_line < offset + 1)) ||
//            (((index_param_type == POSITION) || (index_param_type == POSITION_INTERVAL)) &&
//             ((read_line < offset + 1) || ((read_line > offset + 1) && (type[0] == NB_INTERNODE))))) {
        if (((index_param_type == TIME) || (index_param_type == TIME_INTERVAL) ||
             (index_param_type == POSITION) || (index_param_type == POSITION_INTERVAL)) &&
             (read_line < offset + 1)) {
          status = false;
          error.update(STAT_parsing[STATP_VARIABLE_TYPE]);
        }
      }

      initial_nb_line = line;
    }

    if (status) {
      vector_size = nb_variable;

      if (index_param_type != IMPLICIT_TYPE) {
        if (old_format) {
          nb_variable--;
        }
        else {
          vector_size++;
        }
      }

      nb_sequence = 0;
      lstatus = true;

//      while (buffer.readLine(in_file , true)) {
      while (getline(in_file , buffer)) {
        position = buffer.find('#');
        if (position != string::npos) {
          buffer.erase(position);
        }

        if (!(buffer.empty())) {
          trimmed_buffer = trim_right_copy_if(buffer , is_any_of(" \t"));

          if ((!(trimmed_buffer.empty())) && (trimmed_buffer.find('\\' , trimmed_buffer.length() - 1) == string::npos)) {
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

    // 2nd pass: analysis of the sequence format

    if (status) {
//      in_file.close();
//      in_file.open(path.c_str() , ios::in);

      in_file.clear();
      in_file.seekg(0 , ios::beg);

      offset = (index_param_type == IMPLICIT_TYPE ? 0 : 1);

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
        getline(in_file , buffer);
        line++;

#       ifdef DEBUG
        cout << line << "  " << buffer << endl;
#       endif

      }
      while (line < initial_nb_line);

      max_length = 0;

      switch (index_param_type) {
      case TIME :
        index = -1;
        break;
      case POSITION :
        index = 0;
        break;
      }

      i = 0;

      while (getline(in_file , buffer)) {
        line++;

#       ifdef DEBUG
        cout << line << "  " << buffer << endl;
#       endif

        position = buffer.find('#');
        if (position != string::npos) {
          buffer.erase(position);
        }

        j = 0;
        k = 0;
        line_continue = false;

        tokenizer tok_buffer(buffer , separator);

        for (tokenizer::iterator token = tok_buffer.begin();token != tok_buffer.end();token++) {
          if (line_continue) {
            status = false;
            error.update(STAT_parsing[STATP_FORMAT] , line , j);
            break;
          }

          if ((vector_size > 1) && (j % (vector_size + 1) == vector_size)) {
            if (*token == "\\") {
              line_continue = true;
              length[i]++;
            }

            else {
              if (*token != "|") {
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
            if ((vector_size == 1) && (*token == "\\")) {
              line_continue = true;
            }

            else {
              lstatus = true;

              if (((index_param_type != IMPLICIT_TYPE) && (k == 0)) || (type[k - offset] != REAL_VALUE)) {
/*                try {
                  int_value = stoi(*token);   in C++ 11
                }
                catch(invalid_argument &arg) {
                  lstatus = false;
                } */
                int_value = atoi(token->c_str());

                if ((lstatus) && (((k == 0) && (((index_param_type == TIME) || (index_param_type == POSITION) ||
                        (index_param_type == POSITION_INTERVAL)) && (int_value < 0)) ||
                      ((index_param_type == TIME_INTERVAL) && (int_value <= 0))))) {
//                     ((k == 1) && (type[k - 1] == NB_INTERNODE) && (int_value < 0)))) {
                  lstatus = false;
                }
              }

              else {
/*                try {
                  real_value = stod(*token);   in C++ 11
                }
                catch(invalid_argument &arg) {
                  lstatus = false;
                } */
                real_value = atof(token->c_str());
              }

              if (!lstatus) {
                status = false;
                error.update(STAT_parsing[STATP_DATA_TYPE] , line , j + 1);
              }

              else if (k == 0) {
                switch (index_param_type) {

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
            if (((line_continue) || ((index_param_type != POSITION) && (index_param_type != POSITION_INTERVAL))) &&
                (k != vector_size)) {
              status = false;
              error.update(STAT_parsing[STATP_FORMAT] , line , j);
            }
          }

          if (!line_continue) {
            if ((index_param_type == POSITION) || (index_param_type == POSITION_INTERVAL)) {
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

            switch (index_param_type) {
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
        cout << i << " " << length[i] << " | ";
      }
      cout << endl;
#     endif
    }

    // 3rd pass: sequence copy

    if (status) {
//      in_file.close();
//      in_file.open(path.c_str() , ios::in);

      in_file.clear();
      in_file.seekg(0 , ios::beg);

      seq = new Sequences(nb_sequence , NULL , length , NULL ,
                          index_param_type , nb_variable , type);

      line = 0;

      do {
        getline(in_file , buffer);
        line++;
      }
      while (line < initial_nb_line);

      i = 0;
      j = 0;

      while (getline(in_file , buffer)) {
        position = buffer.find('#');
        if (position != string::npos) {
          buffer.erase(position);
        }

        k = 0;
        m = 0;
        line_continue = false;

        tokenizer tok_buffer(buffer , separator);

        for (tokenizer::iterator token = tok_buffer.begin();token != tok_buffer.end();token++) {
          if ((vector_size > 1) && (m % (vector_size + 1) == vector_size)) {
            if (*token == "\\") {
              line_continue = true;
            }
            k = 0;
            j++;
          }

          else {
            if ((vector_size == 1) && (*token == "\\")) {
              line_continue = true;
            }

            else {
              if ((index_param_type != IMPLICIT_TYPE) && (k == 0)) {
//                seq->index_parameter[i][j] = stoi(*token);   in C++ 11
                seq->index_parameter[i][j] = atoi(token->c_str());
              }
              else if (type[k - offset] != REAL_VALUE) {
//                seq->int_sequence[i][k - offset][j] = stoi(*token);   in C++ 11
                seq->int_sequence[i][k - offset][j] = atoi(token->c_str());
              }
              else {
//                seq->real_sequence[i][k - offset][j] = stod(*token);   in C++ 11
                seq->real_sequence[i][k - offset][j] = atof(token->c_str());
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

      if ((seq->index_param_type == TIME_INTERVAL) || (seq->index_param_type == POSITION_INTERVAL)) {
        seq->index_parameter_computation();
      }

      if (seq->index_parameter) {
        seq->build_index_parameter_frequency_distribution();
      }
//      if ((seq->index_param_type == TIME) || ((seq->index_param_type == POSITION) &&
//          (seq->type[0] != NB_INTERNODE))) {
      if ((seq->index_param_type == TIME) || (seq->index_param_type == POSITION)) {
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


/*--------------------------------------------------------------*/
/**
 *  \brief Writing on a single line of a Sequences object.
 *
 *  \param[in,out] os stream.
 */
/*--------------------------------------------------------------*/

ostream& Sequences::line_write(ostream &os) const

{
  os << nb_sequence << " " << SEQ_label[nb_sequence == 1 ? SEQL_SEQUENCE : SEQL_SEQUENCES] << "   "
     << nb_variable << " " << STAT_word[nb_variable == 1 ? STATW_VARIABLE : STATW_VARIABLES] << "   "
     << SEQ_label[nb_sequence == 1 ? SEQL_LENGTH : SEQL_CUMUL_LENGTH] << ": " << cumul_length;

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a Sequences object.
 *
 *  \param[in,out] os           stream,
 *  \param[in]     exhaustive   flag detail level,
 *  \param[in]     comment_flag flag comment.
 */
/*--------------------------------------------------------------*/

ostream& Sequences::ascii_write(ostream &os , bool exhaustive , bool comment_flag) const

{
  int i , j , k;
  int *int_value , *pint_value;
  double mean , variance , median , lower_quartile , upper_quartile , *real_value , *preal_value;


  if (index_parameter) {
    os << SEQ_word[SEQW_INDEX_PARAMETER] << " : "
       << SEQ_index_parameter_word[index_param_type];
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
    os << (index_param_type == TIME ? SEQ_label[SEQL_TIME] : SEQ_label[SEQL_POSITION])
       << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";
    index_parameter_distribution->ascii_characteristic_print(os , false , comment_flag);

    if (exhaustive) {
      os << "\n";
      if (comment_flag) {
        os << "# ";
      }
      os << "   | " << (index_param_type == TIME ? SEQ_label[SEQL_TIME] : SEQ_label[SEQL_POSITION])
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
    os << (index_param_type == TIME ? SEQ_label[SEQL_TIME_INTERVAL] : SEQ_label[SEQL_POSITION_INTERVAL])
       << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";
    index_interval->ascii_characteristic_print(os , false , comment_flag);

    if (exhaustive) {
      os << "\n";
      if (comment_flag) {
        os << "# ";
      }
      os << "   | " << (index_param_type == TIME ? SEQ_label[SEQL_TIME_INTERVAL] : SEQ_label[SEQL_POSITION_INTERVAL])
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
        mean = mean_computation(i);
        variance = variance_computation(i , mean);

        if (variance > 0.) {
          switch (type[i]) {

          case INT_VALUE : {
            int_value = new int[cumul_length];
            pint_value = int_value;
            for (j = 0;j < nb_sequence;j++) {
              for (k = 0;k < length[j];k++) {
                *pint_value++ = int_sequence[j][i][k];
              }
            }

            lower_quartile = quantile_computation(cumul_length , int_value , 0.25);
            median = quantile_computation(cumul_length , int_value , 0.5);
            upper_quartile = quantile_computation(cumul_length , int_value , 0.75);

            delete [] int_value;
            break;
          }

          case REAL_VALUE : {
            real_value = new double[cumul_length];
            preal_value = real_value;
            for (j = 0;j < nb_sequence;j++) {
              for (k = 0;k < length[j];k++) {
                *preal_value++ = real_sequence[j][i][k];
              }
            }

            lower_quartile = quantile_computation(cumul_length , real_value , 0.25);
            median = quantile_computation(cumul_length , real_value , 0.5);
            upper_quartile = quantile_computation(cumul_length , real_value , 0.75);

            delete [] real_value;
            break;
          }
          }
        }

        else {
          median = mean;
        }

        os << "\n";
        if (comment_flag) {
          os << "# ";
        }
        os << STAT_label[STATL_SAMPLE_SIZE] << ": " << cumul_length << endl;

        if (comment_flag) {
          os << "# ";
        }
        os << STAT_label[STATL_MEAN] << ": " << mean << "   "
           << STAT_label[STATL_MEDIAN] << ": " << median << endl;

        if (comment_flag) {
          os << "# ";
        }
        os << STAT_label[STATL_VARIANCE] << ": " << variance << "   "
           << STAT_label[STATL_STANDARD_DEVIATION] << ": " << sqrt(variance);
        if (variance > 0.) {
          os << "   " << STAT_label[STATL_LOWER_QUARTILE] << ": " << lower_quartile
             << "   " << STAT_label[STATL_UPPER_QUARTILE] << ": " << upper_quartile;
        }
        os << endl;

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


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a Sequences object.
 *
 *  \param[in,out] os         stream,
 *  \param[in]     exhaustive flag detail level.
 */
/*--------------------------------------------------------------*/

ostream& Sequences::ascii_write(ostream &os , bool exhaustive) const

{
  return ascii_write(os , exhaustive , false);
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a Sequences object in a file.
 *
 *  \param[in] error      reference on a StatError object,
 *  \param[in] path       file path,
 *  \param[in] exhaustive flag detail level.
 *
 *  \return               error status.
 */
/*--------------------------------------------------------------*/

bool Sequences::ascii_write(StatError &error , const string path ,
                            bool exhaustive) const

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
    ascii_write(out_file , exhaustive , false);
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of sequences.
 *
 *  \param[in,out] os                          stream,
 *  \param[in]     format                      format (LINE/COLUMN/VECTOR/POSTERIOR_PROBABILITY),
 *  \param[in]     comment_flag                flag comment,
 *  \param[in]     posterior_probability       posterior probabilities of the most probable state sequences,
 *  \param[in]     entropy                     entropies of state sequences,
 *  \param[in]     nb_state_sequence           numbers of state sequences (hidden Markovian models),
 *  \param[in]     posterior_state_probability posterior probabilities of the most probable initial state,
 *  \param[in]     line_nb_character           number of characters per line.
 */
/*--------------------------------------------------------------*/

ostream& Sequences::ascii_print(ostream &os , output_sequence_format format , bool comment_flag ,
                                double *posterior_probability , double *entropy ,
                                double *nb_state_sequence , double *posterior_state_probability ,
                                int line_nb_character) const

{
  int i , j , k , m;


  switch (format) {

  case COLUMN : {
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

      if (index_param_type == POSITION) {
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

  case VECTOR : {
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

      if (index_param_type == POSITION) {
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

  case LINE : {
    int buff , start , width;
    ios_base::fmtflags format_flags;


    format_flags = os.setf(ios::right , ios::adjustfield);

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
          if ((j - start) * (width + 1) > 10000) {
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
          if (index_param_type == POSITION) {
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

    os.setf(format_flags , ios::adjustfield);
    break;
  }

  case ARRAY : {
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

      if (index_param_type == POSITION) {
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

  case POSTERIOR_PROBABILITY : {
    if ((posterior_probability) && (entropy) && (nb_state_sequence)) {
      bool *selected_sequence;
      int index , width[7];
      double max , *divergence;
      ios_base::fmtflags format_flags;


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
      if (posterior_state_probability) {
        width[6] = column_width(nb_sequence , posterior_state_probability) + ASCII_SPACE;
      }

      selected_sequence = new bool[nb_sequence];
      for (i = 0;i < nb_sequence;i++) {
        selected_sequence[i] = false;
      }

      format_flags = os.setf(ios::left , ios::adjustfield);

      os << "\n" << (posterior_state_probability ? 7 : 6) << " " << STAT_word[STATW_VARIABLES] << endl;

      i = 1;
      os << "\n" << STAT_word[STATW_VARIABLE] << " " << i++ << " : " << STAT_variable_word[INT_VALUE] << endl;

      if (posterior_state_probability) {
        os << STAT_word[STATW_VARIABLE] << " " << i++ << " : " << STAT_variable_word[REAL_VALUE] << "   ";
        if (comment_flag) {
          os << "# ";
        }
        os << SEQ_label[SEQL_POSTERIOR_INITIAL_STATE_PROBABILITY] << endl;
      }

      os << STAT_word[STATW_VARIABLE] << " " << i++ << " : " << STAT_variable_word[REAL_VALUE] << "   ";
      if (comment_flag) {
        os << "# ";
      }
      os << SEQ_label[SEQL_POSTERIOR_STATE_SEQUENCE_PROBABILITY] << endl;

      os << STAT_word[STATW_VARIABLE] << " " << i++ << " : " << STAT_variable_word[REAL_VALUE] << "   ";
      if (comment_flag) {
        os << "# ";
      }
      os << SEQ_label[SEQL_STATE_SEQUENCE_ENTROPY] << endl;

      os << STAT_word[STATW_VARIABLE] << " " << i++ << " : " << STAT_variable_word[REAL_VALUE] << "   ";
      if (comment_flag) {
        os << "# ";
      }
      os << SEQ_label[SEQL_STATE_SEQUENCE_DIVERGENCE] << endl;

      os << STAT_word[STATW_VARIABLE] << " " << i++ << " : " << STAT_variable_word[REAL_VALUE] << "   ";
      if (comment_flag) {
        os << "# ";
      }
      os << SEQ_label[SEQL_NB_STATE_SEQUENCE] << endl;

      os << STAT_word[STATW_VARIABLE] << " " << i << " : " << STAT_variable_word[INT_VALUE] << "    ";
      if (comment_flag) {
        os << "# ";
      }
      os << SEQ_label[SEQL_SEQUENCE_LENGTH] << endl;

      os << "\n";
      for (i = 0;i < nb_sequence;i++) {

#       ifdef DEBUG
        for (j = 0;j < length[i];j++) {  // for Fuji/Braeburn GUs
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

        if (posterior_state_probability) {
          for (j = 0;j < nb_sequence;j++) {
            if ((!selected_sequence[j]) && (posterior_state_probability[j] > max)) {
              max = posterior_state_probability[j];
              index = j;
            }
          }
        }

        else {
          for (j = 0;j < nb_sequence;j++) {
/*            if ((!selected_sequence[j]) && (entropy[j] > max)) {
              max = entropy[j]; */
            if ((!selected_sequence[j]) && (posterior_probability[j] > max)) {
              max = posterior_probability[j];
              index = j;
            }
          }
        }

        selected_sequence[index] = true;

        os << setw(width[0]) << i + 1;

        if (posterior_state_probability) {

#         ifdef MESSAGE
          if (posterior_probability[index] > posterior_state_probability[index] + DOUBLE_ERROR) {
            cout << "\n" << SEQ_label[SEQL_SEQUENCE] << " " << identifier[index] << ", "<< SEQ_error[SEQR_POSTERIOR_PROBABILITY_ORDER] << endl;
          }
#         endif

          os << setw(width[6]) << posterior_state_probability[index];
        }

        os << setw(width[1]) << posterior_probability[index]
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

      os.setf(format_flags , ios::adjustfield);
    }
    break;
  }
  }

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a Sequences object.
 *
 *  \param[in,out] os         stream,
 *  \param[in]     format     format (LINE/COLUMN/VECTOR/POSTERIOR_PROBABILITY),
 *  \param[in]     exhaustive flag detail level.
 */
/*--------------------------------------------------------------*/

ostream& Sequences::ascii_data_write(ostream &os , output_sequence_format format ,
                                     bool exhaustive) const

{
  ascii_write(os , exhaustive , false);
  ascii_print(os , format , false);

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a Sequences object in a file.
 *
 *  \param[in] error      reference on a StatError object,
 *  \param[in] path       file path,
 *  \param[in] format     format (LINE/COLUMN/VECTOR/POSTERIOR_PROBABILITY),
 *  \param[in] exhaustive flag detail level.
 *
 *  \return               error status.
 */
/*--------------------------------------------------------------*/

bool Sequences::ascii_data_write(StatError &error , const string path ,
                                 output_sequence_format format , bool exhaustive) const

{
  bool status = false;
  ofstream out_file(path.c_str());


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


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a Sequences object in a file at the spreadsheet format.
 *
 *  \param[in] error reference on a StatError object,
 *  \param[in] path  file path.
 *
 *  \return          error status.
 */
/*--------------------------------------------------------------*/

bool Sequences::spreadsheet_write(StatError &error , const string path) const

{
  bool status;
  int i , j , k;
  int *int_value , *pint_value;
  double mean , variance , median , lower_quartile , upper_quartile , *real_value , *preal_value;
  ofstream out_file(path.c_str());


  error.init();

  if (!out_file) {
    status = false;
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {
    status = true;

    if (index_parameter) {
      out_file << SEQ_word[SEQW_INDEX_PARAMETER] << "\t"
               << SEQ_index_parameter_word[index_param_type];
    }

    if (index_parameter_distribution) {
      out_file << "\t\t" << SEQ_label[SEQL_MIN_INDEX_PARAMETER] << "\t" << index_parameter_distribution->offset
               << "\t\t" << SEQ_label[SEQL_MAX_INDEX_PARAMETER] << "\t" << index_parameter_distribution->nb_value - 1 << endl;

      out_file << "\n" << (index_param_type == TIME ? SEQ_label[SEQL_TIME] : SEQ_label[SEQL_POSITION])
               << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\t";
      index_parameter_distribution->spreadsheet_characteristic_print(out_file);

      out_file << "\n\t" << (index_param_type == TIME ? SEQ_label[SEQL_TIME] : SEQ_label[SEQL_POSITION])
               << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << endl;
      index_parameter_distribution->spreadsheet_print(out_file);
    }

    else {
      out_file << endl;
    }

    if (index_interval) {
      out_file << "\n" << (index_param_type == TIME ? SEQ_label[SEQL_TIME_INTERVAL] : SEQ_label[SEQL_POSITION_INTERVAL])
               << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\t";
      index_interval->spreadsheet_characteristic_print(out_file);

      out_file << "\n\t" << (index_param_type == TIME ? SEQ_label[SEQL_TIME_INTERVAL] : SEQ_label[SEQL_POSITION_INTERVAL])
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
          mean = mean_computation(i);
          variance = variance_computation(i , mean);

          if (variance > 0.) {
            switch (type[i]) {

            case INT_VALUE : {
              int_value = new int[cumul_length];
              pint_value = int_value;
              for (j = 0;j < nb_sequence;j++) {
                for (k = 0;k < length[j];k++) {
                  *pint_value++ = int_sequence[j][i][k];
                }
              }

              lower_quartile = quantile_computation(cumul_length , int_value , 0.25);
              median = quantile_computation(cumul_length , int_value , 0.5);
              upper_quartile = quantile_computation(cumul_length , int_value , 0.75);

              delete [] int_value;
              break;
            }

            case REAL_VALUE : {
              real_value = new double[cumul_length];
              preal_value = real_value;
              for (j = 0;j < nb_sequence;j++) {
                for (k = 0;k < length[j];k++) {
                  *preal_value++ = real_sequence[j][i][k];
                }
              }

              lower_quartile = quantile_computation(cumul_length , real_value , 0.25);
              median = quantile_computation(cumul_length , real_value , 0.5);
              upper_quartile = quantile_computation(cumul_length , real_value , 0.75);

              delete [] real_value;
              break;
            }
            }
          }

          else {
            median = mean;
          }

          out_file << "\n" << STAT_label[STATL_SAMPLE_SIZE] << "\t" << cumul_length << endl;

          out_file << STAT_label[STATL_MEAN] << "\t" << mean << "\t\t"
                   << STAT_label[STATL_MEDIAN] << "\t" << median << endl;

          out_file << STAT_label[STATL_VARIANCE] << "\t" << variance << "\t\t"
                   << STAT_label[STATL_STANDARD_DEVIATION] << "\t" << sqrt(variance);
          if (variance > 0.) {
            out_file << "\t\t" << STAT_label[STATL_LOWER_QUARTILE] << "\t" << lower_quartile
                     << "\t\t" << STAT_label[STATL_UPPER_QUARTILE] << "\t" << upper_quartile;
          }
          out_file << endl;

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


/*--------------------------------------------------------------*/
/**
 *  \brief Plot of a Sequences object using Gnuplot.
 *
 *  \param[in] error  reference on a StatError object,
 *  \param[in] prefix file prefix,
 *  \param[in] title  figure title.
 *
 *  \return           error status.
 */
/*--------------------------------------------------------------*/

bool Sequences::plot_write(StatError &error , const char *prefix ,
                           const char *title) const

{
  bool status;
  int i , j;
  int nb_histo;
  const FrequencyDistribution *phisto[2];
  ostringstream *data_file_name;


  error.init();

  // writing of the data files

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

    // writing of the script files

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
                 << (index_param_type == TIME ? SEQ_label[SEQL_TIME] : SEQ_label[SEQL_POSITION])
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
                 << (index_param_type == TIME ? SEQ_label[SEQL_TIME_INTERVAL] : SEQ_label[SEQL_POSITION_INTERVAL])
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

          out_file << "plot [" << marginal_histogram[j]->min_value - marginal_histogram[j]->bin_width << ":"
                   << marginal_histogram[j]->max_value + marginal_histogram[j]->bin_width << "] [0:"
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


/*--------------------------------------------------------------*/
/**
 *  \brief Plot of a Sequences object.
 *
 *  \return MultiPlotSet object.
 */
/*--------------------------------------------------------------*/

MultiPlotSet* Sequences::get_plotable() const

{
  int i , j;
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

    // index parameter frequency distribution

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
    legend << (index_param_type == TIME ? SEQ_label[SEQL_TIME] : SEQ_label[SEQL_POSITION])
           << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
    plot[i][0].legend = legend.str();

    plot[i][0].style = "impulses";

    index_parameter_distribution->plotable_frequency_write(plot[i][0]);
    i++;
  }

  if (index_interval) {

    // between-index interval frequency distribution

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
    legend << (index_param_type == TIME ? SEQ_label[SEQL_TIME_INTERVAL] : SEQ_label[SEQL_POSITION_INTERVAL])
           << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
    plot[i][0].legend = legend.str();

    plot[i][0].style = "impulses";

    index_interval->plotable_frequency_write(plot[i][0]);
    i++;
  }

  for (j = 0;j < nb_variable;j++) {
    if (marginal_distribution[j]) {

      // marginal frequency distribution

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

      // marginal histogram

      plot[i].xrange = Range(marginal_histogram[j]->min_value - marginal_histogram[j]->bin_width ,
                             marginal_histogram[j]->max_value + marginal_histogram[j]->bin_width);
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

  // sequence length frequency distribution

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


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of sequences at the Gnuplot format.
 *
 *  \param[in] path    file path,
 *  \param[in] ilength sequence length.
 *
 *  \return            error status.
 */
/*--------------------------------------------------------------*/

bool Sequences::plot_print(const char *path , int ilength) const

{
  bool status = false;
  int i , j , k;
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


/*--------------------------------------------------------------*/
/**
 *  \brief Plot of a Sequences object using Gnuplot.
 *
 *  \param[in] error  reference on a StatError object,
 *  \param[in] prefix file prefix,
 *  \param[in] title  figure title.
 *
 *  \return           error status.
 */
/*--------------------------------------------------------------*/

bool Sequences::plot_data_write(StatError &error , const char *prefix ,
                                const char *title) const

{
  bool status;
  int i , j , k;
  int min_index_parameter , max_index_parameter , *pfrequency , *length_nb_sequence;
  double min , max;
  ostringstream *data_file_name;


  error.init();

  if (nb_sequence > PLOT_NB_SEQUENCE) {
    status = false;
    error.update(SEQ_error[SEQR_NB_SEQUENCE]);
  }

  else {

    // writing of the data file

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

      // writing of the script files

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

          if (max == min) {
            max = min + 1;
          }
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
                     << MIN(min , 0) << ":" << (max >= 0. ? max * YSCALE : max * (2. - YSCALE)) << "] ";
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
                     << ":" << (max >= 0. ? max * YSCALE : max * (2. - YSCALE)) << "] ";
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


/*--------------------------------------------------------------*/
/**
 *  \brief Plot of a Sequences object.
 *
 *  \param[in] error reference on a StatError object.
 *
 *  \return          MultiPlotSet object.
 */
/*--------------------------------------------------------------*/

MultiPlotSet* Sequences::get_plotable_data(StatError &error) const

{
  int i , j , k , m , n;
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

        if (max == min) {
          max = min + 1;
        }
        if ((j + 1 < nb_variable) && (type[j + 1] == AUXILIARY)) {
          if (min_value[j + 1] < min) {
            min = min_value[j + 1];
          }
          if (max_value[j + 1] > max) {
            max = max_value[j + 1];
          }
        }

        plot[i].yrange = Range(MIN(min , 0) , (max >= 0. ? max * YSCALE : max * (2. - YSCALE)));
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


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the maximum length of sequences.
 */
/*--------------------------------------------------------------*/

void Sequences::max_length_computation()

{
  int i;


  max_length = length[0];
  for (i = 1;i < nb_sequence;i++) {
    if (length[i] > max_length) {
      max_length = length[i];
    }
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the cumulative length of sequences.
 */
/*--------------------------------------------------------------*/

void Sequences::cumul_length_computation()

{
  int i;


  cumul_length = 0;
  for (i = 0;i < nb_sequence;i++) {
    cumul_length += length[i];
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Construction of the sequence length frequency distribution.
 */
/*--------------------------------------------------------------*/

void Sequences::build_length_frequency_distribution()

{
  int i;


  length_distribution = new FrequencyDistribution(max_length + 1);

  length_distribution->nb_element = nb_sequence;
  for (i = 0;i < nb_sequence;i++) {
    (length_distribution->frequency[length[i]])++;
  }

  length_distribution->nb_value_computation();
  length_distribution->offset_computation();
  length_distribution->max_computation();
  length_distribution->mean_computation();
  length_distribution->variance_computation();
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of index parameters from between-index intervals.
 */
/*--------------------------------------------------------------*/

void Sequences::index_parameter_computation()

{
  if ((index_param_type == TIME_INTERVAL) || (index_param_type == POSITION_INTERVAL)) {
    int i , j;


    switch (index_param_type) {
    case TIME_INTERVAL :
      index_param_type = TIME;
      break;
    case POSITION_INTERVAL :
      index_param_type = POSITION;
      break;
    }

    for (i = 0;i < nb_sequence;i++) {
      for (j = 1;j < (index_param_type == POSITION ? length[i] + 1 : length[i]);j++) {
        index_parameter[i][j] += index_parameter[i][j - 1];
      }
    }
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the minimum value taken by the index parameter.
 */
/*--------------------------------------------------------------*/

int Sequences::min_index_parameter_computation() const

{
  int i;
  int min_index_parameter = I_DEFAULT;


  if (index_parameter) {
    min_index_parameter = index_parameter[0][0];
    for (i = 1;i < nb_sequence;i++) {
      if (index_parameter[i][0] < min_index_parameter) {
        min_index_parameter = index_parameter[i][0];
      }
    }
  }

  return min_index_parameter;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the maximum value taken by the index parameter.
 *
 *  \param[in] last_position flag last position.
 */
/*--------------------------------------------------------------*/

int Sequences::max_index_parameter_computation(bool last_position) const

{
  int i;
  int max_index_parameter = I_DEFAULT;


  if (index_parameter) {
    if ((index_param_type == TIME) || (last_position)) {
      max_index_parameter = index_parameter[0][length[0] - 1];
      for (i = 1;i < nb_sequence;i++) {
        if (index_parameter[i][length[i] - 1] > max_index_parameter) {
          max_index_parameter = index_parameter[i][length[i] - 1];
        }
      }
    }

    else {
      max_index_parameter = index_parameter[0][length[0]];
      for (i = 1;i < nb_sequence;i++) {
        if (index_parameter[i][length[i]] > max_index_parameter) {
          max_index_parameter = index_parameter[i][length[i]];
        }
      }
    }
  }

  return max_index_parameter;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the index parameter frequency distribution.
 */
/*--------------------------------------------------------------*/

void Sequences::build_index_parameter_frequency_distribution()

{
  if (index_parameter) {
    int i , j;


    index_parameter_distribution = new FrequencyDistribution(max_index_parameter_computation() + 1);

    for (i = 0;i < nb_sequence;i++) {
      for (j = 0;j < (index_param_type == POSITION ? length[i] + 1 : length[i]);j++) {
        (index_parameter_distribution->frequency[index_parameter[i][j]])++;
      }
    }

    index_parameter_distribution->offset_computation();
    index_parameter_distribution->nb_element = cumul_length;
    if (index_param_type == POSITION) {
      index_parameter_distribution->nb_element += nb_sequence;
    }
    index_parameter_distribution->max_computation();
    index_parameter_distribution->mean_computation();
    index_parameter_distribution->variance_computation();
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Extraction of the frequency distribution of between-index intervals.
 */
/*--------------------------------------------------------------*/

void Sequences::index_interval_computation()

{
//  if ((index_param_type == TIME) || ((index_param_type == POSITION) &&
//       (type[0] != NB_INTERNODE))) {
  if ((index_param_type == TIME) || (index_param_type == POSITION)) {
    int i , j;


    index_interval = new FrequencyDistribution(max_index_parameter_computation(true) + 1);

    // constitution of the frequency distribution of between-index intervals

    for (i = 0;i < nb_sequence;i++) {
      for (j = 1;j < length[i];j++) {
        (index_interval->frequency[index_parameter[i][j] - index_parameter[i][j - 1]])++;
      }
    }

    index_interval->nb_value_computation();
    index_interval->offset_computation();
    index_interval->nb_element = cumul_length - nb_sequence;
    index_interval->max_computation();
    index_interval->mean_computation();
    index_interval->variance_computation();
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Extraction of the frequency distribution of between-index intervals
 *         of a value of an integer-valued variable.
 *
 *  \param[in] error    reference on a StatError object,
 *  \param[in] variable variable index,
 *  \param[in] value    value.
 *
 *  \return             FrequencyDistribution object.
 */
/*--------------------------------------------------------------*/

FrequencyDistribution* Sequences::value_index_interval_computation(StatError &error , int variable ,
                                                                   int value) const

{
  bool status = true;
  int i , j;
  int previous_index_param , *pindex_param , *pisequence;
  FrequencyDistribution *value_index_interval;


  value_index_interval = NULL;
  error.init();

  if ((index_param_type != TIME) && (index_param_type != POSITION)) {
    status = false;
    error.update(SEQ_error[SEQR_INDEX_PARAMETER_TYPE]);
  }

  if ((variable < 1) || (variable > nb_variable)) {
    status = false;
    error.update(STAT_error[STATR_VARIABLE_INDEX]);
  }

  else {
    variable--;

    if (!marginal_distribution[variable]) {
      status = false;
      ostringstream error_message;
      error_message << STAT_label[STATL_VARIABLE] << " " << variable + 1 << ": "
                    << STAT_error[STATR_MARGINAL_FREQUENCY_DISTRIBUTION];
      error.update((error_message.str()).c_str());
    }

    else if ((value < marginal_distribution[variable]->offset) ||
             (value >= marginal_distribution[variable]->nb_value) ||
             (marginal_distribution[variable]->frequency[value] <= 1)) {
      status = false;
      error.update(SEQ_error[SEQR_VALUE]);
    }
  }

  if (status) {
    value_index_interval = new FrequencyDistribution(max_index_parameter_computation(true) + 1);

    for (i = 0;i < nb_sequence;i++) {
      pindex_param = index_parameter[i];
      pisequence = int_sequence[i][variable];
      previous_index_param = I_DEFAULT;

      for (j = 0;j < length[i];j++) {
        if (*pisequence == value) {
          if (previous_index_param != I_DEFAULT) {
            (value_index_interval->frequency[*pindex_param - previous_index_param])++;
          }
          previous_index_param = *pindex_param;
        }

        pindex_param++;
        pisequence++;
      }
    }

    // computation of the frequency distribution characteristics

    value_index_interval->nb_value_computation();
    value_index_interval->offset_computation();
    value_index_interval->nb_element_computation();

    if (value_index_interval->nb_element == 0) {
      delete value_index_interval;
      value_index_interval = NULL;
      error.update(STAT_error[STATR_EMPTY_SAMPLE]);
    }

    else {
      value_index_interval->max_computation();
      value_index_interval->mean_computation();
      value_index_interval->variance_computation();
    }
  }

  return value_index_interval;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the minimum value taken by a variable.
 *
 *  \param[in] variable variable index.
 */
/*--------------------------------------------------------------*/

void Sequences::min_value_computation(int variable)

{
  int i , j;


  if ((type[variable] != REAL_VALUE) && (type[variable] != AUXILIARY)) {
    min_value[variable] = int_sequence[0][variable][0];

    for (i = 0;i < nb_sequence;i++) {
      for (j = 0;j < length[i];j++) {
        if (int_sequence[i][variable][j] < min_value[variable]) {
          min_value[variable] = int_sequence[i][variable][j];
        }
      }
    }
  }

  else {
    min_value[variable] = real_sequence[0][variable][0];

    for (i = 0;i < nb_sequence;i++) {
      for (j = 0;j < length[i];j++) {
        if (real_sequence[i][variable][j] < min_value[variable]) {
          min_value[variable] = real_sequence[i][variable][j];
        }
      }
    }
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the maximum value taken by a variable.
 *
 *  \param[in] variable variable index.
 */
/*--------------------------------------------------------------*/

void Sequences::max_value_computation(int variable)

{
  int i , j;


  if ((type[variable] != REAL_VALUE) && (type[variable] != AUXILIARY)) {
    max_value[variable] = int_sequence[0][variable][0];

    for (i = 0;i < nb_sequence;i++) {
      for (j = 0;j < length[i];j++) {
        if (int_sequence[i][variable][j] > max_value[variable]) {
          max_value[variable] = int_sequence[i][variable][j];
        }
      }
    }
  }

  else {
    max_value[variable] = real_sequence[0][variable][0];

    for (i = 0;i < nb_sequence;i++) {
      for (j = 0;j < length[i];j++) {
        if (real_sequence[i][variable][j] > max_value[variable]) {
          max_value[variable] = real_sequence[i][variable][j];
        }
      }
    }
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the marginal frequency distribution for
 *         a positive integer-valued variable.
 *
 *  \param[in] variable variable index.
 */
/*--------------------------------------------------------------*/

void Sequences::marginal_frequency_distribution_computation(int variable)

{
  int i , j;


  for (i = 0;i < marginal_distribution[variable]->nb_value;i++) {
    marginal_distribution[variable]->frequency[i] = 0;
  }

  for (i = 0;i < nb_sequence;i++) {
    for (j = 0;j < length[i];j++) {
      (marginal_distribution[variable]->frequency[int_sequence[i][variable][j]])++;
    }
  }

  marginal_distribution[variable]->offset = (int)min_value[variable];
  marginal_distribution[variable]->nb_element_computation();
//  marginal_distribution[variable]->nb_element = cumul_length;
  marginal_distribution[variable]->max_computation();
  marginal_distribution[variable]->mean_computation();
  marginal_distribution[variable]->variance_computation();
}


/*--------------------------------------------------------------*/
/**
 *  \brief Construction of the marginal frequency distribution for
 *         a positive integer-valued variable
 *
 *  \param[in] variable variable index.
 */
/*--------------------------------------------------------------*/

void Sequences::build_marginal_frequency_distribution(int variable)

{
  if (type[variable] != AUXILIARY) {
    if ((type[variable] != REAL_VALUE) && (min_value[variable] >= 0) &&
        (max_value[variable] <= MARGINAL_DISTRIBUTION_MAX_VALUE)) {
      marginal_distribution[variable] = new FrequencyDistribution((int)max_value[variable] + 1);
      marginal_frequency_distribution_computation(variable);
    }

    else {
      build_marginal_histogram(variable);
    }
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Construction of the marginal histogram for a variable.
 *
 *  \param[in] variable   variable index,
 *  \param[in] bin_width  bin width,
 *  \param[in] imin_value minimum value.
 */
/*--------------------------------------------------------------*/

void Sequences::build_marginal_histogram(int variable , double bin_width , double imin_value)

{
  if ((!marginal_histogram[variable]) || (bin_width != marginal_histogram[variable]->bin_width) ||
      (imin_value != D_INF)) {
    int i , j;
    int *pisequence;
    double *prsequence;


    // construction of the histogram

    if (bin_width == D_DEFAULT) {
      bin_width = MAX(::round((max_value[variable] - min_value[variable]) * HISTOGRAM_FREQUENCY / cumul_length) , 1);

#     ifdef MESSAGE
      cout << "\n" << STAT_label[STATL_VARIABLE] << " " << variable + 1 << " - "
           << STAT_label[STATL_BIN_WIDTH] << ": " << bin_width << endl;
//           << " (" << min_value[variable] << ", " << max_value[variable] << ")"
#     endif

    }

    if (imin_value == D_INF) {
      imin_value = floor(min_value[variable] / bin_width) * bin_width;
    }

    if (marginal_histogram[variable]) {
      marginal_histogram[variable]->nb_bin = (int)floor((max_value[variable] - imin_value) / bin_width) + 1;

      delete [] marginal_histogram[variable]->frequency;
      marginal_histogram[variable]->frequency = new int[marginal_histogram[variable]->nb_bin];
    }

    else {
      marginal_histogram[variable] = new Histogram((int)floor((max_value[variable] - imin_value) / bin_width) + 1 , false);

      marginal_histogram[variable]->nb_element = cumul_length;
      marginal_histogram[variable]->type = type[variable];
    }

    marginal_histogram[variable]->bin_width = bin_width;
    marginal_histogram[variable]->min_value = imin_value;
    marginal_histogram[variable]->max_value = ceil(max_value[variable] / bin_width) * bin_width;

    // computation of bin frequencies

    for (i = 0;i < marginal_histogram[variable]->nb_bin;i++) {
      marginal_histogram[variable]->frequency[i] = 0;
    }

    if ((type[variable] != REAL_VALUE) && (type[variable] != AUXILIARY)) {
      for (i = 0;i < nb_sequence;i++) {
        pisequence = int_sequence[i][variable];
        for (j = 0;j < length[i];j++) {
//          (marginal_histogram[variable]->frequency[(int)((*pisequence++ - imin_value) / bin_width)])++;
          (marginal_histogram[variable]->frequency[(int)floor((*pisequence++ - imin_value) / bin_width)])++;
        }
      }
    }

    else {
      for (i = 0;i < nb_sequence;i++) {
        prsequence = real_sequence[i][variable];
        for (j = 0;j < length[i];j++) {
//          (marginal_histogram[variable]->frequency[(int)((*prsequence++ - imin_value) / bin_width)])++;
          (marginal_histogram[variable]->frequency[(int)floor((*prsequence++ - imin_value) / bin_width)])++;
        }
      }
    }

    marginal_histogram[variable]->max_computation();
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Change of the bin width of the marginal histogram.
 *
 *  \param[in] error      reference on a StatError object,
 *  \param[in] variable   variable index,
 *  \param[in] bin_width  bin width,
 *  \param[in] imin_value minimum value.
 *
 *  \return               error status.
 */
/*--------------------------------------------------------------*/

bool Sequences::select_bin_width(StatError &error , int variable ,
                                 double bin_width , double imin_value)

{
  bool status = true;


  error.init();

  if ((variable < 1) || (variable > nb_variable)) {
    status = false;
    error.update(STAT_error[STATR_VARIABLE_INDEX]);
  }

  else {
    variable--;

    if (!marginal_histogram[variable]) {
      status = false;
      error.update(STAT_error[STATR_MARGINAL_HISTOGRAM]);
    }
    if ((bin_width <= 0.) || ((type[variable] != REAL_VALUE) && ((int)bin_width != bin_width))) {
      status = false;
      error.update(STAT_error[STATR_HISTOGRAM_BIN_WIDTH]);
    }
    if ((imin_value != D_INF) && ((imin_value <= min_value[variable] - bin_width) ||
         (imin_value > min_value[variable]) || ((type[variable] != REAL_VALUE) &&
          ((int)imin_value != imin_value)))) {
      status = false;
      error.update(STAT_error[STATR_HISTOGRAM_MIN_VALUE]);
    }
  }

  if (status) {
    build_marginal_histogram(variable , bin_width , imin_value);
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Mean computation for a variable.
 *
 *  \param[in] variable variable index.
 */
/*--------------------------------------------------------------*/

double Sequences::mean_computation(int variable) const

{
  int i , j;
  double mean;


  if (marginal_distribution[variable]) {
    mean = marginal_distribution[variable]->mean;
  }

  else {
    mean = 0.;

    if ((type[variable] != REAL_VALUE) && (type[variable] != AUXILIARY)) {
      for (i = 0;i < nb_sequence;i++) {
        for (j = 0;j < length[i];j++) {
          mean += int_sequence[i][variable][j];
        }
      }
    }

    else {
      for (i = 0;i < nb_sequence;i++) {
        for (j = 0;j < length[i];j++) {
          mean += real_sequence[i][variable][j];
        }
      }
    }

    mean /= cumul_length;
  }

  return mean;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Variance computation for a variable.
 *
 *  \param[in] variable variable index,
 *  \param[in] mean     mean.
 *
 *  \return             variance.
 */
/*--------------------------------------------------------------*/

double Sequences::variance_computation(int variable , double mean) const

{
  int i , j;
  double variance , diff;
  long double square_sum;


  if (marginal_distribution[variable]) {
    variance = marginal_distribution[variable]->variance;
  }

  else {
    if (cumul_length > 1) {
      square_sum = 0.;

      if ((type[variable] != REAL_VALUE) && (type[variable] != AUXILIARY)) {
        for (i = 0;i < nb_sequence;i++) {
          for (j = 0;j < length[i];j++) {
            diff = int_sequence[i][variable][j] - mean;
            square_sum += diff * diff;
          }
        }
      }

      else {
        for (i = 0;i < nb_sequence;i++) {
          for (j = 0;j < length[i];j++) {
            diff = real_sequence[i][variable][j] - mean;
            square_sum += diff * diff;
          }
        }
      }

      variance = square_sum / (cumul_length - 1);
    }

    else {
      variance = 0.;
    }
  }

  return variance;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the mean absolute deviation for a variable.
 *
 *  \param[in] variable variable index,
 *  \param[in] location location measure (e.g. mean or median).
 *
 *  \return             mean absolute deviation.
 */
/*--------------------------------------------------------------*/

double Sequences::mean_absolute_deviation_computation(int variable , double location) const

{
  int i , j;
  double mean_absolute_deviation;


  if (marginal_distribution[variable]) {
    mean_absolute_deviation = marginal_distribution[variable]->mean_absolute_deviation_computation(location);
  }

  else {
    mean_absolute_deviation = 0.;

    if ((type[variable] != REAL_VALUE) && (type[variable] != AUXILIARY)) {
      for (i = 0;i < nb_sequence;i++) {
        for (j = 0;j < length[i];j++) {
          mean_absolute_deviation += fabs(int_sequence[i][variable][j] - location);
        }
      }
    }

    else {
      for (i = 0;i < nb_sequence;i++) {
        for (j = 0;j < length[i];j++) {
          mean_absolute_deviation += fabs(real_sequence[i][variable][j] - location);
        }
      }
    }

    mean_absolute_deviation /= cumul_length;
  }

  return mean_absolute_deviation;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the mean absolute difference for a variable.
 *
 *  \param[in] variable variable index.
 *
 *  \return             mean absolute difference.
 */
/*--------------------------------------------------------------*/

double Sequences::mean_absolute_difference_computation(int variable) const

{
  int i , j , k , m;
  double mean_absolute_difference;


  mean_absolute_difference = 0.;

  if (cumul_length > 1) {
    if ((type[variable] != REAL_VALUE) && (type[variable] != AUXILIARY)) {
      for (i = 0;i < nb_sequence;i++) {
        for (j = 0;j < length[i];j++) {
          for (k = j + 1;k < length[i];k++) {
            mean_absolute_difference += abs(int_sequence[i][variable][j] -
                                            int_sequence[i][variable][k]);
          }
        }

        for (j = i + 1;j < nb_sequence;j++) {
          for (k = 0;k < length[i];k++) {
            for (m = 0;m < length[j];m++) {
              mean_absolute_difference += abs(int_sequence[i][variable][k] -
                                              int_sequence[j][variable][m]);
            }
          }
        }
      }
    }

    else {
      for (i = 0;i < nb_sequence;i++) {
        for (j = 0;j < length[i];j++) {
          for (k = j + 1;k < length[i];k++) {
            mean_absolute_difference += fabs(real_sequence[i][variable][j] -
                                             real_sequence[i][variable][k]);
          }
        }

        for (j = i + 1;j < nb_sequence;j++) {
          for (k = 0;k < length[i];k++) {
            for (m = 0;m < length[j];m++) {
              mean_absolute_difference += fabs(real_sequence[i][variable][k] -
                                               real_sequence[j][variable][m]);
            }
          }
        }
      }
    }

    mean_absolute_difference = 2 * mean_absolute_difference /
                               (cumul_length * (double)(cumul_length - 1));
  }

  return mean_absolute_difference;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the coefficient of skewness for a variable.
 *
 *  \param[in] variable variable index,
 *  \param[in] mean     mean,
 *  \param[in] variance variance.
 *
 *  \return             coefficient of skewness.
 */
/*--------------------------------------------------------------*/

double Sequences::skewness_computation(int variable , double mean , double variance) const

{
  int i , j;
  double skewness , diff;
  long double cube_sum;


  if (marginal_distribution[variable]) {
    skewness = marginal_distribution[variable]->skewness_computation();
  }

  else {
    if ((cumul_length > 2) && (variance > 0.)) {
      cube_sum = 0.;

      if ((type[variable] != REAL_VALUE) && (type[variable] != AUXILIARY)) {
        for (i = 0;i < nb_sequence;i++) {
          for (j = 0;j < length[i];j++) {
            diff = int_sequence[i][variable][j] - mean;
            cube_sum += diff * diff * diff;
          }
        }
      }

      else {
        for (i = 0;i < nb_sequence;i++) {
          for (j = 0;j < length[i];j++) {
            diff = real_sequence[i][variable][j] - mean;
            cube_sum += diff * diff * diff;
          }
        }
      }

      skewness = cube_sum * cumul_length / ((cumul_length - 1) *
                  (double)(cumul_length - 2) * pow(variance , 1.5));
    }

    else {
      skewness = 0.;
    }
  }

  return skewness;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the excess kurtosis for a variable:
 *         excess kurtosis = coefficient of kurtosis - 3.
 *
 *  \param[in] variable variable index,
 *  \param[in] mean     mean,
 *  \param[in] variance variance.
 *
 *  \return             excess kurtosis.
 */
/*--------------------------------------------------------------*/

double Sequences::kurtosis_computation(int variable , double mean , double variance) const

{
  int i , j;
  double kurtosis , diff;
  long double power_sum;


  if (marginal_distribution[variable]) {
    kurtosis = marginal_distribution[variable]->kurtosis_computation();
  }

  else {
    if (variance > 0.) {
      power_sum = 0.;

      if ((type[variable] != REAL_VALUE) && (type[variable] != AUXILIARY)) {
        for (i = 0;i < nb_sequence;i++) {
          for (j = 0;j < length[i];j++) {
            diff = int_sequence[i][variable][j] - mean;
            power_sum += diff * diff * diff * diff;
          }
        }
      }

      else {
        for (i = 0;i < nb_sequence;i++) {
          for (j = 0;j < length[i];j++) {
            diff = real_sequence[i][variable][j] - mean;
            power_sum += diff * diff * diff * diff;
          }
        }
      }

      kurtosis = power_sum / ((cumul_length - 1) * variance * variance) - 3.;
    }

    else {
      kurtosis = -2.;
    }
  }

  return kurtosis;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the mean direction for a circular variable.
 *
 *  \param[in] variable variable index,
 *  \param[in] unit     unit (DEGREE/RADIAN).
 *
 *  \return             mean direction.
 */
/*--------------------------------------------------------------*/

double* Sequences::mean_direction_computation(int variable , angle_unit unit) const

{
  int i , j;
  double *mean_direction;


  mean_direction = new double[4];
//  mean_direction = new double[2];

  mean_direction[0] = 0.;
  mean_direction[1] = 0.;

  switch (type[variable]) {

  case INT_VALUE : {
    for (i = 0;i < nb_sequence;i++) {
      for (j = 0;j < length[i];j++) {
        mean_direction[0] += cos(int_sequence[i][variable][j] * M_PI / 180);
        mean_direction[1] += sin(int_sequence[i][variable][j] * M_PI / 180);
      }
    }
    break;
  }

  case REAL_VALUE : {
    switch (unit) {

    case DEGREE : {
      for (i = 0;i < nb_sequence;i++) {
        for (j = 0;j < length[i];j++) {
          mean_direction[0] += cos(real_sequence[i][variable][j] * M_PI / 180);
          mean_direction[1] += sin(real_sequence[i][variable][j] * M_PI / 180);
        }
      }
      break;
    }

    case RADIAN : {
      for (i = 0;i < nb_sequence;i++) {
        for (j = 0;j < length[i];j++) {
          mean_direction[0] += cos(real_sequence[i][variable][j]);
          mean_direction[1] += sin(real_sequence[i][variable][j]);
        }
      }
      break;
    }
    }
    break;
  }
  }

  mean_direction[0] /= cumul_length;
  mean_direction[1] /= cumul_length;

  mean_direction[2] = sqrt(mean_direction[0] * mean_direction[0] +
                           mean_direction[1] * mean_direction[1]);

  if (mean_direction[2] > 0.) {
    mean_direction[3] = atan(mean_direction[1] / mean_direction[0]);

    if (mean_direction[0] < 0.) {
      mean_direction[3] += M_PI;
    }
    if (unit == DEGREE) {
      mean_direction[3] *= (180 / M_PI);
    }
  }

  else {
    mean_direction[3] = D_DEFAULT;
  }

  return mean_direction;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the root mean square error or the mean absolute error for a variable.
 *
 *  \param[in] error       reference on a StatError object,
 *  \param[in] display     flag for displaying the result,
 *  \param[in] variable    variable index,
 *  \param[in] iidentifier sequence identifier,
 *  \param[in] robust      flag computation of mean absolute error.
 *
 *  \return                root mean square error or mean absolute error.
 */
/*--------------------------------------------------------------*/

bool Sequences::mean_error_computation(StatError &error , bool display , int variable ,
                                       int iidentifier , bool robust) const

{
  bool status = true;
  int i , j;
  int index;
  double mean_error , diff;
  long double mean_squared_error;


  error.init();

  if ((variable < 1) || (variable >= nb_variable)) {
    status = false;
    error.update(STAT_error[STATR_VARIABLE_INDEX]);
  }

  else {
    variable--;

    if ((type[variable] != INT_VALUE) && (type[variable] != REAL_VALUE)) {
      status = false;
      ostringstream correction_message;
      correction_message << STAT_variable_word[INT_VALUE] << " or " << STAT_variable_word[REAL_VALUE];
      error.correction_update(STAT_error[STATR_VARIABLE_TYPE] , (correction_message.str()).c_str());
    }

    if (type[variable + 1] != AUXILIARY) {
      status = false;
      ostringstream error_message;
      error_message << STAT_label[STATL_VARIABLE] << " " << variable + 1 << ": "
                    << STAT_error[STATR_VARIABLE_TYPE];
      error.correction_update((error_message.str()).c_str() , STAT_variable_word[AUXILIARY]);
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

  else {
    index = I_DEFAULT;
  }

  if (status) {
    if (robust) {
      mean_error = 0.;

      switch (type[variable]) {

      case INT_VALUE : {
        for (i = 0;i < nb_sequence;i++) {
          if ((index == I_DEFAULT) || (index == i)) {
            for (j = 0;j < length[i];j++) {
              mean_error += fabs(int_sequence[i][variable][j] - real_sequence[i][variable + 1][j]);
            }
          }
        }
        break;
      }

      case REAL_VALUE : {
        for (i = 0;i < nb_sequence;i++) {
          if ((index == I_DEFAULT) || (index == i)) {
            for (j = 0;j < length[i];j++) {
              mean_error += fabs(real_sequence[i][variable][j] - real_sequence[i][variable + 1][j]);
            }
          }
        }
        break;
      }
      }

      if (index == I_DEFAULT) {
        mean_error /= cumul_length;
      }
      else {
        mean_error /= length[index];
      }
    }

    else {
      mean_squared_error = 0.;

      switch (type[variable]) {

      case INT_VALUE : {
        for (i = 0;i < nb_sequence;i++) {
          if ((index == I_DEFAULT) || (index == i)) {
            for (j = 0;j < length[i];j++) {
              diff = int_sequence[i][variable][j] - real_sequence[i][variable + 1][j];
              mean_squared_error += diff * diff;
            }
          }
        }
        break;
      }

      case REAL_VALUE : {
        for (i = 0;i < nb_sequence;i++) {
          if ((index == I_DEFAULT) || (index == i)) {
            for (j = 0;j < length[i];j++) {
              diff = real_sequence[i][variable][j] - real_sequence[i][variable + 1][j];
              mean_squared_error += diff * diff;
            }
          }
        }
        break;
      }
      }

      if (index == I_DEFAULT) {
        mean_error = sqrtl(mean_squared_error / cumul_length);
      }
      else {
        mean_error = sqrtl(mean_squared_error / length[index]);
      }
    }

    if (display) {
      cout << "\n";
      if (((type[0] != STATE) && (nb_variable > 2)) || ((type[0] == STATE) && (nb_variable > 3))) {
        cout << STAT_label[STATL_VARIABLE] << " " << variable + 1 << "   ";
      }

      if (robust) {
        cout << SEQ_label[SEQL_MEAN_ABSOLUTE_ERROR] << ": " << mean_error << endl;
      }
      else {
        cout << SEQ_label[SEQL_ROOT_MEAN_SQUARE_ERROR] << ": " << mean_error << endl;
      }
    }
  }

  return status;
}


};  // namespace sequence_analysis
