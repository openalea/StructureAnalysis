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
#include "tool/rw_tokenizer.h"
#include "tool/rw_cstring.h"
#include "tool/rw_locale.h"
// #include <rw/vstream.h>
// #include <rw/rwfile.h>
#include "stat_tool/stat_tools.h"
#include "stat_tool/curves.h"
#include "stat_tool/markovian.h"
#include "stat_tool/stat_label.h"
#include "sequences.h"
#include "sequence_label.h"
#include "tool/config.h"

using namespace std;

extern int column_width(int min_value , int max_value);
extern char* label(const char *file_name);



/*--------------------------------------------------------------*
 *
 *  Construction d'un objet Sequences a partir d'un fichier.
 *
 *  arguments : reference sur un objet Format_error, path.
 *
 *--------------------------------------------------------------*/

Sequences* sequences_ascii_read(Format_error &error , const char *path)

{
  RWLocaleSnapshot locale("en");
  RWCString buffer , token;
  size_t position;
  bool status , lstatus;
  register int i , j , k , m;
  int line , read_line , initial_nb_line , max_length , nb_variable = 0 , nb_sequence ,
      index , line_continue , type[SEQUENCE_NB_VARIABLE] , *length;
  long value;
  Sequences *seq;
  ifstream in_file(path);


  seq = 0;
  error.init();

  if (!in_file) {
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {

    // 1ere passe : analyse de la ligne definissant le nombre de variables

    status = true;
    length = 0;
    line = 0;
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

        // test nombre de variables

        case 0 : {
          lstatus = locale.stringToNum(token , &value);
          if (lstatus) {
            if ((value < 1) || (value > SEQUENCE_NB_VARIABLE)) {
              lstatus = false;
            }
            else {
              nb_variable = value;
            }
          }

          if (!lstatus) {
            status = false;
            error.update(STAT_parsing[STATP_NB_VARIABLE] , line , i + 1);
          }
          break;
        }

        // test mot cle VARIABLE(S)

        case 1 : {
          if (token != STAT_word[nb_variable == 1 ? STATW_VARIABLE : STATW_VARIABLES]) {
            status = false;
            error.correction_update(STAT_parsing[STATP_KEY_WORD] ,
                                    STAT_word[nb_variable == 1 ? STATW_VARIABLE : STATW_VARIABLES] , line , i + 1);
          }
          break;
        }
        }

        i++;
      }

      if (i > 0) {
        if (i != 2) {
          status = false;
          error.update(STAT_parsing[STATP_FORMAT] , line);
        }

        read_line++;
        break;
      }
    }

    if (read_line < 1) {
      status = false;
      error.update(STAT_parsing[STATP_FORMAT]);
    }

    // analyse des lignes definissant le type de chaque variable

    if (status) {
      for (i = 0;i < nb_variable;i++) {
        type[i] = I_DEFAULT;
      }

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
            lstatus = locale.stringToNum(token , &value);
            if ((lstatus) && (value != read_line + 1)) {
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
            for (j = 0;j < SEQUENCE_NB_TYPE;j++) {
              if (token == STAT_sequence_word[j]) {
                if ((((j == TIME) || (j == TIME_INTERVAL) || (j == POSITION) || (j == POSITION_INTERVAL)) &&
                     (read_line != 0)) || ((j == NB_INTERNODE) && ((read_line != 1) ||
                      ((read_line == 1) && (type[0] != POSITION) && (type[0] != POSITION_INTERVAL))))) {
                  status = false;
                  error.update(SEQ_parsing[SEQP_VARIABLE_TYPE] , line , i + 1);
                }

                else {
                  if (j == STATE) {
                    j = INT_VALUE;
                  }

                  type[read_line] = j;
                }
                break;
              }
            }

            if (j == SEQUENCE_NB_TYPE) {
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
        if ((((type[0] == TIME) || (type[0] == TIME_INTERVAL)) && (read_line < 2)) ||
            (((type[0] == POSITION) || (type[0] == POSITION_INTERVAL)) &&
             ((read_line < 2) || ((read_line > 2) && (type[1] == NB_INTERNODE))))) {
          status = false;
          error.update(SEQ_parsing[SEQP_VARIABLE_TYPE]);
        }
      }

      initial_nb_line = line;
    }

    if (status) {
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

      length = new int[nb_sequence];
      for (i = 0;i < nb_sequence;i++) {
        if (nb_variable == 1) {
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

      switch (type[0]) {
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

          if ((nb_variable > 1) && (j % (nb_variable + 1) == nb_variable)) {
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
            if ((nb_variable == 1) && (token == "\\")) {
              line_continue = true;
            }

            else {
              lstatus = locale.stringToNum(token , &value);
              if ((lstatus) && ((((type[k] == TIME) || (type[k] == POSITION) ||
                     (type[k] == POSITION_INTERVAL) || (type[k] == NB_INTERNODE)) && (value < 0)) ||
                   ((type[k] == TIME_INTERVAL) && (value <= 0)))) {
                lstatus = false;
              }

              if (!lstatus) {
                status = false;
                error.update(STAT_parsing[STATP_DATA_TYPE] , line , j + 1);
              }

              else {
                switch (type[k]) {

                case TIME : {
                  if (value <= index) {
                    status = false;
                    error.update(SEQ_parsing[SEQP_TIME_INDEX_ORDER] , line , j + 1);
                  }
                  else {
                    index = value;
                  }
                  break;
                }

                case POSITION : {
                  if (value < index) {
                    status = false;
                    error.update(SEQ_parsing[SEQP_POSITION_ORDER] , line , j + 1);
                  }
                  else {
                    index = value;
                  }
                  break;
                }
                }
              }

              if (nb_variable == 1) {
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
          if (nb_variable > 1) {
            if (((line_continue) || ((type[0] != POSITION) && (type[0] != POSITION_INTERVAL))) &&
                (k != nb_variable)) {
              status = false;
              error.update(STAT_parsing[STATP_FORMAT] , line , j);
            }
          }

          if (!line_continue) {
            if ((type[0] == POSITION) || (type[0] == POSITION_INTERVAL)) {
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

            switch (type[0]) {
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

      seq = new Sequences(nb_variable , type , nb_sequence , 0 , length , false);

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
          if ((nb_variable > 1) && (m % (nb_variable + 1) == nb_variable)) {
            if (token == "\\") {
              line_continue = true;
            }
            k = 0;
            j++;
          }

          else {
            if ((nb_variable == 1) && (token == "\\")) {
              line_continue = true;
            }

            else {
              locale.stringToNum(token , &value);
              seq->sequence[i][k][j] = value;

              if (nb_variable == 1) {
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

      if ((seq->type[0] == TIME_INTERVAL) || (seq->type[0] == POSITION_INTERVAL)) {
        seq->index_computation();
      }
      if ((seq->type[0] == TIME) || ((seq->type[0] == POSITION) && (seq->type[1] == INT_VALUE))) {
        seq->index_interval_computation();
      }

      for (i = 0;i < nb_variable;i++) {
        if (seq->type[i] != POSITION) {
          seq->min_value_computation(i);
          seq->max_value_computation(i);
          seq->build_marginal_histogram(i);
        }
      }
    }

    delete [] length;
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Construction d'un objet Sequences a partir d'un fichier.
 *
 *  arguments : reference sur un objet Format_error, path.
 *
 *--------------------------------------------------------------*/

/* Sequences* old_sequences_ascii_read(Format_error &error , const char *path)

{
  RWCString buffer , token;
  RWCRegexp reg("[0-9]+");
  size_t position , start , extent;
  char number[2];
  bool status;
  register int i , j , k;
  int line , nb_sequence , max_length , type[1] , *length;
  Sequences *seq;
  ifstream in_file(path);


  seq = 0;
  error.init();

  if (!in_file) {
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {

    // 1ere passe : recherche du nombre de sequences

    status = true;
    length = 0;
    type[0] = INT_VALUE;
    line = 0;
    nb_sequence = 0;

    while (buffer.readLine(in_file , false)) {
      line++;

#     ifdef DEBUG
      cout << line << "  " << buffer << endl;
#     endif

      position = buffer.first('#');
      if (position != RW_NPOS) {
        buffer.remove(position);
      }
      i = 0;

      RWCTokenizer next(buffer);

      while (!((token = next()).isNull())) {
        i++;
      }

      if (i > 0) {
        if (i != 1) {
          status = false;
          error.update(STAT_parsing[STATP_FORMAT] , line);
        }
        nb_sequence++;
      }
    }

    if (nb_sequence == 0) {
      status = false;
      error.update(STAT_parsing[STATP_FORMAT]);
    }

#   ifdef DEBUG
//    cout << "\nnumber of sequences : " << nb_sequence << endl;
#   endif

    // 2eme passe : analyse du format des sequences

    if (status) {
      in_file.close();
      in_file.open(path , ios::in);

      length = new int[nb_sequence];

      max_length = 0;
      line = 0;
      i = 0;

      while (buffer.readLine(in_file , false)) {
        line++;

        position = buffer.first('#');
        if (position != RW_NPOS) {
          buffer.remove(position);
        }
        j = 0;

        RWCTokenizer next(buffer);

        while (!((token = next()).isNull())) {
          if (j == 0) {
            length[i] = token.length();

            if (length[i] > max_length) {
              max_length = length[i];
            }

            start = token.index(reg , &extent);
            if ((start != 0) || (extent != length[i])) {
              status = false;
              error.update(STAT_parsing[STATP_DATA_TYPE] , line);
            }
          }
          j++;
        }

        if (j > 0) {
          i++;
        }
      }

      if (max_length <= 1) {
        status = false;
        error.update(SEQ_parsing[SEQP_MAX_SEQUENCE_LENGTH]);
      }

    }

    // 3eme passe : copie des sequences

    if (status) {
      in_file.close();
      in_file.open(path , ios::in);

      seq = new Sequences(1 , type , nb_sequence , 0 , length , false);
      i = 0;

      while (buffer.readLine(in_file , false)) {
        position = buffer.first('#');
        if (position != RW_NPOS) {
          buffer.remove(position);
        }
        j = 0;

        RWCTokenizer next(buffer);

        while (!((token = next()).isNull())) {
          if (j == 0) {
            for (k = 0;k < length[i];k++) {
              number[0] = token[k];
              istrstream(number) >> seq->sequence[i][0][k];
            }
          }
          j++;
        }

        if (j > 0) {
          i++;
        }
      }

      seq->min_value_computation(0);
      seq->max_value_computation(0);
      seq->build_marginal_histogram(0);
    }

    delete [] length;
  }

  return seq;
} */


/*--------------------------------------------------------------*
 *
 *  Ecriture sur une ligne d'un objet Sequences.
 *
 *  argument : stream.
 *
 *--------------------------------------------------------------*/

ostream& Sequences::line_write(ostream &os) const

{
  os << nb_variable << " " << STAT_word[nb_variable == 1 ? STATW_VARIABLE : STATW_VARIABLES] << "   "
     << nb_sequence << " " << SEQ_label[nb_sequence == 1 ? SEQL_SEQUENCE : SEQL_SEQUENCES] << "   "
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


  os << nb_variable << " " << STAT_word[nb_variable == 1 ? STATW_VARIABLE : STATW_VARIABLES] << endl;

  for (i = 0;i < nb_variable;i++) {
    os << "\n" << STAT_word[STATW_VARIABLE] << " " << i + 1 << " : "
       << STAT_sequence_word[type[i]];

    if (type[i] != POSITION) {
      os << "   ";
      if (comment_flag) {
        os << "# ";
      }
      os << "(" << STAT_label[STATL_MIN_VALUE] << ": " << min_value[i] << ", "
         << STAT_label[STATL_MAX_VALUE] << ": " << max_value[i] << ")" << endl;

      os << "\n";
      if (comment_flag) {
        os << "# ";
      }
      os << (type[i] == TIME ? SEQ_label[SEQL_TIME] : STAT_label[STATL_VALUE]) << " "
         << STAT_label[STATL_HISTOGRAM] << " - ";

      if (marginal[i]) {
        marginal[i]->ascii_characteristic_print(os , exhaustive , comment_flag);

        if ((marginal[i]->nb_value <= ASCII_NB_VALUE) || (exhaustive)) {
          os << "\n";
          if (comment_flag) {
            os << "# ";
          }
          os << "   | " << (type[i] == TIME ? SEQ_label[SEQL_TIME] : STAT_label[STATL_VALUE]) << " "
             << STAT_label[STATL_HISTOGRAM] << endl;
          marginal[i]->ascii_print(os , comment_flag);
        }
      }

      else if (type[i] == INT_VALUE) {
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
  }

  if (index_interval) {
    os << "\n";
    if (comment_flag) {
      os << "# ";
    }
    os << (type[0] == TIME ? SEQ_label[SEQL_TIME_INTERVAL] : SEQ_label[SEQL_POSITION_INTERVAL])
       << " " << STAT_label[STATL_HISTOGRAM] << " - ";
    index_interval->ascii_characteristic_print(os , false , comment_flag);

    if (exhaustive) {
      os << "\n";
      if (comment_flag) {
        os << "# ";
      }
      os << "   | " << (type[0] == TIME ? SEQ_label[SEQL_TIME_INTERVAL] : SEQ_label[SEQL_POSITION_INTERVAL])
         << " " << STAT_label[STATL_HISTOGRAM] << endl;
      index_interval->ascii_print(os , comment_flag);
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
 *              flag commentaire, nombre de caracteres par ligne.
 *
 *--------------------------------------------------------------*/

ostream& Sequences::ascii_print(ostream &os , char format , bool comment_flag ,
                                int line_nb_character) const

{
  register int i , j , k , m;


  switch (format) {

  case 'c' : {
    for (i = 0;i < nb_sequence;i++) {
      os << "\n";

#ifdef DEBUG
      for (j = 0;j < length[i];j++) {
        for (k = 0;k < nb_variable;k++) {
          os << sequence[i][k][j] << " ";
        }

        if (j < length[i] - 1) {
          if ((os.rdbuf())->out_waiting() > line_nb_character) {
            os << "\\" << endl;
          }

          else {
            if (nb_variable > 1) {
              os << "| ";
            }
          }
        }
      }
#else
      std::ostringstream sos;

      for (j = 0;j < length[i];j++) {
        for (k = 0;k < nb_variable;k++) {
          sos << sequence[i][k][j] << " ";
        }

        if (j < length[i] - 1) {
          if (sos.str().size() > line_nb_character) {
            os << sos.str() << "\\" << endl;
            sos.str("");
          }

          else {
            if (nb_variable > 1) {
              sos << "| ";
            }
          }
        }
      }
      os << sos.str();
#endif

      if (type[0] == POSITION) {
        os << "| " << sequence[i][0][length[i]];
      }

      os << "   ";
      if (comment_flag) {
        os << "# ";
      }
      os << "(" << identifier[i] << ")" << endl;
    }
    break;
  }

  case 'l' : {
    int buff , start , width;
    long old_adjust;


    old_adjust = os.setf(ios::right , ios::adjustfield);

    width = 0;
    for (i = 0;i < nb_variable;i++) {
      buff = column_width(min_value[i] , max_value[i]);
      if (buff > width) {
        width = buff;
      }
    }

    for (i = 0;i < nb_sequence;i++) {
      os << "\n";
      start = 0;

      for (j = 0;j < length[i];j++) {
        os << setw(j == start ? width : width + 1) << sequence[i][0][j];

        if (j < length[i] - 1) {
          if ((j - start) * (width + 1) > 1000) {
//          if ((j - start) * (width + 1) > line_nb_character) {
            os << " \\" << endl;

            for (k = 1;k < nb_variable;k++) {
              os << setw(width) << sequence[i][k][start];
              for (m = start + 1;m <= j;m++) {
                os << setw(width + 1) << sequence[i][k][m];
              }
              os << " \\" << endl;
            }
            start = j + 1;
          }
        }

        else {
          if (type[0] == POSITION) {
            os << setw(width + 1) << sequence[i][0][length[i]];
          }

          for (k = 1;k < nb_variable;k++) {
            os << endl;
            os << setw(width) << sequence[i][k][start];
            for (m = start + 1;m <= j;m++) {
              os << setw(width + 1) << sequence[i][k][m];
            }
          }
        }
      }

      os << "   ";
      if (comment_flag) {
        os << "# ";
      }
      os << "(" << identifier[i] << ")" << endl;
    }

    os.setf((FMTFLAGS)old_adjust , ios::adjustfield);
    break;
  }

  case 'a' : {
    os << "[";
    for (i = 0;i < nb_sequence;i++) {

#ifdef DEBUG
      os << "[";
      for (j = 0;j < length[i];j++) {
        if (nb_variable == 1) {
          os << sequence[i][0][j];
        }

        else {
          os << "[";
          for (k = 0;k < nb_variable - 1;k++) {
            os << sequence[i][k][j] << ",";
          }
          os << sequence[i][nb_variable - 1][j] << "]";
        }

        if (j < length[i] - 1) {
          os << ",";
          if ((os.rdbuf())->out_waiting() > line_nb_character) {
            os << "\\" << endl;
            os << "  ";
          }
        }
      }
#else
      std::ostringstream sos;

      sos << "[";
      for (j = 0;j < length[i];j++) {
        if (nb_variable == 1) {
          sos << sequence[i][0][j];
        }

        else {
          sos << "[";
          for (k = 0;k < nb_variable - 1;k++) {
            sos << sequence[i][k][j] << ",";
          }
          sos << sequence[i][nb_variable - 1][j] << "]";
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
#endif

      if (type[0] == POSITION) {
        os << ",[" << sequence[i][0][length[i]];
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

    out_file << nb_variable << "\t" << STAT_word[nb_variable == 1 ? STATW_VARIABLE : STATW_VARIABLES] << endl;

    for (i = 0;i < nb_variable;i++) {
      out_file << "\n" << STAT_word[STATW_VARIABLE] << "\t" << i + 1 << "\t"
               << STAT_sequence_word[type[i]];

      if (type[i] != POSITION) {
        out_file << "\t\t" << STAT_label[STATL_MIN_VALUE] << "\t" << min_value[i]
                 << "\t\t" << STAT_label[STATL_MAX_VALUE] << "\t" << max_value[i] << endl;

        out_file << "\n" << (type[i] == TIME ? SEQ_label[SEQL_TIME] : STAT_label[STATL_VALUE]) << " "
                 << STAT_label[STATL_HISTOGRAM] << "\t";

        if (marginal[i]) {
          marginal[i]->spreadsheet_characteristic_print(out_file);

          out_file << "\n\t" << (type[i] == TIME ? SEQ_label[SEQL_TIME] : STAT_label[STATL_VALUE]) << " "
                   << STAT_label[STATL_HISTOGRAM] << endl;
          marginal[i]->spreadsheet_print(out_file);
        }

        else if (type[i] == INT_VALUE) {
          out_file << STAT_label[STATL_SAMPLE_SIZE] << "\t" << cumul_length << endl;

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
    }

    if (index_interval) {
      out_file << "\n" << (type[0] == TIME ? SEQ_label[SEQL_TIME_INTERVAL] : SEQ_label[SEQL_POSITION_INTERVAL])
               << " " << STAT_label[STATL_HISTOGRAM] << "\t";
      index_interval->spreadsheet_characteristic_print(out_file);

      out_file << "\n\t" << (type[0] == TIME ? SEQ_label[SEQL_TIME_INTERVAL] : SEQ_label[SEQL_POSITION_INTERVAL])
               << " " << STAT_label[STATL_HISTOGRAM] << endl;
      index_interval->spreadsheet_print(out_file);
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

  phisto = new const Histogram*[nb_variable + 1];

  nb_histo = 0;
  for (i = 0;i < nb_variable;i++) {
    if ((type[i] != POSITION) && (marginal[i])) {
      phisto[nb_histo++] = marginal[i];
    }
  }

  if (index_interval) {
    phisto[nb_histo++] = index_interval;
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

      for (k = 0;k < nb_variable;k++) {
        if ((type[k] != POSITION) && (marginal[k])) {
          if (marginal[k]->nb_value - 1 < TIC_THRESHOLD) {
            out_file << "set xtics 0,1" << endl;
          }
          if ((int)(marginal[k]->max * YSCALE) + 1 < TIC_THRESHOLD) {
            out_file << "set ytics 0,1" << endl;
          }

          out_file << "plot [0:" << MAX(marginal[k]->nb_value - 1 , 1) << "] [0:"
                   << (int)(marginal[k]->max * YSCALE) + 1 << "] \""
                   << label((data_file_name.str()).c_str()) << "\" using " << j++ << " title \""
                   << STAT_label[STATL_VARIABLE] << " " << k + 1 << " - "
                   << (type[k] == TIME ? SEQ_label[SEQL_TIME] : STAT_label[STATL_VALUE]) << " "
                   << STAT_label[STATL_HISTOGRAM] << "\" with impulses" << endl;

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
                 << (type[0] == TIME ? SEQ_label[SEQL_TIME_INTERVAL] : SEQ_label[SEQL_POSITION_INTERVAL])
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
        for (k = 0;k < nb_variable;k++) {
          out_file << sequence[index[j]][k][i] << " ";
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
  int min_index , max_index , *pfrequency , *length_nb_sequence;
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

      if ((type[0] == TIME) || (type[0] == POSITION)) {
        max_index = 0;
        for (i = 0;i < nb_sequence;i++) {
          if (sequence[i][0][length[i] - 1] > max_index) {
            max_index = sequence[i][0][length[i] - 1];
          }
        }

        min_index = max_index;
        for (i = 0;i < nb_sequence;i++) {
          if (sequence[i][0][0] < min_index) {
            min_index = sequence[i][0][0];
          }
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

        out_file << "set border 15 lw 0\n" << "set tics out\n" << "set xtics nomirror\n";

        for (j = ((type[0] == TIME) || (type[0] == POSITION) ? 1 : 0);j < nb_variable;j++) {
          for (k = 0;k < hlength->nb_value;k++) {
            length_nb_sequence[k] = 0;
          }

          out_file << "set title \"";
          if (title) {
            out_file << title;
            if (nb_variable > ((type[0] == TIME) || (type[0] == POSITION) ? 2 : 1)) {
              out_file << " - ";
            }
          }

          if ((type[0] == TIME) || (type[0] == POSITION)) {
            if (nb_variable > 2) {
              out_file << STAT_label[STATL_VARIABLE] << " " << j;
            }
            out_file << "\"\n\n";

            if (max_index - min_index < TIC_THRESHOLD) {
              out_file << "set xtics 0,1" << endl;
            }
            if (max_value[j] - min_value[j] < TIC_THRESHOLD) {
              out_file << "set ytics " << MIN(min_value[j] , 0) << ",1" << endl;
            }

            out_file << "plot [" << min_index << ":" << max_index << "] ["
                     << MIN(min_value[j] , 0) << ":" << MAX(max_value[j] , min_value[j] + 1) << "] ";
            for (k = 0;k < nb_sequence;k++) {
              out_file << "\"" << label((data_file_name[length[k]].str()).c_str()) << "\" using "
                       << length_nb_sequence[length[k]] * nb_variable + 1 << " : "
                       << length_nb_sequence[length[k]] * nb_variable + j + 1;
              if (nb_sequence <= PLOT_TITLE_NB_SEQUENCE) {
                out_file << " title \"" << identifier[k] << "\" with linespoints";
              }
              else {
                out_file << " notitle with linespoints";
              }
              if (k < nb_sequence - 1) {
                out_file << ",\\";
              }
              out_file << endl;
              length_nb_sequence[length[k]]++;
            }

            if (max_index - min_index < TIC_THRESHOLD) {
              out_file << "set xtics autofreq" << endl;
            }
            if (max_value[j] - min_value[j] < TIC_THRESHOLD) {
              out_file << "set ytics autofreq" << endl;
            }
          }

          else {
            if (nb_variable > 1) {
              out_file << STAT_label[STATL_VARIABLE] << " " << j + 1;
            }
            out_file << "\"\n\n";

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
              if (nb_sequence <= PLOT_TITLE_NB_SEQUENCE) {
                out_file << " title \"" << identifier[k] << "\" with linespoints";
              }
              else {
                out_file << " notitle with linespoints";
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

          if ((i == 0) && (j < nb_variable - 1)) {
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
 *  Fonctions pour la persistance.
 *
 *--------------------------------------------------------------*/

/* RWDEFINE_COLLECTABLE(Sequences , STATI_SEQUENCES);


RWspace Sequences::binaryStoreSize() const

{
  register int i , j;
  RWspace size;


  size = sizeof(nb_variable) + sizeof(type) * nb_variable + sizeof(min_value) * nb_variable +
         sizeof(max_value) * nb_variable;
  for (i = 0;i < nb_variable;i++) {
    size += sizeof(true);
    if (marginal[i]) {
      size += marginal[i]->binaryStoreSize();
    }
  }

  size += sizeof(nb_sequence) + sizeof(*identifier) * nb_sequence + sizeof(max_length) +
          sizeof(cumul_length) + sizeof(*length) * nb_sequence + hlength->binaryStoreSize();

  size += sizeof(true);
  if (index_interval) {
    size += index_interval->binaryStoreSize();
  }

  for (i = 0;i < nb_sequence;i++) {
    for (j = 0;j < nb_variable;j++) {
      size += sizeof(***sequence) * (type[j] == POSITION ? length[i] + 1 : length[i]);
    }
  }

  return size;
}


void Sequences::restoreGuts(RWvistream &is)

{
  bool status;
  register int i , j , k;
  int blength , *psequence;


  remove();

  is >> nb_variable;

  type = new int[nb_variable];
  for (i = 0;i < nb_variable;i++) {
    is >> type[i];
  }

  min_value = new int[nb_variable];
  for (i = 0;i < nb_variable;i++) {
    is >> min_value[i];
  }
  max_value = new int[nb_variable];
  for (i = 0;i < nb_variable;i++) {
    is >> max_value[i];
  }

  marginal = new Histogram*[nb_variable];
  for (i = 0;i < nb_variable;i++) {
    is >> status;
    if (status) {
      marginal[i] = new Histogram();
      marginal[i]->restoreGuts(is);
    }
    else {
      marginal[i] = 0;
    }
  }

  is >> nb_sequence;

  identifier = new int[nb_sequence];
  for (i = 0;i < nb_sequence;i++) {
    is >> identifier[i];
  }

  is >> max_length >> cumul_length;

  length = new int[nb_sequence];
  for (i = 0;i < nb_sequence;i++) {
    is >> length[i];
  }

  hlength = new Histogram();
  hlength->restoreGuts(is);

  is >> status;
  if (status) {
    index_interval = new Histogram();
    index_interval->restoreGuts(is);
  }
  else {
    index_interval = 0;
  }

  sequence = new int**[nb_sequence];
  for (i = 0;i < nb_sequence;i++) {
    sequence[i] = new int*[nb_variable];
    for (j = 0;j < nb_variable;j++) {
      blength = (type[j] == POSITION ? length[i] + 1 : length[i]);
      sequence[i][j] = new int[blength];
      psequence = sequence[i][j];
      for (k = 0;k < blength;k++) {
        is >> *psequence++;
      }
    }
  }
}


void Sequences::restoreGuts(RWFile &file)

{
  bool status;
  register int i , j;
  int blength;


  remove();

  file.Read(nb_variable);

  type = new int[nb_variable];
  file.Read(type , nb_variable);

  min_value = new int[nb_variable];
  file.Read(min_value , nb_variable);
  max_value = new int[nb_variable];
  file.Read(max_value , nb_variable);

  marginal = new Histogram*[nb_variable];
  for (i = 0;i < nb_variable;i++) {
    file.Read(status);
    if (status) {
      marginal[i] = new Histogram();
      marginal[i]->restoreGuts(file);
    }
    else {
      marginal[i] = 0;
    }
  }

  file.Read(nb_sequence);

  identifier = new int[nb_sequence];
  file.Read(identifier , nb_sequence);

  file.Read(max_length);
  file.Read(cumul_length);

  length = new int[nb_sequence];
  file.Read(length , nb_sequence);

  hlength = new Histogram();
  hlength->restoreGuts(file);

  file.Read(status);
  if (status) {
    index_interval = new Histogram();
    index_interval->restoreGuts(file);
  }
  else {
    index_interval = 0;
  }

  sequence = new int**[nb_sequence];
  for (i = 0;i < nb_sequence;i++) {
    sequence[i] = new int*[nb_variable];
    for (j = 0;j < nb_variable;j++) {
      blength = (type[j] == POSITION ? length[i] + 1 : length[i]);
      sequence[i][j] = new int[blength];
      file.Read(sequence[i][j] , blength);
    }
  }
}


void Sequences::saveGuts(RWvostream &os) const

{
  register int i , j , k;
  int *psequence;


  os << nb_variable;

  for (i = 0;i < nb_variable;i++) {
    os << type[i];
  }

  for (i = 0;i < nb_variable;i++) {
    os << min_value[i];
  }
  for (i = 0;i < nb_variable;i++) {
    os << max_value[i];
  }

  for (i = 0;i < nb_variable;i++) {
    if (marginal[i]) {
      os << true;
      marginal[i]->saveGuts(os);
    }
    else {
      os << false;
    }
  }

  os << nb_sequence;

  for (i = 0;i < nb_sequence;i++) {
    os << identifier[i];
  }

  os << max_length << cumul_length;

  for (i = 0;i < nb_sequence;i++) {
    os << length[i];
  }
  hlength->saveGuts(os);

  if (index_interval) {
    os << true;
    index_interval->saveGuts(os);
  }
  else {
    os << false;
  }

  for (i = 0;i < nb_sequence;i++) {
    for (j = 0;j < nb_variable;j++) {
      psequence = sequence[i][j];
      for (k = 0;k < (type[j] == POSITION ? length[i] + 1 : length[i]);k++) {
        os << *psequence++;
      }
    }
  }
}


void Sequences::saveGuts(RWFile &file) const

{
  register int i , j;


  file.Write(nb_variable);

  file.Write(type , nb_variable);

  file.Write(min_value , nb_variable);
  file.Write(max_value , nb_variable);

  for (i = 0;i < nb_variable;i++) {
    if (marginal[i]) {
      file.Write(true);
      marginal[i]->saveGuts(file);
    }
    else {
      file.Write(false);
    }
  }

  file.Write(nb_sequence);

  file.Write(identifier , nb_sequence);

  file.Write(max_length);
  file.Write(cumul_length);

  file.Write(length , nb_sequence);
  hlength->saveGuts(file);

  if (index_interval) {
    file.Write(true);
    index_interval->saveGuts(file);
  }
  else {
    file.Write(false);
  }

  for (i = 0;i < nb_sequence;i++) {
    for (j = 0;j < nb_variable;j++) {
      file.Write(sequence[i][j] , (type[j] == POSITION ? length[i] + 1 : length[i]));
    }
  }
} */
