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
 *       $Id: vector_distance.cpp 18023 2015-04-23 07:08:20Z guedon $
 *
 *       Forum for V-Plants developers
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

#include "stat_tools.h"
#include "markovian.h"
#include "vectors.h"
#include "stat_label.h"

using namespace std;


namespace stat_tool {



/*--------------------------------------------------------------*
 *
 *  Constructeur par defaut de la classe VectorDistance.
 *
 *--------------------------------------------------------------*/

VectorDistance::VectorDistance()

{
  nb_variable = 0;
  distance_type = I_DEFAULT;

  variable_type = NULL;
  weight = NULL;
  dispersion = NULL;

  nb_value = NULL;
  symbol_distance = NULL;

  period = NULL;
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe VectorDistance.
 *
 *  arguments : nombre de variables, type et poids de chaque variable,
 *              type de distance (ABSOLUTE_VALUE/QUADRATIC).
 *
 *--------------------------------------------------------------*/

VectorDistance::VectorDistance(int inb_variable , int *ivariable_type ,
                               double *iweight , int idistance_type)

{
  register int i;


  nb_variable = inb_variable;
  distance_type = idistance_type;

  variable_type = new int[nb_variable];
  weight = new double[nb_variable];
  dispersion = new double[nb_variable];
  nb_value = new int[nb_variable];

  for (i = 0;i < nb_variable;i++) {
    variable_type[i] = ivariable_type[i];

    if (iweight) {
      weight[i] = iweight[i];
    }
    else {
      weight[i] = 1. / (double)nb_variable;
    }

    dispersion[i] = D_DEFAULT;
    nb_value[i] = 0;
  }

  symbol_distance = new double**[nb_variable];
  for (i = 0;i < nb_variable;i++) {
    symbol_distance[i] = NULL;
  }

  period = new int[nb_variable];
  for (i = 0;i < nb_variable;i++) {
    period[i] = I_DEFAULT;
  }
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe VectorDistance.
 *
 *  arguments : nombre de variables, type de distance (ABSOLUTE_VALUE/QUADRATIC),
 *              type et poids de chaque variable, nombre de valeurs,
 *              matrice des distances entre symboles, periodes.
 *
 *--------------------------------------------------------------*/

VectorDistance::VectorDistance(int inb_variable , int idistance_type , int *ivariable_type ,
                               double *iweight , int *inb_value , double ***isymbol_distance ,
                               int *iperiod)

{
  register int i , j , k;
  double *pdistance , *cdistance;


  nb_variable = inb_variable;
  distance_type = idistance_type;

  variable_type = new int[nb_variable];
  weight = new double[nb_variable];
  dispersion = new double[nb_variable];
  nb_value = new int[nb_variable];

  for (i = 0;i < nb_variable;i++) {
    variable_type[i] = ivariable_type[i];
    weight[i] = iweight[i];
    dispersion[i] = D_DEFAULT;
    nb_value[i] = inb_value[i];
  }

  symbol_distance = new double**[nb_variable];

  for (i = 0;i < nb_variable;i++) {
    if (isymbol_distance[i]) {
      symbol_distance[i] = new double*[nb_value[i]];
      for (j = 0;j < nb_value[i];j++) {
        symbol_distance[i][j] = new double[nb_value[i]];

        pdistance = symbol_distance[i][j];
        cdistance = isymbol_distance[i][j];
        for (k = 0;k < nb_value[i];k++) {
          *pdistance++ = *cdistance++;
        }
      }
    }

    else {
      symbol_distance[i] = NULL;
    }
  }

  period = new int[nb_variable];
  for (i = 0;i < nb_variable;i++) {
    period[i] = iperiod[i];
  }
}


/*--------------------------------------------------------------*
 *
 *  Copie d'un objet VectorDistance.
 *
 *  argument : reference sur un objet VectorDistance.
 *
 *--------------------------------------------------------------*/

void VectorDistance::copy(const VectorDistance &vector_dist)

{
  register int i , j , k;
  double *pdistance , *cdistance;


  nb_variable = vector_dist.nb_variable;
  distance_type = vector_dist.distance_type;

  variable_type = new int[nb_variable];
  weight = new double[nb_variable];
  dispersion = new double[nb_variable];
  nb_value = new int[nb_variable];

  for (i = 0;i < nb_variable;i++) {
    variable_type[i] = vector_dist.variable_type[i];
    weight[i] = vector_dist.weight[i];
    dispersion[i] = vector_dist.dispersion[i];
    nb_value[i] = vector_dist.nb_value[i];
  }

  symbol_distance = new double**[nb_variable];

  for (i = 0;i < nb_variable;i++) {
    if (vector_dist.symbol_distance[i]) {
      symbol_distance[i] = new double*[nb_value[i]];

      for (j = 0;j < nb_value[i];j++) {
        symbol_distance[i][j] = new double[nb_value[i]];

        pdistance = symbol_distance[i][j];
        cdistance = vector_dist.symbol_distance[i][j];
        for (k = 0;k < nb_value[i];k++) {
          *pdistance++ = *cdistance++;
        }
      }
    }

    else {
      symbol_distance[i] = NULL;
    }
  }

  period = new int[nb_variable];
  for (i = 0;i < nb_variable;i++) {
    period[i] = vector_dist.period[i];
  }
}


/*--------------------------------------------------------------*
 *
 *  Destruction des champs d'un objet VectorDistance.
 *
 *--------------------------------------------------------------*/

void VectorDistance::remove()

{
  register int i , j;


  delete [] variable_type;
  delete [] weight;
  delete [] dispersion;

  if (symbol_distance) {
    for (i = 0;i < nb_variable;i++) {
      if (symbol_distance[i]) {
        for (j = 0;j < nb_value[i];j++) {
          delete [] symbol_distance[i][j];
        }
        delete [] symbol_distance[i];
      }
    }
    delete [] symbol_distance;
  }

  delete [] nb_value;

  delete [] period;
}


/*--------------------------------------------------------------*
 *
 *  Destructeur de la classe VectorDistance.
 *
 *--------------------------------------------------------------*/

VectorDistance::~VectorDistance()

{
  remove();
}


/*--------------------------------------------------------------*
 *
 *  Operateur d'assignement de la classe VectorDistance.
 *
 *  argument : reference sur un objet VectorDistance.
 *
 *--------------------------------------------------------------*/

VectorDistance& VectorDistance::operator=(const VectorDistance &vector_dist)

{
  if (&vector_dist != this) {
    remove();
    copy(vector_dist);
  }

  return *this;
}


/*--------------------------------------------------------------*
 *
 *  Construction d'un objet VectorDistance a partir d'un fichier.
 *
 *  arguments : reference sur un objet StatError, path.
 *
 *--------------------------------------------------------------*/

VectorDistance* vector_distance_ascii_read(StatError &error , const char *path)

{
  RWLocaleSnapshot locale("en");
  RWCString buffer , token;
  size_t position;
  bool status , lstatus , weight_flag;
  register int i , j , k , m;
  int line , read_line , read_type , nb_variable , distance_type = ABSOLUTE_VALUE ,
      variable , symbol , variable_type[VECTOR_NB_VARIABLE] , nb_value[VECTOR_NB_VARIABLE] ,
      period[VECTOR_NB_VARIABLE];
  long value;
  double cumul , weight[VECTOR_NB_VARIABLE] , ***symbol_distance;
  VectorDistance *vector_dist;
  ifstream in_file(path);


  vector_dist = NULL;
  error.init();

  if (!in_file) {
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {
    status = true;
    line = 0;
    read_line = 0;

    // analyse de la ligne definissant le nombre de variables

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
        switch (i) {

        // test nombre de variables

        case 0 : {
          lstatus = locale.stringToNum(token , &value);
          if (lstatus) {
            if ((value < 1) || (value > VECTOR_NB_VARIABLE)) {
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

    if ((status) && (nb_variable > 1)) {

      // analyse de la ligne definissant le type de distance

      while (buffer.readLine(in_file , false)) {
        line++;

#       ifdef DEBUG
        cout << line << "  " << buffer << endl;
#       endif

        position = buffer.first('#');
        if (position != RW_NPOS) {
          buffer.remove(position);
        }
        i = 0;

        RWCTokenizer next(buffer);

        while (!((token = next()).isNull())) {
          switch (i) {

          // test mot cle DISTANCE

          case 0 : {
            if (token != STAT_word[STATW_DISTANCE]) {
              status = false;
              error.correction_update(STAT_parsing[STATP_KEY_WORD] , STAT_word[STATW_DISTANCE] , line , i + 1);
            }
            break;
          }

          // test separateur

          case 1 : {
            if (token != ":") {
              status = false;
              error.update(STAT_parsing[STATP_SEPARATOR] , line , i + 1);
            }
            break;
          }

          // test mot cle correspondant au type de distance

          case 2 : {
            for (j = ABSOLUTE_VALUE;j <= QUADRATIC;j++) {
              if (token == STAT_distance_word[j]) {
                distance_type = j;
                break;
              }
            }

            if (j == QUADRATIC + 1) {
              status = false;
              error.update(STAT_parsing[STATP_KEY_WORD] , line , i + 1);
            }
            break;
          }
          }

          i++;
        }

        if (i > 0) {
          if (i != 3) {
            status = false;
            error.update(STAT_parsing[STATP_FORMAT] , line);
          }
          break;
        }
      }
    }

    if (status) {

      // analyse des lignes definissant le type de chaque variable et
      // les matrices de distances entre symboles

      symbol_distance = new double**[nb_variable];
      for (i = 0;i < nb_variable;i++) {
        variable_type[i] = NUMERIC;
        weight[i] = 1. / (double)nb_variable;
        nb_value[i] = 0;
        symbol_distance[i] = NULL;
        period[i] = I_DEFAULT;
      }

      read_type = 0;
      variable = -1;
      weight_flag = false;
      cumul = 0.;

      while (buffer.readLine(in_file , false)) {
        line++;

#       ifdef DEBUG
        cout << line << "  " << buffer << endl;
#       endif

        position = buffer.first('#');
        if (position != RW_NPOS) {
          buffer.remove(position);
        }
        i = 0;

        RWCTokenizer next(buffer);

        while (!((token = next()).isNull())) {
          if ((i == 0) && (read_type == 0)) {
            if (token == STAT_word[STATW_VARIABLE]) {
              variable++;
            }

            else if ((variable >= 0) && (variable < nb_variable) &&
                     (variable_type[variable] == SYMBOLIC)) {
              lstatus = locale.stringToNum(token , &value);

              if (lstatus) {
                read_type = 1;

                if ((value < 2) || (value > NB_SYMBOL)) {
                  lstatus = false;
                  nb_value[variable] = NB_SYMBOL;
                }
                else {
                  nb_value[variable] = value;
                }
              }

              if (!lstatus) {
                status = false;
                ostringstream correction_message;
                correction_message << "between 2 and " << NB_SYMBOL;
                error.correction_update(STAT_parsing[STATP_NB_SYMBOL] , (correction_message.str()).c_str() ,
                                        line , i + 1);
              }
            }
          }

          switch (read_type) {

          case 0 : {

            switch (i) {

            // test indice de la variable

            case 1 : {
              lstatus = locale.stringToNum(token , &value);
              if ((lstatus) && (value != variable + 1)) {
                lstatus = false;
              }

              if (!lstatus) {
                status = false;
                error.correction_update(STAT_parsing[STATP_VARIABLE_INDEX] , variable + 1 , line , i + 1);
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
              if ((variable >= 0) && (variable < nb_variable)) {
                for (j = SYMBOLIC;j <= CIRCULAR;j++) {
                  if (token == STAT_variable_type_word[j]) {
                    variable_type[variable] = j;
                    break;
                  }
                }

                if (j == CIRCULAR + 1) {
                  status = false;
                  error.update(STAT_parsing[STATP_KEY_WORD] , line , i + 1);
                }
              }
              break;
            }

            // test mot cle WEIGHT

            case 4 : {
              if (token != STAT_word[STATW_WEIGHT]) {
                status = false;
                error.correction_update(STAT_parsing[STATP_KEY_WORD] , STAT_word[STATW_WEIGHT] , line , i + 1);
              }
              break;
            }

            // test separateur

            case 5 : {
              if (token != ":") {
                status = false;
                error.update(STAT_parsing[STATP_SEPARATOR] , line , i + 1);
              }
              break;
            }

            // test valeur du poids

            case 6 : {
              if ((variable >= 0) && (variable < nb_variable)) {
                lstatus = locale.stringToNum(token , weight + variable);
                if (lstatus) {
                  if ((weight[variable] <= 0.) || (weight[variable] > 1. - cumul + DOUBLE_ERROR)) {
                    lstatus = false;
                  }
                  else {
                    cumul += weight[variable];
                  }
                }

                if (!lstatus) {
                  status = false;
                  error.update(STAT_parsing[STATP_WEIGHT_VALUE] , line , i + 1);
                }
              }
              break;
            }
            }

            break;
          }

          case 1 : {

            // test mot cle SYMBOLS

            if ((i == 1) && (token != STAT_word[STATW_SYMBOLS])) {
              status = false;
              error.update(STAT_parsing[STATP_KEY_WORD] , line , i + 1);
            }

            break;
          }

          case 2 : {
            if (i <= MIN(symbol , nb_value[variable] - 1)) {
              lstatus = locale.stringToNum(token , symbol_distance[variable][symbol] + i);
              if (lstatus) {
                if (((i < symbol) && (symbol_distance[variable][symbol][i] <= 0.)) ||
                    ((i == symbol) && (symbol_distance[variable][symbol][i] != 0.))) {
                  lstatus = false;
                }
              }

              if (!lstatus) {
                status = false;
                error.update(STAT_parsing[STATP_LOCAL_DISTANCE] , line , i + 1);
              }
            }

            break;
          }

          case 3 : {

            switch (i) {

            // test mot cle PERIOD

            case 0 : {
              if (token != STAT_word[STATW_PERIOD]) {
                status = false;
                error.correction_update(STAT_parsing[STATP_KEY_WORD] , STAT_word[STATW_PERIOD] , line , i + 1);
              }
              break;
            }

            // test separateur

            case 1 : {
              if (token != ":") {
                status = false;
                error.update(STAT_parsing[STATP_SEPARATOR] , line , i + 1);
              }
              break;
            }

            // test valeur de la periode

            case 2 : {
              if ((variable >= 0) && (variable < nb_variable)) {
                lstatus = locale.stringToNum(token , &value);
                if (lstatus) {
                  if (value < 2) {
                    lstatus = false;
                  }
                  else {
                    period[variable] = value;
                  }
                }

                if (!lstatus) {
                  status = false;
                  error.update(STAT_parsing[STATP_PERIOD_VALUE] , line , i + 1);
                }
              }
              break;
            }
            }

            break;
          }
          }

          i++;
        }

        if (i > 0) {
          switch (read_type) {

          case 0 : {
            if (variable == 0) {
              if (i == 7) {
                weight_flag = true;
              }
              else if (i != 4) {
                status = false;
                error.update(STAT_parsing[STATP_FORMAT] , line);
              }
            }

            else if (((!weight_flag) && (i != 4)) || ((weight_flag) && (i != 7))) {
              status = false;
              error.update(STAT_parsing[STATP_FORMAT] , line);
            }

            if ((variable >= 0) && (variable < nb_variable) &&
                (variable_type[variable] == CIRCULAR)) {
              read_type = 3;
            }
            break;
          }

          case 1 : {
            if (i != 2) {
              status = false;
              error.update(STAT_parsing[STATP_FORMAT] , line);
            }

            read_type = 2;
            symbol = 0;

            symbol_distance[variable] = new double*[nb_value[variable]];
            for (j = 0;j < nb_value[variable];j++) {
              symbol_distance[variable][j] = new double[nb_value[variable]];
            }
            break;
          }

          case 2 : {
            if (i != symbol + 1) {
              status = false;
              error.update(STAT_parsing[STATP_FORMAT] , line);
            }

            symbol++;
            if (symbol == nb_value[variable]) {
              read_type = 0;
            }
            break;
          }

          case 3 : {
            if (i != 3) {
              status = false;
              error.update(STAT_parsing[STATP_FORMAT] , line);
            }

            read_type = 0;
            break;
          }
          }
        }
      }

      if ((read_type != 0) || (variable != nb_variable - 1)) {
        status = false;
        error.update(STAT_parsing[STATP_FORMAT] , line);
      }

      if ((weight_flag) && (cumul < 1. - DOUBLE_ERROR)) {
        status = false;
        error.update(STAT_parsing[STATP_PROBABILITY_SUM]);
      }

      if (status) {
        for (i = 0;i < nb_variable;i++) {
          if (symbol_distance[i]) {
            for (j = 0;j < nb_value[i];j++) {
              for (k = j + 1;k < nb_value[i];k++) {
                symbol_distance[i][j][k] = symbol_distance[i][k][j];
              }
            }

            // verification de l'inegalite triagulaire

            for (j = 0;j < nb_value[i];j++) {
              for (k = j + 1;k < nb_value[i];k++) {
                for (m = 0;m < nb_value[i];m++) {
                  if ((m != j) && (m != k) &&
                      (symbol_distance[i][j][m] + symbol_distance[i][m][k] < symbol_distance[i][j][k])) {
                    status = false;
                    ostringstream error_message;
                    error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                                  << STAT_label[STATL_SYMBOL] << " " << j << ", "
                                  << STAT_label[STATL_SYMBOL] << " " << m << ", "
                                  << STAT_label[STATL_SYMBOL] << " " << k << ": "
                                  << STAT_parsing[STATP_TRIANGLE_INEQUALITY];
                    error.update((error_message.str()).c_str());
                  }
                }
              }
            }
          }
        }
      }

      if (status) {
        vector_dist = new VectorDistance(nb_variable , distance_type , variable_type ,
                                         weight , nb_value , symbol_distance , period);
      }
    }
  }

  return vector_dist;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture sur une ligne d'un objet VectorDistance.
 *
 *  argument : stream.
 *
 *--------------------------------------------------------------*/

ostream& VectorDistance::line_write(ostream &os) const

{
  os << nb_variable << " " << STAT_word[nb_variable == 1 ? STATW_VARIABLE : STATW_VARIABLES];
  if (nb_variable > 1) {
    os << "   " << STAT_word[STATW_DISTANCE] << " : " << STAT_distance_word[distance_type];
  }

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet VectorDistance.
 *
 *  arguments : stream, flag niveau de detail.
 *
 *--------------------------------------------------------------*/

ostream& VectorDistance::ascii_write(ostream &os , bool exhaustive) const

{
  register int i , j , k;
  int buff , width[2];
  long old_adjust;


  old_adjust = os.setf(ios::right , ios::adjustfield);

  os << nb_variable << " " << STAT_word[nb_variable == 1 ? STATW_VARIABLE : STATW_VARIABLES] << endl;
  if (nb_variable > 1) {
    os << "\n" << STAT_word[STATW_DISTANCE] << " : " << STAT_distance_word[distance_type] << endl;
  }

  for (i = 0;i < nb_variable;i++) {
    os << "\n" << STAT_word[STATW_VARIABLE] << " " << i + 1 << " : "
       << STAT_variable_type_word[variable_type[i]];
    if (nb_variable > 1) {
      os << "   " << STAT_word[STATW_WEIGHT] << " : " << weight[i];
    }
    os << endl;

    switch (variable_type[i]) {

    case SYMBOLIC : {
      if (symbol_distance[i]) {
        os << "\n" << nb_value[i] << " " << STAT_word[STATW_SYMBOLS] << endl;

        // calcul des largeurs des colonnes

        width[0] = column_width(nb_value[i]);

        width[1] = 0;
        for (j = 0;j < nb_value[i];j++) {
          buff = column_width(nb_value[i] , symbol_distance[i][j]);
          if (buff > width[1]) {
            width[1] = buff;
          }
        }
        width[1] += ASCII_SPACE;

        // ecriture de la matrice des distances entre symboles

        os << "\n" << setw(width[0] + width[1]) << 0;
        for (j = 1;j < nb_value[i];j++) {
          os << setw(width[1]) << j;
        }
        os << endl;

        for (j = 0;j < nb_value[i];j++) {
          os << setw(width[0]) << j;
          for (k = 0;k < nb_value[i];k++) {
            os << setw(width[1]) << symbol_distance[i][j][k];
          }
          os << endl;
        }
      }
      break;
    }

    case CIRCULAR : {
      os << STAT_word[STATW_PERIOD] << " : " << period[i] << endl;
      break;
    }
    }
  }

  os.setf((FMTFLAGS)old_adjust , ios::adjustfield);

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Fonctions pour la compatibilite avec la classe StatInterface.
 *
 *--------------------------------------------------------------*/

bool VectorDistance::ascii_write(StatError &error , const char *path ,
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
    ascii_write(out_file , exhaustive );
  }

  return status;
}


bool VectorDistance::spreadsheet_write(StatError &error , const char *path) const

{
  return false;
}


bool VectorDistance::plot_write(StatError &error , const char *prefix ,
                                const char *title) const

{
  return false;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la distance maximum pour chaque symbole
 *
 *  argument : indice de la variable.
 *
 *--------------------------------------------------------------*/

double* VectorDistance::max_symbol_distance_computation(int variable) const

{
  register int i , j;
  double *max_symbol_distance = NULL;


  if ((variable_type[variable] == SYMBOLIC) && (symbol_distance[variable])) {
    max_symbol_distance = new double[nb_value[variable]];

    for (i = 0;i < nb_value[variable];i++) {
      max_symbol_distance[i] = 0.;
      for (j = 0;j < nb_value[variable];j++) {
        if (symbol_distance[variable][i][j] > max_symbol_distance[i]) {
          max_symbol_distance[i] = symbol_distance[variable][i][j];
        }
      }
    }
  }

  return max_symbol_distance;
}


/*--------------------------------------------------------------*
 *
 *  Mise a jour de la quantite de standardisation.
 *
 *  arguments : indice de la variable, dispersion.
 *
 *--------------------------------------------------------------*/

void VectorDistance::dispersion_update(int variable , double idispersion) const

{
  dispersion[variable] = idispersion;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la quantite de standardisation.
 *
 *  arguments : indice de la variable, loi marginale empirique, rangs.
 *
 *--------------------------------------------------------------*/

void VectorDistance::dispersion_computation(int variable , const FrequencyDistribution *marginal_distribution ,
                                            double *rank) const

{
  register int i , j;
  int *pfrequency;
  double sum , distance = 1.;


  dispersion[variable] = 0.;

  for (i = marginal_distribution->offset;i < marginal_distribution->nb_value;i++) {
    if (marginal_distribution->frequency[i] > 0) {
      pfrequency = marginal_distribution->frequency + i + 1;
      sum = 0.;
      for (j = i + 1;j < marginal_distribution->nb_value;j++) {
        if (*pfrequency > 0) {
          switch (variable_type[variable]) {

          case SYMBOLIC : {
            if (symbol_distance[variable]) {
              distance = symbol_distance[variable][i][j];
            }
            break;
          }

          case ORDINAL : {
            distance = rank[j] - rank[i];
            break;
          }

          case NUMERIC : {
            distance = j - i;
            break;
          }

          case CIRCULAR : {
            distance = MIN(j - i , i + period[variable] - j);
            break;
          }
          }

          switch (distance_type) {
          case ABSOLUTE_VALUE :
            sum += *pfrequency * distance;
            break;
          case QUADRATIC :
            sum += *pfrequency * distance * distance;
            break;
          }
        }

        pfrequency++;
      }

      dispersion[variable] += marginal_distribution->frequency[i] * sum;
    }
  }

  dispersion[variable] = 2 * dispersion[variable] / (marginal_distribution->nb_element *
                          (double)(marginal_distribution->nb_element - 1));
  if (dispersion[variable] == 0.) {
    dispersion[variable] = 1.;
  }

# ifdef DEBUG
  cout << "\n" << STAT_label[STATL_VARIABLE] << " " << variable << "   dispersion: "
       << dispersion[variable];

  switch (variable_type[variable]) {

  case ORDINAL : {
    double dispersion2 = marginal_distribution->nb_element *
                         ((double)marginal_distribution->nb_element *
                          (double)marginal_distribution->nb_element - 1);

    pfrequency = marginal_distribution->frequency + marginal_distribution->offset;
    for (i = marginal_distribution->offset;i < marginal_distribution->nb_value;i++) {
      if (*pfrequency > 1) {
        dispersion2 -= *pfrequency * ((double)*pfrequency * (double)*pfrequency - 1);
      }
      pfrequency++;
    }

    switch (distance_type) {
    case ABSOLUTE_VALUE :
      dispersion2 /= (3 * marginal_distribution->nb_element *
                      (double)(marginal_distribution->nb_element - 1));
      break;
    case QUADRATIC :
      dispersion2 /= (6 * (marginal_distribution->nb_element - 1));
      break;
    }

    cout << " | " << dispersion2;

    if (distance_type == ABSOLUTE_VALUE) {
      int cumul;
      double previous_rank;

      pfrequency = marginal_distribution->frequency + marginal_distribution->offset;
      cumul = *pfrequency++;
      previous_rank = rank[marginal_distribution->offset];
      dispersion2 = 0.;

      for (i = marginal_distribution->offset + 1;i < marginal_distribution->nb_value;i++) {
        if (*pfrequency > 0) {
          dispersion2 += cumul * (double)(marginal_distribution->nb_element - cumul) *
                         (rank[i] - previous_rank);
          cumul += *pfrequency;
          previous_rank = rank[i];
        }
        pfrequency++;
      }

      dispersion2 = 2 * dispersion2 / (marginal_distribution->nb_element *
                     (double)(marginal_distribution->nb_element - 1));

      cout << " | " << dispersion2;
    }
    break;
  }

  case NUMERIC : {
    if (distance_type == ABSOLUTE_VALUE) {
      int previous_value , cumul;
      double dispersion2 = 0.;

      pfrequency = marginal_distribution->frequency + marginal_distribution->offset;
      cumul = *pfrequency++;
      previous_value = marginal_distribution->offset;

      for (i = marginal_distribution->offset + 1;i < marginal_distribution->nb_value;i++) {
        if (*pfrequency > 0) {
          dispersion2 += cumul * (double)(marginal_distribution->nb_element - cumul) *
                         (i - previous_value);
          cumul += *pfrequency;
          previous_value = i;
        }
        pfrequency++;
      }

      dispersion2 = 2 * dispersion2 / (marginal_distribution->nb_element *
                     (double)(marginal_distribution->nb_element - 1));

      cout << " | " << dispersion2;
    }

    cout << " | ";
    switch (distance_type) {
    case ABSOLUTE_VALUE :
      cout << "sqrt(2) * " << STAT_label[STATL_MEAN_ABSOLUTE_DEVIATION] << ": "
           << sqrt(2.) * marginal_distribution->mean_absolute_deviation_computation();
      break;
    case QUADRATIC :
      cout << 2 * marginal_distribution->variance;
      break;
    }
    break;
  }
  }

  cout << endl;
# endif

}


};  // namespace stat_tool
