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

#include <string>
#include <sstream>
#include <iomanip>

#include <boost/tokenizer.hpp>

#include "vectors.h"
#include "stat_label.h"

using namespace std;
using namespace boost;


namespace stat_tool {



/*--------------------------------------------------------------*/
/**
 *  \brief Default constructor of the VectorDistance class.
 */
/*--------------------------------------------------------------*/

VectorDistance::VectorDistance()

{
  nb_variable = 0;
  distance_type = ABSOLUTE_VALUE;

  var_type = NULL;
  weight = NULL;
  dispersion = NULL;

  nb_value = NULL;
  category_distance = NULL;

  period = NULL;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Constructor of the VectorDistance class.
 *
 *  \param[in] inb_variable   number of variables,
 *  \param[in] ivar_type      variable types,
 *  \param[in] iweight        variable weights,
 *  \param[in] idistance_type distance type (ABSOLUTE_VALUE/QUADRATIC).
 */
/*--------------------------------------------------------------*/

VectorDistance::VectorDistance(int inb_variable , variable_type *ivar_type ,
                               double *iweight , metric idistance_type)

{
  int i;


  nb_variable = inb_variable;
  distance_type = idistance_type;

  var_type = new variable_type[nb_variable];
  weight = new double[nb_variable];
  dispersion = new double[nb_variable];
  nb_value = new int[nb_variable];

  for (i = 0;i < nb_variable;i++) {
    var_type[i] = ivar_type[i];

    if (iweight) {
      weight[i] = iweight[i];
    }
    else {
      weight[i] = 1. / (double)nb_variable;
    }

    dispersion[i] = D_DEFAULT;
    nb_value[i] = 0;
  }

  category_distance = new double**[nb_variable];
  for (i = 0;i < nb_variable;i++) {
    category_distance[i] = NULL;
  }

  period = new int[nb_variable];
  for (i = 0;i < nb_variable;i++) {
    period[i] = I_DEFAULT;
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Constructor of the VectorDistance class.
 *
 *  \param[in] inb_variable       number of variables,
 *  \param[in] idistance_type     distance type (ABSOLUTE_VALUE/QUADRATIC),
 *  \param[in] ivar_type          variable types,
 *  \param[in] iweight            variable weights,
 *  \param[in] inb_value          number of categories (for categorical variables),
 *  \param[in] icategory_distance between-category distance matrices (for categorical variables),
 *  \param[in] iperiod            periods (for circular variables).
 */
/*--------------------------------------------------------------*/

VectorDistance::VectorDistance(int inb_variable , metric idistance_type , variable_type *ivar_type ,
                               double *iweight , int *inb_value , double ***icategory_distance ,
                               int *iperiod)

{
  int i , j , k;
  double *pdistance , *cdistance;


  nb_variable = inb_variable;
  distance_type = idistance_type;

  var_type = new variable_type[nb_variable];
  weight = new double[nb_variable];
  dispersion = new double[nb_variable];
  nb_value = new int[nb_variable];

  for (i = 0;i < nb_variable;i++) {
    var_type[i] = ivar_type[i];
    weight[i] = iweight[i];
    dispersion[i] = D_DEFAULT;
    nb_value[i] = inb_value[i];
  }

  category_distance = new double**[nb_variable];

  for (i = 0;i < nb_variable;i++) {
    if (icategory_distance[i]) {
      category_distance[i] = new double*[nb_value[i]];
      for (j = 0;j < nb_value[i];j++) {
        category_distance[i][j] = new double[nb_value[i]];

        pdistance = category_distance[i][j];
        cdistance = icategory_distance[i][j];
        for (k = 0;k < nb_value[i];k++) {
          *pdistance++ = *cdistance++;
        }
      }
    }

    else {
      category_distance[i] = NULL;
    }
  }

  period = new int[nb_variable];
  for (i = 0;i < nb_variable;i++) {
    period[i] = iperiod[i];
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Copy of a VectorDistance object.
 *
 *  \param[in] vector_dist reference on a VectorDistance object.
 */
/*--------------------------------------------------------------*/

void VectorDistance::copy(const VectorDistance &vector_dist)

{
  int i , j , k;
  double *pdistance , *cdistance;


  nb_variable = vector_dist.nb_variable;
  distance_type = vector_dist.distance_type;

  var_type = new variable_type[nb_variable];
  weight = new double[nb_variable];
  dispersion = new double[nb_variable];
  nb_value = new int[nb_variable];

  for (i = 0;i < nb_variable;i++) {
    var_type[i] = vector_dist.var_type[i];
    weight[i] = vector_dist.weight[i];
    dispersion[i] = vector_dist.dispersion[i];
    nb_value[i] = vector_dist.nb_value[i];
  }

  category_distance = new double**[nb_variable];

  for (i = 0;i < nb_variable;i++) {
    if (vector_dist.category_distance[i]) {
      category_distance[i] = new double*[nb_value[i]];

      for (j = 0;j < nb_value[i];j++) {
        category_distance[i][j] = new double[nb_value[i]];

        pdistance = category_distance[i][j];
        cdistance = vector_dist.category_distance[i][j];
        for (k = 0;k < nb_value[i];k++) {
          *pdistance++ = *cdistance++;
        }
      }
    }

    else {
      category_distance[i] = NULL;
    }
  }

  period = new int[nb_variable];
  for (i = 0;i < nb_variable;i++) {
    period[i] = vector_dist.period[i];
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Destruction of data members of a VectorDistance object.
 */
/*--------------------------------------------------------------*/

void VectorDistance::remove()

{
  int i , j;


  delete [] var_type;
  delete [] weight;
  delete [] dispersion;

  if (category_distance) {
    for (i = 0;i < nb_variable;i++) {
      if (category_distance[i]) {
        for (j = 0;j < nb_value[i];j++) {
          delete [] category_distance[i][j];
        }
        delete [] category_distance[i];
      }
    }
    delete [] category_distance;
  }

  delete [] nb_value;

  delete [] period;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Destructor of the VectorDistance class.
 */
/*--------------------------------------------------------------*/

VectorDistance::~VectorDistance()

{
  remove();
}


/*--------------------------------------------------------------*/
/**
 *  \brief Assignment operator of the VectorDistance class.
 *
 *  \param[in] vector_dist reference on a VectorDistance object.
 *
 *  \return                VectorDistance object.
 */
/*--------------------------------------------------------------*/

VectorDistance& VectorDistance::operator=(const VectorDistance &vector_dist)

{
  if (&vector_dist != this) {
    remove();
    copy(vector_dist);
  }

  return *this;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Construction of a VectorDistance object from a file.
 *
 *  \param[in] error reference on a StatError object,
 *  \param[in] path  file path.
 *
 *  \return          VectorDistance object.
 */
/*--------------------------------------------------------------*/

VectorDistance* VectorDistance::ascii_read(StatError &error , const string path)

{
  string buffer;
  size_t position;
  typedef tokenizer<char_separator<char>> tokenizer;
  char_separator<char> separator(" \t");
  bool status , lstatus , weight_flag;
  int i , j , k , m;
  int line , read_line , read_type , nb_variable , value , variable , category ,
      nb_value[VECTOR_NB_VARIABLE] , period[VECTOR_NB_VARIABLE];
  metric distance_type = ABSOLUTE_VALUE;
  variable_type var_type[VECTOR_NB_VARIABLE];
  double cumul , weight[VECTOR_NB_VARIABLE] , ***category_distance;
  VectorDistance *vector_dist;
  ifstream in_file(path.c_str());


  vector_dist = NULL;
  error.init();

  if (!in_file) {
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {
    status = true;
    line = 0;
    read_line = 0;

    // analysis of the line defining the number of variables

    while (getline(in_file , buffer)) {
      line++;

#     ifdef DEBUG
      cout << line << "  " << buffer << endl;
#     endif

      position = buffer.find('#');
      if (position != string::npos) {
        buffer.erase(position);
      }
      i = 0;

      tokenizer tok_buffer(buffer , separator);

      for (tokenizer::iterator token = tok_buffer.begin();token != tok_buffer.end();token++) {
        switch (i) {

        // test number of variables

        case 0 : {
          lstatus = true;

/*          try {
            value = stoi(*token);   in C++ 11
          }
          catch(invalid_argument &arg) {
            lstatus = false;
          } */
          value = atoi(token->c_str());

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

        // test VARIABLE(S) keyword

        case 1 : {
          if (*token != STAT_word[nb_variable == 1 ? STATW_VARIABLE : STATW_VARIABLES]) {
            status = false;
            error.correction_update(STAT_parsing[STATP_KEYWORD] ,
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

      // analysis of the line defining the distance type

      while (getline(in_file , buffer)) {
        line++;

#       ifdef DEBUG
        cout << line << "  " << buffer << endl;
#       endif

        position = buffer.find('#');
        if (position != string::npos) {
          buffer.erase(position);
        }
        i = 0;

        tokenizer tok_buffer(buffer , separator);

        for (tokenizer::iterator token = tok_buffer.begin();token != tok_buffer.end();token++) {
          switch (i) {

          // test DISTANCE keyword

          case 0 : {
            if (*token != STAT_word[STATW_DISTANCE]) {
              status = false;
              error.correction_update(STAT_parsing[STATP_KEYWORD] , STAT_word[STATW_DISTANCE] , line , i + 1);
            }
            break;
          }

          // test separator

          case 1 : {
            if (*token != ":") {
              status = false;
              error.update(STAT_parsing[STATP_SEPARATOR] , line , i + 1);
            }
            break;
          }

          // test keyword defining the distance type

          case 2 : {
            for (j = ABSOLUTE_VALUE;j <= QUADRATIC;j++) {
              if (*token == STAT_distance_word[j]) {
                distance_type = (metric)j;
                break;
              }
            }

            if (j == QUADRATIC + 1) {
              status = false;
              error.update(STAT_parsing[STATP_KEYWORD] , line , i + 1);
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

      // analysis of the lines defining the type of each variable and
      // the between-category distance matrices

      category_distance = new double**[nb_variable];
      for (i = 0;i < nb_variable;i++) {
        var_type[i] = NUMERIC;
        weight[i] = 1. / (double)nb_variable;
        nb_value[i] = 0;
        category_distance[i] = NULL;
        period[i] = I_DEFAULT;
      }

      read_type = 0;
      variable = -1;
      weight_flag = false;
      cumul = 0.;

      while (getline(in_file , buffer)) {
        line++;

#       ifdef DEBUG
        cout << line << "  " << buffer << endl;
#       endif

        position = buffer.find('#');
        if (position != string::npos) {
          buffer.erase(position);
        }
        i = 0;

        tokenizer tok_buffer(buffer , separator);

        for (tokenizer::iterator token = tok_buffer.begin();token != tok_buffer.end();token++) {
          if ((i == 0) && (read_type == 0)) {
            if (*token == STAT_word[STATW_VARIABLE]) {
              variable++;
            }

            else if ((variable >= 0) && (variable < nb_variable) &&
                     (var_type[variable] == NOMINAL)) {
              lstatus = true;

/*              try {
                value = stoi(*token);   in C++ 11
              }
              catch(invalid_argument &arg) {
                lstatus = false;
              } */
              value = atoi(token->c_str());

              if (lstatus) {
                read_type = 1;

                if ((value < 2) || (value > NB_CATEGORY)) {
                  lstatus = false;
                  nb_value[variable] = NB_CATEGORY;
                }
                else {
                  nb_value[variable] = value;
                }
              }

              if (!lstatus) {
                status = false;
                ostringstream correction_message;
                correction_message << "between 2 and " << NB_CATEGORY;
                error.correction_update(STAT_parsing[STATP_NB_CATEGORY] , (correction_message.str()).c_str() ,
                                        line , i + 1);
              }
            }
          }

          switch (read_type) {

          case 0 : {

            switch (i) {

            // test variable index

            case 1 : {
              lstatus = true;

/*              try {
                value = stoi(*token);   in C++ 11
              }
              catch(invalid_argument &arg) {
                lstatus = false;
              } */
              value = atoi(token->c_str());

              if ((lstatus) && (value != variable + 1)) {
                lstatus = false;
              }

              if (!lstatus) {
                status = false;
                error.correction_update(STAT_parsing[STATP_VARIABLE_INDEX] , variable + 1 , line , i + 1);
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

            // test keyword corresponding to the variable type

            case 3 : {
              if ((variable >= 0) && (variable < nb_variable)) {
                for (j = NOMINAL;j <= CIRCULAR;j++) {
                  if (*token == STAT_variable_type_word[j]) {
                    var_type[variable] = (variable_type)j;
                    break;
                  }
                }

                if (j == CIRCULAR + 1) {
                  status = false;
                  error.update(STAT_parsing[STATP_KEYWORD] , line , i + 1);
                }
              }
              break;
            }

            // test WEIGHT keyword

            case 4 : {
              if (*token != STAT_word[STATW_WEIGHT]) {
                status = false;
                error.correction_update(STAT_parsing[STATP_KEYWORD] , STAT_word[STATW_WEIGHT] , line , i + 1);
              }
              break;
            }

            // test separator

            case 5 : {
              if (*token != ":") {
                status = false;
                error.update(STAT_parsing[STATP_SEPARATOR] , line , i + 1);
              }
              break;
            }

            // test weight value

            case 6 : {
              if ((variable >= 0) && (variable < nb_variable)) {
                lstatus = true;

/*                try {
                  weight[variable] = stod(*token);   in C++ 11
                }
                catch (invalid_argument &arg) {
                  lstatus = false;
                } */
                weight[variable] = atof(token->c_str());

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

            // test CATEGORIES keyword

            if ((i == 1) && (*token != STAT_word[STATW_CATEGORIES])) {
              status = false;
              error.update(STAT_parsing[STATP_KEYWORD] , line , i + 1);
            }

            break;
          }

          case 2 : {
            if (i <= MIN(category , nb_value[variable] - 1)) {
              lstatus = true;

/*              try {
                category_distance[variable][category][i] = stod(*token);   in C++ 11
              }
              catch (invalid_argument &arg) {
                lstatus = false;
              } */
              category_distance[variable][category][i] = atof(token->c_str());

              if (lstatus) {
                if (((i < category) && (category_distance[variable][category][i] <= 0.)) ||
                    ((i == category) && (category_distance[variable][category][i] != 0.))) {
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

            // test PERIOD keyword

            case 0 : {
              if (*token != STAT_word[STATW_PERIOD]) {
                status = false;
                error.correction_update(STAT_parsing[STATP_KEYWORD] , STAT_word[STATW_PERIOD] , line , i + 1);
              }
              break;
            }

            // test separator

            case 1 : {
              if (*token != ":") {
                status = false;
                error.update(STAT_parsing[STATP_SEPARATOR] , line , i + 1);
              }
              break;
            }

            // test period value

            case 2 : {
              if ((variable >= 0) && (variable < nb_variable)) {
                lstatus = true;

/*                try {
                  value = stoi(*token);   in C++ 11
                }
                catch(invalid_argument &arg) {
                  lstatus = false;
                } */
                value = atoi(token->c_str());
 
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

            if ((variable >= 0) && (variable < nb_variable) && (var_type[variable] == CIRCULAR)) {
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
            category = 0;

            category_distance[variable] = new double*[nb_value[variable]];
            for (j = 0;j < nb_value[variable];j++) {
              category_distance[variable][j] = new double[nb_value[variable]];
            }
            break;
          }

          case 2 : {
            if (i != category + 1) {
              status = false;
              error.update(STAT_parsing[STATP_FORMAT] , line);
            }

            category++;
            if (category == nb_value[variable]) {
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
          if (category_distance[i]) {
            for (j = 0;j < nb_value[i];j++) {
              for (k = j + 1;k < nb_value[i];k++) {
                category_distance[i][j][k] = category_distance[i][k][j];
              }
            }

            // checking of the triagular inequality

            for (j = 0;j < nb_value[i];j++) {
              for (k = j + 1;k < nb_value[i];k++) {
                for (m = 0;m < nb_value[i];m++) {
                  if ((m != j) && (m != k) &&
                      (category_distance[i][j][m] + category_distance[i][m][k] < category_distance[i][j][k])) {
                    status = false;
                    ostringstream error_message;
                    error_message << STAT_label[STATL_VARIABLE] << " " << i + 1 << ": "
                                  << STAT_label[STATL_CATEGORY] << " " << j << ", "
                                  << STAT_label[STATL_CATEGORY] << " " << m << ", "
                                  << STAT_label[STATL_CATEGORY] << " " << k << ": "
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
        vector_dist = new VectorDistance(nb_variable , distance_type , var_type ,
                                         weight , nb_value , category_distance , period);
      }
    }
  }

  return vector_dist;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing on a single line of a VectorDistance object.
 *
 *  \param[in,out] os stream.
 */
/*--------------------------------------------------------------*/

ostream& VectorDistance::line_write(ostream &os) const

{
  os << nb_variable << " " << STAT_word[nb_variable == 1 ? STATW_VARIABLE : STATW_VARIABLES];
  if (nb_variable > 1) {
    os << "   " << STAT_word[STATW_DISTANCE] << " : " << STAT_distance_word[distance_type];
  }

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a VectorDistance object.
 *
 *  \param[in,out] os         stream,
 *  \param[in]     exhaustive flag detail level.
 */
/*--------------------------------------------------------------*/

ostream& VectorDistance::ascii_write(ostream &os , bool exhaustive) const

{
  int i , j , k;
  int buff , width[2];
  ios_base::fmtflags format_flags;


  format_flags = os.setf(ios::right , ios::adjustfield);

  os << nb_variable << " " << STAT_word[nb_variable == 1 ? STATW_VARIABLE : STATW_VARIABLES] << endl;
  if (nb_variable > 1) {
    os << "\n" << STAT_word[STATW_DISTANCE] << " : " << STAT_distance_word[distance_type] << endl;
  }

  for (i = 0;i < nb_variable;i++) {
    os << "\n" << STAT_word[STATW_VARIABLE] << " " << i + 1 << " : "
       << STAT_variable_type_word[var_type[i]];
    if (nb_variable > 1) {
      os << "   " << STAT_word[STATW_WEIGHT] << " : " << weight[i];
    }
    os << endl;

    switch (var_type[i]) {

    case NOMINAL : {
      if (category_distance[i]) {
        os << "\n" << nb_value[i] << " " << STAT_word[STATW_CATEGORIES] << endl;

        // computation of the column widths

        width[0] = column_width(nb_value[i]);

        width[1] = 0;
        for (j = 0;j < nb_value[i];j++) {
          buff = column_width(nb_value[i] , category_distance[i][j]);
          if (buff > width[1]) {
            width[1] = buff;
          }
        }
        width[1] += ASCII_SPACE;

        // writing of the between-category distance matrix

        os << "\n" << setw(width[0] + width[1]) << 0;
        for (j = 1;j < nb_value[i];j++) {
          os << setw(width[1]) << j;
        }
        os << endl;

        for (j = 0;j < nb_value[i];j++) {
          os << setw(width[0]) << j;
          for (k = 0;k < nb_value[i];k++) {
            os << setw(width[1]) << category_distance[i][j][k];
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

  os.setf(format_flags , ios::adjustfield);

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Functions for compatibility with the StatInterface class.
 */
/*--------------------------------------------------------------*/

bool VectorDistance::ascii_write(StatError &error , const string path , bool exhaustive) const

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
    ascii_write(out_file , exhaustive );
  }

  return status;
}


bool VectorDistance::spreadsheet_write(StatError &error , const string path) const

{
  return false;
}


bool VectorDistance::plot_write(StatError &error , const char *prefix ,
                                const char *title) const

{
  return false;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the maximum distance for each category
 *
 *  \param[in] variable variable index.
 *
 *  \return             maximum distances.
 */
/*--------------------------------------------------------------*/

double* VectorDistance::max_category_distance_computation(int variable) const

{
  int i , j;
  double *max_category_distance = NULL;


  if ((var_type[variable] == NOMINAL) && (category_distance[variable])) {
    max_category_distance = new double[nb_value[variable]];

    for (i = 0;i < nb_value[variable];i++) {
      max_category_distance[i] = 0.;
      for (j = 0;j < nb_value[variable];j++) {
        if (category_distance[variable][i][j] > max_category_distance[i]) {
          max_category_distance[i] = category_distance[variable][i][j];
        }
      }
    }
  }

  return max_category_distance;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Update of the standardization quantity.
 *
 *  \param[in] variable    variable index,
 *  \param[in] idispersion variable dispersion.
 */
/*--------------------------------------------------------------*/

void VectorDistance::dispersion_update(int variable , double idispersion) const

{
  dispersion[variable] = idispersion;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the standardization quantity.
 *
 *  \param[in] variable              variable index,
 *  \param[in] marginal_distribution marginal frequency distribution,
 *  \param[in] rank                  ranks.
 */
/*--------------------------------------------------------------*/

void VectorDistance::dispersion_computation(int variable , const FrequencyDistribution *marginal_distribution ,
                                            double *rank) const

{
  int i , j;
  int *pfrequency;
  double sum , distance = 1.;


  dispersion[variable] = 0.;

  for (i = marginal_distribution->offset;i < marginal_distribution->nb_value;i++) {
    if (marginal_distribution->frequency[i] > 0) {
      pfrequency = marginal_distribution->frequency + i + 1;
      sum = 0.;
      for (j = i + 1;j < marginal_distribution->nb_value;j++) {
        if (*pfrequency > 0) {
          switch (var_type[variable]) {

          case NOMINAL : {
            if (category_distance[variable]) {
              distance = category_distance[variable][i][j];
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

  switch (var_type[variable]) {

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
           << sqrt(2.) * marginal_distribution->mean_absolute_deviation_computation(marginal_distribution->mean);
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
