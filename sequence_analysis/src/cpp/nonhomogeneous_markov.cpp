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
 *       $Id: nonhomogeneous_markov.cpp 3257 2007-06-06 12:56:12Z dufourko $
 *
 *       Forum for V-Plants developers: amldevlp@cirad.fr
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

#include "nonhomogeneous_markov.h"
#include "sequence_label.h"

using namespace std;
using namespace boost;
using namespace stat_tool;


namespace sequence_analysis {



/*--------------------------------------------------------------*/
/**
 *  \brief Default constructor of the Function class.
 */
/*--------------------------------------------------------------*/

Function::Function()

{
  residual = NULL;
  frequency = NULL;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Constructor of the Function class.
 *
 *  \param[in] iident     identifier,
 *  \param[in] length     length,
 *  \param[in] iparameter parameters.
 */
/*--------------------------------------------------------------*/

Function::Function(parametric_function iident , int length , double *iparameter)
:RegressionKernel(iident , 0 , length - 1)

{
  int i;


  for (i = 0;i < nb_parameter;i++) {
    parameter[i] = iparameter[i];
  }

  residual = NULL;
  frequency = NULL;

  computation();
}


/*--------------------------------------------------------------*/
/**
 *  \brief Constructor of the Function class.
 *
 *  \param[in] iident identifier,
 *  \param[in] length length.
 */
/*--------------------------------------------------------------*/

Function::Function(parametric_function iident , int length)
:RegressionKernel(iident , 0 , length - 1)

{
  int i;


  residual = new double[max_value + 1];
  for (i = 0;i <= max_value;i++) {
    residual[i] = -D_INF;
  }

  frequency = new int[max_value + 1];
  for (i = 0;i <= max_value;i++) {
    frequency[i] = 0;
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Copy of a Function object.
 *
 *  \param[in] function reference on a Function object.
 */
/*--------------------------------------------------------------*/

void Function::copy(const Function &function)

{
  if ((function.residual) && (function.frequency)) {
    int i;


    residual = new double[max_value + 1];
    for (i = 0;i <= max_value;i++) {
      residual[i] = function.residual[i];
    }

    frequency = new int[max_value + 1];
    for (i = 0;i <= max_value;i++) {
      frequency[i] = function.frequency[i];
    }
  }

  else {
    residual = NULL;
    frequency = NULL;
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Constructor by copy of the Function class.
 *
 *  \param[in] function reference on a Function object.
 */
/*--------------------------------------------------------------*/

Function::Function(const Function &function)

{
  RegressionKernel::copy(function);
  copy(function);
}


/*--------------------------------------------------------------*/
/**
 *  \brief Destruction of the data members of a Function object.
 */
/*--------------------------------------------------------------*/

void Function::remove()

{
  delete [] residual;
  delete [] frequency;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Destructor of the Function class.
 */
/*--------------------------------------------------------------*/

Function::~Function()

{
  remove();
}


/*--------------------------------------------------------------*/
/**
 *  \brief Assignment operator of the Function class.
 *
 *  \param[in] function reference on a Function object.
 *
 *  \return             Function object.
 */
/*--------------------------------------------------------------*/

Function& Function::operator=(const Function &function)

{
  if (&function != this) {
    remove();
    RegressionKernel::remove();

    RegressionKernel::copy(function);
    copy(function);
  }

  return *this;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Analysis of the format of a Function object.
 *
 *  \param[in] error   reference on a StatError object,
 *  \param[in] in_file stream,
 *  \param[in] line    reference on the file line index,
 *  \param[in] length  length,
 *  \param[in] min     lower bound.
 *  \param[in] max     upper bound.
 *
 *  \return            Function object.
 */
/*--------------------------------------------------------------*/

Function* Function::parsing(StatError &error , ifstream &in_file , int &line ,
                            int length , double min , double max)

{
  string buffer;
  size_t position;
  typedef tokenizer<char_separator<char>> tokenizer;
  char_separator<char> separator(" \t");
  bool status = true , lstatus;
  int i , j;
  int nb_parameter = 0 , index;
  parametric_function ident = NONPARAMETRIC_FUNCTION;
  double parameter[3];
  Function *function;


  function = NULL;

  while (getline(in_file , buffer)) {
    line++;

#   ifdef DEBUG
    cout << line << "  " << buffer << endl;
#   endif

    position = buffer.find('#');
    if (position != string::npos) {
      buffer.erase(position);
    }
    i = 0;

    tokenizer tok_buffer(buffer , separator);

    for (tokenizer::iterator token = tok_buffer.begin();token != tok_buffer.end();token++) {
      if (i <= 1) {
        switch (i) {

        // test LOGISTIC/MONOMOLECULAR keyword

        case 0 : {
          for (j = LOGISTIC;j <= MONOMOLECULAR;j++) {
            if (*token == STAT_function_word[j]) {
              ident = (parametric_function)j;
              break;
            }
          }

          if (j == MONOMOLECULAR + 1) {
            status = false;
            error.update(STAT_parsing[STATP_KEYWORD] , line , i + 1);
          }
          else {
            nb_parameter = 3;
            parameter[0] = D_DEFAULT;
          }
          break;
        }

        // test FUNCTION keyword

        case 1 : {
          if (*token != STAT_word[STATW_FUNCTION]) {
            status = false;
            error.correction_update(STAT_parsing[STATP_KEYWORD] , STAT_word[STATW_FUNCTION] , line , i + 1);
          }
          break;
        }
        }
      }

      else {
        switch ((i - 2) % 4) {

        // test PARAMETER keyword

        case 0 : {
          if (*token != STAT_word[STATW_PARAMETER]) {
            status = false;
            error.correction_update(STAT_parsing[STATP_KEYWORD] , STAT_word[STATW_PARAMETER] , line , i + 1);
          }
          break;
        }

        // test parameter index

        case 1 : {
          lstatus = true;

/*          try {
            index = stoi(*token);   in C++ 11
          }
          catch(invalid_argument &arg) {
            lstatus = false;
          } */
          index = atoi(token->c_str());

          if ((lstatus) && (index != (i - 2) / 4 + 1)) {
            lstatus = false;
          }

          if (!lstatus) {
            status = false;
            error.correction_update(STAT_parsing[STATP_PARAMETER_INDEX] , (i - 2) / 4 + 1 , line , i + 1);
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

        // test parameter value

        case 3 : {
          if ((i - 2) / 4 < nb_parameter) {
            lstatus = true;

/*            try {
              parameter[(i - 2) / 4] = stod(*token);   in C++ 11
            }
            catch(invalid_argument &arg) {
              lstatus = false;
            } */
            parameter[(i - 2) / 4] = atof(token->c_str());

            if (lstatus) {
              switch (ident) {

              case LOGISTIC : {
                switch ((i - 2) / 4) {

                case 0 : {
                  if ((parameter[0] < min) || (parameter[0] > max)) {
                    lstatus = false;
                  }
                  break;
                }

                case 1 : {
                  if ((parameter[0] != D_DEFAULT) &&
                      ((parameter[0] / (1. + parameter[1]) < min) || (parameter[0] / (1. + parameter[1]) > max))) {
                    lstatus = false;
                  }
                  break;
                }

                case 2 : {
                  if (parameter[2] <= 0.) {
                    lstatus = false;
                  }
                  break;
                }
                }

                break;
              }

              case MONOMOLECULAR : {
                switch ((i - 2) / 4) {

                case 0 : {
                  if ((parameter[0] < min) || (parameter[0] > max)) {
                    lstatus = false;
                  }
                  break;
                }

                case 1 : {
                  if ((parameter[0] != D_DEFAULT) &&
                      ((parameter[0] + parameter[1] < min) || (parameter[0] + parameter[1] > max))) {
                    lstatus = false;
                  }
                  break;
                }

                case 2 : {
                  if (parameter[2] <= 0.) {
                    lstatus = false;
                  }
                  break;
                }
                }

                break;
              }
              }
            }

            if (!lstatus) {
              status = false;
              error.update(STAT_parsing[STATP_PARAMETER_VALUE] , line , i + 1);
            }
          }
          break;
        }
        }
      }

      i++;
    }

    if (i > 0) {
      if (i != 14) {
        status = false;
        error.update(STAT_parsing[STATP_FORMAT] , line);
      }

      break;
    }
  }

  if (ident == NONPARAMETRIC_FUNCTION) {
    status = false;
    error.update(STAT_parsing[STATP_FORMAT] , line);
  }

  if (status) {
    function = new Function(ident , length , parameter);
  }

  return function;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a Function object and the associated Curves object.
 *
 *  \param[in,out] os         stream,
 *  \param[in]     exhaustive flag detail level,
 *  \param[in]     file_flag  flag file,
 *  \param[in]     curves     pointer on a Curves object.
 */
/*--------------------------------------------------------------*/

ostream& Function::ascii_print(ostream &os , bool exhaustive , bool file_flag ,
                               const Curves *curves) const

{
  int i;
  int *pfrequency , width[6];
  double self_transition_mean , residual_mean , residual_standard_deviation ,
         *standard_residual , square_sum[3];
  ios_base::fmtflags format_flags;


  format_flags = os.setf(ios::right , ios::adjustfield);

  ascii_parameter_print(os);
  os << endl;

  if (curves) {
    self_transition_mean = curves->mean_computation(0);
    square_sum[0] = regression_square_sum_computation(self_transition_mean);
    square_sum[1] = residual_square_sum_computation();
    square_sum[2] = curves->total_square_sum_computation(0 , self_transition_mean);

    os << "\n";
    if (file_flag) {
      os << "# ";
    }
    os << STAT_label[STATL_DETERMINATION_COEFF] << ": "
       << 1. - square_sum[1] / square_sum[2] << endl;

    if (file_flag) {
      os << "# ";
    }
    os << regression_df << " " << STAT_label[STATL_REGRESSION] << " " << STAT_label[STATL_FREEDOM_DEGREES] << "   "
       << residual_df << " " << STAT_label[STATL_RESIDUAL] << " " << STAT_label[STATL_FREEDOM_DEGREES] << endl;

#   ifdef DEBUG
    os << "\n" << STAT_label[STATL_REGRESSION] << " " << STAT_label[STATL_SQUARE_SUM] << ": " << square_sum[0]
       << "   " << STAT_label[STATL_RESIDUAL] << " " << STAT_label[STATL_SQUARE_SUM] << ": " << square_sum[1]
       << "   " << STAT_label[STATL_TOTAL] << " " << STAT_label[STATL_SQUARE_SUM] << ": " << square_sum[2] << endl;
#   endif

    // writing of the residual mean and standard deviation

    residual_mean = residual_mean_computation();
    residual_standard_deviation = sqrt(residual_variance_computation(residual_mean));

    os << "\n";
    if (file_flag) {
      os << "# ";
    }
    os << STAT_label[STATL_RESIDUAL] << " " << STAT_label[STATL_MEAN] << ": " << residual_mean << "   "
       << STAT_label[STATL_RESIDUAL] << " " << STAT_label[STATL_STANDARD_DEVIATION] << ": "
       << residual_standard_deviation << endl;
  }

  if (exhaustive) {
    if (curves) {

      // computation of the standardized residuals

      standard_residual = new double[max_value + 1];

      for (i = 0;i <= max_value;i++) {
        if (frequency[i] > 0) {
          standard_residual[i] = residual[i] / residual_standard_deviation;
        }
      }
    }

    // computation of the column widths

    width[0] = column_width(max_value);
    width[2] = column_width(max_value + 1 , point) + ASCII_SPACE;
    if (curves) {
      width[1] = column_width(curves->length , curves->point[0]) + ASCII_SPACE;
      width[3] = column_width(max_value + 1 , residual) + ASCII_SPACE;
      width[4] = column_width(max_value + 1 , standard_residual) + ASCII_SPACE;
      width[5] = column_width(curves->max_frequency_computation()) + ASCII_SPACE;
    }

    // writing of the observed and theoretical self-transition probabilities, of the residuals,
    // the standardized residuals and the frequencies

    os << "\n";
    if (file_flag) {
      os << "# ";
    }
    os << "  ";
    if (curves) {
      os << " | " << SEQ_label[SEQL_OBSERVED] << " " << SEQ_label[SEQL_SELF_TRANSITION];
    }
    os << " | " << SEQ_label[SEQL_THEORETICAL] << " " << SEQ_label[SEQL_SELF_TRANSITION];
    if (curves) {
      os << " | " << STAT_label[STATL_RESIDUAL] << " | " << STAT_label[STATL_STANDARDIZED_RESIDUAL]
         << " | " << STAT_label[STATL_FREQUENCY];
    }
    os << endl;

    for (i = 0;i <= max_value;i++) {
      if (file_flag) {
        os << "# ";
      }
      os << setw(width[0]) << i;

      if (curves) {
        if (frequency[i] > 0) {
          os << setw(width[1]) << curves->point[0][i];
        }
        else {
          os << setw(width[1]) << " ";
        }
      }

      os << setw(width[2]) << point[i];

      if (curves) {
        if (frequency[i] > 0) {
          os << setw(width[3]) << residual[i];
          os << setw(width[4]) << standard_residual[i];
        }
        else {
          os << setw(width[3]) << " ";
          os << setw(width[4]) << " ";
        }

        os << setw(width[5]) << frequency[i];
      }

      os << endl;
    }

    if (curves) {
      delete [] standard_residual;
    }
  }

  os.setf(format_flags , ios::adjustfield);

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a Function object and the associated Curves object
 *         at the spreadsheet format.
 *
 *  \param[in,out] os     stream,
 *  \param[in]     curves pointer on a Curves object.
 */
/*--------------------------------------------------------------*/

ostream& Function::spreadsheet_print(ostream &os , const Curves *curves) const

{
  int i;
  int *pfrequency;
  double self_transition_mean , residual_mean , residual_standard_deviation , square_sum[3];


  os << STAT_function_word[ident] << " " << STAT_word[STATW_FUNCTION];
  for (i = 0;i < nb_parameter;i++) {
    os << "\t\t" << STAT_word[STATW_PARAMETER] << " " << i + 1 << "\t" << parameter[i];
  }
  os << endl;

  if (curves) {
    self_transition_mean = curves->mean_computation(0);
    square_sum[0] = regression_square_sum_computation(self_transition_mean);
    square_sum[1] = residual_square_sum_computation();
    square_sum[2] = curves->total_square_sum_computation(0 , self_transition_mean);

    os << "\n" << STAT_label[STATL_DETERMINATION_COEFF] << "\t"
       << 1. - square_sum[1] / square_sum[2] << endl;

    os << regression_df << "\t" << STAT_label[STATL_REGRESSION] << " " << STAT_label[STATL_FREEDOM_DEGREES] << "\t\t"
       << residual_df << "\t" << STAT_label[STATL_RESIDUAL] << " " << STAT_label[STATL_FREEDOM_DEGREES] << endl;

#   ifdef DEBUG
    os << "\n" << STAT_label[STATL_REGRESSION] << " " << STAT_label[STATL_SQUARE_SUM] << "\t" << square_sum[0]
       << "\t\t" << STAT_label[STATL_RESIDUAL] << " " << STAT_label[STATL_SQUARE_SUM] << "\t" << square_sum[1]
       << "\t\t" << STAT_label[STATL_TOTAL] << " " << STAT_label[STATL_SQUARE_SUM] << "\t" << square_sum[2] << endl;
#   endif

    // writing of the residual mean and standard deviation

    residual_mean = residual_mean_computation();
    residual_standard_deviation = sqrt(residual_variance_computation(residual_mean));

    os << "\n" << STAT_label[STATL_RESIDUAL] << " " << STAT_label[STATL_MEAN] << "\t" << residual_mean
       << "\t\t" << STAT_label[STATL_RESIDUAL] << " " << STAT_label[STATL_STANDARD_DEVIATION] << "\t"
       << residual_standard_deviation << endl;
  }

  // writing of the observed and theoretical self-transition probabilities, of the residuals,
  // the standardized residuals and the frequencies

  os << "\n";
  if (curves) {
    os << "\t" << SEQ_label[SEQL_OBSERVED] << " " << SEQ_label[SEQL_SELF_TRANSITION];
  }
  os << "\t" << SEQ_label[SEQL_THEORETICAL] << " " << SEQ_label[SEQL_SELF_TRANSITION];
  if (curves) {
    os << "\t" << STAT_label[STATL_RESIDUAL] << "\t" << STAT_label[STATL_STANDARDIZED_RESIDUAL]
       << "\t" << SEQ_label[SEQL_ASYMPTOTE] << "\t" << SEQ_label[SEQL_ASYMPTOTE]
       << "\t" << STAT_label[STATL_FREQUENCY];
  }
  os << endl;

  for (i = 0;i <= max_value;i++) {
    os << i;

    if (curves) {
      os << "\t";
      if (frequency[i] > 0) {
        os << curves->point[0][i];
      }
    }

    os << "\t" << point[i];

    if (curves) {
      if (frequency[i] > 0) {
        os << "\t" << residual[i] << "\t" << residual[i] / residual_standard_deviation;
      }
      else {
        os << "\t\t";
      }

      os << "\t" << (1. - point[i]) / residual_standard_deviation
         << "\t" << -point[i] / residual_standard_deviation << "\t" << frequency[i];
    }

    os << endl;
  }

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of the self-transition probabilities and the bounds on
 *         the standardized residuals at the Gnuplot format.
 *
 *  \param[in] path                        file path,
 *  \param[in] residual_standard_deviation residual standard deviation.
 *
 *  \return                                error status.
 */
/*--------------------------------------------------------------*/

bool Function::plot_print(const char *path , double residual_standard_deviation) const

{
  bool status = false;
  int i;
  ofstream out_file(path);


  if (out_file) {
    status = true;

    for (i = 0;i <= max_value - min_value;i++) {
      out_file << point[i];
      if (residual_standard_deviation != D_DEFAULT) {
        out_file << " " << (1. - point[i]) / residual_standard_deviation
                 << " " << -point[i] / residual_standard_deviation;
      }
      out_file << endl;
    }
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Default constructor of the NonhomogeneousMarkov class.
 */
/*--------------------------------------------------------------*/

NonhomogeneousMarkov::NonhomogeneousMarkov()

{
  markov_data = NULL;

  homogeneity = NULL;
  self_transition = NULL;

  process = NULL;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Constructor of the NonhomogeneousMarkov class.
 *
 *  \param[in] inb_state number of states,
 *  \param[in] ident     identifiers of the self-transition probability functions.
 */
/*--------------------------------------------------------------*/

NonhomogeneousMarkov::NonhomogeneousMarkov(int inb_state , parametric_function *ident)
:Chain(ORDINARY , inb_state)

{
  int i;


  markov_data = NULL;

  homogeneity = new bool[nb_state];
  self_transition = new Function*[nb_state];

  for (i = 0;i < nb_state;i++) {
    if ((ident[i] == LOGISTIC) || (ident[i] == MONOMOLECULAR)) {
      homogeneity[i] = false;
    }
    else {
      homogeneity[i] = true;
    }
    self_transition[i] = NULL;
  }

  process = new CategoricalSequenceProcess(nb_state , nb_state);
}


/*--------------------------------------------------------------*/
/**
 *  \brief Constructor of the NonhomogeneousMarkov class.
 *
 *  \param[in] pchain           pointer on a Chain object,
 *  \param[in] pself_transition pointer on Function objects,
 *  \param[in] length           sequence length.
 */
/*--------------------------------------------------------------*/

NonhomogeneousMarkov::NonhomogeneousMarkov(const Chain *pchain , const Function **pself_transition ,
                                           int length)
:Chain(*pchain)

{
  int i;


  markov_data = NULL;

  homogeneity = new bool[nb_state];
  self_transition = new Function*[nb_state];

  for (i = 0;i < nb_state;i++) {
    if (pself_transition[i]) {
      homogeneity[i] = false;
      self_transition[i] = new Function(*pself_transition[i]);
    }

    else {
      homogeneity[i] = true;
      self_transition[i] = NULL;
    }
  }

  process = new CategoricalSequenceProcess(nb_state , nb_state);

  characteristic_computation(length , true);
}


/*--------------------------------------------------------------*/
/**
 *  \brief Copy of a NonhomogeneousMarkov object.
 *
 *  \param[in] markov              reference on a NonhomogeneousMarkov object,
 *  \param[in] data_flag           flag copy of the included NonhomogeneousMarkovData object,
 *  \param[in] characteristic_flag flag copy of the characteristic distributions.
 */
/*--------------------------------------------------------------*/

void NonhomogeneousMarkov::copy(const NonhomogeneousMarkov &markov , bool data_flag ,
                                bool characteristic_flag)

{
  int i;


  if ((data_flag) && (markov.markov_data)) {
    markov_data = new NonhomogeneousMarkovData(*(markov.markov_data) , false);
  }
  else {
    markov_data = NULL;
  }

  homogeneity = new bool[nb_state];
  self_transition = new Function*[nb_state];

  for (i = 0;i < nb_state;i++) {
    homogeneity[i] = markov.homogeneity[i];
    if (homogeneity[i]) {
      self_transition[i] = NULL;
    }
    else {
      self_transition[i] = new Function(*(markov.self_transition[i]));
    }
  }

  process = new CategoricalSequenceProcess(*(markov.process) , CATEGORICAL_SEQUENCE_PROCESS_COPY ,
                                           characteristic_flag);
}


/*--------------------------------------------------------------*/
/**
 *  \brief Destruction of the data members of a NonhomogeneousMarkov object.
 */
/*--------------------------------------------------------------*/

void NonhomogeneousMarkov::remove()

{
  int i;


  delete markov_data;

  if (self_transition) {
    for (i = 0;i < nb_state;i++) {
      if (!homogeneity[i]) {
        delete self_transition[i];
      }
    }
    delete [] self_transition;
  }

  delete [] homogeneity;

  delete process;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Destructor of the NonhomogeneousMarkov class.
 */
/*--------------------------------------------------------------*/

NonhomogeneousMarkov::~NonhomogeneousMarkov()

{
  remove();
}


/*--------------------------------------------------------------*/
/**
 *  \brief Assignment operator of the NonhomogeneousMarkov class.
 *
 *  \param[in] markov reference on a NonhomogeneousMarkov object.
 *
 *  \return           NonhomogeneousMarkov object.
 */
/*--------------------------------------------------------------*/

NonhomogeneousMarkov& NonhomogeneousMarkov::operator=(const NonhomogeneousMarkov &markov)

{
  if (&markov != this) {
    remove();
    Chain::remove();

    Chain::copy(markov);
    copy(markov);
  }

  return *this;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Extraction of a distribution.
 *
 *  \param[in] error     reference on a StatError object,
 *  \param[in] dist_type distribution type,
 *  \param[in] state     state.
 *
 *  \return              DiscreteParametricModel object.
 */
/*--------------------------------------------------------------*/

DiscreteParametricModel* NonhomogeneousMarkov::extract(StatError &error , process_distribution dist_type ,
                                                       int state) const

{
  bool status = true;
  Distribution *pdist;
  DiscreteParametric *pparam;
  DiscreteParametricModel *dist;
  FrequencyDistribution *phisto;


  dist = NULL;
  error.init();

  pdist = NULL;
  pparam = NULL;

  if ((state < 0) || (state >= process->nb_value)) {
    status = false;
    ostringstream error_message;
    error_message << STAT_label[STATL_STATE] << " " << state << " "
                  << STAT_error[STATR_NOT_PRESENT];
    error.update((error_message.str()).c_str());
  }

  else {
    switch (dist_type) {
    case FIRST_OCCURRENCE :
      pdist = process->first_occurrence[state];
      break;
    case RECURRENCE_TIME :
      if (process->recurrence_time) {
        pdist = process->recurrence_time[state];
      }
      break;
    case SOJOURN_TIME :
      if (process->sojourn_time) {
        pparam = process->sojourn_time[state];
      }
      break;
    case NB_RUN :
      pdist = process->nb_run[state];
      break;
    case NB_OCCURRENCE :
      pdist = process->nb_occurrence[state];
      break;
    }

    if ((!pdist) && (!pparam)) {
      status = false;
      error.update(SEQ_error[SEQR_NON_EXISTING_CHARACTERISTIC_DISTRIBUTION]);
    }
  }

  if (status) {
    phisto = NULL;

    if (markov_data) {
      switch (dist_type) {

      case FIRST_OCCURRENCE : {
        phisto = markov_data->characteristics[0]->first_occurrence[state];
        break;
      }

      case RECURRENCE_TIME : {
        if (markov_data->characteristics[0]->recurrence_time[state]->nb_element > 0) {
          phisto = markov_data->characteristics[0]->recurrence_time[state];
        }
        break;
      }

      case SOJOURN_TIME : {
        if (markov_data->characteristics[0]->sojourn_time[state]->nb_element > 0) {
          phisto = markov_data->characteristics[0]->sojourn_time[state];
        }
        break;
      }

      case NB_RUN : {
        phisto = markov_data->characteristics[0]->nb_run[state];
        break;
      }

      case NB_OCCURRENCE : {
        phisto = markov_data->characteristics[0]->nb_occurrence[state];
        break;
      }
      }
    }

    if (pdist) {
      dist = new DiscreteParametricModel(*pdist , phisto);
    }
    else if (pparam) {
      dist = new DiscreteParametricModel(*pparam , phisto);
    }
  }

  return dist;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Construction of a NonhomogeneousMarkov object from a file.
 *
 *  \param[in] error  reference on a StatError object,
 *  \param[in] path   file path,
 *  \param[in] length sequence length.
 *
 *  \return           NonhomogeneousMarkov object.
 */
/*--------------------------------------------------------------*/

NonhomogeneousMarkov* NonhomogeneousMarkov::ascii_read(StatError &error , const string path , int length)

{
  string buffer;
  size_t position;
  typedef tokenizer<char_separator<char>> tokenizer;
  char_separator<char> separator(" \t");
  bool status , lstatus;
  int i , j;
  int line , homogeneity , nb_state , index;
  const Chain *chain;
  const Function **self_transition;
  NonhomogeneousMarkov *markov;
  ifstream in_file(path.c_str());


  markov = NULL;
  error.init();

  if (!in_file) {
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {
    status = true;
    line = 0;

    if (length < 2) {
      status = false;
      error.update(SEQ_error[SEQR_SHORT_SEQUENCE_LENGTH]);
    }
    if (length > MAX_LENGTH) {
      status = false;
      error.update(SEQ_error[SEQR_LONG_SEQUENCE_LENGTH]);
    }

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

        // test NONHOMOGENEOUS_MARKOV_CHAIN keyword

        if (i == 0) {
          if (*token != SEQ_word[SEQW_NONHOMOGENEOUS_MARKOV_CHAIN]) {
            status = false;
            error.update(STAT_parsing[STATP_KEYWORD] , line);
          }
        }

        i++;
      }

      if (i > 0) {
        if (i != 1) {
          status = false;
          error.update(STAT_parsing[STATP_FORMAT] , line);
        }
        break;
      }
    }

    // analysis of the format and reading of the Markov chain

    chain = Chain::parsing(error , in_file , line , ORDINARY);

    if (chain) {
      nb_state = chain->nb_state;
      self_transition = new const Function*[nb_state];
      for (i = 0;i < nb_state;i++) {
        self_transition[i] = NULL;
      }

      // analysis of the format of the self-transition probability functions

      for (i = 0;i < nb_state;i++) {
        homogeneity = I_DEFAULT;

        while (getline(in_file , buffer)) {
          line++;

#         ifdef DEBUG
          cout << line << "  " << buffer << endl;
#         endif

          position = buffer.find('#');
          if (position != string::npos) {
            buffer.erase(position);
          }
          j = 0;

          tokenizer tok_buffer(buffer , separator);

          for (tokenizer::iterator token = tok_buffer.begin();token != tok_buffer.end();token++) {
            switch (j) {

            // test STATE keyword

            case 0 : {
              if (*token != STAT_word[STATW_STATE]) {
                status = false;
                error.correction_update(STAT_parsing[STATP_KEYWORD] , STAT_word[STATW_STATE] , line , j + 1);
              }
              break;
            }

            // test state index

            case 1 : {
              lstatus = true;

/*              try {
                index = stoi(*token);   in C++ 11
              }
              catch(invalid_argument &arg) {
                lstatus = false;
              } */
              index = atoi(token->c_str());

              if ((lstatus) && (index != i)) {
                lstatus = false;
              }

              if (!lstatus) {
                status = false;
                error.correction_update(STAT_parsing[STATP_STATE_INDEX] , i , line , j + 1);
              }
              break;
            }

            // test HOMOGENEOUS/NONHOMOGENEOUS keyword

            case 2 : {
              if (*token == SEQ_word[SEQW_HOMOGENEOUS]) {
                homogeneity = true;
              }
              else {
                if (*token == SEQ_word[SEQW_NONHOMOGENEOUS]) {
                  homogeneity = false;
                }
                else {
                  status = false;
                  error.update(STAT_parsing[STATP_KEYWORD] , line , j + 1);
                }
              }
              break;
            }
            }

            j++;
          }

          if (j > 0) {
            if (j != 3) {
              status = false;
              error.update(STAT_parsing[STATP_FORMAT] , line);
            }

            if (!homogeneity) {
              self_transition[i] = Function::parsing(error , in_file , line , length);
              if (!self_transition[i]) {
                status = false;
              }
            }

            break;
          }
        }

        if (homogeneity == I_DEFAULT) {
          status = false;
          error.update(STAT_parsing[STATP_FORMAT] , line);
        }
      }

      while (getline(in_file , buffer)) {
        line++;

#       ifdef DEBUG
        cout << line << "  " << buffer << endl;
#       endif

        position = buffer.find('#');
        if (position != string::npos) {
          buffer.erase(position);
        }
        if (!(trim_right_copy_if(buffer , is_any_of(" \t")).empty())) {
          status = false;
          error.update(STAT_parsing[STATP_FORMAT] , line);
        }
      }

      if (status) {
        markov = new NonhomogeneousMarkov(chain , self_transition , length);
      }

      delete chain;

      for (i = 0;i < nb_state;i++) {
        delete self_transition[i];
      }
      delete [] self_transition;
    }
  }

  return markov;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing on a single line of a NonhomogeneousMarkov object.
 *
 *  \param[in,out] os stream.
 */
/*--------------------------------------------------------------*/

ostream& NonhomogeneousMarkov::line_write(ostream &os) const

{
  os << nb_state << " " << STAT_word[STATW_STATES];

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a NonhomogeneousMarkov object and the associated data structure.
 *
 *  \param[in,out] os         stream,
 *  \param[in]     seq        pointer on a NonhomogeneousMarkovData object,
 *  \param[in]     exhaustive flag detail level,
 *  \param[in]     file_flag  flag file.
 */
/*--------------------------------------------------------------*/

ostream& NonhomogeneousMarkov::ascii_write(ostream &os , const NonhomogeneousMarkovData *seq ,
                                           bool exhaustive , bool file_flag) const

{
  int i;


  os << SEQ_word[SEQW_NONHOMOGENEOUS_MARKOV_CHAIN] << endl;
 
  // writing of the Markov chain parameters

  ascii_print(os , file_flag);

  // writing of the self-transition probability function parameters

  for (i = 0;i < nb_state;i++) {
    os << "\n" << STAT_word[STATW_STATE] << " " << i << " ";

    if (homogeneity[i]) {
      os << SEQ_word[SEQW_HOMOGENEOUS] << endl;
    }
    else {
      os << SEQ_word[SEQW_NONHOMOGENEOUS] << endl;
      self_transition[i]->ascii_print(os , exhaustive , file_flag ,
                                      (seq ? seq->self_transition[i] : NULL));
    }
  }

  process->ascii_print(os , 0 , NULL , NULL , (seq ? seq->characteristics[0] : NULL) ,
                       exhaustive , file_flag);

  if (seq) {
    int nb_parameter = nb_parameter_computation();
    double information , likelihood;


    // writing of the sequence length frequency distribution

    os << "\n";
    if (file_flag) {
      os << "# ";
    }
    os << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";
    seq->length_distribution->ascii_characteristic_print(os , false , file_flag);

    if (exhaustive) {
      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << "   | " << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << endl;
      seq->length_distribution->ascii_print(os , file_flag);
    }

    os << "\n";
    if (file_flag) {
      os << "# ";
    }
    os << SEQ_label[SEQL_CUMUL_LENGTH] << ": " << seq->cumul_length << endl;

    // writing of the information quantity of the sequences in the i.i.d. case

    information = seq->iid_information_computation();

    os << "\n";
    if (file_flag) {
      os << "# ";
    }
    os << SEQ_label[SEQL_IID_INFORMATION] << ": " << information << " ("
       << information / seq->cumul_length << ")" << endl;

    // writing of the (penalized) log-likelihoods of the model for the sequences

    likelihood = seq->likelihood;

    os << "\n";
    if (file_flag) {
      os << "# ";
    }
    os << STAT_label[STATL_LIKELIHOOD] << ": " << likelihood << "   ("
       << STAT_label[STATL_NORMALIZED] << ": " << likelihood / seq->cumul_length << ")" << endl;

    if ((likelihood != D_INF) && (nb_component == 1)) {
      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << nb_parameter << " " << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS]
         << "   2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[AIC] << "): "
         << 2 * (likelihood - nb_parameter) << endl;

      if (nb_parameter < seq->cumul_length - 1) {
        os << "\n";
        if (file_flag) {
          os << "# ";
        }
        os << nb_parameter << " " << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS]
           << "   2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[AICc] << "): "
           << 2 * (likelihood - (double)(nb_parameter * seq->cumul_length) /
             (double)(seq->cumul_length - nb_parameter - 1)) << endl;
      }

      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << nb_parameter << " " << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS]
         << "   2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[BIC] << "): "
         << 2 * likelihood - nb_parameter * log((double)seq->cumul_length) << endl;
    }
  }

# ifdef DEBUG
  MultiPlotSet *plot_set;

  plot_set = get_plotable(seq);
  delete plot_set;
# endif

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a NonhomogeneousMarkov object.
 *
 *  \param[in,out] os         stream,
 *  \param[in]     exhaustive flag detail level.
 */
/*--------------------------------------------------------------*/

ostream& NonhomogeneousMarkov::ascii_write(ostream &os , bool exhaustive) const

{
  return ascii_write(os , markov_data , exhaustive , false);
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a NonhomogeneousMarkov object in a file.
 *
 *  \param[in] error      reference on a StatError object,
 *  \param[in] path       file path,
 *  \param[in] exhaustive flag detail level.
 *
 *  \return               error status.
 */
/*--------------------------------------------------------------*/

bool NonhomogeneousMarkov::ascii_write(StatError &error , const string path ,
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
    ascii_write(out_file , markov_data , exhaustive , true);
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a NonhomogeneousMarkov object and the associated data structure
 *         in a file at the spreadsheet format.
 *
 *  \param[in,out] os  stream,
 *  \param[in]     seq pointer on a NonhomogeneousMarkovData object.
 */
/*--------------------------------------------------------------*/

ostream& NonhomogeneousMarkov::spreadsheet_write(ostream &os , const NonhomogeneousMarkovData *seq) const

{
  int i;


  os << SEQ_word[SEQW_NONHOMOGENEOUS_MARKOV_CHAIN] << endl;

  // writing of the Markov chain parameters

  spreadsheet_print(os);

  // writing of the self-transition probability function parameters

  for (i = 0;i < nb_state;i++) {
    os << "\n" << STAT_word[STATW_STATE] << "\t" << i << "\t";

    if (homogeneity[i]) {
      os << SEQ_word[SEQW_HOMOGENEOUS] << endl;
    }
    else {
      os << SEQ_word[SEQW_NONHOMOGENEOUS] << endl;
      self_transition[i]->spreadsheet_print(os , (seq ? seq->self_transition[i] : NULL));
    }
  }

  process->spreadsheet_print(os , 0 , NULL , NULL , (seq ? seq->characteristics[0] : NULL));

  if (seq) {
    int nb_parameter = nb_parameter_computation();
    double information , likelihood;


    // writing of the sequence length frequency distribution

    os << "\n" << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\t";
    seq->length_distribution->spreadsheet_characteristic_print(os);

    os << "\n\t" << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << endl;
    seq->length_distribution->spreadsheet_print(os);

    os << "\n" << SEQ_label[SEQL_CUMUL_LENGTH] << "\t" << seq->cumul_length << endl;

    // writing of the information quantity of the sequences in the i.i.d. case

    information = seq->iid_information_computation();

    os << "\n" << SEQ_label[SEQL_IID_INFORMATION] << "\t" << information << "\t"
       << information / seq->cumul_length << endl;

    // writing of the (penalized) log-likelihoods of the model for the sequences

    likelihood = seq->likelihood;

    os << "\n" << STAT_label[STATL_LIKELIHOOD] << "\t" << likelihood << "\t"
       << STAT_label[STATL_NORMALIZED] << "\t" << likelihood / seq->cumul_length << endl;

    if ((likelihood != D_INF) && (nb_component == 1)) {
      os << "\n" << nb_parameter << "\t" << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS] << "\t"
         << "2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[AIC] << ")\t"
         << 2 * (likelihood - nb_parameter) << endl;

      if (nb_parameter < seq->cumul_length - 1) {
        os << "\n" << nb_parameter << "\t" << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS] << "\t"
           << "2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[AICc] << ")\t"
           << 2 * (likelihood - (double)(nb_parameter * seq->cumul_length) /
              (double)(seq->cumul_length - nb_parameter - 1)) << endl;
      }

      os << "\n" << nb_parameter << "\t" << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS] << "\t"
         << "2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[BIC] << ")\t"
         << 2 * likelihood - nb_parameter * log((double)seq->cumul_length) << endl;
    }
  }

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a NonhomogeneousMarkov object in a file at the spreadsheet format.
 *
 *  \param[in] error reference on a StatError object,
 *  \param[in] path  file path.
 *
 *  \return          error status.
 */
/*--------------------------------------------------------------*/

bool NonhomogeneousMarkov::spreadsheet_write(StatError &error , const string path) const

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
    spreadsheet_write(out_file , markov_data);
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Plot of a NonhomogeneousMarkov object and the associated data structure
 *         using Gnuplot.
 *
 *  \param[in] prefix file prefix,
 *  \param[in] title  figure title,
 *  \param[in] seq    pointer on a NonhomogeneousMarkovData object.
 *
 *  \return           error status.
 */
/*--------------------------------------------------------------*/

bool NonhomogeneousMarkov::plot_write(const char *prefix , const char *title ,
                                      const NonhomogeneousMarkovData *seq) const

{
  bool status;
  int i , j;
  int variable , start , *pfrequency , max_frequency[NB_STATE];
  double residual_mean , residual_standard_deviation , *standard_residual , *presidual ,
         min_standard_residual[NB_STATE] , max_standard_residual[NB_STATE];
  ostringstream data_file_name[NB_STATE * 2];


  if (seq) {
    status = process->plot_print(prefix , title , 0 , NULL , NULL ,
                                 seq->characteristics[0] , seq->length_distribution);
  }
  else {
    status = process->plot_print(prefix , title , 0);
  }

  if (status) {

    // writing of the data files

    for (i = 0;i < nb_state;i++) {
      if (!homogeneity[i]) {
        if (seq) {
          max_frequency[i] = seq->self_transition[i]->max_frequency_computation();

          // computation of the standardized residuals

          residual_mean = self_transition[i]->residual_mean_computation();
          residual_standard_deviation = sqrt(self_transition[i]->residual_variance_computation(residual_mean));

          standard_residual = new double[self_transition[i]->max_value + 1];

          pfrequency = self_transition[i]->frequency;
          presidual = self_transition[i]->residual;
          min_standard_residual[i] = 0.;
          max_standard_residual[i] = 0.;

          for (j = 0;j <= self_transition[i]->max_value;j++) {
            if (*pfrequency++ > 0) {
              standard_residual[j] = *presidual / residual_standard_deviation;
              if (standard_residual[j] < min_standard_residual[i]) {
                min_standard_residual[i] = standard_residual[j];
              }
              if (standard_residual[j] > max_standard_residual[i]) {
                max_standard_residual[i] = standard_residual[j];
              }
              presidual++;
            }
          }
        }

        data_file_name[i * 2] << prefix << i * 2 << ".dat";
        self_transition[i]->plot_print((data_file_name[i * 2].str()).c_str() ,
                                       (seq ? residual_standard_deviation : D_DEFAULT));

        if (seq) {
          data_file_name[i * 2 + 1] << prefix << i * 2 + 1 << ".dat";
          seq->self_transition[i]->plot_print_standard_residual((data_file_name[i * 2 + 1].str()).c_str() ,
                                                                standard_residual);
          delete [] standard_residual;
        }
      }
    }

    // writing of the script files

    for (i = 0;i < 2;i++) {
      ostringstream file_name[2];

      switch (i) {
      case 0 :
        file_name[0] << prefix << 0 << 0 << ".plot";
        break;
      case 1 :
        file_name[0] << prefix << 0 << 0 << ".print";
        break;
      }

      ofstream out_file((file_name[0].str()).c_str());

      if (i == 1) {
        out_file << "set terminal postscript" << endl;
        file_name[1] << label(prefix) << 0 << 0 << ".ps";
        out_file << "set output \"" << file_name[1].str() << "\"\n\n";
      }

      out_file << "set border 15 lw 0\n" << "set tics out\n" << "set xtics nomirror\n"
               << "set title";
      if (title) {
        out_file << " \"" << title << "\"";
      }
      out_file << "\n\n";

      start = true;
      for (j = 0;j < nb_state;j++) {
        if (!homogeneity[j]) {
          if (!start) {
            if (i == 0) {
              out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
            }
            out_file << endl;
          }
          else {
            start = false;
          }

          if (self_transition[j]->max_value < TIC_THRESHOLD) {
            out_file << "set xtics 0,1" << endl;
          }

          out_file << "plot [0:" << self_transition[j]->max_value << "] [0:1] ";
          if (seq) {
            out_file << "\""<< label((data_file_name[j * 2 + 1].str()).c_str())
                     << "\" using 1:2 title \"" << STAT_label[STATL_STATE] << " " << j << " - "
                     << SEQ_label[SEQL_OBSERVED] << " " << SEQ_label[SEQL_SELF_TRANSITION]
                     << "\" with points,\\" << endl;
          }
          out_file << "\"" << label((data_file_name[j * 2].str()).c_str())
                   << "\" using 1 title \"" << STAT_label[STATL_STATE] << " " << j << " - "
                   << SEQ_label[SEQL_THEORETICAL] << " " << SEQ_label[SEQL_SELF_TRANSITION]
                   << "\" with linespoints" << endl;

          if (self_transition[j]->max_value < TIC_THRESHOLD) {
            out_file << "set xtics autofreq" << endl;
          }

          if (seq) {
            if (i == 0) {
              out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
            }
            out_file << endl;

            out_file << "set xlabel \"" << SEQ_label[SEQL_INDEX] << "\"" << endl;
            out_file << "set ylabel \"" << STAT_label[STATL_STANDARDIZED_RESIDUAL] << "\"" << endl;

            if (seq->self_transition[j]->length - 1 < TIC_THRESHOLD) {
              out_file << "set xtics 0,1" << endl;
            }

            out_file << "plot [0:" << seq->self_transition[j]->length - 1 << "] ["
                     << min_standard_residual[j] << ":" << max_standard_residual[j] << "] \""
                     << label((data_file_name[j * 2 + 1].str()).c_str())
                     << "\" using 1:3 notitle with points";
            if (((1. - self_transition[j]->point[0]) / residual_standard_deviation <= max_standard_residual[j]) ||
                ((1. - self_transition[j]->point[seq->self_transition[j]->length - 1]) / residual_standard_deviation <=
                 max_standard_residual[j])) {
              out_file << ",\\\n\"" << label((data_file_name[j * 2].str()).c_str())
                       << "\" using 2 title \"" << SEQ_label[SEQL_ASYMPTOTE] << "\" with lines";
            }
            if ((-self_transition[j]->point[0] / residual_standard_deviation >= min_standard_residual[j]) ||
                (-self_transition[j]->point[seq->self_transition[j]->length - 1] / residual_standard_deviation >=
                 min_standard_residual[j])) {
              out_file << ",\\\n\"" << label((data_file_name[j * 2].str()).c_str())
                       << "\" using 3 title \"" << SEQ_label[SEQL_ASYMPTOTE] << "\" with lines";
            }
            out_file << endl;

            out_file << "set xlabel" << endl;
            out_file << "set ylabel" << endl;

            if (seq->self_transition[j]->length - 1 < TIC_THRESHOLD) {
              out_file << "set xtics autofreq" << endl;
            }

            if (i == 0) {
              out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
            }
            out_file << endl;

            out_file << "set xlabel \"" << SEQ_label[SEQL_INDEX] << "\"" << endl;
            out_file << "set ylabel \"" << STAT_label[STATL_FREQUENCY] << "\"" << endl;

            if (seq->self_transition[j]->length - 1 < TIC_THRESHOLD) {
              out_file << "set xtics 0,1" << endl;
            }
            if ((int)(max_frequency[j] * YSCALE) + 1 < TIC_THRESHOLD) {
              out_file << "set ytics 0,1" << endl;
            }

            out_file << "plot [0:" << seq->self_transition[j]->length - 1
                     << "] [0:" << (int)(max_frequency[j] * YSCALE) + 1 << "] \""
                     << label((data_file_name[j * 2 + 1].str()).c_str())
                     << "\" using 1:4 title \"" << STAT_label[STATL_STATE] << " "
                     << j << " - "<< SEQ_label[SEQL_TRANSITION_COUNTS]
                     << "\" with impulses" << endl;

            out_file << "set xlabel" << endl;
            out_file << "set ylabel" << endl;

            if (seq->self_transition[j]->length - 1 < TIC_THRESHOLD) {
              out_file << "set xtics autofreq" << endl;
            }
            if (max_frequency[j] < TIC_THRESHOLD) {
              out_file << "set ytics autofreq" << endl;
            }
          }
        }
      }

      if (i == 1) {
        out_file << "\nset terminal x11" << endl;
      }

      out_file << "\npause 0 \"" << STAT_label[STATL_END] << "\"" << endl;
    }
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Plot of a NonhomogeneousMarkov object using Gnuplot.
 *
 *  \param[in] error  reference on a StatError object,
 *  \param[in] prefix file prefix,
 *  \param[in] title  figure title.
 *
 *  \return           error status.
 */
/*--------------------------------------------------------------*/

bool NonhomogeneousMarkov::plot_write(StatError &error , const char *prefix ,
                                      const char *title) const

{
  bool status = plot_write(prefix , title , markov_data);

  error.init();

  if (!status) {
    error.update(STAT_error[STATR_FILE_PREFIX]);
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Plot of a NonhomogeneousMarkov object and the associated data structure.
 *
 *  \param[in] seq pointer on a NonhomogeneousMarkovData object .
 *
 *  \return        MultiPlotSet object.
 */
/*--------------------------------------------------------------*/

MultiPlotSet* NonhomogeneousMarkov::get_plotable(const NonhomogeneousMarkovData *seq) const

{
  int i , j , k;
  int nb_plot_set , index_length , index , nb_plot , max_frequency , *pfrequency;
  double residual_mean , residual_standard_deviation , min_standard_residual ,
         max_standard_residual , *standard_residual , *presidual;
  ostringstream legend;
  FrequencyDistribution *length_distribution;
  SequenceCharacteristics *characteristics;
  MultiPlotSet *plot_set;


  if (seq) {
    characteristics = seq->characteristics[0];
    length_distribution = seq->length_distribution;
  }
  else {
    characteristics = NULL;
    length_distribution = NULL;
  }

  // computation of the number of plots

  nb_plot_set = 0;

  if ((process->index_value) || (characteristics)) {
    nb_plot_set++;

    if (characteristics) {
      index_length = characteristics->index_value->plot_length_computation();

      if (characteristics->index_value->frequency[index_length - 1] < MAX_FREQUENCY) {
        nb_plot_set++;
      }
      nb_plot_set++;
    }
  }

  if ((process->first_occurrence) || (characteristics)) {
    for (i = 0;i < nb_state;i++) {
      if ((process->first_occurrence) && (process->first_occurrence[i])) {
        nb_plot_set++;
      }
      else if ((characteristics) && (i < characteristics->nb_value) &&
               (characteristics->first_occurrence[i]->nb_element > 0)) {
        nb_plot_set++;
      }
    }
  }

  if ((process->recurrence_time) || (characteristics)) {
    for (i = 0;i < nb_state;i++) {
      if ((process->recurrence_time) && (process->recurrence_time[i])) {
        nb_plot_set++;
      }
      else if ((characteristics) && (i < characteristics->nb_value) &&
               (characteristics->recurrence_time[i]->nb_element > 0)) {
        nb_plot_set++;
      }
    }
  }

  if ((process->sojourn_time) || (characteristics)) {
    for (i = 0;i < nb_state;i++) {
      if ((process->sojourn_time) && (process->sojourn_time[i])) {
        nb_plot_set++;
      }
      else if ((characteristics) && (i < characteristics->nb_value) &&
               (characteristics->sojourn_time[i]->nb_element > 0)) {
        nb_plot_set++;
      }

      if ((characteristics) && (i < characteristics->nb_value) &&
          (characteristics->initial_run) &&
          (characteristics->initial_run[i]->nb_element > 0)) {
        nb_plot_set++;
      }

      if ((characteristics) && (i < characteristics->nb_value) &&
          (characteristics->final_run[i]->nb_element > 0)) {
        nb_plot_set++;
      }
    }
  }

  if ((process->nb_run) || (process->nb_occurrence) ||
      ((characteristics) && (characteristics->nb_run) && (characteristics->nb_occurrence))) {
    for (i = 0;i < nb_state;i++) {
      if (process->nb_run) {
        nb_plot_set++;
      }
      else if ((characteristics) && (i < characteristics->nb_value) &&
               (characteristics->nb_run) && (characteristics->nb_run[i]->nb_element > 0)) {
        nb_plot_set++;
      }

      if (process->nb_occurrence) {
        nb_plot_set++;
      }
      else if ((characteristics) && (i < characteristics->nb_value) &&
               (characteristics->nb_occurrence) &&
               (characteristics->nb_occurrence[i]->nb_element > 0)) {
        nb_plot_set++;
      }
    }

    if ((characteristics) && (characteristics->nb_run) && (characteristics->nb_occurrence)) {
      nb_plot_set++;
    }
  }

  for (i = 0;i < nb_state;i++) {
    if (!homogeneity[i]) {
      nb_plot_set++;
      if (seq) {
        nb_plot_set += 2;
      }
    }
  }

  plot_set = new MultiPlotSet(nb_plot_set , 1);

  MultiPlotSet &plot = *plot_set;

  plot.border = "15 lw 0";

  plot.variable_nb_viewpoint[0] = 1;

  index = 0;
  for (i = 0;i < nb_state;i++) {
    if (!homogeneity[i]) {

      // self-transition in state i probability function

      plot.variable[index] = 0;
      plot.viewpoint[index] = SELF_TRANSITION;

      plot[index].xrange = Range(0 , self_transition[i]->max_value);
      plot[index].yrange = Range(0. , 1.);

      if (self_transition[i]->max_value < TIC_THRESHOLD) {
        plot[index].xtics = 1;
      }

      if (seq) {
        plot[index].resize(2);

        legend.str("");
        legend << STAT_label[STATL_STATE] << " " << i << " - "
               << SEQ_label[SEQL_OBSERVED] << " " << SEQ_label[SEQL_SELF_TRANSITION];
        plot[index][0].legend = legend.str();

        plot[index][0].style = "points";

        seq->self_transition[i]->plotable_write(0 , plot[index][0]);
        j = 1;
      }

      else {
        plot[index].resize(1);
        j = 0;
      }

      legend.str("");
      legend << STAT_label[STATL_STATE] << " " << i << " - "
             << SEQ_label[SEQL_THEORETICAL] << " " << SEQ_label[SEQL_SELF_TRANSITION];
      plot[index][j].legend = legend.str();

      plot[index][j].style = "linespoint";

      self_transition[i]->plotable_write(plot[index][j]);
      index++;

      if (seq) {

        // computation of the standardized residuals

        residual_mean = self_transition[i]->residual_mean_computation();
        residual_standard_deviation = sqrt(self_transition[i]->residual_variance_computation(residual_mean));

        standard_residual = new double[self_transition[i]->max_value + 1];

        pfrequency = self_transition[i]->frequency;
        presidual = self_transition[i]->residual;
        min_standard_residual = 0.;
        max_standard_residual = 0.;

        for (j = 0;j <= self_transition[i]->max_value;j++) {
          if (*pfrequency++ > 0) {
            standard_residual[j] = *presidual / residual_standard_deviation;
            if (standard_residual[j] < min_standard_residual) {
              min_standard_residual = standard_residual[j];
            }
            if (standard_residual[j] > max_standard_residual) {
              max_standard_residual = standard_residual[j];
            }
            presidual++;
          }
        }

        // standardized residuals

        plot.variable[index] = 0;
        plot.viewpoint[index] = SELF_TRANSITION;

        plot[index].xrange = Range(0 , seq->self_transition[i]->length - 1);
        if (seq->self_transition[i]->length - 1 < TIC_THRESHOLD) {
          plot[index].xtics = 1;
        }

        plot[index].yrange = Range(min_standard_residual , max_standard_residual);

        plot[index].xlabel = SEQ_label[SEQL_INDEX];
        plot[index].ylabel = STAT_label[STATL_STANDARDIZED_RESIDUAL];

        nb_plot = 1;
        if (((1. - self_transition[i]->point[0]) / residual_standard_deviation <= max_standard_residual) ||
            ((1. - self_transition[i]->point[seq->self_transition[i]->length - 1]) / residual_standard_deviation <=
             max_standard_residual)) {
          nb_plot++;
        }
        if ((-self_transition[i]->point[0] / residual_standard_deviation >= min_standard_residual) ||
            (-self_transition[i]->point[seq->self_transition[i]->length - 1] / residual_standard_deviation >=
             min_standard_residual)) {
          nb_plot++;
        }
        plot[index].resize(nb_plot);

        plot[index][0].style = "points";

        pfrequency = self_transition[i]->frequency;
        for (j = 0;j <= self_transition[i]->max_value;j++) {
          if (*pfrequency++ > 0) {
            plot[index][0].add_point(j , standard_residual[j]);
          }
        }

        j = 1;
        if (((1. - self_transition[i]->point[0]) / residual_standard_deviation <= max_standard_residual) ||
            ((1. - self_transition[i]->point[seq->self_transition[i]->length - 1]) / residual_standard_deviation <=
             max_standard_residual)) {
          plot[index][j].legend = SEQ_label[SEQL_ASYMPTOTE];

          plot[index][j].style = "lines";

          for (k = 0;k <= self_transition[i]->max_value;k++) {
            plot[index][j].add_point(k , (1. - self_transition[i]->point[k]) / residual_standard_deviation);
          }
          j++;
        }

        if ((-self_transition[i]->point[0] / residual_standard_deviation >= min_standard_residual) ||
            (-self_transition[i]->point[seq->self_transition[i]->length - 1] / residual_standard_deviation >=
             min_standard_residual)) {
          plot[index][j].legend = SEQ_label[SEQL_ASYMPTOTE];

          plot[index][j].style = "lines";

          for (k = 0;k <= self_transition[i]->max_value;k++) {
            plot[index][j].add_point(k , -self_transition[i]->point[k] / residual_standard_deviation);
          }
        }
        index++;

        // self-transition in state i empirical function

        plot.variable[index] = 0;
        plot.viewpoint[index] = SELF_TRANSITION;

        plot[index].xrange = Range(0 , seq->self_transition[i]->length - 1);
        max_frequency = seq->self_transition[i]->max_frequency_computation();
        plot[index].yrange = Range(0 , ceil(max_frequency * YSCALE));

        if (seq->self_transition[i]->length - 1 < TIC_THRESHOLD) {
          plot[index].xtics = 1;
        }
        if (ceil(max_frequency * YSCALE) < TIC_THRESHOLD) {
          plot[index].ytics = 1;
        }

        plot[index].xlabel = SEQ_label[SEQL_INDEX];
        plot[index].ylabel = STAT_label[STATL_FREQUENCY];

        plot[index].resize(1);

        legend.str("");
        legend << STAT_label[STATL_STATE] << " " << i << " - "
               << SEQ_label[SEQL_TRANSITION_COUNTS];
        plot[index][0].legend = legend.str();

        plot[index][0].style = "impulses";

        seq->self_transition[i]->plotable_frequency_write(plot[index][0]);
        index++;

        delete [] standard_residual;
      }
    }
  }

  process->plotable_write(*plot_set , index , 0 , NULL , NULL ,
                          characteristics , length_distribution);

  return plot_set;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Plot of a NonhomogeneousMarkov object.
 */
/*--------------------------------------------------------------*/

MultiPlotSet* NonhomogeneousMarkov::get_plotable() const

{
  return get_plotable(markov_data);
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the number of parameters of a NonhomogeneousMarkov object.
 *
 *  \return number of parameters.
 */
/*--------------------------------------------------------------*/

int NonhomogeneousMarkov::nb_parameter_computation() const

{
  int i;
  int nb_parameter = Chain::nb_parameter_computation();


  for (i = 0;i < nb_state;i++) {
    if (!homogeneity[i]) {
      nb_parameter += self_transition[i]->nb_parameter - 1;
    }
  }

  return nb_parameter;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Default constructor of the NonhomogeneousMarkovData class.
 */
/*--------------------------------------------------------------*/

NonhomogeneousMarkovData::NonhomogeneousMarkovData()

{
  markov = NULL;
  chain_data = NULL;
  likelihood = D_INF;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Constructor of the NonhomogeneousMarkovData class.
 *
 *  \param[in] ilength_distribution sequence length frequency distribution.
 */
/*--------------------------------------------------------------*/

NonhomogeneousMarkovData::NonhomogeneousMarkovData(const FrequencyDistribution &ilength_distribution)
:MarkovianSequences(ilength_distribution , 1 , NULL , false)

{
  markov = NULL;
  chain_data = NULL;
  likelihood = D_INF;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Construction of a NonhomogeneousMarkovData object from
 *         a MarkovianSequences object.
 *
 *  \param[in] seq reference on a MarkovianSequences object.
 */
/*--------------------------------------------------------------*/

NonhomogeneousMarkovData::NonhomogeneousMarkovData(const MarkovianSequences &seq)
:MarkovianSequences(seq)

{
  markov = NULL;
  chain_data = NULL;
  likelihood = D_INF;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Copy of a NonhomogeneousMarkovData object.
 *
 *  \param[in] seq        reference on a NonhomogeneousMarkovData object,
 *  \param[in] model_flag flag copy of the included NonhomogeneousMarkov object.
 */
/*--------------------------------------------------------------*/

void NonhomogeneousMarkovData::copy(const NonhomogeneousMarkovData &seq , bool model_flag)

{
  if ((model_flag) && (seq.markov)) {
    markov = new NonhomogeneousMarkov(*(seq.markov) , false);
  }
  else {
    markov = NULL;
  }

  if (seq.chain_data) {
    chain_data = new ChainData(*(seq.chain_data));
  }
  else {
    chain_data = NULL;
  }

  likelihood = seq.likelihood;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Destructor of the NonhomogeneousMarkovData class.
 */
/*--------------------------------------------------------------*/

NonhomogeneousMarkovData::~NonhomogeneousMarkovData()

{
  delete markov;
  delete chain_data;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Assignment operator of the NonhomogeneousMarkovData class.
 *
 *  \param[in] seq reference on a NonhomogeneousMarkovData object.
 *
 *  \return         NonhomogeneousMarkovData object.
 */
/*--------------------------------------------------------------*/

NonhomogeneousMarkovData& NonhomogeneousMarkovData::operator=(const NonhomogeneousMarkovData &seq)

{
  if (&seq != this) {
    delete markov;
    delete chain_data;

    remove();
    Sequences::remove();

    Sequences::copy(seq);
    MarkovianSequences::copy(seq);
    copy(seq);
  }

  return *this;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Extraction of a frequency distribution.
 *
 *  \param[in] error      reference on a StatError object,
 *  \param[in] histo_type frequency distribution type,
 *  \param[in] state      state.
 *
 *  \return               DiscreteDistributionData object.
 */
/*--------------------------------------------------------------*/

DiscreteDistributionData* NonhomogeneousMarkovData::extract(StatError &error , process_distribution histo_type ,
                                                            int state) const

{
  bool status = true;
  Distribution *pdist;
  DiscreteParametric *pparam;
  FrequencyDistribution *phisto;
  DiscreteDistributionData *histo;


  histo = NULL;
  error.init();

  if ((state < 0) || (state >= marginal_distribution[0]->nb_value)) {
    status = false;
    ostringstream error_message;
    error_message << STAT_label[STATL_STATE] << " " << state << " "
                  << STAT_error[STATR_NOT_PRESENT];
    error.update((error_message.str()).c_str());
  }

  else {
    switch (histo_type) {
    case FIRST_OCCURRENCE :
      phisto = characteristics[0]->first_occurrence[state];
      break;
    case RECURRENCE_TIME :
      phisto = characteristics[0]->recurrence_time[state];
      break;
    case SOJOURN_TIME :
      phisto = characteristics[0]->sojourn_time[state];
      break;
    case FINAL_RUN :
      phisto = characteristics[0]->final_run[state];
      break;
    case NB_RUN :
      phisto = characteristics[0]->nb_run[state];
      break;
    case NB_OCCURRENCE :
      phisto = characteristics[0]->nb_occurrence[state];
      break;
    }

    if (phisto->nb_element == 0) {
      status = false;
      error.update(STAT_error[STATR_EMPTY_SAMPLE]);
    }
  }

  if (status) {
    pdist = NULL;
    pparam = NULL;

    switch (histo_type) {
    case FIRST_OCCURRENCE :
      pdist = markov->process->first_occurrence[state];
      break;
    case RECURRENCE_TIME :
      if (markov->process->recurrence_time) {
        pdist = markov->process->recurrence_time[state];
      }
      break;
    case SOJOURN_TIME :
      if (markov->process->sojourn_time) {
        pparam = markov->process->sojourn_time[state];
      }
      break;
    case NB_RUN :
      pdist = markov->process->nb_run[state];
      break;
    case NB_OCCURRENCE :
      pdist = markov->process->nb_occurrence[state];
      break;
    }

    if (pdist) {
      histo = new DiscreteDistributionData(*phisto , pdist);
    }
    else {
      histo = new DiscreteDistributionData(*phisto , pparam);
    }
  }

  return histo;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Copy of a NonhomogeneousMarkovData object transforming the implicit index parameters in
 *         explicit index parameters.
 *
 *  \param[in] error reference on a StatError object.
 *
 *  \return          NonhomogeneousMarkovData object.
 */
/*--------------------------------------------------------------*/

NonhomogeneousMarkovData* NonhomogeneousMarkovData::explicit_index_parameter(StatError &error) const

{
  NonhomogeneousMarkovData *seq;


  error.init();

  if (index_parameter) {
    seq = NULL;
    error.update(SEQ_error[SEQR_INDEX_PARAMETER_TYPE]);
  }
  else {
    seq = new NonhomogeneousMarkovData(*this , true , EXPLICIT_INDEX_PARAMETER);
  }

  return seq;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Removing of the index parameters.
 *
 *  \param[in] error reference on a StatError object.
 *
 *  \return          NonhomogeneousMarkovData object.
 */
/*--------------------------------------------------------------*/

NonhomogeneousMarkovData* NonhomogeneousMarkovData::remove_index_parameter(StatError &error) const

{
  NonhomogeneousMarkovData *seq;


  error.init();

  if (!index_parameter) {
    seq = NULL;
    error.update(SEQ_error[SEQR_INDEX_PARAMETER_TYPE]);
  }
  else {
    seq = new NonhomogeneousMarkovData(*this , true , REMOVE_INDEX_PARAMETER);
  }

  return seq;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a NonhomogeneousMarkovData object.
 *
 *  \param[in,out] os         stream,
 *  \param[in]     exhaustive flag detail level.
 */
/*--------------------------------------------------------------*/

ostream& NonhomogeneousMarkovData::ascii_write(ostream &os , bool exhaustive) const

{
  if (markov) {
    markov->ascii_write(os , this , exhaustive , false);
  }

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a NonhomogeneousMarkovData object in a file.
 *
 *  \param[in] error      reference on a StatError object,
 *  \param[in] path       file path,
 *  \param[in] exhaustive flag detail level.
 *
 *  \return               error status.
 */
/*--------------------------------------------------------------*/

bool NonhomogeneousMarkovData::ascii_write(StatError &error , const string path ,
                                           bool exhaustive) const

{
  bool status = false;


  if (markov) {
    ofstream out_file(path.c_str());

    error.init();

    if (!out_file) {
      status = false;
      error.update(STAT_error[STATR_FILE_NAME]);
    }

    else {
      status = true;
      markov->ascii_write(out_file , this , exhaustive , true);
    }
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a NonhomogeneousMarkovData object in a file at the spreadsheet format.
 *
 *  \param[in] error reference on a StatError object,
 *  \param[in] path  file path.
 *
 *  \return          error status.
 */
/*--------------------------------------------------------------*/

bool NonhomogeneousMarkovData::spreadsheet_write(StatError &error , const string path) const

{
  bool status = false;


  if (markov) {
    ofstream out_file(path.c_str());

    error.init();

    if (!out_file) {
      status = false;
      error.update(STAT_error[STATR_FILE_NAME]);
    }

    else {
      status = true;
      markov->spreadsheet_write(out_file , this);
    }
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Plot of a NonhomogeneousMarkovData object using Gnuplot.
 *
 *  \param[in] error  reference on a StatError object,
 *  \param[in] prefix file prefix,
 *  \param[in] title  figure title.
 *
 *  \return           error status.
 */
/*--------------------------------------------------------------*/

bool NonhomogeneousMarkovData::plot_write(StatError &error , const char *prefix ,
                                          const char *title) const

{
  bool status = false;


  if (markov) {
    status = markov->plot_write(prefix , title , this);

    error.init();

    if (!status) {
      error.update(STAT_error[STATR_FILE_PREFIX]);
    }
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Plot of a NonhomogeneousMarkovData object.
 *
 *  \return MultiPlotSet object.
 */
/*--------------------------------------------------------------*/

MultiPlotSet* NonhomogeneousMarkovData::get_plotable() const

{
  MultiPlotSet *plot_set;


  if (markov) {
    plot_set = markov->get_plotable(this);
  }
  else {
    plot_set = NULL;
  }

  return plot_set;
}


};  // namespace sequence_analysis
