/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       StructureAnalysis: Identifying patterns in plant architecture and development
 *
 *       Copyright 1995-2018 CIRAD AGAP
 *
 *       File author(s): Yann Guedon (yann.guedon@cirad.fr)
 *
 *       $Source$
 *       $Id$
 *
 *       Forum for StructureAnalysis developers:
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



#include <cmath>

#include <string>
#include <sstream>
#include <iomanip>
#include <fstream>

#include <boost/math/distributions/normal.hpp>

#include "regression.h"
#include "stat_label.h"

using namespace std;
using namespace boost::math;


namespace stat_tool {



/*--------------------------------------------------------------*/
/**
 *  \brief Default constructor of the RegressionKernel class.
 */
/*--------------------------------------------------------------*/

RegressionKernel::RegressionKernel()

{
  ident = NONPARAMETRIC_FUNCTION;

  min_value = 0;
  max_value = 0;

  regression_df = 0.;
  residual_df = 0.;

  nb_parameter = 0;
  parameter = NULL;
//  step = 1.;
  point = NULL;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Constructor of the RegressionKernel class.
 *
 *  \param[in] iident     identifier,
 *  \param[in] imin_value minimum value of the explanatory variable,
 *  \param[in] imax_value maximum value of the explanatory variable.
 */
/*--------------------------------------------------------------*/

// RegressionKernel::RegressionKernel(parametric_function iident , double imin_value ,
//                                    double imax_value , double istep)
RegressionKernel::RegressionKernel(parametric_function iident , int imin_value , int imax_value)

{
  ident = iident;

  if (ident == LINEAR_FUNCTION) {
    min_value = MIN(imin_value , 0);
  }
  else {
    min_value = imin_value;
  }

  max_value = imax_value;

  regression_df = 0.;
  residual_df = 0.;

  switch (ident) {
  case LINEAR_FUNCTION :
    nb_parameter = 2;
    break;
  case LOGISTIC :
    nb_parameter = 3;
    break;
  case MONOMOLECULAR :
    nb_parameter = 3;
    break;
  case NONPARAMETRIC_FUNCTION :
    nb_parameter = 0;
    break;
  }

  if (nb_parameter > 0) {
    parameter = new double[nb_parameter];
  }
  else {
    parameter = NULL;
  }

//  step = istep;

//  point = new double[(int)((max_value - min_value) * step) + 1];
  point = new double[max_value - min_value + 1];
}


/*--------------------------------------------------------------*/
/**
 *  \brief Copy of a RegressionKernel object.
 *
 *  \param[in] regression reference on a RegressionKernel object.
 */
/*--------------------------------------------------------------*/

void RegressionKernel::copy(const RegressionKernel &regression)

{
  int i;


  ident = regression.ident;

  min_value = regression.min_value;
  max_value = regression.max_value;

  regression_df = regression.regression_df;
  residual_df = regression.residual_df;

  nb_parameter = regression.nb_parameter;

  if (regression.parameter) {
    parameter = new double[nb_parameter];
    for (i = 0;i < nb_parameter;i++) {
      parameter[i] = regression.parameter[i];
    }
  }
  else {
    parameter = NULL;
  }

//  step = regression.step;

//  point = new double[(int)((max_value - min_value) * step) + 1];
//  for (i = 0;i <= (int)((max_value - min_value) * step);i++) {
  point = new double[max_value - min_value + 1];
  for (i = 0;i <= max_value - min_value;i++) {
    point[i] = regression.point[i];
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Destruction of the data members of a RegressionKernel object.
 */
/*--------------------------------------------------------------*/

void RegressionKernel::remove()

{
  delete [] parameter;
  delete [] point;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Destructor of the RegressionKernel class.
 */
/*--------------------------------------------------------------*/

RegressionKernel::~RegressionKernel()

{
  remove();
}


/*--------------------------------------------------------------*/
/**
 *  \brief Assignment operator of the RegressionKernel class.
 *
 *  \param[in] regression reference on a RegressionKernel object.
 *
 *  \return               RegressionKernel object.
 */
/*--------------------------------------------------------------*/

RegressionKernel& RegressionKernel::operator=(const RegressionKernel &regression)

{
  if (&regression != this) {
    remove();
    copy(regression);
  }

  return *this;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of the regression function parameters.
 *
 *  \param[in,out] os stream.
 */
/*--------------------------------------------------------------*/

ostream& RegressionKernel::ascii_parameter_print(ostream &os) const

{
  if (ident == LINEAR_FUNCTION) {
    os << STAT_label[STATL_INTERCEPT] << ": " << parameter[0] << "   "
       << STAT_label[STATL_SLOPE] << ": " << parameter[1] << "   ";
  }

  else if ((ident == LOGISTIC) || (ident == MONOMOLECULAR)) {
    int i;


    os << STAT_function_word[ident] << " " << STAT_word[STATW_FUNCTION];
    for (i = 0;i < nb_parameter;i++) {
      os << "   " << STAT_word[STATW_PARAMETER] << " " << i + 1 << " : " << parameter[i];
    }
  }

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Formal writing of the regression function.
 *
 *  \param[in,out] os stream.
 */
/*--------------------------------------------------------------*/

ostream& RegressionKernel::ascii_formal_print(ostream &os) const

{
  switch (ident) {

  case LINEAR_FUNCTION : {
    os << "Y = " << parameter[0];
    if (parameter[1] < 0.) {
      os << " - " << -parameter[1] << " X";
    }
    else {
      os << " + " << parameter[1] << " X";
    }
    break;
  }

  case LOGISTIC : {
    os << "Y = " << parameter[0] << " / (1";
    if (parameter[1] < 0.) {
      os << " - " << -parameter[1];
    }
    else {
      os << " + " << parameter[1];
    }
    os << " EXP(-" << parameter[2] << " X))";
    break;
  }

  case MONOMOLECULAR : {
    os << "Y = " << parameter[0];
    if (parameter[1] < 0.) {
      os << " - " << -parameter[1];
    }
    else {
      os << " + " << parameter[1];
    }
    os << " EXP(-" << parameter[2] << " X)";
    break;
  }
  }

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of the regression function.
 *
 *  \param[in,out] os stream.
 */
/*--------------------------------------------------------------*/

ostream& RegressionKernel::ascii_print(ostream &os) const

{
  int i;
  int buff , width[2];
  double *ppoint;


  width[0] = column_width(min_value , max_value);
  width[1] = column_width(max_value - min_value + 1 , point) + ASCII_SPACE;

  os << "\n" << STAT_label[STATL_EXPLANATORY_VARIABLE] << " | " << STAT_label[STATL_RESPONSE_VARIABLE] << endl;
  ppoint = point;
  for (i = min_value;i <= max_value;i++) {
    os << setw(width[0]) << i;
    os << setw(width[1]) << *ppoint++ << endl;
  }

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of the regression function at the spreadsheet format.
 *
 *  \param[in,out] os stream.
 */
/*--------------------------------------------------------------*/

ostream& RegressionKernel::spreadsheet_print(ostream &os) const

{
  int i;
  double *ppoint;


  os << "\n" << STAT_label[STATL_EXPLANATORY_VARIABLE] << "\t" << STAT_label[STATL_RESPONSE_VARIABLE] << endl;
  ppoint = point;
  for (i = min_value;i <= max_value;i++) {
    os << i << "\t" << *ppoint++ << endl;
  }

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of the regression function (Gnuplot output).
 *
 *  \param[in] path file path.
 *
 *  \return         error status.
 */
/*--------------------------------------------------------------*/

bool RegressionKernel::plot_print(const char *path) const

{
  bool status = false;
  int i;
  double *ppoint;
  ofstream out_file(path);


  if (out_file) {
    status = true;

    if (ident == LINEAR_FUNCTION) {
      out_file << min_value << " " << point[0] << endl;
      out_file << max_value << " " << point[max_value - min_value] << endl;
    }

    else {
      ppoint = point;
      for (i = min_value;i <= max_value;i++) {
        out_file << i << " " << *ppoint++ << endl;
      }
    }
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of the regression function for plotting.
 *
 *  \param[in] plot reference on a SinglePlot object.
 */
/*--------------------------------------------------------------*/

void RegressionKernel::plotable_write(SinglePlot &plot) const

{
  if (ident == LINEAR_FUNCTION) {
    plot.add_point(min_value , point[0]);
    plot.add_point(max_value , point[max_value - min_value]);
  }

  else {
    int i;
    double *ppoint;

    ppoint = point;
    for (i = min_value;i <= max_value;i++) {
      plot.add_point(i , *ppoint++);
    }
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of a parametric regression function.
 */
/*--------------------------------------------------------------*/

void RegressionKernel::computation()

{
  int i;
  double *ppoint;


  ppoint = point;

  switch (ident) {

  case LINEAR_FUNCTION : {
    for (i = min_value;i <= max_value;i++) {
      *ppoint++ = parameter[0] + parameter[1] * i;
    }
    break;
  }

  case LOGISTIC : {
    for (i = min_value;i <= max_value;i++) {
      *ppoint++ = parameter[0] / (1. + parameter[1] * exp(-parameter[2] * i));
    }
    break;
  }

  case MONOMOLECULAR : {
    for (i = min_value;i <= max_value;i++) {
      *ppoint++ = parameter[0] + parameter[1] * exp(-parameter[2] * i);
    }
    break;
  }
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the minimum value of a regression function.
 *
 *  \return minimum value.
 */
/*--------------------------------------------------------------*/

double RegressionKernel::min_computation() const

{
  int i;
  double min , *ppoint;


  ppoint = point;
  min = *ppoint++;

  for (i = (int)min_value + 1;i <= (int)max_value;i++) {
    if (*ppoint < min) {
      min = *ppoint;
    }
    ppoint++;
  }

  return min;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the maximum value of a regression function.
 *
 *  \return maximum value.
 */
/*--------------------------------------------------------------*/

double RegressionKernel::max_computation() const

{
  int i;
  double max , *ppoint;


  ppoint = point;
  max = *ppoint++;

  for (i = (int)min_value + 1;i <= (int)max_value;i++) {
    if (*ppoint > max) {
      max = *ppoint;
    }
    ppoint++;
  }

  return max;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Default constructor of the Regression class.
 */
/*--------------------------------------------------------------*/

Regression::Regression()

{
  vectors = NULL;

  nb_vector = 0;
  residual = NULL;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Construction of a Regression object from a Vectors object.
 *
 *  \param[in] iident               regression function identifier,
 *  \param[in] explanatory_variable explanatory variable index,
 *  \param[in] response_variable    response variable index,
 *  \param[in] vec                  reference on a Vectors object.
 */
/*--------------------------------------------------------------*/

Regression::Regression(parametric_function iident , int explanatory_variable ,
                       int response_variable , const Vectors &vec)
:RegressionKernel(iident , (int)vec.min_value[explanatory_variable] , (int)vec.max_value[explanatory_variable])

{
  int i;


  vectors = vec.select_variable(explanatory_variable , response_variable);

  nb_vector = vectors->nb_vector;

  residual = new double[nb_vector];
  for (i = 0;i < nb_vector;i++) {
    residual[i] = 0.;
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Copy of a Regression object.
 *
 *  \param[in] regression reference on a Regression object.
 */
/*--------------------------------------------------------------*/

void Regression::copy(const Regression &regression)

{
  int i;


  if (regression.vectors) {
    vectors = new Vectors(*(regression.vectors));
  }
  else {
    vectors = NULL;
  }

  nb_vector = regression.nb_vector;

  residual = new double[nb_vector];
  for (i = 0;i < nb_vector;i++) {
    residual[i] = regression.residual[i];
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Constructor by copy of the Regression class.
 *
 *  \param[in] regression reference on a Regression object.
 */
/*--------------------------------------------------------------*/

Regression::Regression(const Regression &regression)

{
  RegressionKernel::copy(regression);
  copy(regression);
}


/*--------------------------------------------------------------*/
/**
 *  \brief Destruction of the data members of a Regression object.
 */
/*--------------------------------------------------------------*/

void Regression::remove()

{
  delete vectors;
  delete [] residual;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Destructor of the Regression class.
 */
/*--------------------------------------------------------------*/

Regression::~Regression()

{
  remove();
}


/*--------------------------------------------------------------*/
/**
 *  \brief Assignment operator of the Regression class.
 *
 *  \param[in] regression reference on a Regression object.
 *
 *  \return               Regression object.
 */
/*--------------------------------------------------------------*/

Regression& Regression::operator=(const Regression &regression)

{
  if (&regression != this) {
    remove();
    RegressionKernel::remove();

    RegressionKernel::copy(regression);
    copy(regression);
  }

  return *this;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing on a single line of a Regression object.
 *
 *  \param[in,out] os stream.
 */
/*--------------------------------------------------------------*/

ostream& Regression::line_write(ostream &os) const

{
  ascii_parameter_print(os);
  os << STAT_label[ident == LINEAR_FUNCTION ? STATL_R_SQUARED : STATL_DETERMINATION_COEFF] << ": "
     << 1. - residual_square_sum_computation() / (vectors->covariance[1][1] * (nb_vector - 1));

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a Regression object.
 *
 *  \param[in,out] os         stream,
 *  \param[in]     exhaustive flag detail level.
 */
/*--------------------------------------------------------------*/

ostream& Regression::ascii_write(ostream &os , bool exhaustive) const

{
  int i;
  int buff , width[5];
  double t_value , residual_mean , residual_standard_deviation , *estimated_response ,
         *standard_residual , square_sum[3] , df[3] , mean_square[3] , standard_deviation[2];
  Test *test;
  ios_base::fmtflags format_flags;


  format_flags = os.setf(ios::right , ios::adjustfield);

  // writing of the marginal distributions

  os << nb_vector << " " << STAT_label[nb_vector == 1 ? STATL_VECTOR : STATL_VECTORS] << endl;

  for (i = 0;i < 2;i++) {
    os << "\n" << STAT_label[i == 0 ? STATL_EXPLANATORY_VARIABLE : STATL_RESPONSE_VARIABLE]
       << "   (" << STAT_label[STATL_MIN_VALUE] << ": " << vectors->min_value[i] << ", "
       << STAT_label[STATL_MAX_VALUE] << ": " << vectors->max_value[i] << ")" << endl;

    os << "\n" << STAT_label[STATL_MARGINAL] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";

    if (vectors->marginal_distribution[i]) {
      vectors->marginal_distribution[i]->ascii_characteristic_print(os , exhaustive);

      if ((vectors->marginal_distribution[i]->nb_value <= ASCII_NB_VALUE) || (exhaustive)) {
        os << "\n   | " << STAT_label[STATL_MARGINAL] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << endl;
        vectors->marginal_distribution[i]->ascii_print(os);
      }
    }

    else {
      os << STAT_label[STATL_SAMPLE_SIZE] << ": " << vectors->nb_vector << endl;

      os << STAT_label[STATL_MEAN] << ": " << vectors->mean[i] << "   "
         << STAT_label[STATL_VARIANCE] << ": " << vectors->covariance[i][i] << "   "
         << STAT_label[STATL_STANDARD_DEVIATION] << ": " << sqrt(vectors->covariance[i][i]) << endl;

      if ((exhaustive) && (vectors->covariance[i][i] > 0.)) {
        os << STAT_label[STATL_SKEWNESS_COEFF] << ": " << vectors->skewness_computation(i) << "   "
           << STAT_label[STATL_KURTOSIS_COEFF] << ": " << vectors->kurtosis_computation(i) << endl;
      }
    }
  }

  if (ident == LINEAR_FUNCTION) {
    os << "\n" << STAT_label[STATL_CORRELATION_COEFF] << ": "
       << vectors->covariance[0][1] / sqrt(vectors->covariance[0][0] * vectors->covariance[1][1]) << endl;
  }

  square_sum[0] = regression_square_sum_computation();
  square_sum[1] = residual_square_sum_computation();
  square_sum[2] = vectors->covariance[1][1] * (nb_vector - 1);

  if (ident == LINEAR_FUNCTION) {

    // writing of the regression function parameters

    df[0] = regression_df;
    df[1] = residual_df;
    df[2] = regression_df + residual_df;

    for (i = 0;i < 3;i++) {
      mean_square[i] = square_sum[i] / df[i];
    }

    standard_deviation[0] = sqrt(mean_square[1] * (1. / (double)nb_vector +
                                 vectors->mean[0] * vectors->mean[0] / ((nb_vector - 1) * vectors->covariance[0][0])));
    standard_deviation[1] = sqrt(mean_square[1] / ((nb_vector - 1) * vectors->covariance[0][0]));

    test = new Test(STUDENT , false , (int)df[1] , I_DEFAULT , D_DEFAULT);

    test->critical_probability = ref_critical_probability[0];
    test->t_value_computation();

//    os << "\n" << STAT_label[STATL_REFERENCE] << " " << STAT_label[STATL_T_VALUE] << ": "
//       << test->value << "   " << STAT_label[STATL_REFERENCE] << " "
//       << STAT_label[STATL_CRITICAL_PROBABILITY] << ": " << test->critical_probability << endl;

    t_value = test->value;
    delete test;

    os << "\n" << STAT_label[STATL_INTERCEPT] << ": " << parameter[0] << " ("
       << STAT_label[STATL_STANDARD_DEVIATION] << ": " << standard_deviation[0] << ")";

    if (standard_deviation[0] > 0) {
      os << "   " << parameter[0] - t_value * standard_deviation[0]
         << " <= " << STAT_label[STATL_INTERCEPT] << " <= "
         << parameter[0] + t_value * standard_deviation[0] << endl;

      test = new Test(STUDENT , false , (int)df[1] , I_DEFAULT , parameter[0] / standard_deviation[0]);
      test->t_critical_probability_computation();

      test->ascii_print(os , false , true);

      delete test;
    }

    else {
      os << endl;
    }

    os << "\n" << STAT_label[STATL_SLOPE] << ": " << parameter[1] << " ("
       << STAT_label[STATL_STANDARD_DEVIATION] << ": " << standard_deviation[1] << ")";

    if (standard_deviation[1] > 0) {
      os << "   " << parameter[1] - t_value * standard_deviation[1]
         << " <= " << STAT_label[STATL_SLOPE] << " <= "
         << parameter[1] + t_value * standard_deviation[1] << endl;

      test = new Test(STUDENT , false , (int)df[1] , I_DEFAULT , parameter[1] / standard_deviation[1]);
      test->t_critical_probability_computation();

      test->ascii_print(os , false , true);

      delete test;
    }

    else {
      os << endl;
    }
  }

  if (exhaustive) {
    ascii_print(os);
  }

  os << "\n" << STAT_label[ident == LINEAR_FUNCTION ? STATL_R_SQUARED : STATL_DETERMINATION_COEFF]
     << ": " << 1. - square_sum[1] / square_sum[2] << endl;

  switch (ident) {

  case NONPARAMETRIC_FUNCTION : {
    os << regression_df << " " << STAT_label[STATL_REGRESSION] << " " << STAT_label[STATL_FREEDOM_DEGREES] << "   "
       << residual_df << " " << STAT_label[STATL_RESIDUAL] << " " << STAT_label[STATL_FREEDOM_DEGREES];

#   ifdef DEBUG
    os << "   (" << regression_df + residual_df << ")";

    os << "\n" << STAT_label[STATL_REGRESSION] << " " << STAT_label[STATL_SQUARE_SUM] << ": " << square_sum[0]
       << "   " << STAT_label[STATL_RESIDUAL] << " " << STAT_label[STATL_SQUARE_SUM] << ": " << square_sum[1]
       << "   " << STAT_label[STATL_TOTAL] << " " << STAT_label[STATL_SQUARE_SUM] << ": " << square_sum[2];
#   endif

    os << endl;
    break;
  }

  case LINEAR_FUNCTION : {

#   ifdef DEBUG
    os << "\ntest: " << vectors->covariance[0][1] / sqrt(vectors->covariance[0][0] * vectors->covariance[1][1])
       << " | " << sqrt(square_sum[0] / square_sum[2]) << endl;
    os << square_sum[0] + square_sum[1] << " | " << square_sum[2] << endl;
#   endif

    // writing of the table of variance analysis

    width[0] = column_width(nb_vector - 1) + ASCII_SPACE;
    width[1] = column_width(3 , square_sum) + ASCII_SPACE;
    width[2] = column_width(3 , mean_square) + ASCII_SPACE;

    os << "\n" << STAT_label[STATL_VARIANCE_ANALYSIS] << endl;
    os << STAT_label[STATL_VARIATION_SOURCE] << " | " << STAT_label[STATL_FREEDOM_DEGREES] << " | "
       << STAT_label[STATL_SQUARE_SUM] << " | " << STAT_label[STATL_MEAN_SQUARE] << endl;
    for (i = 0;i < 3;i++) {
      switch (i) {
      case 0 :
        os << STAT_label[STATL_REGRESSION];
        break;
      case 1 :
        os << STAT_label[STATL_RESIDUAL] << "  ";
        break;
      case 2 :
        os << STAT_label[STATL_TOTAL] << "     ";
        break;
      }

      os << setw(width[0]) << df[i];
      os << setw(width[1]) << square_sum[i];
      os << setw(width[2]) << mean_square[i] << endl;
    }
    os << endl;

    if (mean_square[1] > 0) {
      test = new Test(FISHER , true , (int)df[0] , (int)df[1] , mean_square[0] / mean_square[1]);
      test->F_critical_probability_computation();

      test->ascii_print(os , false , false);

      delete test;
    }
    break;
  }
  }

  // writing of the mean and standard deviation of residuals

  residual_mean = residual_mean_computation();
  residual_standard_deviation = sqrt(residual_variance_computation(residual_mean));

  os << "\n";
  if (ident != LINEAR_FUNCTION) {
    os << STAT_label[STATL_RESIDUAL] << " " << STAT_label[STATL_MEAN] << ": " << residual_mean
       << "   ";
  }
  os << STAT_label[STATL_RESIDUAL] << " " << STAT_label[STATL_STANDARD_DEVIATION] << ": "
     << residual_standard_deviation << endl;

  if (exhaustive) {
    estimated_response = new double[nb_vector];
    for (i = 0;i < nb_vector;i++) {
      estimated_response[i] = point[vectors->int_vector[i][0] - min_value];
    }

    // computation of standardized residuals

    standard_residual = new double[nb_vector];
    for (i = 0;i < nb_vector;i++) {
      standard_residual[i] = residual[i] / residual_standard_deviation;
    }

    // computation of the column width

    width[0] = column_width((int)vectors->min_value[0] , (int)vectors->max_value[0]);

    if (vectors->type[1] == INT_VALUE) {
      width[1] = column_width((int)vectors->min_value[1] , (int)vectors->max_value[1]);
    }

    else {
      width[1] = 0;
      for (i = 0;i < nb_vector;i++) {
        buff = column_width(1 , vectors->real_vector[i] + 1);
        if (buff > width[1]) {
          width[1] = buff;
        }
      }
    }
    width[1] += ASCII_SPACE;

    width[2] = column_width(nb_vector , estimated_response) + ASCII_SPACE;
    width[3] = column_width(nb_vector , residual) + ASCII_SPACE;
    width[4] = column_width(nb_vector , standard_residual) + ASCII_SPACE;

    // writing of observed and theoretical responses, of residuals and of standardized residuals

    os << "\n" << STAT_label[STATL_EXPLANATORY_VARIABLE] << " | " << STAT_label[STATL_RESPONSE_VARIABLE]
       << " | " << STAT_label[STATL_ESTIMATION] << " | " << STAT_label[STATL_RESIDUAL]
       << " | " << STAT_label[STATL_STANDARDIZED_RESIDUAL] << endl;

    for (i = 0;i < nb_vector;i++) {
      os << setw(width[0]) << vectors->int_vector[i][0];
      if (vectors->type[1] == INT_VALUE) {
        os << setw(width[1]) << vectors->int_vector[i][1];
      }
      else {
        os << setw(width[1]) << vectors->real_vector[i][1];
      }
      os << setw(width[2]) << estimated_response[i];
      os << setw(width[3]) << residual[i];
      os << setw(width[4]) << standard_residual[i];
      os << "   (" << vectors->identifier[i] << ")" << endl;
    }

    delete [] estimated_response;
    delete [] standard_residual;
  }

  os.setf(format_flags , ios::adjustfield);

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a Regression object in a file.
 *
 *  \param[in] error      reference on a StatError object,
 *  \param[in] path       file path,
 *  \param[in] exhaustive flag detail level.
 *
 *  \return               error status.
 */
/*--------------------------------------------------------------*/

bool Regression::ascii_write(StatError &error , const string path ,
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
    ascii_write(out_file , exhaustive);
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a Regression object in a file at the spreadsheet format.
 *
 *  \param[in] error reference on a StatError object,
 *  \param[in] path  file path.
 *
 *  \return          error status.
 */
/*--------------------------------------------------------------*/

bool Regression::spreadsheet_write(StatError &error , const string path) const

{
  bool status;
  int i;
  double t_value , residual_mean , residual_standard_deviation , square_sum[3] ,
         df[3] , mean_square[3] , standard_deviation[2];
  Test *test;
  ofstream out_file(path.c_str());


  error.init();

  if (!out_file) {
    status = false;
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {
    status = true;

    // writing of the marginal distributions

    out_file << nb_vector << "\t" << STAT_label[nb_vector == 1 ? STATL_VECTOR : STATL_VECTORS] << endl;

    for (i = 0;i < 2;i++) {
      out_file << "\n" << STAT_label[i == 0 ? STATL_EXPLANATORY_VARIABLE : STATL_RESPONSE_VARIABLE]
               << "\t\t" << STAT_label[STATL_MIN_VALUE] << "\t" << vectors->min_value[i]
               << "\t\t" << STAT_label[STATL_MAX_VALUE] << "\t" << vectors->max_value[i] << endl;

      out_file << "\n" << STAT_label[STATL_MARGINAL] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\t";

      if (vectors->marginal_distribution[i]) {
        vectors->marginal_distribution[i]->spreadsheet_characteristic_print(out_file);

        out_file << "\n\t" << STAT_label[STATL_MARGINAL] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << endl;
        vectors->marginal_distribution[i]->spreadsheet_print(out_file);
      }

      else {
        out_file << STAT_label[STATL_SAMPLE_SIZE] << "\t" << vectors->nb_vector << endl;

        out_file << STAT_label[STATL_MEAN] << "\t" << vectors->mean[i] << "\t\t"
                 << STAT_label[STATL_VARIANCE] << "\t" << vectors->covariance[i][i] << "\t\t"
                 << STAT_label[STATL_STANDARD_DEVIATION] << "\t" << sqrt(vectors->covariance[i][i]) << endl;

        if (vectors->covariance[i][i] > 0.) {
          out_file << STAT_label[STATL_SKEWNESS_COEFF] << "\t" << vectors->skewness_computation(i) << "\t\t"
                   << STAT_label[STATL_KURTOSIS_COEFF] << "\t" << vectors->kurtosis_computation(i) << endl;
        }
      }
    }

    if (ident == LINEAR_FUNCTION) {
      out_file << "\n" << STAT_label[STATL_CORRELATION_COEFF] << "\t"
               << vectors->covariance[0][1] / sqrt(vectors->covariance[0][0] * vectors->covariance[1][1]) << endl;
    }

    square_sum[0] = regression_square_sum_computation();
    square_sum[1] = residual_square_sum_computation();
    square_sum[2] = vectors->covariance[1][1] * (nb_vector - 1);

    if (ident == LINEAR_FUNCTION) {

      // writing of the regression function parameters

      df[0] = regression_df;
      df[1] = residual_df;
      df[2] = regression_df + residual_df;

      for (i = 0;i < 3;i++) {
        mean_square[i] = square_sum[i] / df[i];
      }

      standard_deviation[0] = sqrt(mean_square[1] * (1. / (double)nb_vector +
                                   vectors->mean[0] * vectors->mean[0] / ((nb_vector - 1) * vectors->covariance[0][0])));
      standard_deviation[1] = sqrt(mean_square[1] / ((nb_vector - 1) * vectors->covariance[0][0]));

      test = new Test(STUDENT , false , (int)df[1] , I_DEFAULT , D_DEFAULT);

      test->critical_probability = ref_critical_probability[0];
      test->t_value_computation();

//      out_file << "\n" << STAT_label[STATL_REFERENCE] << " " << STAT_label[STATL_T_VALUE] << "\t"
//               << test->value << "\t\t" << STAT_label[STATL_REFERENCE] << " "
//               << STAT_label[STATL_CRITICAL_PROBABILITY] << "\t" << test->critical_probability << endl;

      t_value = test->value;
      delete test;

      out_file << "\n" << STAT_label[STATL_INTERCEPT] << "\t" << parameter[0] << "\t\t"
               << STAT_label[STATL_STANDARD_DEVIATION] << "\t" << standard_deviation[0];

      if (standard_deviation[0] > 0) {
        out_file << "\t\t" << parameter[0] - t_value * standard_deviation[0]
                 << " <= " << STAT_label[STATL_INTERCEPT] << " <= "
                 << parameter[0] + t_value * standard_deviation[0] << endl;

        test = new Test(STUDENT , false , (int)df[1] , I_DEFAULT , parameter[0] / standard_deviation[0]);
        test->t_critical_probability_computation();

        test->spreadsheet_print(out_file , true);

        delete test;
      }

      else {
        out_file << endl;
      }

      out_file << "\n" << STAT_label[STATL_SLOPE] << "\t" << parameter[1] << "\t\t"
               << STAT_label[STATL_STANDARD_DEVIATION] << "\t" << standard_deviation[1];

      if (standard_deviation[1] > 0) {
        out_file << "\t\t" << parameter[1] - t_value * standard_deviation[1]
                 << " <= " << STAT_label[STATL_SLOPE] << " <= "
                 << parameter[1] + t_value * standard_deviation[1] << endl;

        test = new Test(STUDENT , false , (int)df[1] , I_DEFAULT , parameter[1] / standard_deviation[1]);
        test->t_critical_probability_computation();

        test->spreadsheet_print(out_file , true);

        delete test;
      }

      else {
        out_file << endl;
      }
    }

    spreadsheet_print(out_file);

    out_file << "\n" << STAT_label[ident == LINEAR_FUNCTION ? STATL_R_SQUARED : STATL_DETERMINATION_COEFF]
             << "\t" << 1. - square_sum[1] / square_sum[2] << endl;

    switch (ident) {

    case NONPARAMETRIC_FUNCTION : {
      out_file << regression_df << "\t" << STAT_label[STATL_REGRESSION] << " " << STAT_label[STATL_FREEDOM_DEGREES] << "\t\t"
               << residual_df << "\t" << STAT_label[STATL_RESIDUAL] << " " << STAT_label[STATL_FREEDOM_DEGREES] << endl;
      break;
    }

    case LINEAR_FUNCTION : {

      // writing of the table of variance analysis

      out_file << "\n" << STAT_label[STATL_VARIANCE_ANALYSIS] << endl;
      out_file << STAT_label[STATL_VARIATION_SOURCE] << "\t" << STAT_label[STATL_FREEDOM_DEGREES] << "\t"
               << STAT_label[STATL_SQUARE_SUM] << "\t" << STAT_label[STATL_MEAN_SQUARE] << endl;
      for (i = 0;i < 3;i++) {
        switch (i) {
        case 0 :
          out_file << STAT_label[STATL_REGRESSION];
          break;
        case 1 :
          out_file << STAT_label[STATL_RESIDUAL];
          break;
       case 2 :
          out_file << STAT_label[STATL_TOTAL];
          break;
        }

        out_file << "\t" << df[i] << "\t" << square_sum[i] << "\t" << mean_square[i] << endl;
      }
      out_file << endl;

      if (mean_square[1] > 0) {
        test = new Test(FISHER , true , (int)df[0] , (int)df[1] , mean_square[0] / mean_square[1]);
        test->F_critical_probability_computation();

        test->spreadsheet_print(out_file , false);

        delete test;
      }
      break;
    }
    }

    // writing of the mean and standard deviation of residuals

    residual_mean = residual_mean_computation();
    residual_standard_deviation = sqrt(residual_variance_computation(residual_mean));

    out_file << "\n";
    if (ident != LINEAR_FUNCTION) {
      out_file << STAT_label[STATL_RESIDUAL] << " " << STAT_label[STATL_MEAN] << "\t" << residual_mean
               << "\t\t";
    }
    out_file << STAT_label[STATL_RESIDUAL] << " " << STAT_label[STATL_STANDARD_DEVIATION] << "\t"
             << residual_standard_deviation << endl;

    // writing of observed and theoretical responses, of residuals and of standardized residuals

    out_file << "\n" << STAT_label[STATL_EXPLANATORY_VARIABLE] << "\t" << STAT_label[STATL_RESPONSE_VARIABLE]
             << "\t" << STAT_label[STATL_ESTIMATION] << "\t" << STAT_label[STATL_RESIDUAL]
             << "\t" << STAT_label[STATL_STANDARDIZED_RESIDUAL] << endl;

    for (i = 0;i < nb_vector;i++) {
      out_file << vectors->int_vector[i][0] << "\t";
      if (vectors->type[1] == INT_VALUE) {
        out_file << vectors->int_vector[i][1] << "\t";
      }
      else {
        out_file << vectors->real_vector[i][1] << "\t";
      }
      out_file << point[vectors->int_vector[i][0] - min_value] << "\t" << residual[i] << "\t"
               << residual[i] / residual_standard_deviation << "\t\t" << vectors->identifier[i] << endl;
    }
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Plot of a Regression object using Gnuplot.
 *
 *  \param[in] error  reference on a StatError object,
 *  \param[in] prefix file prefix,
 *  \param[in] title  figure title.
 *
 *  \return           error status.
 */
/*--------------------------------------------------------------*/

bool Regression::plot_write(StatError &error , const char *prefix ,
                            const char *title) const

{
  bool status;
  int i , j , k;
  int **frequency;
  double residual_mean , residual_standard_deviation , min_standard_residual ,
         max_standard_residual , min_response , max_response , threshold ,
         *standard_residual;
  ostringstream data_file_name[2];


  error.init();

  // writing of the data files

  data_file_name[0] << prefix << 0 << ".dat";
  status = plot_print((data_file_name[0].str()).c_str());

  if (!status) {
    error.update(STAT_error[STATR_FILE_PREFIX]);
  }

  else {

    // computation of the standardized residuals

    residual_mean = residual_mean_computation();
    residual_standard_deviation = sqrt(residual_variance_computation(residual_mean));

    standard_residual = new double[nb_vector];
    min_standard_residual = 0.;
    max_standard_residual = 0.;

    for (i = 0;i < nb_vector;i++) {
      standard_residual[i] = residual[i] / residual_standard_deviation;
      if (standard_residual[i] < min_standard_residual) {
        min_standard_residual = standard_residual[i];
      }
      if (standard_residual[i] > max_standard_residual) {
        max_standard_residual = standard_residual[i];
      }
    }

    data_file_name[1] << prefix << 1 << ".dat";
    vectors->plot_print((data_file_name[1].str()).c_str() , standard_residual);
    delete [] standard_residual;

    // writing of the script files

    min_response = MIN(min_computation() , vectors->min_value[1]);
    max_response = MAX(max_computation() , vectors->max_value[1]);

    normal dist;
    threshold = quantile(complement(dist , 0.025));

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

      out_file << "set xlabel \"" << STAT_label[STATL_EXPLANATORY_VARIABLE] << "\"" << endl;
      out_file << "set ylabel \"" << STAT_label[STATL_RESPONSE_VARIABLE] << "\"" << endl;

      if (max_value - min_value < TIC_THRESHOLD) {
        out_file << "set xtics " << MIN(min_value , 0) << ",1" << endl;
      }
      if (max_response - min_response < TIC_THRESHOLD) {
        out_file << "set ytics " << MIN(min_response , 0) << ",1" << endl;
      }

      if (((vectors->marginal_distribution[0]) &&
           (vectors->marginal_distribution[0]->nb_value <= PLOT_NB_VALUE)) &&
          ((vectors->marginal_distribution[1]) &&
           (vectors->marginal_distribution[1]->nb_value <= PLOT_NB_VALUE))) {
        frequency = vectors->joint_frequency_computation(0 , 1);

        for (j = vectors->marginal_distribution[0]->offset;j < vectors->marginal_distribution[0]->nb_value;j++) {
          for (k = vectors->marginal_distribution[1]->offset;k < vectors->marginal_distribution[1]->nb_value;k++) {
            if (frequency[j][k] > 0) {
              out_file << "set label \"" << frequency[j][k] << "\" at " << j << ","
                       << k << endl;
            }
          }
        }

        for (j = 0;j < vectors->marginal_distribution[0]->nb_value;j++) {
          delete [] frequency[j];
        }
        delete [] frequency;
      }

      if ((min_value >= 0) && (max_value - min_value > min_value * PLOT_RANGE_RATIO)) {
        out_file << "plot [" << 0;
      }
      else {
        out_file << "plot [" << min_value;
      }
      out_file << ":" << MAX(max_value , min_value + 1) << "] [";
      if ((min_response >= 0.) && (max_response - min_response > min_response * PLOT_RANGE_RATIO)) {
        out_file << 0;
      }
      else {
        out_file << min_response;
      }
      out_file << ":" << MAX(max_response , min_response + 1) << "] \""
               << label((data_file_name[1].str()).c_str()) << "\" using 1:2 notitle with points,\\" << endl;
      if (ident == NONPARAMETRIC_FUNCTION) {
        out_file << "\"" << label((data_file_name[0].str()).c_str()) << "\" using 1:2 notitle with lines" << endl;
      }
      else {
        out_file << "\"" << label((data_file_name[0].str()).c_str()) << "\" using 1:2 title \"";
        ascii_formal_print(out_file);
        out_file << "\" with lines" << endl;
      }

      out_file << "unset label" << endl;

      out_file << "set xlabel" << endl;
      out_file << "set ylabel" << endl;

      if (max_value - min_value < TIC_THRESHOLD) {
        out_file << "set xtics autofreq" << endl;
      }
      if (max_response - min_response < TIC_THRESHOLD) {
        out_file << "set ytics autofreq" << endl;
      }

      if (i == 0) {
        out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
      }
      out_file << endl;

      out_file << "set xlabel \"" << STAT_label[STATL_EXPLANATORY_VARIABLE] << "\"" << endl;
      out_file << "set ylabel \"" << STAT_label[STATL_STANDARDIZED_RESIDUAL] << "\"" << endl;

      if (max_value - min_value < TIC_THRESHOLD) {
        out_file << "set xtics 0,1" << endl;
      }

      out_file << "plot [" << vectors->min_value[0] << ":" << MAX(vectors->max_value[0] , vectors->min_value[0] + 1)
               << "] [" << MIN(min_standard_residual , -threshold) << ":" << MAX(max_standard_residual , threshold)
               << "] \"" << label((data_file_name[1].str()).c_str()) << "\" using 1:3 notitle with points";
      if (min_standard_residual < -threshold) {
        out_file << ",\\\n" << -threshold << " notitle with lines";
      }
      if (max_standard_residual > threshold) {
        out_file << ",\\\n" << threshold << " notitle with lines";
      }
      out_file << endl;

      out_file << "set xlabel" << endl;
      out_file << "set ylabel" << endl;

      if (max_value - min_value < TIC_THRESHOLD) {
        out_file << "set xtics autofreq" << endl;
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
 *  \brief Plot of a Regression object.
 *
 *  \return MultiPlotSet object.
 */
/*--------------------------------------------------------------*/

MultiPlotSet* Regression::get_plotable() const

{
  int i;
  int xmin , nb_plot;
  double ymin , min_response , max_response , residual_mean , residual_standard_deviation ,
         min_standard_residual , max_standard_residual , threshold , *standard_residual;
  ostringstream title , legend;
  MultiPlotSet *plot_set;


  plot_set = new MultiPlotSet(2);
  MultiPlotSet &plot = *plot_set;

  title.str("");
  title << STAT_label[ident == LINEAR_FUNCTION ? STATL_LINEAR : STATL_NONPARAMETRIC] << " "
        << STAT_label[STATL_REGRESSION];
  plot.title = title.str();

  plot.border = "15 lw 0";

  // regression function and data

  if ((min_value >= 0) && (max_value - min_value > min_value * PLOT_RANGE_RATIO)) {
    xmin = 0.;
  }
  else {
    xmin = min_value;
  }
  plot[0].xrange = Range(xmin , MAX(max_value , min_value + 1));

  min_response = MIN(min_computation() , vectors->min_value[1]);
  max_response = MAX(max_computation() , vectors->max_value[1]);

  if ((min_response >= 0.) && (max_response - min_response > min_response * PLOT_RANGE_RATIO)) {
    ymin = 0.;
  }
  else {
    ymin = min_response;
  }
  plot[0].yrange = Range(ymin , MAX(max_response , min_response + 1));

  if (max_value - min_value < TIC_THRESHOLD) {
    plot[0].xtics = 1;
  }
  if (max_response - min_response < TIC_THRESHOLD) {
    plot[0].ytics = 1;
  }

  plot[0].xlabel = STAT_label[STATL_EXPLANATORY_VARIABLE];
  plot[0].ylabel = STAT_label[STATL_RESPONSE_VARIABLE];

  if (((vectors->marginal_distribution[0]) &&
       (vectors->marginal_distribution[0]->nb_value <= PLOT_NB_VALUE)) &&
      ((vectors->marginal_distribution[1]) &&
       (vectors->marginal_distribution[1]->nb_value <= PLOT_NB_VALUE))) {
    plot[0].resize(3);
  }
  else {
    plot[0].resize(2);
  }

  plot[0][0].style = "points";

  vectors->plotable_write(plot[0][0] , 0 , 1);

  if (ident == LINEAR_FUNCTION) {
    legend.str("");
    ascii_formal_print(legend);
    plot[0][1].legend = legend.str();
  }

  plot[0][1].style = "lines";

  plotable_write(plot[0][1]);

  if (((vectors->marginal_distribution[0]) &&
       (vectors->marginal_distribution[0]->nb_value <= PLOT_NB_VALUE)) &&
      ((vectors->marginal_distribution[1]) &&
       (vectors->marginal_distribution[1]->nb_value <= PLOT_NB_VALUE))) {
    plot[0][2].label = "true";

    vectors->plotable_frequency_write(plot[0][2] , 0 , 1);
  }

  // standardized residuals 

  residual_mean = residual_mean_computation();
  residual_standard_deviation = sqrt(residual_variance_computation(residual_mean));

  standard_residual = new double[nb_vector];
  min_standard_residual = 0.;
  max_standard_residual = 0.;

  for (i = 0;i < nb_vector;i++) {
    standard_residual[i] = residual[i] / residual_standard_deviation;
    if (standard_residual[i] < min_standard_residual) {
      min_standard_residual = standard_residual[i];
    }
    if (standard_residual[i] > max_standard_residual) {
      max_standard_residual = standard_residual[i];
    }
  }

  normal dist;
  threshold = quantile(complement(dist , 0.025));

  plot[1].xrange = Range(vectors->min_value[0] ,
                         MAX(vectors->max_value[0] , vectors->min_value[0] + 1));
  plot[1].yrange = Range(MIN(min_standard_residual , -threshold) ,
                         MAX(max_standard_residual , threshold));

  if (max_value - min_value < TIC_THRESHOLD) {
    plot[1].xtics = 1;
  }

  plot[1].xlabel = STAT_label[STATL_EXPLANATORY_VARIABLE];
  plot[1].ylabel = STAT_label[STATL_STANDARDIZED_RESIDUAL];

  nb_plot = 1;
  if (min_standard_residual < -threshold) {
    nb_plot++;
  }
  if (max_standard_residual > threshold) {
    nb_plot++;
  }
  plot[1].resize(nb_plot);

  plot[1][0].style = "points";

  for (i = 0;i < vectors->nb_vector;i++) {
    plot[1][0].add_point(vectors->int_vector[i][0] , standard_residual[i]);
  }

  i = 1;
  if (min_standard_residual < -threshold) {
    plot[1][i].style = "lines";

    plot[1][i].add_point(min_value , -threshold);
    plot[1][i].add_point(max_value , -threshold);
    i++;
  }

  if (max_standard_residual > threshold) {
    plot[1][i].style = "lines";

    plot[1][i].add_point(min_value , threshold);
    plot[1][i].add_point(max_value , threshold);
  }

  delete [] standard_residual;

  return plot_set;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the variation explained by a regression function.
 *
 *  \return regression square sum.
 */
/*--------------------------------------------------------------*/

double Regression::regression_square_sum_computation() const

{
  int i;
  double regression_square_sum , diff;


  regression_square_sum = 0.;
  for (i = 0;i < nb_vector;i++) {
    diff = point[vectors->int_vector[i][0] - min_value] - vectors->mean[1];
    regression_square_sum += diff * diff;
  }

  return regression_square_sum;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the residuals.
 */
/*--------------------------------------------------------------*/

void Regression::residual_computation()

{
  int i;


  if (vectors->type[1] == INT_VALUE) {
    for (i = 0;i < nb_vector;i++) {
      residual[i] = vectors->int_vector[i][1] - point[vectors->int_vector[i][0] - min_value];
    }
  }
  else {
    for (i = 0;i < nb_vector;i++) {
      residual[i] = vectors->real_vector[i][1] - point[vectors->int_vector[i][0] - min_value];
    }
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the residual mean.
 *
 *  \return residual mean.
 */
/*--------------------------------------------------------------*/

double Regression::residual_mean_computation() const

{
  int i;
  double residual_mean;


  residual_mean = 0.;
  for (i = 0;i < nb_vector;i++) {
    residual_mean += residual[i];
  }
  residual_mean /= nb_vector;

  return residual_mean;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the residual variance.
 *
 *  \param[in] residual_mean residual mean.
 *
 *  \return                  residual variance.
 */
/*--------------------------------------------------------------*/

double Regression::residual_variance_computation(double residual_mean) const

{
  int i;
  double residual_variance = D_DEFAULT , diff;


  if (residual_df > 0.) {
    residual_variance = 0.;
    for (i = 0;i < nb_vector;i++) {
      diff = residual[i] - residual_mean;
      residual_variance += diff * diff;
    }
    residual_variance /= residual_df;
  }

  return residual_variance;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the residual square sum.
 *
 *  \return residual square sum.
 */
/*--------------------------------------------------------------*/

double Regression::residual_square_sum_computation() const

{
  int i;
  double residual_square_sum;


  residual_square_sum = 0.;
  for (i = 0;i < nb_vector;i++) {
    residual_square_sum += residual[i] * residual[i];
  }

  return residual_square_sum;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Linear regression.
 *
 *  \param[in] error                reference on a StatError object,
 *  \param[in] explanatory_variable explanatory variable index,
 *  \param[in] response_variable    response variable index.
 *
 *  \return                         Regression object.
 */
/*--------------------------------------------------------------*/

Regression* Vectors::linear_regression(StatError &error , int explanatory_variable ,
                                       int response_variable) const

{
  bool status = true;
  Regression *regression;


  regression = NULL;
  error.init();

  if (nb_vector < 3) {
    status = false;
    ostringstream correction_message;
    correction_message << STAT_error[STATR_GREATER_THAN] << " " << 2;
    error.correction_update(STAT_error[STATR_NB_VECTOR] , (correction_message.str()).c_str());
  }

  if (explanatory_variable == response_variable) {
    status = false;
    error.update(STAT_error[STATR_VARIABLE_INDICES]);
  }

  if ((explanatory_variable < 1) || (explanatory_variable > nb_variable)) {
    status = false;
    ostringstream error_message;
    error_message << STAT_label[STATL_VARIABLE] << " " << explanatory_variable << ": "
                  << STAT_error[STATR_VARIABLE_INDEX];
    error.update((error_message.str()).c_str());
  }
  else if (type[explanatory_variable - 1] != INT_VALUE) {
    status = false;
    ostringstream error_message;
    error_message << STAT_label[STATL_VARIABLE] << " " << explanatory_variable << ": "
                  << STAT_error[STATR_VARIABLE_TYPE];
    error.correction_update((error_message.str()).c_str() , STAT_variable_word[INT_VALUE]);
  }

  if ((response_variable < 1) || (response_variable > nb_variable)) {
    status = false;
    ostringstream error_message;
    error_message << STAT_label[STATL_VARIABLE] << " " << response_variable << ": "
                  << STAT_error[STATR_VARIABLE_INDEX];
    error.update((error_message.str()).c_str());
  }

  if (status) {
    explanatory_variable--;
    response_variable--;

    regression = new Regression(LINEAR_FUNCTION , explanatory_variable , response_variable , *this);

    regression->regression_df = 1.;
    regression->residual_df = nb_vector - 2;

    // estimation of linear regression parameters (intercept and slope)

    regression->parameter[1] = regression->vectors->covariance[0][1] / regression->vectors->covariance[0][0];
    regression->parameter[0] = regression->vectors->mean[1] - regression->vectors->mean[0] * regression->parameter[1];

    // computation of the regression function and residuals

    regression->computation();
    regression->residual_computation();
  }

  return regression;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Nonparametric moving-average regression.
 *
 *  \param[in] error                reference on a StatError object,
 *  \param[in] explanatory_variable explanatory variable index,
 *  \param[in] response_variable    response variable index,
 *  \param[in] nb_point             filter half width,
 *  \param[in] filter               filter,
 *  \param[in] algorithm            response computation (AVERAGING/LEAST_SQUARES).
 *
 *  \return                         Regression object.
 */
/*--------------------------------------------------------------*/

Regression* Vectors::moving_average(StatError &error , int explanatory_variable ,
                                    int response_variable , int nb_point , double *filter ,
                                    moving_average_method algorithm) const

{
  bool status = true;
  int i , j , k;
  int width , min , max , min_index , max_index , value_index , *index;
  double norm , diff , local_variance , local_covariance , local_mean[2] , *ppoint ,
         *weight , **smoother_matrix;
  Regression *regression;
  Vectors *vec;


  regression = NULL;
  error.init();

  if ((nb_vector < 3) || (nb_vector > REGRESSION_NB_VECTOR)) {
    status = false;
    error.update(STAT_error[STATR_NB_VECTOR]);
  }

  if (explanatory_variable == response_variable) {
    status = false;
    error.update(STAT_error[STATR_VARIABLE_INDICES]);
  }

  if ((explanatory_variable < 1) || (explanatory_variable > nb_variable)) {
    status = false;
    ostringstream error_message;
    error_message << STAT_label[STATL_VARIABLE] << " " << explanatory_variable << ": "
                  << STAT_error[STATR_VARIABLE_INDEX];
    error.update((error_message.str()).c_str());
  }

  else {
    explanatory_variable--;

    if (type[explanatory_variable] != INT_VALUE) {
      status = false;
      ostringstream error_message;
      error_message << STAT_label[STATL_VARIABLE] << " " << explanatory_variable + 1 << ": "
                    << STAT_error[STATR_VARIABLE_TYPE];
      error.correction_update((error_message.str()).c_str() , STAT_variable_word[INT_VALUE]);
    }

    if (!marginal_distribution[explanatory_variable]) {
      status = false;
      ostringstream error_message;
      error_message << STAT_label[STATL_VARIABLE] << " " << explanatory_variable + 1 << " "
                    << STAT_error[STATR_SHIFTED_SCALED];
      error.update((error_message.str()).c_str());
    }

    else {
      for (i = marginal_distribution[explanatory_variable]->offset;i < marginal_distribution[explanatory_variable]->nb_value;i++) {
        if (marginal_distribution[explanatory_variable]->frequency[i] == 0) {
          status = false;
          ostringstream error_message;
          error_message << STAT_label[STATL_VARIABLE] << " " << explanatory_variable + 1 << ": "
                        << STAT_error[STATR_MISSING_VALUE] << " " << i;
          error.update((error_message.str()).c_str());
        }
      }
    }
  }

  if ((response_variable < 1) || (response_variable > nb_variable)) {
    status = false;
    ostringstream error_message;
    error_message << STAT_label[STATL_VARIABLE] << " " << response_variable << ": "
                  << STAT_error[STATR_VARIABLE_INDEX];
    error.update((error_message.str()).c_str());
  }

  if ((nb_point == 0) && (algorithm == LEAST_SQUARES)) {
    status = false;
    error.update(STAT_error[STATR_LEAST_SQUARE_ALGORITHM]);
  }

  if (status) {
    response_variable--;

    regression = new Regression(NONPARAMETRIC_FUNCTION , explanatory_variable , response_variable , *this);
    vec = regression->vectors;

    weight = new double[nb_vector];

    smoother_matrix = new double*[nb_vector];
    for (i = 0;i < nb_vector;i++) {
      smoother_matrix[i] = new double[nb_vector];
      for (j = 0;j < nb_vector;j++) {
        smoother_matrix[i][j] = 0.;
      }
    }

    // Sort of the individuals by increasing explanatory values

    index = vec->order_computation(0);

#   ifdef DEBUG
    cout << "\ntest" << endl;
#   endif

    ppoint = regression->point;
    for (i = regression->min_value;i <= regression->max_value;i++) {

      // computation of bounds

      switch (algorithm) {

      case AVERAGING : {
        width = MIN(i - regression->min_value , regression->max_value - i);
        if (width > nb_point) {
          width = nb_point;
        }
        min = i - width;
        max = i + width;
        break;
      }

      case LEAST_SQUARES : {
        min = MAX(i - nb_point , regression->min_value);
        max = MIN(i + nb_point , regression->max_value);
        break;
      }
      }

      for (j = 0;j < nb_vector;j++) {
        if (vec->int_vector[index[j]][0] == min) {
          min_index = j;
          break;
        }
      }
      for (j = nb_vector - 1;j >= 0;j--) {
        if (vec->int_vector[index[j]][0] == max) {
          max_index = j;
          break;
        }
      }

      for (j = min_index;j < nb_vector;j++) {
        if (vec->int_vector[index[j]][0] == i) {
          value_index = j;
          break;
        }
      }

      switch (algorithm) {

      // estimation by averaging

      case AVERAGING : {
        local_mean[1] = 0.;
        norm = 0.;

        if (vec->type[1] == INT_VALUE) {
          for (j = min_index;j <= max_index;j++) {
            weight[j] = filter[vec->int_vector[index[j]][0] - (i - nb_point)];
            local_mean[1] += weight[j] * vec->int_vector[index[j]][1];
            norm += weight[j];
          }
        }

        else {
          for (j = min_index;j <= max_index;j++) {
            weight[j] = filter[vec->int_vector[index[j]][0] - (i - nb_point)];
            local_mean[1] += weight[j] * vec->real_vector[index[j]][1];
            norm += weight[j];
          }
        }

#       ifdef DEBUG
        cout << local_mean[1] / norm << " (" << min << ", " << max << " | "
             << min_index << ", " << max_index << ")" << endl;
#       endif

        *ppoint++ = local_mean[1] / norm;

        // computation of the smoother matrix

        for (j = 0;j < vec->marginal_distribution[0]->frequency[i];j++) {
          for (k = min_index;k <= max_index;k++) {
            smoother_matrix[value_index + j][k] = weight[k] / norm;
          }
        }
        break;
      }

      // least-square estimation

      case LEAST_SQUARES : {
        local_mean[0] = 0.;
        local_mean[1] = 0.;
        norm = 0.;

        if (vec->type[1] == INT_VALUE) {
          for (j = min_index;j <= max_index;j++) {
            weight[j] = filter[vec->int_vector[index[j]][0] - (i - nb_point)];
            local_mean[0] += weight[j] * vec->int_vector[index[j]][0];
            local_mean[1] += weight[j] * vec->int_vector[index[j]][1];
            norm += weight[j];
          }
        }

        else {
          for (j = min_index;j <= max_index;j++) {
            weight[j] = filter[vec->int_vector[index[j]][0] - (i - nb_point)];
            local_mean[0] += weight[j] * vec->int_vector[index[j]][0];
            local_mean[1] += weight[j] * vec->real_vector[index[j]][1];
            norm += weight[j];
          }
        }

        local_mean[0] /= norm;
        local_mean[1] /= norm;

        local_variance = 0.;
        local_covariance = 0.;

        if (vec->type[1] == INT_VALUE) {
          for (j = min_index;j <= max_index;j++) {
            diff = vec->int_vector[index[j]][0] - local_mean[0];
            local_variance += weight[j] * diff * diff;

            local_covariance += weight[j] * (vec->int_vector[index[j]][0] - local_mean[0]) *
                                (vec->int_vector[index[j]][1] - local_mean[1]);
          }
        }

        else {
          for (j = min_index;j <= max_index;j++) {
            diff = vec->int_vector[index[j]][0] - local_mean[0];
            local_variance += weight[j] * diff * diff;

            local_covariance += weight[j] * (vec->int_vector[index[j]][0] - local_mean[0]) *
                                (vec->real_vector[index[j]][1] - local_mean[1]);
          }
        }

        local_variance /= norm;
        local_covariance /= norm;

#       ifdef DEBUG
        cout << local_mean[1] + (i - local_mean[0]) * local_covariance / local_variance
             << " (" << min << ", " << max << " | " << min_index << ", " << max_index << ")" << endl;
#       endif

        *ppoint++ = local_mean[1] + (i - local_mean[0]) * local_covariance / local_variance;

        // computation of the smoother matrix

        for (j = 0;j < vec->marginal_distribution[0]->frequency[i];j++) {
          for (k = min_index;k <= max_index;k++) {
            smoother_matrix[value_index + j][k] = weight[k] / norm * (1. + (i - local_mean[0]) *
                                                   (vec->int_vector[index[k]][0] - local_mean[0]) / local_variance);
          }
        }
        break;
      }
      }
    }

#   ifdef DEBUG
    cout << endl;
#   endif

    // computation of the equivalent degrees of freedom

    regression->regression_df = 0.;
    regression->residual_df = nb_vector;

    for (i = 0;i < nb_vector;i++) {
      regression->regression_df += smoother_matrix[i][i];
      regression->residual_df -= 2. * smoother_matrix[i][i];
      for (j = 0;j < nb_vector;j++) {
        regression->residual_df += smoother_matrix[i][j] * smoother_matrix[i][j];
      }
    }

    // computation of residuals

    regression->residual_computation();

    delete [] weight;

    for (i = 0;i < nb_vector;i++) {
      delete [] smoother_matrix[i];
    }
    delete [] smoother_matrix;

    delete [] index;
  }

  return regression;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Nonparametric moving-average regression.
 *
 *  \param[in] error                reference on a StatError object,
 *  \param[in] explanatory_variable explanatory variable index,
 *  \param[in] response_variable    response variable index,
 *  \param[in] dist                 symmetric distribution,
 *  \param[in] algorithm            response computation (AVERAGING/LEAST_SQUARES).
 *
 *  \return                         Regression object.
 */
/*--------------------------------------------------------------*/

Regression* Vectors::moving_average(StatError &error , int explanatory_variable ,
                                    int response_variable , const Distribution &dist ,
                                    moving_average_method algorithm) const

{
  bool status = true;
  Regression *regression;


  regression = NULL;
  error.init();

  if ((dist.offset != 0) || ((dist.nb_value - dist.offset) % 2 == 0)) {
    status = false;
    error.correction_update(STAT_error[STATR_NB_VALUE] , STAT_error[STATR_ODD]);
  }
  if (fabs(dist.skewness_computation()) > SKEWNESS_ROUNDNESS) {
    status = false;
    error.update(STAT_error[STATR_NON_SYMMETRICAL_DISTRIBUTION]);
  }
  if (dist.complement > 0.) {
    status = false;
    error.update(STAT_error[STATR_UNPROPER_DISTRIBUTION]);
  }

  if (status) {
    regression = moving_average(error , explanatory_variable , response_variable ,
                                dist.nb_value / 2 , dist.mass , algorithm);
  }

  return regression;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Nonparametric k-nearest-neighbor regression.
 *
 *  \param[in] error                reference on a StatError object,
 *  \param[in] explanatory_variable explanatory variable index,
 *  \param[in] response_variable    response variable index,
 *  \param[in] span                 neighboring proportion with respect to the sample size,
 *  \param[in] weighting            flag neighbor weighting.
 *
 *  \return                         Regression object.
 */
/*--------------------------------------------------------------*/

Regression* Vectors::nearest_neighbor_smoother(StatError &error , int explanatory_variable ,
                                               int response_variable , double span ,
                                               bool weighting) const

{
  bool status = true , greater;
  int i , j , k;
  int nb_neighbor , min_index , max_index , value_index , value , frequency ,
      max_deviation , *index;
  double norm , var , diff , local_variance , local_covariance , local_mean[2] , *ppoint ,
         *weight , **smoother_matrix;
  Regression *regression;
  Vectors *vec;


  regression = NULL;
  error.init();

  if ((nb_vector < 3) || (nb_vector > REGRESSION_NB_VECTOR)) {
    status = false;
    error.update(STAT_error[STATR_NB_VECTOR]);
  }

  if (explanatory_variable == response_variable) {
    status = false;
    error.update(STAT_error[STATR_VARIABLE_INDICES]);
  }

  if ((explanatory_variable < 1) || (explanatory_variable > nb_variable)) {
    status = false;
    ostringstream error_message;
    error_message << STAT_label[STATL_VARIABLE] << " " << explanatory_variable << ": "
                  << STAT_error[STATR_VARIABLE_INDEX];
    error.update((error_message.str()).c_str());
  }
  else if (type[explanatory_variable - 1] != INT_VALUE) {
    status = false;
    ostringstream error_message;
    error_message << STAT_label[STATL_VARIABLE] << " " << explanatory_variable << ": "
                  << STAT_error[STATR_VARIABLE_TYPE];
    error.correction_update((error_message.str()).c_str() , STAT_variable_word[INT_VALUE]);
  }

  if ((response_variable < 1) || (response_variable > nb_variable)) {
    status = false;
    ostringstream error_message;
    error_message << STAT_label[STATL_VARIABLE] << " " << response_variable << ": "
                  << STAT_error[STATR_VARIABLE_INDEX];
    error.update((error_message.str()).c_str());
  }

  if (status) {
    explanatory_variable--;
    response_variable--;

    regression = new Regression(NONPARAMETRIC_FUNCTION , explanatory_variable , response_variable , *this);
    vec = regression->vectors;

    nb_neighbor = MAX((int)::round(nb_vector * span) , 2);
    if (weighting) {
      nb_neighbor++;
    }

    weight = new double[nb_vector];

    smoother_matrix = new double*[nb_vector];
    for (i = 0;i < nb_vector;i++) {
      smoother_matrix[i] = new double[nb_vector];
      for (j = 0;j < nb_vector;j++) {
        smoother_matrix[i][j] = 0.;
      }
    }

    // Sort of the individuals by increasing explanatory values

    index = vec->order_computation(0);

#   ifdef DEBUG
    cout << "\n";
    for (i = 0;i < nb_vector;i++) {
      cout << index[i] << " ";
    }
    cout << endl;
#   endif

#   ifdef DEBUG
    cout << "\ntest - voisinage : " << nb_neighbor << " " << regression->min_value
         << " " << regression->max_value << endl;
#   endif

    ppoint = regression->point;
    for (i = regression->min_value;i <= regression->max_value;i++) {

      // determination of the nearest neighbor

      min_index = I_DEFAULT;
      do {

        // determination of the current nearest value

        if (min_index == I_DEFAULT) {
          min_index = 0;
          while (vec->int_vector[index[min_index]][0] < i) {
            min_index++;
          }

          if ((vec->int_vector[index[min_index]][0] > i) &&
              (i - vec->int_vector[index[min_index - 1]][0] < vec->int_vector[index[min_index]][0] - i)) {
            min_index--;
            greater = false;
          }
          else {
            greater = true;
          }
          max_index = min_index;

          frequency = 0;
          if (vec->int_vector[index[min_index]][0] == i) {
            value_index = min_index;
            j = min_index;
            do {
              j++;
              frequency++;
            }
            while ((j < nb_vector) && (vec->int_vector[index[j]][0] == i));
          }
        }

        else {
          if (min_index == 0) {
            greater = true;
            max_index++;
          }
          else if (max_index == nb_vector - 1) {
            greater = false;
            min_index--;
          }
          else if (i - vec->int_vector[index[min_index - 1]][0] < vec->int_vector[index[max_index + 1]][0] - i) {
            greater = false;
            min_index--;
          }
          else {
            greater = true;
            max_index++;
          }
        }

        // selection of the individuals taking the current nearest value

        if (!greater) {
          value = vec->int_vector[index[min_index]][0];
          do {
            min_index--;
          }
          while ((min_index >= 0) && (vec->int_vector[index[min_index]][0] == value));
          min_index++;
        }

        else {
          value = vec->int_vector[index[max_index]][0];
          do {
            max_index++;
          }
          while ((max_index < nb_vector) && (vec->int_vector[index[max_index]][0] == value));
          max_index--;
        }
      }
      while ((max_index - min_index + 1 < nb_neighbor) ||
             (vec->int_vector[index[max_index]][0] - vec->int_vector[index[min_index]][0] < NEIGHBORHOOD) ||
             (i - vec->int_vector[index[min_index]][0] == vec->int_vector[index[max_index]][0] - i));

      if (weighting) {
        if (!greater) {
          max_deviation = i - vec->int_vector[index[min_index]][0];
        }
        else {
          max_deviation = vec->int_vector[index[max_index]][0] - i;
        }
      }

      // least-square estimation

      local_mean[0] = 0.;
      local_mean[1] = 0.;
      norm = 0.;

      if (vec->type[1] == INT_VALUE) {
        for (j = min_index;j <= max_index;j++) {
          if (weighting) {
            var = (double)(abs(vec->int_vector[index[j]][0] - i)) / (double)max_deviation;
            var = 1. - var * var * var;
            weight[j] = var * var * var;
          }
          else {
            weight[j] = 1.;
          }
          norm += weight[j];

          local_mean[0] += weight[j] * vec->int_vector[index[j]][0];
          local_mean[1] += weight[j] * vec->int_vector[index[j]][1];
        }
      }

      else {
        for (j = min_index;j <= max_index;j++) {
          if (weighting) {
            var = (double)(abs(vec->int_vector[index[j]][0] - i)) / (double)max_deviation;
            var = 1. - var * var * var;
            weight[j] = var * var * var;
          }
          else {
            weight[j] = 1.;
          }
          norm += weight[j];

          local_mean[0] += weight[j] * vec->int_vector[index[j]][0];
          local_mean[1] += weight[j] * vec->real_vector[index[j]][1];
        }
      }

      local_mean[0] /= norm;
      local_mean[1] /= norm;

      local_variance = 0.;
      local_covariance = 0.;

      if (vec->type[1] == INT_VALUE) {
        for (j = min_index;j <= max_index;j++) {
          diff = vec->int_vector[index[j]][0] - local_mean[0];
          local_variance += weight[j] * diff * diff;

          local_covariance += weight[j] * (vec->int_vector[index[j]][0] - local_mean[0]) *
                              (vec->int_vector[index[j]][1] - local_mean[1]);
        }
      }

      else {
        for (j = min_index;j <= max_index;j++) {
          diff = vec->int_vector[index[j]][0] - local_mean[0];
          local_variance += weight[j] * diff * diff;

          local_covariance += weight[j] * (vec->int_vector[index[j]][0] - local_mean[0]) *
                              (vec->real_vector[index[j]][1] - local_mean[1]);
        }
      }

      local_variance /= norm;
      local_covariance /= norm;

#     ifdef DEBUG
      cout << i << " " << frequency << " | ";
      cout << local_mean[0] << " " << local_mean[1] << " " << local_covariance << " " << local_variance << " | ";
      cout << local_mean[1] + (i - local_mean[0]) * local_covariance / local_variance
           << " (" << min_index << ", " << max_index << ")" << endl;
#     endif

      if (local_variance == 0.) {
        status = false;
        break;
      }

      *ppoint++ = local_mean[1] + (i - local_mean[0]) * local_covariance / local_variance;

      // computation of the smoother matrix

      for (j = 0;j < frequency;j++) {
        for (k = min_index;k <= max_index;k++) {
          smoother_matrix[value_index + j][k] = weight[k] / norm * (1. + (i - local_mean[0]) *
                                                 (vec->int_vector[index[k]][0] - local_mean[0]) / local_variance);
        }
      }
    }

    if (status) {

      // computation of the equivalent degrees of freedom

      regression->regression_df = 0.;
      regression->residual_df = nb_vector;

      for (i = 0;i < nb_vector;i++) {
        regression->regression_df += smoother_matrix[i][i];
        regression->residual_df -= 2. * smoother_matrix[i][i];
        for (j = 0;j < nb_vector;j++) {
          regression->residual_df += smoother_matrix[i][j] * smoother_matrix[i][j];
        }
      }

#     ifdef DEBUG
/*      cout << "\nsmoother matrix" << endl;
      for (i = 0;i < nb_vector;i++) {
        cout << vec->int_vector[index[i]][0] << " ";
      }
      cout << "\n\n";
      for (i = 0;i < nb_vector;i++) {
        for (j = 0;j < nb_vector;j++) {
          cout << smoother_matrix[i][j] << " ";
        }
        cout << endl;
      } */
#     endif

      // computation of residuals

      regression->residual_computation();
    }

    else {
      delete regression;
      regression = NULL;
      error.update(STAT_error[STATR_REGRESSION_FAILURE]);
    }

    delete [] weight;

    for (i = 0;i < nb_vector;i++) {
      delete [] smoother_matrix[i];
    }
    delete [] smoother_matrix;

    delete [] index;
  }

  return regression;
}


};  // namespace stat_tool
