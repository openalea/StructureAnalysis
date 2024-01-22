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
 *       $Id: continuous_parametric.cpp 8175 2010-02-18 10:30:35Z guedon $
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

#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>

#include <boost/tokenizer.hpp>

#include <boost/math/distributions/gamma.hpp>
#include <boost/math/distributions/inverse_gaussian.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/special_functions/bessel.hpp>

#include "vectors.h"
#include "stat_label.h"

using namespace std;
using namespace boost;
using namespace boost::math;


namespace stat_tool {



/*--------------------------------------------------------------*/
/**
 *  \brief Constructor of the ContinuousParametric class.
 *
 *  \param[in] iident            identifier,
 *  \param[in] ilocation         location parameter (Gaussian, inverse Gaussian, autoregressive model, von Mises),
 *                               shape parameter (gamma, zero-inflated gamma) or intercept (linear model),
 *  \param[in] idispersion       dispersion parameter (Gaussian, von Mises, linear model, autoregressive model) or
 *                               scale parameter (gamma, zero-inflated gamma, inverse Gaussian),
 *  \param[in] izero_probability zero probability (zero-inflated gamma), slope (linear model) or
 *                               autoregressive coefficient (autoregressive model),
 *  \param[in] iunit             angle unit (von Mises).
 */
/*--------------------------------------------------------------*/

ContinuousParametric::ContinuousParametric(continuous_parametric iident , double ilocation ,
                                           double idispersion , double izero_probability ,
                                           angle_unit iunit)

{
  ident = iident;

  location = ilocation;
  dispersion = idispersion;
  zero_probability = izero_probability;

  min_value = D_DEFAULT;
  max_value = D_DEFAULT;

  unit = iunit;
  slope_standard_deviation = D_DEFAULT;
  sample_size = 0;
  correlation = D_DEFAULT;
  cumul = NULL;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Copy of a ContinuousParametric object.
 *
 *  \param[in] dist reference on a ContinuousParametric object.
 */
/*--------------------------------------------------------------*/

void ContinuousParametric::copy(const ContinuousParametric &dist)

{
  ident = dist.ident;

  location = dist.location;
  dispersion = dist.dispersion;
  zero_probability = dist.zero_probability;

  min_value = dist.min_value;
  max_value = dist.max_value;

  unit = dist.unit;
  slope_standard_deviation = dist.slope_standard_deviation;
  sample_size = dist.sample_size;
  correlation = dist.correlation;

  if (dist.cumul) {
    int i;

    cumul = new double[VON_MISES_NB_STEP];

    for (i = 0;i < VON_MISES_NB_STEP;i++) {
      cumul[i] = dist.cumul[i];
    }
  }

  else {
    cumul = NULL;
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Destructor of the ContinuousParametric class.
 */
/*--------------------------------------------------------------*/

ContinuousParametric::~ContinuousParametric()

{
  delete [] cumul;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Assignment operator of the ContinuousParametric class.
 *
 *  \param[in] dist reference on a ContinuousParametric object.
 *
 *  \return         ContinuousParametric object.
 */
/*--------------------------------------------------------------*/

ContinuousParametric& ContinuousParametric::operator=(const ContinuousParametric &dist)

{
  if (&dist != this) {
    delete [] cumul;
    copy(dist);
  }

  return *this;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Analysis of the format of a ContinuousParametric object.
 *
 *  \param[in] error      reference on a StatError object,
 *  \param[in] in_file    stream,
 *  \param[in] line       reference on the file line index,
 *  \param[in] last_ident identifier of the last distribution in the list.
 *
 *  \return               ContinuousParametric object.
 */
/*--------------------------------------------------------------*/

ContinuousParametric* ContinuousParametric::parsing(StatError &error , ifstream &in_file ,
                                                    int &line , continuous_parametric last_ident)

{
  string buffer;
  size_t position;
  typedef tokenizer<char_separator<char>> tokenizer;
  char_separator<char> separator(" \t");
  bool status = true , lstatus;
  int i , j;
  continuous_parametric ident = GAUSSIAN;
  union {
    double shape;
    double location;
    double intercept;
  };
  union {
    double scale;
    double dispersion;
  };
  union {
    double zero_probability;
    double slope;
    double autoregressive_coeff;
  };
  ContinuousParametric *dist;


  dist = NULL;
  location = D_INF;
  dispersion = D_DEFAULT;
  slope = D_INF;

  while (getline(in_file , buffer)) {
    line++;

#   ifdef DEBUG
    cout << line << " " << buffer << endl;
#   endif

    position = buffer.find('#');
    if (position != string::npos) {
      buffer.erase(position);
    }
    i = 0;

    tokenizer tok_buffer(buffer , separator);

    for (tokenizer::iterator token = tok_buffer.begin();token != tok_buffer.end();token++) {

      // test distribution name

      if (i == 0) {
        for (j = GAMMA;j <= last_ident;j++) {
          if (*token == STAT_continuous_distribution_word[j]) {
            ident = (continuous_parametric)j;
            break;
          }
        }

        if (j == last_ident + 1) {
          status = false;
          error.update(STAT_parsing[STATP_DISTRIBUTION_NAME] , line , i + 1);
        }
      }

      // test parameter name

      else {
        switch ((i - 1) % 3) {

        case 0 : {
          switch ((i - 1) / 3) {

          // 1st parameter: shape parameter (gamma), mean (Gaussian, inverse Gaussian, autoregressive model),
          // mean direction (von Mises), zero probability (zero-inflated gamma), intercept (linear model)

          case 0 : {
            if ((ident == GAMMA) && (*token != STAT_word[STATW_SHAPE])) {
              status = false;
              error.correction_update(STAT_parsing[STATP_PARAMETER_NAME] , STAT_word[STATW_SHAPE] , line , i + 1);
            }

            if (((ident == INVERSE_GAUSSIAN) || (ident == GAUSSIAN) || (ident == AUTOREGRESSIVE_MODEL)) && (*token != STAT_word[STATW_MEAN])) {
              status = false;
              error.correction_update(STAT_parsing[STATP_PARAMETER_NAME] , STAT_word[STATW_MEAN] , line , i + 1);
            }

            if ((ident == VON_MISES) && (*token != STAT_word[STATW_MEAN_DIRECTION])) {
              status = false;
              error.correction_update(STAT_parsing[STATP_PARAMETER_NAME] , STAT_word[STATW_MEAN_DIRECTION] , line , i + 1);
            }

            if ((ident == ZERO_INFLATED_GAMMA) && (*token != STAT_word[STATW_ZERO_PROBABILITY])) {
              status = false;
              error.correction_update(STAT_parsing[STATP_PARAMETER_NAME] , STAT_word[STATW_ZERO_PROBABILITY] , line , i + 1);
            }

            if ((ident == LINEAR_MODEL) && (*token != STAT_word[STATW_INTERCEPT])) {
              status = false;
              error.correction_update(STAT_parsing[STATP_PARAMETER_NAME] , STAT_word[STATW_INTERCEPT] , line , i + 1);
            }
            break;
          }

          // 2nd parameter: scale parameter (gamma, inverse Gaussian), standard deviation (Gaussian),
          // concentration (von Mises), shape parameter (zero-inflated gamma), slope (linear model),
          // autoregressive coefficient (autoregressive model)

          case 1 : {
            if (((ident == GAMMA) || (ident == INVERSE_GAUSSIAN)) && (*token != STAT_word[STATW_SCALE])) {
              status = false;
              error.correction_update(STAT_parsing[STATP_PARAMETER_NAME] , STAT_word[STATW_SCALE] , line , i + 1);
            }

            if ((ident == GAUSSIAN) && (*token != STAT_word[STATW_STANDARD_DEVIATION])) {
              status = false;
              error.correction_update(STAT_parsing[STATP_PARAMETER_NAME] , STAT_word[STATW_STANDARD_DEVIATION] , line , i + 1);
            }

            if ((ident == VON_MISES) && (*token != STAT_word[STATW_CONCENTRATION])) {
              status = false;
              error.correction_update(STAT_parsing[STATP_PARAMETER_NAME] , STAT_word[STATW_CONCENTRATION] , line , i + 1);
            }

            if ((ident == ZERO_INFLATED_GAMMA) && (*token != STAT_word[STATW_SHAPE])) {
              status = false;
              error.correction_update(STAT_parsing[STATP_PARAMETER_NAME] , STAT_word[STATW_SHAPE] , line , i + 1);
            }

            if ((ident == LINEAR_MODEL) && (*token != STAT_word[STATW_SLOPE])) {
              status = false;
              error.correction_update(STAT_parsing[STATP_PARAMETER_NAME] , STAT_word[STATW_SLOPE] , line , i + 1);
            }

            if ((ident == AUTOREGRESSIVE_MODEL) && (*token != STAT_word[STATW_AUTOREGRESSIVE_COEFF])) {
              status = false;
              error.correction_update(STAT_parsing[STATP_PARAMETER_NAME] , STAT_word[STATW_AUTOREGRESSIVE_COEFF] , line , i + 1);
            }
            break;
          }

          // 3rd parameter: scale parameter (zero-inflated gamma), residual standard deviation (linear model, autoregressive model),

          case 2 : {
            if ((ident == ZERO_INFLATED_GAMMA) && (*token != STAT_word[STATW_SCALE])) {
              status = false;
              error.correction_update(STAT_parsing[STATP_PARAMETER_NAME] , STAT_word[STATW_SCALE] , line , i + 1);
            }

            if (((ident == LINEAR_MODEL) || (ident == AUTOREGRESSIVE_MODEL)) && (*token != STAT_word[STATW_STANDARD_DEVIATION])) {
              status = false;
              error.correction_update(STAT_parsing[STATP_PARAMETER_NAME] , STAT_word[STATW_STANDARD_DEVIATION] , line , i + 1);
            }
            break;
          }
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

        // test parameter value

        case 2 : {
          lstatus = true;

          switch ((i - 1) / 3) {

          // 1st parameter: shape parameter (gamma), mean (inverse Gaussian, Gaussian, autoregressive model),
          // mean direction (von Mises), zero probability (zero-inflated gamma), intercept (linear model)

          case 0 : {
            if (ident == GAMMA) {
/*              try {
                shape = stod(*token);   in C++ 11
              }
              catch(invalid_argument &arg) {
                lstatus = false;
              } */
              shape = atof(token->c_str());

              if ((lstatus) && (shape < 0.)) {
                lstatus = false;
              }
            }

            else if (ident == ZERO_INFLATED_GAMMA) {
/*              try {
                zero_probability = stod(*token);   in C++ 11
              }
              catch(invalid_argument &arg) {
                lstatus = false;
              } */
              zero_probability = atof(token->c_str());

              if ((lstatus) && ((zero_probability < 0.) || (zero_probability > 1.))) {
                lstatus = false;
              }
            }

            else if (ident == LINEAR_MODEL) {
/*              try {
                intercept = stod(*token);   in C++ 11
              }
              catch(invalid_argument &arg) {
                lstatus = false;
              } */
              intercept = atof(token->c_str());
            }

            else {
/*              try {
                location = stod(*token);   in C++ 11
              }
              catch(invalid_argument &arg) {
                lstatus = false;
              } */
              location = atof(token->c_str());

              if ((lstatus) && (ident == INVERSE_GAUSSIAN) && (location <= 0.)) {
                lstatus = false;
              }
            }
            break;
          }

          // 2nd parameter: scale parameter (gamma, inverse Gaussian), standard deviation (Gaussian),
          // concentration (von Mises), shape parameter (zero-inflated gamma), slope (linear model),
          // autoregressive coefficient (autoregressive model)

          case 1 : {
            if ((ident == GAMMA) || (ident == INVERSE_GAUSSIAN)) {
/*              try {
                scale = stod(*token);   in C++ 11
              }
              catch(invalid_argument &arg) {
                lstatus = false;
              } */
              scale = atof(token->c_str());

              if ((lstatus) && (scale <= 0.)) {
                lstatus = false;
              }
            }

            else if (ident == ZERO_INFLATED_GAMMA) {
/*              try {
                shape = stod(*token);   in C++ 11
              }
              catch(invalid_argument &arg) {
                lstatus = false;
              } */
              shape = atof(token->c_str());

              if ((lstatus) && (shape <= 0.)) {
                lstatus = false;
              }
            }

            else if (ident == LINEAR_MODEL) {
/*              try {
                slope = stod(*token);   in C++ 11
              }
              catch(invalid_argument &arg) {
                lstatus = false;
              } */
              slope = atof(token->c_str());
            }

            else if (ident == AUTOREGRESSIVE_MODEL) {
/*              try {
                autoregressive_coeff = stod(*token);   in C++ 11
              }
              catch(invalid_argument &arg) {
                lstatus = false;
              } */
              autoregressive_coeff = atof(token->c_str());

              if ((lstatus) && ((autoregressive_coeff < -1.) || (autoregressive_coeff > 1))) {
                lstatus = false;
              }
            }

            else {
/*              try {
                dispersion = stod(*token);   in C++ 11
              }
              catch(invalid_argument &arg) {
                lstatus = false;
              } */
              dispersion = atof(token->c_str());

              if ((lstatus) && (dispersion <= 0.)) {
                lstatus = false;
              }
            }
            break;
          }

          // 3rd parameter: scale parameter (zero-inflated gamma), residual standard deviation (linear model, autoregressive model)

          case 2 : {
            if (ident == ZERO_INFLATED_GAMMA) {
/*              try {
                scale = stod(*token);   in C++ 11
              }
              catch(invalid_argument &arg) {
                lstatus = false;
              } */
              scale = atof(token->c_str());

              if ((lstatus) && (scale <= 0.)) {
                lstatus = false;
              }
            }

            else if ((ident == LINEAR_MODEL) || (ident == AUTOREGRESSIVE_MODEL)) {
/*              try {
                dispersion = stod(*token);   in C++ 11
              }
              catch(invalid_argument &arg) {
                lstatus = false;
              } */
              dispersion = atof(token->c_str());

              if ((lstatus) && (dispersion <= 0.)) {
                lstatus = false;
              }
            }
            break;
          }
          }

          if (!lstatus) {
            status = false;
            error.update(STAT_parsing[STATP_PARAMETER_VALUE] , line , i + 1);
          }
          break;
        }
        }
      }

      i++;
    }

    if (i > 0) {
      if ((((ident == GAMMA) || (ident == INVERSE_GAUSSIAN) || (ident == GAUSSIAN) || (ident == VON_MISES)) && (i != 7)) ||
          (((ident == ZERO_INFLATED_GAMMA) || (ident == LINEAR_MODEL) || (ident == AUTOREGRESSIVE_MODEL)) && (i != 10))) {
        status = false;
        error.update(STAT_parsing[STATP_FORMAT] , line);
      }

      break;
    }
  }

  if (ident == I_DEFAULT) {
    status = false;
    error.update(STAT_parsing[STATP_FORMAT] , line);
  }

  if (status) {
    dist = new ContinuousParametric(ident , location , dispersion , zero_probability ,
                                    (location < 2 * M_PI ? RADIAN : DEGREE));
  }

  return dist;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of the parameters of a continuous distribution.
 *
 *  \param[in,out] os        stream,
 *  \param[in]     file_flag flag comment.
 */
/*--------------------------------------------------------------*/

ostream& ContinuousParametric::ascii_parameter_print(ostream &os , bool file_flag) const

{
  os << STAT_continuous_distribution_word[ident] << "   ";

  switch (ident) {

  case GAMMA : {
    os << STAT_word[STATW_SHAPE] << " : " << shape;
    if (shape > 0.) {
      os << "   " << STAT_word[STATW_SCALE] << " : " << scale;
    }
    break;
  }

  case INVERSE_GAUSSIAN : {
    os << STAT_word[STATW_MEAN] << " : " << location << "   "
       << STAT_word[STATW_SCALE] << " : " << scale;
    break;
  }

  case GAUSSIAN : {
    os << STAT_word[STATW_MEAN] << " : " << location << "   "
       << STAT_word[STATW_STANDARD_DEVIATION] << " : " << dispersion;
    break;
  }

  case VON_MISES : {
    os << STAT_word[STATW_MEAN_DIRECTION] << " : " << location << "   "
       << STAT_word[STATW_CONCENTRATION] << " : " << dispersion;
    break;
  }

  case ZERO_INFLATED_GAMMA : {
    os << STAT_word[STATW_ZERO_PROBABILITY] << " : " << zero_probability;
    if (zero_probability < 1.) {
      os << "   " << STAT_word[STATW_SHAPE] << " : " << shape
         << "   " << STAT_word[STATW_SCALE] << " : " << scale;
    }
    break;
  }

  case LINEAR_MODEL : {
    os << STAT_word[STATW_INTERCEPT] << " : " << intercept << "   "
       << STAT_word[STATW_SLOPE] << " : " << slope;

    if (sample_size > 0.) {
      Test test(STUDENT , false , sample_size , I_DEFAULT , D_DEFAULT);

      test.critical_probability = ref_critical_probability[0];
      test.t_value_computation();

      if ((!file_flag) && (slope_standard_deviation > 0.)) {
        os << "   (" << slope - test.value * slope_standard_deviation << ", "
           << slope + test.value * slope_standard_deviation << ")";
      }

      os << "   " << STAT_word[STATW_STANDARD_DEVIATION] << " : " << dispersion;

      if ((file_flag) && (slope_standard_deviation > 0.)) {
        os << "   # " << slope - test.value * slope_standard_deviation << ", "
           << slope + test.value * slope_standard_deviation;
      }

      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << STAT_label[STATL_CORRELATION_COEFF] << ": " << correlation << "   "
         << STAT_label[STATL_LIMIT_CORRELATION_COEFF] << ": -/+"
         << test.value / sqrt(test.value * test.value + sample_size) << "   "
         << STAT_label[STATL_CRITICAL_PROBABILITY] << ": " << test.critical_probability;
    }

    else {
      os << "   " << STAT_word[STATW_STANDARD_DEVIATION] << " : " << dispersion;
    }
    break;
  }

  case AUTOREGRESSIVE_MODEL : {
    os << STAT_word[STATW_MEAN] << " : " << location << "   "
       << STAT_word[STATW_AUTOREGRESSIVE_COEFF] << " : " << autoregressive_coeff;

    if (sample_size > 0.) {
      normal dist;
      double standard_normal_value = quantile(complement(dist , 0.025)) , standard_error;

      standard_error = standard_normal_value * sqrt((1. - autoregressive_coeff * autoregressive_coeff) / sample_size);

      if (!file_flag) {
        os << "   (" << MAX(autoregressive_coeff - standard_error , -1.) << ", " << MIN(autoregressive_coeff + standard_error , 1.)
           << "),   " << STAT_label[STATL_NULL_AUTOREGRESSIVE_COEFF_95_CONFIDENCE_LIMIT]
           << ": -/+" << standard_normal_value / sqrt(sample_size);
      }

      os << "   " << STAT_word[STATW_STANDARD_DEVIATION] << " : " << dispersion;

      if (file_flag) {
        os << "   # " << MAX(autoregressive_coeff - standard_error , -1.) << ", " << MIN(autoregressive_coeff + standard_error , 1.)
           << ",   " << STAT_label[STATL_NULL_AUTOREGRESSIVE_COEFF_95_CONFIDENCE_LIMIT]
           << ": -/+" << standard_normal_value / sqrt(sample_size);
      }

      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << STAT_label[STATL_DETERMINATION_COEFF] << ": " << determination_coeff;
    }

    else {
      os << "   " << STAT_word[STATW_STANDARD_DEVIATION] << " : " << dispersion;
    }
    break;
  }
  }
  os << endl;

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of the characteristics of a continuous distribution.
 *
 *  \param[in,out] os        stream,
 *  \param[in]     file_flag flag comment.
 */
/*--------------------------------------------------------------*/

ostream& ContinuousParametric::ascii_characteristic_print(ostream &os , bool file_flag) const

{
  if (file_flag) {
    os << "# ";
  }

  switch (ident) {

  case GAMMA : {
    if (shape != 0.) {
      double variance = shape * scale * scale;
      boost::math::gamma_distribution<double> dist(shape , scale);

      os << STAT_label[STATL_MEAN] << ": " << shape * scale << "   "
         << STAT_label[STATL_MEDIAN] << ": " << quantile(dist , 0.5) << "   "
         << STAT_label[STATL_MODE] << ": " << mode(dist) << endl;
      if (file_flag) {
        os << "# ";
      }
      os << STAT_label[STATL_VARIANCE] << ": " << variance << "   "
         << STAT_label[STATL_STANDARD_DEVIATION] << ": " << sqrt(variance) << endl;

      if (file_flag) {
        os << "# ";
      }
      os << STAT_label[STATL_SKEWNESS_COEFF] << ": " << 2 / sqrt(shape) << "   "
         << STAT_label[STATL_KURTOSIS_COEFF] << ": " << 6 / shape << endl;
    }

    else {
      os << STAT_label[STATL_MEAN] << ": " << 0 << "   "
         << STAT_label[STATL_VARIANCE] << ": " << 0 << endl;
    }
    break;
  }

  case INVERSE_GAUSSIAN : {
    double variance = location * location * location / scale;
    inverse_gaussian dist(location , scale);

    os << STAT_label[STATL_MEDIAN] << ": " << quantile(dist , 0.5) << "   "
       << STAT_label[STATL_MODE] << ": " << mode(dist) << endl;
    if (file_flag) {
      os << "# ";
    }
    os << STAT_label[STATL_VARIANCE] << ": " << variance << "   "
       << STAT_label[STATL_STANDARD_DEVIATION] << ": " << sqrt(variance) << endl;

    if (file_flag) {
      os << "# ";
    }
    os << STAT_label[STATL_SKEWNESS_COEFF] << ": " << 3 * sqrt(location / scale) << "   "
       << STAT_label[STATL_KURTOSIS_COEFF] << ": " << 15 * location / scale << endl;
    break;
  }

  case GAUSSIAN : {
    os << STAT_label[STATL_VARIANCE] << ": " << dispersion * dispersion << endl;
    break;
  }

  case VON_MISES : {
    double mean_resultant_length , standard_deviation;

    mean_resultant_length = cyl_bessel_i(1 , dispersion) / cyl_bessel_i(0 , dispersion);

    os << STAT_label[STATL_MEAN_RESULTANT_LENGTH] << ": " << mean_resultant_length;

    if (mean_resultant_length > 0.) {
      standard_deviation = sqrt(-2 * log(mean_resultant_length));
      if (unit == DEGREE) {
        standard_deviation *= (180 / M_PI);
      }

      os << "   " << STAT_label[STATL_CIRCULAR_STANDARD_DEVIATION] << ": "
         << standard_deviation;
    }
    os << endl;

#   ifdef MESSAGE
    if (dispersion > 0.) {
      standard_deviation = 1. / sqrt(fabs(dispersion));
      if (unit == DEGREE) {
        standard_deviation *= (180 / M_PI);
      }

      if (file_flag) {
        os << "# ";
      }
      os << STAT_label[STATL_STANDARD_DEVIATION] << ": " << standard_deviation << endl;
    }
#   endif

    break;
  }

  case ZERO_INFLATED_GAMMA : {
    double mean , variance;

    if (zero_probability == 1.) {
      mean = 0.;
      variance = 0.;
    }
    else {
      mean = (1 - zero_probability) * shape * scale;
      variance = (1 - zero_probability) * shape * scale * scale * (1 + shape) - mean * mean;
    }

    os << STAT_label[STATL_MEAN] << ": " << mean << "   "
       << STAT_label[STATL_VARIANCE] << ": " << variance << "   "
       << STAT_label[STATL_STANDARD_DEVIATION] << ": " << sqrt(variance) << endl;
    break;
  }
  }

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a continuous distribution.
 *
 *  \param[in,out] os         stream,
 *  \param[in]     file_flag  flag comment,
 *  \param[in]     cumul_flag flag on the writing of the cumulative distribution function,
 *  \param[in]     histo1     pointer on an Histogram object,
 *  \param[in]     histo2     pointer on a FrequencyDistribution object.
 */
/*--------------------------------------------------------------*/

ostream& ContinuousParametric::ascii_print(ostream &os , bool file_flag ,
                                           bool cumul_flag , const Histogram *histo1 ,
                                           const FrequencyDistribution *histo2)

{
  if ((histo1) || (histo2)) {
    int i , j;
    int nb_step , width[5];
    double step , histo_scale , mass , value , *frequency , *dist_cumul , *histo_cumul;
    ios_base::fmtflags format_flags;


    format_flags = os.setf(ios::right , ios::adjustfield);

    if (histo1) {
      step = histo1->bin_width;
      histo_scale = histo1->nb_element;
    }
    else {
      step = histo2->min_interval_computation();
      histo_scale = histo2->nb_element;
    }

    if (ident == GAMMA) {
      min_value = 0.;
      if (shape == 0.) {
        max_value = 0.;
      }
      else {
        max_value = shape * scale;
      }
    }

    else if (ident == INVERSE_GAUSSIAN) {
      min_value = 0.;
      max_value = location;
    }

    else if ((ident == GAUSSIAN) || (ident == VON_MISES)) {
      min_value = location;

      if (histo1) {
        switch (histo1->type) {

        case INT_VALUE : {
          if (histo1->min_value < min_value) {
            min_value = histo1->min_value;
          }
          break;
        }

        case REAL_VALUE : {
          if (histo1->min_value - histo1->bin_width / 2 < min_value) {
            min_value = histo1->min_value - histo1->bin_width / 2;
          }
          break;
        }
        }
      }

      if ((histo2) && (histo2->offset < min_value)) {
        min_value = histo2->offset;
      }

      max_value = location;
    }

    else if (ident == ZERO_INFLATED_GAMMA) {
      min_value = 0.;
      if (zero_probability == 1.) {
        max_value = 0.;
      }
      else {
        max_value = shape * scale;
      }
    }

    if (histo1) {
      switch (histo1->type) {

      case INT_VALUE : {
        if (histo1->max_value > max_value) {
          max_value = histo1->max_value;
        }
        break;
      }

      case REAL_VALUE : {
        if (histo1->max_value + histo1->bin_width / 2 > max_value) {
          max_value = histo1->max_value + histo1->bin_width / 2;
        }
        break;
      }
      }
    }

    if ((histo2) && (histo2->nb_value - 1 > max_value)) {
      max_value = histo2->nb_value - 1;
    }

    if (cumul_flag) {
      if (histo1) {
        histo_cumul = histo1->cumul_computation();
      }
      else {
        histo_cumul = histo2->cumul_computation();
      }
    }

    // computation of a discretized continuous distribution 

    switch (ident) {

    case GAMMA : {
      if (shape == 0.) {
        dist_cumul = new double[1];
        frequency = new double[1];

        dist_cumul[0] = 1.;
        frequency[0] = histo_scale;
      }

      else {
        boost::math::gamma_distribution<double> dist(shape , scale);

        value = quantile(complement(dist , GAMMA_TAIL));
        while (max_value < value) {
          max_value += step;
        }

        nb_step = (int)((max_value - min_value) / step) + 1;
        dist_cumul = new double[nb_step];
        frequency = new double[nb_step];

//        value = step;
//        dist_cumul[0] = cdf(dist , value);
        value = step / 2;
        dist_cumul[0] = cdf(dist , value + step / 2);
        frequency[0] = dist_cumul[0] * histo_scale;

        for (i = 1;i < nb_step;i++) {
          value += step;
//          dist_cumul[i] = cdf(dist , value);
          dist_cumul[i] = cdf(dist , value + step / 2);
          frequency[i] = (dist_cumul[i] - dist_cumul[i - 1]) * histo_scale;
        }
      }
      break;
    }

    case INVERSE_GAUSSIAN : {
      inverse_gaussian dist(location , scale);

      value = quantile(dist , 1. - INVERSE_GAUSSIAN_TAIL);
      while (max_value < value) {
        max_value += step;
      }

      nb_step = (int)((max_value - min_value) / step) + 1;
      dist_cumul = new double[nb_step];
      frequency = new double[nb_step];

//      value = step;
//      dist_cumul[0] = cdf(dist , value);
      value = step / 2;
      dist_cumul[0] = cdf(dist , value + step / 2);
      frequency[0] = dist_cumul[0] * histo_scale;

      for (i = 1;i < nb_step;i++) {
        value += step;
//        dist_cumul[i] = cdf(dist , value);
        dist_cumul[i] = cdf(dist , value + step / 2);
        frequency[i] = (dist_cumul[i] - dist_cumul[i - 1]) * histo_scale;
      }
      break;
    }

    case GAUSSIAN : {
      normal dist(location , dispersion);

      value = quantile(dist , GAUSSIAN_TAIL);
      while (min_value > value) {
        min_value -= step;
      }

//      value = quantile(dist , 1. - GAUSSIAN_TAIL);
      value = quantile(complement(dist , GAUSSIAN_TAIL));
      while (max_value < value) {
        max_value += step;
      }

      nb_step = (int)((max_value - min_value) / step) + 1;
      dist_cumul = new double[nb_step];
      frequency = new double[nb_step];

      value = min_value;
      dist_cumul[0] = cdf(dist , value + step / 2);
      frequency[0] = (dist_cumul[0] - cdf(dist , value - step / 2)) * histo_scale;

      for (i = 1;i < nb_step;i++) {
        value += step;
        dist_cumul[i] = cdf(dist , value + step / 2);
        frequency[i] = (dist_cumul[i] - dist_cumul[i - 1]) * histo_scale;
      }
      break;
    }

    case VON_MISES : {
      while (min_value >= step / 2) {
        min_value -= step;
      }

      value = (unit == DEGREE ? 360 : 2 * M_PI) - 3 * step / 2;
      while (max_value <= value) {
        max_value += step;
      }

      nb_step = (int)((max_value - min_value) / step) + 1;
      frequency = new double[nb_step];
      dist_cumul = new double[nb_step];

      value = min_value;
      mass = von_mises_mass_computation(value - step / 2 , value + step / 2);
      frequency[0] = mass * histo_scale;
      dist_cumul[0] = mass;

      for (i = 1;i < nb_step;i++) {
        value += step;
        mass = von_mises_mass_computation(value - step / 2 , value + step / 2);
        frequency[i] = mass * histo_scale;
        dist_cumul[i] = dist_cumul[i - 1] + mass;
      }
      break;
    }

    case ZERO_INFLATED_GAMMA : {
      if (zero_probability == 1.) {
        dist_cumul = new double[1];
        frequency = new double[1];

        dist_cumul[0] = 1.;
        frequency[0] = histo_scale;
      }

      else {
        boost::math::gamma_distribution<double> dist(shape , scale);

        value = quantile(complement(dist , GAMMA_TAIL));
        while (max_value < value) {
          max_value += step;
        }

        nb_step = (int)((max_value - min_value) / step) + 1;
        dist_cumul = new double[nb_step];
        frequency = new double[nb_step];

        value = 0.;
        dist_cumul[0] = zero_probability;
        frequency[0] = zero_probability * histo_scale;

        for (i = 1;i < nb_step;i++) {
          value += step;
          dist_cumul[i] = zero_probability + (1. - zero_probability) * cdf(dist , value);
          frequency[i] = (dist_cumul[i] - dist_cumul[i - 1]) * histo_scale;
        }
      }
      break;
    }
    }

    // computation of the column widths

    width[0] = column_width(min_value , max_value);

    if (histo1) {
      width[1] = column_width(histo1->max) + ASCII_SPACE;
    }
    else {
      width[1] = column_width(histo2->max) + ASCII_SPACE;
    }

    width[2] = column_width(nb_step , frequency) + ASCII_SPACE;

    if (cumul_flag) {
      if (histo1) {
        width[3] = column_width(histo1->nb_bin , histo_cumul) + ASCII_SPACE;
      }
      else {
        width[3] = column_width(histo2->nb_value - histo2->offset ,
                                histo_cumul + histo2->offset) + ASCII_SPACE;
      }

      width[4] = column_width(nb_step , dist_cumul) + ASCII_SPACE;
    }

    value = min_value;
    if (histo1) {
      i = -1;
    }
    else {
      i = histo2->offset - step;
    }

    for (j = 0;j < nb_step;j++) {
      if (file_flag) {
        os << "# ";
      }
      os << setw(width[0]) << value;

      if (histo1) {
        if ((value >= histo1->min_value) && (value <= histo1->max_value)) {
          os << setw(width[1]) << histo1->frequency[++i];
        }
        else {
          os << setw(width[1]) << " ";
        }
      }

      else {
        if ((value >= histo2->offset) && (value < histo2->nb_value)) {
          i += step;
          os << setw(width[1]) << histo2->frequency[i];
        }
        else {
          os << setw(width[1]) << " ";
        }
      }

      os << setw(width[2]) << frequency[j];

      if (cumul_flag) {
        if (((histo1) && (value >= histo1->min_value) && (value <= histo1->max_value)) ||
            ((histo2) && (value >= histo2->offset) && (value < histo2->nb_value))) {
          os << setw(width[3]) << histo_cumul[i];
        }
        else {
          os << setw(width[3]) << " ";
        }

        os << setw(width[4]) << dist_cumul[j];
      }

      os << endl;
      value += step;
    }

    delete [] frequency;
    delete [] dist_cumul;

    if (cumul_flag) {
      delete [] histo_cumul;
    }

    os.setf(format_flags , ios::adjustfield);
  }

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of the parameters of a continuous distribution at the spreadsheet format.
 *
 *  \param[in,out] os stream.
 */
/*--------------------------------------------------------------*/

ostream& ContinuousParametric::spreadsheet_parameter_print(ostream &os) const

{
  os << "\n" << STAT_continuous_distribution_word[ident] << "\t";

  switch (ident) {

  case GAMMA : {
    os << STAT_word[STATW_SHAPE] << "\t" << shape;
    if (shape != 0.) {
      os << "\t\t" << STAT_word[STATW_SCALE] << "\t" << scale;
    }
    break;
  }

  case INVERSE_GAUSSIAN : {
    os << STAT_word[STATW_MEAN] << "\t" << location << "\t\t"
       << STAT_word[STATW_SCALE] << "\t" << scale;
    break;
  }

  case GAUSSIAN : {
    os << STAT_word[STATW_MEAN] << "\t" << location << "\t\t"
       << STAT_word[STATW_STANDARD_DEVIATION] << "\t" << dispersion;
    break;
  }

  case VON_MISES : {
    os << STAT_word[STATW_MEAN_DIRECTION] << "\t" << location << "\t\t"
       << STAT_word[STATW_CONCENTRATION] << "\t" << dispersion;
    break;
  }

  case ZERO_INFLATED_GAMMA : {
    os << STAT_word[STATW_ZERO_PROBABILITY] << "\t" << zero_probability;
    if (zero_probability < 1.) {
      os << "\t\t" << STAT_word[STATW_SHAPE] << "\t" << shape
         << "\t\t" << STAT_word[STATW_SCALE] << "\t" << scale;
    }
    break;
  }

  case LINEAR_MODEL : {
    os << STAT_word[STATW_INTERCEPT] << "\t" << intercept << "\t\t"
       << STAT_word[STATW_SLOPE] << "\t" << slope;

    if (sample_size > 0.) {
      Test test(STUDENT , false , sample_size , I_DEFAULT , D_DEFAULT);

      test.critical_probability = ref_critical_probability[0];
      test.t_value_computation();

      if (slope_standard_deviation > 0.) {
        os << "\t" << slope - test.value * slope_standard_deviation << "\t"
           << slope + test.value * slope_standard_deviation;
      }

      os << "\t\t" << STAT_word[STATW_STANDARD_DEVIATION] << "\t" << dispersion << endl;

      os << STAT_label[STATL_CORRELATION_COEFF] << "\t" << correlation
         << "\t\t" << STAT_label[STATL_LIMIT_CORRELATION_COEFF] << "\t" << "-/+"
         << test.value / sqrt(test.value * test.value + sample_size)
         << "\t" << STAT_label[STATL_CRITICAL_PROBABILITY] << "\t" << test.critical_probability;
    }

    else {
      os << "\t\t" << STAT_word[STATW_STANDARD_DEVIATION] << "\t" << dispersion;
    }
    break;
  }

  case AUTOREGRESSIVE_MODEL : {
    os << STAT_word[STATW_MEAN] << "\t" << location << "\t\t"
       << STAT_word[STATW_AUTOREGRESSIVE_COEFF] << "\t" << autoregressive_coeff;

    if (sample_size > 0.) {
      normal dist;
      double standard_normal_value = quantile(complement(dist , 0.025)) , standard_error;

      standard_error = standard_normal_value * sqrt((1. - autoregressive_coeff * autoregressive_coeff) / sample_size);

      os << "\t" << MAX(autoregressive_coeff - standard_error , -1.) << "\t" << MIN(autoregressive_coeff + standard_error , 1.)
         << "\t" << STAT_label[STATL_NULL_AUTOREGRESSIVE_COEFF_95_CONFIDENCE_LIMIT]
         << "\t-/+" << standard_normal_value / sqrt(sample_size);
    }

    os << "\t\t" << STAT_word[STATW_STANDARD_DEVIATION] << "\t" << dispersion << endl;

    os << STAT_label[STATL_DETERMINATION_COEFF] << "\t" << determination_coeff;
    break;
  }
  }
  os << endl;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of the characteristics of a continuous distribution at the spreadsheet format.
 *
 *  \param[in,out] os stream.
 */
/*--------------------------------------------------------------*/

ostream& ContinuousParametric::spreadsheet_characteristic_print(ostream &os) const

{
  switch (ident) {

  case GAMMA : {
    if (shape != 0.) {
      double variance = shape * scale * scale;
      boost::math::gamma_distribution<double> dist(shape , scale);

      os << STAT_label[STATL_MEAN] << "\t" << shape * scale << "\t\t"
         << STAT_label[STATL_MEDIAN] << "\t" << quantile(dist , 0.5) << "\t\t"
         << STAT_label[STATL_MODE] << "\t" << mode(dist) << endl;
      os << STAT_label[STATL_VARIANCE] << "\t" << variance << "\t\t"
         << STAT_label[STATL_STANDARD_DEVIATION] << "\t" << sqrt(variance) << endl;

      os << STAT_label[STATL_SKEWNESS_COEFF] << "\t" << 2 / sqrt(shape) << "\t\t"
         << STAT_label[STATL_KURTOSIS_COEFF] << "\t" << 6 / shape << endl;
    }

    else {
      os << STAT_label[STATL_MEAN] << "\t" << 0 << "\t\t"
         << STAT_label[STATL_VARIANCE] << "\t" << 0 << endl;
    }
    break;
  }

  case INVERSE_GAUSSIAN : {
    double variance = location * location * location / scale;
    inverse_gaussian dist(location , scale);

    os << STAT_label[STATL_MEDIAN] << "\t" << quantile(dist , 0.5) << "\t\t"
       << STAT_label[STATL_MODE] << "\t" << mode(dist) << endl;
    os << STAT_label[STATL_VARIANCE] << "\t" << variance << "\t\t"
       << STAT_label[STATL_STANDARD_DEVIATION] << "\t" << sqrt(variance) << endl;

    os << STAT_label[STATL_SKEWNESS_COEFF] << "\t" << 3 * sqrt(location / scale) << "\t\t"
       << STAT_label[STATL_KURTOSIS_COEFF] << "\t" << 15 * location / scale << endl;
    break;
  }

  case GAUSSIAN : {
    os << STAT_label[STATL_VARIANCE] << "\t" << dispersion * dispersion << endl;
    break;
  }

  case VON_MISES : {
    double mean_resultant_length , standard_deviation;


    mean_resultant_length = cyl_bessel_i(1 , dispersion) / cyl_bessel_i(0 , dispersion);

    os << STAT_label[STATL_MEAN_RESULTANT_LENGTH] << "\t" << mean_resultant_length;

    if (mean_resultant_length > 0.) {
      standard_deviation = sqrt(-2 * log(mean_resultant_length));
      if (unit == DEGREE) {
        standard_deviation *= (180 / M_PI);
      }

      os << "\t" << STAT_label[STATL_CIRCULAR_STANDARD_DEVIATION] << "\t"
         << standard_deviation;
    }
    os << endl;
    break;
  }

  case ZERO_INFLATED_GAMMA : {
    double mean , variance;

    if (zero_probability == 1.) {
      mean = 0.;
      variance = 0.;
    }
    else {
      mean = (1 - zero_probability) * shape * scale;
      variance = (1 - zero_probability) * shape * scale * scale * (1 + shape) - mean * mean;
    }

    os << STAT_label[STATL_MEAN] << "\t" << mean << "\t\t"
       << STAT_label[STATL_VARIANCE] << "\t" << variance << "\t\t"
       << STAT_label[STATL_STANDARD_DEVIATION] << "\t" << sqrt(variance) << endl;
    break;
  }
  }

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a continuous distribution at the spreadsheet format.
 *
 *  \param[in,out] os         stream,
 *  \param[in]     cumul_flag flag on the writing of the cumulative distribution function,
 *  \param[in]     histo1     pointer on an Histogram object,
 *  \param[in]     histo2     pointer on a FrequencyDistribution object.
 */
/*--------------------------------------------------------------*/

ostream& ContinuousParametric::spreadsheet_print(ostream &os , bool cumul_flag ,
                                                 const Histogram *histo1 ,
                                                 const FrequencyDistribution *histo2)

{
  if ((histo1) || (histo2)) {
    int i , j;
    int nb_step;
    double step , histo_scale , value , mass , *frequency , *dist_cumul , *histo_cumul;


    if (histo1) {
      step = histo1->bin_width;
      histo_scale = histo1->nb_element;
    }
    else {
      step = histo2->min_interval_computation();
      histo_scale = histo2->nb_element;
    }

    if (ident == GAMMA) {
      min_value = 0.;
      if (shape == 0.) {
        max_value = 0.;
      }
      else {
        max_value = shape * scale;
      }
    }

    else if (ident == INVERSE_GAUSSIAN) {
      min_value = 0.;
      max_value = location;
    }

    else if ((ident == GAUSSIAN) || (ident == VON_MISES)) {
      min_value = location;

      if (histo1) {
        switch (histo1->type) {

        case INT_VALUE : {
          if (histo1->min_value < min_value) {
            min_value = histo1->min_value;
          }
          break;
        }

        case REAL_VALUE : {
          if (histo1->min_value - histo1->bin_width / 2 < min_value) {
            min_value = histo1->min_value - histo1->bin_width / 2;
          }
          break;
        }
        }
      }

      if ((histo2) && (histo2->offset < min_value)) {
        min_value = histo2->offset;
      }

      max_value = location;
    }

    else if (ident == ZERO_INFLATED_GAMMA) {
      min_value = 0.;
      if (zero_probability == 1.) {
        max_value = 0.;
      }
      else {
        max_value = shape * scale;
      }
    }

    if (histo1) {
      switch (histo1->type) {

      case INT_VALUE : {
        if (histo1->max_value > max_value) {
          max_value = histo1->max_value;
        }
        break;
      }

      case REAL_VALUE : {
        if (histo1->max_value + histo1->bin_width / 2 > max_value) {
          max_value = histo1->max_value + histo1->bin_width / 2;
        }
        break;
      }
      }
    }

    if ((histo2) && (histo2->nb_value - 1 > max_value)) {
      max_value = histo2->nb_value - 1;
    }

    if (cumul_flag) {
      if (histo1) {
        histo_cumul = histo1->cumul_computation();
      }
      else {
        histo_cumul = histo2->cumul_computation();
      }
    }

    // computation of a discretized continuous distribution

    switch (ident) {

    case GAMMA : {
      if (shape == 0.) {
        dist_cumul = new double[1];
        frequency = new double[1];

        dist_cumul[0] = 1.;
        frequency[0] = histo_scale;
      }

      else {
        boost::math::gamma_distribution<double> dist(shape , scale);

        value = quantile(complement(dist , GAMMA_TAIL));
        while (max_value < value) {
          max_value += step;
        }

        nb_step = (int)((max_value - min_value) / step) + 1;
        dist_cumul = new double[nb_step];
        frequency = new double[nb_step];

//        value = step;
//        dist_cumul[0] = cdf(dist , value);
        value = step / 2;
        dist_cumul[0] = cdf(dist , value + step / 2);
        frequency[0] = dist_cumul[0] * histo_scale;

        for (i = 1;i < nb_step;i++) {
          value += step;
//          dist_cumul[i] = cdf(dist , value);
          dist_cumul[i] = cdf(dist , value + step / 2);
          frequency[i] = (dist_cumul[i] - dist_cumul[i - 1]) * histo_scale;
        }
      }
      break;
    }

    case INVERSE_GAUSSIAN : {
      inverse_gaussian dist(location , scale);

      value = quantile(dist , 1. - INVERSE_GAUSSIAN_TAIL);
      while (max_value < value) {
        max_value += step;
      }

      nb_step = (int)((max_value - min_value) / step) + 1;
      dist_cumul = new double[nb_step];
      frequency = new double[nb_step];

//      value = step;
//      dist_cumul[0] = cdf(dist , value);
      value = step / 2;
      dist_cumul[0] = cdf(dist , value + step / 2);
      frequency[0] = dist_cumul[0] * histo_scale;

      for (i = 1;i < nb_step;i++) {
        value += step;
//        dist_cumul[i] = cdf(dist , value);
        dist_cumul[i] = cdf(dist , value + step / 2);
        frequency[i] = (dist_cumul[i] - dist_cumul[i - 1]) * histo_scale;
      }
      break;
    }

    case GAUSSIAN : {
      normal dist(location , dispersion);

      value = quantile(dist , GAUSSIAN_TAIL);
      while (min_value > value) {
        min_value -= step;
      }

      value = quantile(complement(dist , GAUSSIAN_TAIL));
      while (max_value < value) {
        max_value += step;
      }

      nb_step = (int)((max_value - min_value) / step) + 1;
      dist_cumul = new double[nb_step];
      frequency = new double[nb_step];

      value = min_value;
      dist_cumul[0] = cdf(dist , value + step / 2);
      frequency[0] = (dist_cumul[0] - cdf(dist , value - step / 2)) * histo_scale;

      for (i = 1;i < nb_step;i++) {
        value += step;
        dist_cumul[i] = cdf(dist , value + step / 2);
        frequency[i] = (dist_cumul[i] - dist_cumul[i - 1]) * histo_scale;
      }
      break;
    }

    case VON_MISES : {
      while (min_value >= step / 2) {
        min_value -= step;
      }

      value = (unit == DEGREE ? 360 : 2 * M_PI) - 3 * step / 2;
      while (max_value <= value) {
        max_value += step;
      }

      nb_step = (int)((max_value - min_value) / step) + 1;
      frequency = new double[nb_step];
      dist_cumul = new double[nb_step];

      value = min_value;
      mass = von_mises_mass_computation(value - step / 2 , value + step / 2);
      frequency[0] = mass * histo_scale;
      dist_cumul[0] = mass;

      for (i = 1;i < nb_step;i++) {
        value += step;
        mass = von_mises_mass_computation(value - step / 2 , value + step / 2);
        frequency[i] = mass * histo_scale;
        dist_cumul[i] = dist_cumul[i - 1] + mass;
      }
      break;
    }

    case ZERO_INFLATED_GAMMA : {
      if (zero_probability == 1.) {
        dist_cumul = new double[1];
        frequency = new double[1];

        dist_cumul[0] = 1.;
        frequency[0] = histo_scale;
      }

      else {
        boost::math::gamma_distribution<double> dist(shape , scale);

        value = quantile(complement(dist , GAMMA_TAIL));
        while (max_value < value) {
          max_value += step;
        }

        nb_step = (int)((max_value - min_value) / step) + 1;
        dist_cumul = new double[nb_step];
        frequency = new double[nb_step];

        value = 0.;
        dist_cumul[0] = zero_probability;
        frequency[0] = zero_probability * histo_scale;

        for (i = 1;i < nb_step;i++) {
          value += step;
          dist_cumul[i] = zero_probability + (1. - zero_probability) * cdf(dist , value);
          frequency[i] = (dist_cumul[i] - dist_cumul[i - 1]) * histo_scale;
        }
      }
      break;
    }
    }

    value = min_value;
    if (histo1) {
      i = -1;
    }
    else {
      i = histo2->offset - step;
    }

    for (j = 0;j < nb_step;j++) {
      os << value << "\t";

      if ((histo1) && (value >= histo1->min_value) && (value <= histo1->max_value)) {
        os << histo1->frequency[++i];
      }

      if ((histo2) && (value >= histo2->offset) && (value < histo2->nb_value)) {
        i += step;
        os << histo2->frequency[i];
      }

      os << "\t" << frequency[j];

      if (cumul_flag) {
        os << "\t";
        if (((histo1) && (value >= histo1->min_value) && (value <= histo1->max_value)) ||
            ((histo2) && (value >= histo2->offset) && (value < histo2->nb_value))) {
          os << histo_cumul[i];
        }

        os << "\t" << dist_cumul[j];
      }

      os << endl;
      value += step;
    }

    delete [] frequency;
    delete [] dist_cumul;

    if (cumul_flag) {
      delete [] histo_cumul;
    }
  }

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of the parameters of a continuous distribution at the Gnuplot format.
 *
 *  \param[in,out] os stream.
 */
/*--------------------------------------------------------------*/

ostream& ContinuousParametric::plot_title_print(ostream &os) const

{
  os << STAT_continuous_distribution_letter[ident] << "(";

  if (ident == GAMMA) {
    os << shape;
    if (shape != 0.) {
      os << ", " << scale;
    }
  }

  else if (ident == INVERSE_GAUSSIAN) {
    os << location << ", " << scale;
  }

  else if (ident == ZERO_INFLATED_GAMMA) {
    os << zero_probability;
    if (zero_probability < 1.) {
      os << ", " << shape << ", " << scale;
    }
  }

  else if (ident == LINEAR_MODEL) {
    os << intercept << ", " << slope << ", " << dispersion;
  }

  else if (ident == AUTOREGRESSIVE_MODEL) {
    os << location << ", " << autoregressive_coeff << ", " << dispersion;
  }

  else {
    os << location << ", " << dispersion;
  }
  os << ")";

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a continuous distribution at the Gnuplot format.
 *
 *  \param[in] path   file path,
 *  \param[in] histo1 pointer on an Histogram object,
 *  \param[in] histo2 pointer on a FrequencyDistribution object.
 *
 *  \return           error status.
 */
/*--------------------------------------------------------------*/

bool ContinuousParametric::plot_print(const char *path , const Histogram *histo1 ,
                                      const FrequencyDistribution *histo2)

{
  bool status = false;
  int i;
  int nb_step;
  double histo_scale , step , buff , value , max , *frequency;
  ofstream out_file(path);


  if (out_file) {
    status = true;

    if (histo1) {
      histo_scale = histo1->bin_width * histo1->nb_element;
    }
    else if (histo2) {
      histo_scale = histo2->nb_element;
    }
    else {
      histo_scale = 1.;
    }

    switch (ident) {

    case GAMMA : {
      min_value = 0.;

      if (shape == 0.) {
        max_value = 0.;
        frequency = new double[1];
        frequency[0] = histo_scale;
        max = frequency[0];
      }

      else {
        boost::math::gamma_distribution<double> dist(shape , scale);

        max_value = quantile(complement(dist , GAMMA_TAIL));
        if ((histo1) && (histo1->max_value + histo1->bin_width > max_value)) {
          max_value = histo1->max_value + histo1->bin_width;
        }
        if ((histo2) && (histo2->nb_value - 1)) {
          max_value = histo2->nb_value - 1;
        }

        step = max_value / GAMMA_NB_STEP;
        if (histo1) {
          buff = histo1->bin_width / GAMMA_NB_SUB_STEP;
        }
        else if (histo2) {
          buff = histo2->min_interval_computation() / GAMMA_NB_SUB_STEP;
        }
        if (((histo1) || (histo2)) && (buff < step)) {
          step = buff;
        }

        nb_step = (int)(max_value / step) + 1;
        frequency = new double[nb_step];

        value = 0.;
        max = 0.;
        for (i = 0;i < nb_step;i++) {
          frequency[i] = pdf(dist , value) * histo_scale;
          value += step;
          if (frequency[i] > max) {
            max = frequency[i];
          }
        }
      }
      break;
    }

    case INVERSE_GAUSSIAN : {
      inverse_gaussian dist(location , scale);

      min_value = 0.;
      max_value = quantile(dist , 1. - INVERSE_GAUSSIAN_TAIL);
      if ((histo1) && (histo1->max_value + histo1->bin_width > max_value)) {
        max_value = histo1->max_value + histo1->bin_width;
      }
      if ((histo2) && (histo2->nb_value - 1)) {
        max_value = histo2->nb_value - 1;
      }

      step = max_value / INVERSE_GAUSSIAN_NB_STEP;
      if (histo1) {
        buff = histo1->bin_width / INVERSE_GAUSSIAN_NB_SUB_STEP;
      }
      else if (histo2) {
        buff = histo2->min_interval_computation() / INVERSE_GAUSSIAN_NB_SUB_STEP;
      }
      if (((histo1) || (histo2)) && (buff < step)) {
        step = buff;
      }

      nb_step = (int)(max_value / step) + 1;
      frequency = new double[nb_step];

      value = 0.;
      max = 0.;
      for (i = 0;i < nb_step;i++) {
        frequency[i] = pdf(dist , value) * histo_scale;
        value += step;
        if (frequency[i] > max) {
          max = frequency[i];
        }
      }
      break;
    }

    case GAUSSIAN : {
      normal dist(location , dispersion);

      min_value = quantile(dist , GAUSSIAN_TAIL);
      if ((histo1) && (histo1->min_value - histo1->bin_width < min_value)) {
        min_value = histo1->min_value - histo1->bin_width;
      }
      if ((histo2) && (histo2->offset < min_value)) {
        min_value = histo2->offset;
      }

      max_value = quantile(complement(dist , GAUSSIAN_TAIL));
      if ((histo1) && (histo1->max_value + histo1->bin_width > max_value)) {
        max_value = histo1->max_value + histo1->bin_width;
      }
      if ((histo2) && (histo2->nb_value - 1)) {
        max_value = histo2->nb_value - 1;
      }

      step = (max_value - min_value) / GAUSSIAN_NB_STEP;
      if (histo1) {
        buff = histo1->bin_width / GAUSSIAN_NB_SUB_STEP;
      }
      else if (histo2) {
        buff = histo2->min_interval_computation() / GAUSSIAN_NB_SUB_STEP;
      }
      if (((histo1) || (histo2)) && (buff < step)) {
        step = buff;
      }

/*     computation complet

      nb_step = (int)((max_value - min_value) / step) + 1;
      frequency = new double[nb_step];

      value = min_value;
      for (i = 0;i < nb_step;i++) {
        frequency[i] = pdf(dist , value) * histo_scale;
        value += step;
      }

      max = frequency[(int)round((location - min_value) / step)]; */

      nb_step = (int)(MAX(max_value - location , location - min_value) / step) + 1;
      frequency = new double[nb_step * 2 - 1];

      value = location;
      i = 0;
      while (value < MAX(max_value , 2 * location - min_value)) {
        frequency[i++] = pdf(dist , value) * histo_scale;
        value += step;
      }

      max = frequency[0];

#     ifdef DEBUG
      cout << "\nNB_STEP: " << i << " | " << nb_step << endl; 
#     endif

      nb_step = i;
      buff = frequency[nb_step - 1];
      for (i = 0;i < nb_step - 1;i++) {
        frequency[nb_step - 1 + i] = frequency[i];
      }
      frequency[2 * (nb_step - 1)] = buff;
      for (i = 1;i < nb_step;i++) {
        frequency[nb_step - 1 - i] = frequency[nb_step - 1 + i];
      }

      value = location - (nb_step - 1) * step;

#     ifdef DEBUG
      cout << "\nMIN_VALUE: " << value << " | " << min_value << endl; 
#     endif

      nb_step = 2 * nb_step - 1;
      break;
    }

    case VON_MISES : {
      frequency = new double[VON_MISES_NB_STEP];
      min_value = 0.;

      switch (unit) {

      case DEGREE : {
        max_value = 360;
        step = max_value / VON_MISES_NB_STEP;

        value = 0.;
        for (i = 0;i < VON_MISES_NB_STEP;i++) {
          frequency[i] = exp(dispersion * cos((value - location) * M_PI / 180)) * histo_scale /
                         (360 * cyl_bessel_i(0 , dispersion));
          value += step;
        }
        break;
      }

      case RADIAN : {
        max_value = 2 * M_PI;
        step = max_value / VON_MISES_NB_STEP;

        value = 0.;
        for (i = 0;i < VON_MISES_NB_STEP;i++) {
          frequency[i] = exp(dispersion * cos(value - location)) * histo_scale /
                         (2 * M_PI * cyl_bessel_i(0 , dispersion));
          value += step;
        }
        break;
      }
      }

      max = frequency[(int)round(location / step)];

      nb_step = VON_MISES_NB_STEP;
      value = 0.;
      break;
    }

    case ZERO_INFLATED_GAMMA : {
      min_value = 0.;

      if (zero_probability == 1.) {
        max_value = 0.;
        frequency = new double[1];
        frequency[0] = histo_scale;
        max = frequency[0];
      }

      else {
        boost::math::gamma_distribution<double> dist(shape , scale);

        max_value = quantile(complement(dist , GAMMA_TAIL));
        if ((histo1) && (histo1->max_value + histo1->bin_width > max_value)) {
          max_value = histo1->max_value + histo1->bin_width;
        }
        if ((histo2) && (histo2->nb_value - 1)) {
          max_value = histo2->nb_value - 1;
        }

        step = max_value / GAMMA_NB_STEP;
        if (histo1) {
          buff = histo1->bin_width / GAMMA_NB_SUB_STEP;
        }
        else if (histo2) {
          buff = histo2->min_interval_computation() / GAMMA_NB_SUB_STEP;
        }
        if (((histo1) || (histo2)) && (buff < step)) {
          step = buff;
        }

        nb_step = (int)(max_value / step) + 1;
        frequency = new double[nb_step];

        value = 0.;
        frequency[0] = zero_probability * histo_scale;
        max = frequency[0];

        for (i = 1;i < nb_step;i++) {
          value += step;
          frequency[i] = (1. - zero_probability) * pdf(dist , value) * histo_scale;
          if (frequency[i] > max) {
            max = frequency[i];
          }
        }
      }
      break;
    }
    }

    for (i = 0;i < nb_step;i++) {
      out_file << value << " " << frequency[i] << endl;
      value += step;
    }

    delete [] frequency;
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a q-q plot at the Gnuplot format.
 *
 *  \param[in] path          file path,
 *  \param[in] nb_value      number of values,
 *  \param[in] empirical_cdf pointer on the empirical cumulative distribution function.
 *
 *  \return                  error status.
 */
/*--------------------------------------------------------------*/

bool ContinuousParametric::q_q_plot_print(const char *path , int nb_value ,
                                          double **empirical_cdf) const

{
  bool status = false;
  int i;
  double **qqplot;
  ofstream out_file(path);


  if (out_file) {
    status = true;

    qqplot = q_q_plot_computation(nb_value , empirical_cdf);

    for (i = 0;i < (ident == VON_MISES ? nb_value : nb_value - 1);i++) {
      out_file << qqplot[0][i] << " " << qqplot[1][i] << endl;
    }

    delete [] qqplot[0];
    delete [] qqplot[1];
    delete [] qqplot;
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Plot of a continuous distribution.
 *
 *  \param[in] plot   reference on a SinglePlot object,
 *  \param[in] histo1 pointer on an Histogram object,
 *  \param[in] histo2 pointer on a FrequencyDistribution object.
 */
/*--------------------------------------------------------------*/

void ContinuousParametric::plotable_write(SinglePlot &plot , const Histogram *histo1 ,
                                          const FrequencyDistribution *histo2)

{
  int i;
  int nb_step;
  double histo_scale , step , buff , value , max , *frequency;


  if (histo1) {
    histo_scale = histo1->bin_width * histo1->nb_element;
  }
  else if (histo2) {
    histo_scale = histo2->nb_element;
  }
  else {
    histo_scale = 1.;
  }

  switch (ident) {

  case GAMMA : {
    min_value = 0.;

    if (shape == 0.) {
      max_value = 0.;
      frequency = new double[1];
      frequency[0] = histo_scale;
      max = frequency[0];
    }

    else {
      boost::math::gamma_distribution<double> dist(shape , scale);

      max_value = quantile(complement(dist , GAMMA_TAIL));
      if ((histo1) && (histo1->max_value + histo1->bin_width > max_value)) {
        max_value = histo1->max_value + histo1->bin_width;
      }
      if ((histo2) && (histo2->nb_value - 1)) {
        max_value = histo2->nb_value - 1;
      }

      step = max_value / GAMMA_NB_STEP;
      if (histo1) {
        buff = histo1->bin_width / GAMMA_NB_SUB_STEP;
      }
      else if (histo2) {
        buff = histo2->min_interval_computation() / GAMMA_NB_SUB_STEP;
      }
      if (((histo1) || (histo2)) && (buff < step)) {
        step = buff;
      }

      nb_step = (int)(max_value / step) + 1;
      frequency = new double[nb_step];

      value = 0.;
      max = 0.;
      for (i = 0;i < nb_step;i++) {
        frequency[i] = pdf(dist , value) * histo_scale;
        value += step;

        if (frequency[i] > max) {
          max = frequency[i];
        }
      }
    }
    break;
  }

  case INVERSE_GAUSSIAN : {
    inverse_gaussian dist(location , scale);

    min_value = 0.;
    max_value = quantile(dist , 1. - INVERSE_GAUSSIAN_TAIL);
    if ((histo1) && (histo1->max_value + histo1->bin_width > max_value)) {
      max_value = histo1->max_value + histo1->bin_width;
    }
    if ((histo2) && (histo2->nb_value - 1)) {
      max_value = histo2->nb_value - 1;
    }

    step = max_value / INVERSE_GAUSSIAN_NB_STEP;
    if (histo1) {
      buff = histo1->bin_width / INVERSE_GAUSSIAN_NB_SUB_STEP;
    }
    else if (histo2) {
      buff = histo2->min_interval_computation() / INVERSE_GAUSSIAN_NB_SUB_STEP;
    }
    if (((histo1) || (histo2)) && (buff < step)) {
      step = buff;
    }

    nb_step = (int)(max_value / step) + 1;
    frequency = new double[nb_step];

    value = 0.;
    max = 0.;
    for (i = 0;i < nb_step;i++) {
      frequency[i] = pdf(dist , value) * histo_scale;
      value += step;
      if (frequency[i] > max) {
        max = frequency[i];
      }
    }
    break;
  }

  case GAUSSIAN : {
    normal dist(location , dispersion);

    min_value = quantile(dist , GAUSSIAN_TAIL);
    if ((histo1) && (histo1->min_value - histo1->bin_width < min_value)) {
      min_value = histo1->min_value - histo1->bin_width;
    }
    if ((histo2) && (histo2->offset < min_value)) {
      min_value = histo2->offset;
    }

    max_value = quantile(complement(dist , GAUSSIAN_TAIL));
    if ((histo1) && (histo1->max_value + histo1->bin_width > max_value)) {
      max_value = histo1->max_value + histo1->bin_width;
    }
    if ((histo2) && (histo2->nb_value - 1)) {
      max_value = histo2->nb_value - 1;
    }

    step = (max_value - min_value) / GAUSSIAN_NB_STEP;
    if (histo1) {
      buff = histo1->bin_width / GAUSSIAN_NB_SUB_STEP;
    }
    else if (histo2) {
      buff = histo2->min_interval_computation() / GAUSSIAN_NB_SUB_STEP;
    }
    if (((histo1) || (histo2)) && (buff < step)) {
      step = buff;
    }

    nb_step = (int)(MAX(max_value - location , location - min_value) / step) + 1;
    frequency = new double[nb_step * 2 - 1];

    value = location;
    i = 0;
    while (value < MAX(max_value , 2 * location - min_value)) {
      frequency[i++] = pdf(dist , value) * histo_scale;
      value += step;
    }

    max = frequency[0];

#   ifdef DEBUG
    cout << "\nNB_STEP: " << i << " | " << nb_step << endl; 
#   endif

    nb_step = i;
    buff = frequency[nb_step - 1];
    for (i = 0;i < nb_step - 1;i++) {
      frequency[nb_step - 1 + i] = frequency[i];
    }
    frequency[2 * (nb_step - 1)] = buff;
    for (i = 1;i < nb_step;i++) {
      frequency[nb_step - 1 - i] = frequency[nb_step - 1 + i];
    }

    value = location - (nb_step - 1) * step;

#   ifdef DEBUG
    cout << "\nMIN_VALUE: " << value << " | " << min_value << endl; 
#   endif

    nb_step = 2 * nb_step - 1;
    break;
  }

  case VON_MISES : {
    frequency = new double[VON_MISES_NB_STEP];
    min_value = 0.;

    switch (unit) {

    case DEGREE : {
      max_value = 360;
      step = max_value / VON_MISES_NB_STEP;

      value = 0.;
      for (i = 0;i < VON_MISES_NB_STEP;i++) {
        frequency[i] = exp(dispersion * cos((value - location) * M_PI / 180)) * histo_scale /
                       (360 * cyl_bessel_i(0 , dispersion));
        value += step;
      }
      break;
    }

    case RADIAN : {
      max_value = 2 * M_PI;
      step = max_value / VON_MISES_NB_STEP;

      value = 0.;
      for (i = 0;i < VON_MISES_NB_STEP;i++) {
        frequency[i] = exp(dispersion * cos(value - location)) * histo_scale /
                       (2 * M_PI * cyl_bessel_i(0 , dispersion));
        value += step;
      }
      break;
    }
    }

    max = frequency[(int)round(location / step)];

    nb_step = VON_MISES_NB_STEP;
    value = 0.;
    break;
  }

  case ZERO_INFLATED_GAMMA : {
    min_value = 0.;

    if (zero_probability == 1.) {
      max_value = 0.;
      frequency = new double[1];
      frequency[0] = histo_scale;
      max = frequency[0];
    }

    else {
      boost::math::gamma_distribution<double> dist(shape , scale);

      max_value = quantile(complement(dist , GAMMA_TAIL));
      if ((histo1) && (histo1->max_value + histo1->bin_width > max_value)) {
        max_value = histo1->max_value + histo1->bin_width;
      }
      if ((histo2) && (histo2->nb_value - 1)) {
        max_value = histo2->nb_value - 1;
      }

      step = max_value / GAMMA_NB_STEP;
      if (histo1) {
        buff = histo1->bin_width / GAMMA_NB_SUB_STEP;
      }
      else if (histo2) {
        buff = histo2->min_interval_computation() / GAMMA_NB_SUB_STEP;
      }
      if (((histo1) || (histo2)) && (buff < step)) {
        step = buff;
      }

      nb_step = (int)(max_value / step) + 1;
      frequency = new double[nb_step];

      value = 0.;
      frequency[0] = zero_probability * histo_scale;
      max = frequency[0];

      for (i = 1;i < nb_step;i++) {
        value += step;
        frequency[i] = (1. - zero_probability) * pdf(dist , value) * histo_scale;

        if (frequency[i] > max) {
          max = frequency[i];
        }
      }
    }
    break;
  }
  }

  for (i = 0;i < nb_step;i++) {
    plot.add_point(value , frequency[i]);
    value += step;
  }

  delete [] frequency;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Plot of a q-q plot.
 *
 *  \param[in] plot          reference on a SinglePlot object,
 *  \param[in] nb_value      number of values,
 *  \param[in] empirical_cdf pointer on the empirical cumulative distribution function.
 */
/*--------------------------------------------------------------*/

void ContinuousParametric::q_q_plotable_write(SinglePlot &plot , int nb_value ,
                                              double **empirical_cdf) const

{
  int i;
  double **qqplot;


  qqplot = q_q_plot_computation(nb_value , empirical_cdf);

  for (i = 0;i < (ident == VON_MISES ? nb_value : nb_value - 1);i++) {
    plot.add_point(qqplot[0][i] , qqplot[1][i]);
  }

  delete [] qqplot[0];
  delete [] qqplot[1];
  delete [] qqplot;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the number of parameters of a continuous distribution.
 *
 *  \return number of parameters.
 */
/*--------------------------------------------------------------*/

int ContinuousParametric::nb_parameter_computation() const

{
  int nb_parameter;


  switch (ident) {

  case GAMMA : {
    if (shape == 0.) {
      nb_parameter = 1;
    }
    else {
      nb_parameter = 2;
    }
    break;
  }

  case INVERSE_GAUSSIAN : {
    nb_parameter = 2;
    break;
  }

  case GAUSSIAN : {
    nb_parameter = 2;
    break;
  }

  case VON_MISES : {
    nb_parameter = 2;
    break;
  }

  case ZERO_INFLATED_GAMMA : {
    if (zero_probability == 1.) {
      nb_parameter = 1;
    }
    else if (zero_probability == 0.) {
      nb_parameter = 2;
    }
    else {
      nb_parameter = 3;
    }
    break;
  }

  case LINEAR_MODEL : {
    nb_parameter = 3;
    break;
  }

  case AUTOREGRESSIVE_MODEL : {
    nb_parameter = 3;
    break;
  }

  default : {
    nb_parameter = 0;
    break;
  }
  }

  return nb_parameter;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the density of a von Mises distribution integrated on an interval.
 *
 *  \param[in] inf lower bound,
 *  \param[in] sup upper bound.
 *
 *  \return        integrated density.
 */
/*--------------------------------------------------------------*/

double ContinuousParametric::von_mises_mass_computation(double inf , double sup) const

{
  double step , mass , value;


//  if (dispersion < VON_MISES_CONCENTRATION_THRESHOLD) {
    mass = 0.;

    switch (unit) {

    case DEGREE : {
      step = MIN((sup - inf) / VON_MISES_NB_SUB_STEP , 360. / VON_MISES_NB_STEP);

      for (value = inf + step / 2;value < sup;value += step) {
        mass += exp(dispersion * cos((value - location) * M_PI / 180));
      }
      mass = mass * step / (360 * cyl_bessel_i(0 , dispersion));
      break;
    }

    case RADIAN : {
      step = MIN((sup - inf) / VON_MISES_NB_SUB_STEP , 2 * M_PI / VON_MISES_NB_STEP);

      for (value = inf + step;value < sup;value += step) {
        mass += exp(dispersion * cos(value - location));
      }
      mass = mass * step / (2 * M_PI * cyl_bessel_i(0 , dispersion));
      break;
    }
    }
/*  }

  else {
    normal dist(location , sqrt(1. / dispersion));

    mass = cdf(dist , sup) - cdf(dist , inf);
  } */

  return mass;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the density of a distribution integrated on an interval.
 *
 *  \param[in] inf lower bound,
 *  \param[in] sup upper bound.
 *
 *  \return        integrated density.
 */
/*--------------------------------------------------------------*/

double ContinuousParametric::mass_computation(double inf , double sup) const

{
  double mass;


  switch (ident) {

  case GAMMA : {
    if (inf < 0.) {
      inf = 0.;
    }
    if (sup < inf) {
      sup = inf;
    }

    if (shape == 0.) {
      if (inf == 0.) {
        mass = 1.;
      }
      else {
        mass = 0.;
      }
    }

    else {
      boost::math::gamma_distribution<double> dist(shape , scale);

      if (inf == sup) {
//         mass = pdf(dist , MAX(inf , CONTINUOUS_POSITIVE_INF_BOUND));  bug boost C++
        mass = pdf(dist , inf);

#       ifdef DEBUG
        if ((shape == 1.) && (scale < 0.1) && (inf == 0.)) {
          cout << "TEST: " << scale << ", " << mass << endl;
        }
#       endif

      }

      else {
        mass = cdf(dist , sup) - cdf(dist , inf);
      }
    }
    break;
  }

  case INVERSE_GAUSSIAN : {
    inverse_gaussian dist(location , scale);

    if (inf < 0.) {
      inf = 0.;
    }
    if (sup < inf) {
      sup = inf;
    }

    if (inf == sup) {
      mass = pdf(dist , inf);
    }
    else {
      mass = cdf(dist , sup) - cdf(dist , inf);
    }
    break;
  }

  case GAUSSIAN : {
    normal dist(location , dispersion);

    if (inf == sup) {
      mass = pdf(dist , inf);
    }
    else {
      mass = cdf(dist , sup) - cdf(dist , inf);
    }
    break;
  }

  case VON_MISES : {
    if (inf == sup) {
      switch (unit) {
      case DEGREE :
        mass = exp(dispersion * cos((inf - location) * M_PI / 180)) /
               (360 * cyl_bessel_i(0 , dispersion));
        break;
      case RADIAN :
        mass = exp(dispersion * cos(inf - location)) /
               (2 * M_PI * cyl_bessel_i(0 , dispersion));
        break;
      }
    }

    else {
      mass = von_mises_mass_computation(inf , sup);
    }
    break;
  }

  case ZERO_INFLATED_GAMMA : {
    if (inf < 0.) {
      inf = 0.;
    }
    if (sup < inf) {
      sup = inf;
    }

    if (zero_probability == 1.) {
      if (inf == 0.) {
        mass = zero_probability;
      }
      else {
        mass = 0.;
      }
    }

    else {
      boost::math::gamma_distribution<double> dist(shape , scale);

      if (inf == 0.) {
        mass = zero_probability;
      }

      else {
        if (inf == sup) {
          mass = (1. - zero_probability) * pdf(dist , inf);
        }
        else {
          mass = (1. - zero_probability) * (cdf(dist , sup) - cdf(dist , inf));
        }
      }

      if (mass > 1.) {

#       ifdef MESSAGE
        cout << "WARNING: " << inf << " " << sup << " " << mass << " | "
             << zero_probability << " " << shape << " " << scale << endl;
#       endif
 
      }
    }
    break;
  }

  case LINEAR_MODEL : {
    normal dist(0. , dispersion);

    if (inf == sup) {
      mass = pdf(dist , inf);
    }
    else {
      mass = cdf(dist , sup) - cdf(dist , inf);
    }
    break;
  }

  case AUTOREGRESSIVE_MODEL : {
    normal dist(0. , dispersion);

    if (inf == sup) {
      mass = pdf(dist , inf);
    }
    else {
      mass = cdf(dist , sup) - cdf(dist , inf);
    }
    break;
  }
  }

//  if (mass > 1.) {
  if ((inf == sup) && (mass > 1.)) {
    mass = 1.;
  }

  return mass;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the cumulative distribution function of a von Mises distribution.
 */
/*--------------------------------------------------------------*/

void ContinuousParametric::von_mises_cumul_computation()

{
  int i , j;
  int start;
  double step , value;


  cumul = new double[VON_MISES_NB_STEP];

  switch (unit) {

  case DEGREE : {
    step = 360. / VON_MISES_NB_STEP;

    if (location < 180) {
      value = location + 180;
    }
    else {
      value = location - 180;
    }
    start = (int)round(value * VON_MISES_NB_STEP / 360);

    cumul[start] = exp(dispersion * cos((value - location) * M_PI / 180)) * step /
                   (360 * cyl_bessel_i(0 , dispersion));
    for (i = start + 1;i < VON_MISES_NB_STEP;i++) {
      value += step;
      cumul[i] = cumul[i - 1] + exp(dispersion * cos((value - location) * M_PI / 180)) * step /
                                (360 * cyl_bessel_i(0 , dispersion));
    }

    value += step - 360;
    cumul[0] = cumul[VON_MISES_NB_STEP - 1] + exp(dispersion * cos((value - location) * M_PI / 180)) * step /
                                              (360 * cyl_bessel_i(0 , dispersion));
    for (i = 1;i < start;i++) {
      value += step;
      cumul[i] = cumul[i - 1] + exp(dispersion * cos((value - location) * M_PI / 180)) * step /
                                (360 * cyl_bessel_i(0 , dispersion));
    }
    break;
  }

  case RADIAN : {
    step = 2 * M_PI / VON_MISES_NB_STEP;

    if (location < M_PI) {
      value = location + M_PI;
    }
    else {
      value = location - M_PI;
    }
    start = (int)round(value * VON_MISES_NB_STEP / (2 * M_PI));

    cumul[start] = exp(dispersion * cos(value - location)) * step /
                   (2 * M_PI * cyl_bessel_i(0 , dispersion));
    for (i = start + 1;i < VON_MISES_NB_STEP;i++) {
      value += step;
      cumul[i] = cumul[i - 1] + exp(dispersion * cos(value - location)) * step /
                                (2 * M_PI * cyl_bessel_i(0 , dispersion));
    }

    value += step - 2 * M_PI;
    cumul[0] = cumul[VON_MISES_NB_STEP - 1] + exp(dispersion * cos(value - location)) * step /
                                              (2 * M_PI * cyl_bessel_i(0 , dispersion));
    for (i = 1;i < start;i++) {
      value += step;
      cumul[i] = cumul[i - 1] + exp(dispersion * cos(value - location)) * step /
                                (2 * M_PI * cyl_bessel_i(0 , dispersion));
    }
    break;
  }
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of a q-q plot.
 *
 *  \param[in] nb_value      number of values,
 *  \param[in] empirical_cdf pointer on the empirical cumulative distribution function.
 *
 *  \return                  q-q plot.
 */
/*--------------------------------------------------------------*/

double** ContinuousParametric::q_q_plot_computation(int nb_value ,
                                                    double **empirical_cdf) const

{
  int i , j;
  double step , value , current_cumul , previous_cumul , **qqplot;


  qqplot = new double*[2];
  qqplot[0] = new double[nb_value];
  qqplot[1] = new double[nb_value];

  for (i = 0;i < nb_value;i++) {
    qqplot[0][i] = empirical_cdf[0][i];
  }

  switch (ident) {

  case GAMMA : {
    if (shape == 0.) {
      for (i = 0;i < nb_value - 1;i++) {
        qqplot[1][i] = 0.;
      }
    }

    else {
      boost::math::gamma_distribution<double> dist(shape , scale);

      for (i = 0;i < nb_value - 1;i++) {
        qqplot[1][i] = quantile(dist , empirical_cdf[1][i]);
      }
    }
    break;
  }

  case INVERSE_GAUSSIAN : {
    inverse_gaussian dist(location , scale);

    for (i = 0;i < nb_value - 1;i++) {
      qqplot[1][i] = quantile(dist , empirical_cdf[1][i]);
    }
    break;
  }

  case GAUSSIAN : {
    normal dist(location , dispersion);

    for (i = 0;i < nb_value - 1;i++) {
      qqplot[1][i] = quantile(dist , empirical_cdf[1][i]);
    }
    break;
  }

  case VON_MISES : {
    current_cumul = 0.;
    value = 0.;

    switch (unit) {

    case DEGREE : {
      step = 360. / VON_MISES_NB_STEP;

      for (i = 0;i < nb_value - 1;i++) {
        while (current_cumul < empirical_cdf[1][i]) {
          value += step;
          previous_cumul = current_cumul;
          current_cumul += exp(dispersion * cos((value - step / 2 - location) * M_PI / 180)) * step /
                           (360 * cyl_bessel_i(0 , dispersion));
        }

        qqplot[1][i] = value - step * (current_cumul - empirical_cdf[1][i]) /
                       (current_cumul - previous_cumul);
      }

      qqplot[1][nb_value - 1] = 360;
      break;
    }

    case RADIAN : {
      step = 2 * M_PI / VON_MISES_NB_STEP;

      for (i = 0;i < nb_value - 1;i++) {
        while (current_cumul < empirical_cdf[1][i]) {
          value += step;
          previous_cumul = current_cumul;
          current_cumul += exp(dispersion * cos(value - step / 2 - location)) * step /
                           (2 * M_PI * cyl_bessel_i(0 , dispersion));
        }

        qqplot[1][i] = value - step * (current_cumul - empirical_cdf[1][i]) /
                       (current_cumul - previous_cumul);
      }

      qqplot[1][nb_value - 1] = 2 * M_PI;
      break;
    }
    }
    break;
  }

  case ZERO_INFLATED_GAMMA : {
    if (zero_probability == 1.) {
      for (i = 0;i < nb_value - 1;i++) {
        qqplot[1][i] = 0.;
      }
    }

    else {
      boost::math::gamma_distribution<double> dist(shape , scale);

      i = 0;
      while (empirical_cdf[1][i] <= zero_probability) {
        qqplot[1][i++] = 0.;
      }
      for (j = i;j < nb_value - 1;i++) {
        qqplot[1][j] = quantile(dist , (empirical_cdf[1][j] - zero_probability) / (1. - zero_probability));
      }
    }
    break;
  }
  }

  return qqplot;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the distance between 2 continuous distributions
 *         (sup of the absolute difference between the cumulative distribution functions
 *          in the case of non-crossing cumulative distribution functions; in the general case,
 *          sum of sup on intervals between 2 crossings of cumulative distribution functions).
 *
 *  \param[in] dist reference on a ContinuousParametric object.
 *
 *  \return         distance between 2 continuous distributions
 */
/*--------------------------------------------------------------*/

double ContinuousParametric::sup_norm_distance_computation(ContinuousParametric &dist)

{
  bool crossing;
  int i;
  double min , max , step , value , buff , distance , max_absolute_diff , cumul1[2] , cumul2[2];

# ifdef MESSAGE
  double overlap;
# endif


  switch (ident) {

  case GAMMA : {
    if ((shape == 0.) || (dist.shape == 0.)) {
      if ((shape == 0.) && (dist.shape == 0.)) {
        distance = 0.;
      }
      else {
        distance = 1.;
      }
    }

    else {
      boost::math::gamma_distribution<double> dist1(shape , scale) , dist2(dist.shape , dist.scale);

#     ifdef MESSAGE
      if ((quantile(complement(dist1 , GAMMA_TAIL)) < quantile(dist1 , 1. - GAMMA_TAIL) - DOUBLE_ERROR) ||
          (quantile(complement(dist1 , GAMMA_TAIL)) > quantile(dist1 , 1. - GAMMA_TAIL) + DOUBLE_ERROR)) {
        cout << "\nERROR: ";
        ascii_parameter_print(cout);
        cout << quantile(complement(dist1 , GAMMA_TAIL)) << " | "
             << quantile(dist1 , 1. - GAMMA_TAIL) << endl;
      }
#     endif

      step = MIN(quantile(complement(dist1 , GAMMA_TAIL)) , quantile(complement(dist2 , GAMMA_TAIL))) / GAMMA_NB_STEP;
      distance = 0.;
      value = CONTINUOUS_POSITIVE_INF_BOUND;
      cumul1[0] = cdf(dist1 , value);
      cumul2[0] = cdf(dist2 , value);
      max_absolute_diff = fabs(cumul1[0] - cumul2[0]);
      crossing = false;
      value = 0.;

      for (i = 1;i < GAMMA_NB_STEP;i++) {
        value += step;
        cumul1[1] = cdf(dist1 , value);
        cumul2[1] = cdf(dist2 , value);

        buff = fabs(cumul1[1] - cumul2[1]);
        if (buff > max_absolute_diff) {
          max_absolute_diff = buff;
        }

        if ((!crossing) && (((cumul1[1] > cumul2[1]) && (cumul1[0] <= cumul2[0])) ||
            ((cumul1[1] <= cumul2[1]) && (cumul1[0] > cumul2[0])))) {
          crossing = true;
          distance = max_absolute_diff;
          max_absolute_diff = 0.;
        }

        cumul1[0] = cumul1[1];
        cumul2[0] = cumul2[1];
      }
      distance += max_absolute_diff;

#     ifdef DEBUG
      value = 0.;
      cumul1[0] = 0.;
      cumul2[0] = 0.;
      overlap = 0.;

      for (i = 0;i < GAMMA_NB_STEP;i++) {
        cumul1[1] = cdf(dist1 , value);
        cumul2[1] = cdf(dist2 , value);
        value += step;

        overlap += MIN(cumul1[1] - cumul1[0] , cumul2[1] - cumul2[0]);
        cumul1[0] = cumul1[1];
        cumul2[0] = cumul2[1];
      }

      cout << "\nSup norm distance: " << distance << " " << 1. - overlap << endl;
#     endif

    }
    break;
  }

  case INVERSE_GAUSSIAN : {
    inverse_gaussian dist1(location , scale) , dist2(dist.location , dist.scale);

#   ifdef DEBUG
    if ((quantile(complement(dist1 , INVERSE_GAUSSIAN_TAIL)) < quantile(dist1 , 1. - INVERSE_GAUSSIAN_TAIL) - DOUBLE_ERROR) ||
        (quantile(complement(dist1 , INVERSE_GAUSSIAN_TAIL)) > quantile(dist1 , 1. - INVERSE_GAUSSIAN_TAIL) + DOUBLE_ERROR)) {
      cout << "\nERROR: ";
      ascii_parameter_print(cout);
      cout << quantile(complement(dist1 , INVERSE_GAUSSIAN_TAIL)) << " | "
           << quantile(dist1 , 1. - INVERSE_GAUSSIAN_TAIL) << endl;
    }
#   endif

    step = MIN(quantile(dist1 , 1. - INVERSE_GAUSSIAN_TAIL) , quantile(dist2 , 1. - INVERSE_GAUSSIAN_TAIL)) /
           INVERSE_GAUSSIAN_NB_STEP;
    distance = 0.;
    value = 0.;
    cumul1[0] = 0.;
    cumul2[0] = 0.;
    max_absolute_diff = 0.;
    crossing = false;

    for (i = 1;i < INVERSE_GAUSSIAN_NB_STEP;i++) {
      value += step;
      cumul1[1] = cdf(dist1 , value);
      cumul2[1] = cdf(dist2 , value);

      buff = fabs(cumul1[1] - cumul2[1]);
      if (buff > max_absolute_diff) {
        max_absolute_diff = buff;
      }

      if ((!crossing) && (((cumul1[1] > cumul2[1]) && (cumul1[0] <= cumul2[0])) ||
          ((cumul1[1] <= cumul2[1]) && (cumul1[0] > cumul2[0])))) {
        crossing = true;
        distance = max_absolute_diff;
        max_absolute_diff = 0.;
      }

      cumul1[0] = cumul1[1];
      cumul2[0] = cumul2[1];
    }
    distance += max_absolute_diff;

#   ifdef DEBUG
    value = 0.;
    cumul1[0] = 0.;
    cumul2[0] = 0.;
    overlap = 0.;

    for (i = 1;i < INVERSE_GAUSSIAN_NB_STEP;i++) {
      value += step;
      cumul1[1] = cdf(dist1 , value);
      cumul2[1] = cdf(dist2 , value);

      overlap += MIN(cumul1[1] - cumul1[0] , cumul2[1] - cumul2[0]);
      cumul1[0] = cumul1[1];
      cumul2[0] = cumul2[1];
    }

    cout << "\n" << STAT_label[STATL_SUP_NORM_DISTANCE] << ": " << distance << " " << 1. - overlap << endl;
#   endif

    break;
  }

  case GAUSSIAN : {
    normal dist1(location , dispersion) , dist2(dist.location , dist.dispersion);

    min = MAX(quantile(dist1 , GAUSSIAN_TAIL) , quantile(dist2 , GAUSSIAN_TAIL));
    max = MIN(quantile(complement(dist1 , GAUSSIAN_TAIL)) , quantile(complement(dist2 , GAUSSIAN_TAIL)));

#   ifdef MESSAGE
    if ((quantile(complement(dist1 , GAUSSIAN_TAIL)) < quantile(dist1 , 1. - GAUSSIAN_TAIL) - DOUBLE_ERROR) ||
        (quantile(complement(dist1 , GAUSSIAN_TAIL)) > quantile(dist1 , 1. - GAUSSIAN_TAIL) + DOUBLE_ERROR)) {
      cout << "\nERROR: ";
      ascii_parameter_print(cout);
      cout << quantile(complement(dist1 , GAUSSIAN_TAIL)) << " | "
           << quantile(dist1 , 1. - GAUSSIAN_TAIL) << endl;
    }
#   endif

    step = (max - min) / GAUSSIAN_NB_STEP;
    distance = 0.;
    value = min;
    cumul1[0] = cdf(dist1 , value);
    cumul2[0] = cdf(dist2 , value);
    max_absolute_diff = fabs(cumul1[0] - cumul2[0]);
    crossing = false;

    for (i = 1;i < GAUSSIAN_NB_STEP;i++) {
      value += step;
      cumul1[1] = cdf(dist1 , value);
      cumul2[1] = cdf(dist2 , value);

      buff = fabs(cumul1[1] - cumul2[1]);
      if (buff > max_absolute_diff) {
        max_absolute_diff = buff;
      }

      if ((!crossing) && (((cumul1[1] > cumul2[1]) && (cumul1[0] <= cumul2[0])) ||
          ((cumul1[1] <= cumul2[1]) && (cumul1[0] > cumul2[0])))) {
        crossing = true;
        distance = max_absolute_diff;
        max_absolute_diff = 0.;
      }

      cumul1[0] = cumul1[1];
      cumul2[0] = cumul2[1];
    }
    distance += max_absolute_diff;

#   ifdef MESSAGE
    value = min;
    cumul1[0] = 0.;
    cumul2[0] = 0.;
    overlap = 0.;

    for (i = 0;i < GAUSSIAN_NB_STEP;i++) {
      cumul1[1] = cdf(dist1 , value);
      cumul2[1] = cdf(dist2 , value);
      value += step;

      overlap += MIN(cumul1[1] - cumul1[0] , cumul2[1] - cumul2[0]);
      cumul1[0] = cumul1[1];
      cumul2[0] = cumul2[1];
    }

    cout << "\n" << STAT_label[STATL_SUP_NORM_DISTANCE] << ": " << distance << " " << 1. - overlap << endl;
#   endif

    break;
  }

  case VON_MISES : {   // bug small concentration
    int inf , sup;

    if (!cumul) {
      von_mises_cumul_computation();
    }
    if (!dist.cumul) {
      dist.von_mises_cumul_computation();
    }

    switch (unit) {

    case DEGREE : {
      if (location < dist.location) {
        if (dist.location - location <= 180) {
          min = dist.location - 180;
          max = location + 180;
        }
        else {
          min = location - 180;
          max = dist.location + 180;
        }
      }

      else {
        if (location - dist.location <= 180) {
          min = location - 180;
          max = dist.location + 180;
        }
        else {
          min = dist.location - 180;
          max = location + 180;
        }
      }

      if (min < 0) {
        min += 360;
      }
      if (max >= 360) {
        max -= 360;
      }

      inf = (int)round(min * VON_MISES_NB_STEP / 360);
      sup = (int)round(max * VON_MISES_NB_STEP / 360);
      break;
    }

    case RADIAN : {
      if (location < dist.location) {
        if (dist.location - location <= M_PI) {
          min = dist.location - M_PI;
          max = location + M_PI;
        }
        else {
          min = location - M_PI;
          max = dist.location + M_PI;
        }
      }

      else {
        if (location - dist.location <= M_PI) {
          min = location - M_PI;
          max = dist.location + M_PI;
        }
        else {
          min = dist.location - M_PI;
          max = location + M_PI;
        }
      }

      if (min < 0) {
        min += 2 * M_PI;
      }
      if (max >= 2 * M_PI) {
        max -= 2 * M_PI;
      }

      inf = (int)round(min * VON_MISES_NB_STEP / (2 * M_PI));
      sup = (int)round(max * VON_MISES_NB_STEP / (2 * M_PI));
      break;
    }
    }

#   ifdef DEBUG
    cout << "\nInf: " << inf << ",  Sup: " << sup << endl;
#   endif

    distance = 0.;
    max_absolute_diff = fabs(cumul[inf] - dist.cumul[inf]);
    crossing = false;

    if (inf <= sup) {
      for (i = inf + 1;i < sup;i++) {
        buff = fabs(cumul[i] - dist.cumul[i]);
        if (buff > max_absolute_diff) {
          max_absolute_diff = buff;
        }

        if ((!crossing) && (((cumul[i] > dist.cumul[i]) && (cumul[i - 1] <= dist.cumul[i - 1])) ||
            ((cumul[i] <= dist.cumul[i]) && (cumul[i - 1] > dist.cumul[i - 1])))) {

#         ifdef DEBUG
          cout << "\nCrossing: " << i << endl;
#         endif

          crossing = true;
          distance = max_absolute_diff;
          max_absolute_diff = 0.;
        }
      }
    }

    else {
      for (i = inf + 1;i < VON_MISES_NB_STEP;i++) {
        buff = fabs(cumul[i] - dist.cumul[i]);
        if (buff > max_absolute_diff) {
          max_absolute_diff = buff;
        }

        if ((!crossing) && (((cumul[i] > dist.cumul[i]) && (cumul[i - 1] <= dist.cumul[i - 1])) ||
            ((cumul[i] <= dist.cumul[i]) && (cumul[i - 1] > dist.cumul[i - 1])))) {
          crossing = true;
          distance = max_absolute_diff;
          max_absolute_diff = 0.;
        }
      }

      buff = fabs(cumul[0] - dist.cumul[0]);
      if (buff > max_absolute_diff) {
        max_absolute_diff = buff;
      }

      if ((!crossing) && (((cumul[0] > dist.cumul[0]) && (cumul[VON_MISES_NB_STEP - 1] <= dist.cumul[VON_MISES_NB_STEP - 1])) ||
          ((cumul[0] <= dist.cumul[0]) && (cumul[VON_MISES_NB_STEP - 1] > dist.cumul[VON_MISES_NB_STEP - 1])))) {
        crossing = true;
        distance = max_absolute_diff;
        max_absolute_diff = 0.;
      }

      for (i = 1;i < sup;i++) {
        buff = fabs(cumul[i] - dist.cumul[i]);
        if (buff > max_absolute_diff) {
          max_absolute_diff = buff;
        }

        if ((!crossing) && (((cumul[i] > dist.cumul[i]) && (cumul[i - 1] <= dist.cumul[i - 1])) ||
            ((cumul[i] <= dist.cumul[i]) && (cumul[i - 1] > dist.cumul[i - 1])))) {
          crossing = true;
          distance = max_absolute_diff;
          max_absolute_diff = 0.;
        }
      }
    }

    distance += max_absolute_diff;

#   ifdef MESSAGE
    step = (unit == DEGREE ? 360. : 2 * M_PI) / VON_MISES_NB_STEP;
    value = 0.;
    overlap = 0.;

    for (i = 0;i < VON_MISES_NB_STEP;i++) {
      overlap += MIN(von_mises_mass_computation(value , value + step) , dist.von_mises_mass_computation(value , value + step));
      value += step;
    }

    cout << "\nSup norm distance: " << distance << " " << 1. - overlap << endl;
#   endif

    break;
  }

  case ZERO_INFLATED_GAMMA : {
    cumul1[0] = zero_probability;
    cumul2[0] = dist.zero_probability;
    max_absolute_diff = fabs(cumul1[0] - cumul2[0]);

    if ((zero_probability == 1.) || (dist.zero_probability == 1.)) {
      distance = max_absolute_diff;

#     ifdef MESSAGE
      cout << "\nSup norm distance: " << distance << " " << 1. - MIN(zero_probability , dist.zero_probability) << endl;
#     endif

    }

    else {
      boost::math::gamma_distribution<double> dist1(shape , scale) , dist2(dist.shape , dist.scale);

      step = MIN(quantile(complement(dist1 , GAMMA_TAIL)) , quantile(complement(dist2 , GAMMA_TAIL))) / GAMMA_NB_STEP;
      value = 0.;
      distance = 0.;
      crossing = false;

      for (i = 1;i < GAMMA_NB_STEP;i++) {
        value += step;
        cumul1[1] = zero_probability + (1. - zero_probability) * cdf(dist1 , value);
        cumul2[1] = dist.zero_probability + (1. - dist.zero_probability) * cdf(dist2 , value);

        buff = fabs(cumul1[1] - cumul2[1]);
        if (buff > max_absolute_diff) {
          max_absolute_diff = buff;
        }

        if ((!crossing) && (((cumul1[1] > cumul2[1]) && (cumul1[0] <= cumul2[0])) ||
            ((cumul1[1] <= cumul2[1]) && (cumul1[0] > cumul2[0])))) {
          crossing = true;
          distance = max_absolute_diff;
          max_absolute_diff = 0.;
        }

        cumul1[0] = cumul1[1];
        cumul2[0] = cumul2[1];
      }
      distance += max_absolute_diff;

#     ifdef MESSAGE
      value = 0.;
      cumul1[0] = zero_probability;
      cumul2[0] = dist.zero_probability;
      overlap = MIN(zero_probability , dist.zero_probability);

      for (i = 1;i < GAMMA_NB_STEP;i++) {
        value += step;
        cumul1[1] = zero_probability + (1. - zero_probability) * cdf(dist1 , value);
        cumul2[1] = dist.zero_probability + (1. - dist.zero_probability) * cdf(dist2 , value);

        overlap += MIN(cumul1[1] - cumul1[0] , cumul2[1] - cumul2[0]);
        cumul1[0] = cumul1[1];
        cumul2[0] = cumul2[1];
      }

      cout << "\nSup norm distance: " << distance << " " << 1. - overlap << endl;
#     endif

    }

    break;
  }
  }

  return distance;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the log-likelihood of a continuous distribution for a sample.
 *
 *  \param[in] variable variable index,
 *  \param[in] dist     pointer on a ContinuousParametric object.
 *
 *  \return             log-likelihood.
 */
/*--------------------------------------------------------------*/

double Vectors::likelihood_computation(int variable , const ContinuousParametric &dist) const

{
  int i;
  double mass , likelihood;


  if (marginal_distribution[variable]) {
    likelihood = marginal_distribution[variable]->likelihood_computation(dist , min_interval[variable]);
  }

  else {
    likelihood = 0.;

    if (((dist.ident == GAMMA) || (dist.ident == INVERSE_GAUSSIAN)) && (min_value[variable] < min_interval[variable] / 2)) {
      switch (type[variable]) {

      case INT_VALUE : {
        for (i = 0;i < nb_vector;i++) {
          mass = dist.mass_computation(int_vector[i][variable] , int_vector[i][variable] + min_interval[variable]);
          if (mass > 0.) {
            likelihood += log(mass);
          }
          else {
            likelihood = D_INF;
            break;
          }
        }
        break;
      }

      case REAL_VALUE : {
        for (i = 0;i < nb_vector;i++) {
          mass = dist.mass_computation(real_vector[i][variable] , real_vector[i][variable] + min_interval[variable]);
          if (mass > 0.) {
            likelihood += log(mass);
          }
          else {
            likelihood = D_INF;
            break;
          }
        }
        break;
      }
      }
    }

    else {
      switch (type[variable]) {

      case INT_VALUE : {
        for (i = 0;i < nb_vector;i++) {
          mass = dist.mass_computation(int_vector[i][variable] - min_interval[variable] / 2 , int_vector[i][variable] + min_interval[variable] / 2);
          if (mass > 0.) {
            likelihood += log(mass);
          }
          else {
            likelihood = D_INF;
            break;
          }
        }
        break;
      }

      case REAL_VALUE : {
        for (i = 0;i < nb_vector;i++) {
          mass = dist.mass_computation(real_vector[i][variable] - min_interval[variable] / 2 , real_vector[i][variable] + min_interval[variable] / 2);
          if (mass > 0.) {
            likelihood += log(mass);
          }
          else {
            likelihood = D_INF;
            break;
          }
        }
        break;
      }
      }
    }
  }

  return likelihood;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Estimation of the parameters of a gamma distribution.
 *
 *  \param[in] variable variable index,
 *  \param[in] dist     pointer on a ContinuousParametric object.
 *
 *  \return             log-likelihood of the estimated distribution.
 */
/*--------------------------------------------------------------*/

double Vectors::gamma_estimation(int variable , ContinuousParametric *dist) const

{
  int i;
  int zero_count;
  double buff , log_geometric_mean , diff , likelihood = D_INF;
//  double variance;


  zero_count = 0;

  switch (type[variable]) {

  case INT_VALUE : {
    for (i = 0;i < nb_vector;i++) {
      if (int_vector[i][variable] == 0) {
        zero_count++;
      }
    }
    break;
  }

  case REAL_VALUE : {
    for (i = 0;i < nb_vector;i++) {
      if (real_vector[i][variable] == 0.) {
        zero_count++;
      }
    }
    break;
  }
  }

  if (zero_count / nb_vector > GAMMA_ZERO_FREQUENCY_THRESHOLD) {
    dist->shape = 0;
    dist->scale = D_DEFAULT;
    likelihood = 0.;
  }

  else {
    if (covariance[variable][variable] > 0.) {
/*      if (sqrt(covariance[variable][variable]) < mean[variable] * GAMMA_VARIATION_COEFF_THRESHOLD) {
        variance = mean[variable] * mean[variable] * GAMMA_VARIATION_COEFF_THRESHOLD * GAMMA_VARIATION_COEFF_THRESHOLD;
      }
      else {
        variance = covariance[variable][variable];
      }

      dist->shape = mean[variable] * mean[variable] / variance;
      dist->scale = variance / mean[variable]; */

      // Hawang & Huang (2012), Ann. Inst. Statist. Math. 54(4), 840-847

      buff = mean[variable] * mean[variable] / covariance[variable][variable];
      if (buff > GAMMA_INVERSE_SAMPLE_SIZE_FACTOR / (double)nb_vector) {
        dist->shape = buff - 1. / (double)nb_vector;
      }
      else {
        dist->shape = buff;
      }
/*      if (dist->shape < GAMMA_MIN_SHAPE_PARAMETER) {
        dist->shape = GAMMA_MIN_SHAPE_PARAMETER;
      } */
      dist->scale = mean[variable] / dist->shape;

      if ((dist->shape >= GAMMA_SHAPE_PARAMETER_THRESHOLD) && (nb_vector < GAMMA_FREQUENCY_THRESHOLD)) {
        log_geometric_mean = 0.;

        switch (type[variable]) {

        case INT_VALUE : {
          for (i = 0;i < nb_vector;i++) {
            if (int_vector[i][variable] > 0) {
              log_geometric_mean += log(int_vector[i][variable]);
            }
          }
          break;
        }

        case REAL_VALUE : {
          for (i = 0;i < nb_vector;i++) {
            if (real_vector[i][variable] > 0.) {
              log_geometric_mean += log(real_vector[i][variable]);
            }
          }
          break;
        }
        }

        log_geometric_mean /= (nb_vector - zero_count);
/*        i = 0;   to be reworked

#       ifdef DEBUG
        cout << "\n" << STAT_word[STATW_SHAPE] << " : " << dist->shape << "   "
             << STAT_word[STATW_SCALE] << " : " << dist->scale << endl;
#       endif

        do {
          dist->scale = exp(log_geometric_mean - digamma(dist->shape));
          dist->shape = mean[variable] / dist->scale;
          i++;

#         ifdef DEBUG
          cout << STAT_word[STATW_SHAPE] << " : " << dist->shape << "   "
               << STAT_word[STATW_SCALE] << " : " << dist->scale << endl;
#         endif

        }
        while (i < MIN(GAMMA_ITERATION_FACTOR * iter , GAMMA_MAX_NB_ITERATION)); */

        // approximations Johnson, Kotz & Balakrishnan, Continuous Univariate Distributions, vol. 1, 2nd ed., pp. 361-362

//        dist->shape = mean[variable] / (2 * (mean[variable] - exp(log_geometric_mean))) - 1./12.;
        diff = log(mean[variable]) - log_geometric_mean;
        dist->shape = (1 + sqrt(1 + 4 * diff / 3)) / (4 * diff);
        dist->scale = mean[variable] / dist->shape;
      }

#     ifdef MESSAGE
      cout << "\n" << STAT_word[STATW_SHAPE] << " : " << dist->shape << "   "
           << STAT_word[STATW_SCALE] << " : " << dist->scale << endl;
#     endif

      likelihood = likelihood_computation(variable , *dist);
    }

    else {
      dist->shape = GAMMA_MIN_SHAPE_PARAMETER;
      dist->scale = GAMMA_DEFAULT_SCALE_PARAMETER;
      likelihood = D_INF;
    }
  }

  return likelihood;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Estimation of the parameters of an inverse Gaussian distribution.
 *
 *  \param[in] variable variable index,
 *  \param[in] dist     pointer on a ContinuousParametric object.
 *
 *  \return             log-likelihood of the estimated distribution.
 */
/*--------------------------------------------------------------*/

double Vectors::inverse_gaussian_estimation(int variable , ContinuousParametric *dist) const

{
  int i;
  double inverse_scale , likelihood = D_INF;


  if (mean[variable] > 0.) {
    dist->location = mean[variable];

    inverse_scale = 0.;

    switch (type[variable]) {

    case INT_VALUE : {
      for (i = 0;i < nb_vector;i++) {
        if (int_vector[i][variable] > 0) {
          inverse_scale += 1. / (double)int_vector[i][variable] - 1. / mean[variable];
        }
      }
      break;
    }

    case REAL_VALUE : {
      for (i = 0;i < nb_vector;i++) {
        if (real_vector[i][variable] > 0.) {
          inverse_scale += 1. / real_vector[i][variable] - 1. / mean[variable];
        }
      }
      break;
    }
    }

    if (inverse_scale > 0.) {
      dist->scale = nb_vector / inverse_scale;
      likelihood = likelihood_computation(variable , *dist);
    }
    else {
      dist->scale = D_DEFAULT;
      likelihood = D_INF;
    }
  }

  return likelihood;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Estimation of the parameters of a Gaussian distribution.
 *
 *  \param[in] variable variable index,
 *  \param[in] dist     pointer on a ContinuousParametric object.
 *
 *  \return             log-likelihood of the estimated distribution.
 */
/*--------------------------------------------------------------*/

double Vectors::gaussian_estimation(int variable , ContinuousParametric *dist) const

{
  double likelihood = D_INF;


  if (covariance[variable][variable] > 0.) {
    dist->location = mean[variable];
    dist->dispersion = sqrt(covariance[variable][variable]);

    if (dist->dispersion / dist->location < GAUSSIAN_MIN_VARIATION_COEFF) {
      dist->dispersion = dist->location * GAUSSIAN_MIN_VARIATION_COEFF;
    }

    likelihood = likelihood_computation(variable , *dist);
  }

  return likelihood;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Estimation of the parameters of a von Mises distribution.
 *
 *  \param[in] variable variable index,
 *  \param[in] dist     pointer on a ContinuousParametric object.
 *
 *  \return             log-likelihood of the estimated distribution.
 */
/*--------------------------------------------------------------*/

double Vectors::von_mises_estimation(int variable , ContinuousParametric *dist) const

{
  int i;
  double *mean_direction , likelihood = D_INF;


  mean_direction = new double[4];
  mean_direction[0] = 0.;
  mean_direction[1] = 0.;

  switch (type[variable]) {

  case INT_VALUE : {
    for (i = 0;i < nb_vector;i++) {
      mean_direction[0] += cos(int_vector[i][variable] * M_PI / 180);
      mean_direction[1] += sin(int_vector[i][variable] * M_PI / 180);
    }
    break;
  }

  case REAL_VALUE : {
    for (i = 0;i < nb_vector;i++) {
      switch (dist->unit) {
      case DEGREE :
        mean_direction[0] += cos(real_vector[i][variable] * M_PI / 180);
        mean_direction[1] += sin(real_vector[i][variable] * M_PI / 180);
        break;
      case RADIAN :
        mean_direction[0] += cos(real_vector[i][variable]);
        mean_direction[1] += sin(real_vector[i][variable]);
        break;
      }
    }
    break;
  }
  }

  mean_direction[0] /= nb_vector;
  mean_direction[1] /= nb_vector;

  mean_direction[2] = sqrt(mean_direction[0] * mean_direction[0] + mean_direction[1] * mean_direction[1]);

  if (mean_direction[2] > 0.) {
    mean_direction[3] = atan(mean_direction[1] / mean_direction[0]);

    if (mean_direction[0] < 0.) {
      mean_direction[3] += M_PI;
    }
    else if (mean_direction[1] < 0.) {
      mean_direction[3] += 2 * M_PI;
    }

    if (dist->unit == DEGREE) {
      mean_direction[3] *= (180 / M_PI);
    }

    dist->location = mean_direction[3];
    dist->dispersion = von_mises_concentration_computation(mean_direction[2]);

    likelihood = likelihood_computation(variable , *dist);
  }

  else {
    dist->location = D_INF;
    likelihood = D_INF;
  }

  delete [] mean_direction;

  return likelihood;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Estimation of the parameters of a continuous parametric distribution.
 *
 *  \param[in] variable variable index,
 *  \param[in] ident    identifier of a continuous parametric distribution.
 *
 *  \return             log-likelihood of the estimated distribution.
 */
/*--------------------------------------------------------------*/

double Vectors::continuous_parametric_estimation(int variable , continuous_parametric ident) const

{
  int nb_parameter;
  double likelihood;
  ContinuousParametric *dist;


  dist = new ContinuousParametric(ident);

  switch (ident) {
  case GAMMA :
    likelihood = gamma_estimation(variable , dist);
    break;
  case INVERSE_GAUSSIAN :
    likelihood = inverse_gaussian_estimation(variable , dist);
    break;
  case GAUSSIAN :
    likelihood = gaussian_estimation(variable , dist);
    break;
  case VON_MISES :
    likelihood = von_mises_estimation(variable , dist);
    break;
  }

  if (likelihood != D_INF) {
    nb_parameter = dist->nb_parameter_computation();

#   ifdef MESSAGE
    cout << "\n" << STAT_label[STATL_LIKELIHOOD] << ": " << likelihood << "   ("
         << STAT_label[STATL_NORMALIZED] << ": " << likelihood / nb_vector << ")"
         << "\n" << nb_parameter << " " << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS]
         << "   2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[BIC] << "): "
         << 2 * likelihood - nb_parameter * log((double)nb_vector) << endl;
#   endif

  }

  delete dist;

  return likelihood;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Simulation using a continuous distribution.
 *
 *  \return generated value.
 */
/*--------------------------------------------------------------*/

double ContinuousParametric::simulation()

{
  int i;
  int start;
  double limit , value , step , current_cumul , previous_cumul;


  limit = ((double)rand() / (RAND_MAX + 1.));

  switch (ident) {

  case GAMMA : {
    if (shape == 0.) {
      value = 0.;
    }

    else {
      boost::math::gamma_distribution<double> dist(shape , scale);

      value = quantile(dist , limit);
    }
    break;
  }

  case INVERSE_GAUSSIAN : {
    inverse_gaussian dist(location , scale);

    value = quantile(dist , limit);
    break;
  }

  case GAUSSIAN : {
    normal dist(location , dispersion);

    value = quantile(dist , limit);
    break;
  }

  case VON_MISES : {
    if (!cumul) {
      von_mises_cumul_computation();
    }

    switch (unit) {

    case DEGREE : {
      step = 360. / VON_MISES_NB_STEP;

      if (location < 180) {
        value = location + 180.;
      }
      else {
        value = location - 180.;
      }

      start = (int)round(value * VON_MISES_NB_STEP / 360);
      break;
    }

    case RADIAN : {
      step = 2 * M_PI / VON_MISES_NB_STEP;

      if (location < M_PI) {
        value = location + M_PI;
      }
      else {
        value = location - M_PI;
      }

      start = (int)round(value * VON_MISES_NB_STEP / (2 * M_PI));
      break;
    }
    }

    for (i = start;i < VON_MISES_NB_STEP;i++) {
      if (cumul[i] >= limit) {
        if (i > start) {
          value -= step * (cumul[i] - limit) / (cumul[i] - cumul[i - 1]);
        }
        break;
      }

      else {
        value += step;
      }
    }

    if (i == VON_MISES_NB_STEP) {
      for (i = 0;i < start;i++) {
        if (cumul[i] >= limit) {
          if (i == 0) {
            value -= step * (cumul[i] - limit) / (cumul[i] - cumul[VON_MISES_NB_STEP - 1]);
          }
          else {
            value -= step * (cumul[i] - limit) / (cumul[i] - cumul[i - 1]);
          }
          break;
        }

        else {
          value += step;
        }
      }
    }

    switch (unit) {

    case DEGREE : {
      if (value > 360) {
        value -= 360;
      }
      break;
    }

    case RADIAN : {
      if (value > 2 * M_PI) {
        value -= 2 * M_PI;
      }
      break;
    }
    }
    break;
  }

  case ZERO_INFLATED_GAMMA : {
    if (limit <= zero_probability) {
      value = 0.;
    }

    else if (zero_probability < 1.) {
      boost::math::gamma_distribution<double> dist(shape , scale);

      value = quantile(dist , (limit - zero_probability) / (1. - zero_probability));
    }
    break;
  }

  case LINEAR_MODEL : {
    normal dist(0. , dispersion);

    value = quantile(dist , limit);
    break;
  }

  case AUTOREGRESSIVE_MODEL : {
    normal dist(0. , dispersion);

    value = quantile(dist , limit);
    break;
  }
  }

  return value;
}


};  // namespace stat_tool
