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

#include <boost/tokenizer.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/classification.hpp>

#include "distribution.h"
#include "stat_label.h"

using namespace std;
using namespace boost;


namespace stat_tool {


extern bool cumul_matching_plot_print(const char *path , int nb_cumul , int *offset ,
                                      int *nb_value , double **cumul);



/*--------------------------------------------------------------*/
/**
 *  \brief Initialization of the parameters of a discrete parametric distribution.
 *
 *  \param[in] iinf_bound   lower bound,
 *  \param[in] isup_bound   upper bound (binomial, uniform),
 *  \param[in] iparameter   parameter (Poisson, negative binomial, Poisson geometric),
 *  \param[in] iprobability probability (binomial, negative binomial, Poisson geometric).
 */
/*--------------------------------------------------------------*/

void DiscreteParametric::init(int iinf_bound , int isup_bound ,
                              double iparameter , double iprobability)

{
  inf_bound = iinf_bound;
  sup_bound = isup_bound;
  parameter = iparameter;
  probability = iprobability;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Initialization of the identifier and parameters of a discrete parametric distribution.
 *
 *  \param[in] iident       identifier,
 *  \param[in] iinf_bound   lower bound,
 *  \param[in] isup_bound   upper bound (binomial, uniform),
 *  \param[in] iparameter   parameter (Poisson, negative binomial, Poisson geometric),
 *  \param[in] iprobability probability (binomial, negative binomial, Poisson geometric).
 */
/*--------------------------------------------------------------*/

void DiscreteParametric::init(discrete_parametric iident , int iinf_bound , int isup_bound ,
                              double iparameter , double iprobability)

{
  ident = iident;

  inf_bound = iinf_bound;
  sup_bound = isup_bound;
  parameter = iparameter;
  probability = iprobability;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Constructor of the DiscreteParametric class.
 *
 *  \param[in] inb_value    number of values,
 *  \param[in] iident       identifier,
 *  \param[in] iinf_bound   lower bound,
 *  \param[in] isup_bound   upper bound (binomial, uniform),
 *  \param[in] iparameter   parameter (Poisson, negative binomial, Poisson geometric),
 *  \param[in] iprobability probability (binomial, negative binomial, Poisson geometric).
 */
/*--------------------------------------------------------------*/

DiscreteParametric::DiscreteParametric(int inb_value , discrete_parametric iident ,
                                       int iinf_bound , int isup_bound ,
                                       double iparameter , double iprobability)
:Distribution(inb_value)

{
  ident = iident;

  inf_bound = iinf_bound;
  sup_bound = isup_bound;
  parameter = iparameter;
  probability = iprobability;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Constructor of the DiscreteParametric class.
 *
 *  \param[in] iident          identifier,
 *  \param[in] iinf_bound      lower bound,
 *  \param[in] isup_bound      upper bound (binomial, uniform),
 *  \param[in] iparameter      parameter (Poisson, negative binomial, Poisson geometric),
 *  \param[in] iprobability    probability (binomial, negative binomial, Poisson geometric),
 *  \param[in] cumul_threshold threshold on the cumulative distribution function.
 */
/*--------------------------------------------------------------*/

DiscreteParametric::DiscreteParametric(discrete_parametric iident , int iinf_bound ,
                                       int isup_bound , double iparameter ,
                                       double iprobability , double cumul_threshold)

{
  ident = iident;

  inf_bound = iinf_bound;
  sup_bound = isup_bound;
  parameter = iparameter;
  probability = iprobability;

  nb_value = 0;

  if ((ident == BINOMIAL) || (ident == UNIFORM)) {
    nb_value = sup_bound + 1;
  }

  else if (ident == PRIOR_SEGMENT_LENGTH) {
    nb_value = sequence_length - no_segment + 2;
  }

  else if ((ident == POISSON) || (ident == NEGATIVE_BINOMIAL) || (ident == POISSON_GEOMETRIC)) {
    nb_value = (int)round(inf_bound + (parametric_mean_computation() - inf_bound +
                                       sqrt(parametric_variance_computation())) * 20.);
    if (nb_value == inf_bound) {
      nb_value++;
    }
  }

  Distribution::init(nb_value);

  computation(1 , cumul_threshold);
}


/*--------------------------------------------------------------*/
/**
 *  \brief Construction of a DiscreteParametric object from a Distribution object.
 *
 *  \param[in] dist            reference on a Distribution object,
 *  \param[in] ialloc_nb_value number of allocated values.
 */
/*--------------------------------------------------------------*/

DiscreteParametric::DiscreteParametric(const Distribution &dist , int ialloc_nb_value)
:Distribution(dist , DISTRIBUTION_COPY , ialloc_nb_value)

{
  ident = CATEGORICAL;

  inf_bound = I_DEFAULT;
  sup_bound = I_DEFAULT;
  parameter = D_DEFAULT;
  probability = D_DEFAULT;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Construction of a DiscreteParametric object from a Distribution object
 *         applying a scaling operation.
 *
 *  \param[in] dist          reference on a Distribution object, 
 *  \param[in] scaling_coeff scaling factor.
 */
/*--------------------------------------------------------------*/

DiscreteParametric::DiscreteParametric(const Distribution &dist , double scaling_coeff)
:Distribution(dist , scaling_coeff)

{
  ident = CATEGORICAL;

  inf_bound = I_DEFAULT;
  sup_bound = I_DEFAULT;
  parameter = D_DEFAULT;
  probability = D_DEFAULT;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Construction of an inter-event distribution from
 *         an initial inter-event distribution up- or down-scaling the time scale.
 *
 *  \param[in] dist          reference on an inter-event distribution,
 *  \param[in] scaling_coeff scaling factor.
 */
/*--------------------------------------------------------------*/

DiscreteParametric::DiscreteParametric(const DiscreteParametric &dist , double scaling_coeff)
:Distribution((int)floor(dist.nb_value * scaling_coeff) + 1)

{
  double scaled_mean , scaled_variance , ratio , shifted_mean;


  // scaled_mean = scaling_coeff * dist.mean;
  // scaled_variance = scaling_coeff * scaling_coeff * dist.variance;
  scaled_mean = scaling_coeff * dist.parametric_mean_computation();
  scaled_variance = scaling_coeff * scaling_coeff * dist.parametric_variance_computation();

  inf_bound = (int)floor(dist.inf_bound * scaling_coeff);
  sup_bound = I_DEFAULT;
  parameter = D_DEFAULT;
  probability = D_DEFAULT;

  shifted_mean = scaled_mean - inf_bound;
  ratio = scaled_variance / shifted_mean;

  // case binomial

  if (ratio < 1. - POISSON_RANGE) {
    ident = BINOMIAL;
    sup_bound = (int)ceil(inf_bound + shifted_mean * shifted_mean /
                (shifted_mean - scaled_variance));
    if (sup_bound <= inf_bound) {
      sup_bound = inf_bound + 1;
    }
    probability = shifted_mean / (sup_bound - inf_bound);
  }

  // case negative binomial

  else if (ratio > 1. + POISSON_RANGE) {
    ident = NEGATIVE_BINOMIAL;
    parameter = shifted_mean * shifted_mean / (scaled_variance - shifted_mean);
    probability = shifted_mean / scaled_variance;
  }

  // case Poisson

  else {
    ident = POISSON;
    parameter = shifted_mean;
  }

  computation();
}


/*--------------------------------------------------------------*/
/**
 *  \brief Construction of a DiscreteParametric object from
 *         a FrequencyDistribution object.
 *
 *  \param[in] histo reference on a FrequencyDistribution object.
 */
/*--------------------------------------------------------------*/

DiscreteParametric::DiscreteParametric(const FrequencyDistribution &histo)
:Distribution(histo)

{
  ident = CATEGORICAL;

  inf_bound = I_DEFAULT;
  sup_bound = I_DEFAULT;
  parameter = D_DEFAULT;
  probability = D_DEFAULT;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Copy of a DiscreteParametric object.
 *
 *  \param[in] dist reference on a DiscreteParametric object.
 */
/*--------------------------------------------------------------*/

void DiscreteParametric::copy(const DiscreteParametric &dist)

{
  ident = dist.ident;

  inf_bound = dist.inf_bound;
  sup_bound = dist.sup_bound;
  parameter = dist.parameter;
  probability = dist.probability;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Constructor by copy of the DiscreteParametric class.
 *
 *  \param[in] dist            reference on a DiscreteParametric object,
 *  \param[in] transform       type of transform (DISTRIBUTION_COPY/NORMALIZATION),
 *  \param[in] ialloc_nb_value number of allocated values.
 */
/*--------------------------------------------------------------*/

DiscreteParametric::DiscreteParametric(const DiscreteParametric &dist ,
                                       distribution_transformation transform , int ialloc_nb_value)
:Distribution(dist , transform , ialloc_nb_value)

{
  copy(dist);
}


/*--------------------------------------------------------------*/
/**
 *  \brief Assignment operator of the DiscreteParametric class.
 *
 *  \param[in] dist reference on a DiscreteParametric object.
 *
 *  \return         DiscreteParametric object.
 */
/*--------------------------------------------------------------*/

DiscreteParametric& DiscreteParametric::operator=(const DiscreteParametric &dist)

{
  if (&dist != this) {
    delete [] mass;
    delete [] cumul;

    Distribution::copy(dist);
    copy(dist);
  }

  return *this;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Analysis of the format of a DiscreteParametric object.
 *
 *  \param[in] error           reference on a StatError object,
 *  \param[in] in_file         stream,
 *  \param[in] line            reference on the file line index,
 *  \param[in] last_ident      identifier of the last distribution in the list,
 *  \param[in] cumul_threshold threshold on the cumulative distribution function,
 *  \param[in] min_inf_bound   minimum lower bound.
 *
 *  \return                    DiscreteParametric object.   
 */
/*--------------------------------------------------------------*/

DiscreteParametric* DiscreteParametric::parsing(StatError &error , ifstream &in_file ,
                                                int &line , discrete_parametric last_ident ,
                                                double cumul_threshold , int min_inf_bound)

{
  string buffer;
  size_t position;
  typedef tokenizer<char_separator<char>> tokenizer;
  char_separator<char> separator(" \t");
  bool status = true , lstatus;
  int i , j;
  discrete_parametric ident = CATEGORICAL;
  union {
    int inf_bound;
    int no_segment;
  };
  union {
    int sup_bound = I_DEFAULT;
    int sequence_length;
  };
  double parameter = D_DEFAULT , probability = D_DEFAULT;
  DiscreteParametric *dist;


  dist = NULL;

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

//    for (auto &token : tok_buffer)  {   in C++ 11, entails replacing pointer by reference for token
    for (tokenizer::iterator token = tok_buffer.begin();token != tok_buffer.end();token++) {

      // test distribution name

      if (i == 0) {
        for (j = BINOMIAL;j <= last_ident;j++) {
          if (*token == STAT_discrete_distribution_word[j]) {
            ident = (discrete_parametric)j;
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

          // 1st parameter: lower bound or number of segments (prior segment length)

          case 0 : {
            if (((ident == BINOMIAL) || (ident == POISSON) || (ident == NEGATIVE_BINOMIAL) ||
                 (ident == POISSON_GEOMETRIC) || (ident == UNIFORM)) && (*token != STAT_word[STATW_INF_BOUND])) {
              status = false;
              error.correction_update(STAT_parsing[STATP_PARAMETER_NAME] , STAT_word[STATW_INF_BOUND] , line , i + 1);
            }

            if ((ident == PRIOR_SEGMENT_LENGTH) && (*token != STAT_word[STATW_NO_SEGMENT])) {
              status = false;
              error.correction_update(STAT_parsing[STATP_PARAMETER_NAME] , STAT_word[STATW_NO_SEGMENT] , line , i + 1);
            }
            break;
          }

          // 2nd parameter: upper bound (binomial, uniform), parameter (Poisson, negative binomial, Poisson geometric)
          // or sequence length (prior segment length)

          case 1 : {
            if (((ident == BINOMIAL) || (ident == UNIFORM)) && (*token != STAT_word[STATW_SUP_BOUND])) {
              status = false;
              error.correction_update(STAT_parsing[STATP_PARAMETER_NAME] , STAT_word[STATW_SUP_BOUND] , line , i + 1);
            }

            if (((ident == POISSON) || (ident == NEGATIVE_BINOMIAL) || (ident == POISSON_GEOMETRIC)) &&
                (*token != STAT_word[STATW_PARAMETER])) {
              status = false;
              error.correction_update(STAT_parsing[STATP_PARAMETER_NAME] , STAT_word[STATW_PARAMETER] , line , i + 1);
            }

            if ((ident == PRIOR_SEGMENT_LENGTH) && (*token != STAT_word[STATW_SEQUENCE_LENGTH])) {
              status = false;
              error.correction_update(STAT_parsing[STATP_PARAMETER_NAME] , STAT_word[STATW_SEQUENCE_LENGTH] , line , i + 1);
            }
            break;
          }

          // 3rd parameter: probability (binomial, negative binomial, Poisson geometric)

          case 2 : {
            if (((ident == BINOMIAL) || (ident == NEGATIVE_BINOMIAL) || (ident == POISSON_GEOMETRIC)) &&
                (*token != STAT_word[STATW_PROBABILITY])) {
              status = false;
              error.correction_update(STAT_parsing[STATP_PARAMETER_NAME] , STAT_word[STATW_PROBABILITY] , line , i + 1);
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

          // 1st parameter: lower bound or number of segments (prior segment length)

          case 0 : {
            if ((ident == BINOMIAL) || (ident == POISSON) || (ident == NEGATIVE_BINOMIAL) ||
                (ident == POISSON_GEOMETRIC) || (ident == UNIFORM)) {
/*              try {
                inf_bound = stoi(*token);   in C++ 11
              }
              catch (invalid_argument &arg) {
                lstatus = false;
              } */
              inf_bound = atoi(token->c_str());

              if ((lstatus) && ((inf_bound < min_inf_bound) || (inf_bound > MAX_INF_BOUND))) {
                lstatus = false;
              }
            }

            if (ident == PRIOR_SEGMENT_LENGTH) {
/*              try {
                no_segment = stoi(*token);   in C++ 11
              }
              catch(invalid_argument &arg) {
                lstatus = false;
              } */
              no_segment = atoi(token->c_str());

              if ((lstatus) && (no_segment < 2)) {
                lstatus = false;
              }
            }
            break;
          }

          // 2nd parameter: upper bound (binomial, uniform), parameter (Poisson, negative binomial, Poisson geometric)
          // or sequence length (prior segment length)

          case 1 : {
            if ((ident == BINOMIAL) || (ident == UNIFORM)) {
/*              try {
                sup_bound = stoi(*token);   in C++ 11
              }
              catch(invalid_argument &arg) {
                lstatus = false;
              } */
              sup_bound = atoi(token->c_str());

              if ((lstatus) && (status)) {
                switch (ident) {

                case BINOMIAL : {
                  if ((sup_bound <= inf_bound) || (sup_bound - inf_bound > MAX_DIFF_BOUND)) {
                    lstatus = false;
                  }
                  break;
                }

                case UNIFORM : {
                  if ((sup_bound < inf_bound) || (sup_bound - inf_bound > MAX_DIFF_BOUND)) {
                    lstatus = false;
                  }
                  break;
                }
                }
              }
            }

            if ((ident == POISSON) || (ident == NEGATIVE_BINOMIAL) || (ident == POISSON_GEOMETRIC)) {
/*              try {
                parameter = stod(*token);   in C++ 11
              }
              catch(invalid_argument &arg) {
                lstatus = false;
              } */
              parameter = atof(token->c_str());

              if ((lstatus) && ((parameter <= 0.) || ((ident == POISSON) && (parameter > MAX_MEAN)))) {
                lstatus = false;
              }
            }

            if (ident == PRIOR_SEGMENT_LENGTH) {
/*              try {
                sequence_length = stoi(*token);   in C++ 11
              }
              catch(invalid_argument &arg) {
                lstatus = false;
              } */
              sequence_length = atoi(token->c_str());

              if ((lstatus) && (status) && (sequence_length <= no_segment)) {
                lstatus = false;
              }
            }
            break;
          }

          // 3rd parameter: probability (binomial, negative binomial, Poisson geometric)

          case 2 : {
            if ((ident == BINOMIAL) || (ident == NEGATIVE_BINOMIAL) || (ident == POISSON_GEOMETRIC)) {
/*              try {
                probability = stod(*token);   in C++ 11
              }
              catch(invalid_argument &arg) {
                lstatus = false;
              } */
              probability = atof(token->c_str());

              if (lstatus) {
                switch (ident) {

                case BINOMIAL : {
                  if ((probability < 0.) || (probability > 1.)) {
                    lstatus = false;
                  }
                  break;
                }

                case NEGATIVE_BINOMIAL : {
                  if ((probability <= 0.) || (probability >= 1.)) {
                    lstatus = false;
                  }
                  else if ((status) && (parameter * (1. - probability) / probability > MAX_MEAN)) {
                    lstatus = false;
                  }
                  break;
                }

                case POISSON_GEOMETRIC : {
                  if ((probability <= 0.) || (probability >= 1.)) {
                    lstatus = false;
                  }
                  else if ((status) && ((inf_bound + parameter) / probability > MAX_MEAN)) {
                    lstatus = false;
                  }
                  break;
                }
                }
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
      if ((((ident == BINOMIAL) || (ident == NEGATIVE_BINOMIAL) || (ident == POISSON_GEOMETRIC)) && (i != 10)) ||
          (((ident == POISSON) || (ident == UNIFORM) || (ident == PRIOR_SEGMENT_LENGTH)) && (i != 7))) {
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
    dist = new DiscreteParametric(ident , inf_bound , sup_bound ,
                                  parameter , probability , cumul_threshold);
  }

  return dist;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of the parameters of a discrete distribution.
 *
 *  \param[in,out] os stream.
 */
/*--------------------------------------------------------------*/

ostream& DiscreteParametric::ascii_print(ostream &os) const

{
  os << STAT_discrete_distribution_word[ident];

  if (ident != PRIOR_SEGMENT_LENGTH) {
    if (inf_bound != I_DEFAULT) {
      os << "   " << STAT_word[STATW_INF_BOUND] << " : " << inf_bound;
    }
    if (sup_bound != I_DEFAULT) {
      os << "   " << STAT_word[STATW_SUP_BOUND] << " : " << sup_bound;
    }
    if (parameter != D_DEFAULT) {
      os << "   " << STAT_word[STATW_PARAMETER] << " : " << parameter;
    }
    if (probability != D_DEFAULT) {
      os << "   " << STAT_word[STATW_PROBABILITY] << " : " << probability;
    }
  }

  else {
    os << "   " << STAT_word[STATW_NO_SEGMENT] << " : " << no_segment
       << "   " << STAT_word[STATW_SEQUENCE_LENGTH] << " : " << sequence_length;
  }
  os << endl;

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of the parameters of a discrete distribution at the spreadsheet format.
 *
 *  \param[in,out] os stream.
 */
/*--------------------------------------------------------------*/

ostream& DiscreteParametric::spreadsheet_print(ostream &os) const

{
  os << STAT_discrete_distribution_word[ident];

  if (ident != PRIOR_SEGMENT_LENGTH) {
    if (inf_bound != I_DEFAULT) {
      os << "\t" << STAT_word[STATW_INF_BOUND] << "\t" << inf_bound;
    }
    if (sup_bound != I_DEFAULT) {
      os << "\t" << STAT_word[STATW_SUP_BOUND] << "\t" << sup_bound;
    }
    if (parameter != D_DEFAULT) {
      os << "\t" << STAT_word[STATW_PARAMETER] << "\t" << parameter;
    }
    if (probability != D_DEFAULT) {
      os << "\t" << STAT_word[STATW_PROBABILITY] << "\t" << probability;
    }
  }

  else {
    os << "\t" << STAT_word[STATW_NO_SEGMENT] << "\t" << no_segment
       << "\t" << STAT_word[STATW_SEQUENCE_LENGTH] << "\t" << sequence_length;
  }
  os << endl;

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of the characteristics of a discrete parametric distribution.
 *
 *  \param[in,out] os           stream,
 *  \param[in]     shape        flag on the writing of the shape characteristics,
 *  \param[in]     comment_flag flag comment.
 */
/*--------------------------------------------------------------*/

ostream& DiscreteParametric::ascii_parametric_characteristic_print(ostream &os , bool shape , bool comment_flag) const

{
  if (ident == CATEGORICAL) {
    ascii_characteristic_print(os , shape , comment_flag);
  }

  else {
    double variance = parametric_variance_computation();


    if (comment_flag) {
      os << "# ";
    }
    os << STAT_label[STATL_MEAN] << ": " << parametric_mean_computation() << "   "
       << STAT_label[STATL_MEDIAN] << ": " << quantile_computation() << "   "
       << STAT_label[STATL_MODE] << ": " << mode_computation() << endl;

    if (comment_flag) {
      os << "# ";
    }
    os << STAT_label[STATL_VARIANCE] << ": " << variance << "   "
       << STAT_label[STATL_STANDARD_DEVIATION] << ": " << sqrt(variance);
    if (variance > 0.) {
      os << "   " << STAT_label[STATL_LOWER_QUARTILE] << ": " << quantile_computation(0.25)
         << "   " << STAT_label[STATL_UPPER_QUARTILE] << ": " << quantile_computation(0.75);
    }
    os << endl;

    if ((shape) && (variance > 0.)) {
      if (comment_flag) {
        os << "# ";
      }
      os << STAT_label[STATL_SKEWNESS_COEFF] << ": " << parametric_skewness_computation() << "   "
         << STAT_label[STATL_KURTOSIS_COEFF] << ": " << parametric_kurtosis_computation() << endl;
    }
  }

# ifdef MESSAGE
  if (ident == PRIOR_SEGMENT_LENGTH) {
    ascii_characteristic_print(os , shape , comment_flag);
  }
# endif

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of the characteristics of a discrete parametric distribution at the spreadsheet format.
 *
 *  \param[in,out] os    stream,
 *  \param[in]     shape flag on the writing of the shape characteristics.
 */
/*--------------------------------------------------------------*/

ostream& DiscreteParametric::spreadsheet_parametric_characteristic_print(ostream &os , bool shape) const

{
  if (ident == CATEGORICAL) {
    spreadsheet_characteristic_print(os , shape);
  }

  else {
    double variance = parametric_variance_computation();


    os << STAT_label[STATL_MEAN] << "\t" << parametric_mean_computation() << "\t\t"
       << STAT_label[STATL_MEDIAN] << "\t" << quantile_computation() << "\t\t"
       << STAT_label[STATL_MODE] << "\t" << mode_computation() << endl;

    os << STAT_label[STATL_VARIANCE] << "\t" << variance << "\t\t"
       << STAT_label[STATL_STANDARD_DEVIATION] << "\t" << sqrt(variance);
    if (variance > 0.) {
      os << "\t\t" << STAT_label[STATL_LOWER_QUARTILE] << "\t" << quantile_computation(0.25)
         << "\t\t" << STAT_label[STATL_UPPER_QUARTILE] << "\t" << quantile_computation(0.75);
    }
    os << endl;

    if ((shape) && (variance > 0.)) {
      os << STAT_label[STATL_SKEWNESS_COEFF] << "\t" << parametric_skewness_computation() << "\t\t"
         << STAT_label[STATL_KURTOSIS_COEFF] << "\t" << parametric_kurtosis_computation() << endl;
    }
  }

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of the parameters of a discrete parametric distribution at the Gnuplot format.
 *
 *  \param[in,out] os stream.
 */
/*--------------------------------------------------------------*/

ostream& DiscreteParametric::plot_title_print(ostream &os) const

{
  if (ident != CATEGORICAL) {
    os << " " << STAT_discrete_distribution_letter[ident] << "(";

    if (ident != PRIOR_SEGMENT_LENGTH) {
      os << inf_bound;
      if (sup_bound != I_DEFAULT) {
        os << ", " << sup_bound;
      }
      if (parameter != D_DEFAULT) {
        os << ", " << parameter;
      }
      if (probability != D_DEFAULT) {
        os << ", " << probability;
      }
    }

    else {
      os << no_segment << ", " << sequence_length;
    }
    os << ")";
  }

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Display of a discrete parametric distribution.
 *
 *  \param[in,out] os   stream,
 *  \param[in]     dist reference on a DiscreteParametric object.
 */
/*--------------------------------------------------------------*/

ostream& operator<<(ostream &os , const DiscreteParametric &dist)

{
  streamsize nb_digits;


  nb_digits = os.precision(5);

  os << endl;
  dist.ascii_print(os);
  dist.Distribution::print(os);

  os.precision(nb_digits);

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the number of parameters of a discrete parametric distribution.
 *
 *  \return number of parameters.
 */
/*--------------------------------------------------------------*/

int DiscreteParametric::nb_parameter_computation()

{
  int bnb_parameter;


  switch (ident) {
  case BINOMIAL :
    bnb_parameter = 3;
    break;
  case POISSON :
    bnb_parameter = 2;
    break;
  case NEGATIVE_BINOMIAL :
    bnb_parameter = 3;
    break;
  case POISSON_GEOMETRIC :
    bnb_parameter = 3;
    break;
  case UNIFORM :
    bnb_parameter = 2;
    break;
  case PRIOR_SEGMENT_LENGTH :
    bnb_parameter = 2;
    break;
  default :
    bnb_parameter = 0;
    break;
  }

  return bnb_parameter;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Update of the number of parameters of a discrete parametric distribution.
 */
/*--------------------------------------------------------------*/

void DiscreteParametric::nb_parameter_update()

{
  nb_parameter = nb_parameter_computation();
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the theooretical mean of a discrete parametric distribution.
 *
 *  \return theoretical mean
 */
/*--------------------------------------------------------------*/

double DiscreteParametric::parametric_mean_computation() const

{
  double parametric_mean;


  switch (ident) {
  case BINOMIAL :
    parametric_mean = inf_bound + (sup_bound - inf_bound) * probability;
    break;
  case POISSON :
    parametric_mean = inf_bound + parameter;
    break;
  case NEGATIVE_BINOMIAL :
    parametric_mean = inf_bound + parameter * (1. - probability) / probability;
    break;
  case POISSON_GEOMETRIC :
    parametric_mean = (inf_bound + parameter) / probability;
    break;
  case UNIFORM :
    parametric_mean = (double)(inf_bound + sup_bound) / 2.;
    break;
  case PRIOR_SEGMENT_LENGTH :
    parametric_mean = (double)sequence_length / (double)no_segment;
    break;
  default :
    parametric_mean = mean;
    break;
  }

  return parametric_mean;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the theoretical variance of a discrete parametric distribution.
 *
 *  \return theoretical variance.
 */
/*--------------------------------------------------------------*/

double DiscreteParametric::parametric_variance_computation() const

{
  double parametric_variance;


  switch (ident) {
  case BINOMIAL :
    parametric_variance = (sup_bound - inf_bound) * probability * (1. - probability);
    break;
  case POISSON :
    parametric_variance = parameter;
    break;
  case NEGATIVE_BINOMIAL :
    parametric_variance = parameter * (1. - probability) / (probability * probability);
    break;
  case POISSON_GEOMETRIC :
    parametric_variance = ((inf_bound + parameter) * (1. - probability) + parameter) /
                          (probability * probability);
    break;
  case UNIFORM :
    parametric_variance = (double)((sup_bound - inf_bound + 2) *
                          (sup_bound - inf_bound)) / 12.;
    break;
  case PRIOR_SEGMENT_LENGTH :
    parametric_variance = ((double)sequence_length * (sequence_length - no_segment) * (no_segment - 1)) /
                          ((double)no_segment * no_segment * (no_segment + 1));
    break;
  default :
    parametric_variance = variance;
    break;
  }

  return parametric_variance;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the theoretical coefficient of skewness of a discrete parametric distribution.
 *
 *  \return theoretical coefficient of skewness.
 */
/*--------------------------------------------------------------*/

double DiscreteParametric::parametric_skewness_computation() const

{
  double parametric_skewness;


  switch (ident) {

  case BINOMIAL : {
    if ((probability > 0.) && (probability < 1.)) {
      parametric_skewness = (1. - 2. * probability) /
                            sqrt((sup_bound - inf_bound) * probability * (1. - probability));
    }
    else {
      parametric_skewness = 0.;
    }
    break;
  }

  case POISSON : {
    parametric_skewness = 1. / sqrt(parameter);
    break;
  }

  case NEGATIVE_BINOMIAL : {
    parametric_skewness = (2. - probability) / sqrt(parameter * (1. - probability));
    break;
  }

  case POISSON_GEOMETRIC : {
    parametric_skewness = (parameter * (1. + 3 * (1. - probability)) +
                           (inf_bound + parameter) * (2. - probability) * (1. - probability)) /
                          pow((inf_bound + parameter) * (1. - probability) + parameter , 1.5);
    break;
  }

  case UNIFORM : {
    parametric_skewness = 0.;
    break;
  }

  case PRIOR_SEGMENT_LENGTH : {
    parametric_skewness = (((double)(no_segment - 2) * (2 * sequence_length - no_segment)) / (double)(no_segment + 2)) *
                          sqrt((double)(no_segment + 1) / ((double)sequence_length * (sequence_length - no_segment) * (no_segment - 1)));
    break;
  }

  default : {
    parametric_skewness = skewness_computation();
    break;
  }
  }

  return parametric_skewness;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the theoretical excess kurtosis of a discrete parametric distribution:
 *         excess kurtosis = coefficient of kurtosis - 3.
 *
 *  \return theoretical excess kurtosis
 */
/*--------------------------------------------------------------*/

double DiscreteParametric::parametric_kurtosis_computation() const

{
  double parametric_kurtosis;


  switch (ident) {

  case BINOMIAL : {
    if ((probability > 0.) && (probability < 1.)) {
      parametric_kurtosis = (1. - 6. * probability * (1. - probability)) /
                            ((sup_bound - inf_bound) * probability * (1. - probability));
    }
    else {
      parametric_kurtosis = -D_INF;
    }
    break;
  }

  case POISSON : {
    parametric_kurtosis = 1. / parameter;
    break;
  }

  case NEGATIVE_BINOMIAL : {
    parametric_kurtosis = (probability * probability + 6. * (1. - probability)) /
                          (parameter * (1. - probability));
    break;
  }

  case UNIFORM : {
    parametric_kurtosis = -0.5;
    break;
  }

  case PRIOR_SEGMENT_LENGTH : {
    parametric_kurtosis = ((double)(no_segment + 1) * (-no_segment * no_segment * no_segment * (no_segment - 1) * (no_segment - 6) +
                            2 * no_segment * no_segment * sequence_length * (5 * no_segment * no_segment - 14 * no_segment + 12) +
                            3 * sequence_length * sequence_length * (sequence_length - 2 * no_segment) *
                            (3 * no_segment * no_segment - 7 * no_segment + 6))) /
                          ((double)sequence_length * (sequence_length - no_segment) * (sequence_length - no_segment) * (no_segment - 1) *
                           (no_segment + 2) * (no_segment + 3)) - 3.;;
    break;
  }

  default : {
    parametric_kurtosis = kurtosis_computation();
    break;
  }
  }

  return parametric_kurtosis;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the distance between 2 discrete parametric distributions
 *         (sup of the absolute difference between the cumulative distribution functions
 *          in the case of non-crossing cumulative distribution functions; in the general case,
 *          sum of sup on intervals between 2 crossings of cumulative distribution functions).
 *
 *  \param[in] dist reference on a DiscreteParametric object.
 *
 *  \return         distance between 2 discrete parametric distributions
 */
/*--------------------------------------------------------------*/

double DiscreteParametric::sup_norm_distance_computation(const DiscreteParametric &dist) const

{
  bool crossing;
  int i;
  int inf , sup;
  double buff , distance , max_absolute_diff;


  inf = MAX(offset , dist.offset);
  sup = MIN(nb_value , dist.nb_value) - 1;
  if (sup < inf) {
    distance = 1.;
  }

  else {
    distance = 0.;
    max_absolute_diff = fabs(cumul[inf] - dist.cumul[inf]);
    crossing = false;

    for (i = inf + 1;i <= sup;i++) {
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
    distance += max_absolute_diff;

#   ifdef DEBUG
    double overlap;

    overlap = 0.;
    for (i = inf;i <= sup;i++) {
      overlap += MIN(mass[i] , dist.mass[i]);
    }

    cout << "\n" << STAT_label[STATL_SUP_NORM_DISTANCE] << ": " << distance << " " << 1. - overlap << endl;
#   endif
  }

  return distance;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Construction of a DiscreteParametricModel object from
 *         a FrequencyDistribution object.
 *
 *  \param[in] histo reference on a FrequencyDistribution object.
 */
/*--------------------------------------------------------------*/

DiscreteParametricModel::DiscreteParametricModel(const FrequencyDistribution &histo)
:DiscreteParametric(histo)

{
  frequency_distribution = new DiscreteDistributionData(histo);
}


/*--------------------------------------------------------------*/
/**
 *  \brief Construction of a DiscreteParametricModel object from
 *         a Distribution object and a FrequencyDistribution object.
 *
 *  \param[in] dist  reference on a Distribution object,
 *  \param[in] histo pointer on a FrequencyDistribution object.
 */
/*--------------------------------------------------------------*/

DiscreteParametricModel::DiscreteParametricModel(const Distribution &dist ,
                                                 const FrequencyDistribution *histo)
:DiscreteParametric(dist)

{
  if ((histo) && (histo->nb_element > 0)) {
    frequency_distribution = new DiscreteDistributionData(*histo);
  }
  else {
    frequency_distribution = NULL;
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Construction of a DiscreteParametricModel object from
 *         a DiscreteParametric object and a FrequencyDistribution object.
 *
 *  \param[in] dist  reference on a DiscreteParametric object,
 *  \param[in] histo pointer on a FrequencyDistribution object.
 */
/*--------------------------------------------------------------*/

DiscreteParametricModel::DiscreteParametricModel(const DiscreteParametric &dist ,
                                                 const FrequencyDistribution *histo)
:DiscreteParametric(dist)

{
  if ((histo) && (histo->nb_element > 0)) {
    frequency_distribution = new DiscreteDistributionData(*histo);
  }
  else {
    frequency_distribution = NULL;
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Constructor by copy of the DiscreteParametricModel class.
 *
 *  \param[in] dist      reference on a DiscreteParametricModel object,
 *  \param[in] data_flag flag copy of the DiscreteDistributionData object.
 */
/*--------------------------------------------------------------*/

DiscreteParametricModel::DiscreteParametricModel(const DiscreteParametricModel &dist ,
                                                 bool data_flag)
:DiscreteParametric(dist)

{
  if ((data_flag) && (dist.frequency_distribution)) {
    frequency_distribution = new DiscreteDistributionData(*(dist.frequency_distribution) , false);
  }
  else {
    frequency_distribution = NULL;
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Destructor of the DiscreteParametricModel class.
 */
/*--------------------------------------------------------------*/

DiscreteParametricModel::~DiscreteParametricModel()

{
  if (frequency_distribution) {
    delete frequency_distribution; 
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Assignment operator of the DiscreteParametricModel class.
 *
 *  \param[in] dist reference on a DiscreteParametricModel object.
 *
 *  \return         DiscreteParametricModel object.
 */
/*--------------------------------------------------------------*/

DiscreteParametricModel& DiscreteParametricModel::operator=(const DiscreteParametricModel &dist)

{
  if (&dist != this) {
    delete frequency_distribution;

    delete [] mass;
    delete [] cumul;

    Distribution::copy(dist);
    DiscreteParametric::copy(dist);

    if (dist.frequency_distribution) {
      frequency_distribution = new DiscreteDistributionData(*(dist.frequency_distribution) , false);
    }
    else {
      frequency_distribution = NULL;
    }
  }

  return *this;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Extraction of the DiscreteDistributionData object included in
 *         a DiscreteParametricModel object.
 *
 *  \param[in] error reference on a StatError object.
 *
 *  \return          DiscreteDistributionData object.
 */
/*--------------------------------------------------------------*/

DiscreteDistributionData* DiscreteParametricModel::extract_data(StatError &error) const

{
  DiscreteDistributionData *histo;


  error.init();

  if (!frequency_distribution) {
    histo = NULL;
    error.update(STAT_error[STATR_NON_EXISTING_FREQUENCY_DISTRIBUTION]);
  }

  else {
    histo = new DiscreteDistributionData(*frequency_distribution);
    histo->distribution = new DiscreteParametricModel(*this , false);
  }

  return histo;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Construction of a DiscreteParametricModel object from a file.
 *
 *  \param[in] error           reference on a StatError object,
 *  \param[in] path            file path,
 *  \param[in] cumul_threshold threshold on the cumulative distribution function.
 *
 *  \return                    DiscreteParametricModel object.
 */
/*--------------------------------------------------------------*/

DiscreteParametricModel* DiscreteParametricModel::ascii_read(StatError &error , const string path ,
                                                             double cumul_threshold)

{
  string buffer;
  size_t position;
  bool status;
  int line;
  DiscreteParametric *pdist;
  DiscreteParametricModel *dist;
  ifstream in_file(path.c_str());


  dist = NULL;
  error.init();

  if (!in_file) {
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {
    status = true;
    line = 0;

    pdist = DiscreteParametric::parsing(error , in_file , line ,
                                        UNIFORM , cumul_threshold);

    if (!pdist) {
      status = false;
    }

    while (getline(in_file , buffer)) {
      line++;

#     ifdef DEBUG
      cout << line << " " << buffer << endl;
#     endif

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
      dist = new DiscreteParametricModel(*pdist);
    }

    delete pdist;
  }

  return dist;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing on a single line of a DiscreteParametricModel object.
 *
 *  \param[in,out] os stream.
 */
/*--------------------------------------------------------------*/

ostream& DiscreteParametricModel::line_write(ostream &os) const

{
  os << STAT_discrete_distribution_word[ident];

  if (ident == CATEGORICAL) {
    if ((mean != D_DEFAULT) && (variance != D_DEFAULT)) {
      os << "   " << STAT_label[STATL_MEAN] << ": " << mean
         << "   " << STAT_label[STATL_VARIANCE] << ": " << variance;
    }
  }

  else {
    os << "   " << STAT_label[STATL_MEAN] << ": " << parametric_mean_computation()
       << "   " << STAT_label[STATL_VARIANCE] << ": " << parametric_variance_computation();
  }

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a discrete parametric distribution and
 *         the associated frequency distribution.
 *
 *  \param[in,out] os         stream,
 *  \param[in]     histo      pointer on a DiscreteDistributionData object,
 *  \param[in]     exhaustive flag detail level,
 *  \param[in]     file_flag  flag file.
 */
/*--------------------------------------------------------------*/

ostream& DiscreteParametricModel::ascii_write(ostream &os , const DiscreteDistributionData *histo ,
                                              bool exhaustive , bool file_flag) const

{
  if (ident == CATEGORICAL) {
    file_flag = false;
  }

  ascii_print(os);
  if (complement > 0.) {
    os << STAT_label[STATL_UNPROPER] << " " << STAT_label[STATL_DISTRIBUTION] << " ("
       << STAT_label[STATL_COMPLEMENTARY_PROBABILITY] << ": " << complement << ")" << endl;
  }

  ascii_parametric_characteristic_print(os , true , file_flag);

  if (file_flag) {
    os << "# ";
  }
  os << STAT_label[STATL_MEAN_ABSOLUTE_DEVIATION] << ": " << mean_absolute_deviation_computation(mean);
  if (mean > 0.) {
    os << "   " << STAT_label[STATL_CONCENTRATION_COEFF] << ": " << concentration_computation();
  }
  os << endl;

  if (histo) {
    double likelihood , information;
    Test test(CHI2);


    os << "\n";
    if (file_flag) {
      os << "# ";
    }
    os << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";
    histo->ascii_characteristic_print(os , true , file_flag);

    if (file_flag) {
      os << "# ";
    }
    os << STAT_label[STATL_MEAN_ABSOLUTE_DEVIATION] << ": " << histo->mean_absolute_deviation_computation(histo->mean);
    if (histo->mean > 0.) {
      os << "   " << STAT_label[STATL_CONCENTRATION_COEFF] << ": " << histo->concentration_computation();
    }
    os << endl;

    likelihood = likelihood_computation(*histo);
    information = histo->information_computation();

    os << "\n";
    if (file_flag) {
      os << "# ";
    }
    os << STAT_label[STATL_LIKELIHOOD] << ": " << likelihood << "   ("
       << STAT_label[STATL_NORMALIZED] << ": " << likelihood / histo->nb_element << ")" << endl;

    if (file_flag) {
      os << "# ";
    }
    os << STAT_label[STATL_MAX_LIKELIHOOD] << ": " << information << "   ("
       << STAT_label[STATL_INFORMATION] << ": " << information / histo->nb_element << ")" << endl;

    if (file_flag) {
      os << "# ";
    }
    os << STAT_label[STATL_DEVIANCE] << ": " << 2 * (information - likelihood) << endl;

    chi2_fit(*histo , test);
    os << "\n";
    test.ascii_print(os , file_flag);
  }

  else {
    if (file_flag) {
      os << "# ";
    }
    os << STAT_label[STATL_INFORMATION] << ": " << information_computation() << endl;
  }

  if (exhaustive) {
    os << "\n";
    if (file_flag) {
      os << "# ";
    }
    os << "  ";
    if (histo) {
      os << " | " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
    }
    os << " | " << STAT_label[STATL_DISTRIBUTION];
    if (histo) {
      os << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " "
         << STAT_label[STATL_FUNCTION];
    }
    os << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_DISTRIBUTION] << " "
       << STAT_label[STATL_FUNCTION] << endl;

    Distribution::ascii_print(os , file_flag , true , false , histo);
  }

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a DiscreteParametricModel object.
 *
 *  \param[in,out] os         stream,
 *  \param[in]     exhaustive flag detail level.
 */
/*--------------------------------------------------------------*/

ostream& DiscreteParametricModel::ascii_write(ostream &os , bool exhaustive) const

{
  return ascii_write(os , frequency_distribution , exhaustive , false);
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a DiscreteParametricModel object in a file.
 *
 *  \param[in] error      reference on a StatError object,
 *  \param[in] path       file path,
 *  \param[in] exhaustive flag detail level.
 *
 *  \return               error status.
 */
/*--------------------------------------------------------------*/

bool DiscreteParametricModel::ascii_write(StatError &error , const string path ,
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
    ascii_write(out_file , frequency_distribution , exhaustive , true);
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a discrete parametric distribution and the associated
 *         frequency distribution in a file at the spreadsheet format.
 *
 *  \param[in,out] os    stream,
 *  \param[in]     histo pointer on a DiscreteDistributionData object.
 */
/*--------------------------------------------------------------*/

ostream& DiscreteParametricModel::spreadsheet_write(ostream &os ,
                                                    const DiscreteDistributionData *histo) const

{


  spreadsheet_print(os);
  if (complement > 0.) {
    os << STAT_label[STATL_UNPROPER] << " " << STAT_label[STATL_DISTRIBUTION] << "\t"
       << STAT_label[STATL_COMPLEMENTARY_PROBABILITY] << "\t" << complement << endl;
  }

  spreadsheet_parametric_characteristic_print(os , true);

  os << STAT_label[STATL_MEAN_ABSOLUTE_DEVIATION] << "\t" << mean_absolute_deviation_computation(mean);
  if (mean > 0.) {
    os << "\t\t" << STAT_label[STATL_CONCENTRATION_COEFF] << "\t" << concentration_computation();
  }
  os << endl;

  if (histo) {
    double likelihood , information;
    Test test(CHI2);


    os << "\n" << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\t";
    histo->spreadsheet_characteristic_print(os , true);

    os << STAT_label[STATL_MEAN_ABSOLUTE_DEVIATION] << "\t" << histo->mean_absolute_deviation_computation(histo->mean);
    if (histo->mean > 0.) {
      os << "\t\t" << STAT_label[STATL_CONCENTRATION_COEFF] << "\t" << histo->concentration_computation();
    }
    os << endl;

    likelihood = likelihood_computation(*histo);
    information = histo->information_computation();

    os << "\n" << STAT_label[STATL_LIKELIHOOD] << "\t" << likelihood << "\t"
       << STAT_label[STATL_NORMALIZED] << "\t" << likelihood / histo->nb_element << endl;
    os << STAT_label[STATL_MAX_LIKELIHOOD] << "\t" << information << "\t"
       << STAT_label[STATL_INFORMATION] << "\t" << information / histo->nb_element << endl;
    os << STAT_label[STATL_DEVIANCE] << "\t" << 2 * (information - likelihood) << endl;

    chi2_fit(*histo , test);
    os << "\n";
    test.spreadsheet_print(os);
  }

  else {
    os << STAT_label[STATL_INFORMATION] << "\t" << information_computation() << endl;
  }

  os << "\n";
  if (histo) {
    os << "\t" << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
  }
  os << "\t" << STAT_label[STATL_DISTRIBUTION];
  if (histo) {
    os << "\t" << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " "
       << STAT_label[STATL_FUNCTION];
  }
  os << "\t" << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_DISTRIBUTION] << " "
     << STAT_label[STATL_FUNCTION];
  if (mean > 0.) {
    if ((histo) && (histo->mean > 0.)) {
      os << "\t" << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " " << STAT_label[STATL_CONCENTRATION] << " "
         << STAT_label[STATL_FUNCTION] << "\t" << STAT_label[STATL_DISTRIBUTION] << " ";
    }
    else {
      os << "\t";
    }
    os << STAT_label[STATL_CONCENTRATION] << " " << STAT_label[STATL_FUNCTION];
  }
  os << endl;

  Distribution::spreadsheet_print(os , true , true , false , histo);

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a DiscreteParametricModel object in a file at the spreadsheet format.
 *
 *  \param[in] error reference on a StatError object,
 *  \param[in] path  file path.
 *
 *  \return          error status.
 */
/*--------------------------------------------------------------*/

bool DiscreteParametricModel::spreadsheet_write(StatError &error , const string path) const

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
    spreadsheet_write(out_file , frequency_distribution);
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Plot of a discrete parametric distribution and
 *         the associated frequency distribution using Gnuplot.
 *
 *  \param[in] error  reference on a StatError object,
 *  \param[in] prefix file prefix,
 *  \param[in] title  figure title,
 *  \param[in] histo  pointer on a DiscreteDistributionData object.
 *
 *  \return           error status.
 */
/*--------------------------------------------------------------*/

bool DiscreteParametricModel::plot_write(StatError &error , const char *prefix , const char *title ,
                                         const DiscreteDistributionData *histo) const

{
  bool status;


  if (histo) {
    int i;
    int plot_nb_value , *poffset , *pnb_value;
    double scale , **pcumul , **pconcentration;
    ostringstream *data_file_name;


    error.init();

    // writing of the data files

    data_file_name = new ostringstream[3];

    poffset = new int[2];
    pnb_value = new int[2];
    pcumul = new double*[2];

    pconcentration = new double*[2];
    pconcentration[1] = NULL;

    data_file_name[0] << prefix << 0 << ".dat";

    poffset[0] = histo->offset;
    pnb_value[0] = histo->nb_value;

    // computation of the cumulative distribution function and concentration function of the frequency distribution

    scale = histo->nb_element / (1. - complement);
    pcumul[0] = histo->cumul_computation(scale);
    pconcentration[0] = histo->concentration_function_computation(scale);

    status = histo->plot_print((data_file_name[0].str()).c_str() , pcumul[0] , pconcentration[0]);

    if (status) {
      data_file_name[1] << prefix << 1 << ".dat";

      poffset[1] = offset;
      pnb_value[1] = nb_value;
      plot_nb_value = plot_nb_value_computation(histo);
      pcumul[1] = cumul;

      // computation of the concentration function

      pconcentration[1] = concentration_function_computation();
      plot_print((data_file_name[1].str()).c_str() , pconcentration[1] , scale);

      if ((variance > 0.) && (histo->variance > 0.)) {
        data_file_name[2] << prefix << 2 << ".dat";
        cumul_matching_plot_print((data_file_name[2].str()).c_str() , 2 , poffset , pnb_value , pcumul);
      }
    }

    delete [] poffset;
    delete [] pnb_value;

    delete [] pcumul[0];
    delete [] pcumul;

    delete [] pconcentration[0];
    delete [] pconcentration[1];
    delete [] pconcentration;

    if (status) {

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

        // distribution fit

        if (plot_nb_value - 1 < TIC_THRESHOLD) {
          out_file << "set xtics 0,1" << endl;
        }

        out_file << "plot [0:" << MAX(plot_nb_value - 1 , 1) << "] [0:"
                 << (int)(MAX(histo->max , max * scale) * YSCALE) + 1
                 << "] \"" << label((data_file_name[0].str()).c_str()) << "\" using 1:3 title \""
                 << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\" with impulses,\\" << endl;
        out_file << "\"" << label((data_file_name[1].str()).c_str()) << "\" using 1:2 title \""
                 << STAT_label[STATL_DISTRIBUTION];
        plot_title_print(out_file);
        out_file << "\" with linespoints" << endl;

        if (plot_nb_value - 1 < TIC_THRESHOLD) {
          out_file << "set xtics autofreq" << endl;
        }

        if (variance > 0.) {

          // cumulative distribution functions

          if (i == 0) {
            out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
          }
          out_file << endl;

          if (plot_nb_value - 1 < TIC_THRESHOLD) {
            out_file << "set xtics 0,1" << endl;
          }

          out_file << "plot [0:" << plot_nb_value - 1 << "] [0:" << 1. - complement << "] ";
          if (histo->variance > 0.) {
            out_file << "\"" << label((data_file_name[0].str()).c_str()) << "\" using 1:5 title \""
                     << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " "
                     << STAT_label[STATL_FUNCTION] << "\" with linespoints,\\" << endl;
          }
          out_file << "\"" << label((data_file_name[1].str()).c_str()) << "\" using 1:3 title \""
                   << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_DISTRIBUTION] << " "
                   << STAT_label[STATL_FUNCTION] << "\" with linespoints" << endl;

          if (plot_nb_value - 1 < TIC_THRESHOLD) {
            out_file << "set xtics autofreq" << endl;
          }

          if (histo->variance > 0.) {

            // matching of cumulative distribution functions taking as reference the distribution one

            if (i == 0) {
              out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
            }
            out_file << endl;

            out_file << "set title";
            out_file << " \"";
            if (title) {
              out_file << title << " - ";
            }
            out_file << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_DISTRIBUTION]
                     << " " << STAT_label[STATL_FUNCTION] << " " << STAT_label[STATL_MATCHING] << "\"\n\n";

            out_file << "set grid\n" << "set xtics 0,0.1\n" << "set ytics 0,0.1" << endl;

            out_file << "plot [0:" << 1. - complement << "] [0:" << 1. - complement << "] \""
                     << label((data_file_name[2].str()).c_str())
                     << "\" using 2:1 notitle with linespoints,\\" << endl;
            out_file << "\"" << label((data_file_name[2].str()).c_str())
                     << "\" using 2:2 notitle with lines" << endl;

/*            out_file << "plot [0:" << 1. - complement << "] [0:" << 1. - complement << "] \""
                     << label((data_file_name[2].str()).c_str()) << "\" using 2:1 title \""
                     << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " "
                     << STAT_label[STATL_FUNCTION] << "\" with linespoints,\\" << endl;
            out_file << "\"" << label((data_file_name[2].str()).c_str()) << "\" using 2:2 title \""
                     << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_DISTRIBUTION] << " "
                     << STAT_label[STATL_FUNCTION] << "\" with linespoints" << endl; */

            out_file << "unset grid\n" << "set xtics autofreq\n" << "set ytics autofreq" << endl;
          }

          // concentration curves

          if (i == 0) {
            out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
          }
          out_file << endl;

          out_file << "set title";
          if (title) {
            out_file << " \"" << title << "\"";
          }
          out_file << "\n\n";

          out_file << "set grid\n" << "set xtics 0,0.1\n" << "set ytics 0,0.1" << endl;

          out_file << "plot [0:" << 1. - complement << "] [0:" << 1. - complement << "] ";
          if (histo->variance > 0.) {
            out_file << "\"" << label((data_file_name[0].str()).c_str()) << "\" using 5:6 title \""
                     << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " " << STAT_label[STATL_CONCENTRATION] << " "
                     << STAT_label[STATL_CURVE] << "\" with linespoints,\\" << endl;
          }
          out_file << "\"" << label((data_file_name[1].str()).c_str()) << "\" using 3:4 title \""
                   << STAT_label[STATL_DISTRIBUTION] << " " << STAT_label[STATL_CONCENTRATION] << " "
                   << STAT_label[STATL_CURVE] << "\" with linespoints,\\" << endl;
          out_file << "\"" << label((data_file_name[1].str()).c_str())
                   << "\" using 3:3 notitle with lines" << endl;

          out_file << "unset grid\n" << "set xtics autofreq\n" << "set ytics autofreq" << endl;
        }

        if (i == 1) {
          out_file << "\nset terminal x11" << endl;
        }

        out_file << "\npause 0 \"" << STAT_label[STATL_END] << "\"" << endl;
      }
    }

    delete [] data_file_name;

    if (!status) {
      error.update(STAT_error[STATR_FILE_PREFIX]);
    }
  }

  else {
    status = Distribution::plot_write(error , prefix , 0 , NULL , title);
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Plot of a DiscreteParametricModel object using Gnuplot.
 *
 *  \param[in] error  reference on a StatError object,
 *  \param[in] prefix file prefix,
 *  \param[in] title  figure title.
 *
 *  \return           error status.
 */
/*--------------------------------------------------------------*/

bool DiscreteParametricModel::plot_write(StatError &error , const char *prefix ,
                                         const char *title) const

{
  return plot_write(error , prefix , title , frequency_distribution);
}


/*--------------------------------------------------------------*/
/**
 *  \brief Plot of a discrete parametric distribution and the associated frequency distribution.
 *
 *  \param[in] histo pointer on a DiscreteDistributionData object.
 *
 *  \return          MultiPlotSet object.
 */
/*--------------------------------------------------------------*/

MultiPlotSet* DiscreteParametricModel::get_plotable(const DiscreteDistributionData *histo) const

{
  MultiPlotSet *plot_set;


  if (histo) {
    int i , j;
    int nb_plot_set , plot_nb_value;
    double scale , *pcumul;
    ostringstream legend , title;


    nb_plot_set = 1;
    if (variance > 0.) {
      nb_plot_set += 3;
    }
    if (histo->variance == 0.) {
      nb_plot_set--;
    }
 
    plot_set = new MultiPlotSet(nb_plot_set);
    MultiPlotSet &plot = *plot_set;

    title.str("");
    title << STAT_label[STATL_DISTRIBUTION];
    if (histo) {
      title << " " << STAT_label[STATL_FIT];
    }
    plot.title = title.str();

    plot.border = "15 lw 0";

    // distribution fit

    plot_nb_value = plot_nb_value_computation(histo);

    plot[0].xrange = Range(0 , MAX(plot_nb_value - 1 , 1));
    plot[0].yrange = Range(0 , ceil(MAX(histo->max ,
                                    max * histo->nb_element) * YSCALE));

    if (plot_nb_value - 1 < TIC_THRESHOLD) {
      plot[0].xtics = 1;
    }

    plot[0].resize(2);

    plot[0][0].legend = STAT_label[STATL_FREQUENCY_DISTRIBUTION];

    plot[0][0].style = "impulses";

    histo->plotable_frequency_write(plot[0][0]);

    legend.str("");
    legend << STAT_label[STATL_DISTRIBUTION];
    plot_title_print(legend);
    plot[0][1].legend = legend.str();

    plot[0][1].style = "linespoints";

    plotable_mass_write(plot[0][1] , histo->nb_element);

    if (variance > 0.) {

      // computation of the cumulative frequency distribution function

      scale = histo->nb_element / (1. - complement);
      pcumul = histo->cumul_computation(scale);

      // cumulative distribution functions

      plot[1].xrange = Range(0 , plot_nb_value - 1);
      plot[1].yrange = Range(0. , 1. - complement);

      if (plot_nb_value - 1 < TIC_THRESHOLD) {
        plot[1].xtics = 1;
      }

      plot[1].resize(histo->variance > 0. ? 2 : 1);

      if (histo->variance > 0.) {
        legend.str("");
        legend << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " "
               << STAT_label[STATL_FUNCTION];
        plot[1][0].legend = legend.str();

        plot[1][0].style = "linespoints";

        histo->plotable_cumul_write(plot[1][0] , pcumul);
        i = 1;
      }

      else {
        i = 0;
      }

      legend.str("");
      legend << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_DISTRIBUTION] << " "
             << STAT_label[STATL_FUNCTION];
      plot[1][i].legend = legend.str();

      plot[1][i].style = "linespoints";

      plotable_cumul_write(plot[1][i]);

      // matching of cumulative distribution functions taking as reference the distribution one

      if (histo->variance > 0.) {
        title.str("");
        title << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_DISTRIBUTION] << " "
              << STAT_label[STATL_FUNCTION] << " " << STAT_label[STATL_MATCHING];
        plot[2].title = title.str();

        plot[2].xrange = Range(0. , 1. - complement);
        plot[2].yrange = Range(0. , 1. - complement);

        plot[2].grid = true;

        plot[2].xtics = 0.1;
        plot[2].ytics = 0.1;

        plot[2].resize(2);

/*        legend.str("");
        legend << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " "
               << STAT_label[STATL_FUNCTION];
               plot[2][0].legend = legend.str(); */

        plot[2][0].style = "linespoints";

        histo->plotable_cumul_matching_write(plot[2][0] , offset , nb_value , cumul , pcumul);

/*        legend.str("");
        legend << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_DISTRIBUTION] << " "
               << STAT_label[STATL_FUNCTION];
        plot[2][1].legend = legend.str();

        plot[2][1].style = "linespoints";

        plotable_cumul_matching_write(plot[2][1] , *this); */

        plot[2][1].style = "lines";

        plot[2][1].add_point(0. , 0.);
        plot[2][1].add_point(1. - complement , 1. - complement);
      }

      // concentration curves

      if (histo->variance == 0.) {
        i = 2;
      }
      else {
        i = 3;
      }

      plot[i].xrange = Range(0. , 1. - complement);
      plot[i].yrange = Range(0. , 1. - complement);

      plot[i].grid = true;

      plot[i].xtics = 0.1;
      plot[i].ytics = 0.1;

      plot[i].resize(histo->variance > 0. ? 3 : 2);

      if (histo->variance > 0.) {
        legend.str("");
        legend << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " " << STAT_label[STATL_CONCENTRATION] << " "
               << STAT_label[STATL_CURVE];
        plot[i][0].legend = legend.str();

        plot[i][0].style = "linespoints";

        histo->plotable_concentration_write(plot[i][0] , pcumul , scale);
        j = 1;
      }

      else {
        j = 0;
      }

      legend.str("");
      legend << STAT_label[STATL_DISTRIBUTION] << " " << STAT_label[STATL_CONCENTRATION] << " "
             << STAT_label[STATL_CURVE];
      plot[i][j].legend = legend.str();

      plot[i][j].style = "linespoints";

      plotable_concentration_write(plot[i][j]);
      j++;

      plot[i][j].style = "lines";

      plot[i][j].add_point(0. , 0.);
      plot[i][j].add_point(1. - complement , 1. - complement);
    }
  }

  else {
    StatError error;

    plot_set = Distribution::get_plotable_distributions(error , 0 , NULL);
  }

  return plot_set;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Plot of a DiscreteParametricModel object.
 *
 *  \return MultiPlotSet object.
 */
/*--------------------------------------------------------------*/

MultiPlotSet* DiscreteParametricModel::get_plotable() const

{
  return get_plotable(frequency_distribution);
}


};  // namespace stat_tool
