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
 *       $Id: mixture.cpp 15005 2013-09-30 14:23:12Z guedon $
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

#include "mixture.h"
#include "stat_label.h"

using namespace std;
using namespace boost;


namespace stat_tool {



/*--------------------------------------------------------------*/
/**
 *  \brief Default constructor of the Mixture class.
 */
/*--------------------------------------------------------------*/

Mixture::Mixture()

{
  mixture_data = NULL;

  nb_component = 0;
  weight = NULL;

  nb_output_process = 0;
  categorical_process = NULL;
  discrete_parametric_process = NULL;
  continuous_parametric_process = NULL;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Constructor of the Mixture class.
 *
 *  \param[in] inb_component      number of components,
 *  \param[in] inb_output_process number of observation processes,
 *  \param[in] nb_value           number of observed values per process.
 */
/*--------------------------------------------------------------*/

Mixture::Mixture(int inb_component , int inb_output_process , int *nb_value)

{
  int i;


  mixture_data = NULL;

  nb_component = inb_component;
  weight = new DiscreteParametric(nb_component);

  nb_output_process = inb_output_process;

  categorical_process = new CategoricalProcess*[nb_output_process];
  discrete_parametric_process = new DiscreteParametricProcess*[nb_output_process];
  continuous_parametric_process = new ContinuousParametricProcess*[nb_output_process];

  for (i = 0;i < nb_output_process;i++) {
    if (nb_value[i] == I_DEFAULT) {
      categorical_process[i] = NULL;
      discrete_parametric_process[i] = NULL;
      continuous_parametric_process[i] = new ContinuousParametricProcess(nb_component);
    }

    else if (nb_value[i] <= NB_OUTPUT) {
      categorical_process[i] = new CategoricalProcess(nb_component , nb_value[i] , true);
      discrete_parametric_process[i] = NULL;
      continuous_parametric_process[i] = NULL;
    }

    else {
      categorical_process[i] = NULL;
      discrete_parametric_process[i] = new DiscreteParametricProcess(nb_component , (int)(nb_value[i] * SAMPLE_NB_VALUE_COEFF));
      continuous_parametric_process[i] = NULL;
    }
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Constructor of the Mixture class (univariate Gaussian with evenly spaced means).
 *
 *  \param[in] inb_component      number of components,
 *  \param[in] offset             mixture offset,
 *  \param[in] mean               mean - offset of the 1st component,
 *  \param[in] standard_deviation standard deviation,
 *  \param[in] common_dispersion  common dispersion parameter,
 */
/*--------------------------------------------------------------*/

Mixture::Mixture(int inb_component , double offset , double mean , double standard_deviation ,
                 bool common_dispersion)

{
  int i , j;
  ContinuousParametric **observation;


  mixture_data = NULL;

  nb_component = inb_component;
  weight = new DiscreteParametric(nb_component);

  nb_output_process = 1;

  categorical_process = new CategoricalProcess*[1];
  discrete_parametric_process = new DiscreteParametricProcess*[1];
  continuous_parametric_process = new ContinuousParametricProcess*[1];

  categorical_process[0] = NULL;
  discrete_parametric_process[0] = NULL;

  observation = new ContinuousParametric*[nb_component];

  for (i = 0;i < nb_component - 1;i++) {
    weight->mass[i] = 1 / (double)nb_component;
    observation[i] = new ContinuousParametric(GAUSSIAN , offset + (i + 1) * mean , standard_deviation);
  }
  weight->mass[nb_component - 1] = 1 / (double)nb_component;
  observation[nb_component - 1] = new ContinuousParametric(GAUSSIAN , offset + nb_component * mean / 2 , 10 * standard_deviation);

  continuous_parametric_process[0] = new ContinuousParametricProcess(nb_component , observation);
  continuous_parametric_process[0]->tied_location = true;
  continuous_parametric_process[0]->tied_dispersion = common_dispersion;
  continuous_parametric_process[0]->offset = offset;

  delete [] observation;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Constructor of the Mixture class (univariate gamma, inverse Gaussian or Gaussian with tied parameters).
 *
 *  \param[in] inb_component      number of components,
 *  \param[in] ident              component identifiers,
 *  \param[in] mean               mean of the 1st component,
 *  \param[in] standard_deviation shape parameter (gamma) / scale parameter (inverse Gaussian) /
 *                                standard deviation (Gaussian) of the 1st component,
 *  \param[in] tied_mean          flag tied means,
 *  \param[in] variance_factor    type of link between the variances (CONVOLUTION_FACTOR/SCALING_FACTOR).
 */
/*--------------------------------------------------------------*/

Mixture::Mixture(int inb_component , int ident , double mean , double standard_deviation ,
                 bool tied_mean , tying_rule variance_factor)

{
  int i , j;
  ContinuousParametric **observation;


  mixture_data = NULL;

  nb_component = inb_component;
  weight = new DiscreteParametric(nb_component);

  nb_output_process = 1;

  categorical_process = new CategoricalProcess*[1];
  discrete_parametric_process = new DiscreteParametricProcess*[1];
  continuous_parametric_process = new ContinuousParametricProcess*[1];

  categorical_process[0] = NULL;
  discrete_parametric_process[0] = NULL;

  observation = new ContinuousParametric*[nb_component];

  i = 1;
  for (j = 0;j < nb_component;j++) {
    weight->mass[j] = 1 / (double)nb_component;

    switch (ident) {

    case GAMMA : {
      switch (variance_factor) {
      case CONVOLUTION_FACTOR :
        observation[j] = new ContinuousParametric(GAMMA , i * standard_deviation , mean / standard_deviation);
        break;
      case SCALING_FACTOR :
        observation[j] = new ContinuousParametric(GAMMA , standard_deviation , i * mean / standard_deviation);
        break;
      }
      break;
    }

    case INVERSE_GAUSSIAN : {
      switch (variance_factor) {
      case CONVOLUTION_FACTOR :
        observation[j] = new ContinuousParametric(INVERSE_GAUSSIAN , i * mean ,  i * i * standard_deviation);
        break;
      case SCALING_FACTOR :
        observation[j] = new ContinuousParametric(INVERSE_GAUSSIAN , i * mean , i * standard_deviation);
        break;
      }
      break;
    }

    case GAUSSIAN : {
      switch (variance_factor) {
      case CONVOLUTION_FACTOR :
        observation[j] = new ContinuousParametric(GAUSSIAN , i * mean , sqrt((double)i) * standard_deviation);
        break;
      case SCALING_FACTOR :
        observation[j] = new ContinuousParametric(GAUSSIAN , i * mean , i * standard_deviation);
        break;
      }
      break;
    }
    }
    i *= 2;
  }

  continuous_parametric_process[0] = new ContinuousParametricProcess(nb_component , observation);
  continuous_parametric_process[0]->tied_location = tied_mean;
  continuous_parametric_process[0]->tied_dispersion = true;

  delete [] observation;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Constructor of the Mixture class.
 *
 *  \param[in] iweight                           weight distribution,
 *  \param[in] inb_output_process                number of observation processes,
 *  \param[in] categorical_observation           pointer on CategoricalProcess objects,
 *  \param[in] discrete_parametric_observation   pointer on DiscreteParametricProcess objects,
 *  \param[in] continuous_parametric_observation pointer on ContinuousParametricProcess objects.
 */
/*--------------------------------------------------------------*/

Mixture::Mixture(const DiscreteParametric *iweight , int inb_output_process ,
                 CategoricalProcess **categorical_observation ,
                 DiscreteParametricProcess **discrete_parametric_observation ,
                 ContinuousParametricProcess **continuous_parametric_observation)

{
  int i;


  mixture_data = NULL;

  nb_component = iweight->nb_value;
  weight = new DiscreteParametric(*iweight);

  nb_output_process = inb_output_process;

  categorical_process = new CategoricalProcess*[nb_output_process];
  discrete_parametric_process = new DiscreteParametricProcess*[nb_output_process];
  continuous_parametric_process = new ContinuousParametricProcess*[nb_output_process];

  for (i = 0;i < nb_output_process;i++) {
    if (categorical_observation[i]) {
      categorical_process[i] = new CategoricalProcess(*categorical_observation[i]);
      discrete_parametric_process[i] = NULL;
      continuous_parametric_process[i] = NULL;
    }
    else if (discrete_parametric_observation[i]) {
      categorical_process[i] = NULL;
      discrete_parametric_process[i] = new DiscreteParametricProcess(*discrete_parametric_observation[i]);
      continuous_parametric_process[i] = NULL;
    }
    else {
      categorical_process[i] = NULL;
      discrete_parametric_process[i] = NULL;
      continuous_parametric_process[i] = new ContinuousParametricProcess(*continuous_parametric_observation[i]);
    }
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Copy of a Mixture object.
 *
 *  \param[in] mixt      reference on a Mixture object,
 *  \param[in] data_flag flag copy of the included MixtureData object.
 */
/*--------------------------------------------------------------*/

void Mixture::copy(const Mixture &mixt , bool data_flag)

{
  int i;


  if ((data_flag) && (mixt.mixture_data)) {
    mixture_data = new MixtureData(*(mixt.mixture_data) , false);
  }
  else {
    mixture_data = NULL;
  }

  nb_component = mixt.nb_component;

  weight = new DiscreteParametric(*(mixt.weight));

  nb_output_process = mixt.nb_output_process;

  if (mixt.categorical_process) {
    categorical_process = new CategoricalProcess*[nb_output_process];

    for (i = 0;i < nb_output_process;i++) {
      if (mixt.categorical_process[i]) {
        categorical_process[i] = new CategoricalProcess(*(mixt.categorical_process[i]));
      }
      else {
        categorical_process[i] = NULL;
      }
    }
  }

  else {
    categorical_process = NULL;
  }

  if (mixt.discrete_parametric_process) {
    discrete_parametric_process = new DiscreteParametricProcess*[nb_output_process];

    for (i = 0;i < nb_output_process;i++) {
      if (mixt.discrete_parametric_process[i]) {
        discrete_parametric_process[i] = new DiscreteParametricProcess(*(mixt.discrete_parametric_process[i]));
      }
      else {
        discrete_parametric_process[i] = NULL;
      }
    }
  }

  else {
    discrete_parametric_process = NULL;
  }

  if (mixt.continuous_parametric_process) {
    continuous_parametric_process = new ContinuousParametricProcess*[nb_output_process];

    for (i = 0;i < nb_output_process;i++) {
      if (mixt.continuous_parametric_process[i]) {
        continuous_parametric_process[i] = new ContinuousParametricProcess(*(mixt.continuous_parametric_process[i]));
      }
      else {
        continuous_parametric_process[i] = NULL;
      }
    }
  }

  else {
    continuous_parametric_process = NULL;
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Destruction of the data members of a Mixture object.
 */
/*--------------------------------------------------------------*/

void Mixture::remove()

{
  int i;


  delete mixture_data;

  delete weight;

  if (categorical_process) {
    for (i = 0;i < nb_output_process;i++) {
      delete categorical_process[i];
    }
    delete [] categorical_process;
  }

  if (discrete_parametric_process) {
    for (i = 0;i < nb_output_process;i++) {
      delete discrete_parametric_process[i];
    }
    delete [] discrete_parametric_process;
  }

  if (continuous_parametric_process) {
    for (i = 0;i < nb_output_process;i++) {
      delete continuous_parametric_process[i];
    }
    delete [] continuous_parametric_process;
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Destructor of the Mixture class.
 */
/*--------------------------------------------------------------*/

Mixture::~Mixture()

{
  remove();
}


/*--------------------------------------------------------------*/
/**
 *  \brief Assignment operator of the Mixture class.
 *
 *  \param[in] mixt reference on a Mixture object.
 *
 *  \return         Mixture object.
 */
/*--------------------------------------------------------------*/

Mixture& Mixture::operator=(const Mixture &mixt)

{
  if (&mixt != this) {
    remove();
    copy(mixt);
  }

  return *this;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Extraction of a discrete parametric component.
 *
 *  \param[in] error    reference on a StatError object,
 *  \param[in] variable variable index,
 *  \param[in] index    component index.
 *
 *  \return             DiscreteParametricModel object.
 */
/*--------------------------------------------------------------*/

DiscreteParametricModel* Mixture::extract(StatError &error , int variable , int index) const

{
  bool status = true;
  int hvariable;
  Distribution *pdist;
  DiscreteParametric *pparam;
  DiscreteParametricModel *dist;
  FrequencyDistribution *phisto;
  CategoricalProcess *process;


  dist = NULL;
  error.init();

  if ((variable < 1) || (variable > nb_output_process)) {
    status = false;
    error.update(STAT_error[STATR_OUTPUT_PROCESS_INDEX]);
  }

  else {
    if ((index < 1) || (index > nb_component)) {
      status = false;
      ostringstream error_message;
      error_message << STAT_label[STATL_COMPONENT] << " " << index << " "
                    << STAT_error[STATR_NOT_PRESENT];
      error.update((error_message.str()).c_str());
    }

    else {
      index--;

      if (categorical_process[variable - 1]) {
        pdist = categorical_process[variable - 1]->observation[index];
        pparam = NULL;
      }
      else if (discrete_parametric_process[variable - 1]) {
        pdist = NULL;
        pparam = discrete_parametric_process[variable - 1]->observation[index];
      }
      else {
        status = false;
        ostringstream correction_message;
        correction_message << STAT_label[STATL_CATEGORICAL] << " or "
                           << STAT_label[STATL_DISCRETE_PARAMETRIC];
        error.correction_update(STAT_error[STATR_OUTPUT_PROCESS_TYPE] , (correction_message.str()).c_str());
      }
    }
  }

  if (status) {
    phisto = NULL;

    if (mixture_data) {
      switch (mixture_data->type[0]) {
      case STATE :
        hvariable = variable;
        break;
      case INT_VALUE :
        hvariable = variable - 1;
        break;
      }

      if (hvariable >= 0) {
        if ((mixture_data->observation_distribution) &&
            (mixture_data->observation_distribution[hvariable])) {
          phisto = mixture_data->observation_distribution[hvariable][index];
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
 *  \brief Extraction of the MixtureData object included in a Mixture object.
 *
 *  \param[in] error reference on a StatError object.
 *
 *  \return          MixtureData object.
 */
/*--------------------------------------------------------------*/

MixtureData* Mixture::extract_data(StatError &error) const

{
  bool status = true;
  MixtureData *vec;


  vec = NULL;
  error.init();

  if (!mixture_data) {
    status = false;
    error.update(STAT_error[STATR_NO_DATA]);
  }
  else if (nb_output_process + 1 != mixture_data->nb_variable) {
    status = false;
    error.update(STAT_error[STATR_OPTIMAL_ASSIGNMENT]);
  }

  if (status) {
    vec = new MixtureData(*mixture_data);
    vec->mixture = new Mixture(*this , false);
  }

  return vec;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Application of a threshold on the parameters of a multivariate mixture.
 *
 *  \param[in] min_probability minimum probability.
 *
 *  \return                    Mixture object.
 */
/*--------------------------------------------------------------*/

Mixture* Mixture::thresholding(double min_probability) const

{
  int i;
  Mixture *mixt;


  mixt = new Mixture(*this , false);

  for (i = 0;i < mixt->nb_output_process;i++) {
    if (mixt->categorical_process[i]) {
      mixt->categorical_process[i]->thresholding(min_probability);
    }
  }

  return mixt;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Construction of a Mixture object from a file.
 *
 *  \param[in] error           reference on a StatError object,
 *  \param[in] path            file path,
 *  \param[in] cumul_threshold threshold on the cumulative distribution function.
 *
 *  \return                    Mixture object.
 */
/*--------------------------------------------------------------*/

Mixture* Mixture::ascii_read(StatError &error , const string path , double cumul_threshold)

{
  string buffer;
  size_t position;
  typedef tokenizer<char_separator<char>> tokenizer;
  char_separator<char> separator(" \t");
  bool status , lstatus;
  int i;
  int line , read_line , nb_component , value , nb_output_process , index;
  observation_process obs_type;
  double proba , cumul;
  DiscreteParametric *weight;
  CategoricalProcess **categorical_observation;
  DiscreteParametricProcess **discrete_parametric_observation;
  ContinuousParametricProcess **continuous_parametric_observation;
  Mixture *mixt;
  ifstream in_file(path.c_str());


  mixt = NULL;
  error.init();

  if (!in_file) {
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {
    status = true;
    line = 0;
    nb_component = 0;

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

        // test MIXTURE keyword

        case 0 : {
          if (*token != STAT_word[STATW_MIXTURE]) {
            status = false;
            error.correction_update(STAT_parsing[STATP_KEYWORD] , STAT_word[STATW_MIXTURE] , line , i + 1);
          }
          break;
        }

        // test number of components

        case 1 : {
          lstatus = true;

/*          try {
            nb_component = stoi(*token);   in C++ 11
          }
          catch(invalid_argument &arg) {
            lstatus = false;
          } */
          nb_component = atoi(token->c_str());

          if ((lstatus) && ((nb_component < 2) || (nb_component > MIXTURE_NB_COMPONENT))) {
            lstatus = false;
          }

          if (!lstatus) {
            status = false;
            error.update(STAT_parsing[STATP_NB_COMPONENT] , line , i + 1);
          }
          break;
        }

        // test COMPONENTS keyword

        case 2 : {
          if (*token != STAT_word[STATW_COMPONENTS]) {
            status = false;
            error.correction_update(STAT_parsing[STATP_KEYWORD] , STAT_word[STATW_COMPONENTS] , line , i + 1);
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

    if (nb_component == 0) {
      status = false;
      error.update(STAT_parsing[STATP_FORMAT] , line);
    }

    if (status) {
      weight = new DiscreteParametric(nb_component);

      // format analysis and reading of the weight distribution

      read_line = 0;
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

        if (read_line == 0) {
          for (tokenizer::iterator token = tok_buffer.begin();token != tok_buffer.end();token++) {

            // test WEIGHTS keyword

            if (i == 0) {
              if (*token != STAT_word[STATW_WEIGHTS]) {
                status = false;
                error.correction_update(STAT_parsing[STATP_KEYWORD] , STAT_word[STATW_WEIGHTS] , line);
              }
            }

            i++;
          }

          if (i > 0) {
            if (i != 1) {
              status = false;
              error.update(STAT_parsing[STATP_FORMAT] , line);
            }
          }
        }

        else {
          cumul = 0.;

          for (tokenizer::iterator token = tok_buffer.begin();token != tok_buffer.end();token++) {
            if (i < nb_component) {
              lstatus = true;

/*              try {
                proba = stod(*token);   in C++ 11
              }
              catch (invalid_argument &arg) {
                lstatus = false;
              } */
              proba = atof(token->c_str());

              if (lstatus) {
                if ((proba <= 0.) || (proba > 1. - cumul + DOUBLE_ERROR)) {
                  lstatus = false;
                }

                else {
                  cumul += proba;
                  weight->mass[i] = proba;
                }
              }

              if (!lstatus) {
                status = false;
                error.update(STAT_parsing[STATP_WEIGHT_VALUE] , line , i + 1);
              }
            }

            i++;
          }

          if (i > 0) {
            if (i != nb_component) {
              status = false;
              error.update(STAT_parsing[STATP_FORMAT] , line);
            }

            if (cumul < 1. - DOUBLE_ERROR) {
              status = false;
              error.update(STAT_parsing[STATP_PROBABILITY_SUM] , line);
            }
          }
        }

        if (i > 0) {
          read_line++;
          if (read_line == 2) {
            break;
          }
        }
      }

      // format analysis and reading of the observation distributions

      if (status) {
        nb_output_process = I_DEFAULT;

        categorical_observation = NULL;
        discrete_parametric_observation = NULL;
        continuous_parametric_observation = NULL;

        while (getline(in_file , buffer)) {
          line++;

#         ifdef DEBUG
          cout << line << "  " << buffer << endl;
#         endif

          position = buffer.find('#');
          if (position != string::npos) {
            buffer.erase(position);
          }
          i = 0;

          tokenizer tok_buffer(buffer , separator);

          for (tokenizer::iterator token = tok_buffer.begin();token != tok_buffer.end();token++) {
            switch (i) {

            // test number of observation processes

            case 0 : {
              lstatus = true;

/*              try {
                value = stoi(*token);   in C++ 11
              }
              catch(invalid_argument &arg) {
                lstatus = false;
              } */
              value = atoi(token->c_str());

              if (lstatus) {
                if ((value < 1) || (value > NB_OUTPUT_PROCESS)) {
                  lstatus = false;
                }
                else {
                  nb_output_process = value;
                }
              }

              if (!lstatus) {
                status = false;
                error.update(STAT_parsing[STATP_NB_OUTPUT_PROCESS] , line , i + 1);
              }
              break;
            }

            // test OUTPUT_PROCESS(ES) keyword

            case 1 : {
              if (*token != STAT_word[nb_output_process == 1 ? STATW_OUTPUT_PROCESS : STATW_OUTPUT_PROCESSES]) {
                status = false;
                error.correction_update(STAT_parsing[STATP_KEYWORD] ,
                                        STAT_word[nb_output_process == 1 ? STATW_OUTPUT_PROCESS : STATW_OUTPUT_PROCESSES] , line , i + 1);
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
            break;
          }
        }

        if (nb_output_process == I_DEFAULT) {
          status = false;
          error.update(STAT_parsing[STATP_FORMAT] , line);
        }

        else {
          categorical_observation = new CategoricalProcess*[nb_output_process];
          discrete_parametric_observation = new DiscreteParametricProcess*[nb_output_process];
          continuous_parametric_observation = new ContinuousParametricProcess*[nb_output_process];

          for (i = 0;i < nb_output_process;i++) {
            categorical_observation[i] = NULL;
            discrete_parametric_observation[i] = NULL;
            continuous_parametric_observation[i] = NULL;
          }

          index = 0;

          while (getline(in_file , buffer)) {
            line++;

#           ifdef DEBUG
            cout << line << "  " << buffer << endl;
#           endif

            position = buffer.find('#');
            if (position != string::npos) {
              buffer.erase(position);
            }
            i = 0;

            tokenizer tok_buffer(buffer , separator);

            for (tokenizer::iterator token = tok_buffer.begin();token != tok_buffer.end();token++) {
              switch (i) {

              // test OUTPUT_PROCESS keyword

              case 0 : {
                if (*token != STAT_word[STATW_OUTPUT_PROCESS]) {
                  status = false;
                  error.correction_update(STAT_parsing[STATP_KEYWORD] , STAT_word[STATW_OUTPUT_PROCESS] , line , i + 1);
                }
                break;
              }

              // test observation process index

              case 1 : {
                index++;
                lstatus = true;

/*                try {
                  value = stoi(*token);   in C++ 11
                }
                catch(invalid_argument &arg) {
                  lstatus = false;
                } */
                value = atoi(token->c_str());

                if ((lstatus) && ((value != index) || (value > nb_output_process))) {
                  lstatus = false;
                }

                if (!lstatus) {
                  status = false;
                  error.update(STAT_parsing[STATP_OUTPUT_PROCESS_INDEX] , line , i + 1);
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

              // test CATEGORICAL/DISCRETE_PARAMETRIC/CONTINUOUS_PARAMETRIC keyword

              case 3 : {
                if ((*token == STAT_word[STATW_CATEGORICAL]) ||
                    (*token == STAT_word[STATW_NONPARAMETRIC])) {
                  obs_type = CATEGORICAL_PROCESS;
                }
                else if ((*token == STAT_word[STATW_DISCRETE_PARAMETRIC]) ||
                         (*token == STAT_word[STATW_PARAMETRIC])) {
                  obs_type = DISCRETE_PARAMETRIC;
                }
                else if (*token == STAT_word[STATW_CONTINUOUS_PARAMETRIC]) {
                  obs_type = CONTINUOUS_PARAMETRIC;
                }
                else {
                  obs_type = DEFAULT_PROCESS;
                  status = false;
                  ostringstream correction_message;
                  correction_message << STAT_word[STATW_CATEGORICAL] << " or "
                                     << STAT_word[STATW_DISCRETE_PARAMETRIC] << " or "
                                     << STAT_word[STATW_CONTINUOUS_PARAMETRIC];
                  error.correction_update(STAT_parsing[STATP_KEYWORD] , (correction_message.str()).c_str() , line , i + 1);
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

              switch (obs_type) {

              case CATEGORICAL_PROCESS : {
                categorical_observation[index - 1] = CategoricalProcess::parsing(error , in_file , line ,
                                                                                 nb_component , MIXTURE , true);
                if (!categorical_observation[index - 1]) {
                  status = false;
                }
                break;
              }

              case DISCRETE_PARAMETRIC : {
                discrete_parametric_observation[index - 1] = DiscreteParametricProcess::parsing(error , in_file , line ,
                                                                                                nb_component , MIXTURE ,
                                                                                                cumul_threshold);
                if (!discrete_parametric_observation[index - 1]) {
                  status = false;
                }
                break;
              }

              case CONTINUOUS_PARAMETRIC : {
                continuous_parametric_observation[index - 1] = ContinuousParametricProcess::parsing(error , in_file , line ,
                                                                                                    nb_component , MIXTURE ,
                                                                                                    VON_MISES);
                if (!continuous_parametric_observation[index - 1]) {
                  status = false;
                }
                break;
              }
              }
            }

            if (index == nb_output_process) {
              break;
            }
          }

          if (index < nb_output_process) {
            status = false;
            error.update(STAT_parsing[STATP_FORMAT] , line);
          }

          else {
            while (getline(in_file , buffer)) {
              line++;

#             ifdef DEBUG
              cout << line << " " << buffer << endl;
#             endif

              position = buffer.find('#');
              if (position != string::npos) {
                buffer.erase(position);
              }
              if (!(trim_right_copy_if(buffer , is_any_of(" \t")).empty())) {
                status = false;
                error.update(STAT_parsing[STATP_FORMAT] , line);
              }
            }
          }

          if (status) {
            weight->max_computation();
            weight->cumul_computation();
            mixt = new Mixture(weight , nb_output_process , categorical_observation ,
                               discrete_parametric_observation , continuous_parametric_observation);

#           ifdef DEBUG
            mixt->ascii_write(cout);
#           endif

          }

          delete weight;

          for (i = 0;i < nb_output_process;i++) {
            delete categorical_observation[i];
            delete discrete_parametric_observation[i];
            delete continuous_parametric_observation[i];
          }
          delete [] categorical_observation;
          delete [] discrete_parametric_observation;
          delete [] continuous_parametric_observation;
        }
      }
    }
  }

  return mixt;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing on a single line of a Mixture object.
 *
 *  \param[in,out] os stream.
 */
/*--------------------------------------------------------------*/

ostream& Mixture::line_write(ostream &os) const

{
  os << nb_component << " " << STAT_word[STATW_COMPONENTS];

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a Mixture object and the associated data structure.
 *
 *  \param[in,out] os         stream,
 *  \param[in]     vec        pointer on a MixtureData object,
 *  \param[in]     exhaustive flag detail level,
 *  \param[in]     file_flag  flag file.
 */
/*--------------------------------------------------------------*/

ostream& Mixture::ascii_write(ostream &os , const MixtureData *vec ,
                               bool exhaustive , bool file_flag) const

{
  int i , j , k;
  int buff , width , variable;
  double mean , **distance;
  FrequencyDistribution *marginal_dist = NULL , **observation_dist = NULL;
  Histogram *marginal_histo = NULL , **observation_histo = NULL;
  ios_base::fmtflags format_flags;


  format_flags = os.setf(ios::left , ios::adjustfield);

  os << STAT_word[STATW_MIXTURE] << " " << nb_component << " " << STAT_word[STATW_COMPONENTS] << endl;

  // writing of the weights

  os << "\n" << STAT_word[STATW_WEIGHTS] << endl;

  width = column_width(nb_component , weight->mass) + 1;

  for (i = 0;i < nb_component;i++) {
    os << setw(width) << weight->mass[i];
  }
  os << endl;

  if (weight->mean != D_DEFAULT) {
    if (file_flag) {
      os << "# ";
    }
    os << STAT_label[STATL_MEAN] << ": " << weight->mean << endl;
  }

  os << "\n" << nb_output_process << " "
     << STAT_word[nb_output_process == 1 ? STATW_OUTPUT_PROCESS : STATW_OUTPUT_PROCESSES] << endl;

  // writing of the distributions associated with each observation process

  distance = new double*[nb_component];
  for (i = 0;i < nb_component;i++) {
    distance[i] = new double[nb_component];
  }

  for (i = 0;i < nb_output_process;i++) {
    os << "\n" << STAT_word[STATW_OUTPUT_PROCESS] << " " << i + 1;

    if (categorical_process[i]) {
      os << " : " << STAT_word[STATW_CATEGORICAL];
    }
    else if (discrete_parametric_process[i]) {
      os << " : " << STAT_word[STATW_DISCRETE_PARAMETRIC];
    }
    else {
      os << " : " << STAT_word[STATW_CONTINUOUS_PARAMETRIC];
    }
    os << endl;

    if (vec) {
      switch (vec->type[0]) {
      case STATE :
        variable = i + 1;
        break;
      default :
        variable = i;
        break;
      }

      if (vec->observation_distribution) {
        observation_dist = vec->observation_distribution[variable];
      }
      marginal_dist = vec->marginal_distribution[variable];

      if (vec->observation_histogram) {
        observation_histo = vec->observation_histogram[variable];
      }
      marginal_histo = vec->marginal_histogram[variable];
    }

    if (categorical_process[i]) {
      categorical_process[i]->ascii_print(os , observation_dist , marginal_dist ,
                                          exhaustive , file_flag , MIXTURE);

      for (j = 0;j < nb_component;j++) {
        distance[j][j] = 0.;
        for (k = j + 1;k < nb_component;k++) {
          distance[j][k] = categorical_process[i]->observation[j]->overlap_distance_computation(*(categorical_process[i]->observation[k]));
          distance[k][j] = distance[j][k];
        }
      }
    }

    else if (discrete_parametric_process[i]) {
      discrete_parametric_process[i]->ascii_print(os , observation_dist , marginal_dist ,
                                                  exhaustive , file_flag , MIXTURE);

      for (j = 0;j < nb_component;j++) {
        distance[j][j] = 0.;
        for (k = j + 1;k < nb_component;k++) {
          distance[j][k] = discrete_parametric_process[i]->observation[j]->sup_norm_distance_computation(*(discrete_parametric_process[i]->observation[k]));
          distance[k][j] = distance[j][k];
        }
      }
    }

    else {
      continuous_parametric_process[i]->ascii_print(os , observation_histo , observation_dist ,
                                                    marginal_histo , marginal_dist ,
                                                    exhaustive , file_flag , MIXTURE);

      if ((vec) && (!marginal_dist)) {
        mean = continuous_parametric_process[i]->mean_computation(weight);

        os << "\n";
        if (file_flag) {
          os << "# ";
        }
        os << STAT_label[STATL_MIXTURE] << " - " << STAT_label[STATL_MEAN] << ": " << mean << "  "
           << STAT_label[STATL_STANDARD_DEVIATION] << ": "
           << sqrt(continuous_parametric_process[i]->variance_computation(weight , mean)) << endl;

        if (file_flag) {
          os << "# ";
        }
        os << STAT_label[STATL_VARIABLE] << " " << variable + 1 << "   "
           << STAT_label[STATL_MEAN] << ": " << vec->mean[variable] << "   "
           << STAT_label[STATL_STANDARD_DEVIATION] << ": " << sqrt(vec->covariance[variable][variable]) << endl;
      }

      for (j = 0;j < nb_component;j++) {
        distance[j][j] = 0.;
        for (k = j + 1;k < nb_component;k++) {
          distance[j][k] = continuous_parametric_process[i]->observation[j]->sup_norm_distance_computation(*(continuous_parametric_process[i]->observation[k]));
          distance[k][j] = distance[j][k];
        }
      }
    }

    width = column_width(nb_component , distance[0]);
    for (j = 1;j < nb_component;j++) {
      buff = column_width(nb_component , distance[j]);
      if (buff > width) {
        width = buff;
      }
    }
    width += ASCII_SPACE;

    os.setf(ios::left , ios::adjustfield);

    os << "\n";
    if (file_flag) {
      os << "# ";
    }
//    os << STAT_label[STATL_COMPONENT_DISTANCE] << endl;
    os << STAT_label[STATL_OBSERVATION_DISTRIBUTION_DISTANCE] << endl;

    for (j = 0;j < nb_component;j++) {
      if (file_flag) {
        os << "# ";
      }
      for (k = 0;k < nb_component;k++) {
        if (k != j) {
          os << setw(width) << distance[j][k];
        }
        else {
          os << setw(width) << "_";
        }
      }
      os << endl;
    }
  }

  for (i = 0;i < nb_component;i++) {
    delete [] distance[i];
  }
  delete [] distance;

  if (vec) {
    int nb_parameter = nb_parameter_computation(MIN_PROBABILITY);
    double information;


    os << "\n";
    if (file_flag) {
      os << "# ";
    }
    os << STAT_label[STATL_NB_VECTOR] << ": " << vec->nb_vector << endl;

    // writing of the information quantity of the vectors in the iid case

    for (i = 0;i < vec->nb_variable;i++) {
      if (vec->type[i] == REAL_VALUE) {
        break;
      }
    }

    if (i == vec->nb_variable) {
      information = vec->information_computation();

      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << STAT_label[STATL_INFORMATION] << ": " << information << " ("
         << information / vec->nb_vector << ")" << endl;
    }

    // writing of the log-likelihoods for the observed vectors

    if (vec->restoration_likelihood != D_INF) {
      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << STAT_label[STATL_CLASSIFICATION_LIKELIHOOD] << ": " << vec->restoration_likelihood << "   ("
         << STAT_label[STATL_NORMALIZED] << ": " << vec->restoration_likelihood / vec->nb_vector << ")" << endl;
    }

    if (vec->sample_entropy != D_DEFAULT) {
      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << STAT_label[STATL_CLASSIFICATION_ENTROPY] << ": " << vec->sample_entropy << "   ("
         << STAT_label[STATL_NORMALIZED] << ": " << vec->sample_entropy / vec->nb_vector << ")" << endl;
    }

    if (vec->likelihood != D_INF) {
      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << STAT_label[STATL_LIKELIHOOD] << ": " << vec->likelihood << "   ("
         << STAT_label[STATL_NORMALIZED] << ": " << vec->likelihood / vec->nb_vector << ")" << endl;
    }

    if (vec->likelihood != D_INF) {
      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << nb_parameter << " " << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS]
         << "   2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[AIC] << "): "
         << 2 * (vec->likelihood - nb_parameter) << endl;

      if (nb_parameter < vec->nb_vector - 1) {
        os << "\n";
        if (file_flag) {
          os << "# ";
        }
        os << nb_parameter << " " << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS]
           << "   2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[AICc] << "): "
           << 2 * (vec->likelihood - (double)(nb_parameter * vec->nb_vector) /
             (double)(vec->nb_vector - nb_parameter - 1)) << endl;
      }

      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << nb_parameter << " " << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS]
         << "   2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[BIC] << "): "
         << 2 * vec->likelihood - nb_parameter * log((double)vec->nb_vector) << endl;

      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << nb_parameter << " " << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS]
         << "   2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[ICL] << "): "
//         << 2 * vec->restoration_likelihood - nb_parameter * log((double)vec->nb_vector) << endl;
         << 2 * (vec->likelihood - vec->sample_entropy) - nb_parameter * log((double)vec->nb_vector) << endl;
    }
  }

  os.setf(format_flags , ios::adjustfield);

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a Mixture object.
 *
 *  \param[in,out] os         stream,
 *  \param[in]     exhaustive flag detail level.
 */
/*--------------------------------------------------------------*/

ostream& Mixture::ascii_write(ostream &os , bool exhaustive) const

{
  return ascii_write(os , mixture_data , exhaustive , false);
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a Mixture object in a file.
 *
 *  \param[in] error      reference on a StatError object,
 *  \param[in] path       file path,
 *  \param[in] exhaustive flag detail level.
 *
 *  \return               error status.
 */
/*--------------------------------------------------------------*/

bool Mixture::ascii_write(StatError &error , const string path ,
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
    ascii_write(out_file , mixture_data , exhaustive , true);
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a Mixture object and the associated data structure
 *         in a file at the spreadsheet format.
 *
 *  \param[in,out] os  stream,
 *  \param[in]     vec pointer on a MixtureData object.
 */
/*--------------------------------------------------------------*/

ostream& Mixture::spreadsheet_write(ostream &os , const MixtureData *vec) const

{
  int i , j , k;
  int variable;
  double **distance;
  FrequencyDistribution *marginal_dist = NULL , **observation_dist = NULL;
  Histogram *marginal_histo = NULL , **observation_histo = NULL;


  os << STAT_word[STATW_MIXTURE] << "\t" << nb_component << "\t" << STAT_word[STATW_DISTRIBUTIONS] << endl;

  // writing of the weights

  os << "\n" << STAT_word[STATW_WEIGHTS] << endl;

  for (i = 0;i < nb_component;i++) {
    os << weight->mass[i] << "\t";
  }
  os << endl;

  // writing of the distributions associated with each observation process

  os << "\n" << nb_output_process << "\t"
     << STAT_word[nb_output_process == 1 ? STATW_OUTPUT_PROCESS : STATW_OUTPUT_PROCESSES] << endl;

  distance = new double*[nb_component];
  for (i = 0;i < nb_component;i++) {
    distance[i] = new double[nb_component];
  }

  for (i = 0;i < nb_output_process;i++) {
    os << "\n" << STAT_word[STATW_OUTPUT_PROCESS] << "\t" << i + 1;

    if (categorical_process[i]) {
      os << "\t" << STAT_word[STATW_CATEGORICAL];
    }
    else if (discrete_parametric_process[i]) {
      os << "\t" << STAT_word[STATW_DISCRETE_PARAMETRIC];
    }
    else {
      os << "\t" << STAT_word[STATW_CONTINUOUS_PARAMETRIC];
    }
    os << endl;

    if (vec) {
      switch (vec->type[0]) {
      case STATE :
        variable = i + 1;
        break;
      default :
        variable = i;
        break;
      }

      if (vec->observation_distribution) {
        observation_dist = vec->observation_distribution[variable];
      }
      marginal_dist = vec->marginal_distribution[variable];

      if (vec->observation_histogram) {
        observation_histo = vec->observation_histogram[variable];
      }
      marginal_histo = vec->marginal_histogram[variable];

      if (categorical_process[i]) {
        categorical_process[i]->spreadsheet_print(os , observation_dist , marginal_dist , MIXTURE);

        for (j = 0;j < nb_component;j++) {
          for (k = j + 1;k < nb_component;k++) {
            distance[j][k] = categorical_process[i]->observation[j]->overlap_distance_computation(*(categorical_process[i]->observation[k]));
            distance[k][j] = distance[j][k];
          }
        }
      }

      else if (discrete_parametric_process[i]) {
        discrete_parametric_process[i]->spreadsheet_print(os , observation_dist , marginal_dist , MIXTURE);

        for (j = 0;j < nb_component;j++) {
          for (k = j + 1;k < nb_component;k++) {
            distance[j][k] = discrete_parametric_process[i]->observation[j]->sup_norm_distance_computation(*(discrete_parametric_process[i]->observation[k]));
            distance[k][j] = distance[j][k];
          }
        }
      }

      else {
        continuous_parametric_process[i]->spreadsheet_print(os , observation_histo , observation_dist ,
                                                            marginal_histo , marginal_dist , MIXTURE);

        for (j = 0;j < nb_component;j++) {
          for (k = j + 1;k < nb_component;k++) {
            distance[j][k] = continuous_parametric_process[i]->observation[j]->sup_norm_distance_computation(*(continuous_parametric_process[i]->observation[k]));
            distance[k][j] = distance[j][k];
          }
        }
      }

//      os << "\n" << STAT_label[STATL_COMPONENT_DISTANCE] << endl;
      os << "\n" << STAT_label[STATL_OBSERVATION_DISTRIBUTION_DISTANCE] << endl;

      for (j = 0;j < nb_component;j++) {
        for (k = 0;k < nb_component;k++) {
          if (k != j) {
            os << distance[j][k];
          }
          os << "\t";
        }
        os << endl;
      }
    }
  }

  for (i = 0;i < nb_component;i++) {
    delete [] distance[i];
  }
  delete [] distance;

  if (vec) {
    int nb_parameter = nb_parameter_computation(MIN_PROBABILITY);
    double information;


    os << "\n" << STAT_label[STATL_NB_VECTOR] << "\t" << vec->nb_vector << endl;

    // writing of the information quantity of the vectors in the iid case

    for (i = 0;i < vec->nb_variable;i++) {
      if (vec->type[i] == REAL_VALUE) {
        break;
      }
    }

    if (i == vec->nb_variable) {
      information = vec->information_computation();

      os << "\n" << STAT_label[STATL_INFORMATION] << "\t" << information << "\t"
         << information / vec->nb_vector << endl;
    }

    // writing of log-likelihoods for the observed vectors

    if (vec->restoration_likelihood != D_INF) {
      os << "\n" << STAT_label[STATL_CLASSIFICATION_LIKELIHOOD] << "\t" << vec->restoration_likelihood << "\t"
         << STAT_label[STATL_NORMALIZED] << "\t" << vec->restoration_likelihood / vec->nb_vector << endl;
    }

    if (vec->sample_entropy != D_DEFAULT) {
      os << "\n" << STAT_label[STATL_CLASSIFICATION_ENTROPY] << "\t" << vec->sample_entropy << "\t"
         << STAT_label[STATL_NORMALIZED] << "\t" << vec->sample_entropy / vec->nb_vector << endl;
    }

    if (vec->likelihood != D_INF) {
      os << "\n" << STAT_label[STATL_LIKELIHOOD] << "\t" << vec->likelihood << "\t"
         << STAT_label[STATL_NORMALIZED] << "\t" << vec->likelihood / vec->nb_vector << endl;
    }

    if (vec->likelihood != D_INF) {
      os << "\n" << nb_parameter << "\t" << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS] << "\t"
         << "2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[AIC] << ")\t"
         << 2 * (vec->likelihood - nb_parameter) << endl;

      if (nb_parameter < vec->nb_vector - 1) {
        os << "\n" << nb_parameter << "\t" << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS] << "\t"
           << "2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[AICc] << ")\t"
           << 2 * (vec->likelihood - (double)(nb_parameter * vec->nb_vector) /
              (double)(vec->nb_vector - nb_parameter - 1)) << endl;
      }

      os << "\n" << nb_parameter << "\t" << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS] << "\t"
         << "2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[BIC] << ")\t"
         << 2 * vec->likelihood - nb_parameter * log((double)vec->nb_vector) << endl;

      os << "\n" << nb_parameter << "\t" << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS] << "\t"
         << "2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[ICL] << ")\t"
//         << 2 * vec->restoration_likelihood - nb_parameter * log((double)vec->nb_vector) << endl;
         << 2 * (vec->likelihood - vec->sample_entropy) - nb_parameter * log((double)vec->nb_vector) << endl;
    }
  }

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a Mixture object in a file at the spreadsheet format.
 *
 *  \param[in] error reference on a StatError object,
 *  \param[in] path  file path.
 *
 *  \return          error status.
 */
/*--------------------------------------------------------------*/

bool Mixture::spreadsheet_write(StatError &error , const string path) const

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
    spreadsheet_write(out_file , mixture_data);
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Plot of a Mixture object and the associated data structure using Gnuplot.
 *
 *  \param[in] prefix file prefix,
 *  \param[in] title  figure title,
 *  \param[in] vec    pointer on a MixtureData object.
 *
 *  \return           error status.
 */
/*--------------------------------------------------------------*/

bool Mixture::plot_write(const char *prefix , const char *title ,
                         const MixtureData *vec) const

{
  bool status;
  int i;
  int variable , nb_value = I_DEFAULT;
  double *empirical_cdf[2];
  FrequencyDistribution *marginal_dist = NULL , **observation_dist = NULL;
  Histogram *marginal_histo = NULL , **observation_histo = NULL;
  ostringstream data_file_name;


  // writing of the data file

  data_file_name << prefix << 0 << ".dat";
  status = weight->plot_print((data_file_name.str()).c_str() , (vec ? vec->marginal_distribution[0] : NULL));

  if (status) {

    // writing of the script files

    for (i = 0;i < 2;i++) {
      ostringstream file_name[2];

      switch (i) {
      case 0 :
        file_name[0] << prefix << 0 << ".plot";
        break;
      case 1 :
        file_name[0] << prefix << 0 << ".print";
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

      if (weight->nb_value - 1 < TIC_THRESHOLD) {
        out_file << "set xtics 0,1" << endl;
      }

      if (vec) {
        out_file << "plot [0:" << weight->nb_value - 1 << "] [0:"
                 << (int)(MAX(vec->marginal_distribution[0]->max ,
                              weight->max * vec->marginal_distribution[0]->nb_element) * YSCALE) + 1
                 << "] \"" << label((data_file_name.str()).c_str()) << "\" using 1"
                 << " title \"" << STAT_label[STATL_WEIGHT] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
                 << "\" with impulses,\\" << endl;
        out_file << "\"" << label((data_file_name.str()).c_str()) << "\" using 2"
                 << " title \"" << STAT_label[STATL_WEIGHT] << " " << STAT_label[STATL_DISTRIBUTION]
                 << "\" with linespoints" << endl;
      }

      else {
        out_file << "plot [0:" << weight->nb_value - 1 << "] [0:" << weight->max * YSCALE
                 << "] \"" << label((data_file_name.str()).c_str()) << "\" using 1"
                 << " title \"" << STAT_label[STATL_WEIGHT] << " " << STAT_label[STATL_DISTRIBUTION]
                 << "\" with linespoints" << endl;
      }

      if (weight->nb_value - 1 < TIC_THRESHOLD) {
        out_file << "set xtics autofreq" << endl;
      }

      if (i == 1) {
        out_file << "\nset terminal x11" << endl;
      }

      out_file << "\npause 0 \"" << STAT_label[STATL_END] << "\"" << endl;
    }

    for (i = 0;i < nb_output_process;i++) {
      if (vec) {
        switch (vec->type[0]) {
        case STATE :
          variable = i + 1;
          break;
        default :
          variable = i;
          break;
        }

        if (vec->observation_distribution) {
          observation_dist = vec->observation_distribution[variable];
        }
        marginal_dist = vec->marginal_distribution[variable];

        if (vec->observation_histogram) {
          observation_histo = vec->observation_histogram[variable];
        }
        marginal_histo = vec->marginal_histogram[variable];

        if (continuous_parametric_process[i]) {
          nb_value = vec->cumulative_distribution_function_computation(variable , empirical_cdf);
        }
      }

      if (categorical_process[i]) {
        categorical_process[i]->plot_print(prefix , title , i + 1 , observation_dist ,
                                           marginal_dist , MIXTURE);
      }
      else if (discrete_parametric_process[i]) {
        discrete_parametric_process[i]->plot_print(prefix , title , i + 1 , observation_dist ,
                                                   marginal_dist , MIXTURE);
      }
      else if (continuous_parametric_process[i]) {
        continuous_parametric_process[i]->plot_print(prefix , title , i + 1 ,
                                                     observation_histo , observation_dist ,
                                                     marginal_histo , marginal_dist , nb_value ,
                                                     (vec ? empirical_cdf : NULL) , MIXTURE);
        if (vec) {
          delete [] empirical_cdf[0];
          delete [] empirical_cdf[1];
        }
      }
    }
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Plot of a Mixture object using Gnuplot.
 *
 *  \param[in] error  reference on a StatError object,
 *  \param[in] prefix file prefix,
 *  \param[in] title  figure title.
 *
 *  \return           error status.
 */
/*--------------------------------------------------------------*/

bool Mixture::plot_write(StatError &error , const char *prefix ,
                         const char *title) const

{
  bool status = plot_write(prefix , title , mixture_data);

  error.init();

  if (!status) {
    error.update(STAT_error[STATR_FILE_PREFIX]);
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Plot of a Mixture object and the associated data structure.
 *
 *  \param[in] vec pointer on a MixtureData object.
 *
 *  \return        MultiPlotSet object.
 */
/*--------------------------------------------------------------*/

MultiPlotSet* Mixture::get_plotable(const MixtureData *vec) const

{
  int i , j;
  int nb_plot_set , index , variable;
  FrequencyDistribution *marginal_dist = NULL , **observation_dist = NULL;
  Histogram *marginal_histo = NULL , **observation_histo = NULL;
  MultiPlotSet *plot_set;
  ostringstream title , legend;


  // computation of the number of plots

  nb_plot_set = 1;

  for (i = 0;i < nb_output_process;i++) {
    if (vec) {
      switch (vec->type[0]) {
      case STATE :
        variable = i + 1;
        break;
      default :
        variable = i;
        break;
      }
    }

    if ((vec->observation_distribution) || (vec->observation_histogram)) {
      nb_plot_set += nb_component;
    }
    else {
      nb_plot_set++;
    }

    if ((categorical_process[i]) && (vec->marginal_distribution[variable])) {
      if ((categorical_process[i]->weight) &&
          (categorical_process[i]->mixture)) {
        nb_plot_set++;
      }
    }

    if ((discrete_parametric_process[i]) && (vec->marginal_distribution[variable])) {
      if ((discrete_parametric_process[i]->weight) &&
          (discrete_parametric_process[i]->mixture)) {
        nb_plot_set += 2;
      }
    }

    if ((continuous_parametric_process[i]) && ((vec->marginal_histogram[variable]) ||
         (vec->marginal_distribution[variable]))) {
      if (continuous_parametric_process[i]->weight) {
        nb_plot_set += 2;
      }
    }
  }

  plot_set = new MultiPlotSet(nb_plot_set , nb_output_process + 1);
  MultiPlotSet &plot = *plot_set;

  plot_set->border = "15 lw 0";

  index = 0;
  plot_set->variable_nb_viewpoint[0] = 0;

  plot[0].xrange = Range(0 , weight->nb_value - 1);

  if (vec) {

    // weights

    plot[0].yrange = Range(0 , ceil(MAX(vec->marginal_distribution[0]->max ,
                                    weight->max * vec->marginal_distribution[0]->nb_element)
                                    * YSCALE));

    if (weight->nb_value - 1 < TIC_THRESHOLD) {
      plot[0].xtics = 1;
    }

    plot[0].resize(2);

    legend.str("");
    legend << STAT_label[STATL_WEIGHT] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
    plot[0][0].legend = legend.str();

    plot[0][0].style = "impulses";

    vec->marginal_distribution[0]->plotable_frequency_write(plot[0][0]);

    legend.str("");
    legend << STAT_label[STATL_WEIGHT] << " " << STAT_label[STATL_DISTRIBUTION];
    plot[0][1].legend = legend.str();

    plot[0][1].style = "linespoints";

    weight->plotable_mass_write(plot[0][1] , vec->marginal_distribution[0]->nb_element);
  }

  else {
    plot[0].yrange = Range(0 , weight->max * YSCALE);
    legend.str("");
    legend << STAT_label[STATL_WEIGHT] << " " << STAT_label[STATL_DISTRIBUTION];
    plot[0][0].legend = legend.str();

    plot[0][0].style = "linespoints";

    weight->plotable_mass_write(plot[0][0] , vec->marginal_distribution[0]->nb_element);
  }

  for (i = 0;i < nb_output_process;i++) {
    if (vec) {
      switch (vec->type[0]) {
      case STATE :
        variable = i + 1;
        break;
      default :
        variable = i;
        break;
      }

      if (vec->observation_distribution) {
        observation_dist = vec->observation_distribution[variable];
      }
      marginal_dist = vec->marginal_distribution[variable];

      if (vec->observation_histogram) {
        observation_histo = vec->observation_histogram[variable];
      }
      marginal_histo = vec->marginal_histogram[variable];
    }

    if (categorical_process[i]) {
      plot_set->variable_nb_viewpoint[i] = 0;
      categorical_process[i]->plotable_write(*plot_set , index , i + 1 , observation_dist ,
                                             marginal_dist , MIXTURE);
    }
    else if (discrete_parametric_process[i]) {
      discrete_parametric_process[i]->plotable_write(*plot_set , index , i + 1 , observation_dist ,
                                                     marginal_dist , MIXTURE);
    }
    else {
      continuous_parametric_process[i]->plotable_write(*plot_set , index , i + 1 ,
                                                       observation_histo , observation_dist ,
                                                       marginal_histo , marginal_dist , MIXTURE);
    }
  }

  return plot_set;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Plot of a Mixture object.
 *
 *  \return MultiPlotSet object.
 */
/*--------------------------------------------------------------*/

MultiPlotSet* Mixture::get_plotable() const

{
  return get_plotable(mixture_data);
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the number of parameters of a Mixture object.
 *
 *  \param[in] min_probability minimum probability.
 *
 *  \return                    number of parameters.
 */
/*--------------------------------------------------------------*/

int Mixture::nb_parameter_computation(double min_probability) const

{
  int i;
  int nb_parameter = nb_component - 1;


  for (i = 0;i < nb_output_process;i++) {
    if (categorical_process[i]) {
      nb_parameter += categorical_process[i]->nb_parameter_computation(min_probability);
    }
    else if (discrete_parametric_process[i]) {
      nb_parameter += discrete_parametric_process[i]->nb_parameter_computation();
    }
    else {
      nb_parameter += continuous_parametric_process[i]->nb_parameter_computation();
    }
  }

  return nb_parameter;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Default constructor of the MixtureData class.
 */
/*--------------------------------------------------------------*/

MixtureData::MixtureData()

{
  mixture = NULL;

  observation_distribution = NULL;
  observation_histogram = NULL;

  likelihood = D_INF;
  restoration_likelihood = D_INF;
  sample_entropy = D_DEFAULT;

  posterior_probability = NULL;
  entropy = NULL;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Constructor of the MixtureData class.
 *
 *  \param[in] inb_vector   number of individuals,
 *  \param[in] inb_variable number of variables,
 *  \param[in] itype        variable types,
 *  \param[in] init_flag    flag initialization.
 */
/*--------------------------------------------------------------*/

MixtureData::MixtureData(int inb_vector , int inb_variable ,
                         variable_nature *itype , bool init_flag)
:Vectors(inb_vector , NULL , inb_variable , itype , init_flag)

{
  mixture = NULL;

  observation_distribution = NULL;
  observation_histogram = NULL;

  likelihood = D_INF;
  restoration_likelihood = D_INF;
  sample_entropy = D_DEFAULT;

  posterior_probability = NULL;
  entropy = NULL;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Construction of a MixtureData object from a Vectors object.
 *
 *  \param[in] vec       reference on a Vectors object,
 *  \param[in] transform type of transform(VECTOR_COPY/ADD_COMPONENT_VARIABLE).
 */
/*--------------------------------------------------------------*/

MixtureData::MixtureData(const Vectors &vec , vector_transformation transform)
:Vectors(vec , transform)

{
  mixture = NULL;

  observation_distribution = NULL;
  observation_histogram = NULL;

  likelihood = D_INF;
  restoration_likelihood = D_INF;
  sample_entropy = D_DEFAULT;

  posterior_probability = NULL;
  entropy = NULL;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Copy of a MixtureData object.
 *
 *  \param[in] vec        reference on a MixtureData object,
 *  \param[in] model_flag flag copy of the Mixture object.
 */
/*--------------------------------------------------------------*/

void MixtureData::copy(const MixtureData &vec , bool model_flag)

{
  int i , j;


  if ((model_flag) && (vec.mixture)) {
    mixture = new Mixture(*(vec.mixture) , false);
  }
  else {
    mixture = NULL;
  }

  if (vec.observation_distribution) {
    observation_distribution = new FrequencyDistribution**[nb_variable];
    observation_distribution[0] = NULL;

    for (i = 1;i < nb_variable;i++) {
      if (vec.observation_distribution[i]) {
        observation_distribution[i] = new FrequencyDistribution*[marginal_distribution[0]->nb_value];
        for (j = 0;j < marginal_distribution[0]->nb_value;j++) {
          observation_distribution[i][j] = new FrequencyDistribution(*(vec.observation_distribution[i][j]));
        }
      }

      else {
        observation_distribution[i] = NULL;
      }
    }
  }

  else {
    observation_distribution = NULL;
  }

  if (vec.observation_histogram) {
    observation_histogram = new Histogram**[nb_variable];
    observation_histogram[0] = NULL;

    for (i = 1;i < nb_variable;i++) {
      if (vec.observation_histogram[i]) {
        observation_histogram[i] = new Histogram*[marginal_distribution[0]->nb_value];
        for (j = 0;j < marginal_distribution[0]->nb_value;j++) {
          observation_histogram[i][j] = new Histogram(*(vec.observation_histogram[i][j]));
        }
      }

      else {
        observation_histogram[i] = NULL;
      }
    }
  }

  else {
    observation_histogram = NULL;
  }

  likelihood = vec.likelihood;
  restoration_likelihood = vec.restoration_likelihood;
  sample_entropy = vec.sample_entropy;

  if (vec.posterior_probability) {
    posterior_probability = new double[nb_vector];
    for (i = 0;i < nb_vector;i++) {
      posterior_probability[i] = vec.posterior_probability[i];
    }
  }
  else {
    posterior_probability = NULL;
  }

  if (vec.entropy) {
    entropy = new double[nb_vector];
    for (i = 0;i < nb_vector;i++) {
      entropy[i] = vec.entropy[i];
    }
  }
  else {
    entropy = NULL;
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Destruction of the data members of the MixtureData class.
 */
/*--------------------------------------------------------------*/

void MixtureData::remove()

{
  int i , j;


  delete mixture;

  if (observation_distribution) {
    for (i = 1;i < nb_variable;i++) {
      if (observation_distribution[i]) {
        for (j = 0;j < marginal_distribution[0]->nb_value;j++) {
          delete observation_distribution[i][j];
        }
        delete [] observation_distribution[i];
      }
    }
    delete [] observation_distribution;
  }

  if (observation_histogram) {
    for (i = 1;i < nb_variable;i++) {
      if (observation_histogram[i]) {
        for (j = 0;j < marginal_distribution[0]->nb_value;j++) {
          delete observation_histogram[i][j];
        }
        delete [] observation_histogram[i];
      }
    }
    delete [] observation_histogram;
  }

  delete [] posterior_probability;
  delete [] entropy;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Destructor of the MixtureData class.
 */
/*--------------------------------------------------------------*/

MixtureData::~MixtureData()

{
  remove();
}


/*--------------------------------------------------------------*/
/**
 *  \brief Assignment operator of the MixtureData class.
 *
 *  \param[in] vec  reference on a MixtureData object.
 *
 *  \return         MixtureData object.
 */
/*--------------------------------------------------------------*/

MixtureData& MixtureData::operator=(const MixtureData &vec)

{
  if (&vec != this) {
    remove();
    Vectors::remove();

    Vectors::copy(vec);
    copy(vec);
  }

  return *this;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Initialization of the 1st variable.
 *
 *  \param[in] itype 1st variable type (STATE/INT_VALUE/REAL_VALUE).
 */
/*--------------------------------------------------------------*/

void MixtureData::state_variable_init(variable_nature itype)

{
  int i , j;


  if (itype != type[0]) {
    if (type[0] == STATE) {
      if (observation_distribution) {
        for (i = 1;i < nb_variable;i++) {
          if (observation_distribution[i]) {
            for (j = 0;j < marginal_distribution[0]->nb_value;j++) {
              delete observation_distribution[i][j];
            }
            delete [] observation_distribution[i];
          }
        }
        delete [] observation_distribution;

        observation_distribution = NULL;
      }

      if (observation_histogram) {
        for (i = 1;i < nb_variable;i++) {
          if (observation_histogram[i]) {
            for (j = 0;j < marginal_distribution[0]->nb_value;j++) {
              delete observation_histogram[i][j];
            }
            delete [] observation_histogram[i];
          }
        }
        delete [] observation_histogram;

        observation_histogram = NULL;
      }
    }

    type[0] = itype;
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Extraction of an empirical component.
 *
 *  \param[in] error    reference on a StatError object,
 *  \param[in] variable variable index,
 *  \param[in] index    component index.
 *
 *  \return             DiscreteDistributionData object.
 */
/*--------------------------------------------------------------*/

DiscreteDistributionData* MixtureData::extract(StatError &error , int variable , int index) const

{
  bool status = true;
  Distribution *pdist;
  DiscreteParametric *pparam;
  FrequencyDistribution *phisto;
  DiscreteDistributionData *histo;


  histo = NULL;
  error.init();

  if ((variable < 2) || (variable > nb_variable)) {
    status = false;
    error.update(STAT_error[STATR_VARIABLE_INDEX]);
  }

  else {
    variable--;

    if ((index < 1) || (index > marginal_distribution[0]->nb_value)) {
      status = false;
      error.update(STAT_error[STATR_FREQUENCY_DISTRIBUTION_INDEX]);
    }

    else {
      index--;

      if (!observation_distribution[variable]) {
        status = false;
        ostringstream error_message;
        error_message << STAT_label[STATL_VARIABLE] << " " << variable + 1 << ": "
                      << STAT_error[STATR_VARIABLE_TYPE];
        error.correction_update((error_message.str()).c_str() , STAT_variable_word[INT_VALUE]);
      }

      else {
        phisto = observation_distribution[variable][index];

        if (phisto->nb_element == 0) {
          status = false;
          error.update(STAT_error[STATR_EMPTY_SAMPLE]);
        }
      }
    }
  }

  if (status) {
    pdist = NULL;
    pparam = NULL;

    if (mixture->categorical_process[variable - 1]) {
      pdist = mixture->categorical_process[variable - 1]->observation[index];
    }
    else if (mixture->discrete_parametric_process[variable - 1]) {
      pparam = mixture->discrete_parametric_process[variable - 1]->observation[index];
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
 *  \brief Writing of a MixtureData object.
 *
 *  \param[in,out] os         stream,
 *  \param[in]     exhaustive flag detail level.
 */
/*--------------------------------------------------------------*/

ostream& MixtureData::ascii_write(ostream &os , bool exhaustive) const

{
  if (mixture) {
    mixture->ascii_write(os , this , exhaustive , false);
  }

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a MixtureData object in a file.
 *
 *  \param[in] error       reference on a StatError object,
 *  \param[in] path        file path,
 *  \param[in] exhaustive  flag detail level.
 *
 *  \return                error status.
 */
/*--------------------------------------------------------------*/

bool MixtureData::ascii_write(StatError &error , const string path , bool exhaustive) const

{
  bool status = false;


  if (mixture) {
    ofstream out_file(path.c_str());

    error.init();

    if (!out_file) {
      status = false;
      error.update(STAT_error[STATR_FILE_NAME]);
    }

    else {
      status = true;
      mixture->ascii_write(out_file , this , exhaustive);
    }
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a MixtureData object.
 *
 *  \param[in,out] os         stream,
 *  \param[in]     exhaustive flag detail level.
 */
/*--------------------------------------------------------------*/

ostream& MixtureData::ascii_data_write(ostream &os , bool exhaustive) const

{
  Vectors::ascii_write(os , exhaustive , false);
  ascii_print(os , false , posterior_probability , entropy);

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a MixtureData object in a file.
 *
 *  \param[in] error      reference on a StatError object,
 *  \param[in] path       file path,
 *  \param[in] exhaustive flag detail level.
 *
 *  \return               error status.
 */
/*--------------------------------------------------------------*/

bool MixtureData::ascii_data_write(StatError &error , const string path ,
                                   bool exhaustive) const

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
    Vectors::ascii_write(out_file , exhaustive , true);
    ascii_print(out_file , true , posterior_probability , entropy);
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a MixtureData object in a file at the spreadsheet format.
 *
 *  \param[in] error reference on a StatError object,
 *  \param[in] path  file path.
 *
 *  \return          error status.
 */
/*--------------------------------------------------------------*/

bool MixtureData::spreadsheet_write(StatError &error , const string path) const

{
  bool status = false;


  if (mixture) {
    ofstream out_file(path.c_str());

    error.init();

    if (!out_file) {
      status = false;
      error.update(STAT_error[STATR_FILE_NAME]);
    }

    else {
      status = true;
      mixture->spreadsheet_write(out_file , this);
    }
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Plot of a MixtureData object using Gnuplot.
 *
 *  \param[in] error  reference on a StatError object,
 *  \param[in] prefix file prefix,
 *  \param[in] title  figure title.
 *
 *  \return           error status.
 */
/*--------------------------------------------------------------*/

bool MixtureData::plot_write(StatError &error , const char *prefix ,
                             const char *title) const

{
  bool status = false;


  if (mixture) {
    status = mixture->plot_write(prefix , title , this);

    error.init();

    if (!status) {
      error.update(STAT_error[STATR_FILE_PREFIX]);
    }
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Plot of a MixtureData object.
 *
 *  \return MultiPlotSet object.
 */
/*--------------------------------------------------------------*/

MultiPlotSet* MixtureData::get_plotable() const

{
  MultiPlotSet *plot_set;


  if (mixture) {
    plot_set = mixture->get_plotable(this);
  }
  else {
    plot_set = NULL;
  }

  return plot_set;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the information quantity.
 *
 *  \return information quantity.
 */
/*--------------------------------------------------------------*/

double MixtureData::information_computation() const

{
  int i;
  double information = 0.;


  for (i = (((type[0] != STATE) || (nb_variable == 1)) ? 0 : 1);i < nb_variable;i++) {
    if (marginal_distribution[i]) {
      information += marginal_distribution[i]->information_computation();
    }
    else {
      information = D_INF;
      break;
    }
  }

  return information;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Update of the observation frequency distributions for a variable.
 *
 *  \param[in] variable     variable index,
 *  \param[in] nb_component number of components.
 */
/*--------------------------------------------------------------*/

void MixtureData::observation_frequency_distribution_computation(int variable , int nb_component)

{
  int i , j;


  // initialization of the observation frequency distributions

  for (i = 0;i < nb_component;i++) {
    for (j = 0;j < marginal_distribution[variable]->nb_value;j++) {
      observation_distribution[variable][i]->frequency[j] = 0;
    }
  }

  // update of the observation frequency distributions

  for (i = 0;i < nb_vector;i++) {
    (observation_distribution[variable][int_vector[i][0]]->frequency[int_vector[i][variable]])++;
  }

  // extraction of the characteristics of the observation frequency distributions

  for (i = 0;i < nb_component;i++) {
    observation_distribution[variable][i]->nb_value_computation();
    observation_distribution[variable][i]->offset_computation();
    observation_distribution[variable][i]->nb_element_computation();
    observation_distribution[variable][i]->max_computation();
    observation_distribution[variable][i]->mean_computation();
    observation_distribution[variable][i]->variance_computation();
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Construction of the observation frequency distributions.
 *
 *  \param[in] nb_component number of components.
 */
/*--------------------------------------------------------------*/

void MixtureData::build_observation_frequency_distribution(int nb_component)

{
  if ((nb_variable > 1) && (!observation_distribution)) {
    int i , j;


    observation_distribution = new FrequencyDistribution**[nb_variable];
    observation_distribution[0] = NULL;

    for (i = 1;i < nb_variable;i++) {
      if (marginal_distribution[i]) {
        observation_distribution[i] = new FrequencyDistribution*[nb_component];
        for (j = 0;j < nb_component;j++) {
          observation_distribution[i][j] = new FrequencyDistribution(marginal_distribution[i]->nb_value);
        }

        observation_frequency_distribution_computation(i , nb_component);
      }

      else {
        observation_distribution[i] = NULL;
      }
    }
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Construction of the observation histograms for a variable.
 *
 *  \param[in] variable     variable index,
 *  \param[in] nb_component number of components,
 *  \param[in] bin_width    bin width.
 */
/*--------------------------------------------------------------*/

void MixtureData::build_observation_histogram(int variable , int nb_component , double bin_width)

{
  if ((!observation_histogram[variable]) || (bin_width != observation_histogram[variable][0]->bin_width)) {
    int i , j;
    double imin_value;


    // construction of the observation histograms

    if (bin_width == D_DEFAULT) {
      bin_width = marginal_histogram[variable]->bin_width;
    }
    imin_value = floor(min_value[variable] / bin_width) * bin_width;

    if (observation_histogram[variable]) {
      for (i = 0;i < nb_component;i++) {
        observation_histogram[variable][i]->nb_bin = (int)floor((max_value[variable] - imin_value) / bin_width) + 1;

        delete [] observation_histogram[variable][i]->frequency;
        observation_histogram[variable][i]->frequency = new int[observation_histogram[variable][i]->nb_bin];
      }
    }

    else {
      observation_histogram[variable] = new Histogram*[nb_component];

      for (i = 0;i < nb_component;i++) {
        observation_histogram[variable][i] = new Histogram((int)floor((max_value[variable] - imin_value) / bin_width) + 1 , false);

        observation_histogram[variable][i]->nb_element = marginal_distribution[0]->frequency[i];
        observation_histogram[variable][i]->type = type[variable];
      }

      // computation of the minimum and maximum values for each component

/*      for (i = 0;i < nb_component;i++) {
        observation_histogram[variable][i]->min_value = max_value[variable];
        observation_histogram[variable][i]->max_value = min_value[variable];
      }

      switch (type[variable]) {

      case INT_VALUE : {
        for (i = 0;i < nb_vector;i++) {
          if (int_vector[i][variable] < observation_histogram[variable][int_vector[i][0]]->min_value) {
            observation_histogram[variable][int_vector[i][0]]->min_value = int_vector[i][variable];
          }
          if (int_vector[i][variable] > observation_histogram[variable][int_vector[i][0]]->max_value) {
            observation_histogram[variable][int_vector[i][0]]->max_value = int_vector[i][variable];
          }
        }
        break;
      }

      case REAL_VALUE : {
        for (i = 0;i < nb_vector;i++) {
          if (real_vector[i][variable] < observation_histogram[variable][int_vector[i][0]]->min_value) {
            observation_histogram[variable][int_vector[i][0]]->min_value = real_vector[i][variable];
          }
          if (real_vector[i][variable] > observation_histogram[variable][int_vector[i][0]]->max_value) {
            observation_histogram[variable][int_vector[i][0]]->max_value = real_vector[i][variable];
          }
        }
        break;
      }
      } */
    }

    for (i = 0;i < nb_component;i++) {
      observation_histogram[variable][i]->bin_width = bin_width;
      observation_histogram[variable][i]->min_value = imin_value;
      observation_histogram[variable][i]->max_value = ceil(max_value[variable] / bin_width) * bin_width;
    }

    // computation of frequencies

    for (i = 0;i < nb_component;i++) {
      for (j = 0;j < observation_histogram[variable][i]->nb_bin;j++) {
        observation_histogram[variable][i]->frequency[j] = 0;
      }
    }

    switch (type[variable]) {

    case INT_VALUE : {
      for (i = 0;i < nb_vector;i++) {
//        (observation_histogram[variable][int_vector[i][0]]->frequency[(int)((int_vector[i][variable] - imin_value) / bin_width)])++;
        (observation_histogram[variable][int_vector[i][0]]->frequency[(int)floor((int_vector[i][variable] - imin_value) / bin_width)])++;
      }
      break;
    }

    case REAL_VALUE : {
      for (i = 0;i < nb_vector;i++) {
//        (observation_histogram[variable][int_vector[i][0]]->frequency[(int)((real_vector[i][variable] - imin_value) / bin_width)]++;
        (observation_histogram[variable][int_vector[i][0]]->frequency[(int)floor((real_vector[i][variable] - imin_value) / bin_width)])++;
      }
      break;
    }
    }

    for (i = 0;i < nb_component;i++) {
      observation_histogram[variable][i]->max_computation();
    }
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Construction of the observation histograms.
 *
 *  \param[in] nb_component number of components.
 */
/*--------------------------------------------------------------*/

void MixtureData::build_observation_histogram(int nb_component)

{
  if ((nb_variable > 1) && (!observation_histogram)) {
    int i;


    observation_histogram = new Histogram**[nb_variable];
    observation_histogram[0] = NULL;

    for (i = 1;i < nb_variable;i++) {
      observation_histogram[i] = NULL;
      if (marginal_histogram[i]) {
        build_observation_histogram(i , nb_component , marginal_histogram[i]->bin_width);
      }
    }
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Change of the bin width of the marginal histogram and
 *         the observation histograms for a variable.
 *
 *  \param[in] error      reference on a StatError object,
 *  \param[in] variable   variable index,
 *  \param[in] bin_width  bin width,
 *  \param[in] imin_value minimum value.
 *
 *  \return               error status.
 */
/*--------------------------------------------------------------*/

bool MixtureData::select_bin_width(StatError &error , int variable ,
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

    if ((observation_histogram) && (observation_histogram[variable])) {
      build_observation_histogram(variable , marginal_distribution[0]->nb_value , bin_width);
    }
  }

  return status;
}


};  // namespace stat_tool
