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

#include <sstream>
#include <iomanip>

#include "tool/rw_tokenizer.h"
#include "tool/rw_cstring.h"
#include "tool/rw_locale.h"
#include "tool/config.h"

#include "stat_tools.h"
#include "distribution.h"
#include "markovian.h"
#include "vectors.h"
#include "mixture.h"
#include "stat_label.h"

using namespace std;


namespace stat_tool {



/*--------------------------------------------------------------*
 *
 *  Constructeur par defaut de la classe Mixture.
 *
 *--------------------------------------------------------------*/

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


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Mixture.
 *
 *  arguments : nombre de composantes, nombre de processus d'observation,
 *              nombre de valeurs observees par processus.
 *
 *--------------------------------------------------------------*/

Mixture::Mixture(int inb_component , int inb_output_process , int *nb_value)

{
  register int i;


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


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Mixture (modele gamma ou gaussien univarié avec parametres lies).
 *
 *  arguments : nombre de composantes, identificateur des composantes (GAMMA / GAUSSIAN),
 *              moyenne et parametre de forme de la premiere composante (GAMMA) /
 *              moyenne et ecart-type de la premiere composante (GAUSSIAN), flag moyennes liees,
 *              type de lien entre les variances (CONVOLUTION_FACTOR / SCALING_FACTOR).
 *
 *--------------------------------------------------------------*/

Mixture::Mixture(int inb_component , int ident , double mean , double standard_deviation ,
                 bool tied_mean , int variance_factor)

{
  register int i , j;
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


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Mixture.
 *
 *  arguments : pointeur sur objet DiscreteParametric,
 *              nombre de processus d'observation, pointeurs sur des objets
 *              CategoricalProcess, DiscreteParametricProcess et
 *              ContinuousParametricProcess.
 *
 *--------------------------------------------------------------*/

Mixture::Mixture(const DiscreteParametric *iweight , int inb_output_process ,
                 CategoricalProcess **categorical_observation ,
                 DiscreteParametricProcess **discrete_parametric_observation ,
                 ContinuousParametricProcess **continuous_parametric_observation)

{
  register int i;


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


/*--------------------------------------------------------------*
 *
 *  Copie d'un objet Mixture.
 *
 *  arguments : reference sur un objet Mixture,
 *              flag copie de l'objet MixtureData.
 *
 *--------------------------------------------------------------*/

void Mixture::copy(const Mixture &mixt , bool data_flag)

{
  register int i;


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


/*--------------------------------------------------------------*
 *
 *  Destruction des champs d'un objet Mixture.
 *
 *--------------------------------------------------------------*/

void Mixture::remove()

{
  register int i;


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


/*--------------------------------------------------------------*
 *
 *  Destructeur de la classe Mixture.
 *
 *--------------------------------------------------------------*/

Mixture::~Mixture()

{
  remove();
}


/*--------------------------------------------------------------*
 *
 *  Operateur d'assignement de la classe Mixture.
 *
 *  argument : reference sur un objet Mixture.
 *
 *--------------------------------------------------------------*/

Mixture& Mixture::operator=(const Mixture &mixt)

{
  if (&mixt != this) {
    remove();
    copy(mixt);
  }

  return *this;
}


/*--------------------------------------------------------------*
 *
 *  Extraction d'une loi.
 *
 *  arguments : reference sur un objet StatError, variable, indice de la composante.
 *
 *--------------------------------------------------------------*/

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


/*--------------------------------------------------------------*
 *
 *  Extraction de la partie "donnees" d'un objet Mixture.
 *
 *  argument : reference sur un objet StatError.
 *
 *--------------------------------------------------------------*/

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


/*--------------------------------------------------------------*
 *
 *  Application d'un seuil sur les parametres d'un melange multivarie.
 *
 *  argument : probabilite minimum.
 *
 *--------------------------------------------------------------*/

Mixture* Mixture::thresholding(double min_probability) const

{
  register int i;
  Mixture *mixt;


  mixt = new Mixture(*this , false);

  for (i = 0;i < mixt->nb_output_process;i++) {
    if (mixt->categorical_process[i]) {
      mixt->categorical_process[i]->thresholding(min_probability);
    }
  }

  return mixt;
}


/*--------------------------------------------------------------*
 *
 *  Construction d'un objet Mixture a partir d'un fichier.
 *
 *  arguments : reference sur un objet StatError, path,
 *              seuil sur les fonctions de repartition des lois parametriques.
 *
 *--------------------------------------------------------------*/

Mixture* mixture_ascii_read(StatError &error , const char *path , double cumul_threshold)

{
  RWLocaleSnapshot locale("en");
  RWCString buffer , token;
  size_t position;
  bool status , lstatus;
  register int i;
  int line , read_line , nb_output_process , output_process_type , index;
  long nb_component , value;
  double proba , cumul;
  DiscreteParametric *weight;
  CategoricalProcess **categorical_observation;
  DiscreteParametricProcess **discrete_parametric_observation;
  ContinuousParametricProcess **continuous_parametric_observation;
  Mixture *mixt;
  ifstream in_file(path);


  mixt = NULL;
  error.init();

  if (!in_file) {
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {
    status = true;
    line = 0;
    nb_component = 0;

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

        // test mot cle MIXTURE

        case 0 : {
          if (token != STAT_word[STATW_MIXTURE]) {
            status = false;
            error.correction_update(STAT_parsing[STATP_KEY_WORD] , STAT_word[STATW_MIXTURE] , line , i + 1);
          }
          break;
        }

        // test nombre de composantes

        case 1 : {
          lstatus = locale.stringToNum(token , &nb_component);
          if ((lstatus) && ((nb_component < 2) || (nb_component > MIXTURE_NB_COMPONENT))) {
            lstatus = false;
          }

          if (!lstatus) {
            status = false;
            error.update(STAT_parsing[STATP_NB_COMPONENT] , line , i + 1);
          }
          break;
        }

        // test mot cle COMPONENTS

        case 2 : {
          if (token != STAT_word[STATW_COMPONENTS]) {
            status = false;
            error.correction_update(STAT_parsing[STATP_KEY_WORD] , STAT_word[STATW_COMPONENTS] , line , i + 1);
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

      // analyse du format et lecture de la loi des poids

      read_line = 0;
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

        if (read_line == 0) {
          while (!((token = next()).isNull())) {

            // test mot cle WEIGHTS

            if (i == 0) {
              if (token != STAT_word[STATW_WEIGHTS]) {
                status = false;
                error.correction_update(STAT_parsing[STATP_KEY_WORD] , STAT_word[STATW_WEIGHTS] , line);
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

          while (!((token = next()).isNull())) {
            if (i < nb_component) {
              lstatus = locale.stringToNum(token , &proba);
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

      // analyse du format et lecture des lois d'observation

      if (status) {
        nb_output_process = I_DEFAULT;

        categorical_observation = NULL;
        discrete_parametric_observation = NULL;
        continuous_parametric_observation = NULL;

        while (buffer.readLine(in_file , false)) {
          line++;

#         ifdef DEBUG
          cout << line << "  " << buffer << endl;
#         endif

          position = buffer.first('#');
          if (position != RW_NPOS) {
            buffer.remove(position);
          }
          i = 0;

          RWCTokenizer next(buffer);

          while (!((token = next()).isNull())) {
            switch (i) {

            // test nombre de processus d'observation

            case 0 : {
              lstatus = locale.stringToNum(token , &value);
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

            // test mot cle OUTPUT_PROCESS(ES)

            case 1 : {
              if (token != STAT_word[nb_output_process == 1 ? STATW_OUTPUT_PROCESS : STATW_OUTPUT_PROCESSES]) {
                status = false;
                error.correction_update(STAT_parsing[STATP_KEY_WORD] ,
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

          while (buffer.readLine(in_file , false)) {
            line++;

#           ifdef DEBUG
            cout << line << "  " << buffer << endl;
#           endif

            position = buffer.first('#');
            if (position != RW_NPOS) {
              buffer.remove(position);
            }
            i = 0;

            RWCTokenizer next(buffer);

            while (!((token = next()).isNull())) {
              switch (i) {

              // test mot cle OUTPUT_PROCESS

              case 0 : {
                if (token == STAT_word[STATW_OUTPUT_PROCESS]) {
                  index++;
                }
                else {
                  status = false;
                  error.correction_update(STAT_parsing[STATP_KEY_WORD] , STAT_word[STATW_OUTPUT_PROCESS] , line , i + 1);
                }
                break;
              }

              // test indice du processus d'observation

              case 1 : {
                lstatus = locale.stringToNum(token , &value);
                if ((lstatus) && ((value != index) || (value > nb_output_process))) {
                  lstatus = false;
                }

                if (!lstatus) {
                  status = false;
                  error.update(STAT_parsing[STATP_OUTPUT_PROCESS_INDEX] , line , i + 1);
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

              // test mot cle CATEGORICAL / DISCRETE_PARAMETRIC / CONTINUOUS_PARAMETRIC

              case 3 : {
                if ((token == STAT_word[STATW_CATEGORICAL]) ||
                    (token == STAT_word[STATW_NONPARAMETRIC])) {
                  output_process_type = CATEGORICAL_PROCESS;
                }
                else if ((token == STAT_word[STATW_DISCRETE_PARAMETRIC]) ||
                         (token == STAT_word[STATW_PARAMETRIC])) {
                  output_process_type = DISCRETE_PARAMETRIC;
                }
                else if (token == STAT_word[STATW_CONTINUOUS_PARAMETRIC]) {
                  output_process_type = CONTINUOUS_PARAMETRIC;
                }
                else {
                  output_process_type = CATEGORICAL_PROCESS - 1;
                  status = false;
                  ostringstream correction_message;
                  correction_message << STAT_word[STATW_CATEGORICAL] << " or "
                                     << STAT_word[STATW_DISCRETE_PARAMETRIC] << " or "
                                     << STAT_word[STATW_CONTINUOUS_PARAMETRIC];
                  error.correction_update(STAT_parsing[STATP_KEY_WORD] , (correction_message.str()).c_str() , line , i + 1);
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

              switch (output_process_type) {

              case CATEGORICAL_PROCESS : {
                categorical_observation[index - 1] = categorical_observation_parsing(error , in_file , line ,
                                                                                     nb_component , MIXTURE , true);
                if (!categorical_observation[index - 1]) {
                  status = false;
                }
                break;
              }

              case DISCRETE_PARAMETRIC : {
                discrete_parametric_observation[index - 1] = discrete_observation_parsing(error , in_file , line ,
                                                                                          nb_component , MIXTURE ,
                                                                                          cumul_threshold);
                if (!discrete_parametric_observation[index - 1]) {
                  status = false;
                }
                break;
              }

              case CONTINUOUS_PARAMETRIC : {
                continuous_parametric_observation[index - 1] = continuous_observation_parsing(error , in_file , line ,
                                                                                              nb_component , MIXTURE ,
                                                                                              VON_MISES);
                if (!continuous_parametric_observation[index - 1]) {
                  status = false;
                }
                break;
              }
              }
            }
          }

          if (index != nb_output_process) {
            status = false;
            error.update(STAT_parsing[STATP_FORMAT] , line);
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


/*--------------------------------------------------------------*
 *
 *  Ecriture sur une ligne d'un objet Mixture.
 *
 *  argument : stream.
 *
 *--------------------------------------------------------------*/

ostream& Mixture::line_write(ostream &os) const

{
  os << nb_component << " " << STAT_word[STATW_COMPONENTS];

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Mixture et de la structure de donnees associee.
 *
 *  arguments : stream, pointeur sur un objet MixtureData,
 *              flag niveau de detail, flag fichier.
 *
 *--------------------------------------------------------------*/

ostream& Mixture::ascii_write(ostream &os , const MixtureData *vec ,
                               bool exhaustive , bool file_flag) const

{
  register int i , j , k;
  int buff , width , variable;
  double mean , **distance;
  FrequencyDistribution *marginal_dist = NULL , **observation_dist = NULL;
  Histogram *marginal_histo = NULL , **observation_histo = NULL;
  long old_adjust;


  old_adjust = os.setf(ios::left , ios::adjustfield);

  os << STAT_word[STATW_MIXTURE] << " " << nb_component << " " << STAT_word[STATW_COMPONENTS] << endl;

  // ecriture des poids

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

  // ecriture des lois associees a chaque processus d'observation

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

    // ecriture de la quantite d'information des vecteurs

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

    // ecriture des vraisemblances des vecteurs

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

  os.setf((FMTFLAGS)old_adjust , ios::adjustfield);

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Mixture.
 *
 *  arguments : stream, flag niveau de detail.
 *
 *--------------------------------------------------------------*/

ostream& Mixture::ascii_write(ostream &os , bool exhaustive) const

{
  return ascii_write(os , mixture_data , exhaustive , false);
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Mixture dans un fichier.
 *
 *  arguments : reference sur un objet StatError, path,
 *              flag niveau de detail.
 *
 *--------------------------------------------------------------*/

bool Mixture::ascii_write(StatError &error , const char *path ,
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
    ascii_write(out_file , mixture_data , exhaustive , true);
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Mixture et de la structure de donnees associee
 *  dans un fichier au format tableur.
 *
 *  arguments : stream, pointeur sur les vecteurs observees.
 *
 *--------------------------------------------------------------*/

ostream& Mixture::spreadsheet_write(ostream &os , const MixtureData *vec) const

{
  register int i , j , k;
  int variable;
  double **distance;
  FrequencyDistribution *marginal_dist = NULL , **observation_dist = NULL;
  Histogram *marginal_histo = NULL , **observation_histo = NULL;


  os << STAT_word[STATW_MIXTURE] << "\t" << nb_component << "\t" << STAT_word[STATW_DISTRIBUTIONS] << endl;

  // ecriture des poids

  os << "\n" << STAT_word[STATW_WEIGHTS] << endl;

  for (i = 0;i < nb_component;i++) {
    os << weight->mass[i] << "\t";
  }
  os << endl;

  // ecriture des lois associees a chaque processus d'observation

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

    // ecriture de la quantite d'information des vecteurs dans le cas iid

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

    // ecriture des vraisemblances des vecteurs

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


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Mixture dans un fichier au format tableur.
 *
 *  arguments : reference sur un objet StatError, path.
 *
 *--------------------------------------------------------------*/

bool Mixture::spreadsheet_write(StatError &error , const char *path) const

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
    spreadsheet_write(out_file , mixture_data);
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Sortie Gnuplot d'un objet Mixture et de la structure de donnees associee.
 *
 *  arguments : prefixe des fichiers, titre des figures,
 *              pointeur sur les vecteurs observees.
 *
 *--------------------------------------------------------------*/

bool Mixture::plot_write(const char *prefix , const char *title ,
                         const MixtureData *vec) const

{
  bool status;
  register int i;
  int variable , nb_value = I_DEFAULT;
  double *empirical_cdf[2];
  FrequencyDistribution *marginal_dist = NULL , **observation_dist = NULL;
  Histogram *marginal_histo = NULL , **observation_histo = NULL;
  ostringstream data_file_name;


  // ecriture du fichier de donnees

  data_file_name << prefix << 0 << ".dat";
  status = weight->plot_print((data_file_name.str()).c_str() , (vec ? vec->marginal_distribution[0] : NULL));

  if (status) {

    // ecriture du fichier de commandes et du fichier d'impression

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


/*--------------------------------------------------------------*
 *
 *  Sortie Gnuplot d'un objet Mixture.
 *
 *  arguments : reference sur un objet StatError, prefixe des fichiers,
 *              titre des figures.
 *
 *--------------------------------------------------------------*/

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


/*--------------------------------------------------------------*
 *
 *  Sortie graphique d'un objet Mixture et de la structure
 *  de donnees associee.
 *
 *  argument : pointeur sur les vecteurs observees.
 *
 *--------------------------------------------------------------*/

MultiPlotSet* Mixture::get_plotable(const MixtureData *vec) const

{
  register int i , j;
  int nb_plot_set , index , variable;
  FrequencyDistribution *marginal_dist = NULL , **observation_dist = NULL;
  Histogram *marginal_histo = NULL , **observation_histo = NULL;
  MultiPlotSet *plot_set;
  ostringstream title , legend;


  // calcul du nombre de vues

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

    // 1ere vue : poids

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


/*--------------------------------------------------------------*
 *
 *  Sortie graphique d'un objet Mixture.
 *
 *--------------------------------------------------------------*/

MultiPlotSet* Mixture::get_plotable() const

{
  return get_plotable(mixture_data);
}


/*--------------------------------------------------------------*
 *
 *  Calcul du nombre de parametres independants d'un objet Mixture.
 *
 *  argument : probabilite minimum.
 *
 *--------------------------------------------------------------*/

int Mixture::nb_parameter_computation(double min_probability) const

{
  register int i;
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


/*--------------------------------------------------------------*
 *
 *  Constructeur par defaut de la classe MixtureData.
 *
 *--------------------------------------------------------------*/

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


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe MixtureData.
 *
 *  arguments : nombre d'individus, nombre de variables,
 *              type de chaque variable, flag initialisation.
 *
 *--------------------------------------------------------------*/

MixtureData::MixtureData(int inb_vector , int inb_variable ,
                         int *itype , bool init_flag)
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


/*--------------------------------------------------------------*
 *
 *  Construction d'un objet MixtureData a partir d'un objet Vectors.
 *
 *  arguments : reference sur un objet Vectors, type de transformation
 *              ('c' : copie, 'a' : ajout d'une variable d'etat).
 *
 *--------------------------------------------------------------*/

MixtureData::MixtureData(const Vectors &vec , char transform)
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


/*--------------------------------------------------------------*
 *
 *  Copie d'un objet MixtureData.
 *
 *  arguments : reference sur un objet MixtureData,
 *              flag copie de l'objet Mixture.
 *
 *--------------------------------------------------------------*/

void MixtureData::copy(const MixtureData &vec , bool model_flag)

{
  register int i , j;


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


/*--------------------------------------------------------------*
 *
 *  Destructeur des champs de la classe MixtureData.
 *
 *--------------------------------------------------------------*/

void MixtureData::remove()

{
  register int i , j;


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


/*--------------------------------------------------------------*
 *
 *  Destructeur de la classe MixtureData.
 *
 *--------------------------------------------------------------*/

MixtureData::~MixtureData()

{
  remove();
}


/*--------------------------------------------------------------*
 *
 *  Operateur d'assignement de la classe MixtureData.
 *
 *  argument : reference sur un objet MixtureData.
 *
 *--------------------------------------------------------------*/

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


/*--------------------------------------------------------------*
 *
 *  Initialisation de la 1ere variable.
 *
 *  argument : type de la 1ere variable (STATE / INT_VALUE / REAL_VALUE).
 *
 *--------------------------------------------------------------*/

void MixtureData::state_variable_init(int itype)

{
  register int i , j;


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


/*--------------------------------------------------------------*
 *
 *  Extraction d'une composante empirique.
 *
 *  arguments : reference sur un objet StatError, variable,
 *              indice de la composante.
 *
 *--------------------------------------------------------------*/

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


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet MixtureData.
 *
 *  arguments : stream, flag niveau de detail.
 *
 *--------------------------------------------------------------*/

ostream& MixtureData::ascii_write(ostream &os , bool exhaustive) const

{
  if (mixture) {
    mixture->ascii_write(os , this , exhaustive , false);
  }

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet MixtureData dans un fichier.
 *
 *  arguments : reference sur un objet StatError, path,
 *              flag niveau de detail.
 *
 *--------------------------------------------------------------*/

bool MixtureData::ascii_write(StatError &error , const char *path , bool exhaustive) const

{
  bool status = false;


  if (mixture) {
    ofstream out_file(path);

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


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet MixtureData.
 *
 *  arguments : stream, flag niveau de detail.
 *
 *--------------------------------------------------------------*/

ostream& MixtureData::ascii_data_write(ostream &os , bool exhaustive) const

{
  Vectors::ascii_write(os , exhaustive , false);
  ascii_print(os , false , posterior_probability , entropy);

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet MixtureData dans un fichier.
 *
 *  arguments : reference sur un objet StatError, path,
 *              flag niveau de detail.
 *
 *--------------------------------------------------------------*/

bool MixtureData::ascii_data_write(StatError &error , const char *path ,
                                   bool exhaustive) const

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
    Vectors::ascii_write(out_file , exhaustive , true);
    ascii_print(out_file , true , posterior_probability , entropy);
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet MixtureData dans un fichier au format tableur.
 *
 *  arguments : reference sur un objet StatError, path.
 *
 *--------------------------------------------------------------*/

bool MixtureData::spreadsheet_write(StatError &error , const char *path) const

{
  bool status = false;


  if (mixture) {
    ofstream out_file(path);

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


/*--------------------------------------------------------------*
 *
 *  Sortie Gnuplot d'un objet MixtureData.
 *
 *  arguments : reference sur un objet StatError, prefixe des fichiers,
 *              titre des figures.
 *
 *--------------------------------------------------------------*/

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


/*--------------------------------------------------------------*
 *
 *  Sortie graphique d'un objet MixtureData.
 *
 *--------------------------------------------------------------*/

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


/*--------------------------------------------------------------*
 *
 *  Calcul de la quantite d'information.
 *
 *--------------------------------------------------------------*/

double MixtureData::information_computation() const

{
  register int i;
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


/*--------------------------------------------------------------*
 *
 *  Accumulation des observations (pour une variable donnee).
 *
 *  arguments : indice de la variable, nombre de composantes.
 *
 *--------------------------------------------------------------*/

void MixtureData::observation_frequency_distribution_computation(int variable , int nb_component)

{
  register int i , j;


  // initialisation des lois empiriques

  for (i = 0;i < nb_component;i++) {
    for (j = 0;j < marginal_distribution[variable]->nb_value;j++) {
      observation_distribution[variable][i]->frequency[j] = 0;
    }
  }

  // mise a jour des lois empiriques

  for (i = 0;i < nb_vector;i++) {
    (observation_distribution[variable][int_vector[i][0]]->frequency[int_vector[i][variable]])++;
  }

  // extraction des caracteristiques des lois empiriques

  for (i = 0;i < nb_component;i++) {
    observation_distribution[variable][i]->nb_value_computation();
    observation_distribution[variable][i]->offset_computation();
    observation_distribution[variable][i]->nb_element_computation();
    observation_distribution[variable][i]->max_computation();
    observation_distribution[variable][i]->mean_computation();
    observation_distribution[variable][i]->variance_computation();
  }
}


/*--------------------------------------------------------------*
 *
 *  Construction des lois empiriques d'observation.
 *
 *  argument : nombre d'etats.
 *
 *--------------------------------------------------------------*/

void MixtureData::build_observation_frequency_distribution(int nb_component)

{
  if ((nb_variable > 1) && (!observation_distribution)) {
    register int i , j;


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


/*--------------------------------------------------------------*
 *
 *  Construction des histogrammes d'observation pour une variable.
 *
 *  arguments : indice de la variable, nombre d'etats, pas de regroupement.
 *
 *--------------------------------------------------------------*/

void MixtureData::build_observation_histogram(int variable , int nb_component , double step)

{
  if ((!observation_histogram[variable]) || (step != observation_histogram[variable][0]->step)) {
    register int i , j;
    double imin_value;


    // construction de l'histogramme

    if (step == D_DEFAULT) {
      step = marginal_histogram[variable]->step;
    }
    imin_value = floor(min_value[variable] / step) * step;

    if (observation_histogram[variable]) {
      for (i = 0;i < nb_component;i++) {
        observation_histogram[variable][i]->nb_category = (int)floor((max_value[variable] - imin_value) / step) + 1;

        delete [] observation_histogram[variable][i]->frequency;
        observation_histogram[variable][i]->frequency = new int[observation_histogram[variable][i]->nb_category];
      }
    }

    else {
      observation_histogram[variable] = new Histogram*[nb_component];

      for (i = 0;i < nb_component;i++) {
        observation_histogram[variable][i] = new Histogram((int)floor((max_value[variable] - imin_value) / step) + 1 , false);

        observation_histogram[variable][i]->nb_element = marginal_distribution[0]->frequency[i];
        observation_histogram[variable][i]->type = type[variable];
      }

      // calcul des valeurs minimums et maximums par etat

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
      observation_histogram[variable][i]->step = step;
      observation_histogram[variable][i]->min_value = imin_value;
      observation_histogram[variable][i]->max_value = ceil(max_value[variable] / step) * step;
    }

    // calcul des frequences

    for (i = 0;i < nb_component;i++) {
      for (j = 0;j < observation_histogram[variable][i]->nb_category;j++) {
        observation_histogram[variable][i]->frequency[j] = 0;
      }
    }

    switch (type[variable]) {

    case INT_VALUE : {
      for (i = 0;i < nb_vector;i++) {
//        (observation_histogram[variable][int_vector[i][0]]->frequency[(int)((int_vector[i][variable] - imin_value) / step)])++;
        (observation_histogram[variable][int_vector[i][0]]->frequency[(int)floor((int_vector[i][variable] - imin_value) / step)])++;
      }
      break;
    }

    case REAL_VALUE : {
      for (i = 0;i < nb_vector;i++) {
//        (observation_histogram[variable][int_vector[i][0]]->frequency[(int)((real_vector[i][variable] - imin_value) / step)]++;
        (observation_histogram[variable][int_vector[i][0]]->frequency[(int)floor((real_vector[i][variable] - imin_value) / step)])++;
      }
      break;
    }
    }

    for (i = 0;i < nb_component;i++) {
      observation_histogram[variable][i]->max_computation();
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Construction des histogrammes d'observation.
 *
 *  argument : nombre d'etats.
 *
 *--------------------------------------------------------------*/

void MixtureData::build_observation_histogram(int nb_component)

{
  if ((nb_variable > 1) && (!observation_histogram)) {
    register int i;


    observation_histogram = new Histogram**[nb_variable];
    observation_histogram[0] = NULL;

    for (i = 1;i < nb_variable;i++) {
      observation_histogram[i] = NULL;
      if (marginal_histogram[i]) {
        build_observation_histogram(i , nb_component , marginal_histogram[i]->step);
      }
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Changement du pas de regroupement de l'histogramme marginal et
 *  des histogrammes d'observation pour une variable donnee.
 *
 *  arguments : reference sur un objet StatError, indice de la variable,
 *              pas de regroupement, valeur minimum.
 *
 *--------------------------------------------------------------*/

bool MixtureData::select_step(StatError &error , int variable ,
                              double step , double imin_value)

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
    if ((step <= 0.) || ((type[variable] != REAL_VALUE) && ((int)step != step))) {
      status = false;
      error.update(STAT_error[STATR_HISTOGRAM_STEP]);
    }
    if ((imin_value != D_INF) && ((imin_value <= min_value[variable] - step) ||
         (imin_value > min_value[variable]) || ((type[variable] != REAL_VALUE) &&
          ((int)imin_value != imin_value)))) {
      status = false;
      error.update(STAT_error[STATR_HISTOGRAM_MIN_VALUE]);
    }
  }

  if (status) {
    build_marginal_histogram(variable , step , imin_value);

    if ((observation_histogram) && (observation_histogram[variable])) {
      build_observation_histogram(variable , marginal_distribution[0]->nb_value , step);
    }
  }

  return status;
}


};  // namespace stat_tool
