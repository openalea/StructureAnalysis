/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       V-Plants: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2016 CIRAD/INRA/Inria Virtual Plants
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



#include <string>
#include <sstream>

#include "tool/rw_tokenizer.h"
#include "tool/rw_cstring.h"
#include "tool/rw_locale.h"

#include "stat_tool/stat_label.h"

#include "hidden_variable_order_markov.h"
#include "sequence_label.h"

using namespace std;
using namespace stat_tool;


namespace sequence_analysis {



/*--------------------------------------------------------------*/
/**
 *  \brief Constructor of the VariableOrderMarkov class.
 *
 *  \param[in] pmarkov                           pointer on a VariableOrderMarkovChain object,
 *  \param[in] inb_output_process                number of observation processes,
 *  \param[in] categorical_observation           pointer on CategoricalProcess objects,
 *  \param[in] discrete_parametric_observation   pointer on DiscreteParametricProcess objects,
 *  \param[in] continuous_parametric_observation pointer on ContinuousParametricProcess objects,
 *  \param[in] length                            sequence length.
 */
/*--------------------------------------------------------------*/

VariableOrderMarkov::VariableOrderMarkov(const VariableOrderMarkovChain *pmarkov , int inb_output_process ,
                                         CategoricalProcess **categorical_observation ,
                                         DiscreteParametricProcess **discrete_parametric_observation ,
                                         ContinuousParametricProcess **continuous_parametric_observation ,
                                         int length)

{
  register int i;


  build(*pmarkov);

  nb_iterator = 0;
  markov_data = NULL;

  nb_output_process = inb_output_process;

  categorical_process = new CategoricalSequenceProcess*[nb_output_process];
  discrete_parametric_process = new DiscreteParametricProcess*[nb_output_process];
  continuous_parametric_process = new ContinuousParametricProcess*[nb_output_process];

  for (i = 0;i < nb_output_process;i++) {
    if (categorical_observation[i]) {
      categorical_process[i] = new CategoricalSequenceProcess(*categorical_observation[i]);
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

  characteristic_computation(length , true);
}


/*--------------------------------------------------------------*/
/**
 *  \brief Destructor of the HiddenVariableOrderMarkov class.
 */
/*--------------------------------------------------------------*/

HiddenVariableOrderMarkov::~HiddenVariableOrderMarkov() {}


/*--------------------------------------------------------------*/
/**
 *  \brief Application of a threshold on the probability parameters of
 *         a hidden variable-order Markov chain.
 *
 *  \param[in] min_probability minimum probability.
 *
 *  \return                    HiddenVariableOrderMarkov object.
 */
/*--------------------------------------------------------------*/

HiddenVariableOrderMarkov* HiddenVariableOrderMarkov::thresholding(double min_probability) const

{
  register int i;
  HiddenVariableOrderMarkov *hmarkov;


  hmarkov = new HiddenVariableOrderMarkov(*this , false);
  hmarkov->VariableOrderMarkovChain::thresholding(min_probability);

  for (i = 0;i < hmarkov->nb_output_process;i++) {
    if (hmarkov->categorical_process[i]) {
      hmarkov->categorical_process[i]->thresholding(min_probability);
    }
  }

  return hmarkov;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Construction of a HiddenVariableOrderMarkov object from a file.
 *
 *  \param[in] error           reference on a StatError object,
 *  \param[in] path            file path,
 *  \param[in] length          sequence length,
 *  \param[in] cumul_threshold threshold on the cumulative parametric distribution functions.
 *
 *  \return                    HiddenVariableOrderMarkov object.
 */
/*--------------------------------------------------------------*/

HiddenVariableOrderMarkov* HiddenVariableOrderMarkov::ascii_read(StatError &error ,
                                                                 const string path , int length ,
                                                                 double cumul_threshold)

{
  RWLocaleSnapshot locale("en");
  RWCString buffer , token;
  size_t position;
  process_type type = DEFAULT_TYPE;
  bool status , lstatus;
  register int i;
  int line , nb_output_process , index;
  observation_process obs_type;
  long value;
  const VariableOrderMarkovChain *imarkov;
  CategoricalProcess **categorical_observation;
  DiscreteParametricProcess **discrete_parametric_observation;
  ContinuousParametricProcess **continuous_parametric_observation;
  HiddenVariableOrderMarkov *hmarkov;
  ifstream in_file(path.c_str());


  hmarkov = NULL;
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

        // test (EQUILIBRIUM_)HIDDEN_MARKOV_CHAIN keyword

        if (i == 0) {
          if (token == SEQ_word[SEQW_HIDDEN_MARKOV_CHAIN]) {
            type = ORDINARY;
          }
          else if (token == SEQ_word[SEQW_EQUILIBRIUM_HIDDEN_MARKOV_CHAIN]) {
            type = EQUILIBRIUM;
          }
          else {
            status = false;
            ostringstream correction_message;
            correction_message << SEQ_word[SEQW_HIDDEN_MARKOV_CHAIN] << " or "
                               << SEQ_word[SEQW_EQUILIBRIUM_HIDDEN_MARKOV_CHAIN];
            error.correction_update(STAT_parsing[STATP_KEYWORD] ,
                                    (correction_message.str()).c_str() , line);
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

    if (type != DEFAULT_TYPE) {

      // analysis of the format and reading of the variable-order Markov chain

      imarkov = VariableOrderMarkovChain::parsing(error , in_file , line , type);

      // analysis of the format and reading of the observation distributions

      if (imarkov) {
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

            // test number of observation processes

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

            // test OUTPUT_PROCESS(ES) keyword

            case 1 : {
              if (token != STAT_word[nb_output_process == 1 ? STATW_OUTPUT_PROCESS : STATW_OUTPUT_PROCESSES]) {
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

              // test OUTPUT_PROCESS keyword

              case 0 : {
                if (token == STAT_word[STATW_OUTPUT_PROCESS]) {
                  index++;
                }
                else {
                  status = false;
                  error.correction_update(STAT_parsing[STATP_KEYWORD] , STAT_word[STATW_OUTPUT_PROCESS] , line , i + 1);
                }
                break;
              }

              // test observation process index

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

              // test separator

              case 2 : {
                if (token != ":") {
                  status = false;
                  error.update(STAT_parsing[STATP_SEPARATOR] , line , i + 1);
                }
                break;
              }

              // test CATEGORICAL/DISCRETE_PARAMETRIC/CONTINUOUS_PARAMETRIC keyword

              case 3 : {
                if ((token == STAT_word[STATW_CATEGORICAL]) ||
                    (token == STAT_word[STATW_NONPARAMETRIC])) {
                  obs_type = CATEGORICAL_PROCESS;
                }
                else if ((token == STAT_word[STATW_DISCRETE_PARAMETRIC]) ||
                         (token == STAT_word[STATW_PARAMETRIC])) {
                  obs_type = DISCRETE_PARAMETRIC;
                }
                else if (token == STAT_word[STATW_CONTINUOUS_PARAMETRIC]) {
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
                                                                                 ((Chain*)imarkov)->nb_state ,
                                                                                 HIDDEN_MARKOV , true);
                if (!categorical_observation[index - 1]) {
                  status = false;
                }
                break;
              }

              case DISCRETE_PARAMETRIC : {
                discrete_parametric_observation[index - 1] = DiscreteParametricProcess::parsing(error , in_file , line ,
                                                                                                ((Chain*)imarkov)->nb_state ,
                                                                                                HIDDEN_MARKOV ,
                                                                                                cumul_threshold);
                if (!discrete_parametric_observation[index - 1]) {
                  status = false;
                }
                break;
              }

              case CONTINUOUS_PARAMETRIC : {
                continuous_parametric_observation[index - 1] = ContinuousParametricProcess::parsing(error , in_file , line ,
                                                                                                    ((Chain*)imarkov)->nb_state ,
                                                                                                    HIDDEN_MARKOV ,
                                                                                                    ZERO_INFLATED_GAMMA);
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
            hmarkov = new HiddenVariableOrderMarkov(imarkov , nb_output_process ,
                                                    categorical_observation ,
                                                    discrete_parametric_observation ,
                                                    continuous_parametric_observation , length);

#           ifdef DEBUG
            hmarkov->ascii_write(cout);
#           endif

          }

          delete imarkov;

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

  return hmarkov;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a HiddenVariableOrderMarkov object in a file.
 *
 *  \param[in,out] os         stream,
 *  \param[in]     exhaustive flag detail level.
 */
/*--------------------------------------------------------------*/

ostream& HiddenVariableOrderMarkov::ascii_write(ostream &os , bool exhaustive) const

{
  VariableOrderMarkov::ascii_write(os , markov_data , exhaustive , false , true);

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a HiddenVariableOrderMarkov object in a file.
 *
 *  \param[in] error      reference on a StatError object,
 *  \param[in] path       file path,
 *  \param[in] exhaustive flag detail level.
 *
 *  \return               error status.
 */
/*--------------------------------------------------------------*/

bool HiddenVariableOrderMarkov::ascii_write(StatError &error , const string path ,
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
    VariableOrderMarkov::ascii_write(out_file , markov_data , exhaustive , true , true);
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a HiddenVariableOrderMarkov object in a file
 *         at the spreadsheet format.
 *
 *  \param[in] error reference on a StatError object,
 *  \param[in] path  file path.
 *
 *  \return          error status.
 */
/*--------------------------------------------------------------*/

bool HiddenVariableOrderMarkov::spreadsheet_write(StatError &error ,
                                                  const string path) const

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
    VariableOrderMarkov::spreadsheet_write(out_file , markov_data , true);
  }

  return status;
}


};  // namespace sequence_analysis
