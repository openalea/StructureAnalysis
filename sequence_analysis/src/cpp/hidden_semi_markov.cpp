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

#include "stat_tool/stat_label.h"

#include "hidden_semi_markov.h"
#include "sequence_label.h"

using namespace std;
using namespace boost;
using namespace stat_tool;


namespace sequence_analysis {



/*--------------------------------------------------------------*/
/**
 *  \brief Constructor of the SemiMarkov class.
 *
 *  \param[in] pchain             pointer on a Chain object,
 *  \param[in] poccupancy         pointer on a CategoricalSequenceProcess object,
 *  \param[in] inb_output_process number of observation processes,
 *  \param[in] pobservation       pointer on CategoricalProcess objects,
 *  \param[in] length             sequence length,
 *  \param[in] counting_flag      flag on the computation of the counting distributions.
 */
/*--------------------------------------------------------------*/

SemiMarkov::SemiMarkov(const Chain *pchain , const CategoricalSequenceProcess *poccupancy ,
                       int inb_output_process , CategoricalProcess **pobservation ,
                       int length , bool counting_flag)
:SemiMarkovChain(pchain , poccupancy)

{
  int i;


  nb_iterator = 0;
  semi_markov_data = NULL;

  state_process = new CategoricalSequenceProcess(*poccupancy);

  for (i = 0;i < nb_state;i++) {
    if (transition[i][i] < 1.) {
      state_process->absorption[i] = 0.;
    }
    else {
      state_process->absorption[i] = 1.;
    }
  }

  sojourn_type = new state_sojourn_type[nb_state];
  forward = new Forward*[nb_state];

  for (i = 0;i < nb_state;i++) {
    sojourn_type[i] = (state_process->sojourn_time[i] ? SEMI_MARKOVIAN : MARKOVIAN);

    if ((sojourn_type[i] == SEMI_MARKOVIAN) && (stype[i] == RECURRENT)) {
      forward[i] = new Forward(*(state_process->sojourn_time[i]));
    }
    else {
      forward[i] = NULL;
    }
  }

  if (type == EQUILIBRIUM) {
    for (i = 0;i < nb_state;i++) {
      initial[i] = 1. / (double)nb_state;
    }
    initial_probability_computation();
  }

  nb_output_process = inb_output_process;

  categorical_process = new CategoricalSequenceProcess*[nb_output_process];
  for (i = 0;i < nb_output_process;i++) {
    categorical_process[i] = new CategoricalSequenceProcess(*pobservation[i]);
  }

  if (length > COUNTING_MAX_LENGTH) {
    counting_flag = false;
  }
  characteristic_computation(length , counting_flag);
}


/*--------------------------------------------------------------*/
/**
 *  \brief Constructor of the SemiMarkov class.
 *
 *  \param[in] pchain                            pointer on a Chain object,
 *  \param[in] poccupancy                        pointer on a CategoricalSequenceProcess object,
 *  \param[in] inb_output_process                number of observation processes,
 *  \param[in] categorical_observation           pointer on CategoricalProcess objects,
 *  \param[in] discrete_parametric_observation   pointer on DiscreteParametricProcess objects,
 *  \param[in] continuous_parametric_observation pointer on ContinuousParametricProcess objects,
 *  \param[in] length                            sequence length,
 *  \param[in] counting_flag                     flag on the computation of the counting distributions.
 */
/*--------------------------------------------------------------*/

SemiMarkov::SemiMarkov(const Chain *pchain , const CategoricalSequenceProcess *poccupancy ,
                       int inb_output_process , CategoricalProcess **categorical_observation ,
                       DiscreteParametricProcess **discrete_parametric_observation ,
                       ContinuousParametricProcess **continuous_parametric_observation ,
                       int length , bool counting_flag)
:SemiMarkovChain(pchain , poccupancy)

{
  int i;


  nb_iterator = 0;
  semi_markov_data = NULL;

  state_process = new CategoricalSequenceProcess(*poccupancy);

  for (i = 0;i < nb_state;i++) {
    if (transition[i][i] < 1.) {
      state_process->absorption[i] = 0.;
    }
    else {
      state_process->absorption[i] = 1.;
    }
  }

  sojourn_type = new state_sojourn_type[nb_state];
  forward = new Forward*[nb_state];

  for (i = 0;i < nb_state;i++) {
    sojourn_type[i] = (state_process->sojourn_time[i] ? SEMI_MARKOVIAN : MARKOVIAN);

    if ((sojourn_type[i] == SEMI_MARKOVIAN) && (stype[i] == RECURRENT)) {
      forward[i] = new Forward(*(state_process->sojourn_time[i]));
    }
    else {
      forward[i] = NULL;
    }
  }

  if (type == EQUILIBRIUM) {
    for (i = 0;i < nb_state;i++) {
      initial[i] = 1. / (double)nb_state;
    }
    initial_probability_computation();
  }

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

  if (length > COUNTING_MAX_LENGTH) {
    counting_flag = false;
  }
  characteristic_computation(length , counting_flag);
}


/*--------------------------------------------------------------*/
/**
 *  \brief Destructor of the HiddenSemiMarkov class.
 */
/*--------------------------------------------------------------*/

HiddenSemiMarkov::~HiddenSemiMarkov() {}


/*--------------------------------------------------------------*/
/**
 *  \brief Application of a threshold on the probability parameters of a hidden semi-Markov chain.
 *
 *  \param[in] min_probability minimum probability.
 *
 *  \return                    HiddenSemiMarkov object.
 */
/*--------------------------------------------------------------*/

HiddenSemiMarkov* HiddenSemiMarkov::thresholding(double min_probability) const

{
  int i;
  HiddenSemiMarkov *hsmarkov;


  hsmarkov = new HiddenSemiMarkov(*this , false , false);
  hsmarkov->Chain::thresholding(min_probability , true);

  for (i = 0;i < hsmarkov->nb_output_process;i++) {
    if (hsmarkov->categorical_process[i]) {
      hsmarkov->categorical_process[i]->thresholding(min_probability);
    }
  }

  return hsmarkov;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Construction of a HiddenSemiMarkov object from a file.
 *
 *  \param[in] error           reference on a StatError object,
 *  \param[in] path            file path,
 *  \param[in] length          sequence length,
 *  \param[in] counting_flag   flag on the computation of the counting distributions,
 *  \param[in] cumul_threshold threshold on the cumulative parametric distribution functions,
 *  \param[in] old_format      flag format of the observation processes.
 *
 *  \return                    HiddenSemiMarkov object.
 */
/*--------------------------------------------------------------*/

HiddenSemiMarkov* HiddenSemiMarkov::ascii_read(StatError &error , const string path ,
                                               int length , bool counting_flag ,
                                               double cumul_threshold , bool old_format)

{
  string buffer;
  size_t position;
  typedef tokenizer<char_separator<char>> tokenizer;
  char_separator<char> separator(" \t");
  process_type type = DEFAULT_TYPE;
  bool status , lstatus;
  int i;
  int line , nb_output_process , value , index;
  observation_process obs_type;
  const Chain *chain;
  const CategoricalSequenceProcess *occupancy;
  CategoricalProcess **categorical_observation;
  DiscreteParametricProcess **discrete_parametric_observation;
  ContinuousParametricProcess **continuous_parametric_observation;
  HiddenSemiMarkov *hsmarkov;
  ifstream in_file(path.c_str());


  hsmarkov = NULL;
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

        // test (EQUILIBRIUM_)HIDDEN_SEMI-MARKOV_CHAIN keyword

        if (i == 0) {
          if (*token == SEQ_word[SEQW_HIDDEN_SEMI_MARKOV_CHAIN]) {
            type = ORDINARY;
          }
          else if (*token == SEQ_word[SEQW_EQUILIBRIUM_HIDDEN_SEMI_MARKOV_CHAIN]) {
            type = EQUILIBRIUM;
          }
          else {
            status = false;
            ostringstream correction_message;
            correction_message << SEQ_word[SEQW_HIDDEN_SEMI_MARKOV_CHAIN] << " or "
                               << SEQ_word[SEQW_EQUILIBRIUM_HIDDEN_SEMI_MARKOV_CHAIN];
            error.correction_update(STAT_parsing[STATP_KEYWORD] , (correction_message.str()).c_str() , line);
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

      // analysis of the format and reading of the Markov chain

      chain = Chain::parsing(error , in_file , line , type);

      if (chain) {

        // analysis of the format and reading of the state occupancy distributions

        occupancy = CategoricalSequenceProcess::occupancy_parsing(error , in_file , line ,
                                                                  *chain , cumul_threshold);
        if (!occupancy) {
          status = false;
        }

        // analysis of the format and reading of the observation distributions

        if (old_format) {
          categorical_observation = CategoricalProcess::old_parsing(error , in_file , line ,
                                                                    chain->nb_state , nb_output_process);

          if (categorical_observation) {
            if (status) {
              hsmarkov = new HiddenSemiMarkov(chain , occupancy , nb_output_process ,
                                              categorical_observation , length , counting_flag);
            }

            for (i = 0;i < nb_output_process;i++) {
              delete categorical_observation[i];
            }
            delete [] categorical_observation;
          }
        }

        else {
          nb_output_process = I_DEFAULT;

          categorical_observation = NULL;
          discrete_parametric_observation = NULL;
          continuous_parametric_observation = NULL;

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

              // test number of observation processes

              case 0 : {
                lstatus = true;

/*                try {
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

#             ifdef DEBUG
              cout << line << "  " << buffer << endl;
#             endif

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

/*                  try {
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
                                                                                   chain->nb_state ,
                                                                                   HIDDEN_MARKOV , true);
/*                  categorical_observation[index - 1] = CategoricalProcess::parsing(error , in_file , line , pour les donnees de suivi de croissance manguier
                                                                                   chain->nb_state ,
                                                                                   HIDDEN_MARKOV , false); */
                  if (!categorical_observation[index - 1]) {
                    status = false;
                  }
                  break;
                }

                case DISCRETE_PARAMETRIC : {
                  discrete_parametric_observation[index - 1] = DiscreteParametricProcess::parsing(error , in_file , line ,
                                                                                                  chain->nb_state ,
                                                                                                  HIDDEN_MARKOV ,
                                                                                                  cumul_threshold);
                  if (!discrete_parametric_observation[index - 1]) {
                    status = false;
                  }
                  break;
                }

                case CONTINUOUS_PARAMETRIC : {
                  continuous_parametric_observation[index - 1] = ContinuousParametricProcess::parsing(error , in_file , line ,
                                                                                                      chain->nb_state ,
                                                                                                      HIDDEN_MARKOV ,
                                                                                                      AUTOREGRESSIVE_MODEL);
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

#               ifdef DEBUG
                cout << line << " " << buffer << endl;
#               endif

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
              hsmarkov = new HiddenSemiMarkov(chain , occupancy , nb_output_process ,
                                              categorical_observation ,
                                              discrete_parametric_observation ,
                                              continuous_parametric_observation ,
                                              length , counting_flag);
            }

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

        delete chain;
        delete occupancy;
      }
    }
  }

  return hsmarkov;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a HiddenSemiMarkov object in a file.
 *
 *  \param[in,out] os         stream,
 *  \param[in]     exhaustive flag detail level.
 */
/*--------------------------------------------------------------*/

ostream& HiddenSemiMarkov::ascii_write(ostream &os , bool exhaustive) const

{
  SemiMarkov::ascii_write(os , semi_markov_data , exhaustive ,
                          false , true);

//  os << "\nEnd state: " << end_state() << endl;

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a HiddenSemiMarkov object in a file.
 *
 *  \param[in] error      reference on a StatError object,
 *  \param[in] path       file path,
 *  \param[in] exhaustive flag detail level.
 *
 *  \return               error status.
 */
/*--------------------------------------------------------------*/

bool HiddenSemiMarkov::ascii_write(StatError &error , const string path ,
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
    SemiMarkov::ascii_write(out_file , semi_markov_data , exhaustive ,
                            true , true);
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a HiddenSemiMarkov object in a file at the spreadsheet format.
 *
 *  \param[in] error reference on a StatError object,
 *  \param[in] path  file path.
 *
 *  \return          error status.
 */
/*--------------------------------------------------------------*/

bool HiddenSemiMarkov::spreadsheet_write(StatError &error , const string path) const

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
    SemiMarkov::spreadsheet_write(out_file , semi_markov_data , true);
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Search for an end state.
 *
 *  \return end state index.
 */
/*--------------------------------------------------------------*/

int HiddenSemiMarkov::end_state() const

{
  int i , j , k;
  int end_state = I_DEFAULT , output;


  for (i = nb_state - 1;i >= 0;i--) {
    if (stype[i] == ABSORBING) {

#     ifdef DEBUG
      cout << "\nstate: " << i << " | ";
#     endif

      end_state = i;

      for (j = 0;j < nb_output_process;j++) {
        if (categorical_process[j]) {
          for (k = categorical_process[j]->observation[i]->offset;k < categorical_process[j]->observation[i]->nb_value;k++) {
            if (categorical_process[j]->observation[i]->mass[k] == 1.) {
              output = k;

#             ifdef DEBUG
              cout << "output: " << output << " | ";
#             endif

              break;
            }
          }

          if (k < categorical_process[j]->observation[i]->nb_value) {
            for (k = 0;k < nb_state;k++) {
              if ((k != i) && (output >= categorical_process[j]->observation[k]->offset) &&
                  (output < categorical_process[j]->observation[k]->nb_value) &&
                  (categorical_process[j]->observation[k]->mass[output] > 0.)) {
                end_state = I_DEFAULT;
                break;
              }
            }
            if (end_state == I_DEFAULT) {
              break;
            }
          }

          else {
            end_state = I_DEFAULT;
            break;
          }
        }

        else {
          for (k = discrete_parametric_process[j]->observation[i]->offset;k < discrete_parametric_process[j]->observation[i]->nb_value;k++) {
            if (discrete_parametric_process[j]->observation[i]->mass[k] == 1.) {
              output = k;
              break;
            }
          }

          if (k < discrete_parametric_process[j]->observation[i]->nb_value) {
            for (k = 0;k < nb_state;k++) {
              if ((k != i) && (output >= discrete_parametric_process[j]->observation[k]->offset) &&
                  (output < discrete_parametric_process[j]->observation[k]->nb_value) &&
                  (discrete_parametric_process[j]->observation[k]->mass[output] > 0.)) {
                end_state = I_DEFAULT;
                break;
              }
            }
            if (end_state == I_DEFAULT) {
              break;
            }
          }

          else {
            end_state = I_DEFAULT;
            break;
          }
        }
      }

#     ifdef DEBUG
      cout << "end state: " << end_state << endl;
#     endif

      if (end_state == i) {
        break;
      }
    }
  }

  return end_state;
}


};  // namespace sequence_analysis
