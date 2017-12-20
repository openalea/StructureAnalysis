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



#include <string>
#include <sstream>
#include <iomanip>

#include <boost/tokenizer.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/classification.hpp>

#include "stat_tool/stat_label.h"

#include "semi_markov.h"
#include "sequence_label.h"

using namespace std;
using namespace boost;
using namespace stat_tool;


namespace sequence_analysis {



/*--------------------------------------------------------------*/
/**
 *  \brief Default constructor of the SemiMarkovChain class.
 */
/*--------------------------------------------------------------*/

SemiMarkovChain::SemiMarkovChain()

{
  sojourn_type = NULL;
  state_process = NULL;
  forward = NULL;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Constructor of the SemiMarkovChain class.
 *
 *  \param[in] itype     process type (ORDINARY/EQUILIBRIUM),
 *  \param[in] inb_state number of states.
 */
/*--------------------------------------------------------------*/

SemiMarkovChain::SemiMarkovChain(process_type itype , int inb_state)
:Chain(itype , inb_state)

{
  sojourn_type = NULL;
  state_process = new CategoricalSequenceProcess(nb_state , nb_state , false);
  forward = NULL;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Constructor of the SemiMarkovChain class.
 *
 *  \param[in] pchain     pointer on a Chain object,
 *  \param[in] poccupancy pointer on a CategoricalSequenceProcess object.
 */
/*--------------------------------------------------------------*/

SemiMarkovChain::SemiMarkovChain(const Chain *pchain , const CategoricalSequenceProcess *poccupancy)
:Chain(*pchain)

{
  int i;


  sojourn_type = new state_sojourn_type[nb_state];

  state_process = new CategoricalSequenceProcess(*poccupancy);
  for (i = 0;i < nb_state;i++) {
    if (transition[i][i] < 1.) {
      state_process->absorption[i] = 0.;
    }
    else {
      state_process->absorption[i] = 1.;
    }
  }

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
}


/*--------------------------------------------------------------*/
/**
 *  \brief Copy of a SemiMarkovChain object.
 *
 *  \param[in] smarkov reference on a SemiMarkovChain object,
 *  \param[in] param   parameter (if > 0: number of allocated values for the state occupancy distributions).
 */
/*--------------------------------------------------------------*/

void SemiMarkovChain::copy(const SemiMarkovChain &smarkov , int param)

{
  int i;


  sojourn_type = new state_sojourn_type[nb_state];
  for (i = 0;i < nb_state;i++) {
    sojourn_type[i] = smarkov.sojourn_type[i];
  }

  forward = new Forward*[nb_state];

  for (i = 0;i < nb_state;i++) {
    if (smarkov.forward[i]) {
      forward[i] = new Forward(*(smarkov.forward[i]) , param);
    }
    else {
      forward[i] = NULL;
    }
  }

  switch (param) {
  case I_DEFAULT :
    state_process = new CategoricalSequenceProcess(*(smarkov.state_process));
    break;
  case 0 :
    state_process = new CategoricalSequenceProcess(*(smarkov.state_process) , INIT_OCCUPANCY , I_DEFAULT);
    break;
  default :
    state_process = new CategoricalSequenceProcess(*(smarkov.state_process) , INIT_OCCUPANCY , param);
    break;
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Destruction of the data members of a SemiMarkovChain object.
 */
/*--------------------------------------------------------------*/

void SemiMarkovChain::remove()

{
  int i;


  delete [] sojourn_type;

  delete state_process;

  if (forward) {
    for (i = 0;i < nb_state;i++) {
      delete forward[i];
    }
    delete [] forward;
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Destructor of the SemiMarkovChain class.
 */
/*--------------------------------------------------------------*/

SemiMarkovChain::~SemiMarkovChain()

{
  remove();
}


/*--------------------------------------------------------------*/
/**
 *  \brief Assignment operator of the SemiMarkovChain class.
 *
 *  \param[in] smarkov reference on a SemiMarkovChain object.
 *
 *  \return            SemiMarkovChain object.
 */
/*--------------------------------------------------------------*/

SemiMarkovChain& SemiMarkovChain::operator=(const SemiMarkovChain &smarkov)

{
  if (&smarkov != this) {
    remove();
    Chain::remove();

    Chain::copy(smarkov);
    copy(smarkov);
  }

  return *this;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the number of parameters of a SemiMarkovChain object.
 *
 *  \param[in] min_probability minimum probability.
 *
 *  \return                    number of parameters.
 */
/*--------------------------------------------------------------*/

int SemiMarkovChain::nb_parameter_computation(double min_probability) const

{
  int i;
  int nb_parameter = Chain::nb_parameter_computation(min_probability);


  for (i = 0;i < nb_state;i++) {
    if (sojourn_type[i] == SEMI_MARKOVIAN) {
      nb_parameter += state_process->sojourn_time[i]->nb_parameter_computation();
      if (state_process->sojourn_time[i]->inf_bound == 1) {
        nb_parameter--;
      }
    }
  }

  return nb_parameter;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Default constructor of the SemiMarkov class.
 */
/*--------------------------------------------------------------*/

SemiMarkov::SemiMarkov()

{
  nb_iterator = 0;
  semi_markov_data = NULL;

  nb_output_process = 0;
  categorical_process = NULL;
  discrete_parametric_process = NULL;
  continuous_parametric_process = NULL;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Constructor of the SemiMarkov class.
 *
 *  \param[in] itype              process type (ORDINARY/EQUILIBRIUM),
 *  \param[in] inb_state          number of states,
 *  \param[in] inb_output_process number of observation processes,
 *  \param[in] nb_value           number of observed values for each observation process.
 */
/*--------------------------------------------------------------*/

SemiMarkov::SemiMarkov(process_type itype , int inb_state , int inb_output_process , int *nb_value)
:SemiMarkovChain(itype , inb_state)

{
  int i;


  nb_iterator = 0;
  semi_markov_data = NULL;

  nb_output_process = inb_output_process;

  categorical_process = new CategoricalSequenceProcess*[nb_output_process];
  discrete_parametric_process = new DiscreteParametricProcess*[nb_output_process];
  continuous_parametric_process = new ContinuousParametricProcess*[nb_output_process];

  for (i = 0;i < nb_output_process;i++) {
    if (nb_value[i] == I_DEFAULT) {
      categorical_process[i] = NULL;
      discrete_parametric_process[i] = NULL;
      continuous_parametric_process[i] = new ContinuousParametricProcess(nb_state);
    }

    else if (nb_value[i] <= NB_OUTPUT) {
      categorical_process[i] = new CategoricalSequenceProcess(nb_state , nb_value[i] , true);
      discrete_parametric_process[i] = NULL;
      continuous_parametric_process[i] = NULL;
    }

    else {
      categorical_process[i] = NULL;
      discrete_parametric_process[i] = new DiscreteParametricProcess(nb_state , (int)(nb_value[i] * SAMPLE_NB_VALUE_COEFF));
      continuous_parametric_process[i] = NULL;
    }
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Constructor of the SemiMarkov class.
 *
 *  \param[in] pchain        pointer on a Chain object,
 *  \param[in] poccupancy    pointer on a CategoricalSequenceProcess object,
 *  \param[in] pobservation  pointer on a CategoricalProcess object,
 *  \param[in] length        sequence length,
 *  \param[in] counting_flag flag on the computation of the counting distributions.
 */
/*--------------------------------------------------------------*/

SemiMarkov::SemiMarkov(const Chain *pchain , const CategoricalSequenceProcess *poccupancy ,
                       const CategoricalProcess *pobservation ,
                       int length , bool counting_flag)
:SemiMarkovChain(pchain , poccupancy)

{
  int i;


  nb_iterator = 0;
  semi_markov_data = NULL;

  nb_output_process = (pobservation ? 1 : 0);

  if (nb_output_process == 1) {
    categorical_process = new CategoricalSequenceProcess*[nb_output_process];
    categorical_process[0] = new CategoricalSequenceProcess(*pobservation);
  }
  else {
    categorical_process = NULL;
  }

  discrete_parametric_process = NULL;
  continuous_parametric_process = NULL;

  if (length > COUNTING_MAX_LENGTH) {
    counting_flag = false;
  }
  characteristic_computation(length , counting_flag);
}


/*--------------------------------------------------------------*/
/**
 *  \brief Copy of a SemiMarkov object.
 *
 *  \param[in] smarkov   reference on a SemiMarkov object,
 *  \param[in] data_flag flag copy of the included SemiMarkovData object,
 *  \param[in] param     parameter.
 */
/*--------------------------------------------------------------*/

void SemiMarkov::copy(const SemiMarkov &smarkov , bool data_flag , int param)

{
  int i;


  nb_iterator = 0;

  if ((data_flag) && (smarkov.semi_markov_data)) {
    semi_markov_data = new SemiMarkovData(*(smarkov.semi_markov_data) , false);
  }
  else {
    semi_markov_data = NULL;
  }

  nb_output_process = smarkov.nb_output_process;

  if (smarkov.categorical_process) {
    categorical_process = new CategoricalSequenceProcess*[nb_output_process];

    switch (param) {

    case I_DEFAULT : {
      for (i = 0;i < nb_output_process;i++) {
        if (smarkov.categorical_process[i]) {
          categorical_process[i] = new CategoricalSequenceProcess(*(smarkov.categorical_process[i]));
        }
        else {
          categorical_process[i] = NULL;
        }
      }
      break;
    }

    default : {
      for (i = 0;i < nb_output_process;i++) {
        if (smarkov.categorical_process[i]) {
          categorical_process[i] = new CategoricalSequenceProcess(*(smarkov.categorical_process[i]) ,
                                                                  CATEGORICAL_SEQUENCE_PROCESS_COPY , false);
        }
        else {
          categorical_process[i] = NULL;
        }
      }
      break;
    }
    }
  }

  else {
    categorical_process = NULL;
  }

  if (smarkov.discrete_parametric_process) {
    discrete_parametric_process = new DiscreteParametricProcess*[nb_output_process];

    for (i = 0;i < nb_output_process;i++) {
      if (smarkov.discrete_parametric_process[i]) {
        discrete_parametric_process[i] = new DiscreteParametricProcess(*(smarkov.discrete_parametric_process[i]));
      }
      else {
        discrete_parametric_process[i] = NULL;
      }
    }
  }

  else {
    discrete_parametric_process = NULL;
  }

  if (smarkov.continuous_parametric_process) {
    continuous_parametric_process = new ContinuousParametricProcess*[nb_output_process];

    for (i = 0;i < nb_output_process;i++) {
      if (smarkov.continuous_parametric_process[i]) {
        continuous_parametric_process[i] = new ContinuousParametricProcess(*(smarkov.continuous_parametric_process[i]));
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
 *  \brief Destruction of the data members of a SemiMarkov object.
 */
/*--------------------------------------------------------------*/

void SemiMarkov::remove()

{
  int i;


  delete semi_markov_data;

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
 *  \brief Destructor of the SemiMarkov class.
 */
/*--------------------------------------------------------------*/

SemiMarkov::~SemiMarkov()

{
  remove();
}


/*--------------------------------------------------------------*/
/**
 *  \brief Destruction of a SemiMarkov object taking account of
 *         the number of iterators pointing to it.
 */
/*--------------------------------------------------------------*/

void SemiMarkov::conditional_delete()

{
  if (nb_iterator == 0) {
    delete this;
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Assignment operator of the SemiMarkov class.
 *
 *  \param[in] smarkov reference on a SemiMarkov object.
 *
 *  \return            SemiMarkov object.
 */
/*--------------------------------------------------------------*/

SemiMarkov& SemiMarkov::operator=(const SemiMarkov &smarkov)

{
  if ((&smarkov != this) && (nb_iterator == 0)) {
    remove();
    SemiMarkovChain::remove();
    Chain::remove();

    Chain::copy(smarkov);
    SemiMarkovChain::copy(smarkov);
    copy(smarkov);
  }

  return *this;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Extraction of a distribution.
 *
 *  \param[in] error     reference on a StatError object,
 *  \param[in] dist_type distribution type,
 *  \param[in] variable  variable index,
 *  \param[in] value     state or observation.
 *
 *  \return              DiscreteParametricModel object.
 */
/*--------------------------------------------------------------*/

DiscreteParametricModel* SemiMarkov::extract(StatError &error , process_distribution dist_type ,
                                             int variable , int value) const

{
  bool status = true;
  int hvariable;
  Distribution *pdist;
  DiscreteParametric *pparam;
  DiscreteParametricModel *dist;
  FrequencyDistribution *phisto;
  CategoricalSequenceProcess *process;


  dist = NULL;
  error.init();

  pdist = NULL;
  pparam = NULL;

  if (dist_type == OBSERVATION) {
    if ((variable < 1) || (variable > nb_output_process)) {
      status = false;
      error.update(STAT_error[STATR_OUTPUT_PROCESS_INDEX]);
    }

    else {
      if ((value < 0) || (value >= nb_state)) {
        status = false;
        ostringstream error_message;
        error_message << STAT_label[STATL_STATE] << " " << value << " "
                      << STAT_error[STATR_NOT_PRESENT];
        error.update((error_message.str()).c_str());
      }

      else {
        if (categorical_process[variable - 1]) {
          pdist = categorical_process[variable - 1]->observation[value];
        }
        else if (discrete_parametric_process[variable - 1]) {
          pparam = discrete_parametric_process[variable - 1]->observation[value];
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
  }

  else {
    if ((variable < 0) || (variable > nb_output_process)) {
      status = false;
      error.update(STAT_error[STATR_OUTPUT_PROCESS_INDEX]);
    }

    else {
      if (variable == 0) {
        process = state_process;
      }

      else {
        process = categorical_process[variable - 1];

        if (!process) {
          status = false;
          ostringstream error_message;
          error_message << STAT_label[STATL_VARIABLE] << " " << variable << ": " 
                        << SEQ_error[SEQR_CHARACTERISTICS_NOT_COMPUTED];
          error.update((error_message.str()).c_str());
        }
      }

      if ((process) && ((value < 0) || (value >= process->nb_value))) {
        status = false;
        ostringstream error_message;
        error_message << STAT_label[variable == 0 ? STATL_STATE : STATL_OUTPUT] << " "
                      << value << " " << STAT_error[STATR_NOT_PRESENT];
        error.update((error_message.str()).c_str());
      }

      if (status) {
        switch (dist_type) {
        case FIRST_OCCURRENCE :
          pdist = process->first_occurrence[value];
          break;
        case RECURRENCE_TIME :
          pdist = process->recurrence_time[value];
          break;
        case SOJOURN_TIME :
          pparam = process->sojourn_time[value];
          break;
        case NB_RUN :
          pdist = process->nb_run[value];
          break;
        case NB_OCCURRENCE :
          pdist = process->nb_occurrence[value];
          break;
        }

        if ((!pdist) && (!pparam)) {
          status = false;
          error.update(SEQ_error[SEQR_NON_EXISTING_CHARACTERISTIC_DISTRIBUTION]);
        }
      }
    }
  }

  if (status) {
    phisto = NULL;

    if (semi_markov_data) {
      switch (semi_markov_data->type[0]) {
      case STATE :
        hvariable = variable;
        break;
      case INT_VALUE :
        hvariable = variable - 1;
        break;
      }

      if (hvariable >= 0) {
        switch (dist_type) {

        case OBSERVATION : {
          if ((semi_markov_data->observation_distribution) &&
              (semi_markov_data->observation_distribution[hvariable])) {
            phisto = semi_markov_data->observation_distribution[hvariable][value];
          }
          break;
        }

        case FIRST_OCCURRENCE : {
          phisto = semi_markov_data->characteristics[hvariable]->first_occurrence[value];
          break;
        }

        case RECURRENCE_TIME : {
          if (semi_markov_data->characteristics[hvariable]->recurrence_time[value]->nb_element > 0) {
            phisto = semi_markov_data->characteristics[hvariable]->recurrence_time[value];
          }
          break;
        }

        case SOJOURN_TIME : {
          if (semi_markov_data->characteristics[hvariable]->sojourn_time[value]->nb_element > 0) {
            phisto = semi_markov_data->characteristics[hvariable]->sojourn_time[value];
          }
          break;
        }

        case NB_RUN : {
          phisto = semi_markov_data->characteristics[hvariable]->nb_run[value];
          break;
        }

        case NB_OCCURRENCE : {
          phisto = semi_markov_data->characteristics[hvariable]->nb_occurrence[value];
          break;
        }
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
 *  \brief Extraction of a forward recurrence time distribution.
 *
 *  \param[in] error      reference on a StatError object,
 *  \param[in] state      state,
 *  \param[in] histo_type type of associated frequency distribution.
 *
 *  \return               DiscreteParametricModel object.
 */
/*--------------------------------------------------------------*/

DiscreteParametricModel* SemiMarkov::extract(StatError &error , int state ,
                                             process_distribution histo_type) const

{
  bool status = true;
  Distribution *pdist;
  DiscreteParametricModel *dist;
  FrequencyDistribution *phisto;


  dist = NULL;
  error.init();

  pdist = NULL;

  if ((state < 0) || (state >= nb_state)) {
    status = false;
    ostringstream error_message;
    error_message << STAT_label[STATL_STATE] << " " << state << " "
                  << STAT_error[STATR_NOT_PRESENT];
    error.update((error_message.str()).c_str());
  }

  else {
    pdist = forward[state];

    if (!pdist) {
      status = false;
      error.update(SEQ_error[SEQR_NON_EXISTING_FORWARD_DISTRIBUTION]);
    }

    else {
      phisto = NULL;

      if ((semi_markov_data) && (semi_markov_data->type[0] == STATE)) {
        switch (histo_type) {

        case INITIAL_RUN : {
          if ((semi_markov_data->characteristics[0]->initial_run) &&
              (semi_markov_data->characteristics[0]->initial_run[state]->nb_element > 0)) {
            phisto = semi_markov_data->characteristics[0]->initial_run[state];
          }
          break;
        }

        case FINAL_RUN : {
          if (semi_markov_data->characteristics[0]->final_run[state]->nb_element > 0) {
            phisto = semi_markov_data->characteristics[0]->final_run[state];
          }
          break;
        }
        }
      }

      dist = new DiscreteParametricModel(*pdist , phisto);
    }
  }

  return dist;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Extraction of the SemiMarkovData object included in a SemiMarkov object.
 *
 *  \param[in] error reference on a StatError object.
 *
 *  \return          SemiMarkovData object.
 */
/*--------------------------------------------------------------*/

SemiMarkovData* SemiMarkov::extract_data(StatError &error) const

{
  bool status = true;
  SemiMarkovData *seq;


  seq = NULL;
  error.init();

  if (!semi_markov_data) {
    status = false;
    error.update(STAT_error[STATR_NO_DATA]);
  }
  else if (nb_output_process + 1 != semi_markov_data->nb_variable) {
    status = false;
    error.update(SEQ_error[SEQR_STATE_SEQUENCES]);
  }

  if (status) {
    seq = new SemiMarkovData(*semi_markov_data);
    seq->semi_markov = new SemiMarkov(*this , false);
  }

  return seq;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Application of a threshold on the probability parameters of a semi-Markov chain.
 *
 *  \param[in] min_probability minimum probability.
 *
 *  \return                    SemiMarkov object.
 */
/*--------------------------------------------------------------*/

SemiMarkov* SemiMarkov::thresholding(double min_probability) const

{
  int i;
  SemiMarkov *smarkov;


  smarkov = new SemiMarkov(*this , false , 0);
  smarkov->Chain::thresholding(min_probability , true);

  for (i = 0;i < smarkov->nb_output_process;i++) {
    if (smarkov->categorical_process[i]) {
      smarkov->categorical_process[i]->thresholding(min_probability);
    }
  }

  return smarkov;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Construction of a SemiMarkov object from a file.
 *
 *  \param[in] error           reference on a StatError object,
 *  \param[in] path            file path,
 *  \param[in] length          sequence length,
 *  \param[in] counting_flag   flag on the computation of the counting distributions,
 *  \param[in] cumul_threshold threshold on the state occupancy cumulative distribution functions.
 *
 *  \return                    SemiMarkov object.
 */
/*--------------------------------------------------------------*/

SemiMarkov* SemiMarkov::ascii_read(StatError &error , const string path , int length ,
                                   bool counting_flag , double cumul_threshold)

{
  string buffer;
  size_t position;
  typedef tokenizer<char_separator<char>> tokenizer;
  char_separator<char> separator(" \t");
  process_type type = DEFAULT_TYPE;
  bool status;
  int i;
  int line;
  const Chain *chain;
  const CategoricalSequenceProcess *occupancy;
  const CategoricalProcess *observation;
  SemiMarkov *smarkov;
  ifstream in_file(path.c_str());


  smarkov = NULL;
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

        // test (EQUILIBRIUM_)SEMI-MARKOV_CHAIN keyword

        if (i == 0) {
          if (*token == SEQ_word[SEQW_SEMI_MARKOV_CHAIN]) {
            type = ORDINARY;
          }
          else if (*token == SEQ_word[SEQW_EQUILIBRIUM_SEMI_MARKOV_CHAIN]) {
            type = EQUILIBRIUM;
          }
          else {
            status = false;
            ostringstream correction_message;
            correction_message << SEQ_word[SEQW_SEMI_MARKOV_CHAIN] << " or "
                               << SEQ_word[SEQW_EQUILIBRIUM_SEMI_MARKOV_CHAIN];
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

      // analysis of the format and reading of the Markov chain

      chain = Chain::parsing(error , in_file , line , type);

      if (chain) {

        // analysis of the format and reading of the state occupancy distributions

        occupancy = CategoricalSequenceProcess::occupancy_parsing(error , in_file , line ,
                                                                  *chain , cumul_threshold);
        if (!occupancy) {
          status = false;
        }

        // analysis of the format and reading of the categorical observation distributions

        observation = NULL;

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

            // test OUTPUT_PROCESS keyword

            if (i == 0) {
              if (*token != STAT_word[STATW_OUTPUT_PROCESS]) {
                status = false;
                error.correction_update(STAT_parsing[STATP_KEYWORD] , STAT_word[STATW_OUTPUT_PROCESS] , line);
              }
            }

            i++;
          }

          if (i > 0) {
            if (i != 1) {
              status = false;
              error.update(STAT_parsing[STATP_FORMAT] , line);
            }

            observation = CategoricalProcess::parsing(error , in_file , line , chain->nb_state ,
                                                      HIDDEN_MARKOV , false);
            if (!observation) {
              status = false;
            }

            break;
          }
        }

        while (getline(in_file , buffer)) {
          line++;

#         ifdef DEBUG
          cout << line << "  " << buffer << endl;
#         endif

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
          smarkov = new SemiMarkov(chain , occupancy , observation , length , counting_flag);
        }

        delete chain;
        delete occupancy;
        delete observation;
      }
    }
  }

  return smarkov;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing on a single line of a SemiMarkov object.
 *
 *  \param[in,out] os stream.
 */
/*--------------------------------------------------------------*/

ostream& SemiMarkov::line_write(ostream &os) const

{
  os << nb_state << " " << STAT_word[STATW_STATES];

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a SemiMarkov object and the associated data structure.
 *
 *  \param[in,out] os         stream,
 *  \param[in]     seq        pointer on a SemiMarkovData object,
 *  \param[in]     exhaustive flag detail level,
 *  \param[in]     file_flag  flag file,
 *  \param[in]     hidden     flag hidden model.
 */
/*--------------------------------------------------------------*/

ostream& SemiMarkov::ascii_write(ostream &os , const SemiMarkovData *seq ,
                                 bool exhaustive , bool file_flag , bool hidden) const

{
  int i , j , k;
  int buff , width , variable;
  double **distance;
  FrequencyDistribution *marginal_dist = NULL , **observation_dist = NULL;
  Histogram *marginal_histo = NULL , **observation_histo = NULL;
  SequenceCharacteristics *characteristics = NULL;
  ios_base::fmtflags format_flags;


  format_flags = os.setf(ios::left , ios::adjustfield);

  if (hidden) {
    switch (type) {
    case ORDINARY :
      os << SEQ_word[SEQW_HIDDEN_SEMI_MARKOV_CHAIN] << endl;
      break;
    case EQUILIBRIUM :
      os << SEQ_word[SEQW_EQUILIBRIUM_HIDDEN_SEMI_MARKOV_CHAIN] << endl;
      break;
    }
  }

  else {
    switch (type) {
    case ORDINARY :
      os << SEQ_word[SEQW_SEMI_MARKOV_CHAIN] << endl;
      break;
    case EQUILIBRIUM :
      os << SEQ_word[SEQW_EQUILIBRIUM_SEMI_MARKOV_CHAIN] << endl;
      break;
    }
  }

  // writing of the Markov chain parameters

  ascii_print(os , file_flag);

  // writing of the state occupancy distributions

  if ((seq) && (seq->type[0] == STATE)) {
    characteristics = seq->characteristics[0];
  }
  else {
    characteristics = NULL;
  }

  state_process->ascii_print(os , 0 , NULL , NULL , characteristics ,
                             exhaustive , file_flag , forward);

  if (hidden) {
    for (i = 0;i < nb_output_process;i++) {
      if (discrete_parametric_process[i]) {
        if (discrete_parametric_process[i]->weight) {
          width = column_width(nb_state , discrete_parametric_process[i]->weight->mass);
        }
        else {
          width = 0;
        }
        if (discrete_parametric_process[i]->restoration_weight) {
          buff = column_width(nb_state , discrete_parametric_process[i]->restoration_weight->mass);
          if (buff > width) {
            width = buff;
          }
        }
        width++;

        if (discrete_parametric_process[i]->weight) {
          os << "\n";
          if (file_flag) {
            os << "# ";
          }
          os << STAT_label[STATL_THEORETICAL] << " " << SEQ_label[SEQL_STATE_PROBABILITY] << ": ";

          for (j = 0;j < nb_state;j++) {
            os << setw(width) << discrete_parametric_process[i]->weight->mass[j];
          }
          os << endl;
        }

        if (discrete_parametric_process[i]->restoration_weight) {
          if (file_flag) {
            os << "# ";
          }
          os << STAT_label[STATL_RESTORATION] << " " << SEQ_label[SEQL_STATE_PROBABILITY] << ": ";

          for (j = 0;j < nb_state;j++) {
            os << setw(width) << discrete_parametric_process[i]->restoration_weight->mass[j];
          }
          os << endl;
        }

        break;
      }

      else if (continuous_parametric_process[i]) {
        if (continuous_parametric_process[i]->weight) {
          width = column_width(nb_state , continuous_parametric_process[i]->weight->mass);
        }
        else {
          width = 0;
        }
        if (continuous_parametric_process[i]->restoration_weight) {
          buff = column_width(nb_state , continuous_parametric_process[i]->restoration_weight->mass);
          if (buff > width) {
            width = buff;
          }
        }
        width++;

        if (continuous_parametric_process[i]->weight) {
          os << "\n";
          if (file_flag) {
            os << "# ";
          }
          os << STAT_label[STATL_THEORETICAL] << " " << SEQ_label[SEQL_STATE_PROBABILITY] << ": ";

          for (j = 0;j < nb_state;j++) {
            os << setw(width) << continuous_parametric_process[i]->weight->mass[j];
          }
          os << endl;
        }

        if (continuous_parametric_process[i]->restoration_weight) {
          if (file_flag) {
            os << "# ";
          }
          os << STAT_label[STATL_RESTORATION] << " " << SEQ_label[SEQL_STATE_PROBABILITY] << ": ";

          for (j = 0;j < nb_state;j++) {
            os << setw(width) << continuous_parametric_process[i]->restoration_weight->mass[j];
          }
          os << endl;
        }

        break;
      }
    }

    os << "\n" << nb_output_process << " "
       << STAT_word[nb_output_process == 1 ? STATW_OUTPUT_PROCESS : STATW_OUTPUT_PROCESSES] << endl;
  }

  // writing of the distributions associated with each observation process

  if (hidden) {
    distance = new double*[nb_state];
    for (i = 0;i < nb_state;i++) {
      distance[i] = new double[nb_state];
    }
  }

  for (i = 0;i < nb_output_process;i++) {
    os << "\n" << STAT_word[STATW_OUTPUT_PROCESS];

    if (hidden) {
      os << " " << i + 1;

      if (categorical_process[i]) {
        os << " : " << STAT_word[STATW_CATEGORICAL];
      }
      else if (discrete_parametric_process[i]) {
        os << " : " << STAT_word[STATW_DISCRETE_PARAMETRIC];
      }
      else {
        os << " : " << STAT_word[STATW_CONTINUOUS_PARAMETRIC];
      }
    }
    os << endl;

    if ((continuous_parametric_process[i]) && ((continuous_parametric_process[i]->ident == LINEAR_MODEL) ||
         (continuous_parametric_process[i]->ident == AUTOREGRESSIVE_MODEL))) {
      for (j = 0;j < nb_state;j++) {
        os << "\n" << STAT_word[STATW_STATE] << " " << j << " "
           << STAT_word[STATW_OBSERVATION_MODEL] << endl;
        continuous_parametric_process[i]->observation[j]->ascii_parameter_print(os , file_flag);
      }
    }

    if (seq) {
      switch (seq->type[0]) {
      case STATE :
        variable = i + 1;
        break;
      default :
        variable = i;
        break;
      }

      if ((continuous_parametric_process[i]) && (continuous_parametric_process[i]->ident == AUTOREGRESSIVE_MODEL)) {
        if (seq->type[0] == STATE) {
          seq->autoregressive_model_ascii_print(os , variable , continuous_parametric_process[i] , file_flag);
        }
      }

      else if ((categorical_process[i]) || (discrete_parametric_process[i]) ||
               ((continuous_parametric_process[i]) && (continuous_parametric_process[i]->ident != LINEAR_MODEL))) {
        if (seq->observation_distribution) {
          observation_dist = seq->observation_distribution[variable];
        }
        marginal_dist = seq->marginal_distribution[variable];

        if (seq->observation_histogram) {
          observation_histo = seq->observation_histogram[variable];
        }
        marginal_histo = seq->marginal_histogram[variable];

        characteristics = seq->characteristics[variable];
      }
    }

    if (categorical_process[i]) {
      categorical_process[i]->ascii_print(os , i + 1 , observation_dist , marginal_dist ,
                                          characteristics , exhaustive , file_flag);

      if (hidden) {
        for (j = 0;j < nb_state;j++) {
          distance[j][j] = 0.;

          for (k = j + 1;k < nb_state;k++) {
            if ((transition[j][k] > MIN_PROBABILITY) || (transition[k][j] > MIN_PROBABILITY)) {
              distance[j][k] = categorical_process[i]->observation[j]->overlap_distance_computation(*(categorical_process[i]->observation[k]));
            }
            else {
              distance[j][k] = 1.;
            }

            distance[k][j] = distance[j][k];
          }
        }
      }
    }

    else if (discrete_parametric_process[i]) {
      discrete_parametric_process[i]->ascii_print(os , observation_dist , marginal_dist ,
                                                  exhaustive , file_flag);

      if (hidden) {
        for (j = 0;j < nb_state;j++) {
          distance[j][j] = 0.;

          for (k = j + 1;k < nb_state;k++) {
            if ((transition[j][k] > MIN_PROBABILITY) || (transition[k][j] > MIN_PROBABILITY)) {
              distance[j][k] = discrete_parametric_process[i]->observation[j]->sup_norm_distance_computation(*(discrete_parametric_process[i]->observation[k]));
            }
            else {
              distance[j][k] = 1.;
            }

            distance[k][j] = distance[j][k];
          }
        }
      }
    }

    else if ((continuous_parametric_process[i]->ident != LINEAR_MODEL) &&
             (continuous_parametric_process[i]->ident != AUTOREGRESSIVE_MODEL)) {
      continuous_parametric_process[i]->ascii_print(os , observation_histo , observation_dist ,
                                                    marginal_histo , marginal_dist ,
                                                    exhaustive , file_flag);

      if (hidden) {
        for (j = 0;j < nb_state;j++) {
          distance[j][j] = 0.;

          for (k = j + 1;k < nb_state;k++) {
            if ((transition[j][k] > MIN_PROBABILITY) || (transition[k][j] > MIN_PROBABILITY)) {
              distance[j][k] = continuous_parametric_process[i]->observation[j]->sup_norm_distance_computation(*(continuous_parametric_process[i]->observation[k]));
            }
            else {
              distance[j][k] = 1.;
            }

            distance[k][j] = distance[j][k];
          }
        }
      }
    }

    if ((hidden) && ((categorical_process[i]) || (discrete_parametric_process[i]) ||
         ((continuous_parametric_process[i]) && (continuous_parametric_process[i]->ident != LINEAR_MODEL) &&
          (continuous_parametric_process[i]->ident != AUTOREGRESSIVE_MODEL)))) {
      width = column_width(nb_state , distance[0]);
      for (j = 1;j < nb_state;j++) {
        buff = column_width(nb_state , distance[j]);
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
      os << STAT_label[STATL_CONSECUTIVE_STATE_OBSERVATION_DISTRIBUTION_DISTANCE] << endl;

      for (j = 0;j < nb_state;j++) {
        if (file_flag) {
          os << "# ";
        }
        for (k = 0;k < nb_state;k++) {
          if ((k != j) && (transition[j][k] > MIN_PROBABILITY)) {
            os << setw(width) << distance[j][k];
          }
          else {
            os << setw(width) << "_";
          }
        }
        os << endl;
      }
    }
  }

  if (hidden) {
    for (i = 0;i < nb_state;i++) {
      delete [] distance[i];
    }
    delete [] distance;
  }

  if (seq) {
    int nb_parameter = nb_parameter_computation(hidden ? MIN_PROBABILITY : 0.);
    double information;


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

    // writing of the information quantity of the observed sequences in the i.i.d. case

    for (i = 0;i < seq->nb_variable;i++) {
      if (seq->type[i] == REAL_VALUE) {
        break;
      }
    }

    if (i == seq->nb_variable) {
      information = seq->iid_information_computation();

      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << SEQ_label[SEQL_IID_INFORMATION] << ": " << information << " ("
         << information / seq->cumul_length << ")" << endl;
    }

    // writing of the (penalized) log-likelihoods of the model for sequences

    if (hidden) {
      if (seq->restoration_likelihood != D_INF) {
        os << "\n";
        if (file_flag) {
          os << "# ";
        }
        os << SEQ_label[SEQL_STATE_SEQUENCES_LIKELIHOOD] << ": " << seq->restoration_likelihood << "   ("
           << STAT_label[STATL_NORMALIZED] << ": " << seq->restoration_likelihood / seq->cumul_length << ")" << endl;
      }

      if (seq->sample_entropy != D_DEFAULT) {
        os << "\n";
        if (file_flag) {
          os << "# ";
        }
        os << SEQ_label[SEQL_STATE_SEQUENCE_ENTROPY] << ": " << seq->sample_entropy << "   ("
           << STAT_label[STATL_NORMALIZED] << ": " << seq->sample_entropy / seq->cumul_length << ")" << endl;
      }

      if (seq->likelihood != D_INF) {
        os << "\n";
        if (file_flag) {
          os << "# ";
        }
        os << SEQ_label[SEQL_OBSERVED_SEQUENCES_LIKELIHOOD] << ": " << seq->likelihood << "   ("
           << STAT_label[STATL_NORMALIZED] << ": " << seq->likelihood / seq->cumul_length << ")" << endl;
      }
    }

    else {
      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << STAT_label[STATL_LIKELIHOOD] << ": " << seq->likelihood << "   ("
         << STAT_label[STATL_NORMALIZED] << ": " << seq->likelihood / seq->cumul_length << ")" << endl;
    }

    if (seq->likelihood != D_INF) {
      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << nb_parameter << " " << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS]
         << "   2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[AIC] << "): "
         << 2 * (seq->likelihood - nb_parameter) << endl;

      if (nb_parameter < seq->cumul_length - 1) {
        os << "\n";
        if (file_flag) {
          os << "# ";
        }
        os << nb_parameter << " " << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS]
           << "   2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[AICc] << "): "
           << 2 * (seq->likelihood - (double)(nb_parameter * seq->cumul_length) /
             (double)(seq->cumul_length - nb_parameter - 1)) << endl;
      }

      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << nb_parameter << " " << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS]
         << "   2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[BIC] << "): "
         << 2 * seq->likelihood - nb_parameter * log((double)seq->cumul_length) << endl;

      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << nb_parameter << " " << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS]
         << "   2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[BICc] << "): "
         << 2 * seq->likelihood - penalty_computation(hidden , (hidden ? MIN_PROBABILITY : 0.)) << endl;
    }

//    if ((hidden) && (seq->restoration_likelihood != D_INF)) {
    if ((hidden) && (seq->likelihood != D_INF)) {
      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << nb_parameter << " " << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS]
         << "   2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[ICL] << "): "
//         << 2 * seq->restoration_likelihood - nb_parameter * log((double)seq->cumul_length) << endl;
         << 2 * (seq->likelihood - seq->sample_entropy) - nb_parameter * log((double)seq->cumul_length) << endl;

      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << nb_parameter << " " << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS]
         << "   2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[ICLc] << "): "
//         << 2 * seq->restoration_likelihood - penalty_computation(hidden , MIN_PROBABILITY) << endl;
         << 2 * (seq->likelihood - seq->sample_entropy) - penalty_computation(hidden , MIN_PROBABILITY) << endl;
    }
  }

  os.setf(format_flags , ios::adjustfield);

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a SemiMarkov object.
 *
 *  \param[in,out] os         stream,
 *  \param[in]     exhaustive flag detail level.
 */
/*--------------------------------------------------------------*/

ostream& SemiMarkov::ascii_write(ostream &os , bool exhaustive) const

{
  return ascii_write(os , semi_markov_data , exhaustive , false);
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a SemiMarkov object in a file.
 *
 *  \param[in] error      reference on a StatError object,
 *  \param[in] path       file path,
 *  \param[in] exhaustive flag detail level.
 *
 *  \return               error status.
 */
/*--------------------------------------------------------------*/

bool SemiMarkov::ascii_write(StatError &error , const string path ,
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
    ascii_write(out_file , semi_markov_data , exhaustive , true);
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a SemiMarkov object and the associated data structure
 *         in a file at the spreadsheet format.
 *
 *  \param[in,out] os     stream,
 *  \param[in]     seq    pointer on a SemiMarkovData object,
 *  \param[in]     hidden flag hidden model.
 */
/*--------------------------------------------------------------*/

ostream& SemiMarkov::spreadsheet_write(ostream &os , const SemiMarkovData *seq ,
                                       bool hidden) const

{
  int i , j , k;
  int variable;
  double **distance;
  FrequencyDistribution *marginal_dist = NULL , **observation_dist = NULL;
  Histogram *marginal_histo = NULL , **observation_histo = NULL;
  SequenceCharacteristics *characteristics = NULL;


  if (hidden) {
    switch (type) {
    case ORDINARY :
      os << SEQ_word[SEQW_HIDDEN_SEMI_MARKOV_CHAIN] << endl;
      break;
    case EQUILIBRIUM :
      os << SEQ_word[SEQW_EQUILIBRIUM_HIDDEN_SEMI_MARKOV_CHAIN] << endl;
      break;
    }
  }

  else {
    switch (type) {
    case ORDINARY :
      os << SEQ_word[SEQW_SEMI_MARKOV_CHAIN] << endl;
      break;
    case EQUILIBRIUM :
      os << SEQ_word[SEQW_EQUILIBRIUM_SEMI_MARKOV_CHAIN] << endl;
      break;
    }
  }

  // writing of the Markov chain parameters

  spreadsheet_print(os);

  // writing of the state occupancy distributions

  if ((seq) && (seq->type[0] == STATE)) {
    characteristics = seq->characteristics[0];
  }
  else {
    characteristics = NULL;
  }

  state_process->spreadsheet_print(os , 0 , NULL , NULL , characteristics , forward);

  // writing of the distributions associated with each observation process

  if (hidden) {
    os << "\n" << nb_output_process << "\t"
       << STAT_word[nb_output_process == 1 ? STATW_OUTPUT_PROCESS : STATW_OUTPUT_PROCESSES] << endl;
  }

  if (hidden) {
    distance = new double*[nb_state];
    for (i = 0;i < nb_state;i++) {
      distance[i] = new double[nb_state];
    }
  }

  for (i = 0;i < nb_output_process;i++) {
    os << "\n" << STAT_word[STATW_OUTPUT_PROCESS];

    if (hidden) {
      os << "\t" << i + 1;

      if (categorical_process[i]) {
        os << "\t" << STAT_word[STATW_CATEGORICAL];
      }
      else if (discrete_parametric_process[i]) {
        os << "\t" << STAT_word[STATW_DISCRETE_PARAMETRIC];
      }
      else {
        os << "\t" << STAT_word[STATW_CONTINUOUS_PARAMETRIC];
      }
    }
    os << endl;

    if ((continuous_parametric_process[i]) && ((continuous_parametric_process[i]->ident == LINEAR_MODEL) ||
         (continuous_parametric_process[i]->ident == AUTOREGRESSIVE_MODEL))) {
      for (j = 0;j < nb_state;j++) {
        os << "\n" << STAT_word[STATW_STATE] << " " << j << "\t"
           << STAT_word[STATW_OBSERVATION_MODEL] << endl;
        continuous_parametric_process[i]->observation[j]->spreadsheet_parameter_print(os);
      }
    }

    if (seq) {
      switch (seq->type[0]) {
      case STATE :
        variable = i + 1;
        break;
      default :
        variable = i;
        break;
      }

      if ((continuous_parametric_process[i]) && (continuous_parametric_process[i]->ident == LINEAR_MODEL)) {
        seq->linear_model_spreadsheet_print(os , variable , continuous_parametric_process[i]);
      }
      else if ((continuous_parametric_process[i]) && (continuous_parametric_process[i]->ident == AUTOREGRESSIVE_MODEL)) {
        if (seq->type[0] == STATE) {
          seq->autoregressive_model_spreadsheet_print(os , variable , continuous_parametric_process[i]);
        }
      }

      else {
        if (seq->observation_distribution) {
          observation_dist = seq->observation_distribution[variable];
        }
        marginal_dist = seq->marginal_distribution[variable];

        if (seq->observation_histogram) {
          observation_histo = seq->observation_histogram[variable];
        }
        marginal_histo = seq->marginal_histogram[variable];

        characteristics = seq->characteristics[variable];
      }
    }

    if (categorical_process[i]) {
      categorical_process[i]->spreadsheet_print(os , i + 1 , observation_dist , marginal_dist ,
                                                characteristics);

      if (hidden) {
        for (j = 0;j < nb_state;j++) {
          for (k = j + 1;k < nb_state;k++) {
            if ((transition[j][k] > MIN_PROBABILITY) || (transition[k][j] > MIN_PROBABILITY)) {
              distance[j][k] = categorical_process[i]->observation[j]->overlap_distance_computation(*(categorical_process[i]->observation[k]));
              distance[k][j] = distance[j][k];
            }
          }
        }
      }
    }

    else if (discrete_parametric_process[i]) {
      discrete_parametric_process[i]->spreadsheet_print(os , observation_dist , marginal_dist);

      if (hidden) {
        for (j = 0;j < nb_state;j++) {
          for (k = j + 1;k < nb_state;k++) {
            if ((transition[j][k] > MIN_PROBABILITY) || (transition[k][j] > MIN_PROBABILITY)) {
              distance[j][k] = discrete_parametric_process[i]->observation[j]->sup_norm_distance_computation(*(discrete_parametric_process[i]->observation[k]));
              distance[k][j] = distance[j][k];
            }
          }
        }
      }
    }

    else if ((continuous_parametric_process[i]->ident != LINEAR_MODEL) &&
             (continuous_parametric_process[i]->ident != AUTOREGRESSIVE_MODEL)) {
      continuous_parametric_process[i]->spreadsheet_print(os , observation_histo , observation_dist ,
                                                          marginal_histo , marginal_dist);

      if (hidden) {
        for (j = 0;j < nb_state;j++) {
          for (k = j + 1;k < nb_state;k++) {
            if ((transition[j][k] > MIN_PROBABILITY) || (transition[k][j] > MIN_PROBABILITY)) {
              distance[j][k] = continuous_parametric_process[i]->observation[j]->sup_norm_distance_computation(*(continuous_parametric_process[i]->observation[k]));
              distance[k][j] = distance[j][k];
            }
          }
        }
      }
    }

    if ((hidden) && ((categorical_process[i]) || (discrete_parametric_process[i]) ||
         ((continuous_parametric_process[i]) && (continuous_parametric_process[i]->ident != LINEAR_MODEL) &&
          (continuous_parametric_process[i]->ident != AUTOREGRESSIVE_MODEL)))) {
      os << "\n" << STAT_label[STATL_CONSECUTIVE_STATE_OBSERVATION_DISTRIBUTION_DISTANCE] << endl;

      for (j = 0;j < nb_state;j++) {
        for (k = 0;k < nb_state;k++) {
          if ((k != j) && (transition[j][k] > MIN_PROBABILITY)) {
            os << distance[j][k];
          }
          os << "\t";
        }
        os << endl;
      }
    }
  }

  if (hidden) {
    for (i = 0;i < nb_state;i++) {
      delete [] distance[i];
    }
    delete [] distance;
  }

  if (seq) {
    int nb_parameter = nb_parameter_computation(hidden ? MIN_PROBABILITY : 0.);
    double information;


    // writing of the sequence length frequency distribution

    os << "\n" << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\t";
    seq->length_distribution->spreadsheet_characteristic_print(os);

    os << "\n\t" << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << endl;
    seq->length_distribution->spreadsheet_print(os);

    os << "\n" << SEQ_label[SEQL_CUMUL_LENGTH] << "\t" << seq->cumul_length << endl;

    // writing of the information quantity of the observed sequences in the i.i.d. case

    for (i = 0;i < seq->nb_variable;i++) {
      if (seq->type[i] == REAL_VALUE) {
        break;
      }
    }

    if (i == seq->nb_variable) {
      information = seq->iid_information_computation();

      os << "\n" << SEQ_label[SEQL_IID_INFORMATION] << "\t" << information << "\t"
         << information / seq->cumul_length << endl;
    }

    // writing of the (penalized) log-likelihoods of the model for sequences

    if (hidden) {
      if (seq->restoration_likelihood != D_INF) {
        os << "\n" << SEQ_label[SEQL_STATE_SEQUENCES_LIKELIHOOD] << "\t" << seq->restoration_likelihood << "\t"
           << STAT_label[STATL_NORMALIZED] << "\t" << seq->restoration_likelihood / seq->cumul_length << endl;
      }

      if (seq->sample_entropy != D_DEFAULT) {
        os << "\n" << SEQ_label[SEQL_STATE_SEQUENCE_ENTROPY] << "\t" << seq->sample_entropy << "\t"
           << STAT_label[STATL_NORMALIZED] << "\t" << seq->sample_entropy / seq->cumul_length << endl;
      }

      if (seq->likelihood != D_INF) {
        os << "\n" << SEQ_label[SEQL_OBSERVED_SEQUENCES_LIKELIHOOD] << "\t" << seq->likelihood << "\t"
           << STAT_label[STATL_NORMALIZED] << "\t" << seq->likelihood / seq->cumul_length << endl;
      }
    }

    else {
      os << "\n" << STAT_label[STATL_LIKELIHOOD] << "\t" << seq->likelihood << "\t"
         << STAT_label[STATL_NORMALIZED] << "\t" << seq->likelihood / seq->cumul_length << endl;
    }

    if (seq->likelihood != D_INF) {
      os << "\n" << nb_parameter << "\t" << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS] << "\t"
         << "2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[AIC] << ")\t"
         << 2 * (seq->likelihood - nb_parameter) << endl;

      if (nb_parameter < seq->cumul_length - 1) {
        os << "\n" << nb_parameter << "\t" << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS] << "\t"
           << "2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[AICc] << ")\t"
           << 2 * (seq->likelihood - (double)(nb_parameter * seq->cumul_length) /
              (double)(seq->cumul_length - nb_parameter - 1)) << endl;
      }

      os << "\n" << nb_parameter << "\t" << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS] << "\t"
         << "2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[BIC] << ")\t"
         << 2 * seq->likelihood - nb_parameter * log((double)seq->cumul_length) << endl;

      os << "\n" << nb_parameter << "\t" << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS] << "\t"
         << "2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[BICc] << ")\t"
         << 2 * seq->likelihood - penalty_computation(hidden , (hidden ? MIN_PROBABILITY : 0.)) << endl;
    }

//    if ((hidden) && (seq->restoration_likelihood != D_INF)) {
    if ((hidden) && (seq->likelihood != D_INF)) {
      os << "\n" << nb_parameter << "\t" << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS] << "\t"
         << "2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[ICL] << ")\t"
//         << 2 * seq->restoration_likelihood - nb_parameter * log((double)seq->cumul_length) << endl;
         << 2 * (seq->likelihood - seq->sample_entropy) - nb_parameter * log((double)seq->cumul_length) << endl;

      os << "\n" << nb_parameter << "\t" << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS] << "\t"
         << "2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[ICLc] << ")\t"
//         << 2 * seq->restoration_likelihood - penalty_computation(hidden , MIN_PROBABILITY) << endl;
         << 2 * (seq->likelihood - seq->sample_entropy) - penalty_computation(hidden , MIN_PROBABILITY) << endl;
    }
  }

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a SemiMarkov object in a file at the spreadsheet format.
 *
 *  \param[in] error reference on a StatError object,
 *  \param[in] path  file path.
 *
 *  \return          error status.
 */
/*--------------------------------------------------------------*/

bool SemiMarkov::spreadsheet_write(StatError &error , const string path) const

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
    spreadsheet_write(out_file , semi_markov_data);
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Plot of a SemiMarkov object and the associated data structure using Gnuplot.
 *
 *  \param[in] prefix file prefix,
 *  \param[in] title  figure title,
 *  \param[in] seq    pointer on a SemiMarkovData object.
 *
 *  \return           error status.
 */
/*--------------------------------------------------------------*/

bool SemiMarkov::plot_write(const char *prefix , const char *title ,
                            const SemiMarkovData *seq) const

{
  bool status;
  int i;
  int variable , nb_value = I_DEFAULT;
  double *empirical_cdf[2];
  FrequencyDistribution *length_distribution = NULL , *marginal_dist = NULL , **observation_dist = NULL;
  Histogram *marginal_histo = NULL , **observation_histo = NULL;
  SequenceCharacteristics *characteristics = NULL;


  if ((seq) && (seq->type[0] == STATE)) {
    characteristics = seq->characteristics[0];
    length_distribution = seq->length_distribution;
  }
  else {
    characteristics = NULL;
  }

  status = state_process->plot_print(prefix , title , 0 , NULL , NULL ,
                                     characteristics , length_distribution , forward);

  if (status) {
    if (seq) {
      length_distribution = seq->length_distribution;
    }

    for (i = 0;i < nb_output_process;i++) {
      if (seq) {
        switch (seq->type[0]) {
        case STATE :
          variable = i + 1;
          break;
        default :
          variable = i;
          break;
        }

        if ((continuous_parametric_process[i]) && (continuous_parametric_process[i]->ident == LINEAR_MODEL)) {
          seq->linear_model_plot_print(prefix , title , variable , continuous_parametric_process[i]);
        }
        else if ((continuous_parametric_process[i]) && (continuous_parametric_process[i]->ident == AUTOREGRESSIVE_MODEL)) {
          if (seq->type[0] == STATE) {
            seq->autoregressive_model_plot_print(prefix , title , variable , continuous_parametric_process[i]);
          }
        }

        else {
          if (seq->observation_distribution) {
            observation_dist = seq->observation_distribution[variable];
          }
          marginal_dist = seq->marginal_distribution[variable];

          if (seq->observation_histogram) {
            observation_histo = seq->observation_histogram[variable];
          }
          marginal_histo = seq->marginal_histogram[variable];

          characteristics = seq->characteristics[variable];

          if (continuous_parametric_process[i]) {
            nb_value = seq->cumulative_distribution_function_computation(variable , empirical_cdf);
          }
        }
      }

      if (categorical_process[i]) {
        categorical_process[i]->plot_print(prefix , title , i + 1 , observation_dist ,
                                           marginal_dist , characteristics ,
                                           length_distribution);
      }
      else if (discrete_parametric_process[i]) {
        discrete_parametric_process[i]->plot_print(prefix , title , i + 1 , observation_dist ,
                                                   marginal_dist);
      }
      else if ((continuous_parametric_process[i]->ident != LINEAR_MODEL) &&
               (continuous_parametric_process[i]->ident != AUTOREGRESSIVE_MODEL)) {
        continuous_parametric_process[i]->plot_print(prefix , title , i + 1 ,
                                                     observation_histo , observation_dist ,
                                                     marginal_histo , marginal_dist ,
                                                     nb_value , (seq ? empirical_cdf : NULL));
        if (seq) {
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
 *  \brief Plot of a SemiMarkov object using Gnuplot.
 *
 *  \param[in] error  reference on a StatError object,
 *  \param[in] prefix file prefix,
 *  \param[in] title  figure title.
 *
 *  \return           error status.
 */
/*--------------------------------------------------------------*/

bool SemiMarkov::plot_write(StatError &error , const char *prefix ,
                            const char *title) const

{
  bool status = plot_write(prefix , title , semi_markov_data);

  error.init();

  if (!status) {
    error.update(STAT_error[STATR_FILE_PREFIX]);
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Plot of a SemiMarkov object and the associated data structure.
 *
 *  \param[in] seq pointer on a SemiMarkovData object.
 *
 *  \return        MultiPlotSet object.
 */
/*--------------------------------------------------------------*/

MultiPlotSet* SemiMarkov::get_plotable(const SemiMarkovData *seq) const

{
  int i , j;
  int nb_plot_set , index_length , index , variable;
  FrequencyDistribution *length_distribution = NULL , *marginal_dist = NULL , **observation_dist = NULL;
  Histogram *marginal_histo = NULL , **observation_histo = NULL;
  SequenceCharacteristics *characteristics = NULL;
  MultiPlotSet *plot_set;


  if ((seq) && (seq->type[0] == STATE)) {
    characteristics = seq->characteristics[0];
  }
  else {
    characteristics = NULL;
  }

  // computation of the number of plots

  nb_plot_set = 0;

  if ((state_process->index_value) || (characteristics)) {
    nb_plot_set++;

    if (characteristics) {
      index_length = characteristics->index_value->plot_length_computation();

      if (characteristics->index_value->frequency[index_length - 1] < MAX_FREQUENCY) {
        nb_plot_set++;
      }
      nb_plot_set++;
    }
  }

  if ((state_process->first_occurrence) || (characteristics)) {
    for (i = 0;i < nb_state;i++) {
      if ((state_process->first_occurrence) &&
          (state_process->first_occurrence[i])) {
        nb_plot_set++;
      }
      else if ((characteristics) && (i < characteristics->nb_value) &&
               (characteristics->first_occurrence[i]->nb_element > 0)) {
        nb_plot_set++;
      }
    }
  }

  if ((state_process->recurrence_time) || (characteristics)) {
    for (i = 0;i < nb_state;i++) {
      if ((state_process->recurrence_time) &&
          (state_process->recurrence_time[i])) {
        nb_plot_set++;
      }
      else if ((characteristics) && (i < characteristics->nb_value) &&
               (characteristics->recurrence_time[i]->nb_element > 0)) {
        nb_plot_set++;
      }
    }
  }

  if ((state_process->sojourn_time) || (characteristics)) {
    for (i = 0;i < nb_state;i++) {
      if ((state_process->sojourn_time) &&
          (state_process->sojourn_time[i])) {
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

      if ((forward) && (forward[i])) {
        nb_plot_set++;
      }
      else if ((characteristics) && (i < characteristics->nb_value) &&
               (characteristics->final_run[i]->nb_element > 0)) {
        nb_plot_set++;
      }
    }
  }

  if ((state_process->nb_run) || (state_process->nb_occurrence) ||
      ((characteristics) && (characteristics->nb_run) && (characteristics->nb_occurrence))) {
    for (i = 0;i < nb_state;i++) {
      if (state_process->nb_run) {
        nb_plot_set++;
      }
      else if ((characteristics) && (i < characteristics->nb_value) &&
               (characteristics->nb_run) && (characteristics->nb_run[i]->nb_element > 0)) {
        nb_plot_set++;
      }

      if (state_process->nb_occurrence) {
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

  for (i = 0;i < nb_output_process;i++) {
    if (seq) {
      switch (seq->type[0]) {
      case STATE :
        variable = i + 1;
        break;
      default :
        variable = i;
        break;
      }

      characteristics = seq->characteristics[variable];
    }

    if (categorical_process[i]) {
      if ((categorical_process[i]->index_value) || (characteristics)) {
        nb_plot_set++;

        if (characteristics) {
          index_length = characteristics->index_value->plot_length_computation();

          if (characteristics->index_value->frequency[index_length - 1] < MAX_FREQUENCY) {
            nb_plot_set++;
          }
          nb_plot_set++;
        }
      }

      if ((categorical_process[i]->first_occurrence) || (characteristics)) {
        for (j = 0;j < categorical_process[i]->nb_value;j++) {
          if ((categorical_process[i]->first_occurrence) &&
              (categorical_process[i]->first_occurrence[j])) {
            nb_plot_set++;
          }
          else if ((characteristics) && (j < characteristics->nb_value) &&
                   (characteristics->first_occurrence[j]->nb_element > 0)) {
            nb_plot_set++;
          }
        }
      }

      if ((categorical_process[i]->recurrence_time) || (characteristics)) {
        for (j = 0;j < categorical_process[i]->nb_value;j++) {
          if ((categorical_process[i]->recurrence_time) &&
              (categorical_process[i]->recurrence_time[j])) {
            nb_plot_set++;
          }
          else if ((characteristics) && (i < characteristics->nb_value) &&
                   (characteristics->recurrence_time[j]->nb_element > 0)) {
            nb_plot_set++;
          }
        }
      }

      if ((categorical_process[i]->sojourn_time) || (characteristics)) {
        for (j = 0;j < categorical_process[i]->nb_value;j++) {
          if ((categorical_process[i]->sojourn_time) &&
              (categorical_process[i]->sojourn_time[j])) {
            nb_plot_set++;
          }
          else if ((characteristics) && (i < characteristics->nb_value) &&
                   (characteristics->sojourn_time[j]->nb_element > 0)) {
            nb_plot_set++;
          }

/*          if ((characteristics) && (j < characteristics->nb_value) &&
              (characteristics->initial_run) &&
              (characteristics->initial_run[j]->nb_element > 0)) {
            nb_plot_set++;
          } */

          if ((characteristics) && (j < characteristics->nb_value) &&
              (characteristics->final_run[j]->nb_element > 0)) {
            nb_plot_set++;
          }
        }
      }

      if ((categorical_process[i]->nb_run) || (categorical_process[i]->nb_occurrence) ||
          ((characteristics) && (characteristics->nb_run) && (characteristics->nb_occurrence))) {
        for (j = 0;j < categorical_process[i]->nb_value;j++) {
          if (categorical_process[i]->nb_run) {
            nb_plot_set++;
          }
          else if ((characteristics) && (j < characteristics->nb_value) &&
                   (characteristics->nb_run) && (characteristics->nb_run[j]->nb_element > 0)) {
            nb_plot_set++;
          }

          if (categorical_process[i]->nb_occurrence) {
            nb_plot_set++;
          }
          else if ((characteristics) && (j < characteristics->nb_value) &&
                   (characteristics->nb_occurrence) &&
                   (characteristics->nb_occurrence[j]->nb_element > 0)) {
            nb_plot_set++;
          }
        }

        if ((characteristics) && (characteristics->nb_run) && (characteristics->nb_occurrence)) {
          nb_plot_set++;
        }
      }
    }

    if ((seq->observation_distribution) || (seq->observation_histogram)) {
      nb_plot_set += nb_state;
    }
    else {
      nb_plot_set++;
    }

    if ((categorical_process[i]) && (seq->marginal_distribution[variable])) {
      if ((categorical_process[i]->weight) &&
          (categorical_process[i]->mixture)) {
        nb_plot_set++;
      }
      if ((categorical_process[i]->restoration_weight) &&
          (categorical_process[i]->restoration_mixture)) {
        nb_plot_set++;
      }
    }

    if ((discrete_parametric_process[i]) && (seq->marginal_distribution[variable])) {
      if ((discrete_parametric_process[i]->weight) &&
          (discrete_parametric_process[i]->mixture)) {
        nb_plot_set += 2;
      }
      if ((discrete_parametric_process[i]->restoration_weight) &&
          (discrete_parametric_process[i]->restoration_mixture)) {
        nb_plot_set += 2;
      }
    }

    if ((continuous_parametric_process[i]) && (continuous_parametric_process[i]->ident != LINEAR_MODEL) &&
        (continuous_parametric_process[i]->ident != AUTOREGRESSIVE_MODEL) &&
        ((seq->marginal_histogram[variable]) || (seq->marginal_distribution[variable]))) {
      if (continuous_parametric_process[i]->weight) {
        nb_plot_set += 2;
      }
      if (continuous_parametric_process[i]->restoration_weight) {
        nb_plot_set += 2;
      }
    }

    if ((continuous_parametric_process[i]) && ((continuous_parametric_process[i]->ident == LINEAR_MODEL) ||
         ((continuous_parametric_process[i]->ident == AUTOREGRESSIVE_MODEL) && (seq->type[0] == STATE)))) {
      nb_plot_set += nb_state;
    }
  }

  plot_set = new MultiPlotSet(nb_plot_set , nb_output_process + 1);
  plot_set->border = "15 lw 0";

  if ((seq) && (seq->type[0] == STATE)) {
    characteristics = seq->characteristics[0];
    length_distribution = seq->length_distribution;
  }
  else {
    characteristics = NULL;
  }

  index = 0;
  plot_set->variable_nb_viewpoint[0] = 0;
  state_process->plotable_write(*plot_set , index , 0 , NULL , NULL , characteristics ,
                                length_distribution , forward);

  if (seq) {
    length_distribution = seq->length_distribution;
  }

  for (i = 0;i < nb_output_process;i++) {
    if (seq) {
      switch (seq->type[0]) {
      case STATE :
        variable = i + 1;
        break;
      default :
        variable = i;
        break;
      }

      if ((continuous_parametric_process[i]) && (continuous_parametric_process[i]->ident == LINEAR_MODEL)) {
        seq->linear_model_plotable_write(*plot_set , index , variable , continuous_parametric_process[i]);
      }
      else if ((continuous_parametric_process[i]) && (continuous_parametric_process[i]->ident == AUTOREGRESSIVE_MODEL)) {
        if (seq->type[0] == STATE) {
          seq->autoregressive_model_plotable_write(*plot_set , index , variable , continuous_parametric_process[i]);
        }
      }

      else {
        if (seq->observation_distribution) {
          observation_dist = seq->observation_distribution[variable];
        }
        marginal_dist = seq->marginal_distribution[variable];

        if (seq->observation_histogram) {
          observation_histo = seq->observation_histogram[variable];
        }
        marginal_histo = seq->marginal_histogram[variable];

        characteristics = seq->characteristics[variable];
      }
    }

    if (categorical_process[i]) {
      plot_set->variable_nb_viewpoint[i] = 0;
      categorical_process[i]->plotable_write(*plot_set , index , i + 1 , observation_dist ,
                                             marginal_dist , characteristics ,
                                             length_distribution);
    }
    else if (discrete_parametric_process[i]) {
      discrete_parametric_process[i]->plotable_write(*plot_set , index , i + 1 , observation_dist ,
                                                     marginal_dist);
    }
    else if ((continuous_parametric_process[i]->ident != LINEAR_MODEL) &&
             (continuous_parametric_process[i]->ident != AUTOREGRESSIVE_MODEL)) {
      continuous_parametric_process[i]->plotable_write(*plot_set , index , i + 1 ,
                                                       observation_histo , observation_dist ,
                                                       marginal_histo , marginal_dist);
    }
  }

  return plot_set;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Plot of a SemiMarkov object.
 *
 *  \return MultiPlotSet object.
 */
/*--------------------------------------------------------------*/

MultiPlotSet* SemiMarkov::get_plotable() const

{
  return get_plotable(semi_markov_data);
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the number of parameters of a SemiMarkov object.
 *
 *  \param[in] min_probability minimum probability.
 *
 *  \return                    number of parameters.
 */
/*--------------------------------------------------------------*/

int SemiMarkov::nb_parameter_computation(double min_probability) const

{
  int i;
  int nb_parameter = SemiMarkovChain::nb_parameter_computation(min_probability);


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
 *  \brief Computation of an adaptative penalty.
 *
 *  \param[in] hidden          flag hidden model,
 *  \param[in] min_probability minimum probability.
 *
 *  \return                    adaptative penalty.
 */
/*--------------------------------------------------------------*/

double SemiMarkov::penalty_computation(bool hidden , double min_probability) const

{
  int i , j , k;
  int nb_parameter , sample_size;
  double sum , *memory , *state_marginal;
  double penalty = 0.;


  if (semi_markov_data) {
    if (hidden) {
      memory = memory_computation();

      state_marginal = new double[nb_state];

      switch (type) {

      case ORDINARY : {
        sum = 0.;
        for (i = 0;i < state_process->length->nb_value - 2;i++) {
          sum += (1. - state_process->length->cumul[i + 1]);
        }
        for (i = 0;i < nb_state;i++) {
          memory[i] /= sum;
        }

        for (i = 0;i < nb_state;i++) {
          state_marginal[i] = 0.;
        }
        for (i = 0;i < state_process->length->nb_value - 1;i++) {
          for (j = 0;j < nb_state;j++) {
            state_marginal[j] += state_process->index_value->point[j][i] *
                                 (1. - state_process->length->cumul[i]);
          }
        }

        sum = 0.;
        for (i = 0;i < nb_state;i++) {
          sum += state_marginal[i];
        }
        for (i = 0;i < nb_state;i++) {
          state_marginal[i] /= sum;
        }
        break;
      }

      case EQUILIBRIUM : {
        for (i = 0;i < nb_state;i++) {
          state_marginal[i] = initial[i];
        }
        break;
      }
      }
    }

    for (i = 0;i < nb_state;i++) {
      nb_parameter = 0;
      if (!hidden) {
        sample_size = 0;
      }
      for (j = 0;j < nb_state;j++) {
        if (transition[i][j] > min_probability) {
          nb_parameter++;
          if (!hidden) {
            sample_size += semi_markov_data->chain_data->transition[i][j];
          }
        }
      }

      nb_parameter--;

      if (nb_parameter > 0) {
        if (hidden) {
          if (memory[i] > 0.) {
            penalty += nb_parameter * log(memory[i] * semi_markov_data->cumul_length);
          }
        }

        else {
          if (sample_size > 0) {
            penalty += nb_parameter * log((double)sample_size);
          }
        }
      }
    }

    for (i = 0;i < nb_state;i++) {
      if (sojourn_type[i] == SEMI_MARKOVIAN) {
        nb_parameter = state_process->sojourn_time[i]->nb_parameter_computation();
        if (state_process->sojourn_time[i]->inf_bound == 1) {
          nb_parameter--;
        }

        if (hidden) {
          penalty += nb_parameter * log(state_marginal[i] * semi_markov_data->cumul_length);
        }
        else {
          penalty += nb_parameter *
                     log((double)semi_markov_data->marginal_distribution[0]->frequency[i]);
        }
      }
    }

    for (i = 0;i < nb_output_process;i++) {
      if (categorical_process[i]) {
        for (j = 0;j < nb_state;j++) {
          nb_parameter = 0;
          for (k = 0;k < categorical_process[i]->nb_value;k++) {
            if (categorical_process[i]->observation[j]->mass[k] > min_probability) {
              nb_parameter++;
            }
          }

          nb_parameter--;

          if (nb_parameter > 0) {
            if (hidden) {
              penalty += nb_parameter * log(state_marginal[j] * semi_markov_data->cumul_length);
            }
            else {
              penalty += nb_parameter *
                         log((double)semi_markov_data->marginal_distribution[0]->frequency[j]);
            }
          }
        }
      }

      else if (discrete_parametric_process[i]) {
        for (j = 0;j < nb_state;j++) {
          nb_parameter = discrete_parametric_process[i]->observation[j]->nb_parameter_computation();

          if (hidden) {
            penalty += nb_parameter * log(state_marginal[j] * semi_markov_data->cumul_length);
          }
          else {
            penalty += nb_parameter *
                       log((double)semi_markov_data->marginal_distribution[0]->frequency[j]);
          }
        }
      }

      else {
        for (j = 0;j < nb_state;j++) {
          nb_parameter = continuous_parametric_process[i]->observation[j]->nb_parameter_computation();

          if (hidden) {
            penalty += nb_parameter * log(state_marginal[j] * semi_markov_data->cumul_length);
          }
          else {
            penalty += nb_parameter *
                       log((double)semi_markov_data->marginal_distribution[0]->frequency[j]);
          }
        }
      }
    }

    if (hidden) {
      delete [] memory;
      delete [] state_marginal;
    }
  }

  return penalty;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Default constructor of the SemiMarkovData class.
 */
/*--------------------------------------------------------------*/

SemiMarkovData::SemiMarkovData()

{
  semi_markov = NULL;
  chain_data = NULL;

  likelihood = D_INF;
  restoration_likelihood = D_INF;
  sample_entropy = D_DEFAULT;

  posterior_probability = NULL;
  posterior_state_probability = NULL;
  entropy = NULL;
  nb_state_sequence = NULL;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Constructor of the SemiMarkovData class.
 *
 *  \param[in] ilength_distribution sequence length frequency distribution,
 *  \param[in] inb_variable         number of variables,
 *  \param[in] itype                variable types,
 *  \param[in] init_flag            flag initialization.
 */
/*--------------------------------------------------------------*/

SemiMarkovData::SemiMarkovData(const FrequencyDistribution &ilength_distribution , int inb_variable ,
                               variable_nature *itype , bool init_flag)
:MarkovianSequences(ilength_distribution , inb_variable , itype , init_flag)

{
  semi_markov = NULL;
  chain_data = NULL;

  likelihood = D_INF;
  restoration_likelihood = D_INF;
  sample_entropy = D_DEFAULT;

  posterior_probability = NULL;
  posterior_state_probability = NULL;
  entropy = NULL;
  nb_state_sequence = NULL;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Construction of a SemiMarkovData object from a MarkovianSequences object
 *         adding a state variable.
 *
 *  \param[in] seq reference on a MarkovianSequences object.
 */
/*--------------------------------------------------------------*/

SemiMarkovData::SemiMarkovData(const MarkovianSequences &seq)
:MarkovianSequences(seq , ADD_STATE_VARIABLE , UNCHANGED)

{
  semi_markov = NULL;
  chain_data = NULL;

  likelihood = D_INF;
  restoration_likelihood = D_INF;
  sample_entropy = D_DEFAULT;

  posterior_probability = NULL;
  posterior_state_probability = NULL;
  entropy = NULL;
  nb_state_sequence = NULL;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Construction of a SemiMarkovData object from a MarkovianSequences object.
 *
 *  \param[in] seq              reference on a MarkovianSequences object,
 *  \param[in] transform        type of transform (SEQUENCE_COPY/ADD_STATE_VARIABLE),
 *  \param[in] initial_run_flag addition/removing of the initial run length frequency distributions.
 */
/*--------------------------------------------------------------*/

SemiMarkovData::SemiMarkovData(const MarkovianSequences &seq , sequence_transformation transform ,
                               bool initial_run_flag)
:MarkovianSequences(seq , transform , (initial_run_flag ? ADD_INITIAL_RUN : REMOVE_INITIAL_RUN))

{
  semi_markov = NULL;
  chain_data = NULL;

  likelihood = D_INF;
  restoration_likelihood = D_INF;
  sample_entropy = D_DEFAULT;

  posterior_probability = NULL;
  posterior_state_probability = NULL;
  entropy = NULL;
  nb_state_sequence = NULL;
}


/*--------------------------------------------------------------*/
/**
 *   \brief Copy of a SemiMarkovData object.
 *
 *  \param[in] seq        reference on a SemiMarkovData object,
 *  \param[in] model_flag flag copy of the included SemiMarkov object.
 */
/*--------------------------------------------------------------*/

void SemiMarkovData::copy(const SemiMarkovData &seq , bool model_flag)

{
  int i;


  if ((model_flag) && (seq.semi_markov)) {
    semi_markov = new SemiMarkov(*(seq.semi_markov) , false);
  }
  else {
    semi_markov = NULL;
  }

  if (seq.chain_data) {
    chain_data = new ChainData(*(seq.chain_data));
  }
  else {
    chain_data = NULL;
  }

  likelihood = seq.likelihood;
  restoration_likelihood = seq.restoration_likelihood;
  sample_entropy = seq.sample_entropy;

  if (seq.posterior_probability) {
    posterior_probability = new double[nb_sequence];
    for (i = 0;i < nb_sequence;i++) {
      posterior_probability[i] = seq.posterior_probability[i];
    }
  }
  else {
    posterior_probability = NULL;
  }

  if (seq.posterior_state_probability) {
    posterior_state_probability = new double[nb_sequence];
    for (i = 0;i < nb_sequence;i++) {
      posterior_state_probability[i] = seq.posterior_state_probability[i];
    }
  }
  else {
    posterior_state_probability = NULL;
  }

  if (seq.entropy) {
    entropy = new double[nb_sequence];
    for (i = 0;i < nb_sequence;i++) {
      entropy[i] = seq.entropy[i];
    }
  }
  else {
    entropy = NULL;
  }

  if (seq.nb_state_sequence) {
    nb_state_sequence = new double[nb_sequence];
    for (i = 0;i < nb_sequence;i++) {
      nb_state_sequence[i] = seq.nb_state_sequence[i];
    }
  }
  else {
    nb_state_sequence = NULL;
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Destructor of the SemiMarkovData class.
 */
/*--------------------------------------------------------------*/

SemiMarkovData::~SemiMarkovData()

{
  delete semi_markov;
  delete chain_data;

  delete [] posterior_probability;
  delete [] posterior_state_probability;
  delete [] entropy;
  delete [] nb_state_sequence;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Assignment operator of the SemiMarkovData class.
 *
 *  \param[in] seq reference on a SemiMarkovData object.
 *
 *  \return        SemiMarkovData object.
 */
/*--------------------------------------------------------------*/

SemiMarkovData& SemiMarkovData::operator=(const SemiMarkovData &seq)

{
  if (&seq != this) {
    delete semi_markov;
    delete chain_data;

    delete [] posterior_probability;
    delete [] posterior_state_probability;
    delete [] entropy;
    delete [] nb_state_sequence;

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
 *  \param[in] variable   variable,
 *  \param[in] value      state or observation.
 *
 *  \return               DiscreteDistributionData object.
 */
/*--------------------------------------------------------------*/

DiscreteDistributionData* SemiMarkovData::extract(StatError &error , process_distribution histo_type ,
                                                  int variable , int value) const

{
  bool status = true;
  Distribution *pdist;
  DiscreteParametric *pparam;
  FrequencyDistribution *phisto;
  DiscreteDistributionData *histo;
  CategoricalSequenceProcess *process;


  histo = NULL;
  error.init();

  phisto = NULL;

  if (histo_type == OBSERVATION) {
    if ((variable < 2) || (variable > nb_variable)) {
      status = false;
      error.update(STAT_error[STATR_VARIABLE_INDEX]);
    }

    else {
      variable--;

      if ((value < 0) || (value >= marginal_distribution[0]->nb_value)) {
        status = false;
        ostringstream error_message;
        error_message << STAT_label[STATL_STATE] << " " << value << " "
                      << STAT_error[STATR_NOT_PRESENT];
        error.update((error_message.str()).c_str());
      }

      else {
        if (!observation_distribution[variable]) {
          status = false;
          ostringstream error_message;
          error_message << STAT_label[STATL_VARIABLE] << " " << variable + 1 << ": "
                        << STAT_error[STATR_VARIABLE_TYPE];
          error.correction_update((error_message.str()).c_str() , STAT_variable_word[INT_VALUE]);
        }

        else {
          phisto = observation_distribution[variable][value];

          if (phisto->nb_element == 0) {
            status = false;
            error.update(STAT_error[STATR_EMPTY_SAMPLE]);
          }
        }
      }
    }
  }

  else {
    if ((variable < 1) || (variable > nb_variable)) {
      status = false;
      error.update(STAT_error[STATR_VARIABLE_INDEX]);
    }

    else {
      variable--;

      if (!characteristics[variable]) {
        status = false;
        ostringstream error_message;
        error_message << STAT_label[STATL_VARIABLE] << " " << variable + 1 << ": "
                      << SEQ_error[SEQR_CHARACTERISTICS_NOT_COMPUTED];
        error.update((error_message.str()).c_str());
      }

      else if ((value < 0) || (value >= marginal_distribution[variable]->nb_value)) {
        status = false;
        ostringstream error_message;
        error_message << STAT_label[variable == 0 ? STATL_STATE : STATL_OUTPUT] << " "
                      << value << " " << STAT_error[STATR_NOT_PRESENT];
        error.update((error_message.str()).c_str());
      }

      if (status) {
        switch (histo_type) {

        case FIRST_OCCURRENCE : {
          phisto = characteristics[variable]->first_occurrence[value];
          break;
        }

        case RECURRENCE_TIME : {
          phisto = characteristics[variable]->recurrence_time[value];
          break;
        }

        case SOJOURN_TIME : {
          phisto = characteristics[variable]->sojourn_time[value];
          break;
        }

        case INITIAL_RUN : {
          if (characteristics[variable]->initial_run) {
            phisto = characteristics[variable]->initial_run[value];
          }
          else {
            status = false;
            error.update(STAT_error[STATR_NON_EXISTING_FREQUENCY_DISTRIBUTION]);
          }
          break;
        }

        case FINAL_RUN : {
          phisto = characteristics[variable]->final_run[value];
          break;
        }

        case NB_RUN : {
          phisto = characteristics[variable]->nb_run[value];
          break;
        }

        case NB_OCCURRENCE : {
          phisto = characteristics[variable]->nb_occurrence[value];
          break;
        }
        }

        if ((phisto) && (phisto->nb_element == 0)) {
          status = false;
          error.update(STAT_error[STATR_EMPTY_SAMPLE]);
        }
      }
    }
  }

  if (status) {
    if (variable == 0) {
      process = semi_markov->state_process;
    }
    else {
      process = semi_markov->categorical_process[variable - 1];
    }

    pdist = NULL;
    pparam = NULL;

    switch (histo_type) {

    case OBSERVATION : {
      if (semi_markov->categorical_process[variable - 1]) {
        pdist = semi_markov->categorical_process[variable - 1]->observation[value];
      }
      else if (semi_markov->discrete_parametric_process[variable - 1]) {
        pparam = semi_markov->discrete_parametric_process[variable - 1]->observation[value];
      }
      break;
    }

    case FIRST_OCCURRENCE : {
      pdist = process->first_occurrence[value];
      break;
    }

    case RECURRENCE_TIME : {
      pdist = process->recurrence_time[value];
      break;
    }

    case SOJOURN_TIME : {
      pparam = process->sojourn_time[value];
      break;
    }

    case INITIAL_RUN : {
      if ((variable == 0) && (semi_markov->forward)) {
        pdist = semi_markov->forward[value];
      }
      break;
    }

    case FINAL_RUN : {
      if ((variable == 0) && (semi_markov->forward)) {
        pdist = semi_markov->forward[value];
      }
      break;
    }

    case NB_RUN : {
      pdist = process->nb_run[value];
      break;
    }

    case NB_OCCURRENCE : {
      pdist = process->nb_occurrence[value];
      break;
    }
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
 *  \brief Copy of a SemiMarkovData object transforming the implicit index parameters in
 *         explicit index parameters.
 *
 *  \param[in] error reference on a StatError object.
 *
 *  \return          SemiMarkovData object.
 */
/*--------------------------------------------------------------*/

SemiMarkovData* SemiMarkovData::explicit_index_parameter(StatError &error) const

{
  SemiMarkovData *seq;


  error.init();

  if (index_parameter) {
    seq = NULL;
    error.update(SEQ_error[SEQR_INDEX_PARAMETER_TYPE]);
  }
  else {
    seq = new SemiMarkovData(*this , true , EXPLICIT_INDEX_PARAMETER);
  }

  return seq;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Removing of the index parameters.
 *
 *  \param[in] error reference on a StatError object.
 *
 *  \return          SemiMarkovData object.
 */
/*--------------------------------------------------------------*/

SemiMarkovData* SemiMarkovData::remove_index_parameter(StatError &error) const

{
  SemiMarkovData *seq;


  error.init();

  if (!index_parameter) {
    seq = NULL;
    error.update(SEQ_error[SEQR_INDEX_PARAMETER_TYPE]);
  }
  else {
    seq = new SemiMarkovData(*this , true , REMOVE_INDEX_PARAMETER);
  }

  return seq;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Construction of the auxiliary variables corresponding to
 *         the restored state sequences.
 *
 *  \param[in] error reference on a StatError object.
 *
 *  \return          SemiMarkovData object.
 */
/*--------------------------------------------------------------*/

MarkovianSequences* SemiMarkovData::build_auxiliary_variable(StatError &error) const

{
  bool status = true;
  int i;
  MarkovianSequences *seq;


  seq = NULL;
  error.init();

  if (type[0] != STATE) {
    status = false;
    ostringstream error_message;
    error_message << STAT_label[STATL_VARIABLE] << " 1: "
                  << STAT_error[STATR_VARIABLE_TYPE];
    error.correction_update((error_message.str()).c_str() , STAT_variable_word[STATE]);
  }

  for (i = 0;i < semi_markov->nb_output_process;i++) {
    if (((semi_markov->discrete_parametric_process) && (semi_markov->discrete_parametric_process[i])) ||
        ((semi_markov->continuous_parametric_process) && (semi_markov->continuous_parametric_process[i]))) {
      break;
    }
  }

  if (i == semi_markov->nb_output_process) {
    status = false;
    error.update(SEQ_error[SEQR_PARAMETRIC_PROCESS]);
  }

  if (status) {
    seq = MarkovianSequences::build_auxiliary_variable(semi_markov->discrete_parametric_process ,
                                                       semi_markov->continuous_parametric_process);
  }

  return seq;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Building of residual sequences on the basis of restored state sequences.
 *
 *  \param[in] error reference on a StatError object.
 *
 *  \return          MarkovianSequences object.
 */
/*--------------------------------------------------------------*/

MarkovianSequences* SemiMarkovData::residual_sequences(StatError &error) const

{
  MarkovianSequences *seq;


  error.init();

  if (type[0] != STATE) {
    seq = NULL;

    ostringstream error_message;
    error_message << STAT_label[STATL_VARIABLE] << " 1: "
                  << STAT_error[STATR_VARIABLE_TYPE];
    error.correction_update((error_message.str()).c_str() , STAT_variable_word[STATE]);
  }

  else {
    seq = MarkovianSequences::residual_sequences(semi_markov->categorical_process ,
                                                 semi_markov->discrete_parametric_process ,
                                                 semi_markov->continuous_parametric_process);
  }

  return seq;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a SemiMarkovData object.
 *
 *  \param[in,out] os         stream,
 *  \param[in]     exhaustive flag detail level.
 */
/*--------------------------------------------------------------*/

ostream& SemiMarkovData::ascii_write(ostream &os , bool exhaustive) const

{
  if (semi_markov) {
    semi_markov->ascii_write(os , this , exhaustive , false ,
                             CategoricalSequenceProcess::test_hidden(semi_markov->nb_output_process , semi_markov->categorical_process));
  }

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a SemiMarkovData object in a file.
 *
 *  \param[in] error      reference on a StatError object,
 *  \param[in] path       file path,
 *  \param[in] exhaustive flag detail level.
 *
 *  \return               error status.
 */
/*--------------------------------------------------------------*/

bool SemiMarkovData::ascii_write(StatError &error , const string path ,
                                 bool exhaustive) const

{
  bool status = false;


  if (semi_markov) {
    ofstream out_file(path.c_str());

    error.init();

    if (!out_file) {
      status = false;
      error.update(STAT_error[STATR_FILE_NAME]);
    }

    else {
      status = true;
      semi_markov->ascii_write(out_file , this , exhaustive , true ,
                               CategoricalSequenceProcess::test_hidden(semi_markov->nb_output_process , semi_markov->categorical_process));
    }
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a SemiMarkovData object.
 *
 *  \param[in,out] os         stream,
 *  \param[in]     format     format (line/column),
 *  \param[in]     exhaustive flag detail level.
 *
 *  \return                   error status.
 */
/*--------------------------------------------------------------*/

ostream& SemiMarkovData::ascii_data_write(ostream &os , output_sequence_format format ,
                                          bool exhaustive) const

{
  MarkovianSequences::ascii_write(os , exhaustive , false);
  ascii_print(os , format , false , posterior_probability , entropy , nb_state_sequence , posterior_state_probability);

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a SemiMarkovData object in a file.
 *
 *  \param[in] error      reference on a StatError object,
 *  \param[in] path       file path,
 *  \param[in] format     format (line/column),
 *  \param[in] exhaustive flag detail level.
 *
 *  \return               error status.
 */
/*--------------------------------------------------------------*/

bool SemiMarkovData::ascii_data_write(StatError &error , const string path ,
                                      output_sequence_format format , bool exhaustive) const

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
    if (format != 'a') {
      MarkovianSequences::ascii_write(out_file , exhaustive , true);
    }
    ascii_print(out_file , format , true , posterior_probability , entropy , nb_state_sequence , posterior_state_probability);
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a SemiMarkovData object in a file at the spreadsheet format.
 *
 *  \param[in] error reference on a StatError object,
 *  \param[in] path  file path.
 *
 *  \return          error status.
 */
/*--------------------------------------------------------------*/

bool SemiMarkovData::spreadsheet_write(StatError &error , const string path) const

{
  bool status = false;


  if (semi_markov) {
    ofstream out_file(path.c_str());

    error.init();

    if (!out_file) {
      status = false;
      error.update(STAT_error[STATR_FILE_NAME]);
    }

    else {
      status = true;
      semi_markov->spreadsheet_write(out_file , this ,
                                     CategoricalSequenceProcess::test_hidden(semi_markov->nb_output_process , semi_markov->categorical_process));
    }
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Plot of a SemiMarkovData object using Gnuplot.
 *
 *  \param[in] error  reference on a StatError object,
 *  \param[in] prefix file prefix,
 *  \param[in] title  figure title.
 *
 *  \return           error status.
 */
/*--------------------------------------------------------------*/

bool SemiMarkovData::plot_write(StatError &error , const char *prefix ,
                                const char *title) const

{
  bool status = false;


  if (semi_markov) {
    status = semi_markov->plot_write(prefix , title , this);

    error.init();

    if (!status) {
      error.update(STAT_error[STATR_FILE_PREFIX]);
    }
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Plot of a SemiMarkovData object.
 *
 *  \return MultiPlotSet object.
 */
/*--------------------------------------------------------------*/

MultiPlotSet* SemiMarkovData::get_plotable() const

{
  MultiPlotSet *plot_set;


  if (semi_markov) {
    plot_set = semi_markov->get_plotable(this);
  }
  else {
    plot_set = NULL;
  }

  return plot_set;
}


};  // namespace sequence_analysis
