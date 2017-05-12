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



#include <iomanip>

#include <boost/tokenizer.hpp>

#include "markovian.h"
#include "stat_label.h"

using namespace std;
using namespace boost;


namespace stat_tool {



/*--------------------------------------------------------------*/
/**
 *  \brief Constructor of the CategoricalProcess class.
 *
 *  \param[in] inb_state        number of states,
 *  \param[in] inb_value        number of categories,
 *  \param[in] observation_flag flag on the construction of the observation distributions.
 */
/*--------------------------------------------------------------*/

CategoricalProcess::CategoricalProcess(int inb_state , int inb_value ,
                                       bool observation_flag)

{
  nb_state = inb_state;
  nb_value = inb_value;

  if (observation_flag) {
    int i;


    observation = new Distribution*[nb_state];
    for (i = 0;i < nb_state;i++) {
      observation[i] = new Distribution(nb_value);
    }
  }

  else {
    observation = NULL;
  }

  weight = NULL;
  mixture = NULL;
  restoration_weight = NULL;
  restoration_mixture = NULL;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Constructor of the CategoricalProcess class.
 *
 *  \param[in] inb_state               number of states,
 *  \param[in] inb_value               number of categories,
 *  \param[in] observation_probability pointer on the observation probabilities.
 */
/*--------------------------------------------------------------*/

CategoricalProcess::CategoricalProcess(int inb_state , int inb_value ,
                                       double **observation_probability)

{
  int i , j;


  nb_state = inb_state;
  nb_value = inb_value;

  observation = new Distribution*[nb_state];
  for (i = 0;i < nb_state;i++) {
    observation[i] = new Distribution(nb_value);
    for (j = 0;j < nb_value;j++) {
      observation[i]->mass[j] = observation_probability[i][j];
    }

    observation[i]->cumul_computation();

    observation[i]->max_computation();
//    observation[i]->mean_computation();
//    observation[i]->variance_computation();
  }

  weight = NULL;
  mixture = NULL;
  restoration_weight = NULL;
  restoration_mixture = NULL;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Constructor of the CategoricalProcess class.
 *
 *  \param[in] inb_state    number of states,
 *  \param[in] pobservation pointer on the observation distributions.
 */
/*--------------------------------------------------------------*/

CategoricalProcess::CategoricalProcess(int inb_state , Distribution **pobservation)

{
  int i;

  nb_value = 0;
  nb_state = inb_state;

  observation = new Distribution*[nb_state];
  for (i = 0;i < nb_state;i++) {
    observation[i] = new Distribution(*pobservation[i]);

    observation[i]->cumul_computation();

    observation[i]->max_computation();
    nb_value = max(nb_value, observation[i]->nb_value);
  }

  weight = NULL;
  mixture = NULL;
  restoration_weight = NULL;
  restoration_mixture = NULL;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Copy of a CategoricalProcess object.
 *
 *  \param[in] process reference on a CategoricalProcess object.
 */
/*--------------------------------------------------------------*/

void CategoricalProcess::copy(const CategoricalProcess &process)

{
  int i;


  nb_state = process.nb_state;
  nb_value = process.nb_value;

  if (process.observation) {
    observation = new Distribution*[nb_state];
    for (i = 0;i < nb_state;i++) {
      observation[i] = new Distribution(*(process.observation[i]));
    }
  }

  else {
    observation = NULL;
  }

  if ((process.weight) && (process.mixture)) {
    weight = new Distribution(*(process.weight));
    mixture = new Distribution(*(process.mixture));
  }
  else {
    weight = NULL;
    mixture = NULL;
  }

  if ((process.restoration_weight) && (process.restoration_mixture)) {
    restoration_weight = new Distribution(*(process.restoration_weight));
    restoration_mixture = new Distribution(*(process.restoration_mixture));
  }
  else {
    restoration_weight = NULL;
    restoration_mixture = NULL;
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Destruction of the data members of a CategoricalProcess object.
 */
/*--------------------------------------------------------------*/

void CategoricalProcess::remove()

{
  if (observation) {
    int i;


    for (i = 0;i < nb_state;i++) {
      delete observation[i];
    }
    delete [] observation;

    observation = NULL;
  }

  delete weight;
  weight = NULL;
  delete mixture;
  mixture = NULL;

  delete restoration_weight;
  restoration_weight = NULL;
  delete restoration_mixture;
  restoration_mixture = NULL;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Destructor of the CategoricalProcess class.
 */
/*--------------------------------------------------------------*/

CategoricalProcess::~CategoricalProcess()

{
  remove();
}


/*--------------------------------------------------------------*/
/**
 *  \brief Assignment operator of the CategoricalProcess class.
 *
 *  \param[in] process reference on a CategoricalProcess object.
 *
 *  \return            CategoricalProcess object.
 */
/*--------------------------------------------------------------*/

CategoricalProcess& CategoricalProcess::operator=(const CategoricalProcess &process)

{
  if (&process != this) {
    remove();
    copy(process);
  }

  return *this;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Analysis of the format of categorical observation distributions.
 *
 *  \param[in] error    reference on a StatError object,
 *  \param[in] in_file  stream,
 *  \param[in] line     reference on the file line index,
 *  \param[in] nb_state number of states,
 *  \param[in] model    model type,
 *  \param[in] hidden   flag on the overlap of observation distributions.
 *
 *  \return             CategoricalProcess object.
 */
/*--------------------------------------------------------------*/

CategoricalProcess* CategoricalProcess::parsing(StatError &error , ifstream &in_file , int &line ,
                                                int nb_state , model_type model , bool hidden)

{
  string buffer;
  size_t position;
  typedef tokenizer<char_separator<char>> tokenizer;
  char_separator<char> separator(" \t");
  enum{PSTATE , POUTPUT};
  bool status = true , lstatus , defined_output[NB_OUTPUT];
  int i , j;
  int type , state , output , value , nb_value;
  double proba , cumul , **observation_probability;
  CategoricalProcess *process;


  process = NULL;

  for (i = 0;i < NB_OUTPUT;i++) {
    defined_output[i] = false;
  }

  observation_probability = new double*[nb_state];
  for (i = 0;i < nb_state;i++) {
    observation_probability[i] = new double[NB_OUTPUT];
    for (j = 0;j < NB_OUTPUT;j++) {
      observation_probability[i][j] = 0.;
    }
  }

  state = -1;
  output = NB_OUTPUT;
  cumul = 0.;

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
      switch (i) {

      case 0 : {

        // test COMPONENT/STATE keyword

        if ((*token == STAT_word[STATW_STATE]) || (*token == STAT_word[STATW_COMPONENT])) {
          type = PSTATE;

          switch (model) {

          case MIXTURE : {
            if (*token != STAT_word[STATW_COMPONENT]) {
              status = false;
              error.correction_update(STAT_parsing[STATP_KEYWORD] ,
                                      STAT_word[STATW_COMPONENT] , line , i + 1);
            }
            break;
          }

          case HIDDEN_MARKOV : {
            if (*token != STAT_word[STATW_STATE]) {
              status = false;
              error.correction_update(STAT_parsing[STATP_KEYWORD] ,
                                      STAT_word[STATW_STATE] , line , i + 1);
            }
            break;
          }
          }

          if (state >= 0) {
            if (cumul < 1. - DOUBLE_ERROR) {
              status = false;
              error.update(STAT_parsing[STATP_PROBABILITY_SUM] , line);
            }
          }

          state++;
          output = I_DEFAULT;
          cumul = 0.;
        }

        // test OUTPUT keyword

        else if (*token == STAT_word[STATW_OUTPUT]) {
          type = POUTPUT;
        }

        else {
          status = false;
          error.update(STAT_parsing[STATP_KEYWORD] , line , i + 1);
        }

        break;
      }

      // test component/state index or observation category

      case 1 : {
        lstatus = true;

/*        try {
          value = stoi(*token);   in C++ 11
        }
        catch(invalid_argument &arg) {
          lstatus = false;
        } */
        value = atoi(token->c_str());

        switch (type) {

        case PSTATE : {
          if ((lstatus) && ((value != state) || (value >= nb_state))) {
            lstatus = false;
          }

          if (!lstatus) {
            status = false;

            switch (model) {
            case MIXTURE :
              error.update(STAT_parsing[STATP_COMPONENT_INDEX] , line , i + 1);
              break;
            case HIDDEN_MARKOV :
              error.update(STAT_parsing[STATP_STATE_INDEX] , line , i + 1);
              break;
            }
          }
          break;
        }

        case POUTPUT : {
          if (lstatus) {
            if ((value <= output) || (value >= NB_OUTPUT)) {
              lstatus = false;
            }
            else {
              defined_output[value] = true;
              output = value;
            }
          }

          if (!lstatus) {
            status = false;
            error.update(STAT_parsing[STATP_OUTPUT_INDEX] , line , i + 1);
          }
          break;
        }
        }

        break;
      }

      case 2 : {
        switch (type) {

        // test OBSERVATION_DISTRIBUTION keyword

        case PSTATE : {
          if (*token != STAT_word[STATW_OBSERVATION_DISTRIBUTION]) {
            status = false;
            error.correction_update(STAT_parsing[STATP_KEYWORD] , STAT_word[STATW_OBSERVATION_DISTRIBUTION] , line , i + 1);
          }
          break;
        }

        // test separator

        case POUTPUT : {
          if (*token != ":") {
            status = false;
            error.update(STAT_parsing[STATP_SEPARATOR] , line , i + 1);
          }
          break;
        }
        }

        break;
      }

      // test observation probability value

      case 3 : {
        if (type == POUTPUT) {
          lstatus = true;

/*          try {
            proba = stod(*token);   in C++ 11
          }
          catch (invalid_argument &arg) {
            lstatus = false;
          } */
          proba = atof(token->c_str());

          if (lstatus) {
            if ((proba < 0.) || (proba > 1. - cumul + DOUBLE_ERROR)) {
              lstatus = false;
            }
            else {
              cumul += proba;
              if (status) {
                observation_probability[state][output] = proba;
              }
            }
          }

          if (!lstatus) {
            status = false;
            error.update(STAT_parsing[STATP_OBSERVATION_PROBABILITY] , line , i + 1);
          }
        }
        break;
      }
      }

      i++;
    }

    if (i > 0) {
      if (((type == PSTATE) && (i != 2) && (i != 3)) || ((type == POUTPUT) && (i != 4))) {
        status = false;
        error.update(STAT_parsing[STATP_FORMAT] , line);
      }
    }

    if ((state == nb_state - 1) && (cumul >= 1. - DOUBLE_ERROR)) {
      break;
    }
  }

  if (state != nb_state - 1) {
    status = false;
    error.update(STAT_parsing[STATP_NB_OBSERVATION_DISTRIBUTION]);
  }
  if (cumul < 1. - DOUBLE_ERROR) {
    status = false;
    error.update(STAT_parsing[STATP_PROBABILITY_SUM]);
  }

  i = NB_OUTPUT - 1;
  while (!defined_output[i]) {
    i--;
  }
  nb_value = i + 1;

  // test consecutive categories

  for (i = 0;i < nb_value;i++) {
    if (!defined_output[i]) {
      status = false;
      error.update(STAT_parsing[STATP_NON_CONSECUTIVE_OUTPUTS]);
      break;
    }
  }

  if (status) {
    process = new CategoricalProcess(nb_state , nb_value , observation_probability);

    // test overlap or not of observation distributions

    lstatus = process->test_hidden();

    if ((!hidden) && (lstatus)) {
      status = false;
      error.update(STAT_parsing[STATP_OBSERVATION_DISTRIBUTION_OVERLAP]);
    }

    if ((hidden) && (!lstatus)) {
      status = false;
      error.update(STAT_parsing[STATP_OBSERVATION_DISTRIBUTION_NON_OVERLAP]);
    }

    if (!status) {
      delete process;
      process = NULL;
    }
  }

  for (i = 0;i < nb_state;i++) {
    delete [] observation_probability[i];
  }
  delete [] observation_probability;

  return process;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Analysis of the format of observation distributions
 *         for a collection of categorical observation processes.
 *
 *  \param[in] error             reference on a StatError object,
 *  \param[in] in_file           stream,
 *  \param[in] line              reference on the file line index,
 *  \param[in] nb_state          number of states,
 *  \param[in] nb_output_process number of observation processes.
 *
 *  \return                      categorical observation processes.
 */
/*--------------------------------------------------------------*/

CategoricalProcess** CategoricalProcess::old_parsing(StatError &error , ifstream &in_file , int &line ,
                                                     int nb_state , int &nb_output_process)

{
  string buffer;
  size_t position;
  typedef tokenizer<char_separator<char>> tokenizer;
  char_separator<char> separator(" \t");
  bool status = true , lstatus;
  int i;
  int value , index;
  CategoricalProcess **process;


  nb_output_process = I_DEFAULT;
  process = NULL;

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

      // test OBSERVATION_PROBABILITIES keyword

      if (i == 0) {
        if (*token != STAT_word[STATW_OBSERVATION_PROBABILITIES]) {
          status = false;
          error.correction_update(STAT_parsing[STATP_KEYWORD] , STAT_word[STATW_OBSERVATION_PROBABILITIES] , line);
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
      switch (i) {

      // test number of observation processes

      case 0 : {
        lstatus = true;

/*        try {
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

      // test VARIABLE(S) keyword

      case 1 : {
        if (*token != STAT_word[nb_output_process == 1 ? STATW_VARIABLE : STATW_VARIABLES]) {
          status = false;
          error.correction_update(STAT_parsing[STATP_KEYWORD] ,
                                  STAT_word[nb_output_process == 1 ? STATW_VARIABLE : STATW_VARIABLES] , line , i + 1);
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
    process = new CategoricalProcess*[nb_output_process];
    for (i = 0;i < nb_output_process;i++) {
      process[i] = NULL;
    }

    index = 0;

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

        // test VARIABLE keyword

        case 0 : {
          if (*token == STAT_word[STATW_VARIABLE]) {
            index++;
          }
          else {
            status = false;
            error.correction_update(STAT_parsing[STATP_KEYWORD] , STAT_word[STATW_VARIABLE] , line , i + 1);
          }
          break;
        }

        // test observation process index

        case 1 : {
          lstatus = true;

/*          try {
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
        }

        i++;
      }

      if (i > 0) {
        if (i != 2) {
          status = false;
          error.update(STAT_parsing[STATP_FORMAT] , line);
        }

        process[index - 1] = CategoricalProcess::parsing(error , in_file , line , nb_state ,
                                                         HIDDEN_MARKOV , true);
        if (!process[index - 1]) {
          status = false;
        }
      }
    }

    if (index != nb_output_process) {
      status = false;
      error.update(STAT_parsing[STATP_FORMAT] , line);
    }

    if (!status) {
      for (i = 0;i < nb_output_process;i++) {
        delete process[i];
      }
      delete [] process;

      process = NULL;
    }
  }

  return process;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a CategoricalProcess object.
 *
 *  \param[in,out] os                    stream,
 *  \param[in]     empirical_observation pointer on the observation frequency distributions,
 *  \param[in]     marginal_distribution pointer on the marginal frequency distribution,
 *  \param[in]     exhaustive            flag detail level,
 *  \param[in]     file_flag             flag file,
 *  \param[in]     model                 model type.
 */
/*--------------------------------------------------------------*/

ostream& CategoricalProcess::ascii_print(ostream &os , FrequencyDistribution **empirical_observation ,
                                         FrequencyDistribution *marginal_distribution ,
                                         bool exhaustive , bool file_flag , model_type model) const

{
  if (observation) {
    int i , j;
    int buff , width[2];
    double *pmass , scale[NB_STATE];
    const Distribution *pobservation[NB_STATE];
    ios_base::fmtflags format_flags;


    format_flags = os.setf(ios::left , ios::adjustfield);

    for (i = 0;i < nb_state;i++) {

      // observation distributions

      os << "\n";
      switch (model) {
      case MIXTURE :
        os << STAT_word[STATW_COMPONENT];
        break;
      case HIDDEN_MARKOV :
        os << STAT_word[STATW_STATE];
        break;
      }
      os << " " << i << " " << STAT_word[STATW_OBSERVATION_DISTRIBUTION] << endl;

      pmass = observation[i]->mass + observation[i]->offset;
      for (j = observation[i]->offset;j < observation[i]->nb_value;j++) {
        if (*pmass > 0.) {
          os << STAT_word[STATW_OUTPUT] << " " << j << " : " << *pmass << endl;
        }
        pmass++;
      }

      if ((empirical_observation) && (empirical_observation[i]->nb_element > 0) && (exhaustive)) {
        os << "\n";
        if (file_flag) {
          os << "# ";
        }
        os << "   | ";
        switch (model) {
        case MIXTURE :
          os << STAT_label[STATL_COMPONENT];
          break;
        case HIDDEN_MARKOV :
          os << STAT_label[STATL_STATE];
          break;
        }
        os << " " << i << " " << STAT_label[STATL_OBSERVATION]
           << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " | ";
        switch (model) {
        case MIXTURE :
          os << STAT_label[STATL_COMPONENT];
          break;
        case HIDDEN_MARKOV :
          os << STAT_label[STATL_STATE];
          break;
        }
        os << " " << i << " " << STAT_label[STATL_OBSERVATION]
           << " " << STAT_label[STATL_DISTRIBUTION] << endl;

        observation[i]->ascii_print(os , file_flag , false , false , empirical_observation[i]);
      }
    }

    // computation of the column width

    width[0]= column_width(nb_state - 1) + ASCII_SPACE;

    width[1] = 0;
    for (i = 0;i < nb_state;i++) {
      buff= column_width(observation[i]->nb_value , observation[i]->mass);
      if (buff > width[1]) {
        width[1] = buff;
      }
    }
    width[1] += ASCII_SPACE;

    // observation probability matrix

    os << "\n";
    if (file_flag) {
      os << "# ";
    }
    os << STAT_label[STATL_OBSERVATION_PROBABILITIY_MATRIX] << endl;

    os << "\n";
    if (file_flag) {
      os << "# ";
    }
    os << setw(width[0]) << " ";
    for (i = 0;i < nb_value;i++) {
      os << setw(width[1]) << i;
    }

    for (i = 0;i < nb_state;i++) {
      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << setw(width[0]) << i;
      for (j = 0;j < nb_value;j++) {
        os << setw(width[1]) << observation[i]->mass[j];
      }
    }
    os << endl;

    if (marginal_distribution) {
      double likelihood , information;
      Test test(CHI2);


      if ((weight) && (mixture)) {
        os << "\n";
        if (file_flag) {
          os << "# ";
        }
        os << STAT_label[STATL_THEORETICAL] << " " << STAT_label[STATL_WEIGHTS] << ":";

        for (i = 0;i < nb_state;i++) {
          os << " " << weight->mass[i];
        }
        os << endl;

        likelihood = mixture->likelihood_computation(*marginal_distribution);
        information = marginal_distribution->information_computation();

        os << "\n";
        if (file_flag) {
          os << "# ";
        }
        os << STAT_label[STATL_LIKELIHOOD] << ": " << likelihood << "   ("
           << STAT_label[STATL_NORMALIZED] << ": " << likelihood / marginal_distribution->nb_element << ")" << endl;

        if (file_flag) {
          os << "# ";
        }
        os << STAT_label[STATL_MAX_LIKELIHOOD] << ": " << information << "   ("
           << STAT_label[STATL_INFORMATION] << ": " << information / marginal_distribution->nb_element << ")" << endl;

        if (file_flag) {
          os << "# ";
        }
        os << STAT_label[STATL_DEVIANCE] << ": " << 2 * (information - likelihood) << endl;

        mixture->chi2_fit(*marginal_distribution , test);
        os << "\n";
        test.ascii_print(os , file_flag);

        if (exhaustive) {
          for (i = 0;i < nb_state;i++) {
            pobservation[i] = observation[i];
            scale[i] = weight->mass[i] * marginal_distribution->nb_element;
          }

          os << "\n";
          if (file_flag) {
            os << "# ";
          }
          os << "   | " << STAT_label[STATL_MARGINAL] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
          for (i = 0;i < nb_state;i++) {
            os << " | ";
            switch (model) {
            case MIXTURE :
              os << STAT_label[STATL_COMPONENT];
              break;
            case HIDDEN_MARKOV :
              os << STAT_label[STATL_STATE];
              break;
            }
            os << " " << i << " " << STAT_label[STATL_OBSERVATION] << " " << STAT_label[STATL_DISTRIBUTION];
          }
          os << " | " << STAT_label[STATL_MIXTURE] << " | " << STAT_label[STATL_CUMULATIVE]
             << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " " << STAT_label[STATL_FUNCTION]
             << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_MIXTURE]
             << " " << STAT_label[STATL_FUNCTION] << endl;

          mixture->ascii_print(os , nb_state , pobservation , scale ,
                               file_flag , true , marginal_distribution);
        }
      }

      if ((restoration_weight) && (restoration_mixture)) {
        os << "\n";
        if (file_flag) {
          os << "# ";
        }
        os << STAT_label[STATL_RESTORATION] << " " << STAT_label[STATL_WEIGHTS] << ":";

        for (i = 0;i < nb_state;i++) {
          os << " " << restoration_weight->mass[i];
        }
        os << endl;

        likelihood = restoration_mixture->likelihood_computation(*marginal_distribution);
        information = marginal_distribution->information_computation();

        os << "\n";
        if (file_flag) {
          os << "# ";
        }
        os << STAT_label[STATL_LIKELIHOOD] << ": " << likelihood << "   ("
           << STAT_label[STATL_NORMALIZED] << ": " << likelihood / marginal_distribution->nb_element << ")" << endl;

        if (file_flag) {
          os << "# ";
        }
        os << STAT_label[STATL_MAX_LIKELIHOOD] << ": " << information << "   ("
           << STAT_label[STATL_INFORMATION] << ": " << information / marginal_distribution->nb_element << ")" << endl;

        if (file_flag) {
          os << "# ";
        }
        os << STAT_label[STATL_DEVIANCE] << ": " << 2 * (information - likelihood) << endl;

        restoration_mixture->chi2_fit(*marginal_distribution , test);
        os << "\n";
        test.ascii_print(os , file_flag);

        if (exhaustive) {
          for (i = 0;i < nb_state;i++) {
            pobservation[i] = observation[i];
            scale[i] = restoration_weight->mass[i] * marginal_distribution->nb_element;
          }

          os << "\n";
          if (file_flag) {
            os << "# ";
          }
          os << "   | " << STAT_label[STATL_MARGINAL] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
          for (i = 0;i < nb_state;i++) {
            os << " | ";
            switch (model) {
            case MIXTURE :
              os << STAT_label[STATL_COMPONENT];
              break;
            case HIDDEN_MARKOV :
              os << STAT_label[STATL_STATE];
              break;
            }
            os << " " << i << " " << STAT_label[STATL_OBSERVATION] << " " << STAT_label[STATL_DISTRIBUTION];
          }
          os << " | " << STAT_label[STATL_MIXTURE] << " | " << STAT_label[STATL_CUMULATIVE]
             << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " " << STAT_label[STATL_FUNCTION]
             << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_MIXTURE]
             << " " << STAT_label[STATL_FUNCTION] << endl;

          restoration_mixture->ascii_print(os , nb_state , pobservation , scale ,
                                           file_flag , true , marginal_distribution);
        }
      }
    }

    os.setf(format_flags , ios::adjustfield);
  }

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a CategoricalProcess object at the spreadsheet format.
 *
 *  \param[in,out] os                    stream,
 *  \param[in]     empirical_observation pointer on the observation frequency distributions,
 *  \param[in]     marginal_distribution pointer on the marginal frequency distribution,
 *  \param[in]     model                 model type.
 */
/*--------------------------------------------------------------*/

ostream& CategoricalProcess::spreadsheet_print(ostream &os ,FrequencyDistribution **empirical_observation ,
                                               FrequencyDistribution *marginal_distribution ,
                                               model_type model) const

{
  if (observation) {
    int i , j;
    double *pmass , scale[NB_STATE];
    const Distribution *pobservation[NB_STATE];


    for (i = 0;i < nb_state;i++) {
      os << "\n";
      switch (model) {
      case MIXTURE :
        os << STAT_word[STATW_COMPONENT];
        break;
      case HIDDEN_MARKOV :
        os << STAT_word[STATW_STATE];
        break;
      }
      os << " " << i << "\t" << STAT_word[STATW_OBSERVATION_DISTRIBUTION] << endl;

      pmass= observation[i]->mass+observation[i]->offset;
      for (j = observation[i]->offset;j < observation[i]->nb_value;j++) {
        if (*pmass > 0.) {
          os << STAT_word[STATW_OUTPUT] << "\t" << j << "\t" << *pmass << endl;
        }
        pmass++;
      }

      if ((empirical_observation) && (empirical_observation[i]->nb_element > 0)) {
        os << "\n\t";
        switch (model) {
        case MIXTURE :
          os << STAT_label[STATL_COMPONENT];
          break;
        case HIDDEN_MARKOV :
          os << STAT_label[STATL_STATE];
          break;
        }
        os << " " << i << " " << STAT_label[STATL_OBSERVATION]
           << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\t";
        switch (model) {
        case MIXTURE :
          os << STAT_label[STATL_COMPONENT];
          break;
        case HIDDEN_MARKOV :
          os << STAT_label[STATL_STATE];
          break;
        }
        os << " " << i << " " << STAT_label[STATL_OBSERVATION]
           << " " << STAT_label[STATL_DISTRIBUTION] << endl;

        observation[i]->spreadsheet_print(os , false , false , false , empirical_observation[i]);
      }
    }

    if (marginal_distribution) {
      double likelihood , information;
      Test test(CHI2);


      if ((weight) && (mixture)) {
        os << "\n" << STAT_label[STATL_THEORETICAL] << " " << STAT_label[STATL_WEIGHTS];
        for (i = 0;i < nb_state;i++) {
          os << "\t" << weight->mass[i];
        }
        os << endl;

        likelihood = mixture->likelihood_computation(*marginal_distribution);
        information = marginal_distribution->information_computation();

        os << "\n" << STAT_label[STATL_LIKELIHOOD] << "\t" << likelihood << "\t"
           << STAT_label[STATL_NORMALIZED] << "\t" << likelihood / marginal_distribution->nb_element << endl;
        os << STAT_label[STATL_MAX_LIKELIHOOD] << "\t" << information << "\t"
           << STAT_label[STATL_INFORMATION] << "\t" << information / marginal_distribution->nb_element << endl;
        os << STAT_label[STATL_DEVIANCE] << "\t" << 2 * (information - likelihood) << endl;

        mixture->chi2_fit(*marginal_distribution , test);
        os << "\n";
        test.spreadsheet_print(os);

        for (i = 0;i < nb_state;i++) {
          pobservation[i] = observation[i];
          scale[i] = weight->mass[i] * marginal_distribution->nb_element;
        }

        os << "\n\t" << STAT_label[STATL_MARGINAL] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
        for (i = 0;i < nb_state;i++) {
          os << "\t";
          switch (model) {
          case MIXTURE :
            os << STAT_label[STATL_COMPONENT];
            break;
          case HIDDEN_MARKOV :
            os << STAT_label[STATL_STATE];
            break;
          }
          os << " " << i << " " << STAT_label[STATL_OBSERVATION] << " " << STAT_label[STATL_DISTRIBUTION];
        }
        os << "\t" << STAT_label[STATL_MIXTURE] << "\t" << STAT_label[STATL_CUMULATIVE]
           << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " " << STAT_label[STATL_FUNCTION]
           << "\t" << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_MIXTURE]
           << " " << STAT_label[STATL_FUNCTION] << endl;

        mixture->spreadsheet_print(os , nb_state , pobservation , scale , true ,
                                   marginal_distribution);
      }

      if ((restoration_weight) && (restoration_mixture)) {
        os << "\n" << STAT_label[STATL_RESTORATION] << " " << STAT_label[STATL_WEIGHTS];
        for (i = 0;i < nb_state;i++) {
          os << "\t" << restoration_weight->mass[i];
        }
        os << endl;

        likelihood = restoration_mixture->likelihood_computation(*marginal_distribution);
        information = marginal_distribution->information_computation();

        os << "\n" << STAT_label[STATL_LIKELIHOOD] << "\t" << likelihood << "\t"
           << STAT_label[STATL_NORMALIZED] << "\t" << likelihood / marginal_distribution->nb_element << endl;
        os << STAT_label[STATL_MAX_LIKELIHOOD] << "\t" << information << "\t"
           << STAT_label[STATL_INFORMATION] << "\t" << information / marginal_distribution->nb_element << endl;
        os << STAT_label[STATL_DEVIANCE] << "\t" << 2 * (information - likelihood) << endl;

        restoration_mixture->chi2_fit(*marginal_distribution , test);
        os << "\n";
        test.spreadsheet_print(os);

        for (i = 0;i < nb_state;i++) {
          pobservation[i] = observation[i];
          scale[i] = restoration_weight->mass[i] * marginal_distribution->nb_element;
        }

        os << "\n\t" << STAT_label[STATL_MARGINAL] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
        for (i = 0;i < nb_state;i++) {
          os << "\t";
          switch (model) {
          case MIXTURE :
            os << STAT_label[STATL_COMPONENT];
            break;
          case HIDDEN_MARKOV :
            os << STAT_label[STATL_STATE];
            break;
          }
          os << " " << i << " " << STAT_label[STATL_OBSERVATION] << " " << STAT_label[STATL_DISTRIBUTION];
        }
        os << "\t" << STAT_label[STATL_MIXTURE] << "\t" << STAT_label[STATL_CUMULATIVE]
           << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " " << STAT_label[STATL_FUNCTION]
           << "\t" << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_MIXTURE]
           << " " << STAT_label[STATL_FUNCTION] << endl;

        restoration_mixture->spreadsheet_print(os , nb_state , pobservation , scale , true ,
                                               marginal_distribution);
      }
    }
  }

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Plot of a CategoricalProcess object using Gnuplot.
 *
 *  \param[in] prefix                file prefix,
 *  \param[in] title                 figure title,
 *  \param[in] process               observation process index,
 *  \param[in] empirical_observation pointer on the observation frequency distributions,
 *  \param[in] marginal_distribution pointer on the marginal frequency distribution,
 *  \param[in] model                 model type.
 *
 *  \return                          error status.
 */
/*--------------------------------------------------------------*/

bool CategoricalProcess::plot_print(const char *prefix , const char *title , int process ,
                                    FrequencyDistribution **empirical_observation ,
                                    FrequencyDistribution *marginal_distribution ,
                                    model_type model) const

{
  bool status = false;

  if (observation) {
    int i , j , k , m;
    int nb_dist , nb_histo , *dist_nb_value;
    double *scale;
    const Distribution **pdist;
    const FrequencyDistribution **phisto;
    ostringstream data_file_name;


    // writing of the data files

    pdist = new const Distribution*[3 * nb_state + 2];
    dist_nb_value = new int[3 * nb_state + 2];
    scale = new double[3 * nb_state + 2];
    phisto = new const FrequencyDistribution*[nb_state + 2];

    nb_histo = 0;
    for (i = 0;i < nb_state;i++) {
      pdist[i] = observation[i];
      dist_nb_value[i] = observation[i]->nb_value;

      if ((empirical_observation) && (empirical_observation[i]->nb_element > 0)) {
        phisto[nb_histo++] = empirical_observation[i];
        scale[i] = empirical_observation[i]->nb_element;
      }
      else {
        scale[i] = 1.;
      }
    }

    if (marginal_distribution) {
      nb_dist = nb_state;

      if ((weight) && (mixture)) {
        for (i = 0;i < nb_state;i++) {
          pdist[nb_dist] = observation[i];
          dist_nb_value[nb_dist] = observation[i]->nb_value;
          scale[nb_dist++] = weight->mass[i] * marginal_distribution->nb_element;
        }

        pdist[nb_dist] = mixture;
        dist_nb_value[nb_dist] = mixture->nb_value;
        phisto[nb_histo++] = marginal_distribution;
        scale[nb_dist++] = marginal_distribution->nb_element;
      }

      if ((restoration_weight) && (restoration_mixture)) {
        for (i = 0;i < nb_state;i++) {
          pdist[nb_dist] = observation[i];
          dist_nb_value[nb_dist] = observation[i]->nb_value;
          scale[nb_dist++] = restoration_weight->mass[i] * marginal_distribution->nb_element;
        }

        pdist[nb_dist] = restoration_mixture;
        dist_nb_value[nb_dist] = restoration_mixture->nb_value;
        phisto[nb_histo++] = marginal_distribution;
        scale[nb_dist++] = marginal_distribution->nb_element;
      }
    }

    data_file_name << prefix << process << ".dat";
    status = stat_tool::plot_print((data_file_name.str()).c_str(), nb_dist , pdist ,
                                   scale , dist_nb_value , nb_histo , phisto);

    if (status) {

      // writing of the script files

      for (i = 0;i < 2;i++) {
        ostringstream file_name[2];

        switch (i) {
        case 0 :
          file_name[0] << prefix << process << ".plot";
          break;
        case 1 :
          file_name[0] << prefix << process << ".print";
          break;
        }

        ofstream out_file((file_name[0].str()).c_str());

        if (i == 1) {
          out_file << "set terminal postscript" << endl;
          file_name[1] << label(prefix) << process << 0 << ".ps";
          out_file << "set output \"" << file_name[1].str() << "\"\n\n";
        }

        out_file << "set border 15 lw 0\n" << "set tics out\n" << "set xtics nomirror\n"
                 << "set title";

        if (title) {
          out_file << " \"" << title << " - " << STAT_label[STATL_OUTPUT_PROCESS]
                   << " " << process << "\"";
        }
        out_file << "\n\n";

        j = 0;

        for (k = 0;k < nb_state;k++) {
          if (dist_nb_value[k] - 1 < TIC_THRESHOLD) {
            out_file << "set xtics 0,1" << endl;
          }

          if ((empirical_observation) && (empirical_observation[k]->nb_element > 0)) {
            out_file << "plot [0:" << dist_nb_value[k] - 1 << "] [0:"
                     << (int)(MAX(phisto[j]->max , pdist[k]->max * scale[k]) * YSCALE) + 1
                     << "] \"" << label((data_file_name.str()).c_str()) << "\" using " << j + 1
                     << " title \"";
            switch (model) {
            case MIXTURE :
              out_file << STAT_label[STATL_COMPONENT];
              break;
            case HIDDEN_MARKOV :
              out_file << STAT_label[STATL_STATE];
              break;
            }
            out_file << " " << k << " " << STAT_label[STATL_OBSERVATION]
                     << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\" with impulses,\\" << endl;
            out_file << "\"" << label((data_file_name.str()).c_str()) << "\" using " << nb_histo + k + 1
                     << " title \"";
            switch (model) {
            case MIXTURE :
              out_file << STAT_label[STATL_COMPONENT];
              break;
            case HIDDEN_MARKOV :
              out_file << STAT_label[STATL_STATE];
              break;
            }
            out_file << " " << k << " " << STAT_label[STATL_OBSERVATION]
                     << " " << STAT_label[STATL_DISTRIBUTION] << "\" with linespoints" << endl;
            j++;
          }

          else {
            out_file << "plot [0:" << dist_nb_value[k] - 1 << "] [0:"
                     << MIN(pdist[k]->max * YSCALE , 1.) << "] \""
                     << label((data_file_name.str()).c_str()) << "\" using " << nb_histo + k + 1
                     << " title \"";
            switch (model) {
            case MIXTURE :
              out_file << STAT_label[STATL_COMPONENT];
              break;
            case HIDDEN_MARKOV :
              out_file << STAT_label[STATL_STATE];
              break;
            }
            out_file << " " << k << " " << STAT_label[STATL_OBSERVATION]
                     << " " << STAT_label[STATL_DISTRIBUTION] << "\" with linespoints" << endl;
          }

          if (dist_nb_value[k] - 1 < TIC_THRESHOLD) {
            out_file << "set xtics autofreq" << endl;
          }

          if ((i == 0) && (k < nb_state - 1)) {
             out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
          }
        }

        if (marginal_distribution) {
          k = nb_state;

          if ((weight) && (mixture)) {
            if (i == 0) {
              out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
            }
            out_file << endl;

            out_file << "set title \"";
            if (title) {
              out_file << title << " - ";
            }
            out_file << STAT_label[STATL_OUTPUT_PROCESS] << " " << process << " - "
                     << STAT_label[STATL_THEORETICAL] << " " << STAT_label[STATL_WEIGHTS] << "\"\n\n";

            if (nb_value - 1 < TIC_THRESHOLD) {
              out_file << "set xtics 0,1" << endl;
            }

            out_file << "plot [0:" << nb_value - 1 << "] [0:"
                     << (int)(MAX(marginal_distribution->max , mixture->max * marginal_distribution->nb_element) * YSCALE) + 1
                     << "] \"" << label((data_file_name.str()).c_str()) << "\" using " << j + 1
                     << " title \"" << STAT_label[STATL_MARGINAL] << " "
                     << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\" with impulses,\\" << endl;
            j++;

            for (m = 0;m < nb_state;m++) {
              out_file << "\"" << label((data_file_name.str()).c_str()) << "\" using " << nb_histo + k + 1
                       << " title \"";
              switch (model) {
              case MIXTURE :
                out_file << STAT_label[STATL_COMPONENT];
                break;
              case HIDDEN_MARKOV :
                out_file << STAT_label[STATL_STATE];
                break;
              }
              out_file << " " << m << " " << STAT_label[STATL_OBSERVATION] << " " << STAT_label[STATL_DISTRIBUTION]
                       << "\" with linespoints,\\" << endl;
              k++;
            }

            out_file << "\"" << label((data_file_name.str()).c_str()) << "\" using " << nb_histo + k + 1
                     << " title \"" << STAT_label[STATL_MIXTURE] << "\" with linespoints" << endl;
            k++;

            if (nb_value - 1 < TIC_THRESHOLD) {
              out_file << "set xtics autofreq" << endl;
            }
          }

          if ((restoration_weight) && (restoration_mixture)) {
            if (i == 0) {
              out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
            }
            out_file << endl;

            out_file << "set title \"";
            if (title) {
              out_file << title << " - ";
            }
            out_file << STAT_label[STATL_OUTPUT_PROCESS] << " " << process << " - "
                     << STAT_label[STATL_RESTORATION] << " " << STAT_label[STATL_WEIGHTS] << "\"\n\n";

            if (nb_value - 1 < TIC_THRESHOLD) {
              out_file << "set xtics 0,1" << endl;
            }

            out_file << "plot [0:" << nb_value - 1 << "] [0:"
                     << (int)(MAX(marginal_distribution->max , restoration_mixture->max * marginal_distribution->nb_element) * YSCALE) + 1
                     << "] \"" << label((data_file_name.str()).c_str()) << "\" using " << j + 1
                     << " title \"" << STAT_label[STATL_MARGINAL] << " "
                     << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\" with impulses,\\" << endl;
            j++;

            for (m = 0;m < nb_state;m++) {
              out_file << "\"" << label((data_file_name.str()).c_str()) << "\" using " << nb_histo + k + 1
                       << " title \"";
              switch (model) {
              case MIXTURE :
                out_file << STAT_label[STATL_COMPONENT];
                break;
              case HIDDEN_MARKOV :
                out_file << STAT_label[STATL_STATE];
                break;
              }
              out_file << " " << m << " " << STAT_label[STATL_OBSERVATION] << " " << STAT_label[STATL_DISTRIBUTION]
                       << "\" with linespoints,\\" << endl;
              k++;
            }

            out_file << "\"" << label((data_file_name.str()).c_str()) << "\" using " << nb_histo + k + 1
                     << " title \"" << STAT_label[STATL_MIXTURE] << "\" with linespoints" << endl;
            k++;

            if (nb_value - 1 < TIC_THRESHOLD) {
              out_file << "set xtics autofreq" << endl;
            }
          }
        }

        if (i == 1) {
          out_file << "\nset terminal x11" << endl;
        }

        out_file << "\npause 0 \"" << STAT_label[STATL_END] << "\"" << endl;
      }
    }

    delete [] pdist;
    delete [] dist_nb_value;
    delete [] scale;
    delete [] phisto;
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Plot of a CategoricalProcess object.
 *
 *  \param[in] plot                  reference on a MultiPlotSet object,
 *  \param[in] index                 MultiPlot index,
 *  \param[in] process               observation process index,
 *  \param[in] empirical_observation pointer on the observation frequency distributions,
 *  \param[in] marginal_distribution pointer on the marginal frequency distribution,
 *  \param[in] model                 model type.
 */
/*--------------------------------------------------------------*/

void CategoricalProcess::plotable_write(MultiPlotSet &plot , int &index , int process ,
                                        FrequencyDistribution **empirical_observation ,
                                        FrequencyDistribution *marginal_distribution ,
                                        model_type model) const

{
  if (observation) {
    int i , j;
    double scale , max;
    ostringstream title , legend;


    plot.variable_nb_viewpoint[process] = 1;

    // computation of the number of plots

/*    if (empirical_observation) {
      nb_plot_set = nb_state;
    }
    else {
      nb_plot_set = 1;
    }

    if (marginal_distribution) {
      if ((weight) && (mixture)) {
        nb_plot_set++;
      }
      if ((restoration_weight) && (restoration_mixture)) {
        nb_plot_set++;
      }
    } */

    if (empirical_observation) {
      for (i = 0;i < nb_state;i++) {

        // fit of observation distribution

        plot.variable[index] = process;
        plot.viewpoint[index] = OBSERVATION;

        title.str("");
        title << STAT_label[STATL_OUTPUT_PROCESS] << " " << process;
        plot[index].title = title.str();

        plot[index].xrange = Range(0 , observation[i]->nb_value - 1);
        if (observation[i]->nb_value - 1 < TIC_THRESHOLD) {
          plot[index].xtics = 1;
        }

        if (empirical_observation[i]->nb_element > 0) {
          scale = empirical_observation[i]->nb_element;
          plot[index].yrange = Range(0 , ceil(MAX(empirical_observation[i]->max ,
                                                  observation[i]->max * scale) * YSCALE));

          plot[index].resize(2);

          legend.str("");
          switch (model) {
          case MIXTURE :
            legend << STAT_label[STATL_COMPONENT];
            break;
          case HIDDEN_MARKOV :
            legend << STAT_label[STATL_STATE];
            break;
          }
          legend << " " << i << " " << STAT_label[STATL_OBSERVATION]
                 << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
          plot[index][0].legend = legend.str();

          plot[index][0].style = "impulses";

          empirical_observation[i]->plotable_frequency_write(plot[index][0]);
          j = 1;
        }

        else {
          scale = 1;
          plot[index].yrange = Range(0 , MIN(observation[i]->max * YSCALE , 1.));

          plot[index].resize(1);
          j = 0;
        }

        legend.str("");
        switch (model) {
        case MIXTURE :
          legend << STAT_label[STATL_COMPONENT];
          break;
        case HIDDEN_MARKOV :
          legend << STAT_label[STATL_STATE];
          break;
        }
        legend << " " << i << " " << STAT_label[STATL_OBSERVATION] << " " << STAT_label[STATL_DISTRIBUTION];
        plot[index][j].legend = legend.str();

        plot[index][j].style = "linespoints";

        observation[i]->plotable_mass_write(plot[index][j] , scale);
        index++;
      }
    }

    else {

      // observation distributions

      plot.variable[index] = process;
      plot.viewpoint[index] = OBSERVATION;

      title.str("");
      title << STAT_label[STATL_OUTPUT_PROCESS] << " " << process;
      plot[index].title = title.str();

      plot[index].xrange = Range(0 , nb_value - 1);
      if (nb_value - 1 < TIC_THRESHOLD) {
        plot[index].xtics = 1;
      }

      max = observation[0]->max;
      for (i = 1;i < nb_state;i++) {
        if (observation[i]->max > max) {
          max = observation[i]->max;
        }
      }
      plot[index].yrange = Range(0 , MIN(max * YSCALE , 1.));

      plot[index].resize(nb_state);

      for (i = 0;i < nb_state;i++) {
        legend.str("");
        switch (model) {
        case MIXTURE :
          legend << STAT_label[STATL_COMPONENT];
          break;
        case HIDDEN_MARKOV :
          legend << STAT_label[STATL_STATE];
          break;
        }
        legend << " " << i << " " << STAT_label[STATL_OBSERVATION] << " " << STAT_label[STATL_DISTRIBUTION];
        plot[index][i].legend = legend.str();

        plot[index][i].style = "linespoints";

        observation[i]->plotable_mass_write(plot[index][i]);
      }

      index++;
    }

    if (marginal_distribution) {
      if ((weight) && (mixture)) {

        // fit of the mixture of observation distributions (theoretical weights)

        title.str("");
        title << STAT_label[STATL_OUTPUT_PROCESS] << " " << process << " - "
              << STAT_label[STATL_THEORETICAL] << " " << STAT_label[STATL_WEIGHTS];
        plot[index].title = title.str();

        plot[index].xrange = Range(0 , nb_value - 1);
        if (nb_value - 1 < TIC_THRESHOLD) {
          plot[index].xtics = 1;
        }

        plot[index].yrange = Range(0 , ceil(MAX(marginal_distribution->max ,
                                                mixture->max * marginal_distribution->nb_element) * YSCALE));

        plot[index].resize(nb_state + 2);

        legend.str("");
        legend << STAT_label[STATL_MARGINAL] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
        plot[index][0].legend = legend.str();

        plot[index][0].style = "impulses";

        marginal_distribution->plotable_frequency_write(plot[index][0]);

        for (i = 0;i < nb_state;i++) {
          legend.str("");
          switch (model) {
          case MIXTURE :
            legend << STAT_label[STATL_COMPONENT];
            break;
          case HIDDEN_MARKOV :
            legend << STAT_label[STATL_STATE];
            break;
          }
          legend << " " << i << " " << STAT_label[STATL_OBSERVATION] << " " << STAT_label[STATL_DISTRIBUTION];
          observation[i]->plot_title_print(legend);
          plot[index][i + 1].legend = legend.str();

          plot[index][i + 1].style = "linespoints";

          observation[i]->plotable_mass_write(plot[index][i + 1] ,
                                              weight->mass[i] * marginal_distribution->nb_element);
        }

        plot[index][nb_state + 2].legend = STAT_label[STATL_MIXTURE];

        plot[index][nb_state + 2].style = "linespoints";

        mixture->plotable_mass_write(plot[index][nb_state + 2] ,
                                     marginal_distribution->nb_element);

        index++;
      }

      if ((restoration_weight) && (restoration_mixture)) {

        // fit of the mixture of observation distributions (weights deduced from the restoration)

        title.str("");
        if (process > 0) {
          title << STAT_label[STATL_OUTPUT_PROCESS] << " " << process;
        }
        title << STAT_label[STATL_RESTORATION] << " " << STAT_label[STATL_WEIGHTS];
        plot[index].title = title.str();

        plot[index].xrange = Range(0 , nb_value - 1);
        if (nb_value - 1 < TIC_THRESHOLD) {
          plot[index].xtics = 1;
        }

        plot[index].yrange = Range(0 , ceil(MAX(marginal_distribution->max ,
                                                restoration_mixture->max * marginal_distribution->nb_element) * YSCALE));

        plot[index].resize(nb_state + 2);

        legend.str("");
        legend << STAT_label[STATL_MARGINAL] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
        plot[index][0].legend = legend.str();

        plot[index][0].style = "impulses";

        marginal_distribution->plotable_frequency_write(plot[index][0]);

        for (i = 0;i < nb_state;i++) {
          legend.str("");
          switch (model) {
          case MIXTURE :
            legend << STAT_label[STATL_COMPONENT];
            break;
          case HIDDEN_MARKOV :
            legend << STAT_label[STATL_STATE];
            break;
          }
          legend << " " << i << " " << STAT_label[STATL_OBSERVATION] << " " << STAT_label[STATL_DISTRIBUTION];
          observation[i]->plot_title_print(legend);
          plot[index][i + 1].legend = legend.str();

          plot[index][i + 1].style = "linespoints";

          observation[i]->plotable_mass_write(plot[index][i + 1] ,
                                              restoration_weight->mass[i] * marginal_distribution->nb_element);
        }

        plot[index][nb_state + 2].legend = STAT_label[STATL_MIXTURE];

        plot[index][nb_state + 2].style = "linespoints";

        restoration_mixture->plotable_mass_write(plot[index][nb_state + 2] ,
                                                 marginal_distribution->nb_element);

        index++;
      }
    }
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief test of the overlap between the categories observed in each state.
 *
 *  \return hidden process or not.
 */
/*--------------------------------------------------------------*/

bool CategoricalProcess::test_hidden() const

{
  bool hidden = false;
  int i , j;
  int nb_occurrence;


  for (i = 0;i < nb_value;i++) {
    nb_occurrence = 0;
    for (j = 0;j < nb_state;j++) {
      if (observation[j]->mass[i] > 0.) {
        nb_occurrence++;
      }
    }

    if (nb_occurrence > 1) {
      hidden = true;
      break;
    }
  }

  return hidden;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Application of a threshold on the observation probabilities.
 *
 *  \param[in] min_probability minimum probability.
 */
/*--------------------------------------------------------------*/

void CategoricalProcess::thresholding(double min_probability)

{
  bool stop;
  int i , j;
  int nb_correction;
  double norm , *pmass;


  if (min_probability > THRESHOLDING_FACTOR / (double)nb_value) {
    min_probability = THRESHOLDING_FACTOR / (double)nb_value;
  }

  for (i = 0;i < nb_state;i++) {
    do {
      stop = true;
      pmass = observation[i]->mass;
      nb_correction = 0;
      norm = 0.;

      for (j = 0;j < nb_value;j++) {
        if (*pmass <= min_probability) {
          nb_correction++;
          *pmass++ = min_probability;
        }
        else {
          norm += *pmass++;
        }
      }

      if (nb_correction > 0) {
        pmass = observation[i]->mass;
        for (j = 0;j < nb_value;j++) {
          if (*pmass > min_probability) {
            *pmass *= (1. - nb_correction * min_probability) / norm;
            if (*pmass < min_probability) {
              stop = false;
            }
          }
          pmass++;
        }
      }
    }
    while (!stop);

    if (nb_correction > 0) {
      observation[i]->cumul_computation();
      observation[i]->max_computation();
    }
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Permutation of observation distributions. The permutation validity
 *         should be checked by the calling function.
 *
 *  \param[in] permut permutation.
 */
/*--------------------------------------------------------------*/

void CategoricalProcess::state_permutation(int *permut) const

{
  int i;
  Distribution **pobservation = new Distribution*[nb_state];


  for (i = 0;i < nb_state;i++) {
    pobservation[permut[i]] = observation[i];
  }
  for (i = 0;i < nb_state;i++) {
    observation[i] = pobservation[i];
  }
  delete [] pobservation;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the number of free parameters of a categorical
 *         observation process.
 *
 *  \param[in] min_probability minimum probability.
 *
 *  \return                    number of free parameters.
 */
/*--------------------------------------------------------------*/

int CategoricalProcess::nb_parameter_computation(double min_probability) const

{
  int i , j;
  int nb_parameter = 0;


  for (i = 0;i < nb_state;i++) {
    for (j = 0;j < nb_value;j++) {
      if (observation[i]->mass[j] > min_probability) {
        nb_parameter++;
      }
    }

    nb_parameter--;
  }

  return nb_parameter;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of a mixture of categorical observation distributions.
 *
 *  \param[in] pweight pointer on the weight distribution.
 *
 *  \return            mixture of categorical observation distributions.
 */
/*--------------------------------------------------------------*/

Distribution* CategoricalProcess::mixture_computation(Distribution *pweight)

{
  int i , j;
  Distribution *mixture;


  mixture = new Distribution(nb_value);

  mixture->offset = 0;

  for (i = 0;i < nb_value;i++) {
    mixture->mass[i] = 0.;
    for (j = 0;j < nb_state;j++) {
      mixture->mass[i] += pweight->mass[j] * observation[j]->mass[i];
    }
  }

  mixture->cumul_computation();

  mixture->max_computation();
  mixture->mean_computation();
  mixture->variance_computation();

  return mixture;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Initialization of the observation probabilities.
 */
/*--------------------------------------------------------------*/

void CategoricalProcess::init()

{
  bool *active_output;
  int i , j , k;
  int nb_noise , *state_nb_value;
  double cumul , noise_proba = NOISE_PROBABILITY , *pmass;


  active_output = new bool[nb_value];
  state_nb_value = new int[nb_state];

  cumul = 0.;
  pmass = observation[0]->mass;
  for (i = 0;i < nb_value;i++) {
    cumul += *pmass++;
  }

  // uniform initialization of the observation probabilities

  if (cumul > 1. - DOUBLE_ERROR) {
    for (i = 0;i < nb_state;i++) {
      pmass = observation[i]->mass;
      state_nb_value[i] = 0;
      for (j = 0;j < nb_value;j++) {
        if (*pmass++ > 0.) {
          state_nb_value[i]++;
        }
      }

      pmass = observation[i]->mass;
      for (j = 0;j < nb_value;j++) {
        if (*pmass > 0.) {
          *pmass = 1. / (double)state_nb_value[i];
        }
        pmass++;
      }
    }
  }

  else {
    for (i = 0;i < nb_state;i++) {
      state_nb_value[i] = nb_value;
      pmass = observation[i]->mass;
      for (j = 0;j < nb_value;j++) {
        *pmass++ = 1. / (double)nb_value;
      }
    }
  }

  // perturbation of the observation distributions

  nb_noise = nb_value;

  for (i = 0;i < nb_value;i++) {
    active_output[i] = false;
    for (j = 0;j < nb_state;j++) {
      if ((observation[j]->mass[i] > 0.) && (observation[j]->mass[i] < 1.)) {
        active_output[i] = true;
        nb_noise--;
        break;
      }
    }
  }

  for (i = 0;i < nb_state;i++) {
    if (state_nb_value[i] > 1) {
      pmass = observation[i]->mass;
      for (j = 0;j < nb_value;j++) {
        if ((active_output[j]) && (*pmass > 0.)) {
          active_output[j] = false;
          nb_noise++;
          break;
        }
        pmass++;
      }

      *pmass += noise_proba;
      pmass = observation[i]->mass;
      for (k = 0;k < j;k++) {
        if (*pmass > 0.) {
          *pmass -= noise_proba / (state_nb_value[i] - 1);
        }
        pmass++;
      }
      pmass++;
      for (k = j + 1;k < nb_value;k++) {
        if (*pmass > 0.) {
          *pmass -= noise_proba / (state_nb_value[i] - 1);
        }
        pmass++;
      }

      if ((i < nb_state - 1) && (state_nb_value[i + 1] > 1)) {
        for (j = 0;j < nb_value;j++) {
          if (active_output[j]) {
            break;
          }
        }

        if (j == nb_value) {
          for (j = 0;j < nb_value;j++) {
            for (k = 0;k < nb_state;k++) {
              if ((observation[k]->mass[j] > 0.) && (observation[k]->mass[j] < 1.)) {
                active_output[j] = true;
                nb_noise--;
                break;
              }
            }
          }
          noise_proba += NOISE_PROBABILITY;
        }
      }
    }
  }

  delete [] active_output;
  delete [] state_nb_value;
}


};  // namespace stat_tool
