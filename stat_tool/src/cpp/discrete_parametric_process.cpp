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



#include <sstream>
#include <iomanip>

#include <boost/tokenizer.hpp>

#include "markovian.h"
#include "stat_label.h"

using namespace std;
using namespace boost;


namespace stat_tool {



/*--------------------------------------------------------------*/
/**
 *  \brief Constructor of the DiscreteParametricProcess class.
 *
 *  \param[in] inb_state number of states,
 *  \param[in] inb_value number of values.
 */
/*--------------------------------------------------------------*/

DiscreteParametricProcess::DiscreteParametricProcess(int inb_state , int inb_value)

{
  int i;


  nb_state = inb_state;
  nb_value = inb_value;

  if (nb_state > 0) {
    observation = new DiscreteParametric*[nb_state];
    if (nb_value > 0) {
      for (i = 0;i < nb_state;i++) {
        observation[i] = new DiscreteParametric(UNIFORM , 0 , nb_value , D_DEFAULT , D_DEFAULT);
      }
    }

    else {
      for (i = 0;i < nb_state;i++) {
        observation[i] = NULL;
      }
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
 *  \brief Constructor of the DiscreteParametricProcess class.
 *
 *  \param[in] inb_state    number of states,
 *  \param[in] pobservation pointer on the observation distributions.
 */
/*--------------------------------------------------------------*/

DiscreteParametricProcess::DiscreteParametricProcess(int inb_state , DiscreteParametric **pobservation)

{
  int i;


  nb_state = inb_state;

  nb_value = 0;
  for (i = 0;i < nb_state;i++) {
    if (pobservation[i]->nb_value > nb_value) {
      nb_value = pobservation[i]->nb_value;
    }
  }

  observation = new DiscreteParametric*[nb_state];
  for (i = 0;i < nb_state;i++) {
    observation[i] = new DiscreteParametric(*pobservation[i] , DISTRIBUTION_COPY , nb_value);
  }

  weight = NULL;
  mixture = NULL;
  restoration_weight = NULL;
  restoration_mixture = NULL;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Copy of a DiscreteParametricProcess object.
 *
 *  \param[in] process reference on a DiscreteParametricProcess object.
 */
/*--------------------------------------------------------------*/

void DiscreteParametricProcess::copy(const DiscreteParametricProcess &process)

{
  int i;


  nb_state = process.nb_state;
  nb_value = process.nb_value;

  observation = new DiscreteParametric*[nb_state];
  for (i = 0;i < nb_state;i++) {
    observation[i] = new DiscreteParametric(*(process.observation[i]) , DISTRIBUTION_COPY , nb_value);
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
 *  \brief Destruction of the data members of a DiscreteParametricProcess object.
 */
/*--------------------------------------------------------------*/

void DiscreteParametricProcess::remove()

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
 *  \brief Destructor of the DiscreteParametricProcess class.
 */
/*--------------------------------------------------------------*/

DiscreteParametricProcess::~DiscreteParametricProcess()

{
  remove();
}


/*--------------------------------------------------------------*/
/**
 *  \brief Assignment operator of the DiscreteParametricProcess class.
 *
 *  \param[in] process reference on a DiscreteParametricProcess object.
 *
 *  \return            DiscreteParametricProcess object.
 */
/*--------------------------------------------------------------*/

DiscreteParametricProcess& DiscreteParametricProcess::operator=(const DiscreteParametricProcess &process)

{
  if (&process != this) {
    remove();
    copy(process);
  }

  return *this;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Analysis of the format of discrete parametric observation distributions.
 *
 *  \param[in] error           reference on a StatError object,
 *  \param[in] in_file         stream,
 *  \param[in] line            reference on the file line index,
 *  \param[in] nb_state        number of states,
 *  \param[in] model           model type,
 *  \param[in] cumul_threshold threshold on the cumulative distribution functions.
 *
 *  \return                    DiscreteParametricProcess object.
 */
/*--------------------------------------------------------------*/

DiscreteParametricProcess* DiscreteParametricProcess::parsing(StatError &error , ifstream &in_file ,
                                                              int &line , int nb_state ,
                                                              model_type model , double cumul_threshold)

{
  string buffer;
  size_t position;
  typedef tokenizer<char_separator<char>> tokenizer;
  char_separator<char> separator(" \t");
  bool status = true , lstatus;
  int i , j;
  int index;
  DiscreteParametric **dist;
  DiscreteParametricProcess *process;


  process = NULL;

  dist = new DiscreteParametric*[nb_state];
  for (i = 0;i < nb_state;i++) {
    dist[i] = NULL;
  }

  for (i = 0;i < nb_state;i++) {
    while (getline(in_file , buffer)) {
      line++;

#     ifdef DEBUG
      cout << line << "  " << buffer << endl;
#     endif

      position = buffer.find('#');
      if (position != string::npos) {
        buffer.erase(position);
      }
      j = 0;

      tokenizer tok_buffer(buffer , separator);

      for (tokenizer::iterator token = tok_buffer.begin();token != tok_buffer.end();token++) {
        switch (j) {

        // test COMPONENT/STATE keyword

        case 0 : {
          switch (model) {

          case MIXTURE : {
            if (*token != STAT_word[STATW_COMPONENT]) {
              status = false;
              error.correction_update(STAT_parsing[STATP_KEYWORD] ,
                                      STAT_word[STATW_COMPONENT] , line , j + 1);
            }
            break;
          }

          case HIDDEN_MARKOV : {
            if (*token != STAT_word[STATW_STATE]) {
              status = false;
              error.correction_update(STAT_parsing[STATP_KEYWORD] ,
                                      STAT_word[STATW_STATE] , line , j + 1);
            }
            break;
          }
          }
          break;
        }

        // test component/state index

        case 1 : {
          lstatus = true;

/*          try {
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

            switch (model) {
            case MIXTURE :
              error.correction_update(STAT_parsing[STATP_COMPONENT_INDEX] , i , line , j + 1);
              break;
            case HIDDEN_MARKOV :
              error.correction_update(STAT_parsing[STATP_STATE_INDEX] , i , line , j + 1);
              break;
            }
          }
          break;
        }

        // test OBSERVATION_DISTRIBUTION keyword

        case 2 : {
          if (*token != STAT_word[STATW_OBSERVATION_DISTRIBUTION]) {
            status = false;
            error.correction_update(STAT_parsing[STATP_KEYWORD] , STAT_word[STATW_OBSERVATION_DISTRIBUTION] , line , j + 1);
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

        dist[i] = DiscreteParametric::parsing(error , in_file , line ,
                                              UNIFORM , cumul_threshold);
        if (!dist[i]) {
          status = false;
        }

        break;
      }
    }

    if ((j == 0) && (!dist[i])) {
      status = false;
      error.update(STAT_parsing[STATP_FORMAT] , line);
    }
  }

  if (status) {
    process = new DiscreteParametricProcess(nb_state , dist);
  }

  for (i = 0;i < nb_state;i++) {
    delete dist[i];
  }
  delete [] dist;

  return process;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a DiscreteParametricProcess object.
 *
 *  \param[in,out] os                    stream,
 *  \param[in]     empirical_observation pointer on the observation frequency distributions,
 *  \param[in]     marginal_distribution pointer on the marginal frequency distribution,
 *  \param[in]     exhaustive            flag detail level,
 *  \param[in]     file_flag             flag file,
 *  \param[in]     model                 model type.
 */
/*--------------------------------------------------------------*/

ostream& DiscreteParametricProcess::ascii_print(ostream &os , FrequencyDistribution **empirical_observation ,
                                                FrequencyDistribution *marginal_distribution ,
                                                bool exhaustive , bool file_flag , model_type model) const

{
  int i , j;
  double scale[NB_STATE];
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
    os << " " << i << " " << STAT_word[STATW_OBSERVATION_DISTRIBUTION] << endl;
    observation[i]->ascii_print(os);
    observation[i]->ascii_parametric_characteristic_print(os , false , file_flag);

    if (empirical_observation) {
      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      switch (model) {
      case MIXTURE :
        os << STAT_label[STATL_COMPONENT];
        break;
      case HIDDEN_MARKOV :
        os << STAT_label[STATL_STATE];
        break;
      }
      os << " " << i << " " << STAT_label[STATL_OBSERVATION]
         << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";
      empirical_observation[i]->ascii_characteristic_print(os , false , file_flag);

      if (exhaustive) {
        os << "\n";
        if (file_flag) {
          os << "# ";
        }
        os << "  ";
        if (empirical_observation[i]->nb_element > 0) {
          os << " | ";
          switch (model) {
          case MIXTURE :
            os << STAT_label[STATL_COMPONENT];
            break;
          case HIDDEN_MARKOV :
            os << STAT_label[STATL_STATE];
            break;
          }
          os << " " << i << " " << STAT_label[STATL_OBSERVATION]
             << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
        }
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
        if (empirical_observation[i]->nb_element > 0) {
          os << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " "
             << STAT_label[STATL_FUNCTION];
        }
        os << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_DISTRIBUTION] << " "
           << STAT_label[STATL_FUNCTION] << endl;

        observation[i]->Distribution::ascii_print(os , file_flag , true , false ,
                                                  (empirical_observation[i]->nb_element > 0 ? empirical_observation[i] : NULL));
      }
    }
  }

  if ((!empirical_observation) && (exhaustive)) {
    int buff , width[3];
    ios_base::fmtflags format_flags;

 
    format_flags = os.setf(ios::right , ios::adjustfield);

    // computation of the column widths

    width[0] = column_width(nb_value - 1);

    width[1] = column_width(observation[0]->nb_value , observation[0]->mass);
    for (i = 1;i < nb_state;i++) {
      buff = column_width(observation[i]->nb_value , observation[i]->mass);
      if (buff > width[1]) {
        width[1] = buff;
      }
    }
    width[1] += ASCII_SPACE;

    width[2] = column_width(observation[0]->nb_value , observation[0]->cumul);
    for (i = 1;i < nb_state;i++) {
      buff = column_width(observation[i]->nb_value , observation[i]->cumul);
      if (buff > width[2]) {
        width[2] = buff;
      }
    }
    width[2] += ASCII_SPACE;

    os << "\n";
    if (file_flag) {
      os << "# ";
    }
    os << "  ";
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
      os << " " << i << " "  << STAT_label[STATL_CUMULATIVE]
         << " " << STAT_label[STATL_DISTRIBUTION] << " " << STAT_label[STATL_FUNCTION];
    }
    os << endl;

    for (i = 0;i < nb_value;i++) {
      if (file_flag) {
        os << "# ";
      }
      os << setw(width[0]) << i;

      for (j = 0;j < nb_state;j++) {
        if (i < observation[j]->nb_value) {
          os << setw(width[1]) << observation[j]->mass[i];
        }
        else {
          os << setw(width[1]) << " ";
        }
      }

      for (j = 0;j < nb_state;j++) {
        if (i < observation[j]->nb_value) {
          os << setw(width[2]) << observation[j]->cumul[i];
        }
        else {
          os << setw(width[2]) << " ";
        }
      }
      os << endl;
    }

    os.setf(format_flags , ios::adjustfield);
  }

  if (marginal_distribution) {
    double likelihood , information;
    Test test(CHI2);


    if ((weight) && (mixture)) {
      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << STAT_label[STATL_MIXTURE] << " - "
         << STAT_label[STATL_THEORETICAL] << " " << STAT_label[STATL_WEIGHTS] << ":";

      for (i = 0;i < nb_state;i++) {
        os << " " << weight->mass[i];
      }
      os << endl;

      mixture->ascii_characteristic_print(os , false , file_flag);

      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << STAT_label[STATL_MARGINAL] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";
      marginal_distribution->ascii_characteristic_print(os , false , file_flag);

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
      os << STAT_label[STATL_MIXTURE] << " - " 
         << STAT_label[STATL_RESTORATION] << " " << STAT_label[STATL_WEIGHTS] << ":";

      for (i = 0;i < nb_state;i++) {
        os << " " << restoration_weight->mass[i];
      }
      os << endl;

      restoration_mixture->ascii_characteristic_print(os , false , file_flag);

      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << STAT_label[STATL_MARGINAL] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";
      marginal_distribution->ascii_characteristic_print(os , false , file_flag);

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

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a DiscreteParametricProcess object at the spreadsheet format.
 *
 *  \param[in,out] os                    stream,
 *  \param[in]     empirical_observation pointer on the observation frequency distributions,
 *  \param[in]     marginal_distribution pointer on the marginal frequency distribution,
 *  \param[in]     model                 model type.
 */
/*--------------------------------------------------------------*/

ostream& DiscreteParametricProcess::spreadsheet_print(ostream &os , FrequencyDistribution **empirical_observation ,
                                                      FrequencyDistribution *marginal_distribution ,
                                                      model_type model) const

{
  int i , j;
  double scale[NB_STATE];
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
    os << " " << i << "\t"  << STAT_word[STATW_OBSERVATION_DISTRIBUTION] << endl;
    observation[i]->spreadsheet_print(os);
    observation[i]->spreadsheet_parametric_characteristic_print(os);

    if (empirical_observation) {
      os << "\n";
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
      empirical_observation[i]->spreadsheet_characteristic_print(os);

      os << "\n";
      if (empirical_observation[i]->nb_element > 0) {
        os << "\t";
        switch (model) {
        case MIXTURE :
          os << STAT_label[STATL_COMPONENT];
          break;
        case HIDDEN_MARKOV :
          os << STAT_label[STATL_STATE];
          break;
        }
        os << " " << i << " " << STAT_label[STATL_OBSERVATION]
           << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
      }
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
      if (empirical_observation[i]->nb_element > 0) {
        os << "\t" << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " "
           << STAT_label[STATL_FUNCTION];
      }
      os << "\t" << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_DISTRIBUTION] << " "
         << STAT_label[STATL_FUNCTION] << endl;

      observation[i]->Distribution::spreadsheet_print(os , true , false , false ,
                                                      (empirical_observation[i]->nb_element > 0 ? empirical_observation[i] : NULL));
    }
  }

  if (!empirical_observation) {
    os << "\n";
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
      os << " " << i << " " << STAT_label[STATL_CUMULATIVE]
         << " " << STAT_label[STATL_DISTRIBUTION] << " " << STAT_label[STATL_FUNCTION];
    }
    os << endl;

    for (i = 0;i < nb_value;i++) {
      os << i;

      for (j = 0;j < nb_state;j++) {
       os << "\t";
        if (i < observation[j]->nb_value) {
          os << observation[j]->mass[i];
        }
      }

      for (j = 0;j < nb_state;j++) {
       os << "\t";
        if (i < observation[j]->nb_value) {
          os << observation[j]->cumul[i];
        }
      }
      os << endl;
    }
  }

  if (marginal_distribution) {
    double likelihood , information;
    Test test(CHI2);


    if ((weight) && (mixture)) {
      os << "\n" << STAT_label[STATL_MIXTURE] << "\t"
         << STAT_label[STATL_THEORETICAL] << " " << STAT_label[STATL_WEIGHTS];
      for (i = 0;i < nb_state;i++) {
        os << "\t" << weight->mass[i];
      }
      os << endl;

      mixture->spreadsheet_characteristic_print(os);

      os << "\n" << STAT_label[STATL_MARGINAL] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\t";
      marginal_distribution->spreadsheet_characteristic_print(os);

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
      os << "\n" << STAT_label[STATL_MIXTURE] << "\t"
         << STAT_label[STATL_RESTORATION] << " " << STAT_label[STATL_WEIGHTS];
      for (i = 0;i < nb_state;i++) {
        os << "\t" << restoration_weight->mass[i];
      }
      os << endl;

      restoration_mixture->spreadsheet_characteristic_print(os);

      os << "\n" << STAT_label[STATL_MARGINAL] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\t";
      marginal_distribution->spreadsheet_characteristic_print(os);

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

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Plot of a DiscreteParametricProcess object using Gnuplot.
 *
 *  \param[in] prefix                file prefix,
 *  \param[in] title                 figure title,
 *  \param[in] process               observation process index,
 *  \param[in] empirical_observation pointer on  the observation frequency distributions,
 *  \param[in] marginal_distribution pointer on the marginal frequency distribution,
 *  \param[in] model                 model type.
 *
 *  \return                          error status.
 */
/*--------------------------------------------------------------*/

bool DiscreteParametricProcess::plot_print(const char *prefix , const char *title , int process ,
                                           FrequencyDistribution **empirical_observation ,
                                           FrequencyDistribution *marginal_distribution ,
                                           model_type model) const

{
  bool status;
  int i , j , k , m;
  int nb_dist , *dist_nb_value;
  double max , *scale;
  const Distribution **pdist;
  const FrequencyDistribution **phisto;
  ostringstream data_file_name[3 * NB_STATE + 1];


  // writing of the data files

  nb_dist = 0;
  if (empirical_observation) {
    nb_dist += nb_state;
  }
  if (marginal_distribution) {
    if ((weight) && (mixture)) {
      nb_dist++;
    }
    if ((restoration_weight) && (restoration_mixture)) {
      nb_dist++;
    }
  }

  if (nb_dist > 0) {
    pdist = new const Distribution*[nb_dist];
    scale = new double[nb_dist + 1];
    phisto = new const FrequencyDistribution*[nb_dist];

    if (empirical_observation) {
      for (i = 0;i < nb_state;i++) {
        pdist[i] = observation[i];
        phisto[i] = empirical_observation[i];

        if (empirical_observation[i]->nb_element > 0) {
          scale[i] = phisto[i]->nb_element;
        }
        else {
          scale[i] = 1.;
        }
      }
    }

    else {
      i = 0;
    }

    if (marginal_distribution) {
      if ((weight) && (mixture)) {
        pdist[i] = mixture;
        phisto[i] = marginal_distribution;
        scale[i] = marginal_distribution->nb_element;
        i++;
      }
      if ((restoration_weight) && (restoration_mixture)) {
        pdist[i] = restoration_mixture;
        phisto[i] = marginal_distribution;
        scale[i] = marginal_distribution->nb_element;
      }
    }

    data_file_name[0] << prefix << process << ".dat";
    status = stat_tool::plot_print((data_file_name[0].str()).c_str() , nb_dist , pdist , scale ,
                                   NULL , nb_dist , phisto);
  }

  else {
    data_file_name[1] << prefix << process << 1 << ".dat";
    status = observation[0]->plot_print((data_file_name[1].str()).c_str());
  }

  if (status) {
    i = (nb_dist > 0 ? 1 : 2);
    if (!empirical_observation) {
      for (j = (nb_dist > 0 ? 0 : 1);j < nb_state;j++) {
        data_file_name[i] << prefix << process << i << ".dat";
        observation[j]->plot_print((data_file_name[i].str()).c_str());
        i++;
      }
    }

    if (marginal_distribution) {
      if ((weight) && (mixture)) {
        for (j = 0;j < nb_state;j++) {
          pdist[0] = observation[j];
          scale[nb_dist] = weight->mass[j] * marginal_distribution->nb_element;

          data_file_name[i] << prefix << process << i << ".dat";
          stat_tool::plot_print((data_file_name[i].str()).c_str() , 1 , pdist , scale + nb_dist ,
                                NULL , 0 , NULL);
          i++;
        }
      }

      if ((restoration_weight) && (restoration_mixture)) {
        for (j = 0;j < nb_state;j++) {
          pdist[0] = observation[j];
          scale[nb_dist] = restoration_weight->mass[j] * marginal_distribution->nb_element;

          data_file_name[i] << prefix << process << i << ".dat";
          stat_tool::plot_print((data_file_name[i].str()).c_str() , 1 , pdist , scale + nb_dist ,
                                NULL , 0 , NULL);
          i++;
        }
      }
    }

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
        file_name[1] << label(prefix) << process << ".ps";
        out_file << "set output \"" << file_name[1].str() << "\"\n\n";
      }

      out_file << "set border 15 lw 0\n" << "set tics out\n" << "set xtics nomirror\n"
               << "set title \"";
      if (title) {
        out_file << title << " - ";
      }
      out_file << STAT_label[STATL_OUTPUT_PROCESS] << " " << process << "\"\n\n";

      j = 1;

      if (empirical_observation) {
        for (k = 0;k < nb_state;k++) {
          if (observation[k]->nb_value - 1 < TIC_THRESHOLD) {
            out_file << "set xtics 0,1" << endl;
          }

          if (empirical_observation[k]->nb_element > 0) {
            out_file << "plot [0:" << observation[k]->nb_value - 1 << "] [0:"
                     << (int)(MAX(empirical_observation[k]->max , observation[k]->max * scale[k]) * YSCALE) + 1
                     << "] \"" << label((data_file_name[0].str()).c_str()) << "\" using " << k + 1
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
            out_file << "\"" << label((data_file_name[0].str()).c_str()) << "\" using " << nb_dist + k + 1
                     << " title \"";
            switch (model) {
            case MIXTURE :
              out_file << STAT_label[STATL_COMPONENT];
              break;
            case HIDDEN_MARKOV :
              out_file << STAT_label[STATL_STATE];
              break;
            }
            out_file << " " << k << " " << STAT_label[STATL_OBSERVATION] << " " << STAT_label[STATL_DISTRIBUTION];
            observation[k]->plot_title_print(out_file);
            out_file << "\" with linespoints" << endl;
          }

          else {
            out_file << "plot [0:" << observation[k]->nb_value - 1 << "] [0:"
                     << MIN(observation[k]->max * YSCALE , 1.) << "] \""
                     << label((data_file_name[0].str()).c_str()) << "\" using " << nb_dist + k + 1
                     << " title \"";
            switch (model) {
            case MIXTURE :
              out_file << STAT_label[STATL_COMPONENT];
              break;
            case HIDDEN_MARKOV :
              out_file << STAT_label[STATL_STATE];
              break;
            }
            out_file << " " << k << " " << STAT_label[STATL_OBSERVATION] << " " << STAT_label[STATL_DISTRIBUTION];
            observation[k]->plot_title_print(out_file);
            out_file << "\" with linespoints" << endl;
          }

          if (observation[k]->nb_value - 1 < TIC_THRESHOLD) {
            out_file << "set xtics autofreq" << endl;
          }

          if ((i == 0) && (k < nb_state - 1)) {
            out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
          }
          out_file << endl;
        }
      }

      else {
        if (nb_value - 1 < TIC_THRESHOLD) {
          out_file << "set xtics 0,1" << endl;
        }

        max = observation[0]->max;
        for (k = 1;k < nb_state;k++) {
          if (observation[k]->max > max) {
            max = observation[k]->max;
          }
        }

        out_file << "plot [0:" << nb_value - 1 << "] [0:"
                 << MIN(max * YSCALE , 1.) << "] ";

        for (k = 0;k < nb_state;k++) {
          out_file << "\"" << label((data_file_name[j].str()).c_str()) << "\" title \"";
          switch (model) {
          case MIXTURE :
            out_file << STAT_label[STATL_COMPONENT];
            break;
          case HIDDEN_MARKOV :
            out_file << STAT_label[STATL_STATE];
            break;
          }
          out_file << " " << k << " " << STAT_label[STATL_OBSERVATION] << " " << STAT_label[STATL_DISTRIBUTION];
          observation[k]->plot_title_print(out_file);
          out_file << "\" with linespoints";
          if (k < nb_state - 1) {
            out_file << ",\\";
          }
          out_file << endl;
          j++;
        }

        if (nb_value - 1 < TIC_THRESHOLD) {
          out_file << "set xtics autofreq" << endl;
        }
        k = 0;
      }

      if (marginal_distribution) {
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

          k++;
          out_file << "plot [0:" << nb_value - 1 << "] [0:"
                   << (int)(MAX(marginal_distribution->max , mixture->max * marginal_distribution->nb_element) * YSCALE) + 1
                   << "] \"" << label((data_file_name[0].str()).c_str()) << "\" using " << k
                   << " title \"" << STAT_label[STATL_MARGINAL] << " "
                   << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\" with impulses,\\" << endl;

          for (m = 0;m < nb_state;m++) {
            out_file << "\"" << label((data_file_name[j].str()).c_str()) << "\" title \"";
            switch (model) {
            case MIXTURE :
              out_file << STAT_label[STATL_COMPONENT];
              break;
            case HIDDEN_MARKOV :
              out_file << STAT_label[STATL_STATE];
              break;
            }
            out_file << " " << m << " " << STAT_label[STATL_OBSERVATION] << " " << STAT_label[STATL_DISTRIBUTION];
            observation[m]->plot_title_print(out_file);
            out_file << "\" with linespoints,\\" << endl;
            j++;
          }

          out_file << "\"" << label((data_file_name[0].str()).c_str()) << "\" using " << nb_dist + k
                   << " title \"" << STAT_label[STATL_MIXTURE] << "\" with linespoints" << endl;

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

          k++;
          out_file << "plot [0:" << nb_value - 1 << "] [0:"
                   << (int)(MAX(marginal_distribution->max , restoration_mixture->max * marginal_distribution->nb_element) * YSCALE) + 1
                   << "] \"" << label((data_file_name[0].str()).c_str()) << "\" using " << k
                   << " title \"" << STAT_label[STATL_MARGINAL] << " "
                   << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\" with impulses,\\" << endl;

          for (m = 0;m < nb_state;m++) {
            out_file << "\"" << label((data_file_name[j].str()).c_str()) << "\" title \"";
            switch (model) {
            case MIXTURE :
              out_file << STAT_label[STATL_COMPONENT];
              break;
            case HIDDEN_MARKOV :
              out_file << STAT_label[STATL_STATE];
              break;
            }
            out_file << " " << m << " " << STAT_label[STATL_OBSERVATION] << " " << STAT_label[STATL_DISTRIBUTION];
            observation[m]->plot_title_print(out_file);
            out_file << "\" with linespoints,\\" << endl;
            j++;
          }

          out_file << "\"" << label((data_file_name[0].str()).c_str()) << "\" using " << nb_dist + k
                   << " title \"" << STAT_label[STATL_MIXTURE] << "\" with linespoints" << endl;

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

  if (nb_dist > 0) {
    delete [] pdist;
    delete [] scale;
    delete [] phisto;
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Plot of a DiscreteParametricProcess object.
 *
 *  \param[in] plot                  reference on a MultiPlotSet object,
 *  \param[in] index                 MultiPlot index,
 *  \param[in] process               observation process index,
 *  \param[in] empirical_observation pointer on the observation frequency distributions,
 *  \param[in] marginal_distribution pointer on the marginal frequency distribution,
 *  \param[in] model                 model type.
 */
/*--------------------------------------------------------------*/

void DiscreteParametricProcess::plotable_write(MultiPlotSet &plot , int &index , int process ,
                                               FrequencyDistribution **empirical_observation ,
                                               FrequencyDistribution *marginal_distribution ,
                                               model_type model) const

{
  int i , j;
  double scale , max;
  ostringstream title , legend;


  plot.variable_nb_viewpoint[process] = 1;

  if (empirical_observation) {
    for (i = 0;i < nb_state;i++) {

      // fit of observation distribution

      plot.variable[index] = process;
//      plot.viewpoint[index] = OBSERVATION;

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
        legend << " " << i << " " << STAT_label[STATL_OBSERVATION] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
        plot[index][0].legend = legend.str();

        plot[index][0].style = "impulses";

        empirical_observation[i]->plotable_frequency_write(plot[index][0]);
        j = 1;
      }

      else {
        scale = 1.;
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
      observation[i]->plot_title_print(legend);
      plot[index][j].legend = legend.str();

      plot[index][j].style = "linespoints";

      observation[i]->plotable_mass_write(plot[index][j] , scale);
      index++;
    }
  }

  else {

    // observation distributions

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
      observation[i]->plot_title_print(legend);
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

      // cumulative distribution functions

      title.str("");
      title << STAT_label[STATL_OUTPUT_PROCESS] << " " << process << " - "
            << STAT_label[STATL_RESTORATION] << " " << STAT_label[STATL_WEIGHTS];
      plot[index].title = title.str();

      plot[index].xrange = Range(0 , nb_value);
      plot[index].yrange = Range(0. , 1.);

      if (nb_value - 1 < TIC_THRESHOLD) {
        plot[index].xtics = 1;
      }

      plot[index].resize(2);

      legend.str("");
      legend << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_MARGINAL] << " "
             << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " " << STAT_label[STATL_FUNCTION];
      plot[index][0].legend = legend.str();

      plot[index][0].style = "linespoints";

      marginal_distribution->plotable_cumul_write(plot[index][0]);

      legend.str("");
      legend << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_MIXTURE]
             << " " << STAT_label[STATL_FUNCTION];
      plot[index][1].legend = legend.str();

      plot[index][1].style = "linespoints";

      mixture->plotable_cumul_write(plot[index][1]);

      index++;
    }

    if ((restoration_weight) && (restoration_mixture)) {

      // fit of the mixture of observation distributions (restoration weights)

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

      // cumulative distribution functions

      title.str("");
      if (process > 0) {
        title << STAT_label[STATL_OUTPUT_PROCESS] << " " << process;
      }
      title << STAT_label[STATL_RESTORATION] << " " << STAT_label[STATL_WEIGHTS];
      plot[index].title = title.str();

      plot[index].xrange = Range(0 , nb_value);
      plot[index].yrange = Range(0. , 1.);

      if (nb_value - 1 < TIC_THRESHOLD) {
        plot[index].xtics = 1;
      }

      plot[index].resize(2);

      legend.str("");
      legend << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_MARGINAL] << " "
             << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " " << STAT_label[STATL_FUNCTION];
      plot[index][0].legend = legend.str();

      plot[index][0].style = "linespoints";

      marginal_distribution->plotable_cumul_write(plot[index][0]);

      legend.str("");
      legend << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_MIXTURE]
             << " " << STAT_label[STATL_FUNCTION];
      plot[index][1].legend = legend.str();

      plot[index][1].style = "linespoints";

      restoration_mixture->plotable_cumul_write(plot[index][1]);

      index++;
    }
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the support of the observation process.
 */
/*--------------------------------------------------------------*/

void DiscreteParametricProcess::nb_value_computation()

{
  int i;


  nb_value = observation[0]->nb_value;
  for (i = 0;i < nb_state;i++) {
    if (observation[i]->nb_value > nb_value) {
      nb_value = observation[i]->nb_value;
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

void DiscreteParametricProcess::state_permutation(int *permut) const

{
  int i;
  DiscreteParametric **pobservation;


  pobservation = new DiscreteParametric*[nb_state];

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
 *  \brief Computation of the number of free parameters of a discrete
 *         parametric observation process.
 *
 *  \return number of free parameters.
 */
/*--------------------------------------------------------------*/

int DiscreteParametricProcess::nb_parameter_computation() const

{
  int i;
  int nb_parameter = 0;


  for (i = 0;i < nb_state;i++) {
    nb_parameter += observation[i]->nb_parameter_computation();
    if (observation[i]->inf_bound == 0) {
      nb_parameter--;
    }
  }

  return nb_parameter;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the mean of a mixture of discrete parametric 
 *         observation distributions.
 *
 *  \param[in] pweight pointer on the weight distribution.
 *
 *  \return            mixture mean.
 */
/*--------------------------------------------------------------*/

double DiscreteParametricProcess::mean_computation(Distribution *pweight) const

{
  int i;
  double mean;


  mean = 0.;
  for (i = 0;i < nb_state;i++) {
    mean += pweight->mass[i] * observation[i]->parametric_mean_computation();
  }

  return mean;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the variance of a mixture of discrete parametric 
 *         observation distributions.
 *
 *  \param[in] pweight pointer on the weight distribution,
 *  \param[in] mean    mean.
 *
 *  \return            mixture variance.
 */
/*--------------------------------------------------------------*/

double DiscreteParametricProcess::variance_computation(Distribution *pweight , double mean) const

{
  int i;
  double variance , component_mean;


  if (mean == D_INF) {
    mean = mean_computation(pweight);
  }

  variance = -mean * mean;
  for (i = 0;i < nb_state;i++) {
    component_mean = observation[i]->parametric_mean_computation();
    variance += pweight->mass[i] * (observation[i]->parametric_variance_computation() +
                                    component_mean * component_mean);
  }

  return variance;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of a mixture of discrete parametric observation distributions.
 *
 *  \param[in] pweight pointer on the weight distribution.
 *
 *  \return            mixture of discrete parametric observation distributions.
 */
/*--------------------------------------------------------------*/

Distribution* DiscreteParametricProcess::mixture_computation(Distribution *pweight)

{
  int i , j;
  Distribution *mixture;


  mixture = new Distribution(nb_value);

  mixture->offset = observation[0]->offset;
  for (i = 1;i < nb_state;i++) {
    if (observation[i]->offset < mixture->offset) {
      mixture->offset = observation[i]->offset;
    }
  }

  for (i = 0;i < mixture->offset;i++) {
    mixture->mass[i] = 0.;
  }
  for (i = mixture->offset;i < mixture->nb_value;i++) {
    mixture->mass[i] = 0.;
    for (j = 0;j < nb_state;j++) {
      if ((i >= observation[j]->offset) && (i < observation[j]->nb_value)) {
        mixture->mass[i] += pweight->mass[j] * observation[j]->mass[i];
      }
    }
  }

  mixture->cumul_computation();

  mixture->max_computation();
//  mixture->mean_computation();
//  mixture->variance_computation();
  mixture->mean = mean_computation(pweight);
  mixture->variance = variance_computation(pweight , mixture->mean);

  return mixture;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Initialization of discrete parametric observation distributions
 *         using perturbation.
 */
/*--------------------------------------------------------------*/

void DiscreteParametricProcess::init()

{
  int i , j;
  double noise_proba , *pmass;


  for (i = 0;i < nb_state;i++) {
    noise_proba = NOISE_PROBABILITY * nb_state / nb_value;

    pmass = observation[i]->mass;
    for (j = 0;j < i * nb_value / nb_state;j++) {
      *pmass++ -= noise_proba / (nb_state - 1);
    }
    for (j = i * nb_value / nb_state;j < (i + 1) * nb_value / nb_state;j++) {
      *pmass++ += noise_proba;
    }
    for (j = (i + 1) * nb_value / nb_state;j < nb_value;j++) {
      *pmass++ -= noise_proba / (nb_state - 1);
    }
  }
}


};  // namespace stat_tool
