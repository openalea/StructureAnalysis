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
 *  \brief Constructor of the Chain class.
 *
 *  \param[in] itype     type,
 *  \param[in] inb_state number of states,
 *  \param[in] init_flag flag initialization.
 */
/*--------------------------------------------------------------*/

Chain::Chain(process_type itype , int inb_state , bool init_flag)

{
  type = itype;
  nb_state = inb_state;

  accessibility = NULL;

  nb_component = 0;
  component_nb_state = NULL;
  component = NULL;
  stype = NULL;

  if (nb_state == 0) {
    nb_row = 0;

    initial = NULL;
    transition = NULL;
  }

  else {
    int i , j;


    nb_row = nb_state;

    initial = new double[type == ORDINARY ? nb_state : nb_row];

    if (init_flag) {
      for (i = 0;i < (type == ORDINARY ? nb_state : nb_row);i++) {
        initial[i] = 0.;
      }
    }

    transition = new double*[nb_row];
    for (i = 0;i < nb_row;i++) {
      transition[i] = new double[nb_state];

      if (init_flag) {
        for (j = 0;j < nb_state;j++) {
          transition[i][j] = 0.;
        }
      }
    }
  }

  cumul_initial = NULL;
  cumul_transition = NULL;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Constructor of the Chain class.
 *
 *  \param[in] itype     type,
 *  \param[in] inb_state number of states,
 *  \param[in] inb_row   number of rows of the transition probability matrix,
 *  \param[in] init_flag flag initialization.
 */
/*--------------------------------------------------------------*/

Chain::Chain(process_type itype , int inb_state , int inb_row , bool init_flag)

{
  int i , j;


  type = itype;
  nb_state = inb_state;
  nb_row = inb_row;

  accessibility = NULL;

  nb_component = 0;
  component_nb_state = NULL;
  component = NULL;
  stype = NULL;

  initial = new double[type == ORDINARY ? nb_state : nb_row];

  if (init_flag) {
    for (i = 0;i < (type == ORDINARY ? nb_state : nb_row);i++) {
      initial[i] = 0.;
    }
  }

  transition = new double*[nb_row];
  for (i = 0;i < nb_row;i++) {
    transition[i] = new double[nb_state];

    if (init_flag) {
      for (j = 0;j < nb_state;j++) {
        transition[i][j] = 0.;
      }
    }
  }

  cumul_initial = NULL;
  cumul_transition = NULL;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Copy of Markov chain parameters.
 *
 *  \param[in] chain reference on a Chain object.
 */
/*--------------------------------------------------------------*/

void Chain::parameter_copy(const Chain &chain)

{
  int i , j;


  for (i = 0;i < (type == ORDINARY ? nb_state : nb_row);i++) {
    initial[i] = chain.initial[i];
  }

  if (chain.cumul_initial) {
    for (i = 0;i < (type == ORDINARY ? nb_state : nb_row);i++) {
      cumul_initial[i] = chain.cumul_initial[i];
    }
  }

  for (i = 0;i < nb_row;i++) {
    for (j = 0;j < nb_state;j++) {
      transition[i][j] = chain.transition[i][j];
    }
  }

  if (chain.cumul_transition) {
    for (i = 0;i < nb_row;i++) {
      for (j = 0;j < nb_state;j++) {
        cumul_transition[i][j] = chain.cumul_transition[i][j];
      }
    }
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Copy of a Chain object.
 *
 *  \param[in] chain reference on a Chain object.
 */
/*--------------------------------------------------------------*/

void Chain::copy(const Chain &chain)

{
  int i , j;


  type = chain.type;
  nb_state = chain.nb_state;
  nb_row = chain.nb_row;

  if (chain.nb_component > 0) {
    accessibility = new bool*[nb_state];
    for (i = 0;i < nb_state;i++) {
      accessibility[i] = new bool[nb_state];
      for (j = 0;j < nb_state;j++) {
        accessibility[i][j] = chain.accessibility[i][j];
      }
    }

    nb_component = chain.nb_component;
    component_nb_state = new int[nb_component];
    component = new int*[nb_component];

    for (i = 0;i < nb_component;i++) {
      component_nb_state[i] = chain.component_nb_state[i];
      component[i] = new int[component_nb_state[i]];
      for (j = 0;j < component_nb_state[i];j++) {
        component[i][j] = chain.component[i][j];
      }
    }

    stype = new state_type[nb_state];
    for (i = 0;i < nb_state;i++) {
      stype[i] = chain.stype[i];
    }
  }

  else {
    accessibility = NULL;
    nb_component = 0;
    component_nb_state = NULL;
    component = NULL;
    stype = NULL;
  }

  initial = new double[type == ORDINARY ? nb_state : nb_row];

  if (chain.cumul_initial) {
    cumul_initial = new double[type == ORDINARY ? nb_state : nb_row];
  }
  else {
    cumul_initial = NULL;
  }

  transition = new double*[nb_row];
  for (i = 0;i < nb_row;i++) {
    transition[i] = new double[nb_state];
  }

  if (chain.cumul_transition) {
    cumul_transition = new double*[nb_row];
    for (i = 0;i < nb_row;i++) {
      cumul_transition[i] = new double[nb_state];
    }
  }
  else {
    cumul_transition = NULL;
  }

  parameter_copy(chain);
}


/*--------------------------------------------------------------*/
/**
 *  \brief Destruction of the data members of a Chain object.
 */
/*--------------------------------------------------------------*/

void Chain::remove()

{
  int i;


  if (accessibility) {
    for (i = 0;i < nb_state;i++) {
      delete [] accessibility[i];
    }
    delete [] accessibility;
  }

  delete [] component_nb_state;

  if (component) {
    for (i = 0;i < nb_component;i++) {
      delete [] component[i];
    }
    delete [] component;
  }

  delete [] stype;

  delete [] initial;
  delete [] cumul_initial;

  if (transition) {
    for (i = 0;i < nb_row;i++) {
      delete [] transition[i];
    }
    delete [] transition;
  }

  if (cumul_transition) {
    for (i = 0;i < nb_row;i++) {
      delete [] cumul_transition[i];
    }
    delete [] cumul_transition;
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Destructor of the Chain class.
 */
/*--------------------------------------------------------------*/

Chain::~Chain()

{
  remove();
}


/*--------------------------------------------------------------*/
/**
 *  \brief Assignment operator of the Chain class.
 *
 *  \param[in] chain reference on a Chain object.
 *
 *  \return          Chain object.
 */
/*--------------------------------------------------------------*/

Chain& Chain::operator=(const Chain &chain)

{
  if (&chain != this) {
    remove();
    copy(chain);
  }

  return *this;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Analysis of the format of a Chain object.
 *
 *  \param[in] error   reference on a StatError object,
 *  \param[in] in_file stream,
 *  \param[in] line    reference on the file line index,
 *  \param[in] type    process type (ORDINARY/EQUILIBRIUM).
 *
 *  \return            Chain object.
 */
/*--------------------------------------------------------------*/

Chain* Chain::parsing(StatError &error , ifstream &in_file , int &line , process_type type)

{
  string buffer;
  size_t position;
  typedef tokenizer<char_separator<char>> tokenizer;
  char_separator<char> separator(" \t");
  bool status = true , lstatus , **logic_transition;
  int i;
  int read_line , nb_state = 0 , value;
  double proba , cumul;
  Chain *chain;


  chain = NULL;

  // analysis of the lines defining the number of states and the order (or memory length)

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

      // test number of states

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
          if ((value < 2) || (value > NB_STATE)) {
            lstatus = false;
          }
          else {
            nb_state = value;
          }
        }

        if (!lstatus) {
          status = false;
          error.update(STAT_parsing[STATP_NB_STATE] , line , i + 1);
        }
        break;
      }

      // test STATES keyword

      case 1 : {
        if (*token != STAT_word[STATW_STATES]) {
          status = false;
          error.correction_update(STAT_parsing[STATP_KEYWORD] , STAT_word[STATW_STATES] , line , i + 1);
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

  if (nb_state == 0) {
    status = false;
    error.update(STAT_parsing[STATP_FORMAT] , line);
  }

  if (status) {
    chain = new Chain(type , nb_state , nb_state , false);

    // analysis initial/transition probabilities

    read_line = 0;
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

      if ((read_line == 0) || ((type == ORDINARY) && (read_line == 2))) {
        for (tokenizer::iterator token = tok_buffer.begin();token != tok_buffer.end();token++) {

          // test INITIAL_PROBABILITIES/TRANSITION_PROBABILITIES keyword

          if (i == 0) {
            if ((type == ORDINARY) && (read_line == 0)) {
              if (*token != STAT_word[STATW_INITIAL_PROBABILITIES]) {
                status = false;
                error.correction_update(STAT_parsing[STATP_KEYWORD] , STAT_word[STATW_INITIAL_PROBABILITIES] , line);
              }
            }

            else {
              if (*token != STAT_word[STATW_TRANSITION_PROBABILITIES]) {
                status = false;
                error.correction_update(STAT_parsing[STATP_KEYWORD] , STAT_word[STATW_TRANSITION_PROBABILITIES] , line);
              }
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
          if (i < nb_state) {
            lstatus = true;

/*            try {
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

                if ((type == ORDINARY) && (read_line == 1)) {
                  chain->initial[i] = proba;
                }
                else {
                  chain->transition[read_line - (type == ORDINARY ? 3 : 1)][i] = proba;
                }
              }
            }

            if (!lstatus) {
              status = false;
              if ((type == ORDINARY) && (read_line == 1)) {
                error.update(STAT_parsing[STATP_INITIAL_PROBABILITY] , line , i + 1);
              }
              else {
                error.update(STAT_parsing[STATP_TRANSITION_PROBABILITY] , line , i + 1);
              }
            }
          }

          i++;
        }

        if (i > 0) {
          if (i != nb_state) {
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
        if (((type == ORDINARY) && (read_line == chain->nb_row + 3)) ||
            ((type == EQUILIBRIUM) && (read_line == chain->nb_row + 1))) {
          break;
        }
      }
    }

    if (status) {

      // test accessible states

      status = chain->strongly_connected_component_research(error);

      logic_transition = chain->logic_transition_computation();
      status = chain->strongly_connected_component_research(error , logic_transition);

      if (status) {

        // test irreducibility in the equilibrium process case

        chain->component_computation(logic_transition);
        if ((type == EQUILIBRIUM) && (chain->nb_component > 1)) {
          status = false;
          error.correction_update(STAT_parsing[STATP_CHAIN_STRUCTURE] , STAT_parsing[STATP_IRREDUCIBLE]);
        }
      }

      for (i = 0;i < nb_state;i++) {
        delete [] logic_transition[i];
      }
      delete [] logic_transition;
    }

    if (!status) {
      delete chain;
      chain = NULL;
    }
  }

  return chain;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a Chain object.
 *
 *  \param[in,out] os        stream,
 *  \param[in]     file_flag file flag.
 */
/*--------------------------------------------------------------*/

ostream& Chain::ascii_print(ostream &os , bool file_flag) const

{
  int i , j;
  int buff , width;
  ios_base::fmtflags format_flags;


  format_flags = os.setf(ios::left , ios::adjustfield);

  os << "\n" << nb_state << " " << STAT_word[STATW_STATES] << endl;

  // computation of the column width

  width = column_width(nb_state , initial);
  for (i = 0;i < nb_row;i++) {
    buff = column_width(nb_state , transition[i]);
    if (buff > width) {
      width = buff;
    }
  }
  width += ASCII_SPACE;

  os << "\n";

  switch (type) {

  case ORDINARY : {
    os << STAT_word[STATW_INITIAL_PROBABILITIES] << endl;
    break;
  }

  case EQUILIBRIUM : {
    if (file_flag) {
      os << "# ";
    }
    os << STAT_label[STATL_STATIONARY_PROBABILITIES] << endl;
    break;
  }
  }

  if ((type == EQUILIBRIUM) && (file_flag)) {
    os << "# ";
  }
  for (i = 0;i < nb_state;i++) {
    os << setw(width) << initial[i];
  }
  os << endl;

  os << "\n" << STAT_word[STATW_TRANSITION_PROBABILITIES] << endl;

  for (i = 0;i < nb_row;i++) {
    for (j = 0;j < nb_state;j++) {
      os << setw(width) << transition[i][j];
    }
    os << endl;
  }

  if (nb_component > 0) {
    for (i = 0;i < nb_component;i++) {
      os << "\n";
      if (file_flag) {
        os << "# ";
      }

      switch (stype[component[i][0]]) {
      case TRANSIENT :
        os << STAT_label[STATL_TRANSIENT] << " ";
        break;
      default :
        os << STAT_label[STATL_RECURRENT] << " ";
        break;
      }
      os << STAT_label[STATL_CLASS] << ": " << STAT_label[component_nb_state[i] == 1 ? STATL_STATE : STATL_STATES];

      for (j = 0;j < component_nb_state[i];j++) {
        os << " " << component[i][j];
      }

      if (stype[component[i][0]] == ABSORBING) {
        os << " (" << STAT_label[STATL_ABSORBING] << " " << STAT_label[STATL_STATE] << ")";
      }
    }
    os << endl;
  }

  os.setf(format_flags , ios::adjustfield);

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a Chain object at the spreadsheet format.
 *
 *  \param[in,out] os stream.
 */
/*--------------------------------------------------------------*/

ostream& Chain::spreadsheet_print(ostream &os) const

{
  int i , j;


  os << "\n" << nb_state << "\t" << STAT_word[STATW_STATES] << endl;

  switch (type) {
  case ORDINARY :
    os << "\n" << STAT_word[STATW_INITIAL_PROBABILITIES] << endl;
    break;
  case EQUILIBRIUM :
    os << "\n" << STAT_label[STATL_STATIONARY_PROBABILITIES] << endl;
    break;
  }

  for (i = 0;i < nb_state;i++) {
    os << initial[i] << "\t";
  }
  os << endl;

  os << "\n" << STAT_word[STATW_TRANSITION_PROBABILITIES] << endl;

  for (i = 0;i < nb_row;i++) {
    for (j = 0;j < nb_state;j++) {
      os << transition[i][j] << "\t";
    }
    os << endl;
  }

  if (nb_component > 0) {
    for (i = 0;i < nb_component;i++) {
      switch (stype[component[i][0]]) {
      case TRANSIENT :
        os << "\n"<< STAT_label[STATL_TRANSIENT] << " ";
        break;
      default :
        os << "\n" << STAT_label[STATL_RECURRENT] << " ";
        break;
      }
      os << STAT_label[STATL_CLASS] << "\t" << STAT_label[component_nb_state[i] == 1 ? STATL_STATE : STATL_STATES];

      for (j = 0;j < component_nb_state[i];j++) {
        os << "\t" << component[i][j];
      }

      if (stype[component[i][0]] == ABSORBING) {
        os << " (" << STAT_label[STATL_ABSORBING] << " " << STAT_label[STATL_STATE] << ")";
      }
    }
    os << endl;
  }

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Construction of the cumulative initial and transition distribution functions of a Chain object.
 */
/*--------------------------------------------------------------*/

void Chain::create_cumul()

{
  int i;


  if (!cumul_initial) {
    cumul_initial = new double[type == ORDINARY ? nb_state : nb_row];
  }

  if (!cumul_transition) {
    cumul_transition = new double*[nb_row];
    for (i = 0;i < nb_row;i++) {
      cumul_transition[i] = new double[nb_state];
    }
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Destruction of the cumulative initial and transition distribution functions of a Chain object.
 */
/*--------------------------------------------------------------*/

void Chain::remove_cumul()

{
  int i;


  if (cumul_initial) {
    delete [] cumul_initial;

    cumul_initial = NULL;
  }

  if (cumul_transition) {
    for (i = 0;i < nb_row;i++) {
      delete [] cumul_transition[i];
    }
    delete [] cumul_transition;

    cumul_transition = NULL;
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of cumulative initial and transition distribution functions.
 */
/*--------------------------------------------------------------*/

void Chain::cumul_computation()

{
  if ((cumul_initial) && (cumul_transition)) {
    int i;


    stat_tool::cumul_computation((type == ORDINARY ? nb_state : nb_row) , initial , cumul_initial);

    for (i = 0;i < nb_row;i++) {
      stat_tool::cumul_computation(nb_state , transition[i] , cumul_transition[i]);
    }
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Log transform of Markov chain parameters.
 */
/*--------------------------------------------------------------*/

void Chain::log_computation()

{
  if ((cumul_initial) && (cumul_transition)) {
    int i;


    stat_tool::log_computation((type == ORDINARY ? nb_state : nb_row) , initial , cumul_initial);

    for (i = 0;i < nb_row;i++) {
      stat_tool::log_computation(nb_state , transition[i] , cumul_transition[i]);
    }
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Constructor by copy of the ChainData class.
 *
 *  \param[in] chain_data reference on a ChainData object.
 */
/*--------------------------------------------------------------*/

ChainData::ChainData(const ChainData &chain_data)

{
  copy(chain_data);
}


};  // namespace stat_tool
