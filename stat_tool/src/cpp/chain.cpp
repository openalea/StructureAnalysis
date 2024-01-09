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
 *       $Id: chain.cpp 17973 2015-04-23 06:33:24Z guedon $
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

#include "tool/rw_tokenizer.h"
#include "tool/rw_cstring.h"
#include "tool/rw_locale.h"
#include "tool/config.h"

#include "stat_tools.h"
#include "markovian.h"
#include "stat_label.h"

using namespace std;


namespace stat_tool {



/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Chain.
 *
 *  arguments : type, nombre d'etats, flag initialisation.
 *
 *--------------------------------------------------------------*/

Chain::Chain(char itype , int inb_state , bool init_flag)

{
  type = itype;
  nb_state = inb_state;

  accessibility = NULL;

  nb_component = 0;
  component_nb_state = NULL;
  component = NULL;
  state_type = NULL;

  if (nb_state == 0) {
    nb_row = 0;

    initial = NULL;
    transition = NULL;
  }

  else {
    register int i , j;


    nb_row = nb_state;

    initial = new double[type == 'o' ? nb_state : nb_row];

    if (init_flag) {
      for (i = 0;i < (type == 'o' ? nb_state : nb_row);i++) {
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


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Chain.
 *
 *  arguments : type, nombre d'etats, nombre de lignes de la matrice
 *              des probabilites de transition, flag initialisation.
 *
 *--------------------------------------------------------------*/

Chain::Chain(char itype , int inb_state , int inb_row , bool init_flag)

{
  register int i , j;


  type = itype;
  nb_state = inb_state;
  nb_row = inb_row;

  accessibility = NULL;

  nb_component = 0;
  component_nb_state = NULL;
  component = NULL;
  state_type = NULL;

  initial = new double[type == 'o' ? nb_state : nb_row];

  if (init_flag) {
    for (i = 0;i < (type == 'o' ? nb_state : nb_row);i++) {
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


/*--------------------------------------------------------------*
 *
 *  Copie des parametres d'une chaine de Markov.
 *
 *  argument : reference sur un objet Chain.
 *
 *--------------------------------------------------------------*/

void Chain::parameter_copy(const Chain &chain)

{
  register int i , j;


  for (i = 0;i < (type == 'o' ? nb_state : nb_row);i++) {
    initial[i] = chain.initial[i];
  }

  if (chain.cumul_initial) {
    for (i = 0;i < (type == 'o' ? nb_state : nb_row);i++) {
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


/*--------------------------------------------------------------*
 *
 *  Copie d'un objet Chain.
 *
 *  argument : reference sur un objet Chain.
 *
 *--------------------------------------------------------------*/

void Chain::copy(const Chain &chain)

{
  register int i , j;


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

    state_type = new char[nb_state];
    for (i = 0;i < nb_state;i++) {
      state_type[i] = chain.state_type[i];
    }
  }

  else {
    accessibility = NULL;
    nb_component = 0;
    component_nb_state = NULL;
    component = NULL;
    state_type = NULL;
  }

  initial = new double[type == 'o' ? nb_state : nb_row];

  if (chain.cumul_initial) {
    cumul_initial = new double[type == 'o' ? nb_state : nb_row];
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


/*--------------------------------------------------------------*
 *
 *  Destruction des champs d'un objet Chain.
 *
 *--------------------------------------------------------------*/

void Chain::remove()

{
  register int i;


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

  delete [] state_type;

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


/*--------------------------------------------------------------*
 *
 *  Destructeur de la classe Chain.
 *
 *--------------------------------------------------------------*/

Chain::~Chain()

{
  remove();
}


/*--------------------------------------------------------------*
 *
 *  Operateur d'assignement de la classe Chain.
 *
 *  argument : reference sur un objet Chain.
 *
 *--------------------------------------------------------------*/

Chain& Chain::operator=(const Chain &chain)

{
  if (&chain != this) {
    remove();
    copy(chain);
  }

  return *this;
}


/*--------------------------------------------------------------*
 *
 *  Analyse du format d'un objet Chain.
 *
 *  arguments : reference sur un objet StatError, stream,
 *              reference sur l'indice de la ligne lue, type du processus
 *              ('o' : ordinaire, 'e' : en equilibre).
 *
 *--------------------------------------------------------------*/

Chain* chain_parsing(StatError &error , ifstream &in_file , int &line , char type)

{
  RWLocaleSnapshot locale("en");
  RWCString buffer , token;
  size_t position;
  bool status = true , lstatus , **logic_transition;
  register int i;
  int read_line , nb_state = 0;
  long value;
  double proba , cumul;
  Chain *chain;


  chain = NULL;

  // analyse lignes definissant le nombre d'etats et l'ordre

  while (buffer.readLine(in_file , false)) {
    line++;

#   ifdef DEBUG
    cout << line << "  " << buffer << endl;
#   endif

    position = buffer.first('#');
    if (position != RW_NPOS) {
      buffer.remove(position);
    }
    i = 0;

    RWCTokenizer next(buffer);

    while (!((token = next()).isNull())) {
      switch (i) {

      // test valeur nombre d'etats

      case 0 : {
        lstatus = locale.stringToNum(token , &value);
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

      // test mot cle STATES

      case 1 : {
        if (token != STAT_word[STATW_STATES]) {
          status = false;
          error.correction_update(STAT_parsing[STATP_KEY_WORD] , STAT_word[STATW_STATES] , line , i + 1);
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

    // analyse probabilites initiales / probabilites de transition

    read_line = 0;
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

      if ((read_line == 0) || ((type == 'o') && (read_line == 2))) {
        while (!((token = next()).isNull())) {

          // test mot cle INITIAL_PROBABILITIES / TRANSITION_PROBABILITIES

          if (i == 0) {
            if ((type == 'o') && (read_line == 0)) {
              if (token != STAT_word[STATW_INITIAL_PROBABILITIES]) {
                status = false;
                error.correction_update(STAT_parsing[STATP_KEY_WORD] , STAT_word[STATW_INITIAL_PROBABILITIES] , line);
              }
            }

            else {
              if (token != STAT_word[STATW_TRANSITION_PROBABILITIES]) {
                status = false;
                error.correction_update(STAT_parsing[STATP_KEY_WORD] , STAT_word[STATW_TRANSITION_PROBABILITIES] , line);
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

        while (!((token = next()).isNull())) {
          if (i < nb_state) {
            lstatus = locale.stringToNum(token , &proba);
            if (lstatus) {
              if ((proba < 0.) || (proba > 1. - cumul + DOUBLE_ERROR)) {
                lstatus = false;
              }

              else {
                cumul += proba;

                if ((type == 'o') && (read_line == 1)) {
                  chain->initial[i] = proba;
                }
                else {
                  chain->transition[read_line - (type == 'o' ? 3 : 1)][i] = proba;
                }
              }
            }

            if (!lstatus) {
              status = false;
              if ((type == 'o') && (read_line == 1)) {
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
        if (((type == 'o') && (read_line == chain->nb_row + 3)) ||
            ((type == 'e') && (read_line == chain->nb_row + 1))) {
          break;
        }
      }
    }

    if (status) {

      // test etats atteignables

      status = chain->connex_component_research(error);

      logic_transition = chain->logic_transition_computation();
      status = chain->connex_component_research(error , logic_transition);

      if (status) {

        // test irreductibilite dans le cas en equilibre

        chain->component_computation(logic_transition);
        if ((type == 'e') && (chain->nb_component > 1)) {
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


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Chain.
 *
 *  arguments : stream, flag fichier.
 *
 *--------------------------------------------------------------*/

ostream& Chain::ascii_print(ostream &os , bool file_flag) const

{
  register int i , j;
  int buff , width;
  long old_adjust;


  old_adjust = os.setf(ios::left , ios::adjustfield);

  os << "\n" << nb_state << " " << STAT_word[STATW_STATES] << endl;

  // calcul des largeurs des colonnes

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

  case 'o' : {
    os << STAT_word[STATW_INITIAL_PROBABILITIES] << endl;
    break;
  }

  case 'e' : {
    if (file_flag) {
      os << "# ";
    }
    os << STAT_label[STATL_STATIONARY_PROBABILITIES] << endl;
    break;
  }
  }

  if ((type == 'e') && (file_flag)) {
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

      switch (state_type[component[i][0]]) {
      case 't' :
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

      if (state_type[component[i][0]] == 'a') {
        os << " (" << STAT_label[STATL_ABSORBING] << " " << STAT_label[STATL_STATE] << ")";
      }
    }
    os << endl;
  }

  os.setf((FMTFLAGS)old_adjust , ios::adjustfield);

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Chain au format tableur.
 *
 *  argument : stream.
 *
 *--------------------------------------------------------------*/

ostream& Chain::spreadsheet_print(ostream &os) const

{
  register int i , j;


  os << "\n" << nb_state << "\t" << STAT_word[STATW_STATES] << endl;

  switch (type) {
  case 'o' :
    os << "\n" << STAT_word[STATW_INITIAL_PROBABILITIES] << endl;
    break;
  case 'e' :
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
      switch (state_type[component[i][0]]) {
      case 't' :
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

      if (state_type[component[i][0]] == 'a') {
        os << " (" << STAT_label[STATL_ABSORBING] << " " << STAT_label[STATL_STATE] << ")";
      }
    }
    os << endl;
  }

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Creation des champs de type "cumul" d'un objet Chain.
 *
 *--------------------------------------------------------------*/

void Chain::create_cumul()

{
  register int i;


  if (!cumul_initial) {
    cumul_initial = new double[type == 'o' ? nb_state : nb_row];
  }

  if (!cumul_transition) {
    cumul_transition = new double*[nb_row];
    for (i = 0;i < nb_row;i++) {
      cumul_transition[i] = new double[nb_state];
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Destruction des champs de type "cumul" d'un objet Chain.
 *
 *--------------------------------------------------------------*/

void Chain::remove_cumul()

{
  register int i;


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


/*--------------------------------------------------------------*
 *
 *  Calcul des fonctions de repartition.
 *
 *--------------------------------------------------------------*/

void Chain::cumul_computation()

{
  if ((cumul_initial) && (cumul_transition)) {
    register int i;


    stat_tool::cumul_computation((type == 'o' ? nb_state : nb_row) , initial , cumul_initial);

    for (i = 0;i < nb_row;i++) {
      stat_tool::cumul_computation(nb_state , transition[i] , cumul_transition[i]);
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Transfomation en log des parametres d'une chaine de Markov.
 *
 *--------------------------------------------------------------*/

void Chain::log_computation()

{
  if ((cumul_initial) && (cumul_transition)) {
    register int i;


    stat_tool::log_computation((type == 'o' ? nb_state : nb_row) , initial , cumul_initial);

    for (i = 0;i < nb_row;i++) {
      stat_tool::log_computation(nb_state , transition[i] , cumul_transition[i]);
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Constructeur par copie de la classe ChainData.
 *
 *  argument : reference sur un objet ChainData.
 *
 *--------------------------------------------------------------*/

ChainData::ChainData(const ChainData &chain_data)

{
  copy(chain_data);
}


};  // namespace stat_tool
