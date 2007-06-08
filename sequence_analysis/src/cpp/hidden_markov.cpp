/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       AMAPmod: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2002 UMR Cirad/Inra Modelisation des Plantes
 *
 *       File author(s): Y. Guedon (yann.guedon@cirad.fr)
 *
 *       $Source$
 *       $Id$
 *
 *       Forum for AMAPmod developers: amldevlp@cirad.fr
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



#include "tool/rw_tokenizer.h"
#include "tool/rw_cstring.h"
// #include <rw/vstream.h>
// #include <rw/rwfile.h>
#include "stat_tool/stat_tools.h"
#include "stat_tool/regression.h"
#include "stat_tool/curves.h"
#include "stat_tool/markovian.h"
#include "stat_tool/stat_label.h"
#include "sequences.h"
#include "markov.h"
#include "hidden_markov.h"
#include "sequence_label.h"

using namespace std;



/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Markov.
 *
 *  arguments : pointeurs sur un objet Chain, ordre, nombre de processus d'observation,
 *              pointeurs sur des objets Nonparametric_process, longueur des sequences.
 *
 *--------------------------------------------------------------*/

Markov::Markov(const Chain *pchain , int iorder , int inb_output_process ,
               Nonparametric_process **pobservation , int length)
:Chain(*pchain)

{
  register int i;


  nb_iterator = 0;
  markov_data = 0;

  order = iorder;

  self_row = new int[nb_state];
  self_row_computation();

  homogeneous = true;

  homogeneity = new bool[nb_state];
  self_transition = new Function*[nb_state];

  for (i = 0;i < nb_state;i++) {
    homogeneity[i] = true;
    self_transition[i] = 0;
  }

  nb_output_process = inb_output_process;
  process = new Nonparametric_sequence_process*[nb_output_process + 1];
  process[0] = new Nonparametric_sequence_process(nb_state , nb_state , false);
  for (i = 1;i <= nb_output_process;i++) {
    process[i] = new Nonparametric_sequence_process(*pobservation[i - 1]);
  }

  component_computation();
  characteristic_computation(length , true);
}


/*--------------------------------------------------------------*
 *
 *  Construction d'une chaine de Markov a partir d'une chaine de Markov
 *  initiale en ajoutant un etat.
 *
 *  arguments : reference sur un objet Chain, ordre, etat de reference.
 *
 *--------------------------------------------------------------*/

Chain::Chain(const Chain &chain , int order , int state)

{
  register int i , j , k;
  int row , power[ORDER] , state_index[ORDER];
  double norm;


  type = chain.type;
  nb_state = chain.nb_state + 1;

  nb_row = 1;
  for (i = 0;i < order;i++) {
    nb_row *= nb_state;
  }

  if (chain.nb_component > 0) {
    accessibility = new bool*[nb_state];
    for (i = 0;i <= state;i++) {
      accessibility[i] = new bool[nb_state];
      for (j = 0;j <= state;j++) {
        accessibility[i][j] = chain.accessibility[i][j];
      }
      for (j = state + 1;j < nb_state;j++) {
        accessibility[i][j] = chain.accessibility[i][j - 1];
      }
    }

    for (i = state + 1;i < nb_state;i++) {
      accessibility[i] = new bool[nb_state];
      for (j = 0;j <= state;j++) {
        accessibility[i][j] = chain.accessibility[i - 1][j];
      }
      for (j = state + 1;j < nb_state;j++) {
        accessibility[i][j] = chain.accessibility[i - 1][j - 1];
      }
    }

    nb_component = chain.nb_component;
    component_nb_state = new int[nb_component];
    component = new int*[nb_component];

    for (i = 0;i < nb_component;i++) {
      for (j = 0;j < chain.component_nb_state[i];j++) {
        if (chain.component[i][j] == state) {
          break;
        }
      }

      if (j < chain.component_nb_state[i]) {
        component_nb_state[i] = chain.component_nb_state[i] + 1;
        component[i] = new int[component_nb_state[i]];
        for (k = 0;k <= j;k++) {
          if (chain.component[i][k] <= state) {
            component[i][k] = chain.component[i][k];
          }
          else {
            component[i][k] = chain.component[i][k] + 1;
          }
        }
        component[i][j + 1] = state + 1;
        for (k = j + 2;k < component_nb_state[i];k++) {
          if (chain.component[i][k - 1] <= state) {
            component[i][k] = chain.component[i][k - 1];
          }
          else {
            component[i][k] = chain.component[i][k - 1] + 1;
          }
        }
      }

      else {
        component_nb_state[i] = chain.component_nb_state[i];
        component[i] = new int[component_nb_state[i]];
        for (k = 0;k < component_nb_state[i];k++) {
          if (chain.component[i][k] <= state) {
            component[i][k] = chain.component[i][k];
          }
          else {
            component[i][k] = chain.component[i][k] + 1;
          }
        }
      }
    }

    state_type = new char[nb_state];
    for (i = 0;i <= state;i++) {
      state_type[i] = chain.state_type[i];
    }
    for (i = state + 1;i < nb_state;i++) {
      state_type[i] = chain.state_type[i - 1];
    }
  }

  else {
    accessibility = 0;
    nb_component = 0;
    component_nb_state = 0;
    component = 0;
    state_type = 0;
  }

  initial = new double[nb_state];
  norm = 1. + chain.initial[state];
  for (i = 0;i <= state;i++) {
    initial[i] = chain.initial[i] / norm;
  }
  for (i = state + 1;i < nb_state;i++) {
    initial[i] = chain.initial[i - 1] / norm;
  }

  cumul_initial = 0;

  i = 1;
  for (j = 0;j < order;j++) {
    power[j] = i;
    i *= chain.nb_state;
  }

  for (i = 0;i < order;i++) {
    state_index[i] = 0;
  }

  transition = new double*[nb_row];
  for (i = 0;i < nb_row;i++) {
    transition[i] = new double[nb_state];

    row = 0;
    for (j = 0;j < order;j++) {
      if (state_index[j] <= state) {
        row += state_index[j] * power[j];
      }
      else {
        row += (state_index[j] - 1) * power[j];
      }
    }

    norm = 1. + chain.transition[row][state];

    for (j = 0;j <= state;j++) {
      transition[i][j] = chain.transition[row][j] / norm;
    }
    for (j = state + 1;j < nb_state;j++) {
      transition[i][j] = chain.transition[row][j - 1] / norm;
    }

    // mise a jour des indices des etats

    for (j = 0;j < order;j++) {
      if (state_index[j] < nb_state - 1) {
        state_index[j]++;
        break;
      }
      else {
        state_index[j] = 0;
      }
    }
  }

  cumul_transition = 0;
}


/*--------------------------------------------------------------*
 *
 *  Construction d'un objet Markov a partir d'un objet Markov initial
 *  en ajoutant un etat.
 *
 *  arguments : reference sur un objet Markov,
 *              etat de reference.
 *
 *--------------------------------------------------------------*/

Markov::Markov(const Markov &markov , int state)
:Chain(markov , markov.order , state)

{
  register int i;


  nb_iterator = 0;
  markov_data = 0;

  order = markov.order;
  
  self_row = new int[nb_state];
  self_row_computation();

  homogeneous = true;

  homogeneity = new bool[nb_state];
  self_transition = new Function*[nb_state];

  for (i = 0;i < nb_state;i++) {
    homogeneity[i] = true;
    self_transition[i] = 0;
  }

  nb_output_process = markov.nb_output_process;

  process = new Nonparametric_sequence_process*[nb_output_process + 1];
  process[0] = new Nonparametric_sequence_process(nb_state , nb_state , false);
  for (i = 1;i <= nb_output_process;i++) {
    process[i] = new Nonparametric_sequence_process(*(markov.process[i]) , 's' , state);
  }
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Hidden_markov.
 *
 *  arguments : reference sur un objet Hidden_markov,
 *              probabilite de rester dans un etat.
 *
 *--------------------------------------------------------------*/

Hidden_markov::Hidden_markov(const Hidden_markov &hmarkov , double self_transition)
:Markov(hmarkov , false)

{
  register int i;


  // initialisation des parametres de la chaine de Markov

  init(self_transition);

  // initialisation des lois d'observation

  for (i = 0;i < nb_output_process;i++) {
    process[i + 1]->init();
  }
}


/*--------------------------------------------------------------*
 *
 *  Destructeur de la classe Hidden_markov.
 *
 *--------------------------------------------------------------*/

Hidden_markov::~Hidden_markov() {}


/*--------------------------------------------------------------*
 *
 *  Application d'un seuil sur les parametres d'une chaine de Markov cachee.
 *
 *  argument : probabilite minimum.
 *
 *--------------------------------------------------------------*/

Hidden_markov* Hidden_markov::thresholding(double min_probability) const

{
  register int i;
  Hidden_markov *hmarkov;


  hmarkov = new Hidden_markov(*this , false , false);

  hmarkov->Chain::thresholding(min_probability);
  hmarkov->component_computation();

  for (i = 1;i <= nb_output_process;i++) {
    process[i]->thresholding(min_probability);
  }

  return hmarkov;
}


/*--------------------------------------------------------------*
 *
 *  Construction d'un objet Hidden_markov a partir d'un fichier.
 *
 *  arguments : reference sur un objet Format_error, path,
 *              longueur des sequences, flag format des processus d'observation.
 *
 *--------------------------------------------------------------*/

Hidden_markov* hidden_markov_ascii_read(Format_error &error , const char *path ,
                                        int length , bool old_format)

{
  RWCString buffer , token;
  size_t position;
  bool status;
  register int i;
  int line , order , nb_output_process;
  const Chain *chain;
  Nonparametric_process **observation;
  Hidden_markov *hmarkov;
  ifstream in_file(path);


  hmarkov = 0;
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

        // test mot cle HIDDEN_MARKOV_CHAIN

        if (i == 0) {
          if (token != SEQ_word[SEQW_HIDDEN_MARKOV_CHAIN]) {
            status = false;
            error.correction_update(STAT_parsing[STATP_KEY_WORD] , SEQ_word[SEQW_HIDDEN_MARKOV_CHAIN] , line);
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

    // analyse du format et lecture de la chaine de Markov

    chain = chain_parsing(error , in_file , line , 'o' , true , order);

    if (chain) {

      // analyse du format et lecture des lois d'observation

      switch (old_format) {
      case false :
        observation = observation_parsing(error , in_file , line , chain->nb_state ,
                                          nb_output_process);
        break;
      case true :
        observation = old_observation_parsing(error , in_file , line , chain->nb_state ,
                                              nb_output_process);
        break;
      }

      if (observation) {
        if (status) {
          hmarkov = new Hidden_markov(chain , order , nb_output_process , observation , length);
        }

        for (i = 0;i < nb_output_process;i++) {
          delete observation[i];
        }
        delete [] observation;
      }

      delete chain;
    }
  }

  return hmarkov;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Hidden_markov dans un fichier.
 *
 *  arguments : stream, flag niveau de detail.
 *
 *--------------------------------------------------------------*/

ostream& Hidden_markov::ascii_write(ostream &os , bool exhaustive) const

{
  Markov::ascii_write(os , markov_data , exhaustive , false , true);

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Hidden_markov dans un fichier.
 *
 *  arguments : reference sur un objet Format_error, path,
 *              flag niveau de detail.
 *
 *--------------------------------------------------------------*/

bool Hidden_markov::ascii_write(Format_error &error , const char *path ,
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
    Markov::ascii_write(out_file , markov_data , exhaustive , true , true);
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Hidden_markov dans un fichier au format tableur.
 *
 *  arguments : reference sur un objet Format_error, path.
 *
 *--------------------------------------------------------------*/

bool Hidden_markov::spreadsheet_write(Format_error &error , const char *path) const

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
    Markov::spreadsheet_write(out_file , markov_data , true);
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Fonctions pour la persistance.
 *
 *--------------------------------------------------------------*/

/* RWDEFINE_COLLECTABLE(Hidden_markov , STATI_HIDDEN_MARKOV);


RWspace Hidden_markov::binaryStoreSize() const

{
  return Markov::binaryStoreSize();
}


void Hidden_markov::restoreGuts(RWvistream &is)

{
  Markov::restoreGuts(is);
}


void Hidden_markov::restoreGuts(RWFile &file)

{
  Markov::restoreGuts(file);
}


void Hidden_markov::saveGuts(RWvostream &os) const

{
  Markov::saveGuts(os);
}


void Hidden_markov::saveGuts(RWFile &file) const

{
  Markov::saveGuts(file);
} */
