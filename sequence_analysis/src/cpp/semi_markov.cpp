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



#include <sstream>
#include "tool/rw_tokenizer.h"
#include "tool/rw_cstring.h"
#include "tool/rw_locale.h"
// #include <rw/vstream.h>
// #include <rw/rwfile.h>
#include "stat_tool/stat_tools.h"
#include "stat_tool/distribution.h"
#include "stat_tool/curves.h"
#include "stat_tool/markovian.h"
#include "stat_tool/stat_label.h"
#include "sequences.h"
#include "semi_markov.h"
#include "sequence_label.h"

using namespace std;



/*--------------------------------------------------------------*
 *
 *  Constructeur par defaut de la classe Semi_markov.
 *
 *--------------------------------------------------------------*/

Semi_markov::Semi_markov()

{
  nb_iterator = 0;
  semi_markov_data = 0;

  state_subtype = 0;
  forward = 0;

  nb_output_process = 0;
  nonparametric_process = 0;
  parametric_process = 0;
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Semi_markov.
 *
 *  arguments : type, nombre d'etats, nombre de processus d'observation,
 *              nombre de valeurs observees par processus.
 *
 *--------------------------------------------------------------*/

Semi_markov::Semi_markov(char itype , int inb_state , int inb_output_process , int *nb_value)
:Chain(itype , inb_state)

{
  register int i;


  nb_iterator = 0;
  semi_markov_data = 0;

  state_subtype = 0;
  forward = 0;

  nb_output_process = inb_output_process;

  nonparametric_process = new Nonparametric_sequence_process*[nb_output_process + 1];
  nonparametric_process[0] = new Nonparametric_sequence_process(nb_state , nb_state , false);
  parametric_process = new Parametric_process*[nb_output_process + 1];
  parametric_process[0] = 0;

  for (i = 1;i <= nb_output_process;i++) {
    if (*nb_value <= NB_OUTPUT) {
      nonparametric_process[i] = new Nonparametric_sequence_process(nb_state , *nb_value++ , true);
      parametric_process[i] = 0;
    }
    else {
      nonparametric_process[i] = 0;
      parametric_process[i] = new Parametric_process(nb_state , (int)(*nb_value++ * SAMPLE_NB_VALUE_COEFF));
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Semi_markov.
 *
 *  arguments : pointeur sur un objet Chain, sur un objet Nonparametric_sequence_process et
 *              sur un objet Nonparametric_process, longueur des sequences,
 *              flag sur le calcul des lois de comptage.
 *
 *--------------------------------------------------------------*/

Semi_markov::Semi_markov(const Chain *pchain , const Nonparametric_sequence_process *poccupancy ,
                         const Nonparametric_process *pobservation ,
                         int length , bool counting_flag)
:Chain(*pchain)

{
  register int i;


  nb_iterator = 0;
  semi_markov_data = 0;

  nb_output_process = (pobservation ? 1 : 0);
  nonparametric_process = new Nonparametric_sequence_process*[nb_output_process + 1];
  nonparametric_process[0] = new Nonparametric_sequence_process(*poccupancy);
  if (pobservation) {
    nonparametric_process[1] = new Nonparametric_sequence_process(*pobservation);
  }

  parametric_process = 0;

  for (i = 0;i < nb_state;i++) {
    if (transition[i][i] < 1.) {
      nonparametric_process[0]->absorption[i] = 0.;
    }
    else {
      nonparametric_process[0]->absorption[i] = 1.;
    }
  }

  state_subtype = new int[nb_state];
  forward = new Forward*[nb_state];

  for (i = 0;i < nb_state;i++) {
    state_subtype[i] = (nonparametric_process[0]->sojourn_time[i] ? SEMI_MARKOVIAN : MARKOVIAN);

    if ((state_subtype[i] == SEMI_MARKOVIAN) && (state_type[i] == 'r')) {
      forward[i] = new Forward(*(nonparametric_process[0]->sojourn_time[i]));
    }
    else {
      forward[i] = 0;
    }
  }

  if (type == 'e') {
    for (i = 0;i < nb_state;i++) {
      initial[i] = 1. / (double)nb_state;
    }
    initial_probability_computation();
  }

  if (length > COUNTING_MAX_LENGTH) {
    counting_flag = false;
  }
  characteristic_computation(length , counting_flag);
}


/*--------------------------------------------------------------*
 *
 *  Copie d'un objet Semi_markov.
 *
 *  arguments : reference sur un objet Semi_markov,
 *              flag copie de l'objet Semi_markov_data,
 *              parametre (si strictement positif :
 *              nombre de valeurs allouees pour les lois d'occupation des etats).
 *
 *--------------------------------------------------------------*/

void Semi_markov::copy(const Semi_markov &smarkov , bool data_flag , int param)

{
  register int i;


  nb_iterator = 0;

  if ((data_flag) && (smarkov.semi_markov_data)) {
    semi_markov_data = new Semi_markov_data(*(smarkov.semi_markov_data) , false);
  }
  else {
    semi_markov_data = 0;
  }

  state_subtype = new int[nb_state];
  for (i = 0;i < nb_state;i++) {
    state_subtype[i] = smarkov.state_subtype[i];
  }

  forward = new Forward*[nb_state];

  for (i = 0;i < nb_state;i++) {
    if (smarkov.forward[i]) {
      forward[i] = new Forward(*(smarkov.forward[i]) , param);
    }
    else {
      forward[i] = 0;
    }
  }

  nb_output_process = smarkov.nb_output_process;

  nonparametric_process = new Nonparametric_sequence_process*[nb_output_process + 1];

  if (smarkov.parametric_process) {
    parametric_process = new Parametric_process*[nb_output_process + 1];
    parametric_process[0] = 0;
  }
  else {
    parametric_process = 0;
  }

  switch (param) {

  case I_DEFAULT : {
    for (i = 0;i <= nb_output_process;i++) {
      if (smarkov.nonparametric_process[i]) {
        nonparametric_process[i] = new Nonparametric_sequence_process(*(smarkov.nonparametric_process[i]));
      }
    }
    break;
  }

  case 0 : {
    nonparametric_process[0] = new Nonparametric_sequence_process(*(smarkov.nonparametric_process[0]) ,
                                                                  'o' , I_DEFAULT);
    for (i = 1;i <= nb_output_process;i++) {
      if (smarkov.nonparametric_process[i]) {
        nonparametric_process[i] = new Nonparametric_sequence_process(*(smarkov.nonparametric_process[i]) ,
                                                                      'c' , false);
      }
    }
    break;
  }

  default : {
    nonparametric_process[0] = new Nonparametric_sequence_process(*(smarkov.nonparametric_process[0]) ,
                                                                  'o' , param);
    for (i = 1;i <= nb_output_process;i++) {
      if (smarkov.nonparametric_process[i]) {
        nonparametric_process[i] = new Nonparametric_sequence_process(*(smarkov.nonparametric_process[i]) ,
                                                                      'c' , false);
      }
    }
    break;
  }
  }

  for (i = 1;i <= nb_output_process;i++) {
    if (smarkov.nonparametric_process[i]) {
      if (smarkov.parametric_process) {
        parametric_process[i] = 0;
      }
    }
    else {
      nonparametric_process[i] = 0;
      parametric_process[i] = new Parametric_process(*(smarkov.parametric_process[i]));
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Destruction des champs d'un objet Semi_markov.
 *
 *--------------------------------------------------------------*/

void Semi_markov::remove()

{
  register int i;


  delete semi_markov_data;

  delete [] state_subtype;

  if (forward) {
    for (i = 0;i < nb_state;i++) {
      delete forward[i];
    }
    delete [] forward;
  }

  if (nonparametric_process) {
    for (i = 0;i <= nb_output_process;i++) {
      delete nonparametric_process[i];
    }
    delete [] nonparametric_process;
  }

  if (parametric_process) {
    for (i = 1;i <= nb_output_process;i++) {
      delete parametric_process[i];
    }
    delete [] parametric_process;
  }
}


/*--------------------------------------------------------------*
 *
 *  Destructeur de la classe Semi_markov.
 *
 *--------------------------------------------------------------*/

Semi_markov::~Semi_markov()

{
  remove();
}


/*--------------------------------------------------------------*
 *
 *  Destruction d'un objet Semi_markov en tenant compte du nombre
 *  d'iterateurs pointant dessus.
 *
 *--------------------------------------------------------------*/

void Semi_markov::conditional_delete()

{
  if (nb_iterator == 0) {
    delete this;
  }
}


/*--------------------------------------------------------------*
 *
 *  Operateur d'assignement de la classe Semi_markov.
 *
 *  argument : reference sur un objet Semi_markov.
 *
 *--------------------------------------------------------------*/

Semi_markov& Semi_markov::operator=(const Semi_markov &smarkov)

{
  if ((&smarkov != this) && (nb_iterator == 0)) {
    remove();
    Chain::remove();

    Chain::copy(smarkov);
    copy(smarkov);
  }

  return *this;
}


/*--------------------------------------------------------------*
 *
 *  Extraction d'une loi.
 *
 *  arguments : reference sur un objet Format_error, type de loi,
 *              variable, valeur (etat dans le cas des lois d'observation).
 *
 *--------------------------------------------------------------*/

Parametric_model* Semi_markov::extract(Format_error &error , int type ,
                                       int variable , int value) const

{
  bool status = true;
  int hvariable;
  Distribution *pdist;
  Parametric *pparam;
  Parametric_model *dist;
  Histogram *phisto;


  dist = 0;
  error.init();

  pdist = 0;
  pparam = 0;

  if (type == OBSERVATION) {
    if ((variable < 1) || (variable > nb_output_process)) {
      status = false;
      error.update(SEQ_error[SEQR_OUTPUT_PROCESS_INDEX]);
    }

    else {
      if ((value < 0) || (value >= nb_state)) {
        status = false;
        ostringstream error_message;
        error_message << STAT_label[STATL_STATE] << " " << value << " "
                      << SEQ_error[SEQR_NOT_PRESENT];
        error.update((error_message.str()).c_str());
      }

      else {
        if (nonparametric_process[variable]) {
          pdist = nonparametric_process[variable]->observation[value];
        }
        else {
          pparam = parametric_process[variable]->observation[value];
        }
      }
    }
  }

  else {
    if ((variable < 0) || (variable > nb_output_process)) {
      status = false;
      error.update(SEQ_error[SEQR_OUTPUT_PROCESS_INDEX]);
    }

    else {
      if (!nonparametric_process[variable]) {
        status = false;
        ostringstream error_message;
        error_message << STAT_label[STATL_VARIABLE] << " " << variable << ": " 
                      << SEQ_error[SEQR_CHARACTERISTICS_NOT_COMPUTED];
        error.update((error_message.str()).c_str());
      }

      else if ((value < 0) || (value >= nonparametric_process[variable]->nb_value)) {
        status = false;
        ostringstream error_message;
        error_message << STAT_label[variable == 0 ? STATL_STATE : STATL_OUTPUT] << " "
                      << value << " " << SEQ_error[SEQR_NOT_PRESENT];
        error.update((error_message.str()).c_str());
      }

      if (status) {
        switch (type) {
        case FIRST_OCCURRENCE :
          pdist = nonparametric_process[variable]->first_occurrence[value];
          break;
        case RECURRENCE_TIME :
          pdist = nonparametric_process[variable]->recurrence_time[value];
          break;
        case SOJOURN_TIME :
          pparam = nonparametric_process[variable]->sojourn_time[value];
          break;
        case NB_RUN :
          pdist = nonparametric_process[variable]->nb_run[value];
          break;
        case NB_OCCURRENCE :
          pdist = nonparametric_process[variable]->nb_occurrence[value];
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
    phisto = 0;

    if (semi_markov_data) {
      switch (semi_markov_data->type[0]) {
      case INT_VALUE :
        hvariable = variable - 1;
        break;
      case STATE :
        hvariable = variable;
        break;
      }

      if (hvariable >= 0) {
        switch (type) {

        case OBSERVATION : {
          if ((semi_markov_data->observation) && (semi_markov_data->observation[hvariable])) {
            phisto = semi_markov_data->observation[hvariable][value];
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
      dist = new Parametric_model(*pdist , phisto);
    }
    else if (pparam) {
      dist = new Parametric_model(*pparam , phisto);
    }
  }

  return dist;
}


/*--------------------------------------------------------------*
 *
 *  Extraction d'une loi de l'intervalle de temps residuel.
 *
 *  arguments : reference sur un objet Format_error, etat,
 *              type d'histogramme associe.
 *
 *--------------------------------------------------------------*/

Parametric_model* Semi_markov::extract(Format_error &error , int state , int histogram_type) const

{
  bool status = true;
  Distribution *pdist;
  Parametric_model *dist;
  Histogram *phisto;


  dist = 0;
  error.init();

  pdist = 0;

  if ((state < 0) || (state >= nb_state)) {
    status = false;
    ostringstream error_message;
    error_message << STAT_label[STATL_STATE] << " " << state << " "
                  << SEQ_error[SEQR_NOT_PRESENT];
    error.update((error_message.str()).c_str());
  }

  else {
    pdist = forward[state];

    if (!pdist) {
      status = false;
      error.update(SEQ_error[SEQR_NON_EXISTING_FORWARD_DISTRIBUTION]);
    }

    else {
      phisto = 0;

      if ((semi_markov_data) && (semi_markov_data->type[0] == STATE)) {
        switch (histogram_type) {

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

      dist = new Parametric_model(*pdist , phisto);
    }
  }

  return dist;
}


/*--------------------------------------------------------------*
 *
 *  Extraction de la partie "donnees" d'un objet Semi_markov.
 *
 *  argument : reference sur un objet Format_error.
 *
 *--------------------------------------------------------------*/

Semi_markov_data* Semi_markov::extract_data(Format_error &error) const

{
  bool status = true;
  Semi_markov_data *seq;


  seq = 0;
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
    seq = new Semi_markov_data(*semi_markov_data);
    seq->semi_markov = new Semi_markov(*this , false);
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Application d'un seuil sur les parametres d'une semi-chaine de Markov.
 *
 *  argument : probabilite minimum.
 *
 *--------------------------------------------------------------*/

Semi_markov* Semi_markov::thresholding(double min_probability) const

{
  register int i;
  Semi_markov *smarkov;


  smarkov = new Semi_markov(*this , false , false);

  smarkov->Chain::thresholding(min_probability);
  smarkov->component_computation();

  for (i = 1;i <= smarkov->nb_output_process;i++) {
    if (smarkov->nonparametric_process[i]) {
      smarkov->nonparametric_process[i]->thresholding(min_probability);
    }
  }

  return smarkov;
}


/*--------------------------------------------------------------*
 *
 *  Construction d'un objet Semi_markov a partir d'un fichier.
 *
 *  arguments : reference sur un objet Format_error, path,
 *              longueur des sequences, flag sur le calcul des lois de comptage,
 *              seuil sur les fonctions de repartition des lois d'occupation des etats.
 *
 *--------------------------------------------------------------*/

Semi_markov* semi_markov_ascii_read(Format_error &error , const char *path , int length ,
                                    bool counting_flag , double cumul_threshold)

{
  RWCString buffer , token;
  size_t position;
  char type = 'v';
  bool status;
  register int i;
  int line , order = 1;
  const Chain *chain;
  const Nonparametric_sequence_process *occupancy;
  const Nonparametric_process *observation;
  Semi_markov *smarkov;
  ifstream in_file(path);


  smarkov = 0;
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

        // test mot cle (EQUILIBRIUM) SEMI-MARKOV_CHAIN

        if (i == 0) {
          if (token == SEQ_word[SEQW_SEMI_MARKOV_CHAIN]) {
            type = 'o';
          }
          else if (token == SEQ_word[SEQW_EQUILIBRIUM_SEMI_MARKOV_CHAIN]) {
            type = 'e';
          }
          else {
            status = false;
            ostringstream correction_message;
            correction_message << SEQ_word[SEQW_SEMI_MARKOV_CHAIN] << " or "
                               << SEQ_word[SEQW_EQUILIBRIUM_SEMI_MARKOV_CHAIN];
            error.correction_update(STAT_parsing[STATP_KEY_WORD] ,
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

    if (type != 'v') {

      // analyse du format et lecture de la chaine de Markov

      chain = chain_parsing(error , in_file , line , type , false , order);

      if (chain) {

        // analyse du format et lecture des lois d'occupation de etats

        occupancy = occupancy_parsing(error , in_file , line , *chain , cumul_threshold);
        if (!occupancy) {
          status = false;
        }

        // analyse du format et lecture des lois d'observation

        observation = 0;

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

            // test mot cle OUTPUT_PROCESS

            if (i == 0) {
              if (token != STAT_word[STATW_OUTPUT_PROCESS]) {
                status = false;
                error.correction_update(STAT_parsing[STATP_KEY_WORD] , STAT_word[STATW_OUTPUT_PROCESS] , line);
              }
            }

            i++;
          }

          if (i > 0) {
            if (i != 1) {
              status = false;
              error.update(STAT_parsing[STATP_FORMAT] , line);
            }

            // analyse du format et lecture des lois d'observation

            observation = observation_parsing(error , in_file , line , chain->nb_state , false);
            if (!observation) {
              status = false;
            }

            break;
          }
        }

        while (buffer.readLine(in_file , false)) {
          line++;

#         ifdef DEBUG
          cout << line << "  " << buffer << endl;
#         endif

          position = buffer.first('#');
          if (position != RW_NPOS) {
            buffer.remove(position);
          }
          if (!(buffer.isNull())) {
            status = false;
            error.update(STAT_parsing[STATP_FORMAT] , line);
          }
        }

        if (status) {
          smarkov = new Semi_markov(chain , occupancy , observation , length , counting_flag);
        }

        delete chain;
        delete occupancy;
        delete observation;
      }
    }
  }

  return smarkov;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture sur une ligne d'un objet Semi_markov.
 *
 *  argument : stream.
 *
 *--------------------------------------------------------------*/

ostream& Semi_markov::line_write(ostream &os) const

{
  os << nb_state << " " << STAT_word[STATW_STATES];

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Semi_markov et de la structure de donnees associee.
 *
 *  arguments : stream, pointeur sur un objet Semi_markov_data,
 *              flag niveau de detail, flag fichier, flag semi-Markov cache.
 *
 *--------------------------------------------------------------*/

ostream& Semi_markov::ascii_write(ostream &os , const Semi_markov_data *seq ,
                                  bool exhaustive , bool file_flag , bool hidden) const

{
  register int i;
  int variable;
  Histogram **observation = 0;
  Sequence_characteristics *characteristics = 0;


  switch (hidden) {

  case false : {
    switch (type) {
    case 'o' :
      os << SEQ_word[SEQW_SEMI_MARKOV_CHAIN] << endl;
      break;
    case 'e' :
      os << SEQ_word[SEQW_EQUILIBRIUM_SEMI_MARKOV_CHAIN] << endl;
      break;
    }
    break;
  }

  case true : {
    switch (type) {
    case 'o' :
      os << SEQ_word[SEQW_HIDDEN_SEMI_MARKOV_CHAIN] << endl;
      break;
    case 'e' :
      os << SEQ_word[SEQW_EQUILIBRIUM_HIDDEN_SEMI_MARKOV_CHAIN] << endl;
      break;
    }
    break;
  }
  }

  // ecriture des parametres de la chaine de Markov

  ascii_print(os , file_flag);

  // ecriture des lois d'occupation des etats

  if ((seq) && (seq->type[0] == STATE)) {
    characteristics = seq->characteristics[0];
  }
  nonparametric_process[0]->ascii_print(os , 0 , 0 , characteristics , exhaustive ,
                                        file_flag , forward);

  if (hidden) {
    if ((type == 'o') && (nonparametric_process[0]->index_value)) {
      register int j;
      double sum , *state_marginal;


      state_marginal = new double[nb_state];

      for (i = 0;i < nb_state;i++) {
        state_marginal[i] = 0.;
      }

      for (i = 0;i < nonparametric_process[0]->length->nb_value - 1;i++) {
        for (j = 0;j < nb_state;j++) {
          state_marginal[j] += nonparametric_process[0]->index_value->point[j][i] *
                               (1. - nonparametric_process[0]->length->cumul[i]);
        }
      }

      sum = 0.;
      for (i = 0;i < nb_state;i++) {
        sum += state_marginal[i];
      }

      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << SEQ_label[SEQL_STATE_PROBABILITY] << endl;
      if (file_flag) {
        os << "# ";
      }
      for (i = 0;i < nb_state;i++) {
        os << state_marginal[i] / sum  << "  ";
      }
      os << endl;

      delete state_marginal;
    }

    os << "\n" << nb_output_process << " "
       << STAT_word[nb_output_process == 1 ? STATW_OUTPUT_PROCESS : STATW_OUTPUT_PROCESSES] << endl;
  }

  // ecriture des lois associees a chaque processus d'observation

  for (i = 1;i <= nb_output_process;i++) {
    os << "\n" << STAT_word[STATW_OUTPUT_PROCESS];

    if (hidden) {
      os << " " << i;

      if (nonparametric_process[i]) {
        os << " : " << STAT_word[STATW_NONPARAMETRIC];
      }
      else {
        os << " : " << STAT_word[STATW_PARAMETRIC];
      }
    }
    os << endl;

    if (seq) {
      switch (seq->type[0]) {
      case INT_VALUE :
        variable = i - 1;
        break;
      case STATE :
        variable = i;
        break;
      }

      if (seq->observation) {
        observation = seq->observation[variable];
      }
      if (seq->characteristics[variable]) {
        characteristics = seq->characteristics[variable];
      }
    }

    if (nonparametric_process[i]) {
      nonparametric_process[i]->ascii_print(os , i , observation , characteristics ,
                                            exhaustive , file_flag);
    }
    else {
      parametric_process[i]->ascii_print(os , observation , exhaustive , file_flag);
    }
  }

  if (seq) {
    int nb_parameter = nb_parameter_computation(hidden ? MIN_PROBABILITY : 0.);
    double information , likelihood;


    // ecriture de l'histogramme des longueurs des sequences

    os << "\n";
    if (file_flag) {
      os << "# ";
    }
    os << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_HISTOGRAM] << " - ";
    seq->hlength->ascii_characteristic_print(os , false , file_flag);

    if (exhaustive) {
      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << "   | " << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_HISTOGRAM] << endl;
      seq->hlength->ascii_print(os , file_flag);
    }

    os << "\n";
    if (file_flag) {
      os << "# ";
    }
    os << SEQ_label[SEQL_CUMUL_LENGTH] << ": " << seq->cumul_length << endl;

    // ecriture de la quantite d'information des sequences dans le cas iid

    information = seq->iid_information_computation();

    os << "\n";
    if (file_flag) {
      os << "# ";
    }
    os << SEQ_label[SEQL_IID_INFORMATION] << ": " << information << " ("
       << information / seq->cumul_length << ")" << endl;

    // ecriture des vraisemblances des sequences

    switch (hidden) {

    case false : {
      likelihood = seq->likelihood;

      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << STAT_label[STATL_LIKELIHOOD] << ": " << likelihood << "   ("
         << STAT_label[STATL_NORMALIZED] << ": " << likelihood / seq->cumul_length << ")" << endl;
      break;
    }

    case true : {
      if (seq->likelihood != D_INF) {
        os << "\n";
        if (file_flag) {
          os << "# ";
        }
        os << SEQ_label[SEQL_STATE_SEQUENCES_LIKELIHOOD] << ": " << seq->likelihood << "   ("
           << STAT_label[STATL_NORMALIZED] << ": " << seq->likelihood / seq->cumul_length << ")" << endl;
      }

      likelihood = seq->hidden_likelihood;

      if (likelihood != D_INF) {
        os << "\n";
        if (file_flag) {
          os << "# ";
        }
        os << SEQ_label[SEQL_OBSERVED_SEQUENCES_LIKELIHOOD] << ": " << likelihood << "   ("
           << STAT_label[STATL_NORMALIZED] << ": " << likelihood / seq->cumul_length << ")" << endl;
      }
      break;
    }
    }

    if (likelihood != D_INF) {
      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << nb_parameter << " " << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS]
         << "   2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[AIC] << "): "
         << 2 * (likelihood - nb_parameter) << endl;

      if (nb_parameter < seq->cumul_length - 1) {
        os << "\n";
        if (file_flag) {
          os << "# ";
        }
        os << nb_parameter << " " << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS]
           << "   2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[AICc] << "): "
           << 2 * (likelihood - (double)(nb_parameter * seq->cumul_length) /
             (double)(seq->cumul_length - nb_parameter - 1)) << endl;
      }

      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << nb_parameter << " " << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS]
         << "   2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[BIC] << "): "
         << 2 * likelihood - nb_parameter * log((double)seq->cumul_length) << endl;

      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << nb_parameter << " " << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS]
         << "   2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[BICc] << "): "
         << 2 * likelihood - penalty_computation(hidden , (hidden ? MIN_PROBABILITY : 0.)) << endl;
    }

    if ((hidden) && (seq->likelihood != D_INF)) {
      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << nb_parameter << " " << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS]
         << "   2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[ICL] << "): "
         << 2 * seq->likelihood - nb_parameter * log((double)seq->cumul_length) << endl;

      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << nb_parameter << " " << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS]
         << "   2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[ICLc] << "): "
         << 2 * seq->likelihood - penalty_computation(hidden , MIN_PROBABILITY) << endl;
    }
  }

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Semi_markov.
 *
 *  arguments : stream, flag niveau de detail.
 *
 *--------------------------------------------------------------*/

ostream& Semi_markov::ascii_write(ostream &os , bool exhaustive) const

{
  return ascii_write(os , semi_markov_data , exhaustive , false);
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Semi_markov dans un fichier.
 *
 *  arguments : reference sur un objet Format_error, path,
 *              flag niveau de detail.
 *
 *--------------------------------------------------------------*/

bool Semi_markov::ascii_write(Format_error &error , const char *path ,
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
    ascii_write(out_file , semi_markov_data , exhaustive , true);
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Semi_markov et de la structure de donnees associee
 *  dans un fichier au format tableur.
 *
 *  arguments : stream, pointeur sur les sequences observees,
 *              flag semi-Markov cache.
 *
 *--------------------------------------------------------------*/

ostream& Semi_markov::spreadsheet_write(ostream &os , const Semi_markov_data *seq ,
                                        bool hidden) const

{
  register int i;
  int variable;
  Histogram **observation = 0;
  Sequence_characteristics *characteristics = 0;


  switch (hidden) {

  case false : {
    switch (type) {
    case 'o' :
      os << SEQ_word[SEQW_SEMI_MARKOV_CHAIN] << endl;
      break;
    case 'e' :
      os << SEQ_word[SEQW_EQUILIBRIUM_SEMI_MARKOV_CHAIN] << endl;
      break;
    }
    break;
  }

  case true : {
    switch (type) {
    case 'o' :
      os << SEQ_word[SEQW_HIDDEN_SEMI_MARKOV_CHAIN] << endl;
      break;
    case 'e' :
      os << SEQ_word[SEQW_EQUILIBRIUM_HIDDEN_SEMI_MARKOV_CHAIN] << endl;
      break;
    }
    break;
  }
  }

  // ecriture des parametres de la chaine de Markov

  spreadsheet_print(os);

  // ecriture des lois d'occupation des etats

  if ((seq) && (seq->type[0] == STATE)) {
    characteristics = seq->characteristics[0];
  }
  nonparametric_process[0]->spreadsheet_print(os , 0 , 0 , characteristics , forward);

  // ecriture des lois associees a chaque processus d'observation

  if (hidden) {
    os << "\n" << nb_output_process << "\t"
       << STAT_word[nb_output_process == 1 ? STATW_OUTPUT_PROCESS : STATW_OUTPUT_PROCESSES] << endl;
  }

  for (i = 1;i <= nb_output_process;i++) {
    os << "\n" << STAT_word[STATW_OUTPUT_PROCESS];

    if (hidden) {
      os << "\t" << i;

      if (nonparametric_process[i]) {
        os << "\t" << STAT_word[STATW_NONPARAMETRIC];
      }
      else {
        os << "\t" << STAT_word[STATW_PARAMETRIC];
      }
    }
    os << endl;

    if (seq) {
      switch (seq->type[0]) {
      case INT_VALUE :
        variable = i - 1;
        break;
      case STATE :
        variable = i;
        break;
      }

      if (seq->observation) {
        observation = seq->observation[variable];
      }
      if (seq->characteristics[variable]) {
        characteristics = seq->characteristics[variable];
      }
    }

    if (nonparametric_process[i]) {
      nonparametric_process[i]->spreadsheet_print(os , i , observation , characteristics);
    }
    else {
      parametric_process[i]->spreadsheet_print(os , observation);
    }
  }

  if (seq) {
    int nb_parameter = nb_parameter_computation(hidden ? MIN_PROBABILITY : 0.);
    double information , likelihood;


    // ecriture de l'histogramme des longueurs des sequences

    os << "\n" << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_HISTOGRAM] << "\t";
    seq->hlength->spreadsheet_characteristic_print(os);

    os << "\n\t" << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_HISTOGRAM] << endl;
    seq->hlength->spreadsheet_print(os);

    os << "\n" << SEQ_label[SEQL_CUMUL_LENGTH] << "\t" << seq->cumul_length << endl;

    // ecriture de la quantite d'information des sequences dans le cas iid

    information = seq->iid_information_computation();

    os << "\n" << SEQ_label[SEQL_IID_INFORMATION] << "\t" << information << "\t"
       << information / seq->cumul_length << endl;

    // ecriture des vraisemblances des sequences

    switch (hidden) {

    case false : {
      likelihood = seq->likelihood;

      os << "\n" << STAT_label[STATL_LIKELIHOOD] << "\t" << likelihood << "\t"
         << STAT_label[STATL_NORMALIZED] << "\t" << likelihood / seq->cumul_length << endl;
      break;
    }

    case true : {
      if (seq->likelihood != D_INF) {
        os << "\n" << SEQ_label[SEQL_STATE_SEQUENCES_LIKELIHOOD] << "\t" << seq->likelihood << "\t"
           << STAT_label[STATL_NORMALIZED] << "\t" << seq->likelihood / seq->cumul_length << endl;
      }

      likelihood = seq->hidden_likelihood;

      if (likelihood != D_INF) {
        os << "\n" << SEQ_label[SEQL_OBSERVED_SEQUENCES_LIKELIHOOD] << "\t" << likelihood << "\t"
           << STAT_label[STATL_NORMALIZED] << "\t" << likelihood / seq->cumul_length << endl;
      }
      break;
    }
    }

    if (likelihood != D_INF) {
      os << "\n" << nb_parameter << "\t" << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS] << "\t"
         << "2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[AIC] << ")\t"
         << 2 * (likelihood - nb_parameter) << endl;

      if (nb_parameter < seq->cumul_length - 1) {
        os << "\n" << nb_parameter << "\t" << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS] << "\t"
           << "2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[AICc] << ")\t"
           << 2 * (likelihood - (double)(nb_parameter * seq->cumul_length) /
              (double)(seq->cumul_length - nb_parameter - 1)) << endl;
      }

      os << "\n" << nb_parameter << "\t" << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS] << "\t"
         << "2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[BIC] << ")\t"
         << 2 * likelihood - nb_parameter * log((double)seq->cumul_length) << endl;

      os << "\n" << nb_parameter << "\t" << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS] << "\t"
         << "2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[BICc] << ")\t"
         << 2 * likelihood - penalty_computation(hidden , (hidden ? MIN_PROBABILITY : 0.)) << endl;
    }

    if (seq->likelihood != D_INF) {
      os << "\n" << nb_parameter << "\t" << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS] << "\t"
         << "2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[ICL] << ")\t"
         << 2 * seq->likelihood - nb_parameter * log((double)seq->cumul_length) << endl;

      os << "\n" << nb_parameter << "\t" << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS] << "\t"
         << "2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[BICc] << ")\t"
         << 2 * seq->likelihood - penalty_computation(hidden , MIN_PROBABILITY) << endl;
    }
  }

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Semi_markov dans un fichier au format tableur.
 *
 *  arguments : reference sur un objet Format_error, path.
 *
 *--------------------------------------------------------------*/

bool Semi_markov::spreadsheet_write(Format_error &error , const char *path) const

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
    spreadsheet_write(out_file , semi_markov_data);
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Sortie Gnuplot d'un objet Semi_markov et de la structure
 *  de donnees associee.
 *
 *  arguments : prefixe des fichiers, titre des figures,
 *              pointeur sur les sequences observees.
 *
 *--------------------------------------------------------------*/

bool Semi_markov::plot_write(const char *prefix , const char *title ,
                             const Semi_markov_data *seq) const

{
  bool status;
  register int i;
  int variable;
  Histogram *hlength = 0 , **observation = 0;
  Sequence_characteristics *characteristics = 0;


  if ((seq) && (seq->type[0] == STATE)) {
    characteristics = seq->characteristics[0];
    hlength = seq->hlength;
  }

  status = nonparametric_process[0]->plot_print(prefix , title , 0 , 0 , characteristics ,
                                                hlength , forward);

  if (status) {
    if (seq) {
      hlength = seq->hlength;
    }

    for (i = 1;i <= nb_output_process;i++) {
      if (seq) {
        switch (seq->type[0]) {
        case INT_VALUE :
          variable = i - 1;
          break;
        case STATE :
          variable = i;
          break;
        }

        if (seq->observation) {
          observation = seq->observation[variable];
        }
        if (seq->characteristics[variable]) {
          characteristics = seq->characteristics[variable];
        }
      }

      if (nonparametric_process[i]) {
        nonparametric_process[i]->plot_print(prefix , title , i , observation ,
                                             characteristics , hlength);
      }
      else {
        parametric_process[i]->plot_print(prefix , title , i , observation);
      }
    }
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Sortie Gnuplot d'un objet Semi_markov.
 *
 *  arguments : reference sur un objet Format_error, prefixe des fichiers,
 *              titre des figures.
 *
 *--------------------------------------------------------------*/

bool Semi_markov::plot_write(Format_error &error , const char *prefix ,
                             const char *title) const

{
  bool status = plot_write(prefix , title , semi_markov_data);

  error.init();

  if (!status) {
    error.update(STAT_error[STATR_FILE_PREFIX]);
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Fonctions pour la persistance.
 *
 *--------------------------------------------------------------*/

/* RWDEFINE_COLLECTABLE(Semi_markov , STATI_SEMI_MARKOV);


RWspace Semi_markov::binaryStoreSize() const

{
  register int i;
  RWspace size;


  size = Chain::binaryStoreSize() + sizeof(nb_output_process);

  for (i = 0;i <= nb_output_process;i++) {
    size += nonparametric_process[i]->binaryStoreSize();
  }

  if (semi_markov_data) {
    size += semi_markov_data->recursiveStoreSize();
  }

  return size;
}


void Semi_markov::restoreGuts(RWvistream &is)

{
  register int i;


  remove();

  Chain::restoreGuts(is);

  is >> nb_output_process;

  nonparametric_process = new Nonparametric_sequence_process*[nb_output_process + 1];
  for (i = 0;i <= nb_output_process;i++) {
    nonparametric_process[i] = new Nonparametric_sequence_process();
    nonparametric_process[i]->restoreGuts(is);
  }

  is >> semi_markov_data;
  if (semi_markov_data == RWnilCollectable) {
    semi_markov_data = 0;
  }
}


void Semi_markov::restoreGuts(RWFile &file)

{
  register int i;


  remove();

  Chain::restoreGuts(file);

  file.Read(nb_output_process);

  nonparametric_process = new Nonparametric_sequence_process*[nb_output_process + 1];
  for (i = 0;i <= nb_output_process;i++) {
    nonparametric_process[i] = new Nonparametric_sequence_process();
    nonparametric_process[i]->restoreGuts(file);
  }

  file >> semi_markov_data;
  if (semi_markov_data == RWnilCollectable) {
    semi_markov_data = 0;
  }
}


void Semi_markov::saveGuts(RWvostream &os) const

{
  register int i;


  Chain::saveGuts(os);

  os << nb_output_process;
  for (i = 0;i <= nb_output_process;i++) {
    nonparametric_process[i]->saveGuts(os);
  }

  if (semi_markov_data) {
    os << semi_markov_data;
  }
  else {
    os << RWnilCollectable;
  }
}


void Semi_markov::saveGuts(RWFile &file) const

{
  register int i;


  Chain::saveGuts(file);

  file.Write(nb_output_process);
  for (i = 0;i <= nb_output_process;i++) {
    nonparametric_process[i]->saveGuts(file);
  }

  if (semi_markov_data) {
    file << semi_markov_data;
  }
  else {
    file << RWnilCollectable;
  }
} */


/*--------------------------------------------------------------*
 *
 *  Calcul du nombre de parametres independants d'un objet Semi_markov.
 *
 *  argument : probabilite minimum.
 *
 *--------------------------------------------------------------*/

int Semi_markov::nb_parameter_computation(double min_probability) const

{
  register int i;
  int nb_parameter = Chain::nb_parameter_computation(min_probability);


  for (i = 0;i < nb_state;i++) {
    if (state_subtype[i] == SEMI_MARKOVIAN) {
      nb_parameter += nonparametric_process[0]->sojourn_time[i]->nb_parameter_computation();
      if (nonparametric_process[0]->sojourn_time[i]->inf_bound == 1) {
        nb_parameter--;
      }
    }
  }

  for (i = 1;i <= nb_output_process;i++) {
    if (nonparametric_process[i]) {
      nb_parameter += nonparametric_process[i]->nb_parameter_computation(min_probability);
    }
    else {
      nb_parameter += parametric_process[i]->nb_parameter_computation();
    }
  }

  return nb_parameter;
}


/*--------------------------------------------------------------*
 *
 *  Calcul d'une penalite adaptative.
 *
 *  arguments : flag semi-Markov cache, probabilite minimum.
 *
 *--------------------------------------------------------------*/

double Semi_markov::penalty_computation(bool hidden , double min_probability) const

{
  register int i , j , k;
  int nb_parameter , sample_size;
  double sum , *memory , *state_marginal;
  double penalty = 0.;


  if (semi_markov_data) {
    if (hidden) {
      memory = memory_computation();

      state_marginal = new double[nb_state];

      switch (type) {

      case 'o' : {
        sum = 0.;
        for (i = 0;i < nonparametric_process[0]->length->nb_value - 2;i++) {
          sum += (1. - nonparametric_process[0]->length->cumul[i + 1]);
        }
        for (i = 0;i < nb_state;i++) {
          memory[i] /= sum;
        }

        for (i = 0;i < nb_state;i++) {
          state_marginal[i] = 0.;
        }
        for (i = 0;i < nonparametric_process[0]->length->nb_value - 1;i++) {
          for (j = 0;j < nb_state;j++) {
            state_marginal[j] += nonparametric_process[0]->index_value->point[j][i] *
                                 (1. - nonparametric_process[0]->length->cumul[i]);
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

      case 'e' : {
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
        switch (hidden) {

        case false : {
          if (sample_size > 0) {
            penalty += nb_parameter * log((double)sample_size);
          }
          break;
        }

        case true : {
          if (memory[i] > 0.) {
            penalty += nb_parameter * log(memory[i] * semi_markov_data->cumul_length);
          }
          break;
        }
        }
      }
    }

    for (i = 0;i < nb_state;i++) {
      if (state_subtype[i] == SEMI_MARKOVIAN) {
        nb_parameter = nonparametric_process[0]->sojourn_time[i]->nb_parameter_computation();
        if (nonparametric_process[0]->sojourn_time[i]->inf_bound == 1) {
          nb_parameter--;
        }

        switch (hidden) {
        case false :
          penalty += nb_parameter * log((double)semi_markov_data->marginal[0]->frequency[i]);
          break;
        case true :
          penalty += nb_parameter * log(state_marginal[i] * semi_markov_data->cumul_length);
          break;
        }
      }
    }

    for (i = 1;i <= nb_output_process;i++) {
      if (nonparametric_process[i]) {
        for (j = 0;j < nb_state;j++) {
          nb_parameter = 0;
          for (k = 0;k < nonparametric_process[i]->nb_value;k++) {
            if (nonparametric_process[i]->observation[j]->mass[k] > min_probability) {
              nb_parameter++;
            }
          }

          nb_parameter--;

          if (nb_parameter > 0) {
            switch (hidden) {
              case false :
              penalty += nb_parameter * log((double)semi_markov_data->marginal[0]->frequency[j]);
              break;
            case true :
              penalty += nb_parameter * log(state_marginal[j] * semi_markov_data->cumul_length);
              break;
            }
          }
        }
      }

      else {
        for (j = 0;j < nb_state;j++) {
          nb_parameter = parametric_process[i]->observation[j]->nb_parameter_computation();

          switch (hidden) {
            case false :
            penalty += nb_parameter * log((double)semi_markov_data->marginal[0]->frequency[j]);
            break;
          case true :
            penalty += nb_parameter * log(state_marginal[j] * semi_markov_data->cumul_length);
            break;
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


/*--------------------------------------------------------------*
 *
 *  Constructeur par defaut de la classe Semi_markov_data.
 *
 *--------------------------------------------------------------*/

Semi_markov_data::Semi_markov_data()

{
  semi_markov = 0;
  chain_data = 0;
  likelihood = D_INF;
  hidden_likelihood = D_INF;
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Semi_markov_data.
 *
 *  arguments : nombre de variables, histogramme des longueurs des sequences,
 *              flag initialisation.
 *
 *--------------------------------------------------------------*/

Semi_markov_data::Semi_markov_data(int inb_variable , const Histogram &ihlength ,
                                   bool init_flag)
:Markovian_sequences(inb_variable , ihlength , init_flag)

{
  semi_markov = 0;
  chain_data = 0;
  likelihood = D_INF;
  hidden_likelihood = D_INF;
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Semi_markov_data.
 *
 *  arguments : nombre de sequences, longueurs des sequences, sequences,
 *              pointeur sur un objet Semi_markov.
 *
 *--------------------------------------------------------------*/

Semi_markov_data::Semi_markov_data(int inb_sequence , int *ilength , int ***isequence ,
                                   const Semi_markov &ismarkov)
:Markovian_sequences(ismarkov.nb_output_process + 1 , inb_sequence , ilength , isequence)

{
  type[0] = STATE;
  semi_markov = new Semi_markov(ismarkov , false);

  // extraction des caracteristiques des sequences simulees

  build_transition_count(semi_markov);
  build_observation_histogram();

  semi_markov->characteristic_computation(*this , true);

  // calcul de la vraisemblance

  likelihood = semi_markov->likelihood_computation(*this);
}


/*--------------------------------------------------------------*
 *
 *  Construction d'un objet Semi_markov_data a partir d'un objet Markovian_sequences.
 *
 *  arguments : reference sur un objet Markovian_sequences,
 *              ajout/suppression des histogrammes de temps de sejour initial.
 *
 *--------------------------------------------------------------*/

Semi_markov_data::Semi_markov_data(const Markovian_sequences &seq , bool initial_run_flag)
:Markovian_sequences(seq , 'c' , (initial_run_flag ? ADD_INITIAL_RUN : REMOVE_INITIAL_RUN))

{
  semi_markov = 0;
  chain_data = 0;
  likelihood = D_INF;
  hidden_likelihood = D_INF;
}


/*--------------------------------------------------------------*
 *
 *  Construction d'un objet Semi_markov_data a partir d'un objet Markovian_sequences
 *  avec ajout d'une variable.
 *
 *  arguments : reference sur un objet Markovian_sequences, indice de la variable.
 *
 *--------------------------------------------------------------*/

Semi_markov_data::Semi_markov_data(const Markovian_sequences &seq , int variable)
:Markovian_sequences(seq , 'a' , variable , DEFAULT)

{
  semi_markov = 0;
  chain_data = 0;
  likelihood = D_INF;
  hidden_likelihood = D_INF;
}


/*--------------------------------------------------------------*
 *
 *  Construction d'un objet Semi_markov_data a partir d'un objet Markovian_sequences
 *  avec ajout d'une variable.
 *
 *  arguments : reference sur un objet Markovian_sequences, indice de la variable,
 *              ajout/suppression des histogrammes de temps de sejour initial.
 *
 *--------------------------------------------------------------*/

Semi_markov_data::Semi_markov_data(const Markovian_sequences &seq , int variable , bool initial_run_flag)
:Markovian_sequences(seq , 'a' , variable , (initial_run_flag ? ADD_INITIAL_RUN : REMOVE_INITIAL_RUN))

{
  semi_markov = 0;
  chain_data = 0;
  likelihood = D_INF;
  hidden_likelihood = D_INF;
}


/*--------------------------------------------------------------*
 *
 *  Copie d'un objet Semi_markov_data.
 *
 *  arguments : reference sur un objet Semi_markov_data,
 *              flag copie de l'objet Semi_markov.
 *
 *--------------------------------------------------------------*/

void Semi_markov_data::copy(const Semi_markov_data &seq , bool model_flag)

{
  if ((model_flag) && (seq.semi_markov)) {
    semi_markov = new Semi_markov(*(seq.semi_markov) , false);
  }
  else {
    semi_markov = 0;
  }

  if (seq.chain_data) {
    chain_data = new Chain_data(*(seq.chain_data));
  }
  else {
    chain_data = 0;
  }

  likelihood = seq.likelihood;
  hidden_likelihood = seq.hidden_likelihood;
}


/*--------------------------------------------------------------*
 *
 *  Destructeur de la classe Semi_markov_data.
 *
 *--------------------------------------------------------------*/

Semi_markov_data::~Semi_markov_data()

{
  delete semi_markov;
  delete chain_data;
}


/*--------------------------------------------------------------*
 *
 *  Operateur d'assignement de la classe Semi_markov_data.
 *
 *  argument : reference sur un objet Semi_markov_data.
 *
 *--------------------------------------------------------------*/

Semi_markov_data& Semi_markov_data::operator=(const Semi_markov_data &seq)

{
  if (&seq != this) {
    delete semi_markov;
    delete chain_data;

    remove();
    Sequences::remove();

    Sequences::copy(seq);
    Markovian_sequences::copy(seq);
    copy(seq);
  }

  return *this;
}


/*--------------------------------------------------------------*
 *
 *  Extraction d'un histogramme.
 *
 *  arguments : reference sur un objet Format_error, type d'histogramme,
 *              variable, valeur.
 *
 *--------------------------------------------------------------*/

Distribution_data* Semi_markov_data::extract(Format_error &error , int type ,
                                             int variable , int value) const

{
  bool status = true;
  Distribution *pdist;
  Parametric *pparam;
  Histogram *phisto;
  Distribution_data *histo;


  histo = 0;
  error.init();

  phisto = 0;

  if (type == OBSERVATION) {
    if ((variable < 2) || (variable > nb_variable)) {
      status = false;
      error.update(STAT_error[STATR_VARIABLE_INDEX]);
    }

    else {
      variable--;

      if ((value < 0) || (value >= marginal[0]->nb_value)) {
        status = false;
        ostringstream error_message;
        error_message << STAT_label[STATL_STATE] << " " << value << " "
                      << SEQ_error[SEQR_NOT_PRESENT];
        error.update((error_message.str()).c_str());
      }

      else {
        phisto = observation[variable][value];

        if (phisto->nb_element == 0) {
          status = false;
          error.update(STAT_error[STATR_EMPTY_HISTOGRAM]);
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

      else if ((value < 0) || (value >= marginal[variable]->nb_value)) {
        status = false;
        ostringstream error_message;
        error_message << STAT_label[variable == 0 ? STATL_STATE : STATL_OUTPUT] << " "
                      << value << " " << SEQ_error[SEQR_NOT_PRESENT];
        error.update((error_message.str()).c_str());
      }

      if (status) {
        switch (type) {

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
            error.update(STAT_error[STATR_NON_EXISTING_HISTOGRAM]);
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
          error.update(STAT_error[STATR_EMPTY_HISTOGRAM]);
        }
      }
    }
  }

  if (status) {
    pdist = 0;
    pparam = 0;

    switch (type) {

    case OBSERVATION : {
      if (semi_markov->nonparametric_process[variable]) {
        pdist = semi_markov->nonparametric_process[variable]->observation[value];
      }
      else {
        pparam = semi_markov->parametric_process[variable]->observation[value];
      }
      break;
    }

    case FIRST_OCCURRENCE : {
      pdist = semi_markov->nonparametric_process[variable]->first_occurrence[value];
      break;
    }

    case RECURRENCE_TIME : {
      pdist = semi_markov->nonparametric_process[variable]->recurrence_time[value];
      break;
    }

    case SOJOURN_TIME : {
      pparam = semi_markov->nonparametric_process[variable]->sojourn_time[value];
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
      pdist = semi_markov->nonparametric_process[variable]->nb_run[value];
      break;
    }

    case NB_OCCURRENCE : {
      pdist = semi_markov->nonparametric_process[variable]->nb_occurrence[value];
      break;
    }
    }

    if (pdist) {
      histo = new Distribution_data(*phisto , pdist);
    }
    else {
      histo = new Distribution_data(*phisto , pparam);
    }
  }

  return histo;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Semi_markov_data.
 *
 *  arguments : stream, flag niveau de detail.
 *
 *--------------------------------------------------------------*/

ostream& Semi_markov_data::ascii_write(ostream &os , bool exhaustive) const

{
  if (semi_markov) {
    semi_markov->ascii_write(os , this , exhaustive , false ,
                             ::test_hidden(semi_markov->nb_output_process , semi_markov->nonparametric_process));
  }

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Semi_markov_data dans un fichier.
 *
 *  arguments : reference sur un objet Format_error, path,
 *              flag niveau de detail.
 *
 *--------------------------------------------------------------*/

bool Semi_markov_data::ascii_write(Format_error &error , const char *path ,
                                   bool exhaustive) const

{
  bool status = false;


  if (semi_markov) {
    ofstream out_file(path);

    error.init();

    if (!out_file) {
      status = false;
      error.update(STAT_error[STATR_FILE_NAME]);
    }

    else {
      status = true;
      semi_markov->ascii_write(out_file , this , exhaustive , true ,
                               ::test_hidden(semi_markov->nb_output_process , semi_markov->nonparametric_process));
    }
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Semi_markov_data dans un fichier au format tableur.
 *
 *  arguments : reference sur un objet Format_error, path.
 *
 *--------------------------------------------------------------*/

bool Semi_markov_data::spreadsheet_write(Format_error &error , const char *path) const

{
  bool status = false;


  if (semi_markov) {
    ofstream out_file(path);

    error.init();

    if (!out_file) {
      status = false;
      error.update(STAT_error[STATR_FILE_NAME]);
    }

    else {
      status = true;
      semi_markov->spreadsheet_write(out_file , this ,
                                     ::test_hidden(semi_markov->nb_output_process , semi_markov->nonparametric_process));
    }
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Sortie Gnuplot d'un objet Semi_markov_data.
 *
 *  arguments : reference sur un objet Format_error, prefixe des fichiers,
 *              titre des figures.
 *
 *--------------------------------------------------------------*/

bool Semi_markov_data::plot_write(Format_error &error , const char *prefix ,
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


/*--------------------------------------------------------------*
 *
 *  Fonctions pour la persistance.
 *
 *--------------------------------------------------------------*/

/* RWDEFINE_COLLECTABLE(Semi_markov_data , STATI_SEMI_MARKOV_DATA);


RWspace Semi_markov_data::binaryStoreSize() const

{
  RWspace size = Markovian_sequences::binaryStoreSize() + sizeof(likelihood) + sizeof(hidden_likelihood);

  size += sizeof(true);
  if (chain_data) {
    size += chain_data->binaryStoreSize();
  }

  if (semi_markov) {
    size += semi_markov->recursiveStoreSize();
  }

  return size;
}


void Semi_markov_data::restoreGuts(RWvistream &is)

{
  bool status;


  delete semi_markov;
  delete chain_data;

  Markovian_sequences::restoreGuts(is);

  is >> status;
  if (status) {
    chain_data = new Chain_data();
    chain_data->restoreGuts(is);
  }
  else {
    chain_data = 0;
  }

  is >> likelihood >> hidden_likelihood;

  is >> semi_markov;
  if (semi_markov == RWnilCollectable) {
    semi_markov = 0;
  }
}


void Semi_markov_data::restoreGuts(RWFile &file)

{
  bool status;


  delete semi_markov;
  delete chain_data;

  Markovian_sequences::restoreGuts(file);

  file.Read(status);
  if (status) {
    chain_data = new Chain_data();
    chain_data->restoreGuts(file);
  }
  else {
    chain_data = 0;
  }

  file.Read(likelihood);
  file.Read(hidden_likelihood);

  file >> semi_markov;
  if (semi_markov == RWnilCollectable) {
    semi_markov = 0;
  }
}


void Semi_markov_data::saveGuts(RWvostream &os) const

{
  Markovian_sequences::saveGuts(os);

  if (chain_data) {
    os << true;
    chain_data->saveGuts(os);
  }
  else {
    os << false;
  }

  os << likelihood << hidden_likelihood;

  if (semi_markov) {
    os << semi_markov;
  }
  else {
    os << RWnilCollectable;
  }
}


void Semi_markov_data::saveGuts(RWFile &file) const

{
  Markovian_sequences::saveGuts(file);

  if (chain_data) {
    file.Write(true);
    chain_data->saveGuts(file);
  }
  else {
    file.Write(false);
  }

  file.Write(likelihood);
  file.Write(hidden_likelihood);

  if (semi_markov) {
    file << semi_markov;
  }
  else {
    file << RWnilCollectable;
  }
} */
