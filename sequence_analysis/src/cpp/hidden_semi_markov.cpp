/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       V-Plants: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2010 CIRAD/INRIA Virtual Plants
 *
 *       File author(s): Y. Guedon (yann.guedon@cirad.fr)
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
#include <sstream>
#include "tool/rw_tokenizer.h"
#include "tool/rw_cstring.h"
#include "tool/rw_locale.h"

#include "stat_tool/stat_tools.h"
#include "stat_tool/curves.h"
#include "stat_tool/markovian.h"
#include "stat_tool/stat_label.h"
#include "sequences.h"
#include "semi_markov.h"
#include "hidden_semi_markov.h"
#include "sequence_label.h"

using namespace std;



/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe SemiMarkov.
 *
 *  arguments : pointeur sur un objet Chain et sur un objet NonparametricSequenceProcess,
 *              nombre de processus d'observation, pointeurs sur
 *              des objets NonparametricProcess, longueur des sequences,
 *              flag sur le calcul des lois de comptage.
 *
 *--------------------------------------------------------------*/

SemiMarkov::SemiMarkov(const Chain *pchain , const NonparametricSequenceProcess *poccupancy ,
                       int inb_output_process , NonparametricProcess **pobservation ,
                       int length , bool counting_flag)
:Chain(*pchain)

{
  register int i;


  nb_iterator = 0;
  semi_markov_data = NULL;

  nb_output_process = inb_output_process;

  nonparametric_process = new NonparametricSequenceProcess*[nb_output_process + 1];
  nonparametric_process[0] = new NonparametricSequenceProcess(*poccupancy);
  for (i = 1;i <= nb_output_process;i++) {
    nonparametric_process[i] = new NonparametricSequenceProcess(*pobservation[i - 1]);
  }

  parametric_process = NULL;

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
      forward[i] = NULL;
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
 *  Constructeur de la classe SemiMarkov.
 *
 *  arguments : pointeur sur un objet Chain et sur un objet NonparametricSequenceProcess,
 *              nombre de processus d'observation, pointeurs sur des objets
 *              NonparametricProcess et DiscreteParametricProcess, longueur des sequences,
 *              flag sur le calcul des lois de comptage.
 *
 *--------------------------------------------------------------*/

SemiMarkov::SemiMarkov(const Chain *pchain , const NonparametricSequenceProcess *poccupancy ,
                       int inb_output_process , NonparametricProcess **nonparametric_observation ,
                       DiscreteParametricProcess **parametric_observation , int length , bool counting_flag)
:Chain(*pchain)

{
  register int i;


  nb_iterator = 0;
  semi_markov_data = NULL;

  nb_output_process = inb_output_process;

  nonparametric_process = new NonparametricSequenceProcess*[nb_output_process + 1];
  nonparametric_process[0] = new NonparametricSequenceProcess(*poccupancy);
  parametric_process = new DiscreteParametricProcess*[nb_output_process + 1];
  parametric_process[0] = NULL;

  for (i = 1;i <= nb_output_process;i++) {
    if (nonparametric_observation[i - 1]) {
      nonparametric_process[i] = new NonparametricSequenceProcess(*nonparametric_observation[i - 1]);
      parametric_process[i] = NULL;
    }
    else {
      nonparametric_process[i] = NULL;
      parametric_process[i] = new DiscreteParametricProcess(*parametric_observation[i - 1]);
    }
  }

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
      forward[i] = NULL;
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
 *  Destructeur de la classe HiddenSemiMarkov.
 *
 *--------------------------------------------------------------*/

HiddenSemiMarkov::~HiddenSemiMarkov() {}


/*--------------------------------------------------------------*
 *
 *  Application d'un seuil sur les parametres d'une semi-chaine de Markov cachee.
 *
 *  argument : probabilite minimum.
 *
 *--------------------------------------------------------------*/

HiddenSemiMarkov* HiddenSemiMarkov::thresholding(double min_probability) const

{
  register int i;
  HiddenSemiMarkov *hsmarkov;


  hsmarkov = new HiddenSemiMarkov(*this , false , false);
  hsmarkov->Chain::thresholding(min_probability , true);

  for (i = 1;i <= hsmarkov->nb_output_process;i++) {
    if (hsmarkov->nonparametric_process[i]) {
      hsmarkov->nonparametric_process[i]->thresholding(min_probability);
    }
  }

  return hsmarkov;
}


/*--------------------------------------------------------------*
 *
 *  Construction d'un objet HiddenSemiMarkov a partir d'un fichier.
 *
 *  arguments : reference sur un objet StatError, path, longueur des sequences,
 *              flag sur le calcul des lois de comptage, seuil sur les fonctions de repartition
 *              des lois parametriques, flag format des processus d'observation.
 *
 *--------------------------------------------------------------*/

HiddenSemiMarkov* hidden_semi_markov_ascii_read(StatError &error , const char *path ,
                                                int length , bool counting_flag ,
                                                double cumul_threshold , bool old_format)

{
  RWLocaleSnapshot locale("en");
  RWCString buffer , token;
  size_t position;
  char type = 'v';
  bool status , lstatus , nonparametric;
  register int i;
  int line , nb_output_process , index;
  long value;
  const Chain *chain;
  const NonparametricSequenceProcess *occupancy;
  NonparametricProcess **nonparametric_observation;
  DiscreteParametricProcess **parametric_observation;
  HiddenSemiMarkov *hsmarkov;
  ifstream in_file(path);


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

        // test mot cle (EQUILIBRIUM) HIDDEN_SEMI-MARKOV_CHAIN

        if (i == 0) {
          if (token == SEQ_word[SEQW_HIDDEN_SEMI_MARKOV_CHAIN]) {
            type = 'o';
          }
          else if (token == SEQ_word[SEQW_EQUILIBRIUM_HIDDEN_SEMI_MARKOV_CHAIN]) {
            type = 'e';
          }
          else {
            status = false;
            ostringstream correction_message;
            correction_message << SEQ_word[SEQW_HIDDEN_SEMI_MARKOV_CHAIN] << " or "
                               << SEQ_word[SEQW_EQUILIBRIUM_HIDDEN_SEMI_MARKOV_CHAIN];
            error.correction_update(STAT_parsing[STATP_KEY_WORD] , (correction_message.str()).c_str() , line);
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

      chain = chain_parsing(error , in_file , line , type);

      if (chain) {

        // analyse du format et lecture des lois d'occupation de etats

        occupancy = occupancy_parsing(error , in_file , line , *chain , cumul_threshold);
        if (!occupancy) {
          status = false;
        }

        // analyse du format et lecture des lois d'observation

        switch (old_format) {

        case false : {
          nb_output_process = I_DEFAULT;

          nonparametric_observation = NULL;
          parametric_observation = NULL;

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
            nonparametric_observation = new NonparametricProcess*[nb_output_process];
            parametric_observation = new DiscreteParametricProcess*[nb_output_process];
            for (i = 0;i < nb_output_process;i++) {
              nonparametric_observation[i] = NULL;
              parametric_observation[i] = NULL;
            }

            index = 0;

            while (buffer.readLine(in_file , false)) {
              line++;

#             ifdef DEBUG
              cout << line << "  " << buffer << endl;
#             endif

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
                  nonparametric = true;

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

                // test mot cle NONPARAMETRIC / PARAMETRIC

                case 3 : {
                  if (token == STAT_word[STATW_NONPARAMETRIC]) {
                    nonparametric = true;
                  }
                  else if (token == STAT_word[STATW_PARAMETRIC]) {
                    nonparametric = false;
                  }
                  else {
                    status = false;
                    ostringstream correction_message;
                    correction_message << STAT_word[STATW_NONPARAMETRIC] << " or "
                                       << STAT_word[STATW_PARAMETRIC];
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

                switch (nonparametric) {

                case true : {
                  nonparametric_observation[index - 1] = observation_parsing(error , in_file , line ,
                                                                             chain->nb_state , true);
                  if (!nonparametric_observation[index - 1]) {
                    status = false;
                  }
                  break;
                }

                case false : {
                  parametric_observation[index - 1] = observation_parsing(error , in_file , line ,
                                                                          chain->nb_state ,
                                                                          cumul_threshold);
                  if (!parametric_observation[index - 1]) {
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
              hsmarkov = new HiddenSemiMarkov(chain , occupancy , nb_output_process ,
                                              nonparametric_observation , parametric_observation ,
                                              length , counting_flag);
            }

            for (i = 0;i < nb_output_process;i++) {
              delete nonparametric_observation[i];
              delete parametric_observation[i];
            }
            delete [] nonparametric_observation;
            delete [] parametric_observation;
          }
          break;
        }

        case true : {
          nonparametric_observation = old_observation_parsing(error , in_file , line ,
                                                              chain->nb_state ,
                                                              nb_output_process);

          if (nonparametric_observation) {
            if (status) {
              hsmarkov = new HiddenSemiMarkov(chain , occupancy , nb_output_process ,
                                              nonparametric_observation , length , counting_flag);
            }

            for (i = 0;i < nb_output_process;i++) {
              delete nonparametric_observation[i];
            }
            delete [] nonparametric_observation;
          }
          break;
        }
        }

        delete chain;
        delete occupancy;
      }
    }
  }

  return hsmarkov;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet HiddenSemiMarkov dans un fichier.
 *
 *  arguments : stream, flag niveau de detail.
 *
 *--------------------------------------------------------------*/

ostream& HiddenSemiMarkov::ascii_write(ostream &os , bool exhaustive) const

{
  SemiMarkov::ascii_write(os , semi_markov_data , exhaustive ,
                          false , true);

//  os << "\nEnd state: " << end_state() << endl;

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet HiddenSemiMarkov dans un fichier.
 *
 *  arguments : reference sur un objet StatError, path,
 *              flag niveau de detail.
 *
 *--------------------------------------------------------------*/

bool HiddenSemiMarkov::ascii_write(StatError &error , const char *path ,
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
    SemiMarkov::ascii_write(out_file , semi_markov_data , exhaustive ,
                            true , true);
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet HiddenSemiMarkov dans un fichier au format tableur.
 *
 *  arguments : reference sur un objet StatError, path.
 *
 *--------------------------------------------------------------*/

bool HiddenSemiMarkov::spreadsheet_write(StatError &error , const char *path) const

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
    SemiMarkov::spreadsheet_write(out_file , semi_markov_data , true);
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Recherche etat de fin.
 *
 *--------------------------------------------------------------*/

int HiddenSemiMarkov::end_state() const

{
  register int i , j , k;
  int end_state = I_DEFAULT , output;


  for (i = nb_state - 1;i >= 0;i--) {
    if (state_type[i] == 'a') {

#     ifdef DEBUG
      cout << "\nstate: " << i << " | ";
#     endif

      end_state = i;

      for (j = 1;j <= nb_output_process;j++) {
        if (nonparametric_process[j]) {
          for (k = nonparametric_process[j]->observation[i]->offset;k < nonparametric_process[j]->observation[i]->nb_value;k++) {
            if (nonparametric_process[j]->observation[i]->mass[k] == 1.) {
              output = k;

#             ifdef DEBUG
              cout << "output: " << output << " | ";
#             endif

              break;
            }
          }

          if (k < nonparametric_process[j]->observation[i]->nb_value) {
            for (k = 0;k < nb_state;k++) {
              if ((k != i) && (output >= nonparametric_process[j]->observation[k]->offset) &&
                  (output < nonparametric_process[j]->observation[k]->nb_value) &&
                  (nonparametric_process[j]->observation[k]->mass[output] > 0.)) {
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
          for (k = parametric_process[j]->observation[i]->offset;k < parametric_process[j]->observation[i]->nb_value;k++) {
            if (parametric_process[j]->observation[i]->mass[k] == 1.) {
              output = k;
              break;
            }
          }

          if (k < parametric_process[j]->observation[i]->nb_value) {
            for (k = 0;k < nb_state;k++) {
              if ((k != i) && (output >= parametric_process[j]->observation[k]->offset) &&
                  (output < parametric_process[j]->observation[k]->nb_value) &&
                  (parametric_process[j]->observation[k]->mass[output] > 0.)) {
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
