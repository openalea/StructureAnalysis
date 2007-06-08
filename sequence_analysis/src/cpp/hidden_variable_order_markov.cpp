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
#include "stat_tool/stat_tools.h"
#include "stat_tool/curves.h"
#include "stat_tool/markovian.h"
#include "stat_tool/stat_label.h"
#include "sequences.h"
#include "variable_order_markov.h"
#include "hidden_variable_order_markov.h"
#include "sequence_label.h"

using namespace std;



/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Variable_order_markov.
 *
 *  arguments : pointeur sur un objet Variable_order_markov,
 *              nombre de processus d'observation, pointeurs sur des objets
 *              Nonparametric_process et Parametric_process, longueur des sequences.
 *
 *--------------------------------------------------------------*/

Variable_order_markov::Variable_order_markov(const Variable_order_markov *pmarkov , int inb_output_process ,
                                             Nonparametric_process **nonparametric_observation ,
                                             Parametric_process **parametric_observation , int length)

{
  register int i;
  int nb_terminal;


  memory_tree_completion(*pmarkov);

  switch (type) {

  case 'o' : {
    non_terminal_transition_probability_computation();
    break;
  }

  case 'e' : {
    nb_terminal = (nb_row - 1) * (nb_state - 1) / nb_state + 1;

    for (i = 1;i < nb_row;i++) {
      if (!child[i]) {
        initial[i] = 1. / (double)nb_terminal;
      }
      else {
        initial[i] = 0.;
      }
    }

    initial_probability_computation();
    break;
  }
  }

  nb_output_process = inb_output_process;

  nonparametric_process = new Nonparametric_sequence_process*[nb_output_process + 1];
  nonparametric_process[0] = new Nonparametric_sequence_process(nb_state , nb_state);
  parametric_process = new Parametric_process*[nb_output_process + 1];
  parametric_process[0] = 0;

  for (i = 1;i <= nb_output_process;i++) {
    if (nonparametric_observation[i - 1]) {
      nonparametric_process[i] = new Nonparametric_sequence_process(*nonparametric_observation[i - 1]);
      parametric_process[i] = 0;
    }
    else {
      nonparametric_process[i] = 0;
      parametric_process[i] = new Parametric_process(*parametric_observation[i - 1]);
    }
  }

  characteristic_computation(length , true);
}


/*--------------------------------------------------------------*
 *
 *  Destructeur de la classe Hidden_variable_order_markov.
 *
 *--------------------------------------------------------------*/

Hidden_variable_order_markov::~Hidden_variable_order_markov() {}


/*--------------------------------------------------------------*
 *
 *  Application d'un seuil sur les parametres d'une chaine de Markov cachee.
 *
 *  argument : probabilite minimum.
 *
 *--------------------------------------------------------------*/

Hidden_variable_order_markov* Hidden_variable_order_markov::thresholding(double min_probability) const

{
  register int i;
  Hidden_variable_order_markov *hmarkov;


  hmarkov = new Hidden_variable_order_markov(*this , false);
  hmarkov->threshold_application(min_probability);

  for (i = 1;i <= nb_output_process;i++) {
    if (nonparametric_process[i]) {
      nonparametric_process[i]->thresholding(min_probability);
    }
  }

  return hmarkov;
}


/*--------------------------------------------------------------*
 *
 *  Construction d'un objet Hidden_variable_order_markov a partir d'un fichier.
 *
 *  arguments : reference sur un objet Format_error, path, longueur des sequences,
 *              seuil sur les fonctions de repartition des lois parametriques.
 *
 *--------------------------------------------------------------*/

Hidden_variable_order_markov* hidden_variable_order_markov_ascii_read(Format_error &error ,
                                                                      const char *path , int length ,
                                                                      double cumul_threshold)

{
  RWLocaleSnapshot locale("en");
  RWCString buffer , token;
  size_t position;
  char type = 'v';
  bool status , lstatus , nonparametric;
  register int i;
  int line , nb_output_process , index;
  long value;
  const Variable_order_markov *imarkov;
  Nonparametric_process **nonparametric_observation;
  Parametric_process **parametric_observation;
  Hidden_variable_order_markov *hmarkov;
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

        // test mot cle (EQUILIBRIUM) HIDDEN_MARKOV_CHAIN

        if (i == 0) {
          if (token == SEQ_word[SEQW_HIDDEN_MARKOV_CHAIN]) {
            type = 'o';
          }
          else if (token == SEQ_word[SEQW_EQUILIBRIUM_HIDDEN_MARKOV_CHAIN]) {
            type = 'e';
          }
          else {
            status = false;
            ostringstream correction_message;
            correction_message << SEQ_word[SEQW_HIDDEN_MARKOV_CHAIN] << " or "
                               << SEQ_word[SEQW_EQUILIBRIUM_HIDDEN_MARKOV_CHAIN];
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

      // analyse du format et lecture de la chaine de Markov d'ordre variable

      imarkov = variable_order_markov_parsing(error , in_file , line , type);

      // analyse du format et lecture des lois d'observation

      if (imarkov) {
        nb_output_process = I_DEFAULT;

        nonparametric_observation = 0;
        parametric_observation = 0;

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
          nonparametric_observation = new Nonparametric_process*[nb_output_process];
          parametric_observation = new Parametric_process*[nb_output_process];
          for (i = 0;i < nb_output_process;i++) {
            nonparametric_observation[i] = 0;
            parametric_observation[i] = 0;
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
                                                                           ((Chain*)imarkov)->nb_state , true);
                if (!nonparametric_observation[index - 1]) {
                  status = false;
                }
                break;
              }

              case false : {
                parametric_observation[index - 1] = observation_parsing(error , in_file , line ,
                                                                        ((Chain*)imarkov)->nb_state ,
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
            hmarkov = new Hidden_variable_order_markov(imarkov , nb_output_process ,
                                                       nonparametric_observation ,
                                                       parametric_observation , length);

#           ifdef DEBUG
            hmarkov->ascii_write(cout);
#           endif

          }

          delete imarkov;

          for (i = 0;i < nb_output_process;i++) {
            delete nonparametric_observation[i];
            delete parametric_observation[i];
          }
          delete [] nonparametric_observation;
          delete [] parametric_observation;
        }
      }
    }
  }

  return hmarkov;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Hidden_variable_order_markov dans un fichier.
 *
 *  arguments : stream, flag niveau de detail.
 *
 *--------------------------------------------------------------*/

ostream& Hidden_variable_order_markov::ascii_write(ostream &os , bool exhaustive) const

{
  Variable_order_markov::ascii_write(os , markov_data , exhaustive , false , true);

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Hidden_variable_order_markov dans un fichier.
 *
 *  arguments : reference sur un objet Format_error, path,
 *              flag niveau de detail.
 *
 *--------------------------------------------------------------*/

bool Hidden_variable_order_markov::ascii_write(Format_error &error , const char *path ,
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
    Variable_order_markov::ascii_write(out_file , markov_data , exhaustive , true , true);
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Hidden_variable_order_markov dans un fichier au format tableur.
 *
 *  arguments : reference sur un objet Format_error, path.
 *
 *--------------------------------------------------------------*/

bool Hidden_variable_order_markov::spreadsheet_write(Format_error &error ,
                                                     const char *path) const

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
    Variable_order_markov::spreadsheet_write(out_file , markov_data , true);
  }

  return status;
}
