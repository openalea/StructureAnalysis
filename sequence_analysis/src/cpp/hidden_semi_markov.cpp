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
 *       $Id: hidden_semi_markov.cpp 18050 2015-04-23 09:43:10Z guedon $
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
#include <vector>

#include "tool/rw_tokenizer.h"
#include "tool/rw_cstring.h"
#include "tool/rw_locale.h"

#include "stat_tool/stat_tools.h"
#include "stat_tool/curves.h"
#include "stat_tool/distribution.h"
#include "stat_tool/markovian.h"
#include "stat_tool/vectors.h"
#include "stat_tool/distance_matrix.h"
#include "stat_tool/stat_label.h"

#include "sequences.h"
#include "semi_markov.h"
#include "hidden_semi_markov.h"
#include "sequence_label.h"

using namespace std;
using namespace stat_tool;


namespace sequence_analysis {



/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe SemiMarkov.
 *
 *  arguments : pointeur sur un objet Chain et sur un objet CategoricalSequenceProcess,
 *              nombre de processus d'observation, pointeurs sur
 *              des objets CategoricalProcess, longueur des sequences,
 *              flag sur le calcul des lois de comptage.
 *
 *--------------------------------------------------------------*/

SemiMarkov::SemiMarkov(const Chain *pchain , const CategoricalSequenceProcess *poccupancy ,
                       int inb_output_process , CategoricalProcess **pobservation ,
                       int length , bool counting_flag)
:SemiMarkovChain(pchain , poccupancy)

{
  register int i;


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

  state_subtype = new int[nb_state];
  forward = new Forward*[nb_state];

  for (i = 0;i < nb_state;i++) {
    state_subtype[i] = (state_process->sojourn_time[i] ? SEMI_MARKOVIAN : MARKOVIAN);

    if ((state_subtype[i] == SEMI_MARKOVIAN) && (state_type[i] == 'r')) {
      forward[i] = new Forward(*(state_process->sojourn_time[i]));
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


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe SemiMarkov.
 *
 *  arguments : pointeur sur un objet Chain et sur un objet CategoricalSequenceProcess,
 *              nombre de processus d'observation, pointeurs sur des objets
 *              CategoricalProcess, DiscreteParametricProcess et
 *              ContinuousParametricProcess, longueur des sequences,
 *              flag sur le calcul des lois de comptage.
 *
 *--------------------------------------------------------------*/

SemiMarkov::SemiMarkov(const Chain *pchain , const CategoricalSequenceProcess *poccupancy ,
                       int inb_output_process , CategoricalProcess **categorical_observation ,
                       DiscreteParametricProcess **discrete_parametric_observation ,
                       ContinuousParametricProcess **continuous_parametric_observation ,
                       int length , bool counting_flag)
:SemiMarkovChain(pchain , poccupancy)

{
  register int i;


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

  state_subtype = new int[nb_state];
  forward = new Forward*[nb_state];

  for (i = 0;i < nb_state;i++) {
    state_subtype[i] = (state_process->sojourn_time[i] ? SEMI_MARKOVIAN : MARKOVIAN);

    if ((state_subtype[i] == SEMI_MARKOVIAN) && (state_type[i] == 'r')) {
      forward[i] = new Forward(*(state_process->sojourn_time[i]));
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

  for (i = 0;i < hsmarkov->nb_output_process;i++) {
    if (hsmarkov->categorical_process[i]) {
      hsmarkov->categorical_process[i]->thresholding(min_probability);
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
  bool status , lstatus;
  register int i;
  int line , nb_output_process , output_process_type , index;
  long value;
  const Chain *chain;
  const CategoricalSequenceProcess *occupancy;
  CategoricalProcess **categorical_observation;
  DiscreteParametricProcess **discrete_parametric_observation;
  ContinuousParametricProcess **continuous_parametric_observation;
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

          categorical_observation = NULL;
          discrete_parametric_observation = NULL;
          continuous_parametric_observation = NULL;

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
            categorical_observation = new CategoricalProcess*[nb_output_process];
            discrete_parametric_observation = new DiscreteParametricProcess*[nb_output_process];
            continuous_parametric_observation = new ContinuousParametricProcess*[nb_output_process];

            for (i = 0;i < nb_output_process;i++) {
              categorical_observation[i] = NULL;
              discrete_parametric_observation[i] = NULL;
              continuous_parametric_observation[i] = NULL;
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

                // test mot cle CATEGORICAL / DISCRETE_PARAMETRIC / CONTINUOUS_PARAMETRIC

                case 3 : {
                  if ((token == STAT_word[STATW_CATEGORICAL]) ||
                      (token == STAT_word[STATW_NONPARAMETRIC])) {
                    output_process_type = CATEGORICAL_PROCESS;
                  }
                  else if ((token == STAT_word[STATW_DISCRETE_PARAMETRIC]) ||
                           (token == STAT_word[STATW_PARAMETRIC])) {
                    output_process_type = DISCRETE_PARAMETRIC;
                  }
                  else if (token == STAT_word[STATW_CONTINUOUS_PARAMETRIC]) {
                    output_process_type = CONTINUOUS_PARAMETRIC;
                  }
                  else {
                    output_process_type = CATEGORICAL_PROCESS - 1;
                    status = false;
                    ostringstream correction_message;
                    correction_message << STAT_word[STATW_CATEGORICAL] << " or "
                                       << STAT_word[STATW_DISCRETE_PARAMETRIC] << " or "
                                       << STAT_word[STATW_CONTINUOUS_PARAMETRIC];
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

                switch (output_process_type) {

                case CATEGORICAL_PROCESS : {
                  categorical_observation[index - 1] = categorical_observation_parsing(error , in_file , line ,
                                                                                       chain->nb_state ,
                                                                                       HIDDEN_MARKOV , true);
/*                  categorical_observation[index - 1] = categorical_observation_parsing(error , in_file , line , pour les donnees de suivi de croissance manguier
                                                                                       chain->nb_state ,
                                                                                       HIDDEN_MARKOV , false); */
                  if (!categorical_observation[index - 1]) {
                    status = false;
                  }
                  break;
                }

                case DISCRETE_PARAMETRIC : {
                  discrete_parametric_observation[index - 1] = discrete_observation_parsing(error , in_file , line ,
                                                                                            chain->nb_state ,
                                                                                            HIDDEN_MARKOV ,
                                                                                            cumul_threshold);
                  if (!discrete_parametric_observation[index - 1]) {
                    status = false;
                  }
                  break;
                }

                case CONTINUOUS_PARAMETRIC : {
                  continuous_parametric_observation[index - 1] = continuous_observation_parsing(error , in_file , line ,
                                                                                                chain->nb_state ,
                                                                                                HIDDEN_MARKOV ,
                                                                                                LINEAR_MODEL);
                  if (!continuous_parametric_observation[index - 1]) {
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
          break;
        }

        case true : {
          categorical_observation = old_categorical_observation_parsing(error , in_file , line ,
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





