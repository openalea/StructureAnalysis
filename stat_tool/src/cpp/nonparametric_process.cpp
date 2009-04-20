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


#include <iomanip>

#include "tool/rw_tokenizer.h"
#include "tool/rw_cstring.h"
#include "tool/rw_locale.h"

#include "stat_tools.h"
#include "markovian.h"
#include "stat_label.h"

using namespace std;

extern int column_width(int);
extern int column_width(int nb_value , const double *value , double scale = 1.);
extern char* label(const char*);


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Nonparametric_process.
 *
 *  arguments : nombre d'etats, nombre de valeurs,
 *              flag sur les lois d'observation.
 *
 *--------------------------------------------------------------*/

Nonparametric_process::Nonparametric_process(int inb_state , int inb_value , int observation_flag)

{
  nb_state = inb_state;
  nb_value = inb_value;

  if (observation_flag) {
    register int i;


    observation = new Distribution*[nb_state];
    for (i = 0;i < nb_state;i++) {
      observation[i] = new Distribution(nb_value);
    }
  }

  else {
    observation = 0;
  }
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Nonparametric_process.
 *
 *  arguments : nombre d'etats, nombre de valeurs,
 *              probabilites d'observation.
 *
 *--------------------------------------------------------------*/

Nonparametric_process::Nonparametric_process(int inb_state , int inb_value ,
                                             double **observation_probability)

{
  register int i , j;


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
}
/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Nonparametric_process.
 *
 *  arguments : nombre d'etats, lois d'observation.
 *
 *--------------------------------------------------------------*/

Nonparametric_process::Nonparametric_process(int inb_state , Distribution **pobservation)

{
  register int i;

  nb_value = 0;
  nb_state = inb_state;

  observation = new Distribution*[nb_state];
  for (i = 0;i < nb_state;i++) {
    observation[i] = new Distribution(*pobservation[i]);

    observation[i]->cumul_computation();

    observation[i]->max_computation();
    nb_value = max(nb_value, observation[i]->nb_value);
  }  

}


/*--------------------------------------------------------------*
 *
 *  Copie d'un objet Nonparametric_process.
 *
 *  argument : reference sur un objet Nonparametric_process.
 *
 *--------------------------------------------------------------*/

void Nonparametric_process::copy(const Nonparametric_process &process)

{
  nb_state = process.nb_state;
  nb_value = process.nb_value;

  if (process.observation) {
    register int i;


    observation = new Distribution*[nb_state];
    for (i = 0;i < nb_state;i++) {
      observation[i] = new Distribution(*(process.observation[i]));
    }
  }

  else {
    observation = 0;
  }
}


/*--------------------------------------------------------------*
 *
 *  Copie d'un objet Nonparametric_process avec ajout d'un etat.
 *
 *  arguments : reference sur un objet Nonparametric_process, etat de reference.
 *
 *--------------------------------------------------------------*/

void Nonparametric_process::add_state(const Nonparametric_process &process , int state)

{
  nb_state = process.nb_state + 1;
  nb_value = process.nb_value;

  if (process.observation) {
    register int i;


    observation = new Distribution*[nb_state];
    for (i = 0;i <= state;i++) {
      observation[i] = new Distribution(*(process.observation[i]));
    }
    for (i = state + 1;i < nb_state;i++) {
      observation[i] = new Distribution(*(process.observation[i - 1]));
    }
  }

  else {
    observation = 0;
  }
}


/*--------------------------------------------------------------*
 *
 *  Constructeur par copie de la classe Nonparametric_process.
 *
 *  arguments : reference sur un objet Nonparametric_process,
 *              type de manipulation ('c' : copy, 's' : state), etat de reference.
 *
 *--------------------------------------------------------------*/

Nonparametric_process::Nonparametric_process(const Nonparametric_process &process , char manip , int state)

{
  switch (manip) {
  case 'c' :
    copy(process);
    break;
  case 's' :
    add_state(process , state);
    break;
  }
}


/*--------------------------------------------------------------*
 *
 *  Destruction des champs d'un objet Nonparametric_process.
 *
 *--------------------------------------------------------------*/

void Nonparametric_process::remove()

{
  if (observation) {
    register int i;


    for (i = 0;i < nb_state;i++) {
      delete observation[i];
    }
    delete [] observation;

    observation = 0;
  }
}


/*--------------------------------------------------------------*
 *
 *  Destructeur de la classe Nonparametric_process.
 *
 *--------------------------------------------------------------*/

Nonparametric_process::~Nonparametric_process()

{
  remove();
}


/*--------------------------------------------------------------*
 *
 *  Operateur d'assignement de la classe Nonparametric_process.
 *
 *  argument : reference sur un objet Nonparametric_process.
 *
 *--------------------------------------------------------------*/

Nonparametric_process& Nonparametric_process::operator=(const Nonparametric_process &process)

{
  if (&process != this) {
    remove();
    copy(process);
  }

  return *this;
}


/*--------------------------------------------------------------*
 *
 *  test recouvrement entre les valeurs observees dans chaque etat.
 *
 *--------------------------------------------------------------*/

bool Nonparametric_process::test_hidden() const

{
  bool hidden = false;
  register int i , j;
  int nb_occur;


  for (i = 0;i < nb_value;i++) {
    nb_occur = 0;
    for (j = 0;j < nb_state;j++) {
      if (observation[j]->mass[i] > 0.) {
        nb_occur++;
      }
    }

    if (nb_occur > 1) {
      hidden = true;
      break;
    }
  }

  return hidden;
}


/*--------------------------------------------------------------*
 *
 *  Application d'un seuil sur les probabilites d'observation.
 *
 *  argument : probabilite minimum.
 *
 *--------------------------------------------------------------*/

void Nonparametric_process::thresholding(double min_probability)

{
  bool stop;
  register int i , j;
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


/*--------------------------------------------------------------*
 *
 *  Application d'une permutation des etats. La validite de la 
 *  permutation doit etre verifiee par la procedure appelante.
 *
 *--------------------------------------------------------------*/

void Nonparametric_process::state_permutation(int *perm) const
{
  register int i;
  Distribution **pobservation= new Distribution*[nb_state];
  
  for (i= 0; i < nb_state; i++)
    pobservation[perm[i]]= observation[i];
  for (i= 0; i < nb_state; i++)
    observation[i]= pobservation[i];
  delete [] pobservation;
  pobservation= NULL;

}

/*--------------------------------------------------------------*
 *
 *  Analyse du format des lois d'observation.
 *
 *  arguments : reference sur un objet Format_error, stream,
 *              reference sur l'indice de la ligne lue, nombre d'etats,
 *              nombre de processus d'observation, flag sur le recouvrement.
 *
 *--------------------------------------------------------------*/

Nonparametric_process* observation_parsing(Format_error &error , ifstream &in_file ,
                                           int &line , int nb_state , bool hidden)

{
  RWLocaleSnapshot locale("en");
  RWCString buffer , token;
  size_t position;
  enum{PSTATE , POUTPUT};
  bool status = true , lstatus , defined_output[NB_OUTPUT];
  register int i , j;
  int type , state , output , nb_value;
  long value;
  double proba , cumul , **observation_probability;
  Nonparametric_process *process;


  process = 0;

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

      case 0 : {

        // test mot cle STATE

        if (token == STAT_word[STATW_STATE]) {
          type = PSTATE;

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

        // test mot cle OUTPUT

        else if (token == STAT_word[STATW_OUTPUT]) {
          type = POUTPUT;
        }

        else {
          status = false;
          error.update(STAT_parsing[STATP_KEY_WORD] , line , i + 1);
        }

        break;
      }

      // test indice de l'etat / de l'observation

      case 1 : {
        switch (type) {

        case PSTATE : {
          lstatus = locale.stringToNum(token , &value);
          if ((lstatus) && ((value != state) || (value >= nb_state))) {
            lstatus = false;
          }

          if (!lstatus) {
            status = false;
            error.update(STAT_parsing[STATP_STATE_INDEX] , line , i + 1);
          }
          break;
        }

        case POUTPUT : {
          lstatus = locale.stringToNum(token , &value);
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

        // test mot cle OBSERVATION_DISTRIBUTION

        case PSTATE : {
          if (token != STAT_word[STATW_OBSERVATION_DISTRIBUTION]) {
            status = false;
            error.correction_update(STAT_parsing[STATP_KEY_WORD] , STAT_word[STATW_OBSERVATION_DISTRIBUTION] , line , i + 1);
          }
          break;
        }

        // test separateur

        case POUTPUT : {
          if (token != ":") {
            status = false;
            error.update(STAT_parsing[STATP_SEPARATOR] , line , i + 1);
          }
          break;
        }
        }

        break;
      }

      // test valeur de la probabilite d'observation

      case 3 : {
        if (type == POUTPUT) {
          lstatus = locale.stringToNum(token , &proba);
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

  // test valeurs consecutives

  for (i = 0;i < nb_value;i++) {
    if (!defined_output[i]) {
      status = false;
      error.update(STAT_parsing[STATP_NON_CONSECUTIVE_OUTPUTS]);
      break;
    }
  }

  if (status) {
    process = new Nonparametric_process(nb_state , nb_value , observation_probability);

    // test recouvrement ou pas entre les lois d'observation

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
      process = 0;
    }
  }

  for (i = 0;i < nb_state;i++) {
    delete [] observation_probability[i];
  }
  delete [] observation_probability;

  return process;
}


/*--------------------------------------------------------------*
 *
 *  Analyse du format des lois d'observation pour les
 *  differents processus d'observation.
 *
 *  arguments : reference sur un objet Format_error, stream,
 *              reference sur l'indice de la ligne lue, nombre d'etats,
 *              nombre de processus d'observation.
 *
 *--------------------------------------------------------------*/

Nonparametric_process** observation_parsing(Format_error &error , ifstream &in_file , int &line ,
                                            int nb_state , int &nb_output_process)

{
  RWLocaleSnapshot locale("en");
  RWCString buffer , token;
  size_t position;
  bool status = true , lstatus;
  register int i;
  int index;
  long value;
  Nonparametric_process **process;


  nb_output_process = I_DEFAULT;
  process = 0;

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
    process = new Nonparametric_process*[nb_output_process];
    for (i = 0;i < nb_output_process;i++) {
      process[i] = 0;
    }

    index = 0;

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

        // test mot cle NONPARAMETRIC

        case 3 : {
          if (token != STAT_word[STATW_NONPARAMETRIC]) {
            status = false;
            error.correction_update(STAT_parsing[STATP_KEY_WORD] , STAT_word[STATW_NONPARAMETRIC] , line , i + 1);
          }
          break;
        }
        }

        i++;
      }

      if (i > 0) {
        if ((i != 2) && (i != 4)) {
          status = false;
          error.update(STAT_parsing[STATP_FORMAT] , line);
        }

        process[index - 1] = observation_parsing(error , in_file , line , nb_state , true);
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

      process = 0;
    }
  }

  return process;
}


/*--------------------------------------------------------------*
 *
 *  Analyse du format des lois d'observation pour les
 *  differents processus d'observation.
 *
 *  arguments : reference sur un objet Format_error, stream,
 *              reference sur l'indice de la ligne lue, nombre d'etats,
 *              nombre de processus d'observation.
 *
 *--------------------------------------------------------------*/

Nonparametric_process** old_observation_parsing(Format_error &error , ifstream &in_file , int &line ,
                                                int nb_state , int &nb_output_process)

{
  RWLocaleSnapshot locale("en");
  RWCString buffer , token;
  size_t position;
  bool status = true , lstatus;
  register int i;
  int index;
  long value;
  Nonparametric_process **process;


  nb_output_process = I_DEFAULT;
  process = 0;

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

      // test mot cle OBSERVATION_PROBABILITIES

      if (i == 0) {
        if (token != STAT_word[STATW_OBSERVATION_PROBABILITIES]) {
          status = false;
          error.correction_update(STAT_parsing[STATP_KEY_WORD] , STAT_word[STATW_OBSERVATION_PROBABILITIES] , line);
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

      // test mot cle VARIABLE(S)

      case 1 : {
        if (token != STAT_word[nb_output_process == 1 ? STATW_VARIABLE : STATW_VARIABLES]) {
          status = false;
          error.correction_update(STAT_parsing[STATP_KEY_WORD] ,
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
    process = new Nonparametric_process*[nb_output_process];
    for (i = 0;i < nb_output_process;i++) {
      process[i] = 0;
    }

    index = 0;

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
        switch (i) {

        // test mot cle VARIABLE

        case 0 : {
          if (token == STAT_word[STATW_VARIABLE]) {
            index++;
          }
          else {
            status = false;
            error.correction_update(STAT_parsing[STATP_KEY_WORD] , STAT_word[STATW_VARIABLE] , line , i + 1);
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
        }

        i++;
      }

      if (i > 0) {
        if (i != 2) {
          status = false;
          error.update(STAT_parsing[STATP_FORMAT] , line);
        }

        process[index - 1] = observation_parsing(error , in_file , line , nb_state , true);
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

      process = 0;
    }
  }

  return process;
}


/*--------------------------------------------------------------*
 *
 *  Calcul du nombre de parametres independants d'un processus
 *  d'observation non-parametrique.
 *
 *  argument : probabilite minimum.
 *
 *--------------------------------------------------------------*/

int Nonparametric_process::nb_parameter_computation(double min_probability) const

{
  register int i , j;
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


/*--------------------------------------------------------------*
 *
 *  initialisation des probabilites d'observation.
 *
 *--------------------------------------------------------------*/

void Nonparametric_process::init()

{
  bool *active_output;
  register int i , j , k;
  int nb_noise , *state_nb_value;
  double cumul , noise_proba = NOISE_PROBABILITY , *pmass;


  active_output = new bool[nb_value];
  state_nb_value = new int[nb_state];

  cumul = 0.;
  pmass = observation[0]->mass;
  for (i = 0;i < nb_value;i++) {
    cumul += *pmass++;
  }

  // initialisation uniforme des probabilites d'observation

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

  // perturbation des lois d'observation

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

/*--------------------------------------------------------------*
 *
 *  Affichage des lois d'observation a partir d'un flux de sortie,
 *  des histogrammes correspondant aux observations pour chaque etat,
 *  des flags sur l'affichage exhaustif et dans un fichier
 *
 *--------------------------------------------------------------*/

std::ostream& Nonparametric_process::ascii_print(std::ostream &os, 
						 Histogram **empirical_observation, 
						 bool exhaustive, bool file_flag) const {

  register int i , j;
  int buff , width[2];
  double *pmass;

  if (observation != NULL)
    {
      // affichage des lois d'observation
      for(i= 0; i < nb_state; i++)
      {
         os << "\n" << STAT_word[STATW_STATE] << " " << i << " "
            << STAT_word[STATW_OBSERVATION_DISTRIBUTION] << endl;
         pmass= observation[i]->mass+observation[i]->offset;

         for(j= observation[i]->offset; j < observation[i]->nb_value; j++)
         {
            if (*pmass > 0.)
               os << STAT_word[STATW_OUTPUT] << " " << j << " : " << *pmass << endl;
            pmass++;
         }

	 // affichage des histogrammes correspondants
         if ((empirical_observation != NULL) && (exhaustive))
         {
            os << "\n";
            if (file_flag)
               os << "# ";

            os << "   | " << STAT_label[STATL_STATE] << " " << i << " "
               << STAT_label[STATL_OBSERVATION] << " " << STAT_label[STATL_HISTOGRAM]
               << " | " << STAT_label[STATL_STATE] << " " << i << " "
               << STAT_label[STATL_OBSERVATION] << " " << STAT_label[STATL_DISTRIBUTION] << endl;

            observation[i]->ascii_print(os, file_flag, false, false, empirical_observation[i]);
         }
      }

      // calcul de la largeur de colonne

      width[0]= column_width(nb_state-1);

      width[1]= 0;
      for(i= 0; i < nb_state; i++)
      {
         buff= column_width(observation[i]->nb_value, observation[i]->mass);
         if (buff > width[1])
            width[1] = buff;
      }
      width[1]+= ASCII_SPACE;

      // affichage de la matrice des probabilites d'observation
      os << "\n";
      if (file_flag)
         os << "# ";

      os << STAT_label[STATL_OBSERVATION_PROBABILITIY_MATRIX] << endl;

      os << "\n";
      if (file_flag)
         os << "# ";

      os << setw(width[0]+width[1]) << 0;
      for(i= 1; i < nb_value; i++)
         os << setw(width[1]) << i;

      for(i= 0; i < nb_state; i++)
      {
         os << "\n";
         if (file_flag)
            os << "# ";

         os << setw(width[0]) << i ;
         for(j= 0; j < nb_value; j++)
            os << setw(width[1]) << observation[i]->mass[j];
      }
      os << endl;
    }
  
  return os;

}
/*--------------------------------------------------------------*
 *
 *  Affichage des lois d'observation au format tableur 
 *  a partir d'un flux de sortie et des histogrammes correspondant 
 *  aux observations pour chaque etat.
 *
 *--------------------------------------------------------------*/

std::ostream& Nonparametric_process::spreadsheet_print(std::ostream &os, 
						       Histogram **empirical_observation) const {

  register int i , j;
  int buff , width[2];
  double *pmass;

  if (observation != NULL)
    {
      // affichage des lois d'observation
      for(i= 0; i < nb_state; i++)
      {
         os << "\n" << STAT_word[STATW_STATE] << "\t" << i << "\t"
            << STAT_word[STATW_OBSERVATION_DISTRIBUTION] << endl;
         pmass= observation[i]->mass+observation[i]->offset;

         for(j= observation[i]->offset; j < observation[i]->nb_value; j++)
         {
            if (*pmass > 0.)
               os << STAT_word[STATW_OUTPUT] << "\t" << j << "\t" << *pmass << endl;
            pmass++;
         }

	 // affichage des histogrammes correspondants
         if ((empirical_observation != NULL))
         {
            os << "\n";

            os << "   | " << STAT_label[STATL_STATE] << "\t" << i << "\t"
               << STAT_label[STATL_OBSERVATION] << "\t" << STAT_label[STATL_HISTOGRAM]
               << " | " << STAT_label[STATL_STATE] << "\t" << i << "\t"
               << STAT_label[STATL_OBSERVATION] << "\t" << STAT_label[STATL_DISTRIBUTION] << endl;

            observation[i]->spreadsheet_print(os, false, false, false, empirical_observation[i]);
         }
      }


      // affichage de la matrice des probabilites d'observation
      os << "\n";

      os << STAT_label[STATL_OBSERVATION_PROBABILITIY_MATRIX] << endl;

      os << "\n";

      for(i= 0; i < nb_state; i++) {
        for(j= 0; j < observation[i]->nb_value; j++) 
	  os << observation[i]->mass[j] << "\t";
	os << "\n";
      }

    }
  
  return os;

}

/*--------------------------------------------------------------*
 *
 *  Sortie gnuplot d'un processus non parametrique a partir du
 *  prefixe de nom de fichier, le titre du graphe, la valeur 
 *  de l'etat, et les histogrammes correspondants aux lois
 *  d'observation
 *
 *--------------------------------------------------------------*/

bool Nonparametric_process::plot_print(const char *prefix, const char *title, 
				       int process, Histogram **empirical_observation) const {
				       
   bool status= false, start;
   register int val, i, j= 0, k= 0;
   int nb_histo, nb_dist, histo_index,
       dist_index, plot_index, *dist_nb_value, *index_dist; 
   double *scale;
   const Distribution **pdist;
   const Histogram **phisto;
   ostringstream data_file_name[2];

   // data file printing

   data_file_name[0] << prefix << process << 0 << ".dat";


   // name of the file containing the data for all the graphs
   data_file_name[1] << prefix << process << 1 << ".dat";

   // list of the successive distributions to be printed
   pdist= new const Distribution*[1 * nb_value + nb_state];
   dist_nb_value= new int[1 * nb_value + nb_state];
   // scales for printing distributions and histograms together
   scale= new double[1 * nb_value + nb_state];
   // list of the successive histograms to be printed
   phisto= new const Histogram*[1 + 2 * nb_value + nb_state];
   // correspondance between the histograms and the distributions ?
   index_dist= new int[1 + 2 * nb_value + nb_state];

   nb_histo= 0;
   nb_dist= 0;

   if (observation != NULL)
   {
      for(val= 0; val < nb_state; val++)
      {
         pdist[nb_dist]= observation[val];
         dist_nb_value[nb_dist]= observation[val]->nb_value;

         if (empirical_observation != NULL)
         {
           phisto[nb_histo]= empirical_observation[val];
           index_dist[nb_histo]= nb_dist;
           scale[nb_dist++]= phisto[nb_histo++]->nb_element;
         }
         else
            scale[nb_dist++]= 1.;
      }
   }

   // create the file, named prefix(variable+1)1.dat, containing
   // the values for all the gnuplot graphs

   status= ::plot_print((data_file_name[1].str()).c_str(), nb_dist, pdist,
			scale, dist_nb_value, nb_histo, phisto, index_dist);

   // write the command and printing files
   if (status) {
      // i == 0 -> output on terminal (.plot file)
      // i == 1 -> output into a postscript file  (.plot file)

     histo_index= 0; // or 1?
     dist_index= 0;

      // observation
      if (observation != NULL)
      {
         for(i= 0; i < 2; i++)
         {
            ostringstream file_name[2];

            switch (i)
            {
              case 0 :
                 file_name[0] << prefix << process << 0 << ".plot";
                 break;
              case 1 :
                 file_name[0] << prefix << process << 0 << ".print";
                 break;
            }

            ofstream out_file((file_name[0].str()).c_str());

            if (i == 1)
            {
               out_file << "set terminal postscript" << endl;
               file_name[1] << label(prefix) << process << 0 << ".ps";
               out_file << "set output \"" << file_name[1].str() << "\"\n\n";
            }

            out_file << "set border 15 lw 0\n" << "set tics out\n" << "set xtics nomirror\n"
                     << "set title";

            if (title != NULL)
               out_file << " \"" << title << " - " << STAT_label[STATL_OUTPUT_PROCESS]
                        << "\t" << process << "\"";
            out_file << "\n\n";

            j= histo_index;
            k= dist_index;

            for(val= 0; val < nb_state; val++)
            {
               if (dist_nb_value[k] - 1 < TIC_THRESHOLD)
                  out_file << "set xtics 0,1" << endl;

               if (empirical_observation != NULL)
               {
                  out_file << "plot [0:" << dist_nb_value[k]-1 << "] [0:"
                           << (int)(MAX(phisto[j]->max, pdist[k]->max * scale[k])*YSCALE)+1
                           << "] \"" << label((data_file_name[1].str()).c_str()) << "\" using " << j+1
                           << " title \"" << STAT_label[STATL_STATE] << "\t" << val << "\t"
                           << STAT_label[STATL_OBSERVATION] << "\t" << STAT_label[STATL_HISTOGRAM]
                           << "\" with impulses,\\" << endl;
                  out_file << "\"" << label((data_file_name[1].str()).c_str()) << "\" using " << nb_histo+k+1
                           << " title \"" << STAT_label[STATL_STATE] << "\t" << val << "\t"
                           << STAT_label[STATL_OBSERVATION] << "\t" << STAT_label[STATL_DISTRIBUTION]
                           << "\" with linespoints" << endl;
                  j++;
               }
               else
               {
                  out_file << "plot [0:" << dist_nb_value[k]-1 << "] [0:"
                           << MIN(pdist[k]->max*YSCALE, 1.) << "] \""
                           << label((data_file_name[1].str()).c_str()) << "\" using " << nb_histo+k+1
                           << " title \"" << STAT_label[STATL_STATE] << "\t" << val << "\t"
                           << STAT_label[STATL_OBSERVATION] << "\t" << STAT_label[STATL_DISTRIBUTION]
                           << "\" with linespoints" << endl;
               }

               if (dist_nb_value[k]-1 < TIC_THRESHOLD)
                  out_file << "set xtics autofreq" << endl;
               k++;

               if ((i == 0) && (val < nb_state - 1))
                  out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
               out_file << endl;
            }

            if (i == 1)
               out_file << "\nset terminal x11" << endl;

            out_file << "\npause 0 \"" << STAT_label[STATL_END] << "\"" << endl;
         }
      }

      delete [] pdist;
      delete [] dist_nb_value;
      delete [] scale;
      delete [] phisto;
      delete [] index_dist;
   }

   return status;

}
