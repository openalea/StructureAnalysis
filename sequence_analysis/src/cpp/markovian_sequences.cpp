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
#include <iomanip>
// #include <rw/vstream.h>
// #include <rw/rwfile.h>
#include "stat_tool/stat_tools.h"
#include "stat_tool/distribution.h"
#include "stat_tool/curves.h"
#include "stat_tool/markovian.h"
#include "stat_tool/stat_label.h"
#include "sequences.h"
#include "sequence_label.h"
#include "tool/config.h"

using namespace std;


extern int* select_variable(int nb_variable , int selected_nb_variable ,
                            int *selected_variable , bool keep);

extern int column_width(int value);
extern char* label(const char *file_name);



/*--------------------------------------------------------------*
 *
 *  Constructeur par defaut de la classe Markovian_sequences.
 *
 *--------------------------------------------------------------*/

Markovian_sequences::Markovian_sequences()

{
  self_transition = 0;
  observation = 0;
  characteristics = 0;
}


/*--------------------------------------------------------------*
 *
 *  Initialisation par defaut des champs de la classe Markovian_sequences.
 *
 *--------------------------------------------------------------*/

void Markovian_sequences::init()

{
  register int i;


  self_transition = 0;
  observation = 0;

  characteristics = new Sequence_characteristics*[nb_variable];
  for (i = 0;i < nb_variable;i++) {
    characteristics[i] = 0;
  }
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Markovian_sequences.
 *
 *  arguments : nombre de variables, nombre de sequences,
 *              longueurs des sequences, sequences.
 *
 *--------------------------------------------------------------*/

Markovian_sequences::Markovian_sequences(int inb_variable , int inb_sequence ,
                                         int *ilength , int ***isequence)
:Sequences(inb_variable , inb_sequence , ilength , isequence)

{
  init();
  build_characteristic();
}


/*--------------------------------------------------------------*
 *
 *  Construction d'un objet Markovian_sequences a partir d'un objet Sequences.
 *
 *  argument : reference sur un objet Sequences.
 *
 *--------------------------------------------------------------*/

Markovian_sequences::Markovian_sequences(const Sequences &seq)
:Sequences(seq)

{
  init();
  build_characteristic();
}


/*--------------------------------------------------------------*
 *
 *  Copie d'un objet Markovian_sequences.
 *
 *  arguments : reference sur un objet Markovian_sequences, inversion /
 *              ajout/suppression des histogrammes de temps de sejour initial.
 *
 *--------------------------------------------------------------*/

void Markovian_sequences::copy(const Markovian_sequences &seq , int param)

{
  bool initial_run_flag;
  register int i , j;


  if ((seq.self_transition) && (param != REVERSE)) {
    self_transition = new Self_transition*[marginal[0]->nb_value];
    for (i = 0;i < marginal[0]->nb_value;i++) {
      if (seq.self_transition[i]) {
        self_transition[i] = new Self_transition(*(seq.self_transition[i]));
      }
      else {
        self_transition[i] = 0;
      }
    }
  }

  else {
    self_transition = 0;
  }

  if (seq.observation) {
    observation = new Histogram**[nb_variable];
    observation[0] = 0;

    for (i = 1;i < nb_variable;i++) {
      observation[i] = new Histogram*[marginal[0]->nb_value];
      for (j = 0;j < marginal[0]->nb_value;j++) {
        observation[i][j] = new Histogram(*(seq.observation[i][j]));
      }
    }
  }

  else {
    observation = 0;
  }

  characteristics = new Sequence_characteristics*[nb_variable];

  for (i = 0;i < nb_variable;i++) {
    if (seq.characteristics[i]) {
      if ((param == ADD_INITIAL_RUN) || (param == REMOVE_INITIAL_RUN)) {
        switch (param) {
        case ADD_INITIAL_RUN :
          initial_run_flag = true;
          break;
        case REMOVE_INITIAL_RUN :
          initial_run_flag = false;
          break;
        }

        characteristics[i] = new Sequence_characteristics(*(seq.characteristics[i]) , initial_run_flag);

        if (((seq.characteristics[i]->initial_run) && (!initial_run_flag)) ||
            ((!(seq.characteristics[i]->initial_run)) && (initial_run_flag))) {
           build_sojourn_time_histogram(i , initial_run_flag);
        }
      }

      else if (param == REVERSE) {
        characteristics[i] = new Sequence_characteristics(*(seq.characteristics[i]), 'r');

        build_index_value(i);
        build_first_occurrence_histogram(i);

        if (!(seq.characteristics[i]->initial_run)) {
          build_sojourn_time_histogram(i);
        }
      }

      else {
        characteristics[i] = new Sequence_characteristics(*(seq.characteristics[i]));
      }
    }

    else {
      characteristics[i] = 0;
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Copie d'un objet Markovian_sequences avec ajout d'une variable.
 *
 *  arguments : reference sur un objet Markovian_sequences, indice de la variable,
 *              ajout/suppression des histogrammes de temps de sejour initial.
 *
 *--------------------------------------------------------------*/

void Markovian_sequences::add_variable(const Markovian_sequences &seq ,
                                       int variable , int param)

{
  bool initial_run_flag;
  register int i , j;


  self_transition = 0;
  observation = 0;

  characteristics = new Sequence_characteristics*[nb_variable];

  i = 0;
  for (j = 0;j < nb_variable;j++) {
    if (j != variable) {
      if (seq.characteristics[i]) {
        if ((param == ADD_INITIAL_RUN) || (param == REMOVE_INITIAL_RUN)) {
          switch (param) {
          case ADD_INITIAL_RUN :
            initial_run_flag = true;
            break;
          case REMOVE_INITIAL_RUN :
            initial_run_flag = false;
            break;
          }

          characteristics[j] = new Sequence_characteristics(*(seq.characteristics[i]) , initial_run_flag);

          if (((seq.characteristics[i]->initial_run) && (!initial_run_flag)) ||
              ((!(seq.characteristics[i]->initial_run)) && (initial_run_flag))) {
             build_sojourn_time_histogram(j , initial_run_flag);
          }
        }

        else {
          characteristics[j] = new Sequence_characteristics(*(seq.characteristics[i]));
        }
      }

      else {
        characteristics[j] = 0;
      }

      i++;
    }

    else {
      characteristics[j] = 0;
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Constructeur par copie de la classe Markovian_sequences.
 *
 *  arguments : reference sur un objet Markovian_sequences, type de transformation
 *              ('c' : copie, 'a' : addition d'une variable),
 *              inversion / ajout/suppression des histogrammes de temps de sejour initial /
 *              indice de la variable ('a').
 *
 *--------------------------------------------------------------*/

Markovian_sequences::Markovian_sequences(const Markovian_sequences &seq , char transform ,
                                         int param1 , int param2)

{
  switch (transform) {
  case 'c' :
    Sequences::copy(seq , (param1 == REVERSE ? true : false));
    copy(seq , param1);
    break;
  case 'a' :
    Sequences::add_variable(seq , param1);
    add_variable(seq , param1 , param2);
    break;
  default :
    Sequences::copy(seq);
    copy(seq);
    break;
  }
}


/*--------------------------------------------------------------*
 *
 *  Destructeur des champs de la classe Markovian_sequences.
 *
 *--------------------------------------------------------------*/

void Markovian_sequences::remove()

{
  register int i , j;


  if (self_transition) {
    for (i = 0;i < marginal[0]->nb_value;i++) {
      delete self_transition[i];
    }
    delete [] self_transition;
  }

  if (observation) {
    for (i = 1;i < nb_variable;i++) {
      for (j = 0;j < marginal[0]->nb_value;j++) {
        delete observation[i][j];
      }
      delete [] observation [i];
    }
    delete [] observation;
  }

  if (characteristics) {
    for (i = 0;i < nb_variable;i++) {
      delete characteristics[i];
    }
    delete [] characteristics;
  }  
}


/*--------------------------------------------------------------*
 *
 *  Destructeur de la classe Markovian_sequences.
 *
 *--------------------------------------------------------------*/

Markovian_sequences::~Markovian_sequences()

{
  remove();
}


/*--------------------------------------------------------------*
 *
 *  Operateur d'assignement de la classe Markovian_sequences.
 *
 *  argument : reference sur un objet Markovian_sequences.
 *
 *--------------------------------------------------------------*/

Markovian_sequences& Markovian_sequences::operator=(const Markovian_sequences &seq)

{
  if (&seq != this) {
    remove();
    Sequences::remove();

    Sequences::copy(seq);
    copy(seq);
  }

  return *this;
}


/*--------------------------------------------------------------*
 *
 *  Initialisation de la 1ere variable.
 *
 *  argument : type de la 1ere variable (STATE / INT_VALUE).
 *
 *--------------------------------------------------------------*/

void Markovian_sequences::state_variable_init(int itype)

{
  register int i , j;


  if (itype != type[0]) {
    if (type[0] == STATE) {
      if (self_transition) {
        for (i = 0;i < marginal[0]->nb_value;i++) {
          delete self_transition[i];
        }
        delete [] self_transition;

        self_transition = 0;
      }

      if (observation) {
        for (i = 1;i < nb_variable;i++) {
          for (j = 0;j < marginal[0]->nb_value;j++) {
            delete observation[i][j];
          }
          delete [] observation[i];
        }
        delete [] observation;

        observation = 0;
      }
    }

    type[0] = itype;
  }

  for (i = 1;i < nb_variable;i++) {
    type[i] = INT_VALUE;
  }
}


/*--------------------------------------------------------------*
 *
 *  Extraction d'un histogramme.
 *
 *  arguments : reference sur un objet Format_error, type d'histogramme,
 *              variable, valeur.
 *
 *--------------------------------------------------------------*/

Distribution_data* Markovian_sequences::extract(Format_error &error , int type ,
                                                int variable , int value) const

{
  bool status = true;
  Histogram *phisto;
  Distribution_data *histo;


  histo = 0;
  error.init();

  if ((variable < 1) || (variable > nb_variable)) {
    status = false;
    error.update(STAT_error[STATR_VARIABLE_INDEX]);
  }

  if (status) {
    variable--;

    if (!characteristics[variable]) {
      status = false;
      error.update(SEQ_error[SEQR_CHARACTERISTICS_NOT_COMPUTED]);
    }

    else if ((value < 0) || (value >= marginal[variable]->nb_value)) {
      status = false;
      ostringstream error_message;
      error_message << STAT_label[STATL_VALUE] << " " << value << " "
                    << SEQ_error[SEQR_NOT_PRESENT];
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
          phisto = 0;
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

  if (status) {
    histo = new Distribution_data(*phisto);
  }

  return histo;
}


/*--------------------------------------------------------------*
 *
 *  Fusion d'objets Markovian_sequences.
 *
 *  arguments : reference sur un objet Format_error, nombre d'objets Markovian_sequences,
 *              pointeurs sur les objets Markovian_sequences.
 *
 *--------------------------------------------------------------*/

Markovian_sequences* Markovian_sequences::merge(Format_error &error , int nb_sample ,
                                                const Markovian_sequences **iseq) const

{
  bool status = true;
  register int i , j , k , m , n;
  int inb_sequence , nb_histo , *itype , *ilength , *psequence , *csequence;
  const Histogram **phisto;
  Markovian_sequences *seq;
  const Markovian_sequences **pseq;


  seq = 0;
  error.init();

  for (i = 0;i < nb_sample;i++) {
    if (iseq[i]->nb_variable != nb_variable) {
      status = false;
      ostringstream error_message;
      error_message << STAT_label[STATL_SAMPLE] << " " << i + 2 << ": "
                    << STAT_error[STATR_NB_VARIABLE];
      error.correction_update((error_message.str()).c_str() , nb_variable);
    }
  }

  if (status) {
    nb_sample++;
    pseq = new const Markovian_sequences*[nb_sample];

    pseq[0] = this;
    for (i = 1;i < nb_sample;i++) {
      pseq[i] = iseq[i - 1];
    }

    // calcul du nombre et de la longueur des sequences

    inb_sequence = 0;
    for (i = 0;i < nb_sample;i++) {
      inb_sequence += pseq[i]->nb_sequence;
    }

    itype = new int[nb_variable];

    for (i = 0;i < nb_variable;i++) {
      for (j = 1;j < nb_sample;j++) {
        if (pseq[j]->type[i] != pseq[0]->type[i]) {
          break;
        }
      }
      if (j < nb_sample) {
        itype[i] = INT_VALUE;
      }
      else {
        itype[i] = pseq[0]->type[i];
      }
    }

    ilength = new int[inb_sequence];

    i = 0;
    for (j = 0;j < nb_sample;j++) {
      for (k = 0;k < pseq[j]->nb_sequence;k++) {
        ilength[i++] = pseq[j]->length[k];
      }
    }

    seq = new Markovian_sequences(nb_variable , itype , inb_sequence ,
                                  0 , ilength , false);
    delete [] itype;
    delete [] ilength;

    // copie des sequences

    i = 0;
    for (j = 0;j < nb_sample;j++) {
      for (k = 0;k < pseq[j]->nb_sequence;k++) {
        for (m = 0;m < seq->nb_variable;m++) {
          psequence = seq->sequence[i][m];
          csequence = pseq[j]->sequence[k][m];
          for (n = 0;n < pseq[j]->length[k];n++) {
            *psequence++ = *csequence++;
          }
        }
        i++;
      }
    }

    phisto = new const Histogram*[nb_sample];

    for (i = 0;i < seq->nb_variable;i++) {
      seq->min_value[i] = pseq[0]->min_value[i];
      seq->max_value[i] = pseq[0]->max_value[i];
      for (j = 1;j < nb_sample;j++) {
        if (pseq[j]->min_value[i] < seq->min_value[i]) {
          seq->min_value[i] = pseq[j]->min_value[i];
        }
        if (pseq[j]->max_value[i] > seq->max_value[i]) {
          seq->max_value[i] = pseq[j]->max_value[i];
        }
      }

      for (j = 0;j < nb_sample;j++) {
        phisto[j] = pseq[j]->marginal[i];
      }
      seq->marginal[i] = new Histogram(nb_sample , phisto);

      for (j = 0;j < nb_sample;j++) {
        if (!(pseq[j]->characteristics[i])) {
          break;
        }
      }

      if (j == nb_sample) {
        seq->characteristics[i] = new Sequence_characteristics();

        seq->characteristics[i]->nb_value = seq->marginal[i]->nb_value;

        seq->build_index_value(i);

        seq->characteristics[i]->first_occurrence = new Histogram*[seq->marginal[i]->nb_value];
        seq->characteristics[i]->recurrence_time = new Histogram*[seq->marginal[i]->nb_value];
        seq->characteristics[i]->nb_run = new Histogram*[seq->marginal[i]->nb_value];
        seq->characteristics[i]->nb_occurrence = new Histogram*[seq->marginal[i]->nb_value];

        for (j = 0;j < seq->marginal[i]->nb_value;j++) {
          nb_histo = 0;
          for (k = 0;k < nb_sample;k++) {
            if (j < pseq[k]->marginal[i]->nb_value) {
              phisto[nb_histo++] = pseq[k]->characteristics[i]->first_occurrence[j];
            }
          }
          seq->characteristics[i]->first_occurrence[j] = new Histogram(nb_histo , phisto);

          nb_histo = 0;
          for (k = 0;k < nb_sample;k++) {
            if (j < pseq[k]->marginal[i]->nb_value) {
              phisto[nb_histo++] = pseq[k]->characteristics[i]->recurrence_time[j];
            }
          }
          seq->characteristics[i]->recurrence_time[j] = new Histogram(nb_histo , phisto);

          nb_histo = 0;
          for (k = 0;k < nb_sample;k++) {
            if (j < pseq[k]->marginal[i]->nb_value) {
              phisto[nb_histo++] = pseq[k]->characteristics[i]->nb_run[j];
            }
          }
          seq->characteristics[i]->nb_run[j] = new Histogram(nb_histo , phisto);

          nb_histo = 0;
          for (k = 0;k < nb_sample;k++) {
            if (j < pseq[k]->marginal[i]->nb_value) {
              phisto[nb_histo++] = pseq[k]->characteristics[i]->nb_occurrence[j];
            }
          }
          seq->characteristics[i]->nb_occurrence[j] = new Histogram(nb_histo , phisto);
        }

        for (j = 1;j < nb_sample;j++) {
          if (((pseq[0]->characteristics[i]->initial_run) && (!(pseq[j]->characteristics[i]->initial_run))) ||
              ((!(pseq[0]->characteristics[i]->initial_run)) && (pseq[j]->characteristics[i]->initial_run))) {
            break;
          }
        }

        if (j == nb_sample) {
          seq->characteristics[i]->sojourn_time = new Histogram*[seq->marginal[i]->nb_value];
          if (pseq[0]->characteristics[i]->initial_run) {
            seq->characteristics[i]->initial_run = new Histogram*[seq->marginal[i]->nb_value];
          }
          seq->characteristics[i]->final_run = new Histogram*[seq->marginal[i]->nb_value];

          for (j = 0;j < seq->marginal[i]->nb_value;j++) {
            nb_histo = 0;
            for (k = 0;k < nb_sample;k++) {
              if (j < pseq[k]->marginal[i]->nb_value) {
                phisto[nb_histo++] = pseq[k]->characteristics[i]->sojourn_time[j];
              }
            }
            seq->characteristics[i]->sojourn_time[j] = new Histogram(nb_histo , phisto);

            if (pseq[0]->characteristics[i]->initial_run) {
              nb_histo = 0;
              for (k = 0;k < nb_sample;k++) {
                if (j < pseq[k]->marginal[i]->nb_value) {
                  phisto[nb_histo++] = pseq[k]->characteristics[i]->initial_run[j];
                }
              }
              seq->characteristics[i]->initial_run[j] = new Histogram(nb_histo , phisto);
            }

            nb_histo = 0;
            for (k = 0;k < nb_sample;k++) {
              if (j < pseq[k]->marginal[i]->nb_value) {
                phisto[nb_histo++] = pseq[k]->characteristics[i]->final_run[j];
              }
            }
            seq->characteristics[i]->final_run[j] = new Histogram(nb_histo , phisto);
          }
        }

        else {
          seq->build_sojourn_time_histogram(i , (characteristics[i]->initial_run ? true : false));
        }
      }

      else {
        seq->build_characteristic(i , true , (((characteristics[i]) && (characteristics[i]->initial_run)) ? true : false));
      }
    }

    delete [] pseq;
    delete [] phisto;
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Regroupement des valeurs d'une variable donnee.
 *
 *  arguments : reference sur un objet Format_error, indice de la variable,
 *              pas de regroupement, flag pour ajouter une variable.
 *
 *--------------------------------------------------------------*/

Markovian_sequences* Markovian_sequences::cluster(Format_error &error , int ivariable ,
                                                  int step , bool add_flag) const

{
  bool status = true;
  register int i;
  int variable , offset , *itype;
  Markovian_sequences *seq;


  seq = 0;
  error.init();

  if ((ivariable < 1) || (ivariable > nb_variable)) {
    status = false;
    error.update(STAT_error[STATR_VARIABLE_INDEX]);
  }
  if (step < 1) {
    status = false;
    error.update(STAT_error[STATR_CLUSTERING_STEP]);
  }

  if (status) {
    ivariable--;

    switch (add_flag) {
    case false :
      variable = ivariable;
      offset = 0;
      break;
    case true :
      variable = 0;
      offset = 1;
      break;
    }

    itype = new int[nb_variable + offset];
    for (i = 0;i < nb_variable;i++) {
      itype[i + offset] = type[i];
    }
    itype[variable] = INT_VALUE;

    seq = new Markovian_sequences(nb_variable + offset , itype , nb_sequence ,
                                  identifier , length , false);
    delete [] itype;

    seq->Sequences::cluster(*this , ivariable , step , add_flag);

    for (i = 0;i < seq->nb_variable;i++) {
      if (i == variable) {
        seq->build_characteristic(i , true , (((characteristics[ivariable]) && (characteristics[ivariable]->initial_run)) ? true : false));
      }
      else if (characteristics[i - offset]) {
        seq->characteristics[i] = new Sequence_characteristics(*(characteristics[i - offset]));
      }
    }
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Transcodage des symboles d'une variable donnee.
 *
 *  arguments : reference sur un objet Format_error, indice de la variable,
 *              table de transcodage des symboles, flag pour ajouter une variable.
 *
 *--------------------------------------------------------------*/

Markovian_sequences* Markovian_sequences::transcode(Format_error &error , int ivariable ,
                                                    int *symbol , bool add_flag) const

{
  bool status = true , *presence;
  register int i;
  int variable , offset , min_symbol , max_symbol , *itype;
  Markovian_sequences *seq;


  seq = 0;
  error.init();

  if ((ivariable < 1) || (ivariable > nb_variable)) {
    status = false;
    error.update(STAT_error[STATR_VARIABLE_INDEX]);
  }

  else {
    ivariable--;

    min_symbol = marginal[ivariable]->nb_value;
    max_symbol = 0;

    for (i = 0;i < marginal[ivariable]->nb_value - marginal[ivariable]->offset;i++) {
      if ((symbol[i] < 0) || (symbol[i] >= (add_flag ? marginal[ivariable]->nb_value - 1 : marginal[ivariable]->nb_value))) {
        status = false;
        ostringstream error_message;
        error_message << STAT_label[STATL_SYMBOL] << " " << symbol[i] << " "
                      << STAT_error[STATR_NOT_ALLOWED];
        error.update((error_message.str()).c_str());
      }

      else {
        if (symbol[i] < min_symbol) {
          min_symbol = symbol[i];
        }
        if (symbol[i] > max_symbol) {
          max_symbol = symbol[i];
        }
      }
    }

    if ((min_symbol != 0) || (max_symbol == 0)) {
      status = false;
      error.update(STAT_error[STATR_NB_SYMBOL]);
    }

    if (status) {
      presence = new bool[max_symbol + 1];
      for (i = 0;i <= max_symbol;i++) {
        presence[i] = false;
      }

      for (i = 0;i < marginal[ivariable]->nb_value - marginal[ivariable]->offset;i++) {
        presence[symbol[i]] = true;
      }

      for (i = 0;i <= max_symbol;i++) {
        if (!presence[i]) {
          status = false;
          ostringstream error_message;
          error_message << STAT_error[STATR_MISSING_SYMBOL] << " " << i;
          error.update((error_message.str()).c_str());
        }
      }

      delete [] presence;
    }

    if (status) {
      switch (add_flag) {
      case false :
        variable = ivariable;
        offset = 0;
        break;
      case true :
        variable = 0;
        offset = 1;
        break;
      }

      itype = new int[nb_variable + offset];
      for (i = 0;i < nb_variable;i++) {
        itype[i + offset] = type[i];
      }
      itype[variable] = INT_VALUE;

      seq = new Markovian_sequences(nb_variable + offset , itype , nb_sequence ,
                                    identifier , length , false);
      delete [] itype;

      seq->Sequences::transcode(*this , ivariable , 0 , max_symbol , symbol , add_flag);

      for (i = 0;i < seq->nb_variable;i++) {
        if (i == variable) {
          seq->build_characteristic(i , true , (((characteristics[ivariable]) && (characteristics[ivariable]->initial_run)) ? true : false));
        }
        else if (characteristics[i - offset]) {
          seq->characteristics[i] = new Sequence_characteristics(*(characteristics[i - offset]));
        }
      }
    }
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Transcodage des symboles d'une variable donnee.
 *
 *  arguments : references sur un objet Format_error et
 *              sur un objet Nonparametric_sequence_process.
 *
 *--------------------------------------------------------------*/

Markovian_sequences* Markovian_sequences::transcode(Format_error &error ,
                                                    const Nonparametric_sequence_process *process) const

{
  register int i , j;
  int *symbol;
  Markovian_sequences *seq;


  symbol = new int[process->nb_value];
  for (i = 0;i < process->nb_state;i++) {
    for (j = 0;j < process->nb_value;j++) {
      if (process->observation[i]->mass[j] > 0.) {
        symbol[j] = i;
      }
    }
  }

  seq = transcode(error , 1 , symbol , true);
  delete [] symbol;

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Regroupement des symboles d'une variable donnee.
 *
 *  arguments : reference sur un objet Format_error, indice de la variable,
 *              nombres de classes, bornes pour regrouper les symboles,
 *              flag pour ajouter une variable.
 *
 *--------------------------------------------------------------*/

Markovian_sequences* Markovian_sequences::cluster(Format_error &error , int ivariable , int nb_class ,
                                                  int *ilimit , bool add_flag) const

{
  bool status = true;
  register int i , j , k;
  int variable , offset , *symbol , *limit , *itype;
  Markovian_sequences *seq;


  seq = 0;
  error.init();

  if ((ivariable < 1) || (ivariable > nb_variable)) {
    status = false;
    error.update(STAT_error[STATR_VARIABLE_INDEX]);
  }

  else {
    ivariable--;

    if ((nb_class < 2) || (nb_class >= marginal[ivariable]->nb_value)) {
      status = false;
      error.update(STAT_error[STATR_NB_CLASS]);
    }
  }

  if (status) {
    limit = new int[nb_class + 1];
    limit[0] = marginal[ivariable]->offset;
    for (i = 1;i < nb_class;i++) {
      limit[i] = ilimit[i - 1];
    }
    limit[nb_class] = marginal[ivariable]->nb_value;

    for (i = 1;i <= nb_class;i++) {
      if (limit[i] <= limit[i - 1]) {
        status = false;
        error.update(STAT_error[STATR_CLUSTER_LIMIT]);
      }
    }

    if (status) {
      symbol = new int[marginal[ivariable]->nb_value];

      i = 0;
      for (j = 0;j < nb_class;j++) {
        for (k = limit[j];k < limit[j + 1];k++) {
          symbol[i++] = j;
        }
      }

      switch (add_flag) {
      case false :
        variable = ivariable;
        offset = 0;
        break;
      case true :
        variable = 0;
        offset = 1;
        break;
      }

      itype = new int[nb_variable + offset];
      for (i = 0;i < nb_variable;i++) {
        itype[i + offset] = type[i];
      }
      itype[variable] = INT_VALUE;

      seq = new Markovian_sequences(nb_variable + offset , itype , nb_sequence ,
                                    identifier , length , false);
      delete [] itype;

      seq->Sequences::transcode(*this , ivariable , 0 , nb_class - 1 , symbol , add_flag);

      delete [] symbol;

      for (i = 0;i < seq->nb_variable;i++) {
        if (i == variable) {
          seq->build_characteristic(i , true , (((characteristics[ivariable]) && (characteristics[ivariable]->initial_run)) ? true : false));
        }
        else if (characteristics[i - offset]) {
          seq->characteristics[i] = new Sequence_characteristics(*(characteristics[i - offset]));
        }
      }
    }

    delete [] limit;
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Selection de variables.
 *
 *  arguments : reference sur un objet Format_error, nombre de variables,
 *              indices des variables, flag pour conserver ou rejeter
 *              les variables selectionnees.
 *
 *--------------------------------------------------------------*/

Markovian_sequences* Markovian_sequences::select_variable(Format_error &error , int inb_variable ,
                                                          int *ivariable , bool keep) const

{
  bool status = true , *selected_variable;
  register int i;
  int *variable , *itype;
  Markovian_sequences *seq;


  seq = 0;
  error.init();

  if ((inb_variable < 1) || (inb_variable > (keep ? nb_variable : nb_variable - 1))) {
    status = false;
    error.update(STAT_error[STATR_NB_SELECTED_VARIABLE]);
  }

  if (status) {
    selected_variable = new bool[nb_variable + 1];
    for (i = 1;i <= nb_variable;i++) {
      selected_variable[i] = false;
    }

    for (i = 0;i < inb_variable;i++) {
      if ((ivariable[i] < 1) || (ivariable[i] > nb_variable)) {
        status = false;
        ostringstream error_message;
        error_message << ivariable[i] << ": " << STAT_error[STATR_VARIABLE_INDEX];
        error.update((error_message.str()).c_str());
      }

      else if (selected_variable[ivariable[i]]) {
        status = false;
        ostringstream error_message;
        error_message << STAT_label[STATL_VARIABLE] << " " << ivariable[i] << " "
                      << STAT_error[STATR_ALREADY_SELECTED];
        error.update((error_message.str()).c_str());
      }
      else {
        selected_variable[ivariable[i]] = true;
      }
    }

    delete [] selected_variable;
  }

  if (status) {
    variable = ::select_variable(nb_variable , inb_variable , ivariable , keep);

    itype = new int[keep ? inb_variable : nb_variable - inb_variable];
    for (i = 0;i < (keep ? inb_variable : nb_variable - inb_variable);i++) {
      itype[i] = type[variable[i]];
    }

    seq = new Markovian_sequences((keep ? inb_variable : nb_variable - inb_variable) ,
                                  itype , nb_sequence , identifier , length , false);

    seq->Sequences::select_variable(*this , variable);

    for (i = 0;i < seq->nb_variable;i++) {
      if (characteristics[variable[i]]) {
        seq->characteristics[i] = new Sequence_characteristics(*(characteristics[variable[i]]));
      }
    }

    delete [] variable;
    delete [] itype;
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Suppression de la 1ere variable.
 *
 *--------------------------------------------------------------*/

Markovian_sequences* Markovian_sequences::remove_variable_1() const

{
  register int i;
  int *variable , *itype;
  Markovian_sequences *seq;


  variable = new int[nb_variable - 1];
  itype = new int[nb_variable - 1];
  for (i = 0;i < nb_variable - 1;i++) {
    variable[i] = i + 1;
    itype[i] = type[i + 1];
  }

  seq = new Markovian_sequences(nb_variable - 1 , itype , nb_sequence ,
                                identifier , length , false);

  seq->Sequences::select_variable(*this , variable);

  for (i = 0;i < seq->nb_variable;i++) {
    if (characteristics[i + 1]) {
      seq->characteristics[i] = new Sequence_characteristics(*(characteristics[i + 1]));
    }
  }

  delete [] variable;
  delete [] itype;

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Concatenation des variables d'objets Markovian_sequences.
 *
 *  arguments : reference sur un objet Format_error, nombre d'objets Markovian_sequences,
 *              pointeurs sur les objets Markovian_sequences, echantillon de reference pour les identificateurs.
 *
 *--------------------------------------------------------------*/

Markovian_sequences* Markovian_sequences::merge_variable(Format_error &error , int nb_sample ,
                                                         const Markovian_sequences **iseq , int ref_sample) const

{
  bool status = true;
  register int i , j , k , m;
  int inb_variable , *itype , *iidentifier , *psequence , *csequence;
  Markovian_sequences *seq;
  const Markovian_sequences **pseq;


  seq = 0;
  error.init();

  for (i = 0;i < nb_sample;i++) {
    if (iseq[i]->nb_sequence != nb_sequence) {
      status = false;
      ostringstream error_message;
      error_message << STAT_label[STATL_SAMPLE] << " " << i + 2 << ": "
                    << SEQ_error[SEQR_NB_SEQUENCE];
      error.update((error_message.str()).c_str());
    }

    else {
      for (j = 0;j < nb_sequence;j++) {
        if (iseq[i]->length[j] != length[j]) {
          status = false;
          ostringstream error_message;
          error_message << STAT_label[STATL_SAMPLE] << " " << i + 2 << ": "
                        << SEQ_label[SEQL_SEQUENCE] << " " << j + 1 << ": "
                        << SEQ_error[SEQR_SEQUENCE_LENGTH];
          error.update((error_message.str()).c_str());
        }
      }
    }
  }

  if ((ref_sample != I_DEFAULT) && ((ref_sample < 1) || (ref_sample > nb_sample + 1))) {
    status = false;
    error.update(STAT_error[STATR_SAMPLE_INDEX]);
  }

  if (status) {
    nb_sample++;
    pseq = new const Markovian_sequences*[nb_sample];

    pseq[0] = this;
    inb_variable = nb_variable;
    for (i = 1;i < nb_sample;i++) {
      pseq[i] = iseq[i - 1];
      inb_variable += iseq[i - 1]->nb_variable;
    }

    // calcul des identificateurs des sequences

    if (ref_sample == I_DEFAULT) {
      for (i = 0;i < nb_sequence;i++) {
        for (j = 1;j < nb_sample;j++) {
          if (pseq[j]->identifier[i] != pseq[0]->identifier[i]) {
            break;
          }
        }
        if (j < nb_sample) {
          break;
        }
      }

      if (i < nb_sequence) {
        iidentifier = 0;
      }
      else {
        iidentifier = pseq[0]->identifier;
      }
    }

    else {
      iidentifier = pseq[--ref_sample]->identifier;
    }

    itype = new int[inb_variable];
    inb_variable = 0;
    for (i = 0;i < nb_sample;i++) {
      for (j = 0;j < pseq[i]->nb_variable;j++) {
        itype[inb_variable++] = pseq[i]->type[j];
      }
    }

    seq = new Markovian_sequences(inb_variable , itype , nb_sequence ,
                                  iidentifier , length , false);
    delete [] itype;

    // copie des sequences

    for (i = 0;i < nb_sequence;i++) {
      inb_variable = 0;
      for (j = 0;j < nb_sample;j++) {
        for (k = 0;k < pseq[j]->nb_variable;k++) {
          psequence = seq->sequence[i][inb_variable++];
          csequence = pseq[j]->sequence[i][k];
          for (m = 0;m < length[i];m++) {
            *psequence++ = *csequence++;
          }
        }
      }
    }

    inb_variable = 0;
    for (i = 0;i < nb_sample;i++) {
      for (j = 0;j < pseq[i]->nb_variable;j++) {
        seq->max_value[inb_variable] = pseq[i]->max_value[j];
        seq->marginal[inb_variable] = new Histogram(*(pseq[i]->marginal[j]));

        if (pseq[i]->characteristics[j]) {
          seq->characteristics[inb_variable] = new Sequence_characteristics(*(pseq[i]->characteristics[j]));
        }
        inb_variable++;
      }
    }

    delete [] pseq;
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Copie des sequence avec extraction des temps de sejour initiaux.
 *
 *  argument : reference sur un objet Format_error.
 *
 *--------------------------------------------------------------*/

 Markovian_sequences* Markovian_sequences::initial_run_computation(Format_error &error) const

{
  register int i;
  Markovian_sequences *seq;


  error.init();

  for (i = 0;i < nb_variable;i++) {
    if ((characteristics[i]) && (characteristics[i]->initial_run)) {
      break;
    }
  }

  if (i < nb_variable) {
    seq = 0;
    error.update(SEQ_error[SEQR_INITIAL_RUN_ALREADY_BUILT]);
  }

  else {
    seq = new Markovian_sequences(*this , 'c' , ADD_INITIAL_RUN);
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Ajout d'une serie de vecteurs absorbants a la fin de chaque sequence.
 *
 *  arguments : reference sur un objet Format_error, longueur.
 *
 *--------------------------------------------------------------*/

 Markovian_sequences* Markovian_sequences::add_absorbing_run(Format_error &error , int length_value) const

{
  bool initial_run_flag;
  register int i , j , k;
  int *ilength , *psequence , *csequence;
  Markovian_sequences *seq;


  seq = 0;
  error.init();

  if ((length_value != I_DEFAULT) && ((length_value <= max_length) ||
       (length_value > max_length + MAX_ABSORBING_RUN_LENGTH))) {
    error.update(SEQ_error[SEQR_SEQUENCE_LENGTH]);
  }

  else {
    if (length_value == I_DEFAULT) {
      length_value = max_length + ABSORBING_RUN_LENGTH;
    }

    ilength = new int[nb_sequence];
    for (i = 0;i < nb_sequence;i++) {
      ilength[i] = length_value;
    }

    seq = new Markovian_sequences(nb_variable , nb_sequence , identifier ,
                                  ilength , false);
    delete [] ilength;

    for (i = 0;i < seq->nb_sequence;i++) {
      for (j = 0;j < seq->nb_variable;j++) {
        psequence = seq->sequence[i][j];
        csequence = sequence[i][j];
        for (k = 0;k < length[i];k++) {
          *psequence++ = *csequence++;
        }
        for (k = length[i];k < seq->length[i];k++) {
          *psequence++ = marginal[j]->nb_value;
        }
      }
    }

    for (i = 0;i < seq->nb_variable;i++) {
      seq->min_value[i] = min_value[i];
      seq->max_value[i] = max_value[i] + 1;
      seq->build_marginal_histogram(i);
    }

    initial_run_flag = false;
    for (i = 0;i < nb_variable;i++) {
      if ((characteristics[i]) && (characteristics[i]->initial_run)) {
        initial_run_flag = true;
        break;
      }
    }

    seq->build_characteristic(I_DEFAULT , true , initial_run_flag);
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Segmentation de sequences.
 *
 *  arguments : reference sur un objet Format_error,
 *              longueur des sequences.
 *
 *--------------------------------------------------------------*/

Markovian_sequences* Markovian_sequences::split(Format_error &error , int step) const

{
  register int i , j , k , m;
  int inb_sequence , last_length , nb_segment , *ilength , *psequence , *csequence;
  Markovian_sequences *seq;


  error.init();

  if ((step < 1) || (step > max_length)) {
    seq = 0;
    error.update(SEQ_error[SEQR_SEQUENCE_LENGTH]);
  }

  else {
    ilength = new int[cumul_length / step + nb_sequence];

    inb_sequence = 0;
    for (i = 0;i < nb_sequence;i++) {
      for (j = 0;j < length[i] / step;j++) {
        ilength[inb_sequence++] = step;
      }
      last_length = length[i] % step;
      if (last_length > 0) {
        ilength[inb_sequence++] = last_length;
      }
    }

#   ifdef DEBUG
    cout << "\nTEST: " << inb_sequence << " | " << cumul_length / step + nb_sequence << endl;
#   endif

    seq = new Markovian_sequences(nb_variable , inb_sequence , 0 , ilength , false);
    delete [] ilength;

    // copie des sequences

    inb_sequence = 0;
    for (i = 0;i < nb_sequence;i++) {
      nb_segment = (length[i] % step == 0 ? length[i] / step : length[i] / step + 1);
      for (j = 0;j < nb_variable;j++) {
        csequence = sequence[i][j];
        for (k = 0;k < nb_segment;k++) {
          psequence = seq->sequence[inb_sequence + k][j];
          for (m = 0;m < seq->length[inb_sequence + k];m++) {
            *psequence++ = *csequence++;
          }
        }
      }

      inb_sequence += nb_segment;
    }

    for (i = 0;i < seq->nb_variable;i++) {
      seq->min_value[i] = min_value[i];
      seq->max_value[i] = max_value[i];
      seq->marginal[i] = new Histogram(*(marginal[i]));
    }

    seq->build_characteristic();
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la quantite d'information en faisant l'hypothese
 *  de variables aleatoires independantes et equidistribuees.
 *
 *--------------------------------------------------------------*/

double Markovian_sequences::iid_information_computation() const

{
  register int i;
  double information = 0.;


  for (i = (((type[0] != STATE) || (nb_variable == 1)) ? 0 : 1);i < nb_variable;i++) {
    information += marginal[i]->information_computation();
  }

  return information;
}


/*--------------------------------------------------------------*
 *
 *  Calcul des echantillons correspondant a l'evolution des probabilites
 *  de rester dans un etat donne d'une chaine de Markov non-homogene.
 *
 *  argument : etat.
 *
 *--------------------------------------------------------------*/

void Markovian_sequences::self_transition_computation(int state)

{
  register int i , j;
  int num , denom , *pfrequency;
  double *ppoint;


  pfrequency = self_transition[state]->frequency;
  ppoint = self_transition[state]->point[0];

  for (i = 0;i < max_length - 1;i++) {
    num = 0;
    denom = 0;

    for (j = 0;j < nb_sequence;j++) {
      if (i < length[j] - 1) {
        if (sequence[j][0][i] == state) {
          if (sequence[j][0][i + 1] == state) {
            num++;
          }
          denom++;
        }
      }
    }

    *pfrequency++ = denom;
    if (denom > 0) {
      *ppoint++ = (double)num / (double)denom;
    }
    else {
      *ppoint++ = D_DEFAULT;
    }
  }

# ifdef DEBUG
  pfrequency = self_transition[state]->frequency;
  ppoint = self_transition[state]->point[0];
  double sum = 0.;

  for (i = 0;i < max_length - 1;i++) {
    sum += *pfrequency++ * *ppoint++;
  }

  cout << "\naverage self-transition probability : "
       << sum / self_transition->nb_element_computation() << endl;
# endif

}


/*--------------------------------------------------------------*
 *
 *  Calcul des echantillons correspondant a l'evolution des probabilites
 *  de rester dans un etat donne d'une chaine de Markov non-homogene.
 *
 *--------------------------------------------------------------*/

void Markovian_sequences::self_transition_computation()

{
  if (!self_transition) {
    register int i;


    state_variable_init();
    self_transition = new Self_transition*[marginal[0]->nb_value];

    for (i = 0;i < marginal[0]->nb_value;i++) {
      self_transition[i] = new Self_transition(max_length - 1);
      self_transition_computation(i);
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Calcul des echantillons correspondant a l'evolution des probabilites
 *  de rester dans un etat donne d'une chaine de Markov non-homogene.
 *
 *  argument : homogeneite des etats.
 *
 *--------------------------------------------------------------*/

void Markovian_sequences::self_transition_computation(bool *homogeneity)

{
  if (!self_transition) {
    register int i;


    state_variable_init();
    self_transition = new Self_transition*[marginal[0]->nb_value];

    for (i = 0;i < marginal[0]->nb_value;i++) {
      switch (homogeneity[i]) {
      case false :
        self_transition[i] = new Self_transition(max_length - 1);
        self_transition_computation(i);
        break;
      case true :
        self_transition[i] = 0;
        break;
      }
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Accumulation des observations (pour une variable donnee).
 *
 *  argument : indice de la variable.
 *
 *--------------------------------------------------------------*/

void Markovian_sequences::observation_histogram_computation(int variable)

{
  register int i , j;
  int *pfrequency , *pstate , *poutput;


  // initialisation des histogrammes

  for (i = 0;i < marginal[0]->nb_value;i++) {
    pfrequency = observation[variable][i]->frequency;
    for (j = 0;j < marginal[variable]->nb_value;j++) {
      *pfrequency++ = 0;
    }
  }

  // mise a jour des histogrammes

  for (i = 0;i < nb_sequence;i++) {
    pstate = sequence[i][0];
    poutput = sequence[i][variable];
    for (j = 0;j < length[i];j++) {
      (observation[variable][*pstate++]->frequency[*poutput++])++;
    }
  }

  // extraction des caracteristiques des histogrammes

  for (i = 0;i < marginal[0]->nb_value;i++) {
    if (!characteristics[variable]) {
      observation[variable][i]->nb_value_computation();
    }
    observation[variable][i]->offset_computation();
    observation[variable][i]->nb_element_computation();
    observation[variable][i]->max_computation();

    if (!characteristics[variable]) {
      observation[variable][i]->mean_computation();
      observation[variable][i]->variance_computation();
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Accumulation des observations.
 *
 *--------------------------------------------------------------*/

void Markovian_sequences::observation_histogram_computation()

{
  register int i;


  for (i = 1;i < nb_variable;i++) {
    observation_histogram_computation(i);
  }
}


/*--------------------------------------------------------------*
 *
 *  Construction des histogrammes correspondant aux probabilites
 *  d'observation.
 *
 *  argument : nombre d'etats.
 *
 *--------------------------------------------------------------*/

void Markovian_sequences::create_observation_histogram(int nb_state)

{
  if ((nb_variable > 1) && (!observation)) {
    register int i , j;


    observation = new Histogram**[nb_variable];
    observation[0] = 0;

    for (i = 1;i < nb_variable;i++) {
      observation[i] = new Histogram*[nb_state];
      for (j = 0;j < nb_state;j++) {
        observation[i][j] = new Histogram(marginal[i]->nb_value);
      }
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Construction des histogrammes correspondant aux probabilites
 *  d'observation.
 *
 *--------------------------------------------------------------*/

void Markovian_sequences::build_observation_histogram()

{
  create_observation_histogram(marginal[0]->nb_value);
  observation_histogram_computation();
}


/*--------------------------------------------------------------*
 *
 *  test recouvrement entre les valeurs observees dans chaque etat.
 *
 *  argument : indice de la variable.
 *
 *--------------------------------------------------------------*/

bool Markovian_sequences::test_hidden(int variable) const

{
  bool hidden = false , **occurrence;
  register int i , j;
  int nb_occur , *pstate , *poutput;


  if (variable > 0) {
    occurrence = new bool*[marginal[0]->nb_value];
    for (i = 0;i < marginal[0]->nb_value;i++) {
      occurrence[i] = new bool[marginal[variable]->nb_value];
      for (j = 0;j < marginal[variable]->nb_value;j++) {
        occurrence[i][j] = false;
      }
    }

    for (i = 0;i < nb_sequence;i++) {
      pstate = sequence[i][0];
      poutput = sequence[i][variable];
      for (j = 0;j < length[i];j++) {
        occurrence[*pstate++][*poutput++] = true;
      }
    }

    for (i = 0;i < marginal[variable]->nb_value;i++) {
      nb_occur = 0;
      for (j = 0;j < marginal[0]->nb_value;j++) {
        if (occurrence[j][i]) {
          nb_occur++;
        }
      }

      if (nb_occur > 1) {
        hidden = true;
        break;
      }
    }

    for (i = 0;i < marginal[0]->nb_value;i++) {
      delete [] occurrence[i];
    }
    delete [] occurrence;
  }

  return hidden;
}


/*--------------------------------------------------------------*
 *
 *  Extraction des probabilites de chaque valeur en fonction de l'index
 *  (pour une variable donnee).
 *
 *  argument : indice de la variable.
 *
 *--------------------------------------------------------------*/

void Markovian_sequences::build_index_value(int variable)

{
  register int i , j;
  int total , *frequency;


  // creation d'un objet Curves

  characteristics[variable]->index_value = new Curves(marginal[variable]->nb_value ,
                                                      max_length , true , false);
  frequency = new int[marginal[variable]->nb_value];

  // calcul des probabilites de chaque valeur en fonction de l'index

  for (i = 0;i < max_length;i++) {
    for (j = 0;j < marginal[variable]->nb_value;j++) {
      frequency[j] = 0;
    }

    for (j = 0;j < nb_sequence;j++) {
      if (i < length[j]) {
        frequency[sequence[j][variable][i]]++;
      }
    }

    total = 0;
    for (j = 0;j < marginal[variable]->nb_value;j++) {
      total += frequency[j];
    }
    characteristics[variable]->index_value->frequency[i] = total;
    for (j = 0;j < marginal[variable]->nb_value;j++) {
      characteristics[variable]->index_value->point[j][i] = (double)frequency[j] / (double)total;
    }
  }

  delete [] frequency;
}


/*--------------------------------------------------------------*
 *
 *  Construction des histogrammes du temps avant la premiere occurrence
 *  d'une valeur (pour une variable donnee).
 *
 *  argument : indice de la variable.
 *
 *--------------------------------------------------------------*/

void Markovian_sequences::build_first_occurrence_histogram(int variable)

{
  bool *occurrence;
  register int i , j;
  int nb_value , *psequence;
  Histogram **first_occurrence;


  // creation des histogrammes

  first_occurrence = new Histogram*[marginal[variable]->nb_value];
  for (i = 0;i < marginal[variable]->nb_value;i++) {
    first_occurrence[i] = new Histogram(max_length);
  }

/*  characteristics[variable]->first_occurrence = new Histogram*[marginal[variable]->nb_value];
  for (i = 0;i < marginal[variable]->nb_value;i++) {
    characteristics[variable]->first_occurrence[i] = new Histogram(max_length);
  } */

  // mise a jour des histogrammes

  occurrence = new bool[marginal[variable]->nb_value];

  for (i = 0;i < nb_sequence;i++) {
    nb_value = 0;
    for (j = 0;j < marginal[variable]->nb_value;j++) {
      occurrence[j] = false;
    }

    psequence = sequence[i][variable];
    for (j = 0;j < length[i];j++) {
      if (!occurrence[*psequence]) {
        occurrence[*psequence] = true;
        (first_occurrence[*psequence]->frequency[j])++;
//       (characteristics[variable]->first_occurrence[*psequence]->frequency[j])++;

        nb_value++;
        if (nb_value == marginal[variable]->nb_value) {
          break;
        }
      }

      psequence++;
    }
  }

  delete [] occurrence;

  // extraction des caracteristiques des histogrammes

  for (i = 0;i < marginal[variable]->nb_value;i++) {
    first_occurrence[i]->nb_value_computation();
    first_occurrence[i]->offset_computation();
    first_occurrence[i]->nb_element_computation();
    first_occurrence[i]->max_computation();
    first_occurrence[i]->mean_computation();
    first_occurrence[i]->variance_computation();
  }

/*  for (i = 0;i < marginal[variable]->nb_value;i++) {
    characteristics[variable]->first_occurrence[i]->nb_value_computation();
    characteristics[variable]->first_occurrence[i]->offset_computation();
    characteristics[variable]->first_occurrence[i]->nb_element_computation();
    characteristics[variable]->first_occurrence[i]->max_computation();
    characteristics[variable]->first_occurrence[i]->mean_computation();
    characteristics[variable]->first_occurrence[i]->variance_computation();
  } */

  characteristics[variable]->first_occurrence = new Histogram*[marginal[variable]->nb_value];
  for (i = 0;i < marginal[variable]->nb_value;i++) {
    characteristics[variable]->first_occurrence[i] = new Histogram(*(first_occurrence[i]));
    delete first_occurrence[i];
  }
  delete [] first_occurrence;
}


/*--------------------------------------------------------------*
 *
 *  Construction des histogrammes du temps de retour dans une valeur
 *  (pour une variable donnee).
 *
 *  argument : indice de la variable.
 *
 *--------------------------------------------------------------*/

void Markovian_sequences::build_recurrence_time_histogram(int variable)

{
  register int i , j;
  int *index , *psequence;
  Histogram **recurrence_time;


  // creation des histogrammes

  recurrence_time = new Histogram*[marginal[variable]->nb_value];
  for (i = 0;i < marginal[variable]->nb_value;i++) {
    recurrence_time[i] = new Histogram(max_length);
  }

/*  characteristics[variable]->recurrence_time = new Histogram*[marginal[variable]->nb_value];
  for (i = 0;i < marginal[variable]->nb_value;i++) {
    characteristics[variable]->recurrence_time[i] = new Histogram(max_length);
  } */

  // mise a jour des histogrammes

  index = new int[marginal[variable]->nb_value];

  for (i = 0;i < nb_sequence;i++) {
    for (j = 0;j < marginal[variable]->nb_value;j++) {
      index[j] = I_DEFAULT;
    }
    psequence = sequence[i][variable];

    for (j = 0;j < length[i];j++) {
      if (index[*psequence] != I_DEFAULT) {
        (recurrence_time[*psequence]->frequency[j - index[*psequence]])++;
//        (characteristics[variable]->recurrence_time[*psequence]->frequency[j - index[*psequence]])++;
      }
      index[*psequence++] = j;
    }
  }

  delete [] index;

  // extraction des caracteristiques des histogrammes

  for (i = 0;i < marginal[variable]->nb_value;i++) {
    recurrence_time[i]->nb_value_computation();
    recurrence_time[i]->offset_computation();
    recurrence_time[i]->nb_element_computation();
    recurrence_time[i]->max_computation();
    recurrence_time[i]->mean_computation();
    recurrence_time[i]->variance_computation();
  }

/*  for (i = 0;i < marginal[variable]->nb_value;i++) {
    characteristics[variable]->recurrence_time[i]->nb_value_computation();
    characteristics[variable]->recurrence_time[i]->offset_computation();
    characteristics[variable]->recurrence_time[i]->nb_element_computation();
    characteristics[variable]->recurrence_time[i]->max_computation();
    characteristics[variable]->recurrence_time[i]->mean_computation();
    characteristics[variable]->recurrence_time[i]->variance_computation();
  } */

  characteristics[variable]->recurrence_time = new Histogram*[marginal[variable]->nb_value];
  for (i = 0;i < marginal[variable]->nb_value;i++) {
    characteristics[variable]->recurrence_time[i] = new Histogram(*(recurrence_time[i]));
    delete recurrence_time[i];
  }
  delete [] recurrence_time;
}


/*--------------------------------------------------------------*
 *
 *  Construction des histogrammes du temps de sejour dans une valeur
 *  (pour une variable donnee).
 *
 *  arguments : indice de la variable, flag sur la creation des histogrammes
 *              de temps de sejour initial.
 *
 *--------------------------------------------------------------*/

void Markovian_sequences::build_sojourn_time_histogram(int variable , int initial_run_flag)

/* {
  characteristics[variable]->create_sojourn_time_histogram(max_length , initial_run_flag);
  sojourn_time_histogram_computation(variable);
} */

{
  register int i , j;
  int run_length , *psequence;
  Histogram **sojourn_time , **initial_run , **final_run;


  // creation des histogrammes

  sojourn_time = new Histogram*[marginal[variable]->nb_value];
  for (i = 0;i < marginal[variable]->nb_value;i++) {
    sojourn_time[i] = new Histogram(max_length + 1);
  }

  if (initial_run_flag) {
    initial_run = new Histogram*[marginal[variable]->nb_value];
    for (i = 0;i < marginal[variable]->nb_value;i++) {
      initial_run[i] = new Histogram(max_length + 1);
    }
  }

  final_run = new Histogram*[marginal[variable]->nb_value];
  for (i = 0;i < marginal[variable]->nb_value;i++) {
    final_run[i] = new Histogram(max_length + 1);
  }

  // mise a jour des histogrammes

  for (i = 0;i < nb_sequence;i++) {
    psequence = sequence[i][variable];
    run_length = 1;
    for (j = 1;j < length[i];j++) {
      if (*(psequence + 1) != *psequence) {
        if ((initial_run_flag) && (run_length == j)) {
          (initial_run[*psequence]->frequency[run_length])++;
        }
        else {
          (sojourn_time[*psequence]->frequency[run_length])++;
        }
        run_length = 0;
      }

      run_length++;
      psequence++;
    }

    if ((initial_run_flag) && (run_length == length[i])) {
      (initial_run[*psequence]->frequency[run_length])++;
    }
    (final_run[*psequence]->frequency[run_length])++;
  }

  // extraction des caracteristiques des histogrammes

  for (i = 0;i < marginal[variable]->nb_value;i++) {
    sojourn_time[i]->nb_value_computation();
    sojourn_time[i]->offset_computation();
    sojourn_time[i]->nb_element_computation();
    sojourn_time[i]->max_computation();
    sojourn_time[i]->mean_computation();
    sojourn_time[i]->variance_computation();
  }

  if (initial_run_flag) {
    for (i = 0;i < marginal[variable]->nb_value;i++) {
      initial_run[i]->nb_value_computation();
      initial_run[i]->offset_computation();
      initial_run[i]->nb_element_computation();
      initial_run[i]->max_computation();
      initial_run[i]->mean_computation();
      initial_run[i]->variance_computation();
    }
  }

  for (i = 0;i < marginal[variable]->nb_value;i++) {
    final_run[i]->nb_value_computation();
    final_run[i]->offset_computation();
    final_run[i]->nb_element_computation();
    final_run[i]->max_computation();
    final_run[i]->mean_computation();
    final_run[i]->variance_computation();
  }

  characteristics[variable]->sojourn_time = new Histogram*[marginal[variable]->nb_value];
  for (i = 0;i < marginal[variable]->nb_value;i++) {
    characteristics[variable]->sojourn_time[i] = new Histogram(*(sojourn_time[i]));
    delete sojourn_time[i];
  }
  delete [] sojourn_time;

  if (initial_run_flag) {
    characteristics[variable]->initial_run = new Histogram*[marginal[variable]->nb_value];
    for (i = 0;i < marginal[variable]->nb_value;i++) {
      characteristics[variable]->initial_run[i] = new Histogram(*(initial_run[i]));
      delete initial_run[i];
    }
    delete [] initial_run;
  }

  characteristics[variable]->final_run = new Histogram*[marginal[variable]->nb_value];
  for (i = 0;i < marginal[variable]->nb_value;i++) {
    characteristics[variable]->final_run[i] = new Histogram(*(final_run[i]));
    delete final_run[i];
  }
  delete [] final_run;
}


/*--------------------------------------------------------------*
 *
 *  Accumulation du temps de sejour dans une valeur
 *  (pour une variable donnee).
 *
 *  argument : indice de la variable.
 *
 *--------------------------------------------------------------*/

void Markovian_sequences::sojourn_time_histogram_computation(int variable)

{
  register int i , j;
  int run_length , *pfrequency , *psequence;


  // initialisation des histogrammes

  for (i = 0;i < marginal[variable]->nb_value;i++) {
    characteristics[variable]->sojourn_time[i]->offset = 1;
    characteristics[variable]->sojourn_time[i]->nb_value = characteristics[variable]->sojourn_time[i]->alloc_nb_value;

    pfrequency = characteristics[variable]->sojourn_time[i]->frequency;
    for (j = 0;j < characteristics[variable]->sojourn_time[i]->nb_value;j++) {
      *pfrequency++ = 0;
    }

    if (characteristics[variable]->initial_run) {
      characteristics[variable]->initial_run[i]->offset = 1;
      characteristics[variable]->initial_run[i]->nb_value = characteristics[variable]->initial_run[i]->alloc_nb_value;

      pfrequency = characteristics[variable]->initial_run[i]->frequency;
      for (j = 0;j < characteristics[variable]->initial_run[i]->nb_value;j++) {
        *pfrequency++ = 0;
      }
    }

    characteristics[variable]->final_run[i]->offset = 1;
    characteristics[variable]->final_run[i]->nb_value = characteristics[variable]->final_run[i]->alloc_nb_value;

    pfrequency = characteristics[variable]->final_run[i]->frequency;
    for (j = 0;j < characteristics[variable]->final_run[i]->nb_value;j++) {
      *pfrequency++ = 0;
    }
  }

  // mise a jour des histogrammes

  for (i = 0;i < nb_sequence;i++) {
    psequence = sequence[i][variable];
    run_length = 1;
    for (j = 1;j < length[i];j++) {
      if (*(psequence + 1) != *psequence) {
        if ((characteristics[variable]->initial_run) && (run_length == j)) {
          (characteristics[variable]->initial_run[*psequence]->frequency[run_length])++;
        }
        else {
          (characteristics[variable]->sojourn_time[*psequence]->frequency[run_length])++;
        }
        run_length = 0;
      }

      run_length++;
      psequence++;
    }

    if ((characteristics[variable]->initial_run) && (run_length == length[i])) {
      (characteristics[variable]->initial_run[*psequence]->frequency[run_length])++;
    }
    (characteristics[variable]->final_run[*psequence]->frequency[run_length])++;
  }

  // extraction des caracteristiques des histogrammes

  for (i = 0;i < marginal[variable]->nb_value;i++) {
    characteristics[variable]->sojourn_time[i]->nb_value_computation();
    characteristics[variable]->sojourn_time[i]->offset_computation();
    characteristics[variable]->sojourn_time[i]->nb_element_computation();
    characteristics[variable]->sojourn_time[i]->max_computation();
    characteristics[variable]->sojourn_time[i]->mean_computation();
    characteristics[variable]->sojourn_time[i]->variance_computation();
  }

  if (characteristics[variable]->initial_run) {
    for (i = 0;i < marginal[variable]->nb_value;i++) {
      characteristics[variable]->initial_run[i]->nb_value_computation();
      characteristics[variable]->initial_run[i]->offset_computation();
      characteristics[variable]->initial_run[i]->nb_element_computation();
      characteristics[variable]->initial_run[i]->max_computation();
      characteristics[variable]->initial_run[i]->mean_computation();
      characteristics[variable]->initial_run[i]->variance_computation();
    }
  }

  for (i = 0;i < marginal[variable]->nb_value;i++) {
    characteristics[variable]->final_run[i]->nb_value_computation();
    characteristics[variable]->final_run[i]->offset_computation();
    characteristics[variable]->final_run[i]->nb_element_computation();
    characteristics[variable]->final_run[i]->max_computation();
    characteristics[variable]->final_run[i]->mean_computation();
    characteristics[variable]->final_run[i]->variance_computation();
  }
}


/*--------------------------------------------------------------*
 *
 *  Accumulation du temps de sejour censures dans un etat.
 *
 *  arguments : pointeurs sur les histogrammes des temps de sejour initiaux et finaux et
 *              des longueurs des sequences dans le cas d'un seul etat visite.
 *
 *--------------------------------------------------------------*/

void Markovian_sequences::censored_sojourn_time_histogram_computation(Histogram **initial_run ,
                                                                      Histogram **final_run ,
                                                                      Histogram **single_run) const

{
  register int i , j;
  int *psequence;


  for (i = 0;i < nb_sequence;i++) {
    psequence = sequence[i][0];
    for (j = 1;j < length[i];j++) {
      if (*(psequence + 1) != *psequence) {
        (initial_run[*psequence]->frequency[j])++;
        break;
      }
      psequence++;
    }

    psequence = sequence[i][0] + length[i] - 1;
    if (j == length[i]) {
      (single_run[*psequence]->frequency[length[i]])++;
    }

    else {
      for (j = 1;j < length[i];j++) {
        if (*(psequence - 1) != *psequence) {
          (final_run[*psequence]->frequency[j])++;
          break;
        }
        psequence--;
      }
    }
  }

  // extraction des caracteristiques des histogrammes des temps de sejour censures

  for (i = 0;i < marginal[0]->nb_value;i++) {
    initial_run[i]->nb_value_computation();
    initial_run[i]->offset_computation();
    initial_run[i]->nb_element_computation();

#   ifdef DEBUG
    initial_run[i]->max_computation();
    initial_run[i]->mean_computation();
    initial_run[i]->variance_computation();

    os << initial_run[i];
#   endif

    final_run[i]->nb_value_computation();
    final_run[i]->offset_computation();
    final_run[i]->nb_element_computation();

#   ifdef DEBUG
    final_run[i]->max_computation();
    final_run[i]->mean_computation();
    final_run[i]->variance_computation();

    os << final_run[i];
#   endif

    single_run[i]->nb_value_computation();
    single_run[i]->offset_computation();
    single_run[i]->nb_element_computation();

#   ifdef DEBUG
    single_run[i]->max_computation();
    single_run[i]->mean_computation();
    single_run[i]->variance_computation();

    os << single_run[i];
#   endif

  }
}


/*--------------------------------------------------------------*
 *
 *  Construction des histogrammes du nombre de series d'une valeur par sequence
 *  (pour une variable donnee).
 *
 *  argument : indice de la variable.
 *
 *--------------------------------------------------------------*/

void Markovian_sequences::build_nb_run_histogram(int variable)

{
  register int i , j;
  int *psequence , *count;
  Histogram **nb_run;


  // creation des histogrammes

  nb_run = new Histogram*[marginal[variable]->nb_value];
  for (i = 0;i < marginal[variable]->nb_value;i++) {
    nb_run[i] = new Histogram((max_length % 2 == 0 ?
                               max_length / 2 : max_length / 2 + 1) + 1);
  }

/*  characteristics[variable]->nb_run = new Histogram*[marginal[variable]->nb_value];
  for (i = 0;i < marginal[variable]->nb_value;i++) {
    characteristics[variable]->nb_run[i] = new Histogram((max_length % 2 == 0 ?
                                                          max_length / 2 : max_length / 2 + 1) + 1);
  } */

  // mise a jour des histogrammes

  count = new int[marginal[variable]->nb_value];

  for (i = 0;i < nb_sequence;i++) {
    for (j = 0;j < marginal[variable]->nb_value;j++) {
      count[j] = 0;
    }

    psequence = sequence[i][variable];
    count[*psequence++]++;
    for (j = 1;j < length[i];j++) {
      if (*psequence != *(psequence - 1)) {
        count[*psequence]++;
      }
      psequence++;
    }

    for (j = 0;j < marginal[variable]->nb_value;j++) {
      (nb_run[j]->frequency[count[j]])++;
//      (characteristics[variable]->nb_run[j]->frequency[count[j]])++;
    }
  }

  delete [] count;

  // extraction des caracteristiques des histogrammes

  for (i = 0;i < marginal[variable]->nb_value;i++) {
    nb_run[i]->nb_value_computation();
    nb_run[i]->offset_computation();
    nb_run[i]->nb_element = nb_sequence;
    nb_run[i]->max_computation();
    nb_run[i]->mean_computation();
    nb_run[i]->variance_computation();
  }

/*  for (i = 0;i < marginal[variable]->nb_value;i++) {
    characteristics[variable]->nb_run[i]->nb_value_computation();
    characteristics[variable]->nb_run[i]->offset_computation();
    characteristics[variable]->nb_run[i]->nb_element = nb_sequence;
    characteristics[variable]->nb_run[i]->max_computation();
    characteristics[variable]->nb_run[i]->mean_computation();
    characteristics[variable]->nb_run[i]->variance_computation();
  } */

  characteristics[variable]->nb_run = new Histogram*[marginal[variable]->nb_value];
  for (i = 0;i < marginal[variable]->nb_value;i++) {
    characteristics[variable]->nb_run[i] = new Histogram(*(nb_run[i]));
    delete nb_run[i];
  }
  delete [] nb_run;
}


/*--------------------------------------------------------------*
 *
 *  Construction des histogrammes du nombre d'occurrences
 *  d'une valeur par sequence (pour une variable donnee).
 *
 *  argument : indice de la variable.
 *
 *--------------------------------------------------------------*/

void Markovian_sequences::build_nb_occurrence_histogram(int variable)

{
  register int i , j;
  int *psequence , *count;
  Histogram **nb_occurrence;


  // creation des histogrammes

  nb_occurrence = new Histogram*[marginal[variable]->nb_value];
  for (i = 0;i < marginal[variable]->nb_value;i++) {
    nb_occurrence[i] = new Histogram(max_length + 1);
  }

/*  characteristics[variable]->nb_occurrence = new Histogram*[marginal[variable]->nb_value];
  for (i = 0;i < marginal[variable]->nb_value;i++) {
    characteristics[variable]->nb_occurrence[i] = new Histogram(max_length + 1);
  } */

  // mise a jour des histogrammes

  count = new int[marginal[variable]->nb_value];

  for (i = 0;i < nb_sequence;i++) {
    for (j = 0;j < marginal[variable]->nb_value;j++) {
      count[j] = 0;
    }

    psequence = sequence[i][variable];
    for (j = 0;j < length[i];j++) {
      count[*psequence++]++;
    }

    for (j = 0;j < marginal[variable]->nb_value;j++) {
      (nb_occurrence[j]->frequency[count[j]])++;
//      (characteristics[variable]->nb_occurrence[j]->frequency[count[j]])++;
    }
  }

  delete [] count;

  // extraction des caracteristiques des histogrammes

  for (i = 0;i < marginal[variable]->nb_value;i++) {
    nb_occurrence[i]->nb_value_computation();
    nb_occurrence[i]->offset_computation();
    nb_occurrence[i]->nb_element = nb_sequence;
    nb_occurrence[i]->max_computation();
    nb_occurrence[i]->mean_computation();
    nb_occurrence[i]->variance_computation();
  }

/*  for (i = 0;i < marginal[variable]->nb_value;i++) {
    characteristics[variable]->nb_occurrence[i]->nb_value_computation();
    characteristics[variable]->nb_occurrence[i]->offset_computation();
    characteristics[variable]->nb_occurrence[i]->nb_element = nb_sequence;
    characteristics[variable]->nb_occurrence[i]->max_computation();
    characteristics[variable]->nb_occurrence[i]->mean_computation();
    characteristics[variable]->nb_occurrence[i]->variance_computation();
  } */

  characteristics[variable]->nb_occurrence = new Histogram*[marginal[variable]->nb_value];
  for (i = 0;i < marginal[variable]->nb_value;i++) {
    characteristics[variable]->nb_occurrence[i] = new Histogram(*(nb_occurrence[i]));
    delete nb_occurrence[i];
  }
  delete [] nb_occurrence;
}


/*--------------------------------------------------------------*
 *
 *  Extraction des caracteristiques d'un echantillon de sequences.
 *
 *  Calcul du nombre de valeurs et de l'histogramme des valeurs,
 *  extraction des probabilites des valeurs en fonction de l'index,
 *  construction des histogrammes du temps avant la premiere occurrence d'une valeur,
 *  des histogrammes du temps de retour dans une valeur,
 *  des histogrammes du temps de sejour dans une valeur,
 *  des histogrammes du nombre de series d'une valeur par sequence et
 *  des histogrammes du nombre d'occurrences d'une valeur par sequence.
 *
 *  argument : indice de la variable, flags sur la construction des histogrammes
 *             du temps de sejour et de temps de sejour initial.
 *
 *--------------------------------------------------------------*/

void Markovian_sequences::build_characteristic(int variable , bool sojourn_time_flag ,
                                               bool initial_run_flag)

{
  register int i , j;
  bool build;


  for (i = 0;i < nb_variable;i++) {
    if ((variable == I_DEFAULT) || (i == variable)) {
      build = true;

      if (marginal[i]->nb_value > NB_OUTPUT) {
        build = false;
      }

      else {
        for (j = 0;j < marginal[i]->nb_value;j++) {
          if (marginal[i]->frequency[j] == 0) {
            build = false;
            break;
          }
        }
      }

      if (build) {
        if (sojourn_time_flag) {
          characteristics[i] = new Sequence_characteristics(marginal[i]->nb_value);
        }

        build_index_value(i);

        build_first_occurrence_histogram(i);
        build_recurrence_time_histogram(i);

/*        if (sojourn_time_flag) {
          characteristics[i]->create_sojourn_time_histogram(max_length , initial_run_flag);
        }
        sojourn_time_histogram_computation(i); */

        switch (sojourn_time_flag) {
        case false :
          sojourn_time_histogram_computation(i);
          break;
        case true :
          build_sojourn_time_histogram(i , initial_run_flag);
          break;
        }

        if (max_length <= COUNT_MAX_LENGTH) {
          build_nb_run_histogram(i);
          build_nb_occurrence_histogram(i);
        }
      }
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Comptage des frequences des mots de longueur donnee.
 *
 *  argument : reference sur un objet Format_error, stream, indice de la variable,
 *             longueur des mots, etats de debut et de fin, frequence minimum.
 *
 *--------------------------------------------------------------*/

bool Markovian_sequences::word_count(Format_error &error , ostream &os , int variable ,
                                     int word_length , int begin_state , int end_state ,
                                     int min_frequency) const

{
  bool status = true , *selected_word;
  register int i , j , k;
  int nb_state , nb_word , max_nb_word , value , max_frequency , total_frequency , width ,
      *power , *frequency , *word_value , *psequence , *index , **word;
  double *probability;
  long old_adjust;


  error.init();

  if ((variable < 1) || (variable > nb_variable)) {
    status = false;
    error.update(STAT_error[STATR_VARIABLE_INDEX]);
  }

  else {
    variable--;
    nb_state = marginal[variable]->nb_value - marginal[variable]->offset;
 
    if ((nb_state < 2) || (nb_state > NB_STATE)) {
      status = false;
      error.update(SEQ_error[SEQR_NB_STATE]);
    }

    else {
      if (pow((double)nb_state , word_length) > MAX_NB_WORD) {
        status = false;
        error.update(SEQ_error[SEQR_MAX_NB_WORD]);
      }

      if ((begin_state != I_DEFAULT) && ((begin_state < marginal[variable]->offset) ||
           (begin_state >= marginal[variable]->nb_value) || (marginal[variable]->frequency[begin_state] == 0))) {
        status = false;
        ostringstream error_message;
        error_message << STAT_label[STATL_VARIABLE] << " " << variable + 1 << ": "
                      << STAT_label[STATL_STATE] << " " << begin_state << " "
                      << SEQ_error[SEQR_NOT_PRESENT];
        error.update((error_message.str()).c_str());
      }

      if ((end_state != I_DEFAULT) && ((end_state < marginal[variable]->offset) ||
           (end_state >= marginal[variable]->nb_value) || (marginal[variable]->frequency[end_state] == 0))) {
        status = false;
        ostringstream error_message;
        error_message << STAT_label[STATL_VARIABLE] << " " << variable + 1 << ": "
                      << STAT_label[STATL_STATE] << " " << end_state << " "
                      << SEQ_error[SEQR_NOT_PRESENT];
        error.update((error_message.str()).c_str());
      }
    }
  }

  if (min_frequency < 1) {
    status = false;
    error.update(SEQ_error[SEQR_MIN_FREQUENCY]);
  }

  if (status) {
    power = new int[word_length];

    i = 1;
    for (j = 0;j < word_length - 1;j++) {
      power[j] = i;
      i *= marginal[variable]->nb_value;
    }
    power[word_length - 1] = i;

    // calcul du nombre de mots

    max_nb_word = 0;
    for (i = MAX(hlength->offset , word_length);i < hlength->nb_value;i++) {
      max_nb_word += hlength->frequency[i] * (i - (word_length - 1));
    }
    nb_word = (int)pow((double)marginal[variable]->nb_value , word_length);
    if (nb_word < max_nb_word) {
      max_nb_word = nb_word; 
    }

    frequency = new int[max_nb_word];
    word_value = new int[max_nb_word];
    word = new int*[max_nb_word];

    nb_word = 0;
    for (i = 0;i < nb_sequence;i++) {
      for (j = 0;j < length[i] - (word_length - 1);j++) {
        if (((begin_state == I_DEFAULT) || (sequence[i][variable][j] == begin_state)) &&
            ((end_state == I_DEFAULT) || (sequence[i][variable][j + word_length - 1] == end_state))) {

          // calcul de la valeur du mot

          psequence = sequence[i][variable] + j;
          value = 0;
          for (k = 0;k < word_length;k++) {
            value += *psequence++ * power[k];
          }

          // recherche du mot

          for (k = 0;k < nb_word;k++) {
            if (value == word_value[k]) {
              frequency[k]++;
              break;
            }
          }

          // creation du mot

          if (k == nb_word) {
            frequency[nb_word] = 1;
            word_value[nb_word] = value;
            word[nb_word] = new int[word_length];

            psequence = sequence[i][variable] + j;
            for (k = 0;k < word_length;k++) {
              word[nb_word][k] = *psequence++;
            }
            nb_word++;
          }
        }
      }
    }

    // tri des mots par frequence decroissante

    index = new int[nb_word];
    selected_word = new bool[nb_word];
    probability = new double[nb_word];

    total_frequency = 0;
    for (i = 0;i < nb_word;i++) {
      total_frequency += frequency[i];
      selected_word[i] = false;
    }

    for (i = 0;i < nb_word;i++) {
      max_frequency = 0;
      for (j = 0;j < nb_word;j++) {
        if ((!selected_word[j]) && (frequency[j] > max_frequency)) {
          max_frequency = frequency[j];
          index[i] = j;
        }
      }
      if (frequency[index[i]] < min_frequency) {
        break;
      }

      selected_word[index[i]] = true;
      probability[index[i]] = (double)frequency[index[i]] / (double)total_frequency;

      if (i == 0) {
        width = column_width(max_frequency);
      }
    }

    // sortie

    old_adjust = os.setf(ios::right , ios::adjustfield);

    for (j = 0;j < i;j++) {
      for (k = 0;k < word_length;k++) {
        os << word[index[j]][k] << " ";
      }
      os << "   " << setw(width) << frequency[index[j]]
         << "   " << probability[index[j]];

      if (j == 0) {
        os << "   (" << nb_word << " " << SEQ_label[SEQL_WORDS] << ", "
           << total_frequency << ")";
      }
      os << endl;
    }

    os.setf((FMTFLAGS)old_adjust , ios::adjustfield);

    delete [] power;
    delete [] frequency;
    delete [] word_value;
    delete [] index;
    delete [] selected_word;
    delete [] probability;

    for (i = 0;i < nb_word;i++) {
      delete [] word[i];
    }
    delete [] word;
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Markovian_sequences.
 *
 *  arguments : stream, flag niveau de detail, flag fichier.
 *
 *--------------------------------------------------------------*/

ostream& Markovian_sequences::ascii_write(ostream &os , bool exhaustive , bool comment_flag) const

{
  register int i;


  os << nb_variable << " " << STAT_word[nb_variable == 1 ? STATW_VARIABLE : STATW_VARIABLES] << endl;

  if ((self_transition) && (exhaustive)) {
    for (i = 0;i < marginal[0]->nb_value;i++) {
      if (self_transition[i]) {
        os << "\n";
        if (comment_flag) {
          os << "# ";
        }
        os << "   | " << STAT_label[STATL_STATE] << " " << i << " - "
           << SEQ_label[SEQL_OBSERVED] << " " << SEQ_label[SEQL_SELF_TRANSITION] << endl;

        self_transition[i]->ascii_print(os , comment_flag);
      }
    }
  }

  for (i = 0;i < nb_variable;i++) {
    os << "\n" << STAT_word[STATW_VARIABLE] << " " << i + 1 << " : "
       << STAT_sequence_word[type[i]] << "   ";
    if (comment_flag) {
      os << "# ";
    }
    os << marginal[i]->nb_value << " ";

    switch (type[i]) {
    case STATE :
       os << STAT_label[marginal[i]->nb_value == 1 ? STATL_STATE : STATL_STATES] << endl;
       break;
    case INT_VALUE :
       os << STAT_label[marginal[i]->nb_value == 1 ? STATL_VALUE : STATL_VALUES] << endl;
       break;
    }

    os << "\n";
    if (comment_flag) {
      os << "# ";
    }
    os << STAT_label[type[i] == STATE ? STATL_STATE : STATL_VALUE] << " "
       << STAT_label[STATL_HISTOGRAM] << " - ";
    marginal[i]->ascii_characteristic_print(os , false , comment_flag);

    if ((marginal[i]->nb_value <= ASCII_NB_VALUE) || (exhaustive)) {
      os << "\n";
      if (comment_flag) {
        os << "# ";
      }
      os << "   | " << STAT_label[type[i] == STATE ? STATL_STATE : STATL_VALUE] << " "
         << STAT_label[STATL_HISTOGRAM] << endl;
      marginal[i]->ascii_print(os , comment_flag);
    }

    if (characteristics[i]) {
      characteristics[i]->ascii_print(os , type[i] , *hlength , exhaustive , comment_flag);
    }
  }

  os << "\n";
  if (comment_flag) {
    os << "# ";
  }
  os << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_HISTOGRAM] << " - ";
  hlength->ascii_characteristic_print(os , false , comment_flag);

  if (exhaustive) {
    os << "\n";
    if (comment_flag) {
      os << "# ";
    }
    os << "   | " << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_HISTOGRAM] << endl;
    hlength->ascii_print(os , comment_flag);
  }

  os << "\n";
  if (comment_flag) {
    os << "# ";
  }
  os << SEQ_label[SEQL_CUMUL_LENGTH] << ": " << cumul_length << endl;

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Markovian_sequences.
 *
 *  arguments : stream, flag niveau de detail.
 *
 *--------------------------------------------------------------*/

ostream& Markovian_sequences::ascii_write(ostream &os , bool exhaustive) const

{
  return ascii_write(os , exhaustive , false);
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Markovian_sequences dans un fichier.
 *
 *  arguments : reference sur un objet Format_error, path,
 *              flag niveau de detail.
 *
 *--------------------------------------------------------------*/

bool Markovian_sequences::ascii_write(Format_error &error , const char *path ,
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
    ascii_write(out_file , exhaustive , false);
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Markovian_sequences.
 *
 *  arguments : stream, format (ligne/colonne), flag niveau de detail.
 *
 *--------------------------------------------------------------*/

ostream& Markovian_sequences::ascii_data_write(ostream &os , char format , bool exhaustive) const

{
  ascii_write(os , exhaustive , false);
  ascii_print(os , format , false);

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Markovian_sequences dans un fichier.
 *
 *  arguments : reference sur un objet Format_error, path,
 *              format (ligne/colonne), flag niveau de detail.
 *
 *--------------------------------------------------------------*/

bool Markovian_sequences::ascii_data_write(Format_error &error , const char *path ,
                                           char format , bool exhaustive) const

{
  bool status = false;
  ofstream out_file(path);


  error.init();

  if (!out_file) {
    status = false;
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {
    status = true;
    if (format != 'a') {
      ascii_write(out_file , exhaustive , true);
    }
    ascii_print(out_file , format , true);
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Markovian_sequences dans un fichier au format tableur.
 *
 *  arguments : reference sur un objet Format_error, path.
 *
 *--------------------------------------------------------------*/

bool Markovian_sequences::spreadsheet_write(Format_error &error , const char *path) const

{
  bool status;
  register int i;
  Curves *smoothed_curves;
  ofstream out_file(path);


  error.init();

  if (!out_file) {
    status = false;
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {
    status = true;

    out_file << nb_variable << "\t" << STAT_word[nb_variable == 1 ? STATW_VARIABLE : STATW_VARIABLES] << endl;

    if (self_transition) {
      for (i = 0;i < marginal[0]->nb_value;i++) {
        if (self_transition[i]) {
          out_file << "\n\t" << STAT_label[STATL_STATE] << " " << i << " - "
                   << SEQ_label[SEQL_OBSERVED] << " " << SEQ_label[SEQL_SELF_TRANSITION] << endl;
          self_transition[i]->spreadsheet_print(out_file);

          smoothed_curves = new Curves(*(self_transition[i]) , 's');

          out_file << "\n\t" << STAT_label[STATL_STATE] << " " << i << " - "
                   << SEQ_label[SEQL_SMOOTHED] << " " << SEQ_label[SEQL_OBSERVED] << " " << SEQ_label[SEQL_SELF_TRANSITION] << endl;
          smoothed_curves->spreadsheet_print(out_file);

          delete smoothed_curves;
        }
      }
    }

    for (i = 0;i < nb_variable;i++) {
      out_file << "\n" << STAT_word[STATW_VARIABLE] << "\t" << i + 1 << "\t\t"
               << marginal[i]->nb_value << " ";

      switch (type[i]) {
      case STATE :
         out_file << STAT_label[marginal[i]->nb_value == 1 ? STATL_STATE : STATL_STATES] << endl;
         break;
      case INT_VALUE :
         out_file << STAT_label[marginal[i]->nb_value == 1 ? STATL_VALUE : STATL_VALUES] << endl;
         break;
      }

      out_file << "\n" << STAT_label[type[i] == STATE ? STATL_STATE : STATL_VALUE] << " "
               << STAT_label[STATL_HISTOGRAM] << "\t";
      marginal[i]->spreadsheet_characteristic_print(out_file);

      out_file << "\n\t" << STAT_label[type[i] == STATE ? STATL_STATE : STATL_VALUE] << " "
               << STAT_label[STATL_HISTOGRAM] << endl;
      marginal[i]->spreadsheet_print(out_file);

      if (characteristics[i]) {
        characteristics[i]->spreadsheet_print(out_file , type[i] , *hlength);
      }
    }

    out_file << "\n" << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_HISTOGRAM] << "\t";
    hlength->spreadsheet_characteristic_print(out_file);

    out_file << "\n\t" << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_HISTOGRAM] << endl;
    hlength->spreadsheet_print(out_file);

    out_file << "\n" << SEQ_label[SEQL_CUMUL_LENGTH] << "\t" << cumul_length << endl;
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Sortie Gnuplot d'un objet Markovian_sequences pour une variable donnee
 *  dans le cas d'absence de lois caracteristiques.
 *
 *  arguments : prefixe des fichiers, titre des figures, indice de la variable,
 *              nombre de variables.
 *
 *--------------------------------------------------------------*/

bool Markovian_sequences::plot_print(const char *prefix , const char *title , int variable ,
                                     int nb_variable) const

{
  bool status;
  register int i;
  const Histogram *phisto[1];
  ostringstream data_file_name;


  // ecriture du fichier de donnees

  data_file_name << prefix << variable + 1 << ".dat";

  phisto[0] = hlength;
  status = marginal[variable]->plot_print((data_file_name.str()).c_str() , 1 , phisto);

  if (status) {
    for (i = 0;i < 2;i++) {
      ostringstream file_name[2];

      switch (i) {

      case 0 : {
        if (nb_variable == 1) {
          file_name[0] << prefix << ".plot";
        }
        else {
          file_name[0] << prefix << variable + 1 << ".plot";
        }
        break;
      }

      case 1 : {
        if (nb_variable == 1) {
          file_name[0] << prefix << ".print";
        }
        else {
          file_name[0] << prefix << variable + 1 << ".print";
        }
        break;
      }
      }

      ofstream out_file((file_name[0].str()).c_str());

      if (i == 1) {
        out_file << "set terminal postscript" << endl;

        if (nb_variable == 1) {
          file_name[1] << label(prefix) << ".ps";
        }
        else {
          file_name[1] << label(prefix) << variable + 1 << ".ps";
        }
        out_file << "set output \"" << file_name[1].str() << "\"\n\n";
      }

      out_file << "set border 15 lw 0\n" << "set tics out\n" << "set xtics nomirror\n"
               << "set title";
      if (title) {
        out_file << " \"" << title << "\"";
      }
      out_file << "\n\n";

      if (marginal[variable]->nb_value - 1 < TIC_THRESHOLD) {
        out_file << "set xtics 0,1" << endl;
      }
      if ((int)(marginal[variable]->max * YSCALE) + 1 < TIC_THRESHOLD) {
        out_file << "set ytics 0,1" << endl;
      }

      out_file << "plot [0:" << MAX(marginal[variable]->nb_value - 1 , 1) << "] [0:"
               << (int)(marginal[variable]->max * YSCALE) + 1 << "] \""
               << label((data_file_name.str()).c_str()) << "\" using 1 title \""
               << STAT_label[STATL_VARIABLE] << " " << variable + 1 << " - "
               << STAT_label[type[variable] == STATE ? STATL_STATE : STATL_VALUE] << " "
               << STAT_label[STATL_HISTOGRAM] << "\" with impulses" << endl;

      if (marginal[variable]->nb_value - 1 < TIC_THRESHOLD) {
        out_file << "set xtics autofreq" << endl;
      }
      if ((int)(marginal[variable]->max * YSCALE) + 1 < TIC_THRESHOLD) {
        out_file << "set ytics autofreq" << endl;
      }

      if (i == 0) {
        out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
      }
      out_file << endl;

      if (hlength->nb_value - 1 < TIC_THRESHOLD) {
        out_file << "set xtics 0,1" << endl;
      }
      if ((int)(hlength->max * YSCALE) + 1 < TIC_THRESHOLD) {
        out_file << "set ytics 0,1" << endl;
      }

      out_file << "plot [0:" << hlength->nb_value - 1 << "] [0:"
               << (int)(hlength->max * YSCALE) + 1 << "] \""
               << label((data_file_name.str()).c_str()) << "\" using 2 title \""
               << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_HISTOGRAM]
               << "\" with impulses" << endl;

      if (hlength->nb_value - 1 < TIC_THRESHOLD) {
        out_file << "set xtics autofreq" << endl;
      }
      if ((int)(hlength->max * YSCALE) + 1 < TIC_THRESHOLD) {
        out_file << "set ytics autofreq" << endl;
      }

      if (i == 1) {
        out_file << "\nset terminal x11" << endl;
      }

      out_file << "\npause 0 \"" << STAT_label[STATL_END] << "\"" << endl;
    }
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Sortie Gnuplot d'un objet Markovian_sequences.
 *
 *  arguments : reference sur un objet Format_error, prefixe des fichiers,
 *              titre des figures.
 *
 *--------------------------------------------------------------*/

bool Markovian_sequences::plot_write(Format_error &error , const char *prefix ,
                                     const char *title) const

{
  bool status , start;
  register int i , j;
  int max_frequency[NB_OUTPUT];
  ostringstream data_file_name[NB_OUTPUT];


  error.init();

  if (characteristics[0]) {
    status = characteristics[0]->plot_print(prefix , title , 0 , nb_variable , type[0] , *hlength);
  }
  else {
    status = plot_print(prefix , title , 0 , nb_variable);
  }

  if (status) {
    for (i = 1;i < nb_variable;i++) {
      if (characteristics[i]) {
        characteristics[i]->plot_print(prefix , title , i , nb_variable , type[i] , *hlength);
      }
      else {
        plot_print(prefix , title , i , nb_variable);
      }
    }

    if (self_transition) {
      for (i = 0;i < marginal[0]->nb_value;i++) {
        if (self_transition[i]) {
          max_frequency[i] = self_transition[i]->max_frequency_computation();

          data_file_name[i] << prefix << i << ".dat";
          self_transition[i]->plot_print_0((data_file_name[i].str()).c_str());
        }
      }

      // ecriture des fichiers de commandes et des fichiers d'impression

      for (i = 0;i < 2;i++) {
        ostringstream file_name[2];

        switch (i) {
        case 0 :
          file_name[0] << prefix << 1 << 0 << ".plot";
          break;
        case 1 :
          file_name[0] << prefix << 1 << 0 << ".print";
          break;
        }

        ofstream out_file((file_name[0].str()).c_str());

        if (i == 1) {
          out_file << "set terminal postscript" << endl;
          file_name[1] << label(prefix) << 1 << 0 << ".ps";
          out_file << "set output \"" << file_name[1].str() << "\"\n\n";
        }

        out_file << "set border 15 lw 0\n" << "set tics out\n" << "set xtics nomirror\n"
                 << "set title";
        if ((title) || (nb_variable > 1)) {
          out_file << " \"";
          if (title) {
            out_file << title;
            if (nb_variable > 1) {
              out_file << " - ";
            }
          }
          if (nb_variable > 1) {
            out_file << STAT_label[STATL_VARIABLE] << " " << 1;
          }
          out_file << "\"";
        }
        out_file << "\n\n";

        start = true;
        for (j = 0;j < marginal[0]->nb_value;j++) {
          if (self_transition[j]) {
            if (!start) {
              if (i == 0) {
                out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
              }
              out_file << endl;
            }
            else {
              start = false;
            }

            if (self_transition[j]->length - 1 < TIC_THRESHOLD) {
              out_file << "set xtics 0,1" << endl;
            }

            out_file << "plot [0:" << self_transition[j]->length - 1 << "] [0:1] \""
                     << label((data_file_name[j].str()).c_str()) << "\" using 1:2 title \""
                     << STAT_label[STATL_STATE] << " " << j << " - " << SEQ_label[SEQL_OBSERVED] << " "
                     << SEQ_label[SEQL_SELF_TRANSITION] << "\" with points" << endl;

            if (self_transition[j]->length - 1 < TIC_THRESHOLD) {
              out_file << "set xtics autofreq" << endl;
            }

            if (i == 0) {
              out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
            }
            out_file << endl;

            out_file << "set xlabel \"" << SEQ_label[SEQL_INDEX] << "\"" << endl;
            out_file << "set ylabel \"" << STAT_label[STATL_FREQUENCY] << "\"" << endl;

            if (self_transition[j]->length - 1 < TIC_THRESHOLD) {
              out_file << "set xtics 0,1" << endl;
            }
            if ((int)(max_frequency[j] * YSCALE) + 1 < TIC_THRESHOLD) {
              out_file << "set ytics 0,1" << endl;
            }

            out_file << "plot [0:" << self_transition[j]->length - 1
                     << "] [0:" << (int)(max_frequency[j] * YSCALE) + 1 << "] \""
                     << label((data_file_name[j].str()).c_str())
                     << "\" using 1:3 notitle with impulses" << endl;

            if (self_transition[j]->length - 1 < TIC_THRESHOLD) {
              out_file << "set xtics autofreq" << endl;
            }
            if ((int)(max_frequency[j] * YSCALE) + 1 < TIC_THRESHOLD) {
              out_file << "set ytics autofreq" << endl;
            }

            out_file << "set xlabel" << endl;
            out_file << "set ylabel" << endl;
          }
        }

        if (i == 1) {
          out_file << "\nset terminal x11" << endl;
        }

        out_file << "\npause 0 \"" << STAT_label[STATL_END] << "\"" << endl;
      }
    }
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Fonctions pour la persistance.
 *
 *--------------------------------------------------------------*/

/* RWDEFINE_COLLECTABLE(Markovian_sequences , STATI_MARKOVIAN_SEQUENCES);


RWspace Markovian_sequences::binaryStoreSize() const

{
  register int i , j;
  RWspace size;


  size = Sequences::binaryStoreSize();

  size += sizeof(true);
  if (self_transition) {
    for (i = 0;i < marginal[0]->nb_value;i++) {
      size += sizeof(true);
      if (self_transition[i]) {
        size += self_transition[i]->binaryStoreSize();
      }
    }
  }

  size += sizeof(true);
  if (observation) {
    for (i = 1;i < nb_variable;i++) {
      for (j = 0;j < marginal[0]->nb_value;j++) {
        size += observation[i][j]->binaryStoreSize();
      }
    }
  }

  for (i = 0;i < nb_variable;i++) {
    size += sizeof(true);
    if (characteristics[i]) {
      characteristics[i]->binaryStoreSize();
    }
  }

  return size;
}


void Markovian_sequences::restoreGuts(RWvistream &is)

{
  bool status;
  register int i , j;


  remove();

  Sequences::restoreGuts(is);

  is >> status;
  if (status) {
    self_transition = new Curves*[marginal[0]->nb_value];
    for (i = 0;i < marginal[0]->nb_value;i++) {
      is >> status;
      if (status) {
        self_transition[i] = new Curves();
        self_transition[i]->restoreGuts(is);
      }
      else {
        self_transition[i] = 0;
      }
    }
  }
  else {
    self_transition = 0;
  }

  is >> status;
  if (status) {
    observation = new Histogram**[nb_variable];
    observation[0] = 0;
    for (i = 1;i < nb_variable;i++) {
      observation[i] = new Histogram*[marginal[0]->nb_value];
      for (j = 0;j < marginal[0]->nb_value;j++) {
        observation[i][j] = new Histogram();
        observation[i][j]->restoreGuts(is);
      }
    }
  }
  else {
    observation = 0;
  }

  characteristics = new Sequence_characteristics*[nb_variable];
  for (i = 0;i < nb_variable;i++) {
    is >> status;
    if (status) {
      characteristics[i] = new Sequence_characteristics();
      characteristics[i]->restoreGuts(is);
    }
    else {
      characteristics[i] = 0;
    }
  }
}


void Markovian_sequences::restoreGuts(RWFile &file)

{
  bool status;
  register int i , j;


  remove();

  Sequences::restoreGuts(file);

  file.Read(status);
  if (status) {
    self_transition = new Curves*[marginal[0]->nb_value];
    for (i = 0;i < marginal[0]->nb_value;i++) {
      file.Read(status);
      if (status) {
        self_transition[i] = new Curves();
        self_transition[i]->restoreGuts(file);
      }
      else {
        self_transition[i] = 0;
      }
    }
  }
  else {
    self_transition = 0;
  }

  file.Read(status);
  if (status) {
    observation = new Histogram**[nb_variable];
    observation[0] = 0;
    for (i = 1;i < nb_variable;i++) {
      observation[i] = new Histogram*[marginal[0]->nb_value];
      for (j = 0;j < marginal[0]->nb_value;j++) {
        observation[i][j] = new Histogram();
        observation[i][j]->restoreGuts(file);
      }
    }
  }
  else {
    observation = 0;
  }

  characteristics = new Sequence_characteristics*[nb_variable];
  for (i = 0;i < nb_variable;i++) {
    file.Read(status);
    if (status) {
      characteristics[i] = new Sequence_characteristics();
      characteristics[i]->restoreGuts(file);
    }
    else {
      characteristics[i] = 0;
    }
  }
}


void Markovian_sequences::saveGuts(RWvostream &os) const

{
  register int i , j;


  Sequences::saveGuts(os);

  if (self_transition) {
    os << true;
    for (i = 0;i < marginal[0]->nb_value;i++) {
      if (self_transition[i]) {
        os << true;
        self_transition[i]->saveGuts(os);
      }
      else {
        os << false;
      }
    }
  }
  else {
    os << false;
  }

  if (observation) {
    os << true;
    for (i = 1;i < nb_variable;i++) {
      for (j = 0;j < marginal[0]->nb_value;j++) {
        observation[i][j]->saveGuts(os);
      }
    }
  }
  else {
    os << false;
  }

  for (i = 0;i < nb_variable;i++) {
    if (characteristics[i]) {
      os << true;
      characteristics[i]->saveGuts(os);
    }
    else {
      os << false;
    }
  }
}


void Markovian_sequences::saveGuts(RWFile &file) const

{
  register int i , j;


  Sequences::saveGuts(file);

  if (self_transition) {
    file.Write(true);
    for (i = 0;i < marginal[0]->nb_value;i++) {
      if (self_transition[i]) {
        file.Write(true);
        self_transition[i]->saveGuts(file);
      }
      else {
        file.Write(false);
      }
    }
  }
  else {
    file.Write(false);
  }

  if (observation) {
    file.Write(true);
    for (i = 1;i < nb_variable;i++) {
      for (j = 0;j < marginal[0]->nb_value;j++) {
        observation[i][j]->saveGuts(file);
      }
    }
  }
  else {
    file.Write(false);
  }

  for (i = 0;i < nb_variable;i++) {
    if (characteristics[i]) {
      file.Write(true);
      characteristics[i]->saveGuts(file);
    }
    else {
      file.Write(false);
    }
  }
} */


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Markovian_sequences dans un fichier au format MTG.
 *
 *  arguments : reference sur un objet Format_error, path,
 *              type de chaque variable (NUMERIC/SYMBOLIC).
 *
 *--------------------------------------------------------------*/

bool Markovian_sequences::mtg_write(Format_error &error , const char *path , int *itype) const

{
  bool status;
  register int i , j , k , m;
  ofstream out_file(path);


  error.init();

  if (!out_file) {
    status = false;
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {
    status = true;

    // ecriture de l'entete

    out_file << "CODE:\tFORM-A" << endl;

    out_file << "\nCLASSES:\nSYMBOL\tSCALE\tDECOMPOSITION\tINDEXATION\tDEFINITION" << endl;
    out_file << "$\t0\tFREE\tFREE\tIMPLICIT" << endl;
    out_file << "U\t1\tFREE\tFREE\tIMPLICIT" << endl;
    out_file << "E\t2\tFREE\tFREE\tIMPLICIT" << endl;

    for (i = 0;i < nb_variable;i++) {
      switch (itype[i]) {

      case SYMBOLIC : {
        for (j = 1;j < marginal[i]->nb_value;j++) {
          out_file << (char)('F' + j) << "\t2\tFREE\tFREE\tIMPLICIT" << endl;
        }
        break;
      }

      case NUMERIC : {
        out_file << "F\t2\tFREE\tFREE\tIMPLICIT" << endl;
        break;
      }
      }
    }

    out_file << "\nDESCRIPTION:\nLEFT\tRIGHT\tRELTYPE\tMAX" << endl;
    out_file << "E\tE\t<\t1" << endl;

    for (i = 0;i < nb_variable;i++) {
      switch (itype[i]) {

      case SYMBOLIC : {
        out_file << "E\t";
        for (j = 1;j < marginal[i]->nb_value;j++) {
          out_file << (char)('F' + j);
          if (j < marginal[i]->nb_value - 1) {
            out_file << ",";
          }
        }
        out_file << "\t+\t1" << endl;
        break;
      }

      case NUMERIC : {
        out_file << "E\tF\t+\t?" << endl;
        break;
      }
      }
    }

    out_file << "\nFEATURES:\nNAME\tTYPE" << endl;

    // ecriture du code topologique

    out_file << "\nMTG:\nENTITY-CODE\n" << endl;

    for (i = 0;i < nb_sequence;i++) {
      out_file << "/U" << i + 1 << endl;

      for (j = 0;j < length[i];j++) {
        if (j == 0) {
          out_file << "\t/";
        }
        else {
          out_file << "\t^<";
        }
        out_file << 'E' << j + 1 << endl;

        for (k = 0;k < nb_variable;k++) {
          switch (itype[k]) {

          case SYMBOLIC : {
            if (sequence[i][k][j] > 0) {
              out_file <<"\t\t+" << (char)('F' + sequence[i][k][j]) << 1 << endl;
            }
            break;
          }

          case NUMERIC : {
            for (m = 0;m < sequence[i][k][j];m++) {
              out_file <<"\t\t+F" << m + 1 << endl;
            }
            break;
          }
          }
        }
      }
    }
  }

  return status;
}
