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



#include <sstream>
#include <iomanip>

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
 *  Constructeur par defaut de la classe MarkovianSequences.
 *
 *--------------------------------------------------------------*/

MarkovianSequences::MarkovianSequences()

{
  self_transition = NULL;
  observation = NULL;
  characteristics = NULL;
}


/*--------------------------------------------------------------*
 *
 *  Initialisation par defaut des champs de la classe MarkovianSequences.
 *
 *--------------------------------------------------------------*/

void MarkovianSequences::init()

{
  register int i;


  self_transition = NULL;
  observation = NULL;

  characteristics = new SequenceCharacteristics*[nb_variable];
  for (i = 0;i < nb_variable;i++) {
    characteristics[i] = NULL;
  }
}


/*--------------------------------------------------------------*
 *
 *  Construction d'un objet MarkovianSequences a partir d'un objet Sequences.
 *
 *  argument : reference sur un objet Sequences.
 *
 *--------------------------------------------------------------*/

MarkovianSequences::MarkovianSequences(const Sequences &seq)
:Sequences(seq)

{
  init();
  build_characteristic();
}


/*--------------------------------------------------------------*
 *
 *  Copie d'un objet MarkovianSequences.
 *
 *  arguments : reference sur un objet MarkovianSequences, inversion /
 *              ajout/suppression des lois empiriques de temps de sejour initial.
 *
 *--------------------------------------------------------------*/

void MarkovianSequences::copy(const MarkovianSequences &seq , int param)

{
  bool initial_run_flag;
  register int i , j;


  if ((seq.self_transition) && (param != REVERSE)) {
    self_transition = new SelfTransition*[marginal[0]->nb_value];
    for (i = 0;i < marginal[0]->nb_value;i++) {
      if (seq.self_transition[i]) {
        self_transition[i] = new SelfTransition(*(seq.self_transition[i]));
      }
      else {
        self_transition[i] = NULL;
      }
    }
  }

  else {
    self_transition = NULL;
  }

  if (seq.observation) {
    observation = new FrequencyDistribution**[nb_variable];
    observation[0] = NULL;

    for (i = 1;i < nb_variable;i++) {
      observation[i] = new FrequencyDistribution*[marginal[0]->nb_value];
      for (j = 0;j < marginal[0]->nb_value;j++) {
        observation[i][j] = new FrequencyDistribution(*(seq.observation[i][j]));
      }
    }
  }

  else {
    observation = NULL;
  }

  characteristics = new SequenceCharacteristics*[nb_variable];

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

        characteristics[i] = new SequenceCharacteristics(*(seq.characteristics[i]) , initial_run_flag);

        if (((seq.characteristics[i]->initial_run) && (!initial_run_flag)) ||
            ((!(seq.characteristics[i]->initial_run)) && (initial_run_flag))) {
           build_sojourn_time_frequency_distribution(i , initial_run_flag);
        }
      }

      else if (param == REVERSE) {
        characteristics[i] = new SequenceCharacteristics(*(seq.characteristics[i]), 'r');

        build_index_value(i);
        build_first_occurrence_frequency_distribution(i);

        if (!(seq.characteristics[i]->initial_run)) {
          build_sojourn_time_frequency_distribution(i);
        }
      }

      else {
        characteristics[i] = new SequenceCharacteristics(*(seq.characteristics[i]));
      }
    }

    else {
      characteristics[i] = NULL;
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Copie d'un objet MarkovianSequences avec ajout d'une variable d'etat.
 *
 *  arguments : reference sur un objet MarkovianSequences,
 *              ajout/suppression des lois empiriques de temps de sejour initial.
 *
 *--------------------------------------------------------------*/

void MarkovianSequences::add_state_variable(const MarkovianSequences &seq , int param)

{
  bool initial_run_flag;
  register int i;


  self_transition = NULL;
  observation = NULL;

  characteristics = new SequenceCharacteristics*[nb_variable];
  characteristics[0] = NULL;

  for (i = 0;i < seq.nb_variable;i++) {
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

        characteristics[i + 1] = new SequenceCharacteristics(*(seq.characteristics[i]) , initial_run_flag);

        if (((seq.characteristics[i]->initial_run) && (!initial_run_flag)) ||
            ((!(seq.characteristics[i]->initial_run)) && (initial_run_flag))) {
          build_sojourn_time_frequency_distribution(i + 1 , initial_run_flag);
        }
      }

      else {
        characteristics[i + 1] = new SequenceCharacteristics(*(seq.characteristics[i]));
      }
    }

    else {
      characteristics[i + 1] = NULL;
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Constructeur par copie de la classe MarkovianSequences.
 *
 *  arguments : reference sur un objet MarkovianSequences, type de transformation ('c' : copie,
 *              'a' : addition d'une variable d'etat, 'r' : suppression du parametre d'index),
 *              inversion / ajout/suppression des lois empiriques de temps de sejour initial.
 *
 *--------------------------------------------------------------*/

MarkovianSequences::MarkovianSequences(const MarkovianSequences &seq , char transform ,
                                       int param)

{
  switch (transform) {
  case 'c' :
    Sequences::copy(seq , (param == REVERSE ? true : false));
    copy(seq , param);
    break;
  case 'a' :
    Sequences::add_state_variable(seq);
    add_state_variable(seq , param);
    break;
  case 'r' :
    Sequences::remove_index_parameter(seq);
    copy(seq);
    break;
  default :
    Sequences::copy(seq);
    copy(seq);
    break;
  }
}


/*--------------------------------------------------------------*
 *
 *  Destructeur des champs de la classe MarkovianSequences.
 *
 *--------------------------------------------------------------*/

void MarkovianSequences::remove()

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
 *  Destructeur de la classe MarkovianSequences.
 *
 *--------------------------------------------------------------*/

MarkovianSequences::~MarkovianSequences()

{
  remove();
}


/*--------------------------------------------------------------*
 *
 *  Operateur d'assignement de la classe MarkovianSequences.
 *
 *  argument : reference sur un objet MarkovianSequences.
 *
 *--------------------------------------------------------------*/

MarkovianSequences& MarkovianSequences::operator=(const MarkovianSequences &seq)

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

void MarkovianSequences::state_variable_init(int itype)

{
  register int i , j;


  if (itype != type[0]) {
    if (type[0] == STATE) {
      if (self_transition) {
        for (i = 0;i < marginal[0]->nb_value;i++) {
          delete self_transition[i];
        }
        delete [] self_transition;

        self_transition = NULL;
      }

      if (observation) {
        for (i = 1;i < nb_variable;i++) {
          for (j = 0;j < marginal[0]->nb_value;j++) {
            delete observation[i][j];
          }
          delete [] observation[i];
        }
        delete [] observation;

        observation = NULL;
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
 *  Extraction d'un loi empirique.
 *
 *  arguments : reference sur un objet StatError, type de loi empirique,
 *              variable, valeur.
 *
 *--------------------------------------------------------------*/

DiscreteDistributionData* MarkovianSequences::extract(StatError &error , int type ,
                                                      int variable , int value) const

{
  bool status = true;
  FrequencyDistribution *phisto;
  DiscreteDistributionData *histo;


  histo = NULL;
  error.init();

  if ((variable < 1) || (variable > nb_variable)) {
    status = false;
    error.update(STAT_error[STATR_VARIABLE_INDEX]);
  }

  else {
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
          phisto = NULL;
          status = false;
          error.update(STAT_error[STATR_NON_EXISTING_FREQUENCY_DISTRIBUTION]);
        }
        break;
      }

      case FINAL_RUN : {
        phisto = characteristics[variable]->final_run[value];
        break;
      }

      case NB_RUN : {
        if (characteristics[variable]->nb_run) {
          phisto = characteristics[variable]->nb_run[value];
        }
        else {
          phisto = NULL;
          status = false;
          error.update(STAT_error[STATR_NON_EXISTING_FREQUENCY_DISTRIBUTION]);
        }
        break;
      }

      case NB_OCCURRENCE : {
        if (characteristics[variable]->nb_occurrence) {
          phisto = characteristics[variable]->nb_occurrence[value];
        }
        else {
          phisto = NULL;
          status = false;
          error.update(STAT_error[STATR_NON_EXISTING_FREQUENCY_DISTRIBUTION]);
        }
        break;
      }
      }

      if ((phisto) && (phisto->nb_element == 0)) {
        status = false;
        error.update(STAT_error[STATR_EMPTY_SAMPLE]);
      }
    }
  }

  if (status) {
    histo = new DiscreteDistributionData(*phisto);
  }

  return histo;
}


/*--------------------------------------------------------------*
 *
 *  Fusion d'objets MarkovianSequences.
 *
 *  arguments : reference sur un objet StatError, nombre d'objets MarkovianSequences,
 *              pointeurs sur les objets MarkovianSequences.
 *
 *--------------------------------------------------------------*/

MarkovianSequences* MarkovianSequences::merge(StatError &error , int nb_sample ,
                                              const MarkovianSequences **iseq) const

{
  bool status = true;
  register int i , j , k , m , n;
  int inb_sequence , nb_histo , *ilength , *pindex_param , *cindex_param ,
      *pisequence , *cisequence;
  double *prsequence , *crsequence;
  const FrequencyDistribution **phisto;
  MarkovianSequences *seq;
  const MarkovianSequences **pseq;


  seq = NULL;
  error.init();

  for (i = 0;i < nb_sample;i++) {
    if (iseq[i]->index_parameter_type != index_parameter_type) {
      status = false;
      ostringstream error_message;
      error_message << STAT_label[STATL_SAMPLE] << " " << i + 2 << ": "
                    << SEQ_error[SEQR_INDEX_PARAMETER_TYPE];

      if (index_parameter_type == IMPLICIT_TYPE) {
        error.update((error_message.str()).c_str());
      }
      else {
        error.correction_update((error_message.str()).c_str() , SEQ_index_parameter_word[index_parameter_type]);
      }
    }
  }

  for (i = 0;i < nb_sample;i++) {
    if (iseq[i]->nb_variable != nb_variable) {
      status = false;
      ostringstream error_message;
      error_message << STAT_label[STATL_SAMPLE] << " " << i + 2 << ": "
                    << STAT_error[STATR_NB_VARIABLE];
      error.correction_update((error_message.str()).c_str() , nb_variable);
    }

    else {
      for (j = 0;j < nb_variable;j++) {
        if (iseq[i]->type[j] != type[j]) {
          status = false;
          ostringstream error_message;
          error_message << STAT_label[STATL_SAMPLE] << " " << i + 2 << ": "
                        << STAT_label[STATL_VARIABLE] << " " << j + 1 << ": "
                        << STAT_error[STATR_VARIABLE_TYPE];
          error.correction_update((error_message.str()).c_str() , STAT_variable_word[type[j]]);
        }
      }
    }
  }

  if (status) {
    nb_sample++;
    pseq = new const MarkovianSequences*[nb_sample];

    pseq[0] = this;
    for (i = 1;i < nb_sample;i++) {
      pseq[i] = iseq[i - 1];
    }

    // calcul du nombre et de la longueur des sequences

    inb_sequence = 0;
    for (i = 0;i < nb_sample;i++) {
      inb_sequence += pseq[i]->nb_sequence;
    }

    ilength = new int[inb_sequence];

    i = 0;
    for (j = 0;j < nb_sample;j++) {
      for (k = 0;k < pseq[j]->nb_sequence;k++) {
        ilength[i++] = pseq[j]->length[k];
      }
    }

    seq = new MarkovianSequences(inb_sequence , 0 , ilength , index_parameter_type , 
                                 nb_variable , type);
    delete [] ilength;

    phisto = new const FrequencyDistribution*[nb_sample];

    // copie des sequences

    if (index_parameter_type == TIME) {
      i = 0;
      for (j = 0;j < nb_sample;j++) {
        for (k = 0;k < pseq[j]->nb_sequence;k++) {
          pindex_param = seq->index_parameter[i];
          cindex_param = pseq[j]->index_parameter[k];
          for (m = 0;m < pseq[j]->length[k];m++) {
            *pindex_param++ = *cindex_param++;
          }
          i++;
        }
      }

      for (i = 0;i < nb_sample;i++) {
        phisto[i] = pseq[i]->hindex_parameter;
      }
      seq->hindex_parameter = new FrequencyDistribution(nb_sample , phisto);

      for (i = 0;i < nb_sample;i++) {
        phisto[i] = pseq[i]->index_interval;
      }
      seq->index_interval = new FrequencyDistribution(nb_sample , phisto);
    }

    i = 0;
    for (j = 0;j < nb_sample;j++) {
      for (k = 0;k < pseq[j]->nb_sequence;k++) {
        for (m = 0;m < pseq[j]->nb_variable;m++) {
          if (pseq[j]->type[m] != REAL_VALUE) {
            pisequence = seq->int_sequence[i][m];
            cisequence = pseq[j]->int_sequence[k][m];
            for (n = 0;n < pseq[j]->length[k];n++) {
              *pisequence++ = *cisequence++;
            }
          }

          else {
            prsequence = seq->real_sequence[i][m];
            crsequence = pseq[j]->real_sequence[k][m];
            for (n = 0;n < pseq[j]->length[k];n++) {
              *prsequence++ = *crsequence++;
            }
          }
        }
        i++;
      }
    }

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

      if (seq->type[i] != REAL_VALUE) {
        for (j = 0;j < nb_sample;j++) {
          phisto[j] = pseq[j]->marginal[i];
        }
        seq->marginal[i] = new FrequencyDistribution(nb_sample , phisto);
      }

      for (j = 0;j < nb_sample;j++) {
        if (!(pseq[j]->characteristics[i])) {
          break;
        }
      }

      if (j == nb_sample) {
        seq->characteristics[i] = new SequenceCharacteristics();

        seq->characteristics[i]->nb_value = seq->marginal[i]->nb_value;

        seq->build_index_value(i);

        seq->characteristics[i]->first_occurrence = new FrequencyDistribution*[seq->marginal[i]->nb_value];
        seq->characteristics[i]->recurrence_time = new FrequencyDistribution*[seq->marginal[i]->nb_value];

        for (j = 0;j < seq->marginal[i]->nb_value;j++) {
          nb_histo = 0;
          for (k = 0;k < nb_sample;k++) {
            if (j < pseq[k]->marginal[i]->nb_value) {
              phisto[nb_histo++] = pseq[k]->characteristics[i]->first_occurrence[j];
            }
          }
          seq->characteristics[i]->first_occurrence[j] = new FrequencyDistribution(nb_histo , phisto);

          nb_histo = 0;
          for (k = 0;k < nb_sample;k++) {
            if (j < pseq[k]->marginal[i]->nb_value) {
              phisto[nb_histo++] = pseq[k]->characteristics[i]->recurrence_time[j];
            }
          }
          seq->characteristics[i]->recurrence_time[j] = new FrequencyDistribution(nb_histo , phisto);
        }

        for (j = 1;j < nb_sample;j++) {
          if (((pseq[0]->characteristics[i]->initial_run) && (!(pseq[j]->characteristics[i]->initial_run))) ||
              ((!(pseq[0]->characteristics[i]->initial_run)) && (pseq[j]->characteristics[i]->initial_run))) {
            break;
          }
        }

        if (j == nb_sample) {
          seq->characteristics[i]->sojourn_time = new FrequencyDistribution*[seq->marginal[i]->nb_value];
          if (pseq[0]->characteristics[i]->initial_run) {
            seq->characteristics[i]->initial_run = new FrequencyDistribution*[seq->marginal[i]->nb_value];
          }
          seq->characteristics[i]->final_run = new FrequencyDistribution*[seq->marginal[i]->nb_value];

          for (j = 0;j < seq->marginal[i]->nb_value;j++) {
            nb_histo = 0;
            for (k = 0;k < nb_sample;k++) {
              if (j < pseq[k]->marginal[i]->nb_value) {
                phisto[nb_histo++] = pseq[k]->characteristics[i]->sojourn_time[j];
              }
            }
            seq->characteristics[i]->sojourn_time[j] = new FrequencyDistribution(nb_histo , phisto);

            if (pseq[0]->characteristics[i]->initial_run) {
              nb_histo = 0;
              for (k = 0;k < nb_sample;k++) {
                if (j < pseq[k]->marginal[i]->nb_value) {
                  phisto[nb_histo++] = pseq[k]->characteristics[i]->initial_run[j];
                }
              }
              seq->characteristics[i]->initial_run[j] = new FrequencyDistribution(nb_histo , phisto);
            }

            nb_histo = 0;
            for (k = 0;k < nb_sample;k++) {
              if (j < pseq[k]->marginal[i]->nb_value) {
                phisto[nb_histo++] = pseq[k]->characteristics[i]->final_run[j];
              }
            }
            seq->characteristics[i]->final_run[j] = new FrequencyDistribution(nb_histo , phisto);
          }
        }

        else {
          seq->build_sojourn_time_frequency_distribution(i , (characteristics[i]->initial_run ? true : false));
        }

        for (j = 0;j < nb_sample;j++) {
          if ((!(pseq[j]->characteristics[i]->nb_run)) && (!(pseq[j]->characteristics[i]->nb_occurrence))) {
            break;
          }
        }

        if (j == nb_sample) {
          seq->characteristics[i]->nb_run = new FrequencyDistribution*[seq->marginal[i]->nb_value];
          seq->characteristics[i]->nb_occurrence = new FrequencyDistribution*[seq->marginal[i]->nb_value];

          for (j = 0;j < seq->marginal[i]->nb_value;j++) {
            nb_histo = 0;
            for (k = 0;k < nb_sample;k++) {
              if (j < pseq[k]->marginal[i]->nb_value) {
                phisto[nb_histo++] = pseq[k]->characteristics[i]->nb_run[j];
              }
            }
            seq->characteristics[i]->nb_run[j] = new FrequencyDistribution(nb_histo , phisto);

            nb_histo = 0;
            for (k = 0;k < nb_sample;k++) {
              if (j < pseq[k]->marginal[i]->nb_value) {
                phisto[nb_histo++] = pseq[k]->characteristics[i]->nb_occurrence[j];
              }
            }
            seq->characteristics[i]->nb_occurrence[j] = new FrequencyDistribution(nb_histo , phisto);
          }
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
 *  Regroupement des valeurs d'une variable.
 *
 *  arguments : reference sur un objet StatError, indice de la variable,
 *              pas de regroupement, mode (FLOOR/ROUND/CEIL).
 *
 *--------------------------------------------------------------*/

MarkovianSequences* MarkovianSequences::cluster(StatError &error , int variable ,
                                                int step , int mode) const

{
  bool status = true;
  register int i;
  MarkovianSequences *seq;


  seq = NULL;
  error.init();

  if ((variable < 1) || (variable > nb_variable)) {
    status = false;
    error.update(STAT_error[STATR_VARIABLE_INDEX]);
  }

  else {
    variable--;

    if ((type[variable] != INT_VALUE) && (type[variable] != REAL_VALUE)) {
      status = false;
      ostringstream correction_message;
      correction_message << STAT_variable_word[INT_VALUE] << " or " << STAT_variable_word[REAL_VALUE];
      error.correction_update(STAT_error[STATR_VARIABLE_TYPE] , (correction_message.str()).c_str());
    }
  }

  if (step < 1) {
    status = false;
    error.update(STAT_error[STATR_CLUSTERING_STEP]);
  }

  if (status) {
    seq = new MarkovianSequences(nb_sequence , identifier , length , index_parameter_type ,
                                 nb_variable , type);
    seq->Sequences::cluster(*this , variable , step , mode);

    for (i = 0;i < seq->nb_variable;i++) {
      if (i == variable) {
        seq->build_characteristic(i , true , (((characteristics[i]) && (characteristics[i]->initial_run)) ? true : false));
      }
      else if (characteristics[i]) {
        seq->characteristics[i] = new SequenceCharacteristics(*(characteristics[i]));
      }
    }
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Transcodage des symboles d'une variable entiere.
 *
 *  arguments : reference sur un objet StatError, indice de la variable,
 *              table de transcodage des symboles, flag pour ajouter une variable.
 *
 *--------------------------------------------------------------*/

MarkovianSequences* MarkovianSequences::transcode(StatError &error , int ivariable ,
                                                  int *symbol , bool add_flag) const

{
  bool status = true , *presence;
  register int i;
  int variable , offset , min_symbol , max_symbol , *itype;
  MarkovianSequences *seq;


  seq = NULL;
  error.init();

  if ((ivariable < 1) || (ivariable > nb_variable)) {
    status = false;
    error.update(STAT_error[STATR_VARIABLE_INDEX]);
  }

  else {
    ivariable--;

    if ((type[ivariable] != INT_VALUE) && (type[ivariable] != STATE)) {
      status = false;
      ostringstream correction_message;
      correction_message << STAT_variable_word[INT_VALUE] << " or " << STAT_variable_word[STATE];
      error.correction_update(STAT_error[STATR_VARIABLE_TYPE] , (correction_message.str()).c_str());
    }

    else {
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

      seq = new MarkovianSequences(nb_sequence , identifier , length , index_parameter_type ,
                                   nb_variable + offset , itype);
      delete [] itype;

      seq->Sequences::transcode(*this , ivariable , 0 , max_symbol , symbol , add_flag);

      for (i = 0;i < seq->nb_variable;i++) {
        if (i == variable) {
          seq->build_characteristic(i , true , (((characteristics[ivariable]) && (characteristics[ivariable]->initial_run)) ? true : false));
        }
        else if (characteristics[i - offset]) {
          seq->characteristics[i] = new SequenceCharacteristics(*(characteristics[i - offset]));
        }
      }
    }
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Transcodage des symboles d'une variable entiere.
 *
 *  arguments : references sur un objet StatError et
 *              sur un objet NonparametricSequenceProcess.
 *
 *--------------------------------------------------------------*/

MarkovianSequences* MarkovianSequences::transcode(StatError &error ,
                                                  const NonparametricSequenceProcess *process) const

{
  register int i , j;
  int *symbol;
  MarkovianSequences *seq;


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
 *  Suppression des valeurs non-representees d'une variable entiere.
 *
 *  arguments : reference sur un objet StatError, stream, indice de la variable,
 *              flag pour ajouter une variable.
 *
 *--------------------------------------------------------------*/

MarkovianSequences* MarkovianSequences::consecutive_values(StatError &error , ostream &os ,
                                                           int ivariable , bool add_flag) const

{
  bool status = true;
  register int i , j;
  int variable , offset , max , *symbol , *itype;
  MarkovianSequences *seq;


  seq = NULL;
  error.init();

  if ((ivariable < 1) || (ivariable > nb_variable)) {
    status = false;
    error.update(STAT_error[STATR_VARIABLE_INDEX]);
  }

  else {
    ivariable--;

    if ((type[ivariable] != INT_VALUE) && (type[ivariable] != STATE)) {
      status = false;
      ostringstream correction_message;
      correction_message << STAT_variable_word[INT_VALUE] << " or " << STAT_variable_word[STATE];
      error.correction_update(STAT_error[STATR_VARIABLE_TYPE] , (correction_message.str()).c_str());
    }

    else {
      for (i = 0;i < marginal[ivariable]->nb_value;i++) {
        if (marginal[ivariable]->frequency[i] == 0) {
          break;
        }
      }

      if (i == marginal[ivariable]->nb_value) {
        status = false;
        ostringstream error_message;
        error_message << STAT_label[STATL_VARIABLE] << " " << ivariable + 1 << ": "
                      << SEQ_error[SEQR_CONSECUTIVE_VALUES];
        error.update((error_message.str()).c_str());
      }
    }
  }

  if (status) {

#   ifdef MESSAGE
    os << "\n" << SEQ_label[SEQL_MISSING_VALUE] << ":";
    for (i = 0;i < marginal[ivariable]->nb_value;i++) {
      if (marginal[ivariable]->frequency[i] == 0) {
        os << " " << i;
      }
    }
    os << endl;
#   endif

    symbol = new int[marginal[ivariable]->nb_value - marginal[ivariable]->offset];

//    i = 0;
    i = -1;
    for (j = marginal[ivariable]->offset;j < marginal[ivariable]->nb_value;j++) {
//      symbol[j - marginal[ivariable]->offset] = i;
      if (marginal[ivariable]->frequency[j] > 0) {
        i++;
      }
      symbol[j - marginal[ivariable]->offset] = i;
    }
//    max = i - 1;
    max = i;

#   ifdef DEBUG
    cout << "\nTest :";
    for (i = 0;i < marginal[ivariable]->nb_value - marginal[ivariable]->offset;i++) {
      cout << " " << symbol[i];
    }
    cout << endl;
#   endif

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

    seq = new MarkovianSequences(nb_sequence , identifier , length , index_parameter_type ,
                                 nb_variable + offset , itype);
    delete [] itype;

    seq->Sequences::transcode(*this , ivariable , 0 , max , symbol , add_flag);
    delete [] symbol;

    for (i = 0;i < seq->nb_variable;i++) {
      if (i == variable) {
        seq->build_characteristic(i , true , (((characteristics[ivariable]) && (characteristics[ivariable]->initial_run)) ? true : false));
      }
      else if (characteristics[i - offset]) {
        seq->characteristics[i] = new SequenceCharacteristics(*(characteristics[i - offset]));
      }
    }
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Regroupement des symboles d'une variable entiere.
 *
 *  arguments : reference sur un objet StatError, indice de la variable,
 *              nombres de classes, bornes pour regrouper les symboles,
 *              flag pour ajouter une variable.
 *
 *--------------------------------------------------------------*/

MarkovianSequences* MarkovianSequences::cluster(StatError &error , int ivariable , int nb_class ,
                                                int *ilimit , bool add_flag) const

{
  bool status = true;
  register int i , j , k;
  int variable , offset , *symbol , *limit , *itype;
  MarkovianSequences *seq;


  seq = NULL;
  error.init();

  if ((ivariable < 1) || (ivariable > nb_variable)) {
    status = false;
    error.update(STAT_error[STATR_VARIABLE_INDEX]);
  }

  else {
    ivariable--;

    if ((type[ivariable] != INT_VALUE) && (type[ivariable] != STATE)) {
      status = false;
      ostringstream correction_message;
      correction_message << STAT_variable_word[INT_VALUE] << " or " << STAT_variable_word[STATE];
      error.correction_update(STAT_error[STATR_VARIABLE_TYPE] , (correction_message.str()).c_str());
    }

    else if ((nb_class < 2) || (nb_class >= marginal[ivariable]->nb_value)) {
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

      seq = new MarkovianSequences(nb_sequence , identifier , length , index_parameter_type ,
                                   nb_variable + offset , itype);
      delete [] itype;

      seq->Sequences::transcode(*this , ivariable , 0 , nb_class - 1 , symbol , add_flag);
      delete [] symbol;

      for (i = 0;i < seq->nb_variable;i++) {
        if (i == variable) {
          seq->build_characteristic(i , true , (((characteristics[ivariable]) && (characteristics[ivariable]->initial_run)) ? true : false));
        }
        else if (characteristics[i - offset]) {
          seq->characteristics[i] = new SequenceCharacteristics(*(characteristics[i - offset]));
        }
      }
    }

    delete [] limit;
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Regroupement des valeurs d'une variable reelle.
 *
 *  arguments : reference sur un objet StatError, indice de la variable,
 *              nombre de classes, bornes pour regrouper les valeurs.
 *
 *--------------------------------------------------------------*/

MarkovianSequences* MarkovianSequences::cluster(StatError &error , int variable ,
                                                int nb_class , double *ilimit) const

{
  bool status = true;
  register int i;
  int *itype;
  double *limit;
  MarkovianSequences *seq;


  seq = NULL;
  error.init();

  if ((variable < 1) || (variable > nb_variable)) {
    status = false;
    error.update(STAT_error[STATR_VARIABLE_INDEX]);
  }

  else {
    variable--;

    if (type[variable] != REAL_VALUE) {
      status = false;
      error.correction_update(STAT_error[STATR_VARIABLE_TYPE] , STAT_variable_word[REAL_VALUE]);
    }
  }

  if (nb_class < 2) {
    status = false;
    error.update(STAT_error[STATR_NB_CLASS]);
  }

  if (status) {
    limit = new double[nb_class + 1];
    limit[0] = min_value[variable];
    for (i = 1;i < nb_class;i++) {
      limit[i] = ilimit[i - 1];
    }
    limit[nb_class] = max_value[variable] + DOUBLE_ERROR;

    for (i = 0;i < nb_class;i++) {
      if (limit[i] >= limit[i + 1]) {
        status = false;
        error.update(STAT_error[STATR_CLUSTER_LIMIT]);
      }
    }

    if (status) {
      itype = new int[nb_variable];

      for (i = 0;i < nb_variable;i++) {
        itype[i] = type[i];
      }
      itype[variable] = INT_VALUE;

      seq = new MarkovianSequences(nb_sequence , identifier , length , index_parameter_type ,
                                   nb_variable , itype);
      delete [] itype;

      seq->Sequences::cluster(*this , variable , nb_class , limit);

      for (i = 0;i < seq->nb_variable;i++) {
        if (i == variable) {
          seq->build_characteristic(i);
        }
        else if (characteristics[i]) {
          seq->characteristics[i] = new SequenceCharacteristics(*(characteristics[i]));
        }
      }
    }

    delete [] limit;
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Suppression du parametre d'index.
 *
 *  argument : reference sur un objet StatError.
 *
 *--------------------------------------------------------------*/

MarkovianSequences* MarkovianSequences::remove_index_parameter(StatError &error) const

{
  MarkovianSequences *seq;


  error.init();

  if (!index_parameter) {
    seq = NULL;
    error.update(SEQ_error[SEQR_INDEX_PARAMETER_TYPE]);
  }
  else {
    seq = new MarkovianSequences(*this , 'r');
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Selection de variables.
 *
 *  arguments : reference sur un objet StatError, nombre de variables,
 *              indices des variables, flag pour conserver ou rejeter
 *              les variables selectionnees.
 *
 *--------------------------------------------------------------*/

MarkovianSequences* MarkovianSequences::select_variable(StatError &error , int inb_variable ,
                                                        int *ivariable , bool keep) const

{
  bool status = true , *selected_variable;
  register int i;
  int bnb_variable , *variable , *itype;
  MarkovianSequences *seq;


  seq = NULL;
  error.init();

  if ((inb_variable < 1) || (inb_variable > (keep ? nb_variable : nb_variable - 1))) {
    status = false;
    error.update(STAT_error[STATR_NB_SELECTED_VARIABLE]);
  }

  else {
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

    bnb_variable = (keep ? inb_variable : nb_variable - inb_variable);

    itype = new int[bnb_variable];
    for (i = 0;i < bnb_variable;i++) {
      itype[i] = type[variable[i]];
    }

    seq = new MarkovianSequences(nb_sequence , identifier , length , index_parameter_type ,
                                 bnb_variable , itype);

    seq->Sequences::select_variable(*this , variable);

    for (i = 0;i < seq->nb_variable;i++) {
      if (characteristics[variable[i]]) {
        seq->characteristics[i] = new SequenceCharacteristics(*(characteristics[variable[i]]));
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

MarkovianSequences* MarkovianSequences::remove_variable_1() const

{
  register int i;
  int *variable , *itype;
  MarkovianSequences *seq;


  variable = new int[nb_variable - 1];
  itype = new int[nb_variable - 1];
  for (i = 0;i < nb_variable - 1;i++) {
    variable[i] = i + 1;
    itype[i] = type[i + 1];
  }

  seq = new MarkovianSequences(nb_sequence , identifier , length , index_parameter_type ,
                               nb_variable - 1 , itype);

  seq->Sequences::select_variable(*this , variable);

  for (i = 0;i < seq->nb_variable;i++) {
    if (characteristics[i + 1]) {
      seq->characteristics[i] = new SequenceCharacteristics(*(characteristics[i + 1]));
    }
  }

  delete [] variable;
  delete [] itype;

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Concatenation des variables d'objets MarkovianSequences.
 *
 *  arguments : reference sur un objet StatError, nombre d'objets MarkovianSequences,
 *              pointeurs sur les objets MarkovianSequences, echantillon de reference pour les identificateurs.
 *
 *--------------------------------------------------------------*/

MarkovianSequences* MarkovianSequences::merge_variable(StatError &error , int nb_sample ,
                                                       const MarkovianSequences **iseq , int ref_sample) const

{
  bool status = true;
  register int i , j , k , m;
  int inb_variable , *iidentifier , *itype , *pisequence , *cisequence;
  double *prsequence , *crsequence;
  MarkovianSequences *seq;
  const MarkovianSequences **pseq;


  seq = NULL;
  error.init();

  for (i = 0;i < nb_sample;i++) {
    if ((iseq[i]->index_parameter_type != IMPLICIT_TYPE) &&
        (iseq[i]->index_parameter_type != index_parameter_type)) {
      status = false;
      ostringstream error_message;
      error_message << STAT_label[STATL_SAMPLE] << " " << i + 2 << ": "
                    << SEQ_error[SEQR_INDEX_PARAMETER_TYPE];

      if (index_parameter_type == IMPLICIT_TYPE) {
        error.update((error_message.str()).c_str());
      }
      else {
        error.correction_update((error_message.str()).c_str() , SEQ_index_parameter_word[index_parameter_type]);
      }
    }

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

        else if ((iseq[i]->index_parameter_type == TIME) &&
                 (iseq[i]->index_parameter_type == index_parameter_type)) {
          for (k = 0;k < length[j];k++) {
            if (iseq[i]->index_parameter[j][k] != index_parameter[j][k]) {
              status = false;
              ostringstream error_message;
              error_message << STAT_label[STATL_SAMPLE] << " " << i + 2 << ": "
                            << SEQ_label[SEQL_SEQUENCE] << " " << j + 1 << ": "
                            << SEQ_label[SEQL_INDEX_PARAMETER] << " " << k << ": "
                            << SEQ_error[SEQR_INDEX_PARAMETER];
              error.update((error_message.str()).c_str());
            }
          }
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
    pseq = new const MarkovianSequences*[nb_sample];

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
        itype[inb_variable] = pseq[i]->type[j];
        if ((inb_variable > 0) && (itype[inb_variable] == STATE)) {
          itype[inb_variable] = INT_VALUE;
        }
        inb_variable++;
      }
    }

    seq = new MarkovianSequences(nb_sequence , iidentifier , length , index_parameter_type ,
                                 inb_variable , itype);
    delete [] itype;

    // copie des sequences

    if (hindex_parameter) {
      seq->hindex_parameter = new FrequencyDistribution(*hindex_parameter);
    }
    if (index_interval) {
      seq->index_interval = new FrequencyDistribution(*index_interval);
    }

    if (index_parameter) {
      for (i = 0;i < nb_sequence;i++) {
        for (j = 0;j < length[i];j++) {
          seq->index_parameter[i][j] = index_parameter[i][j];
        }
      }
    }

    for (i = 0;i < nb_sequence;i++) {
      inb_variable = 0;
      for (j = 0;j < nb_sample;j++) {
        for (k = 0;k < pseq[j]->nb_variable;k++) {
          if (seq->type[inb_variable] != REAL_VALUE) {
            pisequence = seq->int_sequence[i][inb_variable++];
            cisequence = pseq[j]->int_sequence[i][k];
            for (m = 0;m < length[i];m++) {
              *pisequence++ = *cisequence++;
            }
          }

          else {
            prsequence = seq->real_sequence[i][inb_variable++];
            crsequence = pseq[j]->real_sequence[i][k];
            for (m = 0;m < length[i];m++) {
              *prsequence++ = *crsequence++;
            }
          }
        }
      }
    }

    inb_variable = 0;
    for (i = 0;i < nb_sample;i++) {
      for (j = 0;j < pseq[i]->nb_variable;j++) {
        seq->min_value[inb_variable] = pseq[i]->min_value[j];
        seq->max_value[inb_variable] = pseq[i]->max_value[j];
        if (pseq[i]->marginal[j]) {
          seq->marginal[inb_variable] = new FrequencyDistribution(*(pseq[i]->marginal[j]));
        }
        if (pseq[i]->characteristics[j]) {
          seq->characteristics[inb_variable] = new SequenceCharacteristics(*(pseq[i]->characteristics[j]));
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
 *  argument : reference sur un objet StatError.
 *
 *--------------------------------------------------------------*/

 MarkovianSequences* MarkovianSequences::initial_run_computation(StatError &error) const

{
  register int i;
  MarkovianSequences *seq;


  error.init();

  for (i = 0;i < nb_variable;i++) {
    if ((characteristics[i]) && (characteristics[i]->initial_run)) {
      break;
    }
  }

  if (i < nb_variable) {
    seq = NULL;
    error.update(SEQ_error[SEQR_INITIAL_RUN_ALREADY_BUILT]);
  }

  else {
    seq = new MarkovianSequences(*this , 'c' , ADD_INITIAL_RUN);
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Ajout d'une serie de vecteurs absorbants a la fin de chaque sequence.
 *
 *  arguments : reference sur un objet StatError, longueur des sequences,
 *              longueur des series finales absorbantes.
 *
 *--------------------------------------------------------------*/

MarkovianSequences* MarkovianSequences::add_absorbing_run(StatError &error ,
                                                          int sequence_length ,
                                                          int run_length) const

{
  bool status = true , initial_run_flag;
  register int i , j , k;
  int end_value , *ilength , *pisequence , *cisequence;
  double *prsequence , *crsequence;
  MarkovianSequences *seq;


  seq = NULL;
  error.init();

  if (index_parameter_type == TIME) {
    status = false;
    error.update(SEQ_error[SEQR_INDEX_PARAMETER_TYPE]);
  }

  if ((sequence_length != I_DEFAULT) && ((sequence_length <= max_length) ||
       (sequence_length > max_length + MAX_ABSORBING_RUN_LENGTH))) {
    status = false;
    error.update(SEQ_error[SEQR_SEQUENCE_LENGTH]);
  }

  if ((run_length != I_DEFAULT) && ((run_length < 1) ||
       (run_length > MAX_ABSORBING_RUN_LENGTH))) {
    status = false;
    error.update(SEQ_error[SEQR_RUN_LENGTH]);
  }

  if (status) {
    ilength = new int[nb_sequence];

    if (run_length == I_DEFAULT) {
     if (sequence_length == I_DEFAULT) {
        sequence_length = max_length + ABSORBING_RUN_LENGTH;
      }

      for (i = 0;i < nb_sequence;i++) {
        ilength[i] = sequence_length;
      }
    }

    else {
      for (i = 0;i < nb_sequence;i++) {
        ilength[i] = length[i] + run_length;
      }
    }

    seq = new MarkovianSequences(nb_sequence , identifier , ilength , index_parameter_type ,
                                 nb_variable , type);
    delete [] ilength;

    for (i = 0;i < seq->nb_sequence;i++) {
      for (j = 0;j < seq->nb_variable;j++) {
        if (seq->type[j] != REAL_VALUE) {
          pisequence = seq->int_sequence[i][j];
          cisequence = int_sequence[i][j];
          for (k = 0;k < length[i];k++) {
            *pisequence++ = *cisequence++;
          }

          if ((min_value[j] > 0) || (max_value[j] < 0)) {
            end_value = 0;
          }
          else {
            end_value = max_value[j] + 1;
          }

          for (k = length[i];k < seq->length[i];k++) {
            *pisequence++ = end_value;
          }

#         ifdef DEBUG

          // pour les donnees Fuji/Braeburn de successions d'UCs

          if (run_length == 1) {
            if ((j == 0) || (j == 1)) {
              pisequence--;
              *pisequence = *(pisequence - 1);
            }
            else if (j > 2) {
              *--pisequence = 0;
            }
          }
#         endif
        }

        else {
          prsequence = seq->real_sequence[i][j];
          crsequence = real_sequence[i][j];
          for (k = 0;k < length[i];k++) {
            *prsequence++ = *crsequence++;
          }

          if ((min_value[j] > 0.) || (max_value[j] < 0.)) {
            end_value = 0;
          }
          else {
            if (ceil(max_value[j]) > max_value[j]) {
              end_value = ceil(max_value[j]);
            }
            else {
              end_value = max_value[j] + 1;
            }
          }

          for (k = length[i];k < seq->length[i];k++) {
            *prsequence++ = end_value;
          }
        }
      }
    }

    for (i = 0;i < seq->nb_variable;i++) {
      seq->min_value_computation(i);
      seq->max_value_computation(i);
      seq->build_marginal_frequency_distribution(i);
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
 *  arguments : reference sur un objet StatError,
 *              longueur des sequences.
 *
 *--------------------------------------------------------------*/

MarkovianSequences* MarkovianSequences::split(StatError &error , int step) const

{
  register int i , j , k , m;
  int inb_sequence , last_length , nb_segment , *ilength , *pindex_param , *cindex_param ,
      *pisequence , *cisequence;
  double *prsequence , *crsequence;
  MarkovianSequences *seq;


  error.init();

  if ((step < 1) || (step > max_length)) {
    seq = NULL;
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

    seq = new MarkovianSequences(inb_sequence , 0 , ilength , index_parameter_type ,
                                 nb_variable , type);
    delete [] ilength;

    // copie des sequences

    inb_sequence = 0;
    for (i = 0;i < nb_sequence;i++) {
      nb_segment = (length[i] % step == 0 ? length[i] / step : length[i] / step + 1);

      if (seq->index_parameter_type == TIME) {
        cindex_param = index_parameter[i];
        for (j = 0;j < nb_segment;j++) {
          pindex_param = seq->index_parameter[inb_sequence + j];
          for (k = 0;k < seq->length[inb_sequence + j];k++) {
            *pindex_param++ = *cindex_param++;
          }
        }
      }

      for (j = 0;j < seq->nb_variable;j++) {
        if (seq->type[j] != REAL_VALUE) {
          cisequence = int_sequence[i][j];
          for (k = 0;k < nb_segment;k++) {
            pisequence = seq->int_sequence[inb_sequence + k][j];
            for (m = 0;m < seq->length[inb_sequence + k];m++) {
              *pisequence++ = *cisequence++;
            }
          }
        }

        else {
          crsequence = real_sequence[i][j];
          for (k = 0;k < nb_segment;k++) {
            prsequence = seq->real_sequence[inb_sequence + k][j];
            for (m = 0;m < seq->length[inb_sequence + k];m++) {
              *prsequence++ = *crsequence++;
            }
          }
        }
      }

      inb_sequence += nb_segment;
    }

    if (seq->index_parameter_type == TIME) {
      seq->hindex_parameter = new FrequencyDistribution(*hindex_parameter);
      seq->index_interval_computation();
    }

    for (i = 0;i < seq->nb_variable;i++) {
      seq->min_value[i] = min_value[i];
      seq->max_value[i] = max_value[i];
      if (marginal[i]) {
        seq->marginal[i] = new FrequencyDistribution(*marginal[i]);
      }
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

double MarkovianSequences::iid_information_computation() const

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

void MarkovianSequences::self_transition_computation(int state)

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
        if (int_sequence[j][0][i] == state) {
          if (int_sequence[j][0][i + 1] == state) {
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

void MarkovianSequences::self_transition_computation()

{
  if (!self_transition) {
    register int i;


    state_variable_init();
    self_transition = new SelfTransition*[marginal[0]->nb_value];

    for (i = 0;i < marginal[0]->nb_value;i++) {
      self_transition[i] = new SelfTransition(max_length - 1);
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

void MarkovianSequences::self_transition_computation(bool *homogeneity)

{
  if (!self_transition) {
    register int i;


    state_variable_init();
    self_transition = new SelfTransition*[marginal[0]->nb_value];

    for (i = 0;i < marginal[0]->nb_value;i++) {
      switch (homogeneity[i]) {
      case false :
        self_transition[i] = new SelfTransition(max_length - 1);
        self_transition_computation(i);
        break;
      case true :
        self_transition[i] = NULL;
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

void MarkovianSequences::observation_frequency_distribution_computation(int variable)

{
  register int i , j;
  int *pfrequency , *pstate , *poutput;


  // initialisation des lois empiriques

  for (i = 0;i < marginal[0]->nb_value;i++) {
    pfrequency = observation[variable][i]->frequency;
    for (j = 0;j < marginal[variable]->nb_value;j++) {
      *pfrequency++ = 0;
    }
  }

  // mise a jour des lois empiriques

  for (i = 0;i < nb_sequence;i++) {
    pstate = int_sequence[i][0];
    poutput = int_sequence[i][variable];
    for (j = 0;j < length[i];j++) {
      (observation[variable][*pstate++]->frequency[*poutput++])++;
    }
  }

  // extraction des caracteristiques des lois empiriques

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

void MarkovianSequences::observation_frequency_distribution_computation()

{
  register int i;


  for (i = 1;i < nb_variable;i++) {
    observation_frequency_distribution_computation(i);
  }
}


/*--------------------------------------------------------------*
 *
 *  Construction des lois empiriques d'observation.
 *
 *  argument : nombre d'etats.
 *
 *--------------------------------------------------------------*/

void MarkovianSequences::create_observation_frequency_distribution(int nb_state)

{
  if ((nb_variable > 1) && (!observation)) {
    register int i , j;


    observation = new FrequencyDistribution**[nb_variable];
    observation[0] = NULL;

    for (i = 1;i < nb_variable;i++) {
      observation[i] = new FrequencyDistribution*[nb_state];
      for (j = 0;j < nb_state;j++) {
        observation[i][j] = new FrequencyDistribution(marginal[i]->nb_value);
      }
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Construction des lois empiriques d'observation.
 *
 *--------------------------------------------------------------*/

void MarkovianSequences::build_observation_frequency_distribution()

{
  create_observation_frequency_distribution(marginal[0]->nb_value);
  observation_frequency_distribution_computation();
}


/*--------------------------------------------------------------*
 *
 *  test recouvrement entre les valeurs observees dans chaque etat.
 *
 *  argument : indice de la variable.
 *
 *--------------------------------------------------------------*/

bool MarkovianSequences::test_hidden(int variable) const

{
  bool hidden = true;

  if ((variable > 0) && (type[variable] == INT_VALUE)) {
    bool **occurrence;
    register int i , j;
    int nb_occurrence , *pstate , *poutput;


    hidden = false;

    occurrence = new bool*[marginal[0]->nb_value];
    for (i = 0;i < marginal[0]->nb_value;i++) {
      occurrence[i] = new bool[marginal[variable]->nb_value];
      for (j = 0;j < marginal[variable]->nb_value;j++) {
        occurrence[i][j] = false;
      }
    }

    for (i = 0;i < nb_sequence;i++) {
      pstate = int_sequence[i][0];
      poutput = int_sequence[i][variable];
      for (j = 0;j < length[i];j++) {
        occurrence[*pstate++][*poutput++] = true;
      }
    }

    for (i = 0;i < marginal[variable]->nb_value;i++) {
      nb_occurrence = 0;
      for (j = 0;j < marginal[0]->nb_value;j++) {
        if (occurrence[j][i]) {
          nb_occurrence++;
        }
      }

      if (nb_occurrence > 1) {
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
 *  Extraction des probabilites de chaque valeur entiere en fonction de l'index
 *  (pour une variable donnee).
 *
 *  argument : indice de la variable.
 *
 *--------------------------------------------------------------*/

void MarkovianSequences::build_index_value(int variable)

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
        frequency[int_sequence[j][variable][i]]++;
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
 *  Construction des lois empiriques du temps avant la 1ere occurrence
 *  d'une valeur entiere (pour une variable donnee).
 *
 *  argument : indice de la variable.
 *
 *--------------------------------------------------------------*/

void MarkovianSequences::build_first_occurrence_frequency_distribution(int variable)

{
  bool *occurrence;
  register int i , j;
  int nb_value , *pisequence;
  FrequencyDistribution **first_occurrence;


  // creation des lois empiriques

  first_occurrence = new FrequencyDistribution*[marginal[variable]->nb_value];
  for (i = 0;i < marginal[variable]->nb_value;i++) {
    first_occurrence[i] = new FrequencyDistribution(max_length);
  }

/*  characteristics[variable]->first_occurrence = new FrequencyDistribution*[marginal[variable]->nb_value];
  for (i = 0;i < marginal[variable]->nb_value;i++) {
    characteristics[variable]->first_occurrence[i] = new FrequencyDistribution(max_length);
  } */

  // mise a jour des lois empiriques

  occurrence = new bool[marginal[variable]->nb_value];

  for (i = 0;i < nb_sequence;i++) {
    nb_value = 0;
    for (j = 0;j < marginal[variable]->nb_value;j++) {
      occurrence[j] = false;
    }

    pisequence = int_sequence[i][variable];
    for (j = 0;j < length[i];j++) {
      if (!occurrence[*pisequence]) {
        occurrence[*pisequence] = true;
        (first_occurrence[*pisequence]->frequency[j])++;
//       (characteristics[variable]->first_occurrence[*pisequence]->frequency[j])++;

        nb_value++;
        if (nb_value == marginal[variable]->nb_value) {
          break;
        }
      }

      pisequence++;
    }
  }

  delete [] occurrence;

  // extraction des caracteristiques des lois empiriques

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

  characteristics[variable]->first_occurrence = new FrequencyDistribution*[marginal[variable]->nb_value];
  for (i = 0;i < marginal[variable]->nb_value;i++) {
    characteristics[variable]->first_occurrence[i] = new FrequencyDistribution(*(first_occurrence[i]));
    delete first_occurrence[i];
  }
  delete [] first_occurrence;
}


/*--------------------------------------------------------------*
 *
 *  Construction des lois empiriques du temps de retour dans une valeur entiere
 *  (pour une variable donnee).
 *
 *  argument : indice de la variable.
 *
 *--------------------------------------------------------------*/

void MarkovianSequences::build_recurrence_time_frequency_distribution(int variable)

{
  register int i , j;
  int *index , *pisequence;
  FrequencyDistribution **recurrence_time;


  // creation des lois empiriques

  recurrence_time = new FrequencyDistribution*[marginal[variable]->nb_value];
  for (i = 0;i < marginal[variable]->nb_value;i++) {
    recurrence_time[i] = new FrequencyDistribution(max_length);
  }

/*  characteristics[variable]->recurrence_time = new FrequencyDistribution*[marginal[variable]->nb_value];
  for (i = 0;i < marginal[variable]->nb_value;i++) {
    characteristics[variable]->recurrence_time[i] = new FrequencyDistribution(max_length);
  } */

  // mise a jour des lois empiriques

  index = new int[marginal[variable]->nb_value];

  for (i = 0;i < nb_sequence;i++) {
    for (j = 0;j < marginal[variable]->nb_value;j++) {
      index[j] = I_DEFAULT;
    }
    pisequence = int_sequence[i][variable];

    for (j = 0;j < length[i];j++) {
      if (index[*pisequence] != I_DEFAULT) {
        (recurrence_time[*pisequence]->frequency[j - index[*pisequence]])++;
//        (characteristics[variable]->recurrence_time[*pisequence]->frequency[j - index[*pisequence]])++;
      }
      index[*pisequence++] = j;
    }
  }

  delete [] index;

  // extraction des caracteristiques des lois empiriques

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

  characteristics[variable]->recurrence_time = new FrequencyDistribution*[marginal[variable]->nb_value];
  for (i = 0;i < marginal[variable]->nb_value;i++) {
    characteristics[variable]->recurrence_time[i] = new FrequencyDistribution(*(recurrence_time[i]));
    delete recurrence_time[i];
  }
  delete [] recurrence_time;
}


/*--------------------------------------------------------------*
 *
 *  Construction des lois empiriques du temps de sejour dans une valeur entiere
 *  (pour une variable donnee).
 *
 *  arguments : indice de la variable, flag sur la creation des lois empiriques
 *              de temps de sejour initial.
 *
 *--------------------------------------------------------------*/

void MarkovianSequences::build_sojourn_time_frequency_distribution(int variable , int initial_run_flag)

/* {
  characteristics[variable]->create_sojourn_time_frequency_distribution(max_length , initial_run_flag);
  sojourn_time_frequency_distribution_computation(variable);
} */

{
  register int i , j;
  int run_length , *pisequence;
  FrequencyDistribution **sojourn_time , **initial_run , **final_run;


  // creation des lois empiriques

  sojourn_time = new FrequencyDistribution*[marginal[variable]->nb_value];
  for (i = 0;i < marginal[variable]->nb_value;i++) {
    sojourn_time[i] = new FrequencyDistribution(max_length + 1);
  }

  if (initial_run_flag) {
    initial_run = new FrequencyDistribution*[marginal[variable]->nb_value];
    for (i = 0;i < marginal[variable]->nb_value;i++) {
      initial_run[i] = new FrequencyDistribution(max_length + 1);
    }
  }

  final_run = new FrequencyDistribution*[marginal[variable]->nb_value];
  for (i = 0;i < marginal[variable]->nb_value;i++) {
    final_run[i] = new FrequencyDistribution(max_length + 1);
  }

  // mise a jour des lois empiriques

  for (i = 0;i < nb_sequence;i++) {
    pisequence = int_sequence[i][variable];
    run_length = 1;
    for (j = 1;j < length[i];j++) {
      if (*(pisequence + 1) != *pisequence) {
        if ((initial_run_flag) && (run_length == j)) {
          (initial_run[*pisequence]->frequency[run_length])++;
        }
        else {
          (sojourn_time[*pisequence]->frequency[run_length])++;
        }
        run_length = 0;
      }

      run_length++;
      pisequence++;
    }

    if ((initial_run_flag) && (run_length == length[i])) {
      (initial_run[*pisequence]->frequency[run_length])++;
    }
    (final_run[*pisequence]->frequency[run_length])++;
  }

  // extraction des caracteristiques des lois empiriques

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

  characteristics[variable]->sojourn_time = new FrequencyDistribution*[marginal[variable]->nb_value];
  for (i = 0;i < marginal[variable]->nb_value;i++) {
    characteristics[variable]->sojourn_time[i] = new FrequencyDistribution(*(sojourn_time[i]));
    delete sojourn_time[i];
  }
  delete [] sojourn_time;

  if (initial_run_flag) {
    characteristics[variable]->initial_run = new FrequencyDistribution*[marginal[variable]->nb_value];
    for (i = 0;i < marginal[variable]->nb_value;i++) {
      characteristics[variable]->initial_run[i] = new FrequencyDistribution(*(initial_run[i]));
      delete initial_run[i];
    }
    delete [] initial_run;
  }

  characteristics[variable]->final_run = new FrequencyDistribution*[marginal[variable]->nb_value];
  for (i = 0;i < marginal[variable]->nb_value;i++) {
    characteristics[variable]->final_run[i] = new FrequencyDistribution(*(final_run[i]));
    delete final_run[i];
  }
  delete [] final_run;
}


/*--------------------------------------------------------------*
 *
 *  Accumulation du temps de sejour dans une valeur entiere
 *  (pour une variable donnee).
 *
 *  argument : indice de la variable.
 *
 *--------------------------------------------------------------*/

void MarkovianSequences::sojourn_time_frequency_distribution_computation(int variable)

{
  register int i , j;
  int run_length , *pfrequency , *pisequence;


  // initialisation des lois empiriques

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

  // mise a jour des lois empiriques

  for (i = 0;i < nb_sequence;i++) {
    pisequence = int_sequence[i][variable];
    run_length = 1;
    for (j = 1;j < length[i];j++) {
      if (*(pisequence + 1) != *pisequence) {
        if ((characteristics[variable]->initial_run) && (run_length == j)) {
          (characteristics[variable]->initial_run[*pisequence]->frequency[run_length])++;
        }
        else {
          (characteristics[variable]->sojourn_time[*pisequence]->frequency[run_length])++;
        }
        run_length = 0;
      }

      run_length++;
      pisequence++;
    }

    if ((characteristics[variable]->initial_run) && (run_length == length[i])) {
      (characteristics[variable]->initial_run[*pisequence]->frequency[run_length])++;
    }
    (characteristics[variable]->final_run[*pisequence]->frequency[run_length])++;
  }

  // extraction des caracteristiques des lois empiriques

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
 *  arguments : pointeurs sur les lois empiriques des temps de sejour initiaux et finaux et
 *              des longueurs des sequences dans le cas d'un seul etat visite.
 *
 *--------------------------------------------------------------*/

void MarkovianSequences::censored_sojourn_time_frequency_distribution_computation(FrequencyDistribution **initial_run ,
                                                                                  FrequencyDistribution **final_run ,
                                                                                  FrequencyDistribution **single_run) const

{
  register int i , j;
  int *pisequence;


  for (i = 0;i < nb_sequence;i++) {
    pisequence = int_sequence[i][0];
    for (j = 1;j < length[i];j++) {
      if (*(pisequence + 1) != *pisequence) {
        (initial_run[*pisequence]->frequency[j])++;
        break;
      }
      pisequence++;
    }

    pisequence = int_sequence[i][0] + length[i] - 1;
    if (j == length[i]) {
      (single_run[*pisequence]->frequency[length[i]])++;
    }

    else {
      for (j = 1;j < length[i];j++) {
        if (*(pisequence - 1) != *pisequence) {
          (final_run[*pisequence]->frequency[j])++;
          break;
        }
        pisequence--;
      }
    }
  }

  // extraction des caracteristiques des lois empiriques des temps de sejour censures

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
 *  Construction des lois empiriques du nombre de series d'une valeur entiere
 *  par sequence (pour une variable donnee).
 *
 *  argument : indice de la variable.
 *
 *--------------------------------------------------------------*/

void MarkovianSequences::build_nb_run_frequency_distribution(int variable)

{
  register int i , j;
  int *pisequence , *count;
  FrequencyDistribution **nb_run;


  // creation des lois empiriques

  nb_run = new FrequencyDistribution*[marginal[variable]->nb_value];
  for (i = 0;i < marginal[variable]->nb_value;i++) {
    nb_run[i] = new FrequencyDistribution((max_length % 2 == 0 ?
                               max_length / 2 : max_length / 2 + 1) + 1);
  }

/*  characteristics[variable]->nb_run = new FrequencyDistribution*[marginal[variable]->nb_value];
  for (i = 0;i < marginal[variable]->nb_value;i++) {
    characteristics[variable]->nb_run[i] = new FrequencyDistribution((max_length % 2 == 0 ?
                                                                     max_length / 2 : max_length / 2 + 1) + 1);
  } */

  // mise a jour des lois empiriques

  count = new int[marginal[variable]->nb_value];

  for (i = 0;i < nb_sequence;i++) {
    for (j = 0;j < marginal[variable]->nb_value;j++) {
      count[j] = 0;
    }

    pisequence = int_sequence[i][variable];
    count[*pisequence++]++;
    for (j = 1;j < length[i];j++) {
      if (*pisequence != *(pisequence - 1)) {
        count[*pisequence]++;
      }
      pisequence++;
    }

    for (j = 0;j < marginal[variable]->nb_value;j++) {
      (nb_run[j]->frequency[count[j]])++;
//      (characteristics[variable]->nb_run[j]->frequency[count[j]])++;
    }
  }

  delete [] count;

  // extraction des caracteristiques des lois empiriques

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

  characteristics[variable]->nb_run = new FrequencyDistribution*[marginal[variable]->nb_value];
  for (i = 0;i < marginal[variable]->nb_value;i++) {
    characteristics[variable]->nb_run[i] = new FrequencyDistribution(*(nb_run[i]));
    delete nb_run[i];
  }
  delete [] nb_run;
}


/*--------------------------------------------------------------*
 *
 *  Construction des lois empiriques du nombre d'occurrences
 *  d'une valeur entiere par sequence (pour une variable donnee).
 *
 *  argument : indice de la variable.
 *
 *--------------------------------------------------------------*/

void MarkovianSequences::build_nb_occurrence_frequency_distribution(int variable)

{
  register int i , j;
  int *pisequence , *count;
  FrequencyDistribution **nb_occurrence;


  // creation des lois empiriques

  nb_occurrence = new FrequencyDistribution*[marginal[variable]->nb_value];
  for (i = 0;i < marginal[variable]->nb_value;i++) {
    nb_occurrence[i] = new FrequencyDistribution(max_length + 1);
  }

/*  characteristics[variable]->nb_occurrence = new FrequencyDistribution*[marginal[variable]->nb_value];
  for (i = 0;i < marginal[variable]->nb_value;i++) {
    characteristics[variable]->nb_occurrence[i] = new FrequencyDistribution(max_length + 1);
  } */

  // mise a jour des lois empiriques

  count = new int[marginal[variable]->nb_value];

  for (i = 0;i < nb_sequence;i++) {
    for (j = 0;j < marginal[variable]->nb_value;j++) {
      count[j] = 0;
    }

    pisequence = int_sequence[i][variable];
    for (j = 0;j < length[i];j++) {
      count[*pisequence++]++;
    }

    for (j = 0;j < marginal[variable]->nb_value;j++) {
      (nb_occurrence[j]->frequency[count[j]])++;
//      (characteristics[variable]->nb_occurrence[j]->frequency[count[j]])++;
    }
  }

  delete [] count;

  // extraction des caracteristiques des lois empiriques

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

  characteristics[variable]->nb_occurrence = new FrequencyDistribution*[marginal[variable]->nb_value];
  for (i = 0;i < marginal[variable]->nb_value;i++) {
    characteristics[variable]->nb_occurrence[i] = new FrequencyDistribution(*(nb_occurrence[i]));
    delete nb_occurrence[i];
  }
  delete [] nb_occurrence;
}


/*--------------------------------------------------------------*
 *
 *  Extraction des caracteristiques d'un echantillon de sequences pour
 *  les variables entieres ayant un petit nombre de valeurs positives consecutives.
 *
 *  Calcul du nombre de valeurs et de la loi marginale empirique des valeurs,
 *  extraction des probabilites des valeurs en fonction de l'index,
 *  construction des lois empiriques du temps avant la 1ere occurrence d'une valeur,
 *  du temps de retour dans une valeur, du temps de sejour dans une valeur,
 *  du nombre de series d'une valeur par sequence et
 *  du nombre d'occurrences d'une valeur par sequence.
 *
 *  argument : indice de la variable, flags sur la construction des lois empiriques
 *             du temps de sejour et de temps de sejour initial.
 *
 *--------------------------------------------------------------*/

void MarkovianSequences::build_characteristic(int variable , bool sojourn_time_flag ,
                                              bool initial_run_flag)

{
  register int i , j;
  bool build;


  for (i = 0;i < nb_variable;i++) {
    if (((variable == I_DEFAULT) || (i == variable)) && (marginal[i])) {
      build = true;

      if (marginal[i]->nb_value > NB_OUTPUT) {
        build = false;
      }

      else if (type[i] != STATE) {
        for (j = 0;j < marginal[i]->nb_value;j++) {
          if (marginal[i]->frequency[j] == 0) {
            build = false;
            break;
          }
        }
      }

      if (build) {
        if (sojourn_time_flag) {
          characteristics[i] = new SequenceCharacteristics(marginal[i]->nb_value);
        }

        build_index_value(i);

        build_first_occurrence_frequency_distribution(i);
        build_recurrence_time_frequency_distribution(i);

/*        if (sojourn_time_flag) {
          characteristics[i]->create_sojourn_time_frequency_distribution(max_length , initial_run_flag);
        }
        sojourn_time_frequency_distribution_computation(i); */

        switch (sojourn_time_flag) {
        case false :
          sojourn_time_frequency_distribution_computation(i);
          break;
        case true :
          build_sojourn_time_frequency_distribution(i , initial_run_flag);
          break;
        }

        if (max_length <= COUNT_MAX_LENGTH) {
          build_nb_run_frequency_distribution(i);
          build_nb_occurrence_frequency_distribution(i);
        }
      }
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Comptage des frequences des mots de longueur donnee.
 *
 *  argument : reference sur un objet StatError, stream, indice de la variable,
 *             longueur des mots, etats de debut et de fin, frequence minimum.
 *
 *--------------------------------------------------------------*/

bool MarkovianSequences::word_count(StatError &error , ostream &os , int variable ,
                                    int word_length , int begin_state , int end_state ,
                                    int min_frequency) const

{
  bool status = true , *selected_word;
  register int i , j , k;
  int nb_state , nb_word , max_nb_word , value , max_frequency , total_frequency , width ,
      *power , *frequency , *word_value , *pisequence , *index , **word;
  double nb_word_bound , *probability;
  long old_adjust;


  error.init();

  if ((variable < 1) || (variable > nb_variable)) {
    status = false;
    error.update(STAT_error[STATR_VARIABLE_INDEX]);
  }

  else {
    variable--;

    if ((type[variable] != INT_VALUE) && (type[variable] != STATE)) {
      status = false;
      ostringstream correction_message;
      correction_message << STAT_variable_word[INT_VALUE] << " or " << STAT_variable_word[STATE];
      error.correction_update(STAT_error[STATR_VARIABLE_TYPE] , (correction_message.str()).c_str());
    }

    else {
      nb_state = marginal[variable]->nb_value - marginal[variable]->offset;
 
      if ((nb_state < 2) || (nb_state > NB_STATE)) {
        status = false;
        error.update(SEQ_error[SEQR_NB_STATE]);
      }

      else {
        max_nb_word = 0;
        for (i = MAX(hlength->offset , word_length);i < hlength->nb_value;i++) {
          max_nb_word += hlength->frequency[i] * (i - (word_length - 1));
        }
        nb_word_bound = pow((double)nb_state , word_length);
        if (nb_word_bound < max_nb_word) {
          max_nb_word = (int)nb_word_bound; 
        }

        if (max_nb_word > MAX_NB_WORD) {
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

    frequency = new int[max_nb_word];
    word_value = new int[max_nb_word];
    word = new int*[max_nb_word];

    nb_word = 0;
    for (i = 0;i < nb_sequence;i++) {
      for (j = 0;j < length[i] - (word_length - 1);j++) {
        if (((begin_state == I_DEFAULT) || (int_sequence[i][variable][j] == begin_state)) &&
            ((end_state == I_DEFAULT) || (int_sequence[i][variable][j + word_length - 1] == end_state))) {

          // calcul de la valeur du mot

          pisequence = int_sequence[i][variable] + j;
          value = 0;
          for (k = 0;k < word_length;k++) {
            value += *pisequence++ * power[k];
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

            pisequence = int_sequence[i][variable] + j;
            for (k = 0;k < word_length;k++) {
              word[nb_word][k] = *pisequence++;
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
 *  Ecriture d'un objet MarkovianSequences.
 *
 *  arguments : stream, flag niveau de detail, flag fichier.
 *
 *--------------------------------------------------------------*/

ostream& MarkovianSequences::ascii_write(ostream &os , bool exhaustive , bool comment_flag) const

{
  register int i;
  double mean , variance;


  if (index_parameter_type == TIME) {
    os << SEQ_word[SEQW_INDEX_PARAMETER] << " : "
       << SEQ_index_parameter_word[index_parameter_type] << "   ";
    if (comment_flag) {
      os << "# ";
    }
    os << "(" << SEQ_label[SEQL_MIN_INDEX_PARAMETER] << ": " << hindex_parameter->offset << ", "
       << SEQ_label[SEQL_MAX_INDEX_PARAMETER] << ": " << hindex_parameter->nb_value - 1 << ")" << endl;

    os << "\n";
    if (comment_flag) {
      os << "# ";
    }
    os << SEQ_label[SEQL_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";
    hindex_parameter->ascii_characteristic_print(os , false , comment_flag);

    if (exhaustive) {
      os << "\n";
      if (comment_flag) {
        os << "# ";
      }
      os << "   | " << SEQ_label[SEQL_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << endl;
      hindex_parameter->ascii_print(os , comment_flag);
    }
    os << endl;
  }

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
       << STAT_variable_word[type[i]] << "   ";
    if (comment_flag) {
      os << "# ";
    }

    if (type[i] == STATE) {
      os << "(" << marginal[i]->nb_value << " "
         << STAT_label[marginal[i]->nb_value == 1 ? STATL_STATE : STATL_STATES] << ")" << endl;
    }
    else {
      os << "(" << STAT_label[STATL_MIN_VALUE] << ": " << min_value[i] << ", "
         << STAT_label[STATL_MAX_VALUE] << ": " << max_value[i] << ")" << endl;
    }

    if (marginal[i]) {
      os << "\n";
      if (comment_flag) {
        os << "# ";
      }
      os << STAT_label[type[i] == STATE ? STATL_STATE : STATL_MARGINAL] << " "
         << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";
      marginal[i]->ascii_characteristic_print(os , false , comment_flag);

      if ((marginal[i]->nb_value <= ASCII_NB_VALUE) || (exhaustive)) {
        os << "\n";
        if (comment_flag) {
          os << "# ";
        }
        os << "   | " << STAT_label[type[i] == STATE ? STATL_STATE : STATL_MARGINAL] << " "
           << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << endl;
        marginal[i]->ascii_print(os , comment_flag);
      }
    }

    else {
      os << "\n";
      if (comment_flag) {
        os << "# ";
      }
      os << STAT_label[STATL_SAMPLE_SIZE] << ": " << cumul_length << endl;

      mean = mean_computation(i);
      variance = variance_computation(i , mean);

      if (comment_flag) {
        os << "# ";
      }
      os << STAT_label[STATL_MEAN] << ": " << mean << "   "
         << STAT_label[STATL_VARIANCE] << ": " << variance << "   "
         << STAT_label[STATL_STANDARD_DEVIATION] << ": " << sqrt(variance) << endl;

      if ((exhaustive) && (variance > 0.)) {
        if (comment_flag) {
          os << "# ";
        }
        os << STAT_label[STATL_SKEWNESS_COEFF] << ": " << skewness_computation(i , mean , variance) << "   "
           << STAT_label[STATL_KURTOSIS_COEFF] << ": " << kurtosis_computation(i , mean , variance) << endl;
      }
    }

    if (characteristics[i]) {
      characteristics[i]->ascii_print(os , type[i] , *hlength , exhaustive , comment_flag);
    }
  }

  os << "\n";
  if (comment_flag) {
    os << "# ";
  }
  os << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";
  hlength->ascii_characteristic_print(os , false , comment_flag);

  if (exhaustive) {
    os << "\n";
    if (comment_flag) {
      os << "# ";
    }
    os << "   | " << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << endl;
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
 *  Ecriture d'un objet MarkovianSequences.
 *
 *  arguments : stream, flag niveau de detail.
 *
 *--------------------------------------------------------------*/

ostream& MarkovianSequences::ascii_write(ostream &os , bool exhaustive) const

{
  return ascii_write(os , exhaustive , false);
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet MarkovianSequences dans un fichier.
 *
 *  arguments : reference sur un objet StatError, path,
 *              flag niveau de detail.
 *
 *--------------------------------------------------------------*/

bool MarkovianSequences::ascii_write(StatError &error , const char *path ,
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
 *  Ecriture d'un objet MarkovianSequences.
 *
 *  arguments : stream, format (ligne/colonne), flag niveau de detail.
 *
 *--------------------------------------------------------------*/

ostream& MarkovianSequences::ascii_data_write(ostream &os , char format , bool exhaustive) const

{
  ascii_write(os , exhaustive , false);
  ascii_print(os , format , false);

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet MarkovianSequences dans un fichier.
 *
 *  arguments : reference sur un objet StatError, path,
 *              format (ligne/colonne), flag niveau de detail.
 *
 *--------------------------------------------------------------*/

bool MarkovianSequences::ascii_data_write(StatError &error , const char *path ,
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
 *  Ecriture d'un objet MarkovianSequences dans un fichier au format tableur.
 *
 *  arguments : reference sur un objet StatError, path.
 *
 *--------------------------------------------------------------*/

bool MarkovianSequences::spreadsheet_write(StatError &error , const char *path) const

{
  bool status;
  register int i;
  double mean , variance;
  Curves *smoothed_curves;
  ofstream out_file(path);


  error.init();

  if (!out_file) {
    status = false;
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {
    status = true;

    if (index_parameter_type == TIME) {
      out_file << SEQ_word[SEQW_INDEX_PARAMETER] << "\t"
               << SEQ_index_parameter_word[index_parameter_type] << "\t\t"
               << SEQ_label[SEQL_MIN_INDEX_PARAMETER] << "\t" << hindex_parameter->offset << "\t\t"
               << SEQ_label[SEQL_MAX_INDEX_PARAMETER] << "\t" << hindex_parameter->nb_value - 1 << endl;

      out_file << "\n" << SEQ_label[SEQL_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\t";
      hindex_parameter->spreadsheet_characteristic_print(out_file);

      out_file << "\n\t" << SEQ_label[SEQL_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << endl;
      hindex_parameter->spreadsheet_print(out_file);
      out_file << endl;
    }

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
      out_file << "\n" << STAT_word[STATW_VARIABLE] << "\t" << i + 1 << "\t"
               << STAT_variable_word[type[i]];

      if (type[i] == STATE) {
        out_file << "\t\t" << marginal[i]->nb_value << "\t"
                 << STAT_label[marginal[i]->nb_value == 1 ? STATL_STATE : STATL_STATES] << endl;
      }
      else {
        out_file << "\t\t" << STAT_label[STATL_MIN_VALUE] << "\t" << min_value[i]
                 << "\t\t" << STAT_label[STATL_MAX_VALUE] << "\t" << max_value[i] << endl;
      }

      if (marginal[i]) {
        out_file << "\n" << STAT_label[type[i] == STATE ? STATL_STATE : STATL_MARGINAL] << " "
                 << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\t";
        marginal[i]->spreadsheet_characteristic_print(out_file);

        out_file << "\n\t" << STAT_label[type[i] == STATE ? STATL_STATE : STATL_MARGINAL] << " "
                 << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << endl;
        marginal[i]->spreadsheet_print(out_file);
      }

      else {
        out_file << "\n" << STAT_label[STATL_SAMPLE_SIZE] << "\t" << cumul_length << endl;

        mean = mean_computation(i);
        variance = variance_computation(i , mean);

        out_file << STAT_label[STATL_MEAN] << "\t" << mean << "\t\t"
                 << STAT_label[STATL_VARIANCE] << "\t" << variance << "\t\t"
                 << STAT_label[STATL_STANDARD_DEVIATION] << "\t" << sqrt(variance) << endl;

        if (variance > 0.) {
          out_file << STAT_label[STATL_SKEWNESS_COEFF] << "\t" << skewness_computation(i , mean , variance) << "\t\t"
                   << STAT_label[STATL_KURTOSIS_COEFF] << "\t" << kurtosis_computation(i , mean , variance) << endl;
        }
      }

      if (characteristics[i]) {
        characteristics[i]->spreadsheet_print(out_file , type[i] , *hlength);
      }
    }

    out_file << "\n" << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\t";
    hlength->spreadsheet_characteristic_print(out_file);

    out_file << "\n\t" << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << endl;
    hlength->spreadsheet_print(out_file);

    out_file << "\n" << SEQ_label[SEQL_CUMUL_LENGTH] << "\t" << cumul_length << endl;
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Sortie Gnuplot d'un objet MarkovianSequences pour une variable donnee
 *  dans le cas d'absence de lois caracteristiques.
 *
 *  arguments : prefixe des fichiers, titre des figures, indice de la variable,
 *              nombre de variables.
 *
 *--------------------------------------------------------------*/

bool MarkovianSequences::plot_print(const char *prefix , const char *title , int variable ,
                                    int nb_variable) const

{
  bool status;
  register int i;
  int nb_histo;
  const FrequencyDistribution *phisto[2];
  ostringstream data_file_name;


  // ecriture du fichier de donnees

  data_file_name << prefix << variable + 1 << ".dat";

  nb_histo = 0;

  phisto[nb_histo++] = hlength;
  if (hindex_parameter) {
    phisto[nb_histo++] = hindex_parameter;
  }

  status = marginal[variable]->plot_print((data_file_name.str()).c_str() , nb_histo , phisto);

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

      if (marginal[variable]) {
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
                 << STAT_label[type[variable] == STATE ? STATL_STATE : STATL_MARGINAL] << " "
                 << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\" with impulses" << endl;

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
      }

      if (hlength->nb_value - 1 < TIC_THRESHOLD) {
        out_file << "set xtics 0,1" << endl;
      }
      if ((int)(hlength->max * YSCALE) + 1 < TIC_THRESHOLD) {
        out_file << "set ytics 0,1" << endl;
      }

      out_file << "plot [0:" << hlength->nb_value - 1 << "] [0:"
               << (int)(hlength->max * YSCALE) + 1 << "] \""
               << label((data_file_name.str()).c_str()) << "\" using 2 title \""
               << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
               << "\" with impulses" << endl;

      if (hlength->nb_value - 1 < TIC_THRESHOLD) {
        out_file << "set xtics autofreq" << endl;
      }
      if ((int)(hlength->max * YSCALE) + 1 < TIC_THRESHOLD) {
        out_file << "set ytics autofreq" << endl;
      }

      if (hindex_parameter) {
        if (i == 0) {
          out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
        }
        out_file << endl;

        if (hindex_parameter->nb_value - 1 < TIC_THRESHOLD) {
          out_file << "set xtics 0,1" << endl;
        }
        if ((int)(hindex_parameter->max * YSCALE) + 1 < TIC_THRESHOLD) {
          out_file << "set ytics 0,1" << endl;
        }

        out_file << "plot [" << hindex_parameter->offset << ":"
                 << hindex_parameter->nb_value - 1 << "] [0:"
                 << (int)(hindex_parameter->max * YSCALE) + 1 << "] \""
                 << label((data_file_name.str()).c_str()) << "\" using 3 title \""
                 << SEQ_label[SEQL_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
                 << "\" with impulses" << endl;

        if (hindex_parameter->nb_value - 1 < TIC_THRESHOLD) {
          out_file << "set xtics autofreq" << endl;
        }
        if ((int)(hindex_parameter->max * YSCALE) + 1 < TIC_THRESHOLD) {
          out_file << "set ytics autofreq" << endl;
        }
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
 *  Sortie Gnuplot d'un objet MarkovianSequences.
 *
 *  arguments : reference sur un objet StatError, prefixe des fichiers,
 *              titre des figures.
 *
 *--------------------------------------------------------------*/

bool MarkovianSequences::plot_write(StatError &error , const char *prefix ,
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
          self_transition[i]->plot_print_standard_residual((data_file_name[i].str()).c_str());
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
                     << "\" using 1:3 title \"" << SEQ_label[SEQL_TRANSITION_COUNTS]
                     << "\" with impulses" << endl;

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
 *  Sortie graphique d'un objet MarkovianSequences pour une variable donnee
 *  dans le cas d'absence de lois caracteristiques.
 *
 *  arguments : indice de la variable, nombre de variables.
 *
 *--------------------------------------------------------------*/

void MarkovianSequences::plotable_write(MultiPlotSet &plot , int &index , int variable) const

{
  ostringstream legend;


  /*  nb_plot_set = 1;
  if (marginal[variable]) {
    nb_plot_set++;
  }
  if (hindex_parameter) {
    nb_plot_set++;
  } */

  plot.variable_nb_viewpoint[variable] = 1;

  if (marginal[variable]) {

    // vue : loi marginale empirique

    plot.variable[index] = variable;

    plot[index].xrange = Range(0 , MAX(marginal[variable]->nb_value - 1 , 1));
    plot[index].yrange = Range(0 , ceil(marginal[variable]->max * YSCALE));

    if (marginal[variable]->nb_value - 1 < TIC_THRESHOLD) {
      plot[index].xtics = 1;
    }
    if (ceil(marginal[variable]->max * YSCALE) < TIC_THRESHOLD) {
      plot[index].ytics = 1;
    }

    plot[index].resize(1);

    legend.str("");
    legend << STAT_label[STATL_VARIABLE] << " " << variable + 1 << " - "
           << STAT_label[STATL_MARGINAL] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
    plot[index][0].legend = legend.str();

    plot[index][0].style = "impulses";

    marginal[variable]->plotable_frequency_write(plot[index][0]);
    index++;
  }

  // vue : loi empirique des longueurs des sequences

  plot.variable[index] = variable;

  plot[index].xrange = Range(0 , hlength->nb_value - 1);
  plot[index].yrange = Range(0 , ceil(hlength->max * YSCALE));

  if (hlength->nb_value - 1 < TIC_THRESHOLD) {
    plot[index].xtics = 1;
  }
  if (ceil(hlength->max * YSCALE) < TIC_THRESHOLD) {
    plot[index].ytics = 1;
  }

  plot[index].resize(1);

  legend.str("");
  legend << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
  plot[index][0].legend = legend.str();

  plot[index][0].style = "impulses";

  hlength->plotable_frequency_write(plot[index][0]);
  index++;

  if (hindex_parameter) {

    // vue : loi empirique des parametres d'index

    plot.variable[index] = variable;

    plot[index].xrange = Range(hindex_parameter->offset , hindex_parameter->nb_value - 1);
    plot[index].yrange = Range(0 , ceil(hindex_parameter->max * YSCALE));

    if (hindex_parameter->nb_value - 1 < TIC_THRESHOLD) {
      plot[index].xtics = 1;
    }
    if (ceil(hindex_parameter->max * YSCALE) < TIC_THRESHOLD) {
      plot[index].ytics = 1;
    }

    plot[index].resize(1);

    legend.str("");
    legend << SEQ_label[SEQL_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
    plot[index][0].legend = legend.str();

    plot[index][0].style = "impulses";

    hindex_parameter->plotable_frequency_write(plot[index][0]);
    index++;
  }
}


/*--------------------------------------------------------------*
 *
 *  Sortie graphique d'un objet MarkovianSequences.
 *
 *--------------------------------------------------------------*/

MultiPlotSet* MarkovianSequences::get_plotable() const

{
  register int i , j;
  int nb_plot_set , index_length , index , max_frequency;
  ostringstream title , legend;
  MultiPlotSet *plot_set;


  // calcul du nombre de vues

  nb_plot_set = 0;

  for (i = 0;i < nb_variable;i++) {
    if (characteristics[i]) {
      index_length = characteristics[i]->index_value->plot_length_computation();

      nb_plot_set += 2;
      if (characteristics[i]->index_value->frequency[index_length - 1] < MAX_FREQUENCY) {
        nb_plot_set++;
      }

      nb_plot_set++;
      for (j = 0;j < characteristics[i]->nb_value;j++) {
        if (characteristics[i]->first_occurrence[j]->nb_element > 0) {
          nb_plot_set++;
        }
      }

      nb_plot_set++;
      for (j = 0;j < characteristics[i]->nb_value;j++) {
        if (characteristics[i]->recurrence_time[j]->nb_element > 0) {
          nb_plot_set++;
        }
      }

      nb_plot_set++;
      for (j = 0;j < characteristics[i]->nb_value;j++) {
        if (characteristics[i]->sojourn_time[j]->nb_element > 0) {
          nb_plot_set++;
        }
        if ((characteristics[i]->initial_run) &&
            (characteristics[i]->initial_run[j]->nb_element > 0)) {
          nb_plot_set++;
        }
        if (characteristics[i]->final_run[j]->nb_element > 0) {
          nb_plot_set++;
        }
      }

      if ((characteristics[i]->nb_run) && (characteristics[i]->nb_occurrence)) {
        nb_plot_set += 3;
        for (j = 0;j < characteristics[i]->nb_value;j++) {
          if ((characteristics[i]->nb_run[j]->nb_element > 0) &&
              (characteristics[i]->nb_occurrence[j]->nb_element > 0)) {
            nb_plot_set += 2;
          }
        }
      }
    }

    else {
      nb_plot_set++;
      if (marginal[i]) {
        nb_plot_set++;
      }
      if (hindex_parameter) {
        nb_plot_set++;
      }
    }
  }

  for (i = 0;i < marginal[0]->nb_value;i++) {
    if (self_transition[i]) {
      nb_plot_set += 2;
    }
  }

  plot_set = new MultiPlotSet(nb_plot_set , nb_variable);

  MultiPlotSet &plot = *plot_set;

  plot.border = "15 lw 0";

  for (i = 0;i < nb_variable;i++) {
    plot.variable_nb_viewpoint[i] = 0;
  }

  index = 0;
  if (self_transition) {
    plot.variable_nb_viewpoint[0]++;

    for (i = 0;i < marginal[0]->nb_value;i++) {
      if (self_transition[i]) {

        // vue : probabilites de rester dans l'etat i indexees

        if (nb_variable > 1) {
          title.str("");
          title << STAT_label[STATL_VARIABLE] << " " << 1;
          plot[index].title = title.str();
        }

        plot[index].xrange = Range(0 , self_transition[i]->length - 1);
        plot[index].yrange = Range(0. , 1.);

        if (self_transition[i]->length - 1 < TIC_THRESHOLD) {
          plot[index].xtics = 1;
        }

        plot[index].resize(1);

        legend.str("");
        legend << STAT_label[STATL_STATE] << " " << i << " - "
               << SEQ_label[SEQL_OBSERVED] << " " << SEQ_label[SEQL_SELF_TRANSITION];
        plot[index][0].legend = legend.str();

        plot[index][0].style = "linespoint";

        self_transition[i]->plotable_write(0 , plot[index][0]);
        index++;

        // vue : loi empirique des comptages de transition indexes

        if (nb_variable > 1) {
          title.str("");
          title << STAT_label[STATL_VARIABLE] << " " << 1;
          plot[index].title = title.str();
        }

        plot[index].xrange = Range(0 , self_transition[i]->length - 1);
        max_frequency = self_transition[i]->max_frequency_computation();
        plot[index].yrange = Range(0 , ceil(max_frequency * YSCALE));

        if (self_transition[i]->length - 1 < TIC_THRESHOLD) {
          plot[index].xtics = 1;
        }
        if (ceil(max_frequency * YSCALE) < TIC_THRESHOLD) {
          plot[index].ytics = 1;
        }

        plot[index].xlabel = SEQ_label[SEQL_INDEX];
        plot[index].ylabel = STAT_label[STATL_FREQUENCY];

        plot[index].resize(1);

        legend.str("");
        legend << STAT_label[STATL_STATE] << " " << i << " - "
               << SEQ_label[SEQL_TRANSITION_COUNTS];
        plot[index][0].legend = legend.str();

        plot[index][0].style = "impulses";

        self_transition[i]->plotable_frequency_write(plot[index][0]);
        index++;
      }
    }
  }

  for (i = 0;i < nb_variable;i++) {
    if (characteristics[i]) {
      characteristics[i]->plotable_write(plot , index , i , type[i] , *hlength);
    }
    else {
      plotable_write(plot , index , i);
    }
  }

  return plot_set;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet MarkovianSequences dans un fichier au format MTG.
 *
 *  arguments : reference sur un objet StatError, path,
 *              type de chaque variable (NUMERIC/SYMBOLIC).
 *
 *--------------------------------------------------------------*/

bool MarkovianSequences::mtg_write(StatError &error , const char *path , int *itype) const

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
            if (int_sequence[i][k][j] > 0) {
              out_file <<"\t\t+" << (char)('F' + int_sequence[i][k][j]) << 1 << endl;
            }
            break;
          }

          case NUMERIC : {
            for (m = 0;m < int_sequence[i][k][j];m++) {
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
