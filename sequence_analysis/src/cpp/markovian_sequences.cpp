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
 *       $Id: markovian_sequences.cpp 18058 2015-04-23 10:47:34Z guedon $
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
#include "stat_tool/curves.h"
#include "stat_tool/distribution.h"
#include "stat_tool/markovian.h"
#include "stat_tool/vectors.h"
#include "stat_tool/distance_matrix.h"
#include "stat_tool/stat_label.h"

#include "sequences.h"
#include "sequence_label.h"
#include "tool/config.h"

using namespace std;
using namespace stat_tool;


namespace sequence_analysis {



/*--------------------------------------------------------------*
 *
 *  Constructeur par defaut de la classe MarkovianSequences.
 *
 *--------------------------------------------------------------*/

MarkovianSequences::MarkovianSequences()

{
  min_interval = NULL;

  self_transition = NULL;
  observation_distribution = NULL;
  observation_histogram = NULL;
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


  min_interval = new double[nb_variable];
  for (i = 0;i < nb_variable;i++) {
    min_interval[i] = 0.;
  }

  self_transition = NULL;
  observation_distribution = NULL;
  observation_histogram = NULL;

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
  register int i;


  init();

  min_interval = new double[nb_variable];
  for (i = 0;i < nb_variable;i++) {
    min_interval_computation(i);
  }

  build_characteristic();
}


/*--------------------------------------------------------------*
 *
 *  Construction d'un objet MarkovianSequences avec ajout de variables auxilliaires.
 *
 *  argument : reference sur un objet MarkovianSequences,
 *             flags sur l'ajout des variables auxilliaires.
 *
 *--------------------------------------------------------------*/

MarkovianSequences::MarkovianSequences(const MarkovianSequences &seq , bool *auxiliary)
:Sequences(seq , auxiliary)

{
  register int i , j;


  min_interval = new double[nb_variable];

  self_transition = NULL;
  observation_distribution = NULL;
  observation_histogram = NULL;

  characteristics = new SequenceCharacteristics*[nb_variable];

  i = 0;
  for (j = 0;j < seq.nb_variable;j++) {
    min_interval[i] = seq.min_interval[j];

    if (seq.characteristics[j]) {
       characteristics[i] = new SequenceCharacteristics(*(seq.characteristics[j]));
    }
    else {
      characteristics[i] = NULL;
    }
    i++;

    if (auxiliary[j]) {
      min_interval[i] = 0.;
      characteristics[i] = NULL;
      i++;
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Copie d'un objet MarkovianSequences.
 *
 *  arguments : reference sur un objet MarkovianSequences,
 *              ajout/suppression des lois empiriques de temps de sejour initial.
 *
 *--------------------------------------------------------------*/

void MarkovianSequences::copy(const MarkovianSequences &seq , int param)

{
  bool initial_run_flag;
  register int i , j;


  min_interval = new double[nb_variable];
  for (i = 0;i < nb_variable;i++) {
    min_interval[i] = seq.min_interval[i];
  }

  if (seq.self_transition) {
    self_transition = new SelfTransition*[marginal_distribution[0]->nb_value];
    for (i = 0;i < marginal_distribution[0]->nb_value;i++) {
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

  if (seq.observation_distribution) {
    observation_distribution = new FrequencyDistribution**[nb_variable];
    observation_distribution[0] = NULL;

    for (i = 1;i < nb_variable;i++) {
      if (seq.observation_distribution[i]) {
        observation_distribution[i] = new FrequencyDistribution*[marginal_distribution[0]->nb_value];
        for (j = 0;j < marginal_distribution[0]->nb_value;j++) {
          observation_distribution[i][j] = new FrequencyDistribution(*(seq.observation_distribution[i][j]));
        }
      }

      else {
        observation_distribution[i] = NULL;
      }
    }
  }

  else {
    observation_distribution = NULL;
  }

  if (seq.observation_histogram) {
    observation_histogram = new Histogram**[nb_variable];
    observation_histogram[0] = NULL;

    for (i = 1;i < nb_variable;i++) {
      if (seq.observation_histogram[i]) {
        observation_histogram[i] = new Histogram*[marginal_distribution[0]->nb_value];
        for (j = 0;j < marginal_distribution[0]->nb_value;j++) {
          observation_histogram[i][j] = new Histogram(*(seq.observation_histogram[i][j]));
        }
      }

      else {
        observation_histogram[i] = NULL;
      }
    }
  }

  else {
    observation_histogram = NULL;
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
 *  Copie d'un objet MarkovianSequences avec inversion du sens de parcours des sequences.
 *
 *  arguments : reference sur un objet MarkovianSequences.
 *
 *--------------------------------------------------------------*/

void MarkovianSequences::reverse(const MarkovianSequences &seq)

{
  register int i , j;


  min_interval = new double[nb_variable];
  for (i = 0;i < nb_variable;i++) {
    min_interval[i] = seq.min_interval[i];
  }

  self_transition = NULL;

  if (seq.observation_distribution) {
    observation_distribution = new FrequencyDistribution**[nb_variable];
    observation_distribution[0] = NULL;

    for (i = 1;i < nb_variable;i++) {
      if (seq.observation_distribution[i]) {
        observation_distribution[i] = new FrequencyDistribution*[marginal_distribution[0]->nb_value];
        for (j = 0;j < marginal_distribution[0]->nb_value;j++) {
          observation_distribution[i][j] = new FrequencyDistribution(*(seq.observation_distribution[i][j]));
        }
      }

      else {
        observation_distribution[i] = NULL;
      }
    }
  }

  else {
    observation_distribution = NULL;
  }

  if (seq.observation_histogram) {
    observation_histogram = new Histogram**[nb_variable];
    observation_histogram[0] = NULL;

    for (i = 1;i < nb_variable;i++) {
      if (seq.observation_histogram[i]) {
        observation_histogram[i] = new Histogram*[marginal_distribution[0]->nb_value];
        for (j = 0;j < marginal_distribution[0]->nb_value;j++) {
          observation_histogram[i][j] = new Histogram(*(seq.observation_histogram[i][j]));
        }
      }

      else {
        observation_histogram[i] = NULL;
      }
    }
  }

  else {
    observation_histogram = NULL;
  }

  characteristics = new SequenceCharacteristics*[nb_variable];

  for (i = 0;i < nb_variable;i++) {
    if (seq.characteristics[i]) {
      characteristics[i] = new SequenceCharacteristics(*(seq.characteristics[i]), 'r');

      build_index_value(i);
      build_first_occurrence_frequency_distribution(i);

      if (!(seq.characteristics[i]->initial_run)) {
        build_sojourn_time_frequency_distribution(i);
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


  min_interval = new double[nb_variable];
  min_interval[0] = 0.;
  for (i = 0;i < seq.nb_variable;i++) {
    min_interval[i + 1] = seq.min_interval[i];
  }

  self_transition = NULL;
  observation_distribution = NULL;
  observation_histogram = NULL;

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
 *              'r' : inversion du sens de parcours, 'a' : ajout d'une variable d'etat,
 *              'm' : suppression du parametre d'index, 'e' : parametre d'index rendu explicite),
 *              ajout/suppression des lois empiriques de temps de sejour initial.
 *
 *--------------------------------------------------------------*/

MarkovianSequences::MarkovianSequences(const MarkovianSequences &seq , char transform , int param)

{
  switch (transform) {
  case 'c' :
    Sequences::copy(seq);
    copy(seq , param);
    break;
  case 'r' :
    Sequences::reverse(seq);
    reverse(seq);
    break;
  case 'a' :
    Sequences::add_state_variable(seq);
    add_state_variable(seq , param);
    break;
  case 'm' :
    Sequences::remove_index_parameter(seq);
    copy(seq);
    break;
  case 'e' :
    Sequences::explicit_index_parameter(seq);
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


  delete [] min_interval;

  if (self_transition) {
    for (i = 0;i < marginal_distribution[0]->nb_value;i++) {
      delete self_transition[i];
    }
    delete [] self_transition;
  }

  if (observation_distribution) {
    for (i = 1;i < nb_variable;i++) {
      if (observation_distribution[i]) {
        for (j = 0;j < marginal_distribution[0]->nb_value;j++) {
          delete observation_distribution[i][j];
        }
        delete [] observation_distribution[i];
      }
    }
    delete [] observation_distribution;
  }

  if (observation_histogram) {
    for (i = 1;i < nb_variable;i++) {
      if (observation_histogram[i]) {
        for (j = 0;j < marginal_distribution[0]->nb_value;j++) {
          delete observation_histogram[i][j];
        }
        delete [] observation_histogram[i];
      }
    }
    delete [] observation_histogram;
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
 *  argument : type de la 1ere variable (STATE / INT_VALUE / REAL_VALUE).
 *
 *--------------------------------------------------------------*/

void MarkovianSequences::state_variable_init(int itype)

{
  register int i , j;


  if (itype != type[0]) {
    if (type[0] == STATE) {
      if (self_transition) {
        for (i = 0;i < marginal_distribution[0]->nb_value;i++) {
          delete self_transition[i];
        }
        delete [] self_transition;

        self_transition = NULL;
      }

      if (observation_distribution) {
        for (i = 1;i < nb_variable;i++) {
          if (observation_distribution[i]) {
            for (j = 0;j < marginal_distribution[0]->nb_value;j++) {
              delete observation_distribution[i][j];
            }
            delete [] observation_distribution[i];
          }
        }
        delete [] observation_distribution;

        observation_distribution = NULL;
      }

      if (observation_histogram) {
        for (i = 1;i < nb_variable;i++) {
          if (observation_histogram[i]) {
            for (j = 0;j < marginal_distribution[0]->nb_value;j++) {
              delete observation_histogram[i][j];
            }
            delete [] observation_histogram[i];
          }
        }
        delete [] observation_histogram;

        observation_histogram = NULL;
      }
    }

    type[0] = itype;
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

    else if ((value < 0) || (value >= marginal_distribution[variable]->nb_value)) {
      status = false;
      ostringstream error_message;
      error_message << STAT_label[STATL_VALUE] << " " << value << " "
                    << STAT_error[STATR_NOT_PRESENT];
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
  register int i , j , k , m , n , p , q;
  int inb_sequence , cumul_nb_sequence , nb_histo , *ilength , *iidentifier ,
      **ivertex_identifier;
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

    // calcul du nombre de sequences

    inb_sequence = 0;
    for (i = 0;i < nb_sample;i++) {
      inb_sequence += pseq[i]->nb_sequence;
    }

    // comparaison des identificateurs des sequences

    iidentifier = new int[inb_sequence];

    cumul_nb_sequence = 0;
    i = 0;
    for (j = 0;j < nb_sample;j++) {
      for (k = 0;k < pseq[j]->nb_sequence;k++) {
        iidentifier[i] = pseq[j]->identifier[k];

        for (m = 0;m < cumul_nb_sequence;m++) {
          if (iidentifier[i] == iidentifier[m]) {
            delete [] iidentifier;
            iidentifier = NULL;
            break;
          }
        }

        if (!iidentifier) {
          break;
        }
        i++;
      }

      if (!iidentifier) {
        break;
      }
      cumul_nb_sequence += pseq[j]->nb_sequence;
    }

    // copie des longueurs des sequences

    ilength = new int[inb_sequence];

    i = 0;
    for (j = 0;j < nb_sample;j++) {
      for (k = 0;k < pseq[j]->nb_sequence;k++) {
        ilength[i++] = pseq[j]->length[k];
      }
    }

    // comparaison des identificateurs des vertex

    for (i = 0;i < nb_sample;i++) {
      if (!(pseq[i]->vertex_identifier)) {
        break;
      }
    }

    if (i == nb_sample) {
      ivertex_identifier = new int*[inb_sequence];

      cumul_nb_sequence = 0;
      i = 0;
      for (j = 0;j < nb_sample;j++) {
        for (k = 0;k < pseq[j]->nb_sequence;k++) {
          ivertex_identifier[i] = new int[pseq[j]->length[k]];
          for (m = 0;m < pseq[j]->length[k];m++) {
            ivertex_identifier[i][m] = pseq[j]->vertex_identifier[k][m];

            for (n = 0;n < cumul_nb_sequence;n++) {
              for (p = 0;p < ilength[n];p++) {
                if (ivertex_identifier[i][m] == ivertex_identifier[n][p]) {
                  for (q = 0;q <= i;q++) {
                    delete [] ivertex_identifier[q];
                  }
                  delete [] ivertex_identifier;
                  ivertex_identifier = NULL;
                  break;
                }
              }

              if (!ivertex_identifier) {
                break;
              }
            }

            if (!ivertex_identifier) {
              break;
            }
          }

          if (!ivertex_identifier) {
            break;
          }
          i++;
        }

        if (!ivertex_identifier) {
          break;
        }
        cumul_nb_sequence += pseq[j]->nb_sequence;
      }
    }

    else {
      ivertex_identifier = NULL;
    }

    seq = new MarkovianSequences(inb_sequence , iidentifier , ilength , ivertex_identifier ,
                                 index_parameter_type , nb_variable , type);
    delete [] iidentifier;
    delete [] ilength;

    if (ivertex_identifier) {
      for (i = 0;i < inb_sequence;i++) {
        delete [] ivertex_identifier[i];
      }
      delete [] ivertex_identifier;
    }

    phisto = new const FrequencyDistribution*[nb_sample];

    // copy of index parameters

    if (index_parameter_type == TIME) {
      i = 0;
      for (j = 0;j < nb_sample;j++) {
        for (k = 0;k < pseq[j]->nb_sequence;k++) {
          for (m = 0;m < pseq[j]->length[k];m++) {
            seq->index_parameter[i][m] = pseq[j]->index_parameter[k][m];
          }
          i++;
        }
      }

      for (i = 0;i < nb_sample;i++) {
        phisto[i] = pseq[i]->index_parameter_distribution;
      }
      seq->index_parameter_distribution = new FrequencyDistribution(nb_sample , phisto);

      for (i = 0;i < nb_sample;i++) {
        phisto[i] = pseq[i]->index_interval;
      }
      seq->index_interval = new FrequencyDistribution(nb_sample , phisto);
    }

    // copie des valeurs

    i = 0;
    for (j = 0;j < nb_sample;j++) {
      for (k = 0;k < pseq[j]->nb_sequence;k++) {
        for (m = 0;m < pseq[j]->nb_variable;m++) {
          if ((pseq[j]->type[m] != REAL_VALUE) && (pseq[j]->type[m] != AUXILIARY)) {
            for (n = 0;n < pseq[j]->length[k];n++) {
              seq->int_sequence[i][m][n] = pseq[j]->int_sequence[k][m][n];
            }
          }

          else {
            for (n = 0;n < pseq[j]->length[k];n++) {
              seq->real_sequence[i][m][n] = pseq[j]->real_sequence[k][m][n];
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

      if (seq->type[i] != AUXILIARY) {
        if (seq->type[i] != REAL_VALUE) {
          for (j = 0;j < nb_sample;j++) {
            phisto[j] = pseq[j]->marginal_distribution[i];
          }
          seq->marginal_distribution[i] = new FrequencyDistribution(nb_sample , phisto);
        }

        else {
          seq->build_marginal_histogram(i);
        }

        seq->min_interval[i] = pseq[0]->min_interval[i];
        for (j = 1;j < nb_sample;j++) {
          if (pseq[j]->min_interval[i] < seq->min_interval[i]) {
            seq->min_interval[i] = pseq[j]->min_interval[i];
          }
        }

        for (j = 0;j < nb_sample;j++) {
          if (!(pseq[j]->characteristics[i])) {
            break;
          }
        }

        if (j == nb_sample) {
          seq->characteristics[i] = new SequenceCharacteristics();

          seq->characteristics[i]->nb_value = seq->marginal_distribution[i]->nb_value;

          seq->build_index_value(i);

          seq->characteristics[i]->first_occurrence = new FrequencyDistribution*[seq->marginal_distribution[i]->nb_value];
          seq->characteristics[i]->recurrence_time = new FrequencyDistribution*[seq->marginal_distribution[i]->nb_value];

          for (j = 0;j < seq->marginal_distribution[i]->nb_value;j++) {
            nb_histo = 0;
            for (k = 0;k < nb_sample;k++) {
              if (j < pseq[k]->marginal_distribution[i]->nb_value) {
                phisto[nb_histo++] = pseq[k]->characteristics[i]->first_occurrence[j];
              }
            }
            seq->characteristics[i]->first_occurrence[j] = new FrequencyDistribution(nb_histo , phisto);

            nb_histo = 0;
            for (k = 0;k < nb_sample;k++) {
              if (j < pseq[k]->marginal_distribution[i]->nb_value) {
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
            seq->characteristics[i]->sojourn_time = new FrequencyDistribution*[seq->marginal_distribution[i]->nb_value];
            if (pseq[0]->characteristics[i]->initial_run) {
              seq->characteristics[i]->initial_run = new FrequencyDistribution*[seq->marginal_distribution[i]->nb_value];
            }
            seq->characteristics[i]->final_run = new FrequencyDistribution*[seq->marginal_distribution[i]->nb_value];

            for (j = 0;j < seq->marginal_distribution[i]->nb_value;j++) {
              nb_histo = 0;
              for (k = 0;k < nb_sample;k++) {
                if (j < pseq[k]->marginal_distribution[i]->nb_value) {
                  phisto[nb_histo++] = pseq[k]->characteristics[i]->sojourn_time[j];
                }
              }
              seq->characteristics[i]->sojourn_time[j] = new FrequencyDistribution(nb_histo , phisto);

              if (pseq[0]->characteristics[i]->initial_run) {
                nb_histo = 0;
                for (k = 0;k < nb_sample;k++) {
                  if (j < pseq[k]->marginal_distribution[i]->nb_value) {
                    phisto[nb_histo++] = pseq[k]->characteristics[i]->initial_run[j];
                  }
                }
                seq->characteristics[i]->initial_run[j] = new FrequencyDistribution(nb_histo , phisto);
              }

              nb_histo = 0;
              for (k = 0;k < nb_sample;k++) {
                if (j < pseq[k]->marginal_distribution[i]->nb_value) {
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
            seq->characteristics[i]->nb_run = new FrequencyDistribution*[seq->marginal_distribution[i]->nb_value];
            seq->characteristics[i]->nb_occurrence = new FrequencyDistribution*[seq->marginal_distribution[i]->nb_value];

            for (j = 0;j < seq->marginal_distribution[i]->nb_value;j++) {
              nb_histo = 0;
              for (k = 0;k < nb_sample;k++) {
                if (j < pseq[k]->marginal_distribution[i]->nb_value) {
                  phisto[nb_histo++] = pseq[k]->characteristics[i]->nb_run[j];
                }
              }
              seq->characteristics[i]->nb_run[j] = new FrequencyDistribution(nb_histo , phisto);

              nb_histo = 0;
              for (k = 0;k < nb_sample;k++) {
                if (j < pseq[k]->marginal_distribution[i]->nb_value) {
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

    if ((type[variable] == INT_VALUE) && (variable + 1 < nb_variable) &&
        (type[variable + 1] == AUXILIARY)) {
      status = false;
      ostringstream error_message;
      error_message << STAT_label[STATL_VARIABLE] << " " << variable + 1 << ": "
                    << STAT_error[STATR_VARIABLE_TYPE];
      error.update((error_message.str()).c_str());
    }
  }

  if (step < 1) {
    status = false;
    error.update(STAT_error[STATR_CLUSTERING_STEP]);
  }

  if (status) {
    seq = new MarkovianSequences(*this , variable , type[variable]);
    seq->Sequences::cluster(*this , variable , step , mode);

    for (i = 0;i < seq->nb_variable;i++) {
      if (i == variable) {
        seq->min_interval_computation(i);
        seq->build_characteristic(i , true , (((characteristics[i]) && (characteristics[i]->initial_run)) ? true : false));
      }

      else {
        seq->min_interval[i] = min_interval[i];
        if (characteristics[i]) {
          seq->characteristics[i] = new SequenceCharacteristics(*(characteristics[i]));
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

    if ((ivariable + 1 < nb_variable) && (type[ivariable + 1] == AUXILIARY)) {
      status = false;
      ostringstream error_message;
      error_message << STAT_label[STATL_VARIABLE] << " " << ivariable + 1 << ": "
                    << STAT_error[STATR_VARIABLE_TYPE];
      error.update((error_message.str()).c_str());
    }

    if (status) {
      min_symbol = marginal_distribution[ivariable]->nb_value;
      max_symbol = 0;

      for (i = 0;i < marginal_distribution[ivariable]->nb_value - marginal_distribution[ivariable]->offset;i++) {
        if ((symbol[i] < 0) || (symbol[i] >= (add_flag ? marginal_distribution[ivariable]->nb_value - 1 : marginal_distribution[ivariable]->nb_value))) {
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

      for (i = 0;i < marginal_distribution[ivariable]->nb_value - marginal_distribution[ivariable]->offset;i++) {
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

      seq = new MarkovianSequences(nb_sequence , identifier , length , vertex_identifier ,
                                   index_parameter_type , nb_variable + offset , itype);
      delete [] itype;

      seq->Sequences::transcode(*this , ivariable , 0 , max_symbol , symbol , add_flag);

      for (i = 0;i < seq->nb_variable;i++) {
        if (i == variable) {
          seq->min_interval_computation(i);
          seq->build_characteristic(i , true , (((characteristics[ivariable]) && (characteristics[ivariable]->initial_run)) ? true : false));
        }

        else {
          seq->min_interval[i] = min_interval[i - offset];
          if (characteristics[i - offset]) {
            seq->characteristics[i] = new SequenceCharacteristics(*(characteristics[i - offset]));
          }
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
 *              sur un objet CategoricalSequenceProcess.
 *
 *--------------------------------------------------------------*/

MarkovianSequences* MarkovianSequences::transcode(StatError &error ,
                                                  const CategoricalSequenceProcess *process) const

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
      for (i = 0;i < marginal_distribution[ivariable]->nb_value;i++) {
        if (marginal_distribution[ivariable]->frequency[i] == 0) {
          break;
        }
      }

      if (i == marginal_distribution[ivariable]->nb_value) {
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
    for (i = 0;i < marginal_distribution[ivariable]->nb_value;i++) {
      if (marginal_distribution[ivariable]->frequency[i] == 0) {
        os << " " << i;
      }
    }
    os << endl;
#   endif

    symbol = new int[marginal_distribution[ivariable]->nb_value - marginal_distribution[ivariable]->offset];

//    i = 0;
    i = -1;
    for (j = marginal_distribution[ivariable]->offset;j < marginal_distribution[ivariable]->nb_value;j++) {
//      symbol[j - marginal_distribution[ivariable]->offset] = i;
      if (marginal_distribution[ivariable]->frequency[j] > 0) {
        i++;
      }
      symbol[j - marginal_distribution[ivariable]->offset] = i;
    }
//    max = i - 1;
    max = i;

#   ifdef DEBUG
    cout << "\nTest :";
    for (i = 0;i < marginal_distribution[ivariable]->nb_value - marginal_distribution[ivariable]->offset;i++) {
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

    seq = new MarkovianSequences(nb_sequence , identifier , length , vertex_identifier ,
                                 index_parameter_type , nb_variable + offset , itype);
    delete [] itype;

    seq->Sequences::transcode(*this , ivariable , 0 , max , symbol , add_flag);
    delete [] symbol;

    for (i = 0;i < seq->nb_variable;i++) {
      if (i == variable) {
        seq->min_interval_computation(i);
        seq->build_characteristic(i , true , (((characteristics[ivariable]) && (characteristics[ivariable]->initial_run)) ? true : false));
      }

      else if (characteristics[i - offset]) {
        seq->min_interval[i] = min_interval[i - offset];
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

    else if ((nb_class < 2) || (nb_class >= marginal_distribution[ivariable]->nb_value)) {
      status = false;
      error.update(STAT_error[STATR_NB_CLASS]);
    }

    if ((ivariable + 1 < nb_variable) && (type[ivariable + 1] == AUXILIARY)) {
      status = false;
      ostringstream error_message;
      error_message << STAT_label[STATL_VARIABLE] << " " << ivariable + 1 << ": "
                    << STAT_error[STATR_VARIABLE_TYPE];
      error.update((error_message.str()).c_str());
    }
  }

  if (status) {
    limit = new int[nb_class + 1];
    limit[0] = marginal_distribution[ivariable]->offset;
    for (i = 1;i < nb_class;i++) {
      limit[i] = ilimit[i - 1];
    }
    limit[nb_class] = marginal_distribution[ivariable]->nb_value;

    for (i = 1;i <= nb_class;i++) {
      if (limit[i] <= limit[i - 1]) {
        status = false;
        error.update(STAT_error[STATR_CLUSTER_LIMIT]);
      }
    }

    if (status) {
      symbol = new int[marginal_distribution[ivariable]->nb_value];

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

      seq = new MarkovianSequences(nb_sequence , identifier , length , vertex_identifier ,
                                   index_parameter_type , nb_variable + offset , itype);
      delete [] itype;

      seq->Sequences::transcode(*this , ivariable , 0 , nb_class - 1 , symbol , add_flag);
      delete [] symbol;

      for (i = 0;i < seq->nb_variable;i++) {
        if (i == variable) {
          seq->min_interval_computation(i);
          seq->build_characteristic(i , true , (((characteristics[ivariable]) && (characteristics[ivariable]->initial_run)) ? true : false));
        }
        else {
          seq->min_interval[i] = min_interval[i - offset];
          if (characteristics[i - offset]) {
            seq->characteristics[i] = new SequenceCharacteristics(*(characteristics[i - offset]));
          }
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

    if ((variable + 1 < nb_variable) && (type[variable + 1] == AUXILIARY)) {
      status = false;
      ostringstream error_message;
      error_message << STAT_label[STATL_VARIABLE] << " " << variable + 1 << ": "
                    << STAT_error[STATR_VARIABLE_TYPE];
      error.update((error_message.str()).c_str());
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
      seq = new MarkovianSequences(*this , variable , type[variable]);
      seq->Sequences::cluster(*this , variable , nb_class , limit);

      for (i = 0;i < seq->nb_variable;i++) {
        if (i == variable) {
          seq->min_interval_computation(i);
          seq->build_characteristic(i);
        }

        else {
          seq->min_interval[i] = min_interval[i];
          if (characteristics[i]) {
            seq->characteristics[i] = new SequenceCharacteristics(*(characteristics[i]));
          }
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
    seq = new MarkovianSequences(*this , 'm');
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Copie d'un objet MarkovianSequences avec transformation du parametre d'index implicite
 *  en parametre d'index explicite.
 *
 *  argument : reference sur un objet StatError.
 *
 *--------------------------------------------------------------*/

MarkovianSequences* MarkovianSequences::explicit_index_parameter(StatError &error) const

{
  MarkovianSequences *seq;


  error.init();

  if (index_parameter) {
    seq = NULL;
    error.update(SEQ_error[SEQR_INDEX_PARAMETER_TYPE]);
  }
  else {
    seq = new MarkovianSequences(*this , 'e');
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

    for (i = 0;i < bnb_variable;i++) {
      if ((type[variable[i]] == AUXILIARY) &&
          ((i == 0) || (variable[i - 1] != variable[i] - 1))) {
        status = false;
        error.update(SEQ_error[SEQR_VARIABLE_INDICES]);
      }
    }

    if (status) {
      itype = new int[bnb_variable];
      for (i = 0;i < bnb_variable;i++) {
        itype[i] = type[variable[i]];
      }

      seq = new MarkovianSequences(nb_sequence , identifier , length , vertex_identifier ,
                                   index_parameter_type , bnb_variable , itype);

      seq->Sequences::select_variable(*this , variable);

      for (i = 0;i < seq->nb_variable;i++) {
        seq->min_interval[i] = min_interval[variable[i]];
        if (characteristics[variable[i]]) {
          seq->characteristics[i] = new SequenceCharacteristics(*(characteristics[variable[i]]));
        }
      }

      delete [] itype;
    }

    delete [] variable;
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

  seq = new MarkovianSequences(nb_sequence , identifier , length , vertex_identifier ,
                               index_parameter_type , nb_variable - 1 , itype);

  seq->Sequences::select_variable(*this , variable);

  for (i = 0;i < seq->nb_variable;i++) {
    seq->min_interval[i] = min_interval[i + 1];
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
 *              pointeurs sur les objets MarkovianSequences,
 *              echantillon de reference pour les identificateurs.
 *
 *--------------------------------------------------------------*/

MarkovianSequences* MarkovianSequences::merge_variable(StatError &error , int nb_sample ,
                                                       const MarkovianSequences **iseq , int ref_sample) const

{
  bool status = true;
  register int i , j , k , m;
  int inb_variable , *iidentifier , *itype , **ivertex_identifier;
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

    // comparaison des identificateurs des sequences

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
        iidentifier = NULL;
      }
      else {
        iidentifier = pseq[0]->identifier;
      }
    }

    else {
      ref_sample--;
      iidentifier = pseq[ref_sample]->identifier;
    }

    // comparaison des identificateurs des vertex

    if (ref_sample == I_DEFAULT) {
      for (i = 0;i < nb_sample;i++) {
        if (!(pseq[i]->vertex_identifier)) {
          break;
        }
      }

      if (i == nb_sample) {
        for (i = 0;i < nb_sequence;i++) {
          for (j = 1;j < nb_sample;j++) {
            for (k = 0;k < pseq[j]->length[i];k++) {
              if (pseq[j]->vertex_identifier[i][k] != pseq[0]->vertex_identifier[i][k]) {
                break;
              }
            }

            if (k < pseq[j]->length[i]) {
              break;
            }
          }

          if (j < nb_sample) {
            break;
          }
        }

        if (i < nb_sequence) {
          ivertex_identifier = NULL;
        }
        else {
          ivertex_identifier = pseq[0]->vertex_identifier;
        }
      }

      else {
        ivertex_identifier = NULL;
      }
    }

    else {
      ivertex_identifier = pseq[ref_sample]->vertex_identifier;
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

    seq = new MarkovianSequences(nb_sequence , iidentifier , length , ivertex_identifier ,
                                 index_parameter_type , inb_variable , itype);
    delete [] itype;

    // copy of index parameters

    if (index_parameter_distribution) {
      seq->index_parameter_distribution = new FrequencyDistribution(*index_parameter_distribution);
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

    // copie des valeurs

    for (i = 0;i < nb_sequence;i++) {
      inb_variable = 0;
      for (j = 0;j < nb_sample;j++) {
        for (k = 0;k < pseq[j]->nb_variable;k++) {
          if ((seq->type[inb_variable] != REAL_VALUE) && (seq->type[inb_variable] != AUXILIARY)) {
            for (m = 0;m < length[i];m++) {
              seq->int_sequence[i][inb_variable][m] = pseq[j]->int_sequence[i][k][m];
            }
          }

          else {
            for (m = 0;m < length[i];m++) {
              seq->real_sequence[i][inb_variable][m] = pseq[j]->real_sequence[i][k][m];
            }
          }

          inb_variable++;
        }
      }
    }

    inb_variable = 0;
    for (i = 0;i < nb_sample;i++) {
      for (j = 0;j < pseq[i]->nb_variable;j++) {
        seq->min_value[inb_variable] = pseq[i]->min_value[j];
        seq->max_value[inb_variable] = pseq[i]->max_value[j];

        if (pseq[i]->marginal_distribution[j]) {
          seq->marginal_distribution[inb_variable] = new FrequencyDistribution(*(pseq[i]->marginal_distribution[j]));
        }
        if (pseq[i]->marginal_histogram[j]) {
          seq->marginal_histogram[inb_variable] = new Histogram(*(pseq[i]->marginal_histogram[j]));
        }

        seq->min_interval[inb_variable] = pseq[i]->min_interval[j];

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
  int end_value , *ilength;
  double mean , variance , *deviation;
  MarkovianSequences *seq;


  seq = NULL;
  error.init();

/*  if (index_parameter_type == TIME) {
    status = false;
    error.update(SEQ_error[SEQR_INDEX_PARAMETER_TYPE]);
  } */

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

    seq = new MarkovianSequences(nb_sequence , identifier , ilength , vertex_identifier ,
                                 index_parameter_type , nb_variable , type , false);
    delete [] ilength;

    // copy of vertex identifiers

    if (vertex_identifier) {
      for (i = 0;i < seq->nb_sequence;i++) {
        for (j = 0;j < length[i];j++) {
          seq->vertex_identifier[i][j] = vertex_identifier[i][j];
        }
        for (j = length[i];j < seq->length[i];j++) {
          seq->vertex_identifier[i][j] = I_DEFAULT;
        }
      }
    }

    // copy of index parameters

    if (index_parameter) {
      for (i = 0;i < seq->nb_sequence;i++) {
        for (j = 0;j < length[i];j++) {
          seq->index_parameter[i][j] = index_parameter[i][j];
        }
        for (j = length[i];j < seq->length[i];j++) {
          seq->index_parameter[i][j] = seq->index_parameter[i][j - 1] + 1;
        }
      }

      seq->build_index_parameter_frequency_distribution();
      if (index_interval) {
        seq->index_interval_computation();
      }
    }

    deviation = new double[nb_variable];

    for (i = 0;i < nb_variable;i++) {
      if ((type[i] == REAL_VALUE) || (type[i] == AUXILIARY)) {
        mean = mean_computation(i);
        variance = variance_computation(i , mean);
        deviation[i] = sqrt(variance) / ABSORBING_RUN_STANDARD_DEVIATION_FACTOR;
      }
    }

    // copie des sequences avec ajout series finales absorbantes

    for (i = 0;i < seq->nb_sequence;i++) {
      for (j = 0;j < seq->nb_variable;j++) {
        if ((seq->type[j] != REAL_VALUE) && (seq->type[j] != AUXILIARY)) {
          for (k = 0;k < length[i];k++) {
            seq->int_sequence[i][j][k] = int_sequence[i][j][k];
          }

          if (min_value[j] > 0) {
            end_value = 0;
          }
          else {
            end_value = max_value[j] + 1;
          }

          for (k = length[i];k < seq->length[i];k++) {
            seq->int_sequence[i][j][k] = end_value;
          }

#         ifdef DEBUG
          if (run_length == 1) {  // pour UCs Fuji/Braeburn
            if ((j == 0) || (j == 1)) {
              seq->int_sequence[i][j][seq->length[i] - 1] = seq->int_sequence[i][j][seq->length[i] - 2];
            }
            else if (j > 2) {
              seq->int_sequence[i][j][seq->length[i] - 1] = 0;
            }
          }
#         endif

        }

        else {
          for (k = 0;k < length[i];k++) {
            seq->real_sequence[i][j][k] = real_sequence[i][j][k];
          }

          if (min_value[j] >= 5 * deviation[j]) {
            for (k = length[i];k < seq->length[i];k++) {
              seq->real_sequence[i][j][k] = min_value[j] - (k % 2 + 4) * deviation[j];
            }
          }

          else {
            for (k = length[i];k < seq->length[i];k++) {
              seq->real_sequence[i][j][k] = max_value[j] + (k % 2 + 4) * deviation[j];
            }
          }
        }
      }
    }

    delete [] deviation;

    for (i = 0;i < seq->nb_variable;i++) {
      seq->min_value_computation(i);
      seq->max_value_computation(i);

      seq->build_marginal_frequency_distribution(i);

      seq->min_interval_computation(i);
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
 *  Construction des variables auxilliaires correspondant a
 *  la restauration des sequences.
 *
 *  arguments : pointeurs sur des objets DiscreteParametricProcess et
 *              sur des objets ContinuousParametricProcess.
 *
 *--------------------------------------------------------------*/

MarkovianSequences* MarkovianSequences::build_auxiliary_variable(DiscreteParametricProcess **discrete_process ,
                                                                 ContinuousParametricProcess **continuous_process) const

{
  bool *auxiliary;
  register int i , j , k , m;
  int *pstate;
  double *pauxiliary , *mean;
  MarkovianSequences *seq;


  auxiliary = new bool[nb_variable];

  auxiliary[0] = false;
  for (i = 1;i < nb_variable;i++) {
    if (((discrete_process) && (discrete_process[i - 1])) ||
        ((continuous_process) && (continuous_process[i - 1]))) {
      auxiliary[i] = true;
    }
    else {
      auxiliary[i] = false;
    }
  }

  seq = new MarkovianSequences(*this , auxiliary);

  i = 0;
  for (j = 1;j < nb_variable;j++) {
    i++;

    if ((discrete_process) && (discrete_process[j - 1])) {
      i++;
      for (k = 0;k < nb_sequence;k++) {
        pstate = seq->int_sequence[k][0];
        pauxiliary = seq->real_sequence[k][i];

        for (m = 0;m < length[k];m++) {
          *pauxiliary++ = discrete_process[j - 1]->observation[*pstate++]->mean;
        }
      }
    }

    else if ((continuous_process) && (continuous_process[j - 1])) {
      i++;

      if ((continuous_process[j - 1]->ident == GAMMA) || (continuous_process[j - 1]->ident == ZERO_INFLATED_GAMMA)) {
        mean = new double [continuous_process[j - 1]->nb_state];

        switch (continuous_process[j - 1]->ident) {

        case GAMMA : {
          for (k = 0;k < continuous_process[j - 1]->nb_state;k++) {
            mean[k] = continuous_process[j - 1]->observation[k]->shape *
                      continuous_process[j - 1]->observation[k]->scale;
          }
          break;
        }

        case ZERO_INFLATED_GAMMA : {
          for (k = 0;k < continuous_process[j - 1]->nb_state;k++) {
            if (continuous_process[j - 1]->observation[k]->zero_probability == 1.) {
              mean[k] = 0.;
            }
            else {
              mean[k] = (1 - continuous_process[j - 1]->observation[k]->zero_probability) *
                        continuous_process[j - 1]->observation[k]->shape *
                        continuous_process[j - 1]->observation[k]->scale;
            }
          }
          break;
        }
        }

        for (k = 0;k < nb_sequence;k++) {
          pstate = seq->int_sequence[k][0];
          pauxiliary = seq->real_sequence[k][i];

          for (m = 0;m < length[k];m++) {
            *pauxiliary++ = mean[*pstate++];
          }
        }

        delete [] mean;
      }

      else if ((continuous_process[j - 1]->ident == GAUSSIAN) || (continuous_process[j - 1]->ident == VON_MISES)) {
        for (k = 0;k < nb_sequence;k++) {
          pstate = seq->int_sequence[k][0];
          pauxiliary = seq->real_sequence[k][i];

          for (m = 0;m < length[k];m++) {
            *pauxiliary++ = continuous_process[j - 1]->observation[*pstate++]->location;
          }
        }
      }

      else if (continuous_process[j - 1]->ident == LINEAR_MODEL) {
        switch (index_parameter_type) {

        case IMPLICIT_TYPE : {
          for (k = 0;k < nb_sequence;k++) {
            pstate = seq->int_sequence[k][0];
            pauxiliary = seq->real_sequence[k][i];

            for (m = 0;m < length[k];m++) {
              *pauxiliary++ = continuous_process[j - 1]->observation[*pstate]->intercept +
                              continuous_process[j - 1]->observation[*pstate]->slope * m;
              pstate++;
            }
          }
          break;
        }

        case TIME : {
          for (k = 0;k < nb_sequence;k++) {
            pstate = seq->int_sequence[k][0];
            pauxiliary = seq->real_sequence[k][i];

            for (m = 0;m < length[k];m++) {
              *pauxiliary++ = continuous_process[j - 1]->observation[*pstate]->intercept +
                              continuous_process[j - 1]->observation[*pstate]->slope * index_parameter[k][m];
              pstate++;
            }
          }
          break;
        }
        }
      }
    }
  }

  for (i = 1;i < seq->nb_variable;i++) {
    if (seq->type[i] == AUXILIARY) {
      seq->min_value_computation(i);
      seq->max_value_computation(i);
    }
  }

  delete [] auxiliary;

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

    seq = new MarkovianSequences(inb_sequence , NULL , ilength , vertex_identifier ,
                                 index_parameter_type , nb_variable , type);
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
      seq->index_parameter_distribution = new FrequencyDistribution(*index_parameter_distribution);
      seq->index_interval_computation();
    }

    for (i = 0;i < seq->nb_variable;i++) {
      seq->min_value[i] = min_value[i];
      seq->max_value[i] = max_value[i];

      if (marginal_distribution[i]) {
        seq->marginal_distribution[i] = new FrequencyDistribution(*marginal_distribution[i]);
      }
      if (marginal_histogram[i]) {
        seq->marginal_histogram[i] = new Histogram(*marginal_histogram[i]);
      }

      seq->min_interval_computation(i);
    }

    seq->build_characteristic();
  }

  return seq;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la fonction de repartition empirique pour une variable.
 *
 *  arguments : indice de la variable, (valeurs, fonction de repartition).
 *
 *--------------------------------------------------------------*/

int MarkovianSequences::cumulative_distribution_function_computation(int variable , double **cdf) const

{
  register int i , j , k;
  int cumul , int_min , int_value , frequency;
  double real_min , real_value;


  if (marginal_distribution[variable]) {
    i = marginal_distribution[variable]->cumulative_distribution_function_computation(cdf);
  }

  else {
    cdf[0] = new double[cumul_length];
    cdf[1] = new double[cumul_length];

    cumul = 0;
    i = 0;

    switch (type[variable]) {

    case INT_VALUE : {
      do {

        // recherche de la valeur minimum courante

        if (cumul == 0) {
          int_value = (int)min_value[variable];
        }

        else {
          int_min = (int)max_value[variable] + 1;
          for (j = 0;j < nb_sequence;j++) {
            for (k = 0;k < length[j];k++) {
              if ((int_sequence[j][variable][k] > int_value) &&
                  (int_sequence[j][variable][k] < int_min)) {
                int_min = int_sequence[j][variable][k];
              }
            }
          }
          int_value = int_min;
        }

        // recherche du nombre de vecteurs prenant pour la variable selectionnee
        // la valeur minimum courante

        frequency = 0;
        for (j = 0;j < nb_sequence;j++) {
          for (k = 0;k < length[j];k++) {
            if (int_sequence[j][variable][k] == int_value) {
              frequency++;
            }
          }
        }

        cdf[0][i] = int_value;
        cdf[1][i] = (cumul + (double)(frequency + 1) / 2.) / (double)cumul_length;
        cumul += frequency;
        i++;
      }
      while (cumul < cumul_length);
      break;
    }

    case REAL_VALUE : {
      do {

        // recherche de la valeur minimum courante

        if (cumul == 0) {
          real_value = min_value[variable];
        }

        else {
          real_min = max_value[variable] + 1;
          for (j = 0;j < nb_sequence;j++) {
            for (k = 0;k < length[j];k++) {
              if ((real_sequence[j][variable][k] > real_value) &&
                  (real_sequence[j][variable][k] < real_min)) {
                real_min = real_sequence[j][variable][k];
              }
            }
          }
          real_value = real_min;
        }

        // recherche du nombre de vecteurs prenant pour la variable selectionnee
        // la valeur minimum courante

        frequency = 0;
        for (j = 0;j < nb_sequence;j++) {
          for (k = 0;k < length[j];k++) {
            if (real_sequence[j][variable][k] == real_value) {
              frequency++;
            }
          }
        }

        cdf[0][i] = real_value;
        cdf[1][i] = (cumul + (double)(frequency + 1) / 2.) / (double)cumul_length;
        cumul += frequency;
        i++;
      }
      while (cumul < cumul_length);
      break;
    }
    }
  }

# ifdef DEBUG
  cout << "\nCumul: ";
  for (j = 0;j < i;j++) {
    cout << cdf[0][j] << " " << cdf[1][j] << " | ";
  }
  cout << endl;
# endif

  return i;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la fonction de repartition empirique pour un etat et une variable.
 *
 *  arguments : indice de la variable, etat, (valeurs, fonction de repartition).
 *
 *--------------------------------------------------------------*/

int MarkovianSequences::cumulative_distribution_function_computation(int variable , int state ,
                                                                     double **cdf) const

{
  register int i , j , k;
  int cumul , int_min , int_value , frequency;
  double real_min , real_value;


  if (observation_distribution[variable]) {
    i = observation_distribution[variable][state]->cumulative_distribution_function_computation(cdf);
  }

  else {
    cdf[0] = new double[marginal_distribution[0]->frequency[state]];
    cdf[1] = new double[marginal_distribution[0]->frequency[state]];

    cumul = 0;
    i = 0;

    switch (type[variable]) {

    case INT_VALUE : {
      do {

        // recherche de la valeur minimum courante

        if (cumul == 0) {
          int_value = (int)min_value[variable];
        }

        else {
          int_min = (int)max_value[variable] + 1;
          for (j = 0;j < nb_sequence;j++) {
            for (k = 0;k < length[j];k++) {
              if ((int_sequence[j][0][k] == state) &&
                  (int_sequence[j][variable][k] > int_value) &&
                  (int_sequence[j][variable][k] < int_min)) {
                int_min = int_sequence[j][variable][k];
              }
            }
          }
          int_value = int_min;
        }

        // recherche du nombre de vecteurs prenant pour la variable selectionnee
        // la valeur minimum courante

        frequency = 0;
        for (j = 0;j < nb_sequence;j++) {
          for (k = 0;k < length[j];k++) {
            if ((int_sequence[j][0][k] == state) &&
                (int_sequence[j][variable][k] == int_value)) {
              frequency++;
            }
          }
        }

        cdf[0][i] = int_value;
        cdf[1][i] = (cumul + (double)(frequency + 1) / 2.) /
                    (double)marginal_distribution[0]->frequency[state];
        cumul += frequency;
        i++;
      }
      while (cumul < marginal_distribution[0]->frequency[state]);
      break;
    }

    case REAL_VALUE : {
      do {

        // recherche de la valeur minimum courante

        if (cumul == 0) {
          real_value = min_value[variable];
        }

        else {
          real_min = max_value[variable] + 1;
          for (j = 0;j < nb_sequence;j++) {
            for (k = 0;k < length[j];k++) {
              if ((int_sequence[j][0][k] == state) &&
                  (real_sequence[j][variable][k] > real_value) &&
                  (real_sequence[j][variable][k] < real_min)) {
                real_min = real_sequence[j][variable][k];
              }
            }
          }
          real_value = real_min;
        }

        // recherche du nombre de vecteurs prenant pour la variable selectionnee
        // la valeur minimum courante

        frequency = 0;
        for (j = 0;j < nb_sequence;j++) {
          for (k = 0;k < length[j];k++) {
            if ((int_sequence[j][0][k] == state) &&
                (real_sequence[j][variable][k] == real_value)) {
              frequency++;
            }
          }
        }

        cdf[0][i] = real_value;
        cdf[1][i] = (cumul + (double)(frequency + 1) / 2.) /
                    (double)marginal_distribution[0]->frequency[state];
        cumul += frequency;
        i++;
      }
      while (cumul < marginal_distribution[0]->frequency[state]);
      break;
    }
    }
  }

  return i;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de l'intervalle minimum entre 2 valeurs pour une variable.
 *
 *  argument : indice de la variable.
 *
 *--------------------------------------------------------------*/

void MarkovianSequences::min_interval_computation(int variable)

{
  if (marginal_distribution[variable]) {
    min_interval[variable] = marginal_distribution[variable]->min_interval_computation();
  }

  else if (type[variable] != AUXILIARY) {
    bool *selected_value;
    register int i , j , k , m;
    int int_min , int_value , nb_value , max_frequency , *frequency , *index;
    double real_min , real_value;


    min_interval[variable] = max_value[variable] - min_value[variable];
    i = 0;

    switch (type[variable]) {

    case INT_VALUE : {
      do {

        // recherche de la valeur minimum courante

        if (i == 0) {
          int_value = (int)min_value[variable];
        }

        else {
          int_min = (int)max_value[variable] + 1;
          for (j = 0;j < nb_sequence;j++) {
            for (k = 0;k < length[j];k++) {
              if ((int_sequence[j][variable][k] > int_value) &&
                  (int_sequence[j][variable][k] < int_min)) {
                int_min = int_sequence[j][variable][k];
              }
            }
          }

          if (int_min - int_value < min_interval[variable]) {
            min_interval[variable] = int_min - int_value;
          }
          int_value = int_min;
        }

        // recherche du nombre de vecteurs prenant pour la variable selectionnee
        // la valeur minimum courante

        for (j = 0;j < nb_sequence;j++) {
          for (k = 0;k < length[j];k++) {
            if (int_sequence[j][variable][k] == int_value) {
              i++;
            }
          }
        }
      }
      while (i < cumul_length);
      break;
    }

    case REAL_VALUE : {
//      double max_interval = 0.;

      frequency = new int[cumul_length];

      j = 0;

      do {

        // recherche de la valeur minimum courante

        if (i == 0) {
          real_value = min_value[variable];
        }

        else {
          real_min = max_value[variable] + 1;
          for (k = 0;k < nb_sequence;k++) {
            for (m = 0;m < length[k];m++) {
              if ((real_sequence[k][variable][m] > real_value) &&
                  (real_sequence[k][variable][m] < real_min)) {
                real_min = real_sequence[k][variable][m];
              }
            }
          }

          if (real_min - real_value < min_interval[variable]) {
            min_interval[variable] = real_min - real_value;
          }
/*          if (real_min - real_value > max_interval) {
            max_interval = real_min - real_value;
          } */
          real_value = real_min;
        }

        // recherche du nombre de vecteurs prenant pour la variable selectionnee
        // la valeur minimum courante

        frequency[j] = 0;
        for (k = 0;k < nb_sequence;k++) {
          for (m = 0;m < length[k];m++) {
            if (real_sequence[k][variable][m] == real_value) {
              i++;
              frequency[j]++;
            }
          }
        }
        j++;
      }
      while (i < cumul_length);

      // recherche de la frequence mediane

      nb_value = j;
      selected_value = new bool[nb_value];
      index = new int [nb_value];

      for (i = 0;i < nb_value;i++) {
        selected_value[i] = false;
      }

      i = 0;
      do {
        max_frequency = 0;
        for (j = 0;j < nb_value;j++) {
          if ((!selected_value[j]) && (frequency[j] > max_frequency)) {
            max_frequency = frequency[j];
          }
        }

        for (j = 0;j < nb_value;j++) {
          if (frequency[j] == max_frequency) {
            selected_value[j] = true;
            index[i++] = j;
          }
        }
      }
      while (i < nb_value);

#     ifdef MESSAGE
      cout << "\n" << STAT_label[STATL_VARIABLE] << " " <<  variable + 1 << ": "
           << min_interval[variable] << " " << frequency[index[nb_value / 2]];
#     endif

      // a finaliser

      if (frequency[index[nb_value / 2]] == 1) {
        min_interval[variable] = 0.;
      }

#     ifdef MESSAGE
      cout << " | " << nb_value << " " << cumul_length << endl;
#     endif

      delete [] frequency;
      delete [] selected_value;
      delete [] index;
      break;
    }
    }
  }
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
    if (marginal_distribution[i]) {
      information += marginal_distribution[i]->information_computation();
    }
    else {
      information = D_INF;
      break;
    }
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
  int num , denom;


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

    self_transition[state]->frequency[i] = denom;
    if (denom > 0) {
      self_transition[state]->point[0][i] = (double)num / (double)denom;
    }
    else {
      self_transition[state]->point[0][i] = D_DEFAULT;
    }
  }

# ifdef DEBUG
  double sum = 0.;

  for (i = 0;i < max_length - 1;i++) {
    sum += self_transition[state]->frequency[i] * self_transition[state]->point[0][i];
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
    self_transition = new SelfTransition*[marginal_distribution[0]->nb_value];

    for (i = 0;i < marginal_distribution[0]->nb_value;i++) {
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
    self_transition = new SelfTransition*[marginal_distribution[0]->nb_value];

    for (i = 0;i < marginal_distribution[0]->nb_value;i++) {
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
 *  Calcul de la distribution marginale des etats a partir de la restauration.
 *
 *--------------------------------------------------------------*/

Distribution* MarkovianSequences::weight_computation() const

{
  register int i;
  Distribution *weight;


  if (type[0] == STATE) {
    weight = new Distribution(marginal_distribution[0]->nb_value);

    for (i = 0;i < marginal_distribution[0]->nb_value;i++) {
      weight->mass[i] = (double)marginal_distribution[0]->frequency[i] /
                        (double)marginal_distribution[0]->nb_element;
    }

    weight->cumul_computation();
    weight->max = (double)marginal_distribution[0]->max /
                  (double)marginal_distribution[0]->nb_element;
  }

  else {
    weight = NULL;
  }

  return weight;
}


/*--------------------------------------------------------------*
 *
 *  Accumulation des observations (pour une variable donnee).
 *
 *  arguments : indice de la variable, nombre d'etats.
 *
 *--------------------------------------------------------------*/

void MarkovianSequences::observation_frequency_distribution_computation(int variable ,
                                                                        int nb_state)

{
  register int i , j;
  int *pstate , *poutput;


  // initialisation des lois empiriques

  for (i = 0;i < nb_state;i++) {
    for (j = 0;j < marginal_distribution[variable]->nb_value;j++) {
      observation_distribution[variable][i]->frequency[j] = 0;
    }
  }

  // mise a jour des lois empiriques

  for (i = 0;i < nb_sequence;i++) {
    pstate = int_sequence[i][0];
    poutput = int_sequence[i][variable];
    for (j = 0;j < length[i];j++) {
      (observation_distribution[variable][*pstate++]->frequency[*poutput++])++;
    }
  }

  // extraction des caracteristiques des lois empiriques

  for (i = 0;i < nb_state;i++) {
    if (!characteristics[variable]) {
      observation_distribution[variable][i]->nb_value_computation();
    }
    observation_distribution[variable][i]->offset_computation();
    observation_distribution[variable][i]->nb_element_computation();
    observation_distribution[variable][i]->max_computation();

    if (!characteristics[variable]) {
      observation_distribution[variable][i]->mean_computation();
      observation_distribution[variable][i]->variance_computation();
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Construction des lois empiriques d'observation.
 *
 *  argument : nombre d'etats.
 *
 *--------------------------------------------------------------*/

void MarkovianSequences::build_observation_frequency_distribution(int nb_state)

{
  if ((nb_variable > 1) && (!observation_distribution)) {
    register int i , j;


    observation_distribution = new FrequencyDistribution**[nb_variable];
    observation_distribution[0] = NULL;

    for (i = 1;i < nb_variable;i++) {
      if (marginal_distribution[i]) {
        observation_distribution[i] = new FrequencyDistribution*[nb_state];
        for (j = 0;j < nb_state;j++) {
          observation_distribution[i][j] = new FrequencyDistribution(marginal_distribution[i]->nb_value);
        }

        observation_frequency_distribution_computation(i , nb_state);
      }

      else {
        observation_distribution[i] = NULL;
      }
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Construction des histogrammes d'observation pour une variable.
 *
 *  arguments : indice de la variable, nombre d'etats, pas de regroupement.
 *
 *--------------------------------------------------------------*/

void MarkovianSequences::build_observation_histogram(int variable , int nb_state , double step)

{
  if ((!observation_histogram[variable]) || (step != observation_histogram[variable][0]->step)) {
    register int i , j;
    int *pstate , *pioutput;
    double imin_value , *proutput;


    // construction de l'histogramme

    if (step == D_DEFAULT) {
      step = marginal_histogram[variable]->step;
    }
    imin_value = floor(min_value[variable] / step) * step;

    if (observation_histogram[variable]) {
      for (i = 0;i < nb_state;i++) {
        observation_histogram[variable][i]->nb_category = (int)floor((max_value[variable] - imin_value) / step) + 1;

        delete [] observation_histogram[variable][i]->frequency;
        observation_histogram[variable][i]->frequency = new int[observation_histogram[variable][i]->nb_category];
      }
    }

    else {
      observation_histogram[variable] = new Histogram*[nb_state];

      for (i = 0;i < nb_state;i++) {
        observation_histogram[variable][i] = new Histogram((int)floor((max_value[variable] - imin_value) / step) + 1 , false);

        observation_histogram[variable][i]->nb_element = marginal_distribution[0]->frequency[i];
        observation_histogram[variable][i]->type = type[variable];
      }

      // calcul des valeurs minimums et maximums par etat

/*      for (i = 0;i < nb_state;i++) {
        observation_histogram[variable][i]->min_value = max_value[variable];
        observation_histogram[variable][i]->max_value = min_value[variable];
      }

      switch (type[variable]) {

      case INT_VALUE : {
        for (i = 0;i < nb_sequence;i++) {
          pstate = int_sequence[i][0];
          pioutput = int_sequence[i][variable];

          for (j = 0;j < length[i];j++) {
            if (*pioutput < observation_histogram[variable][*pstate]->min_value) {
              observation_histogram[variable][*pstate]->min_value = *pioutput;
            }
            if (*pioutput > observation_histogram[variable][*pstate]->max_value) {
              observation_histogram[variable][*pstate]->max_value = *pioutput;
            }
            pstate++;
            pioutput++;
          }
        }
        break;
      }

      case REAL_VALUE : {
        for (i = 0;i < nb_sequence;i++) {
          pstate = int_sequence[i][0];
          proutput = real_sequence[i][variable];

          for (j = 0;j < length[i];j++) {
            if (*proutput < observation_histogram[variable][*pstate]->min_value) {
              observation_histogram[variable][*pstate]->min_value = *proutput;
            }
            if (*proutput > observation_histogram[variable][*pstate]->max_value) {
              observation_histogram[variable][*pstate]->max_value = *proutput;
            }
            pstate++;
            proutput++;
          }
        }
        break;
      }
      } */
    }

    for (i = 0;i < nb_state;i++) {
      observation_histogram[variable][i]->step = step;
      observation_histogram[variable][i]->min_value = imin_value;
      observation_histogram[variable][i]->max_value = ceil(max_value[variable] / step) * step;
    }

    // calcul des frequences

    for (i = 0;i < nb_state;i++) {
      for (j = 0;j < observation_histogram[variable][i]->nb_category;j++) {
        observation_histogram[variable][i]->frequency[j] = 0;
      }
    }

    switch (type[variable]) {

    case INT_VALUE : {
      for (i = 0;i < nb_sequence;i++) {
        pstate = int_sequence[i][0];
        pioutput = int_sequence[i][variable];
        for (j = 0;j < length[i];j++) {
//          (observation_histogram[variable][*pstate++]->frequency[(int)((*pioutput++ - imin_value) / step)])++;
          (observation_histogram[variable][*pstate++]->frequency[(int)floor((*pioutput++ - imin_value) / step)])++;
        }
      }
      break;
    }

    case REAL_VALUE : {
      for (i = 0;i < nb_sequence;i++) {
        pstate = int_sequence[i][0];
        proutput = real_sequence[i][variable];
        for (j = 0;j < length[i];j++) {
//          (observation_histogram[variable][*pstate++]->frequency[(int)((*proutput++ - imin_value) / step)]++;
          (observation_histogram[variable][*pstate++]->frequency[(int)floor((*proutput++ - imin_value) / step)])++;
        }
      }
      break;
    }
    }

    for (i = 0;i < nb_state;i++) {
      observation_histogram[variable][i]->max_computation();
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Construction des histogrammes d'observation.
 *
 *  argument : nombre d'etats.
 *
 *--------------------------------------------------------------*/

void MarkovianSequences::build_observation_histogram(int nb_state)

{
  if ((nb_variable > 1) && (!observation_histogram)) {
    register int i;


    observation_histogram = new Histogram**[nb_variable];
    observation_histogram[0] = NULL;

    for (i = 1;i < nb_variable;i++) {
      observation_histogram[i] = NULL;
      if (marginal_histogram[i]) {
        build_observation_histogram(i , nb_state , marginal_histogram[i]->step);
      }
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Changement du pas de regroupement de l'histogramme marginal et
 *  des histogrammes d'observation pour une variable donnee.
 *
 *  arguments : reference sur un objet StatError, indice de la variable,
 *              pas de regroupement, valeur minimum.
 *
 *--------------------------------------------------------------*/

bool MarkovianSequences::select_step(StatError &error , int variable ,
                                     double step , double imin_value)

{
  bool status = true;


  error.init();

  if ((variable < 1) || (variable > nb_variable)) {
    status = false;
    error.update(STAT_error[STATR_VARIABLE_INDEX]);
  }

  else {
    variable--;

    if (!marginal_histogram[variable]) {
      status = false;
      error.update(STAT_error[STATR_MARGINAL_HISTOGRAM]);
    }
    if ((step <= 0.) || ((type[variable] != REAL_VALUE) &&
         (type[variable] != AUXILIARY) && ((int)step != step))) {
      status = false;
      error.update(STAT_error[STATR_HISTOGRAM_STEP]);
    }
    if ((imin_value != D_INF) && ((imin_value <= min_value[variable] - step) ||
         (imin_value > min_value[variable]) || ((type[variable] != REAL_VALUE) &&
          (type[variable] != AUXILIARY) && ((int)imin_value != imin_value)))) {
      status = false;
      error.update(STAT_error[STATR_HISTOGRAM_MIN_VALUE]);
    }
  }

  if (status) {
    build_marginal_histogram(variable , step , imin_value);

    if ((observation_histogram) && (observation_histogram[variable])) {
      build_observation_histogram(variable , marginal_distribution[0]->nb_value , step);
    }
  }

  return status;
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

    occurrence = new bool*[marginal_distribution[0]->nb_value];
    for (i = 0;i < marginal_distribution[0]->nb_value;i++) {
      occurrence[i] = new bool[marginal_distribution[variable]->nb_value];
      for (j = 0;j < marginal_distribution[variable]->nb_value;j++) {
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

    for (i = 0;i < marginal_distribution[variable]->nb_value;i++) {
      nb_occurrence = 0;
      for (j = 0;j < marginal_distribution[0]->nb_value;j++) {
        if (occurrence[j][i]) {
          nb_occurrence++;
        }
      }

      if (nb_occurrence > 1) {
        hidden = true;
        break;
      }
    }

    for (i = 0;i < marginal_distribution[0]->nb_value;i++) {
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

  characteristics[variable]->index_value = new Curves(marginal_distribution[variable]->nb_value ,
                                                      max_length , true , false , false);
  frequency = new int[marginal_distribution[variable]->nb_value];

  // calcul des probabilites de chaque valeur en fonction de l'index

  for (i = 0;i < max_length;i++) {
    for (j = 0;j < marginal_distribution[variable]->nb_value;j++) {
      frequency[j] = 0;
    }

    for (j = 0;j < nb_sequence;j++) {
      if (i < length[j]) {
        frequency[int_sequence[j][variable][i]]++;
      }
    }

    total = 0;
    for (j = 0;j < marginal_distribution[variable]->nb_value;j++) {
      total += frequency[j];
    }
    characteristics[variable]->index_value->frequency[i] = total;
    for (j = 0;j < marginal_distribution[variable]->nb_value;j++) {
      characteristics[variable]->index_value->point[j][i] = (double)frequency[j] / (double)total;
    }
  }

  delete [] frequency;
}


/*--------------------------------------------------------------*
 *
 *  Extraction des probabilites de chaque valeur entiere en fonction de l'index
 *  (pour une variable donnee).
 *
 *  argument : indice de la variable.
 *
 *--------------------------------------------------------------*/

void MarkovianSequences::build_explicit_index_value(int variable)

{
  register int i , j , k , m;
  int total , *frequency , *index;


  // creation d'un objet Curves

  i = 0;
  for (j = index_parameter_distribution->offset;j < index_parameter_distribution->nb_value;j++) {
    if (index_parameter_distribution->frequency[j] > 0) {
      i++;
    }
  }
  characteristics[variable]->explicit_index_value = new Curves(marginal_distribution[variable]->nb_value ,
                                                               i , true , true , false);

  frequency = new int[marginal_distribution[variable]->nb_value];
  index = new int[nb_sequence];

  // calcul des probabilites de chaque valeur en fonction de l'index explicite

  for (i = 0;i < nb_sequence;i++) {
    index[i] = 0;
  }

  i = 0;
  for (j = index_parameter_distribution->offset;j < index_parameter_distribution->nb_value;j++) {
    if (index_parameter_distribution->frequency[j] > 0) {
      characteristics[variable]->explicit_index_value->index_parameter[i] = j;

      for (k = 0;k < marginal_distribution[variable]->nb_value;k++) {
        frequency[k] = 0;
      }

      for (k = 0;k < nb_sequence;k++) {
        m = 0;
        m = index[k];
        while ((m < length[k] - 1) && (index_parameter[k][m] < j)) {
          m++;
        }
        if (index_parameter[k][m] == j) {
          index[k] = m;
          frequency[int_sequence[k][variable][m]]++;
        }
      }

      total = 0;
      for (k = 0;k < marginal_distribution[variable]->nb_value;k++) {
        total += frequency[k];
      }
      characteristics[variable]->explicit_index_value->frequency[i] = total;
      for (k = 0;k < marginal_distribution[variable]->nb_value;k++) {
        characteristics[variable]->explicit_index_value->point[k][i] = (double)frequency[k] / (double)total;
      }
      i++;
    }
  }

  delete [] frequency;
  delete [] index;
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

  first_occurrence = new FrequencyDistribution*[marginal_distribution[variable]->nb_value];
  for (i = 0;i < marginal_distribution[variable]->nb_value;i++) {
    first_occurrence[i] = new FrequencyDistribution(max_length);
  }

/*  characteristics[variable]->first_occurrence = new FrequencyDistribution*[marginal_distribution[variable]->nb_value];
  for (i = 0;i < marginal_distribution[variable]->nb_value;i++) {
    characteristics[variable]->first_occurrence[i] = new FrequencyDistribution(max_length);
  } */

  // mise a jour des lois empiriques

  occurrence = new bool[marginal_distribution[variable]->nb_value];

  for (i = 0;i < nb_sequence;i++) {
    nb_value = 0;
    for (j = 0;j < marginal_distribution[variable]->nb_value;j++) {
      occurrence[j] = false;
    }

    pisequence = int_sequence[i][variable];
    for (j = 0;j < length[i];j++) {
      if (!occurrence[*pisequence]) {
        occurrence[*pisequence] = true;
        (first_occurrence[*pisequence]->frequency[j])++;
//       (characteristics[variable]->first_occurrence[*pisequence]->frequency[j])++;

        nb_value++;
        if (nb_value == marginal_distribution[variable]->nb_value) {
          break;
        }
      }

      pisequence++;
    }
  }

  delete [] occurrence;

  // extraction des caracteristiques des lois empiriques

  for (i = 0;i < marginal_distribution[variable]->nb_value;i++) {
    first_occurrence[i]->nb_value_computation();
    first_occurrence[i]->offset_computation();
    first_occurrence[i]->nb_element_computation();
    first_occurrence[i]->max_computation();
    first_occurrence[i]->mean_computation();
    first_occurrence[i]->variance_computation();
  }

/*  for (i = 0;i < marginal_distribution[variable]->nb_value;i++) {
    characteristics[variable]->first_occurrence[i]->nb_value_computation();
    characteristics[variable]->first_occurrence[i]->offset_computation();
    characteristics[variable]->first_occurrence[i]->nb_element_computation();
    characteristics[variable]->first_occurrence[i]->max_computation();
    characteristics[variable]->first_occurrence[i]->mean_computation();
    characteristics[variable]->first_occurrence[i]->variance_computation();
  } */

  characteristics[variable]->first_occurrence = new FrequencyDistribution*[marginal_distribution[variable]->nb_value];
  for (i = 0;i < marginal_distribution[variable]->nb_value;i++) {
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

  recurrence_time = new FrequencyDistribution*[marginal_distribution[variable]->nb_value];
  for (i = 0;i < marginal_distribution[variable]->nb_value;i++) {
    recurrence_time[i] = new FrequencyDistribution(max_length);
  }

/*  characteristics[variable]->recurrence_time = new FrequencyDistribution*[marginal_distribution[variable]->nb_value];
  for (i = 0;i < marginal_distribution[variable]->nb_value;i++) {
    characteristics[variable]->recurrence_time[i] = new FrequencyDistribution(max_length);
  } */

  // mise a jour des lois empiriques

  index = new int[marginal_distribution[variable]->nb_value];

  for (i = 0;i < nb_sequence;i++) {
    for (j = 0;j < marginal_distribution[variable]->nb_value;j++) {
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

  for (i = 0;i < marginal_distribution[variable]->nb_value;i++) {
    recurrence_time[i]->nb_value_computation();
    recurrence_time[i]->offset_computation();
    recurrence_time[i]->nb_element_computation();
    recurrence_time[i]->max_computation();
    recurrence_time[i]->mean_computation();
    recurrence_time[i]->variance_computation();
  }

/*  for (i = 0;i < marginal_distribution[variable]->nb_value;i++) {
    characteristics[variable]->recurrence_time[i]->nb_value_computation();
    characteristics[variable]->recurrence_time[i]->offset_computation();
    characteristics[variable]->recurrence_time[i]->nb_element_computation();
    characteristics[variable]->recurrence_time[i]->max_computation();
    characteristics[variable]->recurrence_time[i]->mean_computation();
    characteristics[variable]->recurrence_time[i]->variance_computation();
  } */

  characteristics[variable]->recurrence_time = new FrequencyDistribution*[marginal_distribution[variable]->nb_value];
  for (i = 0;i < marginal_distribution[variable]->nb_value;i++) {
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

  sojourn_time = new FrequencyDistribution*[marginal_distribution[variable]->nb_value];
  for (i = 0;i < marginal_distribution[variable]->nb_value;i++) {
    sojourn_time[i] = new FrequencyDistribution(max_length + 1);
  }

  if (initial_run_flag) {
    initial_run = new FrequencyDistribution*[marginal_distribution[variable]->nb_value];
    for (i = 0;i < marginal_distribution[variable]->nb_value;i++) {
      initial_run[i] = new FrequencyDistribution(max_length + 1);
    }
  }

  final_run = new FrequencyDistribution*[marginal_distribution[variable]->nb_value];
  for (i = 0;i < marginal_distribution[variable]->nb_value;i++) {
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

  for (i = 0;i < marginal_distribution[variable]->nb_value;i++) {
    sojourn_time[i]->nb_value_computation();
    sojourn_time[i]->offset_computation();
    sojourn_time[i]->nb_element_computation();
    sojourn_time[i]->max_computation();
    sojourn_time[i]->mean_computation();
    sojourn_time[i]->variance_computation();
  }

  if (initial_run_flag) {
    for (i = 0;i < marginal_distribution[variable]->nb_value;i++) {
      initial_run[i]->nb_value_computation();
      initial_run[i]->offset_computation();
      initial_run[i]->nb_element_computation();
      initial_run[i]->max_computation();
      initial_run[i]->mean_computation();
      initial_run[i]->variance_computation();
    }
  }

  for (i = 0;i < marginal_distribution[variable]->nb_value;i++) {
    final_run[i]->nb_value_computation();
    final_run[i]->offset_computation();
    final_run[i]->nb_element_computation();
    final_run[i]->max_computation();
    final_run[i]->mean_computation();
    final_run[i]->variance_computation();
  }

  characteristics[variable]->sojourn_time = new FrequencyDistribution*[marginal_distribution[variable]->nb_value];
  for (i = 0;i < marginal_distribution[variable]->nb_value;i++) {
    characteristics[variable]->sojourn_time[i] = new FrequencyDistribution(*(sojourn_time[i]));
    delete sojourn_time[i];
  }
  delete [] sojourn_time;

  if (initial_run_flag) {
    characteristics[variable]->initial_run = new FrequencyDistribution*[marginal_distribution[variable]->nb_value];
    for (i = 0;i < marginal_distribution[variable]->nb_value;i++) {
      characteristics[variable]->initial_run[i] = new FrequencyDistribution(*(initial_run[i]));
      delete initial_run[i];
    }
    delete [] initial_run;
  }

  characteristics[variable]->final_run = new FrequencyDistribution*[marginal_distribution[variable]->nb_value];
  for (i = 0;i < marginal_distribution[variable]->nb_value;i++) {
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
  int run_length , *pisequence;


  // initialisation des lois empiriques

  for (i = 0;i < marginal_distribution[variable]->nb_value;i++) {
    characteristics[variable]->sojourn_time[i]->offset = 1;
    characteristics[variable]->sojourn_time[i]->nb_value = characteristics[variable]->sojourn_time[i]->alloc_nb_value;

    for (j = 0;j < characteristics[variable]->sojourn_time[i]->nb_value;j++) {
      characteristics[variable]->sojourn_time[i]->frequency[j] = 0;
    }

    if (characteristics[variable]->initial_run) {
      characteristics[variable]->initial_run[i]->offset = 1;
      characteristics[variable]->initial_run[i]->nb_value = characteristics[variable]->initial_run[i]->alloc_nb_value;

      for (j = 0;j < characteristics[variable]->initial_run[i]->nb_value;j++) {
        characteristics[variable]->initial_run[i]->frequency[j] = 0;
      }
    }

    characteristics[variable]->final_run[i]->offset = 1;
    characteristics[variable]->final_run[i]->nb_value = characteristics[variable]->final_run[i]->alloc_nb_value;

    for (j = 0;j < characteristics[variable]->final_run[i]->nb_value;j++) {
      characteristics[variable]->final_run[i]->frequency[j] = 0;
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

  for (i = 0;i < marginal_distribution[variable]->nb_value;i++) {
    characteristics[variable]->sojourn_time[i]->nb_value_computation();
    characteristics[variable]->sojourn_time[i]->offset_computation();
    characteristics[variable]->sojourn_time[i]->nb_element_computation();
    characteristics[variable]->sojourn_time[i]->max_computation();
    characteristics[variable]->sojourn_time[i]->mean_computation();
    characteristics[variable]->sojourn_time[i]->variance_computation();
  }

  if (characteristics[variable]->initial_run) {
    for (i = 0;i < marginal_distribution[variable]->nb_value;i++) {
      characteristics[variable]->initial_run[i]->nb_value_computation();
      characteristics[variable]->initial_run[i]->offset_computation();
      characteristics[variable]->initial_run[i]->nb_element_computation();
      characteristics[variable]->initial_run[i]->max_computation();
      characteristics[variable]->initial_run[i]->mean_computation();
      characteristics[variable]->initial_run[i]->variance_computation();
    }
  }

  for (i = 0;i < marginal_distribution[variable]->nb_value;i++) {
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

  for (i = 0;i < marginal_distribution[0]->nb_value;i++) {
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

  nb_run = new FrequencyDistribution*[marginal_distribution[variable]->nb_value];
  for (i = 0;i < marginal_distribution[variable]->nb_value;i++) {
    nb_run[i] = new FrequencyDistribution((max_length % 2 == 0 ?
                               max_length / 2 : max_length / 2 + 1) + 1);
  }

/*  characteristics[variable]->nb_run = new FrequencyDistribution*[marginal_distribution[variable]->nb_value];
  for (i = 0;i < marginal_distribution[variable]->nb_value;i++) {
    characteristics[variable]->nb_run[i] = new FrequencyDistribution((max_length % 2 == 0 ?
                                                                     max_length / 2 : max_length / 2 + 1) + 1);
  } */

  // mise a jour des lois empiriques

  count = new int[marginal_distribution[variable]->nb_value];

  for (i = 0;i < nb_sequence;i++) {
    for (j = 0;j < marginal_distribution[variable]->nb_value;j++) {
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

    for (j = 0;j < marginal_distribution[variable]->nb_value;j++) {
      (nb_run[j]->frequency[count[j]])++;
//      (characteristics[variable]->nb_run[j]->frequency[count[j]])++;
    }
  }

  delete [] count;

  // extraction des caracteristiques des lois empiriques

  for (i = 0;i < marginal_distribution[variable]->nb_value;i++) {
    nb_run[i]->nb_value_computation();
    nb_run[i]->offset_computation();
    nb_run[i]->nb_element = nb_sequence;
    nb_run[i]->max_computation();
    nb_run[i]->mean_computation();
    nb_run[i]->variance_computation();
  }

/*  for (i = 0;i < marginal_distribution[variable]->nb_value;i++) {
    characteristics[variable]->nb_run[i]->nb_value_computation();
    characteristics[variable]->nb_run[i]->offset_computation();
    characteristics[variable]->nb_run[i]->nb_element = nb_sequence;
    characteristics[variable]->nb_run[i]->max_computation();
    characteristics[variable]->nb_run[i]->mean_computation();
    characteristics[variable]->nb_run[i]->variance_computation();
  } */

  characteristics[variable]->nb_run = new FrequencyDistribution*[marginal_distribution[variable]->nb_value];
  for (i = 0;i < marginal_distribution[variable]->nb_value;i++) {
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

  nb_occurrence = new FrequencyDistribution*[marginal_distribution[variable]->nb_value];
  for (i = 0;i < marginal_distribution[variable]->nb_value;i++) {
    nb_occurrence[i] = new FrequencyDistribution(max_length + 1);
  }

/*  characteristics[variable]->nb_occurrence = new FrequencyDistribution*[marginal_distribution[variable]->nb_value];
  for (i = 0;i < marginal_distribution[variable]->nb_value;i++) {
    characteristics[variable]->nb_occurrence[i] = new FrequencyDistribution(max_length + 1);
  } */

  // mise a jour des lois empiriques

  count = new int[marginal_distribution[variable]->nb_value];

  for (i = 0;i < nb_sequence;i++) {
    for (j = 0;j < marginal_distribution[variable]->nb_value;j++) {
      count[j] = 0;
    }

    pisequence = int_sequence[i][variable];
    for (j = 0;j < length[i];j++) {
      count[*pisequence++]++;
    }

    for (j = 0;j < marginal_distribution[variable]->nb_value;j++) {
      (nb_occurrence[j]->frequency[count[j]])++;
//      (characteristics[variable]->nb_occurrence[j]->frequency[count[j]])++;
    }
  }

  delete [] count;

  // extraction des caracteristiques des lois empiriques

  for (i = 0;i < marginal_distribution[variable]->nb_value;i++) {
    nb_occurrence[i]->nb_value_computation();
    nb_occurrence[i]->offset_computation();
    nb_occurrence[i]->nb_element = nb_sequence;
    nb_occurrence[i]->max_computation();
    nb_occurrence[i]->mean_computation();
    nb_occurrence[i]->variance_computation();
  }

/*  for (i = 0;i < marginal_distribution[variable]->nb_value;i++) {
    characteristics[variable]->nb_occurrence[i]->nb_value_computation();
    characteristics[variable]->nb_occurrence[i]->offset_computation();
    characteristics[variable]->nb_occurrence[i]->nb_element = nb_sequence;
    characteristics[variable]->nb_occurrence[i]->max_computation();
    characteristics[variable]->nb_occurrence[i]->mean_computation();
    characteristics[variable]->nb_occurrence[i]->variance_computation();
  } */

  characteristics[variable]->nb_occurrence = new FrequencyDistribution*[marginal_distribution[variable]->nb_value];
  for (i = 0;i < marginal_distribution[variable]->nb_value;i++) {
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
    if (((variable == I_DEFAULT) || (i == variable)) && (marginal_distribution[i])) {
      build = true;

      if (marginal_distribution[i]->nb_value > NB_OUTPUT) {
        build = false;
      }

      else if (type[i] != STATE) {
        for (j = 0;j < marginal_distribution[i]->nb_value;j++) {
          if (marginal_distribution[i]->frequency[j] == 0) {
            build = false;
            break;
          }
        }
      }

      if (build) {
        if (sojourn_time_flag) {
          characteristics[i] = new SequenceCharacteristics(marginal_distribution[i]->nb_value);
        }

        build_index_value(i);
        if (index_parameter) {
          build_explicit_index_value(i);
        }

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
      nb_state = marginal_distribution[variable]->nb_value - marginal_distribution[variable]->offset;

      if ((nb_state < 2) || (nb_state > NB_STATE)) {
        status = false;
        error.update(SEQ_error[SEQR_NB_STATE]);
      }

      else {
        max_nb_word = 0;
        for (i = MAX(length_distribution->offset , word_length);i < length_distribution->nb_value;i++) {
          max_nb_word += length_distribution->frequency[i] * (i - (word_length - 1));
        }
        nb_word_bound = pow((double)nb_state , word_length);
        if (nb_word_bound < max_nb_word) {
          max_nb_word = (int)nb_word_bound;
        }

        if (max_nb_word > MAX_NB_WORD) {
          status = false;
          error.update(SEQ_error[SEQR_MAX_NB_WORD]);
        }

        if ((begin_state != I_DEFAULT) && ((begin_state < marginal_distribution[variable]->offset) ||
             (begin_state >= marginal_distribution[variable]->nb_value) ||
             (marginal_distribution[variable]->frequency[begin_state] == 0))) {
          status = false;
          ostringstream error_message;
          error_message << STAT_label[STATL_VARIABLE] << " " << variable + 1 << ": "
                        << STAT_label[STATL_STATE] << " " << begin_state << " "
                        << STAT_error[STATR_NOT_PRESENT];
          error.update((error_message.str()).c_str());
        }

        if ((end_state != I_DEFAULT) && ((end_state < marginal_distribution[variable]->offset) ||
             (end_state >= marginal_distribution[variable]->nb_value) ||
             (marginal_distribution[variable]->frequency[end_state] == 0))) {
          status = false;
          ostringstream error_message;
          error_message << STAT_label[STATL_VARIABLE] << " " << variable + 1 << ": "
                        << STAT_label[STATL_STATE] << " " << end_state << " "
                        << STAT_error[STATR_NOT_PRESENT];
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
      i *= marginal_distribution[variable]->nb_value;
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
    os << "(" << SEQ_label[SEQL_MIN_INDEX_PARAMETER] << ": " << index_parameter_distribution->offset << ", "
       << SEQ_label[SEQL_MAX_INDEX_PARAMETER] << ": " << index_parameter_distribution->nb_value - 1 << ")" << endl;

    os << "\n";
    if (comment_flag) {
      os << "# ";
    }
    os << SEQ_label[SEQL_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";
    index_parameter_distribution->ascii_characteristic_print(os , false , comment_flag);

    if (exhaustive) {
      os << "\n";
      if (comment_flag) {
        os << "# ";
      }
      os << "   | " << SEQ_label[SEQL_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << endl;
      index_parameter_distribution->ascii_print(os , comment_flag);
    }

    if (index_interval) {
      os << "\n";
      if (comment_flag) {
        os << "# ";
      }
      os << SEQ_label[SEQL_TIME_INTERVAL] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";
      index_interval->ascii_characteristic_print(os , false , comment_flag);

      if (exhaustive) {
        os << "\n";
        if (comment_flag) {
          os << "# ";
        }
        os << "   | " << SEQ_label[SEQL_TIME_INTERVAL] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << endl;
        index_interval->ascii_print(os , comment_flag);
      }
    }

    os << "\n";
  }

  os << nb_variable << " " << STAT_word[nb_variable == 1 ? STATW_VARIABLE : STATW_VARIABLES] << endl;

  if ((self_transition) && (exhaustive)) {
    for (i = 0;i < marginal_distribution[0]->nb_value;i++) {
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
       << STAT_variable_word[type[i]];

    if (type[i] != AUXILIARY) {
      os  << "   ";
      if (comment_flag) {
        os << "# ";
      }

      if (type[i] == STATE) {
        os << "(" << marginal_distribution[i]->nb_value << " "
           << STAT_label[marginal_distribution[i]->nb_value == 1 ? STATL_STATE : STATL_STATES] << ")" << endl;
      }
      else {
        os << "(" << STAT_label[STATL_MIN_VALUE] << ": " << min_value[i] << ", "
           << STAT_label[STATL_MAX_VALUE] << ": " << max_value[i] << ")" << endl;
      }

      if (marginal_distribution[i]) {
        os << "\n";
        if (comment_flag) {
          os << "# ";
        }
        os << STAT_label[type[i] == STATE ? STATL_STATE : STATL_MARGINAL] << " "
           << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";
        marginal_distribution[i]->ascii_characteristic_print(os , false , comment_flag);

        if ((marginal_distribution[i]->nb_value <= ASCII_NB_VALUE) || (exhaustive)) {
          os << "\n";
          if (comment_flag) {
            os << "# ";
          }
          os << "   | " << STAT_label[STATL_FREQUENCY] << endl;
          marginal_distribution[i]->ascii_print(os , comment_flag);
        }

#       ifdef DEBUG
        ContinuousParametric *dist;
;
        marginal_histogram[i] = new Histogram(*marginal_distribution[i]);

        dist = new ContinuousParametric(VON_MISES , 137.5 , 180 * 180 / (50. * 50. * M_PI * M_PI) , DEGREE);
        os << "\n";
        dist->ascii_parameter_print(os);
        dist->ascii_print(os , comment_flag , true , NULL , marginal_distribution[i]);
        os << "\n";
        dist->ascii_print(os , comment_flag , true , marginal_histogram[i]);
        delete dist;

        dist = new ContinuousParametric(GAUSSIAN , 137.5 , 50.);
        os << "\n";
        dist->ascii_parameter_print(os);
        dist->ascii_print(os , comment_flag , true , NULL , marginal_distribution[i]);
//        os << "\n";
//        dist->ascii_print(os , comment_flag , true , marginal_histogram[i]);
        delete dist;

        delete marginal_histogram[i];
        marginal_histogram[i] = NULL;
#       endif

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

        if ((variance > 0.) && (exhaustive)) {
          if (comment_flag) {
            os << "# ";
          }
          os << STAT_label[STATL_SKEWNESS_COEFF] << ": " << skewness_computation(i , mean , variance) << "   "
             << STAT_label[STATL_KURTOSIS_COEFF] << ": " << kurtosis_computation(i , mean , variance) << endl;
        }

        if ((marginal_histogram[i]) && (exhaustive)) {
          os << "\n";
          if (comment_flag) {
            os << "# ";
          }
          os << STAT_label[STATL_MARGINAL] << " " << STAT_label[STATL_HISTOGRAM] << endl;

          os << "\n";
          if (comment_flag) {
            os << "# ";
          }
          os << " " << STAT_label[STATL_VALUE] << "  | " << STAT_label[STATL_FREQUENCY] << endl;
          marginal_histogram[i]->ascii_print(os , comment_flag);
        }
      }

      if (characteristics[i]) {
        characteristics[i]->ascii_print(os , type[i] , *length_distribution , exhaustive , comment_flag);
      }
    }

    else {

#     ifdef MESSAGE
      mean = mean_computation(i);
      variance = variance_computation(i , mean);

      os << "\n";
      if (comment_flag) {
        os << "# ";
      }
      os << STAT_label[STATL_MEAN] << ": " << mean << "   "
         << STAT_label[STATL_VARIANCE] << ": " << variance << "   "
         << STAT_label[STATL_STANDARD_DEVIATION] << ": " << sqrt(variance) << endl;
#     endif

//      os << endl;
    }
  }

  os << "\n";
  if (comment_flag) {
    os << "# ";
  }
  os << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";
  length_distribution->ascii_characteristic_print(os , false , comment_flag);

  if (exhaustive) {
    os << "\n";
    if (comment_flag) {
      os << "# ";
    }
    os << "   | " << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << endl;
    length_distribution->ascii_print(os , comment_flag);
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
               << SEQ_label[SEQL_MIN_INDEX_PARAMETER] << "\t" << index_parameter_distribution->offset << "\t\t"
               << SEQ_label[SEQL_MAX_INDEX_PARAMETER] << "\t" << index_parameter_distribution->nb_value - 1 << endl;

      out_file << "\n" << SEQ_label[SEQL_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\t";
      index_parameter_distribution->spreadsheet_characteristic_print(out_file);

      out_file << "\n\t" << SEQ_label[SEQL_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << endl;
      index_parameter_distribution->spreadsheet_print(out_file);
      out_file << endl;
    }

    out_file << nb_variable << "\t" << STAT_word[nb_variable == 1 ? STATW_VARIABLE : STATW_VARIABLES] << endl;

    if (self_transition) {
      for (i = 0;i < marginal_distribution[0]->nb_value;i++) {
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

      if (type[i] != AUXILIARY) {
        if (type[i] == STATE) {
          out_file << "\t\t" << marginal_distribution[i]->nb_value << "\t"
                   << STAT_label[marginal_distribution[i]->nb_value == 1 ? STATL_STATE : STATL_STATES] << endl;
        }
        else {
          out_file << "\t\t" << STAT_label[STATL_MIN_VALUE] << "\t" << min_value[i]
                   << "\t\t" << STAT_label[STATL_MAX_VALUE] << "\t" << max_value[i] << endl;
        }

        if (marginal_distribution[i]) {
          out_file << "\n" << STAT_label[type[i] == STATE ? STATL_STATE : STATL_MARGINAL] << " "
                   << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\t";
          marginal_distribution[i]->spreadsheet_characteristic_print(out_file);

          out_file << "\n\t" << STAT_label[STATL_FREQUENCY] << endl;
          marginal_distribution[i]->spreadsheet_print(out_file);
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

          if (marginal_histogram[i]) {
            out_file << "\n" << STAT_label[STATL_MARGINAL] << " " << STAT_label[STATL_HISTOGRAM] << endl;
            out_file << "\n" << STAT_label[STATL_VALUE] << "\t" << STAT_label[STATL_FREQUENCY] << endl;
            marginal_histogram[i]->spreadsheet_print(out_file);
          }
        }

        if (characteristics[i]) {
          characteristics[i]->spreadsheet_print(out_file , type[i] , *length_distribution);
        }
      }

      else {
        out_file << endl;
      }
    }

    out_file << "\n" << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\t";
    length_distribution->spreadsheet_characteristic_print(out_file);

    out_file << "\n\t" << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << endl;
    length_distribution->spreadsheet_print(out_file);

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
  const FrequencyDistribution *phisto[1];
  ostringstream data_file_name[2];


  // ecriture du fichier de donnees

  data_file_name[0] << prefix << variable + 1 << 0 << ".dat";

  nb_histo = 0;
  if (index_parameter_distribution) {
    phisto[nb_histo++] = index_parameter_distribution;
  }

  status = length_distribution->plot_print((data_file_name[0].str()).c_str() , nb_histo , phisto);

  if (status) {
    if (marginal_distribution[variable]) {
      data_file_name[1] << prefix << variable + 1 << 1 << ".dat";
      marginal_distribution[variable]->plot_print((data_file_name[1].str()).c_str());
    }
    else if (marginal_histogram[variable]) {
      data_file_name[1] << prefix << variable + 1 << 1 << ".dat";
      marginal_histogram[variable]->plot_print((data_file_name[1].str()).c_str());
    }

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

      if (marginal_distribution[variable]) {
        if (marginal_distribution[variable]->nb_value - 1 < TIC_THRESHOLD) {
          out_file << "set xtics 0,1" << endl;
        }
        if ((int)(marginal_distribution[variable]->max * YSCALE) + 1 < TIC_THRESHOLD) {
          out_file << "set ytics 0,1" << endl;
        }

        out_file << "plot [0:" << MAX(marginal_distribution[variable]->nb_value - 1 , 1) << "] [0:"
                 << (int)(marginal_distribution[variable]->max * YSCALE) + 1 << "] \""
                 << label((data_file_name[1].str()).c_str()) << "\" using 1 title \""
                 << STAT_label[STATL_VARIABLE] << " " << variable + 1 << " - "
                 << STAT_label[type[variable] == STATE ? STATL_STATE : STATL_MARGINAL] << " "
                 << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\" with impulses" << endl;

        if (marginal_distribution[variable]->nb_value - 1 < TIC_THRESHOLD) {
          out_file << "set xtics autofreq" << endl;
        }
        if ((int)(marginal_distribution[variable]->max * YSCALE) + 1 < TIC_THRESHOLD) {
          out_file << "set ytics autofreq" << endl;
        }

        if (i == 0) {
          out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
        }
        out_file << endl;
      }

      else if (marginal_histogram[variable]) {
        if ((int)(marginal_histogram[variable]->max * YSCALE) + 1 < TIC_THRESHOLD) {
          out_file << "set ytics 0,1" << endl;
        }

        out_file << "plot [" << marginal_histogram[variable]->min_value - marginal_histogram[variable]->step << ":"
                 << marginal_histogram[variable]->max_value + marginal_histogram[variable]->step << "] [0:"
                 << (int)(marginal_histogram[variable]->max * YSCALE) + 1 << "] \""
                 << label((data_file_name[1].str()).c_str()) << "\" using 1:2 title \""
                 << STAT_label[STATL_VARIABLE] << " " << variable + 1 << " "
                 << STAT_label[STATL_MARGINAL] << " " << STAT_label[STATL_HISTOGRAM]
                 << "\" with histeps" << endl;

        if ((int)(marginal_histogram[variable]->max * YSCALE) + 1 < TIC_THRESHOLD) {
          out_file << "set ytics autofreq" << endl;
        }

        if (i == 0) {
          out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
        }
        out_file << endl;
      }

      if (length_distribution->nb_value - 1 < TIC_THRESHOLD) {
        out_file << "set xtics 0,1" << endl;
      }
      if ((int)(length_distribution->max * YSCALE) + 1 < TIC_THRESHOLD) {
        out_file << "set ytics 0,1" << endl;
      }

      out_file << "plot [0:" << length_distribution->nb_value - 1 << "] [0:"
               << (int)(length_distribution->max * YSCALE) + 1 << "] \""
               << label((data_file_name[0].str()).c_str()) << "\" using 1 title \""
               << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
               << "\" with impulses" << endl;

      if (length_distribution->nb_value - 1 < TIC_THRESHOLD) {
        out_file << "set xtics autofreq" << endl;
      }
      if ((int)(length_distribution->max * YSCALE) + 1 < TIC_THRESHOLD) {
        out_file << "set ytics autofreq" << endl;
      }

      if (index_parameter_distribution) {
        if (i == 0) {
          out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
        }
        out_file << endl;

        if (index_parameter_distribution->nb_value - 1 < TIC_THRESHOLD) {
          out_file << "set xtics 0,1" << endl;
        }
        if ((int)(index_parameter_distribution->max * YSCALE) + 1 < TIC_THRESHOLD) {
          out_file << "set ytics 0,1" << endl;
        }

        out_file << "plot [" << index_parameter_distribution->offset << ":"
                 << index_parameter_distribution->nb_value - 1 << "] [0:"
                 << (int)(index_parameter_distribution->max * YSCALE) + 1 << "] \""
                 << label((data_file_name[0].str()).c_str()) << "\" using 2 title \""
                 << SEQ_label[SEQL_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
                 << "\" with impulses" << endl;

        if (index_parameter_distribution->nb_value - 1 < TIC_THRESHOLD) {
          out_file << "set xtics autofreq" << endl;
        }
        if ((int)(index_parameter_distribution->max * YSCALE) + 1 < TIC_THRESHOLD) {
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
    status = characteristics[0]->plot_print(prefix , title , 0 , nb_variable , type[0] , *length_distribution);
  }
  else {
    status = plot_print(prefix , title , 0 , nb_variable);
  }

  if (status) {
    for (i = 1;i < nb_variable;i++) {
      if (characteristics[i]) {
        characteristics[i]->plot_print(prefix , title , i , nb_variable , type[i] , *length_distribution);
      }
      else {
        plot_print(prefix , title , i , nb_variable);
      }
    }

    if (self_transition) {
      for (i = 0;i < marginal_distribution[0]->nb_value;i++) {
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
        for (j = 0;j < marginal_distribution[0]->nb_value;j++) {
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
  if ((marginal_distribution[variable]) || (marginal_histogram[variable])) {
    nb_plot_set++;
  }
  if (index_parameter_distribution) {
    nb_plot_set++;
  } */

  plot.variable_nb_viewpoint[variable] = 1;

  if (marginal_distribution[variable]) {

    // vue : loi marginale empirique

    plot.variable[index] = variable;

    plot[index].xrange = Range(0 , MAX(marginal_distribution[variable]->nb_value - 1 , 1));
    plot[index].yrange = Range(0 , ceil(marginal_distribution[variable]->max * YSCALE));

    if (marginal_distribution[variable]->nb_value - 1 < TIC_THRESHOLD) {
      plot[index].xtics = 1;
    }
    if (ceil(marginal_distribution[variable]->max * YSCALE) < TIC_THRESHOLD) {
      plot[index].ytics = 1;
    }

    plot[index].resize(1);

    legend.str("");
    legend << STAT_label[STATL_VARIABLE] << " " << variable + 1 << " - "
           << STAT_label[STATL_MARGINAL] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
    plot[index][0].legend = legend.str();

    plot[index][0].style = "impulses";

    marginal_distribution[variable]->plotable_frequency_write(plot[index][0]);
  }

  else if (marginal_histogram[variable]) {

    // vue : histogramme marginal

    plot[index].xrange = Range(marginal_histogram[variable]->min_value - marginal_histogram[variable]->step ,
                               marginal_histogram[variable]->max_value + marginal_histogram[variable]->step);
    plot[index].yrange = Range(0 , ceil(marginal_histogram[variable]->max * YSCALE));

    if (ceil(marginal_histogram[variable]->max * YSCALE) < TIC_THRESHOLD) {
      plot[index].ytics = 1;
    }

    plot[index].resize(1);

    legend.str("");
    legend << STAT_label[STATL_VARIABLE] << " " << variable + 1 << " "
           << STAT_label[STATL_MARGINAL] << " " << STAT_label[STATL_HISTOGRAM];
    plot[index][0].legend = legend.str();

    plot[index][0].style = "histeps";

    marginal_histogram[variable]->plotable_write(plot[index][0]);
  }

  index++;

  // vue : loi empirique des longueurs des sequences

  plot.variable[index] = variable;

  plot[index].xrange = Range(0 , length_distribution->nb_value - 1);
  plot[index].yrange = Range(0 , ceil(length_distribution->max * YSCALE));

  if (length_distribution->nb_value - 1 < TIC_THRESHOLD) {
    plot[index].xtics = 1;
  }
  if (ceil(length_distribution->max * YSCALE) < TIC_THRESHOLD) {
    plot[index].ytics = 1;
  }

  plot[index].resize(1);

  legend.str("");
  legend << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
  plot[index][0].legend = legend.str();

  plot[index][0].style = "impulses";

  length_distribution->plotable_frequency_write(plot[index][0]);
  index++;

  if (index_parameter_distribution) {

    // vue : loi empirique des parametres d'index

    plot.variable[index] = variable;

    plot[index].xrange = Range(index_parameter_distribution->offset , index_parameter_distribution->nb_value - 1);
    plot[index].yrange = Range(0 , ceil(index_parameter_distribution->max * YSCALE));

    if (index_parameter_distribution->nb_value - 1 < TIC_THRESHOLD) {
      plot[index].xtics = 1;
    }
    if (ceil(index_parameter_distribution->max * YSCALE) < TIC_THRESHOLD) {
      plot[index].ytics = 1;
    }

    plot[index].resize(1);

    legend.str("");
    legend << SEQ_label[SEQL_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
    plot[index][0].legend = legend.str();

    plot[index][0].style = "impulses";

    index_parameter_distribution->plotable_frequency_write(plot[index][0]);
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

    else if (type[i] != AUXILIARY) {
      nb_plot_set++;
      if ((marginal_distribution[i]) || (marginal_histogram[i])) {
        nb_plot_set++;
      }
      if (index_parameter_distribution) {
        nb_plot_set++;
      }
    }
  }

  if (self_transition) {
    for (i = 0;i < marginal_distribution[0]->nb_value;i++) {
      if (self_transition[i]) {
        nb_plot_set += 2;
      }
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

    for (i = 0;i < marginal_distribution[0]->nb_value;i++) {
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
      characteristics[i]->plotable_write(plot , index , i , type[i] , *length_distribution);
    }
    else if (type[i] != AUXILIARY) {
      plotable_write(plot , index , i);
    }
  }

  return plot_set;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un modele lineaire d'observation de type tendance ajuste au format tableur.
 *
 *  arguments : stream, variable, pointeur sur un processus d'observation continu.
 *
 *--------------------------------------------------------------*/

ostream& MarkovianSequences::linear_model_spreadsheet_print(ostream &os , int variable ,
                                                            ContinuousParametricProcess *process) const

{
  bool *used_sequence;
  register int i , j , k , m , n , r;
  int frequency , *index;
  double buff;


  if (type[variable] == INT_VALUE) {
    used_sequence = new bool[nb_sequence];
  }
  if (index_parameter_type == TIME) {
    index = new int[nb_sequence];
  }

  os << "\n" << SEQ_label[SEQL_INDEX] << "\t" << STAT_label[STATL_OBSERVATION];
  for (i = 0;i < process->nb_state;i++) {
    os << "\t" << STAT_label[STATL_STATE] << " " << i << " " << STAT_label[STATL_LINEAR_MODEL];
  }
  if (type[variable] == INT_VALUE) {
    os << "\t" << STAT_label[STATL_FREQUENCY];
  }
  os << endl;

  switch (type[0]) {

  case STATE : {
    for (i = 0;i < process->nb_state;i++) {
      switch (index_parameter_type) {

      case IMPLICIT_TYPE : {
        switch (type[variable]) {

        case INT_VALUE : {
          for (j = 0;j < max_length;j++) {
            for (k = 0;k < nb_sequence;k++) {
              used_sequence[k] = false;
            }
            for (k = 0;k < nb_sequence;k++) {
              if ((j < length[k]) && (int_sequence[k][0][j] == i) && (!used_sequence[k])) {
                used_sequence[k] = true;
                os << j << "\t" << int_sequence[k][variable][j];
                for (m = 0;m <= i;m++) {
                  os << "\t";
                }
                os << process->observation[i]->intercept + process->observation[i]->slope * j;

                frequency = 1;
                for (m = k + 1;m < nb_sequence;m++) {
                  if ((j < length[m]) && (int_sequence[m][0][j] == i) &&
                      (int_sequence[m][variable][j] == int_sequence[k][variable][j])) {
                    used_sequence[m] = true;
                    frequency++;
                  }
                }
                for (m = i;m < process->nb_state;m++) {
                  os << "\t";
                }
                os << frequency << endl;
              }
            }
          }
          break;
        }

        case REAL_VALUE : {
          for (j = 0;j < max_length;j++) {
            for (k = 0;k < nb_sequence;k++) {
              if ((j < length[k]) && (int_sequence[k][0][j] == i)) {
                os << j << "\t" << real_sequence[k][variable][j];
                for (m = 0;m <= i;m++) {
                  os << "\t";
                }
                os << process->observation[i]->intercept + process->observation[i]->slope * j << endl;
              }
            }
          }
          break;
        }
        }
        break;
      }

      case TIME : {
        for (j = 0;j < nb_sequence;j++) {
          index[j] = 0;
        }

        switch (type[variable]) {

        case INT_VALUE : {
          for (j = index_parameter_distribution->offset;j < index_parameter_distribution->nb_value;j++) {
            if (index_parameter_distribution->frequency[j] > 0) {
              for (k = 0;k < nb_sequence;k++) {
                used_sequence[k] = false;
              }
              for (k = 0;k < nb_sequence;k++) {
                m = index[k];
                while ((m < length[k] - 1) && (index_parameter[k][m] < j)) {
                  m++;
                }

                if ((index_parameter[k][m] == j) && (int_sequence[k][0][m] == i) && (!used_sequence[k])) {
                  index[k] = m;
                  used_sequence[k] = true;
                  os << j << "\t" << int_sequence[k][variable][m];
                  for (n = 0;n <= i;n++) {
                    os << "\t";
                  }
                  os << process->observation[i]->intercept + process->observation[i]->slope * j;

                  frequency = 1;
                  for (n = k + 1;n < nb_sequence;n++) {
                    r = index[n];
                    while ((r < length[n] - 1) && (index_parameter[n][r] < j)) {
                      r++;
                    }

                    if ((index_parameter[n][r] == j) && (int_sequence[n][0][r] == i) &&
                        (int_sequence[n][variable][r] == int_sequence[k][variable][m])) {
                      index[n] = r;
                      used_sequence[n] = true;
                      frequency++;
                    }
                  }

                  for (n = i;n < process->nb_state;n++) {
                    os << "\t";
                  }
                  os << frequency << endl;
                }
              }
            }
          }
          break;
        }

        case REAL_VALUE : {
          for (j = index_parameter_distribution->offset;j < index_parameter_distribution->nb_value;j++) {
            if (index_parameter_distribution->frequency[j] > 0) {
              for (k = 0;k < nb_sequence;k++) {
                m = index[k];
                while ((m < length[k] - 1) && (index_parameter[k][m] < j)) {
                  m++;
                }

                if ((index_parameter[k][m] == j) && (int_sequence[k][0][m] == i)) {
                  index[k] = m;
                  os << j << "\t" << real_sequence[k][variable][m];
                  for (n = 0;n <= i;n++) {
                    os << "\t";
                  }
                  os << process->observation[i]->intercept + process->observation[i]->slope * j << endl;
                }
              }
            }
          }
          break;
        }
        }
        break;
      }
      }
    }
    break;
  }

  default : {
    switch (index_parameter_type) {

    case IMPLICIT_TYPE : {
      switch (type[variable]) {

      case INT_VALUE : {
        for (i = 0;i < max_length;i++) {
          for (j = 0;j < nb_sequence;j++) {
            used_sequence[j] = false;
          }
          for (j = 0;j < nb_sequence;j++) {
            if ((i < length[j]) && (!used_sequence[j])) {
              used_sequence[j] = true;
              os << i << "\t" << int_sequence[j][variable][i];
              for (k = 0;k < process->nb_state;k++) {
                os << "\t";
                buff = process->observation[k]->intercept + process->observation[k]->slope * i;
                if ((buff >= min_value[variable]) && (buff <= max_value[variable])) {
                  os << buff;
                }
              }

              frequency = 1;
              for (k = j + 1;k < nb_sequence;k++) {
                if ((i < length[k]) && (int_sequence[k][variable][i] == int_sequence[j][variable][i])) {
                  used_sequence[k] = true;
                  frequency++;
                }
              }
              os << "\t" << frequency << endl;
            }
          }
        }
        break;
      }

      case REAL_VALUE : {
        for (i = 0;i < max_length;i++) {
          for (j = 0;j < nb_sequence;j++) {
            if (i < length[j]) {
              os << i << "\t" << real_sequence[j][variable][i];
              for (k = 0;k < process->nb_state;k++) {
                os << "\t";
                buff = process->observation[k]->intercept + process->observation[k]->slope * i;
                if ((buff >= min_value[variable]) && ( buff <= max_value[variable])) {
                  os << buff;
                }
              }
              os << endl;
            }
          }
        }
        break;
      }
      }
      break;
    }

    case TIME : {
      for (i = 0;i < nb_sequence;i++) {
        index[i] = 0;
      }

      switch (type[variable]) {

      case INT_VALUE : {
        for (i = index_parameter_distribution->offset;i < index_parameter_distribution->nb_value;i++) {
          if (index_parameter_distribution->frequency[i] > 0) {
            for (j = 0;j < nb_sequence;j++) {
              used_sequence[j] = false;
            }
            for (j = 0;j < nb_sequence;j++) {
              k = index[j];
              while ((k < length[j] - 1) && (index_parameter[j][k] < i)) {
                k++;
              }

              if ((index_parameter[j][k] == i) && (!used_sequence[j])) {
                index[j] = k;
                used_sequence[j] = true;
                os << i << "\t" << int_sequence[j][variable][k];
                for (m = 0;m < process->nb_state;m++) {
                  os << "\t";
                  buff = process->observation[m]->intercept + process->observation[m]->slope * i;
                  if ((buff >= min_value[variable]) && (buff <= max_value[variable])) {
                    os << buff;
                  }
                }

                frequency = 1;
                for (m = j + 1;m < nb_sequence;m++) {
                  n = index[m];
                  while ((n < length[m] - 1) && (index_parameter[m][n] < i)) {
                    n++;
                  }

                  if ((index_parameter[m][n] == i) && (int_sequence[m][variable][n] == int_sequence[j][variable][k])) {
                    index[m] = n;
                    used_sequence[m] = true;
                    frequency++;
                  }
                }
                os << "\t" << frequency << endl;
              }
            }
          }
        }
        break;
      }

      case REAL_VALUE : {
        for (i = index_parameter_distribution->offset;i < index_parameter_distribution->nb_value;i++) {
          if (index_parameter_distribution->frequency[i] > 0) {
            for (j = 0;j < nb_sequence;j++) {
              k = index[j];
              while ((k < length[j] - 1) && (index_parameter[j][k] < i)) {
                k++;
              }

              if (index_parameter[j][k] == i) {
                index[j] = k;
                os << i << "\t" << real_sequence[j][variable][k];
                for (m = 0;m < process->nb_state;m++) {
                  os << "\t";
                  buff = process->observation[m]->intercept + process->observation[m]->slope * i;
                  if ((buff >= min_value[variable]) && (buff <= max_value[variable])) {
                    os << buff;
                  }
                }
                os << endl;
              }
            }
          }
        }
        break;
      }
      }
      break;
    }
    }
    break;
  }
  }

  if (type[variable] == INT_VALUE) {
    delete [] used_sequence;
  }
  if (index_parameter_type == TIME) {
    delete [] index;
  }

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Sortie Gnuplot d'un modele lineaire d'observation de type tendance ajuste.
 *
 *  arguments : prefixe des fichiers, titre des figures, variable,
 *              pointeur sur un processus d'observation continu.
 *
 *--------------------------------------------------------------*/

bool MarkovianSequences::linear_model_plot_print(const char *prefix , const char *title , int variable ,
                                                 ContinuousParametricProcess *process) const

{
  bool status = false;
  register int i , j;
  int process_index , *state_min_index_parameter , *state_max_index_parameter , *pstate;
  double buff , *state_min_value , *state_max_value;
  ostringstream data_file_name[NB_STATE * 2 + 1];
  ofstream *out_data_file[NB_STATE + 1];


  // ecriture des fichier de donnees

  state_min_index_parameter = new int[process->nb_state + 1];
  state_max_index_parameter = new int[process->nb_state + 1];

  switch (type[0]) {

  case STATE : {
    process_index = variable;

    switch (index_parameter_type) {

    case IMPLICIT_TYPE : {
      for (i = 0;i < process->nb_state;i++) {
        state_min_index_parameter[i] = max_length - 1;
        state_max_index_parameter[i] = 0;
      }

      for (i = 0;i < nb_sequence;i++) {
        pstate = int_sequence[i][0];

        for (j = 0;j < length[i];j++) {
          if (j < state_min_index_parameter[*pstate]) {
            state_min_index_parameter[*pstate] = j;
          }
          if (j > state_max_index_parameter[*pstate]) {
            state_max_index_parameter[*pstate] = j;
          }
          pstate++;
        }
      }

      state_min_index_parameter[process->nb_state] = 0;
      state_max_index_parameter[process->nb_state] = max_length - 1;
      break;
    }

    case TIME : {
      for (i = 0;i < process->nb_state;i++) {
        state_min_index_parameter[i] = index_parameter_distribution->nb_value - 1;
        state_max_index_parameter[i] = index_parameter_distribution->offset;
      }

      for (i = 0;i < nb_sequence;i++) {
        pstate = int_sequence[i][0];

        for (j = 0;j < length[i];j++) {
          if (index_parameter[i][j] < state_min_index_parameter[*pstate]) {
            state_min_index_parameter[*pstate] = index_parameter[i][j];
          }
          if (index_parameter[i][j] > state_max_index_parameter[*pstate]) {
            state_max_index_parameter[*pstate] = index_parameter[i][j];
          }
          pstate++;
        }
      }

      state_min_index_parameter[process->nb_state] = index_parameter_distribution->offset;
      state_max_index_parameter[process->nb_state] = index_parameter_distribution->nb_value - 1;
      break;
    }
    }

    for (i = 0;i < process->nb_state;i++) {
      if (state_max_index_parameter[i] == state_min_index_parameter[i]) {
        state_max_index_parameter[i]++;
      }

      if (marginal_distribution[0]->frequency[i] == 0) {
        switch (index_parameter_type) {

        case IMPLICIT_TYPE : {
          state_min_index_parameter[i] = 0;
          state_max_index_parameter[i] = max_length - 1;
          break;
        }

        case TIME : {
          state_min_index_parameter[i] = index_parameter_distribution->offset;
          state_max_index_parameter[i] = index_parameter_distribution->nb_value - 1;
          break;
        }
        }

        buff = (min_value[variable] - process->observation[i]->intercept) / process->observation[i]->slope;
        if ((process->observation[i]->slope > 0.) && (buff > state_min_index_parameter[i])) {
          state_min_index_parameter[i] = ceil(buff);
        }
        if ((process->observation[i]->slope < 0.) && (buff < state_max_index_parameter[i])) {
          state_max_index_parameter[i] = floor(buff);
        }

        buff = (max_value[variable] - process->observation[i]->intercept) / process->observation[i]->slope;
        if ((process->observation[i]->slope < 0.) && (buff > state_min_index_parameter[i])) {
          state_min_index_parameter[i] = ceil(buff);
        }
        if ((process->observation[i]->slope > 0.) && (buff < state_max_index_parameter[i])) {
          state_max_index_parameter[i] = floor(buff);
        }
      }
    }

    state_min_value = new double[process->nb_state];
    state_max_value = new double[process->nb_state];
 
    for (i = 0;i < process->nb_state;i++) {
      state_min_value[i] = max_value[variable];
      state_max_value[i] = min_value[variable];
    }

    switch (type[variable]) {

    case INT_VALUE : {
      for (i = 0;i < nb_sequence;i++) {
        pstate = int_sequence[i][0];

        for (j = 0;j < length[i];j++) {
          if (int_sequence[i][variable][j] < state_min_value[*pstate]) {
            state_min_value[*pstate] = int_sequence[i][variable][j];
          }
          if (int_sequence[i][variable][j] > state_max_value[*pstate]) {
            state_max_value[*pstate] = int_sequence[i][variable][j];
          }
          pstate++;
        }
      }
      break;
    }

    case REAL_VALUE : {
      for (i = 0;i < nb_sequence;i++) {
        pstate = int_sequence[i][0];

        for (j = 0;j < length[i];j++) {
          if (real_sequence[i][variable][j] < state_min_value[*pstate]) {
            state_min_value[*pstate] = real_sequence[i][variable][j];
          }
          if (real_sequence[i][variable][j] > state_max_value[*pstate]) {
            state_max_value[*pstate] = real_sequence[i][variable][j];
          }
          pstate++;
        }
      }
      break;
    }
    }
    break;
  }

  default : {
    process_index = variable + 1;

    switch (index_parameter_type) {

    case IMPLICIT_TYPE : {
      for (i = 0;i <= process->nb_state;i++) {
        state_min_index_parameter[i] = 0;
        state_max_index_parameter[i] = max_length - 1;
      }
      break;
    }

    case TIME : {
      for (i = 0;i <= process->nb_state;i++) {
        state_min_index_parameter[i] = index_parameter_distribution->offset;
        state_max_index_parameter[i] = index_parameter_distribution->nb_value - 1;
      }
      break;
    }
    }

    for (i = 0;i < process->nb_state;i++) {
      buff = (min_value[variable] - process->observation[i]->intercept) / process->observation[i]->slope;
      if ((process->observation[i]->slope > 0.) && (buff > state_min_index_parameter[i])) {
        state_min_index_parameter[i] = ceil(buff);
      }
      if ((process->observation[i]->slope < 0.) && (buff < state_max_index_parameter[i])) {
        state_max_index_parameter[i] = floor(buff);
      }

      buff = (max_value[variable] - process->observation[i]->intercept) / process->observation[i]->slope;
      if ((process->observation[i]->slope < 0.) && (buff > state_min_index_parameter[i])) {
        state_min_index_parameter[i] = ceil(buff);
      }
      if ((process->observation[i]->slope > 0.) && (buff < state_max_index_parameter[i])) {
        state_max_index_parameter[i] = floor(buff);
      }

//      cout << "\n" << STAT_label[STATL_STATE] << " " << i << ": " << state_min_index_parameter[i]
//           << ", " << state_max_index_parameter[i] << endl;
    }
    break;
  }
  }

  for (i = 0;i < process->nb_state;i++) {
    data_file_name[i * 2] << prefix << process_index << i * 2 << ".dat";
    out_data_file[0]= new ofstream((data_file_name[i * 2].str()).c_str());

    if (out_data_file[0]) {
      status = true;

      *out_data_file[0] << state_min_index_parameter[i] << " "
                        << process->observation[i]->intercept + process->observation[i]->slope * state_min_index_parameter[i] << endl;
      *out_data_file[0] << state_max_index_parameter[i] << " "
                        << process->observation[i]->intercept + process->observation[i]->slope * state_max_index_parameter[i] << endl;

      out_data_file[0]->close();
      delete out_data_file[0];
    }
  }

  if (type[0] == STATE) {
    for (i = 0;i < process->nb_state;i++) {
      if (marginal_distribution[0]->frequency[i] > 0) {
        data_file_name[i * 2 + 1] << prefix << process_index << i * 2 + 1 << ".dat";
        out_data_file[i] = new ofstream ((data_file_name[i * 2 + 1].str()).c_str());
      }
    }
    switch (index_parameter_type) {

    case IMPLICIT_TYPE : {
      switch (type[variable]) {

      case INT_VALUE : {
        for (i = 0;i < nb_sequence;i++) {
          pstate = int_sequence[i][0];
          for (j = 0;j < length[i];j++) {
            *out_data_file[*pstate++] << j << " " << int_sequence[i][variable][j] << endl;
           }
        }
        break;
      }

      case REAL_VALUE : {
        for (i = 0;i < nb_sequence;i++) {
          pstate = int_sequence[i][0];
          for (j = 0;j < length[i];j++) {
            *out_data_file[*pstate++] << j << " " << real_sequence[i][variable][j] << endl;
          }
        }
        break;
      }
      }
      break;
    }

    case TIME : {
      switch (type[variable]) {

      case INT_VALUE : {
        for (i = 0;i < nb_sequence;i++) {
          pstate = int_sequence[i][0];
          for (j = 0;j < length[i];j++) {
            *out_data_file[*pstate++] << index_parameter[i][j] << " " << int_sequence[i][variable][j] << endl;
          }
        }
        break;
      }

      case REAL_VALUE : {
        for (i = 0;i < nb_sequence;i++) {
          pstate = int_sequence[i][0];
          for (j = 0;j < length[i];j++) {
            *out_data_file[*pstate++] << index_parameter[i][j] << " " << real_sequence[i][variable][j] << endl;
          }
        }
        break;
      }
      }
      break;
    }
    }

    for (i = 0;i < process->nb_state;i++) {
      if (marginal_distribution[0]->frequency[i] > 0) {
        out_data_file[i]->close();
        delete out_data_file[i];
      }
    }
  }

  data_file_name[process->nb_state * 2] << prefix << process_index << process->nb_state * 2 << ".dat";
  out_data_file[process->nb_state] = new ofstream ((data_file_name[process->nb_state * 2].str()).c_str());

  if (out_data_file[process->nb_state]) {
    switch (index_parameter_type) {

    case IMPLICIT_TYPE : {
      switch (type[variable]) {

      case INT_VALUE : {
        for (i = 0;i < nb_sequence;i++) {
          for (j = 0;j < length[i];j++) {
            *out_data_file[process->nb_state] << j << " " << int_sequence[i][variable][j] << endl;
          }
        }
        break;
      }

      case REAL_VALUE : {
        for (i = 0;i < nb_sequence;i++) {
          for (j = 0;j < length[i];j++) {
            *out_data_file[process->nb_state] << j << " " << real_sequence[i][variable][j] << endl;
          }
        }
        break;
      }
      }
      break;
    }

    case TIME : {
      switch (type[variable]) {

      case INT_VALUE : {
        for (i = 0;i < nb_sequence;i++) {
          for (j = 0;j < length[i];j++) {
            *out_data_file[process->nb_state] << index_parameter[i][j] << " " << int_sequence[i][variable][j] << endl;
          }
        }
        break;
      }

      case REAL_VALUE : {
        for (i = 0;i < nb_sequence;i++) {
          for (j = 0;j < length[i];j++) {
            *out_data_file[process->nb_state] << index_parameter[i][j] << " " << real_sequence[i][variable][j] << endl;
          }
        }
        break;
      }
      }
      break;
    }
    }

    out_data_file[process->nb_state]->close();
    delete out_data_file[process->nb_state];
  }

  if (status) {

    // ecriture du fichiers de commandes et du fichiers d'impression

    for (i = 0;i < 2;i++) {
      ostringstream file_name[2];

      switch (i) {
      case 0 :
        file_name[0] << prefix << process_index << ".plot";
        break;
      case 1 :
        file_name[0] << prefix << process_index << ".print";
        break;
      }

      ofstream out_file((file_name[0].str()).c_str());

      if (i == 1) {
        out_file << "set terminal postscript" << endl;
        file_name[1] << label(prefix) << process << ".ps";
        out_file << "set output \"" << file_name[1].str() << "\"\n\n";
      }

      out_file << "set border 15 lw 0\n" << "set tics out\n" << "set xtics nomirror\n"
               << "set title \"";
      if (title) {
        out_file << title << " - ";
      }
      out_file << STAT_label[STATL_OUTPUT_PROCESS] << " " << process_index << "\"\n\n";

      out_file << "set xlabel \"" << SEQ_label[SEQL_INDEX] << "\"" << endl;
      out_file << "set ylabel \"" << STAT_label[STATL_OBSERVATION] << "\"" << endl;

      if (type[0] == STATE) {
        for (j = 0;j < process->nb_state;j++) {
          if (marginal_distribution[0]->frequency[j] > 0) {
            if (state_max_index_parameter[j] - state_min_index_parameter[j] < TIC_THRESHOLD) {
              out_file << "set xtics " << state_min_index_parameter[j] << ",1" << endl;
            }
            if (state_max_value[j] - state_min_value[j] < TIC_THRESHOLD) {
              out_file << "set ytics " << MIN(state_min_value[j] , 0) << ",1" << endl;
            }

            out_file << "plot [" << state_min_index_parameter[j] << ":" << state_max_index_parameter[j] << "] [";
/*            if ((state_min_value[j] >= 0.) && (state_max_value[j] - state_min_value[j] > state_min_value[j] * PLOT_RANGE_RATIO)) {
              out_file << 0;
            }
            else { */
              out_file << state_min_value[j];
//            }
            out_file << ":" << MAX(state_max_value[j] , state_min_value[j] + 1) << "] \""
                     << label((data_file_name[2 * j + 1].str()).c_str()) << "\" using 1:2 notitle with points,\\" << endl;
            out_file << "\"" << label((data_file_name[2 * j].str()).c_str()) << "\" using 1:2 title \""
                     << STAT_label[STATL_STATE] << " " << j << " " << STAT_label[STATL_OBSERVATION] << " "
                     << STAT_label[STATL_MODEL] << "\" with lines" << endl;

            if (state_max_index_parameter[j] - state_min_index_parameter[j] < TIC_THRESHOLD) {
              out_file << "set xtics autofreq" << endl;
            }
            if (state_max_value[j] - state_min_value[j] < TIC_THRESHOLD) {
              out_file << "set ytics autofreq" << endl;
            }

            if (i == 0) {
              out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
            }
            out_file << endl;
          }
        }
      }

      if (state_max_index_parameter[process->nb_state] - state_min_index_parameter[process->nb_state] < TIC_THRESHOLD) {
        out_file << "set xtics " << state_min_index_parameter[process->nb_state] << ",1" << endl;
      }
      if (max_value[variable] - min_value[variable] < TIC_THRESHOLD) {
        out_file << "set ytics " << MIN(min_value[variable] , 0) << ",1" << endl;
      }

      out_file << "plot [" << state_min_index_parameter[process->nb_state] << ":" << state_max_index_parameter[process->nb_state] << "] [";
/*      if ((min_value[variable] >= 0.) && (max_value[variable] - min_value[variable] > min_value[variable] * PLOT_RANGE_RATIO)) {
        out_file << 0;
      }
      else { */
        out_file << min_value[variable];
//      }
      out_file << ":" << MAX(max_value[variable] , min_value[variable] + 1) << "] \""
               << label((data_file_name[2 * process->nb_state].str()).c_str()) << "\" using 1:2 notitle with points,\\" << endl;
      for (j = 0;j < process->nb_state;j++) {
        out_file << "\"" << label((data_file_name[2 * j].str()).c_str()) << "\" using 1:2 title \""
                 << STAT_label[STATL_STATE] << " " << j << " " << STAT_label[STATL_OBSERVATION] << " "
                 << STAT_label[STATL_MODEL] << "\" with lines";
        if (j < process->nb_state - 1) {
          out_file << ",\\";
        }
        out_file << endl;
      }

      if (state_max_index_parameter[process->nb_state] - state_min_index_parameter[process->nb_state] < TIC_THRESHOLD) {
        out_file << "set xtics autofreq" << endl;
      }
      if (max_value[variable] - min_value[variable] < TIC_THRESHOLD) {
        out_file << "set ytics autofreq" << endl;
      }

      out_file << "set xlabel" << endl;
      out_file << "set ylabel" << endl;

      if (i == 1) {
        out_file << "\nset terminal x11" << endl;
      }

      out_file << "\npause 0 \"" << STAT_label[STATL_END] << "\"" << endl;
    }
  }

  delete [] state_min_index_parameter;
  delete [] state_max_index_parameter;
  if (type[0] == STATE) {
    delete [] state_min_value;
    delete [] state_max_value;
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Sortie graphique d'un modele lineaire d'observation de type tendance ajuste.
 *
 *  arguments : prefixe des fichiers, indice du MultiPlot, variable,
 *              pointeur sur un processus d'observation continu.
 *
 *--------------------------------------------------------------*/

void MarkovianSequences::linear_model_plotable_write(MultiPlotSet &plot , int &index , int variable ,
                                                     ContinuousParametricProcess *process) const

{
  bool status = false;
  register int i , j;
  int process_index , plot_offset , *state_min_index_parameter , *state_max_index_parameter , *pstate;
  double buff , *state_min_value , *state_max_value;
  ostringstream title , legend;


  // calcul des bornes

  state_min_index_parameter = new int[process->nb_state + 1];
  state_max_index_parameter = new int[process->nb_state + 1];

  switch (type[0]) {

  case STATE : {
    process_index = variable;
    plot_offset = process->nb_state;

    switch (index_parameter_type) {

    case IMPLICIT_TYPE : {
      for (i = 0;i < process->nb_state;i++) {
        state_min_index_parameter[i] = max_length - 1;
        state_max_index_parameter[i] = 0;
      }

      for (i = 0;i < nb_sequence;i++) {
        pstate = int_sequence[i][0];

        for (j = 0;j < length[i];j++) {
          if (j < state_min_index_parameter[*pstate]) {
            state_min_index_parameter[*pstate] = j;
          }
          if (j > state_max_index_parameter[*pstate]) {
            state_max_index_parameter[*pstate] = j;
          }
          pstate++;
        }
      }

      state_min_index_parameter[process->nb_state] = 0;
      state_max_index_parameter[process->nb_state] = max_length - 1;
      break;
    }

    case TIME : {
      for (i = 0;i < process->nb_state;i++) {
        state_min_index_parameter[i] = index_parameter_distribution->nb_value - 1;
        state_max_index_parameter[i] = index_parameter_distribution->offset;
      }

      for (i = 0;i < nb_sequence;i++) {
        pstate = int_sequence[i][0];

        for (j = 0;j < length[i];j++) {
          if (index_parameter[i][j] < state_min_index_parameter[*pstate]) {
            state_min_index_parameter[*pstate] = index_parameter[i][j];
          }
          if (index_parameter[i][j] > state_max_index_parameter[*pstate]) {
            state_max_index_parameter[*pstate] = index_parameter[i][j];
          }
          pstate++;
        }
      }

      state_min_index_parameter[process->nb_state] = index_parameter_distribution->offset;
      state_max_index_parameter[process->nb_state] = index_parameter_distribution->nb_value - 1;
      break;
    }
    }

    for (i = 0;i < process->nb_state;i++) {
      if (state_max_index_parameter[i] == state_min_index_parameter[i]) {
        state_max_index_parameter[i]++;
      }

      if (marginal_distribution[0]->frequency[i] == 0) {
        switch (index_parameter_type) {

        case IMPLICIT_TYPE : {
          state_min_index_parameter[i] = 0;
          state_max_index_parameter[i] = max_length - 1;
          break;
        }

        case TIME : {
          state_min_index_parameter[i] = index_parameter_distribution->offset;
          state_max_index_parameter[i] = index_parameter_distribution->nb_value - 1;
          break;
        }
        }

        buff = (min_value[variable] - process->observation[i]->intercept) / process->observation[i]->slope;
        if ((process->observation[i]->slope > 0.) && (buff > state_min_index_parameter[i])) {
          state_min_index_parameter[i] = ceil(buff);
        }
        if ((process->observation[i]->slope < 0.) && (buff < state_max_index_parameter[i])) {
          state_max_index_parameter[i] = floor(buff);
        }

        buff = (max_value[variable] - process->observation[i]->intercept) / process->observation[i]->slope;
        if ((process->observation[i]->slope < 0.) && (buff > state_min_index_parameter[i])) {
          state_min_index_parameter[i] = ceil(buff);
        }
        if ((process->observation[i]->slope > 0.) && (buff < state_max_index_parameter[i])) {
          state_max_index_parameter[i] = floor(buff);
        }
      }
    }

    state_min_value = new double[process->nb_state];
    state_max_value = new double[process->nb_state];

    for (i = 0;i < process->nb_state;i++) {
      state_min_value[i] = max_value[variable];
      state_max_value[i] = min_value[variable];
    }

    switch (type[variable]) {

    case INT_VALUE : {
      for (i = 0;i < nb_sequence;i++) {
        pstate = int_sequence[i][0];

        for (j = 0;j < length[i];j++) {
          if (int_sequence[i][variable][j] < state_min_value[*pstate]) {
            state_min_value[*pstate] = int_sequence[i][variable][j];
          }
          if (int_sequence[i][variable][j] > state_max_value[*pstate]) {
            state_max_value[*pstate] = int_sequence[i][variable][j];
          }
          pstate++;
        }
      }
      break;
    }

    case REAL_VALUE : {
      for (i = 0;i < nb_sequence;i++) {
        pstate = int_sequence[i][0];

        for (j = 0;j < length[i];j++) {
          if (real_sequence[i][variable][j] < state_min_value[*pstate]) {
            state_min_value[*pstate] = real_sequence[i][variable][j];
          }
          if (real_sequence[i][variable][j] > state_max_value[*pstate]) {
            state_max_value[*pstate] = real_sequence[i][variable][j];
          }
          pstate++;
        }
      }
      break;
    }
    }
    break;
  }

  default : {
    process_index = variable + 1;
    plot_offset = 0;

    switch (index_parameter_type) {

    case IMPLICIT_TYPE : {
      for (i = 0;i <= process->nb_state;i++) {
        state_min_index_parameter[i] = 0;
        state_max_index_parameter[i] = max_length - 1;
      }
      break;
    }

    case TIME : {
      for (i = 0;i <= process->nb_state;i++) {
        state_min_index_parameter[i] = index_parameter_distribution->offset;
        state_max_index_parameter[i] = index_parameter_distribution->nb_value - 1;
      }
      break;
    }
    }

    for (i = 0;i < process->nb_state;i++) {
      buff = (min_value[variable] - process->observation[i]->intercept) / process->observation[i]->slope;
      if ((process->observation[i]->slope > 0.) && (buff > state_min_index_parameter[i])) {
        state_min_index_parameter[i] = ceil(buff);
      }
      if ((process->observation[i]->slope < 0.) && (buff < state_max_index_parameter[i])) {
        state_max_index_parameter[i] = floor(buff);
      }

      buff = (max_value[variable] - process->observation[i]->intercept) / process->observation[i]->slope;
      if ((process->observation[i]->slope < 0.) && (buff > state_min_index_parameter[i])) {
        state_min_index_parameter[i] = ceil(buff);
      }
      if ((process->observation[i]->slope > 0.) && (buff < state_max_index_parameter[i])) {
        state_max_index_parameter[i] = floor(buff);
      }
    }
    break;
  }
  }

  plot.variable_nb_viewpoint[variable] = 1;

  title.str("");
  title << STAT_label[STATL_OUTPUT_PROCESS] << " " << process_index;

  // vues observations et fonction lineaire par etat

  if (type[0] == STATE) {
    for (i = 0;i < process->nb_state;i++) {
      plot[index + i].title = title.str();

      plot[index + i].xrange = Range(state_min_index_parameter[i] , state_max_index_parameter[i]);
      plot[index + i].yrange = Range(state_min_value[i] , MAX(state_max_value[i] , state_min_value[i] + 1));

      if (state_max_index_parameter[i] - state_min_index_parameter[i] < TIC_THRESHOLD) {
        plot[index + i].xtics = 1;
      }
      if (state_max_value[i] - state_min_value[i] < TIC_THRESHOLD) {
        plot[index + i].ytics = 1;
      }

      plot[index + i].xlabel = SEQ_label[SEQL_INDEX];
      plot[index + i].ylabel = STAT_label[STATL_OBSERVATION];

      plot[index + i].resize(2);

/*      legend.str("");
      legend << STAT_label[STATL_STATE] << " " << i << " " << STAT_label[STATL_OBSERVATION];
      plot[index + i][0].legend = legend.str(); */

      legend.str("");
      legend << STAT_label[STATL_STATE] << " " << i << " " << STAT_label[STATL_OBSERVATION] << " "
             << STAT_label[STATL_MODEL];
      plot[index + i][1].legend = legend.str();

      plot[index + i][0].style = "points";
      plot[index + i][1].style = "lines";
    }

    switch (index_parameter_type) {

    case IMPLICIT_TYPE : {
      switch (type[variable]) {

      case INT_VALUE : {
        for (i = 0;i < nb_sequence;i++) {
          pstate = int_sequence[i][0];
          for (j = 0;j < length[i];j++) {
            plot[index + *pstate++][0].add_point(j , int_sequence[i][variable][j]);
          }
        }
        break;
      }

      case REAL_VALUE : {
        for (i = 0;i < nb_sequence;i++) {
          pstate = int_sequence[i][0];
          for (j = 0;j < length[i];j++) {
            plot[index + *pstate++][0].add_point(j , real_sequence[i][variable][j]);
          }
        }
        break;
      }
      }
      break;
    }

    case TIME : {
      switch (type[variable]) {

      case INT_VALUE : {
        for (i = 0;i < nb_sequence;i++) {
          pstate = int_sequence[i][0];
          for (j = 0;j < length[i];j++) {
            plot[index + *pstate++][0].add_point(index_parameter[i][j] , int_sequence[i][variable][j]);
          }
        }
        break;
      }

      case REAL_VALUE : {
        for (i = 0;i < nb_sequence;i++) {
          pstate = int_sequence[i][0];
          for (j = 0;j < length[i];j++) {
            plot[index + *pstate++][0].add_point(index_parameter[i][j] , real_sequence[i][variable][j]);
          }
        }
        break;
      }
      }
      break;
    }
    }

    for (i = 0;i < process->nb_state;i++) {
      plot[index + i][1].add_point(state_min_index_parameter[i] , process->observation[i]->intercept +
                                   process->observation[i]->slope * state_min_index_parameter[i]);
      plot[index + i][1].add_point(state_max_index_parameter[i] , process->observation[i]->intercept +
                                   process->observation[i]->slope * state_max_index_parameter[i]);
    }
  }

  // vue observations et fonctions lineaires

  plot[index + plot_offset].title = title.str();

  plot[index + plot_offset].xrange = Range(state_min_index_parameter[process->nb_state] ,
                                           state_max_index_parameter[process->nb_state]);
  plot[index + plot_offset].yrange = Range(min_value[variable] , MAX(max_value[variable] , min_value[variable]));

  if (state_max_index_parameter[process->nb_state] - state_min_index_parameter[process->nb_state] < TIC_THRESHOLD) {
    plot[index + plot_offset].xtics = 1;
  }
  if (max_value[variable] - min_value[variable] < TIC_THRESHOLD) {
    plot[index + plot_offset].ytics = 1;
  }

  plot[index + plot_offset].xlabel = SEQ_label[SEQL_INDEX];
  plot[index + plot_offset].ylabel = STAT_label[STATL_OBSERVATION];

  plot[index + plot_offset].resize(process->nb_state + 1);

/*  legend.str("");
  legend << STAT_label[STATL_OBSERVATION];
  plot[index + plot_offset][0].legend = legend.str(); */

  plot[index + plot_offset][0].style = "points";

  for (i = 0;i < process->nb_state;i++) {
    legend.str("");
    legend << STAT_label[STATL_STATE] << " " << i << " " << STAT_label[STATL_OBSERVATION] << " "
           << STAT_label[STATL_MODEL];
    plot[index + plot_offset][i + 1].legend = legend.str();

    plot[index + plot_offset][i + 1].style = "lines";
  }

  switch (index_parameter_type) {

  case IMPLICIT_TYPE : {
    switch (type[variable]) {

    case INT_VALUE : {
      for (i = 0;i < nb_sequence;i++) {
        for (j = 0;j < length[i];j++) {
          plot[index + plot_offset][0].add_point(j , int_sequence[i][variable][j]);
        }
      }
      break;
    }

    case REAL_VALUE : {
      for (i = 0;i < nb_sequence;i++) {
        for (j = 0;j < length[i];j++) {
          plot[index + plot_offset][0].add_point(j , real_sequence[i][variable][j]);
        }
      }
      break;
    }
    }
    break;
  }

  case TIME : {
    switch (type[variable]) {

    case INT_VALUE : {
      for (i = 0;i < nb_sequence;i++) {
        for (j = 0;j < length[i];j++) {
          plot[index + plot_offset][0].add_point(index_parameter[i][j] , int_sequence[i][variable][j]);
        }
      }
      break;
    }

    case REAL_VALUE : {
      for (i = 0;i < nb_sequence;i++) {
        for (j = 0;j < length[i];j++) {
          plot[index + plot_offset][0].add_point(index_parameter[i][j] , real_sequence[i][variable][j]);
        }
      }
      break;
    }
    }
    break;
  }
  }

  for (i = 0;i < process->nb_state;i++) {
    plot[index + plot_offset][i + 1].add_point(state_min_index_parameter[i] , process->observation[i]->intercept +
                                               process->observation[i]->slope * state_min_index_parameter[i]);
    plot[index + plot_offset][i + 1].add_point(state_max_index_parameter[i] , process->observation[i]->intercept +
                                               process->observation[i]->slope * state_max_index_parameter[i]);
  }

  switch (type[0]) {
  case STATE :
    index += process->nb_state + 1;
    break;
  default :
    index++;
    break;
  }

  delete [] state_min_index_parameter;
  delete [] state_max_index_parameter;
  if (type[0] == STATE) {
    delete [] state_min_value;
    delete [] state_max_value;
  }
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
        for (j = 1;j < marginal_distribution[i]->nb_value;j++) {
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
        for (j = 1;j < marginal_distribution[i]->nb_value;j++) {
          out_file << (char)('F' + j);
          if (j < marginal_distribution[i]->nb_value - 1) {
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


};  // namespace sequence_analysis
