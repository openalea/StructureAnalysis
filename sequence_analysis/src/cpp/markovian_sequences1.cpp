/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       V-Plants: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2017 CIRAD/INRA/Inria Virtual Plants
 *
 *       File author(s): Yann Guedon (yann.guedon@cirad.fr)
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



// #include <string>
#include <vector>
#include <sstream>
#include <iomanip>

#include <boost/math/distributions/normal.hpp>

#include "stat_tool/stat_label.h"

// #include "stat_tool/quantile_computation.hpp"

#include "sequences.h"
#include "sequence_label.h"

using namespace std;
using namespace stat_tool;


namespace sequence_analysis {



/*--------------------------------------------------------------*/
/**
 *  \brief Default constructor of the MarkovianSequences class.
 */
/*--------------------------------------------------------------*/

MarkovianSequences::MarkovianSequences()

{
  min_interval = NULL;

  self_transition = NULL;
  observation_distribution = NULL;
  observation_histogram = NULL;
  characteristics = NULL;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Default initialization of the data members of the MarkovianSequences class.
 */
/*--------------------------------------------------------------*/

void MarkovianSequences::init()

{
  int i;


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


/*--------------------------------------------------------------*/
/**
 *  \brief Construction of a MarkovianSequences object from a Sequences object.
 *
 *  \param[in] seq reference on a Sequences object.
 */
/*--------------------------------------------------------------*/

MarkovianSequences::MarkovianSequences(const Sequences &seq)
:Sequences(seq)

{
  int i;


  init();

  min_interval = new double[nb_variable];
  for (i = 0;i < nb_variable;i++) {
    min_interval_computation(i);
  }

  build_characteristic();
}


/*--------------------------------------------------------------*/
/**
 *  \brief Construction of a MarkovianSequences object adding auxiliary variables.
 *
 *  \param[in] seq       reference on a MarkovianSequences object,
 *  \param[in] auxiliary flags on the addition of auxiliary variables.
 */
/*--------------------------------------------------------------*/

MarkovianSequences::MarkovianSequences(const MarkovianSequences &seq , bool *auxiliary)
:Sequences(seq , auxiliary)

{
  int i , j;


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


/*--------------------------------------------------------------*/
/**
 *  \brief Copy of a MarkovianSequences object.
 *
 *  \param[in] seq   reference on a MarkovianSequences object,
 *  \param[in] param addition/removing of the initial run length frequency distributions.
 */
/*--------------------------------------------------------------*/

void MarkovianSequences::copy(const MarkovianSequences &seq , initial_run param)

{
  bool initial_run_flag;
  int i , j;


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


/*--------------------------------------------------------------*/
/**
 *  \brief Copy of a MarkovianSequences object reversing the direction of sequences.
 *
 *  \param[in] seq reference on a MarkovianSequences object.
 */
/*--------------------------------------------------------------*/

void MarkovianSequences::reverse(const MarkovianSequences &seq)

{
  int i , j;


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
      characteristics[i] = new SequenceCharacteristics(*(seq.characteristics[i]) , REVERSE);

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


/*--------------------------------------------------------------*/
/**
 *  \brief Copy of a MarkovianSequences object adding a state variable.
 *
 *  \param[in] seq   reference on a MarkovianSequences object,
 *  \param[in] param addition/removing of the initial run length frequency distributions.
 */
/*--------------------------------------------------------------*/

void MarkovianSequences::add_state_variable(const MarkovianSequences &seq , initial_run param)

{
  bool initial_run_flag;
  int i;


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


/*--------------------------------------------------------------*/
/**
 *  \brief Constructor by copy of the MarkovianSequences class.
 *
 *  \param[in] seq       reference on a MarkovianSequences object,
 *  \param[in] transform type of transform,
 *  \param[in] param     addition/removing of the initial run length frequency distributions.
 */
/*--------------------------------------------------------------*/

MarkovianSequences::MarkovianSequences(const MarkovianSequences &seq , sequence_transformation transform ,
                                       initial_run param)

{
  switch (transform) {
  case SEQUENCE_COPY :
    Sequences::copy(seq);
    copy(seq , param);
    break;
  case REVERSE :
    Sequences::reverse(seq);
    reverse(seq);
    break;
  case ADD_STATE_VARIABLE :
    Sequences::add_state_variable(seq);
    add_state_variable(seq , param);
    break;
  case EXPLICIT_INDEX_PARAMETER :
    Sequences::explicit_index_parameter(seq);
    copy(seq);
    break;
  case REMOVE_INDEX_PARAMETER :
    Sequences::remove_index_parameter(seq);
    copy(seq);
    break;
  default :
    Sequences::copy(seq);
    copy(seq);
    break;
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Destructor of data members of the MarkovianSequences class.
 */
/*--------------------------------------------------------------*/

void MarkovianSequences::remove()

{
  int i , j;


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


/*--------------------------------------------------------------*/
/**
 *  \brief Destructor of the MarkovianSequences class.
 */
/*--------------------------------------------------------------*/

MarkovianSequences::~MarkovianSequences()

{
  remove();
}


/*--------------------------------------------------------------*/
/**
 *  \brief Assignement operator of the MarkovianSequences class.
 *
 *  \param[in] seq reference on a MarkovianSequences object.
 *
 *  \return        MarkovianSequences object.
 */
/*--------------------------------------------------------------*/

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


/*--------------------------------------------------------------*/
/**
 *  \brief Initialization of a state variable (1st variable).
 *
 *  \param[in] itype 1st variable type (STATE/INT_VALUE/REAL_VALUE).
 */
/*--------------------------------------------------------------*/

void MarkovianSequences::state_variable_init(variable_nature itype)

{
  int i , j;


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


/*--------------------------------------------------------------*/
/**
 *  \brief Extraction of a frequency distribution.
 *
 *  \param[in] error    reference on a StatError object,
 *  \param[in] type     frequency distribution type,
 *  \param[in] variable variable index,
 *  \param[in] value    value.
 *
 *  \return             DiscreteDistributionData object.
 */
/*--------------------------------------------------------------*/

DiscreteDistributionData* MarkovianSequences::extract(StatError &error , process_distribution type ,
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


/*--------------------------------------------------------------*/
/**
 *  \brief Merging of MarkovianSequences objects.
 *
 *  \param[in] error     reference on a StatError object,
 *  \param[in] nb_sample number of MarkovianSequences objects,
 *  \param[in] iseq      pointer on the MarkovianSequences objects.
 *
 *  \return              MarkovianSequences object.
 */
/*--------------------------------------------------------------*/

MarkovianSequences* MarkovianSequences::merge(StatError &error , int nb_sample ,
                                              const MarkovianSequences **iseq) const

{
  bool status = true;
  int i , j , k , m , n , p , q;
  int inb_sequence , cumul_nb_sequence , nb_histo , *ilength , *iidentifier ,
      **ivertex_identifier;
  const FrequencyDistribution **phisto;
  MarkovianSequences *seq;
  const MarkovianSequences **pseq;


  seq = NULL;
  error.init();

  for (i = 0;i < nb_sample;i++) {
    if (iseq[i]->index_param_type != index_param_type) {
      status = false;
      ostringstream error_message;
      error_message << STAT_label[STATL_SAMPLE] << " " << i + 2 << ": "
                    << SEQ_error[SEQR_INDEX_PARAMETER_TYPE];

      if (index_param_type == IMPLICIT_TYPE) {
        error.update((error_message.str()).c_str());
      }
      else {
        error.correction_update((error_message.str()).c_str() , SEQ_index_parameter_word[index_param_type]);
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

    // computation of the number of sequences

    inb_sequence = 0;
    for (i = 0;i < nb_sample;i++) {
      inb_sequence += pseq[i]->nb_sequence;
    }

    // comparison of sequence identifiers

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

    // copy of sequence lengths

    ilength = new int[inb_sequence];

    i = 0;
    for (j = 0;j < nb_sample;j++) {
      for (k = 0;k < pseq[j]->nb_sequence;k++) {
        ilength[i++] = pseq[j]->length[k];
      }
    }

    // comparison of vertex identifiers

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
                                 index_param_type , nb_variable , type);
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

    if (index_param_type == TIME) {
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

    // copy of values

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


/*--------------------------------------------------------------*/
/**
 *  \brief Merging of MarkovianSequences objects.
 *
 *  \param[in] error     reference on a StatError object,
 *  \param[in] nb_sample number of MarkovianSequences objects,
 *  \param[in] iseq      pointer on the MarkovianSequences objects.
 *
 *  \return              MarkovianSequences object.
 */
/*--------------------------------------------------------------*/

MarkovianSequences* MarkovianSequences::merge(StatError &error , int nb_sample ,
                                              const vector<MarkovianSequences> iseq) const

{
  int i;
  MarkovianSequences *seq;
  const MarkovianSequences **pseq;


  pseq = new const MarkovianSequences*[nb_sample];
  for (i = 0;i < nb_sample;i++) {
    pseq[i] = new MarkovianSequences(iseq[i]);
  }

  seq = merge(error , nb_sample , pseq);

  for (i = 0;i < nb_sample;i++) {
    delete pseq[i];
  }
  delete [] pseq;

  return seq;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Clustering of values of a variable.
 *
 *  \param[in] error    reference on a StatError object,
 *  \param[in] variable variable index,
 *  \param[in] step     clustering step,
 *  \param[in] mode     mode (FLOOR/ROUND/CEIL).
 *
 *  \return             MarkovianSequences object.
 */
/*--------------------------------------------------------------*/

MarkovianSequences* MarkovianSequences::cluster(StatError &error , int variable ,
                                                int step , rounding mode) const

{
  bool status = true;
  int i;
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


/*--------------------------------------------------------------*/
/**
 *  \brief Transcoding of categories of an integer-valued variable.
 *
 *  \param[in] error        reference on a StatError object,
 *  \param[in] ivariable    variable index,
 *  \param[in] category     transcoding table,
 *  \param[in] add_variable flag for adding a variable.
 *
 *  \return                 MarkovianSequences object.
 */
/*--------------------------------------------------------------*/

MarkovianSequences* MarkovianSequences::transcode(StatError &error , int ivariable ,
                                                  int *category , bool add_variable) const

{
  bool status = true , *presence;
  int i;
  int variable , offset , min_category , max_category;
  variable_nature *itype;
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
      min_category = marginal_distribution[ivariable]->nb_value;
      max_category = 0;

      for (i = 0;i < marginal_distribution[ivariable]->nb_value - marginal_distribution[ivariable]->offset;i++) {
        if ((category[i] < 0) || (category[i] >= (add_variable ? marginal_distribution[ivariable]->nb_value - 1 : marginal_distribution[ivariable]->nb_value))) {
          status = false;
          ostringstream error_message;
          error_message << STAT_label[STATL_CATEGORY] << " " << category[i] << " "
                        << STAT_error[STATR_NOT_ALLOWED];
          error.update((error_message.str()).c_str());
        }

        else {
          if (category[i] < min_category) {
            min_category = category[i];
          }
          if (category[i] > max_category) {
            max_category = category[i];
          }
        }
      }

      if ((min_category != 0) || (max_category == 0)) {
        status = false;
        error.update(STAT_error[STATR_NB_CATEGORY]);
      }
    }

    if (status) {
      presence = new bool[max_category + 1];
      for (i = 0;i <= max_category;i++) {
        presence[i] = false;
      }

      for (i = 0;i < marginal_distribution[ivariable]->nb_value - marginal_distribution[ivariable]->offset;i++) {
        presence[category[i]] = true;
      }

      for (i = 0;i <= max_category;i++) {
        if (!presence[i]) {
          status = false;
          ostringstream error_message;
          error_message << STAT_error[STATR_MISSING_CATEGORY] << " " << i;
          error.update((error_message.str()).c_str());
        }
      }

      delete [] presence;
    }

    if (status) {
      if (add_variable) {
        variable = 0;
        offset = 1;
      }
      else {
        variable = ivariable;
        offset = 0;
      }

      itype = new variable_nature[nb_variable + offset];
      for (i = 0;i < nb_variable;i++) {
        itype[i + offset] = type[i];
      }
      itype[variable] = INT_VALUE;

      seq = new MarkovianSequences(nb_sequence , identifier , length , vertex_identifier ,
                                   index_param_type , nb_variable + offset , itype);
      delete [] itype;

      seq->Sequences::transcode(*this , ivariable , 0 , max_category , category , add_variable);

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


/*--------------------------------------------------------------*/
/**
 *  \brief Transcoding of categories of an integer-valued variable.
 *
 *  \param[in] error        reference on a StatError object,
 *  \param[in] ivariable    variable index,
 *  \param[in] category     transcoding table,
 *  \param[in] add_variable flag for adding a variable.
 *
 *  \return                 MarkovianSequences object.
 */
/*--------------------------------------------------------------*/

MarkovianSequences* MarkovianSequences::transcode(StatError &error , int ivariable ,
                                                  vector<int> category , bool add_variable) const

{
  return transcode(error , ivariable , category.data() , add_variable);
}


/*--------------------------------------------------------------*/
/**
 *  \brief Transcoding of categories of an integer-valued variable.
 *
 *  \param[in] error   reference on a StatError object,
 *  \param[in] process reference on a CategoricalSequenceProcess object.
 *
 *  \return            MarkovianSequences object.
 */
/*--------------------------------------------------------------*/

MarkovianSequences* MarkovianSequences::transcode(StatError &error ,
                                                  const CategoricalSequenceProcess *process) const

{
  int i , j;
  int *category;
  MarkovianSequences *seq;


  category = new int[process->nb_value];
  for (i = 0;i < process->nb_state;i++) {
    for (j = 0;j < process->nb_value;j++) {
      if (process->observation[i]->mass[j] > 0.) {
        category[j] = i;
      }
    }
  }

  seq = transcode(error , 1 , category , true);
  delete [] category;

  return seq;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Removing of the non-represented values of an integer-valued variable.
 *
 *  \param[in] error        reference on a StatError object,
 *  \param[in] display      flag for displaying the non-represented values,
 *  \param[in] ivariable    variable index,
 *  \param[in] add_variable flag for adding a variable.
 *
 *  \return                 MarkovianSequences object.
 */
/*--------------------------------------------------------------*/

MarkovianSequences* MarkovianSequences::consecutive_values(StatError &error , bool display ,
                                                           int ivariable , bool add_variable) const

{
  bool status = true;
  int i , j;
  int variable , offset , max , *category;
  variable_nature *itype;
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
    if (display) {
      cout << "\n" << SEQ_label[SEQL_MISSING_VALUE] << ":";
      for (i = 0;i < marginal_distribution[ivariable]->nb_value;i++) {
        if (marginal_distribution[ivariable]->frequency[i] == 0) {
          cout << " " << i;
        }
      }
      cout << endl;
    }

    category = new int[marginal_distribution[ivariable]->nb_value - marginal_distribution[ivariable]->offset];

//    i = 0;
    i = -1;
    for (j = marginal_distribution[ivariable]->offset;j < marginal_distribution[ivariable]->nb_value;j++) {
//      category[j - marginal_distribution[ivariable]->offset] = i;
      if (marginal_distribution[ivariable]->frequency[j] > 0) {
        i++;
      }
      category[j - marginal_distribution[ivariable]->offset] = i;
    }
//    max = i - 1;
    max = i;

#   ifdef DEBUG
    cout << "\nTest :";
    for (i = 0;i < marginal_distribution[ivariable]->nb_value - marginal_distribution[ivariable]->offset;i++) {
      cout << " " << category[i];
    }
    cout << endl;
#   endif

    if (add_variable) {
      variable = 0;
      offset = 1;
    }
    else {
      variable = ivariable;
      offset = 0;
    }

    itype = new variable_nature[nb_variable + offset];
    for (i = 0;i < nb_variable;i++) {
      itype[i + offset] = type[i];
    }
    itype[variable] = INT_VALUE;

    seq = new MarkovianSequences(nb_sequence , identifier , length , vertex_identifier ,
                                 index_param_type , nb_variable + offset , itype);
    delete [] itype;

    seq->Sequences::transcode(*this , ivariable , 0 , max , category , add_variable);
    delete [] category;

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


/*--------------------------------------------------------------*/
/**
 *  \brief Partitioning of values of an integer-valued variable.
 *
 *  \param[in] error        reference on a StatError object,
 *  \param[in] ivariable    variable index,
 *  \param[in] nb_class     number of classes,
 *  \param[in] ilimit       integer limits between classes (beginning of classes),
 *  \param[in] add_variable flag for adding a variable.
 *
 *  \return                 MarkovianSequences object.
 */
/*--------------------------------------------------------------*/

MarkovianSequences* MarkovianSequences::cluster(StatError &error , int ivariable , int nb_class ,
                                                int *ilimit , bool add_variable) const

{
  bool status = true;
  int i , j , k;
  int variable , offset , *category , *limit;
  variable_nature *itype;
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
      category = new int[marginal_distribution[ivariable]->nb_value];

      i = 0;
      for (j = 0;j < nb_class;j++) {
        for (k = limit[j];k < limit[j + 1];k++) {
          category[i++] = j;
        }
      }

      if (add_variable) {
        variable = 0;
        offset = 1;
      }
      else {
        variable = ivariable;
        offset = 0;
      }

      itype = new variable_nature[nb_variable + offset];
      for (i = 0;i < nb_variable;i++) {
        itype[i + offset] = type[i];
      }
      itype[variable] = INT_VALUE;

      seq = new MarkovianSequences(nb_sequence , identifier , length , vertex_identifier ,
                                   index_param_type , nb_variable + offset , itype);
      delete [] itype;

      seq->Sequences::transcode(*this , ivariable , 0 , nb_class - 1 , category , add_variable);
      delete [] category;

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


/*--------------------------------------------------------------*/
/**
 *  \brief Partitioning of values of an integer-valued variable.
 *
 *  \param[in] error        reference on a StatError object,
 *  \param[in] ivariable    variable index,
 *  \param[in] nb_class     number of classes,
 *  \param[in] ilimit       integer limits between classes (beginning of classes),
 *  \param[in] add_variable flag for adding a variable.
 *
 *  \return                 MarkovianSequences object.
 */
/*--------------------------------------------------------------*/

MarkovianSequences* MarkovianSequences::cluster(StatError &error , int ivariable , int nb_class ,
                                                vector<int> ilimit , bool add_variable) const

{
  return cluster(error , ivariable , nb_class , ilimit.data() , add_variable);
}


/*--------------------------------------------------------------*/
/**
 *  \brief Partitioning of values of a real-valued variable.
 *
 *  \param[in] error    reference on a StatError object,
 *  \param[in] variable variable index,
 *  \param[in] nb_class number of classes,
 *  \param[in] ilimit   real limits between classes (beginning of classes).
 *
 *  \return             MarkovianSequences object.
 */
/*--------------------------------------------------------------*/

MarkovianSequences* MarkovianSequences::cluster(StatError &error , int variable ,
                                                int nb_class , double *ilimit) const

{
  bool status = true;
  int i;
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


/*--------------------------------------------------------------*/
/**
 *  \brief Partitioning of values of a real-valued variable.
 *
 *  \param[in] error    reference on a StatError object,
 *  \param[in] variable variable index,
 *  \param[in] nb_class number of classes,
 *  \param[in] ilimit   real limits between classes (beginning of classes).
 *
 *  \return             MarkovianSequences object.
 */
/*--------------------------------------------------------------*/

MarkovianSequences* MarkovianSequences::cluster(StatError &error , int variable ,
                                                int nb_class , vector<double> ilimit) const

{
  return cluster(error , variable , nb_class , ilimit.data());
}


/*--------------------------------------------------------------*/
/**
 *  \brief Copy of a MarkovianSequences object transforming the implicit index parameters in
 *         explicit index parameters.
 *
 *  \param[in] error reference on a StatError object.
 *
 *  \return          MarkovianSequences object.
 */
/*--------------------------------------------------------------*/

MarkovianSequences* MarkovianSequences::explicit_index_parameter(StatError &error) const

{
  MarkovianSequences *seq;


  error.init();

  if (index_parameter) {
    seq = NULL;
    error.update(SEQ_error[SEQR_INDEX_PARAMETER_TYPE]);
  }
  else {
    seq = new MarkovianSequences(*this , EXPLICIT_INDEX_PARAMETER);
  }

  return seq;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Removing of the index parameters.
 *
 *  \param[in] error reference on a StatError object.
 *
 *  \return          MarkovianSequences object.
 */
/*--------------------------------------------------------------*/

MarkovianSequences* MarkovianSequences::remove_index_parameter(StatError &error) const

{
  MarkovianSequences *seq;


  error.init();

  if (!index_parameter) {
    seq = NULL;
    error.update(SEQ_error[SEQR_INDEX_PARAMETER_TYPE]);
  }
  else {
    seq = new MarkovianSequences(*this , REMOVE_INDEX_PARAMETER);
  }

  return seq;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Selection of variables.
 *
 *  \param[in] error        reference on a StatError object,
 *  \param[in] inb_variable number of variables,
 *  \param[in] ivariable    variable indices,
 *  \param[in] keep         flag for keeping or rejecting the selected variables.
 *
 *  \return                 MarkovianSequences object.
 */
/*--------------------------------------------------------------*/

MarkovianSequences* MarkovianSequences::select_variable(StatError &error , int inb_variable ,
                                                        int *ivariable , bool keep) const

{
  bool status = true , *selected_variable;
  int i;
  int bnb_variable , *variable;
  variable_nature *itype;
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
      itype = new variable_nature[bnb_variable];
      for (i = 0;i < bnb_variable;i++) {
        itype[i] = type[variable[i]];
      }

      seq = new MarkovianSequences(nb_sequence , identifier , length , vertex_identifier ,
                                   index_param_type , bnb_variable , itype);

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


/*--------------------------------------------------------------*/
/**
 *  \brief Selection of variables.
 *
 *  \param[in] error        reference on a StatError object,
 *  \param[in] inb_variable number of variables,
 *  \param[in] ivariable    variable indices,
 *  \param[in] keep         flag for keeping or rejecting the selected variables.
 *
 *  \return                 MarkovianSequences object.
 */
/*--------------------------------------------------------------*/

MarkovianSequences* MarkovianSequences::select_variable(StatError &error , int inb_variable ,
                                                        vector<int> ivariable , bool keep) const

{
  return select_variable(error , inb_variable , ivariable.data() , keep);
}


/*--------------------------------------------------------------*/
/**
 *  \brief Removing of the 1st variable.
 *
 *  \return MarkovianSequences object.
 */
/*--------------------------------------------------------------*/

MarkovianSequences* MarkovianSequences::remove_variable_1() const

{
  int i;
  int *variable;
  variable_nature *itype;
  MarkovianSequences *seq;


  variable = new int[nb_variable - 1];
  itype = new variable_nature[nb_variable - 1];
  for (i = 0;i < nb_variable - 1;i++) {
    variable[i] = i + 1;
    itype[i] = type[i + 1];
  }

  seq = new MarkovianSequences(nb_sequence , identifier , length , vertex_identifier ,
                               index_param_type , nb_variable - 1 , itype);

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


/*--------------------------------------------------------------*/
/**
 *  \brief Merging of variables of MarkovianSequences objects.
 *
 *  \param[in] error      reference on a StatError object,
 *  \param[in] nb_sample  number of MarkovianSequences objects,
 *  \param[in] iseq       pointer on the MarkovianSequences objects,
 *  \param[in] ref_sample reference MarkovianSequences object for the identifiers.
 *
 *  \return               MarkovianSequences object.
 */
/*--------------------------------------------------------------*/

MarkovianSequences* MarkovianSequences::merge_variable(StatError &error , int nb_sample ,
                                                       const MarkovianSequences **iseq , int ref_sample) const

{
  bool status = true;
  int i , j , k , m;
  int inb_variable , *iidentifier , **ivertex_identifier;
  variable_nature *itype;
  MarkovianSequences *seq;
  const MarkovianSequences **pseq;


  seq = NULL;
  error.init();

  for (i = 0;i < nb_sample;i++) {
    if ((iseq[i]->index_param_type != IMPLICIT_TYPE) &&
        (iseq[i]->index_param_type != index_param_type)) {
      status = false;
      ostringstream error_message;
      error_message << STAT_label[STATL_SAMPLE] << " " << i + 2 << ": "
                    << SEQ_error[SEQR_INDEX_PARAMETER_TYPE];

      if (index_param_type == IMPLICIT_TYPE) {
        error.update((error_message.str()).c_str());
      }
      else {
        error.correction_update((error_message.str()).c_str() , SEQ_index_parameter_word[index_param_type]);
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

        else if ((iseq[i]->index_param_type == TIME) &&
                 (iseq[i]->index_param_type == index_param_type)) {
          for (k = 0;k < length[j];k++) {
            if (iseq[i]->index_parameter[j][k] != index_parameter[j][k]) {
              status = false;
              ostringstream error_message;
              error_message << STAT_label[STATL_SAMPLE] << " " << i + 2 << ": "
                            << SEQ_label[SEQL_SEQUENCE] << " " << j + 1 << ": "
                            << SEQ_label[SEQL_INDEX] << " " << k << ": "
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

    // comparison of sequence identifiers

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

    // comparison of vertex identifiers

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

    itype = new variable_nature[inb_variable];
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
                                 index_param_type , inb_variable , itype);
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

    // copy of values

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


/*--------------------------------------------------------------*/
/**
 *  \brief Merging of variables of MarkovianSequences objects.
 *
 *  \param[in] error      reference on a StatError object,
 *  \param[in] nb_sample  number of MarkovianSequences objects,
 *  \param[in] iseq       pointer on the MarkovianSequences objects,
 *  \param[in] ref_sample reference MarkovianSequences object for the identifiers.
 *
 *  \return               MarkovianSequences object.
 */
/*--------------------------------------------------------------*/

MarkovianSequences* MarkovianSequences::merge_variable(StatError &error , int nb_sample ,
                                                       const vector<MarkovianSequences> iseq , int ref_sample) const

{
  int i;
  MarkovianSequences *seq;
  const MarkovianSequences **pseq;


  pseq = new const MarkovianSequences*[nb_sample];
  for (i = 0;i < nb_sample;i++) {
    pseq[i] = new MarkovianSequences(iseq[i]);
  }

  seq = merge_variable(error , nb_sample , pseq , ref_sample);

  for (i = 0;i < nb_sample;i++) {
    delete pseq[i];
  }
  delete [] pseq;

  return seq;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Copy of sequences with extraction of the initial runs.
 *
 *  \param[in] error reference on a StatError object.
 *
 *  \return          MarkovianSequences object.
 */
/*--------------------------------------------------------------*/

 MarkovianSequences* MarkovianSequences::initial_run_computation(StatError &error) const

{
  int i;
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
    seq = new MarkovianSequences(*this , SEQUENCE_COPY , ADD_INITIAL_RUN);
  }

  return seq;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Addition of an absorbing run at the end of each sequence.
 *
 *  \param[in] error           reference on a StatError object,
 *  \param[in] run_length      absorbing run length,
 *  \param[in] sequence_length sequence length,
 *  \param[in] add_variable    flag for adding a binary variable (0: data, 1: end absorbing run).
 *
 *  \return                    MarkovianSequences object.
 */
/*--------------------------------------------------------------*/

MarkovianSequences* MarkovianSequences::add_absorbing_run(StatError &error , int run_length ,
                                                          int sequence_length , bool add_variable) const

{
  bool status = true , initial_run_flag;
  int i , j , k;
  int inb_variable , end_value , *ilength;
  double mean , variance , limit , *standard_deviation;
  variable_nature *itype;
  MarkovianSequences *seq;


  seq = NULL;
  error.init();

/*  if (index_param_type == TIME) {
    status = false;
    error.update(SEQ_error[SEQR_INDEX_PARAMETER_TYPE]);
  } */

  if ((run_length != I_DEFAULT) && ((run_length < 1) ||
       (run_length > MAX_ABSORBING_RUN_LENGTH))) {
    status = false;
    error.update(SEQ_error[SEQR_RUN_LENGTH]);
  }

  if ((sequence_length != I_DEFAULT) && ((sequence_length <= max_length) ||
       (sequence_length > max_length + MAX_ABSORBING_RUN_LENGTH))) {
    status = false;
    error.update(SEQ_error[SEQR_SEQUENCE_LENGTH]);
  }

  if (status) {
    if (add_variable) {
      inb_variable = nb_variable + 1;
    }
    else {
      inb_variable = nb_variable;
    }

    itype = new variable_nature[inb_variable];
    for (i = 0;i < nb_variable;i++) {
      itype[i] = type[i];
    }
    if (add_variable) {
      itype[nb_variable] = INT_VALUE;
    }

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
                                 index_param_type , inb_variable , itype , false);
    delete [] itype;
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

    standard_deviation = new double[nb_variable];

    for (i = 0;i < nb_variable;i++) {
      if ((type[i] == REAL_VALUE) || (type[i] == AUXILIARY)) {
        mean = mean_computation(i);
        variance = variance_computation(i , mean);
        standard_deviation[i] = sqrt(variance) / ABSORBING_RUN_STANDARD_DEVIATION_FACTOR;
      }
    }

    // copy of sequences with addition of an end absorbing run

    for (i = 0;i < seq->nb_sequence;i++) {
      for (j = 0;j < nb_variable;j++) {
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
          if (run_length == 1) {  // for Fuji/Braeburn GUs
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

          // random generation of absorbing run values

          if (min_value[j] >= 10 * standard_deviation[j]) {
            normal dist(min_value[j] - 2 * standard_deviation[j] , standard_deviation[j]);

            for (k = length[i];k < seq->length[i];k++) {
              limit = ((double)rand() / (RAND_MAX + 1.));
              seq->real_sequence[i][j][k] = quantile(dist , limit);
            }
          }

          else {
            normal dist(max_value[j] + 2 * standard_deviation[j] , standard_deviation[j]);

            for (k = length[i];k < seq->length[i];k++) {
              limit = ((double)rand() / (RAND_MAX + 1.));
              seq->real_sequence[i][j][k] = quantile(dist , limit);
            }
          }

/*          if (min_value[j] >= 5 * standard_deviation[j]) {
            for (k = length[i];k < seq->length[i];k++) {
              seq->real_sequence[i][j][k] = min_value[j] - (k % 2 + 4) * standard_deviation[j];
            }
          }

          else {
            for (k = length[i];k < seq->length[i];k++) {
              seq->real_sequence[i][j][k] = max_value[j] + (k % 2 + 4) * standard_deviation[j];
            }
          } */
        }
      }

      // addition of a binary variable

      if (add_variable) {
        for (j = 0;j < length[i];j++) {
          seq->int_sequence[i][nb_variable][j] = 0;
        }
        for (j = length[i];j < seq->length[i];j++) {
          seq->int_sequence[i][nb_variable][j] = 1;
        }
      }
    }

    delete [] standard_deviation;

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


/*--------------------------------------------------------------*/
/**
 *  \brief Construction of auxiliary variables corresponding to
 *         restored state sequences.
 *
 *  \param[in] discrete_process   pointer on DiscreteParametricProcess objects,
 *  \param[in] continuous_process pointer on ContinuousParametricProcess objects.
 *
 *  \return                       MarkovianSequences object.
 */
/*--------------------------------------------------------------*/

MarkovianSequences* MarkovianSequences::build_auxiliary_variable(DiscreteParametricProcess **discrete_process ,
                                                                 ContinuousParametricProcess **continuous_process) const

{
  bool *auxiliary;
  int i , j , k , m;
  int *pstate;
  double *mean;
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
        for (m = 0;m < length[k];m++) {
          seq->real_sequence[k][i][m] = discrete_process[j - 1]->observation[*pstate++]->mean;
        }
      }
    }

    else if ((continuous_process) && (continuous_process[j - 1])) {
      i++;

      if ((continuous_process[j - 1]->ident == GAMMA) || (continuous_process[j - 1]->ident == ZERO_INFLATED_GAMMA) ||
          (continuous_process[j - 1]->ident == GAUSSIAN) || (continuous_process[j - 1]->ident == INVERSE_GAUSSIAN) ||
          (continuous_process[j - 1]->ident == VON_MISES)) {
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

        default : {
          for (k = 0;k < continuous_process[j - 1]->nb_state;k++) {
            mean[k] = continuous_process[j - 1]->observation[k]->location;
          }
          break;
        }
        }

        for (k = 0;k < nb_sequence;k++) {
          pstate = seq->int_sequence[k][0];
          for (m = 0;m < length[k];m++) {
            seq->real_sequence[k][i][m] = mean[*pstate++];
          }
        }

        delete [] mean;
      }

      else if (continuous_process[j - 1]->ident == LINEAR_MODEL) {
        switch (index_param_type) {

        case IMPLICIT_TYPE : {
          for (k = 0;k < nb_sequence;k++) {
            pstate = seq->int_sequence[k][0];
            for (m = 0;m < length[k];m++) {
              seq->real_sequence[k][i][m] = continuous_process[j - 1]->observation[*pstate]->intercept +
                                            continuous_process[j - 1]->observation[*pstate]->slope * m;
              pstate++;
            }
          }
          break;
        }

        case TIME : {
          for (k = 0;k < nb_sequence;k++) {
            pstate = seq->int_sequence[k][0];
            for (m = 0;m < length[k];m++) {
              seq->real_sequence[k][i][m] = continuous_process[j - 1]->observation[*pstate]->intercept +
                                            continuous_process[j - 1]->observation[*pstate]->slope * index_parameter[k][m];
              pstate++;
            }
          }
          break;
        }
        }
      }

      else if (continuous_process[j - 1]->ident == AUTOREGRESSIVE_MODEL) {
        switch (type[i - 1]) {

        case INT_VALUE : {
          for (k = 0;k < nb_sequence;k++) {
            pstate = seq->int_sequence[k][0];
            for (m = 0;m < length[k];m++) {
              if ((m == 0) || (*pstate != *(pstate - 1))) {
                seq->real_sequence[k][i][m] = continuous_process[j - 1]->observation[*pstate]->location;
              }
              else {
                seq->real_sequence[k][i][m] = continuous_process[j - 1]->observation[*pstate]->location +
                                              continuous_process[j - 1]->observation[*pstate]->autoregressive_coeff *
                                              (seq->int_sequence[k][i - 1][m - 1] - continuous_process[j - 1]->observation[*pstate]->location);
              }
              pstate++;
            }
          }
          break;
        }

        case REAL_VALUE : {
          for (k = 0;k < nb_sequence;k++) {
            pstate = seq->int_sequence[k][0];
            for (m = 0;m < length[k];m++) {
              if ((m == 0) || (*pstate != *(pstate - 1))) {
                seq->real_sequence[k][i][m] = continuous_process[j - 1]->observation[*pstate]->location;
              }
              else {
                seq->real_sequence[k][i][m] = continuous_process[j - 1]->observation[*pstate]->location +
                                              continuous_process[j - 1]->observation[*pstate]->autoregressive_coeff *
                                              (seq->real_sequence[k][i - 1][m - 1] - continuous_process[j - 1]->observation[*pstate]->location);
              }
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


/*--------------------------------------------------------------*/
/**
 *  \brief Building of residual sequences on the basis of restored state sequences.
 *
 *  \param[in] categorical_process pointer on CategoricalSequenceProcess objects,
 *  \param[in] discrete_process    pointer on DiscreteParametricProcess objects,
 *  \param[in] continuous_process  pointer on ContinuousParametricProcess objects.
 *
 *  \return                        MarkovianSequences object.
 */
/*--------------------------------------------------------------*/

MarkovianSequences* MarkovianSequences::residual_sequences(CategoricalSequenceProcess **categorical_process ,
                                                           DiscreteParametricProcess **discrete_process ,
                                                           ContinuousParametricProcess **continuous_process) const

{
  int i , j , k;
  int *pstate;
  double *mean;
  variable_nature *itype;
  MarkovianSequences *seq;


  itype = new variable_nature[nb_variable];
  itype[0] = type[0];
  for (i = 1;i < nb_variable;i++) {
    itype[i] = REAL_VALUE;
  }

  seq = new MarkovianSequences(nb_sequence , identifier , length , vertex_identifier ,
                               index_param_type , nb_variable , itype);
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
      for (j = 0;j < (index_param_type == POSITION ? length[i] + 1 : length[i]);j++) {
        seq->index_parameter[i][j] = index_parameter[i][j];
      }
    }
  }

  // copy of restored state sequences

  for (i = 0;i < nb_sequence;i++) {
    for (j = 0;j < length[i];j++) {
      seq->int_sequence[i][0][j] = int_sequence[i][0][j];
    }
  }

  seq->min_value[0] = min_value[0];
  seq->max_value[0] = max_value[0];
  seq->marginal_distribution[0] = new FrequencyDistribution(*marginal_distribution[0]);

  // computation of residual sequences

  for (i = 1;i < nb_variable;i++) {
    if ((categorical_process) && (categorical_process[i - 1])) {
      for (j = 0;j < nb_sequence;j++) {
        pstate = seq->int_sequence[j][0];
        for (k = 0;k < length[j];k++) {
          seq->real_sequence[j][i][k] = int_sequence[j][i][k] - categorical_process[i - 1]->observation[*pstate++]->mean;
        }
      }
    }

    else if ((discrete_process) && (discrete_process[i - 1])) {
      for (j = 0;j < nb_sequence;j++) {
        pstate = seq->int_sequence[j][0];
        for (k = 0;k < length[j];k++) {
          seq->real_sequence[j][i][k] = int_sequence[j][i][k] - discrete_process[i - 1]->observation[*pstate++]->mean;
        }
      }
    }

    else if ((continuous_process) && (continuous_process[i - 1])) {
      if ((continuous_process[i - 1]->ident == GAMMA) || (continuous_process[i - 1]->ident == ZERO_INFLATED_GAMMA) ||
          (continuous_process[i - 1]->ident == GAUSSIAN) || (continuous_process[i - 1]->ident == INVERSE_GAUSSIAN) ||
          (continuous_process[i - 1]->ident == VON_MISES)) {
        mean = new double [continuous_process[i - 1]->nb_state];

        switch (continuous_process[i - 1]->ident) {

        case GAMMA : {
          for (j = 0;j < continuous_process[i - 1]->nb_state;j++) {
            mean[j] = continuous_process[i - 1]->observation[j]->shape *
                      continuous_process[i - 1]->observation[j]->scale;
          }
          break;
        }

        case ZERO_INFLATED_GAMMA : {
          for (j = 0;j < continuous_process[i - 1]->nb_state;j++) {
            if (continuous_process[i - 1]->observation[j]->zero_probability == 1.) {
              mean[j] = 0.;
            }
            else {
              mean[j] = (1 - continuous_process[i - 1]->observation[j]->zero_probability) *
                        continuous_process[i - 1]->observation[j]->shape *
                        continuous_process[i - 1]->observation[j]->scale;
            }
          }
          break;
        }

        default : {
          for (j = 0;j < continuous_process[i - 1]->nb_state;j++) {
            mean[j] = continuous_process[i - 1]->observation[j]->location;
          }
          break;
        }
        }

        switch (type[i]) {

        case INT_VALUE : {
          for (j = 0;j < nb_sequence;j++) {
            pstate = seq->int_sequence[j][0];
            for (k = 0;k < length[j];k++) {
              seq->real_sequence[j][i][k] = int_sequence[j][i][k] - mean[*pstate++];
            }
          }
          break;
        }

        case REAL_VALUE : {
          for (j = 0;j < nb_sequence;j++) {
            pstate = seq->int_sequence[j][0];
            for (k = 0;k < length[j];k++) {
              seq->real_sequence[j][i][k] = real_sequence[j][i][k] - mean[*pstate++];
            }
          }
          break;
        }
        }

        delete [] mean;
      }

      else if (continuous_process[i - 1]->ident == LINEAR_MODEL) {
        switch (index_param_type) {

        case IMPLICIT_TYPE : {
          switch (type[i]) {

          case INT_VALUE : {
            for (j = 0;j < nb_sequence;j++) {
              pstate = seq->int_sequence[j][0];

              for (k = 0;k < length[j];k++) {
                seq->real_sequence[j][i][k] = int_sequence[j][i][k] - (continuous_process[i - 1]->observation[*pstate]->intercept +
                                               continuous_process[i - 1]->observation[*pstate]->slope * k);
                pstate++;
              }
            }
            break;
          }

          case REAL_VALUE : {
            for (j = 0;j < nb_sequence;j++) {
              pstate = seq->int_sequence[j][0];

              for (k = 0;k < length[j];k++) {
                seq->real_sequence[j][i][k] = real_sequence[j][i][k] - (continuous_process[i - 1]->observation[*pstate]->intercept +
                                               continuous_process[i - 1]->observation[*pstate]->slope * k);
                pstate++;
              }
            }
            break;
          }
          }
          break;
        }

        case TIME : {
          switch (type[i]) {

          case INT_VALUE : {
            for (j = 0;j < nb_sequence;j++) {
              pstate = seq->int_sequence[j][0];

              for (k = 0;k < length[j];k++) {
                seq->real_sequence[j][i][k] = int_sequence[j][i][k] - (continuous_process[i - 1]->observation[*pstate]->intercept +
                                               continuous_process[i - 1]->observation[*pstate]->slope * index_parameter[j][k]);
                pstate++;
              }
            }
            break;
          }

          case REAL_VALUE : {
            for (j = 0;j < nb_sequence;j++) {
              pstate = seq->int_sequence[j][0];

              for (k = 0;k < length[j];k++) {
                seq->real_sequence[j][i][k] = real_sequence[j][i][k] - (continuous_process[i - 1]->observation[*pstate]->intercept +
                                               continuous_process[i - 1]->observation[*pstate]->slope * index_parameter[j][k]);
                pstate++;
              }
            }
            break;
          }
          }
          break;
        }
        }
      }

      else if (continuous_process[i - 1]->ident == AUTOREGRESSIVE_MODEL) {
        switch (type[i]) {

        case INT_VALUE : {
          for (j = 0;j < nb_sequence;j++) {
            pstate = seq->int_sequence[j][0];

            for (k = 0;k < length[j];k++) {
              if ((k == 0) || (*pstate != *(pstate - 1))) {
                seq->real_sequence[j][i][k] = int_sequence[j][i][k] - continuous_process[i - 1]->observation[*pstate]->location;
              }
              else {
                seq->real_sequence[j][i][k] = int_sequence[j][i][k] - (continuous_process[i - 1]->observation[*pstate]->location +
                                               continuous_process[i - 1]->observation[*pstate]->autoregressive_coeff *
                                               (int_sequence[j][i][k - 1] - continuous_process[i - 1]->observation[*pstate]->location));
              }
              pstate++;
            }
          }
          break;
        }

        case REAL_VALUE : {
          for (j = 0;j < nb_sequence;j++) {
            pstate = seq->int_sequence[j][0];

            for (k = 0;k < length[j];k++) {
              if ((k == 0) || (*pstate != *(pstate - 1))) {
                seq->real_sequence[j][i][k] = real_sequence[j][i][k] - continuous_process[i - 1]->observation[*pstate]->location;
              }
              else {
                seq->real_sequence[j][i][k] = real_sequence[j][i][k] - (continuous_process[i - 1]->observation[*pstate]->location +
                                               continuous_process[i - 1]->observation[*pstate]->autoregressive_coeff * 
                                               (real_sequence[j][i][k - 1] - continuous_process[i - 1]->observation[*pstate]->location));
              }
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
    seq->min_value_computation(i);
    seq->max_value_computation(i);

    seq->build_marginal_histogram(i);
  }

  return seq;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Split of sequences in equal length segments.
 *
 *  \param[in] error reference on a StatError object,
 *  \param[in] step  sequence lengths.
 *
 *  \return          MarkovianSequences object.
 */
/*--------------------------------------------------------------*/

MarkovianSequences* MarkovianSequences::split(StatError &error , int step) const

{
  int i , j , k , m;
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
                                 index_param_type , nb_variable , type);
    delete [] ilength;

    // copy of sequences

    inb_sequence = 0;
    for (i = 0;i < nb_sequence;i++) {
      nb_segment = (length[i] % step == 0 ? length[i] / step : length[i] / step + 1);

      if (seq->index_param_type == TIME) {
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

    if (seq->index_param_type == TIME) {
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


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the cumulative frequency distribution function for a variable.
 *
 *  \param[in] variable variable index,
 *  \param[in] cdf      (value, cumulative distribution function).
 *
 *  \return             number of values.
 */
/*--------------------------------------------------------------*/

int MarkovianSequences::cumulative_distribution_function_computation(int variable , double **cdf) const

{
  int i , j , k;
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

        // search for the current minimum value

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

        // determination of the number of vectors taken the current minimum value
        // for the selected variable

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

        // search for the current minimum value

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

        // determination of the number of vectors taken the current minimum value
        // for the selected variable

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


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the cumulative frequency distribution function for a variable.
 *
 *  \param[in] variable variable index,
 *  \param[in] state    state,
 *  \param[in] cdf      (value, cumulative distribution function).
 *
 *  \return             number of values.
 */
/*--------------------------------------------------------------*/

int MarkovianSequences::cumulative_distribution_function_computation(int variable , int state ,
                                                                     double **cdf) const

{
  int i , j , k;
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

        // search for the current minimum value

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

        // determination of the number of vectors taken the current minimum value
        // for the selected variable

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

        // search for the current minimum value

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

        // determination of the number of vectors taken the current minimum value
        // for the selected variable

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


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the minimum interval between 2 values for a variable.
 *
 *  \param[in] variable variable index.
 */
/*--------------------------------------------------------------*/

void MarkovianSequences::min_interval_computation(int variable)

{
  if (marginal_distribution[variable]) {
    min_interval[variable] = marginal_distribution[variable]->min_interval_computation();
  }

  else if (type[variable] != AUXILIARY) {
    bool *selected_value;
    int i , j , k , m;
    int int_min , int_value , nb_value , max_frequency , *frequency , *index;
    double real_min , real_value;


    min_interval[variable] = max_value[variable] - min_value[variable];
    i = 0;

    switch (type[variable]) {

    case INT_VALUE : {
      do {

        // search for the current minimum value

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

        // determination of the number of vectors taken the current minimum value
        // for the selected variable

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

        // search for the current minimum value

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

        // determination of the number of vectors taken the current minimum value
        // for the selected variable

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

      // search for the median frequency

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

#     ifdef DEBUG
      cout << "\n" << STAT_label[STATL_VARIABLE] << " " << variable + 1 << ": "
           << min_interval[variable] << " " << frequency[index[nb_value / 2]];
#     endif

      // to be finalized

      if (frequency[index[nb_value / 2]] == 1) {
        min_interval[variable] = 0.;
      }

#     ifdef DEBUG
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


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the information quantity assuming i.i.d. variables.
 */
/*--------------------------------------------------------------*/

double MarkovianSequences::iid_information_computation() const

{
  int i;
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


/*--------------------------------------------------------------*/
/**
 *  \brief Extraction of samples corresponding to the changes in self-transition probability
 *         for a state of a non-homogeneous Markov chain.
 *
 *  \param[in] state state.
 */
/*--------------------------------------------------------------*/

void MarkovianSequences::self_transition_computation(int state)

{
  int i , j;
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


/*--------------------------------------------------------------*/
/**
 *  \brief Extraction of samples corresponding to the changes in self-transition probability
 *         for states of a non-homogeneous Markov chain.
 */
/*--------------------------------------------------------------*/

void MarkovianSequences::self_transition_computation()

{
  if (!self_transition) {
    int i;


    state_variable_init();
    self_transition = new SelfTransition*[marginal_distribution[0]->nb_value];

    for (i = 0;i < marginal_distribution[0]->nb_value;i++) {
      self_transition[i] = new SelfTransition(max_length - 1);
      self_transition_computation(i);
    }
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Extraction of samples corresponding to the changes in self-transition probability
 *         for states of a non-homogeneous Markov chain.
 *
 *  \param[in] homogeneity state homogeneities.
 */
/*--------------------------------------------------------------*/

void MarkovianSequences::self_transition_computation(bool *homogeneity)

{
  if (!self_transition) {
    int i;


    state_variable_init();
    self_transition = new SelfTransition*[marginal_distribution[0]->nb_value];

    for (i = 0;i < marginal_distribution[0]->nb_value;i++) {
      if (homogeneity[i]) {
        self_transition[i] = NULL;
      }
      else {
        self_transition[i] = new SelfTransition(max_length - 1);
        self_transition_computation(i);
      }
    }
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the marginal state frequency distribution from
 *         the restored state sequences.
 *
 *  \return Distribution object.
 */
/*--------------------------------------------------------------*/

Distribution* MarkovianSequences::weight_computation() const

{
  int i;
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


/*--------------------------------------------------------------*/
/**
 *  \brief Update of the observation frequency distributions for a variable.
 *
 *  \param[in] variable variable index,
 *  \param[in] nb_state number of states.
 */
/*--------------------------------------------------------------*/

void MarkovianSequences::observation_frequency_distribution_computation(int variable ,
                                                                        int nb_state)

{
  int i , j;
  int *pstate , *poutput;


  // initialization of the frequency distributions

  for (i = 0;i < nb_state;i++) {
    for (j = 0;j < marginal_distribution[variable]->nb_value;j++) {
      observation_distribution[variable][i]->frequency[j] = 0;
    }
  }

  // update of the frequency distributions

  for (i = 0;i < nb_sequence;i++) {
    pstate = int_sequence[i][0];
    poutput = int_sequence[i][variable];
    for (j = 0;j < length[i];j++) {
      (observation_distribution[variable][*pstate++]->frequency[*poutput++])++;
    }
  }

  // computation of the frequency distribution characteristics

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


/*--------------------------------------------------------------*/
/**
 *  \brief Construction of the observation frequency distributions.
 *
 *  \param[in] nb_state number of states.
 */
/*--------------------------------------------------------------*/

void MarkovianSequences::build_observation_frequency_distribution(int nb_state)

{
  if ((nb_variable > 1) && (!observation_distribution)) {
    int i , j;


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


/*--------------------------------------------------------------*/
/**
 *  \brief Construction of the observation histograms for a variable.
 *
 *  \param[in] variable  variable index,
 *  \param[in] nb_state  number of states, 
 *  \param[in] bin_width bin width.
 */
/*--------------------------------------------------------------*/

void MarkovianSequences::build_observation_histogram(int variable , int nb_state , double bin_width)

{
  if ((!observation_histogram[variable]) || (bin_width != observation_histogram[variable][0]->bin_width)) {
    int i , j;
    int *pstate , *pioutput;
    double imin_value , *proutput;


    // construction of the histograms

    if (bin_width == D_DEFAULT) {
      bin_width = marginal_histogram[variable]->bin_width;
    }
    imin_value = floor(min_value[variable] / bin_width) * bin_width;

    if (observation_histogram[variable]) {
      for (i = 0;i < nb_state;i++) {
        observation_histogram[variable][i]->nb_bin = (int)floor((max_value[variable] - imin_value) / bin_width) + 1;

        delete [] observation_histogram[variable][i]->frequency;
        observation_histogram[variable][i]->frequency = new int[observation_histogram[variable][i]->nb_bin];
      }
    }

    else {
      observation_histogram[variable] = new Histogram*[nb_state];

      for (i = 0;i < nb_state;i++) {
        observation_histogram[variable][i] = new Histogram((int)floor((max_value[variable] - imin_value) / bin_width) + 1 , false);

        observation_histogram[variable][i]->nb_element = marginal_distribution[0]->frequency[i];
        observation_histogram[variable][i]->type = type[variable];
      }

      // computation of the minimum and maximum values for each state

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
      observation_histogram[variable][i]->bin_width = bin_width;
      observation_histogram[variable][i]->min_value = imin_value;
      observation_histogram[variable][i]->max_value = ceil(max_value[variable] / bin_width) * bin_width;
    }

    // update of the histogram bin frequencies

    for (i = 0;i < nb_state;i++) {
      for (j = 0;j < observation_histogram[variable][i]->nb_bin;j++) {
        observation_histogram[variable][i]->frequency[j] = 0;
      }
    }

    switch (type[variable]) {

    case INT_VALUE : {
      for (i = 0;i < nb_sequence;i++) {
        pstate = int_sequence[i][0];
        pioutput = int_sequence[i][variable];
        for (j = 0;j < length[i];j++) {
//          (observation_histogram[variable][*pstate++]->frequency[(int)((*pioutput++ - imin_value) / bin_width)])++;
          (observation_histogram[variable][*pstate++]->frequency[(int)floor((*pioutput++ - imin_value) / bin_width)])++;
        }
      }
      break;
    }

    case REAL_VALUE : {
      for (i = 0;i < nb_sequence;i++) {
        pstate = int_sequence[i][0];
        proutput = real_sequence[i][variable];
        for (j = 0;j < length[i];j++) {
//          (observation_histogram[variable][*pstate++]->frequency[(int)((*proutput++ - imin_value) / bin_width)]++;
          (observation_histogram[variable][*pstate++]->frequency[(int)floor((*proutput++ - imin_value) / bin_width)])++;
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


/*--------------------------------------------------------------*/
/**
 *  \brief Construction of the observation histograms.
 *
 *  \param[in] nb_state number of states.
 */
/*--------------------------------------------------------------*/

void MarkovianSequences::build_observation_histogram(int nb_state)

{
  if ((nb_variable > 1) && (!observation_histogram)) {
    int i;


    observation_histogram = new Histogram**[nb_variable];
    observation_histogram[0] = NULL;

    for (i = 1;i < nb_variable;i++) {
      observation_histogram[i] = NULL;
      if (marginal_histogram[i]) {
        build_observation_histogram(i , nb_state , marginal_histogram[i]->bin_width);
      }
    }
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Change of the bin width of the marginal histogram and of
 *         the observation histograms for a variable.
 *
 *  \param[in] error      reference on a StatError object,
 *  \param[in] variable   variable index,
 *  \param[in] bin_width  bin_width bin width,
 *  \param[in] imin_value minimum value.
 *
 *  \return               error status.
 */
/*--------------------------------------------------------------*/

bool MarkovianSequences::select_bin_width(StatError &error , int variable ,
                                          double bin_width , double imin_value)

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
    if ((bin_width <= 0.) || ((type[variable] != REAL_VALUE) &&
         (type[variable] != AUXILIARY) && ((int)bin_width != bin_width))) {
      status = false;
      error.update(STAT_error[STATR_HISTOGRAM_BIN_WIDTH]);
    }
    if ((imin_value != D_INF) && ((imin_value <= min_value[variable] - bin_width) ||
         (imin_value > min_value[variable]) || ((type[variable] != REAL_VALUE) &&
          (type[variable] != AUXILIARY) && ((int)imin_value != imin_value)))) {
      status = false;
      error.update(STAT_error[STATR_HISTOGRAM_MIN_VALUE]);
    }
  }

  if (status) {
    build_marginal_histogram(variable , bin_width , imin_value);

    if ((observation_histogram) && (observation_histogram[variable])) {
      build_observation_histogram(variable , marginal_distribution[0]->nb_value , bin_width);
    }
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Test of the overlap of values observed in the different states.
 *
 *  \param[in] variable variable index.
 *
 *  \return             state variable hidden or not.
 */
/*--------------------------------------------------------------*/

bool MarkovianSequences::test_hidden(int variable) const

{
  bool hidden = true;

  if ((variable > 0) && (type[variable] == INT_VALUE)) {
    bool **occurrence;
    int i , j;
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


/*--------------------------------------------------------------*/
/**
 *  \brief Extraction of the probabilities of each value as a function of
 *         the index for a variable.
 *
 *  \param[in] variable variable index.
 */
/*--------------------------------------------------------------*/

void MarkovianSequences::build_index_value(int variable)

{
  int i , j;
  int total , *frequency;


  // construction of a Curves object

  characteristics[variable]->index_value = new Curves(marginal_distribution[variable]->nb_value ,
                                                      max_length , true , false , false);
  frequency = new int[marginal_distribution[variable]->nb_value];

  // computation of the probabilities of each value as a function of the index parameter

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


/*--------------------------------------------------------------*/
/**
 *  \brief Extraction of the probabilities of each value as a function of
 *         the explicit index for a variable.
 *
 *  \param[in] variable variable index.
 */
/*--------------------------------------------------------------*/

void MarkovianSequences::build_explicit_index_value(int variable)

{
  int i , j , k , m;
  int total , *frequency , *index;


  // construction of a Curves object

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

  // computation of the probabilities of each value as a function of the explicit index parameter

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


/*--------------------------------------------------------------*/
/**
 *  \brief Construction of the frequency distributions of the times to the 1st occurrence of
 *         each value of an integer variable.
 *
 *  \param[in] variable variable index.
 */
/*--------------------------------------------------------------*/

void MarkovianSequences::build_first_occurrence_frequency_distribution(int variable)

{
  bool *occurrence;
  int i , j;
  int nb_value , *pisequence;
  FrequencyDistribution **first_occurrence;


  // construction of the frequency distributions

  first_occurrence = new FrequencyDistribution*[marginal_distribution[variable]->nb_value];
  for (i = 0;i < marginal_distribution[variable]->nb_value;i++) {
    first_occurrence[i] = new FrequencyDistribution(max_length);
  }

/*  characteristics[variable]->first_occurrence = new FrequencyDistribution*[marginal_distribution[variable]->nb_value];
  for (i = 0;i < marginal_distribution[variable]->nb_value;i++) {
    characteristics[variable]->first_occurrence[i] = new FrequencyDistribution(max_length);
  } */

  // update of the frequency distributions

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

  // computation of the frequency distribution characteristics

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


/*--------------------------------------------------------------*/
/**
 *  \brief Construction of the frequency distributions of recurrence times in
 *         each value of an integer variable.
 *
 *  \param[in] variable variable index.
 */
/*--------------------------------------------------------------*/

void MarkovianSequences::build_recurrence_time_frequency_distribution(int variable)

{
  int i , j;
  int *index , *pisequence;
  FrequencyDistribution **recurrence_time;


  // construction of the frequency distributions

  recurrence_time = new FrequencyDistribution*[marginal_distribution[variable]->nb_value];
  for (i = 0;i < marginal_distribution[variable]->nb_value;i++) {
    recurrence_time[i] = new FrequencyDistribution(max_length);
  }

/*  characteristics[variable]->recurrence_time = new FrequencyDistribution*[marginal_distribution[variable]->nb_value];
  for (i = 0;i < marginal_distribution[variable]->nb_value;i++) {
    characteristics[variable]->recurrence_time[i] = new FrequencyDistribution(max_length);
  } */

  // update of the frequency distributions

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

  // computation of the frequency distribution characteristics

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


/*--------------------------------------------------------------*/
/**
 *  \brief Construction of the frequency distributions of sojourn times in
 *         each value of an integer variable.
 *
 *  \param[in] variable         variable index,
 *  \param[in] initial_run_flag flag on the construction of left-censored sojourn time frequency distributions.
 */
/*--------------------------------------------------------------*/

void MarkovianSequences::build_sojourn_time_frequency_distribution(int variable , int initial_run_flag)

/* {
  characteristics[variable]->create_sojourn_time_frequency_distribution(max_length , initial_run_flag);
  sojourn_time_frequency_distribution_computation(variable);
} */

{
  int i , j;
  int run_length , *pisequence;
  FrequencyDistribution **sojourn_time , **initial_run , **final_run;


  // construction of the frequency distributions

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

  // update of the frequency distributions

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

  // computation of the frequency distribution characteristics

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


/*--------------------------------------------------------------*/
/**
 *  \brief Update of the frequency distributions of sojourn times in
 *         each value of an integer variable.
 *
 *  \param[in] variable variable index.
 */
/*--------------------------------------------------------------*/

void MarkovianSequences::sojourn_time_frequency_distribution_computation(int variable)

{
  int i , j;
  int run_length , *pisequence;


  // initialization of the frequency distributions

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

  // update of the frequency distributions

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

  // computation of the frequency distribution characteristics

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


/*--------------------------------------------------------------*/
/**
 *  \brief Update of the frequency distributions of censored sojourn times in
 *         each value of an integer variable.
 *
 *  \param[in] initial_run pointer on the left-censored sojourn time frequency distributions,
 *  \param[in] final_run   pointer on the right-censored sojourn time frequency distributions,
 *  \param[in] single_run  pointer on the sequence length frequency distributions for the case of
 *                         a single visited state.
 */
/*--------------------------------------------------------------*/

void MarkovianSequences::censored_sojourn_time_frequency_distribution_computation(FrequencyDistribution **initial_run ,
                                                                                  FrequencyDistribution **final_run ,
                                                                                  FrequencyDistribution **single_run) const

{
  int i , j;
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

  // computation of the characteristics of the left- and right-censored sojourn time frequency distributions
  // and single run frequency distributions

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


/*--------------------------------------------------------------*/
/**
 *  \brief Construction of the frequency distributions of the number of runs
 *         per sequence of each value of an integer variable.
 *
 *  \param[in] variable variable index.
 */
/*--------------------------------------------------------------*/

void MarkovianSequences::build_nb_run_frequency_distribution(int variable)

{
  int i , j;
  int *pisequence , *count;
  FrequencyDistribution **nb_run;


  // construction of the frequency distributions

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

  // update of the frequency distributions

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

  // computation of the frequency distribution characteristics

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


/*--------------------------------------------------------------*/
/**
 *  \brief Construction of the frequency distributions of the number of occurrences
 *         per sequence of each value of an integer variable.
 *
 *  \param[in] variable variable index.
 */
/*--------------------------------------------------------------*/

void MarkovianSequences::build_nb_occurrence_frequency_distribution(int variable)

{
  int i , j;
  int *pisequence , *count;
  FrequencyDistribution **nb_occurrence;


  // construction of the frequency distributions

  nb_occurrence = new FrequencyDistribution*[marginal_distribution[variable]->nb_value];
  for (i = 0;i < marginal_distribution[variable]->nb_value;i++) {
    nb_occurrence[i] = new FrequencyDistribution(max_length + 1);
  }

/*  characteristics[variable]->nb_occurrence = new FrequencyDistribution*[marginal_distribution[variable]->nb_value];
  for (i = 0;i < marginal_distribution[variable]->nb_value;i++) {
    characteristics[variable]->nb_occurrence[i] = new FrequencyDistribution(max_length + 1);
  } */

  // update of the frequency distributions

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

  // computation of the frequency distribution characteristics

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


/*--------------------------------------------------------------*/
/**
 *  \brief Extraction of the characteristics of sequences for categorical variables.
 *         Computation of the number of values, construction of the marginal frequency distribution,
 *         of the probabilities of each value as a function of the index parameter,
 *         construction of the frequency distributions of the times to the 1st occurrence of each value,
 *         of the recurrence times in each value, of the sojourn times in each value,
 *         of the number of runs of each value per sequence,
 *         of the number of occurrences of each value per sequence.
 *
 *  \param[in] variable          variable index,
 *  \param[in] sojourn_time_flag flag on the construction of sojourn time frequency distributions
 *  \param[in] initial_run_flag  flag on the construction of left-censored sojourn time frequency distributions.
 */
/*--------------------------------------------------------------*/

void MarkovianSequences::build_characteristic(int variable , bool sojourn_time_flag ,
                                              bool initial_run_flag)

{
  int i , j;
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

        if (sojourn_time_flag) {
          build_sojourn_time_frequency_distribution(i , initial_run_flag);
        }
        else {
          sojourn_time_frequency_distribution_computation(i);
        }

        if (max_length <= COUNTING_FREQUENCY_MAX_LENGTH) {
          build_nb_run_frequency_distribution(i);
          build_nb_occurrence_frequency_distribution(i);
        }
      }
    }
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Count of words of fixed length.
 *
 *  \param[in] error         reference on a StatError object,
 *  \param[in] display       flag for displaying word counts,
 *  \param[in] variable      variable index,
 *  \param[in] word_length   word length,
 *  \param[in] begin_state   begin state,
 *  \param[in] end_state     end state,
 *  \param[in] min_frequency minimum frequency.
 *
 *  \return                  error status.
 */
/*--------------------------------------------------------------*/

bool MarkovianSequences::word_count(StatError &error , bool display , int variable ,
                                    int word_length , int begin_state , int end_state ,
                                    int min_frequency) const

{
  bool status = true , *selected_word;
  int i , j , k;
  int nb_state , nb_word , max_nb_word , value , max_frequency , total_frequency , width ,
      *power , *frequency , *word_value , *pisequence , *index , **word;
  double nb_word_bound , *probability;
  ios_base::fmtflags format_flags;


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

          // computation of the word score

          pisequence = int_sequence[i][variable] + j;
          value = 0;
          for (k = 0;k < word_length;k++) {
            value += *pisequence++ * power[k];
          }

          // word search

          for (k = 0;k < nb_word;k++) {
            if (value == word_value[k]) {
              frequency[k]++;
              break;
            }
          }

          // word construction

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

    // sort of words by decreasing frequencies

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

    // display of word counts

    if (display) {
      format_flags = cout.setf(ios::right , ios::adjustfield);

      for (j = 0;j < i;j++) {
        for (k = 0;k < word_length;k++) {
          cout << word[index[j]][k] << " ";
        }
        cout << "   " << setw(width) << frequency[index[j]]
             << "   " << probability[index[j]];

        if (j == 0) {
          cout << "   (" << nb_word << " " << SEQ_label[SEQL_WORDS] << ", "
               << total_frequency << ")";
        }
        cout << endl;
      }

      cout.setf(format_flags , ios::adjustfield);
    }

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


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of a the autocorrelation function for a given state and
 *         a variable for the assessment of semi-Markov switching autoregressive models.
 *
 *  \param[in] correl    reference on a Correlation object,
 *  \param[in] state     state,
 *  \param[in] variable  variable index.
 */
/*--------------------------------------------------------------*/

void MarkovianSequences::autocorrelation_computation(Correlation &correl , int state , int variable) const

{
  if (marginal_distribution[0]->frequency[state] >= AUTOCORRELATION_MIN_FREQUENCY) {
    int i , j , k;
    int max_lag;
    double mean , variance , diff;


    // computation of mean and variance

    mean = 0.;

    switch (type[variable]) {

    case INT_VALUE : {
      for (i = 0;i < nb_sequence;i++) {
        for (j = 0;j < length[i];j++) {
          if (int_sequence[i][0][j] == state) {
            mean += int_sequence[i][variable][j];
          }
        }
      }
      break;
    }

    case REAL_VALUE : {
      for (i = 0;i < nb_sequence;i++) {
        for (j = 0;j < length[i];j++) {
          if (int_sequence[i][0][j] == state) {
            mean += real_sequence[i][variable][j];
          }
        }
      }
      break;
    }
    }

    mean /= marginal_distribution[0]->frequency[state];

    variance = 0.;

    switch (type[variable]) {

    case INT_VALUE : {
      for (i = 0;i < nb_sequence;i++) {
        for (j = 0;j < length[i];j++) {
          if (int_sequence[i][0][j] == state) {
            diff = int_sequence[i][variable][j] - mean;
            variance += diff * diff;
          }
        }
      }
      break;
    }

    case REAL_VALUE : {
      for (i = 0;i < nb_sequence;i++) {
        for (j = 0;j < length[i];j++) {
          if (int_sequence[i][0][j] == state) {
            diff = real_sequence[i][variable][j] - mean;
            variance += diff * diff;
          }
        }
      }
      break;
    }
    }

    variance /= marginal_distribution[0]->frequency[state];

    // computation of the autocorrelation coefficients

    correl.point[0][0] = 1.;
    correl.frequency[0] = marginal_distribution[0]->frequency[state];

    for (i = 1;i < correl.length;i++) {
      correl.point[0][i] = 0.;
      correl.frequency[i] = 0;

      switch (type[variable]) {

      case INT_VALUE : {
        for (j = 0;j < nb_sequence;j++) {
          for (k = 0;k < i;k++) {
            if (int_sequence[j][0][k] == state) {
              if ((k == 0) || (int_sequence[j][0][k - 1] != state)) {
                max_lag = 0;
              }
              else {
                max_lag++;
              }
            }
          }

          for (k = i;k < length[j];k++) {
            if (int_sequence[j][0][k] == state) {
              if (int_sequence[j][0][k - 1] != state) {
                max_lag = 0;
              }
              else {
                max_lag++;
              }

              if (i <= max_lag) {
                correl.point[0][i] += (int_sequence[j][variable][k] - mean) * (int_sequence[j][variable][k - i] - mean);
                correl.frequency[i]++;
              }
            }
          }
        }
        break;
      }

      case REAL_VALUE : {
        for (j = 0;j < nb_sequence;j++) {
          for (k = 0;k < i;k++) {
            if (int_sequence[j][0][k] == state) {
              if ((k == 0) || (int_sequence[j][0][k - 1] != state)) {
                max_lag = 0;
              }
              else {
                max_lag++;
              }
            }
          }

          for (k = i;k < length[j];k++) {
            if (int_sequence[j][0][k] == state) {
              if (int_sequence[j][0][k - 1] != state) {
                max_lag = 0;
              }
              else {
                max_lag++;
              }

              if (i <= max_lag) {
                correl.point[0][i] += (real_sequence[j][variable][k] - mean) * (real_sequence[j][variable][k - i] - mean);
                correl.frequency[i]++;
              }
            }
          }
        }
        break;
      }
      }

      if ((correl.frequency[i] >= marginal_distribution[0]->frequency[state] * AUTOCORRELATION_FREQUENCY_RATIO) &&
          (correl.frequency[i] >= AUTOCORRELATION_MIN_FREQUENCY)) {
        correl.point[0][i] /= correl.frequency[i] * variance;
      }
      else {
        correl.length = i;
        break;
      }
    }
  }

  else {
    correl.length = 0;
  }
}


};  // namespace sequence_analysis
