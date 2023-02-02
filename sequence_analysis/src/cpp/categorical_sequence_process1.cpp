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



#include "stat_tool/stat_label.h"

#include <boost/tokenizer.hpp>

#include "sequences.h"
#include "sequence_label.h"

using namespace std;
using namespace boost;
using namespace stat_tool;


namespace sequence_analysis {



/*--------------------------------------------------------------*/
/**
 *  \brief Construction of the distributions of a CategoricalSequenceProcess object
 *         except the observation distributions.
 *
 *  \param[in] ilength       sequence length distribution,
 *  \param[in] homogeneity   state homogeneity,
 *  \param[in] counting_flag flag counting distributions.
 */
/*--------------------------------------------------------------*/

void CategoricalSequenceProcess::create_characteristic(const Distribution &ilength ,
                                                       bool *homogeneity ,
                                                       bool counting_flag)

{
  bool homogeneous = true;
  int i;
  int max_length = ilength.nb_value - 1;


  if (length) {
    delete length;
  }
  length = new Distribution(ilength);

  if (index_value) {
    delete index_value;
  }
  index_value = new Curves(nb_value , max_length , false , false , false);

  if (!no_occurrence) {
    no_occurrence = new double[nb_value];
  }

  for (i = 0;i < nb_value;i++) {
    no_occurrence[i] = 0.;
  }

  if (first_occurrence) {
    for (i = 0;i < nb_value;i++) {
      delete first_occurrence[i];
    }
  }
  else {
    first_occurrence = new Distribution*[nb_value];
  }

  for (i = 0;i < nb_value;i++) {
    first_occurrence[i] = new Distribution(NB_VALUE);
  }

  if (!absorption) {
    absorption = new double[nb_value];
  }

  for (i = 0;i < nb_value;i++) {
    absorption[i] = 0.;
  }

  if (sojourn_time) {
    for (i = 0;i < nb_value;i++) {
      delete sojourn_time[i];
    }
  }
  else {
    sojourn_time = new DiscreteParametric*[nb_value];
  }

  for (i = 0;i < nb_value;i++) {
    if (homogeneity[i]) {
      sojourn_time[i] = new DiscreteParametric(NB_VALUE);
    }
    else {
      sojourn_time[i] = NULL;
      homogeneous = false;
    }
  }

  if (homogeneous) {
    if (!leave) {
      leave = new double[nb_value];
    }

    for (i = 0;i < nb_value;i++) {
      leave[i] = 0.;
    }

    if (recurrence_time) {
      for (i = 0;i < nb_value;i++) {
        delete recurrence_time[i];
      }
    }
    else {
      recurrence_time = new Distribution*[nb_value];
    }

    for (i = 0;i < nb_value;i++) {
      recurrence_time[i] = new Distribution(NB_VALUE);
    }
  }

  if (counting_flag) {
    if (nb_run) {
      for (i = 0;i < nb_value;i++) {
        delete nb_run[i];
      }
    }
    else {
      nb_run = new Distribution*[nb_value];
    }

    for (i = 0;i < nb_value;i++) {
      nb_run[i] = new Distribution((max_length % 2 == 0 ?
                                    max_length / 2 : max_length / 2 + 1) + 1);
    }

    if (nb_occurrence) {
      for (i = 0;i < nb_value;i++) {
        delete nb_occurrence[i];
      }
    }
    else {
      nb_occurrence = new Distribution*[nb_value];
    }

    for (i = 0;i < nb_value;i++) {
      nb_occurrence[i] = new Distribution(max_length + 1);
    }
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Construction of the distributions of a CategoricalSequenceProcess object
 *         except the observation distributions.
 *
 *  \param[in] ilength           sequence length distribution,
 *  \param[in] sojourn_time_flag flag sojourn time distributions,
 *  \param[in] counting_flag     flag counting distributions.
 */
/*--------------------------------------------------------------*/

void CategoricalSequenceProcess::create_characteristic(const Distribution &ilength ,
                                                       bool sojourn_time_flag , bool counting_flag)

{
  int i;
  int max_length = ilength.nb_value - 1;


  if (length) {
    delete length;
  }
  length = new Distribution(ilength);

  if (index_value) {
    delete index_value;
  }
  index_value = new Curves(nb_value , max_length , false , false , false);

  if (!no_occurrence) {
    no_occurrence = new double[nb_value];
  }

  for (i = 0;i < nb_value;i++) {
    no_occurrence[i] = 0.;
  }

  if (first_occurrence) {
    for (i = 0;i < nb_value;i++) {
      delete first_occurrence[i];
    }
  }
  else {
    first_occurrence = new Distribution*[nb_value];
  }

  for (i = 0;i < nb_value;i++) {
    first_occurrence[i] = new Distribution(NB_VALUE);
  }

  if (!leave) {
    leave = new double[nb_value];
  }

  for (i = 0;i < nb_value;i++) {
    leave[i] = 0.;
  }

  if (recurrence_time) {
    for (i = 0;i < nb_value;i++) {
      delete recurrence_time[i];
    }
  }
  else {
    recurrence_time = new Distribution*[nb_value];
  }

  for (i = 0;i < nb_value;i++) {
    recurrence_time[i] = new Distribution(NB_VALUE);
  }

  if (sojourn_time_flag) {
    if (!absorption) {
      absorption = new double[nb_value];
    }

    for (i = 0;i < nb_value;i++) {
      absorption[i] = 0.;
    }

    if (sojourn_time) {
      for (i = 0;i < nb_value;i++) {
        delete sojourn_time[i];
      }
    }
    else {
      sojourn_time = new DiscreteParametric*[nb_value];
    }

    for (i = 0;i < nb_value;i++) {
      sojourn_time[i] = new DiscreteParametric(NB_VALUE);
    }
  }

  if (counting_flag) {
    if (nb_run) {
      for (i = 0;i < nb_value;i++) {
        delete nb_run[i];
      }
    }
    else {
      nb_run = new Distribution*[nb_value];
    }

    for (i = 0;i < nb_value;i++) {
      nb_run[i] = new Distribution((max_length % 2 == 0 ?
                                    max_length / 2 : max_length / 2 + 1) + 1);
    }

    if (nb_occurrence) {
      for (i = 0;i < nb_value;i++) {
        delete nb_occurrence[i];
      }
    }
    else {
      nb_occurrence = new Distribution*[nb_value];
    }

    for (i = 0;i < nb_value;i++) {
      nb_occurrence[i] = new Distribution(max_length + 1);
    }
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Constructor of the CategoricalSequenceProcess class.
 *
 *  \param[in] inb_state        number of states,
 *  \param[in] inb_value        number of categories,
 *  \param[in] observation_flag flag observation distributions.
 */
/*--------------------------------------------------------------*/

CategoricalSequenceProcess::CategoricalSequenceProcess(int inb_state , int inb_value ,
                                                       int observation_flag)
:CategoricalProcess(inb_state , inb_value , observation_flag)

{
  length = NULL;
  index_value = NULL;
  no_occurrence = NULL;
  first_occurrence = NULL;
  leave = NULL;
  recurrence_time = NULL;
  absorption = NULL;
  sojourn_time = NULL;
  nb_run = NULL;
  nb_occurrence = NULL;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Constructor of the CategoricalSequenceProcess class.
 *
 *  \param[in] inb_state number of states,
 *  \param[in] occupancy state occupancy distributions.
 */
/*--------------------------------------------------------------*/

CategoricalSequenceProcess::CategoricalSequenceProcess(int inb_state ,
                                                       DiscreteParametric **occupancy)

{
  int i;


  nb_state = inb_state;
  nb_value = inb_state;

  observation = NULL;

  length = NULL;
  index_value = NULL;
  no_occurrence = NULL;
  first_occurrence = NULL;
  leave = NULL;
  recurrence_time = NULL;

  absorption = new double[nb_state];
  sojourn_time = new DiscreteParametric*[nb_state];
  for (i = 0;i < nb_state;i++) {
    if (occupancy[i]) {
      sojourn_time[i] = new DiscreteParametric(*occupancy[i]);
    }
    else {
      sojourn_time[i] = NULL;
    }
  }

  nb_run = NULL;
  nb_occurrence = NULL;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Construction of a CategoricalSequenceProcess object from
 *         a CategoricalProcess object.
 *
 *  \param[in] process reference on a CategoricalProcess object.
 */
/*--------------------------------------------------------------*/

CategoricalSequenceProcess::CategoricalSequenceProcess(const CategoricalProcess &process)
:CategoricalProcess(process)

{
  length = NULL;
  index_value = NULL;
  no_occurrence = NULL;
  first_occurrence = NULL;
  leave = NULL;
  recurrence_time = NULL;
  absorption = NULL;
  sojourn_time = NULL;
  nb_run = NULL;
  nb_occurrence = NULL;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Copy of a CategoricalSequenceProcess object.
 *
 *  \param[in] process             reference on a CategoricalSequenceProcess object,
 *  \param[in] characteristic_flag flag copy of the characteristic distributions.
 */
/*--------------------------------------------------------------*/

void CategoricalSequenceProcess::copy(const CategoricalSequenceProcess &process ,
                                      bool characteristic_flag)

{
  if (characteristic_flag) {
    int i;


    if (process.length) {
      length = new Distribution(*(process.length));
    }
    else {
      length = NULL;
    }

    if (process.index_value) {
      index_value = new Curves(*(process.index_value));
    }
    else {
      index_value = NULL;
    }

    if (process.no_occurrence) {
      no_occurrence = new double[nb_value];
      for (i = 0;i < nb_value;i++) {
        no_occurrence[i] = process.no_occurrence[i];
      }
    }
    else {
      no_occurrence = NULL;
    }

    if (process.first_occurrence) {
      first_occurrence = new Distribution*[nb_value];
      for (i = 0;i < nb_value;i++) {
        first_occurrence[i] = new Distribution(*(process.first_occurrence[i]));
      }
    }
    else {
      first_occurrence = NULL;
    }

    if (process.leave) {
      leave = new double[nb_value];
      for (i = 0;i < nb_value;i++) {
        leave[i] = process.leave[i];
      }
    }
    else {
      leave = NULL;
    }

    if (process.recurrence_time) {
      recurrence_time = new Distribution*[nb_value];
      for (i = 0;i < nb_value;i++) {
        if (process.recurrence_time[i]) {
          recurrence_time[i] = new Distribution(*(process.recurrence_time[i]));
        }
        else {
          recurrence_time[i] = NULL;
        }
      }
    }
    else {
      recurrence_time = NULL;
    }

    if (process.absorption) {
      absorption = new double[nb_value];
      for (i = 0;i < nb_value;i++) {
        absorption[i] = process.absorption[i];
      }
    }
    else {
      absorption = NULL;
    }

    if (process.sojourn_time) {
      sojourn_time = new DiscreteParametric*[nb_value];
      for (i = 0;i < nb_value;i++) {
        if (process.sojourn_time[i]) {
          sojourn_time[i] = new DiscreteParametric(*(process.sojourn_time[i]));
        }
        else {
          sojourn_time[i] = NULL;
        }
      }
    }
    else {
      sojourn_time = NULL;
    }

    if (process.nb_run) {
      nb_run = new Distribution*[nb_value];
      for (i = 0;i < nb_value;i++) {
        nb_run[i] = new Distribution(*(process.nb_run[i]));
      }
    }
    else {
      nb_run = NULL;
    }

    if (process.nb_occurrence) {
      nb_occurrence = new Distribution*[nb_value];
      for (i = 0;i < nb_value;i++) {
        nb_occurrence[i] = new Distribution(*(process.nb_occurrence[i]));
      }
    }
    else {
      nb_occurrence = NULL;
    }
  }

  else {
    length = NULL;
    index_value = NULL;
    no_occurrence = NULL;
    first_occurrence = NULL;
    leave = NULL;
    recurrence_time = NULL;
    absorption = NULL;
    sojourn_time = NULL;
    nb_run = NULL;
    nb_occurrence = NULL;
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Copy of the state occupancy distributions.
 *
 *  \param[in] process            reference on a CategoricalSequenceProcess object,
 *  \param[in] occupancy_nb_value number of allocated values for the state occupancy distributions.
 */
/*--------------------------------------------------------------*/

void CategoricalSequenceProcess::init_occupancy(const CategoricalSequenceProcess &process ,
                                                int occupancy_nb_value)

{
  int i;


  nb_state = process.nb_state;
  nb_value = process.nb_value;

  observation = NULL;

  length = NULL;
  index_value = NULL;
  no_occurrence = NULL;
  first_occurrence = NULL;
  leave = NULL;
  recurrence_time = NULL;

  absorption = new double[nb_value];
  sojourn_time = new DiscreteParametric*[nb_value];
  for (i = 0;i < nb_value;i++) {
    absorption[i] = process.absorption[i];
    if ((process.sojourn_time[i]) && (process.sojourn_time[i]->ident != CATEGORICAL)) {
      sojourn_time[i] = new DiscreteParametric(*(process.sojourn_time[i]) ,
                                               DISTRIBUTION_COPY , occupancy_nb_value);
    }
    else {
      sojourn_time[i] = NULL;
    }
  }

  nb_run = NULL;
  nb_occurrence = NULL;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Constructor by copy of the CategoricalSequenceProcess class.
 *
 *  \param[in] process   reference on a CategoricalSequenceProcess object,
 *  \param[in] transform type of transform (CATEGORICAL_SEQUENCE_PROCESS_COPY/INIT_OCCUPANCY),
 *  \param[in] param     flag on the computation of the characteristic distributions/
 *                       number of allocated values for the state occupancy distributions.
 */
/*--------------------------------------------------------------*/

CategoricalSequenceProcess::CategoricalSequenceProcess(const CategoricalSequenceProcess &process ,
                                                       categorical_sequence_process_transformation transform ,
                                                       int param)

{
  switch (transform) {
  case CATEGORICAL_SEQUENCE_PROCESS_COPY :
    CategoricalProcess::copy(process);
    copy(process , param);
    break;
  case INIT_OCCUPANCY :
    init_occupancy(process , param);
    break;
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Destruction of the data members of a CategoricalSequenceProcess object.
 */
/*--------------------------------------------------------------*/

void CategoricalSequenceProcess::remove()

{
  int i;


  if (length) {
    delete length;

    length = NULL;
  }

  if (index_value) {
    delete index_value;

    index_value = NULL;
  }

  if (no_occurrence) {
    delete [] no_occurrence;

    no_occurrence = NULL;
  }

  if (first_occurrence) {
    for (i = 0;i < nb_value;i++) {
      delete first_occurrence[i];
    }
    delete [] first_occurrence;

    first_occurrence = NULL;
  }

  if (leave) {
    delete [] leave;

    leave = NULL;
  }

  if (recurrence_time) {
    for (i = 0;i < nb_value;i++) {
      delete recurrence_time[i];
    }
    delete [] recurrence_time;

    recurrence_time = NULL;
  }

  if (absorption) {
    delete [] absorption;

    absorption = NULL;
  }

  if (sojourn_time) {
    for (i = 0;i < nb_value;i++) {
      delete sojourn_time[i];
    }
    delete [] sojourn_time;

    sojourn_time = NULL;
  }

  if (nb_run) {
    for (i = 0;i < nb_value;i++) {
      delete nb_run[i];
    }
    delete [] nb_run;

    nb_run = NULL;
  }

  if (nb_occurrence) {
    for (i = 0;i < nb_value;i++) {
      delete nb_occurrence[i];
    }
    delete [] nb_occurrence;

    nb_occurrence = NULL;
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Destructor of the CategoricalSequenceProcess class.
 */
/*--------------------------------------------------------------*/

CategoricalSequenceProcess::~CategoricalSequenceProcess()

{
  remove();
}


/*--------------------------------------------------------------*/
/**
 *  \brief Assignment operator of the CategoricalSequenceProcess class.
 *
 *  \param[in] process reference on a CategoricalSequenceProcess object.
 *
 *  \return            CategoricalSequenceProcess object.
 */
/*--------------------------------------------------------------*/

CategoricalSequenceProcess& CategoricalSequenceProcess::operator=(const CategoricalSequenceProcess &process)

{
  if (&process != this) {
    remove();
    CategoricalProcess::remove();

    CategoricalProcess::copy(process);
    copy(process);
  }

  return *this;
}


/*--------------------------------------------------------------*/
/**
 *  \brief test model hidden.
 *
 *  \param[in] nb_output_process number of observation processes,
 *  \param[in] process           pointer on the observation processes.
 *
 *  \return                      model hidden or not.
 */
/*--------------------------------------------------------------*/

bool CategoricalSequenceProcess::test_hidden(int nb_output_process , CategoricalSequenceProcess **process)

{
  bool hidden = false;
  int i;


  for (i = 0;i < nb_output_process;i++) {
    if (process[i]) {
      hidden = process[i]->CategoricalProcess::test_hidden();
      if (hidden) {
        break;
      }
    }

    else {
      hidden = true;
      break;
    }
  }

  return hidden;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Analysis of the format of state occupancy distributions.
 *
 *  \param[in] error           reference on a StatError object,
 *  \param[in] in_file         stream,
 *  \param[in] line            reference on the file line index,
 *  \param[in] chain           reference on a Chain object,
 *  \param[in] cumul_threshold threshold on the cumulative distribution function.
 *
 *  \return                    CategoricalSequenceProcess object.
 */
/*--------------------------------------------------------------*/

CategoricalSequenceProcess* CategoricalSequenceProcess::occupancy_parsing(StatError &error , ifstream &in_file ,
                                                                          int &line , const Chain &chain ,
                                                                          double cumul_threshold)

{
  string buffer;
  size_t position;
  typedef tokenizer<char_separator<char>> tokenizer;
  char_separator<char> separator(" \t");
  bool status = true , lstatus;
  int i , j;
  int index;
  DiscreteParametric **dist;
  CategoricalSequenceProcess *process;


  process = NULL;

  dist = new DiscreteParametric*[chain.nb_state];
  for (i = 0;i < chain.nb_state;i++) {
    dist[i] = NULL;
  }

  for (i = 0;i < chain.nb_state;i++) {
    if (chain.transition[i][i] == 0.) {
      while (getline(in_file , buffer)) {
        line++;

#       ifdef DEBUG
        cout << line << "  " << buffer << endl;
#       endif

        position = buffer.find('#');
        if (position != string::npos) {
          buffer.erase(position);
        }
        j = 0;

        tokenizer tok_buffer(buffer , separator);

        for (tokenizer::iterator token = tok_buffer.begin();token != tok_buffer.end();token++) {
          switch (j) {

          // test STATE keyword

          case 0 : {
            if (*token != STAT_word[STATW_STATE]) {
              status = false;
              error.correction_update(STAT_parsing[STATP_KEYWORD] , STAT_word[STATW_STATE] , line , j + 1);
            }
            break;
          }

          // test state index

          case 1 : {
            lstatus = true;

/*            try {
              index = stoi(*token);   in C++ 11
            }
            catch(invalid_argument &arg) {
              lstatus = false;
            } */
            index = atoi(token->c_str());

            if ((lstatus) && (index != i)) {
              lstatus = false;
            }

            if (!lstatus) {
              status = false;
              error.correction_update(STAT_parsing[STATP_STATE_INDEX] , i , line , j + 1);
            }
            break;
          }

          // test OCCUPANCY_DISTRIBUTION keyword

          case 2 : {
            if (*token != SEQ_word[SEQW_OCCUPANCY_DISTRIBUTION]) {
              status = false;
              error.correction_update(STAT_parsing[STATP_KEYWORD] , SEQ_word[SEQW_OCCUPANCY_DISTRIBUTION] , line , j + 1);
            }
            break;
          }
          }

          j++;
        }

        if (j > 0) {
          if (j != 3) {
            status = false;
            error.update(STAT_parsing[STATP_FORMAT] , line);
          }

          dist[i] = DiscreteParametric::parsing(error , in_file , line , UNIFORM ,
                                                cumul_threshold , 1);
          if (!dist[i]) {
            status = false;
          }
          else if (dist[i]->mean == 1.) {
            delete dist[i];
            dist[i] = NULL;
          }

          break;
        }
      }

      if ((j == 0) && (!dist[i])) {
        status = false;
        error.update(STAT_parsing[STATP_FORMAT] , line);
      }
    }
  }

  if (status) {
    process = new CategoricalSequenceProcess(chain.nb_state , dist);
  }

  for (i = 0;i < chain.nb_state;i++) {
    delete dist[i];
  }
  delete [] dist;

  return process;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the marginal state distribution for an ordinary process.
 *
 *  \return Distribution object.
 */
/*--------------------------------------------------------------*/

Distribution* CategoricalSequenceProcess::weight_computation() const

{
  int i , j;
  double sum;
  Distribution *weight;


  weight = new Distribution(nb_state);

  for (i = 0;i < nb_state;i++) {
    weight->mass[i] = 0.;
  }

  for (i = 0;i < length->nb_value - 1;i++) {
    for (j = 0;j < nb_state;j++) {
      weight->mass[j] += index_value->point[j][i] * (1. - length->cumul[i]);
    }
  }

  sum = 0.;
  for (i = 0;i < nb_state;i++) {
    sum += weight->mass[i];
  }
  for (i = 0;i < nb_state;i++) {
    weight->mass[i] /= sum;
  }

  weight->cumul_computation();
  weight->max_computation();

  return weight;
}


};  // namespace sequence_analysis
