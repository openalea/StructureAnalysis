/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       StructureAnalysis: Identifying patterns in plant architecture and development
 *
 *       Copyright 1995-2019 CIRAD AGAP
 *
 *       File author(s): Yann Guedon (yann.guedon@cirad.fr)
 *
 *       $Source$
 *       $Id$
 *
 *       Forum for StructureAnalysis developers:
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



#include <cmath>

#include <sstream>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <boost/tokenizer.hpp>

#include "stat_tool/stat_label.h"

#include "renewal.h"
#include "sequence_label.h"

using namespace std;
using namespace boost;
using namespace stat_tool;


namespace sequence_analysis {



/*--------------------------------------------------------------*/
/**
 *  \brief Construction of the observation period frequency distribution and
 *         the number of events frequency distributions on the basis of
 *         triplets {observation period, number of events, frequency}.
 */
/*--------------------------------------------------------------*/

void TimeEvents::build_frequency_distribution()

{
  int i;
  int max_nb_event , *ptime , *pnb_event , *pfrequency;


  // construction of the observation period frequency distribution and
  // the number of events frequency distributions

  htime = new FrequencyDistribution(time[nb_class - 1] + 1);

  hnb_event = new FrequencyDistribution*[time[nb_class - 1] + 1];
  for (i = 0;i <= time[nb_class - 1];i++) {
    hnb_event[i] = NULL;
  }

  ptime = time;
  pnb_event = nb_event;
  pfrequency = frequency;
  max_nb_event = 0;

  for (i = 0;i < nb_class - 1;i++) {
    if (*(ptime + 1) != *ptime) {
      hnb_event[*ptime] = new FrequencyDistribution(*pnb_event + 1);
      if (*pnb_event > max_nb_event) {
        max_nb_event = *pnb_event;
      }
    }
    htime->frequency[*ptime++] += *pfrequency++;
    pnb_event++;
  }

  hnb_event[*ptime] = new FrequencyDistribution(*pnb_event + 1);
  if (*pnb_event > max_nb_event) {
    max_nb_event = *pnb_event;
  }
  htime->frequency[*ptime] += *pfrequency;

  mixture = new FrequencyDistribution(max_nb_event + 1);

  htime->offset_computation();
  htime->nb_element = nb_element;
  htime->max_computation();
  htime->mean_computation();
  htime->variance_computation();

  // update of the number of events frequency distributions

  ptime = time;
  pnb_event = nb_event;
  pfrequency = frequency;

  for (i = 0;i < nb_class;i++) {
    hnb_event[*ptime++]->frequency[*pnb_event] += *pfrequency;
    mixture->frequency[*pnb_event++] += *pfrequency++;
  }

  for (i = htime->offset;i < htime->nb_value;i++) {
    if (htime->frequency[i]) {
      hnb_event[i]->offset_computation();
      hnb_event[i]->nb_element_computation();
      hnb_event[i]->max_computation();
      hnb_event[i]->mean_computation();
      hnb_event[i]->variance_computation();
    }
  }

  mixture->offset_computation();
  mixture->nb_element = nb_element;
  mixture->max_computation();
  mixture->mean_computation();
  mixture->variance_computation();
}


/*--------------------------------------------------------------*/
/**
 *  \brief Construction of the triplets {observation period, number of events, frequency}
 *         on the basis of the observation period frequency distribution and
 *         the number of events frequency distributions.
 */
/*--------------------------------------------------------------*/

void TimeEvents::build_sample()

{
  int i , j;
  int *ptime , *pnb_event , *pfrequency , *hfrequency;


  nb_class = 0;

  for (i = htime->offset;i < htime->nb_value;i++) {
    if (htime->frequency[i] > 0) {
      hfrequency = hnb_event[i]->frequency + hnb_event[i]->offset;
      for (j = hnb_event[i]->offset;j < hnb_event[i]->nb_value;j++) {
        if (*hfrequency++ > 0) {
          nb_class++;
        }
      }
    }
  }

  time = new int[nb_class];
  nb_event = new int[nb_class];
  frequency = new int[nb_class];

  ptime = time;
  pnb_event = nb_event;
  pfrequency = frequency;

  for (i = htime->offset;i < htime->nb_value;i++) {
    if (htime->frequency[i] > 0) {
      hfrequency = hnb_event[i]->frequency + hnb_event[i]->offset;
      for (j = hnb_event[i]->offset;j < hnb_event[i]->nb_value;j++) {
        if (*hfrequency > 0) {
          *ptime++ = i;
          *pnb_event++ = j;
          *pfrequency++ = *hfrequency;
        }
        hfrequency++;
      }
    }
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Construction of a TimeEvents object from pairs {observation period, number of events}.
 *
 *  \param[in] inb_element number of pairs {observation period, number of events},
 *  \param[in] itime       pointer on the observation periods,
 *  \param[in] inb_event   pointer on the numbers of events.
 */
/*--------------------------------------------------------------*/

void TimeEvents::build(int inb_element , int *itime , int *inb_event)

{
  int i , j , k;
  int btime , min_time , max_time , bnb_event , min_nb_event , max_nb_event ,
      nb_selected , *ptime , *pnb_event , *selected_nb_event , *snb_event;


  nb_element = inb_element;

  // computation of the maximum/minimum observation periods and the maximum/minimum numbers of events

  ptime = itime;
  pnb_event = inb_event;
  max_time = 1;
  max_nb_event = 0;

  for (i = 0;i < nb_element;i++) {
    if (*ptime > max_time) {
      max_time = *ptime;
    }
    if (*pnb_event > max_nb_event) {
      max_nb_event = *pnb_event;
    }
    ptime++;
    pnb_event++;
  }

  ptime = itime;
  pnb_event = inb_event;
  min_time = max_time;
  min_nb_event = max_nb_event;

  for (i = 0;i < nb_element;i++) {
    if (*ptime < min_time) {
      min_time = *ptime;
    }
    if (*pnb_event < min_nb_event) {
      min_nb_event = *pnb_event;
    }
    ptime++;
    pnb_event++;
  }

  nb_class = MIN(nb_element , (max_time - min_time + 1) *
                 (max_nb_event - min_nb_event + 1));
  time = new int[nb_class];
  nb_event = new int[nb_class];
  frequency = new int[nb_class];

  selected_nb_event = new int[nb_element];

  // sort of the pairs {observation period, number of events} first by increasing observation period
  // then by increasing number of events

  btime = 0;
  nb_class = 0;
  i = 0;

  do {

    // search for the minimum unselected observation period

    ptime = itime;
    min_time = max_time;
    for (j = 0;j < nb_element;j++) {
      if ((*ptime > btime) && (*ptime < min_time)) {
        min_time = *ptime;
      }
      ptime++;
    }
    btime = min_time;

    // extraction of the pairs corresponding to the selected observation period

    ptime = itime;
    pnb_event = inb_event;
    nb_selected = 0;
    for (j = 0;j < nb_element;j++) {
      if (*ptime == btime) {
        selected_nb_event[nb_selected++] = *pnb_event;
      }
      ptime++;
      pnb_event++;
    }

    // sort of the pairs corresponding to the selected observation period

    bnb_event = -1;
    j = 0;

    do {

      // search for the minimum unselected number of events

      snb_event = selected_nb_event;
      min_nb_event = max_nb_event;
      for (k = 0;k < nb_selected;k++) {
        if ((*snb_event > bnb_event) && (*snb_event < min_nb_event)) {
          min_nb_event = *snb_event;
        }
        snb_event++;
      }
      bnb_event = min_nb_event;

      // constitution of the triplets {observation period, number of events, frequency}

      time[nb_class] = btime;
      nb_event[nb_class] = bnb_event;

      frequency[nb_class] = 0;
      snb_event = selected_nb_event;
      for (k = 0;k < nb_selected;k++) {
        if (*snb_event == bnb_event) {
          i++;
          j++;
          frequency[nb_class]++;
        }
        snb_event++;
      }
      nb_class++;
    }
    while (j < nb_selected);

  }
  while (i < nb_element);

  delete [] selected_nb_event;

  build_frequency_distribution();
}


/*--------------------------------------------------------------*/
/**
 *  \brief Constructor of the TimeEvents class.
 *
 *  \param[in] inb_class number of classes corresponding to fixed {observation period, number of events}.
 */
/*--------------------------------------------------------------*/

TimeEvents::TimeEvents(int inb_class)

{
  nb_element = 0;
  nb_class = inb_class;

  if (nb_class == 0) {
    time = NULL;
    nb_event = NULL;
    frequency = NULL;
  }

  else {
    int i;

    time = new int[nb_class];
    nb_event = new int[nb_class];
    frequency = new int[nb_class];

    for (i = 0;i < nb_class;i++) {
      time[i] = 0;
      nb_event[i] = 0;
      frequency[i] = 0;
    }
  }

  htime = NULL;
  hnb_event = NULL;
  mixture = NULL;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Construction of a TimeEvents object from a FrequencyDistribution object.
 *
 *  \param[in] inb_event reference on a FrequencyDistribution object,
 *  \param[in] itime     observation period.
 */
/*--------------------------------------------------------------*/

TimeEvents::TimeEvents(FrequencyDistribution &inb_event, int itime)

{
  int i;
  int *ptime , *pnb_event , *pfrequency;


  nb_class = 0;
  for (i = inb_event.offset;i < inb_event.nb_value;i++) {
    if (inb_event.frequency[i] > 0) {
      nb_class++;
    }
  }

  time = new int[nb_class];
  nb_event = new int[nb_class];
  frequency = new int[nb_class];

  nb_element = inb_event.nb_element;

  // constitution of the triplets {observation period, number of events, frequency}

  ptime = time;
  pnb_event = nb_event;
  pfrequency = frequency;

  for (i = inb_event.offset;i < inb_event.nb_value;i++) {
    if (inb_event.frequency[i] > 0) {
      *ptime++ = itime;
      *pnb_event++ = i;
      *pfrequency++ = inb_event.frequency[i];
    }
  }

  // construction of the observation period frequency distribution and
  // the number of events frequency distributions

  htime = new FrequencyDistribution(itime + 1);

  htime->frequency[itime] = inb_event.nb_element;
  htime->offset = itime;
  htime->nb_element = inb_event.nb_element;
  htime->max = inb_event.nb_element;
  htime->mean = itime;
  htime->variance = 0.;

  hnb_event = new FrequencyDistribution*[itime + 1];
  for (i = 0;i < itime;i++) {
    hnb_event[i] = NULL;
  }
  hnb_event[itime] = new FrequencyDistribution(inb_event);

  mixture = new FrequencyDistribution(inb_event);
}


/*--------------------------------------------------------------*/
/**
 *  \brief Copy of a TimeEvents object.
 *
 *  \param[in] timev reference on a TimeEvents object.
 */
/*--------------------------------------------------------------*/

void TimeEvents::copy(const TimeEvents &timev)

{
  int i;


  // copy of the triplets {observation period, number of events, frequency}

  nb_element = timev.nb_element;
  nb_class = timev.nb_class;

  time = new int[nb_class];
  nb_event = new int[nb_class];
  frequency = new int[nb_class];

  for (i = 0;i < nb_class;i++) {
    time[i] = timev.time[i];
    nb_event[i] = timev.nb_event[i];
    frequency[i] = timev.frequency[i];
  }

  // copy of the observation period frequency distribution and
  // the number of events frequency distributions

  htime = new FrequencyDistribution(*(timev.htime));

  hnb_event = new FrequencyDistribution*[htime->nb_value];

  for (i = 0;i < htime->offset;i++) {
    hnb_event[i] = NULL;
  }

  for (i = htime->offset;i < htime->nb_value;i++) {
    if (htime->frequency[i] > 0) {
      hnb_event[i] = new FrequencyDistribution(*(timev.hnb_event[i]));
    }
    else {
      hnb_event[i] = NULL;
    }
  }

  mixture = new FrequencyDistribution(*(timev.mixture));
}


/*--------------------------------------------------------------*/
/**
 *  \brief Merging of TimeEvents objects.
 *
 *  \param[in] nb_sample number of TimeEvents objects,
 *  \param[in] ptimev    pointer on the TimeEvents objects.
 */
/*--------------------------------------------------------------*/

void TimeEvents::merge(int nb_sample , const TimeEvents **ptimev)

{
  int i , j;
  int nb_histo;
  const FrequencyDistribution **phisto;


  nb_element = 0;
  for (i = 0;i < nb_sample;i++) {
    nb_element += ptimev[i]->nb_element;
  }

  phisto = new const FrequencyDistribution*[nb_sample];

  // merging of the observation period frequency distributions

  for (i = 0;i < nb_sample;i++) {
    phisto[i] = ptimev[i]->htime;
  }
  htime = new FrequencyDistribution(nb_sample , phisto);

  // merging of the number of events frequency distributions for a given observation period

  hnb_event = new FrequencyDistribution*[htime->nb_value];

  for (i = 0;i < htime->offset;i++) {
    hnb_event[i] = NULL;
  }

  for (i = htime->offset;i < htime->nb_value;i++) {
    if (htime->frequency[i]) {
      nb_histo = 0;
      for (j = 0;j < nb_sample;j++) {
        if ((i < ptimev[j]->htime->nb_value) && (ptimev[j]->hnb_event[i])) {
          phisto[nb_histo++] = ptimev[j]->hnb_event[i];
        }
      }
      hnb_event[i] = new FrequencyDistribution(nb_histo , phisto);
    }

    else {
      hnb_event[i] = NULL;
    }
  }

  for (i = 0;i < nb_sample;i++) {
    phisto[i] = ptimev[i]->mixture;
  }
  mixture = new FrequencyDistribution(nb_sample , phisto);

  delete [] phisto;

  build_sample();
}


/*--------------------------------------------------------------*/
/**
 *  \brief Merging of TimeEvents objects.
 *
 *  \param[in] nb_sample number of TimeEvents objects,
 *  \param[in] itimev    pointer on the TimeEvents objects.
 *
 *  \return              TimeEvents object.
 */
/*--------------------------------------------------------------*/

TimeEvents* TimeEvents::merge(int nb_sample , const vector<TimeEvents> &itimev) const

{
  int i;
  TimeEvents *timev;
  const TimeEvents **ptimev;


  nb_sample++;
  ptimev = new const TimeEvents*[nb_sample];

  ptimev[0] = this;
  for (i = 1;i < nb_sample;i++) {
    ptimev[i] = new TimeEvents(itimev[i - 1]);
  }

  timev = new TimeEvents(nb_sample , ptimev);

  for (i = 1;i < nb_sample;i++) {
    delete ptimev[i];
  }
  delete [] ptimev;

  return timev;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Destruction of the data members of a TimeEvents object.
 */
/*--------------------------------------------------------------*/

void TimeEvents::remove()

{
  int i;


  delete [] time;
  delete [] nb_event;
  delete [] frequency;

  if (hnb_event) {
    for (i = htime->offset;i < htime->nb_value;i++) {
      if (htime->frequency[i] > 0) {
        delete hnb_event[i];
      }
    }
    delete [] hnb_event;
  }

  delete htime;
  delete mixture;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Destructor of the TimeEvents class.
 */
/*--------------------------------------------------------------*/

TimeEvents::~TimeEvents()

{
  remove();
}


/*--------------------------------------------------------------*/
/**
 *  \brief Assignment operator of the TimeEvents class.
 *
 *  \param[in] timev reference on a TimeEvents object.
 *
 *  \return          TimeEvents object.
 */
/*--------------------------------------------------------------*/

TimeEvents& TimeEvents::operator=(const TimeEvents &timev)

{
  if (&timev != this) {
    remove();
    copy(timev);
  }

  return *this;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Extraction of the number of events frequency distribution for a given observation period or
 *         of the mixture the number of events frequency distributions.
 *
 *  \param[in] error      reference on a StatError object,
 *  \param[in] histo_type frequency distribution type (NB_EVENT/NB_EVENT_MIXTURE),
 *  \param[in] itime      observation period.
 *
 *  \return               DiscreteDistributionData object.
 */
/*--------------------------------------------------------------*/

DiscreteDistributionData* TimeEvents::extract(StatError &error , renewal_distribution histo_type ,
                                              int itime) const

{
  DiscreteDistributionData *histo;


  error.init();

  if (histo_type == NB_EVENT) {
    if ((itime < htime->offset) || (itime >= htime->nb_value) || (htime->frequency[itime] == 0)) {
      histo = NULL;
      error.update(SEQ_error[SEQR_OBSERVATION_TIME]);
    }
    else {
      histo = new DiscreteDistributionData(*hnb_event[itime]);
    }
  }

  else if (histo_type == NB_EVENT_MIXTURE) {
    histo = new DiscreteDistributionData(*mixture);
  }

  return histo;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Change of the time unit of a TimeEvents object.
 *
 *  \param[in] error         reference on a StatError object,
 *  \param[in] scaling_coeff scaling factor.
 *
 *  \return                  TimeEvents object.
 */
/*--------------------------------------------------------------*/

TimeEvents* TimeEvents::time_scaling(StatError &error , int scaling_coeff) const

{
  bool status = true;
  int i;
  TimeEvents *timev;


  timev = NULL;
  error.init();

  if (scaling_coeff < 1) {
    status = false;
    error.update(STAT_error[STATR_SCALING_COEFF]);
  }
  if ((htime->nb_value - 1) * scaling_coeff > MAX_TIME) {
    status = false;
    error.update(SEQ_error[SEQR_LONG_OBSERVATION_TIME]);
  }

  if (status) {
    timev = new TimeEvents(nb_class);

    timev->nb_element = nb_element;

    for (i = 0;i < nb_class;i++) {
      timev->time[i] = time[i] * scaling_coeff;
      timev->nb_event[i] = nb_event[i];
      timev->frequency[i] = frequency[i];
    }

    timev->build_frequency_distribution();
  }

  return timev;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Selection of triplets {observation period, number of events, frequency} on
 *         an observation period criterion.
 *
 *  \param[in] error    reference on a StatError object,
 *  \param[in] min_time minimum observation period,
 *  \param[in] max_time maximum observation period.
 *
 *  \return             TimeEvents object.
 */
/*--------------------------------------------------------------*/

TimeEvents* TimeEvents::time_select(StatError &error , int min_time , int max_time) const

{
  bool status = true;
  int i , j;
  int bnb_class;
  TimeEvents *timev;


  timev = NULL;
  error.init();

  if ((min_time < 1) || (min_time > max_time)) {
    status = false;
    error.update(SEQ_error[SEQR_MIN_TIME]);
  }
  if ((max_time < htime->offset) || (max_time < min_time)) {
    status = false;
    error.update(SEQ_error[SEQR_MAX_TIME]);
  }

  if (status) {

    // computation of the number of classes

    bnb_class = 0;
    for (i = 0;i < nb_class;i++) {
      if ((time[i] >= min_time) && (time[i] <= max_time)) {
         bnb_class++;
      }
    }

    // copy of the selected triplets

    timev = new TimeEvents(bnb_class);

    i = 0;
    for (j = 0;j < nb_class;j++) {
      if ((time[j] >= min_time) && (time[j] <= max_time)) {
        timev->time[i] = time[j];
        timev->nb_event[i] = nb_event[j];
        timev->frequency[i] = frequency[j];
        i++;
      }
    }

    timev->nb_element_computation();

    if (timev->nb_element > 0) {
      timev->build_frequency_distribution();
    }

    else {
      delete timev;
      timev = NULL;
      error.update(SEQ_error[SEQR_EMPTY_RENEWAL_DATA]);
    }
  }

  return timev;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Selection of triplets {observation period, number of events, frequency} on
 *         a number of events criterion.
 *
 *  \param[in] error        reference on a StatError object,
 *  \param[in] min_nb_event minimum number of events,
 *  \param[in] max_nb_event maximum number of events.
 */
/*--------------------------------------------------------------*/

TimeEvents* TimeEvents::nb_event_select(StatError &error , int min_nb_event , int max_nb_event) const

{
  bool status = true;
  int i , j;
  int bnb_class;
  TimeEvents *timev;


  timev = NULL;
  error.init();

  if ((min_nb_event < 0) || (min_nb_event > max_nb_event)) {
    status = false;
    error.update(SEQ_error[SEQR_MIN_NB_EVENT]);
  }
  if ((max_nb_event < mixture->offset) || (max_nb_event < min_nb_event)) {
    status = false;
    error.update(SEQ_error[SEQR_MAX_NB_EVENT]);
  }

  if (status) {

    // computation of the number of classes

    bnb_class = 0;
    for (i = 0;i < nb_class;i++) {
      if ((nb_event[i] >= min_nb_event) && (nb_event[i] <= max_nb_event)) {
         bnb_class++;
      }
    }

    // copy of the selected triplets

    timev = new TimeEvents(bnb_class);

    i = 0;
    for (j = 0;j < nb_class;j++) {
      if ((nb_event[j] >= min_nb_event) && (nb_event[j] <= max_nb_event)) {
        timev->time[i] = time[j];
        timev->nb_event[i] = nb_event[j];
        timev->frequency[i] = frequency[j];
        i++;
      }
    }

    timev->nb_element_computation();

    if (timev->nb_element > 0) {
      timev->build_frequency_distribution();
    }

    else {
      delete timev;
      timev = NULL;
      error.update(SEQ_error[SEQR_EMPTY_RENEWAL_DATA]);
    }
  }

  return timev;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Construction of a TimeEvents object from a FrequencyDistribution object.
 *
 *  \param[in] error    reference on a StatError object,
 *  \param[in] nb_event reference on a FrequencyDistribution object,
 *  \param[in] itime    observation period.
 *
 *  \return             TimeEvents object.
 */
/*--------------------------------------------------------------*/

TimeEvents* TimeEvents::build(StatError &error , FrequencyDistribution &nb_event , int itime)

{
  bool status = true;
  TimeEvents *timev;


  timev = NULL;
  error.init();

  if ((nb_event.nb_value == 1) || (itime / (nb_event.nb_value - 1) < MIN_INTER_EVENT)) {
    status = false;
    error.update(SEQ_error[SEQR_SHORT_OBSERVATION_TIME]);
  }
  if (itime > MAX_TIME) {
    status = false;
    error.update(SEQ_error[SEQR_LONG_OBSERVATION_TIME]);
  }

  if (status) {
    timev = new TimeEvents(nb_event , itime);
  }

  return timev;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Construction of a TimeEvents object from pairs {observation period, number of events}.
 *
 *  \param[in] error    reference on a StatError object,
 *  \param[in] itime    pairs {observation period, number of events}.
 *
 *  \return             TimeEvents object.
 */
/*--------------------------------------------------------------*/

TimeEvents* TimeEvents::build(StatError &error , const vector<vector<int> > &time_nb_event)

{
  bool status = true;
  int i;
  int inb_element , *itime , *inb_event;
  TimeEvents *timev;


  timev = NULL;
  error.init();

  if (!time_nb_event.empty()) {
    inb_element = time_nb_event.size();

    for (i = 0;i < inb_element;i++) {
      if (time_nb_event[i].size() != 2) {
        status = false;
        error.update(SEQ_error[SEQR_TIME_NB_EVENT_PAIR] , i);
      }
    }
  }
  else {
    status = false;
    error.update(STAT_error[STATR_EMPTY_SAMPLE]);
  }
  
  if (status) {
    itime = new int[inb_element];
    inb_event = new int[inb_element];

    for (i = 0;i < inb_element;i++) {
      itime[i] = time_nb_event[i][0];
      inb_event[i] = time_nb_event[i][1];
    }

    timev = new TimeEvents();
    timev->build(inb_element , itime , inb_event);

    delete [] itime;
    delete [] inb_event;
  }

  return timev;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Construction of a TimeEvents object from a file whose format is:
 *         one triplet {observation period > 0, number of events >= 0, frequency >= 0} per line.
 *
 *  \param[in] error reference on a StatError object,
 *  \param[in] path  file path.
 *
 *  \return          TimeEvents object.
 */
/*--------------------------------------------------------------*/

TimeEvents* TimeEvents::ascii_read(StatError &error , const string path)

{
  string buffer;
  size_t position;
  typedef tokenizer<char_separator<char>> tokenizer;
  char_separator<char> separator(" \t");
  bool status , lstatus;
  int i , j;
  int line , nb_class , nb_element , value , time , nb_event;
  TimeEvents *timev;
  ifstream in_file(path.c_str());


  timev = NULL;
  error.init();

  if (!in_file) {
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {

    // 1st pass: analysis of each line format and search for
    // the number of triplets {observation period, number of events, frequency}

    status = true;
    line = 0;
    time = 0;
    nb_event = -1;
    nb_class = 0;
    nb_element = 0;

    while (getline(in_file , buffer)) {
      line++;

#     ifdef DEBUG
      cout << line << "  " << buffer << endl;
#     endif

      position = buffer.find('#');
      if (position != string::npos) {
        buffer.erase(position);
      }
      i = 0;

      tokenizer tok_buffer(buffer , separator);

      for (tokenizer::iterator token = tok_buffer.begin();token != tok_buffer.end();token++) {
        if (i <= 2) {
          lstatus = true;

/*          try {
            value = stoi(*token);   in C++ 11
          }
          catch(invalid_argument &arg) {
            lstatus = false;
          } */
          value = atoi(token->c_str());

          // test observation period > 0, number of events >= 0, frequency >= 0

          if ((lstatus) && (((i == 0) && (value <= 0)) ||
               ((i > 0) && (value < 0)))) {
            lstatus = false;
          }

          if (!lstatus) {
            status = false;
            error.update(STAT_parsing[STATP_DATA_TYPE] , line , i + 1);
          }

          else {
            switch (i) {

            // test ordered samples (observation period)

            case 0 : {
              if (value < time) {
                status = false;
                error.update(SEQ_parsing[SEQP_TIME_ORDER] , line , i + 1);
              }
              else if (value > time) {
                time = value;
                nb_event = -1;
              }
              if (value > MAX_TIME) {
                status = false;
                error.update(SEQ_parsing[SEQP_MAX_TIME] , line , i + 1);
              }
              break;
            }

            // test ordered samples (number of events)

            case 1 : {
              if (value <= nb_event) {
                status = false;
                error.update(SEQ_parsing[SEQP_NB_EVENT_ORDER] , line , i + 1);
              }
              else {
                nb_event = value;
              }
              break;
            }

            case 2 : {
              if (value > 0) {
                nb_class++;
                nb_element += value;
              }
              break;
            }
            }
          }
        }

        i++;
      }

      // test 3 items per line

      if ((i > 0) && (i != 3)) {
        status = false;
        error.correction_update(STAT_parsing[STATP_NB_TOKEN] , 3 , line);
      }
    }

    if (nb_element == 0) {
      status = false;
      error.update(STAT_parsing[STATP_EMPTY_SAMPLE]);
    }

    // 2nd pass: data reading

    if (status) {
//      in_file.close();
//      in_file.open(path.c_str() , ios::in);
      in_file.clear();
      in_file.seekg(0,ios::beg);

      timev = new TimeEvents(nb_class);
      timev->nb_element = nb_element;

      i = 0;

      while (getline(in_file , buffer)) {
        position = buffer.find('#');
        if (position != string::npos) {
          buffer.erase(position);
        }
        j = 0;

        tokenizer tok_buffer(buffer , separator);

        for (tokenizer::iterator token = tok_buffer.begin();token != tok_buffer.end();token++) {
          switch (j) {
          case 0 :
//            timev->time[i] = stoi(*token);   in C++ 11
            timev->time[i] = atoi(token->c_str());
            break;
          case 1 :
//            timev->nb_event[i] = stoi(*token);   in C++ 11
            timev->nb_event[i] = atoi(token->c_str());
            break;
          case 2 :
//            timev->frequency[i] = stoi(*token);   in C++ 11
            timev->frequency[i] = atoi(token->c_str());
            break;
          }

          j++;
        }

        if (timev->frequency[i] > 0) {
          i++;
        }
      }

      timev->build_frequency_distribution();
    }
  }

  return timev;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Construction of a TimeEvents object from a file whose format is:
 *         one pair {observation period > 0, number of events >= 0} per line.
 *
 *  \param[in] error reference on a StatError object,
 *  \param[in] path  file path.
 *
 *  \return          TimeEvents object.
 */
/*--------------------------------------------------------------*/

TimeEvents* TimeEvents::old_ascii_read(StatError &error , const string path)

{
  string buffer;
  size_t position;
  typedef tokenizer<char_separator<char>> tokenizer;
  char_separator<char> separator(" \t");
  bool status , lstatus;
  int i , j;
  int line , nb_element , value , *time , *nb_event;
  TimeEvents *timev;
  ifstream in_file(path.c_str());


  timev = NULL;
  error.init();

  if (!in_file) {
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {

    // 1st pass: analysis of each line format and search for
    // the number of pairs {observation period, number of events}

    status = true;
    line = 0;
    nb_element = 0;

    while (getline(in_file , buffer)) {
      line++;

#     ifdef DEBUG
      cout << line << "  " << buffer << endl;
#     endif

      position = buffer.find('#');
      if (position != string::npos) {
        buffer.erase(position);
      }
      i = 0;

      tokenizer tok_buffer(buffer , separator);

      for (tokenizer::iterator token = tok_buffer.begin();token != tok_buffer.end();token++) {
        if (i <= 1) {
          lstatus = true;

/*          try {
            value = stoi(*token);   in C++ 11
          }
          catch(invalid_argument &arg) {
            lstatus = false;
          } */
          value = atoi(token->c_str());

          // test observation period > 0, number of events >= 0

          if ((lstatus) && (((i == 0) && ((value <= 0) || (value > MAX_TIME))) ||
               ((i == 1) && (value < 0)))) {
            lstatus = false;
          }

          if (!lstatus) {
            status = false;
            error.update(STAT_parsing[STATP_DATA_TYPE] , line , i + 1);
          }
        }

        i++;
      }

      // test 2 items per line

      if (i > 0) {
        if (i != 2) {
          status = false;
          error.correction_update(STAT_parsing[STATP_NB_TOKEN] , 2 , line);
        }
        nb_element++;
      }
    }

    if (nb_element == 0) {
      status = false;
      error.update(STAT_parsing[STATP_FORMAT] , line);
    }

    // 2nd pass: data reading

    if (status) {
//      in_file.close();
//      in_file.open(path.c_str() , ios::in);
        in_file.clear();
        in_file.seekg(0,ios::beg);

      time = new int[nb_element];
      nb_event = new int[nb_element];
      i = 0;

      while (getline(in_file , buffer)) {
        position = buffer.find('#');
        if (position != string::npos) {
          buffer.erase(position);
        }
        j = 0;

        tokenizer tok_buffer(buffer , separator);

        for (tokenizer::iterator token = tok_buffer.begin();token != tok_buffer.end();token++) {
          switch (j) {
          case 0 :
//            time[i] = stoi(*token);   in C++ 11
            time[i] = atoi(token->c_str());
            break;
          case 1 :
//            nb_event[i] = stoi(*token);   in C++ 11
            nb_event[i] = atoi(token->c_str());
            break;
          }

          j++;
        }

        i++;
      }

      timev = new TimeEvents(nb_element , time , nb_event);

      delete [] time;
      delete [] nb_event;
    }
  }

  return timev;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing on a single line of a TimeEvents object.
 *
 *  \param[in,out] os stream.
 */
/*--------------------------------------------------------------*/

ostream& TimeEvents::line_write(ostream &os) const

{
  os << STAT_label[STATL_SAMPLE_SIZE] << ": " << nb_element << "   "
     << SEQ_label[SEQL_OBSERVATION_TIME] << " " << STAT_label[STATL_MEAN] << ": " << htime->mean << "   "
     << SEQ_label[SEQL_NB_EVENT] << " " << STAT_label[STATL_MEAN] << ": " << mixture->mean;

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a TimeEvents object.
 *
 *  \param[in,out] os         stream,
 *  \param[in]     exhaustive flag detail level,
 *  \param[in]     type       renewal process type (ORDINARY/EQUILIBRIUM).
 */
/*--------------------------------------------------------------*/

ostream& TimeEvents::ascii_write(ostream &os , bool exhaustive , process_type type) const

{
  int i;


  if ((htime->variance > 0.) && (exhaustive)) {
    os << SEQ_label[SEQL_OBSERVATION_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";
    htime->ascii_characteristic_print(os);

    os << "\n   | " << SEQ_label[SEQL_OBSERVATION_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << endl;
    htime->ascii_print(os);
  }

  for (i = htime->offset;i < htime->nb_value;i++) {
    if (htime->frequency[i] > 0) {
      if (((exhaustive) && (htime->variance > 0.)) || (i > htime->offset)) {
        os << "\n";
      }
      os << SEQ_label[SEQL_NB_EVENT] << " " << SEQ_label[SEQL_DURING] << " "
         << i << " " << SEQ_label[SEQL_TIME_UNIT] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";
      hnb_event[i]->ascii_characteristic_print(os);
      os << STAT_label[STATL_VARIANCE_MEAN_RATIO] << ": "
         << hnb_event[i]->variance / hnb_event[i]->mean << endl;
      if (hnb_event[i]->variance > 0.) {
        os << STAT_label[STATL_SKEWNESS_COEFF] << ": " << hnb_event[i]->skewness_computation() << "   "
           << STAT_label[STATL_KURTOSIS_COEFF] << ": " << hnb_event[i]->kurtosis_computation() << endl;
      }

      switch (type) {

      case ORDINARY : {
        os << "\n" << SEQ_label[SEQL_1_CENSORED_INTER_EVENT] << ": " << hnb_event[i]->nb_element
           << " (" << 1. / (hnb_event[i]->mean + 1.) << ")" << endl;

        os << SEQ_label[SEQL_COMPLETE_INTER_EVENT] << ": " << hnb_event[i]->mean *
              hnb_event[i]->nb_element << " ("
           << hnb_event[i]->mean / (hnb_event[i]->mean + 1.) << ")" << endl;
        break;
      }

      case EQUILIBRIUM : {
        os << "\n" << SEQ_label[SEQL_2_CENSORED_INTER_EVENT] << ": " << hnb_event[i]->frequency[0]
           << " (" << hnb_event[i]->frequency[0] / (hnb_event[i]->nb_element * (hnb_event[i]->mean + 1.))
           << ")" << endl;

        os << SEQ_label[SEQL_1_CENSORED_INTER_EVENT] << ": " << (hnb_event[i]->nb_element -
               hnb_event[i]->frequency[0]) * 2 << " ("
           << (hnb_event[i]->nb_element - hnb_event[i]->frequency[0]) * 2. /
              (hnb_event[i]->nb_element * (hnb_event[i]->mean + 1.)) << ")" << endl;

        os << SEQ_label[SEQL_COMPLETE_INTER_EVENT] << ": " << (hnb_event[i]->mean - 1.) *
              hnb_event[i]->nb_element + hnb_event[i]->frequency[0] << " ("
           << (hnb_event[i]->mean - 1. + (double)hnb_event[i]->frequency[0] /
               (double)hnb_event[i]->nb_element) / (hnb_event[i]->mean + 1.) << ")" << endl;
        break;
      }
      }

      if (exhaustive) {
        os << "\n   | " << SEQ_label[SEQL_NB_EVENT] << " " << SEQ_label[SEQL_DURING] << " "
           << i << " " << SEQ_label[SEQL_TIME_UNIT] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << endl;
        hnb_event[i]->ascii_print(os);
      }
    }
  }

  if ((htime->variance > 0.) && (exhaustive)) {
    os << "\n" << SEQ_label[SEQL_NB_EVENT] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";
    mixture->ascii_characteristic_print(os);

    switch (type) {

    case ORDINARY : {
      os << "\n" << SEQ_label[SEQL_1_CENSORED_INTER_EVENT] << ": " << mixture->nb_element
         << " (" << 1. / (mixture->mean + 1.) << ")" << endl;

      os << SEQ_label[SEQL_COMPLETE_INTER_EVENT] << ": " << mixture->mean *
            mixture->nb_element << " ("
         << mixture->mean / (mixture->mean + 1.) << ")" << endl;
      break;
    }

    case EQUILIBRIUM : {
      os << "\n" << SEQ_label[SEQL_2_CENSORED_INTER_EVENT] << ": " << mixture->frequency[0]
         << " (" << mixture->frequency[0] / (mixture->nb_element * (mixture->mean + 1.))
         << ")" << endl;

      os << SEQ_label[SEQL_1_CENSORED_INTER_EVENT] << ": " << (mixture->nb_element -
             mixture->frequency[0]) * 2 << " ("
         << (mixture->nb_element - mixture->frequency[0]) * 2. /
            (mixture->nb_element * (mixture->mean + 1.)) << ")" << endl;

      os << SEQ_label[SEQL_COMPLETE_INTER_EVENT] << ": " << (mixture->mean - 1.) *
            mixture->nb_element + mixture->frequency[0] << " ("
         << (mixture->mean - 1. + (double)mixture->frequency[0] /
             (double)mixture->nb_element) / (mixture->mean + 1.) << ")" << endl;
      break;
    }
    }

    os << "\n   | " << SEQ_label[SEQL_NB_EVENT] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << endl;
    mixture->ascii_print(os);
  }

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a TimeEvents object.
 *
 *  \param[in,out] os         stream,
 *  \param[in]     exhaustive flag detail level.
 */
/*--------------------------------------------------------------*/

ostream& TimeEvents::ascii_write(ostream &os , bool exhaustive) const

{
  return ascii_write(os , exhaustive , DEFAULT_TYPE);
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a TimeEvents object in a file.
 *
 *  \param[in,out] os         stream,
 *  \param[in]     exhaustive flag detail level,
 *  \param[in]     type       renewal process type (ORDINARY/EQUILIBRIUM).
 */
/*--------------------------------------------------------------*/

ostream& TimeEvents::ascii_file_write(ostream &os , bool exhaustive , process_type type) const

{
  int i;
  int max_frequency , width[3];


  if ((htime->variance > 0.) && (exhaustive)) {
    os << "# " << SEQ_label[SEQL_OBSERVATION_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";
    htime->ascii_characteristic_print(os , false , true);

    os << "\n#    | " << SEQ_label[SEQL_OBSERVATION_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << endl;
    htime->ascii_print(os , true);
  }

  // computation of the column widths

  if (exhaustive) {
    width[0] = column_width(time[nb_class - 1]);
    width[1] = column_width(nb_event[nb_class - 1]) + ASCII_SPACE;

    max_frequency = 0;
    for (i = 0;i < nb_class;i++) {
      if (frequency[i] > max_frequency) {
        max_frequency = frequency[i];
      }
    }

    width[2] = column_width(max_frequency) + ASCII_SPACE;
  }

  for (i = 0;i < nb_class;i++) {
    if ((i == 0) || ((i > 0) && (time[i] > time[i - 1]))) {
      if (((exhaustive) && (htime->variance > 0.)) || (i > 0)) {
        os << "\n";
      }
      if (exhaustive) {
        os << "# ";
      }
      os << SEQ_label[SEQL_NB_EVENT] << " " << SEQ_label[SEQL_DURING] << " "
         << time[i] << " " << SEQ_label[SEQL_TIME_UNIT] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";
      hnb_event[time[i]]->ascii_characteristic_print(os , false , exhaustive);
      if (exhaustive) {
        os << "# ";
      }
      os << STAT_label[STATL_VARIANCE_MEAN_RATIO] << ": "
         << hnb_event[time[i]]->variance / hnb_event[time[i]]->mean << endl;
      if (hnb_event[time[i]]->variance > 0.) {
        if (exhaustive) {
          os << "# ";
        }
        os << STAT_label[STATL_SKEWNESS_COEFF] << ": " << hnb_event[time[i]]->skewness_computation() << "   "
           << STAT_label[STATL_KURTOSIS_COEFF] << ": " << hnb_event[time[i]]->kurtosis_computation() << endl;
      }

      switch (type) {

      case ORDINARY : {
        os << "\n# " << SEQ_label[SEQL_1_CENSORED_INTER_EVENT] << ": " << hnb_event[time[i]]->nb_element
           << " (" << 1. / (hnb_event[time[i]]->mean + 1.) << ")" << endl;

        os << "# " << SEQ_label[SEQL_COMPLETE_INTER_EVENT] << ": " << hnb_event[time[i]]->mean *
              hnb_event[time[i]]->nb_element << " ("
           << hnb_event[time[i]]->mean / (hnb_event[time[i]]->mean + 1.) << ")" << endl;
        break;
      }

      case EQUILIBRIUM : {
        os << "\n# " << SEQ_label[SEQL_2_CENSORED_INTER_EVENT] << ": " << hnb_event[time[i]]->frequency[0]
           << " (" << hnb_event[time[i]]->frequency[0] / (hnb_event[time[i]]->nb_element * (hnb_event[time[i]]->mean + 1.))
           << ")" << endl;

        os << "# " << SEQ_label[SEQL_1_CENSORED_INTER_EVENT] << ": " << (hnb_event[time[i]]->nb_element -
               hnb_event[time[i]]->frequency[0]) * 2 << " ("
           << (hnb_event[time[i]]->nb_element - hnb_event[time[i]]->frequency[0]) * 2. /
              (hnb_event[time[i]]->nb_element * (hnb_event[time[i]]->mean + 1.)) << ")" << endl;

        os << "# " << SEQ_label[SEQL_COMPLETE_INTER_EVENT] << ": " << (hnb_event[time[i]]->mean - 1.) *
              hnb_event[time[i]]->nb_element + hnb_event[time[i]]->frequency[0] << " ("
           << (hnb_event[time[i]]->mean - 1. + (double)hnb_event[time[i]]->frequency[0] /
               (double)hnb_event[time[i]]->nb_element) / (hnb_event[time[i]]->mean + 1.) << ")" << endl;
        break;
      }
      }

      if (exhaustive) {
        os << "\n";
      }
    }

    // writing of the triplets (observation period, number of events, frequency)

    if (exhaustive) {
      os << setw(width[0]) << time[i];
      os << setw(width[1]) << nb_event[i];
      os << setw(width[2]) << frequency[i] << endl;
    }
  }

  if ((htime->variance > 0.) && (exhaustive)) {
    os << "\n# " << SEQ_label[SEQL_NB_EVENT] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";
    mixture->ascii_characteristic_print(os , false , true);

    switch (type) {

    case ORDINARY : {
      os << "\n# " << SEQ_label[SEQL_1_CENSORED_INTER_EVENT] << ": " << mixture->nb_element
         << " (" << 1. / (mixture->mean + 1.) << ")" << endl;

      os << "# " << SEQ_label[SEQL_COMPLETE_INTER_EVENT] << ": " << mixture->mean *
            mixture->nb_element << " ("
         << mixture->mean / (mixture->mean + 1.) << ")" << endl;
      break;
    }

    case EQUILIBRIUM : {
      os << "\n# " << SEQ_label[SEQL_2_CENSORED_INTER_EVENT] << ": " << mixture->frequency[0]
         << " (" << mixture->frequency[0] / (mixture->nb_element * (mixture->mean + 1.))
         << ")" << endl;

      os << "# " << SEQ_label[SEQL_1_CENSORED_INTER_EVENT] << ": " << (mixture->nb_element -
             mixture->frequency[0]) * 2 << " ("
         << (mixture->nb_element - mixture->frequency[0]) * 2. /
            (mixture->nb_element * (mixture->mean + 1.)) << ")" << endl;

      os << "# " << SEQ_label[SEQL_COMPLETE_INTER_EVENT] << ": " << (mixture->mean - 1.) *
            mixture->nb_element + mixture->frequency[0] << " ("
         << (mixture->mean - 1. + (double)mixture->frequency[0] /
             (double)mixture->nb_element) / (mixture->mean + 1.) << ")" << endl;
      break;
    }
    }

    os << "\n#    | " << SEQ_label[SEQL_NB_EVENT] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << endl;
    mixture->ascii_print(os , true);
  }

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a TimeEvents object in a file.
 *
 *  \param[in] error      reference on a StatError object,
 *  \param[in] path       file path,
 *  \param[in] exhaustive flag detail level.
 *
 *  \return               error status.
 */
/*--------------------------------------------------------------*/

bool TimeEvents::ascii_write(StatError &error , const string path ,
                             bool exhaustive) const

{
  bool status;
  ofstream out_file(path.c_str());


  error.init();

  if (!out_file) {
    status = false;
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {
    status = true;
    ascii_file_write(out_file , exhaustive);
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a TimeEvents object at the spreadsheet format.
 *
 *  \param[in,out] os   stream,
 *  \param[in]     type renewal process type (ORDINARY/EQUILIBRIUM).
 */
/*--------------------------------------------------------------*/

ostream& TimeEvents::spreadsheet_write(ostream &os , process_type type) const

{
  int i;


  if (htime->variance > 0.) {
    os << SEQ_label[SEQL_OBSERVATION_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\t";
    htime->spreadsheet_characteristic_print(os);

    os << "\n\t" << SEQ_label[SEQL_OBSERVATION_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << endl;
    htime->spreadsheet_print(os);
  }

  for (i = htime->offset;i < htime->nb_value;i++) {
    if (htime->frequency[i] > 0) {
      if ((htime->variance > 0.) || (i > htime->offset)) {
        os << "\n";
      }
      os << SEQ_label[SEQL_NB_EVENT] << " " << SEQ_label[SEQL_DURING] << " "
         << i << " " << SEQ_label[SEQL_TIME_UNIT] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\t";
      hnb_event[i]->spreadsheet_characteristic_print(os);
      os << STAT_label[STATL_VARIANCE_MEAN_RATIO] << "\t"
         << hnb_event[i]->variance / hnb_event[i]->mean << endl;
      if (hnb_event[i]->variance > 0.) {
        os << STAT_label[STATL_SKEWNESS_COEFF] << "\t" << hnb_event[i]->skewness_computation() << "\t"
           << STAT_label[STATL_KURTOSIS_COEFF] << "\t" << hnb_event[i]->kurtosis_computation() << endl;
      }

      switch (type) {

      case ORDINARY : {
        os << "\n" << SEQ_label[SEQL_1_CENSORED_INTER_EVENT] << "\t" << hnb_event[i]->nb_element
           << "\t" << 1. / (hnb_event[i]->mean + 1.) << endl;

        os << SEQ_label[SEQL_COMPLETE_INTER_EVENT] << "\t" << hnb_event[i]->mean *
              hnb_event[i]->nb_element << "\t"
           << hnb_event[i]->mean / (hnb_event[i]->mean + 1.) << endl;
        break;
      }

      case EQUILIBRIUM : {
        os << "\n" << SEQ_label[SEQL_2_CENSORED_INTER_EVENT] << "\t" << hnb_event[i]->frequency[0] << "\t"
           << hnb_event[i]->frequency[0] / (hnb_event[i]->nb_element * (hnb_event[i]->mean + 1.)) << endl;

        os << SEQ_label[SEQL_1_CENSORED_INTER_EVENT] << "\t" << (hnb_event[i]->nb_element -
               hnb_event[i]->frequency[0]) * 2 << "\t"
           << (hnb_event[i]->nb_element - hnb_event[i]->frequency[0]) * 2. /
              (hnb_event[i]->nb_element * (hnb_event[i]->mean + 1.)) << endl;

        os << SEQ_label[SEQL_COMPLETE_INTER_EVENT] << "\t" << (hnb_event[i]->mean - 1.) *
              hnb_event[i]->nb_element + hnb_event[i]->frequency[0] << "\t"
           << (hnb_event[i]->mean - 1. + (double)hnb_event[i]->frequency[0] /
               (double)hnb_event[i]->nb_element) / (hnb_event[i]->mean + 1.) << endl;
        break;
      }
      }

      os << "\n\t" << SEQ_label[SEQL_NB_EVENT] << " " << SEQ_label[SEQL_DURING] << " "
         << i << " " << SEQ_label[SEQL_TIME_UNIT] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << endl;
      hnb_event[i]->spreadsheet_print(os);
    }
  }

  if (htime->variance > 0.) {
    os << "\n" << SEQ_label[SEQL_NB_EVENT] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\t";
    mixture->spreadsheet_characteristic_print(os);

    switch (type) {

    case ORDINARY : {
      os << "\n" << SEQ_label[SEQL_1_CENSORED_INTER_EVENT] << "\t" << mixture->nb_element
         << "\t" << 1. / (mixture->mean + 1.) << endl;

      os << SEQ_label[SEQL_COMPLETE_INTER_EVENT] << "\t" << mixture->mean *
            mixture->nb_element << "\t"
         << mixture->mean / (mixture->mean + 1.) << endl;
      break;
    }

    case EQUILIBRIUM : {
      os << "\n" << SEQ_label[SEQL_2_CENSORED_INTER_EVENT] << "\t" << mixture->frequency[0] << "\t"
         << mixture->frequency[0] / (mixture->nb_element * (mixture->mean + 1.)) << endl;

      os << SEQ_label[SEQL_1_CENSORED_INTER_EVENT] << "\t" << (mixture->nb_element -
             mixture->frequency[0]) * 2 << "\t"
         << (mixture->nb_element - mixture->frequency[0]) * 2. /
            (mixture->nb_element * (mixture->mean + 1.)) << endl;

      os << SEQ_label[SEQL_COMPLETE_INTER_EVENT] << "\t" << (mixture->mean - 1.) *
            mixture->nb_element + mixture->frequency[0] << "\t"
         << (mixture->mean - 1. + (double)mixture->frequency[0] /
             (double)mixture->nb_element) / (mixture->mean + 1.) << endl;
      break;
    }
    }

    os << "\n\t" << SEQ_label[SEQL_NB_EVENT] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << endl;
    mixture->spreadsheet_print(os);
  }

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a TimeEvents object in a file at the spreadsheet format.
 *
 *  \param[in] error reference on a StatError object,
 *  \param[in] path  file path.
 *
 *  \return          error status.
 */
/*--------------------------------------------------------------*/

bool TimeEvents::spreadsheet_write(StatError &error , const string path) const

{
  bool status;
  ofstream out_file(path.c_str());


  error.init();

  if (!out_file) {
    status = false;
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {
    status = true;
    spreadsheet_write(out_file);
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Plot of a TimeEvents object using Gnuplot.
 *
 *  \param[in] error  reference on a StatError object,
 *  \param[in] prefix file prefix,
 *  \param[in] title  figure title.
 *
 *  \return           error status.
 */
/*--------------------------------------------------------------*/

bool TimeEvents::plot_write(StatError &error , const char *prefix ,
                            const char *title) const

{
  bool status;
  int i , j , k;
  int nb_histo;
  const FrequencyDistribution **phisto;
  ostringstream data_file_name;


  error.init();

  // writing of data file

  data_file_name << prefix << ".dat";

  nb_histo = 0;
  for (i = htime->offset + 1;i < htime->nb_value;i++) {
    if (htime->frequency[i] > 0) {
      nb_histo++;
    }
  }
  if (htime->variance > 0.) {
    nb_histo += 2;
  }

  phisto = new const FrequencyDistribution*[nb_histo];

  nb_histo = 0;
  if (htime->variance > 0.) {
    phisto[nb_histo++] = htime;
  }
  for (i = htime->offset + 1;i < htime->nb_value;i++) {
    if (htime->frequency[i] > 0) {
      phisto[nb_histo++] = hnb_event[i];
    }
  }
  if (htime->variance > 0.) {
    phisto[nb_histo++] = mixture;
  }

  status = hnb_event[htime->offset]->plot_print((data_file_name.str()).c_str() , nb_histo , phisto);

  delete [] phisto;

  if (!status) {
    error.update(STAT_error[STATR_FILE_PREFIX]);
  }

  // writing of script files

  else {
    for (i = 0;i < 2;i++) {
      j = 1;

      ostringstream file_name[2];

      switch (i) {
      case 0 :
        file_name[0] << prefix << ".plot";
        break;
      case 1 :
        file_name[0] << prefix << ".print";
        break;
      }

      ofstream out_file((file_name[0].str()).c_str());

      if (i == 1) {
        out_file << "set terminal postscript" << endl;
        file_name[1] << label(prefix) << ".ps";
        out_file << "set output \"" << file_name[1].str() << "\"\n\n";
      }

      out_file << "set border 15 lw 0\n" << "set tics out\n" << "set xtics nomirror\n"
               << "set title";
      if (title) {
        out_file << " \"" << title << "\"";
      }
      out_file << "\n\n";

      if (htime->variance > 0.) {
        if (htime->nb_value - 1 < TIC_THRESHOLD) {
          out_file << "set xtics 0,1" << endl;
        }
        if ((int)(htime->max * YSCALE) + 1 < TIC_THRESHOLD) {
          out_file << "set ytics 0,1" << endl;
        }

        out_file << "plot [0:" << htime->nb_value - 1 << "] [0:"
                 << (int)(htime->max * YSCALE) + 1 << "] \""
                 << label((data_file_name.str()).c_str()) << "\" using " << j++
                 << " title \"" << SEQ_label[SEQL_OBSERVATION_TIME] << " "
                 << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\" with impulses" << endl;

        if (htime->nb_value - 1 < TIC_THRESHOLD) {
          out_file << "set xtics autofreq" << endl;
        }
        if ((int)(htime->max * YSCALE) + 1 < TIC_THRESHOLD) {
          out_file << "set ytics autofreq" << endl;
        }
      }

      for (k = htime->offset;k < htime->nb_value;k++) {
        if (htime->frequency[k] > 0) {
          if (((htime->variance > 0.) || (k > htime->offset)) && (i == 0)) {
            out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
          }
          out_file << endl;

          if (hnb_event[k]->nb_value - 1 < TIC_THRESHOLD) {
            out_file << "set xtics 0,1" << endl;
          }
          if ((int)(hnb_event[k]->max * YSCALE) + 1 < TIC_THRESHOLD) {
            out_file << "set ytics 0,1" << endl;
          }

          out_file << "plot [0:" << hnb_event[k]->nb_value - 1 << "] [0:"
                   << (int)(hnb_event[k]->max * YSCALE) + 1 << "] \""
                   << label((data_file_name.str()).c_str()) << "\" using " << j++
                   << " title \"" << SEQ_label[SEQL_NB_EVENT] << " "
                   << SEQ_label[SEQL_DURING] << " " << k << " " << SEQ_label[SEQL_TIME_UNIT] << " "
                   << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\" with impulses" << endl;

          if (hnb_event[k]->nb_value - 1 < TIC_THRESHOLD) {
            out_file << "set xtics autofreq" << endl;
          }
          if ((int)(hnb_event[k]->max * YSCALE) + 1 < TIC_THRESHOLD) {
            out_file << "set ytics autofreq" << endl;
          }
        }
      }

      if (htime->variance > 0.) {
        if (i == 0) {
          out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
        }
        out_file << endl;

        if (mixture->nb_value - 1 < TIC_THRESHOLD) {
          out_file << "set xtics 0,1" << endl;
        }
        if ((int)(mixture->max * YSCALE) + 1 < TIC_THRESHOLD) {
          out_file << "set ytics 0,1" << endl;
        }

        out_file << "plot [0:" << mixture->nb_value - 1 << "] [0:"
                 << (int)(mixture->max * YSCALE) + 1 << "] \""
                 << label((data_file_name.str()).c_str()) << "\" using " << j
                 << " title \"" << SEQ_label[SEQL_NB_EVENT] << " "
                 << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\" with impulses" << endl;

        if (mixture->nb_value - 1 < TIC_THRESHOLD) {
          out_file << "set xtics autofreq" << endl;
        }
        if ((int)(mixture->max * YSCALE) + 1 < TIC_THRESHOLD) {
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


/*--------------------------------------------------------------*/
/**
 *  \brief Plot of a TimeEvents object.
 *
 *  \return MultiPlotSet object.
 */
/*--------------------------------------------------------------*/

MultiPlotSet* TimeEvents::get_plotable() const

{
  int i , j , k , m;
  int nb_plot_set , nb_histo , max_nb_value , max_frequency;
  double shift;
  const FrequencyDistribution *phisto[2] , **merged_histo;
  ostringstream legend;
  MultiPlotSet *plot_set;


  nb_plot_set = 1;
  if (htime->variance > 0.) {
    nb_plot_set += 2;
  }

  plot_set = new MultiPlotSet(nb_plot_set);
  MultiPlotSet &plot = *plot_set;

  plot.border = "15 lw 0";

  i = 0;
  if (htime->variance > 0.) {

    // observation period frequency distribution

    plot[i].xrange = Range(0 , htime->nb_value - 1);
    plot[i].yrange = Range(0 , ceil(htime->max * YSCALE));

    if (htime->nb_value - 1 < TIC_THRESHOLD) {
      plot[i].xtics = 1;
    }
    if (ceil(htime->max * YSCALE) < TIC_THRESHOLD) {
      plot[i].ytics = 1;
    }

    plot[i].resize(1);

    legend.str("");
    legend << SEQ_label[SEQL_OBSERVATION_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
    plot[i][0].legend = legend.str();

    plot[i][0].style = "impulses";

    htime->plotable_frequency_write(plot[i][0]);
    i++;
  }

  // number of events frequency distribution for each observation period

  nb_histo = 0;
  max_nb_value = 0;
  max_frequency = 0;

  for (j = htime->offset;j < htime->nb_value;j++) {
    if (htime->frequency[j] > 0) {
      nb_histo++;

      // computation of the maximum number of values and the maximum frequency

      if (hnb_event[j]->nb_value > max_nb_value) {
        max_nb_value = hnb_event[j]->nb_value;
      }
      if (hnb_event[j]->max > max_frequency) {
        max_frequency = hnb_event[j]->max;
      }
    }
  }

  plot[i].xrange = Range(0 , max_nb_value);
  plot[i].yrange = Range(0 , ceil(max_frequency * YSCALE));

  if (max_nb_value < TIC_THRESHOLD) {
    plot[i].xtics = 1;
  }
  if (ceil(max_frequency * YSCALE) < TIC_THRESHOLD) {
    plot[i].ytics = 1;
  }

  plot[i].resize(nb_histo);

  j = 0;
  shift = 0.;

  for (k = htime->offset;k < htime->nb_value;k++) {
    if (htime->frequency[k] > 0) {
      legend.str("");
      legend << SEQ_label[SEQL_NB_EVENT] << " " << SEQ_label[SEQL_DURING] << " " << k << " "
             << SEQ_label[SEQL_TIME_UNIT] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
      plot[i][j].legend = legend.str();

      plot[i][j].style = "impulses";

      for (m = hnb_event[k]->offset;m < hnb_event[k]->nb_value;m++) {
        if (hnb_event[k]->frequency[m] > 0) {
          plot[i][j].add_point(m + shift , hnb_event[k]->frequency[m]);
        }
      }

      if (PLOT_SHIFT * (nb_histo - 1) < PLOT_MAX_SHIFT) {
        shift += PLOT_SHIFT;
      }
      else {
        shift += PLOT_MAX_SHIFT / (nb_histo - 1);
      }

      j++;
    }
  }

  if (htime->variance > 0.) {
    i++;

    // superimposed number of events frequency distributions

    merged_histo = new const FrequencyDistribution*[nb_histo];

    j = nb_histo - 1;
    for (k = htime->nb_value - 1;k >= htime->offset;k--) {
      if (htime->frequency[k] > 0) {
        if (j == nb_histo - 1) {
          merged_histo[j] = new FrequencyDistribution(*hnb_event[k]);
        }

        else {
          phisto[0] = merged_histo[j + 1];
          phisto[1] = hnb_event[k];
          merged_histo[j] = new FrequencyDistribution(2 , phisto);
        }

        j--;
      }
    }

    plot[i].xrange = Range(0 , merged_histo[0]->nb_value - 1);
    plot[i].yrange = Range(0 , ceil(merged_histo[0]->max * YSCALE));

    if (merged_histo[0]->nb_value - 1 < TIC_THRESHOLD) {
      plot[i].xtics = 1;
    }
    if (ceil(merged_histo[0]->max * YSCALE) < TIC_THRESHOLD) {
      plot[i].ytics = 1;
    }

    plot[i].resize(nb_histo);

    j = 0;
    for (k = htime->offset;k < htime->nb_value;k++) {
      if (htime->frequency[k] > 0) {
        legend.str("");
        legend << SEQ_label[SEQL_NB_EVENT] << " " << SEQ_label[SEQL_DURING] << " " << k << " "
               << SEQ_label[SEQL_TIME_UNIT] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
        plot[i][j].legend = legend.str();

        plot[i][j].style = "impulses";

        merged_histo[j]->plotable_frequency_write(plot[i][j]);
        j++;
      }
    }

    for (j = 0;j < nb_histo;j++) {
      delete merged_histo[j];
    }
    delete [] merged_histo;

/*    plot[i].resize(1);

    legend.str("");
    legend << SEQ_label[SEQL_NB_EVENT] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
    plot[i][0].legend = legend.str();

    plot[i][0].style = "impulses";

    mixture->plotable_frequency_write(plot[i][0]); */
  }

  return plot_set;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Search for the minimum mean time interval betweeen 2 events.
 *
 *  \return minimum mean time interval betweeen 2 events.
 */
/*--------------------------------------------------------------*/

double TimeEvents::min_inter_event_computation() const

{
  int i;
  double ratio , min_ratio;


  min_ratio = time[nb_class - 1];
  for (i = 0;i < nb_class;i++) {
    if (nb_event[i] > 0) {
      ratio = (double)time[i] / (double)nb_event[i];
      if (ratio < min_ratio) {
        min_ratio = ratio;
      }
    }
  }

  return min_ratio;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the sample size of a TimeEvents object.
 */
/*--------------------------------------------------------------*/

void TimeEvents::nb_element_computation()

{
  int i;


  nb_element = 0;
  for (i = 0;i < nb_class;i++) {
    nb_element += frequency[i];
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Default constructor of the RenewalData class.
 */
/*--------------------------------------------------------------*/

RenewalData::RenewalData()

{
  renewal = NULL;

  type = EQUILIBRIUM;

  length = NULL;
  sequence = NULL;

  inter_event = NULL;
  within = NULL;
  length_bias = NULL;
  backward = NULL;
  forward = NULL;

  index_event = NULL;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Constructor of the RenewalData class.
 *
 *  \param[in] nb_element sample size,
 *  \param[in] itime      observation period.
 */
/*--------------------------------------------------------------*/

RenewalData::RenewalData(int nb_element , int itime)

{
  renewal = NULL;

  type = EQUILIBRIUM;

  length = new int[nb_element];
  sequence = new int*[nb_element];

  inter_event = NULL;
  within = new FrequencyDistribution(itime);
  length_bias = NULL;
  backward = new FrequencyDistribution(itime);
  forward = new FrequencyDistribution(itime + 1);

  index_event = NULL;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Construction of a RenewalData object from a TimeEvents object.
 *
 *  \param[in] timev reference on a TimeEvents object,
 *  \param[in] itype renewal process type (ORDINARY/EQUILIBRIUM).
 */
/*--------------------------------------------------------------*/

RenewalData::RenewalData(const TimeEvents &timev , process_type itype)
:TimeEvents(timev)

{
  renewal = NULL;

  type = itype;

  length = NULL;
  sequence = NULL;

  inter_event = NULL;
  within = NULL;
  length_bias = NULL;
  backward = NULL;
  forward = NULL;

  index_event = NULL;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Construction of a RenewalData object from a Renewal object.
 *
 *  \param[in] itype renewal process type (ORDINARY/EQUILIBRIUM),
 *  \param[in] renew reference on a Renewal object.
 */
/*--------------------------------------------------------------*/

RenewalData::RenewalData(process_type itype , const Renewal &renew)

{
  renewal = NULL;

  type = itype;

  length = NULL;
  sequence = NULL;

  inter_event = new FrequencyDistribution(*(renew.inter_event));
  within = new FrequencyDistribution(*(renew.inter_event));
  length_bias = new FrequencyDistribution(*(renew.length_bias));
  backward = new FrequencyDistribution(renew.backward->alloc_nb_value);
  forward = new FrequencyDistribution(*(renew.forward));

  index_event = NULL;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Constructor by merging of the RenewalData class.
 *
 *  \param[in] nb_sample number of RenewalData objects,
 *  \param[in] itimev    pointer on the RenewalData objects.
 */
/*--------------------------------------------------------------*/

RenewalData::RenewalData(int nb_sample , const RenewalData **itimev)

{
  int i , j , k;
  const TimeEvents **ptimev;


  ptimev = new const TimeEvents*[nb_sample];

  for (i = 0;i < nb_sample;i++) {
    ptimev[i] = itimev[i];
  }
  TimeEvents::merge(nb_sample , ptimev);

  delete [] ptimev;

  renewal = NULL;

  type = itimev[0]->type;

  length = new int[nb_element];
  sequence = new int*[nb_element];

  i = 0;
  for (j = 0;j < nb_sample;j++) {
    for (k = 0;k < itimev[j]->nb_element;k++) {
      length[i] = itimev[j]->length[k];
      sequence[i] = new int[length[i]];
      i++;
    }
  }

  inter_event = NULL;
  within = NULL;
  length_bias = NULL;
  backward = NULL;
  forward = NULL;

  index_event = NULL;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Copy of a RenewalData object.
 *
 *  \param[in] timev      reference on a RenewalData object,
 *  \param[in] model_flag flag copy of the Renewal object.
 */
/*--------------------------------------------------------------*/

void RenewalData::copy(const RenewalData &timev , bool model_flag)

{
  int i , j;


  if ((model_flag) && (timev.renewal)) {
    renewal = new Renewal(*(timev.renewal) , false);
  }
  else {
    renewal = NULL;
  }

  type = timev.type;

  length = new int[nb_element];
  for (i = 0;i < nb_element;i++) {
    length[i] = timev.length[i];
  }

  sequence = new int*[nb_element];
  for (i = 0;i < nb_element;i++) {
    sequence[i] = new int[length[i]];
    for (j = 0;j < length[i];j++) {
      sequence[i][j] = timev.sequence[i][j];
    }
  }

  if (timev.inter_event) {
    inter_event = new FrequencyDistribution(*(timev.inter_event));
  }
  else {
    inter_event = NULL;
  }
  within = new FrequencyDistribution(*(timev.within));
  if (timev.length_bias) {
    length_bias = new FrequencyDistribution(*(timev.length_bias));
  }
  else {
    length_bias = NULL;
  }
  backward = new FrequencyDistribution(*(timev.backward));
  forward = new FrequencyDistribution(*(timev.forward));

  index_event = new Curves(*(timev.index_event));
}


/*--------------------------------------------------------------*/
/**
 *  \brief Destruction of the data members of a RenewalData object.
 */
/*--------------------------------------------------------------*/

void RenewalData::remove()

{
  delete renewal;

  delete [] length;

  if (sequence) {
    int i;

    for (i = 0;i < nb_element;i++) {
      delete [] sequence[i];
    }
    delete [] sequence;
  }

  delete inter_event;
  delete within;
  delete length_bias;
  delete backward;
  delete forward;

  delete index_event;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Destructor of the RenewalData class.
 */
/*--------------------------------------------------------------*/

RenewalData::~RenewalData()

{
  remove();
}


/*--------------------------------------------------------------*/
/**
 *  \brief Assignment operator of the RenewalData class.
 *
 *  \param[in] timev reference on a RenewalData object.
 *
 *  \return          RenewalData object.
 */
/*--------------------------------------------------------------*/

RenewalData& RenewalData::operator=(const RenewalData &timev)

{
  if (&timev != this) {
    remove();
    TimeEvents::remove();

    TimeEvents::copy(timev);
    copy(timev);
  }

  return *this;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Merging of RenewalData objects.
 *
 *  \param[in] error     reference on a StatError object,
 *  \param[in] nb_sample number of RenewalData objects,
 *  \param[in] itimev    pointer on the RenewalData objects.
 *
 *  \return              RenewalData object.
 */
/*--------------------------------------------------------------*/

RenewalData* RenewalData::merge(StatError &error , int nb_sample ,
                                const RenewalData **itimev) const

{
  bool status = true;
  int i , j , k , m;
  const FrequencyDistribution **phisto;
  RenewalData *timev;
  const RenewalData **ptimev;


  timev = NULL;
  error.init();

  for (i = 0;i < nb_sample;i++) {
    if ((itimev[i]->type != type) || ((itimev[i]->inter_event) && (!inter_event)) ||
        ((!(itimev[i]->inter_event)) && (inter_event))) {
      status = false;
      ostringstream error_message;
      error_message << STAT_label[STATL_SAMPLE] << " " << i + 2 << ": "
                    << SEQ_error[SEQR_INCOMPATIBLE_RENEWAL_DATA];
      error.update((error_message.str()).c_str());
    }
  }

  if (status) {
    nb_sample++;
    ptimev = new const RenewalData*[nb_sample];

    ptimev[0] = this;
    for (i = 1;i < nb_sample;i++) {
      ptimev[i] = itimev[i - 1];
    }

    timev = new RenewalData(nb_sample , ptimev);

    // copy of the sequences of events

    i = 0;
    for (j = 0;j < nb_sample;j++) {
      for (k = 0;k < ptimev[j]->nb_element;k++) {
        for (m = 0;m < ptimev[j]->length[k];m++) {
          timev->sequence[i][m] = ptimev[j]->sequence[k][m];
        }
        i++;
      }
    }

    phisto = new const FrequencyDistribution*[nb_sample];

    if (inter_event) {
      for (i = 0;i < nb_sample;i++) {
        phisto[i] = ptimev[i]->inter_event;
      }
      timev->inter_event = new FrequencyDistribution(nb_sample , phisto);
    }

    for (i = 0;i < nb_sample;i++) {
      phisto[i] = ptimev[i]->within;
    }
    timev->within = new FrequencyDistribution(nb_sample , phisto);

    if (length_bias) {
      for (i = 0;i < nb_sample;i++) {
        phisto[i] = ptimev[i]->length_bias;
      }
      timev->length_bias = new FrequencyDistribution(nb_sample , phisto);
    }

    for (i = 0;i < nb_sample;i++) {
      phisto[i] = ptimev[i]->backward;
    }
    timev->backward = new FrequencyDistribution(nb_sample , phisto);

    for (i = 0;i < nb_sample;i++) {
      phisto[i] = ptimev[i]->forward;
    }
    timev->forward = new FrequencyDistribution(nb_sample , phisto);

    timev->build_index_event(timev->type == ORDINARY ? 0 : 1);

    delete [] phisto;
    delete [] ptimev;
  }

  return timev;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Merging of RenewalData objects.
 *
 *  \param[in] error     reference on a StatError object,
 *  \param[in] nb_sample number of RenewalData objects,
 *  \param[in] itimev    RenewalData objects.
 *
 *  \return              RenewalData object.
 */
/*--------------------------------------------------------------*/

RenewalData* RenewalData::merge(StatError &error , int nb_sample ,
                                const vector<RenewalData> &itimev) const

{
  int i;
  RenewalData *timev;
  const RenewalData **ptimev;


  ptimev = new const RenewalData*[nb_sample];
  for (i = 0;i < nb_sample;i++) {
    ptimev[i] = new RenewalData(itimev[i]);
  }

  timev = merge(error , nb_sample , ptimev);

  for (i = 0;i < nb_sample;i++) {
    delete ptimev[i];
  }
  delete [] ptimev;

  return timev;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Extraction of a frequency distribution.
 *
 *  \param[in] error      reference on a StatError object,
 *  \param[in] histo_type frequency distribution type,
 *  \param[in] itime      observation period.
 *
 *  \return               DiscreteDistributionData object.
 */
/*--------------------------------------------------------------*/

DiscreteDistributionData* RenewalData::extract(StatError &error , renewal_distribution histo_type ,
                                               int itime) const

{
  bool status = true;
  Distribution *pdist;
  DiscreteParametric *pparam;
  FrequencyDistribution *phisto;
  DiscreteDistributionData *histo;


  error.init();

  if (histo_type == NB_EVENT) {
    if ((itime < htime->offset) || (itime >= htime->nb_value) || (htime->frequency[itime] == 0)) {
      histo = NULL;
      error.update(SEQ_error[SEQR_OBSERVATION_TIME]);
    }
    else {
      histo = new DiscreteDistributionData(*hnb_event[itime] ,
                                           (renewal ? renewal->nb_event[itime] : NULL));
    }
  }

  else {
    histo = NULL;

    switch (histo_type) {

    case INTER_EVENT : {
      if (inter_event) {
        phisto = inter_event;
      }
      else {
        status = false;
        ostringstream error_message;
        error_message << SEQ_label[SEQL_INTER_EVENT] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " "
                      << STAT_error[STATR_NOT_PRESENT];
        error.update((error_message.str()).c_str());
      }
      break;
    }

    case WITHIN_OBSERVATION_PERIOD : {
      phisto = within;
      break;
    }

    case LENGTH_BIAS : {
      if (length_bias) {
        phisto = length_bias;
      }
      else {
        status = false;
        ostringstream error_message;
        error_message << SEQ_label[SEQL_LENGTH_BIASED] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " "
                      << STAT_error[STATR_NOT_PRESENT];
        error.update((error_message.str()).c_str());
      }
      break;
    }

    case BACKWARD_RECURRENCE_TIME : {
      phisto = backward;
      break;
    }

    case FORWARD_RECURRENCE_TIME : {
      phisto = forward;
      break;
    }

    case NB_EVENT_MIXTURE : {
      phisto = mixture;
      break;
    }
    }

    if ((status) && (phisto->nb_element == 0)) {
      status = false;
      error.update(STAT_error[STATR_EMPTY_SAMPLE]);
    }

    if (status) {
      pdist = NULL;
      pparam = NULL;

      if (renewal) {
        switch (histo_type) {
        case INTER_EVENT :
          pparam = renewal->inter_event;
          break;
        case WITHIN_OBSERVATION_PERIOD :
          pparam = renewal->inter_event;
          break;
        case LENGTH_BIAS :
          pdist = renewal->length_bias;
          break;
        case BACKWARD_RECURRENCE_TIME :
          pdist = renewal->backward;
          break;
        case FORWARD_RECURRENCE_TIME :
          pdist = renewal->forward;
          break;
        case NB_EVENT_MIXTURE :
          pdist = renewal->mixture;
          break;
        }
      }

      if (pdist) {
        histo = new DiscreteDistributionData(*phisto , pdist);
      }
      else {
        histo = new DiscreteDistributionData(*phisto , pparam);
      }
    }
  }

  return histo;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a RenewalData object.
 *
 *  \param[in,out] os         stream,
 *  \param[in]     exhaustive flag detail level,
 *  \param[in]     file_flag  flag file.
 */
/*--------------------------------------------------------------*/

ostream& RenewalData::ascii_write(ostream &os , bool exhaustive , bool file_flag) const

{
  int i , j;
  int nb_value , max , width[2];
  ios_base::fmtflags format_flags;


  format_flags = os.setf(ios::right , ios::adjustfield);

  // writing of the inter-event frequency distribution,
  // the frequency distribution of time intervals between events within the observation period,
  // the length-biased frequency distribution,
  // the backward and forward recurrence time frequency distributions

  if (inter_event) {
    os << SEQ_label[SEQL_INTER_EVENT] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";
    inter_event->ascii_characteristic_print(os , false , file_flag);
    if (file_flag) {
      os << "# ";
    }
    os << STAT_label[STATL_VARIATION_COEFF] << ": "
       << sqrt(inter_event->variance) / inter_event->mean << endl;
  }

  if ((exhaustive) || (!inter_event)) {
    os << "\n";
    if (file_flag) {
      os << "# ";
    }
    os << STAT_label[STATL_OBSERVATION_INTER_EVENT] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";
    within->ascii_characteristic_print(os , false , file_flag);
  }

  if (exhaustive) {
    if (length_bias) {
      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << SEQ_label[SEQL_LENGTH_BIASED] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";
      length_bias->ascii_characteristic_print(os , false , file_flag);
    }

    os << "\n";
    if (file_flag) {
      os << "# ";
    }
    os << STAT_label[STATL_BACKWARD] << " " << SEQ_label[SEQL_RECURRENCE_TIME] << " "
       << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";
    backward->ascii_characteristic_print(os , false , file_flag);

    os << "\n";
    if (file_flag) {
      os << "# ";
    }
    os << STAT_label[STATL_FORWARD] << " " << SEQ_label[SEQL_RECURRENCE_TIME] << " "
       << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";
    forward->ascii_characteristic_print(os , false , file_flag);

    nb_value = within->nb_value;
    if ((length_bias) && (length_bias->nb_value > nb_value)) {
      nb_value = forward->nb_value;
    }
    if (backward->nb_value > nb_value) {
      nb_value = backward->nb_value;
    }
    if (forward->nb_value > nb_value) {
      nb_value = forward->nb_value;
    }

    width[0] = column_width(nb_value - 1);

    if (inter_event) {
      max = inter_event->max;
    }
    else {
      max = within->max;
    }
    if (backward->max > max) {
      max = backward->max;
    }
    if (forward->max > max) {
      max = forward->max;
    }
    width[1] = column_width(max) + ASCII_SPACE;

    os << "\n";
    if (file_flag) {
      os << "# ";
    }
    os << "   ";
    if (inter_event) {
      os << " | " << SEQ_label[SEQL_INTER_EVENT] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
    }
    os << " | " << STAT_label[STATL_OBSERVATION_INTER_EVENT] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
    if (length_bias) {
      os << " | " << SEQ_label[SEQL_LENGTH_BIASED] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
    }
    os << " | " << STAT_label[STATL_BACKWARD] << " " << SEQ_label[SEQL_RECURRENCE_TIME]
       << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " | " << STAT_label[STATL_FORWARD]
       << " " << SEQ_label[SEQL_RECURRENCE_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];

    for (i = 0;i < nb_value;i++) {
      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << setw(width[0]) << i;

      if (inter_event) {
        if (i < inter_event->nb_value) {
          os << setw(width[1]) << inter_event->frequency[i];
        }
        else {
          os << setw(width[1]) << " ";
        }
      }

      if (i < within->nb_value) {
        os << setw(width[1]) << within->frequency[i];
      }
      else {
        os << setw(width[1]) << " ";
      }

      if (length_bias) {
        if (i < length_bias->nb_value) {
          os << setw(width[1]) << length_bias->frequency[i];
        }
        else {
          os << setw(width[1]) << " ";
        }
      }

      if (i < backward->nb_value) {
        os << setw(width[1]) << backward->frequency[i];
      }
      else {
        os << setw(width[1]) << " ";
      }

      if (i < forward->nb_value) {
        os << setw(width[1]) << forward->frequency[i];
      }
      else {
        os << setw(width[1]) << " ";
      }
    }
    os << endl;
  }

  os << "\n";
  if (file_flag) {
    ascii_file_write(os , exhaustive , type);
  }
  else {
    TimeEvents::ascii_write(os , exhaustive , type);
  }

  // writing of no-event/event probabilities as a function of time

  if (exhaustive) {
    os << "\n";
    if (file_flag) {
      os << "# ";
    }
    os << "   | " << SEQ_label[SEQL_OBSERVED] << " " << SEQ_label[SEQL_NO_EVENT_PROBABILITY]
       << " | " << SEQ_label[SEQL_OBSERVED] << " " << SEQ_label[SEQL_EVENT_PROBABILITY]
       << " | " << STAT_label[STATL_FREQUENCY] << endl;

    index_event->ascii_print(os , file_flag);

    // writing of the sequences of events

    os << "\n";
    if (file_flag) {
      os << "# ";
    }
    os << SEQ_label[SEQL_SEQUENCES] << endl;

    for (i = 0;i < nb_element;i++) {
      os << "\n";
      if (file_flag) {
        os << "# ";
      }

      for (j = 0;j < length[i];j++) {
        if ((j > 0) && ((2 * j) % LINE_NB_CHARACTER == 0)) {
          os << "\\" << endl;
          if (file_flag) {
            os << "# ";
          }
        }

        os << sequence[i][j] << " ";
      }

      os << endl;
    }
  }

  os.setf(format_flags , ios::adjustfield);

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a RenewalData object.
 *
 *  \param[in,out] os         stream,
 *  \param[in]     exhaustive flag detail level.
 */
/*--------------------------------------------------------------*/

ostream& RenewalData::ascii_write(ostream &os , bool exhaustive) const

{
  if (renewal) {
    renewal->ascii_write(os , this , exhaustive , false);
  }
  else {
    ascii_write(os , exhaustive , false);
  }

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a RenewalData object in a file.
 *
 *  \param[in] error      reference on a StatError object,
 *  \param[in] path       file path,
 *  \param[in] exhaustive flag detail level.
 *
 *  \return               error status.
 */
/*--------------------------------------------------------------*/

bool RenewalData::ascii_write(StatError &error , const string path ,
                              bool exhaustive) const

{
  bool status = false;


  ofstream out_file(path.c_str());

  error.init();

  if (!out_file) {
    status = false;
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {
    status = true;

    if (renewal) {
      renewal->ascii_write(out_file , this , exhaustive , true);
    }
    else {
      ascii_write(out_file , exhaustive , true);
    }
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a RenewalData object at the spreadsheet format.
 *
 *  \param[in,out] os stream.
 */
/*--------------------------------------------------------------*/

ostream& RenewalData::spreadsheet_write(ostream &os) const

{
  int i;
  int nb_value;


  // writing of the inter-event frequency distribution,
  // the frequency distribution of time intervals between events within the observation period,
  // the length-biased frequency distribution,
  // the backward and forward recurrence time frequency distributions

  if (inter_event) {
    os << "\n" << SEQ_label[SEQL_INTER_EVENT] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\t";
    inter_event->spreadsheet_characteristic_print(os);
    os << STAT_label[STATL_VARIATION_COEFF] << "\t"
       << sqrt(inter_event->variance) / inter_event->mean << endl;
  }

  os << "\n" << STAT_label[STATL_OBSERVATION_INTER_EVENT] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\t";
  within->spreadsheet_characteristic_print(os);

  if (length_bias) {
    os << "\n" << SEQ_label[SEQL_LENGTH_BIASED] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\t";
    length_bias->spreadsheet_characteristic_print(os);
  }
  os << "\n" << STAT_label[STATL_BACKWARD] << " " << SEQ_label[SEQL_RECURRENCE_TIME] << " "
     << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\t";
  backward->spreadsheet_characteristic_print(os);

  os << "\n" << STAT_label[STATL_FORWARD] << " " << SEQ_label[SEQL_RECURRENCE_TIME] << " "
     << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\t";
  forward->spreadsheet_characteristic_print(os);

  nb_value = within->nb_value;
  if ((length_bias) && (length_bias->nb_value > nb_value)) {
    nb_value = forward->nb_value;
  }
  if (backward->nb_value > nb_value) {
    nb_value = backward->nb_value;
  }
  if (forward->nb_value > nb_value) {
    nb_value = forward->nb_value;
  }

  os << "\n";
  if (inter_event) {
    os << "\t" << SEQ_label[SEQL_INTER_EVENT] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
  }
  os << "\t" << STAT_label[STATL_OBSERVATION_INTER_EVENT] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
  if (length_bias) {
    os << " \t" << SEQ_label[SEQL_LENGTH_BIASED] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
  }
  os << "\t" << STAT_label[STATL_BACKWARD] << " " << SEQ_label[SEQL_RECURRENCE_TIME]
     << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\t" << STAT_label[STATL_FORWARD]
     << " " << SEQ_label[SEQL_RECURRENCE_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];

  for (i = 0;i < nb_value;i++) {
    os << "\n" << i;

    if (inter_event) {
      os << "\t";
      if (i < inter_event->nb_value) {
        os << inter_event->frequency[i];
      }
    }

    os << "\t";
    if (i < within->nb_value) {
      os << within->frequency[i];
    }

    if (length_bias) {
      os << "\t";
      if (i < length_bias->nb_value) {
        os << length_bias->frequency[i];
      }
    }

    os << "\t";
    if (i < backward->nb_value) {
      os << backward->frequency[i];
    }
    os << "\t";
    if (i < forward->nb_value) {
      os << forward->frequency[i];
    }
  }
  os << endl;

  os << "\n";
  TimeEvents::spreadsheet_write(os);

  // writing of no-event/event probabilities as a function of time

  os << "\n\t" << SEQ_label[SEQL_OBSERVED] << " " << SEQ_label[SEQL_NO_EVENT_PROBABILITY]
     << "\t" << SEQ_label[SEQL_OBSERVED] << " " << SEQ_label[SEQL_EVENT_PROBABILITY]
     << "\t" << STAT_label[STATL_FREQUENCY] << endl;

  index_event->spreadsheet_print(os);

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a RenewalData object in a file at the spreadsheet format.
 *
 *  \param[in] error reference on a StatError object,
 *  \param[in] path  file path.
 *
 *  \return          error status.
 */
/*--------------------------------------------------------------*/

bool RenewalData::spreadsheet_write(StatError &error , const string path) const

{
  bool status = false;


  ofstream out_file(path.c_str());

  error.init();

  if (!out_file) {
    status = false;
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {
    status = true;

    if (renewal) {
      renewal->spreadsheet_write(out_file , this);
    }
    else {
      spreadsheet_write(out_file);
    }
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Plot of a RenewalData object using Gnuplot.
 *
 *  \param[in] error  reference on a StatError object,
 *  \param[in] prefix file prefix,
 *  \param[in] title  figure title.
 *
 *  \return           error status.
 */
/*--------------------------------------------------------------*/

bool RenewalData::plot_write(StatError &error , const char *prefix ,
                             const char *title) const

{
  bool status = false;


  error.init();

  if (renewal) {
    status = renewal->plot_write(prefix , title , this);
  }

  else {
    int i , j , k;
    int nb_histo;
    const FrequencyDistribution **phisto;
    ostringstream data_file_name[2];


    // writing of the data files

    data_file_name[0] << prefix << 0 << ".dat";

    nb_histo = 2;
    if (inter_event) {
      nb_histo++;
    }
    if (within->nb_element > 0) {
      nb_histo++;
    }
    if (length_bias) {
      nb_histo++;
    }

    for (i = htime->offset;i < htime->nb_value;i++) {
      if (htime->frequency[i] > 0) {
        nb_histo++;
      }
    }
    if (htime->variance > 0.) {
      nb_histo += 2;
    }

    phisto = new const FrequencyDistribution*[nb_histo];

    nb_histo = 0;
    if (inter_event) {
      phisto[nb_histo++] = inter_event;
    }
    if (within->nb_element > 0) {
      phisto[nb_histo++] = within;
    }
    if (length_bias) {
      phisto[nb_histo++] = length_bias;
    }
    phisto[nb_histo++] = backward;
    phisto[nb_histo++] = forward;

    if (htime->variance > 0.) {
      phisto[nb_histo++] = htime;
    }
    for (i = htime->offset;i < htime->nb_value;i++) {
      if (htime->frequency[i] > 0) {
        phisto[nb_histo++] = hnb_event[i];
      }
    }
    if (htime->variance > 0.) {
      phisto[nb_histo++] = mixture;
    }

    status = phisto[0]->plot_print((data_file_name[0].str()).c_str() , nb_histo - 1 , phisto + 1);

    delete [] phisto;

    if (status) {
      data_file_name[1] << prefix << 1 << ".dat";
      index_event->plot_print((data_file_name[1].str()).c_str());

      // writing of the script files

      for (i = 0;i < 2;i++) {
        j = 1;

        ostringstream file_name[2];

        switch (i) {
        case 0 :
          file_name[0] << prefix << ".plot";
          break;
        case 1 :
          file_name[0] << prefix << ".print";
          break;
        }

        ofstream out_file((file_name[0].str()).c_str());

        if (i == 1) {
          out_file << "set terminal postscript" << endl;
          file_name[1] << label(prefix) << ".ps";
          out_file << "set output \"" << file_name[1].str() << "\"\n\n";
        }

        out_file << "set border 15 lw 0\n" << "set tics out\n" << "set xtics nomirror\n"
                 << "set title";
        if (title) {
          out_file << " \"" << title << "\"";
        }
        out_file << "\n\n";

        if (inter_event) {
          if (inter_event->nb_value - 1 < TIC_THRESHOLD) {
            out_file << "set xtics 0,1" << endl;
          }
          if ((int)(inter_event->max * YSCALE) + 1 < TIC_THRESHOLD) {
            out_file << "set ytics 0,1" << endl;
          }

          out_file << "plot [0:" << inter_event->nb_value - 1 << "] [0:"
                   << (int)(inter_event->max * YSCALE) + 1 << "] \""
                   << label((data_file_name[0].str()).c_str()) << "\" using " << j++
                   << " title \"" << SEQ_label[SEQL_INTER_EVENT] << " "
                   << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\" with impulses" << endl;

          if (inter_event->nb_value - 1 < TIC_THRESHOLD) {
            out_file << "set xtics autofreq" << endl;
          }
          if ((int)(inter_event->max * YSCALE) + 1 < TIC_THRESHOLD) {
            out_file << "set ytics autofreq" << endl;
          }

          if (i == 0) {
            out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
          }
          out_file << endl;
        }

        if (within->nb_element > 0) {
          if (within->nb_value - 1 < TIC_THRESHOLD) {
            out_file << "set xtics 0,1" << endl;
          }
          if ((int)(within->max * YSCALE) + 1 < TIC_THRESHOLD) {
            out_file << "set ytics 0,1" << endl;
          }

          out_file << "plot [0:" << within->nb_value - 1 << "] [0:"
                   << (int)(within->max * YSCALE) + 1 << "] \""
                   << label((data_file_name[0].str()).c_str()) << "\" using " << j++
                   << " title \"" << STAT_label[STATL_OBSERVATION_INTER_EVENT] << " "
                   << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\" with impulses" << endl;

          if (within->nb_value - 1 < TIC_THRESHOLD) {
            out_file << "set xtics autofreq" << endl;
          }
          if ((int)(within->max * YSCALE) + 1 < TIC_THRESHOLD) {
            out_file << "set ytics autofreq" << endl;
          }

          if (i == 0) {
            out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
          }
          out_file << endl;
        }

        if (length_bias) {
          if (length_bias->nb_value - 1 < TIC_THRESHOLD) {
            out_file << "set xtics 0,1" << endl;
          }
          if ((int)(length_bias->max * YSCALE) + 1 < TIC_THRESHOLD) {
            out_file << "set ytics 0,1" << endl;
          }

          out_file << "plot [0:" << length_bias->nb_value - 1 << "] [0:"
                   << (int)(length_bias->max * YSCALE) + 1 << "] \""
                   << label((data_file_name[0].str()).c_str()) << "\" using " << j++
                   << " title \"" << SEQ_label[SEQL_LENGTH_BIASED] << " "
                   << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\" with impulses" << endl;

          if (length_bias->nb_value - 1 < TIC_THRESHOLD) {
            out_file << "set xtics autofreq" << endl;
          }
          if ((int)(length_bias->max * YSCALE) + 1 < TIC_THRESHOLD) {
            out_file << "set ytics autofreq" << endl;
          }

          if (i == 0) {
            out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
          }
          out_file << endl;
        }

        if (backward->nb_value - 1 < TIC_THRESHOLD) {
          out_file << "set xtics 0,1" << endl;
        }
        if ((int)(backward->max * YSCALE) + 1 < TIC_THRESHOLD) {
          out_file << "set ytics 0,1" << endl;
        }

        out_file << "plot [0:" << backward->nb_value - 1 << "] [0:"
                 << (int)(backward->max * YSCALE) + 1 << "] \""
                 << label((data_file_name[0].str()).c_str()) << "\" using " << j++
                 << " title \"" << STAT_label[STATL_BACKWARD] << " "
                 << SEQ_label[SEQL_RECURRENCE_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
                 << "\" with impulses" << endl;

        if (backward->nb_value - 1 < TIC_THRESHOLD) {
          out_file << "set xtics autofreq" << endl;
        }
        if ((int)(backward->max * YSCALE) + 1 < TIC_THRESHOLD) {
          out_file << "set ytics autofreq" << endl;
        }

        if (i == 0) {
          out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
        }
        out_file << endl;

        if (forward->nb_value - 1 < TIC_THRESHOLD) {
          out_file << "set xtics 0,1" << endl;
        }
        if ((int)(forward->max * YSCALE) + 1 < TIC_THRESHOLD) {
          out_file << "set ytics 0,1" << endl;
        }

        out_file << "plot [0:" << forward->nb_value - 1 << "] [0:"
                 << (int)(forward->max * YSCALE) + 1 << "] \""
                 << label((data_file_name[0].str()).c_str()) << "\" using " << j++
                 << " title \"" << STAT_label[STATL_FORWARD] << " "
                 << SEQ_label[SEQL_RECURRENCE_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
                 << "\" with impulses" << endl;

        if (forward->nb_value - 1 < TIC_THRESHOLD) {
          out_file << "set xtics autofreq" << endl;
        }
        if ((int)(forward->max * YSCALE) + 1 < TIC_THRESHOLD) {
          out_file << "set ytics autofreq" << endl;
        }

        if (i == 0) {
          out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
        }
        out_file << endl;

        if (htime->variance > 0.) {
          if (htime->nb_value - 1 < TIC_THRESHOLD) {
            out_file << "set xtics 0,1" << endl;
          }
          if ((int)(htime->max * YSCALE) + 1 < TIC_THRESHOLD) {
            out_file << "set ytics 0,1" << endl;
          }

          out_file << "plot [0:" << htime->nb_value - 1 << "] [0:"
                   << (int)(htime->max * YSCALE) + 1 << "] \""
                   << label((data_file_name[0].str()).c_str()) << "\" using " << j++
                   << " title \"" << SEQ_label[SEQL_OBSERVATION_TIME] << " "
                   << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\" with impulses" << endl;

          if (htime->nb_value - 1 < TIC_THRESHOLD) {
            out_file << "set xtics autofreq" << endl;
          }
          if ((int)(htime->max * YSCALE) + 1 < TIC_THRESHOLD) {
            out_file << "set ytics autofreq" << endl;
          }
        }

        for (k = htime->offset;k < htime->nb_value;k++) {
          if (htime->frequency[k] > 0) {
            if (((htime->variance > 0.) || (k > htime->offset)) && (i == 0)) {
              out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
            }
            out_file << endl;

            if (hnb_event[k]->nb_value - 1 < TIC_THRESHOLD) {
              out_file << "set xtics 0,1" << endl;
            }
            if ((int)(hnb_event[k]->max * YSCALE) + 1 < TIC_THRESHOLD) {
              out_file << "set ytics 0,1" << endl;
            }

            out_file << "plot [0:" << hnb_event[k]->nb_value - 1 << "] [0:"
                     << (int)(hnb_event[k]->max * YSCALE) + 1 << "] \""
                     << label((data_file_name[0].str()).c_str()) << "\" using " << j++
                     << " title \"" << SEQ_label[SEQL_NB_EVENT] << " "
                     << SEQ_label[SEQL_DURING] << " " << k << " " << SEQ_label[SEQL_TIME_UNIT] << " "
                     << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\" with impulses" << endl;

            if (hnb_event[k]->nb_value - 1 < TIC_THRESHOLD) {
              out_file << "set xtics autofreq" << endl;
            }
            if ((int)(hnb_event[k]->max * YSCALE) + 1 < TIC_THRESHOLD) {
              out_file << "set ytics autofreq" << endl;
            }
          }
        }

        if (htime->variance > 0.) {
          if (i == 0) {
            out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
          }
          out_file << endl;

          if (mixture->nb_value - 1 < TIC_THRESHOLD) {
            out_file << "set xtics 0,1" << endl;
          }
          if ((int)(mixture->max * YSCALE) + 1 < TIC_THRESHOLD) {
            out_file << "set ytics 0,1" << endl;
          }

          out_file << "plot [0:" << mixture->nb_value - 1 << "] [0:"
                   << (int)(mixture->max * YSCALE) + 1 << "] \""
                   << label((data_file_name[0].str()).c_str()) << "\" using " << j
                   << " title \"" << SEQ_label[SEQL_NB_EVENT] << " "
                   << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\" with impulses" << endl;

          if (mixture->nb_value - 1 < TIC_THRESHOLD) {
            out_file << "set xtics autofreq" << endl;
          }
          if ((int)(mixture->max * YSCALE) + 1 < TIC_THRESHOLD) {
            out_file << "set ytics autofreq" << endl;
          }
        }

        if (i == 0) {
          out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
        }
        out_file << endl;

        if (index_event->length - 1 < TIC_THRESHOLD) {
          out_file << "set xtics 0,1" << endl;
        }

        out_file << "plot [" << index_event->offset << ":" << index_event->length - 1 << "] [0:1] "
                 << "\"" << label((data_file_name[1].str()).c_str()) << "\" using "
                 << 1 << " title \"" << SEQ_label[SEQL_OBSERVED] << " "
                 << SEQ_label[SEQL_NO_EVENT_PROBABILITY] << " \" with linespoints,\\" << endl;
        out_file << "\"" << label((data_file_name[1].str()).c_str()) << "\" using "
                 << 2 << " title \"" << SEQ_label[SEQL_OBSERVED] << " "
                 << SEQ_label[SEQL_EVENT_PROBABILITY] << " \" with linespoints" << endl;

        if (index_event->length - 1 < TIC_THRESHOLD) {
          out_file << "set xtics autofreq" << endl;
        }

        if (i == 1) {
          out_file << "\nset terminal x11" << endl;
        }

        out_file << "\npause 0 \"" << STAT_label[STATL_END] << "\"" << endl;
      }
    }
  }

  if (!status) {
    error.update(STAT_error[STATR_FILE_PREFIX]);
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Plot of a RenewalData object.
 *
 *  \return MultiPlotSet object.
 */
/*--------------------------------------------------------------*/

MultiPlotSet* RenewalData::get_plotable() const

{
  MultiPlotSet *plot_set;


  if (renewal) {
    plot_set = renewal->get_plotable(this);
  }

  else {
    int i , j , k , m;
    int nb_plot_set , nb_histo , max_nb_value , max_frequency;
    double shift;
    const FrequencyDistribution *phisto[2] , **merged_histo;
    ostringstream title , legend;


    nb_plot_set = 3;
    if (inter_event) {
      nb_plot_set++;
    }
    if ((within->nb_element > 0) || (length_bias)) {
      nb_plot_set++;
    }
    if (htime->variance > 0.) {
      nb_plot_set += 2;
    }

    plot_set = new MultiPlotSet(nb_plot_set);
    MultiPlotSet &plot = *plot_set;

    plot.border = "15 lw 0";

    i = 0;
    if (inter_event) {

      // inter-event frequency distribution

      plot[i].xrange = Range(0 , inter_event->nb_value - 1);
      plot[i].yrange = Range(0 , ceil(inter_event->max * YSCALE));

      if (inter_event->nb_value - 1 < TIC_THRESHOLD) {
        plot[i].xtics = 1;
      }
      if (ceil(inter_event->max * YSCALE) < TIC_THRESHOLD) {
        plot[i].ytics = 1;
      }

      plot[i].resize(within->nb_element > 0 ? 2 : 1);

      legend.str("");
      legend << SEQ_label[SEQL_LENGTH_BIASED] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
      plot[i][0].legend = legend.str();

      plot[i][0].style = "impulses";

      inter_event->plotable_frequency_write(plot[i][0]);

      if (within->nb_element > 0) {
        legend.str("");
        legend << STAT_label[STATL_OBSERVATION_INTER_EVENT] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
        plot[i][1].legend = legend.str();

        plot[i][1].style = "impulses";

        within->plotable_frequency_write(plot[i][1]);
      }

      i++;
    }

    if ((within->nb_element > 0) || (length_bias)) {

      // frequency distribution of time intervals between events within the observation period and
      // length-biased frequency distribution

      nb_histo = 0;
      max_nb_value = 0;
      max_frequency = 0;

      if (within->nb_element > 0) {
        nb_histo++;

        if (within->nb_value > max_nb_value) {
          max_nb_value = within->nb_value;
        }
        if (within->max > max_frequency) {
          max_frequency = within->max;
        }
      }

      if (length_bias) {
        nb_histo++;

        if (length_bias->nb_value > max_nb_value) {
          max_nb_value = length_bias->nb_value;
        }
        if (length_bias->max > max_frequency) {
          max_frequency = length_bias->max;
        }
      }

      plot[i].xrange = Range(0 , max_nb_value);
      plot[i].yrange = Range(0 , ceil(max_frequency * YSCALE));

      if (max_nb_value < TIC_THRESHOLD) {
        plot[i].xtics = 1;
      }
      if (ceil(max_frequency * YSCALE) < TIC_THRESHOLD) {
        plot[i].ytics = 1;
      }

      plot[i].resize(nb_histo);

      j = 0;
      if (within->nb_element > 0) {
        legend.str("");
        legend << STAT_label[STATL_OBSERVATION_INTER_EVENT] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
        plot[i][j].legend = legend.str();

        plot[i][j].style = "impulses";

        within->plotable_frequency_write(plot[i][j]);
        j++;
      }

      if (length_bias) {
        legend.str("");
        legend << SEQ_label[SEQL_LENGTH_BIASED] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
        plot[i][j].legend = legend.str();

        plot[i][j].style = "impulses";

        for (k = length_bias->offset;k < length_bias->nb_value;k++) {
          if (length_bias->frequency[k] > 0) {
            plot[i][j].add_point(k + j * PLOT_SHIFT , length_bias->frequency[k]);
          }
        }
      }

      i++;
    }

    // backward and forward recurrence time frequency distributions

    max_nb_value = MAX(backward->nb_value , forward->nb_value);
    max_frequency = MAX(backward->max , forward->max);

    plot[i].xrange = Range(0 , max_nb_value);
    plot[i].yrange = Range(0 , ceil(max_frequency * YSCALE));

    if (max_nb_value < TIC_THRESHOLD) {
      plot[i].xtics = 1;
    }
    if (ceil(max_frequency * YSCALE) < TIC_THRESHOLD) {
      plot[i].ytics = 1;
    }

    plot[i].resize(2);

    legend.str("");
    legend << STAT_label[STATL_BACKWARD] << " " << SEQ_label[SEQL_RECURRENCE_TIME] << " "
           << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
    plot[i][0].legend = legend.str();

    plot[i][0].style = "impulses";

    backward->plotable_frequency_write(plot[i][0]);

    legend.str("");
    legend << STAT_label[STATL_FORWARD] << " " << SEQ_label[SEQL_RECURRENCE_TIME] << " "
           << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
    plot[i][1].legend = legend.str();

    plot[i][1].style = "impulses";

    for (j = forward->offset;j < forward->nb_value;j++) {
      if (forward->frequency[j] > 0) {
        plot[i][1].add_point(j + PLOT_SHIFT , forward->frequency[j]);
      }
    }
    i++;

    if (htime->variance > 0.) {

      // observation period frequency distribution

      plot[i].xrange = Range(0 , htime->nb_value - 1);
      plot[i].yrange = Range(0 , ceil(htime->max * YSCALE));

      if (htime->nb_value - 1 < TIC_THRESHOLD) {
        plot[i].xtics = 1;
      }
      if (ceil(htime->max * YSCALE) < TIC_THRESHOLD) {
        plot[i].ytics = 1;
      }

      plot[i].resize(1);

      legend.str("");
      legend << SEQ_label[SEQL_OBSERVATION_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
      plot[i][0].legend = legend.str();

      plot[i][0].style = "impulses";

      htime->plotable_frequency_write(plot[i][0]);
      i++;
    }

    // number of events frequency distribution for each observation period

    if (htime->variance > 0.) {
      title.str("");
      title << SEQ_label[SEQL_NB_EVENT] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTIONS];
      plot[i].title = title.str();
    }

    // computation of the maximum number of values and the maximum frequency

    nb_histo = 0;
    max_nb_value = 0;
    max_frequency = 0;

    for (j = htime->offset;j < htime->nb_value;j++) {
      if (htime->frequency[j] > 0) {
        nb_histo++;

        if (hnb_event[j]->nb_value > max_nb_value) {
          max_nb_value = hnb_event[j]->nb_value;
        }
        if (hnb_event[j]->max > max_frequency) {
          max_frequency = hnb_event[j]->max;
        }
      }
    }

    plot[i].xrange = Range(0 , max_nb_value);
    plot[i].yrange = Range(0 , ceil(max_frequency * YSCALE));

    if (max_nb_value < TIC_THRESHOLD) {
      plot[i].xtics = 1;
    }
    if (ceil(max_frequency * YSCALE) < TIC_THRESHOLD) {
      plot[i].ytics = 1;
    }

    plot[i].resize(nb_histo);

    j = 0;
    shift = 0.;

    for (k = htime->offset;k < htime->nb_value;k++) {
      if (htime->frequency[k] > 0) {
        legend.str("");
        if (htime->variance > 0.) {
          legend << k << " " << SEQ_label[SEQL_TIME_UNIT];
        }
        else {
          legend << SEQ_label[SEQL_NB_EVENT] << " " << SEQ_label[SEQL_DURING] << " " << k << " "
                 << SEQ_label[SEQL_TIME_UNIT] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
        }
        plot[i][j].legend = legend.str();

        plot[i][j].style = "impulses";

        for (m = hnb_event[k]->offset;m < hnb_event[k]->nb_value;m++) {
          if (hnb_event[k]->frequency[m] > 0) {
            plot[i][j].add_point(m + shift , hnb_event[k]->frequency[m]);
          }
        }

        if (PLOT_SHIFT * (nb_histo - 1) < PLOT_MAX_SHIFT) {
          shift += PLOT_SHIFT;
        }
        else {
          shift += PLOT_MAX_SHIFT / (nb_histo - 1);
        }

        j++;
      }
    }
    i++;

    if (htime->variance > 0.) {

      // superimposed number of events frequency distributions

      merged_histo = new const FrequencyDistribution*[nb_histo];

      j = nb_histo - 1;
      for (k = htime->nb_value - 1;k >= htime->offset;k--) {
        if (htime->frequency[k] > 0) {
          if (j == nb_histo - 1) {
            merged_histo[j] = new FrequencyDistribution(*hnb_event[k]);
          }

          else {
            phisto[0] = merged_histo[j + 1];
            phisto[1] = hnb_event[k];
            merged_histo[j] = new FrequencyDistribution(2 , phisto);
          }

          j--;
        }
      }

      title.str("");
      title << SEQ_label[SEQL_NB_EVENT] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTIONS];
      plot[i].title = title.str();

      plot[i].xrange = Range(0 , merged_histo[0]->nb_value - 1);
      plot[i].yrange = Range(0 , ceil(merged_histo[0]->max * YSCALE));

      if (merged_histo[0]->nb_value - 1 < TIC_THRESHOLD) {
        plot[i].xtics = 1;
      }
      if (ceil(merged_histo[0]->max * YSCALE) < TIC_THRESHOLD) {
        plot[i].ytics = 1;
      }

      plot[i].resize(nb_histo);

      j = 0;
      for (k = htime->offset;k < htime->nb_value;k++) {
        if (htime->frequency[k] > 0) {
          legend.str("");
          legend << k << " " << SEQ_label[SEQL_TIME_UNIT];
          plot[i][j].legend = legend.str();

          plot[i][j].style = "impulses";

          merged_histo[j]->plotable_frequency_write(plot[i][j]);
          j++;
        }
      }

      for (j = 0;j < nb_histo;j++) {
        delete merged_histo[j];
      }
      delete [] merged_histo;

/*      plot[i].resize(1);

      legend.str("");
      legend << SEQ_label[SEQL_NB_EVENT] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
      plot[i][0].legend = legend.str();

      plot[i][0].style = "impulses";

      mixture->plotable_frequency_write(plot[i][0]); */

      i++;
    }

    // no-event/event probabilities as a function of time

    plot[i].xrange = Range(0 , index_event->length - 1);
    plot[i].yrange = Range(0. , 1.);

    if (index_event->length - 1 < TIC_THRESHOLD) {
      plot[i].xtics = 1;
    }

    plot[i].resize(2);

    legend.str("");
    legend << SEQ_label[SEQL_OBSERVED] << " " << SEQ_label[SEQL_NO_EVENT_PROBABILITY];
    plot[i][0].legend = legend.str();

    plot[i][0].style = "linespoints";

    index_event->plotable_write(0 , plot[i][0]);

    legend.str("");
    legend << SEQ_label[SEQL_OBSERVED] << " " << SEQ_label[SEQL_EVENT_PROBABILITY];
    plot[i][1].legend = legend.str();

    plot[i][1].style = "linespoints";

    index_event->plotable_write(1 , plot[i][1]);
  }

  return plot_set;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Extraction of no-event/event probabilities as a function of time.
 *
 *  \param[in] offset shift (0/1).
 */
/*--------------------------------------------------------------*/

void RenewalData::build_index_event(int offset)

{
  int i , j;
  int frequency[2];


  index_event = new Curves(2 , htime->nb_value , true , false , false);
  index_event->offset = offset;

  for (i = index_event->offset;i < index_event->length;i++) {
    frequency[0] = 0;
    frequency[1] = 0;

    for (j = 0;j < nb_element;j++) {
      if (i - index_event->offset < length[j]) {
        frequency[sequence[j][i - index_event->offset]]++;
      }
    }

    index_event->frequency[i] = frequency[0] + frequency[1];
    index_event->point[0][i] = (double)frequency[0] / (double)index_event->frequency[i];
    index_event->point[1][i] = (double)frequency[1] / (double)index_event->frequency[i];
  }
}


};  // namespace sequence_analysis
