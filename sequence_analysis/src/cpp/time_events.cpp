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
 *       $Id: time_events.cpp 18079 2015-04-23 10:58:57Z guedon $
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
#include <iomanip>
#include <iostream>

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

#include "renewal.h"
#include "sequence_label.h"
#include "tool/config.h"

using namespace std;
using namespace stat_tool;


namespace sequence_analysis {



/*--------------------------------------------------------------*
 *
 *  Construction des lois empiriques du temps et du nombre d'evenements
 *  a partir des echantillons {temps, nombre d'evenements, effectif}.
 *
 *--------------------------------------------------------------*/

void TimeEvents::build_frequency_distribution()

{
  register int i;
  int max_nb_event , *ptime , *pnb_event , *pfrequency;


  // construction de la loi empirique des temps et
  // creation des lois empiriques du nombre d'evenements

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

  // constrution des lois empiriques du nombre d'evenements

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


/*--------------------------------------------------------------*
 *
 *  Construction des echantillons {temps, nombre d'evenements, effectif}
 *  a partir des lois empiriques du temps et du nombre d'evenements.
 *
 *--------------------------------------------------------------*/

void TimeEvents::build_sample()

{
  register int i , j;
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


/*--------------------------------------------------------------*
 *
 *  Construction d'un objet TimeEvents a partir d'echantillons
 *  {temps, nombre d'evenements}.
 *
 *  argument : nombre d'echantillons {temps, nombre d'evenements},
 *             pointeurs sur les temps et sur les nombres d'evenements.
 *
 *--------------------------------------------------------------*/

void TimeEvents::build(int inb_element , int *itime , int *inb_event)

{
  register int i , j , k;
  int btime , min_time , max_time , bnb_event , min_nb_event , max_nb_event ,
      nb_selected , *ptime , *pnb_event , *selected_nb_event , *snb_event;


  nb_element = inb_element;

  // calcul temps maximum/minimum et nombre d'evenements maximum/minimum

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

  // tri des echantillons {temps, nombre d'evenements}
  // par temps croissant puis par nombre d'evenements croissant

  btime = 0;
  nb_class = 0;
  i = 0;

  do {

    // recherche du temps minimum courant

    ptime = itime;
    min_time = max_time;
    for (j = 0;j < nb_element;j++) {
      if ((*ptime > btime) && (*ptime < min_time)) {
        min_time = *ptime;
      }
      ptime++;
    }
    btime = min_time;

    // extraction des echantillons correspondant au temps courant

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

    // tri des echantillons correspondant au temps courant

    bnb_event = -1;
    j = 0;

    do {

      // recherche du nombre d'evenements minimum courant

      snb_event = selected_nb_event;
      min_nb_event = max_nb_event;
      for (k = 0;k < nb_selected;k++) {
        if ((*snb_event > bnb_event) && (*snb_event < min_nb_event)) {
          min_nb_event = *snb_event;
        }
        snb_event++;
      }
      bnb_event = min_nb_event;

      // constitution de l'echantillon {temps, nombre d'evenements, effectif}

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


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe TimeEvents.
 *
 *  argument : nombre de classes.
 *
 *--------------------------------------------------------------*/

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
    register int i;
    int *ptime , *pnb_event , *pfrequency;

    time = new int[nb_class];
    nb_event = new int[nb_class];
    frequency = new int[nb_class];

    ptime = time;
    pnb_event = nb_event;
    pfrequency = frequency;

    for (i = 0;i < nb_class;i++) {
      *ptime++ = 0;
      *pnb_event++ = 0;
      *pfrequency++ = 0;
    }
  }

  htime = NULL;
  hnb_event = NULL;
  mixture = NULL;
}


/*--------------------------------------------------------------*
 *
 *  Construction d'un objet TimeEvents a partir d'un objet FrequencyDistribution.
 *
 *  arguments : reference sur un objet FrequencyDistribution, temps d'observation.
 *
 *--------------------------------------------------------------*/

TimeEvents::TimeEvents(FrequencyDistribution &inb_event, int itime)

{
  register int i;
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

  // constitution des echantillons

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

  // construction des lois empiriques

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


/*--------------------------------------------------------------*
 *
 *  Copie d'un objet TimeEvents.
 *
 *  argument : reference sur un objet TimeEvents.
 *
 *--------------------------------------------------------------*/

void TimeEvents::copy(const TimeEvents &timev)

{
  register int i;
  int *ptime , *pnb_event , *pfrequency , *ttime , *tnb_event , *tfrequency;


  // copie des echantillons

  nb_element = timev.nb_element;
  nb_class = timev.nb_class;

  time = new int[nb_class];
  nb_event = new int[nb_class];
  frequency = new int[nb_class];

  ptime = time;
  ttime = timev.time;
  pnb_event = nb_event;
  tnb_event = timev.nb_event;
  pfrequency = frequency;
  tfrequency = timev.frequency;

  for (i = 0;i < nb_class;i++) {
    *ptime++ = *ttime++;
    *pnb_event++ = *tnb_event++;
    *pfrequency++ = *tfrequency++;
  }

  // copie des lois empiriques

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


/*--------------------------------------------------------------*
 *
 *  Fusion d'objets TimeEvents.
 *
 *  argument : nombre d'objets TimeEvents,
 *             pointeurs sur les objets TimeEvents.
 *
 *--------------------------------------------------------------*/

void TimeEvents::merge(int nb_sample , const TimeEvents **ptimev)

{
  register int i , j;
  int nb_histo;
  const FrequencyDistribution **phisto;


  nb_element = 0;
  for (i = 0;i < nb_sample;i++) {
    nb_element += ptimev[i]->nb_element;
  }

  phisto = new const FrequencyDistribution*[nb_sample];

  // fusion des lois empiriques du temps d'observation

  for (i = 0;i < nb_sample;i++) {
    phisto[i] = ptimev[i]->htime;
  }
  htime = new FrequencyDistribution(nb_sample , phisto);

  // fusion des lois empiriques du nombre d'evenements

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


/*--------------------------------------------------------------*
 *
 *  Destruction des champs d'un objet TimeEvents.
 *
 *--------------------------------------------------------------*/

void TimeEvents::remove()

{
  register int i;


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


/*--------------------------------------------------------------*
 *
 *  Destructeur de la classe TimeEvents.
 *
 *--------------------------------------------------------------*/

TimeEvents::~TimeEvents()

{
  remove();
}


/*--------------------------------------------------------------*
 *
 *  Operateur d'assignement de la classe TimeEvents.
 *
 *  argument : reference sur un objet TimeEvents.
 *
 *--------------------------------------------------------------*/

TimeEvents& TimeEvents::operator=(const TimeEvents &timev)

{
  if (&timev != this) {
    remove();
    copy(timev);
  }

  return *this;
}


/*--------------------------------------------------------------*
 *
 *  Extraction de la loi empirique du nombre d'evenements pour
 *  un temps d'observation donne.
 *
 *  arguments : reference sur un objet StatError, type de loi empirique,
 *              temps d'observation.
 *
 *--------------------------------------------------------------*/

DiscreteDistributionData* TimeEvents::extract(StatError &error , int histo_type ,
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


/*--------------------------------------------------------------*
 *
 *  Changement de l'unite de temps d'un objet TimeEvents.
 *
 *  arguments : reference sur un objet StatError, facteur d'echelle.
 *
 *--------------------------------------------------------------*/

TimeEvents* TimeEvents::time_scaling(StatError &error , int scaling_coeff) const

{
  bool status = true;
  register int i;
  int *ptime , *pnb_event , *pfrequency , *ttime , *tnb_event , *tfrequency;
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

    ptime = timev->time;
    ttime = time;
    pnb_event = timev->nb_event;
    tnb_event = nb_event;
    pfrequency = timev->frequency;
    tfrequency = frequency;

    for (i = 0;i < nb_class;i++) {
      *ptime++ = *ttime++ * scaling_coeff;
      *pnb_event++ = *tnb_event++;
      *pfrequency++ = *tfrequency++;
    }

    timev->build_frequency_distribution();
  }

  return timev;
}


/*--------------------------------------------------------------*
 *
 *  Selection d'echantillons {temps, nombre d'evenements}
 *  sur un critere de temps d'observation.
 *
 *  arguments : reference sur un objet StatError,
 *              bornes sur le temps d'observation.
 *
 *--------------------------------------------------------------*/

TimeEvents* TimeEvents::time_select(StatError &error , int min_time ,
                                    int max_time) const

{
  bool status = true;
  register int i;
  int bnb_class , *ptime , *pnb_event , *pfrequency , *ttime , *tnb_event , *tfrequency;
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

    // calcul du nombre de classes

    ttime = time;
    bnb_class = 0;

    for (i = 0;i < nb_class;i++) {
      if ((*ttime >= min_time) && (*ttime <= max_time)) {
         bnb_class++;
      }
      ttime++;
    }

    // copie des echantillons selectionnes

    timev = new TimeEvents(bnb_class);

    ptime = timev->time;
    ttime = time;
    pnb_event = timev->nb_event;
    tnb_event = nb_event;
    pfrequency = timev->frequency;
    tfrequency = frequency;

    for (i = 0;i < nb_class;i++) {
      if ((*ttime >= min_time) && (*ttime <= max_time)) {
        *ptime++ = *ttime;
        *pnb_event++ = *tnb_event;
        *pfrequency++ = *tfrequency;
      }
      ttime++;
      tnb_event++;
      tfrequency++;
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


/*--------------------------------------------------------------*
 *
 *  Selection d'echantillons {temps, nombre d'evenements}
 *  sur un critere de nombre d'evenements.
 *
 *  arguments : reference sur un objet StatError,
 *              bornes sur le nombre d'evenements.
 *
 *--------------------------------------------------------------*/

TimeEvents* TimeEvents::nb_event_select(StatError &error , int min_nb_event ,
                                        int max_nb_event) const

{
  bool status = true;
  register int i;
  int bnb_class , *ptime , *pnb_event , *pfrequency , *ttime , *tnb_event , *tfrequency;
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

    // calcul du nombre de classes

    tnb_event = nb_event;
    bnb_class = 0;

    for (i = 0;i < nb_class;i++) {
      if ((*tnb_event >= min_nb_event) && (*tnb_event <= max_nb_event)) {
         bnb_class++;
      }
      tnb_event++;
    }

    // copie des echantillons selectionnes

    timev = new TimeEvents(bnb_class);

    ptime = timev->time;
    ttime = time;
    pnb_event = timev->nb_event;
    tnb_event = nb_event;
    pfrequency = timev->frequency;
    tfrequency = frequency;

    for (i = 0;i < nb_class;i++) {
      if ((*tnb_event >= min_nb_event) && (*tnb_event <= max_nb_event)) {
        *ptime++ = *ttime;
        *pnb_event++ = *tnb_event;
        *pfrequency++ = *tfrequency;
      }
      ttime++;
      tnb_event++;
      tfrequency++;
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


/*--------------------------------------------------------------*
 *
 *  Construction d'un objet TimeEvents a partir d'un objet FrequencyDistribution.
 *
 *  arguments : references sur un objet StatError et sur un objet FrequencyDistribution,
 *              temps d'observation.
 *
 *--------------------------------------------------------------*/

TimeEvents* build_time_events(StatError &error , FrequencyDistribution &nb_event , int itime)

{
  bool status = true;
  TimeEvents *timev;


  timev = NULL;
  error.init();

  if (itime / (nb_event.nb_value - 1) < MIN_INTER_EVENT) {
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


/*--------------------------------------------------------------*
 *
 *  Construction d'un objet TimeEvents a partir d'un fichier.
 *  Format : n lignes de la forme (temps > 0) (nombre d'evenements >= 0)
 *           (effectif >= 0).
 *
 *  arguments : reference sur un objet StatError, path.
 *
 *--------------------------------------------------------------*/

TimeEvents* time_events_ascii_read(StatError &error , const char *path)

{
  RWLocaleSnapshot locale("en");
  RWCString buffer , token;
  size_t position;
  bool status , lstatus;
  register int i , j;
  int line , nb_class , nb_element;
  long value , time , nb_event;
  TimeEvents *timev;
  ifstream in_file(path);


  timev = NULL;
  error.init();

  if (!in_file) {
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {

    // 1ere passe : analyse des lignes du fichier et recherche
    // du nombre d'echantillons {temps, nombre d'evenements, effectif}

    status = true;
    line = 0;
    time = 0;
    nb_event = -1;
    nb_class = 0;
    nb_element = 0;

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
        if (i <= 2) {
          lstatus = locale.stringToNum(token , &value);

          // test entier strictement positif (temps) ou
          // entier positif ou nul (nombre d'evenements, effectif)

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

            // test echantillons ordonnes (temps)

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

            // test echantillons ordonnes (nombre d'evenements)

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

      // test trois valeurs par ligne

      if ((i > 0) && (i != 3)) {
        status = false;
        error.correction_update(STAT_parsing[STATP_NB_TOKEN] , 3 , line);
      }
    }

    if (nb_element == 0) {
      status = false;
      error.update(STAT_parsing[STATP_EMPTY_SAMPLE]);
    }

    // 2eme passe : lecture du fichier

    if (status) {
//      in_file.close();
//      in_file.open(path , ios::in);
      in_file.clear();
      in_file.seekg(0,ios::beg);

      timev = new TimeEvents(nb_class);
      timev->nb_element = nb_element;

      i = 0;

      while (buffer.readLine(in_file , false)) {
        position = buffer.first('#');
        if (position != RW_NPOS) {
          buffer.remove(position);
        }
        j = 0;

        RWCTokenizer next(buffer);

        while (!((token = next()).isNull())) {
          locale.stringToNum(token , &value);

          switch (j) {
          case 0 :
            timev->time[i] = value;
            break;
          case 1 :
            timev->nb_event[i] = value;
            break;
          case 2 :
            timev->frequency[i] = value;
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


/*--------------------------------------------------------------*
 *
 *  Construction d'un objet TimeEvents a partir d'un fichier.
 *  Format : n lignes de la forme (temps > 0) (nombre d'evenements >= 0).
 *
 *  arguments : reference sur un objet StatError, path.
 *
 *--------------------------------------------------------------*/

TimeEvents* old_time_events_ascii_read(StatError &error , const char *path)

{
  RWLocaleSnapshot locale("en");
  RWCString buffer , token;
  size_t position;
  bool status , lstatus;
  register int i , j;
  int line , nb_element , *ptime , *pnb_event;
  long value;
  TimeEvents *timev;
  ifstream in_file(path);


  timev = NULL;
  error.init();

  if (!in_file) {
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {

    // 1ere passe : analyse des lignes du fichier et recherche
    // du nombre d'echantillons {temps, nombre d'evenements}

    status = true;
    line = 0;
    nb_element = 0;

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
        if (i <= 1) {
          lstatus = locale.stringToNum(token , &value);

          // test entier strictement positif (temps) ou
          // entier positif ou nul (nombre d'evenements)

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

      // test deux valeurs par ligne

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

    // 2eme passe : lecture du fichier

    if (status) {
//      in_file.close();
//      in_file.open(path , ios::in);
        in_file.clear();
        in_file.seekg(0,ios::beg);

      ptime = new int[nb_element];
      pnb_event = new int[nb_element];
      i = 0;

      while (buffer.readLine(in_file , false)) {
        position = buffer.first('#');
        if (position != RW_NPOS) {
          buffer.remove(position);
        }
        j = 0;

        RWCTokenizer next(buffer);

        while (!((token = next()).isNull())) {
          locale.stringToNum(token , &value);

          switch (j) {
          case 0 :
            ptime[i] = value;
            break;
          case 1 :
            pnb_event[i] = value;
            break;
          }

          j++;
        }

        i++;
      }

      timev = new TimeEvents(nb_element , ptime , pnb_event);

      delete [] ptime;
      delete [] pnb_event;
    }
  }

  return timev;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture sur une ligne d'un objet TimeEvents.
 *
 *  argument : stream.
 *
 *--------------------------------------------------------------*/

ostream& TimeEvents::line_write(ostream &os) const

{
  os << STAT_label[STATL_SAMPLE_SIZE] << ": " << nb_element << "   "
     << SEQ_label[SEQL_OBSERVATION_TIME] << " " << STAT_label[STATL_MEAN] << ": " << htime->mean << "   "
     << SEQ_label[SEQL_NB_EVENT] << " " << STAT_label[STATL_MEAN] << ": " << mixture->mean;

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet TimeEvents.
 *
 *  arguments : stream, flag niveau de detail, type de processus.
 *
 *--------------------------------------------------------------*/

ostream& TimeEvents::ascii_write(ostream &os , bool exhaustive , char type) const

{
  register int i;


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

      case 'o' : {
        os << "\n" << SEQ_label[SEQL_1_CENSORED_INTER_EVENT] << ": " << hnb_event[i]->nb_element
           << " (" << 1. / (hnb_event[i]->mean + 1.) << ")" << endl;

        os << SEQ_label[SEQL_COMPLETE_INTER_EVENT] << ": " << hnb_event[i]->mean *
              hnb_event[i]->nb_element << " ("
           << hnb_event[i]->mean / (hnb_event[i]->mean + 1.) << ")" << endl;
        break;
      }

      case 'e' : {
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

    case 'o' : {
      os << "\n" << SEQ_label[SEQL_1_CENSORED_INTER_EVENT] << ": " << mixture->nb_element
         << " (" << 1. / (mixture->mean + 1.) << ")" << endl;

      os << SEQ_label[SEQL_COMPLETE_INTER_EVENT] << ": " << mixture->mean *
            mixture->nb_element << " ("
         << mixture->mean / (mixture->mean + 1.) << ")" << endl;
      break;
    }

    case 'e' : {
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


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet TimeEvents.
 *
 *  arguments : stream, flag niveau de detail.
 *
 *--------------------------------------------------------------*/

ostream& TimeEvents::ascii_write(ostream &os , bool exhaustive) const

{
  return ascii_write(os , exhaustive , 'v');
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet TimeEvents dans un fichier.
 *
 *  arguments : reference sur un objet StatError, path,
 *              flag niveau de detail, type de processus.
 *
 *--------------------------------------------------------------*/

ostream& TimeEvents::ascii_file_write(ostream &os , bool exhaustive , char type) const

{
  register int i;
  int max_frequency , *pfrequency , width[3];


  if ((htime->variance > 0.) && (exhaustive)) {
    os << "# " << SEQ_label[SEQL_OBSERVATION_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";
    htime->ascii_characteristic_print(os , false , true);

    os << "\n#    | " << SEQ_label[SEQL_OBSERVATION_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << endl;
    htime->ascii_print(os , true);
  }

  // calcul des largeurs des colonnes

  if (exhaustive) {
    width[0] = column_width(time[nb_class - 1]);
    width[1] = column_width(nb_event[nb_class - 1]) + ASCII_SPACE;

    pfrequency = frequency;
    max_frequency = 0;

    for (i = 0;i < nb_class;i++) {
      if (*pfrequency > max_frequency) {
        max_frequency = *pfrequency;
      }
      pfrequency++;
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

      case 'o' : {
        os << "\n# " << SEQ_label[SEQL_1_CENSORED_INTER_EVENT] << ": " << hnb_event[time[i]]->nb_element
           << " (" << 1. / (hnb_event[time[i]]->mean + 1.) << ")" << endl;

        os << "# " << SEQ_label[SEQL_COMPLETE_INTER_EVENT] << ": " << hnb_event[time[i]]->mean *
              hnb_event[time[i]]->nb_element << " ("
           << hnb_event[time[i]]->mean / (hnb_event[time[i]]->mean + 1.) << ")" << endl;
        break;
      }

      case 'e' : {
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

    // ecriture des echantillons (temps, nombre d'evenements)

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

    case 'o' : {
      os << "\n# " << SEQ_label[SEQL_1_CENSORED_INTER_EVENT] << ": " << mixture->nb_element
         << " (" << 1. / (mixture->mean + 1.) << ")" << endl;

      os << "# " << SEQ_label[SEQL_COMPLETE_INTER_EVENT] << ": " << mixture->mean *
            mixture->nb_element << " ("
         << mixture->mean / (mixture->mean + 1.) << ")" << endl;
      break;
    }

    case 'e' : {
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


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet TimeEvents dans un fichier.
 *
 *  arguments : reference sur un objet StatError, path,
 *              flag niveau de detail.
 *
 *--------------------------------------------------------------*/

bool TimeEvents::ascii_write(StatError &error , const char *path ,
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
    ascii_file_write(out_file , exhaustive);
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet TimeEvents au format tableur.
 *
 *  argument : stream, type de processus.
 *
 *--------------------------------------------------------------*/

ostream& TimeEvents::spreadsheet_write(ostream &os , char type) const

{
  register int i;


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

      case 'o' : {
        os << "\n" << SEQ_label[SEQL_1_CENSORED_INTER_EVENT] << "\t" << hnb_event[i]->nb_element
           << "\t" << 1. / (hnb_event[i]->mean + 1.) << endl;

        os << SEQ_label[SEQL_COMPLETE_INTER_EVENT] << "\t" << hnb_event[i]->mean *
              hnb_event[i]->nb_element << "\t"
           << hnb_event[i]->mean / (hnb_event[i]->mean + 1.) << endl;
        break;
      }

      case 'e' : {
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

    case 'o' : {
      os << "\n" << SEQ_label[SEQL_1_CENSORED_INTER_EVENT] << "\t" << mixture->nb_element
         << "\t" << 1. / (mixture->mean + 1.) << endl;

      os << SEQ_label[SEQL_COMPLETE_INTER_EVENT] << "\t" << mixture->mean *
            mixture->nb_element << "\t"
         << mixture->mean / (mixture->mean + 1.) << endl;
      break;
    }

    case 'e' : {
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


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet TimeEvents dans un fichier au format tableur.
 *
 *  arguments : reference sur un objet StatError, path.
 *
 *--------------------------------------------------------------*/

bool TimeEvents::spreadsheet_write(StatError &error , const char *path) const

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
    spreadsheet_write(out_file);
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Sortie Gnuplot d'un objet TimeEvents.
 *
 *  arguments : reference sur un objet StatError, prefixe des fichiers,
 *              titre des figures.
 *
 *--------------------------------------------------------------*/

bool TimeEvents::plot_write(StatError &error , const char *prefix ,
                            const char *title) const

{
  bool status;
  register int i , j , k;
  int nb_histo;
  const FrequencyDistribution **phisto;
  ostringstream data_file_name;


  error.init();

  // ecriture du fichier de donnees

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

  // ecriture du fichier de commandes et du fichier d'impression

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


/*--------------------------------------------------------------*
 *
 *  Sortie graphique d'un objet TimeEvents.
 *
 *--------------------------------------------------------------*/

MultiPlotSet* TimeEvents::get_plotable() const

{
  register int i , j , k , m;
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

    // vue : loi empirique des temps d'observation

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

  // vue : lois de comptage pour chaque temps d'observation

  nb_histo = 0;
  max_nb_value = 0;
  max_frequency = 0;

  for (j = htime->offset;j < htime->nb_value;j++) {
    if (htime->frequency[j] > 0) {
      nb_histo++;

      // calcul du nombre de valeurs maximum et de la frequence maximum

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

    // vue : lois empiriques de comptage superposees

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


/*--------------------------------------------------------------*
 *
 *  Recherche du temps minimum moyen entre 2 evenements.
 *
 *--------------------------------------------------------------*/

double TimeEvents::min_inter_event_computation() const

{
  register int i;
  int *ptime , *pnb_event;
  double ratio , min_ratio;


  ptime = time;
  pnb_event = nb_event;
  min_ratio = time[nb_class - 1];

  for (i = 0;i < nb_class;i++) {
    if (*pnb_event > 0) {
      ratio = (double)*ptime / (double)*pnb_event;
      if (ratio < min_ratio) {
        min_ratio = ratio;
      }
    }

    ptime++;
    pnb_event++;
  }

  return min_ratio;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de l'effectif total d'un objet TimeEvents.
 *
 *--------------------------------------------------------------*/

void TimeEvents::nb_element_computation()

{
  register int i;
  int *pfrequency;


  pfrequency = frequency;
  nb_element = 0;
  for (i = 0;i < nb_class;i++) {
    nb_element += *pfrequency++;
  }
}


/*--------------------------------------------------------------*
 *
 *  Constructeur par defaut de la classe RenewalData.
 *
 *--------------------------------------------------------------*/

RenewalData::RenewalData()

{
  renewal = NULL;

  type = 'v';

  length = NULL;
  sequence = NULL;

  inter_event = NULL;
  within = NULL;
  length_bias = NULL;
  backward = NULL;
  forward = NULL;

  index_event = NULL;
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe RenewalData.
 *
 *  arguments : nombre de fenetres d'observation, temps entre 2 dates d'observation.
 *
 *--------------------------------------------------------------*/

RenewalData::RenewalData(int nb_element , int itime)

{
  renewal = NULL;

  type = 'e';

  length = new int[nb_element];
  sequence = new int*[nb_element];

  inter_event = NULL;
  within = new FrequencyDistribution(itime);
  length_bias = NULL;
  backward = new FrequencyDistribution(itime);
  forward = new FrequencyDistribution(itime + 1);

  index_event = NULL;
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe RenewalData.
 *
 *  arguments : reference sur un objet TimeEvents,
 *              type de processus ('o' : ordinaire, 'e' : en equilibre).
 *
 *--------------------------------------------------------------*/

RenewalData::RenewalData(const TimeEvents &timev , int itype)
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


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe RenewalData.
 *
 *  arguments : type de processus ('o' : ordinaire, 'e' : en equilibre),
 *              reference sur un objet Renewal.
 *
 *--------------------------------------------------------------*/

RenewalData::RenewalData(int itype , const Renewal &renew)

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


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe RenewalData.
 *
 *  argument : nombre d'objets RenewalData,
 *             pointeurs sur les objets RenewalData.
 *
 *--------------------------------------------------------------*/

RenewalData::RenewalData(int nb_sample , const RenewalData **itimev)

{
  register int i , j , k;
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


/*--------------------------------------------------------------*
 *
 *  Copie d'un objet RenewalData.
 *
 *  arguments : reference sur un objet RenewalData,
 *              flag copie de l'objet Renewal.
 *
 *--------------------------------------------------------------*/

void RenewalData::copy(const RenewalData &timev , bool model_flag)

{
  register int i , j;
  int *psequence , *csequence;


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

    psequence = sequence[i];
    csequence = timev.sequence[i];
    for (j = 0;j < length[i];j++) {
      *psequence++ = *csequence++;
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


/*--------------------------------------------------------------*
 *
 *  Destruction des champs d'un objet RenewalData.
 *
 *--------------------------------------------------------------*/

void RenewalData::remove()

{
  delete renewal;

  delete [] length;

  if (sequence) {
    register int i;

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


/*--------------------------------------------------------------*
 *
 *  Destructeur de la classe RenewalData.
 *
 *--------------------------------------------------------------*/

RenewalData::~RenewalData()

{
  remove();
}


/*--------------------------------------------------------------*
 *
 *  Operateur d'assignement de la classe RenewalData.
 *
 *  argument : reference sur un objet RenewalData.
 *
 *--------------------------------------------------------------*/

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


/*--------------------------------------------------------------*
 *
 *  Fusion d'objets RenewalData.
 *
 *  arguments : reference sur un objet StatError, nombre d'objets RenewalData,
 *              pointeurs sur les objets RenewalData.
 *
 *--------------------------------------------------------------*/

RenewalData* RenewalData::merge(StatError &error , int nb_sample ,
                                const RenewalData **itimev) const

{
  bool status = true;
  register int i , j , k , m;
  int *psequence , *csequence;
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

    // copie des sequences

    i = 0;
    for (j = 0;j < nb_sample;j++) {
      for (k = 0;k < ptimev[j]->nb_element;k++) {
        psequence = timev->sequence[i];
        csequence = ptimev[j]->sequence[k];
        for (m = 0;m < ptimev[j]->length[k];m++) {
          *psequence++ = *csequence++;
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

    timev->build_index_event(timev->type == 'o' ? 0 : 1);

    delete [] phisto;
    delete [] ptimev;
  }

  return timev;
}


/*--------------------------------------------------------------*
 *
 *  Extraction d'une loi empirique.
 *
 *  arguments : reference sur un objet StatError, type de loi empirique,
 *              temps d'observation.
 *
 *--------------------------------------------------------------*/

DiscreteDistributionData* RenewalData::extract(StatError &error , int histo_type ,
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


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet RenewalData.
 *
 *  arguments : stream, flag niveau de detail, flag fichier.
 *
 *--------------------------------------------------------------*/

ostream& RenewalData::ascii_write(ostream &os , bool exhaustive , bool file_flag) const

{
  register int i , j;
  int nb_value , max , width[2];
  long old_adjust;


  old_adjust = os.setf(ios::right , ios::adjustfield);

  // ecriture des lois empiriques des intervalles de temps entre 2 evenements,
  // des intervalles de temps a l'interieur de la periode d'observation,
  // des intervalles de temps entre 2 evenements recouvrant une date d'observation,
  // des intervalles de temps apres le dernier evenement et des intervalles de temps residuel

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
    os << SEQ_label[SEQL_OBSERVATION_INTER_EVENT] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";
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
    os << SEQ_label[SEQL_BACKWARD] << " " << SEQ_label[SEQL_RECURRENCE_TIME] << " "
       << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";
    backward->ascii_characteristic_print(os , false , file_flag);

    os << "\n";
    if (file_flag) {
      os << "# ";
    }
    os << SEQ_label[SEQL_FORWARD] << " " << SEQ_label[SEQL_RECURRENCE_TIME] << " "
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
    os << " | " << SEQ_label[SEQL_OBSERVATION_INTER_EVENT] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
    if (length_bias) {
      os << " | " << SEQ_label[SEQL_LENGTH_BIASED] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
    }
    os << " | " << SEQ_label[SEQL_BACKWARD] << " " << SEQ_label[SEQL_RECURRENCE_TIME]
       << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " | " << SEQ_label[SEQL_FORWARD]
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
  switch (file_flag) {
  case false :
    TimeEvents::ascii_write(os , exhaustive , type);
    break;
  case true :
    ascii_file_write(os , exhaustive , type);
    break;
  }

  // ecriture des probabilites de non-evenement/evenement fonction du temps

  if (exhaustive) {
    os << "\n";
    if (file_flag) {
      os << "# ";
    }
    os << "   | " << SEQ_label[SEQL_OBSERVED] << " " << SEQ_label[SEQL_NO_EVENT_PROBABILITY]
       << " | " << SEQ_label[SEQL_OBSERVED] << " " << SEQ_label[SEQL_EVENT_PROBABILITY]
       << " | " << STAT_label[STATL_FREQUENCY] << endl;

    index_event->ascii_print(os , file_flag);

    // ecriture des sequences d'evenements

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

  os.setf((FMTFLAGS)old_adjust , ios::adjustfield);

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet RenewalData.
 *
 *  arguments : stream, flag niveau de detail.
 *
 *--------------------------------------------------------------*/

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


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet RenewalData dans un fichier.
 *
 *  arguments : reference sur un objet StatError, path,
 *              flag niveau de detail.
 *
 *--------------------------------------------------------------*/

bool RenewalData::ascii_write(StatError &error , const char *path ,
                              bool exhaustive) const

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

    if (renewal) {
      renewal->ascii_write(out_file , this , exhaustive , true);
    }
    else {
      ascii_write(out_file , exhaustive , true);
    }
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un processus de renouvellement et de la structure
 *  de donnees associee au format tableur.
 *
 *  argument : stream.
 *
 *--------------------------------------------------------------*/

ostream& RenewalData::spreadsheet_write(ostream &os) const

{
  register int i;
  int nb_value;


  // ecriture des loi empiriques des intervalles de temps entre 2 evenements,
  // des intervalles de temps a l'interieur de la periode d'observation
  // des intervalles de temps entre 2 evenements recouvrant une date d'observation,
  // des intervalles de temps apres le dernier evenement et des intervalles de temps residuel

  if (inter_event) {
    os << "\n" << SEQ_label[SEQL_INTER_EVENT] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\t";
    inter_event->spreadsheet_characteristic_print(os);
    os << STAT_label[STATL_VARIATION_COEFF] << "\t"
       << sqrt(inter_event->variance) / inter_event->mean << endl;
  }

  os << "\n" << SEQ_label[SEQL_OBSERVATION_INTER_EVENT] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\t";
  within->spreadsheet_characteristic_print(os);

  if (length_bias) {
    os << "\n" << SEQ_label[SEQL_LENGTH_BIASED] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\t";
    length_bias->spreadsheet_characteristic_print(os);
  }
  os << "\n" << SEQ_label[SEQL_BACKWARD] << " " << SEQ_label[SEQL_RECURRENCE_TIME] << " "
     << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\t";
  backward->spreadsheet_characteristic_print(os);

  os << "\n" << SEQ_label[SEQL_FORWARD] << " " << SEQ_label[SEQL_RECURRENCE_TIME] << " "
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
  os << "\t" << SEQ_label[SEQL_OBSERVATION_INTER_EVENT] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
  if (length_bias) {
    os << " \t" << SEQ_label[SEQL_LENGTH_BIASED] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
  }
  os << "\t" << SEQ_label[SEQL_BACKWARD] << " " << SEQ_label[SEQL_RECURRENCE_TIME]
     << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\t" << SEQ_label[SEQL_FORWARD]
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

  // ecriture des probabilites de non-evenement/evenement fonction du temps

  os << "\n\t" << SEQ_label[SEQL_OBSERVED] << " " << SEQ_label[SEQL_NO_EVENT_PROBABILITY]
     << "\t" << SEQ_label[SEQL_OBSERVED] << " " << SEQ_label[SEQL_EVENT_PROBABILITY]
     << "\t" << STAT_label[STATL_FREQUENCY] << endl;

  index_event->spreadsheet_print(os);

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet RenewalData dans un fichier au format tableur.
 *
 *  arguments : reference sur un objet StatError, path.
 *
 *--------------------------------------------------------------*/

bool RenewalData::spreadsheet_write(StatError &error , const char *path) const

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

    if (renewal) {
      renewal->spreadsheet_write(out_file , this);
    }
    else {
      spreadsheet_write(out_file);
    }
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Sortie Gnuplot d'un objet RenewalData.
 *
 *  arguments : reference sur un objet StatError, prefixe des fichiers,
 *              titre des figures.
 *
 *--------------------------------------------------------------*/

bool RenewalData::plot_write(StatError &error , const char *prefix ,
                             const char *title) const

{
  bool status = false;


  error.init();

  if (renewal) {
    status = renewal->plot_write(prefix , title , this);
  }

  else {
    register int i , j , k;
    int nb_histo;
    const FrequencyDistribution **phisto;
    ostringstream data_file_name[2];


    // ecriture des fichiers de donnees

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

      // ecriture du fichier de commandes et du fichier d'impression

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
                   << " title \"" << SEQ_label[SEQL_OBSERVATION_INTER_EVENT] << " "
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
                 << " title \"" << SEQ_label[SEQL_BACKWARD] << " "
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
                 << " title \"" << SEQ_label[SEQL_FORWARD] << " "
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


/*--------------------------------------------------------------*
 *
 *  Sortie graphique d'un objet RenewalData.
 *
 *--------------------------------------------------------------*/

MultiPlotSet* RenewalData::get_plotable() const

{
  MultiPlotSet *plot_set;


  if (renewal) {
    plot_set = renewal->get_plotable(this);
  }

  else {
    register int i , j , k , m;
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

      // vue : loi inter-evenement empirique

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
        legend << SEQ_label[SEQL_OBSERVATION_INTER_EVENT] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
        plot[i][1].legend = legend.str();

        plot[i][1].style = "impulses";

        within->plotable_frequency_write(plot[i][1]);
      }

      i++;
    }

    if ((within->nb_element > 0) || (length_bias)) {

      // vue : lois empiriques de l'intervalle de temps entre 2 evenements a l'interieur
      // de la periode d'observation et de l'intervalle de temps entre 2 evenements
      // recouvrant une date d'observation

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
        legend << SEQ_label[SEQL_OBSERVATION_INTER_EVENT] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
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

    // vue : lois empiriques de l'intervalle de temps apres le dernier evenement et
    // de l'intervalle de temps residuel

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
    legend << SEQ_label[SEQL_BACKWARD] << " " << SEQ_label[SEQL_RECURRENCE_TIME] << " "
           << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
    plot[i][0].legend = legend.str();

    plot[i][0].style = "impulses";

    backward->plotable_frequency_write(plot[i][0]);

    legend.str("");
    legend << SEQ_label[SEQL_FORWARD] << " " << SEQ_label[SEQL_RECURRENCE_TIME] << " "
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

      // vue : loi empirique des temps d'observation

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

    // vue : lois de comptage pour chaque temps d'observation

    if (htime->variance > 0.) {
      title.str("");
      title << SEQ_label[SEQL_NB_EVENT] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTIONS];
      plot[i].title = title.str();
    }

    // calcul du nombre de valeurs maximum et de la frequence maximum

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

      // vue : lois empiriques de comptage superposees

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

    // vue : probabilites de non-evenement/evenement fonction du temps

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


/*--------------------------------------------------------------*
 *
 *  Extraction des probabilites de non-evenement/evenement fonction du temps.
 *
 *  argument : decalage (0/1).
 *
 *--------------------------------------------------------------*/

void RenewalData::build_index_event(int offset)

{
  register int i , j;
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
