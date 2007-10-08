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



#include <math.h>
#include <sstream>
#include <iomanip>
#include <iostream>
#include "tool/rw_tokenizer.h"
#include "tool/rw_cstring.h"
#include "tool/rw_locale.h"
// #include <rw/vstream.h>
// #include <rw/rwfile.h>
#include "stat_tool/stat_tools.h"
#include "stat_tool/distribution.h"
#include "stat_tool/curves.h"
#include "stat_tool/stat_label.h"
#include "renewal.h"
#include "sequence_label.h"
#include "tool/config.h"

using namespace std;

extern int column_width(int value);
extern char* label(const char *file_name);



/*--------------------------------------------------------------*
 *
 *  Construction des histogrammes du temps et du nombre d'evenements
 *  a partir des echantillons {temps, nombre d'evenements, effectif}.
 *
 *--------------------------------------------------------------*/

void Time_events::build_histogram()

{
  register int i;
  int max_nb_event , *ptime , *pnb_event , *pfrequency;


  // construction de l'histogramme des temps et
  // creation des histogrammes du nombre d'evenements

  htime = new Histogram(time[nb_class - 1] + 1);

  hnb_event = new Histogram*[time[nb_class - 1] + 1];
  for (i = 0;i <= time[nb_class - 1];i++) {
    hnb_event[i] = 0;
  }

  ptime = time;
  pnb_event = nb_event;
  pfrequency = frequency;
  max_nb_event = 0;

  for (i = 0;i < nb_class - 1;i++) {
    if (*(ptime + 1) != *ptime) {
      hnb_event[*ptime] = new Histogram(*pnb_event + 1);
      if (*pnb_event > max_nb_event) {
        max_nb_event = *pnb_event;
      }
    }
    htime->frequency[*ptime++] += *pfrequency++;
    pnb_event++;
  }

  hnb_event[*ptime] = new Histogram(*pnb_event + 1);
  if (*pnb_event > max_nb_event) {
    max_nb_event = *pnb_event;
  }
  htime->frequency[*ptime] += *pfrequency;

  hmixture = new Histogram(max_nb_event + 1);

  htime->offset_computation();
  htime->nb_element = nb_element;
  htime->max_computation();
  htime->mean_computation();
  htime->variance_computation();

  // constrution des l'histogrammes du nombre d'evenements

  ptime = time;
  pnb_event = nb_event;
  pfrequency = frequency;

  for (i = 0;i < nb_class;i++) {
    hnb_event[*ptime++]->frequency[*pnb_event] += *pfrequency;
    hmixture->frequency[*pnb_event++] += *pfrequency++;
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

  hmixture->offset_computation();
  hmixture->nb_element = nb_element;
  hmixture->max_computation();
  hmixture->mean_computation();
  hmixture->variance_computation();
}


/*--------------------------------------------------------------*
 *
 *  Construction des echantillons {temps, nombre d'evenements, effectif}
 *  a partir des histogrammes du temps et du nombre d'evenements.
 *
 *--------------------------------------------------------------*/

void Time_events::build_sample()

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
 *  Construction d'un objet Time_events a partir d'echantillons
 *  {temps, nombre d'evenements}.
 *
 *  argument : nombre d'echantillons {temps, nombre d'evenements},
 *             pointeurs sur les temps et sur les nombres d'evenements.
 *
 *--------------------------------------------------------------*/

void Time_events::build(int inb_element , int *itime , int *inb_event)

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

  build_histogram();
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Time_events.
 *
 *  argument : nombre de classes.
 *
 *--------------------------------------------------------------*/

Time_events::Time_events(int inb_class)

{
  nb_element = 0;
  nb_class = inb_class;

  if (nb_class == 0) {
    time = 0;
    nb_event = 0;
    frequency = 0;
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

  htime = 0;
  hnb_event = 0;
  hmixture = 0;
}


/*--------------------------------------------------------------*
 *
 *  Copie d'un objet Time_events.
 *
 *  argument : reference sur un objet Time_events.
 *
 *--------------------------------------------------------------*/

void Time_events::copy(const Time_events &timev)

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

  // copie des histogrammes

  htime = new Histogram(*(timev.htime));

  hnb_event = new Histogram*[htime->nb_value];

  for (i = 0;i < htime->offset;i++) {
    hnb_event[i] = 0;
  }

  for (i = htime->offset;i < htime->nb_value;i++) {
    if (htime->frequency[i] > 0) {
      hnb_event[i] = new Histogram(*(timev.hnb_event[i]));
    }
    else {
      hnb_event[i] = 0;
    }
  }

  hmixture = new Histogram(*(timev.hmixture));
}


/*--------------------------------------------------------------*
 *
 *  Fusion d'objets Time_events.
 *
 *  argument : nombre d'objets Time_events,
 *             pointeurs sur les objets Time_events.
 *
 *--------------------------------------------------------------*/

void Time_events::merge(int nb_sample , const Time_events **ptimev)

{
  register int i , j;
  int nb_histo;
  const Histogram **phisto;


  nb_element = 0;
  for (i = 0;i < nb_sample;i++) {
    nb_element += ptimev[i]->nb_element;
  }

  phisto = new const Histogram*[nb_sample];

  // fusion des histogrammes du temps d'observation

  for (i = 0;i < nb_sample;i++) {
    phisto[i] = ptimev[i]->htime;
  }
  htime = new Histogram(nb_sample , phisto);

  // fusion des histogrammes du nombre d'evenements

  hnb_event = new Histogram*[htime->nb_value];

  for (i = 0;i < htime->offset;i++) {
    hnb_event[i] = 0;
  }

  for (i = htime->offset;i < htime->nb_value;i++) {
    if (htime->frequency[i]) {
      nb_histo = 0;
      for (j = 0;j < nb_sample;j++) {
        if ((i < ptimev[j]->htime->nb_value) && (ptimev[j]->hnb_event[i])) {
          phisto[nb_histo++] = ptimev[j]->hnb_event[i];
        }
      }
      hnb_event[i] = new Histogram(nb_histo , phisto);
    }

    else {
      hnb_event[i] = 0;
    }
  }

  for (i = 0;i < nb_sample;i++) {
    phisto[i] = ptimev[i]->hmixture;
  }
  hmixture = new Histogram(nb_sample , phisto);

  delete [] phisto;

  build_sample();
}


/*--------------------------------------------------------------*
 *
 *  Destruction des champs d'un objet Time_events.
 *
 *--------------------------------------------------------------*/

void Time_events::remove()

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
  delete hmixture;
}


/*--------------------------------------------------------------*
 *
 *  Destructeur de la classe Time_events.
 *
 *--------------------------------------------------------------*/

Time_events::~Time_events()

{
  remove();
}


/*--------------------------------------------------------------*
 *
 *  Operateur d'assignement de la classe Time_events.
 *
 *  argument : reference sur un objet Time_events.
 *
 *--------------------------------------------------------------*/

Time_events& Time_events::operator=(const Time_events &timev)

{
  if (&timev != this) {
    remove();
    copy(timev);
  }

  return *this;
}


/*--------------------------------------------------------------*
 *
 *  Extraction de l'histogramme du nombre d'evenements pour
 *  un temps d'observation donne.
 *
 *  arguments : reference sur un objet Format_error, type d'histogramme,
 *              temps d'observation.
 *
 *--------------------------------------------------------------*/

Distribution_data* Time_events::extract(Format_error &error , int histo_type , int itime) const

{
  Distribution_data *histo;


  error.init();

  if (histo_type == NB_EVENT) {
    if ((itime < htime->offset) || (itime >= htime->nb_value) || (htime->frequency[itime] == 0)) {
      histo = 0;
      error.update(SEQ_error[SEQR_OBSERVATION_TIME]);
    }
    else {
      histo = new Distribution_data(*hnb_event[itime]);
    }
  }

  else if (histo_type == MIXTURE) {
    histo = new Distribution_data(*hmixture);
  }

  return histo;
}


/*--------------------------------------------------------------*
 *
 *  Changement de l'unite de temps d'un objet Time_events.
 *
 *  arguments : reference sur un objet Format_error, facteur d'echelle.
 *
 *--------------------------------------------------------------*/

Time_events* Time_events::time_scaling(Format_error &error , int scaling_coeff) const

{
  bool status = true;
  register int i;
  int *ptime , *pnb_event , *pfrequency , *ttime , *tnb_event , *tfrequency;
  Time_events *timev;


  timev = 0;
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
    timev = new Time_events(nb_class);

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

    timev->build_histogram();
  }

  return timev;
}


/*--------------------------------------------------------------*
 *
 *  Selection d'echantillons {temps, nombre d'evenements}
 *  sur un critere de temps d'observation.
 *
 *  arguments : reference sur un objet Format_error,
 *              bornes sur le temps d'observation.
 *
 *--------------------------------------------------------------*/

Time_events* Time_events::time_select(Format_error &error , int min_time ,
                                      int max_time) const

{
  bool status = true;
  register int i;
  int bnb_class , *ptime , *pnb_event , *pfrequency , *ttime , *tnb_event , *tfrequency;
  Time_events *timev;


  timev = 0;
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

    timev = new Time_events(bnb_class);

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
      timev->build_histogram();
    }

    else {
      delete timev;
      timev = 0;
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
 *  arguments : reference sur un objet Format_error,
 *              bornes sur le nombre d'evenements.
 *
 *--------------------------------------------------------------*/

Time_events* Time_events::nb_event_select(Format_error &error , int min_nb_event ,
                                          int max_nb_event) const

{
  bool status = true;
  register int i;
  int bnb_class , *ptime , *pnb_event , *pfrequency , *ttime , *tnb_event , *tfrequency;
  Time_events *timev;


  timev = 0;
  error.init();

  if ((min_nb_event < 0) || (min_nb_event > max_nb_event)) {
    status = false;
    error.update(SEQ_error[SEQR_MIN_NB_EVENT]);
  }
  if ((max_nb_event < hmixture->offset) || (max_nb_event < min_nb_event)) {
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

    timev = new Time_events(bnb_class);

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
      timev->build_histogram();
    }

    else {
      delete timev;
      timev = 0;
      error.update(SEQ_error[SEQR_EMPTY_RENEWAL_DATA]);
    }
  }

  return timev;
}


/*--------------------------------------------------------------*
 *
 *  Construction d'un objet Time_events a partir d'un objet Histogram.
 *
 *  arguments : reference sur un objet Format_error, temps d'observation.
 *
 *--------------------------------------------------------------*/

Time_events* Histogram::build_time_events(Format_error &error , int itime) const

{
  bool status = true;
  register int i;
  int nb_class , *ptime , *pnb_event , *pfrequency , *hfrequency;
  Time_events *timev;


  timev = 0;
  error.init();

  if (itime / (nb_value - 1) < MIN_INTER_EVENT) {
    status = false;
    error.update(SEQ_error[SEQR_SHORT_OBSERVATION_TIME]);
  }
  if (itime > MAX_TIME) {
    status = false;
    error.update(SEQ_error[SEQR_LONG_OBSERVATION_TIME]);
  }

  if (status) {
    hfrequency = frequency + offset;
    nb_class = 0;
    for (i = offset;i < nb_value;i++) {
      if (*hfrequency++ > 0) {
        nb_class++;
      }
    }

    timev = new Time_events(nb_class);

    timev->nb_element = nb_element;

    // constitution des echantillons

    ptime = timev->time;
    pnb_event = timev->nb_event;
    pfrequency = timev->frequency;
    hfrequency = frequency + offset;

    for (i = offset;i < nb_value;i++) {
      if (*hfrequency > 0) {
        *ptime++ = itime;
        *pnb_event++ = i;
        *pfrequency++ = *hfrequency;
      }
      hfrequency++;
    }

    // construction des histogrammes

    timev->htime = new Histogram(itime + 1);

    timev->htime->frequency[itime] = nb_element;
    timev->htime->offset = itime;
    timev->htime->nb_element = nb_element;
    timev->htime->max = nb_element;
    timev->htime->mean = itime;
    timev->htime->variance = 0.;

    timev->hnb_event = new Histogram*[itime + 1];
    for (i = 0;i < itime;i++) {
      timev->hnb_event[i] = 0;
    }
    timev->hnb_event[itime] = new Histogram(*this);

    timev->hmixture = new Histogram(*this);
  }

  return timev;
}


/*--------------------------------------------------------------*
 *
 *  Construction d'un objet Time_events a partir d'un fichier.
 *  Format : n lignes de la forme (temps > 0) (nombre d'evenements >= 0)
 *           (effectif >= 0).
 *
 *  arguments : reference sur un objet Format_error, path.
 *
 *--------------------------------------------------------------*/

Time_events* time_events_ascii_read(Format_error &error , const char *path)

{
  RWLocaleSnapshot locale("en");
  RWCString buffer , token;
  size_t position;
  bool status , lstatus;
  register int i , j;
  int line , nb_class , nb_element;
  long value , time , nb_event;
  Time_events *timev;
  ifstream in_file(path);


  timev = 0;
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

      timev = new Time_events(nb_class);
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

      timev->build_histogram();
    }
  }

  return timev;
}


/*--------------------------------------------------------------*
 *
 *  Construction d'un objet Time_events a partir d'un fichier.
 *  Format : n lignes de la forme (temps > 0) (nombre d'evenements >= 0).
 *
 *  arguments : reference sur un objet Format_error, path.
 *
 *--------------------------------------------------------------*/

Time_events* old_time_events_ascii_read(Format_error &error , const char *path)

{
  RWLocaleSnapshot locale("en");
  RWCString buffer , token;
  size_t position;
  bool status , lstatus;
  register int i , j;
  int line , nb_element , *ptime , *pnb_event;
  long value;
  Time_events *timev;
  ifstream in_file(path);


  timev = 0;
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

      timev = new Time_events(nb_element , ptime , pnb_event);

      delete [] ptime;
      delete [] pnb_event;
    }
  }

  return timev;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture sur une ligne d'un objet Time_events.
 *
 *  argument : stream.
 *
 *--------------------------------------------------------------*/

ostream& Time_events::line_write(ostream &os) const

{
  os << STAT_label[STATL_SAMPLE_SIZE] << ": " << nb_element << "   "
     << SEQ_label[SEQL_OBSERVATION_TIME] << " " << STAT_label[STATL_MEAN] << ": " << htime->mean << "   "
     << SEQ_label[SEQL_NB_EVENT] << " " << STAT_label[STATL_MEAN] << ": " << hmixture->mean;

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Time_events.
 *
 *  arguments : stream, flag niveau de detail, type de processus.
 *
 *--------------------------------------------------------------*/

ostream& Time_events::ascii_write(ostream &os , bool exhaustive , char type) const

{
  register int i;


  if ((htime->variance > 0.) && (exhaustive)) {
    os << SEQ_label[SEQL_OBSERVATION_TIME] << " " << STAT_label[STATL_HISTOGRAM] << " - ";
    htime->ascii_characteristic_print(os);

    os << "\n   | " << SEQ_label[SEQL_OBSERVATION_TIME] << " " << STAT_label[STATL_HISTOGRAM] << endl;
    htime->ascii_print(os);
  }

  for (i = htime->offset;i < htime->nb_value;i++) {
    if (htime->frequency[i] > 0) {
      if (((exhaustive) && (htime->variance > 0.)) || (i > htime->offset)) {
        os << "\n";
      }
      os << SEQ_label[SEQL_NB_EVENT] << " " << SEQ_label[SEQL_DURING] << " "
         << i << " " << SEQ_label[SEQL_TIME_UNIT] << " " << STAT_label[STATL_HISTOGRAM] << " - ";
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
           << i << " " << SEQ_label[SEQL_TIME_UNIT] << " " << STAT_label[STATL_HISTOGRAM] << endl;
        hnb_event[i]->ascii_print(os);
      }
    }
  }

  if ((htime->variance > 0.) && (exhaustive)) {
    os << "\n" << SEQ_label[SEQL_NB_EVENT] << " " << STAT_label[STATL_HISTOGRAM] << " - ";
    hmixture->ascii_characteristic_print(os);

    switch (type) {

    case 'o' : {
      os << "\n" << SEQ_label[SEQL_1_CENSORED_INTER_EVENT] << ": " << hmixture->nb_element
         << " (" << 1. / (hmixture->mean + 1.) << ")" << endl;

      os << SEQ_label[SEQL_COMPLETE_INTER_EVENT] << ": " << hmixture->mean *
            hmixture->nb_element << " ("
         << hmixture->mean / (hmixture->mean + 1.) << ")" << endl;
      break;
    }

    case 'e' : {
      os << "\n" << SEQ_label[SEQL_2_CENSORED_INTER_EVENT] << ": " << hmixture->frequency[0]
         << " (" << hmixture->frequency[0] / (hmixture->nb_element * (hmixture->mean + 1.))
         << ")" << endl;

      os << SEQ_label[SEQL_1_CENSORED_INTER_EVENT] << ": " << (hmixture->nb_element -
             hmixture->frequency[0]) * 2 << " ("
         << (hmixture->nb_element - hmixture->frequency[0]) * 2. /
            (hmixture->nb_element * (hmixture->mean + 1.)) << ")" << endl;

      os << SEQ_label[SEQL_COMPLETE_INTER_EVENT] << ": " << (hmixture->mean - 1.) *
            hmixture->nb_element + hmixture->frequency[0] << " ("
         << (hmixture->mean - 1. + (double)hmixture->frequency[0] /
             (double)hmixture->nb_element) / (hmixture->mean + 1.) << ")" << endl;
      break;
    }
    }

    os << "\n   | " << SEQ_label[SEQL_NB_EVENT] << " " << STAT_label[STATL_HISTOGRAM] << endl;
    hmixture->ascii_print(os);
  }

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Time_events.
 *
 *  arguments : stream, flag niveau de detail.
 *
 *--------------------------------------------------------------*/

ostream& Time_events::ascii_write(ostream &os , bool exhaustive) const

{
  return ascii_write(os , exhaustive , 'v');
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Time_events dans un fichier.
 *
 *  arguments : reference sur un objet Format_error, path,
 *              flag niveau de detail, type de processus.
 *
 *--------------------------------------------------------------*/

ostream& Time_events::ascii_file_write(ostream &os , bool exhaustive , char type) const

{
  register int i;
  int max_frequency , *pfrequency , width[3];


  if ((htime->variance > 0.) && (exhaustive)) {
    os << "# " << SEQ_label[SEQL_OBSERVATION_TIME] << " " << STAT_label[STATL_HISTOGRAM] << " - ";
    htime->ascii_characteristic_print(os , false , true);

    os << "\n#    | " << SEQ_label[SEQL_OBSERVATION_TIME] << " " << STAT_label[STATL_HISTOGRAM] << endl;
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
         << time[i] << " " << SEQ_label[SEQL_TIME_UNIT] << " " << STAT_label[STATL_HISTOGRAM] << " - ";
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
    os << "\n# " << SEQ_label[SEQL_NB_EVENT] << " " << STAT_label[STATL_HISTOGRAM] << " - ";
    hmixture->ascii_characteristic_print(os , false , true);

    switch (type) {

    case 'o' : {
      os << "\n# " << SEQ_label[SEQL_1_CENSORED_INTER_EVENT] << ": " << hmixture->nb_element
         << " (" << 1. / (hmixture->mean + 1.) << ")" << endl;

      os << "# " << SEQ_label[SEQL_COMPLETE_INTER_EVENT] << ": " << hmixture->mean *
            hmixture->nb_element << " ("
         << hmixture->mean / (hmixture->mean + 1.) << ")" << endl;
      break;
    }

    case 'e' : {
      os << "\n# " << SEQ_label[SEQL_2_CENSORED_INTER_EVENT] << ": " << hmixture->frequency[0]
         << " (" << hmixture->frequency[0] / (hmixture->nb_element * (hmixture->mean + 1.))
         << ")" << endl;

      os << "# " << SEQ_label[SEQL_1_CENSORED_INTER_EVENT] << ": " << (hmixture->nb_element -
             hmixture->frequency[0]) * 2 << " ("
         << (hmixture->nb_element - hmixture->frequency[0]) * 2. /
            (hmixture->nb_element * (hmixture->mean + 1.)) << ")" << endl;

      os << "# " << SEQ_label[SEQL_COMPLETE_INTER_EVENT] << ": " << (hmixture->mean - 1.) *
            hmixture->nb_element + hmixture->frequency[0] << " ("
         << (hmixture->mean - 1. + (double)hmixture->frequency[0] /
             (double)hmixture->nb_element) / (hmixture->mean + 1.) << ")" << endl;
      break;
    }
    }

    os << "\n#    | " << SEQ_label[SEQL_NB_EVENT] << " " << STAT_label[STATL_HISTOGRAM] << endl;
    hmixture->ascii_print(os , true);
  }

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Time_events dans un fichier.
 *
 *  arguments : reference sur un objet Format_error, path,
 *              flag niveau de detail.
 *
 *--------------------------------------------------------------*/

bool Time_events::ascii_write(Format_error &error , const char *path ,
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
 *  Ecriture d'un objet Time_events au format tableur.
 *
 *  argument : stream, type de processus.
 *
 *--------------------------------------------------------------*/

ostream& Time_events::spreadsheet_write(ostream &os , char type) const

{
  register int i;


  if (htime->variance > 0.) {
    os << SEQ_label[SEQL_OBSERVATION_TIME] << " " << STAT_label[STATL_HISTOGRAM] << "\t";
    htime->spreadsheet_characteristic_print(os);

    os << "\n\t" << SEQ_label[SEQL_OBSERVATION_TIME] << " " << STAT_label[STATL_HISTOGRAM] << endl;
    htime->spreadsheet_print(os);
  }

  for (i = htime->offset;i < htime->nb_value;i++) {
    if (htime->frequency[i] > 0) {
      if ((htime->variance > 0.) || (i > htime->offset)) {
        os << "\n";
      }
      os << SEQ_label[SEQL_NB_EVENT] << " " << SEQ_label[SEQL_DURING] << " "
         << i << " " << SEQ_label[SEQL_TIME_UNIT] << " " << STAT_label[STATL_HISTOGRAM] << "\t";
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
         << i << " " << SEQ_label[SEQL_TIME_UNIT] << " " << STAT_label[STATL_HISTOGRAM] << endl;
      hnb_event[i]->spreadsheet_print(os);
    }
  }

  if (htime->variance > 0.) {
    os << "\n" << SEQ_label[SEQL_NB_EVENT] << " " << STAT_label[STATL_HISTOGRAM] << "\t";
    hmixture->spreadsheet_characteristic_print(os);

    switch (type) {

    case 'o' : {
      os << "\n" << SEQ_label[SEQL_1_CENSORED_INTER_EVENT] << "\t" << hmixture->nb_element
         << "\t" << 1. / (hmixture->mean + 1.) << endl;

      os << SEQ_label[SEQL_COMPLETE_INTER_EVENT] << "\t" << hmixture->mean *
            hmixture->nb_element << "\t"
         << hmixture->mean / (hmixture->mean + 1.) << endl;
      break;
    }

    case 'e' : {
      os << "\n" << SEQ_label[SEQL_2_CENSORED_INTER_EVENT] << "\t" << hmixture->frequency[0] << "\t"
         << hmixture->frequency[0] / (hmixture->nb_element * (hmixture->mean + 1.)) << endl;

      os << SEQ_label[SEQL_1_CENSORED_INTER_EVENT] << "\t" << (hmixture->nb_element -
             hmixture->frequency[0]) * 2 << "\t"
         << (hmixture->nb_element - hmixture->frequency[0]) * 2. /
            (hmixture->nb_element * (hmixture->mean + 1.)) << endl;

      os << SEQ_label[SEQL_COMPLETE_INTER_EVENT] << "\t" << (hmixture->mean - 1.) *
            hmixture->nb_element + hmixture->frequency[0] << "\t"
         << (hmixture->mean - 1. + (double)hmixture->frequency[0] /
             (double)hmixture->nb_element) / (hmixture->mean + 1.) << endl;
      break;
    }
    }

    os << "\n\t" << SEQ_label[SEQL_NB_EVENT] << " " << STAT_label[STATL_HISTOGRAM] << endl;
    hmixture->spreadsheet_print(os);
  }

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Time_events dans un fichier au format tableur.
 *
 *  arguments : reference sur un objet Format_error, path.
 *
 *--------------------------------------------------------------*/

bool Time_events::spreadsheet_write(Format_error &error , const char *path) const

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
 *  Sortie Gnuplot d'un objet Time_events.
 *
 *  arguments : reference sur un objet Format_error, prefixe des fichiers,
 *              titre des figures.
 *
 *--------------------------------------------------------------*/

bool Time_events::plot_write(Format_error &error , const char *prefix ,
                             const char *title) const

{
  bool status;
  register int i , j , k;
  int nb_histo;
  const Histogram **phisto;
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

  phisto = new const Histogram*[nb_histo];

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
    phisto[nb_histo++] = hmixture;
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
                 << STAT_label[STATL_HISTOGRAM] << "\" with impulses" << endl;

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
                   << STAT_label[STATL_HISTOGRAM] << "\" with impulses" << endl;

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

        if (hmixture->nb_value - 1 < TIC_THRESHOLD) {
          out_file << "set xtics 0,1" << endl;
        }
        if ((int)(hmixture->max * YSCALE) + 1 < TIC_THRESHOLD) {
          out_file << "set ytics 0,1" << endl;
        }

        out_file << "plot [0:" << hmixture->nb_value - 1 << "] [0:"
                 << (int)(hmixture->max * YSCALE) + 1 << "] \""
                 << label((data_file_name.str()).c_str()) << "\" using " << j
                 << " title \"" << SEQ_label[SEQL_NB_EVENT] << " "
                 << STAT_label[STATL_HISTOGRAM] << "\" with impulses" << endl;

        if (hmixture->nb_value - 1 < TIC_THRESHOLD) {
          out_file << "set xtics autofreq" << endl;
        }
        if ((int)(hmixture->max * YSCALE) + 1 < TIC_THRESHOLD) {
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
 *  Recherche du temps minimum moyen entre 2 evenements.
 *
 *--------------------------------------------------------------*/

double Time_events::min_inter_event_computation() const

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
 *  Fonctions pour la persistance.
 *
 *--------------------------------------------------------------*/

/* RWDEFINE_COLLECTABLE(Time_events , STATI_TIME_EVENTS);


RWspace Time_events::binaryStoreSize() const

{
  register int i;
  RWspace size;


  size = sizeof(nb_element) + sizeof(nb_class) +
         sizeof(*time) * nb_class + sizeof(*nb_event) * nb_class +
         sizeof(*frequency) * nb_class + htime->binaryStoreSize();

  for (i = htime->offset;i < htime->nb_value;i++) {
    if (htime->frequency[i] > 0) {
      size += hnb_event[i]->binaryStoreSize();
    }
  }

  size += hmixture->binaryStoreSize();

  return size;
}


void Time_events::restoreGuts(RWvistream &is)

{
  register int i;
  int *ptime , *pnb_event , *pfrequency;


  remove();

  is >> nb_element >> nb_class;

  time = new int[nb_class];
  ptime = time;
  for (i = 0;i < nb_class;i++) {
    is >> *ptime++;
  }

  nb_event = new int[nb_class];
  pnb_event = nb_event;
  for (i = 0;i < nb_class;i++) {
    is >> *pnb_event++;
  }

  frequency = new int[nb_class];
  pfrequency = frequency;
  for (i = 0;i < nb_class;i++) {
    is >> *pfrequency++;
  }

  htime = new Histogram();
  htime->restoreGuts(is);

  hnb_event = new Histogram*[htime->nb_value];

  for (i = 0;i < htime->offset;i++) {
    hnb_event[i] = 0;
  }
  for (i = htime->offset;i < htime->nb_value;i++) {
    if (htime->frequency[i] > 0) {
      hnb_event[i] = new Histogram();
      hnb_event[i]->restoreGuts(is);
    }
    else {
      hnb_event[i] = 0;
    }
  }

  hmixture = new Histogram();
  hmixture->restoreGuts(is);
}


void Time_events::restoreGuts(RWFile &file)

{
  register int i;


  remove();

  file.Read(nb_element);
  file.Read(nb_class);

  time = new int[nb_class];
  file.Read(time , nb_class);
  nb_event = new int[nb_class];
  file.Read(nb_event , nb_class);
  frequency = new int[nb_class];
  file.Read(frequency , nb_class);

  htime = new Histogram();
  htime->restoreGuts(file);

  hnb_event = new Histogram*[htime->nb_value];

  for (i = 0;i < htime->offset;i++) {
    hnb_event[i] = 0;
  }
  for (i = htime->offset;i < htime->nb_value;i++) {
    if (htime->frequency[i] > 0) {
      hnb_event[i] = new Histogram();
      hnb_event[i]->restoreGuts(file);
    }
    else {
      hnb_event[i] = 0;
    }
  }

  hmixture = new Histogram();
  hmixture->restoreGuts(file);
}


void Time_events::saveGuts(RWvostream &os) const

{
  register int i;
  int *ptime , *pnb_event , *pfrequency;


  os << nb_element << nb_class;

  ptime = time;
  for (i = 0;i < nb_class;i++) {
    os << *ptime++;
  }

  pnb_event = nb_event;
  for (i = 0;i < nb_class;i++) {
    os << *pnb_event++;
  }

  pfrequency = frequency;
  for (i = 0;i < nb_class;i++) {
    os << *pfrequency++;
  }

  htime->saveGuts(os);

  for (i = htime->offset;i < htime->nb_value;i++) {
    if (htime->frequency[i] > 0) {
      hnb_event[i]->saveGuts(os);
    }
  }

  hmixture->saveGuts(os);
}


void Time_events::saveGuts(RWFile &file) const

{
  register int i;


  file.Write(nb_element);
  file.Write(nb_class);

  file.Write(time , nb_class);
  file.Write(nb_event , nb_class);
  file.Write(frequency , nb_class);

  htime->saveGuts(file);

  for (i = htime->offset;i < htime->nb_value;i++) {
    if (htime->frequency[i] > 0) {
      hnb_event[i]->saveGuts(file);
    }
  }

  hmixture->saveGuts(file);
} */


/*--------------------------------------------------------------*
 *
 *  Calcul de l'effectif total d'un objet Time_events.
 *
 *--------------------------------------------------------------*/

void Time_events::nb_element_computation()

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
 *  Constructeur par defaut de la classe Renewal_data.
 *
 *--------------------------------------------------------------*/

Renewal_data::Renewal_data()

{
  renewal = 0;

  type = 'v';

  length = 0;
  sequence = 0;

  inter_event = 0;
  within = 0;
  length_bias = 0;
  backward = 0;
  forward = 0;

  index_event = 0;
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Renewal_data.
 *
 *  arguments : nombre de fenetres d'observation, temps entre 2 dates d'observation.
 *
 *--------------------------------------------------------------*/

Renewal_data::Renewal_data(int nb_element , int itime)

{
  renewal = 0;

  type = 'e';

  length = new int[nb_element];
  sequence = new int*[nb_element];

  inter_event = 0;
  within = new Histogram(itime);
  length_bias = 0;
  backward = new Histogram(itime);
  forward = new Histogram(itime + 1);

  index_event = 0;
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Renewal_data.
 *
 *  arguments : reference sur un objet Time_events,
 *              type de processus ('o' : ordinaire, 'e' : en equilibre).
 *
 *--------------------------------------------------------------*/

Renewal_data::Renewal_data(const Time_events &timev , int itype)
:Time_events(timev)

{
  renewal = 0;

  type = itype;

  length = 0;
  sequence = 0;

  inter_event = 0;
  within = 0;
  length_bias = 0;
  backward = 0;
  forward = 0;

  index_event = 0;
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Renewal_data.
 *
 *  arguments : type de processus ('o' : ordinaire, 'e' : en equilibre),
 *              reference sur un objet Renewal.
 *
 *--------------------------------------------------------------*/

Renewal_data::Renewal_data(int itype , const Renewal &renew)

{
  renewal = 0;

  type = itype;

  length = 0;
  sequence = 0;

  inter_event = new Histogram(*(renew.inter_event));
  within = new Histogram(*(renew.inter_event));
  length_bias = new Histogram(*(renew.length_bias));
  backward = new Histogram(renew.backward->alloc_nb_value);
  forward = new Histogram(*(renew.forward));

  index_event = 0;
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Renewal_data.
 *
 *  argument : nombre d'objets Renewal_data,
 *             pointeurs sur les objets Renewal_data.
 *
 *--------------------------------------------------------------*/

Renewal_data::Renewal_data(int nb_sample , const Renewal_data **itimev)

{
  register int i , j , k;
  const Time_events **ptimev;


  ptimev = new const Time_events*[nb_sample];

  for (i = 0;i < nb_sample;i++) {
    ptimev[i] = itimev[i];
  }
  Time_events::merge(nb_sample , ptimev);

  delete [] ptimev;

  renewal = 0;

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

  inter_event = 0;
  within = 0;
  length_bias = 0;
  backward = 0;
  forward = 0;

  index_event = 0;
}


/*--------------------------------------------------------------*
 *
 *  Copie d'un objet Renewal_data.
 *
 *  arguments : reference sur un objet Renewal_data,
 *              flag copie de l'objet Renewal.
 *
 *--------------------------------------------------------------*/

void Renewal_data::copy(const Renewal_data &timev , bool model_flag)

{
  register int i , j;
  int *psequence , *csequence;


  if ((model_flag) && (timev.renewal)) {
    renewal = new Renewal(*(timev.renewal) , false);
  }
  else {
    renewal = 0;
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
    inter_event = new Histogram(*(timev.inter_event));
  }
  else {
    inter_event = 0;
  }
  within = new Histogram(*(timev.within));
  if (timev.length_bias) {
    length_bias = new Histogram(*(timev.length_bias));
  }
  else {
    length_bias = 0;
  }
  backward = new Histogram(*(timev.backward));
  forward = new Histogram(*(timev.forward));

  index_event = new Curves(*(timev.index_event));
}


/*--------------------------------------------------------------*
 *
 *  Destruction des champs d'un objet Renewal_data.
 *
 *--------------------------------------------------------------*/

void Renewal_data::remove()

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
 *  Destructeur de la classe Renewal_data.
 *
 *--------------------------------------------------------------*/

Renewal_data::~Renewal_data()

{
  remove();
}


/*--------------------------------------------------------------*
 *
 *  Operateur d'assignement de la classe Renewal_data.
 *
 *  argument : reference sur un objet Renewal_data.
 *
 *--------------------------------------------------------------*/

Renewal_data& Renewal_data::operator=(const Renewal_data &timev)

{
  if (&timev != this) {
    remove();
    Time_events::remove();

    Time_events::copy(timev);
    copy(timev);
  }

  return *this;
}


/*--------------------------------------------------------------*
 *
 *  Fusion d'objets Renewal_data.
 *
 *  arguments : reference sur un objet Format_error, nombre d'objets Renewal_data,
 *              pointeurs sur les objets Renewal_data.
 *
 *--------------------------------------------------------------*/

Renewal_data* Renewal_data::merge(Format_error &error , int nb_sample , const Renewal_data **itimev) const

{
  bool status = true;
  register int i , j , k , m;
  int *psequence , *csequence;
  const Histogram **phisto;
  Renewal_data *timev;
  const Renewal_data **ptimev;


  timev = 0;
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
    ptimev = new const Renewal_data*[nb_sample];

    ptimev[0] = this;
    for (i = 1;i < nb_sample;i++) {
      ptimev[i] = itimev[i - 1];
    }

    timev = new Renewal_data(nb_sample , ptimev);

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

    phisto = new const Histogram*[nb_sample];

    if (inter_event) {
      for (i = 0;i < nb_sample;i++) {
        phisto[i] = ptimev[i]->inter_event;
      }
      timev->inter_event = new Histogram(nb_sample , phisto);
    }

    for (i = 0;i < nb_sample;i++) {
      phisto[i] = ptimev[i]->within;
    }
    timev->within = new Histogram(nb_sample , phisto);

    if (length_bias) {
      for (i = 0;i < nb_sample;i++) {
        phisto[i] = ptimev[i]->length_bias;
      }
      timev->length_bias = new Histogram(nb_sample , phisto);
    }

    for (i = 0;i < nb_sample;i++) {
      phisto[i] = ptimev[i]->backward;
    }
    timev->backward = new Histogram(nb_sample , phisto);

    for (i = 0;i < nb_sample;i++) {
      phisto[i] = ptimev[i]->forward;
    }
    timev->forward = new Histogram(nb_sample , phisto);

    timev->build_index_event(timev->type == 'o' ? 0 : 1);

    delete [] phisto;
    delete [] ptimev;
  }

  return timev;
}


/*--------------------------------------------------------------*
 *
 *  Extraction d'un histogramme.
 *
 *  arguments : reference sur un objet Format_error, type d'histogramme,
 *              temps d'observation.
 *
 *--------------------------------------------------------------*/

Distribution_data* Renewal_data::extract(Format_error &error , int histo_type , int itime) const

{
  bool status = true;
  Distribution *pdist;
  Parametric *pparam;
  Histogram *phisto;
  Distribution_data *histo;


  error.init();

  if (histo_type == NB_EVENT) {
    if ((itime < htime->offset) || (itime >= htime->nb_value) || (htime->frequency[itime] == 0)) {
      histo = 0;
      error.update(SEQ_error[SEQR_OBSERVATION_TIME]);
    }
    else {
      histo = new Distribution_data(*hnb_event[itime] ,
                                    (renewal ? renewal->nb_event[itime] : 0));
    }
  }

  else {
    histo = 0;

    switch (histo_type) {

    case INTER_EVENT : {
      if (inter_event) {
        phisto = inter_event;
      }
      else {
        status = false;
        ostringstream error_message;
        error_message << SEQ_label[SEQL_INTER_EVENT] << " " << STAT_label[STATL_HISTOGRAM] << " "
                      << SEQ_error[SEQR_NOT_PRESENT];
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
        error_message << SEQ_label[SEQL_LENGTH_BIASED] << " " << STAT_label[STATL_HISTOGRAM] << " "
                      << SEQ_error[SEQR_NOT_PRESENT];
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

    case MIXTURE : {
      phisto = hmixture;
      break;
    }
    }

    if ((status) && (phisto->nb_element == 0)) {
      status = false;
      error.update(STAT_error[STATR_EMPTY_HISTOGRAM]);
    }

    if (status) {
      pdist = 0;
      pparam = 0;

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
        case MIXTURE :
          pdist = renewal->mixture;
          break;
        }
      }

      if (pdist) {
        histo = new Distribution_data(*phisto , pdist);
      }
      else {
        histo = new Distribution_data(*phisto , pparam);
      }
    }
  }

  return histo;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Renewal_data.
 *
 *  arguments : stream, flag niveau de detail, flag fichier.
 *
 *--------------------------------------------------------------*/

ostream& Renewal_data::ascii_write(ostream &os , bool exhaustive , bool file_flag) const

{
  register int i , j;
  int nb_value , max , width[2];
  long old_adjust;


  old_adjust = os.setf(ios::right , ios::adjustfield);

  // ecriture des histogrammes des intervalles de temps entre 2 evenements,
  // des intervalles de temps a l'interieur de la periode d'observation,
  // des intervalles de temps entre 2 evenements recouvrant une date d'observation,
  // des intervalles de temps apres le dernier evenement, des intervalles de temps residuel

  if (inter_event) {
    os << SEQ_label[SEQL_INTER_EVENT] << " " << STAT_label[STATL_HISTOGRAM] << " - ";
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
    os << SEQ_label[SEQL_OBSERVATION_INTER_EVENT] << " " << STAT_label[STATL_HISTOGRAM] << " - ";
    within->ascii_characteristic_print(os , false , file_flag);
  }

  if (exhaustive) {
    if (length_bias) {
      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << SEQ_label[SEQL_LENGTH_BIASED] << " " << STAT_label[STATL_HISTOGRAM] << " - ";
      length_bias->ascii_characteristic_print(os , false , file_flag);
    }

    os << "\n";
    if (file_flag) {
      os << "# ";
    }
    os << SEQ_label[SEQL_BACKWARD] << " " << SEQ_label[SEQL_RECURRENCE_TIME] << " "
       << STAT_label[STATL_HISTOGRAM] << " - ";
    backward->ascii_characteristic_print(os , false , file_flag);

    os << "\n";
    if (file_flag) {
      os << "# ";
    }
    os << SEQ_label[SEQL_FORWARD] << " " << SEQ_label[SEQL_RECURRENCE_TIME] << " "
       << STAT_label[STATL_HISTOGRAM] << " - ";
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
      os << " | " << SEQ_label[SEQL_INTER_EVENT] << " " << STAT_label[STATL_HISTOGRAM];
    }
    os << " | " << SEQ_label[SEQL_OBSERVATION_INTER_EVENT] << " " << STAT_label[STATL_HISTOGRAM];
    if (length_bias) {
      os << " | " << SEQ_label[SEQL_LENGTH_BIASED] << " " << STAT_label[STATL_HISTOGRAM];
    }
    os << " | " << SEQ_label[SEQL_BACKWARD] << " " << SEQ_label[SEQL_RECURRENCE_TIME]
       << " " << STAT_label[STATL_HISTOGRAM] << " | " << SEQ_label[SEQL_FORWARD]
       << " " << SEQ_label[SEQL_RECURRENCE_TIME] << " " << STAT_label[STATL_HISTOGRAM];

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
    Time_events::ascii_write(os , exhaustive , type);
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
 *  Ecriture d'un objet Renewal_data.
 *
 *  arguments : stream, flag niveau de detail.
 *
 *--------------------------------------------------------------*/

ostream& Renewal_data::ascii_write(ostream &os , bool exhaustive) const

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
 *  Ecriture d'un objet Renewal_data dans un fichier.
 *
 *  arguments : reference sur un objet Format_error, path,
 *              flag niveau de detail.
 *
 *--------------------------------------------------------------*/

bool Renewal_data::ascii_write(Format_error &error , const char *path ,
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

ostream& Renewal_data::spreadsheet_write(ostream &os) const

{
  register int i;
  int nb_value;


  // ecriture de l'histogramme des intervalles de temps entre 2 evenements,
  // de l'histogramme des intervalles de temps a l'interieur de la periode d'observation
  // de l'histogramme des intervalles de temps entre 2 evenements recouvrant une date d'observation,
  // de l'histogramme des intervalles de temps apres le dernier evenement
  // de l'histogramme des intervalles de temps residuel

  if (inter_event) {
    os << "\n" << SEQ_label[SEQL_INTER_EVENT] << " " << STAT_label[STATL_HISTOGRAM] << "\t";
    inter_event->spreadsheet_characteristic_print(os);
    os << STAT_label[STATL_VARIATION_COEFF] << "\t"
       << sqrt(inter_event->variance) / inter_event->mean << endl;
  }

  os << "\n" << SEQ_label[SEQL_OBSERVATION_INTER_EVENT] << " " << STAT_label[STATL_HISTOGRAM] << "\t";
  within->spreadsheet_characteristic_print(os);

  if (length_bias) {
    os << "\n" << SEQ_label[SEQL_LENGTH_BIASED] << " " << STAT_label[STATL_HISTOGRAM] << "\t";
    length_bias->spreadsheet_characteristic_print(os);
  }
  os << "\n" << SEQ_label[SEQL_BACKWARD] << " " << SEQ_label[SEQL_RECURRENCE_TIME] << " "
     << STAT_label[STATL_HISTOGRAM] << "\t";
  backward->spreadsheet_characteristic_print(os);

  os << "\n" << SEQ_label[SEQL_FORWARD] << " " << SEQ_label[SEQL_RECURRENCE_TIME] << " "
     << STAT_label[STATL_HISTOGRAM] << "\t";
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
    os << "\t" << SEQ_label[SEQL_INTER_EVENT] << " " << STAT_label[STATL_HISTOGRAM];
  }
  os << "\t" << SEQ_label[SEQL_OBSERVATION_INTER_EVENT] << " " << STAT_label[STATL_HISTOGRAM];
  if (length_bias) {
    os << " \t" << SEQ_label[SEQL_LENGTH_BIASED] << " " << STAT_label[STATL_HISTOGRAM];
  }
  os << "\t" << SEQ_label[SEQL_BACKWARD] << " " << SEQ_label[SEQL_RECURRENCE_TIME]
     << " " << STAT_label[STATL_HISTOGRAM] << "\t" << SEQ_label[SEQL_FORWARD]
     << " " << SEQ_label[SEQL_RECURRENCE_TIME] << " " << STAT_label[STATL_HISTOGRAM];

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
  Time_events::spreadsheet_write(os);

  // ecriture des probabilites de non-evenement/evenement fonction du temps

  os << "\n\t" << SEQ_label[SEQL_OBSERVED] << " " << SEQ_label[SEQL_NO_EVENT_PROBABILITY]
     << "\t" << SEQ_label[SEQL_OBSERVED] << " " << SEQ_label[SEQL_EVENT_PROBABILITY]
     << "\t" << STAT_label[STATL_FREQUENCY] << endl;

  index_event->spreadsheet_print(os);

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Renewal_data dans un fichier au format tableur.
 *
 *  arguments : reference sur un objet Format_error, path.
 *
 *--------------------------------------------------------------*/

bool Renewal_data::spreadsheet_write(Format_error &error , const char *path) const

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
 *  Sortie Gnuplot d'un objet Renewal_data.
 *
 *  arguments : reference sur un objet Format_error, prefixe des fichiers,
 *              titre des figures.
 *
 *--------------------------------------------------------------*/

bool Renewal_data::plot_write(Format_error &error , const char *prefix ,
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
    const Histogram **phisto;
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

    phisto = new const Histogram*[nb_histo];

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
      phisto[nb_histo++] = hmixture;
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
                   << STAT_label[STATL_HISTOGRAM] << "\" with impulses" << endl;

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
                   << STAT_label[STATL_HISTOGRAM] << "\" with impulses" << endl;

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
                   << STAT_label[STATL_HISTOGRAM] << "\" with impulses" << endl;

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
                 << SEQ_label[SEQL_RECURRENCE_TIME] << " " << STAT_label[STATL_HISTOGRAM]
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
                 << SEQ_label[SEQL_RECURRENCE_TIME] << " " << STAT_label[STATL_HISTOGRAM]
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
                   << STAT_label[STATL_HISTOGRAM] << "\" with impulses" << endl;

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
                     << STAT_label[STATL_HISTOGRAM] << "\" with impulses" << endl;

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

          if (hmixture->nb_value - 1 < TIC_THRESHOLD) {
            out_file << "set xtics 0,1" << endl;
          }
          if ((int)(hmixture->max * YSCALE) + 1 < TIC_THRESHOLD) {
            out_file << "set ytics 0,1" << endl;
          }

          out_file << "plot [0:" << hmixture->nb_value - 1 << "] [0:"
                   << (int)(hmixture->max * YSCALE) + 1 << "] \""
                   << label((data_file_name[0].str()).c_str()) << "\" using " << j
                   << " title \"" << SEQ_label[SEQL_NB_EVENT] << " "
                   << STAT_label[STATL_HISTOGRAM] << "\" with impulses" << endl;

          if (hmixture->nb_value - 1 < TIC_THRESHOLD) {
            out_file << "set xtics autofreq" << endl;
          }
          if ((int)(hmixture->max * YSCALE) + 1 < TIC_THRESHOLD) {
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
 *  Fonctions pour la persistance.
 *
 *--------------------------------------------------------------*/

/* RWDEFINE_COLLECTABLE(Renewal_data , STATI_RENEWAL_DATA);


RWspace Renewal_data::binaryStoreSize() const

{
  register int i;
  RWspace size;


  size = Time_events::binaryStoreSize() + sizeof(type) + sizeof(*length) * nb_element;

  for (i = 0;i < nb_element;i++) {
    size += sizeof(**sequence) * length[i];
  }

  if (inter_event) {
    size += inter_event->binaryStoreSize();
  }
  size += 2 * sizeof(true) + within->binaryStoreSize();
  if (length_bias) {
    size += length_bias->binaryStoreSize();
  }
  size += backward->binaryStoreSize() + forward->binaryStoreSize() + index_event->binaryStoreSize();

  if (renewal) {
    size += renewal->recursiveStoreSize();
  }

  return size;
}


void Renewal_data::restoreGuts(RWvistream &is)

{
  bool status;
  register int i , j;
  int *psequence;


  remove();

  Time_events::restoreGuts(is);

  is >> type;

  length = new int[nb_element];
  for (i = 0;i < nb_element;i++) {
    is >> length[i];
  }

  sequence = new int*[nb_element];
  for (i = 0;i < nb_element;i++) {
    sequence[i] = new int[length[i]];
    psequence = sequence[i];
    for (j = 0;j < length[i];j++) {
      is >> *psequence++;
    }
  }

  is >> status;
  if (status) {
    inter_event = new Histogram();
    inter_event->restoreGuts(is);
  }
  else {
    inter_event = 0;
  }

  within = new Histogram();
  within->restoreGuts(is);

  is >> status;
  if (status) {
    length_bias = new Histogram();
    length_bias->restoreGuts(is);
  }
  else {
    length_bias = 0;
  }

  backward = new Histogram();
  backward->restoreGuts(is);
  forward = new Histogram();
  forward->restoreGuts(is);

  index_event = new Curves();
  index_event->restoreGuts(is);

  is >> renewal;
  if (renewal == RWnilCollectable) {
    renewal = 0;
  }
}


void Renewal_data::restoreGuts(RWFile &file)

{
  bool status;
  register int i;


  remove();

  Time_events::restoreGuts(file);

  file.Read(type);

  length = new int[nb_element];
  file.Read(length , nb_element);

  sequence = new int*[nb_element];
  for (i = 0;i < nb_element;i++) {
    sequence[i] = new int[length[i]];
    file.Read(sequence[i] , length[i]);
  }

  file.Read(status);
  if (status) {
    inter_event = new Histogram();
    inter_event->restoreGuts(file);
  }
  else {
    inter_event = 0;
  }

  within = new Histogram();
  within->restoreGuts(file);

  file.Read(status);
  if (status) {
    length_bias = new Histogram();
    length_bias->restoreGuts(file);
  }
  else {
    length_bias = 0;
  }

  backward = new Histogram();
  backward->restoreGuts(file);
  forward = new Histogram();
  forward->restoreGuts(file);

  index_event = new Curves();
  index_event->restoreGuts(file);

  file >> renewal;
  if (renewal == RWnilCollectable) {
    renewal = 0;
  }
}


void Renewal_data::saveGuts(RWvostream &os) const

{
  register int i , j;
  int *psequence;


  Time_events::saveGuts(os);

  os << type;

  for (i = 0;i < nb_element;i++) {
    os << length[i];
  }

  for (i = 0;i < nb_element;i++) {
    psequence = sequence[i];
    for (j = 0;j < length[i];j++) {
      os << *psequence++;
    }
  }

  if (inter_event) {
    os << true;
    inter_event->saveGuts(os);
  }
  else {
    os << false;
  }

  within->saveGuts(os);

  if (length_bias) {
    os << true;
    length_bias->saveGuts(os);
  }
  else {
    os << false;
  }

  backward->saveGuts(os);
  forward->saveGuts(os);

  index_event->saveGuts(os);

  if (renewal) {
    os << renewal;
  }
  else {
    os << RWnilCollectable;
  }
}


void Renewal_data::saveGuts(RWFile &file) const

{
  register int i;


  Time_events::saveGuts(file);

  file.Write(type);

  file.Write(length , nb_element);

  for (i = 0;i < nb_element;i++) {
    file.Write(sequence[i] , length[i]);
  }

  if (inter_event) {
    file.Write(true);
    inter_event->saveGuts(file);
  }
  else {
    file.Write(false);
  }

  within->saveGuts(file);

  if (length_bias) {
    file.Write(true);
    length_bias->saveGuts(file);
  }
  else {
    file.Write(false);
  }

  backward->saveGuts(file);
  forward->saveGuts(file);

  index_event->saveGuts(file);

  if (renewal) {
    file << renewal;
  }
  else {
    file << RWnilCollectable;
  }
} */


/*--------------------------------------------------------------*
 *
 *  Extraction des probabilites de non-evenement/evenement fonction du temps.
 *
 *  argument : decalage (0/1).
 *
 *--------------------------------------------------------------*/

void Renewal_data::build_index_event(int offset)

{
  register int i , j;
  int frequency[2];


  index_event = new Curves(2 , htime->nb_value , true);
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
