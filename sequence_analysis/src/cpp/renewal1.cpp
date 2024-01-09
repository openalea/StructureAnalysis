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
 *       $Id: renewal1.cpp 18063 2015-04-23 10:50:05Z guedon $
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

using namespace std;
using namespace stat_tool;


namespace sequence_analysis {



/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe NbEvent.
 *
 *  arguments : type, temps d'observation, nombre de valeurs,
 *              identificateur et parametres de la loi.
 *
 *--------------------------------------------------------------*/

NbEvent::NbEvent(char itype , int itime , int inb_value , int iident ,
                 int iinf_bound , int isup_bound , double iparameter , double iprobability)
:DiscreteParametric(inb_value , iident , iinf_bound , isup_bound , iparameter , iprobability)

{
  type = itype;
  time = itime;
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe NbEvent.
 *
 *  arguments : type, temps d'observation, reference sur un objet DiscreteParametric.
 *
 *--------------------------------------------------------------*/

NbEvent::NbEvent(char itype , int itime , DiscreteParametric &inter_event)

{
  type = itype;
  time = itime;

  switch (type) {
  case 'o' :
    Distribution::init(time / inter_event.offset + 1);
    break;
  case 'e' :
    Distribution::init((time - 1) / inter_event.offset + 2);
    break;
  }

  ident = inter_event.ident;

  inf_bound = inter_event.inf_bound;
  sup_bound = inter_event.sup_bound;
  parameter = inter_event.parameter;
  probability = inter_event.probability;

  switch (type) {
  case 'o' :
    ordinary_computation(inter_event);
    break;
  case 'e' :
    computation(inter_event);
    break;
  }
}


/*--------------------------------------------------------------*
 *
 *  Constructeur par copie de la classe NbEvent.
 *
 *  arguments : reference sur un objet NbEvent, nombre de valeurs allouees.
 *
 *--------------------------------------------------------------*/

NbEvent::NbEvent(const NbEvent &nb_event , int ialloc_nb_value)
:DiscreteParametric(nb_event , 'c' , ialloc_nb_value)

{
  type = nb_event.type;
  time = nb_event.time;
}


/*--------------------------------------------------------------*
 *
 *  Initialisation du type ('o' : ordinaire, 'e' : en equilibre)
 *  d'un processus de renouvellement.
 *
 *  argument : type.
 *
 *--------------------------------------------------------------*/

void Renewal::type_init(int itype)

{
  if (itype != type) {
    register int i;


    type = itype;

    for (i = time->offset;i < time->nb_value;i++) {
      if (time->mass[i] > 0.) {
        nb_event[i]->type = itype;
      }
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Initialisation des parametres des lois d'un processus de renouvellement.
 *
 *  arguments : parametres d'une loi.
 *
 *--------------------------------------------------------------*/

void Renewal::init(int inf_bound , int sup_bound , double parameter , double probability)

{
  register int i;


  inter_event->init(inf_bound , sup_bound , parameter , probability);
  length_bias->init(inf_bound , sup_bound , parameter , probability);
  backward->init(inf_bound , sup_bound , parameter , probability);
  forward->init(inf_bound , sup_bound , parameter , probability);

  for (i = time->offset;i < time->nb_value;i++) {
    if (time->mass[i] > 0.) {
      nb_event[i]->init(inf_bound , sup_bound , parameter , probability);
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Initialisation de l'identificateur et des parametres des lois
 *  d'un processus de renouvellement.
 *
 *  arguments : identificateur et parametres d'une loi.
 *
 *--------------------------------------------------------------*/

void Renewal::init(int ident , int inf_bound , int sup_bound ,
                   double parameter , double probability)

{
  register int i;


  inter_event->init(ident , inf_bound , sup_bound , parameter , probability);
  length_bias->init(ident , inf_bound , sup_bound , parameter , probability);
  backward->init(ident , inf_bound , sup_bound , parameter , probability);
  forward->init(ident , inf_bound , sup_bound , parameter , probability);

  for (i = time->offset;i < time->nb_value;i++) {
    if (time->mass[i] > 0.) {
      nb_event[i]->init(ident , inf_bound , sup_bound , parameter , probability);
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Constructeur par defaut de la classe Renewal.
 *
 *--------------------------------------------------------------*/

Renewal::Renewal()

{
  nb_iterator = 0;
  renewal_data = NULL;

  type = 'v';

  time = NULL;

  inter_event = NULL;
  length_bias = NULL;
  backward = NULL;
  forward = NULL;

  nb_event_max = 0;
  nevent_time = NULL;

  nb_event = NULL;
  mixture = NULL;

  index_event = NULL;
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Renewal.
 *
 *  arguments : type, reference sur la loi empirique du temps d'observation et
 *              sur la loi inter-evenement.
 *
 *--------------------------------------------------------------*/

Renewal::Renewal(char itype , const FrequencyDistribution &htime ,
                 const DiscreteParametric &iinter_event)

{
  register int i;
  int nb_value;


  nb_iterator = 0;
  renewal_data = NULL;

  type = itype;

  time = new Distribution(htime);

  inter_event = new DiscreteParametric(iinter_event , 'c' , (int)(iinter_event.nb_value * NB_VALUE_COEFF));
  length_bias = new LengthBias(inter_event->alloc_nb_value , inter_event->ident ,
                               inter_event->inf_bound , inter_event->sup_bound ,
                               inter_event->parameter , inter_event->probability);
  backward = new Backward(inter_event->alloc_nb_value - 1 , inter_event->ident ,
                          inter_event->inf_bound , inter_event->sup_bound ,
                          inter_event->parameter , inter_event->probability);
  forward = new Forward(inter_event->alloc_nb_value , inter_event->ident ,
                        inter_event->inf_bound , inter_event->sup_bound ,
                        inter_event->parameter , inter_event->probability);

  nb_event = new NbEvent*[time->nb_value];

  for (i = 0;i < time->offset;i++) {
    nb_event[i] = NULL;
  }
  for (i = time->offset;i < time->nb_value;i++) {
    if (time->mass[i] > 0.) {
      switch (type) {
      case 'o' :
        nb_value = i / inter_event->offset + 1;
        break;
      case 'e' :
        nb_value = (i - 1) / inter_event->offset + 2;
        break;
      }

      nb_event[i] = new NbEvent(type , i , nb_value , inter_event->ident ,
                                inter_event->inf_bound , inter_event->sup_bound ,
                                inter_event->parameter , inter_event->probability);
    }

    else {
      nb_event[i] = NULL;
    }
  }

  mixture = new Distribution(nb_value);

  nb_event_max = nb_value - 1;
  nevent_time = new DiscreteParametric*[nb_event_max + 1];
  for (i = 0;i <= nb_event_max;i++) {
    nevent_time[i] = NULL;
  }

  index_event = new Curves(2 , time->nb_value);
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Renewal.
 *
 *  arguments : type, references sur la loi du temps d'observation et
 *              sur la loi inter-evenement.
 *
 *--------------------------------------------------------------*/

Renewal::Renewal(char itype , const Distribution &itime ,
                 const DiscreteParametric &iinter_event)

{
  register int i;
  int nb_value;


  nb_iterator = 0;
  renewal_data = NULL;

  type = itype;

  time = new Distribution(itime);

  inter_event = new DiscreteParametric(iinter_event , 'n');
  length_bias = new LengthBias(*inter_event);
  backward = new Backward(inter_event->nb_value - 1 , inter_event->ident ,
                          inter_event->inf_bound , inter_event->sup_bound ,
                          inter_event->parameter , inter_event->probability);
  forward = new Forward(*inter_event);

  nb_event = new NbEvent*[time->nb_value];

  for (i = 0;i < time->offset;i++) {
    nb_event[i] = NULL;
  }
  for (i = time->offset;i < time->nb_value;i++) {
    if (time->mass[i] > 0.) {
      switch (type) {
      case 'o' :
        nb_value = i / inter_event->offset + 1;
        break;
      case 'e' :
        nb_value = (i - 1) / inter_event->offset + 2;
        break;
      }

      nb_event[i] = new NbEvent(type , i , nb_value , inter_event->ident ,
                                inter_event->inf_bound , inter_event->sup_bound ,
                                inter_event->parameter , inter_event->probability);
    }

    else {
      nb_event[i] = NULL;
    }
  }

  mixture = new Distribution(nb_value);

  nb_event_max = nb_value - 1;
  nevent_time = new DiscreteParametric*[nb_event_max + 1];
  for (i = 0;i <= nb_event_max;i++) {
    nevent_time[i] = NULL;
  }

  index_event = new Curves(2 , time->nb_value);

  computation(false);
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Renewal.
 *
 *  arguments : references sur un objet RenewalData et sur la loi inter-evenement.
 *
 *--------------------------------------------------------------*/

Renewal::Renewal(const RenewalData &irenewal_data , const DiscreteParametric &iinter_event)

{
  register int i;
  int nb_value;


  nb_iterator = 0;
  renewal_data = new RenewalData(irenewal_data);
  renewal_data->type = 'e';

  type = 'e';

  time = new Distribution(*(renewal_data->htime));

  inter_event = new DiscreteParametric(iinter_event);
  length_bias = new LengthBias(*inter_event);
  backward = new Backward(inter_event->nb_value - 1 , inter_event->ident ,
                          inter_event->inf_bound , inter_event->sup_bound ,
                          inter_event->parameter , inter_event->probability);
  forward = new Forward(*inter_event);

  nb_event = new NbEvent*[time->nb_value];

  for (i = 0;i < time->offset;i++) {
    nb_event[i] = NULL;
  }
  for (i = time->offset;i < time->nb_value;i++) {
    if (time->mass[i] > 0.) {
      switch (type) {
      case 'o' :
        nb_value = i / inter_event->offset + 1;
        break;
      case 'e' :
        nb_value = (i - 1) / inter_event->offset + 2;
        break;
      }

      nb_event[i] = new NbEvent(type , i , nb_value , inter_event->ident ,
                                inter_event->inf_bound , inter_event->sup_bound ,
                                inter_event->parameter , inter_event->probability);
    }

    else {
      nb_event[i] = NULL;
    }
  }

  mixture = new Distribution(nb_value);

  nb_event_max = nb_value - 1;
  nevent_time = new DiscreteParametric*[nb_event_max + 1];
  for (i = 0;i <= nb_event_max;i++) {
    nevent_time[i] = NULL;
  }

  index_event = new Curves(2 , time->nb_value);

  computation(false);
}


/*--------------------------------------------------------------*
 *
 *  Copie d'un objet Renewal.
 *
 *  arguments : reference sur un objet Renewal.
 *              flag copie de l'objet renewal_data.
 *
 *--------------------------------------------------------------*/

void Renewal::copy(const Renewal &renew , bool data_flag)

{
  register int i;


  nb_iterator = 0;

  if ((data_flag) && (renew.renewal_data)) {
    renewal_data = new RenewalData(*(renew.renewal_data));
  }
  else {
    renewal_data = NULL;
  }

  type = renew.type;

  time = new Distribution(*(renew.time));

  inter_event = new DiscreteParametric(*(renew.inter_event));
  length_bias = new LengthBias(*(renew.length_bias));
  backward = new Backward(*(renew.backward) , renew.backward->alloc_nb_value);
  forward = new Forward(*(renew.forward));

  nb_event_max = renew.mixture->nb_value - 1;

  nevent_time = new DiscreteParametric*[nb_event_max + 1];

  nevent_time[0] = NULL;
  for (i = 1;i <= nb_event_max;i++) {
    if (renew.nevent_time[i]) {
      nevent_time[i] = new DiscreteParametric(*(renew.nevent_time[i]));
    }
    else {
      nevent_time[i] = NULL;
    }
  }

  nb_event = new NbEvent*[time->nb_value];

  for (i = 0;i < time->offset;i++) {
    nb_event[i] = NULL;
  }
  for (i = time->offset;i < time->nb_value;i++) {
    if (time->mass[i] > 0.) {
      nb_event[i] = new NbEvent(*(renew.nb_event[i]));
    }
    else {
      nb_event[i] = NULL;
    }
  }

  mixture = new Distribution(*(renew.mixture));

  index_event = new Curves(*(renew.index_event));
}


/*--------------------------------------------------------------*
 *
 *  Destruction des champs d'un objet Renewal.
 *
 *--------------------------------------------------------------*/

void Renewal::remove()

{
  register int i;


  delete renewal_data;

  delete inter_event;
  delete length_bias;
  delete backward;
  delete forward;

  if (nevent_time) {
    for (i = 1;i <= nb_event_max;i++) {
      delete nevent_time[i];
    }
    delete [] nevent_time;
  }

  if (nb_event) {
    for (i = time->offset;i < time->nb_value;i++) {
      if (time->mass[i] > 0.) {
        delete nb_event[i];
      }
    }
    delete [] nb_event;
  }

  delete mixture;

  delete index_event;

  delete time;
}


/*--------------------------------------------------------------*
 *
 *  Destructeur de la classe Renewal.
 *
 *--------------------------------------------------------------*/

Renewal::~Renewal()

{
  remove();
}


/*--------------------------------------------------------------*
 *
 *  Destruction d'un objet Renewal en tenant compte du nombre
 *  d'iterateurs pointant dessus.
 *
 *--------------------------------------------------------------*/

void Renewal::conditional_delete()

{
  if (nb_iterator == 0) {
    delete this;
  }
}


/*--------------------------------------------------------------*
 *
 *  Operateur d'assignement de la classe Renewal.
 *
 *  argument : reference sur un objet Renewal.
 *
 *--------------------------------------------------------------*/

Renewal& Renewal::operator=(const Renewal &renew)

{
  if ((&renew != this) && (nb_iterator == 0)) {
    remove();
    copy(renew);
  }

  return *this;
}


/*--------------------------------------------------------------*
 *
 *  Extraction d'une loi.
 *
 *  arguments : reference sur un objet StatError, type de loi,
 *              temps d'observation.
 *
 *--------------------------------------------------------------*/

DiscreteParametricModel* Renewal::extract(StatError &error , int dist_type , int itime) const

{
  Distribution *pdist;
  DiscreteParametric *pparam;
  DiscreteParametricModel *dist;
  FrequencyDistribution *phisto;


  if (dist_type == NB_EVENT) {
    error.init();

    if ((itime < time->offset) || (itime >= time->nb_value) || (time->mass[itime] == 0.)) {
      dist = NULL;
      error.update(SEQ_error[SEQR_OBSERVATION_TIME]);
    }
    else {
      dist = new DiscreteParametricModel(*((Distribution*)nb_event[itime]) ,
                                         (renewal_data ? renewal_data->hnb_event[itime] : 0));
    }
  }

  else {
    pdist = NULL;
    pparam = NULL;

    switch (dist_type) {
    case INTER_EVENT :
      pparam = inter_event;
      break;
    case LENGTH_BIAS :
      pdist = length_bias;
      break;
    case BACKWARD_RECURRENCE_TIME :
      pdist = backward;
      break;
    case FORWARD_RECURRENCE_TIME :
      pdist = forward;
      break;
    case NB_EVENT_MIXTURE :
      pdist = mixture;
      break;
    }

    phisto = NULL;
    if (renewal_data) {
      switch (dist_type) {

      case INTER_EVENT : {
        if (renewal_data->inter_event) {
          phisto = renewal_data->inter_event;
        }
        break;
      }

      case LENGTH_BIAS : {
        if (renewal_data->length_bias) {
          phisto = renewal_data->length_bias;
        }
        break;
      }

      case BACKWARD_RECURRENCE_TIME : {
        phisto = renewal_data->backward;
        break;
      }

      case FORWARD_RECURRENCE_TIME : {
        phisto = renewal_data->forward;
        break;
      }

      case NB_EVENT_MIXTURE : {
        phisto = renewal_data->mixture;
        break;
      }
      }

      if ((phisto) && (phisto->nb_element == 0)) {
        phisto = NULL;
      }
    }

    if (pdist) {
      dist = new DiscreteParametricModel(*pdist , phisto);
    }
    else if (pparam) {
      dist = new DiscreteParametricModel(*pparam , phisto);
    }
  }

  return dist;
}


/*--------------------------------------------------------------*
 *
 *  Construction d'un objet Renewal a partir d'une loi inter-evenement.
 *
 *  arguments : reference sur un objet StatError, reference sur la loi inter-evenement,
 *              type du processus ('o' : ordinaire, 'e' : en equilibre),
 *              temps d'observation.
 *
 *--------------------------------------------------------------*/

Renewal* renewal_building(StatError &error , const DiscreteParametric &inter_event ,
                          char type , int time)

{
  bool status = true;
  Renewal *renew;


  renew = NULL;
  error.init();

  if (inter_event.offset == 0) {
    status = false;
    error.update(STAT_error[STATR_MIN_VALUE]);
  }
  if (time < MAX(inter_event.offset , 2)) {
    status = false;
    error.update(SEQ_error[SEQR_SHORT_OBSERVATION_TIME]);
  }
  if (time > MAX_TIME) {
    status = false;
    error.update(SEQ_error[SEQR_LONG_OBSERVATION_TIME]);
  }

  if (status) {
    DiscreteParametric dtime(UNIFORM , time , time , D_DEFAULT , D_DEFAULT);
    renew = new Renewal(type , dtime , inter_event);
  }

  return renew;
}


/*--------------------------------------------------------------*
 *
 *  Construction d'un objet Renewal a partir d'un fichier.
 *
 *  arguments : reference sur un objet StatError, path,
 *              type du processus ('o' : ordinaire, 'e' : en equilibre),
 *              temps d'observation, seuil sur la fonction de repartition
 *              de la loi inter-evenement.
 *
 *--------------------------------------------------------------*/

Renewal* renewal_ascii_read(StatError &error , const char *path ,
                            char type , int time , double cumul_threshold)

{
  RWLocaleSnapshot locale("en");
  RWCString buffer , token;
  size_t position;
  bool status;
  int line;
  DiscreteParametric *inter_event;
  Renewal *renew;
  ifstream in_file(path);


  renew = NULL;
  error.init();

  if (!in_file) {
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {
    status = true;
    line = 0;

    inter_event = discrete_parametric_parsing(error , in_file , line ,
                                              NEGATIVE_BINOMIAL , cumul_threshold , 1);

    if (!inter_event) {
      status = false;
    }

    else {
      if (time < MAX(inter_event->offset , 2)) {
        status = false;
        error.update(SEQ_error[SEQR_SHORT_OBSERVATION_TIME]);
      }
      if (time > MAX_TIME) {
        status = false;
        error.update(SEQ_error[SEQR_LONG_OBSERVATION_TIME]);
      }
    }

    while (buffer.readLine(in_file , false)) {
      line++;

#     ifdef DEBUG
      cout << line << " " << buffer << endl;
#     endif

      position = buffer.first('#');
      if (position != RW_NPOS) {
        buffer.remove(position);
      }
      if (!(buffer.isNull())) {
        status = false;
        error.update(STAT_parsing[STATP_FORMAT] , line);
      }
    }

    if (status) {
      DiscreteParametric dtime(UNIFORM , time , time , D_DEFAULT , D_DEFAULT);
      renew = new Renewal(type , dtime , *inter_event);
    }

    delete inter_event;
  }

  return renew;
}


};  // namespace sequence_analysis
