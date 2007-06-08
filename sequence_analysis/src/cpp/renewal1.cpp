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

using namespace std;



/*--------------------------------------------------------------*
 *
 *  Construction d'une loi inter-evenement a partir d'une loi
 *  inter-evenement initiale par dilatation/retraction de l'echelle des temps.
 *
 *  arguments : reference sur une loi inter-evenement, facteur d'echelle.
 *
 *--------------------------------------------------------------*/

Parametric::Parametric(const Parametric &dist , double scale)
:Distribution((int)(dist.nb_value * scale) + 1)

{
  double bmean , bvariance , ratio , shift_mean;


  bmean = scale * dist.mean;
  bvariance = scale * scale * dist.variance;

  inf_bound = (int)(dist.inf_bound * scale);
  sup_bound = I_DEFAULT;
  parameter = D_DEFAULT;
  probability = D_DEFAULT;

  shift_mean = bmean - inf_bound;
  ratio = bvariance / shift_mean;

  // cas binomiale

  if (ratio < 1. - POISSON_RANGE) {
    ident = BINOMIAL;
    sup_bound = (int)ceil(inf_bound + shift_mean * shift_mean / (shift_mean - bvariance));
    if (sup_bound <= inf_bound) {
      sup_bound = inf_bound + 1;
    }
    probability = shift_mean / (sup_bound - inf_bound);
  }

  // cas binomiale negative

  else if (ratio > 1. + POISSON_RANGE) {
    ident = NEGATIVE_BINOMIAL;
    parameter = shift_mean * shift_mean / (bvariance - shift_mean);
    probability = shift_mean / bvariance;
  }

  // cas Poisson

  else {
    ident = POISSON;
    parameter = shift_mean;
  }

  computation();
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Nb_event.
 *
 *  arguments : type, temps d'observation, nombre de valeurs,
 *              identificateur et parametres de la loi.
 *
 *--------------------------------------------------------------*/

Nb_event::Nb_event(char itype , int itime , int inb_value , int iident ,
                   int iinf_bound , int isup_bound , double iparameter , double iprobability)
:Parametric(inb_value , iident , iinf_bound , isup_bound , iparameter , iprobability)

{
  type = itype;
  time = itime;
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Nb_event.
 *
 *  arguments : type, temps d'observation, reference sur un objet Parametric.
 *
 *--------------------------------------------------------------*/

Nb_event::Nb_event(char itype , int itime , Parametric &inter_event)

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
 *  Constructeur par copie de la classe Nb_event.
 *
 *  arguments : reference sur un objet Nb_event, nombre de valeurs allouees.
 *
 *--------------------------------------------------------------*/

Nb_event::Nb_event(const Nb_event &nb_event , int ialloc_nb_value)
:Parametric(nb_event , 'c' , ialloc_nb_value)

{
  type = nb_event.type;
  time = nb_event.time;
}


/*--------------------------------------------------------------*
 *
 *  Fonctions pour la persistance.
 *
 *--------------------------------------------------------------*/

/* RWspace Nb_event::binaryStoreSize(int ialloc_nb_value) const

{
  RWspace size = Parametric::binaryStoreSize(ialloc_nb_value) +
                 sizeof(type) + sizeof(time);

  return size;
}


void Nb_event::restoreGuts(RWvistream &is)

{
  Parametric::restoreGuts(is);

  is >> type >> time;
}


void Nb_event::restoreGuts(RWFile &file)

{
  Parametric::restoreGuts(file);

  file.Read(type);
  file.Read(time);
}


void Nb_event::saveGuts(RWvostream &os , int ialloc_nb_value) const

{
  Parametric::saveGuts(os , ialloc_nb_value);

  os << type << time;
}


void Nb_event::saveGuts(RWFile &file , int ialloc_nb_value) const

{
  Parametric::saveGuts(file , ialloc_nb_value);

  file.Write(type);
  file.Write(time);
} */


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
  renewal_data = 0;

  type = 'v';

  time = 0;

  inter_event = 0;
  length_bias = 0;
  backward = 0;
  forward = 0;

  nb_event_max = 0;
  nevent_time = 0;

  nb_event = 0;
  mixture = 0;

  index_event = 0;
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Renewal.
 *
 *  arguments : type, reference sur l'histogramme du temps d'observation et
 *              sur la loi inter-evenement.
 *
 *--------------------------------------------------------------*/

Renewal::Renewal(char itype , const Histogram &htime , const Parametric &iinter_event)

{
  register int i;
  int nb_value;


  nb_iterator = 0;
  renewal_data = 0;

  type = itype;

  time = new Distribution(htime);

  inter_event = new Parametric(iinter_event , 'c' , (int)(iinter_event.nb_value * NB_VALUE_COEFF));
  length_bias = new Length_bias(inter_event->alloc_nb_value , inter_event->ident ,
                                inter_event->inf_bound , inter_event->sup_bound ,
                                inter_event->parameter , inter_event->probability);
  backward = new Backward(inter_event->alloc_nb_value - 1 , inter_event->ident ,
                          inter_event->inf_bound , inter_event->sup_bound ,
                          inter_event->parameter , inter_event->probability);
  forward = new Forward(inter_event->alloc_nb_value , inter_event->ident ,
                        inter_event->inf_bound , inter_event->sup_bound ,
                        inter_event->parameter , inter_event->probability);

  nb_event = new Nb_event*[time->nb_value];

  for (i = 0;i < time->offset;i++) {
    nb_event[i] = 0;
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

      nb_event[i] = new Nb_event(type , i , nb_value , inter_event->ident ,
                                 inter_event->inf_bound , inter_event->sup_bound ,
                                 inter_event->parameter , inter_event->probability);
    }

    else {
      nb_event[i] = 0;
    }
  }

  mixture = new Distribution(nb_value);

  nb_event_max = nb_value - 1;
  nevent_time = new Parametric*[nb_event_max + 1];
  for (i = 0;i <= nb_event_max;i++) {
    nevent_time[i] = 0;
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
                 const Parametric &iinter_event)

{
  register int i;
  int nb_value;


  nb_iterator = 0;
  renewal_data = 0;

  type = itype;

  time = new Distribution(itime);

  inter_event = new Parametric(iinter_event , 'n');
  length_bias = new Length_bias(*inter_event);
  backward = new Backward(inter_event->nb_value - 1 , inter_event->ident ,
                          inter_event->inf_bound , inter_event->sup_bound ,
                          inter_event->parameter , inter_event->probability);
  forward = new Forward(*inter_event);

  nb_event = new Nb_event*[time->nb_value];

  for (i = 0;i < time->offset;i++) {
    nb_event[i] = 0;
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

      nb_event[i] = new Nb_event(type , i , nb_value , inter_event->ident ,
                                 inter_event->inf_bound , inter_event->sup_bound ,
                                 inter_event->parameter , inter_event->probability);
    }

    else {
      nb_event[i] = 0;
    }
  }

  mixture = new Distribution(nb_value);

  nb_event_max = nb_value - 1;
  nevent_time = new Parametric*[nb_event_max + 1];
  for (i = 0;i <= nb_event_max;i++) {
    nevent_time[i] = 0;
  }

  index_event = new Curves(2 , time->nb_value);

  computation(false);
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Renewal.
 *
 *  arguments : references sur un objet Renewal_data et sur la loi inter-evenement.
 *
 *--------------------------------------------------------------*/

Renewal::Renewal(const Renewal_data &irenewal_data , const Parametric &iinter_event)

{
  register int i;
  int nb_value;


  nb_iterator = 0;
  renewal_data = new Renewal_data(irenewal_data);
  renewal_data->type = 'e';

  type = 'e';

  time = new Distribution(*(renewal_data->htime));

  inter_event = new Parametric(iinter_event);
  length_bias = new Length_bias(*inter_event);
  backward = new Backward(inter_event->nb_value - 1 , inter_event->ident ,
                          inter_event->inf_bound , inter_event->sup_bound ,
                          inter_event->parameter , inter_event->probability);
  forward = new Forward(*inter_event);

  nb_event = new Nb_event*[time->nb_value];

  for (i = 0;i < time->offset;i++) {
    nb_event[i] = 0;
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

      nb_event[i] = new Nb_event(type , i , nb_value , inter_event->ident ,
                                 inter_event->inf_bound , inter_event->sup_bound ,
                                 inter_event->parameter , inter_event->probability);
    }

    else {
      nb_event[i] = 0;
    }
  }

  mixture = new Distribution(nb_value);

  nb_event_max = nb_value - 1;
  nevent_time = new Parametric*[nb_event_max + 1];
  for (i = 0;i <= nb_event_max;i++) {
    nevent_time[i] = 0;
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
    renewal_data = new Renewal_data(*(renew.renewal_data));
  }
  else {
    renewal_data = 0;
  }

  type = renew.type;

  time = new Distribution(*(renew.time));

  inter_event = new Parametric(*(renew.inter_event));
  length_bias = new Length_bias(*(renew.length_bias));
  backward = new Backward(*(renew.backward) , renew.backward->alloc_nb_value);
  forward = new Forward(*(renew.forward));

  nb_event_max = renew.mixture->nb_value - 1;

  nevent_time = new Parametric*[nb_event_max + 1];

  nevent_time[0] = 0;
  for (i = 1;i <= nb_event_max;i++) {
    if (renew.nevent_time[i]) {
      nevent_time[i] = new Parametric(*(renew.nevent_time[i]));
    }
    else {
      nevent_time[i] = 0;
    }
  }

  nb_event = new Nb_event*[time->nb_value];

  for (i = 0;i < time->offset;i++) {
    nb_event[i] = 0;
  }
  for (i = time->offset;i < time->nb_value;i++) {
    if (time->mass[i] > 0.) {
      nb_event[i] = new Nb_event(*(renew.nb_event[i]));
    }
    else {
      nb_event[i] = 0;
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
 *  arguments : reference sur un objet Format_error, type de loi,
 *              temps d'observation.
 *
 *--------------------------------------------------------------*/

Parametric_model* Renewal::extract(Format_error &error , int dist_type , int itime) const

{
  Distribution *pdist;
  Parametric *pparam;
  Parametric_model *dist;
  Histogram *phisto;


  if (dist_type == NB_EVENT) {
    error.init();

    if ((itime < time->offset) || (itime >= time->nb_value) || (time->mass[itime] == 0.)) {
      dist = 0;
      error.update(SEQ_error[SEQR_OBSERVATION_TIME]);
    }
    else {
      dist = new Parametric_model(*((Distribution*)nb_event[itime]) ,
                                  (renewal_data ? renewal_data->hnb_event[itime] : 0));
    }
  }

  else {
    pdist = 0;
    pparam = 0;

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
    case MIXTURE :
      pdist = mixture;
      break;
    }

    phisto = 0;
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

      case MIXTURE : {
        phisto = renewal_data->hmixture;
        break;
      }
      }

      if ((phisto) && (phisto->nb_element == 0)) {
        phisto = 0;
      }
    }

    if (pdist) {
      dist = new Parametric_model(*pdist , phisto);
    }
    else if (pparam) {
      dist = new Parametric_model(*pparam , phisto);
    }
  }

  return dist;
}


/*--------------------------------------------------------------*
 *
 *  Construction d'un objet Renewal a partir d'une loi inter-evenement.
 *
 *  arguments : reference sur un objet Format_error, reference sur la loi inter-evenement,
 *              type du processus ('o' : ordinaire, 'e' : en equilibre),
 *              temps d'observation.
 *
 *--------------------------------------------------------------*/

Renewal* renewal_building(Format_error &error , const Parametric &inter_event ,
                          char type , int time)

{
  bool status = true;
  Renewal *renew;


  renew = 0;
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
    Parametric dtime(UNIFORM , time , time , D_DEFAULT , D_DEFAULT);
    renew = new Renewal(type , dtime , inter_event);
  }

  return renew;
}


/*--------------------------------------------------------------*
 *
 *  Construction d'un objet Renewal a partir d'un fichier.
 *
 *  arguments : reference sur un objet Format_error, path,
 *              type du processus ('o' : ordinaire, 'e' : en equilibre),
 *              temps d'observation, seuil sur la fonction de repartition
 *              de la loi inter-evenement.
 *
 *--------------------------------------------------------------*/

Renewal* renewal_ascii_read(Format_error &error , const char *path ,
                            char type , int time , double cumul_threshold)

{
  RWLocaleSnapshot locale("en");
  RWCString buffer , token;
  size_t position;
  bool status;
  int line;
  Parametric *inter_event;
  Renewal *renew;
  ifstream in_file(path);


  renew = 0;
  error.init();

  if (!in_file) {
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {
    status = true;
    line = 0;

    inter_event = parametric_parsing(error , in_file , line ,
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
      Parametric dtime(UNIFORM , time , time , D_DEFAULT , D_DEFAULT);
      renew = new Renewal(type , dtime , *inter_event);
    }

    delete inter_event;
  }

  return renew;
}
