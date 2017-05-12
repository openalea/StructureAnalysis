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



#include <string>

#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/classification.hpp>

#include "stat_tool/stat_label.h"

#include "renewal.h"
#include "sequence_label.h"

using namespace std;
using namespace boost;
using namespace stat_tool;


namespace sequence_analysis {



/*--------------------------------------------------------------*/
/**
 *  \brief Constructor of the NbEvent class.
 *
 *  \param[in] itype        renewal process type (ORDINARY/EQUILIBRIUM),
 *  \param[in] itime        observation period,
 *  \param[in] inb_value    number of values,
 *  \param[in] iident       distribution identifier,
 *  \param[in] iinf_bound   lower bound,
 *  \param[in] isup_bound   upper bound (binomial, uniform),
 *  \param[in] iparameter   parameter (Poisson, negative binomial),
 *  \param[in] iprobability probability (binomial, negative binomial).
 */
/*--------------------------------------------------------------*/

NbEvent::NbEvent(process_type itype , int itime , int inb_value , discrete_parametric iident ,
                 int iinf_bound , int isup_bound , double iparameter , double iprobability)
:DiscreteParametric(inb_value , iident , iinf_bound , isup_bound , iparameter , iprobability)

{
  type = itype;
  time = itime;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Constructor of the NbEvent class.
 *
 *  \param[in] itype       renewal process type (ORDINARY/EQUILIBRIUM),
 *  \param[in] itime       observation period,
 *  \param[in] inter_event reference on a DiscreteParametric object.
 */
/*--------------------------------------------------------------*/

NbEvent::NbEvent(process_type itype , int itime , DiscreteParametric &inter_event)

{
  type = itype;
  time = itime;

  switch (type) {
  case ORDINARY :
    Distribution::init(time / inter_event.offset + 1);
    break;
  case EQUILIBRIUM :
    Distribution::init((time - 1) / inter_event.offset + 2);
    break;
  }

  ident = inter_event.ident;

  inf_bound = inter_event.inf_bound;
  sup_bound = inter_event.sup_bound;
  parameter = inter_event.parameter;
  probability = inter_event.probability;

  switch (type) {
  case ORDINARY :
    ordinary_computation(inter_event);
    break;
  case EQUILIBRIUM :
    computation(inter_event);
    break;
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Constructor by copy of the NbEvent class.
 *
 *  \param[in] nb_event        reference on a NbEvent object,
 *  \param[in] ialloc_nb_value number of allocated values.
 */
/*--------------------------------------------------------------*/

NbEvent::NbEvent(const NbEvent &nb_event , int ialloc_nb_value)
:DiscreteParametric(nb_event , DISTRIBUTION_COPY , ialloc_nb_value)

{
  type = nb_event.type;
  time = nb_event.time;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Initialization of the renewal process type (ORDINARY/EQUILIBRIUM).
 *
 *  \param[in] itype renewal process type.
 */
/*--------------------------------------------------------------*/

void Renewal::type_init(process_type itype)

{
  if (itype != type) {
    int i;


    type = itype;

    for (i = time->offset;i < time->nb_value;i++) {
      if (time->mass[i] > 0.) {
        nb_event[i]->type = itype;
      }
    }
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Initialization of the inter-event distribution parameters.
 *
 *  \param[in] inf_bound   lower bound,
 *  \param[in] sup_bound   upper bound (binomial, uniform),
 *  \param[in] parameter   parameter (Poisson, negative binomial),
 *  \param[in] probability probability (binomial, negative binomial).
 */
/*--------------------------------------------------------------*/

void Renewal::init(int inf_bound , int sup_bound , double parameter , double probability)

{
  int i;


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


/*--------------------------------------------------------------*/
/**
 *  \brief Initialization of the identifier and parameters of the inter-event distribution.
 *
 *  \param[in] ident       distribution identifier,
 *  \param[in] inf_bound   lower bound,
 *  \param[in] sup_bound   upper bound (binomial, uniform),
 *  \param[in] parameter   parameter (Poisson, negative binomial),
 *  \param[in] probability probability (binomial, negative binomial).
 */
/*--------------------------------------------------------------*/

void Renewal::init(discrete_parametric ident , int inf_bound , int sup_bound ,
                   double parameter , double probability)

{
  int i;


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


/*--------------------------------------------------------------*/
/**
 *  \brief Default constructor of the Renewal class.
 */
/*--------------------------------------------------------------*/

Renewal::Renewal()

{
  nb_iterator = 0;
  renewal_data = NULL;

  type = ORDINARY;

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


/*--------------------------------------------------------------*/
/**
 *  \brief Constructor of the Renewal class.
 *
 *  \param[in] itype        renewal process type (ORDINARY/EQUILIBRIUM),
 *  \param[in] htime        reference on the observation period frequency distribution,
 *  \param[in] iinter_event reference on on the inter-event distribution.
 */
/*--------------------------------------------------------------*/

Renewal::Renewal(process_type itype , const FrequencyDistribution &htime ,
                 const DiscreteParametric &iinter_event)

{
  int i;
  int nb_value;


  nb_iterator = 0;
  renewal_data = NULL;

  type = itype;

  time = new Distribution(htime);

  inter_event = new DiscreteParametric(iinter_event , DISTRIBUTION_COPY ,
                                       (int)(iinter_event.nb_value * NB_VALUE_COEFF));
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
      case ORDINARY :
        nb_value = i / inter_event->offset + 1;
        break;
      case EQUILIBRIUM :
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


/*--------------------------------------------------------------*/
/**
 *  \brief Constructor of the Renewal class.
 *
 *  \param[in] itype        renewal process type,
 *  \param[in] itime        reference on the observation distribution,
 *  \param[in] iinter_event reference on the inter-event distribution.
 */
/*--------------------------------------------------------------*/

Renewal::Renewal(process_type itype , const Distribution &itime ,
                 const DiscreteParametric &iinter_event)

{
  int i;
  int nb_value;


  nb_iterator = 0;
  renewal_data = NULL;

  type = itype;

  time = new Distribution(itime);

  inter_event = new DiscreteParametric(iinter_event , NORMALIZATION);
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
      case ORDINARY :
        nb_value = i / inter_event->offset + 1;
        break;
      case EQUILIBRIUM :
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


/*--------------------------------------------------------------*/
/**
 *  \brief Constructor of the Renewal class.
 *
 *  \param[in] irenewal_data reference on a RenewalData object,
 *  \param[in] iinter_event  reference on the inter-event distribution.
 */
/*--------------------------------------------------------------*/

Renewal::Renewal(const RenewalData &irenewal_data , const DiscreteParametric &iinter_event)

{
  int i;
  int nb_value;


  nb_iterator = 0;
  renewal_data = new RenewalData(irenewal_data);
  renewal_data->type = EQUILIBRIUM;

  type = EQUILIBRIUM;

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
      case ORDINARY :
        nb_value = i / inter_event->offset + 1;
        break;
      case EQUILIBRIUM :
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


/*--------------------------------------------------------------*/
/**
 *  \brief Copy of a Renewal object.
 *
 *  \param[in] renew     reference on a Renewal object,
 *  \param[in] data_flag flag copy of the included RenewalData object.
 */
/*--------------------------------------------------------------*/

void Renewal::copy(const Renewal &renew , bool data_flag)

{
  int i;


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


/*--------------------------------------------------------------*/
/**
 *  \brief Destruction of the data members of a Renewal object.
 */
/*--------------------------------------------------------------*/

void Renewal::remove()

{
  int i;


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


/*--------------------------------------------------------------*/
/**
 *  \brief Destructor of the Renewal class.
 */
/*--------------------------------------------------------------*/

Renewal::~Renewal()

{
  remove();
}


/*--------------------------------------------------------------*/
/**
 *  \brief Destruction of a Renewal object taking account of
 *         the number of iterators pointing to it.
 */
/*--------------------------------------------------------------*/

void Renewal::conditional_delete()

{
  if (nb_iterator == 0) {
    delete this;
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Assignment operator of the Renewal class.
 *
 *  \param[in] renew reference on a Renewal object.
 *
 *  \return          Renewal object.
 */
/*--------------------------------------------------------------*/

Renewal& Renewal::operator=(const Renewal &renew)

{
  if ((&renew != this) && (nb_iterator == 0)) {
    remove();
    copy(renew);
  }

  return *this;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Extraction of a distribution.
 *
 *  \param[in] error     reference on a StatError object,
 *  \param[in] dist_type distribution type,
 *  \param[in] itime     observation period.
 *
 *  \return              DiscreteParametricModel object.
 */
/*--------------------------------------------------------------*/

DiscreteParametricModel* Renewal::extract(StatError &error , renewal_distribution dist_type ,
                                          int itime) const

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


/*--------------------------------------------------------------*/
/**
 *  \brief Construction of a Renewal object on the basis of an inter-event distribution.
 *
 *  \param[in] error       reference on a StatError object,
 *  \param[in] inter_event reference on the inter-event distribution,
 *  \param[in] type        renewal process type (ORDINARY/EQUILIBRIUM),
 *  \param[in] time        observation period.
 *
 *  \return                Renewal object.
 */
/*--------------------------------------------------------------*/

Renewal* Renewal::building(StatError &error , const DiscreteParametric &inter_event ,
                           process_type type , int time)

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


/*--------------------------------------------------------------*/
/**
 *  \brief Construction of a Renewal object from a file.
 *
 *  \param[in] error           reference on a StatError object,
 *  \param[in] path            file path,
 *  \param[in] type            renewal process type (ORDINARY/EQUILIBRIUM),
 *  \param[in] time            observation period,
 *  \param[in] cumul_threshold threshold on the cumulative inter-event distribution function.
 *
 *  \return                    Renewal object.
 */
/*--------------------------------------------------------------*/

Renewal* Renewal::ascii_read(StatError &error , const string path ,
                             process_type type , int time , double cumul_threshold)

{
  string buffer;
  size_t position;
  bool status;
  int line;
  DiscreteParametric *inter_event;
  Renewal *renew;
  ifstream in_file(path.c_str());


  renew = NULL;
  error.init();

  if (!in_file) {
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {
    status = true;
    line = 0;

    inter_event = DiscreteParametric::parsing(error , in_file , line ,
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

    while (getline(in_file , buffer)) {
      line++;

#     ifdef DEBUG
      cout << line << " " << buffer << endl;
#     endif

      position = buffer.find('#');
      if (position != string::npos) {
        buffer.erase(position);
      }
      if (!(trim_right_copy_if(buffer , is_any_of(" \t")).empty())) {
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
