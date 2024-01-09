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
 *       $Id: renewal2.cpp 18064 2015-04-23 10:50:21Z guedon $
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
 *  Ecriture sur une ligne d'un objet Renewal.
 *
 *  argument : stream.
 *
 *--------------------------------------------------------------*/

ostream& Renewal::line_write(ostream &os) const

{
  switch (type) {
  case 'o' :
    os << SEQ_label[SEQL_ORDINARY_RENEWAL] << " - ";
    break;
  case 'e' :
    os << SEQ_label[SEQL_EQUILIBRIUM_RENEWAL] << " - ";
    break;
  }

  os << SEQ_label[SEQL_INTER_EVENT] << " " << STAT_label[STATL_DISTRIBUTION] << "   "
     << STAT_discrete_distribution_word[inter_event->ident] << "   "
     << STAT_label[STATL_MEAN] << ": " << inter_event->mean << "   "
     << STAT_label[STATL_VARIANCE] << ": " << inter_event->variance;

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un processus de renouvellement et de la structure de donnees associee.
 *
 *  arguments : stream, pointeur sur un objet RenewalData,
 *              flag niveau de detail, flag fichier.
 *
 *--------------------------------------------------------------*/

ostream& Renewal::ascii_write(ostream &os , const RenewalData *timev ,
                              bool exhaustive , bool file_flag) const

{
  register int i , j;
  int nb_dist , inf , sup;
  double likelihood , information , *scale;
  const Distribution **pdist;
  Test test(CHI2);


  if (file_flag) {
    os << "# ";
  }
  switch (type) {
  case 'o' :
    os << SEQ_label[SEQL_ORDINARY_RENEWAL] << endl;
    break;
  case 'e' :
    os << SEQ_label[SEQL_EQUILIBRIUM_RENEWAL] << endl;
    break;
  }

  // ecriture de la loi inter-evenement

  os << "\n";
  if (file_flag) {
    os << "# ";
  }
  os << SEQ_label[SEQL_INTER_EVENT] << " " << STAT_label[STATL_DISTRIBUTION] << endl;
  inter_event->ascii_print(os);
  inter_event->ascii_parametric_characteristic_print(os , false , file_flag);
  if (file_flag) {
    os << "# ";
  }
  os << STAT_label[STATL_VARIATION_COEFF] << ": "
     << sqrt(inter_event->variance) / inter_event->mean << endl;

  if ((timev) && (timev->within) && (timev->backward) && (timev->forward)) {
    if (timev->inter_event) {
      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << SEQ_label[SEQL_INTER_EVENT] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";
      timev->inter_event->ascii_characteristic_print(os , false , file_flag);
      if (file_flag) {
        os << "# ";
      }
      os << STAT_label[STATL_VARIATION_COEFF] << ": "
         << sqrt(timev->inter_event->variance) / timev->inter_event->mean << endl;

      if (exhaustive) {
        os << "\n";
        if (file_flag) {
          os << "# ";
        }
        os << "   | " << SEQ_label[SEQL_INTER_EVENT] << " " << STAT_label[STATL_DISTRIBUTION]
           << " | " << STAT_label[STATL_CUMULATIVE] << " " << SEQ_label[SEQL_INTER_EVENT] << " "
           << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " " << STAT_label[STATL_FUNCTION]
           << " | " << STAT_label[STATL_CUMULATIVE] << " " << SEQ_label[SEQL_INTER_EVENT] << " "
           << STAT_label[STATL_DISTRIBUTION] << " " << STAT_label[STATL_FUNCTION] << endl;

        inter_event->Distribution::ascii_print(os , file_flag , true , false , timev->inter_event);
      }
    }

    if (exhaustive) {
      pdist = new const Distribution*[1];
      scale = new double[1];

      pdist[0] = inter_event;

      // ecriture de la loi empirique des intervalles de temps a l'interieur de la periode d'observation

      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << SEQ_label[SEQL_OBSERVATION_INTER_EVENT] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";
      timev->within->ascii_characteristic_print(os , false , file_flag);

      if (timev->within->nb_element > 0) {
        os << "\n";
        if (file_flag) {
          os << "# ";
        }
        os << "   | " << SEQ_label[SEQL_OBSERVATION_INTER_EVENT] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
           << " | " << SEQ_label[SEQL_INTER_EVENT] << " " << STAT_label[STATL_DISTRIBUTION]
           << " | " << STAT_label[STATL_CUMULATIVE] << " " << SEQ_label[SEQL_OBSERVATION_INTER_EVENT] << " "
           << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " " << STAT_label[STATL_FUNCTION]
           << " | " << STAT_label[STATL_CUMULATIVE] << " " << SEQ_label[SEQL_INTER_EVENT] << " "
           << STAT_label[STATL_DISTRIBUTION] << " " << STAT_label[STATL_FUNCTION] << endl;

        inter_event->Distribution::ascii_print(os , file_flag , true , false , timev->within);
      }

      // ecriture de la loi biaisee par la longueur

      if (timev->length_bias) {
        os << "\n";
        if (file_flag) {
          os << "# ";
        }
        os << SEQ_label[SEQL_LENGTH_BIASED] << " " << STAT_label[STATL_DISTRIBUTION] << endl;
        length_bias->ascii_characteristic_print(os , false , file_flag);

        os << "\n";
        if (file_flag) {
          os << "# ";
        }
        os << SEQ_label[SEQL_LENGTH_BIASED] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";
        timev->length_bias->ascii_characteristic_print(os , false , file_flag);

        scale[0] = timev->length_bias->nb_element;

        os << "\n";
        if (file_flag) {
          os << "# ";
        }
        os << "   | " << SEQ_label[SEQL_LENGTH_BIASED] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
           << " | " << SEQ_label[SEQL_LENGTH_BIASED] << " " << STAT_label[STATL_DISTRIBUTION]
           << " | " << SEQ_label[SEQL_INTER_EVENT] << " " << STAT_label[STATL_DISTRIBUTION]
           << " | " << STAT_label[STATL_CUMULATIVE] << " " << SEQ_label[SEQL_LENGTH_BIASED] << " "
           << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " " << STAT_label[STATL_FUNCTION]
           << " | " << STAT_label[STATL_CUMULATIVE] << " " << SEQ_label[SEQL_LENGTH_BIASED] << " "
           << STAT_label[STATL_DISTRIBUTION] << " " << STAT_label[STATL_FUNCTION] << endl;

        length_bias->Distribution::ascii_print(os , 1 , pdist , scale , file_flag , true ,
                                               timev->length_bias , true);
      }

      // ecriture de la loi des intervalles de temps apres le dernier evenement

      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << SEQ_label[SEQL_BACKWARD] << " " << SEQ_label[SEQL_RECURRENCE_TIME] << " "
         << STAT_label[STATL_DISTRIBUTION] << endl;
      backward->ascii_characteristic_print(os , false , file_flag);

      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << SEQ_label[SEQL_BACKWARD] << " " << SEQ_label[SEQL_RECURRENCE_TIME] << " "
         << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";
      timev->backward->ascii_characteristic_print(os , false , file_flag);

      scale[0] = timev->backward->nb_element;

      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << "   | " << SEQ_label[SEQL_BACKWARD] << " " << SEQ_label[SEQL_RECURRENCE_TIME]
         << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " | " << SEQ_label[SEQL_BACKWARD]
         << " " << SEQ_label[SEQL_RECURRENCE_TIME] << " " << STAT_label[STATL_DISTRIBUTION]
         << " | " << SEQ_label[SEQL_INTER_EVENT] << " " << STAT_label[STATL_DISTRIBUTION]
         << " | " << STAT_label[STATL_CUMULATIVE] << " " << SEQ_label[SEQL_BACKWARD]
         << " " << SEQ_label[SEQL_RECURRENCE_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
         << " " << STAT_label[STATL_FUNCTION] << " | " << STAT_label[STATL_CUMULATIVE]
         << " " << SEQ_label[SEQL_BACKWARD] << " " << SEQ_label[SEQL_RECURRENCE_TIME]
         << " " << STAT_label[STATL_DISTRIBUTION] << " " << STAT_label[STATL_FUNCTION] << endl;

      backward->Distribution::ascii_print(os , 1 , pdist , scale , file_flag , true ,
                                          timev->backward , true);

      // ecriture de la loi des intervalles de temps residuel

      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << SEQ_label[SEQL_FORWARD] << " " << SEQ_label[SEQL_RECURRENCE_TIME] << " "
         << STAT_label[STATL_DISTRIBUTION] << endl;
      forward->ascii_characteristic_print(os , false , file_flag);

      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << SEQ_label[SEQL_FORWARD] << " " << SEQ_label[SEQL_RECURRENCE_TIME] << " "
         << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";
      timev->forward->ascii_characteristic_print(os , false , file_flag);

      scale[0] = timev->forward->nb_element;

      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << "   | " << SEQ_label[SEQL_FORWARD] << " " << SEQ_label[SEQL_RECURRENCE_TIME]
         << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " | " << SEQ_label[SEQL_FORWARD]
         << " " << SEQ_label[SEQL_RECURRENCE_TIME]<< " " << STAT_label[STATL_DISTRIBUTION]
         << " | " << SEQ_label[SEQL_INTER_EVENT] << " " << STAT_label[STATL_DISTRIBUTION]
         << " | " << STAT_label[STATL_CUMULATIVE] << " " << SEQ_label[SEQL_FORWARD]
         << " " << SEQ_label[SEQL_RECURRENCE_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
         << " " << STAT_label[STATL_FUNCTION] << " | " << STAT_label[STATL_CUMULATIVE]
         << " " << SEQ_label[SEQL_FORWARD] << " " << SEQ_label[SEQL_RECURRENCE_TIME]
         << " " << STAT_label[STATL_DISTRIBUTION] << " " << STAT_label[STATL_FUNCTION] << endl;

      forward->Distribution::ascii_print(os , 1 , pdist , scale , file_flag , true ,
                                         timev->forward , true);

      delete [] pdist;
      delete [] scale;
    }
  }

  if ((exhaustive) && ((!timev) || (!(timev->length_bias)))) {

    // ecriture de la loi inter-evenement, de la loi biaisee par la longueur,
    // de la loi des intervalles de temps apres le dernier evenement et
    // de la loi des intervalles de temps residuel

    os << "\n";
    if (file_flag) {
      os << "# ";
    }
    os << SEQ_label[SEQL_LENGTH_BIASED] << " " << STAT_label[STATL_DISTRIBUTION] << endl;
    length_bias->ascii_characteristic_print(os , false , file_flag);

    os << "\n";
    if (file_flag) {
      os << "# ";
    }
    os << SEQ_label[SEQL_BACKWARD] << " " << SEQ_label[SEQL_RECURRENCE_TIME] << " "
       << STAT_label[STATL_DISTRIBUTION] << endl;
    backward->ascii_characteristic_print(os , false , file_flag);

    os << "\n";
    if (file_flag) {
      os << "# ";
    }
    os << SEQ_label[SEQL_FORWARD] << " " << SEQ_label[SEQL_RECURRENCE_TIME] << " "
       << STAT_label[STATL_DISTRIBUTION] << endl;
    forward->ascii_characteristic_print(os , false , file_flag);

    pdist = new const Distribution*[3];
    scale = new double[3];

    pdist[0] = length_bias;
    scale[0] = 1.;
    pdist[1] = backward;
    scale[1] = 1.;
    pdist[2] = forward;
    scale[2] = 1.;

    os << "\n";
    if (file_flag) {
      os << "# ";
    }
    os << "   | " << SEQ_label[SEQL_INTER_EVENT] << " " << STAT_label[STATL_DISTRIBUTION]
       << " | " << SEQ_label[SEQL_LENGTH_BIASED] << " " << STAT_label[STATL_DISTRIBUTION]
       << " | " << SEQ_label[SEQL_BACKWARD] << " " << SEQ_label[SEQL_RECURRENCE_TIME]
       << " " << STAT_label[STATL_DISTRIBUTION] << " | " << SEQ_label[SEQL_FORWARD]
       << " " << SEQ_label[SEQL_RECURRENCE_TIME] << " " << STAT_label[STATL_DISTRIBUTION]
       << " | " << STAT_label[STATL_CUMULATIVE] << " " << SEQ_label[SEQL_INTER_EVENT]
       << " " << STAT_label[STATL_DISTRIBUTION] << " " << STAT_label[STATL_FUNCTION] << endl;

    inter_event->Distribution::ascii_print(os , 3 , pdist , scale , file_flag , true ,
                                           NULL , true);

    delete [] pdist;
    delete [] scale;
  }

  if (exhaustive) {

    // ecriture des lois du temps avant le neme evenement

    inf = (timev ? timev->mixture->offset : mixture->offset);
    if (inf < 1) {
      inf = 1;
    }

    sup = (timev ? timev->mixture->nb_value : mixture->nb_value);

    if (inf < sup) {
      pdist = new const Distribution*[sup];
      scale = new double[sup];
      nb_dist = 0;

      for (i = inf + 1;i < sup;i++) {
        pdist[nb_dist] = nevent_time[i];
        scale[nb_dist++] = 1.;
      }

      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << "  ";
      for (i = inf;i < sup;i++) {
        os << " | " << SEQ_label[SEQL_TIME_UP] << " " << i << " " << STAT_label[STATL_DISTRIBUTION];
      }
      os << endl;

      nevent_time[inf]->Distribution::ascii_print(os , nb_dist , pdist , scale , file_flag ,
                                                  false , NULL , true);

      delete [] pdist;
      delete [] scale;
    }
  }

  // ecriture des lois du nombre d'evenements

  for (i = time->offset;i < time->nb_value;i++) {
    if (time->mass[i] > 0.) {
      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << SEQ_label[SEQL_NB_EVENT] << " " << SEQ_label[SEQL_DURING] << " " << i << " "
         << SEQ_label[SEQL_TIME_UNIT] << " " << STAT_label[STATL_DISTRIBUTION] << endl;
      nb_event[i]->ascii_characteristic_print(os , false , file_flag);
      if (file_flag) {
        os << "# ";
      }
      os << STAT_label[STATL_VARIANCE_MEAN_RATIO] << ": " << nb_event[i]->variance / nb_event[i]->mean << endl;
      if (nb_event[i]->variance > 0.) {
        if (file_flag) {
          os << "# ";
        }
        os << STAT_label[STATL_SKEWNESS_COEFF] << ": " << nb_event[i]->skewness_computation() << "   "
           << STAT_label[STATL_KURTOSIS_COEFF] << ": " << nb_event[i]->kurtosis_computation() << endl;
      }

      switch (type) {

      case 'o' : {
        os << "\n";
        if (file_flag) {
          os << "# ";
        }
        os << SEQ_label[SEQL_1_CENSORED_INTER_EVENT] << ": " << 1. / (nb_event[i]->mean + 1.) << endl;

        if (file_flag) {
          os << "# ";
        }
        os << SEQ_label[SEQL_COMPLETE_INTER_EVENT] << ": "
           << nb_event[i]->mean / (nb_event[i]->mean + 1.) << endl;
        break;
      }

      case 'e' : {
        os << "\n";
        if (file_flag) {
          os << "# ";
        }
        os << SEQ_label[SEQL_2_CENSORED_INTER_EVENT] << ": "
           << nb_event[i]->mass[0] / (nb_event[i]->mean + 1.) << endl;

        if (file_flag) {
          os << "# ";
        }
        os << SEQ_label[SEQL_1_CENSORED_INTER_EVENT] << ": "
           << 2. * (1. - nb_event[i]->mass[0]) / (nb_event[i]->mean + 1.) << endl;

        if (file_flag) {
          os << "# ";
        }
        os << SEQ_label[SEQL_COMPLETE_INTER_EVENT] << ": "
           << (nb_event[i]->mean - 1. + nb_event[i]->mass[0]) / (nb_event[i]->mean + 1.) << endl;
        break;
      }
      }

      if (timev) {
        os << "\n";
        if (file_flag) {
          os << "# ";
        }
        os << SEQ_label[SEQL_NB_EVENT] << " " << SEQ_label[SEQL_DURING] << " " << i << " "
           << SEQ_label[SEQL_TIME_UNIT] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";
        timev->hnb_event[i]->ascii_characteristic_print(os , false , file_flag);
        if (timev->hnb_event[i]->mean > 0.) {
          if (file_flag) {
            os << "# ";
          }
          os << STAT_label[STATL_VARIANCE_MEAN_RATIO] << ": "
             << timev->hnb_event[i]->variance / timev->hnb_event[i]->mean << endl;
        }
        if (timev->hnb_event[i]->variance > 0.) {
          if (file_flag) {
            os << "# ";
          }
          os << STAT_label[STATL_SKEWNESS_COEFF] << ": " << timev->hnb_event[i]->skewness_computation() << "   "
             << STAT_label[STATL_KURTOSIS_COEFF] << ": " << timev->hnb_event[i]->kurtosis_computation() << endl;
        }

        switch (type) {

        case 'o' : {
          os << "\n";
          if (file_flag) {
            os << "# ";
          }
          os << SEQ_label[SEQL_1_CENSORED_INTER_EVENT] << ": " << timev->hnb_event[i]->nb_element << " ("
             << 1. / (timev->hnb_event[i]->mean + 1.) << ")" << endl;

          if (file_flag) {
            os << "# ";
          }
          os << SEQ_label[SEQL_COMPLETE_INTER_EVENT] << ": " << timev->hnb_event[i]->mean *
                timev->hnb_event[i]->nb_element << " ("
             << timev->hnb_event[i]->mean / (timev->hnb_event[i]->mean + 1.) << ")" << endl;
          break;
        }

        case 'e' : {
          os << "\n";
          if (file_flag) {
            os << "# ";
          }
          os << SEQ_label[SEQL_2_CENSORED_INTER_EVENT] << ": " << timev->hnb_event[i]->frequency[0] << " ("
             << timev->hnb_event[i]->frequency[0] / (timev->hnb_event[i]->nb_element * (timev->hnb_event[i]->mean + 1.))
             << ")" << endl;

          if (file_flag) {
            os << "# ";
          }
          os << SEQ_label[SEQL_1_CENSORED_INTER_EVENT] << ": " << 2 * (timev->hnb_event[i]->nb_element -
                 timev->hnb_event[i]->frequency[0]) << " ("
             << 2. * (timev->hnb_event[i]->nb_element - timev->hnb_event[i]->frequency[0]) /
                (timev->hnb_event[i]->nb_element * (timev->hnb_event[i]->mean + 1.)) << ")" << endl;

          if (file_flag) {
            os << "# ";
          }
          os << SEQ_label[SEQL_COMPLETE_INTER_EVENT] << ": " << (timev->hnb_event[i]->mean - 1.) *
                timev->hnb_event[i]->nb_element + timev->hnb_event[i]->frequency[0] << " ("
             << (timev->hnb_event[i]->mean - 1. + (double)timev->hnb_event[i]->frequency[0] /
                 (double)timev->hnb_event[i]->nb_element) / (timev->hnb_event[i]->mean + 1.) << ")" << endl;
          break;
        }
        }

        likelihood = nb_event[i]->likelihood_computation(*timev->hnb_event[i]);
        information = timev->hnb_event[i]->information_computation();

        os << "\n";
        if (file_flag) {
          os << "# ";
        }
        os << STAT_label[STATL_LIKELIHOOD] << ": "<< likelihood << "   ("
           << STAT_label[STATL_NORMALIZED] << ": " << likelihood / timev->hnb_event[i]->nb_element << ")" << endl;

        if (file_flag) {
          os << "# ";
        }
        os << STAT_label[STATL_MAX_LIKELIHOOD] << ": " << information << "   ("
           << STAT_label[STATL_INFORMATION] << ": " << information / timev->hnb_event[i]->nb_element << ")" << endl;

        if (file_flag) {
          os << "# ";
        }
        os << STAT_label[STATL_DEVIANCE] << ": " << 2 * (information - likelihood) << endl;

        nb_event[i]->chi2_fit(*(timev->hnb_event[i]) , test);
        os << "\n";
        test.ascii_print(os , file_flag);
      }

      if ((exhaustive) && ((time->variance == 0.) || (timev))) {
        os << "\n";
        if (file_flag) {
          os << "# ";
        }
        os << "  ";
        if (timev) {
          os << " | " << SEQ_label[SEQL_NB_EVENT] << " " << SEQ_label[SEQL_DURING] << " " << i
             << " " << SEQ_label[SEQL_TIME_UNIT] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
        }
        os << " | " << SEQ_label[SEQL_NB_EVENT] << " " << SEQ_label[SEQL_DURING] << " " << i
           << " " << SEQ_label[SEQL_TIME_UNIT] << " " << STAT_label[STATL_DISTRIBUTION];
        if (timev) {
          os << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " "
             << STAT_label[STATL_FUNCTION];
        }
        os << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_DISTRIBUTION] << " "
           << STAT_label[STATL_FUNCTION] << endl;

        nb_event[i]->Distribution::ascii_print(os , file_flag , true , false ,
                                               (timev ? timev->hnb_event[i] : NULL));
      }
    }
  }

  if (time->variance > 0.) {

    // ecriture du melange de lois du nombre d'evenements

    if (exhaustive) {
      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << SEQ_label[SEQL_NB_EVENT_MIXTURE] << endl;
      mixture->ascii_characteristic_print(os , false , file_flag);

      switch (type) {

      case 'o' : {
        os << "\n";
        if (file_flag) {
          os << "# ";
        }
        os << SEQ_label[SEQL_1_CENSORED_INTER_EVENT] << ": " << 1. / (mixture->mean + 1.) << endl;

        if (file_flag) {
          os << "# ";
        }
        os << SEQ_label[SEQL_COMPLETE_INTER_EVENT] << ": "
           << mixture->mean / (mixture->mean + 1.) << endl;
        break;
      }

      case 'e' : {
        os << "\n";
        if (file_flag) {
          os << "# ";
        }
        os << SEQ_label[SEQL_2_CENSORED_INTER_EVENT] << ": "
           << mixture->mass[0] / (mixture->mean + 1.) << endl;

        if (file_flag) {
          os << "# ";
        }
        os << SEQ_label[SEQL_1_CENSORED_INTER_EVENT] << ": "
           << 2. * (1. - mixture->mass[0]) / (mixture->mean + 1.) << endl;

        if (file_flag) {
          os << "# ";
        }
        os << SEQ_label[SEQL_COMPLETE_INTER_EVENT] << ": "
           << (mixture->mean - 1. + mixture->mass[0]) / (mixture->mean + 1.) << endl;
        break;
      }
      }

      if (timev) {
        os << "\n";
        if (file_flag) {
          os << "# ";
        }
        os << SEQ_label[SEQL_NB_EVENT] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";
        timev->mixture->ascii_characteristic_print(os , false , file_flag);

        switch (type) {

        case 'o' : {
          os << "\n";
          if (file_flag) {
            os << "# ";
          }
          os << SEQ_label[SEQL_1_CENSORED_INTER_EVENT] << ": " << timev->mixture->nb_element << " ("
             << 1. / (timev->mixture->mean + 1.) << ")" << endl;

          if (file_flag) {
            os << "# ";
          }
          os << SEQ_label[SEQL_COMPLETE_INTER_EVENT] << ": " << timev->mixture->mean *
                timev->mixture->nb_element << " ("
             << timev->mixture->mean / (timev->mixture->mean + 1.) << ")" << endl;
          break;
        }

        case 'e' : {
          os << "\n";
          if (file_flag) {
            os << "# ";
          }
          os << SEQ_label[SEQL_2_CENSORED_INTER_EVENT] << ": " << timev->mixture->frequency[0] << " ("
             << timev->mixture->frequency[0] / (timev->mixture->nb_element * (timev->mixture->mean + 1.))
             << ")" << endl;

          if (file_flag) {
            os << "# ";
          }
          os << SEQ_label[SEQL_1_CENSORED_INTER_EVENT] << ": " << 2 * (timev->mixture->nb_element -
                 timev->mixture->frequency[0]) << " ("
             << 2. * (timev->mixture->nb_element - timev->mixture->frequency[0]) /
                (timev->mixture->nb_element * (timev->mixture->mean + 1.)) << ")" << endl;

          if (file_flag) {
            os << "# ";
          }
          os << SEQ_label[SEQL_COMPLETE_INTER_EVENT] << ": " << (timev->mixture->mean - 1.) *
                timev->mixture->nb_element + timev->mixture->frequency[0] << " ("
             << (timev->mixture->mean - 1. + (double)timev->mixture->frequency[0] /
                 (double)timev->mixture->nb_element) / (timev->mixture->mean + 1.) << ")" << endl;
          break;
        }
        }

        likelihood = likelihood_computation(*timev);
        information = timev->information_computation();

        os << "\n";
        if (file_flag) {
          os << "# ";
        }
        os << STAT_label[STATL_LIKELIHOOD] << ": " << likelihood << "   ("
           << STAT_label[STATL_NORMALIZED] << ": " << likelihood / timev->nb_element << ")" << endl;
        if (file_flag) {
          os << "# ";
        }
        os << STAT_label[STATL_MAX_LIKELIHOOD] << ": " << information << "   ("
           << STAT_label[STATL_INFORMATION] << ": " << information / timev->nb_element << ")" << endl;

        if (file_flag) {
          os << "# ";
        }
        os << STAT_label[STATL_DEVIANCE] << ": " << 2 * (information - likelihood) << endl;

        mixture->chi2_fit(*(timev->mixture) , test);
        os << "\n";
        test.ascii_print(os , file_flag);
      }

      pdist = new const Distribution*[time->nb_value];
      scale = new double[time->nb_value];
      nb_dist = 0;

      for (i = time->offset;i < time->nb_value;i++) {
        if (time->mass[i] > 0.) {
          pdist[nb_dist] = nb_event[i];

          if (timev) {
            scale[nb_dist++] = timev->nb_element * time->mass[i];
          }
          else {
            scale[nb_dist++] = time->mass[i];
          }
        }
      }

      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << "  ";
      if (timev) {
        os << " | " << SEQ_label[SEQL_NB_EVENT] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
      }
      os << " | " << SEQ_label[SEQL_NB_EVENT_MIXTURE];
      for (i = time->offset;i < time->nb_value;i++) {
        if (time->mass[i] > 0.) {
          os << " | " << SEQ_label[SEQL_NB_EVENT] << " " << SEQ_label[SEQL_DURING] << " " << i
             << " " << SEQ_label[SEQL_TIME_UNIT] << " " << STAT_label[STATL_DISTRIBUTION];
        }
      }
      if (timev) {
        os << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " "
           << STAT_label[STATL_FUNCTION];
      }
      os << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_MIXTURE] << " "
         << STAT_label[STATL_FUNCTION] << endl;

      mixture->ascii_print(os , nb_dist , pdist , scale , file_flag , true ,
                           (timev ? timev->mixture : NULL));

      delete [] pdist;
      delete [] scale;
    }

    // ecriture temps d'observation

    if (timev) {
      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << SEQ_label[SEQL_OBSERVATION_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";
      timev->htime->ascii_characteristic_print(os , false , file_flag);

      if (exhaustive) {
        os << "\n";
        if (file_flag) {
          os << "# ";
        }
        os << "   | " << SEQ_label[SEQL_OBSERVATION_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << endl;
        timev->htime->ascii_print(os , file_flag);
      }
    }
  }

  if (exhaustive) {

    // ecriture des probabilites de non-evenement/evenement fonction du temps

    os << "\n";
    if (file_flag) {
      os << "# ";
    }
    os << "  ";

    if ((timev) && (timev->index_event)) {
      os << " | " << SEQ_label[SEQL_OBSERVED] << " " << SEQ_label[SEQL_NO_EVENT_PROBABILITY];
    }
    os << " | " << SEQ_label[SEQL_THEORETICAL] << " " << SEQ_label[SEQL_NO_EVENT_PROBABILITY];
    if ((timev) && (timev->index_event)) {
      os << " | " << SEQ_label[SEQL_OBSERVED] << " " << SEQ_label[SEQL_EVENT_PROBABILITY];
    }
    os << " | " << SEQ_label[SEQL_THEORETICAL] << " " << SEQ_label[SEQL_EVENT_PROBABILITY];
    if ((timev) && (timev->index_event)) {
      os << " | " << STAT_label[STATL_FREQUENCY];
    }
    os << endl;

    index_event->ascii_print(os , file_flag , (((timev) && (timev->index_event)) ? timev->index_event : NULL));

    // ecriture des sequences d'evenements

    if ((timev) && (timev->sequence)) {
      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << SEQ_label[SEQL_SEQUENCES] << endl;

      for (i = 0;i < timev->nb_element;i++) {
        os << "\n";
        if (file_flag) {
          os << "# ";
        }

        for (j = 0;j < timev->length[i];j++) {
          if ((j > 0) && ((2 * j) % LINE_NB_CHARACTER == 0)) {
            os << "\\" << endl;
            if (file_flag) {
              os << "# ";
            }
          }

          os << timev->sequence[i][j] << " ";
        }

        os << endl;
      }
    }
  }
  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Renewal.
 *
 *  arguments : stream, flag niveau de detail.
 *
 *--------------------------------------------------------------*/

ostream& Renewal::ascii_write(ostream &os , bool exhaustive) const

{
  return ascii_write(os , renewal_data , exhaustive , false);
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Renewal dans un fichier.
 *
 *  arguments : reference sur un objet StatError, path,
 *              flag niveau de detail.
 *
 *--------------------------------------------------------------*/

bool Renewal::ascii_write(StatError &error , const char *path ,
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
    ascii_write(out_file , renewal_data , exhaustive , true);
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un processus de renouvellement et de la structure
 *  de donnees associee au format tableur.
 *
 *  arguments : stream, pointeur sur un objet RenewalData.
 *
 *--------------------------------------------------------------*/

ostream& Renewal::spreadsheet_write(ostream &os , const RenewalData *timev) const

/* {
  os << inter_event->cumul[1] << "\t" << inter_event->cumul[4] << "\t" << inter_event->cumul[12] << "\t"
  os << inter_event->cumul[6] << "\t" << inter_event->cumul[9] << "\t" << inter_event->cumul[12] << "\t"
     << inter_event->mean << "\t" << sqrt(inter_event->variance) << "\t"
     << inter_event->second_difference_norm_computation() << "\t"
     << timev->hnb_event[time->offset]->mean << "\t" << nb_event[time->offset]->mean << "\t"
     << timev->hnb_event[time->offset]->variance << "\t" << nb_event[time->offset]->variance << "\t"
     << 2 * (timev->information_computation()- likelihood_computation(*timev)) << endl;

  return os;
} */

{
  register int i;
  int nb_dist , inf , sup;
  double likelihood , information , *scale;
  const Distribution **pdist;
  Test test(CHI2);


  switch (type) {
  case 'o' :
    os << SEQ_label[SEQL_ORDINARY_RENEWAL] << endl;
    break;
  case 'e' :
    os << SEQ_label[SEQL_EQUILIBRIUM_RENEWAL] << endl;
    break;
  }

  // ecriture de la loi inter-evenement

  os << "\n" << SEQ_label[SEQL_INTER_EVENT] << " " << STAT_label[STATL_DISTRIBUTION] << endl;
  inter_event->spreadsheet_print(os);
  inter_event->spreadsheet_parametric_characteristic_print(os);
  os << STAT_label[STATL_VARIATION_COEFF] << "\t"
     << sqrt(inter_event->variance) / inter_event->mean << endl;

  if ((timev) && (timev->within) && (timev->backward) && (timev->forward)) {
    if (timev->inter_event) {
      os << "\n" << SEQ_label[SEQL_INTER_EVENT] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\t";
      timev->inter_event->spreadsheet_characteristic_print(os);
      os << STAT_label[STATL_VARIATION_COEFF] << "\t"
         << sqrt(timev->inter_event->variance) / timev->inter_event->mean << endl;

      os << "\n\t" << SEQ_label[SEQL_INTER_EVENT] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
         << "\t" << SEQ_label[SEQL_INTER_EVENT] << " " << STAT_label[STATL_DISTRIBUTION]
         << "\t" << STAT_label[STATL_CUMULATIVE] << " " << SEQ_label[SEQL_INTER_EVENT] << " "
         << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " " << STAT_label[STATL_FUNCTION]
         << "\t" << STAT_label[STATL_CUMULATIVE] << " " << SEQ_label[SEQL_INTER_EVENT] << " "
         << STAT_label[STATL_DISTRIBUTION] << " " << STAT_label[STATL_FUNCTION] << endl;

      inter_event->Distribution::spreadsheet_print(os , true , false , false , timev->inter_event);
    }

    pdist = new const Distribution*[1];
    scale = new double[1];

    pdist[0] = inter_event;

    // ecriture de la loi empirique des intervalles de temps a l'interieur de la periode d'observation

    os << "\n" << SEQ_label[SEQL_OBSERVATION_INTER_EVENT] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\t";
    timev->within->spreadsheet_characteristic_print(os);

    if (timev->within->nb_element > 0) {
      os << "\n\t" << SEQ_label[SEQL_OBSERVATION_INTER_EVENT] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
         << "\t" << SEQ_label[SEQL_INTER_EVENT] << " " << STAT_label[STATL_DISTRIBUTION]
         << "\t" << STAT_label[STATL_CUMULATIVE] << " " << SEQ_label[SEQL_OBSERVATION_INTER_EVENT] << " "
         << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " " << STAT_label[STATL_FUNCTION]
         << "\t" << STAT_label[STATL_CUMULATIVE] << " " << SEQ_label[SEQL_INTER_EVENT] << " "
         << STAT_label[STATL_DISTRIBUTION] << " " << STAT_label[STATL_FUNCTION] << endl;

      inter_event->Distribution::spreadsheet_print(os , true , false , false , timev->within);
    }

    // ecriture de la loi biaisee par la longueur

    if (timev->length_bias) {
      os << "\n" << SEQ_label[SEQL_LENGTH_BIASED] << " " << STAT_label[STATL_DISTRIBUTION] << endl;
      length_bias->spreadsheet_characteristic_print(os);
      os << "\n" << SEQ_label[SEQL_LENGTH_BIASED] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\t";
      timev->length_bias->spreadsheet_characteristic_print(os);

      scale[0] = timev->length_bias->nb_element;

      os << "\n\t" << SEQ_label[SEQL_LENGTH_BIASED] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
         << "\t" << SEQ_label[SEQL_LENGTH_BIASED] << " " << STAT_label[STATL_DISTRIBUTION]
         << "\t" << SEQ_label[SEQL_INTER_EVENT] << " " << STAT_label[STATL_DISTRIBUTION]
         << "\t" << STAT_label[STATL_CUMULATIVE] << " " << SEQ_label[SEQL_LENGTH_BIASED] << " "
         << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " " << STAT_label[STATL_FUNCTION]
         << "\t" << STAT_label[STATL_CUMULATIVE] << " " << SEQ_label[SEQL_LENGTH_BIASED] << " "
         << STAT_label[STATL_DISTRIBUTION] << " " << STAT_label[STATL_FUNCTION] << endl;

      length_bias->Distribution::spreadsheet_print(os , 1 , pdist , scale , true ,
                                                   timev->length_bias , true);
    }

    // ecriture de la loi des intervalles de temps apres le dernier evenement

    os << "\n" << SEQ_label[SEQL_BACKWARD] << " " << SEQ_label[SEQL_RECURRENCE_TIME] << " "
       << STAT_label[STATL_DISTRIBUTION] << endl;
    backward->spreadsheet_characteristic_print(os);
    os << "\n" << SEQ_label[SEQL_BACKWARD] << " " << SEQ_label[SEQL_RECURRENCE_TIME] << " "
       << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\t";
    timev->backward->spreadsheet_characteristic_print(os);

    scale[0] = timev->backward->nb_element;

    os << "\n\t" << SEQ_label[SEQL_BACKWARD] << " " << SEQ_label[SEQL_RECURRENCE_TIME]
       << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\t" << SEQ_label[SEQL_BACKWARD]
       << " " << SEQ_label[SEQL_RECURRENCE_TIME] << " " << STAT_label[STATL_DISTRIBUTION]
       << "\t" << SEQ_label[SEQL_INTER_EVENT] << " " << STAT_label[STATL_DISTRIBUTION]
       << "\t" << STAT_label[STATL_CUMULATIVE] << " " << SEQ_label[SEQL_BACKWARD]
       << " " << SEQ_label[SEQL_RECURRENCE_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
       << " " << STAT_label[STATL_FUNCTION] << "\t" << STAT_label[STATL_CUMULATIVE]
       << " " << SEQ_label[SEQL_BACKWARD] << " " << SEQ_label[SEQL_RECURRENCE_TIME]
       << " " << STAT_label[STATL_DISTRIBUTION] << " " << STAT_label[STATL_FUNCTION] << endl;

    backward->Distribution::spreadsheet_print(os , 1 , pdist , scale , true ,
                                              timev->backward , true);

    // ecriture de la loi des intervalles de temps residuel

    os << "\n" << SEQ_label[SEQL_FORWARD] << " " << SEQ_label[SEQL_RECURRENCE_TIME]
       << " " << STAT_label[STATL_DISTRIBUTION] << endl;
    forward->spreadsheet_characteristic_print(os);
    os << "\n" << SEQ_label[SEQL_FORWARD] << " " << SEQ_label[SEQL_RECURRENCE_TIME]
       << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\t";
    timev->forward->spreadsheet_characteristic_print(os);

    scale[0] = timev->forward->nb_element;

    os << "\n\t" << SEQ_label[SEQL_FORWARD] << " " << SEQ_label[SEQL_RECURRENCE_TIME]
       << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\t" << SEQ_label[SEQL_FORWARD]
       << " " << SEQ_label[SEQL_RECURRENCE_TIME] << " " << STAT_label[STATL_DISTRIBUTION]
       << "\t" << SEQ_label[SEQL_INTER_EVENT] << " " << STAT_label[STATL_DISTRIBUTION]
       << "\t" << STAT_label[STATL_CUMULATIVE] << " " << SEQ_label[SEQL_FORWARD]
       << " " << SEQ_label[SEQL_RECURRENCE_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
       << " " << STAT_label[STATL_FUNCTION] << "\t" << STAT_label[STATL_CUMULATIVE]
       << " " << SEQ_label[SEQL_FORWARD] << " " << SEQ_label[SEQL_RECURRENCE_TIME]
       << " " << STAT_label[STATL_DISTRIBUTION] << " " << STAT_label[STATL_FUNCTION] << endl;

    forward->Distribution::spreadsheet_print(os , 1 , pdist , scale , true ,
                                             timev->forward , true);

    delete [] pdist;
    delete [] scale;
  }

  if ((!timev) || (!(timev->length_bias))) {

    // ecriture de la loi inter-evenement, de la loi biaisee par la longueur
    // de la loi des intervalles de temps apres le dernier evenement et
    // de la loi des intervalles de temps residuel

    os << "\n" << SEQ_label[SEQL_LENGTH_BIASED] << " " << STAT_label[STATL_DISTRIBUTION] << endl;
    length_bias->spreadsheet_characteristic_print(os);
    os << "\n" << SEQ_label[SEQL_BACKWARD] << " " << SEQ_label[SEQL_RECURRENCE_TIME] << " "
       << STAT_label[STATL_DISTRIBUTION] << endl;
    backward->spreadsheet_characteristic_print(os);
    os << "\n" << SEQ_label[SEQL_FORWARD] << " " << SEQ_label[SEQL_RECURRENCE_TIME] << " "
       << STAT_label[STATL_DISTRIBUTION] << endl;
    forward->spreadsheet_characteristic_print(os);

    pdist = new const Distribution*[3];
    scale = new double[3];

    pdist[0] = length_bias;
    scale[0] = 1.;
    pdist[1] = backward;
    scale[1] = 1.;
    pdist[2] = forward;
    scale[2] = 1.;

    os << "\n\t" << SEQ_label[SEQL_INTER_EVENT] << " " << STAT_label[STATL_DISTRIBUTION]
       << "\t" << SEQ_label[SEQL_LENGTH_BIASED] << " " << STAT_label[STATL_DISTRIBUTION]
       << "\t" << SEQ_label[SEQL_BACKWARD] << " " << SEQ_label[SEQL_RECURRENCE_TIME]
       << " " << STAT_label[STATL_DISTRIBUTION] << "\t" << SEQ_label[SEQL_FORWARD]
       << " " << SEQ_label[SEQL_RECURRENCE_TIME] << " " << STAT_label[STATL_DISTRIBUTION]
       << "\t" << STAT_label[STATL_CUMULATIVE] << " " << SEQ_label[SEQL_INTER_EVENT]
       << " " << STAT_label[STATL_DISTRIBUTION] << " " << STAT_label[STATL_FUNCTION] << endl;

    inter_event->Distribution::spreadsheet_print(os , 3 , pdist , scale , true , NULL , true);

    delete [] pdist;
    delete [] scale;
  }

  // ecriture des lois du temps avant le neme evenement

  inf = (timev ? timev->mixture->offset : mixture->offset);
  if (inf < 1) {
    inf = 1;
  }

  sup = (timev ? timev->mixture->nb_value : mixture->nb_value);

  if (inf < sup) {
    pdist = new const Distribution*[sup];
    scale = new double[sup];
    nb_dist = 0;

    for (i = inf + 1;i < sup;i++) {
      pdist[nb_dist] = nevent_time[i];
      scale[nb_dist++] = 1.;
    }

    os << "\n";
    for (i = inf;i < sup;i++) {
      os << "\t" << SEQ_label[SEQL_TIME_UP] << " " << i << " " << STAT_label[STATL_DISTRIBUTION];
    }
    os << endl;

    nevent_time[inf]->Distribution::spreadsheet_print(os , nb_dist , pdist , scale ,
                                                      false , NULL , true);

    delete [] pdist;
    delete [] scale;
  }

  // ecriture des lois du nombre d'evenements

  for (i = time->offset;i < time->nb_value;i++) {
    if (time->mass[i] > 0.) {
      os << "\n" << SEQ_label[SEQL_NB_EVENT] << " " << SEQ_label[SEQL_DURING] << " " << i
         << " " << SEQ_label[SEQL_TIME_UNIT] << " " << STAT_label[STATL_DISTRIBUTION] << endl;
      nb_event[i]->spreadsheet_characteristic_print(os);
      os << STAT_label[STATL_VARIANCE_MEAN_RATIO] << "\t" << nb_event[i]->variance / nb_event[i]->mean << endl;
      if (nb_event[i]->variance > 0.) {
        os << STAT_label[STATL_SKEWNESS_COEFF] << "\t" << nb_event[i]->skewness_computation() << "\t"
           << STAT_label[STATL_KURTOSIS_COEFF] << "\t" << nb_event[i]->kurtosis_computation() << endl;
      }

      switch (type) {

      case 'o' : {
        os << "\n" << SEQ_label[SEQL_1_CENSORED_INTER_EVENT] << "\t" << 1. / (nb_event[i]->mean + 1.) << endl;

        os << SEQ_label[SEQL_COMPLETE_INTER_EVENT] << "\t"
           << nb_event[i]->mean / (nb_event[i]->mean + 1.) << endl;
        break;
      }

      case 'e' : {
        os << "\n" << SEQ_label[SEQL_2_CENSORED_INTER_EVENT] << "\t"
           << nb_event[i]->mass[0] / (nb_event[i]->mean + 1.) << endl;

        os << SEQ_label[SEQL_1_CENSORED_INTER_EVENT] << "\t"
           << 2. * (1. - nb_event[i]->mass[0]) / (nb_event[i]->mean + 1.) << endl;

        os << SEQ_label[SEQL_COMPLETE_INTER_EVENT] << "\t"
           << (nb_event[i]->mean - 1. + nb_event[i]->mass[0]) / (nb_event[i]->mean + 1.) << endl;
        break;
      }
      }

      if (timev) {
        os << "\n" << SEQ_label[SEQL_NB_EVENT] << " " << SEQ_label[SEQL_DURING] << " " << i
           << " " << SEQ_label[SEQL_TIME_UNIT] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\t";
        timev->hnb_event[i]->spreadsheet_characteristic_print(os);
        if (timev->hnb_event[i]->mean > 0.) {
          os << STAT_label[STATL_VARIANCE_MEAN_RATIO] << "\t"
             << timev->hnb_event[i]->variance / timev->hnb_event[i]->mean << endl;
        }
        if (timev->hnb_event[i]->variance > 0.) {
          os << STAT_label[STATL_SKEWNESS_COEFF] << "\t" << timev->hnb_event[i]->skewness_computation() << "\t"
             << STAT_label[STATL_KURTOSIS_COEFF] << "\t" << timev->hnb_event[i]->kurtosis_computation() << endl;
        }

        switch (type) {

        case 'o' : {
          os << "\n" << SEQ_label[SEQL_1_CENSORED_INTER_EVENT] << "\t" << timev->hnb_event[i]->nb_element << "\t"
             << 1. / (timev->hnb_event[i]->mean + 1.) << endl;

          os << SEQ_label[SEQL_COMPLETE_INTER_EVENT] << "\t" << timev->hnb_event[i]->mean *
                timev->hnb_event[i]->nb_element << "\t"
             << timev->hnb_event[i]->mean / (timev->hnb_event[i]->mean + 1.) << endl;
          break;
        }

        case 'e' : {
          os << "\n" << SEQ_label[SEQL_2_CENSORED_INTER_EVENT] << "\t" << timev->hnb_event[i]->frequency[0] << "\t"
             << timev->hnb_event[i]->frequency[0] / (timev->hnb_event[i]->nb_element * (timev->hnb_event[i]->mean + 1.)) << endl;

          os << SEQ_label[SEQL_1_CENSORED_INTER_EVENT] << "\t" << 2 * (timev->hnb_event[i]->nb_element -
                 timev->hnb_event[i]->frequency[0]) << "\t"
             << 2. * (timev->hnb_event[i]->nb_element - timev->hnb_event[i]->frequency[0]) /
                (timev->hnb_event[i]->nb_element * (timev->hnb_event[i]->mean + 1.)) << endl;

          os << SEQ_label[SEQL_COMPLETE_INTER_EVENT] << "\t" << (timev->hnb_event[i]->mean - 1.) *
                timev->hnb_event[i]->nb_element + timev->hnb_event[i]->frequency[0] << "\t"
             << (timev->hnb_event[i]->mean - 1. + (double)timev->hnb_event[i]->frequency[0] /
                 (double)timev->hnb_event[i]->nb_element) / (timev->hnb_event[i]->mean + 1.) << endl;
          break;
        }
        }

        likelihood = nb_event[i]->likelihood_computation(*timev->hnb_event[i]);
        information = timev->hnb_event[i]->information_computation();

        os << "\n" << STAT_label[STATL_LIKELIHOOD] << "\t" << likelihood << "\t"
           << STAT_label[STATL_NORMALIZED] << "\t" << likelihood / timev->hnb_event[i]->nb_element << endl;
        os << STAT_label[STATL_MAX_LIKELIHOOD] << "\t" << information << "\t"
           << STAT_label[STATL_INFORMATION] << "\t" << information / timev->hnb_event[i]->nb_element << endl;
        os << STAT_label[STATL_DEVIANCE] << "\t" << 2 * (information - likelihood) << endl;

        nb_event[i]->chi2_fit(*(timev->hnb_event[i]) , test);
        os << "\n";
        test.spreadsheet_print(os);
      }

      if ((time->variance == 0.) || (timev)) {
        os << "\n";
        if (timev) {
          os << "\t" << SEQ_label[SEQL_NB_EVENT] << " " << SEQ_label[SEQL_DURING] << " " << i
             << " " << SEQ_label[SEQL_TIME_UNIT] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
        }
        os << "\t" << SEQ_label[SEQL_NB_EVENT] << " " << SEQ_label[SEQL_DURING] << " " << i
           << " " << SEQ_label[SEQL_TIME_UNIT] << " " << STAT_label[STATL_DISTRIBUTION];
        if (timev) {
          os << "\t" << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " "
             << STAT_label[STATL_FUNCTION];
        }
        os << "\t" << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_DISTRIBUTION] << " "
           << STAT_label[STATL_FUNCTION] << endl;

        nb_event[i]->Distribution::spreadsheet_print(os , true , false , false ,
                                                     (timev ? timev->hnb_event[i] : NULL));
      }
    }
  }

  if (time->variance > 0.) {

    // ecriture du melange de lois du nombre d'evenements

    os << "\n" << SEQ_label[SEQL_NB_EVENT_MIXTURE] << endl;
    mixture->spreadsheet_characteristic_print(os);

    switch (type) {

    case 'o' : {
      os << "\n" << SEQ_label[SEQL_1_CENSORED_INTER_EVENT] << "\t" << 1. / (mixture->mean + 1.) << endl;

      os << SEQ_label[SEQL_COMPLETE_INTER_EVENT] << "\t"
         << mixture->mean / (mixture->mean + 1.) << endl;
      break;
    }

    case 'e' : {
      os << "\n" << SEQ_label[SEQL_2_CENSORED_INTER_EVENT] << "\t"
         << mixture->mass[0] / (mixture->mean + 1.) << endl;

      os << SEQ_label[SEQL_1_CENSORED_INTER_EVENT] << "\t"
         << 2. * (1. - mixture->mass[0]) / (mixture->mean + 1.) << endl;

      os << SEQ_label[SEQL_COMPLETE_INTER_EVENT] << "\t"
         << (mixture->mean - 1. + mixture->mass[0]) / (mixture->mean + 1.) << endl;
      break;
    }
    }

    if (timev) {
      os << "\n" << SEQ_label[SEQL_NB_EVENT] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\t";
      timev->mixture->spreadsheet_characteristic_print(os);

      switch (type) {

      case 'o' : {
        os << "\n" << SEQ_label[SEQL_1_CENSORED_INTER_EVENT] << "\t" << timev->mixture->nb_element << "\t"
           << 1. / (timev->mixture->mean + 1.) << endl;

        os << SEQ_label[SEQL_COMPLETE_INTER_EVENT] << "\t" << timev->mixture->mean *
              timev->mixture->nb_element << "\t"
           << timev->mixture->mean / (timev->mixture->mean + 1.) << endl;
        break;
      }

      case 'e' : {
        os << "\n" << SEQ_label[SEQL_2_CENSORED_INTER_EVENT] << "\t" << timev->mixture->frequency[0] << "\t"
           << timev->mixture->frequency[0] / (timev->mixture->nb_element * (timev->mixture->mean + 1.)) << endl;

        os << SEQ_label[SEQL_1_CENSORED_INTER_EVENT] << "\t" << 2 * (timev->mixture->nb_element -
               timev->mixture->frequency[0]) << "\t"
           << 2. * (timev->mixture->nb_element - timev->mixture->frequency[0]) /
              (timev->mixture->nb_element * (timev->mixture->mean + 1.)) << endl;

        os << SEQ_label[SEQL_COMPLETE_INTER_EVENT] << "\t" << (timev->mixture->mean - 1.) *
              timev->mixture->nb_element + timev->mixture->frequency[0] << "\t"
           << (timev->mixture->mean - 1. + (double)timev->mixture->frequency[0] /
               (double)timev->mixture->nb_element) / (timev->mixture->mean + 1.) << endl;
        break;
      }
      }

      likelihood = likelihood_computation(*timev);
      information = timev->information_computation();

      os << "\n" << STAT_label[STATL_LIKELIHOOD] << "\t" << likelihood << "\t"
         << STAT_label[STATL_NORMALIZED] << "\t" << likelihood / timev->nb_element << endl;
      os << STAT_label[STATL_MAX_LIKELIHOOD] << "\t" << information << "\t"
         << STAT_label[STATL_INFORMATION] << "\t" << information / timev->nb_element << endl;
      os << STAT_label[STATL_DEVIANCE] << "\t" << 2 * (information - likelihood) << endl;

      mixture->chi2_fit(*(timev->mixture) , test);
      os << "\n";
      test.spreadsheet_print(os);
    }

    pdist = new const Distribution*[time->nb_value];
    scale = new double[time->nb_value];
    nb_dist = 0;

    for (i = time->offset;i < time->nb_value;i++) {
      if (time->mass[i] > 0.) {
        pdist[nb_dist] = nb_event[i];

        if (timev) {
          scale[nb_dist++] = timev->nb_element * time->mass[i];
        }
        else {
          scale[nb_dist++] = time->mass[i];
        }
      }
    }

    os << "\n";
    if (timev) {
      os << "\t" << SEQ_label[SEQL_NB_EVENT] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
    }
    os << "\t" << SEQ_label[SEQL_NB_EVENT_MIXTURE];
    for (i = time->offset;i < time->nb_value;i++) {
      if (time->mass[i] > 0.) {
        os << "\t" << SEQ_label[SEQL_NB_EVENT] << " " << SEQ_label[SEQL_DURING] << " " << i
           << " " << SEQ_label[SEQL_TIME_UNIT] << " " << STAT_label[STATL_DISTRIBUTION];
      }
    }
    if (timev) {
      os << "\t" << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " "
         << STAT_label[STATL_FUNCTION];
    }
    os << "\t" << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_MIXTURE] << " "
       << STAT_label[STATL_FUNCTION] << endl;

    mixture->spreadsheet_print(os , nb_dist , pdist , scale , true ,
                               (timev ? timev->mixture : NULL));

    delete [] pdist;
    delete [] scale;

    // ecriture temps d'observation

    if (timev) {
      os << "\n" << SEQ_label[SEQL_OBSERVATION_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\t";
      timev->htime->spreadsheet_characteristic_print(os);

      os << "\n\t" << SEQ_label[SEQL_OBSERVATION_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << endl;
      timev->htime->spreadsheet_print(os);
    }
  }

  // ecriture des probabilites de non-evenement/evenement fonction du temps

  os << "\n";
  if ((timev) && (timev->index_event)) {
    os << "\t" << SEQ_label[SEQL_OBSERVED] << " " << SEQ_label[SEQL_NO_EVENT_PROBABILITY];
  }
  os << "\t" << SEQ_label[SEQL_THEORETICAL] << " " << SEQ_label[SEQL_NO_EVENT_PROBABILITY];
  if ((timev) && (timev->index_event)) {
    os << "\t" << SEQ_label[SEQL_OBSERVED] << " " << SEQ_label[SEQL_EVENT_PROBABILITY];
  }
  os << "\t" << SEQ_label[SEQL_THEORETICAL] << " " << SEQ_label[SEQL_EVENT_PROBABILITY];
  if ((timev) && (timev->index_event)) {
    os << "\t" << STAT_label[STATL_FREQUENCY];
  }
  os << endl;

  index_event->spreadsheet_print(os , (((timev) && (timev->index_event)) ? timev->index_event : NULL));

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Renewal dans un fichier au format tableur.
 *
 *  arguments : reference sur un objet StatError, path.
 *
 *--------------------------------------------------------------*/

bool Renewal::spreadsheet_write(StatError &error , const char *path) const

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
    spreadsheet_write(out_file , renewal_data);
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Sortie Gnuplot d'un processus de renouvellement.
 *
 *  arguments : prefixe des fichiers, titre des figures,
 *              pointeur sur un objet RenewalData.
 *
 *--------------------------------------------------------------*/

bool Renewal::plot_write(const char *prefix , const char *title ,
                         const RenewalData *timev) const

{
  bool status;
  register int i , j , k , m , n;
  int nb_file , nb_dist , nb_histo , nb_time , inf , sup;
  double max , *scale;
  const FrequencyDistribution **phisto;
  const Distribution **pdist;
  ostringstream *data_file_name;


  // ecriture des fichiers de donnees

  nb_file = 2;
  if ((!timev) || (!(timev->length_bias))) {
    nb_file++;
  }

  nb_dist = PLOT_NEVENT_TIME;
  nb_histo = 0;
  if ((timev) && (timev->within) && (timev->backward) && (timev->forward)) {
    if (timev->inter_event) {
      nb_dist++;
      nb_histo++;
    }
    if (timev->within->nb_element > 0) {
      nb_dist++;
      nb_histo++;
    }
    if (timev->length_bias) {
      nb_dist += 2;
      nb_histo++;
    }

    nb_dist += 4;
    nb_histo += 2;
  }

  if ((!timev) || (!(timev->length_bias))) {
    nb_dist += 3;
  }

  nb_time = 0;
  for (i = time->offset;i < time->nb_value;i++) {
    if (time->mass[i] > 0.) {
      nb_time++;
    }
  }

  if (((time->variance == 0.) || (timev)) && (nb_time <= PLOT_NB_TIME)) {
    nb_dist += nb_time;
    if (timev) {
      nb_histo += nb_time;
    }
  }

  if (time->variance > 0.) {
    nb_dist += 1;
    if (timev) {
      nb_histo += 2;
    }

    if (nb_time <= PLOT_NB_TIME) {
      nb_file += nb_time;
    }
  }

  data_file_name = new ostringstream[nb_file];

  data_file_name[0] << prefix << 0 << ".dat";

  pdist = new const Distribution*[nb_dist];
  scale = new double[nb_dist];
  if (timev) {
    phisto = new const FrequencyDistribution*[nb_histo];
  }

  nb_histo = 0;
  nb_dist = 0;

  if ((timev) && (timev->within) && (timev->backward) && (timev->forward)) {
    if (timev->inter_event) {
      phisto[nb_histo++] = timev->inter_event;
      pdist[nb_dist] = inter_event;
      scale[nb_dist++] = timev->inter_event->nb_element;
    }

    if (timev->within->nb_element > 0) {
      phisto[nb_histo++] = timev->within;
      pdist[nb_dist] = inter_event;
      scale[nb_dist++] = timev->within->nb_element;
    }

    if (timev->length_bias) {
      phisto[nb_histo++] = timev->length_bias;
      pdist[nb_dist] = length_bias;
      scale[nb_dist++] = timev->length_bias->nb_element;
      pdist[nb_dist] = inter_event;
      scale[nb_dist++] = timev->length_bias->nb_element;
    }

    phisto[nb_histo++] = timev->backward;
    pdist[nb_dist] = backward;
    scale[nb_dist++] = timev->backward->nb_element;
    pdist[nb_dist] = inter_event;
    scale[nb_dist++] = timev->backward->nb_element;

    phisto[nb_histo++] = timev->forward;
    pdist[nb_dist] = forward;
    scale[nb_dist++] = timev->forward->nb_element;
    pdist[nb_dist] = inter_event;
    scale[nb_dist++] = timev->forward->nb_element;
  }

  if ((!timev) || (!(timev->length_bias))) {
    pdist[nb_dist] = inter_event;
    scale[nb_dist++] = 1.;
    pdist[nb_dist] = length_bias;
    scale[nb_dist++] = 1.;
    pdist[nb_dist] = forward;
    scale[nb_dist++] = 1.;
  }

  inf = (timev ? timev->mixture->offset : mixture->offset);
  if (inf < 1) {
    inf = 1;
  }

  sup = (timev ? timev->mixture->nb_value : mixture->nb_value);
  if (sup - inf > PLOT_NEVENT_TIME) {
    sup = inf + PLOT_NEVENT_TIME;
  }

  for (i = inf;i < sup;i++) {
    pdist[nb_dist] = nevent_time[i];
    scale[nb_dist++] = 1.;
  }

  if ((time->variance > 0.) && (timev)) {
    phisto[nb_histo++] = timev->htime;
  }

  if (((time->variance == 0.) || (timev)) && (nb_time <= PLOT_NB_TIME)) {
    for (i = time->offset;i < time->nb_value;i++) {
      if (time->mass[i] > 0.) {
        pdist[nb_dist] = nb_event[i];

        if (timev) {
          phisto[nb_histo++] = timev->hnb_event[i];
          scale[nb_dist++] = timev->hnb_event[i]->nb_element;
        }
        else {
          scale[nb_dist++] = 1.;
        }
      }
    }
  }

  if (time->variance > 0.) {
    pdist[nb_dist] = mixture;

    if (timev) {
      phisto[nb_histo++] = timev->mixture;
      scale[nb_dist++] = timev->mixture->nb_element;
    }
    else {
      scale[nb_dist++] = 1.;
    }
  }

  status = plot_print((data_file_name[0].str()).c_str() , nb_dist , pdist ,
                      scale , NULL , nb_histo , phisto);

  if (status) {
    i = 1;
    if ((!timev) || (!(timev->length_bias))) {
      data_file_name[i] << prefix << i << ".dat";
      backward->plot_print((data_file_name[i++].str()).c_str());
    }

    if ((time->variance > 0.) && (nb_time <= PLOT_NB_TIME)) {
      for (j = time->offset;j < time->nb_value;j++) {
        if (time->mass[j] > 0.) {
          data_file_name[i] << prefix << i << ".dat";

          pdist[0] = nb_event[j];

          if (timev) {
            scale[0] = time->mass[j] * timev->nb_element;
          }
          else {
            scale[0] = time->mass[j];
          }

          plot_print((data_file_name[i++].str()).c_str() , 1 , pdist ,
                     scale , NULL , 0 , NULL);
        }
      }
    }

    data_file_name[nb_file - 1] << prefix << nb_file - 1 << ".dat";
    index_event->plot_print((data_file_name[nb_file - 1].str()).c_str() , index_event->length ,
                            (((timev) && (timev->index_event)) ? timev->index_event : NULL));

    // ecriture du fichier de commandes et du fichier d'impression

    for (i = 0;i < 2;i++) {
      j = 1;
      k = nb_histo + 1;

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

      if ((timev) && (timev->within) && (timev->backward) && (timev->forward)) {
        if (inter_event->nb_value - 1 < TIC_THRESHOLD) {
          out_file << "set xtics 0,1" << endl;
        }

        if (timev->inter_event) {
          out_file << "plot [0:" << inter_event->nb_value - 1 << "] [0:"
                   << (int)(MAX(timev->inter_event->max , inter_event->max * timev->inter_event->nb_element) * YSCALE) + 1
                   << "] \"" << label((data_file_name[0].str()).c_str()) << "\" using " << j++
                   << " title \"" << SEQ_label[SEQL_INTER_EVENT] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
                   << "\" with impulses,\\" << endl;
          out_file << "\"" << label((data_file_name[0].str()).c_str()) << "\" using " << k++
                   << " title \"" << SEQ_label[SEQL_INTER_EVENT] << " " << STAT_label[STATL_DISTRIBUTION];
          inter_event->plot_title_print(out_file);
          out_file << "\" with linespoints" << endl;

          if (i == 0) {
            out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
          }
          out_file << endl;
        }

        if (timev->within->nb_element > 0) {
          out_file << "plot [0:" << inter_event->nb_value - 1 << "] [0:"
                   << (int)(MAX(timev->within->max , inter_event->max * timev->within->nb_element) * YSCALE) + 1
                   << "] \"" << label((data_file_name[0].str()).c_str()) << "\" using " << j++
                   << " title \"" << SEQ_label[SEQL_OBSERVATION_INTER_EVENT] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
                   << "\" with impulses,\\" << endl;
          out_file << "\"" << label((data_file_name[0].str()).c_str()) << "\" using " << k++
                   << " title \"" << SEQ_label[SEQL_INTER_EVENT] << " " << STAT_label[STATL_DISTRIBUTION];
          inter_event->plot_title_print(out_file);
          out_file << "\" with linespoints" << endl;

          if (i == 0) {
            out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
          }
          out_file << endl;
        }

        if (timev->length_bias) {
          max = MAX(length_bias->max , inter_event->max);

          out_file << "plot [0:" << inter_event->nb_value - 1 << "] [0:"
                   << (int)(MAX(timev->length_bias->max , max * timev->length_bias->nb_element) * YSCALE) + 1
                   << "] \"" << label((data_file_name[0].str()).c_str()) << "\" using " << j++
                   << " title \"" << SEQ_label[SEQL_LENGTH_BIASED] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
                   << "\" with impulses,\\" << endl;
          out_file << "\"" << label((data_file_name[0].str()).c_str()) << "\" using " << k++
                   << " title \"" << SEQ_label[SEQL_LENGTH_BIASED] << " " << STAT_label[STATL_DISTRIBUTION]
                   << "\" with linespoints,\\" << endl;
          out_file << "\"" << label((data_file_name[0].str()).c_str()) << "\" using " << k++
                   << " title \"" << SEQ_label[SEQL_INTER_EVENT] << " " << STAT_label[STATL_DISTRIBUTION];
          inter_event->plot_title_print(out_file);
          out_file << "\" with linespoints" << endl;

          if (i == 0) {
            out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
          }
          out_file << endl;
        }

        max = MAX(backward->max , inter_event->max);

        out_file << "plot [0:" << inter_event->nb_value - 1 << "] [0:"
                 << (int)(MAX(timev->backward->max , max * timev->backward->nb_element) * YSCALE) + 1
                 << "] \"" << label((data_file_name[0].str()).c_str()) << "\" using " << j++
                 << " title \"" << SEQ_label[SEQL_BACKWARD] << " " << SEQ_label[SEQL_RECURRENCE_TIME]
                 << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\" with impulses,\\" << endl;
        out_file << "\"" << label((data_file_name[0].str()).c_str()) << "\" using " << k++
                 << " title \""<< SEQ_label[SEQL_BACKWARD] << " " << SEQ_label[SEQL_RECURRENCE_TIME]
                 << " " << STAT_label[STATL_DISTRIBUTION] << "\" with linespoints,\\" << endl;
        out_file << "\"" << label((data_file_name[0].str()).c_str()) << "\" using " << k++
                 << " title \"" << SEQ_label[SEQL_INTER_EVENT] << " " << STAT_label[STATL_DISTRIBUTION];
        inter_event->plot_title_print(out_file);
        out_file << "\" with linespoints" << endl;

        if (i == 0) {
          out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
        }
        out_file << endl;

        max = MAX(forward->max , inter_event->max);

        out_file << "plot [0:" << inter_event->nb_value - 1 << "] [0:"
                 << (int)(MAX(timev->forward->max , max * timev->forward->nb_element) * YSCALE) + 1
                 << "] \"" << label((data_file_name[0].str()).c_str()) << "\" using " << j++
                 << " title \"" << SEQ_label[SEQL_FORWARD] << " " << SEQ_label[SEQL_RECURRENCE_TIME]
                 << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\" with impulses,\\" << endl;
        out_file << "\"" << label((data_file_name[0].str()).c_str()) << "\" using " << k++
                 << " title \"" << SEQ_label[SEQL_FORWARD] << " " << SEQ_label[SEQL_RECURRENCE_TIME]
                 << " " << STAT_label[STATL_DISTRIBUTION] << "\" with linespoints,\\" << endl;
        out_file << "\"" << label((data_file_name[0].str()).c_str()) << "\" using " << k++
                 << " title \"" << SEQ_label[SEQL_INTER_EVENT] << " " << STAT_label[STATL_DISTRIBUTION];
        inter_event->plot_title_print(out_file);
        out_file << "\" with linespoints" << endl;

        if (inter_event->nb_value - 1 < TIC_THRESHOLD) {
          out_file << "set xtics autofreq" << endl;
        }

        if (!(timev->length_bias)) {
          if (i == 0) {
            out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
          }
          out_file << endl;
        }
      }

      if ((!timev) || (!(timev->length_bias))) {
        if (inter_event->nb_value - 1 < TIC_THRESHOLD) {
          out_file << "set xtics 0,1" << endl;
        }

        max = inter_event->max;
        if (length_bias->max > max) {
          max = length_bias->max;
        }
        if (backward->max > max) {
          max = backward->max;
        }

        out_file << "plot [0:" << inter_event->nb_value - 1 << "] [0:"
                 << MIN(max * YSCALE , 1.) << "] \""
                 << label((data_file_name[0].str()).c_str()) << "\" using " << k++ << " title \""
                 << SEQ_label[SEQL_INTER_EVENT] << " " << STAT_label[STATL_DISTRIBUTION];
        inter_event->plot_title_print(out_file);
        out_file << "\" with linespoints,\\" << endl;
        out_file << "\"" << label((data_file_name[0].str()).c_str()) << "\" using " << k++ << " title \""
                 << SEQ_label[SEQL_LENGTH_BIASED] << " " << STAT_label[STATL_DISTRIBUTION]
                 << "\" with linespoints,\\" << endl;
        out_file << "\"" << label((data_file_name[1].str()).c_str()) << "\" title \""
                 << SEQ_label[SEQL_BACKWARD] << " " << SEQ_label[SEQL_RECURRENCE_TIME]
                 << " " << STAT_label[STATL_DISTRIBUTION] << "\" with linespoints,\\" << endl;
        out_file << "\"" << label((data_file_name[0].str()).c_str()) << "\" using " << k++ << " title \""
                 << SEQ_label[SEQL_FORWARD] << " " << SEQ_label[SEQL_RECURRENCE_TIME]
                 << " " << STAT_label[STATL_DISTRIBUTION] << "\" with linespoints" << endl;

        if (inter_event->nb_value - 1 < TIC_THRESHOLD) {
          out_file << "set xtics autofreq" << endl;
        }
      }

      if (inf < sup) {
        if (i == 0) {
          out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
        }
        out_file << endl;

        if (nevent_time[sup - 1]->nb_value - 1 < TIC_THRESHOLD) {
          out_file << "set xtics 0,1" << endl;
        }

        out_file << "plot [0:" << nevent_time[sup - 1]->nb_value - 1 << "] [0:"
                 << MIN(nevent_time[inf]->max * YSCALE , 1.) << "] ";
        for (m = inf;m < sup;m++) {
          out_file << "\"" << label((data_file_name[0].str()).c_str()) << "\" using " << k++
                   << " title \"" << SEQ_label[SEQL_TIME_UP] << " " << m << " "
                   << STAT_label[STATL_DISTRIBUTION] << "\" with linespoints";
          if (m < sup - 1) {
            out_file << ",\\";
          }
          out_file << endl;
        }

        if (nevent_time[sup - 1]->nb_value - 1 < TIC_THRESHOLD) {
          out_file << "set xtics autofreq" << endl;
        }
      }

      if ((time->variance > 0.) && (timev)) {
        if (i == 0) {
          out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
        }
        out_file << endl;

        if (timev->htime->nb_value - 1 < TIC_THRESHOLD) {
          out_file << "set xtics 0,1" << endl;
        }
        if ((int)(timev->htime->max * YSCALE) + 1 < TIC_THRESHOLD) {
          out_file << "set ytics 0,1" << endl;
        }

        out_file << "plot [0:" << timev->htime->nb_value - 1 << "] [0:"
                 << (int)(timev->htime->max * YSCALE) + 1 << "] \""
                 << label((data_file_name[0].str()).c_str()) << "\" using " << j++ << " title \""
                 << SEQ_label[SEQL_OBSERVATION_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
                 << "\" with impulses" << endl;

        if (timev->htime->nb_value - 1 < TIC_THRESHOLD) {
          out_file << "set xtics autofreq" << endl;
        }
        if ((int)(timev->htime->max * YSCALE) + 1 < TIC_THRESHOLD) {
          out_file << "set ytics autofreq" << endl;
        }
      }

      if (((time->variance == 0.) || (timev)) && (nb_time <= PLOT_NB_TIME)) {
        for (m = time->offset;m < time->nb_value;m++) {
          if (time->mass[m] > 0.) {
            if (i == 0) {
              out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
            }
            out_file << endl;

            if (nb_event[m]->nb_value - 1 < TIC_THRESHOLD) {
              out_file << "set xtics 0,1" << endl;
            }

            if (timev) {
              out_file << "plot [0:" << nb_event[m]->nb_value - 1 << "] [0:"
                       << (int)(MAX(timev->hnb_event[m]->max , nb_event[m]->max * timev->hnb_event[m]->nb_element) * YSCALE) + 1
                       << "] \"" << label((data_file_name[0].str()).c_str()) << "\" using " << j++ << " title \""
                       << SEQ_label[SEQL_NB_EVENT] << " " << SEQ_label[SEQL_DURING] << " " << m << " "
                       << SEQ_label[SEQL_TIME_UNIT] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\" with impulses,\\" << endl;
            }

            else {
              out_file << "plot [0:" << nb_event[m]->nb_value - 1 << "] [0:"
                       << MIN(nb_event[m]->max * YSCALE , 1.) << "] ";
            }

            out_file << "\"" << label((data_file_name[0].str()).c_str()) << "\" using " << k++ << " title \""
                     << SEQ_label[SEQL_NB_EVENT] << " " << SEQ_label[SEQL_DURING] << " " << m << " "
                     << SEQ_label[SEQL_TIME_UNIT] << " " << STAT_label[STATL_DISTRIBUTION] << "\" with linespoints" << endl;

            if (nb_event[m]->nb_value - 1 < TIC_THRESHOLD) {
              out_file << "set xtics autofreq" << endl;
            }
          }
        }
      }

      if (time->variance > 0.) {
        if (i == 0) {
          out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
        }
        out_file << endl;

        if (mixture->nb_value - 1 < TIC_THRESHOLD) {
          out_file << "set xtics 0,1" << endl;
        }

        if (timev) {
          out_file << "plot [0:" << mixture->nb_value - 1 << "] [0:"
                   << (int)(MAX(timev->mixture->max , mixture->max * timev->mixture->nb_element) * YSCALE) + 1
                   << "] \"" << label((data_file_name[0].str()).c_str()) << "\" using " << j++ << " title \""
                   << SEQ_label[SEQL_NB_EVENT] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
                   << "\" with impulses,\\" << endl;
        }

        else {
          out_file << "plot [0:" << mixture->nb_value - 1 << "] [0:"
                   << MIN(mixture->max * YSCALE , 1.) << "] ";
        }

        out_file << "\"" << label((data_file_name[0].str()).c_str()) << "\" using " << k++
                 << " title \"" << SEQ_label[SEQL_NB_EVENT_MIXTURE] << "\" with linespoints";

        if (nb_time <= PLOT_NB_TIME) {
          m = (((!timev) || (!(timev->length_bias))) ? 2 : 1);
          for (n = time->offset;n < time->nb_value;n++) {
            if (time->mass[n] > 0.) {
              out_file << ",\\" << endl;
              out_file << "\"" << label((data_file_name[m++].str()).c_str()) << "\" title \""
                       << SEQ_label[SEQL_NB_EVENT] << " " << SEQ_label[SEQL_DURING] << " " << n << " "
                       << SEQ_label[SEQL_TIME_UNIT] << " " << STAT_label[STATL_DISTRIBUTION] << "\" with linespoints";
            }
          }
        }
        out_file << endl;

        if (mixture->nb_value - 1 < TIC_THRESHOLD) {
          out_file << "set xtics autofreq" << endl;
        }
      }

      if (i == 0) {
        out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
      }
      out_file << endl;

      if (index_event->length - 1 < TIC_THRESHOLD) {
        out_file << "set xtics 0,1" << endl;
      }

      out_file << "plot [" << index_event->offset << ":" << index_event->length - 1 << "] [0:1] ";
      if ((timev) && (timev->index_event)) {
        out_file << "\"" << label((data_file_name[nb_file - 1].str()).c_str()) << "\" using "
                 << 3 << " title \"" << SEQ_label[SEQL_OBSERVED] << " "
                 << SEQ_label[SEQL_NO_EVENT_PROBABILITY] << " \" with linespoints,\\" << endl;
      }
      out_file << "\"" << label((data_file_name[nb_file - 1].str()).c_str()) << "\" using "
               << 1 << " title \"" << SEQ_label[SEQL_THEORETICAL] << " "
               << SEQ_label[SEQL_NO_EVENT_PROBABILITY] << " \" with linespoints,\\" << endl;
      if ((timev) && (timev->index_event)) {
        out_file << "\"" << label((data_file_name[nb_file - 1].str()).c_str()) << "\" using "
                 << 4 << " title \"" << SEQ_label[SEQL_OBSERVED] << " "
                 << SEQ_label[SEQL_EVENT_PROBABILITY] << " \" with linespoints,\\" << endl;
      }
      out_file << "\"" << label((data_file_name[nb_file - 1].str()).c_str()) << "\" using "
               << 2 << " title \"" << SEQ_label[SEQL_THEORETICAL] << " "
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

  delete [] pdist;
  delete [] scale;
  if (timev) {
    delete [] phisto;
  }

  delete [] data_file_name;

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Sortie Gnuplot d'un objet Renewal.
 *
 *  arguments : reference sur un objet StatError, prefixe des fichiers,
 *              titre des figures.
 *
 *--------------------------------------------------------------*/

bool Renewal::plot_write(StatError &error , const char *prefix ,
                         const char *title) const

{
  bool status = plot_write(prefix , title , renewal_data);

  error.init();

  if (!status) {
    error.update(STAT_error[STATR_FILE_PREFIX]);
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Sortie graphique d'un processus de renouvellement.
 *
 *  argument : pointeur sur un objet RenewalData.
 *
 *--------------------------------------------------------------*/

MultiPlotSet* Renewal::get_plotable(const RenewalData *timev) const

{
  register int i , j , k;
  int nb_plot_set , nb_time , inf , sup , scale;
  double max;
  ostringstream legend , title;
  MultiPlotSet *plot_set;


  nb_plot_set = 1;
  if ((timev) && (timev->within) && (timev->backward) && (timev->forward)) {
    nb_plot_set += 2;
    if (timev->inter_event) {
      nb_plot_set++;
    }
    if (timev->within->nb_element > 0) {
      nb_plot_set++;
    }
    if (timev->length_bias) {
      nb_plot_set++;
    }
  }

  if ((!timev) || (!(timev->length_bias))) {
    nb_plot_set++;
  }

  inf = (timev ? timev->mixture->offset : mixture->offset);
  if (inf < 1) {
    inf = 1;
  }

  sup = (timev ? timev->mixture->nb_value : mixture->nb_value);
  if (sup - inf > PLOT_NEVENT_TIME) {
    sup = inf + PLOT_NEVENT_TIME;
  }

  if (inf < sup) {
    nb_plot_set++;
  }

  if ((time->variance > 0.) && (timev)) {
    nb_plot_set++;
  }

  nb_time = 0;
  for (i = time->offset;i < time->nb_value;i++) {
    if (time->mass[i] > 0.) {
      nb_time++;
    }
  }

  if (nb_time <= PLOT_NB_TIME) {
    if (!timev) {
      nb_plot_set++;
    }
    else {
      nb_plot_set += 2 * nb_time;
    }
  }

  if (time->variance > 0.) {
    nb_plot_set++;
  }

  plot_set = new MultiPlotSet(nb_plot_set);
  MultiPlotSet &plot = *plot_set;

  plot.border = "15 lw 0";

  i = 0;
  if ((timev) && (timev->within) && (timev->backward) && (timev->forward)) {
    if (timev->inter_event) {

      // vue : loi inter-evenement ajustee

      plot[i].xrange = Range(0 , inter_event->nb_value - 1);
      plot[i].yrange = Range(0 , ceil(MAX(timev->inter_event->max ,
                                          inter_event->max * timev->inter_event->nb_element) * YSCALE));

      if (inter_event->nb_value - 1 < TIC_THRESHOLD) {
        plot[i].xtics = 1;
      }

      plot[i].resize(2);

      legend.str("");
      legend << SEQ_label[SEQL_INTER_EVENT] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
      plot[i][0].legend = legend.str();

      plot[i][0].style = "impulses";

      timev->inter_event->plotable_frequency_write(plot[i][0]);

      legend.str("");
      legend << SEQ_label[SEQL_INTER_EVENT] << " " << STAT_label[STATL_DISTRIBUTION];
      inter_event->plot_title_print(legend);
      plot[i][1].legend = legend.str();

      plot[i][1].style = "linespoints";

      inter_event->plotable_mass_write(plot[i][1] , timev->inter_event->nb_element);
      i++;
    }

    if (timev->within->nb_element > 0) {

      // vue : loi empirique de l'intervalle de temps entre 2 evenements a l'interieur
      // de la periode d'observation ajustee par la loi inter-evenement

      plot[i].xrange = Range(0 , inter_event->nb_value - 1);
      plot[i].yrange = Range(0 , ceil(MAX(timev->within->max ,
                                          inter_event->max * timev->within->nb_element) * YSCALE));

      if (inter_event->nb_value - 1 < TIC_THRESHOLD) {
        plot[i].xtics = 1;
      }

      plot[i].resize(2);

      legend.str("");
      legend << SEQ_label[SEQL_OBSERVATION_INTER_EVENT] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
      plot[i][0].legend = legend.str();

      plot[i][0].style = "impulses";

      timev->within->plotable_frequency_write(plot[i][0]);

      legend.str("");
      legend << SEQ_label[SEQL_INTER_EVENT] << " " << STAT_label[STATL_DISTRIBUTION];
      inter_event->plot_title_print(legend);
      plot[i][1].legend = legend.str();

      plot[i][1].style = "linespoints";

      inter_event->plotable_mass_write(plot[i][1] , timev->within->nb_element);
      i++;
    }

    if (timev->length_bias) {

      // vue : loi biaisee par la longueur ajustee

      plot[i].xrange = Range(0 , inter_event->nb_value - 1);

      max = MAX(length_bias->max , inter_event->max);
      plot[i].yrange = Range(0 , ceil(MAX(timev->length_bias->max ,
                                          max * timev->length_bias->nb_element) * YSCALE));

      if (inter_event->nb_value - 1 < TIC_THRESHOLD) {
        plot[i].xtics = 1;
      }

      plot[i].resize(3);

      legend.str("");
      legend << SEQ_label[SEQL_LENGTH_BIASED] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
      plot[i][0].legend = legend.str();

      plot[i][0].style = "impulses";

      timev->length_bias->plotable_frequency_write(plot[i][0]);

      legend.str("");
      legend << SEQ_label[SEQL_LENGTH_BIASED] << " " << STAT_label[STATL_DISTRIBUTION];
      plot[i][1].legend = legend.str();

      plot[i][1].style = "linespoints";

      length_bias->plotable_mass_write(plot[i][1] , timev->length_bias->nb_element);

      legend.str("");
      legend << SEQ_label[SEQL_INTER_EVENT] << " " << STAT_label[STATL_DISTRIBUTION];
      inter_event->plot_title_print(legend);
      plot[i][2].legend = legend.str();

      plot[i][2].style = "linespoints";

      inter_event->plotable_mass_write(plot[i][2] , timev->length_bias->nb_element);
      i++;
    }

    // vue : loi de l'intervalle de temps apres le dernier evenement ajustee

    plot[i].xrange = Range(0 , inter_event->nb_value - 1);

    max = MAX(backward->max , inter_event->max);
    plot[i].yrange = Range(0 , ceil(MAX(timev->backward->max ,
                                        max * timev->backward->nb_element) * YSCALE));

    if (inter_event->nb_value - 1 < TIC_THRESHOLD) {
      plot[i].xtics = 1;
    }

    plot[i].resize(3);

    legend.str("");
    legend << SEQ_label[SEQL_BACKWARD] << " " << SEQ_label[SEQL_RECURRENCE_TIME] << " "
           << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
    plot[i][0].legend = legend.str();

    plot[i][0].style = "impulses";

    timev->backward->plotable_frequency_write(plot[i][0]);

    legend.str("");
    legend << SEQ_label[SEQL_BACKWARD] << " " << SEQ_label[SEQL_RECURRENCE_TIME] << " "
           << STAT_label[STATL_DISTRIBUTION];
    plot[i][1].legend = legend.str();

    plot[i][1].style = "linespoints";

    backward->plotable_mass_write(plot[i][1] , timev->backward->nb_element);

    legend.str("");
    legend << SEQ_label[SEQL_INTER_EVENT] << " " << STAT_label[STATL_DISTRIBUTION];
    inter_event->plot_title_print(legend);
    plot[i][2].legend = legend.str();

    plot[i][2].style = "linespoints";

    inter_event->plotable_mass_write(plot[i][2] , timev->backward->nb_element);
    i++;

    // vue : loi de l'intervalle de temps residuel ajustee

    plot[i].xrange = Range(0 , inter_event->nb_value - 1);

    max = MAX(forward->max , inter_event->max);
    plot[i].yrange = Range(0 , ceil(MAX(timev->forward->max ,
                                        max * timev->forward->nb_element) * YSCALE));

    if (inter_event->nb_value - 1 < TIC_THRESHOLD) {
      plot[i].xtics = 1;
    }

    plot[i].resize(3);

    legend.str("");
    legend << SEQ_label[SEQL_FORWARD] << " " << SEQ_label[SEQL_RECURRENCE_TIME] << " "
           << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
    plot[i][0].legend = legend.str();

    plot[i][0].style = "impulses";

    timev->forward->plotable_frequency_write(plot[i][0]);

    legend.str("");
    legend << SEQ_label[SEQL_FORWARD] << " " << SEQ_label[SEQL_RECURRENCE_TIME] << " "
           << STAT_label[STATL_DISTRIBUTION];
    plot[i][1].legend = legend.str();

    plot[i][1].style = "linespoints";

    forward->plotable_mass_write(plot[i][1] , timev->forward->nb_element);

    legend.str("");
    legend << SEQ_label[SEQL_INTER_EVENT] << " " << STAT_label[STATL_DISTRIBUTION];
    inter_event->plot_title_print(legend);
    plot[i][2].legend = legend.str();

    plot[i][2].style = "linespoints";

    inter_event->plotable_mass_write(plot[i][2] , timev->forward->nb_element);
    i++;
  }

  if ((!timev) || (!(timev->length_bias))) {

    // vue : loi inter-evenement, loi biaisee par la longueur, loi de l'intervalle de temps
    // apres le dernier evenement et loi de l'intervalle de temps residuel

    plot[i].xrange = Range(0 , inter_event->nb_value - 1);

    max = inter_event->max;
    if (length_bias->max > max) {
      max = length_bias->max;
    }
    if (backward->max > max) {
      max = backward->max;
    }
    plot[i].yrange = Range(0. , MIN(max * YSCALE , 1.));

    if (inter_event->nb_value - 1 < TIC_THRESHOLD) {
      plot[i].xtics = 1;
    }

    plot[i].resize(4);

    legend.str("");
    legend << SEQ_label[SEQL_INTER_EVENT] << " " << STAT_label[STATL_DISTRIBUTION];
    inter_event->plot_title_print(legend);
    plot[i][0].legend = legend.str();

    plot[i][0].style = "linespoints";

    inter_event->plotable_mass_write(plot[i][0]);

    legend.str("");
    legend << SEQ_label[SEQL_LENGTH_BIASED] << " " << STAT_label[STATL_DISTRIBUTION];
    plot[i][1].legend = legend.str();

    plot[i][1].style = "linespoints";

    length_bias->plotable_mass_write(plot[i][1]);

    legend.str("");
    legend << SEQ_label[SEQL_BACKWARD] << " " << SEQ_label[SEQL_RECURRENCE_TIME] << " "
           << STAT_label[STATL_DISTRIBUTION];
    plot[i][2].legend = legend.str();

    plot[i][2].style = "linespoints";

    backward->plotable_mass_write(plot[i][2]);

    legend.str("");
    legend << SEQ_label[SEQL_FORWARD] << " " << SEQ_label[SEQL_RECURRENCE_TIME] << " "
           << STAT_label[STATL_DISTRIBUTION];
    plot[i][3].legend = legend.str();

    plot[i][3].style = "linespoints";

    forward->plotable_mass_write(plot[i][3]);
    i++;
  }

  if (inf < sup) {

    // vue : lois du temps avant le neme evenement

    plot[i].xrange = Range(0 , nevent_time[sup - 1]->nb_value - 1);
    plot[i].yrange = Range(0. , MIN(nevent_time[inf]->max * YSCALE , 1.));

    if (nevent_time[sup - 1]->nb_value - 1 < TIC_THRESHOLD) {
      plot[i].xtics = 1;
    }

    plot[i].resize(sup - inf);

    for (j = inf;j < sup;j++) {
      legend.str("");
      legend << SEQ_label[SEQL_TIME_UP] << " " << j << " " << STAT_label[STATL_DISTRIBUTION];
      plot[i][j - inf].legend = legend.str();

      plot[i][j - inf].style = "linespoints";

      nevent_time[j]->plotable_mass_write(plot[i][j - inf]);
    }

    i++;
  }

  if ((time->variance > 0.) && (timev)) {

    // vue : loi empirique des temps d'observation

    plot[i].xrange = Range(0 , timev->htime->nb_value - 1);
    plot[i].yrange = Range(0 , ceil(timev->htime->max * YSCALE));

    if (timev->htime->nb_value - 1 < TIC_THRESHOLD) {
      plot[i].xtics = 1;
    }
    if (ceil(timev->htime->max * YSCALE) < TIC_THRESHOLD) {
      plot[i].ytics = 1;
    }

    plot[i].resize(1);

    legend.str("");
    legend << SEQ_label[SEQL_OBSERVATION_TIME] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
    plot[i][0].legend = legend.str();

    plot[i][0].style = "impulses";

    timev->htime->plotable_frequency_write(plot[i][0]);
    i++;
  }

  if (nb_time <= PLOT_NB_TIME) {
    if (!timev) {

      // vue : lois du nombre d'evenements

      plot[i].xrange = Range(0 , nb_event[time->nb_value - 1]->nb_value - 1);
      plot[i].yrange = Range(0. , MIN(nb_event[time->offset]->max * YSCALE , 1.));

      if (nb_event[time->nb_value - 1]->nb_value - 1 < TIC_THRESHOLD) {
        plot[i].xtics = 1;
      }

      plot[i].resize(nb_time);

      j = 0;
      for (k = time->offset;k < time->nb_value;k++) {
        if (time->mass[k] > 0.) {
          legend.str("");
          legend << SEQ_label[SEQL_NB_EVENT] << " " << SEQ_label[SEQL_DURING] << " " << k << " "
                 << SEQ_label[SEQL_TIME_UNIT] << " " << STAT_label[STATL_DISTRIBUTION];
          plot[i][j].legend = legend.str();

          plot[i][j].style = "linepoints";

          nb_event[k]->plotable_mass_write(plot[i][j]);

          j++;
        }
      }

      i++;
    }

    else {
      for (j = time->offset;j < time->nb_value;j++) {
        if (time->mass[j] > 0.) {

          // vue : loi du nombre d'evenements ajustee

          plot[i].xrange = Range(0 , nb_event[j]->nb_value - 1);
          plot[i].yrange = Range(0 , ceil(MAX(timev->hnb_event[j]->max ,
                                              nb_event[j]->max * timev->hnb_event[j]->nb_element) * YSCALE));

          if (nb_event[j]->nb_value - 1 < TIC_THRESHOLD) {
            plot[i].xtics = 1;
          }

          plot[i].resize(2);

          legend.str("");
          legend << SEQ_label[SEQL_NB_EVENT] << " " << SEQ_label[SEQL_DURING] << " " << j << " "
                 << SEQ_label[SEQL_TIME_UNIT] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
          plot[i][0].legend = legend.str();

          plot[i][0].style = "impulses";

          timev->hnb_event[j]->plotable_frequency_write(plot[i][0]);

          legend.str("");
          legend << SEQ_label[SEQL_NB_EVENT] << " " << SEQ_label[SEQL_DURING] << " " << j << " "
                 << SEQ_label[SEQL_TIME_UNIT] << " " << STAT_label[STATL_DISTRIBUTION];
          plot[i][1].legend = legend.str();

          plot[i][1].style = "linespoints";

          nb_event[j]->plotable_mass_write(plot[i][1] , timev->hnb_event[j]->nb_element);
          i++;

          // vue : fonctions de repartition des lois du nombre d'evenements

          title.str("");
          title << SEQ_label[SEQL_NB_EVENT] << " " << SEQ_label[SEQL_DURING] << " "
                << j << " " << SEQ_label[SEQL_TIME_UNIT];
          plot[i].title = title.str();

          plot[i].xrange = Range(0 , nb_event[j]->nb_value - 1);
          plot[i].yrange = Range(0 , 1.);

          if (nb_event[j]->nb_value - 1 < TIC_THRESHOLD) {
            plot[i].xtics = 1;
          }

          plot[i].resize(2);

          legend.str("");
          legend << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " "
                 << STAT_label[STATL_FUNCTION];
          plot[i][0].legend = legend.str();

          plot[i][0].style = "linespoints";

          timev->hnb_event[j]->plotable_cumul_write(plot[i][0]);

          legend.str("");
          legend << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_DISTRIBUTION] << " "
                 << STAT_label[STATL_FUNCTION];
          plot[i][1].legend = legend.str();

          plot[i][1].style = "linespoints";

          nb_event[j]->plotable_cumul_write(plot[i][1]);
          i++;
        }
      }
    }
  }

  if (time->variance > 0.) {
    if (timev) {
      plot[i].yrange = Range(0 , ceil(MAX(timev->mixture->max ,
                                          mixture->max * timev->mixture->nb_element) * YSCALE));
      plot[i].resize(2);

      legend.str("");
      legend << SEQ_label[SEQL_NB_EVENT] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
      plot[i][0].legend = legend.str();

      plot[i][0].style = "impulses";

      timev->mixture->plotable_frequency_write(plot[i][0]);

      scale = timev->mixture->nb_element;
      j = 1;
    }

    else {
      plot[i].yrange = Range(0. , MIN(mixture->max * YSCALE , 1.));
      plot[i].resize(1);

      scale = 1.;
      j = 0;
    }

    plot[i].xrange = Range(0 , mixture->nb_value - 1);

    if (mixture->nb_value - 1 < TIC_THRESHOLD) {
      plot[i].xtics = 1;
    }

    plot[i][j].legend = SEQ_label[SEQL_NB_EVENT_MIXTURE];

    plot[i][j].style = "linespoints";

    mixture->plotable_mass_write(plot[i][j] , scale);
    i++;
  }

  // vue : probabilites de non-evenement/evenement fonction du temps ajustees

  plot[i].xrange = Range(index_event->offset , index_event->length - 1);
  plot[i].yrange = Range(0. , 1.);

  if (index_event->length - 1 < TIC_THRESHOLD) {
    plot[i].xtics = 1;
  }

  if ((timev) && (timev->index_event)) {
    plot[i].resize(4);
  }
  else {
    plot[i].resize(2);
  }

  j = 0;
  if ((timev) && (timev->index_event)) {
    legend.str("");
    legend << SEQ_label[SEQL_OBSERVED] << " " << SEQ_label[SEQL_NO_EVENT_PROBABILITY];
    plot[i][j].legend = legend.str();

    plot[i][j].style = "linespoints";

    timev->index_event->plotable_write(0 , plot[i][j]);
    j++;
  }

  legend.str("");
  legend << SEQ_label[SEQL_THEORETICAL] << " " << SEQ_label[SEQL_NO_EVENT_PROBABILITY];
  plot[i][j].legend = legend.str();

  plot[i][j].style = "linespoints";

  index_event->plotable_write(0 , plot[i][j]);
  j++;

  if ((timev) && (timev->index_event)) {
    legend.str("");
    legend << SEQ_label[SEQL_OBSERVED] << " " << SEQ_label[SEQL_EVENT_PROBABILITY];
    plot[i][j].legend = legend.str();

    plot[i][j].style = "linespoints";

    timev->index_event->plotable_write(1 , plot[i][j]);
    j++;
  }

  legend.str("");
  legend << SEQ_label[SEQL_THEORETICAL] << " " << SEQ_label[SEQL_EVENT_PROBABILITY];
  plot[i][j].legend = legend.str();

  plot[i][j].style = "linespoints";

  index_event->plotable_write(1 , plot[i][j]);

  return plot_set;
}


/*--------------------------------------------------------------*
 *
 *  Sortie graphique d'un processus de renouvellement.
 *
 *--------------------------------------------------------------*/

MultiPlotSet* Renewal::get_plotable() const

{
  return get_plotable(renewal_data);
}


};  // namespace sequence_analysis
