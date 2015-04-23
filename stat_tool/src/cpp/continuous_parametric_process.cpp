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
 *       $Id: continuous_parametric_process.cpp 8176 2010-02-18 10:30:53Z guedon $
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

#include <boost/math/distributions/gamma.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/special_functions/bessel.hpp>

#include "tool/rw_tokenizer.h"
#include "tool/rw_cstring.h"
#include "tool/rw_locale.h"
#include "tool/config.h"

#include "stat_tools.h"
#include "markovian.h"
#include "stat_label.h"

using namespace std;
using namespace boost::math;


namespace stat_tool {


const static double bilateral_tail[7] = {0.05, 0.025, 0.01, 0.005, 0.0025, 0.001, 0.0005};
const static double posterior_threshold[7] = {0.25, 0.1, 0.05, 0.025, 0.01, 0.005, 0.0025};



/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe ContinuousParametricProcess.
 *
 *  argument : nombre d'etats.
 *
 *--------------------------------------------------------------*/

ContinuousParametricProcess::ContinuousParametricProcess(int inb_state)

{
  register int i;


  nb_state = inb_state;
  ident = I_DEFAULT;
  tied_location = false;
  tied_dispersion = false;
  unit = I_DEFAULT;

  if (nb_state > 0) {
    observation = new ContinuousParametric*[nb_state];
    for (i = 0;i < nb_state;i++) {
      observation[i] = NULL;
    }
  }

  else {
    observation = NULL;
  }

  weight = NULL;
  restoration_weight = NULL;
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe ContinuousParametricProcess.
 *
 *  arguments : nombre d'etats, lois d'observation.
 *
 *--------------------------------------------------------------*/

ContinuousParametricProcess::ContinuousParametricProcess(int inb_state , ContinuousParametric **pobservation)

{
  register int i;


  nb_state = inb_state;
  ident = pobservation[0]->ident;
  tied_location = false;
  tied_dispersion = false;
  unit = pobservation[0]->unit;

  observation = new ContinuousParametric*[nb_state];
  for (i = 0;i < nb_state;i++) {
    observation[i] = new ContinuousParametric(*pobservation[i]);
  }

  weight = NULL;
  restoration_weight = NULL;
}


/*--------------------------------------------------------------*
 *
 *  Copie d'un objet ContinuousParametricProcess.
 *
 *  argument : reference sur un objet ContinuousParametricProcess.
 *
 *--------------------------------------------------------------*/

void ContinuousParametricProcess::copy(const ContinuousParametricProcess &process)

{
  register int i;


  nb_state = process.nb_state;
  ident = process.ident;
  tied_location = process.tied_location;
  tied_dispersion = process.tied_dispersion;
  unit = process.unit;

  observation = new ContinuousParametric*[nb_state];
  for (i = 0;i < nb_state;i++) {
    observation[i] = new ContinuousParametric(*(process.observation[i]));
  }

  if (process.weight) {
    weight = new Distribution(*(process.weight));
  }
  else {
    weight = NULL;
  }

  if (process.restoration_weight) {
    restoration_weight = new Distribution(*(process.restoration_weight));
  }
  else {
    restoration_weight = NULL;
  }
}


/*--------------------------------------------------------------*
 *
 *  Destruction des champs d'un objet ContinuousParametricProcess.
 *
 *--------------------------------------------------------------*/

void ContinuousParametricProcess::remove()

{
  if (observation) {
    register int i;


    for (i = 0;i < nb_state;i++) {
      delete observation[i];
    }
    delete [] observation;

    observation = NULL;
  }

  delete weight;
  weight = NULL;

  delete restoration_weight;
  restoration_weight = NULL;
}


/*--------------------------------------------------------------*
 *
 *  Destructeur de la classe ContinuousParametricProcess.
 *
 *--------------------------------------------------------------*/

ContinuousParametricProcess::~ContinuousParametricProcess()

{
  remove();
}


/*--------------------------------------------------------------*
 *
 *  Operateur d'assignement de la classe ContinuousParametricProcess.
 *
 *  argument : reference sur un objet ContinuousParametricProcess.
 *
 *--------------------------------------------------------------*/

ContinuousParametricProcess& ContinuousParametricProcess::operator=(const ContinuousParametricProcess &process)

{
  if (&process != this) {
    remove();
    copy(process);
  }

  return *this;
}


/*--------------------------------------------------------------*
 *
 *  Analyse du format des lois d'observation.
 *
 *  arguments : reference sur un objet StatError, stream,
 *              reference sur l'indice de la ligne lue, nombre d'etats,
 *              type de modele, identificateur de la derniere loi dans la liste.
 *
 *--------------------------------------------------------------*/

ContinuousParametricProcess* continuous_observation_parsing(StatError &error , ifstream &in_file ,
                                                            int &line , int nb_state , int model ,
                                                            int last_ident)

{
  RWLocaleSnapshot locale("en");
  RWCString buffer , token;
  size_t position;
  bool status = true , lstatus;
  register int i , j;
  long index;
  ContinuousParametric **dist;
  ContinuousParametricProcess *process;


  process = NULL;

  dist = new ContinuousParametric*[nb_state];
  for (i = 0;i < nb_state;i++) {
    dist[i] = NULL;
  }

  for (i = 0;i < nb_state;i++) {
    while (buffer.readLine(in_file , false)) {
      line++;

#     ifdef DEBUG
      cout << line << "  " << buffer << endl;
#     endif

      position = buffer.first('#');
      if (position != RW_NPOS) {
        buffer.remove(position);
      }
      j = 0;

      RWCTokenizer next(buffer);

      while (!((token = next()).isNull())) {
        switch (j) {

        // test mot cle COMPONENT / STATE

        case 0 : {
          switch (model) {

          case MIXTURE : {
            if (token != STAT_word[STATW_COMPONENT]) {
              status = false;
              error.correction_update(STAT_parsing[STATP_KEY_WORD] ,
                                      STAT_word[STATW_COMPONENT] , line , j + 1);
            }
            break;
          }

          case HIDDEN_MARKOV : {
            if (token != STAT_word[STATW_STATE]) {
              status = false;
              error.correction_update(STAT_parsing[STATP_KEY_WORD] ,
                                      STAT_word[STATW_STATE] , line , j + 1);
            }
            break;
          }
          }
          break;
        }

        // test indice de la composante / de l'etat

        case 1 : {
          lstatus = locale.stringToNum(token , &index);
          if ((lstatus) && (index != i)) {
            lstatus = false;
          }

          if (!lstatus) {
            status = false;

            switch (model) {
            case MIXTURE :
              error.correction_update(STAT_parsing[STATP_COMPONENT_INDEX] , i , line , j + 1);
              break;
            case HIDDEN_MARKOV :
              error.correction_update(STAT_parsing[STATP_STATE_INDEX] , i , line , j + 1);
              break;
            }
          }
          break;
        }

        // test mot cle OBSERVATION_DISTRIBUTION

        case 2 : {
          if ((token != STAT_word[STATW_OBSERVATION_DISTRIBUTION]) &&
              (token != STAT_word[STATW_OBSERVATION_MODEL])) {
            status = false;
            error.correction_update(STAT_parsing[STATP_KEY_WORD] ,
                                    STAT_word[STATW_OBSERVATION_DISTRIBUTION] , line , j + 1);
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

        switch (model) {
        case MIXTURE :
          dist[i] = continuous_parametric_parsing(error , in_file , line , VON_MISES);
          break;
        case HIDDEN_MARKOV :
          dist[i] = continuous_parametric_parsing(error , in_file , line , last_ident);
          break;
        }
        if (!dist[i]) {
          status = false;
        }

        break;
      }
    }
  }

  if (status) {
    for (i = 1;i < nb_state;i++) {
      if (dist[i]->ident != dist[0]->ident) {
        status = false;
        error.correction_update(STAT_parsing[STATP_OBSERVATION_DISTRIBUTION_TYPE] ,
                                STAT_continuous_distribution_word[dist[0]->ident]);
      }
    }
  }

  if (status) {
    if (dist[0]->ident == VON_MISES) {
      for (i = 0;i < nb_state;i++) {
        if (dist[i]->unit == DEGREE) {
          for (j = 0;j < nb_state;j++) {
            dist[j]->unit = DEGREE;
          }
          break;
        }
      }
    }

    process = new ContinuousParametricProcess(nb_state , dist);
  }

  for (i = 0;i < nb_state;i++) {
    delete dist[i];
  }
  delete [] dist;

  return process;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet ContinuousParametricProcess.
 *
 *  arguments : stream, pointeurs sur les histogrammes d'observation ou
 *              les lois d'observation empiriques et sur l'histogramme marginale ou
 *              la loi marginale empirique, flag niveau de detail,
 *              flag fichier, type de modele.
 *
 *--------------------------------------------------------------*/

ostream& ContinuousParametricProcess::ascii_print(ostream &os , Histogram **observation_histogram ,
                                                  FrequencyDistribution **observation_distribution ,
                                                  Histogram *marginal_histogram ,
                                                  FrequencyDistribution *marginal_distribution ,
                                                  bool exhaustive , bool file_flag , int model) const

{
  register int i , j , k;
  int nb_step , nb_element , buff , width[5];
  long old_adjust;
  double step , value , min_value , max_value , mass , *observation_cumul ,
         **frequency , **cumul;
  gamma_distribution<double> **gamma_dist;
  normal **gaussian_dist;


  if ((marginal_histogram) || (marginal_distribution)) {
    if (marginal_histogram) {
      step = marginal_histogram->step;
    }
    else {
      step = marginal_distribution->min_interval_computation();
    }

    if ((ident == GAMMA) || (ident == ZERO_INFLATED_GAMMA)) {
      if (marginal_histogram) {
        max_value = marginal_histogram->max_value;

/*        switch (marginal_histogram->type) {
        case INT_VALUE :
          max_value = marginal_histogram->max_value;
          break;
        case REAL_VALUE :
          max_value = marginal_histogram->max_value + marginal_histogram->step;
          break;
        } */
      }

      else {
        max_value = marginal_distribution->nb_value - 1;
      }
    }

    else if ((ident == GAUSSIAN) || (ident == VON_MISES)) {
      if (marginal_histogram) {
        min_value = marginal_histogram->min_value;
        max_value = marginal_histogram->max_value;

/*        switch (marginal_histogram->type) {
        case INT_VALUE :
          min_value = marginal_histogram->min_value;
          max_value = marginal_histogram->max_value;
          break;
        case REAL_VALUE :
          min_value = marginal_histogram->min_value - marginal_histogram->step / 2;
          max_value = marginal_histogram->max_value + marginal_histogram->step / 2;
          break;
        } */
      }

      else {
        min_value = marginal_distribution->offset;
        max_value = marginal_distribution->nb_value - 1;
      }
    }

    frequency = new double*[nb_state + 1];
    cumul = new double*[nb_state + 1];

    if ((ident == GAMMA) || (ident == ZERO_INFLATED_GAMMA)) {
      min_value = 0.;
      gamma_dist = new gamma_distribution<double>*[nb_state];

      for (i = 0;i < nb_state;i++) {
        if (((ident == GAMMA) && (observation[i]->shape > 0.)) ||
            ((ident == ZERO_INFLATED_GAMMA) && (observation[i]->zero_probability < 1.))) {
          gamma_dist[i] = new gamma_distribution<double>(observation[i]->shape , observation[i]->scale);

          value = quantile(complement(*gamma_dist[i] , GAMMA_TAIL));
          while (max_value < value) {
            max_value += step;
          }
        }

        else {
          gamma_dist[i] = NULL;
        }
      }

      nb_step = (int)(max_value / step) + 1;
    }

    else if (ident == GAUSSIAN) {
      gaussian_dist = new normal*[nb_state];

      for (i = 0;i < nb_state;i++) {
        gaussian_dist[i] = new normal(observation[i]->location , observation[i]->dispersion);

        value = quantile(*gaussian_dist[i] , GAUSSIAN_TAIL);
        while (min_value > value) {
          min_value -= step;
        }

        value = quantile(complement(*gaussian_dist[i] , GAUSSIAN_TAIL));
        while (max_value < value) {
          max_value += step;
        }
      }

      nb_step = (int)((max_value - min_value) / step) + 1;
    }

    else if (ident == VON_MISES) {
      while (min_value >= step / 2) {
        min_value -= step;
      }

      value = (unit == DEGREE ? 360 : 2 * M_PI) - 3 * step / 2;
      while (max_value <= value) {
        max_value += step;
      }

      nb_step = (int)((max_value - min_value) / step) + 1;
    }

    for (i = 0;i < nb_state;i++) {
      frequency[i] = new double[nb_step];
      cumul[i] = new double[nb_step];
    }

    // calcul des lois d'observation discretisees

    switch (ident) {

    case GAMMA : {
      for (i = 0;i < nb_state;i++) {
        if (observation[i]->shape == 0.) {
          cumul[i][0] = 1.;
          frequency[i][0] = cumul[i][0];

          for (j = 1;j < nb_step;j++) {
            cumul[i][j] = cumul[i][0];
            frequency[i][j] = 0.;
          }
        }

        else {
//          value = step;
//          cumul[i][0] = cdf(*gamma_dist[i] , value);
          value = step / 2;
          cumul[i][0] = cdf(*gamma_dist[i] , value + step / 2);
          frequency[i][0] = cumul[i][0];

          for (j = 1;j < nb_step;j++) {
            value += step;
//            cumul[i][j] = cdf(*gamma_dist[i] , value);
            cumul[i][j] = cdf(*gamma_dist[i] , value + step / 2);
            frequency[i][j] = (cumul[i][j] - cumul[i][j - 1]);
          }
        }
      }
      break;
    }

    case ZERO_INFLATED_GAMMA : {
      for (i = 0;i < nb_state;i++) {
        cumul[i][0] = observation[i]->zero_probability;
        frequency[i][0] = cumul[i][0];

        if (observation[i]->zero_probability == 1.) {
          for (j = 1;j < nb_step;j++) {
            cumul[i][j] = cumul[i][0];
            frequency[i][j] = 0.;
          }
        }

        else {
          value = 0.;
          for (j = 1;j < nb_step;j++) {
            value += step;
            cumul[i][j] = observation[i]->zero_probability + (1. - observation[i]->zero_probability) *
                          cdf(*gamma_dist[i] , value);
            frequency[i][j] = (cumul[i][j] - cumul[i][j - 1]);
          }
        }
      }
      break;
    }

    case GAUSSIAN : {
      for (i = 0;i < nb_state;i++) {
        value = min_value;
        cumul[i][0] = cdf(*gaussian_dist[i] , value + step / 2);
        frequency[i][0] = (cumul[i][0] - cdf(*gaussian_dist[i] , value - step / 2));

        for (j = 1;j < nb_step;j++) {
          value += step;
          cumul[i][j] = cdf(*gaussian_dist[i] , value + step / 2);
          frequency[i][j] = (cumul[i][j] - cumul[i][j - 1]);
        }
      }
      break;
    }

    case VON_MISES : {
      for (i = 0;i < nb_state;i++) {
        value = min_value;
        frequency[i][0] = observation[i]->von_mises_mass_computation(value - step / 2 , value + step / 2);
        cumul[i][0] = frequency[i][0];

        for (j = 1;j < nb_step;j++) {
          value += step;
          mass = observation[i]->von_mises_mass_computation(value - step / 2 , value + step / 2);
          frequency[i][j] = mass;
          cumul[i][j] = cumul[i][j - 1] + mass;
        }
      }
      break;
    }
    }

    width[0] = column_width(min_value , max_value);
  }

  for (i = 0;i < nb_state;i++) {
    os << "\n";
    switch (model) {
    case MIXTURE :
      os << STAT_word[STATW_COMPONENT];
      break;
    case HIDDEN_MARKOV :
      os << STAT_word[STATW_STATE];
      break;
    }
    os << " " << i << " " << STAT_word[STATW_OBSERVATION_DISTRIBUTION] << endl;
    observation[i]->ascii_parameter_print(os);
    observation[i]->ascii_characteristic_print(os , file_flag);

    if (((observation_histogram) || (observation_distribution)) &&
        ((ident != ZERO_INFLATED_GAMMA) || ((ident == ZERO_INFLATED_GAMMA) && (observation[i]->zero_probability < 1.)))) {
      if (observation_distribution) {
        os << "\n";
        if (file_flag) {
          os << "# ";
        }
        switch (model) {
        case MIXTURE :
          os << STAT_label[STATL_COMPONENT];
          break;
        case HIDDEN_MARKOV :
          os << STAT_label[STATL_STATE];
          break;
        }
        os << " " << i << " " << STAT_label[STATL_OBSERVATION] << " "
           << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";

        if ((ident == GAMMA) || (ident == ZERO_INFLATED_GAMMA)) {
          observation_distribution[i]->ascii_characteristic_print(os , true , file_flag);
        }
        else if (ident == GAUSSIAN) {
          observation_distribution[i]->ascii_characteristic_print(os , false , file_flag);
        }
        else if (ident == VON_MISES) {
          observation_distribution[i]->ascii_circular_characteristic_print(os , file_flag);
        }
      }

      if (exhaustive) {
        if (i == 0) {
          old_adjust = os.setf(ios::right , ios::adjustfield);
        }

        if (observation_histogram) {
          nb_element = observation_histogram[i]->nb_element;
        }
        else {
          nb_element = observation_distribution[i]->nb_element;
        }

        if (nb_element > 0) {
          if (observation_histogram) {
            observation_cumul = observation_histogram[i]->cumul_computation();
          }
          else {
            observation_cumul = observation_distribution[i]->cumul_computation();
          }

          // calcul des largeurs des colonnes

          if (observation_histogram) {
            width[1] = column_width(observation_histogram[i]->max) + ASCII_SPACE;
          }
          else {
            width[1] = column_width(observation_distribution[i]->max) + ASCII_SPACE;
          }

          width[2] = column_width(nb_step , frequency[i] , nb_element) + ASCII_SPACE;

          if (observation_histogram) {
            width[3] = column_width(observation_histogram[i]->nb_category , observation_cumul) + ASCII_SPACE;
          }
          else {
            width[3] = column_width(observation_distribution[i]->nb_value - observation_distribution[i]->offset ,
                                    observation_cumul + observation_distribution[i]->offset) + ASCII_SPACE;
          }

          width[4] = column_width(nb_step , cumul[i]) + ASCII_SPACE;

          os << "\n";
          if (file_flag) {
            os << "# ";
          }
          os << "   | ";
          switch (model) {
          case MIXTURE :
            os << STAT_label[STATL_COMPONENT];
            break;
          case HIDDEN_MARKOV :
            os << STAT_label[STATL_STATE];
            break;
          }
          os << " " << i << " " << STAT_label[STATL_OBSERVATION] << " ";
          if (observation_histogram) {
            os << STAT_label[STATL_HISTOGRAM];
          }
          else {
            os << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
          }
          os << " | ";
          switch (model) {
          case MIXTURE :
            os << STAT_label[STATL_COMPONENT];
            break;
          case HIDDEN_MARKOV :
            os << STAT_label[STATL_STATE];
            break;
          }
          os << " " << i << " " << STAT_label[STATL_OBSERVATION] << " " << STAT_label[STATL_DISTRIBUTION]
             << " | " << STAT_label[STATL_CUMULATIVE] << " ";
          if (observation_histogram) {
            os << STAT_label[STATL_HISTOGRAM];
          }
          else {
            os << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
          }
          os << " " << STAT_label[STATL_FUNCTION]
             << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_DISTRIBUTION]
             << " " << STAT_label[STATL_FUNCTION] << endl;

          value = min_value;
          if (observation_histogram) {
            j = -1;
          }
          else {
            j = observation_distribution[i]->offset - step;
          }

          for (k = 0;k < nb_step;k++) {
            os << setw(width[0]) << value;

            if ((observation_histogram) && (value >= observation_histogram[i]->min_value) &&
                (value < observation_histogram[i]->max_value)) {
              os << setw(width[1]) << observation_histogram[i]->frequency[++j];
            }
            else if ((observation_distribution) && (value >= observation_distribution[i]->offset) &&
                     (value < observation_distribution[i]->nb_value)) {
              j += step;
              os << setw(width[1]) << observation_distribution[i]->frequency[j];
            }
            else {
              os << setw(width[1]) << " ";
            }

            os << setw(width[2]) << frequency[i][k] * nb_element;

            if (((observation_histogram) && (value >= observation_histogram[i]->min_value) &&
                 (value < observation_histogram[i]->max_value)) ||
                ((observation_distribution) && (value >= observation_distribution[i]->offset) &&
                 (value < observation_distribution[i]->nb_value))) {
              os << setw(width[3]) << observation_cumul[j];
            }
            else {
              os << setw(width[3]) << " ";
            }

            os << setw(width[4]) << cumul[i][k] << endl;

            value += step;
          }

          delete [] observation_cumul;
        }
      }
    }
  }

  if (((marginal_histogram) || (marginal_distribution)) &&
      ((weight) || (restoration_weight))) {
    int nb_negative_step , offset;
    double mean , variance , likelihood , information , *marginal_cumul;
    Distribution *mixture;
    FrequencyDistribution *clustered_histo;
    Test test(CHI2);


    if (marginal_histogram) {
      nb_element = marginal_histogram->nb_element;
      marginal_cumul = marginal_histogram->cumul_computation();
    }

    else {
      nb_element = marginal_distribution->nb_element;
      marginal_cumul = marginal_distribution->cumul_computation();

      clustered_histo = new FrequencyDistribution(*marginal_distribution , 'c' , (int)step);

      nb_negative_step = 0;
      if (ident == GAUSSIAN) {
        value = min_value;
        while (value < 0.) {
          value += step;
          nb_negative_step++;
        }
      }

      offset = 0;
      value = min_value;
      while (value >= 0.) {
        value -= step;
        offset++;
      }
      if (offset > 0) {
        offset--;
      }

      mixture = new Distribution(nb_step - nb_negative_step + offset);
      mixture->offset = offset;

      information = marginal_distribution->information_computation();
    }

    if (exhaustive) {

     // calcul des largeurs des colonnes

      if (marginal_histogram) {
        width[1] = column_width(marginal_histogram->max) + ASCII_SPACE;
        width[3] = column_width(marginal_histogram->nb_category , marginal_cumul) + ASCII_SPACE;
      }
      else {
        width[1] = column_width(marginal_distribution->max) + ASCII_SPACE;
        width[3] = column_width(marginal_distribution->nb_value - marginal_distribution->offset ,
                                marginal_cumul + marginal_distribution->offset) + ASCII_SPACE;
      }
    }

    frequency[nb_state] = new double[nb_step];
    cumul[nb_state] = new double[nb_step];

    if (weight) {
      if ((marginal_distribution) || (exhaustive)) {
        os << "\n";
        if (file_flag) {
          os << "# ";
        }
        os << STAT_label[STATL_MIXTURE] << " - "
           << STAT_label[STATL_THEORETICAL] << " " << STAT_label[STATL_WEIGHTS] << ":";
      }

      // calcul du melange

      for (i = 0;i < nb_step;i++) {
        frequency[nb_state][i] = 0.;
        cumul[nb_state][i] = 0.;
        for (j = 0;j < nb_state;j++) {
          frequency[nb_state][i] += weight->mass[j] * frequency[j][i];
          cumul[nb_state][i] += weight->mass[j] * cumul[j][i];
        }
//        frequency[nb_state][i] *= nb_element;
      }

      if (marginal_distribution) {
        for (i = nb_negative_step;i < nb_step;i++) {
          mixture->mass[i - nb_negative_step + offset] = frequency[nb_state][i];
        }

#       ifdef DEBUG
        os << "\nTest: " << nb_negative_step << " | " << offset << endl;
        for (i = 0;i < nb_step - nb_negative_step + offset;i++) {
          if (i < clustered_histo->nb_value) {
            os << clustered_histo->frequency[i];
          }
          else {
            os << "  ";
          }
          os << " | " << mixture->mass[i] * nb_element << endl;
        }
#       endif

        for (i = 0;i < nb_state;i++) {
          os << " " << weight->mass[i];
        }
        os << endl;

        if (ident != VON_MISES) {
          mean = mean_computation(weight);
          variance = variance_computation(weight , mean);

          if (file_flag) {
            os << "# ";
          }
          os << STAT_label[STATL_MEAN] << ": " << mean << "   "
             << STAT_label[STATL_VARIANCE] << ": " << variance << "   "
             << STAT_label[STATL_STANDARD_DEVIATION] << ": " << sqrt(variance) << endl;

          os << "\n";
          if (file_flag) {
            os << "# ";
          }
          os << STAT_label[STATL_MARGINAL] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";
          marginal_distribution->ascii_characteristic_print(os , false , file_flag);
        }

        likelihood = mixture->likelihood_computation(*clustered_histo);

        os << "\n";
        if (file_flag) {
          os << "# ";
        }
        os << STAT_label[STATL_LIKELIHOOD] << ": " << likelihood << "   ("
           << STAT_label[STATL_NORMALIZED] << ": " << likelihood / marginal_distribution->nb_element << ")" << endl;

        if (file_flag) {
          os << "# ";
        }
        os << STAT_label[STATL_MAX_LIKELIHOOD] << ": " << information << "   ("
           << STAT_label[STATL_INFORMATION] << ": " << information / marginal_distribution->nb_element << ")" << endl;

        if (file_flag) {
          os << "# ";
        }
        os << STAT_label[STATL_DEVIANCE] << ": " << 2 * (information - likelihood) << endl;

        mixture->chi2_fit(*clustered_histo , test);
        os << "\n";
        test.ascii_print(os , file_flag);
      }

      if (exhaustive) {

       // calcul des largeurs des colonnes

        width[2] = column_width(nb_step , frequency[nb_state] , nb_element);
        for (i = 0;i < nb_state;i++) {
          buff = column_width(nb_step , frequency[i] , weight->mass[i] * nb_element);
          if (buff > width[2]) {
            width[2] = buff;
          }
        }
        width[2] += ASCII_SPACE;

        width[4] = column_width(nb_step , cumul[nb_state]) + ASCII_SPACE;

        os << "\n";
        if (file_flag) {
          os << "# ";
        }
        os << "   | " << STAT_label[STATL_MARGINAL] << " ";
        if (marginal_histogram) {
          os << STAT_label[STATL_HISTOGRAM];
        }
        else {
          os << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
        }
        for (i = 0;i < nb_state;i++) {
          os << " | ";
          switch (model) {
          case MIXTURE :
            os << STAT_label[STATL_COMPONENT];
            break;
          case HIDDEN_MARKOV :
            os << STAT_label[STATL_STATE];
            break;
          }
          os << " " << i << " " << STAT_label[STATL_OBSERVATION] << " " << STAT_label[STATL_DISTRIBUTION];
        }
        os << " | " << STAT_label[STATL_MIXTURE] << " | " << STAT_label[STATL_CUMULATIVE] << " ";
        if (marginal_histogram) {
          os << STAT_label[STATL_HISTOGRAM];
        }
        else {
          os << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
        }
        os << " " << STAT_label[STATL_FUNCTION] << " | " << STAT_label[STATL_CUMULATIVE]
           << " " << STAT_label[STATL_MIXTURE] << " " << STAT_label[STATL_FUNCTION] << endl;

        value = min_value;
        if (marginal_histogram) {
          i = -1;
        }
        else {
          i = marginal_distribution->offset - step;
        }

        for (j = 0;j < nb_step;j++) {
          if (file_flag) {
            os << "# ";
          }
          os << setw(width[0]) << value;

          if (marginal_histogram) {
            if ((value >= marginal_histogram->min_value) &&
                (value < marginal_histogram->max_value)) {
              os << setw(width[1]) << marginal_histogram->frequency[++i];
            }
            else {
              os << setw(width[1]) << " ";
            }
          }

          else {
            if ((value >= marginal_distribution->offset) &&
                (value < marginal_distribution->nb_value)) {
              i += step;
              os << setw(width[1]) << marginal_distribution->frequency[i];
            }
            else {
              os << setw(width[1]) << " ";
            }
          }

          for (k = 0;k < nb_state;k++) {
            os << setw(width[2]) << weight->mass[k] * frequency[k][j] * nb_element;
          }
          os << setw(width[2]) << frequency[nb_state][j] * nb_element;

          if (((marginal_histogram) && (value >= marginal_histogram->min_value) &&
               (value < marginal_histogram->max_value)) ||
              ((marginal_distribution) && (value >= marginal_distribution->offset) &&
               (value < marginal_distribution->nb_value))) {
            os << setw(width[3]) << marginal_cumul[i];
          }
          else {
            os << setw(width[3]) << " ";
          }

          os << setw(width[4]) << cumul[nb_state][j] << endl;

          value += step;
        }
      }
    }

    if (restoration_weight) {
      if ((marginal_distribution) || (exhaustive)) {
        os << "\n";
        if (file_flag) {
          os << "# ";
        }
        os << STAT_label[STATL_MIXTURE] << " - "
           << STAT_label[STATL_RESTORATION] << " " << STAT_label[STATL_WEIGHTS] << ":";
      }

      // calcul du melange

      for (i = 0;i < nb_step;i++) {
        frequency[nb_state][i] = 0.;
        cumul[nb_state][i] = 0.;
        for (j = 0;j < nb_state;j++) {
          frequency[nb_state][i] += restoration_weight->mass[j] * frequency[j][i];
          cumul[nb_state][i] += restoration_weight->mass[j] * cumul[j][i];
        }
//        frequency[nb_state][i] *= nb_element;
      }

      if (marginal_distribution) {
        for (i = nb_negative_step;i < nb_step;i++) {
          mixture->mass[i - nb_negative_step + offset] = frequency[nb_state][i];
        }

        for (i = 0;i < nb_state;i++) {
          os << " " << restoration_weight->mass[i];
        }
        os << endl;

        if (ident != VON_MISES) {
          mean = mean_computation(restoration_weight);
          variance = variance_computation(restoration_weight , mean);

          if (file_flag) {
            os << "# ";
          }
          os << STAT_label[STATL_MEAN] << ": " << mean << "   "
             << STAT_label[STATL_VARIANCE] << ": " << variance << "   "
             << STAT_label[STATL_STANDARD_DEVIATION] << ": " << sqrt(variance) << endl;

          os << "\n";
          if (file_flag) {
            os << "# ";
          }
          os << STAT_label[STATL_MARGINAL] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";
          marginal_distribution->ascii_characteristic_print(os , false , file_flag);
        }

        likelihood = mixture->likelihood_computation(*clustered_histo);

        os << "\n";
        if (file_flag) {
          os << "# ";
        }
        os << STAT_label[STATL_LIKELIHOOD] << ": " << likelihood << "   ("
        << STAT_label[STATL_NORMALIZED] << ": " << likelihood / marginal_distribution->nb_element << ")" << endl;

        if (file_flag) {
          os << "# ";
        }
        os << STAT_label[STATL_MAX_LIKELIHOOD] << ": " << information << "   ("
           << STAT_label[STATL_INFORMATION] << ": " << information / marginal_distribution->nb_element << ")" << endl;

        if (file_flag) {
          os << "# ";
        }
        os << STAT_label[STATL_DEVIANCE] << ": " << 2 * (information - likelihood) << endl;

        mixture->chi2_fit(*clustered_histo , test);
        os << "\n";
        test.ascii_print(os , file_flag);
      }

      if (exhaustive) {

       // calcul des largeurs des colonnes

        width[2] = column_width(nb_step , frequency[nb_state] , nb_element);
        for (i = 0;i < nb_state;i++) {
          buff = column_width(nb_step , frequency[i] , restoration_weight->mass[i] * nb_element);
          if (buff > width[2]) {
            width[2] = buff;
          }
        }
        width[2] += ASCII_SPACE;

        width[4] = column_width(nb_step , cumul[nb_state]) + ASCII_SPACE;

        os << "\n";
        if (file_flag) {
          os << "# ";
        }
        os << "   | " << STAT_label[STATL_MARGINAL] << " ";
        if (marginal_histogram) {
          os << STAT_label[STATL_HISTOGRAM];
        }
        else {
          os << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
        }
        for (i = 0;i < nb_state;i++) {
          os << " | ";
          switch (model) {
          case MIXTURE :
            os << STAT_label[STATL_COMPONENT];
            break;
          case HIDDEN_MARKOV :
            os << STAT_label[STATL_STATE];
            break;
          }
          os << " " << i << " " << STAT_label[STATL_OBSERVATION] << " " << STAT_label[STATL_DISTRIBUTION];
        }
        os << " | " << STAT_label[STATL_MIXTURE] << " | " << STAT_label[STATL_CUMULATIVE] << " ";
        if (marginal_histogram) {
          os << STAT_label[STATL_HISTOGRAM];
        }
        else {
          os << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
        }
        os << " " << STAT_label[STATL_FUNCTION] << " | " << STAT_label[STATL_CUMULATIVE]
           << " " << STAT_label[STATL_MIXTURE] << " " << STAT_label[STATL_FUNCTION] << endl;

        value = min_value;
        if (marginal_histogram) {
          i = -1;
        }
        else {
          i = marginal_distribution->offset - step;
        }

        for (j = 0;j < nb_step;j++) {
          if (file_flag) {
            os << "# ";
          }
          os << setw(width[0]) << value;

          if (marginal_histogram) {
            if ((value >= marginal_histogram->min_value) &&
                (value < marginal_histogram->max_value)) {
              os << setw(width[1]) << marginal_histogram->frequency[++i];
            }
            else {
              os << setw(width[1]) << " ";
            }
          }

          else {
            if ((value >= marginal_distribution->offset) &&
                (value < marginal_distribution->nb_value)) {
              i += step;
              os << setw(width[1]) << marginal_distribution->frequency[i];
            }
            else {
              os << setw(width[1]) << " ";
            }
          }

          for (k = 0;k < nb_state;k++) {
            os << setw(width[2]) << restoration_weight->mass[k] * frequency[k][j] * nb_element;
          }
          os << setw(width[2]) << frequency[nb_state][j] * nb_element;

          if (((marginal_histogram) && (value >= marginal_histogram->min_value) &&
               (value < marginal_histogram->max_value)) ||
              ((marginal_distribution) && (value >= marginal_distribution->offset) &&
               (value < marginal_distribution->nb_value))) {
            os << setw(width[3]) << marginal_cumul[i];
          }
          else {
            os << setw(width[3]) << " ";
          }

          os << setw(width[4]) << cumul[nb_state][j] << endl;

          value += step;
        }
      }
    }

    delete [] marginal_cumul;

    if (marginal_distribution) {
      delete clustered_histo;
      delete mixture;
    }

    delete [] frequency[nb_state];
    delete [] cumul[nb_state];
  }

  if ((marginal_histogram) || (marginal_distribution)) {
    if ((ident == GAMMA) || (ident == ZERO_INFLATED_GAMMA)) {
      for (i = 0;i < nb_state;i++) {
        delete gamma_dist[i];
      }
      delete [] gamma_dist;
    }

    else if (ident ==  GAUSSIAN) {
      for (i = 0;i < nb_state;i++) {
        delete gaussian_dist[i];
      }
      delete [] gaussian_dist;
    }

    for (i = 0;i < nb_state;i++) {
      delete [] frequency[i];
    }
    delete [] frequency;

    for (i = 0;i < nb_state;i++) {
      delete [] cumul[i];
    }
    delete [] cumul;

    os.setf((FMTFLAGS)old_adjust , ios::adjustfield);
  }

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet ContinuousParametricProcess au format tableur.
 *
 *  arguments : stream, pointeurs sur les histogrammes d'observation ou
 *              les lois d'observation empiriques et sur l'histogramme marginale ou
 *              la loi marginale empirique, type de modele.
 *
 *--------------------------------------------------------------*/

ostream& ContinuousParametricProcess::spreadsheet_print(ostream &os , Histogram **observation_histogram ,
                                                        FrequencyDistribution **observation_distribution ,
                                                        Histogram *marginal_histogram ,
                                                        FrequencyDistribution *marginal_distribution ,
                                                        int model) const

{
  register int i , j , k;
  int nb_step , nb_element;
  double step , value , min_value , max_value , mass , *observation_cumul ,
         **frequency , **cumul;
  gamma_distribution<double> **gamma_dist;
  normal **gaussian_dist;


  if ((marginal_histogram) || (marginal_distribution)) {
    if (marginal_histogram) {
      step = marginal_histogram->step;
    }
    else {
      step = marginal_distribution->min_interval_computation();
    }

    if ((ident == GAMMA) || (ident == ZERO_INFLATED_GAMMA)) {
      if (marginal_histogram) {
        max_value = marginal_histogram->max_value;

/*        switch (marginal_histogram->type) {
        case INT_VALUE :
          max_value = marginal_histogram->max_value;
          break;
        case REAL_VALUE :
          max_value = marginal_histogram->max_value + marginal_histogram->step;
          break;
        } */
      }

      else {
        max_value = marginal_distribution->nb_value - 1;
      }
    }

    else if ((ident == GAUSSIAN) || (ident == VON_MISES)) {
      if (marginal_histogram) {
        min_value = marginal_histogram->min_value;
        max_value = marginal_histogram->max_value;

/*        switch (marginal_histogram->type) {

        case INT_VALUE :
          min_value = marginal_histogram->min_value;
          max_value = marginal_histogram->max_value;
          break;
        case REAL_VALUE :
          min_value = marginal_histogram->min_value - marginal_histogram->step / 2;
          max_value = marginal_histogram->max_value + marginal_histogram->step / 2;
          break;
        } */
      }

      else {
        min_value = marginal_distribution->offset;
        max_value = marginal_distribution->nb_value - 1;
      }
    }

    frequency = new double*[nb_state + 1];
    cumul = new double*[nb_state + 1];

    if ((ident == GAMMA) || (ident == ZERO_INFLATED_GAMMA)) {
      min_value = 0.;
      gamma_dist = new gamma_distribution<double>*[nb_state];

      for (i = 0;i < nb_state;i++) {
        if (((ident == GAMMA) && (observation[i]->shape > 0.)) ||
            ((ident == ZERO_INFLATED_GAMMA) && (observation[i]->zero_probability < 1.))) {
          gamma_dist[i] = new gamma_distribution<double>(observation[i]->shape , observation[i]->scale);

          value = quantile(complement(*gamma_dist[i] , GAMMA_TAIL));
          while (max_value < value) {
            max_value += step;
          }
        }

        else {
          gamma_dist[i] = NULL;
        }
      }

      nb_step = (int)(max_value / step) + 1;
    }

    else if (ident == GAUSSIAN) {
      gaussian_dist = new normal*[nb_state];

      for (i = 0;i < nb_state;i++) {
        gaussian_dist[i] = new normal(observation[i]->location , observation[i]->dispersion);

        value = quantile(*gaussian_dist[i] , GAUSSIAN_TAIL);
        while (min_value > value) {
          min_value -= step;
        }

        value = quantile(complement(*gaussian_dist[i] , GAUSSIAN_TAIL));
        while (max_value < value) {
          max_value += step;
        }
      }

      nb_step = (int)((max_value - min_value) / step) + 1;
    }

    else if (ident == VON_MISES) {
      while (min_value >= step / 2) {
        min_value -= step;
      }

      value = (unit == DEGREE ? 360 : 2 * M_PI) - 3 * step / 2;
      while (max_value <= value) {
        max_value += step;
      }

      nb_step = (int)((max_value - min_value) / step) + 1;
    }

    for (i = 0;i < nb_state;i++) {
      frequency[i] = new double[nb_step];
    }
    for (i = 0;i < nb_state;i++) {
      cumul[i] = new double[nb_step];
    }

    // calcul des lois d'observation discretisees

    switch (ident) {

    case GAMMA : {
      for (i = 0;i < nb_state;i++) {
        if (observation[i]->shape == 0.) {
          cumul[i][0] = 1.;
          frequency[i][0] = cumul[i][0];

          for (j = 1;j < nb_step;j++) {
            cumul[i][j] = cumul[i][0];
            frequency[i][j] = 0.;
          }
        }

        else {
//          value = step;
//          cumul[i][0] = cdf(*gamma_dist[i] , value);
          value = step / 2;
          cumul[i][0] = cdf(*gamma_dist[i] , value + step / 2);
          frequency[i][0] = cumul[i][0];

          for (j = 1;j < nb_step;j++) {
            value += step;
//            cumul[i][j] = cdf(*gamma_dist[i] , value);
            cumul[i][j] = cdf(*gamma_dist[i] , value + step / 2);
            frequency[i][j] = (cumul[i][j] - cumul[i][j - 1]);
          }
        }
      }
      break;
    }

    case ZERO_INFLATED_GAMMA : {
      for (i = 0;i < nb_state;i++) {
        cumul[i][0] = observation[i]->zero_probability;
        frequency[i][0] = cumul[i][0];

        if (observation[i]->zero_probability == 1.) {
          for (j = 1;j < nb_step;j++) {
            cumul[i][j] =  cumul[i][0];
            frequency[i][j] = 0.;
          }
        }

        else {
          value = 0.;
          for (j = 1;j < nb_step;j++) {
            value += step;
            cumul[i][j] = observation[i]->zero_probability + (1. - observation[i]->zero_probability) *
                          cdf(*gamma_dist[i] , value);
            frequency[i][j] = (cumul[i][j] - cumul[i][j - 1]);
          }
        }
      }
      break;
    }

    case GAUSSIAN : {
      for (i = 0;i < nb_state;i++) {
        value = min_value;
        cumul[i][0] = cdf(*gaussian_dist[i] , value + step / 2);
        frequency[i][0] = (cumul[i][0] - cdf(*gaussian_dist[i] , value - step / 2));

        for (j = 1;j < nb_step;j++) {
          value += step;
          cumul[i][j] = cdf(*gaussian_dist[i] , value + step / 2);
          frequency[i][j] = (cumul[i][j] - cumul[i][j - 1]);
        }
      }
      break;
    }

    case VON_MISES : {
      for (i = 0;i < nb_state;i++) {
        value = min_value;
        frequency[i][0] = observation[i]->von_mises_mass_computation(value - step / 2 , value + step / 2);
        cumul[i][0] = frequency[i][0];

        for (j = 1;j < nb_step;j++) {
          value += step;
          mass = observation[i]->von_mises_mass_computation(value - step / 2 , value + step / 2);
          frequency[i][j] = mass;
          cumul[i][j] = cumul[i][j - 1] + mass;
        }
      }
      break;
    }
    }
  }

  for (i = 0;i < nb_state;i++) {
    os << "\n";
    switch (model) {
    case MIXTURE :
      os << STAT_word[STATW_COMPONENT];
      break;
    case HIDDEN_MARKOV :
      os << STAT_word[STATW_STATE];
      break;
    }
    os << " " << i << "\t" << STAT_word[STATW_OBSERVATION_DISTRIBUTION] << endl;
    observation[i]->spreadsheet_parameter_print(os);
    observation[i]->spreadsheet_characteristic_print(os);

    if (((observation_histogram) || (observation_distribution)) &&
        ((ident != ZERO_INFLATED_GAMMA) || ((ident == ZERO_INFLATED_GAMMA) && (observation[i]->zero_probability < 1.)))) {
      if (observation_distribution) {
        os << "\n";
        switch (model) {
        case MIXTURE :
          os << STAT_label[STATL_COMPONENT];
          break;
        case HIDDEN_MARKOV :
          os << STAT_label[STATL_STATE];
          break;
        }
        os << " " << i << " " << STAT_label[STATL_OBSERVATION]
           << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\t";

        if ((ident == GAMMA) || (ident == ZERO_INFLATED_GAMMA)) {
          observation_distribution[i]->spreadsheet_characteristic_print(os , true);
        }
        else if (ident == GAUSSIAN) {
          observation_distribution[i]->spreadsheet_characteristic_print(os , false);
        }
        else if (ident == VON_MISES) {
          observation_distribution[i]->spreadsheet_circular_characteristic_print(os);
        }
      }

      if (observation_histogram) {
        nb_element = observation_histogram[i]->nb_element;
      }
      else {
        nb_element = observation_distribution[i]->nb_element;
      }

      if (nb_element > 0) {
        if (observation_histogram) {
          observation_cumul = observation_histogram[i]->cumul_computation();
        }
        else {
          observation_cumul = observation_distribution[i]->cumul_computation();
        }

        os << "\n\t";
        switch (model) {
        case MIXTURE :
          os << STAT_label[STATL_COMPONENT];
          break;
        case HIDDEN_MARKOV :
          os << STAT_label[STATL_STATE];
          break;
        }
        os << " " << i << " " << STAT_label[STATL_OBSERVATION];
        if (observation_histogram) {
          os << " " << STAT_label[STATL_HISTOGRAM];
        }
        else {
          os << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
        }
        os << "\t";
        switch (model) {
        case MIXTURE :
          os << STAT_label[STATL_COMPONENT];
          break;
        case HIDDEN_MARKOV :
          os << STAT_label[STATL_STATE];
          break;
        }
        os << " " << i << " " << STAT_label[STATL_OBSERVATION] << " " << STAT_label[STATL_DISTRIBUTION]
           << "\t" << STAT_label[STATL_CUMULATIVE] << " ";
        if (observation_histogram) {
          os << STAT_label[STATL_HISTOGRAM];
        }
        else {
          os << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
        }
        os  << " " << STAT_label[STATL_FUNCTION]
           << "\t" << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_DISTRIBUTION]
           << " " << STAT_label[STATL_FUNCTION] << endl;

        value = min_value;
        if (observation_histogram) {
          j = -1;
        }
        else {
          j = observation_distribution[i]->offset - step;
        }

        for (k = 0;k < nb_step;k++) {
          os << value << "\t";

          if ((observation_histogram) && (value >= observation_histogram[i]->min_value) &&
              (value < observation_histogram[i]->max_value)) {
            os << observation_histogram[i]->frequency[++j];
          }
          if ((observation_distribution) && (value >= observation_distribution[i]->offset) &&
              (value < observation_distribution[i]->nb_value)) {
            j += step;
            os << observation_distribution[i]->frequency[j];
          }

          os << "\t" << frequency[i][k] * nb_element;

          os << "\t";
          if (((observation_histogram) && (value >= observation_histogram[i]->min_value) &&
               (value < observation_histogram[i]->max_value)) ||
              ((observation_distribution) && (value >= observation_distribution[i]->offset) &&
               (value < observation_distribution[i]->nb_value))) {
            os << observation_cumul[j];
          }

          os << "\t" << cumul[i][k] << endl;

          value += step;
        }

        delete [] observation_cumul;
      }
    }
  }

  if ((!marginal_histogram) && (!marginal_distribution)) {
    os << "\n";
    for (i = 0;i < nb_state;i++) {
      os << "\t";
      switch (model) {
      case MIXTURE :
        os << STAT_label[STATL_COMPONENT];
        break;
      case HIDDEN_MARKOV :
        os << STAT_label[STATL_STATE];
        break;
      }
      os << " " << i << " " << STAT_label[STATL_OBSERVATION] << " " << STAT_label[STATL_DISTRIBUTION];
    }
    os << endl;

    switch (ident) {

    case GAMMA : {
      gamma_dist = new gamma_distribution<double>*[nb_state];
      max_value = 0.;

      for (i = 0;i < nb_state;i++) {
        if (observation[i]->shape > 0.) {
          gamma_dist[i] = new gamma_distribution<double>(observation[i]->shape , observation[i]->scale);

          value = quantile(complement(*gamma_dist[i] , GAMMA_TAIL));
          if (value > max_value) {
            max_value = value;
          }
        }

        else {
          gamma_dist[i] = NULL;
        }
      }

      step = max_value / GAMMA_NB_STEP + 1;
      value = 0.;

      for (i = 0;i < nb_step;i++) {
        os << value;
        for (j = 0;j < nb_state;j++) {
          os << "\t";
          if (observation[i]->shape > 0.) {
            os << pdf(*gamma_dist[j] , value);
          }
        }
        value += step;
        os << endl;
      }

      for (i = 0;i < nb_state;i++) {
        delete gamma_dist[i];
      }
      delete [] gamma_dist;
      break;
    }

    case ZERO_INFLATED_GAMMA : {
      gamma_dist = new gamma_distribution<double>*[nb_state];
      max_value = 0.;

      for (i = 0;i < nb_state;i++) {
        if (observation[i]->zero_probability < 1.) {
          gamma_dist[i] = new gamma_distribution<double>(observation[i]->shape , observation[i]->scale);

          value = quantile(complement(*gamma_dist[i] , GAMMA_TAIL));
          if (value > max_value) {
            max_value = value;
          }
        }

        else {
          gamma_dist[i] = NULL;
        }
      }

      step = max_value / GAMMA_NB_STEP + 1;

      value = 0.;
      os << value;
      for (i = 0;i < nb_state;i++) {
        os << "\t" << observation[i]->zero_probability;
      }
      os << endl;

      for (i = 1;i < nb_step;i++) {
        value += step;
        os << value;
        for (j = 0;j < nb_state;j++) {
          os << "\t";
          if (observation[i]->zero_probability < 1.) {
            os << (1. - observation[j]->zero_probability) * pdf(*gamma_dist[j] , value);
          }
        }
        os << endl;
      }

      for (i = 0;i < nb_state;i++) {
        delete gamma_dist[i];
      }
      delete [] gamma_dist;
      break;
    }

    case GAUSSIAN : {
      gaussian_dist = new normal*[nb_state];

      for (i = 0;i < nb_state;i++) {
        gaussian_dist[i] = new normal(observation[i]->location , observation[i]->dispersion);
      }

      min_value = quantile(*gaussian_dist[0] , GAUSSIAN_TAIL);
      for (i = 1;i < nb_state;i++) {
        value = quantile(*gaussian_dist[i] , GAUSSIAN_TAIL);
        if (value < min_value) {
          min_value = value;
        }
      }

      max_value = quantile(complement(*gaussian_dist[0] , GAUSSIAN_TAIL));
      for (i = 1;i < nb_state;i++) {
        value = quantile(complement(*gaussian_dist[i] , GAUSSIAN_TAIL));
        if (value > max_value) {
          max_value = value;
        }
      }

      step = (max_value - min_value) / GAUSSIAN_NB_STEP + 1;

      value = min_value;
      for (i = 0;i < GAUSSIAN_NB_STEP;i++) {
        os << value;
        for (j = 0;j < nb_state;j++) {
          os << "\t" << pdf(*gaussian_dist[j] , value);
        }
        value += step;
        os << endl;
      }

      for (i = 0;i < nb_state;i++) {
        delete gaussian_dist[i];
      }
      delete [] gaussian_dist;
      break;
    }

    case VON_MISES : {
      value = 0.;

      switch (unit) {

      case DEGREE : {
        step = 360. / VON_MISES_NB_STEP;

        for (i = 0;i < VON_MISES_NB_STEP;i++) {
          os << value;
          for (j = 0;j < nb_state;j++) {
            os << "\t" << exp(observation[j]->dispersion * cos((value - observation[j]->location) * M_PI / 180)) /
                          (360 * cyl_bessel_i(0 , observation[j]->dispersion));
          }
          value += step;
          os << endl;
        }
        break;
      }

      case RADIAN : {
        step = 2 * M_PI / VON_MISES_NB_STEP;

        for (i = 0;i < VON_MISES_NB_STEP;i++) {
          os << value;
          for (j = 0;j < nb_state;j++) {
            os << "\t" << exp(observation[j]->dispersion * cos(value - observation[j]->location)) /
                          (2 * M_PI * cyl_bessel_i(0 , observation[j]->dispersion));
          }
          value += step;
          os << endl;
        }
        break;
      }
      }
      break;
    }
    }
  }

  if (((marginal_histogram) || (marginal_distribution)) &&
      ((weight) || (restoration_weight))) {
    int nb_negative_step , offset;
    double mean , variance , likelihood , information , *marginal_cumul;
    Distribution *mixture;
    FrequencyDistribution *clustered_histo;
    Test test(CHI2);


    if (marginal_histogram) {
      nb_element = marginal_histogram->nb_element;
      marginal_cumul = marginal_histogram->cumul_computation();
    }

    else {
      nb_element = marginal_distribution->nb_element;
      marginal_cumul = marginal_distribution->cumul_computation();

      clustered_histo = new FrequencyDistribution(*marginal_distribution , 'c' , (int)step);

      nb_negative_step = 0;
      if (ident == GAUSSIAN) {
        value = min_value;
        while (value < 0.) {
          value += step;
          nb_negative_step++;
        }
      }

      offset = 0;
      value = min_value;
      while (value >= 0.) {
        value -= step;
        offset++;
      }
      if (offset > 0) {
        offset--;
      }

      mixture = new Distribution(nb_step - nb_negative_step + offset);
      mixture->offset = offset;

      information = marginal_distribution->information_computation();
    }

    frequency[nb_state] = new double[nb_step];
    cumul[nb_state] = new double[nb_step];

    if (weight) {
      os << "\n" << STAT_label[STATL_MIXTURE] << "\t"
         << STAT_label[STATL_THEORETICAL] << " " << STAT_label[STATL_WEIGHTS];
      for (i = 0;i < nb_state;i++) {
        os << "\t" << weight->mass[i];
      }
      os << endl;

      // calcul du melange

      for (i = 0;i < nb_step;i++) {
        frequency[nb_state][i] = 0.;
        cumul[nb_state][i] = 0.;
        for (j = 0;j < nb_state;j++) {
          frequency[nb_state][i] += weight->mass[j] * frequency[j][i];
          cumul[nb_state][i] += weight->mass[j] * cumul[j][i];
        }
//        frequency[nb_state][i] *= nb_element;
      }

      if (marginal_distribution) {
        for (i = nb_negative_step;i < nb_step;i++) {
          mixture->mass[i - nb_negative_step + offset] = frequency[nb_state][i];
        }

        if (ident != VON_MISES) {
          mean = mean_computation(weight);
          variance = variance_computation(weight , mean);

          os << STAT_label[STATL_MEAN] << "\t" << mean << "\t"
             << STAT_label[STATL_VARIANCE] << "\t" << variance << "\t"
             << STAT_label[STATL_STANDARD_DEVIATION] << "\t" << sqrt(variance) << endl;

          os << "\n" << STAT_label[STATL_MARGINAL] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\t";
          marginal_distribution->spreadsheet_characteristic_print(os);
        }

        likelihood = mixture->likelihood_computation(*clustered_histo);

        os << "\n" << STAT_label[STATL_LIKELIHOOD] << "\t" << likelihood << "\t"
           << STAT_label[STATL_NORMALIZED] << "\t" << likelihood / marginal_distribution->nb_element << endl;
        os << STAT_label[STATL_MAX_LIKELIHOOD] << "\t" << information << "\t"
           << STAT_label[STATL_INFORMATION] << "\t" << information / marginal_distribution->nb_element << endl;
        os << STAT_label[STATL_DEVIANCE] << "\t" << 2 * (information - likelihood) << endl;

        mixture->chi2_fit(*clustered_histo , test);
        os << "\n";
        test.spreadsheet_print(os);
      }

      os << "\n\t" << STAT_label[STATL_MARGINAL] << " ";
      if (marginal_histogram) {
        os << STAT_label[STATL_HISTOGRAM];
      }
      else {
        os << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
      }
      for (i = 0;i < nb_state;i++) {
        os << "\t";
        switch (model) {
        case MIXTURE :
          os << STAT_label[STATL_COMPONENT];
          break;
        case HIDDEN_MARKOV :
          os << STAT_label[STATL_STATE];
          break;
        }
        os << " " << i << " " << STAT_label[STATL_OBSERVATION] << " " << STAT_label[STATL_DISTRIBUTION];
      }
      os << "\t" << STAT_label[STATL_MIXTURE] << "\t" << STAT_label[STATL_CUMULATIVE] << " ";
      if (marginal_histogram) {
        os << STAT_label[STATL_HISTOGRAM];
      }
      else {
        os << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
      }
      os << " " << STAT_label[STATL_FUNCTION] << "\t" << STAT_label[STATL_CUMULATIVE]
         << " " << STAT_label[STATL_MIXTURE] << " " << STAT_label[STATL_FUNCTION] << endl;

      value = min_value;
      if (marginal_histogram) {
        i = -1;
      }
      else {
        i = marginal_distribution->offset - step;
      }

      for (j = 0;j < nb_step;j++) {
        os << value << "\t";

        if ((marginal_histogram) && (value >= marginal_histogram->min_value) &&
            (value < marginal_histogram->max_value)) {
          os << marginal_histogram->frequency[++i];
        }
        if ((marginal_distribution) && (value >= marginal_distribution->offset) &&
            (value < marginal_distribution->nb_value)) {
          i += step;
          os << marginal_distribution->frequency[i];
        }

        for (k = 0;k < nb_state;k++) {
          os << "\t" << weight->mass[k] * frequency[k][j] * nb_element;
        }
        os << "\t" << frequency[nb_state][j] * nb_element;

        os << "\t";
        if (((marginal_histogram) &&(value >= marginal_histogram->min_value) &&
             (value < marginal_histogram->max_value)) ||
            ((marginal_distribution) &&(value >= marginal_distribution->offset) &&
             (value < marginal_distribution->nb_value))) {
          os << marginal_cumul[i];
        }

        os << "\t" << cumul[nb_state][j] << endl;

        value += step;
      }
    }

    if (restoration_weight) {
      os << "\n" << STAT_label[STATL_MIXTURE] << "\t"
         << STAT_label[STATL_RESTORATION] << " " << STAT_label[STATL_WEIGHTS];
      for (i = 0;i < nb_state;i++) {
        os << "\t" << restoration_weight->mass[i];
      }
      os << endl;

      // calcul du melange

      for (i = 0;i < nb_step;i++) {
        frequency[nb_state][i] = 0.;
        cumul[nb_state][i] = 0.;
        for (j = 0;j < nb_state;j++) {
          frequency[nb_state][i] += restoration_weight->mass[j] * frequency[j][i];
          cumul[nb_state][i] += restoration_weight->mass[j] * cumul[j][i];
        }
//        frequency[nb_state][i] *= nb_element;
      }

      if (marginal_distribution) {
        for (i = nb_negative_step;i < nb_step;i++) {
          mixture->mass[i - nb_negative_step + offset] = frequency[nb_state][i];
        }

        if (ident != VON_MISES) {
          mean = mean_computation(restoration_weight);
          variance = variance_computation(restoration_weight , mean);

          os << STAT_label[STATL_MEAN] << "\t" << mean << "\t"
             << STAT_label[STATL_VARIANCE] << "\t" << variance << "\t"
             << STAT_label[STATL_STANDARD_DEVIATION] << "\t" << sqrt(variance) << endl;

          os << "\n" << STAT_label[STATL_MARGINAL] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\t";
          marginal_distribution->spreadsheet_characteristic_print(os);
        }

        likelihood = mixture->likelihood_computation(*clustered_histo);

        os << "\n" << STAT_label[STATL_LIKELIHOOD] << "\t" << likelihood << "\t"
           << STAT_label[STATL_NORMALIZED] << "\t" << likelihood / marginal_distribution->nb_element << endl;
        os << STAT_label[STATL_MAX_LIKELIHOOD] << "\t" << information << "\t"
           << STAT_label[STATL_INFORMATION] << "\t" << information / marginal_distribution->nb_element << endl;
        os << STAT_label[STATL_DEVIANCE] << "\t" << 2 * (information - likelihood) << endl;

        mixture->chi2_fit(*clustered_histo , test);
        os << "\n";
        test.spreadsheet_print(os);
      }

      os << "\n\t" << STAT_label[STATL_MARGINAL] << " ";
      if (marginal_histogram) {
        os << STAT_label[STATL_HISTOGRAM];
      }
      else {
        os << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
      }
      for (i = 0;i < nb_state;i++) {
        os << "\t";
        switch (model) {
        case MIXTURE :
          os << STAT_label[STATL_COMPONENT];
          break;
        case HIDDEN_MARKOV :
          os << STAT_label[STATL_STATE];
          break;
        }
        os << " " << i << " " << STAT_label[STATL_OBSERVATION] << " " << STAT_label[STATL_DISTRIBUTION];
      }
      os << "\t" << STAT_label[STATL_MIXTURE] << "\t" << STAT_label[STATL_CUMULATIVE] << " ";
      if (marginal_histogram) {
        os << STAT_label[STATL_HISTOGRAM];
      }
      else {
        os << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
      }
      os << " " << STAT_label[STATL_FUNCTION] << "\t" << STAT_label[STATL_CUMULATIVE]
         << " " << STAT_label[STATL_MIXTURE] << " " << STAT_label[STATL_FUNCTION] << endl;

      value = min_value;
      if (marginal_histogram) {
        i = -1;
      }
      else {
        i = marginal_distribution->offset - step;
      }

      for (j = 0;j < nb_step;j++) {
        os << value << "\t";

        if ((marginal_histogram) && (value >= marginal_histogram->min_value) &&
            (value < marginal_histogram->max_value)) {
          os << marginal_histogram->frequency[++i];
        }
        if ((marginal_distribution) && (value >= marginal_distribution->offset) &&
            (value < marginal_distribution->nb_value)) {
          i += step;
          os << marginal_distribution->frequency[i];
        }

        for (k = 0;k < nb_state;k++) {
          os << "\t" << restoration_weight->mass[k] * frequency[k][j] * nb_element;
        }
        os << "\t" << frequency[nb_state][j] * nb_element;

        os << "\t";
        if (((marginal_histogram) &&(value >= marginal_histogram->min_value) &&
             (value < marginal_histogram->max_value)) ||
            ((marginal_distribution) &&(value >= marginal_distribution->offset) &&
             (value < marginal_distribution->nb_value))) {
          os << marginal_cumul[i];
        }

        os << "\t" << cumul[nb_state][j] << endl;

        value += step;
      }
    }

    delete [] marginal_cumul;

    if (marginal_distribution) {
      delete clustered_histo;
      delete mixture;
    }

    delete [] frequency[nb_state];
    delete [] cumul[nb_state];
  }

  if ((marginal_histogram) || (marginal_distribution)) {
    if ((ident == GAMMA) || (ident == ZERO_INFLATED_GAMMA)) {
      for (i = 0;i < nb_state;i++) {
        delete gamma_dist[i];
      }
      delete [] gamma_dist;
    }

    else if (ident == GAUSSIAN) {
      for (i = 0;i < nb_state;i++) {
        delete gaussian_dist[i];
      }
      delete [] gaussian_dist;
    }

    for (i = 0;i < nb_state;i++) {
      delete [] frequency[i];
    }
    delete [] frequency;

    for (i = 0;i < nb_state;i++) {
      delete [] cumul[i];
    }
    delete [] cumul;
  }

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Calcul d'un q-q plot.
 *
 *  arguments : valeur minimum, pas, nombre de pas, pointeurs sur
 *              la fonction de repartition theorique, nombre de valeurs,
 *              pointeurs sur la fonction de repartition empirique.
 *
 *--------------------------------------------------------------*/

double** q_q_plot_computation(double min_value , double step ,
                              int nb_step , double *theoretical_cdf ,
                              int nb_value , double **empirical_cdf)

{
  register int i , j;
  double value , **qqplot;


  qqplot = new double*[2];
  qqplot[0] = new double[nb_value];
  qqplot[1] = new double[nb_value];

  for (i = 0;i < nb_value;i++) {
    qqplot[0][i] = empirical_cdf[0][i];
  }

  value = min_value;
  i = 0;
  for (j = 0;j < nb_value;j++) {
    while ((i < nb_step) && (theoretical_cdf[i] < empirical_cdf[1][j])) {
      i++;
      value += step;
    }

    if (i == 0) {
      qqplot[1][j] = value - step * (theoretical_cdf[i] - empirical_cdf[1][j]) /
                     theoretical_cdf[i];
    }
    else {
      qqplot[1][j] = value - step * (theoretical_cdf[i] - empirical_cdf[1][j]) /
                     (theoretical_cdf[i] - theoretical_cdf[i - 1]);
    }
  }

  return qqplot;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un q-q plot au format Gnuplot.
 *
 *  arguments : path, nombre de valeurs, pointeur sur le q-q plot.
 *
 *--------------------------------------------------------------*/

bool q_q_plot_print(const char *path , int nb_value , double **qqplot)

{
  bool status = false;
  register int i;
  ofstream out_file(path);


  if (out_file) {
    status = true;

    for (i = 0;i < nb_value;i++) {
      out_file << qqplot[0][i] << " " << qqplot[1][i] << endl;
    }

    delete [] qqplot[0];
    delete [] qqplot[1];
    delete [] qqplot;
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Sortie Gnuplot d'un objet ContinuousParametricProcess.
 *
 *  arguments : prefixe des fichiers, titre des figures,
 *              indice du processus d'observation,
 *              pointeurs sur les histogrammes d'observation ou
 *              les lois d'observation empiriques et sur l'histogramme marginale ou
 *              la loi marginale empirique, nombre de valeurs,
 *              pointeur sur la fonction de repartition empirique,
 *              type de modele.
 *
 *--------------------------------------------------------------*/

bool ContinuousParametricProcess::plot_print(const char *prefix , const char *title ,
                                             int process , Histogram **observation_histogram ,
                                             FrequencyDistribution **observation_distribution ,
                                             Histogram *marginal_histogram ,
                                             FrequencyDistribution *marginal_distribution ,
                                             int nb_value , double **empirical_cdf , int model) const

{
  bool status = false;
  register int i , j , k;
  int min_interval , nb_step , dist_index;
  double value , min_value , max_value , step , buff , max , *scale , *dist_max ,
         **cumul , **frequency , **qqplot;
  gamma_distribution<double> **gamma_dist;
  normal **gaussian_dist;
  ostringstream data_file_name[NB_STATE + 4];


  data_file_name[0] << prefix << process << 0 << ".dat";
  ofstream out_data_file((data_file_name[0].str()).c_str());

  if (out_data_file) {
    status = true;

    frequency = new double*[nb_state + 2];
    frequency[nb_state] = NULL;
    frequency[nb_state + 1] = NULL;

    cumul = new double*[nb_state + 2];
    cumul[nb_state] = NULL;
    cumul[nb_state + 1] = NULL;

    scale = new double[nb_state + 1];
    dist_max = new double[nb_state];

    // calcul des lois d'observation discretisees

/*    if (empirical_cdf) {
      step = empirical_cdf[0][1] - empirical_cdf[0][0];
      for (i = 2;i < nb_value;i++) {
        if (empirical_cdf[0][i] - empirical_cdf[0][i - 1] < step) {
          step = empirical_cdf[0][i] - empirical_cdf[0][i - 1];
        }
      }
    }
    else { */
      step = D_DEFAULT;
//    }

    if (marginal_distribution) {
      min_interval = marginal_distribution->min_interval_computation();
    }

    switch (ident) {

    case GAMMA : {
      gamma_dist = new gamma_distribution<double>*[nb_state];

      min_value = 0.;
      max_value = 0.;

      for (i = 0;i < nb_state;i++) {
        if (observation[i]->shape > 0.) {
          gamma_dist[i] = new gamma_distribution<double>(observation[i]->shape , observation[i]->scale);

          value = quantile(complement(*gamma_dist[i] , GAMMA_TAIL));
          if (value > max_value) {
            max_value = value;
          }
        }

        else {
          gamma_dist[i] = NULL;
        }
      }

      if ((marginal_histogram) && (marginal_histogram->max_value + marginal_histogram->step > max_value)) {
        max_value = marginal_histogram->max_value + marginal_histogram->step;
      }
      if ((marginal_distribution) && (marginal_distribution->nb_value - 1 > max_value)) {
        max_value = marginal_distribution->nb_value - 1;
      }

      if (step == D_DEFAULT) {
        step = max_value / GAMMA_NB_STEP;
      }
      else {
        buff = max_value / GAMMA_NB_STEP;
        if (buff < step) {
          step = buff;
        }
      }

      if (marginal_histogram) {
        buff = marginal_histogram->step / GAMMA_NB_SUB_STEP;
      }
      else if (marginal_distribution) {
        buff = (double)min_interval / GAMMA_NB_SUB_STEP;
      }
      if (((marginal_histogram) || (marginal_distribution)) && (buff < step)) {
        step = buff;
      }

#     ifdef DEBUG
      cout << "\nTest: " << max_value << " | " << step << endl;
#     endif

      nb_step = (int)(max_value / step) + 1;

      for (i = 0;i < nb_state;i++) {
        frequency[i] = new double[nb_step];
        cumul[i] = new double[nb_step];

        dist_max[i] = 0.;
        value = 0.;

        if (observation[i]->shape == 0.) {
          frequency[i][0] = 1.;
          cumul[i][0] = frequency[i][0];

          for (j = 1;j < nb_step;j++) {
            frequency[i][j] = 0.;
            cumul[i][j] = cumul[i][0];
          }
        }

        else {
          for (j = 0;j < nb_step;j++) {
            frequency[i][j] = pdf(*gamma_dist[i] , value);
            if (frequency[i][j] > dist_max[i]) {
              dist_max[i] = frequency[i][j];
            }
            cumul[i][j] = cdf(*gamma_dist[i] , value);
            value += step;
          }
        }
      }
      break;
    }

    case ZERO_INFLATED_GAMMA : {
      gamma_dist = new gamma_distribution<double>*[nb_state];

      min_value = 0.;
      max_value = 0.;

      for (i = 0;i < nb_state;i++) {
        if (observation[i]->zero_probability < 1.) {
          gamma_dist[i] = new gamma_distribution<double>(observation[i]->shape , observation[i]->scale);

          value = quantile(complement(*gamma_dist[i] , GAMMA_TAIL));
          if (value > max_value) {
            max_value = value;
          }
        }

        else {
          gamma_dist[i] = NULL;
        }
      }

      if ((marginal_histogram) && (marginal_histogram->max_value + marginal_histogram->step > max_value)) {
        max_value = marginal_histogram->max_value + marginal_histogram->step;
      }
      if ((marginal_distribution) && (marginal_distribution->nb_value - 1 > max_value)) {
        max_value = marginal_distribution->nb_value - 1;
      }

      if (step == D_DEFAULT) {
        step = max_value / GAMMA_NB_STEP;
      }
      else {
        buff = max_value / GAMMA_NB_STEP;
        if (buff < step) {
          step = buff;
        }
      }

      if (marginal_histogram) {
        buff = marginal_histogram->step / GAMMA_NB_SUB_STEP;
      }
      else if (marginal_distribution) {
        buff = (double)min_interval / GAMMA_NB_SUB_STEP;
      }
      if (((marginal_histogram) || (marginal_distribution)) && (buff < step)) {
        step = buff;
      }

#     ifdef DEBUG
      cout << "\nTest: " << max_value << " | " << step << endl;
#     endif

      nb_step = (int)(max_value / step) + 1;

      for (i = 0;i < nb_state;i++) {
        frequency[i] = new double[nb_step];
        cumul[i] = new double[nb_step];

        value = 0.;
        frequency[i][0] = observation[i]->zero_probability;
        cumul[i][0] = observation[i]->zero_probability;
        dist_max[i] = frequency[i][0];

        if (observation[i]->zero_probability < 1.) {
          for (j = 1;j < nb_step;j++) {
            value += step;
            frequency[i][j] = (1. - observation[i]->zero_probability) * pdf(*gamma_dist[i] , value);
            if (frequency[i][j] > dist_max[i]) {
              dist_max[i] = frequency[i][j];
            }
            cumul[i][j] = observation[i]->zero_probability + (1. - observation[i]->zero_probability) *
                          cdf(*gamma_dist[i] , value);
          }
        }

        else {
          for (j = 1;j < nb_step;j++) {
            frequency[i][j] = 0.;
            cumul[i][j] = cumul[i][0];
          }
        }
      }
      break;
    }

    case GAUSSIAN : {
      gaussian_dist = new normal*[nb_state];

      for (i = 0;i < nb_state;i++) {
        gaussian_dist[i] = new normal(observation[i]->location , observation[i]->dispersion);
      }

      min_value = quantile(*gaussian_dist[0] , GAUSSIAN_TAIL);
      for (i = 1;i < nb_state;i++) {
        value = quantile(*gaussian_dist[i] , GAUSSIAN_TAIL);
        if (value < min_value) {
          min_value = value;
        }
      }
      if ((marginal_histogram) && (marginal_histogram->min_value - marginal_histogram->step < min_value)) {
        min_value = marginal_histogram->min_value - marginal_histogram->step;
      }
      if ((marginal_distribution) && (marginal_distribution->offset < min_value)) {
        min_value = marginal_distribution->offset;
      }

      max_value = quantile(complement(*gaussian_dist[0] , GAUSSIAN_TAIL));
      for (i = 1;i < nb_state;i++) {
        value = quantile(complement(*gaussian_dist[i] , GAUSSIAN_TAIL));
        if (value > max_value) {
          max_value = value;
        }
      }
      if ((marginal_histogram) && (marginal_histogram->max_value + marginal_histogram->step > max_value)) {
        max_value = marginal_histogram->max_value + marginal_histogram->step;
      }
      if ((marginal_distribution) && (marginal_distribution->nb_value - 1 > max_value)) {
        max_value = marginal_distribution->nb_value - 1;
      }

      if (step == D_DEFAULT) {
        step = (max_value - min_value) / GAUSSIAN_NB_STEP;
      }
      else {
        buff = (max_value - min_value) / GAUSSIAN_NB_STEP;
        if (buff < step) {
          step = buff;
        }
      }

      if (marginal_histogram) {
        buff = marginal_histogram->step / GAUSSIAN_NB_SUB_STEP;
      }
      else if (marginal_distribution) {
        buff = (double)min_interval / GAUSSIAN_NB_SUB_STEP;
      }
      if (((marginal_histogram) || (marginal_distribution)) && (buff < step)) {
        step = buff;
      }

#     ifdef DEBUG
      cout << "\nTest: " << min_value << " " << max_value << " | " << endl;
#     endif

      nb_step = (int)((max_value - min_value) / step) + 1;

      for (i = 0;i < nb_state;i++) {
        frequency[i] = new double[nb_step];
        cumul[i] = new double[nb_step];

        value = min_value;
        for (j = 0;j < nb_step;j++) {
          frequency[i][j] = pdf(*gaussian_dist[i] , value);
          cumul[i][j] = cdf(*gaussian_dist[i] , value);
          value += step;
        }

        dist_max[i] = frequency[i][(int)round((observation[i]->location - min_value) / step)];
      }
      break;
    }

    case VON_MISES : {
      min_value = 0.;

      switch (unit) {
      case DEGREE :
        max_value = 360;
        break;
      case RADIAN :
        max_value = 2 * M_PI;
        break;
      }

      if (step == D_DEFAULT) {
        step = max_value / VON_MISES_NB_STEP;
      }
      else {
        buff = max_value / VON_MISES_NB_STEP;
        if (buff < step) {
          step = buff;
        }
      }
      nb_step = max_value / step;

      for (i = 0;i < nb_state;i++) {
        frequency[i] = new double[nb_step];
        cumul[i] = new double[nb_step];

        value = 0.;

        switch (unit) {

        case DEGREE : {
          for (j = 0;j < nb_step;j++) {
            frequency[i][j] = exp(observation[i]->dispersion * cos((value - observation[i]->location) * M_PI / 180)) /
                              (360 * cyl_bessel_i(0 , observation[i]->dispersion));
            value += step;
          }
          break;
        }

        case RADIAN : {
          for (j = 0;j < nb_step;j++) {
            frequency[i][j] = exp(observation[i]->dispersion * cos(value - observation[i]->location)) /
                              (2 * M_PI * cyl_bessel_i(0 , observation[i]->dispersion));
            value += step;
          }
          break;
        }
        }

        dist_max[i] = frequency[i][(int)round(observation[i]->location / step)];

        value = 0.;
        cumul[i][0] = 0.;

        switch (unit) {

        case DEGREE : {
          for (j = 1;j < nb_step;j++) {
            cumul[i][j] = cumul[i][j - 1] + exp(observation[i]->dispersion * cos((value - observation[i]->location) * M_PI / 180)) * step /
                          (360 * cyl_bessel_i(0 , observation[i]->dispersion));
            value += step;
          }
          break;
        }

        case RADIAN : {
          for (j = 1;j < nb_step;j++) {
            cumul[i][j] = cumul[i][j - 1] + exp(observation[i]->dispersion * cos(value - observation[i]->location)) * step /
                          (2 * M_PI * cyl_bessel_i(0 , observation[i]->dispersion));
            value += step;
          }
          break;
        }
        }
      }
      break;
    }
    }

    if ((marginal_histogram) || (marginal_distribution)) {

      // calcul des melanges

      if (weight) {
        frequency[nb_state] = new double[nb_step];

        for (i = 0;i < nb_step;i++) {
          frequency[nb_state][i] = 0.;
          for (j = 0;j < nb_state;j++) {
            frequency[nb_state][i] += weight->mass[j] * frequency[j][i];
          }
        }
      }

      if (restoration_weight) {
        frequency[nb_state + 1] = new double[nb_step];

        for (i = 0;i < nb_step;i++) {
          frequency[nb_state + 1][i] = 0.;
          for (j = 0;j < nb_state;j++) {
            frequency[nb_state + 1][i] += restoration_weight->mass[j] * frequency[j][i];
          }
        }
      }
    }

    // ecriture des fichier de donnees

    value = min_value;
    for (i = 0;i < nb_step;i++) {
      if ((ident == ZERO_INFLATED_GAMMA) && (i == 0)) {
        for (j = 0;j < nb_state;j++) {
          if ((observation_histogram) && (observation_histogram[j]->nb_element > 0)) {
            scale[j] = observation_histogram[j]->nb_element;
          }
          else if ((observation_distribution) && (observation_distribution[j]->nb_element > 0)) {
            scale[j] = observation_distribution[j]->nb_element;
          }
          else {
            scale[j] = 1.;
          }
        }

        if (marginal_histogram) {
          scale[nb_state] = marginal_histogram->nb_element;
        }
        else if (marginal_distribution) {
          scale[nb_state] = marginal_distribution->nb_element;
        }
      }

      else if ((i == 0) || ((ident == ZERO_INFLATED_GAMMA) && (i == 1))) {
        for (j = 0;j < nb_state;j++) {
          if ((observation_histogram) && (observation_histogram[j]->nb_element > 0)) {
            scale[j] = observation_histogram[j]->step * observation_histogram[j]->nb_element;
          }
          else if ((observation_distribution) && (observation_distribution[j]->nb_element > 0)) {
            scale[j] = min_interval * observation_distribution[j]->nb_element;
          }
          else {
            scale[j] = 1.;
          }
        }

        if (marginal_histogram) {
          scale[nb_state] = marginal_histogram->step * marginal_histogram->nb_element;
        }
        else if (marginal_distribution) {
          scale[nb_state] = min_interval * marginal_distribution->nb_element;
        }
      }

      out_data_file << value;
      for (j = 0;j < nb_state;j++) {
        out_data_file << " " << frequency[j][i] * scale[j];
      }

      if ((marginal_histogram) || (marginal_distribution)) {
        if (weight) {
          for (j = 0;j < nb_state;j++) {
            out_data_file << " " << weight->mass[j] * frequency[j][i] * scale[nb_state];
          }
          out_data_file << " " << frequency[nb_state][i] * scale[nb_state];
        }

        if (restoration_weight) {
          for (j = 0;j < nb_state;j++) {
            out_data_file << " " << restoration_weight->mass[j] * frequency[j][i] *
                                    scale[nb_state];
          }
          out_data_file << " " << frequency[nb_state + 1][i] * scale[nb_state];
        }
      }
      out_data_file << endl;

      value += step;
    }

    if (observation_histogram) {
      for (i = 0;i < nb_state;i++) {
        if (observation_histogram[i]->nb_element > 0) {
          data_file_name[i + 1] << prefix << process << i + 1 << ".dat";
          observation_histogram[i]->plot_print((data_file_name[i + 1].str()).c_str());
        }
      }
    }

    if (observation_distribution) {
      for (i = 0;i < nb_state;i++) {
        if (observation_distribution[i]->nb_element > 0) {
          data_file_name[i + 1] << prefix << process << i + 1 << ".dat";
          observation_distribution[i]->plot_print((data_file_name[i + 1].str()).c_str());
        }
      }
    }

    if (marginal_histogram) {
      data_file_name[nb_state + 1] << prefix << process << nb_state + 1 << ".dat";
      marginal_histogram->plot_print((data_file_name[nb_state + 1].str()).c_str());
    }

    if (marginal_distribution) {
      data_file_name[nb_state + 1] << prefix << process << nb_state + 1 << ".dat";
      marginal_distribution->plot_print((data_file_name[nb_state + 1].str()).c_str());
    }

    if (empirical_cdf) {
      if (weight) {
        cumul[nb_state] = new double[nb_step];

        for (j = 0;j < nb_step;j++) {
          cumul[nb_state][j] = 0.;
          for (k = 0;k < nb_state;k++) {
            cumul[nb_state][j] += weight->mass[k] * cumul[k][j];
          }
        }

#       ifdef DEBUG
        cout << "\nCumul: ";
        buff = 0;
        for (j = 0;j < nb_step;j++) {
          buff += frequency[nb_state][j] * step;
          cout << buff << " " << cumul[nb_state][j] << " | ";
        }
        cout << endl;
#       endif

        qqplot = q_q_plot_computation(min_value , step , nb_step , cumul[nb_state] ,
                                      nb_value - 1 , empirical_cdf);

        data_file_name[nb_state + 2] << prefix << process << nb_state + 2 << ".dat";
        q_q_plot_print((data_file_name[nb_state + 2].str()).c_str() , nb_value - 1 , qqplot);
      }

      if (restoration_weight) {
        cumul[nb_state + 1] = new double[nb_step];

        for (j = 0;j < nb_step;j++) {
          cumul[nb_state + 1][j] = 0.;
          for (k = 0;k < nb_state;k++) {
            cumul[nb_state + 1][j] += restoration_weight->mass[k] * cumul[k][j];
          }
        }

        qqplot = q_q_plot_computation(min_value , step , nb_step , cumul[nb_state + 1] ,
                                      nb_value - 1 , empirical_cdf);

        data_file_name[nb_state + 3] << prefix << process << nb_state + 3 << ".dat";
        q_q_plot_print((data_file_name[nb_state + 3].str()).c_str() , nb_value - 1 , qqplot);
      }
    }

    // ecriture du fichiers de commandes et du fichiers d'impression

    for (i = 0;i < 2;i++) {
      ostringstream file_name[2];

      switch (i) {
      case 0 :
        file_name[0] << prefix << process << ".plot";
        break;
      case 1 :
        file_name[0] << prefix << process << ".print";
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
      out_file << STAT_label[STATL_OUTPUT_PROCESS] << " " << process << "\"\n\n";

      if ((observation_histogram) || (observation_distribution)) {
        for (j = 0;j < nb_state;j++) {
          if ((observation_histogram) && (observation_histogram[j]->nb_element > 0)) {
            max = (int)(MAX(observation_histogram[j]->max ,
                            dist_max[j] * scale[j]) * YSCALE) + 1;
          }
          else if ((observation_distribution) && (observation_distribution[j]->nb_element > 0)) {
            max = (int)(MAX(observation_distribution[j]->max ,
                            dist_max[j] * scale[j]) * YSCALE) + 1;
          }
          else {
            max = MIN(dist_max[j] * YSCALE , 1.);
          }

          out_file << "plot [" << min_value << ":" << max_value << "] [0:"
                   << max << "] "; 
          if (((observation_histogram) && (observation_histogram[j]->nb_element > 0)) ||
              ((observation_distribution) && (observation_distribution[j]->nb_element > 0))) {
            out_file << "\"" << label((data_file_name[j + 1].str()).c_str()) << "\"";
            if ((observation_histogram) && (observation_histogram[j]->nb_element > 0)) {
              out_file << " using 1:2";
            }
            out_file << " title \"";
            switch (model) {
            case MIXTURE :
              out_file << STAT_label[STATL_COMPONENT];
              break;
            case HIDDEN_MARKOV :
              out_file << STAT_label[STATL_STATE];
              break;
            }
            out_file << " " << j << " " << STAT_label[STATL_OBSERVATION] << " ";
            if ((observation_histogram) && (observation_histogram[j]->nb_element > 0)) {
              out_file << STAT_label[STATL_HISTOGRAM] << "\" with histeps,\\" << endl;
            }
            else {
              out_file << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
                       << "\" with impulses,\\" << endl;
            }
          }
          out_file << "\"" << label((data_file_name[0].str()).c_str())
                   << "\" using 1:" << j + 2 << " title \"";
          switch (model) {
          case MIXTURE :
            out_file << STAT_label[STATL_COMPONENT];
            break;
          case HIDDEN_MARKOV :
            out_file << STAT_label[STATL_STATE];
            break;
          }
          out_file << " " << j << " " << STAT_label[STATL_OBSERVATION] << " " << STAT_label[STATL_DISTRIBUTION] << " ";
          observation[j]->plot_title_print(out_file);
          out_file << "\" with lines" << endl;

          if ((i == 0) && (j < nb_state - 1)) {
            out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
          }
          out_file << endl;
        }
      }

      else {
        max = dist_max[0];
        for (j = 1;j < nb_state;j++) {
          if (dist_max[j] > max) {
            max = dist_max[j];
          }
        }

        out_file << "plot [" << min_value  << ":" << max_value << "] [0:"
                 << MIN(max * YSCALE , 1.) << "] ";
        for (j = 0;j < nb_state;j++) {
          out_file << "\"" << label((data_file_name[0].str()).c_str()) << "\" using 1:" << j + 2
                   << " title \"";
          switch (model) {
          case MIXTURE :
            out_file << STAT_label[STATL_COMPONENT];
            break;
          case HIDDEN_MARKOV :
            out_file << STAT_label[STATL_STATE];
            break;
          }
          out_file << " " << j << " " << STAT_label[STATL_OBSERVATION] << " " << STAT_label[STATL_DISTRIBUTION] << " ";
          observation[j]->plot_title_print(out_file);
          out_file << "\" with lines";
          if (j < nb_state - 1) {
            out_file << ",\\";
          }
          out_file << endl;
        }
      }

      if ((marginal_histogram) || (marginal_distribution)) {
        if (weight) {
          if (marginal_histogram) {
            max = marginal_histogram->max;
           }
           else {
            max = marginal_distribution->max;
          }

          for (j = 0;j < nb_step;j++) {
            if (frequency[nb_state][i] * scale[nb_state] > max) {
              max = frequency[nb_state][i] * scale[nb_state];
            }
          }

          if (i == 0) {
            out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
          }
          out_file << endl;

          out_file << "set title \"";
          if (title) {
            out_file << title << " - ";
          }
          out_file << STAT_label[STATL_OUTPUT_PROCESS] << " " << process << " - "
                   << STAT_label[STATL_THEORETICAL] << " " << STAT_label[STATL_WEIGHTS] << "\"\n\n";

          out_file << "plot [" << min_value << ":" << max_value << "] [0:"
                   << max * YSCALE << "] \""
                   << label((data_file_name[nb_state + 1].str()).c_str()) << "\"";
          if (marginal_histogram) {
            out_file << " using 1:2";
          }
          out_file << " title \""  << STAT_label[STATL_MARGINAL] << " ";
          if (marginal_histogram) {
            out_file << STAT_label[STATL_HISTOGRAM] << "\" with histeps,\\" << endl;
          }
          else {
            out_file << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
                     << "\" with impulses,\\" << endl;
          }
          for (j = 0;j < nb_state;j++) {
            out_file << "\"" << label((data_file_name[0].str()).c_str())
                     << "\" using 1:" << nb_state + j + 2 << " title \"";
            switch (model) {
            case MIXTURE :
              out_file << STAT_label[STATL_COMPONENT];
              break;
            case HIDDEN_MARKOV :
              out_file << STAT_label[STATL_STATE];
              break;
            }
            out_file << " " << j << " " << STAT_label[STATL_OBSERVATION] << " " << STAT_label[STATL_DISTRIBUTION] << " ";
            observation[j]->plot_title_print(out_file);
            out_file << "\" with lines,\\" << endl;
          }
          out_file << "\"" << label((data_file_name[0].str()).c_str())
                   << "\" using 1:" << 2 * nb_state + 2
                   << " title \"" << STAT_label[STATL_MIXTURE]<< "\" with lines" << endl;

          if (i == 0) {
            out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
          }
          out_file << endl;

          out_file << "set grid" << endl;

          out_file << "plot [" << min_value << ":" << max_value << "] ["
                   << min_value << ":" << max_value << "] \""
                   << label((data_file_name[nb_state + 2].str()).c_str())
                   << "\" using 1:2 notitle with points" << endl;

          out_file << "unset grid" << endl;
        }

        if (restoration_weight) {
          if (marginal_histogram) {
            max = marginal_histogram->max;
           }
           else {
            max = marginal_distribution->max;
          }

          for (j = 0;j < nb_step;j++) {
            if (frequency[nb_state + 1][i] * scale[nb_state] > max) {
              max = frequency[nb_state + 1][i] * scale[nb_state];
            }
          }

          if (i == 0) {
            out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
          }
          out_file << endl;

          out_file << "set title \"";
          if (title) {
            out_file << title << " - ";
          }
          out_file << STAT_label[STATL_OUTPUT_PROCESS] << " " << process << " - "
                   << STAT_label[STATL_RESTORATION] << " " << STAT_label[STATL_WEIGHTS] << "\"\n\n";

          out_file << "plot [" << min_value << ":" << max_value << "] [0:"
                   << max * YSCALE << "] \""
                   << label((data_file_name[nb_state + 1].str()).c_str()) << "\"";
          if (marginal_histogram) {
            out_file << " using 1:2";
          }
          out_file << " title \""  << STAT_label[STATL_MARGINAL] << " ";
          if (marginal_histogram) {
            out_file << STAT_label[STATL_HISTOGRAM] << "\" with histeps,\\" << endl;
          }
          else {
            out_file << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
                     << "\" with impulses,\\" << endl;
          }

          if (weight) {
            dist_index = 2 * nb_state + 3;
          }
          else {
            dist_index = nb_state + 2;
          }

          for (j = 0;j < nb_state;j++) {
            out_file << "\"" << label((data_file_name[0].str()).c_str())
                     << "\" using 1:" << dist_index + j << " title \"";
            switch (model) {
            case MIXTURE :
              out_file << STAT_label[STATL_COMPONENT];
              break;
            case HIDDEN_MARKOV :
              out_file << STAT_label[STATL_STATE];
              break;
            }
            out_file << " " << j << " " << STAT_label[STATL_OBSERVATION] << " " << STAT_label[STATL_DISTRIBUTION] << " ";
            observation[j]->plot_title_print(out_file);
            out_file << "\" with lines,\\" << endl;
          }
          out_file << "\"" << label((data_file_name[0].str()).c_str())
                   << "\" using 1:" << dist_index + nb_state
                   << " title \"" << STAT_label[STATL_MIXTURE]<< "\" with lines" << endl;

          if (i == 0) {
            out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
          }
          out_file << endl;

          out_file << "set grid" << endl;

          out_file << "plot [" << min_value << ":" << max_value << "] ["
                   << min_value << ":" << max_value << "] \""
                   << label((data_file_name[nb_state + 3].str()).c_str())
                   << "\" using 1:2 notitle with points" << endl;

          out_file << "unset grid" << endl;
        }
      }

      if (i == 1) {
        out_file << "\nset terminal x11" << endl;
      }

      out_file << "\npause 0 \"" << STAT_label[STATL_END] << "\"" << endl;
    }

    if ((ident == GAMMA) || (ident == ZERO_INFLATED_GAMMA)) {
      for (i = 0;i < nb_state;i++) {
        delete gamma_dist[i];
      }
      delete [] gamma_dist;
    }

    else if (ident == GAUSSIAN) {
      for (i = 0;i < nb_state;i++) {
        delete gaussian_dist[i];
      }
      delete [] gaussian_dist;
    }

    for (i = 0;i < nb_state + 2;i++) {
      delete [] frequency[i];
    }
    delete [] frequency;

    for (i = 0;i < nb_state + 2;i++) {
      delete [] cumul[i];
    }
    delete [] cumul;

    delete [] scale;
    delete [] dist_max;
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Sortie graphique d'un q-q plot.
 *
 *  arguments : reference sur un objet SinglePlot,
 *              nombre de valeurs, pointeur sur le q-q plot.
 *
 *--------------------------------------------------------------*/

void q_q_plotable_write(SinglePlot &plot , int nb_value , double **qqplot)

{
  register int i;


  for (i = 0;i < nb_value;i++) {
    plot.add_point(qqplot[0][i] , qqplot[1][i]);
  }

  delete [] qqplot[0];
  delete [] qqplot[1];
  delete [] qqplot;
}


/*--------------------------------------------------------------*
 *
 *  Sortie graphique d'un objet ContinuousParametricProcess.
 *
 *  arguments : reference sur un objet MultiPlotSet, indice du MultiPlot,
 *              indice du processus d'observation,
 *              pointeurs sur les histogrammes d'observation ou
 *              les lois d'observation empiriques et sur l'histogramme marginale ou
 *              la loi marginale empirique, nombre de valeurs,
 *              pointeur sur la fonction de repartition empirique, type de modele.
 *
 *--------------------------------------------------------------*/

void ContinuousParametricProcess::plotable_write(MultiPlotSet &plot , int &index , int process ,
                                                 Histogram **observation_histogram ,
                                                 FrequencyDistribution **observation_distribution ,
                                                 Histogram *marginal_histogram ,
                                                 FrequencyDistribution *marginal_distribution ,
                                                 int nb_value , double **empirical_cdf , int model) const

{
  register int i , j , k;
  int min_interval , nb_step;
  double value , min_value , max_value , step , buff , scale , max , *dist_max ,
         **cumul , **frequency , **qqplot;
  gamma_distribution<double> **gamma_dist;
  normal **gaussian_dist;
  ostringstream title , legend;


  frequency = new double*[nb_state + 1];
  cumul = new double*[nb_state + 1];

  dist_max = new double[nb_state];

  // calcul des lois d'observation discretisees

/*  if (empirical_cdf) {
    step = empirical_cdf[0][1] - empirical_cdf[0][0];
    for (i = 2;i < nb_value;i++) {
      if (empirical_cdf[0][i] - empirical_cdf[0][i - 1] < step) {
        step = empirical_cdf[0][i] - empirical_cdf[0][i - 1];
      }
    }
  }
  else { */
    step = D_DEFAULT;
//  }

  if (marginal_distribution) {
    min_interval = marginal_distribution->min_interval_computation();
  }

  switch (ident) {

  case GAMMA : {
    gamma_dist = new gamma_distribution<double>*[nb_state];

    min_value = 0.;
    max_value = 0.;

    for (i = 0;i < nb_state;i++) {
      if (observation[i]->shape > 0.) {
        gamma_dist[i] = new gamma_distribution<double>(observation[i]->shape , observation[i]->scale);

        value = quantile(complement(*gamma_dist[i] , GAMMA_TAIL));
        if (value > max_value) {
          max_value = value;
        }
      }

      else {
        gamma_dist[i] = NULL;
      }
    }

    if ((marginal_histogram) && (marginal_histogram->max_value + marginal_histogram->step > max_value)) {
      max_value = marginal_histogram->max_value + marginal_histogram->step;
    }
    if ((marginal_distribution) && (marginal_distribution->nb_value - 1 > max_value)) {
      max_value = marginal_distribution->nb_value - 1;
    }

    if (step == D_DEFAULT) {
      step = max_value / GAMMA_NB_STEP;
    }
    else {
      buff = max_value / GAMMA_NB_STEP;
      if (buff < step) {
        step = buff;
      }
    }

    if (marginal_histogram) {
      buff = marginal_histogram->step / GAMMA_NB_SUB_STEP;
    }
    else if (marginal_distribution) {
      buff = (double)min_interval / GAMMA_NB_SUB_STEP;
    }
    if (((marginal_histogram) || (marginal_distribution)) && (buff < step)) {
      step = buff;
    }

#   ifdef DEBUG
    cout << "\nTest: " << max_value << " | " << step << endl;
#   endif

    nb_step = (int)(max_value / step) + 1;

    for (i = 0;i < nb_state;i++) {
      frequency[i] = new double[nb_step];
      cumul[i] = new double[nb_step];

      dist_max[i] = 0.;
      value = 0.;

      if (observation[i]->shape == 0.) {
        frequency[i][0] = 1.;
        cumul[i][0] = frequency[i][0];

        for (j = 1;j < nb_step;j++) {
          frequency[i][j] = 0.;
          cumul[i][j] = cumul[i][0];
        }
      }

      else {
        for (j = 0;j < nb_step;j++) {
          frequency[i][j] = pdf(*gamma_dist[i] , value);
          if (frequency[i][j] > dist_max[i]) {
            dist_max[i] = frequency[i][j];
          }
          cumul[i][j] = cdf(*gamma_dist[i] , value);
          value += step;
        }
      }
    }
    break;
  }

  case ZERO_INFLATED_GAMMA : {
    gamma_dist = new gamma_distribution<double>*[nb_state];

    min_value = 0.;
    max_value = 0.;

    for (i = 0;i < nb_state;i++) {
      if (observation[i]->zero_probability < 1.) {
        gamma_dist[i] = new gamma_distribution<double>(observation[i]->shape , observation[i]->scale);

        value = quantile(complement(*gamma_dist[i] , GAMMA_TAIL));
        if (value > max_value) {
          max_value = value;
        }
      }

      else {
        gamma_dist[i] = NULL;
      }
    }

    if ((marginal_histogram) && (marginal_histogram->max_value + marginal_histogram->step > max_value)) {
      max_value = marginal_histogram->max_value + marginal_histogram->step;
    }
    if ((marginal_distribution) && (marginal_distribution->nb_value - 1 > max_value)) {
      max_value = marginal_distribution->nb_value - 1;
    }

    if (step == D_DEFAULT) {
      step = max_value / GAMMA_NB_STEP;
    }
    else {
      buff = max_value / GAMMA_NB_STEP;
      if (buff < step) {
        step = buff;
      }
    }

    if (marginal_histogram) {
      buff = marginal_histogram->step / GAMMA_NB_SUB_STEP;
    }
    else if (marginal_distribution) {
      buff = (double)min_interval / GAMMA_NB_SUB_STEP;
    }
    if (((marginal_histogram) || (marginal_distribution)) && (buff < step)) {
      step = buff;
    }

#   ifdef DEBUG
    cout << "\nTest: " << max_value << " | " << step << endl;
#   endif

    nb_step = (int)(max_value / step) + 1;

    for (i = 0;i < nb_state;i++) {
      frequency[i] = new double[nb_step];
      cumul[i] = new double[nb_step];

      value = 0.;
      frequency[i][0] = observation[i]->zero_probability;
      cumul[i][0] = observation[i]->zero_probability;
      dist_max[i] = frequency[i][0];

      if (observation[i]->zero_probability < 1.) {
        for (j = 1;j < nb_step;j++) {
          value += step;
          frequency[i][j] = (1. - observation[i]->zero_probability) * pdf(*gamma_dist[i] , value);
          if (frequency[i][j] > dist_max[i]) {
            dist_max[i] = frequency[i][j];
          }
          cumul[i][j] = observation[i]->zero_probability + (1. - observation[i]->zero_probability) *
                        cdf(*gamma_dist[i] , value);
        }
      }

      else {
        for (j = 1;j < nb_step;j++) {
          frequency[i][j] = 0.;
          cumul[i][j] = cumul[i][0];
        }
      }
    }
    break;
  }

  case GAUSSIAN : {
    gaussian_dist = new normal*[nb_state];

    for (i = 0;i < nb_state;i++) {
      gaussian_dist[i] = new normal(observation[i]->location , observation[i]->dispersion);
    }

    min_value = quantile(*gaussian_dist[0] , GAUSSIAN_TAIL);
    for (i = 1;i < nb_state;i++) {
      value = quantile(*gaussian_dist[i] , GAUSSIAN_TAIL);
      if (value < min_value) {
        min_value = value;
      }
    }
    if ((marginal_histogram) && (marginal_histogram->min_value - marginal_histogram->step < min_value)) {
      min_value = marginal_histogram->min_value - marginal_histogram->step;
    }
    if ((marginal_distribution) && (marginal_distribution->offset < min_value)) {
      min_value = marginal_distribution->offset;
    }

//    max_value = quantile(*gaussian_dist[0] , 1. - GAUSSIAN_TAIL);
    max_value = quantile(complement(*gaussian_dist[0] , GAUSSIAN_TAIL));
    for (i = 1;i < nb_state;i++) {
      value = quantile(complement(*gaussian_dist[i] , GAUSSIAN_TAIL));
      if (value > max_value) {
        max_value = value;
      }
    }
    if ((marginal_histogram) && (marginal_histogram->max_value + marginal_histogram->step > max_value)) {
      max_value = marginal_histogram->max_value + marginal_histogram->step;
    }
    if ((marginal_distribution) && (marginal_distribution->nb_value - 1 > max_value)) {
      max_value = marginal_distribution->nb_value - 1;
    }

    if (step == D_DEFAULT) {
      step = (max_value - min_value) / GAUSSIAN_NB_STEP;
    }
    else {
      buff = (max_value - min_value) / GAUSSIAN_NB_STEP;
      if (buff < step) {
        step = buff;
      }
    }

    if (marginal_histogram) {
      buff = marginal_histogram->step / GAUSSIAN_NB_SUB_STEP;
    }
    else if (marginal_distribution) {
      buff = (double)min_interval / GAUSSIAN_NB_SUB_STEP;
    }
    if (((marginal_histogram) || (marginal_distribution)) && (buff < step)) {
      step = buff;
    }

    nb_step = (int)((max_value - min_value) / step) + 1;

    for (i = 0;i < nb_state;i++) {
      frequency[i] = new double[nb_step];
      cumul[i] = new double[nb_step];

      value = min_value;
      for (j = 0;j < nb_step;j++) {
        frequency[i][j] = pdf(*gaussian_dist[i] , value);
        cumul[i][j] = cdf(*gaussian_dist[i] , value);
        value += step;
      }

      dist_max[i] = frequency[i][(int)round((observation[i]->location - min_value) / step)];
    }
    break;
  }

  case VON_MISES : {
    min_value = 0.;

    switch (unit) {
    case DEGREE :
      max_value = 360;
      break;
    case RADIAN :
      max_value = 2 * M_PI;
      break;
    }

    if (step == D_DEFAULT) {
      step = max_value / VON_MISES_NB_STEP;
    }
    else {
      buff = max_value / VON_MISES_NB_STEP;
      if (buff < step) {
        step = buff;
      }
    }
    nb_step = max_value / step;

    for (i = 0;i < nb_state;i++) {
      frequency[i] = new double[nb_step];
      cumul[i] = new double[nb_step];

      value = 0.;

      switch (unit) {

      case DEGREE : {
        for (j = 0;j < nb_step;j++) {
          frequency[i][j] = exp(observation[i]->dispersion * cos((value - observation[i]->location) * M_PI / 180)) /
                            (360 * cyl_bessel_i(0 , observation[i]->dispersion));
          value += step;
        }
        break;
      }

      case RADIAN : {
        for (j = 0;j < nb_step;j++) {
          frequency[i][j] = exp(observation[i]->dispersion * cos(value - observation[i]->location)) /
                            (2 * M_PI * cyl_bessel_i(0 , observation[i]->dispersion));
          value += step;
        }
        break;
      }
      }

      dist_max[i] = frequency[i][(int)round(observation[i]->location / step)];

      value = 0.;
      cumul[i][0] = 0.;

      switch (unit) {

      case DEGREE : {
        for (j = 1;j < nb_step;j++) {
          cumul[i][j] = cumul[i][j - 1] + exp(observation[i]->dispersion * cos((value - observation[i]->location) * M_PI / 180)) * step /
                        (360 * cyl_bessel_i(0 , observation[i]->dispersion));
          value += step;
        }
        break;
      }

      case RADIAN : {
        for (j = 1;j < nb_step;j++) {
          cumul[i][j] = cumul[i][j - 1] + exp(observation[i]->dispersion * cos(value - observation[i]->location)) * step /
                       (2 * M_PI * cyl_bessel_i(0 , observation[i]->dispersion));
          value += step;
        }
        break;
      }
      }
    }
    break;
  }
  }

  plot.variable_nb_viewpoint[process] = 1;

  if ((observation_histogram) || (observation_distribution)) {
    for (i = 0;i < nb_state;i++) {

      // vue : ajustement loi d'observation

      plot.variable[index] = process;
//      plot.viewpoint[index] = OBSERVATION;

      title.str("");
      title << STAT_label[STATL_OUTPUT_PROCESS] << " " << process;
      plot[index].title = title.str();

      plot[index].xrange = Range(min_value , max_value);

      if ((observation_histogram) && (observation_histogram[i]->nb_element > 0)) {
        scale = observation_histogram[i]->step * observation_histogram[i]->nb_element;
        max = (int)(MAX(observation_histogram[i]->max ,
                        dist_max[i] * scale) * YSCALE) + 1;
      }
      else if ((observation_distribution) && (observation_distribution[i]->nb_element > 0)) {
        scale = min_interval * observation_distribution[i]->nb_element;
        max = (int)(MAX(observation_distribution[i]->max ,
                        dist_max[i] * scale) * YSCALE) + 1;
      }
      else {
        scale = 1.;
        max = MIN(dist_max[i] * YSCALE , 1.);
      }

      plot[index].yrange = Range(0 , max);

      if (((observation_histogram) && (observation_histogram[j]->nb_element > 0)) ||
          ((observation_distribution) && (observation_distribution[j]->nb_element > 0))) {
        plot[index].resize(2);

        legend.str("");
        switch (model) {
        case MIXTURE :
          legend << STAT_label[STATL_COMPONENT];
          break;
        case HIDDEN_MARKOV :
          legend << STAT_label[STATL_STATE];
          break;
        }
        legend << " " << i << " " << STAT_label[STATL_OBSERVATION] << " ";
        if ((observation_histogram) && (observation_histogram[j]->nb_element > 0)) {
          legend << STAT_label[STATL_HISTOGRAM];
        }
        else {
          legend << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
        }
        plot[index][0].legend = legend.str();

        if ((observation_histogram) && (observation_histogram[j]->nb_element > 0)) {
          plot[index][0].style = "histeps";
          observation_histogram[i]->plotable_write(plot[index][0]);
        }
        else {
          plot[index][0].style = "impulses";
          observation_distribution[i]->plotable_frequency_write(plot[index][0]);
        }

        j = 1;
      }

      else {
        plot[index].resize(1);

        j = 0;
      }

      legend.str("");
      switch (model) {
      case MIXTURE :
        legend << STAT_label[STATL_COMPONENT];
        break;
      case HIDDEN_MARKOV :
        legend << STAT_label[STATL_STATE];
        break;
      }
      legend << " " << i << " " << STAT_label[STATL_OBSERVATION] << " " << STAT_label[STATL_DISTRIBUTION];
      observation[i]->plot_title_print(legend);
      plot[index][j].legend = legend.str();

      plot[index][j].style = "lines";

      value = min_value;
      for (k = 0;k < nb_step;k++) {
        plot[index][j].add_point(value , frequency[i][k] * scale);
        value += step;
      }

      index++;
    }
  }

  else {
    title.str("");
    title << STAT_label[STATL_OUTPUT_PROCESS] << " " << process;
    plot[index].title = title.str();

    plot[index].xrange = Range(min_value , max_value);

    // vue : lois d'observation

    max = dist_max[0];
    for (i = 1;i < nb_state;i++) {
      if (dist_max[i] > max) {
        max = dist_max[i];
      }
    }
    plot[index].yrange = Range(0 , MIN(max * YSCALE , 1.));

    for (i = 0;i < nb_state;i++) {
      legend.str("");
      switch (model) {
      case MIXTURE :
        legend << STAT_label[STATL_COMPONENT];
        break;
      case HIDDEN_MARKOV :
        legend << STAT_label[STATL_STATE];
        break;
      }
      legend << " " << i << " " << STAT_label[STATL_OBSERVATION] << " " << STAT_label[STATL_DISTRIBUTION];
      observation[i]->plot_title_print(legend);
      plot[index][i].legend = legend.str();

      plot[index][i].style = "lines";

      value = min_value;
      for (j = 0;j < nb_step;j++) {
        plot[index][j].add_point(value , frequency[i][j]);
        value += step;
      }
    }
    index++;
  }

  if ((marginal_histogram) || (marginal_distribution)) {
    frequency[nb_state] = new double[nb_step];
    cumul[nb_state] = new double[nb_step];

    if (marginal_histogram) {
      scale = marginal_histogram->step * marginal_histogram->nb_element;
    }
    else {
      scale = min_interval * marginal_distribution->nb_element;
    }

    if (weight) {

      // vue : ajustement melange de lois d'observation - poids theoriques

      for (i = 0;i < nb_step;i++) {
        frequency[nb_state][i] = 0.;
        for (j = 0;j < nb_state;j++) {
          frequency[nb_state][i] += weight->mass[j] * frequency[j][i];
        }
      }

      title.str("");
      title << STAT_label[STATL_OUTPUT_PROCESS] << " " << process << " - "
            << STAT_label[STATL_THEORETICAL] << " " << STAT_label[STATL_WEIGHTS];
      plot[index].title = title.str();

      plot[index].xrange = Range(min_value , max_value);

      if (marginal_histogram) {
        max = marginal_histogram->max;
      }
      else {
        max = marginal_distribution->max;
      }

      for (i = 0;i < nb_step;i++) {
        if (frequency[nb_state][i] * scale > max) {
          max = frequency[nb_state][i] * scale;
        }
      }

      plot[index].yrange = Range(0 , max * YSCALE);

      plot[index].resize(nb_state + 2);

      legend.str("");
      legend << STAT_label[STATL_MARGINAL] << " ";
      if (marginal_histogram) {
        legend << STAT_label[STATL_HISTOGRAM];
      }
      else {
        legend << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
      }
      plot[index][0].legend = legend.str();

      if (marginal_histogram) {
        plot[index][0].style = "histeps";
        marginal_histogram->plotable_write(plot[index][0]);
      }
      else {
        plot[index][0].style = "impulses";
        marginal_distribution->plotable_frequency_write(plot[index][0]);
      }

      for (i = 0;i < nb_state;i++) {
        legend.str("");
        switch (model) {
        case MIXTURE :
          legend << STAT_label[STATL_COMPONENT];
          break;
        case HIDDEN_MARKOV :
          legend << STAT_label[STATL_STATE];
          break;
        }
        legend << " " << i << " " << STAT_label[STATL_OBSERVATION] << " " << STAT_label[STATL_DISTRIBUTION];
        observation[i]->plot_title_print(legend);
        plot[index][i + 1].legend = legend.str();

        plot[index][i + 1].style = "lines";

        value = min_value;
        for (j = 0;j < nb_step;j++) {
          plot[index][i + 1].add_point(value , weight->mass[i] * frequency[i][j] * scale);
          value += step;
        }
      }

      plot[index][nb_state + 2].legend = STAT_label[STATL_MIXTURE];

      plot[index][nb_state + 2].style = "lines";

      value = min_value;
      for (i = 0;i < nb_step;i++) {
        plot[index][nb_state + 2].add_point(value , frequency[nb_state][i] * scale);
        value += step;
      }

      index++;

      // vue : q-q plot

      for (i = 0;i < nb_step;i++) {
        cumul[nb_state][i] = 0.;
        for (j = 0;j < nb_state;j++) {
          cumul[nb_state][i] += weight->mass[j] * cumul[j][i];
        }
      }

#     ifdef DEBUG
      cout << "\nCumul: ";
      buff = 0;
      for (i = 0;i < nb_step;i++) {
        buff += frequency[nb_state][i] * step;
        cout << buff << " " << cumul[nb_state][i] << " | ";
      }
      cout << endl;
#     endif

      title.str("");
      title << STAT_label[STATL_OUTPUT_PROCESS] << " " << process << " - "
            << STAT_label[STATL_THEORETICAL] << " " << STAT_label[STATL_WEIGHTS];
      plot[index].title = title.str();

      plot[index].xrange = Range(min_value , max_value);
      plot[index].yrange = Range(min_value , max_value);

      plot[index].resize(2);

      plot[index][0].style = "lines";
      plot[index][0].add_point(min_value , min_value);
      plot[index][0].add_point(max_value , max_value);

      legend.str("");
      legend << STAT_label[STATL_MIXTURE] << " " << STAT_label[STATL_Q_Q_PLOT];
      plot[index][1].legend = legend.str();

      plot[index][1].style = "points";

      qqplot = q_q_plot_computation(min_value , step , nb_step , cumul[nb_state] ,
                                    nb_value - 1 , empirical_cdf);
      q_q_plotable_write(plot[index][1] , nb_value - 1 , qqplot);
      index++;
    }

    if (restoration_weight) {

      // vue : ajustement melange de lois d'observation - poids deduits de la restauration

      for (i = 0;i < nb_step;i++) {
        frequency[nb_state][i] = 0.;
        for (j = 0;j < nb_state;j++) {
          frequency[nb_state][i] += restoration_weight->mass[j] * frequency[j][i];
        }
      }

      title.str("");
      title << STAT_label[STATL_OUTPUT_PROCESS] << " " << process << " - "
            << STAT_label[STATL_RESTORATION] << " " << STAT_label[STATL_WEIGHTS];
      plot[index].title = title.str();

      plot[index].xrange = Range(min_value , max_value);

      if (marginal_histogram) {
        max = marginal_histogram->max;
       }
       else {
        max = marginal_distribution->max;
      }

      for (i = 0;i < nb_step;i++) {
        if (frequency[nb_state][i] * scale > max) {
          max = frequency[nb_state][i] * scale;
        }
      }

      plot[index].yrange = Range(0 , max * YSCALE);

      plot[index].resize(nb_state + 2);

      legend.str("");
      legend << STAT_label[STATL_MARGINAL] << " ";
      if (marginal_histogram) {
        legend << STAT_label[STATL_HISTOGRAM];
      }
      else {
        legend << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
      }
      plot[index][0].legend = legend.str();

      if (marginal_histogram) {
        plot[index][0].style = "histeps";
        marginal_histogram->plotable_write(plot[index][0]);
      }
      else {
        plot[index][0].style = "impulses";
        marginal_distribution->plotable_frequency_write(plot[index][0]);
      }

      for (i = 0;i < nb_state;i++) {
        legend.str("");
        switch (model) {
        case MIXTURE :
          legend << STAT_label[STATL_COMPONENT];
          break;
        case HIDDEN_MARKOV :
          legend << STAT_label[STATL_STATE];
          break;
        }
        legend << " " << i << " " << STAT_label[STATL_OBSERVATION] << " " << STAT_label[STATL_DISTRIBUTION];
        observation[i]->plot_title_print(legend);
        plot[index][i + 1].legend = legend.str();

        plot[index][i + 1].style = "lines";

        value = min_value;
        for (j = 0;j < nb_step;j++) {
          plot[index][i + 1].add_point(value , restoration_weight->mass[i] * frequency[i][j] * scale);
          value += step;
        }
      }

      plot[index][nb_state + 2].legend = STAT_label[STATL_MIXTURE];

      plot[index][nb_state + 2].style = "lines";

      value = min_value + step / 2;
      for (i = 0;i < nb_step;i++) {
        plot[index][nb_state + 2].add_point(value , frequency[nb_state][i] * scale);
        value += step;
      }

      index++;

      // vue : q-q plot

      for (i = 0;i < nb_step;i++) {
        cumul[nb_state][i] = 0.;
        for (j = 0;j < nb_state;j++) {
          cumul[nb_state][i] += restoration_weight->mass[j] * cumul[j][i];
        }
      }

      title.str("");
      title << STAT_label[STATL_OUTPUT_PROCESS] << " " << process << " - "
            << STAT_label[STATL_RESTORATION] << " " << STAT_label[STATL_WEIGHTS];
      plot[index].title = title.str();

      plot[index].xrange = Range(min_value , max_value);
      plot[index].yrange = Range(min_value , max_value);

      plot[index].resize(2);

      plot[index][0].style = "lines";
      plot[index][0].add_point(min_value , min_value);
      plot[index][0].add_point(max_value , max_value);

      legend.str("");
      legend << STAT_label[STATL_MIXTURE] << " " << STAT_label[STATL_Q_Q_PLOT];
      plot[index][1].legend = legend.str();

      plot[index][1].style = "points";

      qqplot = q_q_plot_computation(min_value , step , nb_step , cumul[nb_state] ,
                                    nb_value - 1 , empirical_cdf);
      q_q_plotable_write(plot[index][1] , nb_value - 1 , qqplot);
      index++;
    }

    delete [] frequency[nb_state];
    delete [] cumul[nb_state];
  }

  if ((ident == GAMMA) || (ident == ZERO_INFLATED_GAMMA)) {
    for (i = 0;i < nb_state;i++) {
      delete gamma_dist[i];
    }
    delete [] gamma_dist;
  }

  else if (ident == GAUSSIAN) {
    for (i = 0;i < nb_state;i++) {
      delete gaussian_dist[i];
    }
    delete [] gaussian_dist;
  }

  for (i = 0;i < nb_state;i++) {
    delete [] frequency[i];
  }
  delete [] frequency;

  for (i = 0;i < nb_state;i++) {
    delete [] cumul[i];
  }
  delete [] cumul;

  delete [] dist_max;
}


/*--------------------------------------------------------------*
 *
 *  Calcul du nombre de parametres independants d'un processus
 *  d'observation continu parametrique.
 *
 *--------------------------------------------------------------*/

int ContinuousParametricProcess::nb_parameter_computation() const

{
  register int i;
  int nb_parameter = 0;


  for (i = 0;i < nb_state;i++) {
    nb_parameter += observation[i]->nb_parameter_computation();
  }

/*  if ((ident == GAUSSIAN) || (ident == VON_MISES)) {
    for (i = 1;i < nb_state;i++) {
      if (observation[i]->dispersion != observation[0]->dispersion) {
        break;
      }
    }

    if (i == nb_state) {
      nb_parameter -= (nb_state - 1);
    }
  } */

  if (tied_location) {
    nb_parameter -= (nb_state - 1);
  }
  if (tied_dispersion) {
    nb_parameter -= (nb_state - 1);
  }

  return nb_parameter;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la moyenne d'un melange de lois d'observation continues parametriques.
 *
 *  argument : loi des poids.
 *
 *--------------------------------------------------------------*/

double ContinuousParametricProcess::mean_computation(Distribution *pweight) const

{
  register int i;
  double mean;


  switch (ident) {

  case GAMMA : {
    mean = 0.;
    for (i = 0;i < nb_state;i++) {
      mean += pweight->mass[i] * observation[i]->shape * observation[i]->scale;
    }
    break;
  }

  case ZERO_INFLATED_GAMMA : {
    mean = 0.;
    for (i = 0;i < nb_state;i++) {
      mean += pweight->mass[i] * (1 - observation[i]->zero_probability) *
              observation[i]->shape * observation[i]->scale;
    }
    break;
  }

  case GAUSSIAN : {
    mean = 0.;
    for (i = 0;i < nb_state;i++) {
      mean += pweight->mass[i] * observation[i]->location;
    }
    break;
  }

  default : {
    mean = D_INF;
    break;
  }
  }

  return mean;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la variance d'un melange de lois d'observation continues parametriques.
 *
 *  arguments : loi des poids, moyenne.
 *
 *--------------------------------------------------------------*/

double ContinuousParametricProcess::variance_computation(Distribution *pweight , double mean) const

{
  register int i;
  double variance;


  if (mean == D_INF) {
    mean = mean_computation(pweight);
  }

  switch (ident) {

  case GAMMA : {
    variance = -mean * mean;
    for (i = 0;i < nb_state;i++) {
      variance += pweight->mass[i] * (observation[i]->shape * observation[i]->scale * observation[i]->scale) *
                  (1 + observation[i]->shape);
    }
    break;
  }

  case ZERO_INFLATED_GAMMA : {
    variance = -mean * mean;
    for (i = 0;i < nb_state;i++) {
      variance += pweight->mass[i] * (1 - observation[i]->zero_probability) *
                  (observation[i]->shape * observation[i]->scale * observation[i]->scale) *
                  (1 + observation[i]->shape);
    }
    break;
  }

  case GAUSSIAN : {
    variance = -mean * mean;
    for (i = 0;i < nb_state;i++) {
      variance += pweight->mass[i] * (observation[i]->dispersion * observation[i]->dispersion +
                                      observation[i]->location * observation[i]->location);
    }
    break;
  }

  default : {
    variance = D_DEFAULT;
    break;
  }
  }

  return variance;
}


/*--------------------------------------------------------------*
 *
 *  Choix de l'unite pour des lois d'observation de von Mises.
 *
 *  argument : unite (DEGREE/RADIAN).
 *
 *--------------------------------------------------------------*/

void ContinuousParametricProcess::select_unit(int iunit)

{
  if (ident == VON_MISES) {
    register int i;


    unit = iunit;
    for (i = 0;i < nb_state;i++) {
      observation[i]->unit = iunit;
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Initialisation des lois d'observation.
 *
 *  arguments : identificateur des lois, valeurs minimum et maximum,
 *              moyenne et variance empiriques de la variable.
 *
 *--------------------------------------------------------------*/

void ContinuousParametricProcess::init(int iident , double min_value , double max_value ,
                                       double mean , double variance)

{
  register int i;


  ident = iident;

  switch (ident) {

  case GAMMA : {
    for (i = 0;i < nb_state;i++) {
      observation[i] = new ContinuousParametric(ident , 1. , mean * (1. + MEAN_SHIFT_COEFF * (i - (double)(nb_state - 1) / 2) / (double)nb_state));
    }
    break;
  }

  case ZERO_INFLATED_GAMMA : {
    for (i = 0;i < nb_state;i++) {
      observation[i] = new ContinuousParametric(ident , 1. , mean * (1. + MEAN_SHIFT_COEFF * (i - (double)(nb_state - 1) / 2) / (double)nb_state) ,
                                                (double)(nb_state - i) / (double)(2 * nb_state));
    }
    break;
  }

  case GAUSSIAN : {
    for (i = 0;i < nb_state;i++) {
      observation[i] = new ContinuousParametric(ident , mean * (1. + MEAN_SHIFT_COEFF * (i - (double)(nb_state - 1) / 2) / (double)nb_state) ,
                                                sqrt(variance));
    }
    break;
  }

  case VON_MISES : {
    for (i = 0;i < nb_state;i++) {
      observation[i] = new ContinuousParametric(ident , mean * (1. + MEAN_SHIFT_COEFF * (i - (double)(nb_state - 1) / 2) / (double)nb_state) ,
                                                 1. / variance);
    }
    break;
  }
  }

# ifdef MESSAGE
  cout << "\nInitializations" << endl;
  if ((ident == GAMMA) || (ident == ZERO_INFLATED_GAMMA)) {
    cout << STAT_word[STATW_SCALE] << " | ";
    for (i = 0;i < nb_state;i++) {
      cout << STAT_label[STATL_STATE] << " " << i << ": " << observation[i]->scale << endl;
    }
  }
  else if ((ident == GAUSSIAN) || (ident == VON_MISES)) {
    for (i = 0;i < nb_state;i++) {
      cout << STAT_label[STATL_STATE] << " " << i << ": " << observation[i]->location << endl;
    }
  }
  cout << endl;
# endif

}


/*--------------------------------------------------------------*
 *
 *  Calcul d'intervalles sur des criteres de quantiles et de probabilites
 *  a posteriori pour des poids supposes egaux sur les lois d'observation.
 *
 *  argument : stream.
 *
 *--------------------------------------------------------------*/

ostream& ContinuousParametricProcess::interval_computation(ostream &os)

{
  register int i , j , k;
  int nb_step , posterior_mode;
  double step , value , min_value , max_value , cumul , norm , max_posterior ,
         posterior_value , **quantile_limit , **posterior , **posterior_limit;
  normal **dist;


  quantile_limit = new double*[2 * nb_state];
  for (i = 0;i < 2 * nb_state;i++) {
    quantile_limit[i] = new double[10];
  }

  posterior = new double*[nb_state];

  posterior_limit = new double*[2 * nb_state];
  for (i = 0;i < 2 * nb_state;i++) {
    posterior_limit[i] = new double[10];
  }

  // calcul des lois a posteriori pour des poids egaux des lois

  switch (ident) {

  case GAUSSIAN : {
    dist = new normal*[nb_state];

    for (i = 0;i < nb_state;i++) {
      dist[i] = new normal(observation[i]->location , observation[i]->dispersion);
      posterior[i] = new double[GAUSSIAN_NB_STEP];
    }

    // calcul des quantiles

    for (i = 0;i < nb_state;i++) {
      for (j = 0;j < 7;j++) {
        quantile_limit[2 * i][j] = quantile(*dist[i] , bilateral_tail[j] / 2.);
        quantile_limit[2 * i + 1][j] = quantile(complement(*dist[i] , bilateral_tail[j] / 2.));
      }
    }

    // calcul des probabilites a posteriori

    min_value = quantile(*dist[0] , GAUSSIAN_TAIL);
    for (i = 1;i < nb_state;i++) {
      value = quantile(*dist[i] , GAUSSIAN_TAIL);
      if (value < min_value) {
        min_value = value;
      }
    }

    max_value = quantile(complement(*dist[0] , GAUSSIAN_TAIL));
    for (i = 1;i < nb_state;i++) {
      value = quantile(complement(*dist[i] , GAUSSIAN_TAIL));
      if (value > max_value) {
        max_value = value;
      }
    }

    step = (max_value - min_value) / GAUSSIAN_NB_STEP + 1;

    value = min_value;
    for (i = 0;i < GAUSSIAN_NB_STEP;i++) {
      norm = 0.;
      for (j = 0;j < nb_state;j++) {
        posterior[j][i] = pdf(*dist[j] , value);
        norm += posterior[j][i];
      }

      for (j = 0;j < nb_state;j++) {
        posterior[j][i] /= norm;
      }
      value += step;
    }

    for (i = 0;i < nb_state;i++) {
      delete dist[i];
    }
    delete [] dist;

    nb_step = GAUSSIAN_NB_STEP;
    break;
  }

  case VON_MISES : {
    switch (unit) {
    case DEGREE :
      step = 360. / VON_MISES_NB_STEP;
      break;
    case RADIAN :
      step = 2 * M_PI / VON_MISES_NB_STEP;
      break;
    }

    // calcul des quantiles

    for (i = 0;i < nb_state;i++) {
      value = observation[i]->location + step / 2;
      cumul = 0.5;
      j = 0;

      do {
        switch (unit) {

        case DEGREE : {
          do {
            value -= step;
            cumul -= exp(observation[i]->dispersion * cos((value - observation[i]->location) * M_PI / 180)) * step /
                     (360 * cyl_bessel_i(0 , observation[i]->dispersion));
          }
          while (cumul > bilateral_tail[j] / 2);

          quantile_limit[2 * i][j] = value;
          quantile_limit[2 * i + 1][j] = 2 * observation[i]->location - value;

          if (quantile_limit[2 * i][j] < 0) {
            quantile_limit[2 * i][j] += 360;
          }
          if (quantile_limit[2 * i + 1][j] >= 360) {
            quantile_limit[2 * i + 1][j] -= 360;
          }
          break;
        }

        case RADIAN : {
          do {
            value -= step;
            cumul -= exp(observation[i]->dispersion * cos(value - observation[i]->location)) * step /
                     (2 * M_PI * cyl_bessel_i(0 , observation[i]->dispersion));
          }
          while (cumul > bilateral_tail[j] / 2);

          quantile_limit[2 * i][j] = value;
          quantile_limit[2 * i + 1][j] = 2 * observation[i]->location - value;
 
          if (quantile_limit[2 * i][j] < 0) {
            quantile_limit[2 * i][j] += 2 * M_PI;
          }
          if (quantile_limit[2 * i + 1][j] >= 2 * M_PI) {
            quantile_limit[2 * i + 1][j] -= 2 * M_PI;
          }
          break;
        }
        }

        j++;
      }
      while (j < 7);
    }

    // calcul des probabilites a posteriori

    for (i = 0;i < nb_state;i++) {
      posterior[i] = new double[VON_MISES_NB_STEP];
    }

    value = 0.;

    switch (unit) {

    case DEGREE : {
      for (i = 0;i < VON_MISES_NB_STEP;i++) {
        norm = 0.;
        for (j = 0;j < nb_state;j++) {
          posterior[j][i] = exp(observation[j]->dispersion * cos((value - observation[j]->location) * M_PI / 180)) /
                            (360 * cyl_bessel_i(0 , observation[j]->dispersion));
          norm += posterior[j][i];
        }

        for (j = 0;j < nb_state;j++) {
          posterior[j][i] /= norm;
        }
        value += step;
      }
      break;
    }

    case RADIAN : {
      for (i = 0;i < VON_MISES_NB_STEP;i++) {
        norm = 0.;
        for (j = 0;j < nb_state;j++) {
          posterior[j][i] = exp(observation[j]->dispersion * cos(value - observation[j]->location)) /
                            (2 * M_PI * cyl_bessel_i(0 , observation[j]->dispersion));
          norm += posterior[j][i];
        }

        for (j = 0;j < nb_state;j++) {
          posterior[j][i] /= norm;
        }
        value += step;
      }
      break;
    }
    }

    nb_step = VON_MISES_NB_STEP;
    break;
  }
  }

  // calcul des seuils sur les probabilites a posteriori

# ifdef MESSAGE
  cout << "\nPosterior mode:";
# endif

  for (i = 0;i < nb_state;i++) {
    max_posterior = 0.;
    posterior_mode = 0;

    switch (ident) {
    case GAUSSIAN :
      value = min_value;
      break;
    case VON_MISES :
      value = 0.;
      break;
    }
    posterior_value = value;

    for (j = 1;j < nb_step;j++) {
      value += step;
      if  (posterior[i][j] > max_posterior) {
        max_posterior = posterior[i][j];
        posterior_mode = j;
        posterior_value = value;
      }
    }

#   ifdef MESSAGE
    cout << " " << posterior_value;
#   endif

    j = posterior_mode;
    value = posterior_value;
    k = 0;
    do {
      do {
        if ((ident == VON_MISES) && (j == 0)) {
          j = nb_step - 1;
          value = j * step;
        }
        else {
          j--;
          value -= step;
        }
      }
      while ((posterior[i][j] > posterior_threshold[k]));

      posterior_limit[2 * i][k++] = value;
    }
    while (k < 7);

    j = posterior_mode;
    value = posterior_value;
    k = 0;
    do {
      do {
        if ((ident == VON_MISES) && (j == nb_step - 1)) {
          j = 0;
          value = 0.;
        }
        else {
          j++;
          value += step;
        }
      }
      while ((posterior[i][j] > posterior_threshold[k]));

      posterior_limit[2 * i + 1][k++] = value;
    }
    while (k < 7);
  }

# ifdef MESSAGE
  cout << endl;
# endif

  // sortie des seuils

  for (i = 0;i < nb_state;i++) {
    os << "\n" << STAT_word[STATW_STATE] << " " << i << " "
       << STAT_word[STATW_OBSERVATION_DISTRIBUTION] << endl;
    observation[i]->ascii_parameter_print(os);
    observation[i]->ascii_characteristic_print(os);

    os << "\nQuantiles" << endl;
    for (j = 0;j < 7;j++) {
      os << bilateral_tail[j] << ": " << quantile_limit[2 * i][j] << ", "
         << quantile_limit[2 * i + 1][j] << " - ";

      if ((ident == GAUSSIAN) || (quantile_limit[2 * i + 1][j] > quantile_limit[2 * i][j])) {
        os << quantile_limit[2 * i + 1][j] - quantile_limit[2 * i][j] << endl;
      }

      else {
        switch (unit) {
        case DEGREE :
          os << quantile_limit[2 * i + 1][j] - (quantile_limit[2 * i][j] - 360) << endl;
          break;
        case RADIAN :
          os << quantile_limit[2 * i + 1][j] - (quantile_limit[2 * i][j] - 2 * M_PI) << endl;
          break;
        }
      }
    }

    os << "\nPosterior probability limits" << endl;
    for (j = 0;j < 7;j++) {
      os << posterior_threshold[j] << ": " << posterior_limit[2 * i][j] << ", "
         << posterior_limit[2 * i + 1][j] << " - ";

      if ((ident == GAUSSIAN) || (posterior_limit[2 * i + 1][j] > posterior_limit[2 * i][j])) {
        os << posterior_limit[2 * i + 1][j] - posterior_limit[2 * i][j] << endl;
      }

      else {
        switch (unit) {
        case DEGREE :
          os << posterior_limit[2 * i + 1][j] - (posterior_limit[2 * i][j] - 360) << endl;
          break;
        case RADIAN :
          os << posterior_limit[2 * i + 1][j] - (posterior_limit[2 * i][j] - 2 * M_PI) << endl;
          break;
        }
      }
    }

    os << endl;
  }

  for (i = 0;i < 2 * nb_state;i++) {
    delete [] quantile_limit[i];
  }
  delete [] quantile_limit;

  for (i = 0;i < nb_state;i++) {
    delete [] posterior[i];
  }
  delete [] posterior;

  for (i = 0;i < 2 * nb_state;i++) {
    delete [] posterior_limit[i];
  }
  delete [] posterior_limit;

  return os;
}


};  // namespace stat_tool
