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
 *       $Id: continuous_parametric.cpp 8175 2010-02-18 10:30:35Z guedon $
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

#include <boost/math/distributions/gamma.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/special_functions/bessel.hpp>

#include "tool/rw_tokenizer.h"
#include "tool/rw_cstring.h"
#include "tool/rw_locale.h"
#include "tool/config.h"

#include "stat_tools.h"
#include "stat_label.h"

using namespace std;
using namespace boost::math;


namespace stat_tool {



/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe ContinuousParametric.
 *
 *  arguments : identificateur, parametres de localisation et de dispersion (Gauss ou von Mises)/
 *              de forme et d'echelle (gamma), probabilite pour 0 (zero-inflated gamma),
 *              intercept, pente et parametre de dispersion (modele lineaire)
 *              unite (loi de von Mises).
 *
 *--------------------------------------------------------------*/

ContinuousParametric::ContinuousParametric(int iident , double ilocation , double idispersion ,
                                           double izero_probability , int iunit)

{
  ident = iident;

  location = ilocation;
  dispersion = idispersion;
  zero_probability = izero_probability;

  min_value = D_DEFAULT;
  max_value = D_DEFAULT;

  unit = iunit;
  slope_standard_deviation = D_DEFAULT;
  sample_size = 0;
  correlation = D_DEFAULT;
  cumul = NULL;
}


/*--------------------------------------------------------------*
 *
 *  Copie d'un objet ContinuousParametric.
 *
 *  argument : reference sur un objet ContinuousParametric.
 *
 *--------------------------------------------------------------*/

void ContinuousParametric::copy(const ContinuousParametric &dist)

{
  ident = dist.ident;

  location = dist.location;
  dispersion = dist.dispersion;
  zero_probability = dist.zero_probability;

  min_value = dist.min_value;
  max_value = dist.max_value;

  unit = dist.unit;
  slope_standard_deviation = dist.slope_standard_deviation;
  sample_size = dist.sample_size;
  correlation = dist.correlation;

  if (dist.cumul) {
    register int i;

    cumul = new double[VON_MISES_NB_STEP];

    for (i = 0;i < VON_MISES_NB_STEP;i++) {
      cumul[i] = dist.cumul[i];
    }
  }

  else {
    cumul = NULL;
  }
}


/*--------------------------------------------------------------*
 *
 *  Destructeur de la classe ContinuousParametric.
 *
 *--------------------------------------------------------------*/

ContinuousParametric::~ContinuousParametric()

{
  delete [] cumul;
}


/*--------------------------------------------------------------*
 *
 *  Operateur d'assignement de la classe ContinuousParametric.
 *
 *  argument : reference sur un objet ContinuousParametric.
 *
 *--------------------------------------------------------------*/

ContinuousParametric& ContinuousParametric::operator=(const ContinuousParametric &dist)

{
  if (&dist != this) {
    delete [] cumul;
    copy(dist);
  }

  return *this;
}


/*--------------------------------------------------------------*
 *
 *  Analyse du format d'un objet ContinuousParametric.
 *
 *  arguments : reference sur un objet StatError, stream,
 *              reference sur l'indice de la ligne lue, identificateur
 *              de la derniere loi dans la liste.
 *
 *--------------------------------------------------------------*/

ContinuousParametric* continuous_parametric_parsing(StatError &error , ifstream &in_file ,
                                                    int &line , int last_ident)

{
  RWLocaleSnapshot locale("en");
  RWCString buffer , token;
  size_t position;
  bool status = true , lstatus;
  register int i , j;
  int ident = I_DEFAULT;
  union {
    double shape;
    double location;
    double intercept;
  };
  union {
    double scale;
    double dispersion;
  };
  union {
    double zero_probability;
    double slope;
  };
  ContinuousParametric *dist;


  dist = NULL;
  location = D_INF;
  dispersion = D_DEFAULT;
  slope = D_INF;

  while (buffer.readLine(in_file , false)) {
    line++;

#   ifdef DEBUG
    cout << line << " " << buffer << endl;
#   endif

    position = buffer.first('#');
    if (position != RW_NPOS) {
      buffer.remove(position);
    }
    i = 0;

    RWCTokenizer next(buffer);

    while (!((token = next()).isNull())) {

      // test nom de la loi

      if (i == 0) {
        for (j = GAMMA;j <= last_ident;j++) {
          if (token == STAT_continuous_distribution_word[j]) {
            ident = j;
            break;
          }
        }

        if (j == last_ident + 1) {
          status = false;
          error.update(STAT_parsing[STATP_DISTRIBUTION_NAME] , line , i + 1);
        }
      }

      // test nom du parametre

      else {
        switch ((i - 1) % 3) {

        case 0 : {
          switch ((i - 1) / 3) {

          // 1er parametre : parametre de forme (gamma), probabilite pour 0 (zero-inflated gamma),
          // moyenne (Gauss), direction moyenne (von Mises), intercept (modele lineaire)

          case 0 : {
            if ((ident == GAMMA) && (token != STAT_word[STATW_SHAPE])) {
              status = false;
              error.correction_update(STAT_parsing[STATP_PARAMETER_NAME] , STAT_word[STATW_SHAPE] , line , i + 1);
            }

            if ((ident == ZERO_INFLATED_GAMMA) && (token != STAT_word[STATW_ZERO_PROBABILITY])) {
              status = false;
              error.correction_update(STAT_parsing[STATP_PARAMETER_NAME] , STAT_word[STATW_ZERO_PROBABILITY] , line , i + 1);
            }

            if ((ident == GAUSSIAN) && (token != STAT_word[STATW_MEAN])) {
              status = false;
              error.correction_update(STAT_parsing[STATP_PARAMETER_NAME] , STAT_word[STATW_MEAN] , line , i + 1);
            }

            if ((ident == VON_MISES) && (token != STAT_word[STATW_MEAN_DIRECTION])) {
              status = false;
              error.correction_update(STAT_parsing[STATP_PARAMETER_NAME] , STAT_word[STATW_MEAN_DIRECTION] , line , i + 1);
            }

            if ((ident == LINEAR_MODEL) && (token != STAT_word[STATW_INTERCEPT])) {
              status = false;
              error.correction_update(STAT_parsing[STATP_PARAMETER_NAME] , STAT_word[STATW_INTERCEPT] , line , i + 1);
            }
            break;
          }

          // 2eme parametre : parametre d'echelle (gamma), parametre de forme (zero-inflated gamma),
          // ecart-type (Gauss), concentration (von Mises), pente (modele lineaire)

          case 1 : {
            if ((ident == GAMMA) && (token != STAT_word[STATW_SCALE])) {
              status = false;
              error.correction_update(STAT_parsing[STATP_PARAMETER_NAME] , STAT_word[STATW_SCALE] , line , i + 1);
            }

            if ((ident == ZERO_INFLATED_GAMMA) && (token != STAT_word[STATW_SHAPE])) {
              status = false;
              error.correction_update(STAT_parsing[STATP_PARAMETER_NAME] , STAT_word[STATW_SHAPE] , line , i + 1);
            }

            if ((ident == GAUSSIAN) && (token != STAT_word[STATW_STANDARD_DEVIATION])) {
              status = false;
              error.correction_update(STAT_parsing[STATP_PARAMETER_NAME] , STAT_word[STATW_STANDARD_DEVIATION] , line , i + 1);
            }

            if ((ident == VON_MISES) && (token != STAT_word[STATW_CONCENTRATION])) {
              status = false;
              error.correction_update(STAT_parsing[STATP_PARAMETER_NAME] , STAT_word[STATW_CONCENTRATION] , line , i + 1);
            }

            if ((ident == LINEAR_MODEL) && (token != STAT_word[STATW_SLOPE])) {
              status = false;
              error.correction_update(STAT_parsing[STATP_PARAMETER_NAME] , STAT_word[STATW_SLOPE] , line , i + 1);
            }
            break;
          }

          // 3eme parametre : parametre d'echelle (zeo_inflated gamma), ecart-type (modele lineaire)

          case 2 : {
            if ((ident == ZERO_INFLATED_GAMMA) && (token != STAT_word[STATW_SCALE])) {
              status = false;
              error.correction_update(STAT_parsing[STATP_PARAMETER_NAME] , STAT_word[STATW_SCALE] , line , i + 1);
            }

            if ((ident == LINEAR_MODEL) && (token != STAT_word[STATW_STANDARD_DEVIATION])) {
              status = false;
              error.correction_update(STAT_parsing[STATP_PARAMETER_NAME] , STAT_word[STATW_STANDARD_DEVIATION] , line , i + 1);
            }
            break;
          }
          }

          break;
        }

        // test separateur

        case 1 : {
          if (token != ":") {
            status = false;
            error.update(STAT_parsing[STATP_SEPARATOR] , line , i + 1);
          }
          break;
        }

        // test valeur du parametre

        case 2 : {
          switch ((i - 1) / 3) {

          // 1er parametre : parametre de forme (gamma), probabilite pour 0 (zero-inflated gamma),
          // moyenne (Gauss), direction moyenne (von Mises)

          case 0 : {
            if (ident == GAMMA) {
              lstatus = locale.stringToNum(token , &shape);
              if ((lstatus) && (shape < 0.)) {
                lstatus = false;
              }
            }

            else if (ident == ZERO_INFLATED_GAMMA) {
              lstatus = locale.stringToNum(token , &zero_probability);
              if ((lstatus) && ((zero_probability < 0.) || (zero_probability > 1.))) {
                lstatus = false;
              }
            }

            else if (ident == LINEAR_MODEL) {
              lstatus = locale.stringToNum(token , &intercept);
            }

            else {
              lstatus = locale.stringToNum(token , &location);
            }
            break;
          }

          // 2eme parametre : parametre d'echelle (gamma), parametre de forme (zero-inflated gamma),
          // ecart-type (Gauss), concentration (von Mises), pente (modele lineaire)

          case 1 : {
            if (ident == GAMMA) {
              lstatus = locale.stringToNum(token , &scale);
              if ((lstatus) && (scale <= 0.)) {
                lstatus = false;
              }
            }

            else if (ident == ZERO_INFLATED_GAMMA) {
              lstatus = locale.stringToNum(token , &shape);
              if ((lstatus) && (shape <= 0.)) {
                lstatus = false;
              }
            }

            else if (ident == LINEAR_MODEL) {
              lstatus = locale.stringToNum(token , &slope);
            }

            else {
              lstatus = locale.stringToNum(token , &dispersion);
              if ((lstatus) && (dispersion <= 0.)) {
                lstatus = false;
              }
            }
            break;
          }

          // 3eme parametre : parametre d'echelle (zero-inflated gamma), ecart-type (modele lineaire)

          case 2 : {
            if (ident == ZERO_INFLATED_GAMMA) {
              lstatus = locale.stringToNum(token , &scale);
              if ((lstatus) && (scale <= 0.)) {
                lstatus = false;
              }
            }

            else if (ident == LINEAR_MODEL) {
              lstatus = locale.stringToNum(token , &dispersion);
              if ((lstatus) && (dispersion <= 0.)) {
                lstatus = false;
              }
            }
            break;
          }
          }

          if (!lstatus) {
            status = false;
            error.update(STAT_parsing[STATP_PARAMETER_VALUE] , line , i + 1);
          }
          break;
        }
        }
      }

      i++;
    }

    if (i > 0) {
      if ((((ident == GAMMA) || (ident == GAUSSIAN) || (ident == VON_MISES)) && (i != 7)) ||
          (((ident == ZERO_INFLATED_GAMMA) || (ident == LINEAR_MODEL)) && (i != 10))) {
        status = false;
        error.update(STAT_parsing[STATP_FORMAT] , line);
      }

      break;
    }
  }

  if (ident == I_DEFAULT) {
    status = false;
    error.update(STAT_parsing[STATP_FORMAT] , line);
  }

  if (status) {
    dist = new ContinuousParametric(ident , location , dispersion , zero_probability ,
                                    (location < 2 * M_PI ? RADIAN : DEGREE));
  }

  return dist;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture des parametres d'une loi continue.
 *
 *  arguments : stream, flag commentaire.
 *
 *--------------------------------------------------------------*/

ostream& ContinuousParametric::ascii_parameter_print(ostream &os , bool file_flag) const

{
  os << STAT_continuous_distribution_word[ident] << "   ";

  switch (ident) {

  case GAMMA : {
    os << STAT_word[STATW_SHAPE] << " : " << shape;
    if (shape > 0.) {
      os << "   " << STAT_word[STATW_SCALE] << " : " << scale;
    }
    break;
  }

  case ZERO_INFLATED_GAMMA : {
    os << STAT_word[STATW_ZERO_PROBABILITY] << " : " << zero_probability;
    if (zero_probability < 1.) {
      os << "   " << STAT_word[STATW_SHAPE] << " : " << shape
         << "   " << STAT_word[STATW_SCALE] << " : " << scale;
    }
    break;
  }

  case GAUSSIAN : {
    os << STAT_word[STATW_MEAN] << " : " << location << "   "
       << STAT_word[STATW_STANDARD_DEVIATION] << " : " << dispersion;
    break;
  }

  case VON_MISES : {
    os << STAT_word[STATW_MEAN_DIRECTION] << " : " << location << "   "
       << STAT_word[STATW_CONCENTRATION] << " : " << dispersion;
    break;
  }

  case LINEAR_MODEL : {
    Test test(STUDENT , false , sample_size , I_DEFAULT , D_DEFAULT);

    os << STAT_word[STATW_INTERCEPT] << " : " << intercept << "   "
       << STAT_word[STATW_SLOPE] << " : " << slope;

    test.critical_probability = ref_critical_probability[0];
    test.t_value_computation();

    if ((!file_flag) && (slope_standard_deviation > 0.)) {
      os << "   (" << slope - test.value * slope_standard_deviation << ", "
         << slope + test.value * slope_standard_deviation << ")";
    }

    os << "   " << STAT_word[STATW_STANDARD_DEVIATION] << " : " << dispersion;

    if ((file_flag) && (slope_standard_deviation > 0.)) {
      os << "   # " << slope - test.value * slope_standard_deviation << ", "
         << slope + test.value * slope_standard_deviation;
    }

    os << "\n";
    if (file_flag) {
      os << "# ";
    }
    os << STAT_label[STATL_CORRELATION_COEFF] << ": " << correlation << "   "
       << STAT_label[STATL_LIMIT_CORRELATION_COEFF] << ": -/+"
       << test.value / sqrt(test.value * test.value + sample_size) << "   "
       << STAT_label[STATL_CRITICAL_PROBABILITY] << ": " << test.critical_probability;
    break;
  }
  }
  os << endl;

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture des caracteristiques d'une loi continue.
 *
 *  arguments : stream, flag commentaire.
 *
 *--------------------------------------------------------------*/

ostream& ContinuousParametric::ascii_characteristic_print(ostream &os , bool file_flag) const

{
  if (file_flag) {
    os << "# ";
  }

  switch (ident) {

  case GAMMA : {
    double variance = shape * scale * scale;

    os << STAT_label[STATL_MEAN] << ": " << shape * scale << "   "
       << STAT_label[STATL_VARIANCE] << ": " << variance << "   "
       << STAT_label[STATL_STANDARD_DEVIATION] << ": " << sqrt(variance) << endl;

    if (shape > 0.) {
      if (file_flag) {
        os << "# ";
      }
      os << STAT_label[STATL_SKEWNESS_COEFF] << ": " << 2 / sqrt(shape) << "   "
         << STAT_label[STATL_KURTOSIS_COEFF] << ": " << 6 / shape << endl;
    }
    break;
  }

  case ZERO_INFLATED_GAMMA : {
    double mean , variance;

    if (zero_probability == 1.) {
      mean = 0.;
      variance = 0.;
    }
    else {
      mean = (1 - zero_probability) * shape * scale;
      variance = (1 - zero_probability) * shape * scale * scale * (1 + shape) - mean * mean;
    }

    os << STAT_label[STATL_MEAN] << ": " <<  mean << "   "
       << STAT_label[STATL_VARIANCE] << ": " << variance << "   "
       << STAT_label[STATL_STANDARD_DEVIATION] << ": " << sqrt(variance) << endl;
    break;
  }

  case GAUSSIAN : {
    os << STAT_label[STATL_VARIANCE] << ": " << dispersion * dispersion << endl;
    break;
  }

  case VON_MISES : {
    double mean_resultant_length , standard_deviation;


    mean_resultant_length = cyl_bessel_i(1 , dispersion) / cyl_bessel_i(0 , dispersion);

    os << STAT_label[STATL_MEAN_RESULTANT_LENGTH] << ": " << mean_resultant_length;

    if (mean_resultant_length > 0.) {
      standard_deviation = sqrt(-2 * log(mean_resultant_length));
      if (unit == DEGREE) {
        standard_deviation *= (180 / M_PI);
      }

      os << "   " << STAT_label[STATL_CIRCULAR_STANDARD_DEVIATION] << ": "
         << standard_deviation;
    }
    os << endl;

#   ifdef MESSAGE
    if (dispersion > 0.) {
      standard_deviation = 1. / sqrt(fabs(dispersion));
      if (unit == DEGREE) {
        standard_deviation *= (180 / M_PI);
      }

      os << STAT_label[STATL_STANDARD_DEVIATION] << ": " <<  standard_deviation << endl;
    }
#   endif

    break;
  }
  }

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'une loi continue.
 *
 *  arguments : stream, flag commentaire, flags sur l'ecriture de la fonction
 *              de repartition, pointeurs sur un objet Histogram et
 *              sur un objet FrequencyDistribution.
 *
 *--------------------------------------------------------------*/

ostream& ContinuousParametric::ascii_print(ostream &os , bool file_flag ,
                                           bool cumul_flag , const Histogram *histo1 ,
                                           const FrequencyDistribution *histo2)

{
  if ((histo1) || (histo2)) {
    register int i , j;
    int nb_step , width[5];
    long old_adjust;
    double step , scale , mass , value , *frequency , *dist_cumul , *histo_cumul;


    old_adjust = os.setf(ios::right , ios::adjustfield);

    if (histo1) {
      step = histo1->step;
      scale = histo1->nb_element;
    }
    else {
      step = histo2->min_interval_computation();
      scale = histo2->nb_element;
    }

    if (ident == GAMMA) {
      min_value = 0.;
      if (shape == 0.) {
        max_value = 0.;
      }
      else {
        max_value = shape * scale;
      }
    }

    else if (ident == ZERO_INFLATED_GAMMA) {
      min_value = 0.;
      if (zero_probability == 1.) {
        max_value = 0.;
      }
      else {
        max_value = shape * scale;
      }
    }

    else if ((ident == GAUSSIAN) || (ident == VON_MISES)) {
      min_value = location;

      if (histo1) {
        switch (histo1->type) {

        case INT_VALUE : {
          if (histo1->min_value < min_value) {
            min_value = histo1->min_value;
          }
          break;
        }

        case REAL_VALUE : {
          if (histo1->min_value - histo1->step / 2 < min_value) {
            min_value = histo1->min_value - histo1->step / 2;
          }
          break;
        }
        }
      }

      if ((histo2) && (histo2->offset < min_value)) {
        min_value = histo2->offset;
      }

      max_value = location;
    }

    if (histo1) {
      switch (histo1->type) {

      case INT_VALUE : {
        if (histo1->max_value > max_value) {
          max_value = histo1->max_value;
        }
        break;
      }

      case REAL_VALUE : {
        if (histo1->max_value + histo1->step / 2 > max_value) {
          max_value = histo1->max_value + histo1->step / 2;
        }
        break;
      }
      }
    }

    if ((histo2) && (histo2->nb_value - 1 > max_value)) {
      max_value = histo2->nb_value - 1;
    }

    if (cumul_flag) {
      if (histo1) {
        histo_cumul = histo1->cumul_computation();
      }
      else {
        histo_cumul = histo2->cumul_computation();
      }
    }

    // calcul d'une loi continue discretisee

    switch (ident) {

    case GAMMA : {
      if (shape == 0.) {
        dist_cumul = new double[1];
        frequency = new double[1];

        dist_cumul[0] = 1.;
        frequency[0] = scale;
      }

      else {
        gamma_distribution<double> dist(shape , scale);

        value = quantile(complement(dist , GAMMA_TAIL));
        while (max_value < value) {
          max_value += step;
        }

        nb_step = (int)((max_value - min_value) / step) + 1;
        dist_cumul = new double[nb_step];
        frequency = new double[nb_step];

//        value = step;
//        dist_cumul[0] = cdf(dist , value);
        value = step / 2;
        dist_cumul[0] = cdf(dist , value + step / 2);
        frequency[0] = dist_cumul[0] * scale;

        for (i = 1;i < nb_step;i++) {
          value += step;
//          dist_cumul[i] = cdf(dist , value);
          dist_cumul[i] = cdf(dist , value + step / 2);
          frequency[i] = (dist_cumul[i] - dist_cumul[i - 1]) * scale;
        }
      }
      break;
    }

    case ZERO_INFLATED_GAMMA : {
      if (zero_probability == 1.) {
        dist_cumul = new double[1];
        frequency = new double[1];

        dist_cumul[0] = 1.;
        frequency[0] = scale;
      }

      else {
        gamma_distribution<double> dist(shape , scale);

        value = quantile(complement(dist , GAMMA_TAIL));
        while (max_value < value) {
          max_value += step;
        }

        nb_step = (int)((max_value - min_value) / step) + 1;
        dist_cumul = new double[nb_step];
        frequency = new double[nb_step];

        value = 0.;
        dist_cumul[0] = zero_probability;
        frequency[0] = zero_probability * scale;

        for (i = 1;i < nb_step;i++) {
          value += step;
          dist_cumul[i] = zero_probability + (1. - zero_probability) * cdf(dist , value);
          frequency[i] = (dist_cumul[i] - dist_cumul[i - 1]) * scale;
        }
      }
      break;
    }

    case GAUSSIAN : {
      normal dist(location , dispersion);

      value = quantile(dist , GAUSSIAN_TAIL);
      while (min_value > value) {
        min_value -= step;
      }

//      value = quantile(dist , 1. - GAUSSIAN_TAIL);
      value = quantile(complement(dist , GAUSSIAN_TAIL));
      while (max_value < value) {
        max_value += step;
      }

      nb_step = (int)((max_value - min_value) / step) + 1;
      dist_cumul = new double[nb_step];
      frequency = new double[nb_step];

      value = min_value;
      dist_cumul[0] = cdf(dist , value + step / 2);
      frequency[0] = (dist_cumul[0] - cdf(dist , value - step / 2)) * scale;

      for (i = 1;i < nb_step;i++) {
        value += step;
        dist_cumul[i] = cdf(dist , value + step / 2);
        frequency[i] = (dist_cumul[i] - dist_cumul[i - 1]) * scale;
      }
      break;
    }

    case VON_MISES : {
      while (min_value >= step / 2) {
        min_value -= step;
      }

      value = (unit == DEGREE ? 360 : 2 * M_PI) - 3 * step / 2;
      while (max_value <= value) {
        max_value += step;
      }

      nb_step = (int)((max_value - min_value) / step) + 1;
      frequency = new double[nb_step];
      dist_cumul = new double[nb_step];

      value = min_value;
      mass = von_mises_mass_computation(value - step / 2 , value + step / 2);
      frequency[0] = mass * scale;
      dist_cumul[0] = mass;

      for (i = 1;i < nb_step;i++) {
        value += step;
        mass = von_mises_mass_computation(value - step / 2 , value + step / 2);
        frequency[i] = mass * scale;
        dist_cumul[i] = dist_cumul[i - 1] + mass;
      }
      break;
    }
    }

    // calcul des largeurs des colonnes

    width[0] = column_width(min_value , max_value);

    if (histo1) {
      width[1] = column_width(histo1->max) + ASCII_SPACE;
    }
    else {
      width[1] = column_width(histo2->max) + ASCII_SPACE;
    }

    width[2] = column_width(nb_step , frequency) + ASCII_SPACE;

    if (cumul_flag) {
      if (histo1) {
        width[3] = column_width(histo1->nb_category , histo_cumul) + ASCII_SPACE;
      }
      else {
        width[3] = column_width(histo2->nb_value - histo2->offset ,
                                histo_cumul + histo2->offset) + ASCII_SPACE;
      }

      width[4] = column_width(nb_step , dist_cumul) + ASCII_SPACE;
    }

    value = min_value;
    if (histo1) {
      i = -1;
    }
    else {
      i = histo2->offset - step;
    }

    for (j = 0;j < nb_step;j++) {
      if (file_flag) {
        os << "# ";
      }
      os << setw(width[0]) << value;

      if (histo1) {
        if ((value >= histo1->min_value) && (value <= histo1->max_value)) {
          os << setw(width[1]) << histo1->frequency[++i];
        }
        else {
          os << setw(width[1]) << " ";
        }
      }

      else {
        if ((value >= histo2->offset) && (value < histo2->nb_value)) {
          i += step;
          os << setw(width[1]) << histo2->frequency[i];
        }
        else {
          os << setw(width[1]) << " ";
        }
      }

      os << setw(width[2]) << frequency[j];

      if (cumul_flag) {
        if (((histo1) && (value >= histo1->min_value) && (value <= histo1->max_value)) ||
            ((histo2) && (value >= histo2->offset) && (value < histo2->nb_value))) {
          os << setw(width[3]) << histo_cumul[i];
        }
        else {
          os << setw(width[3]) << " ";
        }

        os << setw(width[4]) << dist_cumul[j];
      }

      os << endl;
      value += step;
    }

    delete [] frequency;
    delete [] dist_cumul;

    if (cumul_flag) {
      delete [] histo_cumul;
    }

    os.setf((FMTFLAGS)old_adjust , ios::adjustfield);
  }

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture des parametres d'une loi continue au format tableur.
 *
 *  argument : stream.
 *
 *--------------------------------------------------------------*/

ostream& ContinuousParametric::spreadsheet_parameter_print(ostream &os) const

{
  os << "\n" << STAT_continuous_distribution_word[ident] << "\t";

  switch (ident) {

  case GAMMA : {
    os << STAT_word[STATW_SHAPE] << "\t" << shape;
    if (shape > 0.) {
      os << "\t" << STAT_word[STATW_SCALE] << "\t" << scale;
    }
    break;
  }

  case ZERO_INFLATED_GAMMA : {
    os << STAT_word[STATW_ZERO_PROBABILITY] << "\t" << zero_probability;
    if (zero_probability < 1.) {
      os << "\t" << STAT_word[STATW_SHAPE] << "\t" << shape
         << "\t" << STAT_word[STATW_SCALE] << "\t" << scale;
    }
    break;
  }

  case GAUSSIAN : {
    os << STAT_word[STATW_MEAN] << "\t" << location << "\t"
       << STAT_word[STATW_STANDARD_DEVIATION] << "\t" << dispersion;
    break;
  }

  case VON_MISES : {
    os << STAT_word[STATW_MEAN_DIRECTION] << "\t" << location << "\t"
       << STAT_word[STATW_CONCENTRATION] << "\t" << dispersion;
    break;
  }

  case LINEAR_MODEL : {
    Test test(STUDENT , false , sample_size , I_DEFAULT , D_DEFAULT);

    os << STAT_word[STATW_INTERCEPT] << "\t" << intercept << "\t"
       << STAT_word[STATW_SLOPE] << "\t" << slope;

    test.critical_probability = ref_critical_probability[0];
    test.t_value_computation();

    if (slope_standard_deviation > 0.) {
      os << "\t" << slope - test.value * slope_standard_deviation << "\t"
         << slope + test.value * slope_standard_deviation;
    }

    os << "\t" << STAT_word[STATW_STANDARD_DEVIATION] << "\t" << dispersion
       << "\n" << STAT_label[STATL_CORRELATION_COEFF] << "\t" << correlation
       << "\t" << STAT_label[STATL_LIMIT_CORRELATION_COEFF] << "\t" << "-/+"
       << test.value / sqrt(test.value * test.value + sample_size)
       << "\t" << STAT_label[STATL_CRITICAL_PROBABILITY] << "\t" << test.critical_probability;
    break;
  }
  }
  os << endl;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture des caracteristiques d'une loi continue au format tableur.
 *
 *  argument : stream.
 *
 *--------------------------------------------------------------*/

ostream& ContinuousParametric::spreadsheet_characteristic_print(ostream &os) const

{
  switch (ident) {

  case GAMMA : {
    double variance = shape * scale * scale;

    os << STAT_label[STATL_MEAN] << "\t" << shape * scale << "\t"
       << STAT_label[STATL_VARIANCE] << "\t" << variance << "\t"
       << STAT_label[STATL_STANDARD_DEVIATION] << "\t" << sqrt(variance) << endl;

    if (shape > 0.) {
      os << STAT_label[STATL_SKEWNESS_COEFF] << "\t" << 2 / sqrt(shape) << "\t"
         << STAT_label[STATL_KURTOSIS_COEFF] << "\t" << 6 / shape << endl;
    }
    break;
  }

  case ZERO_INFLATED_GAMMA : {
    double mean , variance;

    if (zero_probability == 1.) {
      mean = 0.;
      variance = 0.;
    }
    else {
      mean = (1 - zero_probability) * shape * scale;
      variance = (1 - zero_probability) * shape * scale * scale * (1 + shape) - mean * mean;
    }

    os << STAT_label[STATL_MEAN] << "\t" <<  mean << "\t"
       << STAT_label[STATL_VARIANCE] << "\t" << variance << "\t"
       << STAT_label[STATL_STANDARD_DEVIATION] << "\t" << sqrt(variance) << endl;
    break;
  }

  case GAUSSIAN : {
    os << STAT_label[STATL_VARIANCE] << "\t" << dispersion * dispersion << endl;
    break;
  }

  case VON_MISES : {
    double mean_resultant_length , standard_deviation;


    mean_resultant_length = cyl_bessel_i(1 , dispersion) / cyl_bessel_i(0 , dispersion);

    os << STAT_label[STATL_MEAN_RESULTANT_LENGTH] << "\t" << mean_resultant_length;

    if (mean_resultant_length > 0.) {
      standard_deviation = sqrt(-2 * log(mean_resultant_length));
      if (unit == DEGREE) {
        standard_deviation *= (180 / M_PI);
      }

      os << "\t" << STAT_label[STATL_CIRCULAR_STANDARD_DEVIATION] << "\t"
         << standard_deviation;
    }
    os << endl;
    break;
  }
  }

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'une loi continue au format tableur.
 *
 *  arguments : stream, flag sur l'ecriture de la fonction de repartition,
 *              pointeurs sur un objet Histogram et
 *              sur un objet FrequencyDistribution.
 *
 *--------------------------------------------------------------*/

ostream& ContinuousParametric::spreadsheet_print(ostream &os , bool cumul_flag ,
                                                 const Histogram *histo1 ,
                                                 const FrequencyDistribution *histo2)

{
  if ((histo1) || (histo2)) {
    register int i , j;
    int nb_step;
    double step , scale , value , mass , *frequency , *dist_cumul , *histo_cumul;


    if (histo1) {
      step = histo1->step;
      scale = histo1->nb_element;
    }
    else {
      step = histo2->min_interval_computation();
      scale = histo2->nb_element;
    }

    if (ident == GAMMA) {
      min_value = 0.;
      if (shape == 0.) {
        max_value = 0.;
      }
      else {
        max_value = shape * scale;
      }
    }

    else if (ident == ZERO_INFLATED_GAMMA) {
      min_value = 0.;
      if (zero_probability == 1.) {
        max_value = 0.;
      }
      else {
        max_value = shape * scale;
      }
    }

    else if ((ident == GAUSSIAN) || (ident == VON_MISES)) {
      min_value = location;

      if (histo1) {
        switch (histo1->type) {

        case INT_VALUE : {
          if (histo1->min_value < min_value) {
            min_value = histo1->min_value;
          }
          break;
        }

        case REAL_VALUE : {
          if (histo1->min_value - histo1->step / 2 < min_value) {
            min_value = histo1->min_value - histo1->step / 2;
          }
          break;
        }
        }
      }

      if ((histo2) && (histo2->offset < min_value)) {
        min_value = histo2->offset;
      }

      max_value = location;
    }

    if (histo1) {
      switch (histo1->type) {

      case INT_VALUE : {
        if (histo1->max_value > max_value) {
          max_value = histo1->max_value;
        }
        break;
      }

      case REAL_VALUE : {
        if (histo1->max_value + histo1->step / 2 > max_value) {
          max_value = histo1->max_value + histo1->step / 2;
        }
        break;
      }
      }
    }

    if ((histo2) && (histo2->nb_value - 1 > max_value)) {
      max_value = histo2->nb_value - 1;
    }

    if (cumul_flag) {
      if (histo1) {
        histo_cumul = histo1->cumul_computation();
      }
      else {
        histo_cumul = histo2->cumul_computation();
      }
    }

    // calcul d'une loi continue discretisee

    switch (ident) {

    case GAMMA : {
      if (shape == 0.) {
        dist_cumul = new double[1];
        frequency = new double[1];

        dist_cumul[0] = 1.;
        frequency[0] = scale;
      }

      else {
        gamma_distribution<double> dist(shape , scale);

        value = quantile(complement(dist , GAMMA_TAIL));
        while (max_value < value) {
          max_value += step;
        }

        nb_step = (int)((max_value - min_value) / step) + 1;
        dist_cumul = new double[nb_step];
        frequency = new double[nb_step];

//        value = step;
//        dist_cumul[0] = cdf(dist , value);
        value = step / 2;
        dist_cumul[0] = cdf(dist , value + step / 2);
        frequency[0] = dist_cumul[0] * scale;

        for (i = 1;i < nb_step;i++) {
          value += step;
//          dist_cumul[i] = cdf(dist , value);
          dist_cumul[i] = cdf(dist , value + step / 2);
          frequency[i] = (dist_cumul[i] - dist_cumul[i - 1]) * scale;
        }
      }
      break;
    }

    case ZERO_INFLATED_GAMMA : {
      if (zero_probability == 1.) {
        dist_cumul = new double[1];
        frequency = new double[1];

        dist_cumul[0] = 1.;
        frequency[0] = scale;
      }

      else {
        gamma_distribution<double> dist(shape , scale);

        value = quantile(complement(dist , GAMMA_TAIL));
        while (max_value < value) {
          max_value += step;
        }

        nb_step = (int)((max_value - min_value) / step) + 1;
        dist_cumul = new double[nb_step];
        frequency = new double[nb_step];

        value = 0.;
        dist_cumul[0] = zero_probability;
        frequency[0] = zero_probability * scale;

        for (i = 1;i < nb_step;i++) {
          value += step;
          dist_cumul[i] = zero_probability + (1. - zero_probability) * cdf(dist , value);
          frequency[i] = (dist_cumul[i] - dist_cumul[i - 1]) * scale;
        }
      }
      break;
    }

    case GAUSSIAN : {
      normal dist(location , dispersion);

      value = quantile(dist , GAUSSIAN_TAIL);
      while (min_value > value) {
        min_value -= step;
      }

      value = quantile(complement(dist , GAUSSIAN_TAIL));
      while (max_value < value) {
        max_value += step;
      }

      nb_step = (int)((max_value - min_value) / step) + 1;
      dist_cumul = new double[nb_step];
      frequency = new double[nb_step];

      value = min_value;
      dist_cumul[0] = cdf(dist , value + step / 2);
      frequency[0] = (dist_cumul[0] - cdf(dist , value - step / 2)) * scale;

      for (i = 1;i < nb_step;i++) {
        value += step;
        dist_cumul[i] = cdf(dist , value + step / 2);
        frequency[i] = (dist_cumul[i] - dist_cumul[i - 1]) * scale;
      }
      break;
    }

    case VON_MISES : {
      while (min_value >= step / 2) {
        min_value -= step;
      }

      value = (unit == DEGREE ? 360 : 2 * M_PI) - 3 * step / 2;
      while (max_value <= value) {
        max_value += step;
      }

      nb_step = (int)((max_value - min_value) / step) + 1;
      frequency = new double[nb_step];
      dist_cumul = new double[nb_step];

      value = min_value;
      mass = von_mises_mass_computation(value - step / 2 , value + step / 2);
      frequency[0] = mass * scale;
      dist_cumul[0] = mass;

      for (i = 1;i < nb_step;i++) {
        value += step;
        mass = von_mises_mass_computation(value - step / 2 , value + step / 2);
        frequency[i] = mass * scale;
        dist_cumul[i] = dist_cumul[i - 1] + mass;
      }
      break;
    }
    }

    value = min_value;
    if (histo1) {
      i = -1;
    }
    else {
      i = histo2->offset - step;
    }

    for (j = 0;j < nb_step;j++) {
      os << value << "\t";

      if ((histo1) && (value >= histo1->min_value) && (value <= histo1->max_value)) {
        os << histo1->frequency[++i];
      }

      if ((histo2) && (value >= histo2->offset) && (value < histo2->nb_value)) {
        i += step;
        os << histo2->frequency[i];
      }

      os << "\t" << frequency[j];

      if (cumul_flag) {
        os << "\t";
        if (((histo1) && (value >= histo1->min_value) && (value <= histo1->max_value)) ||
            ((histo2) && (value >= histo2->offset) && (value < histo2->nb_value))) {
          os << histo_cumul[i];
        }

        os << "\t" << dist_cumul[j];
      }

      os << endl;
      value += step;
    }

    delete [] frequency;
    delete [] dist_cumul;

    if (cumul_flag) {
      delete [] histo_cumul;
    }
  }

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture des parametres d'une loi continue au format Gnuplot.
 *
 *  argument : stream.
 *
 *--------------------------------------------------------------*/

ostream& ContinuousParametric::plot_title_print(ostream &os) const

{
  os << STAT_continuous_distribution_letter[ident] << "(";

  if (ident == GAMMA) {
    os << shape;
    if (shape > 0.) {
      os << ", " << scale;
    }
  }
  else if (ident == ZERO_INFLATED_GAMMA) {
    os << zero_probability;
    if (zero_probability < 1.) {
      os << ", " << shape << ", " << scale;
    }
  }
  else if (ident == LINEAR_MODEL) {
    os << intercept << ", " << slope << ", " << dispersion;
  }
  else {
    os << location << ", " << dispersion;
  }
  os << ")";

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'une loi continue au format Gnuplot.
 *
 *  arguments : path, pointeurs sur un objet Histogram et
 *              sur un objet FrequencyDistribution.
 *
 *--------------------------------------------------------------*/

bool ContinuousParametric::plot_print(const char *path , const Histogram *histo1 ,
                                      const FrequencyDistribution *histo2)

{
  bool status = false;
  register int i;
  int nb_step;
  double scale , step , buff , value , max , *frequency;
  ofstream out_file(path);


  if (out_file) {
    status = true;

    if (histo1) {
      scale = histo1->step * histo1->nb_element;
    }
    else if (histo2) {
      scale = histo2->nb_element;
    }
    else {
      scale = 1.;
    }

    switch (ident) {

    case GAMMA : {
      min_value = 0.;

      if (shape == 0.) {
        max_value = 0.;
        frequency = new double[1];
        frequency[0] = scale;
        max = frequency[0];
      }

      else {
        gamma_distribution<double> dist(shape , scale);

        max_value = quantile(complement(dist , GAMMA_TAIL));
        if ((histo1) && (histo1->max_value + histo1->step > max_value)) {
          max_value = histo1->max_value + histo1->step;
        }
        if ((histo2) && (histo2->nb_value - 1)) {
          max_value = histo2->nb_value - 1;
        }

        step = max_value / GAMMA_NB_STEP;
        if (histo1) {
          buff = histo1->step / GAMMA_NB_SUB_STEP;
        }
        else if (histo2) {
          buff = histo2->min_interval_computation() / GAMMA_NB_SUB_STEP;
        }
        if (((histo1) || (histo2)) && (buff < step)) {
          step = buff;
        }

        nb_step = (int)(max_value / step) + 1;
        frequency = new double[nb_step];

        value = 0.;
        max = 0.;
        for (i = 0;i < nb_step;i++) {
          frequency[i] = pdf(dist , value) * scale;
          value += step;
          if (frequency[i] > max) {
            max = frequency[i];
          }
        }
      }
      break;
    }

    case ZERO_INFLATED_GAMMA : {
      min_value = 0.;

      if (zero_probability == 1.) {
        max_value = 0.;
        frequency = new double[1];
        frequency[0] = scale;
        max = frequency[0];
      }

      else {
        gamma_distribution<double> dist(shape , scale);

        max_value = quantile(complement(dist , GAMMA_TAIL));
        if ((histo1) && (histo1->max_value + histo1->step > max_value)) {
          max_value = histo1->max_value + histo1->step;
        }
        if ((histo2) && (histo2->nb_value - 1)) {
          max_value = histo2->nb_value - 1;
        }

        step = max_value / GAMMA_NB_STEP;
        if (histo1) {
          buff = histo1->step / GAMMA_NB_SUB_STEP;
        }
        else if (histo2) {
          buff = histo2->min_interval_computation() / GAMMA_NB_SUB_STEP;
        }
        if (((histo1) || (histo2)) && (buff < step)) {
          step = buff;
        }

        nb_step = (int)(max_value / step) + 1;
        frequency = new double[nb_step];

        value = 0.;
        frequency[0] = zero_probability * scale;
        max = frequency[0];

        for (i = 1;i < nb_step;i++) {
          value += step;
          frequency[i] = (1. - zero_probability) * pdf(dist , value) * scale;
          if (frequency[i] > max) {
            max = frequency[i];
          }
        }
      }
      break;
    }

    case GAUSSIAN : {
      normal dist(location , dispersion);

      min_value = quantile(dist , GAUSSIAN_TAIL);
      if ((histo1) && (histo1->min_value - histo1->step < min_value)) {
        min_value = histo1->min_value - histo1->step;
      }
      if ((histo2) && (histo2->offset < min_value)) {
        min_value = histo2->offset;
      }

      max_value = quantile(complement(dist , GAUSSIAN_TAIL));
      if ((histo1) && (histo1->max_value + histo1->step > max_value)) {
        max_value = histo1->max_value + histo1->step;
      }
      if ((histo2) && (histo2->nb_value - 1)) {
        max_value = histo2->nb_value - 1;
      }

      step = (max_value - min_value) / GAUSSIAN_NB_STEP;
      if (histo1) {
        buff = histo1->step / GAUSSIAN_NB_SUB_STEP;
      }
      else if (histo2) {
        buff = histo2->min_interval_computation() / GAUSSIAN_NB_SUB_STEP;
      }
      if (((histo1) || (histo2)) && (buff < step)) {
        step = buff;
      }

/*     calcul complet

      nb_step = (int)((max_value - min_value) / step) + 1;
      frequency = new double[nb_step];

      value = min_value;
      for (i = 0;i < nb_step;i++) {
        frequency[i] = pdf(dist , value) * scale;
        value += step;
      }

      max = frequency[(int)round((location - min_value) / step)]; */

      nb_step = (int)(MAX(max_value - location , location - min_value) / step) + 1;
      frequency = new double[nb_step * 2 - 1];

      value = location;
      i = 0;
      while (value < MAX(max_value , 2 * location - min_value)) {
        frequency[i++] = pdf(dist , value) * scale;
        value += step;
      }

      max = frequency[0];

#     ifdef DEBUG
      cout << "\nNB_STEP: " << i << " | " << nb_step << endl;
#     endif

      nb_step = i;
      buff = frequency[nb_step - 1];
      for (i = 0;i < nb_step - 1;i++) {
        frequency[nb_step - 1 + i] = frequency[i];
      }
      frequency[2 * (nb_step - 1)] = buff;
      for (i = 1;i < nb_step;i++) {
        frequency[nb_step - 1 - i] = frequency[nb_step - 1 + i];
      }

      value = location - (nb_step - 1) * step;

#     ifdef DEBUG
      cout << "\nMIN_VALUE: " << value << " | " << min_value << endl;
#     endif

      nb_step = 2 * nb_step - 1;
      break;
    }

    case VON_MISES : {
      frequency = new double[VON_MISES_NB_STEP];
      min_value = 0.;

      switch (unit) {

      case DEGREE : {
        max_value = 360;
        step =  max_value / VON_MISES_NB_STEP;

        value = 0.;
        for (i = 0;i < VON_MISES_NB_STEP;i++) {
          frequency[i] = exp(dispersion * cos((value - location) * M_PI / 180)) * scale /
                         (360 * cyl_bessel_i(0 , dispersion));
          value += step;
        }
        break;
      }

      case RADIAN : {
        max_value = 2 * M_PI;
        step =  max_value / VON_MISES_NB_STEP;

        value = 0.;
        for (i = 0;i < VON_MISES_NB_STEP;i++) {
          frequency[i] = exp(dispersion * cos(value - location)) * scale /
                         (2 * M_PI * cyl_bessel_i(0 , dispersion));
          value += step;
        }
        break;
      }
      }

      max = frequency[(int)round(location / step)];

      nb_step = VON_MISES_NB_STEP;
      value = 0.;
      break;
    }
    }

    for (i = 0;i < nb_step;i++) {
      out_file << value << " " << frequency[i] << endl;
      value += step;
    }

    delete [] frequency;
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un q-q plot au format Gnuplot.
 *
 *  arguments : path, nombre de valeurs, pointeur sur la fonction
 *              de repartition empirique.
 *
 *--------------------------------------------------------------*/

bool ContinuousParametric::q_q_plot_print(const char *path , int nb_value ,
                                          double **empirical_cdf) const

{
  bool status = false;
  register int i;
  double **qqplot;
  ofstream out_file(path);


  if (out_file) {
    status = true;

    qqplot = q_q_plot_computation(nb_value , empirical_cdf);

    for (i = 0;i < (ident == VON_MISES ? nb_value : nb_value - 1);i++) {
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
 *  Sortie graphique d'une loi continue.
 *
 *  arguments : reference sur un objet SinglePlot,
 *              pointeurs sur un objet Histogram et
 *              sur un objet FrequencyDistribution.
 *
 *--------------------------------------------------------------*/

void ContinuousParametric::plotable_write(SinglePlot &plot , const Histogram *histo1 ,
                                          const FrequencyDistribution *histo2)

{
  register int i;
  int nb_step;
  double scale , step , buff , value , max , *frequency;


  if (histo1) {
    scale = histo1->step * histo1->nb_element;
  }
  else if (histo2) {
    scale = histo2->nb_element;
  }
  else {
    scale = 1.;
  }

  switch (ident) {

  case GAMMA : {
    min_value = 0.;

    if (shape == 0.) {
      max_value = 0.;
      frequency = new double[1];
      frequency[0] = scale;
      max = frequency[0];
    }

    else {
      gamma_distribution<double> dist(shape , scale);

      max_value = quantile(complement(dist , GAMMA_TAIL));
      if ((histo1) && (histo1->max_value + histo1->step > max_value)) {
        max_value = histo1->max_value + histo1->step;
      }
      if ((histo2) && (histo2->nb_value - 1)) {
        max_value = histo2->nb_value - 1;
      }

      step = max_value / GAMMA_NB_STEP;
      if (histo1) {
        buff = histo1->step / GAMMA_NB_SUB_STEP;
      }
      else if (histo2) {
        buff = histo2->min_interval_computation() / GAMMA_NB_SUB_STEP;
      }
      if (((histo1) || (histo2)) && (buff < step)) {
        step = buff;
      }

      nb_step = (int)(max_value / step) + 1;
      frequency = new double[nb_step];

      value = 0.;
      max = 0.;
      for (i = 0;i < nb_step;i++) {
        frequency[i] = pdf(dist , value) * scale;
        value += step;

        if (frequency[i] > max) {
          max = frequency[i];
        }
      }
    }
    break;
  }

  case ZERO_INFLATED_GAMMA : {
    min_value = 0.;

    if (zero_probability == 1.) {
      max_value = 0.;
      frequency = new double[1];
      frequency[0] = scale;
      max = frequency[0];
    }

    else {
      gamma_distribution<double> dist(shape , scale);

      max_value = quantile(complement(dist , GAMMA_TAIL));
      if ((histo1) && (histo1->max_value + histo1->step > max_value)) {
        max_value = histo1->max_value + histo1->step;
      }
      if ((histo2) && (histo2->nb_value - 1)) {
        max_value = histo2->nb_value - 1;
      }

      step = max_value / GAMMA_NB_STEP;
      if (histo1) {
        buff = histo1->step / GAMMA_NB_SUB_STEP;
      }
      else if (histo2) {
        buff = histo2->min_interval_computation() / GAMMA_NB_SUB_STEP;
      }
      if (((histo1) || (histo2)) && (buff < step)) {
        step = buff;
      }

      nb_step = (int)(max_value / step) + 1;
      frequency = new double[nb_step];

      value = 0.;
      frequency[0] = zero_probability * scale;
      max = frequency[0];

      for (i = 1;i < nb_step;i++) {
        value += step;
        frequency[i] = (1. - zero_probability) * pdf(dist , value) * scale;

        if (frequency[i] > max) {
          max = frequency[i];
        }
      }
    }
    break;
  }

  case GAUSSIAN : {
    normal dist(location , dispersion);

    min_value = quantile(dist , GAUSSIAN_TAIL);
    if ((histo1) && (histo1->min_value - histo1->step < min_value)) {
      min_value = histo1->min_value - histo1->step;
    }
    if ((histo2) && (histo2->offset < min_value)) {
      min_value = histo2->offset;
    }

    max_value = quantile(complement(dist , GAUSSIAN_TAIL));
    if ((histo1) && (histo1->max_value + histo1->step > max_value)) {
      max_value = histo1->max_value + histo1->step;
    }
    if ((histo2) && (histo2->nb_value - 1)) {
      max_value = histo2->nb_value - 1;
    }

    step = (max_value - min_value) / GAUSSIAN_NB_STEP;
    if (histo1) {
      buff = histo1->step / GAUSSIAN_NB_SUB_STEP;
    }
    else if (histo2) {
      buff = histo2->min_interval_computation() / GAUSSIAN_NB_SUB_STEP;
    }
    if (((histo1) || (histo2)) && (buff < step)) {
      step = buff;
    }

    nb_step = (int)(MAX(max_value - location , location - min_value) / step) + 1;
    frequency = new double[nb_step * 2 - 1];

    value = location;
    i = 0;
    while (value < MAX(max_value , 2 * location - min_value)) {
      frequency[i++] = pdf(dist , value) * scale;
      value += step;
    }

    max = frequency[0];

#   ifdef DEBUG
    cout << "\nNB_STEP: " << i << " | " << nb_step << endl;
#   endif

    nb_step = i;
    buff = frequency[nb_step - 1];
    for (i = 0;i < nb_step - 1;i++) {
      frequency[nb_step - 1 + i] = frequency[i];
    }
    frequency[2 * (nb_step - 1)] = buff;
    for (i = 1;i < nb_step;i++) {
      frequency[nb_step - 1 - i] = frequency[nb_step - 1 + i];
    }

    value = location - (nb_step - 1) * step;

#   ifdef DEBUG
    cout << "\nMIN_VALUE: " << value << " | " << min_value << endl;
#   endif

    nb_step = 2 * nb_step - 1;
    break;
  }

  case VON_MISES : {
    frequency = new double[VON_MISES_NB_STEP];
    min_value = 0.;

    switch (unit) {

    case DEGREE : {
      max_value = 360;
      step = max_value / VON_MISES_NB_STEP;

      value = 0.;
      for (i = 0;i < VON_MISES_NB_STEP;i++) {
        frequency[i] = exp(dispersion * cos((value - location) * M_PI / 180)) * scale /
                       (360 * cyl_bessel_i(0 , dispersion));
        value += step;
      }
      break;
    }

    case RADIAN : {
      max_value = 2 * M_PI;
      step = max_value / VON_MISES_NB_STEP;

      value = 0.;
      for (i = 0;i < VON_MISES_NB_STEP;i++) {
        frequency[i] = exp(dispersion * cos(value - location)) * scale /
                       (2 * M_PI * cyl_bessel_i(0 , dispersion));
        value += step;
      }
      break;
    }
    }

    max = frequency[(int)round(location / step)];

    nb_step = VON_MISES_NB_STEP;
    value = 0.;
    break;
  }
  }

  for (i = 0;i < nb_step;i++) {
    plot.add_point(value , frequency[i]);
    value += step;
  }

  delete [] frequency;
}


/*--------------------------------------------------------------*
 *
 *  Sortie graphique d'un q-q plot.
 *
 *  arguments : reference sur un objet SinglePlot, nombre de valeurs,
 *              pointeur sur la fonction de repartition empirique.
 *
 *--------------------------------------------------------------*/

void ContinuousParametric::q_q_plotable_write(SinglePlot &plot , int nb_value ,
                                              double **empirical_cdf) const

{
  register int i;
  double **qqplot;


  qqplot = q_q_plot_computation(nb_value , empirical_cdf);

  for (i = 0;i < (ident == VON_MISES ? nb_value : nb_value - 1);i++) {
    plot.add_point(qqplot[0][i] , qqplot[1][i]);
  }

  delete [] qqplot[0];
  delete [] qqplot[1];
  delete [] qqplot;
}


/*--------------------------------------------------------------*
 *
 *  Calcul du nombre de parametres d'une loi continue.
 *
 *--------------------------------------------------------------*/

int ContinuousParametric::nb_parameter_computation() const

{
  int nb_parameter;


  switch (ident) {

  case GAMMA : {
    if (shape == 0.) {
      nb_parameter = 1;
    }
    else {
      nb_parameter = 2;
    }
    break;
  }

  case ZERO_INFLATED_GAMMA : {
    if (zero_probability == 1.) {
      nb_parameter = 1;
    }
    else if (zero_probability == 0.) {
      nb_parameter = 2;
    }
    else {
      nb_parameter = 3;
    }
    break;
  }

  case GAUSSIAN : {
    nb_parameter = 2;
    break;
  }

  case LINEAR_MODEL : {
    nb_parameter = 3;
    break;
  }

  case VON_MISES : {
    nb_parameter = 2;
    break;
  }

  default : {
    nb_parameter = 0;
    break;
  }
  }

  return nb_parameter;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la densite d'une loi de von Mises integree sur un intervalle.
 *
 *  arguments : bornes de l'intervalle.
 *
 *--------------------------------------------------------------*/

double ContinuousParametric::von_mises_mass_computation(double inf , double sup) const

{
  double step , mass , value;


//  if (dispersion < VON_MISES_CONCENTRATION_THRESHOLD) {
    mass = 0.;

    switch (unit) {

    case DEGREE : {
      step = MIN((sup - inf) / VON_MISES_NB_SUB_STEP , 360. / VON_MISES_NB_STEP);

      for (value = inf + step / 2;value < sup;value += step) {
        mass += exp(dispersion * cos((value - location) * M_PI / 180));
      }
      mass = mass * step / (360 * cyl_bessel_i(0 , dispersion));
      break;
    }

    case RADIAN : {
      step = MIN((sup - inf) / VON_MISES_NB_SUB_STEP , 2 * M_PI / VON_MISES_NB_STEP);

      for (value = inf + step;value < sup;value += step) {
        mass += exp(dispersion * cos(value - location));
      }
      mass = mass * step / (2 * M_PI * cyl_bessel_i(0 , dispersion));
      break;
    }
    }
/*  }

  else {
    normal dist(location , sqrt(1. / dispersion));

    mass = cdf(dist , sup) - cdf(dist , inf);
  } */

  return mass;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la densite integree sur un intervalle.
 *
 *  arguments : bornes de l'intervalle.
 *
 *--------------------------------------------------------------*/

double ContinuousParametric::mass_computation(double inf , double sup) const

{
  double mass;


  switch (ident) {

  case GAMMA : {
    if (inf < 0.) {
      inf = 0.;
    }
    if (sup < inf) {
      sup = inf;
    }

    if (shape == 0.) {
      if (inf == 0.) {
        mass = 1.;
      }
      else {
        mass = 0.;
      }
    }

    else {
      gamma_distribution<double> dist(shape , scale);

      if (inf == sup) {
        mass = pdf(dist , MAX(inf , 1.e-12));  // bug boost C++

#       ifdef DEBUG
        if ((shape == 1.) && (scale < 0.1) && (inf == 0.)) {
          cout << "TEST: " << scale << ", " << mass << endl;
        }
#       endif

      }

      else {
        mass = cdf(dist , sup) - cdf(dist , inf);
      }
    }
    break;
  }

  case ZERO_INFLATED_GAMMA : {
    if (inf < 0.) {
      inf = 0.;
    }
    if (sup < inf) {
      sup = inf;
    }

    if (zero_probability == 1.) {
      if (inf == 0.) {
        mass = zero_probability;
      }
      else {
        mass = 0.;
      }
    }

    else {
      gamma_distribution<double> dist(shape , scale);

      if (inf == 0.) {
        mass = zero_probability;
      }

      else {
        if (inf == sup) {
          mass = (1. - zero_probability) * pdf(dist , inf);
        }
        else {
          mass = (1. - zero_probability) * (cdf(dist , sup) - cdf(dist , inf));
        }
      }

      if (mass > 1.) {

#       ifdef MESSAGE
        cout << "WARNING: " << inf << " " << sup << " " << mass << " | "
             << zero_probability << " " << shape << " " << scale << endl;
#       endif

      }
    }
    break;
  }

  case GAUSSIAN : {
    normal dist(location , dispersion);

    if (inf == sup) {
      mass = pdf(dist , inf);
    }
    else {
      mass = cdf(dist , sup) - cdf(dist , inf);
    }
    break;
  }

  case VON_MISES : {
    if (inf == sup) {
      switch (unit) {
      case DEGREE :
        mass = exp(dispersion * cos((inf - location) * M_PI / 180)) /
               (360 * cyl_bessel_i(0 , dispersion));
        break;
      case RADIAN :
        mass = exp(dispersion * cos(inf - location)) /
               (2 * M_PI * cyl_bessel_i(0 , dispersion));
        break;
      }
    }

    else {
      mass = von_mises_mass_computation(inf , sup);
    }
    break;
  }

  case LINEAR_MODEL : {
    normal dist(0. , dispersion);

    if (inf == sup) {
      mass = pdf(dist , inf);
    }
    else {
      mass = cdf(dist , sup) - cdf(dist , inf);
    }
    break;
  }
  }

//  if (mass > 1.) {
  if ((inf == sup) && (mass > 1.)) {
    mass = 1.;
  }

  return mass;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la fonction de repartition d'une loi de von Mises.
 *
 *--------------------------------------------------------------*/

void ContinuousParametric::von_mises_cumul_computation()

{
  register int i , j;
  int start;
  double step , value;


  cumul = new double[VON_MISES_NB_STEP];

  switch (unit) {

  case DEGREE : {
    step = 360. / VON_MISES_NB_STEP;

    if (location < 180) {
      value = location + 180;
    }
    else {
      value = location - 180;
    }
    start = (int)round(value * VON_MISES_NB_STEP / 360);

    cumul[start] = exp(dispersion * cos((value - location) * M_PI / 180)) * step /
                   (360 * cyl_bessel_i(0 , dispersion));
    for (i = start + 1;i < VON_MISES_NB_STEP;i++) {
      value += step;
      cumul[i] = cumul[i - 1] + exp(dispersion * cos((value - location) * M_PI / 180)) * step /
                                (360 * cyl_bessel_i(0 , dispersion));
    }

    value += step - 360;
    cumul[0] = cumul[VON_MISES_NB_STEP - 1] + exp(dispersion * cos((value - location) * M_PI / 180)) * step /
                                              (360 * cyl_bessel_i(0 , dispersion));
    for (i = 1;i < start;i++) {
      value += step;
      cumul[i] = cumul[i - 1] + exp(dispersion * cos((value - location) * M_PI / 180)) * step /
                                (360 * cyl_bessel_i(0 , dispersion));
    }
    break;
  }

  case RADIAN : {
    step = 2 * M_PI / VON_MISES_NB_STEP;

    if (location < M_PI) {
      value = location + M_PI;
    }
    else {
      value = location - M_PI;
    }
    start = (int)round(value * VON_MISES_NB_STEP / (2 * M_PI));

    cumul[start] = exp(dispersion * cos(value - location)) * step /
                   (2 * M_PI * cyl_bessel_i(0 , dispersion));
    for (i = start + 1;i < VON_MISES_NB_STEP;i++) {
      value += step;
      cumul[i] = cumul[i - 1] + exp(dispersion * cos(value - location)) * step /
                                (2 * M_PI * cyl_bessel_i(0 , dispersion));
    }

    value += step - 2 * M_PI;
    cumul[0] = cumul[VON_MISES_NB_STEP - 1] + exp(dispersion * cos(value - location)) * step /
                                              (2 * M_PI * cyl_bessel_i(0 , dispersion));
    for (i = 1;i < start;i++) {
      value += step;
      cumul[i] = cumul[i - 1] + exp(dispersion * cos(value - location)) * step /
                                (2 * M_PI * cyl_bessel_i(0 , dispersion));
    }
    break;
  }
  }
}


/*--------------------------------------------------------------*
 *
 *  Calcul d'un q-q plot.
 *
 *  arguments : nombre de valeurs, pointeur sur la fonction de repartition empirique.
 *
 *--------------------------------------------------------------*/

double** ContinuousParametric::q_q_plot_computation(int nb_value ,
                                                    double **empirical_cdf) const

{
  register int i , j;
  double step , value , current_cumul , previous_cumul , **qqplot;


  qqplot = new double*[2];
  qqplot[0] = new double[nb_value];
  qqplot[1] = new double[nb_value];

  for (i = 0;i < nb_value;i++) {
    qqplot[0][i] = empirical_cdf[0][i];
  }

  switch (ident) {

  case GAMMA : {
    if (shape == 0.) {
      for (i = 0;i < nb_value - 1;i++) {
        qqplot[1][i] = 0.;
      }
    }

    else {
      gamma_distribution<double> dist(shape , scale);

      for (i = 0;i < nb_value - 1;i++) {
        qqplot[1][i] = quantile(dist , empirical_cdf[1][i]);
      }
    }
    break;
  }

  case ZERO_INFLATED_GAMMA : {
    if (zero_probability == 1.) {
      for (i = 0;i < nb_value - 1;i++) {
        qqplot[1][i] = 0.;
      }
    }

    else {
      gamma_distribution<double> dist(shape , scale);

      i = 0;
      while (empirical_cdf[1][i] <= zero_probability) {
        qqplot[1][i++] = 0.;
      }
      for (j = i;j < nb_value - 1;i++) {
        qqplot[1][j] = quantile(dist , (empirical_cdf[1][j] - zero_probability) / (1. - zero_probability));
      }
    }
    break;
  }

  case GAUSSIAN : {
    normal dist(location , dispersion);

    for (i = 0;i < nb_value - 1;i++) {
      qqplot[1][i] = quantile(dist , empirical_cdf[1][i]);
    }
    break;
  }

  case VON_MISES : {
    current_cumul = 0.;
    value = 0.;

    switch (unit) {

    case DEGREE : {
      step = 360. / VON_MISES_NB_STEP;

      for (i = 0;i < nb_value - 1;i++) {
        while (current_cumul < empirical_cdf[1][i]) {
          value += step;
          previous_cumul = current_cumul;
          current_cumul += exp(dispersion * cos((value - step / 2 - location) * M_PI / 180)) * step /
                           (360 * cyl_bessel_i(0 , dispersion));
        }

        qqplot[1][i] = value - step * (current_cumul - empirical_cdf[1][i]) /
                       (current_cumul - previous_cumul);
      }

      qqplot[1][nb_value - 1] = 360;
      break;
    }

    case RADIAN : {
      step = 2 * M_PI / VON_MISES_NB_STEP;

      for (i = 0;i < nb_value - 1;i++) {
        while (current_cumul < empirical_cdf[1][i]) {
          value += step;
          previous_cumul = current_cumul;
          current_cumul += exp(dispersion * cos(value - step / 2 - location)) * step /
                           (2 * M_PI * cyl_bessel_i(0 , dispersion));
        }

        qqplot[1][i] = value - step * (current_cumul - empirical_cdf[1][i]) /
                       (current_cumul - previous_cumul);
      }

      qqplot[1][nb_value - 1] = 2 * M_PI;
      break;
    }
    }
    break;
  }
  }

  return qqplot;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la distance entre deux lois continues (sup de la difference
 *  absolue des fonctions de repartition dans le cas de fonctions de repartition
 *  ne se croisant pas; sinon somme des sup avant et apres croisement
 *  de la difference des fonctions de repartition).
 *
 *--------------------------------------------------------------*/

double ContinuousParametric::sup_norm_distance_computation(ContinuousParametric &dist)

{
  bool crossing;
  register int i;
  double min , max , step , value , buff , distance , max_absolute_diff , cumul1[2] , cumul2[2], overlap;

  switch (ident) {

  case GAMMA : {
    if ((shape == 0.) || (dist.shape == 0.)) {
      if ((shape == 0.) && (dist.shape == 0.)) {
        distance = 0.;
      }
      else {
        distance = 1.;
      }
    }

    else {
      gamma_distribution<double> dist1(shape , scale) , dist2(dist.shape , dist.scale);

      step = MIN(quantile(complement(dist1 , GAMMA_TAIL)) , quantile(complement(dist2 , GAMMA_TAIL))) / GAMMA_NB_STEP;
      distance = 0.;
      value = 0.;
      cumul1[0] = cdf(dist1 , value);
      cumul2[0] = cdf(dist2 , value);
      max_absolute_diff = fabs(cumul1[0] - cumul2[0]);
      crossing = false;

      for (i = 1;i < GAMMA_NB_STEP;i++) {
        value += step;
        cumul1[1] = cdf(dist1 , value);
        cumul2[1] = cdf(dist2 , value);

        buff = fabs(cumul1[1] - cumul2[1]);
        if (buff > max_absolute_diff) {
          max_absolute_diff = buff;
        }

        if ((!crossing) && (((cumul1[1] > cumul2[1]) && (cumul1[0] <= cumul2[0])) ||
            ((cumul1[1] <= cumul2[1]) && (cumul1[0] > cumul2[0])))) {
          crossing = true;
          distance = max_absolute_diff;
          max_absolute_diff = 0.;
        }

        cumul1[0] = cumul1[1];
        cumul2[0] = cumul2[1];
      }
      distance += max_absolute_diff;

#     ifdef DEBUG
      value = 0.;
      cumul1[0] = 0.;
      cumul2[0] = 0.;
      overlap = 0.;

      for (i = 0;i < GAMMA_NB_STEP;i++) {
        cumul1[1] = cdf(dist1 , value);
        cumul2[1] = cdf(dist2 , value);
        value += step;

        overlap += MIN(cumul1[1] - cumul1[0] , cumul2[1] - cumul2[0]);
        cumul1[0] = cumul1[1];
        cumul2[0] = cumul2[1];
      }

      cout << "\nSup norm distance: " << distance << " " << 1. - overlap << endl;
#     endif

    }
    break;
  }

  case ZERO_INFLATED_GAMMA : {
    cumul1[0] = zero_probability;
    cumul2[0] = dist.zero_probability;
    max_absolute_diff = fabs(cumul1[0] - cumul2[0]);

    if ((zero_probability == 1.) || (dist.zero_probability == 1.)) {
      distance = max_absolute_diff;

#     ifdef MESSAGE
      cout << "\nSup norm distance: " << distance << " " << 1. - MIN(zero_probability , dist.zero_probability) << endl;
#     endif

    }

    else {
      gamma_distribution<double> dist1(shape , scale) , dist2(dist.shape , dist.scale);

      step = MIN(quantile(complement(dist1 , GAMMA_TAIL)) , quantile(complement(dist2 , GAMMA_TAIL))) / GAMMA_NB_STEP;
      value = 0.;
      distance = 0.;
      crossing = false;

      for (i = 1;i < GAMMA_NB_STEP;i++) {
        value += step;
        cumul1[1] = zero_probability + (1. - zero_probability) * cdf(dist1 , value);
        cumul2[1] = dist.zero_probability + (1. - dist.zero_probability) * cdf(dist2 , value);

        buff = fabs(cumul1[1] - cumul2[1]);
        if (buff > max_absolute_diff) {
          max_absolute_diff = buff;
        }

        if ((!crossing) && (((cumul1[1] > cumul2[1]) && (cumul1[0] <= cumul2[0])) ||
            ((cumul1[1] <= cumul2[1]) && (cumul1[0] > cumul2[0])))) {
          crossing = true;
          distance = max_absolute_diff;
          max_absolute_diff = 0.;
        }

        cumul1[0] = cumul1[1];
        cumul2[0] = cumul2[1];
      }
      distance += max_absolute_diff;

#     ifdef MESSAGE
      value = 0.;
      cumul1[0] = zero_probability;
      cumul2[0] = dist.zero_probability;
      overlap = MIN(zero_probability , dist.zero_probability);

      for (i = 1;i < GAMMA_NB_STEP;i++) {
        value += step;
        cumul1[1] = zero_probability + (1. - zero_probability) * cdf(dist1 , value);
        cumul2[1] = dist.zero_probability + (1. - dist.zero_probability) * cdf(dist2 , value);

        overlap += MIN(cumul1[1] - cumul1[0] , cumul2[1] - cumul2[0]);
        cumul1[0] = cumul1[1];
        cumul2[0] = cumul2[1];
      }

      cout << "\nSup norm distance: " << distance << " " << 1. - overlap << endl;
#     endif

    }

    break;
  }

  case GAUSSIAN : {
    normal dist1(location , dispersion) , dist2(dist.location , dist.dispersion);

    min = MAX(quantile(dist1 , GAUSSIAN_TAIL) , quantile(dist2 , GAUSSIAN_TAIL));
    max = MIN(quantile(complement(dist1 , GAUSSIAN_TAIL)) , quantile(complement(dist2 , GAUSSIAN_TAIL)));

    step = (max - min) / GAUSSIAN_NB_STEP;
    distance = 0.;
    value = min;
    cumul1[0] = cdf(dist1 , value);
    cumul2[0] = cdf(dist2 , value);
    max_absolute_diff = fabs(cumul1[0] - cumul2[0]);
    crossing = false;

    for (i = 1;i < GAUSSIAN_NB_STEP;i++) {
      value += step;
      cumul1[1] = cdf(dist1 , value);
      cumul2[1] = cdf(dist2 , value);

      buff = fabs(cumul1[1] - cumul2[1]);
      if (buff > max_absolute_diff) {
        max_absolute_diff = buff;
      }

      if ((!crossing) && (((cumul1[1] > cumul2[1]) && (cumul1[0] <= cumul2[0])) ||
          ((cumul1[1] <= cumul2[1]) && (cumul1[0] > cumul2[0])))) {
        crossing = true;
        distance = max_absolute_diff;
        max_absolute_diff = 0.;
      }

      cumul1[0] = cumul1[1];
      cumul2[0] = cumul2[1];
    }
    distance += max_absolute_diff;

#   ifdef DEBUG
    value = min;
    cumul1[0] = 0.;
    cumul2[0] = 0.;
    overlap = 0.;

    for (i = 0;i < GAUSSIAN_NB_STEP;i++) {
      cumul1[1] = cdf(dist1 , value);
      cumul2[1] = cdf(dist2 , value);
      value += step;

      overlap += MIN(cumul1[1] - cumul1[0] , cumul2[1] - cumul2[0]);
      cumul1[0] = cumul1[1];
      cumul2[0] = cumul2[1];
    }

    cout << "\nSup norm distance: " << distance << " " << 1. - overlap << endl;
#   endif

    break;
  }

  case VON_MISES : {   // bug concentration faible
    int inf , sup;

    if (!cumul) {
      von_mises_cumul_computation();
    }
    if (!dist.cumul) {
      dist.von_mises_cumul_computation();
    }

    switch (unit) {

    case DEGREE : {
      if (location < dist.location) {
        if (dist.location - location <= 180) {
          min = dist.location - 180;
          max = location + 180;
        }
        else {
          min = location - 180;
          max = dist.location + 180;
        }
      }

      else {
        if (location - dist.location <= 180) {
          min = location - 180;
          max = dist.location + 180;
        }
        else {
          min = dist.location - 180;
          max = location + 180;
        }
      }

      if (min < 0) {
        min += 360;
      }
      if (max >= 360) {
        max -= 360;
      }

      inf = (int)round(min * VON_MISES_NB_STEP / 360);
      sup = (int)round(max * VON_MISES_NB_STEP / 360);
      break;
    }

    case RADIAN : {
      if (location < dist.location) {
        if (dist.location - location <= M_PI) {
          min = dist.location - M_PI;
          max = location + M_PI;
        }
        else {
          min = location - M_PI;
          max = dist.location + M_PI;
        }
      }

      else {
        if (location - dist.location <= M_PI) {
          min = location - M_PI;
          max = dist.location + M_PI;
        }
        else {
          min = dist.location - M_PI;
          max = location + M_PI;
        }
      }

      if (min < 0) {
        min += 2 * M_PI;
      }
      if (max >= 2 * M_PI) {
        max -= 2 * M_PI;
      }

      inf = (int)round(min * VON_MISES_NB_STEP / (2 * M_PI));
      sup = (int)round(max * VON_MISES_NB_STEP / (2 * M_PI));
      break;
    }
    }

#   ifdef DEBUG
    cout << "\nInf: " << inf << ",  Sup: " << sup << endl;
#   endif

    distance = 0.;
    max_absolute_diff = fabs(cumul[inf] - dist.cumul[inf]);
    crossing = false;

    if (inf <= sup) {
      for (i = inf + 1;i < sup;i++) {
        buff = fabs(cumul[i] - dist.cumul[i]);
        if (buff > max_absolute_diff) {
          max_absolute_diff = buff;
        }

        if ((!crossing) && (((cumul[i] > dist.cumul[i]) && (cumul[i - 1] <= dist.cumul[i - 1])) ||
            ((cumul[i] <= dist.cumul[i]) && (cumul[i - 1] > dist.cumul[i - 1])))) {

#         ifdef DEBUG
          cout << "\nCrossing: " << i << endl;
#         endif

          crossing = true;
          distance = max_absolute_diff;
          max_absolute_diff = 0.;
        }
      }
    }

    else {
      for (i = inf + 1;i < VON_MISES_NB_STEP;i++) {
        buff = fabs(cumul[i] - dist.cumul[i]);
        if (buff > max_absolute_diff) {
          max_absolute_diff = buff;
        }

        if ((!crossing) && (((cumul[i] > dist.cumul[i]) && (cumul[i - 1] <= dist.cumul[i - 1])) ||
            ((cumul[i] <= dist.cumul[i]) && (cumul[i - 1] > dist.cumul[i - 1])))) {
          crossing = true;
          distance = max_absolute_diff;
          max_absolute_diff = 0.;
        }
      }

      buff = fabs(cumul[0] - dist.cumul[0]);
      if (buff > max_absolute_diff) {
        max_absolute_diff = buff;
      }

      if ((!crossing) && (((cumul[0] > dist.cumul[0]) && (cumul[VON_MISES_NB_STEP - 1] <= dist.cumul[VON_MISES_NB_STEP - 1])) ||
          ((cumul[0] <= dist.cumul[0]) && (cumul[VON_MISES_NB_STEP - 1] > dist.cumul[VON_MISES_NB_STEP - 1])))) {
        crossing = true;
        distance = max_absolute_diff;
        max_absolute_diff = 0.;
      }

      for (i = 1;i < sup;i++) {
        buff = fabs(cumul[i] - dist.cumul[i]);
        if (buff > max_absolute_diff) {
          max_absolute_diff = buff;
        }

        if ((!crossing) && (((cumul[i] > dist.cumul[i]) && (cumul[i - 1] <= dist.cumul[i - 1])) ||
            ((cumul[i] <= dist.cumul[i]) && (cumul[i - 1] > dist.cumul[i - 1])))) {
          crossing = true;
          distance = max_absolute_diff;
          max_absolute_diff = 0.;
        }
      }
    }

    distance += max_absolute_diff;

#   ifdef MESSAGE
    step = (unit == DEGREE ? 360. : 2 * M_PI) / VON_MISES_NB_STEP;
    value = 0.;
    overlap = 0.;

    for (i = 0;i < VON_MISES_NB_STEP;i++) {
      overlap += MIN(von_mises_mass_computation(value , value + step) , dist.von_mises_mass_computation(value , value + step));
      value += step;
    }

    cout << "\nSup norm distance: " << distance << " " << 1. - overlap << endl;
#   endif

    break;
  }
  }

  return distance;
}


/*--------------------------------------------------------------*
 *
 *  Simulation d'une loi continue.
 *
 *--------------------------------------------------------------*/

double ContinuousParametric::simulation()

{
  register int i;
  int start;
  double limit , value , step , current_cumul , previous_cumul;


  limit = ((double)rand() / (RAND_MAX + 1.));

  switch (ident) {

  case GAMMA : {
    if (shape == 0.) {
      value = 0.;
    }

    else {
      gamma_distribution<double> dist(shape , scale);

      value = quantile(dist , limit);
    }
    break;
  }

  case ZERO_INFLATED_GAMMA : {
    if (limit <= zero_probability) {
      value = 0.;
    }

    else if (zero_probability < 1.) {
      gamma_distribution<double> dist(shape , scale);

      value = quantile(dist , (limit - zero_probability) / (1. - zero_probability));
    }
    break;
  }

  case GAUSSIAN : {
    normal dist(location , dispersion);

    value = quantile(dist , limit);
    break;
  }

  case VON_MISES : {
    if (!cumul) {
      von_mises_cumul_computation();
    }

    switch (unit) {

    case DEGREE : {
      step = 360. / VON_MISES_NB_STEP;

      if (location < 180) {
        value = location + 180.;
      }
      else {
        value = location - 180.;
      }

      start = (int)round(value * VON_MISES_NB_STEP / 360);
      break;
    }

    case RADIAN : {
      step = 2 * M_PI / VON_MISES_NB_STEP;

      if (location < M_PI) {
        value = location + M_PI;
      }
      else {
        value = location - M_PI;
      }

      start = (int)round(value * VON_MISES_NB_STEP / (2 * M_PI));
      break;
    }
    }

    for (i = start;i < VON_MISES_NB_STEP;i++) {
      if (cumul[i] >= limit) {
        if (i > start) {
          value -= step * (cumul[i] - limit) / (cumul[i] - cumul[i - 1]);
        }
        break;
      }

      else {
        value += step;
      }
    }

    if (i == VON_MISES_NB_STEP) {
      for (i = 0;i < start;i++) {
        if (cumul[i] >= limit) {
          if (i == 0) {
            value -= step * (cumul[i] - limit) / (cumul[i] - cumul[VON_MISES_NB_STEP - 1]);
          }
          else {
            value -= step * (cumul[i] - limit) / (cumul[i] - cumul[i - 1]);
          }
          break;
        }

        else {
          value += step;
        }
      }
    }

    switch (unit) {

    case DEGREE : {
      if (value > 360) {
        value -= 360;
      }
      break;
    }

    case RADIAN : {
      if (value > 2 * M_PI) {
        value -= 2 * M_PI;
      }
      break;
    }
    }
    break;
  }

  case LINEAR_MODEL : {
    normal dist(0. , dispersion);

    value = quantile(dist , limit);
    break;
  }
  }

  return value;
}


};  // namespace stat_tool
