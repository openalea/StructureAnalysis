/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       V-Plants: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2010 CIRAD/INRIA Virtual Plants
 *
 *       File author(s): Y. Guedon (yann.guedon@cirad.fr)
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


extern int column_width(int value);
extern int column_width(double min_value , double max_value);
extern int column_width(int nb_value , const double *value , double scale = 1.);



/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe ContinuousParametric.
 *
 *  arguments : identificateur, parametres de moyenne et de dispersion,
 *              unite (loi de von Mises).
 *
 *--------------------------------------------------------------*/

ContinuousParametric::ContinuousParametric(int iident , double ilocation ,
                                           double idispersion , int iunit)

{
  ident = iident;

  location = ilocation;
  dispersion = idispersion;

  min_value = D_DEFAULT;
  max_value = D_DEFAULT;

  unit = iunit;
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

  min_value = dist.min_value;
  max_value = dist.max_value;

  unit = dist.unit;

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
 *              reference sur l'indice de la ligne lue.
 *
 *--------------------------------------------------------------*/

ContinuousParametric* continuous_parametric_parsing(StatError &error , ifstream &in_file ,
                                                    int &line)

{
  RWLocaleSnapshot locale("en");
  RWCString buffer , token;
  size_t position;
  bool status = true , lstatus;
  register int i , j;
  int ident = I_DEFAULT;
  double location = D_INF , dispersion = D_DEFAULT;
  ContinuousParametric *dist;


  dist = NULL;

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
        for (j = GAUSSIAN;j <= VON_MISES;j++) {
          if (token == STAT_continuous_distribution_word[j]) {
            ident = j;
            break;
          }
        }

        if (j == VON_MISES + 1) {
          status = false;
          error.update(STAT_parsing[STATP_DISTRIBUTION_NAME] , line , i + 1);
        }
      }

      // test nom du parametre

      else {
        switch ((i - 1) % 3) {

        case 0 : {
          switch ((i - 1) / 3) {

          // 1er parametre : moyenne

          case 0 : {
            if (token != STAT_word[STATW_MEAN]) {
              status = false;
              error.correction_update(STAT_parsing[STATP_PARAMETER_NAME] , STAT_word[STATW_MEAN] , line , i + 1);
            }
            break;
          }

          // 2eme parametre : ecart-type (Gauss) ou concentration (von Mises)

          case 1 : {
            if ((ident == GAUSSIAN) && (token != STAT_word[STATW_STANDARD_DEVIATION])) {
              status = false;
              error.correction_update(STAT_parsing[STATP_PARAMETER_NAME] , STAT_word[STATW_STANDARD_DEVIATION] , line , i + 1);
            }

            if ((ident == VON_MISES) && (token != STAT_word[STATW_CONCENTRATION])) {
              status = false;
              error.correction_update(STAT_parsing[STATP_PARAMETER_NAME] , STAT_word[STATW_CONCENTRATION] , line , i + 1);
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

          // 1er parametre : moyenne

          case 0 : {
            lstatus = locale.stringToNum(token , &location);
            break;
          }

          // 2eme parametre : ecart-type (Gauss) ou concentration (von Mises)

          case 1 : {
            lstatus = locale.stringToNum(token , &dispersion);
            if ((lstatus) && (dispersion <= 0.)) {
              lstatus = false;
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
      if (i != 7) {
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
    dist = new ContinuousParametric(ident , location , dispersion ,
                                    (location < 2 * M_PI ? RADIAN : DEGREE));
  }

  return dist;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture des parametres d'une loi continue.
 *
 *  argument : stream.
 *
 *--------------------------------------------------------------*/

ostream& ContinuousParametric::ascii_parameter_print(ostream &os) const

{
  os << STAT_continuous_distribution_word[ident] << "   "
     << STAT_word[STATW_MEAN] << " : " << location << "   ";

  switch (ident) {
  case GAUSSIAN :
    os << STAT_word[STATW_STANDARD_DEVIATION] << " : " << dispersion;
    break;
  case VON_MISES :
    os << STAT_word[STATW_CONCENTRATION] << " : " << dispersion;
    break;
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
  double variance = variance_computation();


  if (file_flag) {
    os << "# ";
  }
/*  if ((ident != GAUSSIAN) && (ident != VON_MISES)) {
    os << STAT_label[STATL_MEAN] << ": " << mean << "   ";
  } */
  os << STAT_label[STATL_VARIANCE] << ": " << variance;
  if (ident != GAUSSIAN) {
    os << "   " << STAT_label[STATL_STANDARD_DEVIATION] << ": " << sqrt(variance);
  }
  os << endl;

/*  if ((ident != GAUSSIAN) && (ident != VON_MISES) && (variance > 0.)) {
    if (file_flag) {
      os << "# ";
    }
    os << STAT_label[STATL_SKEWNESS_COEFF] << ": " << skewness_computation() << "   "
       << STAT_label[STATL_KURTOSIS_COEFF] << ": " << kurtosis_computation() << endl;
  } */

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
      frequency = new double[nb_step];
      dist_cumul = new double[nb_step];

      value = min_value;
      dist_cumul[0] = cdf(dist , value + step / 2);
      frequency[0] = (dist_cumul[0] - cdf(dist , value - step / 2)) * scale;
      value += step;

      for (i = 1;i < nb_step;i++) {
        dist_cumul[i] = cdf(dist , value + step / 2);
        frequency[i] = (dist_cumul[i] - dist_cumul[i - 1]) * scale;
        value += step;
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
      value += step;

      for (i = 1;i < nb_step;i++) {
        mass = von_mises_mass_computation(value - step / 2 , value + step / 2);
        frequency[i] = mass * scale;
        dist_cumul[i] = dist_cumul[i - 1] + mass;
        value += step;
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
  os << "\n" << STAT_continuous_distribution_word[ident] << "\t"
     << STAT_word[STATW_MEAN] << "\t" << location << "\t";

  switch (ident) {
  case GAUSSIAN :
    os << STAT_word[STATW_STANDARD_DEVIATION] << "\t" << dispersion;
    break;
  case VON_MISES :
    os << STAT_word[STATW_CONCENTRATION] << "\t" << dispersion;
    break;
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
  double variance = variance_computation();


/*  if ((ident != GAUSSIAN) && (ident != VON_MISES)) {
    os << STAT_label[STATL_MEAN] << "\t" << mean_computation << "\t";
  } */
  os << STAT_label[STATL_VARIANCE] << "\t" << variance;
  if (ident != GAUSSIAN) {
    os << "\t" << STAT_label[STATL_STANDARD_DEVIATION] << "\t" << sqrt(variance);
  }
  os << endl;

/*  if ((ident != GAUSSIAN) && (ident != VON_MISES) && (variance > 0.)) {
    os << STAT_label[STATL_SKEWNESS_COEFF] << "\t" << skewness_computation() << "\t"
       << STAT_label[STATL_KURTOSIS_COEFF] << "\t" << kurtosis_computation() << endl;
  } */

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
      frequency = new double[nb_step];
      dist_cumul = new double[nb_step];

      value = min_value;
      dist_cumul[0] = cdf(dist , value + step / 2);
      frequency[0] = (dist_cumul[0] - cdf(dist , value - step / 2)) * scale;
      value += step;

      for (i = 1;i < nb_step;i++) {
        dist_cumul[i] = cdf(dist , value + step / 2);
        frequency[i] = (dist_cumul[i] - dist_cumul[i - 1]) * scale;
        value += step;
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
      value += step;

      for (i = 1;i < nb_step;i++) {
        mass = von_mises_mass_computation(value - step / 2 , value + step / 2);
        frequency[i] = mass * scale;
        dist_cumul[i] = dist_cumul[i - 1] + mass;
        value += step;
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

  if (location != D_INF) {
    os << location;
  }
  if (dispersion != D_DEFAULT) {
    os << ", " << dispersion;
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
        frequency[i] = pdf(dist , value) *  scale;
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
        max_value = 360.;
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
      max_value = 360.;
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
  case GAUSSIAN :
    nb_parameter = 2;
    break;
  case VON_MISES :
    nb_parameter = 2;
    break;
  default :
    nb_parameter = 0;
    break;
  }

  return nb_parameter;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la variance d'une loi continue.
 *
 *--------------------------------------------------------------*/

double ContinuousParametric::variance_computation() const

{
  double variance;


  switch (ident) {

  case GAUSSIAN : {
    variance = dispersion * dispersion;
    break;
  }

  case VON_MISES : {
    switch (unit) {
    case DEGREE :
      variance = 180 * 180 / (dispersion * M_PI * M_PI);
      break;
    case RADIAN :
      variance = 1. / dispersion;
      break;
    }
    break;
  }
  }

  return variance;
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


//  if (dispersion < CONCENTRATION_THRESHOLD) {
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

  case GAUSSIAN : {
    normal dist(location , dispersion);

    mass = cdf(dist , sup) - cdf(dist , inf);
    break;
  }

  case VON_MISES : {
    mass = von_mises_mass_computation(inf , sup);
    break;
  }
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
  register int i;
  double step , value;


  cumul = new double[VON_MISES_NB_STEP];
  value = step / 2;

  switch (unit) {

  case DEGREE : {
    step = 360. / VON_MISES_NB_STEP;

    cumul[0] = 0.;
    for (i = 1;i < VON_MISES_NB_STEP;i++) {
      cumul[i] = cumul[i - 1] + exp(dispersion * cos((value - location) * M_PI / 180)) * step /
                 (360 * cyl_bessel_i(0 , dispersion));
      value += step;
    }
    break;
  }

  case RADIAN : {
    step = 2 * M_PI / VON_MISES_NB_STEP;

    cumul[0] = 0.;
    for (i = 1;i < VON_MISES_NB_STEP;i++) {
      cumul[i] = cumul[i - 1] + exp(dispersion * cos(value - location)) * step /
                 (2 * M_PI * cyl_bessel_i(0 , dispersion));
      value += step;
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
  register int i;
  double step , value , current_cumul , previous_cumul , **qqplot;


  qqplot = new double*[2];
  qqplot[0] = new double[nb_value];
  qqplot[1] = new double[nb_value];

  for (i = 0;i < nb_value;i++) {
    qqplot[0][i] = empirical_cdf[0][i];
  }

  switch (ident) {

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
 *  Simulation d'une loi continue.
 *
 *--------------------------------------------------------------*/

double ContinuousParametric::simulation() const

{
  register int i;
  double limit , value , step , current_cumul , previous_cumul;


  limit = ((double)rand() / (RAND_MAX + 1.));

  switch (ident) {

  case GAUSSIAN : {
    normal dist(location , dispersion);

    value = quantile(dist , limit);
    break;
  }

  case VON_MISES : {
    if (cumul) {
      value = 0.;
      i = 0;

      switch (unit) {
      case DEGREE :
        step = 360. / VON_MISES_NB_STEP;
        break;
      case RADIAN :
        step = 2 * M_PI / VON_MISES_NB_STEP;
        break;
      }

      do {
        value += step;
        i++;
      }
      while (cumul[i] < limit);

      value -= step * (cumul[i] - limit) / (cumul[i] - cumul[i - 1]);
    }

    else {
      current_cumul = 0.;
      value = 0.;

      switch (unit) {

      case DEGREE : {
        step = 360. / VON_MISES_NB_STEP;

        do {
          value += step;
          previous_cumul = current_cumul;
          current_cumul += exp(dispersion * cos((value - step / 2 - location) * M_PI / 180)) * step /
                           (360 * cyl_bessel_i(0 , dispersion));
        }
        while (current_cumul < limit);
        break;
      }

      case RADIAN : {
        step = 2 * M_PI / VON_MISES_NB_STEP;

        do {
          value += step;
          previous_cumul = current_cumul;
          current_cumul += exp(dispersion * cos(value - step / 2 - location)) * step /
                           (2 * M_PI * cyl_bessel_i(0 , dispersion));
        }
        while (current_cumul < limit);
        break;
      }
      }

      value -= step * (current_cumul - limit) / (current_cumul - previous_cumul);
    }
    break;
  }
  }

  return value;
}
