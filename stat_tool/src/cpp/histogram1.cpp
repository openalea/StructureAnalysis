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



#include <limits.h>
#include <math.h>
#include <sstream>
#include <iomanip>

#include "tool/config.h"

// #include <rw/vstream.h>
// #include <rw/rwfile.h>
#include "stat_tools.h"
#include "distribution.h"
#include "curves.h"
#include "stat_label.h"

using namespace std;


extern int column_width(int value);
extern int column_width(int nb_value , const double *value , double scale = 1.);
extern char* label(const char *file_name);



/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Histogram.
 *
 *  arguments : nombre d'echantillons, echantillons.
 *
 *--------------------------------------------------------------*/

Histogram::Histogram(int inb_element , int *pelement)

{
  register int i;
  int *pfrequency;


  nb_element = inb_element;

  nb_value = 0;
  for (i = 0;i < nb_element;i++) {
    if (*pelement > nb_value) {
      nb_value = *pelement;
    }
    pelement++;
  }
  pelement -= nb_element;

  nb_value++;
  alloc_nb_value = nb_value;
  frequency = new int[nb_value];

  pfrequency = frequency;
  for (i = 0;i < nb_value;i++) {
    *pfrequency++ = 0;
  }

  for (i = 0;i < nb_element;i++) {
    frequency[*pelement++]++;
  }

  // calcul des caracteristiques de l'histogramme

  offset_computation();
  max_computation();
  mean_computation();
  variance_computation();
}


/*--------------------------------------------------------------*
 *
 *  Translation d'un histogramme.
 *
 *  arguments : reference sur un objet Histogram, parametre de translation.
 *
 *--------------------------------------------------------------*/

void Histogram::shift(const Histogram &histo , int shift_param)

{
  register int i;
  int *pfrequency , *cfrequency;


  // calcul des caracteristiques de l'histogramme

  nb_element = histo.nb_element;
  nb_value = histo.nb_value + shift_param;
  alloc_nb_value = nb_value;
  offset = histo.offset + shift_param;
  max = histo.max;
  mean = histo.mean + shift_param;
  variance = histo.variance;

  // copie des frequences

  frequency = new int[nb_value];

  pfrequency = frequency;
  for (i = 0;i < offset;i++) {
    *pfrequency++ = 0;
  }
  cfrequency = histo.frequency + histo.offset;
  for (i = offset;i < nb_value;i++) {
    *pfrequency++ = *cfrequency++;
  }
}


/*--------------------------------------------------------------*
 *
 *  Regroupement des valeurs d'un histogramme.
 *
 *  arguments : reference sur un objet Histogram, pas pour le regroupement.
 *
 *--------------------------------------------------------------*/

void Histogram::cluster(const Histogram &histo , int step)

{
  register int i;
  int *pfrequency , *cfrequency;


  nb_element = histo.nb_element;
  nb_value = (histo.nb_value - 1) / step + 1;
  alloc_nb_value = nb_value;
  offset = histo.offset / step;

  // cumul des frequences

  frequency = new int[nb_value];

  pfrequency = frequency - 1;
  cfrequency = histo.frequency;

  for (i = 0;i < histo.nb_value;i++) {
    if (i % step == 0) {
      *++pfrequency = *cfrequency++;
    }
    else {
      *pfrequency += *cfrequency++;
    }
  }

  // calcul des caracteristiques de l'histogramme

  max_computation();
  mean_computation();
  variance_computation();
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Histogram.
 *
 *  arguments : reference sur un objet Histogram, type de transformation
 *              ('s' : translation, 'c' : groupement des valeurs),
 *              pas de regroupement ('c') / parametre de translation ('s').
 *
 *--------------------------------------------------------------*/

Histogram::Histogram(const Histogram &histo , char transform , int param)

{
  switch (transform) {
  case 's' :
    shift(histo , param);
    break;
  case 'c' :
    cluster(histo , param);
    break;
  default :
    copy(histo);
    break;
  }
}


/*--------------------------------------------------------------*
 *
 *  Operateur d'egalite de la classe Histogram.
 *
 *  argument : reference sur un objet Histogram.
 *
 *--------------------------------------------------------------*/

bool Histogram::operator==(const Histogram &histo) const

{
  bool status = true;
  register int i;
  int *pfrequency , *cfrequency;


  if ((offset != histo.offset) || (nb_value != histo.nb_value) ||
      (nb_element != histo.nb_element)) {
    status = false;
  }

  else {
    pfrequency = frequency + offset;
    cfrequency = histo.frequency + offset;
    for (i = offset;i < nb_value;i++) {
      if (*pfrequency++ != *cfrequency++) {
        status = false;
        break;
      }
    }
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Translation d'un histogramme.
 *
 *  arguments : reference sur un objet Format_error, parametre de translation.
 *
 *--------------------------------------------------------------*/

Distribution_data* Histogram::shift(Format_error &error , int shift_param) const

{
  Distribution_data *histo;


  error.init();

  if (shift_param < -offset) {
    histo = 0;
    ostringstream correction_message;
    correction_message << STAT_error[STATR_GREATER_THAN] << " " << -offset;
    error.correction_update(STAT_error[STATR_SHIFT_VALUE] , (correction_message.str()).c_str());
  }

  else {
    histo = new Distribution_data(*this , 's' , shift_param);
  }

  return histo;
}


/*--------------------------------------------------------------*
 *
 *  Regroupement des valeurs d'un histogramme.
 *
 *  arguments : reference sur un objet Format_error, pas pour le regroupement.
 *
 *--------------------------------------------------------------*/

Distribution_data* Histogram::cluster(Format_error &error , int step) const

{
  Distribution_data *histo;


  error.init();

  if (step < 1) {
    histo = 0;
    error.update(STAT_error[STATR_CLUSTERING_STEP]);
  }

  else {
    histo = new Distribution_data(*this , 'c' , step);
  }

  return histo;
}


/*--------------------------------------------------------------*
 *
 *  Regroupement des valeurs d'un histogramme en fonction
 *  de l'augmentation de la quantite d'information.
 *
 *  arguments : reference sur un objet Format_error, proportion de la
 *              quantite d'information de l'histogramme initial, stream.
 *
 *--------------------------------------------------------------*/

Distribution_data* Histogram::cluster(Format_error &error , double ratio , ostream &os) const

{
  bool status = true , stop = false;
  register int i;
  int step = 1 , *pfrequency , *cfrequency;
  double information , reference_information , previous_information;
  Distribution_data *histo , *previous_histo;


  previous_histo = 0;
  error.init();

  if ((ratio < 0.) || (ratio > 1.)) {
    status = false;
    error.update(STAT_error[STATR_INFORMATION_RATIO]);
  }
  if (variance == 0.) {
    status = false;
    error.update(STAT_error[STATR_NULL_INFORMATION]);
  }

  if (status) {
    reference_information = information_computation();
    previous_information = reference_information;
    previous_histo = new Distribution_data(*this);
    histo = new Distribution_data((nb_value - 1) / 2 + 1);
    histo->nb_element = nb_element;

    do {
      step++;

      // regroupement des valeurs

      pfrequency = histo->frequency - 1;
      cfrequency = frequency;

      for (i = 0;i < nb_value;i++) {
        if (i % step == 0) {
          *++pfrequency = *cfrequency++;
        }
        else {
          *pfrequency += *cfrequency++;
        }
      }

      histo->offset = offset / step;
      histo->nb_value = (nb_value - 1) / step + 1;

      // calcul de la quantite d'information

      information = histo->information_computation();

#     ifdef DEBUG
      cout << "\n" << STAT_label[STATL_CLUSTERING_STEP] << ": " << step
           << "   " << STAT_label[STATL_INFORMATION] << ": " << information
           << "   " << STAT_label[STATL_INFORMATION_RATIO] << ": "
           << information / reference_information << endl;
#     endif

      if (information / reference_information < ratio) {
        stop = true;

        if (fabs(information / reference_information - ratio) <
            fabs(previous_information / reference_information - ratio)) {
          previous_information = information;
          delete previous_histo;
          previous_histo = new Distribution_data(*histo);
        }
        else {
          step--;
        }
      }

      else {
        previous_information = information;
        delete previous_histo;
        previous_histo = new Distribution_data(*histo);
      }
    }
    while (!stop);

    delete histo;

#   ifdef MESSAGE
    os << STAT_label[STATL_INFORMATION_RATIO] << ": "
       << previous_information / reference_information << "   "
       << STAT_label[STATL_CLUSTERING_STEP] << ": " << step << endl;
#   endif

    // calcul des caracteristiques de l'histogramme

    previous_histo->max_computation();
    previous_histo->mean_computation();
    previous_histo->variance_computation();
  }

  return previous_histo;
}


/*--------------------------------------------------------------*
 *
 *  Regroupement des valeurs d'un histogramme.
 *
 *  arguments : reference sur un objet Format_error, nombres de classes,
 *              bornes pour regrouper les valeurs.
 *
 *--------------------------------------------------------------*/

Distribution_data* Histogram::cluster(Format_error &error , int nb_class , int *ilimit) const

{
  bool status = true;
  register int i , j;
  int *pfrequency , *cfrequency , *limit;
  Distribution_data *histo;


  histo = 0;
  error.init();

  if ((nb_class < 2) || (nb_class >= nb_value)) {
    status = false;
    error.update(STAT_error[STATR_NB_CLASS]);
  }

  else {
    limit = new int[nb_class + 1];
    limit[0] = 0;
    for (i = 1;i < nb_class;i++) {
      limit[i] = ilimit[i - 1];
    }
    limit[nb_class] = nb_value;

    for (i = 0;i < nb_class;i++) {
      if (limit[i] >= limit[i + 1]) {
        status = false;
        error.update(STAT_error[STATR_CLUSTER_LIMIT]);
      }
    }

    if (status) {
      histo = new Distribution_data(nb_class);

      // regroupement des valeurs

      pfrequency = histo->frequency - 1;
      cfrequency = frequency;

      for (i = 0;i < histo->nb_value;i++) {
        *++pfrequency = *cfrequency++;
        for (j = limit[i] + 1;j < limit[i + 1];j++) {
          *pfrequency += *cfrequency++;
        }
      }

      // calcul des caracteristiques de l'histogramme

      histo->offset_computation();
      histo->nb_element = nb_element;
      histo->max_computation();
      histo->mean_computation();
      histo->variance_computation();
    }

    delete [] limit;
  }

  return histo;
}


/*--------------------------------------------------------------*
 *
 *  Transcodage des symboles.
 *
 *  arguments : reference sur un objet Format_error,
 *              table de transcodage des symboles.
 *
 *--------------------------------------------------------------*/

Distribution_data* Histogram::transcode(Format_error &error , int *symbol) const

{
  bool status = true , *presence;
  register int i;
  int min_symbol , max_symbol , *psymbol , *pfrequency , *cfrequency;
  Distribution_data *histo;


  histo = 0;
  error.init();

  min_symbol = INT_MAX;
  max_symbol = 0;

  for (i = 0;i < nb_value - offset;i++) {
    if (symbol[i] < 0) {
      status = false;
      ostringstream error_message;
      error_message << STAT_label[STATL_SYMBOL] << " " << symbol[i] << " "
                    << STAT_error[STATR_NOT_ALLOWED];
      error.update((error_message.str()).c_str());
    }
    else {
      if (symbol[i] < min_symbol) {
        min_symbol = symbol[i];
      }
      if (symbol[i] > max_symbol) {
        max_symbol = symbol[i];
      }
    }
  }

  if (max_symbol - min_symbol == 0) {
    status = false;
    error.update(STAT_error[STATR_NB_SYMBOL]);
  }

  if (max_symbol - min_symbol > nb_value - 1 - offset) {
    status = false;
    error.update(STAT_error[STATR_NON_CONSECUTIVE_SYMBOLS]);
  }

  if (status) {
    presence = new bool[max_symbol + 1];
    for (i = min_symbol;i <= max_symbol;i++) {
      presence[i] = false;
    }

    for (i = 0;i < nb_value - offset;i++) {
      presence[symbol[i]] = true;
    }

    for (i = min_symbol;i <= max_symbol;i++) {
      if (!presence[i]) {
        status = false;
        ostringstream error_message;
        error_message << STAT_error[STATR_MISSING_SYMBOL] << " " << i;
        error.update((error_message.str()).c_str());
      }
    }

    delete [] presence;
  }

  if (status) {
    histo = new Distribution_data(max_symbol + 1);

    // transcodage des symboles

    pfrequency = histo->frequency;
    for (i = 0;i < histo->nb_value;i++) {
      *pfrequency++ = 0;
    }

    psymbol = symbol;
    cfrequency = frequency + offset;
    for (i = 0;i < nb_value - offset;i++) {
      histo->frequency[*psymbol++] += *cfrequency++;
    }

    // calcul des caracteristiques de l'histogramme

    histo->offset = min_symbol;
    histo->nb_element = nb_element;
    histo->max_computation();
    histo->mean_computation();
    histo->variance_computation();
  }

  return histo;
}


/*--------------------------------------------------------------*
 *
 *  Selection d'individus dans une plage de valeurs.
 *
 *  arguments : reference sur un objet Format_error, bornes sur les valeurs,
 *              flag pour conserver ou rejeter les individus selectionnes.
 *
 *--------------------------------------------------------------*/

Distribution_data* Histogram::value_select(Format_error &error , int min_value ,
                                           int max_value , bool keep) const

{
  bool status = true;
  register int i;
  int *pfrequency , *cfrequency;
  Distribution_data *histo;


  histo = 0;
  error.init();

  if ((min_value < 0) || (min_value >= nb_value) || (min_value > max_value)) {
    status = false;
    error.update(STAT_error[STATR_MIN_VALUE]);
  }
  if ((max_value < offset) || (max_value < min_value)) {
    status = false;
    error.update(STAT_error[STATR_MAX_VALUE]);
  }

  if (status) {
    switch (keep) {

    case false : {
      histo = new Distribution_data(nb_value);

      pfrequency = histo->frequency;
      cfrequency = frequency;
      for (i = 0;i < min_value;i++) {
        *pfrequency++ = *cfrequency++;
      }
      for (i = min_value;i <= MIN(max_value , nb_value - 1);i++) {
        *pfrequency++ = 0;
      }
      if (max_value + 1 < nb_value) {
        cfrequency += max_value - min_value + 1;
        for (i = max_value + 1;i < nb_value;i++) {
          *pfrequency++ = *cfrequency++;
        }
      }

      break;
    }

    case true : {
      histo = new Distribution_data(MIN(max_value + 1 , nb_value));

      // copie des valeurs

      pfrequency = histo->frequency;
      for (i = 0;i < min_value;i++) {
        *pfrequency++ = 0;
      }
      cfrequency = frequency + min_value;
      for (i = min_value;i < histo->nb_value;i++) {
        *pfrequency++ = *cfrequency++;
      }
      break;
    }
    }

    // calcul des caracteristiques de l'histogramme

    histo->nb_value_computation();
    histo->offset_computation();
    histo->nb_element_computation();

    if (histo->nb_element > 0) {
      histo->max_computation();
      histo->mean_computation();
      histo->variance_computation();
    }

    else {
      delete histo;
      histo = 0;
      error.update(STAT_error[STATR_EMPTY_HISTOGRAM]);
    }
  }

  return histo;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un histogramme.
 *
 *  arguments : stream, flag commentaire,
 *              flag sur l'ecriture de la fonction de repartition.
 *
 *--------------------------------------------------------------*/

ostream& Histogram::ascii_print(ostream &os , int comment_flag , bool cumul_flag) const

{
  register int i;
  int width[3];
  long old_adjust;
  double *cumul;


  old_adjust = os.setf(ios::right , ios::adjustfield);

  // calcul des largeurs des colonnes

  width[0] = column_width(nb_value - 1);
  width[1] = column_width(max) + ASCII_SPACE;

  if (cumul_flag) {
    cumul = cumul_computation();
    width[2] = column_width(nb_value , cumul) + ASCII_SPACE;
  }

  // ecriture des frequences et de la fonction de repartition

  for (i = 0;i < nb_value;i++) {
    if (comment_flag == 1) {
      os << "# ";
    }
    os << setw(width[0]) << i;
    os << setw(width[1]) << frequency[i];

    if (cumul_flag) {
      if (comment_flag == 2) {
        os << "  #";
      }
      os << setw(width[2]) << cumul[i];
    }
    os << endl;
  }

  if (cumul_flag) {
    delete [] cumul;
  }

  os.setf((FMTFLAGS)old_adjust , ios::adjustfield);

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Histogram.
 *
 *  arguments : stream, flag niveau de detail, flag fichier.
 *
 *--------------------------------------------------------------*/

ostream& Histogram::ascii_write(ostream &os , bool exhaustive , bool file_flag) const

{
  double information = information_computation();


  if (exhaustive && file_flag) {
    os << "# ";
  }
  os << STAT_label[STATL_HISTOGRAM] << " - ";
  ascii_characteristic_print(os , true , exhaustive && file_flag);

  if (exhaustive && file_flag) {
    os << "# ";
  }
  os << STAT_label[STATL_MEAN_ABSOLUTE_DEVIATION] << ": " << mean_absolute_deviation_computation();
  if (mean > 0.) {
    os << "   " << STAT_label[STATL_CONCENTRATION_COEFF] << ": " << concentration_computation();
  }
  os << endl;

  if (exhaustive && file_flag) {
    os << "# ";
  }
  os << STAT_label[STATL_INFORMATION] << ": " << information << " ("
     << information / nb_element << ")" << endl;

  if (exhaustive) {
    os << "\n";
    if (file_flag) {
      os << "# ";
    }
    os << "   | " << STAT_label[STATL_HISTOGRAM] << " | " << STAT_label[STATL_CUMULATIVE]
       << " " << STAT_label[STATL_HISTOGRAM] << " " << STAT_label[STATL_FUNCTION] << endl;
    ascii_print(os , (file_flag ? 2 : 0) , true);
  }

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Histogram dans un fichier.
 *
 *  arguments : reference sur un objet Format_error, path.
 *
 *--------------------------------------------------------------*/

bool Histogram::ascii_write(Format_error &error , const char *path) const

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
    ascii_write(out_file , true , true);
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture des caracteristiques d'un histogramme au format tableur.
 *
 *  arguments : stream, flag ecriture des parametres de forme.
 *
 *--------------------------------------------------------------*/

ostream& Histogram::spreadsheet_characteristic_print(ostream &os , bool shape) const

{
  os << STAT_label[STATL_SAMPLE_SIZE] << "\t" << nb_element << endl;

  if ((mean != D_DEFAULT) && (variance != D_DEFAULT)) {
    os << STAT_label[STATL_MEAN] << "\t" << mean << "\t\t"
       << STAT_label[STATL_VARIANCE] << "\t" << variance << "\t\t"
       << STAT_label[STATL_STANDARD_DEVIATION] << "\t" << sqrt(variance) << endl;

    if ((shape) && (variance > 0.)) {
      os << STAT_label[STATL_SKEWNESS_COEFF] << "\t" << skewness_computation() << "\t\t"
         << STAT_label[STATL_KURTOSIS_COEFF] << "\t" << kurtosis_computation() << endl;
    }
  }

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un histogramme au format tableur.
 *
 *  arguments : stream, flags sur l'ecriture de la fonction de repartition et
 *              de la fonction de concentration.
 *
 *--------------------------------------------------------------*/

ostream& Histogram::spreadsheet_print(ostream &os , bool cumul_flag , bool concentration_flag) const

{
  register int i;
  double *cumul, *concentration;


  if ((!cumul_flag) || (mean == 0.)) {
    concentration_flag = false;
  }

  if (cumul_flag) {
    cumul = cumul_computation();
  }
  if (concentration_flag) {
    concentration = concentration_function_computation();
  }

  for (i = 0;i < nb_value;i++) {
    os << i << "\t" << frequency[i];
    if (cumul_flag) {
      os << "\t" << cumul[i];
    }
    if (concentration_flag) {
      os << "\t" << concentration[i];
    }
    os << endl;
  }

  if (cumul_flag) {
    delete [] cumul;
  }
  if (concentration_flag) {
    delete [] concentration;
  }

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un histogramme au format Gnuplot.
 *
 *  arguments : path, pointeurs sur les fonctions de repartition et
 *              de concentration, decalage de l'histogramme.
 *
 *--------------------------------------------------------------*/

bool Histogram::plot_print(const char *path , double *cumul ,
                           double *concentration , double shift) const

{
  bool status = false;
  register int i;
  ofstream out_file(path);


  if (out_file) {
    status = true;

    if ((offset == 0) && (variance > 0.)) {
      out_file << -1 << " " << -1 + shift << " " << 0 << " "
               << 0 << " " << 0 << " " << 0 << endl;
    }

    for (i = 0;i < nb_value;i++) {
      out_file << i << " " << i + shift << " " << frequency[i] << " "
               << (double)frequency[i] / (double)nb_element;
      if (variance > 0.) {
        out_file << " " << cumul[i] << " " << concentration[i];
      }
      out_file << endl;
    }

    if (nb_value == 1) {
      out_file << nb_value << " " << nb_value + shift << " "
               << 0 << " " << 0 << endl;
    }
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'une famille histogrammes au format Gnuplot.
 *
 *  arguments : path, nombre d'histogrammes,
 *              pointeurs sur les histogrammes.
 *
 *--------------------------------------------------------------*/

bool Histogram::plot_print(const char *path , int nb_histo ,
                           const Histogram **histo) const

{
  bool status = false;
  register int i , j;
  int plot_nb_value = nb_value;
  ofstream out_file(path);


  if (out_file) {
    status = true;

    for (i = 0;i < nb_histo;i++) {
      if (histo[i]->nb_value > plot_nb_value) {
        plot_nb_value = histo[i]->nb_value;
      }
    }
    if (plot_nb_value < 2) {
      plot_nb_value = 2;
    }

    for (i = 0;i < plot_nb_value;i++) {
      if (i < nb_value) {
        out_file << frequency[i];
      }
      else {
        out_file << 0;
      }

      for (j = 0;j < nb_histo;j++) {
        if (i < histo[j]->nb_value) {
          out_file << " " << histo[j]->frequency[i];
        }
        else {
          out_file << " " << 0;
        }
      }

      out_file << endl;
    }
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Calcul des taux de survie a partir d'un histogramme et
 *  ecriture du resultat.
 *
 *  argument : stream.
 *
 *--------------------------------------------------------------*/

ostream& Histogram::survival_ascii_write(ostream &os) const

{
  Curves *survival_rate;


  os << STAT_label[STATL_HISTOGRAM] << " - ";
  ascii_characteristic_print(os);

  os << "\n   | " << STAT_label[STATL_HISTOGRAM] << " | " << STAT_label[STATL_CUMULATIVE]
     << " " << STAT_label[STATL_HISTOGRAM] << " " << STAT_label[STATL_FUNCTION] << endl;
  ascii_print(os , false , true);

  survival_rate = new Curves(*this);

  os << "\n   | " << STAT_label[STATL_DEATH_PROBABILITY] << " | "
     << STAT_label[STATL_SURVIVAL_PROBABILITY] << " | " << STAT_label[STATL_FREQUENCY] << endl;
  survival_rate->ascii_print(os);

  delete survival_rate;

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Calcul des taux de survie a partir d'un histogramme et
 *  ecriture du resultat dans un fichier.
 *
 *  arguments : reference sur un objet Format_error, path.
 *
 *--------------------------------------------------------------*/

bool Histogram::survival_ascii_write(Format_error &error , const char *path) const

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
    survival_ascii_write(out_file);
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Calcul des taux de survie a partir d'un histogramme et ecriture
 *  du resultat dans un fichier au format tableur.
 *
 *  arguments : reference sur un objet Format_error, path.
 *
 *--------------------------------------------------------------*/

bool Histogram::survival_spreadsheet_write(Format_error &error , const char *path) const

{
  bool status;
  Curves *survival_rate;
  ofstream out_file(path);


  error.init();

  if (!out_file) {
    status = false;
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {
    status = true;

    out_file << STAT_label[STATL_HISTOGRAM] << "\t";
    spreadsheet_characteristic_print(out_file);

    out_file << "\n\t" << STAT_label[STATL_HISTOGRAM] << "\t" << STAT_label[STATL_CUMULATIVE]
             << " " << STAT_label[STATL_HISTOGRAM] << " " << STAT_label[STATL_FUNCTION] << endl;
    spreadsheet_print(out_file , true);

    survival_rate = new Curves(*this);

    out_file << "\n\t" << STAT_label[STATL_DEATH_PROBABILITY] << "\t"
             << STAT_label[STATL_SURVIVAL_PROBABILITY] << "\t" << STAT_label[STATL_FREQUENCY] << endl;
    survival_rate->spreadsheet_print(out_file);

    delete survival_rate;
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un histogramme, de la loi et de la fonction de survie deduites
 *  au format Gnuplot.
 *
 *  arguments : path, pointeur sur la fonction de survie.
 *
 *--------------------------------------------------------------*/

bool Histogram::survival_plot_print(const char *path , double *survivor) const

{
  bool status = false;
  register int i;
  ofstream out_file(path);


  if (out_file) {
    status = true;

    for (i = 0;i < nb_value;i++) {
      out_file << frequency[i] << " " << (double)frequency[i] / (double)nb_element << " "
               << survivor[i] << endl;
    }
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Calcul des taux de survie a partir d'un histogramme et
 *  sortie Gnuplot du resultat.
 *
 *  arguments : reference sur un objet Format_error, prefixe des fichiers,
 *              titre des figures.
 *
 *--------------------------------------------------------------*/

bool Histogram::survival_plot_write(Format_error &error , const char *prefix ,
                                    const char *title) const

{
  bool status;
  register int i;
  double *survivor;
  Curves *survival_rate;
  ostringstream data_file_name[2];


  error.init();

  if (variance == 0.) {
    status = false;
    error.update(STAT_error[STATR_PLOT_NULL_VARIANCE]);
  }

  else {

    // ecriture des fichiers de donnees

    data_file_name[0] << prefix << 0 << ".dat";
    survivor = survivor_function_computation();

    status = survival_plot_print((data_file_name[0].str()).c_str() , survivor);

    delete [] survivor;

    if (!status) {
      error.update(STAT_error[STATR_FILE_PREFIX]);
    }

    else {
      survival_rate = new Curves(*this);

      data_file_name[1] << prefix << 1 << ".dat";
      survival_rate->plot_print((data_file_name[1].str()).c_str());

      // ecriture du fichier de commandes et du fichier d'impression

      for (i = 0;i < 2;i++) {
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

        // histogramme

        if (nb_value - 1 < TIC_THRESHOLD) {
          out_file << "set xtics 0,1" << endl;
        }
        if ((int)(max * YSCALE) + 1 < TIC_THRESHOLD) {
          out_file << "set ytics 0,1" << endl;
        }

        out_file << "plot [0:" << nb_value - 1 << "] [0:"
                 << (int)(max * YSCALE) + 1 << "] \""
                 << label((data_file_name[0].str()).c_str()) << "\" using 1 title \""
                 << STAT_label[STATL_HISTOGRAM] << "\" with impulses" << endl;

        if ((int)(max * YSCALE) + 1 < TIC_THRESHOLD) {
          out_file << "set ytics autofreq" << endl;
        }

        // loi et fonction de survie

        if (i == 0) {
          out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
        }
        out_file << endl;

        out_file << "plot [0:" << nb_value - 1 << "] [0:1] \""
                 << label((data_file_name[0].str()).c_str()) << "\" using 2 title \""
                 << STAT_label[STATL_DISTRIBUTION] << "\" with linespoints,\\" << endl;
        out_file << "\"" << label((data_file_name[0].str()).c_str()) << "\" using 3 title \""
                 << STAT_label[STATL_SURVIVOR] << " " << STAT_label[STATL_FUNCTION]
                 << "\" with linespoints" << endl;

        if (nb_value - 1 < TIC_THRESHOLD) {
          out_file << "set xtics autofreq" << endl;
        }

        // taux de survie

        if (i == 0) {
          out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
        }
        out_file << endl;

        if (survival_rate->length - 1 < TIC_THRESHOLD) {
          out_file << "set xtics 0,1" << endl;
        }

        out_file << "plot [" << survival_rate->offset << ":" << survival_rate->length - 1
                 << "] [0:1] \"" << label((data_file_name[1].str()).c_str()) << "\" using 1 title \""
                 << STAT_label[STATL_DEATH_PROBABILITY] << "\" with linespoints,\\" << endl;
        out_file << "\"" << label((data_file_name[1].str()).c_str()) << "\" using 2 title \""
                 << STAT_label[STATL_SURVIVAL_PROBABILITY] << "\" with linespoints" << endl;

        if (survival_rate->length - 1 < TIC_THRESHOLD) {
          out_file << "set xtics autofreq" << endl;
        }

        if (i == 1) {
          out_file << "\nset terminal x11" << endl;
        }

        out_file << "\npause 0 \"" << STAT_label[STATL_END] << "\"" << endl;
      }

      delete survival_rate;
    }
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un histogramme.
 *
 *  argument : reference sur un objet SinglePlot.
 *
 *--------------------------------------------------------------*/

void Histogram::plotable_frequency_write(SinglePlot &plot) const

{
  register int i;

  for (i = offset;i < nb_value;i++) 
    plot.add_point(i, frequency[i]);
  
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un histogramme.
 *
 *  arguments : reference sur un objet SinglePlot,
 *              facteur d'echelle (valeur par defaut : 1).
 *
 *--------------------------------------------------------------*/

void Histogram::plotable_cumul_write(SinglePlot &plot , double scale) const

{
  register int i;
  double *cumul;

  cumul = cumul_computation(scale);

  for (i = MAX(offset - 1 , 0); i < nb_value; i++) 
    plot.add_point(i, cumul[i]);


  delete [] cumul;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la fonction de repartition deduite d'un histogramme.
 *
 *  argument : facteur d'echelle.
 *
 *--------------------------------------------------------------*/

double* Histogram::cumul_computation(double scale) const

{
  register int i;
  int *pfrequency;
  double *cumul , *pcumul;


  if (scale == D_DEFAULT) {
    scale = nb_element;
  }

  cumul = new double[nb_value];

  pcumul = cumul;
  for (i = 0;i < offset;i++) {
    *pcumul++ = 0.;
  }

  pfrequency = frequency + offset;

  *pcumul++ = *pfrequency++ / scale;
  for (i = offset + 1;i < nb_value;i++) {
    *pcumul = *(pcumul - 1) + *pfrequency++ / scale;
    pcumul++;
  }

  return cumul;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la fonction de survie deduite d'un histogramme.
 *
 *  argument : facteur d'echelle.
 *
 *--------------------------------------------------------------*/

double* Histogram::survivor_function_computation(double scale) const

{
  register int i;
  int *pfrequency;
  double *survivor_function , *psurvivor;


  if (scale == D_DEFAULT) {
    scale = nb_element;
  }

  survivor_function = new double[nb_value];

  psurvivor = survivor_function + nb_value;
  pfrequency = frequency + nb_value;

  *--psurvivor = 0.;
  for (i = nb_value - 2;i >= MAX(offset - 1 , 0);i--) {
    psurvivor--;
    *psurvivor = *(psurvivor + 1) + *--pfrequency / scale;
  }

  for (i = offset - 2;i >= 0;i--) {
    psurvivor--;
    *psurvivor = *(psurvivor + 1);
  }

  return survivor_function;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la fonction de concentration deduite d'un histogramme.
 *
 *  argument : facteur d'echelle.
 *
 *--------------------------------------------------------------*/

double* Histogram::concentration_function_computation(double scale) const

{
  register int i;
  int *pfrequency;
  double norm , *concentration_function , *pconcentration;


  if (mean > 0.) {
    if (scale == D_DEFAULT) {
      scale = nb_element;
    }

    concentration_function = new double[nb_value];

    pconcentration = concentration_function;
    for (i = 0;i < offset;i++) {
      *pconcentration++ = 0.;
    }

    pfrequency = frequency + offset;
    norm = mean * scale;

    *pconcentration++ = *pfrequency++ * offset / norm;
    for (i = offset + 1;i < nb_value;i++) {
      *pconcentration = *(pconcentration - 1) + *pfrequency++ * i / norm;
      pconcentration++;
    }
  }

  else {
    concentration_function = 0;
  }

  return concentration_function;
}


/*--------------------------------------------------------------*
 *
 *  Calcul du coefficient de concentration d'un histogramme.
 *
 *--------------------------------------------------------------*/

double Histogram::concentration_computation() const

{
  register int i;
  int *pfrequency;
  double concentration = D_DEFAULT , *concentration_function;


  if (mean > 0.) {
    concentration_function = concentration_function_computation();
    pfrequency = frequency + offset;

    concentration = *pfrequency++ * concentration_function[offset];
    for (i = offset + 1;i < nb_value;i++) {
      concentration += *pfrequency++ * (concentration_function[i - 1] + concentration_function[i]);
    }

    concentration = 1. - concentration / nb_element;

    delete [] concentration_function;

#   ifdef DEBUG
    int previous_value , cumul;
    double concentration2 = 0.;

    pfrequency = frequency + offset;
    cumul = *pfrequency++;
    previous_value = offset;

    for (i = offset + 1;i < nb_value;i++) {
      if (*pfrequency > 0) {
        concentration2 += cumul * (double)(nb_element - cumul) * (i - previous_value);
        cumul += *pfrequency;
        previous_value = i;
      }
      pfrequency++;
    }

    concentration2 /= (nb_element * (double)nb_element * mean);

    cout << "\n" << STAT_label[STATL_CONCENTRATION_COEFF] << ": " << concentration
         << " | " << concentration2 << endl;
#   endif

  }

  return concentration;
}


/*--------------------------------------------------------------*
 *
 *  Mise a jour d'un histogramme a frequences entieres a partir
 *  d'un histogramme a frequences reelles par arrondi.
 *
 *  arguments : pointeur sur un histogramme a frequences reelles,
 *              effectif de l'histogramme.
 *
 *--------------------------------------------------------------*/

void Histogram::update(const Reestimation<double> *reestim , int inb_element)

{
  register int i , j;
  int index , *pfrequency;
  double scale , sum , max_frequency , *real_frequency , *rfrequency , *cfrequency;


  // copie de l'histogramme reel et mise a l'echelle

  real_frequency = new double[reestim->nb_value];

  rfrequency = real_frequency + reestim->offset;
  cfrequency = reestim->frequency + reestim->offset;
  scale = inb_element / reestim->nb_element;
  for (i = reestim->offset;i < reestim->nb_value;i++) {
    *rfrequency++ = *cfrequency++ * scale;
  }

  // calcul des frequences

  pfrequency = frequency;
  for (i = 0;i < reestim->offset;i++) {
    *pfrequency++ = 0;
  }

  rfrequency = real_frequency + reestim->offset;
  sum = 0.;

  for (i = reestim->offset;i < reestim->nb_value;i++) {
    *pfrequency = (int)*rfrequency;
    *rfrequency -= *pfrequency++;
    if (*rfrequency > 0.) {
      sum += *rfrequency;
    }
    rfrequency++;
  }

  for (i = reestim->nb_value;i < alloc_nb_value;i++) {
    *pfrequency++ = 0;
  }

  // prise en compte des arrondis

  for (i = 0;i < (int)round(sum);i++) {
    rfrequency = real_frequency + reestim->offset;
    max_frequency = 0.;

    for (j = reestim->offset;j < reestim->nb_value;j++) {
      if (*rfrequency > max_frequency) {
        max_frequency = *rfrequency;
        index = j;
      }
      rfrequency++;
    }

    real_frequency[index] = 0.;
    frequency[index]++;
  }

  delete [] real_frequency;

  // calcul des caracteristiques de l'histogramme

  nb_value_computation();
  offset_computation();
  nb_element_computation();
  max_computation();
  mean_computation();
  variance_computation();
}


/*--------------------------------------------------------------*
 *
 *  Calcul d'un histogramme a partir d'un histogramme initial
 *  en changeant l'effectif.
 *
 *  argument : effectif total.
 *
 *--------------------------------------------------------------*/

Histogram* Histogram::frequency_scale(int inb_element) const

{
  register int i , j;
  int index , *pfrequency , *cfrequency;
  double sum , real_max , *real_frequency , *rfrequency;
  Histogram *histo;


  real_frequency = new double[nb_value];
  histo = new Histogram(nb_value);

  cfrequency = frequency + offset;
  rfrequency = real_frequency + offset;
  pfrequency = histo->frequency + offset;
  sum = 0.;

  for (i = offset;i < nb_value;i++) {
    if (*cfrequency > 0) {
      *rfrequency = *cfrequency * (double)inb_element / (double)nb_element;
      *pfrequency = (int)*rfrequency;
      *rfrequency -= *pfrequency;
      if (*rfrequency > 0.) {
        sum += *rfrequency;
      }
    }

    else {
      *rfrequency = 0.;
    }

    cfrequency++;
    rfrequency++;
    pfrequency++;
  }

  // prise en compte des arrondis

  for (i = 0;i < (int)round(sum);i++) {
    rfrequency = real_frequency + offset;
    real_max = 0.;

    for (j = offset;j < nb_value;j++) {
      if (*rfrequency > real_max) {
        real_max = *rfrequency;
        index = j;
      }
      rfrequency++;
    }

    real_frequency[index] = 0.;
    (histo->frequency[index])++;
  }

  delete [] real_frequency;

# ifdef DEBUG
  cout << "\n" << STAT_label[STATL_SAMPLE_SIZE] << " : " << inb_element << " | ";
  histo->nb_element_computation();
  cout << histo->nb_element << endl;
# endif

  // calcul des caracteristiques de l'histogramme

  histo->nb_value_computation();
  histo->offset_computation();
  histo->nb_element = inb_element;
  histo->max_computation();
  histo->mean_computation();
  histo->variance_computation();

  return histo;
}


/*--------------------------------------------------------------*
 *
 *  Calcul des rangs a partir d'un histogramme.
 *
 *--------------------------------------------------------------*/

double* Histogram::rank_computation() const

{
  register int i;
  double *rank;


  rank = new double[nb_value];

  rank[offset] = (double)(1 + frequency[offset]) / 2.;
  for (i = offset + 1;i < nb_value;i++) {
    rank[i] = rank[i - 1] + (double)(frequency[i - 1] + frequency[i]) / 2.;
  }

  return rank;
}
