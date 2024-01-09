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
 *       $Id: frequency_distribution1.cpp 18001 2015-04-23 06:56:57Z guedon $
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



#include <limits.h>
#include <math.h>
#include <sstream>
#include <iomanip>

#include "tool/config.h"

#include "stat_tools.h"
#include "distribution.h"
#include "curves.h"
#include "stat_label.h"

using namespace std;


namespace stat_tool {



/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe FrequencyDistribution.
 *
 *  arguments : nombre d'individus, individus.
 *
 *--------------------------------------------------------------*/

FrequencyDistribution::FrequencyDistribution(int inb_element , int *pelement)

{
  register int i;


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

  for (i = 0;i < nb_value;i++) {
    frequency[i] = 0;
  }
  for (i = 0;i < nb_element;i++) {
    frequency[*pelement++]++;
  }

  // calcul des caracteristiques de la loi empirique

  offset_computation();
  max_computation();
  mean_computation();
  variance_computation();
}


/*--------------------------------------------------------------*
 *
 *  Translation d'une loi empirique.
 *
 *  arguments : reference sur un objet FrequencyDistribution, parametre de translation.
 *
 *--------------------------------------------------------------*/

void FrequencyDistribution::shift(const FrequencyDistribution &histo , int shift_param)

{
  register int i;
  int *cfrequency;


  // calcul des caracteristiques de la loi empirique

  nb_element = histo.nb_element;
  nb_value = histo.nb_value + shift_param;
  alloc_nb_value = nb_value;
  offset = histo.offset + shift_param;
  max = histo.max;
  mean = histo.mean + shift_param;
  variance = histo.variance;

  // copie des frequences

  frequency = new int[nb_value];

  for (i = 0;i < offset;i++) {
    frequency[i] = 0;
  }
  cfrequency = histo.frequency + histo.offset;
  for (i = offset;i < nb_value;i++) {
    frequency[i] = *cfrequency++;
  }
}


/*--------------------------------------------------------------*
 *
 *  Regroupement des valeurs d'une loi empirique.
 *
 *  arguments : reference sur un objet FrequencyDistribution, pas pour le regroupement,
 *              mode (FLOOR/ROUND/CEIL).
 *
 *--------------------------------------------------------------*/

void FrequencyDistribution::cluster(const FrequencyDistribution &histo ,
                                    int step , int mode)

{
  register int i;


  nb_element = histo.nb_element;

  switch (mode) {
  case FLOOR :
    offset = histo.offset / step;
    nb_value = (histo.nb_value - 1) / step + 1;
    break;
  case ROUND :
    offset = (histo.offset + step / 2) / step;
    nb_value = (histo.nb_value - 1 + step / 2) / step + 1;
    break;
  case CEIL :
    offset = (histo.offset + step - 1) / step;
    nb_value = (histo.nb_value + step - 2) / step + 1;
    break;
  }

  alloc_nb_value = nb_value;

  // cumul des frequences

  frequency = new int[nb_value];

  for (i = 0;i < nb_value;i++) {
    frequency[i] = 0;
  }

  switch (mode) {

  case FLOOR : {
    for (i = histo.offset;i < histo.nb_value;i++) {
      frequency[i / step] += histo.frequency[i];
    }
    break;
  }

  case ROUND : {
    for (i = histo.offset;i < histo.nb_value;i++) {
      frequency[(i + step / 2) / step] += histo.frequency[i];
//      frequency[(int)round((double)i / (double)step)] += histo.frequency[i];
    }
    break;
  }

  case CEIL : {
    for (i = histo.offset;i < histo.nb_value;i++) {
      frequency[(i + step - 1) / step] += histo.frequency[i];
//      frequency[(int)ceil((double)i / (double)step)] += histo.frequency[i];
    }
    break;
  }
  }

  // calcul des caracteristiques de la loi empirique

  max_computation();
  mean_computation();
  variance_computation();
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe FrequencyDistribution.
 *
 *  arguments : reference sur un objet FrequencyDistribution, type de transformation
 *              ('s' : translation, 'c' : groupement des valeurs),
 *              pas de regroupement ('c') / parametre de translation ('s'),
 *              mode regroupement (FLOOR/ROUND/CEIL).
 *
 *--------------------------------------------------------------*/

FrequencyDistribution::FrequencyDistribution(const FrequencyDistribution &histo ,
                                             char transform , int param , int mode)

{
  switch (transform) {
  case 's' :
    shift(histo , param);
    break;
  case 'c' :
    cluster(histo , param , mode);
    break;
  default :
    copy(histo);
    break;
  }
}


/*--------------------------------------------------------------*
 *
 *  Operateur d'egalite de la classe FrequencyDistribution.
 *
 *  argument : reference sur un objet FrequencyDistribution.
 *
 *--------------------------------------------------------------*/

bool FrequencyDistribution::operator==(const FrequencyDistribution &histo) const

{
  bool status = true;
  register int i;


  if ((offset != histo.offset) || (nb_value != histo.nb_value) ||
      (nb_element != histo.nb_element)) {
    status = false;
  }

  else {
    for (i = offset;i < nb_value;i++) {
      if (frequency[i] != histo.frequency[i]) {
        status = false;
        break;
      }
    }
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Translation d'une loi empirique.
 *
 *  arguments : reference sur un objet StatError, parametre de translation.
 *
 *--------------------------------------------------------------*/

DiscreteDistributionData* FrequencyDistribution::shift(StatError &error ,
                                                       int shift_param) const

{
  DiscreteDistributionData *histo;


  error.init();

  if (shift_param < -offset) {
    histo = NULL;
    ostringstream correction_message;
    correction_message << STAT_error[STATR_GREATER_THAN] << " " << -offset;
    error.correction_update(STAT_error[STATR_SHIFT_VALUE] , (correction_message.str()).c_str());
  }

  else {
    histo = new DiscreteDistributionData(*this , 's' , shift_param);
  }

  return histo;
}


/*--------------------------------------------------------------*
 *
 *  Regroupement des valeurs d'une loi empirique.
 *
 *  arguments : reference sur un objet StatError, pas pour le regroupement,
 *              mode (FLOOR/ROUND/CEIL).
 *
 *--------------------------------------------------------------*/

DiscreteDistributionData* FrequencyDistribution::cluster(StatError &error ,
                                                         int step , int mode) const

{
  DiscreteDistributionData *histo;


  error.init();

  if (step < 1) {
    histo = NULL;
    error.update(STAT_error[STATR_CLUSTERING_STEP]);
  }

  else {
    histo = new DiscreteDistributionData(*this , 'c' , step , mode);
  }

  return histo;
}


/*--------------------------------------------------------------*
 *
 *  Regroupement des valeurs d'une loi empirique en fonction
 *  de l'augmentation de la quantite d'information.
 *
 *  arguments : reference sur un objet StatError, proportion de la
 *              quantite d'information de la loi empirique initial, stream.
 *
 *--------------------------------------------------------------*/

DiscreteDistributionData* FrequencyDistribution::cluster(StatError &error ,
                                                         double ratio , ostream &os) const

{
  bool status = true , stop = false;
  register int i;
  int step = 1 , *pfrequency , *cfrequency;
  double information , reference_information , previous_information;
  DiscreteDistributionData *histo , *previous_histo;


  previous_histo = NULL;
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
    previous_histo = new DiscreteDistributionData(*this);
    histo = new DiscreteDistributionData((nb_value - 1) / 2 + 1);
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
          previous_histo = new DiscreteDistributionData(*histo);
        }
        else {
          step--;
        }
      }

      else {
        previous_information = information;
        delete previous_histo;
        previous_histo = new DiscreteDistributionData(*histo);
      }
    }
    while (!stop);

    delete histo;

#   ifdef MESSAGE
    os << STAT_label[STATL_INFORMATION_RATIO] << ": "
       << previous_information / reference_information << "   "
       << STAT_label[STATL_CLUSTERING_STEP] << ": " << step << endl;
#   endif

    // calcul des caracteristiques de la loi empirique

    previous_histo->max_computation();
    previous_histo->mean_computation();
    previous_histo->variance_computation();
  }

  return previous_histo;
}


/*--------------------------------------------------------------*
 *
 *  Regroupement des valeurs d'une loi empirique.
 *
 *  arguments : reference sur un objet StatError, nombres de classes,
 *              bornes pour regrouper les valeurs.
 *
 *--------------------------------------------------------------*/

DiscreteDistributionData* FrequencyDistribution::cluster(StatError &error ,
                                                         int nb_class , int *ilimit) const

{
  bool status = true;
  register int i , j;
  int *pfrequency , *cfrequency , *limit;
  DiscreteDistributionData *histo;


  histo = NULL;
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
      histo = new DiscreteDistributionData(nb_class);

      // regroupement des valeurs

      pfrequency = histo->frequency - 1;
      cfrequency = frequency;

      for (i = 0;i < histo->nb_value;i++) {
        *++pfrequency = *cfrequency++;
        for (j = limit[i] + 1;j < limit[i + 1];j++) {
          *pfrequency += *cfrequency++;
        }
      }

      // calcul des caracteristiques de la loi empirique

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
 *  arguments : reference sur un objet StatError,
 *              table de transcodage des symboles.
 *
 *--------------------------------------------------------------*/

DiscreteDistributionData* FrequencyDistribution::transcode(StatError &error ,
                                                           int *symbol) const

{
  bool status = true , *presence;
  register int i;
  int min_symbol , max_symbol , *cfrequency;
  DiscreteDistributionData *histo;


  histo = NULL;
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
    histo = new DiscreteDistributionData(max_symbol + 1);

    // transcodage des symboles

    for (i = 0;i < histo->nb_value;i++) {
      histo->frequency[i] = 0;
    }

    cfrequency = frequency + offset;
    for (i = 0;i < nb_value - offset;i++) {
      histo->frequency[symbol[i]] += *cfrequency++;
    }

    // calcul des caracteristiques de la loi empirique

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
 *  arguments : reference sur un objet StatError, bornes sur les valeurs,
 *              flag pour conserver ou rejeter les individus selectionnes.
 *
 *--------------------------------------------------------------*/

DiscreteDistributionData* FrequencyDistribution::value_select(StatError &error , int min_value ,
                                                              int max_value , bool keep) const

{
  bool status = true;
  register int i;
  DiscreteDistributionData *histo;


  histo = NULL;
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
      histo = new DiscreteDistributionData(nb_value);

      for (i = 0;i < min_value;i++) {
        histo->frequency[i] = frequency[i];
      }
      for (i = min_value;i <= MIN(max_value , nb_value - 1);i++) {
        histo->frequency[i] = 0;
      }
      if (max_value + 1 < nb_value) {
        for (i = max_value + 1;i < nb_value;i++) {
          histo->frequency[i] = frequency[i];
        }
      }

      break;
    }

    case true : {
      histo = new DiscreteDistributionData(MIN(max_value + 1 , nb_value));

      // copie des valeurs

       for (i = 0;i < min_value;i++) {
        histo->frequency[i] = 0;
      }
      for (i = min_value;i < histo->nb_value;i++) {
        histo->frequency[i] = frequency[i];
      }
      break;
    }
    }

    // calcul des caracteristiques de la loi empirique

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
      histo = NULL;
      error.update(STAT_error[STATR_EMPTY_SAMPLE]);
    }
  }

  return histo;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'une loi empirique.
 *
 *  arguments : stream, flag commentaire,
 *              flag sur l'ecriture de la fonction de repartition.
 *
 *--------------------------------------------------------------*/

ostream& FrequencyDistribution::ascii_print(ostream &os , int comment_flag ,
                                            bool cumul_flag) const

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
 *  Ecriture d'un objet FrequencyDistribution.
 *
 *  arguments : stream, flag niveau de detail, flag fichier.
 *
 *--------------------------------------------------------------*/

ostream& FrequencyDistribution::ascii_write(ostream &os , bool exhaustive ,
                                            bool file_flag) const

{
  double information = information_computation();


  if (exhaustive && file_flag) {
    os << "# ";
  }
  os << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";
  ascii_characteristic_print(os , exhaustive , file_flag);

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
    os << "   | " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " | " << STAT_label[STATL_CUMULATIVE]
       << " " << STAT_label[STATL_DISTRIBUTION] << " " << STAT_label[STATL_FUNCTION] << endl;
    ascii_print(os , (file_flag ? 2 : 0) , true);
  }

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet FrequencyDistribution dans un fichier.
 *
 *  arguments : reference sur un objet StatError, path.
 *
 *--------------------------------------------------------------*/

bool FrequencyDistribution::ascii_write(StatError &error , const char *path) const

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
 *  Ecriture des caracteristiques d'une loi empirique au format tableur.
 *
 *  arguments : stream, flag ecriture des parametres de forme.
 *
 *--------------------------------------------------------------*/

ostream& FrequencyDistribution::spreadsheet_characteristic_print(ostream &os , bool shape) const

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
 *  Ecriture des caracteristiques d'une loi empirique
 *  pour une variable circulaire au format tableur.
 *
 *  argument : stream.
 *
 *--------------------------------------------------------------*/

ostream& FrequencyDistribution::spreadsheet_circular_characteristic_print(ostream &os) const

{
  double mean_direction[4];


  os << STAT_label[STATL_SAMPLE_SIZE] << "\t" << nb_element << endl;

  mean_direction_computation(mean_direction);

  os << STAT_label[STATL_MEAN_DIRECTION] << "\t" << mean_direction[3];
  if (mean_direction[2] > 0.) {
    os << "\t" << STAT_label[STATL_MEAN_RESULTANT_LENGTH] << "\t" << mean_direction[2]
       << "\t" << STAT_label[STATL_CIRCULAR_STANDARD_DEVIATION] << "\t"
       << 180 * sqrt(-2 * log(mean_direction[2])) / M_PI;
  }
  os << endl;

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'une loi empirique au format tableur.
 *
 *  arguments : stream, flags sur l'ecriture de la fonction de repartition et
 *              de la fonction de concentration.
 *
 *--------------------------------------------------------------*/

ostream& FrequencyDistribution::spreadsheet_print(ostream &os , bool cumul_flag ,
                                                  bool concentration_flag) const

{
  register int i;
  double *cumul, *concentration;


  if ((!cumul_flag) || (variance == D_DEFAULT) || (variance == 0.)) {
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
 *  Ecriture d'une loi empirique au format Gnuplot.
 *
 *  arguments : path, pointeurs sur les fonctions de repartition et
 *              de concentration, decalage de la loi empirique.
 *
 *--------------------------------------------------------------*/

bool FrequencyDistribution::plot_print(const char *path , double *cumul ,
                                       double *concentration , double shift) const

{
  bool status = false;
  register int i;
  ofstream out_file(path);


  if (out_file) {
    status = true;

    if ((offset == 0) && (((variance == D_DEFAULT) && (nb_value > 1)) ||
         (variance > 0.))) {
      out_file << -1 << " " << -1 + shift << " " << 0 << " "
               << 0 << " " << 0 << " " << 0 << endl;
    }

    for (i = 0;i < nb_value;i++) {
      out_file << i << " " << i + shift << " " << frequency[i] << " "
               << (double)frequency[i] / (double)nb_element;
      if (((variance == D_DEFAULT) && (nb_value > offset + 1)) || (variance > 0.)) {
        out_file << " " << cumul[i];
      }
      if ((variance != D_DEFAULT) && (variance > 0.)) {
        out_file << " " << concentration[i];
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
 *  Ecriture d'une famille de lois empiriques au format Gnuplot.
 *
 *  arguments : path, nombre de lois empiriques,
 *              pointeurs sur les lois empiriques.
 *
 *--------------------------------------------------------------*/

bool FrequencyDistribution::plot_print(const char *path , int nb_histo ,
                                       const FrequencyDistribution **histo) const

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
 *  Ecriture d'une loi empirique.
 *
 *  argument : reference sur un objet SinglePlot.
 *
 *--------------------------------------------------------------*/

void FrequencyDistribution::plotable_frequency_write(SinglePlot &plot) const

{
  register int i;


  for (i = offset;i < nb_value;i++) {
    if (frequency[i] > 0) {
      plot.add_point(i , frequency[i]);
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Ecriture de la loi deduite de la loi empirique.
 *
 *  argument : reference sur un objet SinglePlot.
 *
 *--------------------------------------------------------------*/

void FrequencyDistribution::plotable_mass_write(SinglePlot &plot) const

{
  register int i;


  for (i =  MAX(offset - 1 , 0);i < nb_value;i++) {
    plot.add_point(i , (double)frequency[i] / (double)nb_element);
  }
  if ((double)frequency[nb_value - 1] / (double)nb_element > PLOT_MASS_THRESHOLD) {
    plot.add_point(nb_value , 0.);
  }
}


/*--------------------------------------------------------------*
 *
 *  Ecriture de la fonction de repartition deduite d'une loi empirique.
 *
 *  arguments : reference sur un objet SinglePlot, pointeur sur la fonction de repartition,
 *              facteur d'echelle.
 *
 *--------------------------------------------------------------*/

void FrequencyDistribution::plotable_cumul_write(SinglePlot &plot , double *icumul ,
                                                 double scale) const

{
  register int i;
  double *cumul;


  if (icumul) {
    cumul = icumul;
  }
  else {
    cumul = cumul_computation(scale);
  }

  for (i = MAX(offset - 1 , 0);i < nb_value;i++) {
    plot.add_point(i , cumul[i]);
  }

  if (!icumul) {
    delete [] cumul;
  }
}


/*--------------------------------------------------------------*
 *
 *  Ecriture de la mise en correspondance d'une fonction de repartition avec
 *  la fonction de repartition d'une loi de reference.
 *
 *  arguments : reference sur un objet SinglePlot, bornes et reference
 *              sur la fonction de repartition de reference,
 *              pointeur sur la fonction de repartition.
 *
 *--------------------------------------------------------------*/

void FrequencyDistribution::plotable_cumul_matching_write(SinglePlot &plot , int reference_offset ,
                                                          int  reference_nb_value , double *reference_cumul ,
                                                          double *icumul) const

{
  register int i;
  double *cumul;


  if (icumul) {
    cumul = icumul;
  }
  else {
    cumul = cumul_computation();
  }

  plot.add_point(0. , 0.);
  for (i = MIN(reference_offset , 1);i < offset;i++) {
    plot.add_point(reference_cumul[i] , 0.);
  }
  for (i = offset;i < nb_value;i++) {
    plot.add_point(reference_cumul[i] , cumul[i]);
  }
  for (i = nb_value;i < reference_nb_value;i++) {
    plot.add_point(reference_cumul[i] , cumul[nb_value - 1]);
  }

  if (!icumul) {
    delete [] cumul;
  }
}


/*--------------------------------------------------------------*
 *
 *  Ecriture de la courbe de concentration deduite d'une loi empirique.
 *
 *  argument : reference sur un objet SinglePlot, pointeur sur la fonction de repartition,
 *             facteur d'echelle.
 *
 *--------------------------------------------------------------*/

void FrequencyDistribution::plotable_concentration_write(SinglePlot &plot , double *icumul ,
                                                         double scale) const

{
  register int i;
  double *cumul , *concentration;


  if (icumul) {
    cumul = icumul;
  }
  else {
    cumul = cumul_computation(scale);
  }

  concentration = concentration_function_computation(scale);

  plot.add_point(0. , 0.);
  for (i = offset;i < nb_value;i++) {
    plot.add_point(cumul[i] , concentration[i]);
  }

  if (!icumul) {
    delete [] cumul;
  }
  delete [] concentration;
}


/*--------------------------------------------------------------*
 *
 *  Calcul des taux de survie a partir d'une loi empirique et
 *  ecriture du resultat.
 *
 *  argument : stream.
 *
 *--------------------------------------------------------------*/

ostream& FrequencyDistribution::survival_ascii_write(ostream &os) const

{
  Curves *survival_rate;


  os << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";
  ascii_characteristic_print(os);

  os << "\n   | " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " | " << STAT_label[STATL_CUMULATIVE]
     << " " << STAT_label[STATL_DISTRIBUTION] << " " << STAT_label[STATL_FUNCTION] << endl;
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
 *  Calcul des taux de survie a partir d'une loi empirique et
 *  ecriture du resultat dans un fichier.
 *
 *  arguments : reference sur un objet StatError, path.
 *
 *--------------------------------------------------------------*/

bool FrequencyDistribution::survival_ascii_write(StatError &error , const char *path) const

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
 *  Calcul des taux de survie a partir d'une loi empirique et ecriture
 *  du resultat dans un fichier au format tableur.
 *
 *  arguments : reference sur un objet StatError, path.
 *
 *--------------------------------------------------------------*/

bool FrequencyDistribution::survival_spreadsheet_write(StatError &error , const char *path) const

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

    out_file << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\t";
    spreadsheet_characteristic_print(out_file);

    out_file << "\n\t" << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\t" << STAT_label[STATL_CUMULATIVE]
             << " " << STAT_label[STATL_DISTRIBUTION] << " " << STAT_label[STATL_FUNCTION] << endl;
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
 *  Ecriture d'une loi empirique, de la loi et de la fonction de survie deduites
 *  au format Gnuplot.
 *
 *  arguments : path, pointeur sur la fonction de survie.
 *
 *--------------------------------------------------------------*/

bool FrequencyDistribution::survival_plot_print(const char *path , double *survivor) const

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
 *  Calcul des taux de survie a partir d'une loi empirique et
 *  sortie Gnuplot du resultat.
 *
 *  arguments : reference sur un objet StatError, prefixe des fichiers,
 *              titre des figures.
 *
 *--------------------------------------------------------------*/

bool FrequencyDistribution::survival_plot_write(StatError &error , const char *prefix ,
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

        // loi empirique

        if (nb_value - 1 < TIC_THRESHOLD) {
          out_file << "set xtics 0,1" << endl;
        }
        if ((int)(max * YSCALE) + 1 < TIC_THRESHOLD) {
          out_file << "set ytics 0,1" << endl;
        }

        out_file << "plot [0:" << nb_value - 1 << "] [0:"
                 << (int)(max * YSCALE) + 1 << "] \""
                 << label((data_file_name[0].str()).c_str()) << "\" using 1 title \""
                 << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\" with impulses" << endl;

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
 *  Calcul et ecriture de la fonction de survie d'une loi empirique.
 *
 *  argument : reference sur un objet SinglePlot.
 *
 *--------------------------------------------------------------*/

void FrequencyDistribution::plotable_survivor_write(SinglePlot &plot) const

{
  register int i;
  double *survivor;


  survivor = survivor_function_computation();

  for (i = MAX(offset - 1 , 0);i < nb_value;i++) {
    plot.add_point(i , survivor[i]);
  }

  delete [] survivor;
}


/*--------------------------------------------------------------*
 *
 *  Calcul des taux de survie a partir d'une loi empirique et sortie graphique du resultat.
 *
 *  argument : reference sur un objet StatError.
 *
 *--------------------------------------------------------------*/

MultiPlotSet* FrequencyDistribution::survival_get_plotable(StatError &error) const

{
  MultiPlotSet *plot_set;


  error.init();

  if (variance == 0.) {
    plot_set = NULL;
    error.update(STAT_error[STATR_PLOT_NULL_VARIANCE]);
  }

  else {
    register int i , j;
    int xmax;
    Curves *survival_rate;
    ostringstream legend;


    plot_set = new MultiPlotSet(3);
    MultiPlotSet &plot = *plot_set;

    plot.title = "Survival analysis";
    plot.border = "15 lw 0";

    // 1ere vue : loi empirique

    plot[0].xrange = Range(0 , nb_value - 1);
    plot[0].yrange = Range(0 , ceil(max * YSCALE));

    if (nb_value - 1 < TIC_THRESHOLD) {
      plot[0].xtics = 1;
    }

    plot[0].resize(1);

    plot[0][0].legend = STAT_label[STATL_FREQUENCY_DISTRIBUTION];

    plot[0][0].style = "impulses";

    plotable_frequency_write(plot[0][0]);

    // 2eme vue : loi et fonction de survie

    xmax = nb_value - 1;
    if ((double)frequency[xmax] / (double)nb_element > PLOT_MASS_THRESHOLD) {
      xmax++;
    }
    plot[1].xrange = Range(0 , xmax);

    plot[1].yrange = Range(0. , 1.);

    if (nb_value - 1 < TIC_THRESHOLD) {
      plot[1].xtics = 1;
    }

    plot[1].resize(2);

    plot[1][0].legend = STAT_label[STATL_DISTRIBUTION];

    plot[1][0].style = "linespoints";

    plotable_mass_write(plot[1][0]);

    legend.str("");
    legend << STAT_label[STATL_SURVIVOR] << " " << STAT_label[STATL_FUNCTION];
    plot[1][1].legend = legend.str();

    plot[1][1].style = "linespoints";

    plotable_survivor_write(plot[1][1]);

    // 3eme vue : taux de survie

    survival_rate = new Curves(*this);

    plot[2].xrange = Range(survival_rate->offset , survival_rate->length - 1);
    plot[2].yrange = Range(0. , 1.);

    if (survival_rate->length - 1 < TIC_THRESHOLD) {
      plot[2].xtics = 1;
    }

    plot[2].resize(2);

    plot[2][0].legend = STAT_label[STATL_DEATH_PROBABILITY];
    plot[2][0].style = "linespoints";

    plot[2][1].legend = STAT_label[STATL_SURVIVAL_PROBABILITY];
    plot[2][1].style = "linespoints";

    survival_rate->plotable_write(plot[2]);

    delete survival_rate;
  }

  return plot_set;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la fonction de repartition deduite d'une loi empirique.
 *
 *  argument : facteur d'echelle.
 *
 *--------------------------------------------------------------*/

double* FrequencyDistribution::cumul_computation(double scale) const

{
  register int i;
  double *cumul;


  if (scale == D_DEFAULT) {
    scale = nb_element;
  }

  cumul = new double[nb_value];

  for (i = 0;i < offset;i++) {
    cumul[i] = 0.;
  }
  cumul[offset] = frequency[offset] / scale;
  for (i = offset + 1;i < nb_value;i++) {
    cumul[i] = cumul[i - 1] + frequency[i] / scale;
  }

  return cumul;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la fonction de survie deduite d'une loi empirique.
 *
 *  argument : facteur d'echelle.
 *
 *--------------------------------------------------------------*/

double* FrequencyDistribution::survivor_function_computation(double scale) const

{
  register int i;
  double *survivor_function;


  if (scale == D_DEFAULT) {
    scale = nb_element;
  }

  survivor_function = new double[nb_value];

  survivor_function[nb_value - 1] = 0.;
  for (i = nb_value - 2;i >= MAX(offset - 1 , 0);i--) {
    survivor_function[i] = survivor_function[i + 1] + frequency[i + 1] / scale;
  }
  for (i = offset - 2;i >= 0;i--) {
    survivor_function[i] = survivor_function[i + 1];
  }

  return survivor_function;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la fonction de concentration deduite d'une loi empirique.
 *
 *  argument : facteur d'echelle.
 *
 *--------------------------------------------------------------*/

double* FrequencyDistribution::concentration_function_computation(double scale) const

{
  register int i;
  double norm , *concentration_function;


  if ((variance != D_DEFAULT) && (variance > 0.)) {
    if (scale == D_DEFAULT) {
      scale = nb_element;
    }

    concentration_function = new double[nb_value];

    for (i = 0;i < offset;i++) {
      concentration_function[i] = 0.;
    }
    norm = mean * scale;

    concentration_function[offset] = frequency[offset] * offset / norm;
    for (i = offset + 1;i < nb_value;i++) {
      concentration_function[i] = concentration_function[i - 1] + frequency[i] * i / norm;
    }
  }

  else {
    concentration_function = NULL;
  }

  return concentration_function;
}


/*--------------------------------------------------------------*
 *
 *  Calcul du coefficient de concentration d'une loi empirique.
 *
 *--------------------------------------------------------------*/

double FrequencyDistribution::concentration_computation() const

{
  register int i;
  double concentration = D_DEFAULT , *concentration_function;


  if ((mean > 0.) && (variance > 0.)) {
    concentration_function = concentration_function_computation();

    concentration = frequency[offset] * concentration_function[offset];
    for (i = offset + 1;i < nb_value;i++) {
      concentration += frequency[i] * (concentration_function[i - 1] + concentration_function[i]);
    }
    concentration = 1. - concentration / nb_element;

    delete [] concentration_function;

#   ifdef DEBUG
    int previous_value , cumul;
    double concentration2 = 0.;

    cumul = frequency[offset];
    previous_value = offset;

    for (i = offset + 1;i < nb_value;i++) {
      if (frequency[i] > 0) {
        concentration2 += cumul * (double)(nb_element - cumul) * (i - previous_value);
        cumul += frequency[i];
        previous_value = i;
      }
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
 *  Mise a jour d'une loi empirique a frequences entieres a partir
 *  d'une loi empirique a frequences reelles par arrondi.
 *
 *  arguments : pointeur sur une loi empirique a frequences reelles,
 *              effectif de la loi empirique.
 *
 *--------------------------------------------------------------*/

void FrequencyDistribution::update(const Reestimation<double> *reestim , int inb_element)

{
  register int i , j;
  int index;
  double scale , sum , max_frequency , *real_frequency;


  // copie de la loi empirique reel et mise a l'echelle

  real_frequency = new double[reestim->nb_value];

  scale = inb_element / reestim->nb_element;
  for (i = reestim->offset;i < reestim->nb_value;i++) {
    real_frequency[i] = reestim->frequency[i] * scale;
  }

  // calcul des frequences

  for (i = 0;i < reestim->offset;i++) {
    frequency[i] = 0;
  }

  sum = 0.;
  for (i = reestim->offset;i < reestim->nb_value;i++) {
    frequency[i] = (int)real_frequency[i];
    real_frequency[i] -= frequency[i];
    if (real_frequency[i] > 0.) {
      sum += real_frequency[i];
    }
  }

  for (i = reestim->nb_value;i < alloc_nb_value;i++) {
    frequency[i] = 0;
  }

  // prise en compte des arrondis

  for (i = 0;i < (int)round(sum);i++) {
    max_frequency = 0.;
    for (j = reestim->offset;j < reestim->nb_value;j++) {
      if (real_frequency[j] > max_frequency) {
        max_frequency = real_frequency[j];
        index = j;
      }
    }

    real_frequency[index] = 0.;
    frequency[index]++;
  }

  // calcul des caracteristiques de la loi empirique

  nb_value_computation();
  offset_computation();
  nb_element_computation();
  max_computation();
  mean_computation();
  variance_computation();

  delete [] real_frequency;
}


/*--------------------------------------------------------------*
 *
 *  Calcul d'une loi empirique a partir d'une loi empirique initial
 *  en changeant l'effectif.
 *
 *  argument : effectif total.
 *
 *--------------------------------------------------------------*/

FrequencyDistribution* FrequencyDistribution::frequency_scale(int inb_element) const

{
  register int i , j;
  int index;
  double sum , real_max , *real_frequency;
  FrequencyDistribution *histo;


  real_frequency = new double[nb_value];
  histo = new FrequencyDistribution(nb_value);

  sum = 0.;
  for (i = offset;i < nb_value;i++) {
    if (frequency[i] > 0) {
      real_frequency[i] = frequency[i] * (double)inb_element / (double)nb_element;
      histo->frequency[i] = (int)real_frequency[i];
      real_frequency[i] -= histo->frequency[i];
      if (real_frequency[i] > 0.) {
        sum += real_frequency[i];
      }
    }

    else {
      real_frequency[i] = 0.;
    }
  }

  // prise en compte des arrondis

  for (i = 0;i < (int)round(sum);i++) {
    real_max = 0.;
    for (j = offset;j < nb_value;j++) {
      if (real_frequency[j] > real_max) {
        real_max = real_frequency[j];
        index = j;
      }
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

  // calcul des caracteristiques de la loi empirique

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
 *  Calcul des rangs a partir d'une loi empirique.
 *
 *--------------------------------------------------------------*/

double* FrequencyDistribution::rank_computation() const

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


/*--------------------------------------------------------------*
 *
 *  Calcul de la fonction de repartition empirique.
 *
 *  argument : (valeurs, fonction de repartition).
 *
 *--------------------------------------------------------------*/

int FrequencyDistribution::cumulative_distribution_function_computation(double **cdf) const

{
  register int i , j;
  int buff , cumul;


  buff = MIN(nb_value - offset , nb_element);
  cdf[0] = new double[buff];
  cdf[1] = new double[buff];

  cumul = 0;
  i = 0;
  for (j = offset;j < nb_value;j++) {
    if (frequency[j] > 0) {
      cdf[0][i] = j;
      cdf[1][i] = (cumul + (double)(frequency[j] + 1) / 2.) /  (double)nb_element;
      cumul += frequency[j];
      i++;
    }
  }

  return i;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de l'intervalle minimum entre 2 valeurs.
 *
 *--------------------------------------------------------------*/

int FrequencyDistribution::min_interval_computation() const

{
  register int i;
  int min_interval , previous_value;


  min_interval = nb_value;
  previous_value = offset;

  for (i = offset + 1;i < nb_value;i++) {
    if (frequency[i] > 0) {
      if (i - previous_value < min_interval) {
        min_interval = i - previous_value;
      }
      previous_value = i;
    }
  }

  return min_interval;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la vraisemblance d'une loi continue pour un echantillon.
 *
 *  arguments : reference sur un objet ContinuousParametric,
 *              intervalle minimum entre 2 valeurs.
 *
 *--------------------------------------------------------------*/

double FrequencyDistribution::likelihood_computation(const ContinuousParametric &dist ,
                                                     int min_interval) const

{
  register int i;
  double mass , likelihood = 0.;


  if (nb_element > 0) {
    if (min_interval == I_DEFAULT) {
      min_interval = min_interval_computation();
    }

    for (i = offset;i < nb_value;i++) {
      if (frequency[i] > 0) {
        if ((dist.ident == GAMMA) && (offset < (double)min_interval / 2.)) {
          mass = dist.mass_computation(i , i + (double)min_interval);
        }
        else {
          mass = dist.mass_computation(i - (double)min_interval / 2. , i + (double)min_interval / 2.);
        }

        if (mass > 0.) {
          likelihood += frequency[i] * log(mass);
        }
        else {
          likelihood = D_INF;
          break;
        }
      }
    }
  }

  return likelihood;
}


};  // namespace stat_tool
