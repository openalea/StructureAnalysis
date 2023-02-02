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



#include <limits.h>
#include <math.h>

#include <string>
#include <vector>
#include <sstream>
#include <iomanip>

#include "distribution.h"
#include "curves.h"
#include "stat_label.h"

using namespace std;


namespace stat_tool {



/*--------------------------------------------------------------*/
/**
 *  \brief Constructor of the FrequencyDistribution class.
 *
 *  \param[in] inb_element number of individuals,
 *  \param[in] pelement    individuals.
 */
/*--------------------------------------------------------------*/

FrequencyDistribution::FrequencyDistribution(int inb_element , int *pelement)

{
  int i;


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

  // computation of the frequency distribution characteristics

  offset_computation();
  max_computation();
  mean_computation();
  variance_computation();
}


/*--------------------------------------------------------------*/
/**
 *  \brief Shifting of a frequency distribution.
 *
 *  \param[in] histo       reference on a FrequencyDistribution object,
 *  \param[in] shift_param shifting parameter.
 */
/*--------------------------------------------------------------*/

void FrequencyDistribution::shift(const FrequencyDistribution &histo , int shift_param)

{
  int i;
  int *cfrequency;


  // computation of the frequency distribution characteristics

  nb_element = histo.nb_element;
  nb_value = histo.nb_value + shift_param;
  alloc_nb_value = nb_value;
  offset = histo.offset + shift_param;
  max = histo.max;
  mean = histo.mean + shift_param;
  variance = histo.variance;

  // copy of frequencies

  frequency = new int[nb_value];

  for (i = 0;i < offset;i++) {
    frequency[i] = 0;
  }
  cfrequency = histo.frequency + histo.offset;
  for (i = offset;i < nb_value;i++) {
    frequency[i] = *cfrequency++;
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Clustering of values of a frequency distribution.
 *
 *  \param[in] histo reference on a FrequencyDistribution object,
 *  \param[in] step  clustering step,
 *  \param[in] mode  mode (FLOOR/ROUND/CEIL).
 */
/*--------------------------------------------------------------*/

void FrequencyDistribution::cluster(const FrequencyDistribution &histo ,
                                    int step , rounding mode)

{
  int i;


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

  // cumul of frequencies

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

  // computation of the frequency distribution characteristics

  max_computation();
  mean_computation();
  variance_computation();
}


/*--------------------------------------------------------------*/
/**
 *  \brief Constructor of the FrequencyDistribution class.
 *
 *  \param[in] histo     reference on a FrequencyDistribution object,
 *  \param[in] transform type of transform (SHIFT/CLUSTER),
 *  \param[in] param     shifting parameter (SHIFT) / clustering step (CLUSTER),
 *  \param[in] mode      clustering mode (FLOOR/ROUND/CEIL).
 */
/*--------------------------------------------------------------*/

FrequencyDistribution::FrequencyDistribution(const FrequencyDistribution &histo ,
                                             frequency_distribution_transformation transform ,
                                             int param , rounding mode)

{
  switch (transform) {
  case SHIFT :
    shift(histo , param);
    break;
  case CLUSTER :
    cluster(histo , param , mode);
    break;
  default :
    copy(histo);
    break;
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Equality operator of the FrequencyDistribution class.
 *
 *  \param[in] histo reference on a FrequencyDistribution object.
 *
 *  \return          equality or not of the frequency distributions.
 */
/*--------------------------------------------------------------*/

bool FrequencyDistribution::operator==(const FrequencyDistribution &histo) const

{
  bool status = true;
  int i;


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


/*--------------------------------------------------------------*/
/**
 *  \brief Merging of FrequencyDistribution objects.
 *
 *  \param[in] nb_sample number of FrequencyDistribution objects,
 *  \param[in] ihisto    pointer on the FrequencyDistribution objects.
 *
 *  \return              DiscreteDistributionData object.
 */
/*--------------------------------------------------------------*/

DiscreteDistributionData* FrequencyDistribution::merge(int nb_sample ,
                                                       const vector<FrequencyDistribution> ihisto) const

{
  int i;
  DiscreteDistributionData *histo;
  const FrequencyDistribution **phisto;


  nb_sample++;
  phisto = new const FrequencyDistribution*[nb_sample];

  phisto[0] = this;
  for (i = 1;i < nb_sample;i++) {
    phisto[i] = new FrequencyDistribution(ihisto[i - 1]);
  }

  histo = new DiscreteDistributionData(nb_sample , phisto);

  for (i = 1;i < nb_sample;i++) {
    delete phisto[i];
  }
  delete [] phisto;

  return histo;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Shifting of a frequency distribution.
 *
 *  \param[in] error       reference on a StatError object,
 *  \param[in] shift_param shifting parameter.
 *
 *  \return                DiscreteDistributionData object.
 */
/*--------------------------------------------------------------*/

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
    histo = new DiscreteDistributionData(*this , SHIFT , shift_param);
  }

  return histo;
}


/*--------------------------------------------------------------*/
/**
 *  \brief clustering of values of a frequency distribution.
 *
 *  \param[in] error reference on a StatError object,
 *  \param[in] step  clustering step,
 *  \param[in] mode  mode (FLOOR/ROUND/CEIL).
 *
 *  \return          DiscreteDistributionData object.
 */
/*--------------------------------------------------------------*/

DiscreteDistributionData* FrequencyDistribution::cluster(StatError &error , int step ,
                                                         rounding mode) const

{
  DiscreteDistributionData *histo;


  error.init();

  if (step < 1) {
    histo = NULL;
    error.update(STAT_error[STATR_CLUSTERING_STEP]);
  }

  else {
    histo = new DiscreteDistributionData(*this , CLUSTER , step , mode);
  }

  return histo;
}


/*--------------------------------------------------------------*/
/**
 *  \brief clustering of values of a frequency distribution using
 *         an information quantity criterion.
 *
 *  \param[in] error   reference on a StatError object,
 *  \param[in] ratio   proportion of the information quantity of
 *                     the initial frequency distribution,
 *  \param[in] display flag for displaying the clustering step.
 *
 *  \return            DiscreteDistributionData object.
 */
/*--------------------------------------------------------------*/

DiscreteDistributionData* FrequencyDistribution::cluster(StatError &error , double ratio ,
                                                         bool display) const

{
  bool status = true , stop = false;
  int i;
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

      // clustering of values

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

      // computation of the information quantity

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

    if (display) {
      cout << STAT_label[STATL_INFORMATION_RATIO] << ": "
           << previous_information / reference_information << "   "
           << STAT_label[STATL_CLUSTERING_STEP] << ": " << step << endl;
    }

    // computation of the frequency distribution characteristics

    previous_histo->max_computation();
    previous_histo->mean_computation();
    previous_histo->variance_computation();
  }

  return previous_histo;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Partitioning of values of a frequency distribution.
 *
 *  \param[in] error    reference on a StatError object,
 *  \param[in] nb_class number of classes,
 *  \param[in] ilimit   limits between classes (beginning of classes).
 *
 *  \return             DiscreteDistributionData object.
 */
/*--------------------------------------------------------------*/

DiscreteDistributionData* FrequencyDistribution::cluster(StatError &error , int nb_class ,
                                                         int *ilimit) const

{
  bool status = true;
  int i , j;
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

      // partitioning of values

      pfrequency = histo->frequency - 1;
      cfrequency = frequency;

      for (i = 0;i < histo->nb_value;i++) {
        *++pfrequency = *cfrequency++;
        for (j = limit[i] + 1;j < limit[i + 1];j++) {
          *pfrequency += *cfrequency++;
        }
      }

      // computation of the frequency distribution characteristics

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


/*--------------------------------------------------------------*/
/**
 *  \brief Partitioning of values of a frequency distribution.
 *
 *  \param[in] error    reference on a StatError object,
 *  \param[in] nb_class number of classes,
 *  \param[in] ilimit   limits between classes (beginning of classes).
 *
 *  \return             DiscreteDistributionData object.
 */
/*--------------------------------------------------------------*/

DiscreteDistributionData* FrequencyDistribution::cluster(StatError &error , int nb_class ,
                                                         vector<int> ilimit) const

{
  return cluster(error , nb_class , ilimit.data());
}


/*--------------------------------------------------------------*/
/**
 *  \brief Transcoding of categories.
 *
 *  \param[in] error    reference on a StatError object,
 *  \param[in] category transcoding table.
 *
 *  \return             DiscreteDistributionData object.
 */
/*--------------------------------------------------------------*/

DiscreteDistributionData* FrequencyDistribution::transcode(StatError &error ,
                                                           int *category) const

{
  bool status = true , *presence;
  int i;
  int min_category , max_category , *cfrequency;
  DiscreteDistributionData *histo;


  histo = NULL;
  error.init();

  min_category = INT_MAX;
  max_category = 0;

  for (i = 0;i < nb_value - offset;i++) {
    if (category[i] < 0) {
      status = false;
      ostringstream error_message;
      error_message << STAT_label[STATL_CATEGORY] << " " << category[i] << " "
                    << STAT_error[STATR_NOT_ALLOWED];
      error.update((error_message.str()).c_str());
    }
    else {
      if (category[i] < min_category) {
        min_category = category[i];
      }
      if (category[i] > max_category) {
        max_category = category[i];
      }
    }
  }

  if (max_category - min_category == 0) {
    status = false;
    error.update(STAT_error[STATR_NB_CATEGORY]);
  }

  if (max_category - min_category > nb_value - 1 - offset) {
    status = false;
    error.update(STAT_error[STATR_NON_CONSECUTIVE_CATEGORIES]);
  }

  if (status) {
    presence = new bool[max_category + 1];
    for (i = min_category;i <= max_category;i++) {
      presence[i] = false;
    }

    for (i = 0;i < nb_value - offset;i++) {
      presence[category[i]] = true;
    }

    for (i = min_category;i <= max_category;i++) {
      if (!presence[i]) {
        status = false;
        ostringstream error_message;
        error_message << STAT_error[STATR_MISSING_CATEGORY] << " " << i;
        error.update((error_message.str()).c_str());
      }
    }

    delete [] presence;
  }

  if (status) {
    histo = new DiscreteDistributionData(max_category + 1);

    // transcoding of categories

    for (i = 0;i < histo->nb_value;i++) {
      histo->frequency[i] = 0;
    }

    cfrequency = frequency + offset;
    for (i = 0;i < nb_value - offset;i++) {
      histo->frequency[category[i]] += *cfrequency++;
    }

    // computation of the frequency distribution characteristics

    histo->offset = min_category;
    histo->nb_element = nb_element;
    histo->max_computation();
    histo->mean_computation();
    histo->variance_computation();
  }

  return histo;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Transcoding of categories.
 *
 *  \param[in] error    reference on a StatError object,
 *  \param[in] category transcoding table.
 *
 *  \return             DiscreteDistributionData object.
 */
/*--------------------------------------------------------------*/

DiscreteDistributionData* FrequencyDistribution::transcode(StatError &error ,
                                                           vector<int> category) const

{
  return transcode(error , category.data());
}


/*--------------------------------------------------------------*/
/**
 *  \brief Selection of a range of values.
 *
 *  \param[in] error     reference on a StatError object,
 *  \param[in] min_value lowest value,
 *  \param[in] max_value highest value,
 *  \param[in] keep      flag for keeping or rejecting the selected values.
 *
 *  \return              DiscreteDistributionData object.
 */
/*--------------------------------------------------------------*/

DiscreteDistributionData* FrequencyDistribution::value_select(StatError &error , int min_value ,
                                                              int max_value , bool keep) const

{
  bool status = true;
  int i;
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
    if (keep) {
      histo = new DiscreteDistributionData(MIN(max_value + 1 , nb_value));

      // copy of frequencies

       for (i = 0;i < min_value;i++) {
        histo->frequency[i] = 0;
      }
      for (i = min_value;i < histo->nb_value;i++) {
        histo->frequency[i] = frequency[i];
      }
    }

    else {
      histo = new DiscreteDistributionData(nb_value);

      // copy of frequencies

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
    }

    // computation of the frequency distribution characteristics

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


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a frequency distribution.
 *
 *  \param[in,out] os           stream,
 *  \param[in]     comment_flag flag comment,
 *  \param[in]     cumul_flag   flag on the writing of the cumulative distribution function.
 */
/*--------------------------------------------------------------*/

ostream& FrequencyDistribution::ascii_print(ostream &os , int comment_flag ,
                                            bool cumul_flag) const

{
  int i;
  int width[3];
  double *cumul;
  ios_base::fmtflags format_flags;


  format_flags = os.setf(ios::right , ios::adjustfield);

  // computation of the column widths

  width[0] = column_width(nb_value - 1);
  width[1] = column_width(max) + ASCII_SPACE;

  if (cumul_flag) {
    cumul = cumul_computation();
    width[2] = column_width(nb_value , cumul) + ASCII_SPACE;
  }

  // writing of the frequencies and the cumulative distribution function

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

  os.setf(format_flags , ios::adjustfield);

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a FrequencyDistribution object.
 *
 *  \param[in,out] os         stream,
 *  \param[in]     exhaustive flag detail level,
 *  \param[in]     file_flag  flag file.
 */
/*--------------------------------------------------------------*/

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
  os << STAT_label[STATL_MEAN_ABSOLUTE_DEVIATION] << ": " << mean_absolute_deviation_computation(mean);
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


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a FrequencyDistribution object in a file.
 *
 *  \param[in] error reference on a StatError object,
 *  \param[in] path  file path.
 *
 *  \return          error status.
 */
/*--------------------------------------------------------------*/

bool FrequencyDistribution::ascii_write(StatError &error , const string path) const

{
  bool status;
  ofstream out_file(path.c_str());


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


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of the characteristics of a frequency distribution at the spreadsheet format.
 *
 *  \param[in,out] os    stream,
 *  \param[in]     shape flag on the writing of the shape characteristics.
 */
/*--------------------------------------------------------------*/

ostream& FrequencyDistribution::spreadsheet_characteristic_print(ostream &os , bool shape) const

{
  os << STAT_label[STATL_SAMPLE_SIZE] << "\t" << nb_element << endl;

  if ((mean != D_DEFAULT) && (variance != D_DEFAULT)) {
    os << STAT_label[STATL_MEAN] << "\t" << mean << "\t\t"
       << STAT_label[STATL_MEDIAN] << "\t" << quantile_computation() << "\t\t"
       << STAT_label[STATL_MODE] << "\t" << mode_computation() << endl;

    os << STAT_label[STATL_VARIANCE] << "\t" << variance << "\t\t"
       << STAT_label[STATL_STANDARD_DEVIATION] << "\t" << sqrt(variance);
    if (variance > 0.) {
      os << "\t\t" << STAT_label[STATL_LOWER_QUARTILE] << "\t" << quantile_computation(0.25)
         << "\t\t" << STAT_label[STATL_UPPER_QUARTILE] << "\t" << quantile_computation(0.75);
    }
    os << endl;

    if ((shape) && (variance > 0.)) {
      os << STAT_label[STATL_SKEWNESS_COEFF] << "\t" << skewness_computation() << "\t\t"
         << STAT_label[STATL_KURTOSIS_COEFF] << "\t" << kurtosis_computation() << endl;
    }
  }

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of the characteristics of a frequency distribution of
 *         a circular variable at the spreadsheet format.
 *
 *  \param[in,out] os stream.
 */
/*--------------------------------------------------------------*/

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


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a frequency distribution at the spreadsheet format.
 *
 *  \param[in,out] os                 stream,
 *  \param[in]     cumul_flag         flag on the writing of the cumulative distribution function,
 *  \param[in]     concentration_flag flag on the writing of the concentration function.
 */
/*--------------------------------------------------------------*/

ostream& FrequencyDistribution::spreadsheet_print(ostream &os , bool cumul_flag ,
                                                  bool concentration_flag) const

{
  int i;
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


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a frequency distribution at the Gnuplot format.
 *
 *  \param[in] path          file path,
 *  \param[in] cumul         pointer on the cumulative distribution function,
 *  \param[in] concentration pointer on the concentration function,
 *  \param[in] shift         shift of the frequency distribution.
 *
 *  \return                  error status.
 */
/*--------------------------------------------------------------*/

bool FrequencyDistribution::plot_print(const char *path , double *cumul ,
                                       double *concentration , double shift) const

{
  bool status = false;
  int i;
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


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a family of frequency distributions at the Gnuplot format.
 *
 *  \param[in] path     file path,
 *  \param[in] nb_histo number of frequency distributions,
 *  \param[in] histo    pointer on the frequency distributions.
 *
 *  \return             error status.
 */
/*--------------------------------------------------------------*/

bool FrequencyDistribution::plot_print(const char *path , int nb_histo ,
                                       const FrequencyDistribution **histo) const

{
  bool status = false;
  int i , j;
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


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a frequency distribution.
 *
 *  \param[in] plot reference on a SinglePlot object.
 */
/*--------------------------------------------------------------*/

void FrequencyDistribution::plotable_frequency_write(SinglePlot &plot) const

{
  int i;


  for (i = offset;i < nb_value;i++) {
    if (frequency[i] > 0) {
      plot.add_point(i , frequency[i]);
    }
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of the probability mass function deduced from a frequency distribution.
 *
 *  \param[in] plot reference on a SinglePlot object.
 */
/*--------------------------------------------------------------*/

void FrequencyDistribution::plotable_mass_write(SinglePlot &plot) const

{
  int i;


  for (i = MAX(offset - 1 , 0);i < nb_value;i++) {
    plot.add_point(i , (double)frequency[i] / (double)nb_element);
  }
  if ((double)frequency[nb_value - 1] / (double)nb_element > PLOT_MASS_THRESHOLD) {
    plot.add_point(nb_value , 0.);
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of the cumulative distribution function deduced from a frequency distribution.
 *
 *  \param[in] plot   reference on a SinglePlot object,
 *  \param[in] icumul pointer on the cumulative distribution function,
 *  \param[in] scale  scaling factor.
 */
/*--------------------------------------------------------------*/

void FrequencyDistribution::plotable_cumul_write(SinglePlot &plot , double *icumul ,
                                                 double scale) const

{
  int i;
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


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of the matching of a cumulative distribution function with
 *         a reference cumulative distribution function.
 *
 *  \param[in] plot               reference on a SinglePlot object,
 *  \param[in] reference_offset   reference minimum value,
 *  \param[in] reference_nb_value reference number of values,
 *  \param[in] reference_cumul    pointer on the reference cumulative distribution function,
 *  \param[in] icumul             pointer on the cumulative distribution function.
 */
/*--------------------------------------------------------------*/

void FrequencyDistribution::plotable_cumul_matching_write(SinglePlot &plot , int reference_offset ,
                                                          int  reference_nb_value , double *reference_cumul ,
                                                          double *icumul) const

{
  int i;
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


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of the concentration curve deduced from a frequency distribution.
 *
 *  \param[in] plot   reference on a SinglePlot object,
 *  \param[in] icumul pointer on the cumulative distribution function,
 *  \param[in] scale  scaling factor.
 */
/*--------------------------------------------------------------*/

void FrequencyDistribution::plotable_concentration_write(SinglePlot &plot , double *icumul ,
                                                         double scale) const

{
  int i;
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


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the survival rates from a frequency distribution and
 *         writing of the result.
 *
 *  \param[in,out] os stream.
 */
/*--------------------------------------------------------------*/

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


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the survival rates from a frequency distribution and
 *         writing of the result in a file.
 *
 *  \param[in] error reference on a StatError object,
 *  \param[in] path  file path.
 *
 *  \return          error status.
 */
/*--------------------------------------------------------------*/

bool FrequencyDistribution::survival_ascii_write(StatError &error , const string path) const

{
  bool status;
  ofstream out_file(path.c_str());


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


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the survival rates from a frequency distribution and
 *         writing of the result in a file at the spreadsheet format.
 *
 *  \param[in] error reference on a StatError object,
 *  \param[in] path  file path.
 *
 *  \return          error status.
 */
/*--------------------------------------------------------------*/

bool FrequencyDistribution::survival_spreadsheet_write(StatError &error , const string path) const

{
  bool status;
  Curves *survival_rate;
  ofstream out_file(path.c_str());


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


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a frequency distribution, the deduced probability mass and
 *         survivor functions at the Gnuplot format.
 *
 *  \param[in] path     file path,
 *  \param[in] survivor pointer on the survivor function.
 *
 *  \return             error status.
 */
/*--------------------------------------------------------------*/

bool FrequencyDistribution::survival_plot_print(const char *path , double *survivor) const

{
  bool status = false;
  int i;
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


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the survival rates from a frequency distribution and
 *         plot of the result using Gnuplot.
 *
 *  \param[in] error  reference on a StatError object,
 *  \param[in] prefix file prefix,
 *  \param[in] title  figure title.
 *
 *  \return           error status.
 */
/*--------------------------------------------------------------*/

bool FrequencyDistribution::survival_plot_write(StatError &error , const char *prefix ,
                                                const char *title) const

{
  bool status;
  int i;
  double *survivor;
  Curves *survival_rate;
  ostringstream data_file_name[2];


  error.init();

  if (variance == 0.) {
    status = false;
    error.update(STAT_error[STATR_PLOT_NULL_VARIANCE]);
  }

  else {

    // writing of the data files

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

      // writing of the script files

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

        // frequency distribution

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

        // probability mass and survivor functions

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

        // survival rates

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


/*--------------------------------------------------------------*/
/**
 *  \brief Computation and writing of the survivor function computed from
 *         a frequency distribution.
 *
 *  \param[in] plot reference on a SinglePlot object.
 */
/*--------------------------------------------------------------*/

void FrequencyDistribution::plotable_survivor_write(SinglePlot &plot) const

{
  int i;
  double *survivor;


  survivor = survivor_function_computation();

  for (i = MAX(offset - 1 , 0);i < nb_value;i++) {
    plot.add_point(i , survivor[i]);
  }

  delete [] survivor;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the survival rates from a frequency distribution and
 *         plot of the result.
 *
 *  \param[in] error reference on a StatError object.
 *
 *  \return          plots.
 */
/*--------------------------------------------------------------*/

MultiPlotSet* FrequencyDistribution::survival_get_plotable(StatError &error) const

{
  MultiPlotSet *plot_set;


  error.init();

  if (variance == 0.) {
    plot_set = NULL;
    error.update(STAT_error[STATR_PLOT_NULL_VARIANCE]);
  }

  else {
    int i , j;
    int xmax;
    Curves *survival_rate;
    ostringstream legend;


    plot_set = new MultiPlotSet(3);
    MultiPlotSet &plot = *plot_set;

    plot.title = "Survival analysis";
    plot.border = "15 lw 0";

    // frequency distribution

    plot[0].xrange = Range(0 , nb_value - 1);
    plot[0].yrange = Range(0 , ceil(max * YSCALE));

    if (nb_value - 1 < TIC_THRESHOLD) {
      plot[0].xtics = 1;
    }

    plot[0].resize(1);

    plot[0][0].legend = STAT_label[STATL_FREQUENCY_DISTRIBUTION];

    plot[0][0].style = "impulses";

    plotable_frequency_write(plot[0][0]);

    // probability mass and survivor functions

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

    // survival rates

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


/*--------------------------------------------------------------*/
/**
 *  \brief Computation de la cumulative distribution function deduced from
 *         a frequency distribution.
 *
 *  \param[in] scale scaling factor.
 *
 *  \return          cumulative distribution function.
 */
/*--------------------------------------------------------------*/

double* FrequencyDistribution::cumul_computation(double scale) const

{
  int i;
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


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the survivor function deduced from a frequency distribution.
 *
 *  \param[in] scale scaling factor.
 *
 *  \return          survivor function.
 */
/*--------------------------------------------------------------*/

double* FrequencyDistribution::survivor_function_computation(double scale) const

{
  int i;
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


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the concentration function deduced from a frequency distribution.
 *
 *  \param[in] scale scaling factor.
 *
 *  \return          concentration function.
 */
/*--------------------------------------------------------------*/

double* FrequencyDistribution::concentration_function_computation(double scale) const

{
  int i;
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


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the coefficient of concentration of a frequency distribution.
 *
 *  \return coefficient of concentration.
 */
/*--------------------------------------------------------------*/

double FrequencyDistribution::concentration_computation() const

{
  int i;
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


/*--------------------------------------------------------------*/
/**
 *  \brief Update of an integer frequency distribution from
 *         a real frequency distribution by rounding.
 *
 *  \param[in] reestim     pointer on the real frequency distribution,
 *  \param[in] inb_element cumulative frequencies.
 */
/*--------------------------------------------------------------*/

void FrequencyDistribution::update(const Reestimation<double> *reestim , int inb_element)

{
  int i , j;
  int index;
  double scale , sum , max_frequency , *real_frequency;


  // copy of the real frequency distribution and scaling

  real_frequency = new double[reestim->nb_value];

  scale = inb_element / reestim->nb_element;
  for (i = reestim->offset;i < reestim->nb_value;i++) {
    real_frequency[i] = reestim->frequency[i] * scale;
  }

  // computation of frequencies

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

  // rounding

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

  // computation of the frequency distribution characteristics

  nb_value_computation();
  offset_computation();
  nb_element_computation();
  max_computation();
  mean_computation();
  variance_computation();

  delete [] real_frequency;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of a frequency distribution from an initial frequency distribution
 *         by scaling frequencies.
 *
 *  \param[in] inb_element cumulative frequencies.
 *
 *  \return                scaled frequency distribution.
 */
/*--------------------------------------------------------------*/

FrequencyDistribution* FrequencyDistribution::frequency_scale(int inb_element) const

{
  int i , j;
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

  // rounding

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

  // computation of the frequency distribution characteristics

  histo->nb_value_computation();
  histo->offset_computation();
  histo->nb_element = inb_element;
  histo->max_computation();
  histo->mean_computation();
  histo->variance_computation();

  return histo;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of ranks of values on the basis of a frequency distribution.
 *
 *  \return ranks.
 */
/*--------------------------------------------------------------*/

double* FrequencyDistribution::rank_computation() const

{
  int i;
  double *rank;


  rank = new double[nb_value];

  rank[offset] = (double)(1 + frequency[offset]) / 2.;
  for (i = offset + 1;i < nb_value;i++) {
    rank[i] = rank[i - 1] + (double)(frequency[i - 1] + frequency[i]) / 2.;
  }

  return rank;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the cumulative frequency distribution function.
 *
 *  \param[in] cdf (values, cumulative frequency distribution function).
 *
 *  \return        number of values.
 */
/*--------------------------------------------------------------*/

int FrequencyDistribution::cumulative_distribution_function_computation(double **cdf) const

{
  int i , j;
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


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the minimum interval between 2 values.
 *
 *  \return minimum interval.
 */
/*--------------------------------------------------------------*/

int FrequencyDistribution::min_interval_computation() const

{
  int i;
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


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the log-likelihood of a continuous distribution for a frequency distribution.
 *
 *  \param[in] dist         reference on a ContinuousParametric object,
 *  \param[in] min_interval minimum interval between 2 values.
 *
 *  \return                 log-likelihood.
 */
/*--------------------------------------------------------------*/

double FrequencyDistribution::likelihood_computation(const ContinuousParametric &dist ,
                                                     int min_interval) const

{
  int i;
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
