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



#include <math.h>

#include <sstream>
#include <iomanip>
#include <cstring>

#include "stat_tools.h"
#include "curves.h"
#include "stat_label.h"

using namespace std;


namespace stat_tool {



/*--------------------------------------------------------------*/
/**
 *  \brief Copy of the probability mass function.
 *
 *  \param[in] dist      reference on a Distribution object,
 *  \param[in] inb_value number of values.
 */
/*--------------------------------------------------------------*/

void Distribution::mass_copy(const Distribution &dist , int inb_value)

{
  int i;


  if ((inb_value != I_DEFAULT) && (inb_value < dist.nb_value)) {
    nb_value = inb_value;
  }
  else {
    nb_value = dist.nb_value;
  }
  offset = MIN(dist.offset , nb_value - 1);

  for (i = 0;i < nb_value;i++) {
    mass[i] = dist.mass[i];
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Copy of a Distribution object in the case where the number of
 *         allocated values is the same for the 2 Distribution objects.
 *
 *  \param[in] dist reference on a Distribution object.
 */
/*--------------------------------------------------------------*/

void Distribution::equal_size_copy(const Distribution &dist)

{
  int i;


  nb_value = dist.nb_value;

  offset = dist.offset;
  max = dist.max;
  complement = dist.complement;
  mean = dist.mean;
  variance = dist.variance;
  nb_parameter = dist.nb_parameter;

  for (i = 0;i < nb_value;i++) {
    mass[i] = dist.mass[i];
    cumul[i] = dist.cumul[i];
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Initialization of a Distribution object.
 *
 *  \param[in] inb_value number of values.
 */
/*--------------------------------------------------------------*/

void Distribution::init(int inb_value)

{
  nb_value = inb_value;
  alloc_nb_value = nb_value;

  offset = 0;
  max = 0.;
  complement = 0.;
  mean = D_DEFAULT;
  variance = D_DEFAULT;
  nb_parameter = 0;

  if (nb_value == 0) {
    mass = NULL;
    cumul = NULL;
  }

  else {
    int i;


    mass = new double[nb_value];
    cumul = new double[nb_value];

    for (i = 0;i < nb_value;i++) {
      mass[i] = 0.;
      cumul[i] = 0.;
    }
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Constructor of the Distribution class.
 *
 *  \param[in] inb_value number of values.
 */
/*--------------------------------------------------------------*/

Distribution::Distribution(int inb_value)

{
  init(inb_value);
}


/*--------------------------------------------------------------*/
/**
 *  \brief Constructor of the Distribution class.
 *
 *  \param[in] inb_value number of values,
 *  \param[in] imass     probability mass function.
 */
/*--------------------------------------------------------------*/

Distribution::Distribution(int inb_value , double *imass)

{
  int i;


  nb_value = inb_value;
  alloc_nb_value = nb_value;

  complement = 0.;
  nb_parameter = 0;

  mass = new double[nb_value];
  cumul = new double[nb_value];

  for (i = 0;i < nb_value;i++) {
    mass[i] = imass[i];
  }

  cumul_computation();

  offset_computation();
  max_computation();
  mean_computation();
  variance_computation();
}


/*--------------------------------------------------------------*/
/**
 *  \brief Construction of a Distribution object from an initial
 *         Distribution object applying a scaling operation.
 *
 *  \param[in] dist          reference on a Distribution object,
 *  \param[in] scaling_coeff scaling factor.
 */
/*--------------------------------------------------------------*/

Distribution::Distribution(const Distribution &dist , double scaling_coeff)

{
  int i , j;
  int min , max;


  nb_value = (int)floor(dist.nb_value * scaling_coeff) + 1;
  alloc_nb_value = nb_value;

  offset = (int)floor(dist.offset * scaling_coeff);
  complement = dist.complement;
  nb_parameter = 0;

  mass = new double[nb_value];
  cumul = new double[nb_value];

  for (i = 0;i < nb_value;i++) {
    mass[i] = 0.;
  }

  for (i = dist.offset;i < dist.nb_value;i++) {
    min = (int)floor(i * scaling_coeff);
    max = (int)floor((i + 1) * scaling_coeff);

    if (scaling_coeff < 1.) {
      mass[min] += (max - i * scaling_coeff) * dist.mass[i] / scaling_coeff;
      mass[max] += ((i + 1) * scaling_coeff - max) * dist.mass[i] / scaling_coeff;
    }

    else {
      mass[min] += (min + 1 - i * scaling_coeff) * dist.mass[i] / scaling_coeff;
      for (j = min + 1;j < max;j++) {
        mass[j] += dist.mass[i] / scaling_coeff;
      }
      mass[max] += ((i + 1) * scaling_coeff - max) * dist.mass[i] / scaling_coeff;
    }
  }

  cumul_computation();

  max_computation();
  mean_computation();
  variance_computation();
}


/*--------------------------------------------------------------*/
/**
 *  \brief Construction of a Distribution object from a FrequencyDistribution object.
 *
 *  \param[in] histo reference on a FrequencyDistribution object.
 */
/*--------------------------------------------------------------*/

Distribution::Distribution(const FrequencyDistribution &histo)

{
  nb_value = histo.nb_value;
  alloc_nb_value = nb_value;

  complement = 0.;
  nb_parameter = 0;

  mass = new double[nb_value];
  cumul = new double[nb_value];

  histo.distribution_estimation(this);
}


/*--------------------------------------------------------------*/
/**
 *  \brief Copy of a Distribution object.
 *
 *  \param[in] dist            reference on a Distribution object,
 *  \param[in] ialloc_nb_value number of allocated values.
 */
/*--------------------------------------------------------------*/

void Distribution::copy(const Distribution &dist , int ialloc_nb_value)

{
  int i;


  nb_value = dist.nb_value;
  if (ialloc_nb_value == I_DEFAULT) {
    alloc_nb_value = nb_value;
  }
  else {
    alloc_nb_value = MAX(nb_value , ialloc_nb_value);
  }

  offset = dist.offset;
  max = dist.max;
  complement = dist.complement;
  mean = dist.mean;
  variance = dist.variance;
  nb_parameter = 0;

  mass = new double[alloc_nb_value];
  cumul = new double[alloc_nb_value];

  for (i = 0;i < nb_value;i++) {
    mass[i] = dist.mass[i];
    cumul[i] = dist.cumul[i];
  }
  for (i = nb_value;i < alloc_nb_value;i++) {
    mass[i] = 0.;
    cumul[i] = 0.;
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Copy of a Distribution object with renormalization.
 *
 *  \param[in] dist reference on a Distribution object.
 */
/*--------------------------------------------------------------*/

void Distribution::normalization_copy(const Distribution &dist)

{
  if (dist.complement == 0.) {
    copy(dist);
  }

  else {
    int i;


    nb_value = dist.nb_value;
    alloc_nb_value = nb_value;

    offset = dist.offset;
    complement = 0.;
    nb_parameter = 0;

    mass = new double[nb_value];
    cumul = new double[nb_value];

    for (i = 0;i < nb_value;i++) {
      mass[i] = dist.mass[i] / (1. - dist.complement);
    }

    cumul_computation();

    max = dist.max / (1. - dist.complement);
    mean_computation();
    variance_computation();
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Copy constructor of the Distribution class.
 *
 *  \param[in] dist            reference on a Distribution object,
 *  \param[in] transform       type of transform (DISTRIBUTION_COPY/NORMALIZATION),
 *  \param[in] ialloc_nb_value number of allocated values.
 */
/*--------------------------------------------------------------*/

Distribution::Distribution(const Distribution &dist , distribution_transformation transform ,
                           int ialloc_nb_value)

{
  switch (transform) {
  case DISTRIBUTION_COPY :
    copy(dist , ialloc_nb_value);
    break;
  case NORMALIZATION :
    normalization_copy(dist);
    break;
  default :
    copy(dist);
    break;
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Destructor of the Distribution class.
 */
/*--------------------------------------------------------------*/

Distribution::~Distribution()

{
  delete [] mass;
  delete [] cumul;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Assignment operator of the Distribution class.
 *
 *  \param[in] dist reference on a Distribution object.
 *
 *  \return         Distribution object.
 */
/*--------------------------------------------------------------*/

Distribution& Distribution::operator=(const Distribution &dist)

{
  if (&dist != this) {
    delete [] mass;
    delete [] cumul;

    copy(dist);
  }

  return *this;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Equality operator of the Distribution class.
 *
 *  \param[in] dist reference on a Distribution object.
 *
 *  \return         equality or not of the discrete distributions.
 */
/*--------------------------------------------------------------*/

bool Distribution::operator==(const Distribution &dist) const

{
  bool status = true;
  int i;


  if ((offset != dist.offset) || (nb_value != dist.nb_value)) {
    status = false;
  }

  else {
    for (i = offset;i < nb_value;i++) {
      if (mass[i] != dist.mass[i]) {
        status = false;
        break;
      }
    }
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of the characteristics of a discrete distribution.
 *
 *  \param[in,out] os           stream,
 *  \param[in]     shape        flag on the writing of the shape characteristics,
 *  \param[in]     comment_flag flag comment.
 */
/*--------------------------------------------------------------*/

ostream& Distribution::ascii_characteristic_print(ostream &os , bool shape , bool comment_flag) const

{
  if ((mean != D_DEFAULT) && (variance != D_DEFAULT)) {
    if (comment_flag) {
      os << "# ";
    }
    os << STAT_label[STATL_MEAN] << ": " << mean << "   ";
    if (complement == 0.) {
      os << STAT_label[STATL_MEDIAN] << ": " << quantile_computation() << "   ";
    }
    os << STAT_label[STATL_MODE] << ": " << mode_computation() << endl;

    if (comment_flag) {
      os << "# ";
    }
    os << STAT_label[STATL_VARIANCE] << ": " << variance << "   "
       << STAT_label[STATL_STANDARD_DEVIATION] << ": " << sqrt(variance);
    if ((complement == 0.) && (variance > 0.)) {
      os << "   " << STAT_label[STATL_LOWER_QUARTILE] << ": " << quantile_computation(0.25)
         << "   " << STAT_label[STATL_UPPER_QUARTILE] << ": " << quantile_computation(0.75);
    }
    os << endl;

    if ((shape) && (variance > 0.)) {
      if (comment_flag) {
        os << "# ";
      }
      os << STAT_label[STATL_SKEWNESS_COEFF] << ": " << skewness_computation() << "   "
         << STAT_label[STATL_KURTOSIS_COEFF] << ": " << kurtosis_computation() << endl;
    }
  }

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the column width for an integer value.
 *
 *  \param[in] value value.
 *
 *  \return          column width.
 */
/*--------------------------------------------------------------*/

int column_width(int value)

{
  ostringstream ostring;
  ostring << value;
  return (ostring.str()).size();
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the column width for integer values.
 *
 *  \param[in] min_value minimum value,
 *  \param[in] max_value maximum value.
 *
 *  \return              column width.
 */
/*--------------------------------------------------------------*/

int column_width(int min_value , int max_value)

{
  int width , max_width;


  {
    ostringstream ostring;
    ostring << min_value;
    max_width = (ostring.str()).size();
  }

  {
    ostringstream ostring;
    ostring << max_value;
    width = (ostring.str()).size();
    if (width > max_width) {
      max_width = width;
    }
  }

  return max_width;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the column width for real values.
 *
 *  \param[in] nb_value number of values,
 *  \param[in] value    pointer on the values,
 *  \param[in] scale    scaling factor.
 *
 *  \return             column width.
 */
/*--------------------------------------------------------------*/

int column_width(int nb_value , const double *value , double scale)

{
  int i;
  int width , max_width = 0;


  for (i = 0;i < nb_value;i++) {
    ostringstream ostring;
    ostring << value[i] * scale;
    width = (ostring.str()).size();
    if (width > max_width) {
      max_width = width;
    }
  }

  return max_width;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a discrete distribution and a frequency distribution.
 *
 *  \param[in,out] os            stream,
 *  \param[in]     comment_flag  flag comment,
 *  \param[in]     cumul_flag    flag on the writing of the cumulative distribution functions,
 *  \param[in]     nb_value_flag flag on the computation of the number of values,
 *  \param[in]     histo         pointer on a frequency distribution.
 */
/*--------------------------------------------------------------*/

ostream& Distribution::ascii_print(ostream &os , bool comment_flag , bool cumul_flag ,
                                   bool nb_value_flag , const FrequencyDistribution *histo) const

{
  int i;
  int ascii_nb_value = nb_value , width[5];
  double scale , *histo_cumul , *pcumul;
  ios_base::fmtflags format_flags;


  format_flags = os.setf(ios::right , ios::adjustfield);

  // computation of the scaling factor, the cumulative distribution function deduced from
  // the frequency distribution and the number of values

  if (histo) {
    scale = histo->nb_element / (1. - complement);
    if (cumul_flag) {
      histo_cumul = histo->cumul_computation(scale);
    }
  }

  else {
    scale = 1.;
  }

  if ((nb_value_flag) && (nb_value >= 2)) {
    pcumul = cumul + nb_value - 1;
    while (*--pcumul > 1. - complement - ASCII_ROUNDNESS) {
      ascii_nb_value--;
    }

    if ((histo) && (histo->nb_value > ascii_nb_value)) {
      ascii_nb_value = histo->nb_value;
    }
  }

  // computation of column widths

  width[0] = column_width(ascii_nb_value - 1);
  if (histo) {
    width[1] = column_width(histo->max) + ASCII_SPACE;
  }
  width[2] = column_width(ascii_nb_value , mass , scale) + ASCII_SPACE;
  if (cumul_flag) {
    if (histo) {
      width[3] = column_width(histo->nb_value , histo_cumul , 1.) + ASCII_SPACE;
    }
    width[4] = column_width(ascii_nb_value , cumul , 1.) + ASCII_SPACE;
  }

  // writing of the probability mass function or of the empirical and theoretical frequency distributions

  for (i = 0;i < ascii_nb_value;i++) {
    if (comment_flag) {
      os << "# ";
    }
    os << setw(width[0]) << i;

    if (histo) {
      if (i < histo->nb_value) {
        os << setw(width[1]) << histo->frequency[i];
      }
      else {
        os << setw(width[1]) << " ";
      }
    }

    os << setw(width[2]) << mass[i] * scale;

   // writing of the empirical and theoretical cumulative distribution functions

    if (cumul_flag) {
      if (histo) {
        if (i < histo->nb_value) {
          os << setw(width[3]) << histo_cumul[i];
        }
        else {
          os << setw(width[3]) << " ";
        }
      }

      os << setw(width[4]) << cumul[i];
    }
    os << endl;
  }

  if ((histo) && (cumul_flag)) {
    delete [] histo_cumul;
  }

  os.setf(format_flags , ios::adjustfield);

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a family of discrete distributions and a frequency distribution.
 *
 *  \param[in,out] os           stream,
 *  \param[in]     nb_dist      number of distributions,
 *  \param[in]     dist         pointer on the discrete distributions,
 *  \param[in]     dist_scale   scaling factors,
 *  \param[in]     comment_flag flag comment,
 *  \param[in]     cumul_flag   flag on the writing of the cumulative distribution functions,
 *  \param[in]     histo        pointer on a frequency distribution,
 *  \param[in]     mass_first   flag on the distribution order.
 */
/*--------------------------------------------------------------*/

ostream& Distribution::ascii_print(ostream &os , int nb_dist , const Distribution **dist ,
                                   double *dist_scale , bool comment_flag , bool cumul_flag ,
                                   const FrequencyDistribution *histo , bool mass_first) const

{
  int i , j;
  int ascii_nb_value = nb_value , *width;
  double scale , *histo_cumul;
  ios_base::fmtflags format_flags;


  format_flags = os.setf(ios::right , ios::adjustfield);

  width = new int[nb_dist + 5];

  // computation of the scaling factor, the cumulative distribution function deduced from
  // the frequency distribution and the number of values

  if (histo) {
    scale = histo->nb_element / (1. - complement);
    if (cumul_flag) {
      histo_cumul = histo->cumul_computation(scale);
    }
  }

  else {
    scale = 1.;
  }

  for (i = 0;i < nb_dist;i++) {
    if (dist[i]->nb_value > ascii_nb_value) {
      ascii_nb_value = dist[i]->nb_value;
    }
  }

  // computation of the column widths

  width[0] = column_width(ascii_nb_value - 1);
  if (histo) {
    width[1] = column_width(histo->max) + ASCII_SPACE;
  }
  for (i = 0;i < nb_dist;i++) {
    width[i + 2] = column_width(dist[i]->nb_value , dist[i]->mass , dist_scale[i]) + ASCII_SPACE;
  }
  width[nb_dist + 2] = column_width(nb_value , mass , scale) + ASCII_SPACE;
  if (cumul_flag) {
    if (histo) {
      width[nb_dist + 3] = column_width(histo->nb_value , histo_cumul , 1.) + ASCII_SPACE;
    }
    width[nb_dist + 4] = column_width(nb_value , cumul , 1.) + ASCII_SPACE;
  }

  // writing of the probability mass functions or the empirical and theoretical frequency distributions and
  // the empirical and theoretical cumulative distribution functions

  for (i = 0;i < ascii_nb_value;i++) {
    if (comment_flag) {
      os << "# ";
    }
    os << setw(width[0]) << i;

    if (histo) {
      if (i < histo->nb_value) {
        os << setw(width[1]) << histo->frequency[i];
      }
      else {
        os << setw(width[1]) << " ";
      }
    }

    if (mass_first) {
      if (i < nb_value) {
        os << setw(width[nb_dist + 2]) << mass[i] * scale;
      }
      else {
        os << setw(width[nb_dist + 2]) << " ";
      }
    }

    for (j = 0;j < nb_dist;j++) {
      if (i < dist[j]->nb_value) {
        os << setw(width[j + 2]) << dist[j]->mass[i] * dist_scale[j];
      }
      else {
        os << setw(width[j + 2]) << " ";
      }
    }

    if (!mass_first) {
      if (i < nb_value) {
        os << setw(width[nb_dist + 2]) << mass[i] * scale;
      }
      else {
        os << setw(width[nb_dist + 2]) << " ";
      }
    }

    if (cumul_flag) {
      if (histo) {
        if (i < histo->nb_value) {
          os << setw(width[nb_dist + 3]) << histo_cumul[i];
        }
        else {
          os << setw(width[nb_dist + 3]) << " ";
        }
      }

      if (i < nb_value) {
        os << setw(width[nb_dist + 4]) << cumul[i];
      }
    }
    os << endl;
  }

  delete [] width;
  if ((histo) && (cumul_flag)) {
    delete [] histo_cumul;
  }

  os.setf(format_flags , ios::adjustfield);

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of the characteristics of a discrete distribution at the spreadsheet format.
 *
 *  \param[in,out] os    stream,
 *  \param[in]     shape flag on the writing of the shape characteristics.
 */
/*--------------------------------------------------------------*/

ostream& Distribution::spreadsheet_characteristic_print(ostream &os , bool shape) const

{
  if ((mean != D_DEFAULT) && (variance != D_DEFAULT)) {
    os << STAT_label[STATL_MEAN] << "\t" << mean << "\t\t";
    if (complement == 0.) {
      os << STAT_label[STATL_MEDIAN] << "\t" << quantile_computation() << "\t\t";
    }
    os << STAT_label[STATL_MODE] << "\t" << mode_computation() << endl;

    os << STAT_label[STATL_VARIANCE] << "\t" << variance << "\t\t"
       << STAT_label[STATL_STANDARD_DEVIATION] << "\t" << sqrt(variance);
    if ((complement == 0.) && (variance > 0.)) {
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
 *  \brief Writing of a discrete distribution and a frequency distribution at the spreadsheet format.
 *
 *  \param[in,out] os                 stream,
 *  \param[in]     cumul_flag         flag on the writing of the cumulative distribution functions,
 *  \param[in]     concentration_flag flag on the writing of the concentration functions,
 *  \param[in]     nb_value_flag      flag on the computation of the number of values,
 *  \param[in]     histo              pointer on a frequency distribution.
 */
/*--------------------------------------------------------------*/

ostream& Distribution::spreadsheet_print(ostream &os , bool cumul_flag ,
                                         bool concentration_flag , bool nb_value_flag ,
                                         const FrequencyDistribution *histo) const

{
  int i;
  int spreadsheet_nb_value = nb_value;
  double scale , *pcumul , *histo_cumul , *concentration , *histo_concentration;


  if ((!cumul_flag) || (mean == 0.)) {
    concentration_flag = false;
  }

  // computation of the scaling factor, the cumulative distribution function deduced from
  // the frequency distribution, the concentration functions and the number of values

  if (histo) {
    scale = histo->nb_element / (1. - complement);
    if (cumul_flag) {
      histo_cumul = histo->cumul_computation(scale);
    }
    if ((concentration_flag) && (histo->mean > 0.)) {
      histo_concentration = histo->concentration_function_computation(scale);
    }
  }

  else {
    scale = 1.;
  }

  if (concentration_flag) {
    concentration = concentration_function_computation();
  }

  if ((nb_value_flag) && (nb_value >= 2)) {
    pcumul = cumul + nb_value - 1;
    while (*--pcumul > 1. - complement - SPREADSHEET_ROUNDNESS) {
      spreadsheet_nb_value--;
    }

    if ((histo) && (histo->nb_value > spreadsheet_nb_value)) {
      spreadsheet_nb_value = histo->nb_value;
    }
  }

  // writing of the probability mass function or the empirical and theoretical frequency distributions and
  // the empirical and theoretical cumulative distribution and concentration functions

  for (i = 0;i < spreadsheet_nb_value;i++) {
    os << i;

    if (histo) {
      os << "\t";
      if (i < histo->nb_value) {
        os << histo->frequency[i];
      }
    }

    os << "\t" << mass[i] * scale;

    if (cumul_flag) {
      if (histo) {
        os << "\t";
        if (i < histo->nb_value) {
          os << histo_cumul[i];
        }
      }

      os << "\t" << cumul[i];
    }

    if (concentration_flag) {
      if ((histo) && (histo->mean > 0.)) {
        os << "\t";
        if (i < histo->nb_value) {
          os << histo_concentration[i];
        }
      }

      os << "\t" << concentration[i];
    }
    os << endl;
  }

  if (histo) {
    if (cumul_flag) {
      delete [] histo_cumul;
    }
    if ((concentration_flag) && (histo->mean > 0.)) {
      delete [] histo_concentration;
    }
  }

  if (concentration_flag) {
    delete [] concentration;
  }

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a family of discrete distributions and a frequency distribution
 *         at the spreadsheet format.
 *
 *  \param[in,out] os         stream,
 *  \param[in]     nb_dist    number of distributions,
 *  \param[in]     dist       pointer on the discrete distributions,
 *  \param[in]     dist_scale scaling factors,
 *  \param[in]     cumul_flag flag on the writing of the cumulative distribution functions,
 *  \param[in]     histo      pointer on a frequency distribution,
 *  \param[in]     mass_first flag on the distribution order.
 */
/*--------------------------------------------------------------*/

ostream& Distribution::spreadsheet_print(ostream &os , int nb_dist , const Distribution **dist ,
                                         double *dist_scale , bool cumul_flag ,
                                         const FrequencyDistribution *histo , bool mass_first) const

{
  int i , j;
  int spreadsheet_nb_value = nb_value;
  double scale , *histo_cumul;


  // computation of the scaling factor, the cumulative distribution function deduced from
  // the frequency distribution and the number of values

  if (histo) {
    scale = histo->nb_element / (1. - complement);
    if (cumul_flag) {
      histo_cumul = histo->cumul_computation(scale);
    }
  }

  else {
    scale = 1.;
  }

  for (i = 0;i < nb_dist;i++) {
    if (dist[i]->nb_value > spreadsheet_nb_value) {
      spreadsheet_nb_value = dist[i]->nb_value;
    }
  }

  // writing of the probability mass functions or the empirical and theoretical frequency distributions and
  // the empirical and theoretical cumulative distribution functions

  for (i = 0;i < spreadsheet_nb_value;i++) {
    os << i;

    if (histo) {
      os << "\t";
      if (i < histo->nb_value) {
        os << histo->frequency[i];
      }
    }

    if (mass_first) {
      os << "\t";
      if (i < nb_value) {
        os << mass[i] * scale;
      }
    }

    for (j = 0;j < nb_dist;j++) {
      os << "\t";
      if (i < dist[j]->nb_value) {
        os << dist[j]->mass[i] * dist_scale[j];
      }
    }

    if (!mass_first) {
      os << "\t";
      if (i < nb_value) {
        os << mass[i] * scale;
      }
    }

    if (cumul_flag) {
      if (histo) {
        os << "\t";
        if (i < histo->nb_value) {
          os << histo_cumul[i];
        }
      }

      if (i < nb_value) {
        os << "\t" << cumul[i];
      }
    }
    os << endl;
  }

  if ((histo) && (cumul_flag)) {
    delete [] histo_cumul;
  }

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the number of plotted values (Gnuplot output).
 *
 *  \param[in] histo pointer a frequency distribution.
 *
 *  \return          number of plotted values.
 */
/*--------------------------------------------------------------*/

int Distribution::plot_nb_value_computation(const FrequencyDistribution *histo) const

{
  int plot_nb_value = nb_value;
  double *pcumul;


  if (nb_value >= 2) {
    pcumul = cumul + nb_value - 1;
    while ((*--pcumul > 1. - complement - PLOT_ROUNDNESS) && (plot_nb_value > 2)) {
      plot_nb_value--;
    }

    if ((histo) && (histo->nb_value > plot_nb_value)) {
      plot_nb_value = histo->nb_value;
    }
  }

  return plot_nb_value;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a discrete distribution at the Gnuplot format.
 *
 *  \param[in] path          file path,
 *  \param[in] concentration pointer on the concentration function,
 *  \param[in] scale         scaling factor.
 *
 *  \return                  error status.
 */
/*--------------------------------------------------------------*/

bool Distribution::plot_print(const char *path , double *concentration ,
                              double scale) const

{
  bool status = false;
  int i;
  ofstream out_file(path);


  if (out_file) {
    status = true;

    if ((offset == 0) && (variance > 0.)) {
      out_file << -1 << " " << 0 << " " << 0 << " " << 0 << endl;
    }

    for (i = 0;i < nb_value;i++) {
      out_file << i << " " << mass[i] * scale;
      if (variance > 0.) {
        out_file << " " << cumul[i] << " " << concentration[i];
      }
      out_file << endl;
    }

    if (nb_value == 1) {
      out_file << nb_value << " " << 0 << endl;
    }
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a discrete distribution and a frequency distribution
 *         at the Gnuplot format.
 *
 *  \param[in] path  file path,
 *  \param[in] histo pointer on a frequency distribution.
 *
 *  \return          error status.
 */
/*--------------------------------------------------------------*/

bool Distribution::plot_print(const char *path , const FrequencyDistribution *histo) const

{
  bool status = false;
  int i;
  double scale;
  ofstream out_file(path);


  if (out_file) {
    status = true;

    // computation of the scaling factor

    if (histo) {
      scale = histo->nb_element / (1. - complement);
    }
    else {
      scale = 1.;
    }

    // writing of the probability mass function or the empirical and theoretical frequency distributions

    for (i = 0;i < nb_value;i++) {
      if (histo) {
        if (i < histo->nb_value) {
          out_file << histo->frequency[i] << " ";
        }
        else {
          out_file << 0 << " ";
        }
      }

      out_file << mass[i] * scale << endl;
    }
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a family of a discrete distributions and a family of
 *         frequency distributions at the Gnuplot format.
 *
 *  \param[in] path          file path,
 *  \param[in] nb_dist       number of distributions,
 *  \param[in] dist          pointer on the discrete distributions,
 *  \param[in] scale         scaling factors,
 *  \param[in] dist_nb_value number of values of the distributions,
 *  \param[in] nb_histo      number of frequency distributions,
 *  \param[in] histo         pointer on the frequency distributions.
 *
 *  \return                  error status.
 */
/*--------------------------------------------------------------*/

bool plot_print(const char *path , int nb_dist , const Distribution **dist ,
                double *scale , int *dist_nb_value , int nb_histo ,
                const FrequencyDistribution **histo)

{
  bool status = false;
  int i , j;
  int plot_nb_value = 0;
  ofstream out_file(path);


  if (out_file) {
    status = true;

    // computation of the number of values

    if (histo) {
      for (i = 0;i < nb_histo;i++) {
        if (histo[i]->nb_value > plot_nb_value) {
          plot_nb_value = histo[i]->nb_value;
        }
      }
    }

    if (!dist_nb_value) {
      for (i = 0;i < nb_dist;i++) {
        if (dist[i]->nb_value > plot_nb_value) {
          plot_nb_value = dist[i]->nb_value;
        }
      }
    }

    else {
      for (i = 0;i < nb_dist;i++) {
        if (dist_nb_value[i] > plot_nb_value) {
          plot_nb_value = dist_nb_value[i];
        }
      }
    }

    // writing of the probability mass functions or the empirical and theoretical frequency distributions

    for (i = 0;i < plot_nb_value;i++) {
      if (histo) {
        for (j = 0;j < nb_histo;j++) {
          if (i < histo[j]->nb_value) {
            out_file << histo[j]->frequency[i] << " ";
          }
          else {
            out_file << 0 << " ";
          }
        }
      }

      if (!dist_nb_value) {
        for (j = 0;j < nb_dist;j++) {
          if (i < dist[j]->nb_value) {
            out_file << dist[j]->mass[i] * scale[j] << " ";
          }
          else {
            out_file << 0. << " ";
          }
        }
      }

      else {
        for (j = 0;j < nb_dist;j++) {
          if (i < dist_nb_value[j]) {
            out_file << dist[j]->mass[i] * scale[j] << " ";
          }
          else {
            out_file << 0. << " ";
          }
        }
      }

      out_file << endl;
    }
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of family of cumulative distribution functions for their matching
 *         at the Gnuplot format
 *
 *  \param[in] path     file path,
 *  \param[in] nb_cumul number of cumulative distribution functions,
 *  \param[in] offset   pointer on the lowest values,
 *  \param[in] nb_value pointer on the number of values,
 *  \param[in] cumul    pointer on the cumulative distribution functions.
 *
 *  \return             error status.
 */
/*--------------------------------------------------------------*/

bool cumul_matching_plot_print(const char *path , int nb_cumul , int *offset ,
                               int *nb_value , double **cumul)

{
  bool status = false;
  int i , j;
  int plot_offset , plot_nb_value;
  ofstream out_file(path);


  if (out_file) {
    status = true;

    plot_offset = offset[0];
    plot_nb_value = nb_value[0];
    for (i = 1;i < nb_cumul;i++) {
      if (offset[i] < plot_offset) {
        plot_offset = offset[i];
      }
      if (nb_value[i] > plot_nb_value) {
        plot_nb_value = nb_value[i];
      }
    }

    for (i = 0;i < nb_cumul;i++) {
      out_file << 0 << " ";
    }
    out_file << endl;

    for (i = plot_offset;i < plot_nb_value;i++) {
      for (j = 0;j < nb_cumul;j++) {
        if (i < nb_value[j]) {
          out_file << cumul[j][i] << " ";
        }
        else {
          out_file << cumul[j][nb_value[j] - 1] << " ";
        }
      }
      out_file << endl;
    }
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Extraction of a file name label (Gnuplot output).
 *
 *  \param[in] file_name file name.
 *
 *  \return              label.
 */
/*--------------------------------------------------------------*/

char* label(const char *file_name)

{
  char *pfile_name;


  if (*file_name == '/') {
    pfile_name = (char*)file_name;
  }

  else {
    pfile_name = (char*)strrchr(file_name , '/');
    if (pfile_name) {
      pfile_name++;
    }
    else {
      pfile_name = (char*)file_name;
    }
  }

  return pfile_name;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Plot of family of a discrete distributions using Gnuplot:
 *         - probability mass functions and cumulative distribution functions,
 *         - matching of cumulative distribution functions,
 *         - concentration curves.
 *
 *  \param[in] error   reference on a StatError object,
 *  \param[in] prefix  file prefix,
 *  \param[in] nb_dist number of distributions,
 *  \param[in] idist   pointer on the discrete distributions,
 *  \param[in] title   figure title.
 *
 *  \return            error status.
 */
/*--------------------------------------------------------------*/

bool Distribution::plot_write(StatError &error , const char *prefix , int nb_dist ,
                              const Distribution **idist , const char *title) const

{
  bool status;


  error.init();

  if (nb_dist > PLOT_NB_DISTRIBUTION) {
    status = false;
    error.update(STAT_error[STATR_PLOT_NB_DISTRIBUTION]);
  }

  else {
    bool cumul_concentration_flag;
    int i , j , k;
    int plot_nb_value , max_nb_value , max_range , reference_matching ,
        reference_concentration , *poffset , *pnb_value;
    double max , min_complement , **pcumul , **concentration;
    ostringstream *data_file_name;
    const Distribution **dist;


    nb_dist++;
    dist = new const Distribution*[nb_dist];

    dist[0] = this;
    for (i = 1;i < nb_dist;i++) {
      dist[i] = idist[i - 1];
    }

    // writing of the data files

    data_file_name = new ostringstream[nb_dist + 1];

    poffset = new int[nb_dist];
    pnb_value = new int[nb_dist];
    pcumul = new double*[nb_dist];

    concentration = new double*[nb_dist];
    for (i = 0;i < nb_dist;i++) {
      concentration[i] = NULL;
    }

    max_nb_value = 0;
    max = 0.;
    cumul_concentration_flag = false;
    max_range = 0;
    min_complement = 1.;

    for (i = 0;i < nb_dist;i++) {
      data_file_name[i] << prefix << i << ".dat";

      poffset[i] = dist[i]->offset;
      pnb_value[i] = dist[i]->nb_value;
      pcumul[i] = dist[i]->cumul;

      // computation of the concentration functions

      concentration[i] = dist[i]->concentration_function_computation();

      // computation of the maximum number of values, the largest support and
      // the maximum probability

      plot_nb_value = dist[i]->plot_nb_value_computation();
      if (plot_nb_value > max_nb_value) {
        max_nb_value = plot_nb_value;
      }
      if (dist[i]->max > max) {
        max = dist[i]->max;
      }

      if (dist[i]->variance > 0.) {
        cumul_concentration_flag = true;
        if (dist[i]->nb_value - dist[i]->offset > max_range) {
          max_range = dist[i]->nb_value - dist[i]->offset;
          reference_matching = i;
        }
        if (dist[i]->complement < min_complement) {
          min_complement = dist[i]->complement;
          reference_concentration = i;
        }
      }

      status = dist[i]->plot_print((data_file_name[i].str()).c_str() , concentration[i] , 1.);

      if (!status) {
        break;
      }
    }

    if ((status) && (cumul_concentration_flag) && (nb_dist > 1)) {
      data_file_name[nb_dist] << prefix << nb_dist << ".dat";
      cumul_matching_plot_print((data_file_name[nb_dist].str()).c_str() , nb_dist ,
                                poffset , pnb_value , pcumul);
    }

    delete [] poffset;
    delete [] pnb_value;
    delete [] pcumul;

    for (i = 0;i < nb_dist;i++) {
      delete [] concentration[i];
    }
    delete [] concentration;

    if (status) {

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

        // probability mass functions

        if (MAX(max_nb_value , 2) - 1 < TIC_THRESHOLD) {
          out_file << "set xtics 0,1" << endl;
        }

        out_file << "plot [0:" << MAX(max_nb_value , 2) - 1 << "] [0:"
                 << MIN(max * YSCALE , 1.) << "] ";
        for (j = 0;j < nb_dist;j++) {
          out_file << "\"" << label((data_file_name[j].str()).c_str()) << "\" using 1:2 title \""
                   << STAT_label[STATL_DISTRIBUTION];
          if (nb_dist > 1) {
            out_file << " " << j + 1;
          }
          dist[j]->plot_title_print(out_file);
          out_file << "\" with linespoints";
          if (j < nb_dist - 1) {
            out_file << ",\\";
          }
          out_file << endl;
        }

        if (MAX(max_nb_value , 2) - 1 < TIC_THRESHOLD) {
          out_file << "set xtics autofreq" << endl;
        }

        if (cumul_concentration_flag) {

          // cumulative distribution functions

          if (i == 0) {
            out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
          }
          out_file << endl;

          if (max_nb_value - 1 < TIC_THRESHOLD) {
            out_file << "set xtics 0,1" << endl;
          }

          out_file << "plot [0:" << max_nb_value - 1 << "] [0:" << 1. - min_complement << "] ";
          j = 0;
          for (k = 0;k < nb_dist;k++) {
            if (dist[k]->variance > 0.) {
              if (j > 0) {
                out_file << ",\\\n";
              }
              j++;
              out_file << "\"" << label((data_file_name[k].str()).c_str()) << "\" using 1:3 title \""
                       << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_DISTRIBUTION];
              if (nb_dist > 1) {
                out_file << " " << k + 1;
              }
              out_file << " " << STAT_label[STATL_FUNCTION] << "\" with linespoints";
            }
          }
          out_file << endl;

          if (max_nb_value - 1 < TIC_THRESHOLD) {
            out_file << "set xtics autofreq" << endl;
          }

          // matching of the cumulative distribution functions taking as reference
          // the distribution with the largest support

          if (nb_dist > 1) {
            if (i == 0) {
              out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
            }
            out_file << endl;

            out_file << "set title";
            out_file << " \"";
            if (title) {
              out_file << title << " - ";
            }
            out_file << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_DISTRIBUTION]
                     << " " << STAT_label[STATL_FUNCTION] << " " << STAT_label[STATL_MATCHING] << "\"\n\n";

            out_file << "set grid\n" << "set xtics 0,0.1\n" << "set ytics 0,0.1" << endl;

            out_file << "plot [0:" << 1. - min_complement << "] [0:" << 1. - min_complement << "] ";
            j = 0;
            for (k = 0;k < nb_dist;k++) {
              if (dist[k]->variance > 0.) {
                if (j > 0) {
                  out_file << ",\\\n";
                }
                j++;
                out_file << "\"" << label((data_file_name[nb_dist].str()).c_str()) << "\" using "
                         << reference_matching + 1 << ":" << k + 1 << " title \""
                         << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_DISTRIBUTION] << " "
                         << k + 1 << " " << STAT_label[STATL_FUNCTION] << "\" with linespoints";
              }
            }
            out_file << endl;

            out_file << "unset grid\n" << "set xtics autofreq\n" << "set ytics autofreq" << endl;
          }

          // concentration curves

          if (i == 0) {
            out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
          }
          out_file << endl;

          out_file << "set title";
          if (title) {
            out_file << " \"" << title << "\"";
          }
          out_file << "\n\n";

          out_file << "set grid\n" << "set xtics 0,0.1\n" << "set ytics 0,0.1" << endl;

          out_file << "plot [0:" << 1. - min_complement << "] [0:" << 1. - min_complement << "] ";
          for (j = 0;j < nb_dist;j++) {
            if (dist[j]->variance > 0.) {
              out_file << "\"" << label((data_file_name[j].str()).c_str()) << "\" using 3:4 title \""
                       << STAT_label[STATL_CONCENTRATION] << " " << STAT_label[STATL_CURVE];
              if (nb_dist > 1) {
                out_file << " " << j + 1;
              }
              out_file << "\" with linespoints,\\" << endl;
            }
          }
          out_file << "\""<< label((data_file_name[reference_concentration].str()).c_str())
                   << "\" using 3:3 notitle with lines" << endl;

          out_file << "unset grid\n" << "set xtics autofreq\n" << "set ytics autofreq" << endl;
        }

        if (i == 1) {
          out_file << "\nset terminal x11" << endl;
        }

        out_file << "\npause 0 \"" << STAT_label[STATL_END] << "\"" << endl;
      }
    }

    delete [] dist;
    delete [] data_file_name;

    if (!status) {
      error.update(STAT_error[STATR_FILE_PREFIX]);
    }
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of the probability mass function.
 *
 *  \param[in] plot  reference on a SinglePlot object,
 *  \param[in] scale scaling factor (default value: 1).
 */
/*--------------------------------------------------------------*/

void Distribution::plotable_mass_write(SinglePlot &plot , double scale) const

{
  int i;


  for (i = MAX(offset - 1 , 0);i < nb_value;i++) {
    plot.add_point(i , mass[i] * scale);
  }
  if ((cumul[nb_value - 1] > 1. - DOUBLE_ERROR) &&
      (mass[nb_value - 1] > PLOT_MASS_THRESHOLD)) {
    plot.add_point(nb_value , 0.);
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of the cumulative discrete distribution.
 *
 *  \param[in] plot reference on a SinglePlot object.
 */
/*--------------------------------------------------------------*/

void Distribution::plotable_cumul_write(SinglePlot &plot) const

{
  int i;


  for (i = MAX(offset - 1 , 0);i < nb_value;i++) {
    plot.add_point(i , cumul[i]);
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of the matching of the cumulative distribution function with
 *         the cumulative distribution function of a reference distribution.
 *
 *  \param[in] plot           reference on a SinglePlot object
 *  \param[in] reference_dist reference on the reference distribution.
 */
/*--------------------------------------------------------------*/

void Distribution::plotable_cumul_matching_write(SinglePlot &plot ,
                                                 const Distribution &reference_dist) const

{
  int i;


  plot.add_point(0. , 0.);
  for (i = MIN(reference_dist.offset , 1);i < offset;i++) {
    plot.add_point(reference_dist.cumul[i] , 0.);
  }
  for (i = offset;i < nb_value;i++) {
    plot.add_point(reference_dist.cumul[i] , cumul[i]);
  }
  for (i = nb_value;i < reference_dist.nb_value;i++) {
    plot.add_point(reference_dist.cumul[i] , cumul[nb_value - 1]);
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of the concentration curve deduced from a discrete distribution.
 *
 *  \param[in] plot reference on a SinglePlot object.
 */
/*--------------------------------------------------------------*/

void Distribution::plotable_concentration_write(SinglePlot &plot) const

{
  int i;
  double *concentration;


  concentration = concentration_function_computation();

  plot.add_point(0. , 0.);
  for (i = offset;i < nb_value;i++) {
    plot.add_point(cumul[i] , concentration[i]);
  }

  delete [] concentration;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Plot of family of a discrete distributions:
 *         - probability mass functions and cumulative distribution functions,
 *         - matching of cumulative distribution functions,
 *         - concentration curves.
 *
 *  \param[in] error   reference on a StatError object,
 *  \param[in] nb_dist number of distributions,
 *  \param[in] idist   pointer on the discrete distributions.
 *
 *  \return            MultiPlotSet object.
 */
/*--------------------------------------------------------------*/

MultiPlotSet* Distribution::get_plotable_distributions(StatError &error , int nb_dist ,
                                                       const Distribution **idist) const

{
  MultiPlotSet *plot_set;


  error.init();

  if (nb_dist > PLOT_NB_DISTRIBUTION) {
    plot_set = NULL;
    error.update(STAT_error[STATR_PLOT_NB_DISTRIBUTION]);
  }

  else {
    int i , j , k;
    int plot_nb_value , xmax , max_range , cumul_concentration_nb_dist ,
        nb_plot_set , reference_matching;
    double ymax , min_complement;
    const Distribution **dist;
    ostringstream legend , title;


    nb_dist++;
    dist = new const Distribution*[nb_dist];

    dist[0] = this;
    for (i = 1;i < nb_dist;i++) {
      dist[i] = idist[i - 1];
    }

    xmax = 0;
    ymax = 0.;
    cumul_concentration_nb_dist = 0;
    max_range = 0;
    min_complement = 1.;

    for (i = 0;i < nb_dist;i++) {

      // computation of the maximum number of values, the largest support and
      // the maximum probability

      plot_nb_value = dist[i]->plot_nb_value_computation();
      if (plot_nb_value > xmax) {
        xmax = plot_nb_value;
      }
      if (dist[i]->max > ymax) {
        ymax = dist[i]->max;
      }

      if (dist[i]->variance > 0.) {
        cumul_concentration_nb_dist++;
        if (dist[i]->nb_value - dist[i]->offset > max_range) {
          max_range = dist[i]->nb_value - dist[i]->offset;
          reference_matching = i;
        }
        if (dist[i]->complement < min_complement) {
          min_complement = dist[i]->complement;
        }
      }
    }

    nb_plot_set = 1;
    if (cumul_concentration_nb_dist > 0) {
      nb_plot_set += 3;
      if (nb_dist == 1) {
        nb_plot_set--;
      }
    }

    plot_set = new MultiPlotSet(nb_plot_set);
    MultiPlotSet &plot = *plot_set;

    if (nb_dist == 1) {
      plot.title = STAT_label[STATL_DISTRIBUTION];
    }
    else {
      plot.title = STAT_label[STATL_DISTRIBUTIONS];
    }
    plot.border = "15 lw 0";

    // probability mass functions

    plot[0].xrange = Range(0 , MAX(xmax , 2) - 1);
    plot[0].yrange = Range(0. , MIN(ymax * YSCALE , 1.));

    if (MAX(xmax , 2) - 1 < TIC_THRESHOLD) {
      plot[0].xtics = 1;
    }

    plot[0].resize(nb_dist);

    for (i = 0;i < nb_dist;i++) {
      legend.str("");
      legend << STAT_label[STATL_DISTRIBUTION];
      if (nb_dist > 1) {
        legend << " " << i + 1;
      }
      dist[i]->plot_title_print(legend);
      plot[0][i].legend = legend.str();

      plot[0][i].style = "linespoints";

      dist[i]->plotable_mass_write(plot[0][i]);
    }

    if (cumul_concentration_nb_dist > 0) {

      // cumulative distribution functions

      plot[1].xrange = Range(0 , xmax - 1);
      plot[1].yrange = Range(0. , 1. - min_complement);

      if (xmax - 1 < TIC_THRESHOLD) {
        plot[1].xtics = 1;
      }

      plot[1].resize(cumul_concentration_nb_dist);

      i = 0;
      for (j = 0;j < nb_dist;j++) {
        if (dist[j]->variance > 0.) {
          legend.str("");
          legend << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_DISTRIBUTION];
          if (nb_dist > 1) {
            legend << " " << j + 1;
          }
          legend << " " << STAT_label[STATL_FUNCTION];
          plot[1][i].legend = legend.str();

          plot[1][i].style = "linespoints";

          dist[j]->plotable_cumul_write(plot[1][i]);
          i++;
        }
      }

      // matching of cumulative distribution functions taking as reference
      // the distribution with the largest support

      if (nb_dist > 1) {
        title.str("");
        title << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_DISTRIBUTION]
              << " " << STAT_label[STATL_FUNCTION] << " " << STAT_label[STATL_MATCHING];
        plot[2].title = title.str();

        plot[2].xrange = Range(0. , 1. - min_complement);
        plot[2].yrange = Range(0. , 1. - min_complement);

        plot[2].grid = true;

        plot[2].xtics = 0.1;
        plot[2].ytics = 0.1;

        plot[2].resize(cumul_concentration_nb_dist);

        i = 0;
        for (j = 0;j < nb_dist;j++) {
          if (dist[j]->variance > 0.) {
            legend.str("");
            legend << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_DISTRIBUTION]
                   << " " << j + 1 << " " << STAT_label[STATL_FUNCTION];
            plot[2][i].legend = legend.str();

            plot[2][i].style = "linespoints";

            dist[j]->plotable_cumul_matching_write(plot[2][i] , *dist[reference_matching]);
            i++;
          }
        }

        i = 3;
      }

      else {
        i = 2;
      }

      // concentration curves

      plot[i].xrange = Range(0. , 1. - min_complement);
      plot[i].yrange = Range(0. , 1. - min_complement);

      plot[i].grid = true;

      plot[i].xtics = 0.1;
      plot[i].ytics = 0.1;

      plot[i].resize(cumul_concentration_nb_dist + 1);

      j = 0;
      for (k = 0;k < nb_dist;k++) {
        if (dist[k]->variance > 0.) {
          legend.str("");
          legend << STAT_label[STATL_CONCENTRATION] << " " << STAT_label[STATL_CURVE];
          if (nb_dist > 1) {
            legend << " " << k + 1;
          }
          plot[i][j].legend = legend.str();

          plot[i][j].style = "linespoints";

          dist[k]->plotable_concentration_write(plot[i][j]);
          j++;
        }
      }

      plot[i][j].style = "lines";

      plot[i][j].add_point(0. , 0.);
      plot[i][j].add_point(1. - min_complement , 1. - min_complement);
    }

    delete [] dist;
  }

  return plot_set;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Plot of a discrete distribution.
 *
 *  \return MultiPlotSet object.
 */
/*--------------------------------------------------------------*/

MultiPlotSet* Distribution::get_plotable() const

{
  StatError error;

  return get_plotable_distributions(error , 0 , NULL);
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the survival rates from a discrete distribution and
 *         writing of the result.
 *
 *  \param[in,out] os stream.
 */
/*--------------------------------------------------------------*/

ostream& Distribution::survival_ascii_write(ostream &os) const

{
  Curves *survival_rate;


  ascii_characteristic_print(os);

  os << "\n   | " << STAT_label[STATL_DISTRIBUTION] << " | " << STAT_label[STATL_CUMULATIVE] << " "
     << STAT_label[STATL_DISTRIBUTION] << " " << STAT_label[STATL_FUNCTION] << endl;
  ascii_print(os , false , true , false);

  survival_rate = new Curves(*this);

  os << "\n   | " << STAT_label[STATL_DEATH_PROBABILITY]
     << " | " << STAT_label[STATL_SURVIVAL_PROBABILITY] << endl;
  survival_rate->ascii_print(os);

  delete survival_rate;

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation ofs survival rates from a discrete distribution and
 *         writing of the result in a file.
 *
 *  \param[in] error reference on a StatError object,
 *  \param[in] path  file path.
 *
 *  \return          error status.
 */
/*--------------------------------------------------------------*/

bool Distribution::survival_ascii_write(StatError &error , const string path) const

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
 *  \brief Computation of survival rates from a discrete distribution and
 *         writing of the result in a file at the spreadsheet format.
 *
 *  \param[in] error reference on a StatError object,
 *  \param[in] path  file path.
 *
 *  \return          error status.
 */
/*--------------------------------------------------------------*/

bool Distribution::survival_spreadsheet_write(StatError &error , const string path) const

{
  bool status;
  ofstream out_file(path.c_str());
  Curves *survival_rate;


  error.init();

  if (!out_file) {
    status = false;
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {
    status = true;

    spreadsheet_characteristic_print(out_file);

    out_file << "\n\t" << STAT_label[STATL_DISTRIBUTION] << "\t" << STAT_label[STATL_CUMULATIVE] << " "
             << STAT_label[STATL_DISTRIBUTION] << " " << STAT_label[STATL_FUNCTION] << endl;
    spreadsheet_print(out_file , true , false , false);

    survival_rate = new Curves(*this);

    out_file << "\n\t" << STAT_label[STATL_DEATH_PROBABILITY]
             << "\t" << STAT_label[STATL_SURVIVAL_PROBABILITY] << endl;
    survival_rate->spreadsheet_print(out_file);

    delete survival_rate;
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a discrete distribution and the associated survivor function
 *         at the Gnuplot format.
 *
 *  \param[in] path     file path,
 *  \param[in] survivor pointer on the survivor function.
 *
 *  \return             error status.
 */
/*--------------------------------------------------------------*/

bool Distribution::survival_plot_print(const char *path , double *survivor) const

{
  bool status = false;
  int i;
  ofstream out_file(path);


  if (out_file) {
    status = true;

    for (i = 0;i < nb_value;i++) {
      out_file << mass[i] << " " << survivor[i] << endl;
    }
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the survival rates from a discrete distribution and
 *         plot of the result using Gnuplot.
 *
 *  \param[in] error  reference on a StatError object,
 *  \param[in] prefix file prefix,
 *  \param[in] title  figure title.
 *
 *  \return           error status.
 */
/*--------------------------------------------------------------*/

bool Distribution::survival_plot_write(StatError &error , const char *prefix ,
                                       const char *title) const

{
  bool status;


  error.init();

  if (variance == 0.) {
    status = false;
    error.update(STAT_error[STATR_PLOT_NULL_VARIANCE]);
  }

  else {
    int i;
    double *survivor;
    Curves *survival_rate;
    ostringstream data_file_name[2];


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

        // probability mass function and survivor function

        if (nb_value - 1 < TIC_THRESHOLD) {
          out_file << "set xtics 0,1" << endl;
        }

        out_file << "plot [0:" << nb_value - 1 << "] [0:" << 1. - complement << "] \""
                 << label((data_file_name[0].str()).c_str()) << "\" using 1 title \""
                 << STAT_label[STATL_DISTRIBUTION] << "\" with linespoints,\\"<< endl;
        out_file << "\"" << label((data_file_name[0].str()).c_str()) << "\" using 2 title \""
                 << STAT_label[STATL_SURVIVOR] << " " << STAT_label[STATL_FUNCTION]
                 << "\" with linespoints" << endl;

        if (nb_value - 1 < TIC_THRESHOLD) {
          out_file << "set xtics autofreq" << endl;
        }

        if (i == 0) {
          out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
        }
        out_file << endl;

        // survival rates

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
 *  \brief Computation and writing of the survivor function of a discrete distribution.
 *
 *  \param[in] plot reference on a SinglePlot object.
 */
/*--------------------------------------------------------------*/

void Distribution::plotable_survivor_write(SinglePlot &plot) const

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
 *  \brief Computation of survival rates from a discrete distribution and
 *         plot of the result.
 *
 *  \param[in] error reference on a StatError object.
 *
 *  \return          MultiPlotSet object.
 */
/*--------------------------------------------------------------*/

MultiPlotSet* Distribution::survival_get_plotable(StatError &error) const

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


    plot_set = new MultiPlotSet(2);
    MultiPlotSet &plot = *plot_set;

    plot.title = "Survival analysis";
    plot.border = "15 lw 0";

    // probability mass function and survivor function

    xmax = nb_value - 1;
    if ((cumul[xmax] > 1. - DOUBLE_ERROR) &&
        (mass[xmax] > PLOT_MASS_THRESHOLD)) {
      xmax++;
    }
    plot[0].xrange = Range(0 , xmax);

    plot[0].yrange = Range(0. , 1. - complement);

    if (nb_value - 1 < TIC_THRESHOLD) {
      plot[0].xtics = 1;
    }

    plot[0].resize(2);

    legend.str("");
    legend << STAT_label[STATL_DISTRIBUTION];
    plot_title_print(legend);
    plot[0][0].legend = legend.str();

    plot[0][0].style = "linespoints";

    plotable_mass_write(plot[0][0]);

    legend.str("");
    legend << STAT_label[STATL_SURVIVOR] << " " << STAT_label[STATL_FUNCTION];
    plot[0][1].legend = legend.str();

    plot[0][1].style = "linespoints";

    plotable_survivor_write(plot[0][1]);

    // survival rates

    survival_rate = new Curves(*this);

    if (survival_rate->length - 1 < TIC_THRESHOLD) {
      plot[1].xtics = 1;
    }

    plot[1].xrange = Range(survival_rate->offset , survival_rate->length - 1);
    plot[1].yrange = Range(0. , 1.);

    plot[1].resize(2);

    plot[1][0].legend = STAT_label[STATL_DEATH_PROBABILITY];
    plot[1][0].style = "linespoints";

    plot[1][1].legend = STAT_label[STATL_SURVIVAL_PROBABILITY];
    plot[1][1].style = "linespoints";

    survival_rate->plotable_write(plot[1]);

    delete survival_rate;
  }

  return plot_set;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Display of a discrete distribution.
 *
 *  \param[in,out] os stream.
 */
/*--------------------------------------------------------------*/

ostream& Distribution::print(ostream &os) const

{
  int i;


  ascii_characteristic_print(os);

  os << "offset : " << offset;
  if (max > 0.) {
    os << "   maximum : " << max;
  }
  os << endl;

  os << "probability mass function (" << nb_value << ") : ";
  for (i = 0;i < nb_value;i++) {
    os << mass[i] << " ";
  }
  os << endl;

  os << "cumulative distribution function : ";
  for (i = 0;i < nb_value;i++) {
    os << cumul[i] << " ";
  }
  os << endl;

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Display of a discrete distribution.
 *
 *  \param[in,out] os   stream,
 *  \param[in]     dist reference on a Distribution object.
 */
/*--------------------------------------------------------------*/

ostream& operator<<(ostream &os , const Distribution &dist)

{
  streamsize nb_digits;


  nb_digits = os.precision(5);

  os << endl;
  dist.print(os);

  os.precision(nb_digits);

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the number of possible values from 0.
 */
/*--------------------------------------------------------------*/

void Distribution::nb_value_computation()

{
  double *pmass;


  pmass = mass + alloc_nb_value;
  nb_value = alloc_nb_value;

  while ((*--pmass == 0.) && (nb_value > 0)) {
    nb_value--;
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the number of values of null probability from 0.
 */
/*--------------------------------------------------------------*/

void Distribution::offset_computation()

{
  double *pmass;


  pmass = mass;
  offset = 0;

  while ((*pmass++ == 0.) && (offset < nb_value - 1)) {
    offset++;
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Extraction of the probability at the mode of a discrete distribution.
 */
/*--------------------------------------------------------------*/

void Distribution::max_computation()

{
  int i;


  max = 0.;
  for (i = offset;i < nb_value;i++) {
    if (mass[i] > max) {
      max = mass[i];
    }
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Extraction of the mode of a discrete distribution.
 *
 *  \return mode.
 */
/*--------------------------------------------------------------*/

double Distribution::mode_computation() const

{
  int i;
  double max_mass , mode;


  max_mass = 0.;
  for (i = offset;i < nb_value;i++) {
    if (mass[i] > max_mass) {
      max_mass = mass[i];
      mode = i;
    }
  }

  i = mode;
  while (mass[i + 1] == mass[i]) {
    i++;
  }
  if (i > mode) {
    mode = (i + mode) / 2.;
  }

  return mode;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the mean of a discrete distribution.
 */
/*--------------------------------------------------------------*/

void Distribution::mean_computation()

{
  if (cumul[nb_value - 1] > 0.) {
    int i;


    mean = 0.;
    for (i = offset;i < nb_value;i++) {
      mean += mass[i] * i;
    }
    mean /= cumul[nb_value - 1];
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of a quantile of a discrete distribution.
 *
 *  \param[in] icumul value of the cumulative distribution function.
 *
 *  \return           quantile.
 */
/*--------------------------------------------------------------*/

double Distribution::quantile_computation(double icumul) const

{
  int i;
  double quantile = D_DEFAULT;


  if ((complement == 0.) && (icumul <= cumul[nb_value - 1])) {
    for (i = offset;i < nb_value;i++) {
      if (cumul[i] >= icumul) {
        quantile = i;
        if (cumul[i] == icumul) {
          quantile += 0.5;
        }
        break;
      }
    }
  }

  return quantile;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the variance of a discrete distribution.
 */
/*--------------------------------------------------------------*/

void Distribution::variance_computation()

{
  if (mean != D_DEFAULT) {
    int i;
    double diff;
    long double square_sum;


    square_sum = 0.;
    for (i = offset;i < nb_value;i++) {
      diff = i - mean;
      square_sum += mass[i] * diff * diff;
    }
    variance = square_sum / cumul[nb_value - 1];

    if (variance < 0.) {
      variance = 0.;
    }
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the mean absolute deviation of a discrete distribution.
 *
 *  \param[in] location location measure (e.g. mean or median).
 *
 *  \return              mean absolute deviation.
 */
/*--------------------------------------------------------------*/

double Distribution::mean_absolute_deviation_computation(double location) const

{
  int i;
  double mean_absolute_deviation;


  mean_absolute_deviation = 0.;
  for (i = offset;i < nb_value;i++) {
    mean_absolute_deviation += mass[i] * fabs(i - location);
  }
  mean_absolute_deviation /= cumul[nb_value - 1];

  return mean_absolute_deviation;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the coefficient of skewness of a discrete distribution.
 *
 *  \return coefficient of skewness.
 */
/*--------------------------------------------------------------*/

double Distribution::skewness_computation() const

{
  int i;
  double skewness = D_INF , diff;
  long double cube_sum;


  if ((mean != D_DEFAULT) && (variance != D_DEFAULT)) {
    if (variance > 0.) {
      cube_sum = 0.;
      for (i = offset;i < nb_value;i++) {
        diff = i - mean;
        cube_sum += mass[i] * diff * diff * diff;
      }
      skewness = cube_sum / (cumul[nb_value - 1] * pow(variance , 1.5));
    }

    else {
      skewness = 0.;
    }
  }

  return skewness;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the excess kurtosis of a discrete distribution:
 *         excess kurtosis = coefficient of kurtosis - 3.
 *
 *  \return excess kurtosis.
 */
/*--------------------------------------------------------------*/

double Distribution::kurtosis_computation() const

{
  int i;
  double kurtosis = D_INF , diff;
  long double power_sum;


  if ((mean != D_DEFAULT) && (variance != D_DEFAULT)) {
    if (variance > 0.) {
      power_sum = 0.;
      for (i = offset;i < nb_value;i++) {
        diff = i - mean;
        power_sum += mass[i] * diff * diff * diff * diff;
      }
      kurtosis = power_sum / (cumul[nb_value - 1] * variance * variance) - 3.;
    }

    else {
      kurtosis = -2.;
    }
  }

  return kurtosis;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the information quantity of a discrete distribution.
 *
 *  \return information quantity.
 */
/*--------------------------------------------------------------*/

double Distribution::information_computation() const

{
  int i;
  double information = D_INF;


  if (cumul[nb_value - 1] > 0.) {
    information = 0.;
    for (i = offset;i < nb_value;i++) {
      if (mass[i] > 0.) {
        information += mass[i] * log(mass[i]);
      }
    }

    if (complement > 0.) {
      information -= cumul[nb_value - 1] * log(1. - complement);
    }

    information /= cumul[nb_value - 1];
  }

  return information;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the sum of squared first-order differences.
 *
 *  \return sum of squared first-order differences.
 */
/*--------------------------------------------------------------*/

double Distribution::first_difference_norm_computation() const

{
  int i;
  double first_difference_norm , buff;


  first_difference_norm = mass[offset] * mass[offset];
  for (i = offset;i < nb_value - 1;i++) {
    buff = mass[i + 1] - mass[i];
    first_difference_norm += buff * buff;
  }
  first_difference_norm += mass[nb_value - 1] * mass[nb_value - 1];
  first_difference_norm /= cumul[nb_value - 1];

  return first_difference_norm;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the sum of squared second-order differences.
 *
 *  \return sum of squared second-order differences.
 */
/*--------------------------------------------------------------*/

double Distribution::second_difference_norm_computation() const

{
  int i;
  double second_difference_norm , buff;


  second_difference_norm = mass[offset] * mass[offset];
  buff = -2 * mass[offset] + mass[offset + 1];
  second_difference_norm += buff * buff;

  for (i = offset + 1;i < nb_value - 1;i++) {
    buff = mass[i - 1] - 2 * mass[i] + mass[i + 1];
    second_difference_norm += buff * buff;
  }

  buff = mass[nb_value - 2] - 2 * mass[nb_value - 1];
  second_difference_norm += buff * buff;
  second_difference_norm += mass[nb_value - 1] * mass[nb_value - 1];

  second_difference_norm /= cumul[nb_value - 1];

  return second_difference_norm;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the cumulative distribution function of a discrete distribution.
 *
 *  \param[in] nb_value number of values,
 *  \param[in] pmass    pointer on the probability mass function,
 *  \param[in] pcumul   pointer on the cumulative distribution function.
 */
/*--------------------------------------------------------------*/

void cumul_computation(int nb_value , const double *pmass , double *pcumul)

{
  int i;


  *pcumul = *pmass;
  for (i = 1;i < nb_value;i++) {
    pcumul++;
    *pcumul = *(pcumul - 1) + *++pmass;
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the cumulative distribution function of a discrete distribution.
 */
/*--------------------------------------------------------------*/

void Distribution::cumul_computation()

{
  int i;


  for (i = 0;i < offset;i++) {
    cumul[i] = 0.;
  }
  stat_tool::cumul_computation(nb_value - offset , mass + offset , cumul + offset);
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the survivor function of a discrete distribution.
 *
 *  \return survivor function.
 */
/*--------------------------------------------------------------*/

double* Distribution::survivor_function_computation() const

{
  int i;
  double *survivor_function;


  survivor_function = new double[nb_value];

  for (i = 0;i < offset;i++) {
    survivor_function[i] = 1. - complement;
  }
  for (i = offset;i < nb_value;i++) {
    survivor_function[i] = 1. - complement - cumul[i];
  }

  return survivor_function;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the concentration function of a discrete distribution.
 *
 *  \return concentration function.
 */
/*--------------------------------------------------------------*/

double* Distribution::concentration_function_computation() const

{
  int i;
  double *concentration_function;


  if ((mean > 0.) && (variance > 0.)) {
    concentration_function = new double[nb_value];

    for (i = 0;i < offset;i++) {
      concentration_function[i] = 0.;
    }
    concentration_function[offset] = mass[offset] * offset / mean;
    for (i = offset + 1;i < nb_value;i++) {
      concentration_function[i] = concentration_function[i - 1] + mass[i] * i / mean;
    }
  }

  else {
    concentration_function = NULL;
  }

  return concentration_function;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the coefficient of concentration of a discrete distribution.
 *
 *  \return coefficient of concentration.
 */
/*--------------------------------------------------------------*/

double Distribution::concentration_computation() const

{
  int i;
  double concentration = D_DEFAULT , *concentration_function;


  if ((mean > 0.) && (variance > 0.)) {
    concentration_function = concentration_function_computation();

    concentration = mass[offset] * concentration_function[offset];
    for (i = offset + 1;i < nb_value;i++) {
      concentration += mass[i] * (concentration_function[i - 1] + concentration_function[i]);
    }

    concentration = 1. - concentration / (cumul[nb_value - 1] * concentration_function[nb_value - 1]);

    delete [] concentration_function;

#   ifdef DEBUG
    int previous_value;
    double concentration2 = 0.;

    previous_value = offset;
    for (i = offset + 1;i < nb_value;i++) {
      if (mass[i] > 0.) {
        concentration2 += cumul[i - 1] * (1. - cumul[i - 1]) * (i - previous_value);
//        concentration2 += cumul[i - 1] * (cumul[nb_value - 1] - cumul[i - 1]) * (i - previous_value);
        previous_value = i;
      }
    }

    concentration2 /= mean;

    cout << "\n" << STAT_label[STATL_CONCENTRATION_COEFF] << ": " << concentration << " | " << concentration2 << endl;
#   endif

  }

  return concentration;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the distance between 2 discrete distributions (1 - overlap).
 *
 *  \return overlap distance.
 */
/*--------------------------------------------------------------*/

double Distribution::overlap_distance_computation(const Distribution &dist) const

{
  int i;
  double overlap;


  overlap = 0.;
  for (i = MAX(offset , dist.offset);i < MIN(nb_value , dist.nb_value);i++) {
    overlap += MIN(mass[i] , dist.mass[i]);
  }

  return (1. - overlap);
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the log-probability mass function.
 *
 *  \param[in] nb_value number of values,
 *  \param[in] pmass    pointer on the probability mass function,
 *  \param[in] plog     pointer on the log-probability mass function.
 */
/*--------------------------------------------------------------*/

void log_computation(int nb_value , const double *pmass , double *plog)

{
  int i;


  for (i = 0;i < nb_value;i++) {
    if (*pmass > 0.) {
      *plog = log(*pmass);
    }
    else {
      *plog = D_INF;
    }
    plog++;
    pmass++;
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Ccomputation of the log-probability mass function.
 */
/*--------------------------------------------------------------*/

void Distribution::log_computation()

{
  stat_tool::log_computation(nb_value , mass , cumul);
}


};  // namespace stat_tool
