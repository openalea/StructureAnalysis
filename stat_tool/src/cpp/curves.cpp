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



#include <iomanip>

#include "curves.h"
#include "stat_label.h"

using namespace std;


namespace stat_tool {



/*--------------------------------------------------------------*/
/**
 *  \brief Default constructor of the Curves class.
 */
/*--------------------------------------------------------------*/

Curves::Curves()

{
  nb_curve = 0;
  length = 0;
  offset = 0;
  index_parameter = NULL;
  frequency = NULL;
  point = NULL;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Constructor of the Curves class.
 *
 *  \param[in] inb_curve            number of curves,
 *  \param[in] ilength              curve length,
 *  \param[in] frequency_flag       flag on the frequencies,
 *  \param[in] index_parameter_flag flag on the index parameters,
 *  \param[in] init_flag            flag initialization.
 */
/*--------------------------------------------------------------*/

Curves::Curves(int inb_curve , int ilength , bool frequency_flag ,
               bool index_parameter_flag , bool init_flag)

{
  int i , j;


  nb_curve = inb_curve;
  length = ilength;
  offset = 0;

  if (index_parameter_flag) {
    index_parameter = new int[length];
  }
  else {
    index_parameter = NULL;
  }

  if (frequency_flag) {
    frequency = new int[length];

    if (init_flag) {
      for (i = 0;i < length;i++) {
        frequency[i] = 0;
      }
    }
  }

  else {
    frequency = NULL;
  }

  point = new double*[nb_curve];
  for (i = 0;i < nb_curve;i++) {
    point[i] = new double[length];

    if (init_flag) {
      for (j = 0;j < length;j++) {
        point[i][j] = 0.;
      }
    }
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Construction of a Curves object from a Distribution object.
 *
 *  \param[in] dist reference on a Distribution object.
 */
/*--------------------------------------------------------------*/

Curves::Curves(const Distribution &dist)

{
  int i;
  double norm , *pcumul;


  nb_curve = 2;

  pcumul = dist.cumul + dist.nb_value - 1;
  length = dist.nb_value;
  while ((length > 1) && (*pcumul == *(pcumul - 1))) {
    pcumul--;
    length--;
  }

  offset = 0;
  index_parameter = NULL;
  frequency = NULL;

  point = new double*[2];
  for (i = 0;i < 2;i++) {
    point[i] = new double[length];
  }

  norm = 1. - dist.complement;
  for (i = 0;i < length;i++) {
    if ((dist.mass[i] <= norm) && ((1. - dist.complement - dist.cumul[i]) <= norm)) {
      point[0][i] = dist.mass[i] / norm;
      point[1][i] = (1. - dist.complement - dist.cumul[i]) / norm;
      norm = 1. - dist.complement - dist.cumul[i];
    }

    else {
      break;
    }
  }

  length = i;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Construction of a Curves object from a FrequencyDistribution object.
 *
 *  \param[in] histo reference on a FrequencyDistribution object.
 */
/*--------------------------------------------------------------*/

Curves::Curves(const FrequencyDistribution &histo)

{
  int i;
  int norm;


  nb_curve = 2;
  length = histo.nb_value;
  offset = 0;

  index_parameter = NULL;
  frequency = new int[length];

  point = new double*[2];
  for (i = 0;i < 2;i++) {
    point[i] = new double[length];
  }

  norm = histo.nb_element;
  for (i = 0;i < length;i++) {
    point[0][i] = (double)histo.frequency[i] / (double)norm;
    point[1][i] = (double)(norm - histo.frequency[i]) / (double)norm;
    frequency[i] = norm;
    norm -= histo.frequency[i];
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Copy of a Curves object.
 *
 *  \param[in] curves reference on a Curves object.
 */
/*--------------------------------------------------------------*/

void Curves::copy(const Curves &curves)

{
  int i , j;


  nb_curve = curves.nb_curve;
  length = curves.length;
  offset = curves.offset;

  if (curves.index_parameter) {
    index_parameter = new int[length];

    for (i = 0;i < length;i++) {
      index_parameter[i] = curves.index_parameter[i];
    }
  }

  else {
    index_parameter = NULL;
  }

  if (curves.frequency) {
    frequency = new int[length];

    for (i = 0;i < length;i++) {
      frequency[i] = curves.frequency[i];
    }
  }

  else {
    frequency = NULL;
  }

  point = new double*[nb_curve];
  for (i = 0;i < nb_curve;i++) {
    point[i] = new double[length];

    for (j = 0;j < length;j++) {
      point[i][j] = curves.point[i][j];
    }
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Copy of a Curves object with curve smoothing.
 *
 *  \param[in] curves        reference on a Curves object,
 *  \param[in] max_frequency threshold on the frequencies for the smoothing.
 */
/*--------------------------------------------------------------*/

void Curves::smooth(const Curves &curves , int max_frequency)

{
  int i , j , k;
  int range , buff , min , max , total;
  double sum;


  nb_curve = curves.nb_curve;
  length = curves.length;
  offset = curves.offset;

  if (curves.index_parameter) {
    index_parameter = new int[length];

    for (i = 0;i < length;i++) {
      index_parameter[i] = curves.index_parameter[i];
    }
  }

  else {
    index_parameter = NULL;
  }

  frequency = new int[length];

  point = new double*[nb_curve];
  for (i = 0;i < nb_curve;i++) {
    point[i] = new double[length];
  }

  for (i = 0;i < offset;i++) {
    for (j = 0;j < nb_curve;j++) {
      point[j][i] = 0.;
    }
    frequency[i] = 0;
  }

  for (i = offset;i < length;i++) {
    frequency[i] = curves.frequency[i];

    if (frequency[i] < max_frequency) {
      range = MAX_RANGE;
      if ((frequency[i] > 0) && (max_frequency / frequency[i] < range)) {
        range = max_frequency / frequency[i];
      }

      buff = MIN(range , i);
      buff = MIN(buff , length - 1 - i);
      min = i - buff;
      max = i + buff;

      total = 0;
      for (j = min;j <= max;j++) {
        total += curves.frequency[j];
      }

      if (total > 0) {
        for (j = 0;j < nb_curve;j++) {
          sum = 0.;
          for (k = min;k <= max;k++) {
            sum += curves.point[j][k] * curves.frequency[k];
          }
          point[j][i] = sum / total;
        }
      }

      else {
        for (j = 0;j < nb_curve;j++) {
          point[j][i] = curves.point[j][i];
        }
      }
    }

    else {
      for (j = 0;j < nb_curve;j++) {
        point[j][i] = curves.point[j][i];
      }
    }
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Constructor by copy of the Curves class.
 *
 *  \param[in] curves        reference on a Curves object,
 *  \param[in] transform     type of transform (CURVE_COPY/SMOOTHING),
 *  \param[in] max_frequency threshold on the frequencies for the smoothing.
 */
/*--------------------------------------------------------------*/

Curves::Curves(const Curves &curves , curve_transformation transform , int max_frequency)

{
  if (!(curves.frequency)) {
    transform = CURVE_COPY;
  }

  switch (transform) {
  case SMOOTHING :
    smooth(curves , max_frequency);
    break;
  default :
    copy(curves);
    break;
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Destruction of data members of a Curves object.
 */
/*--------------------------------------------------------------*/

void Curves::remove()

{
  delete [] index_parameter;
  delete [] frequency;

  if (point) {
    int i;

    for (i = 0;i < nb_curve;i++) {
      delete [] point[i];
    }
    delete [] point;
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Destructor of the Curves class.
 */
/*--------------------------------------------------------------*/

Curves::~Curves()

{
  remove();
}


/*--------------------------------------------------------------*/
/**
 *  \brief Assignment operator of the Curves class.
 *
 *  \param[in] curves reference on a Curves object.
 *
 *  \return           Curves object.
 */
/*--------------------------------------------------------------*/

Curves& Curves::operator=(const Curves &curves)

{
  if (&curves != this) {
    remove();
    copy(curves);
  }

  return *this;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of 2 families of curves of the same length.
 *
 *  \param[in,out] os        stream,
 *  \param[in]     file_flag flag file,
 *  \param[in]     curves    pointer on a Curves object.
 */
/*--------------------------------------------------------------*/

ostream& Curves::ascii_print(ostream &os , bool file_flag , const Curves *curves) const

{
  int i , j , k;
  int nb_column = nb_curve + 1 , *width;
  ios_base::fmtflags format_flags;


  format_flags = os.setf(ios::right , ios::adjustfield);

  // computation of the column width

  if (frequency) {
    nb_column++;
  }
  if (curves) {
    nb_column += nb_curve;
    if (curves->frequency) {
      nb_column++;
    }
  }

  width = new int[nb_column];

  if (index_parameter) {
    width[0] = column_width(index_parameter[length - 1]);
  }
  else {
    width[0] = column_width(length - 1);
  }

  i = 1;
  for (j = 0;j < nb_curve;j++) {
    if ((curves) && (j < curves->nb_curve)) {
      width[i++] = column_width(length - offset , curves->point[j] + offset) + ASCII_SPACE;
    }
    width[i++] = column_width(length - offset , point[j] + offset) + ASCII_SPACE;
  }

  if ((curves) && (curves->frequency)) {
    width[i++] = column_width(curves->frequency[offset]) + ASCII_SPACE;
  }
  if (frequency) {
    width[i] = column_width(frequency[offset]) + ASCII_SPACE;
  }

  // writing of curves

  for (i = offset;i < length;i++) {
    if (file_flag) {
      os << "# ";
    }
    if (index_parameter) {
      os << setw(width[0]) << index_parameter[i];
    }
    else {
      os << setw(width[0]) << i;
    }

    j = 1;
    for (k = 0;k < nb_curve;k++) {
      if ((curves) && (k < curves->nb_curve)) {
        os << setw(width[j++]) << curves->point[k][i];
      }
      os << setw(width[j++]) << point[k][i];
    }

    if ((curves) && (curves->frequency)) {
      os << setw(width[j++]) << curves->frequency[i];
    }
    if (frequency) {
      os << setw(width[j]) << frequency[i];
    }
    os << endl;
  }

  delete [] width;

  os.setf(format_flags , ios::adjustfield);

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of 2 families of curves of the same length at the spreadsheet format.
 *
 *  \param[in,out] os     stream,
 *  \param[in]     curves pointer on a Curves object.
 */
/*--------------------------------------------------------------*/

ostream& Curves::spreadsheet_print(ostream &os , const Curves *curves) const

{
  int i , j;


  for (i = offset;i < length;i++) {
    if (index_parameter) {
      os << index_parameter[i];
    }
    else {
      os << i;
    }

    for (j = 0;j < nb_curve;j++) {
      if ((curves) && (j < curves->nb_curve)) {
        os << "\t" << curves->point[j][i];
      }
      os << "\t" << point[j][i];
    }

    if ((curves) && (curves->frequency)) {
      os << "\t" << curves->frequency[i];
    }
    if (frequency) {
      os << "\t" << frequency[i];
    }
    os << endl;
  }

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the curve length to be plotted (Gnuplot output).
 *
 *  \return curve length.
 */
/*--------------------------------------------------------------*/

int Curves::plot_length_computation() const

{
  int plot_length = length , min_frequency , *pfrequency;


  if (frequency) {
    pfrequency = frequency + plot_length;
    min_frequency = MIN(max_frequency_computation() / 2 , PLOT_MIN_FREQUENCY);

    while (*--pfrequency < min_frequency) {
      plot_length--;
    }
  }

  return plot_length;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of 2 families of curves at the Gnuplot format.
 *
 *  \param[in] path     file path,
 *  \param[in] ilength  curve length,
 *  \param[in] curves_0 pointer on a Curves object,
 *  \param[in] curves_1 pointer on a Curves object.
 *
 *  \return             error status.
 */
/*--------------------------------------------------------------*/

bool Curves::plot_print(const char *path , int ilength ,
                        const Curves *curves_0 , const Curves *curves_1) const

{
  bool status = false;
  int i , j;
  ofstream out_file(path);


  if (out_file) {
    status = true;

    if (ilength == I_DEFAULT) {
      ilength = length;
    }

    for (i = 0;i < ilength;i++) {
      if (index_parameter) {
        out_file << index_parameter[i] << " ";
      }
      for (j = 0;j < nb_curve;j++) {
        out_file << point[j][i] << " ";
      }

      if (curves_0) {
        for (j = 0;j < curves_0->nb_curve;j++) {
          out_file << curves_0->point[j][i] << " ";
        }
      }

      if (curves_1) {
        for (j = 0;j < curves_1->nb_curve;j++) {
          out_file << curves_1->point[j][i] << " ";
        }
      }

      out_file << endl;
    }
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a curve.
 *
 *  \param[in] index curve index,
 *  \param[in] plot  reference on a SinglePlot object.
 */
/*--------------------------------------------------------------*/

void Curves::plotable_write(int index , SinglePlot &plot) const

{
  int i;


  if (index_parameter) {
    for (i = offset;i < length;i++) {
      plot.add_point(index_parameter[i] , point[index][i]);
    }
  }

  else {
    for (i = offset;i < length;i++) {
      plot.add_point(i , point[index][i]);
    }
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of the curve family.
 *
 *  \param[in] plot reference on a MultiPlot object.
 */
/*--------------------------------------------------------------*/

void Curves::plotable_write(MultiPlot &plot) const

{
  int i , j;


  if (index_parameter) {
    for (i = offset;i < length;i++) {
      for (j = 0;j < nb_curve;j++) {
        plot[j].add_point(index_parameter[i] , point[j][i]);
      }
    }
  }

  else {
    for (i = offset;i < length;i++) {
      for (j = 0;j < nb_curve;j++) {
        plot[j].add_point(i , point[j][i]);
      }
    }
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of frequencies.
 *
 *  \param[in] plot reference on a SinglePlot object.
 */
/*--------------------------------------------------------------*/

void Curves::plotable_frequency_write(SinglePlot &plot) const

{
  int i;


  if (index_parameter) {
    for (i = offset;i < length;i++) {
      plot.add_point(index_parameter[i] , frequency[i]);
    }
  }

  else {
    for (i = offset;i < length;i++) {
      plot.add_point(i , frequency[i]);
    }
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Display of a Curves object.
 *
 *  \param[in,out] os     stream,
 *  \param[in]     curves reference on a Curves object.
 */
/*--------------------------------------------------------------*/

ostream& operator<<(ostream &os , const Curves &curves)

{
  int i , j;
  streamsize nb_digits;


  nb_digits = os.precision(4);

  os << endl;
  if (curves.nb_curve == 1) {
    os << curves.nb_curve << " curve";
  }
  else {
    os << curves.nb_curve << " curves";
  }
  os << "   length : " << curves.length << endl;

  for (i = curves.offset;i < curves.length;i++) {
    os << "(";
    for (j = 0;j < curves.nb_curve;j++) {
      os << curves.point[j][i] << " ";
    }

    if (curves.frequency) {
      os << "| " << curves.frequency[i];
    }
    os << ") ";
  }
  os << endl;

  os.precision(nb_digits);

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the maximum frequency of a Curves object.
 *
 *  \return maximum frequency.
 */
/*--------------------------------------------------------------*/

int Curves::max_frequency_computation() const

{
  int max_frequency = I_DEFAULT;


  if (frequency) {
    int i;


    max_frequency = 0;
    for (i = offset;i < length;i++) {
      if (frequency[i] > max_frequency) {
        max_frequency = frequency[i];
      }
    }
  }

  return max_frequency;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the cumulative frequency of a Curves object.
 *
 *  \return cumulative frequency.
 */
/*--------------------------------------------------------------*/

int Curves::nb_element_computation() const

{
  int nb_element = I_DEFAULT;


  if (frequency) {
    int i;


    nb_element = 0;
    for (i = offset;i < length;i++) {
      nb_element += frequency[i];
    }
  }

  return nb_element;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the curve mean.
 *
 *  \param[in] index curve index.
 *
 *  \return          mean.
 */
/*--------------------------------------------------------------*/

double Curves::mean_computation(int index) const

{
  int i;
  int nb_element;
  double mean;


  nb_element = 0;
  mean = 0.;
  for (i = offset;i < length;i++) {
    if (frequency[i] > 0) {
      nb_element += frequency[i];
      mean += frequency[i] * point[index][i];
    }
  }
  mean /= nb_element;

  return mean;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of the total variation for a curve.
 *
 *  \param[in] index curve index,
 *  \param[in] mean  mean.
 *
 *  \return          total variation.
 */
/*--------------------------------------------------------------*/

double Curves::total_square_sum_computation(int index , double mean) const

{
  int i;
  double total_square_sum , diff;


  total_square_sum = 0.;
  for (i = offset;i < length;i++) {
    if (frequency[i] > 0) {
      diff = point[index][i] - mean;
      total_square_sum += frequency[i] * diff * diff;
    }
  }

  return total_square_sum;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a curve and the corresponding standardized residuals at the Gnuplot format.
 *
 *  \param[in] path              file path,
 *  \param[in] standard_residual pointer on the standardized residuals.
 *
 *  \return                      error status.
 */
/*--------------------------------------------------------------*/

bool Curves::plot_print_standard_residual(const char *path , double *standard_residual) const

{
  bool status = false;
  int i;
  ofstream out_file(path);


  if (out_file) {
    status = true;

    // writing of observed responses and standardized residuals

    for (i = 0;i < length;i++) {
      if (frequency[i] > 0) {
        out_file << i << " " << point[0][i];
        if (standard_residual) {
          out_file << " " << standard_residual[i];
        }
        out_file << " " << frequency[i] << endl;
      }
    }
  }

  return status;
}


};  // namespace stat_tool
