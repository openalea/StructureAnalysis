/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       V-Plants: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2016 CIRAD/INRA/Inria Virtual Plants
 *
 *       File author(s): Yann Guedon (yann.guedon@cirad.fr)
 *
 *       $Source$
 *       $Id: histogram.cpp 8644 2010-04-14 13:09:28Z guedon $
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

#include "tool/config.h"

#include "stat_tools.h"

using namespace std;


namespace stat_tool {



/*--------------------------------------------------------------*/
/**
 *  \brief Construction of a Histogram object.
 *
 *  \param[in] inb_bin   number of bins,
 *  \param[in] init_flag flag initialization.
 */
/*--------------------------------------------------------------*/

Histogram::Histogram(int inb_bin , bool init_flag)

{
  nb_element = I_DEFAULT;
  nb_bin = inb_bin;
  bin_width = D_DEFAULT;
  max = 0;

  if (nb_bin == 0) {
    frequency = NULL;
  }

  else {
    frequency = new int[nb_bin];

    if (init_flag) {
      register int i;

      for (i = 0;i < nb_bin;i++) {
        frequency[i] = 0;
      }
    }
  }

  type = I_DEFAULT;
  min_value = D_DEFAULT;
  max_value = D_DEFAULT;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Construction of an Histogram object from a FrequencyDistribution object.
 *
 *  \param[in] histo reference on a FrequencyDistribution object.
 */
/*--------------------------------------------------------------*/

Histogram::Histogram(const FrequencyDistribution &histo)

{
  register int i , j;


  nb_element = histo.nb_element;
  max = histo.max;

  bin_width = histo.min_interval_computation();

  nb_bin = (histo.nb_value - 1 - histo.offset) / bin_width + 1;
  frequency = new int[nb_bin];

  i = histo.offset;
  for (j = 0;j < nb_bin;j++) {
    frequency[j] = histo.frequency[i];
    i += bin_width;
  }

  type = INT_VALUE;

  min_value = histo.offset;
  max_value = histo.nb_value - 1;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Copy of an Histogram object.
 *
 *  \param[in] reference on a Histogram object.
 */
/*--------------------------------------------------------------*/

void Histogram::copy(const Histogram &histo)

{
  register int i;


  nb_element = histo.nb_element;
  nb_bin = histo.nb_bin;
  bin_width = histo.bin_width;
  max = histo.max;

  frequency = new int[nb_bin];

  for (i = 0;i < nb_bin;i++) {
    frequency[i] = histo.frequency[i];
  }

  type = histo.type;
  min_value = histo.min_value;
  max_value = histo.max_value;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Destructor of the Histogram class.
 */
/*--------------------------------------------------------------*/

Histogram::~Histogram()

{
  delete [] frequency;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Assignment operator of the Histogram class.
 *
 *  \param[in]  histo reference on a Histogram object,
 *
 *  \param[out] this  Histogram object.
 */
/*--------------------------------------------------------------*/

Histogram& Histogram::operator=(const Histogram &histo)

{
  if (&histo != this) {
    delete [] frequency;

    copy(histo);
  }

  return *this;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of a column width for real values.
 *
 *  \param[in]  min_value minimum value,
 *  \param[in]  max_value maximum value,
 *
 *  \param[out] max_width column width.
 */
/*--------------------------------------------------------------*/

int column_width(double min_value , double max_value)

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
 *  \brief Writing of a histogram.
 *
 *  \param[in,out] os           stream,
 *  \param[in]     comment_flag flag comment.
 */
/*--------------------------------------------------------------*/

ostream& Histogram::ascii_print(ostream &os , bool comment_flag) const

{
  register int i;
  int width[2];
  long old_adjust;
  double first_value , value;


  old_adjust = os.setf(ios::right , ios::adjustfield);

  // computation of column widths

  switch (type) {
  case INT_VALUE :
    first_value = min_value;
    break;
  case REAL_VALUE :
    first_value = min_value + bin_width / 2;
    break;
  }

  width[0] = column_width(first_value , first_value + (nb_bin - 1) * bin_width);
  width[1] = column_width(max) + ASCII_SPACE;

  value = first_value;
  for (i = 0;i < nb_bin;i++) {
//    if (frequency[i] > 0) {
      if (comment_flag) {
        os << "# ";
      }
      os << setw(width[0]) << value;
      os << setw(width[1]) << frequency[i];
      os << endl;
//    }
    value += bin_width;
  }

  os.setf((FMTFLAGS)old_adjust , ios::adjustfield);

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a histogram at the spreadsheet format.
 *
 *  \param[in,out] os stream.
 */
/*--------------------------------------------------------------*/

ostream& Histogram::spreadsheet_print(ostream &os) const

{
  register int i;
  double value;


  switch (type) {
  case INT_VALUE :
    value = min_value;
    break;
  case REAL_VALUE :
    value = min_value + bin_width / 2;
    break;
  }

  for (i = 0;i < nb_bin;i++) {
//    if (frequency[i] > 0) {
      os << value << "\t" << frequency[i] << endl;
//    }
    value += bin_width;
  }

  return os;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a histogram at the Gnuplot format.
 *
 *  \param[in]  path   file path,
 *
 *  \param[out] status error status.
 */
/*--------------------------------------------------------------*/

bool Histogram::plot_print(const char *path) const

{
  bool status = false;
  register int i;
  double value;
  ofstream out_file(path);


  if (out_file) {
    status = true;

    switch (type) {
    case INT_VALUE :
      value = min_value;
      break;
    case REAL_VALUE :
      value = min_value + bin_width / 2;
      break;
    }

    out_file << value - bin_width << " " << 0 << endl;
    for (i = 0;i < nb_bin;i++) {
      out_file << value << " " << frequency[i] << endl;
      value += bin_width;
    }
    out_file << value << " " << 0 << endl;
  }

  return status;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a histogram at the plotable format.
 *
 *  \param[in] plot reference on a SinglePlot object.
 */
/*--------------------------------------------------------------*/

void Histogram::plotable_write(SinglePlot &plot) const

{
  register int i;
  double value;


  switch (type) {
  case INT_VALUE :
    value = min_value;
    break;
  case REAL_VALUE :
    value = min_value + bin_width / 2;
    break;
  }

  plot.add_point(value - bin_width , 0);
  for (i = 0;i < nb_bin;i++) {
    plot.add_point(value , frequency[i]);
    value += bin_width;
  }
  plot.add_point(value , 0);
}


/*--------------------------------------------------------------*/
/**
 *  \brief Determination of the maximum bin frequency for a histogram.
 */
/*--------------------------------------------------------------*/

void Histogram::max_computation()

{
  register int i;


  max = 0;
  for (i = 0;i < nb_bin;i++) {
    if (frequency[i] > max) {
      max = frequency[i];
    }
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Computation of a cumulative distribution function from a histogram.
 *
 *  \param[out] cumul cumulative distribution function.
 */
/*--------------------------------------------------------------*/

double* Histogram::cumul_computation() const

{
  register int i;
  double *cumul;


  cumul = new double[nb_bin];

  cumul[0] = (double)frequency[0] / (double)nb_element;
  for (i = 1;i < nb_bin;i++) {
    cumul[i] = cumul[i - 1] + (double)frequency[i] / (double)nb_element;
  }

  return cumul;
}


};  // namespace stat_tool
