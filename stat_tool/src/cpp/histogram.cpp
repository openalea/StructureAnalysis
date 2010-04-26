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
 *       $Id: distribution.cpp 8644 2010-04-14 13:09:28Z guedon $
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


extern int column_width(int value);



/*--------------------------------------------------------------*
 *
 *  Construction d'un objet Histogram.
 *
 *  arguments : nombre de categories, flag initialisation.
 *
 *--------------------------------------------------------------*/

Histogram::Histogram(int inb_category , bool init_flag)

{
  nb_individual = I_DEFAULT;
  nb_category = inb_category;
  step = D_DEFAULT;
  max = 0;

  if (nb_category == 0) {
    frequency = NULL;
  }

  else {
    frequency = new int[nb_category];

    if (init_flag) {
      register int i;

      for (i = 0;i < nb_category;i++) {
        frequency[i] = 0;
      }
    }
  }

  min_value = D_DEFAULT;
  max_value = D_DEFAULT;
}


/*--------------------------------------------------------------*
 *
 *  Copie d'un objet Histogram.
 *
 *  argument : reference sur un objet Histogram.
 *
 *--------------------------------------------------------------*/

void Histogram::copy(const Histogram &histo)

{
  register int i;


  nb_individual = histo.nb_individual;
  nb_category = histo.nb_category;
  step = histo.step;
  max = histo.max;

  frequency = new int[nb_category];

  for (i = 0;i < nb_category;i++) {
    frequency[i] = histo.frequency[i];
  }

  min_value = histo.min_value;
  max_value = histo.max_value;
}


/*--------------------------------------------------------------*
 *
 *  Constructeur par copie de la classe Histogram.
 *
 *  argument : reference sur un objet Histogram.
 *
 *--------------------------------------------------------------*/

Histogram::Histogram(const Histogram &histo)

{
  copy(histo);
}


/*--------------------------------------------------------------*
 *
 *  Destructeur de la classe Histogram.
 *
 *--------------------------------------------------------------*/

Histogram::~Histogram()

{
  delete [] frequency;
}


/*--------------------------------------------------------------*
 *
 *  Operateur d'assignement de la classe Histogram.
 *
 *  argument : reference sur un objet Histogram.
 *
 *--------------------------------------------------------------*/

Histogram& Histogram::operator=(const Histogram &histo)

{
  if (&histo != this) {
    delete [] frequency;

    copy(histo);
  }

  return *this;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la largeur d'une colonne de reels.
 *
 *  arguments : valeur minimum, valeur maximum.
 *
 *--------------------------------------------------------------*/

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


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un histogramme.
 *
 *  arguments : stream, flag commentaire.
 *
 *--------------------------------------------------------------*/

ostream& Histogram::ascii_print(ostream &os , bool comment_flag) const

{
  register int i;
  int width[2];
  long old_adjust;
  double value;


  old_adjust = os.setf(ios::right , ios::adjustfield);

  // calcul des largeurs des colonnes

  width[0] = column_width(min_value + step / 2 , min_value + (nb_category - 0.5) * step);
  width[1] = column_width(max) + ASCII_SPACE;

  value = min_value + step / 2;
  for (i = 0;i < nb_category;i++) {
//    if (frequency[i] > 0) {
      if (comment_flag) {
        os << "# ";
      }
      os << setw(width[0]) << value;
      os << setw(width[1]) << frequency[i];
      os << endl;
//    }
    value += step;
  }

  os.setf((FMTFLAGS)old_adjust , ios::adjustfield);

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un histogramme au format tableur.
 *
 *  argument : stream.
 *
 *--------------------------------------------------------------*/

ostream& Histogram::spreadsheet_print(ostream &os) const

{
  register int i;
  double value;


  value = min_value + step / 2;
  for (i = 0;i < nb_category;i++) {
//    if (frequency[i] > 0) {
      os << value << "\t" << frequency[i] << endl;
//    }
    value += step;
  }

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un histogramme au format Gnuplot.
 *
 *  argument : path.
 *
 *--------------------------------------------------------------*/

bool Histogram::plot_print(const char *path) const

{
  bool status = false;
  register int i;
  double value;
  ofstream out_file(path);


  if (out_file) {
    status = true;

    value = min_value + step / 2;
    for (i = 0;i < nb_category;i++) {
      out_file << value << " " << frequency[i] << endl;
      value += step;
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

void Histogram::plotable_write(SinglePlot &plot) const

{
  register int i;
  double value;


  value = min_value + step / 2;
  for (i = 0;i < nb_category;i++) {
    plot.add_point(value , frequency[i]);
    value += step;
  }
}


/*--------------------------------------------------------------*
 *
 *  Recherche de la frequence maximum d'un histogramme.
 *
 *--------------------------------------------------------------*/

void Histogram::max_computation()

{
  register int i;


  max = 0;
  for (i = 0;i < nb_category;i++) {
    if (frequency[i] > max) {
      max = frequency[i];
    }
  }
}
