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
 *       $Id: curves.cpp 17987 2015-04-23 06:43:48Z guedon $
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

#include "tool/config.h"

#include "stat_tools.h"
#include "curves.h"
#include "stat_label.h"

using namespace std;


namespace stat_tool {



/*--------------------------------------------------------------*
 *
 *  Constructeur par defaut de la classe Curves.
 *
 *--------------------------------------------------------------*/

Curves::Curves()

{
  nb_curve = 0;
  length = 0;
  offset = 0;
  index_parameter = NULL;
  frequency = NULL;
  point = NULL;
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Curves.
 *
 *  arguments : nombre de courbes, longueur des courbes, flag sur les frequences,
 *              flag sur les parametres d'index, flag initialisation.
 *
 *--------------------------------------------------------------*/

Curves::Curves(int inb_curve , int ilength , bool frequency_flag ,
               bool index_parameter_flag , bool init_flag)

{
  register int i , j;


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


/*--------------------------------------------------------------*
 *
 *  Construction d'un objet Curves a partir d'un objet Distribution.
 *
 *  argument : reference sur un objet Distribution.
 *
 *--------------------------------------------------------------*/

Curves::Curves(const Distribution &dist)

{
  register int i;
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


/*--------------------------------------------------------------*
 *
 *  Construction d'un objet Curves a partir d'un objet FrequencyDistribution.
 *
 *  argument : reference sur un objet FrequencyDistribution.
 *
 *--------------------------------------------------------------*/

Curves::Curves(const FrequencyDistribution &histo)

{
  register int i;
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


/*--------------------------------------------------------------*
 *
 *  Copie d'un objet Curves.
 *
 *  argument : reference sur un objet Curves.
 *
 *--------------------------------------------------------------*/

void Curves::copy(const Curves &curves)

{
  register int i , j;


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


/*--------------------------------------------------------------*
 *
 *  Copie d'un objet Curves avec lissage.
 *
 *  arguments : reference sur un objet Curves,
 *              seuil sur les frequences pour le lissage.
 *
 *--------------------------------------------------------------*/

void Curves::smooth(const Curves &curves , int max_frequency)

{
  register int i , j , k;
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


/*--------------------------------------------------------------*
 *
 *  Constructeur par copie de la classe Curves.
 *
 *  arguments : reference sur un objet Curves,
 *              type de transformation ('c' : copie , 's' : lissage),
 *              seuil sur les frequences pour le lissage.
 *
 *--------------------------------------------------------------*/

Curves::Curves(const Curves &curves , char transform , int max_frequency)

{
  if (!(curves.frequency)) {
    transform = 'c';
  }

  switch (transform) {
  case 's' :
    smooth(curves , max_frequency);
    break;
  default :
    copy(curves);
    break;
  }
}


/*--------------------------------------------------------------*
 *
 *  Destruction des champs d'un objet Curves.
 *
 *--------------------------------------------------------------*/

void Curves::remove()

{
  delete [] index_parameter;
  delete [] frequency;

  if (point) {
    register int i;

    for (i = 0;i < nb_curve;i++) {
      delete [] point[i];
    }
    delete [] point;
  }
}


/*--------------------------------------------------------------*
 *
 *  Destructeur de la classe Curves.
 *
 *--------------------------------------------------------------*/

Curves::~Curves()

{
  remove();
}


/*--------------------------------------------------------------*
 *
 *  Operateur d'assignement de la classe Curves.
 *
 *  argument : reference sur un objet Curves.
 *
 *--------------------------------------------------------------*/

Curves& Curves::operator=(const Curves &curves)

{
  if (&curves != this) {
    remove();
    copy(curves);
  }

  return *this;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture de 2 familles de courbes de meme longueur.
 *
 *  arguments : stream, flag fichier, pointeurs sur un objet Curves.
 *
 *--------------------------------------------------------------*/

ostream& Curves::ascii_print(ostream &os , bool file_flag , const Curves *curves) const

{
  register int i , j , k;
  int nb_column = nb_curve + 1 , *width;
  long old_adjust;


  old_adjust = os.setf(ios::right , ios::adjustfield);

  // calcul des largeurs des colonnes

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

  // ecriture des points

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

  os.setf((FMTFLAGS)old_adjust , ios::adjustfield);

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture de 2 familles de courbes de meme longueur au format tableur.
 *
 *  arguments : stream, pointeurs sur un objet Curves.
 *
 *--------------------------------------------------------------*/

ostream& Curves::spreadsheet_print(ostream &os , const Curves *curves) const

{
  register int i , j;


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


/*--------------------------------------------------------------*
 *
 *  Calcul de la longueur des courbes a afficher (sortie Gnuplot).
 *
 *--------------------------------------------------------------*/

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


/*--------------------------------------------------------------*
 *
 *  Ecriture de 2 familles de courbes au format Gnuplot.
 *
 *  arguments : path, longueur des courbes, pointeur sur 2 objets Curves.
 *
 *--------------------------------------------------------------*/

bool Curves::plot_print(const char *path , int ilength ,
                        const Curves *curves_0 , const Curves *curves_1) const

{
  bool status = false;
  register int i , j;
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


/*--------------------------------------------------------------*
 *
 *  Ecriture d'une courbe.
 *
 *  arguments : indice de la courbe, reference sur un objet SinglePlot.
 *
 *--------------------------------------------------------------*/

void Curves::plotable_write(int index , SinglePlot &plot) const

{
  register int i;


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


/*--------------------------------------------------------------*
 *
 *  Ecriture d'une famille de courbes.
 *
 *  argument : reference sur un objet MultiPlot.
 *
 *--------------------------------------------------------------*/

void Curves::plotable_write(MultiPlot &plot) const

{
  register int i , j;


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


/*--------------------------------------------------------------*
 *
 *  Ecriture des frequences.
 *
 *  argument : reference sur un objet SinglePlot.
 *
 *--------------------------------------------------------------*/

void Curves::plotable_frequency_write(SinglePlot &plot) const

{
  register int i;


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


/*--------------------------------------------------------------*
 *
 *  Visualisation d'un objet Curves.
 *
 *  arguments : stream, reference sur un objet Curves.
 *
 *--------------------------------------------------------------*/

ostream& operator<<(ostream &os , const Curves &curves)

{
  register int i , j;


  os.precision(4);

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

  os.precision(6);

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de l'effectif maximum d'un objet Curves
 *
 *--------------------------------------------------------------*/

int Curves::max_frequency_computation() const

{
  int max_frequency = I_DEFAULT;


  if (frequency) {
    register int i;


    max_frequency = 0;
    for (i = offset;i < length;i++) {
      if (frequency[i] > max_frequency) {
        max_frequency = frequency[i];
      }
    }
  }

  return max_frequency;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de l'effectif total d'un objet Curves
 *
 *--------------------------------------------------------------*/

int Curves::nb_element_computation() const

{
  int nb_element = I_DEFAULT;


  if (frequency) {
    register int i;


    nb_element = 0;
    for (i = offset;i < length;i++) {
      nb_element += frequency[i];
    }
  }

  return nb_element;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la moyenne d'une courbe.
 *
 *  argument : indice de la courbe.
 *
 *--------------------------------------------------------------*/

double Curves::mean_computation(int index) const

{
  register int i;
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


/*--------------------------------------------------------------*
 *
 *  Calcul de la variation totale.
 *
 *  arguments : indice de la courbe, moyenne.
 *
 *--------------------------------------------------------------*/

double Curves::total_square_sum_computation(int index , double mean) const

{
  register int i;
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


};  // namespace stat_tool
