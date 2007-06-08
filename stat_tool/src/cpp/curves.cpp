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



#include <iomanip>
// #include <rw/vstream.h>
// #include <rw/rwfile.h>
#include "stat_tools.h"
#include "curves.h"
#include "stat_label.h"

#include "tool/config.h"

using namespace std;


extern int column_width(int value);
extern int column_width(int nb_value , const double *value , double scale = 1.);



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
  frequency = 0;
  point = 0;
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Curves.
 *
 *  arguments : nombre de courbes, longueur des courbes,
 *              flag sur les frequences, flag initialisation.
 *
 *--------------------------------------------------------------*/

Curves::Curves(int inb_curve , int ilength , bool frequency_flag , bool init_flag)

{
  register int i , j;
  int *pfrequency;
  double *ppoint;


  nb_curve = inb_curve;
  length = ilength;
  offset = 0;

  if (frequency_flag) {
    frequency = new int[length];

    if (init_flag) {
      pfrequency = frequency;
      for (i = 0;i < length;i++) {
        *pfrequency++ = 0;
      }
    }
  }

  else {
    frequency = 0;
  }

  point = new double*[nb_curve];
  for (i = 0;i < nb_curve;i++) {
    point[i] = new double[length];

    if (init_flag) {
      ppoint = point[i];
      for (j = 0;j < length;j++) {
        *ppoint++ = 0.;
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
  double norm , *pmass , *pcumul;


  nb_curve = 2;

  pcumul = dist.cumul + dist.nb_value - 1;
  length = dist.nb_value;
  while ((length > 1) && (*pcumul == *(pcumul - 1))) {
    pcumul--;
    length--;
  }

  offset = 0;
  frequency = 0;

  point = new double*[2];
  for (i = 0;i < 2;i++) {
    point[i] = new double[length];
  }

  pmass = dist.mass;
  pcumul = dist.cumul;
  norm = 1. - dist.complement;

  for (i = 0;i < length;i++) {
    if ((*pmass <= norm) && ((1. - dist.complement - *pcumul) <= norm)) {
      point[0][i] = *pmass++ / norm;
      point[1][i] = (1. - dist.complement - *pcumul) / norm;
      norm = 1. - dist.complement - *pcumul++;
    }

    else {
      break;
    }
  }

  length = i;
}


/*--------------------------------------------------------------*
 *
 *  Construction d'un objet Curves a partir d'un objet Histogram.
 *
 *  argument : reference sur un objet Histogram.
 *
 *--------------------------------------------------------------*/

Curves::Curves(const Histogram &histo)

{
  register int i;
  int norm , *hfrequency;


  nb_curve = 2;
  length = histo.nb_value;
  offset = 0;

  point = new double*[2];
  for (i = 0;i < 2;i++) {
    point[i] = new double[length];
  }
  frequency = new int[length];

  hfrequency = histo.frequency;
  norm = histo.nb_element;

  for (i = 0;i < length;i++) {
    point[0][i] = (double)*hfrequency / (double)norm;
    point[1][i] = (double)(norm - *hfrequency) / (double)norm;
    frequency[i] = norm;
    norm -= *hfrequency++;
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
  int *pfrequency , *cfrequency;
  double *ppoint , *cpoint;


  nb_curve = curves.nb_curve;
  length = curves.length;
  offset = curves.offset;

  if (curves.frequency) {
    frequency = new int[length];

    pfrequency = frequency;
    cfrequency = curves.frequency;
    for (i = 0;i < length;i++) {
      *pfrequency++ = *cfrequency++;
    }
  }

  else {
    frequency = 0;
  }

  point = new double*[nb_curve];
  for (i = 0;i < nb_curve;i++) {
    point[i] = new double[length];

    ppoint = point[i];
    cpoint = curves.point[i];
    for (j = 0;j < length;j++) {
      *ppoint++ = *cpoint++;
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
  int range , buff , min , max , total , *pfrequency;
  double sum , *ppoint;


  nb_curve = curves.nb_curve;
  length = curves.length;
  offset = curves.offset;

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

      pfrequency = curves.frequency + min;
      total = 0;
      for (j = min;j <= max;j++) {
        total += *pfrequency++;
      }

      if (total > 0) {
        for (j = 0;j < nb_curve;j++) {
          ppoint = curves.point[j] + min;
          pfrequency = curves.frequency + min;
          sum = 0.;
          for (k = min;k <= max;k++) {
            sum += *ppoint++ * *pfrequency++;
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

  width[0] = column_width(length - 1);

  i = 1;
  for (j = 0;j < nb_curve;j++) {
    if (curves) {
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
    os << setw(width[0]) << i;

    j = 1;
    for (k = 0;k < nb_curve;k++) {
      if (curves) {
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
    os << i;

    for (j = 0;j < nb_curve;j++) {
      if (curves) {
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
      for (j = 0;j < nb_curve;j++) {
        out_file << point[j][i] << " ";
      }

      if (curves_0) {
        for (j = 0;j < nb_curve;j++) {
          out_file << curves_0->point[j][i] << " ";
        }
      }

      if (curves_1) {
        for (j = 0;j < nb_curve;j++) {
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
 *  Fonctions pour la persistance.
 *
 *--------------------------------------------------------------*/

/* RWspace Curves::binaryStoreSize() const

{
  RWspace size = sizeof(nb_curve) + sizeof(length) + sizeof(offset);


  size += sizeof(true);
  if (frequency) {
    size += sizeof(*frequency) * length;
  }

  size += sizeof(**point) * nb_curve * length;

  return size;
}


void Curves::restoreGuts(RWvistream &is)

{
  bool status;
  register int i , j;


  remove();

  is >> nb_curve >> length >> offset;

  is >> status;
  if (status) {
    frequency = new int[length];
    for (i = 0;i < length;i++) {
      is >> frequency[i];
    }
  }
  else {
    frequency = 0;
  }

  point = new double*[nb_curve];
  for (i = 0;i < nb_curve;i++) {
    point[i] = new double[length];
    for (j = 0;j < length;j++) {
      is >> point[i][j];
    }
  }
}


void Curves::restoreGuts(RWFile &file)

{
  bool status;
  register int i;


  remove();

  file.Read(nb_curve);
  file.Read(length);
  file.Read(offset);

  file.Read(status);
  if (status) {
    frequency = new int[length];
    file.Read(frequency , length);
  }
  else {
    frequency = 0;
  }

  point = new double*[nb_curve];
  for (i = 0;i < nb_curve;i++) {
    point[i] = new double[length];
    file.Read(point[i] , length);
  }
}


void Curves::saveGuts(RWvostream &os) const

{
  register int i , j;


  os << nb_curve << length << offset;

  if (frequency) {
    os << true;
    for (i = 0;i < length;i++) {
      os << frequency[i];
    }
  }
  else {
    os << false;
  }

  for (i = 0;i < nb_curve;i++) {
    for (j = 0;j < length;j++) {
      os << point[i][j];
    }
  }
}


void Curves::saveGuts(RWFile &file) const

{
  register int i;


  file.Write(nb_curve);
  file.Write(length);
  file.Write(offset);

  if (frequency) {
    file.Write(true);
    file.Write(frequency , length);
  }
  else {
    file.Write(false);
  }

  for (i = 0;i < nb_curve;i++) {
    file.Write(point[i] , length);
  }
} */


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
    int *pfrequency;


    pfrequency = frequency + offset;
    max_frequency = 0;
    for (i = offset;i < length;i++) {
      if (*pfrequency > max_frequency) {
        max_frequency = *pfrequency;
      }
      pfrequency++;
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
    int *pfrequency;


    pfrequency = frequency + offset;
    nb_element = 0;
    for (i = offset;i < length;i++) {
      nb_element += *pfrequency++;
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
  int nb_element , *pfrequency;
  double mean , *ppoint;


  pfrequency = frequency + offset;
  ppoint = point[index] + offset;
  nb_element = 0;
  mean = 0.;

  for (i = offset;i < length;i++) {
    if (*pfrequency > 0) {
      nb_element += *pfrequency;
      mean += *pfrequency * *ppoint;
    }
    pfrequency++;
    ppoint++;
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
  int *pfrequency;
  double total_square_sum , diff , *ppoint;


  pfrequency = frequency + offset;
  ppoint = point[index] + offset;
  total_square_sum = 0.;

  for (i = offset;i < length;i++) {
    if (*pfrequency > 0) {
      diff = *ppoint - mean;
      total_square_sum += *pfrequency * diff * diff;
    }
    pfrequency++;
    ppoint++;
  }

  return total_square_sum;
}
