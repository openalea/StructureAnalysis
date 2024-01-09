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
 *       $Id: distribution.cpp 17997 2015-04-23 06:55:51Z guedon $
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
#include "curves.h"
#include "stat_label.h"

using namespace std;


namespace stat_tool {



/*--------------------------------------------------------------*
 *
 *  Copie des probabilites de chaque valeur.
 *
 *  arguments : reference sur un objet Distribution, nombre de valeurs.
 *
 *--------------------------------------------------------------*/

void Distribution::mass_copy(const Distribution &dist , int inb_value)

{
  register int i;


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


/*--------------------------------------------------------------*
 *
 *  Copie d'un objet Distribution dans le cas ou le nombre
 *  de valeurs allouees est le meme pour les deux objets.
 *
 *  arguments : reference sur un objet Distribution.
 *
 *--------------------------------------------------------------*/

void Distribution::equal_size_copy(const Distribution &dist)

{
  register int i;


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


/*--------------------------------------------------------------*
 *
 *  Initialisation d'un objet Distribution.
 *
 *  argument : nombre de valeurs.
 *
 *--------------------------------------------------------------*/

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
    register int i;


    mass = new double[nb_value];
    cumul = new double[nb_value];

    for (i = 0;i < nb_value;i++) {
      mass[i] = 0.;
      cumul[i] = 0.;
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Distribution.
 *
 *  argument : nombre de valeurs.
 *
 *--------------------------------------------------------------*/

Distribution::Distribution(int inb_value)

{
  init(inb_value);
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Distribution.
 *
 *  arguments : nombre de valeurs, masses de probabilite.
 *
 *--------------------------------------------------------------*/

Distribution::Distribution(int inb_value , double *imass)

{
  register int i;


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


/*--------------------------------------------------------------*
 *
 *  Construction d'un objet Distribution a partir d'un objet
 *  Distribution initial avec changement d'echelle.
 *
 *  arguments : reference sur un objet Distribution, facteur d'echelle.
 *
 *--------------------------------------------------------------*/

Distribution::Distribution(const Distribution &dist , double scaling_coeff)

{
  register int i , j;
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


/*--------------------------------------------------------------*
 *
 *  Construction d'un objet Distribution a partir d'un objet FrequencyDistribution.
 *
 *  argument : reference sur un objet FrequencyDistribution.
 *
 *--------------------------------------------------------------*/

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


/*--------------------------------------------------------------*
 *
 *  Copie d'un objet Distribution.
 *
 *  arguments : reference sur un objet Distribution, nombre de valeurs allouees.
 *
 *--------------------------------------------------------------*/

void Distribution::copy(const Distribution &dist , int ialloc_nb_value)

{
  register int i;


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


/*--------------------------------------------------------------*
 *
 *  Copie d'un objet Distribution avec renormalisation.
 *
 *  argument : reference sur un objet Distribution.
 *
 *--------------------------------------------------------------*/

void Distribution::normalization_copy(const Distribution &dist)

{
  if (dist.complement == 0.) {
    copy(dist);
  }

  else {
    register int i;


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


/*--------------------------------------------------------------*
 *
 *  Constructeur par copie de la classe Distribution.
 *
 *  arguments : reference sur un objet Distribution, type de transformation
 *              ('c' : copie , 'n' : copie avec renormalisation),
 *              nombre de valeurs allouees.
 *
 *--------------------------------------------------------------*/

Distribution::Distribution(const Distribution &dist , char transform , int ialloc_nb_value)

{
  switch (transform) {
  case 'c' :
    copy(dist , ialloc_nb_value);
    break;
  case 'n' :
    normalization_copy(dist);
    break;
  default :
    copy(dist);
    break;
  }
}


/*--------------------------------------------------------------*
 *
 *  Destructeur de la classe Distribution.
 *
 *--------------------------------------------------------------*/

Distribution::~Distribution()

{
    delete [] mass;
    delete [] cumul;
}


/*--------------------------------------------------------------*
 *
 *  Operateur d'assignement de la classe Distribution.
 *
 *  argument : reference sur un objet Distribution.
 *
 *--------------------------------------------------------------*/

Distribution& Distribution::operator=(const Distribution &dist)

{
  if (&dist != this) {
    delete [] mass;
    delete [] cumul;

    copy(dist);
  }

  return *this;
}


/*--------------------------------------------------------------*
 *
 *  Operateur d'egalite de la classe Distribution.
 *
 *  argument : reference sur un objet Distribution.
 *
 *--------------------------------------------------------------*/

bool Distribution::operator==(const Distribution &dist) const

{
  bool status = true;
  register int i;


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


/*--------------------------------------------------------------*
 *
 *  Ecriture des caracteristiques d'une loi.
 *
 *  arguments : stream, flag ecriture des parametres de forme, flag commentaire.
 *
 *--------------------------------------------------------------*/

ostream& Distribution::ascii_characteristic_print(ostream &os , bool shape , bool comment_flag) const

{
  if ((mean != D_DEFAULT) && (variance != D_DEFAULT)) {
    if (comment_flag) {
      os << "# ";
    }
    os << STAT_label[STATL_MEAN] << ": " << mean << "   "
       << STAT_label[STATL_VARIANCE] << ": " << variance << "   "
       << STAT_label[STATL_STANDARD_DEVIATION] << ": " << sqrt(variance) << endl;

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


/*--------------------------------------------------------------*
 *
 *  Calcul de la largeur d'une colonne d'entiers.
 *
 *  argument : valeur.
 *
 *--------------------------------------------------------------*/

int column_width(int value)

{
  ostringstream ostring;
  ostring << value;
  return (ostring.str()).size();
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la largeur d'une colonne d'entiers.
 *
 *  arguments : valeur minimum, valeur maximum.
 *
 *--------------------------------------------------------------*/

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


/*--------------------------------------------------------------*
 *
 *  Calcul de la largeur d'une colonne de reels.
 *
 *  arguments : nombre de valeurs, pointeur sur les valeurs reelles,
 *              facteur d'echelle.
 *
 *--------------------------------------------------------------*/

int column_width(int nb_value , const double *value , double scale)

{
  register int i;
  int width , max_width = 0;


  for (i = 0;i < nb_value;i++) {
    ostringstream ostring;
    ostring << *value++ * scale;
    width = (ostring.str()).size();
    if (width > max_width) {
      max_width = width;
    }
  }

  return max_width;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'une loi et d'une loi empirique.
 *
 *  arguments : stream, flag commentaire, flags sur l'ecriture de la fonction
 *              de repartition et sur le calcul du nombre de valeurs,
 *              pointeur sur une loi empirique.
 *
 *--------------------------------------------------------------*/

ostream& Distribution::ascii_print(ostream &os , bool comment_flag , bool cumul_flag ,
                                   bool nb_value_flag , const FrequencyDistribution *histo) const

{
  register int i;
  int ascii_nb_value = nb_value , width[5];
  long old_adjust;
  double scale , *histo_cumul , *pcumul;


  old_adjust = os.setf(ios::right , ios::adjustfield);

  // calcul du facteur d'echelle, de la fonction de repartition deduite de
  // la loi empirique et du nombre de valeurs

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

  // calcul des largeurs des colonnes

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

  // ecriture des probabilites de chaque valeur

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

  os.setf((FMTFLAGS)old_adjust , ios::adjustfield);

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'une famille de lois et d'une loi empirique.
 *
 *  arguments : stream, nombre de lois, pointeurs sur les lois,
 *              facteurs d'echelle, flag commentaire,
 *              flag sur l'ecriture de la fonction de repartition,
 *              pointeur sur une loi empirique, flag sur l'ordre des lois.
 *
 *--------------------------------------------------------------*/

ostream& Distribution::ascii_print(ostream &os , int nb_dist , const Distribution **dist ,
                                   double *dist_scale , bool comment_flag , bool cumul_flag ,
                                   const FrequencyDistribution *histo , bool mass_first) const

{
  register int i , j;
  int ascii_nb_value = nb_value , *width;
  long old_adjust;
  double scale , *histo_cumul;


  old_adjust = os.setf(ios::right , ios::adjustfield);

  width = new int[nb_dist + 5];

  // calcul du facteur d'echelle, de la fonction de repartition deduite de
  // la loi empirique et du nombre de valeurs

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

  // calcul des largeurs des colonnes

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

  // ecriture des probabilites et de la fonction de repartition

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

  os.setf((FMTFLAGS)old_adjust , ios::adjustfield);

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture des caracteristiques d'une loi au format tableur.
 *
 *  arguments : stream, flag ecriture des parametres de forme.
 *
 *--------------------------------------------------------------*/

ostream& Distribution::spreadsheet_characteristic_print(ostream &os , bool shape) const

{
  if ((mean != D_DEFAULT) && (variance != D_DEFAULT)) {
    os << STAT_label[STATL_MEAN] << "\t" << mean << "\t"
       << STAT_label[STATL_VARIANCE] << "\t" << variance << "\t"
       << STAT_label[STATL_STANDARD_DEVIATION] << "\t" << sqrt(variance) << endl;

    if ((shape) && (variance > 0.)) {
      os << STAT_label[STATL_SKEWNESS_COEFF] << "\t" << skewness_computation() << "\t"
         << STAT_label[STATL_KURTOSIS_COEFF] << "\t" << kurtosis_computation() << endl;
    }
  }

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'une loi et d'une loi empirique au format tableur.
 *
 *  arguments : stream, flags sur l'ecriture de la fonction de repartition et
 *              de la fonction de concentration, flag sur le calcul du nombre de valeurs,
 *              pointeur sur une loi empirique.
 *
 *--------------------------------------------------------------*/

ostream& Distribution::spreadsheet_print(ostream &os , bool cumul_flag ,
                                         bool concentration_flag , bool nb_value_flag ,
                                         const FrequencyDistribution *histo) const

{
  register int i;
  int spreadsheet_nb_value = nb_value;
  double scale , *pcumul , *histo_cumul , *concentration , *histo_concentration;


  if ((!cumul_flag) || (mean == 0.)) {
    concentration_flag = false;
  }

  // calcul du facteur d'echelle, de la fonction de repartition deduite de
  // la loi empirique, des fonctions de concentration et du nombre de valeurs

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

  // ecriture des probabilites, des fonctions de repartition et de concentration

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


/*--------------------------------------------------------------*
 *
 *  Ecriture d'une famille de lois et d'une loi empirique au format tableur.
 *
 *  arguments : stream, nombre de lois, pointeurs sur les lois,
 *              facteurs d'echelle, flag sur l'ecriture de la fonction
 *              de repartition, pointeur sur une loi empirique,
 *              flag sur l'ordre des lois.
 *
 *--------------------------------------------------------------*/

ostream& Distribution::spreadsheet_print(ostream &os , int nb_dist , const Distribution **dist ,
                                         double *dist_scale , bool cumul_flag ,
                                         const FrequencyDistribution *histo , bool mass_first) const

{
  register int i , j;
  int spreadsheet_nb_value = nb_value;
  double scale , *histo_cumul;


  // calcul du facteur d'echelle, de la fonction de repartition deduite de
  // la loi empirique et du nombre de valeurs

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

  // ecriture des probabilites de chaque valeur

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


/*--------------------------------------------------------------*
 *
 *  Calcul du nombre de valeurs a afficher (sortie Gnuplot).
 *
 *  argument : pointeur sur une loi empirique.
 *
 *--------------------------------------------------------------*/

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


/*--------------------------------------------------------------*
 *
 *  Ecriture d'une loi au format Gnuplot.
 *
 *  arguments : path, pointeur sur la fonction de concentration,
 *              facteur d'echelle.
 *
 *--------------------------------------------------------------*/

bool Distribution::plot_print(const char *path , double *concentration ,
                              double scale) const

{
  bool status = false;
  register int i;
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


/*--------------------------------------------------------------*
 *
 *  Ecriture d'une loi et d'une loi empirique au format Gnuplot.
 *
 *  arguments : path, pointeur sur une loi empirique.
 *
 *--------------------------------------------------------------*/

bool Distribution::plot_print(const char *path , const FrequencyDistribution *histo) const

{
  bool status = false;
  register int i;
  double scale;
  ofstream out_file(path);


  if (out_file) {
    status = true;

    // calcul du facteur d'echelle

    if (histo) {
      scale = histo->nb_element / (1. - complement);
    }
    else {
      scale = 1.;
    }

    // ecriture des frequences et des probabilites de chaque valeur

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


/*--------------------------------------------------------------*
 *
 *  Ecriture d'une famille de lois et d'une famille de lois empiriques
 *  au format Gnuplot.
 *
 *  arguments : path, nombre de lois, pointeurs sur les lois,
 *              facteurs d'echelle, nombre de valeurs des lois,
 *              nombre de lois empiriques, pointeurs sur les lois empiriques.
 *
 *--------------------------------------------------------------*/

bool plot_print(const char *path , int nb_dist , const Distribution **dist ,
                double *scale , int *dist_nb_value , int nb_histo ,
                const FrequencyDistribution **histo)

{
  bool status = false;
  register int i , j;
  int plot_nb_value = 0;
  ofstream out_file(path);


  if (out_file) {
    status = true;

    // calcul du nombre de valeurs

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

    // ecriture des frequences et des probabilites de chaque valeur

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


/*--------------------------------------------------------------*
 *
 *  Ecriture d'une famille de fonctions de repartition en vue du matching
 *  au format Gnuplot
 *
 *  arguments : path, nombre de fonctions de repartition,
 *              pointeurs sur les plus petites valeurs, les nombres de valeurs et
 *              les fonctions de repartition.
 *
 *--------------------------------------------------------------*/

bool cumul_matching_plot_print(const char *path , int nb_cumul , int *offset ,
                               int *nb_value , double **cumul)

{
  bool status = false;
  register int i , j;
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


/*--------------------------------------------------------------*
 *
 *  Extraction du label d'un nom de fichier.
 *
 *  argument : nom d'un fichier.
 *
 *--------------------------------------------------------------*/

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


/*--------------------------------------------------------------*
 *
 *  Sortie Gnuplot d'une famille de lois :
 *  - lois et fonctions de repartition,
 *  - mise en correspondance des fonctions de repartition,
 *  - courbes de concentration.
 *
 *  arguments : reference sur un objet StatError, prefixe des fichiers,
 *              nombre de lois, pointeurs sur les lois, titre des figures.
 *
 *--------------------------------------------------------------*/

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
    register int i , j , k;
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

    // ecriture des fichiers de donnees

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

      // calcul des fonctions de concentration

      concentration[i] = dist[i]->concentration_function_computation();

      // calcul du nombre de valeurs maximum, de l'etendue maximum et
      // de la probabilite maximum

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

        // lois

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

          // fonctions de repartition

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

          // mise en correspondance des fonctions de repartition en prenant
          // comme reference la loi dont l'etendue est la plus grande

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

          // courbes de concentration

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


/*--------------------------------------------------------------*
 *
 *  Ecriture d'une loi.
 *
 *  arguments : reference sur un objet SinglePlot,
 *              facteur d'echelle (valeur par defaut : 1).
 *
 *--------------------------------------------------------------*/

void Distribution::plotable_mass_write(SinglePlot &plot , double scale) const

{
  register int i;


  for (i = MAX(offset - 1 , 0);i < nb_value;i++) {
    plot.add_point(i , mass[i] * scale);
  }
  if ((cumul[nb_value - 1] > 1. - DOUBLE_ERROR) &&
      (mass[nb_value - 1] > PLOT_MASS_THRESHOLD)) {
    plot.add_point(nb_value , 0.);
  }
}


/*--------------------------------------------------------------*
 *
 *  Ecriture de la fonction de repartition d'une loi.
 *
 *  argument : reference sur un objet SinglePlot.
 *
 *--------------------------------------------------------------*/

void Distribution::plotable_cumul_write(SinglePlot &plot) const

{
  register int i;


  for (i = MAX(offset - 1 , 0);i < nb_value;i++) {
    plot.add_point(i , cumul[i]);
  }
}


/*--------------------------------------------------------------*
 *
 *  Ecriture de la mise en correspondance d'une fonction de repartition avec
 *  la fonction de repartition d'une loi de reference.
 *
 *  argument : reference sur un objet SinglePlot et sur la loi de reference.
 *
 *--------------------------------------------------------------*/

void Distribution::plotable_cumul_matching_write(SinglePlot &plot ,
                                                 const Distribution &reference_dist) const

{
  register int i;


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


/*--------------------------------------------------------------*
 *
 *  Ecriture de la courbe de concentration d'une loi.
 *
 *  argument : reference sur un objet SinglePlot.
 *
 *--------------------------------------------------------------*/

void Distribution::plotable_concentration_write(SinglePlot &plot) const

{
  register int i;
  double *concentration;


  concentration = concentration_function_computation();

  plot.add_point(0. , 0.);
  for (i = offset;i < nb_value;i++) {
    plot.add_point(cumul[i] , concentration[i]);
  }

  delete [] concentration;
}


/*--------------------------------------------------------------*
 *
 *  Sortie graphique d'une famille de lois :
 *  - lois et fonctions de repartition,
 *  - mise en correspondance des fonctions de repartition,
 *  - courbes de concentration.
 *
 *  arguments : reference sur un objet StatError,
 *              nombre de lois, pointeurs sur les lois.
 *
 *--------------------------------------------------------------*/

MultiPlotSet* Distribution::get_plotable_distributions(StatError &error , int nb_dist ,
                                                       const Distribution **idist) const

// usage : plot[vue][distribution]
{
  MultiPlotSet *plot_set;


  error.init();

  if (nb_dist > PLOT_NB_DISTRIBUTION) {
    plot_set = NULL;
    error.update(STAT_error[STATR_PLOT_NB_DISTRIBUTION]);
  }

  else {
    register int i , j , k;
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

      // calcul du nombre de valeurs maximum, de l'etendue maximum et
      // de la probabilite maximum

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

    // 1ere vue : lois

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

      // 2eme vue : fonctions de repartition

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

      // 3eme vue : mise en correspondance des fonctions de repartition en prenant
      // comme reference la loi dont l'etendue est la plus grande

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

      // 4eme vue : courbes de concentration

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


/*--------------------------------------------------------------*
 *
 *  Sortie graphique d'une loi
 *
 *--------------------------------------------------------------*/

MultiPlotSet* Distribution::get_plotable() const

{
  StatError error;

  return get_plotable_distributions(error , 0 , NULL);
}


/*--------------------------------------------------------------*
 *
 *  Calcul des taux de survie a partir d'une loi et ecriture du resultat.
 *
 *  argument : stream.
 *
 *--------------------------------------------------------------*/

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


/*--------------------------------------------------------------*
 *
 *  Calcul des taux de survie a partir d'une loi et
 *  ecriture du resultat dans un fichier.
 *
 *  arguments : reference sur un objet StatError, path.
 *
 *--------------------------------------------------------------*/

bool Distribution::survival_ascii_write(StatError &error , const char *path) const

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
 *  Calcul des taux de survie a partir d'une loi et ecriture
 *  du resultat dans un fichier au format tableur.
 *
 *  arguments : reference sur un objet StatError, path.
 *
 *--------------------------------------------------------------*/

bool Distribution::survival_spreadsheet_write(StatError &error , const char *path) const

{
  bool status;
  ofstream out_file(path);
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


/*--------------------------------------------------------------*
 *
 *  Ecriture d'une loi et de la fonction de survie associee
 *  au format Gnuplot.
 *
 *  arguments : path, pointeur sur la fonction de survie.
 *
 *--------------------------------------------------------------*/

bool Distribution::survival_plot_print(const char *path , double *survivor) const

{
  bool status = false;
  register int i;
  ofstream out_file(path);


  if (out_file) {
    status = true;

    for (i = 0;i < nb_value;i++) {
      out_file << mass[i] << " " << survivor[i] << endl;
    }
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Calcul des taux de survie a partir d'une loi et sortie Gnuplot du resultat.
 *
 *  arguments : reference sur un objet StatError, prefixe des fichiers,
 *              titre des figures.
 *
 *--------------------------------------------------------------*/

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
    register int i;
    double *survivor;
    Curves *survival_rate;
    ostringstream data_file_name[2];


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

        // loi et fonction de survie

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

        // taux de survie

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
 *  Calcul et ecriture de la fonction de survie d'une loi.
 *
 *  argument : reference sur un objet SinglePlot.
 *
 *--------------------------------------------------------------*/

void Distribution::plotable_survivor_write(SinglePlot &plot) const

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
 *  Calcul des taux de survie a partir d'une loi et sortie graphique du resultat.
 *
 *  argument : reference sur un objet StatError.
 *
 *--------------------------------------------------------------*/

MultiPlotSet* Distribution::survival_get_plotable(StatError &error) const

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


    plot_set = new MultiPlotSet(2);
    MultiPlotSet &plot = *plot_set;

    plot.title = "Survival analysis";
    plot.border = "15 lw 0";

    // 1ere vue : loi et fonction de survie

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

    // 2eme vue : taux de survie

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


/*--------------------------------------------------------------*
 *
 *  Visualisation d'une loi.
 *
 *  argument : stream.
 *
 *--------------------------------------------------------------*/

ostream& Distribution::print(ostream &os) const

{
  register int i;


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


/*--------------------------------------------------------------*
 *
 *  Visualisation d'une distribution.
 *
 *  arguments : stream, reference sur un objet Distribution.
 *
 *--------------------------------------------------------------*/

ostream& operator<<(ostream &os , const Distribution &dist)

{
  os.precision(5);

  os << endl;
  dist.print(os);

  os.precision(6);

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Calcul du nombre de valeurs prises par une v.a..
 *
 *--------------------------------------------------------------*/

void Distribution::nb_value_computation()

{
  double *pmass;


  pmass = mass + alloc_nb_value;
  nb_value = alloc_nb_value;

  while ((*--pmass == 0.) && (nb_value > 0)) {
    nb_value--;
  }
}


/*--------------------------------------------------------------*
 *
 *  Calcul du nombre de valeurs de probabilite nulle a partir de 0.
 *
 *--------------------------------------------------------------*/

void Distribution::offset_computation()

{
  double *pmass;


  pmass = mass;
  offset = 0;

  while ((*pmass++ == 0.) && (offset < nb_value - 1)) {
    offset++;
  }
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la valeur de probabilite maximum.
 *
 *--------------------------------------------------------------*/

void Distribution::max_computation()

{
  register int i;


  max = 0.;
  for (i = offset;i < nb_value;i++) {
    if (mass[i] > max) {
      max = mass[i];
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la moyenne d'une loi.
 *
 *--------------------------------------------------------------*/

void Distribution::mean_computation()

{
  if (cumul[nb_value - 1] > 0.) {
    register int i;


    mean = 0.;
    for (i = offset;i < nb_value;i++) {
      mean += mass[i] * i;
    }
    mean /= cumul[nb_value - 1];
  }
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la variance d'une loi.
 *
 *--------------------------------------------------------------*/

void Distribution::variance_computation()

{
  if (mean != D_DEFAULT) {
    register int i;
    double diff;


    variance = 0.;
    for (i = offset;i < nb_value;i++) {
      diff = i - mean;
      variance += mass[i] * diff * diff;
    }
    variance /= cumul[nb_value - 1];

    if (variance < 0.) {
      variance = 0.;
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Calcul de l'ecart absolu moyen d'une loi.
 *
 *--------------------------------------------------------------*/

double Distribution::mean_absolute_deviation_computation() const

{
  register int i;
  double mean_absolute_deviation = D_DEFAULT;


  if (mean != D_DEFAULT) {
    mean_absolute_deviation = 0.;
    for (i = offset;i < nb_value;i++) {
      mean_absolute_deviation += mass[i] * fabs(i - mean);
    }
    mean_absolute_deviation /= cumul[nb_value - 1];

    if (mean_absolute_deviation < 0.) {
      mean_absolute_deviation = 0.;
    }
  }

  return mean_absolute_deviation;
}


/*--------------------------------------------------------------*
 *
 *  Calcul du coefficient d'asymetrie d'une loi.
 *
 *--------------------------------------------------------------*/

double Distribution::skewness_computation() const

{
  register int i;
  double skewness = D_INF , diff;


  if ((mean != D_DEFAULT) && (variance != D_DEFAULT)) {
    skewness = 0.;

    if (variance > 0.) {
      for (i = offset;i < nb_value;i++) {
        diff = i - mean;
        skewness += mass[i] * diff * diff * diff;
      }
      skewness /= (cumul[nb_value - 1] * pow(variance , 1.5));
    }
  }

  return skewness;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de l'exces d'applatissement d'une loi :
 *  exces d'applatissement = coefficient d'applatissement - 3.
 *
 *--------------------------------------------------------------*/

double Distribution::kurtosis_computation() const

{
  register int i;
  double kurtosis = D_INF , diff;


  if ((mean != D_DEFAULT) && (variance != D_DEFAULT)) {
    if (variance == 0.) {
      kurtosis = -2.;
    }

    else {
      kurtosis = 0.;
      for (i = offset;i < nb_value;i++) {
        diff = i - mean;
        kurtosis += mass[i] * diff * diff * diff * diff;
      }
      kurtosis = kurtosis / (cumul[nb_value - 1] * variance * variance) - 3.;
    }
  }

  return kurtosis;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la quantite d'information d'une loi.
 *
 *--------------------------------------------------------------*/

double Distribution::information_computation() const

{
  register int i;
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


/*--------------------------------------------------------------*
 *
 *  Calcul de la somme des carres des differences premieres.
 *
 *--------------------------------------------------------------*/

double Distribution::first_difference_norm_computation() const

{
  register int i;
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


/*--------------------------------------------------------------*
 *
 *  Calcul de la somme des carres des differences secondes.
 *
 *--------------------------------------------------------------*/

double Distribution::second_difference_norm_computation() const

{
  register int i;
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


/*--------------------------------------------------------------*
 *
 *  Calcul de la fonction de repartition d'une loi discrete.
 *
 *  arguments : nombre de valeurs, pointeurs sur les probabilites et
 *              sur la fonction de repartition correspondante.
 *
 *--------------------------------------------------------------*/

void cumul_computation(int nb_value , const double *pmass , double *pcumul)

{
  register int i;


  *pcumul = *pmass;
  for (i = 1;i < nb_value;i++) {
    pcumul++;
    *pcumul = *(pcumul - 1) + *++pmass;
  }
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la fonction de repartition d'une loi discrete.
 *
 *--------------------------------------------------------------*/

void Distribution::cumul_computation()

{
  register int i;


  for (i = 0;i < offset;i++) {
    cumul[i] = 0.;
  }
  stat_tool::cumul_computation(nb_value - offset , mass + offset , cumul + offset);
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la fonction de survie d'une loi discrete.
 *
 *--------------------------------------------------------------*/

double* Distribution::survivor_function_computation() const

{
  register int i;
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


/*--------------------------------------------------------------*
 *
 *  Calcul de la fonction de concentration d'une loi discrete.
 *
 *--------------------------------------------------------------*/

double* Distribution::concentration_function_computation() const

{
  register int i;
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


/*--------------------------------------------------------------*
 *
 *  Calcul du coefficient de concentration d'une loi discrete.
 *
 *--------------------------------------------------------------*/

double Distribution::concentration_computation() const

{
  register int i;
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


/*--------------------------------------------------------------*
 *
 *  Calcul de la distance entre deux lois discretes (1 - recouvrement entre les 2 lois).
 *
 *--------------------------------------------------------------*/

double Distribution::overlap_distance_computation(const Distribution &dist) const

{
  register int i;
  double overlap;


  overlap = 0.;
  for (i = MAX(offset , dist.offset);i < MIN(nb_value , dist.nb_value);i++) {
    overlap += MIN(mass[i] , dist.mass[i]);
  }

  return (1. - overlap);
}


/*--------------------------------------------------------------*
 *
 *  Calcul des logarithmes des probabilites de chaque valeur.
 *
 *  arguments : nombre de valeurs, pointeurs sur les probabilites et
 *              sur la fonction transformee par log correspondante.
 *
 *--------------------------------------------------------------*/

void log_computation(int nb_value , const double *pmass , double *plog)

{
  register int i;


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


/*--------------------------------------------------------------*
 *
 *  Calcul des logarithmes des probabilites de chaque valeur.
 *
 *--------------------------------------------------------------*/

void Distribution::log_computation()

{
  stat_tool::log_computation(nb_value , mass , cumul);
}


};  // namespace stat_tool
