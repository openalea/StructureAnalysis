/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *
 *       Copyright 2005-2009 UMR DAP
 *
 *       File author(s): F. Chaubert-Pereira (chaubert@cirad.fr)
 *
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
#include <iostream>
#include <iomanip>
#include "stat_tool/stat_label.h"
#include "tool/util_math.h"
#include "continuous_histo.h"
#include "stat_tool/stat_tools.h"
#include "tool/config.h"

using namespace std;

/*-------------------------------------------------
 *
 * Constructeur de la classe Continuous_histo
 *
 *-------------------------------------------------*/

Continuous_histo::Continuous_histo()
{
  mean = 0;
  variance = 0;
  frequency = 0;
  nb_element = 0;
  nb_class = 0;
  min_value = 0;
  max_value = 0;
}


/*-------------------------------------------------
 *
 * Constructeur de la classe Continuous_histo
 *
 * arguments : nombre de classes
 *
 *-------------------------------------------------*/

Continuous_histo::Continuous_histo(int inb_class)
{
  nb_element = 0;
  nb_class = inb_class;
  min_value = 0;
  max_value = 0;
  mean = D_DEFAULT;
  variance = D_DEFAULT;

  if (nb_class == 0) {
    frequency = 0;
  }

  else {
    register int i;
    double *pfrequency;

    frequency = new double[nb_class];
    pfrequency = frequency;
    for (i = 0;i < nb_class;i++) {
      *pfrequency++ = 0;
    }
  }
}


/*-------------------------------------------------
 *
 * Constructeur de la classe Continuous_histo
 *
 * arguments : nombre de classes, effectifs par classes
 *
 *-------------------------------------------------*/

Continuous_histo::Continuous_histo(int inb_element, double *pelement)
{
  register int i;
  register int j;
  double tmp1, tmp2;

  double *pfrequency;

  nb_element = inb_element;

  min_value = pelement[0];
  max_value = pelement[0];
  for (i = 1;i < nb_element; i++) {
    if (min_value > pelement[i]) {
      min_value = pelement[i];
    }
    if (max_value < pelement[i]) {
      max_value = pelement[i];
    }
  }

  nb_class_computation();	

  frequency = new double[nb_class];
  pfrequency = frequency;

  for (i = 0;i < nb_class;i++) {
    frequency[i] = 0.;
    for (j = 0;j < nb_element; j++) {
      tmp1 = min_value + i * (max_value - min_value) / nb_class;
      tmp2 = min_value + (i + 1) * (max_value - min_value) / nb_class;
      if ((pelement[j] >= tmp1) && (pelement[j] < tmp2)){
	frequency[i] = frequency[i] + 1.;
      }
      if ((i == (nb_class - 1)) && (pelement[j] == tmp2)) {
	frequency[i] = frequency[i] + 1.; 
      }
    }
  }

  // calcul des caracteristiques de l'histogramme
  mean_computation();
  variance_computation();
}


/*-------------------------------------------------
 *
 * Constructeur par copie de la classe Continuous_histo
 *
 * arguments : reference sur un objet Conitnuous_histo
 *
 *-------------------------------------------------*/

Continuous_histo::Continuous_histo(const Continuous_histo  &cont_histo)
{
  copy(cont_histo);
}


/*-------------------------------------------------
 *
 * Constructeur de la classe Continuous_histo
 *
 * arguments : reference sur un objet Conitnuous_histo
 *
 *-------------------------------------------------*/

Continuous_histo::Continuous_histo(const Continuous_distribution &cont_dist)
{
  register int i;
  double *pdensity;
  pdensity = cont_dist.density;

  nb_element = cont_dist.nb_value;
  nb_class_computation();

  while (( nb_element % nb_class) != 0) {
    nb_class = nb_class+1;
  }
  
  min_value = cont_dist.min_value;
  max_value = cont_dist.max_value;
  
  mean = cont_dist.mean;
  variance = cont_dist.variance;

  frequency = new double[nb_class];

  for (i = 0; i < nb_class; i++) {
    frequency[i] = (pdensity[(i+1) * nb_element / nb_class] + pdensity[i * nb_element/nb_class]) / 2.;
  }
}


/*----------------------------------------------------------------------
 *
 * Constructeur: Passage d'un histogramme des effectifs 
 * à un histogramme des fréquences dont l'aire est égale à X
 * Pour avoir celui des effectifs: 
 * area = nb_element * (max_value-min_value)/nb_class
 *
 * arguments : reference sur un objet ontinuous_histo, nombre
 *
 *----------------------------------------------------------------------*/

Continuous_histo::Continuous_histo(const Continuous_histo &cont_histo, double area) 
{ 
  register int i;
  double step, tmp, sum = 0.;
  double *pfrequency , *cfrequency;

  nb_element = cont_histo.nb_element;
  nb_class = cont_histo.nb_class;
  max_value = cont_histo.max_value;
  min_value = cont_histo.min_value;
  mean = cont_histo.mean;
  variance = cont_histo.variance;

  step = (max_value - min_value) / nb_class;

  frequency = new double[nb_class];

  pfrequency = frequency;
  cfrequency = cont_histo.frequency;
  for (i = 0; i < nb_class;i++) {
    sum = sum + cfrequency[i];
  }

  for (i = 0;i < nb_class;i++) {
    *pfrequency++ = (*cfrequency++ * area) / (sum * step);
  }
}


/*-------------------------------------------------
 *
 * Constructeur de la classe Continuous_histo
 * (changement d'abscisse)
 *
 * arguments : décalage, reference sur un objet Conitnuous_histo
 *
 *-------------------------------------------------*/

Continuous_histo::Continuous_histo(double param, const Continuous_histo &cont_histo)
{
    shift(cont_histo , param);
}


/*-------------------------------------------------
 *
 * Destructeur de la classe Continuous_histo
 *
 *-------------------------------------------------*/

Continuous_histo::~Continuous_histo()
{
  if (frequency) {
    delete [] frequency;
  }
}


/*--------------------------------------------------
 *
 * Calcul du nombre de classes
 *
 *--------------------------------------------------*/

void Continuous_histo::nb_class_computation()
{
  nb_class = int(round(1+(10/3)*log10(nb_element)));
}


/*--------------------------------------------------
 *
 * Calcul de la moyenne d'un histogramme
 *
 *--------------------------------------------------*/

void Continuous_histo::mean_computation()
{
   if (nb_element > 0) {
    register int i;
    double tmp, sum = 0.;
    double *pfrequency;

    pfrequency = frequency ;
    for (i = 0;i < nb_class;i++) {
      sum = sum +  pfrequency[i];
    }
    mean = 0.;

    for (i = 0;i < nb_class;i++) {
      tmp = min_value + i * (max_value - min_value) / nb_class + (max_value - min_value) / (2. * nb_class);
      mean += *pfrequency++ * tmp;
     }

    mean /= sum;
  }
}


/*--------------------------------------------------
 *
 * Calcul de la variance d'un histogramme
 *
 *--------------------------------------------------*/

void Continuous_histo::variance_computation()
{
  variance = 0.;

  if (nb_element > 1) {
    register int i;
    double *pfrequency;
    double diff, tmp, sum = 0.;

    pfrequency = frequency ;
    for (i = 0;i < nb_class;i++) {
      sum = sum + pfrequency[i];
    }
    for (i = 0;i < nb_class;i++) {
      tmp = min_value + i * (max_value - min_value) / nb_class + (max_value - min_value) / (2. * nb_class);
      diff = tmp - mean;
      variance += *pfrequency++ * diff * diff;
    }
    
    variance /= sum;  
  }
}


/*--------------------------------------------------
 *
 * Calcul de l'effectif total de l'histogramme
 *
 *--------------------------------------------------*/
void Continuous_histo::nb_element_computation()
{
  register int i;
  double *pfrequency;

  pfrequency = frequency;
  nb_element = 0;

  for (i = 0;i < nb_class;i++) {
    nb_element += (int)*pfrequency++;
  }
}


/*--------------------------------------------------------------
 *
 * Copie d'un objet Continuous_histo.
 *
 * arguments : reference sur un objet Continuous_histo
 *
 *--------------------------------------------------------------*/

void Continuous_histo::copy(const Continuous_histo &cont_histo)
{
  register int i;
  double *pfrequency , *cfrequency;

  nb_element = cont_histo.nb_element;
  nb_class = cont_histo.nb_class;
  max_value = cont_histo.max_value;
  min_value = cont_histo.min_value;
  mean = cont_histo.mean;
  variance = cont_histo.variance;

  // copie des frequences

  frequency = new double[nb_class];
  pfrequency = frequency;
  cfrequency = cont_histo.frequency;
  for (i = 0;i < nb_class;i++) {
    *pfrequency++ = *cfrequency++;
  }
}


/*--------------------------------------------------------------
 *
 * Operateur d'assignement de la classe Continuous_histo. 
 *
 * arguments : reference sur un objet Continuous_histo
 *
 *--------------------------------------------------------------*/

Continuous_histo& Continuous_histo::operator=(const Continuous_histo &cont_histo)
{
  if (&cont_histo != this) {
    delete [] frequency;

    copy(cont_histo);
  }

  return *this;
}


/*---------------------------------------------------------------
 *
 * Ecriture des caractéristiques d'un objet Continuous_histo
 *
 * arguments : stream, flag de commentaire
 *
 *---------------------------------------------------------------*/

ostream& Continuous_histo::ascii_characteristic_print(ostream &os , bool comment_flag) const
{
  os << STAT_label[STATL_SAMPLE_SIZE] << ": " << nb_element << "  "
     << "classes number: " << nb_class << endl;
  os << "Minimum: " << min_value << "   "
     << "Maximum: " << max_value <<endl;
  if (variance > 0.) {
    if (comment_flag) {
      os << "# ";
    }
    os << STAT_label[STATL_MEAN] << ": " << mean << "   "
       << STAT_label[STATL_VARIANCE] << ": " << variance << "   "
       << STAT_label[STATL_STANDARD_DEVIATION] << ": " << sqrt(variance) << endl;
     }
  else {
    os << "Variance must be positive"<<endl;
  }

  return os;
}


/*-------------------------------------------------------------- 
 *
 * Visualisation d'un objet Continuous_histo.
 *
 * arguments : stream
 *
 *--------------------------------------------------------------*/

ostream& Continuous_histo::print(ostream &os) const
{
  register int i;

  os << endl;
  ascii_characteristic_print(os);

  os << "frequencies ( " << nb_class << " classes, step: " <<((max_value - min_value)/nb_class)<<") : "<<endl;
  for (i = 0;i < nb_class;i++) {
    os << "[ "<< min_value + i * (max_value - min_value) / nb_class <<" ; "
       << min_value + (i+1) * (max_value - min_value) / nb_class << " [: " << frequency[i] << " "<<endl;
  }
  os << endl;

  return os;
}


/*--------------------------------------------------------------
 * 
 *Calcul de la quantite d'information d'un histogramme continu.
 * voir feuille imprimée sur l'entropie
 *
 *--------------------------------------------------------------*/

double Continuous_histo::information_computation() const // does not work
{
  register int i;
  double *pfrequency;
  double sum = 0.;
  double information = D_INF;

  if (nb_element > 0) {
    pfrequency = frequency;
    information = 0.;

    for (i = 0;i < nb_class;i++) {
      sum = sum + pfrequency[i];    
    }

    for (i = 0;i < nb_class;i++) {
      if (*pfrequency > 0) {
        information += *pfrequency * log((double)*pfrequency / sum);
      }
      pfrequency++;
    }
  }
  information += information + log((max_value - min_value / nb_class));
  return information;
}


/*--------------------------------------------------------------
 *
 * Calcul de la log-vraisemblance d'un histogramme continu.
 *
 * arguments : seuil de rejet
 *
 *--------------------------------------------------------------*/

double Continuous_histo::likelihood_computation(double icumul_threshold) const
{
  Continuous_distribution *cont_dist;
  cont_dist=new Continuous_distribution(NORMAL, (*this).mean, (*this).variance, icumul_threshold, (*this).nb_element);
  return cont_dist->likelihood_computation();
}


//--------------------------------------------------------------
// Calcul de la log-vraisemblance d'un histogramme continu.
//--------------------------------------------------------------
// double Continuous_histo::likelihood_computation() const
// {
//   register int i;
//   double *pfrequency, *pfrequency2;
//   double sum=0.;
//   double information = D_INF;
//   Continuous_histo *histo_tmp;
//   Continuous_histo *histo_tmp_freq;
//   histo_tmp = new Continuous_histo(*this, nb_element * (max_value - min_value)/nb_class);
//   histo_tmp_freq = new Continuous_histo(*this, 1.);

//   if (nb_element > 0) {
//     pfrequency = histo_tmp->frequency;
//     pfrequency2 = histo_tmp_freq->frequency;
//     information = 0.;

//     for (i = 0;i < nb_class;i++) {
//       sum = sum + pfrequency[i];    
//     }
//     cout<<"sum: "<< sum<<endl;
//     for (i = 0;i < nb_class;i++) {
//       if (*pfrequency > 0) {
//         information += log((double)*pfrequency * (double)*pfrequency2 /( sum*(max_value - min_value / nb_class)));
// 	  cout<<"information: "<<information<<endl;
// 	cout<<"pfrequency: "<<*pfrequency<<endl;
// 	cout<<"pfrequency2: "<<*pfrequency2<<endl;

//       }
//       pfrequency++;
//       pfrequency2++;
//     }
//   }
//   return information;
// }


/*--------------------------------------------------------------
 *
 * Calcul de la vraisemblance d'un histogramme pour une loi continue donnee.
 *
 * arguments : reference sur un objet Continuous_distribution
 *
 *--------------------------------------------------------------*/

double Continuous_histo::likelihood_computation(Continuous_distribution &cont_dist) const
{
  return cont_dist.likelihood_computation();
}


/*--------------------------------------------------------------
 *
 * Estimation d'une loi continue a partir d'un histogramme.
 *
 * arguments : reference sur un objet Continuous_distribution
 *
 *--------------------------------------------------------------*/

void Continuous_histo::distribution_estimation(Continuous_distribution *cont_dist) const
{
  if (nb_element > 0) {
    register int i;
    double *pfrequency;
    double *pdensity;
    double *pcumul, *qcumul;
    double tmp, step, sum = 0.;

    cont_dist->mean = mean;
    cont_dist->variance = variance;
    cont_dist->nb_value = nb_class + 1;
    cont_dist->min_value = min_value - (max_value - min_value) / (2*nb_class);
    cont_dist->max_value = max_value + (max_value - min_value) / (2*nb_class);

    pdensity = new double[nb_class+2];
    pcumul = new double[nb_class+2];
    qcumul = new double[nb_class+2];

    for (i = 0;i < 1;i++) {
      pdensity[0] = 0.;
      pcumul[0] = 0.;
      qcumul[0] = 0.;
    }
    
    qcumul[nb_class+1] = 1.;
    pdensity[nb_class+1] = 0.;
    pcumul[nb_class+1] = 1.;
    
    step = (max_value - min_value) / nb_class;

    pfrequency = frequency;

    for (i = 0;i < nb_class;i++) {
      sum = sum + pfrequency[i];
    }

    for (i = 1;i <= nb_class;i++) {
      tmp = (double)pfrequency[i-1] / (sum * step);
      pdensity[i] = tmp;
      qcumul[i] = qcumul[i-1] + pdensity[i] * step;
    }
    
    for (i = 1;i <= nb_class; i++) {
      pcumul[i] = (qcumul[i]+ qcumul[i-1]) / 2.;
    }

    cont_dist->density = pdensity;
    cont_dist->cumul = pcumul;
  }
}


/*---------------------------------------------------------
 *
 * Operateur d'egalite de la classe Continuous_histo
 *
 * arguments : reference sur un objet Continuous_histo
 *
 *---------------------------------------------------------*/

bool Continuous_histo::operator==(const Continuous_histo &cont_histo) const
{
  bool status = true;
  register int i;
  double *pfrequency , *cfrequency;


  if ((nb_class != cont_histo.nb_class) || (min_value != cont_histo.min_value) ||
      (nb_element != cont_histo.nb_element) || (max_value != cont_histo.max_value)) {
    status = false;
  }

  else {
    pfrequency = frequency;
    cfrequency = cont_histo.frequency;
    for (i = 0;i < nb_class;i++) {
      if (*pfrequency++ != *cfrequency++) {
        status = false;
        break;
      }
    }
  }

  return status;
}


/*----------------------------------------------------------------
 *
 * Translation d'un histogramme continue.
 *
 * arguments : reference sur un objet Continuous_histo, parammètre de decalage
 *
 * ATTENTION: dans le cas exponentiel, 
 * l'histogramme translaté n'a plus une distribution exponentielle
 *
 *----------------------------------------------------------------*/

void Continuous_histo::shift(const Continuous_histo &cont_histo , double shift_param)
{
  register int i;
  double *pfrequency , *cfrequency;
  
  // calcul des caracteristiques de l'histogramme
  nb_element = cont_histo.nb_element;
  nb_class = cont_histo.nb_class;
  min_value = cont_histo.min_value + shift_param;
  max_value = cont_histo.max_value + shift_param;

  // copie des frequences
  frequency = new double[nb_class];
  
  pfrequency = frequency;
  cfrequency = cont_histo.frequency;
  for (i = 0;i < nb_class;i++) {
    *pfrequency++ = *cfrequency++;
  }
  
  mean_computation();
  variance_computation();
}


/*--------------------------------------------------------------
 *
 * Ecriture d'un histogramme continu.
 *
 * arguments: stream, flag de commentaire, flag de densité et de 
 *            repartition
 *
 * --------------------------------------------------------------*/

ostream& Continuous_histo::ascii_print(ostream &os , int comment_flag , bool cumul_flag, bool density_flag) const
{
  register int i;
  int width[6];
  long old_adjust;
  double *cumul, *density;
  double step;

  step = (max_value - min_value) / nb_class;
  
  old_adjust = os.setf(ios::right , ios::adjustfield);
  
  // calcul des largeurs des colonnes
  width[0] = 3;
  width[1] = 18 + ASCII_SPACE;
  width[2] = 18 + ASCII_SPACE;
  width[3] = 15 + ASCII_SPACE;

  if (density_flag) {
    density = density_computation();
    width[4] =  27 + ASCII_SPACE;
  }
  if (cumul_flag) {
    cumul = cumul_computation();
    width[5] =  30 + ASCII_SPACE;
  }
  
  // ecriture 
  for (i = 0;i < nb_class;i++) {
    if (comment_flag == 1) {
      os << "# ";
    }
    os << setw(width[0]) << i;
    os << setw(width[1]) << min_value + i * step;
    os << setw(width[2]) << min_value + (i + 1) * step;
    os << setw(width[3]) << frequency[i];

    if (density_flag) {
      if (comment_flag == 0) {
        os << "  #";
      }
      os << setw(width[4]) << density[i+1];
    }
    if (cumul_flag) {
      if (comment_flag == 0) {
        os << "  #";
      }
      os << setw(width[5]) << cumul[i+1];
    }
    os << endl;
  }
  
  if (density_flag) {
    delete [] density;
  }
  if (cumul_flag) {
    delete [] cumul;
  }

  os.setf((FMTFLAGS)old_adjust , ios::adjustfield);

  return os;
}


/*--------------------------------------------------------------
 *
 * Ecriture d'un objet Continuous_histo.
 *
 * arguments: stream, flag de sauvegarde
 *
 * --------------------------------------------------------------*/

ostream& Continuous_histo::ascii_write(ostream &os , bool exhaustive , bool file_flag) const
{
  double information = information_computation();

  if (exhaustive && file_flag) {
    os << "# ";
  }
  os << STAT_label[STATL_HISTOGRAM] << " - ";
  ascii_characteristic_print(os , exhaustive && file_flag);

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
    os << "    | " << "class lower bound" << " | " << "class upper bound" << " |      " << STAT_label[STATL_HISTOGRAM] << " | " 
       << "density" << " " << STAT_label[STATL_HISTOGRAM] << " " << STAT_label[STATL_FUNCTION] <<" | " << STAT_label[STATL_CUMULATIVE]
       << " " << STAT_label[STATL_HISTOGRAM] << " " << STAT_label[STATL_FUNCTION] << endl;
    if (file_flag) {
      ascii_print(os , 1  , true, true);
    }
    else {
      ascii_print(os, 0 ,true, true);
    }
  }

  return os;
}


/*--------------------------------------------------------------
 *
 * Ecriture d'un objet Continuous_histo dans un fichier.
 *
 * arguments: adresse
 *
 * --------------------------------------------------------------*/

bool Continuous_histo::ascii_write( const char *path) const
{
  bool status;
  ofstream out_file(path);

  if (!out_file) {
    status = false;
  }
  else {
    status = true;
    ascii_write(out_file , true , true);
  }

  return status;
}


/*--------------------------------------------------------------
 *
 * Ecriture des caractéristiques d'un histogramme continu au format tableur
 *
 * arguments: stream
 *
 * --------------------------------------------------------------*/

ostream& Continuous_histo::spreadsheet_characteristic_print(ostream &os) const
{
  os << STAT_label[STATL_SAMPLE_SIZE] << "\t" << nb_element << "\t"
     << "classes number\t " << nb_class <<endl;
  os << "Minimum" << "\t" << min_value << "\t" 
     << "Maximum" << "\t" << max_value << "\t" <<endl;

  if (variance > 0.) {
    os << STAT_label[STATL_MEAN] << "\t" << mean << "\t\t"
       << STAT_label[STATL_VARIANCE] << "\t" << variance << "\t\t"
       << STAT_label[STATL_STANDARD_DEVIATION] << "\t" << sqrt(variance) << endl;
  }
  else {
    os << "Variance must be positive" <<endl;
  }

  return os;
}


/*--------------------------------------------------------------
 *
 * Ecriture d'un histogramme continu au format tableur.
 *
 * arguments: stream, flag de commentaire, flag de densité et de 
 *            repartition
 *
 * --------------------------------------------------------------*/

ostream& Continuous_histo::spreadsheet_print(ostream &os , bool cumul_flag, bool density_flag) const
{
  register int i;
  double *cumul;
  double *density;
  double step;

  step = (max_value - min_value) / nb_class;

  if (cumul_flag) {
    cumul = cumul_computation();
  }
  if (density_flag) {
    density = density_computation();
  }

  for (i = 0;i < nb_class;i++) {
    os << (min_value + i * step) << "\t" << (min_value + (i + 1) * step) << "\t" << frequency[i];
    if (density_flag) {
      os << "\t" << density[i+1];
    }
    if (cumul_flag) {
      os << "\t" << cumul[i+1];
    }
    os << endl;
  }

  if (density_flag) {
    delete [] density;
  }
  if (cumul_flag) {
    delete [] cumul;
  }

  return os;
}


/*-------------------------------------------------------------------
 *
 * Calcul de la fonction de densité déduite d'un histogramme continu
 *
 *-------------------------------------------------------------------*/

double* Continuous_histo::density_computation( ) const
{
  register int i;
  double *pdensity, *pfrequency;
  double tmp, step, sum = 0.;
  pdensity = new double[nb_class+2];

  pdensity[0] = 0.;
  pdensity[nb_class+1] = 0.;
  
  step = (max_value - min_value) / nb_class;
  
  pfrequency = frequency;

  for (i = 0;i < nb_class;i++) {
    sum = sum + pfrequency[i];
  }
  
  for (i = 1;i <= nb_class;i++) {
    tmp = (double)pfrequency[i-1] / (sum * step);
    pdensity[i] = tmp;
  }

  return pdensity;
}


/*-------------------------------------------------------
 *
 * Calcul de la fonction de répartition 
 * déduite d'un histogramme continu
 *
 *-------------------------------------------------------*/

double* Continuous_histo::cumul_computation( ) const
{
  register int i;
  double *pcumul, *qcumul, *pfrequency;
  double tmp, step, sum = 0.;

  pcumul = new double[nb_class+2];
  qcumul = new double[nb_class+2];
  
  pcumul[0] = 0.;
  qcumul[0] = 0.;
      
  qcumul[nb_class+1] = 1.;
  pcumul[nb_class+1] = 1.;
    
  step = (max_value - min_value) / nb_class;
  
  pfrequency = frequency;

  for (i = 0;i < nb_class;i++) {
    sum = sum + pfrequency[i];
  }

  for (i = 1;i <= nb_class;i++) {
    tmp = (double)pfrequency[i-1] / (sum * step);
    qcumul[i] = qcumul[i-1] + tmp * step;
  }
  
  for (i = 1;i<= nb_class; i++) {
    pcumul[i] = (qcumul[i]+ qcumul[i-1]) / 2.;
  }

  return pcumul;
}
