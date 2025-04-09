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
#include <time.h>
#include "stat_tool/stat_label.h"
#include "stat_tool/stat_tools.h"
#include "tool/config.h"
#include "continuous_distribution.h"

using namespace std;


/*--------------------------------------------------
 *
 * Constructeur de la classe Continuous_distribution
 *
 *--------------------------------------------------*/

Continuous_distribution::Continuous_distribution()
{
  type = 0;
  mean = 0.;
  variance = 1.;
  density = 0;
  cumul = 0;
  nb_value = 0;
  min_value = 0.;
  max_value = 0.;
  cumul_threshold = CUMUL_DEFAULT;
}


/*--------------------------------------------------
 *
 * Constructeur de la classe Continuous_distribution
 *
 * arguments : nombre de valeurs
 *
 *--------------------------------------------------*/

Continuous_distribution::Continuous_distribution (int inb_value)
{
  nb_value = inb_value;

  type = 0;
  mean = 0.;
  variance = 1.;
  min_value = 0.;
  max_value = 0.;
  cumul_threshold = CUMUL_DEFAULT;

  if (nb_value == 0) {
    density = 0;
    cumul = 0;
  }

  else {
    register int i;
    double *pdensity , *pcumul;

    density = new double[nb_value + 1];
    cumul = new double[nb_value + 1];

    pdensity = density;
    pcumul = cumul;

    for (i = 0;i <= nb_value;i++) {
      *pdensity++ = 0.;
      *pcumul++ = 0.;
    }

    pdensity = NULL;
    delete pdensity;

    pcumul = NULL;
    delete pcumul; 

  }
}


/*--------------------------------------------------
 *
 * Constructeur de la classe Continuous_distribution
 *
 * arguments : type de loi, moyenne, variance, 
 *             seuil de rejet, nombre de valeurs
 *
 *--------------------------------------------------*/

Continuous_distribution::Continuous_distribution (int itype, double imean, double ivariance, double icumul_threshold, int inb_value)
{
  register int i;
  double tmp;
  double step;

  type = itype;
  mean = imean;
  variance = ivariance;
  cumul_threshold = icumul_threshold;

  nb_value = inb_value;


  switch(type){
  case NORMAL : 
    if (variance < 0. ) {
      cout<<"ERROR: variance must be positive"<<endl;
      
      density = new double[0];
      cumul = new double[0];
    }
    else {
      qnorm();
      min_value = - (*this).max_value * sqrt(variance) + mean;
      max_value = (*this).max_value * sqrt(variance) + mean;
      step = (max_value - min_value) / nb_value;

      density = new double[nb_value+1];
      cumul = new double[nb_value+1];

      density[0] = dnorm(min_value, mean, variance);
      cumul[0] = pnorm( (min_value - mean) / sqrt(variance) );
      for (i = 1; i <= nb_value ; i++) {
	density[i] = dnorm((min_value + i * step), mean, variance);
	cumul[i] = pnorm(( (min_value + i * step) - mean) / sqrt(variance));
      }
    }
    
    break;
  
  case EXPONENTIAL :
    tmp = mean * mean;
    if ((mean <= 0.) || (variance != tmp)  ) {
      cout<<"ERROR: mean must be positive or variance must be the mean squared"<<endl;
      
      density = new double[0];
      cumul = new double[0];
    }
    else {
      qexp();
      min_value = 0.;
      max_value = (*this).max_value * mean ; 
      step = max_value / nb_value;
      
      density = new double[nb_value+2];
      cumul = new double[nb_value+2];

      density[0] = dexp(min_value, mean);
      cumul[0] = pexp( min_value,mean );
      for (i = 1; i <= nb_value ; i++) {
	density[i] = dexp(min_value + i * step, mean);
	cumul[i] = pexp( min_value + i * step,mean);
      }

    }
    break;
  }
}


/*--------------------------------------------------
 *
 * Constructeur de la classe Continuous_distribution
 *
 * arguments : borne inférieure, borne supérieure, 
 *             type de loi, seuil de rejet, nombre de valeurs
 *
 *--------------------------------------------------*/

Continuous_distribution::Continuous_distribution(double imin_value, double imax_value, int itype, double icumul_threshold, int inb_value)
{
  register int i;
  double step;

  type = itype;
  max_value = imax_value;
  min_value = imin_value;
  cumul_threshold = icumul_threshold; 

  nb_value=inb_value;

  switch(type) {
  case NORMAL:
    mean = (max_value + min_value) / 2; 
    variance = pow(( max_value - mean) / qnorm(cumul_threshold),2);
    step = (max_value - min_value) / nb_value;

    density = new double[nb_value+1];
    cumul = new double[nb_value+1];
    
    density[0] = dnorm(min_value, mean, variance);
    cumul[0] = pnorm( (min_value - mean) / sqrt(variance) );
    for (i = 1; i <= nb_value ; i++) {
      density[i] = dnorm((min_value + i * step), mean, variance);
      cumul[i] = pnorm(( (min_value + i * step) - mean) / sqrt(variance));
    }
    
    break;
    
  case EXPONENTIAL :
    mean = max_value / qexp(cumul_threshold);
    variance = mean * mean;

    step = (max_value - min_value) / nb_value;

    density = new double[nb_value+1];
    cumul = new double[nb_value+1];
    
    density[0] = dexp(min_value, mean);
    cumul[0] = pexp( min_value,mean );
    for (i = 1; i <= nb_value ; i++) {
      density[i] = dexp(min_value + i * step, mean);
      cumul[i] = pexp( min_value + i * step,mean);
    }

    break;
  }
}


/*--------------------------------------------------
 *
 * Constructeur de la classe Continuous_distribution
 *
 * arguments : référence sur un objet de type Continuous_histo
 *
 *--------------------------------------------------*/

Continuous_distribution::Continuous_distribution(const Continuous_histo &cont_histo)
{
  cont_histo.distribution_estimation(this);
}


/*--------------------------------------------------
 *
 * Quantile d'une loi exponentielle avec un taux lambda=1
 *
 * arguments : nombre
 *
 *--------------------------------------------------*/

double qexp (double p)
{
  double tmp;

  tmp = - log(1-p);
  
  return tmp;
}


void Continuous_distribution::qexp()
{
  max_value = ::qexp(1. - cumul_threshold);
}


double Continuous_distribution::qexp(double p)
{
  double val;

  val = ::qexp(1. - p);

  return val;
}


/*--------------------------------------------------
 *
 * Quantile d'une loi normale centrée réduite
 *
 * arguments : nombre
 *
 *--------------------------------------------------*/

double qnorm (double p)
{
  double split=0.42;
  double a0=  2.50662823884;
  double a1=-18.61500062529;
  double a2= 41.39119773534;
  double a3=-25.44106049637;
  double b1= -8.47351093090;
  double b2= 23.08336743743;
  double b3=-21.06224101826;
  double b4=  3.13082909833;
  double c0= -2.78718931138;
  double c1= -2.29796479134;
  double c2=  4.85014127135;
  double c3=  2.32121276858;
  double d1=  3.54388924762;
  double d2=  1.63706781897;
  double q=p-0.5;
  double r, ppnd;

  if(fabs(q) <= split)
    {
      r=q*q;
      ppnd=q*(((a3*r+a2)*r+a1)*r+a0)/((((b4*r+b3)*r+b2)*r+b1)*r+1);
    }
  else
    {
      r=p;
      if(q>0) r=1-p;
      if(r>0) 
	{
	  r=sqrt(-log(r));
	  ppnd=(((c3*r+c2)*r+c1)*r+c0)/((d2*r+d1)*r+1);
	  if(q<0) ppnd=-ppnd;
	}
      else 
	{
	  ppnd=0;
	}
    }
  return ppnd;
}


void Continuous_distribution::qnorm()
{

  if (cumul_threshold < 0.5) {
    max_value = ::qnorm(1. - cumul_threshold / 2);
  }
  else {
    max_value = ::qnorm((1. - cumul_threshold) / 2 + 0.5);
  }

}


double Continuous_distribution::qnorm(double p)
{
  double val;
  if (p < 0.5) {
    val = ::qnorm(1. - p / 2);
  }
  else {
    val = ::qnorm((1. - p) / 2 + 0.5);
  }
  return val ;
}


/*--------------------------------------------------
 *
 * Fonction de répartition d'une loi normale centrée réduite
 *
 * arguments : nombre
 *
 *--------------------------------------------------*/

double Continuous_distribution::pnorm(double x)
{
  double  c = .70710678118654752440;

  if(x >= 0)
    return 0.5*(1+erf(c*x));
  else
    return 0.5*erfc(-c*x);
}


void Continuous_distribution::pnorm()
{
  cumul_threshold = pnorm(min_value) * 2;
}


/*--------------------------------------------------
 *
 * Fonction de répartition d'une loi exponentielle
 *
 * arguments : nombre, moyenne
 *
 *--------------------------------------------------*/

double Continuous_distribution::pexp(double x, double imean)
{
  if(x >= 0)
    return 1. - exp( - x / imean);
  else
    return 0.;
}


void Continuous_distribution::pexp()
{
  cumul_threshold = 1 - pexp(max_value,mean);
}


/*--------------------------------------------------
 *
 * Densité d'une loi normale
 *
 * arguments : nombre, moyenne, variance
 *
 *--------------------------------------------------*/

double Continuous_distribution::dnorm(double d, double imean, double ivariance)
{
  return ( 1. / ( sqrt ( 2. * S_PI * ivariance))) * exp( - ((d-imean)*(d-imean)) / (2 * ivariance)); 
}


/*--------------------------------------------------
 *
 * Densité d'une loi exponentielle avec un taux lambda = 1
 *
 * arguments : nombre, moyenne
 *
 *--------------------------------------------------*/

double Continuous_distribution::dexp(double d, double imean)
{
  return ( 1. / imean) * exp( - d / imean); 
}


/*--------------------------------------------------
 *
 * Approximation de la fonction de repartition 
 * d'une loi normale centree reduite(pnorm)
 *
 * arguments : nombre, moyenne
 *
 *--------------------------------------------------*/

double standard_normal_cumul_density(double value)
{
  register int i;
  double critical_probability, term, var[7];

  var[1] = fabs(value);
  for(i = 2;i <= 6; i++) {
    var[i] = var[i-1] * var[1];
  }

  term = 1. + 0.0498673470 * var[1] + 0.0211410061 * var[2] + 0.0032776263 * var[3] + 
    0.0000380036 * var[4] + 0.0000488906 * var[5] + 0.0000053830 * var[6];

  var[0] = 1.;
  for(i = 0;i < 16; i++) {
    var[0] *= term;
  }
  critical_probability = 0.5 / var[0];

  if(value > 0.){
    critical_probability = 1. - critical_probability;
  }

  return critical_probability;
}


/*--------------------------------------------------
 *
 * Approximation d'une valeur pour une loi normale centree reduite(qnorm)
 *
 * arguments : nombre, moyenne
 *
 *--------------------------------------------------*/

double standard_normal_value_threshold(double p)
{
  register int i;
  double value_threshold, var[3];

  if ( p <= 0.5 ) {
    var[0] = sqrt ( log ( 1. / (p * p) ) );
  }  
  else {
    var[0] = sqrt ( log ( 1. / ((1-p) *(1- p)) ) );
  }

  for (i = 1;i <= 2;i++) {
    var[i] = var[i-1] * var[0];
  }

  value_threshold = var[0] - ( 2.515517 + 0.802853 * var[0] + 0.010328 * var[1]) / 
    (1. + 1.432788 * var[0] + 0.189269 * var[1] + 0.001308 * var[2]);

  if (p < 0.5) {
    value_threshold = - value_threshold;
  }

  return value_threshold;
}


/*---------------------------------------------------------------
 *
 * Destructeur de la classe Continuous_observation
 *
 *---------------------------------------------------------------*/ 

Continuous_distribution::~Continuous_distribution()
{
  if (density) {
    delete[] density;
  }
  if (cumul) {
    delete[] cumul;
  }
}


/*---------------------------------------------------------------
 *
 * Copie d'un objet Continuous_distribution
 *
 * arguments :  reference sur un objet de type Continuous_distribution
 *
 *---------------------------------------------------------------*/

void Continuous_distribution::copy(const Continuous_distribution &cont_dist)
{
  register int i;
  double *pdensity, *pcumul, *cdensity, *ccumul;

  type = cont_dist.type;
  mean = cont_dist.mean;
  variance = cont_dist.variance;
  nb_value = cont_dist.nb_value;
  min_value = cont_dist.min_value;
  max_value = cont_dist.max_value;
  cumul_threshold = cont_dist.cumul_threshold;

  density = new double[nb_value+1];
  cumul = new double[nb_value+1];

  pdensity = density;
  cdensity = cont_dist.density;
  pcumul = cumul;
  ccumul = cont_dist.cumul;

  for (i = 0; i <= nb_value ; i++) {
    *pdensity++ = *cdensity++;
    *pcumul++ = *ccumul++;  
  }

  pdensity = NULL;
  cdensity = NULL;
  pcumul = NULL;
  ccumul = NULL;

  delete pdensity;
  delete cdensity;
  delete pcumul;
  delete ccumul;
}


/*---------------------------------------------------------------
 *
 * Opérateur d'assignementde la classe Continuous_distribution
 *
 * arguments :  reference sur un objet de type Continuous_distribution
 *
 *---------------------------------------------------------------*/

Continuous_distribution& Continuous_distribution::operator=(const Continuous_distribution &cont_dist)
{
  if (&cont_dist != this) {
    delete [] density;
    delete [] cumul;

    copy(cont_dist);
  }

  return *this;
}


/*---------------------------------------------------------------
 *
 * Constructeur par copie de la classe Continuous_distribution
 *
 * arguments :  reference sur un objet de type Continuous_distribution
 *
 *---------------------------------------------------------------*/

Continuous_distribution::Continuous_distribution(const Continuous_distribution &cont_dist)
{
  copy(cont_dist);
}


/*------------------------------------------------------------
 *
 * Calcul de la borne inférieure de la distribution
 *
 *------------------------------------------------------------*/

void Continuous_distribution::min_value_computation()
{
  switch(type){
  case NORMAL : 
    qnorm();
    min_value = - (*this).max_value * sqrt(variance) + mean;
    break;
  
  case EXPONENTIAL :
    qexp();
    min_value = 0.;
    break;
  }
}


/*------------------------------------------------------------
 *
 * Calcul de la borne supérieure de la distribution
 *
 *------------------------------------------------------------*/

void Continuous_distribution::max_value_computation()
{
  switch(type){
  case NORMAL : 
    qnorm();
    max_value = (*this).max_value * sqrt(variance) + mean;
    break;
  
  case EXPONENTIAL :
    qexp();
    max_value = (*this).max_value * mean ; 
    break;
  }
}


/*------------------------------------------------------------
 *
 * Calcul de la valeur de densité maximale
 *
 *------------------------------------------------------------*/

double Continuous_distribution::max_computation()
{
  register int i;
  double max;
  double *pdensity;

  pdensity = density;
  max = 0.;
  for (i = 0;i <= nb_value;i++) {
    if (*pdensity > max) {
      max = *pdensity;
    }
    pdensity++;
  }

  pdensity = NULL;
  delete pdensity;

  return max;
}


/*--------------------------------------------------------------
 *
 * Calcul du nombre de valeurs
 *
 *--------------------------------------------------------------*/

void Continuous_distribution::nb_value_computation()//does not work
{
  register int i = 0;
  int sum = 0;

  nb_value = (*this).nb_value * sqrt((*this).variance)/2.;
}


/*--------------------------------------------------------------
 *
 * Calcul de la moyenne d'une loi continue.(Methode des trapèzes)
 *
 *--------------------------------------------------------------*/

void Continuous_distribution::mean_computation()
{
  register int i;
  double step;

  step = (max_value - min_value) / nb_value;
  
  variance = min_value * density[0] + (min_value + (nb_value * step)) * density[nb_value]; 

  for(i = 1;i < nb_value;i++) {
    mean = mean + density[i] * (min_value + i*step) ;
  }
  mean *= step;
}


/*--------------------------------------------------------------
 *
 * Calcul de la variance d'une loi continue.(Methode des trapèzes)
 *
 *--------------------------------------------------------------*/

void Continuous_distribution::variance_computation()
{
  register int i;
  double step;

  step = (max_value - min_value) / nb_value;
  
  variance = pow(min_value - mean,2) * density[0] + pow(min_value + (nb_value * step) - mean,2) * density[nb_value]; 

  for(i = 1;i < nb_value;i++) {
    variance = variance + density[i] * pow(min_value + i*step - mean,2) ;
  }
  variance *= step;
}


/*--------------------------------------------------------------
 *
 * Operateur d'egalite de la classe Continuous_distribution.
 *
 * arguments : reference sur un objet de la classe Continuous_distribution
 *
 *--------------------------------------------------------------*/

bool Continuous_distribution::operator==(const Continuous_distribution &cont_dist) const
{
  bool status = true;
  register int i;
  double *pdensity , *ddensity;


  if ((min_value != cont_dist.min_value) || (max_value != cont_dist.max_value) || (nb_value != cont_dist.nb_value)) {
    status = false;
  }

  else {
    pdensity = density;
    ddensity = cont_dist.density ;
    for (i = 0;i <= nb_value;i++) {
      if (*pdensity++ != *ddensity++) {
        status = false;
        break;
      }
    }

    pdensity = NULL;
    delete pdensity;

    ddensity = NULL;
    delete ddensity;
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture des caracteristiques d'une loi continue.
 *
 *  arguments : stream, flag ecriture des parametres de forme, flag commentaire.
 *
 *--------------------------------------------------------------*/

ostream& Continuous_distribution::ascii_characteristic_print(ostream &os , bool comment_flag) const
{
  if ((mean != D_DEFAULT) && (variance != D_DEFAULT)) {
    if (comment_flag) {
      os << "# ";
    }
    os << STAT_label[STATL_MEAN] << ": " << mean << "   "
       << STAT_label[STATL_VARIANCE] << ": " << variance << "   "
       << STAT_label[STATL_STANDARD_DEVIATION] << ": " << sqrt(variance) << endl;
  }
  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture des caracteristiques d'une loi continue au format tableur.
 *
 *  arguments : stream, flag ecriture des parametres de forme.
 *
 *--------------------------------------------------------------*/

ostream& Continuous_distribution::spreadsheet_characteristic_print(ostream &os) const
{
  if ((mean != D_DEFAULT) && (variance != D_DEFAULT)) {
    os << STAT_label[STATL_MEAN] << "\t" << mean << "\t"
       << STAT_label[STATL_VARIANCE] << "\t" << variance << "\t"
       << STAT_label[STATL_STANDARD_DEVIATION] << "\t" << sqrt(variance) << endl;
  }
  
  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'une loi au format Gnuplot.(pas testée)
 *
 *  arguments : path, pointeur sur la fonction de concentration,
 *              facteur d'echelle.
 *
 *--------------------------------------------------------------*/

bool Continuous_distribution::plot_print(const char *path , double *concentration ,
					 double scale) const
{
  bool status = false;
  register int i;
  double step;
  ofstream out_file(path);
  
  step = (max_value - min_value) / nb_value;
  if (out_file) {
    status = true;
    
    for (i = 0;i <= nb_value;i++) {
      out_file << (min_value+i*step) << " " << density[i] * scale;
      if (variance > 0.) {
        out_file << " " << cumul[i] << " ";
      }
      out_file << endl;
    }
  }

  return status;
}


//--------------------------------------------------------------
//  Ecriture d'une famille de fonctions de repartition en vue du matching
//  au format Gnuplot
//--------------------------------------------------------------
// bool cumul_matching_plot_print(const char *path , int nb_cumul , double *min_value , double *max_value,
//                                double *step , double **cumul)

// {
//   bool status = false;
//   register int i , j;
//   double plot_min_value, plot_max_value;
//   int plot_nb_value;
//   ofstream out_file(path);
//   int *nb_value;
//   nb_value=new int[nb_cumul];

//   for (i = 0;i < nb_cumul; i++) {
//     nb_value[i]=int(round((max_value[i] - min_value[i]) / step[i]));
//   }
 
//   if (out_file) {
//     status = true;

//     plot_min_value = min_value[0];
//     plot_max_value = max_value[0];
//     plot_nb_value = nb_value[0];
//     for (i = 1;i < nb_cumul;i++) {
//       if (min_value[i] < plot_min_value) {
//         plot_min_value = min_value[i];
//       }
//       if (max_value[i] > plot_max_value) {
// 	plot_max_value = max_value[i];
//       }
//       //if (nb_value[i] > plot_nb_value) {
//       //  plot_nb_value = nb_value[i];
//       //}
//     }

//     for (i = 0;i < nb_cumul;i++) {
//       out_file << 0 << " ";
//     }
//     out_file << endl;

//     for (i = 0;i < plot_nb_value;i++) {
//       for (j = 0;j < nb_cumul;j++) {
//         if (i < nb_value[j]) {
//           out_file << cumul[j][i] << " ";
//         }
//         else {
//           out_file << cumul[j][nb_value[j] - 1] << " ";
//         }
//       }
//       out_file << endl;
//     }
//   }

//   return status;
// }


/*--------------------------------------------------------------*
 *
 * Visualisation d'une loi.
 *
 * argument : stream.
 *
 *--------------------------------------------------------------*/

ostream& Continuous_distribution::print(ostream &os) const
{
  register int i;
  double step;
  double tmp;

  step = (max_value - min_value) / nb_value;

  ascii_characteristic_print(os);

  os << "minimum : " << min_value;
  os << "   maximum : " << max_value;
 
  os << endl;
  if(density){
    os << "density function (nb_value: " << nb_value << ", step: " << step <<") : ";
    for (i = 0;i <= nb_value;i++) {
      os << density[i] << " ";
    }
    os << endl;
  }
  if (cumul) {
    os << "cumulative distribution function (nb_value: " << nb_value << ", step: " << step <<") : ";
    for (i = 0;i <= nb_value;i++) {
      os << cumul[i] << " ";
    }
    os << endl;
  }
  return os;
}


/*--------------------------------------------------------------*
 *
 * Visualisation d'une distribution.
 *
 * arguments : stream, reference sur un objet Continuous_distribution.
 *
 *--------------------------------------------------------------*/

ostream& operator<<(ostream &os , const Continuous_distribution &cont_dist)
{
  os.precision(5);

  os << endl;
  cont_dist.print(os);

  os.precision(6);

  return os;
}


/*--------------------------------------------------------------*
 *
 * Calcul de la fonction de répartition d'une loi continue
 *
 * arguments : nombre de valeurs, seuil de rejet, borne inférieure,
 *              borne supérieure, densité, fonction de répartition.
 *
 *--------------------------------------------------------------*/

void cumul_computation(int inb_value , double icumul_threshold, double imin_value, double imax_value, const double *pdensity , double *pcumul)
{
  register int i;
  double step;

  step = (imax_value - imin_value) / inb_value;

  pcumul[0] = icumul_threshold/2.;

  for (i = 1;i < inb_value;i++) {
    pcumul[i] = pcumul[i-1] + pdensity[i] * step;
  }

  pcumul[inb_value] = 1-icumul_threshold/2.;
}


/*--------------------------------------------------------------
 *
 * Calcul de la fonction de repartition d'une loi continue.
 *
 *--------------------------------------------------------------*/

void Continuous_distribution::cumul_cont_computation()
{
  double *pcumul;

  delete [] cumul;
  cumul = new double[nb_value+1];
  pcumul = cumul;

  ::cumul_computation(nb_value , cumul_threshold, min_value, max_value,  density , pcumul);
}


/*---------------------------------------------------------------
 *
 * Calcul de la fonction de densite d'une loi continue
 *
 *---------------------------------------------------------------*/

void Continuous_distribution::density_cont_computation()
{
  register int i;
  double step;
  double *pdensity;
  step = (max_value - min_value) / nb_value;
  
  delete [] density;
  density = new double[nb_value+1];

  pdensity = density;
  switch (type){
  case NORMAL:
    for ( i = 0; i<= nb_value; i++){
      *pdensity++ = dnorm(min_value+i*step, mean, variance);
    }
    break;
  case EXPONENTIAL:
    for ( i = 0; i<= nb_value; i++){
      *pdensity++ = dexp(min_value+i*step, mean);
    }
    break;
  }

  pdensity = NULL;
  delete pdensity;
}


/*--------------------------------------------------------------
 *
 * Calcul des logarithmes des densités
 *
 * arguments : nombre de valeurs, densités, logarithmes des densités
 *
 *--------------------------------------------------------------*/

void log_computation(int inb_value, const double *pdensity, double *plog)
{
  register int i;

  for (i = 0;i <= inb_value;i++) {
    if (*pdensity > 0.) {
      *plog = log(*pdensity);
    }
    else {
      *plog = D_INF;
    }

    plog++;
    pdensity++;
  }
}
 

/*---------------------------------------------------------------
 *
 * Calcul des logarithmes des densités
 *
 *---------------------------------------------------------------*/

void Continuous_distribution::log_computation()
{
  if(cumul){
    ::log_computation(nb_value, density, cumul);
  }
}


/*---------------------------------------------------------------
 *
 * Tirage aléatoire selon une loi normale
 *
 * arguments : moyenne, variance
 *
 *---------------------------------------------------------------*/

double random_normal(double imean, double ivariance) 
{
  return imean + sqrt(ivariance) * sqrt( -2*log(((float)rand())/RAND_MAX))* cos(2*S_PI*((float)rand())/RAND_MAX);
}

double Continuous_distribution::random_normal()
{
  return ::random_normal(mean,variance);
}


/*---------------------------------------------------------------
 *
 * Tirage aléatoire selon une loi exponentielle
 *
 * arguments : moyenne,
 *
 *---------------------------------------------------------------*/

double random_expo(double imean)
{
  return - imean * log(((float)rand())/RAND_MAX);
}

double Continuous_distribution::random_expo()
{
  return ::random_expo(mean);
}


/*----------------------------------------------------------------
 *
 * Calcul de la log vraisemblance de données continues
 *
 *----------------------------------------------------------------*/

double Continuous_distribution::likelihood_computation()
{
  register int i;
  double tmp = 0.;
  for (i=0;i<= nb_value; i++) {
    tmp += log(density[i]);
  }
  
  return tmp;
}


/*-----------------------------------------------------------------
 *
 * Calcul de la log vraisemblance des données continues 
 * à partir de l'histogramme continu
 *
 * arguments : reference sur un objet Continuous_histo, seuil de rejet
 *
 *-----------------------------------------------------------------*/

double Continuous_distribution::likelihood_computation(const Continuous_histo &cont_histo,
						       double icumul_threshold) const
{ 
  return cont_histo.likelihood_computation(icumul_threshold); 
}


/*--------------------------------------------------------------*
 *
 * Ecriture d'une famille de lois et d'un histogramme.
 * en supposant que la loi et l'histogramme 
 * ont les mêmes bornes (min_value et max_value)
 *
 *  arguments : stream, nombre de lois, pointeurs sur les lois,
 *              facteurs d'echelle, flag commentaire,
 *              flag sur l'ecriture de la fonction de repartition,
 *              pointeur sur un histogramme.
 *
 *--------------------------------------------------------------*/

ostream& Continuous_distribution::ascii_print(ostream &os , bool comment_flag , bool density_flag, bool cumul_flag ,
					      bool nb_class_flag , const Continuous_histo *cont_histo) const
{
  register int i;
  int ascii_nb_value = 0; 
  int width[6];
  long old_adjust;
  double *cont_histo_cumul , *cont_histo_density, *pcumul;
  double step;

  if (cont_histo) {
    ascii_nb_value = cont_histo->nb_class;
    step = (cont_histo->max_value - cont_histo->min_value) / cont_histo->nb_class;
  }
  
  old_adjust = os.setf(ios::right , ios::adjustfield);

  // calcul du facteur d'echelle, de la fonction de repartition deduite de
  // l'histogramme et du nombre de valeurs
  
  if (cont_histo) {
    if (density_flag) {
      cont_histo_density = cont_histo->density_computation();
    }     
    if (cumul_flag) {
      cont_histo_cumul = cont_histo->cumul_computation();
    }
    
    if ((nb_class_flag) && (cont_histo->nb_class >= 2)) {
      pcumul = cumul + cont_histo->nb_class - 1;
      while (*--pcumul > 1. - ASCII_ROUNDNESS) {
	ascii_nb_value--;
      }
      
      if ((cont_histo) && (cont_histo->nb_class > ascii_nb_value)) {
	ascii_nb_value = cont_histo->nb_class;
      }
    }
  }
 
  // calcul des largeurs des colonnes
  
  width[0] = 3;
  if (cont_histo) {
    width[1] = 13 + ASCII_SPACE;
  }
  width[2] = 27 + ASCII_SPACE;
  if (density_flag) {
    width[3] = 18 + ASCII_SPACE;
  }
  width[4] = 30 + ASCII_SPACE;
  if (cumul_flag) {
    width[5] = 20 + ASCII_SPACE;
  }
  
  // ecriture des probabilites de chaque valeur
    for (i = 0;i < ascii_nb_value;i++) {
    if (comment_flag) {
      os << "# ";
    }
    os << setw(width[0]) << i;
    
    if (cont_histo) {
      if (i < cont_histo->nb_class) {
        os << setw(width[1]) << cont_histo->frequency[i];
      }
      else {
        os << setw(width[1]) << " ";
      }
    }
    
    if (density_flag) {
      if (cont_histo) {
        if (i < cont_histo->nb_class) {
          os << setw(width[2]) << cont_histo_density[i+1];
        }
        else {
          os << setw(width[2]) << " ";
        }
      }

      os << setw(width[3]) << (density[i * nb_value/cont_histo->nb_class ] + density[(i+1) * nb_value/cont_histo->nb_class]) / 2.;
    }
    
    if (cumul_flag) {
      if (cont_histo) {
        if (i < cont_histo->nb_class) {
          os << setw(width[4]) << cont_histo_cumul[i+1];
        }
        else {
          os << setw(width[4]) << " ";
        }
      }
      
      os << setw(width[5]) << (cumul[i * nb_value/cont_histo->nb_class ] + cumul[(i+1) * nb_value/cont_histo->nb_class]) / 2.;
    }
    os << endl;
  }
  
  if ((cont_histo) && (density_flag)) {
    delete [] cont_histo_density;
  } 

  if ((cont_histo) && (cumul_flag)) {
    delete [] cont_histo_cumul;
  }
  
  pcumul = NULL;
  delete pcumul;


  os.setf((FMTFLAGS)old_adjust , ios::adjustfield);
  
  return os;
}


/*--------------------------------------------------------------*
 *
 * Ecriture d'un objet Continuous_distribution
 *
 *  arguments : stream, flag commentaire,
 *              flag sur l'ecriture de la fonction de repartition,
 *              pointeur sur un histogramme.
 *
 *--------------------------------------------------------------*/

ostream& Continuous_distribution::ascii_write(ostream &os , bool exhaustive , bool file_flag, 
					      const Continuous_histo *cont_histo, double icumul_threshold) const
{
  double likelihood = (*this).likelihood_computation(*cont_histo, icumul_threshold);

  if (exhaustive && file_flag) {
    os << "# ";
  }
  os << STAT_label[STATL_DISTRIBUTION] << " - ";
  ascii_characteristic_print(os , exhaustive && file_flag);

  if (exhaustive && file_flag) {
    os << "# ";
  }
  os << STAT_label[STATL_LIKELIHOOD] << ": " << likelihood << " ("
     << likelihood / cont_histo->nb_element << ")" << endl;

  if (exhaustive) {
    os << "\n";
    if (file_flag) {
      os << "# ";
    }
    os << "    |    " << STAT_label[STATL_HISTOGRAM] << " | " 
       << "density" << " " << STAT_label[STATL_HISTOGRAM] << " " << STAT_label[STATL_FUNCTION] << " | "
       << "density" << "  " << STAT_label[STATL_FUNCTION] << " | " 
       << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_HISTOGRAM] << " " << STAT_label[STATL_FUNCTION] <<" | "
       << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_FUNCTION] << endl;
    if (file_flag) {
      ascii_print(os , true , true, true, true, cont_histo);
    }
    else {
      ascii_print(os, false ,true, true, true, cont_histo);
    }
  }
  return os;
}


/*--------------------------------------------------------------*
 *
 * Ecriture d'une famille de lois et d'un histogramme.
 * en supposant que la loi et l'histogramme 
 * ont les mêmes bornes (min_value et max_value) au format tableur
 *
 * arguments : stream, nombre de lois, pointeurs sur les lois,
 *              facteurs d'echelle, flag commentaire,
 *              flag sur l'ecriture de la fonction de repartition,
 *              pointeur sur un histogramme.
 *
 *--------------------------------------------------------------*/
ostream& Continuous_distribution::spreadsheet_print(ostream &os , bool cumul_flag , bool density_flag ,
						    bool nb_class_flag , const Continuous_histo *cont_histo) const
{
  register int i;
  int spreadsheet_nb_class = 0;
  double *pcumul , *cont_histo_cumul , *cont_histo_density;
  double step;

  if (cont_histo) {
    step = (cont_histo->max_value - cont_histo->min_value) / cont_histo->nb_class;
    spreadsheet_nb_class = cont_histo->nb_class;
  }

  // calcul de la fonction de repartition deduite de
  // l'histogramme, de la densité déduite de l'histogramme et du nombre de valeurs
  if (cont_histo) {
    if (density_flag) {
      cont_histo_density = cont_histo->density_computation();
    }
    if (cumul_flag) {
      cont_histo_cumul = cont_histo->cumul_computation();
    }
    
    if ((nb_class_flag) && (cont_histo->nb_class >= 2)) {
      pcumul = cumul + cont_histo->nb_class - 1;
      while (*--pcumul > 1. - SPREADSHEET_ROUNDNESS) {
	spreadsheet_nb_class--;
      }
      
      if ((cont_histo) && (cont_histo->nb_class > spreadsheet_nb_class)) {
	spreadsheet_nb_class = cont_histo->nb_class;
      }
    }
  }  

  // ecriture des effectifs, des densités, des fonctions de repartition
  for (i = 0;i < spreadsheet_nb_class;i++) {
    os << i;
    
    if (cont_histo) {
      os << "\t";
      if (i < cont_histo->nb_class) {
	os <<cont_histo->frequency[i];
      }
    }
    
    if (density_flag) {
      if (cont_histo) {
	os << "\t";
	if (i < cont_histo->nb_class) {
	  os << cont_histo_density[i+1];
	}
      }
	
      os << "\t" << (density[i * nb_value/cont_histo->nb_class ] + density[(i+1) * nb_value/cont_histo->nb_class]) / 2.; 
    }
    if (cumul_flag) {
      if (cont_histo) {
	os << "\t";
        if (i < cont_histo->nb_class) {
          os << cont_histo_cumul[i+1];
        }
      }
      
      os << "\t" << (cumul[i * nb_value/cont_histo->nb_class ] + cumul[(i+1) * nb_value/cont_histo->nb_class]) / 2.;
    }
    
    os << endl;
  }
  
  if (cont_histo) {
    if (density_flag) {
      delete [] cont_histo_density;
    }
    if (cumul_flag) {
      delete [] cont_histo_cumul;
    }
  }


  pcumul = NULL;
  delete pcumul;
  
  return os;
}


/*--------------------------------------------------------------*
 *
 * Ecriture d'un objet Continuous_distribution au format tableur
 *
 * arguments : stream, flag commentaire,
 *              flag sur l'ecriture de la fonction de repartition,
 *              pointeur sur un histogramme.
 *
 *--------------------------------------------------------------*/

ostream& Continuous_distribution::spreadsheet_write(ostream &os , bool exhaustive , bool file_flag, 
						    const Continuous_histo *cont_histo, double icumul_threshold) const
{
  double likelihood = likelihood_computation(*cont_histo, icumul_threshold);

  if (exhaustive && file_flag) {
    os << "# ";
  }
  os << STAT_label[STATL_DISTRIBUTION] << " - ";
  ascii_characteristic_print(os , exhaustive && file_flag);

  if (exhaustive && file_flag) {
    os << "# ";
  }
  os << STAT_label[STATL_LIKELIHOOD] << ": " << likelihood << " ("
     << likelihood / cont_histo->nb_element << ")" << endl;

  if (exhaustive) {
    os << "\n";
    if (file_flag) {
      os << "# ";
    }
    os << "    |    " << STAT_label[STATL_HISTOGRAM] << " | " 
       << "density" << " " << STAT_label[STATL_HISTOGRAM] << " " << STAT_label[STATL_FUNCTION] << " | "
       << "density" << "  " << STAT_label[STATL_FUNCTION] << " | " 
       << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_HISTOGRAM] << " " << STAT_label[STATL_FUNCTION] <<" | "
       << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_FUNCTION] << endl;
    if (file_flag) {
      spreadsheet_print(os , true , true, true,  cont_histo);
    }
    else {
      spreadsheet_print(os, false ,true, true,  cont_histo);
    }
  }

  return os;
}
