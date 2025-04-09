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
#include <sstream>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_eigen.h>
#include "tool/rw_tokenizer.h"
#include "tool/rw_cstring.h"
#include "tool/rw_locale.h"
#include "stat_tool/stat_tools.h"
#include "stat_tool/stat_label.h"
#include "continuous_parametric.h"

using namespace std;


extern void convert_matrix_double (double **essai, int nb_row, int nb_col, gsl_matrix* convert_essai);
extern void convert_vector_double (double *essai, int length, gsl_vector* essai_convert);


/*-----------------------------------------------------
 *
 * Constructeur de la classe Continuous_parametric
 *
 *-----------------------------------------------------*/

Continuous_parametric::Continuous_parametric() 
  : Continuous_distribution()
{
  ident = NO_RANDOM;
  nb_covariable = 0;
  nb_random = 0;
  residual_variance = 1.;
  random_variance = 0;
  regression = 0;
}


/*-----------------------------------------------------
 *
 * Constructeur de la classe Continuous_parametric
 *
 * arguments : nombre de valeur, nombre de covariable, nombre 
 *             d'effets aléatoires, variance résiduelle, écart-type
 *             des effets aléatoires, paramètres de régression,
 *             LM/LMM    
 *
 *-----------------------------------------------------*/

Continuous_parametric::Continuous_parametric(int inb_value, int inb_covariable, int inb_random, 
					     double iresidual_variance, double *irandom_variance, 
					     double *iregression, int iident)
  :Continuous_distribution(inb_value)
{ 
  register int i;
  double *pregression, *prandom;
  
  ident = iident;
  nb_covariable = inb_covariable;
  nb_random = inb_random;
  residual_variance = iresidual_variance;

  if (nb_random > 0){
    random_variance = new double[nb_random];
    prandom = random_variance;
    for (i = 0; i < nb_random ; i++) {
      *prandom++ = *irandom_variance++;
    }
  }
  else {
    random_variance = 0;
  }

  regression = new double[nb_covariable];

  pregression = regression;


  for (i = 0; i < nb_covariable ; i++) {
    *pregression++ = *iregression++;
  }

  pregression = NULL;
  delete pregression;

  prandom = NULL;
  delete prandom;

  parametric_variance_computation();
}


/*-----------------------------------------------------
 *
 * Constructeur de la classe Continuous_parametric
 *
 * arguments : nombre de covariable, nombre 
 *             d'effets aléatoires, variance résiduelle, écart-type
 *             des effets aléatoires, paramètres de régression,
 *             LM/LMM    
 *
 *-----------------------------------------------------*/

Continuous_parametric::Continuous_parametric(int inb_covariable, int inb_random, double iresidual_variance, 
					     double *irandom_variance, double *iregression, int iident)
  :Continuous_distribution(NB_VALUE_DEFAULT)
{
  register int i;
  double *pregression, *prandom;
  
  ident = iident;
  nb_covariable = inb_covariable;
  nb_random = inb_random;
  residual_variance = iresidual_variance;
  
  if (nb_random > 0){
    random_variance = new double[nb_random];
    prandom = random_variance;
    for (i = 0; i < nb_random ; i++) {
      *prandom++ = *irandom_variance++;
    }
  }
  else {
    random_variance = 0;
  }
 
  regression = new double[nb_covariable];
  
  pregression = regression;


  for (i = 0; i < nb_covariable ; i++) {
    *pregression++ = *iregression++;
  }

  pregression = NULL;
  delete pregression;

  prandom = NULL;
  delete prandom;

}


/*-----------------------------------------------------
 *
 * Constructeur de la classe Continuous_parametric
 *
 * arguments : moyenne, variance, nombre de covariable, nombre 
 *             d'effets aléatoires, variance résiduelle, écart-type
 *             des effets aléatoires, paramètres de régression,
 *             LM/LMM, nombre de valeurs, seuil de rejet    
 *
 *-----------------------------------------------------*/

Continuous_parametric::Continuous_parametric(double imean, double ivariance,
					     int inb_covariable, int inb_random, 
					     double iresidual_variance, 
					     double *irandom_variance, double *iregression, int iident, int inb_value, 
					     double icumul_threshold)
  :Continuous_distribution(NORMAL, imean, ivariance, icumul_threshold, inb_value)
{
  register int i;
  double *pregression, *prandom;
  
  ident = iident;
  nb_covariable = inb_covariable;
  nb_random = inb_random;
  residual_variance = iresidual_variance;
  
  
  if (nb_random > 0){
    random_variance = new double[nb_random];
    prandom = random_variance;
    for (i = 0; i < nb_random ; i++) {
      *prandom++ = *irandom_variance++;
    }
  }
  else {
    random_variance = 0;
  }

  regression = new double[nb_covariable];
  pregression = regression;
  for (i = 0; i < nb_covariable ; i++) {
    *pregression++ = *iregression++;
  }
  
  pregression = NULL;
  delete pregression;

  prandom = NULL;
  delete prandom;
}


/*-----------------------------------------------------
 *
 * Constructeur de la classe Continuous_parametric
 *
 * arguments : covaribles, effets aléatoires individuels, effets 
 *             aléatoires temporels, nombre de covariable, nombre 
 *             d'effets aléatoires, variance résiduelle, écart-type
 *             des effets aléatoires, paramètres de régression,
 *             LM/LMM, nombre de valeurs, seuil de rejet    
 *
 *-----------------------------------------------------*/

Continuous_parametric::Continuous_parametric(double *covar, double effect, 
					     double year_effect,
					     int inb_covariable, int inb_random, double iresidual_variance, 
					     double *irandom_variance, double *iregression, int iident, int inb_value, 
					     double icumul_threshold)
{
  double imean, ivariance;
  register int i;
  double *pregression, *prandom, *qregression, *qrandom;
  
  ident = iident;
  nb_covariable = inb_covariable;
  nb_random = inb_random;
  residual_variance = iresidual_variance;
    
  if (nb_random > 0){
    random_variance = new double[nb_random];
    prandom = random_variance;
    for (i = 0; i < nb_random ; i++) {
      *prandom++ = *irandom_variance++;
    }
  }
  else {
    random_variance = 0;
  }

  regression = new double[nb_covariable];
  
  pregression = regression;

  qregression = iregression;
  
  for (i = 0; i < nb_covariable ; i++) {
    *pregression++ = *qregression++;
  }
  
  imean = parametric_mean_computation(covar, effect, year_effect);
  ivariance = parametric_variance_computation();
  (*this) = Continuous_parametric(imean, ivariance,  inb_covariable, inb_random, iresidual_variance, 
				  irandom_variance, iregression, iident,  inb_value, icumul_threshold);

  pregression = NULL;
  delete pregression;

  prandom = NULL;
  delete prandom;

  qregression = NULL;
  delete qregression;

  qrandom = NULL;
  delete qrandom;

}


/*-------------------------------------------------------
 *
 * Constructeur de la classe Continuous_parametric
 *
 * arguments : reference sur un objet Continuous_distribution
 *
 *-------------------------------------------------------*/

Continuous_parametric::Continuous_parametric(const Continuous_distribution &cont_dist)
  :Continuous_distribution(cont_dist)
{
  ident = NO_RANDOM;
  nb_covariable = 0;
  nb_random = 0;
  residual_variance = 1.;
  random_variance = 0;
  regression = 0;

}


/*-------------------------------------------------------
 *
 * Constructeur de la classe Continuous_parametric
 *
 * arguments : reference sur un objet Continuous_histo
 *
 *-------------------------------------------------------*/
 
Continuous_parametric::Continuous_parametric(const Continuous_histo &cont_histo)
  :Continuous_distribution(cont_histo)
{
  ident = NO_RANDOM;
  nb_covariable = 0;
  nb_random = 0;
  residual_variance = 1.;
  random_variance = 0;
  regression = 0;
}


/*-------------------------------------------------------
 *
 * Constructeur de la classe Continuous_parametric
 *
 * arguments : reference sur un objet Continuous_parametric
 *
 *-------------------------------------------------------*/

Continuous_parametric::Continuous_parametric(const Continuous_parametric &cont_dist)
  :Continuous_distribution(cont_dist)
{
  copy(cont_dist);
}


/*-----------------------------------------------------
 *
 * Destructeur de la classe Continuous_parametric
 *
 *-----------------------------------------------------*/

Continuous_parametric::~Continuous_parametric()
{
  if (random_variance) {
    delete [] random_variance;
  }

  if (regression) {
    delete [] regression;
  }
}


/*---------------------------------------------------------------
 *
 * Copie d'un objet Continuous_parametric
 *
 * arguments : réference sur un objet Continuous_parametric
 *
 *---------------------------------------------------------------*/

void Continuous_parametric::copy(const Continuous_parametric &cont_dist)
{
  register int i;
  double *pregression, *prandom, *cregression, *crandom;

  ident = cont_dist.ident;
  nb_covariable = cont_dist.nb_covariable;
  nb_random = cont_dist.nb_random;
  residual_variance = cont_dist.residual_variance;


  if (nb_random > 0){
    random_variance = new double[nb_random];
    prandom = random_variance;
    crandom = cont_dist.random_variance;

    for (i = 0; i < nb_random ; i++) {
      *prandom++ = *crandom++;
    }
  }
  else {
    random_variance = 0;
  }

  regression = new double[nb_covariable];

  regression = regression;
  regression = cont_dist.regression;

  pregression = NULL;
  delete pregression;

  prandom = NULL;
  delete prandom;

  cregression = NULL;
  delete cregression;

  crandom = NULL;
  delete crandom;
}


/*--------------------------------------------------------------
 *
 * Operateur d'assignement de la classe Continuous_parametric.
 *
 * argument : reference sur un objet Continuous_parametric
 *
 *--------------------------------------------------------------*/

Continuous_parametric& Continuous_parametric::operator=(const Continuous_parametric &cont_dist)
{
  if (&cont_dist != this) {
    delete [] density;
    delete [] cumul;

    Continuous_distribution::copy(cont_dist);
    copy(cont_dist);
  }

  return *this;
}


/*--------------------------------------------------------------
 *
 * Ecriture des parametres d'une loi continue.
 *
 * argument : stream
 *
 *--------------------------------------------------------------*/

ostream& Continuous_parametric::ascii_print(ostream &os) const
{
  register int i;
  if (!ident)
    os << "NO_RANDOM " << "   " << endl;
  else
    os << "RANDOM " << "   " << endl;

  if (nb_covariable != 0) {
    os << "Number of covariables: " << nb_covariable << "   ";
  }
  if (nb_random != 0) {
    os << "Number of random effects: " << nb_random << "   ";
  }

  os << endl;

  os << "Residual variance: " << residual_variance << "   ";

  os << endl;

  if (random_variance) {
    if (nb_random == 1) {
      for (i = 0; i < nb_random; i++) {
	if ( ident == HETERO_RANDOM){
	  os << "Inter individual heterogeneity parameter " << i << " : " << sqrt(random_variance[i]) << " ; ";
	}
	else{
	  os << "Common environment parameter " << i << " : " << sqrt(random_variance[i]) << " ; ";
	}
      }
    }
    else {
      os << "Inter individual heterogeneity parameter " << i << " : " << sqrt(random_variance[i]) << " ; " << endl;
      os << "Common environment parameter " << i << " : " << sqrt(random_variance[i]) << " ; ";
    }
    os << endl;
  }  

  if (regression) {
    for (i = 0; i < nb_covariable; i++) {
      os << "Regressor " << i << " : " << regression[i] << " ; "<<endl;
    }
    os << endl;
  }  

  os << endl;

  return os;
}


/*--------------------------------------------------------------
 *
 * Ecriture des parametres d'une loi continue au format tableur.
 *
 * argument : stream
 *
 *--------------------------------------------------------------*/

ostream& Continuous_parametric::spreadsheet_print(ostream &os) const
{
  register int i;
  if (!ident)
    os << "NO_RANDOM " << "   " << endl;
  else
    os << "RANDOM " << "   " << endl;

  if (nb_covariable != 0) {
    os << "Number of covariables " << "\t" << nb_covariable << "\t";
  }
  if (nb_random != 0) {
    os << "Number of random effects " <<"\t" << nb_random << "\t";
  }
  os << endl;

  os << "Residual variance " << "\t" << residual_variance << "\t";

  os << endl;

  if (!random_variance) {
    if (nb_random == 1) {
      for (i = 0; i < nb_random; i++) {
	if ( ident == HETERO_RANDOM){
	  os << "Inter individual heterogeneity parameter " << i << "\t" << sqrt(random_variance[i]) << "\t";
	}
	else{
	  os << "Common environment parameter " << i << "\t" << sqrt(random_variance[i]) << "\t";
	}
      }
    }
    else {
      os << "Inter individual heterogeneity parameter " << i << "\t" << sqrt(random_variance[i]) << "\t" << endl;
      os << "Common environment parameter " << i << "\t" << sqrt(random_variance[i]) << "\t";
    }
    os << endl;
  }  

  if (regression) {
    for (i = 0; i < nb_covariable; i++) {
      os << "Regressor " << i << "\t" << regression[i] << "\t" << endl;
    }
    os << endl;
  }  

 

  os << endl;

  return os;
}


/*--------------------------------------------------------------
 *
 * Visualisation d'une loi continue parametrique.
 *
 * argument : stream, reference sur un objet Continuous_parametric
 *
 *--------------------------------------------------------------*/

ostream& operator<<(ostream &os , const Continuous_parametric &cont_dist)
{
  os.precision(5);

  os << endl;
  cont_dist.ascii_print(os);
  cont_dist.Continuous_distribution::print(os);

  os.precision(6);

  return os;
}


/*--------------------------------------------------------------
 *
 * Calcul du nombre de parametres d'une loi.
 *
 *--------------------------------------------------------------*/

int Continuous_parametric::nb_parameter_computation()
{
  int nb_parameter;

  switch (ident) {
  case NO_RANDOM :
    nb_parameter = nb_covariable + 1;
    break;
  case HETERO_RANDOM : 
    nb_parameter = nb_covariable + nb_random + 1;
    break;
  case YEAR_RANDOM : 
    nb_parameter = nb_covariable + nb_random + 1;
    break;
  }

  return nb_parameter;
}


/*--------------------------------------------------------------
 *
 * Calcul de la moyenne d'une loi parametrique gasussienne pour un individu.
 *
 * arguments : covariables, effets aléatoires individuels,
 *             effets aléatoires temporels
 *
 *--------------------------------------------------------------*/

double Continuous_parametric::parametric_mean_computation(double *covar, double effect, 
							  double year_effect) const
{
  register int i;
  double parametric_mean = 0.;

  switch (ident) {
  case NO_RANDOM :
    for (i = 0; i < nb_covariable; i++) {
      parametric_mean = parametric_mean + covar[i] * regression[i];   
    }
    break;
    
  case HETERO_RANDOM :
    for (i = 0; i < nb_covariable; i++) {
      parametric_mean = parametric_mean + covar[i] * regression[i];
    }
    if (effect){
      parametric_mean = parametric_mean + sqrt(random_variance[0]) * effect;
    }
    break;
  
  case YEAR_RANDOM :
    for (i = 0; i < nb_covariable; i++) {
      parametric_mean = parametric_mean + covar[i] * regression[i];
    }
    if (year_effect){
      if (nb_random == 1){
	parametric_mean = parametric_mean + sqrt(random_variance[0]) * year_effect;
      }
      else{
	parametric_mean = parametric_mean + sqrt(random_variance[1]) * year_effect;
      }
    }
    break;
  
  default :
    parametric_mean = mean;
    break;
  }

  return parametric_mean;
}


/*--------------------------------------------------------------
 *
 * Calcul de la moyenne marginale d'une loi parametrique gasussienne pour un individu.
 *
 * arguments : covariables
 *
 *--------------------------------------------------------------*/

double Continuous_parametric::parametric_mean_computation_marginale(double *covar) const
{
  register int i;
  double parametric_mean = 0.;

  switch (ident) {
  case NO_RANDOM :
    for (i = 0; i < nb_covariable; i++) {
      parametric_mean = parametric_mean + covar[i] * regression[i];   
    }
    break;
  case HETERO_RANDOM :
    for (i = 0; i < nb_covariable; i++) {
      parametric_mean = parametric_mean + covar[i] * regression[i];
    }
    break;
  case YEAR_RANDOM :
    for (i = 0; i < nb_covariable; i++) {
      parametric_mean = parametric_mean + covar[i] * regression[i];
    }
    break;
  default :
    parametric_mean = mean;
    break;
  }

  return parametric_mean;
}


/*--------------------------------------------------------------
 *
 * Calcul de la moyenne d'une loi parametrique continue pour n individus.
 * pas implémenté dans le cas avec effets aléatoires
 *
 * arguments : nombre d'individus, covariables
 *
 *--------------------------------------------------------------*/

double* Continuous_parametric::n_parametric_mean_computation(int n, double **covar) const
{
  register int i;
  register int j;
  double sum;
  double* parametric_mean;

  parametric_mean = new double[n];

  switch (ident) {
  case NO_RANDOM :
    for (i = 0; i < n; i++) {
      sum = 0.;
      for(j = 0; j <nb_covariable; j++){
	sum = sum + covar[i][j] * regression[j];    
      }
      parametric_mean[i] = sum;
    }
    break;
  case HETERO_RANDOM :
    break;
  default :   
    for (i = 0; i < n; i++) {
      parametric_mean[i] = mean;
    }
    break;
  }

  return parametric_mean;
}


/*--------------------------------------------------------------
 *
 *  Calcul de la variance d'une loi parametrique continue.
 *
 *--------------------------------------------------------------*/

double Continuous_parametric::parametric_variance_computation() const
{
  double parametric_variance;

  switch (ident) {
  case NO_RANDOM :
    parametric_variance = residual_variance;
    break;
  case HETERO_RANDOM :
    parametric_variance = residual_variance;
    break;
  case YEAR_RANDOM :
    parametric_variance = residual_variance;
    break;
  default :
    parametric_variance = variance;
    break;
  }

  return parametric_variance;
}


/*--------------------------------------------------------------
 *
 *  Calcul de la variance marginale d'une loi parametrique continue.
 *
 *--------------------------------------------------------------*/

double Continuous_parametric::parametric_variance_computation_marginale() const
{
  register int i;
  double parametric_variance;

  switch (ident) {
  case NO_RANDOM :
    parametric_variance = residual_variance;
    break;
  case HETERO_RANDOM :
    parametric_variance = residual_variance;
    for (i = 0; i < nb_random; i++) {
      parametric_variance = parametric_variance + random_variance[i];    
    }
    break;
  case YEAR_RANDOM :
    parametric_variance = residual_variance;
    for (i = 0; i < nb_random; i++) {
      parametric_variance = parametric_variance + random_variance[i];    
    }
    break;
  default :
    parametric_variance = variance;
    break;
  }

  return parametric_variance;
}


/*-----------------------------------------------------------
 *
 * Calcul d'une loi normale
 *
 * arguments : covariables, effets aléatoires individuels,
 *             effets aléatoires temporels, seuil de rejet
 *
 *-----------------------------------------------------------*/

void Continuous_parametric::normal_computation(double *covar, double effect, double year_effect, double icumul_threshold)
{
  register int i;
  double step, proba , *pdensity;

  type = NORMAL;
  cumul_threshold = icumul_threshold;
  nb_value = NB_VALUE_DEFAULT; 

  pdensity = density;

  mean = (*this).parametric_mean_computation(covar, effect, year_effect);
  variance = (*this).parametric_variance_computation(); 
  qnorm();
  min_value = - (*this).max_value * sqrt(variance) + mean;
  max_value = (*this).max_value * sqrt(variance) + mean;
  step = (max_value - min_value) / nb_value;

  for (i = 0;i <= nb_value;i++) {
    proba = ( 1. / ( sqrt ( 2. * S_PI * variance))) *
      exp( - (((min_value + i * step)-mean)*((min_value + i * step)-mean)) / (2 * variance)); 
    
    *pdensity++ = proba;
  }
  cumul_cont_computation();

  pdensity = NULL;
  delete pdensity;

}


/*-----------------------------------------------------------
 *
 * Calcul d'une loi normale marginale
 *
 * arguments : covariables, seuil de rejet
 *
 *-----------------------------------------------------------*/

void Continuous_parametric::normal_computation_marginale(double *covar, double icumul_threshold)
{
  register int i;
  double step, proba , *pdensity;

  type = NORMAL;
  cumul_threshold = icumul_threshold;
  nb_value = NB_VALUE_DEFAULT; 

  pdensity = density;

  mean = (*this).parametric_mean_computation_marginale(covar);
  variance = (*this).parametric_variance_computation_marginale(); 
  qnorm();
  min_value = - (*this).max_value * sqrt(variance) + mean;
  max_value = (*this).max_value * sqrt(variance) + mean;
  step = (max_value - min_value) / nb_value;

  for (i = 0;i <= nb_value;i++) {
    proba = ( 1. / ( sqrt ( 2. * S_PI * variance))) *
      exp( - (((min_value + i * step)-mean)*((min_value + i * step)-mean)) / (2 * variance)); 
    *pdensity++ = proba;
  }

  cumul_cont_computation();


  pdensity = NULL;
  delete pdensity;
}


/*-----------------------------------------------------------
 *
 * Calcul de la densité d'une observation de moyenne X\beta 
 * et de variance 
 *
 * arguments : observation, covariables, effets aléatoires individuels
 *             effets aléatoires temporels
 *
 *-----------------------------------------------------------*/

double Continuous_parametric::density_computation(double observ, double *covar, double effect, double year_effect)
{
  double pdensity;

  mean = (*this).parametric_mean_computation(covar, effect, year_effect);
  variance = (*this).parametric_variance_computation();
  
  switch (type) {
  case NORMAL :
    pdensity = ( 1. / ( sqrt ( 2. * S_PI * variance))) * exp( - ((observ - mean)*(observ - mean)) / (2 * variance)); 
    break;
  case EXPONENTIAL :
    break;
  default : 
    break;
  }
  
  //  mean = (*this).parametric_mean_computation_marginale(covar);
  // variance = (*this).parametric_variance_computation_marginale();

  return pdensity;
}


/*-----------------------------------------------------------
 *
 * Calcul de la densité marginale d'une observation de moyenne X\beta 
 * et de variance 
 *
 * arguments : observation, covariables
 *
 *-----------------------------------------------------------*/

double Continuous_parametric::density_computation_marginale(double observ, double *covar)
{
  double pdensity;

  mean = (*this).parametric_mean_computation_marginale(covar);
  variance = (*this).parametric_variance_computation_marginale();
  
  switch (type) {
  case NORMAL :
    pdensity = ( 1. / ( sqrt ( 2. * S_PI * variance))) * exp( - ((observ - mean)*(observ - mean)) / (2 * variance)); 
    break;
  case EXPONENTIAL :
    break;
  default : 
    break;
  }
  
  return pdensity;
}


/*-----------------------------------------------------------
 *
 * Calcul de la fonction de répartition d'une observation de moyenne X\beta 
 * et de variance 
 *
 * arguments : observation, covariables, effets aléatoires
 *             individuels, effets aléatoires temporels
 *
 *-----------------------------------------------------------*/

double Continuous_parametric::cumul_computation(double observ, double *covar, double effect, double year_effect)
{
  double pcumul;

  mean = (*this).parametric_mean_computation(covar, effect, year_effect);
  variance = (*this).parametric_variance_computation();
  
  switch (type) {
  case NORMAL :
    pcumul = pnorm((observ - mean)/sqrt(variance)); 
    break;
  case EXPONENTIAL :
    break;
  default : 
    break;
  }  

  return pcumul;
}


/*-----------------------------------------------------------
 *
 * Calcul de la fonction de répartition marginale d'une observation de moyenne X\beta 
 * et de variance 
 *
 * arguments : observation, covariables
 *
 *-----------------------------------------------------------*/

double Continuous_parametric::cumul_computation_marginale(double observ, double *covar)
{
  double pcumul;

  mean = (*this).parametric_mean_computation_marginale(covar);
  variance = (*this).parametric_variance_computation_marginale();
  
  switch (type) {
  case NORMAL :
    pcumul = pnorm((observ - mean)/sqrt(variance)); 
    break;
  case EXPONENTIAL :
    break;
  default : 
    break;
  }
  
  return pcumul;
}


/*-----------------------------------------------------------
 *
 * Tirage aléatoire selon une distribution continue paramètrique
 *
 * arguments : covariables
 *
 *-----------------------------------------------------------*/

double Continuous_parametric::simulation(double *covar)
{
  double temp;
  
  mean = (*this).parametric_mean_computation_marginale(covar);
  variance = (*this).parametric_variance_computation();
  
  switch (type) {
  case NORMAL :
    temp = random_normal(); 
    break;
  case EXPONENTIAL :
    temp = random_expo();
    break;
  default : 
    break;
  }
  
  return temp;
}


/*-----------------------------------------------------------
 *
 * Tirage aléatoire selon une distribution continue paramètrique
 *
 * arguments : covariables
 *
 *-----------------------------------------------------------*/
double Continuous_parametric::simulation_marginale(double *covar)
{
  double temp;
  
  mean = (*this).parametric_mean_computation_marginale(covar);
  variance = (*this).parametric_variance_computation_marginale();
  
  switch (type) {
  case NORMAL :
    temp = random_normal(); 
    break;
  case EXPONENTIAL :
    temp = random_expo();
    break;
  default : 
    break;
  }
  
  return temp;
}


/*-------------------------------------------------------------
 *
 * Calcul des paramètres de régression par
 * la méthode des moindres carrés ordinaires (MCO)
 *
 * arguments : nombre d'observations, vecteur des observations, 
 *             matrice des covariables, matrice de pondération
 *
 *-------------------------------------------------------------*/
 
double* Continuous_parametric::ordinary_least_squares(int nb_observ, double *observ, double **covar, double **weight)
{
  int i, s;

  double *parameter;
  parameter = new double[nb_covariable];

  gsl_permutation *p = gsl_permutation_alloc(nb_covariable);
  gsl_matrix *trans = gsl_matrix_calloc(nb_covariable, nb_covariable);      
  gsl_vector *beta = gsl_vector_calloc(nb_covariable);
  gsl_vector *Y = gsl_vector_calloc(nb_observ);
  gsl_matrix *X = gsl_matrix_calloc( nb_observ, nb_covariable);
  gsl_matrix *gamma = gsl_matrix_calloc(nb_observ, nb_observ);
  gsl_matrix *Xgamma = gsl_matrix_calloc(nb_covariable, nb_observ);
  gsl_matrix *XgammaX = gsl_matrix_calloc(nb_covariable, nb_covariable);
  gsl_matrix *inv_trans = gsl_matrix_calloc(nb_covariable, nb_covariable);
  gsl_vector *XgammaY = gsl_vector_calloc(nb_covariable);

  convert_matrix_double(covar, nb_observ, nb_covariable, X);

  if(!weight){
    gsl_matrix_set_identity(gamma);
  }
  else {
    convert_matrix_double(weight, nb_observ, nb_observ, gamma);
  }

  convert_vector_double(observ, nb_observ, Y);
  
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1., X, gamma, 0., Xgamma); // t_Xgamma
  
  
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1., Xgamma, X, 0., XgammaX); // t_XgammaX
  
  gsl_blas_dgemv(CblasNoTrans, 1., Xgamma, Y, 0., XgammaY);//t_XgammaY

	    
  gsl_matrix_memcpy(trans, XgammaX);
  gsl_linalg_LU_decomp(trans, p, &s);
  gsl_linalg_LU_invert(trans, p, inv_trans); //t_XgammaX^{-1}
  
  gsl_blas_dgemv(CblasNoTrans, 1., inv_trans, XgammaY, 0., beta); // t_XgammaX-1} t_XgammaY

  for (i = 0; i < nb_covariable; i++){
    parameter[i] = gsl_vector_get(beta,i);
  }  

  gsl_vector_free(XgammaY);
  gsl_matrix_free(inv_trans);
  gsl_matrix_free(XgammaX);
  gsl_matrix_free(Xgamma);
  gsl_matrix_free(gamma);
  gsl_matrix_free(X);
  gsl_vector_free(Y);
  gsl_vector_free(beta);
  gsl_matrix_free(trans);
  gsl_permutation_free(p);

  return parameter;

} 


/*-------------------------------------------------------------
 *
 * Calcul de la variance résiduelle par 
 * maximum de vraisemblance (cas gaussien)
 *
 * arguments : nombre d'observations, vecteur des observations, 
 *             matrice des covariables, matrice de pondération
 *
 *-------------------------------------------------------------*/
 
double Continuous_parametric::variance_estimation(int nb_observ, double *observ, double **covar, double **weight)
{
  int i, s;

  double parameter;
  double num = 0.;
  double trace = 0.;
    
  gsl_vector *beta = gsl_vector_calloc(nb_covariable);
  gsl_vector *Y = gsl_vector_calloc(nb_observ);
  gsl_matrix *X = gsl_matrix_calloc( nb_observ, nb_covariable);
  gsl_matrix *gamma = gsl_matrix_calloc(nb_observ, nb_observ);
  gsl_vector *Xbeta = gsl_vector_calloc(nb_observ);
  gsl_vector *gammasub = gsl_vector_calloc(nb_observ);

  convert_matrix_double(covar, nb_observ, nb_covariable, X);

  if(!weight){
    gsl_matrix_set_identity(gamma);
  }
  else {
    convert_matrix_double(weight, nb_observ, nb_observ, gamma);
  }
  
  convert_vector_double(observ, nb_observ, Y);

  convert_vector_double(regression, nb_covariable, beta);

  gsl_blas_dgemv(CblasNoTrans, 1., X, beta, 0., Xbeta);//Xbeta

  gsl_vector_sub(Y, Xbeta);//Y-Xbeta
  
  gsl_blas_dgemv(CblasNoTrans, 1., gamma, Y, 0., gammasub); // gamma(Y-Xbeta)
  
  for (i = 0; i < nb_observ; i++){
    num = num + gsl_vector_get(Y,i) * gsl_vector_get(gammasub,i);
  }
 
  for (i = 0; i < nb_observ; i++){
    trace = trace + gsl_matrix_get(gamma,i,i);
  }

  parameter = num/trace;

  gsl_vector_free(gammasub);
  gsl_vector_free(Xbeta);
  gsl_matrix_free(gamma);
  gsl_matrix_free(X);
  gsl_vector_free(Y);
  gsl_vector_free(beta);
 
  return parameter;

} 


/*------------------------------------------------------------
 *
 * Mise à jour de la variance residuelle
 *
 * arguments : variance résiduelle
 *
 *------------------------------------------------------------*/

void Continuous_parametric::Set_residual_variance(const double iresidual_variance) 
{ 
  residual_variance = iresidual_variance; 

  this->variance = parametric_variance_computation();
  this->min_value_computation();
  this->max_value_computation();
  //this->nb_value_computation();
  //  this->density_cont_computation();
//   this->cumul_cont_computation();
}


/*------------------------------------------------------------
 *
 * Mise à jour des écart-types des effets aléatoires
 *
 * arguments : écart-types des effets aléatoires
 *
 *------------------------------------------------------------*/
void Continuous_parametric::Set_random_variance(const double *irandom_variance)
{ 
  register int i;

  double *prandom_variance;
  random_variance = new double[nb_random];

  prandom_variance = random_variance;

  for (i = 0; i < nb_random ; i++) {
    *prandom_variance++ = *irandom_variance++;
  }

  //this->variance = parametric_variance_computation();
  this->min_value_computation();
  this->max_value_computation();
  //this->nb_value_computation();
  //  this->density_cont_computation();
  // this->cumul_cont_computation();

  prandom_variance = NULL;
  delete prandom_variance;

}


/*-------------------------------------------------------------
 *
 * Mise à jour des paramètres de régression fixes
 *
 * arguments : paramètres de régression
 *
 *-------------------------------------------------------------*/

void Continuous_parametric::Set_regression(const double *iregression) 
{ 
  register int i;
  double *pregression;
  regression = new double[nb_covariable];

  pregression = regression;

  for (i = 0; i < nb_covariable ; i++) {
    *pregression++ = *iregression++;
  }

  pregression = NULL;
  delete pregression;

}
