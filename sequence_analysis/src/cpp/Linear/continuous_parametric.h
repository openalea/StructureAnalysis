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


#ifndef CONTINUOUS_PARAMETRIC_H
#define CONTINUOUS_PARAMETRIC_H

#include <sstream>
#include "continuous_distribution.h"

enum {
  NO_RANDOM,
  HETERO_RANDOM,
  YEAR_RANDOM
};


class Continuous_parametric : public Continuous_distribution {

  friend std::ostream& operator<<(std::ostream &os , const Continuous_parametric &cont_dist);
  
 public:
  int nb_covariable;             // number of covariates (without intercept)
  int nb_random;                 // number of random effects
  int ident;                     // type of distribution or parameter (NO_RANDOM/HETERO_RANDOM/YEAR_RANDOM)
  double residual_variance;      // residual variance
  double *random_variance;       // variance of random effects (weigth tau*tau)
  double *regression;            // regression parameter

  Continuous_parametric();
  Continuous_parametric(int inb_value, int inb_covariable, int inb_random, double iresidual_variance, 
			double *irandom_variance, double *iregression, int iident = NO_RANDOM);
  Continuous_parametric(int inb_covariable, int inb_random, double iresidual_variance, 
			double *irandom_variance, double *iregression, int iident = NO_RANDOM);
  Continuous_parametric(double *covar, double effect, double year_effect, 
			int inb_covariable, int inb_random, double iresidual_variance, 
			double *irandom_variance, double *iregression, int iident = NO_RANDOM,
			int inb_value = NB_VALUE_DEFAULT, double icumul_threshold = CUMUL_DEFAULT);
  Continuous_parametric(double imean, double ivariance, int inb_covariable, int inb_random, 
			double iresidual_variance, double *irandom_variance, double *iregression, 
			int iident = NO_RANDOM, int inb_value = NB_VALUE_DEFAULT, double icumul_threshold = CUMUL_DEFAULT);
  Continuous_parametric(const Continuous_distribution &cont_dist);
  Continuous_parametric(const Continuous_histo &cont_histo);
  Continuous_parametric(const Continuous_parametric &cont_dist);
  ~Continuous_parametric();

  void copy(const Continuous_parametric &cont_dist);

  Continuous_parametric& operator=(const Continuous_parametric &cont_dist);
  std::ostream& ascii_print(std::ostream &os) const;
  std::ostream& spreadsheet_print(std::ostream &os) const;

  int nb_parameter_computation();
  double parametric_mean_computation(double *covar, double effect, double year_effect) const;
  double parametric_mean_computation_marginale(double *covar) const;

  double* n_parametric_mean_computation(int n, double **covar) const; //not implemented in random case

  double parametric_variance_computation() const;
  double parametric_variance_computation_marginale() const;

  void normal_computation(double *covar, double effect, double year_effect, double icumul_threshold = CUMUL_DEFAULT);
  void normal_computation_marginale(double *covar, double icumul_threshold = CUMUL_DEFAULT);

  double density_computation(double observ, double *covar, double effect = 0, double year_effect = 0);
  double density_computation_marginale(double observ, double *covar);

  double cumul_computation(double observ, double *covar, double effect, double year_effect);
  double cumul_computation_marginale(double observ, double *covar);

  double simulation(double *covar);
  double simulation_marginale(double *covar);


  double* ordinary_least_squares(int nb_observ, double *observ, double **covar, double **weight = 0); 
  double variance_estimation(int nb_observ, double *observ, double **covar, double **weight = 0);

  // modifications des attributs

  void Set_residual_variance(const double iresidual_variance);
  void Set_random_variance(const double *irandom_variance);
  void Set_regression(const double *iregression);

  // private:
  // protected:
};

#endif
