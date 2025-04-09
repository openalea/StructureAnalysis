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


#ifndef CONTINUOUS_DISTRIBUTION_H
#define CONTINUOUS_DISTRIBUTION_H

#include <fstream>
#include "continuous_histo.h"

const double CUMUL_DEFAULT = 1.e-6; 
const int NB_VALUE_DEFAULT = 100000;
const double S_PI = 3.141592653;



enum {
  NORMAL,
  EXPONENTIAL
};

class Continuous_histo;

class Continuous_distribution { 

  friend std::ostream& operator<<(std::ostream&, const Continuous_distribution&);


  public:
  int type;                     // continuous distribution type (NORMAL/EXPONENTIAL)
  double mean;                  // mean
  double variance;              // variance
  double *density;              // density function
  double *cumul;                // cumulative distribution function
  int nb_value;                 // number of elements
  double min_value;             // lower bound
  double max_value;             // upper bound
  double cumul_threshold;       // unexplained part of the cumulative distribution function de la fonction de repartition rejettee 
                                // (must be between 0. and 1.)

  Continuous_distribution();
  Continuous_distribution (int inb_value);
  Continuous_distribution(int itype, double imean, double ivariance, double icumul_threshold = CUMUL_DEFAULT, 
			  int inb_value = NB_VALUE_DEFAULT );// see default paremeter (standard Gaussian distribution)
  Continuous_distribution(double imin_value,double imax_value, int itype, double icumul_threshold = CUMUL_DEFAULT, 
			  int inb_value = NB_VALUE_DEFAULT);
  Continuous_distribution(const Continuous_distribution &cont_dist);
  Continuous_distribution(const Continuous_histo &cont_histo); 
 
  ~Continuous_distribution();

  void copy(const Continuous_distribution &cont_dist);

  Continuous_distribution& operator=(const Continuous_distribution&); 
  bool operator==(const Continuous_distribution &cont_dist) const;
  bool operator!=(const Continuous_distribution &cont_dist) const {return !(*this == cont_dist); }
  std::ostream& ascii_characteristic_print(std::ostream &os, bool comment_flag = false) const; 
  std::ostream& spreadsheet_characteristic_print(std::ostream &os) const;
  std::ostream& ascii_print(std::ostream &os , bool comment_flag , bool density_flag, bool cumul_flag ,
			    bool nb_class_flag , const Continuous_histo *cont_histo) const;
  std::ostream& ascii_write(std::ostream &os , bool exhaustive , bool file_flag, 
			    const Continuous_histo *cont_histo, double icumul_threshold = CUMUL_DEFAULT) const;
  std::ostream& spreadsheet_print(std::ostream &os , bool cumul_flag , bool density_flag ,
				  bool nb_class_flag , const Continuous_histo *cont_histo) const;
  std::ostream& spreadsheet_write(std::ostream &os , bool exhaustive , bool file_flag, 
				  const Continuous_histo *cont_histo, double icumul_threshold = CUMUL_DEFAULT) const;

  bool plot_print(const char *path , double *concentration ,double scale) const;
  
  std::ostream& print(std::ostream&) const;

  void qnorm();
  void qexp();
  void pnorm();
  void pexp();
  double pnorm(double x);
  double pexp(double x, double imean = 1.);
  double qnorm(double p);
  double qexp(double p);
  double dnorm(double d, double imean = 0., double ivariance = 1.);
  double dexp(double d, double imean = 1.);

  void min_value_computation();
  void max_value_computation();
  double max_computation();

  void nb_value_computation(); //does not work

  void mean_computation();
  void variance_computation();
  void density_cont_computation();
  void cumul_cont_computation();
  void log_computation();
  double random_normal();
  double random_expo();

  double likelihood_computation(); 
  double likelihood_computation(const Continuous_histo &cont_histo, double icumul_threshold = CUMUL_DEFAULT ) const ; // does not work



  //private:

  //protected:

};


#endif
