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

#ifndef CONTINUOUS_HISTO_H
#define CONTINUOUS_HISTO_H

#include <fstream>

#include "continuous_distribution.h"

class Continuous_distribution;

class Continuous_histo { 

  friend std::ostream& operator<<(std::ostream &os, const Continuous_histo &cont_histo)
  { return cont_histo.print(os); }


  public:
  double mean;                            // mean
  double variance;                        // variance
  double *frequency;                      // value frequency in a range 
  int nb_element;                         // total frequency
  int nb_class;                           // number of classes 
  double min_value;                       // lower bound
  double max_value;                       // upper bound

  Continuous_histo();
  Continuous_histo(int inb_class);
  Continuous_histo(int inb_element, double *pelement);
  Continuous_histo(const Continuous_histo &cont_histo);
  Continuous_histo(const Continuous_distribution &cont_dist);
  Continuous_histo(const Continuous_histo &cont_histo, double area);
  Continuous_histo(double param, const Continuous_histo &cont_histo);
  ~Continuous_histo();

  void nb_class_computation();
  void mean_computation();
  void variance_computation();
  void nb_element_computation();

  void copy(const Continuous_histo &cont_histo);

  Continuous_histo& operator=(const Continuous_histo &cont_histo);

  std::ostream& ascii_characteristic_print(std::ostream &os , bool comment_flag=false) const;
  std::ostream& spreadsheet_characteristic_print(std::ostream &os) const;
  std::ostream& spreadsheet_print(std::ostream &os , bool cumul_flag, bool density_flag) const;
  std::ostream& print(std::ostream &os) const;
  std::ostream& ascii_write(std::ostream &os , bool exhaustive , bool file_flag) const;
  std::ostream& ascii_print(std::ostream &os , int comment_flag , bool cumul_flag, bool density_flag) const;
  bool ascii_write(const char *path) const;

  double information_computation() const;// does not work (voir feuille imprimee sur l'entropie)
  double likelihood_computation(double icumul_threshold ) const;
  double likelihood_computation(Continuous_distribution &cont_dist) const; 

  void distribution_estimation(Continuous_distribution *cont_dist) const;

  bool operator==(const Continuous_histo&) const;
  bool operator!=(const Continuous_histo &cont_histo) const {return !(*this == cont_histo); }

  void shift(const Continuous_histo &cont_histo , double shift_param);

  double* density_computation()const;
  double* cumul_computation()const;

  //private:

  //protected:

};


#endif
