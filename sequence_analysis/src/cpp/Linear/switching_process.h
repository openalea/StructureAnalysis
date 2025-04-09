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


#ifndef SWITCHING_PROCESS_H
#define SWITCHING_PROCESS_H

#include <sstream>
#include "continuous_parametric.h"


const double SWITCHING_NOISE_PROBABILITY = 0.05; // disturbance of observation densities 

class Switching_process {

  friend std::ostream& operator<<(std::ostream &os, const Switching_process &sw_process)
  { return sw_process.ascii_print(os, 0, false, false);}

 public:
  int nb_state;                                // number of states
  int nb_value;                                // number of values of output process
  Continuous_parametric **observation;         // set of linear mixed models 

  Switching_process(int inb_state = 0, int inb_value = 0);
  Switching_process(int inb_state, Continuous_parametric **pobservation);
  Switching_process(const Switching_process &process, char manip = 'c', int state = -1);
  ~Switching_process();

  Switching_process& operator=(const Switching_process &process);

  void copy(const Switching_process  &process);
  void add_state(const Switching_process &process, int state);
  void remove();

  std::ostream& ascii_print(std::ostream &os, Continuous_histo **empirical_observation = 0 ,
			    bool exhaustive = true, bool file_flag = true) const;
  std::ostream& spreadsheet_print(std::ostream &os , Continuous_histo **empirical_observation = 0) const;

  void nb_value_computation();
  int nb_parameter_computation() const;

  void init();

  // private:
  // protected:
};

#endif
