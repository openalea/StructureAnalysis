/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       V-Plants: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2017 CIRAD/INRA/Inria Virtual Plants
 *
 *       File author(s): Yann Guedon (yann.guedon@cirad.fr)
 *
 *       $Source$
 *       $Id: quantile_computation.hpp 16084 2014-03-17 14:50:26Z guedon $
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



#ifndef QUANTILE_COMPUTATION_HPP
#define QUANTILE_COMPUTATION_HPP


#include <math.h>
#include <iostream>

using namespace std;


namespace stat_tool {


template <typename Type>
double quantile_computation(int nb_individual , Type *value , double cumul);



/*--------------------------------------------------------------*/
/**
 *  \brief Computation of a quantile for a discrete or continuous sample.
 *
 *  \param[in] nb_individual sample size,
 *  \param[in] value         pointer on the values,
 *  \param[in] cumul         value of the cumulative distribution function.
 *
 *  \return                  quantile.
 */
/*--------------------------------------------------------------*/

template <typename Type>
double quantile_computation(int nb_individual , Type *value , double cumul)

{
  int i , j;
  int *index;
  Type bvalue , min_value , max_value;
  double quantile;


  index = new int[nb_individual];

  max_value = value[0];
  for (i = 1;i < nb_individual;i++) {
    if (value[i] > max_value) {
      max_value = value[i];
    }
  }

  i = 0;

  do {

    // determination of the current minimum value

    if (i == 0) {
      min_value = value[0];
      for (j = 1;j < nb_individual;j++) {
        if (value[j] < min_value) {
          min_value = value[j];
        }
      }
    }

    else {
      bvalue = max_value + 1;
      for (j = 0;j < nb_individual;j++) {
        if ((value[j] > min_value) && (value[j] < bvalue)) {
          bvalue = value[j];
        }
      }
      min_value = bvalue;
    }

    // selection of the individuals taking the current minimum value

    for (j = 0;j < nb_individual;j++) {
      if (value[j] == min_value) {
        index[i++] = j;
      }
    }
  }
  while (i < nb_individual);

  // quantile extraction 

  if (cumul < 1.) {
    i = (int)(nb_individual * cumul);
  }
  else {
    i = nb_individual - 1;
  }
  quantile = value[index[i]];

  delete [] index;

  return quantile;
}


};  // namespace stat_tool



#endif
