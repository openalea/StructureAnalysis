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
 *       $Id: distribution_reestimation.hpp 18000 2015-04-23 06:56:37Z guedon $
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



#ifndef DISTRIBUTION_REESTIMATION_HPP
#define DISTRIBUTION_REESTIMATION_HPP



namespace stat_tool {


template <typename Type>
void reestimation(int nb_value , Type *reestim , double *pmass ,
                  double min_probability , bool null_probability);



/*--------------------------------------------------------------*
 *
 *  Reestimation des probabilites d'une loi discrete.
 *
 *  arguments : nombre de valeurs, pointeurs sur les quantites
 *              de reestimation et sur les probabilites de la loi discrete,
 *              probabilite minimum, flag pour reestimer des probabilites nulles.
 *
 *--------------------------------------------------------------*/

template <typename Type>
void reestimation(int nb_value , Type *reestim , double *pmass ,
                  double min_probability , bool null_probability)

{
  register int i;
  int nb_correction;
  double sum , norm;


  sum = 0.;
  for (i = 0;i < nb_value;i++) {
    sum += *reestim++;
  }
  reestim -= nb_value;

  if (sum > 0.) {
    norm = 0.;
    nb_correction = 0;

    for (i = 0;i < nb_value;i++) {
      if (*pmass > 0.) {
        if (*reestim > min_probability * sum) {
          *pmass = *reestim / sum;
          norm += *pmass;
        }

        else {
          if ((*reestim > 0.) || (!null_probability)) {
            nb_correction++;
            *pmass = min_probability;
          }
          else {
            *pmass = 0.;
          }
        }
      }

      reestim++;
      pmass++;
    }

    if (nb_correction > 0) {
      reestim -= nb_value;
      pmass -= nb_value;

      for (i = 0;i < nb_value;i++) {
        if ((*pmass > 0.) && (*reestim > min_probability * sum)) {
          *pmass *= (1. - nb_correction * min_probability) / norm;
        }
        reestim++;
        pmass++;
      }
    }
  }

  else {
    nb_correction = 0;
    for (i = 0;i < nb_value;i++) {
      if (*pmass++ > 0.) {
        nb_correction++;
      }
    }
    pmass -= nb_value;

    for (i = 0;i < nb_value;i++) {
      if (*pmass > 0.) {
        *pmass = 1. / (double)nb_correction;
      }
      pmass++;
    }
  }
}


};  // namespace stat_tool



#endif
