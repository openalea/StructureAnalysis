/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       AMAPmod: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2000 UMR Cirad/Inra Modelisation des Plantes
 *
 *       File author(s): P.ferraro (pascal.ferraro@cirad.fr)
 *
 *       $Source$
 *       $Id$
 *
 *       Forum for AMAPmod developers    : amldevlp@cirad.fr
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



#ifndef SB_END_SPACE_FREE_MATCHING_U_HEADER
#define SB_END_SPACE_FREE_MATCHING_U_HEADER

#include "matching_unordered.h"


/**
 *\class EndSpaceFreeMatching 
 *\brief Algorithm for evaluating End Space Free Mathcing
 *\par Presentation
 * In order to compute similarity between two tree graphs, we have extended
 * the algorithm proposed by Zhang \cite{Zha93} to define a distance that not consider any space
 * at the beginnig or the end of tree matching.
 *\par Requirements
 * - Two TreeGraphs defined at the same scale;
 * - A NodeCost (method for computing the local distance),
 *\author Pascal Ferraro
 *\date 2005
 */

class EndSpaceFreeMatching_U : public Matching_U
{

  public :
  /* Only beginning space are implemented */
    EndSpaceFreeMatching_U(TreeGraph& , TreeGraph& ,NodeCost& ) ;
    DistanceType  match();
    DistanceType  match_begin_T2();
	DistanceType  match_begin_end_T2();
 
};


#endif

