/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       AMAPmod: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2000 UMR Cirad/Inra Modelisation des Plantes
 *
 *       File author(s): A.ouangraoua (ouangrao@labri.fr)
 *
 *       $Source: /usr/cvsmaster/AMAPmod/src/TreeMatching/POrdered/matching_PO_end_space.h,v $
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


/*matching ne tenant pas compte du nombre de composantes connexes*/

#ifndef SB_MATCHING_PO_END_SPACE_HEADER
#define SB_MATCHING_PO_END_SPACE_HEADER
#include <vector>

#include "matching_PO.h"

typedef std::vector<int> CaseVector;

/**
 *\class Matching_PO_End_Space With Complex
 *\brief Algorithm for comparing two partially ordered tree graph
 *\par Presentation
 * In order to compute a distance between two  partially ordered tree graph with end space free,
 *  we have combined the algorithm to compute edit distance between partially ordered trees and
 * the algorithm to compute end space free edit distance between trees 
 *\par Requirements
 * - Two TreeGraphs defined at the same scale;
 * - A NodeCost (method for computing the local distance),
 *\author Aida ouangraoua
 *\date 2005
 */

typedef  std::vector<std::vector<DistanceType> > DistanceTab;
typedef  std::vector<DistanceType> InDelTab;


class Matching_PO_End_Space : public Matching_PO
{
  
public :
  //Constructor
  Matching_PO_End_Space(TreeGraph& , TreeGraph& ,NodeCost& ) ;
  DistanceType  match();

};


#endif

