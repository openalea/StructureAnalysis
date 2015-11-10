/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       AMAPmod: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2000 UMR Cirad/Inra Modelisation des Plantes
 *
 *       File author(s): A.ouangraoua (ouangrao@labri.fr)
 *
 *       $Source: /usr/cvsmaster/AMAPmod/src/TreeMatching/POrdered/matching_PO_local.h,v $
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

#ifndef SB_MATCHING_PO_LOCAL_HEADER
#define SB_MATCHING_PO_LOCAL_HEADER
#include <vector>

#include "matching_PO.h"
#include "../localmatchpath.h"

typedef std::vector<int> CaseVector;

/**
 *\class Matching_PO_Local With Complex
 *\brief Algorithm for comparing two partially ordered tree graph
 *\par Presentation
 * In order to compute a distance between two  partially ordered tree graphs locally,
 *  we have extended the algorithm to compute edit local similarity between unordered trees and
 * the algorithm to compute edit distance between partially ordered trees 
 *\par Requirements
 * - Two TreeGraphs defined at the same scale;
 * - A NodeCost (method for computing the local distance),
 *\author Aida ouangraoua
 *\date 2005
 */

typedef  std::vector<std::vector<DistanceType> > DistanceTab;
typedef  std::vector<DistanceType> InDelTab;


class Matching_PO_Local : public Matching_PO
{
  
public :
  //Constructor
  Matching_PO_Local(TreeGraph& , TreeGraph& ,NodeCost& ) ;
 //Distance Between Trees
  void distanceBetweenTree(int ,int ) ;
  //Distance Between Forests
  void distanceBetweenForest(int ,int ) ;
  //Distances Betweeen Forest Class
  void distanceBetweenForestClass(int ,int ) ;
  //Distances Betweeen Classes
  void distanceBetweenClass(int ,int ) ;
  
  DistanceType  match();

};


#endif

