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


#ifndef SB_MATCHING_UMINCO_HEADER
#define SB_MATCHING_UMINCO_HEADER

#include "matching_unordered.h"


/**
 *\class Matching_UMINCO
 *\brief Algorithm for comparing two unordered tree graph minimizing the number of connected components.
 *\par Presentation
 * Computes an optimal valid mapping which defines a minimum of  connected components
 *\par Requirements
 * - Two TreeGraphs defined at the same scale;
 * - A NodeCost (method for computing the local distance), 
 *\author Pascal ferraro
 *\date 1999
 */

class Matching_U_MinCo : public Matching_U
{
	
  public : 
  Matching_U_MinCo(){};
    Matching_U_MinCo(TreeGraph& , TreeGraph& ,NodeCost& ) ;
    void make(TreeGraph& , TreeGraph& ,NodeCost& ) ;
    //Destructor
    ~Matching_U_MinCo();
    //Distance Between Trees
    DistanceType distanceBetweenTree(int ,int ); 
    //Distances Betweeen Forest
    DistanceType distanceBetweenForest(int ,int );
    //Operator
    DistanceType  match();
    int getNBC(int,int) const;
    int getNBCRef(int,int) const;

  protected :
    MatchingConnectedTable _nbInputTreeConnectedComponents;
    MatchingConnectedTable _nbReferenceTreeConnectedComponents;
    MatchingConnectedTable _nbInputForestConnectedComponents;
    MatchingConnectedTable _nbReferenceForestConnectedComponents;
    MatchingConnectedTable _InputRootTable;
    MatchingConnectedTable _nbInputRootMapped;
    MatchingConnectedTable _ReferenceRootTable;
    MatchingConnectedTable _nbReferenceRootMapped;
};


#endif

