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


#ifndef SB_MATCHING4_HEADER
#define SB_MATCHING4_HEADER

#include "definitions.h" 
#include "matchpath.h" 
#include "choicetable.h" 
#include "mdtable.h"
#include "treegraph.h"
#include "sequence.h" 
#include "nodecost.h"
#include "wnodecost.h"
#include "mnodecost.h"
#include "mctable.h"
#include "matching.h"


/**
 *\class MatchingMinimizeComponents
 *\brief Algorithm for comparing two unordered tree graph 
 *\par Presentation
 * Computes an optimal valid mapping which defines a minimum of  connected components
 *\par Requirements
 * - Two TreeGraphs defined at the same scale;
 * - A NodeCost (method for computing the local distance), 
 *\author Pascal ferraro
 *\date 1999
 */

class MatchingMinimizeComponents : public Matching
{
	
  public :
    MatchingMinimizeComponents(TreeGraph& , TreeGraph& ,NodeCost& ) ;
    void make(TreeGraph& , TreeGraph& ,NodeCost& ) ;
    //Destructor
    ~MatchingMinimizeComponents();
    //Distance Between Trees
    DistanceType distanceBetweenTree(int ,int ); 
    //Distances Betweeen Forest
    DistanceType distanceBetweenForest(int ,int );
    //Operator
    DistanceType  match();
    int getNBC(int,int) const;
    int getNBCRef(int,int) const;


    int getNBCT1(int,int) const;
    int getNBCRefT2(int,int) const;


  private :
    MatchingConnectedTable _nbTreeConnectedComponents;
    MatchingConnectedTable _nbInputConnectedComponents;
    MatchingConnectedTable _nbReferenceConnectedComponents;
    MatchingConnectedTable _InputRootTable;
    MatchingConnectedTable _nbInputRootMapped;
    MatchingConnectedTable _ReferenceRootTable;
    MatchingConnectedTable _nbReferenceRootMapped;


    MatchingConnectedTable _nbInputConnectedComponentsT1;
    MatchingConnectedTable _nbReferenceConnectedComponentsT2;
    MatchingConnectedTable _InputRootTableT1;
    MatchingConnectedTable _nbInputRootMappedT1;
    MatchingConnectedTable _ReferenceRootTableT2;
    MatchingConnectedTable _nbReferenceRootMappedT2;



    ChoiceTable _choicesT1;
    ChoiceTable _choicesT2;
 

};


#endif

