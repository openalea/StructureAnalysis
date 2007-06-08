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


/*matching ne tenant pas compte du nombre de composantes connexes*/

#ifndef SB_SELFSIMILARMATCHING_HEADER
#define SB_SELFSIMILARMATCHING_HEADER

#include <vector>
#include "definitions.h"
#include "matchpath.h"
#include "choicetable.h"
#include "mdtable.h"
#include "treegraph.h"
#include "sequence.h"
#include "nodecost.h"




/**
 *\class Matching With Complex
 *\brief Algorithm for comparing two unordered tree graph
 *\par Presentation
 * In order to compute a distance between two multiscale tree graph, we have extended
 * the algorithm proposed by Zhang \cite{Zha93}. The distance is computed as the minimum
 * cost of sequence of edit operations needed to transform one tree graph into another one.
 *\par RequirementsSelfSimilar
 * - Two TreeGraphs defined at the same scale;
 * - A NodeCost (method for computing the local distance),
 *\author Pascal ferraro
 *\date 1999
 */


class SelfSimilarMatching
{

  public :
    //Constructor
    SelfSimilarMatching() {};
    SelfSimilarMatching(TreeGraph& , TreeGraph& ,NodeCost& ) ;
    SelfSimilarMatching(TreeGraph& , TreeGraph& ,NodeCost&, int& ) ;
    void make(TreeGraph& , TreeGraph& ,NodeCost& ) ;
    //Destructor
    ~SelfSimilarMatching();

    //Distance Between Trees
    virtual DistanceType distanceBetweenTree(int ,int ) ;
    //Ints Betweeen Forest
    virtual DistanceType distanceBetweenForest(int ,int ) ;



    //Operator
    DistanceType getDBT(int ,int ) const;
    DistanceType getDBF(int ,int ) const;
    DistanceType  match();
    DistanceType getDistanceMatrix(int ,int) const;
    int operator()(int );
    int getNbMin(){return _nb_min;};
    
  protected :
    TreeGraph* T1;
    TreeGraph* T2;
    MatchingDistanceTable _distances;
    NodeCost* ND;
    MatchPath _restrMapp;
    VertexVector _restrMappList;
    DistanceVectorTable _distanceMatrix;
    int _nb_min;
};


#endif

