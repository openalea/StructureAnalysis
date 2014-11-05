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


#ifndef SB_EXT_MATCHING_HEADER
#define SB_EXT_MATCHING_HEADER

#include "definitions.h"
#include "matchpath.h"
#include "emchoicetable.h"
#include "emdtable.h"
#include "treegraph.h"
#include "sequence.h"
#include "nodecost.h"
#include "wnodecost.h"
#include "mnodecost.h"


class ExtendedMatching
{

public :
  //Constructor
  ExtendedMatching() {};
  ExtendedMatching(TreeGraph& , TreeGraph& ,NodeCost& ) ;
  void make(TreeGraph& , TreeGraph& ,NodeCost& ) ;
  //Destructor
  ~ExtendedMatching();
  //Distance Between Trees
  DistanceType distanceBetweenTree(int ,int );
  //Distances Betweeen Forest
  DistanceType distanceBetweenForest(int ,int );
  //Operator
  DistanceType getDBT(int ,int ) const;
  DistanceType getDBF(int ,int ) const;
  DistanceType getDBMF(int ,int ) const;
  DistanceType getDBPF(int ,int ) const;
  DistanceType getDBIMFRF(int ,int ) const;
  DistanceType getDBIFRMF(int ,int ) const;
  DistanceType match();
  void getList(int ,int ,Sequence*);
  void ForestList(int ,int ,Sequence& );
  void pfList(int ,int ,Sequence& );
  void mfList(int ,int ,Sequence& );
  void ifrmfList(int ,int ,Sequence& );
  void imfrfList(int ,int ,Sequence& );
  void TreeList(int ,int ,Sequence& );
  int operator()(int );


private :
  TreeGraph* T1;
  TreeGraph* T2;
  ExtendedMatchingDistanceTable _distances;
  ExtendedMatchingChoiceTable _choices;
  NodeCost* ND;
  MatchPath _restrMappIFAndRMF;
  MatchPath _restrMappIMFAndRF;
  MatchPath _restrMappIPFAndRPF;
  MatchPath _restrMappIMFAndRMF;
  VertexVector _restrMappIFAndRMFList;
  VertexVector _restrMappIMFAndRFList;
  VertexVector _restrMappIPFAndRPFList;
  VertexVector _restrMappIMFAndRMFList;
  int M(int,int);
  AmlBoolean sameComplex(TreeNode* ,TreeNode* );
  int Lat(ChoiceList* L, int vertex);
};


#endif

