/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       AMAPmod: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2000 UMR Cirad/Inra Modelisation des Plantes
 *
 *       File author(s): P.ferraro (pascal.ferraro@labri.fr)
 *
 *       $Source: /usr/cvsmaster/AMAPmod/src/TreeMatching/Ordered/matching_O.cpp,v $
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


#include "matching_O.h"

// -------------
// Constructeur
// -------------

Matching_O::Matching_O(TreeGraph& input,TreeGraph& reference,NodeCost& nodeDistance)
{
  TreeNode* node =  (&input)->getNode(0);
  int scale = (&input)->getMTG()->vscale(node->getVertex());
  _msm = new MS_O_Matching(input, reference ,nodeDistance,scale);
}

// -------------
// Destructeur
// -------------
Matching_O::~Matching_O()
{
  delete (MS_O_Matching*) _msm;
}



//---------------------------------------------------------------------------------------
// Calcul de la distance entre les deux arbres T1 et T2
//---------------------------------------------------------------------------------------
DistanceType  Matching_O::run()
{

  return _msm->match();

}

 DistanceType Matching_O::getInsertionCost(){return _msm->getInsertionCost();}
  DistanceType Matching_O::getDeletionCost(){return _msm->getDeletionCost();}
  DistanceType Matching_O::getMatchingCost(){return _msm->getMatchingCost();}

 Sequence* Matching_O::getSequence(){return _msm->getSequence();} 

void Matching_O::TreeList(int i,int j,Sequence& s){
  Sequence* align = getSequence();
  if (align->getSize())
    {
      align->reset();
      do
	{
	  s.append(align->getCurrent()->getIV(),align->getCurrent()->getRV(),align->getCurrent()->getCost());
	  //cout<<align->getCurrent()->getIV()<<" "<<align->getCurrent()->getRV()<<endl;
	} while(align->next());
      
    } 
}

DistanceType Matching_O::getDBT(int i, int j){
  TreeNode* _node1;
  TreeNode* _node2;
  int n1;
  int n2;
  int nbV1 = _msm->T1->getNbVertex();
  int nbV2 = _msm->T2->getNbVertex();
  
  if(i!=-1){
    _node1 =  _msm->T1->getNode(i);
    n1 = _node1->getNumPostfix();
  } 
  if(j!=-1){ 
    _node2 =  _msm->T2->getNode(j);
    n2 = _node2->getNumPostfix();
  }
  if(i!=-1 && j==-1){
    return  _msm->fDist[_msm->fix(nbV1-1,nbV2-1)][n1][0];
   }
  if(i==-1 && j!=-1){
    return  _msm->fDist[_msm->fix(nbV1-1,nbV2-1)][0][n2];
  }
  
  if(i!=-1 && j!=-1){
    return _msm->getTreeDistance(n1, n2);
  }
   if(i==-1 && j==-1){
     return 0;
   }
   
}
