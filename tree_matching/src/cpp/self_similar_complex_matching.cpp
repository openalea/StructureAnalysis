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


#include "self_similar_complex_matching.h"
#include "dec_matrix.h"


  // -------------
  // Constructeur
  // -------------
SelfSimilarComplexMatching::SelfSimilarComplexMatching(TreeGraph& input,TreeGraph& reference,NodeCost& nodeDistance)
{
  T1=&input;
  T2=&reference;
  ND=&nodeDistance;
  _distances.make(*T1,*T2,nodeDistance);
  _d_l_w.make(*T1,*T2,nodeDistance);
  _d_v_l.make(*T1,*T2,nodeDistance);
  _d_l_l.make(*T1,*T2,nodeDistance);
  _d_v_w.make(*T1,*T2,nodeDistance);
  _restrDistances_v_w.make(*T1,*T2,nodeDistance);
  _restrDistances_l_l.make(*T1,*T2,nodeDistance);
  _restrDistances.make(*T1,*T2,nodeDistance);
 // _choices est un tableau de listes retenant les tentatives successives alignements durant l'algo.
  // c'est donc un tableau de |T1| lignes et |T2| colonnes initialise a 0
  _choices.resize(T1->getNbVertex(),T2->getNbVertex());
  _choices_v_l.resize(T1->getNbVertex(),T2->getNbVertex());
  _choices_l_w.resize(T1->getNbVertex(),T2->getNbVertex());
  _choices_v_w.resize(T1->getNbVertex(),T2->getNbVertex());
  _choices_l_l.resize(T1->getNbVertex(),T2->getNbVertex());
  // constante qui va permettre de calculer l'alignement restreint
  _distanceMatrix = DistanceVectorTable(T1->getNbVertex());
  _restrMapp_v_w.link(I_MAX(T1->getDegree(),T2->getDegree()),*_restrDistances_v_w.getDistanceTable());
  _restrMapp_l_l.link(I_MAX(T1->getDegree(),T2->getDegree()),*_restrDistances_l_l.getDistanceTable());
}



DistanceType SelfSimilarComplexMatching::match()
{
  
  const int size1 = T1->getNbVertex();
  const int size2 = T2->getNbVertex();
  DistanceType D=0;
  if ((!T1->isNull())&&(!T2->isNull()))
    {
      for (int input_vertex=T1->getNbVertex()-1;input_vertex>=0;input_vertex--)
	{
	  _distances.openDistancesVector(input_vertex);
	  _distanceMatrix[input_vertex].resize(size2);
	  _d_v_l.openDistancesVector(input_vertex);
	  _d_l_w.openDistancesVector(input_vertex);
	  _d_v_w.openDistancesVector(input_vertex);
	  _d_l_l.openDistancesVector(input_vertex);
	  _restrDistances_v_w.openDistancesVector(input_vertex);
	  _restrDistances_l_l.openDistancesVector(input_vertex);
	  _restrDistances.openDistancesVector(input_vertex);
	  
	  for (int reference_vertex=T2->getNbVertex()-1;reference_vertex>=0;reference_vertex--)
	    {
		  //std::cerr<<"Distance Between Forest "<<input_vertex<<"  -  "<<reference_vertex<<endl;
	      distanceBetweenForest(input_vertex,reference_vertex);
	      DistanceType d = distanceBetweenTree(input_vertex,reference_vertex);
	      (_distanceMatrix[input_vertex])[reference_vertex]=d;
	    }
	  if (int(100. - 100*input_vertex/size1)%10 == 0)
		  std::cerr << "\x0d" << "Already computed : "<<int(100. - 100*input_vertex/size1)<<" %"<< endl;
	  
	  for (int i=1;i<=T1->getNbChild(input_vertex);i++)
	    {
	      _distances.closeDistancesVector(T1->child(input_vertex,i));
	      _d_v_l.closeDistancesVector(T1->child(input_vertex,i));
	      _d_l_w.closeDistancesVector(T1->child(input_vertex,i));
	      _d_v_w.closeDistancesVector(T1->child(input_vertex,i));
	      _d_l_l.closeDistancesVector(T1->child(input_vertex,i));
	      _restrDistances_v_w.closeDistancesVector(T1->child(input_vertex,i));
	      _restrDistances_l_l.closeDistancesVector(T1->child(input_vertex,i));
	      _restrDistances.closeDistancesVector(T1->child(input_vertex,i));
	    }
	}

      D=getDBT(0,0);
      
    }
  else
    {
      if (T1->isNull())
	{
	  if (!T2->isNull()) {D=_distances.referenceTreeFromEmpty(0);}
	}
      else
	{
	  D=_distances.inputTreeToEmpty(0);
	}
    }
  return(D);
}
// --------------------------------------------
// Renvoie les distances entre arbres et forets
// --------------------------------------------
DistanceType SelfSimilarComplexMatching::getDistanceMatrix(int iv,int rv) 
{
  return((_distanceMatrix[iv])[rv]);
}

