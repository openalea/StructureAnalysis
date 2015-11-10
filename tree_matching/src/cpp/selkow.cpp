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


#include "selkow.h"

  // -------------
  // Constructeur
  // -------------
Selkow::Selkow(TreeGraph& input,TreeGraph& reference,NodeCost& nodeDistance)
{
  T1=&input;
  T2=&reference;
  ND=&nodeDistance;
  _distances.make(*T1,*T2,nodeDistance);
//  _insertCost.make(*T1,*T2,nodeDistance);
//  _deleteCost.make(*T1,*T2,nodeDistance);
//  _substitutionCost.make(*T1,*T2,nodeDistance);
//  _sumNbCaseVector = CaseVector(3,0);
//  _nbCaseVector = CaseVector(7,0);
  _distanceMatrix = DistanceVectorTable(T1->getNbVertex());
  _choices.resize(T1->getNbVertex(),T2->getNbVertex());
  // constante qui va permettre de calculer l'alignement restreint
  _restrMapp.link(I_MAX(T1->getDegree(),T2->getDegree()),*_distances.getDistanceTable());
}
  //--------------
  // Destructeur
  // -------------
Selkow::~Selkow()
{
  int size1 = T1->getNbVertex();
  int size2 = T2->getNbVertex();
//    for (int iv=0;iv<=size1-1;iv++)
//      {
//        delete (DistanceType*) _distanceMatrix[iv];
//      }
//    delete (DistanceType**) _distanceMatrix;
}

// ----------------------------------------------------------------------------------
// Calcule la distance entre les deux arbres T1[input_vertex] et T2[reference_vertex]
// ----------------------------------------------------------------------------------
DistanceType Selkow::distanceBetweenTree(int input_vertex,int reference_vertex)
{
  // On stocke dans ni et nj le nombre d'enfants de T1[i] et T2[j]
  int ni=T1->getNbChild(input_vertex);
  int nj=T2->getNbChild(reference_vertex);
  DistanceType MIN;
//  DistanceType insert,del,sub;
  int child_input,child_ref;
  child_input = T1->child(input_vertex,1);
  child_ref = T2->child(reference_vertex,1);
  MIN = getDBF(child_input,child_ref)+_distances.getCCost(input_vertex,reference_vertex);
  //cerr<<input_vertex<<" - "<<reference_vertex<<" c_cost = "<<_distances.getCCost(input_vertex,reference_vertex)<<endl;

  _choices.createList(input_vertex,reference_vertex);
  _choices.putFirst(input_vertex,reference_vertex,3);
  _choices.putLast(input_vertex,reference_vertex,reference_vertex);
//  insert=getInBF(child_input,child_ref);
//  del=getDeBF(child_input,child_ref);
//  sub=getSuBF(child_input,child_ref)+_distances.getCCost(input_vertex,reference_vertex);

  // On range dans le tableau des distances, la distance entre les arbres de racines input_vertex et
  // reference_vertex.
  //  _choices.putFirst(input_vertex,reference_vertex,MTC);
//  _insertCost.putDBT(input_vertex,reference_vertex,insert);
//  _deleteCost.putDBT(input_vertex,reference_vertex,del);
//  _substitutionCost.putDBT(input_vertex,reference_vertex,sub);
  
  _distances.putDBT(input_vertex,reference_vertex,MIN);
  return(MIN);
}



// -------------------------------
// Calcule la distance entre deux
// forets
// -------------------------------
DistanceType Selkow::distanceBetweenForest(int input_vertex,int reference_vertex)
{
// ni et nj representent le nombre de forets a comparees
//  cerr<<endl<<"Distance Between Forest : "<<input_vertex<<" - "<<reference_vertex<<" = ";
  int ni=T1->getNbChild(input_vertex);
  int nj=T2->getNbChild(reference_vertex);
  DistanceType cost1,cost2,cost3;
//DistanceType insert,del,sub;
  DistanceType DIST;
  int MFC=0;

  DIST=MAXDIST;

  int input_right_brother=T1->rightBrother(input_vertex);

  
  cost1=getDBF(input_right_brother,reference_vertex)+getDBT(input_vertex,EMPTY_TREE);
  if (cost1<=DIST) {
    DIST=cost1;
    MFC=1;
  }

  int reference_right_brother = T2->rightBrother(reference_vertex);
  
  cost2=getDBF(input_vertex,reference_right_brother)+getDBT(EMPTY_TREE,reference_vertex);

  if (cost2<=DIST) {
    DIST=cost2;
    MFC=2;
  }
   cost3=getDBF(input_right_brother,reference_right_brother)+getDBT(input_vertex,reference_vertex);

  if (cost3<=DIST) { 
    DIST=cost3; 
    MFC=3; 
  }
  //-------------------------------
  //We maintain the matching lists
  // On maintient les listes d'alignement
  //-------------------------------
//   _choices.createList(input_vertex,reference_vertex);
  _choices.putLast(input_vertex,reference_vertex,MFC);
  
  
  switch(MFC)
    {
    case 1 :
      {
	_choices.putLast(input_vertex,reference_vertex,input_right_brother);
//	del=getDBT(input_right_brother,EMPTY_TREE)+getDeBF(input_right_brother,reference_vertex);
//	insert=getInBF(input_right_brother,reference_vertex);
//	sub=getSuBF(input_right_brother,reference_vertex);
      }
      break;
    case 2 :
      {
	_choices.putLast(input_vertex,reference_vertex,reference_right_brother);
//	del=getDeBF(input_vertex,reference_right_brother);
//	insert=getDBT(EMPTY_TREE,reference_right_brother)+getInBF(input_vertex,reference_right_brother);
//	sub=getSuBF(input_vertex,reference_right_brother);
      }
      break;
    case 3 :
      {
	_choices.putLast(input_vertex,reference_vertex,reference_vertex);
//	sub=getSuBT(input_vertex,reference_vertex)+getSuBF(input_right_brother,reference_right_brother);
//	del=getDeBT(input_vertex,reference_vertex)+getDeBF(input_right_brother,reference_right_brother);
//	insert=getInBT(input_vertex,reference_vertex)+getInBF(input_right_brother,reference_right_brother);
      }
      break;
    default :   break;
    }
  
//  _deleteCost.putDBF(input_vertex,reference_vertex,del);
//  _insertCost.putDBF(input_vertex,reference_vertex,insert);
//  _substitutionCost.putDBF(input_vertex,reference_vertex,sub);
  _distances.putDBF(input_vertex,reference_vertex,DIST);
  return(DIST);
}


DistanceType Selkow::match()
{
  const int size1 = T1->getNbVertex();
  const int size2 = T2->getNbVertex();
  int pourcent = 0;
  DistanceType D=0;
   _distanceMatrix.resize(size1);
//  cerr << "Already computed : 0% matched ...    " << endl;
  if ((!T1->isNull())&&(!T2->isNull()))
    {
      for (int input_vertex=size1-1;input_vertex>=0;input_vertex--)
        {
            _distances.openDistancesVector(input_vertex);
//	    _insertCost.openDistancesVector(input_vertex);
// 	    _deleteCost.openDistancesVector(input_vertex);
// 	    _substitutionCost.openDistancesVector(input_vertex);
	   _distanceMatrix[input_vertex].resize(size2);
          for (int reference_vertex=size2-1;reference_vertex>=0;reference_vertex--)
            {
              DistanceType d =distanceBetweenTree(input_vertex,reference_vertex);
			  //cerr<<input_vertex<<" - "<<reference_vertex<<" = "<<d<<endl;
	      if ((input_vertex!=0)&&(reference_vertex!=0))
		distanceBetweenForest(input_vertex,reference_vertex);
	      (_distanceMatrix[input_vertex])[reference_vertex]=d;
            }
//		  if (int(100.- (100.*input_vertex)/size1)%10 == 0){
//			  pourcent += 10;
//	          cerr << "Already computed : "<< int(100.- (100.*input_vertex)/size1) <<"% "<<" matched ... " << endl;
//		  }
          for (int i=1;i<=T1->getNbChild(input_vertex);i++)
            {
	      _distances.closeDistancesVector(T1->child(input_vertex,i));
//	      _insertCost.closeDistancesVector(T1->child(input_vertex,i));
//	      _deleteCost.closeDistancesVector(T1->child(input_vertex,i));
//	      _substitutionCost.closeDistancesVector(T1->child(input_vertex,i));
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
  cerr<< "\x0d"<<endl;
  return(D);
}


void Selkow::TreeList(int input_vertex,int reference_vertex,Sequence& sequence)
{
  if ((!T1->isNull())&&(!T2->isNull()))
    {
      ChoiceList* L=_choices.getList(input_vertex,reference_vertex);
      int tree_choice=L->front();
      sequence.append(input_vertex,reference_vertex,_distances.getCCost(input_vertex,reference_vertex));
      int i_node=T1->child(input_vertex,1);
      int r_node=T2->child(reference_vertex,1);
      if ((r_node!=-1)&&(i_node!=-1)) ForestList(i_node,r_node,sequence);
    }
}

void Selkow::ForestList(int input_vertex,int reference_vertex,Sequence& sequence)
{
  ChoiceList* L=_choices.getList(input_vertex,reference_vertex);
  int forest_choice=Lat(L,2);
  switch(forest_choice)
    {
    case 1: ForestList(Lat(L,3),reference_vertex,sequence);break;
    case 2: ForestList(input_vertex,Lat(L,3),sequence);break;
    case 3: {
      TreeList(input_vertex,reference_vertex,sequence);
      int i_node=T1->rightBrother(input_vertex);
      int r_node=T2->rightBrother(reference_vertex);
      if ((r_node!=-1)&&(i_node!=-1)) ForestList(i_node,r_node,sequence);
    }break;
    default : break;
    }
}

DistanceType Selkow::getDBF(int input_vertex,int reference_vertex) 
{
  return(_distances.getDBOrderedF(input_vertex,reference_vertex));
}

DistanceType Selkow::getInBF(int input_vertex,int reference_vertex) 
{
  return(_insertCost.getDBOrderedF(input_vertex,reference_vertex));
}

DistanceType Selkow::getDeBF(int input_vertex,int reference_vertex) 
{
  return(_deleteCost.getDBOrderedF(input_vertex,reference_vertex));
}

DistanceType Selkow::getSuBF(int input_vertex,int reference_vertex) 
{
  return(_substitutionCost.getDBOrderedF(input_vertex,reference_vertex));
}
