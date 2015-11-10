/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       AMAPmod: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2000 UMR Cirad/Inra Modelisation des Plantes
 *
 *       File author(s): P.ferraro (pascal.ferraro@cirad.fr)
 *
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


#include "local_matching_u.h"


  // -------------
  // Constructeur
  // -------------
LocalMatching_U::LocalMatching_U(TreeGraph& input,TreeGraph& reference,NodeCost& nodeDistance):
  Matching_U(input,reference,nodeDistance)
{
  _localrestrMapp.link(I_MAX(T1->getDegree(),T2->getDegree()),*_distances.getDistanceTable());
  _maxScore = 0.;
  _maxInput = 0;
  _maxRef   = 0;
  _init = 1;
}

// -------------
// Destructeur
// -------------
LocalMatching_U::~LocalMatching_U()
{
  int size1 = T1->getNbVertex();
  int size2 = T2->getNbVertex();
}

// ----------------------------------------------------------------------------------
// Calcule la distance entre les deux arbres T1[input_vertex] et T2[reference_vertex]
// ----------------------------------------------------------------------------------
DistanceType LocalMatching_U::distanceBetweenTree(int input_vertex,int reference_vertex)
{
  // On stocke dans ni et nj le nombre d'enfants de T1[i] et T2[j]
  int ni=T1->getNbChild(input_vertex);
  int nj=T2->getNbChild(reference_vertex);

  DistanceType cost1,cost2,cost3,dist1,dist2;
  DistanceType max,MAX=0.;
  int im=0,jm=0,MTC=0;
  int i;
  
  //cout<<"Distance Between subtrees "<<input_vertex<<" - "<<reference_vertex<<endl;

  //----------------------------------------------------------------------
  //Case 1 : We search the reference_tree as a subtree of the input_tree
  //----------------------------------------------------------------------
  max=0.;
  for (i=1;i<=ni;i++)
    {
      // On cherche parmi tous les fils de input_vertex celui dont l'arbre est le plus ressemblant a T2
      int input_child=T1->child(input_vertex,i);
      dist1=getDBT(input_child,reference_vertex)+_distances.getDCost(input_child);
	  // On conserve la plus grande distance
      if (dist1>max) 
	{ 
	  max=dist1; 
	  im=input_child; 
	}
    }
  cost1=max;
  // On conserve le cout maximum
  if (cost1>MAX) 
    { 
      MAX=cost1; 
      MTC=1; 
    }

  //--------------------------------------------------------------------
  //Case2 : We search the input_tree as a subtree of the reference_tree
  //--------------------------------------------------------------------
  max=0;
  for (i=1;i<=nj;i++)
    {
      // On recherche parmi tous les fils de T2 celui qui ressemble le plus a T1
      int reference_child=T2->child(reference_vertex,i);
      dist2=getDBT(input_vertex,reference_child)+_distances.getICost(reference_child);
	  if (dist2>max) 
	{ 
	  max=dist2;
	  jm=reference_child; 
	}
    }
  cost2=max;
  // On conserve le cout s'il est superieur au precedent
  if (cost2>MAX) 
    { 
      MAX=cost2; 
      MTC=2; 
    }

  //----------------------------------------------------------------------------------
  //Case3 : We evaluate the matching between the input_forest and the reference_forest
  //----------------------------------------------------------------------------------
  // Le cout est celui de l'alignement des deux forets
  // plus celui de l'echange de T1(i) en T2(j)
  cost3=getDBF(input_vertex,reference_vertex);
  cost3=cost3+_distances.getCCost(input_vertex,reference_vertex);
  //cout<<"cost 3 = "<<cost3<<endl;
  if (cost3>MAX) 
    {
      MAX=cost3; 
      MTC=3;  
    }

  //-----------------------------------
  // We maintain the matching lists
  //-----------------------------------
  switch (MTC)
    {
    case 1 :{
      _choices.putFirst(input_vertex,reference_vertex,im);
      _choices.putLast(input_vertex,reference_vertex,EMPTY_NODE);
    }break;
    case 2 :{
      _choices.putFirst(input_vertex,reference_vertex,jm);
      _choices.putLast(input_vertex,reference_vertex,M(input_vertex,jm));
    }break;
    case 3 :{
      _choices.putFirst(input_vertex,reference_vertex,EMPTY_NODE);
      _choices.putLast(input_vertex,reference_vertex,reference_vertex);
    }break;
    default :   assert(0);break;
    }
  _choices.putFirst(input_vertex,reference_vertex,MTC);
  // On range dans le tableau des distances, la distance entre les arbres de racines input_vertex et
  // reference_vertex.
  _distances.putDBT(input_vertex,reference_vertex,MAX);
  // We remind the maximum score
  if (MAX>_maxScore)
    {
      _maxScore = MAX;
      _maxInput = input_vertex;
      _maxRef   = reference_vertex;
    }

  return(MAX);

}



// -------------------------------
// Calcule la distance entre deux
// forets
// -------------------------------
DistanceType LocalMatching_U::distanceBetweenForest(int input_vertex,int reference_vertex)
{
// ni et nj representent le nombre de forets a comparees
  int ni=T1->getNbChild(input_vertex);
  int nj=T2->getNbChild(reference_vertex);
  DistanceType cost1,cost2,cost3,dist1,dist2;
  DistanceType max,DIST;
  int im=0,jm=0,MFC=0;
  int i;

  DIST=0.;
  //  cout<<"Distance Between subforests "<<input_vertex<<" - "<<reference_vertex<<endl;

  //------------------------------------------------------------------------
  //Case 1 : We search the reference_forest as a subtree of the input_forest
  // On met en correspondance une sous-foret d'un arbre de F1 avec la foret F2
  //------------------------------------------------------------------------
  max=0.;
  for (i=1;i<=ni;i++)
    {
      int input_child=T1->child(input_vertex,i);
      dist1=getDBF(input_child,reference_vertex)+_distances.getDCost(input_child);
      if (dist1>max) 
	{ 
	  max=dist1;
	  im=input_child;
	}
    }
  cost1=max;


  if (cost1>=DIST) 
    {
      DIST=cost1;
      MFC=1;
    }

  //------------------------------------------------------------------------
  //Case 2 : We search the input_forest as a subtree of the reference_forest
  // On met en correspondance une sous-foret d'un arbre de F2 avec la foret F1
  //------------------------------------------------------------------------
  max=0.;
  for (i=1;i<=nj;i++)
    {
      int reference_child=T2->child(reference_vertex,i);
      dist2=getDBF(input_vertex,reference_child)+_distances.getICost(reference_child);
      if (dist2>max) 
	{ 
	  max=dist2;
	  jm=reference_child;
	}
    }
  cost2=max;

  
  if (cost2>=DIST) 
    {
      DIST=cost2;
      MFC=2;
    }

  //---------------------------------------------------------------------------------------------
  //Case 3 : We evaluate the restricted mapping between the input_forest and the reference_forest
  // On evalue l'alignement restreint entre les deux forets
  //---------------------------------------------------------------------------------------------

  // On fabrique le graphe de flot necessaire a la resolution du probleme
  NodeList* input_list=new NodeList();
  NodeList* reference_list=new NodeList();
  for (int s1=1;s1<=ni;s1++) 
    { 
      input_list->push_back(T1->child(input_vertex,s1)); 
    }
  for (int s2=1;s2<=nj;s2++) 
    { 
      reference_list->push_back(T2->child(reference_vertex,s2));
    }
  _localrestrMapp.make(*input_list,*reference_list);
  _restrMappList.resize(ni+nj+3,EMPTY_NODE);

  // THE INPUT FOREST IS EMPTY_TREE
  // All the reference vertices are paired with empty
  // Si la foret initiale est vide, il faut inserer toutes les arbres de
  // la foret de reference et tous les noeuds de references sont associes
  // avec le noeud vide
  if (ni==0)
    {
      _restrMappList[1]=2;
      for (i=1;i<=nj;i++) 
	{ 
	  _restrMappList[i+1]=1; 
	}
      //cost3=getDBF(EMPTY_TREE,reference_vertex);
      //score to EMPTY Tree is null
      cost3 = 0.;
    }
  else
    {
      // THE REFERENCE FOREST IS EMPTY_TREE
      // All the input vertices are paired with empty
      // Si c'est l'arbre de reference qui est vide,
      // il faut supprimer la foret initiale et tous les
      // noeuds de cette foret seront associer avec un
      // noeud vide
      if (nj==0)
        {
          _restrMappList[2]=1;
          for (i=1;i<=ni;i++) 
	    { 
	      _restrMappList[i]=ni+1; 
	    }
          //cost3=getDBF(input_vertex,EMPTY_TREE);
	  cost3 = 0.;
        }
      else
        {
          //BOTH FOREST ARE NOT EMPTY_TREE
          // A retricted mapping must be calculated
          // Sinon on resout le probleme de flot maximum de cout minimum
          cost3=_localrestrMapp.maxCostFlow(_restrMappList);
	  //cout<<"cost 3 = "<<cost3<<endl;

        }
    }

 

  if (cost3>=DIST) 
    { 
      DIST=cost3; 
      MFC=3; 
    }

  //-------------------------------
  //We maintain the matching lists
  // On maintient les listes d'alignement
  //-------------------------------
  _choices.createList(input_vertex,reference_vertex);
  _choices.putFirst(input_vertex,reference_vertex,MFC);


  switch(MFC)
    {
    case 1 :
      {
        _choices.putLast(input_vertex,reference_vertex,im);
      }
      break;
    case 2 :
      {
        _choices.putLast(input_vertex,reference_vertex,jm);
      }
      break;
      case 3 :
        {
	  for (int i=1;i<=T1->getNbChild(input_vertex);i++)
	    {
	      _choices.putLast(input_vertex,reference_vertex,_localrestrMapp.who(_restrMappList[i]));
	    }
        }
        break;
    default :   break;
    }


  delete (NodeList*) input_list;
  delete (NodeList*) reference_list;
  _distances.putDBF(input_vertex,reference_vertex,DIST);
  return(DIST);

}



int LocalMatching_U::getMaxInput()
{
  return (_maxInput);
}

int LocalMatching_U::getMaxRef()
{
  return (_maxRef);
}

void LocalMatching_U::TreeList(int input_vertex,int reference_vertex,Sequence& sequence)
{
  if ((!T1->isNull())&&(!T2->isNull()))
    {
      if (_init) {
	_init = 0;
	TreeList(_maxInput,_maxRef,sequence);
      }
      else {
	ChoiceList* L=_choices.getList(input_vertex,reference_vertex);
	int tree_choice=L->front();
	switch(tree_choice)
	  {
	  case 1:
	    {
	      TreeList(Lat(L,1),reference_vertex,sequence);
	    }
	    break;
	  case 2:
	    {
            TreeList(input_vertex,Lat(L,1),sequence);
	    }
	    break;
	  case 3: {
	    sequence.append(input_vertex,reference_vertex,_distances.getCCost(input_vertex,reference_vertex));
	    ForestList(input_vertex,reference_vertex,sequence);
	  }break;
	  default : break;
	  }
      }
    }
}


DistanceType LocalMatching_U::match()
{
  const int size1 = T1->getNbVertex();
  const int size2 = T2->getNbVertex();
  DistanceType D=0;
  if ((!T1->isNull())&&(!T2->isNull()))
    {
      for (int input_vertex=size1-1;input_vertex>=0;input_vertex--)
        {
          _distances.openDistancesVector(input_vertex);
          for (int reference_vertex=size2-1;reference_vertex>=0;reference_vertex--)
            {
              distanceBetweenForest(input_vertex,reference_vertex);
              DistanceType d =distanceBetweenTree(input_vertex,reference_vertex);
            }
	  if (int(100. - 100*input_vertex/size1)%10 == 0)
	    cerr << "\x0d" << "Already computed : "<<int(100. - 100*input_vertex/size1) <<"% " <<" matched ... " << flush;
          for (int i=1;i<=T1->getNbChild(input_vertex);i++)
            {
              _distances.closeDistancesVector(T1->child(input_vertex,i));
            }
        }
      D=getDBT(_maxInput,_maxRef);
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
DistanceType LocalMatching_U::getDistanceMatrix(int iv,int rv) const
{
  return(_distanceMatrix[iv][rv]);
}

DistanceType LocalMatching_U::getDBT(int input_vertex,int reference_vertex) const
{
  if ((input_vertex == -1) || (reference_vertex==-1))
    return (0.);
  else
    return(_distances.getDBT(input_vertex,reference_vertex));
}

DistanceType LocalMatching_U::getDBF(int input_vertex,int reference_vertex) const
{
  if ((input_vertex == -1) || (reference_vertex==-1))
    return (0.);
  else
    return(_distances.getDBF(input_vertex,reference_vertex));
}









