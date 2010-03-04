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


#include "matching_unordered.h"

  // -------------
  // Constructeur
  // -------------
Matching_U::Matching_U(TreeGraph& input,TreeGraph& reference,NodeCost& nodeDistance):
  T1(&input),
  T2(&reference),
  ND(&nodeDistance),
  _spaceOptimization(0),
  _sumNbCaseVector(3,0),
  _nbCaseVector(7,0),
  _distanceMatrix(T1->getNbVertex())
{
  _i_v = 0;
  _r_v = 0;
  _distances.make(*T1,*T2,nodeDistance);
  _insertCost.make(*T1,*T2,nodeDistance);
  _deleteCost.make(*T1,*T2,nodeDistance);
  _substitutionCost.make(*T1,*T2,nodeDistance);
  // _choices est un tableau de listes retenant les tentatives successives alignements durant l'algo.
  // c'est donc un tableau de |T1| lignes et |T2| colonnes initialise a 0
  _choices.resize(T1->getNbVertex(),T2->getNbVertex());
  // constante qui va permettre de calculer l'alignement restreint
  _restrMapp.link(I_MAX(T1->getDegree(),T2->getDegree()),*_distances.getDistanceTable());
}
// -------------
Matching_U::Matching_U(TreeGraph& input,TreeGraph& reference,NodeCost& nodeDistance,int& spaceOpti):
  T1(&input),
  T2(&reference),
  ND(&nodeDistance),
  _spaceOptimization(1),
  _distanceMatrix(T1->getNbVertex())
{
  _distances.make(*T1,*T2,nodeDistance);
  _restrMapp.link(I_MAX(T1->getDegree(),T2->getDegree()),*_distances.getDistanceTable());
}
// Destructeur
  // -------------
Matching_U::~Matching_U()
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
DistanceType Matching_U::distanceBetweenTree(int input_vertex,int reference_vertex)
{
  // On stocke dans ni et nj le nombre d'enfants de T1[i] et T2[j]
  int ni=T1->getNbChild(input_vertex);
  int nj=T2->getNbChild(reference_vertex);

  DistanceType cost1,cost2,cost3,dist1,dist2,insert,del,sub;
  DistanceType min,MIN=2*MAXDIST;
  int im=0,jm=0,MTC=0;
  int i;

  //----------------------------------------------------------------------
  //Case 1 : We search the reference_tree as a subtree of the input_tree
  //         On cherche a mettre en correspondance l'arbre de reference
  //         avec un sous arbre de l'arbre initial, il faut donc effacer
  //         T1 moins le sous arbre qui ressemble le plus a T2
  //----------------------------------------------------------------------
  min=MAXDIST;
  // cout de l'effacement de l'arbre initial
  cost1=getDBT(input_vertex,EMPTY_TREE);
  for (i=1;i<=ni;i++)
    {
      // On cherche parmi tous les fils de input_vertex celui dont l'arbre est le plus ressemblant a T2
      int input_child=T1->child(input_vertex,i);
      // la distance est donc le passage de T1[iam] en T2[j] - l'effacement de T[iam] qui a ete
      // compte precedemment
      dist1=getDBT(input_child,reference_vertex)-getDBT(input_child,EMPTY_TREE);
      // On conserve la plus petite distance
      if (dist1<min) { min=dist1; im=input_child; }
    }
  cost1=cost1+min;
  // On conserve le cout minimum
  if (cost1<MIN) { MIN=cost1; MTC=1; }

  //--------------------------------------------------------------------
  //Case2 : We search the input_tree as a subtree of the reference_tree
  //        On cherche a mettre en correspondance T1 et un sous arbre de
  //        T2, il faut donc inserer T2 dans T1 moins l'arbre qui
  //        ressemble le plus a T2 qu'on transforme
  //--------------------------------------------------------------------
  min=MAXDIST;
  // cout de l'insertion de l'arbre T2
  cost2=getDBT(EMPTY_TREE,reference_vertex);
  for (i=1;i<=nj;i++)
    {
      // On recherche parmi tous les fils de T2 celui qui ressemble le plus a T1
      int reference_child=T2->child(reference_vertex,i);
      dist2=getDBT(input_vertex,reference_child)-getDBT(EMPTY_TREE,reference_child);
      if (dist2<min) { min=dist2;jm=reference_child; };
    }
  cost2=cost2+min;
  // On conserve le cout s'il est inferieur au precedent
  if (cost2<MIN) { MIN=cost2; MTC=2; }

  //----------------------------------------------------------------------------------
  //Case3 : We evaluate the matching between the input_forest and the reference_forest
  // On evalue la mise en correspondance des arbres des deux forets issues de T1 et T2
  //----------------------------------------------------------------------------------
  // Le cout est celui de l'alignement des deux forets
  // plus celui de l'echange de T1(i) en T2(j)
  cost3=getDBF(input_vertex,reference_vertex);
  cost3=cost3+_distances.getCCost(input_vertex,reference_vertex);
  
  
  // On conserve le cout s'il est inferieur au precedent
  if (cost3<MIN) { MIN=cost3; MTC=3;  }

  //-----------------------------------
  // We maintain the matching lists
  // mise a jour des listes d'alignement
  //-----------------------------------
  if (!_spaceOptimization){
    switch (MTC)
      {
      case 1 :{
	_choices.putFirst(input_vertex,reference_vertex,im);
	_choices.putLast(input_vertex,reference_vertex,-1);
	insert=getInBT(im,reference_vertex);
	del=getDBT(input_vertex,EMPTY_TREE)-getDBT(im,EMPTY_TREE)+getDeBT(im,reference_vertex);
	sub=getSuBT(im,reference_vertex);
      }break;
      case 2 :{
	_choices.putFirst(input_vertex,reference_vertex,jm);
	_choices.putLast(input_vertex,reference_vertex,M(input_vertex,jm));
	insert=getDBT(EMPTY_TREE,reference_vertex)-getDBT(EMPTY_TREE,jm)+getInBT(input_vertex,jm);
	del=getDeBT(input_vertex,jm);
	sub=getSuBT(input_vertex,jm);
      }break;
      case 3 :{
	_choices.putFirst(input_vertex,reference_vertex,-1);
	_choices.putLast(input_vertex,reference_vertex,reference_vertex);
	insert=getInBF(input_vertex,reference_vertex);
	del=getDeBF(input_vertex,reference_vertex);
	sub=getSuBF(input_vertex,reference_vertex)+_distances.getCCost(input_vertex,reference_vertex);
      }break;
      default :   assert(0);break;
      }
    _choices.putFirst(input_vertex,reference_vertex,MTC);
    // On range dans le tableau des distances, la distance entre les arbres de racines input_vertex et
    // reference_vertex.
    _insertCost.putDBT(input_vertex,reference_vertex,insert);
    _deleteCost.putDBT(input_vertex,reference_vertex,del);
    _substitutionCost.putDBT(input_vertex,reference_vertex,sub);
  
    
    // On calcule les differents cas rencontres
    
    if ((cost1==cost2) && (cost1==cost3))
      { _nbCaseVector[0]++;}
    if ((cost1==cost2) && (cost1<cost3))
      { _nbCaseVector[1]++;}
    if ((cost1==cost3) && (cost1<cost2))
      { _nbCaseVector[2]++;}
    if ((cost2==cost3) && (cost2<cost1))
      { _nbCaseVector[3]++;}
    if ((cost1<cost3) && (cost1<cost2))
      { _nbCaseVector[4]++;}
    if ((cost2<cost3) && (cost2<cost1))
      { _nbCaseVector[5]++;}
    if ((cost3<cost1) && (cost3<cost2))
      { _nbCaseVector[6]++;}
    
    _sumNbCaseVector[MTC-1]++;
  }
  _distances.putDBT(input_vertex,reference_vertex,MIN);
  //  cerr<<input_vertex<<" - "<<reference_vertex<<" => "<<"cost1 = "<<cost1<<" - cost 2 = "<<cost2<<" - cost 3 = "<<cost3<<endl;
  return(MIN);
}



// -------------------------------
// Calcule la distance entre deux
// forets
// -------------------------------
DistanceType Matching_U::distanceBetweenForest(int input_vertex,int reference_vertex)
{
// ni et nj representent le nombre de forets a comparees
  int ni=T1->getNbChild(input_vertex);
  int nj=T2->getNbChild(reference_vertex);
  DistanceType cost1,cost2,cost3,dist1,dist2,insert,del,sub;
  DistanceType min,DIST;
  int im=0,jm=0,MFC=0;
  int i;

  DIST=MAXDIST;

  //------------------------------------------------------------------------
  //Case 1 : We search the reference_forest as a subtree of the input_forest
  // On met en correspondance une sous-foret d'un arbre de F1 avec la foret F2
  //------------------------------------------------------------------------
  min=MAXDIST;
  cost1=getDBF(input_vertex,EMPTY_TREE);
  for (i=1;i<=ni;i++)
    {
      int input_child=T1->child(input_vertex,i);
      dist1=getDBF(input_child,reference_vertex)-getDBF(input_child,EMPTY_TREE);
      if (dist1<min) { min=dist1;im=input_child;}
    }
  cost1=cost1+min;

  if (cost1<=DIST) {DIST=cost1;MFC=1;}

  //------------------------------------------------------------------------
  //Case 2 : We search the input_forest as a subtree of the reference_forest
  // On met en correspondance une sous-foret d'un arbre de F2 avec la foret F1
  //------------------------------------------------------------------------
  min=MAXDIST;
  cost2=getDBF(EMPTY_TREE,reference_vertex);
  for (i=1;i<=nj;i++)
    {
      int reference_child=T2->child(reference_vertex,i);
      dist2=getDBF(input_vertex,reference_child)-getDBF(EMPTY_TREE,reference_child);
      if (dist2<min) { min=dist2;jm=reference_child;}
    }
  cost2=cost2+min;

  if (cost2<=DIST) {DIST=cost2;MFC=2;}

  //---------------------------------------------------------------------------------------------
  //Case 3 : We evaluate the restricted mapping between the input_forest and the reference_forest
  // On evalue l'alignement restreint entre les deux forets
  //---------------------------------------------------------------------------------------------

  // On fabrique le graphe de flot necessaire a la resolution du probleme
  NodeList* input_list=new NodeList();
  NodeList* reference_list=new NodeList();
  for (int s1=1;s1<=ni;s1++) { input_list->push_back(T1->child(input_vertex,s1)); }
  for (int s2=1;s2<=nj;s2++) { reference_list->push_back(T2->child(reference_vertex,s2)); }
  _restrMapp.make(*input_list,*reference_list);
  _restrMappList.resize(ni+nj+3,EMPTY_NODE);

  // THE INPUT FOREST IS EMPTY_TREE
  // All the reference vertices are paired with empty
  // Si la foret initiale est vide, il faut inserer toutes les arbres de
  // la foret de reference et tous les noeuds de references sont associes
  // avec le noeud vide
  if (ni==0)
    {
      _restrMappList[1]=2;
      for (i=1;i<=nj;i++) { _restrMappList[i+1]=1; }
      cost3=getDBF(EMPTY_TREE,reference_vertex);
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
          for (i=1;i<=ni;i++) { _restrMappList[i]=ni+1; }
          cost3=getDBF(input_vertex,EMPTY_TREE);
        }
      else
        {
          //BOTH FOREST ARE NOT EMPTY_TREE
          // A retricted mapping must be calculated
          // Sinon on resout le probleme de flot maximum de cout minimum
          cost3=_restrMapp.minCostFlow(_restrMappList);

        }
    }


  if (cost3<DIST) { DIST=cost3; MFC=3; }

  if (!_spaceOptimization){
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
	  del=getDBF(input_vertex,EMPTY_TREE)-getDBF(im,EMPTY_TREE)+getDeBF(im,reference_vertex);
	  insert=getInBF(im,reference_vertex);
	  sub=getSuBF(im,reference_vertex);
	}
	break;
      case 2 :
	{
	  _choices.putLast(input_vertex,reference_vertex,jm);
	  del=getDeBF(input_vertex,jm);
	  insert=getDBF(EMPTY_TREE,reference_vertex)-getDBF(EMPTY_TREE,jm)+getInBF(input_vertex,jm);
	  sub=getSuBF(input_vertex,jm);
	}
	break;
      case 3 :
        {
          if (ni==0)
            {
              insert=getDBF(EMPTY_TREE,reference_vertex);
              del=0;
              sub=0;
            }
          else if (nj==0)
            {
              insert=0;
              del=getDBF(input_vertex,EMPTY_TREE);
              sub=0;
            }
          else
	    {
	      sub=0;
	      del=0;
	      for (int i=1;i<=T1->getNbChild(input_vertex);i++)
		{
		  _choices.putLast(input_vertex,reference_vertex,_restrMapp.who(_restrMappList[i]));
		  if (_restrMapp.who(_restrMappList[i])!=EMPTY_NODE)
		    {
		      sub=sub+getSuBT(T1->child(input_vertex,i),_restrMapp.who(_restrMappList[i]));
		      del=del+getDeBT(T1->child(input_vertex,i),_restrMapp.who(_restrMappList[i]));
		    }
		  else
		    {
		      del=del+getDBT(T1->child(input_vertex,i),EMPTY_TREE);
		    }
		}
	      insert=DIST-sub-del;
	    }
        }
        break;
      default :   break;
      }

    _deleteCost.putDBF(input_vertex,reference_vertex,del);
    _insertCost.putDBF(input_vertex,reference_vertex,insert);
    _substitutionCost.putDBF(input_vertex,reference_vertex,sub);
  }
  delete (NodeList*) input_list;
  delete (NodeList*) reference_list;
  _distances.putDBF(input_vertex,reference_vertex,DIST);
  return(DIST);

}



void Matching_U::getList(int input_vertex, int reference_vertex, Sequence* sequence)
{
  TreeList(input_vertex,reference_vertex,*sequence);
}

int Matching_U::Lat(ChoiceList* L, int vertex)
{
  ChoiceList::iterator begin;
  begin = L->begin();
  for (int i=0;i<vertex;i++)
    begin++;
  return(*begin);
}



void Matching_U::TreeList(int input_vertex,int reference_vertex,Sequence& sequence)
{
  if ((!T1->isNull())&&(!T2->isNull()))
    {
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

void Matching_U::ForestList(int input_vertex,int reference_vertex,Sequence& sequence)
{
  ChoiceList* L=_choices.getList(input_vertex,reference_vertex);
  int forest_choice=Lat(L,2);
  switch(forest_choice)
    {
    case 1: ForestList(Lat(L,3),reference_vertex,sequence);break;
    case 2: ForestList(input_vertex,Lat(L,3),sequence);break;
    case 3: {
      for (int i=1;i<=T1->getNbChild(input_vertex);i++)
        {
          int i_node=T1->child(input_vertex,i);
          int r_node=Lat(L,2+i);
          if (r_node!=-1) TreeList(i_node,r_node,sequence);
        }
    }break;
    default : break;
    }
}

DistanceType Matching_U::match()
{
  cerr<<"Unordered Matching"<<endl;
  const int size1 = T1->getNbVertex();
  const int size2 = T2->getNbVertex();
  DistanceType D=0;
  _distanceMatrix.resize(size1);
  //cerr << "" << "Already computed :     ";
  if ((!T1->isNull())&&(!T2->isNull()))
    {
      for (int input_vertex=size1-1;input_vertex>=0;input_vertex--)
	{
	  _distances.openDistancesVector(input_vertex);
	  // cout<<"dist "<<input_vertex<<endl;			if (int(100. - 100*input_vertex/size1)%10 == 0)
			  cerr << "\x0d" << " ok Already computed : "<<int(100. - 100*input_vertex/size1) <<"% " <<" matched ... " << flush;

	  if (!_spaceOptimization){
	    _insertCost.openDistancesVector(input_vertex);
	    //cout<<"ins "<<input_vertex<<endl;
	    _deleteCost.openDistancesVector(input_vertex);
	    //cout<<"del "<<input_vertex<<endl;
	    _substitutionCost.openDistancesVector(input_vertex);
	    //cout<<"sub "<<input_vertex<<endl;
	  }
	  _distanceMatrix[input_vertex].resize(size2);
	  for (int reference_vertex=size2-1;reference_vertex>=0;reference_vertex--)
	    {
	      distanceBetweenForest(input_vertex,reference_vertex);
	      DistanceType d =distanceBetweenTree(input_vertex,reference_vertex);
	      (_distanceMatrix[input_vertex])[reference_vertex]=d;
	    }
	  if (int(100. - 100*input_vertex/size1)%10 == 0)
	    cerr << "\x0d" << "ok Already computed : "<<int(100. - 100*input_vertex/size1) <<"% " <<" matched ... " << flush;
	  for (int i=1;i<=T1->getNbChild(input_vertex);i++)
	    {
	      _distances.closeDistancesVector(T1->child(input_vertex,i));
	      if (!_spaceOptimization){
		_insertCost.closeDistancesVector(T1->child(input_vertex,i));
		_deleteCost.closeDistancesVector(T1->child(input_vertex,i));
		_substitutionCost.closeDistancesVector(T1->child(input_vertex,i));
	      }
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
  //cerr<< "\x0d"<<endl;
  return(D);
}

DistanceType  Matching_U::getInsertCost()
{
  return(getInBT(0,0));
}

DistanceType  Matching_U::getDeleteCost()
{
  return(getDeBT(0,0));
}

DistanceType  Matching_U::getSubstitutionCost()
{
  return(getSuBT(0,0));
}

// --------------------------------------------
// Renvoie les distances entre arbres et forets
// --------------------------------------------
DistanceType Matching_U::getDistanceMatrix(int iv,int rv) const
{
  return(_distanceMatrix[iv][rv]);
}

DistanceType Matching_U::getDBT(int input_vertex,int reference_vertex) const
{
  return(_distances.getDBT(input_vertex,reference_vertex));
}

DistanceType Matching_U::getDBF(int input_vertex,int reference_vertex) const
{
  return(_distances.getDBF(input_vertex,reference_vertex));
}
DistanceType Matching_U::getInBT(int input_vertex,int reference_vertex) const
{
  return(_insertCost.getDBT(input_vertex,reference_vertex));
}

DistanceType Matching_U::getInBF(int input_vertex,int reference_vertex) const
{
  return(_insertCost.getDBF(input_vertex,reference_vertex));
}

DistanceType Matching_U::getDeBT(int input_vertex,int reference_vertex) const
{
  return(_deleteCost.getDBT(input_vertex,reference_vertex));
}



DistanceType Matching_U::getDeBF(int input_vertex,int reference_vertex) const
{
  return(_deleteCost.getDBF(input_vertex,reference_vertex));
}

DistanceType Matching_U::getSuBF(int input_vertex,int reference_vertex) const
{
  return(_substitutionCost.getDBF(input_vertex,reference_vertex));
}

DistanceType Matching_U::getSuBT(int input_vertex,int reference_vertex) const
{
  return(_substitutionCost.getDBT(input_vertex,reference_vertex));
}


// renvoie le dernier element de la liste de la case node du tableau maintenant les listes d'alignement
int Matching_U::M(int i_node,int r_node)
{
  return(_choices.getList(i_node,r_node)->back());
}









