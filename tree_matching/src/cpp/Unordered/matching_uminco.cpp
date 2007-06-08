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

#include "matching_uminco.h"

  // -------------
  // Constructeur
  // -------------
Matching_U_MinCo::Matching_U_MinCo(TreeGraph& input,TreeGraph& reference,NodeCost& nodeDistance):
Matching_U(input,reference,nodeDistance)
{
  _nbInputTreeConnectedComponents.make(*T1,*T2);
  _nbInputForestConnectedComponents.make(*T1,*T2);
  _nbReferenceTreeConnectedComponents.make(*T1,*T2);
  _nbReferenceForestConnectedComponents.make(*T1,*T2);
  _InputRootTable.make(*T1,*T2);
  _nbInputRootMapped.make(*T1,*T2);
  _ReferenceRootTable.make(*T1,*T2);
  _nbReferenceRootMapped.make(*T1,*T2);
}
  // -------------
  // Destructeur
  // -------------
Matching_U_MinCo::~Matching_U_MinCo()
{
}

  // ----------------------------------------------------------------------------------
  // Calcule la distance entre les deux arbres T1[input_vertex] et T2[reference_vertex]
  // ----------------------------------------------------------------------------------
DistanceType Matching_U_MinCo::distanceBetweenTree(int input_vertex,int reference_vertex)
{
  // On stocke dans ni et nj le nombre d'enfants de T1[i] et T2[j]
  int ni=T1->getNbChild(input_vertex);
  int nj=T2->getNbChild(reference_vertex);
  DistanceType cost1,cost2,cost3,dist1,dist2;
  DistanceType min,MIN=2*MAXDIST;
  int im=-1,jm=-1,MTC=0;
  int i;
  int nbcon,NBCON1=MAXINT,NBCON2=MAXINT,NBCON=MAXINT;
  int nbconref,NBCONREF1=MAXINT,NBCONREF2=MAXINT,NBCONREF=MAXINT;

//----------------------------------------------------------------------------------
//Case3 : We evaluate the matching between the input_forest and the reference_forest
// On evalue la mise en correspondance des arbres des deux forets issues de T1 et T2
//----------------------------------------------------------------------------------
// Le cout est celui de l'alignement des deux forets
// plus celui de l'echange de T1(i) en T2(j)


  cost3=getDBF(input_vertex,reference_vertex);
  cost3=cost3+_distances.getCCost(input_vertex,reference_vertex);

  nbcon=_nbInputForestConnectedComponents.getNBC(input_vertex,reference_vertex)
    -_nbInputRootMapped.getNBC(input_vertex,reference_vertex)
    +1;

  nbconref=_nbReferenceForestConnectedComponents.getNBC(input_vertex,reference_vertex)
    -_nbReferenceRootMapped.getNBC(input_vertex,reference_vertex)
    +1;

  // On conserve le cout s'il est inferieur au precedent
  if (cost3<MIN) 
    { 
      MIN=cost3; 
      MTC=3; 
      NBCON=nbcon;
      NBCONREF=nbconref;
    }
  


//----------------------------------------------------------------------
//Case 1 : We search the reference_tree as a subtree of the input_tree
//         On cherche a mettre en correspondance l'arbre de reference
//         avec un sous arbre de l'arbre initial, il faut donc effacer
//         T1 moins le sous arbre qui ressemble le plus a T2
//----------------------------------------------------------------------
  min=MAXDIST;
// cout de l'effacement de l'arbre initial
  cost1=getDBT(input_vertex,EMPTY_TREE);

  // nbcon et nbconref permettent de compter le nbre de composantes
  // connexes des arbres initial et final
  nbcon=0;
  nbconref=0;
  for (i=1;i<=ni;i++)
    {
      // On cherche parmi tous les fils de input_vertex celui dont l'arbre est le plus ressemblant a T2
      int input_child=T1->child(input_vertex,i);

      // la distance est donc le passage de T1[iam] en T2[j] - l'effacement de T[iam] qui a ete
      // compte precedemment
      dist1=getDBT(input_child,reference_vertex)-getDBT(input_child,EMPTY_TREE);

      //On compte le nbre de composantes connexes des deux arbres
      nbcon = getNBC(input_child,reference_vertex);
      nbconref = getNBCRef(input_child,reference_vertex);

      // On conserve la plus petite distance
      if (dist1<min)
	{ 
	  min=dist1; 
	  im=input_child; 
	  NBCON1=nbcon;
	  NBCONREF1=nbconref;
	}

      // Si la distance est egale, on conserve celui qui minimise le nbre
      // de composantes connexes de T1+T2
      if (dist1==min)
	{
	  if (nbcon+nbconref<NBCON1+nbconref<NBCONREF1)
	    {
	      im=input_child; 
	      NBCON1=nbcon;
	      NBCONREF1=nbconref;
	    }
	}
    }
  cost1=cost1+min;
  // On conserve le cout minimum
  if (cost1<MIN)
    {
      MIN=cost1;
      MTC=1; 
      NBCON=NBCON1;
      NBCONREF=NBCONREF1;
    }
  if (cost1==MIN)
    {
      if (NBCON1+NBCONREF1<NBCON+NBCONREF)
	{
	  MTC=1;
	  NBCON=NBCON1;
	  NBCONREF=NBCONREF1;
	}
    }

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

      // Calcul du nbre de c. c.   
      nbcon=getNBC(input_vertex,reference_child);
      nbconref=getNBCRef(input_vertex,reference_child);
      
      if (dist2<min)
	{ 
	  min=dist2;
	  jm=reference_child;
	  NBCON2=nbcon;
	  NBCONREF2=nbconref;
	}
      if (dist2==min)
	{
	  if (nbcon+nbconref<NBCONREF2+NBCON2)
	    {
	      jm=reference_child;
	      NBCON2=nbcon;
	      NBCONREF2=nbconref;
	    }
	}
    }           
  cost2=cost2+min; 
  // On conserve le cout s'il est inferieur au precedent
  if (cost2<MIN)
    {
      MIN=cost2; 
      MTC=2;
      NBCON=NBCON2;
      NBCONREF=NBCONREF2;
    }
  if (cost2==MIN)
    {

      if (NBCON2+NBCONREF2<NBCONREF+NBCON)
	{
	  MTC=2;
	  NBCON=NBCON2;
	  NBCONREF=NBCONREF2;
	}
    }

//-----------------------------------
// We maintain the matching lists
// mise a jour des listes d'alignement
//-----------------------------------
  switch (MTC)
    {
    case 1 :{
      _choices.putFirst(input_vertex,reference_vertex,im); 
      _choices.putLast(input_vertex,reference_vertex,-1);
      _InputRootTable.putNBC(input_vertex,reference_vertex,0);
      _ReferenceRootTable.putNBC(input_vertex,reference_vertex,1);
    };break;
    case 2 :{
      _choices.putFirst(input_vertex,reference_vertex,jm);
      _choices.putLast(input_vertex,reference_vertex,M(input_vertex,jm));
      _InputRootTable.putNBC(input_vertex,reference_vertex,1);
      _ReferenceRootTable.putNBC(input_vertex,reference_vertex,0);
    };break;
    case 3 :{
      _choices.putFirst(input_vertex,reference_vertex,-1);
      _choices.putLast(input_vertex,reference_vertex,reference_vertex);
      _InputRootTable.putNBC(input_vertex,reference_vertex,1);
      _ReferenceRootTable.putNBC(input_vertex,reference_vertex,1);
    };break;
    default : 	assert(0);break;
    }
  _choices.putFirst(input_vertex,reference_vertex,MTC);
  // On range dans le tableau des distances, la distance entre les arbres de racines input_vertex et
  // reference_vertex.
  _nbInputTreeConnectedComponents.putNBC(input_vertex,reference_vertex,NBCON);
  _nbReferenceTreeConnectedComponents.putNBC(input_vertex,reference_vertex,NBCONREF);
  _distances.putDBT(input_vertex,reference_vertex,MIN);

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
  return(MIN);

}



// -------------------------------
// Calcule la distance entre deux
// forets
// -------------------------------
DistanceType Matching_U_MinCo::distanceBetweenForest(int input_vertex,int reference_vertex)
{
// ni et nj representent le nombre de forets a comparees
  int ni=T1->getNbChild(input_vertex);
  int nj=T2->getNbChild(reference_vertex);
  DistanceType cost1,cost2,cost3,dist1,dist2;
  DistanceType min,DIST;
  int im=-1,jm=-1,MFC=0;
  int i;int cf,CF1=MAXINT,CF2=MAXINT,CF3=MAXINT,CF=MAXINT; 
  int sa,SA2=0,SA=0;
  int sar,SAR1=0,SAR=0;
//cf est le nbre de composantes connexes de la foret initiale
  // il est egal a la somme des composantes connexes de
  //chaque arbre
  int cfref,CFREF1=MAXINT,CFREF2=MAXINT,CFREF3=MAXINT,CFREF=MAXINT;
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
      cf=_nbInputForestConnectedComponents.getNBC(input_child,reference_vertex);
      cfref=_nbReferenceForestConnectedComponents.getNBC(input_child,reference_vertex);

      sar = _nbReferenceRootMapped.getNBC(input_child,reference_vertex);

      if (dist1<min) 
	{ 
	  min=dist1;
	  im=input_child;
	  CF1=cf;
	  CFREF1=cfref;
	  SAR1=sar;
	}
      if (dist1==min)
	{
	  if (cf+cfref-sar<CF1+cfref-SAR1)
	    {
	      im=input_child;
	      CF1=cf;
	      CFREF1=cfref;
	      SAR1=sar;
	    }
	}
    }
  cost1=cost1+min;

  if (cost1<DIST) 
    {
      DIST=cost1;
      MFC=1;
      CF=CF1;
      CFREF=CFREF1;
      SAR=SAR1;
      _nbInputRootMapped.putNBC(input_vertex,reference_vertex,0);
      _nbReferenceRootMapped.putNBC(input_vertex,reference_vertex,SAR1);
    }

  if (cost1==DIST) 
    {
      if (CF1-SAR1+CFREF1<CF+CFREF-SAR)
	{
	  MFC=1;
	  CF=CF1;
	  CFREF=CFREF1;
	  SAR=SAR1;
	  _nbInputRootMapped.putNBC(input_vertex,reference_vertex,0);
	  _nbReferenceRootMapped.putNBC(input_vertex,reference_vertex,SAR1);
	}
    }

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
      cf=_nbInputForestConnectedComponents.getNBC(input_vertex,reference_child);
      cfref=_nbReferenceForestConnectedComponents.getNBC(input_vertex,reference_child);
      sa=_nbInputRootMapped.getNBC(input_vertex,reference_child);

      if (dist2<min) 
	{ 
	  min=dist2;
	  jm=reference_child;
	  CF2=cf;
	  CFREF2=cfref;
	  SA2=sa;
	}
      if (dist2==min)
	{
	  if ((cf+cfref-sa<CF2-SA2+CFREF2))
	    {
	      jm=reference_child;
	      CF2=cf;
	      CFREF2=cfref;
	      SA2=sa;
	    }
	}
    }
  cost2=cost2+min;

  if (cost2<DIST)
    {
      DIST=cost2;
      MFC=2;
      CF=CF2;
      CFREF=CFREF2;
      _nbInputRootMapped.putNBC(input_vertex,reference_vertex,SA2);
      _nbInputRootMapped.putNBC(input_vertex,reference_vertex,0);
      SA=SA2;
    }
  if (cost2==DIST) 
    {
      if ((CF2-SA2+CFREF2<CF-SA+CFREF))
	{
	  MFC=2;
	  CF=CF2;
	  CFREF=CFREF2;
	  SA=SA2;
	  _nbInputRootMapped.putNBC(input_vertex,reference_vertex,SA2);
	  _nbReferenceRootMapped.putNBC(input_vertex,reference_vertex,0);
	}
    }
//---------------------------------------------------------------------------------------------
//Case 3 : We evaluate the restricted mapping between the input_forest and the reference_forest
// On evalue l'alignement restreint entre les deux forets
//---------------------------------------------------------------------------------------------

// On fabrique le graphe de flot necessaire a la resolution du probleme
  NodeList* input_list=new NodeList();
  NodeList* reference_list=new NodeList();
  int nbInputMappedChild=0;
  int nbReferenceMappedChild=0;
  for (int s1=1;s1<=ni;s1++) { input_list->push_back(T1->child(input_vertex,s1)); };
  for (int s2=1;s2<=nj;s2++) { reference_list->push_back(T2->child(reference_vertex,s2)); };
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
      for (i=1;i<=nj;i++) { _restrMappList[i+1]=1; };
      cost3=getDBF(EMPTY_TREE,reference_vertex);
      CF3=0;
      CFREF3=0;
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
	  for (i=1;i<=ni;i++) { _restrMappList[i]=ni+1; };
	  cost3=getDBF(input_vertex,EMPTY_TREE); 
	  CF3=0;
	  CFREF3=0;
	}
      else
	{
	  //BOTH FOREST ARE NOT EMPTY_TREE
// A retricted mapping must be calculated
// Sinon on resout le probleme de flot maximum de cout minimum
	  cost3=_restrMapp.minCostFlow(_restrMappList);
	  cf=0;
	  cfref=0;
	  for (int i=1;i<=T1->getNbChild(input_vertex);i++)
	    {
	      int input_child=T1->child(input_vertex,i);
	      if (_restrMapp.who(_restrMappList[i])!=EMPTY_NODE)
		{
		  cf=cf+getNBC(input_child,_restrMapp.who(_restrMappList[i]));
		  cfref=cfref+getNBCRef(input_child,_restrMapp.who(_restrMappList[i]));
		  
		  nbInputMappedChild=nbInputMappedChild 
		    + _InputRootTable.getNBC(input_child,_restrMapp.who(_restrMappList[i]));

		  nbReferenceMappedChild=nbReferenceMappedChild 
		    + _ReferenceRootTable.getNBC(input_child,_restrMapp.who(_restrMappList[i]));
		}
	    }
	  CFREF3=cfref;
	  CF3=cf;
	}
    }

  if (cost3<DIST)
    { 
      DIST=cost3; 
      MFC=3;
      CF=CF3;
      CFREF=CFREF3;
      _nbInputRootMapped.putNBC(input_vertex,reference_vertex,nbInputMappedChild);
      _nbReferenceRootMapped.putNBC(input_vertex,reference_vertex,nbReferenceMappedChild);
      
    }
  if (cost3==DIST)
    { 
      if (CF3+CFREF3-nbInputMappedChild-nbReferenceMappedChild<=CF+CFREF-SA-SAR)
      MFC=3;
      CF=CF3;
      CFREF=CFREF3;
      _nbInputRootMapped.putNBC(input_vertex,reference_vertex,nbInputMappedChild);
      _nbReferenceRootMapped.putNBC(input_vertex,reference_vertex,nbReferenceMappedChild);
      
    }

//-------------------------------
//We maintain the matching lists
// On maintient les listes d'alignement
//-------------------------------
  _choices.createList(input_vertex,reference_vertex);
  _choices.putFirst(input_vertex,reference_vertex,MFC);

  _nbInputForestConnectedComponents.putNBC(input_vertex,reference_vertex,CF);
  _nbReferenceForestConnectedComponents.putNBC(input_vertex,reference_vertex,CFREF);
	
  switch(MFC)
    {
    case 1 : _choices.putLast(input_vertex,reference_vertex,im);break;
    case 2 : _choices.putLast(input_vertex,reference_vertex,jm);break;
    case 3 : {
      for (int i=1;i<=T1->getNbChild(input_vertex);i++)
	{
	  _choices.putLast(input_vertex,reference_vertex,_restrMapp.who(_restrMappList[i]));
	};
    }; break;
    default : 	break;
    };
  
	
  delete (NodeList*) input_list;
  delete (NodeList*) reference_list;
  _distances.putDBF(input_vertex,reference_vertex,DIST);
  return(DIST);

}

DistanceType Matching_U_MinCo::match()
{
  cerr<<"Minimum connected Components"<<endl;
  const int size1 = T1->getNbVertex();
  const int size2 = T2->getNbVertex();
  DistanceType D=0;
  if ((!T1->isNull())&&(!T2->isNull()))
    {
      for (int input_vertex=T1->getNbVertex()-1;input_vertex>=0;input_vertex--)
	{
	  _distances.openDistancesVector(input_vertex);
	  _InputRootTable.openIntVector(input_vertex);
	  _ReferenceRootTable.openIntVector(input_vertex);
	  _nbInputTreeConnectedComponents.openIntVector(input_vertex);
	  _nbInputForestConnectedComponents.openIntVector(input_vertex);
	  _nbReferenceTreeConnectedComponents.openIntVector(input_vertex);
	  _nbReferenceForestConnectedComponents.openIntVector(input_vertex);
	  _nbInputRootMapped.openIntVector(input_vertex);
	  _nbReferenceRootMapped.openIntVector(input_vertex);
	  //On initialise le tableau _InputRootTable a False, aucun noeud n'a d'images
	  //au debut de l'algo.
	  for (int j=0;j<T2->getNbVertex();j++)
	    {
	      _InputRootTable.putNBC(input_vertex,j,0);
	      _ReferenceRootTable.putNBC(input_vertex,j,0);
	    }

	  for (int reference_vertex=T2->getNbVertex()-1;reference_vertex>=0;reference_vertex--)
	    {
	      distanceBetweenForest(input_vertex,reference_vertex);
	      distanceBetweenTree(input_vertex,reference_vertex);
	    }
	  if (int(100. - 100*input_vertex/size1)%10 == 0)
	    cerr << "\x0d" << "Already computed : "<<int(100. - 100*input_vertex/size1) <<"% " <<" matched ... " << flush;
	  
	  for (int i=1;i<=T1->getNbChild(input_vertex);i++)
	    {
	      _distances.closeDistancesVector(T1->child(input_vertex,i));
	      _nbInputTreeConnectedComponents.closeIntVector(T1->child(input_vertex,i));
	      _nbInputForestConnectedComponents.closeIntVector(T1->child(input_vertex,i));
	      _InputRootTable.closeIntVector(T1->child(input_vertex,i));
	      _nbInputRootMapped.closeIntVector(T1->child(input_vertex,i));
	      _nbReferenceTreeConnectedComponents.closeIntVector(T1->child(input_vertex,i));
	      _nbReferenceForestConnectedComponents.closeIntVector(T1->child(input_vertex,i));
	      _ReferenceRootTable.closeIntVector(T1->child(input_vertex,i));
	      _nbReferenceRootMapped.closeIntVector(T1->child(input_vertex,i));

	    }
	};
      
      D=getDBT(0,0);
      /*for (int compt=0;compt<3;compt++)
	{
	  cout<<"Cas "<<compt+1<<" : "<<_sumNbCaseVector[compt]<<endl;
	}
      int tot=0;
      for (int compt=0;compt<7;compt++)
	{
	  tot=tot+_nbCaseVector[compt];
	}
      cout<<" Total = "<< tot<<endl<<endl;
      for (int compt=0;compt<7;compt++)
	{
	  cout<<"Pourcentage Cas "<<compt+1<<" : "<<double((_nbCaseVector[compt]))/double(tot)*100.<<endl;
	}
      cout<<endl;
      cout<<"nbre de composantes connexes sur T1 : "<<getNBC(0,0)<<endl;
      cout<<"nbre de composantes connexes sur T2 : "<<getNBCRef(0,0)<<endl;
      */
      
    }
  else
    {
      if (T1->isNull())
	{
	  if (!T2->isNull()) {D=_distances.referenceTreeFromEmpty(0);};
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

int Matching_U_MinCo::getNBC(int input_vertex,int reference_vertex) const
{
  return(_nbInputTreeConnectedComponents.getNBC(input_vertex,reference_vertex));
}

int Matching_U_MinCo::getNBCRef(int input_vertex,int reference_vertex) const
{
  return(_nbReferenceTreeConnectedComponents.getNBC(input_vertex,reference_vertex));
}








